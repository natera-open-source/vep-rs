#!/usr/bin/env python3
"""Hold the documented public interface of the workspace's crates to the record under `tests/public-interface/`.

The record has one file a package, `<package>.txt`, which this script writes and nobody edits by hand. A
file opens with one comment line, then states the package's library target and its cargo features, then
holds one line for each documented public item of the library, tab-separated and sorted by path, kind and
signature:

    library <the library's crate name>, or `none` for a package with no library target
    feature <name>  <what the feature enables, sorted and blank-separated, or `-`>
    <kind>  <path>  <signature>

The feature named `default` is the package's default feature set.

A documented public item is one the API documentation shows: it is reached from the crate root through
public modules and `pub use`, and neither it nor a module, type, trait, variant, impl block or `pub use` on
the way to it carries `#[doc(hidden)]`. An item under `#[doc(hidden)]` has no line, so an item that moves
there reads as removed. The kinds, with what the signature holds:

    module                  `mod`
    struct, union           generic parameters and where clause; whether the fields are named or
                            positional and whether one is private; `#[non_exhaustive]`, `#[repr]`
    enum                    generic parameters and where clause; `#[non_exhaustive]`
    variant                 its shape and an explicit discriminant
    field                   its type; a positional field is named by its index
    trait, trait alias      `unsafe`, `auto`, generic parameters, supertraits, where clause and whether a
                            `dyn` of it can be written
    function, method        `const`, `async`, `unsafe`, the ABI, generic parameters, each parameter with
                            its name and type, the return type and the where clause; a method of a trait
                            says whether the trait provides a body; a method of an impl block that has
                            generic parameters or a where clause names the block
    constant, static        its type; for a static also `mut` and, in an `extern` block, `unsafe`
    associated constant     its type; in a trait, whether the trait gives a value
    associated type         in a trait, its bounds and default
    type alias              generic parameters and the type it stands for
    macro                   its rules as rustdoc prints them, each pattern without its body
    re-export               the path a `pub use` names, for an item that has a documented path of its own
                            or belongs to another crate; the crate a `pub extern crate` names
    implementation          one line for each trait a type implements by hand or by a derive, with the
                            associated types and constants the block sets; listed at the type, or at a
                            trait of this crate where the block is written for a type of another crate
                            or for a type alias
    auto traits             one line a type: which of `Send`, `Sync`, `Unpin`, `UnwindSafe` and
                            `RefUnwindSafe` the compiler derives for it (`!` before one it does not),
                            each with its conditions

Every type in a signature is written by the path of its definition, or by its shortest documented path
where it belongs to a crate of this workspace, so a changed `use` line in a source file changes no line.
`#[deprecated]` is held on every item, and so is a `#[cfg(..)]` the item itself carries: the listing is made
with every feature on, so an item that moves behind a feature keeps its line and gains that mark. Not held: the value of a constant, the body of a function, the
implicit discriminants and the order of an enum's variants, blanket implementations, the marker trait a
derived `PartialEq` adds, the impl blocks rustdoc's output leaves out (the methods written for a `dyn` of a
trait, and another crate's trait implemented for such a `dyn` or through a type alias), and the meaning of a
feature beyond the names it enables.

Beside the listings stands `outside-crates.tsv`, kept by hand: one line for each crate outside the workspace
and the standard library whose types a documented signature may name, `<crate><tab><reason>`, the crate named
as rustdoc names the one that defines the type. A signature that names a type of any other crate is a
finding, and so is a line of the list whose crate no signature names.

How the listing is made: `cargo doc --workspace --lib --no-deps --locked --all-features` writes rustdoc's
JSON output for every library, with hidden items included so that an item re-exported out of a hidden
module is read whole. That output is unstable: on a stable toolchain rustdoc gives it only under
`RUSTC_BOOTSTRAP=1`, which this script sets for the rustdoc processes alone, through a wrapper it hands
cargo as `RUSTDOC`, so that no build script and no dependency is compiled under it and the build shares
every dependency with an ordinary `cargo doc`. The script reads format 57, the format of Rust 1.97.1, and
names any other format it is given on standard error; a shape it does not know ends the run with exit 2
rather than leave an item out.

Every finding is one line on standard error, `ERROR: [public-interface] <package>: <what changed>`, followed
by one line that gives the command that writes the record again; the exit status is then 1. A record that
holds prints `check_public_interface: OK` on standard output.

usage: check_public_interface.py [--update] [--offline] [--root DIR] [--record DIR]

`--update` writes the record again from the tree (the diff is then part of the change), and leaves
`outside-crates.tsv` as it is. `--offline` is handed to cargo. `--root` is the workspace to read (default:
the tree this script stands in) and `--record` its record (default: `tests/public-interface` under the
root).

No dependency beyond the standard library and cargo.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RECORD = Path("tests") / "public-interface"
OUTSIDE_CRATES = "outside-crates.tsv"
TESTED_FORMAT = 57
STANDARD_LIBRARY = frozenset(("std", "core", "alloc", "proc_macro", "test"))
LIBRARY_KINDS = frozenset(("lib", "rlib", "dylib", "cdylib", "staticlib", "proc-macro"))
HIDDEN_ATTRIBUTE = "#[doc(hidden)]"
# The auto traits a caller's code can ask of a type, by the path of their definition.
AUTO_TRAITS = {
    "core::marker::Send": "Send",
    "core::marker::Sync": "Sync",
    "core::marker::Unpin": "Unpin",
    "core::panic::unwind_safe::UnwindSafe": "UnwindSafe",
    "core::panic::unwind_safe::RefUnwindSafe": "RefUnwindSafe",
}
# rustdoc writes a `#[cfg(..)]` as the compiler's own record of it, each part with the place it was read at.
CFG_ATTRIBUTE = "#[attr = CfgTrace(["
CFG_PLACE = re.compile(r",\s*(?:span:\s*)?[^\s,()\[\]{}]+:\d+:\d+: \d+:\d+ \(#\d+\)")
CFG_NAME_VALUE = re.compile(r'NameValue \{ name: "([^"]*)", value: Some\("([^"]*)"\) \}')
CFG_NAME = re.compile(r'NameValue \{ name: "([^"]*)", value: None \}')
# Implemented by the compiler beside a derived `PartialEq`; no stable code can name it.
UNLISTED_TRAITS = frozenset(("core::marker::StructuralPartialEq",))
UPDATE_COMMAND = "python3 scripts/check_public_interface.py --update"
HEADER = (
    "# Written by scripts/check_public_interface.py --update, which reads the tree; "
    "the script's opening comment gives the format."
)
RUSTDOC_WRAPPER = '#!/bin/sh\nRUSTC_BOOTSTRAP=1 exec "${PUBLIC_INTERFACE_RUSTDOC:-rustdoc}" "$@"\n'
RUSTDOC_FLAGS = ("-Z", "unstable-options", "--output-format", "json", "--document-hidden-items")

# The kind a module's direct child is listed under, by the key of its `inner`.
MODULE_CHILD_KINDS = {
    "module": "module",
    "struct": "struct",
    "enum": "enum",
    "union": "union",
    "trait": "trait",
    "trait_alias": "trait alias",
    "function": "function",
    "constant": "constant",
    "static": "static",
    "type_alias": "type alias",
    "macro": "macro",
    "proc_macro": "macro",
    "extern_crate": "re-export",
}

# The kinds of item whose `impls` list the impl blocks written for the type.
OWNER_KINDS = frozenset(("struct", "enum", "union"))

# The kind a member of a trait or of an inherent impl block is listed under.
MEMBER_KINDS = {
    "function": "method",
    "assoc_const": "associated constant",
    "assoc_type": "associated type",
}

TYPE, VALUE, MACRO = "type", "value", "macro"
ALL_NAMESPACES = frozenset((TYPE, VALUE, MACRO))
NAMESPACES = {
    "module": frozenset((TYPE,)),
    "enum": frozenset((TYPE,)),
    "union": frozenset((TYPE,)),
    "trait": frozenset((TYPE,)),
    "trait_alias": frozenset((TYPE,)),
    "type_alias": frozenset((TYPE,)),
    "extern_crate": frozenset((TYPE,)),
    "struct": frozenset((TYPE, VALUE)),
    "variant": frozenset((TYPE, VALUE)),
    "function": frozenset((VALUE,)),
    "constant": frozenset((VALUE,)),
    "static": frozenset((VALUE,)),
    "macro": frozenset((MACRO,)),
    "proc_macro": frozenset((MACRO,)),
}


class Unreadable(Exception):
    """The tree or rustdoc's output holds something this script cannot read."""


def inner_kind(item: dict) -> str:
    """The key of an item's `inner`, which names what the item is; a kind with no content is the name alone."""
    inner = item["inner"]
    return inner if isinstance(inner, str) else next(iter(inner))


def carries_hidden(item: dict) -> bool:
    """Whether the item itself carries `#[doc(hidden)]`."""
    for attribute in item.get("attrs", ()):
        if isinstance(attribute, dict) and attribute.get("other") == HIDDEN_ATTRIBUTE:
            return True
    return False


class Entry:
    """One listed item: its kind and path, whether it is hidden, and what its signature is read from."""

    __slots__ = ("block", "blocks", "hidden", "in_trait", "item", "kind", "path")

    def __init__(self, kind, path, hidden, item, block=None, in_trait=False, blocks=()):
        self.kind = kind
        self.path = path
        self.hidden = hidden
        self.item = item
        # The inherent impl block a method or an associated constant stands in.
        self.block = block
        self.in_trait = in_trait
        # The synthetic impl blocks of an `auto traits` line.
        self.blocks = blocks


class Walker:
    """Walks one rustdoc JSON document from the crate root and collects the items a caller can name."""

    def __init__(self, document: dict):
        self.index = document["index"]
        root = self.index.get(str(document["root"]))
        if root is None or "module" not in root["inner"]:
            raise Unreadable("rustdoc's output names no crate root module")
        self.root = root
        # Item id -> one hidden flag for each path of its own, through public modules alone.
        self.own_paths: dict[int, list[bool]] = {}
        self.entries: list[Entry] = []
        self.recording_own_paths = True
        # Ids of the re-exported items and globbed modules being listed, so that a `pub use` that leads
        # back to one of them ends the walk there.
        self.open_targets: list[int] = []
        self.listed_blocks: set[int] = set()

    def run(self) -> list[Entry]:
        """The entries of the crate, hidden ones included."""
        self.recording_own_paths = True
        self.list_item(self.root, self.root["name"], False, "module")
        self.recording_own_paths = False
        self.entries = []
        self.list_item(self.root, self.root["name"], False, "module")
        return self.entries

    def item(self, item_id) -> dict | None:
        """The item an id names, or None when the document does not hold it."""
        if item_id is None:
            return None
        return self.index.get(str(item_id))

    def add(self, kind, path, hidden, item, **context) -> None:
        """Record one entry or, on the first pass, one path of the item's own."""
        if self.recording_own_paths:
            if kind not in ("implementation", "auto traits"):
                self.own_paths.setdefault(item["id"], []).append(hidden)
            return
        self.entries.append(Entry(kind, path, hidden, item, **context))

    def list_item(self, item, path, inherited_hidden, kind) -> None:
        """List one item at `path`, with every member it holds."""
        hidden = inherited_hidden or carries_hidden(item)
        self.add(kind, path, hidden, item)
        what = inner_kind(item)
        inner = item["inner"][what]
        if what == "module":
            self.list_module_children(inner["items"], path, hidden)
        elif what == "struct":
            self.list_fields(self.field_ids(inner["kind"]), path, hidden, public_only=True)
            self.list_blocks(item, inner["impls"], path, hidden)
        elif what == "union":
            self.list_fields(inner["fields"], path, hidden, public_only=True)
            self.list_blocks(item, inner["impls"], path, hidden)
        elif what == "enum":
            for variant_id in inner["variants"]:
                variant = self.item(variant_id)
                if variant is not None:
                    self.list_item(variant, f"{path}::{variant['name']}", hidden, "variant")
            self.list_blocks(item, inner["impls"], path, hidden)
        elif what == "trait":
            self.list_members(inner["items"], path, hidden, public_only=False, in_trait=True)
            self.list_implementations_for_other_types(inner["implementations"], path, hidden)
        elif what == "variant":
            self.list_fields(self.field_ids(inner["kind"]), path, hidden, public_only=False)

    @staticmethod
    def field_ids(kind) -> list:
        """The field ids of a struct's or a variant's `kind`."""
        if isinstance(kind, str):
            return []
        if "tuple" in kind:
            return kind["tuple"]
        if "plain" in kind:
            return kind["plain"]["fields"]
        if "struct" in kind:
            return kind["struct"]["fields"]
        raise Unreadable(f"unknown struct or variant shape: {sorted(kind)}")

    def list_fields(self, field_ids, path, hidden, public_only) -> None:
        """List fields; rustdoc writes null for a field it leaves out."""
        for field_id in field_ids:
            field = self.item(field_id)
            if field is None:
                continue
            if public_only and field["visibility"] != "public":
                continue
            self.add("field", f"{path}::{field['name']}", hidden or carries_hidden(field), field)

    def list_blocks(self, owner, block_ids, path, hidden) -> None:
        """List what a type's impl blocks give it: members, implemented traits and auto traits."""
        synthetic = []
        for block_id in block_ids:
            block = self.item(block_id)
            if block is None:
                continue
            inner = block["inner"]["impl"]
            block_hidden = hidden or carries_hidden(block)
            if inner["trait"] is None:
                self.list_members(inner["items"], path, block_hidden, public_only=True, block=block)
            elif inner["is_synthetic"]:
                synthetic.append(block)
            elif inner["blanket_impl"] is None:
                self.listed_blocks.add(block["id"])
                self.add("implementation", path, block_hidden, block)
        if synthetic:
            self.add("auto traits", path, hidden, owner, blocks=tuple(synthetic))

    def list_implementations_for_other_types(self, block_ids, path, hidden) -> None:
        """List a trait's implementations for types that are no struct, enum or union of this crate."""
        for block_id in block_ids:
            block = self.item(block_id)
            if block is None or block.get("crate_id", 0) != 0:
                continue
            inner = block["inner"]["impl"]
            if inner["is_synthetic"] or inner["blanket_impl"] is not None:
                continue
            target = inner["for"]
            if "resolved_path" in target:
                owner = self.item(target["resolved_path"]["id"])
                # A block written for a type alias stands in the `impls` of no type.
                if owner is not None and owner.get("crate_id", 0) == 0 and inner_kind(owner) in OWNER_KINDS:
                    continue
            self.add("implementation", path, hidden or carries_hidden(block), block)

    def list_members(self, member_ids, path, hidden, public_only, block=None, in_trait=False) -> None:
        """List the members of a trait or of an inherent impl block."""
        for member_id in member_ids:
            member = self.item(member_id)
            if member is None:
                continue
            if public_only and member["visibility"] != "public":
                continue
            what = inner_kind(member)
            if what not in MEMBER_KINDS:
                raise Unreadable(f"unknown member kind {what!r} at {path}::{member['name']}")
            self.add(
                MEMBER_KINDS[what],
                f"{path}::{member['name']}",
                hidden or carries_hidden(member),
                member,
                block=block,
                in_trait=in_trait,
            )

    def list_module_children(self, child_ids, path, hidden) -> None:
        """List what the module at `path` holds."""
        children = [child for child in map(self.item, child_ids) if child is not None]
        for child in children:
            if child["visibility"] != "public":
                continue
            if inner_kind(child) == "use":
                self.list_use(child, path, hidden, children, None)
            else:
                self.list_named(child, f"{path}::{child['name']}", hidden)

    def list_named(self, item, path, hidden) -> None:
        """List an item a module holds, or one a `pub use` names, at `path`."""
        what = inner_kind(item)
        if what == "variant":
            self.list_item(item, path, hidden, "variant")
            return
        if what not in MODULE_CHILD_KINDS:
            raise Unreadable(f"unknown item kind {what!r} at {path}")
        self.list_item(item, path, hidden, MODULE_CHILD_KINDS[what])

    def list_use(self, use, path, inherited_hidden, siblings, skip) -> None:
        """List what one `pub use` adds to the module at `path`.

        `siblings` are the items declared beside the `pub use`; `skip` says which names an enclosing glob
        leaves out, or is None outside a glob.
        """
        if self.recording_own_paths:
            return
        inner = use["inner"]["use"]
        hidden = inherited_hidden or carries_hidden(use)
        target = self.item(inner["id"])
        if target is not None and target.get("crate_id", 0) != 0:
            target = None
        if inner["is_glob"]:
            self.list_glob(use, target, path, hidden, siblings, skip)
            return
        name = inner["name"]
        namespaces = NAMESPACES.get(inner_kind(target), ALL_NAMESPACES) if target else ALL_NAMESPACES
        if skip is not None and skip(name, namespaces):
            return
        self.list_reexported(use, target, f"{path}::{name}", hidden)

    def list_reexported(self, use, target, path, hidden) -> None:
        """List `target` as the `pub use` item `use` names it at `path`.

        An item with a documented path of its own, and an item of another crate, is one `re-export`
        entry. An item with none (it stands in a private module, or in a hidden one) is listed as the item
        itself at the path the `pub use` gives it, with its members, where rustdoc documents it.
        """
        if target is None or target["id"] in self.open_targets:
            self.add("re-export", path, hidden, use)
            return
        own = self.own_paths.get(target["id"])
        target_hidden = carries_hidden(target)
        if own is not None and (target_hidden or not all(own)):
            self.add("re-export", path, hidden or target_hidden, use)
            return
        self.open_targets.append(target["id"])
        try:
            self.list_named(target, path, hidden)
        finally:
            self.open_targets.pop()

    def list_glob(self, use, target, path, hidden, siblings, outer_skip) -> None:
        """List what the glob `use` of `target` adds to the module at `path`."""
        if target is None:
            self.add("re-export", f"{path}::*", hidden, use)
            return
        if target["id"] in self.open_targets:
            return
        declared = self.declared_names(siblings)

        def skip(name, namespaces):
            if any((name, namespace) in declared for namespace in namespaces):
                return True
            return outer_skip is not None and outer_skip(name, namespaces)

        what = inner_kind(target)
        if what == "module":
            member_ids = target["inner"]["module"]["items"]
        elif what == "enum":
            member_ids = target["inner"]["enum"]["variants"]
        else:
            raise Unreadable(f"a glob of a {what} at {path}")
        members = [member for member in map(self.item, member_ids) if member is not None]
        self.open_targets.append(target["id"])
        try:
            for member in members:
                member_kind = inner_kind(member)
                if member_kind != "variant" and member["visibility"] != "public":
                    continue
                if member_kind == "use":
                    self.list_use(member, path, hidden, members, skip)
                elif not skip(member["name"], NAMESPACES.get(member_kind, ALL_NAMESPACES)):
                    self.list_reexported(use, member, f"{path}::{member['name']}", hidden)
        finally:
            self.open_targets.pop()

    def declared_names(self, siblings) -> set:
        """The (name, namespace) pairs a module declares without a glob."""
        declared = set()
        for sibling in siblings:
            what = inner_kind(sibling)
            if what == "use":
                inner = sibling["inner"]["use"]
                if inner["is_glob"]:
                    continue
                target = self.item(inner["id"])
                namespaces = NAMESPACES.get(inner_kind(target), ALL_NAMESPACES) if target else ALL_NAMESPACES
                name = inner["name"]
            else:
                namespaces = NAMESPACES.get(what, ALL_NAMESPACES)
                name = sibling["name"]
            declared.update((name, namespace) for namespace in namespaces)
        return declared


def one_line(text: str) -> str:
    """The text with every run of white space as one blank."""
    return " ".join(str(text).split())


def cfg_mark(attribute: str) -> str:
    """The `#[cfg(..)]` an item carries, as its source writes the conditions."""
    text = CFG_PLACE.sub("", one_line(attribute[len(CFG_ATTRIBUTE) : -len("])]")]))
    text = CFG_NAME.sub(r"\1", CFG_NAME_VALUE.sub(r'\1 = "\2"', text))
    text = text.replace("All([", "all(").replace("Any([", "any(").replace("Not(", "not(").replace("])", ")")
    if "(#" in text or "{" in text or "[" in text:
        raise Unreadable(f"a cfg attribute in a shape this script does not read: {attribute}")
    return f"#[cfg({text})]"


class Signatures:
    """Writes the signature of each entry of one crate, and notes the crates each signature names."""

    def __init__(self, document: dict, documented_paths: dict[str, str]):
        self.index = document["index"]
        self.paths = document["paths"]
        self.crates = document["external_crates"]
        self.crate_name = self.index[str(document["root"])]["name"]
        # The definition path of each documented item of the workspace -> its shortest documented path.
        self.documented_paths = documented_paths
        self.named: set[str] = set()

    # Paths

    def definition(self, item_id) -> tuple[str, str] | None:
        """The crate and the definition path of the item an id names, as rustdoc's `paths` table has them."""
        entry = self.paths.get(str(item_id))
        if entry is None:
            return None
        crate_id = entry["crate_id"]
        crate = self.crate_name if crate_id == 0 else self.crates[str(crate_id)]["name"]
        return crate, "::".join(entry["path"])

    def named_path(self, item_id, written: str) -> str:
        """The path a signature gives the item an id names, noting its crate; `written` where none is known."""
        definition = self.definition(item_id)
        if definition is None:
            return written
        crate, path = definition
        self.named.add(crate)
        return self.documented_paths.get(path, path)

    def path(self, path: dict) -> str:
        """A path to a type or a trait, with its generic arguments."""
        return self.named_path(path["id"], path["path"]) + self.generic_args(path.get("args"))

    # Types

    def type(self, value) -> str:
        """A type as a signature writes it."""
        if value == "infer":
            return "_"
        if not isinstance(value, dict) or len(value) != 1:
            raise Unreadable(f"unknown type shape: {value!r}")
        what, inner = next(iter(value.items()))
        if what == "resolved_path":
            return self.path(inner)
        if what == "generic":
            return inner
        if what == "primitive":
            return inner
        if what == "tuple":
            parts = [self.type(part) for part in inner]
            return "(" + ", ".join(parts) + ("," if len(parts) == 1 else "") + ")"
        if what == "slice":
            return f"[{self.type(inner)}]"
        if what == "array":
            return f"[{self.type(inner['type'])}; {one_line(inner['len'])}]"
        if what == "borrowed_ref":
            lifetime = f"{inner['lifetime']} " if inner["lifetime"] else ""
            mutable = "mut " if inner["is_mutable"] else ""
            return f"&{lifetime}{mutable}{self.type(inner['type'])}"
        if what == "raw_pointer":
            return ("*mut " if inner["is_mutable"] else "*const ") + self.type(inner["type"])
        if what == "impl_trait":
            return "impl " + self.bounds(inner)
        if what == "dyn_trait":
            parts = [self.poly_trait(part) for part in inner["traits"]]
            if inner["lifetime"]:
                parts.append(inner["lifetime"])
            return "dyn " + " + ".join(parts)
        if what == "function_pointer":
            binder = self.binder(inner["generic_params"])
            return binder + self.header(inner["header"]) + "fn" + self.parameters_and_output(inner["sig"])
        if what == "qualified_path":
            name = inner["name"] + self.generic_args(inner.get("args"))
            self_type = self.type(inner["self_type"])
            if inner["trait"] is None or not inner["trait"].get("path"):
                return f"{self_type}::{name}"
            return f"<{self_type} as {self.path(inner['trait'])}>::{name}"
        if what == "pat":
            return f"{self.type(inner['type'])} is {one_line(inner['__pat_unstable_do_not_use'])}"
        raise Unreadable(f"unknown type shape: {what!r}")

    def generic_args(self, args) -> str:
        """The generic arguments of a path."""
        if args is None:
            return ""
        if args == "return_type_notation":
            return "(..)"
        what, inner = next(iter(args.items()))
        if what == "angle_bracketed":
            parts = [self.generic_arg(arg) for arg in inner["args"]]
            parts += [self.constraint(constraint) for constraint in inner["constraints"]]
            return f"<{', '.join(parts)}>" if parts else ""
        if what == "parenthesized":
            output = f" -> {self.type(inner['output'])}" if inner["output"] is not None else ""
            return "(" + ", ".join(self.type(part) for part in inner["inputs"]) + ")" + output
        raise Unreadable(f"unknown generic arguments shape: {what!r}")

    def generic_arg(self, arg) -> str:
        if arg == "infer":
            return "_"
        what, inner = next(iter(arg.items()))
        if what == "lifetime":
            return inner
        if what == "type":
            return self.type(inner)
        if what == "const":
            return self.constant_expression(inner)
        raise Unreadable(f"unknown generic argument shape: {what!r}")

    @staticmethod
    def constant_expression(constant: dict) -> str:
        return one_line(constant["value"] if constant.get("value") is not None else constant["expr"])

    def term(self, term) -> str:
        what, inner = next(iter(term.items()))
        if what == "type":
            return self.type(inner)
        if what == "constant":
            return self.constant_expression(inner)
        raise Unreadable(f"unknown term shape: {what!r}")

    def constraint(self, constraint: dict) -> str:
        name = constraint["name"] + self.generic_args(constraint.get("args"))
        what, inner = next(iter(constraint["binding"].items()))
        if what == "equality":
            return f"{name} = {self.term(inner)}"
        if what == "constraint":
            return f"{name}: {self.bounds(inner)}"
        raise Unreadable(f"unknown associated item constraint shape: {what!r}")

    # Generic parameters, bounds and where clauses

    def binder(self, generic_params) -> str:
        """The `for<..>` in front of a bound or a function pointer."""
        if not generic_params:
            return ""
        return "for<" + ", ".join(self.generic_param(param) for param in generic_params) + "> "

    def poly_trait(self, poly: dict) -> str:
        return self.binder(poly["generic_params"]) + self.path(poly["trait"])

    def bound(self, bound) -> str:
        what, inner = next(iter(bound.items()))
        if what == "trait_bound":
            modifier = {"none": "", "maybe": "?", "maybe_const": "~const "}.get(inner["modifier"])
            if modifier is None:
                raise Unreadable(f"unknown trait bound modifier: {inner['modifier']!r}")
            return modifier + self.poly_trait(inner)
        if what == "outlives":
            return inner
        if what == "use":
            return (
                "use<"
                + ", ".join(one_line(next(iter(arg.values())) if isinstance(arg, dict) else arg) for arg in inner)
                + ">"
            )
        raise Unreadable(f"unknown bound shape: {what!r}")

    def bounds(self, bounds) -> str:
        return " + ".join(self.bound(bound) for bound in bounds)

    def generic_param(self, param: dict) -> str | None:
        """One generic parameter, or None for the parameter an `impl Trait` argument stands for."""
        what, inner = next(iter(param["kind"].items()))
        if what == "lifetime":
            outlives = f": {' + '.join(inner['outlives'])}" if inner["outlives"] else ""
            return param["name"] + outlives
        if what == "type":
            if inner["is_synthetic"]:
                return None
            text = param["name"]
            if inner["bounds"]:
                text += ": " + self.bounds(inner["bounds"])
            if inner["default"] is not None:
                text += " = " + self.type(inner["default"])
            return text
        if what == "const":
            text = f"const {param['name']}: {self.type(inner['type'])}"
            if inner["default"] is not None:
                text += " = " + one_line(inner["default"])
            return text
        raise Unreadable(f"unknown generic parameter shape: {what!r}")

    def generics(self, generics: dict) -> tuple[str, str]:
        """The `<..>` after a name and the ` where ..` after a signature."""
        params = [text for text in map(self.generic_param, generics["params"]) if text is not None]
        predicates = [self.where_predicate(predicate) for predicate in generics["where_predicates"]]
        return (
            f"<{', '.join(params)}>" if params else "",
            f" where {', '.join(predicates)}" if predicates else "",
        )

    def where_predicate(self, predicate) -> str:
        what, inner = next(iter(predicate.items()))
        if what == "bound_predicate":
            return self.binder(inner["generic_params"]) + f"{self.type(inner['type'])}: {self.bounds(inner['bounds'])}"
        if what == "lifetime_predicate":
            return f"{inner['lifetime']}: {' + '.join(inner['outlives'])}"
        if what == "eq_predicate":
            return f"{self.type(inner['lhs'])} = {self.term(inner['rhs'])}"
        raise Unreadable(f"unknown where predicate shape: {what!r}")

    # Functions

    @staticmethod
    def header(header: dict) -> str:
        words = [
            word
            for word, flag in (("const", "is_const"), ("async", "is_async"), ("unsafe", "is_unsafe"))
            if header[flag]
        ]
        abi = header["abi"]
        if abi != "Rust":
            if isinstance(abi, dict):
                name, detail = next(iter(abi.items()))
                if isinstance(detail, dict) and detail.get("unwind"):
                    name += "-unwind"
                elif isinstance(detail, str):
                    name = detail
            else:
                name = abi
            words.append(f'extern "{name}"')
        return "".join(word + " " for word in words)

    def parameter(self, name: str, value) -> str:
        if name == "self":
            if value == {"generic": "Self"}:
                return "self"
            if (
                isinstance(value, dict)
                and "borrowed_ref" in value
                and value["borrowed_ref"]["type"] == {"generic": "Self"}
            ):
                reference = value["borrowed_ref"]
                lifetime = f"{reference['lifetime']} " if reference["lifetime"] else ""
                return f"&{lifetime}{'mut ' if reference['is_mutable'] else ''}self"
        return f"{name}: {self.type(value)}"

    def parameters_and_output(self, sig: dict) -> str:
        parts = [self.parameter(name, value) for name, value in sig["inputs"]]
        if sig["is_c_variadic"]:
            parts.append("...")
        output = f" -> {self.type(sig['output'])}" if sig["output"] is not None else ""
        return "(" + ", ".join(parts) + ")" + output

    def function(self, item: dict) -> str:
        inner = item["inner"]["function"]
        params, where = self.generics(inner["generics"])
        return self.header(inner["header"]) + "fn" + params + self.parameters_and_output(inner["sig"]) + where

    # Items

    @staticmethod
    def marks(item: dict) -> str:
        """The attributes of an item that a caller's code depends on."""
        marks = []
        for attribute in item.get("attrs", ()):
            if attribute == "non_exhaustive":
                marks.append("#[non_exhaustive]")
            elif attribute == "no_mangle":
                marks.append("#[no_mangle]")
            elif attribute == "macro_export":
                marks.append("#[macro_export]")
            elif isinstance(attribute, dict):
                what, inner = next(iter(attribute.items()))
                if what == "repr":
                    parts = [f"{key}={inner[key]}" for key in sorted(inner) if inner[key] is not None]
                    marks.append(f"#[repr({', '.join(parts)})]")
                elif what in ("export_name", "link_section"):
                    marks.append(f"#[{what} = {json.dumps(inner)}]")
                elif what == "target_feature":
                    marks.append(f"#[target_feature({', '.join(sorted(inner['enable']))})]")
                elif what == "other" and isinstance(inner, str) and inner.startswith(CFG_ATTRIBUTE):
                    marks.append(cfg_mark(inner))
        if item.get("deprecation") is not None:
            marks.append("#[deprecated]")
        return "".join(f"; {mark}" for mark in sorted(marks))

    def fields_shape(self, kind, of_variant: bool) -> str:
        """The shape of a struct's or a variant's fields; a variant's fields are never private."""
        if isinstance(kind, str):
            return "no fields"
        if "tuple" in kind:
            fields = kind["tuple"]
            text = f"{len(fields)} positional field" + ("" if len(fields) == 1 else "s")
            private = 0 if of_variant else sum(1 for field_id in fields if self.is_private_field(field_id))
            return text + (f", {private} private" if private else "")
        named = kind.get("plain") or kind.get("struct")
        private = not of_variant and (
            named["has_stripped_fields"] or any(self.is_private_field(field_id) for field_id in named["fields"])
        )
        return "named fields" + (", some private" if private else "")

    def is_private_field(self, field_id) -> bool:
        field = self.index.get(str(field_id)) if field_id is not None else None
        return field is None or field["visibility"] != "public"

    def block_header(self, block: dict) -> str:
        """The header of an inherent impl block that narrows what its members apply to, or nothing."""
        inner = block["inner"]["impl"]
        params, where = self.generics(inner["generics"])
        target = inner["for"]
        plain = "resolved_path" in target and not self.generic_args(target["resolved_path"].get("args"))
        if not params and not where and plain:
            return ""
        return f"; in impl{params} {self.type(target)}{where}"

    def implementation(self, block: dict) -> str:
        inner = block["inner"]["impl"]
        params, where = self.generics(inner["generics"])
        text = ("unsafe " if inner["is_unsafe"] else "") + f"impl{params} "
        text += ("!" if inner["is_negative"] else "") + self.path(inner["trait"])
        text += f" for {self.type(inner['for'])}{where}"
        members = []
        for member_id in inner["items"]:
            member = self.index.get(str(member_id))
            if member is None:
                continue
            if "assoc_type" in member["inner"] and member["inner"]["assoc_type"]["type"] is not None:
                members.append(f"type {member['name']} = {self.type(member['inner']['assoc_type']['type'])}")
            elif "assoc_const" in member["inner"]:
                members.append(f"const {member['name']}: {self.type(member['inner']['assoc_const']['type'])}")
        if members:
            text += " { " + "; ".join(sorted(members)) + " }"
        return text

    def auto_traits(self, blocks) -> str:
        parts = []
        for block in blocks:
            inner = block["inner"]["impl"]
            definition = self.definition(inner["trait"]["id"])
            name = AUTO_TRAITS.get(definition[1]) if definition else None
            if name is None:
                continue
            _params, where = self.generics(inner["generics"])
            parts.append((name, ("!" if inner["is_negative"] else "") + name + where))
        return ", ".join(text for _name, text in sorted(parts))

    def reexport(self, use: dict) -> str:
        if "extern_crate" in use["inner"]:
            # rustdoc gives the crate's own name as `rename` where the item binds it `as` another name.
            inner = use["inner"]["extern_crate"]
            crate = inner["rename"] or inner["name"]
            self.named.add(self.crate_name if crate == "self" else crate)
            return f"pub extern crate {crate}"
        inner = use["inner"]["use"]
        target = self.named_path(inner["id"], inner["source"]) if inner["id"] is not None else inner["source"]
        return f"pub use {target}" + ("::*" if inner["is_glob"] else "")

    def documents_trait(self, block: dict) -> bool:
        """Whether the trait an impl block implements is one the API documentation shows."""
        definition = self.definition(block["inner"]["impl"]["trait"]["id"])
        if definition is None:
            return False
        crate, path = definition
        if path in UNLISTED_TRAITS:
            return False
        return path in self.documented_paths or not self.in_workspace(crate)

    def in_workspace(self, crate: str) -> bool:
        return crate == self.crate_name or crate in self.workspace

    workspace: frozenset = frozenset()

    def signature(self, entry: Entry) -> str:
        """The signature of one entry; `self.named` then holds the crates it names."""
        self.named = set()
        item = entry.item
        if entry.kind == "implementation":
            return self.implementation(item) + self.marks(item)
        if entry.kind == "auto traits":
            return self.auto_traits(entry.blocks)
        if entry.kind == "re-export":
            return self.reexport(item) + self.marks(item)
        what = inner_kind(item)
        inner = item["inner"][what]
        if what == "module":
            text = "mod"
        elif what == "struct":
            params, where = self.generics(inner["generics"])
            text = f"struct{params}{where}; {self.fields_shape(inner['kind'], of_variant=False)}"
        elif what == "union":
            params, where = self.generics(inner["generics"])
            private = inner["has_stripped_fields"] or any(
                self.is_private_field(field_id) for field_id in inner["fields"]
            )
            text = f"union{params}{where}; named fields" + (", some private" if private else "")
        elif what == "enum":
            params, where = self.generics(inner["generics"])
            text = f"enum{params}{where}"
        elif what == "variant":
            text = self.fields_shape(inner["kind"], of_variant=True)
            if inner["discriminant"] is not None:
                text += f"; = {self.constant_expression(inner['discriminant'])}"
        elif what == "struct_field":
            text = self.type(inner)
        elif what == "trait":
            params, where = self.generics(inner["generics"])
            words = ("unsafe " if inner["is_unsafe"] else "") + ("auto " if inner["is_auto"] else "")
            supertraits = f": {self.bounds(inner['bounds'])}" if inner["bounds"] else ""
            text = f"{words}trait{params}{supertraits}{where}; " + ("dyn" if inner["is_dyn_compatible"] else "no dyn")
        elif what == "trait_alias":
            params, where = self.generics(inner["generics"])
            text = f"trait{params} = {self.bounds(inner['params'])}{where}"
        elif what == "function":
            text = self.function(item)
            if entry.in_trait:
                text += "; provided" if inner["has_body"] else "; required"
        elif what == "constant":
            text = f"const: {self.type(inner['type'])}"
        elif what == "static":
            words = ("unsafe " if inner["is_unsafe"] else "") + "static" + (" mut" if inner["is_mutable"] else "")
            text = f"{words}: {self.type(inner['type'])}"
        elif what == "assoc_const":
            text = f"const: {self.type(inner['type'])}"
            if entry.in_trait:
                text += "; provided" if inner["value"] is not None else "; required"
        elif what == "assoc_type":
            params, where = self.generics(inner["generics"])
            text = f"type{params}"
            if inner["bounds"]:
                text += ": " + self.bounds(inner["bounds"])
            if inner["type"] is not None:
                text += " = " + self.type(inner["type"])
            text += where
        elif what == "type_alias":
            params, where = self.generics(inner["generics"])
            text = f"type{params} = {self.type(inner['type'])}{where}"
        elif what == "macro":
            text = one_line(inner)
        elif what == "proc_macro":
            text = f"proc macro {inner['kind']}" + "".join(f" {helper}" for helper in inner["helpers"])
        else:
            raise Unreadable(f"unknown item kind {what!r} at {entry.path}")
        if entry.block is not None:
            text += self.block_header(entry.block)
        return one_line(text + self.marks(item))


def fail(code: int, *findings: str) -> None:
    for finding in findings:
        print(f"ERROR: [public-interface] {finding}", file=sys.stderr)
    sys.exit(code)


def cargo_metadata(root: Path, offline: bool) -> dict:
    command = ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked"]
    if offline:
        command.append("--offline")
    try:
        done = subprocess.run(command, cwd=root, check=False, capture_output=True, text=True)
    except FileNotFoundError:
        fail(2, "cargo is not on PATH")
    if done.returncode != 0:
        fail(2, f"cargo metadata failed in {root}: {one_line(done.stderr)}")
    return json.loads(done.stdout)


def workspace_packages(metadata: dict) -> list[dict]:
    """The packages of the workspace, by name, each with the crate name of its library target or None."""
    members = set(metadata["workspace_members"])
    packages = []
    for package in metadata["packages"]:
        if package["id"] not in members:
            continue
        libraries = [
            target["name"].replace("-", "_") for target in package["targets"] if LIBRARY_KINDS & set(target["kind"])
        ]
        packages.append(
            {"name": package["name"], "library": libraries[0] if libraries else None, "features": package["features"]}
        )
    return sorted(packages, key=lambda package: package["name"])


def document_libraries(root: Path, metadata: dict, libraries: list[str], offline: bool) -> dict[str, dict]:
    """rustdoc's JSON output of every library of the workspace, by crate name."""
    if not libraries:
        return {}
    command = ["cargo", "doc", "--workspace", "--lib", "--no-deps", "--locked", "--all-features"]
    if offline:
        command.append("--offline")
    with tempfile.TemporaryDirectory(prefix="public-interface-") as scratch:
        wrapper = Path(scratch) / "rustdoc-json"
        wrapper.write_text(RUSTDOC_WRAPPER, encoding="utf-8")
        wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        environment = dict(os.environ)
        environment["PUBLIC_INTERFACE_RUSTDOC"] = environment.get("RUSTDOC", "rustdoc")
        environment["RUSTDOC"] = str(wrapper)
        environment.pop("RUSTDOCFLAGS", None)
        environment["CARGO_ENCODED_RUSTDOCFLAGS"] = "\x1f".join(RUSTDOC_FLAGS)
        done = subprocess.run(command, cwd=root, env=environment, check=False, capture_output=True, text=True)
    if done.returncode != 0:
        sys.stderr.write(done.stderr)
        fail(2, f"`{' '.join(command)}` failed in {root}")
    documents = {}
    for library in libraries:
        file = Path(metadata["target_directory"]) / "doc" / f"{library}.json"
        try:
            with open(file, encoding="utf-8") as handle:
                documents[library] = json.load(handle)
        except (OSError, ValueError) as error:
            fail(2, f"{file}: {error}")
    return documents


def package_lines(package: dict) -> list[str]:
    """The lines of a record that state the package's library target and features."""
    lines = [f"library\t{package['library'] or 'none'}"]
    for name in sorted(package["features"]):
        lines.append(f"feature\t{name}\t{' '.join(sorted(package['features'][name])) or '-'}")
    return lines


def listings(packages: list[dict], documents: dict[str, dict]) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """The item lines of each package, and for each crate outside the workspace the items that name it."""
    walked = {}
    documented_paths: dict[str, str] = {}
    for library, document in documents.items():
        entries = Walker(document).run()
        walked[library] = entries
        local = {key: "::".join(value["path"]) for key, value in document["paths"].items() if value["crate_id"] == 0}
        for entry in entries:
            if entry.hidden or entry.kind in ("re-export", "implementation", "auto traits"):
                continue
            definition = local.get(str(entry.item["id"]))
            if definition is None:
                continue
            known = documented_paths.get(definition)
            if known is None or (entry.path.count("::"), entry.path) < (known.count("::"), known):
                documented_paths[definition] = entry.path
    workspace = frozenset(documents)
    items: dict[str, list[str]] = {}
    outside: dict[str, list[str]] = {}
    for package in packages:
        library = package["library"]
        if library is None:
            items[package["name"]] = []
            continue
        signatures = Signatures(documents[library], documented_paths)
        signatures.workspace = workspace
        rows = set()
        for entry in walked[library]:
            if entry.hidden:
                continue
            if entry.kind == "implementation" and not signatures.documents_trait(entry.item):
                continue
            try:
                signature = signatures.signature(entry)
            except (KeyError, TypeError, StopIteration) as error:
                raise Unreadable(
                    f"{entry.kind} {entry.path}: rustdoc's output is not in a shape this script reads ({error!r})"
                )
            rows.add((entry.path, entry.kind, signature))
            for crate in signatures.named - STANDARD_LIBRARY:
                if not signatures.in_workspace(crate):
                    outside.setdefault(crate, []).append(f"{package['name']}: {entry.kind} {entry.path}")
        items[package["name"]] = [f"{kind}\t{path}\t{signature}" for path, kind, signature in sorted(rows)]
    return items, outside


def parse(text: str) -> list[tuple[str, str, str]]:
    """The (kind, path, signature) of each line of a record; a `library` line has no path."""
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "library":
            rows.append(("library", "", "\t".join(fields[1:])))
        else:
            rows.append((fields[0], fields[1] if len(fields) > 1 else "", "\t".join(fields[2:])))
    return rows


def differences(package: str, recorded: str, current: str) -> list[str]:
    """One finding for each line that the record and the tree do not share."""
    before: dict[tuple[str, str], list[str]] = {}
    after: dict[tuple[str, str], list[str]] = {}
    for rows, into in ((parse(recorded), before), (parse(current), after)):
        for kind, path, signature in rows:
            into.setdefault((kind, path), []).append(signature)
    findings = []
    for key in sorted(set(before) | set(after), key=lambda key: (key[1], key[0])):
        name = " ".join(part for part in key if part)
        old, new = before.get(key, []), after.get(key, [])
        gone = sorted(signature for signature in set(old) if old.count(signature) > new.count(signature))
        come = sorted(signature for signature in set(new) if new.count(signature) > old.count(signature))
        if len(gone) == 1 and len(come) == 1:
            findings.append(f"{package}: changed {name}: `{gone[0]}` became `{come[0]}`")
            continue
        findings += [f"{package}: removed {name}: `{signature}`" for signature in gone]
        findings += [f"{package}: added {name}: `{signature}`" for signature in come]
    return findings


def outside_crate_findings(record: Path, shown: str, outside: dict[str, list[str]]) -> list[str]:
    """Findings on the crates outside the workspace that documented signatures name, against the list."""
    listed: dict[str, str] = {}
    file = record / OUTSIDE_CRATES
    if file.is_file():
        for number, line in enumerate(file.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip() or line.startswith("#"):
                continue
            crate, _, reason = line.partition("\t")
            if not reason.strip():
                return [f"{shown}/{OUTSIDE_CRATES} line {number}: a line is `<crate><tab><reason>`"]
            listed[crate.strip()] = reason
    findings = []
    for crate in sorted(set(outside) - set(listed)):
        for place in sorted(set(outside[crate])):
            findings.append(
                f"{place} names a type of the crate `{crate}`, which {shown}/{OUTSIDE_CRATES} does not list"
            )
    for crate in sorted(set(listed) - set(outside)):
        findings.append(
            f"{shown}/{OUTSIDE_CRATES} lists `{crate}`, which no documented signature names: remove its line"
        )
    return findings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Hold the documented public interface of the workspace's crates to its record."
    )
    parser.add_argument("--update", action="store_true", help="write the record again from the tree")
    parser.add_argument("--offline", action="store_true", help="hand --offline to cargo")
    parser.add_argument("--root", type=Path, default=ROOT, help="the workspace to read")
    parser.add_argument(
        "--record", type=Path, help="the record's directory (default: tests/public-interface under the root)"
    )
    arguments = parser.parse_args(argv)
    root = arguments.root.resolve()
    record = (arguments.record or root / RECORD).resolve()
    try:
        shown = record.relative_to(root).as_posix()
    except ValueError:
        shown = str(record)

    metadata = cargo_metadata(root, arguments.offline)
    packages = workspace_packages(metadata)
    libraries = [package["library"] for package in packages if package["library"]]
    documents = document_libraries(root, metadata, libraries, arguments.offline)
    formats = sorted({document.get("format_version") for document in documents.values()} - {TESTED_FORMAT})
    if formats:
        print(
            f"NOTE: [public-interface] rustdoc wrote JSON format {', '.join(map(str, formats))}; "
            f"this script reads format {TESTED_FORMAT}, the format of the toolchain .github/workflows/ci.yml pins",
            file=sys.stderr,
        )
    try:
        items, outside = listings(packages, documents)
    except Unreadable as error:
        fail(2, str(error))
    except (AttributeError, LookupError, StopIteration, TypeError, ValueError) as error:
        fail(2, f"rustdoc's output is not in a shape this script reads ({error!r})")
    current = {
        package["name"]: "\n".join([HEADER, *package_lines(package), *items[package["name"]]]) + "\n"
        for package in packages
    }

    findings = []
    stale = sorted(file for file in record.glob("*.txt") if file.stem not in current) if record.is_dir() else []
    if arguments.update:
        record.mkdir(parents=True, exist_ok=True)
        written = 0
        for name, text in current.items():
            file = record / f"{name}.txt"
            if not file.is_file() or file.read_text(encoding="utf-8") != text:
                file.write_text(text, encoding="utf-8")
                written += 1
        for file in stale:
            file.unlink()
        print(f"check_public_interface: {len(current)} listings under {shown}, {written} written, {len(stale)} removed")
    else:
        for name, text in current.items():
            file = record / f"{name}.txt"
            if not file.is_file():
                findings.append(f"{name}: {shown}/{name}.txt does not exist, so nothing holds the package's interface")
                continue
            recorded = file.read_text(encoding="utf-8")
            if recorded != text:
                found = differences(name, recorded, text)
                findings += found or [f"{name}: {shown}/{name}.txt is not the file this script writes"]
        findings += [f"{file.stem}: {shown}/{file.name} is the record of no package of the workspace" for file in stale]
    listed = outside_crate_findings(record, shown, outside)
    if findings or listed:
        for finding in findings + listed:
            print(f"ERROR: [public-interface] {finding}", file=sys.stderr)
        closing = []
        if findings:
            closing.append(
                f"a deliberate change writes the record again with `{UPDATE_COMMAND}` and commits the diff of {shown}"
            )
        if listed:
            closing.append(f"{shown}/{OUTSIDE_CRATES} is edited by hand, one line a crate with its reason")
        print(f"ERROR: [public-interface] {'; '.join(closing)}", file=sys.stderr)
        return 1
    if not arguments.update:
        print("check_public_interface: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
