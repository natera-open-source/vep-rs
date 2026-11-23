"""Tests for check_public_interface.py on a made-up workspace the tests write: a library `shapes`, a
package `tool` with a binary alone, and beside the workspace a crate `outsider` that `shapes` depends on by
path. Each test copies that tree with its record, plants one change and reads the findings the check
prints. The last tests read the repository's own record and workflow. They need cargo
on the toolchain `.github/workflows/ci.yml` pins and no network."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "check_public_interface.py"
ROOT = HERE.parent
PREFIX = "ERROR: [public-interface] "
RECORD = "tests/public-interface"
UPDATE = "python3 scripts/check_public_interface.py --update"
CLOSING = f"a deliberate change writes the record again with `{UPDATE}` and commits the diff of {RECORD}"
LIST_CLOSING = f"{RECORD}/outside-crates.tsv is edited by hand, one line a crate with its reason"

WORKSPACE_MANIFEST = """\
[workspace]
resolver = "2"
members = ["shapes", "tool"]
"""

SHAPES_MANIFEST = """\
[package]
name = "shapes"
version = "0.0.0"
edition = "2021"

[features]
default = ["fast"]
fast = []

[dependencies]
outsider = { path = "../../outsider" }
"""

SHAPES_LIBRARY = """\
//! A made-up library.

use std::path::Path;

/// Doubles a value.
pub fn measure(value: u32) -> u32 {
    helper(value)
}

fn helper(value: u32) -> u32 {
    value * 2
}

/// Whether a path has a name.
pub fn named(path: &Path) -> bool {
    path.file_name().is_some()
}

/// How many sides a shape has.
pub fn count(shape: &Shape) -> u8 {
    shape.sides
}

/// Listed in the record.
pub fn listed() {}

/// A shape.
#[derive(Clone, Debug)]
pub struct Shape {
    /// How many sides it has.
    pub sides: u8,
}

/// A colour.
pub enum Colour {
    /// Red.
    Red,
    /// Green.
    Green,
}

#[doc(hidden)]
pub mod workshop {
    /// A tool, offered at the crate root.
    pub struct Tool {
        /// Its weight.
        pub weight: u8,
    }

    /// Stays out of the documentation.
    pub fn sharpen() {}
}

pub use workshop::Tool;
"""

TOOL_MANIFEST = """\
[package]
name = "tool"
version = "0.0.0"
edition = "2021"
"""

OUTSIDER_MANIFEST = """\
[package]
name = "outsider"
version = "0.0.0"
edition = "2021"
"""

OUTSIDER_LIBRARY = """\
//! A crate outside the workspace.

/// A thing of another crate.
pub struct Thing;
"""


def environment() -> dict[str, str]:
    """The environment of a run: nothing that moves cargo's output directory or rustdoc's flags."""
    dropped = ("CARGO_TARGET_DIR", "CARGO_BUILD_TARGET", "RUSTDOC", "RUSTDOCFLAGS", "CARGO_ENCODED_RUSTDOCFLAGS")
    return {name: value for name, value in os.environ.items() if name not in dropped}


def check(tree: Path, *arguments: str) -> subprocess.CompletedProcess:
    command = [sys.executable, str(SCRIPT), "--root", str(tree / "workspace"), "--offline", *arguments]
    return subprocess.run(command, check=False, capture_output=True, text=True, env=environment())


def findings(done: subprocess.CompletedProcess) -> list[str]:
    assert all(line.startswith(PREFIX) for line in done.stderr.splitlines()), done.stderr
    return [line[len(PREFIX) :] for line in done.stderr.splitlines()]


@pytest.fixture(scope="module")
def recorded(tmp_path_factory) -> Path:
    """The made-up tree with its record written, which every test copies."""
    tree = tmp_path_factory.mktemp("recorded")
    files = {
        "workspace/Cargo.toml": WORKSPACE_MANIFEST,
        "workspace/shapes/Cargo.toml": SHAPES_MANIFEST,
        "workspace/shapes/src/lib.rs": SHAPES_LIBRARY,
        "workspace/tool/Cargo.toml": TOOL_MANIFEST,
        "workspace/tool/src/main.rs": "fn main() {}\n",
        "outsider/Cargo.toml": OUTSIDER_MANIFEST,
        "outsider/src/lib.rs": OUTSIDER_LIBRARY,
    }
    for name, text in files.items():
        file = tree / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(text, encoding="utf-8")
    subprocess.run(
        ["cargo", "generate-lockfile", "--offline"],
        cwd=tree / "workspace",
        check=True,
        capture_output=True,
        env=environment(),
    )
    done = check(tree, "--update")
    assert done.returncode == 0, done.stderr
    assert done.stdout == f"check_public_interface: 2 listings under {RECORD}, 2 written, 0 removed\n"
    return tree


@pytest.fixture
def tree(recorded: Path, tmp_path: Path) -> Path:
    """A copy of the recorded tree, without its build output."""
    copy = tmp_path / "tree"
    shutil.copytree(recorded, copy, ignore=shutil.ignore_patterns("target"))
    return copy


def edit(tree: Path, name: str, old: str, new: str) -> None:
    file = tree / "workspace" / name
    text = file.read_text(encoding="utf-8")
    assert text.count(old) == 1, f"{name} holds {old!r} {text.count(old)} times"
    file.write_text(text.replace(old, new), encoding="utf-8")


def refused(tree: Path, *expected: str, closing: str = CLOSING) -> None:
    done = check(tree)
    assert done.returncode == 1, done.stderr
    assert done.stdout == ""
    assert findings(done) == [*expected, closing]


def passes(tree: Path) -> None:
    done = check(tree)
    assert (done.returncode, done.stdout, done.stderr) == (0, "check_public_interface: OK\n", "")


def test_the_recorded_tree_passes_and_its_record_states_each_package(recorded):
    passes(recorded)
    record = recorded / "workspace" / RECORD
    assert sorted(file.name for file in record.iterdir()) == ["shapes.txt", "tool.txt"]
    shapes = (record / "shapes.txt").read_text(encoding="utf-8").splitlines()
    assert shapes[0].startswith("# Written by scripts/check_public_interface.py --update")
    assert shapes[1:4] == ["library\tshapes", "feature\tdefault\tfast", "feature\tfast\t-"]
    assert "function\tshapes::measure\tfn(value: u32) -> u32" in shapes
    assert "function\tshapes::named\tfn(path: &std::path::Path) -> bool" in shapes
    assert "struct\tshapes::Shape\tstruct; named fields" in shapes
    assert "field\tshapes::Shape::sides\tu8" in shapes
    assert "implementation\tshapes::Shape\timpl core::clone::Clone for shapes::Shape" in shapes
    assert "auto traits\tshapes::Shape\tRefUnwindSafe, Send, Sync, Unpin, UnwindSafe" in shapes
    assert "variant\tshapes::Colour::Green\tno fields" in shapes
    assert shapes[4:] == sorted(
        shapes[4:], key=lambda line: (line.split("\t")[1], line.split("\t")[0], line.split("\t")[2])
    )
    assert (record / "tool.txt").read_text(encoding="utf-8").splitlines()[1:] == ["library\tnone"]


def test_an_item_offered_out_of_a_hidden_module_is_listed_whole_and_the_module_is_not(recorded):
    shapes = (recorded / "workspace" / RECORD / "shapes.txt").read_text(encoding="utf-8").splitlines()
    assert "struct\tshapes::Tool\tstruct; named fields" in shapes
    assert "field\tshapes::Tool::weight\tu8" in shapes
    assert not any("workshop" in line for line in shapes)


def test_an_added_public_function_is_refused_by_name(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.",
        "/// Added.\npub fn added() -> u8 {\n    1\n}\n\n/// Listed in the record.",
    )
    refused(tree, "shapes: added function shapes::added: `fn() -> u8`")


def test_a_removed_public_function_is_refused_by_name(tree):
    edit(tree, "shapes/src/lib.rs", "/// Listed in the record.\npub fn listed() {}\n", "")
    refused(tree, "shapes: removed function shapes::listed: `fn()`")


def test_a_changed_parameter_type_is_refused_with_both_signatures(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "pub fn measure(value: u32) -> u32 {\n    helper(value)",
        "pub fn measure(value: u64) -> u32 {\n    helper(value as u32)",
    )
    refused(tree, "shapes: changed function shapes::measure: `fn(value: u32) -> u32` became `fn(value: u64) -> u32`")


def test_a_changed_return_type_is_refused_with_both_signatures(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "pub fn measure(value: u32) -> u32 {\n    helper(value)",
        "pub fn measure(value: u32) -> u64 {\n    helper(value).into()",
    )
    refused(tree, "shapes: changed function shapes::measure: `fn(value: u32) -> u32` became `fn(value: u32) -> u64`")


def test_a_new_variant_of_an_enum_is_refused_by_name(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "    /// Green.\n    Green,\n",
        "    /// Green.\n    Green,\n    /// Blue.\n    Blue,\n",
    )
    refused(tree, "shapes: added variant shapes::Colour::Blue: `no fields`")


def test_a_new_public_field_is_refused_by_name(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "    pub sides: u8,\n",
        "    pub sides: u8,\n    /// How many corners it has.\n    pub corners: u8,\n",
    )
    refused(tree, "shapes: added field shapes::Shape::corners: `u8`")


def test_a_new_private_field_changes_the_line_of_its_struct(tree):
    edit(tree, "shapes/src/lib.rs", "    pub sides: u8,\n", "    pub sides: u8,\n    corners: u8,\n")
    refused(
        tree, "shapes: changed struct shapes::Shape: `struct; named fields` became `struct; named fields, some private`"
    )


def test_an_item_moved_under_doc_hidden_reads_as_removed(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.\npub fn listed",
        "/// Listed in the record.\n#[doc(hidden)]\npub fn listed",
    )
    refused(tree, "shapes: removed function shapes::listed: `fn()`")


def test_an_item_moved_behind_a_feature_changes_its_line(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.\npub fn listed",
        '/// Listed in the record.\n#[cfg(all(feature = "fast", not(windows)))]\npub fn listed',
    )
    refused(
        tree,
        'shapes: changed function shapes::listed: `fn()` became `fn(); #[cfg(all(feature = "fast", not(windows)))]`',
    )


def test_an_implementation_and_a_re_export_behind_a_feature_carry_the_mark(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.",
        '#[cfg(feature = "fast")]\nimpl Default for Colour {\n    fn default() -> Self {\n        Colour::Red\n    }\n}\n\n'
        '#[cfg(feature = "fast")]\npub use std::path::Path as Route;\n\n/// Listed in the record.',
    )
    refused(
        tree,
        "shapes: added implementation shapes::Colour: "
        '`impl core::default::Default for shapes::Colour; #[cfg(feature = "fast")]`',
        'shapes: added re-export shapes::Route: `pub use std::path::Path; #[cfg(feature = "fast")]`',
    )


def test_a_trait_implemented_through_a_type_alias_is_listed_at_the_trait(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.",
        "/// Can be drawn.\npub trait Draw {}\n\n/// Another name for a shape.\npub type Sketch = Shape;\n\n"
        "impl Draw for Sketch {}\n\n/// Listed in the record.",
    )
    refused(
        tree,
        "shapes: added implementation shapes::Draw: `impl shapes::Draw for shapes::Sketch`",
        "shapes: added trait shapes::Draw: `trait; dyn`",
        "shapes: added type alias shapes::Sketch: `type = shapes::Shape`",
    )


def test_a_static_of_an_extern_block_is_listed_as_unsafe(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.",
        'extern "C" {\n    /// Defined by another library.\n    pub static EDGE: u8;\n}\n\n/// Listed in the record.',
    )
    refused(tree, "shapes: added static shapes::EDGE: `unsafe static: u8`")


def test_a_dropped_derive_and_a_lost_auto_trait_are_refused(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "#[derive(Clone, Debug)]\npub struct Shape {\n",
        "#[derive(Debug)]\npub struct Shape {\n    _owner: std::rc::Rc<u8>,\n",
    )
    refused(
        tree,
        "shapes: changed auto traits shapes::Shape: `RefUnwindSafe, Send, Sync, Unpin, UnwindSafe` became "
        "`RefUnwindSafe, !Send, !Sync, Unpin, UnwindSafe`",
        "shapes: removed implementation shapes::Shape: `impl core::clone::Clone for shapes::Shape`",
        "shapes: changed struct shapes::Shape: `struct; named fields` became `struct; named fields, some private`",
    )


def test_a_new_cargo_feature_is_refused_by_name(tree):
    edit(tree, "shapes/Cargo.toml", "fast = []\n", "fast = []\nextra = []\n")
    refused(tree, "shapes: added feature extra: `-`")


def test_a_feature_dropped_from_the_default_set_is_refused(tree):
    edit(tree, "shapes/Cargo.toml", 'default = ["fast"]', "default = []")
    refused(tree, "shapes: changed feature default: `fast` became `-`")


def test_a_library_target_added_to_a_package_that_had_none_is_refused(tree):
    (tree / "workspace" / "tool" / "src" / "lib.rs").write_text(
        "//! A library.\n\n/// Offered.\npub fn offered() {}\n", encoding="utf-8"
    )
    refused(
        tree,
        "tool: changed library: `none` became `tool`",
        "tool: added module tool: `mod`",
        "tool: added function tool::offered: `fn()`",
    )


def test_a_signature_that_names_a_crate_outside_the_list_is_refused_until_the_list_names_it(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.",
        "/// Borrows a thing.\npub fn borrowed(_thing: &outsider::Thing) {}\n\n/// Listed in the record.",
    )
    added = "shapes: added function shapes::borrowed: `fn(_thing: &outsider::Thing)`"
    unlisted = f"shapes: function shapes::borrowed names a type of the crate `outsider`, which {RECORD}/outside-crates.tsv does not list"
    refused(tree, added, unlisted, closing=f"{CLOSING}; {LIST_CLOSING}")
    # The update command writes the listing and leaves the list to a hand.
    done = check(tree, "--update")
    assert done.returncode == 1
    assert findings(done) == [unlisted, LIST_CLOSING]
    refused(tree, unlisted, closing=LIST_CLOSING)
    (tree / "workspace" / RECORD / "outside-crates.tsv").write_text(
        "# crate\treason\noutsider\tits thing is borrowed\n", encoding="utf-8"
    )
    passes(tree)


def test_a_public_extern_crate_item_names_its_crate_and_is_held_to_the_list(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "/// Listed in the record.",
        "pub extern crate outsider as other;\npub extern crate self as me;\n\n/// Listed in the record.",
    )
    refused(
        tree,
        "shapes: added re-export shapes::me: `pub extern crate self`",
        "shapes: added re-export shapes::other: `pub extern crate outsider`",
        f"shapes: re-export shapes::other names a type of the crate `outsider`, which {RECORD}/outside-crates.tsv does not list",
        closing=f"{CLOSING}; {LIST_CLOSING}",
    )


def test_a_listed_crate_that_no_signature_names_is_refused(tree):
    (tree / "workspace" / RECORD / "outside-crates.tsv").write_text(
        "outsider\tits thing is borrowed\n", encoding="utf-8"
    )
    refused(
        tree,
        f"{RECORD}/outside-crates.tsv lists `outsider`, which no documented signature names: remove its line",
        closing=LIST_CLOSING,
    )


def test_a_changed_doc_comment_alone_passes(tree):
    edit(tree, "shapes/src/lib.rs", "/// Doubles a value.", "/// Returns twice the value it is given.")
    passes(tree)


def test_a_changed_private_item_alone_passes(tree):
    edit(
        tree,
        "shapes/src/lib.rs",
        "fn helper(value: u32) -> u32 {\n    value * 2",
        "fn helper(value: u32) -> u32 {\n    value + value",
    )
    edit(
        tree,
        "shapes/src/lib.rs",
        "    /// Stays out of the documentation.\n    pub fn sharpen() {}",
        "    /// Stays out of the documentation.\n    pub fn sharpen(_edge: u8) {}",
    )
    passes(tree)


def test_a_type_written_by_another_path_to_the_same_item_passes(tree):
    edit(tree, "shapes/src/lib.rs", "use std::path::Path;\n", "")
    edit(tree, "shapes/src/lib.rs", "pub fn named(path: &Path)", "pub fn named(path: &std::path::Path)")
    # rustdoc gives a path inside the crate as the source writes it, so the listing alone makes these one.
    edit(tree, "shapes/src/lib.rs", "pub fn count(shape: &Shape)", "pub fn count(shape: &crate::Shape)")
    passes(tree)


def test_the_update_command_writes_the_record_a_changed_tree_then_passes_on(tree):
    edit(tree, "shapes/src/lib.rs", "/// Listed in the record.\npub fn listed() {}\n", "")
    (tree / "workspace" / RECORD / "gone.txt").write_text("library\tnone\n", encoding="utf-8")
    refused(
        tree,
        "shapes: removed function shapes::listed: `fn()`",
        f"gone: {RECORD}/gone.txt is the record of no package of the workspace",
    )
    done = check(tree, "--update")
    assert (done.returncode, done.stderr) == (0, "")
    assert done.stdout == f"check_public_interface: 2 listings under {RECORD}, 1 written, 1 removed\n"
    passes(tree)


def test_a_package_without_a_record_is_refused(tree):
    (tree / "workspace" / RECORD / "tool.txt").unlink()
    refused(tree, f"tool: {RECORD}/tool.txt does not exist, so nothing holds the package's interface")


def test_a_library_that_does_not_compile_ends_the_run_with_exit_2_and_no_finding(tree):
    edit(tree, "shapes/src/lib.rs", "pub fn listed() {}", "pub fn listed() -> Missing {}")
    done = check(tree)
    assert done.returncode == 2
    assert done.stdout == ""
    last = done.stderr.splitlines()[-1]
    assert last.startswith(
        f"{PREFIX}`cargo doc --workspace --lib --no-deps --locked --all-features --offline` failed in "
    )
    assert "added" not in last and "removed" not in last and UPDATE not in done.stderr


def unreadable(tree: Path, tmp_path: Path, rewrite: str, expected: str) -> None:
    """The check ends with exit 2 and the one line `expected` where a stand-in for rustdoc documents the
    library and then applies the sed expression `rewrite` to the JSON output."""
    stand_in = tmp_path / "rustdoc"
    stand_in.write_text(
        '#!/bin/sh\nrustdoc "$@" || exit $?\n'
        f"sed '{rewrite}' target/doc/shapes.json > target/doc/shapes.tmp\n"
        "mv target/doc/shapes.tmp target/doc/shapes.json\n",
        encoding="utf-8",
    )
    stand_in.chmod(0o755)
    command = [sys.executable, str(SCRIPT), "--root", str(tree / "workspace"), "--offline"]
    done = subprocess.run(
        command, check=False, capture_output=True, text=True, env={**environment(), "RUSTDOC": str(stand_in)}
    )
    assert done.returncode == 2, done.stderr
    assert done.stdout == ""
    assert findings(done) == [expected]


def test_rustdoc_output_in_a_shape_the_script_does_not_read_ends_the_run_with_exit_2_and_no_finding(tree, tmp_path):
    unreadable(
        tree,
        tmp_path,
        's/"impls":/"blocks":/g',
        "rustdoc's output is not in a shape this script reads (KeyError('impls'))",
    )


def test_an_item_kind_the_script_does_not_know_ends_the_run_with_exit_2_and_is_named(tree, tmp_path):
    # rustdoc writes a kind that has no content, such as an extern type, as its name alone.
    unreadable(
        tree,
        tmp_path,
        's/"inner":{"module":{"is_crate":false,[^}]*}}/"inner":"extern_type"/',
        "unknown item kind 'extern_type' at shapes::workshop",
    )


def test_two_runs_write_the_same_bytes(tree, tmp_path):
    written = []
    for name in ("first", "second"):
        record = tmp_path / name
        done = check(tree, "--update", "--record", str(record))
        assert done.returncode == 0, done.stderr
        written.append({file.name: file.read_bytes() for file in sorted(record.iterdir())})
        shutil.rmtree(tree / "workspace" / "target")
    assert written[0] == written[1]
    assert written[0]["shapes.txt"] == (tree / "workspace" / RECORD / "shapes.txt").read_bytes()


def workspace_packages() -> dict[str, dict]:
    command = ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked", "--offline"]
    metadata = json.loads(subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True).stdout)
    members = set(metadata["workspace_members"])
    return {package["name"]: package for package in metadata["packages"] if package["id"] in members}


def test_the_repository_s_record_has_one_listing_for_each_package_and_states_its_features():
    packages = workspace_packages()
    record = ROOT / RECORD
    assert sorted(file.stem for file in record.glob("*.txt")) == sorted(packages)
    for name, package in packages.items():
        lines = (record / f"{name}.txt").read_text(encoding="utf-8").splitlines()
        stated = {line.split("\t")[1]: line.split("\t")[2] for line in lines if line.startswith("feature\t")}
        assert stated == {feature: " ".join(sorted(enables)) or "-" for feature, enables in package["features"].items()}
        has_library = any("lib" in target["kind"] for target in package["targets"])
        assert (lines[1] != "library\tnone") == has_library, name


def test_the_converter_s_record_says_it_has_no_library():
    lines = (ROOT / RECORD / "vep-cache-converter.txt").read_text(encoding="utf-8").splitlines()
    assert lines[1] == "library\tnone"
    assert not [line for line in lines[2:] if not line.startswith("feature\t")]


def test_every_listed_outside_crate_carries_a_reason():
    lines = (ROOT / RECORD / "outside-crates.tsv").read_text(encoding="utf-8").splitlines()
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
    assert rows and all(len(row) == 2 and row[0].strip() and len(row[1].split()) >= 4 for row in rows)
    assert len({row[0] for row in rows}) == len(rows)


def test_the_workflow_runs_the_check_offline_after_the_documentation_build():
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    step = "run: python3 scripts/check_public_interface.py --offline"
    assert workflow.count(step) == 1
    condition = workflow[: workflow.index(step)].rstrip().splitlines()[-1]
    assert "matrix.name == 'Module Docs'" in condition
    assert workflow.index("run: ${{ matrix.command }}") < workflow.index(step)
