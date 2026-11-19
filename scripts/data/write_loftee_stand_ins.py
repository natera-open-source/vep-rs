#!/usr/bin/env python3
"""Write stand-ins for the files of LoFTEE's maxEntScan directory in a checkout made without it.

usage: write_loftee_stand_ins.py CHECKOUT [--dry-run]

LoFTEE's Perl module loads scripts from its `maxEntScan/` directory and calls subroutines those
scripts define. In a checkout made without that directory, this script writes there the fewest
files with which the module loads and runs with its splice predictions off. It holds no MaxEntScan
content and writes none: everything it writes comes from the call sites in the checkout's own
`*.pm` and `*.pl` files, which it reads outside that directory, line by line, without their
comments and POD:

- a `require` of a path under the directory names a file to write;
- a call written `&name`, where no file read defines `name` and no `use` line imports it by name,
  names a subroutine to define. A call written without the ampersand is not read, since the same
  spelling calls a function that a `use` line imports without naming it;
- a quoted path under the directory on the line of such a call names a file that subroutine would
  read, and no file is written for it.

One stand-in is written per required path. The stand-in for the first required path, in the order
of the files' names and then of their lines, defines every subroutine, and each other stand-in
defines none: Perl finds a subroutine by its name, whichever file defines it. A subroutine whose
every call sits in a constructor (`sub new`) returns an empty list, which the array or the hash of
its caller holds. Any other subroutine ends the process with exit status 1 and one line on standard
error, `ERROR: [loftee_stand_in] <name> was called: ...`, which says that MaxEntScan scoring is not
available in the stand-in; it never returns a score. The line goes to a duplicate of standard error
taken when the stand-in is loaded, and the process exits, because Ensembl VEP traps a plugin's
`die` and reopens standard error in a forked child. Each stand-in opens with a comment that says it
is a stand-in, which file it replaces and that it holds no MaxEntScan content.

A file already in the directory is the script's own when its name is one the script would write
and its size is the size the script would write, and it is left as it is. The script opens no file
that is already in the directory.

`--dry-run` prints each file with the subroutines it would define, and writes nothing.

Refused, with nothing written: a CHECKOUT that is not a directory; a checkout where no Perl file
requires a path under the directory; a quoted path under the directory on a line that calls no
subroutine the stand-ins define, a file the module would read itself; an entry of the directory
that is not the script's own.

Exit 0 with one line naming the files written and the files left; exit 1 with one
`ERROR: [write_loftee_stand_ins] ...` line for a refusal.
"""

from __future__ import annotations

import argparse
import os
import re
import stat
import sys
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import NamedTuple

COMPONENT = "write_loftee_stand_ins"
STAND_IN_COMPONENT = "loftee_stand_in"
DIRECTORY = "maxEntScan"
CONSTRUCTOR = "new"
SOURCE_SUFFIXES = (".pm", ".pl")
SKIPPED_DIRECTORIES = (DIRECTORY, ".git")

POD_START_RE = re.compile(r"^=[A-Za-z]")
POD_END_RE = re.compile(r"^=cut\b")
SOURCE_END_RE = re.compile(r"^__(END|DATA)__\s*$")
REQUIRE_RE = re.compile(r"""\brequire\s*\(?\s*(["'])""" + DIRECTORY + r"""/([\w.\-/]+)\1""")
PATH_RE = re.compile(r"""(["'])(?:[^"']*/)?(""" + DIRECTORY + r"""(?:/[^"']*)?)\1""")
USE_RE = re.compile(r"^\s*use\s+[A-Za-z_][\w:]*\s+([^;]*);")
NAME_RE = re.compile(r"[A-Za-z_]\w*")
# A declaration, a brace, or an ampersand call; a backslash before the ampersand takes a reference.
TOKEN_RE = re.compile(r"\bsub\s+([A-Za-z_]\w*)|([{}])|(?<![&\w$])(\\?)&([A-Za-z_]\w*)(?![\w:])")

STOP_REASON = (
    "MaxEntScan scoring is not available in the stand-in for {path}, which holds no MaxEntScan content. "
    "Run LoFTEE with its splice predictions off, or with its own " + DIRECTORY + " directory."
)


class Site(NamedTuple):
    file: str
    line: int

    def __str__(self) -> str:
        return f"{self.file} line {self.line}"


@dataclass
class Scan:
    """The call sites of a checkout: the first site that requires each path under the directory;
    per `&name`, each site with whether it sits in a constructor; the names the files define and
    import; and each quoted path under the directory with the `&` names called on its line."""

    requires: dict[str, Site] = field(default_factory=dict)
    calls: dict[str, list[tuple[Site, bool]]] = field(default_factory=dict)
    defined: set[str] = field(default_factory=set)
    imported: set[str] = field(default_factory=set)
    paths: list[tuple[Site, str, list[str]]] = field(default_factory=list)


@dataclass
class StandIn:
    """One file to write: its path in the checkout, the site that requires it, its subroutines as
    (name, whether it returns an empty list, its call sites) and its bytes."""

    path: str
    site: Site
    subroutines: list[tuple[str, bool, list[Site]]]
    text: bytes


def code_of(line: str) -> tuple[str, str]:
    """(the line without its comment; the same with the text inside its quotes blanked)."""
    kept: list[str] = []
    bare: list[str] = []
    quote = ""
    escaped = False
    for index, char in enumerate(line):
        if quote:
            kept.append(char)
            closing = char == quote and not escaped
            escaped = char == "\\" and not escaped
            bare.append(char if closing else " ")
            if closing:
                quote = ""
            continue
        # `$#name` is the last index of an array, not a comment.
        if char == "#" and (index == 0 or line[index - 1] != "$"):
            break
        if char in "\"'":
            quote = char
        kept.append(char)
        bare.append(char)
    return "".join(kept), "".join(bare)


def code_lines(text: str) -> list[tuple[int, str, str]]:
    """(line number, code, code with quoted text blanked) of every line outside POD, up to an
    `__END__` or `__DATA__` line."""
    lines: list[tuple[int, str, str]] = []
    in_pod = False
    for lineno, line in enumerate(text.splitlines(), start=1):
        if in_pod:
            in_pod = not POD_END_RE.match(line)
        elif POD_START_RE.match(line):
            in_pod = True
        elif SOURCE_END_RE.match(line):
            break
        else:
            lines.append((lineno, *code_of(line)))
    return lines


def read_sources(checkout: Path) -> list[tuple[str, str]]:
    """(path in the checkout, text) of every Perl file outside the directory, in the order of the paths."""
    sources: list[tuple[str, str]] = []
    for root, directories, names in os.walk(checkout):
        if Path(root) == checkout:
            directories[:] = [name for name in directories if name not in SKIPPED_DIRECTORIES]
        for name in names:
            if name.endswith(SOURCE_SUFFIXES):
                path = Path(root) / name
                with open(path, encoding="utf-8", errors="replace") as handle:
                    sources.append((path.relative_to(checkout).as_posix(), handle.read()))
    return sorted(sources)


def scan_sources(sources: list[tuple[str, str]]) -> Scan:
    scan = Scan()
    for file, text in sources:
        enclosing, depth = "", 0
        for lineno, code, bare in code_lines(text):
            site = Site(file, lineno)
            use = USE_RE.match(code)
            if use:
                scan.imported.update(NAME_RE.findall(use.group(1)))
            called: list[str] = []
            for declared, brace, reference, name in TOKEN_RE.findall(bare):
                if declared:
                    scan.defined.add(declared)
                    enclosing, depth = declared, 0
                elif brace == "{":
                    depth += 1
                elif brace == "}":
                    depth -= 1
                    if depth <= 0:
                        enclosing, depth = "", 0
                else:
                    in_constructor = enclosing == CONSTRUCTOR and depth > 0 and not reference
                    scan.calls.setdefault(name, []).append((site, in_constructor))
                    called.append(name)
            required: set[int] = set()
            for match in REQUIRE_RE.finditer(code):
                if ".." not in PurePosixPath(match.group(2)).parts:
                    scan.requires.setdefault(match.group(2), site)
                    required.add(match.start(1))
            for match in PATH_RE.finditer(code):
                if match.start(1) not in required:
                    scan.paths.append((site, match.group(2), called))
    return scan


def stand_in_text(path: str, site: Site, subroutines: list[tuple[str, bool, list[Site]]], first: str) -> bytes:
    lines = [
        f"# This file is a stand-in written by {Path(__file__).name}.",
        f"# It replaces {path}, which {site} requires.",
        "# It holds no MaxEntScan content: no model, no table and no line of MaxEntScan's code.",
    ]
    if path != first:
        lines.append(f"# It defines no subroutine: the stand-in for {first} defines them all.")
    elif subroutines:
        lines += [
            "# A subroutine that the module's constructor alone calls returns an empty list. Any other",
            "# ends the process, since it has no score to return.",
        ]
    if any(not empty for _name, empty, _sites in subroutines):
        reason = STOP_REASON.format(path=path)
        lines += [
            "",
            "# Ensembl VEP traps a plugin's die and reopens standard error in a forked child, so the",
            "# message goes to a duplicate of standard error taken at load and the process exits.",
            "my $stand_in_stderr;",
            "open($stand_in_stderr, '>&', \\*STDERR) or $stand_in_stderr = \\*STDERR;",
            "my $stand_in_stop = sub {",
            "    my ($name) = @_;",
            f'    syswrite($stand_in_stderr, "ERROR: [{STAND_IN_COMPONENT}] $name was called: {reason}\\n");',
            "    exit 1;",
            "};",
        ]
    if subroutines:
        lines.append("")
    for name, empty, _sites in subroutines:
        body = "return ();" if empty else f"$stand_in_stop->('{name}');"
        lines.append(f"sub {name} {{ {body} }}")
    lines += ["", "1;"]
    return "".join(f"{line}\n" for line in lines).encode("utf-8")


def plan(scan: Scan) -> tuple[list[StandIn], list[tuple[Site, str]], str]:
    """(the stand-ins to write, the first defining every subroutine; the quoted paths no stand-in
    reads, each with its site; the reason nothing can be written, empty when something can)."""
    if not scan.requires:
        return [], [], f"no Perl file requires a path under {DIRECTORY}/: not a LoFTEE checkout"
    missing = {
        name: sites
        for name, sites in scan.calls.items()
        if name not in scan.defined and name not in scan.imported
    }
    for site, path, called in scan.paths:
        if not any(name in missing for name in called):
            reason = (
                f"{site} reads {path} outside a call to a subroutine the stand-ins define; "
                "a stand-in cannot supply that file"
            )
            return [], [], reason
    subroutines = [
        (name, all(in_constructor for _site, in_constructor in sites), [site for site, _in in sites])
        for name, sites in missing.items()
    ]
    ordered = sorted(scan.requires.items(), key=lambda item: item[1])
    first = f"{DIRECTORY}/{ordered[0][0]}"
    stand_ins = []
    for relative, site in ordered:
        path = f"{DIRECTORY}/{relative}"
        own = subroutines if path == first else []
        stand_ins.append(StandIn(path, site, own, stand_in_text(path, site, own, first)))
    return stand_ins, [(site, path) for site, path, _called in scan.paths], ""


def survey(checkout: Path, stand_ins: list[StandIn]) -> tuple[set[str], list[str]]:
    """(the stand-ins already in the directory, by name and size; every other entry there). The
    entries are listed and measured, never opened."""
    directory = checkout / DIRECTORY
    if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
        return set(), [DIRECTORY]
    sizes = {stand_in.path: len(stand_in.text) for stand_in in stand_ins}
    parents = {parent.as_posix() for path in sizes for parent in PurePosixPath(path).parents}
    own: set[str] = set()
    foreign: list[str] = []
    for root, directories, names in os.walk(directory):
        for name in directories:
            path = Path(root) / name
            relative = path.relative_to(checkout).as_posix()
            if path.is_symlink() or relative not in parents:
                foreign.append(relative)
        for name in names:
            path = Path(root) / name
            relative = path.relative_to(checkout).as_posix()
            status = os.lstat(path)
            if stat.S_ISREG(status.st_mode) and sizes.get(relative) == status.st_size:
                own.add(relative)
            else:
                foreign.append(relative)
    return own, sorted(foreign)


def describe(stand_ins: list[StandIn], own: set[str], unread: list[tuple[Site, str]]) -> list[str]:
    """The lines of a dry run."""
    lines: list[str] = []
    for stand_in in stand_ins:
        verb = "would leave" if stand_in.path in own else "would write"
        lines.append(f"{stand_in.path}: {verb} {len(stand_in.text)} bytes; replaces the file {stand_in.site} requires")
        for name, empty, sites in stand_in.subroutines:
            where = ", ".join(str(site) for site in sites)
            what = "returns an empty list to the constructor" if empty else "ends the run"
            lines.append(f"  {name}: {what} (called at {where})")
        if not stand_in.subroutines:
            lines.append("  no subroutine")
    for site, path in unread:
        lines.append(f"{path}: named at {site}, read by no stand-in, not written")
    return lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("checkout", metavar="CHECKOUT", help=f"a LoFTEE checkout made without its {DIRECTORY} directory")
    parser.add_argument("--dry-run", action="store_true", help="print the files and their subroutines; write nothing")
    args = parser.parse_args(argv)
    checkout = Path(args.checkout)
    if not checkout.is_dir():
        print(f"ERROR: [{COMPONENT}] {checkout} is not a directory", file=sys.stderr)
        return 1
    stand_ins, unread, reason = plan(scan_sources(read_sources(checkout)))
    if reason:
        print(f"ERROR: [{COMPONENT}] {checkout}: {reason}", file=sys.stderr)
        return 1
    own, foreign = survey(checkout, stand_ins)
    if foreign:
        print(
            f"ERROR: [{COMPONENT}] {checkout} holds {foreign[0]}, which this script did not write "
            f"({len(foreign)} such under {DIRECTORY}); nothing is written",
            file=sys.stderr,
        )
        return 1
    if args.dry_run:
        print("\n".join(describe(stand_ins, own, unread)))
        return 0
    written = [stand_in for stand_in in stand_ins if stand_in.path not in own]
    for stand_in in written:
        path = checkout / stand_in.path
        path.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation fails on a file that appeared after the survey and never opens it.
        with open(path, "xb") as handle:
            handle.write(stand_in.text)
    wrote = ", ".join(f"{stand_in.path} ({len(stand_in.subroutines)} subroutines)" for stand_in in written)
    print(f"LoFTEE stand-ins: wrote {wrote or 'none'}; left {', '.join(sorted(own)) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
