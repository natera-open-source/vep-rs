#!/usr/bin/env python3
"""Hold the list `CRATES` of `.github/workflows/crates-io.yml` to the workspace it publishes.

The workflow uploads the crates in the order that list gives them, and crates.io resolves every dependency
of an uploaded crate against the registry. So the list has to name exactly the workspace crates crates.io
takes, each after every workspace crate it depends on. A normal or a build dependency counts, and so does a
development dependency that carries a version requirement, which `cargo publish` keeps in the uploaded
manifest (it drops one that names a path alone).

The workspace is read from `cargo metadata --format-version 1 --no-deps --locked`. `--no-deps` resolves
nothing, so cargo itself does not notice there a version that moved without its pins: the workspace is
therefore also held to carry one version, every dependency of one of its crates on another to read `=` that
version, and a normal or build dependency to carry a version at all.

Every finding is one line on standard error, `ERROR: [crate_publish_order] <finding>`, and the exit status
is 1 when there is one. A list that holds prints one line on standard output.

usage: crate_publish_order.py [--manifest-path Cargo.toml] [--workflow crates-io.yml] [--metadata FILE]

The manifest and the workflow default to those of the tree this script stands in. `--metadata FILE` reads
the `cargo metadata` JSON from FILE instead of running cargo.

No dependency beyond the standard library.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

CRATES_IO = "crates-io"
ROOT = Path(__file__).resolve().parent.parent.parent
WORKFLOW = ROOT / ".github" / "workflows" / "crates-io.yml"


def fail(*findings: str) -> None:
    for finding in findings:
        print(f"ERROR: [crate_publish_order] {finding}", file=sys.stderr)
    sys.exit(1)


def cargo_metadata(manifest_path: Path) -> dict:
    command = ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked", "--manifest-path", str(manifest_path)]
    try:
        done = subprocess.run(command, check=False, capture_output=True, text=True)
    except FileNotFoundError:
        fail("cargo is not on PATH")
    if done.returncode != 0:
        fail(f"cargo metadata failed: {' '.join(done.stderr.split())}")
    return json.loads(done.stdout)


def listed_crates(workflow: Path) -> list[str]:
    """The names of the workflow's `CRATES`, in its order.

    The value is read where it stands as a key of the top-level `env:` block (two spaces of indentation): a
    plain scalar on the key's line, or a folded or literal block scalar on the more indented lines after it.
    """
    try:
        lines = workflow.read_text(encoding="utf-8").splitlines()
    except (OSError, ValueError) as exc:
        fail(f"cannot read {workflow}: {exc}")
    at = [n for n, line in enumerate(lines) if re.match(r"  CRATES:(\s|$)", line)]
    if len(at) != 1:
        fail(f"{workflow} does not state CRATES once as a key of its `env:` block")
    value = lines[at[0]].split(":", 1)[1].strip()
    if re.fullmatch(r"[>|][+-]?", value):
        block = []
        for line in lines[at[0] + 1:]:
            if line.strip() and not line.startswith("    "):
                break
            block.append(line)
        value = " ".join(block)
    names = value.strip("\"'").split()
    if not names or not all(re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in names):
        fail(f"{workflow} states CRATES in a shape this script does not read: {value.strip()!r}")
    return names


def published(package: dict) -> bool:
    """`publish` is null for a crate any registry takes, else the list of registries that take it."""
    return package["publish"] is None or CRATES_IO in package["publish"]


def findings(metadata: dict, listed: list[str]) -> tuple[list[str], str]:
    """Each way `listed` or the workspace's pins would fail a release, and the workspace version."""
    member_ids = set(metadata["workspace_members"])
    members = {p["name"]: p for p in metadata["packages"] if p["id"] in member_ids}
    crates = {name: p for name, p in members.items() if published(p)}
    if not crates:
        return ["the workspace publishes no crate"], ""
    problems: list[str] = []

    versions = sorted({p["version"] for p in crates.values()})
    version = versions[0]
    if len(versions) > 1:
        problems.append("the workspace crates carry more than one version: " + ", ".join(
            f"{name} {p['version']}" for name, p in sorted(crates.items())))

    twice = sorted({name for name in listed if listed.count(name) > 1})
    if twice:
        problems.append("CRATES lists more than once: " + ", ".join(twice))
    unlisted = sorted(set(crates) - set(listed))
    if unlisted:
        problems.append("CRATES leaves out a crate the workspace publishes: " + ", ".join(unlisted))
    unknown = sorted(set(listed) - set(crates))
    if unknown:
        problems.append("CRATES lists a name the workspace does not publish: " + ", ".join(unknown))

    place = {name: n for n, name in reversed(list(enumerate(listed)))}
    for name, package in sorted(crates.items()):
        for dep in package["dependencies"]:
            # A workspace crate is a dependency by path; a registry crate of the same name is not one.
            if dep["name"] not in members or not dep.get("path"):
                continue
            other, kind, req = dep["name"], dep["kind"] or "normal", dep["req"]
            if req == "*":
                if kind != "dev":
                    problems.append(f"{name} names {other} by path alone: a published {kind} dependency needs a version")
                continue
            if other not in crates:
                problems.append(f"{name} depends on {other}, which is not published to crates.io")
                continue
            if len(versions) == 1 and req != f"={version}":
                problems.append(f"{name} requires {other} `{req}`, not `={version}`")
            if name in place and other in place and place[name] <= place[other]:
                problems.append(f"CRATES lists {name} before {other}, which it depends on")
    return problems, version


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-path", type=Path, default=ROOT / "Cargo.toml", help="the workspace's Cargo.toml")
    ap.add_argument("--workflow", type=Path, default=WORKFLOW, help="the workflow file that states CRATES")
    ap.add_argument("--metadata", type=Path, help="read the `cargo metadata` JSON here instead of running cargo")
    a = ap.parse_args()
    listed = listed_crates(a.workflow)
    if a.metadata:
        try:
            metadata = json.loads(a.metadata.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            fail(f"cannot read {a.metadata}: {exc}")
    else:
        metadata = cargo_metadata(a.manifest_path)
    problems, version = findings(metadata, listed)
    if problems:
        fail(*problems)
    print(f"CRATES lists the workspace's {len(listed)} crates, each after its dependencies, each pin ={version}")


if __name__ == "__main__":
    main()
