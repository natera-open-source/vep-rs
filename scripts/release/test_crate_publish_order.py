"""Tests for crate_publish_order.py on `cargo metadata` documents and workflow texts the tests build, one
package per crate with the fields the script reads, and on the repository's own tree: its list `CRATES`
through cargo, the workspace version against the pins of `[workspace.dependencies]` as `Cargo.toml` states
them, and the toolchain of the workflow that publishes against the release workflow's."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "crate_publish_order.py"
ROOT = HERE.parent.parent
WORKFLOWS = ROOT / ".github" / "workflows"
PREFIX = "ERROR: [crate_publish_order] "


def dep(name: str, req: str = "=1.2.3", kind: str | None = None, path: bool = True) -> dict:
    """One entry of a package's `dependencies`: a workspace crate by path, or a registry crate."""
    entry = {"name": name, "req": req, "kind": kind}
    if path:
        entry["path"] = f"/w/crates/{name}"
    else:
        entry["source"] = "registry+https://github.com/rust-lang/crates.io-index"
    return entry


def package(name: str, deps: list[dict] | None = None, version: str = "1.2.3", publish: list[str] | None = None) -> dict:
    return {"name": name, "version": version, "id": f"path+file:///w/crates/{name}#{version}",
            "publish": publish, "dependencies": deps or []}


def workflow(crates: str) -> str:
    """A workflow text whose `env:` block states CRATES as a folded block scalar, as the tree's does."""
    return f"name: made up\nenv:\n  RUST_TOOLCHAIN: 1.2.3\n  CRATES: >-\n    {crates}\n\njobs:\n  verify:\n    steps: []\n"


def run(
    tmp_path: Path, packages: list[dict], crates: str | None = None, *, text: str | None = None,
    outside: tuple[dict, ...] = (),
) -> subprocess.CompletedProcess[str]:
    """The script over a metadata document whose workspace members are `packages` and a workflow that lists
    `crates` (the packages in their own order when not given); `text` is a whole workflow text instead."""
    document = {"packages": [*packages, *outside], "workspace_members": [p["id"] for p in packages]}
    (tmp_path / "metadata.json").write_text(json.dumps(document))
    listed = " ".join(p["name"] for p in packages) if crates is None else crates
    (tmp_path / "crates-io.yml").write_text(workflow(listed) if text is None else text)
    command = [sys.executable, str(SCRIPT), "--metadata", str(tmp_path / "metadata.json"),
               "--workflow", str(tmp_path / "crates-io.yml")]
    return subprocess.run(command, check=False, capture_output=True, text=True)


def passes(out: subprocess.CompletedProcess[str], count: int) -> None:
    assert out.returncode == 0 and out.stderr == "", out
    assert out.stdout == f"CRATES lists the workspace's {count} crates, each after its dependencies, each pin =1.2.3\n"


def refusals(out: subprocess.CompletedProcess[str]) -> list[str]:
    """The findings of a refused run, one per line of standard error."""
    assert out.returncode == 1 and out.stdout == "", out
    lines = out.stderr.splitlines()
    assert lines and all(line.startswith(PREFIX) for line in lines), out.stderr
    return [line.removeprefix(PREFIX) for line in lines]


APP = [
    package("app", [dep("engine"), dep("model")]),
    package("engine", [dep("model"), dep("reader")]),
    package("reader", [dep("model")]),
    package("model"),
]


def test_a_list_that_names_each_crate_after_its_dependencies_holds(tmp_path: Path) -> None:
    passes(run(tmp_path, APP, "model reader engine app"), 4)


def test_a_swapped_pair_is_refused_by_the_crate_and_the_dependency_it_precedes(tmp_path: Path) -> None:
    assert refusals(run(tmp_path, APP, "model engine reader app")) == ["CRATES lists engine before reader, which it depends on"]
    assert refusals(run(tmp_path, APP, "app engine reader model")) == [
        "CRATES lists app before engine, which it depends on",
        "CRATES lists app before model, which it depends on",
        "CRATES lists engine before model, which it depends on",
        "CRATES lists engine before reader, which it depends on",
        "CRATES lists reader before model, which it depends on",
    ]


def test_a_crate_the_list_leaves_out_is_refused(tmp_path: Path) -> None:
    assert refusals(run(tmp_path, APP, "model reader app")) == ["CRATES leaves out a crate the workspace publishes: engine"]


def test_a_name_the_workspace_does_not_publish_is_refused(tmp_path: Path) -> None:
    out = run(tmp_path, APP, "model reader engine app extra", outside=(package("extra"),))
    assert refusals(out) == ["CRATES lists a name the workspace does not publish: extra"]


def test_a_name_listed_twice_is_refused(tmp_path: Path) -> None:
    assert refusals(run(tmp_path, APP, "model reader engine app reader")) == ["CRATES lists more than once: reader"]


def test_every_finding_is_a_line_of_its_own(tmp_path: Path) -> None:
    packages = [package("a", [dep("b", req="^1.2.3")]), package("b"), package("c")]
    assert refusals(run(tmp_path, packages, "a b extra")) == [
        "CRATES leaves out a crate the workspace publishes: c",
        "CRATES lists a name the workspace does not publish: extra",
        "a requires b `^1.2.3`, not `=1.2.3`",
        "CRATES lists a before b, which it depends on",
    ]


def test_the_findings_of_two_crates_come_in_name_order_whatever_the_document_order(tmp_path: Path) -> None:
    packages = [package("b", [dep("z")]), package("a", [dep("z")]), package("z")]
    assert refusals(run(tmp_path, packages, "b a z")) == [
        "CRATES lists a before z, which it depends on",
        "CRATES lists b before z, which it depends on",
    ]


def test_a_build_and_a_versioned_development_dependency_order_a_crate_and_a_path_only_one_does_not(
    tmp_path: Path,
) -> None:
    packages = [
        package("a", [dep("x", kind="build")]),
        package("b", [dep("y", kind="dev")]),
        package("c", [dep("z", req="*", kind="dev")]),
        package("x"),
        package("y"),
        package("z"),
    ]
    passes(run(tmp_path, packages, "c x a y b z"), 6)
    assert refusals(run(tmp_path, packages, "a b c x y z")) == [
        "CRATES lists a before x, which it depends on",
        "CRATES lists b before y, which it depends on",
    ]


def test_a_registry_crate_named_as_a_workspace_crate_orders_nothing(tmp_path: Path) -> None:
    packages = [package("a", [dep("b", req="^9", path=False), dep("serde", req="^1", path=False)]), package("b", [dep("a")])]
    passes(run(tmp_path, packages, "a b"), 2)


def test_a_package_that_is_no_workspace_member_is_not_expected_in_the_list(tmp_path: Path) -> None:
    packages = [package("a", [dep("b")]), package("b")]
    passes(run(tmp_path, packages, "b a", outside=(package("serde", version="1.0.0"),)), 2)


def test_workspace_crates_at_two_versions_are_refused(tmp_path: Path) -> None:
    packages = [package("b", version="1.2.4"), package("a", [dep("b")])]
    assert refusals(run(tmp_path, packages)) == ["the workspace crates carry more than one version: a 1.2.3, b 1.2.4"]


def test_a_pin_that_is_not_the_exact_workspace_version_is_refused(tmp_path: Path) -> None:
    stale = [package("b"), package("a", [dep("b", req="=1.2.2")])]
    assert refusals(run(tmp_path, stale)) == ["a requires b `=1.2.2`, not `=1.2.3`"]
    loose = [package("b"), package("a", [dep("b", req="^1.2.3")]), package("c", [dep("b", req="=1.2.3", kind="dev")])]
    assert refusals(run(tmp_path, loose)) == ["a requires b `^1.2.3`, not `=1.2.3`"]
    bare = [package("b"), package("a", [dep("b", req="1.2.3")])]
    assert refusals(run(tmp_path, bare)) == ["a requires b `1.2.3`, not `=1.2.3`"]
    for kind in ("build", "dev"):
        stale = [package("b"), package("a", [dep("b", req="=1.2.2", kind=kind)])]
        assert refusals(run(tmp_path, stale)) == ["a requires b `=1.2.2`, not `=1.2.3`"]


def test_a_normal_or_build_dependency_by_path_alone_is_refused(tmp_path: Path) -> None:
    packages = [package("b"), package("a", [dep("b", req="*")]), package("c", [dep("b", req="*", kind="build")])]
    assert refusals(run(tmp_path, packages)) == [
        "a names b by path alone: a published normal dependency needs a version",
        "c names b by path alone: a published build dependency needs a version",
    ]


def test_a_cycle_through_a_versioned_development_dependency_is_refused_in_either_order(tmp_path: Path) -> None:
    packages = [package("a", [dep("b")]), package("b", [dep("a", kind="dev")]), package("c")]
    assert refusals(run(tmp_path, packages, "a b c")) == ["CRATES lists a before b, which it depends on"]
    assert refusals(run(tmp_path, packages, "b a c")) == ["CRATES lists b before a, which it depends on"]


def test_a_crate_kept_off_crates_io_is_not_listed_and_a_versioned_dependency_on_it_is_refused(tmp_path: Path) -> None:
    packages = [package("a"), package("tool", [dep("a")], publish=[]), package("b", publish=["crates-io"])]
    passes(run(tmp_path, packages, "a b"), 2)
    assert refusals(run(tmp_path, packages, "a tool b")) == ["CRATES lists a name the workspace does not publish: tool"]
    packages.append(package("c", [dep("tool")]))
    assert refusals(run(tmp_path, packages, "a b c")) == ["c depends on tool, which is not published to crates.io"]
    assert refusals(run(tmp_path, [package("tool", publish=["internal"])], "tool")) == ["the workspace publishes no crate"]


def test_the_list_is_read_as_a_block_or_on_the_key_line_and_refused_in_any_other_shape(tmp_path: Path) -> None:
    head, tail = "name: made up\nenv:\n  RUST_TOOLCHAIN: 1.2.3\n", "\njobs:\n  verify:\n    steps: []\n"
    for stated in (
        "  CRATES: >-\n    model reader\n    engine app\n",
        "  CRATES: |\n    model reader engine\n    app\n  OTHER: app model\n",
        "  CRATES: model reader engine app\n",
        '  CRATES: "model reader engine app"\n',
    ):
        passes(run(tmp_path, APP, text=head + stated + tail), 4)
    path = tmp_path / "crates-io.yml"
    for stated, finding in (
        ("", f"{path} does not state CRATES once as a key of its `env:` block"),
        ("  CRATES: model reader\n  CRATES: engine app\n", f"{path} does not state CRATES once as a key of its `env:` block"),
        ("      CRATES: model reader engine app\n", f"{path} does not state CRATES once as a key of its `env:` block"),
        ("  CRATES: >-\n", f"{path} states CRATES in a shape this script does not read: ''"),
        ("  CRATES: [model, reader, engine, app]\n",
         f"{path} states CRATES in a shape this script does not read: '[model, reader, engine, app]'"),
    ):
        assert refusals(run(tmp_path, APP, text=head + stated + tail)) == [finding]


def test_a_missing_workflow_or_an_unreadable_metadata_file_is_one_error_line(tmp_path: Path) -> None:
    (tmp_path / "broken.json").write_text("{")
    (tmp_path / "crates-io.yml").write_text(workflow("a"))
    out = subprocess.run([sys.executable, str(SCRIPT), "--metadata", str(tmp_path / "broken.json"),
                          "--workflow", str(tmp_path / "crates-io.yml")], check=False, capture_output=True, text=True)
    (finding,) = refusals(out)
    assert finding.startswith(f"cannot read {tmp_path / 'broken.json'}: ")
    out = subprocess.run([sys.executable, str(SCRIPT), "--metadata", str(tmp_path / "broken.json"),
                          "--workflow", str(tmp_path / "absent.yml")], check=False, capture_output=True, text=True)
    (finding,) = refusals(out)
    assert finding.startswith(f"cannot read {tmp_path / 'absent.yml'}: ")
    (tmp_path / "bytes.yml").write_bytes(b"env:\n  CRATES: a \xff\n")
    out = subprocess.run([sys.executable, str(SCRIPT), "--metadata", str(tmp_path / "broken.json"),
                          "--workflow", str(tmp_path / "bytes.yml")], check=False, capture_output=True, text=True)
    (finding,) = refusals(out)
    assert finding.startswith(f"cannot read {tmp_path / 'bytes.yml'}: ")


def test_the_repository_lists_every_workspace_crate_after_its_dependencies(tmp_path: Path) -> None:
    # Run outside the repository with no option, so the tree is found from the script's own place.
    out = subprocess.run([sys.executable, str(SCRIPT)], check=False, capture_output=True, text=True, cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    manifest = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    version = manifest["workspace"]["package"]["version"]
    members = sorted(Path(member).name for member in manifest["workspace"]["members"])
    assert out.stdout == (
        f"CRATES lists the workspace's {len(members)} crates, each after its dependencies, each pin ={version}\n")
    # The list as the workflow states it, read here a second way.
    text = (WORKFLOWS / "crates-io.yml").read_text(encoding="utf-8")
    block = re.search(r"^  CRATES: >-\n((?:    \S.*\n)+)", text, re.MULTILINE)
    assert block, "crates-io.yml does not state CRATES as a folded block of its `env:` block"
    order = block.group(1).split()
    assert sorted(order) == members
    assert order[0] == "vep-core"
    at = {name: n for n, name in enumerate(order)}
    for crate in order:
        crate_manifest = tomllib.loads((ROOT / "crates" / crate / "Cargo.toml").read_text(encoding="utf-8"))
        for table in ("dependencies", "build-dependencies", "dev-dependencies"):
            for name, entry in crate_manifest.get(table, {}).items():
                # An entry that inherits the workspace's pin or states a version stays in the published manifest.
                if name in at and isinstance(entry, dict) and (entry.get("workspace") or "version" in entry):
                    assert at[name] < at[crate], f"{crate} is listed before its {table} entry {name}"


def test_every_workspace_crate_pin_is_the_exact_workspace_version() -> None:
    manifest = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    version = manifest["workspace"]["package"]["version"]
    pins = {name: entry["version"] for name, entry in manifest["workspace"]["dependencies"].items()
            if isinstance(entry, dict) and "path" in entry}
    assert pins, "no workspace crate is pinned under [workspace.dependencies]"
    assert pins == dict.fromkeys(pins, f"={version}")
    for crate in (Path(member).name for member in manifest["workspace"]["members"]):
        crate_manifest = tomllib.loads((ROOT / "crates" / crate / "Cargo.toml").read_text(encoding="utf-8"))
        assert crate_manifest["package"]["version"] == {"workspace": True}, crate
        for table in ("dependencies", "build-dependencies", "dev-dependencies"):
            for name, entry in crate_manifest.get(table, {}).items():
                if name in pins:
                    # A development dependency that names a path alone carries no pin: `cargo publish` leaves
                    # it out of the uploaded manifest, as vep-cache-builder's on vep-cli, which depends on it.
                    if table == "dev-dependencies" and isinstance(entry, dict) and set(entry) == {"path"}:
                        continue
                    assert isinstance(entry, dict) and entry.get("workspace") is True, f"{crate} {table} {name}"


def test_the_publishing_workflow_installs_the_release_workflow_toolchain() -> None:
    stated = {}
    for name in ("crates-io.yml", "release.yml"):
        text = (WORKFLOWS / name).read_text(encoding="utf-8")
        values = re.findall(r"^  RUST_TOOLCHAIN: (\S+)$", text, re.MULTILINE)
        assert len(values) == 1, f"{name} states RUST_TOOLCHAIN {len(values)} times in its `env:` block"
        assert re.fullmatch(r"\d+\.\d+\.\d+", values[0]), f"{name} states RUST_TOOLCHAIN as {values[0]}"
        assert text.count("toolchain: ${{ env.RUST_TOOLCHAIN }}") == text.count("dtolnay/rust-toolchain@"), name
        stated[name] = values[0]
    assert stated["crates-io.yml"] == stated["release.yml"]
