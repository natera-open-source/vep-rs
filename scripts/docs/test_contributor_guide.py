"""The contributor guide quotes `.github/workflows/ci.yml` as the workflow stands.

The guide at the repository root caches the gate a pull request must pass: CI's pinned toolchain,
the minimum supported Rust version, the `--locked` flag, the clippy feature flag, the API
documentation gate and the per-file test floors of the Comparator Tests job. The workflow is the authority, so these tests read it as
text (the standard library alone reads it) and hold every value the guide quotes to it.
They also hold the workflow's own Comparator Tests list to one shape, every floored file run and
every run file floored, and this file's floor to the number of tests it carries. A Cargo root under
`tools/` is outside every `--workspace` command, so they hold the workflow to an audit of each
root's lockfile and a build of each root on the minimum supported Rust version."""

from __future__ import annotations

import datetime
import inspect
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CI = ROOT / ".github" / "workflows" / "ci.yml"
GUIDE = ROOT / "AGENTS.md"
CONTRIBUTING = ROOT / "CONTRIBUTING.md"
MANIFESTS = (ROOT / "Cargo.toml", *sorted((ROOT / "tools").glob("*/Cargo.toml")))
THIS_FILE = "scripts/docs/test_contributor_guide.py"

TOOLCHAIN = re.compile(r'^\s*toolchain:\s*"?([0-9][0-9.]*)"?\s*$', re.MULTILINE)
RUST_VERSION = re.compile(r'^rust-version\s*=\s*"([^"]+)"', re.MULTILINE)
# A floor entry of the Comparator Tests list, and the same shape wherever the guide quotes one.
FLOOR_ENTRY = re.compile(r'"(scripts/\S+\.py) (\d+)"')
# A harness test file's floor, `check <file> <floor>`, in the workflow and in a guide that quotes one.
HARNESS_FLOOR = re.compile(r"\bcheck (scripts/\S+\.sh) (\d+)\b")
MATRIX_COMMAND = re.compile(r"^\s*command:\s*(.+?)\s*$", re.MULTILINE)
MATRIX_COMPARATOR = re.compile(r"^\s*comparator:\s*(.+?)\s*$", re.MULTILINE)
MANIFEST_PATH = re.compile(r"--manifest-path (tools/[^/\s]+)/Cargo\.toml")
PIP_INSTALL = re.compile(r"pip install --quiet (.+?)\s*$", re.MULTILINE)
COLLECT_PATTERN = re.compile(r"grep -cE '(\^scripts/\.\*::)'")
OK_PATTERN = re.compile(r"grep -c '(\^  ok )'")
COMMENT_LINE = re.compile(r"^\s*#.*$", re.MULTILINE)
# A lockfile audit of one tool root, whatever flags stand before `--file`.
AUDIT_FILE = re.compile(r"^cargo audit (?:\S+ )*--file (tools/[^/\s]+)/Cargo\.lock$")
# An advisory one audit passes over, and the allowlist that dates it for the workspace.
AUDIT_IGNORE = re.compile(r"--ignore (RUSTSEC-\d{4}-\d{4}) ")
ADVISORY_SCRIPT = "scripts/check_advisory_reachability.sh"
API_DOCS_SCRIPT = "scripts/check_api_docs.sh"
ADVISORY_ALLOWLIST = ROOT / "scripts" / "advisory-reachability-allow.txt"
RUN_LINE = re.compile(r"^\s*(?:- )?run: (\S.*?)\s*$", re.MULTILINE)
# A block scalar's indicator, `|` or `>` with its chomping and indentation marks and the comment
# YAML allows after them: the lines under it are the command, so the indicator is none.
BLOCK_MARK = re.compile(r"[|>][-+0-9]*(?:\s+#.*)?")


def ci_text() -> str:
    return CI.read_text(encoding="utf-8")


def ci_code() -> str:
    """The workflow with its comment lines removed: a value quoted in a comment is not one the job
    runs, so every pattern that reads a command or a floor reads this text."""
    return COMMENT_LINE.sub("", ci_text())


def guide_text() -> str:
    return GUIDE.read_text(encoding="utf-8")


def rust_version() -> str:
    """The one `rust-version` every manifest declares."""
    declared = {RUST_VERSION.search(m.read_text(encoding="utf-8")).group(1) for m in MANIFESTS}
    assert len(declared) == 1, f"the manifests declare several rust-versions: {sorted(declared)}"
    return declared.pop()


def toolchains() -> tuple[str, str]:
    """The pinned toolchain and the MSRV toolchain the workflow's `toolchain:` inputs name."""
    msrv = rust_version()
    values = set(TOOLCHAIN.findall(ci_code()))
    assert msrv in values, f"no toolchain input names the declared rust-version {msrv}: {sorted(values)}"
    pinned = values - {msrv}
    assert len(pinned) == 1, f"the workflow pins several toolchains beside the MSRV: {sorted(pinned)}"
    return pinned.pop(), msrv


def comparator_block() -> str:
    """The Comparator Tests job's text, comment lines removed, from its key to the next job's."""
    text = ci_code()
    start = text.index("\n  comparator-tests:\n")
    end = text.index("\n  harness-tests:\n", start)
    return text[start:end]


def floor_list() -> dict[str, int]:
    """The `tests=( ... )` entries of the Comparator Tests job, file to floor."""
    block = comparator_block()
    array = block[block.index("tests=(") : block.index(")", block.index("tests=("))]
    entries = FLOOR_ENTRY.findall(array)
    files = [f for f, _ in entries]
    assert len(files) == len(set(files)), "a file is floored twice"
    return {f: int(n) for f, n in entries}


def run_list() -> list[str]:
    """The files an explicit `python3 -m pytest -q \\` continuation names, empty when the job runs
    the floor list's own array (`"${files[@]}"`)."""
    block = comparator_block()
    match = re.search(r"python3 -m pytest -q \\\n((?:\s*scripts/\S+ ?\\?\n)+)", block)
    if match is None:
        assert 'files+=("$file")' in block, "the floored array is not gathered into the run"
        assert 'python3 -m pytest -q "${files[@]}"' in block, (
            "the job neither runs the floored array nor an explicit list"
        )
        return []
    return re.findall(r"(scripts/\S+?)\s*\\?\n", match.group(1))


def matrix_command(name: str) -> str:
    """The `command:` of the check matrix entry named `name`."""
    text = ci_code()
    entry = text[text.index(f"- name: {name}\n") :]
    return MATRIX_COMMAND.search(entry).group(1)


def matrix_comparator(name: str) -> str:
    """The `comparator:` of the check matrix entry named `name`, its command over the tool roots."""
    text = ci_code()
    entry = text[text.index(f"- name: {name}\n") :]
    return MATRIX_COMPARATOR.search(entry).group(1)


def msrv_build_commands() -> list[str]:
    """The `cargo build` lines of the MSRV job, whether `run:` opens its step or follows the step's `if:`."""
    text = ci_code()
    start = text.index("\n  msrv:\n")
    end = text.index("\n  build-release:\n", start)
    builds = re.findall(r"^\s*(?:- )?run: (cargo build .+?)\s*$", text[start:end], re.MULTILINE)
    assert builds, "the MSRV job runs no cargo build"
    return builds


def msrv_loop_manifests() -> list[str]:
    """The manifests whose `rust-version` the MSRV job holds to its toolchain: the files of the
    tree each word of the loop's list names, a `*` standing for one directory as the shell reads it."""
    text = ci_code()
    start = text.index("\n  msrv:\n")
    end = text.index("\n  build-release:\n", start)
    loop = re.search(r"^\s*for manifest in (.+?); do\s*$", text[start:end], re.MULTILINE)
    assert loop, "the MSRV job loops over no manifest"
    return [path.relative_to(ROOT).as_posix() for word in loop.group(1).split() for path in sorted(ROOT.glob(word))]


def allowlisted_advisories() -> dict[str, datetime.date]:
    """The advisories of the dated allowlist, each with the day its entry expires, read as
    `scripts/check_advisory_reachability.sh` reads the file."""
    entries = {}
    for raw in ADVISORY_ALLOWLIST.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            advisory, expiry, _reason = (part.strip() for part in line.split("|"))
            entries[advisory] = datetime.date.fromisoformat(expiry)
    return entries


def one_line_commands(job: str) -> list[str]:
    """The commands the steps of `job` give on their `run:` lines; a step whose `run:` opens a
    block of several lines gives none."""
    return [command for command in RUN_LINE.findall(job) if not BLOCK_MARK.fullmatch(command)]


def advisory_commands() -> list[str]:
    """The one-line `run:` commands of the Advisory Reachability job."""
    text = ci_code()
    start = text.index("\n  advisory-reachability:\n")
    end = text.index("\n  msrv:\n", start)
    return one_line_commands(text[start:end])


def tool_roots() -> list[str]:
    """Every directory under `tools/` that carries a lockfile of its own, as `tools/<name>`."""
    roots = sorted(f"tools/{lock.parent.name}" for lock in (ROOT / "tools").glob("*/Cargo.lock"))
    assert roots, "no directory under tools/ carries a Cargo.lock"
    return roots


def test_the_workflow_pins_one_toolchain_beside_the_msrv():
    pinned, msrv = toolchains()
    assert pinned != msrv
    assert msrv == rust_version()


def test_the_guide_quotes_the_pinned_toolchain_and_the_msrv():
    pinned, msrv = toolchains()
    guide = guide_text()
    assert f"Rust {pinned}" in guide
    assert f"Rust {msrv}" in guide
    assert f"rustup toolchain install {pinned} {msrv}" in guide


def test_the_guide_names_only_the_workflow_s_toolchain_versions():
    """A version the guide spells `Rust <n>` or `+<n>` is one the workflow names."""
    pinned, msrv = toolchains()
    version = r"(1\.[0-9]+(?:\.[0-9]+)?)"
    named = set(re.findall(r"\bRust " + version, guide_text()))
    named |= set(re.findall(r"cargo \+" + version + " ", guide_text()))
    assert named == {pinned, msrv}, f"the guide names {sorted(named)}, the workflow {pinned} and {msrv}"


def test_ci_tests_and_builds_locked():
    assert "--locked" in matrix_command("Test")
    for command in msrv_build_commands():
        assert "--locked" in command, command


def test_the_guide_quotes_the_test_command_with_locked():
    """The guide's `cargo test` line is the Test job's command, with `--no-fail-fast` as its one addition."""
    test = matrix_command("Test")
    assert "--locked" in test
    quoted = re.findall(r"^(cargo test --workspace .+)$", guide_text(), re.MULTILINE)
    assert quoted, "the guide quotes no workspace cargo test line"
    for line in quoted:
        assert line.removesuffix(" --no-fail-fast") == test, f"{line!r} is not the Test job's {test!r}"


def test_the_guide_quotes_the_clippy_command_with_its_feature():
    clippy = matrix_command("Clippy")
    assert "--features vep-cli/bench" in clippy
    assert clippy in guide_text()
    assert clippy in CONTRIBUTING.read_text(encoding="utf-8")


def test_the_guides_quote_the_api_documentation_gate():
    """The Module Docs entry runs the API documentation gate over the workspace and over the tool
    roots, and `AGENTS.md` and `CONTRIBUTING.md` quote both forms."""
    workspace, tools = f"{API_DOCS_SCRIPT} workspace", f"{API_DOCS_SCRIPT} tools"
    assert workspace in matrix_command("Module Docs").split(" && ")
    assert tools in matrix_comparator("Module Docs").split(" && ")
    for page in (GUIDE, CONTRIBUTING):
        text = page.read_text(encoding="utf-8")
        for form in (workspace, tools):
            assert f"`{form}`" in text or f"\n{form}\n" in text, f"{page.name} does not quote {form!r}"


def test_the_guide_quotes_the_msrv_build_command():
    """The guide's `cargo +<msrv> build` line is the MSRV job's workspace build, whole."""
    _, msrv = toolchains()
    workspace = [c for c in msrv_build_commands() if "--workspace" in c]
    assert len(workspace) == 1, workspace
    quoted = re.findall(rf"^cargo \+{re.escape(msrv)} (build .+)$", guide_text(), re.MULTILINE)
    assert quoted == [workspace[0].removeprefix("cargo ")], f"the guide quotes {quoted}, the workflow {workspace[0]!r}"


def test_the_guide_names_every_tool_cargo_root_ci_builds():
    roots = set(MANIFEST_PATH.findall(ci_code()))
    assert roots, "the workflow builds no tool root"
    guide = guide_text()
    for root in roots:
        assert f"`{root}`" in guide, root


def test_the_workflow_audits_and_builds_every_tool_cargo_root():
    """Each tool root's lockfile is audited, the root is built on the MSRV and its manifest's
    `rust-version` is held to that toolchain; the MSRV job's workspace build carries the feature
    that compiles the `bench_plugins` binary."""
    audited = {m.group(1) for m in map(AUDIT_FILE.match, advisory_commands()) if m}
    builds = msrv_build_commands()
    built = {root for command in builds for root in MANIFEST_PATH.findall(command)}
    looped = msrv_loop_manifests()
    assert "Cargo.toml" in looped, f"the rust-version loop leaves out the workspace manifest: {looped}"
    for root in tool_roots():
        assert root in audited, f"no Advisory Reachability step audits {root}/Cargo.lock"
        assert root in built, f"no MSRV step builds {root}"
        assert f"{root}/Cargo.toml" in looped, f"the rust-version loop leaves out {root}/Cargo.toml"
    workspace = [c for c in builds if "--workspace" in c]
    assert len(workspace) == 1, workspace
    assert "--features vep-cli/bench" in workspace[0], workspace[0]


def test_the_guides_quote_the_advisory_job_s_two_commands():
    """`AGENTS.md` and `CONTRIBUTING.md` quote the reachability script the job runs over the
    workspace's lockfile and the one audit command it runs over every tool root's, the root
    written `<name>`. A row differs from that command only by an `--ignore` of an advisory the
    dated allowlist carries for the workspace in an entry that has not expired."""
    commands = advisory_commands()
    assert ADVISORY_SCRIPT in commands, f"the Advisory Reachability job does not run {ADVISORY_SCRIPT}"
    audits = [c for c in commands if AUDIT_FILE.match(c)]
    forms = {re.sub(r"tools/[^/\s]+/", "tools/<name>/", AUDIT_IGNORE.sub("", c)) for c in audits}
    assert len(forms) == 1, f"the job audits the tool lockfiles by {len(forms)} commands, not one: {sorted(forms)}"
    form = forms.pop()
    plain = {AUDIT_IGNORE.sub("", c) for c in audits}
    for root in tool_roots():
        assert form.replace("tools/<name>", root) in plain, f"{root}/Cargo.lock is not audited as {form!r}"
    allowed = allowlisted_advisories()
    for command in audits:
        for advisory in AUDIT_IGNORE.findall(command):
            assert advisory in allowed, f"{ADVISORY_ALLOWLIST.name} does not carry {advisory}: {command!r}"
            assert allowed[advisory] >= datetime.date.today(), f"the entry of {advisory} expired {allowed[advisory]}"
    for page in (GUIDE, CONTRIBUTING):
        text = page.read_text(encoding="utf-8")
        assert f"`{ADVISORY_SCRIPT}`" in text, f"{page.name} does not quote {ADVISORY_SCRIPT}"
        assert f"`{form}`" in text, f"{page.name} does not quote {form!r}"


def test_the_advisory_reader_returns_commands_alone():
    """A step whose `run:` opens a block gives the reader no entry, so the job's list is its
    one-line commands: the reachability script once and one lockfile audit for each tool root."""
    job = (
        "    steps:\n"
        "      - run: |\n"
        "          set -euo pipefail\n"
        "          cargo audit --version\n"
        "      - if: needs.changes.outputs.rust != 'false'\n"
        "        run: >-\n"
        "          cargo audit\n"
        "          --version\n"
        "      - run: | # a comment after the indicator\n"
        "          cargo audit --version\n"
        "      - run: scripts/one.sh\n"
        "      - if: needs.changes.outputs.rust != 'false'\n"
        "        run: cargo audit --file tools/one/Cargo.lock\n"
    )
    assert one_line_commands(job) == ["scripts/one.sh", "cargo audit --file tools/one/Cargo.lock"]
    commands = advisory_commands()
    for command in commands:
        assert not BLOCK_MARK.fullmatch(command), f"{command!r} opens a block and is no command"
    assert commands.count(ADVISORY_SCRIPT) == 1, commands
    audited = sorted(m.group(1) for m in map(AUDIT_FILE.match, commands) if m)
    assert audited == tool_roots(), f"the job audits {audited}, the tree holds {tool_roots()}"


def test_the_manifest_list_is_the_workspace_and_every_tool_root():
    """The manifests held to one `rust-version` are the workspace's and that of each tool root, so
    a root with a lockfile and no manifest, or a manifest and no lockfile, fails here."""
    expected = [ROOT / "Cargo.toml", *(ROOT / root / "Cargo.toml" for root in tool_roots())]
    assert list(MANIFESTS) == expected, f"{[m.relative_to(ROOT).as_posix() for m in MANIFESTS]}"


def test_the_guides_name_the_advisory_allow_list():
    """`AGENTS.md` and `CONTRIBUTING.md` quote the allowlist of dated entries, and it is the file
    the reachability script reads unless `ALLOWLIST` names another."""
    name = ADVISORY_ALLOWLIST.relative_to(ROOT).as_posix()
    script = (ROOT / ADVISORY_SCRIPT).read_text(encoding="utf-8")
    assert f'ALLOWLIST="${{ALLOWLIST:-{name}}}"' in script, f"{ADVISORY_SCRIPT} does not default to {name}"
    for page in (GUIDE, CONTRIBUTING):
        assert f"`{name}`" in page.read_text(encoding="utf-8"), f"{page.name} does not quote {name}"


def test_the_guide_quotes_the_python_install_line():
    installs = set(PIP_INSTALL.findall(ci_code()))
    assert len(installs) == 1, f"the workflow installs several Python sets: {sorted(installs)}"
    assert f"pip install {installs.pop()}" in guide_text()


def test_the_floor_list_feeds_the_run():
    floored = set(floor_list())
    run = run_list()
    if run:
        assert set(run) == floored
        assert len(run) == len(set(run)), "a file is run twice"


def test_every_floor_the_guide_quotes_is_the_workflow_s():
    floors = floor_list()
    quoted = FLOOR_ENTRY.findall(guide_text())
    assert quoted, "the guide quotes no floor entry"
    for file, floor in quoted:
        assert file in floors, f"the guide floors {file}, which the workflow does not list"
        assert int(floor) == floors[file], f"{file}: the guide says {floor}, the workflow {floors[file]}"
    harness = dict(HARNESS_FLOOR.findall(ci_code()))
    assert harness, "the Harness Tests job floors no file"
    for file, floor in HARNESS_FLOOR.findall(guide_text()):
        assert harness.get(file) == floor, f"{file}: the guide says {floor}, the workflow {harness.get(file)}"


def test_the_guide_quotes_this_file_s_floor():
    floors = floor_list()
    assert THIS_FILE in floors
    assert f'"{THIS_FILE} {floors[THIS_FILE]}"' in guide_text()


def test_this_file_s_floor_is_its_test_count():
    module = sys.modules[__name__]
    tests = [n for n, f in inspect.getmembers(module, inspect.isfunction) if n.startswith("test_")]
    assert floor_list()[THIS_FILE] == len(tests)


def test_the_guide_quotes_the_workflow_s_floor_count_patterns():
    """The two floor-measuring commands the guide quotes are the ones the jobs run: the pytest
    collected-count pattern and the bash ok-line pattern."""
    collected = set(COLLECT_PATTERN.findall(ci_code()))
    assert len(collected) == 1, sorted(collected)
    assert f"grep -cE '{collected.pop()}'" in guide_text()
    ok = set(OK_PATTERN.findall(ci_code()))
    assert len(ok) == 1, sorted(ok)
    assert f"grep -c '{ok.pop()}'" in guide_text()
