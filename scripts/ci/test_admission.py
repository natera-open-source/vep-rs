"""Tests for the admission check, on the recorded shapes of this repository's pull requests (the
fixtures beside this file: numbers, logins, names, addresses and model names are placeholders, hunk
text is reduced to its tags, the test regions of the changed Rust files are recorded from the files
themselves), on the constructed shapes of the intended-divergences section, and on constructed
inputs. Every check fails on its recorded or constructed shape and passes once the shape is
corrected; the comment is edited in place; a stale CI run judges nothing; each exemption holds the
close back; the API writes of a run are the expected ones and no other. The two tests that read the
repository's own pull request template and issue forms are skipped only when `ADMISSION_RELEASE_TREE=1`
is set, for a run of this file over a tree whose template and forms predate the rule."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
import urllib.error
from pathlib import Path

import pytest
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import admission as adm

RELEASE_TREE = pytest.mark.skipif(
    os.environ.get("ADMISSION_RELEASE_TREE") == "1",
    reason="the release tree's template and forms predate the rule until the default branch merges in",
)
FIXTURE_DIR = HERE / "fixtures" / "admission"
FIXTURES = sorted(FIXTURE_DIR.glob("*.json"))
FIXTURE_IDS = [path.stem for path in FIXTURES]
FIXTURE_COUNT = 27
INDEX_SAMPLE = (FIXTURE_DIR / "index-sample-main.md").read_text(encoding="utf-8")
INDEX_SAMPLE_REGROUPED = (FIXTURE_DIR / "index-sample-regrouped.md").read_text(encoding="utf-8")
CHECKER = ROOT / adm.PARITY_CHECKER
TEMPLATE = ROOT / ".github" / "pull_request_template.md"
FORMS_DIR = ROOT / ".github" / "ISSUE_TEMPLATE"
MATRIX_HEADER = "item\tperl_reference\ttarget\tstate\ttest1\ttest2\tinterim\n"
# The plugins page as the parity checker reads it: the support table heading over an empty table,
# which asks nothing of a matrix row that is not a plugin's.
PLUGINS_PAGE = (
    "# Plugins\n\n## Support table\n\n"
    "| Plugin | Data file and provider | Terms | Assemblies | Status | Reason |\n"
    "| --- | --- | --- | --- | --- | --- |\n"
)
MAINTAINERS = {"first-maintainer", "second-maintainer"}
REPO = "example/vep-rs"

# The misses each fixture carries, by check name. The ten recorded pull requests were opened before
# the divergences rule, so their descriptions lack its section and that miss is reported without a
# close (`opened_before_divergences_rule` in the fixture); the constructed `divergences-*` shapes
# were opened after it.
RECORDED = {
    "bot-dependency-bump": {"base", "signoff", "declaration", "body", "divergences"},
    "fork-many-commits-no-declaration": {"declaration", "divergences"},
    "fork-performance-change": {"declaration", "divergences"},
    "maintainer-docs-with-declaration": {"divergences"},
    "maintainer-flag-bam": {"divergences"},
    "maintainer-flag-gff-gtf": {"divergences"},
    "maintainer-guide-added": {"divergences"},
    "maintainer-guide-refined": {"divergences"},
    "maintainer-local-gate-fixes": {"changelog", "divergences"},
    "maintainer-manifests-merge-commit": {"divergences"},
    "divergences-none-clean": set(),
    "divergences-none-edits-a-page": {"divergences"},
    "divergences-none-edits-a-linked-page": {"divergences"},
    "divergences-touches-letter": set(),
    "divergences-touches-difference-name": set(),
    "divergences-touches-page-path": set(),
    "divergences-touches-unknown-entry": {"divergences"},
    "divergences-touches-without-rationale": {"divergences"},
    "divergences-touches-no-entry-line": {"divergences"},
    "divergences-none-and-touches": {"divergences"},
    "divergences-none-and-proposes": {"divergences"},
    "divergences-nothing-checked": {"divergences"},
    "divergences-no-section": {"divergences"},
    "divergences-regrouped-index": set(),
    "divergences-proposes-with-description": set(),
    "divergences-proposes-without-description": {"divergences"},
    "divergences-touches-and-proposes": set(),
}
NONE_SECTION = f"\n## {adm.DIVERGENCE_SECTION}\n- [x] {adm.NONE_STATEMENT}\n"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def index_text(fixture: dict, key: str) -> str | None:
    name = fixture.get(key)
    return (FIXTURE_DIR / name).read_text(encoding="utf-8") if name else None


def regions_of(fixture: dict) -> adm.RegionsFor:
    """The recorded test regions of a fixture's Rust files, as the job reads them from git."""
    table = {}
    for entry in fixture["files"]:
        for side, regions in (entry.get("regions") or {}).items():
            path = entry.get("previous_filename") if side == "old" else None
            table[(side, path or entry["filename"])] = [tuple(r) for r in regions]
    return lambda side, path: table.get((side, path), [])


def findings_of(
    fixture: dict, labels: set[str] | None = None, parity: adm.Finding | None = None
) -> list[adm.Finding]:
    pr = fixture["pull_request"]
    return adm.static_findings(
        fixture["contains_base_tip"],
        pr["base"]["ref"],
        fixture["commits"],
        fixture["files"],
        labels or {label["name"] for label in pr["labels"]},
        pr["body"],
        regions_of(fixture),
        fixture.get("head_changelog"),
        parity,
        adm.declared_models(fixture["commits"]),
        base_index=index_text(fixture, "base_index_file"),
        head_index=index_text(fixture, "head_index_file"),
        divergences_close=not fixture.get("opened_before_divergences_rule", False),
    )


def evaluate(
    fixture: dict, labels: set[str] | None = None, parity: adm.Finding | None = None
) -> set[str]:
    return {finding.check for finding in findings_of(fixture, labels, parity)}


def write_tree(root: Path, rows: str = "", tests: dict[str, str] | None = None) -> None:
    """A tree the parity checker reads: the matrix, an args.rs without flags, a plugins page whose
    support table has no rows, the named test files."""
    (root / "tests" / "parity").mkdir(parents=True, exist_ok=True)
    (root / "tests" / "parity" / "matrix.tsv").write_text(MATRIX_HEADER + rows)
    (root / "crates" / "vep-cli" / "src").mkdir(parents=True, exist_ok=True)
    (root / "crates" / "vep-cli" / "src" / "args.rs").write_text("pub struct Args {}\n")
    (root / "docs").mkdir(exist_ok=True)
    (root / "docs" / "plugins.md").write_text(PLUGINS_PAGE)
    for path, text in (tests or {}).items():
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        (root / path).write_text(text)


ROW = "foo\tConfig.pm\terror\tlanded:error\tcrates/vep-cli/src/x.rs::refuses_foo\t\t\n"
TEST_FILE = "#[cfg(test)]\nmod tests {\n    #[test]\n    fn refuses_foo() {}\n}\n"


# --- the Rust scanner ----------------------------------------------------------------------------


def test_regions_cover_a_cfg_test_module_and_a_test_fn_outside_it():
    text = (
        "fn code() {}\n"  # 1
        "\n"
        "#[test]\n"  # 3
        "fn alone() {\n"
        "    assert!(true);\n"
        "}\n"  # 6
        "\n"
        "#[cfg(test)]\n"  # 8
        "mod tests {\n"
        "    use super::*;\n"
        "    #[test]\n"
        "    fn inner() {}\n"
        "}\n"  # 13
        "fn more_code() {}\n"
    )
    assert adm.rust_test_regions(text) == [(3, 6), (8, 13)]


@pytest.mark.parametrize(
    "attr",
    [
        "#[tokio::test]",
        "#[rstest]",
        "#[test_log::test(tokio::test)]",
        "#[test_case(1, 2)]",
        "#[cfg(all(test, unix))]",
    ],
)
def test_regions_recognise_the_test_attribute_shapes(attr):
    text = f"{attr}\nfn f() {{\n}}\nfn g() {{}}\n"
    assert adm.rust_test_regions(text) == [(1, 3)]


@pytest.mark.parametrize(
    "attr",
    [
        "#[cfg(not(test))]",
        "#[cfg_attr(test, derive(Debug))]",
        "#[derive(Debug)]",
        "#[should_panic]",
    ],
)
def test_regions_leave_non_test_attributes_alone(attr):
    text = f"{attr}\nfn f() {{\n}}\n"
    assert adm.rust_test_regions(text) == []


def test_regions_survive_braces_in_strings_chars_comments_and_lifetimes():
    text = (
        "fn code<'a>(s: &'a str) -> &'a str {\n"  # 1
        '    let _ = "{ not a brace }";\n'
        '    let _ = r#"{ " }"#;\n'
        "    let _ = '{';\n"
        "    let _ = '\\'';\n"
        "    // } not a brace\n"
        "    /* } /* nested } */ */\n"
        "    s\n"
        "}\n"  # 9
        "#[cfg(test)]\n"  # 10
        "mod tests {\n"
        "    fn t<'b>() -> &'b str { \"}\" }\n"
        "}\n"  # 13
    )
    assert adm.rust_test_regions(text) == [(10, 13)]


def test_regions_keep_a_pending_attribute_across_commas_in_a_signature():
    text = "#[cfg(test)]\nfn helper(a: &str, b: HashMap<String, u8>) -> Vec<(u8, u8)> {\n    vec![]\n}\nfn code() {}\n"
    assert adm.rust_test_regions(text) == [(1, 4)]


@pytest.mark.parametrize(
    "attr, is_test",
    [
        ("#[cfg(test)]", True),
        ("#[cfg(all(test, feature = \"x\"))]", True),
        ("#[cfg(all(unix, test))]", True),
        ("#[cfg(any(test, not(test)))]", False),
        ("#[cfg(not(all(test, feature = \"never\")))]", False),
        ("#[cfg(feature = \"test-helpers\")]", False),
        ("#[cfg(not(test))]", False),
        ("#[cfg(all(any(test, unix)))]", False),
    ],
)
def test_regions_only_a_cfg_true_in_a_test_build_marks_a_region(attr, is_test):
    text = f"pub fn keep() {{}}\n{attr}\npub fn item() {{}}\n"
    assert adm.rust_test_regions(text) == ([(2, 3)] if is_test else [])
    assert adm.cfg_is_test(attr[2:-1]) is is_test


def test_cfg_parser_bounds_its_depth_and_survives_hostile_predicates():
    deep = "cfg(" + "all(" * 500 + "test" + ")" * 500 + ")"
    assert adm.cfg_is_test(deep) is False
    nested_ok = "cfg(" + "all(" * 3 + "test" + ")" * 3 + ")"
    assert adm.cfg_is_test(nested_ok) is True
    at_cap = "cfg(" + "all(" * adm.CFG_DEPTH_MAX + "test" + ")" * adm.CFG_DEPTH_MAX + ")"
    assert adm.cfg_is_test(at_cap) is True
    past_cap = "cfg(" + "all(" * (adm.CFG_DEPTH_MAX + 1) + "test" + ")" * (adm.CFG_DEPTH_MAX + 1) + ")"
    assert adm.cfg_is_test(past_cap) is False
    for hostile in ("cfg(all(test", "cfg(all(feature = \"a,b\", test))", "cfg(all(\"unterminated, test))", "cfg()"):
        adm.cfg_is_test(hostile)
    assert adm.cfg_is_test("cfg(all(feature = \"a,b\", test))") is True


def test_regions_single_line_item_and_inner_attribute():
    assert adm.rust_test_regions("#[cfg(test)]\nuse std::fmt;\nfn code() {}\n") == [(1, 2)]
    assert adm.rust_test_regions("#![cfg(test)]\nfn a() {}\nfn b() {}\n") == [(1, 3)]


def test_regions_field_attribute_does_not_open_a_region_on_the_next_item():
    text = (
        "struct S {\n    #[cfg(test)]\n    x: u8,\n    y: u8,\n}\nfn code() {\n    let z = 1;\n}\n"
    )
    assert adm.rust_test_regions(text) == [(2, 3)]


def test_regions_cover_a_field_or_variant_attribute_alone():
    """A `#[cfg(test)]` on a struct field, a struct literal's field or an enum variant marks that
    field or variant alone: through its comma, or its last token when no comma follows."""
    text = (
        "struct Lookups {\n"  # 1
        "    rules: Rules,\n"
        "    #[cfg(test)]\n"  # 3
        "    rows: Vec<CachedVariation>,\n"  # 4
        "    #[cfg(test)]\n"  # 5
        "    pub kept: &'static str\n"  # 6
        "}\n"
        "impl Lookups {\n"
        "    fn open() -> Result<Lookups, Error> {\n"
        "        Ok(Lookups {\n"
        "            rules: Rules::default(),\n"
        "            #[cfg(test)]\n"  # 12
        "            rows: Vec::new(),\n"  # 13
        "            #[cfg(test)]\n"  # 14
        '            kept: "",\n'  # 15
        "        })\n"
        "    }\n"
        "}\n"
        "enum Kind {\n"
        "    #[cfg(test)]\n"  # 20
        "    Probe(fn(u8), [u8; 2]),\n"  # 21
        "    #[cfg(test)]\n"  # 22
        "    Pair {\n"
        "        a: u8,\n"
        "    },\n"  # 25
        "    Real,\n"
        "}\n"
        "fn code() {}\n"
    )
    assert adm.rust_test_regions(text) == [(3, 4), (5, 6), (12, 13), (14, 15), (20, 21), (22, 25)]


def test_regions_field_attribute_spans_generics_arguments_and_trailing_comments():
    """A comma inside a field's generic brackets does not end it; an attribute inside a call's
    parentheses ends with the argument; a comment after the last field is not its last token, and a
    string that is ends the field on the string's last line."""
    text = (
        "struct S {\n"  # 1
        "    #[cfg(test)]\n"  # 2
        "    map: HashMap<\n"
        "        String,\n"
        "        Vec<u8>,\n"
        "    >,\n"  # 6
        "    #[cfg(test)]\n"  # 7
        "    last: Box<dyn Fn(u8) -> u8>\n"  # 8
        "    // the body's end\n"
        "}\n"
        "fn f() {\n"
        "    foo(x, #[cfg(test)] b);\n"  # 12
        "    if c {\n"
        "        real();\n"
        "    }\n"
        "    let fast = a < b;\n"
        "}\n"
    )
    assert adm.rust_test_regions(text) == [(2, 6), (7, 8), (12, 12)]
    text = 'fn f() {\n    L {\n        #[cfg(test)]\n        kept: "multi\nline"\n    }\n}\n'
    assert adm.rust_test_regions(text) == [(3, 5)]
    text = "struct S {\n    #[cfg(test)]\n    m: Mapö<\n        String,\n        u8,\n    >,\n}\n"
    assert adm.rust_test_regions(text) == [(2, 6)]


def test_regions_item_words_inside_parentheses_or_as_names_open_no_item():
    """A tuple field or variant element whose type opens with an item word (`fn(u8)`, `*const u8`,
    `&'static str`) and a field, variant or binding spelled with a weak keyword (`union`) are
    targets that end at their comma, so the attribute never swallows the next one; an item inside a
    closure body passed as an argument is still an item."""
    text = "struct S(u8, #[cfg(test)] fn(u8));\n"
    assert adm.rust_test_regions(text) == [(1, 1)]
    text = "struct S(#[cfg(test)] fn(u8), #[cfg(test)] *const u8,);\nfn code() {\n}\n"
    assert adm.rust_test_regions(text) == [(1, 1), (1, 1)]
    text = "enum E {\n    Probe(\n        u8,\n        #[cfg(test)] &'static str,\n    ),\n}\n"
    assert adm.rust_test_regions(text) == [(4, 4)]
    text = "struct S {\n    #[cfg(test)]\n    union: u8,\n    y: u8,\n}\n"
    assert adm.rust_test_regions(text) == [(2, 3)]
    text = "let s = S {\n    #[cfg(test)]\n    union,\n    #[cfg(test)]\n    macro_rules,\n};\n"
    assert adm.rust_test_regions(text) == [(2, 3), (4, 5)]
    text = "match k {\n    #[cfg(test)]\n    union => 1,\n    _ => 2,\n}\n"
    assert adm.rust_test_regions(text) == [(2, 3)]
    text = "#[cfg(test)]\nunion U {\n    a: u8,\n}\n#[cfg(test)]\nmacro_rules! m {\n    () => {};\n}\n"
    assert adm.rust_test_regions(text) == [(1, 4), (5, 8)]
    text = "#[cfg(test)]\nunion Ü<T>\nwhere\n    T: Copy,\n{\n    a: T,\n}\n"
    assert adm.rust_test_regions(text) == [(1, 7)]
    text = (
        "fn f() {\n"
        "    spawn(move || {\n"
        "        #[cfg(test)]\n"  # 3
        "        fn g<T>()\n"
        "        where\n"
        "            T: A,\n"
        "        {\n"
        "        }\n"  # 8
        "        real();\n"
        "    });\n"
        "}\n"
    )
    assert adm.rust_test_regions(text) == [(3, 8)]


def test_regions_hold_every_test_attribute_of_the_repository_sources():
    """Every `#[cfg(test)]` and `#[test]` line of the crates lies in a region, and a region starts at an
    attribute line."""
    checked = 0
    for path in sorted(ROOT.glob("crates/*/src/**/*.rs")):
        text = path.read_text(encoding="utf-8", errors="replace")
        regions = adm.rust_test_regions(text)
        lines = text.splitlines()
        for number, line in enumerate(lines, 1):
            if line.lstrip().startswith(("#[cfg(test)]", "#[test]")):
                assert adm.in_regions(number, regions), f"{path}:{number}"
                checked += 1
        for start, _end in regions:
            assert lines[start - 1].lstrip().startswith("#["), f"{path}:{start}"
    assert checked > 100


# --- the diff ------------------------------------------------------------------------------------


def test_changed_lines_number_removed_lines_in_the_old_file_and_added_in_the_new():
    patch = "@@ -10,3 +10,4 @@ fn x() {\n a\n-b\n+c\n+d\n e\n@@ -0,0 +1,2 @@\n+only\n+new\n\\ No newline at end of file\n"
    removed, added = adm.changed_lines(patch)
    assert removed == [(11, "b")]
    assert added == [(11, "c"), (12, "d"), (1, "only"), (2, "new")]


@pytest.mark.parametrize(
    "path",
    [
        "crates/vep-cli/tests/golden.rs",
        "crates/vep-cli/tests/fixtures/regulatory/README.md",
        "tools/vep-fields-compare/tests/compare.rs",
        "scripts/parity/test_check_parity_matrix.py",
        "scripts/concordance/test_run_clone_measurement.sh",
        "tests/golden/116/GRCh38/manifest.json",
        "tests/golden/116/GRCh38/expected/default.txt",
        "tests/golden/116/GRCh38/expected_hgvs/vcf.txt",
    ],
)
def test_test_paths(path):
    assert adm.is_test_path(path)


@pytest.mark.parametrize(
    "path",
    [
        "crates/vep-cli/src/args.rs",
        "scripts/parity/check_parity_matrix.py",
        "scripts/data/setup.sh",
        "scripts/conftest.py",
        "docs/plugins.md",
        "tests/parity/matrix.tsv",
        "tests/golden/116/GRCh38/variants.vcf",
        "tests/golden/116/GRCh38/provenance.json",
        "tests/fixtures/plugins/scores.txt",
    ],
)
def test_non_test_paths(path):
    assert not adm.is_test_path(path)


def rust_entry(patch: str, **extra) -> dict:
    return {
        "filename": "crates/vep-core/src/lib.rs",
        "status": "modified",
        "additions": 1,
        "deletions": 1,
        "patch": patch,
        **extra,
    }


def test_classify_rust_lines_by_region():
    regions = lambda side, path: [(20, 40)]
    only_code = rust_entry("@@ -5,2 +5,2 @@\n-a\n+b\n")
    only_test = rust_entry("@@ -25,2 +25,2 @@\n-a\n+b\n")
    both = rust_entry("@@ -5,1 +5,1 @@\n-a\n+b\n@@ -30,1 +30,1 @@\n-c\n+d\n")
    assert adm.classify_file(only_code, regions) == (True, False)
    assert adm.classify_file(only_test, regions) == (False, True)
    assert adm.classify_file(both, regions) == (True, True)


def test_classify_removed_test_lines_are_neither_code_nor_a_test_change():
    regions = lambda side, path: [(1, 10)] if side == "old" else []
    removed_test = rust_entry("@@ -3,2 +3,0 @@\n-a\n-b\n", deletions=2, additions=0)
    assert adm.classify_file(removed_test, regions) == (False, False)
    removed_code = rust_entry("@@ -30,2 +30,0 @@\n-a\n-b\n", deletions=2, additions=0)
    assert adm.classify_file(removed_code, regions) == (True, False)


def test_classify_without_a_patch_or_with_a_rename_only():
    no_patch = {
        "filename": "crates/vep-core/src/lib.rs",
        "status": "modified",
        "additions": 400,
        "deletions": 2,
    }
    assert adm.classify_file(no_patch, lambda side, path: []) == (True, False)
    rename = {
        "filename": "crates/vep-core/src/b.rs",
        "previous_filename": "crates/vep-core/src/a.rs",
        "status": "renamed",
        "additions": 0,
        "deletions": 0,
    }
    assert adm.classify_file(rename, lambda side, path: []) == (False, False)


def test_classify_a_file_moved_out_of_a_test_directory_is_code_and_no_test_change():
    none = lambda side, path: []
    moved = {
        "filename": "crates/vep-core/src/b.rs",
        "previous_filename": "crates/vep-core/tests/a.rs",
        "status": "renamed",
        "additions": 2,
        "deletions": 0,
        "patch": "@@ -1,0 +1,2 @@\n+x\n+y\n",
    }
    assert adm.classify_file(moved, none) == (True, False)
    assert adm.check_tests_included([moved], set(), none) is not None
    moved_no_patch = dict(moved)
    del moved_no_patch["patch"]
    assert adm.classify_file(moved_no_patch, none) == (True, False)
    into_tests = {
        "filename": "crates/vep-core/tests/a.rs",
        "previous_filename": "crates/vep-core/src/a.rs",
        "status": "renamed",
        "additions": 0,
        "deletions": 0,
    }
    assert adm.classify_file(into_tests, none) == (False, True)


BASE_RS = "// license\n//! doc\n\npub fn keep() -> u32 { 1 }\n"


def text_entry(path: str, old: str, new: str, status: str = "modified") -> dict:
    """An API file entry built from two texts, the patch as the API renders it (hunks only)."""
    import difflib

    lines = list(difflib.unified_diff(old.splitlines(), new.splitlines(), lineterm="", n=3))
    patch = "\n".join(lines[2:])
    return {
        "filename": path,
        "status": status,
        "additions": sum(1 for line in patch.splitlines() if line.startswith("+")),
        "deletions": sum(1 for line in patch.splitlines() if line.startswith("-")),
        "patch": patch,
    }


def regions_from(texts: dict[tuple[str, str], str]) -> adm.RegionsFor:
    return lambda side, path: adm.rust_test_regions(texts.get((side, path), "")) if texts.get(
        (side, path)
    ) else []


def test_tests_included_a_deleted_test_or_a_removed_test_fn_alone_does_not_count():
    code_path = "crates/vep-core/src/lib.rs"
    code_new = BASE_RS + "\npub fn behaviour_change() -> u32 { 3 }\n"
    deleted_test = {
        "filename": "crates/vep-core/tests/behaviour.rs",
        "status": "removed",
        "additions": 0,
        "deletions": 2,
        "patch": "@@ -1,2 +0,0 @@\n-#[test]\n-fn pins_behaviour() { assert_eq!(1, 1); }",
    }
    regions = regions_from({("new", code_path): code_new, ("old", code_path): BASE_RS})
    files = [text_entry(code_path, BASE_RS, code_new), deleted_test]
    assert adm.check_tests_included(files, set(), regions) is not None
    src_old = BASE_RS + "\n#[cfg(test)]\nmod tests {\n    #[test]\n    fn a() {}\n    #[test]\n    fn b() {}\n}\n"
    src_new = BASE_RS + "\npub fn behaviour_change() -> u32 { 3 }\n\n#[cfg(test)]\nmod tests {\n    #[test]\n    fn a() {}\n}\n"
    regions = regions_from({("new", code_path): src_new, ("old", code_path): src_old})
    assert adm.check_tests_included([text_entry(code_path, src_old, src_new)], set(), regions) is not None
    src_added = src_old.replace("    fn b() {}\n", "    fn b() {}\n    #[test]\n    fn c() {}\n")
    src_added = src_added.replace("pub fn keep() -> u32 { 1 }", "pub fn keep() -> u32 { 2 }")
    regions = regions_from({("new", code_path): src_added, ("old", code_path): src_old})
    assert adm.check_tests_included([text_entry(code_path, src_old, src_added)], set(), regions) is None


def test_tests_included_production_code_under_a_release_cfg_needs_a_test():
    path = "crates/vep-core/src/lib.rs"
    for attr in (
        "#[cfg(any(test, not(test)))]",
        "#[cfg(not(all(test, feature = \"never\")))]",
        "#[cfg(feature = \"test-helpers\")]",
    ):
        new = BASE_RS + f"\n{attr}\npub fn production_change() -> u32 {{ 2 }}\n"
        regions = regions_from({("new", path): new, ("old", path): BASE_RS})
        files = [text_entry(path, BASE_RS, new)]
        assert adm.classify_file(files[0], regions) == (True, False), attr
        assert adm.check_tests_included(files, set(), regions) is not None, attr
        assert adm.is_crate_change(files, regions) is True, attr


def test_benches_and_examples_are_crate_code_with_no_test_asked():
    none = lambda side, path: []
    bench = {"filename": "crates/vep-core/benches/speed.rs", "status": "modified", "additions": 3, "deletions": 1}
    example = {"filename": "crates/vep-cli/examples/run.rs", "status": "added", "additions": 9, "deletions": 0}
    assert adm.check_tests_included([bench, example], set(), none) is None
    assert adm.is_crate_change([bench], none) is True
    bench_test = {
        "filename": "crates/vep-core/benches/speed.rs",
        "status": "modified",
        "additions": 2,
        "deletions": 0,
        "patch": "@@ -1,0 +1,2 @@\n+#[test]\n+fn t() {}\n",
    }
    everywhere = lambda side, path: [(1, 2)]
    assert adm.classify_file(bench_test, everywhere) == (True, False)
    code = {"filename": "crates/vep-core/src/lib.rs", "status": "modified", "additions": 1, "deletions": 0}
    assert adm.check_tests_included([code, bench_test], set(), everywhere) is not None
    fixture = {"filename": "scripts/ci/fixtures/admission/shape.json", "status": "modified", "additions": 1, "deletions": 1}
    testdata = {"filename": "scripts/concordance/testdata/reference_release/expected_115_2/a.txt", "status": "modified", "additions": 1, "deletions": 0}
    assert not adm.is_test_path(fixture["filename"]) and not adm.is_test_path(testdata["filename"])
    script = {"filename": "scripts/golden/build.py", "status": "modified", "additions": 1, "deletions": 0}
    assert adm.check_tests_included([script, fixture, testdata], set(), none) is not None


def test_workflow_files_are_a_miss_for_edits_renames_and_removals():
    ci = text_entry(
        ".github/workflows/ci.yml",
        "jobs:\n  test:\n    name: Test\n    steps:\n      - run: cargo test\n",
        "jobs:\n  test:\n    name: Test\n    steps:\n      - run: 'true'\n",
    )
    finding = adm.check_workflow_files([ci])
    assert finding is not None and finding.closes and "`.github/workflows/ci.yml`" in finding.text
    assert "maintainers'" in finding.text
    renamed_out = {
        "filename": "docs/ci.yml",
        "previous_filename": ".github/workflows/ci.yml",
        "status": "renamed",
        "additions": 0,
        "deletions": 0,
    }
    assert adm.check_workflow_files([renamed_out]) is not None
    renamed_in = {
        "filename": ".github/workflows/extra.yml",
        "previous_filename": "docs/extra.yml",
        "status": "renamed",
        "additions": 0,
        "deletions": 0,
    }
    assert adm.check_workflow_files([renamed_in]) is not None
    removed = {"filename": ".github/workflows/dco.yml", "status": "removed", "additions": 0, "deletions": 20}
    assert adm.check_workflow_files([removed]) is not None
    docs = {"filename": "docs/plugins.md", "status": "modified", "additions": 1, "deletions": 0}
    assert adm.check_workflow_files([docs]) is None


LIB_RS = "pub fn f() {}\n\n#[cfg(test)]\nmod helpers;\n\nmod real;\n"


def test_tree_regions_cover_a_module_file_declared_under_cfg_test():
    texts = {
        "crates/x/src/lib.rs": LIB_RS,
        "crates/x/src/helpers.rs": "pub fn make() -> u8 {\n    1\n}\n",
        "crates/x/src/real.rs": "pub fn g() {}\n#[test]\nfn t() {}\n",
        "crates/x/src/sub/mod.rs": "#[cfg(test)]\npub(crate) mod fixtures;\n",
        "crates/x/src/sub/fixtures.rs": "pub const A: u8 = 1;\npub const B: u8 = 2;\n",
    }
    read = texts.get
    assert adm.tree_regions(read, "crates/x/src/helpers.rs") == [(1, 3)]
    assert adm.tree_regions(read, "crates/x/src/real.rs") == [(2, 3)]
    assert adm.tree_regions(read, "crates/x/src/sub/fixtures.rs") == [(1, 2)]
    assert adm.tree_regions(read, "crates/x/src/missing.rs") == []
    assert adm.module_parents("crates/x/src/sub/mod.rs") == (
        "sub",
        ["crates/x/src/mod.rs", "crates/x/src/lib.rs", "crates/x/src/main.rs", "crates/x/src.rs"],
    )


def test_regions_variant_attribute_does_not_extend_to_the_next_variant():
    text = "pub enum E {\n    #[cfg(test)]\n    Probe,\n    Real {\n        x: u8,\n    },\n}\n"
    assert adm.rust_test_regions(text) == [(2, 3)]
    text = "struct S {\n    #[cfg(test)]\n    x: HashMap<String, u8>,\n    y: u8,\n}\nfn f() {\n}\n"
    assert adm.rust_test_regions(text) == [(2, 3)]
    text = "pub enum E {\n    #[cfg(test)]\n    Probe(&'static str),\n    Real {\n        x: u8,\n    },\n}\n"
    assert adm.rust_test_regions(text) == [(2, 3)]
    text = "pub enum E {\n    #[cfg(test)]\n    Probe(fn(u8), [u8; 2]),\n    Real {\n        x: u8,\n    },\n}\n"
    assert adm.rust_test_regions(text) == [(2, 3)]
    text = "#[test]\nfn f(a: (u8, u8), b: [u8; 2]) -> Vec<(u8, u8)>\nwhere\n    T: Into<u8>,\n{\n}\nfn g() {}\n"
    assert adm.rust_test_regions(text) == [(1, 6)]


def test_classify_scripts_docs_and_test_paths():
    none = lambda side, path: []
    entry = {
        "filename": "scripts/golden/build.py",
        "status": "modified",
        "additions": 1,
        "deletions": 0,
    }
    assert adm.classify_file(entry, none) == (True, False)
    entry = {"filename": "scripts/data/stage.sh", "status": "added", "additions": 9, "deletions": 0}
    assert adm.classify_file(entry, none) == (True, False)
    entry = {
        "filename": "scripts/golden/test_build.py",
        "status": "modified",
        "additions": 1,
        "deletions": 0,
    }
    assert adm.classify_file(entry, none) == (False, True)
    entry = {"filename": "docs/plugins.md", "status": "modified", "additions": 1, "deletions": 0}
    assert adm.classify_file(entry, none) == (False, False)
    entry = {
        "filename": "crates/vep-cli/tests/cli.rs",
        "status": "modified",
        "additions": 1,
        "deletions": 0,
    }
    assert adm.classify_file(entry, none) == (False, True)


def test_tests_included_label_exempts():
    files = [
        {
            "filename": "scripts/golden/build.py",
            "status": "modified",
            "additions": 1,
            "deletions": 0,
        }
    ]
    none = lambda side, path: []
    assert adm.check_tests_included(files, set(), none) is not None
    assert adm.check_tests_included(files, {adm.NO_TESTS_LABEL}, none) is None
    assert adm.check_tests_included([], set(), none) is None


# --- the commits ---------------------------------------------------------------------------------


def commit(sha: str, message: str, email: str = "a@example.com", parents: int = 1) -> dict:
    return {
        "sha": sha,
        "parents": [{"sha": "p"}] * parents,
        "commit": {"message": message, "author": {"name": "A", "email": email}},
    }


def test_signoff_needs_the_author_name_and_address_and_skips_merges():
    signed = commit("1" * 40, "Subject\n\nSigned-off-by: A <a@example.com>\n")
    other_address = commit(
        "2" * 40,
        "Subject\n\nSigned-off-by: Bot <support@example.com>\n",
        email="1+bot[bot]@users.noreply.example.com",
    )
    unsigned = commit("3" * 40, "Subject\n")
    other_name = commit("5" * 40, "Subject\n\nSigned-off-by: Someone Else <a@example.com>\n")
    merge = commit("4" * 40, "Merge main into topic\n", parents=2)
    assert adm.check_signoff([signed, merge]) is None
    finding = adm.check_signoff([signed, other_address, unsigned, other_name])
    assert finding is not None and "`2222222`" in finding.text and "`3333333`" in finding.text
    assert "`5555555`" in finding.text and "`1111111`" not in finding.text
    assert "author's name and address" in finding.text
    assert finding.closes is False and "new pull request from a signed branch" in finding.text


def test_signoff_lists_at_most_five_commits():
    finding = adm.check_signoff([commit(f"{n}" * 40, "Subject\n") for n in range(1, 9)])
    assert finding is not None and "and 3 more" in finding.text


def test_declaration_names_models_and_normalises_no_llm():
    commits = [
        commit(
            "1" * 40,
            "Subject\n\nCo-Authored-By: Example Model 2 <model@example.com>\nSigned-off-by: A <a@example.com>\n",
        ),
        commit("2" * 40, f"Subject\n\nco-authored-by:   {adm.NO_MODEL_DECLARATION.lower()}  \n"),
        commit("5" * 40, f"Subject\n\nCo-Authored-By: {'  '.join(adm.NO_MODEL_DECLARATION.split())}\n"),
        commit("3" * 40, "Subject\n\nCo-Authored-By: Example Model 2 <other@example.com>\n"),
    ]
    assert adm.declared_models(commits) == ["Example Model 2", "No LLM"]
    assert adm.check_declaration(adm.declared_models(commits)) is None
    signed_only = [commit("4" * 40, "Subject\n\nSigned-off-by: A <a@example.com>\n")]
    missing = adm.check_declaration(adm.declared_models(signed_only))
    assert missing is not None and missing.closes is True
    assert "one new signed commit" in missing.text and "do not amend" in missing.text


def test_label_models_come_from_the_newest_declaring_commit_and_stop_at_two_names():
    head = commit("9" * 40, "Subject\n\nCo-Authored-By: Example Model 1 <m@example.com>\n")
    earlier = commit("8" * 40, "Subject\n\nCo-Authored-By: Example Model 2 <m@example.com>\n")
    assert adm.label_models([earlier, head]) == ["Example Model 1"]
    assert adm.label_models([head, earlier]) == ["Example Model 2"]
    assert adm.label_models([]) == []
    assert adm.label_models([commit("7" * 40, "Subject\n")]) == []
    merge_at_head = commit("4" * 40, "Merge main into topic\n", parents=2)
    assert adm.label_models([earlier, head, merge_at_head]) == ["Example Model 1"]
    two = commit(
        "6" * 40,
        "Subject\n\nCo-Authored-By: Example Model 1 <m@example.com>\nCo-Authored-By: No LLM\n",
    )
    assert adm.label_models([two]) == ["Example Model 1", "No LLM"]
    many = [
        commit(f"{n:040x}", f"c\n\nCo-Authored-By: Jane Doe {n} <j{n}@example.com>\n")
        for n in range(1, 5)
    ]
    assert adm.label_models(many) is None
    assert adm.declared_models(many) == [f"Jane Doe {n}" for n in range(1, 5)]


def test_model_label_is_bounded_and_clean():
    assert adm.model_label("Example Model 2") == "model: Example Model 2"
    assert len(adm.model_label("x" * 80)) == adm.LABEL_MAX
    assert adm.model_label("a\x00b\nc") == "model: abc"
    assert (
        adm.model_label("@someone <!-- /admission:static -->")
        == "model: someone -- /admissionstatic --"
    )
    assert adm.model_label("<>") is None


# --- the description and the changelog -----------------------------------------------------------


TEMPLATE_BODY = (
    "## Summary\n<!-- One or two sentences. -->\nAdds a flag.\n\n"
    "## Verification\n<!-- What was run. -->\ncargo test passes.\n\n"
    "## Upgrade impact\n<!-- None, or the change. -->\nNone.\n\n- [x] Every commit is signed off\n"
)


def test_body_sections_present_and_filled():
    assert adm.check_body_sections(TEMPLATE_BODY) is None
    assert (
        adm.check_body_sections("### summary\nx\n## VERIFICATION\ny\n## Upgrade Impact\nz") is None
    )
    bullets = "## Summary\n- adds a flag\n## Verification\n* ran the tests\n## Upgrade impact\n1. none\n"
    assert adm.check_body_sections(bullets) is None
    fenced = "```\n## Summary\nx\n## Verification\nx\n## Upgrade impact\nx\n```\n"
    finding = adm.check_body_sections(fenced)
    assert finding is not None and all(f"`## {name}`" in finding.text for name in adm.SECTIONS)
    output_under_heading = (
        "## Summary\nAdds a flag.\n## Verification\n```\n## not a heading\ncargo test: ok\n```\n"
        "## Upgrade impact\nNone.\n"
    )
    assert adm.check_body_sections(output_under_heading) is None
    four_backticks = (
        "## Summary\nAdds a flag.\n## Verification\n````text\n```\n## not a heading\n````\n"
        "## Upgrade impact\nNone.\n"
    )
    assert adm.check_body_sections(four_backticks) is None
    unterminated = "## Summary\nx\n## Verification\n```\ny\n## Upgrade impact\nz\n"
    finding = adm.check_body_sections(unterminated)
    assert finding is not None and "`## Upgrade impact`" in finding.text


def test_body_sections_on_the_repository_template():
    """The template's checklist sits under its last heading and is not that section's text: the
    unfilled template misses all three sections, and filling each heading's line clears them."""
    template = (ROOT / ".github" / "pull_request_template.md").read_text(encoding="utf-8")
    unfilled = adm.check_body_sections(template)
    assert unfilled is not None
    assert all(f"`## {name}`" in unfilled.text for name in adm.SECTIONS)
    filled = template
    for name in adm.SECTIONS:
        filled = filled.replace(f"## {name}\n", f"## {name}\nFilled in.\n", 1)
    assert adm.check_body_sections(filled) is None
    checklist_only = "## Summary\n- [x] Every commit is signed off\n## Verification\ny\n## Upgrade impact\nz\n"
    finding = adm.check_body_sections(checklist_only)
    assert finding is not None and "`## Summary`" in finding.text


@pytest.mark.parametrize(
    "body, missing",
    [
        ("## Summary\nx\n## Verification\ny\n", "`## Upgrade impact`"),
        (
            "## Summary\n<!-- only the template comment -->\n## Verification\ny\n## Upgrade impact\nz\n",
            "`## Summary`",
        ),
        ("Bumps a dependency.\n", "`## Summary`, `## Verification`, `## Upgrade impact`"),
        (None, "`## Summary`"),
    ],
)
def test_body_sections_missing_or_empty(body, missing):
    finding = adm.check_body_sections(body)
    assert finding is not None and missing in finding.text


CHANGELOG = "# Changelog\n\n## [Unreleased]\n\n### Added\n\n- Entry (#0)\n\n## [0.1.0] - 0000-00-00\n\n### Added\n\n- Old (#0)\n"


def test_unreleased_range():
    assert adm.unreleased_range(CHANGELOG) == (3, 8)
    assert adm.unreleased_range("# Changelog\n\n## [0.1.0] - 0000-00-00\n") is None
    assert adm.unreleased_range("## [Unreleased]\n\n- a\n") == (1, 3)


NONE = lambda side, path: []


def crate_files(changelog_patch: str | None, status: str = "modified") -> list[dict]:
    files = [
        {
            "filename": "crates/vep-core/src/lib.rs",
            "status": "modified",
            "additions": 1,
            "deletions": 0,
        }
    ]
    if changelog_patch is not None:
        files.append(
            {
                "filename": "CHANGELOG.md",
                "status": status,
                "additions": 1,
                "deletions": 0,
                "patch": changelog_patch,
            }
        )
    return files


def test_changelog_line_under_unreleased_satisfies_a_crate_change():
    assert (
        adm.check_changelog(crate_files("@@ -6,2 +6,3 @@\n \n+- Entry (#0)\n \n"), CHANGELOG, NONE)
        is None
    )


@pytest.mark.parametrize(
    "patch",
    [None, "@@ -12,1 +12,2 @@\n- Old (#0)\n+- Newer old (#0)\n", "@@ -6,1 +6,2 @@\n \n+\n"],
)
def test_changelog_missing_outside_unreleased_or_blank_is_a_miss(patch):
    assert adm.check_changelog(crate_files(patch), CHANGELOG, NONE) is not None


def test_changelog_without_hunks_passes_only_when_lines_changed():
    big = {"filename": "CHANGELOG.md", "status": "modified", "additions": 400, "deletions": 0}
    assert adm.check_changelog(crate_files(None) + [big], CHANGELOG, NONE) is None
    mode_only = {"filename": "CHANGELOG.md", "status": "modified", "additions": 0, "deletions": 0}
    assert adm.check_changelog(crate_files(None) + [mode_only], CHANGELOG, NONE) is not None


def test_changelog_not_required_without_a_crate_code_change():
    files = [
        {
            "filename": "crates/vep-cli/tests/cli.rs",
            "status": "modified",
            "additions": 3,
            "deletions": 0,
        },
        {"filename": "docs/plugins.md", "status": "modified", "additions": 1, "deletions": 0},
        {
            "filename": "crates/vep-core/README.md",
            "status": "modified",
            "additions": 3,
            "deletions": 0,
        },
    ]
    assert adm.check_changelog(files, CHANGELOG, NONE) is None
    assert adm.check_changelog(crate_files(None), None, NONE) is not None
    no_unreleased = "# No unreleased\n"
    assert (
        adm.check_changelog(
            crate_files("@@ -1,1 +1,2 @@\n x\n+- Entry (#0)\n"), no_unreleased, NONE
        )
        is not None
    )
    manifest = [
        {
            "filename": "crates/vep-core/Cargo.toml",
            "status": "modified",
            "additions": 1,
            "deletions": 1,
        }
    ]
    assert adm.check_changelog(manifest, CHANGELOG, NONE) is not None


RELEASE_CHANGELOG = (
    "# Changelog\n\n## [Unreleased]\n\n## [0.2.0] - 0000-00-00\n\n### Added\n\n- Entry (#0)\n\n"
    "## [0.1.0] - 0000-00-00\n\n### Added\n\n- Old (#0)\n"
)


def test_changelog_release_heading_inserted_above_existing_entries_satisfies_it():
    head = (
        "# Changelog\n\n## [Unreleased]\n\n## [0.2.0] - 0000-00-00\n\n### Added\n\n- Entry (#0)\n"
    )
    patch = "@@ -3,3 +3,5 @@\n ## [Unreleased]\n \n+## [0.2.0] - 0000-00-00\n+\n ### Added\n"
    manifest = {
        "filename": "crates/vep-core/Cargo.toml",
        "status": "modified",
        "additions": 1,
        "deletions": 1,
    }
    changelog = {
        "filename": "CHANGELOG.md",
        "status": "modified",
        "additions": 2,
        "deletions": 0,
        "patch": patch,
    }
    assert adm.check_changelog([manifest, changelog], head, NONE) is None
    empty = "# Changelog\n\n## [Unreleased]\n\n## [0.2.0] - 0000-00-00\n"
    assert adm.check_changelog([manifest, changelog], empty, NONE) is not None


def test_changelog_line_under_a_release_heading_the_change_adds_satisfies_it():
    patch = (
        "@@ -3,6 +3,6 @@\n ## [Unreleased]\n \n-### Added\n-\n-- Entry (#0)\n"
        "+## [0.2.0] - 0000-00-00\n+\n+### Added\n+\n+- Entry (#0)\n \n"
    )
    manifest = {
        "filename": "crates/vep-core/Cargo.toml",
        "status": "modified",
        "additions": 1,
        "deletions": 1,
    }
    changelog = {
        "filename": "CHANGELOG.md",
        "status": "modified",
        "additions": 5,
        "deletions": 3,
        "patch": patch,
    }
    assert adm.check_changelog([manifest, changelog], RELEASE_CHANGELOG, NONE) is None
    heading_only = "@@ -3,2 +3,4 @@\n ## [Unreleased]\n \n+## [0.2.0] - 0000-00-00\n+\n"
    changelog = {
        "filename": "CHANGELOG.md",
        "status": "modified",
        "additions": 2,
        "deletions": 0,
        "patch": heading_only,
    }
    bare = "# Changelog\n\n## [Unreleased]\n\n## [0.2.0] - 0000-00-00\n\n## [0.1.0] - 0000-00-00\n\n- Old (#0)\n"
    assert adm.check_changelog([manifest, changelog], bare, NONE) is not None


# --- the intended divergences --------------------------------------------------------------------


SECTION_HEAD = f"## {adm.DIVERGENCE_SECTION}\n"
STATEMENTS = {
    "none": adm.NONE_STATEMENT,
    "touches": "This change touches these existing entries, one line each:",
    "proposes": "This change proposes a new intended divergence, described below:",
}


def section(checked: set[str], lines: dict[str, list[str]] | None = None) -> str:
    """A divergences section with the named statements checked and the given lines under each."""
    out = [SECTION_HEAD.rstrip("\n")]
    for kind, text in STATEMENTS.items():
        out.append(f"- [{'x' if kind in checked else ' '}] {text}")
        out.extend(f"  {line}" for line in (lines or {}).get(kind, []))
    return "\n".join(out) + "\n"


def divergences(body: str, files: list[dict] | None = None, head: str | None = None, **kw):
    return adm.check_divergences(TEMPLATE_BODY + body, files or [], INDEX_SAMPLE, head, **kw)


def test_index_names_come_from_rows_that_link_a_page_on_both_index_shapes():
    """A defect letter, a difference's name as the first cell writes it, an identifier of the
    regrouped index, and every link target with and without the `docs/` prefix and any fragment;
    not a header, a number, the literal `none`, a row without a page link, or a page linked outside
    the tables."""
    index = adm.read_index(INDEX_SAMPLE)
    for token in (
        "B",
        "a",
        "JSON escaping",
        "A bare M",
        "A chromosome the cache lacks",
        "divergences/defect-b.md",
        "docs/divergences/defect-b.md",
        "cache-setup.md#chromosome-names",
        "cache-setup.md",
        "docs/cache-setup.md",
        "plugins/spliceai.md#reverse-strand-indels",
        "Docs/Divergences/Defect-B.md",
    ):
        assert adm.entry_name(token) in index.names, token
    for token in ("none", "1", "3", "Defect", "Difference", "Defect A", "divergences/reproducing.md", ""):
        assert adm.entry_name(token) not in index.names, token
    assert index.pages == {
        "docs/divergences/defect-a.md",
        "docs/divergences/defect-b.md",
        "docs/divergences/cnv-tr.md",
        "docs/divergences/defect-m.md",
        "docs/divergences/json-escaping.md",
        "docs/divergences/bare-m-chromosome.md",
        "docs/cache-setup.md",
        "docs/plugins/spliceai.md",
    }
    regrouped = adm.read_index(INDEX_SAMPLE_REGROUPED)
    for token in ("start_cooccurrence_swap", "spliceai_reverse_strand_indel", "divergences/defect-c.md"):
        assert adm.entry_name(token) in regrouped.names, token
    assert adm.entry_name("B") not in regrouped.names and adm.entry_name("3") not in regrouped.names
    both = index | regrouped
    assert adm.entry_name("B") in both.names and adm.entry_name("start_cooccurrence_swap") in both.names
    assert adm.read_index(None) == adm.Index(frozenset(), frozenset())


def test_index_names_on_the_repository_index():
    """Every entry of the repository's own index resolves by its first cell and by its page."""
    index = adm.read_index((ROOT / adm.INDEX_PATH).read_text(encoding="utf-8"))
    assert len(index.names) > 50 and len(index.pages) > 20
    assert adm.entry_name("none") not in index.names
    assert all(not name.isdigit() for name in index.names)
    assert adm.INDEX_PATH not in index.pages
    assert all(page.startswith("docs/") and page.endswith(".md") for page in index.pages)


def test_table_cells_split_on_unescaped_pipes_and_entry_name_normalises():
    assert adm.table_cells("| a | `b \\| c` | [d](e.md) |") == ["a", "`b | c`", "[d](e.md)"]
    assert adm.table_cells("a | b") == ["a", "b"]
    assert adm.is_separator_row("| --- | :---: |") and not adm.is_separator_row("| a | b |")
    assert adm.entry_name("[A bare `M`](divergences/bare-m-chromosome.md)") == "a bare m"
    assert adm.entry_name("  JSON   escaping ") == "json escaping"
    assert adm.entry_name("[`start_cooccurrence_swap`](divergences/defect-b.md)") == "start_cooccurrence_swap"


def test_statements_are_read_by_their_words_with_their_lines():
    lines = [
        "intro text",
        "- [x] This change touches none of the intended divergences.",
        "* [ ] It touches these existing entries:",
        "  - `B`: changes it, because the reasoning assumes a start codon.",
        "  plain continuation",
        "- [X] A new intended divergence is proposed below:",
        "  The description.",
        "- [ ] An unrelated checklist line",
        "  not attached to a statement",
    ]
    statements = adm.read_statements(lines)
    assert set(statements) == {"none", "touches", "proposes"}
    assert statements["none"] == (True, [])
    assert statements["touches"][0] is False
    assert adm.entry_lines(statements["touches"][1]) == [
        ("B", "changes it, because the reasoning assumes a start codon.")
    ]
    assert statements["proposes"] == (True, ["  The description."])
    assert adm.entry_lines(["`docs/divergences/defect-b.md` - removed", "1. `A`: kept", "text", "- `C`"]) == [
        ("docs/divergences/defect-b.md", "removed"),
        ("A", "kept"),
        ("C", ""),
    ]


def test_divergences_none_passes_and_is_refused_on_an_artefact_edit():
    assert divergences(section({"none"})) is None
    for path in (
        adm.INDEX_PATH,
        "docs/divergences/defect-b.md",
        "docs/cache-setup.md",
        "docs/plugins/spliceai.md",
        "scripts/divergences/reproductions/x.toml",
    ):
        finding = divergences(section({"none"}), [{"filename": path, "status": "modified"}])
        assert finding is not None and finding.check == "divergences" and finding.closes
        assert f"`{path}`" in finding.text and "under the second statement" in finding.text, path
    renamed = {"filename": "docs/other.md", "previous_filename": "docs/divergences/defect-b.md"}
    assert divergences(section({"none"}), [renamed]) is not None
    for path in (
        "docs/divergences.md",
        "docs/plugins.md",
        "scripts/concordance/compare_vep_outputs.py",
        "tests/golden/116/GRCh38/manifest.json",
        "crates/vep-effects/src/lib.rs",
    ):
        assert divergences(section({"none"}), [{"filename": path, "status": "modified"}]) is None, path


def test_divergences_exclusivity_and_missing_shapes():
    missing = divergences("")
    assert missing is not None and "lacks the `## Intended divergences` section" in missing.text
    assert divergences(SECTION_HEAD + "text without statements\n") is not None
    nothing = divergences(section(set()))
    assert nothing is not None and "No statement" in nothing.text
    for other in ("touches", "proposes"):
        both = divergences(section({"none", other}, {"touches": ["- `B`: kept, because x."], "proposes": ["desc"]}))
        assert both is not None and "excludes the other two" in both.text
    assert adm.check_divergences(None, [], None, None) is not None
    soft = divergences(section(set()), closes=False)
    assert soft is not None and soft.closes is False


def test_divergences_entries_resolve_against_base_and_head_with_a_rationale():
    ok = {"touches": ["- `B`: changes it, because the index's reasoning is wrong for this record."]}
    assert divergences(section({"touches"}, ok)) is None
    by_name = {"touches": ["- `JSON escaping`: removes it, because the writer changed."]}
    assert divergences(section({"touches"}, by_name)) is None
    by_path = {"touches": ["`docs/divergences/defect-b.md`: keeps it, because the class holds."]}
    assert divergences(section({"touches"}, by_path)) is None
    by_fragment = {"touches": ["- `cache-setup.md#chromosome-names`: changes it, because the synonym rule moved."]}
    assert divergences(section({"touches"}, by_fragment)) is None
    head_only = {"touches": ["- `start_cooccurrence_swap`: changes it, because the shape misses a record."]}
    assert divergences(section({"touches"}, head_only)) is not None
    assert divergences(section({"touches"}, head_only), head=INDEX_SAMPLE_REGROUPED) is None
    unknown = divergences(section({"touches"}, {"touches": ["- `Q`: removed, because x.", "- `1`: kept, because y.", "- `none`: kept, because z."]}))
    assert unknown is not None and "`Q`, `1`, `none`" in unknown.text and "first column" in unknown.text
    bare = divergences(section({"touches"}, {"touches": ["- `B`", "- `A`:"]}))
    assert bare is not None and "`B`, `A`" in bare.text and "says nothing after the entry" in bare.text
    none_listed = divergences(section({"touches"}, {"touches": ["Changed defect B as described in the summary."]}))
    assert none_listed is not None and "no entry line follows it" in none_listed.text
    many = divergences(section({"touches"}, {"touches": [f"- `Z{n}`: x, because y." for n in range(7)]}))
    assert many is not None and "and 2 more" in many.text


def test_divergences_proposal_needs_a_description_and_may_join_the_entries():
    assert divergences(section({"proposes"}, {"proposes": ["The input, both rows, the reason."]})) is None
    empty = divergences(section({"proposes"}, {"proposes": ["", "   "]}))
    assert empty is not None and "no description follows it" in empty.text
    both = {"touches": ["- `B`: removes it, because x."], "proposes": ["A new one, described."]}
    assert divergences(section({"touches", "proposes"}, both)) is None
    half = {"touches": ["- `B`: removes it, because x."]}
    assert divergences(section({"touches", "proposes"}, half)) is not None


def test_divergences_strips_html_comments_and_fenced_headings():
    commented = section({"touches"}) + "<!-- - `B`: an example line in a comment. -->\n"
    finding = divergences(commented)
    assert finding is not None and "no entry line follows it" in finding.text
    hidden = "```\n" + section({"none"}) + "```\n"
    assert divergences(hidden) is not None
    in_comment = "<!--\n" + section({"none"}) + "-->\n" + section(set())
    assert divergences(in_comment) is not None


def test_fenced_lines_under_a_statement_are_neither_statements_nor_entries():
    """A template fragment quoted in a fenced block under the section checks nothing; an
    entry-shaped line inside a fence names no entry; fenced text still counts as a description."""
    quoted_template = section(set()) + "```\n- [x] This change touches none of the intended divergences.\n```\n"
    finding = divergences(quoted_template)
    assert finding is not None and "No statement" in finding.text
    fenced_entry = {"touches": ["- `B`: changes it, because x.", "```", "- `Q`: not a citation", "```"]}
    assert divergences(section({"touches"}, fenced_entry)) is None
    only_fenced = {"touches": ["```", "- `B`: changes it", "```"]}
    finding = divergences(section({"touches"}, only_fenced))
    assert finding is not None and "no entry line follows it" in finding.text
    fenced_description = {"proposes": ["```", "the record and both rows", "```"]}
    assert divergences(section({"proposes"}, fenced_description)) is None
    assert adm.unfenced(["a", "```", "b", "```", "c", "~~~", "d"]) == ["a", "c"]
    assert adm.fenced_lines(["```x`", "- [ ] y"]) == [("```x`", False), ("- [ ] y", False)]


def test_a_line_opening_with_an_inline_code_span_is_not_a_fence():
    """A backtick run followed by text holding another backtick is an inline code span, so the
    headings after it still open sections and fields; a real fence still hides them."""
    assert adm.opens_fence("```") == "```" and adm.opens_fence("~~~ text") == "~~~"
    assert adm.opens_fence("````rust") == "````" and adm.opens_fence("``` cargo test `x`") is None
    assert adm.opens_fence("```ERROR: [cache] shard missing```") is None
    assert adm.opens_fence("text") is None
    body = (
        "## Summary\nAdds a flag.\n## Verification\n```cargo test -p vep-cli```\n"
        "## Upgrade impact\nNone.\n" + section({"none"})
    )
    assert adm.check_body_sections(body) is None
    assert adm.check_divergences(body, [], INDEX_SAMPLE, None) is None
    issue = (
        "### What happened\n\n```ERROR: [cache] shard missing```\nand then it stopped\n\n"
        f"### {adm.ISSUE_CHOICE_LABEL}\n\n{adm.ISSUE_CHOICES['none']}\n\n"
        f"### {adm.ISSUE_ENTRIES_LABEL}\n\n{adm.NO_RESPONSE}\n"
    )
    fields = adm.form_fields(issue)
    assert fields[adm.ISSUE_CHOICE_LABEL.lower()] == adm.ISSUE_CHOICES["none"]
    assert fields["what happened"].startswith("```ERROR") and fields["what happened"].endswith("stopped")
    assert adm.check_issue_divergences(issue, adm.read_index(INDEX_SAMPLE)) is None


def test_an_entry_whose_name_holds_a_backtick_is_cited_in_a_longer_run_or_without_it():
    index = adm.read_index(INDEX_SAMPLE)
    assert adm.entry_name("A bare M") in index.names
    assert adm.entry_lines(["- ``A bare `M` ``: keeps it, because x.", "- `A bare M`: keeps it."]) == [
        ("A bare `M`", "keeps it, because x."),
        ("A bare M", "keeps it."),
    ]
    for line in ("- ``A bare `M` ``: keeps it, because x.", "- `A bare M`: keeps it, because x."):
        assert divergences(section({"touches"}, {"touches": [line]})) is None, line
    assert adm.entry_lines(["- `B`: `x` and `y`"]) == [("B", "`x` and `y`")]
    assert adm.entry_lines(["- ``: nothing", "- ` `: nothing", "- ```A`: x"]) == []
    # A fence line, bare or with an info string, is no entry line, so a fenced block quoted under
    # the second statement leaves a valid citation valid.
    assert adm.entry_lines(["```", "```text", "````", "~~~"]) == []
    fenced = {"touches": ["- `B`: changes it, because x.", "```text", "the two rows, quoted", "```"]}
    assert divergences(section({"touches"}, fenced)) is None
    issue = issue_body(adm.ISSUE_CHOICES["touches"], "`B`: should go, because x.\n```\nrows\n```")
    assert adm.check_issue_divergences(issue, index) is None
    # A description is at most 65,536 characters; a line of backticks, or of spaces before a last
    # character, is read in bounded time.
    hostile = "- " + "`" * 65_536
    started = time.perf_counter()
    adm.entry_lines(
        [hostile, hostile + "x", hostile + "x" + "`" * 65_535 + "x", "- `a`: x" + " " * 65_536 + "x"]
    )
    assert time.perf_counter() - started < 5


def test_entry_remedies_name_the_change_for_a_pull_request_and_the_report_for_an_issue():
    index = adm.read_index(INDEX_SAMPLE)
    bare = divergences(section({"touches"}, {"touches": ["- `B`"]}))
    assert bare is not None and f"add {adm.CHANGE_REMEDY} and why" in bare.text
    none_listed = divergences(section({"touches"}))
    assert none_listed is not None and f"`<entry>`: {adm.CHANGE_REMEDY}, because ...``" in none_listed.text
    bare_issue = adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["touches"], "`B`"), index)
    assert bare_issue is not None and f"add {adm.REPORT_REMEDY} and why" in bare_issue
    assert adm.CHANGE_REMEDY not in bare_issue
    none_issue = adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["touches"]), index)
    assert none_issue is not None and f"`<entry>`: {adm.REPORT_REMEDY}, because ...``" in none_issue


@RELEASE_TREE
def test_divergences_section_on_the_repository_template():
    """The template's fourth section carries the three statements, none checked, with its example
    entry line inside a comment; checking the first statement alone passes."""
    template = TEMPLATE.read_text(encoding="utf-8")
    lines = adm.body_sections(template)[adm.DIVERGENCE_SECTION.lower()]
    statements = adm.read_statements(lines)
    assert set(statements) == {"none", "touches", "proposes"}
    assert not any(checked for checked, _lines in statements.values())
    assert adm.entry_lines(lines) == []
    unfilled = adm.check_divergences(template, [], INDEX_SAMPLE, None)
    assert unfilled is not None and "No statement" in unfilled.text
    none_line = next(line for line in lines if adm.CHECKBOX_RE.match(line) and "none" in line)
    assert adm.check_divergences(template.replace(none_line, none_line.replace("[ ]", "[x]"), 1), [], INDEX_SAMPLE, None) is None
    assert adm.check_body_sections(template) is not None
    assert adm.SECTIONS == ("Summary", "Verification", "Upgrade impact")


@RELEASE_TREE
def test_issue_forms_carry_the_divergence_choice():
    """Each of the ten forms asks the choice as a required dropdown with the three statements,
    then the entries field, after its last required field; blank issues stay off."""
    forms = sorted(path for path in FORMS_DIR.glob("*.yml") if path.name != "config.yml")
    assert len(forms) == 10
    bodies = {form: yaml.safe_load(form.read_text(encoding="utf-8"))["body"] for form in forms}
    wordings: set[tuple[str, str]] = set()
    for form, fields in bodies.items():
        ids = [field.get("id") for field in fields]
        at = ids.index("divergences")
        dropdown, textarea = fields[at], fields[at + 1]
        assert dropdown["type"] == "dropdown", form.name
        assert dropdown["attributes"]["label"] == adm.ISSUE_CHOICE_LABEL, form.name
        assert dropdown["attributes"]["options"] == list(adm.ISSUE_CHOICES.values()), form.name
        assert dropdown["validations"]["required"] is True, form.name
        assert textarea["type"] == "textarea" and textarea["id"] == "divergence-entries", form.name
        assert textarea["attributes"]["label"] == adm.ISSUE_ENTRIES_LABEL, form.name
        assert "render" not in textarea["attributes"], form.name
        assert not any(
            (field.get("validations") or {}).get("required") for field in fields[at + 1 :]
        ), form.name
        wordings.add((dropdown["attributes"]["description"], textarea["attributes"]["description"]))
    assert len(wordings) == 1
    config = yaml.safe_load((FORMS_DIR / "config.yml").read_text(encoding="utf-8"))
    assert config["blank_issues_enabled"] is False


def issue_body(choice: str | None, entries: str = adm.NO_RESPONSE) -> str:
    head = "### Version\n\nvep 0.0.0\n\n"
    if choice is None:
        return head + "### Command\n\nvep\n"
    return (
        f"{head}### {adm.ISSUE_CHOICE_LABEL}\n\n{choice}\n\n"
        f"### {adm.ISSUE_ENTRIES_LABEL}\n\n{entries}\n"
    )


def test_form_fields_and_the_issue_choice():
    fields = adm.form_fields(issue_body(adm.ISSUE_CHOICES["none"]))
    assert fields["version"] == "vep 0.0.0"
    assert fields[adm.ISSUE_CHOICE_LABEL.lower()] == adm.ISSUE_CHOICES["none"]
    assert fields[adm.ISSUE_ENTRIES_LABEL.lower()] == ""
    assert adm.form_fields(None) == {}
    quoted = issue_body(
        adm.ISSUE_CHOICES["none"],
        f"The form I filled in:\n```\n### {adm.ISSUE_CHOICE_LABEL}\n\nsomething else\n```\nend",
    )
    fields = adm.form_fields(quoted)
    assert fields[adm.ISSUE_CHOICE_LABEL.lower()] == adm.ISSUE_CHOICES["none"]
    assert fields[adm.ISSUE_ENTRIES_LABEL.lower()].startswith("The form I filled in:\n```")
    assert fields[adm.ISSUE_ENTRIES_LABEL.lower()].endswith("```\nend")
    index = adm.read_index(INDEX_SAMPLE)
    assert adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["none"]), index) is None
    assert adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["none"], "stray text"), index) is None
    assert (
        adm.check_issue_divergences(
            issue_body(adm.ISSUE_CHOICES["touches"], "`B`: should go, because the reasoning is wrong here.\n`JSON escaping`: keep"),
            index,
        )
        is None
    )
    unknown = adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["touches"], "`Q`: remove it"), index)
    assert unknown is not None and "`Q`" in unknown
    empty = adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["touches"]), index)
    assert empty is not None and "names existing entries and no entry line follows it" in empty
    assert adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["proposes"], "A record, both rows, the reason."), index) is None
    missing = adm.check_issue_divergences(issue_body(adm.ISSUE_CHOICES["proposes"]), index)
    assert missing is not None and "field is empty" in missing
    absent = adm.check_issue_divergences(issue_body(None), index)
    assert absent is not None and f"no `{adm.ISSUE_CHOICE_LABEL}` choice" in absent
    assert adm.check_issue_divergences("", index) is not None
    reworded = adm.check_issue_divergences(issue_body("Touches none of them"), index)
    assert reworded is not None and "no `Intended divergences` choice" in reworded


def issue_event(body: str | None, labels=(), state="open", **extra) -> dict:
    issue = {
        "number": 7,
        "state": state,
        "body": body,
        "labels": [{"name": name} for name in labels],
        "user": {"login": "someone", "type": "User"},
    }
    issue.update(extra)
    return {"action": "opened", "issue": issue}


def issue_repo(tmp_path: Path) -> Path:
    base = repo_dir(tmp_path)
    (base / "docs").mkdir(exist_ok=True)
    (base / adm.INDEX_PATH).write_text(INDEX_SAMPLE, encoding="utf-8")
    return base


def test_run_issue_labels_and_comments_a_miss_and_never_closes(tmp_path, capsys):
    gh = FakeGitHub()
    code = adm.run_issue(issue_event(issue_body(adm.ISSUE_CHOICES["touches"], "`Q`: drop it")), gh, issue_repo(tmp_path))
    assert code == 0
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("POST", f"/repos/{REPO}/issues/7/comments"),
        ("POST", f"/repos/{REPO}/issues/7/labels"),
    ]
    assert writes[0][2]["body"].startswith(adm.ISSUE_MARKER) and "`Q`" in writes[0][2]["body"]
    assert "Edit the issue and the check runs again." in writes[0][2]["body"]
    assert writes[1][2] == {"labels": [adm.NEEDS_DIVERGENCE_LABEL]}
    assert all("state" not in (w[2] if len(w) > 2 else {}) for w in writes)
    assert "ERROR: [admission]" in capsys.readouterr().err
    gh = FakeGitHub()
    assert adm.run_issue(issue_event(issue_body(None)), gh, issue_repo(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [
        ("POST", f"/repos/{REPO}/issues/7/comments"),
        ("POST", f"/repos/{REPO}/issues/7/labels"),
    ]


def test_run_issue_edits_its_comment_in_place_and_clears_a_fixed_issue(tmp_path):
    existing = {"id": 3, "body": f"{adm.ISSUE_MARKER}\nold", "user": {"login": "github-actions[bot]", "type": "Bot"}}
    gh = FakeGitHub(comments=[existing])
    body = issue_body(adm.ISSUE_CHOICES["proposes"])
    assert adm.run_issue(issue_event(body, labels=[adm.NEEDS_DIVERGENCE_LABEL]), gh, issue_repo(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [("PATCH", f"/repos/{REPO}/issues/comments/3")]
    assert "field is empty" in gh.writes()[0][2]["body"]
    gh = FakeGitHub(comments=[existing])
    fixed = issue_body(adm.ISSUE_CHOICES["proposes"], "The record, both rows and the reason.")
    assert adm.run_issue(issue_event(fixed, labels=[adm.NEEDS_DIVERGENCE_LABEL]), gh, issue_repo(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [
        ("PATCH", f"/repos/{REPO}/issues/comments/3"),
        ("DELETE", f"/repos/{REPO}/issues/7/labels/needs-divergence-check"),
    ]
    assert "passes" in gh.writes()[0][2]["body"]


def test_run_issue_writes_nothing_for_a_clean_closed_or_pull_request_issue(tmp_path):
    clean = issue_body(adm.ISSUE_CHOICES["none"])
    gh = FakeGitHub()
    assert adm.run_issue(issue_event(clean), gh, issue_repo(tmp_path)) == 0 and gh.writes() == []
    gh = FakeGitHub()
    assert adm.run_issue(issue_event(issue_body(None), state="closed"), gh, issue_repo(tmp_path)) == 0
    assert gh.writes() == []
    gh = FakeGitHub()
    assert adm.run_issue(issue_event(issue_body(None), pull_request={"url": "x"}), gh, issue_repo(tmp_path)) == 0
    assert gh.calls == []
    gh = FakeGitHub()
    without_index = repo_dir(tmp_path / "without-index")
    assert adm.run_issue(issue_event(clean), gh, without_index) == 0 and gh.writes() == []


def test_run_issue_acts_on_the_labels_and_state_as_re_read(tmp_path):
    """A run queued behind another sees a stale payload: the label it removes or adds, and the
    state it respects, are the issue's as re-read before the writes."""
    fixed = issue_body(adm.ISSUE_CHOICES["proposes"], "The record, both rows and the reason.")
    live = {"number": 7, "state": "open", "labels": [{"name": adm.NEEDS_DIVERGENCE_LABEL}]}
    existing = {"id": 3, "body": f"{adm.ISSUE_MARKER}\nold", "user": {"login": "github-actions[bot]", "type": "Bot"}}
    gh = FakeGitHub(comments=[existing], issue=live)
    assert adm.run_issue(issue_event(fixed, labels=[]), gh, issue_repo(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [
        ("PATCH", f"/repos/{REPO}/issues/comments/3"),
        ("DELETE", f"/repos/{REPO}/issues/7/labels/needs-divergence-check"),
    ]
    gh = FakeGitHub(issue={"number": 7, "state": "open", "labels": [{"name": adm.NEEDS_DIVERGENCE_LABEL}]})
    assert adm.run_issue(issue_event(issue_body(None), labels=[]), gh, issue_repo(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [("POST", f"/repos/{REPO}/issues/7/comments")]
    gh = FakeGitHub(issue={"number": 7, "state": "closed", "labels": []})
    assert adm.run_issue(issue_event(issue_body(None)), gh, issue_repo(tmp_path)) == 0
    assert gh.writes() == []


def test_run_issue_skips_an_issue_opened_before_the_rule_landed(tmp_path, capsys):
    rule = [{"commit": {"committer": {"date": "1999-01-05T00:00:00Z"}}}]
    before = issue_event(issue_body(None), created_at="1999-01-04T00:00:00Z")
    gh = FakeGitHub(rule_landed=rule)
    assert adm.run_issue(dict(before, repository={"default_branch": "trunk"}), gh, issue_repo(tmp_path)) == 0
    assert gh.writes() == [] and "opened before the check" in capsys.readouterr().out
    landing_reads = [call for call in gh.calls if call[1].endswith(f"/repos/{REPO}/commits")]
    assert [call[2] for call in landing_reads] == [{"path": adm.DIVERGENCE_LANDING_PATH, "sha": "trunk"}]
    after = issue_event(issue_body(None), created_at="1999-01-06T00:00:00Z")
    gh = FakeGitHub(rule_landed=rule)
    assert adm.run_issue(after, gh, issue_repo(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [
        ("POST", f"/repos/{REPO}/issues/7/comments"),
        ("POST", f"/repos/{REPO}/issues/7/labels"),
    ]
    gh = FakeGitHub(fail_landed=True)
    assert adm.run_issue(after, gh, issue_repo(tmp_path)) == 0 and gh.writes() == []


# --- the parity matrix ---------------------------------------------------------------------------


def test_parity_resolves_on_a_consistent_tree(tmp_path):
    write_tree(tmp_path, ROW, {"crates/vep-cli/src/x.rs": TEST_FILE})
    assert adm.check_parity(tmp_path, CHECKER) is None


def test_parity_unresolved_cell_is_a_finding_naming_it(tmp_path):
    write_tree(
        tmp_path, ROW, {"crates/vep-cli/src/x.rs": TEST_FILE.replace("refuses_foo", "refuses_bar")}
    )
    finding = adm.check_parity(tmp_path, CHECKER)
    assert (
        finding is not None
        and "refuses_foo" in finding.text
        and str(adm.PARITY_CHECKER) in finding.text
    )


def test_parity_renders_contributor_text_inside_code_spans(tmp_path):
    row = ROW.replace("refuses_foo", "refuses_@someone_`x`")
    write_tree(tmp_path, row, {"crates/vep-cli/src/x.rs": TEST_FILE})
    finding = adm.check_parity(tmp_path, CHECKER)
    assert finding is not None and "`x`" not in finding.text and "@someone_'x'" in finding.text


def test_parity_checker_crash_is_a_closing_finding_naming_its_last_line(tmp_path):
    crashing = tmp_path / "checker.py"
    crashing.write_text("import sys\nprint('boom @x `y`', file=sys.stderr)\nsys.exit(3)\n")
    finding = adm.check_parity(tmp_path, crashing)
    assert finding is not None and finding.closes and finding.check == "parity"
    assert "could not read this change" in finding.text and "boom @x 'y'" in finding.text


def test_parity_checker_killed_silent_or_unimportable_is_an_environment_error(tmp_path):
    silent = tmp_path / "silent.py"
    silent.write_text("import sys\nsys.exit(1)\n")
    with pytest.raises(adm.AdmissionError, match="no output"):
        adm.check_parity(tmp_path, silent)
    killed = tmp_path / "killed.py"
    killed.write_text("import os, signal\nos.kill(os.getpid(), signal.SIGKILL)\n")
    with pytest.raises(adm.AdmissionError, match="exit -9"):
        adm.check_parity(tmp_path, killed)
    unimportable = tmp_path / "unimportable.py"
    unimportable.write_text("import module_that_does_not_exist_anywhere\n")
    with pytest.raises(adm.AdmissionError, match="ModuleNotFoundError"):
        adm.check_parity(tmp_path, unimportable)


def test_parity_checker_timeout_is_an_environment_error(tmp_path, monkeypatch):
    monkeypatch.setattr(adm, "PARITY_TIMEOUT_S", 1)
    hanging = tmp_path / "checker.py"
    hanging.write_text("import time\ntime.sleep(5)\n")
    with pytest.raises(adm.AdmissionError, match="did not finish"):
        adm.check_parity(tmp_path, hanging)


@pytest.mark.parametrize(
    "target, tail",
    [
        ("tests/parity/matrix.tsv", b"\xff\n"),
        ("crates/vep-cli/src/args.rs", b"\n// \xff\n"),
        ("tests/parity/matrix.tsv", b"x\t" + b"a" * 200_000 + b"\n"),
    ],
    ids=["invalid-utf8-matrix", "invalid-utf8-args", "oversized-cell"],
)
def test_parity_contributed_bytes_never_abort_the_run(tmp_path, target, tail):
    write_tree(tmp_path, ROW, {"crates/vep-cli/src/x.rs": TEST_FILE})
    path = tmp_path / target
    path.write_bytes(path.read_bytes() + tail)
    finding = adm.check_parity(tmp_path, CHECKER)
    assert finding is None or finding.check == "parity"


# --- CI ------------------------------------------------------------------------------------------


def job(name: str, conclusion: str) -> dict:
    return {
        "name": name,
        "conclusion": conclusion,
        "html_url": f"https://example.invalid/jobs/{name.replace(' ', '-')}",
    }


def test_check_ci_reports_failed_required_jobs_only():
    jobs = [
        job("Test", "failure"),
        job("Clippy", "success"),
        job("MSRV", "timed_out"),
        job("Build and check links", "failure"),
        job("Format", "cancelled"),
    ]
    findings = adm.check_ci(jobs, set(adm.REQUIRED_CHECKS))
    assert [f.check for f in findings] == ["ci", "ci"]
    assert (
        "`Test`" in findings[0].text
        and "jobs/Test" in findings[0].text
        and "`MSRV`" in findings[1].text
    )


class FakeGitHub:
    """The API surface the modes use, with canned answers and a record of every write."""

    def __init__(
        self,
        *,
        commits=(),
        files=(),
        comments=(),
        rules=None,
        pulls_by_head=(),
        pull=None,
        jobs=(),
        fail_rules=False,
        landed=(),
        rule_landed=None,
        fail_landed=False,
        issue=None,
    ):
        self.repo = REPO
        self.commits, self.files, self.comments = list(commits), list(files), list(comments)
        self.rules, self.pulls_by_head, self.pull, self.issue = rules, list(pulls_by_head), pull, issue
        self.jobs, self.fail_rules, self.landed = list(jobs), fail_rules, list(landed)
        # The commits that landed the divergences rule's fixture; the workflow's own when not given.
        self.rule_landed = list(landed if rule_landed is None else rule_landed)
        self.fail_landed = fail_landed
        self.calls: list[tuple] = []

    def paginate(self, path, params=None):
        self.calls.append(("GET", path) if params is None else ("GET", path, dict(params)))
        if path.endswith(f"/repos/{REPO}/commits"):
            if self.fail_landed:
                raise adm.AdmissionError("503")
            if (params or {}).get("path") == adm.DIVERGENCE_LANDING_PATH:
                return self.rule_landed
            return self.landed
        if path.endswith("/commits"):
            return self.commits
        if path.endswith("/files"):
            return self.files
        if path.endswith("/comments"):
            return self.comments
        if "/actions/runs/" in path:
            return self.jobs
        return []

    def get(self, path, params=None):
        self.calls.append(("GET", path))
        if "/rules/branches/" in path:
            if self.fail_rules:
                raise adm.AdmissionError("403")
            return self.rules
        if path.endswith("/pulls"):
            return self.pulls_by_head
        if "/pulls/" in path:
            return self.pull
        parts = path.split("/")
        if len(parts) > 2 and parts[-2] == "issues" and parts[-1].isdigit():
            return self.issue
        return None

    def post(self, path, body):
        self.calls.append(("POST", path, body))
        return {"id": 77, "body": body.get("body")}

    def patch(self, path, body):
        self.calls.append(("PATCH", path, body))
        return body

    def delete(self, path):
        self.calls.append(("DELETE", path))

    def writes(self):
        return [call for call in self.calls if call[0] != "GET"]


RULES = [
    {
        "type": "required_status_checks",
        "parameters": {"required_status_checks": [{"context": "Test"}, {"context": "Admission"}]},
    },
    {"type": "deletion"},
]


class PagedGitHub(adm.GitHub):
    """The real client over canned pages: a wrapped jobs page, then a bare list page."""

    def __init__(self):
        super().__init__(REPO, None)
        self.pages = [
            (
                200,
                {"Link": '<https://api.example/next>; rel="next"'},
                {"total_count": 3, "jobs": [{"name": "Test"}, {"name": "Clippy"}]},
            ),
            (200, {}, [{"name": "MSRV"}]),
        ]
        self.urls = []

    def request(self, method, path, body=None, params=None):
        self.urls.append((path, params))
        return self.pages.pop(0)


def test_paginate_unwraps_the_jobs_envelope_and_follows_the_next_link():
    gh = PagedGitHub()
    items = gh.paginate("/repos/x/y/actions/runs/1/jobs")
    assert [item["name"] for item in items] == ["Test", "Clippy", "MSRV"]
    assert gh.urls[0] == ("/repos/x/y/actions/runs/1/jobs", {"per_page": 100})
    assert gh.urls[1] == ("https://api.example/next", None)


def test_request_retries_a_read_once_and_a_write_never(monkeypatch):
    attempts = []

    def failing(req, timeout=0):
        attempts.append(req.get_method())
        raise urllib.error.HTTPError(req.full_url, 503, "unavailable", {}, None)

    monkeypatch.setattr(adm.urllib.request, "urlopen", failing)
    monkeypatch.setattr(adm.time, "sleep", lambda s: None)
    gh = adm.GitHub(REPO, "token")
    with pytest.raises(adm.AdmissionError):
        gh.get("/x")
    with pytest.raises(adm.AdmissionError):
        gh.post("/x", {"a": 1})
    assert attempts == ["GET", "GET", "POST"]


def test_required_checks_from_the_ruleset_with_the_fallback():
    assert adm.required_checks(FakeGitHub(rules=RULES), "main") == {"Test", "Admission"}
    assert adm.required_checks(FakeGitHub(rules=[]), "main") == set(adm.REQUIRED_CHECKS)
    assert adm.required_checks(FakeGitHub(fail_rules=True), "main") == set(adm.REQUIRED_CHECKS)


# --- the verdict ---------------------------------------------------------------------------------


def test_maintainers_from_the_repository_codeowners():
    logins = adm.maintainers_from_codeowners((ROOT / adm.CODEOWNERS).read_text())
    assert len(logins) == 2 and all(login == login.lower() for login in logins)


def test_maintainers_from_codeowners_shapes():
    text = "# comment\n*   @Abc @org/team docs@example.com\ndocs/ @beta\n"
    assert adm.maintainers_from_codeowners(text) == {"abc"}
    assert adm.maintainers_from_codeowners("docs/ @beta\n") == set()


FINDING = [adm.Finding("x", "x")]


def test_close_footer_exemptions():
    user = {"login": "someone", "type": "User"}
    maintainer = {"login": "First-Maintainer", "type": "User"}
    bot = {"login": "dependency-bot[bot]", "type": "Bot"}
    assert adm.is_bot({"login": "app/dependabot", "type": "User"}) and adm.is_bot(bot)
    assert not adm.is_bot(user)
    soft = [adm.Finding("signoff", "x", closes=False)]
    assert adm.close_footer([], user, set(), False, MAINTAINERS) == (False, None)
    assert adm.close_footer(FINDING, maintainer, set(), False, MAINTAINERS)[0] is False
    assert (
        adm.close_footer(FINDING, bot, set(), False, MAINTAINERS)[1]
        == "No automatic close: the author is a bot."
    )
    assert adm.close_footer(FINDING, user, {adm.OVERRIDE_LABEL}, False, MAINTAINERS)[0] is False
    assert (
        adm.close_footer(FINDING, user, set(), True, MAINTAINERS)[1]
        == "No automatic close: the pull request is a draft."
    )
    old_pr = adm.close_footer(FINDING, user, set(), False, MAINTAINERS, grandfathered=True)
    assert old_pr[0] is False and "before the admission check" in old_pr[1]
    softly = adm.close_footer(soft, user, set(), False, MAINTAINERS)
    assert softly[0] is False and "DCO check holds the merge" in softly[1]
    close, footer = adm.close_footer(FINDING + soft, user, set(), False, MAINTAINERS)
    assert close is True and footer is not None and "reopen" in footer


def test_predates_and_landed_at():
    assert adm.predates("1999-01-01T10:00:00Z", "1999-01-05T09:00:00Z") is True
    assert adm.predates("1999-01-05T09:00:01Z", "1999-01-05T09:00:00Z") is False
    assert adm.predates(None, "1999-01-05T09:00:00Z") is False
    assert adm.predates("1999-01-01T10:00:00Z", None) is False
    assert adm.landed_at(FakeGitHub(), "main") is None
    dates = ["1999-01-06T00:00:00Z", "1999-01-05T00:00:00Z"]
    gh = FakeGitHub(landed=[{"commit": {"committer": {"date": d}}} for d in dates])
    assert adm.landed_at(gh, "main") == "1999-01-05T00:00:00Z"


# --- the comment ---------------------------------------------------------------------------------


def test_code_span_strips_backticks_newlines_and_length():
    assert adm.code("a `b`\nc  @d") == "`a 'b' c @d`"
    assert adm.code("") == "`-`"
    assert len(adm.code("x" * 500)) == 122


def test_comment_blocks_round_trip_and_drop_a_stale_head():
    head = "a" * 40
    findings = [adm.Finding("base", "The branch does not contain the tip of `main`: merge it in.")]
    static = adm.static_block(findings, head, "https://example.invalid/policy")
    body = adm.render_comment(
        {"static": (head, static, True)}, head, "The pull request is closed until these are fixed."
    )
    assert body.startswith(adm.MARKER) and "closes=1" in body
    blocks = adm.parse_blocks(body)
    assert set(blocks) == {"static"} and blocks["static"][0] == head and blocks["static"][2] is True
    assert adm.block_has_findings(blocks["static"])
    ci = adm.ci_block(
        [adm.Finding("ci", "The CI job `Test` failed (u): fix it and push.")], "b" * 40
    )
    blocks["ci"] = ("b" * 40, ci, True)
    rendered = adm.render_comment(blocks, head, None)
    assert "After CI" not in rendered and "does not contain" in rendered
    passing = adm.render_comment(
        {"static": (head, adm.static_block([], head, "u"), False)}, head, None
    )
    assert "passes on head" in passing and not adm.block_has_findings(
        adm.parse_blocks(passing)["static"]
    )
    assert adm.parse_blocks(passing)["static"][2] is False


# --- the recorded pull requests ------------------------------------------------------------------


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_records_its_verdict(path):
    assert evaluate(load(path)) == RECORDED[path.stem]


def test_every_fixture_has_a_recorded_verdict():
    assert set(RECORDED) == set(FIXTURE_IDS) and len(FIXTURES) == FIXTURE_COUNT


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_divergences_miss_closes_only_after_the_rule_landed(path):
    """A pull request opened before the divergences rule carries the miss as comment-only; one
    opened after it is closed for the miss."""
    fixture = load(path)
    before = fixture.get("opened_before_divergences_rule", False)
    assert before is (not path.stem.startswith("divergences-"))
    divergences = [f for f in findings_of(fixture) if f.check == "divergences"]
    assert ("divergences" in RECORDED[path.stem]) is bool(divergences)
    for finding in divergences:
        assert finding.closes is (not before)


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_divergences_check_flips(path):
    """A description without the section misses; the same description with the none statement
    checked passes unless the change edits a divergence artefact, and then naming the entry passes."""
    fixture = load(path)
    pr = fixture["pull_request"]
    pr["body"] = (pr["body"] or "").split(f"## {adm.DIVERGENCE_SECTION}")[0]
    assert "divergences" in evaluate(fixture)
    pr["body"] += NONE_SECTION
    index = adm.read_index(index_text(fixture, "base_index_file")) | adm.read_index(
        index_text(fixture, "head_index_file")
    )
    artefacts = adm.divergence_artefacts(fixture["files"], index)
    assert ("divergences" in evaluate(fixture)) is bool(artefacts)
    if artefacts:
        pr["body"] = pr["body"].replace(
            f"- [x] {adm.NONE_STATEMENT}",
            "- [x] This change proposes a new intended divergence, described below:\n"
            "  The input, both engines' rows and why vep-rs's row is the one to keep.",
        )
        assert "divergences" not in evaluate(fixture)


def strip_trailers(fixture: dict, key: str) -> None:
    for entry in fixture["commits"]:
        message = entry["commit"]["message"]
        entry["commit"]["message"] = (
            "\n".join(line for line in message.splitlines() if not line.lower().startswith(key))
            + "\n"
        )


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_base_check_flips(path):
    fixture = load(path)
    recorded = "base" in RECORDED[path.stem]
    assert ("base" in evaluate(fixture)) is recorded
    fixture["contains_base_tip"] = not fixture["contains_base_tip"]
    assert ("base" in evaluate(fixture)) is not recorded


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_signoff_check_flips(path):
    fixture = load(path)
    if "signoff" in RECORDED[path.stem]:
        for entry in fixture["commits"]:
            author = entry["commit"]["author"]
            entry["commit"]["message"] += f"Signed-off-by: {author['name']} <{author['email']}>\n"
        assert "signoff" not in evaluate(fixture)
    else:
        strip_trailers(fixture, "signed-off-by:")
        assert "signoff" in evaluate(fixture)


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_declaration_check_flips(path):
    fixture = load(path)
    if "declaration" in RECORDED[path.stem]:
        fixture["commits"][-1]["commit"]["message"] += "Co-Authored-By: No LLM\n"
        assert "declaration" not in evaluate(fixture)
        assert adm.declared_models(fixture["commits"]) == ["No LLM"]
    else:
        strip_trailers(fixture, "co-authored-by:")
        assert "declaration" in evaluate(fixture)


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_body_check_flips(path):
    fixture = load(path)
    if "body" in RECORDED[path.stem]:
        fixture["pull_request"]["body"] = TEMPLATE_BODY
        assert "body" not in evaluate(fixture)
    else:
        fixture["pull_request"]["body"] = fixture["pull_request"]["body"].replace(
            "## Verification", "## Checks"
        )
        assert "body" in evaluate(fixture)


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_tests_check_flips_where_code_changes(path):
    fixture = load(path)
    none = lambda side, path: []
    changes_code = any(
        adm.classify_file(entry, regions_of(fixture))[0] for entry in fixture["files"]
    )
    assert "tests" not in RECORDED[path.stem]
    if not changes_code:
        return
    fixture["files"] = [
        entry for entry in fixture["files"] if not adm.is_test_path(entry["filename"])
    ]
    for entry in fixture["files"]:
        entry.pop("regions", None)
    findings = adm.static_findings(
        True,
        "main",
        fixture["commits"],
        fixture["files"],
        set(),
        TEMPLATE_BODY,
        none,
        None,
        None,
        ["No LLM"],
    )
    assert "tests" in {finding.check for finding in findings}
    labelled = adm.static_findings(
        True,
        "main",
        fixture["commits"],
        fixture["files"],
        {adm.NO_TESTS_LABEL},
        TEMPLATE_BODY,
        none,
        None,
        None,
        ["No LLM"],
    )
    assert "tests" not in {finding.check for finding in labelled}


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_changelog_check_flips_where_a_crate_changes(path):
    fixture = load(path)
    crate_change = adm.is_crate_change(fixture["files"], regions_of(fixture))
    if "changelog" in RECORDED[path.stem]:
        fixture["files"].append(
            {
                "filename": "CHANGELOG.md",
                "status": "modified",
                "additions": 1,
                "deletions": 0,
                "patch": "@@ -6,2 +6,3 @@\n \n+- Entry (#0)\n \n",
            }
        )
        fixture["head_changelog"] = CHANGELOG
        assert "changelog" not in evaluate(fixture)
    elif crate_change:
        fixture["files"] = [
            entry for entry in fixture["files"] if entry["filename"] != "CHANGELOG.md"
        ]
        assert "changelog" in evaluate(fixture)


@pytest.mark.parametrize("path", FIXTURES, ids=FIXTURE_IDS)
def test_fixture_close_decision(path):
    fixture = load(path)
    pr = fixture["pull_request"]
    findings = findings_of(fixture) or [adm.Finding("x", "x")]
    close, footer = adm.close_footer(findings, pr["user"], set(), pr["draft"], MAINTAINERS)
    contributor = pr["user"]["login"] not in MAINTAINERS and pr["user"]["type"] != "Bot"
    assert close is (contributor and any(f.closes for f in findings))
    if contributor and not close and all(f.check == "divergences" for f in findings):
        assert footer is not None and "opened before that section was asked" in footer


# --- the modes -----------------------------------------------------------------------------------


class FakeGit:
    def __init__(
        self,
        head: str,
        contains: bool = True,
        fetched: str | None = None,
        texts: dict | None = None,
        commits: list | None = None,
    ):
        self.head, self.contains, self.fetched, self.texts = (
            head,
            contains,
            fetched or head,
            texts or {},
        )
        self.commits = list(commits or [])

    def fetch_pull_head(self, number):
        return self.fetched

    def log_commits(self, base, head):
        return self.commits

    def rev_parse(self, rev):
        return "b" * 40

    def is_ancestor(self, ancestor, descendant):
        return self.contains

    def merge_base(self, a, b):
        return "m" * 40

    def show(self, rev, path):
        return self.texts.get(path)

    def export_tree(self, rev, dest):
        write_tree(dest, ROW, {"crates/vep-cli/src/x.rs": TEST_FILE})


def repo_dir(tmp_path: Path) -> Path:
    base = tmp_path / "base"
    (base / "scripts" / "parity").mkdir(parents=True, exist_ok=True)
    shutil.copy(CHECKER, base / adm.PARITY_CHECKER)
    (base / ".github").mkdir(exist_ok=True)
    (base / adm.CODEOWNERS).write_text("* @first-maintainer @second-maintainer\n")
    return base


def event_for(fixture: dict, action: str = "opened") -> dict:
    return {"action": action, "pull_request": fixture["pull_request"]}


def test_run_event_comments_labels_and_closes_a_contributor_miss(tmp_path, capsys):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    fixture["commits"][-1]["commit"]["message"] += (
        "Co-Authored-By: Example Model 3 <model@example.com>\n"
    )
    fixture["contains_base_tip"] = False
    gh = FakeGitHub(files=fixture["files"])
    git = FakeGit(fixture["pull_request"]["head"]["sha"], contains=False, commits=fixture["commits"])
    code = adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path))
    assert code == 1
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("POST", f"/repos/{REPO}/issues/0/comments"),
        ("PATCH", f"/repos/{REPO}/pulls/0"),
        ("POST", f"/repos/{REPO}/issues/0/labels"),
    ]
    comment = writes[0][2]["body"]
    assert comment.startswith(adm.MARKER) and "does not contain the tip of `main`" in comment
    assert "closed until" in comment and "declares the model" not in comment
    assert writes[1][2] == {"state": "closed"}
    assert writes[2][2] == {"labels": ["model: Example Model 3"]}
    assert "ERROR: [admission]" in capsys.readouterr().err


def test_run_event_edits_the_existing_comment_and_keeps_the_ci_block_of_the_same_head(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    head = fixture["pull_request"]["head"]["sha"]
    ci = adm.ci_block([adm.Finding("ci", "The CI job `Test` failed (u): fix it and push.")], head)
    old = adm.render_comment(
        {"ci": (head, ci, True), "static": ("0" * 40, "stale", True)}, head, None
    )
    existing = {"id": 5, "body": old, "user": {"login": "github-actions[bot]", "type": "Bot"}}
    gh = FakeGitHub(files=fixture["files"], comments=[existing])
    git = FakeGit(head, commits=fixture["commits"])
    code = adm.run_event(event_for(fixture, "synchronize"), gh, git, repo_dir(tmp_path))
    assert code == 1
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("PATCH", f"/repos/{REPO}/issues/comments/5"),
        ("PATCH", f"/repos/{REPO}/pulls/0"),
    ]
    body = writes[0][2]["body"]
    assert "After CI" in body and "stale" not in body and "declares the model" in body
    assert body.count(adm.MARKER) == 1


def test_sync_model_labels_skips_a_name_with_nothing_left():
    gh = FakeGitHub()
    adm.sync_model_labels(gh, 0, set(), ["<>"])
    assert gh.writes() == []
    adm.sync_model_labels(gh, 0, set(), ["<>", "Example Model 1"])
    assert gh.writes() == [
        ("POST", f"/repos/{REPO}/issues/0/labels", {"labels": ["model: Example Model 1"]})
    ]


def test_run_event_clean_pull_request_writes_nothing_but_the_label(tmp_path):
    fixture = load(FIXTURE_DIR / "maintainer-flag-bam.json")
    fixture["pull_request"]["body"] += NONE_SECTION
    gh = FakeGitHub(files=fixture["files"])
    git = FakeGit(
        fixture["pull_request"]["head"]["sha"],
        texts={"CHANGELOG.md": fixture["head_changelog"]},
        commits=fixture["commits"],
    )
    code = adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path))
    assert code == 0
    assert [w[:2] for w in gh.writes()] == [("POST", f"/repos/{REPO}/issues/0/labels")]


def test_run_event_maintainer_miss_is_commented_not_closed(tmp_path):
    fixture = load(FIXTURE_DIR / "maintainer-flag-bam.json")
    strip_trailers(fixture, "co-authored-by:")
    fixture["pull_request"]["labels"] = [{"name": "model: Example Model 1"}]
    gh = FakeGitHub(files=fixture["files"])
    git = FakeGit(
        fixture["pull_request"]["head"]["sha"],
        texts={"CHANGELOG.md": fixture["head_changelog"]},
        commits=fixture["commits"],
    )
    code = adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path))
    assert code == 0
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("POST", f"/repos/{REPO}/issues/0/comments"),
        ("DELETE", f"/repos/{REPO}/issues/0/labels/model%3A%20Example%20Model%201"),
    ]
    assert "No automatic close: the author is a maintainer." in writes[0][2]["body"]


def test_run_event_override_label_and_draft_are_commented_not_closed(tmp_path):
    codes = []
    for mutate in (
        lambda pr: pr["labels"].append({"name": adm.OVERRIDE_LABEL}),
        lambda pr: pr.update(draft=True),
    ):
        fixture = load(FIXTURE_DIR / "fork-performance-change.json")
        mutate(fixture["pull_request"])
        gh = FakeGitHub(files=fixture["files"])
        head = fixture["pull_request"]["head"]["sha"]
        git = FakeGit(head, commits=fixture["commits"])
        codes.append(adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)))
        assert all(w[1] != f"/repos/{REPO}/pulls/0" for w in gh.writes())
        assert "No automatic close" in gh.writes()[0][2]["body"]
    assert codes == [0, 0]


def test_run_event_signoff_miss_alone_is_commented_not_closed(tmp_path):
    fixture = load(FIXTURE_DIR / "maintainer-docs-with-declaration.json")
    fixture["pull_request"]["user"] = {"login": "someone", "type": "User"}
    fixture["pull_request"]["body"] += NONE_SECTION
    strip_trailers(fixture, "signed-off-by:")
    gh = FakeGitHub(files=fixture["files"])
    head = fixture["pull_request"]["head"]["sha"]
    git = FakeGit(head, commits=fixture["commits"])
    code = adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path))
    assert code == 0
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("POST", f"/repos/{REPO}/issues/0/comments"),
        ("POST", f"/repos/{REPO}/issues/0/labels"),
    ]
    body = writes[0][2]["body"]
    assert "carries no `Signed-off-by`" in body and "DCO check holds the merge" in body
    assert "closes=0" in body


def test_run_event_divergences_miss_closes_only_a_pull_request_opened_after_the_rule(tmp_path):
    """The rule has its own landing date: a pull request opened between the admission check's
    landing and the rule's gets the divergences miss as a comment; one opened after the rule's
    landing is closed for it."""
    fixture = load(FIXTURE_DIR / "maintainer-docs-with-declaration.json")
    pr = fixture["pull_request"]
    pr["user"] = {"login": "someone", "type": "User"}
    pr["state"] = "open"
    head = pr["head"]["sha"]
    workflow = [{"commit": {"committer": {"date": "1999-01-03T00:00:00Z"}}}]
    rule = [{"commit": {"committer": {"date": "1999-01-05T00:00:00Z"}}}]
    pr["created_at"] = "1999-01-04T00:00:00Z"
    gh = FakeGitHub(files=fixture["files"], landed=workflow, rule_landed=rule)
    git = FakeGit(head, commits=fixture["commits"])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("POST", f"/repos/{REPO}/issues/0/comments"),
        ("POST", f"/repos/{REPO}/issues/0/labels"),
    ]
    body = writes[0][2]["body"]
    assert f"lacks the `## {adm.DIVERGENCE_SECTION}` section" in body
    assert "opened before that section was asked" in body and "closes=0" in body
    pr["created_at"] = "1999-01-06T00:00:00Z"
    gh = FakeGitHub(files=fixture["files"], landed=workflow, rule_landed=rule)
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 1
    writes = gh.writes()
    assert ("PATCH", f"/repos/{REPO}/pulls/0", {"state": "closed"}) in writes
    assert "closes=1" in writes[0][2]["body"]


def test_run_event_pull_request_older_than_the_workflow_is_commented_not_closed(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    head = fixture["pull_request"]["head"]["sha"]
    landed = [{"commit": {"committer": {"date": "1999-01-03T00:00:00Z"}}}]
    fixture["pull_request"]["created_at"] = "1999-01-01T00:00:00Z"
    git = FakeGit(head, commits=fixture["commits"])
    gh = FakeGitHub(files=fixture["files"], landed=landed)
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [("POST", f"/repos/{REPO}/issues/0/comments")]
    assert "opened before the admission check" in gh.writes()[0][2]["body"]
    fixture["pull_request"]["created_at"] = "1999-01-04T00:00:00Z"
    gh = FakeGitHub(files=fixture["files"], landed=landed)
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 1


def test_run_event_landing_date_read_failure_never_closes(tmp_path, capsys):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    head = fixture["pull_request"]["head"]["sha"]
    gh = FakeGitHub(files=fixture["files"], fail_landed=True)
    git = FakeGit(head, commits=fixture["commits"])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [("POST", f"/repos/{REPO}/issues/0/comments")]
    assert "landing date was not read" in capsys.readouterr().err


def test_run_event_on_a_closed_pull_request_sends_no_close_and_exits_zero(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    fixture["pull_request"]["state"] = "closed"
    head = fixture["pull_request"]["head"]["sha"]
    gh = FakeGitHub(files=fixture["files"])
    git = FakeGit(head, commits=fixture["commits"])
    assert adm.run_event(event_for(fixture, "edited"), gh, git, repo_dir(tmp_path)) == 0
    assert [w[:2] for w in gh.writes()] == [("POST", f"/repos/{REPO}/issues/0/comments")]


def test_run_event_closes_a_contributor_who_edits_a_workflow_and_spares_a_maintainer(tmp_path):
    fixture = load(FIXTURE_DIR / "maintainer-docs-with-declaration.json")
    ci = text_entry(
        ".github/workflows/ci.yml",
        "jobs:\n  test:\n    name: Test\n    steps:\n      - run: cargo test\n",
        "jobs:\n  test:\n    name: Test\n    steps:\n      - run: 'true'\n",
    )
    head = fixture["pull_request"]["head"]["sha"]
    fixture["pull_request"]["state"] = "open"
    gh = FakeGitHub(files=[ci])
    git = FakeGit(head, commits=fixture["commits"])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    assert all(w[1] != f"/repos/{REPO}/pulls/0" for w in gh.writes())
    fixture["pull_request"]["user"] = {"login": "someone", "type": "User"}
    gh = FakeGitHub(files=[ci])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 1
    writes = gh.writes()
    assert ("PATCH", f"/repos/{REPO}/pulls/0", {"state": "closed"}) in writes
    assert "workflow changes are the maintainers'" in writes[0][2]["body"]
    fixture["pull_request"]["labels"] = [{"name": adm.OVERRIDE_LABEL}]
    gh = FakeGitHub(files=[ci])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    assert all(w[1] != f"/repos/{REPO}/pulls/0" for w in gh.writes())
    body = gh.writes()[0][2]["body"]
    assert "workflow changes are the maintainers'" in body and adm.OVERRIDE_LABEL in body


def test_run_event_labels_from_the_head_commit_and_skips_many_declarations(tmp_path):
    fixture = load(FIXTURE_DIR / "maintainer-docs-with-declaration.json")
    head = fixture["pull_request"]["head"]["sha"]
    commits = fixture["commits"] + [
        commit(f"{n:040x}", f"c\n\nSigned-off-by: A <a@example.com>\nCo-Authored-By: Model {n}\n")
        for n in range(1, 4)
    ]
    for entry in commits:
        entry["commit"]["author"] = fixture["commits"][0]["commit"]["author"]
        entry["commit"]["message"] = entry["commit"]["message"].replace(
            "A <a@example.com>",
            f"{entry['commit']['author']['name']} <{entry['commit']['author']['email']}>",
        )
    gh = FakeGitHub(files=fixture["files"])
    git = FakeGit(head, commits=commits)
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    writes = gh.writes()
    assert all(not w[1].endswith("/labels") for w in writes)
    assert any("label is not set" in (w[2] or {}).get("body", "") for w in writes if len(w) > 2)
    gh = FakeGitHub(files=fixture["files"])
    git = FakeGit(head, commits=commits[-2:])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    labels = [w for w in gh.writes() if w[1].endswith("/labels")]
    assert labels == [("POST", f"/repos/{REPO}/issues/0/labels", {"labels": ["model: Model 3"]})]


def test_run_event_writes_nothing_when_the_head_moved_before_the_write(tmp_path, capsys):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    head = fixture["pull_request"]["head"]["sha"]
    moved = dict(fixture["pull_request"], head={"sha": "e" * 40})
    gh = FakeGitHub(files=fixture["files"], pull=moved)
    declaring = [dict(c) for c in fixture["commits"]]
    declaring[-1]["commit"] = dict(
        declaring[-1]["commit"],
        message=declaring[-1]["commit"]["message"] + "Co-Authored-By: Example Model 1\n",
    )
    git = FakeGit(head, commits=declaring)
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 0
    assert gh.writes() == []
    assert "head moved past" in capsys.readouterr().out


def test_run_event_syncs_the_labels_of_the_pull_request_as_re_read(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    head = fixture["pull_request"]["head"]["sha"]
    re_read = dict(fixture["pull_request"], labels=[{"name": "model: Earlier Model"}])
    gh = FakeGitHub(files=fixture["files"], pull=re_read)
    declaring = [dict(c) for c in fixture["commits"]]
    declaring[-1]["commit"] = dict(
        declaring[-1]["commit"],
        message=declaring[-1]["commit"]["message"] + "Co-Authored-By: Example Model 1\n",
    )
    git = FakeGit(head, commits=declaring)
    adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path))
    label_writes = [w for w in gh.writes() if "/labels" in w[1]]
    assert label_writes == [
        ("DELETE", f"/repos/{REPO}/issues/0/labels/model%3A%20Earlier%20Model"),
        ("POST", f"/repos/{REPO}/issues/0/labels", {"labels": ["model: Example Model 1"]}),
    ]


def test_run_event_keeps_the_other_block_written_after_the_judgement(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    head = fixture["pull_request"]["head"]["sha"]
    ci = adm.ci_block([adm.Finding("ci", "The CI job `Test` failed (u): fix it and push.")], head)
    other = {
        "id": 5,
        "body": adm.render_comment({"ci": (head, ci, True)}, head, None),
        "user": {"login": "github-actions[bot]", "type": "Bot"},
    }
    gh = FakeGitHub(files=fixture["files"], comments=[other])
    git = FakeGit(head, commits=fixture["commits"])
    assert adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path)) == 1
    body = gh.writes()[0][2]["body"]
    assert "admission:ci" in body and "admission:static" in body and body.count(adm.MARKER) == 1


def test_run_dispatch_judges_the_numbered_pull_request(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    gh = FakeGitHub(files=fixture["files"], pull=fixture["pull_request"])
    head = fixture["pull_request"]["head"]["sha"]
    git = FakeGit(head, commits=fixture["commits"])
    assert adm.run_dispatch(0, gh, git, repo_dir(tmp_path)) == 1
    assert ("GET", f"/repos/{REPO}/pulls/0") in gh.calls
    assert [w[:2] for w in gh.writes()] == [
        ("POST", f"/repos/{REPO}/issues/0/comments"),
        ("PATCH", f"/repos/{REPO}/pulls/0"),
    ]


def test_run_event_head_moved_judges_nothing(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    gh = FakeGitHub(files=fixture["files"])
    code = adm.run_event(
        event_for(fixture),
        gh,
        FakeGit(fixture["pull_request"]["head"]["sha"], fetched="c" * 40, commits=fixture["commits"]),
        repo_dir(tmp_path),
    )
    assert code == 0 and gh.calls == []


def test_run_event_parity_miss_names_the_cell(tmp_path):
    fixture = load(FIXTURE_DIR / "maintainer-flag-bam.json")

    class BrokenTreeGit(FakeGit):
        def export_tree(self, rev, dest):
            write_tree(
                dest, ROW, {"crates/vep-cli/src/x.rs": TEST_FILE.replace("refuses_foo", "gone")}
            )

    gh = FakeGitHub(files=fixture["files"])
    git = BrokenTreeGit(
        fixture["pull_request"]["head"]["sha"],
        texts={"CHANGELOG.md": fixture["head_changelog"]},
        commits=fixture["commits"],
    )
    code = adm.run_event(event_for(fixture), gh, git, repo_dir(tmp_path))
    body = gh.writes()[0][2]["body"]
    assert code == 0 and "refuses_foo" in body and "the author is a maintainer" in body


def run_payload(
    head: str, event: str = "pull_request", conclusion: str = "failure", pulls: list | None = None
) -> dict:
    return {
        "workflow_run": {
            "id": 1,
            "event": event,
            "conclusion": conclusion,
            "head_sha": head,
            "head_branch": "topic",
            "head_repository": {"owner": {"login": "contributor"}},
            "pull_requests": pulls or [],
        }
    }


def test_run_workflow_run_closes_a_contributor_whose_required_check_failed(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    pr = fixture["pull_request"]
    head = pr["head"]["sha"]
    gh = FakeGitHub(
        pulls_by_head=[pr],
        rules=RULES,
        jobs=[job("Test", "failure"), job("Build and check links", "failure")],
    )
    code = adm.run_workflow_run(run_payload(head), gh, repo_dir(tmp_path))
    assert code == 1
    writes = gh.writes()
    assert [w[:2] for w in writes] == [
        ("POST", f"/repos/{REPO}/issues/0/comments"),
        ("PATCH", f"/repos/{REPO}/pulls/0"),
    ]
    body = writes[0][2]["body"]
    assert "`Test`" in body and "check links" not in body and "closed until" in body


def test_run_workflow_run_stale_head_or_other_event_judges_nothing(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    pr = fixture["pull_request"]
    gh = FakeGitHub(pulls_by_head=[pr], rules=RULES, jobs=[job("Test", "failure")])
    assert adm.run_workflow_run(run_payload("d" * 40), gh, repo_dir(tmp_path)) == 0
    assert (
        adm.run_workflow_run(run_payload(pr["head"]["sha"], event="push"), gh, repo_dir(tmp_path))
        == 0
    )
    for conclusion in ("cancelled", "action_required", "startup_failure"):
        payload = run_payload(pr["head"]["sha"], conclusion=conclusion)
        assert adm.run_workflow_run(payload, gh, repo_dir(tmp_path)) == 0
    assert gh.writes() == []


def test_run_workflow_run_success_refreshes_an_existing_comment_only(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    pr = fixture["pull_request"]
    head = pr["head"]["sha"]
    gh = FakeGitHub(pulls_by_head=[pr], rules=RULES, jobs=[job("Test", "success")])
    assert (
        adm.run_workflow_run(run_payload(head, conclusion="success"), gh, repo_dir(tmp_path)) == 0
    )
    assert gh.writes() == []
    old = adm.render_comment({"ci": (head, "After CI on head `x`:\n- old", True)}, head, None)
    existing = {"id": 9, "body": old, "user": {"type": "Bot", "login": "github-actions[bot]"}}
    gh = FakeGitHub(
        pulls_by_head=[pr], rules=RULES, jobs=[job("Test", "success")], comments=[existing]
    )
    assert (
        adm.run_workflow_run(run_payload(head, conclusion="success"), gh, repo_dir(tmp_path)) == 0
    )
    assert [w[:2] for w in gh.writes()] == [("PATCH", f"/repos/{REPO}/issues/comments/9")]
    assert "Every required CI check passes" in gh.writes()[0][2]["body"]


def test_resolve_pull_request_prints_the_open_pull_request_number_or_nothing(capsys):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    pr = fixture["pull_request"]
    gh = FakeGitHub(pulls_by_head=[pr])
    assert adm.resolve_pull_request(run_payload(pr["head"]["sha"]), gh) == 0
    assert capsys.readouterr().out.strip() == str(pr["number"])
    assert adm.resolve_pull_request(run_payload("d" * 40), gh) == 0
    assert capsys.readouterr().out == ""
    assert adm.resolve_pull_request(run_payload(pr["head"]["sha"], event="push"), gh) == 0
    assert capsys.readouterr().out == ""


def test_run_workflow_run_skips_a_closed_pull_request_the_run_lists(tmp_path, capsys):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    pr = dict(fixture["pull_request"], state="closed")
    gh = FakeGitHub(pull=pr, pulls_by_head=[], rules=RULES, jobs=[job("Test", "failure")])
    code = adm.run_workflow_run(
        run_payload(pr["head"]["sha"], pulls=[{"number": 0}]), gh, repo_dir(tmp_path)
    )
    assert code == 0 and gh.writes() == []
    assert "no open pull request has head" in capsys.readouterr().out


def test_run_workflow_run_finds_the_pull_request_by_commit_when_the_run_lists_none(tmp_path):
    fixture = load(FIXTURE_DIR / "fork-performance-change.json")
    pr = fixture["pull_request"]
    gh = FakeGitHub(pulls_by_head=[pr], rules=RULES, jobs=[job("Test", "failure")], pull=pr)
    code = adm.run_workflow_run(
        run_payload(pr["head"]["sha"], pulls=[{"number": 0}]), gh, repo_dir(tmp_path)
    )
    assert code == 1 and any(
        call[1] == f"/repos/{REPO}/pulls/0" and call[0] == "GET" for call in gh.calls
    )


# --- the local mode on a scratch repository ------------------------------------------------------


def git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True
    ).stdout


def scratch_repo(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "dev@example.com")
    git(repo, "config", "user.name", "Dev")
    (repo / "scripts" / "parity").mkdir(parents=True)
    shutil.copy(CHECKER, repo / adm.PARITY_CHECKER)
    write_tree(repo, ROW, {"crates/vep-cli/src/x.rs": TEST_FILE})
    (repo / "crates" / "vep-core" / "src").mkdir(parents=True)
    (repo / "crates" / "vep-core" / "src" / "lib.rs").write_text(
        "pub fn f() -> u8 {\n    1\n}\n\n#[cfg(test)]\nmod tests {\n    #[test]\n    fn f_is_one() {\n        assert_eq!(super::f(), 1);\n    }\n}\n"
    )
    (repo / "CHANGELOG.md").write_text(CHANGELOG)
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "Base")
    git(repo, "update-ref", "refs/remotes/origin/main", "HEAD")
    git(repo, "checkout", "-q", "-b", "topic")
    return repo


def test_run_local_passes_a_complete_branch_and_names_each_miss(tmp_path, capsys):
    repo = scratch_repo(tmp_path)
    lib = repo / "crates" / "vep-core" / "src" / "lib.rs"
    lib.write_text(lib.read_text().replace("    1\n", "    2\n").replace("f(), 1", "f(), 2"))
    (repo / "CHANGELOG.md").write_text(
        CHANGELOG.replace("- Entry (#0)\n", "- Entry (#0)\n- `f` returns two.\n")
    )
    git(repo, "commit", "-q", "-s", "-a", "-m", "Return two\n\nCo-Authored-By: No LLM")
    assert adm.run_local(repo, "origin/main", None) == 0
    out = capsys.readouterr().out
    assert "admission: OK (1 commits, 2 files" in out and "model No LLM" in out

    # A second branch: a code-only change, a script without a test, the author's address not the
    # one signed off, no declaration, no changelog line.
    git(repo, "checkout", "-q", "-b", "topic-two", "origin/main")
    lib.write_text(lib.read_text().replace("    1\n", "    3\n"))
    (repo / "scripts" / "new.sh").write_text("echo hi\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "--author=Other <other@example.com>", "-m", "Return three")
    assert adm.run_local(repo, "origin/main", None) == 1
    captured = capsys.readouterr()
    err = captured.err
    assert "no test gains a line" in err
    assert "carries no `Signed-off-by` trailer with its author's name and address" in err
    assert "declares the model" in err
    assert "gains no line under `## [Unreleased]`" in err
    assert "does not contain the tip" not in err
    assert "fix them before pushing" in captured.out


def test_run_local_reports_a_branch_behind_the_base(tmp_path, capsys):
    repo = scratch_repo(tmp_path)
    (repo / "docs.md").write_text("docs\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "Docs\n\nCo-Authored-By: No LLM")
    git(repo, "checkout", "-q", "main")
    (repo / "other.md").write_text("other\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "Advance main")
    git(repo, "update-ref", "refs/remotes/origin/main", "HEAD")
    git(repo, "checkout", "-q", "topic")
    assert adm.run_local(repo, "origin/main", None) == 1
    err = capsys.readouterr().err
    assert "does not contain the tip of `main`" in err
    assert "git fetch upstream && git merge --signoff upstream/main" in err
    with pytest.raises(adm.AdmissionError, match="run `git fetch upstream` first"):
        adm.run_local(repo, "upstream/main", None)


def test_log_commits_survives_separator_bytes_in_a_message_and_orders_parents_first(tmp_path):
    repo = scratch_repo(tmp_path)
    (repo / "a.md").write_text("a\n")
    git(repo, "add", "-A")
    hostile = "Subject\x1e with a record separator\n\nBody\x00free\x1e\n\nSigned-off-by: Dev <dev@example.com>\nCo-Authored-By: No LLM\n"
    subprocess.run(
        ["git", "-C", str(repo), "commit", "-q", "-F", "-"],
        input=hostile.replace("\x00", ""),
        text=True,
        check=True,
    )
    (repo / "b.md").write_text("b\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "Second")
    commits = adm.Git(repo).log_commits("origin/main", "HEAD")
    assert commits[0]["commit"]["message"].startswith("Subject\x1e with a record separator\n")
    assert commits[1]["commit"]["message"].startswith("Second\n")
    assert len(commits) == 2
    assert commits[0]["commit"]["author"] == {"name": "Dev", "email": "dev@example.com"}
    assert adm.declared_models(commits) == ["No LLM"]
    assert adm.label_models(commits) == ["No LLM"]
    assert adm.Git(repo).log_commits("HEAD", "HEAD") == []


def test_log_commits_reports_an_object_the_store_lost_as_an_environment_error(tmp_path, monkeypatch):
    repo = scratch_repo(tmp_path)
    (repo / "a.md").write_text("a\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "One")
    sha = git(repo, "rev-parse", "HEAD").strip()
    real_run = subprocess.run

    def pruned(args, **kwargs):
        if "cat-file" in args:
            return subprocess.CompletedProcess(args, 0, stdout=f"{sha} missing\n".encode(), stderr=b"")
        return real_run(args, **kwargs)

    monkeypatch.setattr(adm.subprocess, "run", pruned)
    with pytest.raises(adm.AdmissionError, match="is not a commit object: missing"):
        adm.Git(repo).log_commits("origin/main", "HEAD")


def test_run_local_warns_on_a_workflow_edit_without_failing_on_it(tmp_path, capsys):
    repo = scratch_repo(tmp_path)
    (repo / ".github" / "workflows").mkdir(parents=True)
    (repo / ".github" / "workflows" / "ci.yml").write_text("on: push\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "Workflow\n\nCo-Authored-By: No LLM")
    assert adm.run_local(repo, "origin/main", None) == 0
    captured = capsys.readouterr()
    assert "WARNING: [admission]" in captured.err and "workflow changes are the maintainers'" in captured.err
    assert adm.OVERRIDE_LABEL in captured.err and "admission: OK" in captured.out


def test_run_local_checks_the_divergences_section_of_a_body_file(tmp_path, capsys):
    """A bare local run certifies none; a saved description is held to the rule, and the index at
    the base and the head resolves the entries it names, the artefact test included."""
    repo = scratch_repo(tmp_path)
    (repo / "docs").mkdir(exist_ok=True)
    (repo / adm.INDEX_PATH).write_text(INDEX_SAMPLE, encoding="utf-8")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "Add the index\n\nCo-Authored-By: No LLM")
    git(repo, "update-ref", "refs/remotes/origin/main", "HEAD")
    (repo / "docs" / "divergences").mkdir()
    (repo / "docs" / "divergences" / "defect-b.md").write_text("# Defect B\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "Add the page\n\nCo-Authored-By: No LLM")
    assert adm.run_local(repo, "origin/main", None) == 1
    err = capsys.readouterr().err
    assert "`docs/divergences/defect-b.md`" in err and "under the second statement" in err
    body = tmp_path / "body.md"
    body.write_text(TEMPLATE_BODY + section({"touches"}, {"touches": ["- `B`: changes it, because the record's start codon is intact."]}))
    assert adm.run_local(repo, "origin/main", body) == 0
    assert "admission: OK" in capsys.readouterr().out
    body.write_text(TEMPLATE_BODY + section({"none", "touches"}, {"touches": ["- `B`: kept, because x."]}))
    assert adm.run_local(repo, "origin/main", body) == 1
    assert "excludes the other two" in capsys.readouterr().err
    body.write_text(TEMPLATE_BODY + section({"touches"}, {"touches": ["- `Q`: kept, because x."]}))
    assert adm.run_local(repo, "origin/main", body) == 1
    assert "`Q`" in capsys.readouterr().err
    body.write_text(TEMPLATE_BODY)
    assert adm.run_local(repo, "origin/main", body) == 1
    assert "lacks the `## Intended divergences` section" in capsys.readouterr().err


def test_export_tree_writes_links_as_files_and_skips_filters(tmp_path):
    repo = scratch_repo(tmp_path)
    (repo / "data.bin").write_bytes(b"\x00\x01")
    (repo / ".gitattributes").write_text("*.bin filter=lfs diff=lfs merge=lfs -text\n")
    (repo / "link.rs").symlink_to("crates/vep-core/src/lib.rs")
    # A machine with `git lfs install` stores the pointer at `git add`: keep the raw bytes in the tree.
    git(repo, "config", "filter.lfs.clean", "")
    git(repo, "config", "filter.lfs.process", "")
    git(repo, "config", "filter.lfs.required", "false")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-s", "-m", "Attributes and a link")
    git(repo, "config", "filter.lfs.smudge", "false")
    git(repo, "config", "filter.lfs.required", "true")
    dest = tmp_path / "export"
    adm.Git(repo).export_tree("HEAD", dest)
    assert (dest / "data.bin").read_bytes() == b"\x00\x01"
    assert not (dest / "link.rs").is_symlink()
    assert (dest / "link.rs").read_text() == "crates/vep-core/src/lib.rs"
    assert (dest / "tests" / "parity" / "matrix.tsv").read_text().startswith(MATRIX_HEADER)
    assert not (repo / "data.bin.lock").exists()


def test_main_usage_and_missing_repository(tmp_path, capsys):
    with pytest.raises(SystemExit):
        adm.main([])
    event = tmp_path / "event.json"
    event.write_text("{}")
    assert adm.main(["--event", str(event), "--repo", "", "--repo-dir", str(tmp_path)]) == 2
    assert "GITHUB_REPOSITORY" in capsys.readouterr().err
