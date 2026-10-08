"""Tests for check_parity_matrix.py on a small fixture tree: one test per rule, plus the exit codes."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "check_parity_matrix.py"
sys.path.insert(0, str(HERE))
import check_parity_matrix as cpm

RUST_TEST = "crates/vep-cli/tests/parity_matrix.rs"

# Every field shape clap accepts that the reader must follow to its flag.
ARGS_RS = """
#[derive(Parser, Debug)]
pub struct Args {
    /// Input file path
    #[arg(short = 'i', long = "input_file", default_value = "STDIN")]
    pub input_file: String,

    /// Row shape, a value list whose brackets must not end the attribute early
    #[arg(
        long = "parquet_shape",
        default_value = "nested",
        value_parser = ["nested", "flat"]
    )]
    pub parquet_shape: String,

    /// A bare `long` takes the field name
    #[arg(long)]
    pub symbol: bool,

    /// A bare `long` on a two-word field is kebab-case on the command line
    #[arg(long, default_value_t = 0)]
    pub max_af: usize,

    /// Aliases are not flags of their own
    #[arg(long = "force_overwrite", visible_alias = "force")]
    pub force_overwrite: bool,

    /// A short-only option has no long name
    #[arg(short = 'z')]
    pub zed: bool,

    #[arg(long)]
    hidden_private: bool,

    #[arg(long)]
    pub(crate) crate_vis: bool,

    #[arg(long)]
    // a line comment between the attribute and the field
    pub after_line_comment: bool,

    #[arg(long)]
    #[doc = "see [link]"]
    pub after_doc_brackets: bool,

    #[arg(long) ]
    pub space_before_bracket: bool,

    #[arg(long)]
    pub r#type: bool,

    #[arg(long = r"raw_long")]
    pub raw_long: bool,

    #[clap(long)]
    pub legacy_attr: bool,

    #[arg (long)]
    pub spaced_attr: bool,

    #[arg(id = "renamed_id", long)]
    pub original_field: bool,
}
"""
FLAGS = [
    "input_file",
    "parquet_shape",
    "symbol",
    "max-af",
    "force_overwrite",
    "hidden-private",
    "crate-vis",
    "after-line-comment",
    "after-doc-brackets",
    "space-before-bracket",
    "type",
    "raw_long",
    "legacy-attr",
    "spaced-attr",
    "renamed_id",
]

# Declarations the test-attribute rule must tell apart.
FIELDS_RS = """
// fn commented_out() {}
/* fn in_a_block_comment() {} */
const S: &str = "fn in_a_string";
const OPENER: &str = "U/*";
const QUOTE: char = '"';
pub fn production_helper() -> u8 { 1 }

#[cfg(test)]
mod tests {
    #[test]
    fn symbol_column_matches() {}

    #[tokio::test]
    async fn async_case() {}

    #[rstest]
    #[case(1)]
    fn rstest_case(#[case] n: u8) {}

    #[test]
    #[ignore]
    pub(crate) fn ignored_case() {}
}
const CLOSER: &str = r#"*/ and // inside a raw string"#;
"""

HEADER = "\t".join(cpm.COLUMNS) + "\n"


def row(
    item,
    target="implemented",
    state="open",
    test1="",
    test2="",
    interim="",
    perl="Config.pm",
):
    return "\t".join([item, perl, target, state, test1, test2, interim]) + "\n"


BASE_ROWS = [
    row("input_file"),
    row("parquet_shape", target="removed"),
    row("symbol"),
    row("max_af"),
    row("force_overwrite"),
    row("rest_service_surface", target="excluded", perl="-"),
] + [row(flag.replace("-", "_")) for flag in FLAGS[5:]]


def rows_with(item, **fields) -> list[str]:
    """BASE_ROWS with the row for `item` replaced by `row(item, **fields)`."""
    return [row(item, **fields) if r.startswith(item + "\t") else r for r in BASE_ROWS]


def build_tree(
    tmp_path: Path,
    rows=None,
    *,
    covers=("symbol",),
    fns=("symbol_column_matches",),
    with_corpus=True,
    args_rs=ARGS_RS,
) -> Path:
    root = tmp_path / "repo"
    (root / "crates" / "vep-cli" / "src").mkdir(parents=True)
    (root / "crates" / "vep-cli" / "tests").mkdir(parents=True)
    (root / "tests" / "parity").mkdir(parents=True)
    (root / "crates" / "vep-cli" / "src" / "args.rs").write_text(args_rs)
    body = "\n".join(f"#[test]\nfn {fn}() {{}}\n" for fn in fns)
    (root / "crates" / "vep-cli" / "src" / "fields.rs").write_text(body)
    (root / "crates" / "vep-cli" / "tests" / "golden.rs").write_text(
        "#[test]\nfn golden_corpora_match() {}\n"
    )
    if with_corpus:
        corpus = root / "tests" / "golden" / "116" / "GRCh38"
        corpus.mkdir(parents=True)
        (corpus / "manifest.json").write_text(
            json.dumps({"covers": list(covers)})
        )
    (root / "tests" / "parity" / "matrix.tsv").write_text(
        HEADER + "".join(BASE_ROWS if rows is None else rows)
    )
    return root


def errors_of(root: Path) -> list[str]:
    errors, _ = cpm.check(root)
    return errors


TEST1 = "crates/vep-cli/src/fields.rs::symbol_column_matches"
TEST2 = "tests/golden/116/GRCh38"


def test_clean_tree_passes_with_summary(tmp_path):
    root = build_tree(tmp_path)
    errors, summary = cpm.check(root)
    assert errors == []
    assert (
        summary
        == "parity: 16 rows, 15 flags covered, 0 landed, 0 interim entries resolve"
    )


def test_scan_args_follows_every_field_shape():
    assert cpm.scan_args(ARGS_RS) == (FLAGS, [])


@pytest.mark.parametrize(
    "snippet,error",
    [
        (
            "    #[arg(long)]\n    fn not_a_field() {}\n",
            "crates/vep-cli/src/args.rs:5: #[arg(...)] is not followed by a field this check can read",
        ),
        (
            "    #[arg(long)\n    pub unclosed: bool,\n",
            "crates/vep-cli/src/args.rs:5: attribute is not closed by `]`",
        ),
        (
            "    #[command(flatten)]\n    pub extra: Flattened,\n",
            f"crates/vep-cli/src/args.rs:5: #[command(flatten is not read by this check; {RUST_TEST} covers the flags it brings",
        ),
        (
            "    #[clap(subcommand)]\n    pub action: Action,\n",
            f"crates/vep-cli/src/args.rs:5: #[clap(subcommand is not read by this check; {RUST_TEST} covers the flags it brings",
        ),
        (
            '    #[command(rename_all = "snake_case")]\n    pub inner: Inner,\n',
            f"crates/vep-cli/src/args.rs:5: rename_all is not read by this check; {RUST_TEST} covers the flags it brings",
        ),
    ],
)
def test_scan_args_refuses_what_it_cannot_read(tmp_path, snippet, error):
    source = (
        "#[derive(Parser)]\npub struct Args {\n    #[arg(long)]\n    pub symbol: bool,\n"
        + snippet
        + "}\n"
    )
    assert cpm.scan_args(source) == (["symbol"], [error])
    root = build_tree(tmp_path, [row("symbol")], args_rs=source)
    assert errors_of(root) == [error]


def test_flag_without_a_row_fails(tmp_path):
    rows = [r for r in BASE_ROWS if not r.startswith("max_af\t")]
    errors = errors_of(build_tree(tmp_path, rows))
    assert errors == [
        "crates/vep-cli/src/args.rs: flag --max-af has no row 'max_af' in tests/parity/matrix.tsv"
    ]


def test_landed_row_needs_test1(tmp_path):
    rows = rows_with("symbol", state="landed:implemented", test2=TEST2)
    errors = errors_of(build_tree(tmp_path, rows))
    assert errors == ["tests/parity/matrix.tsv:4: symbol: landed with no test1"]


def test_landed_implemented_row_needs_test2(tmp_path):
    rows = rows_with("symbol", state="landed:implemented", test1=TEST1)
    errors = errors_of(build_tree(tmp_path, rows))
    assert errors == [
        "tests/parity/matrix.tsv:4: symbol: landed implemented with no test2"
    ]


@pytest.mark.parametrize("target", ["accepted_noop", "error"])
def test_landed_row_of_other_targets_needs_test1_only(tmp_path, target):
    rows = rows_with("symbol", target=target, state=f"landed:{target}", test1=TEST1)
    assert errors_of(build_tree(tmp_path, rows)) == []


def test_landed_removed_row_needs_its_flag_gone(tmp_path):
    rows = rows_with(
        "parquet_shape", target="removed", state="landed:removed", test1=TEST1
    )
    assert errors_of(build_tree(tmp_path, rows)) == [
        "tests/parity/matrix.tsv:3: parquet_shape: landed as removed while crates/vep-cli/src/args.rs still declares the flag"
    ]
    gone = ARGS_RS.replace('long = "parquet_shape"', 'long = "parquet_layout"')
    rows = rows + [row("parquet_layout")]
    assert errors_of(build_tree(tmp_path / "gone", rows, args_rs=gone)) == []


def test_landed_documented_subset_row_with_both_tests_passes(tmp_path):
    rows = rows_with(
        "symbol",
        target="documented_subset",
        state="landed:documented_subset",
        test1=TEST1,
        test2=TEST2,
    )
    errors, summary = cpm.check(build_tree(tmp_path, rows))
    assert errors == []
    assert summary.endswith("1 landed, 0 interim entries resolve")


def test_landed_row_has_no_interim(tmp_path):
    rows = rows_with(
        "symbol",
        target="error",
        state="landed:error",
        test1=TEST1,
        interim="crates/vep-cli/tests/golden.rs::golden_corpora_match",
    )
    assert errors_of(build_tree(tmp_path, rows)) == [
        "tests/parity/matrix.tsv:4: symbol: a landed row has no interim"
    ]


def test_test1_must_name_an_existing_function(tmp_path):
    rows = rows_with(
        "symbol",
        test1="crates/vep-cli/src/fields.rs::no_such_fn;crates/vep-cli/src/gone.rs::x;docs/a.rs::y;crates/../tests/parity/x.rs::outside",
    )
    errors = errors_of(build_tree(tmp_path, rows))
    assert errors == [
        "tests/parity/matrix.tsv:4: symbol: test1 crates/vep-cli/src/fields.rs::no_such_fn crates/vep-cli/src/fields.rs declares no test fn no_such_fn",
        "tests/parity/matrix.tsv:4: symbol: test1 crates/vep-cli/src/gone.rs::x crates/vep-cli/src/gone.rs does not exist",
        "tests/parity/matrix.tsv:4: symbol: test1 docs/a.rs::y is not <path>.rs::<fn> under crates/",
        "tests/parity/matrix.tsv:4: symbol: test1 crates/../tests/parity/x.rs::outside has a `..` segment",
    ]


@pytest.mark.parametrize(
    "fn,declared",
    [
        ("symbol_column_matches", True),
        ("async_case", True),
        ("rstest_case", True),
        ("ignored_case", True),
        ("commented_out", False),
        ("in_a_block_comment", False),
        ("in_a_string", False),
        ("production_helper", False),
    ],
)
def test_test1_needs_a_test_attribute(tmp_path, fn, declared):
    assert cpm.declares_test_fn(FIELDS_RS, fn) is declared
    root = build_tree(
        tmp_path, rows_with("symbol", test1=f"crates/vep-cli/src/fields.rs::{fn}")
    )
    (root / "crates" / "vep-cli" / "src" / "fields.rs").write_text(FIELDS_RS)
    errors = errors_of(root)
    if declared:
        assert errors == []
    else:
        assert errors == [
            f"tests/parity/matrix.tsv:4: symbol: test1 crates/vep-cli/src/fields.rs::{fn} crates/vep-cli/src/fields.rs declares no test fn {fn}"
        ]


def test_several_test1_entries_all_resolve(tmp_path):
    rows = rows_with(
        "symbol",
        test1=f"{TEST1}; crates/vep-cli/tests/golden.rs::golden_corpora_match",
    )
    assert errors_of(build_tree(tmp_path, rows)) == []


def test_test2_corpus_must_exist_with_a_manifest_covering_the_row(tmp_path):
    rows = rows_with("symbol", test2="tests/golden/116/GRCh37")
    assert errors_of(build_tree(tmp_path, rows)) == [
        "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh37 tests/golden/116/GRCh37 does not exist"
    ]
    root = build_tree(tmp_path / "b", rows)
    (root / "tests" / "golden" / "116" / "GRCh37").mkdir()
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh37 tests/golden/116/GRCh37 has no manifest.json"
    ]
    (root / "tests" / "golden" / "116" / "GRCh37" / "manifest.json").write_text("{")
    assert errors_of(root)[0].startswith(
        "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh37 tests/golden/116/GRCh37/manifest.json is not JSON"
    )
    (root / "tests" / "golden" / "116" / "GRCh37" / "manifest.json").write_text(
        json.dumps({"exemplars": []})
    )
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh37 tests/golden/116/GRCh37/manifest.json has no covers list"
    ]


def test_test2_corpus_must_list_the_row_in_covers(tmp_path):
    rows = rows_with("symbol", test2=TEST2)
    errors = errors_of(build_tree(tmp_path, rows, covers=("max_af",)))
    assert errors == [
        "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh38 tests/golden/116/GRCh38/manifest.json covers does not name symbol"
    ]


def test_test2_must_be_a_corpus_path(tmp_path):
    rows = rows_with("symbol", test2=TEST1)
    errors = errors_of(build_tree(tmp_path, rows))
    assert errors == [
        f"tests/parity/matrix.tsv:4: symbol: test2 {TEST1} is not tests/golden/<set>/<corpus>"
    ]


def test_interim_resolves_to_a_test_function(tmp_path):
    rows = rows_with(
        "symbol", interim="crates/vep-cli/src/fields.rs::refused_until_landed"
    )
    root = build_tree(tmp_path, rows)
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: interim crates/vep-cli/src/fields.rs::refused_until_landed crates/vep-cli/src/fields.rs declares no test fn refused_until_landed"
    ]
    root2 = build_tree(
        tmp_path / "ok", rows, fns=("symbol_column_matches", "refused_until_landed")
    )
    errors, summary = cpm.check(root2)
    assert errors == []
    assert summary.endswith("0 landed, 1 interim entries resolve")


def test_landed_at_another_target_fails(tmp_path):
    rows = rows_with(
        "symbol", state="landed:documented_subset", test1=TEST1, test2=TEST2
    )
    assert errors_of(build_tree(tmp_path, rows)) == [
        "tests/parity/matrix.tsv:4: symbol: landed at 'documented_subset', its target is 'implemented'"
    ]
    rows = rows_with("symbol", state="landed:implemented", test1=TEST1, test2=TEST2)
    assert errors_of(build_tree(tmp_path / "same", rows)) == []


def test_bare_landed_is_an_error(tmp_path):
    rows = rows_with("symbol", state="landed", test1=TEST1, test2=TEST2)
    assert errors_of(build_tree(tmp_path, rows)) == [
        "tests/parity/matrix.tsv:4: symbol: state 'landed' must name the target it landed at: landed:<target>"
    ]


def test_excluded_row_never_lands(tmp_path):
    rows = rows_with(
        "rest_service_surface",
        target="excluded",
        state="landed:excluded",
        test1="crates/vep-cli/tests/golden.rs::golden_corpora_match",
        perl="-",
    )
    assert errors_of(build_tree(tmp_path, rows)) == [
        "tests/parity/matrix.tsv:7: rest_service_surface: an excluded row never lands"
    ]


def test_header_and_vocabulary_are_checked(tmp_path):
    root = build_tree(tmp_path)
    (root / "tests" / "parity" / "matrix.tsv").write_text("item\ttarget\n")
    assert errors_of(root) == [
        "tests/parity/matrix.tsv: header is ['item', 'target'], expected "
        + str(cpm.COLUMNS)
    ]
    (root / "tests" / "parity" / "matrix.tsv").write_text(
        HEADER
        + row("input_file", target="ported")
        + row("input_file", perl="")
        + row("", state="done")
        + "short\trow\n"
    )
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:5: 2 columns, expected 7",
        "tests/parity/matrix.tsv:2: input_file: target 'ported' is not one of implemented, documented_subset, accepted_noop, error, removed, excluded",
        "tests/parity/matrix.tsv:3: item 'input_file' repeats line 2",
        "tests/parity/matrix.tsv:3: input_file: empty perl_reference (`-` names none)",
        "tests/parity/matrix.tsv:4: empty item",
        "tests/parity/matrix.tsv:4: : state 'done' is not open or landed:<target>",
    ]


def test_missing_files_are_reported(tmp_path):
    root = build_tree(tmp_path)
    (root / "tests" / "parity" / "matrix.tsv").unlink()
    assert errors_of(root) == ["tests/parity/matrix.tsv is missing"]
    root = build_tree(tmp_path / "b")
    (root / "crates" / "vep-cli" / "src" / "args.rs").unlink()
    assert errors_of(root) == ["crates/vep-cli/src/args.rs is missing"]


def run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_exit_codes_and_output(tmp_path):
    root = build_tree(tmp_path)
    ok = run_cli("--root", str(root))
    assert ok.returncode == 0
    assert (
        ok.stdout
        == "parity: 16 rows, 15 flags covered, 0 landed, 0 interim entries resolve\n"
    )
    assert ok.stderr == ""
    (root / "tests" / "parity" / "matrix.tsv").write_text(
        HEADER + "".join(rows_with("symbol", state="landed:implemented"))
    )
    bad = run_cli("--root", str(root))
    assert bad.returncode == 1
    assert bad.stdout == ""
    assert bad.stderr == (
        "ERROR: [parity] tests/parity/matrix.tsv:4: symbol: landed with no test1\n"
        "ERROR: [parity] tests/parity/matrix.tsv:4: symbol: landed implemented with no test2\n"
    )


def test_repository_matrix_is_clean():
    repo = HERE.parent.parent
    errors, summary = cpm.check(repo)
    assert errors == []
    assert summary.startswith("parity: ")
