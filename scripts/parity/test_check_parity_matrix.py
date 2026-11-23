"""Tests for check_parity_matrix.py on a small fixture tree: one test per rule, plus the exit codes."""

from __future__ import annotations

import json
import re
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


# The reference list of the fixture tree, in byte order: two scripts, modules of the three
# repositories and the BioPerl package.
REFERENCE_FILES = [
    "Bio::Tools::CodonTable",
    "Config.pm",
    "INSTALL.pl",
    "IO/Parser/BaseVCF4.pm",
    "Parser/VCF.pm",
    "Variation/Utils/Sequence.pm",
    "vep",
]
REFERENCE_LIST = "tests/parity/reference_files.txt"


def names_no_file(cell: str, where: str = "tests/parity/matrix.tsv:4: symbol") -> str:
    """The error for a `perl_reference` that is no line of the list and no plugin's module, on
    the `symbol` row unless `where` names another."""
    return (
        f"{where}: perl_reference {cell} names no file of "
        f"{REFERENCE_LIST} and no module of a plugin of the support table"
    )


def rows_with(item, **fields) -> list[str]:
    """BASE_ROWS with the row for `item` replaced by `row(item, **fields)`."""
    return [row(item, **fields) if r.startswith(item + "\t") else r for r in BASE_ROWS]


TABLE_HEADER = (
    "## Support table\n\nOne row per Ensembl VEP plugin.\n\n"
    "| Plugin | Data file and provider | Terms | Assemblies | Status | Reason |\n"
    "| --- | --- | --- | --- | --- | --- |\n"
)
REAL = "supported, verified on real data"
FIXTURES = "supported on the user's licensed file, verified on synthetic fixtures"
NOT_PORTED = "not ported"


def table_row(name, status=REAL, link=None, cells=None) -> str:
    """A support table row; `link` wraps the name in a Markdown link, `cells` replaces the four
    descriptive cells (data, terms, assemblies, reason)."""
    plugin = f"[{name}]({link})" if link else name
    data, terms, assemblies, reason = cells or ("`x.tsv.gz`; a lab", "CC0", "GRCh38", "open data")
    return f"| {plugin} | {data} | {terms} | {assemblies} | {status} | {reason} |\n"


def write_plugins_page(root: Path, table, before=TABLE_HEADER, after="\n## Usage\n") -> None:
    (root / "docs").mkdir(exist_ok=True)
    (root / "docs" / "plugins.md").write_text("# Plugins\n\n" + before + "".join(table) + after)


def build_tree(
    tmp_path: Path,
    rows=None,
    *,
    covers=("symbol",),
    fns=("symbol_column_matches",),
    with_corpus=True,
    args_rs=ARGS_RS,
    table=(),
    reference_files=REFERENCE_FILES,
) -> Path:
    root = tmp_path / "repo"
    (root / "crates" / "vep-cli" / "src").mkdir(parents=True)
    (root / "crates" / "vep-cli" / "tests").mkdir(parents=True)
    (root / "tests" / "parity").mkdir(parents=True)
    write_plugins_page(root, table)
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
    (root / REFERENCE_LIST).write_text("".join(f"{name}\n" for name in reference_files))
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
        == "parity: 16 rows, 15 flags covered, 0 landed, 0 interim entries resolve, 0 support table rows"
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


# The pin file of a fixture tree: helpers, then one test for each way a test reads a public
# text, then one for each shape that reads none.
PIN_FILE = "crates/vep-cli/tests/documented_subsets.rs"
PIN_FILE_RS = r"""
mod common;

fn page_text(path: &str) -> String {
    std::fs::read_to_string(format!("{}/../../docs/{path}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

fn assert_page_states(path: &str, phrase: &str) {
    assert!(page_text(path).contains(phrase), "docs/{path} lacks {phrase:?}");
}

fn readme() -> String {
    std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md")).unwrap()
}

fn first_line() -> String {
    readme().lines().next().unwrap().to_string()
}

fn turns() -> usize {
    turns_back()
}

fn turns_back() -> usize {
    turns()
}

fn help_of(flag: &str) -> String {
    let command = Args::command();
    let argument = command.get_arguments().find(|a| a.get_long() == Some(flag)).unwrap();
    argument.get_long_help().or(argument.get_help()).unwrap().to_string()
}

#[test]
fn names_a_page() {
    assert_page_states("cache-setup.md", "is left out");
}

#[test]
fn includes_a_page() {
    assert!(include_str!("../../../docs/output/known-variants.md").contains("is left out"));
}

#[test]
fn calls_a_reader() {
    assert!(readme().contains("is left out"));
}

#[test]
fn calls_a_reader_through_another() {
    assert!(first_line().contains("vep-rs"));
}

#[test]
fn asks_an_argument_for_its_help() {
    assert!(help_of("symbol").contains("is left out"));
}

#[test]
fn asks_an_argument_for_its_short_help() {
    let command = Args::command();
    let argument = command.get_arguments().next().unwrap();
    assert!(argument.get_help().unwrap().to_string().contains("is left out"));
}

#[test]
fn renders_the_help() {
    assert!(Args::command().render_help().to_string().contains("is left out"));
}

#[test]
fn renders_the_long_help() {
    assert!(Args::command().render_long_help().to_string().contains("is left out"));
}

#[test]
fn behaviour_alone() {
    assert_eq!(consequence_of("21\t5\t.\tA\tG"), "missense_variant");
}

#[test]
fn page_named_in_a_comment() {
    // docs/cache-setup.md states the exclusion.
    /* "cache-setup.md" */
    assert_eq!(1 + 1, 2);
}

#[test]
fn page_named_outside_a_literal() {
    let row = Row::default();
    assert!(row.md.is_empty());
}

#[test]
fn checksum_file_named_in_a_literal() {
    assert!(std::path::Path::new("variants.md5").exists());
}

#[test]
fn reader_of_another_file() {
    assert!(common::readme().contains("is left out"));
}

#[test]
fn reader_called_on_a_value() {
    assert!(Pages::default().readme().is_empty());
}

#[test]
fn reader_named_and_never_called() {
    assert_eq!("readme".len(), 6);
}

#[test]
fn help_named_in_a_literal() {
    assert_eq!("get_long_help".len(), 13);
}

#[test]
fn calls_two_functions_that_call_each_other() {
    assert_eq!(turns(), 0);
}

#[test]
fn page_named_inside_a_message() {
    assert_eq!(refusal(), "convert the cache instead (docs/cache-conversion.md), or build it");
}

#[test]
fn page_named_at_the_end_of_a_message() {
    assert_eq!(refusal(), "the file is not read, see docs/plugins.md");
}

#[test]
fn page_quoted_inside_a_raw_literal() {
    assert_eq!(refusal(), r#"the page "docs/plugins.md" states it"#);
}

fn names_a_page_without_a_test_attribute() {
    assert_page_states("cache-setup.md", "is left out");
}
"""
READS_A_PUBLIC_TEXT = [
    "names_a_page",
    "includes_a_page",
    "calls_a_reader",
    "calls_a_reader_through_another",
    "asks_an_argument_for_its_help",
    "asks_an_argument_for_its_short_help",
    "renders_the_help",
    "renders_the_long_help",
]
READS_NO_PUBLIC_TEXT = [
    "behaviour_alone",
    "page_named_in_a_comment",
    "page_named_outside_a_literal",
    "checksum_file_named_in_a_literal",
    "reader_of_another_file",
    "reader_called_on_a_value",
    "reader_named_and_never_called",
    "help_named_in_a_literal",
    "calls_two_functions_that_call_each_other",
    "page_named_inside_a_message",
    "page_named_at_the_end_of_a_message",
    "page_quoted_inside_a_raw_literal",
]
SUBSET_WITHOUT_ITS_TEXT = (
    "tests/parity/matrix.tsv:4: symbol: landed documented_subset with no test1 that reads a "
    "public text: a string literal that is a Markdown file's path or a help text asked of clap"
)


def stated(fn: str) -> str:
    """The `test1` entry of a test of the fixture pin file."""
    return f"{PIN_FILE}::{fn}"


def write_pin_file(root: Path) -> None:
    (root / PIN_FILE).write_text(PIN_FILE_RS)


def subset_rows(test1: str, state: str = "landed:documented_subset") -> list[str]:
    return rows_with("symbol", target="documented_subset", state=state, test1=test1, test2=TEST2)


def test_landed_documented_subset_row_with_both_tests_passes(tmp_path):
    root = build_tree(tmp_path, subset_rows(f"{TEST1};{stated('names_a_page')}"))
    write_pin_file(root)
    errors, summary = cpm.check(root)
    assert errors == []
    assert summary.endswith("1 landed, 0 interim entries resolve, 0 support table rows")


def test_landed_documented_subset_row_needs_a_test_that_reads_a_public_text(tmp_path):
    # A behaviour test resolves and states no exclusion.
    root = build_tree(tmp_path, subset_rows(TEST1))
    assert errors_of(root) == [SUBSET_WITHOUT_ITS_TEXT]


@pytest.mark.parametrize("fn", READS_A_PUBLIC_TEXT)
def test_a_test_that_reads_a_public_text_lands_a_documented_subset_row(tmp_path, fn):
    root = build_tree(tmp_path, subset_rows(f"{TEST1};{stated(fn)}"))
    write_pin_file(root)
    assert errors_of(root) == []
    assert cpm.reads_public_text(cpm.rust_source(PIN_FILE_RS), fn)


@pytest.mark.parametrize("fn", READS_NO_PUBLIC_TEXT)
def test_a_test_that_reads_no_public_text_lands_no_documented_subset_row(tmp_path, fn):
    # Each stands in the pin file: the file a test stands in shows nothing of what it reads.
    root = build_tree(tmp_path, subset_rows(stated(fn)))
    write_pin_file(root)
    assert errors_of(root) == [SUBSET_WITHOUT_ITS_TEXT]
    assert not cpm.reads_public_text(cpm.rust_source(PIN_FILE_RS), fn)


def test_a_public_text_is_read_by_a_unit_test_beside_the_code(tmp_path):
    entry = "crates/vep-cli/src/fields.rs::exclusion_is_stated"
    root = build_tree(tmp_path, subset_rows(entry))
    (root / "crates" / "vep-cli" / "src" / "fields.rs").write_text(
        'fn page() -> String {\n    std::fs::read_to_string("docs/output.md").unwrap()\n}\n\n'
        "#[cfg(test)]\nmod tests {\n    use super::*;\n\n    #[test]\n"
        '    fn exclusion_is_stated() {\n        assert!(page().contains("is left out"));\n    }\n}\n'
    )
    assert errors_of(root) == []


def test_only_a_landed_documented_subset_row_is_asked_for_a_public_text(tmp_path):
    assert errors_of(build_tree(tmp_path / "open", subset_rows(TEST1, state="open"))) == []
    implemented = rows_with("symbol", state="landed:implemented", test1=TEST1, test2=TEST2)
    assert errors_of(build_tree(tmp_path / "implemented", implemented)) == []


def test_a_test1_entry_that_does_not_resolve_reads_no_public_text(tmp_path):
    missing = "crates/vep-cli/tests/absent.rs::exclusion_is_stated"
    helper = stated("names_a_page_without_a_test_attribute")
    root = build_tree(tmp_path, subset_rows(f"{missing};{helper}"))
    write_pin_file(root)
    assert errors_of(root) == [
        SUBSET_WITHOUT_ITS_TEXT,
        f"tests/parity/matrix.tsv:4: symbol: test1 {missing} crates/vep-cli/tests/absent.rs does not exist",
        (
            f"tests/parity/matrix.tsv:4: symbol: test1 {helper} {PIN_FILE} declares no test fn "
            "names_a_page_without_a_test_attribute"
        ),
    ]
    # A reader that is no test reads nothing for a row either.
    assert not cpm.reads_public_text(cpm.rust_source(PIN_FILE_RS), "readme")


def test_function_bodies_are_read_from_the_brace_after_the_signature():
    source = (
        "trait Reader {\n    fn declared_alone(&self) -> [u8; 2];\n}\n"
        "fn with_brackets(x: [u8; 2], f: impl Fn(u8) -> u8) -> Vec<(u8, u8)>\nwhere\n    u8: Copy,\n"
        '{\n    let open = "{";\n    if x[0] > 0 { inner() }\n}\n'
        "fn inner() {}\n"
    )
    code = cpm.strip_comments(source)
    bodies = cpm.function_bodies(code)
    assert sorted(bodies) == ["inner", "with_brackets"]
    ((_, start, end),) = bodies["with_brackets"]
    assert code[start:end] == '{\n    let open = " ";\n    if x[0] > 0 { inner() }\n}'
    ((_, start, end),) = bodies["inner"]
    assert code[start:end] == "{}"


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
        "tests/parity/matrix.tsv:4: symbol: test1 docs/a.rs::y is not <path>.rs::<fn> under crates/, scripts/<dir>/test_<name>.py::<fn>, scripts/[<dir>/]test_<name>.sh::<check> or tools/<tool>/<path>.rs::<fn>",
        "tests/parity/matrix.tsv:4: symbol: test1 crates/../tests/parity/x.rs::outside has a `..` segment",
    ]


# A pytest module with every declaration shape the scripts resolver must tell apart.
TEST_PY = '''"""def test_in_docstring(): pass"""
import pytest

# def test_commented_out(): pass
NOTE = "def test_in_a_string(): pass"


def helper(): pass


@pytest.mark.parametrize("n", [1])
def test_decorated(n): pass


async def test_async(): pass


class TestGroup:
    def test_in_class(self): pass

    class Inner:
        def test_nested_twice(self): pass


def outer():
    def test_nested_in_function(): pass
'''
SCRIPT_TEST = "scripts/parity/test_check_parity_matrix.py"


def with_script_test(root: Path, source: str = TEST_PY) -> Path:
    (root / "scripts" / "parity").mkdir(parents=True)
    (root / "scripts" / "parity" / "test_check_parity_matrix.py").write_text(source)
    return root


@pytest.mark.parametrize(
    "fn,declared",
    [
        ("test_decorated", True),
        ("test_async", True),
        ("test_in_class", True),
        ("test_nested_twice", False),
        ("test_nested_in_function", False),
        ("test_commented_out", False),
        ("test_in_a_string", False),
        ("test_in_docstring", False),
    ],
)
def test_scripts_test_entry_resolves_to_a_pytest_function(tmp_path, fn, declared):
    """A `test1` entry may name a pytest function of a `scripts/**/test_*.py` module: one declared
    at module level or directly inside a class, never in a comment, a string or a nested function."""
    assert cpm.declares_pytest_fn(TEST_PY, fn) == (declared, None)
    rows = rows_with("symbol", test1=f"{SCRIPT_TEST}::{fn}", perl="-")
    root = with_script_test(build_tree(tmp_path, rows))
    errors = errors_of(root)
    if declared:
        assert errors == []
    else:
        assert errors == [
            f"tests/parity/matrix.tsv:4: symbol: test1 {SCRIPT_TEST}::{fn} {SCRIPT_TEST} declares no def {fn}"
        ]


def test_scripts_test_entry_shapes(tmp_path):
    """The scripts shape is `scripts/<dir>/test_<name>.py::test_<fn>` with no `..` segment, in an
    existing module that parses; `interim` keeps the crate shape."""
    rows = rows_with(
        "symbol",
        test1=(
            "scripts/parity/check_parity_matrix.py::test_x;scripts/test_top.py::test_x;"
            "scripts/parity/test_check_parity_matrix.py::helper;scripts/parity/test_gone.py::test_x;"
            "scripts/../scripts/parity/test_check_parity_matrix.py::test_x"
        ),
        interim=f"{SCRIPT_TEST}::test_decorated",
        perl="-",
    )
    root = with_script_test(build_tree(tmp_path, rows))
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: test1 scripts/parity/check_parity_matrix.py::test_x is not <path>.rs::<fn> under crates/, scripts/<dir>/test_<name>.py::<fn>, scripts/[<dir>/]test_<name>.sh::<check> or tools/<tool>/<path>.rs::<fn>",
        "tests/parity/matrix.tsv:4: symbol: test1 scripts/test_top.py::test_x is not <path>.rs::<fn> under crates/, scripts/<dir>/test_<name>.py::<fn>, scripts/[<dir>/]test_<name>.sh::<check> or tools/<tool>/<path>.rs::<fn>",
        f"tests/parity/matrix.tsv:4: symbol: test1 {SCRIPT_TEST}::helper is not <path>.rs::<fn> under crates/, scripts/<dir>/test_<name>.py::<fn>, scripts/[<dir>/]test_<name>.sh::<check> or tools/<tool>/<path>.rs::<fn>",
        "tests/parity/matrix.tsv:4: symbol: test1 scripts/parity/test_gone.py::test_x scripts/parity/test_gone.py does not exist",
        "tests/parity/matrix.tsv:4: symbol: test1 scripts/../scripts/parity/test_check_parity_matrix.py::test_x has a `..` segment",
        f"tests/parity/matrix.tsv:4: symbol: interim {SCRIPT_TEST}::test_decorated is not <path>.rs::<fn> under crates/",
    ]
    (root / "scripts" / "parity" / "test_check_parity_matrix.py").write_text("def test_x(:\n")
    rows = rows_with("symbol", test1=f"{SCRIPT_TEST}::test_x", perl="-")
    (root / "tests" / "parity" / "matrix.tsv").write_text(HEADER + "".join(rows))
    assert errors_of(root) == [
        f"tests/parity/matrix.tsv:4: symbol: test1 {SCRIPT_TEST}::test_x {SCRIPT_TEST} does not parse (invalid syntax, line 1)"
    ]


def test_process_row_perl_reference_is_dash(tmp_path):
    """A row proven by a scripts test is a process row and names no Ensembl VEP module."""
    rows = rows_with("symbol", test1=f"{SCRIPT_TEST}::test_decorated")
    root = with_script_test(build_tree(tmp_path, rows))
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: a process row names no Ensembl VEP module; perl_reference must be -"
    ]
    rows = rows_with("symbol", test1=f"{SCRIPT_TEST}::test_decorated;{TEST1}", perl="-")
    root = with_script_test(build_tree(tmp_path / "dash", rows))
    assert errors_of(root) == []
    # The crate shape alone keeps the row an Ensembl VEP row, whose perl_reference names a module.
    assert errors_of(build_tree(tmp_path / "crate", rows_with("symbol", test1=TEST1))) == []


TOOL_TEST = "tools/vep-any-compare/tests/reports.rs"


def with_tool_test(root: Path) -> Path:
    (root / "tools" / "vep-any-compare" / "tests").mkdir(parents=True)
    (root / TOOL_TEST).write_text("#[test]\nfn the_verdict_is_red() {}\n\nfn helper() {}\n")
    return root


def test_tools_test_entry_resolves_to_a_test_fn_and_makes_a_process_row(tmp_path):
    """A `test1` entry may name a test of a tool's own Cargo root: it resolves as a crate entry
    does, and its row is a process row, whose perl_reference is `-`."""
    rows = rows_with("symbol", test1=f"{TOOL_TEST}::the_verdict_is_red", perl="-")
    assert errors_of(with_tool_test(build_tree(tmp_path, rows))) == []
    rows = rows_with("symbol", test1=f"{TOOL_TEST}::the_verdict_is_red")
    assert errors_of(with_tool_test(build_tree(tmp_path / "module", rows))) == [
        "tests/parity/matrix.tsv:4: symbol: a process row names no Ensembl VEP module; perl_reference must be -"
    ]


def test_tools_test_entry_shapes(tmp_path):
    """The tools shape is `tools/<tool>/<path>.rs::<fn>` with no `..` segment, naming a function
    under a test attribute in an existing file; `interim` keeps the crate shape."""
    rows = rows_with(
        "symbol",
        test1=(
            f"{TOOL_TEST}::helper;tools/vep-any-compare/tests/gone.rs::x;tools/reports.rs::x;"
            "tools/vep-any-compare/../x.rs::outside"
        ),
        interim=f"{TOOL_TEST}::the_verdict_is_red",
        perl="-",
    )
    assert errors_of(with_tool_test(build_tree(tmp_path, rows))) == [
        f"tests/parity/matrix.tsv:4: symbol: test1 {TOOL_TEST}::helper {TOOL_TEST} declares no test fn helper",
        "tests/parity/matrix.tsv:4: symbol: test1 tools/vep-any-compare/tests/gone.rs::x tools/vep-any-compare/tests/gone.rs does not exist",
        "tests/parity/matrix.tsv:4: symbol: test1 tools/reports.rs::x is not <path>.rs::<fn> under crates/, scripts/<dir>/test_<name>.py::<fn>, scripts/[<dir>/]test_<name>.sh::<check> or tools/<tool>/<path>.rs::<fn>",
        "tests/parity/matrix.tsv:4: symbol: test1 tools/vep-any-compare/../x.rs::outside has a `..` segment",
        f"tests/parity/matrix.tsv:4: symbol: interim {TOOL_TEST}::the_verdict_is_red is not <path>.rs::<fn> under crates/",
    ]


SHELL_SUITE = "scripts/data/test_layout.sh"

# A suite in the shape the tree's shell suites have: `ok` prints a passed check's line, and the
# helpers `same`, `refused` and `refused_twice` hand their first argument on to it.
SUITE_SH = r"""#!/usr/bin/env bash
PASS=0
FAIL=0
ok() {
    PASS=$((PASS + 1))
    echo "  ok   $1"
}
bad() {
    FAIL=$((FAIL + 1))
    echo "  FAIL $1"
}
same() { # same <name> <expected> <found>
    if [ "$2" = "$3" ]; then ok "$1: equal"; else bad "$1: expected [$2], found [$3]"; fi
}
refused() {
    local name="$1" want="$2"
    shift 2
    if ! "$@" 2>/dev/null; then
        ok "$name"
    else
        bad "$name: ran"
    fi
}
refused_twice() { refused "$1" "$2" false; }
note() { echo "$1"; }
quiet() {
    # ok "$1"
    echo "$1"
}
second() { ok "$2"; }

bash -n "$0" && ok "the suite parses" || bad "the suite has a syntax error"
if true; then ok 'a single-quoted name'; else bad 'a single-quoted name'; fi
same "the layout is the committed one" a a
refused "a missing file is refused" x false
refused_twice "a helper of a helper" x
# ok "a commented-out check"
note "a line a helper prints"
quiet "a helper whose call of ok is a comment"
echo "an echoed line"
second x "a second argument"
ok "a longer name than the entry"
ok "counted $PASS checks"
"""


def with_shell_suite(root: Path, source: str = SUITE_SH) -> Path:
    (root / "scripts" / "data").mkdir(parents=True)
    (root / SHELL_SUITE).write_text(source)
    (root / "scripts" / "test_top.sh").write_text('ok() { echo "  ok   $1"; }\nok "at the top"\n')
    return root


@pytest.mark.parametrize(
    ("name", "declared"),
    [
        ("the suite parses", True),
        ("a single-quoted name", True),
        ("the layout is the committed one", True),
        ("a missing file is refused", True),
        ("a helper of a helper", True),
        ("the suite has a syntax error", False),
        ("a commented-out check", False),
        ("a line a helper prints", False),
        ("a helper whose call of ok is a comment", False),
        ("an echoed line", False),
        ("a second argument", False),
        ("a longer name", False),
        ("name than the entry", False),
    ],
)
def test_shell_suite_entry_resolves_to_a_named_check(tmp_path, name, declared):
    """A `test1` entry may name a check of a shell suite under `scripts/`: the text the suite
    gives, whole and quoted, as the first argument of `ok` or of a function of its own that hands
    its first argument on to `ok`, directly or through another such function. The text of a `bad`
    call alone, a comment line of the suite or of a function's body, a call of any other function,
    a later argument and a part of a longer name declare nothing."""
    assert cpm.shell_check_functions(SUITE_SH) == {"ok", "same", "refused", "refused_twice"}
    assert cpm.declares_shell_check(SUITE_SH, name) is declared
    rows = rows_with("symbol", test1=f"{SHELL_SUITE}::{name}", perl="-")
    errors = errors_of(with_shell_suite(build_tree(tmp_path, rows)))
    if declared:
        assert errors == []
    else:
        reason = f"{SHELL_SUITE} names no check {name!r} in a call of a check function"
        assert errors == [
            f"tests/parity/matrix.tsv:4: symbol: test1 {SHELL_SUITE}::{name} {reason}"
        ]


# Shapes of shell text that are no call of a check function, beside the ones that are: a quoted
# text, an argument of another command, a comment after code, a here-document, a script of the
# same name; and the ends of a body: a one-line definition with a comment after it, a closing
# brace that is redirected, a function inside a function, a definition line whose comment ends
# in a brace. A quoted text of several lines is no code on any of its lines, a closing brace
# inside one ends no body, a brace inside a quoted text on a definition's line closes nothing,
# a call after a body's closing brace is outside the body, a here-document's word may stand
# after a backslash or hold a hyphen, a call may follow assignments, and the text of a `$( )`
# inside double quotes is code.
SHAPES_SH = (
    r"""#!/usr/bin/env bash
ok() { echo "  ok   $1"; }
nope() { echo "not ok $1"; }
shown() { echo ok "$1"; }
later() { true # ok "$1"
}
loud() { ok "$1"; } # one line, with a comment after it
plain() { echo "$1"; }
redirected() {
    ok "$1"
} >&2
after_redirect() { echo "$1"; }
outer() {
    inner() { ok "$1"; }
    echo "$1"
}
commented_brace() { # prints {x}
    ok "$1"
}
function keyword_form {
    ok "$1"
}
local_script() {
    ./ok "$1"
}
quoted_brace() {
    awk '
BEGIN { n = 1
}
' /dev/null
    ok "$1"
}
brace_in_text() { echo "}"
    ok "$1"
}
tailed() {
    echo "$1"
}; ok "$1"
cat >"$TMP/stub.sh" <<'STUB'
stubbed() { ok "$1"; }
ok "a call inside a here-document"
STUB
cat <<\PLAIN
ok "a call inside a here-document whose word follows a backslash"
PLAIN
cat <<'END-OF-TEXT'
ok "a call inside a here-document whose word holds a hyphen"
END-OF-TEXT
awk '
    { ok "a call inside a quoted program of several lines" }
' /dev/null
text="two lines
ok 'a call inside a quoted text of several lines'
"
stub='
quoted_stub() { ok "$1"; }
'; ok "a call on the line that closes a quoted text"
"""
    + 'cat <<-TABBED\n\tok "a call inside an indented here-document"\n\tTABBED\n'
    + r"""nope "the text of a failure printer"
shown "an argument of echo"
later "a call after a comment sign in a body"
loud "a one-line helper with a comment after it"
plain "the helper after the one-line helper"
redirected "a helper whose closing brace is redirected"
after_redirect "the helper after the redirected one"
outer "an outer function that hands nothing on"
inner "an inner function"
commented_brace "a helper whose definition line ends in a commented brace"
keyword_form "a helper defined with the function keyword"
local_script "a script of the same name"
stubbed "a function defined inside a here-document"
quoted_brace "a helper whose body holds a closing brace in a quoted program"
brace_in_text "a helper whose definition line holds a closing brace in a quoted text"
tailed "a helper whose closing brace a call follows"
quoted_stub "a function defined inside a quoted text"
LC_ALL=C loud "a call of a helper after an assignment"
LC_ALL=C TZ="U T C" ok "a call after two assignments"
echo LC_ALL=C ok "an assignment that is an argument of echo"
echo "$(ok "a call inside a substitution in double quotes")"
mark="$(printf '%s' '"')"; ok "a call after a double quote inside a substitution"
echo ok "an echoed call"
true # ok "a call behind a comment sign"
ok "a call after the here-documents"
"""
)


@pytest.mark.parametrize(
    ("name", "declared"),
    [
        ("a one-line helper with a comment after it", True),
        ("a helper whose closing brace is redirected", True),
        ("an inner function", True),
        ("a helper whose definition line ends in a commented brace", True),
        ("a helper defined with the function keyword", True),
        ("a call after the here-documents", True),
        ("a helper whose body holds a closing brace in a quoted program", True),
        ("a helper whose definition line holds a closing brace in a quoted text", True),
        ("a call on the line that closes a quoted text", True),
        ("a call of a helper after an assignment", True),
        ("a call after two assignments", True),
        ("a call inside a substitution in double quotes", True),
        ("a call after a double quote inside a substitution", True),
        ("the text of a failure printer", False),
        ("an argument of echo", False),
        ("a call after a comment sign in a body", False),
        ("the helper after the one-line helper", False),
        ("the helper after the redirected one", False),
        ("an outer function that hands nothing on", False),
        ("a script of the same name", False),
        ("a function defined inside a here-document", False),
        ("a call inside a here-document", False),
        ("a call inside an indented here-document", False),
        ("an echoed call", False),
        ("a call behind a comment sign", False),
        ("a helper whose closing brace a call follows", False),
        ("a function defined inside a quoted text", False),
        ("a call inside a here-document whose word follows a backslash", False),
        ("a call inside a here-document whose word holds a hyphen", False),
        ("a call inside a quoted program of several lines", False),
        ("a call inside a quoted text of several lines", False),
        ("an assignment that is an argument of echo", False),
    ],
)
def test_shell_suite_is_read_as_code(name, declared):
    """A call of a check function is a command outside quotes, comments and here-documents, and
    a function's body ends where shell ends it: a failure printer whose text holds the word, an
    argument of another command, a script of the same name and a here-document's lines call
    nothing, and a helper that follows a one-line definition, a redirected closing brace or an
    inner function takes nothing from the body before it. A quote stays open over the lines it
    spans, so a call, a definition and a closing brace inside it are text, and a call on the
    line that closes it is one; the text of a `$( )` inside double quotes is code with quotes
    of its own; a call after a body's closing brace is no part of the body; and assignments
    before a command's name leave it a command."""
    assert cpm.shell_check_functions(SHAPES_SH) == {
        "ok",
        "loud",
        "redirected",
        "inner",
        "commented_brace",
        "keyword_form",
        "quoted_brace",
        "brace_in_text",
    }
    assert cpm.declares_shell_check(SHAPES_SH, name) is declared


def test_shell_suite_that_ends_inside_a_quote_is_read_line_by_line():
    """A quote stays open over the lines after it only in a suite whose quotes pair to its end:
    where the reader is left inside a quote at the file's end (an escaped quote inside
    `$'...'` reads as an opening one), every quote ends with its line, so a check after that
    line resolves."""
    source = (
        'ok() { echo "  ok   $1"; }\n'
        "apostrophe=$'\\''\n"
        'ok "a call after a quote the reader cannot pair"\n'
    )
    assert cpm.shell_code_lines(source) == cpm.shell_code_lines(source, spanning=False)
    assert cpm.shell_code_lines(source)[2] == 'ok "a call after a quote the reader cannot pair"'
    assert cpm.declares_shell_check(source, "a call after a quote the reader cannot pair")


def test_shell_suite_entry_shapes(tmp_path):
    """The shell shape is `scripts/[<dir>/]test_<name>.sh::<check>` with no `..` segment, in an
    existing suite; a name the suite computes is no entry, a suite that defines no `ok` names no
    check, and `test2` and `interim` keep their own shapes."""
    shapes = (
        "is not <path>.rs::<fn> under crates/, scripts/<dir>/test_<name>.py::<fn>, "
        "scripts/[<dir>/]test_<name>.sh::<check> or tools/<tool>/<path>.rs::<fn>"
    )
    rows = rows_with(
        "symbol",
        test1=(
            f"{SHELL_SUITE}::counted $PASS checks;scripts/data/layout.sh::the suite parses;"
            "docs/test_layout.sh::the suite parses;scripts/data/test_gone.sh::the suite parses;"
            "scripts/../scripts/data/test_layout.sh::the suite parses;"
            "scripts/test_top.sh::at the top;scripts/data/test_plain.sh::the suite parses"
        ),
        test2=f"{SHELL_SUITE}::the suite parses",
        interim=f"{SHELL_SUITE}::the suite parses",
        perl="-",
    )
    root = with_shell_suite(build_tree(tmp_path, rows))
    (root / "scripts" / "data" / "test_plain.sh").write_text('echo "the suite parses"\n')
    where = "tests/parity/matrix.tsv:4: symbol:"
    no_ok = (
        "scripts/data/test_plain.sh names no check 'the suite parses' in a call of a check function"
    )
    assert errors_of(root) == [
        f"{where} test1 {SHELL_SUITE}::counted $PASS checks {shapes}",
        f"{where} test1 scripts/data/layout.sh::the suite parses {shapes}",
        f"{where} test1 docs/test_layout.sh::the suite parses {shapes}",
        f"{where} test1 scripts/data/test_gone.sh::the suite parses scripts/data/test_gone.sh does not exist",
        f"{where} test1 scripts/../scripts/data/test_layout.sh::the suite parses has a `..` segment",
        f"{where} test1 scripts/data/test_plain.sh::the suite parses {no_ok}",
        f"{where} test2 {SHELL_SUITE}::the suite parses is not tests/golden/<set>/<corpus>",
        f"{where} interim {SHELL_SUITE}::the suite parses is not <path>.rs::<fn> under crates/",
    ]


def test_shell_suite_entry_makes_a_process_row(tmp_path):
    """A row proven by a check of a shell suite is a process row, as one proven by a pytest
    function is: its perl_reference is `-`, beside a crate entry too."""
    entry = f"{SHELL_SUITE}::the suite parses"
    rows = rows_with("symbol", test1=entry)
    assert errors_of(with_shell_suite(build_tree(tmp_path, rows))) == [
        "tests/parity/matrix.tsv:4: symbol: a process row names no Ensembl VEP module; perl_reference must be -"
    ]
    rows = rows_with("symbol", test1=f"{entry};{TEST1}", perl="-")
    assert errors_of(with_shell_suite(build_tree(tmp_path / "dash", rows))) == []


def test_process_row_test2_is_a_covering_corpus(tmp_path):
    """A landed process row needs its test2 like any implemented row, and the corpus resolves
    through its `covers` as for any other row."""
    rows = rows_with("symbol", state="landed:implemented", test1=f"{SCRIPT_TEST}::test_decorated", perl="-")
    root = with_script_test(build_tree(tmp_path, rows))
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: landed implemented with no test2"
    ]
    rows = rows_with(
        "symbol", state="landed:implemented", test1=f"{SCRIPT_TEST}::test_decorated", test2=TEST2, perl="-"
    )
    root = with_script_test(build_tree(tmp_path / "ok", rows))
    errors, summary = cpm.check(root)
    assert errors == []
    assert summary.endswith("1 landed, 0 interim entries resolve, 0 support table rows")
    root = with_script_test(build_tree(tmp_path / "other", rows, covers=("max_af",)))
    assert errors_of(root) == [
        "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh38 tests/golden/116/GRCh38/manifest.json covers does not name symbol"
    ]


def test_covers_entry_must_be_a_matrix_row(tmp_path):
    """Every name in a `covers` list reached through a test2 cell is a row: the manifest's own
    list and each flag set's, reported once per corpus however many rows name it."""
    rows = rows_with("symbol", test2=TEST2)
    rows = [row("max_af", test2=TEST2) if r.startswith("max_af\t") else r for r in rows]
    root = build_tree(tmp_path, rows, covers=("symbol", "max_af", "retired_row"))
    assert errors_of(root) == [
        "tests/golden/116/GRCh38/manifest.json covers names retired_row, which is not a row"
    ]
    manifest = root / "tests" / "golden" / "116" / "GRCh38" / "manifest.json"
    manifest.write_text(json.dumps({
        "covers": ["symbol"],
        "flag_sets": [{"name": "af", "flags": ["--max_af"], "covers": ["max_af", "ghost"]}],
    }))
    assert errors_of(root) == [
        "tests/golden/116/GRCh38/manifest.json flag set af covers names ghost, which is not a row"
    ]
    # A value that is not a string is reported, not raised on.
    manifest.write_text(json.dumps({"covers": ["symbol", "max_af", ["max_af"], 7], "flag_sets": []}))
    assert errors_of(root) == [
        "tests/golden/116/GRCh38/manifest.json covers names ['max_af'], which is not a row",
        "tests/golden/116/GRCh38/manifest.json covers names 7, which is not a row",
    ]
    # A corpus no test2 cell reaches is not read.
    (root / "tests" / "parity" / "matrix.tsv").write_text(HEADER + "".join(BASE_ROWS))
    assert errors_of(root) == []


def test_test2_resolves_through_a_flag_sets_covers(tmp_path):
    """A row covered by a flag set of the corpus resolves, and a malformed `flag_sets` is a reason."""
    rows = rows_with("symbol", test2=TEST2)
    root = build_tree(tmp_path, rows, covers=("max_af",))
    manifest = root / "tests" / "golden" / "116" / "GRCh38" / "manifest.json"
    manifest.write_text(json.dumps({"covers": ["max_af"], "flag_sets": [{"name": "sym", "covers": ["symbol"]}]}))
    assert errors_of(root) == []
    prefix = "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh38 tests/golden/116/GRCh38/manifest.json"
    for flag_sets, reason in [
        ({"a": 1}, "has a flag_sets that is not a list"),
        ([{"covers": ["symbol"]}], "has a flag_sets entry without a name"),
        ([{"name": "sym", "covers": "symbol"}], "flag set sym has a covers that is not a list"),
    ]:
        manifest.write_text(json.dumps({"covers": ["max_af"], "flag_sets": flag_sets}))
        assert errors_of(root) == [f"{prefix} {reason}"]


# Every literal and comment shape the stripper must read: each line pairs a source with what
# survives (comments gone, literal bodies blanked to their length, delimiters and code kept).
STRIP_CASES = [
    ('let s = "fn in_a_string";', 'let s = "              ";'),
    ('let s = "a\\"b"; // tail', 'let s = "    "; '),
    ('let s = r#"*/ and // "#; x', 'let s = r#"          "#; x'),
    ('let s = r##"q"# fn ghost() {}"##;', 'let s = r##"                 "##;'),
    ('let b = b"/* bytes */";', 'let b = b"           ";'),
    ("let c = '\\''; let d = '\"'; let e = 'x';", "let c = '  '; let d = ' '; let e = ' ';"),
    ("fn f<'a>(x: &'a str) {}", "fn f<'a>(x: &'a str) {}"),
    ("a /* one /* two */ three */ b", "a  b"),
    ("a // it's \"quoted\" /* open\nb", "a \nb"),
    ('let s = "two\nlines";', 'let s = "   \n     ";'),
]


@pytest.mark.parametrize("source,stripped", STRIP_CASES)
def test_strip_comments_blanks_literals_and_drops_comments(source, stripped):
    assert cpm.strip_comments(source) == stripped


# What survives of each source of STRIP_CASES when the literals are kept: the comments alone go.
KEPT_LITERALS = [
    'let s = "fn in_a_string";',
    'let s = "a\\"b"; ',
    'let s = r#"*/ and // "#; x',
    'let s = r##"q"# fn ghost() {}"##;',
    'let b = b"/* bytes */";',
    "let c = '\\''; let d = '\"'; let e = 'x';",
    "fn f<'a>(x: &'a str) {}",
    "a  b",
    "a \nb",
    'let s = "two\nlines";',
]


@pytest.mark.parametrize(
    "source,stripped,kept", [(*case, kept) for case, kept in zip(STRIP_CASES, KEPT_LITERALS)]
)
def test_strip_comments_keeps_literals_index_for_index(source, stripped, kept):
    assert cpm.strip_comments(source, keep_literals=True) == kept
    assert len(kept) == len(stripped)


def test_a_test_declaration_inside_a_literal_declares_nothing():
    source = 'const T: &str = "#[test] fn ghost() {}";\nconst R: &str = r#"#[test]\nfn raw_ghost() {}"#;\n#[test]\nfn real() {}\n'
    assert cpm.declares_test_fn(source, "real")
    assert not cpm.declares_test_fn(source, "ghost")
    assert not cpm.declares_test_fn(source, "raw_ghost")


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
    for not_an_object in ("[]", '"x"', "null"):
        (root / "tests" / "golden" / "116" / "GRCh37" / "manifest.json").write_text(not_an_object)
        assert errors_of(root) == [
            "tests/parity/matrix.tsv:4: symbol: test2 tests/golden/116/GRCh37 tests/golden/116/GRCh37/manifest.json is not a JSON object"
        ]
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
    assert summary.endswith("0 landed, 1 interim entries resolve, 0 support table rows")


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


@pytest.mark.parametrize(
    "cell",
    [
        "-",
        "Config.pm",
        "Config.pm::read_config_from_file",
        "Parser/VCF.pm",
        "Parser/VCF.pm::create_VariationFeatures",
        "Variation/Utils/Sequence.pm",
        "IO/Parser/BaseVCF4.pm::get_IDs",
        "vep",
        "INSTALL.pl",
        "Bio::Tools::CodonTable",
    ],
)
def test_reference_in_a_form_the_readme_gives_passes(tmp_path, cell):
    """A cell is `-`, a line of the reference list (a script, a module by its path, the BioPerl
    package), or such a line with one identifier after `::`."""
    assert errors_of(build_tree(tmp_path, rows_with("symbol", perl=cell))) == []


@pytest.mark.parametrize(
    "cell",
    [
        "VCF.pm",
        "VCF.pm::create_VariationFeatures",
        "Utils/Sequence.pm",
        "Sequence.pm",
        "BaseVCF4.pm",
        "modules/Bio/EnsEMBL/VEP/Config.pm",
        "config.pm",
        "filter_vep",
        "Bio::Tools::SeqStats",
        "Nonesuch.pm",
    ],
)
def test_reference_that_names_no_listed_file_fails(tmp_path, cell):
    """A module named by its last components, by its whole path in the repository or in another
    letter case, a script or a package the list lacks and a `.pm` name that is no module of a
    plugin of the support table are refused, each at its row."""
    assert errors_of(build_tree(tmp_path, rows_with("symbol", perl=cell))) == [names_no_file(cell)]


@pytest.mark.parametrize("sub", ["", "read config", "2nd", "a::b", "status_msg()"])
def test_reference_names_one_identifier_after_the_file(tmp_path, sub):
    """What follows `::` after a listed file, a sub or a parameter, is one identifier."""
    cell = f"Config.pm::{sub}"
    expected = (
        names_no_file(cell)
        if sub == "a::b"
        else f"tests/parity/matrix.tsv:4: symbol: perl_reference {cell}: {sub!r} after `::` is not one identifier"
    )
    assert errors_of(build_tree(tmp_path, rows_with("symbol", perl=cell))) == [expected]


def test_reference_list_is_read_whole(tmp_path):
    """A tree without the list, and a list with a line that is not one name or that does not sort
    after the line before it (a repeated line among them), stop the check at the list."""
    root = build_tree(tmp_path)
    (root / REFERENCE_LIST).unlink()
    assert errors_of(root) == [f"{REFERENCE_LIST} is missing"]
    lines = ["Config.pm", "", "Parser/ VCF.pm", "Parser/VCF.pm", "Parser/VCF.pm", "INSTALL.pl"]
    root = build_tree(tmp_path / "malformed", reference_files=lines)
    assert errors_of(root) == [
        f"{REFERENCE_LIST}:2: '' is not one name",
        f"{REFERENCE_LIST}:3: 'Parser/ VCF.pm' is not one name",
        f"{REFERENCE_LIST}:5: Parser/VCF.pm does not sort after Parser/VCF.pm",
        f"{REFERENCE_LIST}:6: INSTALL.pl does not sort after Parser/VCF.pm",
    ]


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
        == "parity: 16 rows, 15 flags covered, 0 landed, 0 interim entries resolve, 0 support table rows\n"
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


# A matrix with three plugin rows: two named `<item>.pm`, LoFTEE named `LoF.pm` as the module is.
PLUGIN_ROWS = BASE_ROWS + [
    row("AlphaMissense", perl="AlphaMissense.pm", state="landed:implemented", test1=TEST1, test2=TEST2),
    row("AVADA", target="error", perl="AVADA.pm", state="landed:error", test1=TEST1),
    row("LoFTEE", perl="LoF.pm", state="landed:implemented", test1=TEST1, test2=TEST2),
]
PLUGIN_COVERS = ("symbol", "AlphaMissense", "LoFTEE")
PLUGIN_TABLE = [
    table_row("AlphaMissense", REAL, link="plugins/annotation-store.md#alphamissense"),
    table_row("AVADA", NOT_PORTED, cells=("`avada.vcf.gz`", "not stated", "GRCh37", "one GRCh37 file")),
    table_row("LoFTEE", REAL, link="plugins/loftee.md"),
]
PAGE = "docs/plugins.md"


def test_support_table_rows_match_matrix_targets(tmp_path):
    """A table whose every row names a matrix row with the target its status asserts is clean, a
    linked name read as its text, and the summary counts the table's rows."""
    errors, summary = cpm.check(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=PLUGIN_TABLE))
    assert errors == []
    assert summary.endswith("3 landed, 0 interim entries resolve, 3 support table rows")
    table, reasons = cpm.read_support_table(tmp_path / "repo")
    assert reasons == []
    assert [(entry["name"], entry["status"], entry["lineno"]) for entry in table] == [
        ("AlphaMissense", REAL, 9),
        ("AVADA", NOT_PORTED, 10),
        ("LoFTEE", REAL, 11),
    ]


def test_support_table_row_without_matrix_row_fails(tmp_path):
    table = PLUGIN_TABLE[:2] + [table_row("Ghost", FIXTURES, link="ghost.md")] + PLUGIN_TABLE[2:]
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == [
        f"{PAGE}: support table row Ghost has no matrix row"
    ]


def test_supported_status_requires_implemented(tmp_path):
    table = [PLUGIN_TABLE[0], table_row("AVADA", REAL), PLUGIN_TABLE[2]]
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == [
        f"{PAGE}: support table row AVADA status {REAL} does not match target error"
    ]


def test_synthetic_fixtures_status_requires_implemented(tmp_path):
    table = [PLUGIN_TABLE[0], table_row("AVADA", FIXTURES), PLUGIN_TABLE[2]]
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == [
        f"{PAGE}: support table row AVADA status {FIXTURES} does not match target error"
    ]
    # The synthetic-fixtures phrase asserts `implemented` like the real-data one.
    table = [table_row("AlphaMissense", FIXTURES), PLUGIN_TABLE[1], PLUGIN_TABLE[2]]
    assert errors_of(build_tree(tmp_path / "ok", PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == []


def test_not_ported_status_requires_error(tmp_path):
    table = [table_row("AlphaMissense", NOT_PORTED), PLUGIN_TABLE[1], PLUGIN_TABLE[2]]
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == [
        f"{PAGE}: support table row AlphaMissense status not ported does not match target implemented"
    ]
    # A landed `documented_subset` target matches no phrase either.
    rows = [
        row("AVADA", target="documented_subset", perl="AVADA.pm", state="landed:documented_subset", test1=stated("names_a_page"), test2=TEST2)
        if r.startswith("AVADA\t")
        else r
        for r in PLUGIN_ROWS
    ]
    root = build_tree(tmp_path / "subset", rows, covers=PLUGIN_COVERS + ("AVADA",), table=PLUGIN_TABLE)
    write_pin_file(root)
    assert errors_of(root) == [
        f"{PAGE}: support table row AVADA status not ported does not match target documented_subset"
    ]


def test_plugin_matrix_row_absent_from_table_fails(tmp_path):
    """A matrix row whose perl_reference is `<item>.pm` needs a table row; a row named after any
    other module (`Config.pm`, the flags) does not."""
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=PLUGIN_TABLE[:1] + PLUGIN_TABLE[2:])) == [
        "tests/parity/matrix.tsv:19: AVADA is a plugin row (perl_reference AVADA.pm) absent from the support table"
    ]
    assert errors_of(build_tree(tmp_path / "flags", BASE_ROWS)) == []


def test_loftee_row_resolves_through_the_table(tmp_path):
    """LoFTEE's matrix row names `LoF.pm`, so the `.pm` rule never reaches it: its table row resolves
    through the item name, and a table without it raises no absence error."""
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == []
    assert errors_of(build_tree(tmp_path / "absent", PLUGIN_ROWS, covers=PLUGIN_COVERS, table=PLUGIN_TABLE[:2])) == []
    rows = [row("LoFTEE", target="error", perl="LoF.pm", state="landed:error", test1=TEST1) if r.startswith("LoFTEE\t") else r for r in PLUGIN_ROWS]
    assert errors_of(build_tree(tmp_path / "target", rows, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == [
        f"{PAGE}: support table row LoFTEE status {REAL} does not match target error"
    ]


def test_reference_to_a_plugin_module_resolves_through_the_table(tmp_path):
    """A row that is no plugin row may name the module of a plugin of the support table, `<Name>.pm`
    or LoFTEE's `LoF.pm`, with a sub after it too. A `.pm` name that is no such module is refused:
    a plugin's name in another letter case, `LoFTEE.pm`, which is no file, on any row, and `LoF.pm`
    where the table holds no LoFTEE. A row's own `<item>.pm` with a sub after it is refused when
    the item is no plugin of the table: the table rule reads the bare cell alone."""
    for index, cell in enumerate(["AVADA.pm", "AlphaMissense.pm::get_header_info", "LoF.pm", "LoF.pm::run"]):
        rows = [row("symbol", perl=cell) if r.startswith("symbol\t") else r for r in PLUGIN_ROWS]
        assert errors_of(build_tree(tmp_path / str(index), rows, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == []
    for cell in ["Avada.pm", "LoFTEE.pm"]:
        rows = [row("symbol", perl=cell) if r.startswith("symbol\t") else r for r in PLUGIN_ROWS]
        assert errors_of(build_tree(tmp_path / cell, rows, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == [
            names_no_file(cell)
        ]
    rows = [r.replace("\tLoF.pm\t", "\tLoFTEE.pm\t") if r.startswith("LoFTEE\t") else r for r in PLUGIN_ROWS]
    assert errors_of(build_tree(tmp_path / "loftee", rows, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == [
        names_no_file("LoFTEE.pm", "tests/parity/matrix.tsv:20: LoFTEE")
    ]
    rows = [row("symbol", perl="LoF.pm") if r.startswith("symbol\t") else r for r in PLUGIN_ROWS]
    assert errors_of(build_tree(tmp_path / "absent", rows, covers=PLUGIN_COVERS, table=PLUGIN_TABLE[:2])) == [
        names_no_file("LoF.pm")
    ]
    rows = [row("symbol", perl="symbol.pm::run") if r.startswith("symbol\t") else r for r in PLUGIN_ROWS]
    assert errors_of(build_tree(tmp_path / "own", rows, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == [
        names_no_file("symbol.pm::run")
    ]


def test_malformed_table_row_fails(tmp_path):
    """A row without six cells, a repeated name and an empty name are malformed, each at its line;
    a `\\|` inside a cell is a pipe, not a cell boundary; the matrix rules do not run on a
    malformed table."""
    malformed = f"{PAGE}: support table malformed at line"
    table = [
        PLUGIN_TABLE[0],
        "| AVADA | `avada.vcf.gz` | not stated | GRCh37 | not ported |\n",
        PLUGIN_TABLE[2],
        table_row("AlphaMissense", REAL),
        table_row("", REAL),
    ]
    root = build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)
    assert errors_of(root) == [
        f"{malformed} 10: 5 cells, not 6",
        f"{malformed} 12: AlphaMissense repeats line 9",
        f"{malformed} 13: empty plugin name",
    ]
    piped = table_row("AVADA", NOT_PORTED, cells=("`cover=0\\|1`", "not stated", "GRCh37", "one file"))
    root = build_tree(tmp_path / "pipe", PLUGIN_ROWS, covers=PLUGIN_COVERS, table=[PLUGIN_TABLE[0], piped, PLUGIN_TABLE[2]])
    assert errors_of(root) == []
    assert cpm.table_cells("| a | `b\\|c` | d |") == ["a", "`b|c`", "d"]
    # The table's shape: the page, the heading, a table under it, six header cells, a separator.
    root = build_tree(tmp_path / "shape", PLUGIN_ROWS, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)
    (root / "docs" / "plugins.md").unlink()
    assert errors_of(root) == [f"{PAGE} is missing"]
    write_plugins_page(root, PLUGIN_TABLE, before="## Plugins\n\n| Plugin | --- |\n")
    assert errors_of(root) == [f"{PAGE}: no '## Support table' heading"]
    write_plugins_page(root, [], before="## Support table\n\nProse only.\n\n")
    assert errors_of(root) == [f"{PAGE}: no table under '## Support table'"]
    write_plugins_page(root, PLUGIN_TABLE, before="## Support table\n\n| Plugin | Status |\n| --- | --- |\n")
    assert errors_of(root) == [f"{malformed} 5: the header has 2 cells, not 6"]
    header = TABLE_HEADER.splitlines()[-2] + "\n"
    write_plugins_page(root, PLUGIN_TABLE, before="## Support table\n\n" + header)
    assert errors_of(root) == [f"{malformed} 6: no separator row"]
    write_plugins_page(root, PLUGIN_TABLE, before="## Support table\n\n" + header + "| --- | --- |\n")
    assert errors_of(root) == [f"{malformed} 6: the separator has 2 cells, not 6"]


def test_unknown_status_fails(tmp_path):
    table = [PLUGIN_TABLE[0], table_row("AVADA", "ported"), PLUGIN_TABLE[2]]
    phrases = f"'{REAL}', \"{FIXTURES}\", 'not ported'"
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == [
        f"{PAGE}: support table malformed at line 10: status 'ported' is not one of {phrases}"
    ]
    # The matrix rule asserts no target for it: called directly, the row is neither a match nor a mismatch.
    rows, _ = cpm.read_matrix(tmp_path / "repo")
    parsed, _ = cpm.read_support_table(tmp_path / "repo")
    assert [entry["status"] for entry in parsed] == [REAL, "ported", REAL]
    assert cpm.check_support_table(parsed, rows) == []


def test_supported_status_requires_a_landed_row(tmp_path):
    """A supported status is a claim about shipped code: a matrix row still open, whatever its target,
    cannot carry one, while `not ported` fits any open row, a port decided but not landed included."""
    alphamissense_open = [row("AlphaMissense", perl="AlphaMissense.pm") if r.startswith("AlphaMissense\t") else r for r in PLUGIN_ROWS]
    root = build_tree(tmp_path, alphamissense_open, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)
    assert errors_of(root) == [
        f"{PAGE}: support table row AlphaMissense status {REAL} names an open matrix row"
    ]
    table = [table_row("AlphaMissense", FIXTURES), PLUGIN_TABLE[1], PLUGIN_TABLE[2]]
    root = build_tree(tmp_path / "fixtures", alphamissense_open, covers=PLUGIN_COVERS, table=table)
    assert errors_of(root) == [
        f"{PAGE}: support table row AlphaMissense status {FIXTURES} names an open matrix row"
    ]
    avada_open = [row("AVADA", target="error", perl="AVADA.pm") if r.startswith("AVADA\t") else r for r in PLUGIN_ROWS]
    assert errors_of(build_tree(tmp_path / "open-error", avada_open, covers=PLUGIN_COVERS, table=PLUGIN_TABLE)) == []
    table = [table_row("AlphaMissense", NOT_PORTED), PLUGIN_TABLE[1], PLUGIN_TABLE[2]]
    assert errors_of(build_tree(tmp_path / "open-port", alphamissense_open, covers=PLUGIN_COVERS, table=table)) == []


def test_table_rows_keep_alphabetical_order(tmp_path):
    """The rows follow the plugin name, case ignored, so a lower-case `avocado` sits between `AVADA`
    and `LoFTEE`; a row out of place is named with the row before it."""
    table = [PLUGIN_TABLE[1], PLUGIN_TABLE[0], PLUGIN_TABLE[2]]
    assert errors_of(build_tree(tmp_path, PLUGIN_ROWS, covers=PLUGIN_COVERS, table=table)) == [
        f"{PAGE}: support table row AlphaMissense is out of alphabetical order after AVADA"
    ]
    rows = PLUGIN_ROWS + [row("avocado", perl="avocado.pm", target="error", state="landed:error", test1=TEST1)]
    table = PLUGIN_TABLE[:2] + [table_row("avocado", NOT_PORTED)] + PLUGIN_TABLE[2:]
    assert errors_of(build_tree(tmp_path / "casefold", rows, covers=PLUGIN_COVERS, table=table)) == []


README_PLUGINS_ROW_RE = re.compile(r"^\| \*\*Plugins\*\* +\| (\d+) of the (\d+) Ensembl VEP plugins built in", re.MULTILINE)


def test_repository_table_is_clean():
    """The repository's support table parses, every plugin row of the matrix and every table row
    agree, which `cpm.check` proves in full on the tree, and the README's plugin count is the
    table's."""
    repo = HERE.parent.parent
    table, reasons = cpm.read_support_table(repo)
    assert reasons == []
    rows, _ = cpm.read_matrix(repo)
    assert cpm.check_support_table(table, rows) == []
    assert {entry["name"] for entry in table} >= {
        row["item"] for row in rows if row["perl_reference"] == f"{row['item']}.pm"
    } | {"LoFTEE"}
    supported = sum(1 for entry in table if entry["status"] != NOT_PORTED)
    match = README_PLUGINS_ROW_RE.search((repo / "README.md").read_text(encoding="utf-8"))
    assert match is not None, "README.md has no Plugins row counting the built-in plugins"
    assert (int(match.group(1)), int(match.group(2))) == (supported, len(table))


def test_repository_matrix_is_clean():
    repo = HERE.parent.parent
    errors, summary = cpm.check(repo)
    assert errors == []
    assert summary.startswith("parity: ")
