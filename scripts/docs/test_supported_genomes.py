"""Tests for the supported-genomes page's writer, on the committed tables and on copies the tests change: the
committed page is what the committed tables give, a table that gains, loses or changes a row fails the check, the
page's rows and counts are the tables' own, and the species page's corpus rows are the corpora the tree holds."""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import render_supported_genomes as genomes

PAGE = ROOT / genomes.PAGE
SPECIES_PAGE = ROOT / genomes.SPECIES_PAGE
TABLE_DIR = ROOT / genomes.TABLE_DIR
NEW_ROW = ["zz_made_up_genome", "MadeUp1", "ensembl", "", "63", "1", "20", "exact", "0" * 64]
EXCLUDED_OUTCOME = "excluded:made-up-reason"
EXCLUDED_MEANING = "The cache was left out for a made-up reason."
EXCLUDED_SENTENCE = "A cache that was excluded was never compared."
DEFAULT_FLAGS_SENTENCE = "Both engines ran under default flags."
GIVEN_REFERENCE_SENTENCE = (
    "On the human RefSeq and merged caches, whose transcripts are edited from an alignment file, both were also "
    "given `--use_given_ref` ([Genome checks](genome-checks.md))."
)


def table_rows(name: str) -> list[dict[str, str]]:
    """The rows of one committed table file, read without the writer's own reader."""
    opener = gzip.open if name.endswith(".gz") else open
    with opener(TABLE_DIR / name, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def all_rows() -> list[dict[str, str]]:
    return table_rows("genomes.tsv") + table_rows("bacteria.tsv.gz")


def page_sections() -> dict[str, list[list[str]]]:
    """Level-two heading of the committed page -> the cells of each table row under it, header rows left out."""
    sections: dict[str, list[list[str]]] = {}
    heading = ""
    for line in PAGE.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            heading = line[3:]
            sections[heading] = []
        elif line.startswith("| ") and heading:
            cells = [cell.strip() for cell in line.strip("|").split("|")]
            if not set(cells[0]) <= {"-", " "}:
                sections[heading].append(cells)
    return {name: rows[1:] for name, rows in sections.items() if rows}


def plain(cell: str) -> str:
    return cell.replace("`", "")


def count(cell: str) -> int:
    return int(cell.replace(",", ""))


@pytest.fixture
def tree(tmp_path: Path) -> Path:
    """A copy of the files the writer reads, the corpora linked, for a test to change."""
    (tmp_path / "docs").mkdir()
    shutil.copy(PAGE, tmp_path / genomes.PAGE)
    shutil.copy(SPECIES_PAGE, tmp_path / genomes.SPECIES_PAGE)
    shutil.copytree(TABLE_DIR, tmp_path / genomes.TABLE_DIR)
    (tmp_path / "tests" / "golden").symlink_to(ROOT / "tests" / "golden", target_is_directory=True)
    return tmp_path


def rewrite(tree: Path, name: str, change) -> None:
    """Pass a table file's lines, header included, through `change` and write them back."""
    path = tree / genomes.TABLE_DIR / name
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        lines = handle.read().splitlines()
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(change(lines)) + "\n")


def recount(tree: Path, set_name: str, by: int) -> None:
    def change(lines: list[str]) -> list[str]:
        out = []
        for line in lines:
            cells = line.split("\t")
            if cells[0] == set_name:
                cells[3] = str(int(cells[3]) + by)
            out.append("\t".join(cells))
        return out

    rewrite(tree, "sets.tsv", change)


def add_row(tree: Path, name: str, set_name: str) -> None:
    """Append a made-up row of `set_name`, under the set's own cache version where `sets.tsv` lists the set."""
    row = list(NEW_ROW)
    row[3] = set_name
    listed = (tree / genomes.TABLE_DIR / "sets.tsv").read_text(encoding="utf-8").splitlines()
    row[4] = next((line.split("\t")[2] for line in listed if line.split("\t")[0] == set_name), row[4])
    rewrite(tree, name, lambda lines: [*lines, "\t".join(row)])
    recount(tree, set_name, 1)


def exclude_one(tree: Path, set_name: str) -> list[str]:
    """Declare a made-up excluded outcome in the copy and give it, with no digest, as the genome job writes a cache
    it excludes, to the first `exact` row `genomes.tsv` holds of `set_name`. Returns the row's cells as they stood."""
    held: list[list[str]] = []

    def change(lines: list[str]) -> list[str]:
        out = []
        for line in lines:
            cells = line.split("\t")
            if not held and cells[3] == set_name and cells[7] == "exact":
                held.append(list(cells))
                cells[7], cells[8] = EXCLUDED_OUTCOME, "-"
            out.append("\t".join(cells))
        return out

    rewrite(tree, "outcomes.tsv", lambda lines: [*lines, f"{EXCLUDED_OUTCOME}\t{EXCLUDED_MEANING}"])
    rewrite(tree, "genomes.tsv", change)
    assert len(held) == 1
    return held[0]


def compare_every_row(tree: Path) -> None:
    """Give every excluded row of the copy's two table files the outcome `exact` and a digest."""

    def change(lines: list[str]) -> list[str]:
        out = []
        for line in lines:
            cells = line.split("\t")
            if cells[7].startswith("excluded:"):
                cells[7], cells[8] = "exact", "0" * 64
            out.append("\t".join(cells))
        return out

    for name in ("genomes.tsv", "bacteria.tsv.gz"):
        rewrite(tree, name, change)


def four_counts(page: str) -> list[tuple[str, int, int]]:
    """(label, exact, compared) for each figure of a page's line of four counts."""
    line = next(line for line in page.splitlines() if line.startswith("Exact of compared: "))
    return [
        (label, count(exact), count(compared))
        for label, exact, compared in re.findall(r"([a-z][a-z ]+?) ([\d,]+) of ([\d,]+)[;.]", line.split(": ", 1)[1])
    ]


def check(tree: Path, capsys: pytest.CaptureFixture[str]) -> tuple[int, str]:
    code = genomes.main(["--repo", str(tree), "--check"])
    return code, capsys.readouterr().err


def test_the_committed_page_is_what_the_committed_tables_give(capsys):
    assert check(ROOT, capsys) == (0, "")


def test_an_unchanged_copy_passes_the_check(tree, capsys):
    assert check(tree, capsys) == (0, "")


def test_the_same_tree_gives_the_same_bytes():
    assert genomes.render(ROOT).encode("utf-8") == genomes.render(ROOT).encode("utf-8") == PAGE.read_bytes()


@pytest.mark.parametrize("set_name", ["plants", "pangenome"])
def test_a_table_row_without_a_page_row_fails_the_check(tree, capsys, set_name):
    add_row(tree, "genomes.tsv", set_name)
    code, err = check(tree, capsys)
    assert code == 1
    assert err.startswith("ERROR: [render_supported_genomes] docs/supported-genomes.md is not what the tables give: line ")
    assert err.count("\n") == 1


@pytest.mark.parametrize(
    ("name", "set_name"),
    [("bacteria.tsv.gz", "bacteria_0_collection"), ("genomes.tsv", "fungi_ascomycota1_collection")],
)
def test_a_row_added_to_a_summarised_set_fails_the_check_on_its_count(tree, capsys, name, set_name):
    add_row(tree, name, set_name)
    code, err = check(tree, capsys)
    assert code == 1
    assert "is not what the tables give" in err


def test_a_page_row_without_a_table_row_fails_the_check(tree, capsys):
    removed = []

    def change(lines: list[str]) -> list[str]:
        kept = [line for line in lines if not line.startswith("zea_mays\t")]
        removed.append(len(lines) - len(kept))
        return kept

    rewrite(tree, "genomes.tsv", change)
    assert removed == [1]
    recount(tree, "plants", -1)
    code, err = check(tree, capsys)
    assert code == 1
    assert "is not what the tables give" in err
    assert genomes.main(["--repo", str(tree)]) == 0
    assert "`zea_mays`" not in (tree / genomes.PAGE).read_text(encoding="utf-8")


def test_a_row_whose_outcome_changes_fails_the_check(tree, capsys):
    def change(lines: list[str]) -> list[str]:
        out = []
        for line in lines:
            cells = line.split("\t")
            if cells[0] == "zea_mays":
                cells[7], cells[8] = "differs:unclassified", "-"
            out.append("\t".join(cells))
        return out

    rewrite(tree, "genomes.tsv", change)
    code, err = check(tree, capsys)
    assert code == 1
    assert "is not what the tables give" in err


def test_a_page_that_is_missing_fails_the_check(tree, capsys):
    (tree / genomes.PAGE).unlink()
    code, err = check(tree, capsys)
    assert code == 1
    assert "line 1" in err


def crlf_line_ends(page: bytes) -> bytes:
    return page.replace(b"\n", b"\r\n")


def no_last_line_end(page: bytes) -> bytes:
    return page[:-1]


def a_byte_that_is_not_utf8(page: bytes) -> bytes:
    return page.replace(b"# Supported genomes", b"# Supported \xff genomes", 1)


@pytest.mark.parametrize(
    ("change", "words"),
    [
        (crlf_line_ends, "line 1: the page's line end is not the one the tables give"),
        (no_last_line_end, "the page's line end is not the one the tables give"),
        (a_byte_that_is_not_utf8, "line 1: the page has `# Supported \\xff genomes`, the tables give `# Supported genomes`"),
    ],
)
def test_a_page_of_other_bytes_fails_the_check_in_one_line(tree, capsys, change, words):
    page = tree / genomes.PAGE
    changed = change(page.read_bytes())
    page.write_bytes(changed)
    code, err = check(tree, capsys)
    assert code == 1
    assert err.startswith("ERROR: [render_supported_genomes] docs/supported-genomes.md is not what the tables give: line ")
    assert err.count("\n") == 1 and "\r" not in err
    assert words in err
    assert page.read_bytes() == changed


def test_the_check_names_the_last_line_of_a_page_without_its_line_end(tree, capsys):
    page = tree / genomes.PAGE
    held = page.read_bytes()
    lines = held.count(b"\n")
    page.write_bytes(held[:-1])
    assert f"line {lines}: the page's line end" in check(tree, capsys)[1]


def test_a_page_that_cannot_be_written_is_refused_in_one_line(tree, capsys):
    page = tree / genomes.PAGE
    page.unlink()
    page.mkdir()
    code = genomes.main(["--repo", str(tree)])
    err = capsys.readouterr().err
    assert code == 1
    assert err.startswith("ERROR: [render_supported_genomes] ") and err.count("\n") == 1
    assert page.is_dir()


def test_writing_restores_a_page_the_check_refused(tree, capsys):
    add_row(tree, "genomes.tsv", "plants")
    assert check(tree, capsys)[0] == 1
    assert genomes.main(["--repo", str(tree)]) == 0
    assert check(tree, capsys) == (0, "")
    assert "| `zz_made_up_genome` | MadeUp1, `ensembl` | `exact` | - |" in (tree / genomes.PAGE).read_text(encoding="utf-8")


def set_not_listed(tree: Path) -> None:
    rewrite(tree, "sets.tsv", lambda lines: [line for line in lines if not line.startswith("plants\t")])


def set_counted_otherwise(tree: Path) -> None:
    recount(tree, "plants", 1)


def outcome_not_declared(tree: Path) -> None:
    """Give the first row of `genomes.tsv` the outcome `differs:unclassified`, with no digest, and take that
    outcome's line out of `outcomes.tsv`: the copy holds a row whose outcome no line declares, whatever outcomes the
    committed rows read."""

    def change(lines: list[str]) -> list[str]:
        cells = lines[1].split("\t")
        cells[7], cells[8] = "differs:unclassified", "-"
        return [lines[0], "\t".join(cells), *lines[2:]]

    rewrite(tree, "genomes.tsv", change)
    rewrite(tree, "outcomes.tsv", lambda lines: [line for line in lines if not line.startswith("differs:unclassified\t")])


def header_of_other_columns(tree: Path) -> None:
    rewrite(tree, "genomes.tsv", lambda lines: [lines[0].replace("outcome", "result"), *lines[1:]])


def row_of_another_width(tree: Path) -> None:
    rewrite(tree, "genomes.tsv", lambda lines: [*lines, "one\tcell\tshort"])


def set_of_no_known_kind(tree: Path) -> None:
    add_row(tree, "genomes.tsv", "archaea")
    rewrite(tree, "sets.tsv", lambda lines: [*lines, "archaea\t63\t63\t1"])


def table_file_missing(tree: Path) -> None:
    (tree / genomes.TABLE_DIR / "bacteria.tsv.gz").unlink()


def set_without_a_row(tree: Path) -> None:
    rewrite(tree, "genomes.tsv", lambda lines: [line for line in lines if line.split("\t")[3] != "pangenome"])


def row_of_another_cache_version(tree: Path) -> None:
    def change(lines: list[str]) -> list[str]:
        out = []
        for line in lines:
            cells = line.split("\t")
            if cells[0] == "zea_mays":
                cells[4] = "62"
            out.append("\t".join(cells))
        return out

    rewrite(tree, "genomes.tsv", change)


def side_table_blank_line(tree: Path) -> None:
    rewrite(tree, "sets.tsv", lambda lines: [*lines, ""])


def side_table_short_row(tree: Path) -> None:
    rewrite(tree, "outcomes.tsv", lambda lines: [*lines, "differs:made-up"])


def set_listed_twice(tree: Path) -> None:
    rewrite(tree, "sets.tsv", lambda lines: [*lines, "plants\t64\t64\t267"])


def gzip_cut_short(tree: Path) -> None:
    path = tree / genomes.TABLE_DIR / "bacteria.tsv.gz"
    held = path.read_bytes()
    path.write_bytes(held[: len(held) // 2])


def gzip_damaged(tree: Path) -> None:
    """A gzip header before bytes that open no deflate block."""
    (tree / genomes.TABLE_DIR / "bacteria.tsv.gz").write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff" + b"\xff" * 32)


def table_not_utf8(tree: Path) -> None:
    with open(tree / genomes.TABLE_DIR / "genomes.tsv", "ab") as handle:
        handle.write(b"\xff\n")


def corpora_of(tree: Path, manifest: str) -> None:
    """Replace the linked corpora by one made-up corpus whose manifest is the text given."""
    (tree / "tests" / "golden").unlink()
    corpus = tree / genomes.CORPORA / "made-up"
    corpus.mkdir(parents=True)
    (corpus / "manifest.json").write_text(manifest, encoding="utf-8")


def manifest_not_json(tree: Path) -> None:
    corpora_of(tree, "{")


def manifest_no_object(tree: Path) -> None:
    corpora_of(tree, "[]")


def manifest_without_assembly(tree: Path) -> None:
    corpora_of(tree, '{"assembly": 7, "records": []}')


def manifest_species_no_text(tree: Path) -> None:
    corpora_of(tree, '{"assembly": "MadeUp1", "species": ["made_up"], "records": []}')


def manifest_flags_no_list(tree: Path) -> None:
    corpora_of(tree, '{"assembly": "MadeUp1", "flags": 7, "records": []}')


@pytest.mark.parametrize(
    ("damage", "words"),
    [
        (set_not_listed, "set plants: tests/genomes/sets.tsv does not list it"),
        (set_counted_otherwise, "set plants: 267 rows where tests/genomes/sets.tsv counts 268"),
        (outcome_not_declared, "outcome differs:unclassified: tests/genomes/outcomes.tsv does not declare it"),
        (header_of_other_columns, "genomes.tsv: the header is not the columns"),
        (row_of_another_width, "3 cells where the header has 9"),
        (set_of_no_known_kind, "set archaea: no title in SPECIES_SETS"),
        (table_file_missing, "bacteria.tsv.gz: no such file"),
        (set_without_a_row, "set pangenome: no row where tests/genomes/sets.tsv counts "),
        (row_of_another_cache_version, "cache version 62 where tests/genomes/sets.tsv gives set plants 63"),
        (side_table_blank_line, ": 0 cells where the header has 4"),
        (side_table_short_row, ": 1 cells where the header has 2"),
        (set_listed_twice, ": a second row of set plants"),
        (gzip_cut_short, "bacteria.tsv.gz: "),
        (gzip_damaged, "bacteria.tsv.gz: Error "),
        (table_not_utf8, "genomes.tsv: 'utf-8' codec can't decode byte 0xff"),
        (manifest_not_json, "made-up/manifest.json: "),
        (manifest_no_object, "made-up/manifest.json: no JSON object"),
        (manifest_without_assembly, "made-up/manifest.json: names no assembly"),
        (manifest_species_no_text, "made-up/manifest.json: names no species"),
        (manifest_flags_no_list, "made-up/manifest.json: `flags` is no list"),
    ],
)
def test_a_table_the_writer_does_not_render_is_refused_in_one_line(tree, capsys, damage, words):
    before = (tree / genomes.PAGE).read_bytes()
    damage(tree)
    code = genomes.main(["--repo", str(tree)])
    err = capsys.readouterr().err
    assert code == 1
    assert err.startswith("ERROR: [render_supported_genomes] ") and err.count("\n") == 1
    assert words in err
    assert (tree / genomes.PAGE).read_bytes() == before


def test_every_cache_of_a_listed_set_has_one_page_row_and_no_other():
    listed = [row for row in table_rows("genomes.tsv") if not row["set"].endswith("_collection")]
    expected = Counter((row["genome"], f"{row['assembly']}, {row['transcripts']}", row["outcome"]) for row in listed)
    sections = page_sections()
    found = Counter(
        (plain(cells[0]), plain(cells[1]), plain(cells[2]))
        for rows in sections.values()
        for cells in rows
        if len(cells) == 4 and (cells[3] == "-" or cells[3].startswith(("yes", "[yes]")))
    )
    assert found == expected
    assert sum(found.values()) == len(listed) > 1000


def test_a_corpus_cell_says_yes_for_the_caches_the_corpora_annotate_and_for_no_other():
    corpora = set(genomes.read_corpora(ROOT))
    annotated = {
        (row["genome"], f"{row['assembly']}, {row['transcripts']}")
        for row in table_rows("genomes.tsv")
        if not row["set"].endswith("_collection") and (row["genome"], row["assembly"], row["transcripts"]) in corpora
    }
    stated = {
        (plain(cells[0]), plain(cells[1]))
        for rows in page_sections().values()
        for cells in rows
        if len(cells) == 4 and cells[3].startswith(("yes", "[yes]"))
    }
    assert stated == annotated
    assert len(annotated) > 30


def test_a_median_between_two_counts_carries_the_separators_too():
    assert genomes.median_text([1000, 2000]) == "1,500"
    assert genomes.median_text([30, 31]) == "30.5"
    assert genomes.median_text([1000, 2001]) == "1,500.5"


def test_every_collection_has_one_line_that_holds_its_count():
    by_set = Counter(row["set"] for row in all_rows() if row["set"].endswith("_collection"))
    exact = Counter(row["set"] for row in all_rows() if row["set"].endswith("_collection") and row["outcome"] == "exact")
    lines = {}
    for heading, rows in page_sections().items():
        if heading.endswith("collections"):
            for cells in rows:
                assert plain(cells[0]) not in lines
                lines[plain(cells[0])] = cells
    assert {name: count(cells[1]) for name, cells in lines.items()} == dict(by_set)
    for name, cells in lines.items():
        assert cells[2].split(";")[0] == f"exact {exact[name]:,}", name
    assert sum(by_set.values()) + len(table_rows("genomes.tsv")) - sum(
        1 for row in table_rows("genomes.tsv") if row["set"].endswith("_collection")
    ) == len(all_rows())


def test_the_line_of_four_counts_is_exact_of_compared_and_leaves_excluded_rows_out(tree):
    rows = all_rows()

    def figures(keep) -> tuple[int, int]:
        compared = [row for row in rows if keep(row["set"]) and not row["outcome"].startswith("excluded:")]
        return sum(row["outcome"] == "exact" for row in compared), len(compared)

    species = {"vertebrates", "metazoa", "plants", "fungi", "protists"}
    expected = [
        ("species caches", *figures(lambda name: name in species)),
        ("bacterial and archaeal genomes", *figures(lambda name: name.startswith("bacteria_"))),
        (
            "fungal and protist collection genomes",
            *figures(lambda name: name.endswith("_collection") and name.startswith(("fungi_", "protists_"))),
        ),
        ("human pangenome assemblies", *figures(lambda name: name == "pangenome")),
    ]
    assert four_counts(PAGE.read_text(encoding="utf-8")) == expected
    excluded = sum(row["outcome"].startswith("excluded:") for row in rows)
    assert sum(compared for _, _, compared in expected) == len(rows) - excluded
    exclude_one(tree, "plants")
    assert genomes.main(["--repo", str(tree)]) == 0
    (label, exact, compared), others = expected[0], expected[1:]
    written = (tree / genomes.PAGE).read_text(encoding="utf-8")
    assert four_counts(written) == [(label, exact - 1, compared - 1), *others]


def test_the_sentences_on_what_was_compared_set_excluded_caches_apart_only_where_a_row_is_excluded(tree):
    each = "Both engines annotated a small generated input for each"
    once = f"was compared once on the generated records its row counts. {DEFAULT_FLAGS_SENTENCE} {GIVEN_REFERENCE_SENTENCE}"

    def plants(page: str) -> str:
        return next(line for line in page.splitlines() if line.startswith(f"| {genomes.SPECIES_SETS['plants'][0]} "))

    compare_every_row(tree)
    without = genomes.render(tree)
    assert f"{each}.\n" in without and f" Every cache on this page {once}" in without
    assert EXCLUDED_SENTENCE not in without and "merged transcripts.\n" in without
    assert plants(without).endswith("| 0 |")
    genome, assembly, transcripts = exclude_one(tree, "plants")[:3]
    with_one = genomes.render(tree)
    assert f"{each} compared.\n" in with_one and f" Every cache on this page that was not excluded {once}" in with_one
    assert f"merged transcripts. {EXCLUDED_SENTENCE}\n" in with_one
    assert 80 <= len(with_one.splitlines()[2]) <= 160
    assert f"| `{EXCLUDED_OUTCOME}` | 1 | {EXCLUDED_MEANING} |\n" in with_one
    assert plants(with_one).endswith(f"| `{EXCLUDED_OUTCOME}` 1 |")
    assert f"| `{genome}` | {assembly}, `{transcripts}` | `{EXCLUDED_OUTCOME}` | " in with_one


def edited_human_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """The rows of human's own RefSeq and merged caches, the caches Ensembl edits from an alignment file."""
    return [row for row in rows if row["genome"] == genomes.HUMAN and row["transcripts"] in ("merged", "refseq")]


def test_the_page_names_the_option_both_engines_were_given_on_the_caches_edited_from_an_alignment_file():
    text = PAGE.read_text(encoding="utf-8")
    paragraph = text.splitlines()[4]
    assert f" {DEFAULT_FLAGS_SENTENCE} {GIVEN_REFERENCE_SENTENCE} No genome other than human is timed." in paragraph
    assert "under default flags, on the generated records" not in text
    assert text.count("--use_given_ref") == 1
    # The caches the sentence names are rows of the page, each compared.
    named = edited_human_rows(table_rows("genomes.tsv"))
    assert sorted((row["assembly"], row["transcripts"]) for row in named) == [
        ("GRCh37", "merged"),
        ("GRCh37", "refseq"),
        ("GRCh38", "merged"),
        ("GRCh38", "refseq"),
    ]
    assert not any(row["outcome"].startswith("excluded:") for row in named)
    for row in named:
        assert f"| `{row['genome']}` | {row['assembly']}, `{row['transcripts']}` | `{row['outcome']}` | " in text
    # The job gives both engines the option on a cache that names an alignment file and on no other.
    job = (ROOT / "scripts" / "genomes" / "check_genome.sh").read_text(encoding="utf-8")
    sets_the_option = (
        'REFERENCE_FLAG=""\n'
        'case "$(alignment_file_of "$JSON_CACHE/info.json")" in\n'
        "'' | 0) ;;\n"
        '*) REFERENCE_FLAG="--use_given_ref" ;;\n'
        "esac\n"
    )
    assert job.count(sets_the_option) == 1 and job.count("REFERENCE_FLAG=") == 2
    assert job.count('ENSEMBL_VEP_COMMAND="$ENSEMBL_VEP_COMMAND $REFERENCE_FLAG"') == 1
    assert job.count("$SET_FLAG $REFERENCE_FLAG") == 1
    # The page the sentence links says which caches are edited so.
    checks = (ROOT / "docs" / "genome-checks.md").read_text(encoding="utf-8")
    assert "it gives both engines `--use_given_ref` on such a cache and on no other." in checks
    assert "Ensembl's human RefSeq and merged caches are edited so." in checks


def test_the_option_is_named_only_where_a_cache_edited_from_an_alignment_file_was_compared(tree):
    timed = "No genome other than human is timed."
    table = tree / genomes.TABLE_DIR / "genomes.tsv"

    def exclude(keep: str | None):
        """Give every row of human's RefSeq and merged caches but the first of the transcript set `keep` a made-up
        excluded outcome."""
        seen = []

        def change(lines: list[str]) -> list[str]:
            out = []
            for line in lines:
                cells = line.split("\t")
                if cells[0] == genomes.HUMAN and cells[2] in ("merged", "refseq"):
                    seen.append(cells[2])
                    if cells[2] != keep or seen.count(keep) > 1:
                        cells[7], cells[8] = EXCLUDED_OUTCOME, "-"
                out.append("\t".join(cells))
            return out

        return change

    rewrite(tree, "outcomes.tsv", lambda lines: [*lines, f"{EXCLUDED_OUTCOME}\t{EXCLUDED_MEANING}"])
    held = table.read_bytes()
    for keep in ("merged", "refseq"):
        table.write_bytes(held)
        rewrite(tree, "genomes.tsv", exclude(keep))
        one_compared = genomes.render(tree)
        assert f" {DEFAULT_FLAGS_SENTENCE} {GIVEN_REFERENCE_SENTENCE} {timed}" in one_compared, keep
        assert [row["transcripts"] for row in edited_human_rows(genomes.compared_rows(genomes.Tables(tree).rows))] == [keep]
    rewrite(tree, "genomes.tsv", exclude(None))
    none_compared = genomes.render(tree)
    assert f" {DEFAULT_FLAGS_SENTENCE} {timed}" in none_compared
    assert "--use_given_ref" not in none_compared and "alignment file" not in none_compared
    assert len(edited_human_rows(genomes.compared_rows(genomes.Tables(tree).rows))) == 0


def test_the_counts_mode_prints_the_figures_the_page_states(capsys):
    before = PAGE.read_bytes()
    assert genomes.main(["--repo", str(ROOT), "--counts"]) == 0
    printed = [line.split("\t") for line in capsys.readouterr().out.splitlines()]
    text = PAGE.read_text(encoding="utf-8")
    stated = next(line for line in text.splitlines() if line.startswith("Exact of compared: "))
    assert [len(cells) for cells in printed] == [3, 3, 3, 3, 2]
    for label, exact, compared in printed[:4]:
        assert f" {label} {int(exact):,} of {int(compared):,}" in stated
    assert printed[4][0] == "species with a corpus, human apart"
    assert f"The {int(printed[4][1]):,} other species with a committed corpus" in text
    assert PAGE.read_bytes() == before


def test_the_set_summary_counts_every_cache_once_by_outcome():
    rows = all_rows()
    summary = page_sections()["Sets"]
    assert sum(count(re.search(r"\(([\d,]+)\)$", cells[1])[1]) for cells in summary) == len(rows)
    by_outcome: Counter[str] = Counter()
    for cells in summary:
        for part in filter(None, (part.strip() for part in f"{cells[2]};{cells[3]}".split(";"))):
            if part != "0":
                name, figure = plain(part).rsplit(" ", 1)
                by_outcome[name] += count(figure)
    assert by_outcome == Counter(row["outcome"] for row in rows)
    stated = {plain(cells[0]): count(cells[1]) for cells in page_sections()["Outcomes"]}
    assert stated == dict(by_outcome)


def test_the_records_table_gives_each_sets_least_and_greatest():
    rows = [row for row in all_rows() if not row["outcome"].startswith("excluded:")]
    records = page_sections()["Records compared"]
    assert min(count(cells[1]) for cells in records) == min(int(row["records"]) for row in rows)
    assert max(count(cells[3]) for cells in records) == max(int(row["records"]) for row in rows)
    bacteria = [int(row["records"]) for row in table_rows("bacteria.tsv.gz")]
    line = next(cells for cells in records if cells[0].startswith("Bacterial"))
    assert (count(line[1]), count(line[3])) == (min(bacteria), max(bacteria))
    under_ten = sum(value < 10 for value in bacteria)
    assert f"{under_ten:,} of {len(bacteria):,} compared caches had fewer than ten records." in PAGE.read_text(
        encoding="utf-8"
    )


def test_the_title_and_the_opening_paragraph_fit_the_site():
    lines = PAGE.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# Supported genomes" and lines[1] == ""
    assert 80 <= len(lines[2]) <= 160
    assert lines[3] == ""


def species_page_corpus_rows() -> list[str]:
    return re.findall(r"^\|.*\| `tests/golden/116/([^`]+)` \|$", SPECIES_PAGE.read_text(encoding="utf-8"), re.MULTILINE)


def test_the_species_page_tables_every_corpus_of_another_species_once():
    corpora = genomes.read_corpora(ROOT)
    others = {name for (genome, _, _), names in corpora.items() if not genomes.is_human(genome) for name in names}
    assemblies = {
        name
        for (genome, _, _), names in corpora.items()
        if genomes.is_human(genome) and genome != genomes.HUMAN
        for name in names
    }
    tabled = species_page_corpus_rows()
    assert len(tabled) == len(set(tabled))
    assert others <= set(tabled) <= others | assemblies
    assert len(others) > 30


def test_the_species_page_counts_the_species_the_corpora_give():
    corpora = genomes.read_corpora(ROOT)
    species = {genome for genome, _, _ in corpora if not genomes.is_human(genome)}
    for stated in (SPECIES_PAGE, ROOT / "README.md"):
        text = stated.read_text(encoding="utf-8")
        assert f"and on {len(species)} other species" in text, stated
        assert f"Each of the {len(species)} has a golden corpus of its own" in text, stated
    assert f"The {len(species)} genomes exercise tables " in SPECIES_PAGE.read_text(encoding="utf-8")
    assert f"The {len(species)} other species with a committed corpus" in PAGE.read_text(encoding="utf-8")


def test_every_corpus_link_lands_on_a_heading_of_the_species_page():
    headings = {
        genomes.heading_anchor(line[3:])
        for line in SPECIES_PAGE.read_text(encoding="utf-8").splitlines()
        if line.startswith("## ")
    }
    linked = set(re.findall(r"\[yes\]\(cache-species\.md#([^)]+)\)", PAGE.read_text(encoding="utf-8")))
    assert linked and linked <= headings
    assert genomes.heading_anchor(genomes.SPECIES_TABLE_HEADING) in headings


def test_a_pangenome_assembly_with_a_corpus_is_human_and_links_its_own_section(tree):
    annotated = set(genomes.read_corpora(ROOT))
    table_row = next(
        row
        for row in table_rows("genomes.tsv")
        if row["set"] == "pangenome" and (row["genome"], row["assembly"], row["transcripts"]) not in annotated
    )
    before = genomes.render(tree)
    counts = genomes.counts_report(tree)
    golden = tree / "tests" / "golden"
    golden.unlink()
    (tree / genomes.CORPORA).mkdir(parents=True)
    for corpus in sorted((ROOT / genomes.CORPORA).iterdir()):
        (tree / genomes.CORPORA / corpus.name).symlink_to(corpus, target_is_directory=True)
    made_up = tree / genomes.CORPORA / "made-up-assembly"
    made_up.mkdir()
    manifest = {"species": table_row["genome"], "assembly": table_row["assembly"], "records": [{}]}
    (made_up / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    assert genomes.is_human(table_row["genome"]) and table_row["genome"] != genomes.HUMAN
    after = genomes.render(tree)
    assert genomes.counts_report(tree) == counts
    anchor = genomes.heading_anchor(genomes.PANGENOME_HEADING)
    cells = f"| `{table_row['genome']}` | {table_row['assembly']}, `ensembl` | `{table_row['outcome']}` |"
    assert f"{cells} - |" in before and f"{cells} [yes](cache-species.md#{anchor}) |" in after
    changed = [(old, new) for old, new in zip(before.splitlines(), after.splitlines()) if old != new]
    assert len(changed) == 1
    species_page = tree / genomes.SPECIES_PAGE
    text = species_page.read_text(encoding="utf-8")
    assert f"## {genomes.PANGENOME_HEADING}\n" in text
    species_page.write_text(text.replace(f"## {genomes.PANGENOME_HEADING}\n", "## Another heading\n"), encoding="utf-8")
    assert f"{cells} yes |" in genomes.render(tree)

