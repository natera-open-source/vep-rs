"""Tests for parquet_to_vep_tab.py: the inverse of the writer's cell contract per type family,
on files written here under the writer's schema; the canonical form; every input shape; and
the golden corpora round-tripped through the `vep` binary when one is built."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pa = pytest.importorskip("pyarrow")
import pyarrow.parquet as pq  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parquet_to_vep_tab as adapter  # noqa: E402

REPO = HERE.parents[1]
CORPORA = REPO / "tests" / "golden" / "116"
ADAPTER = HERE / "parquet_to_vep_tab.py"

KEY_FIELDS = [
    ("contig_index", pa.int32()),
    ("chrom", pa.string()),
    ("pos", pa.int64()),
    ("end", pa.int64()),
    ("ref", pa.string()),
    ("alt", pa.string()),
]
POSITIONS = pa.struct([("start", pa.int32()), ("end", pa.int32()), ("total", pa.int32())])
NULL = adapter.NULL_TEXT


def positions(start, end, total):
    return {"start": start, "end": end, "total": total}


# One tab column per type family and spelling: (column, Arrow type, stored values, the tab
# cells those values print). A column shorter than the longest is NULL on the other rows.
CELLS = [
    ("Feature", pa.string(), ["ENST00000456328", "", None], ["ENST00000456328", "", NULL]),
    ("Allele", pa.string(), ["-", "A"], ["-", "A"]),
    ("Location", pa.string(), ["21:25000100-25000101"], ["21:25000100-25000101"]),
    ("STRAND", pa.int32(), [-1, 1, None], ["-1", "1", NULL]),
    ("DISTANCE", pa.int32(), [4821], ["4821"]),
    ("MAX_AF", pa.float64(), [0.095, 1.0, 9.8e-05, 0.0, None], ["0.095", "1", "0.000098", "0", NULL]),
    ("MOTIF_SCORE_CHANGE", pa.float64(), [-0.012, 0.1], ["-0.012", "0.1"]),
    ("CANONICAL", pa.bool_(), [True, False, None], ["YES", NULL, NULL]),
    ("PICK", pa.bool_(), [True, False], ["1", NULL]),
    ("GENCODE_PRIMARY", pa.bool_(), [True, False], ["1", NULL]),
    ("HIGH_INF_POS", pa.bool_(), [True, False, None], ["Y", "N", NULL]),
    ("AF", pa.list_(pa.float64()), [[0.0002], [0.2766, 0.095], [None, 0.1]], ["0.0002", "0.2766,0.095", "-,0.1"]),
    ("gnomADe_AF", pa.list_(pa.float64()), [[0.095, 0.1]], ["0.095,0.1"]),
    ("gnomADg_NFE_AF", pa.list_(pa.float64()), [[0.0], None], ["0", NULL]),
    ("SOMATIC", pa.list_(pa.bool_()), [[False, True], [None]], ["0,1", NULL]),
    ("PHENO", pa.list_(pa.bool_()), [[True]], ["1"]),
    ("MANE", pa.list_(pa.string()), [["MANE_Select", "MANE_Plus_Clinical"]], ["MANE_Select,MANE_Plus_Clinical"]),
    (
        "Consequence",
        pa.list_(pa.string()),
        [["missense_variant", "splice_region_variant"], ["intergenic_variant"]],
        ["missense_variant,splice_region_variant", "intergenic_variant"],
    ),
    ("CLIN_SIG", pa.list_(pa.string()), [["benign;likely_benign", "pathogenic"]], ["benign;likely_benign,pathogenic"]),
    ("PUBMED", pa.list_(pa.string()), [["12345", "67890"], [None, "1"]], ["12345,67890", "-,1"]),
    ("EXON", POSITIONS, [positions(3, 3, 12), positions(3, 4, 12)], ["3/12", "3-4/12"]),
    ("INTRON", POSITIONS, [positions(7, 7, 11)], ["7/11"]),
    (
        "cDNA_position",
        POSITIONS,
        [positions(44, 44, None), positions(44, 45, None), positions(44, 44, 1200), positions(None, None, 1200), positions(44, 45, 1200)],
        ["44", "44-45", "44/1200", "-/1200", "44-45/1200"],
    ),
    ("CDS_position", POSITIONS, [positions(146, None, None)], ["146-?"]),
    ("Protein_position", POSITIONS, [positions(None, 52, None), positions(49, None, 300), positions(None, None, None)], ["?-52", "49-?/300", NULL]),
    # A plugin column is a string; a BOOLEAN column of no listed spelling prints a flag.
    ("CADD_PHRED", pa.string(), ["23.4"], ["23.4"]),
    ("FLAGGED", pa.bool_(), [True, False], ["1", "0"]),
]
ROWS = max(len(values) for _, _, values, _ in CELLS)


def table(columns, rows: int, chrom: str = "21") -> pa.Table:
    """A table of the writer's schema: the key columns, then ``columns`` as (name, type,
    values) padded with NULL to ``rows`` rows."""
    arrays = [
        pa.array(list(range(rows)), pa.int32()),
        pa.array([chrom] * rows, pa.string()),
        pa.array([100 + i for i in range(rows)], pa.int64()),
        pa.array([100 + i for i in range(rows)], pa.int64()),
        pa.array(["A"] * rows, pa.string()),
        pa.array(["G"] * rows, pa.string()),
    ]
    names = [name for name, _ in KEY_FIELDS]
    for name, kind, values in columns:
        arrays.append(pa.array(list(values) + [None] * (rows - len(values)), kind))
        names.append(name)
    return pa.table(arrays, names=names)


def write(path: Path, columns, rows: int, chrom: str = "21") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table(columns, rows, chrom).replace_schema_metadata({"vep_rs_version": "0"}), path)


def run(*args: str) -> str:
    proc = subprocess.run(
        [sys.executable, str(ADAPTER), *args], capture_output=True, text=True, check=False
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def fixture(tmp_path: Path) -> Path:
    path = tmp_path / "cells.parquet"
    write(path, [(name, kind, values) for name, kind, values, _ in CELLS], ROWS)
    return path


def test_every_typed_cell_prints_its_tab_text(tmp_path: Path):
    lines = run(str(fixture(tmp_path))).splitlines()
    assert len(lines) == 1 + ROWS
    for row, line in enumerate(lines[1:]):
        cells = line.split("\t")
        assert len(cells) == len(CELLS)
        for (column, _, _, texts), cell in zip(CELLS, cells):
            expected = texts[row] if row < len(texts) else NULL
            assert cell == expected, f"{column} row {row}"


def test_the_header_is_the_tab_columns_after_the_key_columns(tmp_path: Path):
    header = run(str(fixture(tmp_path))).splitlines()[0]
    assert header == "#" + "\t".join(name for name, _, _, _ in CELLS)
    assert not any(key in header.split("\t") for key, _ in KEY_FIELDS)


def test_cell_text_prints_each_value_kind():
    assert adapter.cell_text("Feature", None) == NULL
    assert adapter.cell_text("Allele", "-") == "-"
    assert adapter.cell_text("Feature", "") == ""
    assert adapter.cell_text("STRAND", -1) == "-1"
    assert adapter.cell_text("MAX_AF", 0.095) == "0.095"
    assert adapter.cell_text("CANONICAL", True) == "YES"
    assert adapter.cell_text("CANONICAL", False) == NULL
    assert adapter.cell_text("HIGH_INF_POS", False) == "N"
    assert adapter.cell_text("SOMATIC", [False, True, None]) == "0,1,-"
    assert adapter.cell_text("AF", [0.0002, None]) == "0.0002,-"
    assert adapter.cell_text("Consequence", ["a", None]) == "a,-"
    assert adapter.cell_text("EXON", positions(3, 4, 12)) == "3-4/12"
    assert adapter.cell_text("cDNA_position", positions(None, None, None)) == NULL


def test_shortest_decimal_prints_the_writers_text():
    cases = [
        (0.095, "0.095"),
        (1.0, "1"),
        (9.8e-05, "0.000098"),
        (0.0, "0"),
        (-0.0, "-0"),
        (1e16, "10000000000000000"),
        (1e21, "1000000000000000000000"),
        (0.1 + 0.2, "0.30000000000000004"),
        (5e-324, "0." + "0" * 323 + "5"),
        (float("nan"), "NaN"),
        (float("inf"), "inf"),
        (float("-inf"), "-inf"),
        # Two shortest decimals equally close: the one of larger magnitude.
        (-126430982921026.62, "-126430982921026.63"),
        (126430982921026.62, "126430982921026.63"),
    ]
    for value, text in cases:
        assert adapter.shortest_decimal(value) == text, repr(value)
        if text not in ("NaN", "inf", "-inf"):
            assert float(text) == value, text


def test_parse_double_follows_the_contracts_number_grammar():
    for token in ("0.0950", "1", "-0.012", "9.8e-05", "1E5", ".5", "3.", "+1", "inf", "-Infinity", "nan"):
        assert adapter.parse_double(token) is not None, token
    for token in ("n/a", "1_000", " 1", "0x10", "1e", "-", "", "1,2"):
        assert adapter.parse_double(token) is None, token


def test_canonical_tab_cell_rewrites_a_double_and_an_int32_cell():
    assert adapter.canonical_tab_cell("MAX_AF", "0.0950") == "0.095"
    assert adapter.canonical_tab_cell("MAX_AF", "1.0") == "1"
    assert adapter.canonical_tab_cell("MAX_AF", "9.8e-05") == "0.000098"
    assert adapter.canonical_tab_cell("MAX_AF", NULL) == NULL
    assert adapter.canonical_tab_cell("MAX_AF", "n/a") == "n/a"
    assert adapter.canonical_tab_cell("MOTIF_SCORE_CHANGE", "0.100") == "0.1"
    assert adapter.canonical_tab_cell("gnomADe_AF", "0.0950,0.1,1") == "0.095,0.1,1"
    assert adapter.canonical_tab_cell("Feature", "0.0950") == "0.0950"
    # An INT32 cell is stored as its integer, and a text that is no INT32 as NULL.
    assert adapter.canonical_tab_cell("TSL", "007") == "7"
    assert adapter.canonical_tab_cell("TSL", "00") == "0"
    assert adapter.canonical_tab_cell("STRAND", "-1") == "-1"
    assert adapter.canonical_tab_cell("TSL", "2147483648") == NULL
    assert adapter.canonical_tab_cell("TSL", "1.5") == NULL
    assert adapter.canonical_tab_cell("TSL", "0" * 5000 + "7") == "7"
    assert adapter.canonical_tab_cell("TSL", "9" * 5000) == NULL
    assert adapter.canonical_tab_cell("TSL", " 7") == NULL
    assert adapter.canonical_tab_cell("TSL", NULL) == NULL
    assert adapter.canonical_tab_cell("Feature", "007") == "007"
    assert len(adapter.DOUBLE_LIST_COLUMNS) == 28
    assert {"AF", "AFR_AF", "gnomADe_AF", "gnomADe_REMAINING_AF", "gnomADg_AF", "gnomADg_SAS_AF"} <= adapter.DOUBLE_LIST_COLUMNS
    assert "MAX_AF_POPS" not in adapter.DOUBLE_LIST_COLUMNS


def test_the_parts_of_a_directory_come_out_in_path_order(tmp_path: Path):
    out = tmp_path / "out.parquet"
    columns = [("Feature", pa.string(), ["second"])]
    write(out / "part-00001.parquet", columns, 1)
    write(out / "part-00000.parquet", [("Feature", pa.string(), ["first"])], 1)
    (out / "_vep_rs_manifest.json").write_text("{}", encoding="utf-8")
    assert run(str(out)).splitlines() == ["#Feature", "first", "second"]


def test_a_hive_directory_is_read_file_by_file_with_chrom_from_the_files(tmp_path: Path):
    out = tmp_path / "out.parquet"
    write(out / "chrom=21" / "part-00000.parquet", [("Location", pa.string(), ["21:100"])], 1, chrom="21")
    write(out / "chrom=X" / "part-00000.parquet", [("Location", pa.string(), ["X:100", "X:101"])], 2, chrom="X")
    assert run(str(out)).splitlines() == ["#Location", "21:100", "X:100", "X:101"]
    # The stored `chrom` is what a reader of the files gets, whatever the directory says.
    assert pq.read_table(out / "chrom=X" / "part-00000.parquet").column("chrom").to_pylist() == ["X", "X"]


def test_canonical_sorts_the_rows_and_reprints_every_double_cell(tmp_path: Path):
    parquet = tmp_path / "out.parquet"
    columns = [
        ("Feature", pa.string(), ["b", "a", "c"]),
        ("MAX_AF", pa.float64(), [0.095, 9.8e-05, None]),
        ("AF", pa.list_(pa.float64()), [[0.2766, 0.095], [0.1], None]),
    ]
    write(parquet, columns, 3)
    header = "#Feature\tMAX_AF\tAF"
    canonical = run("--canonical", str(parquet))
    assert canonical.splitlines() == [header, "a\t0.000098\t0.1", "b\t0.095\t0.2766,0.095", "c\t-\t-"]
    plain = run(str(parquet)).splitlines()
    assert sorted(plain[1:]) == canonical.splitlines()[1:]

    tab = tmp_path / "out.tab"
    tab.write_text(
        "## ENSEMBL VARIANT EFFECT PREDICTOR\n## Output produced at ...\n"
        f"{header}\nc\t-\t-\nb\t0.0950\t0.2766,0.0950\na\t9.8e-05\t0.100\n",
        encoding="utf-8",
    )
    assert run("--canonical", str(tab)) == canonical
    gz = tmp_path / "out.tab.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as fh:
        fh.write(tab.read_text(encoding="utf-8"))
    assert run("--canonical", str(gz)) == canonical

    again = tmp_path / "again.tab"
    again.write_text(canonical, encoding="utf-8")
    assert run("--canonical", str(again)) == canonical

    # A line ends at a newline alone: a cell holding another separator character stays
    # one row, and a CRLF file reads like an LF one.
    separators = tmp_path / "separators.tab"
    separators.write_text(f"{header}\r\nb c\t0.0950\t0.1\r\n", encoding="utf-8")
    assert run("--canonical", str(separators)) == f"{header}\nb c\t0.095\t0.1\n"


def test_the_output_option_writes_a_file(tmp_path: Path):
    parquet = fixture(tmp_path)
    out = tmp_path / "out.tab"
    assert run("--input", str(parquet), "--output", str(out)) == ""
    assert out.read_text(encoding="utf-8") == run(str(parquet))


def test_errors_are_reported_on_stderr_with_exit_2(tmp_path: Path, capsys):
    def fails(*args: str) -> str:
        assert adapter.main(list(args)) == 2
        err = capsys.readouterr().err
        assert err.startswith(f"ERROR: [{adapter.COMPONENT}] ")
        return err

    assert "no such file" in fails(str(tmp_path / "absent.parquet"))
    assert "no such file" in fails("--canonical", str(tmp_path / "absent.tab"))
    empty = tmp_path / "empty"
    empty.mkdir()
    assert "no .parquet file" in fails(str(empty))
    mixed = tmp_path / "mixed.parquet"
    write(mixed / "part-00000.parquet", [("Feature", pa.string(), ["a"])], 1)
    write(mixed / "part-00001.parquet", [("Gene", pa.string(), ["g"])], 1)
    assert "tab columns differ" in fails(str(mixed))
    headless = tmp_path / "headless.tab"
    headless.write_text("a\t1\n", encoding="utf-8")
    assert "no # header line" in fails("--canonical", str(headless))
    keys_only = tmp_path / "keys.parquet"
    pq.write_table(table([], 1), keys_only)
    assert "no tab columns" in fails(str(keys_only))
    # A regular file without the Parquet magic is refused before pyarrow opens it; a file
    # with the magic but no footer, a non-UTF-8 tab file and a truncated gzip get the ERROR
    # line too, not a traceback.
    plain = tmp_path / "plain.tab"
    plain.write_text("#Feature\ta\n", encoding="utf-8")
    assert "not a Parquet file" in fails(str(plain))
    corrupt = tmp_path / "corrupt.parquet"
    corrupt.write_bytes(adapter.PARQUET_MAGIC + b"\0" * 8)
    fails(str(corrupt))
    fails("--canonical", str(corrupt))
    binary = tmp_path / "binary.tab"
    binary.write_bytes(b"\xff\xfe\x00")
    fails("--canonical", str(binary))
    truncated = tmp_path / "truncated.tab.gz"
    truncated.write_bytes(gzip.compress(b"#Feature\ta\n" * 100)[:-8])
    fails("--canonical", str(truncated))
    if hasattr(os, "mkfifo"):
        fifo = tmp_path / "fifo"
        os.mkfifo(fifo)
        assert "neither a regular file nor a directory" in fails(str(fifo))
        assert "neither a regular file nor a directory" in fails("--canonical", str(fifo))


def vep_binary() -> Path | None:
    """The built `vep` binary: ``VEP_BINARY``, else the release or debug build under the
    cargo target directory."""
    named = os.environ.get("VEP_BINARY")
    if named:
        return Path(named)
    target = Path(os.environ.get("CARGO_TARGET_DIR", REPO / "target"))
    for profile in ("release", "debug"):
        candidate = target / profile / "vep"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    return None


@pytest.mark.parametrize("corpus", ["GRCh38", "GRCh38-everything", "GRCh38-minor-allele-columns"])
def test_a_golden_corpus_round_trips_through_the_writer(tmp_path: Path, corpus: str):
    vep = vep_binary()
    if vep is None:
        pytest.skip("no vep binary: set VEP_BINARY or build one under the cargo target directory")
    directory = CORPORA / corpus
    if not (directory / "manifest.json").is_file():
        pytest.skip(f"no golden corpus at {directory}")
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    common = [
        str(vep), "--offline", "--json_cache", str(directory / "json_cache"), "--species", "homo_sapiens",
        "--assembly", manifest["assembly"], "--buffer_size", "5000", "--force_overwrite", "--no_stats",
        "--quiet", "-i", str(directory / "variants.vcf"), *(manifest.get("flags") or []),
    ]
    tab = tmp_path / "out.tab"
    parquet = tmp_path / "out.parquet"
    subprocess.run([*common, "--tab", "-o", str(tab)], check=True)
    subprocess.run([*common, "--output_format", "parquet", "-o", str(parquet)], check=True)
    from_parquet = run("--canonical", str(parquet))
    from_tab = run("--canonical", str(tab))
    assert from_parquet == from_tab
    rows = sum(1 for line in tab.read_text(encoding="utf-8").splitlines() if not line.startswith("#"))
    assert rows > 0
    assert len(from_parquet.splitlines()) == 1 + rows
