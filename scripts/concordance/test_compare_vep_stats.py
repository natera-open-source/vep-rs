"""Tests for compare_vep_stats.py: the parser on the committed Ensembl VEP 116.2 statistics
report of tests/golden/116/GRCh38-statistics, the masks, the set rule of the chart
sections, the tolerance, the verdict and the report's shape. The candidate texts are the
committed report with one change each, plus a candidate-shaped report whose volatile rows
and tie orders differ as a second engine's would."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_vep_formats as cvf  # noqa: E402
import compare_vep_stats as cvs  # noqa: E402

CORPUS = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-statistics"
REFERENCE = CORPUS / "expected" / "default.summary.log"
PREDICTIONS = CORPUS / "expected_predictions" / "default.summary.log"

# Stats.pm:806-831 writes the three tabular sections, then the charts of
# generate_chart_data (Stats.pm:505-609) in its order; the fixture run carries no
# prediction flag and no FASTA, so it has neither the two prediction charts nor a
# distribution chart.
FIXTURE_SECTIONS = [
    "VEP run statistics",
    "Data version",
    "General statistics",
    "Variant classes",
    "Consequences (most severe)",
    "Consequences (all)",
    "Coding consequences",
    "Variants by chromosome",
    "Position in protein",
]

# The volatile rows as a second run of another engine writes them: its own command line
# (no `--database 0`, a `--json_cache` directory), UTC times, its own paths.
CANDIDATE_VOLATILE_ROWS = {
    "Annotation sources": "Cache: /home/user/vep-rs/tests/golden/116/GRCh38-statistics/json_cache",
    "Command line options": "vep --assembly GRCh38 --buffer_size 5000 --force_overwrite "
    "--input_file [PATH]/variants.vcf --json_cache [PATH]/json_cache --offline "
    "--output_file [PATH]/default.txt --quiet --species homo_sapiens --stats_text",
    "Start time": "2030-01-01 00:00:00",
    "End time": "2030-01-01 00:00:01",
    "Run time": "1 seconds",
    "Input file": "/home/user/vep-rs/tests/golden/116/GRCh38-statistics/variants.vcf",
    "Output file": "/data/out/GRCh38-statistics-default.txt",
}


@pytest.fixture(scope="module")
def reference() -> str:
    return REFERENCE.read_text(encoding="utf-8")


def with_row(text: str, section: str, key: str, value: str) -> str:
    """The text with the value of row ``key`` of ``[section]`` replaced."""
    out = []
    current = None
    replaced = False
    for line in text.splitlines():
        header = cvs.SECTION_RE.match(line)
        if header:
            current = header.group(1)
        elif current == section and line.split("\t", 1)[0] == key:
            line = f"{key}\t{value}"
            replaced = True
        out.append(line)
    assert replaced, f"no row {key!r} in [{section}]"
    return "\n".join(out) + "\n"


def section_lines(text: str, section: str) -> list[str]:
    lines = []
    current = None
    for line in text.splitlines():
        header = cvs.SECTION_RE.match(line)
        if header:
            current = header.group(1)
        elif current == section and line:
            lines.append(line)
    return lines


def with_section_lines(text: str, section: str, rows: list[str]) -> str:
    """The text with the rows of ``[section]`` replaced by ``rows``."""
    out = []
    current = None
    for line in text.splitlines():
        header = cvs.SECTION_RE.match(line)
        if header:
            current = header.group(1)
            out.append(line)
            if current == section:
                out.extend(rows)
            continue
        if current == section and line:
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def candidate_shaped(text: str) -> str:
    """The reference as a second engine would write it: the volatile rows its own, the
    equal-count ties of [Variant classes] in key order."""
    for key, value in CANDIDATE_VOLATILE_ROWS.items():
        text = with_row(text, cvs.RUN_STATISTICS, key, value)
    classes = section_lines(text, "Variant classes")
    assert classes[1:3] == ["sequence_alteration\t7", "indel\t7"], classes
    return with_section_lines(text, "Variant classes", [classes[0], classes[2], classes[1], *classes[3:]])


def test_sections_parse_in_order(reference):
    sections = cvs.parse_report(reference)
    assert [s.title for s in sections] == FIXTURE_SECTIONS
    assert [s.is_chart for s in sections] == [False] * 3 + [True] * 6
    run_stats = sections[0]
    assert run_stats.ordered_keys() == [
        "VEP version (API)",
        "Annotation sources",
        "Species",
        "Command line options",
        "Start time",
        "End time",
        "Run time",
        "Input file",
        "Output file",
    ]
    assert dict(run_stats.rows)["VEP version (API)"] == " 116.2 (116)"
    assert dict(sections[2].rows)["Lines of input read"] == "340"
    assert dict(sections[3].rows) == {
        "duplication": "1",
        "sequence_alteration": "7",
        "indel": "7",
        "substitution": "9",
        "deletion": "35",
        "insertion": "37",
        "SNV": "243",
    }
    # A row without a tab (a [General statistics] value undefined) keeps a None value;
    # text before the first header is not a statistics report.
    assert cvs.parse_report("[General statistics]\nOverlapped genes\n")[0].rows == [("Overlapped genes", None)]
    with pytest.raises(cvs.ComparatorError, match="precedes the first"):
        cvs.parse_report("Lines of input read\t3\n[General statistics]\n")


def test_masked_rows_never_differ(reference):
    candidate = candidate_shaped(reference)
    assert candidate != reference
    report = cvs.compare_reports(reference, candidate)
    assert report["verdict"] == "green", report["reasons"]
    assert report["differences"] == {}
    assert report["compared"] == {"sections": 9, "rows": 9 + 13 + 7 + 7 + 22 + 26 + 9 + 2 + 10}
    assert set(cvs.MASKED_ROWS) <= set(report["masked_rows"])
    # The path mask leaves the structure of Annotation sources to compare: a second source
    # on one side alone is a difference, the directory alone never is.
    masked = with_row(reference, cvs.RUN_STATISTICS, "Annotation sources", "Cache: [PATH]/116_GRCh38")
    assert cvs.compare_reports(reference, masked)["verdict"] == "green"
    custom = with_row(
        reference, cvs.RUN_STATISTICS, "Annotation sources",
        "Cache: /work/cache/homo_sapiens/116_GRCh38; Custom: /work/corpus/custom/sites.vcf.gz (vcf)",
    )
    report = cvs.compare_reports(reference, custom)
    assert report["verdict"] == "red"
    assert report["differences"]["VEP run statistics"] == [{
        "kind": "value",
        "row": "Annotation sources",
        "expected": "Cache: [PATH]",
        "actual": "Cache: [PATH]; Custom: [PATH] (vcf)",
    }]
    # The rows that are not masked compare verbatim.
    species = with_row(reference, cvs.RUN_STATISTICS, "Species", "mus_musculus")
    assert cvs.compare_reports(reference, species)["differences"]["VEP run statistics"] == [
        {"kind": "value", "row": "Species", "expected": "homo_sapiens", "actual": "mus_musculus"}
    ]
    version = with_row(reference, cvs.RUN_STATISTICS, "VEP version (API)", "116.2 (116)")
    assert cvs.compare_reports(reference, version)["verdict"] == "red"
    # --mask-row sets one more label aside, in whichever section it occurs.
    genes = with_row(reference, cvs.GENERAL_STATISTICS, "Overlapped genes", "97")
    assert cvs.compare_reports(reference, genes)["verdict"] == "red"
    assert cvs.compare_reports(reference, genes, mask_rows=("Overlapped genes",))["verdict"] == "green"


def test_chart_sections_compare_as_sets(reference):
    rows = section_lines(reference, "Consequences (most severe)")
    swapped = with_section_lines(reference, "Consequences (most severe)", [rows[1], rows[0], *rows[2:]])
    assert swapped != reference
    report = cvs.compare_reports(reference, swapped)
    assert report["verdict"] == "green", report["reasons"]
    # Rows of a tabular section keep their order: Stats.pm writes them in a fixed sequence.
    run_rows = section_lines(reference, cvs.RUN_STATISTICS)
    reordered = with_section_lines(reference, cvs.RUN_STATISTICS, [run_rows[2], run_rows[1], run_rows[0], *run_rows[3:]])
    report = cvs.compare_reports(reference, reordered)
    assert report["verdict"] == "red"
    assert report["differences"]["VEP run statistics"] == [{
        "kind": "order",
        "expected": ["VEP version (API)", "Annotation sources", "Species"] + [k for k, _ in cvs.parse_report(reference)[0].rows[3:]],
        "actual": ["Species", "Annotation sources", "VEP version (API)"] + [k for k, _ in cvs.parse_report(reference)[0].rows[3:]],
    }]
    # A key on one side alone is a set difference, whatever its count.
    extra = with_section_lines(reference, "Variant classes", [*section_lines(reference, "Variant classes"), "inversion\t0"])
    report = cvs.compare_reports(reference, extra)
    assert report["differences"]["Variant classes"] == [
        {"kind": "unexpected", "row": "inversion", "expected": None, "actual": 0}
    ]
    assert report["reasons"] == ["[Variant classes]: 1 difference"]
    # A key printed twice is reported, never merged.
    doubled = with_section_lines(reference, "Variant classes", [*section_lines(reference, "Variant classes"), "SNV\t1"])
    report = cvs.compare_reports(reference, doubled)
    assert report["differences"]["Variant classes"] == [{"kind": "duplicate", "row": "SNV", "side": "actual"}]
    # In a tabular section the first value of a repeated key is the one compared, and the
    # repeat adds no order difference.
    rows = section_lines(reference, cvs.GENERAL_STATISTICS)
    repeated = with_section_lines(reference, cvs.GENERAL_STATISTICS, [*rows[:5], "Overlapped genes\t1", *rows[5:]])
    report = cvs.compare_reports(reference, repeated)
    assert report["differences"]["General statistics"] == [
        {"kind": "duplicate", "row": "Overlapped genes", "side": "actual"}
    ]


def test_tolerance_applies_per_chart_section(reference):
    off_by_one = with_row(reference, "Consequences (all)", "intron_variant", "964")
    red = cvs.compare_reports(reference, off_by_one)
    assert red["verdict"] == "red"
    assert red["differences"]["Consequences (all)"] == [
        {"kind": "count", "row": "intron_variant", "expected": 963, "actual": 964}
    ]
    green = cvs.compare_reports(reference, off_by_one, tolerance=1)
    assert green["verdict"] == "green", green["reasons"]
    assert green["tolerance"] == 1
    # The tolerance is a count bound within each chart section, not a budget across them:
    # two charts each off by one pass at 1, and one chart off by two fails at 1.
    two_charts = with_row(off_by_one, "Coding consequences", "missense_variant", "362")
    assert cvs.compare_reports(reference, two_charts, tolerance=1)["verdict"] == "green"
    off_by_two = with_row(reference, "Consequences (all)", "intron_variant", "965")
    assert cvs.compare_reports(reference, off_by_two, tolerance=1)["verdict"] == "red"
    assert cvs.compare_reports(reference, off_by_two, tolerance=2)["verdict"] == "green"
    # Under a tolerance a key on one side alone counts as 0 on the other.
    rare = with_section_lines(reference, "Variant classes", [*section_lines(reference, "Variant classes"), "inversion\t1"])
    assert cvs.compare_reports(reference, rare)["verdict"] == "red"
    assert cvs.compare_reports(reference, rare, tolerance=1)["verdict"] == "green"
    # The tabular sections stay exact under any tolerance.
    lines = with_row(reference, cvs.GENERAL_STATISTICS, "Lines of input read", "341")
    assert cvs.compare_reports(reference, lines, tolerance=5)["verdict"] == "red"
    with pytest.raises(cvs.ComparatorError, match="negative"):
        cvs.compare_reports(reference, reference, tolerance=-1)


def test_verdict_is_red_on_any_count_difference(reference):
    changed = with_row(reference, cvs.GENERAL_STATISTICS, "Overlapped transcripts", "842")
    report = cvs.compare_reports(reference, changed)
    assert report["verdict"] == "red"
    assert report["reasons"] == ["[General statistics]: 1 difference"]
    assert report["differences"] == {
        "General statistics": [
            {"kind": "value", "row": "Overlapped transcripts", "expected": "843", "actual": "842"}
        ]
    }
    # A [Data version] row, a missing row and a missing section are each red.
    assert cvs.compare_reports(reference, with_row(reference, cvs.DATA_VERSION, "dbSNP", "157"))["verdict"] == "red"
    rows = section_lines(reference, cvs.GENERAL_STATISTICS)
    fewer = with_section_lines(reference, cvs.GENERAL_STATISTICS, rows[:-1])
    report = cvs.compare_reports(reference, fewer)
    assert report["differences"]["General statistics"] == [
        {"kind": "missing", "row": "Overlapped regulatory features", "expected": "-"}
    ]
    predictions = PREDICTIONS.read_text(encoding="utf-8")
    report = cvs.compare_reports(reference, predictions)
    assert report["verdict"] == "red"
    assert report["reasons"] == [
        "section [SIFT summary] only in the actual file",
        "section [PolyPhen summary] only in the actual file",
    ]
    assert report["sections"]["actual"][7:9] == ["SIFT summary", "PolyPhen summary"]
    assert report["sections"]["order_agrees"] is True
    # A run with a FASTA writes the distribution charts between [Variants by chromosome]
    # and [Position in protein]; the same sections in another order are red.
    titles = [s.title for s in cvs.parse_report(reference)]
    moved = reference.replace("\n[Position in protein]\n", "\n[Position in protein moved]\n")
    moved = moved.replace("\n[Variants by chromosome]\n", "\n[Position in protein]\n").replace(
        "\n[Position in protein moved]\n", "\n[Variants by chromosome]\n"
    )
    report = cvs.compare_reports(reference, moved)
    assert report["sections"]["order_agrees"] is False
    assert report["reasons"][0] == "the shared sections come in another order in the actual file"
    assert report["sections"]["expected"] == titles
    # A repeated section title is reported once; on either side its first copy is the one
    # compared and the rows compared count it once.
    twice = reference + "\n[Variant classes]\nSNV\t243\n"
    report = cvs.compare_reports(reference, twice)
    assert report["reasons"] == ["section [Variant classes] appears twice in the actual file"]
    assert report["sections"]["order_agrees"] is True
    report = cvs.compare_reports(twice, reference)
    assert report["reasons"] == ["section [Variant classes] appears twice in the expected file"]
    assert report["compared"] == cvs.compare_reports(reference, reference)["compared"]
    # Two empty reports agree on nothing.
    assert cvs.compare_reports("", "")["reasons"] == [
        "the expected file carries no [section] header",
        "the actual file carries no [section] header",
    ]


def test_report_shape_matches_compare_vep_formats(reference, tmp_path, capsys):
    empty = tmp_path / "empty_warnings.txt"
    empty.write_text("", encoding="utf-8")
    formats_report = cvf.compare_warnings(empty, empty, None, None)
    report = cvs.compare_reports(reference, reference)
    shared = {"mode", "format", "verdict", "reasons"}
    assert shared <= set(formats_report)
    assert shared <= set(report)
    assert report["mode"] == report["format"] == "stats"
    assert set(report) == shared | {"tolerance", "masked_rows", "sections", "compared", "differences"}
    assert report["verdict"] in {"green", "red"} and isinstance(report["reasons"], list)
    # The command line: the paths under `expected` and `actual`, the JSON report written,
    # `<label>: <verdict>` on standard output, exit 0 on green, 1 on red, 2 on an error.
    actual = tmp_path / "candidate_summary.txt"
    actual.write_text(candidate_shaped(reference), encoding="utf-8")
    written = tmp_path / "out" / "stats.json"
    assert cvs.main(["--expected", str(REFERENCE), "--actual", str(actual), "--report", str(written)]) == 0
    assert capsys.readouterr().out == "stats: green\n"
    on_disk = json.loads(written.read_text(encoding="utf-8"))
    assert on_disk["expected"] == str(REFERENCE) and on_disk["actual"] == str(actual)
    assert on_disk["verdict"] == "green" and on_disk["differences"] == {}
    actual.write_text(with_row(reference, cvs.GENERAL_STATISTICS, "Variants processed", "338"), encoding="utf-8")
    assert cvs.main(["--expected", str(REFERENCE), "--actual", str(actual)]) == 1
    assert capsys.readouterr().out == "stats: red\n  [General statistics]: 1 difference\n"
    assert cvs.main([
        "--expected", str(REFERENCE), "--actual", str(actual), "--mask-row", "Variants processed",
    ]) == 0
    capsys.readouterr()
    assert cvs.main(["--expected", str(REFERENCE), "--actual", str(tmp_path / "none.txt")]) == 2
    assert capsys.readouterr().err.startswith("ERROR: [compare_vep_stats] cannot read ")
    assert cvs.main(["--expected", str(REFERENCE), "--actual", str(actual), "--tolerance", "-1"]) == 2
    assert capsys.readouterr().err == "ERROR: [compare_vep_stats] --tolerance -1 is negative\n"
    # A gzipped report reads like a plain one; a corrupt one is an error, not a traceback.
    import gzip

    zipped = tmp_path / "ref.summary.log.gz"
    with gzip.open(zipped, "wt", encoding="utf-8") as fh:
        fh.write(reference)
    assert cvs.compare_files(zipped, REFERENCE)["verdict"] == "green"
    corrupt = tmp_path / "corrupt.summary.log.gz"
    corrupt.write_bytes(zipped.read_bytes()[:40] + b"not deflate data" * 8)
    assert cvs.main(["--expected", str(corrupt), "--actual", str(REFERENCE)]) == 2
    assert capsys.readouterr().err.startswith("ERROR: [compare_vep_stats] cannot read ")
