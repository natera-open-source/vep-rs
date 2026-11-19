"""Tests for plugin_gate_rows.py: the counts at the tree's own support table, the three outputs held to
each other, and each refusal on a small planted tree."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "plugin_gate_rows.py"
sys.path.insert(0, str(HERE))
import plugin_gate_rows as pgr

REAL = "supported, verified on real data"
FIXTURES = "supported on the user's licensed file, verified on synthetic fixtures"
NOT_PORTED = "not ported"
PREFIX = "ERROR: [plugin_gate_rows] "
HEADER = "| Plugin | Data file and provider | Terms | Assemblies | Status | Reason |"
SEPARATOR = "| --- | --- | --- | --- | --- | --- |"
LABELS = [
    "removed from or absent from Ensembl VEP's release 116 collection, or resting on a retired resource",
    "needing a live database, a web service or an external program",
    "data for GRCh37 only",
    "specialised data, rarely used",
    "writing files or reports, or rewriting or filtering output",
    "commercial or restricted licence",
    "covered by something vep-rs has",
]


def table_row(name: str, assemblies: str, status: str, data: str = "a file", reason: str = "a reason") -> str:
    return f"| {name} | {data} | terms | {assemblies} | {status} | {reason} |"


# A clean planted tree: one plugin of each kind of data and one not ported.
ROWS = [
    table_row("Alpha", "GRCh37, GRCh38", REAL),
    table_row("Beta", "GRCh38", FIXTURES),
    table_row("Gamma", "GRCh37", NOT_PORTED),
    table_row("TranscriptAnnotator", "GRCh38", REAL),
]
CORPORA = [
    ("GRCh38-plugin-beta", ["Beta", "api-params-parsing"], "GRCh38"),
    ("GRCh38-plugin-transcript-annotator", ["TranscriptAnnotator"], "GRCh38"),
]
GROUPS = [("Gamma", "retired")]


def plant(root: Path, rows: list[str] = ROWS, corpora: list[tuple] = CORPORA, groups: list[tuple] = GROUPS) -> Path:
    page = root / "docs" / "plugins.md"
    page.parent.mkdir(parents=True, exist_ok=True)
    table = "\n".join([HEADER, SEPARATOR, *rows])
    page.write_text(f"# Plugins\n\nProse.\n\n## Support table\n\nProse.\n\n{table}\n\nAfter the table.\n", encoding="utf-8")
    for corpus, covers, assembly in corpora:
        directory = root / "tests" / "golden" / "116" / corpus
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "manifest.json").write_text(json.dumps({"covers": covers, "assembly": assembly}), encoding="utf-8")
    groups_file = root / pgr.GROUPS_FILE
    groups_file.parent.mkdir(parents=True, exist_ok=True)
    groups_file.write_text("".join(f"{plugin}\t{group}\n" for plugin, group in [("plugin", "group"), *groups]), encoding="utf-8")
    return root


def write_all(root: Path, out: Path) -> tuple[int, Path, Path, Path]:
    out.mkdir(parents=True, exist_ok=True)
    rows, not_ported, record = out / "rows.tsv", out / "not_ported.json", out / "record_input.json"
    code = pgr.main(["--root", str(root), "--rows", str(rows), "--not-ported", str(not_ported), "--record-input", str(record)])
    return code, rows, not_ported, record


def refusals(root: Path, tmp_path: Path, capsys) -> list[str]:
    """The refusal lines of a run that must exit 1, write nothing and print only refusals."""
    code, *written = write_all(root, tmp_path / "out")
    captured = capsys.readouterr()
    lines = captured.err.splitlines()
    assert code == 1
    assert lines and all(line.startswith(PREFIX) for line in lines)
    assert captured.out == ""
    assert not any(path.exists() for path in written)
    return [line.removeprefix(PREFIX) for line in lines]


def tsv_rows(path: Path) -> tuple[str, list[dict]]:
    """(the digest of the comment line, the rows as dicts) of a `--rows` file."""
    comment, header, *lines = path.read_text(encoding="utf-8").splitlines()
    assert comment.startswith("# support_table_sha256: ")
    assert header.split("\t") == ["plugin", "assembly", "scored_on", "dataset_or_corpus", "label"]
    return comment.removeprefix("# support_table_sha256: "), [dict(zip(header.split("\t"), line.split("\t"))) for line in lines]


def test_the_trees_table_gives_57_rows_over_31_plugins() -> None:
    rows, _block, _digest, errors = pgr.derive(ROOT)
    assert errors == []
    assert len(rows) == 57
    assert len({row["plugin"] for row in rows}) == 31
    kinds = [row["scored_on"] for row in rows]
    assert (kinds.count("clinvar"), kinds.count("committed_table"), kinds.count("synthetic_fixture")) == (35, 2, 20)


def test_the_trees_rows_name_the_two_clinvar_datasets_and_22_committed_corpora() -> None:
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    clinvar = {(row["assembly"], row["dataset_or_corpus"]) for row in rows if row["scored_on"] == "clinvar"}
    assert clinvar == {("GRCh37", "clinvar_grch37"), ("GRCh38", "clinvar_grch38")}
    corpora = [row["dataset_or_corpus"] for row in rows if row["scored_on"] != "clinvar"]
    assert len(corpora) == len(set(corpora)) == 22
    for row in rows:
        if row["scored_on"] != "clinvar":
            manifest = json.loads((ROOT / "tests/golden/116" / row["dataset_or_corpus"] / "manifest.json").read_text())
            assert row["plugin"] in manifest["covers"] and manifest["assembly"] == row["assembly"]


def test_the_trees_rows_follow_the_table_then_the_assembly() -> None:
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    table, _errors = pgr.parity.read_support_table(ROOT)
    supported = [entry["name"] for entry in table if entry["status"] != NOT_PORTED]
    assert list(dict.fromkeys(row["plugin"] for row in rows)) == supported
    for plugin in supported:
        assemblies = [row["assembly"] for row in rows if row["plugin"] == plugin]
        assert assemblies == sorted(assemblies)


def test_the_trees_cadd_structural_variant_corpus_is_a_row_under_its_label() -> None:
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    cadd = [(row["assembly"], row["dataset_or_corpus"], row["label"]) for row in rows if row["plugin"] == "CADD"]
    assert cadd == [
        ("GRCh37", "GRCh37-plugin-cadd", "CADD"),
        ("GRCh38", "GRCh38-plugin-cadd", "CADD"),
        ("GRCh38", "GRCh38-plugin-cadd-sv", "CADD (structural variants)"),
    ]
    assert [row["label"] for row in rows if row["label"] != row["plugin"]] == ["CADD (structural variants)"]


def test_the_trees_transcript_annotator_is_scored_on_its_committed_table() -> None:
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    own = [(row["assembly"], row["scored_on"], row["dataset_or_corpus"]) for row in rows if row["plugin"] == "TranscriptAnnotator"]
    assert own == [
        ("GRCh37", "committed_table", "GRCh37-plugin-transcript-annotator"),
        ("GRCh38", "committed_table", "GRCh38-plugin-transcript-annotator"),
    ]
    assert {row["plugin"] for row in rows if row["scored_on"] == "committed_table"} == {"TranscriptAnnotator"}


def test_the_trees_groups_are_56_plugins_in_seven_groups_of_87() -> None:
    _rows, block, _digest, errors = pgr.derive(ROOT)
    assert errors == []
    assert block["of"] == 87
    assert [group["label"] for group in block["groups"]] == LABELS
    assert [group["count"] for group in block["groups"]] == [16, 13, 7, 7, 5, 4, 4]
    listed = [plugin for group in block["groups"] for plugin in group["plugins"]]
    assert len(listed) == len(set(listed)) == 56
    assert all(group["count"] == len(group["plugins"]) for group in block["groups"])


def test_the_trees_groups_hold_every_not_ported_plugin_once_in_the_tables_order() -> None:
    _rows, block, _digest, _errors = pgr.derive(ROOT)
    table, _errors = pgr.parity.read_support_table(ROOT)
    not_ported = [entry["name"] for entry in table if entry["status"] == NOT_PORTED]
    assert sorted(plugin for group in block["groups"] for plugin in group["plugins"]) == sorted(not_ported)
    for group in block["groups"]:
        assert group["plugins"] == [name for name in not_ported if name in group["plugins"]]
    first = block["groups"][0]["plugins"]
    assert (first[0], first[-1]) == ("AVADA", "SubsetVCF")


def test_the_trees_two_rows_with_an_escaped_pipe_keep_their_assemblies() -> None:
    page = (ROOT / "docs/plugins.md").read_text(encoding="utf-8").splitlines()
    table, _errors = pgr.parity.read_support_table(ROOT)
    escaped = [entry["name"] for entry in table if "\\|" in page[entry["lineno"] - 1]]
    assert escaped == ["DosageSensitivity", "LOEUF"]
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    for plugin in escaped:
        assert [row["assembly"] for row in rows if row["plugin"] == plugin] == ["GRCh37", "GRCh38"]


def test_an_escaped_pipe_before_the_assemblies_cell_does_not_move_the_cell(tmp_path: Path) -> None:
    rows = [table_row("Alpha", "GRCh38", REAL, data="`mode=a\\|b` of a file", reason="`cover=0\\|1`"), *ROWS[1:]]
    derived, _block, _digest, errors = pgr.derive(plant(tmp_path, rows))
    assert errors == []
    assert [row["assembly"] for row in derived if row["plugin"] == "Alpha"] == ["GRCh38"]


def test_an_assemblies_cell_is_read_without_its_parenthesised_parts() -> None:
    assert pgr.assemblies_of("GRCh37, GRCh38") == ["GRCh37", "GRCh38"]
    assert pgr.assemblies_of("GRCh37, GRCh38 (CADD-SV on GRCh38)") == ["GRCh37", "GRCh38"]
    assert pgr.assemblies_of("GRCh37 (a note); GRCh38 by a note on GRCh37") == ["GRCh37", "GRCh38"]
    assert pgr.assemblies_of("GRCh37 (GRCh38 by remapping)") == ["GRCh37"]
    assert pgr.assemblies_of("GRCh38 (the matches table)") == ["GRCh38"]
    assert pgr.assemblies_of("any; keyed on transcript id") == []
    assert pgr.assemblies_of("none") == []


def test_the_trees_parenthesised_cells_give_the_assemblies_outside_the_parentheses() -> None:
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    assert [row["assembly"] for row in rows if row["plugin"] == "BayesDel"] == ["GRCh37", "GRCh38"]
    assert [row["assembly"] for row in rows if row["plugin"] == "CADD"] == ["GRCh37", "GRCh38", "GRCh38"]


def test_an_assembly_named_only_inside_parentheses_gives_no_row(tmp_path: Path) -> None:
    rows = [ROWS[0], table_row("Beta", "GRCh38 (GRCh37 by remapping)", FIXTURES), *ROWS[2:]]
    derived, _block, _digest, errors = pgr.derive(plant(tmp_path, rows))
    assert errors == []
    assert [row["assembly"] for row in derived if row["plugin"] == "Beta"] == ["GRCh38"]


def test_the_planted_tree_gives_its_rows_and_groups(tmp_path: Path) -> None:
    rows, block, _digest, errors = pgr.derive(plant(tmp_path))
    assert errors == []
    assert [tuple(row.values()) for row in rows] == [
        ("Alpha", "GRCh37", "clinvar", "clinvar_grch37", "Alpha"),
        ("Alpha", "GRCh38", "clinvar", "clinvar_grch38", "Alpha"),
        ("Beta", "GRCh38", "synthetic_fixture", "GRCh38-plugin-beta", "Beta"),
        ("TranscriptAnnotator", "GRCh38", "committed_table", "GRCh38-plugin-transcript-annotator", "TranscriptAnnotator"),
    ]
    assert block["of"] == 4
    assert block["groups"][0] == {"label": LABELS[0], "count": 1, "plugins": ["Gamma"]}
    assert [group["count"] for group in block["groups"][1:]] == [0] * 6


def test_a_corpus_that_stacks_several_plugins_or_names_none_is_no_row(tmp_path: Path) -> None:
    corpora = [
        *CORPORA,
        ("GRCh38-plugin-stack", ["Alpha", "Beta"], "GRCh38"),
        ("GRCh38-plugin-other", ["api-params-parsing"], "GRCh38"),
        ("GRCh38-plugin-no-covers", None, "GRCh38"),
    ]
    rows, _block, _digest, errors = pgr.derive(plant(tmp_path, corpora=corpora))
    assert errors == []
    assert [row["dataset_or_corpus"] for row in rows if row["scored_on"] != "clinvar"] == [
        "GRCh38-plugin-beta",
        "GRCh38-plugin-transcript-annotator",
    ]


def test_a_second_corpus_with_a_fixed_label_is_a_row_after_the_plugins_own(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(pgr.CORPUS_LABELS, "GRCh38-plugin-beta-wide", "Beta (wide records)")
    corpora = [("GRCh38-plugin-beta-wide", ["Beta"], "GRCh38"), *CORPORA]
    rows, _block, _digest, errors = pgr.derive(plant(tmp_path, corpora=corpora))
    assert errors == []
    assert [(row["dataset_or_corpus"], row["label"]) for row in rows if row["plugin"] == "Beta"] == [
        ("GRCh38-plugin-beta", "Beta"),
        ("GRCh38-plugin-beta-wide", "Beta (wide records)"),
    ]


def test_every_label_is_the_plugins_name_or_opens_with_it() -> None:
    rows, _block, _digest, _errors = pgr.derive(ROOT)
    for row in rows:
        assert row["label"] == row["plugin"] or row["label"].startswith(f"{row['plugin']} (")


def test_the_three_outputs_carry_the_same_rows_digest_and_groups(tmp_path: Path, capsys) -> None:
    code, rows_path, not_ported_path, record_path = write_all(ROOT, tmp_path / "out")
    assert code == 0
    assert capsys.readouterr().out == (
        "plugin gate rows: 57 rows over 31 plugins (35 clinvar, 20 synthetic_fixture, 2 committed_table); "
        "56 of 87 plugins not ported, in 7 groups\n"
    )
    digest, rows = tsv_rows(rows_path)
    not_ported = json.loads(not_ported_path.read_text(encoding="utf-8"))
    record = json.loads(record_path.read_text(encoding="utf-8"))
    assert list(record) == ["support_table_sha256", "rows", "not_ported"]
    assert record["support_table_sha256"] == digest
    assert record["not_ported"] == not_ported
    assert list(not_ported) == ["of", "groups"]
    assert all(list(group) == ["label", "count", "plugins"] for group in not_ported["groups"])
    assert len(record["rows"]) == len(rows) == 57
    for written, row in zip(record["rows"], rows):
        assert list(written) == ["plugin", "assembly", "scored_on", "dataset_or_corpus"]
        assert written["plugin"] == row["label"]
        for key in ("assembly", "scored_on", "dataset_or_corpus"):
            assert written[key] == row[key]
    labelled = [written["plugin"] for written, row in zip(record["rows"], rows) if row["label"] != row["plugin"]]
    assert labelled == ["CADD (structural variants)"]


def test_the_digest_is_of_the_tables_lines_alone(tmp_path: Path) -> None:
    root = plant(tmp_path)
    _rows, _block, digest, _errors = pgr.derive(root)
    table = "".join(f"{line}\n" for line in [HEADER, SEPARATOR, *ROWS])
    assert digest == hashlib.sha256(table.encode("utf-8")).hexdigest()
    page = root / "docs" / "plugins.md"
    page.write_text(page.read_text(encoding="utf-8").replace("After the table.", "Other prose."), encoding="utf-8")
    assert pgr.derive(root)[2] == digest
    page.write_text(page.read_text(encoding="utf-8").replace("a reason", "another reason", 1), encoding="utf-8")
    assert pgr.derive(root)[2] != digest


def test_the_trees_digest_is_of_its_89_table_lines() -> None:
    lines = (ROOT / "docs/plugins.md").read_text(encoding="utf-8").splitlines()
    start = lines.index(HEADER)
    table = lines[start : start + 2 + 87]
    assert table[-1].startswith("| VARITY ") and not lines[start + 2 + 87].startswith("|")
    assert pgr.derive(ROOT)[2] == hashlib.sha256("".join(f"{line}\n" for line in table).encode("utf-8")).hexdigest()


def test_a_rerun_writes_byte_identical_files(tmp_path: Path) -> None:
    runs = []
    for name in ("first", "second"):
        out = tmp_path / name
        out.mkdir()
        paths = [out / "rows.tsv", out / "not_ported.json", out / "record_input.json"]
        done = subprocess.run(
            [sys.executable, str(SCRIPT), "--root", str(ROOT), "--rows", str(paths[0]), "--not-ported", str(paths[1]), "--record-input", str(paths[2])],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            check=False,
        )
        assert (done.returncode, done.stderr) == (0, "")
        runs.append([path.read_bytes() for path in paths] + [done.stdout.encode()])
    assert runs[0] == runs[1]
    assert all(data.endswith(b"\n") for data in runs[0])


def test_a_run_without_an_output_option_checks_and_writes_nothing(tmp_path: Path) -> None:
    command = [sys.executable, str(SCRIPT), "--root", str(plant(tmp_path / "tree"))]
    done = subprocess.run(command, capture_output=True, text=True, cwd=tmp_path, check=False)
    assert (done.returncode, done.stderr) == (0, "")
    assert done.stdout == "plugin gate rows: 4 rows over 3 plugins (2 clinvar, 1 synthetic_fixture, 1 committed_table); 1 of 4 plugins not ported, in 7 groups\n"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["tree"]


def test_a_not_ported_plugin_with_no_group_is_refused(tmp_path: Path, capsys) -> None:
    rows = [*ROWS, table_row("Zeta", "GRCh37", NOT_PORTED)]
    assert refusals(plant(tmp_path, rows), tmp_path, capsys) == [
        "docs/plugins.md:15: Zeta is not ported and has no group in scripts/parity/not_ported_groups.tsv"
    ]


def test_a_not_ported_plugin_with_two_groups_is_refused(tmp_path: Path, capsys) -> None:
    groups = [("Gamma", "retired"), ("Gamma", "covered")]
    assert refusals(plant(tmp_path, groups=groups), tmp_path, capsys) == [
        "scripts/parity/not_ported_groups.tsv: Gamma has 2 groups (lines 2, 3); a plugin has one"
    ]


def test_a_not_ported_plugin_listed_twice_in_one_group_is_refused(tmp_path: Path, capsys) -> None:
    groups = [("Gamma", "retired"), ("Gamma", "retired")]
    assert refusals(plant(tmp_path, groups=groups), tmp_path, capsys) == [
        "scripts/parity/not_ported_groups.tsv: Gamma has 2 groups (lines 2, 3); a plugin has one"
    ]


def test_a_group_row_for_a_supported_plugin_is_refused(tmp_path: Path, capsys) -> None:
    groups = [*GROUPS, ("Alpha", "covered")]
    assert refusals(plant(tmp_path, groups=groups), tmp_path, capsys) == [
        "scripts/parity/not_ported_groups.tsv:3: Alpha is not a plugin the support table lists as not ported"
    ]


def test_a_group_row_for_a_plugin_the_table_lacks_is_refused(tmp_path: Path, capsys) -> None:
    groups = [*GROUPS, ("Omega", "licence")]
    assert refusals(plant(tmp_path, groups=groups), tmp_path, capsys) == [
        "scripts/parity/not_ported_groups.tsv:3: Omega is not a plugin the support table lists as not ported"
    ]


def test_a_supported_plugin_with_no_row_is_refused(tmp_path: Path, capsys) -> None:
    rows = [table_row("Alpha", "any; keyed on gene", REAL), *ROWS[1:]]
    assert refusals(plant(tmp_path, rows), tmp_path, capsys) == [
        "docs/plugins.md:11: Alpha is supported and has no row: its Assemblies cell names no assembly"
    ]


def test_a_status_the_table_reader_refuses_is_refused(tmp_path: Path, capsys) -> None:
    rows = [table_row("Alpha", "GRCh38", "supported"), *ROWS[1:]]
    lines = refusals(plant(tmp_path, rows), tmp_path, capsys)
    assert len(lines) == 1
    assert lines[0].startswith("docs/plugins.md: support table malformed at line 11: status 'supported' is not one of ")


def test_a_status_with_no_kind_of_data_is_refused(tmp_path: Path, capsys, monkeypatch) -> None:
    status = "supported, verified on a table of its own"
    monkeypatch.setitem(pgr.parity.STATUS_TARGETS, status, "implemented")
    rows = [table_row("Alpha", "GRCh38", status), *ROWS[1:]]
    assert refusals(plant(tmp_path, rows), tmp_path, capsys) == [
        f"docs/plugins.md:11: Alpha has the status {status!r}, which names no kind of data here"
    ]


def test_a_fixture_plugin_and_assembly_with_no_committed_corpus_is_refused(tmp_path: Path, capsys) -> None:
    rows = [ROWS[0], table_row("Beta", "GRCh37, GRCh38", FIXTURES), *ROWS[2:]]
    assert refusals(plant(tmp_path, rows), tmp_path, capsys) == [
        "docs/plugins.md:12: Beta is scored on corpora, and tests/golden holds none of it on GRCh37"
    ]


def test_a_committed_table_plugin_with_no_committed_corpus_is_refused(tmp_path: Path, capsys) -> None:
    assert refusals(plant(tmp_path, corpora=CORPORA[:1]), tmp_path, capsys) == [
        "docs/plugins.md:14: TranscriptAnnotator is scored on corpora, and tests/golden holds none of it on GRCh38"
    ]


def test_two_corpora_of_one_plugin_and_assembly_under_one_label_are_refused(tmp_path: Path, capsys) -> None:
    corpora = [*CORPORA, ("GRCh38-plugin-beta-wide", ["Beta"], "GRCh38")]
    assert refusals(plant(tmp_path, corpora=corpora), tmp_path, capsys) == [
        (
            "docs/plugins.md:12: Beta has the corpora GRCh38-plugin-beta and GRCh38-plugin-beta-wide on GRCh38 under "
            "the one label 'Beta'; a second corpus needs a label of its own in CORPUS_LABELS"
        )
    ]


def test_a_real_data_plugin_on_an_assembly_without_a_clinvar_dataset_is_refused(tmp_path: Path, capsys) -> None:
    rows = [table_row("Alpha", "GRCh36, GRCh38", REAL), *ROWS[1:]]
    assert refusals(plant(tmp_path, rows), tmp_path, capsys) == [
        "docs/plugins.md:11: Alpha is scored on ClinVar, and GRCh36 has no genome-wide ClinVar dataset"
    ]


def test_a_manifest_that_cannot_be_read_is_refused(tmp_path: Path, capsys) -> None:
    root = plant(tmp_path)
    (root / "tests/golden/116/GRCh38-plugin-beta/manifest.json").write_text("{", encoding="utf-8")
    lines = refusals(root, tmp_path, capsys)
    assert lines[0].startswith("tests/golden/116/GRCh38-plugin-beta/manifest.json is not JSON")
    assert lines[1:] == ["docs/plugins.md:12: Beta is scored on corpora, and tests/golden holds none of it on GRCh38"]


def test_a_plugin_corpus_whose_manifest_names_no_assembly_is_refused(tmp_path: Path, capsys) -> None:
    corpora = [("GRCh38-plugin-beta", ["Beta"], None), CORPORA[1]]
    lines = refusals(plant(tmp_path, corpora=corpora), tmp_path, capsys)
    assert lines[0] == "tests/golden/116/GRCh38-plugin-beta/manifest.json names no assembly"


def test_a_group_row_that_is_not_two_cells_is_refused(tmp_path: Path, capsys) -> None:
    root = plant(tmp_path)
    (root / pgr.GROUPS_FILE).write_text("plugin\tgroup\nGamma\n", encoding="utf-8")
    assert refusals(root, tmp_path, capsys) == ["scripts/parity/not_ported_groups.tsv:2: not a plugin and a group"]


def test_a_group_row_naming_no_group_is_refused(tmp_path: Path, capsys) -> None:
    lines = refusals(plant(tmp_path, groups=[("Gamma", "obsolete")]), tmp_path, capsys)
    assert lines == [
        (
            "scripts/parity/not_ported_groups.tsv:2: Gamma names the group 'obsolete', not one of "
            "retired, live_service, grch37_only, specialised, writes_or_filters, licence, covered"
        )
    ]


def test_a_groups_file_that_is_missing_or_has_another_header_is_refused(tmp_path: Path, capsys) -> None:
    root = plant(tmp_path)
    (root / pgr.GROUPS_FILE).write_text("name\tgroup\nGamma\tretired\n", encoding="utf-8")
    assert refusals(root, tmp_path, capsys) == ["scripts/parity/not_ported_groups.tsv:1: the header is not ['plugin', 'group']"]
    (root / pgr.GROUPS_FILE).unlink()
    assert refusals(root, tmp_path, capsys) == ["scripts/parity/not_ported_groups.tsv is missing"]


def test_a_table_with_no_rows_is_refused(tmp_path: Path, capsys) -> None:
    assert refusals(plant(tmp_path, rows=[]), tmp_path, capsys) == ["docs/plugins.md: the support table has no rows"]


def test_every_refusal_of_a_run_is_printed(tmp_path: Path, capsys) -> None:
    rows = [table_row("Alpha", "none", REAL), table_row("Beta", "GRCh37", FIXTURES), *ROWS[2:], table_row("Zeta", "GRCh37", NOT_PORTED)]
    assert refusals(plant(tmp_path, rows), tmp_path, capsys) == [
        "docs/plugins.md:11: Alpha is supported and has no row: its Assemblies cell names no assembly",
        "docs/plugins.md:12: Beta is scored on corpora, and tests/golden holds none of it on GRCh37",
        "docs/plugins.md:15: Zeta is not ported and has no group in scripts/parity/not_ported_groups.tsv",
    ]


def test_the_trees_data_file_lists_each_not_ported_plugin_once() -> None:
    lines = (ROOT / pgr.GROUPS_FILE).read_text(encoding="utf-8").splitlines()
    assert lines[0] == "plugin\tgroup"
    plugins = [line.split("\t")[0] for line in lines[1:]]
    assert len(plugins) == len(set(plugins)) == 56
    assert plugins == sorted(plugins, key=str.casefold)
