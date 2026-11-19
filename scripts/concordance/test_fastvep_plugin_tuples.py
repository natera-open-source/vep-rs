"""Tests for fastvep_plugin_tuples.py: the plugin tuples it writes from the
supplementary columns of fastVEP's tab output, its field map and its refusals.

Every supplementary cell is fabricated in the layout its column documents, except the
lines of `BINARY_LINES`, which the binary wrote over the repository's synthetic fixtures;
no line comes from a source file. The script runs as a subprocess, the way the harness
runs it.

Run with: ``python3 -m pytest -q scripts/concordance/test_fastvep_plugin_tuples.py``
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple
from urllib.parse import quote

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "fastvep_plugin_tuples.py"
NORMALIZER = HERE / "fastvep_normalize.py"
FIELD_MAP = HERE / "fastvep_field_map.tsv"
GOLDEN = HERE.parents[1] / "tests" / "golden" / "116"
BUILT_IN = [
    "#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "Extra",
]
MAP_HEADER = (
    "fastvep_column\tpart\tvep_field\tsource\tkind\tnote\tlist_separator\treference_set"
)
CUSTOM_SOURCES = {"ClinVar", "gnomAD"}
ERROR = "ERROR: [fastvep_plugin_tuples] "
QUIET = (
    "values outside --fields: 0; parts not written word(score): 0; "
    "cells with several entries: 0; unmapped columns: none"
)


class Run(NamedTuple):
    returncode: int
    stdout: str
    stderr: str
    tuple_file: Path

    def tuples(self) -> list[tuple[str, ...]]:
        text = self.tuple_file.read_text(encoding="utf-8")
        return [tuple(line.split("\t")) for line in text.split("\n")[:-1]]


def _row(location: str, allele: str, feature: str, *cells: str, gene: str = "ENSG1") -> str:
    return "\t".join(
        ["v", location, allele, gene, feature, "Transcript", "missense_variant", "10", "10",
         "4", "A/T", "gcc/Acc", "-", "IMPACT=MODERATE", *cells]
    )


def _table(columns: list[str], rows: list[str]) -> str:
    return "\t".join([*BUILT_IN, *columns]) + "\n" + "".join(row + "\n" for row in rows)


def _script(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True, text=True, encoding="utf-8", check=False,
    )


def _run_text(tmp_path: Path, text: str, *options: str) -> Run:
    raw = tmp_path / "raw.tab"
    raw.write_text(text, encoding="utf-8")
    tuple_file = tmp_path / "tuples.tsv"
    proc = _script(str(raw), str(tuple_file), *options)
    return Run(proc.returncode, proc.stdout, proc.stderr, tuple_file)


def _run(tmp_path: Path, columns: list[str], rows: list[str], *options: str) -> Run:
    return _run_text(tmp_path, _table(columns, rows), *options)


def _one_cell(tmp_path: Path, column: str, cell: str, *options: str) -> Run:
    return _run(tmp_path, [column], [_row("21:1000", "G", "ENST1", cell)], *options)


def _at(pairs: list[tuple[str, str]], location: str = "21:1000", allele: str = "G",
        feature: str = "ENST1") -> list[tuple[str, ...]]:
    return [(location, allele, feature, field, value) for field, value in pairs]


def _map_rows() -> list[dict[str, str]]:
    lines = FIELD_MAP.read_text(encoding="utf-8").splitlines()
    names = lines[0].split("\t")
    return [dict(zip(names, line.split("\t"), strict=True)) for line in lines[1:]]


# One layout per source


GNOMAD_FIELDS = [
    "gnomAD_AF", "gnomAD_AC", "gnomAD_AN", "gnomAD_nhomalt", "gnomAD_AF_afr", "gnomAD_AF_amr",
    "gnomAD_AF_asj", "gnomAD_AF_eas", "gnomAD_AF_fin", "gnomAD_AF_mid", "gnomAD_AF_nfe",
    "gnomAD_AF_oth", "gnomAD_AF_remaining", "gnomAD_AF_sas",
]
GNOMAD_VALUES = [
    "0.0123", "37", "3008", "1", "0.02", "1.5e-5", "0", "0.001", "0.25", "", "0.0004", "0.5",
    "0.125", "3e-7",
]
LAYOUTS = {
    "alphamissense": (
        "FV_ALPHAMISSENSE",
        "G|0.9876|likely_pathogenic",
        [("am_pathogenicity", "0.9876"), ("am_class", "likely_pathogenic")],
    ),
    "revel": ("FV_REVEL", "G|0.5", [("REVEL", "0.5")]),
    "spliceai": (
        "SpliceAI",
        "G|GENE1|0.01|0.00|0.85|0.00|-12|3|0|27",
        [("SpliceAI_pred", "GENE1|0.01|0.00|0.85|0.00|-12|3|0|27")],
    ),
    "primateai": ("FV_PRIMATEAI", "G|0.7312", [("PrimateAI", "0.7312")]),
    "dbnsfp": (
        "FV_DBNSFP",
        "G|deleterious(0.012)|probably_damaging(0.998)",
        [("SIFT_score", "0.012"), ("Polyphen2_HDIV_score", "0.998")],
    ),
    "gene_constraint": (
        "FV_GNOMAD_GENE",
        "GENE1|0.9981|0.2345|3.21|-0.45",
        [("pLI_gene_value", "0.9981"), ("LOEUF", "0.2345")],
    ),
    "clinvar": (
        "FV_CLINVAR",
        (
            "G|Pathogenic&risk_factor|criteria_provided%2C_multiple_submitters%2C_no_conflicts"
            "|Condition_one&Condition_two|single_nucleotide_variant|SO:0001483"
        ),
        [
            ("ClinVar_CLNSIG", "Pathogenic|risk_factor"),
            ("ClinVar_CLNREVSTAT", "criteria_provided,_multiple_submitters,_no_conflicts"),
            ("ClinVar_CLNDN", "Condition_one|Condition_two"),
            ("ClinVar_CLNVC", "single_nucleotide_variant"),
            ("ClinVar_CLNVCSO", "SO:0001483"),
        ],
    ),
    "gnomad": (
        "FV_GNOMAD",
        "G|" + "|".join(GNOMAD_VALUES),
        [pair for pair in zip(GNOMAD_FIELDS, GNOMAD_VALUES, strict=True) if pair[1]],
    ),
}


@pytest.mark.parametrize("layout", sorted(LAYOUTS))
def test_layout_becomes_one_tuple_per_populated_part(tmp_path: Path, layout: str) -> None:
    column, cell, pairs = LAYOUTS[layout]
    run = _one_cell(tmp_path, column, cell)
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at(pairs)


def test_gnomad_layout_keeps_the_printed_number_text_and_skips_its_empty_part(
    tmp_path: Path,
) -> None:
    column, cell, _ = LAYOUTS["gnomad"]
    values = {row[3]: row[4] for row in _one_cell(tmp_path, column, cell).tuples()}
    assert values["gnomAD_AF_amr"] == "1.5e-5" and values["gnomAD_AF_asj"] == "0"
    assert "gnomAD_AF_mid" not in values and len(values) == 13


def test_gene_constraint_parts_without_a_field_are_dropped_and_counted(tmp_path: Path) -> None:
    column, cell, _ = LAYOUTS["gene_constraint"]
    run = _one_cell(tmp_path, column, cell)
    assert run.stdout == ""
    assert run.stderr == (
        "fastvep_plugin_tuples: wrote 2 tuples from 1 rows (0 too short); "
        "by field: pLI_gene_value=1, LOEUF=1; mapped columns: FV_GNOMAD_GENE; "
        "dropped parts with no Ensembl VEP field: 2 (FV_GNOMAD_GENE:3=1, FV_GNOMAD_GENE:4=1); "
        + QUIET
        + "\n"
    )


# Lines as the binary wrote them


# Lines of fastVEP 0.5.0's tab output over the repository's synthetic fixtures: the
# header lines whole, then rows with an empty cell, a lowercase leading part, a deleted
# allele and no feature. The tuples are those a plain split of each cell gives.
BINARY_LINES = {
    'spliceai': (
        [
            '## fastVEP output',
            '## COLUMN=<ID=SpliceAI,Description="SpliceAI annotations. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">',
            '#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSpliceAI',
        ],
        [
            'snv_intergenic\t21:15073000\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-',
            'snv_transcript_without_symbol\t21:16719000\tG\tENSG00000237735\tENST00000753406\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-1\t-\tG|RNU1-98P|0.92|0.57|0.95|0.05|-36|39|-15|38',
            'del_two_genes\t21:17620001-17620002\t-\tENSG00000154639\tENST00001128578\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tC|BTG3-AS1|0.93|0.55|0.23|0.87|-49|-4|-20|24',
            'snv_lowercase\t21:18341500\tT\tENSG00000154646\tENST00000284885\tTranscript\tmissense_variant\t1640\t1477\t493\tD/N\tGat/Aat\t-\tMODERATE\t-\t-1\t-\tt|TMPRSS15|0.28|0.39|0.51|0.41|43|-5|27|30',
        ],
        [
            ('21:16719000', 'G', 'ENST00000753406', 'SpliceAI_pred', 'RNU1-98P|0.92|0.57|0.95|0.05|-36|39|-15|38'),
            ('21:17620001-17620002', '-', 'ENST00001128578', 'SpliceAI_pred', 'BTG3-AS1|0.93|0.55|0.23|0.87|-49|-4|-20|24'),
            ('21:18341500', 'T', 'ENST00000284885', 'SpliceAI_pred', 'TMPRSS15|0.28|0.39|0.51|0.41|43|-5|27|30'),
        ],
    ),
    'primateai': (
        [
            '## fastVEP output',
            '## COLUMN=<ID=FV_PRIMATEAI,Description="fastVEP PrimateAI annotations. Format: ALLELE|SCORE">',
            '#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tFV_PRIMATEAI',
        ],
        [
            'intergenic_snv\t21:19268555\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\tT|0.0017',
            'synonymous_fwd_lowercase\t21:30396095\tT\tENSG00000198390\tENST00000355459\tTranscript\tsynonymous_variant\t65\t9\t3\tY\ttaC/taT\t-\tLOW\t-\t1\t-\tt|0.2914',
            'coding_mnv\t21:30396104-30396105\tAG\tENSG00000198390\tENST00000355459\tTranscript\tstop_gained\t74-75\t18-19\t6-7\tCS/*A\ttgCTct/tgAGct\t-\tHIGH\t-\t1\t-\t-',
        ],
        [
            ('21:19268555', 'T', '-', 'PrimateAI', '0.0017'),
            ('21:30396095', 'T', 'ENST00000355459', 'PrimateAI', '0.2914'),
        ],
    ),
    'dbnsfp': (
        [
            '## fastVEP output',
            '## COLUMN=<ID=FV_DBNSFP,Description="fastVEP dbNSFP annotations. Format: ALLELE|SIFT|POLYPHEN">',
            '#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tFV_DBNSFP',
        ],
        [
            'intergenic_snv\t21:19268555\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-',
            'intron_snv\t21:20204701\tC\tENSG00000233236\tENST00000746617\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-1\t-\tC|tolerated(0.959)|probably_damaging(0.137)',
            'missense_rev_two_changes\t21:29067663\tC\tENSG00000156261\tENST00001085360\tTranscript\tmissense_variant,NMD_transcript_variant\t346\t274\t92\tM/V\tAtg/Gtg\t-\tMODERATE\t-\t-1\t-\tC|tolerated(0.022)|benign(0.807)',
            'synonymous_fwd_lowercase\t21:30396095\tT\tENSG00000198390\tENST00000355459\tTranscript\tsynonymous_variant\t65\t9\t3\tY\ttaC/taT\t-\tLOW\t-\t1\t-\tt|deleterious(0.819)|benign(0.087)',
        ],
        [
            ('21:20204701', 'C', 'ENST00000746617', 'SIFT_score', '0.959'),
            ('21:20204701', 'C', 'ENST00000746617', 'Polyphen2_HDIV_score', '0.137'),
            ('21:29067663', 'C', 'ENST00001085360', 'SIFT_score', '0.022'),
            ('21:29067663', 'C', 'ENST00001085360', 'Polyphen2_HDIV_score', '0.807'),
            ('21:30396095', 'T', 'ENST00000355459', 'SIFT_score', '0.819'),
            ('21:30396095', 'T', 'ENST00000355459', 'Polyphen2_HDIV_score', '0.087'),
        ],
    ),
}


@pytest.mark.parametrize("source", sorted(BINARY_LINES))
def test_lines_the_binary_wrote_give_the_tuples_of_a_plain_split(
    tmp_path: Path, source: str
) -> None:
    header, rows, expected = BINARY_LINES[source]
    run = _run_text(tmp_path, "".join(line + "\n" for line in [*header, *rows]))
    assert run.returncode == 0, run.stderr
    assert run.tuples() == expected
    assert f" from {len(rows)} rows (0 too short); " in run.stderr


# Escapes and lists


def test_every_escape_decodes_and_only_a_bare_ampersand_splits_a_list(tmp_path: Path) -> None:
    cell = (
        "G|Likely%20benign%3B%3D%25%22|criteria_provided%2C_single_submitter"
        "|One%7Ctwo&three%26four|a&b|SO:0001483"
    )
    run = _one_cell(tmp_path, "FV_CLINVAR", cell)
    assert run.tuples() == _at(
        [
            ("ClinVar_CLNSIG", 'Likely benign;=%"'),
            ("ClinVar_CLNREVSTAT", "criteria_provided,_single_submitter"),
            ("ClinVar_CLNDN", "One|two|three&four"),
            ("ClinVar_CLNVC", "a&b"),
            ("ClinVar_CLNVCSO", "SO:0001483"),
        ]
    )


def test_a_tab_or_line_breaking_character_stays_encoded_so_a_tuple_stays_one_line(
    tmp_path: Path,
) -> None:
    breaks = "\t\n\r\x0b\x0c\x1c\x1d\x1e\x85" + chr(0x2028) + chr(0x2029)
    escapes = "".join(quote(char) + "x" for char in breaks)
    assert escapes.startswith("%09x%0Ax%0Dx%0Bx") and escapes.count("%") == 16
    run = _one_cell(tmp_path, "FV_CLINVAR", f"G|a{escapes}|raw\x0bbreak|||")
    text = run.tuple_file.read_text(encoding="utf-8")
    assert text == (
        f"21:1000\tG\tENST1\tClinVar_CLNSIG\ta{escapes}\n"
        "21:1000\tG\tENST1\tClinVar_CLNREVSTAT\traw%0Bbreak\n"
    )
    assert len(text.splitlines()) == 2


def test_an_escaped_space_decodes_to_the_space_the_default_format_writes(tmp_path: Path) -> None:
    # Ensembl VEP's default format keeps the spaces of a `--custom` value and a literal
    # `%20` as they are; its VCF format turns a run of spaces into one `_`.
    run = _one_cell(tmp_path, "FV_CLINVAR", "G||||two%20words|a%20%20b%2520c")
    assert run.tuples() == _at([("ClinVar_CLNVC", "two words"), ("ClinVar_CLNVCSO", "a  b%20c")])


def test_an_escaped_percent_sign_decodes_once(tmp_path: Path) -> None:
    run = _one_cell(tmp_path, "FV_CLINVAR", "G||||50%2520off|")
    assert run.tuples() == _at([("ClinVar_CLNVC", "50%20off")])


def test_values_are_utf8_whatever_the_locale_and_a_byte_that_is_not_passes_through(
    tmp_path: Path,
) -> None:
    row = _row("21:1000", "G", "ENST1", "G|caf%C3%A9|BYTE|esc%FFaped||") + "\n"
    raw = tmp_path / "raw.tab"
    raw.write_bytes(
        _table(["FV_CLINVAR"], []).encode("utf-8")
        + row.encode("utf-8").replace(b"BYTE", b"byte\xff")
    )
    tuple_file = tmp_path / "tuples.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(raw), str(tuple_file)],
        env={**os.environ, "LC_ALL": "C", "PYTHONUTF8": "0", "PYTHONCOERCECLOCALE": "0"},
        capture_output=True, check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert tuple_file.read_bytes() == (
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tcaf\xc3\xa9\n"
        b"21:1000\tG\tENST1\tClinVar_CLNREVSTAT\tbyte\xff\n"
        b"21:1000\tG\tENST1\tClinVar_CLNDN\tesc\xffaped\n"
    )


def test_a_raw_byte_and_an_escape_read_as_one_utf8_text(tmp_path: Path) -> None:
    rows = (
        _row("21:1000", "G", "ENST1", "G|oneRAW%85two||||x%C3%A9|,G|||||xRAWER%A9|") + "\n"
    ).encode("utf-8").replace(b"RAWER", b"\xc3").replace(b"RAW", b"\xc2")
    raw = tmp_path / "raw.tab"
    raw.write_bytes(_table(["FV_CLINVAR"], []).encode("utf-8") + rows)
    tuple_file = tmp_path / "tuples.tsv"
    proc = _script(str(raw), str(tuple_file))
    assert proc.returncode == 0, proc.stderr
    assert tuple_file.read_bytes() == (
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tone%C2%85two\n"
        b"21:1000\tG\tENST1\tClinVar_CLNVCSO\tx\xc3\xa9\n"
    )


def test_a_list_takes_the_separator_its_map_row_names(tmp_path: Path) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_text(
        MAP_HEADER + "\nFV_X\t1\tX_LIST\tX\ttext\tfabricated\t,\t-\n"
        "FV_X\t2\tX_PLAIN\tX\ttext\tfabricated\t-\t-\n",
        encoding="utf-8",
    )
    run = _one_cell(tmp_path, "FV_X", "G|a&b&c|d&e", "--field-map", str(custom))
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at([("X_LIST", "a,b,c"), ("X_PLAIN", "d&e")])


def test_a_joined_range_and_a_score_half_are_percent_decoded(tmp_path: Path) -> None:
    rows = [
        _row(
            "21:1000", "G", "ENST1",
            "G|GENE%2DA|0.01|0.00|0.85|0.00|-12|3|0|27",
            "G|probably%5Fdamaging(0%2E998)|",
        )
    ]
    run = _run(tmp_path, ["SpliceAI", "FV_DBNSFP"], rows)
    assert run.tuples() == _at(
        [
            ("SpliceAI_pred", "GENE-A|0.01|0.00|0.85|0.00|-12|3|0|27"),
            ("SIFT_score", "0.998"),
        ]
    )


# Empty columns and parts


def test_an_empty_column_writes_an_empty_tuple_file(tmp_path: Path) -> None:
    rows = [_row("21:1000", "G", "ENST1", "-"), _row("21:1000", "G", "ENST2", "")]
    run = _run(tmp_path, ["FV_REVEL"], rows)
    assert run.returncode == 0, run.stderr
    assert run.tuple_file.read_bytes() == b""
    assert run.stderr == (
        "fastvep_plugin_tuples: wrote 0 tuples from 2 rows (0 too short); "
        "by field: REVEL=0; mapped columns: FV_REVEL; "
        "dropped parts with no Ensembl VEP field: 0; " + QUIET + "\n"
    )


def test_a_header_without_rows_writes_an_empty_tuple_file(tmp_path: Path) -> None:
    run = _run(tmp_path, ["FV_REVEL"], [])
    assert run.returncode == 0, run.stderr
    assert run.tuple_file.read_bytes() == b""
    assert run.stderr.startswith(
        "fastvep_plugin_tuples: wrote 0 tuples from 0 rows (0 too short); by field: REVEL=0; "
    )


@pytest.mark.parametrize(
    ("column", "cell", "pairs"),
    [
        ("FV_GNOMAD_GENE", "GENE1|0.9981|||", [("pLI_gene_value", "0.9981")]),
        ("FV_GNOMAD_GENE", "GENE1||0.2345||", [("LOEUF", "0.2345")]),
        ("FV_ALPHAMISSENSE", "G||likely_benign", [("am_class", "likely_benign")]),
        ("FV_ALPHAMISSENSE", "G|0.5", [("am_pathogenicity", "0.5")]),
        ("FV_ALPHAMISSENSE", "G", []),
        ("SpliceAI", "G|||||||||", []),
    ],
    ids=["loeuf-empty", "pli-empty", "first-empty", "short-entry", "allele-alone", "all-empty"],
)
def test_an_empty_or_missing_part_writes_no_tuple(
    tmp_path: Path, column: str, cell: str, pairs: list[tuple[str, str]]
) -> None:
    run = _one_cell(tmp_path, column, cell)
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at(pairs)
    assert "dropped parts with no Ensembl VEP field: 0;" in run.stderr


# Rows, features and genes


def test_a_gene_level_column_is_read_from_each_row_and_never_carried_to_another_gene(
    tmp_path: Path,
) -> None:
    gene_one = "GENE1|0.9981|0.2345|3.21|-0.45"
    rows = [
        _row("21:1000", "G", "ENST1", "G|0.5", gene_one, gene="ENSG1"),
        _row("21:1000", "G", "ENST2", "G|0.5", "-", gene="ENSG2"),
        _row("21:1000", "G", "ENST3", "G|0.5", "GENE3|0.0001|1.9|0.1|0.2", gene="ENSG3"),
        _row("21:1000", "G", "ENST4", "G|0.5", gene_one, gene="ENSG1"),
    ]
    run = _run(tmp_path, ["FV_REVEL", "FV_GNOMAD_GENE"], rows)
    assert run.returncode == 0, run.stderr
    own = [("REVEL", "0.5"), ("pLI_gene_value", "0.9981"), ("LOEUF", "0.2345")]
    assert run.tuples() == [
        *_at(own, feature="ENST1"),
        *_at([("REVEL", "0.5")], feature="ENST2"),
        *_at([("REVEL", "0.5"), ("pLI_gene_value", "0.0001"), ("LOEUF", "1.9")], feature="ENST3"),
        *_at(own, feature="ENST4"),
    ]
    assert "(FV_GNOMAD_GENE:3=3, FV_GNOMAD_GENE:4=3)" in run.stderr


def test_an_allele_value_is_written_on_every_feature_row_and_on_an_intergenic_row(
    tmp_path: Path,
) -> None:
    rows = [
        _row("21:2000", "T", "ENST1", "T|0.25"),
        _row("21:2000", "T", "ENST2", "T|0.25"),
        _row("21:3000", "C", "-", "C|0.75", gene="-"),
    ]
    run = _run(tmp_path, ["FV_REVEL"], rows)
    assert run.tuples() == [
        ("21:2000", "T", "ENST1", "REVEL", "0.25"),
        ("21:2000", "T", "ENST2", "REVEL", "0.25"),
        ("21:3000", "C", "-", "REVEL", "0.75"),
    ]
    assert run.stderr.startswith(
        "fastvep_plugin_tuples: wrote 3 tuples from 3 rows (0 too short); by field: REVEL=3; "
    )


def test_a_tuple_carries_the_allele_of_its_row_not_the_leading_part_of_the_cell(
    tmp_path: Path,
) -> None:
    rows = [
        _row("21:4000", "A", "ENST1", "A|0.1"),
        _row("21:4000", "C", "ENST1", "C|0.9"),
        _row("21:5001-5002", "-", "ENST1", "T|0.3"),
    ]
    run = _run(tmp_path, ["FV_REVEL"], rows)
    assert run.tuples() == [
        ("21:4000", "A", "ENST1", "REVEL", "0.1"),
        ("21:4000", "C", "ENST1", "REVEL", "0.9"),
        ("21:5001-5002", "-", "ENST1", "REVEL", "0.3"),
    ]


def test_a_row_cut_short_of_a_supplementary_column_writes_no_tuple(tmp_path: Path) -> None:
    text = _table(["FV_REVEL"], [_row("21:1000", "G", "ENST1", "G|0.5")])
    text += _row("21:2000", "T", "ENST2") + "\n" + "v\t21:3000\tC\n" + "\n"
    run = _run_text(tmp_path, text)
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at([("REVEL", "0.5")])
    assert run.stderr.startswith("fastvep_plugin_tuples: wrote 1 tuples from 4 rows (3 too short); ")


def test_a_carriage_return_inside_a_cell_stays_in_its_value_and_never_ends_the_row(
    tmp_path: Path,
) -> None:
    rows = [
        _row("21:1000", "G", "ENST1", "G|one\rtwo|||a\r|", "G|0.5"),
        _row("21:2000", "T", "ENST2", "-", "T|0.25"),
    ]
    run = _run(tmp_path, ["FV_CLINVAR", "FV_REVEL"], rows)
    assert run.returncode == 0, run.stderr
    assert run.tuple_file.read_bytes() == (
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tone%0Dtwo\n"
        b"21:1000\tG\tENST1\tClinVar_CLNVC\ta%0D\n"
        b"21:1000\tG\tENST1\tREVEL\t0.5\n"
        b"21:2000\tT\tENST2\tREVEL\t0.25\n"
    )
    assert " from 2 rows (0 too short); " in run.stderr


def test_a_carriage_return_before_the_line_feed_ends_the_row_with_it(tmp_path: Path) -> None:
    rows = [
        _row("21:1000", "G", "ENST1", "G|0.5", "G|one\rtwo||||"),
        _row("21:2000", "T", "ENST2", "T|0.25", "-"),
    ]
    text = _table(["FV_REVEL", "FV_CLINVAR"], rows).replace("\n", "\r\n")
    run = _run_text(tmp_path, text)
    assert run.returncode == 0, run.stderr
    assert run.tuple_file.read_bytes() == (
        b"21:1000\tG\tENST1\tREVEL\t0.5\n"
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tone%0Dtwo\n"
        b"21:2000\tT\tENST2\tREVEL\t0.25\n"
    )
    assert " from 2 rows (0 too short); " in run.stderr
    assert "; mapped columns: FV_REVEL, FV_CLINVAR; " in run.stderr


def test_a_carriage_return_in_a_key_cell_is_written_encoded_so_the_tuple_stays_one_line(
    tmp_path: Path,
) -> None:
    run = _run(tmp_path, ["FV_REVEL"], [_row("21:\r1000", "-\r", "ENST\r1", "G|0.5")])
    assert run.returncode == 0, run.stderr
    assert run.tuple_file.read_bytes() == b"21:%0D1000\t-%0D\tENST%0D1\tREVEL\t0.5\n"
    assert " from 1 rows (0 too short); " in run.stderr


def test_several_entries_in_one_cell_each_write_their_parts_and_are_counted(
    tmp_path: Path,
) -> None:
    run = _one_cell(tmp_path, "FV_GNOMAD_GENE", "GENE1|0.9981|||,GENE1||0.2345||")
    assert run.tuples() == _at([("pLI_gene_value", "0.9981"), ("LOEUF", "0.2345")])
    assert "cells with several entries: 1;" in run.stderr


def test_equal_values_of_several_entries_are_written_once(tmp_path: Path) -> None:
    run = _one_cell(tmp_path, "FV_GNOMAD_GENE", "GENE1|0.5|0.25||,GENE1|0.5|0.75||")
    assert run.tuples() == _at(
        [("pLI_gene_value", "0.5"), ("LOEUF", "0.25"), ("LOEUF", "0.75")]
    )
    assert run.stderr.startswith(
        "fastvep_plugin_tuples: wrote 3 tuples from 1 rows (0 too short); "
        "by field: pLI_gene_value=1, LOEUF=2; "
    )


@pytest.mark.parametrize(
    ("cell", "pairs", "unshaped"),
    [
        ("G|deleterious|benign(0.02)", [("Polyphen2_HDIV_score", "0.02")], 1),
        ("G|(0.012)|tolerated()", [("SIFT_score", "0.012")], 0),
        ("G|a(b)(0.5)|", [("SIFT_score", "0.5")], 0),
    ],
    ids=["no-parentheses", "empty-half", "last-parenthesis"],
)
def test_a_word_score_part_is_split_only_when_it_has_that_shape(
    tmp_path: Path, cell: str, pairs: list[tuple[str, str]], unshaped: int
) -> None:
    run = _one_cell(tmp_path, "FV_DBNSFP", cell)
    assert run.tuples() == _at(pairs)
    assert f"parts not written word(score): {unshaped};" in run.stderr


def test_both_halves_of_a_word_score_part_split_at_its_last_parenthesis_and_decode(
    tmp_path: Path,
) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_text(
        MAP_HEADER + "\nFV_X\t1.word\tX_WORD\tX\ttext\tfabricated\t-\t-\n"
        "FV_X\t1.score\tX_SCORE\tX\tnumber\tfabricated\t-\t-\n",
        encoding="utf-8",
    )
    run = _one_cell(tmp_path, "FV_X", "G|a%5Fb(c)(0%2E5)", "--field-map", str(custom))
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at([("X_WORD", "a_b(c)"), ("X_SCORE", "0.5")])


def test_two_half_rows_count_one_part_that_is_not_word_score_once(tmp_path: Path) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_text(
        MAP_HEADER + "\nFV_X\t1.word\tX_WORD\tX\ttext\tfabricated\t-\t-\n"
        "FV_X\t1.score\tX_SCORE\tX\tnumber\tfabricated\t-\t-\n",
        encoding="utf-8",
    )
    run = _one_cell(tmp_path, "FV_X", "G|unsplit", "--field-map", str(custom))
    assert run.returncode == 0, run.stderr
    assert run.tuples() == []
    assert "; parts not written word(score): 1; " in run.stderr


def test_a_part_read_whole_and_by_half_writes_both_and_drops_nothing(tmp_path: Path) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_text(
        MAP_HEADER + "\nFV_X\t1.word\tX_WORD\tX\ttext\tfabricated\t-\t-\n"
        "FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t-\n",
        encoding="utf-8",
    )
    run = _one_cell(tmp_path, "FV_X", "G|a(b)", "--field-map", str(custom))
    assert run.tuples() == _at([("X_WORD", "a"), ("X_ONE", "a(b)")])
    assert "; dropped parts with no Ensembl VEP field: 0; " in run.stderr


def test_a_row_without_a_feature_cell_is_too_short_wherever_the_mapped_column_is(
    tmp_path: Path,
) -> None:
    text = "#FV_REVEL\tLocation\tAllele\tGene\tFeature\nG|0.5\t21:1000\tG\tENSG1\tENST1\nG|0.5\t21:2000\n"
    run = _run_text(tmp_path, text)
    assert run.tuples() == _at([("REVEL", "0.5")])
    assert run.stderr.startswith("fastvep_plugin_tuples: wrote 1 tuples from 2 rows (1 too short); ")


def test_a_word_score_half_without_a_row_is_dropped_and_counted(tmp_path: Path) -> None:
    column, cell, _ = LAYOUTS["dbnsfp"]
    run = _one_cell(tmp_path, column, cell)
    assert (
        "; dropped parts with no Ensembl VEP field: 2 (FV_DBNSFP:1.word=1, FV_DBNSFP:2.word=1); "
        in run.stderr
    )


def test_dropped_parts_are_listed_in_part_order(tmp_path: Path) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_text(
        MAP_HEADER + "\nFV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t-\n", encoding="utf-8"
    )
    cell = "G|" + "|".join("abcdefghijk")
    run = _one_cell(tmp_path, "FV_X", cell, "--field-map", str(custom))
    assert run.tuples() == _at([("X_ONE", "a")])
    assert "(FV_X:2=1, FV_X:3=1, " in run.stderr and "FV_X:9=1, FV_X:10=1, FV_X:11=1); " in run.stderr


def test_a_kept_field_that_writes_no_tuple_reads_zero_in_the_summary(tmp_path: Path) -> None:
    rows = [_row("21:1000", "G", "ENST1", "-", "-", "GENE1|0.9981|0.2345||")]
    run = _run(tmp_path, ["FV_ALPHAMISSENSE", "FV_REVEL", "FV_GNOMAD_GENE"], rows)
    assert run.returncode == 0, run.stderr
    assert (
        "; by field: am_pathogenicity=0, am_class=0, REVEL=0, pLI_gene_value=1, LOEUF=1; "
        in run.stderr
    )


def test_a_cell_repeated_on_the_next_row_is_counted_again_in_every_counter(
    tmp_path: Path,
) -> None:
    row = _row("21:1000", "G", "ENST1", "GENE1|0.5|0.25|1|2,GENE1|0.5|0.75||", "G|deleterious|")
    run = _run(
        tmp_path, ["FV_GNOMAD_GENE", "FV_DBNSFP"], [row, row.replace("ENST1", "ENST2")],
        "--fields", "pLI_gene_value,SIFT_score",
    )
    assert run.tuples() == [
        ("21:1000", "G", "ENST1", "pLI_gene_value", "0.5"),
        ("21:1000", "G", "ENST2", "pLI_gene_value", "0.5"),
    ]
    assert run.stderr == (
        "fastvep_plugin_tuples: wrote 2 tuples from 2 rows (0 too short); "
        "by field: pLI_gene_value=2, SIFT_score=0; mapped columns: FV_GNOMAD_GENE, FV_DBNSFP; "
        "dropped parts with no Ensembl VEP field: 4 (FV_GNOMAD_GENE:3=2, FV_GNOMAD_GENE:4=2); "
        "values outside --fields: 4; parts not written word(score): 2; "
        "cells with several entries: 2; unmapped columns: none\n"
    )


# Column header


def test_the_header_is_the_last_single_hash_line_before_the_first_row(tmp_path: Path) -> None:
    text = (
        "#fabricated comment\n"
        + _table(["FV_PRIMATEAI", "FV_REVEL"], [])
        + "## FV_REVEL\ta meta line, never the header\n"
        + _row("21:1000", "G", "ENST1", "-", "G|0.5")
        + "\n#"
        + _row("21:1500", "A", "ENST9", "A|0.1", "A|0.9")
        + "\n"
        + _row("21:2000", "T", "ENST2", "T|0.7312", "-")
        + "\n"
    )
    run = _run_text(tmp_path, text)
    assert run.returncode == 0, run.stderr
    assert run.tuples() == [
        ("21:1000", "G", "ENST1", "REVEL", "0.5"),
        ("21:2000", "T", "ENST2", "PrimateAI", "0.7312"),
    ]


def test_a_source_column_without_a_map_row_is_named_in_the_summary(tmp_path: Path) -> None:
    rows = [_row("21:1000", "G", "ENST1", "G|1.234", "G|0.5")]
    run = _run(tmp_path, ["FV_PHYLOP", "FV_REVEL"], rows)
    assert run.tuples() == _at([("REVEL", "0.5")])
    assert run.stderr.endswith("unmapped columns: FV_PHYLOP\n")


@pytest.mark.parametrize(
    ("columns", "options", "named"),
    [
        (["FV_REVEL", "FV_REVEL"], [], "FV_REVEL (columns 15 and 16)"),
        (
            ["FV_PHYLOP", "FV_REVEL", "FV_PHYLOP", "FV_GNOMAD_GENE", "FV_REVEL",
             "FV_GNOMAD_GENE", "FV_REVEL"],
            [],
            "FV_REVEL (columns 16, 19 and 21), FV_GNOMAD_GENE (columns 18 and 20)",
        ),
        (
            ["FV_GNOMAD_GENE", "FV_REVEL", "FV_REVEL"],
            ["--fields", "pLI_gene_value"],
            "FV_REVEL (columns 16 and 17)",
        ),
    ],
    ids=["one-column-twice", "two-columns-and-an-unmapped-one", "a-column-fields-leaves-out"],
)
def test_a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions(
    tmp_path: Path, columns: list[str], options: list[str], named: str
) -> None:
    (tmp_path / "tuples.tsv").write_bytes(b"an earlier run\n")
    row = _row("21:1000", "G", "ENST1", *[f"G|0.{place}" for place in range(len(columns))])
    run = _run(tmp_path, columns, [row], *options)
    assert run.returncode == 1
    assert run.stderr == f"{ERROR}the header names a mapped column more than once: {named}\n"
    assert run.tuple_file.read_bytes() == b"an earlier run\n"


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            _table([], [_row("21:1000", "G", "ENST1")]),
            "the field map names none of the header's columns: Uploaded_variation, Location,",
        ),
        (_table(["FV_PHYLOP"], []), "the field map names none of the header's columns: "),
        (_row("21:1000", "G", "ENST1", "G|0.5") + "\n", "the input has no column header line"),
        ("", "the input has no column header line"),
    ],
    ids=["no-mapped-column", "no-mapped-column-no-rows", "no-header", "empty-input"],
)
def test_an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing(
    tmp_path: Path, text: str, message: str
) -> None:
    (tmp_path / "tuples.tsv").write_bytes(b"an earlier run\n")
    run = _run_text(tmp_path, text)
    assert run.returncode == 1
    assert run.stderr.startswith(ERROR + message) and run.stderr.count("\n") == 1
    assert run.tuple_file.read_bytes() == b"an earlier run\n"


# Options and files


def test_fields_keeps_the_named_fields_and_counts_the_rest(tmp_path: Path) -> None:
    column, cell, _ = LAYOUTS["gene_constraint"]
    run = _one_cell(tmp_path, column, cell, "--fields", "pLI_gene_value")
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at([("pLI_gene_value", "0.9981")])
    assert "; by field: pLI_gene_value=1; " in run.stderr
    assert "values outside --fields: 1;" in run.stderr
    assert "dropped parts with no Ensembl VEP field: 2 " in run.stderr


def test_fields_refuses_a_name_the_map_lacks(tmp_path: Path) -> None:
    run = _one_cell(tmp_path, "FV_REVEL", "G|0.5", "--fields", "REVEL,pLI_gene")
    assert run.returncode == 1
    assert run.stderr == f"{ERROR}--fields names no row of {FIELD_MAP}: pLI_gene\n"
    assert not run.tuple_file.exists()


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ("", "--fields names no field"),
        ("REVEL,,am_class", "--fields names no field beside a comma: REVEL,,am_class"),
        ("REVEL,", "--fields names no field beside a comma: REVEL,"),
    ],
    ids=["empty-value", "between-two-commas", "after-the-last-comma"],
)
def test_fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma(
    tmp_path: Path, fields: str, message: str
) -> None:
    run = _one_cell(tmp_path, "FV_REVEL", "G|0.5", "--fields", fields)
    assert run.returncode == 1
    assert run.stderr == f"{ERROR}{message}\n"
    assert not run.tuple_file.exists()


@pytest.mark.parametrize(
    ("fields", "absent"),
    [
        ("REVEL", "REVEL (FV_REVEL)"),
        ("pLI_gene_value,am_class,REVEL", "REVEL (FV_REVEL), am_class (FV_ALPHAMISSENSE)"),
    ],
    ids=["every-named-field", "one-of-several"],
)
def test_fields_refuses_a_named_field_whose_column_the_header_lacks(
    tmp_path: Path, fields: str, absent: str
) -> None:
    column, cell, _ = LAYOUTS["gene_constraint"]
    run = _one_cell(tmp_path, column, cell, "--fields", fields)
    assert run.returncode == 1
    assert run.stderr == (
        f"{ERROR}the header lacks the column of a field --fields names: {absent}\n"
    )
    assert not run.tuple_file.exists()


@pytest.mark.parametrize(
    "options",
    [["--fields"], ["--field-map"], ["--fields", "REVEL", "--verbose"], ["a-third-path"]],
    ids=["fields-without-a-value", "map-without-a-value", "unknown-option", "third-path"],
)
def test_a_malformed_command_line_is_a_usage_error_and_writes_nothing(
    tmp_path: Path, options: list[str]
) -> None:
    run = _one_cell(tmp_path, "FV_REVEL", "G|0.5", *options)
    assert run.returncode == 2
    assert run.stderr.startswith("usage: ")
    assert not run.tuple_file.exists()


def test_an_input_that_cannot_be_read_is_refused(tmp_path: Path) -> None:
    absent = tmp_path / "absent.tab"
    proc = _script(str(absent), str(tmp_path / "tuples.tsv"))
    assert proc.returncode == 1
    assert proc.stderr == f"{ERROR}cannot read the input {absent}: No such file or directory\n"
    assert not (tmp_path / "tuples.tsv").exists()


def test_a_tuple_file_that_cannot_be_opened_is_refused(tmp_path: Path) -> None:
    raw = tmp_path / "raw.tab"
    raw.write_text(_table(["FV_REVEL"], [_row("21:1000", "G", "ENST1", "G|0.5")]), "utf-8")
    absent = tmp_path / "absent" / "tuples.tsv"
    proc = _script(str(raw), str(absent))
    assert proc.returncode == 1
    assert proc.stderr == (
        f"{ERROR}cannot write the tuple file {absent}: No such file or directory\n"
    )


def test_an_earlier_tuple_file_is_replaced_by_a_run_that_is_not_refused(tmp_path: Path) -> None:
    (tmp_path / "tuples.tsv").write_bytes(b"an earlier run\n")
    run = _one_cell(tmp_path, "FV_REVEL", "G|0.5")
    assert run.returncode == 0, run.stderr
    assert run.tuples() == _at([("REVEL", "0.5")])


@pytest.mark.skipif(not Path("/dev/full").exists(), reason="needs a device that refuses writes")
def test_a_write_that_fails_ends_in_an_error_that_is_not_a_refusal(tmp_path: Path) -> None:
    raw = tmp_path / "raw.tab"
    raw.write_text(_table(["FV_REVEL"], [_row("21:1000", "G", "ENST1", "G|0.5")]), "utf-8")
    proc = _script(str(raw), "/dev/full")
    assert proc.returncode != 0
    assert "OSError" in proc.stderr and ERROR not in proc.stderr


def test_a_tuple_file_that_is_the_input_is_refused_and_the_input_is_kept(tmp_path: Path) -> None:
    text = _table(["FV_REVEL"], [_row("21:1000", "G", "ENST1", "G|0.5")])
    raw = tmp_path / "raw.tab"
    raw.write_text(text, encoding="utf-8")
    proc = _script(str(raw), str(raw))
    assert proc.returncode == 1
    assert proc.stderr == f"{ERROR}the tuple file {raw} is the input\n"
    assert raw.read_text(encoding="utf-8") == text


def test_the_normalizer_output_gives_the_same_tuples_as_the_raw_output(tmp_path: Path) -> None:
    columns = ["FV_REVEL", "FV_GNOMAD_GENE"]
    rows = [
        _row("21:100", "G", "-", "G|0.5", "-", gene="-"),
        _row("21:200", "T", "ENST1", "T|0.25", "GENE1|0.9981|0.2345|3.21|-0.45"),
    ]
    text = _table(columns, rows)
    raw = tmp_path / "raw.tab"
    raw.write_text(text, encoding="utf-8")
    normalized = tmp_path / "normalized.tab"
    subprocess.run([sys.executable, str(NORMALIZER), str(raw), str(normalized)], check=True)
    assert normalized.read_text(encoding="utf-8") == text.replace(
        "\t-\t-\tTranscript\t", "\t-\t-\t-\t"
    )
    from_raw = _script(str(raw), str(tmp_path / "raw.tsv"))
    from_normalized = _script(str(normalized), str(tmp_path / "normalized.tsv"))
    assert from_raw.returncode == 0 and from_normalized.returncode == 0
    assert (tmp_path / "raw.tsv").read_bytes() == (tmp_path / "normalized.tsv").read_bytes()
    assert from_raw.stderr == from_normalized.stderr
    assert from_raw.stderr.startswith("fastvep_plugin_tuples: wrote 4 tuples from 2 rows ")


# The field map


GOOD_ROW = "FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t-"
BAD_MAPS = {
    "header": ("fastvep_column\tpart\n" + GOOD_ROW, ":1: the header is not fastvep_column part"),
    "column-count": (
        MAP_HEADER + "\nFV_X\t1\tX_ONE\tX\ttext\tfabricated\t-",
        ":2: 7 columns, not 8",
    ),
    "part": (
        MAP_HEADER + "\n" + GOOD_ROW.replace("\t1\t", "\tfirst\t"),
        ":2: part 'first' is not N, N-M, N.word or N.score",
    ),
    "part-zero": (
        MAP_HEADER + "\n" + GOOD_ROW.replace("\t1\t", "\t0\t"),
        ":2: part '0' must count from 1",
    ),
    "falling-range": (
        MAP_HEADER + "\n" + GOOD_ROW.replace("\t1\t", "\t3-2\t"),
        ":2: part '3-2' must count from 1, and a range must rise",
    ),
    "kind": (
        MAP_HEADER + "\n" + GOOD_ROW.replace("\ttext\t", "\tfloat\t"),
        ":2: kind 'float' is not text or number",
    ),
    "field-twice": (
        MAP_HEADER + "\n" + GOOD_ROW + "\n" + GOOD_ROW.replace("\t1\t", "\t2\t"),
        ":3: X_ONE is mapped twice",
    ),
    "separator-on-a-range": (
        MAP_HEADER + "\nFV_X\t1-2\tX_ONE\tX\ttext\tfabricated\t|\t-",
        ":2: a list separator needs a part written N, not '1-2'",
    ),
    "separator-on-a-half": (
        MAP_HEADER + "\nFV_X\t1.word\tX_ONE\tX\ttext\tfabricated\t|\t-",
        ":2: a list separator needs a part written N, not '1.word'",
    ),
    "empty-separator": (
        MAP_HEADER + "\nFV_X\t1\tX_ONE\tX\ttext\tfabricated\t\t-",
        ":2: an empty fastvep_column, vep_field, list_separator or reference_set",
    ),
    "empty-reference-set": (
        MAP_HEADER + "\nFV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t",
        ":2: an empty fastvep_column, vep_field, list_separator or reference_set",
    ),
    "blank-line": ("", ":1: the header is not fastvep_column part"),
}


@pytest.mark.parametrize("name", sorted(BAD_MAPS))
def test_a_malformed_field_map_is_refused_with_its_line(tmp_path: Path, name: str) -> None:
    text, message = BAD_MAPS[name]
    custom = tmp_path / "map.tsv"
    custom.write_text(text + "\n", encoding="utf-8")
    run = _one_cell(tmp_path, "FV_X", "G|a", "--field-map", str(custom))
    assert run.returncode == 1
    assert run.stderr.startswith(f"{ERROR}{custom}{message}")
    assert run.stderr.count("\n") == 1
    assert not run.tuple_file.exists()


def test_a_missing_field_map_is_refused(tmp_path: Path) -> None:
    absent = tmp_path / "absent.tsv"
    run = _one_cell(tmp_path, "FV_X", "G|a", "--field-map", str(absent))
    assert run.returncode == 1
    assert run.stderr == (
        f"{ERROR}cannot read the field map {absent}: No such file or directory\n"
    )


def test_a_field_map_that_is_not_utf8_text_is_refused(tmp_path: Path) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_bytes(b"\xff\xfe\n")
    run = _one_cell(tmp_path, "FV_X", "G|a", "--field-map", str(custom))
    assert run.returncode == 1
    assert run.stderr == f"{ERROR}the field map {custom} is not UTF-8 text\n"


def test_an_empty_field_map_is_refused(tmp_path: Path) -> None:
    custom = tmp_path / "map.tsv"
    custom.write_bytes(b"")
    run = _one_cell(tmp_path, "FV_X", "G|a", "--field-map", str(custom))
    assert run.returncode == 1
    assert run.stderr.startswith(f"{ERROR}{custom}:1: the header is not fastvep_column part")


def test_the_shipped_map_covers_each_source_column_and_no_constraint_z_score() -> None:
    rows = _map_rows()
    assert {row["fastvep_column"] for row in rows} == {layout[0] for layout in LAYOUTS.values()}
    assert {row["kind"] for row in rows} == {"text", "number"}
    constraint = [row["part"] for row in rows if row["fastvep_column"] == "FV_GNOMAD_GENE"]
    assert constraint == ["1", "2"]


def _header_fields(reference: Path) -> set[str]:
    names = set()
    with gzip.open(reference, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            if " : " in line:
                names.add(line[3:].split(" : ", 1)[0])
    return names


def _reference_fields(source: str) -> set[str]:
    names = set()
    for path in GOLDEN.glob(f"GRCh38-plugin-{source.lower()}*/expected*/default.txt.gz"):
        names |= _header_fields(path)
    return names


def test_every_reference_set_of_the_map_declares_its_field_on_both_assemblies() -> None:
    rows = [row for row in _map_rows() if row["reference_set"] != "-"]
    assert {row["source"] for row in rows} == {"SpliceAI", "PrimateAI", "dbNSFP"}
    for row in rows:
        for assembly in ("GRCh37", "GRCh38"):
            corpus = f"{assembly}-plugin-{row['source'].lower()}*"
            (reference,) = GOLDEN.glob(f"{corpus}/{row['reference_set']}/default.txt.gz")
            assert row["vep_field"] in _header_fields(reference), (assembly, row["vep_field"])


def test_every_plugin_field_of_the_map_is_a_field_of_a_committed_reference() -> None:
    rows = [row for row in _map_rows() if row["source"] not in CUSTOM_SOURCES]
    assert len({row["source"] for row in rows}) == 7
    for row in rows:
        assert row["vep_field"] in _reference_fields(row["source"]), row["vep_field"]


def test_every_custom_field_of_the_map_is_its_short_name_and_a_committed_info_key() -> None:
    keys = set()
    for assembly in ("GRCh37", "GRCh38"):
        text = (GOLDEN / assembly / "variants.vcf").read_text(encoding="utf-8")
        for line in text.splitlines():
            if not line.startswith("#"):
                keys.update(item.split("=", 1)[0] for item in line.split("\t")[7].split(";"))
    rows = [row for row in _map_rows() if row["source"] in CUSTOM_SOURCES]
    assert len(rows) == 19
    for row in rows:
        short_name, _, key = row["vep_field"].partition("_")
        assert short_name == row["source"] and key in keys, row["vep_field"]
