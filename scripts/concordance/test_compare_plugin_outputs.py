"""Tests for compare_plugin_outputs.py's end-to-end paths.

The comparator's `--json-output` branch aggregates per-plugin counters that the
per-file results carry; a counter the aggregate does not initialise raises KeyError on
the first file and the run dies after printing a passing report. These tests drive the
script as a subprocess, the way the harness does, on two-record fixtures.

The registry tests hold `PLUGIN_FIELDS` to the Ensembl VEP references of the golden
corpora and to the support table of docs/plugins.md: a field a plugin writes that no
entry names is a field the comparison never counts.
"""

from __future__ import annotations

import functools
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts" / "parity"))
import check_parity_matrix  # noqa: E402
import compare_plugin_outputs as cpo  # noqa: E402

GOLDEN = ROOT / "tests" / "golden" / "116"
SCRIPT = Path(__file__).with_name("compare_plugin_outputs.py")
HEADER = "\t".join(
    [
        "#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
        "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "Extra",
    ]
)


def _row(loc: str, allele: str, feature: str, extra: str) -> str:
    return "\t".join(
        ["v", loc, allele, "ENSG1", feature, "Transcript", "missense_variant", "10", "10",
         "4", "A/T", "gcc/Acc", "-", extra]
    )


def _write(path: Path, rows: list[str]) -> Path:
    path.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return path


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True, check=False
    )


@pytest.fixture
def pair(tmp_path: Path) -> tuple[Path, Path]:
    perl = _write(
        tmp_path / "perl.txt",
        [_row("21:1000", "T", "ENST1", "IMPACT=MODERATE;REVEL=0.512"),
         _row("21:2000", "G", "ENST2", "IMPACT=MODERATE;REVEL=0.100")],
    )
    rust = _write(
        tmp_path / "rust.txt",
        [_row("21:1000", "T", "ENST1", "IMPACT=MODERATE;REVEL=0.512"),
         _row("21:2000", "G", "ENST2", "IMPACT=MODERATE;REVEL=0.100")],
    )
    return perl, rust


def test_json_output_aggregates_carry_annotated_counts(pair: tuple[Path, Path], tmp_path: Path) -> None:
    perl, rust = pair
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out))
    assert proc.returncode == 0, proc.stderr
    data = json.loads(out.read_text())
    agg = data["aggregates"]["REVEL"]
    assert agg["perl_annotated"] == 2 and agg["rust_annotated"] == 2
    assert agg["match"] == 2 and agg["mismatch"] == 0
    assert agg["pass"] is True


def test_json_output_flags_a_mismatch_and_a_rust_only_record(tmp_path: Path) -> None:
    perl = _write(tmp_path / "perl.txt",
                  [_row("21:1000", "T", "ENST1", "REVEL=0.512")])
    rust = _write(tmp_path / "rust.txt",
                  [_row("21:1000", "T", "ENST1", "REVEL=0.900"),
                   _row("21:3000", "C", "ENST3", "REVEL=0.300")])
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out))
    assert proc.returncode != 0
    agg = json.loads(out.read_text())["aggregates"]["REVEL"]
    assert agg["mismatch"] == 1 and agg["rust_only"] == 1 and agg["pass"] is False


def test_unknown_plugin_name_is_refused(pair: tuple[Path, Path]) -> None:
    perl, rust = pair
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVELL")
    assert proc.returncode == 2
    assert "unknown plugin" in proc.stderr


# Twelve records whose REVEL value differs between the two sides on every record: two
# more mismatches than the default sample keeps.
PERL_ROWS = [_row(f"21:{1000 + i}", "T", f"ENST{i}", "REVEL=0.100") for i in range(12)]
RUST_ROWS = [_row(f"21:{1000 + i}", "T", f"ENST{i}", "REVEL=0.900") for i in range(12)]


@pytest.mark.parametrize(
    ("value", "kept"),
    [("3", 3), ("11", 11), ("0", 0)],
    ids=["under_the_default", "over_the_default", "zero_keeps_none"],
)
def test_max_mismatches_caps_the_json_sample(tmp_path: Path, value: str, kept: int) -> None:
    perl = _write(tmp_path / "perl.txt", PERL_ROWS)
    rust = _write(tmp_path / "rust.txt", RUST_ROWS)
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out), "--max-mismatches", value)
    assert proc.returncode == 1
    data = json.loads(out.read_text())
    entry = data["results"]["perl"]["REVEL"]
    assert [m["key"][0] for m in entry["mismatches_sample"]] == [
        f"21:{1000 + i}" for i in range(kept)
    ]
    assert entry["total_mismatches"] == 12 and entry["stats"]["mismatch"] == 12
    agg = data["aggregates"]["REVEL"]
    assert agg["mismatch"] == 12 and agg["pass"] is False


def test_negative_max_mismatches_is_refused_by_name(pair: tuple[Path, Path], tmp_path: Path) -> None:
    perl, rust = pair
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out), "--max-mismatches", "-1")
    assert proc.returncode == 2
    assert "argument --max-mismatches" in proc.stderr
    assert not out.exists()


def test_max_mismatches_caps_each_file_pair_in_directory_mode(tmp_path: Path) -> None:
    perl_dir, rust_dir = tmp_path / "perl", tmp_path / "rust"
    perl_dir.mkdir()
    rust_dir.mkdir()
    _write(perl_dir / "first.txt", PERL_ROWS)
    _write(rust_dir / "first.txt", RUST_ROWS)
    _write(perl_dir / "second.txt", PERL_ROWS[:5])
    _write(rust_dir / "second.txt", RUST_ROWS[:5])
    out = tmp_path / "result.json"
    proc = _run("--perl-dir", str(perl_dir), "--rust-dir", str(rust_dir), "--plugins", "REVEL",
                "--json-output", str(out), "--max-mismatches", "4")
    assert proc.returncode == 1
    data = json.loads(out.read_text())
    first, second = data["results"]["first"]["REVEL"], data["results"]["second"]["REVEL"]
    assert len(first["mismatches_sample"]) == 4 and first["total_mismatches"] == 12
    assert len(second["mismatches_sample"]) == 4 and second["total_mismatches"] == 5
    assert data["aggregates"]["REVEL"]["mismatch"] == 17


def test_printed_report_shows_one_mismatch_under_max_mismatches_1(tmp_path: Path) -> None:
    perl = _write(tmp_path / "perl.txt", PERL_ROWS)
    rust = _write(tmp_path / "rust.txt", RUST_ROWS)
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--max-mismatches", "1")
    assert proc.returncode == 1
    assert "Mismatch: 12" in proc.stdout
    shown = [line.strip() for line in proc.stdout.splitlines() if line.endswith(": mismatch")]
    assert shown == ["21:1000 T ENST0: mismatch"]


def test_control_absent_max_mismatches_keeps_ten(tmp_path: Path) -> None:
    """Control: the default sample of ten holds whether or not the option is read."""
    perl = _write(tmp_path / "perl.txt", PERL_ROWS)
    rust = _write(tmp_path / "rust.txt", RUST_ROWS)
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out))
    assert proc.returncode == 1
    entry = json.loads(out.read_text())["results"]["perl"]["REVEL"]
    assert len(entry["mismatches_sample"]) == 10 and entry["total_mismatches"] == 12


def test_gnomadc_fields_are_counted_under_the_names_ensembl_vep_writes(tmp_path: Path) -> None:
    extra = "IMPACT=MODIFIER;gnomAD_10x_cov=0.0332;gnomAD_mean_cov=2.1875"
    other = _row("21:3000", "C", "ENST3", "IMPACT=MODERATE;REVEL=0.100")
    perl = _write(tmp_path / "perl.txt",
                  [_row("21:1000", "T", "ENST1", extra), _row("21:2000", "G", "ENST2", extra),
                   other])
    differing = extra.replace("0.0332", "0.5000") + ";gnomAD_total_DP_cov=7"
    rust = _write(tmp_path / "rust.txt",
                  [_row("21:1000", "T", "ENST1", extra), _row("21:2000", "G", "ENST2", differing),
                   other])
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "gnomADc",
                "--json-output", str(out))
    assert proc.returncode != 0
    data = json.loads(out.read_text())
    agg = data["aggregates"]["gnomADc"]
    assert agg["perl_annotated"] == 2 and agg["rust_annotated"] == 2
    assert agg["match"] == 1 and agg["mismatch"] == 1
    assert data["results"]["perl"]["gnomADc"]["fields"] == [
        "gnomAD_mean_cov", "gnomAD_total_DP_cov", "gnomAD_10x_cov"]


def _write_vcf(path: Path, entries: list[tuple[str, str]]) -> Path:
    header = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
              'Format: Allele|Feature|REVEL|gnomAD_10x_cov">')
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    rows = [f"21\t{pos}\t.\tC\tT\t.\t.\tCSQ={csq}" for pos, csq in entries]
    path.write_text("\n".join([header, columns, *rows]) + "\n")
    return path


def test_vcf_outputs_are_compared_on_the_csq_fields_a_plugin_writes(tmp_path: Path) -> None:
    unannotated = [("3000", "T|ENST3||"), ("4000", "T|ENST4|.|.")]
    perl = _write_vcf(tmp_path / "perl.vcf", [("1000", "T|ENST1|0.512|0.0332"),
                                              ("2000", "T|ENST2||0.0332"), *unannotated])
    rust = _write_vcf(tmp_path / "rust.vcf", [("1000", "T|ENST1|0.512|0.0332"),
                                              ("2000", "T|ENST2||0.5000"), *unannotated])
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust),
                "--plugins", "REVEL,gnomADc", "--json-output", str(out))
    assert proc.returncode != 0
    aggregates = json.loads(out.read_text())["aggregates"]
    revel, gnomadc = aggregates["REVEL"], aggregates["gnomADc"]
    assert revel["perl_annotated"] == 1 and revel["rust_annotated"] == 1 and revel["pass"] is True
    assert gnomadc["perl_annotated"] == 2 and gnomadc["rust_annotated"] == 2
    assert gnomadc["match"] == 1 and gnomadc["mismatch"] == 1


def test_a_name_two_requested_plugins_write_counts_for_one_of_them(tmp_path: Path) -> None:
    rows = [_row("21:1000", "T", "ENST1", "IMPACT=MODERATE;SIFT_pred=T;REVEL_score=0.512"),
            _row("21:2000", "G", "ENST2", "IMPACT=MODERATE;SIFT_pred=D")]
    perl, rust = _write(tmp_path / "perl.txt", rows), _write(tmp_path / "rust.txt", rows)
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust),
                "--plugins", "TranscriptAnnotator,dbNSFP", "--json-output", str(out))
    assert proc.returncode != 0
    data = json.loads(out.read_text())
    dbnsfp, annotator = data["aggregates"]["dbNSFP"], data["aggregates"]["TranscriptAnnotator"]
    assert dbnsfp["perl_annotated"] == 2 and dbnsfp["rust_annotated"] == 2
    assert dbnsfp["pass"] is True
    assert annotator["perl_annotated"] == 0 and annotator["pass"] is False
    assert data["results"]["perl"]["dbNSFP"]["fields"] == ["REVEL_score", "SIFT_pred"]
    assert data["results"]["perl"]["dbNSFP"]["shared_fields"] == ["SIFT_pred"]
    assert data["results"]["perl"]["TranscriptAnnotator"]["fields"] == [
        "transcript_annotator_<column>", "Comment"]
    assert data["results"]["perl"]["TranscriptAnnotator"]["shared_fields"] == ["SIFT_pred"]
    assert proc.stdout.count("Also entries of another requested plugin: SIFT_pred\n") == 2

    proc = _run("--perl-output", str(perl), "--rust-output", str(rust),
                "--plugins", "TranscriptAnnotator", "--json-output", str(out))
    assert proc.returncode == 0
    data = json.loads(out.read_text())
    assert data["aggregates"]["TranscriptAnnotator"]["perl_annotated"] == 2
    assert data["results"]["perl"]["TranscriptAnnotator"]["shared_fields"] == []
    assert "Also entries of another requested plugin" not in proc.stdout


def test_plugin_of_field_answers_with_the_plugin_of_the_run_that_writes_a_field() -> None:
    every = cpo.ALL_PLUGIN_NAMES
    assert cpo.plugin_of_field("gnomAD_30x_cov", ["gnomADc"]) == "gnomADc"
    assert cpo.plugin_of_field("gnomAD_30x_cov", ["REVEL", "gnomADMt"]) is None
    assert cpo.plugin_of_field("gnomAD_median_approx_cov", every) == "gnomADc"
    assert cpo.plugin_of_field("gnomAD_x_cov", every) is None
    assert cpo.plugin_of_field("gnomAD_meanx_cov", every) is None
    assert cpo.plugin_of_field("gnomAD_30x_cov_het", every) is None
    assert cpo.plugin_of_field("am_gnomAD_30x_cov", every) is None
    assert cpo.plugin_of_field("gnomADc_mean", every) is None
    assert cpo.plugin_of_field("gnomAD_AN", every) == "gnomADMt"
    assert cpo.plugin_of_field("gnomAD_exomes_AF", every) == "dbNSFP"
    assert cpo.plugin_of_field("dbNSFP_SIFT_score", every) == "dbNSFP"
    assert cpo.plugin_of_field("dbNSFP_", every) is None
    assert cpo.plugin_of_field("transcript_annotator_Comment", every) == "TranscriptAnnotator"
    assert cpo.plugin_of_field("IMPACT", every) is None
    assert cpo.plugin_of_field("HGVSc", [*every, "HGVS"]) is None
    assert cpo.plugin_of_field("LoF_filter", ["LoF"]) == "LoFTEE"
    assert cpo.plugin_of_field("LoF_filter", ["LoFTEE"]) == "LoFTEE"
    assert cpo.plugin_of_field("LoF_filter", ["REVEL", "LoF", "LoFTEE", "LoF"]) == "LoFTEE"
    assert cpo.plugin_of_field("REVEL", ("REVEL", "REVEL")) == "REVEL"
    assert cpo.plugin_of_field("REVEL", {"CADD", "REVEL"}) == "REVEL"
    assert cpo.plugin_of_field("LoF_filter", frozenset({"LoF"})) == "LoFTEE"


def test_a_name_two_plugins_write_is_settled_by_the_run_and_never_guessed() -> None:
    for name in ("SIFT_pred", "SIFT_score"):
        assert cpo.plugin_of_field(name, ["dbNSFP"]) == "dbNSFP"
        assert cpo.plugin_of_field(name, ["REVEL", "TranscriptAnnotator"]) == "TranscriptAnnotator"
        assert cpo.plugin_of_field(name, ["dbNSFP", "TranscriptAnnotator"]) is None
        assert cpo.plugin_of_field(name, ("TranscriptAnnotator", "dbNSFP")) is None
        assert cpo.plugin_of_field(name, cpo.ALL_PLUGIN_NAMES) is None
        assert cpo.plugin_of_field(name, ["REVEL"]) is None
        assert cpo.plugin_of_field(name, []) is None
    with pytest.raises(TypeError):
        cpo.plugin_of_field("SIFT_score")
    for not_a_list in ("dbNSFP", iter(["TranscriptAnnotator", "dbNSFP"]), {"dbNSFP": 1}, None):
        for name in ("SIFT_score", "IMPACT"):
            with pytest.raises(TypeError):
                cpo.plugin_of_field(name, not_a_list)


def test_field_matcher_reads_several_pattern_entries() -> None:
    claims = cpo.field_matcher(["exact", "a_<N>", "b_<name>_c"])
    assert all(map(claims, ["exact", "a_12", "b_x_c"]))
    assert not any(map(claims, ["exactly", "a_", "a_1x", "a_12b_x_c", "b__c", ""]))


class Reference(NamedTuple):
    """A default-format reference of the golden corpora."""

    name: str  # its path under the corpus root
    plugins: list[str] | None  # the plugins its command line names; None without a command line
    flags: tuple[str, ...]  # the words of its command line
    declared: tuple[tuple[str, str], ...]  # each Extra key its header declares, and its description
    keys: frozenset[str]  # the declared keys and the keys its rows carry


@functools.cache
def _references() -> list[Reference]:
    found = []
    for path in sorted(GOLDEN.glob("*/expected*/default.txt.gz")):
        plugins, flags, declared, keys, declaring = None, (), [], set(), False
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("## VEP command-line:"):
                    plugins = re.findall(r"--plugin\s+([^,\s]+)", line)
                    flags = tuple(line.split(":", 1)[1].split())
                    declaring = False
                elif line.startswith("## Extra column keys:"):
                    declaring = True
                elif line.startswith("##"):
                    if declaring and " : " in line:
                        key, description = line[3:].rstrip("\n").split(" : ", 1)
                        declared.append((key, description))
                        keys.add(key)
                elif not line.startswith("#"):
                    cells = line.rstrip("\n").split("\t")
                    if len(cells) >= 14:
                        keys.update(cpo.parse_extra_field(cells[13]))
        name = str(path.relative_to(GOLDEN))
        found.append(Reference(name, plugins, flags, tuple(declared), frozenset(keys)))
    return found


def _core_keys() -> set[str]:
    """The Extra keys of the references whose command line loads no plugin."""
    return set().union(*(ref.keys for ref in _references() if ref.plugins == []))


def _core_descriptions() -> dict[str, set[str]]:
    """Each core key with the descriptions the references that load no plugin give it."""
    described: dict[str, set[str]] = {}
    for ref in _references():
        for key, description in ref.declared if ref.plugins == [] else ():
            described.setdefault(key, set()).add(description)
    return described


def test_every_plugin_key_of_a_golden_reference_is_registered() -> None:
    core = _core_keys()
    assert core, "no reference without a plugin was read"
    unregistered = []
    for ref in _references():
        if ref.plugins is None:
            assert "-plugin-" not in ref.name, f"{ref.name} has no command line naming its plugins"
            continue
        unregistered += [f"{ref.name}: {key}" for key in sorted(ref.keys - core)
                         if cpo.plugin_of_field(key, ref.plugins, ref.flags) is None]
    assert not unregistered, "not one loaded plugin has an entry for: " + "; ".join(unregistered)
    matchers = [cpo.field_matcher(cpo.PLUGIN_FIELDS[plugin]) for plugin in cpo.ALL_PLUGIN_NAMES]
    claimed = [key for key in sorted(core) if any(claims(key) for claims in matchers)]
    assert not claimed, "an entry names a key a run without plugins writes: " + ", ".join(claimed)


def test_core_named_columns_are_the_ones_the_golden_references_declare() -> None:
    """A core key that a plugin's reference declares under a description no reference without a
    plugin gives it is a plugin's column under that name: the lookup answers the plugin for it
    in that reference's run."""
    described = _core_descriptions()
    columns: dict[str, set[str]] = {}
    for ref in _references():
        for key, description in ref.declared if ref.plugins else ():
            if description not in described.get(key, {description}):
                assert len(ref.plugins) == 1, f"{ref.name} declares {key} among several plugins"
                plugin = cpo.PLUGIN_ALIASES.get(ref.plugins[0], ref.plugins[0])
                columns.setdefault(plugin, set()).add(key)
                assert cpo.plugin_of_field(key, ref.plugins, ref.flags) == plugin
    assert {plugin: sorted(keys) for plugin, keys in columns.items()} == {
        plugin: sorted(flags_of) for plugin, flags_of in cpo.CORE_NAMED_COLUMNS.items()}


def test_the_flags_of_a_core_named_column_are_the_ones_that_turn_its_core_key_on() -> None:
    described = _core_descriptions()
    with_command_line = [ref for ref in _references() if ref.plugins is not None]
    for flags_of in cpo.CORE_NAMED_COLUMNS.values():
        for key, flags in flags_of.items():
            for ref in with_command_line:
                is_core_key = any(k == key and d in described[key] for k, d in ref.declared)
                assert is_core_key == any(flag in ref.flags for flag in flags), ref.name
            for flag in flags:
                assert any(flag in ref.flags for ref in with_command_line), flag


def test_a_run_that_loads_the_plugin_and_turns_the_core_key_on_is_refused() -> None:
    assert cpo.plugin_of_field("APPRIS", ["dbNSFP"]) == "dbNSFP"
    assert cpo.plugin_of_field("TSL", ["REVEL", "dbNSFP"], ["--offline", "--appris"]) == "dbNSFP"
    assert cpo.plugin_of_field("APPRIS", ["dbNSFP"], ("--tsl",)) == "dbNSFP"
    assert cpo.plugin_of_field("APPRIS", ["REVEL"], ["--everything"]) is None
    assert cpo.plugin_of_field("SIFT_score", ["dbNSFP"], ["--everything"]) == "dbNSFP"
    for key, flag in (("APPRIS", "--appris"), ("APPRIS", "--everything"),
                      ("TSL", "--tsl"), ("TSL", "--everything")):
        with pytest.raises(ValueError, match=f"dbNSFP and sets {flag}: {key} is both"):
            cpo.plugin_of_field(key, ["dbNSFP"], ["--offline", flag])
    for not_a_list in ("--appris", iter(["--appris"]), {"--appris": 1}):
        with pytest.raises(TypeError):
            cpo.plugin_of_field("APPRIS", ["dbNSFP"], not_a_list)


def test_the_comparison_leaves_a_core_named_column_to_the_core_key(tmp_path: Path) -> None:
    rows = [_row("21:1000", "T", "ENST1", "IMPACT=MODERATE;APPRIS=principal;TSL=1")]
    perl, rust = _write(tmp_path / "perl.txt", rows), _write(tmp_path / "rust.txt", rows)
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "dbNSFP",
                "--json-output", str(out))
    assert proc.returncode != 0
    agg = json.loads(out.read_text())["aggregates"]["dbNSFP"]
    assert agg["perl_annotated"] == 0 and agg["rust_annotated"] == 0 and agg["pass"] is False


def test_every_registry_entry_is_read_from_a_golden_reference() -> None:
    written: dict[str, set[str]] = {}
    for ref in _references():
        for plugin in ref.plugins or ():
            written.setdefault(cpo.PLUGIN_ALIASES.get(plugin, plugin), set()).update(ref.keys)
    unread = [f"{plugin}: {entry}"
              for plugin in cpo.ALL_PLUGIN_NAMES
              for entry in cpo.PLUGIN_FIELDS[plugin]
              if not any(map(cpo.field_matcher([entry]), written.get(plugin, ())))]
    assert not unread, "no reference that loads the plugin writes: " + "; ".join(unread)


def test_only_the_bare_transcript_annotator_columns_have_two_plugins() -> None:
    matchers = [cpo.field_matcher(cpo.PLUGIN_FIELDS[plugin]) for plugin in cpo.ALL_PLUGIN_NAMES]
    candidates = set().union(*(ref.keys for ref in _references()))
    candidates.update(*(cpo.PLUGIN_FIELDS[plugin] for plugin in cpo.ALL_PLUGIN_NAMES))
    shared = {key for key in candidates if sum(claims(key) for claims in matchers) > 1}
    assert shared == {"SIFT_pred", "SIFT_score"}


def test_the_registry_names_exactly_the_supported_plugins_of_the_support_table() -> None:
    rows, errors = check_parity_matrix.read_support_table(ROOT)
    assert not errors, errors
    supported = {row["name"] for row in rows if row["status"].startswith("supported")}
    assert supported == set(cpo.ALL_PLUGIN_NAMES), (
        "PLUGIN_FIELDS of compare_plugin_outputs.py has one entry per supported plugin of the "
        "support table, naming the Extra keys the plugin's golden references carry")
