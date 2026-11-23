"""reproduce_plugin_concordance.py: its dry runs, its refusals, the rows it prints for one corpus and
the workflow step that counts the tests that score."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shlex
import stat
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts/docs"))
sys.path.insert(0, str(ROOT / "scripts/release"))
import render_release_figure as figure
import render_release_notes as notes
import reproduce_plugin_concordance as driver

CORPUS = "GRCh38-plugin-transcript-annotator"
EXPECTED_ROWS = HERE / "testdata/reproduce_plugin_concordance" / f"{CORPUS}.base.jsonl"
# The workflow whose Test job runs this file's tests that score, and the words by which its step
# selects them.
WORKFLOW = ROOT / ".github/workflows/ci.yml"
COUNTED_RUN = "scripts/concordance/test_reproduce_plugin_concordance.py -k"
ERROR_RE = re.compile(r"^ERROR: \[reproduce_plugin_concordance\] .+$")
# The tree's build of tools/vep-plugin-pairs that the score step runs. A test that scores skips
# without one; CI runs those tests after the build and holds them to have run.
BUILT_TOOL = next((path for path in driver.TOOL_BUILDS if path.is_file()), None)
scores_through_the_built_tool = pytest.mark.skipif(BUILT_TOOL is None, reason="tools/vep-plugin-pairs is not built")
# A stand-in for the vep-rs binary: it writes the file STAND_IN_OUTPUT names to its `-o` path.
STAND_IN = """#!/usr/bin/env python3
import gzip, os, sys
out = sys.argv[sys.argv.index("-o") + 1]
with gzip.open(os.environ["STAND_IN_OUTPUT"], "rb") as source, open(out, "wb") as target:
    target.write(source.read())
"""


def named_inputs(tmp_path: Path) -> list[str]:
    """Every named input a caller gives, as places under one made-up directory."""
    words = []
    for name, option in driver.INPUT_OPTIONS.items():
        words += [option, "4" if name == "fork" else str(tmp_path / "inputs" / name)]
    return words


def dry_run(run: dict, tmp_path: Path, capsys, *extra: str) -> tuple[int, list[str], str]:
    argv = ["--run", run.get("family") or run["corpus"], "--work-dir", str(tmp_path / "work"), "--dry-run", *extra]
    argv += ["--set", run["set"]] if "set" in run else []
    given = named_inputs(tmp_path)
    if "family" not in run:
        at = given.index("--input")
        del given[at : at + 2]
    code = driver.main(argv + given)
    captured = capsys.readouterr()
    return code, captured.out.splitlines(), captured.err


def stand_in_engine(tmp_path: Path, monkeypatch, output: Path) -> Path:
    binary = tmp_path / "stand-in-vep"
    binary.write_text(STAND_IN, encoding="utf-8")
    binary.chmod(binary.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("STAND_IN_OUTPUT", str(output))
    return binary


def refusal(capsys, *argv: str) -> str:
    assert driver.main(list(argv)) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    (line,) = captured.err.splitlines()
    assert ERROR_RE.match(line), line
    return line


def test_the_dry_run_of_every_run_prints_its_commands_with_no_input_left_unnamed(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    assert len(runs) >= 76
    for run in runs:
        code, lines, err = dry_run(run, tmp_path, capsys, "--reference", str(tmp_path / "reference.txt"))
        assert (code, err) == (0, ""), driver.run_key(run)
        assert len(lines) == 2, driver.run_key(run)
        score = lines[1].split()
        assert Path(score[0]) in driver.TOOL_BUILDS and score[1] == "score", lines[1]
        assert f" --class-table {HERE / 'plugin_class_table.tsv'} " in lines[1]
        assert "vep-plugin-pairs" not in lines[0]
        for line in lines:
            assert not driver.PLACEHOLDER_RE.search(line), line
            assert "<" not in line.replace("2>", "").replace(" > ", " "), line


def test_a_dry_run_writes_nothing(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    for run in runs[:6]:
        dry_run(run, tmp_path, capsys, "--reference", str(tmp_path / "reference.txt"))
    assert not (tmp_path / "work").exists()


def test_a_dry_run_prints_the_lines_the_run_executes_and_its_score_line_states_the_binarys_digest(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    binary, tool = tmp_path / "engine", tmp_path / "pairing-tool"
    binary.write_bytes(b"#!/bin/sh\n")
    tool.write_bytes(b"")
    run = ["--run", CORPUS, "--set", "base", "--steps", "vep-rs,score", "--work-dir", str(tmp_path / "work"),
           "--pairing-tool", str(tool)]
    assert driver.main([*run, "--vep", str(binary), "--dry-run"]) == 0
    printed = capsys.readouterr().out.splitlines()
    assert len(printed) == 2
    words = shlex.split(printed[1])
    assert words.count("--engine-binary") == 1
    assert words[words.index("--engine-binary") + 1] == "vep-rs=a8076d3d28d21e02012b20eaf7dbf75409a6277134439025f282e368e3305abf"
    executed = []

    def record(command: dict, step: str, _key: str) -> None:
        executed.append(driver.shell_line(command))
        if step == "score":
            report = Path(command["argv"][command["argv"].index("--report") + 1])
            report.write_text('{"scores": []}', encoding="utf-8")

    monkeypatch.setattr(driver, "execute", record)
    assert driver.main([*run, "--vep", str(binary)]) == 0, capsys.readouterr().err
    assert executed == printed
    # A binary that does not stand has no digest to state, in a dry run as in a run.
    assert driver.main([*run, "--vep", str(tmp_path / "no-such-binary"), "--dry-run"]) == 0
    assert "--engine-binary" not in capsys.readouterr().out


def test_a_binary_that_stands_and_cannot_be_read_is_refused_in_one_line_in_a_dry_run_as_in_a_run(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    binary = tmp_path / "engine"
    binary.write_bytes(b"")

    def unreadable(path: Path) -> str:
        raise PermissionError(13, "Permission denied", str(path))

    monkeypatch.setattr(driver, "sha256_of", unreadable)
    score = ["--run", CORPUS, "--set", "base", "--steps", "score", "--vep", str(binary), "--work-dir", str(tmp_path / "work")]
    for mode in (["--dry-run"], []):
        line = refusal(capsys, *score, *mode)
        assert line == f"ERROR: [reproduce_plugin_concordance] the binary {binary} cannot be read: Permission denied"
    assert not (tmp_path / "work").exists()


def test_a_reference_by_path_prints_no_container_command(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    for run in runs:
        _code, lines, _err = dry_run(run, tmp_path, capsys, "--reference", str(tmp_path / "reference.txt"))
        assert not any("docker" in line or "generate_reference_output.sh" in line for line in lines), driver.run_key(run)
        assert f"--reference {tmp_path / 'reference.txt'}" in lines[1]


def test_a_whole_run_runs_ensembl_vep_from_the_pinned_digest_with_its_recorded_words(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    for run in (run for run in runs if "family" in run):
        code, lines, _err = dry_run(run, tmp_path, capsys, "--run-ensembl-vep", "--tool-dir", str(tmp_path / "tools"))
        assert code == 0 and len(lines) == 3, run["family"]
        container = lines[0].split(" > ")[0].split()
        assert container[:3] == ["docker", "run", "--rm"]
        assert driver.pinned_image() in container
        assert container[container.index("--fork") + 1] == "16"
        inside = container[container.index("vep") :]
        assert inside[:5] == ["vep", "-i", "/work/input.vcf", "-o", "/work/perl_out/output.txt"]
        assert not any(str(tmp_path) in word for word in inside), lines[0]
        plugins = [inside[i + 1] for i, word in enumerate(inside) if word == "--plugin"]
        ours = [run["vep_rs"][i + 1] for i, word in enumerate(run["vep_rs"]) if word == "--plugin"]
        assert plugins == [word.replace("{plugin_data}", "/work/plugin_data").replace("{loftee}", "/work/loftee") for word in ours]
        assert str(tmp_path / "work" / run["family"] / "ensembl-vep" / "output.txt") in lines[2]
    loftee = dry_run(next(run for run in runs if run.get("family") == "loftee_grch38"), tmp_path, capsys,
                     "--run-ensembl-vep", "--tool-dir", str(tmp_path / "tools"))[1][0]
    assert f"{tmp_path / 'tools'}:/work/tools/0:ro" in loftee and "/work/tools/0 vep" in loftee
    assert "/work/loftee:ro" in loftee and "loftee_path:/work/loftee" in loftee


def test_the_fork_count_and_the_environment_of_ensembl_vep_follow_the_callers_words(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    gwas = next(run for run in runs if run.get("family") == "gwas")
    _code, lines, _err = dry_run(gwas, tmp_path, capsys, "--run-ensembl-vep", "--ensembl-vep-fork", "4",
                                 "--perl-env", "PERL_HASH_SEED=0", "--perl-env", "PERL_PERTURB_KEYS=0")
    container = lines[0].split()
    assert container[container.index("--fork") + 1] == "4"
    assert "-e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0" in lines[0]
    assert lines[1].split()[lines[1].split().index("--fork") + 1] == "4"
    assert ":/work/dir " in lines[0] and "--dir /work/dir" in lines[0]


def test_a_corpus_set_runs_ensembl_vep_through_the_generator_on_a_copy(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    run = next(run for run in runs if run.get("corpus") == "GRCh38-plugin-bayesdel-fixture" and run["set"] == "bayesdel")
    code, lines, _err = dry_run(run, tmp_path, capsys, "--run-ensembl-vep")
    copy = tmp_path / "work" / driver.run_key(run) / "corpus"
    assert code == 0 and len(lines) == 4
    assert lines[0].split()[:2] == ["cp", "-R"] and lines[0].endswith(str(copy))
    assert f"generate_reference_output.sh --corpus {copy} " in lines[1] and "--flag-set bayesdel" in lines[1]
    assert f"--reference {copy}/expected_bayesdel/default.txt.gz" in lines[3]
    base = next(run for run in runs if run.get("corpus") == "GRCh38-plugin-bayesdel-fixture" and run["set"] == "base")
    assert "--flag-set" not in dry_run(base, tmp_path, capsys, "--run-ensembl-vep")[1][1]


def test_one_line_of_the_generator_pins_the_image_and_this_script_states_no_digest() -> None:
    generator = driver.GENERATOR.read_text(encoding="utf-8")
    assert len(re.findall(r"sha256:[0-9a-f]{64}", generator)) == 1
    assert not re.search(r"[0-9a-f]{64}", Path(driver.__file__).read_text(encoding="utf-8"))
    assert re.fullmatch(r"ensemblorg/ensembl-vep@sha256:[0-9a-f]{64}", driver.pinned_image())


def test_each_refusal_ends_in_its_one_line_and_exit_1(tmp_path: Path, capsys) -> None:
    work = ["--work-dir", str(tmp_path / "work")]
    assert "holds no entry for the run nosuch" in refusal(capsys, "--run", "nosuch", *work)
    assert "holds no entry for the row GWAS on GRCh37" in refusal(capsys, "--row", "GWAS", "--assembly", "GRCh37", *work)
    assert "--row takes --assembly" in refusal(capsys, "--row", "GWAS", *work)
    assert "name one of --run, --row or --list" in refusal(capsys, *work)
    assert "pass --work-dir" in refusal(capsys, "--run", "gwas")
    assert "pass --input" in refusal(capsys, "--run", "gwas", *work)
    assert "pass --reference or --run-ensembl-vep" in refusal(capsys, "--run", "gwas", "--input", "in.vcf", *work)
    assert "the steps are" in refusal(capsys, "--run", "gwas", "--steps", "score,all", *work)
    assert "pass one" in refusal(capsys, "--run", CORPUS, "--reference", "r.txt", "--run-ensembl-vep", *work)
    assert "select one run with --set" in refusal(capsys, "--run", CORPUS, "--report", "r.json", *work)
    assert "--input is for a whole run" in refusal(capsys, "--run", CORPUS, "--input", "in.vcf", *work)
    assert "--perl-env takes NAME=VALUE" in refusal(capsys, "--run", CORPUS, "--perl-env", "SEED", *work)
    line = refusal(capsys, "--run", "gwas", "--input", "in.vcf", "--reference", "r.txt", "--dry-run", *work)
    assert "names {vep_rs_cache}: pass --vep-rs-cache" in line
    line = refusal(capsys, "--run", "GRCh38-plugin-maxentscan-fixture", "--set", "base", "--run-ensembl-vep",
                   "--generated", "g", "--ensembl-vep-cache", "c", "--plugin-modules", "m", "--dry-run", *work)
    assert "its reference is the corpus's committed expected output" in line
    assert not (tmp_path / "work").exists()


def test_a_missing_data_file_is_refused_by_name_before_anything_runs(tmp_path: Path, capsys, monkeypatch) -> None:
    marker = tmp_path / "ran"
    binary = tmp_path / "vep"
    binary.write_text(f"#!/bin/sh\ntouch {marker}\n", encoding="utf-8")
    binary.chmod(0o755)
    places = {name: tmp_path / name for name in ("cache", "plugin_data")}
    for place in places.values():
        place.mkdir()
    for name in ("in.vcf", "ref.fa", "ref.fa.fai", "reference.txt"):
        (tmp_path / name).write_text("", encoding="utf-8")
    line = refusal(
        capsys, "--run", "gwas", "--work-dir", str(tmp_path / "work"), "--vep", str(binary), "--pairing-tool", str(binary),
        "--input", str(tmp_path / "in.vcf"), "--reference", str(tmp_path / "reference.txt"),
        "--vep-rs-cache", str(places["cache"]), "--fasta", str(tmp_path / "ref.fa"), "--plugin-data", str(places["plugin_data"]),
    )
    assert "reads 1 file(s) that are missing" in line
    assert line.endswith(str(places["plugin_data"] / "assembly_agnostic/gwas/gwas_catalog_associations.tsv"))
    assert not marker.exists() and not (tmp_path / "work").exists()


def test_a_missing_index_is_refused_by_name(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    run = next(run for run in runs if run.get("family") == "loftee_grch38")
    args = driver.parse(["--run", "loftee_grch38", "--work-dir", str(tmp_path / "work"), "--reference", "r.txt",
                         *named_inputs(tmp_path)])
    missing = driver.Plan(run, args, _inputs).missing(("check", "vep-rs"))
    data = tmp_path / "inputs/plugin_data/grch38/loftee/grch38"
    assert str(data / "human_ancestor.fa") in missing and str(data / "human_ancestor.fa.fai") in missing


@scores_through_the_built_tool
def test_the_rows_the_built_tool_scores_for_one_corpus_equal_the_committed_ones(tmp_path: Path, capsys, monkeypatch) -> None:
    reference = ROOT / "tests/golden/116" / CORPUS / "expected/default.txt.gz"
    binary = stand_in_engine(tmp_path, monkeypatch, reference)
    report = tmp_path / "report.json"
    rows_out = tmp_path / "rows.jsonl"
    code = driver.main(["--run", CORPUS, "--set", "base", "--vep", str(binary), "--work-dir", str(tmp_path / "work"),
                        "--report", str(report), "--rows-out", str(rows_out)])
    captured = capsys.readouterr()
    assert code == 0, captured.err
    assert captured.out == EXPECTED_ROWS.read_text(encoding="utf-8") == rows_out.read_text(encoding="utf-8")
    assert captured.err.startswith(f"TranscriptAnnotator GRCh38 on {CORPUS}: raw F1 1.000000000, adjusted F1 1.000000000;")
    scored = json.loads(report.read_text(encoding="utf-8"))
    assert scored["class_table"]["path"] == str(HERE / "plugin_class_table.tsv")
    assert scored["run_command"][0] == str(binary)
    assert not (tmp_path / "work" / CORPUS / "base" / "scratch").exists()


def test_the_readers_of_the_plugin_record_take_the_printed_rows() -> None:
    entries = [json.loads(line) for line in EXPECTED_ROWS.read_text(encoding="utf-8").splitlines()]
    assert [entry["projection"] for entry in entries] == ["feature", "allele"]
    record = {"file": "printed rows", "data": {"real_data": [], "synthetic_fixtures": [], "committed_table": entries}}
    index = figure.plugin_entries(record)
    for entry in entries:
        (found,) = index[("committed_table", "TranscriptAnnotator", "GRCh38", CORPUS, "vep-rs", entry["projection"])]
        where = found[1]
        assert figure.f1(*figure.summed_counts([found], figure.RAW_COUNTS), where) == entry["raw_f1"]
        adjusted = figure.summed_counts([found], figure.ADJUSTED_COUNTS)
        assert figure.released_adjusted_f1(*adjusted, where, "vep-rs") == entry["adjusted_f1"]
    table, scored = notes.plugin_table({"rows": [entries[0]]})
    assert scored == 1 and "| TranscriptAnnotator |" in table[-1] and table[-1].count("1.000000000") == 2


def test_a_printed_row_carries_the_keys_of_a_record_entry_in_the_records_order() -> None:
    fixture = json.loads((ROOT / "scripts/docs/testdata/plugin_record_fixture.json").read_text(encoding="utf-8"))
    printed = json.loads(EXPECTED_ROWS.read_text(encoding="utf-8").splitlines()[0])
    wanted = [key for key in fixture["committed_table"][0] if key != "run_manifest_sha256"]
    assert list(printed) == [*wanted, "raw_f1", "adjusted_f1"]
    whole = driver.record_entries(
        {"family": "gwas", "rows": [["GWAS", "GRCh38", "clinvar", "clinvar_grch38"]]},
        {"scores": [{"engine": "vep-rs", "source": "GWAS", "projection": "feature", "perl": 9, "rust": 8, "intersection": 7,
                     "adj_perl": 7, "adj_rust": 7, "adj_intersection": 7, "masked": {"a_class": 3}, "residual_rows": 0,
                     "text_only_differences": 0}]},
    )
    assert list(whole[0]) == [*(key for key in fixture["real_data"][0] if key != "run_manifest_sha256"), "raw_f1", "adjusted_f1"]
    assert (whole[0]["raw_f1"], whole[0]["adjusted_f1"]) == ("0.823529412", "1.000000000")


@pytest.mark.parametrize("counts", [(7, 9, 8), (0, 3, 0), (5, 5, 5), (10**9 - 1, 10**9, 10**9)])
def test_the_f1_is_the_record_readers_on_counts_both_print(counts: tuple[int, int, int]) -> None:
    assert driver.f1_nine(*counts, exact=False) == figure.f1(*counts, "a row")
    assert driver.f1_nine(*counts, exact=True) == figure.released_adjusted_f1(*counts, "a row", "vep-rs")
    assert driver.f1_nine(0, 0, 0, exact=False) is None


def test_unequal_counts_that_round_to_one_read_below_one_in_the_driver_and_are_refused_by_the_record_reader() -> None:
    matched, reference, engine = counts = (10**12, 10**12 + 1, 10**12)
    assert driver.f1_nine(*counts, exact=False) == figure.f1(*counts, "a row") == "1.000000000"
    assert driver.f1_nine(*counts, exact=True) == "0.999999999"
    with pytest.raises(SystemExit) as refused:
        figure.released_adjusted_f1(*counts, "a row", "vep-rs")
    assert refused.value.code == (
        "ERROR: [render_release_figure] a row adjusted counts are unequal and their F1 rounds to 1 at nine places "
        f"(reference {reference}, vep-rs {engine}, matched {matched}); the released version's adjusted F1 is exactly 1 "
        "or the record is not a shipped state"
    )


def test_a_label_in_parentheses_is_a_row_of_its_plugin() -> None:
    runs, _inputs = driver.load_runs()
    run = next(run for run in runs if run.get("corpus") == "GRCh38-plugin-cadd-sv-fixture" and run["set"] == driver.BASE_SET)
    assert run["rows"] == [["CADD (structural variants)", "GRCh38", "synthetic_fixture", "GRCh38-plugin-cadd-sv-fixture"]]
    score = {"engine": "vep-rs", "source": "CADD", "projection": "feature", "perl": 1, "rust": 1, "intersection": 1,
             "adj_perl": 1, "adj_rust": 1, "adj_intersection": 1, "masked": {}, "residual_rows": 0, "text_only_differences": 0}
    (entry,) = driver.record_entries(run, {"scores": [score, {**score, "engine": "another"}]})
    assert (entry["plugin"], entry["expected_set"]) == ("CADD (structural variants)", "base")


def test_a_score_with_a_declared_class_not_applied_is_refused_in_one_line(tmp_path: Path, capsys, monkeypatch) -> None:
    score = {"engine": "vep-rs", "source": "GWAS", "projection": "feature", "perl": 9, "rust": 8, "intersection": 7,
             "masked": {"a_class": 1}, "not_applied": {"b_class": 2, "c_class": 1}, "text_only_differences": 0, "residual_rows": 3,
             "shared_keys": 7}
    corpus = {"corpus": "GRCh38-plugin-gwas", "set": "base", "rows": [["GWAS", "GRCh38", "synthetic_fixture", "GRCh38-plugin-gwas"]]}
    with pytest.raises(driver.Refusal) as refused:
        driver.record_entries(corpus, {"scores": [score]})
    told = "carries no adjusted count: declared and not applied: b_class (2), c_class (1)"
    assert str(refused.value) == f"the score of GWAS on GRCh38-plugin-gwas/base {told}"
    whole = {"family": "gwas", "rows": [["GWAS", "GRCh38", "clinvar", "clinvar_grch38"]]}
    with pytest.raises(driver.Refusal) as refused:
        driver.record_entries(whole, {"scores": [score]}, "a_population_dataset")
    assert str(refused.value) == f"the score of GWAS on a_population_dataset (gwas) {told}"
    scored = {**score, "source": "TranscriptAnnotator"}
    monkeypatch.setattr(driver, "run_plan", lambda plan, steps: driver.record_entries(plan.run, {"scores": [scored]}, plan.dataset))
    line = refusal(capsys, "--run", CORPUS, "--set", "base", "--vep", str(tmp_path / "vep"), "--work-dir", str(tmp_path / "work"))
    assert line == f"ERROR: [reproduce_plugin_concordance] the score of TranscriptAnnotator on {CORPUS}/base {told}"


def test_a_score_with_its_adjusted_counts_writes_its_row() -> None:
    score = {"engine": "vep-rs", "source": "GWAS", "projection": "feature", "perl": 9, "rust": 8, "intersection": 7,
             "adj_perl": 7, "adj_rust": 7, "adj_intersection": 7, "masked": {"a_class": 3}, "not_applied": {},
             "text_only_differences": 0, "residual_rows": 0, "shared_keys": 7}
    corpus = {"corpus": "GRCh38-plugin-gwas", "set": "base", "rows": [["GWAS", "GRCh38", "synthetic_fixture", "GRCh38-plugin-gwas"]]}
    (entry,) = driver.record_entries(corpus, {"scores": [score]})
    assert json.dumps(entry) == (
        '{"plugin": "GWAS", "assembly": "GRCh38", "corpus": "GRCh38-plugin-gwas", "expected_set": "base", "engine": "vep-rs", '
        '"projection": "feature", "scored_on": "synthetic_fixture", "perl": 9, "rust": 8, "intersection": 7, "adj_perl": 7, '
        '"adj_rust": 7, "adj_intersection": 7, "masked": {"a_class": 3}, "residual_rows": 0, "text_only_differences": 0, '
        '"raw_f1": "0.823529412", "adjusted_f1": "1.000000000"}'
    )


def test_the_score_step_runs_the_named_tool_else_the_trees_release_build_else_its_debug_build(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    builds = (tmp_path / "release/vep-plugin-pairs", tmp_path / "debug/vep-plugin-pairs")
    monkeypatch.setattr(driver, "TOOL_BUILDS", builds)
    argv = ["--run", CORPUS, "--set", "base", "--work-dir", str(tmp_path / "work"), "--dry-run"]

    def score_words(*extra: str) -> list[str]:
        assert driver.main([*argv, *extra]) == 0
        return capsys.readouterr().out.splitlines()[-1].split()

    # A dry run checks no file: with no build it names the release build's place.
    assert score_words()[:2] == [str(builds[0]), "score"]
    for build in reversed(builds):
        build.parent.mkdir()
        build.write_text("", encoding="utf-8")
        assert score_words()[:2] == [str(build), "score"]
    named = tmp_path / "named-tool"
    words = score_words("--pairing-tool", str(named), "--pairing-threads", "3", "--class-table", str(tmp_path / "classes.tsv"))
    assert words[:2] == [str(named), "score"] and "--pairing-tool" not in words
    assert words[words.index("--pairing-threads") + 1] == "3"
    assert words[words.index("--class-table") + 1] == str(tmp_path / "classes.tsv") and words.count("--class-table") == 1
    assert "--pairing-threads" not in score_words()


def test_a_score_with_no_tool_built_is_refused_before_any_step_with_the_command_that_builds_it(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    marker = tmp_path / "ran"
    binary = tmp_path / "vep"
    binary.write_text(f"#!/bin/sh\ntouch {marker}\n", encoding="utf-8")
    binary.chmod(0o755)
    unbuilt = tmp_path / "release/vep-plugin-pairs"
    monkeypatch.setattr(driver, "TOOL_BUILDS", (unbuilt, tmp_path / "debug/vep-plugin-pairs"))
    argv = ["--run", CORPUS, "--set", "base", "--vep", str(binary), "--work-dir", str(tmp_path / "work")]
    build = f"'cargo build --release --locked --manifest-path {ROOT / 'tools/vep-plugin-pairs/Cargo.toml'}'"
    told = f"which is not a file: build the tool with {build} or name a built vep-plugin-pairs with --pairing-tool"
    line = refusal(capsys, *argv)
    assert line == f"ERROR: [reproduce_plugin_concordance] the score step of {CORPUS}/base runs {unbuilt}, {told}"
    line = refusal(capsys, *argv, "--pairing-tool", str(tmp_path / "no-such-tool"))
    assert line == f"ERROR: [reproduce_plugin_concordance] the score step of {CORPUS}/base runs {tmp_path / 'no-such-tool'}, {told}"
    assert not marker.exists() and not (tmp_path / "work").exists()
    # A run that names no score step needs no tool.
    assert driver.main([*argv, "--steps", "check,vep-rs"]) == 0
    assert marker.exists()


@scores_through_the_built_tool
def test_the_score_step_alone_runs_the_built_tool_on_an_engine_output_by_path(tmp_path: Path, capsys, monkeypatch) -> None:
    reference = ROOT / "tests/golden/116" / CORPUS / "expected_cols/default.txt.gz"
    output = tmp_path / "held.txt"
    stand_in_engine(tmp_path, monkeypatch, reference)
    line = refusal(capsys, "--run", CORPUS, "--set", "cols", "--steps", "check,score", "--engine-output", str(output),
                   "--work-dir", str(tmp_path / "work"))
    assert line.endswith(str(output))
    output.write_bytes(gzip.decompress(reference.read_bytes()))
    code = driver.main(["--run", CORPUS, "--set", "cols", "--steps", "check,score", "--engine-output", str(output),
                        "--work-dir", str(tmp_path / "work"), "--vep", str(tmp_path / "no-such-binary")])
    captured = capsys.readouterr()
    assert code == 0, captured.err
    assert [json.loads(line)["expected_set"] for line in captured.out.splitlines()] == ["cols", "cols"]


def test_a_step_that_ends_non_zero_is_refused_with_its_error_file(tmp_path: Path, capsys) -> None:
    binary = tmp_path / "vep"
    binary.write_text("#!/bin/sh\necho no >&2\nexit 3\n", encoding="utf-8")
    binary.chmod(0o755)
    line = refusal(capsys, "--run", CORPUS, "--set", "base", "--vep", str(binary), "--pairing-tool", str(binary),
                   "--work-dir", str(tmp_path / "work"))
    assert f"the vep-rs step of {CORPUS}/base ended 3" in line
    assert (tmp_path / "work" / CORPUS / "base" / "vep-rs.stderr").read_text(encoding="utf-8") == "no\n"


@scores_through_the_built_tool
def test_a_relative_path_to_the_binary_is_read_from_the_working_directory_and_scored_by_the_built_tool(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    reference = ROOT / "tests/golden/116" / CORPUS / "expected/default.txt.gz"
    binary = stand_in_engine(tmp_path, monkeypatch, reference)
    monkeypatch.chdir(tmp_path)
    code = driver.main(["--run", CORPUS, "--set", "base", "--vep", f"./{binary.name}", "--work-dir", "work"])
    captured = capsys.readouterr()
    assert code == 0, captured.err
    assert captured.out == EXPECTED_ROWS.read_text(encoding="utf-8")
    report = json.loads((tmp_path / "work" / CORPUS / "base" / "report.json").read_text(encoding="utf-8"))
    assert report["run_command"][0] == str(binary)


def counting_step() -> list[str]:
    """The commands of the workflow step that runs this file's tests that score and counts them, as
    the workflow writes them under the step's `run: |`."""
    lines = WORKFLOW.read_text(encoding="utf-8").splitlines()
    named = [at for at, line in enumerate(lines) if COUNTED_RUN in line and not line.lstrip().startswith("#")]
    assert len(named) == 1, f"{len(named)} lines of {WORKFLOW.name} run this file's tests that score"
    opened = max(at for at in range(named[0]) if re.fullmatch(r"\s*(?:- )?run: \|", lines[at]))
    indent = len(lines[opened + 1]) - len(lines[opened + 1].lstrip())
    block: list[str] = []
    for line in lines[opened + 1 :]:
        if line.strip() and len(line) - len(line.lstrip()) < indent:
            break
        block.append(line[indent:])
    return block


def two_runs_of_the_counting_step(tmp_path: Path, first: str, second: str) -> tuple[int, int, str]:
    """Runs the workflow's step twice on one machine and returns both exit statuses and the first
    run's standard error. Each run has a `RUNNER_TEMP` of its own and a `python3` of its own, which
    prints the summary line given for that run in place of pytest and the committed rows in place
    of the driver. The second run starts once the first has written its summary and ends before
    the first counts it, the order in which a shared scratch file gives the first run the second
    run's summary."""
    block = counting_step()
    summary_written = next(at for at, line in enumerate(block) if "-m pytest" in line)

    def run_of(name: str, summary: str, commands: list[str]) -> list[str]:
        home = tmp_path / name
        (home / "bin").mkdir(parents=True)
        (home / "runner-temp").mkdir()
        stand_in = home / "bin/python3"
        stand_in.write_text(
            "#!/bin/sh\n"
            'case " $* " in\n'
            f"  *' -m pytest '*) printf '%s\\n' {shlex.quote(summary)} ;;\n"
            f"  *) cat {shlex.quote(str(EXPECTED_ROWS))} ;;\n"
            "esac\n",
            encoding="utf-8",
        )
        stand_in.chmod(0o755)
        script = home / "step.sh"
        script.write_text("\n".join(commands) + "\n", encoding="utf-8")
        path = f"{home / 'bin'}{os.pathsep}{os.environ['PATH']}"
        return ["env", f"RUNNER_TEMP={home / 'runner-temp'}", f"PATH={path}", "bash", str(script)]

    second_status = tmp_path / "second.status"
    second_run = shlex.join(run_of("second", second, block))
    meanwhile = f"status=0; {second_run} || status=$?; echo \"$status\" > {shlex.quote(str(second_status))}"
    first_run = run_of("first", first, [*block[: summary_written + 1], meanwhile, *block[summary_written + 1 :]])
    done = subprocess.run(first_run, cwd=ROOT, capture_output=True, text=True, check=False)
    return done.returncode, int(second_status.read_text(encoding="utf-8")), done.stderr


def test_the_workflow_step_that_counts_the_scoring_tests_reads_the_summary_of_its_own_run(tmp_path: Path) -> None:
    ran = "3 passed in 0.01s"
    skipped = "3 skipped in 0.01s"
    # A run whose three tests skipped fails its count beside a run whose three tests ran.
    first, second, stderr = two_runs_of_the_counting_step(tmp_path / "skipped-beside-ran", skipped, ran)
    assert (first, second) == (1, 0), stderr
    # A run whose three tests ran passes its count beside a run whose three tests skipped.
    first, second, stderr = two_runs_of_the_counting_step(tmp_path / "ran-beside-skipped", ran, skipped)
    assert (first, second) == (0, 1), stderr


def test_a_program_that_cannot_start_is_refused_in_one_line(tmp_path: Path, capsys) -> None:
    binary = tmp_path / "vep"
    binary.write_text("not a program\n", encoding="utf-8")
    line = refusal(capsys, "--run", CORPUS, "--set", "base", "--vep", str(binary), "--steps", "vep-rs",
                   "--work-dir", str(tmp_path / "work"))
    assert f"the vep-rs step of {CORPUS}/base cannot start {binary}: " in line


def test_a_whole_genome_cell_names_its_own_dataset_in_the_printed_rows(tmp_path: Path, capsys) -> None:
    run = {"family": "gwas", "dataset": "clinvar_grch38", "rows": [["GWAS", "GRCh38", "clinvar", "clinvar_grch38"]]}
    score = {"engine": "vep-rs", "source": "GWAS", "projection": "feature", "perl": 1, "rust": 1, "intersection": 1,
             "adj_perl": 1, "adj_rust": 1, "adj_intersection": 1, "masked": {}, "residual_rows": 0, "text_only_differences": 0}
    assert driver.record_entries(run, {"scores": [score]})[0]["dataset"] == "clinvar_grch38"
    assert driver.record_entries(run, {"scores": [score]}, "a_population_dataset")[0]["dataset"] == "a_population_dataset"
    runs, inputs = driver.load_runs()
    args = driver.parse(["--run", "gwas", "--dataset", "a_population_dataset", "--work-dir", str(tmp_path), "--reference", "r.txt",
                         *named_inputs(tmp_path)])
    assert driver.Plan(next(run for run in runs if run.get("family") == "gwas"), args, inputs).dataset == "a_population_dataset"
    line = refusal(capsys, "--run", CORPUS, "--dataset", "a_population_dataset", "--work-dir", str(tmp_path / "work"))
    assert "--dataset names a whole run's input" in line


def test_list_prints_every_run_with_its_rows(capsys) -> None:
    assert driver.main(["--list"]) == 0
    lines = capsys.readouterr().out.splitlines()
    runs, _inputs = driver.load_runs()
    assert len(lines) == len(runs)
    assert "gwas\tclinvar\tGWAS GRCh38" in lines
    assert f"{CORPUS}/base\tcommitted_table\tTranscriptAnnotator GRCh38" in lines


def test_a_row_selects_the_run_or_the_sets_it_stands_on(tmp_path: Path, capsys) -> None:
    code = driver.main(["--row", "BayesDel", "--assembly", "GRCh38", "--work-dir", str(tmp_path / "work"), "--dry-run"])
    lines = capsys.readouterr().out.splitlines()
    assert code == 0 and len(lines) == 4
    assert sum("expected_bayesdel/default.txt.gz" in line for line in lines) == 1
    code = driver.main(["--row", "GWAS", "--assembly", "GRCh38", "--work-dir", str(tmp_path / "work"), "--dry-run",
                        "--reference", "r.txt", *named_inputs(tmp_path)])
    assert code == 0 and "--plugins GWAS" in capsys.readouterr().out


SWEEP_STEP = "utrannotator-sweep"
SWEEP_RUNS = ("plugins_grch37", "plugins_grch38")
# The names Docker takes for a container.
DOCKER_NAME_RE = re.compile(r"[a-zA-Z0-9][a-zA-Z0-9_.-]+")
# A stand-in for docker: it appends its words to the file STAND_IN_LOG names and, unless
# STAND_IN_EXIT names a status to end on, writes the seed of the run to the file its `-o` word
# names, under the directory mounted at /work/out. Where STAND_IN_RUNNING names a directory, every
# container it started stays running there under its `--name`, and a second container of a name
# that stands there ends 125, as Docker ends one whose name is in use.
STAND_IN_DOCKER = """#!/bin/sh
printf '%s\\n' "$*" >> "$STAND_IN_LOG"
if [ "${STAND_IN_EXIT:-0}" != 0 ]; then
    echo no >&2
    exit "$STAND_IN_EXIT"
fi
while [ $# -gt 0 ]; do
    case "$1" in
    *:/work/out) out="${1%:/work/out}" ;;
    PERL_HASH_SEED=*) seed="${1#PERL_HASH_SEED=}" ;;
    -o) name="${2##*/}" ;;
    --name) container="$2" ;;
    esac
    shift
done
if [ -n "${STAND_IN_RUNNING:-}" ] && ! mkdir "$STAND_IN_RUNNING/$container" 2>/dev/null; then
    echo "the container name $container is already in use" >&2
    exit 125
fi
printf '%s\\n' "$seed" > "$out/$name"
"""


def sweep_seeds(run: dict) -> list[str]:
    """The seeds a run's `score` words name, in their order: the text between `sweep-seed` and `.txt`."""
    return [word.split("/sweep-seed")[1].removesuffix(".txt") for word in run["score"] if "/sweep-seed" in word]


def combined_run(family: str) -> dict:
    return next(run for run in driver.load_runs()[0] if run.get("family") == family)


def stand_in_docker(tmp_path: Path, monkeypatch) -> Path:
    """A stand-in docker first on PATH; the file its calls are logged in."""
    programs = tmp_path / "programs"
    programs.mkdir()
    (programs / "docker").write_text(STAND_IN_DOCKER, encoding="utf-8")
    (programs / "docker").chmod(0o755)
    log = tmp_path / "docker.log"
    log.write_text("", encoding="utf-8")
    monkeypatch.setenv("PATH", f"{programs}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("STAND_IN_LOG", str(log))
    return log


def sweep_argv(tmp_path: Path, family: str, *extra: str) -> list[str]:
    return ["--run", family, "--steps", SWEEP_STEP, "--work-dir", str(tmp_path / "work"),
            "--reference", str(tmp_path / "reference.txt"), *named_inputs(tmp_path), *extra]


def test_the_sweep_step_stands_in_the_steps_order_and_the_usage_text_and_runs_only_when_named(tmp_path: Path, capsys) -> None:
    assert driver.STEPS == ("check", "ensembl-vep", SWEEP_STEP, "vep-rs", "score")
    assert driver.parse([]).steps == "check,ensembl-vep,vep-rs,score"
    usage = driver.__doc__
    assert "[--steps check,ensembl-vep,utrannotator-sweep,vep-rs,score]" in usage
    assert "- `utrannotator-sweep`, only where `--steps` names it:" in usage
    assert "[--sweep-input FILE]" in usage and "[--overwrite-sweep]" in usage
    line = refusal(capsys, "--run", "gwas", "--steps", "score,all", "--work-dir", str(tmp_path / "work"))
    assert line.endswith("the steps are check, ensembl-vep, utrannotator-sweep, vep-rs, score")
    for family in SWEEP_RUNS:
        _code, lines, _err = dry_run(combined_run(family), tmp_path, capsys, "--run-ensembl-vep")
        assert len(lines) == 3 and not any("PERL_HASH_SEED" in line for line in lines), family


def test_the_sweeps_dry_run_prints_one_seeded_command_a_seed_in_the_order_of_the_score_words(tmp_path: Path, capsys) -> None:
    inputs = tmp_path / "inputs"
    for family in SWEEP_RUNS:
        run = combined_run(family)
        seeds, assembly = sweep_seeds(run), run["assembly"]
        assert len(seeds) == len(set(seeds)) == 5, family
        code, lines, err = dry_run(run, tmp_path, capsys, "--steps", SWEEP_STEP, "--reference", str(tmp_path / "reference.txt"))
        assert (code, err) == (0, "") and len(lines) == 5, family
        for seed, line in zip(seeds, lines, strict=True):
            words = line.split(" > ")[0].split()
            written = inputs / "utrannotator_sweep" / f"sweep-seed{seed}.txt"
            place = hashlib.sha256(os.fsencode(written)).hexdigest()[:12]
            assert words == [
                "docker", "run", "--rm", "--name", f"seeded-{seed}-{assembly}-{place}", "--network", "none",
                "--user", f"{os.getuid()}:{os.getgid()}",
                "-e", "HOME=/tmp", "-e", f"PERL_HASH_SEED={seed}", "-e", "PERL_PERTURB_KEYS=0",
                "-v", f"{inputs / 'ensembl_vep_cache'}:/work/cache:ro",
                "-v", f"{inputs / 'input'}:/work/in/input.vcf:ro",
                "-v", f"{inputs / 'utrannotator_sweep'}:/work/out",
                "-v", f"{inputs}:/work/ref:ro",
                "-v", f"{inputs / 'plugin_data' / assembly.lower() / 'utrannotator'}:/work/plugin_data:ro",
                "-w", "/work/out", driver.pinned_image(),
                "vep", "-i", "/work/in/input.vcf", "-o", f"/work/out/sweep-seed{seed}.txt",
                "--offline", "--cache", "--dir_cache", "/work/cache", "--cache_version", "116",
                "--species", "homo_sapiens", "--assembly", assembly, "--format", "vcf", "--buffer_size", "5000",
                "--fork", "1", "--fasta", "/work/ref/fasta", "--no_stats", "--quiet", "--force_overwrite", "--hgvs",
                "--plugin", f"UTRAnnotator,file=/work/plugin_data/uORF_5UTR_{assembly}_PUBLIC.txt",
                "--dir_plugins", "/plugins",
            ], line
            assert line.endswith(f" > {tmp_path / 'work' / family / f'sweep-seed{seed}.stdout'}"
                                 f" 2> {tmp_path / 'work' / family / f'sweep-seed{seed}.stderr'}")


def test_the_seeds_are_read_from_the_score_words_and_stand_nowhere_in_the_driver(tmp_path: Path, capsys) -> None:
    source = Path(driver.__file__).read_text(encoding="utf-8")
    _runs, inputs = driver.load_runs()
    for family in SWEEP_RUNS:
        for seed in sweep_seeds(combined_run(family)):
            assert f"sweep-seed{seed}" not in source and f"PERL_HASH_SEED={seed}" not in source
    made_up = dict(combined_run("plugins_grch38"))
    made_up["score"] = ["--plugins", "UTRAnnotator", "--utrannotator-sweep", "{utrannotator_sweep}/sweep-seed0042.txt",
                        "--input", "{input}", "--utrannotator-sweep", "{utrannotator_sweep}/sweep-seed7.txt"]
    args = driver.parse(sweep_argv(tmp_path, "plugins_grch38"))
    commands = driver.Plan(made_up, args, inputs).sweep_commands()
    assert [command["argv"][command["argv"].index("-o") + 1] for command in commands] == [
        "/work/out/sweep-seed0042.txt", "/work/out/sweep-seed7.txt"]
    assert ["PERL_HASH_SEED=0042" in command["argv"] for command in commands] == [True, False]
    assert [command["writes"] for command in commands] == [
        tmp_path / "inputs/utrannotator_sweep" / name for name in ("sweep-seed0042.txt", "sweep-seed7.txt")]


def test_a_run_whose_score_words_name_no_seeded_output_gets_no_sweep_command_and_one_line(tmp_path: Path, capsys) -> None:
    runs, _inputs = driver.load_runs()
    corpus = next(run for run in runs if run.get("corpus") == CORPUS and run["set"] == driver.BASE_SET)
    for run in (next(run for run in runs if run.get("family") == "gwas"), corpus):
        code, lines, err = dry_run(run, tmp_path, capsys, "--steps", SWEEP_STEP, "--reference", str(tmp_path / "reference.txt"))
        assert (code, lines) == (0, []), driver.run_key(run)
        assert err == f"{driver.run_key(run)}: the run's score words name no seeded output; the utrannotator-sweep step writes nothing\n"


def test_the_sweep_step_writes_nothing_in_a_dry_run(tmp_path: Path, capsys, monkeypatch) -> None:
    log = stand_in_docker(tmp_path, monkeypatch)
    assert driver.main([*sweep_argv(tmp_path, "plugins_grch37"), "--dry-run"]) == 0
    assert len(capsys.readouterr().out.splitlines()) == 5
    assert log.read_text(encoding="utf-8") == ""
    assert not (tmp_path / "work").exists() and not (tmp_path / "inputs").exists()


def test_the_sweep_step_writes_the_five_outputs_the_score_words_name(tmp_path: Path, capsys, monkeypatch) -> None:
    log = stand_in_docker(tmp_path, monkeypatch)
    seeds = sweep_seeds(combined_run("plugins_grch38"))
    assert driver.main(sweep_argv(tmp_path, "plugins_grch38")) == 0
    captured = capsys.readouterr()
    assert (captured.out, captured.err) == ("", "")
    sweep = tmp_path / "inputs/utrannotator_sweep"
    assert sorted(path.name for path in sweep.iterdir()) == sorted(f"sweep-seed{seed}.txt" for seed in seeds)
    assert [(sweep / f"sweep-seed{seed}.txt").read_text(encoding="utf-8") for seed in seeds] == [f"{seed}\n" for seed in seeds]
    calls = log.read_text(encoding="utf-8").splitlines()
    assert driver.main([*sweep_argv(tmp_path, "plugins_grch38"), "--dry-run"]) == 0
    printed = [line.split(" > ")[0].removeprefix("docker ") for line in capsys.readouterr().out.splitlines()]
    assert calls == printed and len(calls) == 5


def test_a_sweep_output_that_stands_is_not_written_over_without_the_callers_word(tmp_path: Path, capsys, monkeypatch) -> None:
    log = stand_in_docker(tmp_path, monkeypatch)
    seeds = sweep_seeds(combined_run("plugins_grch37"))
    held = tmp_path / "inputs/utrannotator_sweep" / f"sweep-seed{seeds[2]}.txt"
    held.parent.mkdir(parents=True)
    held.write_text("held\n", encoding="utf-8")
    line = refusal(capsys, *sweep_argv(tmp_path, "plugins_grch37"))
    assert line.endswith(f"the utrannotator-sweep step of plugins_grch37 would write over {held}: pass --overwrite-sweep")
    assert log.read_text(encoding="utf-8") == "" and held.read_text(encoding="utf-8") == "held\n"
    assert [path.name for path in held.parent.iterdir()] == [held.name]
    assert driver.main([*sweep_argv(tmp_path, "plugins_grch37"), "--overwrite-sweep"]) == 0
    assert held.read_text(encoding="utf-8") == f"{seeds[2]}\n" and len(log.read_text(encoding="utf-8").splitlines()) == 5


def test_a_seeded_run_that_ends_non_zero_is_refused_with_its_error_file(tmp_path: Path, capsys, monkeypatch) -> None:
    log = stand_in_docker(tmp_path, monkeypatch)
    monkeypatch.setenv("STAND_IN_EXIT", "3")
    seeds = sweep_seeds(combined_run("plugins_grch37"))
    line = refusal(capsys, *sweep_argv(tmp_path, "plugins_grch37"))
    error_file = tmp_path / "work/plugins_grch37" / f"sweep-seed{seeds[0]}.stderr"
    assert line.endswith(f"the utrannotator-sweep step of plugins_grch37 ended 3; its standard error is in {error_file}")
    assert error_file.read_text(encoding="utf-8") == "no\n"
    assert len(log.read_text(encoding="utf-8").splitlines()) == 1
    assert list((tmp_path / "inputs/utrannotator_sweep").iterdir()) == []


def test_a_seeded_container_is_named_for_its_seed_its_assembly_and_the_file_it_writes() -> None:
    name = driver.seeded_container("11", "GRCh38", Path("/data/sweep/sweep-seed11.txt"))
    assert name == "seeded-11-GRCh38-6bd82590f93f"
    assert driver.seeded_container("11", "GRCh38", Path("/data/sweep/sweep-seed11.txt")) == name
    assert driver.seeded_container("11", "GRCh38", Path("/data/other/sweep-seed11.txt")) != name


def test_two_sweeps_run_side_by_side_into_two_directories_and_refuse_each_other_on_one(
    tmp_path: Path, capsys, monkeypatch
) -> None:
    log = stand_in_docker(tmp_path, monkeypatch)
    running = tmp_path / "running"
    running.mkdir()
    monkeypatch.setenv("STAND_IN_RUNNING", str(running))
    seeds = sweep_seeds(combined_run("plugins_grch38"))

    def sweep_into(place: Path, work: Path, *extra: str) -> list[str]:
        return sweep_argv(tmp_path, "plugins_grch38", "--utrannotator-sweep", str(place), "--work-dir", str(work), *extra)

    places = [tmp_path / "one/sweep", tmp_path / "two/sweep"]
    for place in places:
        assert driver.main(sweep_into(place, place.parent / "work")) == 0, capsys.readouterr().err
        assert sorted(path.name for path in place.iterdir()) == sorted(f"sweep-seed{seed}.txt" for seed in seeds)
    names = [path.name for path in running.iterdir()]
    assert len(names) == 10 and all(DOCKER_NAME_RE.fullmatch(name) for name in names)
    # A third sweep into the first directory, whose five containers still run.
    line = refusal(capsys, *sweep_into(places[0], tmp_path / "three/work", "--overwrite-sweep"))
    error_file = tmp_path / "three/work/plugins_grch38" / f"sweep-seed{seeds[0]}.stderr"
    assert line.endswith(f"the utrannotator-sweep step of plugins_grch38 ended 125; its standard error is in {error_file}")
    assert len(log.read_text(encoding="utf-8").splitlines()) == 11


def test_the_sweep_reads_the_callers_records_and_the_runs_input_where_none_is_named(tmp_path: Path, capsys) -> None:
    run = combined_run("plugins_grch37")
    reference = ("--steps", SWEEP_STEP, "--reference", str(tmp_path / "reference.txt"))
    _code, lines, _err = dry_run(run, tmp_path, capsys, *reference)
    assert all(f" -v {tmp_path / 'inputs/input'}:/work/in/input.vcf:ro " in line for line in lines) and len(lines) == 5
    records = tmp_path / "records" / "residual-records.vcf.gz"
    _code, lines, _err = dry_run(run, tmp_path, capsys, *reference, "--sweep-input", str(records))
    assert len(lines) == 5
    for line in lines:
        assert f" -v {records}:/work/in/input.vcf.gz:ro " in line and " vep -i /work/in/input.vcf.gz -o /work/out/sweep-seed" in line
        assert str(tmp_path / "inputs/input") not in line
    usage = " ".join(driver.__doc__.split())
    assert "the whole input, or the records behind the UTRAnnotator residual rows of a first score" in usage


def test_the_check_step_names_what_the_sweep_reads_and_not_the_outputs_it_writes(tmp_path: Path) -> None:
    _runs, inputs = driver.load_runs()
    run = combined_run("plugins_grch38")
    records = tmp_path / "records.vcf"
    args = driver.parse([*sweep_argv(tmp_path, "plugins_grch38"), "--sweep-input", str(records)])
    plan = driver.Plan(run, args, inputs)
    sweep = [str(tmp_path / "inputs/utrannotator_sweep" / f"sweep-seed{seed}.txt") for seed in sweep_seeds(run)]
    missing = plan.missing(("check", SWEEP_STEP, "score"))
    assert not set(sweep) & set(missing)
    for path in (records, tmp_path / "inputs/fasta", tmp_path / "inputs/fasta.fai",
                 tmp_path / "inputs/ensembl_vep_cache/homo_sapiens/116_GRCh38",
                 tmp_path / "inputs/plugin_data/grch38/utrannotator/uORF_5UTR_GRCh38_PUBLIC.txt"):
        assert str(path) in missing, path
    assert set(sweep) <= set(plan.missing(("check", "score")))
    assert str(records) not in plan.missing(("check", "score"))


def test_both_reads_of_the_cache_directory_follow_the_one_species_and_cache_version(tmp_path: Path, monkeypatch) -> None:
    runs, inputs = driver.load_runs()
    for run in runs:
        words = run["ensembl_vep"]
        named = (words[words.index("--species") + 1], words[words.index("--cache_version") + 1])
        assert named == (driver.SPECIES, driver.CACHE_VERSION), driver.run_key(run)
    args = driver.parse(["--run", "plugins_grch38", "--work-dir", str(tmp_path / "work"), "--run-ensembl-vep",
                         *named_inputs(tmp_path)])
    plan = driver.Plan(combined_run("plugins_grch38"), args, inputs)
    cache = tmp_path / "inputs/ensembl_vep_cache"

    def reads() -> tuple[Path, list[str]]:
        """The directory a seeded run reads, and what the check of the `ensembl-vep` step looks for under the cache."""
        looked_for = [path for path in plan.missing(("check", "ensembl-vep")) if path.startswith(f"{cache}{os.sep}")]
        return plan.sweep_commands()[0]["reads"][0], looked_for

    stated = cache / "homo_sapiens" / "116_GRCh38"
    assert reads() == (stated, [str(stated)])
    monkeypatch.setattr(driver, "SPECIES", "arabidopsis_thaliana")
    monkeypatch.setattr(driver, "CACHE_VERSION", "63")
    moved = cache / "arabidopsis_thaliana" / "63_GRCh38"
    assert reads() == (moved, [str(moved)])
    # The species directory and the version prefix stand in the two definitions and nowhere else.
    source = Path(driver.__file__).read_text(encoding="utf-8").splitlines()
    assert [line for line in source if "homo_sapiens" in line] == ['SPECIES = "homo_sapiens"']
    assert [line for line in source if "CACHE_VERSION =" in line] == ['CACHE_VERSION = "116"']
    assert not [line for line in source if "116_" in line]


def inputs_without(tmp_path: Path, *left_out: str) -> list[str]:
    """The named inputs a caller gives, but for the options left out."""
    words = named_inputs(tmp_path)
    for option in left_out:
        at = words.index(option)
        del words[at : at + 2]
    return words


def test_a_sweep_only_call_with_its_own_records_is_asked_for_neither_input_nor_reference(tmp_path: Path, capsys) -> None:
    records = tmp_path / "records" / "residual-records.vcf.gz"
    for family in SWEEP_RUNS:
        argv = ["--run", family, "--work-dir", str(tmp_path / "work"), "--dry-run", "--steps", SWEEP_STEP,
                "--sweep-input", str(records), *inputs_without(tmp_path, "--input")]
        assert driver.main(argv) == 0, family
        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        assert captured.err == "" and len(lines) == 5, family
        assert all(f" -v {records}:/work/in/input.vcf.gz:ro " in line for line in lines), family
        assert not any(str(tmp_path / "inputs/input") in line for line in lines), family
    assert not (tmp_path / "work").exists()


def test_a_sweep_only_call_with_no_records_is_refused_in_one_line(tmp_path: Path, capsys) -> None:
    argv = ["--run", "plugins_grch38", "--work-dir", str(tmp_path / "work"), "--dry-run", "--steps", SWEEP_STEP,
            *inputs_without(tmp_path, "--input")]
    line = refusal(capsys, *argv)
    assert line.endswith(f"the {SWEEP_STEP} step of plugins_grch38 reads a file of records: pass --sweep-input or --input")
    # With the run's input and no reference word the step reads the input, as the usage text says.
    assert driver.main([*argv, "--input", str(tmp_path / "inputs/input")]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 5 and all(f" -v {tmp_path / 'inputs/input'}:/work/in/input.vcf:ro " in line for line in lines)


def test_a_call_is_asked_for_the_input_and_the_reference_only_where_a_named_step_reads_them(tmp_path: Path, capsys) -> None:
    work = ["--work-dir", str(tmp_path / "work"), "--dry-run"]
    given = inputs_without(tmp_path, "--input")
    with_input = ["--input", str(tmp_path / "inputs/input")]
    # vep-rs reads the input and no reference; the score reads both; a whole run reads both.
    assert driver.main(["--run", "gwas", "--steps", "check,vep-rs", *work, *given, *with_input]) == 0
    assert len(capsys.readouterr().out.splitlines()) == 1
    assert "pass --input" in refusal(capsys, "--run", "gwas", "--steps", "vep-rs", *work, *given)
    assert "pass --input" in refusal(capsys, "--run", "gwas", "--steps", "score", "--reference", "r.txt", *work, *given)
    assert "pass --input" in refusal(capsys, "--run", "gwas", "--steps", "ensembl-vep", "--run-ensembl-vep", *work, *given)
    line = refusal(capsys, "--run", "gwas", "--steps", "score", *work, *given, *with_input)
    assert "pass --reference or --run-ensembl-vep" in line
    assert "pass --input" in refusal(capsys, "--run", "gwas", *work, *given)
    assert "pass --reference or --run-ensembl-vep" in refusal(capsys, "--run", "gwas", *work, *given, *with_input)
    # The check step alone reads names of files the other named steps read: none here.
    assert driver.main(["--run", "gwas", "--steps", "check", *work, *given]) == 0
    assert capsys.readouterr().out == ""
    usage = " ".join(driver.__doc__.split())
    assert "asked for `--input` only where a named step reads it" in usage
    assert "for `--reference` or `--run-ensembl-vep` only where `score` is named" in usage
