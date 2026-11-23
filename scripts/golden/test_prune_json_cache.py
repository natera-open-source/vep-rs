"""Tests for the summary `prune_json_cache.py` prints, on a synthetic cache, and for the figures
of that summary a committed corpus states of its own slice."""
from __future__ import annotations

import importlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
pjc = importlib.import_module("prune_json_cache")
GOLDEN = HERE.parent.parent / "tests" / "golden"


def test_prune_summary_counts_the_records_written_and_each_transcript_once(tmp_path: Path, capsys):
    """`transcripts_per_chromosome` counts the transcript records written, one for each shard a
    transcript is stored in, and `distinct_transcripts_per_chromosome` counts a transcript once
    however many shards hold it, told apart as the runtime's loader tells them: by source cache
    and `dbID` (a number or its digits), or by stable id and version where a record has no `dbID`."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    crossing = {"stable_id": "ENST_B", "version": "2", "start": "990000", "end": "1010000", "dbID": "2", "source_cache": "Ensembl"}
    crossing_as_number = {**crossing, "dbID": 2}
    # A merged cache numbers its two sources apart, so the same dbID under RefSeq is another transcript.
    other_source = {"stable_id": "NM_B", "start": "995000", "end": "996000", "dbID": "2", "source_cache": "RefSeq"}
    unnumbered = {"stable_id": "ENST_N", "version": 1, "start": "999000", "end": "1001000"}
    next_version = {**unnumbered, "version": 2}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([crossing, other_source, unnumbered, next_version]))
    (cache / "transcripts" / "21" / "1000001-2000000.json").write_text(json.dumps([crossing_as_number, unnumbered, next_version]))
    (cache / "transcripts" / "22").mkdir()
    (cache / "transcripts" / "22" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_D", "start": "5", "end": "9", "dbID": "2"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t995500\t.\tA\tG\t.\t.\t.\n22\t7\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    # Seven records on 21 in two shards are four transcripts: the one with a dbID, the RefSeq model
    # of the same dbID, and the two versions of the one without.
    assert capsys.readouterr().out == (
        '{"transcripts_per_chromosome": {"21": 7, "22": 1}, "distinct_transcripts_per_chromosome": {"21": 4, "22": 1}, '
        '"shards": 3, "hgnc_donors": 0, "manifest_files": 4, "spans_by_chromosome": {"21": 1, "22": 1}}\n')
    assert len(json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())) == 4
    assert len(json.loads((out / "transcripts" / "21" / "1000001-2000000.json").read_text())) == 3
    # The runtime reads a record's number from a JSON number or a string, takes null or an empty
    # string for none, a fraction for its whole part and a number below zero for zero. Ten records
    # in two shards are five transcripts: two whose dbID is an empty string, the one whose dbID is
    # spelled three ways, the one whose version is empty or absent, and the two numbered below
    # zero, which the runtime loads as one.
    spelled = tmp_path / "spelled"
    (spelled / "transcripts" / "21").mkdir(parents=True)
    span = {"start": "999000", "end": "1001000"}
    no_id = {"stable_id": "ENST_E", "dbID": "", **span}
    padded = {"stable_id": "ENST_P", "dbID": "07", **span}
    no_version = {"stable_id": "ENST_V", "dbID": None, "version": "", **span}
    below_zero = {"stable_id": "ENST_X", "dbID": -1, **span}
    (spelled / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps(
        [no_id, {**no_id, "stable_id": "ENST_F"}, padded, no_version, below_zero]))
    (spelled / "transcripts" / "21" / "1000001-2000000.json").write_text(json.dumps(
        [no_id, {**padded, "dbID": 7}, {**padded, "dbID": 7.9}, {"stable_id": "ENST_V", **span}, {**below_zero, "stable_id": "ENST_Y", "dbID": -2}]))
    assert pjc.prune(spelled, tmp_path / "spelled_out", {"21": [(1000000, 1000000)]}, 5000) == {
        "transcripts_per_chromosome": {"21": 10}, "distinct_transcripts_per_chromosome": {"21": 5}, "shards": 2, "hgnc_donors": 0}


def test_the_tool_and_its_readme_say_what_each_count_holds():
    docstring = " ".join(pjc.__doc__.split())
    assert "counts, per chromosome, the transcript records written (`transcripts_per_chromosome`:" in docstring
    assert "the transcripts kept and the spans searched" not in docstring
    readme = " ".join((HERE.parent / "README.md").read_text(encoding="utf-8").split())
    assert ("Its `transcripts_per_chromosome` counts the transcript records written, one for each shard a transcript occupies, "
            "and its `distinct_transcripts_per_chromosome` counts each transcript once") in readme


def test_a_corpus_record_states_the_figures_the_tool_prints_for_its_slice(tmp_path: Path):
    """A `provenance.json` that states `transcripts_per_chromosome` under `pruned_cache` states
    what the tool prints under that key when the corpus's committed slice is cut again from
    itself, the transcript records, and a `distinct_transcripts_per_chromosome` beside it states
    the tool's figure of that name. The slice is `json_cache/`, or the one leaf of an Ensembl
    layout under `cache_root/`."""
    differing = []
    held = 0
    for record in sorted(GOLDEN.glob("*/*/provenance.json")):
        stated = json.loads(record.read_text(encoding="utf-8")).get("pruned_cache")
        if not isinstance(stated, dict) or "transcripts_per_chromosome" not in stated:
            continue
        corpus = record.parent
        name = str(corpus.relative_to(GOLDEN))
        slices = [found.parent for found in (*corpus.glob("json_cache/transcripts"), *corpus.glob("cache_root/*/*/transcripts"))]
        assert len(slices) == 1, f"{name}: one committed slice expected, found {slices}"
        # One span over every sequence of the slice keeps every record it holds.
        whole = {sequence.name: [(-(1 << 62), 1 << 62)] for sequence in (slices[0] / "transcripts").iterdir() if sequence.is_dir()}
        cut_again = tmp_path / "cut_again"
        printed = pjc.prune(slices[0], cut_again, whole, 0)
        shutil.rmtree(cut_again, ignore_errors=True)
        for key in ("transcripts_per_chromosome", "distinct_transcripts_per_chromosome"):
            if key in stated and stated[key] != printed[key]:
                differing.append(f"{name}: {key} states {stated[key]}, the tool prints {printed[key]}")
        held += 1
    assert not differing, "\n".join(differing)
    assert held > 0, f"no corpus under {GOLDEN} states the figure"


def test_the_cat_and_salmon_records_state_both_figures_of_their_slices():
    """A transcript stored in two region files is two transcript records, so the cat's slice
    holds 52 more transcript records than transcripts on `F1` and the salmon's 41 more on `8`:
    each corpus record states both figures."""
    expected = {
        "F.catus_Fca126_mat1.0-cat": {
            "transcripts_per_chromosome": {"F1": 682, "JAFEKA010000047.1": 18},
            "distinct_transcripts_per_chromosome": {"F1": 630, "JAFEKA010000047.1": 18},
        },
        "Ssal_v3.1-salmon": {
            "transcripts_per_chromosome": {"8": 778, "CAJNNT020003670.1": 73},
            "distinct_transcripts_per_chromosome": {"8": 737, "CAJNNT020003670.1": 73},
        },
    }
    for corpus, figures in expected.items():
        stated = json.loads((GOLDEN / "116" / corpus / "provenance.json").read_text(encoding="utf-8"))["pruned_cache"]
        assert {key: stated.get(key) for key in figures} == figures, corpus
