"""Tests for reevaluate.py: the schema check, the readers of each compare kind, the verdict, the page check, and a run of the two fixture reproductions through fake engines."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import stat
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_divergence_order as cdo
import reevaluate as rv

TESTDATA = HERE / "testdata"
FIXTURES = TESTDATA / "reproductions"
IMAGE = "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da"
UNOBSERVED_BASIS = "well-formed under VCF 4.2; not seen in the measured datasets"
SIGHTING_BASIS = "seen in a structural-variant file of a population dataset: 4,158 pairs"


def repro(name="x", **overrides) -> rv.Reproduction:
    data = {
        "entry": name,
        "page": "docs/intended-divergences.md",
        "verified_against": "116.2",
        "image": IMAGE,
        "input_class": "unobserved",
        "input_class_basis": UNOBSERVED_BASIS,
        "input": {"text": "21\t100\trec\tG\tA\t.\t.\t.\n"},
        "cache": {"assembly": "GRCh38", "json_cache": "tests/golden/116/GRCh38/json_cache"},
        "run": {"flags": ["--hgvs"]},
        "compare": {"kind": "field", "field": "HGVSp", "where": {"record": "rec"}},
        "expected": {"ensembl_vep": ["a"], "vep_rs": ["b"]},
    }
    for key, value in overrides.items():
        if value is None:
            data.pop(key, None)
        elif isinstance(value, dict) and isinstance(data.get(key), dict):
            data[key] = {**data[key], **value}
            for k, v in value.items():
                if v is None:
                    data[key].pop(k)
        else:
            data[key] = value
    return rv.Reproduction(name, Path(f"/nowhere/{name}.toml"), data)


def engine_run(tmp_path: Path, text: str, name: str = "default.txt", exit_code: int = 0, stderr: str = "", stdout: str = "", warnings: str | None = None, stats: str | None = None) -> rv.EngineRun:
    out = tmp_path / name
    out.write_text(text)
    if warnings is not None:
        out.with_name(out.name + "_warnings.txt").write_text(warnings)
    if stats is not None:
        out.with_name(out.name + "_summary.txt").write_text(stats)
    return rv.EngineRun("vep-rs", ["vep"], exit_code, out, stdout, stderr)


TABLE = (
    "## ENSEMBL VARIANT EFFECT PREDICTOR v116.2\n"
    "## Output produced at 2020-01-01\n"
    "## Extra column keys:\n"
    "## HGVSp : HGVS protein sequence name\n"
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    "rec\t21:100\tA\tENSG1\tENST1\tTranscript\tmissense_variant,splice_region_variant\t10\t10\t4\tL/P\tctA/ccA\t-\tIMPACT=MODERATE;HGVSp=ENSP1:p.Leu4Pro;DISTANCE=0\n"
    "rec\t21:100\tA\tENSG1\tENST2\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
    "other\t21:300\tC\tENSG3\tENST3\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=4000\n"
)


def test_validate_accepts_the_committed_shape_and_names_each_violation():
    assert rv.validate(repro()) == []
    assert rv.validate(repro(entry="other")) == ["entry is 'other', the file is x.toml"]
    assert "image is not pinned by digest" in rv.validate(repro(image="ensemblorg/ensembl-vep:release_116.2"))
    assert "input needs exactly one of text and path" in rv.validate(repro(input={"text": "a", "path": "b"}))
    assert "input needs exactly one of text and path" in rv.validate(repro(input={"text": None}))
    assert "cache.assembly is missing" in rv.validate(repro(cache={"assembly": None}))
    assert "cache needs json_cache or gxf" in rv.validate(repro(cache={"json_cache": None}))
    assert rv.validate(repro(cache={"json_cache": None, "gxf": "{corpus}/genes.gtf.gz"})) == []
    assert "cache.source is not ensembl, merged or refseq" in rv.validate(repro(cache={"source": "other"}))
    assert "unknown key cache.prepare" in rv.validate(repro(cache={"prepare": ["x"]}))
    assert "unknown key colour" in rv.validate(repro(colour="x"))
    assert "compare.kind 'flow' is not one of" in rv.validate(repro(compare={"kind": "flow"}))[0]
    assert "compare.kind field needs compare.field" in rv.validate(repro(compare={"field": None}))
    assert "compare.kind header needs compare.pattern" in rv.validate(repro(compare={"kind": "header", "field": None}))
    assert "compare.kind json needs compare.path" in rv.validate(repro(compare={"kind": "json", "field": None}))
    assert "compare.kind stats needs compare.section" in rv.validate(repro(compare={"kind": "stats", "field": None}))
    assert "compare.where has keys colour that compare.kind field does not read" in rv.validate(repro(compare={"where": {"colour": "x"}}))
    assert "compare.where has keys id that compare.kind rows does not read" in rv.validate(repro(compare={"kind": "rows", "field": None, "where": {"id": "x"}}))
    assert "compare.where has keys gene that compare.kind vcf does not read" in rv.validate(repro(compare={"kind": "vcf", "where": {"gene": "x"}}))
    assert "compare.where has keys record that compare.kind exit does not read" in rv.validate(repro(compare={"kind": "exit", "field": None}))
    assert rv.validate(repro(compare={"kind": "json", "field": None, "path": "x", "where": {"input": "l", "id": "r"}})) == []
    assert "compare.where has keys feature that compare.kind vcf does not read" in rv.validate(repro(compare={"kind": "vcf", "field": None, "line": True, "where": {"id": "r", "feature": "ENST1"}}))
    assert "expected.ensembl_vep equals expected.vep_rs: the entry records no divergence" in rv.validate(repro(expected={"ensembl_vep": "a", "vep_rs": ["a"]}))
    assert "expected needs ensembl_vep and vep_rs" in rv.validate(repro(expected={"ensembl_vep": "a", "vep_rs": None}))
    assert "run.format 'csv' is not one of" in rv.validate(repro(run={"format": "csv"}))[0]
    assert "run.tool 'grep' is not one of" in rv.validate(repro(run={"tool": "grep"}))[0]
    assert "a tool run needs run.args" in rv.validate(repro(run={"tool": "filter"}))
    assert rv.validate(repro(run={"tool": "filter", "args": ["-i", "{input}"]}, cache={"assembly": None, "json_cache": None})) == []
    assert "prepare is not a list of strings" in rv.validate(repro(prepare=[1]))
    assert rv.validate(repro(prepare=5)) == ["prepare is not a list of strings"]
    assert rv.validate(repro(input={"format": ["vcf"]})) == ["input.format is not a string"]
    assert rv.validate(repro(run={"tool": ["filter"]}, compare={"kind": ["rows"]})) == ["run.tool is not a string", "compare.kind is not a string"]
    assert "compare.kind field needs compare.field" in rv.validate(repro(compare={"field": None, "line": True}))
    assert rv.validate(repro(compare={"kind": "vcf", "field": None, "line": True, "where": {"id": "r"}})) == []
    assert "image is missing" in rv.validate(repro(image=5))
    assert "run.flags is not a list of strings" in rv.validate(repro(run={"flags": "--hgvs"}))
    assert "run.flags is not a list of strings" in rv.validate(repro(run={"flags": ["--distance", 500]}))
    assert "cache.chromosomes is not a list of strings" in rv.validate(repro(cache={"chromosomes": [21]}))
    assert "prepare is not a list of strings" in rv.validate(repro(prepare="true"))
    assert "expected.ensembl_vep_exit is not an integer" in rv.validate(repro(expected={"ensembl_vep": "a", "vep_rs": "b", "ensembl_vep_exit": "2"}))
    assert "cache.json_cache is not a string" in rv.validate(repro(cache={"json_cache": 5}))
    assert rv.validate(repro(cache="x")) == ["cache is not a table"]
    assert rv.validate(repro(compare=3)) == ["compare is not a table"]
    assert "input.text is not a string" in rv.validate(repro(input={"text": 5}))
    assert rv.validate(repro(expected={"ensembl_vep": "a", "vep_rs": "b", "ensembl_vep_exit": 2})) == []
    assert "compare.pattern is not a regular expression" in rv.validate(repro(compare={"kind": "lines", "field": None, "where": None, "pattern": "("}))[0]
    assert "compare.where.pattern is not a regular expression" in rv.validate(repro(compare={"where": {"pattern": "[" }}))[0]
    assert "compare.mask is not a regular expression" in rv.validate(repro(compare={"mask": ["ok", "(?P<x"]}))[0]
    assert "compare.where is not a table" in rv.validate(repro(compare={"where": "rec"}))


def test_validate_holds_the_input_class_to_three_values_and_its_basis_to_the_class():
    """Every reproduction names one of the three input classes and one line of basis, in the form its class asks for."""
    missing = ["input_class is missing", "input_class_basis is missing"]
    assert rv.validate(repro(input_class=None)) == missing[:1]
    assert rv.validate(repro(input_class_basis=None)) == missing[1:]
    assert rv.validate(repro(input_class=None, input_class_basis=None)) == missing
    assert rv.validate(repro(input_class="", input_class_basis=" ")) == missing
    assert rv.validate(repro(input_class=3, input_class_basis=["a"])) == missing
    for fourth in ("synthetic", "Malformed", "tier 3"):
        assert rv.validate(repro(input_class=fourth)) == [f"input_class {fourth!r} is not one of observed, unobserved, malformed"]
    for parted in ("\n", "\r", "\u2028"):
        assert rv.validate(repro(input_class_basis=f"well-formed under VCF 4.2;{parted}not seen in the measured datasets")) == ["input_class_basis is not one line"]
    assert rv.validate(repro(input_class_basis=UNOBSERVED_BASIS + "\n")) == ["input_class_basis is not one line"]
    unobserved = "input_class_basis of an unobserved entry does not read `well-formed under <format or option contract>; not seen in the measured datasets`"
    for wrong in ("allowed by VCF 4.2", "well-formed under VCF 4.2", "allowed by VCF 4.2; not seen in the measured datasets", "well-formed under ; not seen in the measured datasets", "well-formed under \t ; not seen in the measured datasets"):
        assert rv.validate(repro(input_class_basis=wrong)) == [unobserved], wrong
    assert rv.validate(repro(input_class_basis=UNOBSERVED_BASIS + "; classed by its most reachable shape, a bracket mate: a mate position that is no number breaks the VCF 4.2 specification, section 5.4")) == []
    observed = "input_class_basis of an observed entry is neither a tier 1 or tier 2 `Tier and count` cell nor `seen in <the real data its page names>: <the figure its page measures there>`"
    for right in ("1 · 7.3 per million · no flag needed", "2 · 156.8 per million · 0 flag sets, 28 corpora", SIGHTING_BASIS, "seen in a file (region 21:1-100) of a population dataset: 4,158 pairs: two of them on one record"):
        assert rv.validate(repro(input_class="observed", input_class_basis=right)) == [], right
    # The words `seen in`, the data name and the figure the page measures there, parted by the first colon a space follows.
    unframed = ("a population dataset: 4,158 pairs", "Seen in a population dataset: 4,158 pairs", "seen in a population dataset", "seen in : 4,158 pairs", "seen in : 4,158 pairs: two of them on one record", "seen in :4,158 pairs: two of them on one record", "seen in  a population dataset: 4,158 pairs", "seen in a population dataset:  4,158 pairs", "seen in a population dataset: ", "seen in a population dataset:4,158 pairs", "seen in a population dataset : 4,158 pairs", "seen in a population dataset : 4,158 pairs: two of them on one record")
    for wrong in ("3 · none on the cuts · 9 flag sets, 2 corpora", "4,158 pairs in a population dataset", *unframed):
        assert rv.validate(repro(input_class="observed", input_class_basis=wrong)) == [observed], wrong
    assert rv.SIGHTING.match("seen in a: b: c").groupdict() == {"where": "a", "what": "b: c"}
    assert rv.validate(repro(input_class="malformed", input_class_basis='POS is an Integer (VCF 4.2 specification, section 1.4.1 "Fixed fields"): the POS is text that is no number')) == []
    # A malformed entry's basis is the rule its input breaks, so the basis of another class left under it is refused.
    for other in (UNOBSERVED_BASIS, "well-formed under VCF 4.2", SIGHTING_BASIS, "3 · none on the cuts · 9 flag sets, 2 corpora", "2 · 156.8 per million · 0 flag sets, 28 corpora"):
        assert rv.validate(repro(input_class="malformed", input_class_basis=other)) == ["input_class_basis of a malformed entry reads as the basis of another class, not as the rule the input breaks"], other


def test_rows_field_and_count_read_the_table_by_where(tmp_path):
    run = engine_run(tmp_path, TABLE)
    rows = repro(compare={"kind": "rows", "field": None, "where": {"record": "rec"}})
    assert rv.extract(rows, run) == ["21:100 | A | ENST1 | missense_variant,splice_region_variant", "21:100 | A | ENST2 | intron_variant"]
    columns = repro(compare={"kind": "rows", "field": None, "where": {"record": "rec", "feature": "ENST1"}, "columns": ["Feature", "Codons", "DISTANCE", "SIFT"]})
    assert rv.extract(columns, run) == ["ENST1 | ctA/ccA | 0 | -"]
    field = repro(compare={"kind": "field", "field": "HGVSp", "where": {"record": "rec"}})
    assert rv.extract(field, run) == ["ENSP1:p.Leu4Pro", "-"]
    keyed = repro(compare={"kind": "field", "field": "HGVSp", "where": {"record": "rec"}, "with_key": True})
    assert rv.extract(keyed, run) == ["ENST1=ENSP1:p.Leu4Pro", "ENST2=-"]
    by_terms = repro(compare={"kind": "field", "field": "Protein_position", "where": {"consequence": "splice_region_variant,missense_variant"}})
    assert rv.extract(by_terms, run) == ["4"]
    alternatives = repro(compare={"kind": "count", "field": None, "where": {"feature": ["ENST1", "ENST3"]}})
    assert rv.extract(alternatives, run) == ["2"]
    pattern = repro(compare={"kind": "count", "field": None, "where": {"pattern": "DISTANCE=4000"}})
    assert rv.extract(pattern, run) == ["1"]
    nothing = repro(compare={"kind": "rows", "field": None, "where": {"record": "absent"}})
    assert rv.extract(nothing, run) == []


def test_header_lines_warnings_streams_exit_and_stats_read_their_files(tmp_path):
    run = engine_run(tmp_path, TABLE, exit_code=2, stderr="WARNING: x\nERROR: [input] bad at /host/path/in.vcf\n", stdout="done\n", warnings="WARNING: line 3 : skipped /host/in.vcf\nWARNING: line 4 : kept\n", stats="[General statistics]\nLines of input read\t2\n\n[Position in protein]\n00-10%\t1\n10-20%\t2\n\n[Other]\nx\t1\n")
    assert rv.extract(repro(compare={"kind": "header", "field": None, "pattern": "^## HGVSp"}), run) == ["## HGVSp : HGVS protein sequence name"]
    assert rv.extract(repro(compare={"kind": "lines", "field": None, "pattern": "ENST3"}), run) == [TABLE.splitlines()[7]]
    assert rv.extract(repro(compare={"kind": "warnings", "field": None, "pattern": "skipped"}), run) == ["WARNING: line 3 : skipped /host/in.vcf"]
    assert rv.extract(repro(compare={"kind": "stderr", "field": None, "pattern": "^ERROR"}), run) == ["ERROR: [input] bad at /host/path/in.vcf"]
    assert rv.extract(repro(compare={"kind": "stdout", "field": None, "pattern": "done"}), run) == ["done"]
    assert rv.extract(repro(compare={"kind": "exit", "field": None}), run) == ["2"]
    assert rv.extract(repro(compare={"kind": "stats", "field": None, "section": "[Position in protein]"}), run) == ["00-10%\t1", "10-20%\t2"]
    assert rv.extract(repro(compare={"kind": "stats", "field": None, "section": "[Missing]"}), run) == []
    raw = tmp_path / "raw.txt"
    raw.write_bytes(b"#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\nr\xe9c\t21:1\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\n")
    raw_run = rv.EngineRun("vep-rs", ["vep"], 0, raw, "", "")
    assert rv.extract(repro(compare={"kind": "field", "field": "Uploaded_variation", "where": None}), raw_run) == ["r\\xe9c"]
    zipped = tmp_path / "zipped.bgz"
    with gzip.open(zipped, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\n22\t5\n21\t7\n")
    assert rv.input_chromosomes(repro(), zipped) == ["22", "21"]
    no_files = rv.EngineRun("vep-rs", ["vep"], 1, tmp_path / "missing.txt", "", "")
    assert rv.extract(repro(compare={"kind": "warnings", "field": None, "pattern": "."}), no_files) == []
    assert rv.extract(repro(compare={"kind": "rows", "field": None, "where": None}), no_files) == []


def test_json_reads_a_path_over_the_matching_records(tmp_path):
    records = [
        {"id": "rec", "input": "21\t100\trec\tG\tA", "seq_region_name": "21", "start": 100, "most_severe_consequence": "missense_variant", "transcript_consequences": [{"transcript_id": "ENST1", "hgvsp": "ENSP1:p.Leu4Pro", "consequence_terms": ["missense_variant"]}, {"transcript_id": "ENST2", "consequence_terms": ["intron_variant"]}]},
        {"id": "other", "input": "21\t300\tother\tG\tC", "seq_region_name": "21", "start": 300, "most_severe_consequence": "downstream_gene_variant"},
    ]
    run = engine_run(tmp_path, "\n".join(json.dumps(r) for r in records) + "\n", name="json.jsonl")
    assert rv.extract(repro(compare={"kind": "json", "field": None, "path": "transcript_consequences[transcript_id=ENST1].hgvsp", "where": {"id": "rec"}}), run) == ["ENSP1:p.Leu4Pro"]
    assert rv.extract(repro(compare={"kind": "json", "field": None, "path": "transcript_consequences[].consequence_terms", "where": {"record": "rec"}}), run) == ['["missense_variant"]', '["intron_variant"]']
    assert rv.extract(repro(compare={"kind": "json", "field": None, "path": "most_severe_consequence", "where": None}), run) == ["missense_variant", "downstream_gene_variant"]
    assert rv.extract(repro(compare={"kind": "json", "field": None, "path": "most_severe_consequence", "where": {"location": "21:300"}}), run) == ["downstream_gene_variant"]
    assert rv.extract(repro(compare={"kind": "json", "field": None, "path": "start", "where": {"input": "21\t300\tother\tG\tC"}}), run) == ["300"]
    assert rv.extract(repro(compare={"kind": "json", "field": None, "path": "transcript_consequences[].hgvsp", "where": {"id": "other"}}), run) == []
    assert rv.json_path({"a": [{"k": "1", "v": 1}, {"k": "2", "v": 2}]}, "a[k=2].v") == [2]
    assert rv.json_path({"a": {"b": None}}, "a.b") == [None]
    assert rv.json_value(None) == "null"


def test_vcf_reads_csq_fields_by_record_feature_and_allele(tmp_path):
    vcf = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|HGVSp">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trec\tG\tA,C\t.\t.\tDP=4;CSQ=A|missense_variant|ENST1|ENSP1:p.Leu4Pro,A|intron_variant|ENST2|,C|missense_variant|ENST1|ENSP1:p.Leu4Arg\n"
        "21\t300\tother\tG\tC\t.\t.\t.\n"
    )
    run = engine_run(tmp_path, vcf, name="vcf.vcf")
    assert rv.extract(repro(compare={"kind": "vcf", "field": "HGVSp", "where": {"id": "rec", "feature": "ENST1"}}), run) == ["ENSP1:p.Leu4Pro", "ENSP1:p.Leu4Arg"]
    assert rv.extract(repro(compare={"kind": "vcf", "field": "HGVSp", "where": {"id": "rec", "feature": "ENST1", "allele": "C"}}), run) == ["ENSP1:p.Leu4Arg"]
    assert rv.extract(repro(compare={"kind": "vcf", "field": "Consequence", "where": {"location": "21:100", "allele": "A"}}), run) == ["missense_variant", "intron_variant"]
    assert rv.extract(repro(compare={"kind": "vcf", "field": "CSQ", "where": {"id": "other"}}), run) == ["-"]
    assert rv.extract(repro(compare={"kind": "vcf", "field": None, "line": True, "where": {"pos": "300"}}), run) == ["21\t300\tother\tG\tC\t.\t.\t."]
    assert rv.extract(repro(compare={"kind": "vcf", "field": "DP", "info": "DP", "where": {"id": "rec"}}), run) == ["4"]


def test_normalise_masks_paths_and_patterns_and_sorts_unless_told_otherwise():
    values = ["b at /opt/vep/src/x.pm line 3", "a pid 1234"]
    assert rv.normalise(values, {}) == ["a pid 1234", "b at <path> line 3"]
    assert rv.normalise(values, {"sorted": False, "mask_paths": False}) == values
    assert rv.normalise(values, {"mask": [r"\d+"]}) == ["a pid <masked>", "b at <path> line <masked>"]
    assert rv.normalise(["ENSP1:p.Leu4Pro"], {}) == ["ENSP1:p.Leu4Pro"]
    assert rv.normalise(["a description "], {}) == ["a description"]
    assert rv.normalise(["-/A", "C//G", "*/Q", "21:100 | -/A | ENST1 | x"], {}) == ["*/Q", "-/A", "21:100 | -/A | ENST1 | x", "C//G"]
    assert rv.normalise(["gzip: /work/in/x.bgz: not in gzip format"], {}) == ["gzip: <path> not in gzip format"]
    assert rv.normalise(["a description "], {"strip": False}) == ["a description "]
    assert rv.normalise(["21:100 | A | ENST1 | x"], {}) == ["21:100 | A | ENST1 | x"]


def test_verdict_is_closed_present_or_changed():
    expected = {"ensembl_vep": ["0.00"], "vep_rs": ["-"]}
    assert rv.verdict(expected, {}, ["0.00"], ["-"]) == "still-present"
    assert rv.verdict(expected, {}, ["-"], ["-"]) == "closed"
    assert rv.verdict(expected, {}, ["0.00"], ["0.34"]) == "changed"
    assert rv.verdict(expected, {}, ["0.01"], ["-"]) == "changed"
    assert rv.verdict({"ensembl_vep": "x", "vep_rs": "y"}, {}, ["x"], ["y"]) == "still-present"


def fake_engine(tmp_path: Path, name: str, outputs: Path) -> Path:
    """A stand-in engine: copies `<outputs>/<name>/<output file name>` (and its warnings file) to the `-o` path."""
    script = tmp_path / f"fake-{name}"
    script.write_text(
        "#!/usr/bin/env python3\n"
        "import shutil, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "out = Path(args[args.index('-o') + 1])\n"
        f"src = Path({str(outputs / name)!r}) / out.name\n"
        "shutil.copy(src, out)\n"
        "if src.with_name(src.name + '_warnings.txt').exists():\n"
        "    shutil.copy(src.with_name(src.name + '_warnings.txt'), out.with_name(out.name + '_warnings.txt'))\n"
        "Path(str(out) + '.args').write_text('\\n'.join(args))\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def perl_cache_root(tmp_path: Path) -> Path:
    root = tmp_path / "perl"
    (root / "homo_sapiens" / "116_GRCh38" / "21").mkdir(parents=True)
    return root


def test_the_two_fixture_reproductions_run_through_fake_engines(tmp_path, capsys):
    ensembl = fake_engine(tmp_path, "ensembl-vep", TESTDATA / "outputs")
    vep_rs = fake_engine(tmp_path, "vep-rs", TESTDATA / "outputs")
    root = perl_cache_root(tmp_path)
    results = tmp_path / "results.json"
    args = ["--ensembl-vep", str(ensembl), "--vep", str(vep_rs), "--perl-cache-dir", str(root), "--reproductions", str(FIXTURES), "--work-dir", str(tmp_path / "work"), "--json", str(results)]
    assert rv.main(args) == 0
    out = capsys.readouterr().out
    assert "still_present_fixture\tstill-present" in out
    assert "closed_fixture\tclosed\tEnsembl VEP ['ENSP2:p.Leu4%3D'], vep-rs ['ENSP2:p.Leu4%3D']" in out
    assert "| still-present | 1 |" in out and "| closed | 1 |" in out and "| total | 2 |" in out
    recorded = {r["entry"]: r for r in json.loads(results.read_text())}
    assert recorded["closed_fixture"]["verdict"] == "closed"
    assert recorded["still_present_fixture"]["vep_rs"] == ["21:100 | A | ENST1 | downstream_gene_variant | 5000"]
    assert recorded["still_present_fixture"]["ensembl_vep"] == []
    # The fixed arguments reached both engines, with vep-rs's cache and Ensembl VEP's --dir_cache root.
    ensembl_args = (tmp_path / "work" / "still_present_fixture" / "out" / "ensembl-vep" / "default.txt.args").read_text().split("\n")
    vep_rs_args = (tmp_path / "work" / "still_present_fixture" / "out" / "vep-rs" / "default.txt.args").read_text().split("\n")
    assert ensembl_args[ensembl_args.index("--dir_cache") + 1] == str(root)
    assert "--cache_version" in ensembl_args and "--fork" in ensembl_args and "--minimal" in ensembl_args and "--no_intergenic" in ensembl_args
    assert "--json_cache" in vep_rs_args and "--fork" not in vep_rs_args and "--no_stats" in vep_rs_args
    assert vep_rs_args[vep_rs_args.index("--format") + 1] == "vcf"
    written = (tmp_path / "work" / "still_present_fixture" / "stage" / "input" / "record.vcf").read_text()
    assert written.startswith("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t99\tfx1\t")
    assert rv.main([*args, "--expect-all-present"]) == 1
    assert rv.main([*args, "--entry", "still_present_fixture", "--expect-all-present"]) == 0
    assert rv.main([*args, "--entry", "unknown_fixture"]) == 1
    assert "no reproduction for unknown_fixture" in capsys.readouterr().err


def test_a_local_ensembl_vep_runs_its_sibling_tools(tmp_path):
    """`--ensembl-vep <dir>/vep` with a tool run calls `<dir>/filter_vep` or `<dir>/variant_recoder`; vep-rs runs `vep filter` and `vep recode`."""
    tools = tmp_path / "tools"
    tools.mkdir()
    for name in ("vep", "filter_vep", "variant_recoder"):
        script = tools / name
        script.write_text("#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\nargs = sys.argv[1:]\nPath(args[args.index('-o') + 1]).write_text('')\n")
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
    for tool, program in (("filter", "filter_vep"), ("recode", "variant_recoder")):
        r = repro(run={"tool": tool, "args": ["-i", "{input}", "-o", "{out}"]}, cache={"assembly": None, "json_cache": None})
        paths = rv.Paths(r, tmp_path / f"stage-{tool}", tmp_path / f"out-{tool}", None, container=False)
        (tmp_path / f"stage-{tool}").mkdir()
        options = argparse.Namespace(vep=tools / "vep", ensembl_vep=str(tools / "vep"), image=None, timeout=60)
        inp = rv.write_input(r, paths)
        ensembl = rv.run_engine("ensembl-vep", r, paths, options, None, inp, "default")
        assert ensembl.command[0] == str(tools / program) and ensembl.exit_code == 0
        vep_rs = rv.run_engine("vep-rs", r, paths, options, None, inp, "default")
        assert vep_rs.command[:2] == [str(tools / "vep"), tool]
    r = repro(run={"tool": "filter", "args": ["-i", "{input}", "-o", "{out}"]}, cache={"assembly": None, "json_cache": None})
    paths = rv.Paths(r, tmp_path / "stage-perl", tmp_path / "out-perl", None, container=False)
    (tmp_path / "stage-perl").mkdir()
    options = argparse.Namespace(vep=tools / "vep", ensembl_vep=f"/usr/bin/env {tools / 'vep'}", image=None, timeout=60)
    run = rv.run_engine("ensembl-vep", r, paths, options, None, rv.write_input(r, paths), "default")
    assert run.command[:2] == ["/usr/bin/env", str(tools / "filter_vep")]
    # A tool directory goes first on the local engine's PATH, alone when PATH is unset.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "samtools").write_text("")
    (tools / "filter_vep").write_text(f"#!{sys.executable}\nimport os, sys\nfrom pathlib import Path\nargs = sys.argv[1:]\nPath(args[args.index('-o') + 1]).write_text(os.environ['PATH'])\n")
    r = repro(run={"tool": "filter", "args": ["-i", "{input}", "-o", "{out}"], "tools": ["samtools"]}, cache={"assembly": None, "json_cache": None})
    options = argparse.Namespace(vep=tools / "vep", ensembl_vep=str(tools / "vep"), image=None, timeout=60, tool_dirs=[bin_dir])
    old_path = os.environ["PATH"]
    for environment_path, expected in ((old_path, f"{bin_dir}:{old_path}"), (None, str(bin_dir))):
        if environment_path is None:
            del os.environ["PATH"]
        try:
            paths = rv.Paths(r, tmp_path / "stage-path", tmp_path / "out-path", None, container=False)
            (tmp_path / "stage-path").mkdir(exist_ok=True)
            run = rv.run_engine("ensembl-vep", r, paths, options, None, rv.write_input(r, paths), "default")
        finally:
            os.environ["PATH"] = old_path
        assert run.exit_code == 0 and run.output.read_text() == expected


def test_a_failed_engine_run_is_an_error_unless_its_exit_is_recorded(tmp_path):
    """An engine that exits non-zero or writes no output is a failed run, never the recorded empty value; a recorded exit status makes that status the expectation."""
    ok = rv.EngineRun("ensembl-vep", ["vep"], 0, tmp_path / "default.txt", "", "")
    (tmp_path / "default.txt").write_text("#Uploaded_variation\n")
    failed = rv.EngineRun("ensembl-vep", ["vep"], 2, tmp_path / "missing.txt", "", "ERROR: no cache\n")
    rows = repro(compare={"kind": "rows", "field": None})
    assert rv.run_failure(rows, "ensembl-vep", ok) is None
    assert rv.run_failure(rows, "ensembl-vep", failed) == "ensembl-vep exited 2: ERROR: no cache"
    silent = rv.EngineRun("vep-rs", ["vep"], 0, tmp_path / "missing.txt", "", "")
    assert rv.run_failure(rows, "vep-rs", silent) == "vep-rs wrote no output file"
    recorded = repro(compare={"kind": "rows", "field": None}, expected={"ensembl_vep": [], "vep_rs": ["x"], "ensembl_vep_exit": 2})
    assert rv.run_failure(recorded, "ensembl-vep", failed) is None
    assert rv.run_failure(recorded, "ensembl-vep", ok) == "ensembl-vep exited 0, the entry records 2: "
    zero = repro(compare={"kind": "rows", "field": None}, expected={"ensembl_vep": [], "vep_rs": ["x"], "ensembl_vep_exit": 0})
    assert rv.run_failure(zero, "ensembl-vep", silent.__class__("ensembl-vep", ["vep"], 0, tmp_path / "missing.txt", "", "")) == "ensembl-vep wrote no output file"
    for kind in ("exit", "stderr", "stdout"):
        assert rv.run_failure(repro(compare={"kind": kind, "field": None, "pattern": "x"}), "ensembl-vep", failed) is None
    warnings = repro(compare={"kind": "warnings", "field": None, "pattern": "x"})
    assert rv.run_failure(warnings, "vep-rs", silent) is None
    assert rv.run_failure(warnings, "vep-rs", failed) == "ensembl-vep exited 2: ERROR: no cache".replace("ensembl-vep", "vep-rs")


def test_an_engine_that_dies_reads_as_an_error_in_a_run(tmp_path, capsys):
    """A stand-in Ensembl VEP that exits 2 without output turns the still-present fixture into an error, not a verdict."""
    dying = tmp_path / "dying"
    dying.write_text("#!/usr/bin/env python3\nimport sys\nsys.stderr.write('ERROR: [cache] no transcript shards\\n')\nsys.exit(2)\n")
    dying.chmod(dying.stat().st_mode | stat.S_IEXEC)
    vep_rs = fake_engine(tmp_path, "vep-rs", TESTDATA / "outputs")
    root = perl_cache_root(tmp_path)
    args = ["--ensembl-vep", str(dying), "--vep", str(vep_rs), "--perl-cache-dir", str(root), "--reproductions", str(FIXTURES), "--work-dir", str(tmp_path / "work"), "--entry", "still_present_fixture", "--expect-all-present"]
    assert rv.main(args) == 1
    out = capsys.readouterr().out
    assert "still_present_fixture\terror\tensembl-vep exited 2: ERROR: [cache] no transcript shards" in out
    assert "| error | 1 |" in out
    empty = tmp_path / "empty"
    empty.mkdir()
    assert rv.main(["--ensembl-vep", str(dying), "--vep", str(vep_rs), "--perl-cache-dir", str(root), "--reproductions", str(empty), "--expect-all-present"]) == 1
    assert "no reproduction under" in capsys.readouterr().err
    (empty / "bad.toml").write_text("entry = [\n")
    assert rv.main(["--list", "--reproductions", str(empty)]) == 1
    assert "ERROR: [reevaluate]" in capsys.readouterr().err
    (empty / "bad.toml").unlink()
    (empty / "latin.toml").write_bytes(b'entry = "caf\xe9"\n')
    assert rv.main(["--list", "--reproductions", str(empty)]) == 1
    assert "latin.toml: " in capsys.readouterr().err
    (empty / "latin.toml").unlink()
    text = (FIXTURES / "still_present_fixture.toml").read_text().replace('entry = "still_present_fixture"', 'entry = "number_fixture"').replace('json_cache = "{stage}/json_cache"', "json_cache = 5")
    (empty / "number_fixture.toml").write_text(text)
    assert rv.main(["--ensembl-vep", str(dying), "--vep", str(vep_rs), "--perl-cache-dir", str(root), "--reproductions", str(empty), "--expect-all-present"]) == 1
    assert "number_fixture\terror\tcache.json_cache is not a string" in capsys.readouterr().out


def test_a_cache_the_roots_lack_and_a_failing_prepare_step_are_errors(tmp_path, capsys):
    ensembl = fake_engine(tmp_path, "ensembl-vep", TESTDATA / "outputs")
    vep_rs = fake_engine(tmp_path, "vep-rs", TESTDATA / "outputs")
    empty_root = tmp_path / "empty"
    empty_root.mkdir()
    args = ["--ensembl-vep", str(ensembl), "--vep", str(vep_rs), "--perl-cache-dir", str(empty_root), "--reproductions", str(FIXTURES), "--work-dir", str(tmp_path / "work"), "--entry", "closed_fixture"]
    assert rv.main(args) == 0
    assert "closed_fixture\terror\tno Ensembl VEP cache homo_sapiens/116_GRCh38 with chromosomes 21" in capsys.readouterr().out
    broken = tmp_path / "broken"
    broken.mkdir()
    text = (FIXTURES / "closed_fixture.toml").read_text().replace('entry = "closed_fixture"', 'entry = "broken_fixture"\nprepare = ["false"]')
    (broken / "broken_fixture.toml").write_text(text)
    root = perl_cache_root(tmp_path)
    assert rv.main(["--ensembl-vep", str(ensembl), "--vep", str(vep_rs), "--perl-cache-dir", str(root), "--reproductions", str(broken), "--work-dir", str(tmp_path / "work2"), "--expect-all-present"]) == 1
    assert "broken_fixture\terror\tprepare step failed (1): false" in capsys.readouterr().out


def test_check_holds_the_page_and_the_reproductions_to_one_reference(tmp_path, capsys):
    header = cdo.ENTRY_TABLE_HEADER + "\n| --- | --- | --- | --- | --- |\n"

    def page(rows):
        body = "".join(f"| [`{k}`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | {v} | [P](x.md) |\n" for k, v in rows)
        return "# Index\n\n## Section\n\n" + header + body + "\n## How to read the counts\n\nprose\n"

    page_path = tmp_path / "page.md"
    page_path.write_text(page([("still_present_fixture", "116.2"), ("closed_fixture", "116.2")]))
    assert rv.check(page_path, "116.2", FIXTURES) == []
    assert rv.main(["--check", "--page", str(page_path), "--reproductions", str(FIXTURES)]) == 0
    assert "2 entries have a reproduction verified against 116.2" in capsys.readouterr().out
    page_path.write_text(page([("still_present_fixture", "116.2"), ("closed_fixture", "116.1"), ("orphan_entry", "116.2")]))
    errors = rv.check(page_path, "116.2", FIXTURES)
    assert errors == ["closed_fixture: the page says verified against '116.1', the reference is 116.2", "orphan_entry: no reproduction under " + str(FIXTURES)]
    page_path.write_text(page([("still_present_fixture", "116.2")]))
    assert rv.check(page_path, "116.2", FIXTURES) == ["closed_fixture: a reproduction with no row on the page"]
    page_path.write_text(page([("still_present_fixture", "116.2"), ("closed_fixture", "116.2")]))
    assert rv.check(page_path, "116.3", FIXTURES) == [
        "still_present_fixture: verified against 116.2, the reference is 116.3",
        "still_present_fixture: the page says verified against '116.2', the reference is 116.3",
        "closed_fixture: verified against 116.2, the reference is 116.3",
        "closed_fixture: the page says verified against '116.2', the reference is 116.3",
    ]
    assert rv.main(["--check", "--page", str(page_path), "--reproductions", str(FIXTURES), "--reference", "116.3"]) == 1
    assert "ERROR: [reevaluate]" in capsys.readouterr().err
    assert rv.main(["--check", "--page", str(tmp_path / "missing.md"), "--reproductions", str(FIXTURES)]) == 1


def write_record(directory: Path, name: str, input_class: str | None, basis: str | None) -> None:
    """A copy of the still-present fixture under another name with the given input class and basis; None leaves the field out."""
    fields = {"entry": name, "input_class": input_class, "input_class_basis": basis}
    lines = []
    for line in (FIXTURES / "still_present_fixture.toml").read_text().splitlines():
        key = next((k for k in fields if line.startswith(f"{k} = ")), None)
        if key is None:
            lines.append(line)
        elif fields[key]:
            lines.append(f'{key} = "{fields[key]}"')
    (directory / f"{name}.toml").write_text("\n".join(lines) + "\n")


def test_check_holds_each_input_class_to_the_pages_tier_and_section(tmp_path, capsys):
    """`--check` refuses a reproduction without its input class or its basis and one naming a fourth class, and holds each class to the page: an observed entry's basis is its cell at tier 1 or 2 and a `seen in` line at tier 3, a malformed entry stands in the malformed-input section and no other entry does."""
    header = cdo.ENTRY_TABLE_HEADER + "\n| --- | --- | --- | --- | --- |\n"
    seen, unseen = "2 · 98.6 per million · 15 flag sets, 20 corpora", "3 · none on the cuts · 1 flag set, 0 corpora"
    page_path, registry = tmp_path / "page.md", tmp_path / "registry"
    registry.mkdir()
    check_args = ["--check", "--page", str(page_path), "--reproductions", str(registry)]

    def table(rows):
        return header + "".join(f"| [`{key}`](x.md) | s | {cell} | 116.2 | [P](x.md) |\n" for key, cell in rows)

    def write_page(main, malformed):
        text = "# Index\n\n## Section\n\n" + table(main)
        if malformed:
            text += f"\n## {cdo.MALFORMED_SECTION}\n\n### Section\n\n" + table(malformed)
        page_path.write_text(text + "\n## How to read the counts\n\nprose\n")

    def line(key):
        return next(n for n, text in enumerate(page_path.read_text().splitlines(), 1) if f"`{key}`" in text)

    write_record(registry, "seen_entry", "observed", seen)
    write_record(registry, "unseen_entry", "unobserved", UNOBSERVED_BASIS)
    write_record(registry, "broken_entry", "malformed", "POS is an Integer (VCF 4.2 specification, section 1.4.1, POS): the POS is text that is no number")
    write_page([("seen_entry", seen), ("unseen_entry", unseen)], [("broken_entry", unseen)])
    assert rv.check(page_path, "116.2", registry) == []
    assert rv.main(check_args) == 0
    assert "3 entries have a reproduction verified against 116.2: 1 observed, 1 unobserved, 1 malformed" in capsys.readouterr().out

    # A cell that states no tier is not read as tier 3: the basis is held to it as written.
    write_page([("seen_entry", "tier 2"), ("unseen_entry", unseen)], [("broken_entry", unseen)])
    assert rv.check(page_path, "116.2", registry) == [f"seen_entry: input_class_basis is {seen!r}, the page's `Tier and count` cell is 'tier 2'"]
    write_page([("seen_entry", seen), ("unseen_entry", unseen)], [("broken_entry", unseen)])

    write_record(registry, "unseen_entry", None, None)
    assert rv.check(page_path, "116.2", registry) == ["unseen_entry: input_class is missing", "unseen_entry: input_class_basis is missing"]
    assert rv.main(check_args) == 1
    assert "ERROR: [reevaluate] unseen_entry: input_class is missing" in capsys.readouterr().err
    write_record(registry, "unseen_entry", "unobserved", None)
    assert rv.check(page_path, "116.2", registry) == ["unseen_entry: input_class_basis is missing"]
    write_record(registry, "unseen_entry", "synthetic", UNOBSERVED_BASIS)
    assert rv.check(page_path, "116.2", registry) == ["unseen_entry: input_class 'synthetic' is not one of observed, unobserved, malformed"]
    assert rv.main(check_args) == 1
    assert "is not one of observed, unobserved, malformed" in capsys.readouterr().err
    write_record(registry, "unseen_entry", "unobserved", UNOBSERVED_BASIS)

    stale = "2 · 98.6 per million · 14 flag sets, 20 corpora"
    write_record(registry, "seen_entry", "observed", stale)
    assert rv.check(page_path, "116.2", registry) == [f"seen_entry: input_class_basis is {stale!r}, the page's `Tier and count` cell is {seen!r}"]
    write_record(registry, "seen_entry", "unobserved", UNOBSERVED_BASIS)
    assert rv.check(page_path, "116.2", registry) == [f"line {line('seen_entry')}: seen_entry is tier 2 with input_class 'unobserved': a tier 1 or tier 2 entry is observed"]
    write_record(registry, "seen_entry", "synthetic", seen)
    assert rv.check(page_path, "116.2", registry) == ["seen_entry: input_class 'synthetic' is not one of observed, unobserved, malformed"]
    write_record(registry, "seen_entry", "observed", seen)

    # A tier 3 entry whose page records its shape in real data is observed, and its basis names that data in place of its cell.
    write_record(registry, "unseen_entry", "observed", SIGHTING_BASIS)
    assert rv.check(page_path, "116.2", registry) == []
    assert rv.main(check_args) == 0
    assert "3 entries have a reproduction verified against 116.2: 2 observed, 0 unobserved, 1 malformed" in capsys.readouterr().out
    write_record(registry, "unseen_entry", "observed", seen)
    assert rv.check(page_path, "116.2", registry) == [f"unseen_entry: input_class_basis is {seen!r}: an observed entry of tier 3 reads `seen in <the real data its page names>: <the figure its page measures there>`"]
    write_record(registry, "unseen_entry", "observed", None)
    assert rv.check(page_path, "116.2", registry) == ["unseen_entry: input_class_basis is missing"]
    write_record(registry, "unseen_entry", "observed", unseen)
    assert rv.check(page_path, "116.2", registry) == ["unseen_entry: input_class_basis of an observed entry is neither a tier 1 or tier 2 `Tier and count` cell nor `seen in <the real data its page names>: <the figure its page measures there>`"]
    write_record(registry, "seen_entry", "observed", SIGHTING_BASIS)
    write_record(registry, "unseen_entry", "observed", SIGHTING_BASIS)
    assert rv.check(page_path, "116.2", registry) == [f"seen_entry: input_class_basis is {SIGHTING_BASIS!r}, the page's `Tier and count` cell is {seen!r}"]
    write_record(registry, "seen_entry", "observed", seen)

    write_page([("seen_entry", seen)], [("broken_entry", unseen), ("unseen_entry", unseen)])
    assert rv.check(page_path, "116.2", registry) == [f"line {line('unseen_entry')}: unseen_entry stands under '{cdo.MALFORMED_SECTION}' with input_class 'observed'"]
    write_record(registry, "unseen_entry", "unobserved", UNOBSERVED_BASIS)
    assert rv.check(page_path, "116.2", registry) == [f"line {line('unseen_entry')}: unseen_entry stands under '{cdo.MALFORMED_SECTION}' with input_class 'unobserved'"]
    write_page([("seen_entry", seen), ("broken_entry", unseen), ("unseen_entry", unseen)], [])
    assert rv.check(page_path, "116.2", registry) == [f"line {line('broken_entry')}: broken_entry is malformed and stands under 'Section', outside '{cdo.MALFORMED_SECTION}'"]


def test_list_prints_every_reproduction(capsys):
    assert rv.main(["--list", "--reproductions", str(FIXTURES)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out == ["closed_fixture\t116.2\tdocs/intended-divergences.md", "still_present_fixture\t116.2\tdocs/intended-divergences.md"]


def test_running_needs_an_engine_and_a_binary(capsys, tmp_path):
    with pytest.raises(SystemExit):
        rv.main(["--vep", "/bin/true"])
    assert "--vep and --image" in capsys.readouterr().err
    assert rv.main(["--vep", "/bin/true", "--ensembl-vep", "/bin/true", "--extra-dir", "loftee", "--reproductions", str(FIXTURES)]) == 1
    assert "--extra-dir takes NAME=PATH" in capsys.readouterr().err
    assert rv.main(["--vep", "/bin/true", "--ensembl-vep", "/bin/true", "--extra-dir", f"Loftee={tmp_path}", "--reproductions", str(FIXTURES)]) == 1
    assert rv.main(["--vep", "/bin/true", "--ensembl-vep", "/bin/true", "--extra-dir", f"loftee={tmp_path / 'absent'}", "--reproductions", str(FIXTURES)]) == 1
    assert rv.main(["--vep", "/bin/true", "--ensembl-vep", "/bin/true", "--tool-dir", str(tmp_path / "absent"), "--reproductions", str(FIXTURES)]) == 1
    assert "--tool-dir" in capsys.readouterr().err


def test_write_input_keeps_a_header_the_text_brings_and_adds_samples(tmp_path):
    paths = rv.Paths(repro(input={"text": "##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n"}), tmp_path / "s1", tmp_path / "o1", None)
    assert rv.write_input(paths.repro, paths).read_text() == "##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n"
    paths = rv.Paths(repro(input={"text": "21\t1\tr\tA\tG\t.\t.\t.\tGT\t0/1\t1/1", "samples": ["S1", "S2"], "header": ["##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"]}), tmp_path / "s2", tmp_path / "o2", None)
    text = rv.write_input(paths.repro, paths).read_text()
    assert text.splitlines()[1] == "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    assert text.splitlines()[2] == rv.VCF_HEADER + "\tFORMAT\tS1\tS2"
    assert rv.input_chromosomes(paths.repro, paths.stage / "input" / "record.vcf") == ["21"]
    ensembl = repro(input={"text": "21 1 1 A/G +", "format": "ensembl"})
    paths = rv.Paths(ensembl, tmp_path / "s3", tmp_path / "o3", None)
    assert rv.write_input(ensembl, paths).name == "record.txt"
    assert rv.input_chromosomes(ensembl, paths.stage / "input" / "record.txt") == []
    assert rv.input_chromosomes(repro(cache={"chromosomes": ["MT"]}), tmp_path / "none") == ["MT"]


def test_find_perl_cache_takes_the_first_root_holding_the_chromosomes(tmp_path):
    first = tmp_path / "first"
    (first / "homo_sapiens" / "116_GRCh38" / "21").mkdir(parents=True)
    second = tmp_path / "second"
    (second / "homo_sapiens" / "116_GRCh38" / "21").mkdir(parents=True)
    (second / "homo_sapiens" / "116_GRCh38" / "MT").mkdir()
    (second / "homo_sapiens_merged" / "116_GRCh38" / "21").mkdir(parents=True)
    assert rv.find_perl_cache(repro(), [first, second], ["21"]) == first / "homo_sapiens" / "116_GRCh38"
    assert rv.find_perl_cache(repro(), [first, second], ["21", "MT"]) == second / "homo_sapiens" / "116_GRCh38"
    assert rv.find_perl_cache(repro(cache={"source": "merged"}), [first, second], ["21"]) == second / "homo_sapiens_merged" / "116_GRCh38"
    with pytest.raises(rv.ReproductionError, match="no Ensembl VEP cache homo_sapiens_refseq/116_GRCh38"):
        rv.find_perl_cache(repro(cache={"source": "refseq"}), [first, second], ["21"])
    # A reproduction that prepares its own cache may run with no root cache at all.
    staged_repro = repro(cache={"perl_cache": "{stage}/perl-cache", "source": "refseq"})
    assert rv.find_perl_cache(staged_repro, [first, second], ["21"]) is None
    paths = rv.Paths(staged_repro, tmp_path / "stage", tmp_path / "out", None)
    with pytest.raises(rv.ReproductionError, match="holds no homo_sapiens_refseq/116_GRCh38"):
        rv.prepared_perl_cache(staged_repro, paths)
    staged = tmp_path / "stage" / "perl-cache" / "homo_sapiens_refseq" / "116_GRCh38"
    staged.mkdir(parents=True)
    assert rv.prepared_perl_cache(staged_repro, paths) == staged
    assert rv.cache_relative(repro(cache={"species": "mus_musculus", "assembly": "GRCm39", "version": "116"})) == "mus_musculus/116_GRCm39"


def test_a_prepared_cache_is_read_after_the_prepare_steps(tmp_path, capsys):
    """`cache.perl_cache` names a root the prepare steps write, so it is resolved after them, with the roots' cache in PERL_CACHE; IMAGE is empty without `--image`."""
    ensembl = fake_engine(tmp_path, "ensembl-vep", TESTDATA / "outputs")
    vep_rs = fake_engine(tmp_path, "vep-rs", TESTDATA / "outputs")
    root = perl_cache_root(tmp_path)
    repros = tmp_path / "repros"
    repros.mkdir()
    text = (FIXTURES / "still_present_fixture.toml").read_text()
    step = 'mkdir -p "$STAGE/perl-cache/homo_sapiens" && cp -al "$PERL_CACHE" "$STAGE/perl-cache/homo_sapiens/116_GRCh38" && echo "$PERL_CACHE" > "$STAGE/root.txt" && printf "%s" "${IMAGE-unset}" > "$STAGE/image.txt"'
    text = text.replace('entry = "still_present_fixture"', f"entry = \"prepared_fixture\"\nprepare = ['{step}']")
    text = text.replace('json_cache = "{stage}/json_cache"', 'json_cache = "{stage}/json_cache"\nperl_cache = "{stage}/perl-cache"')
    (repros / "prepared_fixture.toml").write_text(text)
    work = tmp_path / "work"
    assert rv.main(["--ensembl-vep", str(ensembl), "--vep", str(vep_rs), "--perl-cache-dir", str(root), "--reproductions", str(repros), "--work-dir", str(work), "--expect-all-present"]) == 0
    assert "prepared_fixture\tstill-present" in capsys.readouterr().out
    assert (work / "prepared_fixture" / "stage" / "root.txt").read_text().strip() == str(root / "homo_sapiens" / "116_GRCh38")
    assert (work / "prepared_fixture" / "stage" / "image.txt").read_text() == ""
    args = (work / "prepared_fixture" / "out" / "ensembl-vep" / "default.txt.args").read_text().split("\n")
    assert args[args.index("--dir_cache") + 1] == str(work / "prepared_fixture" / "stage" / "perl-cache")


def test_substitute_fills_placeholders_per_engine_and_mounts_what_a_word_names(tmp_path):
    r = repro(cache={"json_cache": "tests/golden/116/GRCh38-custom/json_cache"})
    paths = rv.Paths(r, tmp_path / "stage", tmp_path / "out", tmp_path / "plugins")
    corpus = rv.REPO / "tests" / "golden" / "116" / "GRCh38-custom"
    inp, out = tmp_path / "in.vcf", tmp_path / "out" / "default.txt"
    assert paths.substitute("file={corpus}/custom/sites.vcf.gz", "vep-rs", inp, out) == f"file={corpus}/custom/sites.vcf.gz"
    assert paths.substitute("file={corpus}/custom/sites.vcf.gz", "ensembl-vep", inp, out) == "file=/work/corpus/custom/sites.vcf.gz"
    assert paths.mounts == {"/work/corpus": (corpus, False)}
    assert paths.substitute("{plugins}", "vep-rs", inp, out) == str(corpus / "plugin_data")
    assert paths.substitute("{plugins}", "ensembl-vep", inp, out) == "/work/plugins"
    assert paths.substitute("{stage}/x", "ensembl-vep", inp, out) == "/work/stage/x"
    assert paths.substitute("{input} {out}", "ensembl-vep", inp, out) == "/work/in/in.vcf /work/out/default.txt"
    assert paths.substitute("--plain", "ensembl-vep", inp, out) == "--plain"
    assert set(paths.mounts) == {"/work/corpus", "/work/plugins", "/work/stage"}
    assert paths.mounts["/work/stage"] == (tmp_path / "stage", True)
    with pytest.raises(rv.ReproductionError, match="pass --extra-dir loftee=<path>"):
        paths.substitute("{loftee}/LoF.pm", "ensembl-vep", inp, out)
    extras = rv.Paths(r, tmp_path / "stage", tmp_path / "out", tmp_path / "plugins", extras={"loftee": tmp_path / "loftee"})
    assert extras.substitute("loftee_path:{loftee}", "vep-rs", inp, out) == f"loftee_path:{tmp_path / 'loftee'}"
    assert extras.substitute("loftee_path:{loftee}", "ensembl-vep", inp, out) == "loftee_path:/work/extra/loftee"
    assert extras.mounts["/work/extra/loftee"] == (tmp_path / "loftee", False)
    without_corpus = rv.Paths(repro(run={"tool": "filter", "args": ["x"]}, cache={"json_cache": None, "assembly": None}), tmp_path / "stage", tmp_path / "out", None)
    with pytest.raises(rv.ReproductionError, match="names no directory"):
        without_corpus.substitute("{corpus}/x", "vep-rs", inp, out)
    with pytest.raises(rv.ReproductionError, match="needs --plugins-dir"):
        without_corpus.substitute("{plugins}", "ensembl-vep", inp, out)


def test_tool_directories_come_from_the_tool_dir_option(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    (tools / "samtools").write_text("")
    options = argparse.Namespace(tool_dirs=[tmp_path / "empty", tools])
    (tmp_path / "empty").mkdir()
    assert rv.tool_directories(repro(run={"tools": ["samtools"]}), options) == [tools]
    assert rv.tool_directories(repro(), options) == []
    with pytest.raises(rv.ReproductionError, match="needs the tool bcftools: pass --tool-dir"):
        rv.tool_directories(repro(run={"tools": ["bcftools"]}), options)
    assert "run.tools is not a list of strings" in rv.validate(repro(run={"tools": "samtools"}))


def test_docker_command_mounts_read_only_unless_writable_and_passes_stdin():
    paths = rv.Paths(repro(), Path("/stage"), Path("/out"), None)
    paths.mount(Path("/host/cache"), "/work/cache/homo_sapiens/116_GRCh38")
    paths.mount(Path("/host/out"), "/work/out", writable=True)
    cmd = rv.docker_command("img@sha256:0", paths, "vep", ["-i", "x"], {"PERL_HASH_SEED": "0"}, stdin=True, name="reevaluate-x-1")
    assert cmd[cmd.index("--name") + 1] == "reevaluate-x-1"
    assert cmd[:3] == ["docker", "run", "--rm"]
    assert "-i" in cmd and "-v" in cmd
    assert "/host/cache:/work/cache/homo_sapiens/116_GRCh38:ro" in cmd and "/host/out:/work/out" in cmd
    assert cmd[-4:] == ["img@sha256:0", "vep", "-i", "x"]
    assert "PERL_HASH_SEED=0" in cmd
    assert "-i" not in rv.docker_command("img@sha256:0", paths, "vep", [], {}, stdin=False)[:8]


def test_engine_args_follow_the_fixed_argument_lists(tmp_path):
    r = repro(cache={"fasta": "tests/golden/116/GRCh37-hgvs/reference.fa.gz", "assembly": "GRCh37", "json_cache": "tests/golden/116/GRCh37-hgvs/json_cache"}, run={"flags": ["--hgvs", "--plugin", "pLI"], "stats": True, "buffer_size": 10})
    paths = rv.Paths(r, tmp_path / "stage", tmp_path / "out", tmp_path / "plugins")
    inp, out = tmp_path / "in.vcf", tmp_path / "out" / "tab.txt"
    perl = tmp_path / "perl" / "homo_sapiens" / "116_GRCh37"
    vep_rs = rv.engine_args(r, "vep-rs", paths, inp, out, perl, "tab")
    assert vep_rs[:4] == ["-i", str(inp), "-o", str(out)]
    assert "--json_cache" in vep_rs and "--fork" not in vep_rs and "--stats_text" in vep_rs and "--no_stats" not in vep_rs
    assert vep_rs[vep_rs.index("--buffer_size") + 1] == "10" and vep_rs[-1] == "--tab"
    assert vep_rs[vep_rs.index("--fasta") + 1] == str(rv.REPO / "tests/golden/116/GRCh37-hgvs/reference.fa.gz")
    ensembl = rv.engine_args(r, "ensembl-vep", paths, inp, out, perl, "tab")
    assert ensembl[:4] == ["-i", "/work/in/in.vcf", "-o", "/work/out/tab.txt"]
    assert ensembl[ensembl.index("--dir_cache") + 1] == "/work/cache" and ensembl[ensembl.index("--cache_version") + 1] == "116"
    assert ensembl[ensembl.index("--fork") + 1] == "1" and ensembl[ensembl.index("--dir_plugins") + 1] == "/work/plugins"
    assert ensembl[ensembl.index("--fasta") + 1] == "/work/ref/reference.fa"
    assert (tmp_path / "stage" / "ref" / "reference.fa").exists() and (tmp_path / "stage" / "ref" / "reference.fa.fai").exists()
    assert "/work/cache/homo_sapiens/116_GRCh37" in paths.mounts and paths.mounts["/work/cache/homo_sapiens/116_GRCh37"] == (perl, False)
    grch37 = repro(run={"flags": ["--grch37"]})
    assert "--assembly" not in rv.engine_args(grch37, "vep-rs", paths, inp, out, perl, "default")
    unnamed = repro(run={"assembly_flag": False})
    assert "--assembly" not in rv.engine_args(unnamed, "ensembl-vep", paths, inp, out, perl, "default")
    guess = repro(input={"text": "21:100", "format": "guess"}, run={"format_flag": False})
    args = rv.engine_args(guess, "vep-rs", paths, inp, out, perl, "json")
    assert "--format" not in args and "--json" not in args
    database = repro(run={"perl_mode": "database"})
    db_args = rv.engine_args(database, "ensembl-vep", paths, inp, out, None, "default")
    assert db_args[db_args.index("--database") + 1] == "--db_version" and "--offline" not in db_args and "--dir_cache" not in db_args
    assert "--json_cache" in rv.engine_args(database, "vep-rs", paths, inp, out, None, "default")
    assert "run.perl_mode is not offline or database" in rv.validate(repro(run={"perl_mode": "remote"}))
    bare = rv.Paths(repro(run={"flags": ["--plugin", "pLI"]}), tmp_path / "stage", tmp_path / "out", None)
    with pytest.raises(rv.ReproductionError, match="needs --plugins-dir"):
        rv.engine_args(repro(run={"flags": ["--plugin", "pLI"]}), "ensembl-vep", bare, inp, out, perl, "default")
    own_dir = rv.engine_args(repro(run={"flags": ["--dir_plugins", "{stage}/mods", "--plugin", "LoF"]}), "ensembl-vep", bare, inp, out, perl, "default")
    assert own_dir.count("--dir_plugins") == 1
    gxf = repro(cache={"json_cache": "tests/golden/116/GRCh38-gtf-builder/json_cache", "gxf": "{corpus}/chr21_genes.gtf.gz"})
    assert rv.engine_args(gxf, "ensembl-vep", paths, inp, out, None, "default")[4:6] == ["--gtf", "/work/corpus/chr21_genes.gtf.gz"]
    assert "--json_cache" in rv.engine_args(gxf, "vep-rs", paths, inp, out, None, "default")


def test_the_committed_reproductions_are_valid_and_cover_the_page():
    """Every reproduction beside this test validates, is verified against the reference, and matches a row of the index page."""
    repros = rv.load_reproductions()
    assert repros, "no reproduction committed"
    for r in repros:
        assert rv.validate(r) == [], r.name
    assert rv.check(rv.PAGE, "116.2") == []
