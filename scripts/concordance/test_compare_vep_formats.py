"""Tests for compare_vep_formats.py: the four parsers on the committed Ensembl VEP 116.2
outputs of tests/golden/116/GRCh38-release116, the comparison rules, self-parity and the
warnings mode. The expected values are read off the committed files; the constants below
are the records the reference annotated, its consequence rows and the records it skipped."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_vep_formats as cvf  # noqa: E402

CORPUS = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-release116"
EXPECTED = CORPUS / "expected"
RECORDS = 36
ROWS = 977
SKIPPED_RECORDS = {12, 13, 15}


@pytest.fixture(scope="module")
def texts() -> dict[str, str]:
    return {
        "default": cvf.read_text(EXPECTED / "default.txt.gz"),
        "tab": cvf.read_text(EXPECTED / "tab.txt.gz"),
        "vcf": cvf.read_text(EXPECTED / "vcf.vcf.gz"),
        "json": cvf.read_text(EXPECTED / "json.jsonl.gz"),
    }


@pytest.fixture(scope="module")
def inputs() -> list[str]:
    return cvf.input_records(CORPUS / "variants.vcf")


@pytest.fixture(scope="module")
def manifest() -> dict:
    return json.loads((CORPUS / "manifest.json").read_text(encoding="utf-8"))


def data_lines(text: str) -> list[str]:
    return [line for line in cvf.text_lines(text) if not line.startswith("#")]


def header_lines(text: str) -> list[str]:
    return [line for line in cvf.text_lines(text) if line.startswith("#")]


def rebuild(text: str, rows: list[str]) -> str:
    return "\n".join(header_lines(text) + rows) + "\n"


def no_documented() -> cvf.Documented:
    return cvf.Documented()


def test_parse_default_reads_the_committed_reference(texts):
    entries = cvf.parse_default(texts["default"])
    assert len(entries) == ROWS
    first = entries[0]
    assert first.key == "21:11000001|N.|-|-"
    assert first.location == "21:11000001"
    assert first.consequence == "intergenic_variant"
    assert first.fields["Uploaded_variation"] == "bnd_intergenic_pair_a"
    assert first.fields["Extra"] == "IMPACT=MODIFIER"
    assert first.fields["IMPACT"] == "MODIFIER"
    ablation = next(e for e in entries if e.fields["Uploaded_variation"] == "del_10_6_mb" and e.feature == "ENST00000262354")
    assert ablation.consequence == "transcript_ablation"
    assert ablation.fields["Extra"] == "IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00"
    assert ablation.fields["OverlapBP"] == "23932"
    assert ablation.fields["OverlapPC"] == "100.00"
    assert ablation.fields["STRAND"] == "-1"
    assert len({e.key for e in entries}) == ROWS


def test_parse_tab_reads_the_committed_reference(texts):
    tab = cvf.parse_tab(texts["tab"])
    default = cvf.parse_default(texts["default"])
    assert len(tab) == ROWS
    assert [e.key for e in tab] == [e.key for e in default]
    assert "Extra" not in tab[0].fields
    assert set(tab[0].fields) == set(header_lines(texts["tab"])[-1][1:].split("\t"))
    assert all(t.fields["IMPACT"] == d.fields["IMPACT"] for t, d in zip(tab, default))


def test_parse_vcf_reads_the_committed_reference(texts, inputs):
    entries = cvf.parse_vcf(texts["vcf"], inputs)
    records = [e for e in entries if e.key.endswith("|record")]
    csq = [e for e in entries if not e.key.endswith("|record")]
    assert len(records) == RECORDS
    assert len(csq) == ROWS
    assert {e.record for e in records} == set(range(39)) - SKIPPED_RECORDS
    ins = next(e for e in records if e.fields["ID"] == "multi_ins_ins")
    assert ins.key == "rec11|record"
    assert ins.fields["INFO_without_CSQ"] == "END=30165780"
    assert ins.fields["CSQ_present"] == "yes"
    chunks = [e for e in csq if e.record == 11]
    assert [e.key for e in chunks] == [
        "rec11|insertion|ENST00000286808|Transcript",
        "rec11|insertion|ENST00000286808|Transcript#2",
    ]
    assert chunks[0].consequence == "feature_elongation&3_prime_UTR_variant"
    assert chunks[0].fields["EXON"] == "1/1"
    assert chunks[0].fields["cDNA_position"] == "1025-1026"
    assert chunks[0].fields["BIOTYPE"] == "protein_coding"
    assert chunks[0].fields["HGVSc"] == ""


def test_parse_json_reads_the_committed_reference(texts, inputs):
    entries = cvf.parse_json(texts["json"], inputs)
    records = [e for e in entries if e.key.endswith("|record")]
    consequences = [e for e in entries if not e.key.endswith("|record")]
    assert len(records) == RECORDS
    assert len(consequences) == ROWS
    ins = next(e for e in records if e.key == "rec11|record")
    assert ins.fields["id"] == "multi_ins_ins"
    assert ins.fields["most_severe_consequence"] == "feature_elongation"
    assert ins.fields["start"] == "30165781"
    assert ins.fields["end"] == "30165780"
    assert ins.fields["allele_string"] == "<INS>/<INS>"
    assert "input" not in ins.fields
    assert "transcript_consequences" not in ins.fields
    items = [e for e in consequences if e.record == 11]
    assert [e.key for e in items] == [
        "rec11|insertion|ENST00000286808|transcript",
        "rec11|insertion|ENST00000286808|transcript#2",
    ]
    assert items[0].consequence == "feature_elongation,3_prime_UTR_variant"
    assert items[0].fields["cdna_start"] == "1025"
    assert items[0].fields["cdna_end"] == "1026"
    assert items[0].fields["strand"] == "-1"
    vcf_keys = {(e.record, e.allele, e.feature) for e in cvf.parse_vcf(texts["vcf"], inputs) if e.record is not None}
    assert {(e.record, e.allele, e.feature) for e in entries} == vcf_keys


def test_parse_json_sets_aside_the_engine_identity_key(texts, inputs):
    """vep-rs writes ``JSON_IDENTITY_KEY`` on every JSON object, the input-only object of a
    skipped record included; the reference writes none, so the key is set aside and the
    entries read the same as the reference's."""
    reference = cvf.parse_json(texts["json"], inputs)
    identity = {cvf.JSON_IDENTITY_KEY: "0.5.0"}
    json_lines = [json.dumps({**json.loads(line), **identity}) for line in cvf.text_lines(texts["json"]) if line]
    json_lines.append(json.dumps({"input": "21\t1\tskipped\tA\tG\t.\t.\t.", **identity}))
    entries = cvf.parse_json("\n".join(json_lines) + "\n", inputs)
    annotated, skipped = entries[: len(reference)], entries[len(reference) :]
    assert [(e.key, e.fields) for e in annotated] == [(e.key, e.fields) for e in reference]
    assert all(cvf.JSON_IDENTITY_KEY not in e.fields for e in entries)
    assert [e.key for e in skipped] == [f"obj{RECORDS}|record"]
    assert skipped[0].fields == {}
    assert cvf.compare_entries(reference, annotated, no_documented())["verdict"] == "green"


def test_json_scalar_renders_like_the_rust_comparator():
    assert cvf.json_scalar("x") == "x"
    assert cvf.json_scalar(["a", "b"]) == "a,b"
    assert cvf.json_scalar(1) == "1"
    assert cvf.json_scalar(True) == "true"
    assert cvf.json_scalar(None) == "null"
    assert cvf.json_scalar(0.00001) == "1e-5"
    assert cvf.json_scalar(1e16) == "1e16"
    assert cvf.json_scalar(0.5) == "0.5"
    assert cvf.json_scalar({"name": "PF1", "db": "Pfam"}) == '{"db":"Pfam","name":"PF1"}'
    assert cvf.json_scalar([{"id": "rs1", "af": 0.25}]) == '{"af":0.25,"id":"rs1"}'


def test_ga4gh_vrs_compares_as_a_parsed_object():
    """Ensembl VEP writes the ``GA4GH_VRS`` object with its keys in hash order, so two rows whose
    objects agree compare equal whatever the order of their keys, in the default format's
    ``Extra`` and as a tab column alike; a value that is not JSON is compared as written."""
    perl = ('{"type":"Allele","state":{"sequence":"C","type":"SequenceState"},'
            '"location":{"type":"SequenceLocation","sequence_id":"refseq:NC_000021.9",'
            '"interval":{"type":"SimpleInterval","end":10538885,"start":"10538884"}}}')
    rust = ('{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval",'
            '"start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},'
            '"state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}')
    assert cvf.canonical_json_text(perl) == cvf.canonical_json_text(rust)
    assert cvf.canonical_json_text("21:10538884:A:C") == "21:10538884:A:C"
    assert cvf.canonical_json_text("-") == "-"
    assert cvf.canonical_json_text("1e5") == "1e5"
    header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    row = "rs1\t21:10538885\tC\tG1\tT1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;GA4GH_VRS={};STRAND=1\n"
    perl_entry = cvf.parse_default(header + row.format(perl))[0]
    rust_entry = cvf.parse_default(header + row.format(rust))[0]
    assert perl_entry.fields["GA4GH_VRS"] == rust_entry.fields["GA4GH_VRS"] == cvf.canonical_json_text(perl)
    assert perl_entry.fields["Extra"] == rust_entry.fields["Extra"]
    assert perl_entry.fields["STRAND"] == "1"
    # A pair whose value is not an object, and a bare key, are rebuilt as written.
    dash_row = row.format("-")
    bare_key_row = row.replace("GA4GH_VRS={};", "GA4GH_VRS;")
    dash_entry, bare_key_entry = cvf.parse_default(header + dash_row + bare_key_row)
    assert dash_entry.fields["Extra"] == "IMPACT=MODIFIER;GA4GH_VRS=-;STRAND=1"
    assert dash_entry.fields["GA4GH_VRS"] == "-"
    assert bare_key_entry.fields["Extra"] == "IMPACT=MODIFIER;GA4GH_VRS;STRAND=1"
    assert bare_key_entry.fields["GA4GH_VRS"] == ""
    tab_header = "#Uploaded_variation\tLocation\tAllele\tGA4GH_VRS\n"
    tab_row = "rs1\t21:10538885\tC\t{}\n"
    perl_tab = cvf.parse_tab(tab_header + tab_row.format(perl))[0]
    rust_tab = cvf.parse_tab(tab_header + tab_row.format(rust))[0]
    assert perl_tab.fields["GA4GH_VRS"] == rust_tab.fields["GA4GH_VRS"]
    # Self-parity: the JSON format nests the object with both bounds numbers, where the text
    # formats carry the start as a string; the VCF format is refused with --ga4gh_vrs, so its
    # CSQ declares no such column.
    nested = json.loads(rust)
    nested["location"]["interval"]["start"] = 10538884
    assert cvf.numberified_json_text(perl) == cvf.json_compact(nested)
    assert cvf.numberified_json_text("21:10538884:A:C") == "21:10538884:A:C"
    # A superscript digit is isdigit() but not a decimal int() accepts; it stays a string.
    assert cvf.numberified_json_text('{"start":"\u00b2"}') == '{"start":"\u00b2"}'
    csq = ("Allele", "Consequence", "IMPACT", "Gene", "Feature_type", "Feature", "STRAND")
    vcf = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: '
        + "|".join(csq) + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        + "21\t10538885\trs1\tA\tC\t.\t.\tCSQ=C|intron_variant|MODIFIER|G1|Transcript|T1|1\n"
    )
    obj = {
        "input": "21\t10538885\trs1\tA\tC\t.\t.\t.", "id": "rs1", "seq_region_name": "21",
        "start": 10538885, "end": 10538885, "strand": 1, "allele_string": "A/C",
        "most_severe_consequence": "intron_variant",
        "transcript_consequences": [{
            "variant_allele": "C", "gene_id": "G1", "transcript_id": "T1", "consequence_terms": ["intron_variant"],
            "impact": "MODIFIER", "strand": 1, "ga4gh_vrs": nested,
        }],
    }
    tab = (
        tab_header.rstrip("\n") + "\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\tSTRAND\n"
        + tab_row.format(rust).rstrip("\n") + "\tG1\tT1\tTranscript\tintron_variant\tMODIFIER\t1\n"
    )
    report = cvf.SelfParity(header + row.format(rust), tab, vcf, json.dumps(obj) + "\n").run()
    assert report["verdict"] == "green", report["reasons"]
    assert "GA4GH_VRS" in report["pairs"]["json"]["columns_compared"]
    nested["location"]["interval"]["end"] = 10538886
    report = cvf.SelfParity(header + row.format(rust), tab, vcf, json.dumps(obj) + "\n").run()
    assert sorted(report["pairs"]["json"]["differing"]) == ["GA4GH_VRS"]


def test_dedup_keys_counts_exact_duplicates():
    entries = [cvf.Entry("k", None, None, "", "", "", {}) for _ in range(3)]
    entries.append(cvf.Entry("other", None, None, "", "", "", {}))
    cvf.dedup_keys(entries)
    assert [e.key for e in entries] == ["k", "k#2", "k#3", "other"]


def test_malformed_lines_are_errors(inputs):
    with pytest.raises(cvf.ComparatorError):
        cvf.parse_default("#A\tB\n1\t2\t3\n")
    with pytest.raises(cvf.ComparatorError):
        cvf.parse_vcf("##INFO=<ID=CSQ,Description=\"Format: A|B\">\n1\t2\t.\tA\tT\t.\t.\tCSQ=x\n", inputs)
    with pytest.raises(cvf.ComparatorError):
        cvf.parse_json("[1, 2]\n", inputs)


def test_verdict_is_red_on_any_difference(texts):
    expected = cvf.parse_default(texts["default"])
    rows = data_lines(texts["default"])
    assert cvf.compare_entries(expected, cvf.parse_default(texts["default"]), no_documented())["verdict"] == "green"
    mutated = rows[:]
    mutated[0] = mutated[0].replace("IMPACT=MODIFIER", "IMPACT=LOW")
    report = cvf.compare_entries(expected, cvf.parse_default(rebuild(texts["default"], mutated)), no_documented())
    assert report["verdict"] == "red"
    assert report["one_sided"] == {"expected_only": 0, "actual_only": 0}
    assert report["columns"]["IMPACT"]["differing"] == 1
    assert report["columns"]["Extra"]["differing"] == 1
    assert report["columns"]["IMPACT"]["examples"] == {
        "intergenic_variant": ["21:11000001|N.|-|-: expected 'MODIFIER', got 'LOW'"]
    }
    assert "field IMPACT: 1 differing entries" in report["reasons"]


def test_one_sided_counts(texts):
    expected = cvf.parse_default(texts["default"])
    rows = data_lines(texts["default"])
    dropped = cvf.compare_entries(expected, cvf.parse_default(rebuild(texts["default"], rows[1:])), no_documented())
    assert dropped["one_sided"] == {"expected_only": 1, "actual_only": 0}
    assert dropped["missing"] == {"intergenic_variant": {"count": 1, "examples": ["21:11000001|N.|-|-"]}}
    extra_row = rows[0].replace("21:11000001", "21:11000002")
    added = cvf.compare_entries(expected, cvf.parse_default(rebuild(texts["default"], rows + [extra_row])), no_documented())
    assert added["one_sided"] == {"expected_only": 0, "actual_only": 1}
    assert added["unexpected"]["intergenic_variant"]["examples"] == ["21:11000002|N.|-|-"]
    assert added["verdict"] == "red"
    assert added["counts"]["expected_entries"] == ROWS
    assert added["counts"]["actual_entries"] == ROWS + 1


def test_mask_application(texts, manifest):
    """The manifest documents the two mate-side rows of record 1 under
    breakend_mate_local_read: ENST00000334068 as upstream_gene_variant where the reference
    prints intergenic_variant, ENST00000390690 as coding_sequence_variant,feature_truncation
    where it prints feature_truncation."""
    expected = cvf.parse_default(texts["default"])
    rows = data_lines(texts["default"])
    mutated = rows[:]
    target = next(i for i, r in enumerate(rows) if r.startswith("bnd_intergenic_with_genic_mate\t") and "ENST00000334068" in r)
    mutated[target] = mutated[target].replace("\tintergenic_variant\t", "\tupstream_gene_variant\t")
    other = next(i for i, r in enumerate(rows) if r.startswith("bnd_intergenic_with_genic_mate\t") and "ENST00000390690" in r)
    mutated[other] = mutated[other].replace("\tfeature_truncation\t", "\tfeature_truncation,coding_sequence_variant\t")
    actual = cvf.parse_default(rebuild(texts["default"], mutated))

    unmasked = cvf.compare_entries(expected, actual, cvf.Documented.from_manifest(manifest, set()))
    assert unmasked["verdict"] == "red"
    assert unmasked["columns"]["Consequence"]["differing"] == 2

    masked = cvf.Documented.from_manifest(manifest, {"breakend_mate_local_read"})
    assert masked.applied == {"breakend_mate_local_read": 2}
    assert masked.available == {"breakend_mate_local_read": 2, "reference_skipped_record": 5}
    report = cvf.compare_entries(expected, actual, masked)
    assert report["verdict"] == "green"
    assert report["counts"]["documented_divergences_skipped"] == 2
    assert report["counts"]["compared"] == ROWS - 2

    mutated[target] = mutated[target].replace("\tupstream_gene_variant\t", "\tdownstream_gene_variant\t")
    stale = cvf.compare_entries(expected, cvf.parse_default(rebuild(texts["default"], mutated)), masked)
    assert stale["verdict"] == "red"
    assert len(stale["stale_divergences"]) == 1
    assert "documented divergences not observed: 1" in stale["reasons"]


def test_reference_skipped_record_rows(texts, manifest, inputs):
    """The manifest documents five vep-rs-only keys on records the reference skipped; the
    candidate's consequence rows on them are set aside only under the class's mask, and its
    record-level VCF line on such a record is counted, not compared."""
    skipped_row = (
        "multi_ins_del\t21:30165801-30166000\tinsertion\tENSG00000156282\tENST00000286808\tTranscript\t"
        "feature_elongation,coding_sequence_variant,3_prime_UTR_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1"
    )
    rows = data_lines(texts["default"]) + [skipped_row]
    expected = cvf.parse_default(texts["default"])
    actual = cvf.parse_default(rebuild(texts["default"], rows))
    assert cvf.compare_entries(expected, actual, cvf.Documented.from_manifest(manifest, set()))["one_sided"]["actual_only"] == 1
    masked = cvf.Documented.from_manifest(manifest, {"reference_skipped_record"})
    assert cvf.compare_entries(expected, actual, masked)["verdict"] == "green"

    vcf_rows = data_lines(texts["vcf"]) + [inputs[12]]
    vcf_actual = cvf.parse_vcf(rebuild(texts["vcf"], vcf_rows), inputs)
    vcf_expected = cvf.parse_vcf(texts["vcf"], inputs)
    report = cvf.compare_entries(vcf_expected, vcf_actual, masked)
    assert report["verdict"] == "green"
    assert report["counts"]["record_entries_on_reference_skipped_records"] == 1
    assert cvf.compare_entries(vcf_expected, vcf_actual, no_documented())["one_sided"]["actual_only"] == 1


def test_field_divergence_is_held_to_the_documented_value(texts, manifest):
    documented = dict(manifest)
    documented["field_divergences"] = [{
        "location": "21:11000001", "allele": "N.", "feature": "-", "feature_type": "-",
        "record_indices": [0], "field": "IMPACT", "vep": "MODIFIER", "vep_rs": "LOW",
        "expected_divergence": "impact_example",
    }]
    expected = cvf.parse_default(texts["default"])
    rows = data_lines(texts["default"])
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=LOW")
    actual = cvf.parse_default(rebuild(texts["default"], rows))
    masked = cvf.Documented.from_manifest(documented, {"impact_example"})
    assert cvf.compare_entries(expected, actual, masked)["verdict"] == "green"
    rows[0] = rows[0].replace("IMPACT=LOW", "IMPACT=HIGH")
    stale = cvf.compare_entries(expected, cvf.parse_default(rebuild(texts["default"], rows)), masked)
    assert stale["verdict"] == "red"
    assert stale["stale_divergences"][0].startswith("21:11000001|N.|-|-: field IMPACT documented as 'LOW'")


def test_floor_failure(texts):
    expected = cvf.parse_default(texts["default"])
    actual = cvf.parse_default(texts["default"])
    met = cvf.compare_entries(expected, actual, no_documented(), floors={"IMPACT": ROWS, "OverlapBP": 1})
    assert met["verdict"] == "green"
    assert met["columns"]["IMPACT"] == {
        "expected_populated": ROWS, "actual_populated": ROWS, "differing": 0, "pairs_accepted": 0,
        "floor": ROWS, "floor_met": True,
    }
    empty = cvf.compare_entries(expected, actual, no_documented(), floors={"HGVSc": 1})
    assert empty["verdict"] == "red"
    assert empty["reasons"] == ["floor HGVSc=1: expected populated 0, actual populated 0"]
    short = cvf.compare_entries(expected, actual, no_documented(), floors={"IMPACT": ROWS + 1})
    assert short["verdict"] == "red"
    assert short["columns"]["IMPACT"]["floor_met"] is False


def test_field_pairs_accept_a_value_pair(texts):
    expected = cvf.parse_default(texts["default"])
    rows = data_lines(texts["default"])
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=LOW")
    actual = cvf.parse_default(rebuild(texts["default"], rows))
    pairs = cvf.parse_field_pairs(["IMPACT=MODIFIER=>LOW;Extra=IMPACT=MODIFIER=>IMPACT=LOW"])
    assert pairs == {"IMPACT": {("MODIFIER", "LOW")}, "Extra": {("IMPACT=MODIFIER", "IMPACT=LOW")}}
    accepted = cvf.compare_entries(expected, actual, no_documented(), field_pairs=pairs)
    assert accepted["verdict"] == "green"
    assert accepted["columns"]["IMPACT"]["pairs_accepted"] == 1
    other = cvf.compare_entries(expected, actual, no_documented(), field_pairs={"IMPACT": {("MODIFIER", "HIGH")}})
    assert other["verdict"] == "red"


def test_order_difference_is_red_except_for_parquet(texts):
    expected = cvf.parse_default(texts["default"])
    rows = data_lines(texts["default"])
    swapped = [rows[1], rows[0]] + rows[2:]
    actual = cvf.parse_default(rebuild(texts["default"], swapped))
    ordered = cvf.compare_entries(expected, actual, no_documented())
    assert ordered["verdict"] == "red"
    assert ordered["order"] == {"first_divergence": 0, "expected": "21:11000001|N.|-|-", "actual": "21:11000001|N[21:11500000[|-|-"}
    unordered = cvf.compare_entries(expected, actual, no_documented(), check_order=False)
    assert unordered["verdict"] == "green"
    assert unordered["order"] is None


def test_parquet_actual_goes_through_the_adapter(texts, monkeypatch, tmp_path):
    rows = data_lines(texts["tab"])
    reordered = rebuild(texts["tab"], list(reversed(rows)))
    calls = []

    def fake_adapter(path):
        calls.append(path)
        return reordered

    monkeypatch.setattr(cvf, "parquet_as_tab", fake_adapter)
    expected_file = tmp_path / "tab.txt"
    expected_file.write_text(texts["tab"], encoding="utf-8")
    report = cvf.compare_files(
        "parquet", expected_file, tmp_path / "out.parquet", None, no_documented(), set(), {}, {}, 4
    )
    assert calls == [tmp_path / "out.parquet"]
    assert report["verdict"] == "green"
    assert report["format"] == "parquet"
    assert report["counts"]["compared"] == ROWS


def test_vcf_and_json_formats_need_the_input(inputs):
    for fmt, name in (("vcf", "vcf.vcf"), ("json", "json.jsonl")):
        with pytest.raises(cvf.ComparatorError):
            cvf.parse_format(fmt, EXPECTED / f"{name}.gz", None)
    assert len(cvf.parse_format("vcf", EXPECTED / "vcf.vcf.gz", inputs)) == RECORDS + ROWS


def self_parity(texts, **override) -> dict:
    kwargs = {k: texts[k] for k in ("default", "tab", "vcf", "json")}
    kwargs.update(override)
    return cvf.SelfParity(kwargs["default"], kwargs["tab"], kwargs["vcf"], kwargs["json"]).run()


def test_self_parity_green_on_the_committed_formats(texts):
    report = self_parity(texts)
    assert report["verdict"] == "green", report["reasons"]
    assert report["records"] == RECORDS
    assert report["default_entries"] == ROWS
    for name in ("tab", "vcf", "json"):
        pair = report["pairs"][name]
        assert pair["shared_entries"] == ROWS
        assert pair["one_sided"] == {"default_only": 0, f"{name}_only": 0}
        assert pair["differing"] == {}
    assert "OverlapBP" in report["pairs"]["json"]["columns_compared"]
    assert "OverlapBP" not in report["pairs"]["tab"]["columns_compared"]
    assert report["not_compared"]["vcf_implied"] == ["BIOTYPE", "EXON", "HGNC_ID", "INTRON", "SYMBOL", "SYMBOL_SOURCE"]
    assert report["not_compared"]["tab_undeclared"] == ["OverlapBP", "OverlapPC"]
    assert report["not_compared"]["vcf_undeclared"] == ["OverlapBP", "OverlapPC"]
    assert report["not_compared"]["default_record_level"] == []


def test_self_parity_red_on_a_mutated_copy(texts):
    tab_rows = data_lines(texts["tab"])
    tab_rows[0] = tab_rows[0].replace("\tMODIFIER\t", "\tLOW\t")
    report = self_parity(texts, tab=rebuild(texts["tab"], tab_rows))
    assert report["verdict"] == "red"
    assert report["reasons"] == ["tab: column IMPACT differs on 1 entries"]
    assert report["pairs"]["tab"]["differing"]["IMPACT"]["count"] == 1

    json_lines = cvf.text_lines(texts["json"])
    obj = json.loads(json_lines[0])
    obj["intergenic_consequences"].pop()
    json_lines[0] = json.dumps(obj)
    report = self_parity(texts, json="\n".join(json_lines) + "\n")
    assert report["verdict"] == "red"
    assert report["pairs"]["json"]["one_sided"] == {"default_only": 1, "json_only": 0}

    vcf_rows = data_lines(texts["vcf"])
    vcf_rows[0] = vcf_rows[0].replace("|intergenic_variant|MODIFIER|", "|downstream_gene_variant|MODIFIER|", 1)
    report = self_parity(texts, vcf=rebuild(texts["vcf"], vcf_rows))
    assert report["verdict"] == "red"
    assert report["pairs"]["vcf"]["differing"]["Consequence"]["count"] == 1


def test_self_parity_checks_most_severe_consequence(texts):
    json_lines = cvf.text_lines(texts["json"])
    obj = json.loads(json_lines[0])
    obj["most_severe_consequence"] = "stop_gained"
    json_lines[0] = json.dumps(obj)
    report = self_parity(texts, json="\n".join(json_lines) + "\n")
    assert report["verdict"] == "red"
    assert report["record_checks"]["problems"] == [
        "bnd_intergenic_pair_a|21:11000001: most_severe_consequence 'stop_gained', the object's terms give 'intergenic_variant'"
    ]


def test_self_parity_accepts_a_name_carrying_a_bar(texts):
    renamed = {k: v.replace("bnd_intergenic_pair_a", "bnd|pair_a") for k, v in texts.items()}
    report = self_parity(renamed)
    assert report["verdict"] == "green", report["reasons"]
    assert report["pairs"]["json"]["shared_entries"] == ROWS


def test_self_parity_undoes_the_json_sift_polyphen_and_domains_splits():
    """A synthetic record, since no committed corpus runs --sift, --polyphen or --domains:
    the JSON prediction/score pair and the db/name objects must match the default
    format's ``deleterious(0.01)`` and ``Pfam:PF1,PANTHER:PTHR1`` spellings."""
    fixed = ("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
             "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation")
    values = ("rs1", "1:100", "G", "ENSG1", "ENST1", "Transcript", "missense_variant", "10", "5", "2", "A/G", "gCa/gGa", "rs1")
    extra = {"IMPACT": "MODERATE", "STRAND": "1", "SIFT": "deleterious(0.01)",
             "PolyPhen": "probably_damaging(0.95)", "DOMAINS": "Pfam:PF1,PANTHER:PTHR1"}
    default = "#" + "\t".join(fixed + ("Extra",)) + "\n" + "\t".join(values + (";".join(f"{k}={v}" for k, v in extra.items()),)) + "\n"
    tab = "#" + "\t".join(fixed + tuple(extra)) + "\n" + "\t".join(values + tuple(extra.values())) + "\n"
    csq = ("Allele", "Consequence", "IMPACT", "Gene", "Feature_type", "Feature", "cDNA_position", "CDS_position",
           "Protein_position", "Amino_acids", "Codons", "Existing_variation", "STRAND", "SIFT", "PolyPhen", "DOMAINS")
    row = dict(zip(fixed, values)) | extra
    vcf = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: '
        + "|".join(csq) + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        + "1\t100\trs1\tA\tG\t.\t.\tCSQ=" + "|".join(cvf.vcf_escape(row[c]) for c in csq) + "\n"
    )
    obj = {
        "input": "1\t100\trs1\tA\tG\t.\t.\t.", "id": "rs1", "seq_region_name": "1", "start": 100, "end": 100,
        "strand": 1, "allele_string": "A/G", "most_severe_consequence": "missense_variant",
        "colocated_variants": [{"id": "rs1"}],
        "transcript_consequences": [{
            "variant_allele": "G", "gene_id": "ENSG1", "transcript_id": "ENST1", "consequence_terms": ["missense_variant"],
            "impact": "MODERATE", "strand": 1, "cdna_start": 10, "cdna_end": 10, "cds_start": 5, "cds_end": 5,
            "protein_start": 2, "protein_end": 2, "amino_acids": "A/G", "codons": "gCa/gGa",
            "sift_prediction": "deleterious", "sift_score": 0.01,
            "polyphen_prediction": "probably_damaging", "polyphen_score": 0.95,
            "domains": [{"db": "Pfam", "name": "PF1"}, {"db": "PANTHER", "name": "PTHR1"}],
        }],
    }
    report = cvf.SelfParity(default, tab, vcf, json.dumps(obj) + "\n").run()
    assert report["verdict"] == "green", report["reasons"]
    assert {"SIFT", "PolyPhen", "DOMAINS"} <= set(report["pairs"]["json"]["columns_compared"])
    obj["transcript_consequences"][0]["domains"].pop()
    obj["transcript_consequences"][0]["sift_score"] = 0.02
    report = cvf.SelfParity(default, tab, vcf, json.dumps(obj) + "\n").run()
    assert sorted(report["pairs"]["json"]["differing"]) == ["DOMAINS", "SIFT"]


def test_perl_numify_follows_perl_numification():
    assert cvf.perl_numify("30425949") == "30425949"
    assert cvf.perl_numify("+30425949") == "30425949"
    assert cvf.perl_numify("30_425_949") == "30"
    assert cvf.perl_numify(".") == "0"
    assert cvf.perl_numify("0x1D04A1D") == "0"
    assert cvf.perl_numify("30425949.5") == "30425949.5"
    assert cvf.perl_numify("1e999") == "Inf"


def test_uploaded_variation_and_location_follow_the_output_factory():
    named = {"id": "rs1", "seq_region_name": "21", "start": 100, "end": 100, "allele_string": "A/G",
             "input": "21\t100\trs1\tA\tG\t.\t.\t."}
    assert cvf.uploaded_variation(named) == "rs1"
    assert cvf.record_location(named) == "21:100"
    unnamed = {"id": ".", "seq_region_name": "21", "start": 10795228, "end": 10795228, "allele_string": "C/G/T",
               "input": "21\t10795228\t.\tC\tG,T\t.\t.\t."}
    assert cvf.uploaded_variation(unnamed) == "21_10795228_C/G/T"
    insertion = {"id": ".", "seq_region_name": "21", "start": 9489855, "end": 9489854, "allele_string": "-/AA/AC",
                 "input": "21\t9489854\t.\tT\tTAA,TAC\t.\t.\t."}
    assert cvf.uploaded_variation(insertion) == "21_9489854_T/TAA/TAC"
    assert cvf.record_location(insertion) == "21:9489854-9489855"
    assert cvf.format_coords(5, None) == "5-?"
    assert cvf.format_coords(None, 7) == "?-7"


def test_feature_kind_vcf_escape_and_most_severe():
    assert cvf.feature_kind("Transcript") == "transcript"
    assert cvf.feature_kind("RegulatoryFeature") == "regulatory_feature"
    assert cvf.feature_kind("MotifFeature") == "motif_feature"
    assert cvf.feature_kind("-") == "intergenic"
    assert cvf.feature_kind("regulatory_feature") == "regulatory_feature"
    assert cvf.vcf_escape("a,b;c d|e") == "a&b%3Bc_d&e"
    assert cvf.reverse_complement("ACGT") == "ACGT"
    assert cvf.reverse_complement("AAC") == "GTT"
    assert cvf.reverse_complement("deletion") == "deletion"
    assert cvf.most_severe({"intron_variant", "splice_region_variant", "?"}) == "splice_region_variant"
    assert cvf.most_severe({"?"}) == "?"
    assert cvf.values_equal("100.00", "100")
    assert not cvf.values_equal("1/1", "1")


def test_warnings_equal_after_masking(tmp_path):
    expected = tmp_path / "ref_warnings.txt"
    actual = tmp_path / "cand_warnings.txt"
    expected.write_text(
        "WARNING: line 13 skipped (21 30165800 multi_ins_del N <INS>,<DEL> . . END=3016...): "
        "Invalid allele string N/<INS>,<DEL> or possible parsing error\n"
        "WARNING: line 16 skipped (21 30165850 del_without_end N <DEL> . . SVTYPE=DEL): start > end+1 : (START=30165851, END=30165850)\n"
        "WARNING: Unable to open /data/ref/one.fa for reading\n",
        encoding="utf-8",
    )
    actual.write_text(
        "WARNING: line 16 skipped (21\t30165850\tdel_without_end\tN\t<DEL>): start > end+1 : (START=30165851, END=30165850)\n"
        "WARNING: line 13 skipped (): Invalid allele string N/<INS>,<DEL> or possible parsing error\n"
        "WARNING: Unable to open ./two.fa for reading\n",
        encoding="utf-8",
    )
    report = cvf.compare_warnings(expected, actual, None, None)
    assert report["verdict"] == "green"
    assert report["warnings"]["shared"] == 3
    assert cvf.parse_warning_lines(expected.read_text()) == {
        (13, "Invalid allele string N/<INS>,<DEL> or possible parsing error"),
        (16, "start > end+1 : (START=30165851, END=30165850)"),
        (None, "Unable to open <path> for reading"),
    }
    actual.write_text("WARNING: line 13 skipped (x): Invalid allele string N/<INS>,<DEL> or possible parsing error\n", encoding="utf-8")
    report = cvf.compare_warnings(expected, actual, None, None)
    assert report["verdict"] == "red"
    assert report["warnings"]["expected_only"] == [
        {"line": 16, "description": "start > end+1 : (START=30165851, END=30165850)"},
        {"line": None, "description": "Unable to open <path> for reading"},
    ]
    assert report["reasons"] == ["warnings: 2 entries only in the expected file"]


def test_skipped_variants_files_compare_by_line_and_description(tmp_path):
    expected = tmp_path / "ref_skipped.txt"
    actual = tmp_path / "cand_skipped.txt"
    expected.write_text(
        "[VEP skipped the following variants from /data/in/variants.vcf]\n"
        "Line 13   \t21\t30165800\tmulti_ins_del\tN\t<INS>,<DEL>\tInvalid allele string N/<INS>,<DEL> or possible parsing error\n"
        "Line 16   \t21 30165850 del_without_end N <DEL> . . SVTYPE=DEL   \tstart > end+1 : (START=30165851, END=30165850)\n",
        encoding="utf-8",
    )
    actual.write_text(
        "[VEP skipped the following variants from variants.vcf]\n"
        "Line 13\t21 30165800 multi_ins_del N <INS>,<DEL> . . END=30166000\tInvalid allele string N/<INS>,<DEL> or possible parsing error\n"
        "Line 16\tx\tstart > end+1 : (START=30165851, END=30165850)\n",
        encoding="utf-8",
    )
    warnings = tmp_path / "w.txt"
    warnings.write_text("", encoding="utf-8")
    report = cvf.compare_warnings(warnings, warnings, expected, actual)
    assert report["verdict"] == "green"
    assert report["skipped"]["shared"] == 2
    with pytest.raises(cvf.ComparatorError):
        cvf.compare_warnings(warnings, warnings, expected, None)


def test_cli_writes_the_report_and_exits_by_verdict(tmp_path, capsys):
    report = tmp_path / "r.json"
    code = cvf.main([
        "--format", "default", "--expected", str(EXPECTED / "default.txt.gz"),
        "--actual", str(EXPECTED / "default.txt.gz"), "--report", str(report),
    ])
    assert code == 0
    written = json.loads(report.read_text(encoding="utf-8"))
    assert written["verdict"] == "green"
    assert written["mode"] == "compare"
    assert written["format"] == "default"
    assert capsys.readouterr().out == "default: green\n"

    mutated = tmp_path / "mutated.txt"
    text = cvf.read_text(EXPECTED / "default.txt.gz")
    mutated.write_text(text.replace("IMPACT=MODIFIER", "IMPACT=LOW", 1), encoding="utf-8")
    code = cvf.main([
        "--format", "default", "--expected", str(EXPECTED / "default.txt.gz"),
        "--actual", str(mutated), "--report", str(report), "--floors", "IMPACT=1",
    ])
    assert code == 1
    assert json.loads(report.read_text(encoding="utf-8"))["verdict"] == "red"
    out = capsys.readouterr().out
    assert out.startswith("default: red\n")
    assert "  field IMPACT: 1 differing entries" in out


def test_cli_self_parity_and_warnings_modes(tmp_path):
    report = tmp_path / "sp.json"
    code = cvf.main([
        "--self-parity", "--default", str(EXPECTED / "default.txt.gz"), "--tab", str(EXPECTED / "tab.txt.gz"),
        "--vcf", str(EXPECTED / "vcf.vcf.gz"), "--json", str(EXPECTED / "json.jsonl.gz"),
        "--input", str(CORPUS / "variants.vcf"), "--report", str(report),
    ])
    assert code == 0
    written = json.loads(report.read_text(encoding="utf-8"))
    assert written["mode"] == "self-parity"
    assert written["verdict"] == "green"
    assert written["files"]["input"] == str(CORPUS / "variants.vcf")
    empty = tmp_path / "w.txt"
    empty.write_text("", encoding="utf-8")
    assert cvf.main(["--warnings", "--expected", str(empty), "--actual", str(empty), "--report", str(report)]) == 0
    assert json.loads(report.read_text(encoding="utf-8"))["format"] == "warnings"


def test_cli_errors_exit_2_with_one_error_line(tmp_path, capsys):
    code = cvf.main([
        "--format", "vcf", "--expected", str(EXPECTED / "vcf.vcf.gz"), "--actual", str(EXPECTED / "vcf.vcf.gz"),
    ])
    assert code == 2
    err = capsys.readouterr().err
    assert err == "ERROR: [compare_vep_formats] --format vcf needs --input, the input VCF\n"
    assert cvf.main(["--self-parity", "--default", str(EXPECTED / "default.txt.gz")]) == 2
    assert "needs --tab, --vcf, --json" in capsys.readouterr().err
    assert cvf.main(["--format", "default", "--expected", str(tmp_path / "none.txt"), "--actual", str(tmp_path / "none.txt")]) == 2
    assert capsys.readouterr().err.startswith("ERROR: [compare_vep_formats] cannot read ")
    assert cvf.main([]) == 2


def test_field_pairs_reach_the_extra_column(texts):
    """A pair accepted on an Extra key also settles the raw Extra text that carries it; an
    Extra whose keys are only reordered stays a difference, as every key still agrees."""
    expected = cvf.parse_default(texts["default"])
    rows = [r.replace("OverlapPC=100.00", "OverlapPC=100") for r in data_lines(texts["default"])]
    actual = cvf.parse_default(rebuild(texts["default"], rows))
    changed = sum(1 for r in data_lines(texts["default"]) if "OverlapPC=100.00" in r)
    assert changed > 0
    report = cvf.compare_entries(expected, actual, no_documented(), field_pairs={"OverlapPC": {("100.00", "100")}})
    assert report["verdict"] == "green", report["reasons"]
    assert report["columns"]["OverlapPC"]["pairs_accepted"] == changed
    assert report["columns"]["Extra"]["pairs_accepted"] == changed
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=LOW")
    report = cvf.compare_entries(
        expected, cvf.parse_default(rebuild(texts["default"], rows)), no_documented(), field_pairs={"OverlapPC": {("100.00", "100")}}
    )
    assert report["columns"]["Extra"]["differing"] == 1
    assert report["columns"]["IMPACT"]["differing"] == 1
    ablation_extra = "IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00"
    reordered = [r.replace(ablation_extra, "STRAND=-1;IMPACT=HIGH;OverlapBP=23932;OverlapPC=100.00") for r in data_lines(texts["default"])]
    assert sum(1 for r in data_lines(texts["default"]) if r.endswith("\t" + ablation_extra)) == 1
    report = cvf.compare_entries(expected, cvf.parse_default(rebuild(texts["default"], reordered)), no_documented())
    assert report["verdict"] == "red"
    assert report["columns"]["Extra"]["differing"] == 1
    assert report["columns"]["Extra"]["pairs_accepted"] == 0
    assert all(report["columns"][c]["differing"] == 0 for c in ("IMPACT", "STRAND", "OverlapBP", "OverlapPC"))
    pairs = {"OverlapPC": {("100.00", "100")}}
    assert cvf.extra_pairs_accepted("IMPACT=HIGH;OverlapPC=100.00", "IMPACT=HIGH;OverlapPC=100", pairs)
    assert not cvf.extra_pairs_accepted("IMPACT=HIGH;OverlapPC=100.00", "OverlapPC=100;IMPACT=HIGH", pairs)


def test_unreadable_inputs_and_an_unwritable_report_exit_2(tmp_path, capsys):
    whole = (EXPECTED / "default.txt.gz").read_bytes()
    truncated = tmp_path / "truncated.txt.gz"
    truncated.write_bytes(whole[: len(whole) // 2])
    corrupt = tmp_path / "corrupt.txt.gz"
    corrupt.write_bytes(whole[:2000] + bytes(b ^ 0xFF for b in whole[2000:2400]) + whole[2400:])
    binary = tmp_path / "binary.txt"
    binary.write_bytes(whole)
    for actual in (truncated, corrupt, binary):
        code = cvf.main(["--format", "default", "--expected", str(EXPECTED / "default.txt.gz"), "--actual", str(actual)])
        assert code == 2
        err = capsys.readouterr().err
        assert err.startswith(f"ERROR: [compare_vep_formats] cannot read {actual}: ")
        assert err.count("\n") == 1
    code = cvf.main([
        "--format", "default", "--expected", str(EXPECTED / "default.txt.gz"), "--actual", str(EXPECTED / "default.txt.gz"),
        "--report", str(binary / "sub" / "report.json"),
    ])
    assert code == 2
    assert capsys.readouterr().err.startswith("ERROR: [compare_vep_formats] cannot write ")


def test_self_parity_keeps_the_dash_of_used_ref_and_given_ref_and_skips_an_input_only_object(texts):
    """JSON.pm:327 keeps `-` for Allele, USED_REF and GIVEN_REF where every other `-` is
    dropped; JSON.pm:229-233 writes an object carrying only `input` for a skipped record."""
    rows = data_lines(texts["default"])
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;USED_REF=-;GIVEN_REF=-")
    default = rebuild(texts["default"], rows)
    tab_header = header_lines(texts["tab"])[-1] + "\tUSED_REF\tGIVEN_REF"
    tab_rows = [r + ("\t-\t-") for r in data_lines(texts["tab"])]
    tab = "\n".join(header_lines(texts["tab"])[:-1] + [tab_header] + tab_rows) + "\n"
    json_lines = cvf.text_lines(texts["json"])
    obj = json.loads(json_lines[0])
    obj["intergenic_consequences"][0]["used_ref"] = "-"
    obj["intergenic_consequences"][0]["given_ref"] = "-"
    json_lines[0] = json.dumps(obj)
    skipped_line = "21\t99\tskipped_by_size\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=99999999"
    json_lines.append(json.dumps({"input": skipped_line}))
    vcf = rebuild(texts["vcf"], data_lines(texts["vcf"]) + [skipped_line])
    report = cvf.SelfParity(default, tab, vcf, "\n".join(json_lines) + "\n").run()
    assert report["verdict"] == "green", report["reasons"]
    assert {"USED_REF", "GIVEN_REF"} <= set(report["pairs"]["json"]["columns_compared"])
    assert report["records"] == RECORDS + 1
    del obj["intergenic_consequences"][0]["given_ref"]
    json_lines[0] = json.dumps(obj)
    report = cvf.SelfParity(default, tab, vcf, "\n".join(json_lines) + "\n").run()
    assert report["pairs"]["json"]["differing"]["GIVEN_REF"]["count"] == 1


def test_floors_masks_and_field_pairs_parse():
    assert cvf.parse_floors("HGVSc=30, IMPACT=977") == {"HGVSc": 30, "IMPACT": 977}
    assert cvf.parse_floors(None) == {}
    with pytest.raises(cvf.ComparatorError):
        cvf.parse_floors("HGVSc=many")
    with pytest.raises(cvf.ComparatorError):
        cvf.parse_floors("HGVSc")
    assert cvf.parse_masks("a, b,,c") == {"a", "b", "c"}
    assert cvf.parse_masks(None) == set()
    assert cvf.parse_field_pairs(["HGVSp=p.Arg97%3D=>p.Arg97="]) == {"HGVSp": {("p.Arg97%3D", "p.Arg97=")}}
    with pytest.raises(cvf.ComparatorError):
        cvf.parse_field_pairs(["HGVSp=a:b"])
