"""Tests for compare_vep_formats.py: the four parsers on the committed Ensembl VEP 116.2
outputs of tests/golden/116/GRCh38-release116, the comparison rules, self-parity and the
warnings mode. The expected values are read off the committed files; the constants below
are the records the reference annotated, its consequence rows and the records it skipped."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_vep_formats as cvf  # noqa: E402

GOLDEN = HERE.parent.parent / "tests" / "golden" / "116"
CORPUS = GOLDEN / "GRCh38-release116"
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


def test_json_objects_name_their_record_by_the_first_five_columns(tmp_path):
    """Parser/VCF.pm:464 stores a QUAL of ``0`` as ``.``, so the ``input`` of the committed
    GRCh38-vcf-edge-cases JSON object for ``sv_qual0`` (record 40, QUAL ``0``) is not the
    input line. The object still keys to its record, as a VCF line does, so a candidate
    output carrying the input-only object of the reference-skipped record 26 (the manifest's
    documented ``reference_skipped_record`` key) before it compares green under the class
    and shifts no later object onto its ordinal."""
    corpus = GOLDEN / "GRCh38-vcf-edge-cases"
    inputs = cvf.input_records(corpus / "variants.vcf")
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    text = cvf.read_text(corpus / "expected" / "json.jsonl.gz")
    objects = cvf.json_objects(text)
    assert inputs[40].split("\t")[2:6] == ["sv_qual0", "N", "<DEL>", "0"]
    echo = next(obj["input"] for _, obj in objects if obj["id"] == "sv_qual0")
    assert echo == inputs[40].replace("\t0\tPASS\t", "\t.\tPASS\t") != inputs[40]
    entries = cvf.parse_json(text, inputs)
    assert [e.key for e in entries if e.record == 40] == ["rec40|record", "rec40|deletion|ENST00000286808|transcript"]
    assert all(e.record is not None for e in entries)

    (extra,) = manifest["vep_rs_only_tuples"]
    assert extra["expected_divergence"] == "reference_skipped_record" and extra["record_indices"] == [26]
    identity = {cvf.JSON_IDENTITY_KEY: "0.5.0"}
    records = [e for e in entries if e.key.endswith("|record")]
    position = next(n for (n, _), e in zip(objects, records) if e.record > 26)
    lines = [json.dumps({**obj, **identity}) for _, obj in objects]
    lines.insert(position, json.dumps({"input": inputs[26], **identity}))
    expected_file, actual_file = tmp_path / "expected.jsonl", tmp_path / "actual.jsonl"
    expected_file.write_text(text, encoding="utf-8")
    actual_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    masked = {"reference_skipped_record"}
    report = cvf.compare_files(
        "json", expected_file, actual_file, inputs, cvf.Documented.from_manifest(manifest, masked), masked, {}, {}, 4
    )
    assert report["verdict"] == "green", report["reasons"]
    assert report["counts"]["expected_entries"] == 98
    assert report["counts"]["actual_entries"] == 99
    assert report["counts"]["compared"] == 98
    assert report["counts"]["record_entries_on_reference_skipped_records"] == 1
    assert report["one_sided"] == {"expected_only": 0, "actual_only": 0}
    unmasked = cvf.compare_files(
        "json", expected_file, actual_file, inputs, cvf.Documented.from_manifest(manifest, set()), set(), {}, {}, 4
    )
    assert unmasked["verdict"] == "red"
    assert unmasked["one_sided"] == {"expected_only": 0, "actual_only": 1}
    assert [x for group in unmasked["unexpected"].values() for x in group["examples"]] == ["rec26|record"]


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
    assert masked.available == {"breakend_mate_local_read": 2, "reference_skipped_record": 4}
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


def test_regulatory_feature_rows_mask_drops_the_reference_rows(texts, inputs, tmp_path):
    """Under ``regulatory_feature_rows`` the reference's RegulatoryFeature and MotifFeature
    rows go in every format and a JSON record's most_severe_consequence is recomputed from
    the rows that remain; without the mask the same rows are missing in the candidate, and
    on a reference without such rows the class is unused."""
    regulatory_row = (
        "bnd_intergenic_pair_a\t21:11000001\tN.\t-\tENSR00000347123\tRegulatoryFeature\tregulatory_region_variant"
        "\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter"
    )
    motif_row = (
        "bnd_intergenic_pair_a\t21:11000001\tN.\t-\tENSM00000577720\tMotifFeature\tTF_binding_site_variant"
        "\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;MOTIF_NAME=ENSPFM0086;MOTIF_POS=3;HIGH_INF_POS=N;TRANSCRIPTION_FACTORS=EGR1"
    )
    rows = data_lines(texts["default"])
    reference = cvf.parse_default(rebuild(texts["default"], rows[:1] + [regulatory_row, motif_row] + rows[1:]))
    actual = cvf.parse_default(texts["default"])
    assert cvf.compare_entries(reference, actual, no_documented())["one_sided"]["expected_only"] == 2
    kept = cvf.without_regulatory_rows(reference)
    assert [e.key for e in kept] == [e.key for e in actual]
    assert cvf.compare_entries(kept, actual, no_documented())["verdict"] == "green"

    vcf_rows = data_lines(texts["vcf"])
    csq_header = next(line for line in header_lines(texts["vcf"]) if line.startswith("##INFO=<ID=CSQ,"))
    csq_fields = cvf.csq_format_fields(csq_header)
    values = {"Allele": "N.", "Consequence": "regulatory_region_variant", "IMPACT": "MODIFIER",
              "Feature_type": "RegulatoryFeature", "Feature": "ENSR00000347123", "BIOTYPE": "promoter"}
    chunk = "|".join(values.get(name, "") for name in csq_fields)
    vcf_rows[0] = vcf_rows[0].replace(";CSQ=", f";CSQ={chunk},")
    vcf_reference = cvf.parse_vcf(rebuild(texts["vcf"], vcf_rows), inputs)
    vcf_actual = cvf.parse_vcf(texts["vcf"], inputs)
    assert cvf.compare_entries(vcf_reference, vcf_actual, no_documented())["one_sided"]["expected_only"] == 1
    assert cvf.compare_entries(cvf.without_regulatory_rows(vcf_reference), vcf_actual, no_documented())["verdict"] == "green"

    json_lines = [line for line in cvf.text_lines(texts["json"]) if line]
    first = json.loads(json_lines[0])
    first["regulatory_feature_consequences"] = [{
        "variant_allele": "N.", "regulatory_feature_id": "ENSR00000347123",
        "consequence_terms": ["regulatory_region_variant"], "impact": "MODIFIER", "biotype": "promoter",
    }]
    first["most_severe_consequence"] = "regulatory_region_variant"
    json_reference = cvf.parse_json("\n".join([json.dumps(first)] + json_lines[1:]) + "\n", inputs)
    json_actual = cvf.parse_json(texts["json"], inputs)
    unmasked = cvf.compare_entries(json_reference, json_actual, no_documented())
    assert unmasked["one_sided"]["expected_only"] == 1
    assert unmasked["columns"]["most_severe_consequence"]["differing"] == 1
    kept = cvf.without_regulatory_rows(json_reference)
    assert kept[0].fields["most_severe_consequence"] == "intergenic_variant"
    assert "regulatory_feature_consequences" not in kept[0].fields
    assert cvf.compare_entries(kept, json_actual, no_documented())["verdict"] == "green"
    # A record that lost no row keeps the reference's term, so a wrong one stays red.
    second = json.loads(json_lines[1])
    second["most_severe_consequence"] = "stop_gained"
    kept = cvf.without_regulatory_rows(cvf.parse_json("\n".join([json.dumps(first), json.dumps(second)] + json_lines[2:]) + "\n", inputs))
    assert cvf.compare_entries(kept, json_actual, no_documented())["columns"]["most_severe_consequence"]["differing"] == 1

    expected_warnings = tmp_path / "ref_warnings.txt"
    expected_warnings.write_text("WARNING: line 2 skipped (21 1 x A G . . .): Invalid allele string\n", encoding="utf-8")
    actual_warnings = tmp_path / "cand_warnings.txt"
    notice = f"WARNING: --everything implies --regulatory, which has no effect: {cvf.REGULATORY_NOTICE}\n"
    actual_warnings.write_text(notice + "WARNING: line 2 skipped (): Invalid allele string\n", encoding="utf-8")
    assert cvf.compare_warnings(expected_warnings, actual_warnings, None, None)["verdict"] == "red"
    masked_report = cvf.compare_warnings(expected_warnings, actual_warnings, None, None, {cvf.REGULATORY_FEATURE_ROWS})
    assert masked_report["verdict"] == "green" and masked_report["warnings"]["shared"] == 1
    actual_warnings.write_text("WARNING: line 2 skipped (): Invalid allele string\n", encoding="utf-8")
    missing = cvf.compare_warnings(expected_warnings, actual_warnings, None, None, {cvf.REGULATORY_FEATURE_ROWS})
    assert missing["verdict"] == "red" and missing["reasons"] == ["warnings: the regulatory notice is missing from the actual file"]

    expected_path = tmp_path / "reference.txt"
    expected_path.write_text(rebuild(texts["default"], rows[:1] + [regulatory_row] + rows[1:]), encoding="utf-8")
    actual_path = tmp_path / "actual.txt"
    actual_path.write_text(texts["default"], encoding="utf-8")
    report = cvf.compare_files(
        "default", expected_path, actual_path, None, no_documented(), {cvf.REGULATORY_FEATURE_ROWS}, {}, {}, 4
    )
    assert report["verdict"] == "green"
    assert report["masks"] == {
        "requested": ["regulatory_feature_rows"], "applied": {"regulatory_feature_rows": 1}, "unused": [], "available": {},
    }
    report = cvf.compare_files(
        "default", actual_path, actual_path, None, no_documented(), {cvf.REGULATORY_FEATURE_ROWS}, {}, {}, 4
    )
    assert report["verdict"] == "green"
    assert report["masks"]["applied"] == {} and report["masks"]["unused"] == ["regulatory_feature_rows"]


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
        "nearest_accepted": 0, "order_only": 0, "individual_allele_order": 0, "floor": ROWS, "floor_met": True,
    }
    assert met["floors_not_applicable"] == []
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


def test_unordered_compares_the_shared_entries_without_their_order(tmp_path, capsys):
    """The GRCh38-pick-summary manifest marks the three flag sets whose rows Ensembl VEP
    emits from a hash ``unordered`` (OutputFactory.pm:593 and :605 ``keys %by_allele``, :855
    ``keys %by_gene``). A candidate --per_gene output with the first record's two gene rows
    the other way round is red to the comparator and green under ``--unordered``, which
    turns off the order check alone: a dropped row stays red."""
    corpus = GOLDEN / "GRCh38-pick-summary"
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    unordered_sets = [f["name"] for f in manifest["flag_sets"] if f.get("unordered")]
    assert unordered_sets == ["pick-allele", "per-gene", "pick-allele-gene"]
    expected_file = corpus / "expected_per-gene" / "default.txt.gz"
    text = cvf.read_text(expected_file)
    rows = data_lines(text)
    assert [r.split("\t")[4] for r in rows[:2]] == ["ENST00000815440", "ENST00000624291"]
    swapped = tmp_path / "swapped.txt"
    swapped.write_text(rebuild(text, [rows[1], rows[0]] + rows[2:]), encoding="utf-8")
    report = tmp_path / "r.json"
    args = ["--format", "default", "--expected", str(expected_file), "--actual", str(swapped), "--report", str(report)]
    assert cvf.main(args) == 1
    assert capsys.readouterr().out == "default: red\n  order differs at shared position 0\n"
    written = json.loads(report.read_text(encoding="utf-8"))
    assert written["unordered"] is False
    assert written["order"]["expected"] == "21:8857289|C|ENST00000815440|Transcript"
    assert cvf.main(args + ["--unordered"]) == 0
    assert capsys.readouterr().out == "default: green\n"
    written = json.loads(report.read_text(encoding="utf-8"))
    assert written["unordered"] is True
    assert written["order"] is None
    assert written["counts"]["compared"] == len(rows) == 70
    dropped = tmp_path / "dropped.txt"
    dropped.write_text(rebuild(text, [rows[1]] + rows[2:]), encoding="utf-8")
    assert cvf.main(["--format", "default", "--expected", str(expected_file), "--actual", str(dropped), "--unordered"]) == 1
    assert capsys.readouterr().out == "default: red\n  missing in the candidate: 1 entries\n"


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


def synthetic_self_parity(extra: dict[str, str], consequence: dict) -> dict:
    """The self-parity report of one synthetic missense record: the default format's Extra
    carries IMPACT, STRAND and then ``extra``, the tab and VCF formats spell the same values, and
    the JSON transcript consequence carries ``consequence`` beside the fixed columns' members."""
    fixed = ("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
             "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation")
    values = ("rs1", "1:100", "G", "ENSG1", "ENST1", "Transcript", "missense_variant", "10", "5", "2", "A/G", "gCa/gGa", "rs1")
    extra_column = {"IMPACT": "MODERATE", "STRAND": "1", **extra}
    default = "#" + "\t".join(fixed + ("Extra",)) + "\n" + "\t".join(values + (";".join(f"{k}={v}" for k, v in extra_column.items()),)) + "\n"
    tab = "#" + "\t".join(fixed + tuple(extra_column)) + "\n" + "\t".join(values + tuple(extra_column.values())) + "\n"
    csq = ("Allele", "Consequence", "IMPACT", "Gene", "Feature_type", "Feature", "cDNA_position", "CDS_position",
           "Protein_position", "Amino_acids", "Codons", "Existing_variation", "STRAND", *extra)
    row = dict(zip(fixed, values)) | extra_column
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
            **consequence,
        }],
    }
    return cvf.SelfParity(default, tab, vcf, json.dumps(obj) + "\n").run()


def test_self_parity_undoes_the_json_sift_polyphen_and_domains_splits():
    """A synthetic record beside the committed --sift, --polyphen and --domains corpora: the
    JSON prediction/score pair and the db/name objects must match the default format's
    ``deleterious(0.01)`` and ``Pfam:PF1,PANTHER:PTHR1`` spellings."""
    extra = {"SIFT": "deleterious(0.01)", "PolyPhen": "probably_damaging(0.95)", "DOMAINS": "Pfam:PF1,PANTHER:PTHR1"}
    consequence = {
        "sift_prediction": "deleterious", "sift_score": 0.01,
        "polyphen_prediction": "probably_damaging", "polyphen_score": 0.95,
        "domains": [{"db": "Pfam", "name": "PF1"}, {"db": "PANTHER", "name": "PTHR1"}],
    }
    report = synthetic_self_parity(extra, consequence)
    assert report["verdict"] == "green", report["reasons"]
    assert {"SIFT", "PolyPhen", "DOMAINS"} <= set(report["pairs"]["json"]["columns_compared"])
    report = synthetic_self_parity(extra, {**consequence, "sift_score": 0.02, "domains": [{"db": "Pfam", "name": "PF1"}]})
    assert sorted(report["pairs"]["json"]["differing"]) == ["DOMAINS", "SIFT"]


def test_self_parity_domains_model_repeats_the_last_match_and_numbers_a_digit_name():
    """Ensembl VEP's JSON re-split (`(\\w+):(\\w+)`, OutputFactory/JSON.pm:363-371) keeps `$1` and
    `$2` from the list's last successful match, so an entry the pattern does not match repeats the
    previous object and an opening one is dropped; `numberify` then makes a digit-only name a
    number. `Gene3D:1.10.510.10,Low_complexity_(Seg):seg,Pfam:PF1` is therefore
    `{Gene3D 1}, {Gene3D 1}, {Pfam PF1}` and `Low_complexity_(Seg):seg,Pfam:PF1` is `{Pfam PF1}`.
    A capture of `0` fails `push ... if $1 && $2`, so `Gene3D:0.10.510.10,Low_complexity_(Seg):seg`
    is empty and `Pfam:PF1,Gene3D:0.10.510.10,Low_complexity_(Seg):seg` is `{Pfam PF1}` alone."""
    repeated = [{"db": "Gene3D", "name": 1}, {"db": "Gene3D", "name": 1}, {"db": "Pfam", "name": "PF1"}]
    report = synthetic_self_parity({"DOMAINS": "Gene3D:1.10.510.10,Low_complexity_(Seg):seg,Pfam:PF1"}, {"domains": repeated})
    assert report["verdict"] == "green", report["reasons"]
    report = synthetic_self_parity({"DOMAINS": "Low_complexity_(Seg):seg,Pfam:PF1"}, {"domains": [{"db": "Pfam", "name": "PF1"}]})
    assert report["verdict"] == "green", report["reasons"]
    report = synthetic_self_parity(
        {"DOMAINS": "Gene3D:1.10.510.10,Low_complexity_(Seg):seg,Pfam:PF1"},
        {"domains": [{"db": "Gene3D", "name": "1"}, {"db": "Pfam", "name": "PF1"}]},
    )
    assert sorted(report["pairs"]["json"]["differing"]) == ["DOMAINS"]
    report = synthetic_self_parity({"DOMAINS": "Gene3D:0.10.510.10,Low_complexity_(Seg):seg"}, {"domains": []})
    assert report["verdict"] == "green", report["reasons"]
    report = synthetic_self_parity(
        {"DOMAINS": "Pfam:PF1,Gene3D:0.10.510.10,Low_complexity_(Seg):seg"}, {"domains": [{"db": "Pfam", "name": "PF1"}]}
    )
    assert report["verdict"] == "green", report["reasons"]
    report = synthetic_self_parity(
        {"DOMAINS": "Gene3D:0.10.510.10,Low_complexity_(Seg):seg"},
        {"domains": [{"db": "Gene3D", "name": 0}, {"db": "Gene3D", "name": 0}]},
    )
    assert sorted(report["pairs"]["json"]["differing"]) == ["DOMAINS"]
    assert cvf.perl_numeric_or_text("007") == 7 and cvf.perl_numeric_or_text("PF1") == "PF1"
    # Python's str.isdigit accepts a superscript digit that Perl's looks_like_number refuses.
    assert cvf.perl_numeric_or_text("\u00b2") == "\u00b2"


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


NEAREST_CORPUS = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-identifiers"


def test_nearest_sites_agree_on_the_winning_sites_not_the_ids():
    """Ensembl VEP names one transcript at a winning site; the candidate names every transcript
    there. Transcripts A and B of gene G1 share the site 1000; R of G2 starts at 5000; F of G3
    at 9000."""
    sites = cvf.NearestSites({"A": {1000}, "B": {1000}, "R": {5000}, "F": {9000},
                              "G1": {1000}, "G2": {5000}, "G3": {9000}, "SYM1": {1000}, "SYM3": {9000}})
    assert sites.agree("A", "A,B")
    assert sites.agree("B", "A,B")
    assert sites.agree("A", "A&B")
    assert sites.agree("R,F", "R,F")
    assert sites.agree("G1", "G1")
    assert not sites.agree("A", "B"), "the reference's transcript must be among the candidate's"
    assert not sites.agree("A", "A,R"), "a candidate transcript at another site is a difference"
    assert not sites.agree("R", "R,F"), "the candidate adding the other side's site is a difference"
    assert not sites.agree("A", "A,X"), "a transcript the cache lacks starts nowhere"
    assert not sites.agree("X", "X"), "a transcript the cache lacks never agrees, even with itself"
    assert sites.agree("-", "-") and sites.agree(cvf.ABSENT, cvf.ABSENT)
    assert not sites.agree("A", "-") and not sites.agree("-", "A")
    assert not sites.agree("-", cvf.ABSENT), "an empty value and a missing field stay a difference"
    assert cvf.nearest_list("A&B") == ["A", "B"] and cvf.nearest_list("-") == []
    assert cvf.extra_nearest_accepted("IMPACT=HIGH;NEAREST=A", "IMPACT=HIGH;NEAREST=A,B", sites)
    assert not cvf.extra_nearest_accepted("IMPACT=HIGH;NEAREST=A", "NEAREST=A,B;IMPACT=HIGH", sites)
    assert not cvf.extra_nearest_accepted("IMPACT=HIGH;NEAREST=A", "IMPACT=LOW;NEAREST=A,B", sites)
    assert not cvf.extra_nearest_accepted("IMPACT=HIGH", "IMPACT=HIGH", sites)


def test_nearest_sites_from_cache_settle_the_identifiers_corpus():
    """The nearest-transcript set of the identifiers corpus: every record where Ensembl VEP's
    one transcript is one of several at the winning site differs by id and agrees by site;
    the gene and symbol sets agree verbatim."""
    sites = cvf.NearestSites.from_cache(NEAREST_CORPUS / "json_cache")
    assert sites.sites["ENST00000389125"] == {29940022}
    assert sites.sites["ENST00000399914"] == {29940022}
    assert sites.sites["GRIK1"] >= {29940022} and sites.sites["ENSG00000171189"] >= {29940022}
    assert sites.agree("ENST00000399914", "ENST00000389124,ENST00000389125,ENST00000399907,ENST00000399909,"
                                           "ENST00000399913,ENST00000399914,ENST00001131791")
    text = cvf.read_text(NEAREST_CORPUS / "expected_nearest-transcript" / "default.txt.gz")
    expected = cvf.parse_default(text)
    rows = data_lines(text)
    widened = [r.replace("NEAREST=ENST00000399914", "NEAREST=ENST00000389125,ENST00000399914") for r in rows]
    changed = sum(1 for r in rows if "NEAREST=ENST00000399914" in r)
    assert changed > 0
    actual = cvf.parse_default(rebuild(text, widened))
    plain = cvf.compare_entries(expected, actual, no_documented())
    assert plain["verdict"] == "red" and plain["columns"]["NEAREST"]["differing"] == changed
    settled = cvf.compare_entries(expected, actual, no_documented(), nearest_sites=sites)
    assert settled["verdict"] == "green", settled["reasons"]
    assert settled["columns"]["NEAREST"]["nearest_accepted"] == changed
    assert settled["columns"]["Extra"]["nearest_accepted"] == changed
    wrong = [r.replace("NEAREST=ENST00000399914", "NEAREST=ENST00000400546") for r in rows]
    report = cvf.compare_entries(expected, cvf.parse_default(rebuild(text, wrong)), no_documented(), nearest_sites=sites)
    assert report["verdict"] == "red" and report["columns"]["NEAREST"]["differing"] == changed
    with pytest.raises(cvf.ComparatorError):
        cvf.NearestSites.from_cache(NEAREST_CORPUS / "expected")


def one_shard_cache(root: Path, shard: object) -> Path:
    (root / "transcripts" / "1").mkdir(parents=True)
    (root / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(shard), encoding="utf-8")
    return root


def test_nearest_sites_from_cache_rejects_a_shard_of_another_shape(tmp_path, capsys):
    """A shard that parses as neither a transcript array nor a range-keyed map of arrays, or
    that holds an element other than an object or ``null``, is a read failure: the comparator's
    error line naming the shard and exit 2, not a traceback."""
    other_shapes = ("x", 7, None, [7], ["x"], [[]], {"1-1000000": "x"}, {"1-1000000": {"a": 1}})
    for n, shard in enumerate(other_shapes):
        with pytest.raises(cvf.ComparatorError, match=r"1-1000000\.json"):
            cvf.NearestSites.from_cache(one_shard_cache(tmp_path / f"shape{n}", shard))
    assert cvf.NearestSites.from_cache(one_shard_cache(tmp_path / "empty_array", [])).sites == {}
    assert cvf.NearestSites.from_cache(one_shard_cache(tmp_path / "empty_map", {})).sites == {}
    code = cvf.main([
        "--format", "default", "--expected", str(EXPECTED / "default.txt.gz"), "--actual", str(EXPECTED / "default.txt.gz"),
        "--nearest-sites", str(one_shard_cache(tmp_path / "cli", "x")),
    ])
    assert code == 2
    err = capsys.readouterr().err
    assert err.startswith("ERROR: [compare_vep_formats] ")
    assert "1-1000000.json" in err
    assert err.count("\n") == 1


def test_nearest_sites_from_cache_drops_nulls_and_rejects_a_coordinate_that_is_not_an_integer(tmp_path):
    """A ``null`` is a placeholder in either shard shape and is dropped. A protein-coding
    transcript's site is its start on the forward strand and its end on the reverse, an integer
    or a string of digits; one that is missing or of another kind is a read failure, while
    another biotype's coordinates are never read."""
    plus = {"biotype": "protein_coding", "strand": 1, "start": 10, "end": 20, "stable_id": "P", "gene_stable_id": "G", "gene_symbol": "S"}
    minus = {"biotype": "protein_coding", "strand": -1, "start": "30", "end": "40", "stable_id": "M", "gene_stable_id": "G"}
    other = {"biotype": "lncRNA", "strand": 1, "stable_id": "L"}
    sites = {"P": {10}, "M": {40}, "G": {10, 40}, "S": {10}}
    assert cvf.NearestSites.from_cache(one_shard_cache(tmp_path / "array", [None, plus, None, minus, other])).sites == sites
    wrapped = {"1-1000000": [None, plus], "1000001-2000000": [minus, None, other]}
    assert cvf.NearestSites.from_cache(one_shard_cache(tmp_path / "map", wrapped)).sites == sites
    without_an_integer_site = (
        {"biotype": "protein_coding", "strand": 1, "end": 20, "stable_id": "P"},
        {"biotype": "protein_coding", "strand": -1, "start": 10, "stable_id": "P"},
        {"biotype": "protein_coding", "strand": 1, "start": "abc", "end": 20, "stable_id": "P"},
        {"biotype": "protein_coding", "strand": 1, "start": 12.5, "end": 20, "stable_id": "P"},
        {"biotype": "protein_coding", "strand": 1, "start": True, "end": 20, "stable_id": "P"},
        {"biotype": "protein_coding", "strand": 1, "start": None, "end": 20, "stable_id": "P"},
    )
    for n, transcript in enumerate(without_an_integer_site):
        with pytest.raises(cvf.ComparatorError, match=r"1-1000000\.json: transcript 'P' has no integer (start|end)$"):
            cvf.NearestSites.from_cache(one_shard_cache(tmp_path / f"bad{n}", [transcript]))


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


def test_self_parity_keys_json_objects_by_the_first_five_columns_of_the_input(texts, inputs):
    """A structural variant whose QUAL is ``0`` is echoed with ``.`` in the JSON ``input`` and
    written back with ``.`` in the VCF output (Parser/VCF.pm:464), so under ``--input`` the
    object and the VCF line of record 11 name the same record only through their first five
    columns; the committed formats, read with such an input, stay in parity."""
    qual0 = list(inputs)
    assert qual0[11].split("\t")[2:6] == ["multi_ins_ins", "N", "<INS>,<INS>", "."]
    qual0[11] = qual0[11].replace("\t<INS>,<INS>\t.\t", "\t<INS>,<INS>\t0\t", 1)
    echo = next(json.loads(line)["input"] for line in cvf.text_lines(texts["json"]) if '"multi_ins_ins"' in line)
    assert echo == inputs[11] != qual0[11]
    report = cvf.SelfParity(texts["default"], texts["tab"], texts["vcf"], texts["json"], inputs=qual0).run()
    assert report["verdict"] == "green", report["reasons"]
    assert report["records"] == RECORDS


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


PROTEIN_FEATURES = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-protein-features"
FREQFILTER = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-freqfilter"
# Rows of the protein-features reference whose DOMAINS carries two or more items that differ,
# so reversing the items changes the text: 118 of the 119 multi-item rows in the default
# format, and 116 JSON consequence objects once the JSON re-split has mapped the items.
REORDERABLE_DOMAINS_ROWS = 118
REORDERABLE_DOMAINS_OBJECTS = 116


def reversed_lists(text: str, column: str) -> str:
    """``text`` with the items of every ``column=`` list of the Extra text reversed."""

    def flip(m):
        return column + "=" + ",".join(reversed(m.group(1).split(",")))

    return re.sub(rf"{column}=([^;\t\n]+)", flip, text)


def test_order_only_class_accepts_domains_in_any_order_in_every_format():
    """golden-GRCh38-protein-features under ``--masks order_only``: a candidate printing the
    protein features of a row in another order (a native cache orders them by feature id, the
    Ensembl cache by its array) agrees as a multiset, counted per column as ``order_only``;
    without the class the same rows are differing fields, and a reorder of any other column,
    or a repeated item, stays red under it."""
    inputs = cvf.input_records(PROTEIN_FEATURES / "variants.vcf")
    default = cvf.read_text(PROTEIN_FEATURES / "expected" / "default.txt.gz")
    expected = cvf.parse_default(default)
    actual = cvf.parse_default(reversed_lists(default, "DOMAINS"))
    plain = cvf.compare_entries(expected, actual, no_documented())
    assert plain["verdict"] == "red"
    assert plain["columns"]["DOMAINS"]["differing"] == REORDERABLE_DOMAINS_ROWS
    assert plain["columns"]["Extra"]["differing"] == REORDERABLE_DOMAINS_ROWS
    documented = no_documented()
    accepted = cvf.compare_entries(expected, actual, documented, order_only=cvf.order_only_fields("default"))
    assert accepted["verdict"] == "green"
    assert accepted["columns"]["DOMAINS"]["differing"] == 0
    assert accepted["columns"]["DOMAINS"]["order_only"] == REORDERABLE_DOMAINS_ROWS
    assert accepted["columns"]["Extra"]["order_only"] == REORDERABLE_DOMAINS_ROWS
    assert documented.applied == {cvf.ORDER_ONLY: 2 * REORDERABLE_DOMAINS_ROWS}
    example = accepted["columns"]["DOMAINS"]["order_only_examples"][0]
    assert example.startswith("21:14216193|G|ENST00000400577|Transcript: expected 'PDB-ENSP_mappings:2ywk.A,")
    assert "got 'AFDB-ENSP_mappings:AF-P57052-F1,PDB-ENSP_mappings:2ywk.A'" in example

    # The tab format's DOMAINS column, the same rule.
    tab = cvf.read_text(PROTEIN_FEATURES / "expected" / "tab.txt.gz")
    domains_column = cvf.header_columns(tab).index("DOMAINS")
    tab_rows = []
    for row in data_lines(tab):
        cells = row.split("\t")
        cells[domains_column] = ",".join(reversed(cells[domains_column].split(",")))
        tab_rows.append("\t".join(cells))
    tab_report = cvf.compare_entries(
        cvf.parse_tab(tab), cvf.parse_tab(rebuild(tab, tab_rows)), no_documented(),
        order_only=cvf.order_only_fields("tab"),
    )
    assert tab_report["verdict"] == "green"
    assert tab_report["columns"]["DOMAINS"]["order_only"] == REORDERABLE_DOMAINS_ROWS

    # The VCF format joins the items with & (VCF.pm:391).
    vcf = cvf.read_text(PROTEIN_FEATURES / "expected" / "vcf.vcf.gz")
    vcf_rows = data_lines(vcf)
    flipped = []
    for row in vcf_rows:
        prefix, marker, csq = row.partition(";CSQ=")
        if not csq:
            prefix, marker, csq = row.partition("CSQ=")
        chunks = []
        for chunk in csq.split(","):
            vals = chunk.split("|")
            vals[-1] = "&".join(reversed(vals[-1].split("&")))
            chunks.append("|".join(vals))
        flipped.append(prefix + marker + ",".join(chunks))
    vcf_expected = cvf.parse_vcf(vcf, inputs)
    vcf_actual = cvf.parse_vcf(rebuild(vcf, flipped), inputs)
    assert cvf.compare_entries(vcf_expected, vcf_actual, no_documented())["columns"]["DOMAINS"]["differing"] == REORDERABLE_DOMAINS_ROWS
    vcf_report = cvf.compare_entries(vcf_expected, vcf_actual, no_documented(), order_only=cvf.order_only_fields("vcf"))
    assert vcf_report["verdict"] == "green"
    assert vcf_report["columns"]["DOMAINS"]["order_only"] == REORDERABLE_DOMAINS_ROWS

    # The JSON format's domains array (JSON.pm:363-371), the objects compared as a multiset.
    json_text = cvf.read_text(PROTEIN_FEATURES / "expected" / "json.jsonl.gz")
    objects = []
    for line in cvf.text_lines(json_text):
        if not line:
            continue
        obj = json.loads(line)
        for tc in obj.get("transcript_consequences", []):
            if "domains" in tc:
                tc["domains"] = list(reversed(tc["domains"]))
        objects.append(json.dumps(obj))
    json_expected = cvf.parse_json(json_text, inputs)
    json_actual = cvf.parse_json("\n".join(objects) + "\n", inputs)
    assert cvf.compare_entries(json_expected, json_actual, no_documented())["columns"]["domains"]["differing"] == REORDERABLE_DOMAINS_OBJECTS
    json_report = cvf.compare_entries(json_expected, json_actual, no_documented(), order_only=cvf.order_only_fields("json"))
    assert json_report["verdict"] == "green"
    assert json_report["columns"]["domains"]["order_only"] == REORDERABLE_DOMAINS_OBJECTS
    assert cvf.order_only_fields("json") == {"domains": cvf.JSON_ARRAY, "zyg": ","}

    # Another column's order, or a repeated item, is not an order-only agreement.
    rows = data_lines(default)
    target = next(i for i, r in enumerate(rows) if "DOMAINS=PDB-ENSP_mappings:2ywk.A,AFDB-ENSP_mappings:AF-P57052-F1" in r)
    rows[target] = rows[target].replace(
        "DOMAINS=PDB-ENSP_mappings:2ywk.A,AFDB-ENSP_mappings:AF-P57052-F1",
        "DOMAINS=PDB-ENSP_mappings:2ywk.A,PDB-ENSP_mappings:2ywk.A,AFDB-ENSP_mappings:AF-P57052-F1",
    )
    repeated = cvf.compare_entries(expected, cvf.parse_default(rebuild(default, rows)), no_documented(), order_only=cvf.order_only_fields("default"))
    assert repeated["verdict"] == "red"
    assert repeated["columns"]["DOMAINS"]["differing"] == 1
    rows = data_lines(default)
    rows[target] = rows[target].replace("\tmissense_variant\t", "\tsplice_region_variant,missense_variant\t", 1)
    other = cvf.compare_entries(expected, cvf.parse_default(rebuild(default, rows)), no_documented(), order_only=cvf.order_only_fields("default"))
    assert other["verdict"] == "red"
    assert other["columns"]["Consequence"]["differing"] == 1
    assert other["columns"]["Consequence"]["order_only"] == 0


def test_order_only_class_covers_freqs_and_is_reported_through_compare_files(tmp_path):
    """FREQS is printed in hash order (OutputFactory.pm:1152-1153), so a two-allele record's
    entries agree as a multiset under the class; ``compare_files`` applies the class from
    ``--masks`` and reports it applied or unused, with the other formats untouched."""
    default = cvf.read_text(FREQFILTER / "expected_freqfilter" / "default.txt.gz")
    assert "FREQS=AF:C:0.0007987" in default
    reference = default.replace("FREQS=AF:C:0.0007987", "FREQS=AF:C:0.0007987,AF:T:0.5")
    candidate = default.replace("FREQS=AF:C:0.0007987", "FREQS=AF:T:0.5,AF:C:0.0007987")
    reference_path = tmp_path / "reference.txt"
    reference_path.write_text(reference, encoding="utf-8")
    candidate_path = tmp_path / "candidate.txt"
    candidate_path.write_text(candidate, encoding="utf-8")
    rows = default.count("FREQS=AF:C:0.0007987")

    plain = cvf.compare_files("default", reference_path, candidate_path, None, no_documented(), set(), {}, {}, 4)
    assert plain["verdict"] == "red"
    assert plain["columns"]["FREQS"]["differing"] == rows
    documented = no_documented()
    report = cvf.compare_files("default", reference_path, candidate_path, None, documented, {cvf.ORDER_ONLY}, {}, {}, 4)
    assert report["verdict"] == "green"
    assert report["columns"]["FREQS"]["order_only"] == rows
    assert report["columns"]["Extra"]["order_only"] == rows
    assert report["masks"] == {
        "requested": [cvf.ORDER_ONLY], "applied": {cvf.ORDER_ONLY: 2 * rows}, "unused": [], "available": {},
    }
    # The class is unused on a candidate that prints every list in the reference's order.
    same = cvf.compare_files("default", reference_path, reference_path, None, no_documented(), {cvf.ORDER_ONLY}, {}, {}, 4)
    assert same["verdict"] == "green"
    assert same["masks"]["applied"] == {} and same["masks"]["unused"] == [cvf.ORDER_ONLY]
    assert same["columns"]["FREQS"]["order_only"] == 0
    # The record-level FREQS has no JSON key, so the JSON format compares no such field.
    assert "FREQS" in cvf.RECORD_LEVEL_COLUMNS
    assert cvf.same_items("a,b", "b,a", ",") and not cvf.same_items("a,b", "a,a,b", ",")
    assert cvf.same_items("-", "-", ",") and not cvf.same_items("-", "a", ",")
    assert cvf.same_items('{"db":"x","name":"1"},{"db":"y","name":"2"}', '{"db":"y","name":"2"},{"db":"x","name":"1"}', cvf.JSON_ARRAY)
    assert not cvf.same_items("a,b", "b,a", cvf.JSON_ARRAY)


INDIVIDUAL_ZYG = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-individual-zyg"
# Rows of the individual-zyg reference whose ZYG list names two or more samples, so reversing
# the list changes the text: 1,530 of the 1,534 rows carrying ZYG in the default format, and
# 103 of the 104 record objects carrying zyg in the JSON format.
REORDERABLE_ZYG_ROWS = 1530
REORDERABLE_ZYG_OBJECTS = 103


def test_order_only_class_covers_the_zyg_list_in_every_format():
    """``--individual_zyg`` prints one ``sample:ZYG`` entry per named sample in hash order
    (OutputFactory.pm:983-987), so a candidate naming the samples in another order agrees as
    a multiset under the class in every format: the default and tab ``ZYG`` column with the
    ``Extra`` text carrying it, the ``&``-joined CSQ field, and the record object's ``zyg``
    array, which the JSON format carries once at the top level, so ``ZYG`` is a record-level
    column that self-parity sets aside against the consequence entries and that has no JSON
    floor; a changed zygosity or a sample fewer stays red."""
    expected_dir = INDIVIDUAL_ZYG / "expected_individual-zyg"
    inputs = cvf.input_records(INDIVIDUAL_ZYG / "variants.vcf")
    default = cvf.read_text(expected_dir / "default.txt.gz")
    expected = cvf.parse_default(default)
    actual = cvf.parse_default(reversed_lists(default, "ZYG"))
    plain = cvf.compare_entries(expected, actual, no_documented())
    assert plain["verdict"] == "red"
    assert plain["columns"]["ZYG"]["differing"] == REORDERABLE_ZYG_ROWS
    assert plain["columns"]["Extra"]["differing"] == REORDERABLE_ZYG_ROWS
    documented = no_documented()
    accepted = cvf.compare_entries(expected, actual, documented, order_only=cvf.order_only_fields("default"))
    assert accepted["verdict"] == "green"
    assert accepted["columns"]["ZYG"]["differing"] == 0
    assert accepted["columns"]["ZYG"]["order_only"] == REORDERABLE_ZYG_ROWS
    assert accepted["columns"]["Extra"]["order_only"] == REORDERABLE_ZYG_ROWS
    assert documented.applied == {cvf.ORDER_ONLY: 2 * REORDERABLE_ZYG_ROWS}
    example = accepted["columns"]["ZYG"]["order_only_examples"][0]
    assert example.startswith("21:9614117-9614200|-|ENST00000842019|Transcript: expected 'S7:HOMREF,S4:HOMREF,")
    assert "got 'S1:HOMREF,S5:HOMREF,S6:HOMREF,S8:HOMREF,S3:HOMREF,S2:HOMREF,S4:HOMREF,S7:HOMREF'" in example

    # The tab format's ZYG column, the same rule.
    tab = cvf.read_text(expected_dir / "tab.txt.gz")
    zyg_column = cvf.header_columns(tab).index("ZYG")
    tab_rows = []
    for row in data_lines(tab):
        cells = row.split("\t")
        cells[zyg_column] = ",".join(reversed(cells[zyg_column].split(",")))
        tab_rows.append("\t".join(cells))
    tab_report = cvf.compare_entries(
        cvf.parse_tab(tab), cvf.parse_tab(rebuild(tab, tab_rows)), no_documented(),
        order_only=cvf.order_only_fields("tab"),
    )
    assert tab_report["verdict"] == "green"
    assert tab_report["columns"]["ZYG"]["order_only"] == REORDERABLE_ZYG_ROWS

    # The VCF format joins the items with & (VCF.pm:391); the sample columns follow INFO.
    vcf = cvf.read_text(expected_dir / "vcf.vcf.gz")
    csq_header = next(line for line in cvf.text_lines(vcf) if line.startswith("##INFO=<ID=CSQ"))
    zyg_field = cvf.csq_format_fields(csq_header).index("ZYG")
    flipped = []
    for row in data_lines(vcf):
        cols = row.split("\t")
        info = []
        for kv in cols[7].split(";"):
            if kv.startswith("CSQ="):
                chunks = []
                for chunk in kv[4:].split(","):
                    vals = chunk.split("|")
                    vals[zyg_field] = "&".join(reversed(vals[zyg_field].split("&")))
                    chunks.append("|".join(vals))
                kv = "CSQ=" + ",".join(chunks)
            info.append(kv)
        cols[7] = ";".join(info)
        flipped.append("\t".join(cols))
    vcf_expected = cvf.parse_vcf(vcf, inputs)
    vcf_actual = cvf.parse_vcf(rebuild(vcf, flipped), inputs)
    assert cvf.compare_entries(vcf_expected, vcf_actual, no_documented())["columns"]["ZYG"]["differing"] == REORDERABLE_ZYG_ROWS
    vcf_report = cvf.compare_entries(vcf_expected, vcf_actual, no_documented(), order_only=cvf.order_only_fields("vcf"))
    assert vcf_report["verdict"] == "green"
    assert vcf_report["columns"]["ZYG"]["order_only"] == REORDERABLE_ZYG_ROWS

    # The JSON format's zyg array, once per record object, its strings compared as a multiset.
    json_text = cvf.read_text(expected_dir / "json.jsonl.gz")
    objects = []
    for line in cvf.text_lines(json_text):
        if not line:
            continue
        obj = json.loads(line)
        if "zyg" in obj:
            obj["zyg"] = list(reversed(obj["zyg"]))
        objects.append(json.dumps(obj))
    json_expected = cvf.parse_json(json_text, inputs)
    json_actual = cvf.parse_json("\n".join(objects) + "\n", inputs)
    assert cvf.compare_entries(json_expected, json_actual, no_documented())["columns"]["zyg"]["differing"] == REORDERABLE_ZYG_OBJECTS
    json_report = cvf.compare_entries(json_expected, json_actual, no_documented(), order_only=cvf.order_only_fields("json"))
    assert json_report["verdict"] == "green"
    assert json_report["columns"]["zyg"]["order_only"] == REORDERABLE_ZYG_OBJECTS

    # A changed zygosity, or a sample fewer, is not an order-only agreement.
    full_list = "ZYG=S7:HOMREF,S4:HOMREF,S2:HOMREF,S3:HOMREF,S8:HOMREF,S6:HOMREF,S5:HOMREF,S1:HOMREF"
    rows = data_lines(default)
    target = next(i for i, r in enumerate(rows) if full_list in r)
    rows[target] = rows[target].replace(full_list, full_list.replace("S1:HOMREF", "S1:HET"), 1)
    changed = cvf.compare_entries(expected, cvf.parse_default(rebuild(default, rows)), no_documented(), order_only=cvf.order_only_fields("default"))
    assert changed["verdict"] == "red"
    assert changed["columns"]["ZYG"]["differing"] == 1
    assert changed["columns"]["ZYG"]["order_only"] == 0
    rows = data_lines(default)
    rows[target] = rows[target].replace(full_list, full_list.replace(",S1:HOMREF", ""), 1)
    fewer = cvf.compare_entries(expected, cvf.parse_default(rebuild(default, rows)), no_documented(), order_only=cvf.order_only_fields("default"))
    assert fewer["verdict"] == "red"
    assert fewer["columns"]["ZYG"]["differing"] == 1

    # ZYG is a record-level column: no JSON floor, and the self-parity of one engine's formats
    # checks it once per record, the rows' list against the object's zyg array as multisets,
    # instead of against the JSON consequence entries, which carry no zyg key.
    assert "ZYG" in cvf.RECORD_LEVEL_COLUMNS
    assert cvf.json_floor_fields("ZYG") is None
    report = cvf.SelfParity(default, tab, vcf, json_text, inputs=inputs).run()
    assert report["verdict"] == "green", report["reasons"]
    assert "ZYG" not in report["not_compared"]["default_record_level"]
    assert "ZYG" in report["pairs"]["vcf"]["columns_compared"]
    assert "ZYG" in report["pairs"]["tab"]["columns_compared"]
    reordered = cvf.SelfParity(default, tab, vcf, "\n".join(objects) + "\n", inputs=inputs).run()
    assert reordered["verdict"] == "green", reordered["reasons"]
    changed = "\n".join(objects).replace('"S1:HOMREF"', '"S1:HET"', 1) + "\n"
    report = cvf.SelfParity(default, tab, vcf, changed, inputs=inputs).run()
    assert report["verdict"] == "red"
    assert len(report["record_checks"]["problems"]) == 1
    assert report["record_checks"]["problems"][0].endswith("against the objects' zyg ['S1:HET,S2:HOMREF,S3:HOMREF,S4:HOMREF,S5:HOMREF,S6:HOMREF,S7:HOMREF,S8:HOMREF']")


def test_self_parity_holds_zyg_to_the_rows_only_for_objects_with_a_consequence():
    """Under ``--individual`` a reference-homozygous individual's object carries its ``zyg``
    and no consequence, and the default format writes no row for it (OutputFactory.pm:511),
    so the record check holds only the objects that carry a consequence to the rows."""
    individual = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-individual"
    expected_dir = individual / "expected_individual-all"
    report = cvf.SelfParity(
        cvf.read_text(expected_dir / "default.txt.gz"), cvf.read_text(expected_dir / "tab.txt.gz"),
        cvf.read_text(expected_dir / "vcf.vcf.gz"), cvf.read_text(expected_dir / "json.jsonl.gz"),
        inputs=cvf.input_records(individual / "variants.vcf"),
    ).run()
    assert report["record_checks"]["problems"] == []


INDIVIDUAL = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-individual"
# The hand-built hb_multi record (G, ALT A,T) and its sample S2, genotype 1/2: the copy's allele
# string is G/A/T and each of its four transcripts has one A row and one T row, in that order.
TWO_ALLELE_ROWS = 8


def test_individual_allele_order_class_accepts_a_copys_alleles_in_any_order(tmp_path):
    """Under ``--individual`` Ensembl VEP lists the alternate alleles of an individual's copy in
    hash order (Parser/VCF.pm:639), so a candidate listing hb_multi S2's T before its A agrees
    under the class in every format: each transcript's two rows in the other order, the JSON
    record's ``allele_string``, and under ``--allele_number --uploaded_allele`` the rows'
    ``ALLELE_NUM`` (one permutation of the alternates' numbers over every row of the copy) and
    ``UPLOADED_ALLELE``; each agreement is counted as ``individual_allele_order`` and the
    reordered entries under ``counts``. The same rows without ``IND`` stay red, as do a changed
    allele, a numbering that is not a permutation, a renumbering of one transcript alone, a
    reference allele numbered other than 0 and a reference row moved after the alternates;
    ``order_only`` alone accepts none of it."""

    def swapped_pairs(items: list) -> list:
        return [x for i in range(0, len(items), 2) for x in (items[i + 1], items[i])]

    def s2_rows(rows: list[str]) -> list[int]:
        return [i for i, r in enumerate(rows) if r.startswith("hb_multi\t") and "IND=S2;" in r]

    def with_class() -> cvf.Documented:
        return cvf.Documented.from_manifest(None, {cvf.INDIVIDUAL_ALLELE_ORDER})

    expected_dir = INDIVIDUAL / "expected_individual-all"
    inputs = cvf.input_records(INDIVIDUAL / "variants.vcf")
    default = cvf.read_text(expected_dir / "default.txt.gz")
    rows = data_lines(default)
    targets = s2_rows(rows)
    assert len(targets) == TWO_ALLELE_ROWS
    assert [rows[i].split("\t")[2] for i in targets] == ["A", "T"] * 4
    swapped = list(rows)
    for first, second in zip(targets[0::2], targets[1::2]):
        swapped[first], swapped[second] = rows[second], rows[first]
    expected = cvf.parse_default(default)
    actual = cvf.parse_default(rebuild(default, swapped))
    plain = cvf.compare_entries(expected, actual, no_documented())
    assert plain["verdict"] == "red"
    assert plain["reasons"] == [f"order differs at shared position {targets[0]}"]
    assert plain["order"]["expected"] == "21:34716340|A|ENST00000349499|Transcript#2"
    assert plain["order"]["actual"] == "21:34716340|T|ENST00000349499|Transcript"
    order_only = cvf.compare_entries(expected, actual, no_documented(), order_only=cvf.order_only_fields("default"))
    assert order_only["verdict"] == "red"
    documented = with_class()
    accepted = cvf.compare_entries(expected, actual, documented)
    assert accepted["verdict"] == "green"
    assert accepted["order"] is None
    assert accepted["counts"]["individual_entries_reordered"] == TWO_ALLELE_ROWS
    assert documented.applied == {cvf.INDIVIDUAL_ALLELE_ORDER: TWO_ALLELE_ROWS}
    assert all(col["individual_allele_order"] == 0 for col in accepted["columns"].values())

    # The same eight rows without IND are a record of no individual: the order stays red.
    expected_eight = cvf.parse_default(rebuild(default, [rows[i].replace("IND=S2;", "") for i in targets]))
    actual_eight = cvf.parse_default(rebuild(default, [swapped[i].replace("IND=S2;", "") for i in targets]))
    documented = with_class()
    scoped = cvf.compare_entries(expected_eight, actual_eight, documented)
    assert scoped["verdict"] == "red"
    assert scoped["reasons"] == ["order differs at shared position 0"]
    assert scoped["counts"]["individual_entries_reordered"] == 0
    assert documented.applied == {}

    # Through compare_files the class is reported applied, and unused on a candidate in the
    # reference's order.
    reference_path = tmp_path / "reference.txt"
    reference_path.write_text(default, encoding="utf-8")
    candidate_path = tmp_path / "candidate.txt"
    candidate_path.write_text(rebuild(default, swapped), encoding="utf-8")
    report = cvf.compare_files("default", reference_path, candidate_path, None, no_documented(), {cvf.INDIVIDUAL_ALLELE_ORDER}, {}, {}, 4)
    assert report["verdict"] == "green"
    assert report["masks"] == {
        "requested": [cvf.INDIVIDUAL_ALLELE_ORDER], "applied": {cvf.INDIVIDUAL_ALLELE_ORDER: TWO_ALLELE_ROWS}, "unused": [], "available": {},
    }
    same = cvf.compare_files("default", reference_path, reference_path, None, no_documented(), {cvf.INDIVIDUAL_ALLELE_ORDER}, {}, {}, 4)
    assert same["verdict"] == "green"
    assert same["masks"]["applied"] == {} and same["masks"]["unused"] == [cvf.INDIVIDUAL_ALLELE_ORDER]

    # The JSON record object's allele_string (JSON.pm:241) and its consequence entries.
    json_text = cvf.read_text(expected_dir / "json.jsonl.gz")
    objects = []
    for line in cvf.text_lines(json_text):
        if not line:
            continue
        obj = json.loads(line)
        if obj.get("id") == "hb_multi" and obj.get("ind") == "S2":
            assert obj["allele_string"] == "G/A/T"
            obj["allele_string"] = "G/T/A"
            assert [tc["variant_allele"] for tc in obj["transcript_consequences"]] == ["A", "T"] * 4
            obj["transcript_consequences"] = swapped_pairs(obj["transcript_consequences"])
        objects.append(json.dumps(obj))
    json_expected = cvf.parse_json(json_text, inputs)
    json_actual = cvf.parse_json("\n".join(objects) + "\n", inputs)
    plain = cvf.compare_entries(json_expected, json_actual, no_documented())
    assert plain["verdict"] == "red"
    assert plain["columns"]["allele_string"]["differing"] == 1
    documented = with_class()
    json_report = cvf.compare_entries(json_expected, json_actual, documented)
    assert json_report["verdict"] == "green"
    assert json_report["columns"]["allele_string"]["individual_allele_order"] == 1
    assert json_report["columns"]["allele_string"]["individual_allele_order_examples"] == ["obj434|record: expected 'G/A/T', got 'G/T/A'"]
    assert json_report["counts"]["individual_entries_reordered"] == TWO_ALLELE_ROWS
    assert documented.applied == {cvf.INDIVIDUAL_ALLELE_ORDER: TWO_ALLELE_ROWS + 1}
    changed = [o.replace('"allele_string": "G/T/A"', '"allele_string": "G/T/C"') for o in objects]
    changed_report = cvf.compare_entries(json_expected, cvf.parse_json("\n".join(changed) + "\n", inputs), with_class())
    assert changed_report["verdict"] == "red"
    assert changed_report["reasons"] == ["field allele_string: 1 differing entries"]

    # The VCF format's CSQ entries of the copy's line, each transcript's pair in the other order.
    vcf = cvf.read_text(expected_dir / "vcf.vcf.gz")
    flipped = []
    for row in data_lines(vcf):
        cols = row.split("\t")
        if cols[2] == "hb_multi" and "|S2|" in cols[7]:
            prefix, marker, csq = cols[7].partition("CSQ=")
            chunks = csq.split(",")
            assert [chunk.split("|")[0] for chunk in chunks] == ["A", "T"] * 4
            cols[7] = prefix + marker + ",".join(swapped_pairs(chunks))
        flipped.append("\t".join(cols))
    vcf_expected = cvf.parse_vcf(vcf, inputs)
    vcf_actual = cvf.parse_vcf(rebuild(vcf, flipped), inputs)
    assert cvf.compare_entries(vcf_expected, vcf_actual, no_documented())["verdict"] == "red"
    documented = with_class()
    vcf_report = cvf.compare_entries(vcf_expected, vcf_actual, documented)
    assert vcf_report["verdict"] == "green"
    assert vcf_report["counts"]["individual_entries_reordered"] == TWO_ALLELE_ROWS
    assert documented.applied == {cvf.INDIVIDUAL_ALLELE_ORDER: TWO_ALLELE_ROWS}

    # --process_ref_homs --allele_number --uploaded_allele: the copy's rows are G, A and T per
    # transcript, numbered 0, 1 and 2 in the allele string's order (OutputFactory.pm:1354),
    # UPLOADED_ALLELE the string itself (OutputFactory.pm:1360); the other order numbers T 1
    # and A 2 and prints G/T/A on every row of the copy, the G row included.
    fields_default = cvf.read_text(INDIVIDUAL / "expected_individual-refhoms-fields" / "default.txt.gz")
    rows = data_lines(fields_default)
    targets = s2_rows(rows)
    assert [rows[i].split("\t")[2] for i in targets] == ["G", "A", "T"] * 4
    renumbered = list(rows)
    for i in targets:
        row = rows[i].replace("UPLOADED_ALLELE=G/A/T;", "UPLOADED_ALLELE=G/T/A;")
        if "\tA\t" in row:
            row = row.replace("ALLELE_NUM=1;", "ALLELE_NUM=2;")
        else:
            row = row.replace("ALLELE_NUM=2;", "ALLELE_NUM=1;")
        renumbered[i] = row
    for g in targets[0::3]:
        renumbered[g + 1], renumbered[g + 2] = renumbered[g + 2], renumbered[g + 1]
    fields_expected = cvf.parse_default(fields_default)
    fields_actual = cvf.parse_default(rebuild(fields_default, renumbered))
    plain = cvf.compare_entries(fields_expected, fields_actual, no_documented())
    assert plain["verdict"] == "red"
    assert plain["columns"]["ALLELE_NUM"]["differing"] == TWO_ALLELE_ROWS
    assert plain["columns"]["UPLOADED_ALLELE"]["differing"] == len(targets)
    documented = with_class()
    fields_report = cvf.compare_entries(fields_expected, fields_actual, documented)
    assert fields_report["verdict"] == "green"
    assert fields_report["columns"]["ALLELE_NUM"]["individual_allele_order"] == TWO_ALLELE_ROWS
    assert fields_report["columns"]["UPLOADED_ALLELE"]["individual_allele_order"] == len(targets)
    assert fields_report["columns"]["Extra"]["individual_allele_order"] == len(targets)
    assert fields_report["counts"]["individual_entries_reordered"] == TWO_ALLELE_ROWS
    # The eight reordered rows, their eight ALLELE_NUM values, and UPLOADED_ALLELE with the
    # Extra text carrying it on all twelve rows of the copy.
    assert documented.applied == {cvf.INDIVIDUAL_ALLELE_ORDER: 2 * TWO_ALLELE_ROWS + 2 * len(targets)}
    assert fields_report["columns"]["ALLELE_NUM"]["individual_allele_order_examples"][:2] == [
        "21:34716340|A|ENST00000349499|Transcript#2: expected '1', got '2'",
        "21:34716340|T|ENST00000349499|Transcript: expected '2', got '1'",
    ]
    # Numbering both alternate alleles 1 is not a permutation of the copy's numbers.
    doubled = [r.replace("ALLELE_NUM=2;", "ALLELE_NUM=1;") if i in targets else r for i, r in enumerate(renumbered)]
    doubled_report = cvf.compare_entries(fields_expected, cvf.parse_default(rebuild(fields_default, doubled)), with_class())
    assert doubled_report["verdict"] == "red"
    assert doubled_report["reasons"] == ["field ALLELE_NUM: 4 differing entries", "field Extra: 4 differing entries"]
    assert doubled_report["columns"]["ALLELE_NUM"]["individual_allele_order"] == 0
    # A renumbering of one transcript alone is no permutation applied to the copy: each allele
    # carries one number across the copy's rows, its index in the copy's one allele string
    # (OutputFactory.pm:1354).
    partial = list(rows)
    partial[targets[1]] = rows[targets[1]].replace("ALLELE_NUM=1;", "ALLELE_NUM=2;")
    partial[targets[2]] = rows[targets[2]].replace("ALLELE_NUM=2;", "ALLELE_NUM=1;")
    partial_report = cvf.compare_entries(fields_expected, cvf.parse_default(rebuild(fields_default, partial)), with_class())
    assert partial_report["verdict"] == "red"
    assert partial_report["columns"]["ALLELE_NUM"]["differing"] == 2
    assert partial_report["columns"]["ALLELE_NUM"]["individual_allele_order"] == 0
    # The reference allele keeps 0 and its row keeps its place (Parser/VCF.pm:639 writes the
    # reference first): G numbered 1 with A numbered 0, or the G row after the alternates, is red.
    shifted_zero = list(rows)
    for g in targets[0::3]:
        shifted_zero[g] = rows[g].replace("ALLELE_NUM=0;", "ALLELE_NUM=1;")
        shifted_zero[g + 1] = rows[g + 1].replace("ALLELE_NUM=1;", "ALLELE_NUM=0;")
    zero_report = cvf.compare_entries(fields_expected, cvf.parse_default(rebuild(fields_default, shifted_zero)), with_class())
    assert zero_report["verdict"] == "red"
    assert zero_report["columns"]["ALLELE_NUM"]["differing"] == TWO_ALLELE_ROWS
    assert zero_report["columns"]["ALLELE_NUM"]["individual_allele_order"] == 0
    reference_last = list(rows)
    for g in targets[0::3]:
        reference_last[g], reference_last[g + 1], reference_last[g + 2] = rows[g + 1], rows[g + 2], rows[g]
    last_report = cvf.compare_entries(fields_expected, cvf.parse_default(rebuild(fields_default, reference_last)), with_class())
    assert last_report["verdict"] == "red"
    assert last_report["reasons"] == [f"order differs at shared position {targets[0]}"]

    # A JSON consequence entry names the reference allele as ref_allele: without allele_num the
    # reference row is still known by it, so moved after the alternates it reads red while the
    # alternates' own order stays free.
    def without_numbers(obj: dict) -> dict:
        for tc in obj.get("transcript_consequences", []):
            tc.pop("allele_num", None)
        return obj

    def reordered(obj: dict, order: tuple[int, int, int]) -> str:
        tcs = obj["transcript_consequences"]
        return json.dumps({**obj, "transcript_consequences": [tcs[i + j] for i in range(0, len(tcs), 3) for j in order]})

    fields_json = cvf.read_text(INDIVIDUAL / "expected_individual-refhoms-fields" / "json.jsonl.gz")
    plain_objects, reference_last_objects, alternates_objects = [], [], []
    for line in cvf.text_lines(fields_json):
        if not line:
            continue
        obj = without_numbers(json.loads(line))
        text = json.dumps(obj)
        plain_objects.append(text)
        if obj.get("id") == "hb_multi" and obj.get("ind") == "S2":
            tcs = obj["transcript_consequences"]
            assert [tc["variant_allele"] for tc in tcs] == ["G", "A", "T"] * 4
            assert all(tc["ref_allele"] == "G" for tc in tcs)
            reference_last_objects.append(reordered(obj, (1, 2, 0)))
            alternates_objects.append(reordered(obj, (0, 2, 1)))
        else:
            reference_last_objects.append(text)
            alternates_objects.append(text)
    json_fields_expected = cvf.parse_json("\n".join(plain_objects) + "\n", inputs)
    reference_last_json = cvf.compare_entries(
        json_fields_expected, cvf.parse_json("\n".join(reference_last_objects) + "\n", inputs), with_class()
    )
    assert reference_last_json["verdict"] == "red"
    assert reference_last_json["reasons"][0].startswith("order differs at shared position")
    alternates_json = cvf.compare_entries(
        json_fields_expected, cvf.parse_json("\n".join(alternates_objects) + "\n", inputs), with_class()
    )
    assert alternates_json["verdict"] == "green"
    assert alternates_json["counts"]["individual_entries_reordered"] == TWO_ALLELE_ROWS
    # A number missing on one side is no permutation: the row's own value must be its allele's
    # number in the copy.
    one_missing = []
    for line in cvf.text_lines(fields_json):
        if not line:
            continue
        obj = json.loads(line)
        if obj.get("id") == "hb_multi" and obj.get("ind") == "S2":
            assert obj["transcript_consequences"][0].pop("allele_num") == 0
        one_missing.append(json.dumps(obj))
    one_missing_report = cvf.compare_entries(
        cvf.parse_json(fields_json, inputs), cvf.parse_json("\n".join(one_missing) + "\n", inputs), with_class()
    )
    assert one_missing_report["verdict"] == "red"
    assert one_missing_report["reasons"] == ["field allele_num: 1 differing entries"]
    assert one_missing_report["columns"]["allele_num"]["individual_allele_order"] == 0
    # An insertion's reference row prints REF_ALLELE=-, its own allele: without a number it is
    # still the reference row, and moved after the alternates it reads red.
    indel_targets = [i for i, r in enumerate(rows) if r.startswith("hb_multi_indel\t") and "IND=S4;" in r]
    assert [rows[i].split("\t")[2] for i in indel_targets][:3] == ["-", "A", "AA"]
    unnumbered = [re.sub(r"ALLELE_NUM=\d+;", "", r) if i in indel_targets else r for i, r in enumerate(rows)]
    indel_last = list(unnumbered)
    for g in indel_targets[0::3]:
        indel_last[g], indel_last[g + 1], indel_last[g + 2] = unnumbered[g + 1], unnumbered[g + 2], unnumbered[g]
    indel_expected = cvf.parse_default(rebuild(fields_default, unnumbered))
    indel_report = cvf.compare_entries(indel_expected, cvf.parse_default(rebuild(fields_default, indel_last)), with_class())
    assert indel_report["verdict"] == "red"
    assert indel_report["reasons"] == [f"order differs at shared position {indel_targets[0]}"]


def test_json_floors_count_under_the_json_keys(tmp_path):
    """A floor names the default-format column in every format: in JSON it is counted under
    the column's JSON key (DOMAINS under domains, ENSP under protein_id, SIFT and PolyPhen
    under either of their two keys, so ``--sift s`` counts the bare score), and a
    record-level column, ``Feature_type`` (the consequence list's name, no key of an entry)
    or ``Extra`` (the text formats' container) has no row count there and is reported as not
    applicable instead of red; the text formats count the column itself."""
    assert cvf.json_floor_fields("DOMAINS") == ("domains",)
    assert cvf.json_floor_fields("ENSP") == ("protein_id",)
    assert cvf.json_floor_fields("SYMBOL") == ("gene_symbol",)
    assert cvf.json_floor_fields("SIFT") == ("sift_prediction", "sift_score")
    assert cvf.json_floor_fields("PolyPhen") == ("polyphen_prediction", "polyphen_score")
    assert cvf.json_floor_fields("Protein_position") == ("protein_start", "protein_end")
    assert cvf.json_floor_fields("Feature") == ("transcript_id", "regulatory_feature_id", "motif_feature_id")
    for column in ("Existing_variation", "AF", "gnomADe_NFE_AF", "MAX_AF", "FREQS", "CLIN_SIG", "Feature_type", "Extra"):
        assert cvf.json_floor_fields(column) is None

    inputs = cvf.input_records(PROTEIN_FEATURES / "variants.vcf")
    floors = {"DOMAINS": 157, "SIFT": 219, "PolyPhen": 163, "ENSP": 451, "Existing_variation": 1, "Feature_type": 1}
    json_path = PROTEIN_FEATURES / "expected" / "json.jsonl.gz"
    report = cvf.compare_files("json", json_path, json_path, inputs, no_documented(), set(), floors, {}, 4)
    assert report["verdict"] == "green"
    assert report["reasons"] == []
    assert report["floors_not_applicable"] == ["Existing_variation", "Feature_type"]
    assert report["columns"]["Existing_variation"] == {
        "expected_populated": 0, "actual_populated": 0, "differing": 0, "pairs_accepted": 0,
        "nearest_accepted": 0, "order_only": 0, "individual_allele_order": 0, "floor": 1, "floor_met": None,
        "floor_note": cvf.FLOOR_NOT_APPLICABLE,
    }
    assert report["columns"]["Feature_type"]["floor_met"] is None
    assert report["columns"]["Feature_type"]["floor_note"] == cvf.FLOOR_NOT_APPLICABLE
    for column, fields in (("DOMAINS", ["domains"]), ("SIFT", ["sift_prediction", "sift_score"]),
                           ("PolyPhen", ["polyphen_prediction", "polyphen_score"]), ("ENSP", ["protein_id"])):
        col = report["columns"][column]
        assert col["floor_met"] is True
        assert col["floor_fields"] == fields
        assert col["expected_populated"] == col["actual_populated"] == floors[column]
    # The same counts the default format reaches under the column names themselves.
    default_path = PROTEIN_FEATURES / "expected" / "default.txt.gz"
    text_report = cvf.compare_files("default", default_path, default_path, None, no_documented(), set(), floors, {}, 4)
    assert text_report["floors_not_applicable"] == []
    assert text_report["reasons"] == ["floor Existing_variation=1: expected populated 0, actual populated 0"]
    assert text_report["columns"]["Feature_type"]["floor_met"] is True
    assert "floor_fields" not in text_report["columns"]["DOMAINS"]
    assert text_report["columns"]["DOMAINS"]["expected_populated"] == 157

    # A floor the JSON counts fall short of is red, as in the text formats.
    short = cvf.compare_files("json", json_path, json_path, inputs, no_documented(), set(), {"DOMAINS": 158}, {}, 4)
    assert short["verdict"] == "red"
    assert short["reasons"] == ["floor DOMAINS=158: expected populated 157, actual populated 157"]

    # --sift s prints the bare score, which the JSON format carries as sift_score alone.
    score_corpus = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-sift-score"
    score_inputs = cvf.input_records(score_corpus / "variants.vcf")
    score_json = score_corpus / "expected" / "json.jsonl.gz"
    score_report = cvf.compare_files("json", score_json, score_json, score_inputs, no_documented(), set(), {"SIFT": 84, "PolyPhen": 58}, {}, 4)
    assert score_report["verdict"] == "green"
    assert score_report["columns"]["SIFT"]["expected_populated"] == 84
    assert score_report["columns"]["sift_score"]["expected_populated"] == 84
    assert "sift_prediction" not in score_report["columns"]

    # Through the command line, the floors of a fixture row on its JSON format.
    out = tmp_path / "json.json"
    code = cvf.main([
        "--format", "json", "--expected", str(json_path), "--actual", str(json_path), "--input",
        str(PROTEIN_FEATURES / "variants.vcf"), "--floors", "DOMAINS=157,SIFT=219,PolyPhen=163,ENSP=451,Existing_variation=1",
        "--masks", cvf.ORDER_ONLY, "--report", str(out),
    ])
    assert code == 0
    written = json.loads(out.read_text(encoding="utf-8"))
    assert written["floors_not_applicable"] == ["Existing_variation"]
    assert written["masks"]["unused"] == [cvf.ORDER_ONLY]


def json_floor_counts_match_the_default(corpus: Path) -> dict[str, tuple[int, int]]:
    """Per default-format column of a corpus reference with a JSON floor mapping, the
    default format's populated count beside the JSON format's floor count."""
    inputs = cvf.input_records(corpus / "variants.vcf")
    default = cvf.parse_default(cvf.read_text(corpus / "expected" / "default.txt.gz"))
    objects = cvf.parse_json(cvf.read_text(corpus / "expected" / "json.jsonl.gz"), inputs)
    default_counts = cvf.populated_counts(default)
    out = {}
    for column in sorted({name for e in default for name in e.fields}):
        fields = cvf.json_floor_fields(column)
        if fields is not None:
            out[column] = (default_counts.get(column, 0), cvf.populated_under(objects, fields))
    return out


def test_json_floor_counts_equal_the_default_counts_on_every_mapped_column():
    """On the GRCh38 and GRCh38-everything references of Ensembl VEP 116.2, every default
    column with a JSON floor mapping counts the same under its JSON fields as in the default
    format: the position keys survive a ``?`` start (JSON.pm:340-341 keep the end), the
    record objects' own ``strand`` is not counted, a ``domains`` array left empty by the
    re-split still marks a populated row, and a present key is a populated column."""
    corpus = HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38"
    grch38 = json_floor_counts_match_the_default(corpus)
    assert grch38["cDNA_position"] == (8879, 8879)
    assert grch38["CDS_position"] == (5907, 5907)
    assert grch38["STRAND"] == (20621, 20621)
    mismatched = {column: counts for column, counts in grch38.items() if counts[0] != counts[1]}
    assert mismatched == {}
    # A floor spelled with the JSON key itself is counted the same way, over the
    # consequence entries alone: the record objects carry their own strand.
    inputs = cvf.input_records(corpus / "variants.vcf")
    json_path = corpus / "expected" / "json.jsonl.gz"
    report = cvf.compare_files("json", json_path, json_path, inputs, no_documented(), set(), {"strand": 20621}, {}, 4)
    assert report["columns"]["strand"]["floor_fields"] == ["strand"]
    assert report["columns"]["strand"]["expected_populated"] == 20621 and report["verdict"] == "green"
    everything = json_floor_counts_match_the_default(HERE.parent.parent / "tests" / "golden" / "116" / "GRCh38-everything")
    assert len(everything) == 40
    assert everything["DOMAINS"] == (2055, 2055)
    assert everything["Protein_position"] == (2734, 2734)
    mismatched = {column: counts for column, counts in everything.items() if counts[0] != counts[1]}
    assert mismatched == {}


def test_self_parity_sets_freqs_aside_as_a_record_level_column():
    """FREQS is a record-level column the JSON format never carries (JSON.pm:493-495), so
    the self-parity comparison of one engine's formats lists it under the columns not
    compared against JSON instead of reporting the missing key as a difference."""
    expected = FREQFILTER / "expected_freqfilter"
    report = cvf.SelfParity(
        cvf.read_text(expected / "default.txt.gz"),
        cvf.read_text(expected / "tab.txt.gz"),
        cvf.read_text(expected / "vcf.vcf.gz"),
        cvf.read_text(expected / "json.jsonl.gz"),
        inputs=cvf.input_records(FREQFILTER / "variants.vcf"),
    ).run()
    assert report["verdict"] == "green", report["reasons"]
    assert "FREQS" in report["not_compared"]["default_record_level"]
    assert "FREQS" in report["pairs"]["vcf"]["columns_compared"]
