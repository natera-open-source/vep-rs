"""Tests of the corpus fixture scripts beside this file on synthetic inputs: `prune_json_cache.py`,
`make_custom_fixtures.py` and `extract_reference.py`."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import extract_reference as xr  # noqa: E402
import make_custom_fixtures as mcf  # noqa: E402
import prune_json_cache as pjc  # noqa: E402


CUSTOM_FIXTURE_VCF = (
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "21\t100\trs1\tA\tG\t.\tPASS\t.\n"  # exact
    "21\t200\trs2\tTCA\tT\t.\tPASS\t.\n"  # overlap
    "21\t300\trs3\tC\tA,T\t.\tPASS\t.\n"  # miss
    "21\t300\trs3b\tC\tG\t.\tPASS\t.\n"  # the locus is written once
    "21\t400\trs4\tA\tAT\t.\tPASS\t.\n"  # multi_allelic
    "21\t500\trs5\tG\tC\t.\tPASS\t.\n"  # exact again
    "21\t599\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=900\n"  # sv_span, inside 600..899
    "21\t1000\tbnd1\tA\tA[22:5000[\t.\tPASS\tSVTYPE=BND\n"  # no span
    "22\t50\trs6\tT\tC\t.\tPASS\t.\n"  # overlap, the next chromosome
)


def test_make_custom_fixtures_records_follow_the_class_cycle():
    """A corpus's distinct sequence loci take the classes exact, overlap, miss, multi_allelic in
    turn (a second variant at a locus gets no record), an overlap record keeps POS and REF and
    changes the first ALT's last base (`TCA>T` to `TCA>A`, so `type=exact` sees another allele
    and `type=overlap` the same span: File/VCF.pm:385-420 at release/116.3), a multi-allelic
    record lists the first ALT then that allele with one AF each, a structural variant with END
    gets one single-base record inside its span at a position no sequence variant holds, a
    breakend none; records sort by chromosome then position for tabix."""
    records = [line.split("\t") for line in CUSTOM_FIXTURE_VCF.splitlines() if not line.startswith("#")]
    sites, counts = mcf.derive_sites(records)
    assert counts == {"exact": 2, "overlap": 2, "miss": 1, "multi_allelic": 1, "sv_span": 1}
    assert [s.vcf_line() for s in sites] == [
        "21\t100\texact_1\tA\tG\t.\tPASS\tAF=0.001;CLASS=exact",
        "21\t200\toverlap_1\tTCA\tA\t.\tPASS\tAF=0.002;CLASS=overlap",
        "21\t400\tmulti_1\tA\tAT,AA\t.\tPASS\tAF=0.003,0.503;CLASS=multi_allelic",
        "21\t500\texact_2\tG\tC\t.\tPASS\tAF=0.004;CLASS=exact",
        "21\t600\tsv_1\tN\tA\t.\tPASS\tAF=0.006;CLASS=sv_span",
        "22\t50\toverlap_2\tT\tA\t.\tPASS\tAF=0.005;CLASS=overlap",
    ]
    # The BED score is the record's ordinal, the AF in thousandths.
    assert [s.bed_line() for s in sites][:2] == ["21\t99\t100\texact_1\t1\t+", "21\t199\t202\toverlap_1\t2\t+"]
    assert [s.bedgraph_line() for s in sites][1] == "21\t199\t200\t0.002"
    # The UCSC copies spell the chromosome as UCSC does and score the second record with its
    # ordinal in hundredths minus one (phyloP) and one minus its AF (phastCons).
    assert sites[1].bed_line(mcf.ucsc_chrom(sites[1].chrom)) == "chr21\t199\t202\toverlap_1\t2\t+"
    assert sites[1].bedgraph_line("chr21", sites[1].phylop()) == "chr21\t199\t200\t-0.98"
    assert sites[1].bedgraph_line("chr21", sites[1].phastcons()) == "chr21\t199\t200\t0.998"
    assert mcf.ucsc_chrom("21") == "chr21" and mcf.ucsc_chrom("MT") == "chrM" and mcf.ucsc_chrom("chr1") == "chr1"
    # A promoter window five bases either side of every record, every third without a gene id,
    # the fifth with a second window twenty bases either side, and the enhancer over the first
    # record, sorted by chromosome and start; a window starts at base 1 at the earliest.
    assert mcf.promoter_lines(sites) == [
        "21\t.\tpromoter_window\t95\t105\t.\t+\t.\tID=PW1;gene_id=ENSG90000000001",
        "21\t.\tenhancer\t95\t105\t.\t+\t.\tID=EH1",
        "21\t.\tpromoter_window\t195\t207\t.\t+\t.\tID=PW2;gene_id=ENSG90000000002",
        "21\t.\tpromoter_window\t395\t405\t.\t+\t.\tID=PW3",
        "21\t.\tpromoter_window\t495\t505\t.\t+\t.\tID=PW4;gene_id=ENSG90000000004",
        "21\t.\tpromoter_window\t580\t620\t.\t+\t.\tID=PW5w;gene_id=ENSG90000000005",
        "21\t.\tpromoter_window\t595\t605\t.\t+\t.\tID=PW5;gene_id=ENSG90000000005",
        "22\t.\tpromoter_window\t45\t55\t.\t+\t.\tID=PW6",
    ]
    early = [mcf.Site("1", 3, "exact_1", "A", ("G",), ("0.001",), "exact")]
    assert mcf.promoter_lines(early)[0] == "1\t.\tpromoter_window\t1\t8\t.\t+\t.\tID=PW1;gene_id=ENSG90000000001"
    # A single-base SNV's other allele avoids REF as well as the original ALT.
    assert mcf.other_allele("G", "A") == "C" and mcf.other_allele("A", "G") == "C" and mcf.other_allele("C", "A") == "G"
    assert mcf.other_allele("T", "TCA") == "A" and mcf.other_allele("AT", "A") == "AA"
    assert mcf.structural_end(10, "N", "<DEL>", "SVTYPE=DEL;END=20") == 20
    assert mcf.structural_end(10, "N", "<DEL>", "SVTYPE=DEL;SVLEN=-5") == 15
    assert mcf.structural_end(10, "ACGTACGT", "<CN0>", ".") == 17
    assert mcf.structural_end(10, "N", "N[22:5[", "SVTYPE=BND") is None
    assert mcf.structural_end(10, "N", ".N", "") is None
    assert mcf.contig_length("GRCh38", "chr21", 5) == 46709983 and mcf.contig_length("GRCh37", "21", 5) == 48129895
    assert mcf.contig_length("GRCh38", "NC_000021.9", 5000) == 5000


def test_make_custom_fixtures_writes_indexed_files(tmp_path: Path, capsys):
    """The tool writes the VCF, BED and GFF sources as BGZF (the gzip magic with the extra-field
    flag, `1f 8b 08 04`) with their tabix indexes beside them, the bigWig with its magic,
    records the tool versions and every file's sha256 under provenance.json's custom_fixtures
    key while keeping the record's other keys, and refuses to overwrite without --force."""
    missing = [tool for tool in mcf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    (corpus / "variants.vcf").write_text(CUSTOM_FIXTURE_VCF, encoding="utf-8")
    (corpus / "manifest.json").write_text(json.dumps({"assembly": "GRCh38"}), encoding="utf-8")
    (corpus / "provenance.json").write_text(json.dumps({"vep_image": "x", "input": "hand-built records"}), encoding="utf-8")
    assert mcf.main(["--corpus", str(corpus)]) == 0
    assert capsys.readouterr().out.strip().endswith("loci per class: exact 2, overlap 2, miss 1, multi_allelic 1, sv_span 1")
    out = corpus / "custom"
    for name in ("sites.vcf.gz", "regions.bed.gz", "promoters.gff.gz", "clinvar_source.vcf.gz", "regions_chr.bed.gz"):
        assert (out / name).read_bytes()[:4] == b"\x1f\x8b\x08\x04", name
        assert (out / f"{name}.tbi").is_file(), name
        assert (out / f"{name}.tbi").read_bytes()[:2] == b"\x1f\x8b", name
    with gzip.open(out / "clinvar_source.vcf.gz", "rt", encoding="utf-8") as fh:
        clinvar_lines = fh.read().splitlines()
    assert clinvar_lines[:2] == ["##fileformat=VCFv4.2", "##source=ClinVar"]
    with gzip.open(out / "regions_chr.bed.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines()[:2] == ["chr21\t99\t100\texact_1\t1\t+", "chr21\t199\t202\toverlap_1\t2\t+"]
    for name in (mcf.PHYLOP_BW, mcf.PHASTCONS_BW):
        assert (out / name).read_bytes()[:4] == b"\x26\xfc\x8f\x88", name
    with gzip.open(out / "promoters.gff.gz", "rt", encoding="utf-8") as fh:
        gff_lines = fh.read().splitlines()
    assert gff_lines[:2] == ["##gff-version 3", "21\t.\tpromoter_window\t95\t105\t.\t+\t.\tID=PW1;gene_id=ENSG90000000001"]
    assert len(gff_lines) == 9
    assert (out / "scores.bw").read_bytes()[:4] == b"\x26\xfc\x8f\x88"
    with gzip.open(out / "sites.vcf.gz", "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2" and "##contig=<ID=21>" in lines and "##contig=<ID=22>" in lines
    assert [line.split("\t")[2] for line in lines if not line.startswith("#")] == [
        "exact_1", "overlap_1", "multi_1", "exact_2", "sv_1", "overlap_2",
    ]
    with gzip.open(out / "regions.bed.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines()[0] == "21\t99\t100\texact_1\t1\t+"
    assert sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()) == sorted(mcf.FILES), "no scratch file survives"
    record = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    assert record["vep_image"] == "x" and record["input"] == "hand-built records"
    fixtures = record["custom_fixtures"]
    assert fixtures["generator"] == "scripts/golden/make_custom_fixtures.py"
    assert fixtures["variants_vcf_sha256"] == hashlib.sha256(CUSTOM_FIXTURE_VCF.encode()).hexdigest()
    assert fixtures["assembly"] == "GRCh38"
    assert fixtures["loci"] == {"exact": 2, "overlap": 2, "miss": 1, "multi_allelic": 1, "sv_span": 1}
    assert sorted(fixtures["files"]) == sorted(f"custom/{name}" for name in mcf.FILES)
    for name, digest in fixtures["files"].items():
        assert hashlib.sha256((corpus / name).read_bytes()).hexdigest() == digest, name
    assert sorted(fixtures["tools"]) == ["bgzip", "bigtools", "tabix"]
    assert all(version.split()[-1][0].isdigit() for version in fixtures["tools"].values()), fixtures["tools"]
    with pytest.raises(SystemExit, match="pass --force"):
        mcf.main(["--corpus", str(corpus)])
    assert mcf.main(["--corpus", str(corpus), "--force"]) == 0
    assert json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["custom_fixtures"]["files"] == fixtures["files"], \
        "a regeneration from the same input writes the same bytes"


def test_extract_reference_writes_the_contig_and_its_index(tmp_path: Path):
    src = tmp_path / "genome.fa"
    src.write_text(">1 dna:chromosome\nACGTACGTAC\nGGGG\n>21 dna:chromosome chromosome:GRCh37:21:1:14:1 REF\nAACCGGTTAA\nCCGN\n>22\nTTTT\n", encoding="utf-8")
    out = tmp_path / "reference.fa"
    rc = xr.main([str(src), "--contig", "21", "--out", str(out)])
    assert rc == 0
    assert out.read_text(encoding="utf-8") == ">21\nAACCGGTTAA\nCCGN\n"
    assert Path(str(out) + ".fai").read_text(encoding="utf-8") == "21\t14\t4\t10\t11\n"
    with pytest.raises(SystemExit):
        xr.main([str(src), "--contig", "MT", "--out", str(tmp_path / "none.fa")])


def test_prune_keeps_reachable_transcripts_in_every_shard(tmp_path: Path):
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    t_a = {"stable_id": "ENST_A", "start": "10000", "end": "20000", "dbID": "1"}
    t_b = {"stable_id": "ENST_B", "start": "990000", "end": "1010000", "dbID": "2"}  # straddles the 1M boundary
    t_c = {"stable_id": "ENST_C", "start": "500000", "end": "500500", "dbID": "3"}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([t_a, t_b, t_c]))
    (cache / "transcripts" / "21" / "1000001-2000000.json").write_text(json.dumps([t_b]))
    (cache / "transcripts" / "22").mkdir()
    (cache / "transcripts" / "22" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_D", "start": "5", "end": "9", "dbID": "4"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t24000\t.\tA\tG\t.\t.\t.\n21\t995000\t.\tT\t<DEL>\t.\t.\tEND=996000\n")
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh37", "--cache-version", "115"])
    assert rc == 0
    s1 = json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())
    s2 = json.loads((out / "transcripts" / "21" / "1000001-2000000.json").read_text())
    assert [t["stable_id"] for t in s1] == ["ENST_A", "ENST_B"]  # A within 5 kb flank, C unreachable
    assert [t["stable_id"] for t in s2] == ["ENST_B"]  # the straddling copy survives in both shards
    assert not (out / "transcripts" / "22").exists()
    info = json.loads((out / "info.json").read_text())
    assert info["assembly"] == "GRCh37" and info["cache_version"] == 115
    for name in ("1-1000000.json.gz", "1-1000000.json"):
        assert name.removesuffix(".gz") == (name[: -len(".gz")] if name.endswith(".gz") else name) == "1-1000000.json"  # noqa: FURB188


def test_prune_keeps_the_same_symbol_hgnc_donors_of_an_unnamed_model(tmp_path: Path, capsys):
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    reached = {"stable_id": "NR_A", "start": "10000", "end": "20000", "dbID": "1", "source_cache": "RefSeq",
               "gene_symbol": "ANKRD30BP2", "gene_symbol_source": "EntrezGene"}
    donor = {"stable_id": "ENST_DONOR", "start": "500000", "end": "500500", "dbID": "2", "source_cache": "Ensembl",
             "gene_symbol": "ANKRD30BP2", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:16620"}
    same_symbol_unnamed = {"stable_id": "ENST_UNNAMED", "start": "600000", "end": "600500", "dbID": "3",
                           "source_cache": "Ensembl", "gene_symbol": "ANKRD30BP2"}
    other_symbol = {"stable_id": "ENST_OTHER", "start": "700000", "end": "700500", "dbID": "4", "source_cache": "Ensembl",
                    "gene_symbol": "MRPS6", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:14051"}
    named_reached = {"stable_id": "ENST_B", "start": "1200000", "end": "1210000", "dbID": "5", "source_cache": "Ensembl",
                     "gene_symbol": "MRPS6", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:14051"}
    unnamed_reached = {"stable_id": "NM_B", "start": "1200000", "end": "1210000", "dbID": "6", "source_cache": "RefSeq",
                       "gene_symbol": "MRPS6", "gene_symbol_source": "EntrezGene"}
    spare_donor = {"stable_id": "ENST_SPARE", "start": "1500000", "end": "1500500", "dbID": "7", "source_cache": "Ensembl",
                   "gene_symbol": "MRPS6", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:14051"}
    far_donor = {"stable_id": "ENST_FAR", "start": "1600000", "end": "1600500", "dbID": "8", "source_cache": "Ensembl",
                 "gene_symbol": "ANKRD30BP2", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:16620"}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([reached, donor, same_symbol_unnamed, other_symbol]))
    (cache / "transcripts" / "21" / "1000001-2000000.json").write_text(json.dumps([named_reached, unnamed_reached, spare_donor, far_donor]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t24000\t.\tA\tG\t.\t.\t.\n21\t1205000\t.\tT\tC\t.\t.\t.\n")
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38", "--cache-version", "116"])
    assert rc == 0
    s1 = json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())
    s2 = json.loads((out / "transcripts" / "21" / "1000001-2000000.json").read_text())
    # The donor stays in the shard's order; a same-symbol model without an id and a model of
    # another symbol are not donors; a symbol a reached model already carries with its id pulls
    # nothing in, and a donor is read from the unnamed model's own shard alone.
    assert [t["stable_id"] for t in s1] == ["NR_A", "ENST_DONOR"]
    assert [t["stable_id"] for t in s2] == ["ENST_B", "NM_B"]
    assert json.loads(capsys.readouterr().out)["hgnc_donors"] == 1


def test_prune_regulatory_keeps_the_reachable_features_and_declares_them(tmp_path: Path, capsys):
    """`--regulatory` writes the regulatory and motif features within the flank of a record
    span, in shards named as in the source, the binding matrices the kept motifs name, and
    `info.json` declares `regulatory` with the schema and cell-type list only when a shard was
    written; without the flag the tree is left out and the flag reads false."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116,
                                                 "regulatory": True, "regulatory_schema": 1, "cell_types": ["A549", "HUVEC"]}))
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "10000", "end": "20000", "dbID": "1"}]))
    (cache / "regulatory" / "21").mkdir(parents=True)
    near = {"stable_id": "ENSR_NEAR", "start": 26000, "end": 27000, "feature_type": "RegulatoryFeature", "biotype": "enhancer", "activity": "AI"}
    far = {"stable_id": "ENSR_FAR", "start": 500000, "end": 500500, "feature_type": "RegulatoryFeature", "biotype": "promoter", "activity": "IA"}
    motif = {"stable_id": "ENSM_IN", "start": 23990, "end": 24010, "strand": 1, "feature_type": "MotifFeature", "matrix": "ENSPFM1"}
    far_motif = {"stable_id": "ENSM_FAR", "start": 600000, "end": 600010, "strand": 1, "feature_type": "MotifFeature", "matrix": "ENSPFM2"}
    (cache / "regulatory" / "21" / "1-1000000.json").write_text(json.dumps([near, far, motif, far_motif]))
    (cache / "regulatory" / "binding_matrices.json").write_text(json.dumps(
        {"ENSPFM1": {"transcription_factors": ["TF"], "frequencies": []}, "ENSPFM2": {"transcription_factors": ["TG"], "frequencies": []}}))
    with gzip.open(cache / "regulatory" / "21" / "1000001-2000000.json.gz", "wt") as fh:
        fh.write(json.dumps([{"stable_id": "ENSR_OTHER", "start": 1500000, "end": 1500100, "strand": 0, "feature_type": "RegulatoryFeature"}]))
    (cache / "regulatory" / "22").mkdir()
    (cache / "regulatory" / "22" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENSR_22", "start": 5, "end": 9, "strand": 0, "feature_type": "RegulatoryFeature"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t24000\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--regulatory", "--gzip"])
    assert rc == 0
    summary = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert summary["regulatory_features_per_chromosome"] == {"21": 2} and summary["regulatory_shards"] == 1
    assert summary["binding_matrices"] == 1
    with gzip.open(out / "regulatory" / "21" / "1-1000000.json.gz", "rt") as fh:
        kept = json.load(fh)
    assert [f["stable_id"] for f in kept] == ["ENSR_NEAR", "ENSM_IN"]  # the far feature is beyond the 5 kb flank
    assert kept[1]["matrix"] == "ENSPFM1"
    assert json.loads((out / "regulatory" / "binding_matrices.json").read_text()) == {"ENSPFM1": {"transcription_factors": ["TF"], "frequencies": []}}
    assert not (out / "regulatory" / "21" / "1000001-2000000.json.gz").exists()
    assert not (out / "regulatory" / "22").exists()
    info = json.loads((out / "info.json").read_text())
    assert info["regulatory"] is True and info["regulatory_schema"] == 1 and info["cell_types"] == ["A549", "HUVEC"]
    manifest = (out / "CHECKSUMS.sha256").read_text()
    assert "regulatory/21/1-1000000.json.gz" in manifest and "regulatory/binding_matrices.json" in manifest

    # A record reaching the far feature alone keeps that shard and an empty table: the runtime
    # opens the table of every regulatory tree.
    far_vcf = tmp_path / "far.vcf"
    far_vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t500100\t.\tA\tG\t.\t.\t.\n")
    far_out = tmp_path / "far_out"
    assert pjc.main(["--cache", str(cache), "--out", str(far_out), "--vcf", str(far_vcf), "--regulatory"]) == 0
    assert [f["stable_id"] for f in json.loads((far_out / "regulatory" / "21" / "1-1000000.json").read_text())] == ["ENSR_FAR"]
    assert json.loads((far_out / "regulatory" / "binding_matrices.json").read_text()) == {}

    # A --perl-info run carries the schema and cell-type list from the source's own info.json,
    # since Ensembl VEP's info.txt has neither; a source lacking them is refused.
    info_txt = tmp_path / "info.txt"
    info_txt.write_text("species\thomo_sapiens\nassembly\tGRCh38\nregulatory\t1\n")
    perl_out = tmp_path / "perl_out"
    assert pjc.main(["--cache", str(cache), "--out", str(perl_out), "--vcf", str(vcf), "--regulatory",
                     "--perl-info", str(info_txt), "--cache-version", "116"]) == 0
    perl_info = json.loads((perl_out / "info.json").read_text())
    assert perl_info["regulatory_schema"] == 1 and perl_info["cell_types"] == ["A549", "HUVEC"]
    bare = tmp_path / "bare_cache"
    shutil.copytree(cache, bare)
    (bare / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116, "regulatory": True}))
    assert pjc.main(["--cache", str(bare), "--out", str(tmp_path / "bare_out"), "--vcf", str(vcf), "--regulatory"]) == 2
    assert "lacks regulatory_schema, cell_types" in capsys.readouterr().err

    plain = tmp_path / "plain"
    rc = pjc.main(["--cache", str(cache), "--out", str(plain), "--vcf", str(vcf)])
    assert rc == 0
    assert not (plain / "regulatory").exists()
    info = json.loads((plain / "info.json").read_text())
    assert info["regulatory"] is False and "regulatory_schema" not in info and "cell_types" not in info

    rc = pjc.main(["--cache", str(tmp_path / "none"), "--out", str(tmp_path / "x"), "--vcf", str(vcf), "--regulatory",
                   "--assembly", "GRCh38", "--cache-version", "116"])
    assert rc == 2


def test_prune_writes_the_checksum_manifest_last(tmp_path: Path, capsys):
    """The output ends with `CHECKSUMS.sha256`: one `sha256sum` text-mode line per file of the
    output, the manifest itself excepted, paths relative to the output in byte order, so
    `sha256sum --check` accepts it; `--info-only` refreshes it after rewriting `info.json`, and
    `--manifest-only` writes it alone over an assembled cache directory."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(
        json.dumps([{"stable_id": "ENST_A", "start": "10000", "end": "20000", "dbID": "1"}]))
    (cache / "chr_synonyms.txt").write_text("21\tchr21\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t12000\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38",
                     "--cache-version", "116", "--gzip"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["manifest_files"] == 3
    checksums = out / "CHECKSUMS.sha256"
    lines = checksums.read_text().splitlines()
    assert [line[66:] for line in lines] == ["chr_synonyms.txt", "info.json", "transcripts/21/1-1000000.json.gz"]
    for line in lines:
        digest, path = line[:64], line[66:]
        assert line[64:66] == "  "
        assert hashlib.sha256((out / path).read_bytes()).hexdigest() == digest
    if shutil.which("sha256sum"):
        check = subprocess.run(["sha256sum", "--check", "--strict", "CHECKSUMS.sha256"], cwd=out, capture_output=True, text=True, check=False)
        assert check.returncode == 0, check.stdout + check.stderr
    before = checksums.read_text()
    assert pjc.main(["--info-only", "--out", str(out), "--assembly", "GRCh38", "--cache-version", "117"]) == 0
    after = checksums.read_text()
    assert after != before
    assert hashlib.sha256((out / "info.json").read_bytes()).hexdigest() in after
    (out / "transcripts" / "21" / "1000001-2000000.json").write_text("[]")
    assert pjc.main(["--manifest-only", "--out", str(out)]) == 0
    assert json.loads(capsys.readouterr().out) == {"manifest_files": 4}
    assert "transcripts/21/1000001-2000000.json" in checksums.read_text()
    assert pjc.main(["--manifest-only", "--out", str(tmp_path / "not-a-cache")]) == 2
    assert "holds no info.json" in capsys.readouterr().err


def test_manifest_lists_files_through_directory_links_in_path_byte_order(tmp_path: Path):
    """The writer lists a file reached through a symbolic link to a directory, as the Rust writer
    does, and orders the lines by the bytes of the paths alone: a name that is a proper prefix of
    another whose next byte sorts below a newline comes first, as it does in the Rust writer's
    byte order. A tree nested deeper than `MANIFEST_MAX_DEPTH`, a link to an ancestor, is refused
    with the tool's error instead of walked without end."""
    out = tmp_path / "out"
    (out / "transcripts" / "21").mkdir(parents=True)
    (out / "info.json").write_text("{}\n", encoding="utf-8")
    (out / "info.json\tcopy").write_text("{}\n", encoding="utf-8")
    (out / "transcripts" / "21" / "1-1000000.json").write_text("[]", encoding="utf-8")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (elsewhere / "1-1000000_var.gz").write_bytes(b"gz")
    (out / "variations").mkdir()
    os.symlink(elsewhere, out / "variations" / "21", target_is_directory=True)
    assert pjc.write_manifest(out) == 4
    lines = (out / "CHECKSUMS.sha256").read_text(encoding="utf-8").splitlines()
    assert [line[66:] for line in lines] == [
        "info.json", "info.json\tcopy", "transcripts/21/1-1000000.json", "variations/21/1-1000000_var.gz",
    ]
    assert lines[3][:64] == hashlib.sha256(b"gz").hexdigest()
    os.symlink(out, out / "transcripts" / "loop", target_is_directory=True)
    with pytest.raises(SystemExit) as refusal:
        pjc.write_manifest(out)
    assert str(refusal.value).startswith("ERROR: [prune_json_cache] ")
    assert str(refusal.value).endswith(f": directories nest more than {pjc.MANIFEST_MAX_DEPTH} levels deep")
    (out / "transcripts" / "loop").unlink()
    os.symlink(out / "absent", out / "dangling")
    with pytest.raises(SystemExit) as refusal:
        pjc.write_manifest(out)
    assert str(refusal.value) == f"ERROR: [prune_json_cache] {out / 'dangling'}: a symbolic link to nothing"


def test_prune_writes_contig_order_and_the_kept_chromosomes_synonyms(tmp_path: Path, capsys):
    """`contig_order` lists the output's own chromosome directories in the converter's order
    (autosomes by number, X, Y, MT, then the rest by name), whatever the source advertised, and
    `chr_synonyms.txt` keeps the lines naming a chromosome the output holds."""
    cache = tmp_path / "cache"
    for chrom in ["21", "X", "3", "KI270728.1", "MT", "22"]:
        (cache / "transcripts" / chrom).mkdir(parents=True)
        (cache / "transcripts" / chrom / "1-1000000.json").write_text(
            json.dumps([{"stable_id": f"ENST_{chrom}", "start": "100", "end": "150", "dbID": "1"}]))
    (cache / "info.json").write_text(json.dumps({
        "species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116, "source_versions": {},
        "contig_order": ["1", "2", "3", "21", "22", "X", "Y", "MT", "KI270728.1"]}))
    (cache / "chr_synonyms.txt").write_text(
        "21\tNC_000021.9\n21\tchr21\nX\tchrX\n22\tNC_000022.11\nMT\tchrM\nKI270728.1\tchr21_KI270728v1_random\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "X\t120\t.\tA\tG\t.\t.\t.\n21\t120\t.\tA\tG\t.\t.\t.\n3\t120\t.\tA\tG\t.\t.\t.\n"
                   "KI270728.1\t120\t.\tA\tG\t.\t.\t.\nMT\t120\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf)]) == 0
    summary = json.loads(capsys.readouterr().out)
    info = json.loads((out / "info.json").read_text())
    assert info["contig_order"] == ["3", "21", "X", "MT", "KI270728.1"]
    assert (out / "chr_synonyms.txt").read_text().splitlines() == [
        "21\tNC_000021.9", "21\tchr21", "X\tchrX", "MT\tchrM", "KI270728.1\tchr21_KI270728v1_random"]
    assert summary["synonym_lines"] == 5
    # A source without the file writes none, and an info-only output names no contig.
    (cache / "chr_synonyms.txt").unlink()
    bare = tmp_path / "bare"
    assert pjc.main(["--cache", str(cache), "--out", str(bare), "--vcf", str(vcf)]) == 0
    capsys.readouterr()
    assert not (bare / "chr_synonyms.txt").exists()
    info_only = tmp_path / "info_only"
    assert pjc.main(["--info-only", "--out", str(info_only), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    assert "contig_order" not in json.loads((info_only / "info.json").read_text())


def test_prune_resolves_input_chromosome_names_as_the_runtime_does(tmp_path: Path, capsys):
    """Records on `NC_000021.9` (a synonym the cache's file links), `chr21` and `CHR21` (the
    `chr` prefix removed, in any case, with no synonym line naming them), `Y` (the prefix added,
    for a directory the cache names `chrY`) and `M` (read as `MT`) reach the cache's `21`,
    `chrY` and `MT` directories; a name that resolves to nothing keeps nothing."""
    cache = tmp_path / "cache"
    for chrom in ["21", "chrY", "MT"]:
        (cache / "transcripts" / chrom).mkdir(parents=True)
        (cache / "transcripts" / chrom / "1-1000000.json").write_text(
            json.dumps([{"stable_id": f"ENST_{chrom}", "start": "100", "end": "150", "dbID": "1"},
                        {"stable_id": f"ENST_{chrom}_far", "start": "900000", "end": "900100", "dbID": "2"}]))
    (cache / "chr_synonyms.txt").write_text("21\tNC_000021.9\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "NC_000021.9\t120\t.\tA\tG\t.\t.\t.\nchr21\t130\t.\tA\tG\t.\t.\t.\nCHR21\t140\t.\tA\tG\t.\t.\t.\n"
                   "Y\t120\t.\tA\tG\t.\t.\t.\nM\t120\t.\tA\tG\t.\t.\t.\nNC_000022.11\t120\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["spans_by_chromosome"] == {"21": 3, "MT": 1, "NC_000022.11": 1, "chrY": 1}
    assert summary["transcripts_per_chromosome"] == {"21": 1, "MT": 1, "chrY": 1}
    assert [t["stable_id"] for t in json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())] == ["ENST_21"]


def test_prune_resolves_names_through_the_runs_synonyms_file(tmp_path: Path, capsys):
    """Under `--synonyms` the records resolve through that file alone, as Ensembl VEP reads a
    `--synonyms` file in place of the cache's `chr_synonyms.txt`: `MT` reaches the cache's
    `mitochondrion_genome` through the run's line and the cache's own `NC_000021.9` line does
    not resolve; the copied `chr_synonyms.txt` stays the cache's; a missing file is refused."""
    cache = tmp_path / "cache"
    for chrom in ["21", "mitochondrion_genome"]:
        (cache / "transcripts" / chrom).mkdir(parents=True)
        (cache / "transcripts" / chrom / "1-1000000.json").write_text(
            json.dumps([{"stable_id": f"T_{chrom}", "start": "100", "end": "150", "dbID": "1"}]))
    (cache / "chr_synonyms.txt").write_text("21\tNC_000021.9\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "MT\t120\t.\tA\tG\t.\t.\t.\nNC_000021.9\t120\t.\tA\tG\t.\t.\t.\n21\t130\t.\tA\tG\t.\t.\t.\n")
    synonyms = tmp_path / "synonyms.txt"
    synonyms.write_text("mitochondrion_genome\tMT\n")
    common = ["--cache", str(cache), "--vcf", str(vcf), "--assembly", "BDGP6.54", "--cache-version", "116"]
    assert pjc.main([*common, "--out", str(tmp_path / "plain")]) == 0
    plain = json.loads(capsys.readouterr().out)
    assert plain["spans_by_chromosome"] == {"21": 2, "MT": 1}
    assert plain["transcripts_per_chromosome"] == {"21": 1}
    out = tmp_path / "out"
    assert pjc.main([*common, "--out", str(out), "--synonyms", str(synonyms)]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["spans_by_chromosome"] == {"21": 1, "NC_000021.9": 1, "mitochondrion_genome": 1}
    assert summary["transcripts_per_chromosome"] == {"21": 1, "mitochondrion_genome": 1}
    assert (out / "chr_synonyms.txt").read_text() == "21\tNC_000021.9\n"
    assert pjc.main([*common, "--out", str(tmp_path / "out2"), "--synonyms", str(tmp_path / "absent.txt")]) == 2
    assert not (tmp_path / "out2" / "transcripts").exists()


def test_breakend_mate_positions_count_as_spans(tmp_path: Path, capsys):
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\t.\tN\tN]13:5000]\t.\t.\tSVTYPE=BND\n"
        "21\t200\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=22;POS2=700;END2=701\n"
        "21\t300\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;POS2=900;END2=901\n"
    )
    spans = pjc.record_spans(vcf)
    assert spans["13"] == [(5000, 5000)]
    assert spans["22"] == [(700, 701)]
    assert (300, 901) in spans["21"], "an intra-chromosomal pair keeps the span between the breakends"
    # The summary counts those spans, not records: three records, one mate span each, plus
    # the pair span of the same-chromosome breakend.
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out"), "--vcf", str(vcf), "--assembly", "GRCh37", "--cache-version", "116"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["spans_by_chromosome"] == {"13": 1, "21": 5, "22": 1}
    assert summary["transcripts_per_chromosome"] == {"21": 1} and summary["shards"] == 1


def test_prune_reads_a_bgz_vcf_and_refuses_a_missing_one(tmp_path: Path, capsys):
    """A `.bgz` input (BGZF is gzip) yields the spans of the same records stored plain or as
    `.gz`, and a `--vcf` path that is not a file is refused by name with status 2."""
    records = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "21\t24000\tdel1\tT\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=25000\n"
        "21\t300\tbnd1\tN\tN]22:5000]\t.\tPASS\tSVTYPE=BND\n"
    ).encode("utf-8")
    plain = tmp_path / "variants.vcf"
    plain.write_bytes(records)
    gz = tmp_path / "variants.vcf.gz"
    with gzip.open(gz, "wb") as handle:
        handle.write(records)
    bgz = tmp_path / "variants.vcf.bgz"
    with bgz.open("wb") as handle:
        writer = pjc.BgzfWriter(handle)
        writer.write(records)
        writer.close()
    spans = pjc.record_spans(plain)
    assert spans == {"21": [(100, 100), (24000, 25000), (300, 300)], "22": [(5000, 5000)]}
    assert pjc.record_spans(gz) == spans and pjc.record_spans(bgz) == spans
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    common = ["--cache", str(cache), "--out", str(tmp_path / "out"), "--assembly", "GRCh38", "--cache-version", "116"]
    assert pjc.main([*common, "--vcf", str(bgz)]) == 0
    assert json.loads(capsys.readouterr().out)["spans_by_chromosome"] == {c: len(s) for c, s in sorted(spans.items())}
    missing = tmp_path / "absent.vcf.bgz"
    assert pjc.main([*common, "--vcf", str(bgz), "--vcf", str(missing)]) == 2
    assert capsys.readouterr().err == f"ERROR: [prune_json_cache] --vcf {missing} is not a file\n"
    assert pjc.main([*common, "--vcf", str(tmp_path)]) == 2
    assert capsys.readouterr().err == f"ERROR: [prune_json_cache] --vcf {tmp_path} is not a file\n"


def test_info_json_from_perl_info_txt(tmp_path: Path, capsys):
    info_txt = tmp_path / "info.txt"
    info_txt.write_text(
        "species\thomo_sapiens\nassembly\tGRCh37\nsift\tb\npolyphen\tb\n"
        "source_gencode\tGENCODE 19\nsource_assembly\tGRCh37.p13\nsource_1000genomes\tphase3\n"
        "variation_cols\tchr,variation_name,failed\nvar_type\ttabix\nregulatory\t1\nbuild\t-\n"
        "cell_types\tA549_(m,_58_y),B_cell_(f)\n"
    )
    out = tmp_path / "out"
    rc = pjc.main(["--info-only", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"])
    assert rc == 0
    info = json.loads((out / "info.json").read_text())
    assert info["species"] == "homo_sapiens" and info["assembly"] == "GRCh37" and info["cache_version"] == 115
    assert info["source_versions"] == {"gencode": "GENCODE 19", "assembly": "GRCh37.p13", "1000genomes": "phase3"}
    # No variation slice was written, so the keys that advertise one are dropped.
    assert "variation_cols" not in info and "var_type" not in info
    assert "cell_types" not in info
    # The flag follows the output's own tree, not the source's line: no shard, no regulatory data.
    assert info["regulatory"] is False and "build" not in info
    assert not (out / "transcripts").exists()
    (out / "regulatory" / "21").mkdir(parents=True)
    (out / "regulatory" / "21" / "1-1000000.json").write_text("[]")
    # A tree needs its schema and cell-type list, which no info.txt carries: a refresh keeps
    # the ones the metadata being replaced holds and refuses when neither side has them.
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"]) == 2
    assert "carries regulatory_schema, cell_types" in capsys.readouterr().err
    (out / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh37", "cache_version": 115,
                                               "regulatory": True, "regulatory_schema": 1, "cell_types": ["A549"]}))
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"]) == 0
    refreshed = json.loads((out / "info.json").read_text())
    assert refreshed["regulatory"] is True and refreshed["regulatory_schema"] == 1 and refreshed["cell_types"] == ["A549"]
    # Each missing key is filled alone: an info.txt carrying the schema keeps it as the whole
    # number the runtime reads, and the cell-type list comes from the metadata being replaced.
    schema_txt = tmp_path / "schema_info.txt"
    schema_txt.write_text(info_txt.read_text() + "regulatory_schema\t1\n")
    (out / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh37", "cache_version": 115,
                                               "regulatory": True, "cell_types": ["A549"]}))
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(schema_txt), "--cache-version", "115"]) == 0
    mixed = json.loads((out / "info.json").read_text())
    assert mixed["regulatory"] is True and mixed["regulatory_schema"] == 1 and mixed["cell_types"] == ["A549"]
    assert isinstance(mixed["regulatory_schema"], int)
    bad_schema_txt = tmp_path / "bad_schema_info.txt"
    bad_schema_txt.write_text(info_txt.read_text() + "regulatory_schema\tone\n")
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(bad_schema_txt), "--cache-version", "115"]) == 2
    assert capsys.readouterr().err == f"ERROR: [prune_json_cache] --perl-info {bad_schema_txt}: regulatory_schema 'one' is not a whole number\n"
    # --regulatory names a slice to cut, which an --info-only run does not.
    assert pjc.main(["--info-only", "--regulatory", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"]) == 2


VAR_COLS = ["chr", "variation_name", "failed", "somatic", "start", "end", "allele_string", "strand", "clin_sig",
            "phenotype_or_disease", "clin_sig_allele", "clin_sig_ref_allele", "clinical_impact", "pubmed",
            "var_synonyms", "AF", "gnomADg"]


def _all_vars_line(chrom: str, name: str, start: int, alleles: str = "A/G") -> bytes:
    fields = [chrom, name, ".", ".", str(start), ".", alleles] + ["."] * (len(VAR_COLS) - 7)
    return ("\t".join(fields) + "\n").encode()


def _write_all_vars(path: Path, lines: list[bytes]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        writer = pjc.BgzfWriter(raw)
        for line in lines:
            writer.write(line)
        writer.close()


def _bgzf_blocks(path: Path) -> dict[int, bytes]:
    """Decompressed bgzf blocks by their file offset."""
    data = path.read_bytes()
    blocks, offset = {}, 0
    while offset < len(data):
        assert data[offset:offset + 4] == b"\x1f\x8b\x08\x04", "a bgzf member with the extra field"
        block_size = int.from_bytes(data[offset + 16:offset + 18], "little") + 1
        payload = data[offset + 18:offset + block_size - 8]
        blocks[offset] = zlib.decompress(payload, -15)
        offset += block_size
    return blocks


def _read_chunk(blocks: dict[int, bytes], begin: int, end: int) -> bytes:
    offsets = sorted(blocks)
    out = b""
    for block_offset in offsets:
        if block_offset < begin >> 16 or block_offset > end >> 16:
            continue
        block = blocks[block_offset]
        lo = (begin & 0xFFFF) if block_offset == begin >> 16 else 0
        hi = (end & 0xFFFF) if block_offset == end >> 16 else len(block)
        out += block[lo:hi]
    return out


def _parse_csi(path: Path) -> tuple[dict, dict[int, list[tuple[int, int, int]]]]:
    """The header fields and, for the single reference, bin -> [(loffset, begin, end), ...]."""
    data = gzip.open(path, "rb").read()
    assert data[:4] == b"CSI\x01"
    min_shift, depth, l_aux = struct.unpack("<iii", data[4:16])
    aux = data[16:16 + l_aux]
    fmt, seq_col, beg_col, end_col, meta, skip, l_nm = struct.unpack("<7i", aux[:28])
    names = aux[28:28 + l_nm].split(b"\0")[:-1]
    offset = 16 + l_aux
    n_ref, = struct.unpack("<i", data[offset:offset + 4])
    offset += 4
    header = {"min_shift": min_shift, "depth": depth, "format": fmt, "seq_col": seq_col, "beg_col": beg_col,
              "end_col": end_col, "meta": meta, "skip": skip, "names": names, "n_ref": n_ref}
    bins: dict[int, list[tuple[int, int, int]]] = {}
    for _ in range(n_ref):
        n_bin, = struct.unpack("<i", data[offset:offset + 4])
        offset += 4
        for _ in range(n_bin):
            b, loffset, n_chunk = struct.unpack("<IQi", data[offset:offset + 16])
            offset += 16
            for _ in range(n_chunk):
                begin, end = struct.unpack("<QQ", data[offset:offset + 16])
                offset += 16
                bins.setdefault(b, []).append((loffset, begin, end))
    return header, bins


def _reg2bins(beg: int, end: int, min_shift: int, depth: int) -> list[int]:
    found = []
    for level in range(depth + 1):
        shift = min_shift + 3 * (depth - level)
        first = ((1 << (3 * level)) - 1) // 7
        found.extend(range(first + (beg >> shift), first + ((end - 1) >> shift) + 1))
    return found


def _csi_query(gz: Path, start: int) -> list[str]:
    """The record names an index-driven reader returns for a one-base query at `start`."""
    header, bins = _parse_csi(Path(str(gz) + ".csi"))
    blocks = _bgzf_blocks(gz)
    names = []
    for b in _reg2bins(start - 1, start, header["min_shift"], header["depth"]):
        for _loffset, begin, end in bins.get(b, []):
            for line in _read_chunk(blocks, begin, end).split(b"\n"):
                fields = line.split(b"\t")
                if len(fields) > 4 and int(fields[header["beg_col"] - 1]) == start:
                    names.append(fields[1].decode())
    return names


def test_prune_variation_slice_reindexes(tmp_path: Path):
    src = tmp_path / "perl-cache"
    lines = [_all_vars_line("21", "rs1", 100), _all_vars_line("21", "rs2", 150), _all_vars_line("21", "somatic1", 200, "COSMIC_MUTATION"),
             _all_vars_line("21", "rs3", 200), _all_vars_line("21", "rs3b", 201), _all_vars_line("21", "inherited1", 201, "HGMD_MUTATION"),
             _all_vars_line("21", "rs4", 202),
             _all_vars_line("21", "rs5", 5_000_000), _all_vars_line("21", "rs6", 5_000_001), _all_vars_line("21", "rs7", 5_000_003),
             _all_vars_line("21", "rs8", 9_000_000)]
    # A long record spans two bgzf blocks, as a record with many synonyms does.
    long_record = ("\t".join(["21", "rs9", ".", ".", "9000001", ".", "A/T"] + ["."] * 7 + ["x" * 70_000] + ["."] * 2) + "\n").encode()
    lines.append(long_record)
    lines.append(_all_vars_line("21", "rs10", 9_000_002))
    _write_all_vars(src / "21" / "all_vars.gz", lines)
    _write_all_vars(src / "22" / "all_vars.gz", [_all_vars_line("22", "rs22", 100)])
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    info_txt = tmp_path / "info.txt"
    info_txt.write_text("species\thomo_sapiens\nassembly\tGRCh38\nvariation_cols\t" + ",".join(VAR_COLS) + "\nvar_type\ttabix\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t200\t.\tA\tG\t.\t.\t.\n"
        "21\t5000002\t.\tA\tG\t.\t.\t.\n"
        "21\t9000001\t.\tA\tT\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116",
                   "--variations", str(src), "--exclude-variation", r"^somatic\d", "--exclude-variation", r"^inherited\d"])
    assert rc == 0
    gz = out / "variations" / "21" / "all_vars.gz"
    kept = [line.split(b"\t")[1].decode() for line in gzip.open(gz, "rb").read().splitlines()]
    # Starts within one base of a record span survive (199..201, 5000001..5000003, 9000000..9000002),
    # the two records the patterns name are dropped, chromosome 22 has no record span.
    assert kept == ["rs3", "rs3b", "rs6", "rs7", "rs8", "rs9", "rs10"]
    assert not (out / "variations" / "22").exists()
    header, bins = _parse_csi(Path(str(gz) + ".csi"))
    assert header["names"] == [b"21"] and header["n_ref"] == 1
    assert (header["seq_col"], header["beg_col"], header["end_col"]) == (1, 5, 5), "the index names the start column as both bounds"
    assert (header["format"], header["meta"], header["skip"]) == (0, ord("#"), 0)
    assert (header["min_shift"], header["depth"]) == (pjc.CSI_MIN_SHIFT, pjc.CSI_DEPTH)
    meta_bin = pjc.bin_first(pjc.CSI_DEPTH + 1) + 1
    assert bins[meta_bin][1][1:] == (7, 0), "the pseudo-bin counts the records written"
    # The first bin's `loffset` is the first record's virtual offset, 0: htslib's `hts_itr_query`
    # starts each chunk at the `loffset` of the bin covering the query start, so a later offset
    # hides the first record from a query at its position.
    first_bin = pjc.bin_first(pjc.CSI_DEPTH) + (199 >> pjc.CSI_MIN_SHIFT)
    assert bins[first_bin][0][0] == 0
    for start, expected in [(200, ["rs3"]), (201, ["rs3b"]), (202, []), (5_000_001, ["rs6"]), (5_000_003, ["rs7"]), (9_000_001, ["rs9"]),
                            (9_000_002, ["rs10"]), (100, []), (150, [])]:
        assert _csi_query(gz, start) == expected, start
    info = json.loads((out / "info.json").read_text())
    assert info["variation_cols"] == VAR_COLS and info["var_type"] == "tabix"
    # Without --variations the same run writes no slice and no variation keys.
    out2 = tmp_path / "out2"
    assert pjc.main(["--cache", str(cache), "--out", str(out2), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116"]) == 0
    info2 = json.loads((out2 / "info.json").read_text())
    assert "variation_cols" not in info2 and "var_type" not in info2 and not (out2 / "variations").exists()
    # A --variations directory without a file for any record chromosome, and a negative flank,
    # are refused rather than written as an empty slice.
    empty = tmp_path / "empty-variations"
    (empty / "22").mkdir(parents=True)
    out3 = tmp_path / "out3"
    assert pjc.main(["--cache", str(cache), "--out", str(out3), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116",
                     "--variations", str(empty)]) == 2
    assert not (out3 / "info.json").exists() and not (out3 / "transcripts").exists()
    assert pjc.main(["--cache", str(cache), "--out", str(out3), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116",
                     "--variations", str(src), "--variation-flank", "-1"]) == 2


def test_prune_variation_slice_takes_the_columns_of_the_variation_cache(tmp_path: Path):
    """A slice cut from a native cache beside an Ensembl cache's metadata declares the native
    cache's columns: `--variation-info` replaces `variation_cols` and `var_type` of the info being
    written, and the index is built on the start column of those columns. The option without
    `--variations`, either option under `--info-only`, an info.json without a `variation_cols`
    list, or one whose top level is not an object, is refused; so is a slice whose records carry
    a field count other than the declared columns, or whose file is not gzip or does not inflate,
    before anything is written, while a chromosome file with no record is written as the empty
    indexed slice the Ensembl VEP dumper writes, not refused."""
    native = tmp_path / "native"
    native_cols = VAR_COLS[:15] + ["AF"]
    line = ("\t".join(["21", "rs1", ".", ".", "200", ".", "A/G"] + ["."] * (len(native_cols) - 7)) + "\n").encode()
    _write_all_vars(native / "variations" / "21" / "all_vars.gz", [line])
    (native / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116,
                                                   "variation_cols": native_cols, "var_type": "tabix", "source_versions": {}}))
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    info_txt = tmp_path / "info.txt"
    info_txt.write_text("species\thomo_sapiens\nassembly\tGRCh38\nsource_dbSNP\t156\nvariation_cols\t" + ",".join(VAR_COLS) + "\nvar_type\ttabix\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t200\t.\tA\tG\t.\t.\t.\n")
    common = ["--cache", str(cache), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116"]
    native_info = ["--variation-info", str(native / "info.json")]
    out = tmp_path / "out"
    assert pjc.main([*common, "--out", str(out), "--variations", str(native / "variations"), *native_info]) == 0
    info = json.loads((out / "info.json").read_text())
    assert info["variation_cols"] == native_cols and info["var_type"] == "tabix"
    assert info["source_versions"] == {"dbSNP": "156"}, "the rest of the metadata is still the Ensembl cache's"
    gz = out / "variations" / "21" / "all_vars.gz"
    assert _csi_query(gz, 200) == ["rs1"]
    assert pjc.main([*common, "--out", str(tmp_path / "out2"), *native_info]) == 2
    (tmp_path / "bare.json").write_text(json.dumps({"species": "homo_sapiens"}))
    assert pjc.main([*common, "--out", str(tmp_path / "out3"), "--variations", str(native / "variations"), "--variation-info", str(tmp_path / "bare.json")]) == 2
    assert not (tmp_path / "out3" / "info.json").exists()
    (tmp_path / "list.json").write_text(json.dumps([]))
    assert pjc.main([*common, "--out", str(tmp_path / "out4"), "--variations", str(native / "variations"), "--variation-info", str(tmp_path / "list.json")]) == 2
    assert not (tmp_path / "out4" / "info.json").exists()
    # The Ensembl cache's columns over the native files: the record width disagrees, so the
    # declaration is refused rather than written over records it does not describe.
    assert pjc.main([*common, "--out", str(tmp_path / "out5"), "--variations", str(native / "variations")]) == 2
    assert not (tmp_path / "out5").exists() or not any((tmp_path / "out5").iterdir())
    assert pjc.main(["--info-only", "--out", str(tmp_path / "out6"), "--perl-info", str(info_txt), "--cache-version", "116", *native_info]) == 2
    assert not (tmp_path / "out6" / "info.json").exists()
    assert pjc.main(["--info-only", "--out", str(tmp_path / "out7"), "--perl-info", str(info_txt), "--cache-version", "116",
                     "--variations", str(native / "variations"), *native_info]) == 2
    assert not (tmp_path / "out7" / "info.json").exists()
    # A chromosome file with no record yields the empty bgzf file and an index naming the sequence
    # with no bin, and the metadata keeps the variation keys; a file that is not gzip is refused.
    empty = tmp_path / "empty"
    _write_all_vars(empty / "variations" / "21" / "all_vars.gz", [])
    assert pjc.main([*common, "--out", str(tmp_path / "out8"), "--variations", str(empty / "variations"), *native_info]) == 0
    empty_slice = tmp_path / "out8" / "variations" / "21" / "all_vars.gz"
    assert empty_slice.read_bytes() == pjc.BGZF_EOF
    header, bins = _parse_csi(Path(str(empty_slice) + ".csi"))
    assert header["names"] == [b"21"] and header["n_ref"] == 1 and bins == {}
    assert _csi_query(empty_slice, 200) == []
    info8 = json.loads((tmp_path / "out8" / "info.json").read_text())
    assert info8["variation_cols"] == native_cols and info8["var_type"] == "tabix"
    plain = tmp_path / "plain"
    (plain / "variations" / "21").mkdir(parents=True)
    (plain / "variations" / "21" / "all_vars.gz").write_bytes(b"21\trs1\t.\t.\t200\n")
    assert pjc.main([*common, "--out", str(tmp_path / "out9"), "--variations", str(plain / "variations"), *native_info]) == 2
    assert not (tmp_path / "out9" / "transcripts").exists()
    # A gzip header over a corrupt deflate stream raises zlib.error, not OSError or EOFError.
    corrupt = tmp_path / "corrupt"
    (corrupt / "variations" / "21").mkdir(parents=True)
    (corrupt / "variations" / "21" / "all_vars.gz").write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03" + b"\xff\xff\xff\xff\x00garbage" * 4)
    assert pjc.main([*common, "--out", str(tmp_path / "out10"), "--variations", str(corrupt / "variations"), *native_info]) == 2
    assert not (tmp_path / "out10" / "transcripts").exists()


def test_bgzf_writer_blocks_and_virtual_offsets(tmp_path: Path):
    path = tmp_path / "x.gz"
    with path.open("wb") as raw:
        writer = pjc.BgzfWriter(raw)
        assert writer.tell() == 0
        writer.write(b"a" * 1000)
        assert writer.tell() == 1000
        writer.write(b"b" * pjc.BGZF_BLOCK_BYTES)
        first_block_size = writer.block_offset
        assert first_block_size > 0 and writer.tell() == (first_block_size << 16) | 1000
        writer.close()
    blocks = _bgzf_blocks(path)
    assert [len(b) for b in blocks.values()] == [pjc.BGZF_BLOCK_BYTES, 1000, 0]
    assert gzip.open(path, "rb").read() == b"a" * 1000 + b"b" * pjc.BGZF_BLOCK_BYTES
    assert path.read_bytes().endswith(pjc.BGZF_EOF)


def test_reg2bin_follows_the_binning_scheme():
    # Bins at the finest level start at (8^8 - 1) / 7 and cover 16 kb windows; a span across a
    # window boundary moves up one level.
    finest = pjc.bin_first(pjc.CSI_DEPTH)
    assert pjc.reg2bin(0, 1) == finest
    assert pjc.reg2bin(16_384, 16_385) == finest + 1
    assert pjc.reg2bin(16_383, 16_385) == pjc.bin_first(pjc.CSI_DEPTH - 1)
    assert pjc.bin_level(0) == 0 and pjc.bin_level(finest) == pjc.CSI_DEPTH
    assert pjc.bin_bottom(finest + 3) == 3 and pjc.bin_bottom(pjc.bin_first(pjc.CSI_DEPTH - 1) + 1) == 8
    assert pjc.variation_windows([(100, 100), (103, 110), (300, 300)], 1) == [(99, 111), (299, 301)]
    assert pjc.variation_windows([(1, 5)], 3) == [(1, 8)]


def test_gzip_option_writes_gz_shards_readable_as_source(tmp_path: Path):
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "10000", "end": "20000"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t12000\t.\tA\tG\t.\t.\t.\n")
    out1 = tmp_path / "out1"
    assert pjc.main(["--cache", str(cache), "--out", str(out1), "--vcf", str(vcf), "--gzip", "--assembly", "GRCh37", "--cache-version", "115"]) == 0
    gz = out1 / "transcripts" / "21" / "1-1000000.json.gz"
    assert gz.is_file()
    assert pjc.read_shard(gz)[0]["stable_id"] == "ENST_A"
    # A gzipped cache is itself a valid pruning source.
    out2 = tmp_path / "out2"
    assert pjc.main(["--cache", str(out1), "--out", str(out2), "--vcf", str(vcf), "--assembly", "GRCh37", "--cache-version", "115"]) == 0
    assert json.loads((out2 / "transcripts" / "21" / "1-1000000.json").read_text())[0]["stable_id"] == "ENST_A"


def test_peptide_is_kept_and_cannot_be_dropped(tmp_path: Path):
    """The reference peptide's sequence edits are read from the cached peptide, so a pruned
    cache keeps it and `--drop peptide` is refused."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    tx = {"stable_id": "ENST_A", "start": "10000", "end": "20000",
          "variation_effect_feature_cache": {"peptide": "MLW", "sorted_exons": [], "translateable_seq": "GTGTTATGG"}}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([tx]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t12000\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--drop", "sorted_exons", "--assembly", "GRCh37", "--cache-version", "115"]) == 0
    kept = json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())[0]["variation_effect_feature_cache"]
    assert kept["peptide"] == "MLW" and "sorted_exons" not in kept
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out2"), "--vcf", str(vcf), "--drop", "peptide", "--assembly", "GRCh37", "--cache-version", "115"]) == 2


def test_drop_except_keeps_the_dropped_keys_on_the_named_transcripts(tmp_path: Path, capsys):
    """`--drop-except <file>` names the transcripts that keep every `--drop` key, so a corpus
    carries protein features and prediction matrices for the transcripts its records test and
    for no other; the option needs a `--drop` key to apply to, and a file that cannot be read
    or decoded is refused with the path and the reason, not a traceback."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    features = [{"start": 1, "end": 9, "hseqname": "PF00080", "analysis": "Pfam"}]
    matrices = {"sift": {"analysis": "sift", "predictions_data": "H4sI"}}
    vefc = {"peptide": "MLW", "sorted_exons": [], "protein_features": features, "protein_function_predictions": matrices}
    transcripts = [
        {"stable_id": "ENST_KEPT", "start": "10000", "end": "20000", "variation_effect_feature_cache": vefc},
        {"stable_id": "ENST_DROPPED", "start": "12000", "end": "22000", "variation_effect_feature_cache": vefc},
    ]
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps(transcripts))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t15000\t.\tA\tG\t.\t.\t.\n")
    keep = tmp_path / "keep.txt"
    keep.write_text("ENST_KEPT\n\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--drop", "sorted_exons",
                     "--drop", "protein_features", "--drop", "protein_function_predictions",
                     "--drop-except", str(keep), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    assert json.loads(capsys.readouterr().out)["drop_except"] == 1
    by_id = {t["stable_id"]: t["variation_effect_feature_cache"]
             for t in json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())}
    assert by_id["ENST_KEPT"]["protein_features"] == features
    assert by_id["ENST_KEPT"]["protein_function_predictions"] == matrices
    assert by_id["ENST_KEPT"]["sorted_exons"] == []
    assert set(by_id["ENST_DROPPED"]) == {"peptide"}
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out2"), "--vcf", str(vcf),
                     "--drop-except", str(keep), "--assembly", "GRCh38", "--cache-version", "116"]) == 2
    assert "--drop-except needs at least one --drop key" in capsys.readouterr().err
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out3"), "--vcf", str(vcf), "--drop", "sorted_exons",
                     "--drop-except", str(tmp_path / "missing.txt"), "--assembly", "GRCh38", "--cache-version", "116"]) == 2
    assert f"--drop-except {tmp_path / 'missing.txt'}: " in capsys.readouterr().err
    not_utf8 = tmp_path / "latin1.txt"
    not_utf8.write_bytes(b"\xffENST_KEPT\n")
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out4"), "--vcf", str(vcf), "--drop", "sorted_exons",
                     "--drop-except", str(not_utf8), "--assembly", "GRCh38", "--cache-version", "116"]) == 2
    assert f"--drop-except {not_utf8}: 'utf-8' codec can't decode" in capsys.readouterr().err
