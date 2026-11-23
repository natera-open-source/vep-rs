"""Tests for the synthetic plugin data file generator on hand-built inputs."""
from __future__ import annotations

import gzip
import json
import shutil
import struct
import sys
import zlib
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_plugin_fixture as mpf  # noqa: E402

needs_tools = pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix write the fixture",
)

VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "21\t100\tmissense\tA\tG\t.\t.\t.",
    "21\t200\tmulti\tC\tT,A\t.\t.\t.",
    "21\t300\tdel\tAT\tA\t.\t.\t.",
    "21\t400\tins\tA\tAT\t.\t.\t.",
    "21\t500\tintron\tg\ta\t.\t.\t.",
    "21\t600\tmissed\tG\tA\t.\t.\t.",
    "MT\t3308\tmito\tT\tC\t.\t.\t.",
]) + "\n"

# A default-format output of VCF: two transcripts at 21:100 with the same change, one with
# another change and one intronic; the multi-allelic site with one change per allele; a
# stop gained written `*`; the deletion at one position under the `-` allele and the
# insertion as a range; the intron site with no coding row and the mitochondrial site with a
# missense row.
REFERENCE = "\n".join([
    "## ENSEMBL VARIANT EFFECT PREDICTOR v116.2",
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
    "missense\t21:100\tG\tENSG1\tENST000002\tTranscript\tmissense_variant\t10\t10\t4\tR/G\tAga/Gga\t-\tIMPACT=MODERATE;STRAND=1",
    "missense\t21:100\tG\tENSG1\tENST000001\tTranscript\tmissense_variant\t10\t10\t4\tR/G\tAga/Gga\t-\tIMPACT=MODERATE;STRAND=1",
    "missense\t21:100\tG\tENSG1\tENST000003\tTranscript\tmissense_variant\t13\t13\t5\tK/E\tAag/Gag\t-\tIMPACT=MODERATE;STRAND=1",
    "missense\t21:100\tG\tENSG1\tENST000004\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1",
    "missense\t21:100\tG\tENSR1\tENSR1\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    "multi\t21:200\tT\tENSG2\tENST000005\tTranscript\tstop_gained\t20\t20\t7\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1",
    "multi\t21:200\tA\tENSG2\tENST000005\tTranscript\tsynonymous_variant\t20\t20\t7\tQ\tCag/Aag\t-\tIMPACT=LOW;STRAND=1",
    "del\t21:301\t-\tENSG2\tENST000005\tTranscript\tframeshift_variant\t30\t30\t10\tX\tAcg/cg\t-\tIMPACT=HIGH;STRAND=1",
    "ins\t21:400-401\tT\tENSG2\tENST000005\tTranscript\tframeshift_variant\t40-41\t40-41\t14\tX\tAcg/ATcg\t-\tIMPACT=HIGH;STRAND=1",
    "intron\t21:500\tA\tENSG2\tENST000005\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1",
    "missed\t21:600\tA\tENSG2\tENST000005\tTranscript\tmissense_variant\t60\t60\t20\tA/T\tGcc/Acc\t-\tIMPACT=MODERATE;STRAND=1",
    "mito\tMT:3308\tC\tENSG3\tENST000006\tTranscript\tmissense_variant\t2\t2\t1\tM/T\taTa/aCa\t-\tIMPACT=MODERATE;STRAND=1",
]) + "\n"


def write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(VCF, encoding="utf-8")
    reference = tmp_path / "default.txt.gz"
    with gzip.open(reference, "wt", encoding="utf-8") as handle:
        handle.write(REFERENCE)
    return vcf, reference


def rows_of(path: Path) -> list[list[str]]:
    """The data rows of a written fixture, the header left out, cells split on tabs with the
    line end stripped."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        lines = handle.read().split("\n")
    rows = [line.rstrip("\r").split("\t") for line in lines if line]
    return rows[1:]


def header_of(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return handle.readline()


def tbi_columns(path: Path) -> tuple[int, int, int, str]:
    """The 1-based sequence, start and end columns and the comment prefix of a `.tbi`."""
    raw = zlib.decompress(Path(f"{path}.tbi").read_bytes(), 16 + zlib.MAX_WBITS)
    magic, _n_ref, _format, col_seq, col_beg, col_end, meta = struct.unpack("<4siiiiii", raw[:28])
    assert magic == b"TBI\x01"
    return col_seq, col_beg, col_end, chr(meta)


def run(tmp_path: Path, layout: str, assembly: str, name: str, *extra: str) -> Path:
    vcf, reference = write_inputs(tmp_path)
    out = tmp_path / "plugin_data" / name
    mpf.main(["--layout", layout, "--assembly", assembly, "--seed", "38", "--vcf", str(vcf), "--reference", str(reference), "--out", str(out), "--force", *extra])
    return out


def test_header_order_per_layout():
    """The dbnsfp header the generator writes holds 458 columns,
    opening with the coordinates, alleles and amino acids and carrying the hg19 pair
    at columns 8 and 9; the dbscsnv header it writes holds 18 columns with the hg38
    pair at columns 5 and 6 and the two scores last."""
    assert len(mpf.DBNSFP_COLUMNS) == 458
    assert mpf.DBNSFP_COLUMNS[:6] == ["#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt"]
    assert mpf.DBNSFP_COLUMNS[7:9] == ["hg19_chr", "hg19_pos(1-based)"]
    assert mpf.DBNSFP_COLUMNS[13:15] == ["Ensembl_geneid", "Ensembl_transcriptid"]
    assert mpf.DBNSFP_COLUMNS.index("SIFT_score") == 37
    assert mpf.DBNSFP_COLUMNS.index("Polyphen2_HDIV_score") == 43
    assert mpf.DBNSFP_COLUMNS.index("Interpro_domain") == 448
    assert mpf.DBNSFP_COLUMNS[-1] == "Geuvadis_eQTL_target_gene"
    assert len(set(mpf.DBNSFP_COLUMNS)) == 458
    assert mpf.DBSCSNV_COLUMNS == [
        "chr", "pos", "ref", "alt", "hg38_chr", "hg38_pos", "RefSeq?", "Ensembl?", "RefSeq_region",
        "RefSeq_gene", "RefSeq_functional_consequence", "RefSeq_id_c.change_p.change", "Ensembl_region",
        "Ensembl_gene", "Ensembl_functional_consequence", "Ensembl_id_c.change_p.change", "ada_score", "rf_score",
    ]
    for column in (*mpf.TRANSCRIPT_FILLERS, *mpf.SINGLE_FILLERS):
        assert column in mpf.DBNSFP_COLUMNS, column


def test_tabix_index_column_per_layout():
    """Each layout is indexed as Ensembl VEP's plugin documents it (dbNSFP.pm:62-74,
    dbscSNV.pm:49-59): dbNSFP on
    columns 1 and 2 for GRCh38 and on 8 and 9 for GRCh37, dbscSNV on 1 and 2 for GRCh37 and
    on 5 and 6 for GRCh38 with `c` as the comment character. The generator ends a dbscSNV
    line in CRLF and a dbNSFP line in LF. A name that is no layout is refused."""
    assert mpf.preparation("dbnsfp", "GRCh38").tabix == ("-s", "1", "-b", "2", "-e", "2")
    assert mpf.preparation("dbnsfp", "GRCh37").tabix == ("-s", "8", "-b", "9", "-e", "9")
    assert mpf.preparation("dbnsfp", "GRCh37").line_end == "\n"
    assert mpf.preparation("dbscsnv", "GRCh37").tabix == ("-s", "1", "-b", "2", "-e", "2", "-c", "c")
    assert mpf.preparation("dbscsnv", "GRCh38").tabix == ("-s", "5", "-b", "6", "-e", "6", "-c", "c")
    assert mpf.preparation("dbscsnv", "GRCh38").line_end == "\r\n"
    with pytest.raises(SystemExit, match="unknown layout cadd_sv"):
        mpf.preparation("cadd_sv", "GRCh38")


def test_sites_and_reference_are_read_as_the_plugins_match(tmp_path: Path):
    """Only single-base records are sites (one per alternate base of a multi-allelic record,
    upper-cased), and the output's transcript rows are keyed by chromosome, position and
    allele with the deletion's `-` allele row and the insertion's range location left out; a
    plain output reads as a gzipped one does."""
    vcf, _ = write_inputs(tmp_path)
    sites = mpf.read_sites(vcf)
    assert [(s.chrom, s.pos, s.ref, s.alts) for s in sites] == [
        ("21", 100, "A", ("G",)),
        ("21", 200, "C", ("T", "A")),
        ("21", 500, "G", ("A",)),
        ("21", 600, "G", ("A",)),
        ("MT", 3308, "T", ("C",)),
    ]
    reference = tmp_path / "default.txt"
    reference.write_text(REFERENCE, encoding="utf-8")
    rows = mpf.read_reference(reference)
    assert sorted(rows) == [("21", 100, "G"), ("21", 200, "A"), ("21", 200, "T"), ("21", 500, "A"), ("21", 600, "A"), ("MT", 3308, "C")]
    assert [r.feature for r in rows[("21", 100, "G")]] == ["ENST000002", "ENST000001", "ENST000003", "ENST000004"]
    assert rows[("21", 200, "T")][0].amino_acids == "Q/*"


@needs_tools
def test_dbnsfp_rows_follow_the_output_and_the_preparation(tmp_path: Path):
    """At an input allele with a coding change the file carries one row per distinct change
    listing that change's transcripts sorted, the stop codon written `X`; the other two bases
    get one fabricated row each listing every transcript of the site; the GRCh38 file carries
    the corpus coordinates in columns 1 and 2 and a shifted position in the hg19 columns, `MT`
    is spelled `M`, and the planted row sits inside the detection region with `;` in every
    per-transcript column."""
    out = run(tmp_path, "dbnsfp", "GRCh38", "dbNSFP4.9a_grch38_synthetic.gz")
    assert header_of(out) == "\t".join(mpf.DBNSFP_COLUMNS) + "\n"
    rows = rows_of(out)
    assert tbi_columns(out) == (1, 2, 2, "#")
    columns = {name: i for i, name in enumerate(mpf.DBNSFP_COLUMNS)}
    by_key = {(r[0], r[1], r[2], r[3], r[4], r[5]): r for r in rows}
    own = by_key[("21", "100", "A", "G", "R", "G")]
    assert own[columns["Ensembl_transcriptid"]] == "ENST000001;ENST000002"
    assert own[columns["Ensembl_geneid"]] == "ENSG1;ENSG1"
    assert own[columns["aapos"]] == "4;4"
    assert own[columns["hg19_chr"]] == "21" and own[columns["hg19_pos(1-based)"]] == str(100 + mpf.OTHER_ASSEMBLY_SHIFT)
    assert len(own[columns["SIFT_score"]].split(";")) == 2
    other_change = by_key[("21", "100", "A", "G", "K", "E")]
    assert other_change[columns["Ensembl_transcriptid"]] == "ENST000003"
    stop = by_key[("21", "200", "C", "T", "Q", "X")]
    assert stop[columns["Ensembl_transcriptid"]] == "ENST000005"
    fabricated = [r for r in rows if r[0] == "21" and r[1] == "100" and r[3] != "G"]
    assert sorted(r[3] for r in fabricated) == ["C", "T"]
    assert all(r[columns["Ensembl_transcriptid"]] == "ENST000001;ENST000002;ENST000003;ENST000004" for r in fabricated)
    assert all(r[4] != r[5] for r in fabricated)
    mito = [r for r in rows if r[0] == "M"]
    assert len(mito) == 3 and {r[1] for r in mito} == {"3308"}
    assert by_key[("M", "3308", "T", "C", "M", "T")][columns["Ensembl_transcriptid"]] == "ENST000006"
    planted = [r for r in rows if r[0] == "1"]
    assert len(planted) == 1
    region = mpf.DETECTION_REGION
    assert region[1] <= int(planted[0][1]) <= region[2]
    for column in mpf.TRANSCRIPT_FILLERS:
        assert ";" in planted[0][columns[column]], column
    assert [r[0] for r in rows] == sorted((r[0] for r in rows), key=str)
    for chrom in ("1", "21", "M"):
        positions = [int(r[1]) for r in rows if r[0] == chrom]
        assert positions == sorted(positions)
    assert Path(f"{out}.tbi").is_file()
    spec = json.loads(Path(f"{out}.spec.json").read_text(encoding="utf-8"))
    assert spec["rows"] == {"own_allele": 5, "fabricated": 11, "planted": 1, "miss": 0}
    assert spec["preparation"]["index"] == "tabix -s 1 -b 2 -e 2"
    assert spec["inputs"]["vcf"] == {"variants.vcf": mpf.sha256(tmp_path / "variants.vcf")}


@needs_tools
def test_dbnsfp_grch37_file_carries_the_corpus_in_the_hg19_columns(tmp_path: Path):
    """The GRCh37 preparation puts the corpus chromosome and position in `hg19_chr` and
    `hg19_pos(1-based)`, indexes those columns and sorts the rows by them, with the hg38
    position shifted."""
    out = run(tmp_path, "dbnsfp", "GRCh37", "dbNSFP4.9a_grch37_synthetic.gz")
    assert tbi_columns(out) == (8, 9, 9, "#")
    rows = rows_of(out)
    own = [r for r in rows if r[7] == "21" and r[8] == "100" and r[3] == "G"]
    assert len(own) == 2
    assert {r[0] for r in own} == {"21"} and {r[1] for r in own} == {str(100 + mpf.OTHER_ASSEMBLY_SHIFT)}
    planted = [r for r in rows if r[7] == "1"]
    assert len(planted) == 1 and planted[0][8] == str(mpf.PLANTED_ROW["pos"]) and planted[0][1] == planted[0][8]


@needs_tools
def test_dbscsnv_rows_per_assembly(tmp_path: Path):
    """The GRCh38 file carries the corpus coordinates in `hg38_chr`/`hg38_pos`, the GRCh37
    file in `chr`/`pos`; both carry three rows per site, the bare header, CRLF line ends and
    the site's Ensembl gene ids, and no row at a `--miss` site."""
    out38 = run(tmp_path, "dbscsnv", "GRCh38", "dbscSNV1.1_synthetic_GRCh38.txt.gz", "--miss", "21:600")
    assert header_of(out38) == "\t".join(mpf.DBSCSNV_COLUMNS) + "\r\n"
    assert tbi_columns(out38) == (5, 6, 6, "c")
    rows = rows_of(out38)
    with gzip.open(out38, "rb") as handle:
        assert handle.read().count(b"\r\n") == len(rows) + 1
    at_100 = [r for r in rows if r[4] == "21" and r[5] == "100"]
    assert sorted(r[3] for r in at_100) == ["C", "G", "T"]
    assert all(r[0] == "21" and r[1] == str(100 + mpf.OTHER_ASSEMBLY_SHIFT) for r in at_100)
    assert all(r[13] == "ENSG1" for r in at_100)
    assert not any(r[5] == "600" for r in rows)
    assert len(rows) == 4 * 3
    out37 = run(tmp_path, "dbscsnv", "GRCh37", "dbscSNV1.1_synthetic_GRCh37.txt.gz")
    assert tbi_columns(out37) == (1, 2, 2, "c")
    rows37 = rows_of(out37)
    at_100 = [r for r in rows37 if r[0] == "21" and r[1] == "100"]
    assert len(at_100) == 3 and all(r[5] == str(100 + mpf.OTHER_ASSEMBLY_SHIFT) for r in at_100)
    assert len(rows37) == 5 * 3
    spec = json.loads(Path(f"{out37}.spec.json").read_text(encoding="utf-8"))
    assert spec["rows"] == {"snv": 15, "miss": 0} and spec["preparation"]["line_end"] == "CRLF"


@needs_tools
def test_seed_is_reproducible(tmp_path: Path):
    """The same inputs and seed write the same bytes; another seed writes other values at the
    same rows."""
    first = run(tmp_path, "dbnsfp", "GRCh38", "a_synthetic.gz")
    first_bytes = gzip.decompress(first.read_bytes())
    second = run(tmp_path, "dbnsfp", "GRCh38", "b_synthetic.gz")
    assert gzip.decompress(second.read_bytes()) == first_bytes
    vcf, reference = write_inputs(tmp_path)
    other = tmp_path / "plugin_data" / "c_synthetic.gz"
    mpf.main(["--layout", "dbnsfp", "--assembly", "GRCh38", "--seed", "39", "--vcf", str(vcf), "--reference", str(reference), "--out", str(other), "--force"])
    other_rows = rows_of(other)
    first_rows = rows_of(first)
    assert [r[:4] for r in other_rows] == [r[:4] for r in first_rows]
    assert other_rows != first_rows


@needs_tools
def test_provenance_records_the_fixture_relative_to_the_corpus(tmp_path: Path):
    """With --provenance the corpus record gains a `plugin_fixtures` entry naming the file,
    its index and its spec with their digests, paths relative to the corpus, every other key
    kept."""
    vcf, reference = write_inputs(tmp_path)
    provenance = tmp_path / "provenance.json"
    provenance.write_text(json.dumps({"vep_image": "x", "plugin_fixtures": {"other.gz": {"seed": 1}}}), encoding="utf-8")
    out = tmp_path / "plugin_data" / "dbscSNV1.1_synthetic_GRCh38.txt.gz"
    mpf.main(["--layout", "dbscsnv", "--assembly", "GRCh38", "--seed", "7", "--vcf", str(vcf), "--reference", str(reference), "--out", str(out), "--provenance", str(provenance)])
    record = json.loads(provenance.read_text(encoding="utf-8"))
    assert record["vep_image"] == "x" and record["plugin_fixtures"]["other.gz"] == {"seed": 1}
    entry = record["plugin_fixtures"]["dbscSNV1.1_synthetic_GRCh38.txt.gz"]
    assert entry["generator"] == mpf.GENERATOR and entry["seed"] == 7 and entry["layout"] == "dbscsnv"
    assert set(entry["files"]) == {
        "plugin_data/dbscSNV1.1_synthetic_GRCh38.txt.gz",
        "plugin_data/dbscSNV1.1_synthetic_GRCh38.txt.gz.tbi",
        "plugin_data/dbscSNV1.1_synthetic_GRCh38.txt.gz.spec.json",
    }
    assert entry["files"]["plugin_data/dbscSNV1.1_synthetic_GRCh38.txt.gz"] == mpf.sha256(out)
    assert "inputs" not in entry
    spec = json.loads(Path(f"{out}.spec.json").read_text(encoding="utf-8"))
    assert spec["inputs"]["reference"] == {"default.txt.gz": mpf.sha256(reference)}
    with pytest.raises(SystemExit, match="exists; pass --force"):
        mpf.main(["--layout", "dbscsnv", "--assembly", "GRCh38", "--seed", "7", "--vcf", str(vcf), "--reference", str(reference), "--out", str(out)])


def test_help_gives_the_data_terms_as_the_plugins_state_them(capsys):
    """`--help` words each file's terms from the plugin's module (dbNSFP.pm:54-55 names the two
    branches of its files; dbscSNV.pm names no terms below its licence header, lines 1-16)."""
    with pytest.raises(SystemExit) as stop:
        mpf.main(["--help"])
    assert stop.value.code == 0
    text = " ".join(capsys.readouterr().out.split())
    assert (
        "Ensembl VEP's dbNSFP plugin documents two branches of its files, one for academic and one "
        "for commercial use (dbNSFP.pm:54-55), and its dbscSNV plugin states no terms for its file, "
        "its licence header covering the module alone (dbscSNV.pm:1-16). Neither data file is "
        "redistributed here, so their plugins are verified on made-up fixtures"
    ) in text
    for claim in ("free for academic use", "paid agreement", "publishes no licence text"):
        assert claim not in text


def test_refusals(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """A --layout that names no layout, a malformed --miss, a missing input and an output
    outside the provenance's corpus are refused before anything is written."""
    vcf, reference = write_inputs(tmp_path)
    base = ["--assembly", "GRCh38", "--seed", "1", "--vcf", str(vcf), "--reference", str(reference), "--out", str(tmp_path / "out.gz")]
    with pytest.raises(SystemExit) as refused:
        mpf.main(["--layout", "spliceai_snv", *base])
    assert refused.value.code == 2
    assert "argument --layout: invalid choice: 'spliceai_snv'" in capsys.readouterr().err
    with pytest.raises(SystemExit, match="--miss names a site as chr:pos"):
        mpf.main(["--layout", "dbnsfp", *base, "--miss", "21-100"])
    with pytest.raises(SystemExit, match="input file not found"):
        mpf.main(["--layout", "dbnsfp", "--assembly", "GRCh38", "--seed", "1", "--vcf", str(tmp_path / "none.vcf"), "--reference", str(reference), "--out", str(tmp_path / "out.gz")])
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    with pytest.raises(SystemExit, match="lies outside the corpus"):
        mpf.main(["--layout", "dbnsfp", *base, "--provenance", str(elsewhere / "provenance.json")])
    assert not (tmp_path / "out.gz").exists()
