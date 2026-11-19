"""Tests for the synthetic plugin score file generator on hand-built inputs."""
from __future__ import annotations

import gzip
import hashlib
import json
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plugin_score_fixtures as psf  # noqa: E402

VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "21\t100\tsnv\tA\tG\t.\t.\t.",
    "21\t200\tmulti\tC\tT,A\t.\t.\t.",
    "21\t300\tdel\tAT\tA\t.\t.\t.",
    "21\t400\tins\tA\tAT\t.\t.\t.",
    "21\t500\tmnv\tAG\tCT\t.\t.\t.",
    "21\t600\tsv\tG\t<DEL>\t.\t.\tEND=700",
    "21\t700\tstar\tG\t*\t.\t.\t.",
    "21\t800\tlower\tt\tc\t.\t.\t.",
    "21\t900\tmissed\tG\tA\t.\t.\t.",
    "22\t50\tother_chr\tT\tG\t.\t.\t.",
]) + "\n"


def write_corpus(tmp_path: Path) -> Path:
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    (corpus / "variants.vcf").write_text(VCF, encoding="utf-8")
    return corpus


def test_alleles_are_read_as_the_plugins_compare_them(tmp_path: Path):
    """A shared leading base is dropped and the position moves past it for a deletion but
    stays on the anchor for an insertion; an empty allele is `-`; a longer substitution is
    kept whole; a symbolic or `*` allele writes nothing; alleles are upper-cased."""
    corpus = write_corpus(tmp_path)
    alleles = psf.read_alleles(corpus / "variants.vcf")
    assert [(a.chrom, a.pos, a.ref, a.alt) for a in alleles] == [
        ("21", 100, "A", "G"),
        ("21", 200, "C", "T"),
        ("21", 200, "C", "A"),
        ("21", 301, "T", "-"),
        ("21", 400, "-", "T"),
        ("21", 500, "AG", "CT"),
        ("21", 800, "T", "C"),
        ("21", 900, "G", "A"),
        ("22", 50, "T", "G"),
    ]
    assert [a.is_snv for a in alleles] == [True, True, True, False, False, False, True, True, True]


def test_alleles_follow_both_of_ensembl_veps_trimming_stages():
    """A record with one alternate allele of another length loses its shared leading base
    and then every shared leading and trailing base, so a homopolymer insertion is `-` and
    the inserted bases past the repeat, queried from the repeat's last base, and an interior
    deletion is the deleted base alone; alleles of one length are never trimmed; a record
    with several alternate alleles loses the leading base only when every one of them
    shares it, and nothing more; bases are compared as written and upper-cased afterwards;
    an alternate allele equal to the reference, and a record carrying a symbolic allele,
    write nothing."""
    records = [
        ("TG", ["TGG"], [("-", "G", 101)]),
        ("AA", ["AC"], [("AA", "AC", 100)]),
        ("GAC", ["GC"], [("A", "-", 101)]),
        ("C", ["T", "CA"], [("C", "T", 100), ("C", "CA", 100)]),
        ("A", ["AT", "AG"], [("-", "T", 100), ("-", "G", 100)]),
        ("TTGA", ["t"], [("TTGA", "T", 100)]),
        ("TTGA", ["Ttgatc"], [("TGA", "TGATC", 101)]),
        ("TTGA", ["TTGATC"], [("-", "TC", 103)]),
        ("TTGA", ["T", "TTGA"], [("TGA", "-", 101)]),
        ("AT", ["ATAT"], [("-", "AT", 101)]),
        ("CT", ["C", "*"], [("T", "-", 101)]),
        ("C", ["T", "<DEL>"], []),
        ("C", ["c"], []),
    ]
    for ref, alts, expected in records:
        got = [(a.ref, a.alt, a.pos) for a in psf.trimmed("21", 100, ref, alts)]
        assert got == expected, f"{ref} > {','.join(alts)}"


def test_snv_sites_get_the_three_other_bases_and_misses_get_none(tmp_path: Path):
    """Every single-nucleotide site gets one row per alternate base other than its
    reference, in position order with the chromosomes in natural order; a `--miss` site gets
    none; the clinpred and primateai shapes write no row for an insertion, deletion or
    longer substitution."""
    corpus = write_corpus(tmp_path)
    alleles = psf.read_alleles(corpus / "variants.vcf")
    rows, counts = psf.rows_for(alleles, "clinpred", 7, "bare", {"21:900"})
    cells = [row.split("\t")[:4] for row in rows]
    assert cells == [
        ["21", "100", "A", "C"], ["21", "100", "A", "G"], ["21", "100", "A", "T"],
        ["21", "200", "C", "A"], ["21", "200", "C", "G"], ["21", "200", "C", "T"],
        ["21", "800", "T", "A"], ["21", "800", "T", "C"], ["21", "800", "T", "G"],
        ["22", "50", "T", "A"], ["22", "50", "T", "C"], ["22", "50", "T", "G"],
    ]
    assert counts == {"snv": 12, "non_snv": 0, "miss": 1}
    assert all(len(row.split("\t")) == 5 for row in rows)


def test_bayesdel_adds_the_other_alleles_in_the_compared_representation(tmp_path: Path):
    """The bayesdel shape adds the deletion at its moved position with `-` for the alternate
    allele, the insertion at its anchor with `-` for the reference, and the two-base
    substitution whole, each once, sorted among the single-base rows."""
    corpus = write_corpus(tmp_path)
    alleles = psf.read_alleles(corpus / "variants.vcf")
    rows, counts = psf.rows_for(alleles, "bayesdel", 7, "bare", set())
    cells = [row.split("\t") for row in rows]
    non_snv = [row[:4] for row in cells if row[2] == "-" or row[3] == "-" or len(row[2]) > 1]
    assert non_snv == [["21", "301", "T", "-"], ["21", "400", "-", "T"], ["21", "500", "AG", "CT"]]
    assert counts == {"snv": 15, "non_snv": 3, "miss": 0}
    positions = [(psf.natural_key(row[0]), int(row[1])) for row in cells]
    assert positions == sorted(positions)
    scores = [float(row[4]) for row in cells]
    assert all(psf.BAYESDEL_RANGE[0] <= s <= psf.BAYESDEL_RANGE[1] for s in scores)


def test_scores_depend_on_the_seed_and_the_row_alone(tmp_path: Path):
    """The same seed gives the same rows whatever else the file holds; another seed moves
    every score; the two decoy rows of a site carry scores other than the picked row's."""
    corpus = write_corpus(tmp_path)
    alleles = psf.read_alleles(corpus / "variants.vcf")
    rows_all, _ = psf.rows_for(alleles, "primateai", 3, "chr", set())
    rows_fewer, _ = psf.rows_for(alleles, "primateai", 3, "chr", {"22:50"})
    assert rows_fewer == [row for row in rows_all if not row.startswith("chr22\t")]
    rows_other_seed, _ = psf.rows_for(alleles, "primateai", 4, "chr", set())
    assert [r.split("\t")[-1] for r in rows_all] != [r.split("\t")[-1] for r in rows_other_seed]
    site_100 = [r.split("\t")[-1] for r in rows_all if r.startswith("chr21\t100\t")]
    assert len(site_100) == 3 and len(set(site_100)) == 3
    assert all(r.startswith("chr") for r in rows_all)


def test_primateai_rows_carry_the_eleven_header_columns(tmp_path: Path):
    corpus = write_corpus(tmp_path)
    alleles = psf.read_alleles(corpus / "variants.vcf")
    rows, _ = psf.rows_for(alleles, "primateai", 1, "chr", set())
    header = psf.PRIMATEAI_HEADER.lstrip("#").split("\t")
    assert header == ["chr", "pos", "ref", "alt", "refAA", "altAA", "strand_1pos_0neg", "trinucleotide_context", "UCSC_gene", "ExAC_coverage", "primateDL_score"]
    for row in rows:
        cells = row.split("\t")
        assert len(cells) == len(header)
        assert cells[4] in psf.AMINO_ACIDS and cells[5] in psf.AMINO_ACIDS and cells[4] != cells[5]
        assert cells[6] in ("0", "1")
        assert len(cells[7]) == 3 and cells[7][1] == cells[2]
        assert cells[8].startswith("NM_")
        assert 0.0 <= float(cells[10]) <= 1.0
    assert len(psf.PRIMATEAI_PREAMBLE) == 11 and all(line.startswith("#") for line in psf.PRIMATEAI_PREAMBLE)


def test_main_writes_the_indexed_file_and_records_it(tmp_path: Path, capsys):
    """The tool writes the file as BGZF (the gzip magic with the extra-field flag) with its
    tabix index beside it under data/, the header on the line the shape names, records the
    seed, the misses, both digests and the tool versions under provenance.json's
    plugin_fixtures key beside the record's other keys and another fixture's entry, writes
    the same bytes twice from one seed, and refuses to overwrite without --force."""
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus = write_corpus(tmp_path)
    (corpus / "provenance.json").write_text(json.dumps({"vep_image": "x", "plugin_fixtures": {"other_synthetic.tsv.gz": {"seed": 1}}}), encoding="utf-8")
    name = "PrimateAI_scores_v0.2_GRCh38_sorted_synthetic.tsv.bgz"
    args = ["--corpus", str(corpus), "--shape", "primateai", "--name", name, "--seed", "38", "--chr-style", "chr", "--miss", "21:900"]
    assert psf.main(args) == 0
    assert capsys.readouterr().out.strip() == f"wrote data/{name} and its .tbi: 12 rows (snv 12, non_snv 0, missed sites 1), seed 38"
    out = corpus / "data" / name
    assert out.read_bytes()[:4] == b"\x1f\x8b\x08\x04"
    assert Path(f"{out}.tbi").read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[11] == psf.PRIMATEAI_HEADER and all(line.startswith("#") for line in lines[:12])
    assert lines[12].startswith("chr21\t100\tA\tC\t")
    assert sorted(p.name for p in (corpus / "data").iterdir()) == [name, f"{name}.tbi"], "no scratch file survives"
    first = out.read_bytes()
    record = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    assert record["vep_image"] == "x" and record["plugin_fixtures"]["other_synthetic.tsv.gz"] == {"seed": 1}
    entry = record["plugin_fixtures"][name]
    assert entry["generator"] == "scripts/fixtures/plugin_score_fixtures.py"
    assert entry["shape"] == "primateai" and entry["seed"] == 38 and entry["chr_style"] == "chr"
    assert entry["misses"] == ["21:900"]
    assert entry["variants_vcf_sha256"] == hashlib.sha256(VCF.encode()).hexdigest()
    assert entry["rows"] == {"snv": 12, "non_snv": 0, "miss": 1}
    assert entry["files"] == {f"data/{name}": hashlib.sha256(first).hexdigest(), f"data/{name}.tbi": hashlib.sha256(Path(f"{out}.tbi").read_bytes()).hexdigest()}
    assert sorted(entry["tools"]) == ["bgzip", "tabix"]
    assert all(version.split()[-1][0].isdigit() for version in entry["tools"].values()), entry["tools"]
    with pytest.raises(SystemExit, match="pass --force"):
        psf.main(args)
    assert psf.main([*args, "--force"]) == 0
    assert out.read_bytes() == first
    queried = subprocess.run(["tabix", str(out), "chr21:200-200"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[3] for line in queried] == ["A", "G", "T"]


MITOCHONDRIAL_RECORDS = [
    "MT\t40\tmt_snv\tC\tT\t.\t.\t.",
    "MT\t60\tmt_del\tGA\tG\t.\t.\t.",
    "chrM\t50\tmt_other_spelling\tT\tC\t.\t.\t.",
    "MT\t70\tmt_missed\tA\tG\t.\t.\t.",
]
MITOCHONDRIAL_SPELLINGS = ("MT\t", "chrM\t")


def write_mitochondrial_corpus(tmp_path: Path) -> Path:
    corpus = write_corpus(tmp_path)
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as fh:
        fh.write("".join(f"{line}\n" for line in MITOCHONDRIAL_RECORDS))
    return corpus


def test_the_mitochondrial_contig_is_written_under_the_given_name(tmp_path: Path):
    """Under a mitochondrial name every spelling of that contig is written as the name and
    its rows merge in position order behind the numbered contigs; every other row, and
    every cell of a renamed row after its first, are the bytes the input's spelling gives,
    so no score moves; the contig's own name changes nothing; one position under two
    spellings is refused."""
    corpus = write_mitochondrial_corpus(tmp_path)
    alleles = psf.read_alleles(corpus / "variants.vcf")
    rows, counts = psf.rows_for(alleles, "bayesdel", 7, "bare", {"MT:70"})
    assert counts == {"snv": 21, "non_snv": 4, "miss": 1}
    named = psf.with_mitochondrial_name(rows, "M")
    spelt = [row for row in rows if row.startswith(MITOCHONDRIAL_SPELLINGS)]
    renamed = [row for row in named if row.startswith("M\t")]
    assert len(spelt) == 7 and len(renamed) == 7
    assert [row for row in named if not row.startswith("M\t")] == [row for row in rows if row not in spelt]
    assert sorted(row.split("\t", 1)[1] for row in renamed) == sorted(row.split("\t", 1)[1] for row in spelt)
    assert [int(row.split("\t")[1]) for row in renamed] == [40, 40, 40, 50, 50, 50, 61]
    assert named[-7:] == renamed, "the contig's rows stand together, behind the numbered contigs"
    own = [row for row in rows if row.startswith("MT\t")]
    assert psf.with_mitochondrial_name(own, "MT") == own
    assert [row.split("\t")[0] for row in psf.with_mitochondrial_name(own, "chrMT")] == ["chrMT"] * 4
    with pytest.raises(SystemExit, match="position 50 of the mitochondrial contig has rows as MT and as chrM"):
        psf.with_mitochondrial_name(["MT\t50\tT\tC\t0.1", "chrM\t50\tT\tC\t0.2"], "M")


def test_header_names_stand_in_the_header_line_of_a_shape_read_by_place():
    """The names are written one per column behind the layout's `#`, for the bayesdel and
    clinpred shapes alone; another count than the layout's, an empty name, a name holding
    white space or opening with `#`, and a shape whose header is read by name or has
    several lines are refused."""
    names = psf.header_names_for("bayesdel", "Chr,Start,Ref,Alt,Score")
    assert names == ["Chr", "Start", "Ref", "Alt", "Score"]
    assert psf.headers_for("bayesdel", None, names) == ["#Chr\tStart\tRef\tAlt\tScore"]
    assert psf.headers_for("clinpred", None, psf.header_names_for("clinpred", "a,b,c,d,e")) == ["#a\tb\tc\td\te"]
    assert psf.headers_for("bayesdel", None) == [psf.BAYESDEL_HEADER]
    with pytest.raises(SystemExit, match="4 given, the bayesdel header has 5 columns"):
        psf.header_names_for("bayesdel", "Chr,Pos,Ref,Alt")
    for given in ("Chr,,Ref,Alt,Score", "Chr,Po s,Ref,Alt,Score", "#Chr,Pos,Ref,Alt,Score"):
        with pytest.raises(SystemExit, match="a name is not empty, holds no white space"):
            psf.header_names_for("bayesdel", given)
    for shape in ("primateai", "splicevault", "cadd_snv", "cadd_sv", "spliceai_snv"):
        with pytest.raises(SystemExit, match="--shape bayesdel or clinpred and with no other shape"):
            psf.header_names_for(shape, "a,b,c,d,e")


def test_main_writes_the_named_contig_and_header_and_records_both(tmp_path: Path, capsys):
    """The tool writes the mitochondrial rows under the given name, where tabix lists the
    contig and finds its rows, and the given names in the header line, and records both
    under the fixture's entry; without the two options the entry carries neither key, and
    the rows of the other contigs are the same bytes either way."""
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus = write_mitochondrial_corpus(tmp_path)
    name = "BayesDel_170824_addAF_all_scores_synthetic.txt.gz"
    out = corpus / "data" / name
    args = ["--corpus", str(corpus), "--shape", "bayesdel", "--name", name, "--seed", "5", "--miss", "MT:70"]
    assert psf.main(args) == 0
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        plain = fh.read().splitlines()
    assert plain[0] == psf.BAYESDEL_HEADER
    plain_entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert "mitochondrial_name" not in plain_entry and "header_names" not in plain_entry
    capsys.readouterr()
    header_names = "Chr,pOs,Ref,Alt,bAYESdEL_ADDaf"
    assert psf.main([*args, "--mitochondrial-name", "M", "--header-names", header_names, "--force"]) == 0
    assert capsys.readouterr().out.strip() == f"wrote data/{name} and its .tbi: 25 rows (snv 21, non_snv 4, missed sites 1), seed 5"
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "#Chr\tpOs\tRef\tAlt\tbAYESdEL_ADDaf"
    assert [line for line in lines[1:] if not line.startswith("M\t")] == [line for line in plain[1:] if not line.startswith(MITOCHONDRIAL_SPELLINGS)]
    assert sorted(line.split("\t", 1)[1] for line in lines if line.startswith("M\t")) == sorted(line.split("\t", 1)[1] for line in plain if line.startswith(MITOCHONDRIAL_SPELLINGS))
    contigs = subprocess.run(["tabix", "-l", str(out)], capture_output=True, text=True, check=True).stdout.split()
    assert contigs == ["21", "22", "M"]
    queried = subprocess.run(["tabix", str(out), "M:40-61"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[1] for line in queried] == ["40", "40", "40", "50", "50", "50", "61"]
    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert entry["mitochondrial_name"] == "M" and entry["header_names"] == header_names.split(",")
    assert entry["rows"] == plain_entry["rows"] and entry["seed"] == 5 and entry["misses"] == ["MT:70"]
    assert {key: value for key, value in entry.items() if key not in ("mitochondrial_name", "header_names", "files")} == {key: value for key, value in plain_entry.items() if key != "files"}


TRANSCRIPTS = "\n".join([
    "# site\ttranscript",
    "21:100\tENST00000300305",
    "21:100\tENST00000000001",
    "21:100\tENST00000300305",
    "21:300\tENST00000300305",
    "21:400\tENST00000346798",
    "21:800\tENST00000300305",
    "21:900\tENST00000300305",
    "",
]) + "\n"


def write_splicevault_corpus(tmp_path: Path) -> tuple[Path, Path]:
    corpus = write_corpus(tmp_path)
    (corpus / "data").mkdir()
    transcripts = corpus / "data" / "splicevault_transcripts.tsv"
    transcripts.write_text(TRANSCRIPTS, encoding="utf-8")
    return corpus, transcripts


def test_records_keep_the_alleles_as_the_vcf_writes_them(tmp_path: Path):
    """The splicevault shape reads the anchored alleles as written and upper-cased, one per
    alternate allele, and nothing for a symbolic, `*` or reference-equal allele."""
    corpus = write_corpus(tmp_path)
    records = psf.read_records(corpus / "variants.vcf")
    assert [(r.chrom, r.pos, r.ref, r.alt) for r in records] == [
        ("21", 100, "A", "G"),
        ("21", 200, "C", "T"),
        ("21", 200, "C", "A"),
        ("21", 300, "AT", "A"),
        ("21", 400, "A", "AT"),
        ("21", 500, "AG", "CT"),
        ("21", 800, "T", "C"),
        ("21", 900, "G", "A"),
        ("22", 50, "T", "G"),
    ]
    assert [r.is_snv for r in records] == [True, True, True, False, False, False, True, True, True]


def test_transcripts_are_listed_once_per_site_in_file_order(tmp_path: Path):
    corpus, transcripts = write_splicevault_corpus(tmp_path)
    sites = psf.read_transcripts(transcripts)
    assert sites["21:100"] == ["ENST00000300305", "ENST00000000001"]
    assert sites["21:400"] == ["ENST00000346798"]
    bad = corpus / "bad.tsv"
    bad.write_text("21:100 ENST1\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="chr:pos<TAB>transcript_id"):
        psf.read_transcripts(bad)
    bad.write_text("# nothing\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="no site listed"):
        psf.read_transcripts(bad)


def test_splicevault_rows_follow_the_shape_layout(tmp_path: Path):
    """A single-nucleotide site gets one row per alternate base per listed transcript, an
    insertion or deletion one row per transcript with its anchored alleles, an unlisted or
    missed site none; rows sort by position, then reference, alternate allele and transcript
    id; every row has twelve columns whose events are `Top<rank>:<type>;<impact>;<percent>;
    <frame>` groups and whose frameshift fraction counts them; a percent in the thousands
    carries a space; rows depend on the seed and the row's identity alone."""
    corpus, transcripts = write_splicevault_corpus(tmp_path)
    records = psf.read_records(corpus / "variants.vcf")
    sites = psf.read_transcripts(transcripts)
    rows, counts = psf.splicevault_rows(records, sites, 6, "chr", {"21:900"})
    cells = [row.split("\t") for row in rows]
    assert [row[:5] for row in cells] == [
        ["chr21", "100", "A", "C", "ENST00000000001"], ["chr21", "100", "A", "C", "ENST00000300305"],
        ["chr21", "100", "A", "G", "ENST00000000001"], ["chr21", "100", "A", "G", "ENST00000300305"],
        ["chr21", "100", "A", "T", "ENST00000000001"], ["chr21", "100", "A", "T", "ENST00000300305"],
        ["chr21", "300", "AT", "A", "ENST00000300305"],
        ["chr21", "400", "A", "AT", "ENST00000346798"],
        ["chr21", "800", "T", "A", "ENST00000300305"], ["chr21", "800", "T", "C", "ENST00000300305"],
        ["chr21", "800", "T", "G", "ENST00000300305"],
    ]
    assert counts == {"snv": 9, "non_snv": 2, "miss": 1}
    event = r"Top[1-4]:(ES;\d+(-\d+)?|C[DA];[+-]\d+);(\d+( \d{3})?%|0\.\d{3}%);(Frameshift|inFrame)"
    for row in cells:
        assert len(row) == 12
        assert row[5] in psf.SPLICEVAULT_SITE_TYPES
        assert row[6].startswith("chr21:")
        assert 0.2 <= float(row[7]) <= 0.99
        events = row[9].split("|")
        assert all(re.fullmatch(event, e) for e in events), row[9]
        assert [e.split(":")[0] for e in events] == [f"Top{i}" for i in range(1, len(events) + 1)]
        assert row[8] == f"Frameshift:{sum(e.endswith('Frameshift') for e in events)}/{len(events)}"
        assert row[10].isdigit() and row[11].isdigit()
    percents = [psf.splicevault_percent(random.Random(i)) for i in range(60)]
    assert any(" " in p for p in percents) and any(p.startswith("0.") for p in percents)
    again, _ = psf.splicevault_rows(records, sites, 6, "chr", set())
    assert [r for r in again if not r.startswith("chr21\t900\t")] == rows
    other_seed, _ = psf.splicevault_rows(records, sites, 7, "chr", {"21:900"})
    assert other_seed != rows
    assert cells[0][5:] != cells[1][5:], "the decoy row and the transcript's row differ"


def test_main_writes_the_splicevault_fixture_and_records_its_transcripts(tmp_path: Path, capsys):
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus, transcripts = write_splicevault_corpus(tmp_path)
    name = "SpliceVault_data_GRCh38_synthetic.tsv.gz"
    args = ["--corpus", str(corpus), "--shape", "splicevault", "--name", name, "--seed", "6", "--chr-style", "chr", "--transcripts", str(transcripts)]
    assert psf.main(args) == 0
    assert capsys.readouterr().out.strip() == f"wrote data/{name} and its .tbi: 14 rows (snv 12, non_snv 2, missed sites 0), seed 6"
    out = corpus / "data" / name
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == psf.SPLICEVAULT_HEADER
    assert lines[1].startswith("chr21\t100\tA\tC\tENST00000000001\t")
    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert entry["shape"] == "splicevault" and entry["rows"] == {"snv": 12, "non_snv": 2, "miss": 0}
    assert entry["transcripts"] == {"data/splicevault_transcripts.tsv": hashlib.sha256(TRANSCRIPTS.encode()).hexdigest()}
    queried = subprocess.run(["tabix", str(out), "chr21:400-400"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[:5] for line in queried] == [["chr21", "400", "A", "AT", "ENST00000346798"]]
    with pytest.raises(SystemExit, match="--transcripts is given with --shape splicevault"):
        psf.main(["--corpus", str(corpus), "--shape", "splicevault", "--name", name, "--seed", "6", "--force"])
    with pytest.raises(SystemExit, match="--transcripts is given with --shape splicevault"):
        psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "6", "--transcripts", str(transcripts)])
    with pytest.raises(SystemExit, match="file not found"):
        psf.main(["--corpus", str(corpus), "--shape", "splicevault", "--name", name, "--seed", "6", "--force", "--transcripts", str(corpus / "none.tsv")])


SYMBOLS = "\n".join([
    "# site\tsymbol",
    "21:100\tSOD1",
    "21:100\tOTHER",
    "21:100\tSOD1",
    "21:200\tGENE2",
    "21:300\tDELG",
    "21:400\tINSG",
    "21:500\tMNVG",
    "21:800\tLOW",
    "21:900\tMISSED",
    "",
]) + "\n"


def write_spliceai_corpus(tmp_path: Path) -> tuple[Path, Path]:
    corpus = write_corpus(tmp_path)
    (corpus / "plugin_data").mkdir()
    symbols = corpus / "plugin_data" / "spliceai_symbols.tsv"
    symbols.write_text(SYMBOLS, encoding="utf-8")
    return corpus, symbols


def test_spliceai_rows_follow_the_vcf_layout(tmp_path: Path):
    """The snv shape gives a single-nucleotide site one row per alternate base per listed
    symbol and skips every other allele; the indel shape gives an insertion, a deletion or
    a longer substitution one row per symbol with its alleles as the VCF writes them and
    skips the single-base sites; an unlisted or missed site gets none; every row has the
    eight VCF columns, `.` for ID, QUAL and FILTER, and an INFO of the allele, the symbol,
    four two-decimal delta scores and four integer delta positions; rows sort by position,
    then reference, alternate allele and symbol, and depend on the seed and the row alone."""
    corpus, symbols_path = write_spliceai_corpus(tmp_path)
    records = psf.read_records(corpus / "variants.vcf")
    symbols = psf.read_symbols(symbols_path)
    assert symbols["21:100"] == ["SOD1", "OTHER"]
    snv_rows, snv_counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", {"21:900"})
    snv_cells = [row.split("\t") for row in snv_rows]
    assert [row[:5] for row in snv_cells] == [
        ["21", "100", ".", "A", "C"], ["21", "100", ".", "A", "C"],
        ["21", "100", ".", "A", "G"], ["21", "100", ".", "A", "G"],
        ["21", "100", ".", "A", "T"], ["21", "100", ".", "A", "T"],
        ["21", "200", ".", "C", "A"], ["21", "200", ".", "C", "G"], ["21", "200", ".", "C", "T"],
        ["21", "800", ".", "T", "A"], ["21", "800", ".", "T", "C"], ["21", "800", ".", "T", "G"],
    ]
    assert snv_counts == {"snv": 12, "non_snv": 0, "miss": 1}
    assert [row[7].split("|")[1] for row in snv_cells[:2]] == ["OTHER", "SOD1"]
    info = r"SpliceAI=[ACGT]+\|[A-Z0-9]+(\|(0\.\d{2}|1\.00)){4}(\|-?\d{1,2}){4}"
    for row in snv_cells:
        assert len(row) == 8 and row[2] == "." and row[5] == "." and row[6] == "."
        assert re.fullmatch(info, row[7]), row[7]
        assert row[7].split("|")[0] == f"SpliceAI={row[4]}"
    indel_rows, indel_counts = psf.spliceai_rows(records, symbols, "spliceai_indel", 9, "chr", set())
    assert [row.split("\t")[:5] for row in indel_rows] == [
        ["chr21", "300", ".", "AT", "A"],
        ["chr21", "400", ".", "A", "AT"],
        ["chr21", "500", ".", "AG", "CT"],
    ]
    assert [row.split("\t")[7].split("|")[:2] for row in indel_rows] == [
        ["SpliceAI=A", "DELG"], ["SpliceAI=AT", "INSG"], ["SpliceAI=CT", "MNVG"],
    ]
    assert indel_counts == {"snv": 0, "non_snv": 3, "miss": 0}
    again, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set())
    assert [r for r in again if not r.startswith("21\t900\t")] == snv_rows
    other_seed, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 10, "bare", {"21:900"})
    assert other_seed != snv_rows
    assert snv_cells[0][7].split("|")[2:] != snv_cells[1][7].split("|")[2:], "the two symbols' rows differ"


def test_main_writes_the_spliceai_fixture_under_its_data_dir_and_records_its_symbols(tmp_path: Path, capsys):
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus, symbols = write_spliceai_corpus(tmp_path)
    name = "spliceai_scores.raw.indel.hg38.synthetic.vcf.gz"
    args = ["--corpus", str(corpus), "--shape", "spliceai_indel", "--name", name, "--seed", "9", "--data-dir", "plugin_data", "--symbols", str(symbols)]
    assert psf.main(args) == 0
    assert capsys.readouterr().out.strip() == f"wrote plugin_data/{name} and its .tbi: 3 rows (snv 0, non_snv 3, missed sites 0), seed 9"
    out = corpus / "plugin_data" / name
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[:3] == psf.SPLICEAI_HEADER
    assert "SpliceAIv1.3 synthetic fixture" in lines[1]
    assert lines[3].startswith("21\t300\t.\tAT\tA\t.\t.\tSpliceAI=A|DELG|")
    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert entry["shape"] == "spliceai_indel" and entry["rows"] == {"snv": 0, "non_snv": 3, "miss": 0}
    assert entry["symbols"] == {"plugin_data/spliceai_symbols.tsv": hashlib.sha256(SYMBOLS.encode()).hexdigest()}
    assert set(entry["files"]) == {f"plugin_data/{name}", f"plugin_data/{name}.tbi"}
    # Under the VCF preset the deletion's record spans its reference allele, so a query
    # for the deleted base alone finds it, as Ensembl VEP's plugin queries a deletion.
    queried = subprocess.run(["tabix", str(out), "21:301-301"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[:5] for line in queried] == [["21", "300", ".", "AT", "A"]]
    with pytest.raises(SystemExit, match="--symbols is given with --shape spliceai_snv or spliceai_indel"):
        psf.main(["--corpus", str(corpus), "--shape", "spliceai_snv", "--name", name, "--seed", "9", "--force"])
    with pytest.raises(SystemExit, match="--symbols is given with --shape spliceai_snv or spliceai_indel"):
        psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "9", "--symbols", str(symbols)])
    with pytest.raises(SystemExit, match="--symbols file not found"):
        psf.main(["--corpus", str(corpus), "--shape", "spliceai_snv", "--name", name, "--seed", "9", "--force", "--symbols", str(corpus / "none.tsv")])
    for escape in ("../x", ".."):
        with pytest.raises(SystemExit, match="--dir is a directory name inside the corpus"):
            psf.main(["--corpus", str(corpus), "--shape", "spliceai_snv", "--name", name, "--seed", "9", "--data-dir", escape, "--symbols", str(symbols)])


def test_spliceai_sites_take_the_other_forms_of_the_layout(tmp_path: Path):
    """A site under `--join-entries` has each allele's records joined into records of
    several entries, every entry in one or N a record, in the records' order; one under
    `--no-information` one record per allele whose INFO is `.`; one under
    `--join-alternates` one record listing its alternate bases with every entry; one under
    `--any-base-reference` its records with `N` as the reference. Every entry is the value
    of the record it comes from, a site under no option keeps its records, the class
    counts the rows written and each form its records. A site without a record, one whose
    alleles have a single record, and a call outside the options' shapes are refused."""
    corpus, _ = write_spliceai_corpus(tmp_path)
    records = psf.read_records(corpus / "variants.vcf")
    symbols = {"21:100": ["SOD1", "OTHER", "THIRD"], "21:200": ["GENE2"], "21:800": ["LOW"], "22:50": ["FAR"], "21:300": ["DELG", "DELH"]}
    plain, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set())
    value = {tuple(cells[:5]) + (cells[7].split("|")[1],): cells[7].removeprefix("SpliceAI=") for cells in (row.split("\t") for row in plain)}

    forms = psf.record_forms_for("spliceai_snv", ["21:100"], ["21:200"], ["21:800"], ["22:50"])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), forms)
    cells = [row.split("\t") for row in rows]
    assert [row[:7] for row in cells] == [
        ["21", "100", ".", "A", "C", ".", "."], ["21", "100", ".", "A", "G", ".", "."], ["21", "100", ".", "A", "T", ".", "."],
        ["21", "200", ".", "C", "A", ".", "."], ["21", "200", ".", "C", "G", ".", "."], ["21", "200", ".", "C", "T", ".", "."],
        ["21", "800", ".", "T", "A,C,G", ".", "."],
        ["22", "50", ".", "N", "A", ".", "."], ["22", "50", ".", "N", "C", ".", "."], ["22", "50", ".", "N", "G", ".", "."],
    ]
    assert counts == {"snv": 10, "non_snv": 0, "miss": 0, "join_entries": 3, "no_information": 3, "join_alternates": 1, "any_base_reference": 3}
    for row in cells[:3]:
        entries = row[7].removeprefix("SpliceAI=").split(",")
        assert [entry.split("|")[:2] for entry in entries] == [[row[4], "OTHER"], [row[4], "SOD1"], [row[4], "THIRD"]]
        assert entries == [value[("21", "100", ".", "A", row[4], symbol)] for symbol in ("OTHER", "SOD1", "THIRD")]
    assert [row[7] for row in cells[3:6]] == [".", ".", "."]
    assert cells[6][7] == "SpliceAI=" + ",".join(value[("21", "800", ".", "T", alt, "LOW")] for alt in "ACG")
    assert [row[7] for row in cells[7:]] == ["SpliceAI=" + value[("22", "50", ".", "T", alt, "FAR")] for alt in "ACG"]

    by_two = psf.record_forms_for("spliceai_snv", ["21:100=2"], [], [], [])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), by_two)
    at_site = [row.split("\t") for row in rows if row.startswith("21\t100\t")]
    assert [(row[4], row[7].count(",") + 1) for row in at_site] == [("C", 2), ("C", 1), ("G", 2), ("G", 1), ("T", 2), ("T", 1)]
    assert counts == {"snv": 15, "non_snv": 0, "miss": 0, "join_entries": 6}
    assert [row for row in rows if not row.startswith("21\t100\t")] == [row for row in plain if not row.startswith("21\t100\t")]

    indel = psf.record_forms_for("spliceai_indel", ["21:300"], [], [], [])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_indel", 9, "bare", set(), indel)
    assert [row.split("\t")[7].count(",") + 1 for row in rows] == [2] and counts == {"snv": 0, "non_snv": 1, "miss": 0, "join_entries": 1}

    assert psf.record_forms_for("spliceai_snv", [], [], [], []) is None
    with pytest.raises(SystemExit, match="--join-entries names 21:900, a site with no record in this file"):
        psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", ["21:900"], [], [], []))
    with pytest.raises(SystemExit, match="--join-entries names 21:200, whose alleles have one record each"):
        psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", ["21:200"], [], [], []))
    with pytest.raises(SystemExit, match="are given with --shape spliceai_snv or spliceai_indel and with no other shape"):
        psf.record_forms_for("clinpred", [], ["21:100"], [], [])
    with pytest.raises(SystemExit, match="--join-alternates is given with --shape spliceai_snv"):
        psf.record_forms_for("spliceai_indel", [], [], ["21:300"], [])
    with pytest.raises(SystemExit, match="--any-base-reference is given with --shape spliceai_snv"):
        psf.record_forms_for("spliceai_indel", [], [], [], ["21:300"])
    with pytest.raises(SystemExit, match=r"--join-entries names a site as chr:pos\[=N\], got '21:100=two'"):
        psf.record_forms_for("spliceai_snv", ["21:100=two"], [], [], [])
    with pytest.raises(SystemExit, match="--join-entries 21:100=1: a record of several entries holds two or more"):
        psf.record_forms_for("spliceai_snv", ["21:100=1"], [], [], [])
    with pytest.raises(SystemExit, match="21:100 stands under --join-entries and --no-information; a site takes one form"):
        psf.record_forms_for("spliceai_snv", ["21:100"], ["21:100"], [], [])


def test_spliceai_lines_stand_as_many_times_as_their_site_is_told(tmp_path: Path):
    """Under `--repeat-lines` every record line of a site stands N times, the copies next
    to each other: the plain lines of a site under no other option, and the lines of the
    form a site stands under beside it. Every other record of the file stands, the class
    counts the lines the file holds, a form the lines it holds in that form and
    `repeat_lines` the lines of its sites. A malformed word, a count below two, a site named
    twice, a site without a record and a call outside the spliceai shapes are refused."""
    corpus, _ = write_spliceai_corpus(tmp_path)
    records = psf.read_records(corpus / "variants.vcf")
    symbols = {"21:100": ["SOD1", "OTHER", "THIRD"], "21:200": ["GENE2"], "21:800": ["LOW"], "21:300": ["DELG", "DELH"]}
    plain, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set())
    joined, joined_counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", ["21:100"], [], [], []))
    assert joined_counts == {"snv": 9, "non_snv": 0, "miss": 0, "join_entries": 3}

    forms = psf.record_forms_for("spliceai_snv", ["21:100"], [], [], [], ["21:100=3", "21:200=2"])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), forms)
    assert rows[:9] == [line for line in joined[:3] for _ in range(3)]
    assert all(row.split("\t")[7].count(",") == 2 for row in rows[:9])
    assert rows[9:15] == [line for line in plain if line.startswith("21\t200\t") for _ in range(2)]
    assert rows[15:] == [line for line in plain if line.startswith("21\t800\t")]
    assert counts == {"snv": 18, "non_snv": 0, "miss": 0, "join_entries": 9, "repeat_lines": 15}

    alone = psf.record_forms_for("spliceai_snv", [], [], [], [], ["21:100=2"])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), alone)
    assert rows == [line for line in plain if line.startswith("21\t100\t") for _ in range(2)] + [line for line in plain if not line.startswith("21\t100\t")]
    assert counts == {"snv": 24, "non_snv": 0, "miss": 0, "repeat_lines": 18}
    again, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", [], [], [], [], ["21:100=2"]))
    assert again == rows

    indel = psf.record_forms_for("spliceai_indel", ["21:300"], [], [], [], ["21:300=2"])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_indel", 9, "bare", set(), indel)
    assert len(rows) == 2 and rows[0] == rows[1] and rows[0].split("\t")[7].count(",") == 1
    assert counts == {"snv": 0, "non_snv": 2, "miss": 0, "join_entries": 2, "repeat_lines": 2}

    for word in ("21:100", "21:100=", "21:100=two", "21=2", "21:100=2=3"):
        with pytest.raises(SystemExit, match=r"--repeat-lines names a site and the times its lines stand as chr:pos=N, got "):
            psf.record_forms_for("spliceai_snv", [], [], [], [], [word])
    with pytest.raises(SystemExit, match="--repeat-lines 21:100=1: a line that stands more than once stands twice or more"):
        psf.record_forms_for("spliceai_snv", [], [], [], [], ["21:100=1"])
    with pytest.raises(SystemExit, match="--repeat-lines names 21:100 twice; a site's lines stand one number of times"):
        psf.record_forms_for("spliceai_snv", [], [], [], [], ["21:100=2", "21:100=3"])
    with pytest.raises(SystemExit, match="are given with --shape spliceai_snv or spliceai_indel and with no other shape"):
        psf.record_forms_for("clinpred", [], [], [], [], ["21:100=2"])
    with pytest.raises(SystemExit, match="--repeat-lines names 21:900, a site with no record in this file"):
        psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", [], [], [], [], ["21:900=2"]))


def test_main_writes_the_repeated_lines_and_records_their_sites(tmp_path: Path, capsys):
    """A file written under `--repeat-lines` holds each line of the named sites that many
    times, tabix reads every copy, and the words are recorded with the counts of the lines
    the file holds."""
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus, symbols = write_spliceai_corpus(tmp_path)
    name = "repeated.synthetic.vcf.gz"
    base = ["--corpus", str(corpus), "--shape", "spliceai_snv", "--seed", "9", "--data-dir", "plugin_data", "--symbols", str(symbols), "--name", name]
    assert psf.main([*base, "--join-entries", "21:100", "--repeat-lines", "21:900=3", "--repeat-lines", "21:100=2"]) == 0
    assert capsys.readouterr().out.strip() == (
        f"wrote plugin_data/{name} and its .tbi: 21 rows (snv 21, non_snv 0, missed sites 0, join_entries 6, repeat_lines 15), seed 9"
    )
    out = corpus / "plugin_data" / name
    queried = subprocess.run(["tabix", str(out), "21:100-100"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [(line.split("\t")[4], line.split("\t")[7].count(",") + 1) for line in queried] == [("C", 2), ("C", 2), ("G", 2), ("G", 2), ("T", 2), ("T", 2)]
    assert queried[0] == queried[1] and queried[2] == queried[3] and queried[4] == queried[5]
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        missed = [line for line in fh if line.startswith("21\t900\t")]
    assert len(missed) == 9 and len(set(missed)) == 3
    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert entry["record_forms"] == {"join_entries": ["21:100"], "repeat_lines": ["21:100=2", "21:900=3"]}
    assert entry["rows"] == {"snv": 21, "non_snv": 0, "miss": 0, "join_entries": 6, "repeat_lines": 15}


def test_spliceai_lines_stand_again_with_scores_drawn_again(tmp_path: Path, capsys):
    """Under `--rescore-lines` every record line of a site stands N times, the copies next to
    each other: the first as its form writes it, each further copy naming the same alleles
    and symbols with scores of the layout that depend on the seed, the entry's record and
    the copy's place alone. Every other record of the file stands, a form counts the lines
    the file holds in it and `rescore_lines` the lines of its sites, and the words are
    recorded. A malformed word, a count below two, a site named twice, a site under
    `--repeat-lines` too, a site without a record and a call outside the spliceai shapes
    are refused."""
    corpus, symbols_path = write_spliceai_corpus(tmp_path)
    records = psf.read_records(corpus / "variants.vcf")
    symbols = {"21:100": ["SOD1", "OTHER"], "21:200": ["GENE2"], "21:800": ["LOW"]}
    plain, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set())
    joined, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", ["21:100"], [], [], []))

    forms = psf.record_forms_for("spliceai_snv", ["21:100"], [], [], [], ["21:800=2"], ["21:100=2", "21:200=3"])
    rows, counts = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), forms)
    assert counts == {"snv": 21, "non_snv": 0, "miss": 0, "join_entries": 6, "repeat_lines": 6, "rescore_lines": 15}
    assert rows[0:6:2] == joined[:3]
    entry = r"[ACGT]\|[A-Z0-9]+(\|(0\.\d{2}|1\.00)){4}(\|-?\d{1,2}){4}"
    for first, second in zip(rows[0:6:2], rows[1:6:2]):
        assert first.split("\t")[:7] == second.split("\t")[:7]
        ours, theirs = (line.split("\t")[7].removeprefix("SpliceAI=").split(",") for line in (first, second))
        assert [item.split("|")[:2] for item in ours] == [item.split("|")[:2] for item in theirs]
        assert all(re.fullmatch(entry, item) for item in theirs) and all(a != b for a, b in zip(ours, theirs))
    at_200 = [line for line in plain if line.startswith("21\t200\t")]
    assert rows[6:15:3] == at_200 and len(set(rows[6:15])) == 9
    assert [line.split("\t")[:5] for line in rows[6:15]] == [line.split("\t")[:5] for line in at_200 for _ in range(3)]
    assert rows[15:] == [line for line in plain if line.startswith("21\t800\t") for _ in range(2)]

    again, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", ["21:100"], [], [], [], ["21:800=2"], ["21:100=2", "21:200=3"]))
    assert again == rows
    by_two, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", [], [], [], [], [], ["21:200=2"]))
    assert by_two[6:12] == [line for index, line in enumerate(rows[6:15]) if index % 3 != 2], "a copy's scores do not move with the count"
    other_seed, _ = psf.spliceai_rows(records, symbols, "spliceai_snv", 10, "bare", set(), psf.record_forms_for("spliceai_snv", [], [], [], [], [], ["21:200=2"]))
    assert other_seed[7] != by_two[7]

    for word in ("21:100", "21:100=", "21:100=two", "21=2"):
        with pytest.raises(SystemExit, match=r"--rescore-lines names a site and the times its lines stand as chr:pos=N, got "):
            psf.record_forms_for("spliceai_snv", [], [], [], [], [], [word])
    with pytest.raises(SystemExit, match="--rescore-lines 21:100=1: a line that stands more than once stands twice or more"):
        psf.record_forms_for("spliceai_snv", [], [], [], [], [], ["21:100=1"])
    with pytest.raises(SystemExit, match="--rescore-lines names 21:100 twice; a site's lines stand one number of times"):
        psf.record_forms_for("spliceai_snv", [], [], [], [], [], ["21:100=2", "21:100=3"])
    with pytest.raises(SystemExit, match="21:100 stands under --repeat-lines and --rescore-lines; the copies of a site's lines are written one way"):
        psf.record_forms_for("spliceai_snv", [], [], [], [], ["21:100=2"], ["21:100=2"])
    with pytest.raises(SystemExit, match="are given with --shape spliceai_snv or spliceai_indel and with no other shape"):
        psf.record_forms_for("clinpred", [], [], [], [], [], ["21:100=2"])
    with pytest.raises(SystemExit, match="--rescore-lines names 21:900, a site with no record in this file"):
        psf.spliceai_rows(records, symbols, "spliceai_snv", 9, "bare", set(), psf.record_forms_for("spliceai_snv", [], [], [], [], [], ["21:900=2"]))

    if any(shutil.which(tool) is None for tool in psf.TOOLS):
        return
    name = "rescored.synthetic.vcf.gz"
    base = ["--corpus", str(corpus), "--shape", "spliceai_snv", "--seed", "9", "--data-dir", "plugin_data", "--symbols", str(symbols_path), "--name", name]
    assert psf.main([*base, "--join-entries", "21:100", "--rescore-lines", "21:100=2"]) == 0
    assert capsys.readouterr().out.strip() == (
        f"wrote plugin_data/{name} and its .tbi: 15 rows (snv 15, non_snv 0, missed sites 0, join_entries 6, rescore_lines 6), seed 9"
    )
    queried = subprocess.run(["tabix", str(corpus / "plugin_data" / name), "21:100-100"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert queried == rows[:6]
    recorded = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert recorded["record_forms"] == {"join_entries": ["21:100"], "rescore_lines": ["21:100=2"]}
    assert recorded["rows"] == {"snv": 15, "non_snv": 0, "miss": 0, "join_entries": 6, "rescore_lines": 6}


def bgzf_members(data: bytes) -> list[tuple[int, int]]:
    """Each bgzf block of `data` as its compressed offset and the size it inflates to, read
    from the block size of its header's `BC` subfield and from its `ISIZE` trailer."""
    members, at = [], 0
    while at < len(data):
        assert data[at : at + 4] == b"\x1f\x8b\x08\x04" and data[at + 12 : at + 14] == b"BC"
        size = int.from_bytes(data[at + 16 : at + 18], "little") + 1
        members.append((at, int.from_bytes(data[at + size - 4 : at + size], "little")))
        at += size
    return members


def test_main_writes_two_joined_streams_and_the_record_forms_and_records_them(tmp_path: Path, capsys):
    """Under `--second-stream-at` the file is two bgzip streams joined: the first stream's
    end marker, a block that holds no bytes, stands before the block the named row opens,
    the lines are those of the file written as one stream, and tabix reads the rows on
    either side of the marker. The option and the sites of the record forms are recorded."""
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus, symbols = write_spliceai_corpus(tmp_path)
    one, two = "one_stream.synthetic.vcf.gz", "two_streams.synthetic.vcf.gz"
    base = ["--corpus", str(corpus), "--shape", "spliceai_snv", "--seed", "9", "--data-dir", "plugin_data", "--symbols", str(symbols)]
    assert psf.main([*base, "--name", one]) == 0
    assert psf.main([*base, "--name", two, "--second-stream-at", "21:200"]) == 0
    capsys.readouterr()
    data = (corpus / "plugin_data" / two).read_bytes()
    members = bgzf_members(data)
    empty = [index for index, (_, inflated) in enumerate(members) if inflated == 0]
    assert len(empty) == 2 and empty[1] == len(members) - 1 and 0 < empty[0] < empty[1] - 1, members
    assert [inflated == 0 for _, inflated in bgzf_members((corpus / "plugin_data" / one).read_bytes())].count(True) == 1
    with gzip.open(corpus / "plugin_data" / one, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    with gzip.open(corpus / "plugin_data" / two, "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines() == lines
    second = gzip.decompress(data[members[empty[0] + 1][0] : members[empty[1]][0]]).decode().splitlines()
    assert second == [line for line in lines if line.startswith(("21\t200\t", "21\t800\t", "21\t900\t"))]
    for region, rows in (("21:100-100", 6), ("21:200-200", 3), ("21:800-800", 3)):
        queried = subprocess.run(["tabix", str(corpus / "plugin_data" / two), region], capture_output=True, text=True, check=True).stdout.splitlines()
        assert len(queried) == rows, region
    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert fixtures[two]["second_stream_at"] == "21:200" and "second_stream_at" not in fixtures[one]
    assert "record_forms" not in fixtures[two]

    forms = "forms.synthetic.vcf.gz"
    assert psf.main([*base, "--name", forms, "--join-entries", "21:100", "--no-information", "21:200", "--any-base-reference", "21:800"]) == 0
    assert capsys.readouterr().out.strip() == (
        f"wrote plugin_data/{forms} and its .tbi: 12 rows (snv 12, non_snv 0, missed sites 0, join_entries 3, no_information 3, any_base_reference 3), seed 9"
    )
    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][forms]
    assert entry["record_forms"] == {"join_entries": ["21:100"], "no_information": ["21:200"], "any_base_reference": ["21:800"]}
    assert entry["rows"] == {"snv": 12, "non_snv": 0, "miss": 0, "join_entries": 3, "no_information": 3, "any_base_reference": 3}
    with pytest.raises(SystemExit, match="--second-stream-at names 21:150, where the file has no row"):
        psf.main([*base, "--name", "x.synthetic.vcf.gz", "--second-stream-at", "21:150"])
    with pytest.raises(SystemExit, match="--second-stream-at names a site as chr:pos, got '21'"):
        psf.main([*base, "--name", "x.synthetic.vcf.gz", "--second-stream-at", "21"])


def test_main_refuses_a_malformed_call(tmp_path: Path):
    corpus = write_corpus(tmp_path)
    with pytest.raises(SystemExit, match="file name, not a path"):
        psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "data/x_synthetic.gz", "--seed", "1"])
    with pytest.raises(SystemExit, match="chr:pos"):
        psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "1", "--miss", "21"])
    with pytest.raises(SystemExit, match="no variants.vcf"):
        psf.main(["--corpus", str(tmp_path), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "1"])
    with pytest.raises(SystemExit, match="--header-names is given with --shape bayesdel or clinpred"):
        psf.main(["--corpus", str(corpus), "--shape", "primateai", "--name", "x_synthetic.gz", "--seed", "1", "--header-names", "a,b"])
    with pytest.raises(SystemExit, match="1 given, the clinpred header has 5 columns"):
        psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "1", "--header-names", ""])
    (corpus / "variants.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n21\t600\tsv\tG\t<DEL>\n", encoding="utf-8")
    if all(shutil.which(tool) for tool in psf.TOOLS):
        with pytest.raises(SystemExit, match="no site to score"):
            psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "1"])


CADD_VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "21\t100\tsnv\tA\tG\t.\t.\t.",
    "21\t200\tmulti\tC\tT,A\t.\t.\t.",
    "21\t300\tdel\tAT\tA\t.\t.\t.",
    "21\t400\tins\tA\tAT\t.\t.\t.",
    "21\t500\tmnv\tAG\tCT\t.\t.\t.",
    "21\t600\tshared\tACG\tATG\t.\t.\t.",
    "21\t900\tmissed\tG\tA\t.\t.\t.",
    "21\t5000\tdel_sv\tT\t<DEL>\t.\t.\tSVTYPE=DEL;END=6000",
    "21\t7000\tins_sv\tT\t<INS>\t.\t.\tSVTYPE=INS",
    "21\t8000\tdup_sv\tT\t<DUP>\t.\t.\tSVTYPE=DUP;SVLEN=500",
    "21\t9000\tinv_sv\tT\t<INV>\t.\t.\tSVTYPE=INV;END=9500",
    "21\t9600\ttandem_sv\tT\t<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;END=9900",
    "21\t9950\tbare_del_sv\tT\t<DEL>\t.\t.\tSVTYPE=DEL",
    "21\t12000\tmissed_sv\tT\t<DEL>\t.\t.\tSVTYPE=DEL;END=12500",
]) + "\n"


def write_cadd_corpus(tmp_path: Path) -> Path:
    corpus = tmp_path / "cadd"
    corpus.mkdir()
    (corpus / "variants.vcf").write_text(CADD_VCF, encoding="utf-8")
    return corpus


def test_cadd_snv_rows_score_the_single_base_change_after_trimming(tmp_path: Path):
    """A site is an allele's one changed base once the common bases at both ends are
    trimmed: the substitution `ACG/ATG` at 600 is the `C/T` site at 601, scored beside the
    two other bases; a deletion, an insertion and a two-base substitution are no site; a
    missed site gets no row; the rows carry the six header columns."""
    corpus = write_cadd_corpus(tmp_path)
    rows, counts = psf.cadd_snv_rows(psf.read_alleles(corpus / "variants.vcf"), 17, "bare", {"21:900"})
    cells = [row.split("\t") for row in rows]
    assert [row[:4] for row in cells] == [
        ["21", "100", "A", "C"], ["21", "100", "A", "G"], ["21", "100", "A", "T"],
        ["21", "200", "C", "A"], ["21", "200", "C", "G"], ["21", "200", "C", "T"],
        ["21", "601", "C", "A"], ["21", "601", "C", "G"], ["21", "601", "C", "T"],
    ]
    assert counts == {"snv": 9, "non_snv": 0, "miss": 1}
    for row in cells:
        assert len(row) == 6
        assert psf.CADD_RAW_RANGE[0] <= float(row[4]) <= psf.CADD_RAW_RANGE[1]
        assert psf.CADD_PHRED_RANGE[0] <= float(row[5]) <= psf.CADD_PHRED_RANGE[1]
    assert psf.trimmed_both_ends("ACG", "ATG") == ("C", "T", 1)
    assert psf.trimmed_both_ends("TTGA", "T") == ("TGA", "-", 1)
    assert psf.trimmed_both_ends("-", "T") == ("-", "T", 0)


def test_cadd_indel_rows_keep_the_anchor_and_add_a_longer_decoy(tmp_path: Path):
    """An insertion and a deletion are written as the VCF writes them, anchor base kept,
    each beside a decoy one base longer at the same anchor; substitutions and single-base
    changes get no row; a missed site (its VCF position) gets none."""
    corpus = write_cadd_corpus(tmp_path)
    rows, counts = psf.cadd_indel_rows(psf.read_records(corpus / "variants.vcf"), 17, "bare", set())
    assert [row.split("\t")[:4] for row in rows] == [
        ["21", "300", "AT", "A"], ["21", "300", "ATT", "A"],
        ["21", "400", "A", "AT"], ["21", "400", "A", "ATT"],
    ]
    assert counts == {"snv": 0, "non_snv": 4, "miss": 0}
    rows, counts = psf.cadd_indel_rows(psf.read_records(corpus / "variants.vcf"), 17, "chr", {"21:300"})
    assert [row.split("\t")[:4] for row in rows] == [["chr21", "400", "A", "AT"], ["chr21", "400", "A", "ATT"]]
    assert counts == {"snv": 0, "non_snv": 2, "miss": 1}


def test_cadd_sv_rows_follow_the_shape_layout(tmp_path: Path):
    """A `<DEL>`, `<INS>` or `<DUP>` record gets its own row, with `Start` the VCF position
    and `End` from `END`, from `SVLEN`, or the next base for an insertion, beside a decoy
    of another type over the same span and one of its type one base along; an inversion, a
    tandem duplication, a deletion with neither `END` nor `SVLEN` and a missed record get
    none; rows sort by start and carry the seven header columns."""
    corpus = write_cadd_corpus(tmp_path)
    records = psf.read_structural(corpus / "variants.vcf")
    assert [(r.chrom, r.start, r.end, r.kind) for r in records] == [
        ("21", 5000, 6000, "DEL"), ("21", 7000, 7001, "INS"), ("21", 8000, 8500, "DUP"), ("21", 12000, 12500, "DEL"),
    ]
    rows, counts = psf.cadd_sv_rows(records, 17, "bare", {"21:12000"})
    cells = [row.split("\t") for row in rows]
    assert [row[:4] for row in cells] == [
        ["21", "5000", "6000", "DEL"], ["21", "5000", "6000", "DUP"], ["21", "5001", "6000", "DEL"],
        ["21", "7000", "7001", "DEL"], ["21", "7000", "7001", "INS"], ["21", "7001", "7001", "INS"],
        ["21", "8000", "8500", "DEL"], ["21", "8000", "8500", "DUP"], ["21", "8001", "8500", "DUP"],
    ]
    assert counts == {"sv": 3, "decoy": 6, "miss": 1}
    for row in cells:
        assert len(row) == 7
        assert row[4] == f"fixture_{row[3]}_{row[0]}:{row[1]}-{row[2]}"
        assert psf.CADD_SV_PHRED_RANGE[0] <= float(row[5]) <= psf.CADD_SV_PHRED_RANGE[1]
        assert psf.CADD_SV_RAW_RANGE[0] <= float(row[6]) <= psf.CADD_SV_RAW_RANGE[1]
    assert cells[0][5:] != cells[1][5:], "the decoy of another type carries other scores"
    again, _ = psf.cadd_sv_rows(records, 17, "bare", {"21:12000"})
    assert again == rows
    other_seed, _ = psf.cadd_sv_rows(records, 18, "bare", {"21:12000"})
    assert other_seed != rows


def test_cadd_headers_name_the_build():
    assert psf.headers_for("cadd_snv", "GRCh37") == [
        "##CADD GRCh37-v1.7 layout: a synthetic score file for a vep-rs golden corpus; "
        "the scores are fabricated from a seed and carry no meaning",
        "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED",
    ]
    assert psf.headers_for("cadd_indel", "GRCh38")[0].startswith("##CADD GRCh38-v1.7 ")
    for shape in ("cadd_snv", "cadd_indel"):
        for assembly in psf.CADD_ASSEMBLIES:
            first = psf.headers_for(shape, assembly)[0]
            assert first.startswith(f"##CADD {assembly}-v1.7 layout: a synthetic score file")
            assert "(c)" not in first and "rights reserved" not in first.lower()
    assert psf.headers_for("cadd_sv", "GRCh38") == [
        "## CADD-SV GRCh38 layout: a synthetic score file for a vep-rs golden corpus; the scores are fabricated from a seed and carry no meaning",
        "#Chrom\tStart\tEnd\tType\tName\tCADD-SV_PHRED-score\tCADD-SV_Raw-score",
    ]
    assert psf.headers_for("clinpred", None) == [psf.CLINPRED_HEADER]


def test_main_writes_the_cadd_fixtures_and_their_indexes(tmp_path: Path, capsys):
    """The SNV file's first line names the build and says the file is synthetic, and the
    file is indexed on its position; the CADD-SV file is indexed on a zero-based start and
    an exclusive end, so a tabix query of the base after a deletion's
    start returns its rows; provenance records the assembly; the assembly option is refused
    where it has no meaning."""
    missing = [tool for tool in psf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus = write_cadd_corpus(tmp_path)
    snv = ["--corpus", str(corpus), "--shape", "cadd_snv", "--name", "whole_genome_SNVs_synthetic.tsv.gz", "--seed", "17", "--assembly", "GRCh37", "--miss", "21:900"]
    assert psf.main(snv) == 0
    assert capsys.readouterr().out.strip() == "wrote data/whole_genome_SNVs_synthetic.tsv.gz and its .tbi: 9 rows (snv 9, non_snv 0, missed sites 1), seed 17"
    with gzip.open(corpus / "data" / "whole_genome_SNVs_synthetic.tsv.gz", "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[0].startswith("##CADD GRCh37-v1.7 layout: a synthetic score file") and lines[1] == psf.CADD_HEADER
    queried = subprocess.run(["tabix", str(corpus / "data" / "whole_genome_SNVs_synthetic.tsv.gz"), "21:601-601"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[:4] for line in queried] == [["21", "601", "C", "A"], ["21", "601", "C", "G"], ["21", "601", "C", "T"]]

    sv = ["--corpus", str(corpus), "--shape", "cadd_sv", "--name", "1000G_phase3_SVs_synthetic.tsv.gz", "--seed", "17", "--assembly", "GRCh38"]
    assert psf.main(sv) == 0
    assert capsys.readouterr().out.strip() == "wrote data/1000G_phase3_SVs_synthetic.tsv.gz and its .tbi: 12 rows (sv 4, decoy 8, missed sites 0), seed 17"
    out = corpus / "data" / "1000G_phase3_SVs_synthetic.tsv.gz"
    with gzip.open(out, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[:2] == [psf.CADD_SV_MARKER, psf.CADD_SV_HEADER]
    queried = subprocess.run(["tabix", str(out), "21:5001-5001"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[:4] for line in queried] == [["21", "5000", "6000", "DEL"], ["21", "5000", "6000", "DUP"]]
    queried = subprocess.run(["tabix", str(out), "21:5002-5002"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[:4] for line in queried] == [["21", "5000", "6000", "DEL"], ["21", "5000", "6000", "DUP"], ["21", "5001", "6000", "DEL"]]
    queried = subprocess.run(["tabix", str(out), "21:5000-5000"], capture_output=True, text=True, check=True).stdout.splitlines()
    assert queried == [], "a zero-based start of 5000 covers base 5001, not 5000"
    record = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert record["whole_genome_SNVs_synthetic.tsv.gz"]["assembly"] == "GRCh37"
    assert record["1000G_phase3_SVs_synthetic.tsv.gz"]["shape"] == "cadd_sv" and record["1000G_phase3_SVs_synthetic.tsv.gz"]["rows"] == {"sv": 4, "decoy": 8, "miss": 0}

    with pytest.raises(SystemExit, match="--assembly is given with a cadd shape"):
        psf.main(["--corpus", str(corpus), "--shape", "cadd_indel", "--name", "indels_synthetic.tsv.gz", "--seed", "1"])
    with pytest.raises(SystemExit, match="--assembly is given with a cadd shape"):
        psf.main(["--corpus", str(corpus), "--shape", "clinpred", "--name", "x_synthetic.gz", "--seed", "1", "--assembly", "GRCh38"])
    with pytest.raises(SystemExit, match="GRCh38 alone"):
        psf.main(["--corpus", str(corpus), "--shape", "cadd_sv", "--name", "sv_synthetic.tsv.gz", "--seed", "1", "--assembly", "GRCh37", "--force"])
