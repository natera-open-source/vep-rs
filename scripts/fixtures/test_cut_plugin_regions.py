"""Tests for cut_plugin_regions.py on a hand-built contig, cache and bigWig."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cut_plugin_regions as cut

TOOLS = ("bgzip", "bigtools")

# A 300-base contig: one transcript at 101-200 with exons 101-130 and 171-200, one
# record at 115 and one at 250, outside the transcript.
CONTIG = ("ACGT" * 75)[:300]


def gerp_value(base: int) -> float:
    """The GERP value the test bigWig carries at 0-based `base`."""
    return (base % 7) - 2.5


def write_fasta(path: Path, name: str, sequence: str, line_bases: int = 60) -> None:
    lines = [sequence[i : i + line_bases] for i in range(0, len(sequence), line_bases)]
    path.write_text(f">{name}\n" + "".join(f"{l}\n" for l in lines), encoding="utf-8")
    Path(f"{path}.fai").write_text(f"{name}\t{len(sequence)}\t{len(name) + 2}\t{line_bases}\t{line_bases + 1}\n", encoding="utf-8")


def build_corpus(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    corpus = tmp_path / "corpus"
    (corpus / "json_cache" / "transcripts" / "21").mkdir(parents=True)
    (corpus / "variants.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t115\tinside\tA\tG\t.\t.\t.\n21\t250\toutside\tC\tT\t.\t.\t.\n",
        encoding="utf-8",
    )
    transcript = {
        "stable_id": "ENST_TEST",
        "start": "101",
        "end": "200",
        "strand": 1,
        "exons": [{"start": 101, "end": 130}, {"start": 171, "end": 200}],
    }
    with gzip.open(corpus / "json_cache" / "transcripts" / "21" / "1-1000.json.gz", "wt", encoding="utf-8") as handle:
        json.dump([transcript], handle)
    reference = tmp_path / "genome.fa"
    write_fasta(reference, "21", CONTIG)
    ancestor = tmp_path / "ancestor.fa"
    write_fasta(ancestor, "21", CONTIG.lower().replace("a", "t"))
    bedgraph = tmp_path / "gerp.bedgraph"
    bedgraph.write_text("".join(f"21\t{i}\t{i + 1}\t{gerp_value(i)}\n" for i in range(300)), encoding="utf-8")
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("21\t300\n", encoding="utf-8")
    gerp = tmp_path / "gerp.bw"
    subprocess.run(["bigtools", "bedgraphtobigwig", str(bedgraph), str(sizes), str(gerp)], check=True)
    return corpus, reference, ancestor, gerp


def decompressed(path: Path) -> str:
    return gzip.open(path, "rt", encoding="utf-8").read()


def sequence_of(fasta_text: str) -> str:
    return "".join(l for l in fasta_text.splitlines() if not l.startswith(">"))


def bedgraph_of(bigwig: Path, tmp_path: Path) -> list[tuple[int, int, float]]:
    out = tmp_path / "out.bedgraph"
    subprocess.run(["bigtools", "bigwigtobedgraph", str(bigwig), str(out)], check=True)
    rows = []
    for line in out.read_text(encoding="utf-8").splitlines():
        chrom, start, end, value = line.split("\t")
        assert chrom == "21"
        rows.append((int(start), int(end), float(value)))
    return rows


@pytest.fixture
def tools():
    missing = [tool for tool in TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")


def test_cut_masks_outside_the_regions_and_records_the_fixtures(tmp_path: Path, tools):
    corpus, reference, ancestor, gerp = build_corpus(tmp_path)
    (corpus / "provenance.json").write_text(json.dumps({"input": "two records"}) + "\n", encoding="utf-8")
    rc = cut.main(["--corpus", str(corpus), "--reference", str(reference), "--ancestor", str(ancestor), "--gerp", str(gerp),
                   "--flank", "10", "--record-flank", "5", "--exon-flank", "2"])
    assert rc == 0

    # The reference keeps the transcript span plus the flank (91-210) and each record
    # plus its flank (110-120, 245-255); every other base is N, the contig its full length.
    sequence = sequence_of(decompressed(corpus / "reference.fa.gz"))
    assert len(sequence) == 300
    kept = set(range(91, 211)) | set(range(110, 121)) | set(range(245, 256))
    for position in range(1, 301):
        base = sequence[position - 1]
        assert base == (CONTIG[position - 1] if position in kept else "N"), position
    assert (corpus / "reference.fa.fai").read_text(encoding="utf-8") == "21\t300\t4\t60\t61\n"
    assert (corpus / "reference.fa.gz.fai").read_text(encoding="utf-8") == "21\t300\t4\t60\t61\n"
    assert (corpus / "reference.fa.gz.gzi").is_file()
    assert not (corpus / "reference.fa").exists()

    # The ancestral FASTA keeps the records' neighbourhoods alone.
    ancestral = sequence_of(decompressed(corpus / "plugin_data" / "human_ancestor.fa.gz"))
    assert len(ancestral) == 300
    record_kept = set(range(110, 121)) | set(range(245, 256))
    source = CONTIG.lower().replace("a", "t")
    for position in range(1, 301):
        assert ancestral[position - 1] == (source[position - 1] if position in record_kept else "N"), position
    assert (corpus / "plugin_data" / "human_ancestor.fa.gz.gzi").is_file()
    assert (corpus / "plugin_data" / "human_ancestor.fa.fai").is_file()

    # The bigWig carries the exons plus the flank (99-132, 169-202 in 1-based bases) and
    # no zoom level.
    rows = bedgraph_of(corpus / "plugin_data" / "gerp.bw", tmp_path)
    covered = set()
    for start, end, value in rows:
        for base in range(start, end):
            covered.add(base + 1)
            assert value == gerp_value(base)
    assert covered == set(range(99, 133)) | set(range(169, 203))
    info = subprocess.run(["bigtools", "bigwiginfo", str(corpus / "plugin_data" / "gerp.bw")], check=True, capture_output=True, text=True).stdout
    assert "zoomLevels: 0" in info

    record = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    assert record["input"] == "two records"
    fixtures = record["plugin_fixtures"]
    assert fixtures["generator"] == "scripts/fixtures/cut_plugin_regions.py"
    assert (fixtures["flank"], fixtures["record_flank"], fixtures["exon_flank"]) == (10, 5, 2)
    assert sorted(fixtures["tools"]) == ["bgzip", "bigtools"]
    assert set(fixtures["files"]) == {
        "reference.fa.gz", "reference.fa.gz.fai", "reference.fa.gz.gzi", "reference.fa.fai",
        "plugin_data/human_ancestor.fa.gz", "plugin_data/human_ancestor.fa.gz.fai", "plugin_data/human_ancestor.fa.gz.gzi",
        "plugin_data/human_ancestor.fa.fai", "plugin_data/gerp.bw",
    }
    for name, digest in fixtures["files"].items():
        assert hashlib.sha256((corpus / name).read_bytes()).hexdigest() == digest, name


def test_existing_fixtures_need_force(tmp_path: Path, tools):
    corpus, reference, ancestor, gerp = build_corpus(tmp_path)
    args = ["--corpus", str(corpus), "--reference", str(reference), "--ancestor", str(ancestor), "--gerp", str(gerp)]
    assert cut.main(args) == 0
    with pytest.raises(SystemExit, match="pass --force to overwrite the fixtures"):
        cut.main(args)
    assert cut.main([*args, "--force"]) == 0


def test_gerp_alone_is_cut_around_the_records(tmp_path: Path, tools):
    corpus, _reference, _ancestor, gerp = build_corpus(tmp_path)
    # A record spelt with the `chr` prefix reaches the cache's chromosome 21.
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("chr21\t180\tprefixed\tA\tC\t.\t.\t.\n")
    rc = cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--gerp-regions", "records", "--record-flank", "5"])
    assert rc == 0

    # The bigWig carries each record plus its flank (110-120, 175-185, 245-255 in 1-based
    # bases), the contig at the length the source bigWig declares, and nothing else is written.
    rows = bedgraph_of(corpus / "plugin_data" / "gerp.bw", tmp_path)
    covered = set()
    for start, end, value in rows:
        for base in range(start, end):
            covered.add(base + 1)
            assert value == gerp_value(base)
    assert covered == set(range(110, 121)) | set(range(175, 186)) | set(range(245, 256))
    info = subprocess.run(["bigtools", "bigwiginfo", "--chroms", str(corpus / "plugin_data" / "gerp.bw")], check=True, capture_output=True, text=True).stdout
    assert "zoomLevels: 0" in info
    assert "21 0 300" in info
    assert not (corpus / "reference.fa.gz").exists()
    assert not (corpus / "plugin_data" / "human_ancestor.fa.gz").exists()

    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert fixtures["gerp_regions"] == "records"
    assert set(fixtures["files"]) == {"plugin_data/gerp.bw"}


def test_gerp_name_sets_the_bigwig_file_name(tmp_path: Path, tools):
    corpus, _reference, _ancestor, gerp = build_corpus(tmp_path)
    name = "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
    rc = cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--gerp-name", name])
    assert rc == 0

    # The cut bigWig is written under the given name alone, with the same per-base
    # values, and the provenance record names it by that name.
    assert not (corpus / "plugin_data" / "gerp.bw").exists()
    rows = bedgraph_of(corpus / "plugin_data" / name, tmp_path)
    assert rows
    for start, end, value in rows:
        for base in range(start, end):
            assert value == gerp_value(base)
    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert set(fixtures["files"]) == {f"plugin_data/{name}"}

    for bad in ("sub/gerp.bw", "", "."):
        with pytest.raises(SystemExit, match="is not a file name"):
            cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--gerp-name", bad, "--force"])


def test_nothing_to_cut_is_refused(tmp_path: Path, tools):
    corpus, _reference, _ancestor, _gerp = build_corpus(tmp_path)
    with pytest.raises(SystemExit, match="nothing to cut"):
        cut.main(["--corpus", str(corpus)])


def test_a_record_chromosome_the_cache_lacks_is_refused(tmp_path: Path, tools):
    corpus, reference, ancestor, gerp = build_corpus(tmp_path)
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("22\t10\tother\tA\tC\t.\t.\t.\n")
    with pytest.raises(SystemExit, match="names chromosome 22, which the cache does not hold"):
        cut.main(["--corpus", str(corpus), "--reference", str(reference), "--ancestor", str(ancestor), "--gerp", str(gerp)])


def test_fetch_reads_across_line_breaks():
    fasta_text = ">21\n" + "".join(f"{CONTIG[i:i + 60]}\n" for i in range(0, 300, 60))
    path = Path(__file__).resolve().parent / "_fetch_probe.fa"
    try:
        path.write_text(fasta_text, encoding="utf-8")
        entry = (300, 4, 60, 61)
        assert cut.fetch(path, entry, 59, 62).decode() == CONTIG[58:62]
        assert cut.fetch(path, entry, 1, 300).decode() == CONTIG
        assert cut.fetch(path, entry, 290, 400).decode() == CONTIG[289:]
        assert cut.fetch(path, entry, 5, 4) == b""
    finally:
        path.unlink(missing_ok=True)


def test_fetch_reads_a_crlf_fasta_whole():
    # A FASTA with two-byte line ends indexes 62 bytes per 60-base line; the read
    # budget follows that width, so a span over many lines comes back complete.
    fasta_text = ">21\r\n" + "".join(f"{CONTIG[i:i + 60]}\r\n" for i in range(0, 300, 60))
    path = Path(__file__).resolve().parent / "_fetch_crlf_probe.fa"
    try:
        path.write_bytes(fasta_text.encode())
        entry = (300, 5, 60, 62)
        assert cut.fetch(path, entry, 1, 300).decode() == CONTIG
        assert cut.fetch(path, entry, 59, 62).decode() == CONTIG[58:62]
    finally:
        path.unlink(missing_ok=True)


def test_count_bases_stops_at_the_contig_ends():
    # A flank past either end of the contig is written as nothing, so it is not counted.
    assert cut.count_bases({"21": [(-10, 5), (295, 310)]}, {"21": 300}) == 11
    assert cut.count_bases({"21": [(10, 20), (15, 30)]}, {"21": 300}) == 21
    assert cut.count_bases({}, {"21": 300}) == 0


def test_merge_joins_touching_regions():
    assert cut.merge([(10, 20), (21, 30), (40, 50), (45, 60), (5, 8)]) == [(5, 8), (10, 30), (40, 60)]
