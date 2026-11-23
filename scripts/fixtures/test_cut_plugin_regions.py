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


def build_fasta_corpus(tmp_path: Path) -> tuple[Path, Path, Path]:
    """The corpus, the reference and the ancestral FASTA, written with no external tool."""
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
    return corpus, reference, ancestor


def build_corpus(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    corpus, reference, ancestor = build_fasta_corpus(tmp_path)
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


@pytest.fixture
def bgzip_alone(tmp_path: Path, monkeypatch):
    """A PATH whose one directory holds `bgzip` and no other program."""
    bgzip = shutil.which("bgzip")
    if bgzip is None:
        pytest.skip("bgzip not on PATH")
    programs = tmp_path / "programs"
    programs.mkdir()
    (programs / "bgzip").symlink_to(bgzip)
    monkeypatch.setenv("PATH", str(programs))


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


def test_a_cut_keeps_the_entries_other_generators_recorded(tmp_path: Path, bgzip_alone):
    corpus, reference, _ancestor = build_fasta_corpus(tmp_path)
    # A score-table generator records its file under `plugin_fixtures.<file name>`.
    earlier = {"generator": "scripts/fixtures/plugin_score_fixtures.py", "seed": 37, "files": {"plugin_data/scores.synthetic.txt.gz": "0" * 64}}
    later = {"generator": "scripts/fixtures/make_plugin_fixture.py", "seed": 7, "files": {"plugin_data/table.synthetic.gz": "1" * 64}}
    provenance = corpus / "provenance.json"
    provenance.write_text(json.dumps({"input": "two records", "plugin_fixtures": {"scores.synthetic.txt.gz": earlier}}) + "\n", encoding="utf-8")
    reported = subprocess.run(["bgzip", "--version"], check=True, capture_output=True, text=True).stdout.splitlines()[0]
    written = sorted(["reference.fa.gz", "reference.fa.gz.fai", "reference.fa.gz.gzi", "reference.fa.fai"])

    def cut_keys(flank: int) -> dict:
        return {
            "generator": "scripts/fixtures/cut_plugin_regions.py",
            "flank": flank,
            "record_flank": 5,
            "exon_flank": 50,
            "gerp_regions": "exons",
            "tools": {"bgzip": reported},
            "files": {name: hashlib.sha256((corpus / name).read_bytes()).hexdigest() for name in written},
        }

    # The first cut adds its keys after the entry the record held and drops nothing.
    args = ["--corpus", str(corpus), "--reference", str(reference), "--record-flank", "5"]
    assert cut.main([*args, "--flank", "10"]) == 0
    expected = {"input": "two records", "plugin_fixtures": {"scores.synthetic.txt.gz": earlier, **cut_keys(10)}}
    assert provenance.read_text(encoding="utf-8") == json.dumps(expected, indent=1) + "\n"

    # A forced cut rewrites its own keys where they stand, between an entry recorded
    # before the first cut and one recorded after it.
    record = json.loads(provenance.read_text(encoding="utf-8"))
    record["plugin_fixtures"]["table.synthetic.gz"] = later
    provenance.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    assert cut.main([*args, "--flank", "20", "--force"]) == 0
    expected = {"input": "two records", "plugin_fixtures": {"scores.synthetic.txt.gz": earlier, **cut_keys(20), "table.synthetic.gz": later}}
    assert provenance.read_text(encoding="utf-8") == json.dumps(expected, indent=1) + "\n"


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
    assert list(fixtures["tools"]) == ["bigtools"]


def test_an_ancestor_alone_is_cut_without_bigtools_and_counts_no_reference_base(tmp_path: Path, bgzip_alone, capsys):
    corpus, _reference, ancestor = build_fasta_corpus(tmp_path)
    assert shutil.which("bigtools") is None
    rc = cut.main(["--corpus", str(corpus), "--ancestor", str(ancestor), "--flank", "10", "--record-flank", "5"])
    assert rc == 0

    # The ancestral FASTA keeps each record plus its flank (110-120, 245-255) and no other
    # base of the transcript span, the contig at the length its own index gives.
    source = CONTIG.lower().replace("a", "t")
    record_kept = set(range(110, 121)) | set(range(245, 256))
    sequence = sequence_of(decompressed(corpus / "plugin_data" / "human_ancestor.fa.gz"))
    assert sequence == "".join(source[position - 1] if position in record_kept else "N" for position in range(1, 301))
    assert (corpus / "plugin_data" / "human_ancestor.fa.fai").read_text(encoding="utf-8") == "21\t300\t4\t60\t61\n"
    assert (corpus / "plugin_data" / "human_ancestor.fa.gz.fai").read_text(encoding="utf-8") == "21\t300\t4\t60\t61\n"

    # The corpus gains the four files of the ancestral FASTA and its record, and no reference.
    names = ["human_ancestor.fa.gz", "human_ancestor.fa.gz.fai", "human_ancestor.fa.gz.gzi", "human_ancestor.fa.fai"]
    written = [f"plugin_data/{name}" for name in names]
    assert {path.name for path in corpus.iterdir()} == {"json_cache", "variants.vcf", "provenance.json", "plugin_data"}
    assert {path.name for path in (corpus / "plugin_data").iterdir()} == set(names)
    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert set(fixtures["files"]) == set(written)
    reported = subprocess.run(["bgzip", "--version"], check=True, capture_output=True, text=True).stdout.splitlines()[0]
    assert fixtures["tools"] == {"bgzip": reported}

    # The summary counts no base of a reference the cut does not write.
    summary = {"transcripts": 1, "records": 2, "reference_bases_kept": 0, "gerp_bases": 0, "files": written}
    assert capsys.readouterr().out == json.dumps(summary) + "\n"


def test_each_written_file_keeps_the_contig_lengths_the_help_states(tmp_path: Path, tools, capsys):
    # The three sources give the one contig three lengths: the reference 290 bases, the
    # ancestral FASTA 280 and the bigWig 300.
    corpus, reference, ancestor, gerp = build_corpus(tmp_path)
    write_fasta(reference, "21", CONTIG[:290])
    write_fasta(ancestor, "21", CONTIG.lower().replace("a", "t")[:280])
    flanks = ["--flank", "10", "--record-flank", "5", "--exon-flank", "2"]

    def declared(bigwig: Path) -> str:
        return subprocess.run(["bigtools", "bigwiginfo", "--chroms", str(bigwig)], check=True, capture_output=True, text=True).stdout

    # With all three given, each masked FASTA keeps the length of its own index and the cut
    # bigWig declares the reference's.
    sources = ["--reference", str(reference), "--ancestor", str(ancestor), "--gerp", str(gerp)]
    assert cut.main(["--corpus", str(corpus), *sources, *flanks]) == 0
    assert (corpus / "reference.fa.fai").read_text(encoding="utf-8") == "21\t290\t4\t60\t61\n"
    assert (corpus / "plugin_data" / "human_ancestor.fa.fai").read_text(encoding="utf-8") == "21\t280\t4\t60\t61\n"
    assert len(sequence_of(decompressed(corpus / "reference.fa.gz"))) == 290
    assert len(sequence_of(decompressed(corpus / "plugin_data" / "human_ancestor.fa.gz"))) == 280
    assert "21 0 290" in declared(corpus / "plugin_data" / "gerp.bw")

    # Without the reference the cut bigWig declares the ancestral FASTA's length; cut alone
    # it declares its source's, which the test of that cut holds.
    assert cut.main(["--corpus", str(corpus), *sources[2:], *flanks, "--force"]) == 0
    assert "21 0 280" in declared(corpus / "plugin_data" / "gerp.bw")

    capsys.readouterr()
    with pytest.raises(SystemExit):
        cut.main(["--help"])
    printed = " ".join(capsys.readouterr().out.split())
    assert "Each masked FASTA keeps the contig lengths of its own index." in printed
    assert (
        "The cut bigWig declares the contig lengths of the reference's index, else of the ancestral "
        "FASTA's index, else the lengths the source bigWig declares."
    ) in printed
    assert "The contig lengths come from" not in printed


def test_a_reference_alone_is_cut_without_bigtools_and_writes_no_plugin_data(tmp_path: Path, bgzip_alone):
    corpus, reference, _ancestor = build_fasta_corpus(tmp_path)
    assert shutil.which("bigtools") is None
    rc = cut.main(["--corpus", str(corpus), "--reference", str(reference), "--flank", "10", "--record-flank", "5"])
    assert rc == 0

    # The reference keeps the transcript span plus the flank (91-210, which holds the record
    # at 115 and its flank) and the record at 250 plus its flank (245-255).
    sequence = sequence_of(decompressed(corpus / "reference.fa.gz"))
    kept = set(range(91, 211)) | set(range(245, 256))
    assert sequence == "".join(CONTIG[position - 1] if position in kept else "N" for position in range(1, 301))

    # The corpus gains the four files of the reference and its record, and no directory.
    written = {"reference.fa.gz", "reference.fa.gz.fai", "reference.fa.gz.gzi", "reference.fa.fai"}
    assert {path.name for path in corpus.iterdir()} == {"json_cache", "variants.vcf", "provenance.json", *written}
    assert not (corpus / "plugin_data").exists()
    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert set(fixtures["files"]) == written
    reported = subprocess.run(["bgzip", "--version"], check=True, capture_output=True, text=True).stdout.splitlines()[0]
    assert fixtures["tools"] == {"bgzip": reported}


def leave_alone_on_path(tool: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """PATH becomes one directory that holds a link to `tool` and nothing else."""
    only = tmp_path / "path"
    only.mkdir()
    (only / tool).symlink_to(Path(shutil.which(tool)).resolve())
    monkeypatch.setenv("PATH", str(only))


def test_a_bigwig_only_cut_asks_no_bgzip_version(tmp_path: Path, tools, monkeypatch):
    corpus, _reference, _ancestor, gerp = build_corpus(tmp_path)
    leave_alone_on_path("bigtools", tmp_path, monkeypatch)
    assert shutil.which("bgzip") is None
    rc = cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--gerp-regions", "records", "--record-flank", "5"])
    assert rc == 0

    # A cut that writes the bigWig alone runs bigtools alone, and the record names the
    # tools the cut ran, each with the version it reports.
    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert set(fixtures["files"]) == {"plugin_data/gerp.bw"}
    assert list(fixtures["tools"]) == ["bigtools"]
    assert fixtures["tools"]["bigtools"].split()[-1][0].isdigit(), fixtures["tools"]


def test_a_fasta_only_cut_asks_no_bigtools_version(tmp_path: Path, tools, monkeypatch):
    corpus, reference, ancestor, _gerp = build_corpus(tmp_path)
    leave_alone_on_path("bgzip", tmp_path, monkeypatch)
    assert shutil.which("bigtools") is None
    rc = cut.main(["--corpus", str(corpus), "--reference", str(reference), "--ancestor", str(ancestor)])
    assert rc == 0

    # A cut that writes the two FASTA files alone runs bgzip alone.
    fixtures = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"]
    assert not (corpus / "plugin_data" / "gerp.bw").exists()
    assert set(fixtures["files"]) == {
        "reference.fa.gz", "reference.fa.gz.fai", "reference.fa.gz.gzi", "reference.fa.fai",
        "plugin_data/human_ancestor.fa.gz", "plugin_data/human_ancestor.fa.gz.fai", "plugin_data/human_ancestor.fa.gz.gzi",
        "plugin_data/human_ancestor.fa.fai",
    }
    assert list(fixtures["tools"]) == ["bgzip"]
    assert fixtures["tools"]["bgzip"].split()[-1][0].isdigit(), fixtures["tools"]


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


def test_a_record_named_by_a_synonym_is_cut_on_the_cache_chromosome(tmp_path: Path, bgzip_alone):
    corpus, reference, _ancestor = build_fasta_corpus(tmp_path)
    # The lines link a name to the cache's 21 from the first column and from the second, link
    # `shared` to 22 before 21, and hold a blank line and a line of one name, which link nothing.
    (corpus / "json_cache" / "chr_synonyms.txt").write_text(
        "synonym_of_21\t21\n\nalone\n21\tlater_name\nshared\t22\nshared\t21\nelsewhere\t22\n", encoding="utf-8"
    )
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("synonym_of_21\t280\tby_synonym\tA\tC\t.\t.\t.\n")
        handle.write("later_name\t20\tby_later_name\tA\tC\t.\t.\t.\n")
        handle.write("shared\t50\tby_shared_name\tA\tC\t.\t.\t.\n")
    rc = cut.main(["--corpus", str(corpus), "--reference", str(reference), "--flank", "10", "--record-flank", "5"])
    assert rc == 0

    # The records at 280, 20 and 50 keep their flanks (275-285, 15-25, 45-55) on the contig
    # the cache names 21, the one contig of the cut; a synonym of a chromosome the cache lacks
    # resolves nothing.
    text = decompressed(corpus / "reference.fa.gz")
    assert [line for line in text.splitlines() if line.startswith(">")] == [">21"]
    kept = set(range(15, 26)) | set(range(45, 56)) | set(range(91, 211)) | set(range(245, 256)) | set(range(275, 286))
    assert sequence_of(text) == "".join(CONTIG[position - 1] if position in kept else "N" for position in range(1, 301))
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("elsewhere\t10\tother\tA\tC\t.\t.\t.\n")
    with pytest.raises(SystemExit, match="names chromosome elsewhere, which the cache does not hold"):
        cut.main(["--corpus", str(corpus), "--reference", str(reference), "--force"])


# A 120-base contig no transcript of the corpus's cache lies on.
SCAFFOLD = ("TGCA" * 30)[:120]


def add_contig(path: Path, name: str, sequence: str, line_bases: int = 60) -> None:
    """Appends a contig to a FASTA `write_fasta` wrote, and its line to the index."""
    offset = path.stat().st_size + len(name) + 2
    lines = [sequence[i : i + line_bases] for i in range(0, len(sequence), line_bases)]
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f">{name}\n" + "".join(f"{l}\n" for l in lines))
    with Path(f"{path}.fai").open("a", encoding="utf-8") as handle:
        handle.write(f"{name}\t{len(sequence)}\t{offset}\t{line_bases}\t{line_bases + 1}\n")


def test_a_record_on_a_sequence_without_a_transcript_keeps_its_region(tmp_path: Path, bgzip_alone, capsys):
    corpus, reference, _ancestor = build_fasta_corpus(tmp_path)
    add_contig(reference, "unplaced", SCAFFOLD)
    # The cache holds transcripts on 21 alone, and its synonyms file links `other_name` to
    # the reference's `unplaced`. One record names that contig as the reference does, one by
    # the synonym and one with the `chr` prefix.
    (corpus / "json_cache" / "chr_synonyms.txt").write_text("other_name\tunplaced\n", encoding="utf-8")
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("unplaced\t60\tas_written\tA\tC\t.\t.\t.\n")
        handle.write("other_name\t100\tby_synonym\tA\tC\t.\t.\t.\n")
        handle.write("chrunplaced\t20\tprefixed\tA\tC\t.\t.\t.\n")
    rc = cut.main(["--corpus", str(corpus), "--reference", str(reference), "--flank", "10", "--record-flank", "5"])
    assert rc == 0

    # The cut holds both contigs at their full length: 21 as a cut of the cache's chromosome
    # (91-210 and 245-255), `unplaced` with each of its records plus the flank (15-25, 55-65
    # and 95-105) and every other base N.
    text = decompressed(corpus / "reference.fa.gz")
    assert [line for line in text.splitlines() if line.startswith(">")] == [">21", ">unplaced"]
    first, second = text.split(">unplaced\n")
    kept = set(range(91, 211)) | set(range(245, 256))
    assert sequence_of(first) == "".join(CONTIG[position - 1] if position in kept else "N" for position in range(1, 301))
    kept = set(range(15, 26)) | set(range(55, 66)) | set(range(95, 106))
    assert sequence_of(second) == "".join(SCAFFOLD[position - 1] if position in kept else "N" for position in range(1, 121))
    index = "21\t300\t4\t60\t61\nunplaced\t120\t319\t60\t61\n"
    assert (corpus / "reference.fa.fai").read_text(encoding="utf-8") == index
    assert (corpus / "reference.fa.gz.fai").read_text(encoding="utf-8") == index

    # The run's one line counts the five records and the bases kept on both contigs (131 and 33).
    files = ["reference.fa.gz", "reference.fa.gz.fai", "reference.fa.gz.gzi", "reference.fa.fai"]
    assert capsys.readouterr().out == json.dumps(
        {"transcripts": 1, "records": 5, "reference_bases_kept": 164, "gerp_bases": 0, "files": files}
    ) + "\n"

    # A record on a sequence that neither the cache nor the reference holds is still refused.
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("nowhere\t10\tother\tA\tC\t.\t.\t.\n")
    with pytest.raises(SystemExit) as refusal:
        cut.main(["--corpus", str(corpus), "--reference", str(reference), "--force"])
    assert str(refusal.value) == (
        "ERROR: [cut_plugin_regions] variants.vcf names chromosome nowhere, "
        "which the cache does not hold and the reference has no contig for"
    )


def test_a_record_named_by_a_synonym_of_the_reference_contig_is_cut_on_the_cache_chromosome(tmp_path: Path, bgzip_alone):
    corpus, _reference, _ancestor = build_fasta_corpus(tmp_path)
    # The reference names `chr21` the contig the cache's 21 is cut from, and the cache's
    # synonyms file links `alias` to that name alone.
    reference = tmp_path / "prefixed.fa"
    write_fasta(reference, "chr21", CONTIG)
    (corpus / "json_cache" / "chr_synonyms.txt").write_text("alias\tchr21\n", encoding="utf-8")
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("alias\t280\tby_alias\tA\tC\t.\t.\t.\n")
    rc = cut.main(["--corpus", str(corpus), "--reference", str(reference), "--flank", "10", "--record-flank", "5"])
    assert rc == 0

    # The record's flank (275-285) is kept on the one contig of the cut, which carries the
    # cache's name: the reference's contig is written once.
    text = decompressed(corpus / "reference.fa.gz")
    assert [line for line in text.splitlines() if line.startswith(">")] == [">21"]
    kept = set(range(91, 211)) | set(range(245, 256)) | set(range(275, 286))
    assert sequence_of(text) == "".join(CONTIG[position - 1] if position in kept else "N" for position in range(1, 301))
    assert (corpus / "reference.fa.fai").read_text(encoding="utf-8") == "21\t300\t4\t60\t61\n"


def test_a_bigwig_alone_is_cut_around_a_record_on_a_sequence_without_a_transcript(tmp_path: Path, tools):
    corpus, _reference, _ancestor = build_fasta_corpus(tmp_path)
    bedgraph = tmp_path / "gerp.bedgraph"
    bedgraph.write_text(
        "".join(f"21\t{i}\t{i + 1}\t{gerp_value(i)}\n" for i in range(300))
        + "".join(f"unplaced\t{i}\t{i + 1}\t{gerp_value(i)}\n" for i in range(120)),
        encoding="utf-8",
    )
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("21\t300\nunplaced\t120\n", encoding="utf-8")
    gerp = tmp_path / "gerp.bw"
    subprocess.run(["bigtools", "bedgraphtobigwig", str(bedgraph), str(sizes), str(gerp)], check=True)
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("unplaced\t60\tas_written\tA\tC\t.\t.\t.\n")
        handle.write("unplaced\t3\tnear_the_start\tA\tC\t.\t.\t.\n")

    def rows_of(bigwig: Path) -> list[tuple[str, int]]:
        out = tmp_path / "cut.bedgraph"
        subprocess.run(["bigtools", "bigwigtobedgraph", str(bigwig), str(out)], check=True)
        rows = []
        for line in out.read_text(encoding="utf-8").splitlines():
            chrom, start, end, value = line.split("\t")
            for base in range(int(start), int(end)):
                assert float(value) == gerp_value(base)
                rows.append((chrom, base + 1))
        return rows

    # Cut around the records, the bigWig carries the records of `unplaced` plus their flank
    # (55-65, and 1-8 where the flank passes the contig's first base) beside those of 21
    # (110-120, 245-255), each contig at its declared length.
    rc = cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--gerp-regions", "records", "--record-flank", "5"])
    assert rc == 0
    cut_bigwig = corpus / "plugin_data" / "gerp.bw"
    expected = [("21", base) for base in [*range(110, 121), *range(245, 256)]]
    expected += [("unplaced", base) for base in [*range(1, 9), *range(55, 66)]]
    assert sorted(rows_of(cut_bigwig)) == sorted(expected)
    info = subprocess.run(["bigtools", "bigwiginfo", "--chroms", str(cut_bigwig)], check=True, capture_output=True, text=True).stdout
    assert "21 0 300" in info
    assert "unplaced 1 120" in info

    # Cut to the exons, the records of `unplaced` add no region, and the cut holds 21 alone.
    rc = cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--exon-flank", "2", "--force"])
    assert rc == 0
    assert {chrom for chrom, _base in rows_of(cut_bigwig)} == {"21"}

    # A record on a sequence the bigWig does not declare either is refused.
    with (corpus / "variants.vcf").open("a", encoding="utf-8") as handle:
        handle.write("nowhere\t10\tother\tA\tC\t.\t.\t.\n")
    with pytest.raises(SystemExit) as refusal:
        cut.main(["--corpus", str(corpus), "--gerp", str(gerp), "--force"])
    assert str(refusal.value) == (
        "ERROR: [cut_plugin_regions] variants.vcf names chromosome nowhere, "
        "which the cache does not hold and the bigWig has no contig for"
    )


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
