"""The species is an input of the harness: feature_oracle.sh reads it from the case table's
`species` column, generate_reference_output.sh from a corpus manifest's `species`, and both
default to homo_sapiens. The tests drive the two scripts under --dry-run, over fabricated
data trees and over every committed corpus of another species, and read the commands they
plan; nothing runs an engine. Those corpora are found by their manifests, a `species` other
than homo_sapiens with the parity matrix row of its genome under `covers`, so the tests hold
each one the tree carries, but for a corpus whose `flags` hold `--refseq` or `--merged` and
whose `covers` names no such row.

Run: python3 -m pytest scripts/concordance/test_species_harness.py
"""

from __future__ import annotations

import gzip
import json
import shlex
import shutil
import subprocess
from pathlib import Path

import pytest
from reference_release import release_from_provenance

HERE = Path(__file__).resolve().parent
ORACLE = HERE / "feature_oracle.sh"
GENERATOR = HERE / "generate_reference_output.sh"
GOLDEN = HERE.parents[1] / "tests" / "golden" / "116"
HUMAN = "homo_sapiens"
# A parity matrix row of a genome carries this prefix. So does the row of the species
# parameter itself, which a corpus of any species may cover and which names no genome.
ROW_PREFIX = "species-"
HARNESS_ROW = "species-harness-parameter"
# The flags of a corpus cut from a RefSeq or a merged cache: its cache is the slice of the
# species directory suffixed `_refseq` or `_merged`, which
# crates/vep-cli/tests/transcript_sets.rs reads.
SET_FLAGS = ("--refseq", "--merged")
IMAGE = "ensemblorg/ensembl-vep@sha256:" + "a" * 64
COLUMNS = ["case", "tier", "input", "cache", "fasta", "flags", "formats", "perl_image", "perl_mode",
           "reference_release", "masks", "floors", "field_pairs", "species", "unordered", "reference", "refusal",
           "tool"]
HUMAN_PRIMARIES = [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]


def run(script: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(script), *args], check=False, capture_output=True, text=True)


def write_vcf(path: Path, contig: str) -> None:
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    f"{contig}\t1000\t.\tA\tG\t.\tPASS\t.\n", encoding="utf-8")


def stub_vep(path: Path) -> Path:
    path.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def commands_after(marker: str, text: str) -> list[str]:
    """The command lines following every `DRY-RUN <marker>` line of a dry run."""
    lines = text.splitlines()
    return [lines[i + 1] for i, line in enumerate(lines) if line.startswith(f"DRY-RUN {marker}")]


def test_feature_oracle_dry_run_names_the_case_species_for_both_engines(tmp_path: Path, monkeypatch):
    """A yeast row plans both engines under --species saccharomyces_cerevisiae with the Ensembl
    slice under the species directory; a row whose cell is - plans under homo_sapiens."""
    # The runner stops, a dry run included, where its format comparator is not an executable;
    # a dry run starts none, so the plan is read with a stand-in named.
    monkeypatch.setenv("VEP_FORMATS_COMPARATOR", str(stub_vep(tmp_path / "vep-formats-compare")))
    data = tmp_path / "data"
    yeast_slice = data / "caches/perl/vep-cache/saccharomyces_cerevisiae/116_R64-1-1"
    human_slice = data / "caches/perl/vep-cache/homo_sapiens/116_GRCh38"
    for d in (yeast_slice, human_slice, data / "caches/vep-rs/yeast/transcripts", data / "caches/vep-rs/grch38/transcripts"):
        d.mkdir(parents=True)
    write_vcf(data / "yeast.vcf", "IV")
    write_vcf(data / "human.vcf", "21")
    vep = stub_vep(tmp_path / "vep")
    rows = [
        ["c-yeast", "fixture", "yeast.vcf", "116_R64-1-1=caches/vep-rs/yeast", "-", "-", "default", IMAGE,
         "offline", "116.2", "-", "-", "-", "saccharomyces_cerevisiae", "-", "-", "-", "-"],
        ["c-human", "fixture", "human.vcf", "116_GRCh38", "-", "-", "default", IMAGE,
         "offline", "116.2", "-", "-", "-", "-", "-", "-", "-", "-"],
    ]
    table = tmp_path / "cases.tsv"
    table.write_text("\n".join("\t".join(r) for r in [COLUMNS, *rows]) + "\n", encoding="utf-8")
    out = run(ORACLE, "--cases", str(table), "--data-dir", str(data), "--out", str(tmp_path / "never"),
              "--candidates", str(vep), "--candidate-pins", "stub", "--dry-run")
    assert out.returncode == 0, out.stderr
    perl, = commands_after("c-yeast perl default", out.stdout)
    candidate, = commands_after("c-yeast stub default", out.stdout)
    for command in (perl, candidate):
        words = shlex.split(command)
        assert words[words.index("--species") + 1] == "saccharomyces_cerevisiae", command
        assert words[words.index("--assembly") + 1] == "R64-1-1", command
    assert "--cache_version 116 --species saccharomyces_cerevisiae" in perl
    assert f"--json_cache {data}/caches/vep-rs/yeast --species saccharomyces_cerevisiae" in candidate
    slice_lines = [line for line in out.stdout.splitlines() if line.startswith("DRY-RUN c-yeast cache ")]
    assert slice_lines == [f"DRY-RUN c-yeast cache species=saccharomyces_cerevisiae perl={yeast_slice} candidate={data}/caches/vep-rs/yeast"]
    assert str(yeast_slice).endswith("/saccharomyces_cerevisiae/116_R64-1-1")
    for command in (*commands_after("c-human perl default", out.stdout), *commands_after("c-human stub default", out.stdout)):
        assert "--species homo_sapiens --assembly GRCh38" in command, command
    assert f"DRY-RUN c-human cache species=homo_sapiens perl={human_slice} " in out.stdout
    # The slice is looked for under the species directory: a yeast row whose slice is absent
    # fails preflight by that path.
    yeast_slice.rmdir()
    out = run(ORACLE, "--cases", str(table), "--data-dir", str(data), "--out", str(tmp_path / "never"),
              "--candidates", str(vep), "--candidate-pins", "stub", "--dry-run", "--only", "c-yeast")
    assert out.returncode != 0
    assert f"c-yeast: Ensembl VEP cache slice not found: {yeast_slice}" in out.stderr


def _corpus(tmp_path: Path, manifest: dict) -> Path:
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    write_vcf(corpus / "variants.vcf", "1")
    with gzip.open(corpus / "reference.fa.gz", "wt") as fh:
        fh.write(">1\nACGT\n")
    (corpus / "reference.fa.fai").write_text("1\t4\t3\t4\t5\n", encoding="utf-8")
    (corpus / "manifest.json").write_text(json.dumps({"flags": ["--protein"], "fasta": "reference.fa.gz", **manifest}),
                                          encoding="utf-8")
    return corpus


def test_generate_reference_output_dry_run_reads_the_manifest_species_and_fasta(tmp_path: Path):
    """An arabidopsis corpus plans its four commands under --species arabidopsis_thaliana
    --assembly TAIR10 --cache_version 63 with the manifest's FASTA, read from the slice under
    the species directory; a manifest without species plans under homo_sapiens. A corpus's
    FASTA is the manifest's, as the golden tests pass it, so --fasta-for is a dataset option."""
    cache = tmp_path / "cache"
    (cache / "arabidopsis_thaliana/63_TAIR10/1").mkdir(parents=True)
    for c in HUMAN_PRIMARIES:
        (cache / "homo_sapiens/116_GRCh38" / c).mkdir(parents=True)
    corpus = _corpus(tmp_path, {"species": "arabidopsis_thaliana", "assembly": "TAIR10", "release": "63"})
    out = run(GENERATOR, "--corpus", str(corpus), "--perl-cache-dir", str(cache), "--dry-run")
    assert out.returncode == 0, out.stderr
    commands = [line for line in out.stdout.splitlines() if line.startswith("docker run ")]
    assert len(commands) == 4
    for command in commands:
        words = shlex.split(command)
        assert words[words.index("--species") + 1] == "arabidopsis_thaliana", command
        assert words[words.index("--assembly") + 1] == "TAIR10", command
        assert words[words.index("--cache_version") + 1] == "63", command
        assert words[words.index("--fasta") + 1] == "/work/ref.fa", command
        mounts = [words[i + 1] for i, w in enumerate(words) if w == "-v"]
        assert any(m.startswith(f"{corpus}/.generate.dry-run/reference.fa:/work/ref.fa:") for m in mounts), mounts
    assert f"perl cache check (arabidopsis_thaliana TAIR10): 1 chromosome directories under {cache}/arabidopsis_thaliana/63_TAIR10" in out.stdout
    out = run(GENERATOR, "--corpus", str(corpus), "--perl-cache-dir", str(cache), "--dry-run",
              "--fasta-for", "arabidopsis_thaliana:TAIR10=/x.fa")
    assert out.returncode != 0 and "drop --fasta-for" in out.stderr

    shutil.rmtree(corpus)
    corpus = _corpus(tmp_path, {"assembly": "GRCh38", "release": "116"})
    out = run(GENERATOR, "--corpus", str(corpus), "--perl-cache-dir", str(cache), "--dry-run")
    assert out.returncode == 0, out.stderr
    commands = [line for line in out.stdout.splitlines() if line.startswith("docker run ")]
    assert len(commands) == 4
    for command in commands:
        assert "--cache_version 116 --species homo_sapiens --assembly GRCh38 " in command, command
    assert f"all 25 primary contigs present under {cache}/homo_sapiens/116_GRCh38" in out.stdout


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _transcript_set(manifest: dict) -> str | None:
    """The transcript set a manifest's `flags` name, `refseq` or `merged`; none for a corpus
    cut from the cache of Ensembl's own transcripts."""
    named = [flag.removeprefix("--") for flag in manifest.get("flags") or [] if flag in SET_FLAGS]
    if len(named) > 1:
        raise ValueError(f"the flags name {len(named)} transcript sets: {', '.join(named)}")
    return named[0] if named else None


def _cache_slice(manifest: dict) -> str:
    """The slice of an Ensembl VEP cache root a manifest's run reads: the species directory,
    suffixed `_refseq` or `_merged` where the flags name that set, then `<release>_<assembly>`."""
    transcript_set = _transcript_set(manifest)
    directory = manifest["species"] + (f"_{transcript_set}" if transcript_set else "")
    return f"{directory}/{manifest['release']}_{manifest['assembly']}"


def _species_corpora(golden: Path) -> dict[str, list[Path]]:
    """The corpora under `golden` whose manifest names a species other than homo_sapiens, by
    the parity matrix row of its genome that each one's own `covers` names: a name beginning
    `species-` other than HARNESS_ROW. Such a corpus that names no row of its genome is
    refused by its name, unless its manifest's `flags` hold one of SET_FLAGS: a transcript-set
    corpus covers the rows of its flag and joins no row here."""
    corpora: dict[str, list[Path]] = {}
    for manifest in sorted(golden.glob("*/manifest.json")):
        data = _load(manifest)
        species = data.get("species") or HUMAN
        if species == HUMAN:
            continue
        rows = [row for row in data.get("covers") or [] if row.startswith(ROW_PREFIX) and row != HARNESS_ROW]
        if not rows and any(flag in SET_FLAGS for flag in data.get("flags") or []):
            continue
        if not rows:
            raise ValueError(f"{manifest.parent.name}: the manifest names the species {species} and its covers holds "
                             f"no row of that genome, a name beginning {ROW_PREFIX} other than {HARNESS_ROW}")
        for row in rows:
            corpora.setdefault(row, []).append(manifest.parent)
    return corpora


def test_species_corpora_name_their_species_cache_fasta_and_reference(tmp_path: Path):
    """Each genome's row is covered by one committed corpus. Its manifest names a species, an
    assembly and a cache version; the cache it carries and the record of its Ensembl VEP
    reference name the same three, at the reference release of the human corpus; and its
    reference generation is planned under them with the corpus's own input and FASTA, read
    from the slice under the species directory. Where the manifest's flags name a transcript
    set, that directory is the one suffixed with the set, the corpus's cache names the set and
    each planned command carries the flag."""
    human_release = release_from_provenance(_load(GOLDEN / "GRCh38" / "provenance.json"))
    rows = _species_corpora(GOLDEN)
    assert rows, f"{GOLDEN} holds no corpus of a species other than {HUMAN}"
    for row, found in rows.items():
        assert len(found) == 1, f"{row}: covered by {[c.name for c in found]}, not by one corpus"
        corpus, = found
        manifest = _load(corpus / "manifest.json")
        species, assembly, release = (manifest.get(k) for k in ("species", "assembly", "release"))
        assert assembly, f"{row}: {corpus.name} names the assembly {assembly!r}"
        assert release and release.isdigit(), f"{row}: release is {release!r}"
        fasta = manifest.get("fasta")
        assert fasta and fasta.endswith(".gz"), f"{row}: fasta is {fasta!r}"
        plain_fasta = fasta.removesuffix(".gz")
        transcript_set = _transcript_set(manifest)
        cache_slice = _cache_slice(manifest)

        info = _load(corpus / "json_cache" / "info.json")
        assert (info["species"], info["assembly"], info["cache_version"]) == (species, assembly, int(release)), \
            f"{row}: {corpus.name}/json_cache is the cache of {info['species']} {info['assembly']} {info['cache_version']}"
        if transcript_set:
            assert info.get("source_type") == transcript_set, \
                f"{row}: {corpus.name}/json_cache names the transcript set {info.get('source_type')!r}, the flags {transcript_set!r}"

        record = _load(corpus / "provenance.json")
        assert release_from_provenance(record) == human_release, f"{row}: {corpus.name} records another reference release"
        assert record["perl_cache"] == cache_slice, f"{row}: perl_cache is {record['perl_cache']!r}"
        assert record["cache_version"] == int(release), f"{row}: cache_version is {record['cache_version']!r}"
        recorded = shlex.split(record["vep_command"])
        expected = {"--species": species, "--assembly": assembly, "--cache_version": release}
        for flag, value in {**expected, "--fasta": plain_fasta}.items():
            assert recorded[recorded.index(flag) + 1] == value, f"{row}: recorded {record['vep_command']}"

        cache = tmp_path / row
        (cache / cache_slice / "1").mkdir(parents=True)
        out = run(GENERATOR, "--corpus", str(corpus), "--perl-cache-dir", str(cache), "--dry-run", "--force")
        assert out.returncode == 0, f"{row}: {out.stderr}"
        assert f"perl cache check ({species} {assembly}): 1 chromosome directories under {cache}/{cache_slice};" in out.stdout
        commands = [line for line in out.stdout.splitlines() if line.startswith("docker run ")]
        assert len(commands) == 4, f"{row}: {len(commands)} commands planned"
        for command in commands:
            words = shlex.split(command)
            for flag, value in {**expected, "--fasta": "/work/ref.fa"}.items():
                assert words[words.index(flag) + 1] == value, f"{row}: {command}"
            if transcript_set:
                assert f"--{transcript_set}" in words, f"{row}: {command}"
            mounts = [words[i + 1] for i, w in enumerate(words) if w == "-v"]
            assert f"{corpus}/variants.vcf:/work/input.vcf:ro" in mounts, f"{row}: {mounts}"
            assert f"{corpus}/.generate.dry-run/{plain_fasta}:/work/ref.fa:ro" in mounts, f"{row}: {mounts}"


def _manifest_only(golden: Path, name: str, manifest: dict) -> Path:
    """A corpus directory under `golden` that holds a manifest and nothing else."""
    corpus = golden / name
    corpus.mkdir(parents=True)
    (corpus / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return corpus


def test_a_corpus_of_a_further_species_joins_the_rows_by_its_manifest_alone(tmp_path: Path):
    """Beside the manifests of the committed corpora of other species, a corpus of a species
    the tree carries no corpus of, the Japanese medaka, is read under the row of its genome,
    and under that row alone though it covers the row of the species parameter too. A human
    corpus, by default or by name, adds no row whatever it covers."""
    committed = _species_corpora(GOLDEN)
    for corpus in {corpus for found in committed.values() for corpus in found}:
        (tmp_path / corpus.name).mkdir()
        shutil.copy(corpus / "manifest.json", tmp_path / corpus.name)
    medaka = _manifest_only(tmp_path, "ASM223467v1-medaka", {
        "species": "oryzias_latipes", "assembly": "ASM223467v1", "release": "116",
        "covers": ["genetic-code-tables", "species-medaka-asm223467v1", HARNESS_ROW]})
    _manifest_only(tmp_path, "GRCh38-by-default", {"assembly": "GRCh38", "covers": [HARNESS_ROW, "species-human-grch38"]})
    _manifest_only(tmp_path, "GRCh38-by-name", {"species": HUMAN, "assembly": "GRCh38", "covers": ["species-human-grch38"]})
    assert _species_corpora(tmp_path) == {
        **{row: [tmp_path / corpus.name for corpus in found] for row, found in committed.items()},
        "species-medaka-asm223467v1": [medaka],
    }


def test_a_corpus_of_another_species_that_covers_no_row_of_its_genome_is_refused_by_name(tmp_path: Path):
    """A manifest that names a species other than homo_sapiens and whose `covers` is absent,
    holds no name beginning `species-`, or holds the row of the species parameter alone, stops
    the derivation with the corpus's name."""
    for case, covers in enumerate((None, ["genetic-code-tables"], [HARNESS_ROW, "genetic-code-tables"])):
        manifest = {"species": "danio_rerio", "assembly": "GRCz11", "release": "116"}
        if covers is not None:
            manifest["covers"] = covers
        _manifest_only(tmp_path / str(case), "GRCz11-zebrafish", manifest)
        with pytest.raises(ValueError) as refused:
            _species_corpora(tmp_path / str(case))
        assert str(refused.value) == (
            "GRCz11-zebrafish: the manifest names the species danio_rerio and its covers holds no row of that genome, "
            "a name beginning species- other than species-harness-parameter"
        ), covers


def test_a_transcript_set_corpus_of_another_species_joins_no_row_unless_it_names_one(tmp_path: Path):
    """A manifest that names a species other than homo_sapiens and whose `flags` hold
    `--refseq` or `--merged` is read without a row of its genome: it joins no row, whether its
    `covers` holds the row of the species parameter or not. One that does name a row of its
    genome joins that row, set flag or not, and the same `covers` without the flag is refused."""
    mouse = {"species": "mus_musculus", "assembly": "GRCm39", "release": "116"}
    _manifest_only(tmp_path, "GRCm39-mouse-refseq", {**mouse, "flags": ["--refseq"], "covers": ["refseq"]})
    _manifest_only(tmp_path, "GRCm39-mouse-merged", {**mouse, "flags": ["--merged", "--hgvs"], "covers": ["merged", HARNESS_ROW]})
    cattle = _manifest_only(tmp_path, "ARS-UCD2.0-cattle-merged", {
        "species": "bos_taurus", "assembly": "ARS-UCD2.0", "release": "116", "flags": ["--merged"],
        "covers": ["merged", "species-cattle-ars-ucd2"]})
    assert _species_corpora(tmp_path) == {"species-cattle-ars-ucd2": [cattle]}

    _manifest_only(tmp_path / "plain", "GRCm39-mouse-refseq", {**mouse, "flags": ["--hgvs"], "covers": ["refseq"]})
    with pytest.raises(ValueError) as refused:
        _species_corpora(tmp_path / "plain")
    assert str(refused.value).startswith("GRCm39-mouse-refseq: the manifest names the species mus_musculus and ")


def test_a_species_corpus_of_a_transcript_set_reads_the_suffixed_species_directory(tmp_path: Path):
    """A manifest whose `flags` hold `--merged` or `--refseq` reads the slice under the species
    directory suffixed with the set, and one without either the species directory itself. The
    reference generation of a merged corpus is planned from that suffixed slice, each command
    under `--merged` and `--species` the species, and stops where the cache root holds the
    slice of Ensembl's own transcripts alone."""
    cattle = {"species": "bos_taurus", "assembly": "ARS-UCD2.0", "release": "116"}
    assert _cache_slice({**cattle, "flags": ["--merged", "--hgvs"]}) == "bos_taurus_merged/116_ARS-UCD2.0"
    assert _cache_slice({**cattle, "flags": ["--refseq"]}) == "bos_taurus_refseq/116_ARS-UCD2.0"
    assert _cache_slice({**cattle, "flags": ["--hgvs"]}) == "bos_taurus/116_ARS-UCD2.0"
    assert _cache_slice(cattle) == "bos_taurus/116_ARS-UCD2.0"
    with pytest.raises(ValueError) as refused:
        _cache_slice({**cattle, "flags": ["--refseq", "--merged"]})
    assert str(refused.value) == "the flags name 2 transcript sets: refseq, merged"

    manifest = {**cattle, "flags": ["--merged"]}
    corpus = _corpus(tmp_path, manifest)
    cache = tmp_path / "cache"
    (cache / "bos_taurus/116_ARS-UCD2.0/1").mkdir(parents=True)
    out = run(GENERATOR, "--corpus", str(corpus), "--perl-cache-dir", str(cache), "--dry-run")
    assert out.returncode != 0
    assert f"{cache}/{_cache_slice(manifest)}" in out.stderr, out.stderr
    (cache / _cache_slice(manifest) / "1").mkdir(parents=True)
    out = run(GENERATOR, "--corpus", str(corpus), "--perl-cache-dir", str(cache), "--dry-run")
    assert out.returncode == 0, out.stderr
    assert f"perl cache check (bos_taurus ARS-UCD2.0): 1 chromosome directories under {cache}/{_cache_slice(manifest)};" in out.stdout
    commands = [line for line in out.stdout.splitlines() if line.startswith("docker run ")]
    assert len(commands) == 4
    for command in commands:
        words = shlex.split(command)
        assert "--merged" in words, command
        assert words[words.index("--species") + 1] == "bos_taurus", command
