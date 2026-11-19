"""The species is an input of the harness: feature_oracle.sh reads it from the case table's
`species` column, generate_reference_output.sh from a corpus manifest's `species`, and both
default to homo_sapiens. The tests drive the two scripts under --dry-run, over fabricated
data trees and over the committed corpora of the four other species, and read the commands
they plan; nothing runs an engine.

Run: python3 -m pytest scripts/concordance/test_species_harness.py
"""

from __future__ import annotations

import gzip
import json
import shlex
import shutil
import subprocess
from pathlib import Path

from reference_release import release_from_provenance

HERE = Path(__file__).resolve().parent
ORACLE = HERE / "feature_oracle.sh"
GENERATOR = HERE / "generate_reference_output.sh"
GOLDEN = HERE.parents[1] / "tests" / "golden" / "116"
# The parity matrix row of each of the four other species, with the species and assembly
# its corpus is cut for.
SPECIES_ROWS = {
    "species-arabidopsis-tair10": ("arabidopsis_thaliana", "TAIR10"),
    "species-fly-bdgp6": ("drosophila_melanogaster", "BDGP6.54"),
    "species-mouse-grcm39": ("mus_musculus", "GRCm39"),
    "species-yeast-r64": ("saccharomyces_cerevisiae", "R64-1-1"),
}
IMAGE = "ensemblorg/ensembl-vep@sha256:" + "a" * 64
COLUMNS = ["case", "tier", "input", "cache", "fasta", "flags", "formats", "perl_image", "perl_mode",
           "reference_release", "masks", "floors", "field_pairs", "species", "unordered", "reference", "refusal"]
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


def test_feature_oracle_dry_run_names_the_case_species_for_both_engines(tmp_path: Path):
    """A yeast row plans both engines under --species saccharomyces_cerevisiae with the Ensembl
    slice under the species directory; a row whose cell is - plans under homo_sapiens."""
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
         "offline", "116.2", "-", "-", "-", "saccharomyces_cerevisiae", "-", "-", "-"],
        ["c-human", "fixture", "human.vcf", "116_GRCh38", "-", "-", "default", IMAGE,
         "offline", "116.2", "-", "-", "-", "-", "-", "-", "-"],
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
    FASTA is the manifest's, as the golden tests pass it, so --fasta-for is a suite option."""
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


def _species_corpora() -> dict[str, list[Path]]:
    """The committed corpora whose manifest's own `covers` names each row of SPECIES_ROWS."""
    corpora: dict[str, list[Path]] = {row: [] for row in SPECIES_ROWS}
    for manifest in sorted(GOLDEN.glob("*/manifest.json")):
        for row in _load(manifest).get("covers") or []:
            if row in corpora:
                corpora[row].append(manifest.parent)
    return corpora


def test_species_corpora_name_their_species_cache_fasta_and_reference(tmp_path: Path):
    """Each row of SPECIES_ROWS is covered by one committed corpus. Its manifest names the
    row's species and assembly and a cache version; the cache it carries and the record of
    its Ensembl VEP reference name the same three, at the reference release of the human
    corpus; and its reference generation is planned under them with the corpus's own input
    and FASTA, read from the slice under the species directory."""
    human_release = release_from_provenance(_load(GOLDEN / "GRCh38" / "provenance.json"))
    for row, found in _species_corpora().items():
        assert len(found) == 1, f"{row}: covered by {[c.name for c in found]}, not by one corpus"
        corpus, = found
        manifest = _load(corpus / "manifest.json")
        species, assembly, release = (manifest.get(k) for k in ("species", "assembly", "release"))
        assert (species, assembly) == SPECIES_ROWS[row], f"{row}: {corpus.name} names {species!r} {assembly!r}"
        assert release and release.isdigit(), f"{row}: release is {release!r}"
        fasta = manifest.get("fasta")
        assert fasta and fasta.endswith(".gz"), f"{row}: fasta is {fasta!r}"
        plain_fasta = fasta.removesuffix(".gz")
        cache_slice = f"{species}/{release}_{assembly}"

        info = _load(corpus / "json_cache" / "info.json")
        assert (info["species"], info["assembly"], info["cache_version"]) == (species, assembly, int(release)), \
            f"{row}: {corpus.name}/json_cache is the cache of {info['species']} {info['assembly']} {info['cache_version']}"

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
            mounts = [words[i + 1] for i, w in enumerate(words) if w == "-v"]
            assert f"{corpus}/variants.vcf:/work/input.vcf:ro" in mounts, f"{row}: {mounts}"
            assert f"{corpus}/.generate.dry-run/{plain_fasta}:/work/ref.fa:ro" in mounts, f"{row}: {mounts}"

