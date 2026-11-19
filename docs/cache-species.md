# Other species

Any species with an Ensembl or Ensembl Genomes cache converts and annotates the same way as human ([Cache setup](cache-setup.md), [Cache conversion](cache-conversion.md)).

In each such cache:

- the species directory is the Latin name;
- each transcript's genetic code is its record's `codon_table`;
- the chromosome names are the cache's.

An Ensembl Genomes cache carries that division's own release number, so `--cache_version` names it rather than the Ensembl release.

vep-rs is verified by concordance against Ensembl VEP 116.2 on human GRCh37 and GRCh38, the whole-genome datasets of the current release, and on mouse GRCm39, fruit fly BDGP6.54, Arabidopsis thaliana TAIR10 and yeast R64-1-1.

Each of the four has a golden corpus of its own under `tests/golden/116/`, cut from a conversion of the cache Ensembl ships. A larger comparison dataset stands behind each corpus. That dataset is one chromosome with the organelle genomes for mouse, the fruit fly and Arabidopsis, and every chromosome for yeast. On each of them vep-rs matches Ensembl VEP exactly except for the documented divergences ([Intended divergences](intended-divergences.md)).

`--sift` and `--polyphen` refuse on a cache without their data, and `--af_gnomad` on a cache that declares variation columns without the gnomAD exome column; every other flag that reads data Ensembl ships for human alone runs there as [the flag table](cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone) states.

| Species and assembly | Cache release and `--cache_version` | Genetic codes | Corpus |
| --- | --- | --- | --- |
| `mus_musculus` GRCm39 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/GRCm39-mouse` |
| `drosophila_melanogaster` BDGP6.54 | Ensembl 116; `--cache_version 116`, the default | `mitochondrion_genome` table 5; every other region table 1 | `tests/golden/116/BDGP6.54-fly` |
| `arabidopsis_thaliana` TAIR10 | Ensembl Plants 63; `--cache_version 63` | `Pt` table 11; `Mt` table 1, which the cache states; the nuclear chromosomes table 1 | `tests/golden/116/TAIR10-arabidopsis` |
| `saccharomyces_cerevisiae` R64-1-1 | Ensembl 116; `--cache_version 116`, the default | `Mito` table 3; every other region table 1 | `tests/golden/116/R64-1-1-yeast` |

The genetic codes are the tables the cache names per sequence region. They are read from the core database's `codon_table` attribute by Ensembl's dumper and carried by the conversion. The translation tables themselves are `Bio::Tools::CodonTable`'s ([The codon table of each sequence region](cache-setup.md#the-codon-table-of-each-sequence-region)). The four genomes exercise tables 1, 2, 3, 5 and 11.

Any other species in Ensembl release 116 or Ensembl Genomes release 63 with a GFF3 and a peptide FASTA, or with an Ensembl VEP cache to convert, builds a cache ([Build a cache](cache-setup.md#build-a-cache)). Such a species annotates under the codon tables the cache names, without a concordance measurement. A GFF3 and a FASTA state no genetic code. A build from them names an organelle genome's table with `--codon-table`, or reads it from the release's core dumps.

## Arabidopsis thaliana

Arabidopsis thaliana at Ensembl Plants release 63, whose plastid `Pt` translates under genetic code table 11 and whose mitochondrion `Mt` under table 1:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$VEP_CACHE/arabidopsis_thaliana/63_TAIR10:/input/63_TAIR10:ro" \
  -v /path/to/arabidopsis_cache:/output -v "$PWD/scripts/data:/scripts:ro" \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl /input/63_TAIR10 /output

# Annotate; the cache version is the Plants release
vep -i arabidopsis.vcf -o arabidopsis.txt --offline \
  --json_cache /path/to/arabidopsis_cache \
  --species arabidopsis_thaliana --assembly TAIR10 --cache_version 63
```

The cache names the chromosomes `1` to `5`, `Mt` and `Pt`. Its synonyms resolve the `Chr1` to `Chr5` spellings and the GenBank and RefSeq accessions, so an input spelled any of those ways annotates. Among the accessions are `CP002687.1` and `NC_003075.7` for chromosome 4, `Y08501.2` and `NC_001284.2` for `Mt`, and `AP000423.1` and `NC_000932.1` for `Pt`. Under `--dir_cache` the cache sits at `arabidopsis_thaliana/63_TAIR10`, which `--cache_version 63` selects. `--json_cache` names the directory itself and needs no version.

The Arabidopsis cache carries SIFT matrices, no PolyPhen matrices and variants on the five nuclear chromosomes alone. So `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and `--check_existing` prints empty columns on `Mt` and `Pt`. The corpus `tests/golden/116/TAIR10-arabidopsis` compares both engines on chromosome 4, `Mt` and `Pt` under the default flags, `--everything` and `--hgvs`. The plastid's start codon shapes are among its records.

A cache built with `vep-cache-builder` from the release's GFF3 and protein FASTA sets aside the FASTA records of four mitochondrial transcripts. The transcripts are `ATMG00410.1`, `ATMG00730.1`, `ATMG01080.1` and `ATMG01170.1`. Each record's length is neither its transcript's whole codons nor one fewer. With `--genome-fasta` the build translates the four peptides from the genome; without it the four transcripts carry no peptide.

## Drosophila melanogaster

The fruit fly at release 116, whose cache names the mitochondrion `mitochondrion_genome` under genetic code table 5 while its variant sets say `MT`:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$VEP_CACHE/drosophila_melanogaster/116_BDGP6.54:/input/116_BDGP6.54:ro" \
  -v /path/to/fly_cache:/output -v "$PWD/scripts/data:/scripts:ro" \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl /input/116_BDGP6.54 /output

# Annotate; the synonyms file links the input's MT to the cache's name
printf 'mitochondrion_genome\tMT\n' > synonyms.txt
vep -i fly.vcf -o fly.txt --offline --json_cache /path/to/fly_cache \
  --species drosophila_melanogaster --assembly BDGP6.54 \
  --synonyms synonyms.txt
```

The cache names the chromosome arms `2L`, `2R`, `3L` and `3R`, the chromosomes `4`, `X` and `Y`, `mitochondrion_genome`, `rDNA` and 20 scaffolds (`211000022278760`, `Unmapped_Scaffold_8_D1580_D1567`). Its synonyms resolve the `chr`-prefixed names, `chrM` and the GenBank and RefSeq accessions (`AE014134.6` and `NT_033779.5` for `2L`, `KJ947872.2` and `NC_024511.2` for the mitochondrion). `MT` is the name Ensembl's variant files give the fly's mitochondrion, and it is neither a chromosome nor a synonym of the cache. So an input spelling it that way is skipped in both engines unless a `--synonyms` file links the two names. The flag replaces the cache's own table, so the file also lists any other synonym the input needs.

The fly cache has no SIFT or PolyPhen matrices and no variant. So `--everything` runs with both disabled, as in Ensembl VEP, and `--check_existing` prints empty columns. The corpus `tests/golden/116/BDGP6.54-fly` compares both engines on 2L, the mitochondrion, 4, Y, rDNA and seven of the scaffolds under the default flags, `--everything` and `--hgvs`.

## Mus musculus

The mouse at release 116, whose cache carries SIFT matrices, a regulatory build and the mitochondrion as `MT` under genetic code table 2:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$VEP_CACHE/mus_musculus/116_GRCm39:/input/116_GRCm39:ro" \
  -v /path/to/mouse_cache:/output -v "$PWD/scripts/data:/scripts:ro" \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl /input/116_GRCm39 /output

# Annotate with SIFT predictions and known variants
vep -i mouse.vcf -o mouse.txt --offline --json_cache /path/to/mouse_cache \
  --species mus_musculus --assembly GRCm39 \
  --fasta Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --sift b --check_existing
```

The cache names the chromosomes `1` to `19`, `X`, `Y` and `MT` and 28 unplaced scaffolds by accession (`GL456210.1`, `JH584296.1`). Its synonyms resolve the `chr`-prefixed names, `chrM` and the RefSeq and GenBank accessions (`NC_000085.7` and `CM001012.3` for chromosome 19, `NC_005089.1` for `MT`). Under `--dir_cache` the cache sits at `mus_musculus/116_GRCm39`, which the default `--cache_version` selects.

The mouse cache has SIFT matrices and no PolyPhen matrix: `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and SIFT populated.

Some data is in the mouse cache alone, so the flags that read it print values here and empty cells on the other three ([the flag table](cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone)):

| Data | Caches of the four that carry it |
| --- | --- |
| The TSL, CCDS and GENCODE attributes, gene phenotypes and variant citations | the mouse cache alone |
| APPRIS | the mouse and fruit fly caches |
| UniProt accessions and protein domains | all four |
| Variant synonyms | the mouse and Arabidopsis caches |

The corpus `tests/golden/116/GRCm39-mouse` compares both engines on chromosome 19, the mitochondrion and the unplaced scaffolds under the default flags, `--everything` and `--hgvs --sift b`.

## Saccharomyces cerevisiae

Yeast at release 116, whose cache names the chromosomes `I` to `XVI` and `Mito`, the mitochondrion under genetic code table 3, and carries SIFT matrices and the SGRP variants:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$VEP_CACHE/saccharomyces_cerevisiae/116_R64-1-1:/input/116_R64-1-1:ro" \
  -v /path/to/yeast_cache:/output -v "$PWD/scripts/data:/scripts:ro" \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl /input/116_R64-1-1 /output

# Annotate; the SGRP identifiers come back through --check_existing
vep -i yeast.vcf -o yeast.txt --offline --json_cache /path/to/yeast_cache \
  --species saccharomyces_cerevisiae --assembly R64-1-1 \
  --check_existing --sift b
```

The cache's synonyms are the GenBank accessions alone (`BK006935.2` for `I`, `AJ011856.1` for `Mito`). A `chrI` spelling is reached through the `chr` rule. A mitochondrion spelled `MT`, `M` or `chrM` is skipped unless a `--synonyms` file links it to `Mito`. Under `--dir_cache` the cache sits at `saccharomyces_cerevisiae/116_R64-1-1`, which the default `--cache_version` selects.

Under table 3 a `CTN` codon reads as threonine, `TGA` as tryptophan and `ATA` as a start codon. So a mitochondrial change the nuclear code calls missense can be synonymous, and a stop codon changed to `TGA` is a stop lost. The cache has no PolyPhen matrix: `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and the SIFT column populated. The corpus `tests/golden/116/R64-1-1-yeast` compares both engines on the mitochondrion and twelve nuclear chromosomes under the default flags, `--everything` and `--hgvs`.

## Where next

- [Cache setup](cache-setup.md): building a cache, naming it, the reference FASTA and the chromosome names.
- [Cache conversion](cache-conversion.md): what a converted Ensembl VEP cache carries.
- [The flag table](cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone): what each flag that reads data Ensembl ships for human alone does on another species.
- [Intended divergences](intended-divergences.md): where the two engines differ, and why.
