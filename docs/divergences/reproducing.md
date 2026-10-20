# Reproducing these examples

Every example on the defect pages is one or two VCF records, from a public dataset or hand-written and printed in full in its section, run through both engines with their default consequence settings.

## The engines

Every command on the defect pages and on this one writes both binaries as `vep`. vep-rs's `vep` is a release binary or a source build placed on PATH (`cargo build --release` leaves it at `target/release/vep`). Ensembl VEP's `vep` is the name it has inside the `ensemblorg/ensembl-vep:release_116.2` image (or a local install of release 116.2). The text beside each command names the engine it runs.

## The inputs

All are fetched and laid out by `scripts/data/download_real_world_vcfs.sh`, which pins the ClinVar release:

- ClinVar's GRCh37 and GRCh38 VCFs (`clinvar_grch37.vcf.gz`, `clinvar_grch38.vcf.gz`);
- gnomAD-SV v2.1 (GRCh37, `gnomad_sv_v2.1_chr21.vcf.gz`) and v4.1 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), chromosome 21;
- the 1000 Genomes structural-variant releases, chromosome 21 (`1kg_sv_chr21.vcf.gz` per assembly: Phase 3 for GRCh37, high-coverage for GRCh38);
- ClinVar's structural variants, chromosome 21 (`clinvar_sv_chr21.vcf.gz`);
- the synthetic sets committed under `tests/sv_validation/<assembly>/`, of which the defect pages use `02_mnp_complex`, `05_symbolic_del_ins`, `07_cnv_repeat`, `09_breakends`, `10_special_alleles`, `11_multi_allelic`, `12_complex_imprecise` and `13_vcf45_features` (each `.vcf.gz`).

Contig names are Ensembl's (`21`, not `chr21`; `scripts/README.md` describes the canonicalisation).

## The transcripts

The transcript sets are Ensembl release 116's: Ensembl VEP reads its own cache (`homo_sapiens/116_GRCh37` and `homo_sapiens/116_GRCh38`), and vep-rs reads the JSON conversion of the same cache (`scripts/data/storable_to_json.pl`, described in [`docs/cache-setup.md`](../cache-setup.md), laid out as `json_cache/homo_sapiens/116_GRCh37` and `json_cache/homo_sapiens/116_GRCh38`), so both engines see identical transcripts. Transcript coordinates quoted on the defect pages are read from that cache. The worked examples of classes 13 to 22 run vep-rs on the pruned slice of the GRCh38 cache committed with the golden corpus, `tests/golden/116/GRCh38-release116/json_cache`, which holds the transcripts the corpus's records reach.

## The reference genome

Both engines take the assembly's Ensembl primary-assembly FASTA on every example, the SNP/indel ones included: `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl release 116) and `Homo_sapiens.GRCh37.75.dna.primary_assembly.fa` (Ensembl's GRCh37 primary assembly, first served under release 75 and unchanged since), each with its `.fai` index. The reference bases quoted on the defect pages come from them; this command gives the 674 bases of [class 3's](defect-b.md#3-start_lost-and-start_retained_variant-emitted-together) explicit spelling:

```
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 21:26965974-26966647
```

## Writing a probe input

Every probe input on the defect pages is written as its five leading VCF columns, CHROM POS ID REF ALT, tab-separated. To run one, put it in a file after the two header lines `##fileformat=VCFv4.2` and `#CHROM POS ID REF ALT QUAL FILTER INFO` (tab-separated). Use `.` in QUAL and FILTER and, in INFO, `.` or the keys the class shows (`END=`, `SVTYPE=` and the like); `--format vcf` needs all eight columns. The `alone` inputs are such a file with one record, the `plus one SNV` inputs the same file with the SNV line added.

## The commands

vep-rs:

```
vep -i record.vcf -o vep-rs.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/116_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38 \
  --fork 16 --buffer_size 5000 --force --no_stats
```

Ensembl VEP 116.2, with the flags the reference outputs were generated with:

```
docker run --rm -v "$PWD":/work -v /path/to/vep_cache:/cache -v /path/to/fasta:/fasta ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/record.vcf -o /work/vep.txt --offline --cache \
  --dir_cache /cache --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --format vcf \
  --fasta /fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --buffer_size 5000 --fork 4 --no_stats --force_overwrite
```

Substitute `GRCh37` and its cache and FASTA for the GRCh37 examples.

## The order-dependent class

The reference outputs behind the counts used `--fork 4` for the structural-variant sets and `--fork 16` for the SNP/indel datasets. Two classes depend on the other records of a buffer: [class 23](defect-p.md#23-a-bracket-breakends-own-end-rows-lost-to-the-first-transcript-that-reaches-it) on the order in which the buffer's records ask for their regions, and [class 14](defect-g.md#14-a-coordinate-at-or-past-2147478648-a-breakend-mate-a-span-end-or-the-records-own-position-turns-the-references-buffer-intergenic) on the records loaded beside the one whose coordinate reaches 2,147,478,648, which the turn makes intergenic too. Their pages say which inputs give which rows. No neighbouring record, `--fork` or `--buffer_size` value changes the outcome of any other class.

Flags that do change an outcome are named in the sections that use them:

- `--max_sv_size` (default 10,000,000 on both engines; `-1` lifts the cap on both) decides [class 5's](defect-c.md#5-a-record-ensembl-vep-drops-before-annotation) oversize route, and [class 12](cnv-tr.md#12-cnvtr-read-symbolically-not-expanded) lowers it on Ensembl VEP alone;
- `--dont_skip` lets Ensembl VEP annotate a record it would drop ([5.4](defect-c.md#54-why-ensembl-vep-writes-nothing));
- `--fasta`, which the reference protocol passes, decides whether Ensembl VEP expands a `<CNV:TR>` record (classes 12, [19](defect-l.md#19-an-empty-or-unchanged-cnvtr-expansion-dropped-without-a-row) and [22](defect-o.md#22-a-real-cnvtr-allele-dropped-beside-an-rn0-allele));
- `--distance` moves the coordinate at which class 14's turn happens ([14.6](defect-g.md#146-try-it-yourself)).

Both engines write their warnings to `<output>_warnings.txt` beside the output (`vep.txt_warnings.txt` for the reference runs), and vep-rs echoes each on standard error unless `--quiet`.

## The comparators

The comparators are `scripts/concordance/compare_vep_outputs.py` (SNP/indel) and `scripts/validation/compare_sv_concordance.py` (structural variants); their per-class counters are named in each class's last section, and `scripts/README.md` describes the full measurement.
