# Reproducing the plugin concordance

How to score a supported plugin against Ensembl VEP from this repository: the
data each row needs, both engines' commands and the scorer.

vep-rs matches Ensembl VEP exactly except for the documented divergences. For a
plugin, a release measures that on the plugin's populated fields: one tuple of
Location, Allele, Feature, field and value per field of an output row. The
counts of both engines' tuples and of the shared ones give a raw F1. The
adjusted F1 sets aside the tuples a documented class of
[`docs/intended-divergences.md`](../intended-divergences.md) accounts for.
`scripts/concordance/reproduce_plugin_concordance.py` runs both engines and the
scorer for one row and prints those counts.

## Choose a row and build

1. **Find the row.** `scripts/parity/plugin_gate_rows.py --rows rows.tsv`
   derives, from the support table of [Plugins](../plugins.md), what every
   supported plugin is scored on. There are three kinds of rows:

   | `scored_on` | Plugins | Scored on |
   | --- | --- | --- |
   | `clinvar` | supported, verified on real data | the genome-wide ClinVar dataset of each assembly |
   | `synthetic_fixture` | supported on the user's licensed file | the committed corpora of the plugin, whose data files are synthetic |
   | `committed_table` | TranscriptAnnotator | its committed corpora, whose table is authored for them |

   `scripts/concordance/plugin_gate_runs.toml` lists the runs behind the rows.
   A genome-wide row stands on one of five runs, named `plugins_grch37`,
   `plugins_grch38`, `gwas`, `loftee_grch38` and `loftee_grch37`. A corpus row
   stands on one
   run per set of its corpus: `base`, the manifest's own run, and each of the
   manifest's flag sets. In that file every run states the words vep-rs takes,
   the words Ensembl VEP took for the reference and the scorer's words. The
   driver lists each run with the rows it scores:

   ```bash
   python3 scripts/concordance/reproduce_plugin_concordance.py --list
   ```

   `--row` takes a row's `label` and `--assembly` its assembly. The label is
   the plugin's name, except `CADD (structural variants)` for the second CADD
   corpus and `LoFTEE (synthetic conservation data)` for LoFTEE on GRCh37. That
   row's variants are the dataset's real records, and its conservation data
   are synthetic.

2. **Build the two programs.** The driver runs a `vep` binary and the scorer,
   the `score` command of `tools/vep-plugin-pairs`, whose release build it
   finds in the tree. It needs Python 3.11 or later.

   ```bash
   cargo build --release --locked
   cargo build --release --locked \
     --manifest-path tools/vep-plugin-pairs/Cargo.toml
   ```

## Score a corpus row

3. **Run a corpus row.** A `synthetic_fixture` row and a `committed_table` row
   need nothing outside the repository. The corpus carries its input, its
   pruned cache, its data files and Ensembl VEP's output of every set. Print
   the commands first, then run them:

   ```bash
   python3 scripts/concordance/reproduce_plugin_concordance.py \
     --row BayesDel --assembly GRCh38 --work-dir work --dry-run
   python3 scripts/concordance/reproduce_plugin_concordance.py \
     --row BayesDel --assembly GRCh38 --work-dir work \
     --vep target/release/vep
   ```

   The driver runs `vep` and the scorer once per set. A set that loads no
   plugin, as the `base` set of the BayesDel corpus, populates no field of the
   plugin: it prints no line and is named on standard error.

   The plugins whose status in the support table is
   `supported on the user's licensed file, verified on synthetic fixtures`
   are reproduced on these fixtures: AVI, BayesDel, CADD, ClinPred, dbNSFP,
   dbscSNV, MaxEntScan, PrimateAI, SpliceAI and SpliceVault. Their fixture
   files hold fabricated values in the layout each plugin's module reads, and
   each corpus's `provenance.json` records how they were written.
   TranscriptAnnotator reads a table of the user's own, so its rows are scored
   on the tables committed with its two corpora
   (`--row TranscriptAnnotator` with either assembly).

   MaxEntScan's corpora read a synthetic matrix directory the tree does not
   commit. The `generated` block of the corpus manifest names the command that
   writes it, into `fordownload_synthetic` under a directory that `--generated` names to
   the driver:

   ```bash
   python3 scripts/fixtures/maxentscan_fixture_matrices.py \
     --seed 1 --out work/generated/fordownload_synthetic
   python3 scripts/concordance/reproduce_plugin_concordance.py \
     --row MaxEntScan --assembly GRCh38 --work-dir work \
     --vep target/release/vep --generated work/generated
   ```

## Stage the data for a genome-wide row

4. **Stage the data of a genome-wide row.** A `clinvar` row needs the ClinVar
   input of its assembly, the cache, the reference FASTA and the plugin's data
   files. Steps 1 and 3 of the runbook
   [Reproducing the published concordance](../../scripts/README.md#reproducing-the-published-concordance)
   give the input (`clinvar_grch37` or `clinvar_grch38`, canonical contig
   names), the FASTA and the JSON cache. `$DATA` in the commands below is the
   tree that holds this data.

   The plugin data directory holds each file under a scope directory,
   `grch37`, `grch38` or `assembly_agnostic`, at
   the path of the table below, the one `scripts/data/plugin_data_layout.sh`
   names for the file.

   The support table names each source and its terms.
   `scripts/data/setup_plugin_data.sh` prints the download and the
   preparation of each file
   ([Preparing data files](data-files.md#the-setup_plugin_datash-script)). It
   stages one assembly under one directory, so give it the scope directory,
   and place the files of the `assembly_agnostic` rows under that scope:

   ```bash
   scripts/data/setup_plugin_data.sh --assembly GRCh38 \
     --source upstream --data-dir $DATA/plugin_data/grch38
   ```

   | Plugin | Scope | File under the scope |
   | --- | --- | --- |
   | AlphaMissense | assembly | `alphamissense/alphamissense.tsv.gz`, `.tbi` |
   | AncestralAllele | assembly | `ancestral_allele/homo_sapiens_ancestor.fa.gz`, `.fai`, `.gzi` |
   | Conservation | `grch38` | `loftee/grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw` |
   | DosageSensitivity | `assembly_agnostic` | `dosagesensitivity/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz` |
   | gnomADc | assembly | `gnomadc/gnomad_coverage.tsv.gz`, `.tbi` |
   | gnomADMt | `assembly_agnostic` | `gnomadmt/gnomad.genomes.v3.1.sites.chrM.vcf.bgz`, `.tbi` |
   | GWAS | `assembly_agnostic` | `gwas/gwas_catalog_associations.tsv` |
   | LOEUF | assembly | `loeuf/loeuf_dataset.tsv.gz` on GRCh37, `loeuf/loeuf_dataset_grch38.tsv.gz` on GRCh38, `.tbi` |
   | LoFTEE | `grch37` | under `loftee_synthetic_conservation/`, synthetic conservation data and the ancestral FASTA: `GERP_scores.synthetic.txt.gz`, `.tbi`; `phylocsf_gerp.synthetic.sql`; `human_ancestor.fa.gz`, `.fai`, `.gzi` |
   | LoFTEE | `grch38` | `loftee/grch38/human_ancestor.fa`, `.fai`; the GERP bigWig Conservation reads |
   | LoFtool | `assembly_agnostic` | `loftool/LoFtool_scores.txt` |
   | OpenTargets | `grch38` | `opentargets/open_targets_vep_full.tsv.gz`, `.tbi` |
   | pLI | `assembly_agnostic` | `pli/pLI_values.txt` |
   | REVEL | assembly | `revel/revel.tsv.gz`, `.tbi` |
   | UTRAnnotator | assembly | `utrannotator/uORF_5UTR_GRCh37_PUBLIC.txt` or `uORF_5UTR_GRCh38_PUBLIC.txt` |

   Blosum62, Downstream, HGVSIntronOffset, NearestExonJB, NMD and TSSDistance
   read no data file. The driver's first action is a check, by name, that every
   file the run names stands in the directory with its index. A missing file
   ends the run before either engine starts.

   The `grch37` row of LoFTEE is the run `loftee_grch37`. Its variants are the
   dataset's real records and its conservation data are synthetic: a made-up
   GERP table and a made-up conservation database. No source publishes a
   made-up file, so `scripts/data/setup_plugin_data.sh` prints no download for
   one and reports the run's six files in a summary line of their own.
   `tools/vep-loftee-synthetic-conservation` writes the table, its index and
   the database from the GRCh37 JSON cache, every value a function of the seed
   and the row's key
   ([what the files hold](loftee-grch37.md#how-grch37-is-compared-with-ensembl-vep)):

   ```bash
   TOOL=tools/vep-loftee-synthetic-conservation
   RUN=$DATA/plugin_data/grch37/loftee_synthetic_conservation
   cargo build --release --locked --manifest-path $TOOL/Cargo.toml
   $TOOL/target/release/vep-loftee-synthetic-conservation \
     --cache $DATA/caches/vep-rs/grch37 --seed 1 --out-dir $RUN
   ```

   The run's ancestral FASTA is the AncestralAllele file of the table above,
   written under the bare sequence names that the plugin's `samtools` read
   asks for (the `master` branch's `LoF.pm` lines 552-553). It is block-gzipped
   and indexed with `samtools faidx`:

   ```bash
   ANC=$DATA/plugin_data/grch37/ancestral_allele/homo_sapiens_ancestor.fa.gz
   bgzip -dc $ANC \
     | sed -E 's/^>ANCESTOR_for_[a-z]+:GRCh37:([^:]+):.*$/>\1/' \
     | bgzip > $RUN/human_ancestor.fa.gz
   samtools faidx $RUN/human_ancestor.fa.gz
   ```

## Prepare the Ensembl VEP reference

5. **Write LoFTEE's stand-in files.** LoFTEE's module requires scripts from
   its `maxEntScan/` directory, which holds MaxEntScan's files, and hands
   paths under it to their subroutines. The reference runs with the module's
   splice predictions off, where none of those scores is computed, so the
   directory is left out of the checkout and stand-in files take its place.
   The tag below is on LoFTEE's `grch38` branch:

   ```bash
   git clone --filter=blob:none --no-checkout \
     https://github.com/konradjk/loftee loftee
   git -C loftee sparse-checkout set --no-cone '/*' '!/maxEntScan/'
   git -C loftee checkout v1.0.4_GRCh38
   python3 scripts/data/write_loftee_stand_ins.py loftee
   ```

   A GRCh37 reference runs the plugin's `master` branch. Its checkout is made
   the same way, at the commit that the `provenance.json` of the GRCh37 LoFTEE
   corpus records:

   ```bash
   git clone --filter=blob:none --no-checkout \
     https://github.com/konradjk/loftee loftee-master
   git -C loftee-master sparse-checkout set --no-cone '/*' '!/maxEntScan/'
   git -C loftee-master checkout c7a87dffe3b96a729d299ca3ccb1c25db7b44e30
   python3 scripts/data/write_loftee_stand_ins.py loftee-master
   ```

   The writer reads the module's own Perl files for the scripts they `require`
   from that directory and the subroutines they call that no file of the module
   defines; `--dry-run` prints both. It writes one stand-in per script.

   - A subroutine the module's constructor calls returns an empty list.
   - A subroutine reached only with splice predictions on ends the run with
     a message, so no row carries a made-up score.

   No MaxEntScan file is read at any point. The writer refuses a
   directory that holds a file it did not write. In that checkout `git status`
   lists the stand-ins as modified files, and a `git diff` would fetch the
   upstream files of the same names.

   LoFTEE also calls `samtools` for its ancestral-allele filter, which the
   Ensembl VEP image lacks. `--tool-dir` names a directory holding a `samtools`
   binary that runs inside the image, and the driver puts that directory first
   on the container's `PATH` ([LoFTEE](loftee.md#4-the-ancestral-allele-filter-withheld-from-every-allele-of-a-multi-allelic-record)).

6. **Get Ensembl VEP's output.** The reference of a run is Ensembl VEP 116.3's
   default-format output of the same input under the run's `ensembl_vep`
   words. `--reference` names a stored output.

   `--run-ensembl-vep` writes one with Docker, from the image
   `scripts/concordance/generate_reference_output.sh` pins,
   `ensemblorg/ensembl-vep:release_116.3` at the digest
   `sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3`.
   It needs the Ensembl VEP cache of release 116, whose download that script's
   header describes. The output lands at
   `<work-dir>/<run>/ensembl-vep/output.txt`. One output serves every later
   scoring of the same input, so the output of a genome-wide run is worth
   keeping.

   A corpus set runs through that script on a copy of the corpus in the work
   directory. A corpus whose plugin module the image does not carry also takes
   `--plugin-modules`, the directory of Ensembl VEP's plugin modules, a
   `VEP_plugins` checkout at `release/116`; the driver's error line asks for
   it. MaxEntScan's sets are scored against their committed output alone,
   since that script cannot name their generated directory.

## Score a genome-wide row and a whole-genome cell

7. **Run a genome-wide row.** Print the commands, then run them. `$DATA` is
   the data tree of step 4:

   ```bash
   VCF=$DATA/vcf/grch38/all_variants_canonical
   python3 scripts/concordance/reproduce_plugin_concordance.py \
     --row GWAS --assembly GRCh38 --work-dir work \
     --vep target/release/vep --fork 16 \
     --input $VCF/clinvar_grch38.canonical.vcf.gz \
     --vep-rs-cache $DATA/caches/vep-rs/grch38 \
     --fasta $DATA/reference/grch38/genome.fa \
     --plugin-data $DATA/plugin_data \
     --ensembl-vep-cache $DATA/caches/perl/vep-cache \
     --run-ensembl-vep \
     --pairing-tool tools/vep-plugin-pairs/target/release/vep-plugin-pairs
   ```

   Add `--dry-run` to print the three commands without running them: the
   container, `vep` and the scorer. `--steps` runs a part where the outputs of
   the rest are held, for instance `--steps check,score` on an engine output
   already written. The LoFTEE row adds `--loftee loftee` and `--tool-dir`
   with the `samtools` directory of step 5. The GRCh37 row,
   `LoFTEE (synthetic conservation data)`, adds the same two options with the
   `master` checkout, `--loftee loftee-master`.

   The rows of `plugins_grch37` and `plugins_grch38` run together, all plugins
   of the run loaded at once, and add `--utrannotator-sweep`. That directory
   holds five files the run's `score` words name, `sweep-seed<seed>.txt`. Each
   is Ensembl VEP's default-format output of the input's records with
   `--hgvs` and the UTRAnnotator plugin, written with `PERL_HASH_SEED` set to
   the seed of its name and `PERL_PERTURB_KEYS=0`. The class
   [`utrannotator_consequence_precedence`](../divergences/utrannotator-consequence-precedence.md)
   reads them.

   The same `score` words name the run's REVEL table after `--revel-table`,
   which the class
   [`revel_other_reference_amino_acid_rows`](../divergences/revel-other-reference-amino-acid-rows.md)
   reads.

   The `score` words of a SpliceAI run name its corpus's two synthetic score
   files, each after `--spliceai-table`, which the class
   [`spliceai_multi_entry_cell`](spliceai.md#two-genes-in-one-record) reads.

   The step `utrannotator-sweep` writes the five files. It runs only where
   `--steps` names it:

   ```bash
   python3 scripts/concordance/reproduce_plugin_concordance.py \
     --run plugins_grch38 --work-dir work --steps utrannotator-sweep \
     --input $VCF/clinvar_grch38.canonical.vcf.gz \
     --fasta $DATA/reference/grch38/genome.fa \
     --plugin-data $DATA/plugin_data \
     --ensembl-vep-cache $DATA/caches/perl/vep-cache \
     --utrannotator-sweep work/sweep
   ```

   Each seed runs in a container with no network, in one process, with
   UTRAnnotator alone of the plugins. The runs read the run's `--input` by
   default, and `--sweep-input` names another file of records. For a release's
   own scores the seeded runs read the records behind the rows that a first
   score, made without the five files, leaves as UTRAnnotator residuals. A file
   that stands is written over only under `--overwrite-sweep`.

8. **Run a whole-genome cell.** The same driver scores a plugin group on a
   population dataset: name the run, give one chromosome's file as its input
   and name the dataset for the printed lines. The engines' words, the data
   files and the scorer's words are the run's own. The reference is Ensembl VEP's output
   of that file under the run's words, stored and named with `--reference` as
   below, or written with `--run-ensembl-vep` as in step 7. The sweep directory
   holds the five seeded outputs for that file's records.

   ```bash
   VCF=$DATA/vcf/grch38/all_variants_canonical
   python3 scripts/concordance/reproduce_plugin_concordance.py \
     --run plugins_grch38 --work-dir work --dry-run \
     --vep target/release/vep --fork 16 \
     --input $VCF/gnomad_genomes_v4.1_chr21.canonical.vcf.gz \
     --dataset gnomad_genomes_v4.1_chr21 \
     --reference work/chr21.ensembl-vep.txt \
     --vep-rs-cache $DATA/caches/vep-rs/grch38 \
     --fasta $DATA/reference/grch38/genome.fa \
     --plugin-data $DATA/plugin_data \
     --utrannotator-sweep work/sweep \
     --pairing-tool tools/vep-plugin-pairs/target/release/vep-plugin-pairs
   ```

## Read the scores

9. **Read the result.** Standard output carries one JSON line per plugin and
   projection, in the form of an entry of a release's plugin record. The
   `feature` projection counts the tuples as defined above, and the `allele`
   projection counts them without the Feature.

   A line's keys name the row:

   - `plugin`
   - `assembly`
   - `dataset` and `family` or `corpus` and `expected_set`
   - `engine`
   - `projection`
   - `scored_on`

   Its counts are:

   - the raw `perl`, `rust` and `intersection`
   - the adjusted `adj_perl`, `adj_rust` and `adj_intersection`
   - `masked` by class
   - `residual_rows`
   - `text_only_differences`

   Each line adds `raw_f1` and `adjusted_f1`, computed from those counts to
   nine decimal places. Standard
   error carries one line per row, summed over the sets of a corpus. The
   scorer's own report and its residual files stay in the work directory.

   A release states its figures in a plugin record under
   `docs/concordance-provenance/`, a file named
   `<date>-plugins-v<version>.json` whose `record_type` is `plugins`. Its lists
   `real_data`, `population_data`, `synthetic_fixtures` and `committed_table`
   hold one entry per plugin, assembly, dataset or corpus set and projection,
   in the form above, and its `rows` list repeats the derived rows. A record entry also carries
   `run_manifest_sha256`, the digest of the release's own run manifest, which
   the driver's lines leave out. The expected figures of a row are that
   record's entries for it.

What differs from one run to the next, and what does not:

- Ensembl VEP's NearestExonJB and UTRAnnotator plugins print the items of a
  value in Perl's hash order, which changes between processes. The classes
  [`nearest_exon_jb_item_order`](../divergences/nearest-exon-jb-item-order.md)
  and
  [`utrannotator_annotation_pair_order`](../divergences/utrannotator-annotation-pair-order.md)
  compare such a value as a set of items. `--perl-env PERL_HASH_SEED=0` with
  `--perl-env PERL_PERTURB_KEYS=0` fixes the order for a reference written
  with `--run-ensembl-vep`.
- The header lines of an output state the time of the run and the paths it
  read, and the statistics report states its run time. Wall time depends on
  the machine and enters no figure.
- Against one reference output, every tuple count and both F1 values are the
  same on every run of the same binary, cache, FASTA, input and data files.

## Where next

- [Plugins](../plugins.md): the support table, with each plugin's source, terms
  and status
- [Preparing data files](data-files.md#the-setup_plugin_datash-script): the
  script that prints the download and the preparation of each file
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the
  documented classes the adjusted F1 sets aside
- [Reproducing the published concordance](../../scripts/README.md#reproducing-the-published-concordance):
  the runbook that gives the input, the FASTA and the JSON cache
