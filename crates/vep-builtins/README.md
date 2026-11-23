# vep-builtins

Built-in VEP plugin implementations compiled directly into the `vep` binary, and
the `--custom` annotation sources that run beside them
([Custom annotation sources](#custom-annotation-sources)).

## Architecture

Instead of compiling each plugin as a separate `cdylib` shared library (the C FFI
approach in `vep-plugin`), built-in plugins operate on native Rust types with no
JSON serialization or FFI overhead.

```
vep-cli (runner.rs)
  |
  +-- vep-builtins (this crate)
  |     |
  |     +-- BuiltinPlugin trait
  |     |     +-- feature_types()          <-- the row kinds the plugin writes on (Ensembl VEP's feature_types)
  |     |     +-- annotates_structural_variants()  <-- false unless the Ensembl VEP module declares the kind
  |     |     +-- prefetch() + annotate()  <-- two-phase: parallel I/O, then sequential writes
  |     |     +-- run_batch()              <-- single-phase fallback, once per buffer
  |     |     +-- run()                    <-- per-transcript-consequence fallback
  |     |     +-- run_naming_empty_strings()  <-- run() with the fields whose empty value is the empty string
  |     |
  |     +-- AnnotationStore trait (one query interface over two backends)
  |     |     +-- TabixAnnotator: noodles::tabix + CSI binning index over `.tsv.gz`
  |     |     +-- BinaryAnnotator: mmap-backed `.vpd` data + `.vpdi` index, chosen
  |     |     |   by `open_backend` when a `.vpd` sits beside the `.tsv.gz`
  |     |     |   (`vep-cache-converter convert-plugin-data` writes it) and its
  |     |     |   index records the columns the plugin reads and the SHA256 and
  |     |     |   length the table and the `.vpd` have, and ends with the SHA256
  |     |     |   of its own bytes; a store refused on any of them leaves the
  |     |     |   table read through its `.tbi`, with a warning
  |     |     +-- Batch region merging (1000 bp gap)
  |     |     +-- Allele matching & normalization
  |     |
  |     +-- BuiltinRegistry
  |     |     +-- Name-based plugin resolution
  |     |     +-- Returns the names no built-in claims to the caller
  |     |
  |     +-- custom (the --custom sources, run on each batch after the plugins)
  |           +-- grammar: one CustomSpec per --custom string
  |           +-- source: CustomSources, opened before the input is read
  |           +-- vcf, bed, gff_overlay (tabix-indexed), bigwig
  |
  +-- vep-plugin (C FFI dylib loader, for the names no built-in claims)
```

## Plugins

Support is per reference genome, as in the
[plugin support table](../../docs/plugins.md#support-table). The Mechanism
column names how the plugin reads its data and matches a variant.

| Plugin        | Mechanism                                                                                      | GRCh37 | GRCh38   | Output fields                                                                                                                                 |
| ------------- | ---------------------------------------------------------------------------------------------- | ------ | -------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| CADD          | Annotation store, position + allele match; a structural variant by class and start (`sv=`)     | Yes    | Yes      | `CADD_PHRED`, `CADD_RAW`                                                                                                                      |
| REVEL         | Annotation store, position + allele match, the record's `aaref` and `aaalt` equal to the row's reference and alternate amino acids, transcript match unless `no_match`; missense only | Yes    | Yes      | `REVEL`                                                                                                                                       |
| gnomADc       | Annotation store, position only; the mean over the record's span                               | Yes    | Yes      | `gnomAD_<key>_cov` (`mean`, `median` or `median_approx` and `total_DP`, `1x` to `100x`; the keys chosen by the file's column count)            |
| AlphaMissense | Annotation store, position + allele match, the record's `protein_variant` letters equal to the row's amino acids, transcript match under `transcript_match=1`; missense only | Yes    | Yes      | `am_pathogenicity`, `am_class` by default; any column under `cols=`                                                                           |
| dbscSNV       | Annotation store, position + alternate allele of a single-base substitution; the variant's scores on every transcript row | Yes    | Yes      | `ada_score`, `rf_score`                                                                                                                       |
| LoFTEE        | Consequence filter over the transcript model, `--fasta` splice motifs, GERP                    | Yes    | Yes      | `LoF`, `LoF_filter`, `LoF_flags`, `LoF_info`                                                                                                  |
| LoFtool       | Flat file loaded at init, gene symbol lookup                                                   | Yes    | Yes      | `LoFtool`                                                                                                                                     |
| pLI           | Flat file loaded at init, gene symbol lookup, or transcript stable id lookup in `transcript` mode | Yes | Yes | `pLI_gene_value`; `pLI_transcript_value` in `transcript` mode                                                                                 |
| DosageSensitivity | Flat file loaded at init, gene symbol lookup on transcript rows of sequence and structural variants; under `cover=1` only a row whose variant spans the whole transcript | Yes | Yes | `pHaplo`, `pTriplo`                                                                                                                           |
| MaxEntScan    | Matrix directory loaded at init, `--fasta` windows at the transcript's splice sites and around the variant | Yes | Yes  | `MaxEntScan_ref`, `MaxEntScan_alt`, `MaxEntScan_diff`; with `SWA` `MES-SWA_donor_*` and `MES-SWA_acceptor_*`; with `NCSS` `MES-NCSS_*`; `verbose` adds the sequences |
| GWAS          | Flat file loaded at init, position + risk allele match                                         | No     | Yes      | `GWAS_accessions`, `GWAS_associated_gene`, `GWAS_beta_coef`, `GWAS_odds_ratio`, `GWAS_p_value`, `GWAS_pmid`, `GWAS_risk_allele`, `GWAS_study` |
| SpliceAI      | Annotation store, position + allele + gene symbol match                                        | SNV and indel files | SNV and indel files | `SpliceAI_pred` (default and tab); `SpliceAI_pred_*` (nine fields, VCF or `split_output`); `SpliceAI_cutoff` under `cutoff` |
| AVI           | Annotation store, position + allele match; single-base substitutions only                      | No     | Yes      | `AVI_PHRED`; `AVI_RAW` under `raw=1`; `AVI_LINK` under `links=1`                                                                              |
| AncestralAllele | Indexed FASTA of ancestral sequences, the record's span read as written                      | Yes    | Yes      | `AA`                                                                                                                                          |
| dbNSFP        | Annotation store, position + alternate allele, the amino acid change by default (`pep_match`), the transcript's entry of a column split on `;` under `transcript_match=1`; on rows of the default consequence set or `consequence=` | Yes    | Yes      | The selected columns under their own names (`dbNSFP_` prefixed under `output_prefix=1`), or every column with `ALL`; the user's licensed file |
| BayesDel      | Annotation store, position + reference and alternate allele match; the record's score on every row, intergenic rows included | Yes    | Yes      | `BayesDel`                                                                                                                |
| ClinPred      | Annotation store, position + alternate allele of a single-base substitution; missense, stop-gained, stop-lost and start-lost rows | Yes    | Yes      | `ClinPred`                                                                                                                |
| PrimateAI     | Annotation store, position + alternate allele of a single-base substitution (the reference column is not read); the score on every row, intergenic rows included; a file name carrying an assembly must match the run | Yes    | Yes      | `PrimateAI`                                                                                                               |
| OpenTargets   | Annotation store, position + allele match; one record per association, the selected columns `,`-joined across them (one object per record in JSON) | No     | Yes      | `OpenTargets_<column>`: `gwasLocusToGeneScore`, `gwasGeneId`, `gwasDiseases`, `qtlGeneId` and `qtlBiosampleName` by default; any header column under `cols=` |
| SpliceVault   | Annotation store, position + allele match; the row naming the transcript under `transcript_match=1` (the default), the first allele match under `transcript_match=0` | Yes    | Yes      | `SpliceVault_top_events`, `SpliceVault_site_sample_count`, `SpliceVault_site_max_depth`, `SpliceVault_out_of_frame_events`, `SpliceVault_site_pos`, `SpliceVault_site_type`, `SpliceVault_SpliceAI_delta` |
| TranscriptAnnotator | Annotation store, position + allele match + the transcript's stable id or a RefSeq cross-reference; the user's own table | Yes    | Yes      | `<prefix>_<column>` per selected column (`cols=`), the prefix the file's basename unless `prefix=` names another       |
| gnomADMt      | Annotation store over the gnomAD chrM sites VCF, position + allele match; the record's INFO fields on every row, intergenic rows included | No     | Yes      | `gnomAD_<field>`: `hap_AC_hom`, `hap_AC_het`, `hap_AF_hom`, `hap_AF_het` and `hap_AN` by default; any INFO field under `fields=` |
| LOEUF         | Annotation store over gene spans; the record of the row's transcript under `match_by=transcript`, the gene's minimum under `match_by=gene`; transcript rows other than upstream and downstream ones | Yes    | Yes      | `LOEUF`                                                                                                                   |
| Conservation  | bigWig of per-base GERP scores, the record's span (an insertion's two flanks) averaged, or its maximum under `MAX` | No     | Yes      | `Conservation`                                                                                                            |
| Blosum62      | Embedded matrix, the row's amino-acid change                                                   | Yes    | Yes      | `BLOSUM62`                                                                                                                                    |
| Downstream    | Transcript model: the changed coding sequence translated from the frameshift                   | Yes    | Yes      | `DownstreamProtein`, `ProteinLengthChange`                                                                                                    |
| HGVSIntronOffset | The HGVSc coordinates' intron offsets, under `--hgvs`                                       | Yes    | Yes      | `HGVS_IntronStartOffset`, `HGVS_IntronEndOffset`                                                                                              |
| NearestExonJB | Transcript model: the overlapped exons' junction distances                                     | Yes    | Yes      | `NearestExonJB`                                                                                                                               |
| TSSDistance   | Transcript span: the distance to the transcription start site                                  | Yes    | Yes      | `TSSDistance`                                                                                                                                 |
| NMD           | Transcript model: a truncating variant's position in the exon model against the four escape rules of nonsense-mediated decay | Yes    | Yes      | `NMD`                                                                                                                                         |
| UTRAnnotator  | Transcript model and `--fasta`: the 5' UTR with the allele spliced in, scanned for upstream open reading frame (uORF) effects; the uORF evidence table under `file=` | Yes    | Yes      | `5UTR_consequence`, `5UTR_annotation`, `Existing_uORFs`, `Existing_OutOfFrame_oORFs`, `Existing_InFrame_oORFs`                             |

LoFtool and pLI are keyed by gene symbol, so one data file serves both assemblies;
given no file, each reads the copy bundled beside Ensembl VEP's plugin modules
from the `--dir_plugins` directory.
LoFTEE reads GERP from a tabix TSV on GRCh37 (`gerp_tabix=`) and from a bigWig on
GRCh38 (`gerp_bigwig=`), the form each branch of Ensembl VEP's plugin reads
(`LoF.pm` lines 138-141 of `master`, line 130 of `grch38`), or, when no GERP
word is given, from the path Ensembl VEP's plugin derives (`loftee_path` followed
by the GRCh38 bigWig's upstream name; a fixed path on GRCh37). A GRCh37 run also
passes `run_splice_predictions=0`: Ensembl VEP's plugin for that assembly runs its
splice-prediction layer by default, which is not ported, and a run that leaves it
on is refused. Under `conservation_file=` LoFTEE reads the SQLite database of
Ensembl VEP's plugin: the PhyloCSF items on both assemblies and, on GRCh37, the
GERP total of each whole exon. Conservation reads any bigWig of per-base scores; Ensembl Compara
publishes one for GRCh38 alone, so no GRCh37 file is staged for it.

A built-in plugin writes on the row kinds its `feature_types()` names and, through
`annotates_structural_variants()`, on a structural variant's rows only when its
Ensembl VEP module lists `StructuralVariationFeature` or `BaseVariationFeature`
in `variant_feature_types`. CADD, DosageSensitivity, TSSDistance and UTRAnnotator
do; every other built-in skips those rows, as Ensembl VEP does
([Rows a plugin writes on](../../docs/plugins.md#rows-a-plugin-writes-on)). The
default `run_batch()` applies both gates, and a plugin that overrides
`prefetch()`, `annotate()` or `run_batch()` skips a structural variant itself,
through `annotates_variant()` or its own test of the record, unless it
annotates both kinds.

## Batch Query Strategy

`TabixAnnotator` issues one tabix query per merged region rather than one per
variant:

```text
per variant:
  tabix query(chr, pos)           # one index lookup per variant

per buffer (TabixAnnotator):
  sort the buffer's regions by start within each chromosome
  merge overlapping regions (bridging gaps of up to 1000 bp)
  tabix query per merged region   # one index lookup per merged region
  index the records by position
  for each variant:
    in-memory range lookup        # no I/O
```

`BinaryAnnotator` answers the same merged-region queries by binary search over a
memory-mapped `.vpd` file, with no decompression.

## Usage

Plugins are loaded via the standard `--plugin` CLI flag:

```bash
vep -i input.vcf \
    --assembly GRCh38 \
    --plugin CADD,snv=/path/to/cadd_snv.tsv.gz,indels=/path/to/cadd_indels.tsv.gz \
    --plugin gnomADc,/path/to/gnomad_coverage.tsv.gz \
    --plugin REVEL,file=/path/to/revel.tsv.gz \
    --offline --cache
```

Pass `--assembly` whenever plugins are enabled; without the flag the cache's own
assembly serves. CADD and AlphaMissense use it to reject a data file whose own
build marker contradicts the run; REVEL, dbNSFP and dbscSNV use it to pick the
coordinate column they read and refuse the file indexed
for the other assembly; LoFTEE uses it to select the END_TRUNC cutoff and the
default GERP source. When neither the flag nor the cache names one, REVEL refuses
to load, dbscSNV too unless its second parameter names the file's assembly, dbNSFP
reads the GRCh38 column, LoFTEE keeps the GRCh38 cutoff and default, and the file
checks of CADD and AlphaMissense are skipped. gnomADc carries no in-file build
marker and is not checked; AVI checks none either: a GRCh37 run is warned,
not refused.

Under `--refseq` or `--merged`, five plugins that match rows to Ensembl
transcript ids warn that the match may fail: REVEL unless `no_match`,
AlphaMissense and dbNSFP under `transcript_match=1`, SpliceVault unless
`transcript_match=0`, and pLI in `transcript` mode. The warning is printed on
standard error alone, under `--quiet` too, and never in the warnings file, as
Ensembl VEP prints a plugin's bare `warn` (`traits::print_on_stderr` and
`warn_on_stderr`); dbNSFP's column warning, its file-name version warning and
its two per-row `ERROR:` lines for a malformed transcript list take the same
channel. Every other diagnostic of the crate is a `tracing` event.

The runner resolves `--plugin Name,...` arguments against the built-in registry
first. A name no built-in claims is loaded as a C-ABI dynamic library through
`vep-plugin`.

## Custom annotation sources

The `custom` module reads the `--custom` sources: a tabix-indexed VCF, BED or
GFF (General Feature Format) file, or a bigWig, whose records annotate the input
variants under a short name
([Custom annotation](../../docs/custom-annotation.md)). `custom::grammar`
parses each `--custom` string into a `CustomSpec` as Ensembl VEP reads it and
refuses a malformed string with Ensembl VEP's own text, a remote file, and a GFF
or GTF (Gene Transfer Format) transcript model, which vep-rs annotates against
as a cache built from the file. `custom::source` opens every source before the
input is read (`CustomSources::open`) and runs the set on each batch after the
plugins (`CustomSources::run_batch`): a `prefetch` per source in parallel, then
a sequential `annotate`, under the overlap predicates, the record list and the
summary statistics every format shares. `custom::vcf`, `custom::bed` and
`custom::gff_overlay` read the tabix-indexed formats, the last one the
`gff_type=gencode_promoter` promoter windows, while a GFF in a `gff_type` Ensembl
VEP does not know opens and annotates nothing; `custom::bigwig` reads a bigWig.
Each source leaves its records and statistics on the variant as
`vep_core::custom_annotation::CustomAnnotations`, which the `vep-io` field table
and JSON formatter print.

A `COSV` record named like one an earlier source listed under the same short
name lends that record the fields it lacks and is not listed again. The lending
source counts the match, so a source can match a record and list none of its own
([Records and statistics per input
record](../../docs/custom-annotation.md#records-and-statistics-per-input-record)).

Each `--custom` source of the `custom` module lists the chromosomes its file
names (`CustomSource::chromosomes`: a tabix index's sequence names, or a bigWig's
chromosomes with a leading `chr` removed), read once when the source opens. The
run admits a record on a chromosome the cache lacks when a source lists it
([Chromosomes](../../docs/custom-annotation.md#chromosomes)).
