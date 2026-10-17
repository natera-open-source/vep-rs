# CLI Reference

The `vep` binary accepts the Ensembl VEP command line: every flag is Ensembl VEP's, does what Ensembl VEP's does, or is refused with an error naming what it would do, and the few flags vep-rs adds are marked as its own. A flag the error table below governs (refused, or accepted without effect) carries that table's text as its help.

This page is generated from the command-line definition and the option tables by the `cli_reference` test of `crates/vep-cli`; edit the sources, then run `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference` to rewrite it.

```bash
# Show all flags
vep --help

# Short form
vep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache
```

Worked examples are in [examples.md](examples.md); the output formats and their parity with Ensembl VEP in [output-formats.md](output-formats.md).

## Flags

Grouped by section as `--help` groups them. A repeated single-value flag keeps its last value, as Getopt::Long does; `-o` accumulates and the last wins with a warning, as the `vep` script does (vep:32-33). Long option names are matched exactly as written.

### Input

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--input_file`, `-i` | String | `STDIN` | Input file path (use "STDIN" for standard input) |
| `--input_data` | String | -- | Input data string (accepted for VEP compatibility; not implemented, input is read from --input_file) |
| `--format` | `ensembl`, `vcf`, `hgvs`, `id`, `spdi`, `region`, `guess` | `guess` | Input file format ("guess" = auto-detect) |
| `--gp` | bool | false | Refused with an error: --gp is not supported: it reads 1000 Genomes pilot coordinates from INFO/GP, which no current VCF carries |
| `--ambiguous_hgvs` | bool | false | Refused with an error: --ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g.) only |
| `--phased` | bool | false | Accepted without effect: Ensembl VEP marks the genotypes as phased and nothing reads the mark |
| `--no_check_variants_order` | bool | false | Accepted without effect: vep-rs never checks the input order: every record is annotated where it stands, and an input whose chromosomes interleave keeps more transcript indexes in memory |

### Output

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--output_file`, `-o` | String (repeatable) | `variant_effect_output.txt` | Output file path (use "STDOUT" for standard output, in any case); the last of several wins |
| `--output_format` | `vep`, `vcf`, `json`, `tab`, `parquet` | `vep` | Output format. `parquet` writes a TSV intermediate and post-processes it to Parquet through a `duckdb` CLI subprocess; the `duckdb` binary must be on PATH |
| `--parquet_shape` | `nested`, `flat` | `nested` | Parquet row shape (only effective when --output_format=parquet). nested (default): one row per variant, CSQ subfields are LIST<T>. flat: one row per (variant × transcript-consequence) tuple, CSQ subfields are scalar |
| `--keep_intermediate` (hidden) | bool | false | Keep the Parquet intermediate TSV after DuckDB finalize succeeds. Debug aid: normally the intermediate is deleted once the Parquet directory is written |
| `--parquet_row_group_size` (hidden) | u64 | `122880` | Rows per Parquet row group (only effective when --output_format=parquet). Rounded up to a multiple of 2,048, DuckDB's vector size, which is also the dictionary size limit. Test aid for producing several row groups from a small input |
| `--vcf` | bool | false | Shorthand for --output_format vcf; also switches on --symbol, --biotype and --numbers |
| `--json` | bool | false | Shorthand for --output_format json |
| `--tab` | bool | false | Shorthand for --output_format tab |
| `--compress_output` | `gzip`, `bgzip` | -- | Compress output file (accepted for VEP compatibility; not implemented, output is written uncompressed) |
| `--force_overwrite`, `--force` | bool | false | Overwrite output file if it already exists |
| `--no_headers` | bool | false | Suppress header lines in output |
| `--mask_header_cache_path` | bool | false | Print the cache directory in the output headers as [PATH]/<leaf> |
| `--fields` | String | -- | Accepted for VEP compatibility; not implemented (the columns follow the field flags) |
| `--vcf_info_field` | String | `CSQ` | VCF INFO field name for VEP annotation |
| `--keep_ann` | bool | false | Accepted without effect: Ensembl VEP declares it as a synonym of --keep_csq and reads nothing from it |
| `--solr` | bool | false | Accepted without effect: Ensembl VEP declares it and writes no Solr output |
| `--convert` | String | -- | Refused with an error: --convert has been deprecated |

### Cache / Offline

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--cache` | bool | false | Accepted for VEP compatibility; no effect (the annotation cache is --json_cache) |
| `--offline` | bool | false | Accepted without effect: vep-rs always runs offline and never contacts a database |
| `--dir` | String | `~/.vep` | Top-level VEP data directory (accepted for VEP compatibility; nothing is loaded from it, it only heads the cache path in the output header) |
| `--dir_cache` | String | -- | Cache directory (defaults to --dir; accepted for VEP compatibility, nothing is loaded from it, it only names the cache path in the output header) |
| `--cache_version` | u32 | -- | Cache version named in the output header's cache path (defaults to the VEP release; selects no cache) |
| `--json_cache` | String | -- | Path to the JSON transcript cache directory (built by vep-cache-builder or converted from a Perl VEP cache) |
| `--fasta`, `--fa` | String | -- | FASTA file for reference sequence lookup |
| `--shift_3prime` | bool | false | Enable 3' shifting of insertions/deletions (requires --fasta) |
| `--sereal` | bool | false | Accepted without effect: the cache's serialiser is read from its own metadata |

### Species / Assembly

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--species`, `-s` | String | `homo_sapiens` | Species named in the output header's cache path (the annotated data is the --json_cache contents) |
| `--assembly`, `-a` | String | -- | Genome assembly (e.g., GRCh38) |
| `--grch37` | bool | false | Annotate against GRCh37: sets --assembly GRCh37, over an explicit --assembly, as Ensembl VEP's option set does |

### Database

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--database` | bool | false | Refused with an error: --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with --json_cache, or --cache with --dir_cache |
| `--host` | String | -- | Refused with an error: --host is not supported: it names the Ensembl database host, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache |
| `--port` | String | -- | Refused with an error: --port is not supported: it names the Ensembl database port, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache |
| `--user`, `-u` | String | -- | Refused with an error: --user is not supported: it names the Ensembl database user, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache |
| `--password`, `--pass` | String | -- | Refused with an error: --password is not supported: it is the Ensembl database password, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache |
| `--registry` | String | -- | Refused with an error: --registry is not supported: it names an Ensembl registry file of database connections, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache |
| `--db_version` | String | -- | Refused with an error: --db_version is not supported: it selects an Ensembl database version; --cache_version selects a cache |
| `--genomes` | bool | false | Refused with an error: --genomes is not supported: it selects the Ensembl Genomes database host; vep-rs reads a local cache only |
| `--is_multispecies` | String | -- | Refused with an error: --is_multispecies is not supported: it describes a multi-species Ensembl database; vep-rs reads a local cache only |

### Runtime

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--buffer_size` | usize | `5000` | Number of variants to process per batch |
| `--max_sv_size` | i64 | `10000000` | A structural variant spanning more than this many bases is oversize; -1 lifts the limit, as in Ensembl VEP 116. A <CNV:TR> record whose repeat sequence is longer than this stays symbolic rather than reading a direction of change (under -1 every one does), as in Ensembl VEP 116. Ensembl VEP 116 drops an oversize record before annotation and writes nothing for it in any format; vep-rs annotates it: rows in the default, tab and parquet outputs, its VCF line without consequences, no JSON object (docs/output-formats.md and docs/intended-divergences.md describe the divergence) |
| `--fork`, `--threads` | usize | `0` | Number of parallel annotation worker threads (Ensembl VEP forks processes; vep-rs runs threads and is parallel by default). Default 0 = auto: uses the number of available logical CPUs, capped at 32 (wall time stops improving past that and rises on wider hosts). Set to 1 to force single-threaded annotation, or any positive value to pin the worker count; user-supplied values are not capped |
| `--decompression_threads` | usize | `0` | Number of BGZF decompression worker threads for .gz/.bgz VCF inputs. Default 0 = auto: uses the number of available logical CPUs, capped at 10 (measured sweet-spot: throughput plateaus around 4-10 threads and is flat past that). User-supplied positive values are not capped: pass any integer to pin the worker count, e.g. `--decompression_threads=1` for single-threaded |
| `--transcript_index_impl` (hidden) | `bin`, `coitrees`, `sorted` | `bin` | Transcript overlap index implementation (hidden, for A/B benchmarking). Options: bin (UCSC binning, default), coitrees (implicit interval tree), sorted (sorted array + suffix_max_end). All three must produce identical overlap sets; any difference between implementations indicates a bug |
| `--distance` | String | `5000` | Upstream/downstream distance in bp (or "up,down" for asymmetric) |

### Output Control

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--everything`, `-e` | bool | false | Shorthand for Ensembl VEP's --everything set: --sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype and --mirna, each overriding an explicit value as in Ensembl VEP. Every column of Ensembl VEP's --everything except the regulatory columns: --regulatory is switched on and annotates nothing, so no RegulatoryFeature or MotifFeature row is written |
| `--symbol` | bool | false | Add gene symbol to output |
| `--biotype` | bool | false | Add transcript biotype to output |
| `--canonical` | bool | false | Flag canonical transcripts |
| `--mane` | bool | false | Add MANE Select and MANE Plus Clinical flags |
| `--mane_select` | bool | false | Add MANE Select flag only |
| `--tsl` | bool | false | Add Transcript Support Level |
| `--appris` | bool | false | Add APPRIS isoform annotation |
| `--ccds` | bool | false | Add CCDS transcript ID |
| `--protein` | bool | false | Add Ensembl protein ID |
| `--protein_version` | bool | false | Add the translation version to the Ensembl protein ID (ENSP00000286808.3); needs --protein |
| `--uniprot` | bool | false | Add UniProt cross-references |
| `--xref_refseq` | bool | false | Add RefSeq transcript cross-references |
| `--hgvs` | bool | false | Add HGVS coding and protein notations |
| `--hgvsg` | bool | false | Add the HGVSg column (accepted for VEP compatibility; the value is not computed and stays empty) |
| `--sift` | `s`, `p`, `b` | -- | SIFT predictions (p = prediction, s = score, b = both) |
| `--polyphen` | `s`, `p`, `b` | -- | PolyPhen predictions (p = prediction, s = score, b = both) |
| `--condel` | String | -- | Accepted without effect: Ensembl VEP declares it and reads nothing from it |
| `--numbers` | bool | false | Add exon/intron numbers |
| `--domains` | bool | false | Add overlapping protein domain information |
| `--regulatory` | bool | false | Accepted for VEP compatibility; regulatory-feature annotation is not implemented and the flag warns |
| `--regulatory_gff` | String | -- | Refused with an error: --regulatory_gff is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache, and regulatory annotation is not supported |
| `--extended_promoters` | bool | false | Refused with an error: --extended_promoters is not supported: it would widen the promoters read with --regulatory_gff to their extended bounds, and regulatory annotation is not supported |
| `--variant_class` | bool | false | Add variant class from SO |
| `--check_existing` | bool | false | Check for co-located known variants |
| `--coding_only` | bool | false | Accepted for VEP compatibility; not implemented |
| `--no_intergenic` | bool | false | Accepted for VEP compatibility; not implemented |
| `--hgnc` | bool | false | Accepted without effect: Ensembl VEP declares it and reads nothing from it |

### Allele Frequencies

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--af` | bool | false | Add global allele frequency from 1000 Genomes Phase 3 |
| `--af_1kg` | bool | false | Add continental allele frequencies from 1000 Genomes |
| `--af_gnomade` | bool | false | Add gnomAD exome allele frequencies |
| `--af_gnomadg` | bool | false | Add gnomAD genome allele frequencies |
| `--max_af` | bool | false | Add highest allele frequency across populations |
| `--pubmed` | bool | false | Add PubMed IDs for co-located variants |
| `--gene_phenotype` | bool | false | Add gene-phenotype associations |

### Filtering / Picking

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--pick` | bool | false | Pick one consequence per variant |
| `--pick_allele` | bool | false | Pick one consequence per variant allele |
| `--per_gene` | bool | false | Pick one consequence per gene |
| `--pick_allele_gene` | bool | false | Pick one consequence per allele per gene |
| `--most_severe` | bool | false | Report only the single most severe consequence |
| `--summary` | bool | false | Report a summary of all consequences |
| `--flag_pick` | bool | false | Flag (but don't remove) non-picked consequences |
| `--flag_pick_allele` | bool | false | Flag non-picked consequences per allele |
| `--flag_pick_allele_gene` | bool | false | Flag non-picked consequences per allele per gene |
| `--pick_order` | `mane_select`, `mane_plus_clinical`, `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `rank`, `length`, `ensembl`, `refseq` (comma-separated) | -- | Accepted for VEP compatibility; not implemented (the pick order is fixed to VEP's default) |
| `--allele_number` | bool | false | Add allele number from input VCF |
| `--show_ref_allele` | bool | false | Show reference allele in output |
| `--minimal` | bool | false | Convert alleles to minimal representation; switches on --allele_number as in Ensembl VEP |
| `--total_length` | bool | false | Accepted for VEP compatibility; not implemented |

### Plugins

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--plugin` | String (repeatable) | -- | Load plugin(s) - can be specified multiple times |
| `--custom` | String (repeatable) | -- | Custom annotation source(s) (accepted for VEP compatibility; not implemented, sources are not loaded) |
| `--custom_suppress_filter` | bool | false | Refused with an error: --custom_suppress_filter is not supported: it would keep a --custom VCF source's FILTER column out of the output unless the --custom source's fields argument names it, and custom annotation is not supported |
| `--dir_plugins` | String | -- | Plugin directory |

### Transcript Sets

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--refseq` | bool | false | Suffixes the species directory of the output header's cache path with _refseq when --json_cache is absent (accepted for VEP compatibility; selects no transcript set, which is the --json_cache contents) |
| `--merged` | bool | false | Suffixes the species directory of the output header's cache path with _merged when --json_cache is absent (accepted for VEP compatibility; selects no transcript set, which is the --json_cache contents) |
| `--gencode_basic` | bool | false | Accepted for VEP compatibility; not implemented (the transcript set is the --json_cache contents) |
| `--use_given_ref` | bool | false | Accepted without effect: vep-rs always uses the reference allele the input gives |

### Statistics

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--no_stats` | bool | false | Accepted for VEP compatibility; not implemented (no statistics file is written, with or without it) |
| `--stats_file` | String | -- | Accepted for VEP compatibility; not implemented (no statistics file is written) |
| `--stats_text` | bool | false | Accepted for VEP compatibility; not implemented (no statistics file is written) |

### Verbosity

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--verbose`, `-v` | bool | false | Verbose output (debug-level logging) |
| `--quiet`, `-q` | bool | false | Suppress status/progress messages |

### Miscellaneous

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--dont_skip` | bool | false | Abort the run on a variant that fails parsing instead of skipping it (VEP keeps the variant) |
| `--allow_non_variant` | bool | false | Process non-variant lines in VCF (e.g., reference-only calls) |
| `--mirna` | bool | false | Add the miRNA column (accepted for VEP compatibility; the value is not computed and stays empty) |
| `--lrg` | bool | false | Refused with an error: --lrg is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode) |
| `--tmpdir` | String | -- | Accepted without effect: no module of Ensembl VEP reads it, and vep-rs writes its Parquet intermediate beside the output file |

`--help` prints this list and exits 0, as does a command line with no argument at all (vep:29); `--version` prints `vep <version>` and exits 0.

## Option rules

Ensembl VEP's option rules are applied as data, from `Config.pm:327-645` at release/116.2, in this order after parsing: a flag the error table refuses is refused first; then the option sets switch on the flags they imply, in the order below, a set's value overwriting an explicit one (`apply_option_sets`, `Config.pm:816-884`); then the checks of `check_config` (`Config.pm:944-985`): valid values, incompatible pairs, requirements, deprecations, each refused with Ensembl VEP's own message. A rule naming a flag vep-rs does not declare cannot fire and is listed so the table reads as Ensembl VEP's.

### Option sets (`Config.pm:351-593`)

| Flag | Switches on |
| --- | --- |
| `--everything` | `--sift b`, `--polyphen b`, `--ccds`, `--hgvs`, `--symbol`, `--numbers`, `--domains`, `--regulatory`, `--canonical`, `--protein`, `--biotype`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--uniprot`, `--mane`, `--tsl`, `--appris`, `--variant_class`, `--gene_phenotype`, `--mirna` |
| `--hgvs` | `--hgvsc`, `--hgvsp` |
| `--genomes` | `--host mysql-eg-publicsql.ebi.ac.uk`, `--port 4157` |
| `--grch37` | `--port 3337`, `--assembly GRCh37` |
| `--refseq` | `--core_type otherfeatures` |
| `--individual` | `--allow_non_variant` |
| `--individual_zyg` | `--allow_non_variant` |
| `--cell_type` | `--regulatory` |
| `--regulatory_gff` | `--regulatory` |
| `--filter_common` | `--check_frequency` |
| `--check_frequency`, `--af`, `--af_1kg`, `--af_gnomad`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed` | `--check_existing` |
| `--show_cache_info` | `--offline`, `--format ensembl` |
| `--offline` | `--cache`, `--database 0` |
| `--cache` | `--no_slice_cache` |
| `--convert` | `--no_stats` |
| `--json` | `--rest`, `--output_format json` |
| `--rest` | `--no_escape` |
| `--vcf` | `--output_format vcf`, `--symbol`, `--biotype`, `--numbers` |
| `--tab` | `--output_format tab` |
| `--humdiv` | `--polyphen_analysis humdiv` |
| `--minimal` | `--allele_number` |
| `--gff` | `--custom file=%gff%,format=gff` |
| `--gtf` | `--custom file=%gtf%,format=gtf` |
| `--bigwig` | `--custom %bigwig%,,bigwig,exact` |
| `--phyloP` | `--custom %ucsc_data_root%%ucsc_assembly%/phyloP%phyloP%way/%ucsc_assembly%.phyloP%phyloP%way.bw,phlyoP%phyloP%way,bigwig,exact` |
| `--phastCons` | `--custom %ucsc_data_root%%ucsc_assembly%/phastCons%phastCons%way/%ucsc_assembly%.phastCons%phastCons%way.bw,phastCons%phastCons%way,bigwig,exact` |
| `--exclude_predicted` | `--transcript_filter not stable_id match ^X._` |
| `--bam` | `--use_transcript_ref`, `--bam_edited` |
| `--mane` | `--mane_select`, `--mane_plus_clinical` |
| `--use_given_ref` | `--use_transcript_ref 0` |

`--everything` therefore reaches `--hgvsc` and `--hgvsp` through `--hgvs`, `--mane_select` and `--mane_plus_clinical` through `--mane`, and `--check_existing` through the frequency flags; it does not switch on `--allele_number` or `--hgvsg`. `--regulatory` is switched on and annotates nothing, so every column of Ensembl VEP's `--everything` is produced except the regulatory columns, and no RegulatoryFeature or MotifFeature row is written.

### Incompatible flags (`Config.pm:615-632`)

`Can't use --<flag> and --<other> together`, naming the first pair in this table's order (Ensembl VEP iterates a hash there, so which pair it names varies from run to run); Ensembl VEP skips the check under `--safe`.

| Flag | Cannot be combined with |
| --- | --- |
| `--most_severe` | `--biotype`, `--no_intergenic`, `--protein`, `--symbol`, `--sift`, `--polyphen`, `--coding_only`, `--ccds`, `--mane`, `--mane_select`, `--canonical`, `--xref_refseq`, `--numbers`, `--domains`, `--tsl`, `--appris`, `--uniprot`, `--pick`, `--flag_pick`, `--pick_allele`, `--flag_pick_allele`, `--vcf`, `--summary` |
| `--summary` | `--biotype`, `--no_intergenic`, `--protein`, `--symbol`, `--sift`, `--polyphen`, `--coding_only`, `--ccds`, `--mane`, `--mane_select`, `--canonical`, `--xref_refseq`, `--numbers`, `--domains`, `--tsl`, `--appris`, `--uniprot`, `--pick`, `--flag_pick`, `--pick_allele`, `--flag_pick_allele`, `--vcf`, `--most_severe` |
| `--database` | `--af_1kg`, `--af_gnomad`, `--af_gnomade`, `--af_gnomadg`, `--af_exac`, `--max_af`, `--pubmed`, `--var_synonyms`, `--offline`, `--cache` |
| `--af_gnomade` | `--af_gnomad` |
| `--quiet` | `--verbose` |
| `--refseq` | `--gencode_basic`, `--merged` |
| `--gencode_primary` | `--gencode_basic`, `--refseq` |
| `--json` | `--vcf`, `--tab` |
| `--vcf` | `--json`, `--tab` |
| `--tab` | `--vcf`, `--json` |
| `--individual` | `--minimal`, `--individual_zyg` |
| `--check_ref` | `--lookup_ref` |
| `--check_svs` | `--offline` |
| `--ga4gh_vrs` | `--vcf` |
| `--shift_hgvs` | `--shift_3prime`, `--shift_genomic` |
| `--lrg` | `--offline` |

### Required flags (`Config.pm:607-612`)

`You must set --<required> to use --<flag>`.

| Flag | Needs |
| --- | --- |
| `--original` | `--filters` |
| `--phyloP` | `--ucsc_assembly` |
| `--phastCons` | `--ucsc_assembly` |
| `--custom_multi_allelic` | `--custom` |

### Valid values (`Config.pm:596-604`)

A declared flag carries its list as a value enumeration, so a value outside it is an argument syntax error (exit 2) listing the possible values.

| Flag | Values |
| --- | --- |
| `--format` | `ensembl`, `vcf`, `hgvs`, `id`, `spdi`, `region`, `guess` |
| `--terms` | `SO`, `display`, `NCBI` |
| `--sift` | `s`, `p`, `b` |
| `--polyphen` | `s`, `p`, `b` |
| `--pick_order` | `mane_select`, `mane_plus_clinical`, `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `rank`, `length`, `ensembl`, `refseq` |
| `--nearest` | `transcript`, `gene`, `symbol` |
| `--compress_output` | `gzip`, `bgzip` |

### Deprecated flags (`Config.pm:636-638`)

- `--convert`: `--convert has been deprecated`

### Lists (`Config.pm:327-346`)

- Comma-separated values: `--individual`, `--individual_zyg`, `--cell_type`, `--pick_order`, `--fields`, `--chr`, `--distance`, `--dont_export`.
- Repeatable, each occurrence appending: `--plugin`, `--custom`, `--phyloP`, `--phastCons`, `--transcript_filter`.
- Input formats Ensembl VEP exempts from its sorted-input check (`Config.pm:640-645`): `id`, `hgvs`, `spdi`, `region`. vep-rs performs no such check in any format: every record is annotated in input order, and an input whose chromosomes interleave keeps more transcript indexes in memory.

### Rules of vep-rs's own

- `ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT`

### Where `Config.pm` is not followed

- An output file matching `stdout` in any case does not force `--quiet` or drop `--verbose` (`Config.pm:933-937`): the status lines go to standard error, so the annotation on standard output is clean either way, and `--quiet --verbose` stays refused with `-o STDOUT` as without.
- A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not refused (`Config.pm:993-1017`): vep-rs reads a local cache in every run, so those flags select nothing and the run proceeds on the cache it is given.
- A `%name%` in an option-set value is replaced even when the flag's value is `0` or empty (`Config.pm:867` leaves the reference in place then); no declared flag carries such a value.
- A VCF output selected with `--output_format vcf` fills `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` as `--vcf` does; Ensembl VEP's `vcf` option set applies to the `--vcf` flag alone and leaves the CSQ fields empty under `--output_format vcf`.

## Flags refused or accepted without effect

The error table: every flag vep-rs refuses, and the flags accepted without effect whose only documentation is this table, with what Ensembl VEP does with each (cited to its source at release/116.2). A refused flag exits 1 with the message below under the `cli` component; an accepted one parses and changes nothing. A flag whose description above reads `Accepted for VEP compatibility` parses as well, its description saying what, if anything, it changes. The help text of each entry is this table's.

| Flag | Treatment | Message or reason | Ensembl VEP |
| --- | --- | --- | --- |
| `--database` | refused (`cli`) | `ERROR: [cli] --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:237 declares it; AnnotationSourceAdaptor.pm:146-200 opens the database |
| `--host` | refused (`cli`) | `ERROR: [cli] --host is not supported: it names the Ensembl database host, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:88, default ensembldb.ensembl.org at :287 |
| `--port` | refused (`cli`) | `ERROR: [cli] --port is not supported: it names the Ensembl database port, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:89, default 3306 at :289, 3337 under --grch37 at :401 |
| `--user` | refused (`cli`) | `ERROR: [cli] --user is not supported: it names the Ensembl database user, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:90 (user\|u=s), default anonymous at :288 |
| `--password` | refused (`cli`) | `ERROR: [cli] --password is not supported: it is the Ensembl database password, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:91 (password\|pass=s) |
| `--registry` | refused (`cli`) | `ERROR: [cli] --registry is not supported: it names an Ensembl registry file of database connections, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:87 |
| `--db_version` | refused (`cli`) | `ERROR: [cli] --db_version is not supported: it selects an Ensembl database version; --cache_version selects a cache` | Config.pm:92; CacheDir.pm:323 falls back to it for the cache version |
| `--genomes` | refused (`cli`) | `ERROR: [cli] --genomes is not supported: it selects the Ensembl Genomes database host; vep-rs reads a local cache only` | Config.pm:95; the option set at :390-396 sets the host and port |
| `--is_multispecies` | refused (`cli`) | `ERROR: [cli] --is_multispecies is not supported: it describes a multi-species Ensembl database; vep-rs reads a local cache only` | Config.pm:101; BaseVEP.pm:282 reads it for the registry species |
| `--lrg` | refused (`cli`) | `ERROR: [cli] --lrg is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode)` | Config.pm:221, incompatible with --offline at :631; Runner.pm:766 refuses it offline |
| `--convert` | refused (`cli`) | `ERROR: [cli] --convert has been deprecated` | Config.pm:636-638 lists it as deprecated and :977-985 refuses it with this text |
| `--gp` | refused (`cli`) | `ERROR: [cli] --gp is not supported: it reads 1000 Genomes pilot coordinates from INFO/GP, which no current VCF carries` | Config.pm:120; Parser/VCF.pm:271-285 takes chr:pos from INFO/GP and skips lines without it |
| `--ambiguous_hgvs` | refused (`cli`) | `ERROR: [cli] --ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g.) only` | Config.pm:202; Parser/HGVS.pm:178-186, 224 resolves the protein notation through the database |
| `--af_exac` (not declared; caught before parsing) | refused (`cli`) | `ERROR: [cli] --af_exac was removed in Ensembl VEP release 115: use --af_gnomade` | Config.pm:618 names it among --database's incompatibilities while @VEP_PARAMS (:73-276) does not declare it, so Ensembl VEP refuses it as an unexpected parameter (vep:30) |
| `--regulatory_gff` | refused (`cli`) | `ERROR: [cli] --regulatory_gff is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache, and regulatory annotation is not supported` | Config.pm:210; the option set at :434-439 switches on --regulatory |
| `--extended_promoters` | refused (`cli`) | `ERROR: [cli] --extended_promoters is not supported: it would widen the promoters read with --regulatory_gff to their extended bounds, and regulatory annotation is not supported` | Config.pm:211 |
| `--custom_suppress_filter` | refused (`cli`) | `ERROR: [cli] --custom_suppress_filter is not supported: it would keep a --custom VCF source's FILTER column out of the output unless the --custom source's fields argument names it, and custom annotation is not supported` | Config.pm:268 |
| `--condel` | accepted without effect | Ensembl VEP declares it and reads nothing from it | Config.pm:206 (condel=s), read by no module |
| `--hgnc` | accepted without effect | Ensembl VEP declares it and reads nothing from it | Config.pm:189, read by no module |
| `--keep_ann` | accepted without effect | Ensembl VEP declares it as a synonym of --keep_csq and reads nothing from it | Config.pm:220, read by no module |
| `--phased` | accepted without effect | Ensembl VEP marks the genotypes as phased and nothing reads the mark | Config.pm:146; Parser/VCF.pm:107, 654, 734 set phased on the variant and nothing reads it |
| `--sereal` | accepted without effect | the cache's serialiser is read from its own metadata | Config.pm:249; CacheDir.pm:147 takes the serialiser from the cache's info file |
| `--solr` | accepted without effect | Ensembl VEP declares it and writes no Solr output | Config.pm:215, read by no module |
| `--tmpdir` | accepted without effect | no module of Ensembl VEP reads it, and vep-rs writes its Parquet intermediate beside the output file | Config.pm:259 (tmpdir=s), default /tmp at :304, read by no module |
| `--offline` | accepted without effect | vep-rs always runs offline and never contacts a database | Config.pm:245; the option set at :464-469 switches on --cache and off --database |
| `--use_given_ref` | accepted without effect | vep-rs always uses the reference allele the input gives | Config.pm:255; the option set at :588-592 switches off use_transcript_ref, which CacheDir.pm:412-416 enables for a BAM-edited cache |
| `--no_check_variants_order` | accepted without effect | vep-rs never checks the input order: every record is annotated where it stands, and an input whose chromosomes interleave keeps more transcript indexes in memory | Config.pm:83; InputBuffer.pm:189-258 aborts after 100 records more than 100 bp behind their predecessor |

## Configuration sources

The command line is the one source of flags.

Environment variables vep-rs reads, which Ensembl VEP has no counterpart for:

- `VEP_FASTA_STORAGE`: how the FASTA is held: `mmap` (default) maps the file, `memory` reads it into RAM; any other value is an error.

## Exit status and errors

- `0`: the run completed; also `--help`, `--version` and an empty command line.
- `1`: a refusal or a failed run. One line `ERROR: [<component>] <message>` is written to standard error, where `<message>` is Ensembl VEP's text where it has one, without its own `ERROR: ` prefix. The components:
  - `cli`: an unsupported flag, or a combination the option rules refuse
  - `output`: the output file or format
  - `run`: any other failure of the annotation run, with its cause chain
- `2`: an argument syntax error, in clap's words with its `error:` prefix: an unknown flag, a missing value, a value outside an enumeration.

`--help` and `--version` win over a flag the error table refuses before parsing wherever they stand on the command line; over a flag the parser does not know they win only when they come first, since clap reads left to right (Ensembl VEP's Getopt::Long runs with `pass_through`, Config.pm:62, so any unknown flag reaches the usage there, vep:29).

A panic aborts the process with the panic's own status. Ensembl VEP's own exit codes (255 on `die`, 2 on a failed open) are not copied.
