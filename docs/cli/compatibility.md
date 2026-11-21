# Compatibility with Ensembl VEP

This page gathers the Ensembl VEP behaviour vep-rs refuses or accepts without effect, and how the command line is read. It covers the database flags, the error table, where `Config.pm` is not followed, the configuration sources and their priority, the option spellings, and the exit status and error contract.

## Database flags

vep-rs reads a local cache and never contacts a database. Each database flag is declared so that the command line names it, and refused with the error its description shows.

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--database` | bool | false | Refused with an error: `--database` is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with `--json_cache`, or `--cache` with `--dir_cache` |
| `--host` | String | none | Refused with an error: `--host` is not supported: it names the Ensembl database host, and vep-rs reads a local cache only; name the cache with `--json_cache`, or `--cache` with `--dir_cache` |
| `--port` | String | none | Refused with an error: `--port` is not supported: it names the Ensembl database port, and vep-rs reads a local cache only; name the cache with `--json_cache`, or `--cache` with `--dir_cache` |
| `--user`, `-u` | String | none | Refused with an error: `--user` is not supported: it names the Ensembl database user, and vep-rs reads a local cache only; name the cache with `--json_cache`, or `--cache` with `--dir_cache` |
| `--password`, `--pass` | String | none | Refused with an error: `--password` is not supported: it is the Ensembl database password, and vep-rs reads a local cache only; name the cache with `--json_cache`, or `--cache` with `--dir_cache` |
| `--registry` | String | none | Refused with an error: `--registry` is not supported: it names an Ensembl registry file of database connections, and vep-rs reads a local cache only; name the cache with `--json_cache`, or `--cache` with `--dir_cache` |
| `--db_version` | String | none | Refused with an error: `--db_version` is not supported: it selects an Ensembl database version; `--cache_version` selects a cache |
| `--genomes` | bool | false | Refused with an error: `--genomes` is not supported: it selects the Ensembl Genomes database host; vep-rs reads a local cache only |
| `--is_multispecies` | String | none | Refused with an error: `--is_multispecies` is not supported: it describes a multi-species Ensembl database; vep-rs reads a local cache only |

## Consequence terms never assigned

One term is never assigned: the root `sequence_variant` has no predicate, and Ensembl VEP's consequence engine never assigns it either. The other 40 of Ensembl VEP's 41 Sequence Ontology terms are assigned; the six regulatory and motif terms need `--regulatory` (or `--everything`) and a cache converted with its regulatory features ([Regulatory and motif features](../regulatory.md)).

## Input formats

The input formats `--format` reads, how `guess` detects one, and the formats refused are stated in its row under [Input](input-output.md#input).

## Refused flags

Every flag vep-rs refuses, with the line it prints and what Ensembl VEP does with the flag (cited to its source at release/116.2). A refused flag exits 1 with the message below under the `cli` component, and the help text of a declared one is this table's message.

| Flag | Error | Ensembl VEP |
| --- | --- | --- |
| `--database` | `ERROR: [cli] --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:237 declares it; AnnotationSourceAdaptor.pm:146-200 opens the database |
| `--host` | `ERROR: [cli] --host is not supported: it names the Ensembl database host, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:88, default ensembldb.ensembl.org at :287 |
| `--port` | `ERROR: [cli] --port is not supported: it names the Ensembl database port, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:89, default 3306 at :289, 3337 under `--grch37` at :401 |
| `--user` | `ERROR: [cli] --user is not supported: it names the Ensembl database user, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:90 (user\|u=s), default anonymous at :288 |
| `--password` | `ERROR: [cli] --password is not supported: it is the Ensembl database password, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:91 (password\|pass=s) |
| `--registry` | `ERROR: [cli] --registry is not supported: it names an Ensembl registry file of database connections, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache` | Config.pm:87 |
| `--db_version` | `ERROR: [cli] --db_version is not supported: it selects an Ensembl database version; --cache_version selects a cache` | Config.pm:92; CacheDir.pm:323 falls back to it for the cache version |
| `--genomes` | `ERROR: [cli] --genomes is not supported: it selects the Ensembl Genomes database host; vep-rs reads a local cache only` | Config.pm:95; the option set at :390-396 sets the host and port |
| `--is_multispecies` | `ERROR: [cli] --is_multispecies is not supported: it describes a multi-species Ensembl database; vep-rs reads a local cache only` | Config.pm:101; BaseVEP.pm:282 reads it for the registry species |
| `--lrg` | `ERROR: [cli] --lrg is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode)` | Config.pm:221, incompatible with `--offline` at :631; Runner.pm:766 refuses it offline |
| `--convert` | `ERROR: [cli] --convert has been deprecated` | Config.pm:636-638 lists it as deprecated and :977-985 refuses it with this text when its value is true; a value Perl reads as false, `0` or the empty value, leaves the flag inactive (:902-906) and the run proceeds, there and here. Behind such a value, a true one from a lower source is refused there (:1101, 1158) and not read here |
| `--gp` | `ERROR: [cli] --gp is not supported: it reads a GRCh37 chr:pos from the INFO/GP field that some 1000 Genomes VCFs carry` | Config.pm:120; Parser/VCF.pm:271-285 takes chr:pos from INFO/GP and skips lines without it |
| `--ambiguous_hgvs` | `ERROR: [cli] --ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g., m.), transcript (c., n.) and protein (p.) notation, each resolved to one variant, a protein change with several equally short nucleotide paths skipped with Ensembl VEP's message` | Config.pm:202; Parser/HGVS.pm:178-186, 224 resolves the protein notation to every variant through the database; VariationFeatureAdaptor.pm:1956 refuses an ambiguous one without it |
| `--af_exac` (not declared; caught before parsing) | `ERROR: [cli] --af_exac was removed in Ensembl VEP release 107: use --af_gnomade` | Config.pm:618 names it among `--database`'s incompatibilities while @VEP\_PARAMS (:73-276) does not declare it, so Ensembl VEP refuses it as an unexpected parameter (vep:30) |
| `--gff` | `ERROR: [cli] --gff is not supported: Ensembl VEP reads a GFF3 as its transcript set through --custom file=<gff>,format=gff; build a cache from the file and the genome, then annotate against it: vep-cache-builder --species <species> --assembly <assembly> --release <N> --gtf <file> --genome-fasta <genome.fa> --output-dir <dir> (--gff3 <file> for a GFF3), then vep -i <in> -o <out> --offline --json_cache <dir> --fasta <genome.fa>; the page docs/cache-from-gtf-or-gff3.md states the route` | Config.pm:260 (gff=s); the option set at :531-536 sets `--custom` file=%gff%,format=gff; AnnotationSource/File/GFF.pm:71-116 lists the feature types read |
| `--gtf` | `ERROR: [cli] --gtf is not supported: Ensembl VEP reads a GTF as its transcript set through --custom file=<gtf>,format=gtf; build a cache from the file and the genome, then annotate against it: vep-cache-builder --species <species> --assembly <assembly> --release <N> --gtf <file> --genome-fasta <genome.fa> --output-dir <dir> (--gff3 <file> for a GFF3), then vep -i <in> -o <out> --offline --json_cache <dir> --fasta <genome.fa>; the page docs/cache-from-gtf-or-gff3.md states the route` | Config.pm:261 (gtf=s); the option set at :538-543 sets `--custom` file=%gtf%,format=gtf; AnnotationSource/File/GTF.pm:78-85 lists the feature types read |
| `--bam` | `ERROR: [cli] --bam is not supported: Ensembl VEP applies a BAM's alignments to RefSeq transcript models at run time; the RefSeq cache Ensembl distributes already carries those edits, which vep-rs reads from a converted homo_sapiens_refseq cache under --refseq; a GFF or GTF transcript set of your own is read without the BAM's edits; build a cache from the file and the genome, then annotate against it: vep-cache-builder --species <species> --assembly <assembly> --release <N> --gtf <file> --genome-fasta <genome.fa> --output-dir <dir> (--gff3 <file> for a GFF3), then vep -i <in> -o <out> --offline --json_cache <dir> --fasta <genome.fa>; the page docs/cache-from-gtf-or-gff3.md states the route` | Config.pm:253 (bam=s); the option set at :573-579 switches on use\_transcript\_ref and bam\_edited; AnnotationType/Transcript.pm:421-631 applies the alignments to each transcript as SeqEdits; CacheDir.pm:412-416 reads a BAM-edited cache's info |

## Accepted without effect

Every flag in the table below parses and changes nothing. The table is the error table's: each flag's help is its reason, and the Ensembl VEP column cites what Ensembl VEP does with it at release/116.2.

| Flag | Reason | Ensembl VEP |
| --- | --- | --- |
| `--condel` | Ensembl VEP declares it and reads nothing from it | Config.pm:206 (condel=s), read by no module |
| `--debug` | Ensembl VEP reads it only to print the plugin directories and modules it loads; the annotation output is unchanged | Config.pm:275; Runner.pm:920, 941 |
| `--hgnc` | Ensembl VEP declares it and reads nothing from it | Config.pm:189, read by no module |
| `--keep_ann` | Ensembl VEP declares it as a synonym of `--keep_csq` and reads nothing from it | Config.pm:220, read by no module |
| `--phased` | Ensembl VEP marks the genotypes as phased and nothing reads the mark | Config.pm:146; Parser/VCF.pm:107, 654, 734 set phased on the variant and nothing reads it |
| `--sereal` | the cache's serialiser is read from its own metadata | Config.pm:249; CacheDir.pm:147 takes the serialiser from the cache's info file |
| `--solr` | Ensembl VEP declares it and writes no Solr output | Config.pm:215, read by no module |
| `--tmpdir` | no module of Ensembl VEP reads it, and it places no file vep-rs writes | Config.pm:259 (tmpdir=s), default /tmp at :304, read by no module |
| `--offline` | vep-rs always runs offline and never contacts a database | Config.pm:245; the option set at :464-469 switches on `--cache` and off `--database` |

## Where `Config.pm` is not followed

- An output file matching `stdout` in any case does not force `--quiet` or drop `--verbose` (`Config.pm:933-937`). The status lines go to standard error, so the annotation on standard output is clean either way, and `--quiet --verbose` stays refused with `-o STDOUT` as without. `--quiet` silences those lines and nothing else: the statistics report is written in full under it, as Ensembl VEP keys its counters on `--no_stats` alone (`Stats.pm:230, 277`).
- A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not refused (`Config.pm:993-1017`). vep-rs reads a local cache in every run, so those flags select nothing and the run proceeds on the cache it is given.
- A `%name%` in an option-set value is replaced even when the flag's value is `0` or empty (`Config.pm:867` leaves the reference in place then). No declared flag carries such a value.
- A VCF output selected with `--output_format vcf` fills `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` as `--vcf` does. Ensembl VEP's `vcf` option set applies to the `--vcf` flag alone and leaves the CSQ fields empty under `--output_format vcf`.
- A value from the `--config` file, `vep.ini` or the environment never overwrites one a higher source gives, even a `0`, an empty value or a bare `--cache`. Ensembl VEP's `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a value Perl reads as false, so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`.
- A value from a file or the environment that the parser refuses (`fork many`, an empty `format`) is an argument syntax error (exit 2). Ensembl VEP's `is_valid_param` drops the line (`Config.pm:1046-1047`). It keeps a line whose value is itself an option Getopt::Long reads (`shift_3prime --offline`), the value stored as text, which Perl reads as true.
- `--buffer_size 0` annotates every record, where Ensembl VEP's buffer admits a record only while it holds fewer than the value (`InputBuffer.pm:191, 212`), so there a zero writes the header alone. A value with a minus sign or an underscore (`-1`, `1_0`) is an argument syntax error (exit 2); Getopt::Long reads both as integers, `-1` admitting no record and `1_0` read as 10.
- A byte outside UTF-8 in a value from a file or the environment becomes the replacement character, where Ensembl VEP passes the byte. A `--config` value is a path as written: `-`, surrounding whitespace, a `<`-led or `>`-led value and one opening or closing with `|` name a file so called. Ensembl VEP opens the value with Perl's two-argument `open` (`Config.pm:1075`). That reads `-` as standard input, trims the whitespace, opens a `<`-led or `+<`-led value as the file after the mode, opens a `>`-led value for writing, emptied, or appended to after `>>`, and runs a `|`-led or `|`-closed value as a command.
- A flag read from a file or the environment acts exactly as if typed, refusals included. So `show_cache_info` from a file or the environment prints the cache information and exits 0. Ensembl VEP's driver reads the typed flag alone (vep:39), and a lower source's value runs the annotation there instead. The status lines Ensembl VEP prints under `--verbose` while reading the sources (`Read configuration from <file>`, `Set '<flag>=<value>' from environment variable <VAR>`, `Config.pm:1107, 1145-1164`) are not printed.
- A leading `~` in `--config`, `--dir`, `--dir_cache`, `--cache <dir>`, `--json_cache`, `--synonyms`, `--check_svs` and a file's `dir` is expanded to the home directory. Ensembl VEP leaves the expansion to the shell. A `~/` the shell did not expand (quoted, or from a file or the environment) in `--config`, `--dir`, `--dir_cache`, `--cache`, `--synonyms` or `--full_cache_dir` (`--json_cache` here) is passed on as typed there and names no configuration file, cache directory or synonyms file. The plugins directory under `--dir` alone is globbed there (`Runner.pm:908`). `--check_svs` takes no value there (`Config.pm:125`).
- `--freq_gt_lt` takes `gt` or `lt` and `--freq_filter` `exclude` or `include`; any other word is an argument syntax error (exit 2). Ensembl VEP declares both as `=s` (`Config.pm:137, 139`) and checks neither: an unknown direction passes no allele (`BaseCacheVariation.pm:467-468`), and any filter but `include` reads as `exclude` (`:319`).
- `--shift_hgvs` takes `0` or `1`; any other value is an argument syntax error (exit 2). Ensembl VEP's `shift_hgvs=i` (`Config.pm:201`) takes any integer and `BaseRunner.pm:498` matches it against `/(0|1)/` as a substring. So `2` leaves the default shift on, and `10` sets the shifting variables to 10, which shifts HGVSc and HGVSp (tested against `0`) and not HGVSg (tested against `1`).
- `--max_sv_size` takes `-1` or a size of 0 or more; a value below `-1` is an argument syntax error (exit 2), a negative size being no size. Ensembl VEP's `max_sv_size=i` (`Config.pm:154`) takes any integer and `Parser.pm:491` spares `-1` alone, so `-2` makes every record oversize and the run drops them all, writing the header alone.
- An integer written with an underscore on `--failed`, `--clin_sig_allele`, `--fork`, `--cache_version` or `--max_sv_size` (`--fork 1_0`) is an argument syntax error (exit 2). Getopt::Long drops the underscores of an `=i` value (`Config.pm:119, 147, 151, 154, 239`), so `1_0` reads as 10 there and `--cache_version 11_6` selects release 116. `--shift_3prime` and `--shift_genomic` read an integer as Getopt::Long spells one, the underscores included.
- `-V` is the version, where Getopt::Long's `ignore_case` reads it as `-v`, and an upper-case short names no flag here. A short with its value attached (`-ofile.txt`) and a cluster of switches (`-qv`) are clap's readings, where Getopt::Long refuses the token. A `+` before an option name, which Getopt::Long accepts, is not an option here.
- A flag written with `=` and nothing after it (`--tmpdir=`, `--convert=`) carries the empty value, as `--tmpdir ''` does. Getopt::Long leaves such a token of a flag that requires a value unread, and Ensembl VEP refuses it as an extra parameter (vep:30).
- A trailing comma or a comma-joined value of `--sift`, `--polyphen`, `--terms`, `--nearest`, `--compress_output` or `--format` (`--sift b,`) is an argument syntax error (exit 2). Ensembl VEP's `check_config` splits such a value on commas before its `%VALID` test (`Config.pm:946`) and reads the first token by prefix (`OutputFactory.pm:1808`), so `--sift b,` annotates there (`SIFT=tolerated(0.19)` on one transcript of `21 14216193 C G`). A `0` or empty value of such a flag (`--sift 0`, `--sift ''`, `--polyphen 0`, `--nearest 0`) is an argument syntax error too (exit 2). Ensembl VEP reads a value Perl reads as false as the flag unset (`Config.pm:904-905, 945`) and annotates without the column (`OutputFactory.pm:1805`, `AnnotationType/Transcript.pm:170`).
- The `-o` value names the output file as given. Ensembl VEP opens it with Perl's two-argument `open` (`BaseRunner.pm:226`), which strips surrounding whitespace and reads a leading `>` as the append mode, so `-o ' lead.txt'` writes `lead.txt` there and ` lead.txt` here, and `-o '>out.txt'` appends to `out.txt` there and writes `>out.txt` here.
- The `--synonyms` value names its file as given too, `-`, `0`, the empty value, a value with surrounding whitespace and a `>`-led value included. Ensembl VEP's two-argument `open` (`BaseVEP.pm:644`) reads `-` as standard input, trims the whitespace, and opens a `>`-led value for writing, emptying the file it names and reading an empty table. A value Perl reads as false leaves the flag unset there, so the cache's own file is read (`CacheDir.pm:254`).
- A `--fasta` value Perl reads as false, `0` or the empty value, leaves that flag unset there likewise, so the cache directory's FASTA is auto-detected or the run proceeds without one (`CacheDir.pm:235`, `BaseVEP.pm:521`). vep-rs opens the value as a path and refuses it with `Specified FASTA file/directory <value> not found` (exit 1).
- The `--warning_file` and `--skipped_variants_file` values name their files as given likewise, `-`, surrounding whitespace and a `>`-led value included. Ensembl VEP opens each for writing with the two-argument `open` (`BaseVEP.pm:894`, `BaseRunner.pm:344`), which sends `-` to standard output, trims the whitespace and appends a `>`-led value to the file it names. So `--warning_file -` prints the warnings on standard output there and writes a file named `-` here, `--warning_file ' sp.txt'` writes `sp.txt` there and ` sp.txt` here, and `--skipped_variants_file '>skipped.txt'` appends to `skipped.txt` there and writes `>skipped.txt` here.
- A `--skipped_variants_file` value Perl reads as false, `0` or the empty value, leaves the flag unset there (`Runner.pm:176`), so no file is written. vep-rs writes a file named `0`, and the empty value fails the run as it ends, the output written, with `Could not write to skipped variants file` (exit 1). A false `--warning_file` value is the default warnings file in both engines (`BaseVEP.pm:890`).

## The statistics report

The reports `--stats_text` and `--stats_html` write, the HTML one by default, have Ensembl VEP's sections, rows and charts (`Stats.pm` `dump_text` and `dump_html`), with these differences:

- The HTML report has the same sections, ids, titles, tables and sort orders as Ensembl VEP's (`Stats.pm` `dump_html`), rendered as inline SVG without Ensembl VEP's remote chart loader and images. The page is one file that loads no script, no chart library and no image, so it opens without a network connection. Each pie, bar and area chart is drawn from the same counts, with Ensembl VEP's colours where its `COLOUR_KEYS` names the key (`Constants.pm:299-346`); a click on a chart selects nothing, where the chart library highlights the table row.
- The `Start time` and `End time` rows are UTC, the clock of the `## Output produced at` header line, where Ensembl VEP prints local time (`Utils.pm:467-488`); `Run time` is whole seconds in both.
- The `Command line options` row is the `## VEP command-line:` header's payload, so it names `--json_cache` where Ensembl VEP names `--dir_cache` and carries no `--database 0` (`vep:35`).
- A chart's rows are ordered by count or chromosome as Ensembl VEP orders them (`Stats.pm:760-790`), and two rows its order leaves tied (equal counts; `21` and `chr21`) print in key order, where Ensembl VEP prints them in its hash order, which changes between runs.
- The `Distribution of variants on chromosome` charts need the chromosome lengths, which vep-rs reads from the index of the `--fasta` file alone, one chart per chromosome spelled as the FASTA names it. Ensembl VEP reads them from its database when one is connected and from the FASTA otherwise (`BaseVEP.pm:550-580`). Without a FASTA neither engine writes a distribution chart.
- The report is checked for an existing file before the output and the skipped-variants file are, as Ensembl VEP checks it (`Runner.pm:144-148`, before `:202`). The three files are created only once every check passes, so a refused output leaves no empty report behind, where Ensembl VEP leaves one.

## Configuration sources

A flag of `vep`, and an option of `vep recode`, is read from four sources, the first that names it winning, as Ensembl VEP reads them for `vep` and `variant_recoder` alike (`Config.pm:693-708`):

1. the command line;
2. the `--config` file;
3. `<dir>/vep.ini`, where `<dir>` is `--dir` (the command line's, else the `--config` file's) and `~/.vep` by default;
4. the environment variable `VEP_<FLAG>` (`VEP_DIR_PLUGINS` sets `--dir_plugins`).

A repeatable flag (`--plugin`, `--custom`, `--phyloP`, `--phastCons`, `--transcript_filter`) collects every source's values, the command line's first. No source behind the command line can set `--config`, `--help`, `--version`.

### The `--config` file and `vep.ini` (`Config.pm:1068-1110`)

- One flag per line: the flag's long name without its dashes (one leading dash is allowed), a blank, then its value. A line starting with `#`, or blank, is skipped.
- A switch is on with the value `1` and off with `0`, except `var_synonyms`, whose `0` switches `check_existing` on, with the notice `INFO: Colocated variants will be provided when using --var_synonyms`, and writes no synonyms, as Ensembl VEP tests that flag with `defined` (Runner.pm:117). `cache` is on with `1`, a directory or an empty quoted value and off with `0`. A key alone on its line is dropped. Ensembl VEP drops a scalar or switch key alone on its line too (`Config.pm:1101`) and keeps an empty value for a repeatable key alone on its line (`Config.pm:1094-1097`): `plugin` warns `Failed to compile plugin` there, and `custom` and `transcript_filter` are fatal.
- A flag taking a value takes the first word after its name, or, for a repeatable flag, the rest of the line. Blanks inside double quotes are part of a value and the quotes are removed. A literal `___SPACE___` in a value stays as written. Ensembl VEP's reader uses that text as its placeholder for a blank inside double quotes and turns a literal one into a blank too (`Config.pm:1082, 1089`), so a file named with it is unreachable there.
- `shift_3prime` and `shift_genomic` are on with any integer but `0`, read as text: `00`, `+0` and `-0`, which are 0 and off on the command line, are on from a file or the environment. Ensembl VEP keeps a lower source's value as the text it reads (`Config.pm:1101, 1158`) and tests it for truth, where Getopt::Long stores the command line's value as a number.
- A key is the long name exactly. A prefix, another case or a name that is no option of Ensembl VEP's or of this version names nothing, and the line is dropped. So is an alias, unless the alias is an option of Ensembl VEP's own (`full_cache_dir` names `--json_cache`). Ensembl VEP drops an unknown name too. A prefix or another case of one of its options passes `is_valid_param` (Getopt::Long's `ignore_case` and `auto_abbrev`) and is stored under the written key, which nothing reads: the annotation is unchanged, but the `## VEP command-line:` header echoes the key (`--PROTEIN_VERSION` for `PROTEIN_VERSION 1`, `Config.pm:734-745`).
- A key of a flag this version refuses (`database 1`, `af_exac 1`) reaches that refusal (exit 1). The exception is a switch set to `0` or the empty value: off there, it adds nothing here, and `af_exac 0` adds nothing either. `convert 0` is read as `--convert 0` is, and the run proceeds.
- Of several lines naming one flag, the first whose value is not `0` or empty counts.
- A `--config` file that cannot be opened is `ERROR: [config] Could not open config file "<file>"`. A `--config` value of `0` or the empty value is that error with the empty name, before any file is looked for (`Config.pm:1073`). A `--config` path that is a directory reads as nothing, as Perl's `open` of a directory yields no line (`:1075-1077`). `vep.ini` is read only when it is a file.

```
# flags of a run
cache 1
offline 1
dir_cache /data/vep
species homo_sapiens
assembly GRCh38
plugin LoFtool,/data/LoFtool_scores.txt
```

### Environment variables (`Config.pm:1127-1168`)

`VEP_<FLAG>=<value>` sets `--<flag>`, the name lower-cased, when the value is not empty. The name is read as a configuration file's key is: a name that is no option of Ensembl VEP's or of this version is ignored. One of a flag this version refuses reaches that refusal (`VEP_DATABASE=1` is the `--database` error, exit 1). The variables are read in name order.

Environment variables vep-rs reads outside this mechanism, which Ensembl VEP has no counterpart for:

| Variable | Meaning |
| --- | --- |
| `VEP_FASTA_STORAGE` | how the FASTA is held: `mmap` (default) maps the file and inflates a bgzipped one block by block as queries need them, `memory` reads it into RAM, a bgzipped one inflated whole; any other value is an error |

## Option spellings (`Config.pm:62`)

- Getopt::Long's defaults apply to the command line: a long option is read in any case (`--FORMAT vcf`), abbreviated to any prefix that names one option (`--form vcf`), and after one dash as well as two (`-cache`, `-form vcf`).
- A value follows its option or is attached with `=`, and the token after an option that requires a value is that value whatever it looks like (`--input_file -cache` names the file `-cache`).
- A prefix is resolved against this version's flags first. Where a flag vep-rs adds makes it ambiguous, it is resolved as Getopt::Long resolves it against Ensembl VEP's `VEP_PARAMS` list (`Config.pm:73-276`): `--js` is `--json` though `--json_cache` also starts with it, `--pa` `--password`, `--ve` `--verbose`.
- A prefix that names one flag of this version resolves to it (`--bio` is `--biotype`). A prefix of several flags that Getopt::Long resolves no better (`--fo`), or a name no option has, is an argument syntax error (exit 2).
- A single-letter short is clap's and is matched by case: `-v` is `--verbose`, `-h` the help and `-V` the version, `-qv` the cluster `-q -v`, and `-ofile.txt` the value `file.txt` of `-o`, a spelling Getopt::Long refuses. Ensembl VEP refuses `-h`, a prefix of eight of its options, and has no `--version`.

Of the spellings Getopt::Long resolves to an option this version declares, these resolve differently here:

- none: every such spelling names the same option.

## Exit status and errors

- `0`: the run completed; also `--help`, `--version` and an empty command line.
- `1`: a refusal or a failed run. One line `ERROR: [<component>] <message>` is written to standard error, where `<message>` is Ensembl VEP's text where it has one, without its own `ERROR: ` prefix; a refused `--custom` string carries Ensembl VEP's second line, `LINE: --custom <string>`, where Ensembl VEP prints it. The components are listed below.
- `2`: an argument syntax error, in clap's words with its `error:` prefix: an unknown flag, a prefix of several flags, a missing value, a value outside an enumeration.
- A compressed input cut inside a gzip member is a `run` failure, exit 1 with nothing past the header written, where Ensembl VEP writes the records before the cut and exits 0 ([the divergence page](../divergences/compressed-input.md)).
- A directory named as the input is a `run` failure, exit 1 with no output written, where Ensembl VEP passes the existence test (`Parser.pm:239`), leaves the directory uncompressed (`:242` tests `-f`), opens it, reads nothing from it and exits 0 with the header alone written.
- A line of the run that standard error refuses (the banner, a status, progress or summary line, a warning's echo, the `ERROR:` line) is dropped and never changes the exit status. A completed run exits 0 on a full device or a pipe without a reader, and a refusal exits 1. Ensembl VEP completes on a full device too; on a pipe without a reader it completes when it has nothing to write there and is killed by `SIGPIPE` (exit 141, no output) at its first warning (`BaseVEP.pm:821`).
- An unwritable `--warning_file` is an `output` failure, exit 1, raised after the batch that warned, with that batch's rows written and each warning echoed on standard error unless `--quiet`; under `--show_cache_info` it is raised before the print, with nothing printed. Ensembl VEP dies inside the first warning (`BaseVEP.pm:894`), with the header alone written and no warning echoed.
- The process ends as soon as its outputs are written and closed: the loaded cache is left to the operating system, which reclaims it in one step, instead of being freed one allocation at a time. An allocator asked to report at exit (`MIMALLOC_SHOW_STATS`, `MIMALLOC_VERBOSE`) gets the ordinary exit, so its report prints.

| Component | Meaning |
| --- | --- |
| `cli` | an unsupported flag, or a combination the option rules refuse |
| `input` | the input: a file that does not exist, standard input without `--format`, a first data line of no known format, BCF input, CAid input (looked up in the Ensembl variation database, which vep-rs does not reach), an identifier list against a cache without variation files, a region, SPDI or HGVS input without `--fasta`, or a region input with `--check_ref` |
| `output` | the output file or format, the warnings file or the skipped-variants file |
| `parquet` | the Parquet writer: input out of chromosome and position order without `--no_check_variants_order`, or an output directory holding files vep-rs did not write |
| `plugin` | a `--plugin` naming an Ensembl VEP plugin vep-rs does not provide, or a Perl module (a `.pm` file by path or under the plugins directories), with or without `--safe`; under `--safe`, a plugin that cannot be loaded or a plugins directory that does not exist |
| `cache` | the cache: nothing resolves at the configured path, its layout or metadata is unusable, a transcript or regulatory shard cannot be read (a transcript an HGVS input names included), or a regulatory chromosome directory holds files but no shard, its assembly disagrees with `--assembly` or with its own directory name, or a variation file cannot be opened or read |
| `fasta` | the reference FASTA: the file or directory named cannot be read as an indexed FASTA, or a flag that reads it (`--hgvs`, `--hgvsg`, `--check_ref`, `--lookup_ref`, `--use_transcript_ref`, `--shift_3prime 1`) was given and none was named or found in the cache directory |
| `filter` | a `--transcript_filter` Ensembl VEP would refuse, or whose regular expression is outside the supported dialect |
| `check_svs` | the `--check_svs` file: it is not a file, has no `.csi` or `.tbi` index beside it, or cannot be read |
| `custom` | a `--custom` source: its string is malformed (Ensembl VEP's text), it names what vep-rs does not read (a remote URL, or a GFF or GTF source in transcript mode), or its file cannot be opened or read |
| `config` | the `--config` file, or a `vep.ini` that is a file, cannot be read; or a configuration value the parser accepts but the run cannot use, such as a `--distance` part that is not an integer |
| `run` | any other failure of the annotation run, with its cause chain |

`--help` and `--version` win over every other token of the command line, wherever they stand before a `--` (vep:29 before vep:30 under Getopt::Long's `pass_through`, Config.pm:62); `--help=1` and a `--help` behind `--` are refused. A `--help` or `--version` that is the token after a flag taking a string is that value (`-i --help` names the file) on both engines. After a flag taking a number it is the help or the version, as it is after a flag of this version's own, which Ensembl VEP leaves in `@ARGV` (`--pipeline_depth --help`).

A panic aborts the process with the panic's own status. Ensembl VEP's own exit codes (255 on `die`, 2 on a failed open) are not copied.
