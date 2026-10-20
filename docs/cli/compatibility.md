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

Six regulatory-feature terms need a regulatory-feature cache that is not loaded: `TFBS_ablation`, `TFBS_amplification`, `TF_binding_site_variant`, `regulatory_region_ablation`, `regulatory_region_amplification` and `regulatory_region_variant`. `--regulatory`, alone or through `--everything`, is accepted and warns that regulatory-feature annotation is not implemented.

The root term `sequence_variant` has no predicate, and Ensembl VEP's consequence engine never assigns it either. The other 34 of Ensembl VEP's 41 Sequence Ontology terms are assigned.

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
| `--convert` | `ERROR: [cli] --convert has been deprecated` | Config.pm:636-638 lists it as deprecated and :977-985 refuses it with this text |
| `--gp` | `ERROR: [cli] --gp is not supported: it reads 1000 Genomes pilot coordinates from INFO/GP, which no current VCF carries` | Config.pm:120; Parser/VCF.pm:271-285 takes chr:pos from INFO/GP and skips lines without it |
| `--ambiguous_hgvs` | `ERROR: [cli] --ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g., m.), transcript (c., n.) and protein (p.) notation, each resolved to one variant, a protein change with several equally short nucleotide paths skipped with Ensembl VEP's message` | Config.pm:202; Parser/HGVS.pm:178-186, 224 resolves the protein notation to every variant through the database; VariationFeatureAdaptor.pm:1956 refuses an ambiguous one without it |
| `--af_exac` (not declared; caught before parsing) | `ERROR: [cli] --af_exac was removed in Ensembl VEP release 115: use --af_gnomade` | Config.pm:618 names it among `--database`'s incompatibilities while @VEP\_PARAMS (:73-276) does not declare it, so Ensembl VEP refuses it as an unexpected parameter (vep:30) |
| `--regulatory_gff` | `ERROR: [cli] --regulatory_gff is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache, and regulatory annotation is not supported` | Config.pm:210; the option set at :434-439 switches on `--regulatory` |
| `--extended_promoters` | `ERROR: [cli] --extended_promoters is not supported: it would widen the promoters read with --regulatory_gff to their extended bounds, and regulatory annotation is not supported` | Config.pm:211 |
| `--custom_suppress_filter` | `ERROR: [cli] --custom_suppress_filter is not supported: it would keep a --custom VCF source's FILTER column out of the output unless the --custom source's fields argument names it, and custom annotation is not supported` | Config.pm:268 |

## Accepted without effect

Every flag below parses and changes nothing. The first table is the error table's: each flag's help is its reason, and the Ensembl VEP column cites what Ensembl VEP does with it at release/116.2.

| Flag | Reason | Ensembl VEP |
| --- | --- | --- |
| `--condel` | Ensembl VEP declares it and reads nothing from it | Config.pm:206 (condel=s), read by no module |
| `--hgnc` | Ensembl VEP declares it and reads nothing from it | Config.pm:189, read by no module |
| `--keep_ann` | Ensembl VEP declares it as a synonym of `--keep_csq` and reads nothing from it | Config.pm:220, read by no module |
| `--phased` | Ensembl VEP marks the genotypes as phased and nothing reads the mark | Config.pm:146; Parser/VCF.pm:107, 654, 734 set phased on the variant and nothing reads it |
| `--sereal` | the cache's serialiser is read from its own metadata | Config.pm:249; CacheDir.pm:147 takes the serialiser from the cache's info file |
| `--solr` | Ensembl VEP declares it and writes no Solr output | Config.pm:215, read by no module |
| `--tmpdir` | no module of Ensembl VEP reads it, and vep-rs writes its Parquet intermediate beside the output file | Config.pm:259 (tmpdir=s), default /tmp at :304, read by no module |
| `--offline` | vep-rs always runs offline and never contacts a database | Config.pm:245; the option set at :464-469 switches on `--cache` and off `--database` |
| `--use_given_ref` | vep-rs always uses the reference allele the input gives | Config.pm:255; the option set at :588-592 switches off use\_transcript\_ref, which CacheDir.pm:412-416 enables for a BAM-edited cache |
| `--no_check_variants_order` | vep-rs never checks the input order: every record is annotated where it stands, and an input whose chromosomes interleave keeps more transcript indexes in memory | Config.pm:83; InputBuffer.pm:189-258 aborts after 100 records more than 100 bp behind their predecessor |

The flags below are declared with the description their group page shows, which says what, if anything, each changes; the second column is what Ensembl VEP does with the flag.

| Flag | Ensembl VEP |
| --- | --- |
| `--regulatory` | Reports overlaps with regulatory regions and transcription factor binding sites |
| `--custom` | Adds annotation from a tabix-indexed or bigWig file |
| `--no_stats` | Writes no statistics file |
| `--stats_file` | Names the statistics file |
| `--stats_text` | Writes the statistics as plain text instead of HTML |

## Where `Config.pm` is not followed

- An output file matching `stdout` in any case does not force `--quiet` or drop `--verbose` (`Config.pm:933-937`). The status lines go to standard error, so the annotation on standard output is clean either way, and `--quiet --verbose` stays refused with `-o STDOUT` as without.
- A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not refused (`Config.pm:993-1017`). vep-rs reads a local cache in every run, so those flags select nothing and the run proceeds on the cache it is given.
- A `%name%` in an option-set value is replaced even when the flag's value is `0` or empty (`Config.pm:867` leaves the reference in place then). No declared flag carries such a value.
- A VCF output selected with `--output_format vcf` fills `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` as `--vcf` does. Ensembl VEP's `vcf` option set applies to the `--vcf` flag alone and leaves the CSQ fields empty under `--output_format vcf`.
- A value from the `--config` file, `vep.ini` or the environment never overwrites one a higher source gives, even a `0`, an empty value or a bare `--cache`. Ensembl VEP's `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a value Perl reads as false, so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`.
- A value from a file or the environment that the parser refuses (`fork many`, an empty `format`) is an argument syntax error (exit 2). Ensembl VEP's `is_valid_param` drops the line (`Config.pm:1046-1047`).
- A `--config` path that is a directory, or cannot be read for any reason, is the `config` error; Ensembl VEP's `open` succeeds on a directory and reads nothing from it. A `vep.ini` that is a directory is read as nothing here too. A byte outside UTF-8 in a value becomes the replacement character, where Ensembl VEP passes the byte.
- A flag read from a file or the environment acts exactly as if typed, refusals included. The status lines Ensembl VEP prints under `--verbose` while reading the sources (`Read configuration from <file>`, `Set '<flag>=<value>' from environment variable <VAR>`, `Config.pm:1107, 1145-1164`) are not printed.
- A leading `~` in `--config`, `--dir` and a file's `dir` is expanded to the home directory, which Ensembl VEP leaves to the shell.
- `-V` is the version, where Getopt::Long's `ignore_case` reads it as `-v`, and an upper-case short names no flag here. A short with its value attached (`-ofile.txt`) and a cluster of switches (`-qv`) are clap's readings, where Getopt::Long refuses the token. A `+` before an option name, which Getopt::Long accepts, is not an option here.

## Configuration sources

A flag is read from four sources, the first that names it winning, as Ensembl VEP reads them (`Config.pm:693-708`):

1. the command line;
2. the `--config` file;
3. `<dir>/vep.ini`, where `<dir>` is `--dir` (the command line's, else the `--config` file's) and `~/.vep` by default;
4. the environment variable `VEP_<FLAG>` (`VEP_DIR_PLUGINS` sets `--dir_plugins`).

A repeatable flag (`--plugin`, `--custom`, `--phyloP`, `--phastCons`, `--transcript_filter`) collects every source's values, the command line's first. No source behind the command line can set `--config`, `--help`, `--version`.

### The `--config` file and `vep.ini` (`Config.pm:1068-1110`)

- One flag per line: the flag's long name without its dashes (one leading dash is allowed), a blank, then its value. A line starting with `#`, or blank, is skipped.
- A switch is on with the value `1` and off with `0`. `cache` is on with `1`, a directory or an empty quoted value and off with `0`. A key alone on its line is dropped.
- A flag taking a value takes the first word after its name, or, for a repeatable flag, the rest of the line. Blanks inside double quotes are part of a value and the quotes are removed.
- A key is the long name exactly. A prefix, another case or a name that is no option of Ensembl VEP's or of this version names nothing, and the line is dropped as Ensembl VEP drops it. So is an alias, unless the alias is an option of Ensembl VEP's own (`full_cache_dir` names `--json_cache`).
- A key of an Ensembl VEP option this version lacks (`stats_html 1`) is refused as the flag on the command line is (exit 2). The exception is a switch of Ensembl VEP's set to `0`: off there, it adds nothing here.
- Of several lines naming one flag, the first whose value is not `0` or empty counts.
- A `--config` file that cannot be opened is `ERROR: [config] Could not open config file "<file>"`; `vep.ini` is read only when it is a file.

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

`VEP_<FLAG>=<value>` sets `--<flag>`, the name lower-cased, when the value is not empty. The name is read as a configuration file's key is: a name that is no option of Ensembl VEP's or of this version is ignored, and one of an Ensembl VEP option this version lacks is refused as the flag is (exit 2). The variables are read in name order.

Environment variables vep-rs reads outside this mechanism, which Ensembl VEP has no counterpart for:

| Variable | Meaning |
| --- | --- |
| `VEP_FASTA_STORAGE` | how the FASTA is held: `mmap` (default) maps the file and inflates a bgzipped one block by block as queries need them, `memory` reads it into RAM, a bgzipped one inflated whole; any other value is an error |

## Option spellings (`Config.pm:62`)

- Getopt::Long's defaults apply to the command line: a long option is read in any case (`--FORMAT vcf`), abbreviated to any prefix that names one option (`--form vcf`), and after one dash as well as two (`-cache`, `-form vcf`).
- A value follows its option or is attached with `=`, and the token after an option that requires a value is that value whatever it looks like (`--input_file -cache` names the file `-cache`).
- A prefix is resolved against this version's flags first. Where a flag vep-rs adds makes it ambiguous, it is resolved as Getopt::Long resolves it against Ensembl VEP's `VEP_PARAMS` list (`Config.pm:73-276`): `--js` is `--json` though `--json_cache` also starts with it, `--pa` `--password`, `--ve` `--verbose`.
- An Ensembl VEP option this version lacks is refused by its full name in every spelling (`-stats_html` and `--stats_h` are refused as `--stats_html`, exit 2), never read as a cluster of shorts.
- A prefix that names one flag of this version resolves to it even where Ensembl VEP's larger set leaves it ambiguous (`--sy` is `--symbol`). A prefix of several flags that Getopt::Long resolves no better (`--fo`), or a name no option has, is an argument syntax error (exit 2).
- A single-letter short is clap's and is matched by case: `-v` is `--verbose`, `-V` the version, `-qv` the cluster `-q -v`, and `-ofile.txt` the value `file.txt` of `-o`, a spelling Getopt::Long refuses.

Of the spellings Getopt::Long resolves to an option this version declares, these resolve differently here:

- none: every such spelling names the same option.

## Exit status and errors

- `0`: the run completed; also `--help`, `--version` and an empty command line.
- `1`: a refusal or a failed run. One line `ERROR: [<component>] <message>` is written to standard error, where `<message>` is Ensembl VEP's text where it has one, without its own `ERROR: ` prefix. The components are listed below.
- `2`: an argument syntax error, in clap's words with its `error:` prefix: an unknown flag, a prefix of several flags, a missing value, a value outside an enumeration.

| Component | Meaning |
| --- | --- |
| `cli` | an unsupported flag, or a combination the option rules refuse |
| `input` | the input: a file that does not exist, standard input without `--format`, a first data line of no known format, BCF input, CAid input (which Ensembl VEP resolves through its database), an identifier list against a cache without variation files, or a region or SPDI input without `--fasta` |
| `output` | the output file or format, the warnings file or the skipped-variants file |
| `plugin` | under `--safe`, a plugin that cannot be loaded or a plugins directory that does not exist |
| `cache` | the cache: nothing resolves at the configured path, its layout or metadata is unusable, its assembly disagrees with `--assembly` or with its own directory name, or a variation file cannot be opened or read |
| `fasta` | the reference FASTA: the file or directory named cannot be read as an indexed FASTA, or a flag that reads it (`--hgvs`, `--hgvsg`, `--shift_3prime 1`) was given and none was named or found in the cache directory |
| `filter` | a `--transcript_filter` Ensembl VEP would refuse, or whose regular expression is outside the supported dialect |
| `check_svs` | the `--check_svs` file: it is not a file, has no `.csi` or `.tbi` index beside it, or cannot be read |
| `config` | the `--config` file, or a `vep.ini` that is a file, cannot be read; or a configuration value the parser accepts but the run cannot use, such as a `--distance` part that is not an integer |
| `run` | any other failure of the annotation run, with its cause chain |

`--help` and `--version` win over a flag the error table refuses before parsing, wherever they stand on the command line. Over a flag the parser does not know they win only when they come first, since clap reads left to right. Ensembl VEP's Getopt::Long runs with `pass_through` (Config.pm:62), so there `--help` wins over an unknown flag wherever it stands (vep:29) and an unknown flag alone is refused (vep:30).

A panic aborts the process with the panic's own status. Ensembl VEP's own exit codes (255 on `die`, 2 on a failed open) are not copied.
