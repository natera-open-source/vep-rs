// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CLI reference pages (`docs/cli-reference.md` and the group pages under
//! `docs/cli/`) are rendered from the command-line definition, the option-rule
//! tables and the error table, so the pages, `--help` and the parser cannot
//! disagree. The test asserts every committed page is its rendering;
//! `VEP_UPDATE_CLI_REFERENCE=1` rewrites them instead.

use std::any::TypeId;
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use clap::{Arg, ArgAction, Command, CommandFactory};
use vep_cli::args::Args;
use vep_cli::argv::{
    parser_options, perl_options, prepass, resolve, resolve_spelling, spell_as_declared,
    unresolved_perl_spellings, COMMAND_LINE_ONLY, DEFAULT_DIR, ENV_PREFIX, INI_FILE,
    VERSION_SHORT_BESIDE_ARGUMENTS,
};
use vep_cli::config::{Treatment, COMPONENTS, ENVIRONMENT, ERROR_TABLE};
use vep_cli::filter_cli::FilterArgs;
use vep_cli::option_rules::{
    ALLOW_MULTIPLE, DEPRECATED, INCOMPATIBLE, INCOMPATIBLE_LIFTED, LIST_FLAGS, OPTION_SETS,
    OWN_RULES, REQUIRES, UNSORTABLE, VALID,
};
use vep_cli::recode_cli::RecodeArgs;

const LANDING: &str = "docs/cli-reference.md";

/// The landing page's note on `--cell_type`, the one list flag whose value is
/// not split on every comma: the cache's epigenome names may hold one.
const CELL_TYPE_LIST_NOTE: &str = " `--cell_type` is read against the cache's epigenome names by \
longest match, since a name may hold a comma \
([`--cell_type` names with a comma](divergences/cell-type-comma-names.md)).";

const COMPATIBILITY: &str = "docs/cli/compatibility.md";

/// The page of the companion tool `vep filter`, rendered from its own command
/// line ([`FilterArgs`]); it is not a group of `vep`'s flags, so the landing
/// page's table does not list it.
const FILTER_PAGE: &str = "docs/cli/filter.md";

const FILTER_INTRO: &str = "`vep filter` filters an annotated output line by line, as Ensembl VEP's \
`filter_vep` does, over a VCF or tab output of `vep` or of Ensembl VEP. It takes the same options, \
read as Getopt::Long reads them: after one dash or two, in any case, and abbreviated to any prefix \
that names one option. A repeated option keeps its last value, and every `--filter` given is kept. \
The short forms are `-i`, `-o`, `-f`, `-s`, `-c`, `-l` and `-y`.\n\n\
```bash\n\
# The lines with a moderate-impact missense entry\n\
vep filter -i annotated.vcf \\\n  --filter \"Consequence is missense_variant and IMPACT is MODERATE\"\n\n\
# Only the passing entries, kept in the CSQ field\n\
vep filter -i annotated.vcf --only_matched -f \"SYMBOL is BRCA2\" -o brca2.vcf\n\n\
# The field names the input's header block offers\n\
vep filter -i annotated.tab --list\n\
```\n\n\
A filter is `<field> <operator> <value>`. Several filters are joined by `and` and `or`, evaluated \
left to right with no precedence, grouped by parentheses and inverted by `not`; repeated `--filter` \
flags are joined by `and`. The operators are `is` (`eq`, `=`), `ne` (`!=`), `gt` (`>`), `lt` (`<`), \
`gte` (`>=`), `lte` (`<=`), `exists` (`ex`, `defined`), `nex`, `match` (`matches`, `regex`, `re`), \
`nre`, `in` and `is_child`. `match` is a case-insensitive, unanchored regular expression; `in` takes \
a comma-separated list or a file with one value per line.\n\n\
A field alone tests that it is defined, `ne` without a value is `nex`, and a value starting with \
`#` names another field. A value shaped `text(number)` or `text:number` compares through its \
number against a numeric value, and against an `in` list when the list's first written entry is \
a number, where `filter_vep` picks a hash-order entry \
([the class](../divergences/transcript-filter-in-list-numeric-by-hash-order.md)). Write an `in` \
list's entries in one shape, `HGNC:2559,HGNC:1132` or `2559,1132`, and both tools keep the same \
lines. A value holding several entries joined by `,` or `&` after a digit passes when any entry \
passes.\n\n\
A field name resolves to a key of the line exactly, then in lower \
case, upper case, with a leading `_`, then as the unique case-insensitive prefix of one key; a \
prefix of several keys ends the run with an error.\n\n\
A VCF line's fields are its columns, its INFO keys and, for each entry of the consequence field, \
that entry's fields; an INFO flag without `=` is defined and empty. The line passes when any entry \
passes. Every INFO item under the field's name is read, so an output written with \
[`--keep_csq`](input-output.md#output) is filtered on its earlier annotation and its new one, and \
`--only_matched` keeps each item's passing entries. A tab line's fields are its columns, with an \
`Extra` column exploded into its pairs; a word without `=` among them is undefined. An \
empty value is undefined, and so is `-` under every field but `Allele`; an INFO value on a line \
without a consequence entry stands as written.\n\n\
The header block is written through, and a line passes unchanged unless `--only_matched` or \
`--soft_filter` rewrites it. The input is read as bytes and every line, header or data, is written \
as the bytes it was read from, so a line that is not UTF-8 is filtered as `filter_vep` filters it, \
each byte outside UTF-8 one character in a comparison.\n\n\
A command-line word is read as bytes too, so a value, a field name or a file name outside UTF-8 \
compares byte for byte.\n\n\
`match` runs its pattern over the bytes of the value as Perl runs it over a byte string: `.` is one \
byte, and `\\w`, `\\s`, `\\d` and the case fold are ASCII. A `\\p{}` class, a `\\x{}`, `\\o{}` or octal \
escape above 255 or a `\\N{U+}` escape switches the pattern to Unicode rules over the bytes read as \
Latin-1, and a `u`, `a` or `aa` flag switches its own group, as in Perl. Under those rules the case \
fold is full: `ss` matches `ß` and `\\x{DF}` matches `ss`.\n\n\
`--list` \
sorts the names by byte. An error naming several keys of a line lists them as the line holds them, \
the columns, the INFO keys, then the consequence fields; `filter_vep` lists them in its hash's \
order, which changes from run to run.\n\n\
`--ontology` rewrites `Consequence is <term>` to `Consequence is_child <term>`, which passes an \
entry whose term is the parent or lies below it in the Sequence Ontology. The answer comes from \
the built-in table of consequence terms rather than an Ensembl database, and the connection flags \
are refused by name. Where `vep filter` reads a line differently from `filter_vep` is stated on \
[Abbreviated fields, `is_child` and the regex \
dialect](../divergences/transcript-filter-fields.md#vep-filter).\n\n\
Every error is one line `ERROR: [filter] <message>` on standard error and exit 1; a command-line \
syntax error is exit 2.\n\n\
A line that standard error refuses, on a full device or a pipe without a reader, is dropped and \
changes no exit status. A run that ends with the warning of lines without the consequence field \
exits 0 with its output whole, and a refusal exits 1. `filter_vep` completes on a full device \
too, a refusal still ending non-zero; on a pipe without a reader it is killed by `SIGPIPE` (exit \
141) at its warning or its error, with none of its output written.\n\n";

/// The `vep filter` page: its introduction and one flag table per `--help`
/// heading of its command line.
fn render_filter() -> String {
    render_tool("vep filter", FILTER_INTRO, FilterArgs::command())
}

/// The page of the companion tool `vep recode`, rendered from its own command
/// line ([`RecodeArgs`]), like the `vep filter` page.
const RECODE_PAGE: &str = "docs/cli/recode.md";

const RECODE_INTRO: &str = "`vep recode` translates a variant between its identifiers and notations, as \
Ensembl VEP's `variant_recoder` does, from the cache alone. The input is an identifier, an HGVS, SPDI, \
VCF, region or Ensembl-format record, one per line, from `--input_file` or `--input_data` (`--id`). It \
takes the same options, read as Getopt::Long reads them: after one dash or two, in any case, and \
abbreviated to any prefix that names one option. A repeated option keeps its last value, and every \
`--transcript_filter` given is kept. The short forms are `-i`, `-s`, `-a` and `-u`.\n\n\
`--config` names a file of `<name> <value>` lines. An option is read from the four sources `vep` \
reads its own from, the first that names it winning: the command line, the `--config` file, \
`<dir>/vep.ini`, where `<dir>` is the file's `dir` or else `~/.vep`, then the environment variable \
`VEP_<OPTION>`. A repeatable option collects every source's values in that order, and a name no \
option of this page has is dropped.\n\n\
The options `variant_recoder` settles before it reads the sources are set on the command line \
alone, as there. `fields`, `buffer_size`, `failed`, `merged`, `lrg` and `ambiguous_hgvs` keep the \
values it presets. \
`vcf_string`, `var_synonyms`, `mane_select` and `ga4gh_vrs` add their output fields before the \
sources are read, so a source's value adds none. `pretty` and `vcf_string` are names its reader \
does not know. Every other option of this page but `--config` may be set from a source. A \
`vep` option this page lacks is dropped from every source, where `variant_recoder` reads it \
and applies it ([the class](../divergences/recode-config-file.md)).\n\n\
```bash\n\
# Every notation of an identifier\n\
vep recode --id rs699 --json_cache /path/to/cache\n\n\
# A protein notation, with its VCF form and its MANE Select transcript\n\
vep recode --id \"ENSP00000284967.6:p.Ala331Thr\" --vcf_string --mane_select \\\n  --pretty --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38\n\n\
# Two fields alone, from a file of SPDI notations\n\
vep recode -i variants.spdi.txt --fields hgvsg,hgvsc --json_cache /path/to/cache\n\
```\n\n\
The output is one JSON array on standard output, one object per input line; `--pretty` indents it. \
Each object is keyed by the input's alleles. An allele's object holds `input` and, for each field \
named, the unique values found in that allele's consequence entries. The default fields are `id`, \
`hgvsg`, `hgvsc`, `hgvsp` and `spdi`; `--fields` may name any key of the JSON output instead \
(`transcript_id`, `consequence_terms`).\n\n\
`--vcf_string` adds the allele's VCF form, `chromosome-position-reference-alternate`. `--var_synonyms` \
adds the known synonyms of the variant, one `source: id,id` string per source. `--mane_select` adds \
the `hgvsg`, `hgvsc` and `hgvsp` triple of each MANE Select transcript's entry and switches those \
three fields on. `--ga4gh_vrs` adds a GA4GH VRS Allele object for each `spdi`, which it switches \
on.\n\n\
The translation is an annotation run under a fixed configuration: `--json --hgvsg_use_accession \
--distance 0`, plus the HGVS and SPDI flags the fields need. HGVSg and SPDI name a chromosome by \
its `NC_` accession, and a scaffold or patch by its RefSeq name, where the cache's synonyms have \
one. The VCF form is built from the record's \
alleles for every kind of input, the base before an indel read from the reference, as `to_VCF_record` \
builds it. An input whose record has no consequence entry is `null`.\n\n\
The identifiers are the co-located variants' names per allele, a COSMIC or HGMD identifier listed \
under every allele. A known variant is co-located as the recoder's database query co-locates it: \
every row of the record's region whose span overlaps the record and whose alleles match, however the \
two are spelled, the variants failed in dbSNP included (`--failed 1`). So a deletion written from the \
last base of a repeat finds the row dbSNP writes from the first base, and an insertion finds a row \
only when it spans both of the bases around it.\n\n\
For an identifier input whose stored alleles are an unminimised indel, such as `GGG/GG`, \
`variant_recoder` writes the transcript and protein notations and the MANE Select triple without \
shifting the indel 3' along the transcript. It shifts them for every other input of the same \
variant. The stored class of the database row, `indel`, outlives the minimisation of its alleles, \
and `_return_3prime` shifts an insertion or a deletion alone.\n\n\
`vep recode` shifts them for every input, so those notations of such an identifier differ from the \
recoder's by the length of the repeat. The comparator's class `identifier_indel_transcript_notations` \
leaves them out of a comparison.\n\n\
In `variant_recoder`, a multi-allelic known variant's identifier reaches an allele of the record that \
spells none of its alternates only through the minimised form of its last alternate, so a deletion \
allele written on its own gets no identifier from it. `vep recode` minimises each alternate and names \
the variant on the allele it spells.\n\n\
The warnings raised while an input is read or annotated are carried in the output, under `warnings` \
on the next object written. `--buffer_size`, 1 by default, is the number of inputs read before an \
object is written, as in `variant_recoder`. A run that writes no object is one object of its warnings \
alone, and `--warning_file` writes them to a file as well. A `--warning_file` that cannot be \
written is refused before the input is annotated; a pipe or a device is first met by the write, \
once the input is annotated.\n\n\
The synonyms and the co-located variants are the cache's variation data, where `variant_recoder` \
reads the Ensembl database. The cache's synonyms are the ClinVar, OMIM, UniProt, PharmGKB, PhenCode \
and COSMIC identifiers and never the merged dbSNP identifiers (`Archive dbSNP`) the database adds; \
the comparator's class `archive_dbsnp_synonyms` leaves those out. A run that writes `var_synonyms` \
says so in its `INFO:` line. `variant_recoder --var_synonyms` \
stops with an error on a VCF input and on the identifier of a failed variant, where `vep recode` \
writes the cache's synonyms of the variant the ID column or the identifier names.\n\n\
A RefSeq notation (`NM_`, `NP_`, `XM_`, `XP_`) is written from a merged cache, and no variant is \
mapped to an LRG, as `--lrg` is refused. Every run that completes says the latter in one `INFO:` \
line on standard error, which no flag withholds, and the former where its cache is read as the \
Ensembl set. An input given on an LRG region the cache holds keeps its own notation.\n\n\
A protein notation resolves to one variant, as `--ambiguous_hgvs` is refused. A ClinGen allele \
identifier (`CA...`) is not an input. `ga4gh_vrs` holds a GA4GH VRS Allele object only for a \
chromosome whose synonyms name an `NC_` accession, and is an empty list for any other. \
The database connection flags (`--host`, `--port`, `--user`, `--password`, `--registry`, \
`--db_version`, `--genomes`, `--is_multispecies`) and `--bam` are refused with an error naming the \
flag. `--mane_select` is for human alone.\n\n\
`variant_recoder` reads the merged transcript set on every run, and `vep recode` reads it wherever \
the cache named holds it, with `--merged` or without. Such a cache is a `--json_cache` directory \
whose `info.json` declares the merged set, or a JSON cache in the `_merged` species directory under \
the cache root. Beside the Ensembl set's cache the run would read, it is the cache of the same \
`<cache_version>_<assembly>` directory name, so the run keeps its cache version and assembly; under \
a root with no directory of the Ensembl set, it is the one `--merged` resolves.\n\n\
An Ensembl VEP Storable cache is not one, and a cache of the Ensembl set alone is read as it is. \
Where the directory of the Ensembl set is there and holds no cache the run can read, the run is \
refused as it is without a merged cache, and `--merged` reads the merged one.\n\n\
`--refseq` is refused with `Can't use --refseq and --merged together`, as `variant_recoder` refuses \
it, and with `--gencode_basic` the error names that flag instead. A cache whose `info.json` \
declares the RefSeq set is not an input: `--refseq` alone reads it, and a run that names it is \
refused under the `cache` component. A `--json_cache` directory of RefSeq models converted \
without the `source_type` key is checked against nothing and read as the Ensembl set, its RefSeq \
models annotated. A merged cache carries the RefSeq models beside the Ensembl ones.\n\n\
Every error is one line `ERROR: [<component>] <message>` on standard error and exit 1, the \
component naming what refused the run (`recode`, `config`, `input`, `output`, `cache`, `fasta` or \
`filter`). A command-line syntax error is exit 2, a value outside the list of `--format`, \
`--pick_order` or `--shift_hgvs` among them.\n\n\
A `--transcript_filter` that names `_source_cache` or `dbID` over a cache built from a GTF or \
GFF3, whose transcripts carry neither, draws one `INFO:` line for each key, as under `vep` \
([the two keys](../divergences/transcript-filter-fields.md#a-gtf-or-gff3-as-the-transcript-source)). \
A run that completes writes them before its own line, and no flag withholds them.\n\n\
The annotation run's notices of what it left out of what it was given to read are written the \
same way. A FASTA file passed over for want of its index draws one `INFO:` line \
([Reference FASTA](../cache-setup.md#reference-fasta)). One line names each chromosome on \
which an entry the run writes is without an `hgvsg` or an `hgvsc` that takes the FASTA's \
sequence to write \
([the notations there](../cache-setup.md#a-contig-the-fasta-lacks)). One more counts the \
insertion and deletion alleles whose `hgvsc`, written as a field or in a `mane_select` triple, \
is on a transcript edited away from the genome \
([such a model](../output/refseq-and-merged.md#hgvsc-on-an-edited-model)). A run that completes \
writes these before its own line as well, in the order `vep` prints them, and no flag withholds \
them.\n\n\
A line that standard error refuses, on a full device or a pipe without a reader, is dropped. A \
refusal exits 1 without its `ERROR:` line, and a run that completes exits 0 without its `INFO:` \
lines, all it writes there, its \
warnings going into the output. `variant_recoder` completes on such a stream too, writing nothing \
there; its refusal ends non-zero, with another exit status than on a readable stream or killed by \
`SIGPIPE` (exit 141).\n\n";

/// The `vep recode` page.
fn render_recode() -> String {
    render_tool("vep recode", RECODE_INTRO, RecodeArgs::command())
}

/// What an abbreviated option of a companion tool rests on, the paragraph
/// after the one that offers the abbreviation: [`spell_as_declared`] resolves
/// a prefix against the options the tool's parser declares and nothing else,
/// so an option of the same prefix leaves it resolving to neither, a name
/// written with two dashes, whatever its own, for clap to refuse, and an
/// option named as the prefix takes it.
const TOOL_ABBREVIATIONS: &str = "An abbreviation is resolved against the options of the release \
that reads it, so a later release that adds an option can give a prefix to that option or make it \
name two options, which is refused (exit 2). A script that must keep working \
writes each option's full name.\n\n";

/// A companion tool's page: its title, its introduction, with
/// [`TOOL_ABBREVIATIONS`] after its opening paragraph, and one flag table per
/// `--help` heading of its command line.
fn render_tool(title: &str, intro: &str, mut command: Command) -> String {
    command.build();
    let sections = flag_sections(&command);
    let mut out = String::new();
    writeln!(out, "# {title}\n").unwrap();
    let (opening, rest) = intro
        .split_once("\n\n")
        .unwrap_or_else(|| panic!("the introduction of {title} opens with a paragraph"));
    writeln!(out, "{opening}\n").unwrap();
    out.push_str(TOOL_ABBREVIATIONS);
    out.push_str(rest);
    for (heading, _) in &sections {
        writeln!(out, "## {heading}\n").unwrap();
        flag_table(&mut out, heading, &sections);
    }
    out
}

const REGENERATE: &str = "VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference";

/// One group page: its file under `docs/`, its title, the landing table's
/// summary of it, its hand-written introduction, the `--help` headings it
/// carries, each as a section of the page, and the hand-written text that
/// follows a section's flag table, keyed by the section's title.
struct Group {
    file: &'static str,
    title: &'static str,
    covers: &'static str,
    intro: &'static str,
    sections: &'static [(&'static str, &'static str)],
    trailers: &'static [(&'static str, &'static str)],
}

/// The text under "Species and assembly": how each flag that reads data
/// Ensembl ships for human alone behaves on a cache of a species without it,
/// beside what the Ensembl release 116 mouse cache (GRCm39) holds, read from its
/// `info.txt`, its chromosome 19 variation file and the chromosome 19 transcript
/// shard at 5,000,001-6,000,000 (1,222 transcripts).
const SPECIES_FLAGS: &str = "### Flags that read data Ensembl ships for human alone\n\n\
Ensembl ships SIFT and PolyPhen predictions, population frequencies, the MANE, APPRIS, TSL and \
CCDS attributes and the GENCODE sets for human, and some of them for mouse. A flag that reads \
such data behaves on a cache of a species that lacks it as Ensembl VEP behaves: refused for SIFT \
and PolyPhen, which `--everything` switches off, and, on a cache that declares variation \
columns, for `--af_gnomad`, and otherwise as the table states. \
A column flag declares its columns in the header and leaves every cell empty, a transcript-set \
flag drops every transcript, a frequency filter whose population is accepted keeps or drops every \
record and `--regulatory` writes no regulatory row.\n\n\
The table lists each flag, its behaviour on a cache without the data, and what the Ensembl \
release 116 mouse cache (GRCm39) holds, read from its metadata, its chromosome 19 variants and \
one chromosome 19 transcript shard of 1,222 transcripts.\n\n\
The chromosome 19 variant counts are recorded in \
[the cache's count record](../concordance-provenance/mouse-grcm39-chromosome-19-variation-counts.json).\n\n\
| Flag | On a cache without the data | Mouse at release 116 |\n\
| --- | --- | --- |\n\
| `--sift` | Refused before the input is read: `ERROR: [cache] SIFT not available`; under \
`--everything` or `--safe` switched off with `INFO: disabling SIFT`, so no SIFT column is written | \
SIFT is present (`sift b`): predictions print |\n\
| `--polyphen` | As `--sift`, naming PolyPhen | No PolyPhen: refused, or switched off under \
`--everything` or `--safe` |\n\
| `--af_gnomad` | Refused: `ERROR: [cache] gnomAD data is not available in this cache`, since \
the cache declares no `gnomADe` column; `--everything` skips the check. On a cache that declares \
no variation columns the flag is accepted, as Ensembl VEP accepts it: the gnomAD exome columns \
are declared and every cell is empty | No `gnomADe` column: refused |\n\
| `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af` | Columns declared, every cell \
empty: a frequency is written only from a co-located variant that carries it | No frequency \
column: every cell empty |\n\
| `--check_frequency`, `--filter_common` | No record has frequency data, so `--freq_filter \
exclude` keeps every record and `include` drops every record; `--freq_pop` accepts `1KG_ALL` \
alone, any other population is refused as not declared. On a cache that declares no variation \
columns every `--freq_pop` name is accepted, as Ensembl VEP checks the name only where the cache \
declares variation columns. On a cache that declares no variation columns the filter is inactive, \
as in Ensembl VEP: `--freq_filter include` and `exclude` both keep every record, and `FREQS` is \
empty on every row | The same |\n\
| `--pubmed`, `--var_synonyms` | Columns declared; a cell is filled only where the co-located \
variant carries a citation or a synonym | Both columns declared and filled where a variant carries \
one: on chromosome 19, 186 of 2,809,961 variants carry a citation and 19,293 a synonym |\n\
| `--mane`, `--mane_select` | Columns declared, every cell empty | No MANE attribute: every cell \
empty |\n\
| `--appris`, `--tsl`, `--ccds` | Columns declared, every cell empty | Present: of the shard's \
1,222 transcripts, 205 carry APPRIS, 163 a TSL and 117 a CCDS |\n\
| `--gene_phenotype`, `--uniprot`, `--domains` | Columns declared, every cell empty | Present: in \
the shard, 865 transcripts belong to a gene with a phenotype, 287 carry a UniParc accession and 548 \
a protein domain |\n\
| `--gencode_basic` | Every transcript is dropped as the cache loads and every record is \
`intergenic_variant`; the warnings file gets one line, `WARNING: no transcript of the cache \
carries the gencode_basic attribute; every transcript is filtered` | The attribute is present \
(GENCODE M39; 612 of the shard's 1,222 transcripts): the set applies |\n\
| `--gencode_primary`, `--flag_gencode_primary` | Ensembl VEP's warning that the option is only \
available for human on the GRCh38 assembly; `--gencode_primary` then drops every transcript, with \
the warning above naming `gencode_primary`, and `--flag_gencode_primary` writes an empty \
`GENCODE_PRIMARY` column | The attribute is present (136 of the shard's 1,222 transcripts): the \
assembly warning prints and the set or the column applies all the same |\n\
| `--regulatory` | No RegulatoryFeature or MotifFeature row; one warning says the cache holds no \
regulatory data | Regulatory data is present (83 cell types) |\n\n";

/// The text under "Filtering and picking": how an `in` list of
/// `--transcript_filter` that mixes a number and other text is read, and the
/// shape that keeps both engines on the same rows.
const TRANSCRIPT_FILTER_LISTS: &str =
    "An `in` list of `--transcript_filter` that mixes a number and \
other text, `gene_hgnc_id in HGNC:2559,1132`, reads a `text:number` or `text(number)` value \
through its number when its first written entry is a number and whole otherwise \
([the class](../divergences/transcript-filter-in-list-numeric-by-hash-order.md)). Ensembl VEP \
decides by an entry it picks from the list's hash, so write the entries in one shape, \
`HGNC:2559,HGNC:1132` or `2559,1132`, and both engines keep the same rows.\n\n";

const GROUPS: &[Group] = &[
    Group {
        file: "docs/cli/input-output.md",
        title: "Input and output",
        covers: "The input file, its format and how it is read, the output file, its format and \
                 compression, the header lines, the warnings and skipped-variants files.",
        intro: "These flags name the input file and how it is read, and the output file, its \
                format, its compression and its header lines. The formats themselves are described \
                under [Output formats](../output-formats.md).\n\n",
        sections: &[("Input", "Input"), ("Output", "Output")],
        trailers: &[],
    },
    Group {
        file: "docs/cli/cache-species-runtime.md",
        title: "Cache, species and runtime",
        covers: "The cache directory and its layout, the reference FASTA, the 3' shifting of \
                 insertions and deletions, species and assembly, the batch size, thread counts \
                 and distances.",
        intro: "These flags name the annotation cache and the reference FASTA, the species and \
                assembly, and the batch size, thread counts and distances of a run. The cache is \
                found as Ensembl VEP finds it, under `--dir_cache`, `--species` and \
                `--cache_version`, or named directly with `--json_cache`. Building or converting \
                a cache is described under [Cache setup](../cache-setup.md), and the genomes \
                beside human, with each cache's release, genetic codes and chromosome names, \
                under [Other species](../cache-species.md). The flags that read data Ensembl \
                ships for human alone are listed under [Species and \
                assembly](#flags-that-read-data-ensembl-ships-for-human-alone).\n\n\
                The flags under [3' shifting](#3-shifting) move an insertion or deletion to its \
                most 3' position on a transcript before it is annotated or described in HGVS, \
                and report how far it moved.\n\n",
        sections: &[
            ("Cache / Offline", "Cache and offline"),
            ("3' Shifting", "3' shifting"),
            ("Species / Assembly", "Species and assembly"),
            ("Runtime", "Runtime"),
        ],
        trailers: &[("Species and assembly", SPECIES_FLAGS)],
    },
    Group {
        file: "docs/cli/output-control.md",
        title: "Output control",
        covers: "The field flags, `--everything` among them, that choose which annotation \
                 columns appear, the flags that report co-located known variants, and the two \
                 flags that drop rows.",
        intro: "These flags choose which annotation fields appear in the output, report the \
                known variants co-located with each input variant, and drop the rows outside a \
                coding sequence or outside every transcript. `--everything` switches on exactly \
                Ensembl VEP's `--everything` set, listed under \
                [option sets](../cli-reference.md#option-sets-configpm351-593).\n\n",
        sections: &[
            ("Annotation Fields", "Annotation fields"),
            ("Co-located Variants", "Co-located variants"),
            ("Dropped Rows", "Dropped rows"),
        ],
        trailers: &[],
    },
    Group {
        file: "docs/cli/allele-frequencies.md",
        title: "Allele frequencies and co-located variants",
        covers: "Co-located known variants, population frequencies, the frequency filter, PubMed \
                 identifiers and gene-phenotype associations.",
        intro: "These flags report the population allele frequencies, PubMed identifiers and \
                gene-phenotype associations of co-located known variants, and filter the input by \
                frequency. `--check_existing`, which reports the co-located variants themselves, \
                is listed under [Output control](output-control.md#co-located-variants); the \
                frequency flags `--everything` switches on are listed under \
                [option sets](../cli-reference.md#option-sets-configpm351-593).\n\n",
        sections: &[("Allele Frequencies", "Flags")],
        trailers: &[],
    },
    Group {
        file: "docs/cli/filtering-picking.md",
        title: "Filtering, picking and transcript sets",
        covers: "The pick and flag-pick families, the transcript filter, the allele and position \
                 columns, the transcript-set flags and the reference allele a row reads.",
        intro: "These flags reduce the rows reported per record, allele or gene, flag the ones a \
                pick would drop and filter transcripts by their attributes. The flags under \
                [Alleles and positions](#alleles-and-positions) add the allele columns, trim the \
                alleles and add the total length to the positions. The transcript-set flags \
                select the cache read and the transcript set annotated: the GENCODE basic or \
                primary set, and on a RefSeq or merged cache the RefSeq accessions unless \
                `--all_refseq` keeps every model. The flags under \
                [Reference allele](#reference-allele) choose whether a row's reference allele is \
                the transcript's or the input's.\n\n",
        sections: &[
            ("Filtering / Picking", "Filtering and picking"),
            ("Alleles / Positions", "Alleles and positions"),
            ("Transcript Sets", "Transcript sets"),
            ("Reference Allele", "Reference allele"),
        ],
        trailers: &[("Filtering and picking", TRANSCRIPT_FILTER_LISTS)],
    },
    Group {
        file: "docs/cli/plugins.md",
        title: "Plugins",
        covers: "Loading built-in and dynamic-library plugins, and the custom annotation \
                 sources with their bigWig and conservation-score shortcuts.",
        intro: "These flags load the built-in plugins and dynamic-library plugins, and name the \
                directory the latter are searched in. The built-in plugins, their parameters and \
                their data files are described under [Plugins](../plugins.md); the C ABI a \
                dynamic-library plugin implements is described under \
                [Writing a plugin](../developers/writing-a-plugin.md). A `--plugin` naming an \
                Ensembl VEP plugin vep-rs does not provide ends the run with an error naming the \
                reason, with or without `--safe`, as does a Perl module (`Name.pm` by path, or \
                `Name` with `Name.pm` under the plugins directory). The \
                [support table](../plugins.md) lists every Ensembl VEP plugin and its status.\n\n\
                The flags under [Custom annotation](#custom-annotation) name the `--custom` \
                sources, their bigWig and conservation-score shortcuts, and how a VCF source's \
                values are written. The sources themselves are described under \
                [Custom annotation](../custom-annotation.md).\n\n",
        sections: &[
            ("Plugins", "Loading plugins"),
            ("Custom Annotation", "Custom annotation"),
        ],
        trailers: &[],
    },
    Group {
        file: "docs/cli/statistics-verbosity-misc.md",
        title: "Statistics, verbosity and record handling",
        covers: "The statistics flags, logging, and the two flags that keep a record a run would \
                 otherwise skip.",
        intro: "These flags cover the statistics report written beside the output, the logging \
                level, and the two flags that keep a record a run would otherwise skip.\n\n",
        sections: &[
            ("Statistics", "Statistics"),
            ("Verbosity", "Verbosity"),
            ("Record Handling", "Record handling"),
        ],
        trailers: &[],
    },
    Group {
        file: COMPATIBILITY,
        title: "Compatibility with Ensembl VEP",
        covers: "The database flags, the flags refused or accepted without effect, where \
                 `Config.pm` is not followed, the configuration sources, the option spellings \
                 and the exit status.",
        intro: "This page gathers the Ensembl VEP behaviour vep-rs refuses or accepts without \
                effect, and how the command line is read. It covers the database flags, the \
                error table, where `Config.pm` is not followed, the configuration sources and \
                their priority, the option spellings, and the exit status and error \
                contract.\n\n",
        sections: &[("Database", "Database flags")],
        trailers: &[],
    },
];

/// The declared flags whose help reads `Accepted for VEP compatibility`, with
/// what Ensembl VEP does with each. The rendering checks the list against the
/// help texts both ways, so a flag that gains an implementation leaves this
/// table and a new compatibility flag enters it; the table is empty while no
/// declared flag reads so.
const COMPATIBILITY_FLAGS: &[(&str, &str)] = &[];

const COMPATIBILITY_HELP: &str = "accepted for vep compatibility";

/// The repository path of a page, from this crate's directory.
fn repo_path(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(file)
}

/// A page's link target from the landing page: its path under `docs/`.
fn site_link(file: &str) -> &str {
    file.strip_prefix("docs/")
        .unwrap_or_else(|| panic!("{file} is not a page under docs/"))
}

/// A Markdown table cell: newlines and runs of blanks collapsed, a pipe escaped
/// everywhere, and outside a code span a bare `--flag` put in code (the site's
/// smart punctuation would turn its dashes into an en dash) and `<`, `>` and
/// `_` escaped, so that `<leaf>` is not read as an HTML tag and
/// `[_refseq|_merged]` does not open an emphasis.
fn cell(text: &str) -> String {
    let joined = text.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut out = String::with_capacity(joined.len() + 16);
    for (i, segment) in joined.split('`').enumerate() {
        if i > 0 {
            out.push('`');
        }
        let in_code_span = i % 2 == 1;
        if in_code_span {
            out.push_str(&segment.replace('|', "\\|"));
        } else {
            out.push_str(&prose_cell(segment));
        }
    }
    out
}

/// Prose outside a code span: each bare `--flag` (with an attached `=value`)
/// in code, the rest escaped.
fn prose_cell(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(pos) = rest.find("--") {
        let (before, from) = rest.split_at(pos);
        out.push_str(&escape_prose(before));
        let at_word_start = before.chars().last().is_none_or(|c| !c.is_alphanumeric());
        let name_len = flag_token_len(&from[2..]);
        if at_word_start && name_len > 0 {
            let mut end = 2 + name_len;
            if from[end..].starts_with('=') {
                end += 1 + flag_token_len(&from[end + 1..]);
            }
            out.push_str(&code(&from[..end]));
            rest = &from[end..];
        } else {
            out.push_str("--");
            rest = &from[2..];
        }
    }
    out.push_str(&escape_prose(rest));
    out
}

/// The length in bytes of the flag-shaped word that opens `text`: ASCII
/// letters, digits and underscores, with a hyphen counted only between two of
/// them, so `--core-dumps` is one word and a trailing dash is not.
fn flag_token_len(text: &str) -> usize {
    let bytes = text.as_bytes();
    let is_word_byte = |byte: u8| byte.is_ascii_alphanumeric() || byte == b'_';
    let mut len = 0;
    while len < bytes.len() {
        let byte = bytes[len];
        let inner_hyphen =
            byte == b'-' && len > 0 && bytes.get(len + 1).is_some_and(|next| is_word_byte(*next));
        if !(is_word_byte(byte) || inner_hyphen) {
            break;
        }
        len += 1;
    }
    len
}

fn escape_prose(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        if matches!(c, '|' | '<' | '>' | '_') {
            out.push('\\');
        }
        out.push(c);
    }
    out
}

fn code(text: &str) -> String {
    format!("`{text}`")
}

/// A long flag name in code: `--name`.
fn flag_in_code(name: &str) -> String {
    code(&format!("--{name}"))
}

fn flags_in_code(flags: &[&str]) -> String {
    flags
        .iter()
        .map(|f| flag_in_code(f))
        .collect::<Vec<_>>()
        .join(", ")
}

/// `flags`, a list of Ensembl VEP's option set (`Config.pm:327-346`), split into
/// those the command line declares, which take their values as the list says,
/// and those it does not, which are refused as options this version lacks.
fn declared_and_refused(
    command: &Command,
    flags: &[&'static str],
) -> (Vec<&'static str>, Vec<&'static str>) {
    let declared: Vec<&str> = command.get_arguments().filter_map(Arg::get_long).collect();
    flags
        .iter()
        .copied()
        .partition(|flag| declared.contains(flag))
}

/// The refused flags of a list named with what the parser does with them,
/// `shape` being the list's property in Ensembl VEP (`repeatable`,
/// `comma-separated`), or nothing when every flag of the list is declared. A
/// flag `@VEP_PARAMS` declares is refused on the command line and in a
/// configuration file; one it does not declare (`dont_export`, an option of the
/// `haplo` script) is unknown on the command line and dropped from a
/// configuration file, as Ensembl VEP's `is_valid_param` drops it
/// (`Config.pm:1046-1047`).
fn refused_in_list(refused: &[&str], shape: &str) -> String {
    let perl = perl_options();
    let (lacked, undeclared): (Vec<&str>, Vec<&str>) = refused
        .iter()
        .copied()
        .partition(|flag| perl.iter().any(|option| option.canonical == *flag));
    let mut out = String::new();
    if let Some((names, verb)) = named_options(&lacked) {
        write!(
            out,
            " {names}, {shape} in Ensembl VEP, {verb} this version lacks, refused on the command \
             line and in a configuration file (exit 2)."
        )
        .unwrap();
    }
    if let Some((names, verb)) = named_options(&undeclared) {
        write!(
            out,
            " {names}, {shape} in Ensembl VEP's list, {verb} neither its `vep` script \
             (`Config.pm:73-276`) nor this version declares: unknown on the command line (exit \
             2) and dropped from a configuration file, as in Ensembl VEP."
        )
        .unwrap();
    }
    out
}

/// `flags` in code as prose (`` `--a` and `--b` ``) with the verb that agrees
/// with their number, or nothing for no flag.
fn named_options(flags: &[&str]) -> Option<(String, &'static str)> {
    let (last, rest) = flags.split_last()?;
    Some(if rest.is_empty() {
        (flag_in_code(last), "is an option")
    } else {
        (
            format!("{} and {}", flags_in_code(rest), flag_in_code(last)),
            "are options",
        )
    })
}

/// What follows the names of a flag `--help` leaves out.
const HIDDEN_MARK: &str = " (hidden)";

/// What the mark means, the close of the sentences a page of `vep`'s flags
/// carries under its first table with a marked row. The parser declares a
/// hidden flag like any other, so the pre-pass reads it from the command line,
/// from a configuration file's key and from its `VEP_<FLAG>` variable.
const HIDDEN_NOTE: &str = "are left out of `--help` and are kept for tuning and diagnosis. The \
name, the values, the default and the effect of such a flag may change, and the flag may be \
removed, in any release, however it is given: on the command line, in a configuration file or in \
a `VEP_<FLAG>` variable ([configuration sources](compatibility.md#configuration-sources)).\n\n";

/// The sentences on the mark for a page with `marked` marked rows.
fn hidden_note(marked: usize) -> String {
    format!("Flags marked `(hidden)`, {marked} on this page, {HIDDEN_NOTE}")
}

/// The rows of a `--help` heading's table that carry the mark.
fn marked_rows(heading: &str, sections: &[FlagSection]) -> usize {
    sections
        .iter()
        .filter(|(h, _)| h == heading)
        .flat_map(|(_, rows)| rows)
        .filter(|row| row[0].ends_with(HIDDEN_MARK))
        .count()
}

/// `--long`, then `-s`, then the visible aliases, each in code.
fn flag_names(arg: &Arg) -> String {
    let mut names = Vec::new();
    if let Some(long) = arg.get_long() {
        names.push(flag_in_code(long));
    }
    if let Some(short) = arg.get_short() {
        names.push(code(&format!("-{short}")));
    }
    for alias in arg.get_visible_aliases().into_iter().flatten() {
        names.push(flag_in_code(alias));
    }
    for alias in arg.get_visible_short_aliases().into_iter().flatten() {
        names.push(code(&format!("-{alias}")));
    }
    let mut out = names.join(", ");
    if arg.is_hide_set() {
        out.push_str(HIDDEN_MARK);
    }
    out
}

/// The value type as the reference names it: `bool` for a switch, the possible
/// values when the flag is an enumeration, the Rust type otherwise, with
/// `repeatable` for a flag that accumulates and `comma-separated` for a
/// `LIST_FLAGS` entry, whose one value is split on commas.
fn flag_type(arg: &Arg) -> String {
    if matches!(arg.get_action(), ArgAction::SetTrue) {
        return "bool".to_string();
    }
    let possible = arg.get_possible_values();
    let base = if !possible.is_empty() {
        possible
            .iter()
            .map(|p| code(p.get_name()))
            .collect::<Vec<_>>()
            .join(", ")
    } else {
        let id = arg.get_value_parser().type_id();
        let name = if id == TypeId::of::<String>() {
            "String"
        } else if id == TypeId::of::<usize>() {
            "usize"
        } else if id == TypeId::of::<u32>() {
            "u32"
        } else if id == TypeId::of::<u64>() {
            "u64"
        } else if id == TypeId::of::<i64>() {
            "i64"
        } else if id == TypeId::of::<f64>() {
            "f64"
        } else {
            panic!(
                "--{}: a value type the reference does not name",
                arg.get_id()
            )
        };
        name.to_string()
    };
    let is_list_flag = arg
        .get_long()
        .is_some_and(|long| LIST_FLAGS.contains(&long));
    // `--cell_type` is a list whose names may hold a comma: read by longest
    // match against the cache's names, never split on every comma.
    let qualifier = if arg.get_long() == Some("cell_type") {
        Some("list of the cache's epigenome names")
    } else if is_list_flag {
        Some("comma-separated")
    } else if matches!(arg.get_action(), ArgAction::Append) {
        Some("repeatable")
    } else {
        None
    };
    match qualifier {
        Some(qualifier) => format!("{base} ({qualifier})"),
        None => base,
    }
}

/// `false` for a switch, `none` for no default, a number as written, an empty
/// default as `""`, any other default in code.
fn flag_default(arg: &Arg) -> String {
    if matches!(arg.get_action(), ArgAction::SetTrue) {
        return "false".to_string();
    }
    let defaults: Vec<String> = arg
        .get_default_values()
        .iter()
        .map(|v| {
            let value = v.to_string_lossy();
            if value.is_empty() {
                code("\"\"")
            } else if value.parse::<f64>().is_ok() {
                value.into_owned()
            } else {
                code(&value)
            }
        })
        .collect();
    if defaults.is_empty() {
        "none".to_string()
    } else {
        defaults.join(", ")
    }
}

fn flag_help(arg: &Arg) -> String {
    arg.get_long_help()
        .or(arg.get_help())
        .map(|h| h.to_string())
        .unwrap_or_default()
}

fn table(out: &mut String, header: &[&str], rows: &[Vec<String>]) {
    writeln!(out, "| {} |", header.join(" | ")).unwrap();
    writeln!(
        out,
        "|{}",
        header.iter().map(|_| " --- |").collect::<String>()
    )
    .unwrap();
    for row in rows {
        writeln!(out, "| {} |", row.join(" | ")).unwrap();
    }
    out.push('\n');
}

/// One `--help` heading and the table row of each flag under it.
type FlagSection = (String, Vec<Vec<String>>);

/// The flag rows of every `--help` heading, in the command's order.
fn flag_sections(command: &Command) -> Vec<FlagSection> {
    let mut sections: Vec<FlagSection> = Vec::new();
    for arg in command.get_arguments() {
        let Some(long) = arg.get_long() else {
            continue;
        };
        // clap's own `--help` and `--version` carry no heading; a declared
        // `--version` (the database release of `vep filter`) is a flag like
        // any other.
        if long == "help" || (long == "version" && arg.get_help_heading().is_none()) {
            continue;
        }
        let heading = arg
            .get_help_heading()
            .unwrap_or_else(|| panic!("--{long} has no help heading"))
            .to_string();
        let row = vec![
            flag_names(arg),
            flag_type(arg),
            flag_default(arg),
            cell(&flag_help(arg)),
        ];
        match sections.iter_mut().find(|(h, _)| *h == heading) {
            Some((_, rows)) => rows.push(row),
            None => sections.push((heading, vec![row])),
        }
    }
    sections
}

/// The flag table of one `--help` heading. A heading `GROUPS` names that no
/// flag carries is a stale `GROUPS`.
fn flag_table(out: &mut String, heading: &str, sections: &[FlagSection]) {
    let rows = sections
        .iter()
        .find(|(h, _)| h == heading)
        .map(|(_, rows)| rows)
        .unwrap_or_else(|| {
            panic!(
                "the --help heading {heading:?} named in GROUPS has no flags: a heading renamed \
                 in args.rs is renamed in GROUPS too"
            )
        });
    table(out, &["Flag", "Type", "Default", "Description"], rows);
}

/// A group page: its title, its introduction and one flag table per `--help`
/// heading it carries, with [`hidden_note`] under the first table that marks
/// a flag.
fn render_group(group: &Group, sections: &[FlagSection]) -> String {
    let mut out = String::new();
    writeln!(out, "# {}\n", group.title).unwrap();
    out.push_str(group.intro);
    let marked: usize = group
        .sections
        .iter()
        .map(|(heading, _)| marked_rows(heading, sections))
        .sum();
    let mut noted = false;
    for (heading, section) in group.sections {
        writeln!(out, "## {section}\n").unwrap();
        flag_table(&mut out, heading, sections);
        if !noted && marked_rows(heading, sections) > 0 {
            out.push_str(&hidden_note(marked));
            noted = true;
        }
        if let Some((_, trailer)) = group.trailers.iter().find(|(title, _)| title == section) {
            out.push_str(trailer);
        }
    }
    out
}

/// The landing page's note under the option-set table: the table prints
/// [`OPTION_SETS`], Ensembl VEP's text, whose `--phyloP` row spells the
/// source's short name as `Config.pm:555` does, while the run names the source
/// as its file is named.
const PHYLOP_ROW_NOTE: &str = "The `--phyloP` row is Ensembl VEP's own text, the short name \
`phlyoP%phyloP%way` included; vep-rs names that source `phyloP<N>way` \
([The key of a `--phyloP` source](divergences/phylop-short-name.md)).\n\n";

/// The landing page: what the command line is, the group pages, how a flag is
/// read, and Ensembl VEP's option rules applied as data.
fn render_landing(command: &Command) -> String {
    let mut out = String::new();
    out.push_str("# CLI reference\n\n");
    out.push_str(
        "The `vep` binary accepts the Ensembl VEP command line. Every flag but the few vep-rs \
         adds is Ensembl VEP's and does what Ensembl VEP's does, diverges from it, is accepted \
         without effect, or is refused with an error naming what it would do. This reference and \
         [Intended divergences](intended-divergences.md) document such divergences; the flags \
         vep-rs adds are described beside Ensembl VEP's. A flag the error table governs (refused, \
         or accepted without effect) carries that table's text as its help.\n\n",
    );
    writeln!(
        out,
        "These pages are generated from the command-line definition and the option tables by \
         the `cli_reference` test of `crates/vep-cli`. Edit the sources, then run `{REGENERATE}` \
         to rewrite them.\n"
    )
    .unwrap();
    out.push_str("```bash\n# Show all flags\nvep --help\n\n# Short form\nvep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache\n```\n\n");
    out.push_str(
        "Worked examples are on [Examples](examples.md); the output formats and their parity \
         with Ensembl VEP on [Output formats](output-formats.md); the divergences from Ensembl \
         VEP that a flag's description names on [Intended divergences](intended-divergences.md).\n\n",
    );

    out.push_str("## Option groups\n\n");
    out.push_str(
        "Every flag has a row on one of these pages, grouped by the `--help` heading the command \
         line declares for it: its names, its value type, its default and its description. `--help` \
         itself leaves out the flags marked `(hidden)`. A default of `false` is a switch that \
         is off; `none` marks a flag with no declared default, and where a value applies in its \
         absence the description names it.\n\n",
    );
    let rows: Vec<Vec<String>> = GROUPS
        .iter()
        .map(|g| {
            vec![
                format!("[{}]({})", g.title, site_link(g.file)),
                cell(g.covers),
            ]
        })
        .collect();
    table(&mut out, &["Group", "Covers"], &rows);
    writeln!(
        out,
        "A repeated single-value flag keeps its last value, as Getopt::Long does. `-o` \
         accumulates and the last wins with a warning, as the `vep` script does (vep:32-33). A \
         long option is also read in any case, after one dash and as a unique prefix, which is \
         resolved against the flags of the release that reads it ([option \
         spellings]({0}#option-spellings-configpm62)), and from the `--config` file, `vep.ini` \
         and the `VEP_*` variables ([configuration sources]({0}#configuration-sources)).\n",
        site_link(COMPATIBILITY)
    )
    .unwrap();
    out.push_str(
        "`--help` prints the flag list and exits 0, as does a command line with no argument at \
         all (vep:29); `--version` prints `vep <version>` and exits 0.\n\n",
    );

    out.push_str("## Option rules\n\n");
    out.push_str(
        "Ensembl VEP's option rules are applied as data, from `Config.pm:327-645` at \
         release/116.3, in this order after parsing. A flag the error table refuses is refused \
         first. Then the option sets switch on the flags they imply, in the order below, a set's \
         value overwriting an explicit one (`apply_option_sets`, `Config.pm:816-884`). Then the \
         checks of `check_config` (`Config.pm:944-985`) run: valid values, incompatible pairs, \
         requirements, deprecations, each refused with Ensembl VEP's own message. A rule naming \
         a flag vep-rs does not declare cannot fire and is listed so the table reads as Ensembl \
         VEP's.\n\n",
    );
    out.push_str("### Option sets (`Config.pm:351-593`)\n\n");
    let rows: Vec<Vec<String>> = OPTION_SETS
        .iter()
        .map(|set| {
            vec![
                flags_in_code(set.flags),
                set.set
                    .iter()
                    .map(|(k, v)| {
                        if *v == "1" {
                            flag_in_code(k)
                        } else {
                            code(&format!("--{k} {v}"))
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(", "),
            ]
        })
        .collect();
    table(&mut out, &["Flag", "Switches on"], &rows);
    out.push_str(
        "`--everything` therefore reaches `--hgvsc` and `--hgvsp` through `--hgvs`, \
         `--mane_select` and `--mane_plus_clinical` through `--mane`, and `--check_existing` \
         through the frequency flags; it does not switch on `--allele_number` or `--hgvsg`. \
         `--regulatory` is switched on and emits RegulatoryFeature and MotifFeature rows when \
         the cache carries regulatory features, with every motif field where it carries the \
         binding matrices' counts; `CELL_TYPE` needs `--cell_type`, which is not in Ensembl \
         VEP's set either.\n\n",
    );
    out.push_str(PHYLOP_ROW_NOTE);
    out.push_str("### Incompatible flags (`Config.pm:615-632`)\n\n");
    out.push_str(
        "`Can't use --<flag> and --<other> together`, naming the first pair in this table's \
         order (Ensembl VEP iterates a hash there, so which pair it names varies from run to \
         run); Ensembl VEP skips the check under `--safe`.\n\n",
    );
    let rows: Vec<Vec<String>> = INCOMPATIBLE
        .iter()
        .map(|(flag, others)| vec![flag_in_code(flag), flags_in_code(others)])
        .collect();
    table(&mut out, &["Flag", "Cannot be combined with"], &rows);
    if !INCOMPATIBLE_LIFTED.is_empty() {
        out.push_str("Pairs Ensembl VEP refuses that are accepted here:\n\n");
        let rows: Vec<Vec<String>> = INCOMPATIBLE_LIFTED
            .iter()
            .map(|pair| {
                vec![
                    flag_in_code(pair.flag),
                    flag_in_code(pair.other),
                    cell(pair.reason),
                ]
            })
            .collect();
        table(
            &mut out,
            &["Flag", "Beside", "Why the pair is accepted"],
            &rows,
        );
    }
    out.push_str("### Required flags (`Config.pm:607-612`)\n\n");
    out.push_str("`You must set --<required> to use --<flag>`.\n\n");
    let rows: Vec<Vec<String>> = REQUIRES
        .iter()
        .map(|(flag, required)| vec![flag_in_code(flag), flags_in_code(required)])
        .collect();
    table(&mut out, &["Flag", "Needs"], &rows);
    out.push_str("### Valid values (`Config.pm:596-604`)\n\n");
    out.push_str(
        "A declared flag carries its list as a value enumeration, so a value outside it is an \
         argument syntax error (exit 2) listing the possible values.\n\n",
    );
    let rows: Vec<Vec<String>> = VALID
        .iter()
        .map(|(flag, values)| {
            vec![
                flag_in_code(flag),
                values
                    .iter()
                    .map(|v| code(v))
                    .collect::<Vec<_>>()
                    .join(", "),
            ]
        })
        .collect();
    table(&mut out, &["Flag", "Values"], &rows);
    out.push_str("### Deprecated flags (`Config.pm:636-638`)\n\n");
    for (flag, replacement) in DEPRECATED {
        match replacement {
            Some(new) => writeln!(
                out,
                "- `--{flag}`: `--{flag} has been deprecated - please use --{new} instead`"
            )
            .unwrap(),
            None => writeln!(out, "- `--{flag}`: `--{flag} has been deprecated`").unwrap(),
        }
    }
    out.push('\n');
    out.push_str("### Lists (`Config.pm:327-346`)\n\n");
    let (split, refused) = declared_and_refused(command, LIST_FLAGS);
    writeln!(
        out,
        "- Comma-separated values: {}.{}{}",
        flags_in_code(&split),
        if split.contains(&"cell_type") {
            CELL_TYPE_LIST_NOTE
        } else {
            ""
        },
        refused_in_list(&refused, "comma-separated")
    )
    .unwrap();
    let (collected, refused) = declared_and_refused(command, ALLOW_MULTIPLE);
    writeln!(
        out,
        "- Repeatable, each occurrence appending: {}.{}",
        flags_in_code(&collected),
        refused_in_list(&refused, "repeatable")
    )
    .unwrap();
    writeln!(
        out,
        "- Input formats Ensembl VEP exempts from its sorted-input check (`Config.pm:640-645`): \
         {}. vep-rs performs no such check in any format: every record is annotated in input \
         order, and an input whose chromosomes interleave keeps more transcript indexes in \
         memory.\n",
        UNSORTABLE
            .iter()
            .map(|f| code(f))
            .collect::<Vec<_>>()
            .join(", ")
    )
    .unwrap();
    out.push_str("### Rules of vep-rs's own\n\n");
    out.push_str(
        "Checked after Ensembl VEP's rules on the same state, for behaviour Ensembl VEP has no \
         flag for; each is refused with the line below.\n\n",
    );
    let rows: Vec<Vec<String>> = OWN_RULES
        .iter()
        .map(|rule| {
            vec![
                code(rule.component),
                code(&format!("ERROR: [{}] {}", rule.component, rule.message)),
            ]
        })
        .collect();
    table(&mut out, &["Component", "Error"], &rows);
    out
}

/// What an abbreviated flag of `vep` rests on, an item of the compatibility
/// page's list of option spellings: [`resolve_spelling`] resolves a prefix
/// against the names the parser declares, aliases and hidden flags among them,
/// and only then against Ensembl VEP's list. A flag named as the prefix takes
/// it, and a flag of the same prefix leaves it to Ensembl VEP's list, which
/// resolves it or leaves the token for clap to refuse.
const VEP_ABBREVIATIONS: &str = "- An abbreviation is resolved against the flags of the release \
that reads it, its aliases and its hidden flags included. A later release that adds a flag can \
give a prefix to that flag, or make a prefix that Ensembl VEP's list does not resolve name two \
flags, which is refused (exit 2); a script that must keep working writes each flag's full \
name.\n";

/// The compatibility page: the database flags, the error table, where
/// `Config.pm` is not followed, the configuration sources, the option
/// spellings and the exit status.
fn render_compatibility(group: &Group, command: &Command, sections: &[FlagSection]) -> String {
    let mut out = String::new();
    writeln!(out, "# {}\n", group.title).unwrap();
    out.push_str(group.intro);
    let [(heading, section)] = group.sections else {
        panic!("the compatibility page carries the database flags alone");
    };
    writeln!(out, "## {section}\n").unwrap();
    out.push_str(
        "vep-rs reads a local cache and never contacts a database. Each database flag is \
         declared so that the command line names it, and refused with the error its description \
         shows.\n\n",
    );
    flag_table(&mut out, heading, sections);
    let marked = marked_rows(heading, sections);
    if marked > 0 {
        out.push_str(&hidden_note(marked));
    }

    out.push_str("## Consequence terms never assigned\n\n");
    out.push_str(
        "One term is never assigned: the root `sequence_variant` has no predicate, and Ensembl \
         VEP's consequence engine never assigns it either. The other 40 of Ensembl VEP's 41 \
         Sequence Ontology terms are assigned; the six regulatory and motif terms need \
         `--regulatory` (or `--everything`) and a cache that carries regulatory features \
         ([Regulatory and motif features](../regulatory.md)).\n\n",
    );

    out.push_str("## Input formats\n\n");
    out.push_str(
        "The input formats `--format` reads, how `guess` detects one, and the formats refused \
         are stated in its row under [Input](input-output.md#input).\n\n",
    );

    out.push_str("## Refused flags\n\n");
    out.push_str(
        "Every flag vep-rs refuses, with the line it prints and what Ensembl VEP does with the \
         flag (cited to its source at release/116.3). A refused flag exits 1 with the message \
         below under the `cli` component, and the help text of a declared one is this table's \
         message. The `-V` row is the exception: `-V` is the version's short, and its help text \
         is the version's.\n\n",
    );
    let mut rows: Vec<Vec<String>> = ERROR_TABLE
        .iter()
        .filter_map(|entry| match entry.treatment {
            Treatment::Refused { component, message } => Some(vec![
                error_table_flag(entry.flag, entry.declared),
                code(&format!("ERROR: [{component}] {message}")),
                cell(entry.ensembl_vep),
            ]),
            Treatment::AcceptedWithoutEffect { .. } => None,
        })
        .collect();
    rows.push(vec![
        "`-V` beside other arguments (caught before parsing)".to_string(),
        code(&format!("ERROR: [cli] {VERSION_SHORT_BESIDE_ARGUMENTS}")),
        cell(
            "Config.pm:159 declares verbose|v, and Getopt::Long's `ignore_case`, left on at :62, \
             reads `-V` as that option, so the run proceeds verbosely",
        ),
    ]);
    table(&mut out, &["Flag", "Error", "Ensembl VEP"], &rows);

    out.push_str("## Accepted without effect\n\n");
    out.push_str(
        "Every flag in the table below parses and changes no annotation. The table is the error \
         table's: each flag's help is its reason, and the Ensembl VEP column cites what Ensembl \
         VEP does with it at release/116.3.\n\n",
    );
    let rows: Vec<Vec<String>> = ERROR_TABLE
        .iter()
        .filter_map(|entry| match entry.treatment {
            Treatment::AcceptedWithoutEffect { reason } => Some(vec![
                error_table_flag(entry.flag, entry.declared),
                cell(reason),
                cell(entry.ensembl_vep),
            ]),
            Treatment::Refused { .. } => None,
        })
        .collect();
    table(&mut out, &["Flag", "Reason", "Ensembl VEP"], &rows);
    let declared = compatibility_flags(command);
    if !declared.is_empty() {
        out.push_str(
            "The flags below are declared with the description their group page shows, which \
             says what, if anything, each changes; the second column is what Ensembl VEP does \
             with the flag.\n\n",
        );
        let rows: Vec<Vec<String>> = declared
            .iter()
            .map(|(flag, ensembl_vep)| vec![flag_in_code(flag), cell(ensembl_vep)])
            .collect();
        table(&mut out, &["Flag", "Ensembl VEP"], &rows);
    }

    out.push_str("## Where `Config.pm` is not followed\n\n");
    out.push_str(
        "- An output file matching `stdout` in any case does not force `--quiet` or drop \
         `--verbose` (`Config.pm:934-937`). The status lines go to standard error, so the \
         annotation on standard output is clean either way, and `--quiet --verbose` stays \
         refused with `-o STDOUT` as without. `--quiet` silences those lines and nothing else: \
         the statistics report is written in full under it, as Ensembl VEP keys its counters \
         on `--no_stats` alone (`Stats.pm:230, 277`).\n\
         - A run naming none of `--database`, `--cache`, `--offline` and `--custom` is not \
         refused (`Config.pm:993-1017`). vep-rs reads a local cache in every run, so those flags \
         select nothing and the run proceeds on the cache it is given.\n\
         - A `%name%` in an option-set value is replaced even when the flag's value is `0` or \
         empty (`Config.pm:867` leaves the reference in place then). No declared flag carries \
         such a value.\n\
         - A VCF output selected with `--output_format vcf` fills `SYMBOL`, `BIOTYPE`, `EXON` \
         and `INTRON` as `--vcf` does. Ensembl VEP's `vcf` option set applies to the `--vcf` \
         flag alone and leaves the CSQ fields empty under `--output_format vcf`.\n\
         - A value from the `--config` file, `vep.ini` or the environment never overwrites one \
         a higher source gives, even a `0`, an empty value or a bare `--cache`. Ensembl VEP's \
         `||=` (`Config.pm:1101, 1158`) lets a lower source overwrite a value Perl reads as \
         false, so `--shift_hgvs 0` on its command line yields to a file's `shift_hgvs 1`.\n\
         - A value from a file or the environment that the parser refuses (`fork many`, an \
         empty `format`) is an argument syntax error (exit 2). Ensembl VEP's `is_valid_param` \
         drops the line (`Config.pm:1046-1047`). It keeps a line whose value is itself an \
         option Getopt::Long reads (`shift_3prime --offline`), the value stored as text, which \
         Perl reads as true.\n\
         - `--buffer_size 0` annotates every record, where Ensembl VEP's buffer admits a \
         record only while it holds fewer than the value (`InputBuffer.pm:191, 212`), so \
         there a zero writes the header alone. A value with a minus sign or an underscore \
         (`-1`, `1_0`) is an argument syntax error (exit 2); Getopt::Long reads both as \
         integers, `-1` admitting no record and `1_0` read as 10.\n\
         - A byte outside UTF-8 in a value from a file or the environment becomes the \
         replacement character, where Ensembl VEP passes the byte. A `--config` value is a \
         path as written: `-`, surrounding whitespace, a `<`-led or `>`-led value and one \
         opening or closing with `|` name a file so called. Ensembl VEP opens the value with \
         Perl's two-argument `open` (`Config.pm:1075`). That reads `-` as standard input, trims \
         the whitespace, opens a `<`-led or `+<`-led value as the file after the mode, opens a \
         `>`-led value for writing, emptied, or appended to after `>>`, and runs a `|`-led or \
         `|`-closed value as a command.\n\
         - A flag read from a file or the environment acts exactly as if typed, refusals \
         included. So `show_cache_info` from a file or the environment prints the cache \
         information and exits 0. Ensembl VEP's driver reads the typed flag alone (vep:39), \
         and a lower source's value runs the annotation there instead. The status lines \
         Ensembl VEP prints under `--verbose` while reading the sources (`Read configuration \
         from <file>`, `Set '<flag>=<value>' from environment variable <VAR>`, `Config.pm:1107, \
         1145-1164`) are not printed.\n\
         - A leading `~` in `--config`, `--dir`, `--dir_cache`, `--cache <dir>`, \
         `--json_cache`, `--synonyms`, `--check_svs` and a file's `dir` is expanded to the \
         home directory. Ensembl VEP leaves the expansion to the shell. A `~/` the shell did \
         not expand (quoted, or from a file or the environment) in `--config`, `--dir`, \
         `--dir_cache`, `--cache`, `--synonyms` or `--full_cache_dir` (`--json_cache` here) is \
         passed on as typed there and names no configuration file, cache directory or \
         synonyms file. The plugins directory under `--dir` alone is globbed there \
         (`Runner.pm:908`). `--check_svs` takes no value there (`Config.pm:125`).\n\
         - `--freq_gt_lt` takes `gt` or `lt` and `--freq_filter` `exclude` or `include`; any \
         other word is an argument syntax error (exit 2). Ensembl VEP declares both as `=s` \
         (`Config.pm:137, 139`) and checks neither: an unknown direction passes no allele \
         (`BaseCacheVariation.pm:468-469`), and any filter but `include` reads as `exclude` \
         (`:320`).\n\
         - `--shift_hgvs` takes `0` or `1`; any other value is an argument syntax error \
         (exit 2). Ensembl VEP's `shift_hgvs=i` (`Config.pm:201`) takes any integer and \
         `BaseRunner.pm:498` matches it against `/(0|1)/` as a substring. So `2` leaves the \
         default shift on, and `10` sets the shifting variables to 10, which shifts HGVSc and \
         HGVSp (tested against `0`) and not HGVSg (tested against `1`).\n\
         - `--max_sv_size` takes `-1` or a size of 0 or more; a value below `-1` is an argument \
         syntax error (exit 2), a negative size being no size. Ensembl VEP's `max_sv_size=i` \
         (`Config.pm:154`) takes any integer and `Parser.pm:491` spares `-1` alone, so `-2` \
         makes every record oversize and the run drops them all, writing the header alone.\n\
         - An integer written with an underscore on `--failed`, `--clin_sig_allele`, `--fork`, \
         `--cache_version` or `--max_sv_size` (`--fork 1_0`) is an argument syntax error \
         (exit 2). Getopt::Long drops the underscores of an `=i` value \
         (`Config.pm:119, 147, 151, 154, 239`), so `1_0` reads as 10 there and \
         `--cache_version 11_6` selects release 116. `--shift_3prime` and `--shift_genomic` \
         read an integer as Getopt::Long spells one, the underscores included.\n\
         - `-V` is the version, where Getopt::Long's `ignore_case` reads it as `-v`. Beside any \
         other token `-V` is refused (exit 1) with a line naming `-v` and `--version`, unless \
         `--help` or `--version` is given too.\n\
         - `-A`, `-E`, `-I`, `-O`, `-S`, `-T` and `-U` are argument syntax errors (exit 2), \
         where `ignore_case` reads each as its lower-case short. Four upper-case letters name \
         an option, as they do in Ensembl VEP: `-J` is `--json`, `-Q` `--quiet`, `-W` \
         `--warning_file` and `-X` `--xref_refseq`. Any other upper-case letter is refused by \
         both (vep:30).\n\
         - A short with its value attached (`-ofile.txt`) and a cluster of switches (`-qv`) are \
         clap's readings, where Getopt::Long refuses the token. A `+` before an option name, \
         which Getopt::Long accepts, is not an option here.\n\
         - A flag written with `=` and nothing after it (`--tmpdir=`, `--convert=`) carries \
         the empty value, as `--tmpdir ''` does. Getopt::Long leaves such a token of a flag \
         that requires a value unread, and Ensembl VEP refuses it as an extra parameter \
         (vep:30).\n\
         - A trailing comma or a comma-joined value of `--sift`, `--polyphen`, `--terms`, \
         `--nearest`, `--compress_output` or `--format` (`--sift b,`) is an argument syntax \
         error (exit 2). Ensembl VEP's `check_config` splits such a value on commas before its \
         `%VALID` test (`Config.pm:946`), and `--sift` and `--polyphen` read its first letter \
         (`OutputFactory.pm:1815-1817`), so `--sift b,` annotates there \
         (`SIFT=tolerated(0.19)` on one \
         transcript of `21 14216193 C G`). A `0` or empty value of such a flag (`--sift 0`, \
         `--sift ''`, `--polyphen 0`, `--nearest 0`) is an argument syntax error too (exit 2). \
         For `--sift`, `--polyphen` and `--nearest` Ensembl VEP reads it as the flag unset \
         (`Config.pm:904-905, 945`) and annotates without the column (`OutputFactory.pm:1814`, \
         `AnnotationType/Transcript.pm:170`).\n\
         - The `-o` value names the output file as given. Ensembl VEP opens it with Perl's \
         two-argument `open` (`BaseRunner.pm:226`), which strips surrounding whitespace and \
         reads a leading `>` as the append mode, so `-o ' lead.txt'` writes `lead.txt` there \
         and ` lead.txt` here, and `-o '>out.txt'` appends to `out.txt` there and writes \
         `>out.txt` here.\n\
         - The `--synonyms` value names its file as given too, `-`, `0`, the empty value, a \
         value with surrounding whitespace and a `>`-led value included. Ensembl VEP's \
         two-argument `open` (`BaseVEP.pm:644`) reads `-` as standard input, trims the \
         whitespace, and opens a `>`-led value for writing, emptying the file it names and \
         reading an empty table. A value Perl reads as false leaves the flag unset there, so \
         the cache's own file is read (`CacheDir.pm:254`).\n\
         - A `--fasta` value Perl reads as false, `0` or the empty value, leaves that flag unset \
         there likewise, so the cache directory's FASTA is auto-detected or the run proceeds \
         without one (`CacheDir.pm:235`, `BaseVEP.pm:521`). vep-rs opens the value as a path \
         and refuses it with `Specified FASTA file/directory <value> not found` (exit 1).\n\
         - The `--warning_file` and `--skipped_variants_file` values name their files as \
         given likewise, `-`, surrounding whitespace and a `>`-led value included. Ensembl \
         VEP opens each for writing with the two-argument `open` (`BaseVEP.pm:894`, \
         `BaseRunner.pm:344`), which sends `-` to standard output, trims the whitespace and \
         appends a `>`-led value to the file it names. So `--warning_file -` prints the \
         warnings on standard output there and writes a file named `-` here, \
         `--warning_file ' sp.txt'` writes `sp.txt` there and ` sp.txt` here, and \
         `--skipped_variants_file '>skipped.txt'` appends to `skipped.txt` there and writes \
         `>skipped.txt` here.\n\
         - A `--skipped_variants_file` value Perl reads as false, `0` or the empty value, \
         leaves the flag unset there (`Runner.pm:176`), so no file is written. vep-rs writes a \
         file named `0`, and the empty value is refused before the first record is read, with \
         `Could not write to skipped variants file` (exit 1). A false `--warning_file` value is \
         the default warnings file in both engines (`BaseVEP.pm:890`).\n\n",
    );

    out.push_str("## The statistics report\n\n");
    out.push_str(
        "The reports `--stats_text` and `--stats_html` write, the HTML one by default, have \
         Ensembl VEP's sections, rows and charts (`Stats.pm` `dump_text` and `dump_html`), with \
         these differences:\n\n\
         - The HTML report has the same sections, ids, titles, tables and sort orders as Ensembl \
         VEP's (`Stats.pm` `dump_html`), rendered as inline SVG without Ensembl VEP's remote \
         chart loader and images. The page is one file that loads no script, no chart library \
         and no image, so it opens without a network connection. Each pie, bar and area chart \
         is drawn from the same counts, with Ensembl VEP's colours where its `COLOUR_KEYS` names \
         the key (`Constants.pm:299-346`); a click on a chart selects nothing, where the chart \
         library highlights the table row.\n\
         - Under each chart's heading the HTML report carries the chart's counts in an HTML \
         comment. The comment holds the statement Ensembl VEP's page hands its chart library \
         for that chart, `var <id> = draw<Type>(…, \
         google.visualization.arrayToDataTable([…]), …);` (`Stats.pm:955-963`), written from \
         the counts the drawing and the table show. A program that collects reports and reads \
         a chart's counts from that line of Ensembl VEP's page, as MultiQC's module for \
         Ensembl VEP does, finds the line here. Nothing in the page runs the statement, and \
         the page calls no chart library.\n\
         - The `Start time` and `End time` rows are UTC, the clock of the `## Output produced \
         at` header line, where Ensembl VEP prints local time (`Utils.pm:467-488`); `Run time` \
         is whole seconds in both.\n\
         - The `Command line options` row is the `## VEP command-line:` header's payload, so it \
         names `--json_cache` where Ensembl VEP names `--dir_cache` and carries no `--database \
         0` (`vep:35`).\n\
         - A chart's rows are ordered by count or chromosome as Ensembl VEP orders them \
         (`Stats.pm:760-790`), and two rows its order leaves tied (equal counts; `21` and \
         `chr21`) print in key order, where Ensembl VEP prints them in its hash order, which \
         changes between runs.\n\
         - The `Distribution of variants on chromosome` charts need the chromosome lengths, \
         which vep-rs reads from the index of the run's reference FASTA alone, one chart per \
         chromosome spelled as the FASTA names it. Ensembl VEP reads them from its database \
         when one is connected and from the FASTA otherwise (`BaseVEP.pm:550-580`). Without a \
         FASTA vep-rs writes no distribution chart, and Ensembl VEP writes one only from a \
         connected database. A vep-rs run that writes a report without a FASTA says so in \
         one `INFO:` line on standard error, which `--quiet` withholds.\n\
         - The report is checked for an existing file before the output and the \
         skipped-variants file are, as Ensembl VEP checks it (`Runner.pm:144-148`, before \
         `:202`).\n\
         - Before its first record the run tests that it can write the report, the output, \
         a warnings file that `--warning_file` names, then the skipped-variants file, a pipe \
         or a device excepted. The first that fails is refused by name, and no report, \
         output or skipped-variants file is created or emptied before every check and test \
         has passed; a warning drawn earlier is written to the warnings file. So such a \
         refusal leaves those earlier files as they were and no empty report, where Ensembl \
         VEP leaves one behind a refused output.\n\n",
    );

    out.push_str("## Configuration sources\n\n");
    out.push_str(
        "A flag of `vep`, and an option of `vep recode`, is read from four sources, the first \
         that names it winning, as Ensembl VEP reads them for `vep` and `variant_recoder` alike \
         (`Config.pm:693-708`):\n\n",
    );
    writeln!(
        out,
        "1. the command line;\n\
         2. the `--config` file;\n\
         3. `<dir>/{INI_FILE}`, where `<dir>` is `--dir` (the command line's, else the \
         `--config` file's) and `{DEFAULT_DIR}` by default;\n\
         4. the environment variable `{ENV_PREFIX}<FLAG>` (`{ENV_PREFIX}DIR_PLUGINS` sets \
         `--dir_plugins`).\n"
    )
    .unwrap();
    let (collected, refused) = declared_and_refused(command, ALLOW_MULTIPLE);
    writeln!(
        out,
        "A repeatable flag ({}) collects every source's values, the command line's first.{} No \
         source behind the command line can set {}.\n",
        flags_in_code(&collected),
        refused_in_list(&refused, "repeatable"),
        flags_in_code(COMMAND_LINE_ONLY)
    )
    .unwrap();
    writeln!(
        out,
        "### The `--config` file and `{INI_FILE}` (`Config.pm:1068-1110`)\n"
    )
    .unwrap();
    writeln!(
        out,
        "- One flag per line: the flag's long name without its dashes (one leading dash is \
         allowed), a blank, then its value. A line starting with `#`, or blank, is skipped.\n\
         - A switch is on with the value `1` and off with `0`, except `var_synonyms`, whose \
         `0` switches `check_existing` on, with the notice `INFO: Colocated variants \
         will be provided when using --var_synonyms`, and writes no synonyms, as Ensembl VEP \
         tests that flag with `defined` (Runner.pm:117). `cache` is on with `1`, a \
         directory or an empty quoted value and off with `0`. A key alone on its line is \
         dropped. Ensembl VEP drops a scalar or switch key alone on its line too \
         (`Config.pm:1101`) and keeps an empty value for a repeatable key alone on its line \
         (`Config.pm:1094-1097`): `plugin` warns `Failed to compile plugin` there, and `custom` \
         and `transcript_filter` are fatal.\n\
         - A flag taking a value takes the first word after its name, or, for a repeatable \
         flag, the rest of the line. Blanks inside double quotes are part of a value and the \
         quotes are removed. A literal `___SPACE___` in a value stays as written. Ensembl \
         VEP's reader uses that text as its placeholder for a blank inside double quotes and \
         turns a literal one into a blank too (`Config.pm:1082, 1089`), so a file named with \
         it is unreachable there.\n\
         - `shift_3prime` and `shift_genomic` are on with any integer but `0`, read as text: \
         `00`, `+0` and `-0`, which are 0 and off on the command line, are on from a file or \
         the environment. Ensembl VEP keeps a lower source's value as the text it reads \
         (`Config.pm:1101, 1158`) and tests it for truth, where Getopt::Long stores the \
         command line's value as a number.\n\
         - A key is the long name exactly. A prefix, another case or a name that is no option \
         of Ensembl VEP's or of this version names nothing, and the line is dropped. So is an \
         alias, unless the alias is an option of Ensembl VEP's own (`full_cache_dir` names \
         `--json_cache`). Ensembl VEP drops an unknown name too. A prefix or another case of \
         one of its options passes `is_valid_param` (Getopt::Long's `ignore_case` and \
         `auto_abbrev`) and is stored under the written key, which nothing reads: the \
         annotation is unchanged, but the `## VEP command-line:` header echoes the key \
         (`--PROTEIN_VERSION` for `PROTEIN_VERSION 1`, `Config.pm:734-745`).\n\
         - A key of a flag this version refuses (`database 1`, `af_exac 1`) reaches that \
         refusal (exit 1). The exception is a switch set to `0` or the empty value: off \
         there, it adds nothing here, and `af_exac 0` adds nothing either. `convert 0` is \
         read as `--convert 0` is, and the run proceeds.\n\
         - Of several lines naming one flag, the first whose value is not `0` or empty counts.\n\
         - A `--config` file that cannot be opened is `ERROR: [config] Could not open config \
         file \"<file>\"`. A `--config` value of `0` or the empty value is that error with the \
         empty name, before any file is looked for (`Config.pm:1073`). A `--config` path that \
         is a directory reads as nothing, as Perl's `open` of a directory yields no line \
         (`:1075-1077`). `{INI_FILE}` is read only when it is a file.\n\n\
         ```\n# flags of a run\ncache 1\noffline 1\ndir_cache /data/vep\nspecies homo_sapiens\n\
         assembly GRCh38\nplugin LoFtool,/data/LoFtool_scores.txt\n```\n"
    )
    .unwrap();
    out.push_str("### Environment variables (`Config.pm:1127-1168`)\n\n");
    writeln!(
        out,
        "`{ENV_PREFIX}<FLAG>=<value>` sets `--<flag>`, the name lower-cased, when the value is \
         not empty. The name is read as a configuration file's key is: a name that is no option \
         of Ensembl VEP's or of this version is ignored. One of a flag this version refuses \
         reaches that refusal (`{ENV_PREFIX}DATABASE=1` is the `--database` error, exit 1). \
         The variables are read in name order.\n"
    )
    .unwrap();
    writeln!(
        out,
        "No variable sets `--phyloP` or `--phastCons`, the two flags whose names hold a capital \
         letter: `{ENV_PREFIX}PHYLOP=100` is read as `phylop`, which names nothing, and the run \
         adds no source and prints no message. Ensembl VEP reads the variable the same way: it \
         stores the value under `phylop`, which its option set never reads, and echoes \
         `--phylop 100` in the `## VEP command-line:` header (`Config.pm:1140, 1158`, `:553`). \
         The command line and a configuration file's `phyloP` key set the flag, and \
         `{ENV_PREFIX}CUSTOM=<bigWig>,phyloP100way,bigwig,exact` adds from the environment the \
         source `--phyloP 100` names.\n"
    )
    .unwrap();
    out.push_str(
        "Environment variables vep-rs reads outside this mechanism, which Ensembl VEP has no \
         counterpart for:\n\n",
    );
    let rows: Vec<Vec<String>> = ENVIRONMENT
        .iter()
        .map(|(name, help)| vec![code(name), cell(help)])
        .collect();
    table(&mut out, &["Variable", "Meaning"], &rows);

    out.push_str("## Option spellings (`Config.pm:62`)\n\n");
    out.push_str(
        "- Getopt::Long's defaults apply to the command line: a long option is read in any case \
         (`--FORMAT vcf`), abbreviated to any prefix that names one option (`--form vcf`), and \
         after one dash as well as two (`-cache`, `-form vcf`).\n\
         - A value follows its option or is attached with `=`, and the token after an option \
         that requires a value is that value whatever it looks like (`--input_file -cache` \
         names the file `-cache`).\n\
         - A prefix is resolved against this version's flags first. Where a flag vep-rs adds \
         makes it ambiguous, it is resolved as Getopt::Long resolves it against Ensembl VEP's \
         `VEP_PARAMS` list (`Config.pm:73-276`): `--js` is `--json` though `--json_cache` also \
         starts with it, `--pa` `--password`, `--ve` `--verbose`.\n\
         - A prefix that names one flag of this version resolves to it (`--bio` is `--biotype`). \
         A prefix of several flags that Getopt::Long resolves no better (`--fo`), or a name no \
         option has, is an argument syntax error (exit 2).\n",
    );
    out.push_str(VEP_ABBREVIATIONS);
    out.push_str(
        "- A single-letter short is clap's and is matched by case: `-v` is `--verbose`, `-h` \
         the help and `-V` the version, `-qv` the cluster `-q -v`, and `-ofile.txt` the value \
         `file.txt` of `-o`, a spelling Getopt::Long refuses. Ensembl VEP refuses `-h`, a \
         prefix of eight of its options, and has no `--version`.\n\n",
    );
    out.push_str(
        "Of the spellings Getopt::Long resolves to an option this version declares, these \
         resolve differently here:\n\n",
    );
    let unresolved = unresolved_perl_spellings(command);
    if unresolved.is_empty() {
        out.push_str("- none: every such spelling names the same option.\n\n");
    } else {
        let rows: Vec<Vec<String>> = unresolved
            .iter()
            .map(|u| {
                vec![
                    flag_in_code(&u.spelling),
                    flag_in_code(&u.perl),
                    if u.candidates.is_empty() {
                        "no flag".to_string()
                    } else {
                        flags_in_code(&u.candidates.iter().map(String::as_str).collect::<Vec<_>>())
                    },
                ]
            })
            .collect();
        table(
            &mut out,
            &[
                "Spelling",
                "Ensembl VEP reads",
                "vep-rs flags it could name",
            ],
            &rows,
        );
    }

    out.push_str("## Exit status and errors\n\n");
    out.push_str(
        "- `0`: the run completed; also `--help`, `--version` and an empty command line.\n\
         - `1`: a refusal or a failed run. `ERROR: [<component>] <message>` is written to \
         standard error, on one line unless the message carries a second. `<message>` carries \
         Ensembl VEP's text, without its own `ERROR: ` prefix, where vep-rs quotes that text, \
         and is in vep-rs's own words elsewhere: for example, the `--af_exac` refusal names the \
         flag to use, where Ensembl VEP reports an unexpected parameter. The components are \
         listed below.\n\
         - `2`: an argument syntax error, in clap's words with its `error:` prefix: an unknown \
         flag, a prefix of several flags, a missing value, a value outside an enumeration.\n\
         - Two refusals carry Ensembl VEP's text in parentheses after vep-rs's own words, as \
         `(Ensembl VEP: <text>)`: the `--lrg` refusal and the refusal of a GFF or GTF \
         transcript source without a FASTA. Alone, either text would point to a database \
         connection, which vep-rs never makes.\n\
         - Three refusals carry a second line, as Ensembl VEP's own messages do. A refused \
         `--custom` string carries `LINE: --custom <string>` where Ensembl VEP prints it. An \
         unsupported `--regulatory_gff` key carries `Supported keys: ` and the keys the flag \
         reads, and a `--cell_type` name the cache's regulatory data lacks carries the cache's \
         cell types.\n\
         - A compressed input cut inside a gzip member is a `run` failure, exit 1 with nothing \
         past the header written, where Ensembl VEP writes the records before the cut and \
         exits 0 ([the divergence page](../divergences/compressed-input.md)).\n\
         - A directory named as the input is a `run` failure, exit 1 with no output written, \
         where Ensembl VEP passes the existence test (`Parser.pm:239`), leaves the directory \
         uncompressed (`:242` tests `-f`), opens it, reads nothing from it and exits 0 with the \
         header alone written.\n\
         - A line of the run that standard error refuses (the banner, a status, progress or \
         summary line, a warning's echo, the `ERROR:` line) is dropped and never changes the \
         exit status. So is a line a plugin prints itself, and every line of \
         [`vep filter`](filter.md), [`vep recode`](recode.md), `vep-cache-builder` and \
         `vep-cache-converter`. A completed run exits 0 on a full device or a pipe without a \
         reader, and a refusal exits 1. Ensembl VEP completes on a full device too; on a pipe \
         without a reader it completes when it has nothing to write there and is killed by \
         `SIGPIPE` (exit 141, no output) at its first warning (`BaseVEP.pm:821`).\n\
         - A `--warning_file` that cannot be written is an `output` failure, exit 1, raised \
         before the first record is read, with nothing written, whether or not the run would \
         warn; under `--show_cache_info`, once a warning is drawn and before the print. A \
         default warnings file, a pipe or a device that cannot be written fails only at the \
         warning that writes it: the run ends once any batch in hand is written, and a \
         Parquet run then publishes nothing. Ensembl VEP opens every warnings file at its \
         first warning and dies where that fails (`BaseVEP.pm:894`), with any buffer in \
         hand unwritten and no warning echoed.\n\
         - The process ends as soon as its outputs are written and closed: the loaded cache is \
         left to the operating system, which reclaims it in one step, instead of being freed \
         one allocation at a time. An allocator asked to report at exit (`MIMALLOC_SHOW_STATS`, \
         `MIMALLOC_VERBOSE`) gets the ordinary exit, so its report prints.\n\n",
    );
    let rows: Vec<Vec<String>> = COMPONENTS
        .iter()
        .map(|(component, meaning)| vec![code(component), cell(meaning)])
        .collect();
    table(&mut out, &["Component", "Meaning"], &rows);
    out.push_str(
        "`--help` and `--version` win over every other token of the command line, wherever they \
         stand before a `--` (vep:29 before vep:30 under Getopt::Long's `pass_through`, \
         Config.pm:62); `--help=1` and a `--help` behind `--` are refused. A `--help` or \
         `--version` that is the token after a flag taking a string is that value (`-i \
         --help` names the file) on both engines. After a flag taking a number it is the help \
         or the version, as it is after a flag of this version's own, which Ensembl VEP \
         leaves in `@ARGV` (`--pipeline_depth --help`).\n\n",
    );
    out.push_str(
        "A panic aborts the process with the panic's own status. Ensembl VEP's own exit codes \
         (255 on `die`, 2 on a failed open) are not copied.\n",
    );
    out
}

fn error_table_flag(flag: &str, declared: bool) -> String {
    if declared {
        flag_in_code(flag)
    } else {
        format!("`--{flag}` (not declared; caught before parsing)")
    }
}

/// `COMPATIBILITY_FLAGS` checked against the help texts: every declared flag
/// whose help reads `Accepted for VEP compatibility` is listed, and every
/// listed flag still reads so.
fn compatibility_flags(command: &Command) -> &'static [(&'static str, &'static str)] {
    let from_help: Vec<&str> = command
        .get_arguments()
        .filter(|a| {
            flag_help(a)
                .to_ascii_lowercase()
                .contains(COMPATIBILITY_HELP)
        })
        .filter_map(|a| a.get_long())
        .collect();
    for flag in &from_help {
        assert!(
            COMPATIBILITY_FLAGS.iter().any(|(f, _)| f == flag),
            "--{flag} is accepted for compatibility: add it to COMPATIBILITY_FLAGS with what \
             Ensembl VEP does with it"
        );
    }
    for (flag, _) in COMPATIBILITY_FLAGS {
        assert!(
            from_help.contains(flag),
            "--{flag}'s help does not read as accepted for compatibility: remove it from \
             COMPATIBILITY_FLAGS"
        );
    }
    COMPATIBILITY_FLAGS
}

/// Every page as `(file under docs/, text)`: the landing page, the group
/// pages and the companion tools' pages. A `--help` heading no group carries
/// is a stale `GROUPS`.
fn render(command: &Command) -> Vec<(&'static str, String)> {
    let sections = flag_sections(command);
    let mut pages = vec![(LANDING, render_landing(command))];
    for group in GROUPS {
        let text = if group.file == COMPATIBILITY {
            render_compatibility(group, command, &sections)
        } else {
            render_group(group, &sections)
        };
        pages.push((group.file, text));
    }
    pages.push((FILTER_PAGE, render_filter()));
    pages.push((RECODE_PAGE, render_recode()));
    for (heading, _) in &sections {
        assert!(
            GROUPS
                .iter()
                .any(|g| g.sections.iter().any(|(h, _)| h == heading)),
            "the --help heading {heading:?} belongs to no group page: add it to GROUPS"
        );
    }
    pages
}

#[test]
fn cli_reference_is_current() {
    let mut command = Args::command_with_error_table();
    command.build();
    let update = std::env::var_os("VEP_UPDATE_CLI_REFERENCE").is_some();
    for (file, rendered) in render(&command) {
        let path = repo_path(file);
        if update {
            fs::write(&path, &rendered).unwrap();
            continue;
        }
        let committed = fs::read_to_string(&path).unwrap_or_default();
        if committed != rendered {
            let first = committed
                .lines()
                .zip(rendered.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(committed.lines().count().min(rendered.lines().count()));
            panic!(
                "{file} is not the rendering of the command line and the tables: first \
                 difference at line {}\n  committed: {:?}\n  rendered:  {:?}\nrun {REGENERATE}",
                first + 1,
                committed.lines().nth(first).unwrap_or(""),
                rendered.lines().nth(first).unwrap_or("")
            );
        }
    }
}

/// The `--help` heading `command` lists `flag` under.
fn heading_of<'a>(command: &'a Command, flag: &str) -> &'a str {
    command
        .get_arguments()
        .find(|arg| arg.get_long() == Some(flag))
        .and_then(Arg::get_help_heading)
        .unwrap_or_else(|| panic!("--{flag} has no help heading"))
}

/// The headings `--help` prints over its flag lists, in its order: clap
/// prints a heading where it first meets a flag that carries it. `Options`,
/// which holds clap's own `--help` and `--version`, is left out.
fn help_headings(command: &mut Command) -> Vec<String> {
    command
        .render_help()
        .to_string()
        .lines()
        .filter(|line| !line.starts_with(' '))
        .filter_map(|line| line.strip_suffix(':'))
        .filter(|heading| *heading != "Options")
        .map(str::to_string)
        .collect()
}

/// A heading's anchor on the site, as mdBook derives it: letters and digits
/// in lower case, `_` and `-` kept, each blank a `-`, every other character
/// dropped.
fn anchor(heading: &str) -> String {
    heading
        .chars()
        .filter_map(|c| {
            if c.is_alphanumeric() || c == '_' || c == '-' {
                Some(c.to_ascii_lowercase())
            } else if c.is_whitespace() {
                Some('-')
            } else {
                None
            }
        })
        .collect()
}

/// Each `--help` heading of `vep` names what its flags do. The shift flags,
/// the co-located variant flags, the two flags that drop rows, the allele and
/// position columns, the reference allele flags, the custom annotation flags
/// and the two flags that keep a record a run would skip each stand under a
/// heading of their own; `--mirna`, `--lrg` and `--tmpdir` stand under the
/// heading of their kind. `--help` prints the headings in this order, and
/// none reads `Miscellaneous`.
#[test]
fn flags_of_one_kind_share_a_help_heading() {
    let mut command = Args::command_with_error_table();
    command.build();
    let headings: [(&str, &[&str]); 11] = [
        (
            "3' Shifting",
            &[
                "shift_3prime",
                "shift_hgvs",
                "shift_genomic",
                "shift_length",
            ],
        ),
        (
            "Co-located Variants",
            &[
                "check_existing",
                "no_check_alleles",
                "exclude_null_alleles",
                "failed",
                "clin_sig_allele",
                "clinvar_somatic_classification",
                "var_synonyms",
                "check_svs",
            ],
        ),
        ("Dropped Rows", &["coding_only", "no_intergenic"]),
        (
            "Alleles / Positions",
            &[
                "allele_number",
                "show_ref_allele",
                "uploaded_allele",
                "ambiguity",
                "minimal",
                "total_length",
            ],
        ),
        (
            "Reference Allele",
            &["use_transcript_ref", "use_given_ref", "bam"],
        ),
        ("Plugins", &["plugin", "dir_plugins", "safe", "debug"]),
        (
            "Custom Annotation",
            &[
                "custom",
                "bigwig",
                "phyloP",
                "phastCons",
                "ucsc_assembly",
                "ucsc_data_root",
                "custom_multi_allelic",
                "custom_suppress_filter",
                "gff",
                "gtf",
            ],
        ),
        ("Record Handling", &["dont_skip", "allow_non_variant"]),
        (
            "Annotation Fields",
            &["symbol", "hgvs", "regulatory", "mirna"],
        ),
        ("Database", &["database", "lrg"]),
        ("Runtime", &["buffer_size", "tmpdir"]),
    ];
    for (heading, flags) in headings {
        for flag in flags {
            assert_eq!(heading_of(&command, flag), heading, "--{flag}");
        }
    }
    assert_eq!(
        help_headings(&mut command),
        [
            "Input",
            "Output",
            "Cache / Offline",
            "3' Shifting",
            "Species / Assembly",
            "Database",
            "Runtime",
            "Annotation Fields",
            "Co-located Variants",
            "Dropped Rows",
            "Allele Frequencies",
            "Filtering / Picking",
            "Alleles / Positions",
            "Plugins",
            "Custom Annotation",
            "Transcript Sets",
            "Reference Allele",
            "Statistics",
            "Verbosity",
            "Record Handling",
        ]
    );
}

/// `vep recode` lists the six flags that tune the annotation behind the
/// translation and select no transcript under a heading of their own, printed
/// after `Transcript selection`.
#[test]
fn recode_flags_that_select_no_transcript_have_their_own_heading() {
    let mut command = RecodeArgs::command();
    command.build();
    for flag in [
        "minimal",
        "buffer_size",
        "failed",
        "shift_hgvs",
        "use_transcript_ref",
        "use_given_ref",
    ] {
        assert_eq!(heading_of(&command, flag), "Annotation run", "--{flag}");
    }
    for flag in ["merged", "transcript_filter", "pick", "pick_order"] {
        assert_eq!(
            heading_of(&command, flag),
            "Transcript selection",
            "--{flag}"
        );
    }
    assert_eq!(
        help_headings(&mut command),
        [
            "Input",
            "Output",
            "Cache, species and assembly",
            "Transcript selection",
            "Annotation run",
            "Database",
        ]
    );
}

/// Each page carries one section per kind of flag, in this order, and a link
/// from a page's text to a section of the same page names a section the page
/// has, so a lead that names its sections cannot outlive one of them.
#[test]
fn pages_carry_a_section_per_kind_of_flag() {
    let mut command = Args::command_with_error_table();
    command.build();
    let pages: BTreeMap<&str, String> = render(&command).into_iter().collect();
    let sections_of = |file: &str| -> Vec<&str> {
        pages[file]
            .lines()
            .filter_map(|line| line.strip_prefix("## "))
            .collect()
    };
    assert_eq!(
        sections_of("docs/cli/cache-species-runtime.md"),
        [
            "Cache and offline",
            "3' shifting",
            "Species and assembly",
            "Runtime"
        ]
    );
    assert_eq!(
        sections_of("docs/cli/output-control.md"),
        ["Annotation fields", "Co-located variants", "Dropped rows"]
    );
    assert_eq!(
        sections_of("docs/cli/filtering-picking.md"),
        [
            "Filtering and picking",
            "Alleles and positions",
            "Transcript sets",
            "Reference allele"
        ]
    );
    assert_eq!(
        sections_of("docs/cli/plugins.md"),
        ["Loading plugins", "Custom annotation"]
    );
    assert_eq!(
        sections_of("docs/cli/statistics-verbosity-misc.md"),
        ["Statistics", "Verbosity", "Record handling"]
    );
    assert_eq!(
        sections_of(RECODE_PAGE),
        [
            "Input",
            "Output",
            "Cache, species and assembly",
            "Transcript selection",
            "Annotation run",
            "Database",
        ]
    );
    assert!(
        pages["docs/cli/allele-frequencies.md"]
            .contains("[Output control](output-control.md#co-located-variants)"),
        "{}",
        pages["docs/cli/allele-frequencies.md"]
    );
    for (file, text) in &pages {
        let anchors: Vec<String> = text
            .lines()
            .filter(|line| line.starts_with('#'))
            .map(|line| anchor(line.trim_start_matches('#').trim()))
            .collect();
        for link in text.split("](#").skip(1) {
            let target = link.split(')').next().unwrap_or_default();
            assert!(
                anchors.iter().any(|a| a == target),
                "{file} links #{target}, a section it does not have"
            );
        }
    }
}

/// The option-set table prints Ensembl VEP's text, so its `--phyloP` row
/// carries the source's short name as `Config.pm:555` spells it. The sentence
/// under the table says whose text the row is and links the page of that
/// difference.
#[test]
fn option_set_table_says_whose_spelling_its_phylop_row_is() {
    let mut command = Args::command_with_error_table();
    command.build();
    let landing = render_landing(&command);
    let row = landing
        .find("| `--phyloP` | `--custom ")
        .expect("the option-set table has a --phyloP row");
    let row_text = landing[row..].lines().next().unwrap_or_default();
    assert!(
        row_text.contains(",phlyoP%phyloP%way,bigwig,exact`"),
        "{row_text}"
    );
    let note = landing
        .find(PHYLOP_ROW_NOTE)
        .expect("the landing page carries the note on the --phyloP row");
    let next_table = landing
        .find("### Incompatible flags")
        .expect("the incompatible-flags table follows the option sets");
    assert!(row < note && note < next_table);
    assert!(PHYLOP_ROW_NOTE.contains("(divergences/phylop-short-name.md)"));
    assert!(repo_path("docs/divergences/phylop-short-name.md").is_file());
}

/// Every table-backed flag's help is the table's text, verbatim, so `--help`
/// and the page cannot drift from the refusal.
#[test]
fn help_of_table_backed_flags_is_the_tables() {
    let mut command = Args::command_with_error_table();
    command.build();
    let help: BTreeMap<&str, String> = command
        .get_arguments()
        .filter_map(|a| a.get_long().map(|l| (l, flag_help(a))))
        .collect();
    for entry in ERROR_TABLE.iter().filter(|e| e.declared) {
        assert_eq!(
            help.get(entry.flag),
            Some(&entry.help()),
            "--{}",
            entry.flag
        );
    }
}

/// The compatibility page's sentences on an upper-case letter are the parser's
/// readings. A letter reaches an option where the pre-pass spells it as one:
/// the four that are a prefix of one option of `@VEP_PARAMS`. Every other but
/// `-V` is an unknown argument, seven of them the upper-case form of a short
/// that `@VEP_PARAMS` declares, which Getopt::Long's `ignore_case`
/// (`Config.pm:62`) reads as that short. `-V` keeps the version's help text,
/// which the refusals table's lead-in says of its row.
#[test]
fn upper_case_letters_are_read_as_the_compatibility_page_says() {
    let mut command = Args::command_with_error_table();
    command.build();
    let dir = tempfile::tempdir().unwrap();
    let perl = perl_options();
    let perl_shorts: Vec<String> = perl
        .iter()
        .flat_map(|o| o.names.iter())
        .filter(|n| n.len() == 1)
        .cloned()
        .collect();
    let mut named = Vec::new();
    let mut refused_shorts = Vec::new();
    let mut refused_by_both = Vec::new();
    for letter in ('A'..='Z').filter(|l| *l != 'V') {
        let home = (OsString::from("HOME"), dir.path().as_os_str().to_owned());
        let argv = prepass(vec!["vep".into(), format!("-{letter}").into()], [home]).unwrap();
        let spelled = argv[1].to_str().unwrap();
        let long = spelled
            .strip_prefix("--")
            .filter(|long| command.get_arguments().any(|a| a.get_long() == Some(*long)));
        let in_perl =
            vep_cli::argv::resolve(&letter.to_string(), &perl).map(|o| o.canonical.as_str());
        if let Some(long) = long {
            assert_eq!(in_perl, Some(long), "-{letter}");
            named.push((letter, long.to_string()));
            continue;
        }
        let error = command.clone().try_get_matches_from(&argv).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::UnknownArgument,
            "{argv:?}"
        );
        if perl_shorts.contains(&letter.to_ascii_lowercase().to_string()) {
            refused_shorts.push(letter);
        } else {
            assert_eq!(in_perl, None, "-{letter}");
            refused_by_both.push(letter);
        }
    }
    let named: Vec<(char, &str)> = named.iter().map(|(l, long)| (*l, long.as_str())).collect();
    assert_eq!(
        named,
        [
            ('J', "json"),
            ('Q', "quiet"),
            ('W', "warning_file"),
            ('X', "xref_refseq")
        ]
    );
    assert_eq!(refused_shorts, ['A', 'E', 'I', 'O', 'S', 'T', 'U']);
    assert_eq!(refused_by_both.len(), 14, "{refused_by_both:?}");
    let version = command
        .get_arguments()
        .find(|a| a.get_short() == Some('V'))
        .expect("the version's short");
    assert_eq!(version.get_long(), Some("version"));
    assert_eq!(flag_help(version), "Print version");
    let pages: BTreeMap<&str, String> = render(&command).into_iter().collect();
    for sentence in [
        "- `-A`, `-E`, `-I`, `-O`, `-S`, `-T` and `-U` are argument syntax errors (exit 2), where \
         `ignore_case` reads each as its lower-case short. Four upper-case letters name an \
         option, as they do in Ensembl VEP: `-J` is `--json`, `-Q` `--quiet`, `-W` \
         `--warning_file` and `-X` `--xref_refseq`. Any other upper-case letter is refused by \
         both (vep:30).\n",
        "and the help text of a declared one is this table's message. The `-V` row is the \
         exception: `-V` is the version's short, and its help text is the version's.\n\n",
    ] {
        assert!(
            pages[COMPATIBILITY].contains(sentence),
            "{}",
            pages[COMPATIBILITY]
        );
    }
    assert!(!pages[COMPATIBILITY].contains("an upper-case short names no flag here"));
}

/// The repeatable and comma-separated flags the pages list as collected or
/// split are the ones the command line declares: every repeatable flag is
/// declared, so the landing page's list and the compatibility page's
/// configuration-sources sentence name none as refused, and a configuration
/// file's `phyloP` or `cell_type` line reaches the parser as the flag. The one
/// comma-separated flag this version lacks, `--dont_export`, is refused as an
/// unknown argument and named as refused instead; `@VEP_PARAMS` does not
/// declare it, so the file reader drops its line. `--cell_type` is listed with
/// the note that its value is read by longest match against the cache's names.
#[test]
fn repeatable_flags_list_the_declared_ones_and_name_the_refused() {
    let mut command = Args::command_with_error_table();
    command.build();
    let (collected, refused) = declared_and_refused(&command, ALLOW_MULTIPLE);
    assert_eq!(
        collected,
        [
            "plugin",
            "custom",
            "phyloP",
            "phastCons",
            "transcript_filter"
        ]
    );
    assert!(refused.is_empty(), "{refused:?}");
    let (split, refused_lists) = declared_and_refused(&command, LIST_FLAGS);
    assert_eq!(
        split,
        [
            "individual",
            "individual_zyg",
            "cell_type",
            "pick_order",
            "fields",
            "chr",
            "distance"
        ]
    );
    assert_eq!(refused_lists, ["dont_export"]);
    let unknown_argument = |argv: &[&str]| {
        let error = command.clone().try_get_matches_from(argv).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::UnknownArgument,
            "{argv:?}"
        );
    };
    for flag in &refused_lists {
        unknown_argument(&["vep", &format!("--{flag}"), "100"]);
    }
    let dir = tempfile::tempdir().unwrap();
    let config = dir.path().join("config.txt");
    fs::write(
        &config,
        "cell_type a,b\ndont_export a\nphyloP 100\nphastCons 100\n",
    )
    .unwrap();
    let home = (OsString::from("HOME"), dir.path().as_os_str().to_owned());
    let from_file = prepass(
        vec!["vep".into(), "--config".into(), config.into_os_string()],
        [home],
    )
    .unwrap();
    let from_file: Vec<&str> = from_file.iter().map(|t| t.to_str().unwrap()).collect();
    for token in ["--cell_type=a,b", "--phyloP=100", "--phastCons=100"] {
        assert!(from_file.contains(&token), "{from_file:?}");
        assert!(
            command.clone().try_get_matches_from(["vep", token]).is_ok(),
            "{token}"
        );
    }
    assert!(
        !from_file.iter().any(|t| t.contains("dont_export")),
        "{from_file:?}"
    );
    assert_eq!(refused_in_list(&[], "repeatable"), "");
    assert_eq!(
        refused_in_list(&["phyloP"], "repeatable"),
        " `--phyloP`, repeatable in Ensembl VEP, is an option this version lacks, refused on the \
         command line and in a configuration file (exit 2)."
    );
    let undeclared_sentence = " `--dont_export`, comma-separated in Ensembl VEP's list, is an \
                               option neither its `vep` script (`Config.pm:73-276`) nor this \
                               version declares: unknown on the command line (exit 2) and dropped \
                               from a configuration file, as in Ensembl VEP.";
    assert_eq!(
        refused_in_list(&["dont_export"], "comma-separated"),
        undeclared_sentence
    );
    let pages: BTreeMap<&str, String> = render(&command).into_iter().collect();
    let lists_line = format!(
        "- Comma-separated values: `--individual`, `--individual_zyg`, `--cell_type`, \
         `--pick_order`, `--fields`, `--chr`, `--distance`.{CELL_TYPE_LIST_NOTE}\
         {undeclared_sentence}\n"
    );
    assert!(pages[LANDING].contains(&lists_line), "{}", pages[LANDING]);
    let landing_line =
        "- Repeatable, each occurrence appending: `--plugin`, `--custom`, `--phyloP`, \
                        `--phastCons`, `--transcript_filter`.\n";
    assert!(pages[LANDING].contains(landing_line), "{}", pages[LANDING]);
    let sources_sentence = "A repeatable flag (`--plugin`, `--custom`, `--phyloP`, `--phastCons`, \
                            `--transcript_filter`) collects every source's values, the command \
                            line's first. No source behind the command line can set ";
    assert!(
        pages[COMPATIBILITY].contains(sources_sentence),
        "{}",
        pages[COMPATIBILITY]
    );
    for page in pages.values() {
        assert!(!page.contains("option this version lacks"), "{page}");
    }
}

/// The long names of the flags `command` hides from `--help`.
fn hidden_flags(command: &Command) -> Vec<&str> {
    command
        .get_arguments()
        .filter(|arg| arg.is_hide_set())
        .filter_map(Arg::get_long)
        .collect()
}

/// A page's table rows whose first cell, the flag's names, ends with the mark.
fn marked_lines(page: &str) -> Vec<&str> {
    page.lines()
        .filter(|line| {
            line.strip_prefix("| `--")
                .and_then(|row| row.split(" | ").next())
                .is_some_and(|names| names.ends_with(HIDDEN_MARK))
        })
        .collect()
}

/// Whether a help text lists `flag`: one of its option lines, which clap
/// indents by at most six columns, opens with the long name after any short
/// form.
fn help_lists(help: &str, flag: &str) -> bool {
    let name = format!("--{flag}");
    help.lines().any(|line| {
        let option = line.trim_start();
        let after_short = match option.split_once(", ") {
            Some((short, rest)) if short.len() == 2 && short.starts_with('-') => rest,
            _ => option,
        };
        line.len() - option.len() <= 6
            && after_short
                .strip_prefix(&name)
                .is_some_and(|rest| !rest.starts_with(|c: char| c.is_alphanumeric() || c == '_'))
    })
}

/// `--help` prints every flag but the ones the parser hides, and the
/// reference lists each hidden flag on one row that carries the mark; no
/// other row carries it. A page with a marked row says once, under the table
/// that holds its first one, what the mark means, with its own count of such
/// rows, and a page without one does not. The companion tools hide no flag,
/// so their pages carry no such sentence.
#[test]
fn hidden_flags_carry_the_mark_and_their_page_says_once_what_it_means() {
    let mut command = Args::command_with_error_table();
    command.build();
    let hidden = hidden_flags(&command);
    for help in [
        command.clone().render_help().to_string(),
        command.clone().render_long_help().to_string(),
    ] {
        for flag in command.get_arguments().filter_map(Arg::get_long) {
            assert_eq!(
                help_lists(&help, flag),
                !hidden.contains(&flag),
                "--{flag} in --help"
            );
        }
    }
    for mut tool in [FilterArgs::command(), RecodeArgs::command()] {
        tool.build();
        assert!(hidden_flags(&tool).is_empty(), "{}", tool.get_name());
    }

    let pages = render(&command);
    let mut marked = Vec::new();
    for (file, text) in &pages {
        let rows = marked_lines(text);
        let notes = text.matches(HIDDEN_NOTE).count();
        if rows.is_empty() {
            assert_eq!(notes, 0, "{file} marks no flag");
        } else {
            assert_eq!(notes, 1, "{file} marks {} flags", rows.len());
            let note = hidden_note(rows.len());
            let table_end = text.find(rows[0]).unwrap();
            let after_table = text[table_end..].find("|\n\n").unwrap() + table_end + 3;
            assert!(text[after_table..].starts_with(&note), "{file}: {note}");
        }
        marked.extend(rows);
    }
    assert_eq!(marked.len(), hidden.len(), "{marked:?}");
    for flag in &hidden {
        let row = format!("| `--{flag}`");
        assert_eq!(
            marked.iter().filter(|line| line.starts_with(&row)).count(),
            1,
            "--{flag}"
        );
    }
}

/// The sentence on hidden flags names three ways a run gives one: the command
/// line, a configuration file and a `VEP_<FLAG>` variable. The pre-pass spells
/// each hidden flag from all three, as it spells any declared flag.
#[test]
fn a_hidden_flag_is_read_from_each_source_its_sentence_names() {
    let mut command = Args::command_with_error_table();
    command.build();
    let dir = tempfile::tempdir().unwrap();
    let home = || (OsString::from("HOME"), dir.path().as_os_str().to_owned());
    let names = |argv: Vec<OsString>, flag: &str| {
        let option = format!("--{flag}");
        let with_value = format!("--{flag}=1");
        argv.iter().any(|t| *t == *option || *t == *with_value)
    };
    for flag in hidden_flags(&command) {
        let upper = flag.to_ascii_uppercase();
        // A name in another case is written as declared only when the pre-pass
        // resolves it; any other name stays as typed.
        let typed = prepass(vec!["vep".into(), format!("--{upper}").into()], [home()]).unwrap();
        assert!(names(typed, flag), "--{upper} on the command line");

        let config = dir.path().join(format!("{flag}.txt"));
        fs::write(&config, format!("{flag} 1\n")).unwrap();
        let from_file = prepass(
            vec!["vep".into(), "--config".into(), config.into_os_string()],
            [home()],
        )
        .unwrap();
        assert!(names(from_file, flag), "{flag} in a configuration file");

        let variable = format!("{ENV_PREFIX}{}", flag.to_ascii_uppercase());
        let from_environment = prepass(
            vec!["vep".into()],
            [home(), (OsString::from(&variable), OsString::from("1"))],
        )
        .unwrap();
        assert!(names(from_environment, flag), "{variable}");
    }
}

/// `command` with one more switch, named `long`, built.
fn with_switch(command: &Command, long: &'static str) -> Command {
    let mut added = command
        .clone()
        .arg(Arg::new(long).long(long).action(ArgAction::SetTrue));
    added.build();
    added
}

/// The three pages that offer abbreviated option names each say once, where
/// they offer them, what an abbreviation rests on: the compatibility page in
/// its list of option spellings, the companion tools' pages in the paragraph
/// after the opening one. No other page carries either text.
#[test]
fn the_pages_that_offer_abbreviations_say_what_they_rest_on() {
    let mut command = Args::command_with_error_table();
    command.build();
    let pages: BTreeMap<&str, String> = render(&command).into_iter().collect();
    for (file, text) in &pages {
        let of_vep = usize::from(*file == COMPATIBILITY);
        let of_a_tool = usize::from([FILTER_PAGE, RECODE_PAGE].contains(file));
        assert_eq!(text.matches(VEP_ABBREVIATIONS).count(), of_vep, "{file}");
        assert_eq!(
            text.matches(TOOL_ABBREVIATIONS).count(),
            of_a_tool,
            "{file}"
        );
    }
    let spellings = pages[COMPATIBILITY]
        .split("\n## ")
        .find(|section| section.starts_with("Option spellings"))
        .expect("the compatibility page has a section on option spellings");
    assert!(spellings.contains("abbreviated to any prefix that names one option"));
    assert!(spellings.contains(VEP_ABBREVIATIONS), "{spellings}");
    for file in [FILTER_PAGE, RECODE_PAGE] {
        let mut paragraphs = pages[file].split("\n\n").skip(1);
        let opening = paragraphs.next().unwrap_or_default();
        assert!(
            opening.contains("abbreviated to any prefix that names one option"),
            "{file}: {opening}"
        );
        assert_eq!(
            paragraphs.next(),
            Some(TOOL_ABBREVIATIONS.trim_end()),
            "{file}"
        );
    }
}

/// What the abbreviation sentences say of the readers. `vep` resolves a
/// prefix against the names its parser declares, aliases and hidden flags
/// among them, then against Ensembl VEP's list: one more flag under the same
/// prefix leaves a prefix that list does not resolve to no flag, a flag named
/// as the prefix takes it, and a prefix that list resolves keeps its flag
/// beside one more flag under it. A companion tool resolves against its own
/// options alone: one more option under the same prefix leaves the name,
/// written with two dashes after one dash as after two, for clap to refuse
/// with exit 2, and an option named as the prefix takes it.
#[test]
fn an_added_option_takes_an_abbreviation_or_leaves_it_to_none() {
    let mut command = Args::command_with_error_table();
    command.build();
    let perl = perl_options();
    let reached = |command: &Command, prefix: &str| {
        let known = parser_options(command);
        resolve_spelling(prefix, &known, &perl).map(|option| option.canonical.clone())
    };
    let known = parser_options(&command);
    for arg in command.get_arguments() {
        let Some(long) = arg.get_long() else {
            continue;
        };
        let names = &known.iter().find(|o| o.canonical == long).unwrap().names;
        for alias in arg.get_all_aliases().into_iter().flatten() {
            assert!(names.iter().any(|n| n == alias), "--{long}: {alias}");
        }
    }
    for flag in hidden_flags(&command) {
        assert_eq!(reached(&command, flag).as_deref(), Some(flag));
    }

    assert_eq!(
        reached(&command, "parquet").as_deref(),
        Some("parquet_partition")
    );
    assert!(resolve("parquet", &perl).is_none());
    let beside = with_switch(&command, "parquet_other");
    assert_eq!(reached(&beside, "parquet"), None);
    let named = with_switch(&command, "parquet");
    assert_eq!(reached(&named, "parquet").as_deref(), Some("parquet"));
    assert_eq!(reached(&command, "bio").as_deref(), Some("biotype"));
    let beside = with_switch(&command, "biotype_other");
    assert_eq!(reached(&beside, "bio").as_deref(), Some("biotype"));

    for (tool, prefix, option, beside, program) in [
        (
            FilterArgs::command(),
            "lim",
            "limit",
            "lim_other",
            "vep filter",
        ),
        (
            RecodeArgs::command(),
            "mane",
            "mane_select",
            "mane_other",
            "vep recode",
        ),
    ] {
        let mut tool = tool;
        tool.build();
        let typed = [OsString::from(format!("--{prefix}"))];
        let spelled = |command: &Command| spell_as_declared(&typed, command);
        assert_eq!(spelled(&tool), [OsString::from(format!("--{option}"))]);

        let tied = with_switch(&tool, beside);
        assert_eq!(spelled(&tied), typed, "{program}");
        let one_dash = [OsString::from(format!("-{prefix}"))];
        assert_eq!(spell_as_declared(&one_dash, &tied), typed, "{program}");
        let refusal = tied
            .clone()
            .try_get_matches_from(std::iter::once(OsString::from(program)).chain(spelled(&tied)))
            .unwrap_err();
        assert_eq!(refusal.kind(), clap::error::ErrorKind::UnknownArgument);
        assert_eq!(refusal.exit_code(), 2, "{program}");

        let named = with_switch(&tool, prefix);
        let taken = resolve(prefix, &parser_options(&named)).map(|o| o.canonical.clone());
        assert_eq!(taken.as_deref(), Some(prefix), "{program}");
    }
}
