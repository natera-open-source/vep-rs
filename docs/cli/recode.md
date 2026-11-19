# vep recode

`vep recode` translates a variant between its identifiers and notations, as Ensembl VEP's `variant_recoder` does, from the cache alone. The input is an identifier, an HGVS, SPDI, VCF, region or Ensembl-format record, one per line, from `--input_file` or `--input_data` (`--id`). It takes the same options, read as Getopt::Long reads them: after one dash or two, in any case, and abbreviated to any prefix that names one option. A repeated option keeps its last value, and every `--transcript_filter` given is kept. The short forms are `-i`, `-s`, `-a` and `-u`.

`--config` names a file of `<name> <value>` lines. An option is read from the four sources `vep` reads its own from, the first that names it winning: the command line, the `--config` file, `<dir>/vep.ini`, where `<dir>` is the file's `dir` or else `~/.vep`, then the environment variable `VEP_<OPTION>`. A repeatable option collects every source's values in that order, and a name no option of this page has is dropped.

The options `variant_recoder` settles before it reads the sources are set on the command line alone, as there. `fields`, `buffer_size`, `failed`, `merged`, `lrg` and `ambiguous_hgvs` keep the values it presets. `vcf_string`, `var_synonyms`, `mane_select` and `ga4gh_vrs` add their output fields before the sources are read, so a source's value adds none. `pretty` and `vcf_string` are names its reader does not know. Every other option of this page but `--config` may be set from a source. A `vep` option this page lacks is dropped from every source, where `variant_recoder` reads it and applies it ([the class](../divergences/recode-config-file.md)).

```bash
# Every notation of an identifier
vep recode --id rs699 --json_cache /path/to/cache

# A protein notation, with its VCF form and its MANE Select transcript
vep recode --id "ENSP00000284967.6:p.Ala331Thr" --vcf_string --mane_select \
  --pretty --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38

# Two fields alone, from a file of SPDI notations
vep recode -i variants.spdi.txt --fields hgvsg,hgvsc --json_cache /path/to/cache
```

The output is one JSON array on standard output, one object per input line; `--pretty` indents it. Each object is keyed by the input's alleles. An allele's object holds `input` and, for each field named, the unique values found in that allele's consequence entries. The default fields are `id`, `hgvsg`, `hgvsc`, `hgvsp` and `spdi`; `--fields` may name any key of the JSON output instead (`transcript_id`, `consequence_terms`).

`--vcf_string` adds the allele's VCF form, `chromosome-position-reference-alternate`. `--var_synonyms` adds the known synonyms of the variant, one `source: id,id` string per source. `--mane_select` adds the `hgvsg`, `hgvsc` and `hgvsp` triple of each MANE Select transcript's entry and switches those three fields on. `--ga4gh_vrs` adds a GA4GH VRS Allele object for each `spdi`, which it switches on.

The translation is an annotation run under a fixed configuration: `--json --hgvsg_use_accession --distance 0`, plus the HGVS and SPDI flags the fields need. HGVSg and SPDI name a chromosome by its `NC_` accession, and a scaffold or patch by its RefSeq name, where the cache's synonyms have one. The VCF form is built from the record's alleles for every kind of input, the base before an indel read from the reference, as `to_VCF_record` builds it. An input whose record has no consequence entry is `null`.

The identifiers are the co-located variants' names per allele, a COSMIC or HGMD identifier listed under every allele. A known variant is co-located as the recoder's database query co-locates it: every row of the record's region whose span overlaps the record and whose alleles match, however the two are spelled, the variants failed in dbSNP included (`--failed 1`). So a deletion written from the last base of a repeat finds the row dbSNP writes from the first base, and an insertion finds a row only when it spans both of the bases around it.

For an identifier input whose stored alleles are an unminimised indel, such as `GGG/GG`, `variant_recoder` writes the transcript and protein notations and the MANE Select triple without shifting the indel 3' along the transcript. It shifts them for every other input of the same variant. The stored class of the database row, `indel`, outlives the minimisation of its alleles, and `_return_3prime` shifts an insertion or a deletion alone.

`vep recode` shifts them for every input, so those notations of such an identifier differ from the recoder's by the length of the repeat. The comparator's class `identifier_indel_transcript_notations` leaves them out of a comparison.

In `variant_recoder`, a multi-allelic known variant's identifier reaches an allele of the record that spells none of its alternates only through the minimised form of its last alternate, so a deletion allele written on its own gets no identifier from it. `vep recode` minimises each alternate and names the variant on the allele it spells.

The warnings raised while an input is read or annotated are carried in the output, under `warnings` on the next object written. `--buffer_size`, 1 by default, is the number of inputs read before an object is written, as in `variant_recoder`. A run that writes no object is one object of its warnings alone, and `--warning_file` writes them to a file as well.

The synonyms and the co-located variants are the cache's variation data, where `variant_recoder` reads the Ensembl database. The cache's synonyms are the ClinVar, OMIM, UniProt, PharmGKB, PhenCode and COSMIC identifiers and never the merged dbSNP identifiers (`Archive dbSNP`) the database adds; the comparator's class `archive_dbsnp_synonyms` leaves those out. `variant_recoder --var_synonyms` stops with an error on a VCF input and on the identifier of a failed variant, where `vep recode` writes the cache's synonyms of the variant the ID column or the identifier names.

A RefSeq notation (`NM_`, `NP_`, `XM_`, `XP_`) appears only with a merged cache. An LRG notation never appears, as `--lrg` is refused, and a protein notation resolves to one variant, as `--ambiguous_hgvs` is refused. A ClinGen allele identifier (`CA...`) is not an input. `ga4gh_vrs` holds a GA4GH VRS Allele object only for a chromosome whose synonyms name an `NC_` accession, and is an empty list for any other. The database connection flags (`--host`, `--port`, `--user`, `--password`, `--registry`, `--db_version`, `--genomes`, `--is_multispecies`) and `--bam` are refused with an error naming the flag. `--mane_select` is for human alone.

`variant_recoder` reads the merged transcript set on every run, and `vep recode` reads it wherever the cache named holds it, with `--merged` or without. Such a cache is a `--json_cache` directory whose `info.json` declares the merged set, or a JSON cache in the `_merged` species directory under the cache root. Beside the Ensembl set's cache the run would read, it is the cache of the same `<cache_version>_<assembly>` directory name, so the run keeps its cache version and assembly; under a root with no directory of the Ensembl set, it is the one `--merged` resolves.

An Ensembl VEP Storable cache is not one, and a cache of the Ensembl set alone is read as it is. Where the directory of the Ensembl set is there and holds no cache the run can read, the run is refused as it is without a merged cache, and `--merged` reads the merged one.

`--refseq` is refused with `Can't use --refseq and --merged together`, as `variant_recoder` refuses it, and with `--gencode_basic` the error names that flag instead. A RefSeq cache, which is read under `--refseq` alone, is not an input; a merged cache carries the RefSeq models beside the Ensembl ones.

Every error is one line `ERROR: [<component>] <message>` on standard error and exit 1, the component naming what refused the run (`recode`, `config`, `input`, `cache`, `fasta` or `filter`); a command-line syntax error is exit 2.

## Input

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--config` | String | none | A configuration file, one option per line as \<name\> \<value\> without the dashes (cache\_version 116, pick 1); an option the command line also gives keeps the command line's value |
| `--input_file`, `-i` | String | none | The input file, one variant per line; a gzip or BGZF file is decompressed |
| `--input_data`, `--id` | String | none | The input as a string, read in place of `--input_file`; one variant per line |
| `--format` | String | none | The input format: ensembl, vcf, hgvs, id, spdi or region, detected from the first line when absent |
| `--delimiter` | String | none | The column delimiter of an Ensembl-format or VCF input |

## Output

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--pretty` | bool | false | Indent the JSON output, three spaces per level, as Ensembl VEP prints it |
| `--fields` | String (comma-separated) | none | The fields of each allele's object, comma-separated; default id,hgvsg,hgvsc,hgvsp,spdi. A field named here switches on what computes it: hgvsg, hgvsc, hgvsp and spdi need the FASTA, id the co-located variants. Any key of the annotation's JSON output may be named |
| `--vcf_string` | bool | false | Add vcf\_string: each allele's VCF form as chromosome-position-reference-alternate |
| `--var_synonyms` | bool | false | Add var\_synonyms: the known synonyms of the variant, one "source: id,id" string per source, from the cache's co-located variants |
| `--mane_select` | bool | false | Add mane\_select: the hgvsg, hgvsc and hgvsp of the MANE Select transcript, which switches those three fields on; human only |
| `--ga4gh_vrs` | bool | false | Add ga4gh\_vrs: a GA4GH VRS Allele object for each spdi, which it switches on; an empty list for a chromosome whose synonyms name no NC\_ accession |
| `--warning_file` | String | none | Also write the warnings, which the output carries, to this file |
| `--debug` | bool | false | Accepted as variant\_recoder accepts it: there it prints the plugin directories and modules loaded, of which vep recode has none, so the output is unchanged |

## Cache, species and assembly

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--json_cache` | String | none | The cache directory itself, the one holding info.json and transcripts/, read in place of the `--dir_cache` layout |
| `--dir_cache` | String | none | The cache root: the cache read is \<dir\_cache\>/\<species\>/\<cache\_version\>\_\<assembly\>, the species directory suffixed \_merged when the merged set is read, as Ensembl VEP lays it out |
| `--cache_version` | u32 | none | The cache version, the prefix of the cache directory under the species directory |
| `--fasta` | String | none | The reference FASTA, needed for the HGVS and SPDI fields; a FASTA found in the cache directory is read without the flag |
| `--synonyms` | String | none | A chromosome synonyms file, read in place of the cache's own |
| `--species`, `-s` | String | none | The species, as the cache directory names it; default homo\_sapiens |
| `--assembly`, `-a` | String | none | The assembly, checked against the cache's |
| `--grch37` | bool | false | Select the GRCh37 assembly |

## Transcript selection

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--refseq` | bool | false | Refused: the recoder reads the merged transcript set on every run, so the flag meets "Can't use `--refseq` and `--merged` together" |
| `--merged` | bool | false | Read the merged transcript set, Ensembl and RefSeq models together, the species directory suffixed \_merged; a cache that holds the set is read without the flag, as variant\_recoder presets it |
| `--all_refseq` | bool | false | Keep every model of a merged cache, not the RefSeq accessions alone |
| `--gencode_basic` | bool | false | Use the GENCODE basic transcript set alone |
| `--transcript_filter` | String (repeatable) | none | A transcript filter, as vep's `--transcript_filter` takes it; repeatable |
| `--exclude_predicted` | bool | false | Drop RefSeq's predicted transcript models, the XM\_ and XR\_ accessions |
| `--minimal` | bool | false | Convert each input allele to its minimal representation before annotation |
| `--pick` | bool | false | Keep one consequence per variant |
| `--pick_allele` | bool | false | Keep one consequence per variant allele |
| `--per_gene` | bool | false | Keep one consequence per gene |
| `--pick_allele_gene` | bool | false | Keep one consequence per allele per gene |
| `--pick_order` | String (comma-separated) | none | The criteria the pick flags compare, in order, comma-separated |
| `--buffer_size` | usize | none | The number of variants annotated per batch; default 1 |
| `--failed` | i64 | none | The highest failed value a known variant may carry and still be reported; default 1 |
| `--shift_hgvs` | String | none | Describe an insertion or deletion at its most 3' position in HGVS (1, the default) or where the input placed it (0) |
| `--use_transcript_ref` | bool | false | Read each allele's reference from the transcript instead of the input |
| `--use_given_ref` | bool | false | Keep the input's reference allele on a cache that would otherwise switch `--use_transcript_ref` on |

## Database

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--registry` | String | none | Refused: variant\_recoder reads an Ensembl database; vep recode reads a local cache alone |
| `--host` | String | none | Refused: the database host |
| `--port` | String | none | Refused: the database port |
| `--user`, `-u` | String | none | Refused: the database user |
| `--password`, `--pass` | String | none | Refused: the database password |
| `--db_version` | String | none | Refused: the Ensembl database version; `--cache_version` selects a cache |
| `--genomes` | bool | false | Refused: the Ensembl Genomes database host |
| `--is_multispecies` | String | none | Refused: a multi-species Ensembl database |
| `--lrg` | bool | false | Refused: LRG mapping needs the Ensembl core database |
| `--ambiguous_hgvs` | bool | false | Refused: resolving a protein notation to every matching variant needs the Ensembl database |
| `--bam` | String | none | Refused: an alignment file that would edit the cache's transcript models |

