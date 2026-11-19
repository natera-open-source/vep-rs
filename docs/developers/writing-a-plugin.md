# Writing a plugin

Write a vep-rs plugin as a built-in compiled into the binary, the route for
canonical Ensembl VEP plugins, or as a shared library that ships separately.

This page describes:

- the two plugin systems
- the C ABI a shared library exports and the JSON it receives
- how a built-in plugin is added to the `vep-builtins` crate
- how a plugin's fields are compared against Ensembl VEP's

## Architecture

### Built-in vs dynamic plugins

vep-rs has two plugin systems:

**Built-in plugins** (`vep-builtins` crate) are compiled into the binary. They
operate on native Rust types with no serialization overhead. This is the
recommended approach for canonical VEP plugins.

**Dynamic plugins** (`vep-plugin` crate) are compiled as separate `cdylib`
shared libraries loaded at runtime via C FFI. Communication uses JSON
serialization. This is the route for extensions that ship separately from
vep-rs. The host side is `crates/vep-cli/src/dylib_plugins.rs`.

### Worked example

`crates/vep-plugin/examples/echo_plugin.rs` is a complete plugin. It echoes the
consequence's `transcript_id` as `ECHO_FEATURE` and its first parameter, when
one was given, as `ECHO_PARAM`. It declares both in its header info and lists
`Transcript` as its one feature type. It reads the one field it needs with a
string scan, so it adds no dependency.

Build it as the shared library the loader opens, then pass it by path or by
name. In the commands, `input.vcf` is a VCF of your own and
`/path/to/json_cache` is the JSON cache directory that
[Cache setup](../cache-setup.md) describes:

```bash
cargo build --example echo_plugin -p vep-plugin
# -> target/debug/examples/libecho_plugin.so (Linux)
#    or libecho_plugin.dylib (macOS)

vep -i input.vcf -o output.txt --offline --cache \
    --json_cache /path/to/json_cache --assembly GRCh38 \
    --plugin target/debug/examples/libecho_plugin.so,hello

vep -i input.vcf -o output.txt --offline --cache \
    --json_cache /path/to/json_cache --assembly GRCh38 \
    --dir_plugins target/debug/examples --plugin echo_plugin,hello
```

Every transcript line then carries this in its Extra column, and the header
describes the two fields after the flag fields' lines:

```text
## ECHO_FEATURE : The transcript id the consequence names
## ECHO_PARAM : The first plugin parameter
...
ECHO_FEATURE=<transcript_id>;ECHO_PARAM=hello
```

Under `--tab`, `--vcf` or `--output_format parquet` the two fields are the last
two columns. An intergenic line carries neither field, since the plugin listed
`Transcript` alone.

`crates/vep-plugin/examples/echo_all_rows_plugin.rs` is the other shape: its
`vep_plugin_feature_types` lists all four kinds, so it runs on every row. It
writes the kind of row it was given as `ECHO_ROW`: the consequence's
`feature_type`, or `intergenic` where the consequence document is `null`.

`crates/vep-plugin/examples/faulty_plugin.rs` fails on purpose: its first `run`
returns `[1]`, its second NULL, and every later call `{"FAULTY_OK":"1"}`. A run
with it loaded writes the two warnings of
[a failed `run` call](#when-a-run-call-fails) once each and annotates every
other transcript row. That is how the host's handling of a defective plugin is
tested.

`cargo test -p vep-plugin` builds the three examples and runs
`tests/echo_roundtrip.rs`. That test loads the echo plugins through the same
loader, by path and by name, and reads their header info and feature types.

### How a plugin is found

When `--plugin Name,...` is specified, the runner checks the built-in registry
first. A name no built-in claims is resolved as a shared library.

A token written as a file, one holding a path separator or a suffix
(`./libecho.so`, `libecho.so`), is the library file. A bare file name is opened
from the working directory. A bare name is the library's stem, found as
`libName.so` and then `libName.dylib` under `--dir_plugins`, each a regular
file. The working directory is not searched for it, as Ensembl VEP's `use`
takes only a `Name.pm` module from its library path. These flags load
`/opt/vep/plugins/libecho_plugin.so`:

```bash
--dir_plugins /opt/vep/plugins --plugin echo_plugin
```

The string a plugin returns from `vep_plugin_name` appears only in log and error
messages.

`--dir_plugins` defaults to `<--dir>/Plugins`, as in Ensembl VEP. Its value is
expanded as Ensembl VEP's `glob` expands it. `*`, `?`, `[...]`, `{a,b}`, `~`
and `~user` have their shell meanings, and a backslash quotes the character
after it.

Every match that is a directory is searched, the last match first, and the
first library found is the plugin. When `--plugin` names a plugin and a match
is not a directory, the run writes this message to the warnings file. Under
`--safe` it is an error:

```text
Plugins directory '<dir>' not found.
```

A plugin receives its `--plugin` parameters and nothing else of the run. Ensembl
VEP hands every plugin the whole configuration hash as well (`BaseVepPlugin.pm`
`config`, lines 155-168), which a plugin reads to learn the assembly, the cache
directory or any other flag. The ABI has no such call, so a shared library that
needs a run setting takes it as a parameter.

### What a failed load writes

A name that resolves to neither a built-in nor a library is settled in this
order, the first three rows of the table:

| Case | What the run writes | Without `--safe` | Under `--safe` |
| --- | --- | --- | --- |
| An Ensembl VEP plugin vep-rs does not provide | `ERROR: [plugin] Name is an Ensembl VEP plugin vep-rs does not provide: <reason>`, a line that points at the [support table](../plugins.md) | The run ends | The run ends |
| A Perl module, refused by name whether given as `Name.pm` by path or as `Name` with `Name.pm` under the plugins directory | `ERROR: [plugin] <path>/Name.pm is an Ensembl VEP Perl plugin; vep-rs runs no Perl plugin, see docs/plugins.md` | The run ends | The run ends |
| Any other name | `Failed to load plugin Name: Name is not a built-in plugin and no libName.so or libName.dylib is under <dir>`, in the warnings file | The run continues without the plugin, which is Ensembl VEP's default | An error, as Ensembl VEP's is |
| A file that exists but fails a load check (listed below) | `Failed to compile plugin Name: <reason>`, in the warnings file | The run continues without it, as Ensembl VEP's does when a plugin module fails to compile | An error, as Ensembl VEP's is |
| A plugin whose `vep_plugin_create` returns NULL or whose `vep_plugin_init` returns non-zero; a built-in plugin whose parameters it cannot use (a scores file that does not open) | `Failed to instantiate plugin Name: <reason>`, in the warnings file | The run continues without it, as Ensembl VEP's does when a plugin's constructor dies | An error, as Ensembl VEP's is |

A file fails a load check, and is dropped with a warning, when it:

- does not open as a shared library
- lacks one of the nine required symbols
- reports an ABI version other than 2
- returns header info or feature types the host cannot read (below)

A library under `--dir_plugins` loads before the support table is read, so a
plugin you build under an Ensembl VEP plugin's name is used.

### The C ABI (version 2)

A plugin exports these ten symbols with these signatures
(`crates/vep-plugin/src/abi.rs`); nine are required and `vep_plugin_feature_types`
is optional.

#### `vep_plugin_abi_version`

```rust
extern "C" fn() -> u32
```

Returns `2`; any other value is a load error. A library built against version 1
is refused with `incompatible plugin ABI version: expected 2, got 1`. Version 1
had no header symbol and skipped `init` when no parameter was given.

#### `vep_plugin_create`

```rust
extern "C" fn() -> *mut c_void
```

Returns an opaque instance pointer; NULL is a failed instantiation (above).

#### `vep_plugin_destroy`

```rust
extern "C" fn(*mut c_void)
```

Frees the instance; called when the host drops the plugin.

#### `vep_plugin_name`

```rust
extern "C" fn(*const c_void) -> *const c_char
```

Returns a NUL-terminated name; the pointer must stay valid for the instance's
lifetime.

#### `vep_plugin_version`

```rust
extern "C" fn(*const c_void) -> *const c_char
```

Returns a NUL-terminated version string, under the same lifetime rule.

#### `vep_plugin_init`

```rust
extern "C" fn(*mut c_void, *const *const c_char, usize) -> i32
```

Receives the `--plugin` parameters after the name as NUL-terminated strings;
returns 0 on success. Called exactly once per instance, on every load, with a
count of 0 when no parameter was given. The call comes before
`vep_plugin_header_info`, `vep_plugin_feature_types` and any `vep_plugin_run`.
The fields a plugin declares may therefore depend on its parameters, as Ensembl
VEP's `get_header_info` may (`BaseVepPlugin.pm` `new`, which receives the
parameters, runs before it).

#### `vep_plugin_header_info`

```rust
extern "C" fn(*const c_void) -> *mut c_char
```

Returns a NUL-terminated JSON array of `[key, description]` pairs: the fields the
plugin writes and the description the headers print for each. It is `[]` or
`null` for none. The host frees the string with `vep_plugin_free`. A NULL
return or text that is not such an array is a load error.

The keys may come in any order. The host sorts each plugin's keys and lists
plugins in load order, as Ensembl VEP's `get_plugin_headers` does
(`OutputFactory.pm` lines 2407-2421). The built-in plugins come before the
shared libraries whatever the `--plugin` order; Ensembl VEP, with one kind of
plugin, follows the command line throughout.

#### `vep_plugin_feature_types` (optional)

```rust
extern "C" fn(*const c_void) -> *mut c_char
```

Returns a NUL-terminated JSON array of the kinds of row the plugin runs on. The
kinds are `"Transcript"`, `"RegulatoryFeature"`, `"MotifFeature"` and
`"Intergenic"`, the values Ensembl VEP's `feature_types` takes
(`BaseVepPlugin.pm` lines 220-223 and 229-233). The host frees the string with
`vep_plugin_free`; any other name is a load error.

A library that does not export the symbol runs on transcript rows alone,
Ensembl VEP's default for a plugin that does not override `feature_types`
(`BaseVepPlugin.pm` line 99); the `echo_all_rows_plugin` example lists all four
kinds, so it runs on every row, the intergenic rows included.

A plugin that lists `"RegulatoryFeature"` or `"MotifFeature"` switches
regulatory annotation on for the run, with the status line
`Fetching regulatory features for plugin: <name>`. Ensembl VEP does the same
for such a plugin (`Runner.pm` lines 998 to 1002).

A value a plugin writes prints on the rows of the kinds it lists. A plugin
listing `"Intergenic"` writes on a record's intergenic row. Its field stays off
the record's regulatory and motif rows, a field its header info never declared
included.

The feature type is the one interest check a shared library has. Ensembl VEP
also asks each plugin whether it takes the variant's kind, `VariationFeature` or
`StructuralVariationFeature` (`check_variant_feature_type`, `OutputFactory.pm`
line 2220). The ABI has no such declaration, so a shared library sees sequence
and structural variants alike and reads `is_structural` in the variant document
when it wants one kind.

#### `vep_plugin_run`

```rust
extern "C" fn(*const c_void, *const c_char, *const c_char) -> *mut c_char
```

Takes `(instance, consequence JSON, variant JSON)` and returns a NUL-terminated
JSON object that the host frees with `vep_plugin_free`. On an intergenic row the
consequence JSON is `null`. NULL signals an error.

#### `vep_plugin_free`

```rust
extern "C" fn(*mut c_char)
```

Frees a string returned by `vep_plugin_run`, `vep_plugin_header_info` or
`vep_plugin_feature_types`.

### What `run` receives

After the built-in plugins have run on a batch, the host calls `vep_plugin_run`
once per row a plugin runs on. It does so for every loaded shared library in
`--plugin` order. The rows are each entry of a variant's
`transcript_consequences` whose `feature_type` the plugin listed (the
`Transcript` entries when it exports no `vep_plugin_feature_types`).

A plugin that listed `Intergenic` also runs on the intergenic row of a variant
overlapping no transcript, with or without regulatory or motif entries. There
the consequence document is `null`. A plugin that listed `Transcript` alone
never sees a regulatory, motif or intergenic row. This is the
`check_feature_type` gate of `OutputFactory.pm` `run_plugins` (lines
2220-2223), where a feature-less row is checked as `Intergenic`.

The two JSON documents are the serde serialization of `TranscriptConsequence`
and `InputVariant`.

#### The consequence document

The consequence document is written by the `Serialize`
implementation in `crates/vep-core/src/consequence.rs`. Its keys are the
`TRANSCRIPT_CONSEQUENCE_KEYS` table there, in that order, each marked whether it
is written only when populated; the test
`serialize_emits_the_derive_key_set_in_order` pins the list. A consequence holds
a shared handle on its transcript and reads the transcript's values through
accessors of the key's name, so the document is flat. These have no key:

- the handle itself
- `picked`
- `protein_id_display`
- `sift_term`
- `polyphen_term`
- `protein_start`
- `hgvs_intron_offsets`
- `exon_ids`
- `transcript_sequences`
- `id_rank`
- `record_span_in_cds`

Every field a plugin will find:

| Field                                              | Meaning                                                                                      | Example                                            |
| -------------------------------------------------- | -------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| `transcript_id`                                    | Transcript stable ID                                                                         | `"ENST00000366667"`                                |
| `gene_id`                                          | Gene stable ID                                                                               | `"ENSG00000142168"`                                |
| `gene_symbol`, `gene_symbol_source`                | Gene symbol and its source                                                                   | `"SOD1"`, `"HGNC"`                                 |
| `hgnc_id`                                          | HGNC ID                                                                                      |                                                    |
| `consequences`                                     | Array of the `Consequence` enum's variant names, not the SO terms                            | `["MissenseVariant"]`, `["StopGained"]`            |
| `impact`                                           | Impact of the most severe consequence                                                        | `"HIGH"`, `"MODERATE"`, `"LOW"`, `"MODIFIER"`      |
| `biotype`                                          | Transcript biotype                                                                           | `"protein_coding"`                                 |
| `canonical`                                        | Whether the transcript is canonical                                                          | `true`                                             |
| `cdna_position`, `cds_position`, `protein_position` | Positions, as strings                                                                        | `"428"`, `"301-303"`                               |
| `cdna_total`, `cds_total`, `protein_total`         | The lengths `--total_length` appends to those positions; see below                           |                                                    |
| `own_span_positions`                               | On an allele of a `--minimal` split, the three positions of the allele's own trimmed span, where the columns above carry the record span's; present only then | `{"cdna": "439", "cds": "87", "protein": "29"}` |
| `amino_acids`                                      | Amino acid change                                                                            | `"A/V"`                                            |
| `codons`                                           | Codon change                                                                                 | `"gCc/gTc"`                                        |
| `protein_id`                                       | Ensembl protein ID                                                                           |                                                    |
| `distance`                                         | Distance to the transcript, for upstream and downstream consequences                         | `4123`                                             |
| `strand`                                           | Transcript strand                                                                            | `1`, `-1`                                          |
| `exon`, `intron`                                   | Exon or intron number over the total; a variant overlapping two or more features carries the range `first-last/total` | `"3/5"`, `"2/4"`, `"6-7/12"`                       |
| `hgvsc`, `hgvsp`                                   | HGVS coding and protein notations, as the output prints them before escaping; `--hgvsp_use_prediction` and `--remove_hgvsp_version` shape `hgvsp` here too |  |
| `hgvs_offset`                                      | The signed number of bases the HGVS notations were shifted along the transcript; `null` when nothing shifted |                                   |
| `shift_length`                                     | The bases `--shift_3prime 1` moved the allele along the transcript; `0` otherwise            | `0`                                                |
| `shifted_allele`                                   | Under `--shift_3prime 1`, the inserted bases rotated along the shift, the `Allele` a coding row prints; present only on such a row | `"GAA"`                                            |
| `used_ref`                                         | Under `--use_transcript_ref`, the reference allele the row was analysed against, on the variant's strand: the transcript's inside its exons, the input's elsewhere; present only on a sequence variant's row under the flag | `"A"`                             |
| `sift`, `polyphen`                                 | Prediction and score                                                                         | `"deleterious(0.01)"`                              |
| `domains`                                          | Array of `[source, id]` pairs                                                                | `[["Pfam", "PF00080"]]`                            |
| `feature_type`                                     | The feature kind                                                                             | `"Transcript"`, `"RegulatoryFeature"`, `"MotifFeature"` |
| `flags`                                            | Transcript flags, then `PICK` on a row `--flag_pick` picked                                  | `["cds_start_NF"]`                                 |
| `tsl`                                              | Transcript Support Level                                                                     | `1`                                                |
| `mane_select`, `mane_plus_clinical`                | MANE transcript IDs, when the transcript carries the tag                                     |                                                    |
| `appris`                                           | APPRIS annotation                                                                            | `"principal1"`                                     |
| `ccds`                                             | CCDS cross-reference                                                                         |                                                    |
| `gencode_primary`                                  | `true` on a transcript of the GENCODE primary set                                            | `true`                                             |
| `swissprot`, `trembl`, `uniparc`, `uniprot_isoform` | UniProt cross-references: Swiss-Prot, TrEMBL, UniParc and isoform accessions                 |                                                    |
| `refseq`                                           | RefSeq cross-reference                                                                       |                                                    |
| `transcript_version`, `gene_version`               | The transcript's and the gene's stable-id versions                                           |                                                    |
| `gene_phenotype`                                   | `true` when the gene has a phenotype association                                             | `true`                                             |
| `mirna`                                            | The `miRNA_stem` and `miRNA_loop` terms the variant overlaps on a miRNA transcript; empty elsewhere | `["miRNA_stem"]`                            |
| `cdna_length`                                      | The transcript's spliced length                                                              |                                                    |
| `cds_length`                                       | The translateable sequence's length; present for a translated transcript                     |                                                    |
| `source_cache`                                     | `"Ensembl"` or `"RefSeq"`; present for a merged cache's transcript                           | `"Ensembl"`                                        |
| `bam_edit_status`                                  | `"ok"` or `"failed"`; present for a transcript of a BAM-edited cache whose model the alignment corrected or failed to | `"ok"`                                   |
| `feature_start`, `feature_end`                     | Genomic span of the feature, 1-based inclusive; both `0` when unknown                        |                                                    |
| `plugin_data`                                      | The fields the built-in plugins have already written to this consequence; present only when non-empty. `SpliceVault_top_events` holds its events `\|`-separated, as Ensembl VEP's plugin splits the table's events column (`SpliceVault.pm` line 216), and the text formats join them with `,` or `&` | `{"REVEL": "0.931"}`                      |
| `loftee_ctx`                                       | The exon/intron model; present when a built-in that reads it is loaded (LoFTEE, MaxEntScan, NMD, NearestExonJB on every row, UTRAnnotator on `5_prime_UTR_variant` rows) |                                                    |
| `motif`                                            | A motif row's `name` (the binding matrix's stable id), `pos` (the variant's position in the matrix, absent outside it) and `factors`; present on a sequence variant's `MotifFeature` row alone | `{"name": "ENSPFM0290", "pos": 7, "factors": ["HNF4A"]}` |
| `cell_type`                                        | Under `--cell_type`, the `CELL_TYPE` list of a `RegulatoryFeature` or `MotifFeature` row: one object per listed epigenome the feature has an activity for, `name` and `activity` (`ACTIVE`, `INACTIVE`, `POISED`, `REPRESSED` or `NA`); empty on a motif row; present on those rows alone | `[{"name": "HUVEC", "activity": "REPRESSED"}]` |

`cdna_total`, `cds_total` and `protein_total` are the lengths `--total_length`
appends to the three positions. Each is present only when its rule holds,
whatever the flag. The cDNA total is present wherever the variant overlaps an
exon. The coding and peptide totals are present where the variant lies in the
coding region and the cache carries the sequence.

#### The variant document

The variant document is the serde derive of `InputVariant`
(`crates/vep-core/src/variant.rs`), where a JSON key is the Rust field name. A
field marked `#[serde(skip)]` (`raw_input_bytes`, `split_first_feature_rowless`,
`reference_row_shift_allele`) never appears, and a field marked
`skip_serializing_if` is left out when it is `None`, empty or zero. In both
documents every other absent optional value is `null`.

The sequence-variant fields:

| Field                                                                 | Meaning                                                                                     | Example                                                   |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `chr`                                                                 | Chromosome, normalized                                                                      | `"21"`                                                    |
| `original_chr`                                                        | The name the record prints: as input, or `"MT"` for a bare `M` that the annotation sources (the cache and any `--custom` or `--regulatory_gff` file) know only as `MT` | `"chr21"`, `"MT"`                                         |
| `input_chr`                                                           | The input's `"M"` behind such a rename; absent otherwise                                    | `"M"`                                                     |
| `start`, `end`                                                        | 1-based inclusive; `end` is `start - 1` for an insertion                                    | 9068417                                                 |
| `strand`                                                              | Strand                                                                                      | `"Forward"`, `"Reverse"`                                  |
| `allele_string`                                                       | Alleles in VEP form                                                                         | `"A/G"`, `"AC/-"`, `"-/TT"`                               |
| `written_case_allele_string`                                          | The record's allele string after the anchor chop and any `--minimal` trim, in the case the input wrote it; present only when a base is lower case | `"ca/-"`                                     |
| `ref_allele`, `alt_alleles`                                           | Serialized as arrays of byte values, so read `allele_string` for the sequences              | `[65]`                                                    |
| `allele_index`                                                        | Which ALT of the record this allele is                                                      | `0`                                                       |
| `allele_num`                                                          | The ALLELE_NUM a breakpoint record prints; see below                                        | `0`                                                       |
| `id`                                                                  | The record's name; see below                                                                | `"rs1234"`, `null`                                        |
| `variant_class`                                                       | Variant class                                                                               | `"Snv"`, `"Insertion"`, `"Deletion"`, `"StructuralDeletion"` |
| `transcript_consequences`                                             | The variant's whole list, so the variant document repeats every consequence                 |                                                           |
| `most_severe_consequence`                                             | The most severe consequence across all transcripts                                          | `"MissenseVariant"`                                       |
| `colocated_variants`                                                  | Co-located known variants                                                                   |                                                           |
| `hgvsg`                                                               | The allele's HGVS genomic notation; present under `--hgvsg` when the allele has one         |                                                           |
| `existing_variation`                                                  | The co-located known variants' IDs                                                          |                                                           |
| `plugin_data`                                                         | Variant-level built-in fields; present only when non-empty                                  | `{"CADD_PHRED": "23.1"}`                                  |
| `plugin_records`                                                      | The records a built-in groups under its name for the JSON output, one object per matched data line, a missing cell `null`; present only when non-empty | `{"OpenTargets": [{"OpenTargets_gwasGeneId": "ENSG00000154734"}]}` |
| `custom_annotations`                                                  | The `--custom` fields; absent at plugin time, since the custom sources run after the shared libraries |                                                           |
| `minimised`, `original_allele_string`, `original_start`, `original_end` | Whether the record was minimised under `--minimal`, and the alleles and coordinates before it; on an allele of a split multi-allelic record the record's span |                                                           |
| `rejoin_key_prefix`                                                   | On an allele after the first of a `--minimal` split, the trimmed allele string and `_` the output sorts its rows by; present only then | `"-/A_"`                                                 |
| `intergenic_shift`                                                    | The bases `--shift_3prime 1` moved an insertion or deletion that overlaps no transcript; present only then |                                             |
| `raw_input`, `input_record`                                           | The input line as [the JSON output's `input`](../output/json.md#the-record-object) stores it (trailing empty columns dropped, a structural record's QUAL, FILTER and INFO filled with `.`), and the ordinal of the record it came from | |
| `summary_terms`                                                       | The `--summary` or `--most_severe` row's terms; present on the one allele of a record that carries that row |                                            |
| `intergenic_pick`                                                     | Whether a `--flag_pick` variant marked the allele's intergenic row                          | `true`                                                    |
| `uploaded_allele_string`, `record_allele_string_multi`, `input_line_name` | The `Uploaded_variation` helpers; see below                                             |                                                           |

`allele_num` is the ALLELE_NUM a breakpoint record prints when it is not
`allele_index + 1`. Its values are the record's own breakend 0, its mates
from 1, a bare `<BND>`'s allele 2 behind the reference base. It is absent for
every other record.

`id` is the first `;`-joined ID of a VCF record as written (`.` and `0`
included), the whole ID column of an Ensembl-format line, or the identifier of
an identifier-list input (`--format id`). It is `null` when the column or its
first `;`-element is empty, so `;rs2` gives `null`.

`input_line_name` is the `CHROM_POS_REF/ALT` name built from the input line,
present for VCF and Ensembl-format input.

The structural-variant fields:

| Field                             | Meaning                                                                                                              | Example                                     |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| `is_structural`                   | Whether the allele is a structural variant                                                                           | `true`                                      |
| `sv_end`, `sv_type`, `sv_len`     | INFO END, SVTYPE and SVLEN                                                                                           | 9100000, `"DEL"`, `-31583`                |
| `ci_pos`, `ci_end`                | CIPOS and CIEND confidence intervals                                                                                 | `[-10, 10]`                                 |
| `mate_id`                         | Only on a record that names a mate, a bracket or `<BND>` record: MATEID, PARID or a stand-in; never on a single breakend |                                          |
| `mate_chr`, `mate_pos`            | The mate's chromosome and position                                                                                   | `"21"`, 9150000                           |
| `is_single_breakend`              | Whether the allele is a single breakend                                                                              | `true`                                      |
| `record_so_term`                  | The whole record's class term; see below                                                                             | `"Alu_insertion"`, `"copy_number_variation"` |

`record_so_term` is the whole record's class term as VEP's `class_SO_term` reads
it from every ALT joined by `/`. It is the Allele every structural row prints,
and is absent for a sequence variant.

### What `run` returns, and where it lands

The returned string must be a JSON object. Each entry becomes one field on that
row. A string value is stored verbatim and any other value as its JSON text
(`0.42`, `null`). A key already present, from a built-in or an earlier shared
library, is overwritten. On a consequence row the fields are the consequence's;
on an intergenic row they are the variant's, which that row prints.

The fields appear in the Extra column of the default output format and in the
per-consequence objects of `--json`. There, as for every other field, the key is
lowercased and a numeric-looking value is written as a number.

Every key the plugin declared in `vep_plugin_header_info` is also a column of
`--tab`, an entry of the CSQ Format of `--vcf` and a column of
`--output_format parquet`. Each gets a `## KEY : description` line in the
default and `--tab` headers and a `##KEY=description` line in the VCF header. A
built-in plugin's fields get the same lines (`OutputFactory/VEP_output.pm`,
`Tab.pm` and `VCF.pm`). An undeclared key the plugin writes reaches the Extra
column and `--json` only.

A key a plugin returns that names a column or field of the output (`Codons`,
`TSL`, `SYMBOL`) replaces that value on the row in every format, as Ensembl
VEP's `$line_hash->{$key} = $plugin_results->{$key}` (`OutputFactory.pm` line
2231) does after the core fields set theirs. A declared flag field the plugin
did not write keeps the row's own value under its flag and is empty without it.
A Parquet column keeps the type its name gives it
([Parquet output](../output/parquet.md)), so a value outside that type, a list
in the INT32 `TSL` column, is stored NULL. Keys carrying the plugin's name, as
the canonical plugins' do (`CADD_PHRED`, `SpliceAI_pred`), collide with nothing.

#### When a `run` call fails

A `run` call that returns NULL or text that is not JSON costs the row that
plugin's fields alone and writes the warning `Plugin '<name>' went wrong:
<reason>`. A JSON value that is not an object writes `Plugin '<name>' did not
return a hashref, output ignored!` and is dropped the same way. The row is kept,
the plugins after it still run on it, and the run goes on, as Ensembl VEP keeps
the line when a plugin returns something other than a hash or dies
(`OutputFactory.pm` lines 2226-2245). A warning text is written once per run,
so a plugin that fails the same way on every row draws one line.

Ensembl VEP reads a plugin's `undef` return as a request to drop the line from
the output (`OutputFactory.pm` lines 2238-2240, the `BaseVepFilterPlugin`
mechanism). The ABI has no filter return, so a shared library cannot remove a
row, and a defect in one costs its own fields rather than the run.

#### No `finish` call and no per-line hook

The ABI has no `finish` call. Ensembl VEP calls each plugin's `finish` method, when
it has one, after the last variant (`Runner.pm` lines 168-173). A shared library
that holds a file to flush does so in `vep_plugin_destroy`, which the host calls
on every loaded instance when the run ends. Ensembl VEP's per-line hook, the
`$line_hash` of the whole output line that `run` receives as its second argument
(`BaseVepPlugin.pm` lines 329-370), has no counterpart either. A shared library
reads the consequence and variant documents above, which carry the annotation the
line is rendered from, and writes fields, never the rendered line.

Two parts of that line are absent from the documents. The `--custom` fields are
absent because the custom sources run after the shared libraries. The fields
another shared library wrote on the same row are absent because the documents
are serialised once per row before the row's plugins run. Ensembl VEP's
`$line_hash` holds both (`OutputFactory.pm` lines 1003-1006 and 391-398).

### Batch query strategy

Tabix-based plugins use the `TabixAnnotator` engine (or the `BinaryAnnotator`
when the data file's `.vpd` store is accepted). The engine batches I/O across
the entire variant buffer (`--buffer_size`, default 5000):

1. Collect query regions from all variants in the batch
2. Group the regions by chromosome and sort each chromosome's regions by position
3. Merge overlapping regions, bridging gaps of up to 1000 bp
4. Issue one tabix query per merged region rather than one per variant
5. Index returned records by position
6. Match each variant against the indexed records

Ensembl VEP's `BaseVepTabixPlugin` runs one tabix query per variant whenever a
plugin sets no region expansion, as its CADD and gnomADc plugins do. The batched
strategy issues one per merged region. Its run-time contribution is not
separately measured; the published wall times had no plugin enabled.

### Plugin execution in the pipeline

Plugins execute after consequence assignment:

1. Input batch parsed
2. Consequences calculated (parallel if `--fork > 1`)
3. Co-located variants matched
4. **`registry.run_batch(&mut variants)`**: every built-in plugin's `prefetch()`
   runs in parallel across plugins, then each plugin's `annotate()` (or
   `run_batch()`, for a plugin without `prefetch`) runs sequentially on the full
   batch
5. Shared libraries run on every row whose kind they listed, the intergenic
   rows included for a plugin that listed `Intergenic`
6. The `--custom` sources match their records and add their fields
7. Output written

#### When a built-in plugin fails

A built-in plugin whose `prefetch()`, `annotate()` or `run_batch()` fails on a
batch of two or more variants is retried variant by variant. The plugins after
it still run. A variant whose own run succeeds keeps the plugin's fields; one
whose `prefetch()` or `annotate()` fails loses them, unless the plugin declares
`annotate_fails_per_row()`. One whose `run_batch()`, or such an `annotate()`,
fails keeps every row but the failed ones. Past the first variant failing on
its own, the retry therefore only counts.

`run_batch()` runs every row and reports the first failure after the loop.
Ensembl VEP's per-line `eval` of each plugin's `run` (`OutputFactory.pm` lines
2225-2245) likewise costs a dying line its fields alone.

The run ends with one `Plugin '<Name>' went wrong: <reason>; <footprint>; <N>
variants affected` line per failed plugin, from its first failing batch. After a
failed `prefetch()` or `annotate()` the footprint is `the plugin annotated
nothing on the variants it failed on, the first at <location>`. After a failed
`run_batch()` or per-row `annotate()` it is `the plugin annotated every row but
the ones it failed on, the first of them on the variant at <location>`.

#### Where a built-in plugin writes its fields

Variant-level plugins (CADD, gnomADc, GWAS, OpenTargets) write to
`variant.plugin_data`. A plugin whose variant matches several data lines
(OpenTargets) writes each field `,`-joined across them for the text formats. It
also puts the same lines, one object each, under its name in
`variant.plugin_records`. The JSON output writes that as an array in place of
the joined fields.

Transcript-level plugins write to `consequence.plugin_data`. REVEL, SpliceAI,
dbNSFP, dbscSNV and AlphaMissense do so in their batch `annotate()` after a
tabix `prefetch()`. LoFtool, pLI and LoFTEE do so through the per-consequence
`run()` method.

LoFTEE is the only consequence-_filtering_ plugin. It inspects the transcript
exon/intron model (`LofteeContext` on the consequence, populated when a loaded
built-in reads it). It does not look up a precomputed score by position.

An empty value in either map is Ensembl VEP's undefined value. The default
format writes `KEY=`, the tab format `-` and the JSON format no key.

A plugin whose source defines an empty string, such as a table cell left empty,
names the field in `consequence.plugin_empty_strings`. The consequence document
never carries that list. The tab format then writes an empty cell and the JSON
format `""`. The `GWAS_*` fields, `DownstreamProtein` and `AA` are empty strings
by name.

## Adding a new plugin

1. Create a new module in `crates/vep-builtins/src/plugins/my_plugin.rs`
2. Implement the `BuiltinPlugin` trait:

```rust
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::variant::InputVariant;
use vep_core::consequence::TranscriptConsequence;
use indexmap::IndexMap;

pub struct MyPlugin { /* fields */ }

impl BuiltinPlugin for MyPlugin {
    fn name(&self) -> &str { "MyPlugin" }
    fn header_info(&self) -> Vec<(String, String)> {
        vec![("MY_FIELD".into(), "Description".into())]
    }
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        // Parse params, open data files
        Ok(())
    }
    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Return key-value pairs for the Extra column
        Ok(IndexMap::new())
    }
    // Tabix plugins add prefetch() and annotate() for batch I/O; see below
}
```

3. Add the module to `crates/vep-builtins/src/plugins/mod.rs`
4. Register the name in `crates/vep-builtins/src/registry.rs` (`create_builtin()`)
5. Add tests

### For tabix-based plugins

Open the data file with `open_best_store_with` and
`TabixAnnotator::keeping_trailing_whitespace` (tabix, or the `.vpd` store when
accepted) and implement the two-phase protocol the built-in tabix plugins use.
`prefetch()` takes a shared reference to the batch, so the registry runs it for
every plugin in parallel, and `annotate()` applies the fetched records
sequentially. `run_batch()` is the single-phase fallback for a plugin that
returns `None` from `prefetch()`. A plugin that overrides any of the three
applies `annotates_variant()` to each variant it reads, as the default
`run_batch()` does.

The trait's `annotates_structural_variants()` is false unless the plugin's
Ensembl VEP module lists `StructuralVariationFeature` or
`BaseVariationFeature` in its `variant_feature_types`. A gated plugin writes
nothing on a structural variant's rows.

`open_best_store_with` applies a closure to the tabix reader and leaves a `.vpd`
store as built. `TabixAnnotator::keeping_trailing_whitespace` removes each
record line's line end alone, as Ensembl VEP's tabix iterator does, so an empty
or whitespace trailing cell stays a cell. Every store-backed built-in passes it,
gnomADc included, so a row reads alike through the table and a store. The
trimming reader (`TabixAnnotator::open_plugin_table`) suits a plugin reading
through tabix alone that never prints a row's last cell (LOEUF).

```rust
use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::tabix::{BatchQueryResult, TabixAnnotator, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

fn prefetch(
    &self,
    variants: &[InputVariant],
) -> Result<Option<PrefetchData>, PluginError> {
    let Some(store) = &self.store else { return Ok(None) };
    let regions: Vec<(&str, u64, u64)> = variants.iter()
        .filter(|v| self.annotates_variant(v))
        .map(|v| (v.chr.as_str(), v.start, v.end))
        .collect();
    Ok(Some(PrefetchData::new(store.query_batch(&regions)?)))
}

fn annotate(
    &self,
    variants: &mut [InputVariant],
    data: PrefetchData,
) -> Result<(), PluginError> {
    let batch: BatchQueryResult = data.downcast()?;
    for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
        let records = batch.lookup(&variant.chr, variant.start, variant.end);
        // Match alleles, then write to variant or consequence plugin_data
    }
    Ok(())
}
```

## Plugin concordance testing

`scripts/concordance/run_concordance.sh` runs Ensembl VEP (in Docker) and vep-rs on
the same inputs with the same plugins and compares each plugin's fields. Every
run needs four inputs, and it exits 1 before running anything when one is
missing:

- `--perl-cache-dir`, the directory that contains `homo_sapiens/<version>_<assembly>/`
- `--json-cache-dir`, a vep-rs JSON cache directory containing `info.json`
- `--fasta`, an indexed reference FASTA given to both engines; `--no-fasta` is
  the explicit opt-out and is recorded in `provenance.json`
- a benchmark directory, `--benchmark-dir` or `$VEP_BENCHMARK_DIR`

The plugin flags:

| Flag                          | Default                                | Description                                                                                                     |
| ----------------------------- | -------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| `--plugins <list>`            | (none)                                 | Comma-separated plugin names to test                                                                            |
| `--plugin-data-dir <path>`    | `~/.vep/plugin_data/<assembly>/<tier>` | Plugin data directory (lower-case assembly, from `--assembly`); a path that names the other assembly is refused |
| `--plugin-tier <quick\|full>` | `quick`                                | Subdirectory used by the default data directory                                                                 |

Two data tiers share the layout of `scripts/data/plugin_data_layout.sh`, the
open data files the setup script stages. `quick` holds chromosome-21 subsets of
them, enough for the chr21 benchmark VCFs that
`scripts/data/download_real_world_vcfs.sh` seeds. A genome-wide input
therefore receives plugin values only on its chr21 records. `full` holds the
genome-wide files.

CADD, SpliceAI, AVI, dbNSFP and dbscSNV read no file of the layout. The harness
takes your own file from `CADD_SNV_FILE`, `CADD_INDEL_FILE`,
`SPLICEAI_SNV_FILE`, `SPLICEAI_INDEL_FILE`, `AVI_FILE`, `DBNSFP_FILE` or
`DBSCSNV_FILE`, and leaves the plugin out where its variable is unset.

This command stages either tier, as
[Preparing data files](../plugins/data-files.md) describes. The script is
`scripts/data/setup_plugin_data.sh`, `<asm>` is the assembly, such as `GRCh38`,
and `<tier>` is `quick` or `full`:

```bash
setup_plugin_data.sh --assembly <asm> --tier <tier>
```

### Running the harness

An example harness run, in smoke mode on the `quick` tier:

```bash
scripts/concordance/run_concordance.sh \
  --assembly GRCh38 \
  --perl-cache-dir /path/to/vep_cache \
  --json-cache-dir ~/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --benchmark-dir .vep/vcf/grch38/all_variants \
  --plugins REVEL,gnomADc \
  --plugin-data-dir ~/.vep/plugin_data/grch38/quick \
  --mode smoke
```

With `--plugins` set, the harness:

- maps each plugin to a `--plugin Name,params` flag for both engines, with file
  paths from the layout table (Ensembl VEP's remapped under the container mount
  `/work/plugin_data`). A file of your own is read under
  `/work/user_plugin_data/<variable>`.
- mounts the data directory into the Ensembl VEP container, plus a
  `VEP_plugins/` checkout when one exists at `$VEP_PLUGINS_DIR` (default: a
  `VEP_plugins/` sibling of the repository). Otherwise the
  `ensemblorg/ensembl-vep` image's own `/plugins` copies serve.
- keys its Ensembl VEP output cache on the plugin list, Ensembl VEP's plugin
  argument strings, the data directory and a size-and-mtime fingerprint of its
  files. Cached outputs from another plugin configuration or parameter are
  therefore not reused.
- after the consequence comparison, runs `compare_plugin_outputs.py` per plugin
  and writes `reports/plugin_concordance.json` (per plugin: `match`, `mismatch`,
  unrounded `concordance_pct`, plus `tolerance`, `min_annotated`, `file_count`
  and per-file `results`).

A failing plugin makes the harness exit non-zero. Unless `--no-perf` is passed,
a baseline vep-rs run without plugins is also timed on the first prepared input.
`reports/perf/plugin_timing.json` records the plugin overhead and the wall time
of Ensembl VEP with plugins over that of vep-rs with plugins.

`compare_plugin_outputs.py` also runs standalone on two default-format outputs:

```bash
python3 scripts/concordance/compare_plugin_outputs.py \
    --perl-output perl_vep_with_cadd.txt \
    --rust-output rust_vep_with_cadd.txt \
    --plugins CADD \
    --fields CADD_PHRED,CADD_RAW
```

A release publishes a raw and an adjusted F1 for each supported plugin and
assembly, measured as
[Reproducing the plugin concordance](../plugins/reproducing-concordance.md)
describes. A plugin comparison here is a check on your own data and binary, not
a reproduction of those figures.

### Comparing LoFTEE

A docker-mode run with LoFTEE also needs `LOFTEE_PERL_PATH` set to a checkout of
[konradjk/loftee](https://github.com/konradjk/loftee) (branch `master` for
GRCh37, `grch38` for GRCh38). The harness mounts that checkout into the Ensembl
VEP container. `LOFTEE_DATA_DIR` overrides the default data location,
`<plugin-data-dir>/loftee/<assembly>/`.

The harness passes `run_splice_predictions:0` to Ensembl VEP's `LoF` on both
assemblies and `run_splice_predictions=0` to vep-rs on GRCh37. The reason is
that `master` runs the splice-prediction layer by default and vep-rs does not
implement it.

The harness hands vep-rs the plain `human_ancestor.fa` and Ensembl VEP's `LoF`
plugin the bgzipped `human_ancestor.fa.gz` (vep-rs reads either form).
`setup_plugin_data.sh` reports the directory READY only when the plain file, the
compressed copy and all three indexes are present.

The `END_TRUNC` cutoff differs in sign between the two upstream branches (`-58`
on `grch38`, `+180` on `master`). vep-rs derives it from `--assembly`, so pass
that flag rather than `gerp_end_trunc_cutoff=`. Compare `LoF` and `LoF_filter`
only (`--fields LoF,LoF_filter --tolerance 0`): `LoF_info` and `LoF_flags` match
on GRCh38 without bigWig zoom levels; with them `GERP_DIST` differs.

Three things to know when reading a LoFTEE comparison:

- On GRCh38, Ensembl VEP's LoFTEE plugin reads the GERP bigWig through
  [lossy zoom-level summaries](../plugins/loftee.md#5-gerp-totals-read-from-the-bigwigs-zoom-levels)
  and vep-rs sums exact per-base values, so a borderline `END_TRUNC` call can
  flip.
- On GRCh37, Ensembl VEP's `master` plugin takes whole intermediate exons from
  the `gerp_exons` table of its conservation database, which the harness does
  not pass. A record whose exon walk needs a whole-exon total therefore dies on
  the Ensembl VEP side and is counted as vep-rs-only.
- Ensembl VEP's `ANC_ALLELE` check shells out to `samtools faidx`, which the
  `ensemblorg/ensembl-vep` image does not ship. That filter is therefore never
  applied on the Ensembl VEP side in Docker.

## Where next

- [Plugins](../plugins.md): the user-facing list of plugins.
- [Preparing data files](../plugins/data-files.md): the file format every tabix
  plugin expects and the script that stages the files.
- [JSON output](../output/json.md): the keys of each object of the JSON output.
- [Parquet output](../output/parquet.md): the typed schema of the Parquet
  output.
