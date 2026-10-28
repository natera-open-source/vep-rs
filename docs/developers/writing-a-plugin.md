# Writing a plugin

This page describes the two plugin systems, the C ABI a shared-library plugin
exports and the JSON it receives, how a built-in plugin is added to the
`vep-builtins` crate, and how a plugin's fields are compared against Ensembl
VEP's. The user-facing list of plugins is on the [Plugins](../plugins.md) page.

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

When `--plugin Name,...` is specified, the runner checks the built-in registry
first. A name no built-in claims is resolved as a dylib: as a file path if one
exists at `Name`, otherwise as `libName.so` and then `libName.dylib` under
`--dir_plugins`. The name is the library file name, so these flags load
`/opt/vep/plugins/libecho_plugin.so`:

```bash
--dir_plugins /opt/vep/plugins --plugin echo_plugin
```

The string a plugin returns from `vep_plugin_name` appears only in log and error
messages.

`--dir_plugins` defaults to `<--dir>/Plugins`, as in Ensembl VEP. When
`--plugin` names a plugin and the directory does not exist, the run writes this
message to the warnings file; under `--safe` it is an error:

```text
Plugins directory '<dir>' not found.
```

A name that resolves to neither a built-in nor a library is written to the
warnings file as `Failed to load plugin Name: Name is not a built-in plugin and
no libName.so or libName.dylib is under <dir>`. The run continues without it,
which is Ensembl VEP's default. Under `--safe` it is an error, as Ensembl VEP's
is.

A file that exists but does not open as a shared library, lacks one of the eight
required symbols or reports an ABI version other than 1 is dropped with a warning. The warnings
file records it as `Failed to compile plugin Name: <reason>` and the run continues
without it, as Ensembl VEP's does when a plugin module fails to compile. Under
`--safe` it is an error, as Ensembl VEP's is.

A plugin whose `vep_plugin_create` returns NULL or whose `vep_plugin_init`
returns non-zero is written to the warnings file as `Failed to instantiate
plugin Name: <reason>` and the run continues without it, as Ensembl VEP's does
when a plugin's constructor dies. A built-in plugin whose parameters it cannot
use (a scores file that does not open) is dropped the same way. Under `--safe`
either is an error, as Ensembl VEP's is.

### The C ABI (version 1)

A plugin exports these eight symbols with these signatures
(`crates/vep-plugin/src/abi.rs`), and may export a ninth, `vep_plugin_feature_types`.

#### `vep_plugin_abi_version`

```rust
extern "C" fn() -> u32
```

Returns `1`; any other value is a load error.

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
returns 0 on success. Called only when at least one parameter was given.

#### `vep_plugin_run`

```rust
extern "C" fn(*const c_void, *const c_char, *const c_char) -> *mut c_char
```

Takes `(instance, consequence JSON, variant JSON)` and returns a NUL-terminated
JSON object that the host frees with `vep_plugin_free`. NULL signals an error.

#### `vep_plugin_free`

```rust
extern "C" fn(*mut c_char)
```

Frees a string returned by `vep_plugin_run`.

#### `vep_plugin_feature_types` (optional)

```rust
extern "C" fn(*const c_void) -> *const c_char
```

Returns the feature types the plugin annotates as a NUL-terminated,
comma-separated list of the words Ensembl VEP's `feature_types` takes:
`Transcript`, `RegulatoryFeature`, `MotifFeature`, `Feature` (the three
together) and `Intergenic`. The host reads it once after `vep_plugin_init`, so
the list may depend on the parameters, and the pointer must stay valid for the
instance's lifetime; any other word fails the load. A plugin without the symbol
annotates `Transcript` rows alone, Ensembl VEP's default. A plugin that names
`RegulatoryFeature` or `MotifFeature` switches regulatory annotation on for the
run, as Ensembl VEP does for such a plugin.

### What `run` receives

After the built-in plugins have run on a batch, the host calls `vep_plugin_run`
once per entry in each variant's `transcript_consequences` whose `feature_type`
the plugin declares, for every loaded dylib in `--plugin` order. A plugin
declares `Transcript` unless `vep_plugin_feature_types` says otherwise, as
Ensembl VEP runs a plugin on the feature classes its `feature_types` names. An intergenic row has
no entry, so `Intergenic` in the list gives a plugin no row; Ensembl VEP runs a
plugin that declares it on intergenic rows.

The two JSON documents are the serde serialization of `TranscriptConsequence`
(`crates/vep-core/src/consequence.rs`) and `InputVariant`
(`crates/vep-core/src/variant.rs`), and those two structs with their serde
attributes are the authoritative schema. No field is renamed, so a JSON key is
the Rust field name. A field marked `#[serde(skip)]` (`id_rank`,
`record_span_in_cds`, `raw_input_bytes`, `split_first_feature_rowless`,
`reference_row_shift_allele`) never appears, a field marked
`skip_serializing_if` is left out when it is `None`, empty or zero, and every
other absent optional value is `null`.

**The consequence document.** Every field a plugin will find:

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
| `exon`, `intron`                                   | Exon or intron number over the total                                                         | `"3/5"`, `"2/4"`                                   |
| `hgvsc`, `hgvsp`                                   | HGVS coding and protein notations, as the output prints them before escaping; `--hgvsp_use_prediction` and `--remove_hgvsp_version` shape `hgvsp` here too |  |
| `hgvs_offset`                                      | The signed number of bases the HGVS notations were shifted along the transcript; `null` when nothing shifted |                                   |
| `shift_length`                                     | The bases `--shift_3prime 1` moved the allele along the transcript; `0` otherwise            | `0`                                                |
| `shifted_allele`                                   | Under `--shift_3prime 1`, the inserted bases rotated along the shift, the `Allele` a coding row prints; present only on such a row | `"GAA"`                                            |
| `used_ref`                                         | Under `--use_transcript_ref`, the reference allele the row was analysed against, on the variant's strand: the transcript's inside its exons, the input's elsewhere; present only on a sequence variant's row under the flag | `"A"`                             |
| `sift`, `polyphen`                                 | Prediction and score                                                                         | `"deleterious(0.01)"`                              |
| `domains`                                          | Array of `[source, id]` pairs                                                                | `[["Pfam", "PF00080"]]`                            |
| `feature_type`                                     | The feature kind                                                                             | `"Transcript"`, `"RegulatoryFeature"`, `"MotifFeature"` |
| `flags`                                            | Transcript flags                                                                             | `["cds_start_NF"]`                                 |
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
| `feature_start`, `feature_end`                     | Genomic span of the feature, 1-based inclusive; both `0` when unknown                        |                                                    |
| `plugin_data`                                      | The fields the built-in plugins have already written to this consequence; present only when non-empty | `{"REVEL_score": "0.931"}`                |
| `loftee_ctx`                                       | The exon/intron model; present only when LoFTEE is enabled                                   |                                                    |
| `motif`                                            | A motif row's `name` (the binding matrix's stable id), `pos` (the variant's position in the matrix, absent outside it) and `factors`; present on a sequence variant's `MotifFeature` row alone | `{"name": "ENSPFM0290", "pos": 7, "factors": ["HNF4A"]}` |

`cdna_total`, `cds_total` and `protein_total` are the lengths `--total_length`
appends to the three positions. Each is present only when its rule holds,
whatever the flag: the cDNA total wherever the variant overlaps an exon, the
coding and peptide totals where the variant lies in the coding region and the
cache carries the sequence.

**The variant document.** The sequence-variant fields:

| Field                                                                 | Meaning                                                                                     | Example                                                   |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `chr`                                                                 | Chromosome, normalized                                                                      | `"21"`                                                    |
| `original_chr`                                                        | The name the record prints: as input, or `"MT"` for a bare `M` the cache knows only as `MT` | `"chr21"`, `"MT"`                                         |
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
| `minimised`, `original_allele_string`, `original_start`, `original_end` | Whether the record was minimised under `--minimal`, and the alleles and coordinates before it; on an allele of a split multi-allelic record the record's span |                                                           |
| `rejoin_key_prefix`                                                   | On an allele after the first of a `--minimal` split, the trimmed allele string and `_` the output sorts its rows by; present only then | `"-/A_"`                                                 |
| `intergenic_shift`                                                    | The bases `--shift_3prime 1` moved an insertion or deletion that overlaps no transcript; present only then |                                             |
| `raw_input`, `input_record`                                           | The input line as [the JSON output's `input`](../output/json.md#the-record-object) stores it (trailing empty columns dropped, a structural record's QUAL, FILTER and INFO filled with `.`), and the ordinal of the record it came from | |
| `summary_terms`                                                       | The `--summary` or `--most_severe` row's terms; present on the one allele of a record that carries that row |                                            |
| `intergenic_pick`                                                     | Whether a `--flag_pick` variant marked the allele's intergenic row                          | `true`                                                    |
| `uploaded_allele_string`, `record_allele_string_multi`, `input_line_name` | The `Uploaded_variation` helpers; see below                                             |                                                           |

`allele_num` is the ALLELE_NUM a breakpoint record prints when it is not
`allele_index + 1`: the record's own breakend 0, its mates from 1, a bare
`<BND>`'s allele 2 behind the reference base. It is absent for every other
record.

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
consequence: a string value is stored verbatim, any other value as its JSON text
(`0.42`, `null`), and a key already present, from a built-in or an earlier
dylib, is overwritten. The fields appear in the Extra column of the default
output format and in the per-consequence objects of `--json`, where, as for
every other field, the key is lowercased and a numeric-looking value is written
as a number.

The ABI has no header-declaration symbol, so a dylib field gets no `##`
description line in the header and no column in `--tab`, `--vcf` (CSQ) or
`--output_format parquet` output; those writers take their plugin columns from
the built-in plugins' header info. A `run` call that returns NULL, invalid JSON,
or a JSON value that is not an object aborts the run.

### Worked example

`crates/vep-plugin/examples/echo_plugin.rs` is a complete plugin: it echoes the
consequence's `transcript_id` as `ECHO_FEATURE` and its first parameter, when
one was given, as `ECHO_PARAM`. It reads the one field it needs with a string
scan, so it adds no dependency. Build it as the shared library the loader opens,
then pass it by path or by name:

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

Every consequence line then carries this in its Extra column:

```text
ECHO_FEATURE=<transcript_id>;ECHO_PARAM=hello
```

`cargo test -p vep-plugin` builds the example and runs
`tests/echo_roundtrip.rs`, which loads it through the same loader, by path and
by name.

### Batch query strategy

Tabix-based plugins use the `TabixAnnotator` engine (or the `BinaryAnnotator`
when a `.vpd` store sits beside the data file), which batches I/O across the
entire variant buffer (`--buffer_size`, default 5000):

1. Collect query regions from all variants in the batch
2. Group the regions by chromosome and sort each chromosome's regions by position
3. Merge overlapping regions, bridging gaps of up to 1000 bp
4. Issue one tabix query per merged region rather than one per variant
5. Index returned records by position
6. Match each variant against the indexed records

Ensembl VEP's `BaseVepTabixPlugin` runs one tabix query per variant whenever a
plugin sets no region expansion, as its CADD and gnomADc plugins do. The batched
strategy replaces that with one query per merged region. Its contribution to run
time is not separately measured: the published wall-time measurements were taken
with no plugin enabled.

### Plugin execution in the pipeline

Plugins execute after consequence assignment:

1. Input batch parsed
2. Consequences calculated (parallel if `--fork > 1`)
3. Co-located variants matched
4. **`registry.run_batch(&mut variants)`**: every built-in plugin's `prefetch()`
   runs in parallel across plugins, then each plugin's `annotate()` (or
   `run_batch()`, for a plugin without `prefetch`) runs sequentially on the full
   batch
5. Dynamic (dylib) plugins run on every row whose feature type they declare
   (transcript rows alone without `vep_plugin_feature_types`)
6. Output written

Variant-level plugins (CADD, gnomADc, GWAS) write to `variant.plugin_data`.
Transcript-level plugins write to `consequence.plugin_data`: REVEL, SpliceAI,
dbNSFP, dbscSNV and AlphaMissense do so in their batch `annotate()` after a
tabix `prefetch()`, and LoFtool, pLI and LoFTEE
through the per-consequence `run()` method. LoFTEE is the only
consequence-_filtering_ plugin: it inspects the transcript exon/intron model
(carried on the consequence via `LofteeContext`, populated only when LoFTEE is
active) rather than looking up a precomputed score by position.

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

Open the data file with `open_best_store` (tabix, or the `.vpd` store when
present) and implement the two-phase protocol the built-in tabix plugins use.
`prefetch()` takes a shared reference to the batch, so the registry runs it for
every plugin in parallel, and `annotate()` applies the fetched records
sequentially. `run_batch()` is the single-phase fallback for a plugin that
returns `None` from `prefetch()`.

```rust
use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

fn prefetch(
    &self,
    variants: &[InputVariant],
) -> Result<Option<PrefetchData>, PluginError> {
    let Some(store) = &self.store else { return Ok(None) };
    let regions: Vec<(&str, u64, u64)> = variants.iter()
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
    for variant in variants.iter_mut() {
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

Two data tiers share the layout of `scripts/data/plugin_data_layout.sh`. `quick`
holds chromosome-21 slices of the full files, enough for the chr21 benchmark
VCFs that `scripts/data/download_real_world_vcfs.sh` seeds, so a genome-wide
input receives plugin values only on its chr21 records. Where the upstream host
serves a working remote tabix index, `tabix <url> 21 | bgzip` produces the slice
without a full download. `full` holds the genome-wide files. This command stages
either, as [Preparing data files](../plugins/data-files.md) describes:

```bash
setup_plugin_data.sh --assembly <asm> --tier <tier>
```

An example harness run, in smoke mode on the `quick` tier:

```bash
scripts/concordance/run_concordance.sh \
  --assembly GRCh38 \
  --perl-cache-dir /path/to/vep_cache \
  --json-cache-dir ~/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --benchmark-dir .vep/vcf/grch38/all_variants \
  --plugins CADD,REVEL,gnomADc \
  --plugin-data-dir ~/.vep/plugin_data/grch38/quick \
  --mode smoke
```

With `--plugins` set, the harness:

- maps each plugin to a `--plugin Name,params` flag for both engines, with file
  paths from the layout table (Perl's remapped under the container mount
  `/work/plugin_data`);
- mounts the data directory into the Perl container, plus a `VEP_plugins/`
  checkout when one exists at `$VEP_PLUGINS_DIR` (default: a `VEP_plugins/`
  sibling of the repository; otherwise the `ensemblorg/ensembl-vep` image's own
  `/plugins` copies serve);
- keys its Perl output cache on the plugin list, the Perl plugin argument
  strings, the data directory and a size-and-mtime fingerprint of its files, so
  cached outputs from another plugin configuration or parameter are not reused;
- after the consequence comparison, runs `compare_plugin_outputs.py` per plugin
  and writes `reports/plugin_concordance.json` (per plugin: `match`, `mismatch`,
  unrounded `concordance_pct`, plus `tolerance`, `min_annotated`, `file_count`
  and per-file `results`).

A failing plugin makes the harness exit non-zero. Unless `--no-perf` is passed,
a baseline vep-rs run without plugins is also timed on the first prepared input,
and `reports/perf/plugin_timing.json` records the plugin overhead and the
Perl-with-plugins over vep-rs-with-plugins wall time.

`compare_plugin_outputs.py` also runs standalone on two default-format outputs:

```bash
python3 scripts/concordance/compare_plugin_outputs.py \
    --perl-output perl_vep_with_cadd.txt \
    --rust-output rust_vep_with_cadd.txt \
    --plugins CADD \
    --fields CADD_PHRED,CADD_RAW
```

**LoFTEE.** A docker-mode run with LoFTEE also needs `LOFTEE_PERL_PATH` set to a
checkout of [konradjk/loftee](https://github.com/konradjk/loftee) (branch
`master` for GRCh37, `grch38` for GRCh38), which the harness mounts into the
Perl container; `LOFTEE_DATA_DIR` overrides the default data location,
`<plugin-data-dir>/loftee/<assembly>/`. The harness passes
`run_splice_predictions:0` to Perl's `LoF` on both assemblies, because `master`
(GRCh37) runs the splice-prediction layer by default and vep-rs does not
implement it.

The harness hands vep-rs the plain `human_ancestor.fa` and Perl's `LoF` plugin
the bgzipped `human_ancestor.fa.gz` (vep-rs reads either form), and
`setup_plugin_data.sh` reports the directory READY only when the plain file, the
compressed copy and all three indexes are present. The `END_TRUNC` cutoff
differs in sign between the two upstream branches (`-58` on `grch38`, `+180` on
`master`); vep-rs derives it from `--assembly`, so pass that flag rather than
`gerp_end_trunc_cutoff=`. Compare `LoF` and `LoF_filter` only
(`--fields LoF,LoF_filter --tolerance 0`): `LoF_info` and `LoF_flags` are
informational and differ in formatting between the two implementations.

Three things to know when reading a LoFTEE comparison:

- On GRCh38, Perl's plugin reads the GERP bigWig through lossy zoom-level
  summaries while vep-rs sums exact per-base values, so a borderline variant can
  land on either side of the `END_TRUNC` cutoff.
- On GRCh37, Perl's `master` plugin takes whole intermediate exons from the
  `gerp_exons` table of its conservation database, which the harness does not
  pass, so a record whose exon walk needs a whole-exon total dies on the Perl
  side and is counted as vep-rs-only.
- Perl's `ANC_ALLELE` check shells out to `samtools faidx`, which the
  `ensemblorg/ensembl-vep` image does not ship, so that filter is never applied
  on the Perl side in Docker.

No per-plugin agreement figure is published: the released concordance results
are core-annotation F1 with no plugin enabled, and a plugin comparison here is a
check on your own data and binary, not a reproduction of a published number.
