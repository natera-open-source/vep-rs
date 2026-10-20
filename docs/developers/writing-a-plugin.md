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

A name that resolves to neither is reported with the message below, and the run
continues without it, which is Ensembl VEP's default without `--safe`:

```text
plugin 'Name' not found as a built-in or as a dylib under --dir_plugins; skipped
```

A file that exists but does not open as a shared library, lacks one of the eight
symbols, reports an ABI version other than 1, or returns non-zero from
`vep_plugin_init` aborts the run before annotation starts.

### The C ABI (version 1)

A plugin exports these eight symbols with these signatures
(`crates/vep-plugin/src/abi.rs`).

#### `vep_plugin_abi_version`

```rust
extern "C" fn() -> u32
```

Returns `1`; any other value is a load error.

#### `vep_plugin_create`

```rust
extern "C" fn() -> *mut c_void
```

Returns an opaque instance pointer; NULL is a load error.

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

### What `run` receives

After the built-in plugins have run on a batch, the host calls `vep_plugin_run`
once per entry in each variant's `transcript_consequences` (transcript,
regulatory-feature and motif-feature consequences alike), for every loaded dylib
in `--plugin` order. A variant with no such entries, such as an intergenic one,
is never passed to a plugin.

The two JSON documents are the serde serialization of `TranscriptConsequence`
(`crates/vep-core/src/consequence.rs`) and `InputVariant`
(`crates/vep-core/src/variant.rs`). The Rust field names are the JSON keys, and
those two structs are the authoritative schema. Absent optional values are
`null`.

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
| `amino_acids`                                      | Amino acid change                                                                            | `"A/V"`                                            |
| `codons`                                           | Codon change                                                                                 | `"gCc/gTc"`                                        |
| `protein_id`                                       | Ensembl protein ID                                                                           |                                                    |
| `distance`                                         | Distance to the transcript, for upstream and downstream consequences                         | `4123`                                             |
| `strand`                                           | Transcript strand                                                                            | `1`, `-1`                                          |
| `exon`, `intron`                                   | Exon or intron number over the total                                                         | `"3/5"`, `"2/4"`                                   |
| `hgvsc`, `hgvsp`                                   | HGVS coding and protein notations                                                            |                                                    |
| `sift`, `polyphen`                                 | Prediction and score                                                                         | `"deleterious(0.01)"`                              |
| `domains`                                          | Array of `[source, id]` pairs                                                                | `[["Pfam", "PF00080"]]`                            |
| `feature_type`                                     | The feature kind                                                                             | `"Transcript"`, `"RegulatoryFeature"`, `"MotifFeature"` |
| `flags`                                            | Transcript flags                                                                             | `["cds_start_NF"]`                                 |
| `tsl`                                              | Transcript Support Level                                                                     | `1`                                                |
| `mane_select`, `mane_plus_clinical`                | MANE transcript IDs, when the transcript carries the tag                                     |                                                    |
| `appris`                                           | APPRIS annotation                                                                            | `"principal1"`                                     |
| `ccds`, `swissprot`, `trembl`, `refseq`            | Cross-references                                                                             |                                                    |
| `feature_start`, `feature_end`                     | Genomic span of the feature, 1-based inclusive; both `0` when unknown                        |                                                    |
| `plugin_data`                                      | The fields the built-in plugins have already written to this consequence; present only when non-empty | `{"REVEL_score": "0.931"}`                |
| `loftee_ctx`                                       | The exon/intron model; present only when LoFTEE is enabled                                   |                                                    |

**The variant document.** The sequence-variant fields:

| Field                                                                 | Meaning                                                                                     | Example                                                   |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `chr`                                                                 | Chromosome, normalized                                                                      | `"21"`                                                    |
| `original_chr`                                                        | Chromosome as input                                                                         | `"chr21"`                                                 |
| `start`, `end`                                                        | 1-based inclusive; `end` is `start - 1` for an insertion                                    | 9068417                                                 |
| `strand`                                                              | Strand                                                                                      | `"Forward"`, `"Reverse"`                                  |
| `allele_string`                                                       | Alleles in VEP form                                                                         | `"A/G"`, `"AC/-"`, `"-/TT"`                               |
| `ref_allele`, `alt_alleles`                                           | Serialized as arrays of byte values, so read `allele_string` for the sequences              | `[65]`                                                    |
| `allele_index`                                                        | Which ALT of the record this allele is                                                      | `0`                                                       |
| `allele_num`                                                          | The ALLELE_NUM a breakpoint record prints; see below                                        | `0`                                                       |
| `id`                                                                  | The record's name; see below                                                                | `"rs1234"`, `null`                                        |
| `variant_class`                                                       | Variant class                                                                               | `"Snv"`, `"Insertion"`, `"Deletion"`, `"StructuralDeletion"` |
| `hgvsg`                                                               | Genomic HGVS notation under `--hgvsg`; absent otherwise                                     | `"21:g.25_26del"`                                         |
| `transcript_consequences`                                             | The variant's whole list, so the variant document repeats every consequence                 |                                                           |
| `most_severe_consequence`                                             | The most severe consequence across all transcripts                                          | `"MissenseVariant"`                                       |
| `colocated_variants`, `existing_variation`                            | Co-located known variants and their IDs                                                     |                                                           |
| `plugin_data`                                                         | Variant-level built-in fields; present only when non-empty                                  | `{"CADD_PHRED": "23.1"}`                                  |
| `minimised`, `original_allele_string`, `original_start`, `original_end` | Whether a multi-allelic record was minimised, and the alleles and coordinates before it     |                                                           |
| `raw_input`, `input_record`                                           | The input line, and the ordinal of the record it came from                                  |                                                           |
| `uploaded_allele_string`, `record_allele_string_multi`, `input_line_name` | The `Uploaded_variation` helpers; see below                                             |                                                           |

`allele_num` is the ALLELE_NUM a breakpoint record prints when it is not
`allele_index + 1`: the record's own breakend 0, its mates from 1, a bare
`<BND>`'s allele 2 behind the reference base. It is absent for every other
record.

`id` is the first `;`-joined ID of a VCF record as written, `.` and `0`
included, or the whole ID column of an Ensembl-format line. It is `null` when
the column or its first `;`-element is empty, so `;rs2` gives `null`.

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

Perl VEP's `BaseVepTabixPlugin` runs one tabix query per variant whenever a
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
5. Dynamic (dylib) plugins run on every transcript consequence
6. Output written

Variant-level plugins (CADD, gnomADc, AlphaMissense, GWAS) write to
`variant.plugin_data`. Transcript-level plugins write to
`consequence.plugin_data`: REVEL, SpliceAI, dbNSFP and dbscSNV do so in their
batch `annotate()` after a tabix `prefetch()`, and LoFtool, pLI and LoFTEE
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

`scripts/concordance/run_concordance.sh` runs Perl VEP (in Docker) and vep-rs on
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

vep-rs reads the plain `human_ancestor.fa`, so decompress the upstream bgzipped
file and re-index it (`bgzip -d`, then `samtools faidx`). The harness hands
Perl's `LoF` plugin the bgzipped `human_ancestor.fa.gz`, and
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
