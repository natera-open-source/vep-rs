# Embedding vep-rs as a library

The contract of the `vep-cli` crate's library surface says what a program
that links vep-rs may rely on, how it adds the crates and how the surface is
versioned.

The [Writing a plugin](writing-a-plugin.md) page covers the plugin interfaces.
This page covers a program that holds an annotator and annotates records of
its own, as a service does. Ensembl VEP serves that shape through
`Runner::run_rest`, which annotates one input and returns the JSON output's
objects. The vep-rs surface returns the same objects from a long-lived
annotator.

## Adding the crates

The release workflow publishes the workspace's crates to crates.io when a
maintainer publishes a release page, each under the release's version, so a
release's version appears there minutes after its page. A program adds `vep-cli` from crates.io
at a published version, or as a git dependency pinned to a release's tag, the
form for a release crates.io does not carry:

```toml
[dependencies]
vep-cli = "=X.Y.Z"
```

```toml
[dependencies.vep-cli]
git = "https://github.com/natera-open-source/vep-rs"
tag = "vX.Y.Z"
```

`X.Y.Z` is a released version listed on the [Changelog](../changelog.md) page.
Every release carries a tag `vX.Y.Z`, the binary's `--version` prints the
same number, and neither a tag nor a published version ever moves. A pinned
program therefore builds the same sources until its pin changes.

`vep-cli` re-exports what the embedding path needs, so most programs add it
alone. `vep-core`, `vep-io`, `vep-effects`, `vep-builtins`, `vep-plugin` and
`vep-fasta` are added the same way, at the same version, when a program names
one of their types directly.

## The surface

The surface is the list of `pub use` re-exports at the top of
`crates/vep-cli/src/lib.rs`. A program that uses one of them is relying on this
page; everything else in `vep_cli` is a module of the binary, reachable at its
module path but not bound by this page.

- From `vep_cli::annotator`:
  - `Annotator`: one loaded cache and its settings, built once and shared.
    `Clone` is cheap, since every part of it sits behind an `Arc`, so clones
    serve many threads without a wrapper of their own.
  - `AnnotationResources`, `PredictionConfig`, `build_thread_pool`: the parts
    an annotator is built from by hand, for a program that loads them itself.
  - `RawRow`, `validate_raw_row`, `row_to_input_variant`: the simple-allele
    input, one chromosome, 1-based position, reference and alternate allele
    per row.
  - `AnnotationError`, `AnnotationStats`: the per-row error and the counts an
    annotation returns.
- From `vep_cli::config`: `Config`, the binary's resolved command line.
- From `vep_cli::pick` and `vep_cli::row_filter`: `apply_to_batch`,
  `apply_to_record`, `FilterConfig`, `PickCriterion` and `RowFilters`, the
  `--pick` family and the `--coding_only` and `--no_intergenic` filters. The
  binary applies them to an annotated record before it writes it.
- From `vep_cli::transcript_index`: `LazyTranscriptIndexes`,
  `MaterializedTranscriptIndexes`, `TranscriptSet`, the transcript indexes an
  annotator reads.
- From `vep-core`: `InputVariant`, the record an annotation fills, and
  `Consequence`, `Impact`, `TranscriptConsequence`, the types of what it
  holds.
- From `vep-io`: `JsonOutputFormatter`, the writer of the `--json` output.

## The embedding path

The path is the binary's own, without its reader and writer:

1. `Config::from_argv([...])` takes the binary's command line, the binary name
   first and without the input and output files. It resolves that line as the
   binary does: the option sets, the refusals and the settings of every flag.
   It reads the flags as typed. The `--config` file, `vep.ini` and the `VEP_*`
   variables are read by the step the binary runs before it,
   `vep_cli::argv::prepass`, a module path.

   ```rust
   let config = Config::from_argv([
       "vep", "--offline", "--json_cache", "/caches/116_GRCh38",
       "--species", "homo_sapiens", "--assembly", "GRCh38",
   ])?;
   ```

2. `Annotator::from_config(&config)` resolves the cache the configuration
   names; a configuration with no cache is an error, not an empty index. It
   opens the FASTA when one is named and settles the HGVS and prediction flags
   against the cache, as the binary does. It loads the plugins, then every
   chromosome's transcripts up front, so each request is a lookup. The
   annotator holds a thread pool sized by `--fork`;
   `from_config_with_pool` takes one the program built.

3. The input: `try_annotate_rows(&rows)` takes `RawRow` values, one allele
   per row, and returns a result per row at the same index. A result is the
   annotated `InputVariant`, or `AnnotationError::InvalidInput` for a row the
   simple shape does not carry:

   - a position of 0
   - an empty chromosome or allele
   - `ALT` written `.`
   - a multi-allelic, symbolic or breakend allele

   For a whole VCF record a program builds the record's `InputVariant` values
   itself and calls `annotate_batch(&mut record)`, or `try_annotate_batch`
   for a result per allele. The reader the binary parses a data line with is
   `vep_cli::vcf_parser::parse_vcf_line`, a module path. Beside it,
   `parse_vcf_line_with` takes the run's own `--max_sv_size`, `--minimal` and
   FASTA where the binary's differ from the defaults.

4. The output: `JsonOutputFormatter::new(options, plugin_fields, assembly)`
   writes a record as `vep --json` writes it.

   - `options` selects the optional field groups. Its `Default` is the field
     set of a run with no flag beyond the cache, the species and the assembly.
   - `plugin_fields` is `annotator.plugin_header_fields()`.
   - `assembly` is the name written as the record's `assembly_name`. The
     binary writes the cache's own assembly, else the one in the name of its
     directory under `--dir_cache`, else the `--assembly` value. `None` leaves
     the key out.

   `format_record_value(&alleles)` returns the record's object as a
   `serde_json::Value`, over every allele parsed from one line.
   `serde_json::to_string` of it, through the program's own `serde_json`
   dependency, is the line the binary writes.

### Optional field groups and escaping

A program that switches on an optional field group names the group twice. The
flag goes on the command line (`--symbol`, `--hgvs`, `--sift b`), so the
annotator computes the field. The matching option goes in the formatter's
options, so the field is written. The options type is
`vep_io::output::fields::FieldOptions`, named through a `vep-io` dependency at
the same version.

A BAM-edited cache switches `--use_transcript_ref` on in the annotator's copy
of the configuration, unless the run gives `--use_given_ref`. The program
therefore sets the options' `use_transcript_ref` and `bam_edited` from
`annotator.use_transcript_ref()` and `annotator.bam_edited()`, not from its
own `Config`.

`vep --json` implies `--no_escape`, as Ensembl VEP's `json` option set switches
on `rest` and `rest` switches on `no_escape` (`Config.pm` lines 485-498). The
binary therefore writes a synonymous HGVSp change as `=`.
`FieldOptions::no_escape` is `false` by default. A program writing the binary's
bytes sets it to `true`, or its `--hgvs` output carries `%3D` where the binary
writes `=`.

## What the library leaves to the program

The library has no skip step. The binary validates each record before
annotation and reports the records it drops in its warnings file; a program
applies the checks it wants. One check matters for the output's shape: the
binary writes no JSON object for a record whose span exceeds `--max_sv_size`
(10,000,000 bases unless the flag says otherwise), so a program that wants the
binary's output leaves such a record out.

The binary skips a record on a chromosome the cache does not hold, with a
warning, unless a `--custom` or `--regulatory_gff` file lists it. Such a
record comes back from the library with no consequence and no
`most_severe_consequence`.

The binary's validation also upper-cases a lower-case allele before
annotation, keeping the record's name as written. The library annotates a
record as the program hands it over, so a lower-case allele writes different
bytes.

A program that builds its own transcript index, where `Annotator::from_config`
would load one, hands it the reference sequence with
`LazyTranscriptIndexes::with_reference`. The binary pairs a RefSeq model's stop
codon that an `_rna_edit` completes with the bases the genome holds of it. An
index without the reference sequence keeps the cache's pairs. A record at
such a stop base then reads otherwise than the binary reads it.

The binary's notices (`INFO:` lines) are `tracing` events at the `INFO` level
and its warnings `tracing` events at the `WARN` level. Both go to the
subscriber the program installs. The warnings the binary writes to its
warnings file without a `tracing` event go to standard error instead:

- the plugin loader's warnings
- the warning of `--check_existing` on a cache that holds no variation files
- `Unable to assign consequence type`, the line `--summary` and
  `--most_severe` write for a record with no row

The statistics file, the warnings file and the output headers are the binary's
alone.

The warning the SpliceVault, AlphaMissense, REVEL, pLI and dbNSFP plugins
print under `--refseq` or `--merged`, that transcript matching may fail, goes
to standard error from the plugin itself, in a program as in the binary.

## What is proven

`crates/vep-cli/tests/library_embedding.rs` runs this path over every record of
the `tests/golden/116/GRCh37` dataset, the records past `--max_sv_size` left
out, and compares the lines with the file `vep --json` writes from the same
command line: the two are equal byte for byte. The reader is the binary's own,
so the comparison proves the annotator and the formatter.

## Versioning

While the version is 0.y.z, a release that changes any crate's public API,
this surface included, raises the minor number; a patch release keeps every
public item's name and signature. A pin moves between patch releases of one
minor line without a code change; moving it across minor lines starts with
the [Changelog](../changelog.md) page.
