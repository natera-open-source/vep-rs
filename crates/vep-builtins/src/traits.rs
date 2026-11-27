// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin trait definition for built-in VEP plugins.
//!
//! Built-in plugins operate on native Rust types without JSON serialization
//! or C FFI overhead. The plugin execution model has three layers:
//!
//! - `prefetch` + `annotate`:
//!   Two-phase execution that separates I/O (parallelizable) from annotation
//!   (sequential). Tabix plugins implement these to enable parallel prefetch
//!   across independent data files.
//! - `run_batch`: Single-phase fallback for plugins
//!   that don't benefit from parallel I/O (e.g., in-memory lookups).
//! - `run`: Per-transcript-consequence annotation for
//!   simple per-transcript plugins.

use std::sync::Arc;

use indexmap::IndexMap;
use rayon::prelude::*;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;

use crate::pool_pass::{about_record, map_in_order};

/// Errors from builtin plugin execution.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PluginError {
    #[error("plugin initialization failed: {0}")]
    Init(String),
    #[error("plugin execution failed: {0}")]
    Run(String),
    #[error("tabix I/O error: {0}")]
    Io(#[from] std::io::Error),
}

impl PluginError {
    /// The error's own text, without the words its variant prints before it.
    pub fn detail(&self) -> String {
        match self {
            Self::Init(text) | Self::Run(text) => text.clone(),
            Self::Io(error) => error.to_string(),
        }
    }

    /// The same error with `note` after its text, `; ` between, in the same
    /// variant.
    pub(crate) fn with_note(self, note: &str) -> Self {
        match self {
            Self::Init(text) => Self::Init(format!("{text}; {note}")),
            Self::Run(text) => Self::Run(format!("{text}; {note}")),
            Self::Io(error) => Self::Io(std::io::Error::new(
                error.kind(),
                format!("{error}; {note}"),
            )),
        }
    }
}

/// The output format of the run, as Ensembl VEP's plugins read it from
/// `$self->{config}->{output_format}`: `vep` by default (Config.pm's default),
/// `tab`, `vcf` or `json` under the format flags. SpliceAI prints one field or
/// nine by it (SpliceAI.pm:173-175, 208, 275-278) and AlphaMissense nests its
/// fields under a group key in JSON (AlphaMissense.pm:281-282). Parquet is a
/// vep-rs format with no Ensembl VEP counterpart: one column per field, so a
/// plugin shapes its fields for it as for the tab format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
#[doc(hidden)]
pub enum OutputFormat {
    #[default]
    Vep,
    Tab,
    Vcf,
    Json,
    Parquet,
}

impl OutputFormat {
    /// The format named by its canonical writer name (`vep`, `tab`, `vcf`,
    /// `json`, `parquet`); any other name is the default format.
    pub fn from_name(name: &str) -> Self {
        match name {
            "tab" => Self::Tab,
            "vcf" => Self::Vcf,
            "json" => Self::Json,
            "parquet" => Self::Parquet,
            _ => Self::Vep,
        }
    }
}

/// The run's settings Ensembl VEP's plugins read from `$self->{config}` while
/// they construct: the output format, the `--refseq` / `--merged` transcript
/// set (REVEL.pm:144-146, AlphaMissense.pm:188-190, pLI.pm:187-189,
/// SpliceVault.pm:129-131 and dbNSFP.pm:617-619 warn under either that
/// transcript matching may fail) and
/// the plugins directories, in search order (the last `--dir_plugins` match
/// first, as `unshift @INC` leaves them, Runner.pm:921), where LoFtool.pm:74-78
/// and pLI.pm:115-119 find their bundled score files beside the module that
/// loaded (`$INC{'<Plugin>.pm'}`), and whether the run computes the HGVS
/// notations (`--hgvs` with a reference FASTA), which a plugin reading the
/// HGVSc coordinates needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
#[doc(hidden)]
pub struct RunContext<'a> {
    pub output_format: OutputFormat,
    pub refseq: bool,
    pub merged: bool,
    pub hgvs: bool,
    pub plugin_dirs: &'a [String],
}

/// Prints `line` on standard error alone, the channel a bare `warn` in an
/// Ensembl VEP plugin reaches: printed under `--quiet` too and never in the
/// warnings file, which carries the run's own warnings and is compared with
/// Ensembl VEP's. A line the stream refuses is dropped once the process has
/// called [`vep_core::stderr::drop_refused_writes`], as the `vep` binary does;
/// in any other process the line is written as `eprintln!` writes it.
pub(crate) fn print_on_stderr(line: &str) {
    vep_core::stderr::write(format_args!("{line}\n"));
}

/// [`print_on_stderr`] of `WARNING: <text>`, the shape of most such lines.
pub(crate) fn warn_on_stderr(text: &str) {
    print_on_stderr(&format!("WARNING: {text}"));
}

/// Type-erased container for data returned by [`BuiltinPlugin::prefetch`].
///
/// Each plugin stores its own concrete type (e.g., `BatchQueryResult`) inside
/// the `Box<dyn Any>` and downcasts it in [`BuiltinPlugin::annotate`]. This
/// keeps the trait free of backend-specific types.
#[non_exhaustive]
#[doc(hidden)]
pub struct PrefetchData(pub Box<dyn std::any::Any + Send>);

impl PrefetchData {
    /// Create a new `PrefetchData` wrapping a concrete value.
    pub fn new<T: std::any::Any + Send + 'static>(data: T) -> Self {
        Self(Box::new(data))
    }

    /// Downcast the inner data to a concrete type.
    ///
    /// Returns `Err(PluginError::Run)` if the downcast fails, which indicates
    /// a defect in the plugin (prefetch/annotate type mismatch).
    pub fn downcast<T: 'static>(self) -> Result<T, PluginError> {
        self.0
            .downcast::<T>()
            .map(|b| *b)
            .map_err(|_| PluginError::Run("prefetch data type mismatch".into()))
    }
}

pub(crate) mod sealed {
    /// The supertrait of `BuiltinPlugin` that code outside the crate cannot
    /// name: a type implements `BuiltinPlugin` only where it implements this.
    pub trait Sealed {}
}

/// Trait that all built-in VEP plugins implement.
///
/// Plugins are compiled directly into the binary and operate on native Rust
/// types. There are three execution strategies (the registry picks the best
/// one automatically):
///
/// 1. **Two-phase (prefetch + annotate):** The registry calls [`prefetch`]
///    on all plugins in parallel (via rayon), then calls [`annotate`]
///    sequentially to write results into variants. Tabix plugins use this
///    to overlap I/O from independent data files.
///
/// 2. **Single-phase (run_batch):** Fallback for plugins that return `None`
///    from [`prefetch`]. Called sequentially with `&mut` access to variants.
///
/// 3. **Per-transcript (run):** Default implementation of [`run_batch`]
///    dispatches to [`run`] for each transcript consequence.
///
/// The trait cannot be implemented outside this crate, because the registry
/// builds plugins by name alone and takes none from a caller.
///
/// [`prefetch`]: BuiltinPlugin::prefetch
/// [`annotate`]: BuiltinPlugin::annotate
/// [`run_batch`]: BuiltinPlugin::run_batch
/// [`run`]: BuiltinPlugin::run
#[doc(hidden)]
pub trait BuiltinPlugin: sealed::Sealed + Send + Sync {
    /// Plugin name matching the `--plugin` CLI argument (e.g., "CADD").
    fn name(&self) -> &str;

    /// Header info for output columns: `(field_name, description)` pairs.
    fn header_info(&self) -> Vec<(String, String)>;

    /// The kinds of row the plugin runs on, Ensembl VEP's `feature_types`
    /// (`BaseVepPlugin.pm:99` defaults it to `Transcript`, `BaseVepPlugin.pm:229-233`
    /// reads it): a plugin writes nothing on a row of another kind
    /// (`OutputFactory.pm:2232`), and a value it writes for the whole variant
    /// prints on the rows of these kinds alone. A plugin whose Ensembl VEP module
    /// declares `['Feature', 'Intergenic']` returns [`RowKind::ALL`].
    fn feature_types(&self) -> &[RowKind] {
        &RowKind::DEFAULT
    }

    /// Whether the plugin runs on `consequence`'s row
    /// ([`feature_types`](Self::feature_types)).
    fn annotates_row(&self, consequence: &TranscriptConsequence) -> bool {
        self.feature_types()
            .contains(&RowKind::from(consequence.feature_type))
    }

    /// Whether the plugin runs on structural variants as well as sequence
    /// variants, Ensembl VEP's `variant_feature_types` (`BaseVepPlugin.pm:100`
    /// defaults it to `['VariationFeature']`; `OutputFactory.pm:2229` checks it
    /// before the feature kind, through the `isa` test at
    /// `BaseVepPlugin.pm:290`). A StructuralVariationFeature is a
    /// BaseVariationFeature and not a VariationFeature
    /// (`StructuralVariationFeature.pm:87`, `VariationFeature.pm:118`), so a
    /// plugin with the default writes nothing on a structural variant's rows.
    /// A plugin whose Ensembl VEP module lists `StructuralVariationFeature` or
    /// `BaseVariationFeature` returns `true`.
    fn annotates_structural_variants(&self) -> bool {
        false
    }

    /// Whether the plugin runs on `variant`'s rows
    /// ([`annotates_structural_variants`](Self::annotates_structural_variants)):
    /// every plugin runs on a sequence variant. A plugin that overrides
    /// [`prefetch`](Self::prefetch), [`annotate`](Self::annotate) or
    /// [`run_batch`](Self::run_batch) applies it to each variant it reads.
    fn annotates_variant(&self, variant: &InputVariant) -> bool {
        self.annotates_structural_variants() || !variant.is_structural
    }

    /// Initialize the plugin with CLI parameters.
    ///
    /// Parameters come from `--plugin Name,param1,param2,...` after the name.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError>;

    /// The URL a `--plugin` word names as a data file under the plugin's own
    /// parameter grammar, read by the registry before [`init`](Self::init),
    /// which refuses the plugin with `remote data files are not supported`
    /// naming it. The default reading, `default_remote_data_file`, is the
    /// value of `key=value` or the whole positional word, the grammar of Ensembl
    /// VEP's `params_to_hash` and of the plugins that read by position. A plugin
    /// whose module reads a word another way overrides it (LoFTEE's `key:value`
    /// names no URL, since its module drops a word of three or more colon
    /// fields), so a URL is refused exactly where the module would open it.
    fn remote_data_file<'a>(&self, param: &'a str) -> Option<&'a str> {
        default_remote_data_file(param)
    }

    /// Inject the pipeline's shared reference FASTA, if one was provided via
    /// `--fasta`. Called by the registry after [`init`](BuiltinPlugin::init).
    ///
    /// Default is a no-op; the plugins that read the reference (LoFTEE for
    /// intron boundary motifs, MaxEntScan for splice-site windows, UTRAnnotator
    /// for the 5' UTR and coding sequence) use it, and a plugin that reads a
    /// table through an annotation store hands it to the store, which compares
    /// the table's mitochondrial rows with it on GRCh37
    /// ([`crate::tabix::MitochondrialFrame`]). Sharing the already-loaded
    /// `Arc` avoids opening the FASTA twice.
    fn set_reference_fasta(&mut self, _fasta: Option<Arc<vep_fasta::IndexedFasta>>) {}

    /// Whether the plugin cannot annotate without the reference FASTA. The
    /// registry refuses to load such a plugin when the run has none, before
    /// [`init`](Self::init) runs, so the answer cannot depend on a parameter.
    /// MaxEntScan and UTRAnnotator declare it: Ensembl VEP's MaxEntScan dies
    /// offline without `--fasta` (MaxEntScan.pm:94), while its UTRAnnotator
    /// prints a notice and does not die (UTRAnnotator.pm:88). A plugin that
    /// degrades gracefully (LoFTEE) keeps the default.
    fn requires_reference_fasta(&self) -> bool {
        false
    }

    /// Whether the plugin cannot annotate without the HGVS notations. The
    /// registry refuses to load such a plugin when the run computes none,
    /// before [`init`](Self::init) runs, as HGVSIntronOffset reads its offsets
    /// from the HGVSc coordinates (HGVSIntronOffset.pm:62-65) that `--hgvs`
    /// builds.
    fn requires_hgvs(&self) -> bool {
        false
    }

    /// Inject the run's target assembly, parsed from the run's assembly name
    /// (the cache's own, else the value of `--assembly`). Called by the
    /// registry **before** [`init`](BuiltinPlugin::init), so `init` can use it
    /// while opening data files.
    ///
    /// Default is a no-op. Plugins that need it fall into two groups:
    ///
    /// - **Dual-coordinate files**: REVEL ships one file holding both
    ///   `hg19_pos` and `grch38_pos`, indexed on a different column per
    ///   assembly. Without the assembly the plugin cannot know which column
    ///   holds the coordinates tabix indexed, and reads the wrong one.
    /// - **Per-assembly readers**: LoFTEE reads its GERP scores with `tabix` on
    ///   GRCh37 and from a bigWig on GRCh38, as each branch of Ensembl VEP's
    ///   plugin does (`LoF.pm` lines 138 to 141 of `master`, line 130 of
    ///   `grch38`), with different END_TRUNC cutoffs.
    ///
    /// Plugins keyed on gene symbol (LoFtool, pLI) genuinely ignore this: one file
    /// serves both assemblies.
    ///
    /// Several plugins whose databases are pre-split per assembly override it, to
    /// *validate* that the supplied file matches the requested build rather than
    /// annotating from the wrong genome. A pre-split file is not self-evidently
    /// correct: the two copies share a column layout, so pointing a run at the wrong
    /// one resolves contigs, matches nothing, and exits 0. Each uses whichever
    /// in-file signal that database carries: CADD the assembly's name in any
    /// header line (CADD.pm:151-157, checked when the plugin opens its file),
    /// dbscSNV its coordinate column names, AlphaMissense its per-row `genome`
    /// column (`tabix::first_data_row_field`).
    ///
    /// gnomADc deliberately does not: its two files are different
    /// releases (v2.1 genomes on GRCh37, v4.0 exomes on GRCh38) whose headers the
    /// plugin already reads to discover its output fields, and neither carries a
    /// build marker. There is no in-file signal to validate against, so a check here
    /// would have to infer the build from a contig-name convention, exactly the
    /// heuristic that would break on the next release.
    fn set_assembly(&mut self, _assembly: Option<vep_core::assembly::Assembly>) {}

    /// Inject the run's output format. Called by the registry **before**
    /// [`init`](BuiltinPlugin::init), as Ensembl VEP's plugins read
    /// `$self->{config}->{output_format}` in their constructors.
    ///
    /// Default is a no-op. A plugin whose Ensembl VEP module shapes its fields
    /// by format overrides it: SpliceAI writes one pipe-joined field in the
    /// default and tab formats and nine fields under VCF (SpliceAI.pm:173-175,
    /// 208, 275-278), and AlphaMissense and SpliceAI group their fields under
    /// the plugin's name in JSON (AlphaMissense.pm:281-282, SpliceAI.pm:316).
    fn set_output_format(&mut self, _format: OutputFormat) {}

    /// Inject the run's transcript set: whether `--refseq` or `--merged` was
    /// given. Called by the registry **before** [`init`](BuiltinPlugin::init).
    ///
    /// Default is a no-op. The plugins that match records to transcripts by
    /// Ensembl stable id warn under either flag that the match may fail
    /// (REVEL.pm:144-146, AlphaMissense.pm:188-190, pLI.pm:187-189 and
    /// SpliceVault.pm:129-131, each a bare `warn`, and dbNSFP.pm:617-619; all
    /// five printed here with `warn_on_stderr`), since RefSeq ids never match
    /// an Ensembl id.
    fn set_transcript_set(&mut self, _refseq: bool, _merged: bool) {}

    /// Inject the run's plugins directories, in search order (the last
    /// `--dir_plugins` match first). Called by the registry **before**
    /// [`init`](BuiltinPlugin::init).
    ///
    /// Default is a no-op. A plugin whose Ensembl VEP module reads a bundled
    /// data file from the directory the module loaded from (`$INC{'LoFtool.pm'}`,
    /// LoFtool.pm:74-78; pLI.pm:115-119) overrides it and, given no file
    /// parameter, reads the first directory that holds the file.
    fn set_plugin_dirs(&mut self, _dirs: &[String]) {}

    /// I/O phase: fetch annotation data for a batch of variants.
    ///
    /// This method takes a **shared** reference to the variant buffer, so the
    /// registry can call it on multiple plugins in parallel. The returned
    /// [`PrefetchData`] is an opaque container that will be passed to
    /// [`annotate`](BuiltinPlugin::annotate) in the sequential phase.
    ///
    /// Return `Ok(None)` to opt out of two-phase execution; the registry
    /// will fall back to calling [`run_batch`](BuiltinPlugin::run_batch).
    fn prefetch(&self, _variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        Ok(None)
    }

    /// Whether [`prefetch`](BuiltinPlugin::prefetch) reads the variants'
    /// consequence rows, so it can only run once the batch is annotated. The
    /// default, `false`, says the prefetch reads the records' coordinates,
    /// alleles and class alone, and the registry may run it before the batch
    /// is annotated, overlapping its reads with the annotation of the batch
    /// before ([`BuiltinRegistry::prefetch_ahead`](crate::BuiltinRegistry::prefetch_ahead)).
    /// A plugin that queries only the records with a row of a given consequence
    /// (ClinPred reads the rows it scores) returns `true`.
    fn prefetch_reads_consequences(&self) -> bool {
        false
    }

    /// Annotation phase: apply prefetched data to variants.
    ///
    /// Called sequentially after all [`prefetch`](BuiltinPlugin::prefetch)
    /// calls complete. The `data` argument is the value returned from
    /// `prefetch` for this plugin. An `annotate` that fails writes none of the
    /// plugin's fields on the variant it fails on, as a failed `prefetch`
    /// leaves none: the registry's warning states that the plugin annotated
    /// nothing on the variants it failed on. A plugin that returns `true` from
    /// [`annotate_fails_per_row`](BuiltinPlugin::annotate_fails_per_row)
    /// instead writes every row it does not fail on, and the warning states
    /// that footprint.
    fn annotate(
        &self,
        _variants: &mut [InputVariant],
        _data: PrefetchData,
    ) -> Result<(), PluginError> {
        Ok(())
    }

    /// Whether a failure of [`annotate`](BuiltinPlugin::annotate) costs the
    /// rows it failed on alone, as a [`run_batch`](BuiltinPlugin::run_batch)
    /// failure does, rather than every field of the variant. A plugin whose
    /// Ensembl VEP module dies per output line, so that the lines matched
    /// before the failing row keep its fields, writes every row it does not
    /// fail on and returns `true`; the registry's warning then states that
    /// footprint.
    fn annotate_fails_per_row(&self) -> bool {
        false
    }

    /// Batch-annotate a buffer of variants (single-phase fallback).
    ///
    /// Called once per buffer when the plugin does not implement the two-phase
    /// prefetch/annotate protocol. Tabix plugins should prefer implementing
    /// [`prefetch`](BuiltinPlugin::prefetch) + [`annotate`](BuiltinPlugin::annotate)
    /// instead.
    ///
    /// Default implementation calls [`run`](BuiltinPlugin::run), through
    /// [`run_naming_empty_strings`](BuiltinPlugin::run_naming_empty_strings),
    /// for each variant/consequence pair of the variants the plugin annotates
    /// ([`annotates_variant`](BuiltinPlugin::annotates_variant)), every row
    /// of the batch, and reports the first failure after the loop, so a row
    /// whose `run` fails loses the plugin's fields alone, as the row does under
    /// Ensembl VEP, whose `run_plugins` evals each row's plugin call on its own
    /// (`OutputFactory.pm` lines 2234-2253). An override that fails keeps the
    /// same footprint: the registry's warning states that the rows the plugin
    /// failed on alone lack its fields.
    ///
    /// A warning a row's `run` raises goes on under the input record of the
    /// row's variant, unless it names a record itself, so a run that holds a
    /// batch's warnings reports it with that record
    /// (`vep_core::warning_hold`).
    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        let mut first_failure = None;
        for variant in variants.iter_mut() {
            if let Err(error) = run_rows(self, variant) {
                first_failure.get_or_insert(error);
            }
        }
        match first_failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    /// Annotate a single transcript consequence (per-transcript plugins).
    ///
    /// Returns key-value pairs to add to the Extra column. Return an empty
    /// map if the plugin has no annotation for this consequence.
    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        Ok(IndexMap::new())
    }

    /// [`run`](BuiltinPlugin::run) with the keys of its result whose empty
    /// value is Ensembl VEP's empty string, a table cell the row has and
    /// leaves empty, where every other empty value is its undefined value.
    /// The default [`run_batch`](BuiltinPlugin::run_batch) names those keys
    /// in the row's `plugin_empty_strings`, as a plugin's own batch route
    /// does, and a caller that writes the fields itself names them through
    /// [`TranscriptConsequence::set_plugin_value`]. The default names none.
    fn run_naming_empty_strings(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<(IndexMap<String, String>, Vec<String>), PluginError> {
        Ok((self.run(consequence, variant)?, Vec::new()))
    }
}

/// The default [`BuiltinPlugin::run_batch`] over one variant: every row the
/// plugin annotates gets the fields its `run` returns, and the first failure
/// comes back after the variant's last row. The warnings the rows raise go on
/// under the variant's input record (`about_record`).
fn run_rows<P: BuiltinPlugin + ?Sized>(
    plugin: &P,
    variant: &mut InputVariant,
) -> Result<(), PluginError> {
    if !plugin.annotates_variant(variant) {
        return Ok(());
    }
    let mut first_failure = None;
    about_record(variant.input_record, || {
        for tc_idx in 0..variant.transcript_consequences.len() {
            if !plugin.annotates_row(&variant.transcript_consequences[tc_idx]) {
                continue;
            }
            let row = &variant.transcript_consequences[tc_idx];
            match plugin.run_naming_empty_strings(row, variant) {
                Ok((fields, empty_strings)) => {
                    let row = &mut variant.transcript_consequences[tc_idx];
                    for (key, value) in fields {
                        let empty_string = empty_strings.contains(&key);
                        row.set_plugin_value(key, value, empty_string);
                    }
                }
                Err(error) => {
                    first_failure.get_or_insert(error);
                }
            }
        }
    });
    first_failure.map_or(Ok(()), Err)
}

/// The default [`BuiltinPlugin::run_batch`] with the variants spread over the
/// current rayon pool, for a plugin whose `run` reads shared data alone: the
/// fields, the warnings and the failure returned are the serial default's,
/// the first failure in variant order (`map_in_order`).
pub(crate) fn run_batch_in_parallel<P: BuiltinPlugin + ?Sized>(
    plugin: &P,
    variants: &mut [InputVariant],
) -> Result<(), PluginError> {
    map_in_order(variants.par_iter_mut(), |variant| run_rows(plugin, variant))
        .into_iter()
        .find(Result::is_err)
        .unwrap_or(Ok(()))
}

/// The default reading of [`BuiltinPlugin::remote_data_file`]: the URL
/// (`scheme://...`) a word names as a whole positional word, query and all,
/// or as the value of `key=value`, which `params_to_hash` cuts at the value's
/// first `=` (`perl_key_value`: `file=https://h/f?a=b` names `https://h/f?a`,
/// the file Ensembl VEP's plugin opens).
pub(crate) fn default_remote_data_file(param: &str) -> Option<&str> {
    use crate::custom::grammar::is_remote;
    if is_remote(param) {
        return Some(param);
    }
    crate::params::perl_key_value(param)
        .map(|(_, value)| value)
        .filter(|value| is_remote(value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::FeatureType;

    /// A per-consequence plugin whose `run` marks every row it is handed.
    struct MarkerPlugin;

    impl crate::traits::sealed::Sealed for MarkerPlugin {}

    impl BuiltinPlugin for MarkerPlugin {
        fn name(&self) -> &str {
            "Marker"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            let mut result = IndexMap::new();
            result.insert("Marker".to_string(), "1".to_string());
            Ok(result)
        }
    }

    fn row(feature_type: FeatureType) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.feature_type = feature_type;
        consequence
    }

    /// A plugin that names no kinds runs on transcript rows alone, Ensembl
    /// VEP's default `feature_types` (`BaseVepPlugin.pm:99`).
    #[test]
    fn annotates_row_accepts_transcript_rows_alone() {
        assert_eq!(MarkerPlugin.feature_types(), [RowKind::Transcript]);
        assert!(MarkerPlugin.annotates_row(&row(FeatureType::Transcript)));
        assert!(!MarkerPlugin.annotates_row(&row(FeatureType::RegulatoryFeature)));
        assert!(!MarkerPlugin.annotates_row(&row(FeatureType::MotifFeature)));
    }

    /// Ensembl VEP hands a row to a plugin only when the plugin's
    /// `feature_types` names the row's feature class (OutputFactory.pm:2232,
    /// BaseVepPlugin.pm `check_feature_type`), so the per-consequence `run`
    /// path writes on the transcript row of a record and on neither of its
    /// regulatory rows.
    #[test]
    fn default_run_batch_writes_on_transcript_rows_alone() {
        let mut variants = vec![InputVariant::new(
            "21".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        variants[0].transcript_consequences.extend([
            row(FeatureType::Transcript),
            row(FeatureType::RegulatoryFeature),
            row(FeatureType::MotifFeature),
        ]);

        MarkerPlugin.run_batch(&mut variants).unwrap();

        let written: Vec<bool> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.contains_key("Marker"))
            .collect();
        assert_eq!(written, [true, false, false]);
        assert!(variants[0].plugin_data.is_empty());
    }

    /// A per-consequence plugin whose `run` raises two warnings on every row,
    /// one naming no record and one naming record 99, as the layer that writes
    /// a run's warnings gives them to the hold of the thread raising them.
    struct WarningRowPlugin;

    impl crate::traits::sealed::Sealed for WarningRowPlugin {}

    impl BuiltinPlugin for WarningRowPlugin {
        fn name(&self) -> &str {
            "WarningRow"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            let _ = vep_core::warning_hold::offer(0, format!("row at {}", variant.start));
            let _ = vep_core::warning_hold::offer(99, format!("names 99 at {}", variant.start));
            Ok(IndexMap::new())
        }
    }

    /// The default `run_batch` hands on what a row's `run` raises under the
    /// input record of the row's variant, in variant order, and leaves a
    /// warning that names a record under that record.
    #[test]
    fn default_run_batch_names_the_variants_record_on_a_rows_warnings() {
        use vep_core::warning_hold::{holding, Held};
        let mut variants: Vec<InputVariant> = [(100u64, 7u64), (200, 8)]
            .into_iter()
            .map(|(position, input_record)| {
                let mut variant = InputVariant::new(
                    "21".into(),
                    position,
                    position,
                    b"A".to_vec(),
                    b"G".to_vec(),
                );
                variant.input_record = input_record;
                variant
                    .transcript_consequences
                    .push(row(FeatureType::Transcript));
                variant
            })
            .collect();

        let (ran, raised) = holding(|| WarningRowPlugin.run_batch(&mut variants));
        ran.unwrap();

        let held = |input_record: u64, text: &str| Held {
            input_record,
            text: text.to_string(),
        };
        assert_eq!(
            raised,
            [
                held(7, "row at 100"),
                held(99, "names 99 at 100"),
                held(8, "row at 200"),
                held(99, "names 99 at 200"),
            ]
        );
    }

    /// A per-consequence plugin whose `run` raises a warning on every row and
    /// fails on the rows of a variant at a multiple of three; it writes the
    /// variant's position on every other row. The earlier the variant, the
    /// longer its rows take, so threads end in the reverse of variant order.
    struct SlowFailingPlugin;

    impl crate::traits::sealed::Sealed for SlowFailingPlugin {}

    impl BuiltinPlugin for SlowFailingPlugin {
        fn name(&self) -> &str {
            "SlowFailing"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            std::thread::sleep(std::time::Duration::from_micros(50 * (40 - variant.start)));
            let _ = vep_core::warning_hold::offer(0, format!("row at {}", variant.start));
            if variant.start % 3 == 0 {
                return Err(PluginError::Run(format!("fails at {}", variant.start)));
            }
            Ok(IndexMap::from([(
                "Position".to_string(),
                variant.start.to_string(),
            )]))
        }
    }

    /// The parallel default writes the fields, hands on the warnings and
    /// returns the failure the serial default does: the first failure in
    /// variant order, whatever thread ran which variant.
    #[test]
    fn the_parallel_default_run_batch_matches_the_serial_one() {
        use vep_core::warning_hold::holding;
        let batch = || -> Vec<InputVariant> {
            (1..40u64)
                .map(|position| {
                    let mut variant = InputVariant::new(
                        "21".into(),
                        position,
                        position,
                        b"A".to_vec(),
                        b"G".to_vec(),
                    );
                    variant.input_record = position;
                    variant
                        .transcript_consequences
                        .extend([row(FeatureType::Transcript), row(FeatureType::Transcript)]);
                    variant
                })
                .collect()
        };
        let mut serial = batch();
        let (serial_result, serial_raised) = holding(|| SlowFailingPlugin.run_batch(&mut serial));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap();
        let mut parallel = batch();
        let (parallel_result, parallel_raised) =
            pool.install(|| holding(|| run_batch_in_parallel(&SlowFailingPlugin, &mut parallel)));
        assert_eq!(
            serial_result.unwrap_err().to_string(),
            "plugin execution failed: fails at 3"
        );
        assert_eq!(
            parallel_result.unwrap_err().to_string(),
            "plugin execution failed: fails at 3"
        );
        assert_eq!(parallel_raised, serial_raised);
        assert_eq!(parallel_raised.len(), 2 * 39);
        let fields = |variants: &[InputVariant]| -> Vec<IndexMap<String, String>> {
            variants
                .iter()
                .flat_map(|v| {
                    v.transcript_consequences
                        .iter()
                        .map(|tc| tc.plugin_data.clone())
                })
                .collect()
        };
        assert_eq!(fields(&parallel), fields(&serial));
    }

    /// The value the default `run_batch` writes on a row stands over an
    /// earlier plugin's under the same key (OutputFactory.pm:2239-2241), so the
    /// earlier value's mark in the row's `plugin_empty_strings` goes with it,
    /// and a mark under another key stays.
    #[test]
    fn default_run_batch_replaces_an_earlier_value_and_its_mark() {
        let mut variants = vec![InputVariant::new(
            "21".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        let mut earlier = row(FeatureType::Transcript);
        earlier.set_plugin_value("Marker".to_string(), String::new(), true);
        earlier.set_plugin_value("Other".to_string(), String::new(), true);
        assert_eq!(earlier.plugin_empty_strings, ["Marker", "Other"]);
        variants[0].transcript_consequences.push(earlier);

        MarkerPlugin.run_batch(&mut variants).unwrap();

        let written = &variants[0].transcript_consequences[0];
        assert_eq!(written.plugin_data["Marker"], "1");
        assert_eq!(written.plugin_empty_strings, ["Other"]);
    }

    /// A per-consequence plugin whose `run` hands back two empty values and a
    /// filled one, and which names `Cell` and `Filled` as empty strings.
    struct EmptyCellPlugin;

    impl crate::traits::sealed::Sealed for EmptyCellPlugin {}

    impl BuiltinPlugin for EmptyCellPlugin {
        fn name(&self) -> &str {
            "EmptyCell"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            Ok(IndexMap::from([
                ("Undefined".to_string(), String::new()),
                ("Cell".to_string(), String::new()),
                ("Filled".to_string(), "1".to_string()),
            ]))
        }
        fn run_naming_empty_strings(
            &self,
            consequence: &TranscriptConsequence,
            variant: &InputVariant,
        ) -> Result<(IndexMap<String, String>, Vec<String>), PluginError> {
            let named = ["Cell", "Filled"].map(String::from).to_vec();
            Ok((self.run(consequence, variant)?, named))
        }
    }

    /// The default `run_batch` names in the row's `plugin_empty_strings` the
    /// keys `run_naming_empty_strings` names whose value is empty: a named key
    /// with a value and an empty value that is not named stay unmarked, and a
    /// plugin that names none marks none.
    #[test]
    fn default_run_batch_names_the_empty_strings_a_plugin_names() {
        let record = || {
            let mut variant =
                InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
            variant
                .transcript_consequences
                .push(row(FeatureType::Transcript));
            vec![variant]
        };

        let mut variants = record();
        EmptyCellPlugin.run_batch(&mut variants).unwrap();
        let written = &variants[0].transcript_consequences[0];
        assert_eq!(
            written.plugin_data,
            IndexMap::from([
                ("Undefined".to_string(), String::new()),
                ("Cell".to_string(), String::new()),
                ("Filled".to_string(), "1".to_string()),
            ])
        );
        assert_eq!(written.plugin_empty_strings, ["Cell"]);

        let mut variants = record();
        MarkerPlugin.run_batch(&mut variants).unwrap();
        let written = &variants[0].transcript_consequences[0];
        assert_eq!(written.plugin_data["Marker"], "1");
        assert!(written.plugin_empty_strings.is_empty());
    }

    /// A per-consequence plugin whose Ensembl VEP module would list
    /// `StructuralVariationFeature` among its `variant_feature_types`.
    struct StructuralMarkerPlugin;

    impl crate::traits::sealed::Sealed for StructuralMarkerPlugin {}

    impl BuiltinPlugin for StructuralMarkerPlugin {
        fn name(&self) -> &str {
            "StructuralMarker"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn annotates_structural_variants(&self) -> bool {
            true
        }
        fn run(
            &self,
            consequence: &TranscriptConsequence,
            variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            MarkerPlugin.run(consequence, variant)
        }
    }

    /// A deletion written as the symbolic allele `<DEL>` over `start..=end`,
    /// with one transcript row.
    fn structural_deletion(start: u64, end: u64) -> InputVariant {
        let mut variant =
            InputVariant::new("21".into(), start, end, b"N".to_vec(), b"<DEL>".to_vec());
        variant.is_structural = true;
        variant
            .transcript_consequences
            .push(row(FeatureType::Transcript));
        variant
    }

    /// Ensembl VEP 116.3 runs a plugin on a line only when the plugin's
    /// `variant_feature_types` names the variant's class
    /// (`OutputFactory.pm:2229` `check_variant_feature_type(ref($vf))`), and
    /// the default is `['VariationFeature']` alone (`BaseVepPlugin.pm:100`), a
    /// class a StructuralVariationFeature is not
    /// (`StructuralVariationFeature.pm:87`): a plugin that names no kinds
    /// writes nothing on a structural variant's transcript row and still
    /// writes on the sequence variant beside it.
    #[test]
    fn default_run_batch_skips_structural_variants() {
        assert!(!MarkerPlugin.annotates_structural_variants());
        let mut sequence = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        sequence
            .transcript_consequences
            .push(row(FeatureType::Transcript));
        let mut variants = vec![structural_deletion(100, 500), sequence];
        assert!(!MarkerPlugin.annotates_variant(&variants[0]));
        assert!(MarkerPlugin.annotates_variant(&variants[1]));

        MarkerPlugin.run_batch(&mut variants).unwrap();

        assert!(variants[0].transcript_consequences[0]
            .plugin_data
            .is_empty());
        assert!(variants[1].transcript_consequences[0]
            .plugin_data
            .contains_key("Marker"));
    }

    /// A plugin whose Ensembl VEP module lists `StructuralVariationFeature`
    /// (CADD.pm:186-188, DosageSensitivity.pm:113-115, UTRAnnotator.pm:80-82)
    /// or `BaseVariationFeature` (TSSDistance.pm:73-75) passes the `isa` test
    /// of `BaseVepPlugin.pm:290` for a StructuralVariationFeature, so the
    /// default `run_batch` hands it the structural variant's rows.
    #[test]
    fn a_plugin_declaring_structural_variants_runs_on_them() {
        let mut variants = vec![structural_deletion(100, 500)];
        assert!(StructuralMarkerPlugin.annotates_variant(&variants[0]));

        StructuralMarkerPlugin.run_batch(&mut variants).unwrap();

        assert_eq!(
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("Marker")
                .map(String::as_str),
            Some("1")
        );
    }

    /// A per-consequence plugin whose `run` fails, with the row's `exon` text
    /// as the error, on a row that carries one, and marks every other row.
    struct ExonFailsPlugin;

    impl crate::traits::sealed::Sealed for ExonFailsPlugin {}

    impl BuiltinPlugin for ExonFailsPlugin {
        fn name(&self) -> &str {
            "ExonFails"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run(
            &self,
            consequence: &TranscriptConsequence,
            variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            match &consequence.exon {
                Some(exon) => Err(PluginError::Run(exon.clone())),
                None => MarkerPlugin.run(consequence, variant),
            }
        }
    }

    /// Ensembl VEP's `run_plugins` evals each row's plugin call on its own
    /// (`OutputFactory.pm` lines 2234-2253), so a row whose call dies loses the
    /// plugin's fields alone: the default `run_batch` runs every row, writes
    /// on the rows before and after a failing one and on the next variant's,
    /// and returns the first failure.
    #[test]
    fn default_run_batch_runs_every_row_and_returns_the_first_failure() {
        let failing = |exon: &str| {
            let mut consequence = row(FeatureType::Transcript);
            consequence.exon = Some(exon.to_string());
            consequence
        };
        let mut first = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        first.transcript_consequences.extend([
            row(FeatureType::Transcript),
            failing("2/5"),
            row(FeatureType::Transcript),
            failing("3/5"),
        ]);
        let mut second = InputVariant::new("21".into(), 200, 200, b"C".to_vec(), b"T".to_vec());
        second
            .transcript_consequences
            .push(row(FeatureType::Transcript));
        let mut variants = vec![first, second];

        let error = ExonFailsPlugin
            .run_batch(&mut variants)
            .expect_err("the first failing row's error");
        assert!(
            matches!(&error, PluginError::Run(text) if text == "2/5"),
            "{error}"
        );

        let marked = |variant: &InputVariant| -> Vec<bool> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.plugin_data.contains_key("Marker"))
                .collect()
        };
        assert_eq!(marked(&variants[0]), [true, false, true, false]);
        assert_eq!(marked(&variants[1]), [true]);
    }
}
