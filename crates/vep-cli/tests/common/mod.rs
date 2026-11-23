// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Shared helpers for the golden-corpus and format-parity tests: corpus
//! discovery, running the `vep` binary on a corpus, parsing each output format
//! into comparable entries, and rendering every mismatch in one report.
//!
//! A corpus lives at `tests/golden/<release>/<name>/` with `variants.vcf`,
//! `manifest.json`, its cache and `expected/{default.txt,tab.txt,vcf.vcf,
//! json.jsonl}[.gz]`; `<name>` is the assembly, optionally suffixed
//! (`GRCh37-hgvs`), the manifest's `assembly` is what the run passes to
//! `--assembly` and its `species` (`homo_sapiens` when absent) what it passes
//! to `--species`. The cache is `json_cache/`, passed with `--json_cache`, unless
//! the manifest's `cache_layout` is `ensembl`: then it sits at
//! `cache_root/<species>/<release>_<assembly>/` (the species directory suffixed
//! `_refseq` or `_merged` when the manifest's flags carry `--refseq` or
//! `--merged`) and every run passes
//! `--cache --dir_cache <corpus>/cache_root --cache_version <release>`, so the
//! binary resolves the directory as Ensembl VEP does. The manifest may also
//! carry `flags` (extra command-line flags, such as `--hgvs`; a word may carry
//! `{corpus}`, the corpus directory, so a flag names a file beside the manifest
//! without an absolute path, as a `runs[].argv` token may, and `{plugins}`, the
//! plugin directory each engine reads, here the corpus's `plugin_data` and in
//! the reference run the directory of Ensembl VEP's plugin modules;
//! [`corpus_flags`])
//! and `fasta` (a reference beside the manifest, bgzipped with its `.fai` and
//! `.gzi`, passed to `--fasta` as it is); a word may also carry `{generated}`,
//! the directory the manifest's `generated` block writes before the run: its
//! `command`, run from the repository root, and the `sha256` of every file it
//! writes, for a fixture the corpus reproduces from a seed instead of committing
//! ([`generated_dir`]). With `fasta_in_cache` true the reference and its
//! indexes are placed inside a staged copy of the cache instead, which the run
//! names as its cache with no `--fasta`, so the run relies on the binary's
//! detection of an indexed `.fa` or `.fa.gz` in the cache directory. Its
//! `divergences` list names the (Location, Allele, Feature) keys whose
//! consequence terms are documented to differ between VEP and vep-rs (with the
//! corpus record ordinals they belong to), and `vep_rs_only_tuples` the keys
//! only vep-rs emits, each classified; both are compared against their
//! documented shape rather than for equality. An item whose
//! `vep_rs_consequence_sets` hold VEP's own set beside another documents a key
//! several records share: a row on which vep-rs writes VEP's set is compared
//! for equality like any other, and the rows that differ are held to the
//! item's sets ([`compare_entries`]). A `divergences` item's
//! `vep_rs_consequence_sets` are vep-rs's sets in a run without a reference
//! sequence; an item may state its sets in a run that has one under
//! `vep_rs_consequence_sets_with_reference_sequence`, and such a run is held to
//! those: where they are VEP's own set the key is compared for equality like
//! any other ([`Documented::for_run`], [`has_reference_sequence`]). The class
//! `reference_skipped_record` names a record the reference dropped before
//! annotation (its `reference_rows` is 0 and `reference_warning` quotes the
//! reason): the reference writes nothing for it, so the record-level VCF line
//! and JSON object vep-rs writes for it are counted, not compared, while every
//! consequence entry on it must still be a documented vep-rs-only key.
//! `field_divergences` names keys whose terms agree but whose named field (an
//! HGVS string) is documented to differ, with the value vep-rs prints; the key
//! carries the row's consequence set, so a record sharing a Location and Allele
//! with the documented one keeps its own rows compared exactly, and every other
//! field of a documented key compares exactly. An item that states `rows_alike`
//! documents a key and consequence set several records share: a row of it on
//! which both engines print one value of the field is compared for equality like
//! any other, and the rows that differ are held to the item's value
//! ([`Documented::field_shares_its_key`]); the rows found alike are held to the
//! count the item states. A run under `--nearest transcript`
//! compares its `NEAREST` values as the coordinates of the named transcripts'
//! start sites ([`nearest_as_coordinates`]). `masks` lists the mask classes
//! applied to every run of the corpus before the comparison: under
//! [`REGULATORY_FEATURE_ROWS`] the reference loses its regulatory rows
//! ([`apply_masks`]), under [`ORDER_ONLY`] the `DOMAINS`, `FREQS`, `ZYG` and
//! `NearestExonJB` items, the `CLIN_SIG` lists and the `5UTR_annotation` pair
//! groups of both sides compare as a multiset ([`list_columns_as_multisets`]), under
//! [`UNLOADED_MODULE_PLUGIN`] the reference's run-time warnings about a plugin
//! module vep-rs cannot load stand against vep-rs's load warning for it
//! ([`compare_warnings`]), and under [`CONSERVATION_SHORT_NAME_SPELLING`] and
//! [`PC_DESCRIPTION_TEXT`] the reference text is read with a `--phyloP` source's
//! key spelt `phyloP<N>way` and with vep-rs's description of a `PC` field
//! ([`masked_reference_text`]). A class of [`CORRECTED_BEHAVIOUR_CLASSES`]
//! changes no comparison here: its rows are the documented items of the three
//! lists, and naming it in `masks` states which corrected behaviour they belong
//! to.
//!
//! A corpus is run once per [`Run`]: the manifest's own run (its `flags`,
//! `variants.vcf`, every format, `expected/`) and once per `flag_sets` entry,
//! which names its own `flags` (the whole list for that run), `input` (a file in
//! the corpus, `variants.vcf` by default), `stdin` (the input piped on standard
//! input with `-i STDIN`), `formats` (the formats compared, all four by default),
//! `expected_dir` (`expected_<name>` by default), `unordered` (the entries
//! compared without their order, for a set whose rows the reference emits in
//! hash order; entries sharing a key then pair by what they hold,
//! [`number_shared_keys_by_content`]) and `perl_mode` (the manifest's own, `offline` when it names
//! none; `database` for a reference Ensembl VEP wrote connected to its
//! database, whose default and tab headers compare by their column
//! descriptions alone and whose VCF `##VEP=` line by the tokens it names; `gtf`
//! or `gff` for one it wrote reading its transcripts from the file given to
//! `--gtf` or `--gff` while vep-rs reads the cache built from that file, whose
//! headers compare the same way once the source's `SOURCE` column and header
//! lines are set aside, [`GXF_PERL_MODES`]). An entry carrying its own
//! `divergences` documents its run; one without inherits the manifest's lists.
//! An entry carrying its own `masks` lists the classes of its run; one without
//! inherits the manifest's.
//!
//! A manifest may also carry `runs`, command lines executed after the default
//! run whose expectations are about the process rather than the annotation: each
//! [`ManifestRun`] names `argv` (with the placeholders [`expand_argv`] lists),
//! `env`, `stdin` and an `expect` block of exit status, output file compared as
//! the default run is, text standard output starts with or contains or standard
//! error contains, and the warnings and skipped-variants files compared as
//! entry sets ([`warning_entries`]).
//!
//! Every run's warnings file is compared with the warnings the reference wrote
//! for the corpus ([`compare_warnings`]): the manifest's `records[].reference_warning`
//! strings, each at its record's ordinal, and the unnumbered texts of a
//! `reference_warnings` list (the flag set's own, else the manifest's;
//! [`run_reference_warnings`]), with paths masked and the quoted input line
//! left out. The line vep-rs writes when a run asks for known variants from a
//! cache without variation data has no counterpart in the reference and is set
//! aside on a corpus whose cache holds none.
//!
//! The VCF and JSON formats carry no `Location`, so their entries are keyed by
//! the corpus record ordinal, found by matching each output record to its
//! input line by the first five columns (of the VCF line, of the JSON object's
//! `input`); the VCF output of an input that is not a VCF carries no such line,
//! so its records are the input records in turn, one per line.

#![allow(dead_code)]

use std::borrow::Borrow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use indexmap::IndexMap;

use vep_core::cache_manifest::{CacheManifest, MANIFEST_FILE};
use vep_core::consequence::Consequence;

/// One corpus directory.
pub struct Corpus {
    pub release: String,
    /// The directory name: the assembly, optionally suffixed.
    pub name: String,
    /// The assembly passed to `--assembly`: the manifest's `assembly`, else the name.
    pub assembly: String,
    /// The species passed to `--species`, and the cache directory's name under the
    /// Ensembl layout: the manifest's `species`, else [`DEFAULT_SPECIES`].
    pub species: String,
    pub dir: PathBuf,
    pub manifest: serde_json::Value,
}

/// The species of a manifest that names none.
pub const DEFAULT_SPECIES: &str = "homo_sapiens";

impl Corpus {
    /// Whether the cache sits in Ensembl VEP's layout under `cache_root/`
    /// (`cache_layout: "ensembl"`) rather than at `json_cache/`.
    pub fn ensembl_cache_layout(&self) -> bool {
        match self.manifest.get("cache_layout") {
            None => false,
            Some(serde_json::Value::String(layout)) if layout == "ensembl" => true,
            Some(other) => panic!(
                "{}/{}: cache_layout is {other}, not \"ensembl\"",
                self.release, self.name
            ),
        }
    }

    /// The cache directory itself: `json_cache/`, or the leaf of the Ensembl layout
    /// under [`Self::species_dir`].
    pub fn cache_dir(&self) -> PathBuf {
        if self.ensembl_cache_layout() {
            self.cache_root()
                .join(self.species_dir())
                .join(format!("{}_{}", self.release, self.assembly))
        } else {
            self.dir.join("json_cache")
        }
    }

    /// The transcript-set flags of the manifest's `flags`, `--refseq` and
    /// `--merged`, in that order. A corpus cut from a RefSeq or merged cache
    /// declares the set in its `info.json`, which the binary holds a run's flags
    /// to, so every run of the corpus carries them.
    pub fn transcript_set_flags(&self) -> Vec<String> {
        let flags = strings(&self.manifest["flags"]);
        ["--refseq", "--merged"]
            .into_iter()
            .filter(|flag| flags.iter().any(|f| f == flag))
            .map(str::to_string)
            .collect()
    }

    /// The arguments a run of a RefSeq or merged corpus needs beyond its cache:
    /// [`Self::transcript_set_flags`] and, since a BAM-edited cache switches
    /// `--use_transcript_ref` on and reads the reference, `--fasta` naming the
    /// corpus FASTA; empty for every other corpus.
    pub fn transcript_set_args(&self) -> Vec<String> {
        let mut args = self.transcript_set_flags();
        if !args.is_empty() {
            if let Some(fasta) = self.manifest["fasta"].as_str() {
                args.push("--fasta".to_string());
                args.push(self.dir.join(fasta).display().to_string());
            }
        }
        args
    }

    /// The species directory of the Ensembl layout: the species, with `_refseq`
    /// and then `_merged` appended for each of `--refseq` and `--merged` the
    /// manifest's flags carry, the directory the binary resolves from the same
    /// flags as Ensembl VEP does (`CacheDir.pm:312`).
    pub fn species_dir(&self) -> String {
        let mut dir = self.species.clone();
        for flag in self.transcript_set_flags() {
            dir.push('_');
            dir.push_str(flag.trim_start_matches('-'));
        }
        dir
    }

    /// The cache root of the Ensembl layout, `cache_root/`.
    pub fn cache_root(&self) -> PathBuf {
        self.dir.join("cache_root")
    }

    /// Whether the run relies on the binary finding the reference inside the
    /// cache directory (`fasta_in_cache: true`): the manifest's `fasta` is then
    /// staged into a copy of the cache and no `--fasta` is passed.
    pub fn fasta_in_cache(&self) -> bool {
        match self.manifest.get("fasta_in_cache") {
            None => false,
            Some(serde_json::Value::Bool(flag)) => {
                assert!(
                    !*flag || self.manifest["fasta"].is_string(),
                    "{}/{}: fasta_in_cache needs a fasta",
                    self.release,
                    self.name
                );
                *flag
            }
            Some(other) => panic!(
                "{}/{}: fasta_in_cache is {other}, not a boolean",
                self.release, self.name
            ),
        }
    }

    /// The arguments naming the committed cache: `--json_cache <dir>`, or under
    /// the Ensembl layout `--cache --dir_cache <root> --cache_version <release>`,
    /// so the run resolves `<root>/<species>/<release>_<assembly>` itself.
    pub fn cache_args(&self) -> Vec<String> {
        if self.ensembl_cache_layout() {
            vec![
                "--cache".to_string(),
                "--dir_cache".to_string(),
                self.cache_root().display().to_string(),
                "--cache_version".to_string(),
                self.release.clone(),
            ]
        } else {
            vec![
                "--json_cache".to_string(),
                self.cache_dir().display().to_string(),
            ]
        }
    }

    /// The arguments naming the cache a run under `out_dir` reads: the staged
    /// copy holding the reference for a `fasta_in_cache` corpus, else
    /// [`Self::cache_args`].
    pub fn run_cache_args(&self, out_dir: &Path) -> Vec<String> {
        if self.fasta_in_cache() {
            vec![
                "--json_cache".to_string(),
                staged_cache(self, out_dir).display().to_string(),
            ]
        } else {
            self.cache_args()
        }
    }

    /// The cache directory a run under `out_dir` names: the staged copy holding
    /// the reference for a `fasta_in_cache` corpus, else the committed directory.
    pub fn run_cache_dir(&self, out_dir: &Path) -> PathBuf {
        if self.fasta_in_cache() {
            staged_cache(self, out_dir)
        } else {
            self.cache_dir()
        }
    }

    /// The mask classes the manifest's `masks` list declares, each a known class.
    pub fn masks(&self) -> Vec<String> {
        mask_classes(
            &format!("{}/{}", self.release, self.name),
            &self.manifest["masks"],
        )
    }
}

/// The mask class of a corpus whose cache holds no regulatory shards: the
/// reference's `RegulatoryFeature` and `MotifFeature` rows are set aside.
pub const REGULATORY_FEATURE_ROWS: &str = "regulatory_feature_rows";

/// The mask class of a run whose list-valued columns carry the same items in
/// another order. Ensembl VEP prints a transcript's protein features in the
/// order of its cache's array (`Variation/BaseTranscriptVariation.pm:634-643`),
/// which the core API fills by analysis in Perl hash order with the database's
/// rows inside each analysis, while vep-cache-builder writes them by feature
/// identifier, so `DOMAINS` agrees as a set and not in order; `FREQS` lists the
/// frequencies that filtered a record in hash order (`OutputFactory.pm:1152-1153`),
/// as does the `ZYG` list of `--individual_zyg` (`OutputFactory.pm:983`), which
/// vep-rs prints in the order the run names the samples; the UTRAnnotator
/// plugin's `5UTR_annotation` writes its `key=value` pairs and its groups in
/// hash order (UTRAnnotator.pm:282-288), which vep-rs fixes, as the
/// NearestExonJB plugin lists its exon items (NearestExonJB.pm:267-271); and
/// `CLIN_SIG` under `--clin_sig_allele 1` joins the different term lists of a
/// row's known variants in hash order (`OutputFactory.pm:1120`, `1143-1144`),
/// which vep-rs writes in the order of the known variants. Under the class the
/// items of the four item columns, the lists of `CLIN_SIG` and the pair groups
/// of `5UTR_annotation` in the default, tab and VCF formats and the JSON
/// `domains` and `zyg` arrays and `NearestExonJB` string compare as a multiset
/// ([`list_columns_as_multisets`]); every other field compares exactly.
pub const ORDER_ONLY: &str = "order_only";

/// The mask class of a corpus whose flags name a plugin module of its own that
/// Ensembl VEP loads from the plugins directory and vep-rs, which runs no such
/// module, cannot. The reference's run-time warnings about that plugin,
/// `Plugin '<Name>' did not return a hashref, output ignored!` and `Plugin
/// '<Name>' went wrong: <reason>` (OutputFactory.pm:2244 and :2253), are set
/// aside, and the one line vep-rs writes about it, `Failed to load plugin
/// <Name>: ...`, is required in their place and set aside too; every other
/// warning compares as usual ([`compare_warnings`]).
pub const UNLOADED_MODULE_PLUGIN: &str = "unloaded_module_plugin";

/// The mask class of a run under `--phyloP`: Ensembl VEP's option set names
/// each source `phlyoP<N>way` while the file it reads is `phyloP<N>way`
/// (`Config.pm:555`), so the misspelt key heads the reference's header line,
/// tab column, `Extra` item, CSQ field and JSON `custom_annotations` member
/// where vep-rs writes `phyloP<N>way`. Under the class the reference text is
/// read with the key spelt as vep-rs spells it ([`masked_reference_text`]).
pub const CONSERVATION_SHORT_NAME_SPELLING: &str = "conservation_short_name_spelling";

/// The mask class of a run whose `--custom` source asks for the `PC` field:
/// Ensembl VEP describes it (`AnnotationSourceAdaptor.pm:395-397`) as the
/// share of the input variant covered under `reciprocal=0` and as a minimum
/// reciprocal overlap under `reciprocal=1`, while `File.pm:556-564` computes
/// the share of the annotation record the variant covers without `reciprocal`
/// and the share of the variant the record covers with it. vep-rs's header
/// names the share computed; under the class the reference's
/// `## <key>_PC : <text> from <file>` and `##INFO=<ID=<key>_PC,...>` lines
/// are read with that text ([`masked_reference_text`]).
pub const PC_DESCRIPTION_TEXT: &str = "pc_description_text";

/// The classes of the behaviours vep-rs corrects where Ensembl VEP 116.3 is
/// wrong, each a page or section of `docs/intended-divergences.md`. A corpus
/// documents the rows of such a class in its `divergences`, `vep_rs_only_tuples`
/// or `field_divergences` under `expected_divergence`, which the harness applies
/// whatever the class, and names the class in `masks` so the comparison states
/// what it sets aside; the format comparator's built-in classes of the same
/// names recognise the shapes without a manifest.
const CORRECTED_BEHAVIOUR_CLASSES: &[&str] = &[
    "start_cooccurrence_swap",
    "splice_family_swap",
    "start_codon_removed_without_five_prime_utr",
    "start_codon_kept_without_five_prime_utr",
    "gxf_reverse_strand_utr_sequence",
    "conservation_max_text_order",
    "gwas_haplotype_cross_product",
    "alphamissense_two_digit_transcript_version",
    "cadd_structural_allele_match",
    "splicevault_sibling_allele_rows",
    "splicevault_json_later_rows",
    "transcript_annotator_sibling_allele_rows",
    "utrannotator_gained_start_window",
    "loeuf_zero_minimum",
    "maxentscan_overlapping_windows",
    "loftee_gerp_walk_without_conservation_database",
    "dosage_sensitivity_regulatory_row_warnings",
    "hgvs_shift_shared_across_alleles",
    "hgvs_shift_sequence_edge_window",
    "refseq_edited_model_hgvsc_offset",
    "refseq_edited_model_intron_span_coding_reference",
    "cnv_tr_symbolic_loss_named_a_gain",
    "organelle_peptide_paths_translation_table",
    "hgvs_shifted_indel_written_span_gate",
    "hgvs_shifted_deletion_clamped_span",
    "hgvsp_shifted_indel_written_span_start_test",
    "hgvsp_start_codon_test_past_the_shifted_allele",
    "hgvs_clipped_reference_against_whole_alternate",
    "hgvs_added_bases_left_at_the_span_end",
    "hgvs_insertion_left_where_the_shift_stops",
    "hgvs_added_bases_written_outside_the_transcript",
    "total_length_json_range_end",
    "allele_number_repeated_or_reference_equal_alt",
];

/// The mask classes that change how the harness compares; a manifest may also
/// name a class of [`CORRECTED_BEHAVIOUR_CLASSES`].
const MASK_CLASSES: &[&str] = &[
    REGULATORY_FEATURE_ROWS,
    ORDER_ONLY,
    UNLOADED_MODULE_PLUGIN,
    CONSERVATION_SHORT_NAME_SPELLING,
    PC_DESCRIPTION_TEXT,
    FILE_CACHE_WRITES_NO_SOURCE_FIELD,
];

/// The misspelt key's parts: the prefix Ensembl VEP's set writes and vep-rs's
/// spelling of it; the width and `way` follow each.
const MISSPELT_CONSERVATION_PREFIX: &str = "phlyoP";
const CONSERVATION_PREFIX: &str = "phyloP";

/// Ensembl VEP's two `PC` descriptions, each paired with vep-rs's text for the
/// share that setting of `reciprocal` computes.
const PC_DESCRIPTIONS: [(&str, &str); 2] = [
    (
        "Percentage of input variant covered by reference variant",
        "Percentage of the annotation record covered by the input variant",
    ),
    (
        "Percentage of minimum reciprocal overlap between input variant and reference variant",
        "Percentage of the input variant covered by the annotation record",
    ),
];

/// The reference output of one format as the run's text mask classes read it,
/// applied before the header comparison and the parse. Under
/// [`CONSERVATION_SHORT_NAME_SPELLING`] every `phlyoP<N>way` (a digit run
/// between the prefix and `way`) becomes `phyloP<N>way`, in the header lines
/// and the rows alike; under [`PC_DESCRIPTION_TEXT`] every header line
/// describing a `<key>_PC` field with one of Ensembl VEP's two texts carries
/// vep-rs's text for the share the field reports. Without either class the
/// text is returned as it is.
pub fn masked_reference_text(masks: &[String], text: String) -> String {
    let mut text = text;
    if masks_name(masks, CONSERVATION_SHORT_NAME_SPELLING) {
        text = respelt_conservation_keys(&text);
    }
    if masks_name(masks, PC_DESCRIPTION_TEXT) {
        text = with_pc_descriptions(&text);
    }
    text
}

/// `text` with every `phlyoP<N>way` key spelt `phyloP<N>way`; a `phlyoP` not
/// followed by a digit run and `way` stays as written.
fn respelt_conservation_keys(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find(MISSPELT_CONSERVATION_PREFIX) {
        let after = &rest[at + MISSPELT_CONSERVATION_PREFIX.len()..];
        let digits = after.len() - after.trim_start_matches(|c: char| c.is_ascii_digit()).len();
        let is_key = digits > 0 && after[digits..].starts_with("way");
        out.push_str(&rest[..at]);
        out.push_str(if is_key {
            CONSERVATION_PREFIX
        } else {
            MISSPELT_CONSERVATION_PREFIX
        });
        rest = after;
    }
    out.push_str(rest);
    out
}

/// `text` with every `## <key>_PC : ...` and `##INFO=<ID=<key>_PC,...>` header
/// line carrying vep-rs's `PC` description in place of Ensembl VEP's; the key a
/// header line describes is read as [`without_transcript_source`] reads it.
fn with_pc_descriptions(text: &str) -> String {
    fn described_key(line: &str) -> Option<&str> {
        if let Some(rest) = line.strip_prefix("## ") {
            rest.split_once(" : ").map(|(key, _)| key)
        } else {
            line.strip_prefix("##INFO=<ID=")
                .and_then(|rest| rest.split(',').next())
        }
    }
    let mut out: Vec<String> = Vec::new();
    for line in text.split_inclusive('\n') {
        let is_pc_line = described_key(line.trim_end()).is_some_and(|key| key.ends_with("_PC"));
        let replaced = PC_DESCRIPTIONS
            .iter()
            .find(|(perl_text, _)| is_pc_line && line.contains(perl_text))
            .map(|(perl_text, vep_rs_text)| line.replacen(perl_text, vep_rs_text, 1));
        out.push(replaced.unwrap_or_else(|| line.to_string()));
    }
    out.concat()
}

/// The classes a `masks` list names, each one of [`MASK_CLASSES`] or
/// [`CORRECTED_BEHAVIOUR_CLASSES`].
fn mask_classes(label: &str, masks: &serde_json::Value) -> Vec<String> {
    let masks = strings(masks);
    for class in &masks {
        assert!(
            MASK_CLASSES.contains(&class.as_str())
                || CORRECTED_BEHAVIOUR_CLASSES.contains(&class.as_str()),
            "{label}: unknown mask class {class}"
        );
    }
    masks
}

/// Whether a `masks` list names `class`.
fn masks_name(masks: &[String], class: &str) -> bool {
    masks.iter().any(|m| m == class)
}

/// The `Feature_type` values of the rows [`REGULATORY_FEATURE_ROWS`] drops; in
/// the JSON format the entries of the `regulatory_feature_consequences` and
/// `motif_feature_consequences` arrays.
const REGULATORY_FEATURE_TYPES: [&str; 2] = ["RegulatoryFeature", "MotifFeature"];
/// The warnings-file phrase vep-rs writes when the cache holds no regulatory
/// data; [`REGULATORY_FEATURE_ROWS`] requires the line.
const REGULATORY_NOTICE: &str = "this cache holds no regulatory data";
/// The warnings-file phrase vep-rs writes when a run asks for known variants
/// from a cache that holds no variation data. Ensembl VEP attaches no variation
/// source to such a cache and writes nothing about it (`CacheDir.pm:179`), so
/// [`compare_warnings`] sets the line aside on a corpus whose cache holds no
/// `variations` directory.
const VARIATION_NOTICE: &str = "[cache] no variation data under";

/// Applies a run's mask classes to the reference entries of one format.
/// Under [`REGULATORY_FEATURE_ROWS`] the default and tab rows and the VCF CSQ
/// entries whose `Feature_type` is regulatory go; a JSON record that carried a
/// regulatory array loses it and has its `most_severe_consequence` recomputed
/// from the consequence entries that remain, as vep-rs computes it from the rows
/// it writes, so a record Ensembl VEP annotates against a regulatory feature and
/// no transcript compares as its intergenic row; every other record keeps the
/// reference's value. [`ORDER_ONLY`] rewrites both sides, so it is applied by
/// [`mask_entries`] through [`list_columns_as_multisets`].
pub fn apply_masks(masks: &[String], format: &str, expected: &mut Vec<Entry>) {
    if !masks_name(masks, REGULATORY_FEATURE_ROWS) {
        return;
    }
    if format != "json" {
        expected.retain(|e| {
            !e.fields
                .get("Feature_type")
                .is_some_and(|t| REGULATORY_FEATURE_TYPES.contains(&t.as_str()))
        });
        return;
    }
    let record_of = |key: &str| key.split('|').next().unwrap_or_default().to_string();
    let regulatory_entry = |e: &Entry| e.key.ends_with("|regulatory") || e.key.ends_with("|motif");
    let mut remaining: HashMap<String, Vec<Consequence>> = HashMap::new();
    for entry in expected
        .iter()
        .filter(|e| !e.key.ends_with("|record") && !regulatory_entry(e))
    {
        let terms = remaining.entry(record_of(&entry.key)).or_default();
        for term in entry
            .consequence
            .split([',', '&'])
            .filter(|t| !t.is_empty())
        {
            terms.push(
                Consequence::from_so_term(term).unwrap_or_else(|| {
                    panic!("{}: {term} is not an SO consequence term", entry.key)
                }),
            );
        }
    }
    let dropped_records: HashSet<String> = expected
        .iter()
        .filter(|e| regulatory_entry(e))
        .map(|e| record_of(&e.key))
        .collect();
    expected.retain(|e| !regulatory_entry(e));
    for entry in expected.iter_mut().filter(|e| e.key.ends_with("|record")) {
        let dropped = dropped_records.contains(&record_of(&entry.key));
        // A record that lost no row keeps Ensembl VEP's own value under comparison.
        if dropped && entry.fields.contains_key("most_severe_consequence") {
            let terms = remaining.remove(&record_of(&entry.key)).unwrap_or_default();
            let most_severe = Consequence::most_severe(&terms)
                .map(|c| c.so_term().to_string())
                .unwrap_or_else(|| "?".to_string());
            entry
                .fields
                .insert("most_severe_consequence".to_string(), most_severe);
        }
    }
}

/// The columns [`ORDER_ONLY`] compares as multisets in the text formats; the
/// JSON format carries `DOMAINS` as a consequence's `domains` array, `ZYG` as
/// the record's `zyg` array and `NearestExonJB` as the consequence's
/// `,`-joined `nearestexonjb` string, and has no joined `CLIN_SIG`: each known
/// variant carries its own terms there.
const ORDER_ONLY_COLUMNS: [&str; 5] = ["DOMAINS", "FREQS", "ZYG", "NearestExonJB", "CLIN_SIG"];

/// The separator of the term lists of `CLIN_SIG` in the text format `format`:
/// the `;` Ensembl VEP joins them with (`OutputFactory.pm:1144`), which the
/// default and VCF formats write as `%3B` (`OutputFactory/VEP_output.pm:118`,
/// `OutputFactory/VCF.pm:402`). The lists are compared whole, so `CLIN_SIG`
/// is split at this separator alone.
fn clin_sig_list_separator(format: &str) -> &'static str {
    if format == "tab" {
        ";"
    } else {
        "%3B"
    }
}

/// The column [`ORDER_ONLY`] compares as a multiset of `key=value` groups: the
/// UTRAnnotator plugin writes one `:`-joined group of pairs per affected
/// upstream ORF, the pairs and the groups in Ensembl VEP's hash order
/// (UTRAnnotator.pm:282-288 at VEP_plugins release/116) and in vep-rs's fixed
/// order; the JSON format nests the groups as objects, which compare by key.
const ORDER_ONLY_PAIR_GROUP_COLUMN: &str = "5UTR_annotation";

/// Rewrites every [`ORDER_ONLY_COLUMNS`] value of one side into its items in
/// sorted order: the default and tab formats' `,`-joined items
/// (`Utils.pm:185-198`), the VCF's `&`-joined items (`OutputFactory/VCF.pm:391`),
/// the lists of `CLIN_SIG` (`clin_sig_list_separator`),
/// the JSON `domains` array's objects by their compact text, the JSON `zyg`
/// array's strings and the JSON `nearestexonjb` string's `,`-joined items; the
/// [`ORDER_ONLY_PAIR_GROUP_COLUMN`]'s groups each into
/// their sorted pairs, then into sorted groups. Applied to both sides under
/// [`ORDER_ONLY`], it makes two lists of the same items in any order compare
/// equal while lists differing in an item, or in how often one occurs, still
/// differ. The default format's raw `Extra` text is rewritten with each field;
/// an absent or `-` value stays as written.
pub fn list_columns_as_multisets(entries: &mut [Entry], format: &str) {
    if format == "json" {
        for entry in entries {
            // `domains` is an array of objects, each kept as its compact text by
            // `json_scalar`, so a JSON parse splits them where `,` would not;
            // `zyg` is an array of strings, joined by `,`.
            if let Some(value) = entry.fields.get("domains").cloned() {
                if let Ok(serde_json::Value::Array(objects)) =
                    serde_json::from_str::<serde_json::Value>(&format!("[{value}]"))
                {
                    let sorted = sorted_items(objects.iter().map(json_scalar), ",");
                    entry.fields.insert("domains".to_string(), sorted);
                }
            }
            for field in ["zyg", "nearestexonjb"] {
                if let Some(value) = entry.fields.get(field).cloned() {
                    let sorted = sorted_items(value.split(','), ",");
                    entry.fields.insert(field.to_string(), sorted);
                }
            }
        }
        return;
    }
    let separator = if format == "vcf" { "&" } else { "," };
    for entry in entries {
        for column in ORDER_ONLY_COLUMNS
            .iter()
            .copied()
            .chain([ORDER_ONLY_PAIR_GROUP_COLUMN])
        {
            let Some(value) = entry.fields.get(column).cloned() else {
                continue;
            };
            if value.is_empty() || value == "-" {
                continue;
            }
            let sorted = if column == ORDER_ONLY_PAIR_GROUP_COLUMN {
                sorted_items(
                    value
                        .split(separator)
                        .map(|group| sorted_items(group.split(':'), ":")),
                    separator,
                )
            } else if column == "CLIN_SIG" {
                let lists = clin_sig_list_separator(format);
                sorted_items(value.split(lists), lists)
            } else {
                sorted_items(value.split(separator), separator)
            };
            if let Some(extra) = entry.fields.get_mut("Extra") {
                *extra = extra.replace(&format!("{column}={value}"), &format!("{column}={sorted}"));
            }
            entry.fields.insert(column.to_string(), sorted);
        }
    }
}

/// `items` in sorted order, joined by `separator`.
fn sorted_items<T: Borrow<str> + Ord>(items: impl Iterator<Item = T>, separator: &str) -> String {
    let mut items: Vec<T> = items.collect();
    items.sort_unstable();
    items.join(separator)
}

/// Applies a run's mask classes to the parsed entries of one format before
/// they are compared: the reference loses what [`apply_masks`] drops, and under
/// [`ORDER_ONLY`] both sides have their list columns rewritten by
/// [`list_columns_as_multisets`].
fn mask_entries(masks: &[String], format: &str, expected: &mut Vec<Entry>, actual: &mut [Entry]) {
    apply_masks(masks, format, expected);
    if masks_name(masks, ORDER_ONLY) {
        list_columns_as_multisets(expected, format);
        list_columns_as_multisets(actual, format);
    }
}

/// A copy of the corpus cache in `out_dir` with the reference and its indexes
/// inside it under their own names (made once; later calls reuse it), for a
/// corpus whose run relies on the cache directory's FASTA detection. The copy's
/// checksum manifest is rewritten over the added files, as a user who places a
/// FASTA in a cache directory lists it.
fn staged_cache(corpus: &Corpus, out_dir: &Path) -> PathBuf {
    let target = out_dir.join(format!("{}-cache", corpus.name));
    if !target.is_dir() {
        // Built under another name and renamed into place, so a copy a failed
        // run left half-written is never taken for a whole one.
        let building = out_dir.join(format!("{}-cache.building", corpus.name));
        let _ = fs::remove_dir_all(&building);
        copy_dir(&corpus.cache_dir(), &building);
        let fasta = committed_fasta(corpus);
        link_or_copy_indexed_fasta(&fasta, &building.join(fasta.file_name().unwrap()));
        if building.join(MANIFEST_FILE).is_file() {
            CacheManifest::write(&building).unwrap();
        }
        fs::rename(&building, &target).unwrap();
    }
    target
}

/// Copies a directory tree.
fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let destination = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &destination);
        } else {
            fs::copy(entry.path(), destination).unwrap();
        }
    }
}

/// The reference FASTA committed beside the manifest, the file its `fasta` names.
fn committed_fasta(corpus: &Corpus) -> PathBuf {
    let name = corpus.manifest["fasta"]
        .as_str()
        .unwrap_or_else(|| panic!("{}/{}: a fasta is needed here", corpus.release, corpus.name));
    corpus.dir.join(name)
}

/// Places `fasta` and the `.fai` and `.gzi` beside it at `target` and its
/// suffixes, as hard links where the file system allows and copies otherwise.
fn link_or_copy_indexed_fasta(fasta: &Path, target: &Path) {
    for suffix in ["", ".fai", ".gzi"] {
        let from = PathBuf::from(format!("{}{suffix}", fasta.display()));
        if !from.is_file() {
            continue;
        }
        let to = PathBuf::from(format!("{}{suffix}", target.display()));
        // A target already there may be a hard link of the source, which
        // `fs::copy` truncates before it reads.
        if fs::hard_link(&from, &to).is_err() && !to.is_file() {
            fs::copy(&from, &to).unwrap();
        }
    }
}

/// Every corpus under `tests/golden/`, sorted by (release, assembly).
pub fn corpora() -> Vec<Corpus> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden");
    let mut out = Vec::new();
    let Ok(releases) = fs::read_dir(&root) else {
        return out;
    };
    for rel in releases.flatten() {
        if !rel.path().is_dir() {
            continue;
        }
        for asm in fs::read_dir(rel.path()).unwrap().flatten() {
            let dir = asm.path();
            if !dir.join("variants.vcf").is_file() {
                continue;
            }
            let manifest: serde_json::Value =
                serde_json::from_str(&fs::read_to_string(dir.join("manifest.json")).unwrap())
                    .unwrap();
            let name = asm.file_name().to_string_lossy().into_owned();
            let assembly = manifest["assembly"]
                .as_str()
                .map(str::to_string)
                .unwrap_or_else(|| name.clone());
            let species = manifest["species"]
                .as_str()
                .map_or_else(|| DEFAULT_SPECIES.to_string(), str::to_string);
            out.push(Corpus {
                release: rel.file_name().to_string_lossy().into_owned(),
                name,
                assembly,
                species,
                dir,
                manifest,
            });
        }
    }
    out.sort_by(|a, b| (&a.release, &a.name).cmp(&(&b.release, &b.name)));
    out
}

/// Reads `<dir>/expected/<name>` or its `.gz`, whichever exists.
pub fn read_expected(dir: &Path, name: &str) -> String {
    read_maybe_gz(&dir.join("expected").join(name))
}

/// Reads `path` or `path.gz`, whichever exists.
fn read_maybe_gz(path: &Path) -> String {
    if path.is_file() {
        return fs::read_to_string(path).unwrap();
    }
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let file = fs::File::open(&gz).unwrap_or_else(|e| panic!("{}: {e}", gz.display()));
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// The data lines of an input file (a `.gz` or `.bgz` is decompressed), in
/// order; a line that is not UTF-8 is read with its invalid bytes replaced, as
/// the binary reads it for its text outputs.
fn data_lines(path: &Path) -> Vec<String> {
    input_text(path)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// An input file's text, decompressed when its name says it is gzipped.
fn input_text(path: &Path) -> String {
    let bytes = if path.extension().is_some_and(|e| e == "gz" || e == "bgz") {
        let mut decompressed = Vec::new();
        flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap())
            .read_to_end(&mut decompressed)
            .unwrap();
        decompressed
    } else {
        fs::read(path).unwrap()
    };
    String::from_utf8_lossy(&bytes).into_owned()
}

/// Whether the input file is a VCF: one carries a `##fileformat=VCF` or a
/// `#CHROM` header line, which no other input format does.
fn input_is_vcf(path: &Path) -> bool {
    input_text(path)
        .lines()
        .take_while(|l| l.starts_with('#'))
        .any(|l| l.starts_with("##fileformat=VCF") || l.starts_with("#CHROM"))
}

/// The output formats a run compares when its entry names none.
const ALL_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// One run of a corpus: the manifest's own run (named `default`) or one of its
/// `flag_sets` entries, each with its own flags, input, formats and expected
/// outputs.
pub struct Run<'a> {
    pub corpus: &'a Corpus,
    pub name: String,
    pub flags: Vec<String>,
    /// The input file, relative to the corpus directory.
    pub input: String,
    /// Whether the input is piped on standard input, with `-i STDIN`.
    pub stdin: bool,
    pub formats: Vec<String>,
    /// The directory of the expected outputs, relative to the corpus directory.
    pub expected_dir: String,
    /// The reference's rows within a record come out in Ensembl VEP's hash
    /// order (`--per_gene` emits its gene buckets and `--pick_allele` its allele
    /// buckets from a Perl hash), so the run compares the entries and their
    /// fields, not the order of the entries.
    pub unordered: bool,
    /// How the reference was generated: `offline` (Ensembl VEP on its cache,
    /// the default), `database` (Ensembl VEP connected to the Ensembl
    /// database, for a column it writes only with one) or `gtf` or `gff`
    /// (Ensembl VEP reading its transcripts from the file given to `--gtf` or
    /// `--gff`, [`GXF_PERL_MODES`]). A database reference's preamble describes
    /// a connection rather than a cache, so its default and tab headers compare
    /// by their column descriptions and column line alone; a GTF or GFF
    /// reference's names no cache either and compares the same way, after its
    /// transcript source's `SOURCE` column and header lines are set aside
    /// ([`without_transcript_source`]).
    pub perl_mode: String,
    /// The mask classes applied to the run: the flag set's own `masks`, else
    /// the manifest's ([`Corpus::masks`]).
    pub masks: Vec<String>,
    /// The manifest object whose divergence lists apply to this run.
    documented: &'a serde_json::Value,
    /// The run's own manifest entry: the manifest for its own run, the flag
    /// set's object for a flag set.
    entry: &'a serde_json::Value,
}

impl Run<'_> {
    /// `<release>/<corpus>` for the manifest's own run, `<release>/<corpus> [<set>]` for a flag set.
    pub fn label(&self) -> String {
        if self.name == "default" {
            format!("{}/{}", self.corpus.release, self.corpus.name)
        } else {
            format!(
                "{}/{} [{}]",
                self.corpus.release, self.corpus.name, self.name
            )
        }
    }

    pub fn compares(&self, format: &str) -> bool {
        self.formats.iter().any(|f| f == format)
    }

    /// Whether the run's flags carry `flag`, as a word or as `flag=value`.
    pub fn has_flag(&self, flag: &str) -> bool {
        self.flags
            .iter()
            .any(|f| f == flag || f.starts_with(&format!("{flag}=")))
    }

    /// The value of a `flag value` or `flag=value` pair in the run's flags.
    fn flag_value(&self, flag: &str) -> Option<&str> {
        let prefix = format!("{flag}=");
        self.flags.iter().enumerate().find_map(|(i, f)| {
            if let Some(value) = f.strip_prefix(&prefix) {
                Some(value)
            } else if f == flag {
                self.flags.get(i + 1).map(String::as_str)
            } else {
                None
            }
        })
    }

    /// The INFO field the run's VCF output writes its annotation under:
    /// `--vcf_info_field`, else `CSQ`.
    pub fn vcf_info_field(&self) -> &str {
        self.flag_value("--vcf_info_field").unwrap_or("CSQ")
    }

    /// Whether the run reads a VCF: `--format vcf` says so, another named
    /// format (`ensembl`, `region`, `hgvs`) says not, and without a `--format`
    /// word the input file decides, as Ensembl VEP detects the format itself
    /// ([`input_is_vcf`]). A VCF output carries a copy of a VCF input's line
    /// and of no other format's.
    pub fn reads_vcf(&self) -> bool {
        match self.flag_value("--format") {
            Some(format) => format == "vcf",
            None => input_is_vcf(&self.corpus.dir.join(&self.input)),
        }
    }

    /// The expected output `name` (or its `.gz`) in the run's expected directory.
    pub fn expected(&self, name: &str) -> String {
        read_maybe_gz(&self.corpus.dir.join(&self.expected_dir).join(name))
    }

    /// The run's input data lines, in order.
    pub fn input_records(&self) -> Vec<String> {
        data_lines(&self.corpus.dir.join(&self.input))
    }

    /// The run's documented keys, each `divergences` item read for the run's
    /// setting: with or without a reference sequence ([`Self::has_reference_sequence`]).
    pub fn documented(&self, out_dir: &Path) -> Documented {
        Documented::for_run(self.documented, self.has_reference_sequence(out_dir))
    }

    /// Whether the binary has a reference sequence in this run, by the words
    /// the run gives it ([`base_args`]) and the cache directory it names.
    pub fn has_reference_sequence(&self, out_dir: &Path) -> bool {
        has_reference_sequence(
            &base_args(self, out_dir),
            &self.corpus.run_cache_dir(out_dir),
            self.corpus,
        )
    }

    /// Whether the run's `NEAREST` values compare as transcription start
    /// coordinates: under `--nearest transcript` Ensembl VEP names one transcript
    /// per winning coordinate, whichever its interval tree held first, and vep-rs
    /// names every transcript at the coordinate, so the two agree on the
    /// coordinates and not on the ids ([`nearest_as_coordinates`]).
    pub fn nearest_by_coordinate(&self) -> bool {
        self.flag_value("--nearest") == Some("transcript")
    }

    /// Whether the run masks the reference's regulatory rows ([`REGULATORY_FEATURE_ROWS`]).
    pub fn masks_regulatory_rows(&self) -> bool {
        masks_name(&self.masks, REGULATORY_FEATURE_ROWS)
    }

    /// Whether the run's reference ran a plugin module vep-rs cannot load
    /// ([`UNLOADED_MODULE_PLUGIN`]).
    pub fn masks_unloaded_module_plugin(&self) -> bool {
        masks_name(&self.masks, UNLOADED_MODULE_PLUGIN)
    }

    /// Whether the run compares the columns of [`ORDER_ONLY_COLUMNS`] and of
    /// [`ORDER_ONLY_PAIR_GROUP_COLUMN`] as multisets ([`ORDER_ONLY`]).
    pub fn order_only(&self) -> bool {
        masks_name(&self.masks, ORDER_ONLY)
    }

    /// Whether the reference read its transcripts from a GTF or GFF file
    /// ([`GXF_PERL_MODES`]).
    pub fn transcript_source_reference(&self) -> bool {
        GXF_PERL_MODES.contains(&self.perl_mode.as_str())
    }

    /// Whether the reference's preamble names no cache (a database or a GTF or
    /// GFF reference), so the headers compare by their column descriptions.
    pub fn headers_by_descriptions(&self) -> bool {
        self.perl_mode == "database" || self.transcript_source_reference()
    }
}

/// The `perl_mode` values of a reference Ensembl VEP wrote reading its
/// transcripts from a GTF or GFF annotation source (`--gtf <file>`, `--gff
/// <file>`: `Config.pm:532-541` turns each into a `--custom` source of that
/// format) instead of a cache, while vep-rs annotates against the cache
/// `vep-cache-builder` built from the same file. Such a reference stamps every
/// transcript it reads with the source's short name and writes it as `SOURCE`,
/// with the source's own column and the header lines describing both, which a
/// cache built from the file never carries: the documented intended divergence
/// `file_cache_writes_no_source_field`
/// (`docs/divergences/file-cache-writes-no-source-field.md`).
/// [`without_transcript_source`] sets them aside before the comparison.
pub const GXF_PERL_MODES: [&str; 2] = ["gtf", "gff"];

/// Every `perl_mode` a manifest or flag set may name.
const PERL_MODES: [&str; 4] = ["offline", "database", "gtf", "gff"];

/// The mask class of a run whose reference read its transcripts from a GTF or
/// GFF ([`GXF_PERL_MODES`]): the `SOURCE` column, the source's own column and
/// the header lines describing both, which [`without_transcript_source`] sets
/// aside. A corpus with such a reference names the class in `masks`, so the
/// comparison states what it sets aside ([`compare_run`]).
pub const FILE_CACHE_WRITES_NO_SOURCE_FIELD: &str = "file_cache_writes_no_source_field";

/// The column the transcript source stamps on every row it annotates.
const SOURCE_COLUMN: &str = "SOURCE";

/// A GTF or GFF reference output of one format without its transcript source,
/// the class `file_cache_writes_no_source_field`
/// (`docs/divergences/file-cache-writes-no-source-field.md`): the `SOURCE`
/// column (the default format's `SOURCE=<file>` Extra item, the tab column, the
/// CSQ field and the JSON consequence's `source` key), the source's own
/// annotation column, named by its short name and empty on every row since a
/// transcript source annotates no record (the tab column and the CSQ field),
/// and the header lines describing both (`## SOURCE : Source of transcript`,
/// the `## <file> : <path>` line naming the source file, the tab column
/// headers, the two names in the CSQ format list and the VCF's
/// `##INFO=<ID=<file>,...>` line). The source's line is recognised by its text
/// ending in the key it describes, as Ensembl VEP writes a custom source's
/// header (`OutputFactory.pm:2444-2479`, `## <short name> : <path>`), which no
/// other header line does. In the VCF the fields dropped are those the last
/// `##INFO=<ID=<info_field>,...>` line carrying a `Format: ` clause names, and
/// the entry rewritten is a line's last one under the field, as
/// [`parse_vcf_info`] reads them; a header line of the field without the clause
/// and an earlier entry (one `--keep_csq` kept) stay as written.
pub fn without_transcript_source(format: &str, text: &str, info_field: &str) -> String {
    fn describes_itself(key: &str, description: &str) -> bool {
        description
            .trim_end_matches(['"', '>'])
            .rsplit('/')
            .next()
            .is_some_and(|base| base == key)
    }
    /// The entries of an INFO annotation value (`,`-joined, their fields
    /// `|`-joined) without the fields at `dropped`.
    fn without_fields(entries: &str, dropped: &[usize]) -> String {
        entries
            .split(',')
            .map(|entry| {
                let fields: Vec<&str> = entry.split('|').collect();
                without_positions(&fields, dropped).join("|")
            })
            .collect::<Vec<_>>()
            .join(",")
    }
    let mut dropped: BTreeSet<String> = BTreeSet::from([SOURCE_COLUMN.to_string()]);
    for line in text.lines().filter(|l| l.starts_with("##")) {
        let (key, description) = if let Some(rest) = line.strip_prefix("## ") {
            match rest.split_once(" : ") {
                Some(pair) => pair,
                None => continue,
            }
        } else if let Some(rest) = line.strip_prefix("##INFO=<ID=") {
            let key = rest.split(',').next().unwrap_or_default();
            let description = rest.rsplit("Description=\"").next().unwrap_or_default();
            (key, description)
        } else {
            continue;
        };
        if describes_itself(key, description) {
            dropped.insert(key.to_string());
        }
    }
    let is_dropped = |name: &str| dropped.contains(name);
    let without_items = |extra: &str| -> String {
        let kept: Vec<&str> = extra
            .split(';')
            .filter(|item| !item.split_once('=').is_some_and(|(key, _)| is_dropped(key)))
            .collect();
        if kept.is_empty() {
            "-".to_string()
        } else {
            kept.join(";")
        }
    };
    let mut out: Vec<String> = Vec::new();
    match format {
        "default" | "tab" => {
            let mut dropped_columns: Vec<usize> = Vec::new();
            for line in text.lines() {
                if let Some(rest) = line.strip_prefix("## ") {
                    if rest
                        .split_once(" : ")
                        .is_some_and(|(key, _)| is_dropped(key))
                    {
                        continue;
                    }
                    out.push(line.to_string());
                } else if line.starts_with('#') {
                    let columns: Vec<&str> = line.split('\t').collect();
                    dropped_columns = (0..columns.len())
                        .filter(|&i| is_dropped(columns[i]))
                        .collect();
                    out.push(without_positions(&columns, &dropped_columns).join("\t"));
                } else if dropped_columns.is_empty() {
                    // The default format's column line names no dropped column: its
                    // SOURCE is an item of the Extra column.
                    let mut columns: Vec<String> = line.split('\t').map(str::to_string).collect();
                    if let Some(extra) = columns.last_mut() {
                        *extra = without_items(extra);
                    }
                    out.push(columns.join("\t"));
                } else {
                    let columns: Vec<&str> = line.split('\t').collect();
                    out.push(without_positions(&columns, &dropped_columns).join("\t"));
                }
            }
        }
        "vcf" => {
            let mut dropped_fields: Vec<usize> = Vec::new();
            let info_prefix = format!("##INFO=<ID={info_field},");
            let entry_prefix = format!("{info_field}=");
            for line in text.lines() {
                if let Some((head, fields)) = line
                    .strip_prefix(&info_prefix)
                    .and_then(|rest| rest.split_once("Format: "))
                {
                    let names: Vec<&str> = fields.trim_end_matches("\">").split('|').collect();
                    dropped_fields = (0..names.len()).filter(|&i| is_dropped(names[i])).collect();
                    let kept = without_positions(&names, &dropped_fields).join("|");
                    out.push(format!("{info_prefix}{head}Format: {kept}\">"));
                } else if let Some(rest) = line.strip_prefix("##INFO=<ID=") {
                    if is_dropped(rest.split(',').next().unwrap_or_default()) {
                        continue;
                    }
                    out.push(line.to_string());
                } else if line.starts_with('#') || dropped_fields.is_empty() {
                    out.push(line.to_string());
                } else {
                    let mut columns: Vec<String> = line.split('\t').map(str::to_string).collect();
                    if let Some(info) = columns.get_mut(7) {
                        let mut items: Vec<String> = info.split(';').map(str::to_string).collect();
                        if let Some(at) = items
                            .iter()
                            .rposition(|item| item.starts_with(&entry_prefix))
                        {
                            let entries = &items[at][entry_prefix.len()..];
                            items[at] = format!(
                                "{entry_prefix}{}",
                                without_fields(entries, &dropped_fields)
                            );
                            *info = items.join(";");
                        }
                    }
                    out.push(columns.join("\t"));
                }
            }
        }
        _ => {
            for line in text.lines() {
                if line.trim().is_empty() {
                    continue;
                }
                let mut object: serde_json::Value = serde_json::from_str(line)
                    .unwrap_or_else(|e| panic!("a JSON reference line does not parse: {e}"));
                if let Some(consequences) = object
                    .get_mut("transcript_consequences")
                    .and_then(|v| v.as_array_mut())
                {
                    for consequence in consequences.iter_mut() {
                        if let Some(map) = consequence.as_object_mut() {
                            map.remove("source");
                        }
                    }
                }
                out.push(object.to_string());
            }
        }
    }
    let mut joined = out.join("\n");
    if text.ends_with('\n') {
        joined.push('\n');
    }
    joined
}

/// `items` without the elements at `positions`.
fn without_positions<T: Clone>(items: &[T], positions: &[usize]) -> Vec<T> {
    items
        .iter()
        .enumerate()
        .filter(|(i, _)| !positions.contains(i))
        .map(|(_, item)| item.clone())
        .collect()
}

/// Rewrites every `NEAREST` value (the default, tab and VCF field, the JSON
/// record's `nearest`) from transcript ids to the coordinates of their start
/// sites, each coordinate once in order; an id the cache's protein-coding
/// transcripts lack stays as written, so it still mismatches. The default
/// format's raw `Extra` text is rewritten with it.
pub fn nearest_as_coordinates(entries: &mut [Entry], sites: &HashMap<String, u64>) {
    let resolve = |value: &str| -> String {
        let mut out: Vec<String> = Vec::new();
        for id in value.split(['&', ',']).filter(|s| !s.is_empty()) {
            let coordinate = sites
                .get(id)
                .map(u64::to_string)
                .unwrap_or_else(|| id.to_string());
            if out.last() != Some(&coordinate) {
                out.push(coordinate);
            }
        }
        out.join(",")
    };
    for entry in entries {
        for key in ["NEAREST", "nearest"] {
            let Some(value) = entry.fields.get(key).cloned() else {
                continue;
            };
            let resolved = resolve(&value);
            if let Some(extra) = entry.fields.get_mut("Extra") {
                *extra = extra.replace(&format!("NEAREST={value}"), &format!("NEAREST={resolved}"));
            }
            entry.fields.insert(key.to_string(), resolved);
        }
    }
}

/// The transcription start site of every protein-coding transcript of the
/// corpus cache, by stable id: the start on the forward strand, the end on the
/// reverse (`AnnotationType/Transcript.pm:749`).
pub fn protein_coding_start_sites(corpus: &Corpus) -> HashMap<String, u64> {
    fn integer(v: &serde_json::Value) -> Option<u64> {
        v.as_u64().or_else(|| v.as_str()?.parse().ok())
    }
    let mut sites = HashMap::new();
    let transcripts = corpus.cache_dir().join("transcripts");
    for chromosome in fs::read_dir(&transcripts)
        .unwrap_or_else(|e| panic!("{}: {e}", transcripts.display()))
        .flatten()
    {
        for shard in fs::read_dir(chromosome.path()).unwrap().flatten() {
            let path = shard.path();
            let plain = path.to_str().unwrap().trim_end_matches(".gz");
            let text = read_maybe_gz(Path::new(plain));
            let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
            // A shard is a bare array of transcripts or a map of such arrays.
            let records: Vec<&serde_json::Value> = match &parsed {
                serde_json::Value::Array(items) => items.iter().collect(),
                serde_json::Value::Object(map) => map
                    .values()
                    .filter_map(|v| v.as_array())
                    .flatten()
                    .filter(|v| !v.is_null())
                    .collect(),
                other => panic!("{}: {other}", path.display()),
            };
            for t in records {
                if t["biotype"].as_str() != Some("protein_coding") {
                    continue;
                }
                let forward = t["strand"].as_i64() == Some(1) || t["strand"].as_str() == Some("1");
                let site = integer(&t[if forward { "start" } else { "end" }]).unwrap();
                sites.insert(t["stable_id"].as_str().unwrap().to_string(), site);
            }
        }
    }
    sites
}

fn strings(value: &serde_json::Value) -> Vec<String> {
    value
        .as_array()
        .into_iter()
        .flatten()
        .map(|v| v.as_str().expect("manifest lists are strings").to_string())
        .collect()
}

/// Holds a `reference_warnings` list, the manifest's or a flag set's (`label`
/// names which), to the shape the corpus tool holds the key to: strings, none
/// empty or blank, so no run expects a warning with no text.
fn assert_warning_texts(label: &str, list: &serde_json::Value) {
    let texts: Vec<&str> = list
        .as_array()
        .into_iter()
        .flatten()
        .map(|v| {
            v.as_str()
                .unwrap_or_else(|| panic!("{label}: reference_warnings is not of its type"))
        })
        .collect();
    assert!(
        texts.iter().all(|t| !t.trim().is_empty()),
        "{label}: reference_warnings has an empty text"
    );
}

/// The placeholder a flag word or a `runs[].argv` token carries for the corpus directory.
const CORPUS_PLACEHOLDER: &str = "{corpus}";
/// The placeholder for the plugin directory each engine reads: the corpus's
/// `plugin_data` here, Ensembl VEP's plugin modules in the reference run.
const PLUGINS_PLACEHOLDER: &str = "{plugins}";

/// A flag word or `argv` token with `{corpus}` and `{plugins}` replaced.
fn expand_corpus_placeholders(corpus: &Corpus, word: &str) -> String {
    word.replace(CORPUS_PLACEHOLDER, &corpus.dir.display().to_string())
        .replace(
            PLUGINS_PLACEHOLDER,
            &corpus.dir.join("plugin_data").display().to_string(),
        )
}

/// A `flags` list of the manifest or of a flag set, with `{corpus}` and
/// `{plugins}` inside any word replaced ([`expand_corpus_placeholders`]).
fn corpus_flags(corpus: &Corpus, value: &serde_json::Value) -> Vec<String> {
    strings(value)
        .iter()
        .map(|flag| expand_corpus_placeholders(corpus, flag))
        .collect()
}

/// The `perl_mode` of a manifest or flag set object, one of [`PERL_MODES`];
/// `fallback` when the object names none.
fn perl_mode_of(label: &str, object: &serde_json::Value, fallback: &str) -> String {
    let mode = match object.get("perl_mode") {
        None => fallback,
        Some(serde_json::Value::String(mode)) => mode.as_str(),
        Some(_) => panic!("{label}: perl_mode is not of its type"),
    };
    assert!(
        PERL_MODES.contains(&mode),
        "{label}: perl_mode is one of {}, not {mode:?}",
        PERL_MODES.join(", ")
    );
    mode.to_string()
}

/// The manifest's own run. Its `perl_mode` is the manifest's (`offline` when
/// it names none), which a flag set inherits unless it names its own.
pub fn default_run(corpus: &Corpus) -> Run<'_> {
    let label = format!("{}/{}", corpus.release, corpus.name);
    Run {
        corpus,
        name: "default".to_string(),
        flags: corpus_flags(corpus, &corpus.manifest["flags"]),
        input: "variants.vcf".to_string(),
        stdin: false,
        formats: ALL_FORMATS.iter().map(|f| f.to_string()).collect(),
        expected_dir: "expected".to_string(),
        unordered: false,
        perl_mode: perl_mode_of(&label, &corpus.manifest, "offline"),
        masks: corpus.masks(),
        documented: &corpus.manifest,
        entry: &corpus.manifest,
    }
}

/// Whether a manifest path stays inside the corpus directory: non-empty,
/// relative, with no `..` segment.
fn inside_corpus(path: &str) -> bool {
    let parsed = Path::new(path);
    !path.is_empty()
        && !parsed.is_absolute()
        && !parsed
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir))
}

/// Whether a flag set name is one the corpus tool writes: a letter or digit,
/// then letters, digits, `_` and `-`, and not `default`, the manifest's own run.
fn plain_name(name: &str) -> bool {
    let mut chars = name.chars();
    name != "default"
        && chars.next().is_some_and(|c| c.is_ascii_alphanumeric())
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
}

/// The manifest's own run followed by one run per `flag_sets` entry. An entry
/// the corpus tool would refuse (a name that is `default`, repeated or not a
/// plain word, no format to compare, an empty path or one outside the corpus,
/// a stdin run whose flags name no `--format`) fails here, so a hand-written
/// entry cannot certify its `covers` with a run that compares nothing; so does
/// a manifest or an entry whose `masks` is not a list, or whose
/// `reference_warnings` is not a list of texts or holds an empty one
/// ([`assert_warning_texts`]), before any run.
pub fn runs(corpus: &Corpus) -> Vec<Run<'_>> {
    let label = format!("{}/{}", corpus.release, corpus.name);
    for key in ["masks", "reference_warnings"] {
        assert!(
            corpus.manifest.get(key).is_none() || corpus.manifest[key].is_array(),
            "{label}: {key} is not of its type"
        );
    }
    assert_warning_texts(&label, &corpus.manifest["reference_warnings"]);
    let mut out = vec![default_run(corpus)];
    for set in corpus.manifest["flag_sets"]
        .as_array()
        .into_iter()
        .flatten()
    {
        let name = set["name"]
            .as_str()
            .unwrap_or_else(|| panic!("{label}: a flag set has a name"))
            .to_string();
        assert!(
            plain_name(&name),
            "{label}: flag set name {name:?} is not a plain word other than default"
        );
        assert!(
            out.iter().all(|r| r.name != name),
            "{label}: flag set {name} is named twice"
        );
        // A key that is present reads as its one type: a wrong type never falls
        // back to the default, which would run a set other than the one written.
        let typed = |key: &str, ok: bool| {
            assert!(
                set.get(key).is_none() || ok,
                "{label} flag set {name}: {key} is not of its type"
            );
        };
        typed("formats", set["formats"].is_array());
        typed("flags", set["flags"].is_array());
        typed("input", set["input"].is_string());
        typed("expected_dir", set["expected_dir"].is_string());
        typed("stdin", set["stdin"].is_boolean());
        typed("unordered", set["unordered"].is_boolean());
        typed("perl_mode", set["perl_mode"].is_string());
        typed("masks", set["masks"].is_array());
        typed("reference_warnings", set["reference_warnings"].is_array());
        let set_label = format!("{label} flag set {name}");
        assert_warning_texts(&set_label, &set["reference_warnings"]);
        let masks = if set.get("masks").is_some() {
            mask_classes(&set_label, &set["masks"])
        } else {
            corpus.masks()
        };
        let perl_mode = perl_mode_of(&set_label, set, &out[0].perl_mode);
        let formats = if set["formats"].is_array() {
            strings(&set["formats"])
        } else {
            ALL_FORMATS.iter().map(|f| f.to_string()).collect()
        };
        assert!(
            !formats.is_empty(),
            "{label} flag set {name}: no format to compare"
        );
        for format in &formats {
            assert!(
                ALL_FORMATS.contains(&format.as_str()),
                "{label} flag set {name}: unknown format {format}"
            );
        }
        let input = set["input"].as_str().unwrap_or("variants.vcf").to_string();
        let expected_dir = set["expected_dir"]
            .as_str()
            .map(str::to_string)
            .unwrap_or_else(|| format!("expected_{name}"));
        for path in [&input, &expected_dir] {
            assert!(
                inside_corpus(path),
                "{label} flag set {name}: {path:?} is not a path inside the corpus"
            );
        }
        let flags = corpus_flags(corpus, &set["flags"]);
        let stdin = set["stdin"].as_bool().unwrap_or(false);
        assert!(
            !stdin
                || flags
                    .iter()
                    .any(|f| f == "--format" || f.starts_with("--format=")),
            "{label} flag set {name}: a stdin run names its input format with --format"
        );
        out.push(Run {
            corpus,
            name,
            flags,
            input,
            stdin,
            formats,
            expected_dir,
            unordered: set["unordered"].as_bool().unwrap_or(false),
            perl_mode,
            masks,
            documented: if set.get("divergences").is_some() {
                set
            } else {
                &corpus.manifest
            },
            entry: set,
        });
    }
    out
}

/// The corpus's reference FASTA to pass as `--fasta`: the file the manifest's
/// `fasta` names beside it, with its indexes; `None` when the manifest names
/// none or the reference is staged inside the cache ([`Corpus::fasta_in_cache`]).
fn corpus_fasta(corpus: &Corpus) -> Option<PathBuf> {
    if corpus.fasta_in_cache() || !corpus.manifest["fasta"].is_string() {
        return None;
    }
    Some(committed_fasta(corpus))
}

/// A directory under `out_dir` holding the corpus FASTA and its indexes under a
/// name that carries the species and the assembly, as the files of an Ensembl
/// FASTA directory do, so a run can name it with `--fasta_dir`. Built once per
/// corpus; later calls reuse it.
fn staged_fasta_dir(corpus: &Corpus, out_dir: &Path) -> PathBuf {
    let dir = out_dir.join(format!("{}-fasta_dir", corpus.name));
    let fasta = committed_fasta(corpus);
    let name = fasta.file_name().unwrap().to_string_lossy().into_owned();
    let staged = dir.join(format!("{}.{}.{name}", corpus.species, corpus.assembly));
    if !staged.is_file() {
        fs::create_dir_all(&dir).unwrap();
        link_or_copy_indexed_fasta(&fasta, &staged);
    }
    dir
}

/// The placeholder a flag word or a `runs[].argv` token carries for the
/// directory the manifest's `generated` step writes ([`generated_dir`]).
const GENERATED_PLACEHOLDER: &str = "{generated}";

/// The directory the corpus's `generated` step writes, built once per corpus
/// under `out_dir` and reused by every later run in it. The manifest's
/// `generated` block names the `command` (its words, `{generated}` and
/// `{corpus}` expanded, run from the repository root) and the `sha256` of every
/// file it must write, each path relative to the directory: a fixture a corpus
/// does not commit (a matrix set a seeded generator reproduces) is written here
/// before the run and checked against the digests the reference was generated
/// over, so a generator whose output moved fails here, never as a comparison
/// against a reference built from the old output.
fn generated_dir(corpus: &Corpus, out_dir: &Path) -> PathBuf {
    let label = format!("{}/{}", corpus.release, corpus.name);
    let block = &corpus.manifest["generated"];
    assert!(
        block.is_object(),
        "{label}: a flag names {GENERATED_PLACEHOLDER} but the manifest has no `generated` block"
    );
    let dir = out_dir.join(format!("{}-generated", corpus.name));
    let marker = dir.join(".generated");
    if marker.is_file() {
        return dir;
    }
    fs::create_dir_all(&dir).unwrap();
    let words: Vec<String> = strings(&block["command"])
        .into_iter()
        .map(|word| {
            word.replace(GENERATED_PLACEHOLDER, &dir.display().to_string())
                .replace(CORPUS_PLACEHOLDER, &corpus.dir.display().to_string())
        })
        .collect();
    assert!(
        !words.is_empty(),
        "{label}: the `generated` block names no command"
    );
    let repository = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = Command::new(&words[0])
        .args(&words[1..])
        .current_dir(repository)
        .stdin(Stdio::null())
        .output()
        .unwrap_or_else(|e| panic!("{label}: `{}` does not start: {e}", words.join(" ")));
    assert!(
        output.status.success(),
        "{label}: `{}` failed:\n{}",
        words.join(" "),
        String::from_utf8_lossy(&output.stderr)
    );
    let digests = block["sha256"]
        .as_object()
        .unwrap_or_else(|| panic!("{label}: the `generated` block has no `sha256` map"));
    assert!(
        !digests.is_empty(),
        "{label}: the `generated` block's `sha256` map is empty"
    );
    for (name, expected) in digests {
        let path = dir.join(name);
        let digest = vep_core::cache_manifest::file_digest(&path)
            .unwrap_or_else(|e| panic!("{label}: generated file {name} cannot be read: {e}"));
        let actual: String = digest.iter().map(|byte| format!("{byte:02x}")).collect();
        assert_eq!(
            Some(actual.as_str()),
            expected.as_str(),
            "{label}: generated file {name} does not hash as the manifest records"
        );
    }
    fs::write(marker, b"").unwrap();
    dir
}

/// `words` with `{generated}` inside any word replaced by the corpus's
/// [`generated_dir`], built on the first word that names it.
fn expand_generated(corpus: &Corpus, words: &[String], out_dir: &Path) -> Vec<String> {
    let mut generated: Option<String> = None;
    words
        .iter()
        .map(|word| {
            if !word.contains(GENERATED_PLACEHOLDER) {
                return word.clone();
            }
            let dir = generated
                .get_or_insert_with(|| generated_dir(corpus, out_dir).display().to_string());
            word.replace(GENERATED_PLACEHOLDER, dir)
        })
        .collect()
}

/// The arguments every run of a corpus starts from, after `-i <input>`: the
/// corpus cache (the staged copy under `out_dir` for a `fasta_in_cache` corpus),
/// the fixed arguments, the run's `flags` (a word naming `{generated}` reads the
/// directory the manifest's `generated` step writes under `out_dir`) and the
/// corpus `fasta`. A run whose flags carry `--grch37` names its assembly through
/// that flag alone, so the option set, not `--assembly`, selects it.
fn base_args(run: &Run, out_dir: &Path) -> Vec<String> {
    let corpus = run.corpus;
    let mut args = vec!["--offline".to_string()];
    args.extend(corpus.run_cache_args(out_dir));
    args.extend(["--species".to_string(), corpus.species.clone()]);
    if !run.flags.iter().any(|f| f == "--grch37") {
        args.extend(["--assembly".to_string(), corpus.assembly.clone()]);
    }
    args.extend(
        [
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ]
        .map(str::to_string),
    );
    args.extend(expand_generated(corpus, &run.flags, out_dir));
    if let Some(fasta) = corpus_fasta(corpus) {
        args.extend(["--fasta".to_string(), fasta.display().to_string()]);
    }
    args
}

/// A cache under `dir` with no transcript on chromosome 21, the regulatory
/// fixture shard (`fixtures/regulatory/1-1000000.json`: an enhancer at 100-120
/// and two motifs inside it) on chromosome 21, and a corrupt regulatory shard
/// on chromosome 22, which the command line's lazy loader reads only for a
/// record of that chromosome.
pub fn regulatory_cache(dir: &Path) -> PathBuf {
    regulatory_cache_with(dir, include_bytes!("../fixtures/regulatory/1-1000000.json"))
}

/// The `info.json` of a cache holding regulatory shards in the fixtures' shape:
/// the regulatory marker, the shard schema and the fixtures' cell-type list,
/// which every fixture feature's activity string is aligned to.
pub fn regulatory_info_json() -> String {
    let cell_types: serde_json::Value =
        serde_json::from_slice(include_bytes!("../fixtures/regulatory/cell_types.json")).unwrap();
    serde_json::json!({
        "assembly": "GRCh38",
        "regulatory": true,
        "regulatory_schema": 1,
        "cell_types": cell_types,
    })
    .to_string()
}

/// Writes the fixtures' binding matrices table into `cache`'s regulatory tree.
pub fn write_regulatory_matrices(cache: &Path) {
    fs::create_dir_all(cache.join("regulatory")).unwrap();
    fs::write(
        cache.join("regulatory/binding_matrices.json"),
        include_bytes!("../fixtures/regulatory/binding_matrices.json"),
    )
    .unwrap();
}

/// [`regulatory_cache`] over the slice of the converted chromosome 21 cache
/// (`fixtures/regulatory/chr21-slice.json`, described in the fixtures' README):
/// seven features as the converter wrote them.
pub fn regulatory_cache_of_real_features(dir: &Path) -> PathBuf {
    regulatory_cache_with(
        dir,
        include_bytes!("../fixtures/regulatory/chr21-slice.json"),
    )
}

fn regulatory_cache_with(dir: &Path, shard: &[u8]) -> PathBuf {
    let cache = dir.join("cache");
    fs::create_dir_all(cache.join("transcripts/21")).unwrap();
    fs::create_dir_all(cache.join("regulatory/21")).unwrap();
    fs::create_dir_all(cache.join("regulatory/22")).unwrap();
    fs::write(cache.join("transcripts/21/1-1000000.json"), "[]").unwrap();
    fs::write(cache.join("regulatory/21/1-50000000.json"), shard).unwrap();
    fs::write(cache.join("regulatory/22/1-1000000.json"), "[").unwrap();
    write_regulatory_matrices(&cache);
    fs::write(cache.join("info.json"), regulatory_info_json()).unwrap();
    cache
}

/// The binary with no `VEP_*` variable and no `vep.ini` of the machine's own
/// reaching it: every `VEP_*` variable is removed and `HOME` is a directory
/// under `out_dir` holding an empty `.vep/Plugins`, so a run reads only the
/// flags the test gives and its default plugins directory exists, as the one
/// the reference's Ensembl VEP image names does, so a `--plugin` run warns of
/// no missing directory.
pub fn isolated_vep(out_dir: &Path) -> Command {
    let home = out_dir.join("home");
    fs::create_dir_all(home.join(".vep").join("Plugins")).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    for (key, _) in std::env::vars_os() {
        if key.to_string_lossy().starts_with("VEP_") {
            cmd.env_remove(&key);
        }
    }
    cmd.env("HOME", home);
    cmd
}

/// [`isolated_vep`] with its `HOME` under the test target directory, for a test
/// that spawns the binary outside any output directory.
pub fn isolated_vep_anywhere() -> Command {
    isolated_vep(Path::new(env!("CARGO_TARGET_TMPDIR")))
}

/// The device that refuses every write, where the platform has one.
const FULL_DEVICE: &str = "/dev/full";

/// The streams that take no bytes, by name, for a test that puts standard
/// error on each: the device that refuses every write, where the platform has
/// one, and a pipe whose reader has gone.
pub fn unwritable_streams() -> &'static [&'static str] {
    if Path::new(FULL_DEVICE).exists() {
        &[FULL_DEVICE, "a pipe whose reader has gone"]
    } else {
        &["a pipe whose reader has gone"]
    }
}

/// The stream [`unwritable_streams`] names as `stream`, opened for a child
/// process to write.
pub fn unwritable(stream: &str) -> Stdio {
    if stream == FULL_DEVICE {
        return fs::File::options().write(true).open(stream).unwrap().into();
    }
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

/// The `vep` invocation a run starts from: the run's input (`-i STDIN` when it
/// is piped) and the corpus cache, the fixed arguments, then the run's `flags`
/// and the corpus `fasta`; `out_dir` holds the staged cache of a
/// `fasta_in_cache` corpus. The caller adds the output path and format and runs
/// it through [`run_output`], which pipes the input when the run says so. The
/// process runs in the corpus directory, so a relative path inside a flag (a
/// `--transcript_filter` list file beside the manifest) names a corpus file.
/// `VEP_DIR_PLUGINS` names the corpus directory, as the reference's container
/// names an existing plugin directory through that variable, so a run under
/// `--plugin` warns of no missing directory under the isolated home.
pub fn vep_command(run: &Run, out_dir: &Path) -> Command {
    let corpus = run.corpus;
    let mut cmd = isolated_vep(out_dir);
    cmd.current_dir(&corpus.dir);
    cmd.env("VEP_DIR_PLUGINS", &corpus.dir);
    cmd.arg("-i");
    if run.stdin {
        cmd.arg("STDIN").stdin(Stdio::piped());
    } else {
        cmd.arg(corpus.dir.join(&run.input));
    }
    cmd.args(base_args(run, out_dir));
    cmd
}

/// Runs `cmd` to completion with `input` written on its standard input, its
/// standard output and error captured as [`Command::output`] captures them.
fn output_with_stdin(cmd: &mut Command, input: Vec<u8>) -> Output {
    let mut child = cmd
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("vep binary runs");
    let mut stdin = child.stdin.take().expect("stdin is piped");
    let writer = std::thread::spawn(move || stdin.write_all(&input));
    let output = child.wait_with_output().expect("vep binary runs");
    let written = writer.join().unwrap();
    // A child that exited before reading its input leaves the writer a broken
    // pipe; the caller reports that child's stderr, so the write is checked
    // only for a child that succeeded.
    if output.status.success() {
        written.expect("input written to vep");
    }
    output
}

/// Runs a [`vep_command`] to completion, feeding the run's input file on
/// standard input when the run is a stdin run.
pub fn run_output(run: &Run, cmd: &mut Command) -> Output {
    if !run.stdin {
        return cmd.output().expect("vep binary runs");
    }
    let input = fs::read(run.corpus.dir.join(&run.input)).unwrap();
    output_with_stdin(cmd, input)
}

/// The output path [`run_vep`] writes for a run in one format; the run's
/// warnings file sits beside it as `<path>_warnings.txt`.
pub fn output_path(run: &Run, format: &str, out_dir: &Path) -> PathBuf {
    out_dir.join(format!("{}-{}-{format}.out", run.corpus.name, run.name))
}

/// Runs the `vep` binary for the run in one output format and returns the
/// output file's text; a run under `--compress_output` is decompressed first
/// (gzip or BGZF, both multi-member gzip streams), since the compressed bytes
/// are the encoder's and the text is what the reference compares. `format` is
/// `default`, `tab`, `vcf` or `json`.
pub fn run_vep(run: &Run, format: &str, out_dir: &Path, extra: &[&str]) -> String {
    let out = output_path(run, format, out_dir);
    let mut cmd = vep_command(run, out_dir);
    cmd.arg("-o").arg(&out);
    match format {
        "default" => {}
        "tab" | "vcf" | "json" => {
            cmd.arg(format!("--{format}"));
        }
        other => panic!("unknown format {other}"),
    }
    cmd.args(extra);
    let output = run_output(run, &mut cmd);
    assert!(
        output.status.success(),
        "vep --{format} failed on {}:\n{}",
        run.label(),
        String::from_utf8_lossy(&output.stderr)
    );
    if run.has_flag("--compress_output") {
        let bytes = fs::read(&out).unwrap();
        assert!(
            bytes.starts_with(&[0x1f, 0x8b]),
            "{}: --compress_output wrote no gzip stream",
            run.label()
        );
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .unwrap();
        return text;
    }
    fs::read_to_string(&out).unwrap()
}

/// The corpus input's data lines (`variants.vcf`), in order.
pub fn input_records(corpus: &Corpus) -> Vec<String> {
    data_lines(&corpus.dir.join("variants.vcf"))
}

// --------------------------------------------------------------------------
// Entries: one comparable unit per (record, allele, feature).

/// A comparable unit of output: one row of the default or tab format, one CSQ
/// entry of the VCF format, one consequence object of the JSON format, or one
/// record-level object.
#[derive(Debug, Clone)]
pub struct Entry {
    /// Stable identity within the format: the (Location or input record
    /// ordinal, Allele, Feature, Feature_type) comparison key plus an occurrence
    /// counter for exact duplicates.
    pub key: String,
    /// Corpus record ordinal when the format makes it known.
    pub record: Option<usize>,
    /// `Location` when the format carries it.
    pub location: Option<String>,
    pub allele: String,
    pub feature: String,
    /// Consequence terms as the format prints them (`,` or `&` joined).
    pub consequence: String,
    pub fields: BTreeMap<String, String>,
}

fn dedup_keys(entries: &mut [Entry]) {
    let mut seen: HashMap<String, usize> = HashMap::new();
    for e in entries.iter_mut() {
        let n = seen.entry(e.key.clone()).or_insert(0);
        *n += 1;
        if *n > 1 {
            e.key = format!("{}#{}", e.key, n);
        }
    }
}

/// Numbers the entries that share a key by what they hold, where
/// [`dedup_keys`] numbered them by where they stand: of the entries under one
/// key, the one whose consequence and fields sort first keeps the key and the
/// rest take `#2`, `#3` in that order. An unordered run pairs its entries
/// this way, so two rows of one feature that the engines write in different
/// orders (a transcript the cache and a gene-model file both hold, which
/// Ensembl VEP lists in an order that differs from run to run) meet their own
/// counterparts.
pub fn number_shared_keys_by_content(entries: &mut [Entry]) {
    fn base(key: &str) -> &str {
        match key.rsplit_once('#') {
            Some((base, count))
                if !count.is_empty() && count.bytes().all(|b| b.is_ascii_digit()) =>
            {
                base
            }
            _ => key,
        }
    }
    let mut sharing: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (at, entry) in entries.iter().enumerate() {
        sharing
            .entry(base(&entry.key).to_string())
            .or_default()
            .push(at);
    }
    for (key, mut places) in sharing {
        if places.len() < 2 {
            continue;
        }
        places.sort_by(|&a, &b| {
            (&entries[a].consequence, &entries[a].fields)
                .cmp(&(&entries[b].consequence, &entries[b].fields))
        });
        for (rank, at) in places.into_iter().enumerate() {
            entries[at].key = match rank {
                0 => key.clone(),
                _ => format!("{key}#{}", rank + 1),
            };
        }
    }
}

/// A JSON object written with its keys sorted, or `text` as it is when it is
/// not JSON. Ensembl VEP encodes the `GA4GH_VRS` object with its keys in hash
/// order, which differs from row to row, so the text formats compare the parsed
/// object, as the JSON format's nested object does.
fn canonical_json(text: &str) -> String {
    serde_json::from_str::<serde_json::Value>(text)
        .map(|value| value.to_string())
        .unwrap_or_else(|_| text.to_string())
}

/// Parses the default format: 13 named columns plus `Extra`, whose `key=value`
/// pairs become fields of their own (and `Extra` keeps the raw text, its
/// `GA4GH_VRS` value in [`canonical_json`] form).
pub fn parse_default(text: &str) -> Vec<Entry> {
    let mut cols: Vec<String> = Vec::new();
    let mut entries = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('#') {
            if !line.starts_with("##") {
                cols = h.split('\t').map(str::to_string).collect();
            }
            continue;
        }
        let vals: Vec<&str> = line.split('\t').collect();
        assert_eq!(vals.len(), cols.len(), "column count on line: {line}");
        let mut fields: BTreeMap<String, String> = cols
            .iter()
            .cloned()
            .zip(vals.iter().map(|v| v.to_string()))
            .collect();
        if let Some(vrs) = fields.get_mut("GA4GH_VRS") {
            *vrs = canonical_json(vrs);
        }
        if let Some(extra) = fields.get("Extra").cloned() {
            if extra != "-" {
                let mut extra_text = extra.clone();
                for kv in extra.split(';') {
                    let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
                    let mut v = v.to_string();
                    if k == "GA4GH_VRS" {
                        v = canonical_json(&v);
                        extra_text = extra_text.replace(kv, &format!("{k}={v}"));
                    }
                    fields.insert(k.to_string(), v);
                }
                fields.insert("Extra".to_string(), extra_text);
            }
        }
        let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
        entries.push(Entry {
            key: format!(
                "{}|{}|{}|{}",
                get("Location"),
                get("Allele"),
                get("Feature"),
                get("Feature_type")
            ),
            record: None,
            location: Some(get("Location")),
            allele: get("Allele"),
            feature: get("Feature"),
            consequence: get("Consequence"),
            fields,
        });
    }
    dedup_keys(&mut entries);
    entries
}

/// Parses the tab format (every column named in the `#` line).
pub fn parse_tab(text: &str) -> Vec<Entry> {
    parse_default(text)
}

/// The corpus input records, for mapping an output record back to its ordinal:
/// by identity where the output format carries one, so that a format which
/// omits a skipped record (VEP writes no VCF line and no JSON object for one)
/// keeps every later record aligned, identical input lines matched in order
/// of appearance; in turn where it carries none (the VCF output of an input
/// that is not a VCF), which holds only while the output has one record per
/// input record, so [`InputIndex::records_in_turn`] is checked against the
/// line count once the output is read.
struct InputIndex {
    unmatched: HashMap<String, VecDeque<usize>>,
    /// The next ordinal and the record count of an index pairing output
    /// records with the input records in turn, whatever identity they carry.
    in_order: Option<(usize, usize)>,
}

impl InputIndex {
    /// `key` reduces an input line to the identity the output format carries.
    fn new(inputs: &[String], key: impl Fn(&str) -> String) -> Self {
        let mut unmatched: HashMap<String, VecDeque<usize>> = HashMap::new();
        for (i, line) in inputs.iter().enumerate() {
            unmatched.entry(key(line)).or_default().push_back(i);
        }
        InputIndex {
            unmatched,
            in_order: None,
        }
    }

    /// An index pairing each output record with the next of `records` input
    /// records in turn, for an input whose lines the output format carries no
    /// identity of: the VCF output of an Ensembl-format, region or HGVS input,
    /// where both engines write one line per input record in input order.
    fn in_order(records: usize) -> Self {
        InputIndex {
            unmatched: HashMap::new(),
            in_order: Some((0, records)),
        }
    }

    /// The ordinal of the next unmatched input record with this identity; in
    /// turn, the next record, none once every record is paired.
    fn record(&mut self, key: &str) -> Option<usize> {
        if let Some((next, records)) = self.in_order.as_mut() {
            if *next >= *records {
                return None;
            }
            *next += 1;
            return Some(*next - 1);
        }
        self.unmatched.get_mut(key)?.pop_front()
    }

    /// In turn, the input record count the output's must equal, since a
    /// shorter or longer output is a misalignment the pairing cannot tell from
    /// an exact match. `None` by identity.
    fn records_in_turn(&self) -> Option<usize> {
        self.in_order.map(|(_, records)| records)
    }
}

/// A VCF data line's identity: CHROM, POS, ID, REF and ALT, which VEP and
/// vep-rs both write back unchanged.
pub fn vcf_line_key(line: &str) -> String {
    line.split('\t').take(5).collect::<Vec<_>>().join("\t")
}

/// Parses the VCF format: one entry per CSQ chunk, keyed by the input record
/// ordinal (the line matched to `inputs` by CHROM, POS, ID, REF and ALT),
/// `Allele` and `Feature`, plus one record-level entry per line carrying the
/// eight fixed columns with the CSQ removed from INFO. A line matching no input
/// record is keyed by its output ordinal.
pub fn parse_vcf(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_vcf_info(text, inputs, "CSQ")
}

/// [`parse_vcf`] for an output whose annotation is under `info_field`
/// (`--vcf_info_field`). The chunks parsed are the line's last entry under that
/// key, whose subfields the last `##INFO=<ID=<field>,` header line names (the
/// annotation this run wrote); an earlier entry kept by `--keep_csq`, and the
/// earlier header lines, stay in `INFO_without_CSQ` and the header comparison.
pub fn parse_vcf_info(text: &str, inputs: &[String], info_field: &str) -> Vec<Entry> {
    parse_vcf_entries(text, InputIndex::new(inputs, vcf_line_key), info_field)
}

/// [`parse_vcf_info`] for the VCF output of an input that is not a VCF: its
/// lines carry no input line to match, so each is the next of `inputs` in
/// turn, and the output must hold one line per input record.
pub fn parse_vcf_info_in_order(text: &str, inputs: &[String], info_field: &str) -> Vec<Entry> {
    parse_vcf_entries(text, InputIndex::in_order(inputs.len()), info_field)
}

fn parse_vcf_entries(text: &str, mut index: InputIndex, info_field: &str) -> Vec<Entry> {
    let info_header = format!("##INFO=<ID={info_field},");
    let info_key = format!("{info_field}=");
    let mut csq_fields: Vec<String> = Vec::new();
    let mut entries = Vec::new();
    let mut lines_seen = 0usize;
    for line in text.lines() {
        if line.starts_with(&info_header) {
            let fmt = line.split("Format: ").nth(1).unwrap();
            let fmt = fmt.trim_end_matches("\">");
            csq_fields = fmt.split('|').map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(cols.len() >= 8, "VCF line has fewer than 8 columns: {line}");
        let record = index.record(&vcf_line_key(line));
        let rec_label = record
            .map(|r| format!("rec{r}"))
            .unwrap_or(format!("line{lines_seen}"));
        lines_seen += 1;
        let kvs: Vec<&str> = cols[7].split(';').collect();
        let annotation = kvs.iter().rposition(|kv| kv.starts_with(&info_key));
        let csq = annotation.map(|i| kvs[i][info_key.len()..].to_string());
        let info_rest: Vec<String> = kvs
            .iter()
            .enumerate()
            .filter(|(i, _)| Some(*i) != annotation)
            .map(|(_, kv)| kv.to_string())
            .collect();
        let mut fields = BTreeMap::new();
        for (name, val) in ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
            .iter()
            .zip(cols.iter())
        {
            fields.insert(name.to_string(), val.to_string());
        }
        fields.insert("INFO_without_CSQ".to_string(), info_rest.join(";"));
        fields.insert(
            "CSQ_present".to_string(),
            if csq.is_some() { "yes" } else { "no" }.to_string(),
        );
        entries.push(Entry {
            key: format!("{rec_label}|record"),
            record,
            location: None,
            allele: String::new(),
            feature: String::new(),
            consequence: String::new(),
            fields,
        });
        if let Some(csq) = csq {
            for chunk in csq.split(',') {
                let vals: Vec<&str> = chunk.split('|').collect();
                assert_eq!(
                    vals.len(),
                    csq_fields.len(),
                    "CSQ chunk width on {rec_label}: {chunk}"
                );
                let fields: BTreeMap<String, String> = csq_fields
                    .iter()
                    .cloned()
                    .zip(vals.iter().map(|v| v.to_string()))
                    .collect();
                let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
                entries.push(Entry {
                    key: format!(
                        "{rec_label}|{}|{}|{}",
                        get("Allele"),
                        get("Feature"),
                        get("Feature_type")
                    ),
                    record,
                    location: None,
                    allele: get("Allele"),
                    feature: get("Feature"),
                    consequence: get("Consequence"),
                    fields,
                });
            }
        }
    }
    if let Some(records) = index.records_in_turn() {
        assert_eq!(
            lines_seen, records,
            "an output of an input that is not a VCF pairs its lines with the \
             input records in turn, so it needs one line per record: {lines_seen} \
             data lines against {records} records"
        );
    }
    dedup_keys(&mut entries);
    entries
}

fn json_scalar(v: &serde_json::Value) -> String {
    match v {
        serde_json::Value::String(s) => s.clone(),
        serde_json::Value::Array(items) => {
            items.iter().map(json_scalar).collect::<Vec<_>>().join(",")
        }
        other => other.to_string(),
    }
}

/// The JSON arrays of consequence objects, each an entry of its own: the array,
/// the entry kind and the key naming the feature (none for the intergenic row).
const JSON_CONSEQUENCE_ARRAYS: [(&str, &str, Option<&str>); 4] = [
    (
        "transcript_consequences",
        "transcript",
        Some("transcript_id"),
    ),
    (
        "regulatory_feature_consequences",
        "regulatory",
        Some("regulatory_feature_id"),
    ),
    (
        "motif_feature_consequences",
        "motif",
        Some("motif_feature_id"),
    ),
    ("intergenic_consequences", "intergenic", None),
];

/// Parses the JSON format: one record-level entry per object (top-level
/// scalars, keyed by the corpus record the `input` line names, the vep-rs
/// identity key set aside for [`compare_json_identity`]) and one entry per
/// element of the four consequence arrays ([`JSON_CONSEQUENCE_ARRAYS`]), the
/// regulatory and motif ones keyed by their feature id as the transcript ones
/// are, so a documented field divergence on a motif row applies to its JSON
/// object. An object names its record by the first five columns of its
/// `input`, as a VCF line does: Ensembl VEP stores a structural record's line
/// with QUAL, FILTER and INFO filled (Parser/VCF.pm:464-466), so the whole
/// `input` may differ from the line as written. An object whose `input` matches
/// no record is keyed by its output ordinal.
pub fn parse_json(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut index = InputIndex::new(inputs, vcf_line_key);
    let mut entries = Vec::new();
    for (n, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let obj: serde_json::Value = serde_json::from_str(line).unwrap();
        let input = obj.get("input").and_then(|v| v.as_str()).unwrap_or("");
        let record = index.record(&vcf_line_key(input));
        let rec_label = record
            .map(|r| format!("rec{r}"))
            .unwrap_or(format!("obj{n}"));
        let mut fields = BTreeMap::new();
        for (k, v) in obj.as_object().unwrap() {
            if JSON_CONSEQUENCE_ARRAYS
                .iter()
                .any(|(array, _, _)| array == k)
                || k == "input"
                || k == "vep_rs_version"
            {
                continue;
            }
            fields.insert(k.clone(), json_scalar(v));
        }
        entries.push(Entry {
            key: format!("{rec_label}|record"),
            record,
            location: None,
            allele: String::new(),
            feature: String::new(),
            consequence: String::new(),
            fields,
        });
        for (list, kind, feature_key) in JSON_CONSEQUENCE_ARRAYS {
            let Some(items) = obj.get(list).and_then(|v| v.as_array()) else {
                continue;
            };
            for item in items {
                let fields: BTreeMap<String, String> = item
                    .as_object()
                    .unwrap()
                    .iter()
                    .map(|(k, v)| (k.clone(), json_scalar(v)))
                    .collect();
                let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
                let feature = feature_key.map(get).unwrap_or_default();
                entries.push(Entry {
                    key: format!("{rec_label}|{}|{}|{kind}", get("variant_allele"), feature),
                    record,
                    location: None,
                    allele: get("variant_allele"),
                    feature,
                    consequence: get("consequence_terms"),
                    fields,
                });
            }
        }
    }
    dedup_keys(&mut entries);
    entries
}

/// Whether a run with the argument words `words` over the cache directory
/// `cache_dir` has a reference sequence, as the binary settles it
/// (`fasta_source` of `runner.rs`): the file `--fasta_dir` holds for the
/// corpus's species and assembly, else a `--fasta` word (or `--fa`), else,
/// unless `--no_fasta` is given, the first indexed `.fa` or `.fa.gz` of the
/// cache directory, which a `fasta_in_cache` corpus stages there.
pub fn has_reference_sequence(words: &[String], cache_dir: &Path, corpus: &Corpus) -> bool {
    let named = |flag: &str| {
        words
            .iter()
            .any(|w| w == flag || w.starts_with(&format!("{flag}=")))
    };
    let value = |flag: &str| {
        let prefix = format!("{flag}=");
        words.iter().enumerate().find_map(|(i, w)| {
            if let Some(value) = w.strip_prefix(&prefix) {
                Some(value)
            } else if w == flag {
                words.get(i + 1).map(String::as_str)
            } else {
                None
            }
        })
    };
    let in_fasta_dir = value("--fasta_dir").is_some_and(|dir| {
        matches!(
            vep_fasta::fasta_dir_match(Path::new(dir), &corpus.species, Some(&corpus.assembly)),
            Ok(Some(_))
        )
    });
    in_fasta_dir
        || named("--fasta")
        || named("--fa")
        || (!named("--no_fasta") && vep_fasta::first_indexed_fasta(cache_dir).is_some())
}

/// The key under which a `divergences` item states vep-rs's consequence sets
/// in a run that has a reference sequence.
pub const WITH_REFERENCE_SETS: &str = "vep_rs_consequence_sets_with_reference_sequence";

/// The sets a `divergences` item states under [`WITH_REFERENCE_SETS`], `None`
/// for an item without the key. The key is a non-empty list of distinct sets
/// in sorted order, each spelt as [`consequence_set`] spells it and other than
/// the item's `vep_rs_consequence_sets` as a whole; a list that holds VEP's set
/// (`vep_consequence_set`) holds nothing else. Any other shape fails the test
/// that reads it, whatever the run's setting.
fn with_reference_sets(item: &serde_json::Value) -> Option<Vec<String>> {
    let value = item.get(WITH_REFERENCE_SETS)?;
    let key = format!(
        "{} {} {}",
        item["location"].as_str().unwrap_or(""),
        item["allele"].as_str().unwrap_or(""),
        item["feature"].as_str().unwrap_or("")
    );
    fn refuse(key: &str, reason: &str) -> ! {
        panic!("the item on {key}: {WITH_REFERENCE_SETS} {reason}")
    }
    let sets: Vec<String> = match value.as_array() {
        Some(list) if !list.is_empty() => list
            .iter()
            .map(|s| match s.as_str() {
                Some(s) if !s.is_empty() && consequence_set(s) == s => s.to_string(),
                _ => refuse(
                    &key,
                    "holds an item that is not a sorted, comma-joined set of terms",
                ),
            })
            .collect(),
        _ => refuse(&key, "is not a non-empty list"),
    };
    if !sets.windows(2).all(|pair| pair[0] < pair[1]) {
        refuse(&key, "is not a sorted list of distinct sets");
    }
    let without: Option<Vec<&str>> = item["vep_rs_consequence_sets"]
        .as_array()
        .map(|list| list.iter().filter_map(|s| s.as_str()).collect());
    if without.is_some_and(|without| without == sets) {
        refuse(
            &key,
            "repeats vep_rs_consequence_sets and states no second outcome",
        );
    }
    let vep = consequence_set(item["vep_consequence_set"].as_str().unwrap_or(""));
    if sets.len() > 1 && sets.contains(&vep) {
        refuse(&key, "names VEP's set beside another");
    }
    Some(sets)
}

// --------------------------------------------------------------------------
// Documented divergences.

/// The keys whose consequence terms are documented to differ, from the
/// manifest, addressable by Location or by record ordinal, and the records the
/// reference skipped.
pub struct Documented {
    /// (Location, Allele, Feature) -> documented vep-rs consequence sets.
    by_location: HashMap<(String, String, String), Vec<String>>,
    /// (record ordinal, Allele, Feature) -> documented vep-rs consequence sets.
    by_record: HashMap<(usize, String, String), Vec<String>>,
    /// Keys only vep-rs emits.
    extra_by_location: BTreeSet<(String, String, String)>,
    extra_by_record: BTreeSet<(usize, String, String)>,
    /// (Location, Allele, Feature, consequence set) -> field name -> the value
    /// vep-rs is documented to print where VEP prints another. The consequence
    /// set tells apart two records that share a Location and Allele (an
    /// insertion and a substitution of the same bases), one of whose rows on
    /// the feature carries the field and the other not.
    fields_by_location: HashMap<(String, String, String, String), BTreeMap<String, String>>,
    fields_by_record: HashMap<(usize, String, String, String), BTreeMap<String, String>>,
    /// The documented fields, under the same two keys, of an item that states
    /// `rows_alike`: several records write a row on its key, and on some of
    /// those rows both engines print one value. Each field names its item's
    /// place in `shared_items`.
    shared_fields_by_location: HashMap<(String, String, String, String), BTreeMap<String, usize>>,
    shared_fields_by_record: HashMap<(usize, String, String, String), BTreeMap<String, usize>>,
    /// The items that state `rows_alike`, in the manifest's order.
    shared_items: Vec<SharedField>,
    /// Records the reference wrote nothing for: those a `vep_rs_only_tuples`
    /// entry of the class [`REFERENCE_SKIPPED_RECORD`] names.
    skipped_records: BTreeSet<usize>,
}

/// The class a `vep_rs_only_tuples` entry carries when every record it names
/// is one the reference dropped before annotation.
pub const REFERENCE_SKIPPED_RECORD: &str = "reference_skipped_record";

/// The key of a `field_divergences` item that counts the rows of its key on
/// which both engines print one value of the field: the rows of other records
/// that share the key with the documented one.
pub const ROWS_ALIKE: &str = "rows_alike";

/// A `field_divergences` item that states [`ROWS_ALIKE`]: its key as a report
/// names it, its field and the count it states.
struct SharedField {
    key: String,
    field: String,
    rows_alike: u64,
}

impl Documented {
    /// The documented keys of a run without a reference sequence.
    pub fn from_manifest(manifest: &serde_json::Value) -> Self {
        Self::for_run(manifest, false)
    }

    /// The documented keys of a run, `with_reference` saying whether the binary
    /// has a reference sequence in it. A `divergences` item documents its
    /// `vep_rs_consequence_sets`; in a run with a reference sequence an item that
    /// carries [`WITH_REFERENCE_SETS`] documents those sets instead, and where
    /// they are VEP's own set it documents nothing: its key is compared like any
    /// other.
    pub fn for_run(manifest: &serde_json::Value, with_reference: bool) -> Self {
        let mut d = Documented {
            by_location: HashMap::new(),
            by_record: HashMap::new(),
            extra_by_location: BTreeSet::new(),
            extra_by_record: BTreeSet::new(),
            fields_by_location: HashMap::new(),
            fields_by_record: HashMap::new(),
            shared_fields_by_location: HashMap::new(),
            shared_fields_by_record: HashMap::new(),
            shared_items: Vec::new(),
            skipped_records: BTreeSet::new(),
        };
        let s = |v: &serde_json::Value, k: &str| {
            v.get(k).and_then(|x| x.as_str()).unwrap_or("").to_string()
        };
        for item in manifest["divergences"].as_array().into_iter().flatten() {
            let stated = with_reference_sets(item);
            let sets: Vec<String> = match stated {
                Some(stated) if with_reference => {
                    let vep = consequence_set(item["vep_consequence_set"].as_str().unwrap_or(""));
                    if stated == [vep] {
                        continue;
                    }
                    stated
                }
                _ => item["vep_rs_consequence_sets"]
                    .as_array()
                    .into_iter()
                    .flatten()
                    .filter_map(|x| x.as_str().map(str::to_string))
                    .collect(),
            };
            d.by_location
                .entry((s(item, "location"), s(item, "allele"), s(item, "feature")))
                .or_default()
                .extend(sets.iter().cloned());
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    d.by_record
                        .entry((r as usize, s(item, "allele"), s(item, "feature")))
                        .or_default()
                        .extend(sets.iter().cloned());
                }
            }
        }
        for item in manifest["vep_rs_only_tuples"]
            .as_array()
            .into_iter()
            .flatten()
        {
            d.extra_by_location.insert((
                s(item, "location"),
                s(item, "allele"),
                s(item, "feature"),
            ));
            let reference_skipped =
                item["expected_divergence"].as_str() == Some(REFERENCE_SKIPPED_RECORD);
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    d.extra_by_record
                        .insert((r as usize, s(item, "allele"), s(item, "feature")));
                    if reference_skipped {
                        d.skipped_records.insert(r as usize);
                    }
                }
            }
        }
        for item in manifest["field_divergences"]
            .as_array()
            .into_iter()
            .flatten()
        {
            let (field, value) = (s(item, "field"), s(item, "vep_rs"));
            let terms = consequence_set(&s(item, "consequence_set"));
            let location_key = (
                s(item, "location"),
                s(item, "allele"),
                s(item, "feature"),
                terms.clone(),
            );
            let shared = item.get(ROWS_ALIKE).map(|count| {
                d.shared_items.push(SharedField {
                    key: format!("{}|{}|{}", location_key.0, location_key.1, location_key.2),
                    field: field.clone(),
                    rows_alike: count.as_u64().unwrap_or(0),
                });
                d.shared_items.len() - 1
            });
            if let Some(place) = shared {
                d.shared_fields_by_location
                    .entry(location_key.clone())
                    .or_default()
                    .insert(field.clone(), place);
            }
            d.fields_by_location
                .entry(location_key)
                .or_default()
                .insert(field.clone(), value.clone());
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    let record_key = (
                        r as usize,
                        s(item, "allele"),
                        s(item, "feature"),
                        terms.clone(),
                    );
                    if let Some(place) = shared {
                        d.shared_fields_by_record
                            .entry(record_key.clone())
                            .or_default()
                            .insert(field.clone(), place);
                    }
                    d.fields_by_record
                        .entry(record_key)
                        .or_default()
                        .insert(field.clone(), value.clone());
                }
            }
        }
        d
    }

    /// The documented field values on an entry's key: those under its Location
    /// when it has one and that key is documented, else those under its record.
    fn divergent_fields(&self, e: &Entry) -> Option<&BTreeMap<String, String>> {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        let terms = consequence_set(&e.consequence);
        if let Some(loc) = &e.location {
            if let Some(m) = self.fields_by_location.get(&(
                loc.clone(),
                allele.clone(),
                feature.clone(),
                terms.clone(),
            )) {
                return Some(m);
            }
        }
        let record = e.record?;
        self.fields_by_record.get(&(record, allele, feature, terms))
    }

    /// The documented vep-rs value of `field` for an entry, when the field is
    /// documented to differ on its key. `field` is matched as the format spells
    /// it: the manifest names the default-format key (`HGVSp`), the JSON format
    /// its lower-case form.
    pub fn divergent_field(&self, e: &Entry, field: &str) -> Option<&String> {
        fn lookup<'m>(m: &'m BTreeMap<String, String>, field: &str) -> Option<&'m String> {
            m.get(field).or_else(|| {
                m.iter()
                    .find(|(k, _)| k.eq_ignore_ascii_case(field))
                    .map(|(_, v)| v)
            })
        }
        self.divergent_fields(e).and_then(|m| lookup(m, field))
    }

    /// The place among the items that state [`ROWS_ALIKE`] of the one that
    /// documents `field` on the entry's key, the key looked up as
    /// [`Documented::divergent_field`] looks it up and the field in any letter
    /// case; `None` where no such item documents it. On such a key a row whose
    /// two values agree is the row of another record than the documented one.
    pub fn field_shares_its_key(&self, e: &Entry, field: &str) -> Option<usize> {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        let terms = consequence_set(&e.consequence);
        let named = |fields: &BTreeMap<String, usize>| {
            fields
                .iter()
                .find(|(name, _)| name.eq_ignore_ascii_case(field))
                .map(|(_, place)| *place)
        };
        if let Some(loc) = &e.location {
            let key = (loc.clone(), allele.clone(), feature.clone(), terms.clone());
            if self.fields_by_location.contains_key(&key) {
                return self.shared_fields_by_location.get(&key).and_then(named);
            }
        }
        self.shared_fields_by_record
            .get(&(e.record?, allele, feature, terms))
            .and_then(named)
    }

    /// One line for each item that states [`ROWS_ALIKE`], some row of whose key
    /// a comparison read, and whose count is not the rows it found alike:
    /// `rows` holds, by the item's place, the rows read and those alike.
    fn stale_row_counts(&self, rows: &HashMap<usize, (u64, u64)>) -> Vec<String> {
        let mut lines: Vec<(usize, String)> = rows
            .iter()
            .filter_map(|(place, (_, alike))| {
                let item = &self.shared_items[*place];
                (*alike != item.rows_alike).then(|| {
                    (
                        *place,
                        format!(
                            "{}: field {} documents {} rows alike, {} are",
                            item.key, item.field, item.rows_alike, alike
                        ),
                    )
                })
            })
            .collect();
        lines.sort();
        lines.into_iter().map(|(_, line)| line).collect()
    }

    /// Whether any field of the entry is documented to differ.
    pub fn has_divergent_fields(&self, e: &Entry) -> bool {
        self.divergent_fields(e).is_some()
    }

    /// Documented vep-rs consequence sets for an entry, if its key is documented.
    pub fn divergent(&self, e: &Entry) -> Option<&Vec<String>> {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(loc) = &e.location {
            if let Some(v) = self
                .by_location
                .get(&(loc.clone(), allele.clone(), feature.clone()))
            {
                return Some(v);
            }
        }
        if let Some(r) = e.record {
            if let Some(v) = self.by_record.get(&(r, allele, feature)) {
                return Some(v);
            }
        }
        None
    }

    pub fn extra(&self, e: &Entry) -> bool {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(loc) = &e.location {
            if self
                .extra_by_location
                .contains(&(loc.clone(), allele.clone(), feature.clone()))
            {
                return true;
            }
        }
        if let Some(r) = e.record {
            if self.extra_by_record.contains(&(r, allele, feature)) {
                return true;
            }
        }
        false
    }

    /// Whether any documented divergence or vep-rs-only key belongs to the
    /// record, in which case the record's own summary values (VEP's
    /// `most_severe_consequence`) may legitimately differ.
    pub fn record_divergent(&self, record: usize) -> bool {
        self.by_record.keys().any(|(r, _, _)| *r == record)
            || self.extra_by_record.iter().any(|(r, _, _)| *r == record)
    }

    /// Whether the reference dropped the record before annotation, leaving no
    /// record-level line or object to compare vep-rs's against.
    pub fn reference_skipped(&self, record: usize) -> bool {
        self.skipped_records.contains(&record)
    }

    pub fn is_empty(&self) -> bool {
        self.by_location.is_empty() && self.by_record.is_empty()
    }
}

/// The default format's Extra values and the VCF's CSQ values percent-encode
/// `=` and `;`; a documented field value is compared with both undone.
fn decode_extra(s: &str) -> String {
    s.replace("%3D", "=").replace("%3B", ";")
}

/// The default format prints an absent value as `-`; the VCF and JSON formats
/// leave it empty. Keys compare on the `-` form.
fn dash(s: &str) -> String {
    if s.is_empty() {
        "-".to_string()
    } else {
        s.to_string()
    }
}

/// The field value of an entry that lacks the field.
const ABSENT: &str = "<absent>";

/// A documented field value compares across formats on the default format's
/// spelling of absence and of a comma: `-`, an empty VCF value and a key the
/// default or JSON format leaves out are one value, and the `&` the VCF format
/// writes for every comma of a value (`OutputFactory/VCF.pm:401`), between the
/// items of a list and inside a name alike, is the default format's `,`.
fn documented_form(s: &str) -> String {
    if s == ABSENT {
        "-".to_string()
    } else {
        dash(&decode_extra(s).replace('&', ","))
    }
}

/// Whether a printed field value is the documented one: the two
/// [`documented_form`]s equal, or the printed value equal to the documented
/// one with its trailing fractional zeros dropped, the form the JSON format
/// gives a numeric plugin value (`1311092.0000` is written `1311092`, `0.2010`
/// is `0.201`). No other spelling of the number is admitted.
fn is_documented_value(actual: &str, documented: &str) -> bool {
    let (actual, documented) = (documented_form(actual), documented_form(documented));
    actual == documented || actual == without_trailing_zeros(&documented)
}

/// A decimal number (`-12.3400`) with the trailing zeros of its fractional
/// part dropped, and the point when nothing follows it; any other text, a
/// dotted name or an exponent form among them, is returned as it is.
fn without_trailing_zeros(text: &str) -> &str {
    let digits = text.strip_prefix('-').unwrap_or(text);
    let Some((whole, fraction)) = digits.split_once('.') else {
        return text;
    };
    let decimal = |part: &str| !part.is_empty() && part.bytes().all(|b| b.is_ascii_digit());
    if !decimal(whole) || !decimal(fraction) {
        return text;
    }
    text.trim_end_matches('0').trim_end_matches('.')
}

/// Consequence terms normalised to a sorted, comma-joined set for comparison
/// with the manifest's documented sets.
pub fn consequence_set(s: &str) -> String {
    let mut terms: Vec<&str> = s.split(['&', ',']).filter(|t| !t.is_empty()).collect();
    terms.sort_unstable();
    terms.dedup();
    terms.join(",")
}

// --------------------------------------------------------------------------
// Comparison.

const EXAMPLES_PER_GROUP: usize = 4;

/// Whether vep-rs writes VEP's own consequence set on a row of a documented key
/// whose documented `sets` hold that set beside another. Several records share
/// such a key (every `<CNV:TR>` record of one position writes one Location and
/// Allele), the documented divergence is the row of another of them, and this
/// row is compared like any other. A key documented with VEP's own set alone
/// is set aside whole: its rows agree in their terms and differ in a field.
fn agrees_on_a_shared_key(e: &Entry, a: Option<&Entry>, sets: &[String]) -> bool {
    let own = consequence_set(&e.consequence);
    let stated = |held: bool| sets.iter().any(|set| (consequence_set(set) == own) == held);
    a.is_some_and(|a| consequence_set(&a.consequence) == own) && stated(true) && stated(false)
}

/// Compares expected against actual entries and returns a report of every
/// difference, or an empty string when the two agree. Under `ordered` the keys
/// both sides share must also come in the same order.
pub fn compare_entries(
    expected: &[Entry],
    actual: &[Entry],
    documented: &Documented,
    ordered: bool,
) -> String {
    let actual_by_key: HashMap<&str, &Entry> = actual.iter().map(|e| (e.key.as_str(), e)).collect();
    let expected_by_key: HashMap<&str, &Entry> =
        expected.iter().map(|e| (e.key.as_str(), e)).collect();

    let mut missing: Vec<&Entry> = Vec::new();
    let mut unexpected: Vec<&Entry> = Vec::new();
    // field -> expected consequence set -> examples; plus a total per field.
    let mut mismatches: BTreeMap<String, BTreeMap<String, Vec<String>>> = BTreeMap::new();
    let mut mismatch_totals: BTreeMap<String, usize> = BTreeMap::new();
    let mut stale_divergences: Vec<String> = Vec::new();
    // By the place of an item that states `rows_alike`: the rows of its key read
    // and those on which both engines print one value.
    let mut shared_rows: HashMap<usize, (u64, u64)> = HashMap::new();
    let mut compared = 0usize;
    let mut skipped_documented = 0usize;
    let mut annotated_without_reference = 0usize;

    for e in expected {
        let a = actual_by_key.get(e.key.as_str()).copied();
        let documented_sets = documented
            .divergent(e)
            .filter(|sets| !agrees_on_a_shared_key(e, a, sets));
        if let Some(sets) = documented_sets {
            skipped_documented += 1;
            match a {
                None if sets.is_empty() => {}
                None => stale_divergences.push(format!(
                    "{}: documented vep-rs sets {:?}, but vep-rs emits no row",
                    e.key, sets
                )),
                Some(a) => {
                    let got = consequence_set(&a.consequence);
                    // The JSON record object prints no consequence set of its own, so
                    // a record-level row another format documents for the record (the
                    // one row of a `--summary` run) has nothing to observe on it.
                    let record_object =
                        e.key.ends_with("|record") && e.consequence.is_empty() && got.is_empty();
                    if !record_object && !sets.iter().any(|s| consequence_set(s) == got) {
                        stale_divergences.push(format!(
                            "{}: documented vep-rs sets {:?}, vep-rs now emits {:?} (VEP: {:?})",
                            e.key, sets, a.consequence, e.consequence
                        ));
                    }
                }
            }
            continue;
        }
        let Some(a) = a else {
            missing.push(e);
            continue;
        };
        compared += 1;
        let mut names: BTreeSet<&String> = e.fields.keys().collect();
        names.extend(a.fields.keys());
        let record_summary_differs =
            e.key.ends_with("|record") && e.record.is_some_and(|r| documented.record_divergent(r));
        // The default format's `Extra` is the raw text of the fields compared
        // one by one below, so it differs wherever one of them is documented to.
        let extra_differs = documented.has_divergent_fields(e);
        for name in names {
            if record_summary_differs && name == "most_severe_consequence" {
                continue;
            }
            if extra_differs && name == "Extra" {
                continue;
            }
            let ev = e.fields.get(name).map(String::as_str).unwrap_or(ABSENT);
            let av = a.fields.get(name).map(String::as_str).unwrap_or(ABSENT);
            if let Some(documented_value) = documented.divergent_field(e, name) {
                // The field is documented to differ: vep-rs must print the documented
                // value (compared after the Extra/CSQ percent-encoding of `=` and `;`
                // is undone), or the documentation is stale. On a key whose item
                // states `rows_alike`, a row on which both engines print one value
                // is another record's and agrees like any other row; the rows found
                // so are counted against the item's number.
                if let Some(place) = documented.field_shares_its_key(e, name) {
                    let (read, alike) = shared_rows.entry(place).or_insert((0, 0));
                    *read += 1;
                    if ev == av {
                        *alike += 1;
                        continue;
                    }
                }
                if !is_documented_value(av, documented_value) {
                    stale_divergences.push(format!(
                        "{}: field {name} documented as {:?}, vep-rs now prints {:?} (VEP: {:?})",
                        e.key, documented_value, av, ev
                    ));
                }
                continue;
            }
            if ev != av {
                *mismatch_totals.entry(name.clone()).or_default() += 1;
                let group = mismatches
                    .entry(name.clone())
                    .or_default()
                    .entry(e.consequence.clone())
                    .or_default();
                if group.len() < EXAMPLES_PER_GROUP {
                    group.push(format!("{}: expected {:?}, got {:?}", e.key, ev, av));
                }
            }
        }
    }
    stale_divergences.extend(documented.stale_row_counts(&shared_rows));
    for a in actual {
        if expected_by_key.contains_key(a.key.as_str()) || documented.extra(a) {
            continue;
        }
        // Only the record-level line or object is without a counterpart; a
        // consequence entry on a skipped record is a vep-rs-only key, documented
        // or unexpected like any other.
        if a.key.ends_with("|record") && a.record.is_some_and(|r| documented.reference_skipped(r)) {
            annotated_without_reference += 1;
            continue;
        }
        unexpected.push(a);
    }

    // Order: the keys both sides share, in each side's order.
    let shared_expected: Vec<&str> = expected
        .iter()
        .filter(|e| actual_by_key.contains_key(e.key.as_str()) && documented.divergent(e).is_none())
        .map(|e| e.key.as_str())
        .collect();
    let shared_actual: Vec<&str> = actual
        .iter()
        .filter(|a| {
            expected_by_key.contains_key(a.key.as_str()) && documented.divergent(a).is_none()
        })
        .map(|a| a.key.as_str())
        .collect();
    let mut order_note = String::new();
    if ordered && shared_expected != shared_actual {
        let first = shared_expected
            .iter()
            .zip(shared_actual.iter())
            .position(|(e, a)| e != a)
            .unwrap_or(0);
        order_note = format!(
            "ORDER differs: first divergence at shared position {first}: expected {:?}, got {:?}\n",
            shared_expected.get(first),
            shared_actual.get(first)
        );
    }

    if missing.is_empty()
        && unexpected.is_empty()
        && mismatches.is_empty()
        && stale_divergences.is_empty()
        && order_note.is_empty()
    {
        return String::new();
    }

    let mut report = format!(
        "expected entries {}, actual {}, compared {}, documented divergences skipped {}\n",
        expected.len(),
        actual.len(),
        compared,
        skipped_documented
    );
    if annotated_without_reference > 0 {
        report.push_str(&format!(
            "record-level entries on records the reference skipped: {annotated_without_reference}\n"
        ));
    }
    report.push_str(&order_note);
    if !missing.is_empty() {
        report.push_str(&format!("MISSING in vep-rs: {} entries\n", missing.len()));
        let mut by_cons: BTreeMap<&str, Vec<&Entry>> = BTreeMap::new();
        for e in &missing {
            by_cons.entry(e.consequence.as_str()).or_default().push(e);
        }
        for (cons, items) in by_cons {
            report.push_str(&format!("  [{cons}] {}\n", items.len()));
            for e in items.iter().take(EXAMPLES_PER_GROUP) {
                report.push_str(&format!("    {}\n", e.key));
            }
        }
    }
    if !unexpected.is_empty() {
        report.push_str(&format!(
            "UNEXPECTED in vep-rs: {} entries\n",
            unexpected.len()
        ));
        let mut by_cons: BTreeMap<&str, Vec<&Entry>> = BTreeMap::new();
        for e in &unexpected {
            by_cons.entry(e.consequence.as_str()).or_default().push(e);
        }
        for (cons, items) in by_cons {
            report.push_str(&format!("  [{cons}] {}\n", items.len()));
            for e in items.iter().take(EXAMPLES_PER_GROUP) {
                report.push_str(&format!("    {}\n", e.key));
            }
        }
    }
    for (field, groups) in &mismatches {
        report.push_str(&format!(
            "FIELD {field}: {} mismatching entries in {} consequence sets\n",
            mismatch_totals[field],
            groups.len()
        ));
        for (cons, examples) in groups {
            report.push_str(&format!("  [{cons}]\n"));
            for ex in examples {
                report.push_str(&format!("    {ex}\n"));
            }
        }
    }
    if !stale_divergences.is_empty() {
        report.push_str(&format!(
            "DOCUMENTED DIVERGENCES NOT OBSERVED IN THIS RUN: {} (regenerate the manifest with vep-golden-corpus classify)\n",
            stale_divergences.len()
        ));
        for s in stale_divergences.iter().take(20) {
            report.push_str(&format!("    {s}\n"));
        }
    }
    report
}

// --------------------------------------------------------------------------
// Header comparison.

/// Header lines of the default and tab formats that VEP prints with values no
/// second implementation can reproduce: they compare by prefix.
const VOLATILE_PREFIXES: [&str; 3] = [
    "## Output produced at ",
    "## Using cache in ",
    "## VEP command-line: ",
];

/// The engine identity vep-rs writes beside Ensembl VEP's header lines, which
/// the reference lacks: `## vep-rs <version>` in the default and tab headers,
/// `##vep-rs=<version>` in the VCF header and `vep_rs_version` in every JSON
/// object. The binary under test is this crate's version, so the comparators
/// set aside exactly one line of exactly this text and nothing else.
const VEP_RS_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Removes the identity line from a vep-rs header, or describes its absence.
fn take_identity_line(actual: &mut Vec<&str>, identity: &str) -> Result<(), String> {
    let Some(at) = actual.iter().position(|line| *line == identity) else {
        return Err(format!("vep-rs writes no {identity:?} header line"));
    };
    actual.remove(at);
    Ok(())
}

/// The JSON format's identity check: every object carries `vep_rs_version`
/// naming this version; returns the first object that does not, or an empty
/// string.
pub fn compare_json_identity(actual: &str) -> String {
    for (n, line) in actual.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let obj: serde_json::Value = serde_json::from_str(line).unwrap();
        let found = obj.get("vep_rs_version").and_then(|v| v.as_str());
        if found != Some(VEP_RS_VERSION) {
            return format!(
                "JSON object {n}: vep_rs_version is {found:?}, expected {VEP_RS_VERSION:?}"
            );
        }
    }
    String::new()
}

/// Version lines of the Ensembl Perl API components, which vep-rs has no
/// counterpart for, and the `## Connected to <database> on <host>` line of a
/// reference generated with a database connection (a `--check_svs` run), which
/// vep-rs, reading a local cache only, never writes; they may be absent.
fn is_api_component_line(line: &str) -> bool {
    [
        "ensembl",
        "ensembl-compara",
        "ensembl-funcgen",
        "ensembl-io",
        "ensembl-variation",
    ]
    .iter()
    .any(|c| line.starts_with(&format!("## {c} version ")))
        || line.starts_with("## Connected to ")
}

/// The `## Using API version <n>, DB version <m>` line as an offline run writes
/// it: a reference generated with a database connection names the database's
/// version where an offline run, and vep-rs, write `?`.
fn with_unknown_db_version(line: &str) -> String {
    let api_and_db_versions = line
        .strip_prefix("## Using API version ")
        .and_then(|versions| versions.split_once(", DB version "));
    match api_and_db_versions {
        Some((api, _db)) => format!("## Using API version {api}, DB version ?"),
        None => line.to_string(),
    }
}

/// Compares the `#`-prefixed header lines of a default or tab output, the
/// vep-rs identity line set aside; returns a description of the first
/// difference, or an empty string.
pub fn compare_tab_headers(expected: &str, actual: &str) -> String {
    let exp: Vec<&str> = expected.lines().filter(|l| l.starts_with('#')).collect();
    let mut act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
    if let Err(missing) = take_identity_line(&mut act, &format!("## vep-rs {VEP_RS_VERSION}")) {
        return missing;
    }
    let mut ai = 0usize;
    for (i, e) in exp.iter().enumerate() {
        if is_api_component_line(e) {
            if act.get(ai).is_some_and(|a| *a == *e) {
                ai += 1;
            }
            continue;
        }
        let Some(a) = act.get(ai) else {
            return format!("header line {} missing in vep-rs: {e}", i + 1);
        };
        if let Some(prefix) = VOLATILE_PREFIXES.iter().find(|p| e.starts_with(**p)) {
            if !a.starts_with(prefix) {
                return format!(
                    "header line {}: expected prefix {prefix:?}, got {a:?}",
                    i + 1
                );
            }
        } else if description_masked(a) != description_masked(&with_unknown_db_version(e)) {
            return format!("header line {}: expected {e:?}, got {a:?}", i + 1);
        }
        ai += 1;
    }
    if ai != act.len() {
        return format!(
            "vep-rs writes {} extra header line(s), first: {:?}",
            act.len() - ai,
            act[ai]
        );
    }
    String::new()
}

/// A header line with the paths in a plugin's field description masked
/// ([`mask_paths`]): a plugin that names its data file in the description
/// (`## AA : Ancestral allele from <file>`, `##AA=Ancestral allele from <file>`)
/// names the path each engine was given, the reference's inside its container
/// and vep-rs's on the host, so the two compare as the warnings files do; so
/// does a source version line whose value is a file (`## regulatory_gff version
/// <file>`, the first file of `--regulatory_gff`). Every other header line is
/// returned as it is.
fn description_masked(line: &str) -> String {
    let plugin_description = line.contains(" : ")
        || (line.starts_with("## ") && line.contains(" version "))
        || (line.starts_with("##")
            && !line.starts_with("##INFO=")
            && !line.starts_with("##fileformat=")
            && !line.starts_with("##VEP")
            && line.contains('='));
    if plugin_description {
        mask_paths(line)
    } else {
        line.to_string()
    }
}

/// Compares the column descriptions (`## <field> : <text>`, with the two lines
/// that introduce them) and the column line of a default or tab output, the
/// vep-rs identity line set aside, and nothing else of the preamble: a reference
/// generated against the Ensembl database names its connection and API version
/// where a cache run names its cache and source versions. Returns a description
/// of the first difference, or an empty string.
pub fn compare_column_descriptions(expected: &str, actual: &str) -> String {
    fn described(line: &str) -> bool {
        !line.starts_with("##")
            || line == "## Column descriptions:"
            || line == "## Extra column keys:"
            || line.contains(" : ")
    }
    let exp: Vec<&str> = expected
        .lines()
        .filter(|l| l.starts_with('#') && described(l))
        .collect();
    let mut act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
    if let Err(missing) = take_identity_line(&mut act, &format!("## vep-rs {VEP_RS_VERSION}")) {
        return missing;
    }
    act.retain(|l| described(l));
    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        if description_masked(a) != description_masked(e) {
            return format!("header line {}: expected {e:?}, got {a:?}", i + 1);
        }
    }
    if exp.len() != act.len() {
        return format!(
            "header column lines: expected {}, vep-rs writes {}",
            exp.len(),
            act.len()
        );
    }
    String::new()
}

/// Compares VCF meta lines: `##fileformat` and `##INFO` exact, the `##VEP=`
/// line with its `time`, `cache` and `db` tokens removed, so its API module
/// tokens compare exactly and in their order, the command line by prefix, the
/// `##vep-rs=` identity line set aside; the `#CHROM` line exact.
pub fn compare_vcf_headers(expected: &str, actual: &str) -> String {
    compare_vcf_headers_of(expected, actual, false)
}

/// [`compare_vcf_headers`] for a reference Ensembl VEP wrote connected to its
/// database, whose `##VEP=` line names the assembly and gene build alone where a
/// cache run names its cache's source versions too: the reference's tokens must
/// each appear with their value in the actual line, which may carry more.
pub fn compare_database_vcf_headers(expected: &str, actual: &str) -> String {
    compare_vcf_headers_of(expected, actual, true)
}

fn compare_vcf_headers_of(expected: &str, actual: &str, database: bool) -> String {
    /// The `key=value` tokens of a `##VEP=` line, a quoted value kept whole,
    /// without the three whose value is the run's own: `time`, `cache` and `db`.
    fn vep_line_tokens(line: &str) -> Vec<String> {
        // Tokens are `key=value` or `key="quoted value"`; quoted values may hold spaces.
        let mut tokens: Vec<String> = Vec::new();
        let mut rest = line;
        while !rest.is_empty() {
            let rest_trim = rest.trim_start();
            let eq = match rest_trim.find('=') {
                Some(i) => i,
                None => {
                    tokens.push(rest_trim.to_string());
                    break;
                }
            };
            let key = &rest_trim[..eq];
            let after = &rest_trim[eq + 1..];
            let (value, remainder) = if let Some(q) = after.strip_prefix('"') {
                let close = q.find('"').unwrap_or(q.len());
                (&q[..close], &q[(close + 1).min(q.len())..])
            } else {
                let sp = after.find(' ').unwrap_or(after.len());
                (&after[..sp], &after[sp..])
            };
            // `db` names the database a reference generated with a connection
            // read; an offline run writes no such token, nor does vep-rs. A
            // value that is a path (`regulatory_gff`, the first file of
            // `--regulatory_gff`) names each engine's own copy of the file, the
            // reference's inside its container, so it compares masked.
            if !matches!(key, "time" | "cache" | "db") {
                tokens.push(format!("{key}={}", mask_paths(value)));
            }
            rest = remainder;
        }
        tokens
    }
    let exp: Vec<&str> = expected.lines().filter(|l| l.starts_with('#')).collect();
    let mut act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
    if let Err(missing) = take_identity_line(&mut act, &format!("##vep-rs={VEP_RS_VERSION}")) {
        return missing;
    }
    if exp.len() != act.len() {
        return format!(
            "header line count: expected {}, got {}\nexpected:\n{}\nactual:\n{}",
            exp.len(),
            act.len(),
            exp.join("\n"),
            act.join("\n")
        );
    }
    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        if e.starts_with("##VEP=") {
            let (ne, na) = (vep_line_tokens(e), vep_line_tokens(a));
            let agree = if database {
                ne.iter().all(|token| na.contains(token))
            } else {
                ne == na
            };
            if !agree {
                return format!(
                    "header line {}: expected {:?}, got {:?}",
                    i + 1,
                    ne.join(" "),
                    na.join(" ")
                );
            }
        } else if e.starts_with("##VEP-command-line=") {
            if !a.starts_with("##VEP-command-line=") {
                return format!("header line {}: expected a command line, got {a:?}", i + 1);
            }
        } else if description_masked(a) != description_masked(e) {
            return format!("header line {}: expected {e:?}, got {a:?}", i + 1);
        }
    }
    String::new()
}

/// Compares two outputs written under `--no_headers` line by line: neither may
/// carry a `#` line, and the data lines must be the same in the same order.
/// Returns the first difference, or an empty string. A headed output's entry
/// comparison reads every column of every row, so on a run without documented
/// divergences the two comparisons agree; this one needs no column line.
pub fn compare_headerless(expected: &str, actual: &str) -> String {
    for (side, text) in [("reference", expected), ("vep-rs", actual)] {
        if let Some(line) = text.lines().find(|l| l.starts_with('#')) {
            return format!("{side} writes a header line under --no_headers: {line:?}");
        }
    }
    let exp: Vec<&str> = expected.lines().collect();
    let act: Vec<&str> = actual.lines().collect();
    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        if e != a {
            return format!("line {}: expected {e:?}, got {a:?}", i + 1);
        }
    }
    if exp.len() != act.len() {
        return format!("line count: expected {}, got {}", exp.len(), act.len());
    }
    String::new()
}

/// Runs the `vep` binary for the run in one output format and compares the
/// output with the run's expected file for it: the header report (for JSON,
/// which has no header, the identity-key report), then the entry report, each
/// empty when the two agree; a run under `--no_headers` compares its text
/// outputs line by line, since no column line names their columns. `format` is
/// `default`, `tab`, `vcf` or `json`.
pub fn compare_run(run: &Run, format: &str, out_dir: &Path) -> (String, String) {
    let expected = run.expected(match format {
        "default" => "default.txt",
        "tab" => "tab.txt",
        "vcf" => "vcf.vcf",
        "json" => "json.jsonl",
        other => panic!("unknown format {other}"),
    });
    let expected = if run.transcript_source_reference() {
        assert!(
            masks_name(&run.masks, FILE_CACHE_WRITES_NO_SOURCE_FIELD),
            "{}: a reference read from a GTF or GFF names {FILE_CACHE_WRITES_NO_SOURCE_FIELD} in masks",
            run.label()
        );
        without_transcript_source(format, &expected, run.vcf_info_field())
    } else {
        expected
    };
    let expected = masked_reference_text(&run.masks, expected);
    let actual = run_vep(run, format, out_dir, &[]);
    let documented = run.documented(out_dir);
    if format != "json" && run.has_flag("--no_headers") {
        assert!(
            documented.is_empty(),
            "{}: a --no_headers run compares line by line and documents no divergence",
            run.label()
        );
        return (String::new(), compare_headerless(&expected, &actual));
    }
    let inputs = run.input_records();
    let info_field = run.vcf_info_field();
    let parse = |text: &str| -> Vec<Entry> {
        match format {
            "default" | "tab" => parse_default(text),
            "vcf" if run.reads_vcf() => parse_vcf_info(text, &inputs, info_field),
            "vcf" => parse_vcf_info_in_order(text, &inputs, info_field),
            _ => parse_json(text, &inputs),
        }
    };
    let mut expected_entries = parse(&expected);
    let mut actual_entries = parse(&actual);
    mask_entries(
        &run.masks,
        format,
        &mut expected_entries,
        &mut actual_entries,
    );
    if run.unordered {
        number_shared_keys_by_content(&mut expected_entries);
        number_shared_keys_by_content(&mut actual_entries);
    }
    if run.nearest_by_coordinate() {
        let sites = protein_coding_start_sites(run.corpus);
        nearest_as_coordinates(&mut expected_entries, &sites);
        nearest_as_coordinates(&mut actual_entries, &sites);
    }
    let headers = match format {
        "default" | "tab" if run.headers_by_descriptions() => {
            compare_column_descriptions(&expected, &actual)
        }
        "default" | "tab" => compare_tab_headers(&expected, &actual),
        "vcf" if run.headers_by_descriptions() => compare_database_vcf_headers(&expected, &actual),
        "vcf" => compare_vcf_headers(&expected, &actual),
        _ => compare_json_identity(&actual),
    };
    (
        headers,
        compare_entries(
            &expected_entries,
            &actual_entries,
            &documented,
            !run.unordered,
        ),
    )
}

// --------------------------------------------------------------------------
// Warnings and skipped-variants files.

/// One entry of a warnings or skipped-variants file: the input line number the
/// entry names, if any, and its description with every path masked.
pub type WarningEntry = (Option<u64>, String);

/// A description with every space-delimited token that is a path (opening with
/// `/`, `~/`, `./` or `../`) replaced by `<path>`, so a message naming a file
/// compares across machines.
fn mask_paths(description: &str) -> String {
    description
        .split(' ')
        .map(|token| {
            let is_path = token.starts_with('/')
                || token.starts_with("~/")
                || token.starts_with("./")
                || token.starts_with("../");
            if is_path {
                "<path>"
            } else {
                token
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
        .trim()
        .to_string()
}

/// The plugin a run-time plugin warning of `run_plugins` names: `Plugin '<Name>'
/// did not return a hashref, output ignored!` (OutputFactory.pm:2244) or
/// `Plugin '<Name>' went wrong: <reason>` (`:2253`); `None` for any other text.
fn plugin_run_warning_name(description: &str) -> Option<String> {
    let rest = description.strip_prefix("Plugin '")?;
    let (name, tail) = rest.split_once('\'')?;
    (tail == " did not return a hashref, output ignored!" || tail.starts_with(" went wrong: "))
        .then(|| name.to_string())
}

/// The entries of a warnings file or a skipped-variants file, one per line
/// ([`warning_entry`]); a blank line and the header line of a skipped-variants
/// file are left out.
pub fn warning_entries(text: &str) -> BTreeSet<WarningEntry> {
    text.lines()
        .map(str::trim_end)
        .filter(|line| {
            !line.trim().is_empty() && !line.starts_with("[VEP skipped the following variants from")
        })
        .map(warning_entry)
        .collect()
}

/// The entry of one line, in the line shapes of BaseVEP.pm:802-863 and
/// BaseRunner.pm:364-380: `WARNING: line N skipped (<line>): <description>` (or
/// `variant skipped`), `Line N\t<line>\t<description>`, and any other `WARNING:
/// <description>`, whose trailing `line N` is its number. The quoted input line
/// is not compared, and a line of any other shape is an entry of its own.
fn warning_entry(line: &str) -> WarningEntry {
    if let Some(rest) = line.strip_prefix("Line ") {
        let parts: Vec<&str> = rest.splitn(3, '\t').collect();
        if let [number, _quoted, description] = parts[..] {
            return (number.trim().parse().ok(), mask_paths(description));
        }
    }
    let Some(rest) = line.strip_prefix("WARNING:") else {
        return (None, mask_paths(line));
    };
    let rest = rest.trim_start();
    let skipped = if let Some(after_line) = rest.strip_prefix("line ") {
        after_line
            .split_once(" skipped")
            .map(|(digits, tail)| (digits.parse::<u64>().ok(), tail))
    } else {
        rest.strip_prefix("variant skipped")
            .map(|tail| (None, tail))
    };
    if let Some((number, tail)) = skipped {
        let description = match tail.trim_start().strip_prefix('(') {
            Some(quoted) => quoted.split_once("): ").map_or(quoted, |(_, d)| d),
            None => tail.trim_start(),
        };
        return (number, mask_paths(description));
    }
    match rest.rsplit_once(" line ") {
        Some((head, digits))
            if digits.bytes().all(|b| b.is_ascii_digit()) && !digits.is_empty() =>
        {
            (
                digits.parse::<u64>().ok(),
                mask_paths(head.trim_end_matches(" on")),
            )
        }
        _ => (None, mask_paths(rest)),
    }
}

/// Compares two warnings (or skipped-variants) files as entry sets and returns
/// a report of every entry on one side alone, or an empty string.
pub fn compare_warning_sets(expected: &str, actual: &str) -> String {
    let expected = warning_entries(expected);
    let actual = warning_entries(actual);
    let mut report = String::new();
    for (number, description) in expected.difference(&actual) {
        report.push_str(&format!(
            "  MISSING in vep-rs: line {number:?}: {description}\n"
        ));
    }
    for (number, description) in actual.difference(&expected) {
        report.push_str(&format!(
            "  UNEXPECTED in vep-rs: line {number:?}: {description}\n"
        ));
    }
    report
}

// --------------------------------------------------------------------------
// Statistics reports (`--stats_text`).

/// The three sections a statistics report opens with (Stats.pm:806-823), whose
/// rows compare in order; every later section is a chart.
const SUMMARY_TABULAR_SECTIONS: [&str; 3] =
    ["VEP run statistics", "Data version", "General statistics"];

/// The `[VEP run statistics]` rows whose values no second run reproduces
/// (Stats.pm:645-663): the engine's own command line, times and paths. They
/// compare by key alone, as `tools/vep-stats-compare` compares them.
const VOLATILE_SUMMARY_ROWS: [&str; 7] = [
    "Command line options",
    "Start time",
    "End time",
    "Run time",
    "Input file",
    "Input data",
    "Output file",
];

/// The row naming the cache directory and the custom sources (Stats.pm:666-677),
/// whose paths are masked before the values compare.
const SUMMARY_PATH_ROW: &str = "Annotation sources";

/// One `[title]` section of a statistics report: its title and its `(key, value)`
/// rows in file order.
type SummarySection = (String, Vec<(String, String)>);

/// The sections of a statistics report (Stats.pm:806-831): each `[title]` line
/// opens one, and every non-blank line under it is a `key<TAB>value` row (a line
/// without a tab is a key with an empty value). Text before the first section
/// header is reported as a section of its own, titled by the line.
fn summary_sections(text: &str) -> Vec<SummarySection> {
    let mut sections: Vec<SummarySection> = Vec::new();
    for line in text.lines().map(str::trim_end) {
        if line.is_empty() {
            continue;
        }
        if let Some(title) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
            sections.push((title.to_string(), Vec::new()));
            continue;
        }
        if sections.is_empty() {
            sections.push((format!("text before the first section: {line}"), Vec::new()));
            continue;
        }
        let (key, value) = line.split_once('\t').unwrap_or((line, ""));
        sections
            .last_mut()
            .unwrap()
            .1
            .push((key.to_string(), value.to_string()));
    }
    sections
}

/// A value of `[VEP run statistics]` as it compares: a volatile row's value set
/// aside, every path of `Annotation sources` (a token carrying `/`, up to a blank
/// or the `;` joining two sources) replaced by `[PATH]`, any other value as is.
fn masked_summary_value(key: &str, value: &str) -> String {
    if VOLATILE_SUMMARY_ROWS.contains(&key) {
        return "<masked>".to_string();
    }
    if key != SUMMARY_PATH_ROW {
        return value.to_string();
    }
    value
        .split(' ')
        .map(|token| {
            let (body, tail) = match token.strip_suffix(';') {
                Some(body) => (body, ";"),
                None => (token, ""),
            };
            if body.contains('/') {
                format!("[PATH]{tail}")
            } else {
                token.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Compares two `--stats_text` reports as `tools/vep-stats-compare` does, its
/// rules for the two documented differences excepted, and
/// returns a report of every difference, or an empty string: the same
/// sections in the same order; a tabular section row by row, the same keys in
/// the same order and every value equal once [`masked_summary_value`] has set the
/// volatile values aside; a chart section as a set of `(key, count)` pairs, since
/// `sort_keys` (Stats.pm:760-790) breaks ties in hash order; a key printed twice
/// in a section is a difference on its side. The tool holds the
/// `Novel / existing variants` row of a run under a frequency check to
/// `stats_existing_variants_check_frequency`, and a key of
/// `[Consequences (all)]` that names no consequence type, in a run under a
/// `--custom` source keyed `Consequence`, to
/// `stats_consequences_custom_key_named_consequence`; here the row's two values
/// and the key compare as written, like any other.
pub fn compare_summary_sections(expected: &str, actual: &str) -> String {
    let expected = summary_sections(expected);
    let actual = summary_sections(actual);
    let mut report = String::new();
    let titles = |sections: &[SummarySection]| -> Vec<String> {
        sections.iter().map(|(title, _)| title.clone()).collect()
    };
    let (expected_titles, actual_titles) = (titles(&expected), titles(&actual));
    for (side, titles) in [
        ("the reference", &expected_titles),
        ("vep-rs", &actual_titles),
    ] {
        let mut seen = BTreeSet::new();
        for title in titles.iter().filter(|t| !seen.insert(*t)) {
            report.push_str(&format!("  section [{title}] appears twice in {side}\n"));
        }
    }
    for title in expected_titles
        .iter()
        .filter(|t| !actual_titles.contains(t))
    {
        report.push_str(&format!("  MISSING in vep-rs: section [{title}]\n"));
    }
    for title in actual_titles
        .iter()
        .filter(|t| !expected_titles.contains(t))
    {
        report.push_str(&format!("  UNEXPECTED in vep-rs: section [{title}]\n"));
    }
    // The shared titles once each, in each file's order: a repeated title is reported
    // above, not again as an order difference.
    let shared = |titles: &[String], other: &[String]| -> Vec<String> {
        let mut seen = BTreeSet::new();
        titles
            .iter()
            .filter(|t| other.contains(t) && seen.insert((*t).clone()))
            .cloned()
            .collect()
    };
    let shared_expected = shared(&expected_titles, &actual_titles);
    let shared_actual = shared(&actual_titles, &expected_titles);
    if shared_expected != shared_actual {
        report.push_str(&format!(
            "  sections in another order: expected {shared_expected:?}, vep-rs {shared_actual:?}\n"
        ));
    }
    for (title, expected_rows) in &expected {
        let Some((_, actual_rows)) = actual.iter().find(|(t, _)| t == title) else {
            continue;
        };
        // A key printed twice keeps its first value, as `tools/vep-stats-compare` keeps it.
        let mut keyed = |rows: &[(String, String)], side: &str| -> IndexMap<String, String> {
            let mut map = IndexMap::new();
            for (key, value) in rows {
                if map.contains_key(key) {
                    report.push_str(&format!("  [{title}]: {key} printed twice in {side}\n"));
                } else {
                    map.insert(key.clone(), value.clone());
                }
            }
            map
        };
        let expected_keyed = keyed(expected_rows, "the reference");
        let actual_keyed = keyed(actual_rows, "vep-rs");
        if SUMMARY_TABULAR_SECTIONS.contains(&title.as_str()) {
            for (key, value) in expected_rows
                .iter()
                .filter(|(k, _)| !actual_keyed.contains_key(k))
            {
                report.push_str(&format!("  [{title}]: MISSING in vep-rs: {key}\t{value}\n"));
            }
            for (key, value) in actual_rows
                .iter()
                .filter(|(k, _)| !expected_keyed.contains_key(k))
            {
                report.push_str(&format!(
                    "  [{title}]: UNEXPECTED in vep-rs: {key}\t{value}\n"
                ));
            }
            let order = |rows: &IndexMap<String, String>, other: &IndexMap<String, String>| {
                rows.keys()
                    .filter(|key| other.contains_key(*key))
                    .cloned()
                    .collect::<Vec<String>>()
            };
            let expected_order = order(&expected_keyed, &actual_keyed);
            let actual_order = order(&actual_keyed, &expected_keyed);
            if expected_order != actual_order {
                report.push_str(&format!(
                    "  [{title}]: rows in another order: expected {expected_order:?}, vep-rs {actual_order:?}\n"
                ));
            }
            for key in expected_order {
                let expected_value = masked_summary_value(&key, &expected_keyed[&key]);
                let actual_value = masked_summary_value(&key, &actual_keyed[&key]);
                if expected_value != actual_value {
                    report.push_str(&format!(
                        "  [{title}]: {key}: expected {expected_value:?}, vep-rs {actual_value:?}\n"
                    ));
                }
            }
        } else {
            let pairs = |rows: &[(String, String)]| -> BTreeSet<(String, String)> {
                rows.iter().cloned().collect()
            };
            let (expected_pairs, actual_pairs) = (pairs(expected_rows), pairs(actual_rows));
            for (key, count) in expected_pairs.difference(&actual_pairs) {
                report.push_str(&format!("  [{title}]: MISSING in vep-rs: {key}\t{count}\n"));
            }
            for (key, count) in actual_pairs.difference(&expected_pairs) {
                report.push_str(&format!(
                    "  [{title}]: UNEXPECTED in vep-rs: {key}\t{count}\n"
                ));
            }
        }
    }
    report
}

/// A report's path from the `--stats_file` value, or from the default root
/// `<output>_summary` (BaseRunner.pm:254-271): an extension equal to the
/// report's (`txt` or `html`) is kept, the other of the two is swapped for it,
/// and any other extension, or none, gets the report's appended. The extension
/// is whatever follows the last dot, so a dot in a directory name counts when
/// the file name has none, as it does there.
fn stats_file_of(root: &str, extension: &str) -> PathBuf {
    let name = match root.rsplit_once('.') {
        Some((stem, ext)) if !stem.is_empty() && !ext.is_empty() => {
            match ext.to_ascii_lowercase().as_str() {
                ext if ext == extension => root.to_string(),
                "txt" | "html" => format!("{stem}.{extension}"),
                _ => format!("{root}.{extension}"),
            }
        }
        _ => format!("{root}.{extension}"),
    };
    PathBuf::from(name)
}

/// The text report's path from its root ([`stats_file_of`] with `txt`).
pub fn stats_text_file(root: &str) -> PathBuf {
    stats_file_of(root, "txt")
}

/// The HTML report's path from its root ([`stats_file_of`] with `html`).
pub fn stats_html_file(root: &str) -> PathBuf {
    stats_file_of(root, "html")
}

/// The root the statistics reports of a manifest run are named from: the
/// `--stats_file` (or `--sf`) value, else `<output>_summary`, where the output is
/// the `-o` value as spelled for a run writing to standard output
/// (`STDOUT_summary` in the working directory, since BaseRunner.pm:254 reads
/// `output_file` as given) and the `{out}` path otherwise.
fn stats_root(argv: &[String], out: &Path, out_dir: &Path) -> String {
    let named = ["--stats_file", "--sf"]
        .iter()
        .find_map(|flag| value_after(argv, flag));
    match named {
        Some(root) => root.to_string(),
        None => {
            let spelled = ["-o", "--output_file"]
                .iter()
                .find_map(|flag| value_after(argv, flag));
            match spelled {
                Some(name) if name.eq_ignore_ascii_case("stdout") => out_dir
                    .join(format!("{name}_summary"))
                    .display()
                    .to_string(),
                _ => format!("{}_summary", out.display()),
            }
        }
    }
}

/// The text statistics report a manifest run writes ([`stats_root`] through
/// [`stats_text_file`]).
pub fn stats_report_path(argv: &[String], out: &Path, out_dir: &Path) -> PathBuf {
    stats_text_file(&stats_root(argv, out, out_dir))
}

/// The HTML statistics report a manifest run writes ([`stats_root`] through
/// [`stats_html_file`]).
pub fn stats_html_report_path(argv: &[String], out: &Path, out_dir: &Path) -> PathBuf {
    stats_html_file(&stats_root(argv, out, out_dir))
}

/// The element ids every HTML statistics report carries (Stats.pm:875-896 and
/// the four charts of every run, :898-909, each with its chart element and its
/// table).
pub const HTML_REPORT_IDS: [&str; 16] = [
    "masthead",
    "run_stats",
    "data_version",
    "gen_stats",
    "var_class",
    "var_class_pie",
    "var_class_table",
    "var_cons",
    "var_cons_pie",
    "var_cons_table",
    "consequences",
    "consequences_pie",
    "consequences_table",
    "coding",
    "coding_pie",
    "coding_table",
];

/// What `html` lacks of a self-contained HTML statistics report: the page
/// opening, the title, every id of [`HTML_REPORT_IDS`]; and what it must not
/// hold: a script or a remote (`http`) resource. Empty when it is one.
pub fn check_html_report(html: &str) -> String {
    let mut report = String::new();
    if !html.starts_with("<html>") {
        report.push_str("  the page does not open with <html>\n");
    }
    if !html.contains("<title>VEP summary</title>") {
        report.push_str("  no <title>VEP summary</title>\n");
    }
    for id in HTML_REPORT_IDS {
        if !html.contains(&format!(" id=\"{id}\"")) {
            report.push_str(&format!("  no element with id {id}\n"));
        }
    }
    if html.contains("<script") {
        report.push_str("  the page holds a <script element\n");
    }
    if html.contains("http://") || html.contains("https://") {
        report.push_str("  the page names a remote resource\n");
    }
    report
}

/// The warnings Ensembl VEP wrote for the run's records, from the manifest: one
/// entry per record carrying `reference_warning`, at the record's 1-based
/// ordinal, which is the line number Ensembl VEP reports (header lines are not
/// counted). Those it wrote about no one record are
/// [`run_reference_warnings`]'s.
pub fn reference_warnings(manifest: &serde_json::Value) -> BTreeSet<WarningEntry> {
    manifest["records"]
        .as_array()
        .into_iter()
        .flatten()
        .enumerate()
        .filter_map(|(i, record)| {
            record["reference_warning"]
                .as_str()
                .map(|w| (Some(i as u64 + 1), mask_paths(w)))
        })
        .collect()
}

/// The warnings Ensembl VEP wrote for a run about no one record, each text of
/// its `reference_warnings` list as an unnumbered entry: the flag set's own list
/// when its entry carries one, else the manifest's ([`runs`] holds both to
/// lists of non-empty texts). A tandem repeat or tandem duplication record
/// read without a FASTA draws one such warning (Parser/VCF.pm:521-524).
pub fn run_reference_warnings(run: &Run) -> BTreeSet<WarningEntry> {
    let list = run
        .entry
        .get("reference_warnings")
        .unwrap_or(&run.corpus.manifest["reference_warnings"]);
    strings(list)
        .iter()
        .map(|text| (None, mask_paths(text)))
        .collect()
}

/// The warnings file a run leaves beside its output: `<output>_warnings.txt`.
fn warnings_file_beside(output: &Path) -> PathBuf {
    PathBuf::from(format!("{}_warnings.txt", output.display()))
}

/// Compares the warnings file the run's default-format output left beside it
/// (`run_vep` in the `default` format must have run first) with the reference's
/// warnings, which are the manifest's `reference_warning` entries, the run's
/// `reference_warnings` list and, when the run's expected directory holds
/// `default.txt_warnings.txt`, that file, whose entries must agree with the
/// manifest's. On a corpus whose cache holds no `variations` directory the line
/// vep-rs writes about that ([`VARIATION_NOTICE`]) is no part of the comparison.
/// Returns a report, or an empty string.
pub fn compare_warnings(run: &Run, out_dir: &Path) -> String {
    let path = warnings_file_beside(&output_path(run, "default", out_dir));
    let mut actual = warning_entries(&fs::read_to_string(&path).unwrap_or_default());
    if !run.corpus.cache_dir().join("variations").is_dir() {
        actual.retain(|(_, description)| !description.contains(VARIATION_NOTICE));
    }
    let mut expected = reference_warnings(&run.corpus.manifest);
    expected.extend(run_reference_warnings(run));
    let mut report = String::new();
    if run.masks_regulatory_rows() {
        let before = actual.len();
        actual.retain(|(_, description)| !description.contains(REGULATORY_NOTICE));
        if actual.len() == before {
            report.push_str(&format!(
                "{}: the regulatory notice ({REGULATORY_NOTICE}) is missing from {}\n",
                run.label(),
                path.display()
            ));
        }
    }
    let stored = run
        .corpus
        .dir
        .join(&run.expected_dir)
        .join("default.txt_warnings.txt");
    if stored.is_file() {
        let stored_entries = warning_entries(&fs::read_to_string(&stored).unwrap());
        if stored_entries != expected {
            report.push_str(&format!(
                "{}: {} disagrees with the manifest's reference_warning entries\n",
                run.label(),
                stored.display()
            ));
        }
    }
    if run.masks_unloaded_module_plugin() {
        let unloaded: BTreeSet<String> = expected
            .iter()
            .filter_map(|(_, description)| plugin_run_warning_name(description))
            .collect();
        if unloaded.is_empty() {
            report.push_str(&format!(
                "{}: {UNLOADED_MODULE_PLUGIN} is named but the reference warns about no plugin\n",
                run.label()
            ));
        }
        expected.retain(|(_, description)| plugin_run_warning_name(description).is_none());
        for name in unloaded {
            let load_failure = format!("Failed to load plugin {name}: ");
            let before = actual.len();
            actual.retain(|(_, description)| !description.starts_with(&load_failure));
            let removed = before - actual.len();
            if removed != 1 {
                report.push_str(&format!(
                    "{}: {removed} lines of {} report that the plugin module {name} failed to \
                     load, where one was expected\n",
                    run.label(),
                    path.display()
                ));
            }
        }
    }
    for (number, description) in expected.difference(&actual) {
        report.push_str(&format!(
            "{}: reference warning missing from {}: line {number:?}: {description}\n",
            run.label(),
            path.display()
        ));
    }
    for (number, description) in actual.difference(&expected) {
        report.push_str(&format!(
            "{}: warning not in the reference: line {number:?}: {description}\n",
            run.label()
        ));
    }
    report
}

/// The value following `flag` in an argument list.
fn value_after<'a>(argv: &'a [String], flag: &str) -> Option<&'a str> {
    argv.iter()
        .position(|t| t == flag)
        .and_then(|i| argv.get(i + 1))
        .map(String::as_str)
}

// --------------------------------------------------------------------------
// Manifest runs: command lines with process-level expectations.

/// What a manifest run expects of the process.
pub struct Expect {
    pub exit: i32,
    /// An expected output file (under `expected/`) the run's output is compared
    /// with as the default run's is: the file the run wrote at `{out}`,
    /// decompressed when the command line names `--compress_output`, or its
    /// standard output when the command line names no `{out}`.
    pub output: Option<String>,
    /// The format `output` is compared in: `default`, `tab`, `vcf` or `json`.
    pub format: String,
    pub stdout_starts_with: Option<String>,
    pub stdout_contains: Option<String>,
    pub stderr_contains: Option<String>,
    /// An expected warnings file (under `expected/`) compared with the run's
    /// warnings file, the one `--warning_file` names or `{out}_warnings.txt`, as
    /// the set of (line number, description) its lines carry, paths masked.
    pub warnings_file: Option<String>,
    /// An expected skipped-variants file (under `expected/`) compared the same
    /// way with the file the run's `--skipped_variants_file` names.
    pub skipped_file: Option<String>,
    /// An expected statistics report (under `expected/`, the `--stats_text`
    /// output of the reference run) compared with the run's, the file its
    /// `--stats_file` names or `{out}_summary.txt`, section by section
    /// ([`compare_summary_sections`]).
    pub stats_file: Option<String>,
    /// Whether the run writes the HTML statistics report, the file its
    /// `--stats_file` names through the HTML name or `{out}_summary.html`
    /// (`STDOUT_summary.html` as spelled, in the working directory), checked as
    /// [`check_html_report`] checks one.
    pub stats_html: bool,
}

/// One entry of a manifest's `runs` list.
pub struct ManifestRun {
    pub name: String,
    /// The command line after the binary name, before placeholder expansion.
    pub argv: Vec<String>,
    pub env: BTreeMap<String, String>,
    /// A file in the corpus whose bytes are piped on standard input.
    pub stdin: Option<String>,
    pub expect: Expect,
}

/// The keys an `expect` block may carry.
const EXPECT_KEYS: [&str; 10] = [
    "exit",
    "output",
    "format",
    "stdout_starts_with",
    "stdout_contains",
    "stderr_contains",
    "warnings_file",
    "skipped_file",
    "stats_file",
    "stats_html",
];
const RUN_KEYS: [&str; 5] = ["name", "argv", "env", "stdin", "expect"];

/// The manifest's `runs`, each held to one shape: a plain unique name, `argv` a
/// list of strings, `env` an object of strings, `stdin` a file inside the
/// corpus, and an `expect` block with an integer `exit` and at least one of
/// `output`, `warnings_file`, `skipped_file` and `stats_file` (files under
/// `expected/`), `stdout_starts_with`, `stdout_contains` and `stderr_contains`,
/// each a non-empty string, or `stats_html` true, so an entry that expects
/// nothing cannot certify the rows the manifest covers; a `skipped_file` needs
/// `--skipped_variants_file` in the argv, a `stats_file` needs `--stats_text`
/// there and no `--no_stats` or `{golden}`, which carries it ([`base_args`]),
/// since a run under `--no_stats` writes no report, and `stats_html` needs no
/// `--no_stats` or `{golden}` and `--stats_html` beside any `--stats_text`,
/// which alone turns the HTML report off. An unknown key anywhere is refused
/// rather than ignored.
pub fn manifest_runs(corpus: &Corpus) -> Vec<ManifestRun> {
    let label = format!("{}/{}", corpus.release, corpus.name);
    let mut out: Vec<ManifestRun> = Vec::new();
    for run in corpus.manifest["runs"].as_array().into_iter().flatten() {
        let object = run
            .as_object()
            .unwrap_or_else(|| panic!("{label}: a run is an object"));
        for key in object.keys() {
            assert!(
                RUN_KEYS.contains(&key.as_str()),
                "{label}: run key {key:?} is not one of {RUN_KEYS:?}"
            );
        }
        let name = run["name"]
            .as_str()
            .unwrap_or_else(|| panic!("{label}: a run has a name"))
            .to_string();
        assert!(
            plain_name(&name),
            "{label}: run name {name:?} is not a plain word other than default"
        );
        assert!(
            out.iter().all(|r| r.name != name),
            "{label}: run {name} is named twice"
        );
        assert!(
            run["argv"].is_array(),
            "{label} run {name}: argv is a list of strings"
        );
        let argv = strings(&run["argv"]);
        let mut env = BTreeMap::new();
        if let Some(map) = run.get("env") {
            let map = map
                .as_object()
                .unwrap_or_else(|| panic!("{label} run {name}: env is an object"));
            for (k, v) in map {
                let v = v
                    .as_str()
                    .unwrap_or_else(|| panic!("{label} run {name}: env {k} is a string"));
                env.insert(k.clone(), v.to_string());
            }
        }
        let stdin = run.get("stdin").map(|s| {
            let s = s
                .as_str()
                .unwrap_or_else(|| panic!("{label} run {name}: stdin is a file name"));
            assert!(
                inside_corpus(s),
                "{label} run {name}: stdin {s:?} is not a path inside the corpus"
            );
            s.to_string()
        });
        let expect = run["expect"]
            .as_object()
            .unwrap_or_else(|| panic!("{label} run {name}: expect is an object"));
        for key in expect.keys() {
            assert!(
                EXPECT_KEYS.contains(&key.as_str()),
                "{label} run {name}: expect key {key:?} is not one of {EXPECT_KEYS:?}"
            );
        }
        let exit = expect
            .get("exit")
            .and_then(|e| e.as_i64())
            .unwrap_or_else(|| panic!("{label} run {name}: expect.exit is an integer"));
        let text = |key: &str| -> Option<String> {
            match expect.get(key) {
                None | Some(serde_json::Value::Null) => None,
                Some(v) => {
                    let v = v
                        .as_str()
                        .unwrap_or_else(|| panic!("{label} run {name}: expect.{key} is a string"));
                    assert!(!v.is_empty(), "{label} run {name}: expect.{key} is empty");
                    Some(v.to_string())
                }
            }
        };
        let expected_file = |key: &str| -> Option<String> {
            let file = text(key);
            if let Some(file) = &file {
                assert!(
                    inside_corpus(file),
                    "{label} run {name}: {key} {file:?} is not a path inside the corpus"
                );
            }
            file
        };
        let output = expected_file("output");
        let warnings_file = expected_file("warnings_file");
        let skipped_file = expected_file("skipped_file");
        assert!(
            skipped_file.is_none() || argv.iter().any(|t| t == "--skipped_variants_file"),
            "{label} run {name}: skipped_file needs --skipped_variants_file in argv"
        );
        let stats_file = expected_file("stats_file");
        assert!(
            stats_file.is_none()
                || (argv.iter().any(|t| t == "--stats_text")
                    && !argv.iter().any(|t| t == "--no_stats" || t == "{golden}")),
            "{label} run {name}: stats_file needs --stats_text and no --no_stats or {{golden}} in argv"
        );
        let stats_html = match expect.get("stats_html") {
            None | Some(serde_json::Value::Null) => false,
            Some(v) => v
                .as_bool()
                .unwrap_or_else(|| panic!("{label} run {name}: expect.stats_html is a boolean")),
        };
        assert!(
            !stats_html
                || (!argv.iter().any(|t| t == "--no_stats" || t == "{golden}")
                    && (argv.iter().any(|t| t == "--stats_html")
                        || !argv.iter().any(|t| t == "--stats_text"))),
            "{label} run {name}: stats_html needs no --no_stats or {{golden}} in argv, and --stats_html beside --stats_text"
        );
        let format = text("format").unwrap_or_else(|| "default".to_string());
        assert!(
            ALL_FORMATS.contains(&format.as_str()),
            "{label} run {name}: unknown format {format}"
        );
        let stdout_starts_with = text("stdout_starts_with");
        let stdout_contains = text("stdout_contains");
        let stderr_contains = text("stderr_contains");
        assert!(
            output.is_some()
                || stdout_starts_with.is_some()
                || stdout_contains.is_some()
                || stderr_contains.is_some()
                || warnings_file.is_some()
                || skipped_file.is_some()
                || stats_file.is_some()
                || stats_html,
            "{label} run {name}: expects nothing beyond the exit status"
        );
        out.push(ManifestRun {
            name,
            argv,
            env,
            stdin,
            expect: Expect {
                exit: exit as i32,
                output,
                format,
                stdout_starts_with,
                stdout_contains,
                stderr_contains,
                warnings_file,
                skipped_file,
                stats_file,
                stats_html,
            },
        });
    }
    out
}

/// Expands a run's `argv`: the token `{golden}` becomes the default run's
/// arguments (`-i <input>` and [`base_args`]), `{flags}` the manifest's own
/// flags; inside a token `{input}` is the corpus input path, `{cache}` the
/// corpus cache directory itself, `{cache_root}` the root of the Ensembl layout,
/// `{corpus}` the corpus directory, `{plugins}` its `plugin_data` directory,
/// `{fasta_dir}` the [`staged_fasta_dir`] under `out_dir`, `{generated}` the
/// [`generated_dir`] under `out_dir`, `{assembly}` the assembly and `{out}` the
/// output path the caller chose.
pub fn expand_argv(corpus: &Corpus, argv: &[String], out: &Path, out_dir: &Path) -> Vec<String> {
    let default = default_run(corpus);
    let input = corpus.dir.join("variants.vcf").display().to_string();
    let cache = corpus.run_cache_dir(out_dir).display().to_string();
    let cache_root = corpus.cache_root().display().to_string();
    let out = out.display().to_string();
    let argv = expand_generated(corpus, argv, out_dir);
    let mut expanded = Vec::new();
    for token in &argv {
        match token.as_str() {
            "{golden}" => {
                expanded.push("-i".to_string());
                expanded.push(input.clone());
                expanded.extend(base_args(&default, out_dir));
            }
            "{flags}" => expanded.extend(default.flags.iter().cloned()),
            _ => {
                let token = if token.contains("{fasta_dir}") {
                    token.replace(
                        "{fasta_dir}",
                        &staged_fasta_dir(corpus, out_dir).display().to_string(),
                    )
                } else {
                    token.clone()
                };
                expanded.push(expand_corpus_placeholders(
                    corpus,
                    &token
                        .replace("{input}", &input)
                        .replace("{cache_root}", &cache_root)
                        .replace("{cache}", &cache)
                        .replace("{assembly}", &corpus.assembly)
                        .replace("{out}", &out),
                ))
            }
        }
    }
    expanded
}

/// Executes a manifest run in `out_dir` and returns a report of every
/// expectation it missed, or an empty string.
pub fn execute_manifest_run(corpus: &Corpus, run: &ManifestRun, out_dir: &Path) -> String {
    let out = out_dir.join(format!("{}-run-{}.out", corpus.name, run.name));
    let argv = expand_argv(corpus, &run.argv, &out, out_dir);
    let mut cmd = isolated_vep(out_dir);
    cmd.args(&argv)
        .envs(&run.env)
        .current_dir(out_dir)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let output = match &run.stdin {
        Some(file) => {
            let input = fs::read(corpus.dir.join(file)).unwrap();
            output_with_stdin(&mut cmd, input)
        }
        None => {
            cmd.stdin(Stdio::null());
            cmd.output().expect("vep binary runs")
        }
    };
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    let mut report = String::new();
    let label = format!("{}/{} run {}", corpus.release, corpus.name, run.name);
    if output.status.code() != Some(run.expect.exit) {
        report.push_str(&format!(
            "{label}: exit {:?}, expected {}\nstderr:\n{stderr}",
            output.status.code(),
            run.expect.exit
        ));
    }
    if let Some(prefix) = &run.expect.stdout_starts_with {
        if !stdout.starts_with(prefix) {
            report.push_str(&format!(
                "{label}: stdout does not start with {prefix:?}:\n{}\n",
                stdout.lines().take(3).collect::<Vec<_>>().join("\n")
            ));
        }
    }
    if let Some(needle) = &run.expect.stdout_contains {
        if !stdout.contains(needle) {
            report.push_str(&format!(
                "{label}: stdout does not contain {needle:?}:\n{}\n",
                stdout.lines().take(3).collect::<Vec<_>>().join("\n")
            ));
        }
    }
    if let Some(needle) = &run.expect.stderr_contains {
        if !stderr.contains(needle) {
            report.push_str(&format!(
                "{label}: stderr does not contain {needle:?}:\n{stderr}"
            ));
        }
    }
    if let Some(expected_name) = &run.expect.warnings_file {
        let path = value_after(&argv, "--warning_file")
            .map(PathBuf::from)
            .unwrap_or_else(|| warnings_file_beside(&out));
        let actual = fs::read_to_string(&path).unwrap_or_default();
        let differences = compare_warning_sets(&read_expected(&corpus.dir, expected_name), &actual);
        if !differences.is_empty() {
            report.push_str(&format!(
                "{label}: warnings file {} differs from {expected_name}\n{differences}",
                path.display()
            ));
        }
    }
    if let Some(expected_name) = &run.expect.skipped_file {
        let path = value_after(&argv, "--skipped_variants_file")
            .map(PathBuf::from)
            .expect("manifest_runs requires --skipped_variants_file in the argv");
        let actual = fs::read_to_string(&path).unwrap_or_else(|e| {
            report.push_str(&format!(
                "{label}: no skipped-variants file at {}: {e}\n",
                path.display()
            ));
            String::new()
        });
        let differences = compare_warning_sets(&read_expected(&corpus.dir, expected_name), &actual);
        if !differences.is_empty() {
            report.push_str(&format!(
                "{label}: skipped-variants file {} differs from {expected_name}\n{differences}",
                path.display()
            ));
        }
    }
    if let Some(expected_name) = &run.expect.stats_file {
        let path = stats_report_path(&argv, &out, out_dir);
        let actual = fs::read_to_string(&path).unwrap_or_else(|e| {
            report.push_str(&format!(
                "{label}: no statistics report at {}: {e}\n",
                path.display()
            ));
            String::new()
        });
        let differences =
            compare_summary_sections(&read_expected(&corpus.dir, expected_name), &actual);
        if !differences.is_empty() {
            report.push_str(&format!(
                "{label}: statistics report {} differs from {expected_name}\n{differences}",
                path.display()
            ));
        }
    }
    if run.expect.stats_html {
        let path = stats_html_report_path(&argv, &out, out_dir);
        match fs::read_to_string(&path) {
            Ok(html) => {
                let lacks = check_html_report(&html);
                if !lacks.is_empty() {
                    report.push_str(&format!(
                        "{label}: HTML statistics report {} is not self-contained\n{lacks}",
                        path.display()
                    ));
                }
            }
            Err(e) => report.push_str(&format!(
                "{label}: no HTML statistics report at {}: {e}\n",
                path.display()
            )),
        }
    }
    if let Some(expected_name) = &run.expect.output {
        let actual = if run.argv.iter().any(|t| t.contains("{out}")) {
            match fs::read(&out) {
                Err(e) => {
                    report.push_str(&format!("{label}: no output at {}: {e}\n", out.display()));
                    String::new()
                }
                Ok(bytes)
                    if argv.iter().any(|t| {
                        t == "--compress_output" || t.starts_with("--compress_output=")
                    }) =>
                {
                    let mut text = String::new();
                    if !bytes.starts_with(&[0x1f, 0x8b]) {
                        report.push_str(&format!(
                            "{label}: --compress_output wrote no gzip stream\n"
                        ));
                    } else if let Err(e) =
                        flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut text)
                    {
                        report.push_str(&format!(
                            "{label}: {} does not decompress: {e}\n",
                            out.display()
                        ));
                    }
                    text
                }
                Ok(bytes) => String::from_utf8(bytes).unwrap_or_else(|e| {
                    report.push_str(&format!("{label}: {} is not UTF-8: {e}\n", out.display()));
                    String::new()
                }),
            }
        } else {
            stdout.clone()
        };
        let expected = read_expected(&corpus.dir, expected_name);
        let documented = Documented::for_run(
            &corpus.manifest,
            has_reference_sequence(&argv, &corpus.run_cache_dir(out_dir), corpus),
        );
        let inputs = input_records(corpus);
        let format = run.expect.format.as_str();
        let parse = |text: &str| -> Vec<Entry> {
            match format {
                "default" | "tab" => parse_default(text),
                "vcf" => parse_vcf(text, &inputs),
                _ => parse_json(text, &inputs),
            }
        };
        let mut expected_entries = parse(&expected);
        let mut actual_entries = parse(&actual);
        mask_entries(
            &corpus.masks(),
            format,
            &mut expected_entries,
            &mut actual_entries,
        );
        let header = match format {
            "default" | "tab" => compare_tab_headers(&expected, &actual),
            "vcf" => compare_vcf_headers(&expected, &actual),
            _ => compare_json_identity(&actual),
        };
        let body = compare_entries(&expected_entries, &actual_entries, &documented, true);
        if !header.is_empty() || !body.is_empty() {
            report.push_str(&format!(
                "{label}: output differs from {expected_name}\n{header}{body}"
            ));
        }
    }
    report
}

/// Writes `head` and then `records`, one line each, bgzip-compressed at
/// `<dir>/<name>` with no flush between lines, and a tabix index describing
/// `header` at `<dir>/<name>.tbi`. The index holds each record under the name
/// in its first column, at the one-based span `span` reads from its columns;
/// the `head` lines belong to no chunk.
pub fn indexed_lines(
    dir: &Path,
    name: &str,
    header: noodles::csi::binning_index::index::Header,
    head: &[&str],
    records: &[&str],
    span: impl Fn(&[&str]) -> (usize, usize),
) -> PathBuf {
    use noodles::core::Position;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::{bgzf, tabix};

    let path = dir.join(name);
    let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(header);
    for line in head {
        writeln!(writer, "{line}").unwrap();
    }
    for record in records {
        let start = writer.virtual_position();
        writeln!(writer, "{record}").unwrap();
        let end = writer.virtual_position();
        let fields: Vec<&str> = record.split('\t').collect();
        let (first, last) = span(&fields);
        indexer
            .add_record(
                fields[0],
                Position::try_from(first).unwrap(),
                Position::try_from(last).unwrap(),
                Chunk::new(start, end),
            )
            .unwrap();
    }
    writer.finish().unwrap();
    let mut index_path = path.as_os_str().to_owned();
    index_path.push(".tbi");
    tabix::fs::write(index_path, &indexer.build()).unwrap();
    path
}

/// [`indexed_lines`] with no line ahead of the records; `columns` names the
/// zero-based start and end columns and the offset that makes the start
/// one-based.
pub fn indexed_file(
    dir: &Path,
    name: &str,
    header: noodles::csi::binning_index::index::Header,
    columns: (usize, usize, usize),
    records: &[&str],
) -> PathBuf {
    let (start_column, end_column, start_offset) = columns;
    let position = |column: &str| column.parse::<usize>().unwrap();
    indexed_lines(dir, name, header, &[], records, |fields| {
        (
            position(fields[start_column]) + start_offset,
            position(fields[end_column]),
        )
    })
}

/// A bgzip-compressed BED of `records` (tab-separated, zero-based starts) under
/// `dir` with its tabix index beside it, as `bgzip` and `tabix -p bed` write
/// them.
pub fn indexed_bed(dir: &Path, name: &str, records: &[&str]) -> PathBuf {
    use noodles::csi::binning_index::index::header;
    indexed_file(
        dir,
        name,
        header::Builder::bed().build(),
        (1, 2, 1),
        records,
    )
}

/// A bgzip-compressed GFF3 of `records` under `dir` with its tabix index
/// beside it, as `bgzip` and `tabix -p gff` write them.
pub fn indexed_gff(dir: &Path, name: &str, records: &[&str]) -> PathBuf {
    use noodles::csi::binning_index::index::header;
    indexed_file(
        dir,
        name,
        header::Builder::gff().build(),
        (3, 4, 0),
        records,
    )
}

/// A bgzipped VCF at `<dir>/<name>` with its tabix index beside it: the
/// `header` lines, the `#CHROM` line, then the one `record`, indexed over the
/// bases of its `REF`.
pub fn indexed_vcf(dir: &Path, name: &str, header: &[&str], record: &str) -> PathBuf {
    use noodles::csi::binning_index::index::header::Builder;
    let mut head = header.to_vec();
    head.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    indexed_lines(
        dir,
        name,
        Builder::vcf().build(),
        &head,
        &[record],
        |fields| {
            let pos: usize = fields[1].parse().unwrap();
            (pos, pos + fields[3].len() - 1)
        },
    )
}

/// Sum of file sizes under a directory.
pub fn dir_size(dir: &Path) -> u64 {
    let mut total = 0;
    for entry in fs::read_dir(dir).unwrap().flatten() {
        let p = entry.path();
        total += if p.is_dir() {
            dir_size(&p)
        } else {
            p.metadata().unwrap().len()
        };
    }
    total
}
