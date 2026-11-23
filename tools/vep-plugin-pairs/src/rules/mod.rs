// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class rules: what a documented class accepts of a reference value and an engine value.
//!
//! A class is consulted with the field, the two values, the reference's row and the
//! engine's, the rows of the reference's input record and what the run handed the class's
//! plugin (`Consult`). `BuiltIn` holds the tool's own classes: it tests that the class reads
//! the field (`fields`) and hands the consult to the class's rule.

pub mod cache;
pub mod coding_end;
/// What stands at the two ends of a cache transcript's coding sequence, and the coding
/// terms vep-rs owes on a row there: the sources of `tools/vep-formats-compare`, compiled
/// here so that both tools read a record and derive a row by one rule.
#[path = "../../../vep-formats-compare/src/coding_ends.rs"]
pub mod coding_ends;
pub mod downstream;
pub mod entry;
pub mod exon_spans;
pub mod fields;
pub mod gwas;
pub mod hgvs_offsets;
pub mod loftee;
pub mod maxentscan;
pub mod mitochondrial;
pub mod nearest_exon;
pub mod record;
pub mod revel;
pub mod shared;
pub mod sibling;
pub mod skipped_record;
pub mod spliceai;
pub mod tabix;
pub mod term_pair;
pub mod utrannotator;

pub use entry::{
    blank, consequence_set, decode_extra, field_named, field_of, parse_default, uploaded_alleles,
    Entry, Fields, Row, ABSENT,
};
pub use record::RecordRows;

use crate::pytext::{self, Json};
use std::collections::{BTreeMap, BTreeSet};
use std::io;
use std::sync::Arc;

/// A refusal: the message of one error line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Refusal(pub String);

pub type Res<T> = Result<T, Refusal>;

impl From<io::Error> for Refusal {
    fn from(error: io::Error) -> Refusal {
        Refusal(error.to_string())
    }
}

/// What a run hands the classes of a plugin; `None` where the run gave none.
#[derive(Default)]
pub struct Reads {
    pub gwas_pairs: Option<gwas::GwasCatalogPairs>,
    pub revel_table: Option<revel::RevelTable>,
    pub utrannotator_sweep: Option<utrannotator::UtrannotatorSeedSweep>,
    pub downstream_sequences: Option<Arc<downstream::DownstreamSequences>>,
    pub splice_windows: Option<maxentscan::SpliceWindows>,
    pub spliceai_table: Option<spliceai::SpliceAiTable>,
    pub exon_spans: Option<exon_spans::ExonSpans>,
}

impl Reads {
    /// Reads the shards of a JSON cache on `threads` threads from the next chromosome a
    /// class asks for. What the readers hold and refuse is the same on any number of
    /// threads.
    pub fn read_on(&self, threads: usize) {
        if let Some(sequences) = &self.downstream_sequences {
            sequences.read_on(threads);
        }
        if let Some(windows) = &self.splice_windows {
            windows.read_on(threads);
        }
        if let Some(table) = &self.spliceai_table {
            table.read_on(threads);
        }
        if let Some(spans) = &self.exon_spans {
            spans.read_on(threads);
        }
    }
}

/// One question to a class: whether it accepts the differing field `name` of the
/// reference's row `e`.
pub struct Consult<'a> {
    pub klass: &'a str,
    /// The reference's row.
    pub e: &'a Row,
    /// The field.
    pub name: &'a str,
    /// The reference's value, `ABSENT` where it wrote none.
    pub expected: &'a str,
    /// The engine's value, `ABSENT` where it wrote none.
    pub actual: &'a str,
    /// The engine's row.
    pub a: Option<&'a Row>,
    /// The reference's rows of `e`'s own input record, which the caller vouches are one
    /// record's.
    pub record_rows: Option<&'a [Row]>,
    /// The engine's rows of the same input record.
    pub candidate_rows: Option<&'a [Row]>,
    /// The classes applied with `klass`; empty where none is handed.
    pub declared: &'a BTreeSet<String>,
    /// The assembly the run's command line names; `None` where the run gave no command
    /// line or its words name none.
    pub assembly: Option<&'a str>,
    /// The run's source of the field: the plugin among those the run loaded that the
    /// registry names for it, else the custom source whose prefix the field starts with.
    /// `None` where the caller names none.
    pub source: Option<&'a str>,
    pub gwas_pairs: Option<&'a gwas::GwasCatalogPairs>,
    pub revel_table: Option<&'a revel::RevelTable>,
    pub utrannotator_sweep: Option<&'a utrannotator::UtrannotatorSeedSweep>,
    pub downstream_sequences: Option<&'a downstream::DownstreamSequences>,
    pub splice_windows: Option<&'a maxentscan::SpliceWindows>,
    pub spliceai_table: Option<&'a spliceai::SpliceAiTable>,
    pub exon_spans: Option<&'a exon_spans::ExonSpans>,
}

/// The classes a scorer consults, from any of its threads.
pub trait Rules: Sync {
    /// Whether `klass` is a class of these rules.
    fn knows(&self, klass: &str) -> bool;
    /// Whether the class accepts the consult's two values of its field.
    fn accepts(&self, consult: &Consult) -> Result<bool, Refusal>;
    /// Whether `klass` is handed an empty value a row writes as written.
    fn transcript_reading(&self, klass: &str) -> bool;
    /// Whether `klass` is handed the engine's rows of the record.
    fn sibling_allele(&self, klass: &str) -> bool;
    /// Whether `klass` reads the reference bases a row's allele stands in for.
    fn record_reference(&self, klass: &str) -> bool;
}

/// The class under which a set-valued column agrees when its items agree as a multiset:
/// Ensembl VEP prints such a list in the order of a hash's keys, which carries no meaning.
pub const ORDER_ONLY: &str = "order_only";
/// The columns `order_only` reads as lists.
pub const ORDER_ONLY_COLUMNS: [&str; 4] = ["DOMAINS", "FREQS", "ZYG", "NearestExonJB"];
/// The column `order_only` reads as groups of `key=value` pairs joined by `:`, the pairs in
/// hash order within each group (UTRAnnotator.pm:282-288).
pub const ORDER_ONLY_PAIR_GROUP_COLUMN: &str = "5UTR_annotation";
/// The marker `order_only_fields` gives a JSON array field in place of a list separator.
pub const JSON_ARRAY: &str = "json";
/// The prefix `order_only_fields` gives the pair-group column's separator.
pub const PAIR_GROUPS: &str = "pairs:";
/// The JSON field of each order-only column with its separator.
pub const ORDER_ONLY_JSON_FIELDS: [(&str, &str); 3] = [
    ("domains", JSON_ARRAY),
    ("zyg", ","),
    ("nearestexonjb", ","),
];

pub const GWAS_HAPLOTYPE_CROSS_PRODUCT: &str = "gwas_haplotype_cross_product";
pub const GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE: &str =
    "gwas_interaction_first_snp_allele_whitespace";
pub const SPLICEVAULT_SIBLING_ALLELE_ROWS: &str = "splicevault_sibling_allele_rows";
pub const ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION: &str =
    "alphamissense_two_digit_transcript_version";
pub const SPLICEVAULT_JSON_LATER_ROWS: &str = "splicevault_json_later_rows";
pub const UTRANNOTATOR_ANNOTATION_PAIR_ORDER: &str = "utrannotator_annotation_pair_order";
pub const NEAREST_EXON_JB_ITEM_ORDER: &str = "nearest_exon_jb_item_order";
pub const UTRANNOTATOR_GAINED_START_WINDOW: &str = "utrannotator_gained_start_window";
pub const UTRANNOTATOR_CONSEQUENCE_PRECEDENCE: &str = "utrannotator_consequence_precedence";
pub const UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT: &str = "utrannotator_structural_allele_text";
pub const UTRANNOTATOR_LOST_START_AFTER_INSERTION: &str = "utrannotator_lost_start_after_insertion";
pub const UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE: &str =
    "utrannotator_gained_start_remade_in_place";
pub const UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY: &str =
    "utrannotator_frameshift_at_stop_boundary";
pub const UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON: &str =
    "utrannotator_lost_stop_written_before_stop_codon";
pub const LOEUF_ZERO_MINIMUM: &str = "loeuf_zero_minimum";
pub const MAXENTSCAN_OVERLAPPING_WINDOWS: &str = "maxentscan_overlapping_windows";
pub const ALPHAMISSENSE_SIBLING_ALLELE_ROWS: &str = "alphamissense_sibling_allele_rows";
pub const REVEL_SIBLING_ALLELE_ROWS: &str = "revel_sibling_allele_rows";
pub const REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS: &str = "revel_other_reference_amino_acid_rows";
pub const GNOMADMT_INFO_VALUE_SECOND_EQUALS: &str = "gnomadmt_info_value_second_equals";
pub const LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD: &str = "loftee_anc_allele_multiallelic_record";
pub const LOFTEE_FIRST_ALT_DONOR_FILTER: &str = "loftee_first_alt_donor_filter";
pub const LOFTEE_GERP_ZOOM_LEVEL_SUMMARY: &str = "loftee_gerp_zoom_level_summary";
pub const LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE: &str =
    "loftee_gerp_walk_without_conservation_database";
pub const SPLICE_FAMILY_SWAP: &str = "splice_family_swap";
pub const NEAREST_EXON_JB_STRETCH_PRIMING: &str = "nearest_exon_jb_stretch_priming";
pub const NEAREST_EXON_JB_SHIFTED_SPAN: &str = "nearest_exon_jb_shifted_span";
pub const NEAREST_EXON_JB_ID_LESS_EXON_FOLD: &str = "nearest_exon_jb_id_less_exon_fold";
pub const DOWNSTREAM_REVERSE_STRAND_ROTATION: &str = "downstream_reverse_strand_rotation";
pub const SPLICEAI_CHOPPED_MULTIALLELIC: &str = "spliceai_chopped_multiallelic";
pub const SPLICEAI_MULTI_ENTRY_CELL: &str = "spliceai_multi_entry_cell";
pub const PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING: &str = "plugin_contig_mitochondrial_spelling";
pub const DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M: &str = "dbnsfp_mitochondrial_record_named_m";
pub const HGVS_INTRON_OFFSET_ABUTTING_INSERTION: &str = "hgvs_intron_offset_abutting_insertion";
pub const HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE: &str = "hgvs_shifted_indel_written_span_gate";
pub const HGVS_SHIFTED_DELETION_CLAMPED_SPAN: &str = "hgvs_shifted_deletion_clamped_span";
pub const HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE: &str =
    "hgvs_clipped_reference_against_whole_alternate";
pub const HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END: &str = "hgvs_added_bases_left_at_the_span_end";
pub const HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS: &str =
    "hgvs_insertion_left_where_the_shift_stops";
pub const HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT: &str =
    "hgvs_added_bases_written_outside_the_transcript";
pub const BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES: &str =
    "blosum62_hgvs_shift_shared_across_alleles";
pub const DOWNSTREAM_UTR_FRAME_PAST_CODING_END: &str = "downstream_utr_frame_past_coding_end";
pub const DOWNSTREAM_START_PAST_OWN_CODON: &str = "downstream_start_past_own_codon";
pub const PLUGIN_HGVS_SHIFTED_CDS_COORDINATES: &str = "plugin_hgvs_shifted_cds_coordinates";
pub const PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN: &str =
    "plugin_value_follows_splice_term_past_the_record_span";
pub const PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR: &str =
    "plugin_value_follows_stop_codon_term_pair";
/// A field the engine alone writes on a row of a record the reference dropped before
/// annotation. A scorer settles it from the two outputs' rows; no consult accepts it.
pub const PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD: &str =
    "plugin_value_on_reference_skipped_record";

/// The items of a list-valued field: the `separator`-joined parts of a text format's
/// value, or the elements of a JSON array's inner text (`JSON_ARRAY`), each as its compact
/// text; under a `PAIR_GROUPS` separator the parts are groups of `:`-joined pairs, each
/// group read with its pairs sorted; `None` when the value is not such a list.
pub fn list_items(value: &str, separator: &str) -> Option<Vec<String>> {
    if let Some(separator) = separator.strip_prefix(PAIR_GROUPS) {
        if separator.is_empty() {
            return None;
        }
        return Some(
            value
                .split(separator)
                .map(|group| {
                    let mut pairs: Vec<&str> = group.split(':').collect();
                    pairs.sort_unstable();
                    pairs.join(":")
                })
                .collect(),
        );
    }
    if separator.is_empty() {
        return None;
    }
    if separator != JSON_ARRAY {
        return Some(value.split(separator).map(str::to_string).collect());
    }
    match pytext::json_loads(&format!("[{value}]")) {
        Ok(Json::List(items)) => Some(items.iter().map(pytext::json_compact).collect()),
        _ => None,
    }
}

/// Whether two list-valued field texts carry the same items the same number of times, in
/// any order; an empty or absent value agrees with itself alone.
pub fn same_items(expected: &str, actual: &str, separator: &str) -> bool {
    if blank(expected) || blank(actual) {
        return expected == actual;
    }
    let (Some(mut expected_items), Some(mut actual_items)) = (
        list_items(expected, separator),
        list_items(actual, separator),
    ) else {
        return false;
    };
    expected_items.sort_unstable();
    actual_items.sort_unstable();
    expected_items == actual_items
}

/// The fields the `order_only` class compares as multisets in the format `fmt`, each with
/// its list separator: `,` in the default and tab formats (`convert_arrayref`,
/// Utils.pm:185-198), `&` in VCF (VCF.pm:391), and in JSON the arrays of
/// `ORDER_ONLY_JSON_FIELDS`.
pub fn order_only_fields(fmt: &str) -> BTreeMap<String, String> {
    if fmt == "json" {
        return ORDER_ONLY_JSON_FIELDS
            .iter()
            .map(|(field, separator)| (field.to_string(), separator.to_string()))
            .collect();
    }
    let separator = if fmt == "vcf" { "&" } else { "," };
    let mut fields: BTreeMap<String, String> = ORDER_ONLY_COLUMNS
        .iter()
        .map(|column| (column.to_string(), separator.to_string()))
        .collect();
    fields.insert(
        ORDER_ONLY_PAIR_GROUP_COLUMN.to_string(),
        format!("{PAIR_GROUPS}{separator}"),
    );
    fields
}

/// The refusal of a class this build holds no rule for.
pub fn no_rule(klass: &str) -> Refusal {
    Refusal(format!("{klass}: no rule of this build"))
}

/// The tool's own classes: those of `data/field_classes.tsv`.
pub struct BuiltIn;

impl Rules for BuiltIn {
    fn knows(&self, klass: &str) -> bool {
        fields::class_fields(klass).is_some()
    }

    /// Whether the class accepts the differing field: the field is one the class reads and
    /// the two values have the class's shape. A value that is blank is handed to the rule
    /// as `-`. A class without a rule in this build is refused.
    fn accepts(&self, c: &Consult) -> Result<bool, Refusal> {
        let read = fields::class_fields(c.klass).ok_or_else(|| no_rule(c.klass))?;
        // dbNSFP writes a column under the column's own name (dbNSFP.pm:257-261), which no
        // row of the field table holds: such a field is the class's where the run gives it
        // to the plugin.
        let run_column = c.klass == DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M
            && c.source == Some(mitochondrial::DBNSFP_PLUGIN);
        if !field_named(c.name, read) && !run_column {
            return Ok(false);
        }
        let ev = if blank(c.expected) { "-" } else { c.expected };
        let av = if blank(c.actual) { "-" } else { c.actual };
        match c.klass {
            GWAS_HAPLOTYPE_CROSS_PRODUCT => Ok(gwas::haplotype_cross_product(c, ev, av)),
            GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE => {
                Ok(gwas::interaction_first_snp_allele_whitespace(c, ev, av))
            }
            SPLICEVAULT_SIBLING_ALLELE_ROWS => Ok(!blank(ev)
                && c.record_rows.is_some_and(|rows| {
                    sibling::sibling_allele_value(c.e, read, rows, false, c.candidate_rows)
                })),
            ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION | SPLICEVAULT_JSON_LATER_ROWS => {
                Ok(blank(ev) && !blank(av))
            }
            UTRANNOTATOR_ANNOTATION_PAIR_ORDER => Ok(utrannotator::reordered_pair_groups(ev, av)),
            NEAREST_EXON_JB_ITEM_ORDER => Ok(nearest_exon::reordered_items(ev, av)),
            UTRANNOTATOR_GAINED_START_WINDOW => Ok(utrannotator::gained_start_window(c)),
            UTRANNOTATOR_CONSEQUENCE_PRECEDENCE => Ok(utrannotator::consequence_precedence(c)),
            UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT => {
                Ok(utrannotator::structural_allele_text(c, ev, av))
            }
            UTRANNOTATOR_LOST_START_AFTER_INSERTION => {
                Ok(utrannotator::lost_start_after_insertion(c))
            }
            UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE => {
                Ok(utrannotator::gained_start_remade_in_place(c))
            }
            UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY => {
                Ok(utrannotator::frameshift_at_stop_boundary(c))
            }
            UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON => {
                Ok(utrannotator::lost_stop_written_before_stop_codon(c))
            }
            LOEUF_ZERO_MINIMUM => Ok(blank(ev) && av == "0"),
            MAXENTSCAN_OVERLAPPING_WINDOWS => {
                maxentscan::overlapping_windows(c.e, c.name, ev, av, c.splice_windows)
            }
            ALPHAMISSENSE_SIBLING_ALLELE_ROWS | REVEL_SIBLING_ALLELE_ROWS => Ok(!blank(ev)
                && c.record_rows.is_some_and(|rows| {
                    sibling::sibling_allele_value(c.e, read, rows, true, c.candidate_rows)
                })),
            REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS => revel::revel_other_reference_amino_acid_row(
                c.e,
                c.a.map(|a| &**a),
                ev,
                av,
                c.record_rows,
                c.revel_table,
            ),
            GNOMADMT_INFO_VALUE_SECOND_EQUALS => {
                let written = decode_extra(av);
                let (kept, sign, _) = pytext::partition(&written, "=");
                Ok(!blank(ev) && !sign.is_empty() && decode_extra(ev) == kept)
            }
            LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD => {
                loftee::anc_allele_multiallelic_record(c.e, c.a.map(|a| &**a), c.record_rows)
            }
            LOFTEE_FIRST_ALT_DONOR_FILTER => Ok(loftee::first_alt_donor_filter(
                c.e,
                c.a.map(|a| &**a),
                c.record_rows,
            )),
            LOFTEE_GERP_ZOOM_LEVEL_SUMMARY => Ok(loftee::gerp_zoom_level_summary(
                c.e,
                c.a.map(|a| &**a),
                c.name,
            )),
            LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE => Ok(
                loftee::gerp_walk_without_conservation_database(c.e, c.a.map(|a| &**a), ev, av),
            ),
            SPLICE_FAMILY_SWAP => term_pair::splice_family_swap(c, ev, av),
            NEAREST_EXON_JB_STRETCH_PRIMING => Ok(nearest_exon::stretch_priming(c.e, ev, av)),
            NEAREST_EXON_JB_SHIFTED_SPAN => nearest_exon::shifted_span(c.e, ev, av),
            NEAREST_EXON_JB_ID_LESS_EXON_FOLD => Ok(nearest_exon::id_less_exon_fold(c.e, ev, av)),
            DOWNSTREAM_REVERSE_STRAND_ROTATION => downstream::reverse_strand_rotation(c.e, ev, av),
            SPLICEAI_CHOPPED_MULTIALLELIC => spliceai::chopped_multiallelic(c, ev, av),
            SPLICEAI_MULTI_ENTRY_CELL => spliceai::multi_entry_cell(c, ev, av),
            PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING => Ok(
                mitochondrial::contig_mitochondrial_spelling(c.e, c.name, ev, av, c.record_rows),
            ),
            DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M => Ok(mitochondrial::dbnsfp_record_named_m(
                c.e,
                ev,
                av,
                c.record_rows,
                c.assembly,
            )),
            HGVS_INTRON_OFFSET_ABUTTING_INSERTION => {
                Ok(hgvs_offsets::intron_offset_abutting_insertion(
                    c.e,
                    c.a.map(|a| &**a),
                    c.name,
                    ev,
                    av,
                ))
            }
            HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE => hgvs_offsets::shifted_indel_written_span_gate(
                c.e,
                c.a.map(|a| &**a),
                c.name,
                ev,
                av,
            ),
            HGVS_SHIFTED_DELETION_CLAMPED_SPAN => {
                hgvs_offsets::shifted_deletion_clamped_span_field(
                    c.e,
                    c.a.map(|a| &**a),
                    c.name,
                    ev,
                    av,
                )
            }
            HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE => {
                hgvs_offsets::clipped_reference_field(c.e, c.a.map(|a| &**a), c.name, ev, av)
            }
            HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END => {
                hgvs_offsets::added_bases_left_at_the_span_end(c.e, c.a.map(|a| &**a), c.name)
            }
            HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS => {
                hgvs_offsets::insertion_left_where_the_shift_stops(
                    c.e,
                    c.a.map(|a| &**a),
                    c.name,
                    c.exon_spans,
                )
            }
            HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT => {
                hgvs_offsets::added_bases_written_outside_the_transcript(
                    c.e,
                    c.a.map(|a| &**a),
                    c.name,
                    ev,
                    av,
                )
            }
            BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES => {
                downstream::blosum62_shift_shared_across_alleles(c, ev, av)
            }
            DOWNSTREAM_UTR_FRAME_PAST_CODING_END
            | DOWNSTREAM_START_PAST_OWN_CODON
            | PLUGIN_HGVS_SHIFTED_CDS_COORDINATES => downstream::transcript_reading(c),
            PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN => {
                term_pair::follows_splice_term_past_the_record_span(c, ev, av)
            }
            PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR => {
                term_pair::follows_stop_codon_term_pair(c, ev, av)
            }
            PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD => Ok(false),
            _ => Err(no_rule(c.klass)),
        }
    }

    fn transcript_reading(&self, klass: &str) -> bool {
        fields::transcript_reading(klass)
    }

    fn sibling_allele(&self, klass: &str) -> bool {
        fields::sibling_allele(klass)
    }

    /// The classes of `#record_reference`, and the classes that read a transcript's
    /// sequence: on a row whose two Consequence sets differ those ask whether a class of a
    /// coding end takes the pair (`downstream::classes::downstream_class`), which reads the
    /// count of reference bases the row's allele stands in for.
    fn record_reference(&self, klass: &str) -> bool {
        fields::record_reference(klass) || fields::transcript_reading(klass)
    }
}
