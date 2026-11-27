// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The built-in field classes: what a rule reads and the dispatch by class name.

pub mod added_bases;
pub mod allele_number;
pub mod clipped;
pub mod coding_end;
pub mod coding_end_pairs;
pub mod cut_lists;
pub mod downstream;
pub mod edited_model;
pub mod file_transcript;
pub mod following;
pub mod frameshift_at_a_held_stop;
pub mod frameshift_count;
pub mod frameshift_no_changed_residue;
pub mod gwas;
pub mod hgvsp_stop;
pub mod inserted_span;
pub mod intron_ends;
pub mod kept_start;
pub mod length_edits;
pub mod loftee;
pub mod nearest_exon;
pub mod partial_codon_models;
pub mod record_alleles;
pub mod registry;
pub mod removed_start;
pub mod revel;
pub mod sequence_edges;
pub mod settled_stop;
pub mod shifted_indel;
pub mod shifted_start_test;
pub mod splice_windows;
pub mod spliceai;
pub mod start_test_past_the_allele;
pub mod structural;
pub mod term_pairs;
pub mod total_length;
pub mod transcript_edges;
pub mod transcript_model;
pub mod transcript_reference;
pub mod utrannotator;

use crate::classes as class;
use crate::coding_ends::terms::{CodingEndClass, LandedRule};
use crate::compare::CONSEQUENCE_FIELDS;
use crate::consts::{
    BAM_EDIT_FIELDS, CDNA_START_FIELDS, CDS_POSITION_FIELDS, DISTANCE_FIELDS, FLAGS_FIELDS,
    HGVSC_FIELDS, HGVS_INTRON_OFFSET_FIELDS, HGVS_INTRON_START_OFFSET_FIELDS, HGVS_OFFSET_FIELDS,
    INTEGER_RE_PATTERN, LOF_INFO_FIELDS, NMD_FIELDS, NON_ZERO_INTEGER_RE_PATTERN,
    SHIFT_LENGTH_FIELDS, STRAND_FIELDS, STRUCTURAL_INSERTION_WORD, UTRANNOTATOR_ANNOTATION_FIELDS,
};
use crate::entry::Entry;
use crate::error::Result;
use crate::fasta::ReferenceSequence;
use crate::fields::{as_number, entry_kind, field_named, OrganelleCodonTables};
use crate::lists::{reordered_items, reordered_pair_groups};
use crate::manifest::decode_extra;
use crate::nearest::{record_span, NearestSites};
use crate::pyfmt;
use crate::rules::added_bases::{
    added_bases_left_at_the_span_end, added_bases_written_outside_the_transcript,
};
use crate::rules::allele_number::{allele_number_place_pair, WRITTEN_ALLELE_LIST_FIELDS};
use crate::rules::clipped::{
    clipped_notation_intron_offsets, clipped_reference_against_whole_alternate,
    insertion_left_where_the_shift_stops,
};
use crate::rules::coding_end::{
    codons_read_past_the_coding_end, residues_read_past_the_coding_end,
};
use crate::rules::coding_end_pairs::{
    asked_at, coding_end_following_cell, frameshift_count_at_the_first_residue,
    stands_on_the_record, AtACodingEnd,
};
use crate::rules::downstream::{
    blosum62_later_allele_shared_shift, downstream_class_at, field_of, hgnc_id_copy_pair,
    later_allele_shared_shift, nmd_reads_shifted_pair, organelle_peptide_path, positive_shift,
    unrotated_shifted_insertion, DownstreamSequences, ONE_CODING_CHANGE_FIELDS,
};
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::edited_model::{
    allele_readings, cdna_pair_spans, edited_model_clipped_alleles,
    edited_model_seam_insertion_replacement, edited_model_shifted_insertion_bases,
    shifted_deletion_clamped_span, shifted_deletion_start_offsets,
};
use crate::rules::frameshift_count::frameshift_count_past_a_stop_read_as_a_residue;
use crate::rules::frameshift_no_changed_residue::frameshift_with_no_changed_residue;
use crate::rules::gwas::{gwas_annotated, gwas_pair_named, gwas_site, GwasCatalogPairs};
use crate::rules::hgvsp_stop::{
    extension_counted_from_a_stop_the_peptide_holds, stop_codon_protein_comparison,
    stop_named_on_incomplete_coding_sequence,
    stop_unchanged_where_inserted_residues_keep_the_protein,
    stop_unchanged_where_the_reference_writes_none,
};
use crate::rules::inserted_span::inserted_span_codons_cell;
use crate::rules::intron_ends::{
    hgvsp_intron_end_insertion, hgvsp_intron_start_insertion, IntronEnds,
};
use crate::rules::length_edits::failed_edit_offset_cell;
use crate::rules::loftee::{
    cdna_span_closed, gerp_totals, gerp_walk_the_candidate_alone_finished, lof_filter_item_side,
    nagnag_site_first_transcript_strand_row, notation_indel_length,
};
use crate::rules::nearest_exon::nearest_exon_items;
use crate::rules::record_alleles::{
    first_alternate, record_alleles, record_chromosome, record_indel_row, sibling_rows_of,
};
use crate::rules::registry::registry_class;
use crate::rules::revel::{revel_other_reference_amino_acid_row, sibling_allele_value, RevelTable};
use crate::rules::sequence_edges::SequenceEdges;
use crate::rules::shifted_indel::{
    intron_offset_abutting_insertion, shifted_indel_written_span_gate,
    shifted_indel_written_span_protein_gate, shifted_indel_written_span_start_test,
    start_codon_test_past_the_shifted_allele_rows,
};
use crate::rules::splice_windows::SpliceWindows;
use crate::rules::spliceai::{spliceai_multi_entry_cell_row, SpliceAiTable};
use crate::rules::structural::StructuralVariants;
use crate::rules::term_pairs::{
    acceptor_site_offsets, clamped_insertion_position, deleted_span_coding_reference_cell,
    freqs_differ_by_zero_read_as_na, intron_span_coding_reference_cell, shifted_insertion_end,
    splice_term_past_the_record_span_pair_at, symbolic_loss_named_a_gain_pair,
    term_pair_following_cell, STOP_CODON_PLUGIN_VALUE_PAIRS, TERM_PAIR_CLASSES,
};
use crate::rules::total_length::total_length_range_end_pair;
use crate::rules::transcript_edges::TranscriptEdges;
use crate::rules::transcript_reference::failed_edit_reference_cell;
use crate::rules::utrannotator::{
    annotation_groups, frameshift_at_stop_boundary_row, gained_start_remade_in_place_row,
    gained_start_window_row, lost_start_after_insertion_row,
    lost_stop_written_before_stop_codon_row, structural_allele_text_row,
    structural_copy_number_loss_row, structural_deletion_span_row, structural_insertion_length_row,
    structural_inversion_span_row, utrannotator_precedence_row, without_reference_row,
    UtrannotatorSeedSweep, STRUCTURAL_COPY_NUMBER_WORD, STRUCTURAL_INVERSION_WORD,
    UTRANNOTATOR_CONSEQUENCE_FIELDS,
};
use crate::shared::{blank, consequence_set, feature_kind, insertion_row, ABSENT, FRAMESHIFT_TERM};
use std::collections::BTreeSet;
use std::sync::LazyLock;

pub const TER_COUNT_RE_PATTERN: &str = r"^(.*(?:fs|ext)Ter)(\d+|\?)(?:\n?\z)";
pub static TER_COUNT_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(TER_COUNT_RE_PATTERN).expect("the pattern compiles"));

pub const GWAS_STUDY_FIELDS: &[&str] = &["GWAS_study"];

pub const MAXENTSCAN_WINDOW_TEXT_FIELDS: &[&str] = &["MaxEntScan_ref_seq", "MaxEntScan_alt_seq"];

pub const MAXENTSCAN_DONOR_TEXT_RE_PATTERN: &str = r"[ACGT]{9}";
pub static MAXENTSCAN_DONOR_TEXT_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(MAXENTSCAN_DONOR_TEXT_RE_PATTERN).expect("the pattern compiles")
});

pub const MAXENTSCAN_ACCEPTOR_TEXT_RE_PATTERN: &str = r"[ACGT]{23}";
pub static MAXENTSCAN_ACCEPTOR_TEXT_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(MAXENTSCAN_ACCEPTOR_TEXT_RE_PATTERN).expect("the pattern compiles")
});

pub const GNOMAD_PLUGIN_PREFIX: &str = "gnomad_";

pub const ANC_ALLELE_FILTER: &str = "ANC_ALLELE";

pub const LOF_SUBSTITUTION_TERMS: &[&str] = &[
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
];

pub const GC_TO_GT_DONOR_FILTER: &str = "GC_TO_GT_DONOR";

pub const LOF_DONOR_TERM: &str = "splice_donor_variant";

pub const LOF_TRUNCATING_TERMS: &[&str] = &["frameshift_variant", "stop_gained"];

pub const END_TRUNC_FILTER: &str = "END_TRUNC";

pub const GERP_END_TRUNC_CUTOFF: f64 = -58.0;

pub const GERP_EXACT_BELOW_BASES: i64 = 80;

pub const LAST_EXON_REACH: i64 = 50;

pub const FRAMESHIFT_INTRON_STRETCH: i64 = 12;

pub const MITOCHONDRIAL_SPELLINGS: &[&str] = &["MT", "chrM", "chrMT"];

/// The record spellings dbNSFP.pm:313 names `M` to the plugin's table.
pub const DBNSFP_NAMED_M_SPELLINGS: &[&str] = &["MT", "chrMT"];

pub const GRCH37_NAMES: &[&str] = &["grch37", "hg19", "37"];

/// The `gnomAD_` fields are two plugins': gnomADc's end in `_cov`, the others are gnomADMt's.
pub const GNOMADC_FIELD_SUFFIX: &str = "_cov";

pub const CDNA_COORDINATE_RE_PATTERN: &str = r"\*?-?\d+(?:[+-]\d+)?";
pub static CDNA_COORDINATE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(CDNA_COORDINATE_RE_PATTERN).expect("the pattern compiles"));

pub const INTRON_TERM: &str = "intron_variant";

pub const SPLICE_ACCEPTOR_TERM: &str = "splice_acceptor_variant";

pub const SPLICE_SITE_TERMS: &[&str] = &["splice_acceptor_variant", "splice_donor_variant"];

/// The side tables a comparison reads beside the entries; `None` where the run names none.
#[derive(Debug, Default)]
pub struct Side {
    pub gwas_pairs: Option<GwasCatalogPairs>,
    pub utrannotator_sweep: Option<UtrannotatorSeedSweep>,
    pub codon_tables: Option<OrganelleCodonTables>,
    pub partial_codon_models: Option<partial_codon_models::PartialCodonModels>,
    pub length_edits: Option<length_edits::LengthEdits>,
    /// The reference sequence the candidate was run on, which
    /// `stop_codon_completed_by_polyadenylation` reads past a transcript's last base.
    pub reference_sequence: Option<ReferenceSequence>,
    pub downstream_sequences: Option<DownstreamSequences>,
    pub transcript_edges: Option<TranscriptEdges>,
    pub intron_ends: Option<IntronEnds>,
    pub revel_table: Option<RevelTable>,
    pub spliceai_table: Option<SpliceAiTable>,
    pub splice_windows: Option<SpliceWindows>,
    pub sequence_edges: Option<SequenceEdges>,
    pub nearest_sites: Option<NearestSites>,
    pub structural_variants: Option<StructuralVariants>,
}

impl Side {
    /// Asks `rule` with what the side tables hold of the coding ends of the transcript of
    /// the reference entry `e` (`coding_end_pairs::asked_at`).
    pub fn asked_at<T>(
        &self,
        e: &Entry,
        rule: impl FnOnce(&AtACodingEnd<'_>) -> Result<T>,
    ) -> Result<T> {
        asked_at(
            self.partial_codon_models.as_ref(),
            self.reference_sequence.as_ref(),
            e,
            rule,
        )
    }
}

/// What `field_class_accepts` reads beside the two values.
#[derive(Debug)]
pub struct RuleContext<'a> {
    pub record_rows: Option<&'a [&'a Entry]>,
    pub candidate_rows: Option<&'a [&'a Entry]>,
    pub declared: &'a [&'a str],
    pub assembly: Option<&'a str>,
    pub dbnsfp_columns: &'a [String],
    pub side: &'a Side,
}

/// Whether the built-in field class `klass` accepts the differing field `name` of the
/// entry `e`: the field is one the class reads (`FIELD_CLASSES`) and the two values have
/// the class's shape; `a` is the candidate's entry, read by a class whose shape spans two
/// fields of the row, and `record_rows` the reference's entries of `e`'s own input record
/// (every allele's rows, at least those on its feature), which the caller vouches are one
/// record's: a class that reads a record's alleles accepts nothing without them.
/// `loftee_nagnag_site_first_transcript_strand` reads the record's rows on other features
/// among them and accepts nothing without them.
/// `gwas_pairs` is the catalog the two GWAS classes read and accept nothing without.
/// `revel_table` is the table `revel_other_reference_amino_acid_rows` reads and accepts
/// nothing without.
/// `spliceai_table` is the score files `spliceai_multi_entry_cell` reads and accepts
/// nothing without; it counts the rows the class is asked about and refuses under them.
/// `assembly` is the run's assembly as its manifest names it, which
/// `dbnsfp_mitochondrial_record_named_m` reads: without one the class accepts nothing.
/// `dbnsfp_columns` is the columns the manifest's run hands the dbNSFP plugin
/// (`Documented::dbnsfp_columns`), which the same class reads as its fields beside those of
/// `FIELD_CLASSES`.
/// `codon_tables` is the candidate cache's codon tables, without which
/// `organelle_peptide_paths_translation_table` accepts nothing.
/// `partial_codon_models` is the coding ends of the candidate cache's transcripts, without
/// which `hgvsp_stop_named_on_incomplete_coding_sequence` reads the row alone, the five
/// classes of a coding end take nothing, `hgvsp_frameshift_count_at_the_first_residue`
/// accepts nothing, and a class of a stop codon or a start codon that
/// stands reads the two rows alone; beside it such a class takes a pair only where the
/// row's transcript shows its codon (`coding_end_pairs::stands_beside_a_cache`).
/// `length_edits` is the candidate cache's edits that change the length of a transcript's
/// cDNA, without which `refseq_failed_edit_model_hgvsc_offset` accepts nothing.
/// `downstream_sequences` is the candidate cache's transcripts, without which the classes
/// of `TRANSCRIPT_READING_CLASSES` and `hgvsp_frameshift_count_past_a_stop_read_as_a_residue`
/// accept nothing and `hgvsp_stop_codon_protein_comparison` leaves more than one residue
/// before an inserted stop a difference and reads an extension's count from the row alone.
/// `stop_codon_settled_before_unknown_bases` reads from them the protein notation and the
/// HGVS_OFFSET of an insertion written in an intron and carried into the coding sequence
/// (`settled_stop::carried_insertion_cell`), and without them its rows with a `Codons` cell
/// alone.
/// `transcript_edges` is the candidate cache's exon spans, without which
/// `hgvs_clipped_reference_against_whole_alternate` is its open-side test alone,
/// `hgvs_insertion_left_where_the_shift_stops` reads a `dup` by its count of bases alone and
/// `hgvsp_shifted_indel_written_span_gate` holds the exon's first coding base to the bound
/// its notation gives.
/// `candidate_rows` is the candidate's entries of the same input record, without which the
/// three sibling-allele classes accept nothing.
/// `utrannotator_sweep` is the seeded runs `utrannotator_consequence_precedence` reads and
/// accepts nothing without, `utrannotator_gained_start_remade_in_place` reads on a row
/// whose reference writes another term than the gained start, and
/// `utrannotator_lost_stop_written_before_stop_codon` reads on a row whose candidate alone
/// reports the lost stop. `intron_ends` is the candidate
/// cache's intron ends, which
/// `hgvsp_intron_end_insertion` and `hgvsp_intron_start_insertion` read and accept nothing
/// without.
/// `splice_windows` is the candidate cache's introns, which
/// `maxentscan_overlapping_windows` reads and accepts nothing without.
/// `declared` is the classes named with `klass`: a Downstream value two causes change is
/// accepted only where it names the second cause's class too (`downstream_class`).
/// A class of `TERM_PAIR_FOLLOWING_CLASSES` accepts a cell of a row its own predicate
/// takes, where `term_pair_following_cell` reads the two values from the two sets of terms,
/// or, for a class of a coding end, where the cell is the one the pair gives
/// (`coding_end_pairs::coding_end_following_cell`); it reads `a` and accepts nothing without it.
///
/// `a` is `None` where the caller holds no candidate entry; the record's rows, the declared
/// columns, the assembly, the dbNSFP columns and the side tables are read from `context`.
pub fn field_class_accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    a: Option<&Entry>,
    context: &RuleContext,
) -> Result<bool> {
    let tables = class::tables();
    let Some(patterns) = tables.field_classes_of(klass) else {
        crate::uncaught!("{}", pyfmt::repr_str(klass));
    };
    // dbNSFP writes a column under the column's own name (dbNSFP.pm:257-261), which the
    // run's words state and no table of the comparison holds. The JSON format writes that
    // name in lower case; a name that differs from it in case otherwise is another field.
    let run_column = klass == class::DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M
        && context
            .dbnsfp_columns
            .iter()
            .any(|column| column == name || pyfmt::lower(column) == name);
    if !field_named(name, patterns) && !run_column {
        return Ok(false);
    }
    let ev = if blank(expected_value) {
        "-"
    } else {
        expected_value
    };
    let av = if blank(actual_value) {
        "-"
    } else {
        actual_value
    };
    let record_rows = context.record_rows;
    let side = context.side;
    if tables.in_term_pair_following_classes(klass) {
        let Some(a) = a else {
            return Ok(false);
        };
        let Some((_, pair)) = TERM_PAIR_CLASSES.iter().find(|(held, _)| *held == klass) else {
            crate::uncaught!("{}", pyfmt::repr_str(klass));
        };
        return side.asked_at(e, |at| {
            // The following classes that read a cell otherwise than `term_pair_following_cell`
            // does, beside the classes of a coding end: the residues and the protein notation
            // of a stop codon settled before unknown bases, and the protein notation of a start
            // codon removed without a 5' UTR and of a start codon kept without a 5' UTR. A row
            // the first class's predicate does not take has no `Codons` cell the class reads,
            // and is read beside the candidate cache's transcripts.
            if klass == class::STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES {
                if pair(e, a, at)? {
                    return Ok(settled_stop::stop_codon_settled_before_unknown_bases_cell(
                        e, a, name, ev, av,
                    ));
                }
                let Some(sequences) = side.downstream_sequences.as_ref() else {
                    return Ok(false);
                };
                return settled_stop::carried_insertion_cell(e, a, name, ev, av, sequences);
            }
            if !pair(e, a, at)? {
                return Ok(false);
            }
            if klass == class::START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR {
                return Ok(
                    removed_start::start_codon_removed_without_five_prime_utr_cell(
                        e, a, name, ev, av,
                    ),
                );
            }
            if klass == class::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR {
                return Ok(kept_start::start_codon_kept_without_five_prime_utr_cell(
                    e, a, name, ev, av,
                ));
            }
            // This class reads its cells otherwise as well, from the two values as the
            // entries carry them: an empty plugin value is a field the candidate writes,
            // which a `-` is not.
            if klass == class::FRAMESHIFT_INSIDE_A_HELD_STOP_CODON {
                return Ok(
                    frameshift_at_a_held_stop::frameshift_inside_a_held_stop_codon_cell(
                        e,
                        a,
                        name,
                        expected_value,
                        actual_value,
                    ),
                );
            }
            Ok(match CodingEndClass::of_name(klass) {
                Some(class) => coding_end_following_cell(class, e, a, name, ev, av, at),
                None => term_pair_following_cell(e, a, name, ev, av),
            })
        });
    }
    if klass == class::GXF_REVERSE_STRAND_UTR_SEQUENCE {
        return Ok(file_transcript::file_transcript_stop_count_pair(e, ev, av));
    }
    if klass == class::ALTERNATE_CODON_PAST_THE_CODING_END {
        if field_named(name, CODONS_FIELDS) {
            // A row names its first residue in `Protein_position`, up to a `-` or the `/`
            // of a total length; a JSON record in `protein_start`.
            let first_residue = match e.fields.get("Protein_position") {
                Some(position) => pyfmt::partition(pyfmt::partition(position, "/").0, "-").0,
                None => e.fields.get("protein_start").map_or("", String::as_str),
            };
            return Ok(codons_read_past_the_coding_end(first_residue, ev, av));
        }
        return Ok(a.is_some_and(|a| {
            residues_read_past_the_coding_end(
                &field_of(Some(e), CODONS_FIELDS),
                &field_of(Some(a), CODONS_FIELDS),
                ev,
                av,
            )
        }));
    }
    if klass == class::ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE {
        let (Some(codon_tables), Some(a)) = (side.codon_tables.as_ref(), a) else {
            return Ok(false);
        };
        return Ok(organelle_peptide_path(e, a, ev, av) && codon_tables.table(e)?.is_some());
    }
    if klass == class::CONSERVATION_MAX_TEXT_ORDER {
        return Ok(match (as_number(ev), as_number(av)) {
            (Some(reference), Some(candidate)) => reference < candidate && candidate < 0.0,
            _ => false,
        });
    }
    if klass == class::GWAS_HAPLOTYPE_CROSS_PRODUCT
        || klass == class::GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE
    {
        let Some(pairs) = side.gwas_pairs.as_ref() else {
            return Ok(false);
        };
        let Some((chromosome, position)) = gwas_site(e, record_rows) else {
            return Ok(false);
        };
        let site = (chromosome.as_str(), position.as_str());
        if klass == class::GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE {
            if !blank(ev) || gwas_annotated(e) {
                return Ok(false);
            }
            // The plugin writes a field the catalog row leaves empty as a defined empty
            // string (`GWAS_odds_ratio=`), which the reference's unannotated row has no key
            // for: the class's beside the candidate's filled fields, as the cross-product
            // class reads it.
            if blank(av) && !(actual_value.is_empty() && a.is_some_and(gwas_annotated)) {
                return Ok(false);
            }
            let studies = field_of(a, GWAS_STUDY_FIELDS);
            return Ok(gwas_pair_named(e, &studies, site, &pairs.candidate_only));
        }
        let studies = field_of(Some(e), GWAS_STUDY_FIELDS);
        if !gwas_pair_named(e, &studies, site, &pairs.reference_only) {
            return Ok(false);
        }
        // The plugin writes a field the catalog row leaves empty as a defined empty string
        // (`GWAS_beta_coef=`), which a row the candidate never annotated has no key for;
        // such a pair is the class's when the reference row carries any `GWAS_` value. The
        // empty string is the reference's: a key the candidate alone writes stays a
        // difference.
        if !blank(ev) || blank(av) {
            return Ok(!blank(ev)
                || (expected_value.is_empty() && blank(actual_value) && gwas_annotated(e)));
        }
        // A value the candidate alone writes is its own catalog row's where it reports
        // another study than the reference, whose row leaves the field empty.
        let own_study = field_of(a, GWAS_STUDY_FIELDS);
        return Ok(gwas_annotated(e) && !blank(&own_study) && own_study != studies);
    }
    if klass == class::SPLICEVAULT_SIBLING_ALLELE_ROWS {
        let Some(rows) = record_rows.filter(|_| !blank(ev)) else {
            return Ok(false);
        };
        return Ok(sibling_allele_value(
            e,
            patterns,
            rows,
            false,
            context.candidate_rows,
        ));
    }
    if klass == class::TABIX_ALLELE_CASE_INSENSITIVE_MATCH
        && field_named(name, &[GNOMAD_PLUGIN_PREFIX])
        && pyfmt::lower(name).ends_with(GNOMADC_FIELD_SUFFIX)
    {
        // gnomADc.pm:163 reads its table's rows by position and matches no allele.
        return Ok(false);
    }
    if klass == class::TABIX_ALLELE_CASE_INSENSITIVE_MATCH
        && field_named(name, &[gwas::GWAS_PREFIX])
        && blank(ev)
        && actual_value.is_empty()
    {
        // The GWAS plugin writes a field the catalog row leaves empty as a defined empty
        // string (`GWAS_odds_ratio=`), which the reference's unmatched row has no key for:
        // the class's beside the candidate's filled fields.
        return Ok(a.is_some_and(gwas_annotated) && !gwas_annotated(e));
    }
    if klass == class::ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION
        || klass == class::SPLICEVAULT_JSON_LATER_ROWS
        || klass == class::TABIX_ALLELE_CASE_INSENSITIVE_MATCH
    {
        return Ok(blank(ev) && !blank(av));
    }
    if klass == class::CADD_STRUCTURAL_ALLELE_MATCH {
        return Ok(blank(ev)
            && !blank(av)
            && e.allele != "-"
            && !pyfmt::upper(&e.allele)
                .chars()
                .all(|base| "ACGTN".contains(base)));
    }
    if klass == class::UTRANNOTATOR_ANNOTATION_PAIR_ORDER {
        return Ok(reordered_pair_groups(ev, av));
    }
    if klass == class::NEAREST_EXON_JB_ITEM_ORDER {
        return Ok(reordered_items(ev, av));
    }
    if klass == class::UTRANNOTATOR_GAINED_START_WINDOW {
        return Ok(a.is_some_and(|a| gained_start_window_row(e, a, record_rows)));
    }
    if klass == class::UTRANNOTATOR_CONSEQUENCE_PRECEDENCE {
        let (Some(a), Some(sweep)) = (a, side.utrannotator_sweep.as_ref()) else {
            return Ok(false);
        };
        if !utrannotator_precedence_row(e, a) {
            return Ok(false);
        }
        // An annotation the two entries share is no differing field of its own, so the
        // consequence is accepted only where the runs write it beside both terms.
        let expected = field_of(Some(e), UTRANNOTATOR_ANNOTATION_FIELDS);
        let actual = field_of(Some(a), UTRANNOTATOR_ANNOTATION_FIELDS);
        let shared = expected != ABSENT
            && actual != ABSENT
            && annotation_groups(&expected) == annotation_groups(&actual);
        return Ok(sweep.writes(
            e,
            a,
            shared || field_named(name, UTRANNOTATOR_ANNOTATION_FIELDS),
        ));
    }
    if klass == class::UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT {
        let Some(a) = a else {
            return Ok(false);
        };
        if blank(av) {
            return Ok(!blank(ev) && structural_allele_text_row(e, a));
        }
        if e.allele == STRUCTURAL_INSERTION_WORD {
            return Ok(structural_insertion_length_row(e, a));
        }
        if e.allele == STRUCTURAL_INVERSION_WORD {
            return Ok(structural_inversion_span_row(e, a, record_rows));
        }
        if e.allele == STRUCTURAL_COPY_NUMBER_WORD {
            return Ok(structural_copy_number_loss_row(e, a, record_rows));
        }
        return Ok(structural_deletion_span_row(e, a));
    }
    if klass == class::UTRANNOTATOR_LOST_START_AFTER_INSERTION {
        return Ok(a.is_some_and(|a| lost_start_after_insertion_row(e, a, record_rows)));
    }
    if klass == class::UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE {
        let sweep = side.utrannotator_sweep.as_ref();
        return Ok(a.is_some_and(|a| gained_start_remade_in_place_row(e, a, record_rows, sweep)));
    }
    if klass == class::UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY {
        // Two entries that both report the frameshift write one term: the consequence is the
        // class's only where one of them fills neither field.
        if field_named(name, UTRANNOTATOR_CONSEQUENCE_FIELDS) && !(blank(ev) || blank(av)) {
            return Ok(false);
        }
        return Ok(a.is_some_and(|a| frameshift_at_stop_boundary_row(e, a)));
    }
    if klass == class::UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON {
        let sweep = side.utrannotator_sweep.as_ref();
        return Ok(
            a.is_some_and(|a| lost_stop_written_before_stop_codon_row(e, a, record_rows, sweep))
        );
    }
    if klass == class::UTRANNOTATOR_WITHOUT_REFERENCE {
        return Ok(a.is_some_and(|a| without_reference_row(e, a)));
    }
    if klass == class::LOEUF_ZERO_MINIMUM {
        return Ok(blank(ev) && av == "0");
    }
    if klass == class::MAXENTSCAN_OVERLAPPING_WINDOWS {
        let Some(windows) = side.splice_windows.as_ref() else {
            return Ok(false);
        };
        if !windows.holds(e)? {
            return Ok(false);
        }
        if field_named(name, MAXENTSCAN_WINDOW_TEXT_FIELDS) {
            return Ok(MAXENTSCAN_DONOR_TEXT_WHOLE_RE.is_match(ev)
                && MAXENTSCAN_ACCEPTOR_TEXT_WHOLE_RE.is_match(av));
        }
        return Ok(as_number(ev).is_some() && as_number(av).is_some());
    }
    if klass == class::ALPHAMISSENSE_SIBLING_ALLELE_ROWS
        || klass == class::REVEL_SIBLING_ALLELE_ROWS
    {
        let Some(rows) = record_rows.filter(|_| !blank(ev)) else {
            return Ok(false);
        };
        return Ok(sibling_allele_value(
            e,
            patterns,
            rows,
            true,
            context.candidate_rows,
        ));
    }
    if klass == class::REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS {
        return revel_other_reference_amino_acid_row(
            e,
            a,
            ev,
            av,
            record_rows,
            side.revel_table.as_ref(),
        );
    }
    if klass == class::GNOMADMT_INFO_VALUE_SECOND_EQUALS {
        let decoded = decode_extra(av);
        let (kept, sign, _) = pyfmt::partition(&decoded, "=");
        return Ok(!blank(ev) && !sign.is_empty() && decode_extra(ev) == kept);
    }
    if klass == class::LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD {
        let Some(rows) = record_rows else {
            return Ok(false);
        };
        if lof_filter_item_side(e, a, ANC_ALLELE_FILTER).as_deref() != Some("candidate") {
            return Ok(false);
        }
        if !carries_one_of(&consequence_set(&e.consequence), LOF_SUBSTITUTION_TERMS) {
            return Ok(false);
        }
        if !matches!(e.allele.as_str(), "A" | "C" | "G" | "T") || record_indel_row(e, rows)? {
            return Ok(false);
        }
        return Ok(!sibling_rows_of(e, rows).is_empty());
    }
    if klass == class::LOFTEE_FIRST_ALT_DONOR_FILTER {
        let Some(rows) = record_rows else {
            return Ok(false);
        };
        if lof_filter_item_side(e, a, GC_TO_GT_DONOR_FILTER).is_none() {
            return Ok(false);
        }
        if !carries_one_of(&consequence_set(&e.consequence), &[LOF_DONOR_TERM])
            || sibling_rows_of(e, rows).is_empty()
        {
            return Ok(false);
        }
        // The first alternate's own row reads its own base in both engines, and rows that
        // state no allele order do not tell it from a later allele's.
        return Ok(first_alternate(e, rows).is_some_and(|first| first != e.allele));
    }
    if klass == class::LOFTEE_GERP_ZOOM_LEVEL_SUMMARY {
        let Some((reference_total, candidate_total, items)) = gerp_totals(e, a) else {
            return Ok(false);
        };
        if !carries_one_of(&consequence_set(&e.consequence), LOF_TRUNCATING_TERMS) {
            return Ok(false);
        }
        let item = |key: &str| as_number(items.get(key).map_or("", String::as_str));
        // `BP_DIST` is the interval's bases less one: 79 walks the 80-base interval.
        let Some(walked) = item("BP_DIST") else {
            return Ok(false);
        };
        if walked + 1.0 < GERP_EXACT_BELOW_BASES as f64 {
            return Ok(false);
        }
        if field_named(name, LOF_INFO_FIELDS) {
            return Ok(true);
        }
        let (Some(carrying), Some(reach)) = (
            lof_filter_item_side(e, a, END_TRUNC_FILTER),
            item("DIST_FROM_LAST_EXON"),
        ) else {
            return Ok(false);
        };
        if reach > LAST_EXON_REACH as f64 {
            return Ok(false);
        }
        let (lower, higher) = if carrying == "reference" {
            (reference_total, candidate_total)
        } else {
            (candidate_total, reference_total)
        };
        return Ok(lower <= GERP_END_TRUNC_CUTOFF && GERP_END_TRUNC_CUTOFF < higher);
    }
    if klass == class::LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE {
        let Some(a) = a.filter(|_| blank(ev) && !blank(av)) else {
            return Ok(false);
        };
        // LoF.pm:180 and :251 start the walk on a row that holds either term, and the
        // candidate's plugin reads its own row's terms.
        let truncating =
            |row: &Entry| carries_one_of(&consequence_set(&row.consequence), LOF_TRUNCATING_TERMS);
        return Ok(truncating(e) && truncating(a) && gerp_walk_the_candidate_alone_finished(e, a));
    }
    if klass == class::LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND {
        let (Some(a), Some(rows)) = (a, record_rows) else {
            return Ok(false);
        };
        return Ok(nagnag_site_first_transcript_strand_row(e, a, rows));
    }
    if klass == class::SPLICE_FAMILY_SWAP {
        let Some(a) = a.filter(|_| !blank(ev) && blank(av)) else {
            return Ok(false);
        };
        // Downstream.pm:95-99 write on a frameshift row none of whose terms names a splice
        // site.
        let terms = consequence_set(&e.consequence);
        if !carries_one_of(&terms, &[FRAMESHIFT_TERM])
            || terms.split(',').any(|term| term.contains("splice"))
        {
            return Ok(false);
        }
        let pair = registry_class(
            &terms,
            &consequence_set(&a.consequence),
            &e.allele,
            &[class::SPLICE_FAMILY_SWAP],
        );
        return Ok(pair.as_deref() == Some(class::SPLICE_FAMILY_SWAP));
    }
    if klass == class::NEAREST_EXON_JB_STRETCH_PRIMING {
        let distance = field_of(Some(e), DISTANCE_FIELDS);
        return Ok(blank(ev)
            && !blank(av)
            && pyfmt::is_digit(&distance)
            && whole_number(&distance) <= i128::from(FRAMESHIFT_INTRON_STRETCH));
    }
    if klass == class::NEAREST_EXON_JB_SHIFTED_SPAN {
        if !positive_shift(e)? || feature_kind(&entry_kind(e)) != "transcript" {
            return Ok(false);
        }
        if blank(ev) || blank(av) {
            return Ok(blank(ev) != blank(av));
        }
        let (Some(expected_items), Some(actual_items)) =
            (nearest_exon_items(ev), nearest_exon_items(av))
        else {
            return Ok(false);
        };
        let exons = |items: &[(String, String, String, String)]| -> BTreeSet<String> {
            items.iter().map(|(exon_id, ..)| exon_id.clone()).collect()
        };
        return Ok(exons(&expected_items) != exons(&actual_items));
    }
    if klass == class::NEAREST_START_SITE_ONE_BASE_ABOVE {
        // The class accepts nothing without the candidate cache's start sites. A row is
        // placed by its Location, a JSON record's entry by the span its own keys give.
        let Some(sites) = side.nearest_sites.as_ref() else {
            return Ok(false);
        };
        let field = |name: &str| e.fields.get(name).map(String::as_str);
        return Ok(
            match e.location.as_deref().or(e.allele_location.as_deref()) {
                Some(location) => sites.site_above_start_pair(location, ev, av),
                None => record_span(field("seq_region_name"), field("start"), field("end"))
                    .is_some_and(|(sequence, low, high)| {
                        sites.site_above_start_span(&sequence, low, high, ev, av)
                    }),
            },
        );
    }
    if klass == class::NEAREST_EXON_JB_ID_LESS_EXON_FOLD {
        let (Some(expected_items), Some(actual_items)) =
            (nearest_exon_items(ev), nearest_exon_items(av))
        else {
            return Ok(false);
        };
        let [(_, folded_distance, folded_side, _)] = &expected_items[..] else {
            return Ok(false);
        };
        if expected_items
            .iter()
            .chain(&actual_items)
            .any(|(exon_id, ..)| !exon_id.is_empty())
        {
            return Ok(false);
        }
        if !actual_items
            .iter()
            .any(|(_, distance, side, _)| distance == folded_distance && side == folded_side)
        {
            return Ok(false);
        }
        // One candidate item beside another length is the fold of two overlapped exons: the
        // row begins in one exon and ends in another, across both splice sites between
        // them.
        let terms = consequence_set(&e.consequence);
        return Ok(actual_items.len() > 1
            || (SPLICE_SITE_TERMS
                .iter()
                .all(|term| carries_one_of(&terms, &[term]))
                && cdna_span_closed(e)));
    }
    if klass == class::DOWNSTREAM_REVERSE_STRAND_ROTATION {
        if blank(ev)
            || blank(av)
            || !carries_one_of(&consequence_set(&e.consequence), &[FRAMESHIFT_TERM])
        {
            return Ok(false);
        }
        if !any_field_reads(e, STRAND_FIELDS, "-1") {
            return Ok(false);
        }
        let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
        let length = notation_indel_length(&field_of(Some(e), HGVSC_FIELDS));
        let Some(length) =
            length.filter(|length| INTEGER_WHOLE_RE.is_match(&offset) && *length > 0)
        else {
            return Ok(false);
        };
        return Ok(shifted_past_and_off_a_multiple(&offset, length));
    }
    if klass == class::SHIFTED_ALLELE_ROTATION_BY_STRAND {
        return Ok(a.is_some_and(|a| {
            unrotated_shifted_insertion(e, a, ev, av, &ONE_CODING_CHANGE_FIELDS)
        }));
    }
    if klass == class::SPLICEAI_CHOPPED_MULTIALLELIC {
        let Some(rows) = record_rows.filter(|_| blank(ev) && !blank(av)) else {
            return Ok(false);
        };
        let sequence =
            !e.allele.is_empty() && e.allele.bytes().all(|base| b"ACGTN".contains(&base));
        if !sequence || record_indel_row(e, rows)? {
            return Ok(false);
        }
        return Ok(record_alleles(e, rows).contains("-"));
    }
    if klass == class::SPLICEAI_MULTI_ENTRY_CELL {
        return spliceai_multi_entry_cell_row(e, a, name, ev, av, side.spliceai_table.as_ref());
    }
    if klass == class::PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING {
        if !blank(ev) || blank(av) {
            return Ok(false);
        }
        if field_named(name, &[GNOMAD_PLUGIN_PREFIX])
            && !pyfmt::lower(name).ends_with(GNOMADC_FIELD_SUFFIX)
        {
            return Ok(false);
        }
        return Ok(record_chromosome(e, record_rows)
            .is_some_and(|chromosome| MITOCHONDRIAL_SPELLINGS.contains(&chromosome.as_str())));
    }
    if klass == class::DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M {
        if blank(ev) || !blank(av) {
            return Ok(false);
        }
        let on_grch37 = context.assembly.is_some_and(|assembly| {
            GRCH37_NAMES.contains(&pyfmt::lower(pyfmt::strip(assembly)).as_str())
        });
        return Ok(on_grch37
            && record_chromosome(e, record_rows).is_some_and(|chromosome| {
                DBNSFP_NAMED_M_SPELLINGS.contains(&chromosome.as_str())
            }));
    }
    if klass == class::REFSEQ_EDITED_MODEL_HGVSC_OFFSET {
        if blank(ev) || blank(av) || !any_field_reads(e, BAM_EDIT_FIELDS, "OK") {
            return Ok(false);
        }
        let (expected_name, _, expected_body) = pyfmt::partition(ev, ":");
        let (actual_name, _, actual_body) = pyfmt::partition(av, ":");
        if expected_name != actual_name
            || CDNA_COORDINATE_RE.replace_all(expected_body, "#")
                != CDNA_COORDINATE_RE.replace_all(actual_body, "#")
        {
            return Ok(false);
        }
        let same_lengths = cdna_pair_spans(expected_body)
            .into_iter()
            .zip(cdna_pair_spans(actual_body))
            .all(|spans| match spans {
                (Some(expected), Some(actual)) => expected == actual,
                _ => true,
            });
        return Ok(same_lengths
            || edited_model::deletion_across_lacked_bases(e, expected_body, actual_body));
    }
    if klass == class::REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES {
        return Ok(edited_model_shifted_insertion_bases(e, a, ev, av));
    }
    if klass == class::REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC {
        return Ok(edited_model_seam_insertion_replacement(e, a, ev, av));
    }
    if klass == class::REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC {
        return edited_model_clipped_alleles(e, a, ev, av, false);
    }
    if klass == class::REFSEQ_EDITED_MODEL_HGVS_OFFSET {
        if !blank(ev) || !any_field_reads(e, BAM_EDIT_FIELDS, "OK") {
            return Ok(false);
        }
        let hgvsc = first_field(e, HGVSC_FIELDS);
        let coding = CDS_POSITION_FIELDS
            .iter()
            .any(|field| !blank(e.fields.get(*field).map_or(ABSENT, String::as_str)));
        return Ok(coding && hgvsc.ends_with("dup") && NON_ZERO_INTEGER_WHOLE_RE.is_match(av));
    }
    if klass == class::CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN {
        // The row predicate reads the Consequence whole. The IMPACT differs only where the
        // candidate writes no term in the place of `feature_elongation`, each side's being
        // the one its terms give.
        let Some(a) = a else {
            return Ok(false);
        };
        return Ok(symbolic_loss_named_a_gain_pair(e, a)?
            && (field_named(name, CONSEQUENCE_FIELDS)
                || term_pair_following_cell(e, a, name, ev, av)));
    }
    if klass == class::REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE {
        return Ok(a.is_some_and(|a| intron_span_coding_reference_cell(e, a, name, ev, av)));
    }
    if klass == class::REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS {
        return Ok(a.is_some_and(|a| deleted_span_coding_reference_cell(e, a, name, ev, av)));
    }
    if klass == class::REFSEQ_EDITED_MODEL_INSERTED_SPAN_CODONS {
        let (Some(a), Some(sequences)) = (a, side.downstream_sequences.as_ref()) else {
            return Ok(false);
        };
        return inserted_span_codons_cell(e, a, name, ev, av, sequences);
    }
    if klass == class::HGVS_INTRON_OFFSET_ABUTTING_INSERTION {
        return Ok(intron_offset_abutting_insertion(e, a, name, ev, av));
    }
    if klass == class::HGVS_SHIFTED_INSERTION_CLAMPED_END {
        let offset = first_field(e, HGVS_OFFSET_FIELDS);
        if blank(ev) || blank(av) || blank(offset) || !insertion_row(e) || !positive_shift(e)? {
            return Ok(false);
        }
        // The shift carried the insertion from before the first base to after base k, so
        // its cDNA position starts at the shift's length.
        let shift = first_field(e, SHIFT_LENGTH_FIELDS);
        let cdna_start = pyfmt::partition(first_field(e, CDNA_START_FIELDS), "-").0;
        if cdna_start != shift {
            return Ok(false);
        }
        let (expected_name, _, expected_body) = pyfmt::partition(ev, ":");
        let (actual_name, _, actual_body) = pyfmt::partition(av, ":");
        let clamped = clamped_insertion_position(expected_body);
        let shifted = shifted_insertion_end(actual_body);
        let (
            Some((clamped_numbering, clamped_position, clamped_bases, stripped_count)),
            Some((shifted_numbering, shifted_position, shifted_bases, shifted_count)),
        ) = (clamped, shifted)
        else {
            return Ok(false);
        };
        if expected_name != actual_name
            || (clamped_numbering, clamped_position) != (shifted_numbering, shifted_position)
        {
            return Ok(false);
        }
        let allele_length = pyfmt::len(&e.allele) as i64;
        if shifted_count != allele_length {
            return Ok(false);
        }
        let written = pyfmt::len(&clamped_bases) as i64;
        if written.saturating_add(stripped_count) != allele_length {
            // On a model that carries an edit the reference's `ins` has the row's allele
            // written back whole over the bases the clip left (:1508-1514).
            return Ok(stripped_count == 1
                && written == allele_length
                && shifted_bases.is_none()
                && any_field_reads(e, BAM_EDIT_FIELDS, "OK")
                && allele_readings(e).contains(&clamped_bases));
        }
        return Ok(shifted_bases.is_none_or(|bases| bases.starts_with(&clamped_bases)));
    }
    if klass == class::HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE {
        return shifted_indel_written_span_gate(e, a, name, ev, av);
    }
    if klass == class::HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE {
        return shifted_indel_written_span_protein_gate(
            e,
            a,
            ev,
            av,
            side.transcript_edges.as_ref(),
        );
    }
    if klass == class::HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST {
        return shifted_indel_written_span_start_test(e, a, ev, av);
    }
    if klass == class::HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE {
        return start_codon_test_past_the_shifted_allele_rows(e, a, ev, av);
    }
    if klass == class::HGVS_SHIFTED_DELETION_CLAMPED_SPAN {
        if field_named(name, HGVS_INTRON_START_OFFSET_FIELDS) {
            let Some(a) = a else {
                return Ok(false);
            };
            return shifted_deletion_start_offsets(e, a, ev, av);
        }
        return shifted_deletion_clamped_span(e, ev, av);
    }
    if klass == class::HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE {
        let edges = side.transcript_edges.as_ref();
        if field_named(name, HGVS_INTRON_OFFSET_FIELDS) {
            return clipped_notation_intron_offsets(e, a, name, ev, av, edges);
        }
        return clipped_reference_against_whole_alternate(e, ev, av, edges);
    }
    if klass == class::HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END {
        return added_bases_left_at_the_span_end(e, a, name);
    }
    if klass == class::HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS {
        let edges = side.transcript_edges.as_ref();
        return insertion_left_where_the_shift_stops(e, a, name, edges);
    }
    if klass == class::HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT {
        return added_bases_written_outside_the_transcript(e, a, name, ev, av);
    }
    if klass == class::HGVSP_STOP_CODON_PROTEIN_COMPARISON {
        let flags = field_of(Some(e), FLAGS_FIELDS);
        if stop_unchanged_where_the_reference_writes_none(e, &flags, ev, av)
            || stop_codon_protein_comparison(&e.consequence, &flags, ev, av)?
        {
            return Ok(true);
        }
        let Some(sequences) = side.downstream_sequences.as_ref() else {
            return Ok(false);
        };
        if stop_unchanged_where_inserted_residues_keep_the_protein(e, &flags, ev, av, sequences)? {
            return Ok(true);
        }
        return extension_counted_from_a_stop_the_peptide_holds(e, ev, av, sequences);
    }
    if klass == class::HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE {
        let rule = LandedRule::HgvspStopNamedOnIncompleteCodingSequence;
        if !side.asked_at(e, |at| Ok(stands_on_the_record(rule, at)))? {
            return Ok(false);
        }
        let flags = field_of(Some(e), FLAGS_FIELDS);
        if stop_named_on_incomplete_coding_sequence(&e.consequence, &flags, ev, av) {
            return Ok(true);
        }
        return match side.partial_codon_models.as_ref() {
            Some(models) => {
                hgvsp_stop::stop_named_on_a_transcript_of_a_partial_last_codon(e, ev, av, models)
            }
            None => Ok(false),
        };
    }
    if klass == class::HGVSP_FRAMESHIFT_COUNT_PAST_A_STOP_READ_AS_A_RESIDUE {
        let (Some(a), Some(sequences)) = (a, side.downstream_sequences.as_ref()) else {
            return Ok(false);
        };
        return frameshift_count_past_a_stop_read_as_a_residue(e, a, ev, av, sequences);
    }
    if klass == class::HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE {
        let candidate_terms = a.map(|a| a.consequence.as_str());
        return Ok(frameshift_with_no_changed_residue(
            &e.consequence,
            candidate_terms,
            ev,
            av,
        ));
    }
    if klass == class::HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE {
        let Some(a) = a else {
            return Ok(false);
        };
        let models = side.partial_codon_models.as_ref();
        return frameshift_count_at_the_first_residue(e, a, ev, av, models);
    }
    if klass == class::HGVSP_INTRON_END_INSERTION {
        return hgvsp_intron_end_insertion(e, a, ev, av, side.intron_ends.as_ref());
    }
    if klass == class::HGVSP_INTRON_START_INSERTION {
        return hgvsp_intron_start_insertion(e, a, ev, av, side.intron_ends.as_ref());
    }
    if klass == class::SHIFT_INTRON_TERM_UNSHIFTED_SPAN {
        if blank(ev) || blank(av) || !positive_shift(e)? {
            return Ok(false);
        }
        let terms_of = |value: &str| -> BTreeSet<String> {
            consequence_set(value)
                .split(',')
                .map(str::to_string)
                .collect()
        };
        let (expected_terms, actual_terms) = (terms_of(ev), terms_of(av));
        let mut with_intron = actual_terms.clone();
        with_intron.insert(INTRON_TERM.to_string());
        if actual_terms.contains(INTRON_TERM) || expected_terms != with_intron {
            return Ok(false);
        }
        let hgvsc = first_field(e, HGVSC_FIELDS);
        if blank(hgvsc) {
            return Ok(insertion_row(e) && actual_terms.contains(SPLICE_ACCEPTOR_TERM));
        }
        let Some(offsets) = acceptor_site_offsets(pyfmt::partition(hgvsc, ":").2) else {
            return Ok(false);
        };
        if !offsets.is_empty() {
            return Ok(true);
        }
        // A notation without an intron offset is numbered over the exons alone, so a
        // deletion over a whole intron reads as exonic too; such a span covers both splice
        // sites, and the row carries intron_variant from its shifted span on both engines.
        return Ok(!SPLICE_SITE_TERMS
            .iter()
            .any(|term| actual_terms.contains(*term)));
    }
    if klass == class::FREQUENCY_FILTER_ZERO_VALUE {
        return Ok(freqs_differ_by_zero_read_as_na(ev, av));
    }
    if klass == class::REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET {
        let (Some(a), Some(edits)) = (a, side.length_edits.as_ref()) else {
            return Ok(false);
        };
        return failed_edit_offset_cell(e, a, ev, av, edits, context.declared);
    }
    if klass == class::USE_TRANSCRIPT_REF_CACHE_CDNA {
        return Ok(a.is_some_and(|a| failed_edit_reference_cell(e, a, name, ev, av)));
    }
    if klass == class::HGNC_ID_COPY {
        return Ok(hgnc_id_copy_pair(e, a, ev, av));
    }
    if klass == class::HGVS_SHIFT_SHARED_ACROSS_ALLELES {
        let (Some(a), Some(sequences)) = (a, side.downstream_sequences.as_ref()) else {
            return Ok(false);
        };
        return later_allele_shared_shift(e, a, name, record_rows, sequences);
    }
    if klass == class::BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES {
        let (Some(a), Some(sequences)) = (a, side.downstream_sequences.as_ref()) else {
            return Ok(false);
        };
        return blosum62_later_allele_shared_shift(e, a, ev, av, record_rows, sequences);
    }
    if tables.in_transcript_reading_classes(klass) {
        let (Some(a), Some(sequences)) = (a, side.downstream_sequences.as_ref()) else {
            return Ok(false);
        };
        if field_named(name, NMD_FIELDS) {
            return nmd_reads_shifted_pair(
                e,
                a,
                expected_value,
                actual_value,
                sequences,
                record_rows,
            );
        }
        // The values as the entries carry them: an empty one is the plugin's, a `-` or an
        // absent one is no value.
        let decided = side.asked_at(e, |at| {
            downstream_class_at(
                e,
                a,
                name,
                expected_value,
                actual_value,
                sequences,
                record_rows,
                at,
            )
        })?;
        return Ok(decided.is_some_and(|(first, second)| {
            first == klass
                && second.is_none_or(|second| context.declared.contains(&second.as_str()))
        }));
    }
    if klass == class::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN {
        let Some(a) = a else {
            return Ok(false);
        };
        if !side.asked_at(e, |at| splice_term_past_the_record_span_pair_at(e, a, at))? {
            return Ok(false);
        }
        return Ok(term_pair_following_cell(e, a, name, ev, av));
    }
    if klass == class::PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR {
        let Some(a) = a else {
            return Ok(false);
        };
        let taken = side.asked_at(e, |at| {
            for pair in STOP_CODON_PLUGIN_VALUE_PAIRS {
                if pair(e, a, at)? {
                    return Ok(true);
                }
            }
            Ok(false)
        })?;
        return Ok(taken && term_pair_following_cell(e, a, name, ev, av));
    }
    if klass == class::TOTAL_LENGTH_JSON_RANGE_END {
        return Ok(a.is_some_and(|a| total_length_range_end_pair(e, a, name, ev, av)));
    }
    if klass == class::ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT {
        return Ok(a.is_some_and(|a| {
            a.allele == e.allele
                && allele_number_place_pair(
                    written_allele_list(e),
                    written_allele_list(a),
                    &e.allele,
                    ev,
                    av,
                )
        }));
    }
    Ok(false)
}

/// `pattern` held to the whole text, as `fullmatch` holds it.
pub(crate) fn whole_text(pattern: &str) -> regex::Regex {
    regex::Regex::new(&format!(r"\A(?:{pattern})\z")).expect("the pattern compiles")
}

static MAXENTSCAN_DONOR_TEXT_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(MAXENTSCAN_DONOR_TEXT_RE_PATTERN));
static MAXENTSCAN_ACCEPTOR_TEXT_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(MAXENTSCAN_ACCEPTOR_TEXT_RE_PATTERN));
static INTEGER_WHOLE_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(INTEGER_RE_PATTERN));
static NON_ZERO_INTEGER_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));

/// `int(text)` as its decimal text; a text the interpreter does not convert ends the run
/// as its `ValueError` does.
pub(crate) fn whole_number_text(text: &str) -> String {
    pyfmt::int_checked(text).unwrap_or_else(|error| crate::uncaught!("{error}"))
}

/// `int(text)` as a whole number of any length `i128` holds, the nearest one past it; a
/// text the interpreter does not convert ends the run as `whole_number_text` ends it.
pub(crate) fn whole_number(text: &str) -> i128 {
    let decimal = whole_number_text(text);
    decimal.parse().unwrap_or(if decimal.starts_with('-') {
        i128::MIN
    } else {
        i128::MAX
    })
}

/// Whether `text` is one or more of the ASCII `letters` and nothing else.
pub(crate) fn spelt_in(text: &str, letters: &[u8]) -> bool {
    !text.is_empty() && text.bytes().all(|byte| letters.contains(&byte))
}

/// A whole number as `i64`, the nearest one where it lies past either end.
pub(crate) fn nearest_i64(number: i128) -> i64 {
    i64::try_from(number).unwrap_or(if number < 0 { i64::MIN } else { i64::MAX })
}

/// The value of the first of the fields `names` the entry carries, `ABSENT` for none.
pub(crate) fn first_field<'a>(e: &'a Entry, names: &[&str]) -> &'a str {
    names
        .iter()
        .find_map(|name| e.fields.get(*name))
        .map_or(ABSENT, String::as_str)
}

/// The record's written allele list as the entry carries it, under the name its format gives
/// the field; `None` for an entry that carries none.
fn written_allele_list(e: &Entry) -> Option<&str> {
    WRITTEN_ALLELE_LIST_FIELDS
        .iter()
        .find_map(|name| e.fields.get(*name))
        .map(String::as_str)
}

/// Whether one of the fields `names` of the entry reads `value`.
fn any_field_reads(e: &Entry, names: &[&str], value: &str) -> bool {
    names
        .iter()
        .any(|name| e.fields.get(*name).is_some_and(|held| held == value))
}

/// Whether a set of terms, joined by `,` as `consequence_set` joins it, holds one of `terms`.
fn carries_one_of(joined: &str, terms: &[&str]) -> bool {
    joined.split(',').any(|held| terms.contains(&held))
}

/// Whether the whole number `offset` writes is, its sign set aside, greater than the
/// positive `length` and no multiple of it; the number may be of any length `int` converts.
fn shifted_past_and_off_a_multiple(offset: &str, length: i64) -> bool {
    let decimal = match pyfmt::int_checked(offset) {
        Ok(decimal) => decimal,
        Err(error) => crate::uncaught!("{error}"),
    };
    let digits = decimal.trim_start_matches('-');
    let length = u128::from(length.unsigned_abs());
    // Twenty-one digits are past every `i64`; a shorter run fits `u128`.
    let greater = digits.len() > 20 || digits.parse::<u128>().is_ok_and(|shift| shift > length);
    let remainder = digits.bytes().fold(0u128, |remainder, digit| {
        (remainder * 10 + u128::from(digit - b'0')) % length
    });
    greater && remainder != 0
}
