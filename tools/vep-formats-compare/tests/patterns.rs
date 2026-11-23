// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The patterns the tool holds in public statics compile. A pattern in a private static is
//! out of a test file's reach and is compiled by the functions of its module that match it.

use std::sync::LazyLock;

#[test]
fn every_pattern_in_a_public_static_compiles() {
    LazyLock::force(&vep_formats_compare::compare::MOTIF_SCORE_TEXT);
    LazyLock::force(&vep_formats_compare::parse::recoder::INDEL_NOTATION);
    LazyLock::force(&vep_formats_compare::consts::RECORD_LEVEL_COLUMN_RE);
    LazyLock::force(&vep_formats_compare::consts::MISSPELT_CONSERVATION_KEY);
    LazyLock::force(&vep_formats_compare::rules::TER_COUNT_RE);
    LazyLock::force(&vep_formats_compare::rules::gwas::GWAS_LIST_SPLIT_RE);
    LazyLock::force(&vep_formats_compare::rules::gwas::GWAS_RS_ID_RE);
    LazyLock::force(&vep_formats_compare::rules::utrannotator::CAP_DISTANCE_RE);
    LazyLock::force(&vep_formats_compare::rules::utrannotator::UTR_HGVSC_INSERTION_RE);
    LazyLock::force(&vep_formats_compare::rules::utrannotator::UTR_HGVSC_DUPLICATION_RE);
    LazyLock::force(&vep_formats_compare::rules::MAXENTSCAN_DONOR_TEXT_RE);
    LazyLock::force(&vep_formats_compare::rules::MAXENTSCAN_ACCEPTOR_TEXT_RE);
    LazyLock::force(&vep_formats_compare::rules::nearest_exon::NEAREST_EXON_ITEM_RE);
    LazyLock::force(&vep_formats_compare::rules::nearest_exon::NEAREST_EXON_LIST_RE);
    LazyLock::force(&vep_formats_compare::warnings::FORKED_CHILD_RELAY_RE);
    LazyLock::force(&vep_formats_compare::rules::CDNA_COORDINATE_RE);
    LazyLock::force(&vep_formats_compare::rules::edited_model::CDNA_COORDINATE_PAIR_RE);
    LazyLock::force(&vep_formats_compare::consts::NON_ZERO_INTEGER_RE);
    LazyLock::force(&vep_formats_compare::rules::edited_model::EDITED_INSERTION_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::term_pairs::CLAMPED_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::term_pairs::INSERTION_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::term_pairs::SHIFTED_DUPLICATION_RE);
    LazyLock::force(&vep_formats_compare::consts::INTEGER_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::PROTEIN_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::STOP_CODON_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::STOP_CODON_FLANK_INSERTION_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::STOP_CODON_RANGE_END_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::CODING_COORDINATE_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::ACCEPTOR_INTRON_COORDINATE_RE);
    LazyLock::force(&vep_formats_compare::rules::shifted_indel::DONOR_INTRON_COORDINATE_RE);
    LazyLock::force(&vep_formats_compare::consts::DELETION_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::clipped::COORDINATE_PARTS_AS_READ_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::CLIPPED_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::rules::added_bases::ONE_PAST_COORDINATE_RE);
    LazyLock::force(&vep_formats_compare::rules::clipped::ONE_BEFORE_COORDINATE_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_UNCHANGED_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_EXTENSION_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_REPLACED_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_INSERTION_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_RANGE_REPLACED_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_FRAMESHIFT_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_DELETED_EXTENSION_RE);
    LazyLock::force(&vep_formats_compare::rules::hgvsp_stop::STOP_NAMED_ALONE_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::INTRON_END_COPY_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::INTRON_END_DUPLICATION_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::INTRON_END_INSERTION_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::INTRON_START_INSERTION_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::VERSION_SUFFIX_RE);
    LazyLock::force(&vep_formats_compare::rules::intron_ends::CACHE_SHARD_NAME_RE);
    LazyLock::force(&vep_formats_compare::consts::INDEL_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::consts::INTRON_OFFSET_RE);
    LazyLock::force(&vep_formats_compare::rules::downstream::HGNC_ID_RE);
    LazyLock::force(&vep_formats_compare::rules::downstream::PEPTIDE_PAIR_RE);
    LazyLock::force(&vep_formats_compare::rules::edited_model::CLIPPED_ALLELES_NOTATION_RE);
    LazyLock::force(&vep_formats_compare::consts::POSITION_RE);
    LazyLock::force(&vep_formats_compare::self_parity::POSITION_PART_RE);
    LazyLock::force(&vep_formats_compare::self_parity::TOOL_SCORE_RE);
    LazyLock::force(&vep_formats_compare::self_parity::DOMAIN_RE);
    LazyLock::force(&vep_formats_compare::self_parity::PERL_NUMBER_RE);
    LazyLock::force(&vep_formats_compare::warnings::WARNING_LINE_RE);
    LazyLock::force(&vep_formats_compare::warnings::WARNING_OTHER_RE);
    LazyLock::force(&vep_formats_compare::warnings::SKIPPED_LINE_RE);
    LazyLock::force(&vep_formats_compare::warnings::PATH_RE);
    LazyLock::force(&vep_formats_compare::warnings::UNDEFINED_ALLELE_WARNING);
    LazyLock::force(&vep_formats_compare::warnings::EXCEPTION_RULE);
    LazyLock::force(&vep_formats_compare::shared::UPLOADED_ALLELES_RE);
    LazyLock::force(&vep_formats_compare::shared::UPLOADED_POSITION_RE);
    LazyLock::force(&vep_formats_compare::shared::_VEP_RS_POS_NUMBER);
    LazyLock::force(&vep_formats_compare::rules::structural::_BRACKET_MATE_RE);
    LazyLock::force(&vep_formats_compare::rules::structural::_VCF44_SYMBOLIC);
    LazyLock::force(&vep_formats_compare::rules::structural::_PERL_NUMBER);
    LazyLock::force(&vep_formats_compare::rules::structural::_PERL_NAN);
    LazyLock::force(&vep_formats_compare::rules::structural::_PERL_INFINITY);
}
