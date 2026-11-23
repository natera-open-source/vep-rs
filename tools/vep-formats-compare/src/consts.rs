// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Constants several modules read.

use std::sync::LazyLock;

/// Membership only.
pub const REGULATORY_KINDS: &[&str] = &["motif_feature", "regulatory_feature"];

pub const STRAND_FIELDS: &[&str] = &["STRAND", "strand"];

/// The name `--run` gives the manifest's own run, as the golden harness names it; a flag set
/// is named by its `name`, which the harness holds to another word.
pub const OWN_RUN: &str = "default";

/// The optional key of a `divergences` item: the candidate's consequence sets on the key in
/// a run that has a reference sequence, `vep_rs_consequence_sets` being those of a run
/// without one.
pub const WITH_REFERENCE_SETS: &str = "vep_rs_consequence_sets_with_reference_sequence";

/// The optional key of a `field_divergences` item: the count of the rows of its key and
/// consequence set on which both engines print one value of the field, the rows of other
/// records that share the key with the documented one.
pub const ROWS_ALIKE: &str = "rows_alike";

/// The value of `--reference-sequence` for a candidate run that had a reference sequence.
pub const REFERENCE_PRESENT: &str = "present";

pub const EXAMPLES_PER_GROUP: i64 = 4;

/// OutputFactory/JSON.pm:98-112 at release/116.3.
pub const JSON_RENAME: &[(&str, &str)] = &[
    ("consequence", "consequence_terms"),
    ("gene", "gene_id"),
    ("allele", "variant_allele"),
    ("symbol", "gene_symbol"),
    ("symbol_source", "gene_symbol_source"),
    ("overlapbp", "bp_overlap"),
    ("overlappc", "percentage_overlap"),
    ("refseq", "refseq_transcript_ids"),
    ("ensp", "protein_id"),
    ("chr", "seq_region_name"),
    ("variation_name", "id"),
    ("sv", "colocated_structural_variants"),
    ("clinical_impact", "clinvar_somatic_classification"),
];

pub const JSON_CONSEQUENCE_LISTS: &[(&str, &str)] = &[
    ("transcript_consequences", "transcript"),
    ("intergenic_consequences", "intergenic"),
    ("regulatory_feature_consequences", "regulatory_feature"),
    ("motif_feature_consequences", "motif_feature"),
];

/// Columns the default format prints on every row of a record while the JSON format
/// carries them once per record, in colocated_variants or at the top level, or not at all:
/// every key VariationFeature_to_output_hash sets (OutputFactory.pm:883-1023: SV, VARIANT_CLASS,
/// IND, ZYG, MINIMISED, AMBIGUITY, NEAREST and CHECK_REF) is copied once to the record object
/// (JSON.pm:262), never into a consequence entry; the JSON factory's add_colocated_variant_info
/// returns its hash unchanged (JSON.pm:493-495), so FREQS (OutputFactory.pm:1149-1163) has no
/// JSON key, and JSON.pm:419-426 copies the AF keys alone out of the per-allele hash that
/// carries MAX_AF and MAX_AF_POPS.
///
/// Membership only.
pub const RECORD_LEVEL_COLUMNS: &[&str] = &[
    "AMBIGUITY",
    "CHECK_REF",
    "CLINVAR_SOMATIC_CLASSIFICATION",
    "CLIN_SIG",
    "Existing_variation",
    "FREQS",
    "IND",
    "Location",
    "MAX_AF",
    "MAX_AF_POPS",
    "MINIMISED",
    "NEAREST",
    "PHENO",
    "PUBMED",
    "SOMATIC",
    "SV",
    "Uploaded_variation",
    "VARIANT_CLASS",
    "VAR_SYNONYMS",
    "ZYG",
];

pub const RECORD_LEVEL_COLUMN_RE_PATTERN: &str = r"(^|_)AF(?:\n?\z)";
pub static RECORD_LEVEL_COLUMN_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(RECORD_LEVEL_COLUMN_RE_PATTERN).expect("the pattern compiles")
});

pub const MISSPELT_CONSERVATION_KEY_PATTERN: &str = r"phlyoP(\d+)way";
pub static MISSPELT_CONSERVATION_KEY: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(MISSPELT_CONSERVATION_KEY_PATTERN).expect("the pattern compiles")
});

pub const CONSERVATION_KEY_SPELLING: &str = "phyloP\\1way";

pub const HGVSP_FIELDS: &[&str] = &["HGVSp", "hgvsp"];

pub const CONSERVATION_FIELDS: &[&str] = &["Conservation", "conservation"];

pub const CADD_PREFIX: &str = "cadd_";

pub const SPLICEVAULT_PREFIXES: &[&str] = &["splicevault_", "splicevault"];

pub const SPLICEVAULT_JSON_FIELD: &str = "splicevault";

pub const UTRANNOTATOR_PREFIX: &str = "5utr_";

pub const CDNA_START_FIELDS: &[&str] = &["cDNA_position", "cdna_start"];

pub const CDNA_END_FIELD: &str = "cdna_end";

pub const UTRANNOTATOR_ANNOTATION_FIELDS: &[&str] = &["5UTR_annotation", "5utr_annotation"];

/// A plain insertion's class word; a mobile element's insertion has a word of its own. A
/// stated length that is no multiple of three shifts the frame whatever is inserted
/// (UTRAnnotator.pm:889, :923), so the candidate writes the frameshift alone on such a row:
/// each effect with the six keys that read no inserted base, and `alt_type` as `uORF`
/// where a stop codon past the point settles it.
pub const STRUCTURAL_INSERTION_WORD: &str = "insertion";

pub const LOEUF_FIELDS: &[&str] = &["LOEUF", "loeuf"];

pub const MAXENTSCAN_PREFIX: &str = "maxentscan_";

pub const AMINO_ACIDS_FIELDS: &[&str] = &["Amino_acids"];

pub const REVEL_FIELDS: &[&str] = &["REVEL"];

pub const LOF_INFO_FIELDS: &[&str] = &["LoF_info"];

pub const NEAREST_EXON_JB_FIELDS: &[&str] = &["NearestExonJB"];

pub const DISTANCE_FIELDS: &[&str] = &["DISTANCE"];

pub const SPLICEAI_PREFIXES: &[&str] = &["SpliceAI_", "spliceai"];

pub const TABIX_PLUGIN_FIELDS: &[&str] = &[
    "am_",
    "alphamissense",
    "AVI_",
    "BayesDel",
    "cadd_",
    "ClinPred",
    "ada_score",
    "rf_score",
    "LOEUF",
    "loeuf",
    "OpenTargets_",
    "opentargets",
    "PrimateAI",
    "REVEL",
    "SpliceAI_",
    "spliceai",
    "splicevault_",
    "splicevault",
    "gnomad_",
];

pub const HGVSC_FIELDS: &[&str] = &["HGVSc", "hgvsc"];

pub const BAM_EDIT_FIELDS: &[&str] = &["BAM_EDIT", "bam_edit"];

pub const REFSEQ_OFFSET_FIELDS: &[&str] = &["REFSEQ_OFFSET", "refseq_offset"];

pub const HGVS_OFFSET_FIELDS: &[&str] = &["HGVS_OFFSET", "hgvs_offset"];

pub const CDS_POSITION_FIELDS: &[&str] = &["CDS_position", "cds_start", "cds_end"];

pub const NON_ZERO_INTEGER_RE_PATTERN: &str = r"-?[1-9]\d*";
pub static NON_ZERO_INTEGER_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(NON_ZERO_INTEGER_RE_PATTERN).expect("the pattern compiles"));

pub const SHIFT_LENGTH_FIELDS: &[&str] = &["SHIFT_LENGTH", "shift_length"];

pub const HGVS_INTRON_START_OFFSET_FIELDS: &[&str] = &["HGVS_IntronStartOffset"];

pub const HGVS_INTRON_OFFSET_FIELDS: &[&str] = &["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"];

pub const OUTSIDE_TRANSCRIPT_TERMS: &[&str] = &["downstream_gene_variant", "upstream_gene_variant"];

pub const INTEGER_RE_PATTERN: &str = r"-?\d+";
pub static INTEGER_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(INTEGER_RE_PATTERN).expect("the pattern compiles"));

pub const PROTEIN_POSITION_FIELDS: &[&str] = &["Protein_position", "protein_start", "protein_end"];

pub const SHIFTED_DELETION_SPAN_FIELDS: &[&str] = &["HGVSc", "hgvsc", "HGVS_IntronStartOffset"];

pub const DELETION_NOTATION_RE_PATTERN: &str =
    r"([cn])\.(\*?-?\d+(?:[+-]\d+)?)(?:_(\*?-?\d+(?:[+-]\d+)?))?del";
pub static DELETION_NOTATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(DELETION_NOTATION_RE_PATTERN).expect("the pattern compiles")
});

pub const CLIPPED_REFERENCE_FIELDS: &[&str] = &[
    "HGVSc",
    "hgvsc",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];

pub const ADDED_BASES_FIELDS: &[&str] = &[
    "HGVSc",
    "hgvsc",
    "HGVS_OFFSET",
    "hgvs_offset",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];

pub const SHIFT_BOUND_FIELDS: &[&str] = &[
    "HGVSc",
    "hgvsc",
    "HGVS_OFFSET",
    "hgvs_offset",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
    "HGVSp",
    "hgvsp",
];

pub const FLAGS_FIELDS: &[&str] = &["FLAGS", "flags"];

pub const INDEL_NOTATION_RE_PATTERN: &str =
    r"[cn]\.(\*?-?\d+(?:[+-]\d+)?)(?:_(\*?-?\d+(?:[+-]\d+)?))?(dup|del|ins[ACGT]+)";
pub static INDEL_NOTATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(INDEL_NOTATION_RE_PATTERN).expect("the pattern compiles"));

/// The intron offset a coordinate ends with: the sign follows the position's digits, so the
/// leading minus of a 5' UTR position (`c.-3`) is none.
pub const INTRON_OFFSET_RE_PATTERN: &str = r"(?<=\d)[+-]\d+(?=\n?\z)";
pub static INTRON_OFFSET_RE: LazyLock<fancy_regex::Regex> = LazyLock::new(|| {
    fancy_regex::Regex::new(INTRON_OFFSET_RE_PATTERN).expect("the pattern compiles")
});

pub const CDNA_POSITION_FIELDS: &[&str] = &["cDNA_position", "cdna_start", "cdna_end"];

pub const USED_REF_FIELDS: &[&str] = &["USED_REF", "used_ref"];

pub const FREQS_FIELDS: &[&str] = &["FREQS", "freqs"];

pub const HGNC_ID_FIELDS: &[&str] = &["HGNC_ID", "hgnc_id"];

pub const DOWNSTREAM_PROTEIN_FIELDS: &[&str] = &["DownstreamProtein", "downstreamprotein"];

pub const NMD_FIELDS: &[&str] = &["NMD", "nmd"];

pub const DOSAGE_SENSITIVITY_FIELDS: &[&str] = &["pHaplo", "pTriplo"];

pub const UTRANNOTATOR_EXISTING_FIELDS: &[&str] = &[
    "Existing_InFrame_oORFs",
    "Existing_OutOfFrame_oORFs",
    "Existing_uORFs",
];

/// `hgvs_shift_shared_across_alleles` beside the candidate JSON cache: three columns of a
/// later allele's row are decided without a manifest. The alleles of a record share one
/// TranscriptVariation per transcript, which keeps the CDS pair and the translation pair once
/// read (BaseTranscriptVariation.pm:252-267, :365-380). The first allele's protein notation
/// reads both again at its 3'-shifted position and leaves them there
/// (TranscriptVariationAllele.pm:1724-1733, :844-858, :2404-2410), and reads the shared
/// reference allele's peptide again at those positions (:1780-1783, :1793). A later allele's
/// row prints the kept pairs (OutputFactory.pm:1746, :1750) and that peptide before its own
/// (:1742, TranscriptVariationAllele.pm:618-630). Its own peptide stands at the record's own
/// position: a run that keeps statistics reads every allele's consequences, and with them its
/// peptide, before the record's first row (OutputFactory.pm:347, :1022; Stats.pm:241;
/// VariationFeature.pm:1000-1002; Utils/VariationEffect.pm:785). The shape
/// (`later_allele_shared_shift`): a differing CDS_position, Protein_position or Amino_acids
/// on the row of a second or later alternate allele of its record (`later_alternate`) whose
/// first allele's reference row on the feature carries a non-zero HGVS_OFFSET
/// (`first_allele_shift`). The candidate prints the record's own CDS pair and the reference
/// that span moved by the first allele's offset, both mapped through the cache's mapper pairs
/// (`shifted_pair_readings`). Each side's Protein_position is the codons of its pair, and
/// the two Amino_acids share their alternate part, each side's reference part the coding
/// sequence's translation over its codons under the transcript's table. A first allele's row,
/// a first allele without an offset, a move of another size, a candidate value that is not the
/// record's own, any other column, a row whose record's alleles no entry states in order and a
/// transcript the cache holds without its mapper pairs stay differences; without the cache the
/// class decides nothing from the two outputs.
pub const SHARED_SHIFT_FIELDS: &[&str] = &[
    "CDS_position",
    "cds_start",
    "cds_end",
    "Protein_position",
    "protein_start",
    "protein_end",
    "Amino_acids",
];

/// ensembl-variation Utils/Constants.pm %OVERLAP_CONSEQUENCES at release/116, SO term to rank.
pub const CONSEQUENCE_RANK: &[(&str, i64)] = &[
    ("transcript_ablation", 1),
    ("splice_acceptor_variant", 2),
    ("splice_donor_variant", 3),
    ("stop_gained", 4),
    ("frameshift_variant", 5),
    ("stop_lost", 6),
    ("start_lost", 7),
    ("transcript_amplification", 8),
    ("feature_elongation", 9),
    ("feature_truncation", 10),
    ("inframe_insertion", 11),
    ("inframe_deletion", 12),
    ("missense_variant", 13),
    ("protein_altering_variant", 14),
    ("splice_donor_5th_base_variant", 15),
    ("splice_region_variant", 16),
    ("splice_donor_region_variant", 17),
    ("splice_polypyrimidine_tract_variant", 18),
    ("incomplete_terminal_codon_variant", 19),
    ("start_retained_variant", 20),
    ("stop_retained_variant", 21),
    ("synonymous_variant", 22),
    ("coding_sequence_variant", 23),
    ("mature_miRNA_variant", 24),
    ("5_prime_UTR_variant", 25),
    ("3_prime_UTR_variant", 26),
    ("non_coding_transcript_exon_variant", 27),
    ("intron_variant", 28),
    ("NMD_transcript_variant", 29),
    ("non_coding_transcript_variant", 30),
    ("coding_transcript_variant", 31),
    ("upstream_gene_variant", 32),
    ("downstream_gene_variant", 33),
    ("TFBS_ablation", 34),
    ("TFBS_amplification", 35),
    ("TF_binding_site_variant", 36),
    ("regulatory_region_ablation", 37),
    ("regulatory_region_amplification", 38),
    ("regulatory_region_variant", 39),
    ("intergenic_variant", 40),
    ("sequence_variant", 41),
];

/// The fields whose text is a JSON object, compared in canonical form.
pub const JSON_OBJECT_FIELDS: &[&str] = &["GA4GH_VRS"];

pub const TERM_PAIR_FOLLOWING_FIELDS: &[&str] = &[
    "IMPACT",
    "impact",
    "DownstreamProtein",
    "ProteinLengthChange",
    "NMD",
    "LoF",
    "LoF_info",
    "LoF_filter",
    "LoF_flags",
    "HGVSp",
    "hgvsp",
];

pub const POSITION_RE_PATTERN: &str = r"(?i)([\p{L}\p{N}_]+?)_pos[iİı]t[iİı]on(?:\n?\z)";
pub static POSITION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(POSITION_RE_PATTERN).expect("the pattern compiles"));

pub const COMPLEMENT: &[(i64, i64)] = &[
    (65, 84),
    (67, 71),
    (71, 67),
    (84, 65),
    (78, 78),
    (97, 116),
    (99, 103),
    (103, 99),
    (116, 97),
    (110, 110),
];
