// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variant types: input variants and co-located known variants.

use std::borrow::Cow;
use std::sync::{Arc, OnceLock};

use crate::consequence::{Consequence, ConsequenceList, TranscriptConsequence};
use crate::coordinate::{normalize_chromosome, Strand};
use crate::sv_routing::{perl_true, sv_so_term_and_class, sv_type_abbreviation, sv_type_string};

/// Classification of variant types based on Sequence Ontology (SO) variant classes.
///
/// Covers both small sequence-level variants (SNV, insertion, deletion, indel,
/// substitution) and structural variants (symbolic VCF alleles like `<DEL>`,
/// `<DUP>`, `<INV>`, BND breakpoints, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum VariantClass {
    Snv,
    Insertion,
    Deletion,
    Indel,
    Substitution,
    SequenceAlteration,

    /// `<DEL>`: large structural deletion
    StructuralDeletion,
    /// `<INS>`: large structural insertion
    StructuralInsertion,
    /// `<DUP>`: duplication
    Duplication,
    /// `<DUP:TANDEM>`: tandem duplication
    TandemDuplication,
    /// `<INV>`: inversion
    Inversion,
    /// `<CNV>`, `<CN0>` to `<CN9>`: copy-number variation
    CopyNumberVariation,
    /// `<INS:ME:*>`: mobile element insertion
    MobileElementInsertion,
    /// `<DEL:ME:*>`: mobile element deletion
    MobileElementDeletion,
    /// BND notation (bracket/dot syntax): chromosome breakpoint
    Translocation,
    /// `<CNV:TR>`: tandem repeat
    TandemRepeat,
    /// `<CPX>` or other unrecognised symbolic SV types
    ComplexStructural,
    /// `<NON_REF>`: a gVCF reference-confidence block. Not a variant allele, so it
    /// is annotated against no feature and a record carrying it is written as one
    /// `intergenic_variant` row spanning the block.
    ReferenceBlock,
}

impl VariantClass {
    /// Return the Sequence Ontology term for this variant class.
    ///
    /// SO terms match Ensembl VEP's `%SO_TERMS` mapping exactly.
    pub fn so_term(self) -> &'static str {
        match self {
            VariantClass::Snv => "SNV",
            VariantClass::Insertion => "insertion",
            VariantClass::Deletion => "deletion",
            VariantClass::Indel => "indel",
            VariantClass::Substitution => "substitution",
            VariantClass::SequenceAlteration => "sequence_alteration",
            // Structural variant terms follow Ensembl VEP `%SO_TERMS`.
            VariantClass::StructuralDeletion => "deletion",
            VariantClass::StructuralInsertion => "insertion",
            VariantClass::Duplication => "duplication",
            VariantClass::TandemDuplication => "tandem_duplication",
            VariantClass::Inversion => "inversion",
            VariantClass::CopyNumberVariation => "copy_number_variation",
            VariantClass::MobileElementInsertion => "mobile_element_insertion",
            VariantClass::MobileElementDeletion => "mobile_element_deletion",
            VariantClass::Translocation => "chromosome_breakpoint",
            VariantClass::TandemRepeat => "tandem_repeat",
            VariantClass::ComplexStructural => "complex_structural_alteration",
            // VEP's `class_SO_term` for a type its table lacks is the type as written.
            VariantClass::ReferenceBlock => "<NON_REF>",
        }
    }

    /// Whether this class is a symbolic structural allele annotated from POS+1 (true for every SV class).
    ///
    /// In VCF, POS is the anchor base preceding the structural event.
    /// Ensembl VEP uses POS+1 as the start coordinate for all SV classes:
    /// span-based SVs (DEL, DUP, INV, CNV, CPX, TR) extend from POS+1 to END,
    /// while point-based SVs (INS, ME, BND) use POS+1 as a single position.
    pub fn is_span_type(self) -> bool {
        matches!(
            self,
            VariantClass::StructuralDeletion
                | VariantClass::StructuralInsertion
                | VariantClass::Duplication
                | VariantClass::TandemDuplication
                | VariantClass::Inversion
                | VariantClass::CopyNumberVariation
                | VariantClass::MobileElementInsertion
                | VariantClass::MobileElementDeletion
                | VariantClass::Translocation
                | VariantClass::ComplexStructural
                | VariantClass::TandemRepeat
                | VariantClass::ReferenceBlock
        )
    }

    /// Whether this variant class represents a structural variant (symbolic allele).
    pub fn is_structural(self) -> bool {
        matches!(
            self,
            VariantClass::StructuralDeletion
                | VariantClass::StructuralInsertion
                | VariantClass::Duplication
                | VariantClass::TandemDuplication
                | VariantClass::Inversion
                | VariantClass::CopyNumberVariation
                | VariantClass::MobileElementInsertion
                | VariantClass::MobileElementDeletion
                | VariantClass::Translocation
                | VariantClass::TandemRepeat
                | VariantClass::ComplexStructural
                | VariantClass::ReferenceBlock
        )
    }
}

/// A single-allele input variant to be annotated.
///
/// Represents one allele of one VCF record after multi-allelic splitting.
/// Coordinates are 1-based inclusive (VEP convention). For insertions, the
/// ref allele is `"-"` and `end < start`. Carries annotation results
/// (transcript consequences, co-located variants, plugin data) that are
/// populated during the pipeline.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct InputVariant {
    /// Chromosome name, normalized for cache lookup (e.g., "21", "X", "MT").
    pub chr: String,
    /// The chromosome name the record prints (Location, `seq_region_name`): the
    /// input's spelling (`chr21`, `chrX`, `chrM`), which Ensembl VEP keeps, or
    /// `MT` for an input `M` that reaches the cache only through the rename of
    /// Parser.pm:794-796.
    pub original_chr: String,
    /// The input's spelling where `original_chr` carries the renamed `MT`:
    /// Ensembl VEP names a record (Parser.pm:462) before it renames `M` to `MT`,
    /// so the names built from the input keep `M`. `None` when the two agree.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub input_chr: Option<Arc<str>>,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (inclusive).
    pub end: u64,
    /// Strand.
    pub strand: Strand,
    /// Reference allele (e.g., b"A", b"AC", b"-" for insertions).
    pub ref_allele: Vec<u8>,
    /// Alternate alleles.
    pub alt_alleles: Vec<Vec<u8>>,
    /// Allele string in VEP format (e.g., "A/G", "AC/-", "-/TT").
    pub allele_string: String,
    /// The specific alternate allele this variant represents
    /// (for multi-allelic sites split into separate records).
    pub allele_index: usize,
    /// The ALLELE_NUM Ensembl VEP prints when it is not `allele_index + 1`: the
    /// position of the allele in a breakpoint record's allele string, whose reference
    /// base is part 1 (`N/<BND>` prints 2; StructuralVariationOverlap.pm:64-77), with 0
    /// for the record's own breakend and 1 upward for its mates (`:89-101`); or the
    /// allele's number among the alleles Ensembl VEP kept, counted from 1, of a
    /// sequence record whose earlier ALT was dropped as a repeat or the reference
    /// (VariationFeatureOverlap.pm:459-491; `T` of `C G,G,T` prints 2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allele_num: Option<u32>,
    /// The record's name as Ensembl VEP's parser stores it before `validate_vf`: the
    /// first `;`-joined ID of a VCF record ([`first_vcf_id`]) or the whole ID column
    /// of an Ensembl-format line, `.` and `0` kept as written; `None` when the column
    /// is absent or empty.
    pub id: Option<String>,
    /// Variant class.
    pub variant_class: VariantClass,
    /// Consequence annotations computed during annotation phase.
    pub transcript_consequences: Vec<TranscriptConsequence>,
    /// The allele overlaps no transcript and its rows are regulatory and motif
    /// rows alone, so it carries an `intergenic_variant` row beside them: Perl's
    /// intergenic row is decided by transcripts alone
    /// (VariationFeature.pm:826-839; StructuralVariationFeature.pm:674-690).
    #[serde(default)]
    pub intergenic_overlap: bool,
    /// Most severe consequence across all transcripts.
    pub most_severe_consequence: Option<Consequence>,
    /// The known variants co-located with the whole input record, in the cache's
    /// order with a name listed once: Ensembl VEP matches the record's variation
    /// feature, every allele at once, and keys the matches on it
    /// (`AnnotationType/Variation.pm` `compare_existing`), so an allele of a
    /// multi-allelic record carries the matches of its siblings too, each telling
    /// through [`ColocatedVariant::applies_to_allele`] which rows it prints on.
    pub colocated_variants: Vec<ColocatedVariant>,
    /// The HGVS genomic notation of this allele (`HGVSg`), computed under
    /// `--hgvsg` with a reference FASTA; printed on every row of the allele.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hgvsg: Option<String>,
    /// The `NEAREST` values under `--nearest`: the identifiers of the
    /// protein-coding transcription start sites nearest the record, each once,
    /// printed on every row of the record; `None` without the flag, empty on a
    /// chromosome with no such site.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest: Option<Vec<String>>,
    /// The SPDI notation of this allele (`SPDI`), computed under `--spdi` or
    /// `--ga4gh_vrs` from the record's alleles and span as written; a row prints
    /// it when its `Allele` is the allele it describes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spdi: Option<String>,
    /// The GA4GH VRS allele of this allele (`GA4GH_VRS`), computed under
    /// `--ga4gh_vrs` from the SPDI named by the chromosome's RefSeq accession;
    /// `None` when the synonyms table pairs the chromosome with no `NC_` name.
    /// Printed on the rows that print the SPDI.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ga4gh_vrs: Option<crate::spdi::VrsAllele>,
    /// The notations of the record's other alternate alleles, for a row whose
    /// `Allele` a genomic shift rotated onto one of them (`--shift_genomic 1`);
    /// empty for a bi-allelic record.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sibling_spdi: Vec<crate::spdi::SiblingSpdi>,
    /// Ensembl VEP's `minimised` mark, which `MINIMISED` prints under `--minimal`:
    /// set on a bi-allelic record `minimise_alleles` handled (Parser.pm:1011, so
    /// under `--minimal` on every bi-allelic sequence record, trimmed or not, and
    /// without the flag on a bi-allelic indel) and on every allele of a
    /// multi-allelic record whose `--minimal` split trimmed an allele
    /// (InputBuffer.pm:696).
    pub minimised: bool,
    /// Original allele string before minimisation.
    pub original_allele_string: Option<String>,
    /// The record's start before its `--minimal` split, for every allele of a
    /// multi-allelic record the split trimmed: `rejoin_variants_in_InputBuffer`
    /// (OutputFactory.pm:2291-2293) restores the record's own string and span
    /// before the output hash is built, so Location, the JSON `start` and `end`
    /// and the generated name read the record's span while the alleles and their
    /// consequences stay the trimmed ones. `None` for any other record.
    pub original_start: Option<u64>,
    /// The record's end before its `--minimal` split ([`Self::original_start`]).
    pub original_end: Option<u64>,
    /// The prefix Ensembl VEP's rejoin puts on this allele's transcript keys when
    /// it merges the alleles of a `--minimal` split back into the record: an
    /// allele after the first is stored under `<its allele string>_<transcript
    /// key>` (OutputFactory.pm:2316) and the rows print in the sorted key order
    /// (VariationFeature.pm:666), so `-/T_ENST...` sorts before the first
    /// allele's bare `ENST...` and `G/T_ENST...` after it. `None` for the first
    /// allele and for every record the split left alone.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rejoin_key_prefix: Option<Arc<str>>,
    /// The first allele's feature of this `--minimal` split record writes no row
    /// before the rejoin: every row it had was filtered on the feature's own span,
    /// or it has no alternate allele and no reference row. Ensembl VEP caches
    /// `ref_allele_string` on the first row it builds (OutputFactory.pm:1357,
    /// called at :2289 before the rejoin restores the record's string;
    /// VariationFeature.pm:1533-1539), so with no such row every row of the
    /// record prints the record's own reference as `REF_ALLELE`. Set by the row
    /// filters; `false` for every other record. Not part of the variant's JSON
    /// form.
    #[serde(skip)]
    pub split_first_feature_rowless: bool,
    /// The feature this allele, or the reference row leading it, stands alone as
    /// once `--check_frequency` dropped an allele of its batch: its position among
    /// the alternate alleles of the `--minimal` split record. Ensembl VEP's drop
    /// resets the buffer (BaseCacheVariation.pm:365-369, InputBuffer.pm:535-537),
    /// which cancels the rejoin of every split record in it (OutputFactory.pm:294),
    /// so each surviving feature is written as a record of its own, on its own
    /// span with its own string, and [`Self::same_record`] keeps two features of
    /// one line apart by it. `None` for every other allele. Not part of the
    /// variant's JSON form.
    #[serde(skip)]
    pub standalone_feature: Option<u32>,
    /// Raw input line (for passthrough in some output formats), one buffer shared by
    /// every allele of the record.
    pub raw_input: Option<Arc<str>>,
    /// The bytes of a captured VCF line that is not UTF-8, which the VCF output
    /// writes as read, where [`Self::raw_input`] holds the line with each invalid
    /// sequence replaced; `None` for a line that is UTF-8. One buffer shared by
    /// every allele of the record. Not part of the variant's JSON form, which
    /// carries the text of `raw_input`.
    #[serde(skip)]
    pub raw_input_bytes: Option<Arc<[u8]>>,
    /// The names of [`Self::colocated_variants`] that print on this allele's rows,
    /// in Ensembl VEP's order: non-somatic variants first, then by name prefix
    /// ([`colocated_in_output_order`]), keeping those matched on this allele or on
    /// their span alone (`OutputFactory.pm:1069-1086` at release/116.2).
    pub existing_variation: Vec<String>,
    /// The names of the known structural variants overlapping the variant under
    /// `--check_svs`, sorted and unique as Ensembl VEP prints them in `SV`
    /// (`OutputFactory.pm:962-964`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub overlapping_svs: Vec<String>,
    /// Plugin annotation data (variant-level).
    /// Keyed by field name (e.g., "CADD_PHRED"), value is the annotation string.
    /// Variant-level plugin data is included in the Extra column for every
    /// transcript consequence line.
    #[serde(default, skip_serializing_if = "indexmap::IndexMap::is_empty")]
    pub plugin_data: indexmap::IndexMap<String, String>,

    /// For structural variants: END position from VCF INFO field (1-based inclusive).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sv_end: Option<u64>,
    /// For structural variants: SVTYPE from VCF INFO field.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sv_type: Option<String>,
    /// For structural variants: SVLEN from VCF INFO field.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sv_len: Option<i64>,
    /// For imprecise SVs: CIPOS confidence interval (two i64 values).
    /// First value is typically negative (leftward expansion), second positive (rightward).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ci_pos: Option<(i64, i64)>,
    /// For imprecise SVs: CIEND confidence interval (two i64 values).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ci_end: Option<(i64, i64)>,
    /// Whether this is a structural variant (symbolic allele).
    #[serde(default)]
    pub is_structural: bool,
    /// For BND: the mate's id, INFO/MATEID or PARID when the record carries one, else
    /// a stand-in derived from the record. `Some` on every breakend of a record that
    /// has a mate and `None` on a record whose only ALT is a single breakend, which is
    /// how `display_allele()` and the breakend predicates tell the record's own
    /// breakend (`N.`, `.N`) from one whose Allele is `chromosome_breakpoint`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mate_id: Option<String>,
    /// For BND: mate chromosome parsed from the BND ALT notation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mate_chr: Option<String>,
    /// For BND: mate position parsed from the BND ALT notation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mate_pos: Option<u64>,
    /// For BND: whether this is the single-breakend allele form (e.g., "T.").
    #[serde(default)]
    pub is_single_breakend: bool,
    /// Ensembl VEP's `vep_skip` mark: a structural variant of a type its `%SO_TERMS`
    /// table lacks (`<CPX>`, `<NON_REF>`) or a deletion with neither `END` nor
    /// `SVLEN`, or whose start lies past its end (Parser/VCF.pm:477-481, 555-558).
    /// A marked record VEP keeps is carried in its VCF line without CSQ and as an
    /// `input`-only JSON object; the default and tab formats still list the rows.
    #[serde(default)]
    pub vep_skip: bool,
    /// A record whose span, end minus start, exceeds `--max_sv_size`. VEP carries
    /// its VCF line without CSQ and drops it from the JSON output altogether
    /// (`--json` sets the REST mode in which the parser skips it).
    #[serde(default)]
    pub oversize_sv: bool,
    /// The record of a line none of whose ALTs is an allele to annotate: every ALT a
    /// `*` (the VCF marker for a base an upstream deletion covers) or a `<DEL:*>`,
    /// equal to the reference or a repeat of an earlier ALT, or the `.` of a
    /// non-variant line. Ensembl VEP builds the feature over the reference bases with
    /// the allele string of every ALT (`N/*`) and no allele to annotate, writes no row
    /// in the default and tab formats, carries the VCF line without CSQ and writes the
    /// JSON object with `most_severe_consequence` `?`. No feature annotates it;
    /// `--nearest` still sets [`Self::nearest`].
    #[serde(default)]
    pub no_allele: bool,
    /// Optional annotation target chromosome override.
    ///
    /// Used for paired BND alleles, where the consequence should be calculated
    /// against the mate breakpoint while retaining the original local location
    /// in output.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_chr: Option<String>,
    /// Optional annotation target start position override.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_start: Option<u64>,
    /// Optional annotation target end position override.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_end: Option<u64>,
    /// Ordinal of the input record this allele came from (1-based line count of
    /// the input). Every allele split from one record shares it, which is how the
    /// writers reassemble a record's rows in VEP's order and emit one VCF line or
    /// one JSON object per record.
    #[serde(default)]
    pub input_record: u64,
    /// The allele string VEP prints in `Uploaded_variation` for an ID-less record:
    /// the raw REF/ALT of a bi-allelic record whose alleles differ in length (VEP
    /// minimises those and reports the untrimmed alleles), the anchor-trimmed
    /// alleles of every other record with all ALTs joined by `/`, or, for a
    /// structural variant, its ALTs joined by `/`, a breakpoint record's carrying the
    /// mate INFO/CHR2 names (`N[21:100[`). For a multi-allelic record whose
    /// `--minimal` split trimmed an allele, the record's chopped working string,
    /// which `split_variants` writes over the raw one (InputBuffer.pm:694) and
    /// `UPLOADED_ALLELE` then prints. One buffer shared by every allele of the
    /// record. `None` falls back to `allele_string`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uploaded_allele_string: Option<Arc<str>>,
    /// The record's allele string as VEP works with it, when it is not the allele's
    /// own: for a multi-allelic VCF record every allele with the shared first base
    /// chopped (`TT/T/-`), for a multi-allelic Ensembl-format line the whole allele
    /// column upper-cased (`TT/-/T`; Parser/VEP_input.pm:156, Parser.pm:554), for a
    /// breakpoint record the reference base joined to its ALT side
    /// (`N/N[21:100[/N]21:200]`) unless that side opens or closes with `.`, a native
    /// single breakend VEP keeps bare (Parser/VCF.pm:469 and 499). VEP's JSON
    /// `allele_string` prints this unless `--check_frequency` dropped an allele
    /// ([`Self::working_allele_string`]), while Uploaded_variation prints the raw
    /// one. One buffer shared by every allele of the record.
    #[serde(default)]
    pub record_allele_string_multi: Option<Arc<str>>,
    /// The record's working allele string in the case the input wrote it, when a
    /// base is lower case: the chopped and minimised alleles as they stand when
    /// Parser.pm:462 names a record whose ID Perl reads as false, before
    /// `validate_vf` upper-cases them (Parser.pm:554), and the alleles a record
    /// `--dont_skip` keeps after a refusal ahead of that line is annotated with.
    /// Read while it is the case variant of [`Self::record_allele_string`], so a
    /// copy whose string was rebuilt reads its own; `None` for an upper-case
    /// record. One buffer shared by every allele of the record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub written_case_allele_string: Option<Arc<str>>,
    /// The name VEP builds from the input line for a record its ID column cannot
    /// name (OutputFactory.pm `VariationFeature_to_output_hash`, the `_line`
    /// branches): `CHROM_POS_REF/ALT1/ALT2` from a VCF line's own columns, POS
    /// and alleles as written; `chr_start_alleles` from an Ensembl-format line,
    /// start the lower of its two coordinates. `None` for an input format VEP
    /// names from the variant itself. One buffer shared by every allele of the record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub input_line_name: Option<Arc<str>>,
    /// The SO term of the whole record's structural class, VEP's `class_SO_term`:
    /// the `%SO_TERMS` entry for the type abbreviation Parser.pm `get_SO_term`
    /// derives from every ALT joined by `/`, so a mobile-element record carries its
    /// subtype (`Alu_insertion`) and every ALT of a multi-allelic record carries the
    /// one term. The Allele column of every structural row prints it
    /// (OutputFactory.pm `BaseStructuralVariationOverlapAllele_to_output_hash`), and
    /// the predicates follow it rather than the allele's own form. For a sequence
    /// variant resolved from a variant identifier, the class the variation
    /// database stores for it (`VariationFeature.pm:1117`, read before any
    /// minimisation), and for an allele of a `--minimal` split of a record without
    /// one the class of the first allele's trimmed pair
    /// ([`crate::minimal::apply_minimal_split`]), which [`Self::class_so_term`]
    /// prints. `None` for any other sequence variant and for a structural type the
    /// table lacks.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub record_so_term: Option<String>,
    /// Bases a 3' shift along the genome moves this insertion or deletion when it
    /// overlaps no transcript, towards higher coordinates (`--shift_3prime 1`),
    /// zero when the flanking bases let it move nowhere; `None` when the run does
    /// not shift, the allele is not an indel, or it overlaps a transcript (its
    /// rows then carry each transcript's own shift).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub intergenic_shift: Option<u64>,
    /// The record's `--summary` or `--most_severe` row, carried by one allele of
    /// the record: the unique consequence terms of the rows left after the row
    /// filters and the pick family, in rank order (`--summary`), or the first of
    /// them alone (`--most_severe`); `OutputFactory.pm:426-458` at release/116.2.
    /// The allele then has no transcript consequences, its
    /// `most_severe_consequence` is the first term, and the other alleles of the
    /// record carry no row.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub summary_terms: Option<ConsequenceList>,
    /// The allele's intergenic row is the one a `--flag_pick`, `--flag_pick_allele`
    /// or `--flag_pick_allele_gene` run picked, so it prints `PICK=1`: Ensembl VEP
    /// flags the intergenic allele of a record that overlaps no transcript as it
    /// flags a transcript allele (`OutputFactory.pm:1363` at release/116.2).
    #[serde(default)]
    pub intergenic_pick: bool,
    /// The `FREQS` entries of the record under `--check_frequency`, one
    /// `POPULATION:ALLELE:frequency` per alternate allele of the record as
    /// uploaded, the alleles the filter dropped included, `NA` for an allele
    /// with no frequency (`BaseCacheVariation.pm:434-440`,
    /// `OutputFactory.pm:1148-1163` at release/116.2); empty when no
    /// co-located variant carried frequency data for the population. Every
    /// allele of the record carries the same list.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub freqs: Vec<String>,
    /// The record's allele string after `--check_frequency` dropped alleles:
    /// the reference and the kept alternates (`BaseCacheVariation.pm:355`),
    /// which the JSON `allele_string` and `AMBIGUITY` read in place of the
    /// uploaded one, and `UPLOADED_ALLELE` where Ensembl VEP prints the
    /// feature's string ([`Self::uploaded_allele`]), as Ensembl VEP rewrites
    /// the variation feature's own string. `None` when every allele was kept.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filtered_allele_string: Option<Arc<str>>,
    /// The record's reference allele did not match the reference sequence under
    /// `--check_ref`, and `--dont_skip` kept the record (Parser.pm:625 sets
    /// `check_ref_failed`); every row of the record prints `CHECK_REF` as `failed`
    /// (OutputFactory.pm:1020).
    #[serde(default)]
    pub check_ref_failed: bool,
    /// The sample this allele was copied for under `--individual`: Ensembl VEP
    /// copies the record's feature once per named individual
    /// (`Parser/VCF.pm:629-664`) and prints the name as `IND`
    /// (`OutputFactory.pm:971`). `None` for a record read without the flag or
    /// for a structural record, which is annotated before the per-individual
    /// step (`:244`). Every allele of one individual's copy carries the same
    /// name, and [`Self::same_record`] keeps the copies apart.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub individual: Option<Arc<str>>,
    /// The individual's genotype as Ensembl VEP stores it (`Parser/VCF.pm:657`):
    /// each called allele of the GT field, `.` left out, as the chopped allele
    /// string it stands for, the reference and a `*` included. `ZYG` counts the
    /// distinct values ([`Self::zygosity`]).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub genotype: Vec<String>,
    /// The position of the individual's copy in the `--individual` list, 0 for the
    /// first: Ensembl VEP copies the record once per listed name, a repeated name
    /// included (`Parser/VCF.pm:622`), and [`Self::same_record`] keeps two copies
    /// of one name apart by it.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub copy_ordinal: u32,
    /// The individual carries no alternate allele of the record
    /// (`Parser/VCF.pm:643`): `ZYG` ends in `REF`, and the copy is the reference
    /// alone, the record's reference row under `--process_ref_homs` or a feature
    /// with no allele to annotate otherwise. Also the record read under
    /// `--individual_zyg` for one named sample carrying no alternate allele, whose
    /// allele string Ensembl VEP sets to the reference twice (`:727-730`), so the
    /// row is the reference row and prints the reference's SPDI as such a copy does.
    #[serde(default)]
    pub hom_ref: bool,
    /// The named samples' zygosities of the record under `--individual_zyg`
    /// (`Parser/VCF.pm:683-743`): one entry per sample whose genotype is called,
    /// each name once, in the order the flag names them (the header's under
    /// `all`), which every format prints as the `ZYG` list (Ensembl VEP lists
    /// them in hash order, `OutputFactory.pm:983`). `None` for a record read
    /// without the flag or for a structural record, which is annotated before the
    /// per-sample step (`Parser/VCF.pm:244`); empty when no named sample's
    /// genotype is called, so the record prints no `ZYG`. One buffer shared by
    /// every allele of the record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_zygosities: Option<Arc<ZygosityList>>,
    /// The record's reference-allele row under `--process_ref_homs`: the allele is
    /// the reference itself, numbered 0 and listed before the alternates of every
    /// feature (`OutputFactory.pm:509-516`, `VariationFeatureOverlap.pm:469-474`,
    /// `BaseVariationFeatureOverlap.pm:363-373`), and no HGVS notation is built for
    /// it (`TranscriptVariationAllele.pm:1400`).
    #[serde(default)]
    pub reference_allele_row: bool,
    /// The alternate allele whose 3' shift places a reference row under
    /// `--shift_3prime 1`: the record's last alternate. `reset_shifted_positions`
    /// (OutputFactory.pm:306-323) leaves the transcript variation the record's
    /// alleles share at the coordinates `clear_shifting_variables`
    /// (TranscriptVariationAllele.pm:502-528) sets for the last allele with a
    /// shift, and the reference row prints those coordinates and reads its codon
    /// at them. `None` for every other allele. Not part of the variant's JSON
    /// form.
    #[serde(skip)]
    pub reference_row_shift_allele: Option<Arc<str>>,
}

/// An empty placeholder: no allocation, every field at its zero value. Lets a
/// worker take a finished variant out of a batch (`std::mem::take`) and drop it
/// where it was annotated, instead of leaving the drop to another thread.
impl Default for InputVariant {
    fn default() -> Self {
        Self {
            chr: String::new(),
            original_chr: String::new(),
            input_chr: None,
            start: 0,
            end: 0,
            strand: Strand::Forward,
            ref_allele: Vec::new(),
            alt_alleles: Vec::new(),
            allele_string: String::new(),
            allele_index: 0,
            allele_num: None,
            id: None,
            variant_class: VariantClass::Snv,
            transcript_consequences: Vec::new(),
            intergenic_overlap: false,
            most_severe_consequence: None,
            colocated_variants: Vec::new(),
            hgvsg: None,
            nearest: None,
            spdi: None,
            ga4gh_vrs: None,
            sibling_spdi: Vec::new(),
            minimised: false,
            original_allele_string: None,
            original_start: None,
            original_end: None,
            rejoin_key_prefix: None,
            split_first_feature_rowless: false,
            standalone_feature: None,
            raw_input: None,
            raw_input_bytes: None,
            existing_variation: Vec::new(),
            overlapping_svs: Vec::new(),
            plugin_data: indexmap::IndexMap::new(),
            sv_end: None,
            sv_type: None,
            sv_len: None,
            ci_pos: None,
            ci_end: None,
            is_structural: false,
            mate_id: None,
            mate_chr: None,
            mate_pos: None,
            is_single_breakend: false,
            vep_skip: false,
            oversize_sv: false,
            no_allele: false,
            annotation_chr: None,
            annotation_start: None,
            annotation_end: None,
            input_record: 0,
            uploaded_allele_string: None,
            record_allele_string_multi: None,
            written_case_allele_string: None,
            input_line_name: None,
            record_so_term: None,
            intergenic_shift: None,
            summary_terms: None,
            intergenic_pick: false,
            freqs: Vec::new(),
            filtered_allele_string: None,
            check_ref_failed: false,
            individual: None,
            genotype: Vec::new(),
            copy_ordinal: 0,
            hom_ref: false,
            sample_zygosities: None,
            reference_allele_row: false,
            reference_row_shift_allele: None,
        }
    }
}

impl InputVariant {
    /// Create a new input variant.
    ///
    /// The chromosome name is normalized at construction time (stripping `chr`
    /// prefix and converting `M` to `MT`) so that downstream matching against
    /// the Ensembl cache works regardless of the input naming convention.
    /// The original chromosome name is preserved for output (matching Ensembl VEP
    /// behavior of keeping the input naming convention in Location output).
    pub fn new(
        chr: String,
        start: u64,
        end: u64,
        ref_allele: Vec<u8>,
        alt_allele: Vec<u8>,
    ) -> Self {
        let allele_string = format!(
            "{}/{}",
            String::from_utf8_lossy(&ref_allele),
            String::from_utf8_lossy(&alt_allele),
        );
        let variant_class = classify_variant(&ref_allele, &alt_allele);
        let normalized_chr = normalize_chromosome(&chr);
        Self {
            original_chr: chr,
            input_chr: None,
            chr: normalized_chr,
            start,
            end,
            strand: Strand::Forward,
            ref_allele,
            alt_alleles: vec![alt_allele],
            allele_string,
            allele_index: 0,
            allele_num: None,
            id: None,
            variant_class,
            transcript_consequences: Vec::new(),
            intergenic_overlap: false,
            most_severe_consequence: None,
            colocated_variants: Vec::new(),
            hgvsg: None,
            nearest: None,
            spdi: None,
            ga4gh_vrs: None,
            sibling_spdi: Vec::new(),
            minimised: false,
            original_allele_string: None,
            original_start: None,
            original_end: None,
            rejoin_key_prefix: None,
            split_first_feature_rowless: false,
            standalone_feature: None,
            raw_input: None,
            raw_input_bytes: None,
            existing_variation: Vec::new(),
            overlapping_svs: Vec::new(),
            plugin_data: indexmap::IndexMap::new(),
            sv_end: None,
            sv_type: None,
            sv_len: None,
            ci_pos: None,
            ci_end: None,
            is_structural: false,
            mate_id: None,
            mate_chr: None,
            mate_pos: None,
            is_single_breakend: false,
            vep_skip: false,
            oversize_sv: false,
            no_allele: false,
            record_allele_string_multi: None,
            written_case_allele_string: None,
            annotation_chr: None,
            annotation_start: None,
            annotation_end: None,
            input_record: 0,
            uploaded_allele_string: None,
            input_line_name: None,
            record_so_term: None,
            intergenic_shift: None,
            summary_terms: None,
            intergenic_pick: false,
            freqs: Vec::new(),
            filtered_allele_string: None,
            check_ref_failed: false,
            individual: None,
            genotype: Vec::new(),
            copy_ordinal: 0,
            hom_ref: false,
            sample_zygosities: None,
            reference_allele_row: false,
            reference_row_shift_allele: None,
        }
    }

    /// A record with no allele to annotate, built over its reference alone: a
    /// non-variant VCF line (ALT `.`) or an identifier whose mapping names no
    /// alternate allele. The reference bases from `pos`, ending at
    /// `pos + length(REF) - 1` as Ensembl VEP's feature does (`Parser/VCF.pm:254`, so
    /// an empty REF ends before it starts), the allele string REF and every ALT of
    /// `alts` joined by `/`, and [`Self::no_allele`] set, so no feature annotates it
    /// and every writer carries it as Ensembl VEP carries a feature with no allele.
    pub fn without_allele<'a>(
        chr: String,
        pos: u64,
        ref_allele: &'a str,
        alts: impl IntoIterator<Item = &'a str>,
    ) -> Self {
        let end = (pos + ref_allele.len() as u64).saturating_sub(1);
        let mut variant = Self::new(chr, pos, end, ref_allele.as_bytes().to_vec(), b"*".to_vec());
        variant.allele_string = std::iter::once(ref_allele)
            .chain(alts)
            .collect::<Vec<_>>()
            .join("/");
        variant.no_allele = true;
        variant
    }

    /// Whether two alleles were split from one input record for the same
    /// individual: both carry the same record ordinal, or, when neither carries
    /// one (alleles built from a line alone), they share the line's identity; and
    /// the same [`Self::individual`] at the same [`Self::copy_ordinal`], since
    /// each individual's copy of a record is a feature of its own in Ensembl VEP,
    /// picked and written apart from the others; and the same
    /// [`Self::standalone_feature`], a split feature written on its own being a
    /// record of its own. The pipeline lays a batch out in record order, so the
    /// alleles of one record are adjacent and
    /// `chunk_by_mut(InputVariant::same_record)` yields one chunk per record, or
    /// per individual's copy of it.
    pub fn same_record(&self, other: &InputVariant) -> bool {
        if self.individual != other.individual
            || self.copy_ordinal != other.copy_ordinal
            || self.standalone_feature != other.standalone_feature
        {
            return false;
        }
        if self.input_record != 0 || other.input_record != 0 {
            return self.input_record == other.input_record;
        }
        self.original_chr == other.original_chr
            && self.start == other.start
            && self.id == other.id
            && self.ref_allele == other.ref_allele
            && self.uploaded_allele_string == other.uploaded_allele_string
    }

    /// The `ZYG` value of an individual's copy (`OutputFactory.pm:974-976`):
    /// the [`Zygosity`] of its genotype; `None` for an allele that is no
    /// individual's copy.
    pub fn zygosity(&self) -> Option<String> {
        self.individual.as_ref()?;
        Some(Zygosity::of(&self.genotype, self.hom_ref).to_string())
    }

    /// The record's reference-allele row (`--process_ref_homs`): the allele is
    /// the reference over the record's span, numbered 0, and the row carries the
    /// record-level fields of `self`, one allele of the record. Ensembl VEP builds
    /// it from the reference allele of the feature's string, so it takes the
    /// reference of the record's working string, which an individual's copy sets
    /// to its own (`record_allele_string_multi`), or, for an allele of a
    /// `--minimal` split, the allele's own trimmed reference over its trimmed
    /// span, the split feature's string and coordinates (InputBuffer.pm:674-676),
    /// with the record's span as [`Self::original_start`] and `original_end` and
    /// the allele's [`Self::rejoin_key_prefix`], so the row sorts with the
    /// allele's.
    pub fn reference_allele_row(&self) -> InputVariant {
        let record_allele_string = self.record_allele_string();
        let reference = if self.original_start.is_some() {
            self.ref_allele.clone()
        } else {
            record_allele_string
                .split('/')
                .next()
                .unwrap_or("-")
                .as_bytes()
                .to_vec()
        };
        let mut row = InputVariant::new(
            self.original_chr.clone(),
            self.start,
            self.end,
            reference.clone(),
            reference,
        );
        row.chr.clone_from(&self.chr);
        row.input_chr.clone_from(&self.input_chr);
        row.strand = self.strand;
        // The predicates read an allele of the reference's length as a
        // substitution (BaseVariationFeatureOverlapAllele.pm:377-378) and the
        // record's span as the feature's, so an insertion's reference row keeps
        // the insertion span's class for the span-based predicates.
        if self.variant_class == VariantClass::Insertion {
            row.variant_class = VariantClass::Insertion;
        }
        row.allele_num = Some(0);
        row.id.clone_from(&self.id);
        row.minimised = self.minimised;
        row.original_allele_string
            .clone_from(&self.original_allele_string);
        row.original_start = self.original_start;
        row.original_end = self.original_end;
        row.rejoin_key_prefix.clone_from(&self.rejoin_key_prefix);
        row.raw_input.clone_from(&self.raw_input);
        row.raw_input_bytes.clone_from(&self.raw_input_bytes);
        row.uploaded_allele_string
            .clone_from(&self.uploaded_allele_string);
        row.record_allele_string_multi = Some(record_allele_string.into());
        row.written_case_allele_string
            .clone_from(&self.written_case_allele_string);
        row.input_line_name.clone_from(&self.input_line_name);
        row.record_so_term.clone_from(&self.record_so_term);
        row.input_record = self.input_record;
        row.check_ref_failed = self.check_ref_failed;
        row.individual.clone_from(&self.individual);
        row.copy_ordinal = self.copy_ordinal;
        row.genotype.clone_from(&self.genotype);
        row.hom_ref = self.hom_ref;
        row.sample_zygosities.clone_from(&self.sample_zygosities);
        row.reference_allele_row = true;
        row.reference_row_shift_allele = (!self.no_allele)
            .then(|| Arc::from(String::from_utf8_lossy(self.alt_allele()).as_ref()));
        row
    }

    /// Location string in VEP format ("chr:start" or "chr:start-end"), under
    /// the chromosome name the record prints ([`Self::original_chr`]: the
    /// input's spelling, "chr21:100-200" for chr-prefixed input, or "MT" for a
    /// renamed bare "M").
    pub fn location(&self) -> String {
        let mut out = Vec::new();
        self.write_location(&mut out);
        written_text(out)
    }

    /// Appends [`Self::location`] to `out` without an intermediate `String`.
    pub fn write_location(&self, out: &mut Vec<u8>) {
        // Ensembl VEP's `rejoin_variants_in_InputBuffer` restores a multi-allelic
        // record's shared window coordinates for output.
        let s = self.original_start.unwrap_or(self.start);
        let e = self.original_end.unwrap_or(self.end);
        let lo = s.min(e);
        let hi = s.max(e);
        out.extend_from_slice(self.original_chr.as_bytes());
        out.push(b':');
        push_decimal(out, lo);
        if lo != hi {
            out.push(b'-');
            push_decimal(out, hi);
        }
    }

    /// The identifier that names the record in `Uploaded_variation`: [`Self::id`]
    /// when [`id_names_record`] admits it.
    fn explicit_id(&self) -> Option<&str> {
        self.id.as_deref().filter(|id| id_names_record(id))
    }

    /// The chromosome spelling the record's generated names are built from: the
    /// input's, which Ensembl VEP reads before it renames `M` (Parser.pm:462).
    fn name_chr(&self) -> &str {
        self.input_chr.as_deref().unwrap_or(&self.original_chr)
    }

    /// `Uploaded_variation` as VEP prints it (OutputFactory.pm
    /// `VariationFeature_to_output_hash`): the identifier `explicit_id` admits unless
    /// it is the `chr_start_alleles` name VEP gives an ID-less record; else, for a
    /// record whose ID Perl reads as false and whose string changed after
    /// Parser.pm:462 named it, renamed from `M` or rewritten by `--check_frequency`,
    /// that name ([`Self::variation_name`]); else the name built from the input
    /// line ([`Self::input_line_name`]); else, for an input format VEP names from
    /// the variant, `chr_start_alleles` where start is the anchor-trimmed start and
    /// the alleles are [`Self::uploaded_allele_string`].
    pub fn uploaded_variation(&self) -> String {
        let mut out = Vec::new();
        self.write_uploaded_variation(&mut out);
        written_text(out)
    }

    /// Appends [`Self::uploaded_variation`] to `out` without an intermediate
    /// `String`.
    pub fn write_uploaded_variation(&self, out: &mut Vec<u8>) {
        if let Some(id) = self.explicit_id().filter(|id| !self.is_generated_name(id)) {
            out.extend_from_slice(id.as_bytes());
            return;
        }
        // A record whose ID column Perl reads as false is named before the
        // `M` rename and the `--check_frequency` rewrite (Parser.pm:462) and
        // prints that name when it differs from the one built at output
        // (OutputFactory.pm:887-900), as either makes it; a `.` ID falls through
        // to the line's name there as here.
        if (self.input_chr.is_some() || self.filtered_allele_string.is_some())
            && !self.id.as_deref().is_some_and(perl_true)
        {
            out.extend_from_slice(self.variation_name().as_bytes());
            return;
        }
        if let Some(name) = &self.input_line_name {
            out.extend_from_slice(name.as_bytes());
            return;
        }
        let start = self.original_start.unwrap_or(self.start);
        let alleles = self
            .uploaded_allele_string
            .as_deref()
            .or(self.original_allele_string.as_deref())
            .unwrap_or(&self.allele_string);
        out.extend_from_slice(self.name_chr().as_bytes());
        out.push(b'_');
        push_decimal(out, start);
        out.push(b'_');
        out.extend_from_slice(alleles.as_bytes());
    }

    /// Whether `id` is the `chr_start_alleles` name VEP gives a record whose ID
    /// column is empty (Parser.pm `validate_vf`). The name is built as
    /// `VariationFeature_to_output_hash` builds the one it compares the ID with
    /// (OutputFactory.pm:887-891), from the record's chromosome and allele string
    /// as they stand at output, `MT` for a renamed `M`, so [`Self::original_chr`]
    /// rather than [`Self::name_chr`], and [`Self::working_allele_string`]. The
    /// prefix test keeps an identifier such as `rs123` from building the name on
    /// every row.
    fn is_generated_name(&self, id: &str) -> bool {
        let chr = self.original_chr.as_str();
        if !id.starts_with(chr) {
            return false;
        }
        // The start as the output reads it: `rejoin_variants_in_InputBuffer`
        // (OutputFactory.pm:294) restores a minimised record's own start before the
        // output hash is built.
        let start = self.original_start.unwrap_or(self.start);
        id == format!("{}_{}_{}", chr, start, self.working_allele_string())
    }

    /// The allele string of the whole input record as VEP holds it: the symbolic
    /// ALT(s) of a structural variant, [`Self::record_allele_string_multi`] for a
    /// multi-allelic or breakpoint record, or the minimised `allele_string` of a
    /// bi-allelic record.
    pub fn record_allele_string(&self) -> String {
        if let Some(multi) = &self.record_allele_string_multi {
            return multi.to_string();
        }
        match self.uploaded_allele_string.as_deref() {
            Some(s) if self.is_structural => s.to_string(),
            _ => self.allele_string.clone(),
        }
    }

    /// The record's allele string as the run works with it after the input was
    /// read: [`Self::filtered_allele_string`] once `--check_frequency` dropped
    /// alleles, else [`Self::record_allele_string`]. Ensembl VEP's variation
    /// feature holds this string; its `variation_name`, built before the
    /// filter runs, keeps the uploaded one.
    pub fn working_allele_string(&self) -> String {
        match &self.filtered_allele_string {
            Some(filtered) => filtered.to_string(),
            None => self.record_allele_string(),
        }
    }

    /// VEP's `variation_name` after `validate_vf`, the `id` of its JSON output
    /// (OutputFactory/JSON.pm:236): [`Self::id`] when Perl reads it as true, so a
    /// `.` stays `.`; else `chr_start_alleles` built from
    /// [`Self::record_allele_string`] in the case the input wrote it
    /// ([`Self::written_case_allele_string`]), the name Parser.pm:462 `||=` gives
    /// a record whose ID is absent, empty or `0`.
    pub fn variation_name(&self) -> String {
        if let Some(id) = self.id.as_deref().filter(|id| perl_true(id)) {
            return id.to_string();
        }
        let start = self.original_start.unwrap_or(self.start);
        let record_alleles = self.record_allele_string();
        let alleles = match self.written_case_allele_string.as_deref() {
            Some(written) if written.eq_ignore_ascii_case(&record_alleles) => written.to_string(),
            _ => record_alleles,
        };
        format!("{}_{}_{}", self.name_chr(), start, alleles)
    }

    /// The specific alternate allele for this variant record.
    pub fn alt_allele(&self) -> &[u8] {
        self.alt_alleles
            .first()
            .map(|a| a.as_slice())
            .unwrap_or(b"-")
    }

    /// Allele string for the specific alt allele in display format.
    ///
    /// For a structural variant, the record's SO term ([`Self::record_so_term`], the
    /// `class_SO_term` every row of the record prints in Ensembl VEP), or for a
    /// variant built without one the term of its class, read from its own ALT for
    /// a mobile element. A bracket breakend prints its string as written
    /// (`A]21:33033339]`), the own breakend of a record that names a mate its `A.`
    /// or `.A`, and a translocation allele without a breakend string, a single
    /// breakend on its own or a bare `<BND>`, the class term `chromosome_breakpoint`:
    /// Ensembl VEP's `$vfoa->{breakend}->{string} || $svf->class_SO_term`
    /// (`OutputFactory.pm:2029`). A type the term table lacks prints the type
    /// string Ensembl VEP keeps as the `class_SO_term` ([`Self::structural_type_string`]:
    /// `<CPX>` for a bare `<CPX>`, `CPX` under `SVTYPE=CPX`).
    /// For `<NON_REF>` alleles, returns the raw allele string since Ensembl VEP outputs
    /// the literal `<NON_REF>` (it has no SO term mapping).
    /// For small variants, returns the actual allele bases.
    pub fn display_allele(&self) -> String {
        let mut out = Vec::new();
        self.write_display_allele(&mut out);
        written_text(out)
    }

    /// Ensembl VEP's type string for a structural record (Parser/VCF.pm:469-475): the
    /// ALTs joined by `/`, or INFO/SVTYPE when that is true and the joined string
    /// opens outside the VCF 4.4 forms ([`sv_type_string`]). A type `%SO_TERMS`
    /// lacks keeps it as the record's `class_SO_term` (`:477-481`), so `<CPX>` prints
    /// `<CPX>`, `<CPX>` with `SVTYPE=CPX` prints `CPX` and with `SVTYPE=FOO` `FOO`.
    fn structural_type_string(&self) -> String {
        let joined = match &self.uploaded_allele_string {
            Some(alts) => alts.to_string(),
            None => String::from_utf8_lossy(self.alt_allele()).into_owned(),
        };
        sv_type_string(&joined, self.sv_type.as_deref())
    }

    /// Appends [`Self::display_allele`] to `out` without an intermediate
    /// `String`.
    pub fn write_display_allele(&self, out: &mut Vec<u8>) {
        let raw = String::from_utf8_lossy(self.alt_allele());
        let type_string = (self.variant_class == VariantClass::ComplexStructural)
            .then(|| self.structural_type_string());
        let allele: &str = if self.is_structural || self.variant_class.is_structural() {
            if self.variant_class == VariantClass::Translocation {
                if self.mate_id.is_some() {
                    &raw
                } else {
                    self.variant_class.so_term()
                }
            } else if let Some(term) = &self.record_so_term {
                term
            } else if self.variant_class == VariantClass::ReferenceBlock {
                // Ensembl VEP outputs the allele as written (no SO term mapping).
                &raw
            } else if let Some(type_string) = &type_string {
                type_string
            } else if self.variant_class == VariantClass::MobileElementInsertion
                || self.variant_class == VariantClass::MobileElementDeletion
            {
                // The ALT's own type read as the record's would be (`<INS:ME:ALU>` is
                // `Alu_insertion`, `<DEL:ME>` `mobile_element_deletion`).
                sv_so_term_and_class(&sv_type_abbreviation(&raw, None))
                    .map_or(self.variant_class.so_term(), |(term, _)| term)
            } else {
                self.variant_class.so_term()
            }
        } else {
            // The allele as the row's overlap allele holds it, `variation_feature_seq`
            // (OutputFactory.pm:1347): for an allele of a `--minimal` split that is
            // the trimmed allele, while Location prints the record's span.
            &raw
        };
        out.extend_from_slice(allele.as_bytes());
    }

    /// The `VARIANT_CLASS` value of the record, Ensembl VEP's `class_SO_term`: a
    /// structural record's own term ([`Self::record_so_term`], `Alu_insertion` for
    /// `<INS:ME:ALU>`) or, for a type the term table lacks, the string its `Allele`
    /// column prints ([`Self::display_allele`], `<CPX>` for a bare `<CPX>` and `CPX`
    /// under `SVTYPE=CPX`), since Ensembl VEP prints the one `class_SO_term` in both
    /// columns of a structural row; a
    /// breakend row, whose `Allele` is its breakend string, prints the class term
    /// `chromosome_breakpoint`. A sequence record prints the class its record
    /// carries ([`Self::record_so_term`], the one a variation database stores for
    /// a resolved identifier) or else [`so_variation_class`] of
    /// [`Self::working_allele_string`] upper-cased, the allele string the whole
    /// record works with after `validate_vf` (`Parser.pm:554`) and after
    /// `--check_frequency` dropped alleles (BaseCacheVariation.pm:333, :355 rewrite
    /// the feature's string before any source reads its class, the filtering
    /// source being sorted first at AnnotationSourceAdaptor.pm:98, and
    /// VariationFeature.pm:1118 classes that string), so every allele of a
    /// multi-allelic record prints the one class (`sequence_alteration` for
    /// `TG/T/-`, `indel` for `T/A/TT`, `SNV` for `A/C/G`), a record whose `-`
    /// allele the filter dropped prints the kept pair's class, and a record written
    /// in lower case prints the class of its upper-cased alleles.
    pub fn class_so_term(&self) -> String {
        if self.is_structural || self.variant_class.is_structural() {
            return match &self.record_so_term {
                Some(term) => term.clone(),
                None if self.variant_class == VariantClass::Translocation => {
                    self.variant_class.so_term().to_string()
                }
                None => self.display_allele(),
            };
        }
        match &self.record_so_term {
            Some(term) => term.clone(),
            None => {
                so_variation_class(&self.working_allele_string().to_ascii_uppercase()).to_string()
            }
        }
    }

    /// The `AMBIGUITY` value of the record: [`ambiguity_code`] of
    /// [`Self::working_allele_string`], the allele string Ensembl VEP's
    /// `VariationFeature::ambig_code` reads, and `None` for a structural variant,
    /// whose feature class carries no code (OutputFactory.pm
    /// `VariationFeature_to_output_hash`).
    pub fn ambiguity_code(&self) -> Option<u8> {
        if self.is_structural {
            return None;
        }
        ambiguity_code(&self.working_allele_string())
    }

    /// The `UPLOADED_ALLELE` value of the record: the allele string as uploaded
    /// ([`Self::uploaded_allele_string`], else the allele's own) where Ensembl VEP
    /// keeps it as `original_allele_string` ([`keeps_uploaded_spelling`] under the
    /// run's `--minimal` setting), else the feature's working string, upper-cased
    /// ([`Self::working_allele_string`]: the string `--check_frequency` rewrote or
    /// `--lookup_ref` gave its reference), as `OutputFactory.pm:1360` prints
    /// `original_allele_string || allele_string`; `None` for a structural variant,
    /// whose rows carry neither it nor `REF_ALLELE`
    /// (`BaseStructuralVariationOverlapAllele_to_output_hash` sets `Allele`,
    /// `ALLELE_NUM`, `Consequence`, `IMPACT` and `PICK` alone).
    pub fn uploaded_allele(&self, minimal: bool) -> Option<String> {
        if self.is_structural {
            return None;
        }
        let raw = self
            .uploaded_allele_string
            .as_deref()
            .unwrap_or(&self.allele_string);
        if keeps_uploaded_spelling(raw, minimal) {
            Some(raw.to_string())
        } else {
            Some(self.working_allele_string().to_ascii_uppercase())
        }
    }

    /// [`Self::location`] after a 3' shift of `shift` bases along `strand`
    /// (`--shift_genomic 1`): a transcript row moves by `shift_length * strand`
    /// (OutputFactory.pm:1705), an intergenic row towards higher coordinates
    /// (:1997).
    pub fn location_shifted(&self, shift: u64, strand: i8) -> String {
        let moved = |pos: u64| {
            if strand < 0 {
                pos.saturating_sub(shift)
            } else {
                pos + shift
            }
        };
        let s = moved(self.original_start.unwrap_or(self.start));
        let e = moved(self.original_end.unwrap_or(self.end));
        let (lo, hi) = (s.min(e), s.max(e));
        if lo == hi {
            format!("{}:{}", self.original_chr, lo)
        } else {
            format!("{}:{}-{}", self.original_chr, lo, hi)
        }
    }

    /// The variant read against `reference` in place of its reference allele,
    /// as `--use_transcript_ref` reads each row against the transcript's
    /// (`_raw_allele_hashes`, VariationFeatureOverlap.pm:418-425: the input's
    /// reference allele leaves the allele list and the feature's takes the
    /// reference allele's object). A reference row is that allele, so its
    /// alternate allele follows.
    pub fn against_reference_allele(&self, reference: &[u8]) -> InputVariant {
        let mut reread = self.clone();
        reread.ref_allele = reference.to_vec();
        if reread.reference_allele_row {
            reread.alt_alleles = vec![reference.to_vec()];
        }
        reread.allele_string = std::iter::once(&reread.ref_allele)
            .chain(&reread.alt_alleles)
            .map(|a| String::from_utf8_lossy(a).into_owned())
            .collect::<Vec<_>>()
            .join("/");
        reread
    }

    /// The variant with its alleles on the forward strand: an Ensembl-format
    /// record gives them on the strand its column names, and Ensembl VEP reads
    /// each against a feature as `feature_seq`, reverse-complemented when the
    /// two strands differ (VariationFeatureOverlapAllele.pm:253-256), so a copy
    /// on the forward strand analyses the same way against every transcript. A
    /// forward-strand variant is returned as it is, as is a structural variant,
    /// whose alleles are class words that no `feature_seq` reads. Only an allele
    /// `seq_is_dna` admits (`:315-321`: the nucleotide and IUPAC letters and `-`)
    /// is reverse-complemented, with Ensembl's `reverse_comp`
    /// ([`crate::codon::reverse_complement`]); any other allele is kept as given,
    /// and a `-` allele stays `-`.
    pub fn on_forward_strand(&self) -> InputVariant {
        if !self.alleles_are_read_through_the_complement() {
            return self.clone();
        }
        let mut forward = self.clone();
        forward.strand = Strand::Forward;
        forward.ref_allele = complement_dna_allele(&self.ref_allele);
        forward.alt_alleles = self
            .alt_alleles
            .iter()
            .map(|a| complement_dna_allele(a))
            .collect();
        forward.allele_string = std::iter::once(&forward.ref_allele)
            .chain(&forward.alt_alleles)
            .map(|a| String::from_utf8_lossy(a).into_owned())
            .collect::<Vec<_>>()
            .join("/");
        forward
    }

    /// The reference allele and this row's alternate allele on the forward
    /// strand, the orientation a data file keyed on the chromosome's strand
    /// compares them in: a plugin reverse-complements the alleles of a record
    /// given on the reverse strand before it looks a file row up
    /// (`get_matched_variant_alleles`, Utils/Sequence.pm:1182 and 1196, and the
    /// plugins' own `reverse_comp`). Borrowed as given unless
    /// [`Self::on_forward_strand`] would complement them.
    pub fn forward_strand_alleles(&self) -> (Cow<'_, [u8]>, Cow<'_, [u8]>) {
        if self.alleles_are_read_through_the_complement() {
            (
                Cow::Owned(complement_dna_allele(&self.ref_allele)),
                Cow::Owned(complement_dna_allele(self.alt_allele())),
            )
        } else {
            (
                Cow::Borrowed(&self.ref_allele),
                Cow::Borrowed(self.alt_allele()),
            )
        }
    }

    /// Whether [`Self::on_forward_strand`] complements the alleles: given on the
    /// reverse strand and not a structural variant, whose alleles are class words
    /// that no `feature_seq` reads.
    fn alleles_are_read_through_the_complement(&self) -> bool {
        self.strand == Strand::Reverse && !self.is_structural && !self.variant_class.is_structural()
    }

    /// [`Self::display_allele`] after a 3' shift of `shift` bases along `strand`
    /// (`--shift_genomic 1`, OutputFactory.pm:1349-1352): the inserted bases as
    /// given, rotated along the move, the `hgvs_allele_string` `perform_shift`
    /// keeps (TranscriptVariationAllele.pm:300-301, 332-341: the rotation runs the
    /// other way for an allele given on the reverse strand); every other allele is
    /// unchanged.
    pub fn display_allele_shifted(&self, shift: u64, strand: i8) -> String {
        let inserted = self.alt_allele();
        if self.variant_class != VariantClass::Insertion || inserted.is_empty() || inserted == b"-"
        {
            return self.display_allele();
        }
        let mut rotated = inserted.to_vec();
        let steps = (shift % rotated.len() as u64) as usize;
        let moves_towards_higher = strand >= 0;
        let given_on_forward = self.strand == Strand::Forward;
        if moves_towards_higher == given_on_forward {
            rotated.rotate_left(steps);
        } else {
            rotated.rotate_right(steps);
        }
        String::from_utf8_lossy(&rotated).to_string()
    }
}

/// The `skip_serializing_if` predicate of [`InputVariant::copy_ordinal`].
fn is_zero(value: &u32) -> bool {
    *value == 0
}

/// An allele read on the other strand: Ensembl's `reverse_comp`
/// ([`crate::codon::reverse_complement`]) of an allele `seq_is_dna` admits
/// (VariationFeatureOverlapAllele.pm:315-321: the nucleotide and IUPAC letters
/// and `-`, so a `-` allele stays `-`); any other allele as given.
fn complement_dna_allele(allele: &[u8]) -> Vec<u8> {
    if crate::codon::seq_is_dna(allele) {
        crate::codon::reverse_complement(allele)
    } else {
        allele.to_vec()
    }
}

/// Appends `n` in decimal, byte for byte as `u64`'s `Display` writes it, without
/// the formatting machinery on a path that runs once per row; the buffer holds
/// the 20 digits of `u64::MAX`.
pub fn push_decimal(out: &mut Vec<u8>, mut n: u64) {
    let mut digits = [0u8; 20];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    out.extend_from_slice(&digits[start..]);
}

/// The text a byte writer appended: whole `str` values and ASCII, so UTF-8.
fn written_text(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).expect("byte writers append str slices and ASCII")
}

/// The Sequence Ontology class of an allele string, as ensembl-variation
/// `Utils/Sequence.pm` `SO_variation_class` derives it with a trusted reference
/// allele: single upper-case letters throughout are an `SNV`; `)` followed by a
/// digit is a `tandem_repeat`; a `-` reference with sequence or `INS` alleles an
/// `insertion` and with one `DEL` allele a `deletion`; a sequence reference with
/// `-` or `DEL` alleles alone a `deletion`, and with sequence alleles a
/// `substitution` when every allele has its length, else an `indel`; a `DEL`
/// reference with `-` alleles a `deletion`; a lone `INS` or `DEL` allele an
/// `insertion` or `deletion`; anything else, a `*` allele or deletions of two
/// lengths included, a `sequence_alteration`.
pub fn so_variation_class(allele_string: &str) -> &'static str {
    let is_sequence = |a: &str| !a.is_empty() && a.bytes().all(|b| b.is_ascii_uppercase());
    let parts: Vec<&str> = allele_string.split('/').collect();
    if parts.len() > 1 && parts.iter().all(|a| a.len() == 1 && is_sequence(a)) {
        return "SNV";
    }
    let tandem_repeat = allele_string
        .as_bytes()
        .windows(2)
        .any(|w| w[0] == b')' && w[1].is_ascii_digit());
    if tandem_repeat {
        return "tandem_repeat";
    }
    let mut class = "sequence_alteration";
    match parts.as_slice() {
        [reference, alts @ ..] if !alts.is_empty() => {
            if *reference == "-" {
                if alts.iter().all(|a| is_sequence(a) || a.contains("INS")) {
                    class = "insertion";
                } else if alts.len() == 1 && alts[0].contains("DEL") {
                    class = "deletion";
                }
            } else if is_sequence(reference) {
                if alts.iter().all(|a| a.contains('-') || a.contains("DEL")) {
                    class = "deletion";
                } else if alts.iter().all(|a| is_sequence(a)) {
                    class = if alts.iter().all(|a| a.len() == reference.len()) {
                        "substitution"
                    } else {
                        "indel"
                    };
                }
            } else if reference.contains("DEL") && alts.iter().all(|a| a.contains('-')) {
                class = "deletion";
            }
        }
        [single] => {
            if single.contains("INS") {
                class = "insertion";
            } else if single.contains("DEL") {
                class = "deletion";
            }
        }
        _ => {}
    }
    class
}

/// The IUPAC ambiguity code of an allele string, as ensembl-variation
/// `Utils/Sequence.pm` `ambiguity_code` derives it: the alleles are split on `/`,
/// `|` or `\`, an `N` stands for the four bases, and the distinct alleles,
/// sorted as strings, are joined and upper-cased into the key of a sixteen-entry
/// table (`AC` is `M`, `ACGT` is `N`, `-` alone keeps `-`). The key is built
/// from whole alleles, so `G/AC` reads as `ACG` and prints `V` while `A/AT`
/// reads as `AAT` and prints nothing; every anchor-trimmed indel (`-/T`, `TA/-`)
/// and every record carrying a `*` has no code.
pub fn ambiguity_code(allele_string: &str) -> Option<u8> {
    let mut alleles: Vec<&str> = Vec::new();
    for allele in allele_string.split(['/', '|', '\\']) {
        if allele == "N" {
            alleles.extend(["A", "C", "G", "T"]);
        } else {
            alleles.push(allele);
        }
    }
    alleles.sort_unstable();
    alleles.dedup();
    match alleles.concat().to_ascii_uppercase().as_str() {
        "AC" => Some(b'M'),
        "ACG" => Some(b'V'),
        "ACGT" => Some(b'N'),
        "ACT" => Some(b'H'),
        "AG" => Some(b'R'),
        "AGT" => Some(b'D'),
        "AT" => Some(b'W'),
        "CG" => Some(b'S'),
        "CGT" => Some(b'B'),
        "CT" => Some(b'Y'),
        "GT" => Some(b'K'),
        "C" => Some(b'C'),
        "A" => Some(b'A'),
        "T" => Some(b'T'),
        "G" => Some(b'G'),
        "-" => Some(b'-'),
        _ => None,
    }
}

/// The `ZYG` value of one genotype (`OutputFactory.pm:974-976` and `:984-985` at
/// release/116.2): `HET` when the called alleles hold more than one distinct
/// value, else `HOM` (a genotype with no called allele included), with `REF`
/// appended when the individual carries no alternate allele of the record.
/// Serialised as the printed label, so a serialised [`InputVariant`] names the
/// zygosity as every output format does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum Zygosity {
    Het,
    Hom,
    HetRef,
    HomRef,
}

impl Zygosity {
    /// The value for `genotype`, the called alleles as strings, and whether the
    /// individual carries no alternate allele.
    pub fn of<A: AsRef<str>>(genotype: &[A], hom_ref: bool) -> Self {
        let heterozygous = match genotype.split_first() {
            Some((first, rest)) => rest.iter().any(|allele| allele.as_ref() != first.as_ref()),
            None => false,
        };
        match (heterozygous, hom_ref) {
            (true, false) => Zygosity::Het,
            (false, false) => Zygosity::Hom,
            (true, true) => Zygosity::HetRef,
            (false, true) => Zygosity::HomRef,
        }
    }

    /// The value as the output prints it.
    pub fn label(self) -> &'static str {
        match self {
            Zygosity::Het => "HET",
            Zygosity::Hom => "HOM",
            Zygosity::HetRef => "HETREF",
            Zygosity::HomRef => "HOMREF",
        }
    }
}

impl std::fmt::Display for Zygosity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.label())
    }
}

/// One named sample's zygosity of a record under `--individual_zyg`
/// (`Parser/VCF.pm:723-737` stores the genotype and the reference-homozygous mark
/// per sample; `OutputFactory.pm:981-987` prints them), one entry of the `ZYG`
/// list, written `sample:ZYG`.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SampleZygosity {
    /// The sample's name, as the `#CHROM` header spells it.
    pub sample: Arc<str>,
    pub zygosity: Zygosity,
}

impl std::fmt::Display for SampleZygosity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.sample, self.zygosity)
    }
}

/// The `ZYG` list of a record under `--individual_zyg`, shared by every row of the
/// record: its entries, and the list as each text format joins it, rendered once
/// per separator since a many-sample record writes the same list on every row.
#[derive(Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct ZygosityList {
    pub entries: Vec<SampleZygosity>,
    #[serde(skip)]
    comma_joined: OnceLock<String>,
    #[serde(skip)]
    ampersand_joined: OnceLock<String>,
}

impl ZygosityList {
    pub fn new(entries: Vec<SampleZygosity>) -> Self {
        Self {
            entries,
            comma_joined: OnceLock::new(),
            ampersand_joined: OnceLock::new(),
        }
    }

    /// The entries joined by `separator`: `,` as the default and tab formats
    /// join a list, `&` as the VCF's CSQ does, each rendered once.
    pub fn joined(&self, separator: &str) -> Cow<'_, str> {
        let render = || {
            let mut out = String::new();
            for (index, entry) in self.entries.iter().enumerate() {
                if index > 0 {
                    out.push_str(separator);
                }
                out.push_str(&entry.sample);
                out.push(':');
                out.push_str(entry.zygosity.label());
            }
            out
        };
        match separator {
            "," => Cow::Borrowed(self.comma_joined.get_or_init(render)),
            "&" => Cow::Borrowed(self.ampersand_joined.get_or_init(render)),
            _ => Cow::Owned(render()),
        }
    }
}

impl PartialEq for ZygosityList {
    fn eq(&self, other: &Self) -> bool {
        self.entries == other.entries
    }
}

impl Eq for ZygosityList {}

impl Clone for ZygosityList {
    fn clone(&self) -> Self {
        Self::new(self.entries.clone())
    }
}

/// Whether `UPLOADED_ALLELE` prints the record's allele string as uploaded
/// (Ensembl VEP's `original_allele_string` is set) rather than the feature's
/// working string: the run is `--minimal` (`minimise_alleles` sets it for every
/// record), or the alleles differ in length with no `-` among them (`Parser.pm`
/// `post_process_vfs` ends its length test at the first `-`); `validate_vf`
/// upper-cases the working string every other record prints
/// ([`InputVariant::uploaded_allele`]).
pub fn keeps_uploaded_spelling(allele_string: &str, minimal: bool) -> bool {
    let mut parts = allele_string.split('/');
    let reference_length = parts.next().map_or(0, str::len);
    minimal
        || (!allele_string.contains('-') && parts.any(|allele| allele.len() != reference_length))
}

#[cfg(test)]
mod strand_tests {
    use super::*;

    /// An Ensembl-format allele on the reverse strand reads on the forward strand
    /// as its reverse complement, a `-` staying `-`; a forward-strand variant is
    /// its own copy.
    #[test]
    fn on_forward_strand_reverse_complements_the_alleles_of_a_reverse_strand_input() {
        let mut given = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"TGA".to_vec());
        given.strand = Strand::Reverse;
        let forward = given.on_forward_strand();
        assert_eq!(forward.strand, Strand::Forward);
        assert_eq!(forward.ref_allele, b"-");
        assert_eq!(forward.alt_allele(), b"TCA");
        assert_eq!(forward.allele_string, "-/TCA");
        assert_eq!(forward.variant_class, given.variant_class);
        let snv = InputVariant::new("21".into(), 100, 100, b"G".to_vec(), b"T".to_vec());
        assert_eq!(snv.on_forward_strand().allele_string, "G/T");
    }

    /// `feature_seq` flips only an allele `seq_is_dna` admits
    /// (VariationFeatureOverlapAllele.pm:253, 315-321), with Ensembl's
    /// `reverse_comp`: an IUPAC letter takes its complement (`R` to `Y`), a
    /// structural-variant word or a lower-case run of bases is kept or
    /// complemented as written, and a `*` or symbolic allele is left as given.
    /// A structural variant is not read through `feature_seq` at all, so its
    /// copy is itself, a word the letter set would admit (`CNV`) included.
    #[test]
    fn on_forward_strand_flips_dna_alleles_only_with_the_iupac_complement() {
        let flip = |reference: &[u8], alt: &[u8]| {
            let mut given =
                InputVariant::new("21".into(), 100, 100, reference.to_vec(), alt.to_vec());
            given.strand = Strand::Reverse;
            given.on_forward_strand().allele_string
        };
        assert_eq!(flip(b"G", b"R"), "C/Y");
        assert_eq!(flip(b"ag", b"tn"), "ct/na");
        assert_eq!(flip(b"DEL", b"DUP"), "DEL/DUP");
        assert_eq!(flip(b"G", b"<DEL>"), "C/<DEL>");
        assert_eq!(flip(b"G", b"*"), "C/*");
        let mut sv = InputVariant::new("21".into(), 100, 900, b"CNV".to_vec(), b"CNV".to_vec());
        sv.strand = Strand::Reverse;
        sv.is_structural = true;
        let copy = sv.on_forward_strand();
        assert_eq!(
            (copy.strand, copy.allele_string.as_str()),
            (Strand::Reverse, "CNV/CNV")
        );
    }

    /// The alleles a data file is matched against are the record's on the forward
    /// strand: `T/C` given on the reverse strand reads `A/G`, an insertion's `-`
    /// stays `-`, a forward-strand record and a structural variant give their
    /// alleles as written.
    #[test]
    fn forward_strand_alleles_complement_a_reverse_strand_record_only() {
        let alleles = |variant: &InputVariant| {
            let (reference, alt) = variant.forward_strand_alleles();
            (reference.into_owned(), alt.into_owned())
        };
        let mut reverse = InputVariant::new("21".into(), 100, 100, b"T".to_vec(), b"C".to_vec());
        reverse.strand = Strand::Reverse;
        assert_eq!(alleles(&reverse), (b"A".to_vec(), b"G".to_vec()));
        let mut insertion =
            InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"TGA".to_vec());
        insertion.strand = Strand::Reverse;
        assert_eq!(alleles(&insertion), (b"-".to_vec(), b"TCA".to_vec()));
        let forward = InputVariant::new("21".into(), 100, 100, b"T".to_vec(), b"C".to_vec());
        assert_eq!(alleles(&forward), (b"T".to_vec(), b"C".to_vec()));
        let mut sv = InputVariant::new("21".into(), 100, 900, b"CNV".to_vec(), b"CNV".to_vec());
        sv.strand = Strand::Reverse;
        sv.is_structural = true;
        assert_eq!(alleles(&sv), (b"CNV".to_vec(), b"CNV".to_vec()));
    }

    /// `perform_shift` rotates the allele as given the other way when it was
    /// given on the reverse strand (TranscriptVariationAllele.pm:300-301,
    /// 332-341): `-/TCA` moved 10 bases along a forward transcript prints `CAT`,
    /// its reverse-strand spelling `-/TGA` prints `ATG`.
    #[test]
    fn display_allele_shifted_rotates_a_reverse_strand_allele_the_other_way() {
        let forward = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"TCA".to_vec());
        assert_eq!(forward.display_allele_shifted(10, 1), "CAT");
        let mut reverse = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"TGA".to_vec());
        reverse.strand = Strand::Reverse;
        assert_eq!(reverse.display_allele_shifted(10, 1), "ATG");
        assert_eq!(reverse.display_allele_shifted(10, -1), "GAT");
    }

    /// `_raw_allele_hashes` under `-use_feature_ref` replaces only the reference
    /// allele's object (VariationFeatureOverlap.pm:418-425): the alternates, the
    /// strand, the class and the record-level fields stay the input's, the
    /// allele string reads the new reference, and a reference row, being the
    /// reference allele itself, becomes the new one.
    #[test]
    fn against_reference_allele_swaps_the_reference_and_keeps_the_rest() {
        let mut given = InputVariant::new("21".into(), 100, 101, b"AC".to_vec(), b"-".to_vec());
        given.strand = Strand::Reverse;
        given.id = Some("rs1".into());
        given.alt_alleles.push(b"GT".to_vec());
        let reread = given.against_reference_allele(b"AT");
        assert_eq!(reread.ref_allele, b"AT");
        assert_eq!(reread.alt_alleles, given.alt_alleles);
        assert_eq!(reread.allele_string, "AT/-/GT");
        assert_eq!(
            (
                reread.strand,
                reread.variant_class,
                reread.id.as_deref(),
                reread.start,
                reread.end
            ),
            (Strand::Reverse, given.variant_class, Some("rs1"), 100, 101)
        );
        let reference_row = given.reference_allele_row().against_reference_allele(b"AT");
        assert_eq!(reference_row.alt_alleles, vec![b"AT".to_vec()]);
        assert_eq!(reference_row.allele_string, "AT/AT");
        assert!(reference_row.reference_allele_row);
    }
}

/// The name Ensembl VEP's VCF parser gives a record from its ID field: the first of the
/// `;`-joined IDs (ensembl-io BaseVCF4.pm:436-440 splits the field, Parser/VCF.pm:347
/// takes `$ids->[0]`), `.` and `0` as written; `None` for an empty field or first
/// element, which Perl's `split` and `||=` read alike as no name. An Ensembl-format ID
/// column is not split (Parser/VEP_input.pm:127 takes the column whole).
pub fn first_vcf_id(field: &str) -> Option<&str> {
    field.split(';').next().filter(|id| !id.is_empty())
}

/// Whether a stored name ([`InputVariant::id`]) names its record in `Uploaded_variation`
/// and may stand in for a breakend's mate: Parser.pm:462 `||=` replaces a Perl-false
/// name (`0`) with the generated one, and OutputFactory.pm:900 rejects `.`.
pub fn id_names_record(id: &str) -> bool {
    perl_true(id) && id != "."
}

/// Classify a variant based on ref and alt alleles.
pub fn classify_variant(ref_allele: &[u8], alt_allele: &[u8]) -> VariantClass {
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    if ref_is_dash && !alt_is_dash {
        VariantClass::Insertion
    } else if !ref_is_dash && alt_is_dash {
        VariantClass::Deletion
    } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
        VariantClass::Snv
    } else if ref_allele.len() == alt_allele.len() {
        VariantClass::Substitution
    } else {
        VariantClass::Indel
    }
}

/// A known variant from the variation cache that co-locates with an input variant.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ColocatedVariant {
    /// The chromosome column of the start-indexed variation file, the file's own
    /// chromosome; `None` for a record of the per-region layout, which has none.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chr: Option<String>,
    /// Variant identifier (e.g., "rs699").
    pub id: String,
    /// Start position.
    pub start: u64,
    /// End position.
    pub end: u64,
    /// Allele string.
    pub allele_string: Option<String>,
    /// Strand.
    pub strand: i8,
    /// Whether this is a somatic variant.
    pub somatic: bool,
    /// Clinical significance terms.
    pub clin_sig: Vec<String>,
    /// Phenotype/disease association flag.
    pub phenotype: bool,
    /// PubMed IDs.
    pub pubmed: Vec<String>,
    /// The input alleles the record's alleles were matched to, as
    /// `Utils/Sequence.pm` `get_matched_variant_alleles` returns them (ensembl-variation
    /// release/116). Empty for a record matched on its span alone, a null-allele record
    /// or any record under `--no_check_alleles`, which prints on every allele's row.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub matched_alleles: Vec<MatchedAllele>,
    /// The record's frequency columns as the cache stores them, by population
    /// column name: `allele:frequency` pairs joined by `,`.
    #[serde(default, skip_serializing_if = "std::collections::HashMap::is_empty")]
    pub frequency_strings: std::collections::HashMap<String, String>,
    /// The record's ClinVar allele and somatic fields, when it carries any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clinvar: Option<Box<crate::variation::ClinVarDetail>>,
    /// The record's synonyms as the cache stores them: `source::id,id` groups
    /// joined by `--`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub var_synonyms: Option<String>,
}

impl ColocatedVariant {
    /// Whether the record prints on the rows of `allele` (the row's `Allele`
    /// column): a record matched on its span alone prints on every row, a record
    /// matched by allele only on the rows of an allele it matched
    /// (`OutputFactory.pm:1081-1083` at release/116.2).
    pub fn applies_to_allele(&self, allele: &str) -> bool {
        self.matched_alleles.is_empty()
            || self
                .matched_alleles
                .iter()
                .any(|matched| matched.input_allele == allele)
    }

    /// The record allele matched to the input allele `allele`, the first match
    /// when the record's alleles matched it more than once (`OutputFactory.pm:1200`).
    pub fn record_allele_for(&self, allele: &str) -> Option<&str> {
        self.matched_alleles
            .iter()
            .find(|matched| matched.input_allele == allele)
            .map(|matched| matched.record_allele.as_str())
    }

    /// Ensembl VEP's rank of the record by the first two letters of its name,
    /// case-insensitively: dbSNP `rs` 1, HGMD `cm`, `ci` and `cd` 2, COSMIC `co`
    /// 3, anything else 100 (`OutputFactory.pm:1057-1065`).
    pub fn prefix_rank(&self) -> u8 {
        let prefix: String = self
            .id
            .chars()
            .take(2)
            .flat_map(char::to_lowercase)
            .collect();
        match prefix.as_str() {
            "rs" => 1,
            "cm" | "ci" | "cd" => 2,
            "co" => 3,
            _ => 100,
        }
    }
}

/// One matched allele pair of an input variant and a known variant: the
/// `a_allele`, `a_index`, `b_allele` and `b_index` of `get_matched_variant_alleles`
/// (`Utils/Sequence.pm:1232-1238`, ensembl-variation release/116).
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct MatchedAllele {
    /// The input allele as the input's allele string spells it.
    pub input_allele: String,
    /// The input allele's position among the input alternate alleles.
    pub input_index: usize,
    /// The known variant's allele as its allele string spells it.
    pub record_allele: String,
    /// The known variant's allele's position among its alternate alleles.
    pub record_index: usize,
}

/// The co-located variants in the order Ensembl VEP reports them: non-somatic
/// variants first, then by [`ColocatedVariant::prefix_rank`], ties in the given
/// order (`OutputFactory.pm:1069-1078`; Perl's `sort` is a stable merge sort).
pub fn colocated_in_output_order(variants: &[ColocatedVariant]) -> Vec<&ColocatedVariant> {
    let mut ordered: Vec<&ColocatedVariant> = variants.iter().collect();
    ordered.sort_by_key(|cv| (cv.somatic, cv.prefix_rank()));
    ordered
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `chr` prefix is removed for the cache lookup and kept for the output, as
    /// `_have_chr` accepts the name through its `chr` rule without rewriting it
    /// (Parser.pm:788-792).
    #[test]
    fn chr_prefix_stripped_for_lookup() {
        let v = InputVariant::new("chr21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(v.chr, "21");
        assert_eq!(v.original_chr, "chr21");
        assert_eq!(v.location(), "chr21:100");
        assert_eq!(v.variation_name(), "chr21_100_A/G");
    }

    /// A record renamed from `M` to `MT` (Parser.pm:794-796, applied by the
    /// pipeline once the cache's names are known) prints `MT` in its Location and
    /// keeps `M` in the names Ensembl VEP built before the rename (Parser.pm:462):
    /// the generated `Uploaded_variation` and the JSON `id`.
    #[test]
    fn bare_m_prints_mt() {
        let mut v = InputVariant::new("M".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(v.chr, "MT");
        v.input_chr = Some("M".into());
        v.original_chr = "MT".to_string();
        assert_eq!(v.location(), "MT:100");
        assert_eq!(v.uploaded_variation(), "M_100_A/G");
        assert_eq!(v.variation_name(), "M_100_A/G");
        // The name an ID is compared with is built at output from the renamed
        // chromosome (OutputFactory.pm:887-891, :900): an ID spelled `M_...`
        // differs from it and prints; one spelled `MT_...` equals it, so the line
        // name stands in for it.
        v.id = Some("M_100_A/G".to_string());
        v.input_line_name = Some("M_100_A/G/T".into());
        assert_eq!(v.uploaded_variation(), "M_100_A/G");
        v.id = Some("MT_100_A/G".to_string());
        assert_eq!(v.uploaded_variation(), "M_100_A/G/T");
        // An ID column Perl reads as false, empty or `0`, takes the name of
        // Parser.pm:462, which differs from the output-time name after the rename
        // and so prints over the line name; a `.` ID takes the line name.
        v.id = None;
        assert_eq!(v.uploaded_variation(), "M_100_A/G");
        v.id = Some("0".to_string());
        assert_eq!(v.uploaded_variation(), "M_100_A/G");
        v.id = Some(".".to_string());
        assert_eq!(v.uploaded_variation(), "M_100_A/G/T");
    }

    /// `chrM` is found through the cache's synonyms or its `chr` rule and keeps
    /// its spelling everywhere, as every name but a bare `M` does.
    #[test]
    fn chrm_keeps_input_name() {
        let v = InputVariant::new("chrM".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(v.chr, "MT");
        assert_eq!(v.original_chr, "chrM");
        assert!(v.input_chr.is_none());
        assert_eq!(v.location(), "chrM:100");
        assert_eq!(v.variation_name(), "chrM_100_A/G");
    }

    #[test]
    fn test_classify_variant() {
        assert_eq!(classify_variant(b"A", b"G"), VariantClass::Snv);
        assert_eq!(classify_variant(b"-", b"TT"), VariantClass::Insertion);
        assert_eq!(classify_variant(b"AC", b"-"), VariantClass::Deletion);
        assert_eq!(classify_variant(b"AC", b"GT"), VariantClass::Substitution);
        assert_eq!(classify_variant(b"A", b"TT"), VariantClass::Indel);
    }

    #[test]
    fn test_variant_location() {
        let snp = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert_eq!(snp.location(), "21:25000100");

        let del = InputVariant::new(
            "21".into(),
            25000100,
            25000105,
            b"ACGTAC".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(del.location(), "21:25000100-25000105");

        let ins = InputVariant::new(
            "21".into(),
            25000101,
            25000100,
            b"-".to_vec(),
            b"TT".to_vec(),
        );
        assert_eq!(ins.location(), "21:25000100-25000101");
    }

    /// `Utils/Sequence.pm` `ambiguity_code` (ensembl-variation release/116,
    /// lines 115-130): the fifteen single-base sets and `-` of its `%ambig` table,
    /// keyed by the distinct alleles joined, with `N` standing for the four bases;
    /// every other key, so every anchor-trimmed indel (`-/T`, `TA/-`) and any
    /// record carrying a `*` allele, has no code. `VariationFeature::ambig_code`
    /// reads the record's whole allele string, so a multi-allelic record takes
    /// the code of all its alleles, a bi-allelic indel its minimised string, and
    /// a structural variant never has one (`OutputFactory.pm`
    /// `VariationFeature_to_output_hash` asks a `VariationFeature` alone).
    #[test]
    fn ambiguity_code_follows_perl_table() {
        let table = [
            ("A/C", b'M'),
            ("A/C/G", b'V'),
            ("A/C/G/T", b'N'),
            ("A/C/T", b'H'),
            ("A/G", b'R'),
            ("A/G/T", b'D'),
            ("A/T", b'W'),
            ("C/G", b'S'),
            ("C/G/T", b'B'),
            ("C/T", b'Y'),
            ("G/T", b'K'),
            ("C", b'C'),
            ("A", b'A'),
            ("T", b'T'),
            ("G", b'G'),
            ("-", b'-'),
        ];
        for (alleles, code) in table {
            assert_eq!(ambiguity_code(alleles), Some(code), "{alleles}");
        }
        // Order, repeats and the separator do not matter; `N` expands.
        assert_eq!(ambiguity_code("T/A"), Some(b'W'));
        assert_eq!(ambiguity_code("A/G/A"), Some(b'R'));
        assert_eq!(ambiguity_code("A|G"), Some(b'R'));
        assert_eq!(ambiguity_code("N/A"), Some(b'N'));
        assert_eq!(ambiguity_code("a/c"), Some(b'M'));
        // The key is the joined allele strings, so a multi-base allele whose
        // letters complete a set has a code (Ensembl VEP 116.2 prints these).
        assert_eq!(ambiguity_code("G/AC"), Some(b'V'));
        assert_eq!(ambiguity_code("T/AC"), Some(b'H'));
        assert_eq!(ambiguity_code("T/AG"), Some(b'D'));
        assert_eq!(ambiguity_code("C/GT"), Some(b'B'));
        assert_eq!(ambiguity_code("A/CGT"), Some(b'N'));
        assert_eq!(ambiguity_code("A/GT"), Some(b'D'));
        assert_eq!(ambiguity_code("C/A/GT"), Some(b'N'));
        for none in [
            "A/AT", "G/AT", "-/T", "T/-", "AG/A/-", "A/*", "A/G/*", "AT/GC", "<DEL>", "",
        ] {
            assert_eq!(ambiguity_code(none), None, "{none}");
        }

        let snv = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(snv.ambiguity_code(), Some(b'R'));
        let mut multi = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        multi.record_allele_string_multi = Some("A/G/T".into());
        assert_eq!(multi.ambiguity_code(), Some(b'D'));
        let ins = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"T".to_vec());
        assert_eq!(ins.ambiguity_code(), None);
        let mut sv = InputVariant::new("21".into(), 100, 200, b"N".to_vec(), b"<DEL>".to_vec());
        sv.is_structural = true;
        sv.uploaded_allele_string = Some("<DEL>".into());
        assert_eq!(sv.ambiguity_code(), None);
        let mut sv_like_snv =
            InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        sv_like_snv.is_structural = true;
        assert_eq!(sv_like_snv.ambiguity_code(), None);
    }

    /// `Utils/Sequence.pm` `SO_variation_class` (ensembl-variation release/116,
    /// lines 335-412) on the allele strings its own comments name and on the
    /// record strings a multi-allelic VCF line becomes: the class is the whole
    /// record's, so two deletions of different length are a
    /// `sequence_alteration`, a mixed SNV and insertion an `indel`, and a `*`
    /// allele defeats every pattern (`OutputFactory.pm`
    /// `VariationFeature_to_output_hash` writes `class_SO_term`).
    #[test]
    fn so_variation_class_follows_perl() {
        for (alleles, class) in [
            ("A/T", "SNV"),
            ("A/T/G", "SNV"),
            ("N/A", "SNV"),
            ("(CAG)8/(CAG)9", "tandem_repeat"),
            ("-/ATT", "insertion"),
            ("-/(LARGEINSERTION)", "insertion"),
            ("-/(LARGEDELETION)", "deletion"),
            ("A/-", "deletion"),
            ("TA/-", "deletion"),
            ("A/(LARGEDELETION)", "deletion"),
            ("AA/TT", "substitution"),
            ("AA/TTTT", "indel"),
            ("(LARGEDELETION)/-", "deletion"),
            ("(LARGEINSERTION)", "insertion"),
            ("(308 BP DELETION)", "deletion"),
            ("A/C/*", "sequence_alteration"),
            ("TG/T/-", "sequence_alteration"),
            ("T/A/TT", "indel"),
            ("-/T/-", "sequence_alteration"),
            ("a/t", "sequence_alteration"),
        ] {
            assert_eq!(so_variation_class(alleles), class, "{alleles}");
        }
        let snv = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(snv.class_so_term(), "SNV");
        let lower = InputVariant::new("21".into(), 100, 100, b"a".to_vec(), b"g".to_vec());
        assert_eq!(lower.class_so_term(), "SNV");
        let lower_del = InputVariant::new("21".into(), 100, 100, b"t".to_vec(), b"-".to_vec());
        assert_eq!(lower_del.class_so_term(), "deletion");
        let mut multi = InputVariant::new("21".into(), 101, 102, b"TG".to_vec(), b"T".to_vec());
        multi.record_allele_string_multi = Some("TG/T/-".into());
        assert_eq!(multi.class_so_term(), "sequence_alteration");
        // A resolved identifier prints the class its record carries, the stored
        // one, where its minimised alleles would read as a deletion.
        let mut stored = InputVariant::new("21".into(), 102, 102, b"G".to_vec(), b"-".to_vec());
        stored.record_so_term = Some("indel".into());
        assert_eq!(stored.class_so_term(), "indel");
        let mut sv = InputVariant::new("21".into(), 100, 200, b"N".to_vec(), b"<DEL>".to_vec());
        sv.is_structural = true;
        sv.variant_class = VariantClass::StructuralDeletion;
        sv.uploaded_allele_string = Some("<DEL>".into());
        assert_eq!(sv.class_so_term(), "deletion");
        // A structural record prints the term its own ALT names, as the Allele
        // column does: Ensembl VEP writes `Alu_insertion`, not the class term.
        let mut alu = InputVariant::new(
            "21".into(),
            101,
            100,
            b"N".to_vec(),
            b"<INS:ME:ALU>".to_vec(),
        );
        alu.is_structural = true;
        alu.variant_class = VariantClass::MobileElementInsertion;
        alu.record_so_term = Some("Alu_insertion".into());
        assert_eq!(alu.class_so_term(), "Alu_insertion");
        // A type the term table lacks has no record term; its row prints Perl's
        // type string in the Allele and VARIANT_CLASS columns (Parser/VCF.pm:469-481):
        // the ALT as written, or INFO/SVTYPE when that names the type.
        let mut cpx = InputVariant::new("21".into(), 101, 300, b"N".to_vec(), b"<CPX>".to_vec());
        cpx.is_structural = true;
        cpx.variant_class = VariantClass::ComplexStructural;
        cpx.vep_skip = true;
        cpx.uploaded_allele_string = Some("<CPX>".into());
        assert_eq!(cpx.class_so_term(), "<CPX>");
        assert_eq!(cpx.class_so_term(), cpx.display_allele());
        cpx.sv_type = Some("CPX".into());
        assert_eq!(cpx.class_so_term(), "CPX");
        assert_eq!(cpx.display_allele(), "CPX");
        cpx.sv_type = Some("FOO".into());
        assert_eq!(cpx.class_so_term(), "FOO");
        let mut ctx = cpx.clone();
        ctx.alt_alleles = vec![b"<CTX>".to_vec()];
        ctx.uploaded_allele_string = Some("<CTX>".into());
        ctx.sv_type = None;
        assert_eq!(ctx.class_so_term(), "<CTX>");
        // A breakend row's Allele is its breakend string; its class is the term.
        let mut bnd =
            InputVariant::new("21".into(), 100, 100, b"N".to_vec(), b"N[21:200[".to_vec());
        bnd.is_structural = true;
        bnd.variant_class = VariantClass::Translocation;
        bnd.mate_id = Some("bnd_1".into());
        assert_eq!(bnd.class_so_term(), "chromosome_breakpoint");
        assert_eq!(bnd.display_allele(), "N[21:200[");
    }

    /// `UPLOADED_ALLELE` is the raw string of a record whose alleles differ in
    /// length with no `-` among them, or of every record under `--minimal`, and
    /// the upper-cased working string otherwise (`OutputFactory.pm`
    /// `VariationFeatureOverlapAllele_to_output_hash` reads
    /// `original_allele_string || allele_string`; `Parser.pm` `post_process_vfs`
    /// sets the former when a length differs, ending its test at the first `-`,
    /// `minimise_alleles` sets it for every record, and `validate_vf` upper-cases
    /// the latter), and a structural row has none. An Ensembl-format column
    /// written `at/-` or `tt/-/t` therefore prints upper-cased, as `t/-` does, and
    /// a record whose reference `--lookup_ref` replaced prints the replaced string
    /// (`A/A/G` for `T A,G` over a reference `A`, as Ensembl VEP 116.2 prints).
    #[test]
    fn uploaded_allele_follows_perls_two_sources() {
        for (raw, kept) in [
            ("a/g", false),
            ("A/C/*", false),
            ("at/a", true),
            ("TAG/TA/T", true),
            ("AT/A/*", true),
            ("t/-", false),
            ("at/-", false),
            ("-/at", false),
            ("tt/-/t", false),
        ] {
            assert_eq!(keeps_uploaded_spelling(raw, false), kept, "{raw}");
            assert!(keeps_uploaded_spelling(raw, true), "{raw} under --minimal");
        }
        let mut snv = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        snv.uploaded_allele_string = Some("a/g".into());
        assert_eq!(snv.uploaded_allele(false).as_deref(), Some("A/G"));
        assert_eq!(snv.uploaded_allele(true).as_deref(), Some("a/g"));
        let plain = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(plain.uploaded_allele(false).as_deref(), Some("A/G"));
        let mut looked_up = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        looked_up.allele_index = 1;
        looked_up.uploaded_allele_string = Some("T/A/G".into());
        looked_up.record_allele_string_multi = Some("A/A/G".into());
        assert_eq!(looked_up.uploaded_allele(false).as_deref(), Some("A/A/G"));
        let mut sv = InputVariant::new("21".into(), 100, 200, b"N".to_vec(), b"<DEL>".to_vec());
        sv.is_structural = true;
        sv.uploaded_allele_string = Some("<DEL>".into());
        assert_eq!(sv.uploaded_allele(false), None);
    }

    /// A record `--check_frequency` rewrote works with its kept alleles: the
    /// JSON `allele_string`, `AMBIGUITY`, `VARIANT_CLASS` and, where Ensembl VEP
    /// prints the feature's string, `UPLOADED_ALLELE` read `T/A` for
    /// `21 7759128 T A,G` after the `G` allele is dropped, while the names built
    /// before the filter ran (`variation_name`, the line name) keep the uploaded
    /// one; `21 44195140 G A,GT` keeps its uploaded `UPLOADED_ALLELE`, as a
    /// non-minimised indel does, and prints `indel`. The class is the rewritten
    /// string's because the filtering source annotates first
    /// (AnnotationSourceAdaptor.pm:98) and `class_SO_term` reads the feature's
    /// string (VariationFeature.pm:1118): Ensembl VEP 116.2 prints
    /// `VARIANT_CLASS=SNV` on every row of `21 31148973 GG G,GC` once the exclude
    /// filter drops its common `-` allele, `sequence_alteration` without the
    /// filter. A record whose ID Perl reads as false is named by Parser.pm:462
    /// before the filter, and `Uploaded_variation` prints that name once the
    /// rewrite makes the output's `chr_start_alleles` differ from it
    /// (OutputFactory.pm:887-900); an explicit ID equal to the output's name
    /// yields the line name, as without the filter.
    #[test]
    fn a_frequency_filtered_record_prints_its_kept_alleles() {
        let mut snv =
            InputVariant::new("21".into(), 7759128, 7759128, b"T".to_vec(), b"A".to_vec());
        snv.uploaded_allele_string = Some("T/A/G".into());
        snv.record_allele_string_multi = Some("T/A/G".into());
        snv.input_line_name = Some("21_7759128_T/A/G".into());
        assert_eq!(snv.ambiguity_code(), Some(b'D'));
        snv.filtered_allele_string = Some("T/A".into());
        assert_eq!(snv.working_allele_string(), "T/A");
        assert_eq!(snv.record_allele_string(), "T/A/G");
        assert_eq!(snv.uploaded_allele(false).as_deref(), Some("T/A"));
        assert_eq!(snv.uploaded_allele(true).as_deref(), Some("T/A/G"));
        assert_eq!(snv.ambiguity_code(), Some(b'W'));
        assert_eq!(snv.variation_name(), "21_7759128_T/A/G");
        assert_eq!(snv.class_so_term(), "SNV");
        snv.id = Some(".".into());
        assert_eq!(snv.uploaded_variation(), "21_7759128_T/A/G");

        let mut indel = InputVariant::new(
            "21".into(),
            44195140,
            44195140,
            b"G".to_vec(),
            b"GT".to_vec(),
        );
        indel.uploaded_allele_string = Some("G/A/GT".into());
        indel.record_allele_string_multi = Some("G/A/GT".into());
        indel.filtered_allele_string = Some("G/GT".into());
        assert_eq!(indel.working_allele_string(), "G/GT");
        assert_eq!(indel.uploaded_allele(false).as_deref(), Some("G/A/GT"));
        assert_eq!(indel.class_so_term(), "indel");

        // `21 100 T A,TT` with the `TT` allele dropped: the class is the kept
        // `T/A`'s `SNV`, not the uploaded string's `indel`.
        let mut mixed = InputVariant::new("21".into(), 100, 100, b"T".to_vec(), b"A".to_vec());
        mixed.uploaded_allele_string = Some("T/A/TT".into());
        mixed.record_allele_string_multi = Some("T/A/TT".into());
        assert_eq!(mixed.class_so_term(), "indel");
        mixed.filtered_allele_string = Some("T/A".into());
        assert_eq!(mixed.uploaded_allele(false).as_deref(), Some("T/A/TT"));
        assert_eq!(mixed.class_so_term(), "SNV");

        // `21 31148973 GG G,GC`, anchor-chopped to `G/-/C` at 31148974, after
        // the common `-` allele is dropped.
        let mut del_plus_snv = InputVariant::new(
            "21".into(),
            31148974,
            31148974,
            b"G".to_vec(),
            b"C".to_vec(),
        );
        del_plus_snv.record_allele_string_multi = Some("G/-/C".into());
        assert_eq!(del_plus_snv.class_so_term(), "sequence_alteration");
        del_plus_snv.filtered_allele_string = Some("G/C".into());
        assert_eq!(del_plus_snv.class_so_term(), "SNV");

        // `21 100 <ID> AT A,ATT`, anchor-chopped to `T/-/TT` at 101, after the
        // `TT` allele is dropped.
        let mut chopped = InputVariant::new("21".into(), 101, 101, b"T".to_vec(), b"-".to_vec());
        chopped.uploaded_allele_string = Some("AT/A/ATT".into());
        chopped.record_allele_string_multi = Some("T/-/TT".into());
        chopped.input_line_name = Some("21_100_AT/A/ATT".into());
        for id in [None, Some("."), Some("0"), Some("21_101_T/-/TT")] {
            chopped.id = id.map(str::to_string);
            assert_eq!(chopped.uploaded_variation(), "21_100_AT/A/ATT", "ID {id:?}");
        }
        chopped.filtered_allele_string = Some("T/-".into());
        assert_eq!(chopped.uploaded_allele(false).as_deref(), Some("AT/A/ATT"));
        for (id, name) in [
            (None, "21_101_T/-/TT"),
            (Some("0"), "21_101_T/-/TT"),
            (Some("."), "21_100_AT/A/ATT"),
            (Some("21_101_T/-/TT"), "21_101_T/-/TT"),
            (Some("21_101_T/-"), "21_100_AT/A/ATT"),
            (Some("rs1"), "rs1"),
        ] {
            chopped.id = id.map(str::to_string);
            assert_eq!(chopped.uploaded_variation(), name, "ID {id:?}");
        }
        chopped.id = None;
        assert_eq!(chopped.variation_name(), "21_101_T/-/TT");
    }

    /// The three byte writers append exactly the bytes of the `String` columns,
    /// on every branch each column has: the location of a point, a span and an
    /// insertion; a record named by its ID, by its input line and by its
    /// alleles at the trimmed start; the allele of a sequence variant, an
    /// allele of a `--minimal` split (the trimmed allele under the record's
    /// span), every structural class term and a breakend with and without a
    /// mate.
    #[test]
    fn byte_writers_match_the_string_columns() {
        let mut named = InputVariant::new("chr21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        named.id = Some("rs123".into());
        let del = InputVariant::new(
            "21".into(),
            25000100,
            25000105,
            b"ACGTAC".to_vec(),
            b"-".to_vec(),
        );
        let mut ins = InputVariant::new(
            "21".into(),
            25000101,
            25000100,
            b"-".to_vec(),
            b"TT".to_vec(),
        );
        ins.input_line_name = Some("21_25000100_A/ATT".into());
        let mut multi = InputVariant::new("21".into(), 101, 101, b"AG".to_vec(), b"A".to_vec());
        multi.original_start = Some(100);
        multi.original_end = Some(102);
        multi.original_allele_string = Some("TAG/TA".into());
        multi.uploaded_allele_string = Some("TAG/TA/T".into());
        let mut cpx = InputVariant::new("21".into(), 100, 300, b"N".to_vec(), b"<CPX:X>".to_vec());
        cpx.variant_class = VariantClass::ComplexStructural;
        cpx.is_structural = true;
        let mut alu = make_me_variant(b"<INS:ME:ALU>");
        alu.record_so_term = Some("Alu_insertion".into());
        let bare_me = make_me_variant(b"<INS:ME:SVA>");
        let mut non_ref =
            InputVariant::new("21".into(), 100, 300, b"N".to_vec(), b"<NON_REF>".to_vec());
        non_ref.variant_class = VariantClass::ReferenceBlock;
        non_ref.is_structural = true;
        let mut paired =
            InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"A[21:200[".to_vec());
        paired.variant_class = VariantClass::Translocation;
        paired.is_structural = true;
        paired.mate_id = Some("bnd_1".into());
        let mut single = paired.clone();
        single.mate_id = None;
        let mut dup = InputVariant::new("21".into(), 100, 300, b"N".to_vec(), b"<DUP>".to_vec());
        dup.variant_class = VariantClass::Duplication;
        dup.is_structural = true;
        dup.uploaded_allele_string = Some("<DUP>".into());

        let expected = [
            (&named, "rs123", "chr21:100", "G"),
            (&del, "21_25000100_ACGTAC/-", "21:25000100-25000105", "-"),
            (&ins, "21_25000100_A/ATT", "21:25000100-25000101", "TT"),
            (&multi, "21_100_TAG/TA/T", "21:100-102", "A"),
            (&cpx, "21_100_N/<CPX:X>", "21:100-300", "<CPX:X>"),
            (&alu, "1_100_N/<INS:ME:ALU>", "1:100-200", "Alu_insertion"),
            (
                &bare_me,
                "1_100_N/<INS:ME:SVA>",
                "1:100-200",
                "SVA_insertion",
            ),
            (&non_ref, "21_100_N/<NON_REF>", "21:100-300", "<NON_REF>"),
            (&paired, "1_100_N/A[21:200[", "1:100", "A[21:200["),
            (
                &single,
                "1_100_N/A[21:200[",
                "1:100",
                "chromosome_breakpoint",
            ),
            (&dup, "21_100_<DUP>", "21:100-300", "duplication"),
        ];
        for (v, uploaded, location, allele) in expected {
            let mut out = Vec::new();
            v.write_uploaded_variation(&mut out);
            assert_eq!(out, uploaded.as_bytes(), "{}", v.allele_string);
            assert_eq!(v.uploaded_variation(), uploaded);
            out.clear();
            v.write_location(&mut out);
            assert_eq!(out, location.as_bytes(), "{}", v.allele_string);
            assert_eq!(v.location(), location);
            out.clear();
            v.write_display_allele(&mut out);
            assert_eq!(out, allele.as_bytes(), "{}", v.allele_string);
            assert_eq!(v.display_allele(), allele);
        }
    }

    /// Zero, every value to 1000, each power of ten with its neighbours and the
    /// 20-digit `u64::MAX` write as `Display` does.
    #[test]
    fn push_decimal_writes_every_digit() {
        let mut cases: Vec<u64> = (0..=1000).collect();
        let mut power_of_ten: u64 = 1;
        for _ in 0..19 {
            power_of_ten = power_of_ten.saturating_mul(10);
            cases.extend_from_slice(&[power_of_ten - 1, power_of_ten, power_of_ten + 1]);
        }
        cases.extend_from_slice(&[25000100, 12_345_678_901_234_567_890, u64::MAX - 1, u64::MAX]);
        for n in cases {
            let mut out = b"x".to_vec();
            push_decimal(&mut out, n);
            assert_eq!(out, format!("x{n}").into_bytes());
        }
    }

    #[test]
    fn test_variant_uploaded_variation() {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(v.uploaded_variation(), "21_100_A/G");
        assert_eq!(v.variation_name(), "21_100_A/G");

        // An Ensembl-format ID column is kept whole (Parser/VEP_input.pm:127).
        v.id = Some("rs123;rs456".into());
        assert_eq!(v.uploaded_variation(), "rs123;rs456");
        assert_eq!(v.variation_name(), "rs123;rs456");

        v.id = Some(".".into());
        assert_eq!(v.uploaded_variation(), "21_100_A/G");
        assert_eq!(v.variation_name(), ".");
    }

    /// Ensembl VEP names a record by its line when its stored name is `.`
    /// (OutputFactory.pm:900), absent or `0` (Parser.pm:462 `||=` replaces a Perl-false
    /// name with `chr_start_alleles`, which OutputFactory.pm:900 then rejects as the
    /// generated name), and its JSON `id` (OutputFactory/JSON.pm:236) is that stored name
    /// after the `||=`. On `21 100 . AT A`, minimised to `T/-` at 101, the names `.`,
    /// none and `0` all yield the line name `21_100_AT/A`; the JSON `id` is `.` for `.`
    /// and the generated `21_101_T/-` for none and `0`.
    #[test]
    fn a_perl_false_or_dot_name_yields_the_line_name_and_the_generated_json_id() {
        let mut v = InputVariant::new("21".into(), 101, 101, b"T".to_vec(), b"-".to_vec());
        v.uploaded_allele_string = Some("AT/A".into());
        v.input_line_name = Some("21_100_AT/A".into());
        v.minimised = true;
        for (id, json_id) in [
            (Some("."), "."),
            (None, "21_101_T/-"),
            (Some("0"), "21_101_T/-"),
        ] {
            v.id = id.map(str::to_string);
            assert_eq!(v.uploaded_variation(), "21_100_AT/A", "ID {id:?}");
            assert_eq!(v.variation_name(), json_id, "ID {id:?}");
        }
    }

    /// `first_vcf_id` reads the ID field as ensembl-io BaseVCF4.pm:436-440 and
    /// Parser/VCF.pm:347 do: the first `;`-joined element, `.` and `0` kept, an empty
    /// field or element no name.
    #[test]
    fn first_vcf_id_takes_the_first_element_as_written() {
        assert_eq!(first_vcf_id("rs2;."), Some("rs2"));
        assert_eq!(first_vcf_id(".;rs2"), Some("."));
        assert_eq!(first_vcf_id("0;rs2"), Some("0"));
        assert_eq!(first_vcf_id("0"), Some("0"));
        assert_eq!(first_vcf_id("."), Some("."));
        assert_eq!(first_vcf_id(";rs2"), None);
        assert_eq!(first_vcf_id(""), None);
    }

    #[test]
    fn uploaded_variation_names_a_minimised_indel_by_its_raw_alleles() {
        // `21 100 . AT A`: VEP minimises to T/- at 101 and names the record by the raw
        // alleles at the trimmed start; its JSON id keeps the VCF's `.`.
        let mut v = InputVariant::new("21".into(), 101, 101, b"T".to_vec(), b"-".to_vec());
        v.id = Some(".".into());
        v.uploaded_allele_string = Some("AT/A".into());
        v.minimised = true;
        assert_eq!(v.uploaded_variation(), "21_101_AT/A");
        assert_eq!(v.variation_name(), ".");
        assert_eq!(v.record_allele_string(), "T/-");
        // A record from a format without an ID column is named by its alleles.
        v.id = None;
        v.uploaded_allele_string = None;
        assert_eq!(v.variation_name(), "21_101_T/-");
    }

    #[test]
    fn multi_allelic_and_structural_records_share_one_name() {
        // `21 100 . TAG TA,T`: named by the raw alleles at the chopped start; the
        // JSON allele_string is the chopped working string.
        let mut v = InputVariant::new("21".into(), 101, 101, b"AG".to_vec(), b"A".to_vec());
        v.id = Some(".".into());
        v.uploaded_allele_string = Some("TAG/TA/T".into());
        v.record_allele_string_multi = Some("AG/A/-".into());
        assert_eq!(v.uploaded_variation(), "21_101_TAG/TA/T");
        assert_eq!(v.variation_name(), ".");
        assert_eq!(v.record_allele_string(), "AG/A/-");

        let mut sv = InputVariant::new(
            "21".into(),
            46222496,
            46224000,
            b"T".to_vec(),
            b"<CN0>".to_vec(),
        );
        sv.is_structural = true;
        sv.uploaded_allele_string = Some("<CN0>".into());
        assert_eq!(sv.uploaded_variation(), "21_46222496_<CN0>");
        assert_eq!(sv.record_allele_string(), "<CN0>");
    }

    /// The reference row is the reference of the record's working string, numbered
    /// 0 over the record's span, and shares the record's `minimised` mark
    /// (`VariationFeature_to_output_hash` sets `MINIMISED` from the feature) and
    /// line; for an allele of a `--minimal` split it is the allele's own trimmed
    /// reference over the trimmed span (InputBuffer.pm:674-676), with the record's
    /// span restored for Location and the allele's rejoin key.
    #[test]
    fn reference_row_shares_the_records_fields_and_a_split_alleles_span() {
        let mut minimised = InputVariant::new("21".into(), 101, 102, b"AA".to_vec(), b"-".to_vec());
        minimised.uploaded_allele_string = Some("AAA/A".into());
        minimised.original_allele_string = Some("AAA/A".into());
        minimised.minimised = true;
        minimised.raw_input = Some("21\t100\t.\tAAA\tA".into());
        minimised.raw_input_bytes = Some(Arc::from(b"21\t100\t.\tAAA\tA".as_slice()));
        let row = minimised.reference_allele_row();
        assert!(row.reference_allele_row && row.minimised);
        assert_eq!(
            (row.ref_allele.as_slice(), row.alt_allele()),
            (b"AA".as_slice(), b"AA".as_slice())
        );
        assert_eq!((row.start, row.end, row.allele_num), (101, 102, Some(0)));
        assert_eq!(row.original_allele_string.as_deref(), Some("AAA/A"));
        assert_eq!(row.raw_input_bytes, minimised.raw_input_bytes);
        assert_eq!(row.record_allele_string(), "AA/-");

        // `21 100 . AAA AGG,AAG` under `--minimal`: the second allele trims to `A/G`
        // at 102 and its rows sort under `A/G_`.
        let mut second = InputVariant::new("21".into(), 100, 102, b"AAA".to_vec(), b"AAG".to_vec());
        second.allele_index = 1;
        second.record_allele_string_multi = Some("AAA/AGG/AAG".into());
        assert_eq!(
            crate::minimal::apply_minimal_split(&mut second, 1),
            crate::minimal::MinimalSplit::Trimmed
        );
        let row = second.reference_allele_row();
        assert_eq!(
            (row.ref_allele.as_slice(), row.alt_allele()),
            (b"A".as_slice(), b"A".as_slice())
        );
        assert_eq!((row.start, row.end), (102, 102));
        assert_eq!(
            (row.original_start, row.original_end),
            (Some(100), Some(102))
        );
        assert_eq!(row.rejoin_key_prefix.as_deref(), Some("A/G_"));
        assert!(row.minimised);
        assert_eq!(row.location(), "21:100-102");
        assert_eq!(row.record_allele_string(), "AAA/AGG/AAG");
        assert_eq!(row.uploaded_allele(true).as_deref(), Some("AAA/AGG/AAG"));
        assert_eq!(row.class_so_term(), "substitution");
    }

    /// The reference row names the alternate whose 3' shift places it under
    /// `--shift_3prime 1`, the allele it is built from, so the shift reads the
    /// row at that allele's shifted span; the row of a feature with no alternate
    /// allele names none and is read where it stands.
    #[test]
    fn reference_row_names_the_alternate_it_shifts_with() {
        let insertion = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"GAA".to_vec());
        let row = insertion.reference_allele_row();
        assert_eq!(row.reference_row_shift_allele.as_deref(), Some("GAA"));
        assert_eq!(
            (row.ref_allele.as_slice(), row.alt_allele()),
            (b"-".as_slice(), b"-".as_slice())
        );
        let without_allele = InputVariant::without_allele("21".into(), 100, "A", ["*"]);
        assert!(without_allele.no_allele);
        assert_eq!(
            without_allele
                .reference_allele_row()
                .reference_row_shift_allele,
            None
        );
    }

    #[test]
    fn test_sv_so_terms_match_ensembl_vep() {
        // SO terms must match Ensembl VEP's `%SO_TERMS` exactly.
        assert_eq!(VariantClass::StructuralDeletion.so_term(), "deletion");
        assert_eq!(VariantClass::StructuralInsertion.so_term(), "insertion");
        assert_eq!(VariantClass::Duplication.so_term(), "duplication");
        assert_eq!(
            VariantClass::TandemDuplication.so_term(),
            "tandem_duplication"
        );
        assert_eq!(VariantClass::Inversion.so_term(), "inversion");
        assert_eq!(
            VariantClass::CopyNumberVariation.so_term(),
            "copy_number_variation"
        );
        assert_eq!(
            VariantClass::MobileElementInsertion.so_term(),
            "mobile_element_insertion"
        );
        assert_eq!(
            VariantClass::Translocation.so_term(),
            "chromosome_breakpoint"
        );
        assert_eq!(VariantClass::TandemRepeat.so_term(), "tandem_repeat");
        assert_eq!(
            VariantClass::ComplexStructural.so_term(),
            "complex_structural_alteration"
        );
        // A type absent from `%SO_TERMS` keeps the type string as its class_SO_term.
        assert_eq!(VariantClass::ReferenceBlock.so_term(), "<NON_REF>");
    }

    #[test]
    fn test_is_structural_true_for_sv_types() {
        let sv_types = [
            VariantClass::StructuralDeletion,
            VariantClass::StructuralInsertion,
            VariantClass::Duplication,
            VariantClass::TandemDuplication,
            VariantClass::Inversion,
            VariantClass::CopyNumberVariation,
            VariantClass::MobileElementInsertion,
            VariantClass::MobileElementDeletion,
            VariantClass::Translocation,
            VariantClass::TandemRepeat,
            VariantClass::ComplexStructural,
            VariantClass::ReferenceBlock,
        ];
        for vc in sv_types {
            assert!(vc.is_structural(), "{vc:?} should be structural");
            assert!(vc.is_span_type(), "{vc:?} starts at POS+1");
        }
    }

    #[test]
    fn test_is_structural_false_for_small_variants() {
        let small_types = [
            VariantClass::Snv,
            VariantClass::Insertion,
            VariantClass::Deletion,
            VariantClass::Indel,
            VariantClass::Substitution,
            VariantClass::SequenceAlteration,
        ];
        for vc in small_types {
            assert!(!vc.is_structural(), "{vc:?} should not be structural");
        }
    }

    #[test]
    fn test_input_variant_sv_fields_default() {
        let v = InputVariant::new("1".into(), 100, 500, b"N".to_vec(), b"-".to_vec());
        assert_eq!(v.sv_end, None);
        assert_eq!(v.sv_type, None);
        assert_eq!(v.sv_len, None);
        assert!(!v.is_structural);
        assert_eq!(v.mate_id, None);
    }

    #[test]
    fn test_input_variant_with_sv_fields() {
        let mut v = InputVariant::new("1".into(), 100, 500, b"N".to_vec(), b"-".to_vec());
        v.variant_class = VariantClass::StructuralDeletion;
        v.sv_end = Some(500);
        v.sv_type = Some("DEL".into());
        v.sv_len = Some(-400);
        v.is_structural = true;
        v.mate_id = None;

        assert_eq!(v.variant_class.so_term(), "deletion");
        assert!(v.variant_class.is_structural());
        assert_eq!(v.sv_end, Some(500));
        assert_eq!(v.sv_type.as_deref(), Some("DEL"));
        assert_eq!(v.sv_len, Some(-400));
        assert!(v.is_structural);
    }

    #[test]
    fn test_display_allele_native_single_breakend_uses_so_term() {
        let mut v = InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"A.".to_vec());
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        assert_eq!(v.display_allele(), "chromosome_breakpoint");
    }

    /// An allele of a record with parsed breakends carries `mate_id` and prints its text;
    /// a bracket Perl's regex refuses (`N[21:abc[`) leaves the record without breakends,
    /// and Ensembl VEP 116.2 prints `chromosome_breakpoint` for it.
    #[test]
    fn test_display_allele_paired_breakend_variants_keep_raw_allele() {
        let mut paired =
            InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"A[21:200[".to_vec());
        paired.variant_class = VariantClass::Translocation;
        paired.is_structural = true;
        paired.mate_id = Some("bnd_1".into());
        assert_eq!(paired.display_allele(), "A[21:200[");
        let mut refused = paired.clone();
        refused.alt_alleles = vec![b"N[21:abc[".to_vec()];
        refused.mate_id = None;
        refused.is_single_breakend = true;
        assert_eq!(refused.display_allele(), "chromosome_breakpoint");

        let mut derived_single =
            InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"A.".to_vec());
        derived_single.variant_class = VariantClass::Translocation;
        derived_single.is_structural = true;
        derived_single.is_single_breakend = true;
        derived_single.mate_id = Some("bnd_2".into());
        assert_eq!(derived_single.display_allele(), "A.");
    }

    /// Parser/VCF.pm `create_StructuralVariationFeatures` gives a bracket
    /// breakend record the one allele string `REF/ALT` (`N/N[21:9412900[`), which
    /// OutputFactory.pm `VariationFeature_to_output_hash` builds the record's
    /// generated name from; the parser stores it as `record_allele_string_multi`, so the
    /// record's own breakend (`N.`) answers with the same string as its mate, and a
    /// native single breakend (`.N`) stays bare.
    #[test]
    fn derived_own_breakend_shares_the_record_allele_string() {
        let mut own =
            InputVariant::new("21".into(), 9412501, 9412501, b"N".to_vec(), b"N.".to_vec());
        own.variant_class = VariantClass::Translocation;
        own.is_structural = true;
        own.is_single_breakend = true;
        own.mate_id = Some("21:9412501".into());
        own.uploaded_allele_string = Some("N[21:9412900[".into());
        own.record_allele_string_multi = Some("N/N[21:9412900[".into());
        assert_eq!(own.record_allele_string(), "N/N[21:9412900[");

        let mut mate = own.clone();
        mate.alt_alleles = vec![b"N[21:9412900[".to_vec()];
        mate.allele_string = "N/N[21:9412900[".into();
        mate.is_single_breakend = false;
        assert_eq!(mate.record_allele_string(), "N/N[21:9412900[");

        let mut native =
            InputVariant::new("21".into(), 9412701, 9412701, b"N".to_vec(), b".N".to_vec());
        native.variant_class = VariantClass::Translocation;
        native.is_structural = true;
        native.is_single_breakend = true;
        native.uploaded_allele_string = Some(".N".into());
        assert_eq!(native.record_allele_string(), ".N");
    }

    /// OutputFactory.pm `BaseStructuralVariationOverlapAllele_to_output_hash` prints
    /// `$svf->class_SO_term`, the record's term, on every row: a row carrying the
    /// record's `Alu_insertion` prints it whatever its own ALT reads, and a `<DEL>`
    /// of a `copy_number_variation` record prints the record's term.
    #[test]
    fn record_so_term_labels_every_structural_row() {
        let mut v = make_me_variant(b"<INS:ME:LINE1>");
        v.record_so_term = Some("Alu_insertion".into());
        assert_eq!(v.display_allele(), "Alu_insertion");

        let mut v = InputVariant::new("21".into(), 100, 300, b"N".to_vec(), b"<DEL>".to_vec());
        v.variant_class = VariantClass::CopyNumberVariation;
        v.is_structural = true;
        v.record_so_term = Some("copy_number_variation".into());
        assert_eq!(v.display_allele(), "copy_number_variation");
    }

    #[test]
    fn test_input_variant_bnd_with_mate_id() {
        let mut v = InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"-".to_vec());
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.mate_id = Some("bnd_partner_1".into());

        assert_eq!(v.variant_class.so_term(), "chromosome_breakpoint");
        assert!(v.variant_class.is_structural());
        assert_eq!(v.mate_id.as_deref(), Some("bnd_partner_1"));
    }

    fn make_me_variant(alt: &[u8]) -> InputVariant {
        let mut v = InputVariant::new("1".into(), 100, 200, b"N".to_vec(), alt.to_vec());
        v.variant_class = VariantClass::MobileElementInsertion;
        v.is_structural = true;
        v
    }

    #[test]
    fn test_display_allele_me_alu_insertion() {
        let v = make_me_variant(b"<INS:ME:ALU>");
        assert_eq!(v.display_allele(), "Alu_insertion");
    }

    #[test]
    fn test_display_allele_me_line1_insertion() {
        let v = make_me_variant(b"<INS:ME:LINE1>");
        assert_eq!(v.display_allele(), "LINE1_insertion");
    }

    #[test]
    fn test_display_allele_me_l1_insertion() {
        let v = make_me_variant(b"<INS:ME:L1>");
        assert_eq!(v.display_allele(), "LINE1_insertion");
    }

    #[test]
    fn test_display_allele_me_line_insertion() {
        // Ensembl VEP only maps L1→LINE1, not LINE→LINE1. <INS:ME:LINE> falls
        // through to generic mobile_element_insertion.
        let v = make_me_variant(b"<INS:ME:LINE>");
        assert_eq!(v.display_allele(), "mobile_element_insertion");
    }

    #[test]
    fn test_display_allele_me_sva_insertion() {
        let v = make_me_variant(b"<INS:ME:SVA>");
        assert_eq!(v.display_allele(), "SVA_insertion");
    }

    #[test]
    fn test_display_allele_me_herv_insertion() {
        let v = make_me_variant(b"<INS:ME:HERV>");
        assert_eq!(v.display_allele(), "HERV_insertion");
    }

    #[test]
    fn test_display_allele_me_generic_no_subtype() {
        let v = make_me_variant(b"<INS:ME>");
        assert_eq!(v.display_allele(), "mobile_element_insertion");
    }

    #[test]
    fn test_display_allele_me_unknown_subtype() {
        let v = make_me_variant(b"<INS:ME:UNKNOWN>");
        assert_eq!(v.display_allele(), "mobile_element_insertion");
    }

    #[test]
    fn test_display_allele_me_del_alu() {
        let mut v = make_me_variant(b"<DEL:ME:ALU>");
        v.variant_class = VariantClass::MobileElementDeletion;
        assert_eq!(v.display_allele(), "Alu_deletion");
    }

    #[test]
    fn test_display_allele_me_del_generic() {
        let mut v = make_me_variant(b"<DEL:ME>");
        v.variant_class = VariantClass::MobileElementDeletion;
        assert_eq!(v.display_allele(), "mobile_element_deletion");
    }

    #[test]
    fn test_display_allele_me_del_line() {
        let mut v = make_me_variant(b"<DEL:ME:LINE>");
        v.variant_class = VariantClass::MobileElementDeletion;
        assert_eq!(v.display_allele(), "mobile_element_deletion");
    }
}
