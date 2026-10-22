// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variant types: input variants and co-located known variants.

use std::borrow::Cow;
use std::sync::Arc;

use crate::codon::reverse_complement;
use crate::consequence::{Consequence, TranscriptConsequence};
use crate::coordinate::{normalize_chromosome, Strand};
use crate::sv_routing::{perl_true, sv_so_term_and_class, sv_type_abbreviation};

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
    /// SO terms match Perl VEP's `%SO_TERMS` mapping exactly.
    pub fn so_term(self) -> &'static str {
        match self {
            VariantClass::Snv => "SNV",
            VariantClass::Insertion => "insertion",
            VariantClass::Deletion => "deletion",
            VariantClass::Indel => "indel",
            VariantClass::Substitution => "substitution",
            VariantClass::SequenceAlteration => "sequence_alteration",
            // Structural variant terms follow Perl VEP `%SO_TERMS`.
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
    /// Perl VEP uses POS+1 as the start coordinate for all SV classes:
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
    /// Original chromosome name from input VCF (e.g., "chr21", "chrX", "chrM").
    /// Used in output to preserve the input naming convention, matching Perl VEP behavior.
    pub original_chr: String,
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
    /// for the record's own breakend and 1 upward for its mates (`:89-101`).
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
    /// Most severe consequence across all transcripts.
    pub most_severe_consequence: Option<Consequence>,
    /// Co-located known variants (from variation cache).
    pub colocated_variants: Vec<ColocatedVariant>,
    /// Whether the variant was minimised from a multi-allelic site.
    pub minimised: bool,
    /// Original allele string before minimisation.
    pub original_allele_string: Option<String>,
    /// Original start coordinate before per-allele secondary minimisation.
    /// Used for output parity with Perl VEP, which reverts multi-allelic
    /// records to their shared window coordinates after consequence calculation.
    pub original_start: Option<u64>,
    /// Original end coordinate before per-allele secondary minimisation.
    pub original_end: Option<u64>,
    /// Raw input line (for passthrough in some output formats), one buffer shared by
    /// every allele of the record.
    pub raw_input: Option<Arc<str>>,
    /// Existing variation IDs (from co-located variants).
    pub existing_variation: Vec<String>,
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
    /// For a tandem repeat (`<CNV:TR>`): the alternate allele's length in bases, RUC
    /// times the repeat unit's length when the record carries INFO/RUC, else INFO/RB
    /// rounded down to whole units. `None` when the record states no length or any
    /// repeat sequence of the record is longer than `--max_sv_size`: the record then
    /// describes no direction of change and is read symbolically.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tr_alt_bases: Option<u64>,
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
    /// A structural variant of a type VEP's `%SO_TERMS` table lacks (`<CPX>`,
    /// `<NON_REF>`). VEP carries its VCF line without CSQ and emits only its
    /// `input` in JSON; the default and tab formats still list the rows.
    #[serde(default)]
    pub vep_skip: bool,
    /// A structural variant wider than `--max_sv_size`. VEP carries its VCF
    /// line without CSQ and drops it from the JSON output altogether (`--json`
    /// sets the REST mode in which the parser skips it).
    #[serde(default)]
    pub oversize_sv: bool,
    /// The record of a line whose every ALT is `*`, the VCF marker for a base an
    /// upstream deletion covers: Ensembl VEP builds the feature over the reference
    /// bases with the allele string `N/*` and no allele to annotate, writes no row in
    /// the default and tab formats, carries the VCF line without CSQ and writes the
    /// JSON object with `most_severe_consequence` `?`. Nothing annotates it.
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
    /// mate INFO/CHR2 names (`N[21:100[`). One buffer shared by every allele of the
    /// record. `None` falls back to `allele_string`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uploaded_allele_string: Option<Arc<str>>,
    /// The record's allele string as VEP works with it, when it is not the allele's
    /// own: for a multi-allelic record every allele with the shared first base chopped
    /// (`TT/T/-`), for a breakpoint record the reference base joined to its ALT side
    /// (`N/N[21:100[/N]21:200]`) unless that side opens or closes with `.`, a native
    /// single breakend VEP keeps bare (Parser/VCF.pm:469 and 499). VEP's JSON
    /// `allele_string` prints this while Uploaded_variation prints the raw one. One
    /// buffer shared by every allele of the record.
    #[serde(default)]
    pub record_allele_string_multi: Option<Arc<str>>,
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
    /// the predicates follow it rather than the allele's own form. `None` for a
    /// sequence variant and for a structural type the table lacks.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub record_so_term: Option<String>,
}

/// An empty placeholder: no allocation, every field at its zero value. Lets a
/// worker take a finished variant out of a batch (`std::mem::take`) and drop it
/// where it was annotated, instead of leaving the drop to another thread.
impl Default for InputVariant {
    fn default() -> Self {
        Self {
            chr: String::new(),
            original_chr: String::new(),
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
            most_severe_consequence: None,
            colocated_variants: Vec::new(),
            minimised: false,
            original_allele_string: None,
            original_start: None,
            original_end: None,
            raw_input: None,
            existing_variation: Vec::new(),
            plugin_data: indexmap::IndexMap::new(),
            sv_end: None,
            sv_type: None,
            sv_len: None,
            tr_alt_bases: None,
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
            input_line_name: None,
            record_so_term: None,
        }
    }
}

impl InputVariant {
    /// Create a new input variant.
    ///
    /// The chromosome name is normalized at construction time (stripping `chr`
    /// prefix and converting `M` to `MT`) so that downstream matching against
    /// the Ensembl cache works regardless of the input naming convention.
    /// The original chromosome name is preserved for output (matching Perl VEP
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
            most_severe_consequence: None,
            colocated_variants: Vec::new(),
            minimised: false,
            original_allele_string: None,
            original_start: None,
            original_end: None,
            raw_input: None,
            existing_variation: Vec::new(),
            plugin_data: indexmap::IndexMap::new(),
            sv_end: None,
            sv_type: None,
            sv_len: None,
            tr_alt_bases: None,
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
            annotation_chr: None,
            annotation_start: None,
            annotation_end: None,
            input_record: 0,
            uploaded_allele_string: None,
            input_line_name: None,
            record_so_term: None,
        }
    }

    /// The record of a line whose every ALT is `*` (`N *`, `N *,*`): the reference
    /// bases from `pos`, the allele string REF and every ALT joined by `/`, and
    /// [`Self::no_allele`] set, so nothing annotates it and every writer carries it as
    /// Ensembl VEP carries its `N/*` feature.
    pub fn without_allele<'a>(
        chr: String,
        pos: u64,
        ref_allele: &'a str,
        alts: impl IntoIterator<Item = &'a str>,
    ) -> Self {
        let end = pos + ref_allele.len().saturating_sub(1) as u64;
        let mut variant = Self::new(chr, pos, end, ref_allele.as_bytes().to_vec(), b"*".to_vec());
        variant.allele_string = std::iter::once(ref_allele)
            .chain(alts)
            .collect::<Vec<_>>()
            .join("/");
        variant.no_allele = true;
        variant
    }

    /// Location string in VEP format ("chr:start" or "chr:start-end").
    /// Uses the original chromosome name from the input VCF to preserve
    /// the naming convention (e.g., "chr21:100-200" for chr-prefixed input).
    pub fn location(&self) -> String {
        // Perl VEP's `rejoin_variants_in_InputBuffer` restores a multi-allelic
        // record's shared window coordinates for output.
        let s = self.original_start.unwrap_or(self.start);
        let e = self.original_end.unwrap_or(self.end);
        let lo = s.min(e);
        let hi = s.max(e);
        if lo == hi {
            format!("{}:{}", self.original_chr, lo)
        } else {
            format!("{}:{}-{}", self.original_chr, lo, hi)
        }
    }

    /// The identifier that names the record in `Uploaded_variation`: [`Self::id`]
    /// when [`id_names_record`] admits it.
    fn explicit_id(&self) -> Option<&str> {
        self.id.as_deref().filter(|id| id_names_record(id))
    }

    /// `Uploaded_variation` as VEP prints it (OutputFactory.pm
    /// `VariationFeature_to_output_hash`): the identifier `explicit_id` admits unless
    /// it is the `chr_start_alleles` name VEP gives an ID-less record; else the name built
    /// from the input line ([`Self::input_line_name`]); else, for an input format
    /// VEP names from the variant, `chr_start_alleles` where start is the
    /// anchor-trimmed start and the alleles are [`Self::uploaded_allele_string`].
    pub fn uploaded_variation(&self) -> String {
        if let Some(id) = self.explicit_id().filter(|id| !self.is_generated_name(id)) {
            return id.to_string();
        }
        if let Some(name) = &self.input_line_name {
            return name.to_string();
        }
        let start = self.original_start.unwrap_or(self.start);
        let alleles = self
            .uploaded_allele_string
            .as_deref()
            .or(self.original_allele_string.as_deref())
            .unwrap_or(&self.allele_string);
        format!("{}_{}_{}", self.original_chr, start, alleles)
    }

    /// Whether `id` is the `chr_start_alleles` name VEP gives a record whose ID
    /// column is empty (Parser.pm `validate_vf`). The prefix test keeps an
    /// identifier such as `rs123` from building the name on every row.
    fn is_generated_name(&self, id: &str) -> bool {
        if !id.starts_with(self.original_chr.as_str()) {
            return false;
        }
        id == format!(
            "{}_{}_{}",
            self.original_chr,
            self.start,
            self.record_allele_string()
        )
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

    /// VEP's `variation_name` after `validate_vf`, the `id` of its JSON output
    /// (OutputFactory/JSON.pm:236): [`Self::id`] when Perl reads it as true, so a
    /// `.` stays `.`; else `chr_start_alleles` built from
    /// [`Self::record_allele_string`], the name Parser.pm:462 `||=` gives a record
    /// whose ID is absent, empty or `0`.
    pub fn variation_name(&self) -> String {
        if let Some(id) = self.id.as_deref().filter(|id| perl_true(id)) {
            return id.to_string();
        }
        let start = self.original_start.unwrap_or(self.start);
        format!(
            "{}_{}_{}",
            self.original_chr,
            start,
            self.record_allele_string()
        )
    }

    /// The specific alternate allele for this variant record.
    pub fn alt_allele(&self) -> &[u8] {
        self.alt_alleles
            .first()
            .map(|a| a.as_slice())
            .unwrap_or(b"-")
    }

    /// Whether the record writes its alleles as sequence on strand -1.
    ///
    /// Ensembl-format and region lines may. Ensembl VEP keeps such alleles as written
    /// in the record (Allele column, `allele_string`, `strand`) and reverse-complements
    /// them wherever it reads them on the forward strand or a transcript's:
    /// `VariationFeatureOverlapAllele::feature_seq` (VariationFeatureOverlapAllele.pm:253,
    /// sequence alleles only), the plugins' allele matching
    /// (`get_matched_variant_alleles`, Utils/Sequence.pm:1182) and the VCF writer
    /// (`to_VCF_record`, VariationFeature.pm:2376).
    pub fn has_reverse_strand_alleles(&self) -> bool {
        self.strand == Strand::Reverse && !self.is_structural && !self.variant_class.is_structural()
    }

    /// The reference and this record's alternate allele on the forward strand
    /// ([`Self::has_reverse_strand_alleles`]).
    pub fn forward_strand_alleles(&self) -> (Cow<'_, [u8]>, Cow<'_, [u8]>) {
        if self.has_reverse_strand_alleles() {
            (
                Cow::Owned(reverse_complement(&self.ref_allele)),
                Cow::Owned(reverse_complement(self.alt_allele())),
            )
        } else {
            (
                Cow::Borrowed(&self.ref_allele),
                Cow::Borrowed(self.alt_allele()),
            )
        }
    }

    /// This record with its alleles on the forward strand, the orientation the
    /// consequence engine reads them in; `None` unless
    /// [`Self::has_reverse_strand_alleles`].
    pub fn to_forward_strand(&self) -> Option<InputVariant> {
        if !self.has_reverse_strand_alleles() {
            return None;
        }
        let mut forward = self.clone();
        forward.strand = Strand::Forward;
        forward.ref_allele = reverse_complement(&self.ref_allele);
        for alt in &mut forward.alt_alleles {
            *alt = reverse_complement(alt);
        }
        forward.allele_string = self
            .allele_string
            .split('/')
            .map(|a| String::from_utf8_lossy(&reverse_complement(a.as_bytes())).into_owned())
            .collect::<Vec<_>>()
            .join("/");
        Some(forward)
    }

    /// Allele string for the specific alt allele in display format.
    ///
    /// For a structural variant, the record's SO term ([`Self::record_so_term`], the
    /// `class_SO_term` every row of the record prints in Perl VEP), or for a
    /// variant built without one the term of its class, read from its own ALT for
    /// a mobile element. A bracket breakend prints its string as written
    /// (`A]21:33033339]`), the own breakend of a record that names a mate its `A.`
    /// or `.A`, and a translocation allele without a breakend string, a single
    /// breakend on its own or a bare `<BND>`, the class term `chromosome_breakpoint`:
    /// Perl VEP's `$vfoa->{breakend}->{string} || $svf->class_SO_term`
    /// (`OutputFactory.pm:2029`).
    /// For `<NON_REF>` alleles, returns the raw allele string since Perl VEP outputs
    /// the literal `<NON_REF>` (it has no SO term mapping).
    /// For small variants, returns the actual allele bases.
    pub fn display_allele(&self) -> String {
        if self.is_structural || self.variant_class.is_structural() {
            let raw = String::from_utf8_lossy(self.alt_allele()).to_string();
            if self.variant_class == VariantClass::Translocation {
                if self.mate_id.is_some() {
                    raw
                } else {
                    self.variant_class.so_term().to_string()
                }
            } else if let Some(term) = &self.record_so_term {
                term.clone()
            } else if self.variant_class == VariantClass::ReferenceBlock {
                // Perl VEP outputs the allele as written (no SO term mapping).
                raw
            } else if self.variant_class == VariantClass::ComplexStructural {
                // <CPX>: Perl VEP strips angle brackets to get "CPX" as class_SO_term,
                // then outputs that as the allele string.
                raw.trim_start_matches('<')
                    .trim_end_matches('>')
                    .split(':')
                    .next()
                    .unwrap_or("CPX")
                    .to_string()
            } else if self.variant_class == VariantClass::MobileElementInsertion
                || self.variant_class == VariantClass::MobileElementDeletion
            {
                // The ALT's own type read as the record's would be (`<INS:ME:ALU>` is
                // `Alu_insertion`, `<DEL:ME>` `mobile_element_deletion`).
                sv_so_term_and_class(&sv_type_abbreviation(&raw, None)).map_or_else(
                    || self.variant_class.so_term().to_string(),
                    |(term, _)| term.to_string(),
                )
            } else {
                self.variant_class.so_term().to_string()
            }
        } else if self.original_start.is_some() {
            // Perl VEP outputs the pre-minimised alt of a multi-allelic record;
            // `original_start` is set only for those records.
            self.original_allele_string
                .as_deref()
                .and_then(|s| s.split('/').nth(1))
                .unwrap_or("-")
                .to_string()
        } else {
            String::from_utf8_lossy(self.alt_allele()).to_string()
        }
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
    /// Population allele frequencies.
    pub frequencies: std::collections::HashMap<String, f64>,
    /// PubMed IDs.
    pub pubmed: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn test_sv_so_terms_match_perl_vep() {
        // SO terms must match Perl VEP %SO_TERMS exactly.
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
        // Perl VEP only maps L1→LINE1, not LINE→LINE1. <INS:ME:LINE> falls
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
