// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variant types: input variants and co-located known variants.

use crate::consequence::{Consequence, TranscriptConsequence};
use crate::coordinate::{normalize_chromosome, Strand};

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
    /// Original input identifier (e.g., VCF ID field, or "chr:pos:ref:alt").
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
    /// Raw input line (for passthrough in some output formats).
    pub raw_input: Option<String>,
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
    /// For a tandem repeat (`<CNV:TR>`): the alternate allele's length in bases, read from
    /// INFO/RB, else RUC times the repeat unit's length. `None` when the record carries
    /// neither, so the record describes no direction of change.
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
    /// For BND: MATEID from VCF INFO field.
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
    /// alleles of every other record with all ALTs joined by `/`, or the symbolic
    /// ALT alone for a structural variant. `None` falls back to `allele_string`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uploaded_allele_string: Option<String>,
    /// The record's allele string as VEP works with it, for a multi-allelic
    /// record: every allele with the shared first base chopped (`TT/T/-`). VEP's
    /// JSON `allele_string` prints this while Uploaded_variation prints the raw one.
    #[serde(default)]
    pub record_allele_string_multi: Option<String>,
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
            record_allele_string_multi: None,
            annotation_chr: None,
            annotation_start: None,
            annotation_end: None,
            input_record: 0,
            uploaded_allele_string: None,
        }
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

    /// The explicit input identifier, if the record carried one (VCF ID other
    /// than `.`; the first of several `;`-separated IDs).
    fn explicit_id(&self) -> Option<&str> {
        self.id
            .as_deref()
            .filter(|id| !id.is_empty() && *id != ".")
            .map(|id| id.split(';').next().unwrap_or(id))
    }

    /// `Uploaded_variation` as VEP prints it: the input ID when there is one, else
    /// `chr_start_alleles` where start is the anchor-trimmed start and the alleles
    /// are [`Self::uploaded_allele_string`] (raw REF/ALT for a minimised bi-allelic
    /// record, the record's full allele string otherwise), joined by `/`.
    pub fn uploaded_variation(&self) -> String {
        if let Some(id) = self.explicit_id() {
            return id.to_string();
        }
        let start = self.original_start.unwrap_or(self.start);
        let alleles = self
            .uploaded_allele_string
            .as_deref()
            .or(self.original_allele_string.as_deref())
            .unwrap_or(&self.allele_string);
        format!("{}_{}_{}", self.original_chr, start, alleles)
    }

    /// The allele string of the whole input record as VEP holds it: the symbolic
    /// ALT(s) of a structural variant, every anchor-trimmed allele of a
    /// multi-allelic record joined by `/`, or the minimised `allele_string` of a
    /// bi-allelic record.
    pub fn record_allele_string(&self) -> String {
        // A breakend record's string is REF joined to the paired form
        // (`N/N[21:100[`), except a native single breakend (`N.`), which VEP
        // keeps bare (ensembl-vep `Parser/VCF.pm` `create_StructuralVariationFeatures`).
        if self.variant_class == VariantClass::Translocation && !self.is_single_breakend {
            return self.allele_string.clone();
        }
        if let Some(multi) = &self.record_allele_string_multi {
            return multi.clone();
        }
        match self.uploaded_allele_string.as_deref() {
            Some(s) if self.is_structural => s.to_string(),
            _ => self.allele_string.clone(),
        }
    }

    /// VEP's `variation_name`, the `id` of its JSON output: the first input ID
    /// when there is one; the literal `.` for a VCF record whose ID column is
    /// `.` (VEP keeps the dot as the name); else `chr_start_alleles` built from
    /// [`Self::record_allele_string`].
    pub fn variation_name(&self) -> String {
        if let Some(id) = self.explicit_id() {
            return id.to_string();
        }
        if self.uploaded_allele_string.is_some() {
            return ".".to_string();
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

    /// Allele string for the specific alt allele in display format.
    ///
    /// For most structural variants, returns the SO term (e.g., "deletion", "duplication")
    /// matching Perl VEP behavior. For BND/translocation variants, returns the raw allele
    /// string (e.g., "A]21:33033339]") since Perl VEP outputs the literal breakend notation.
    /// For `<NON_REF>` alleles, returns the raw allele string since Perl VEP outputs
    /// the literal `<NON_REF>` (it has no SO term mapping).
    /// For small variants, returns the actual allele bases.
    pub fn display_allele(&self) -> String {
        if self.is_structural || self.variant_class.is_structural() {
            let raw = String::from_utf8_lossy(self.alt_allele()).to_string();
            if self.variant_class == VariantClass::Translocation {
                if self.is_single_breakend && self.mate_id.is_none() {
                    self.variant_class.so_term().to_string()
                } else {
                    // Paired BNDs keep the literal allele string, as Perl VEP does.
                    raw
                }
            } else if raw.eq_ignore_ascii_case("<NON_REF>") {
                // <NON_REF>: Perl VEP outputs the raw allele (no SO term mapping)
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
                // Perl VEP outputs ME-subtype-specific alleles from the ALT field:
                //   <INS:ME:ALU> → Alu_insertion, <DEL:ME:SVA> → SVA_deletion, etc.
                me_display_allele(self.alt_allele())
                    .unwrap_or_else(|| self.variant_class.so_term().to_string())
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

/// Extract a Perl VEP-style display allele for mobile element variants.
///
/// Perl VEP parses the ME subtype from the symbolic ALT allele and produces
/// subtype-specific names:
///   `<INS:ME:ALU>` → `Alu_insertion`
///   `<DEL:ME:SVA>` → `SVA_deletion`
///   `<INS:ME:LINE1>` → `LINE1_insertion`
///   `<INS:ME>` → `mobile_element_insertion`
fn me_display_allele(raw_alt: &[u8]) -> Option<String> {
    let s = std::str::from_utf8(raw_alt).ok()?;
    let inner = s.strip_prefix('<')?.strip_suffix('>')?;
    let parts: Vec<&str> = inner.split(':').collect();
    if parts.len() < 2 {
        return None;
    }
    let operation = if parts[0].eq_ignore_ascii_case("DEL") {
        "deletion"
    } else {
        "insertion"
    };

    if parts.len() >= 3 {
        let element = parts[2].to_ascii_uppercase();
        let canonical = match element.as_str() {
            "ALU" => Some("Alu"),
            "L1" | "LINE1" => Some("LINE1"),
            "SVA" => Some("SVA"),
            "HERV" => Some("HERV"),
            _ => None,
        };
        if let Some(name) = canonical {
            return Some(format!("{}_{}", name, operation));
        }
    }
    Some(format!("mobile_element_{}", operation))
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

        v.id = Some("rs123;rs456".into());
        assert_eq!(v.uploaded_variation(), "rs123");
        assert_eq!(v.variation_name(), "rs123");

        v.id = Some(".".into());
        assert_eq!(v.uploaded_variation(), "21_100_A/G");
    }

    #[test]
    fn uploaded_variation_names_a_minimised_indel_by_its_raw_alleles() {
        // `21 100 . AT A`: VEP minimises to T/- at 101 and names the record by the raw
        // alleles at the trimmed start; its JSON id keeps the VCF's `.`.
        let mut v = InputVariant::new("21".into(), 101, 101, b"T".to_vec(), b"-".to_vec());
        v.uploaded_allele_string = Some("AT/A".into());
        v.minimised = true;
        assert_eq!(v.uploaded_variation(), "21_101_AT/A");
        assert_eq!(v.variation_name(), ".");
        assert_eq!(v.record_allele_string(), "T/-");
        // A record from a format without an ID column is named by its alleles.
        v.uploaded_allele_string = None;
        assert_eq!(v.variation_name(), "21_101_T/-");
    }

    #[test]
    fn multi_allelic_and_structural_records_share_one_name() {
        // `21 100 . TAG TA,T`: named by the raw alleles at the chopped start; the
        // JSON allele_string is the chopped working string.
        let mut v = InputVariant::new("21".into(), 101, 101, b"AG".to_vec(), b"A".to_vec());
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
        ];
        for vc in sv_types {
            assert!(vc.is_structural(), "{vc:?} should be structural");
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

    #[test]
    fn test_display_allele_paired_breakend_variants_keep_raw_allele() {
        let mut paired =
            InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"A[21:200[".to_vec());
        paired.variant_class = VariantClass::Translocation;
        paired.is_structural = true;
        assert_eq!(paired.display_allele(), "A[21:200[");

        let mut derived_single =
            InputVariant::new("1".into(), 100, 100, b"N".to_vec(), b"A.".to_vec());
        derived_single.variant_class = VariantClass::Translocation;
        derived_single.is_structural = true;
        derived_single.is_single_breakend = true;
        derived_single.mate_id = Some("bnd_2".into());
        assert_eq!(derived_single.display_allele(), "A.");
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
