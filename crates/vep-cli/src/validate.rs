// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's record validation, `validate_vf` (Parser.pm:457-649 at
//! release/116.2), over the alleles vep-rs parsed from one input record.
//!
//! Ensembl VEP builds one feature per input line and runs the checks on it in a
//! fixed order, stopping at the first that fails: the `--chr` list, the
//! coordinates, the `--max_sv_size` bound, the chromosome against the cache, the
//! allele string (upper-cased first), an insertion's coordinates, the reference
//! allele's length against the span, then `--check_ref` or `--lookup_ref`
//! against the reference sequence. A failing record is reported through the
//! warning sink and dropped, or kept as it stands under `--dont_skip`
//! (Parser.pm:189, Parser/VCF.pm:199), the checks after the failing one never
//! run on it. vep-rs holds a record as one [`InputVariant`] per alternate
//! allele, so the verdict is taken on the first allele and what the checks
//! change (the upper-casing, a looked-up reference, the failed check mark) is
//! carried to the record's other alleles.

use std::borrow::Cow;
use std::collections::HashSet;

use vep_core::codon::reverse_complement;
use vep_core::coordinate::{Strand, MAX_COORD};
use vep_core::sv_routing::{perl_split, perl_true, sv_type_abbreviation, sv_type_string};
use vep_core::sv_span::deletion_looks_incomplete;
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{classify_variant, InputVariant};
use vep_core::variation::{perl_number_string, perl_numeric};
use vep_fasta::IndexedFasta;

use crate::pipeline::{Reported, SkipSink};
use crate::transcript_index::LazyTranscriptIndexes;
use crate::vcf_parser::RecordInfo;

/// The reference-sequence checks of a run and the `--chr` list.
pub(crate) struct ReferenceChecks<'a> {
    /// The chromosome names `--chr` admits; `None` without the flag.
    pub chr: Option<&'a ChromosomeList>,
    /// `--check_ref`: a reference allele that differs from the reference
    /// sequence fails the record.
    pub check_ref: bool,
    /// `--lookup_ref`: the reference allele is replaced by the reference
    /// sequence over the record's span; read only when `check_ref` is off, as
    /// `validate_vf`'s `elsif` reads it (Parser.pm:635).
    pub lookup_ref: bool,
    /// The run's reference FASTA, which both checks read.
    pub fasta: Option<&'a IndexedFasta>,
    /// The input format's parser names every record itself (Parser/Region.pm:198,
    /// Parser/HGVS.pm:237, Parser/SPDI.pm:193), so Parser.pm:462 builds no name
    /// from the alleles.
    pub parser_names_records: bool,
}

impl ReferenceChecks<'_> {
    /// Whether the `--chr` list admits a chromosome as the input spells it; every
    /// name without the flag.
    pub(crate) fn admits(&self, name: &str) -> bool {
        self.chr.is_none_or(|admitted| admitted.admits(name))
    }
}

/// No list and neither check: every record passes these three.
#[cfg(test)]
pub(crate) const NO_CHECKS: ReferenceChecks<'static> = ReferenceChecks {
    chr: None,
    check_ref: false,
    lookup_ref: false,
    fasta: None,
    parser_names_records: false,
};

/// The chromosome names a `--chr` value admits: each comma-separated token as
/// written (Ensembl VEP matches the raw input name against the token,
/// Parser.pm:470-472), and a token of two non-negative integers joined by `-`
/// as every integer between them, so `1-5,10,MT` admits `1` to `5`. Ensembl
/// VEP admits the token `1-5` alone.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(crate) struct ChromosomeList {
    names: HashSet<String>,
    ranges: Vec<(u64, u64)>,
}

impl ChromosomeList {
    pub(crate) fn new(tokens: &[String]) -> Self {
        // `str::parse::<u64>` takes a leading `+`, so a bound is tested as a run
        // of decimal digits first and `+1-3` stays a name.
        let bound = |text: &str| -> Option<u64> {
            if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
                return None;
            }
            text.parse().ok()
        };
        let mut list = ChromosomeList::default();
        for token in tokens.iter().filter(|t| !t.is_empty()) {
            let bounds = token
                .split_once('-')
                .and_then(|(from, to)| Some((bound(from)?, bound(to)?)));
            match bounds {
                Some((from, to)) if from <= to => list.ranges.push((from, to)),
                _ => {
                    list.names.insert(token.clone());
                }
            }
        }
        list
    }

    /// Whether a record whose chromosome the input spells `name` is annotated: a
    /// name as written, or the plain decimal spelling of an integer in a range.
    pub(crate) fn admits(&self, name: &str) -> bool {
        self.names.contains(name)
            || name.parse::<u64>().is_ok_and(|n| {
                n.to_string() == name
                    && self
                        .ranges
                        .iter()
                        .any(|&(from, to)| (from..=to).contains(&n))
            })
    }
}

/// Perl's `looks_like_number`, the test of `validate_vf`'s coordinate check
/// (Parser.pm:475): a decimal number with an optional sign, fraction and
/// exponent, or `Inf`, `Infinity` or `NaN` in any case, blanks around it allowed.
pub(crate) fn looks_like_number(text: &str) -> bool {
    let text = text.trim();
    let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
    if ["inf", "infinity", "nan"].contains(&unsigned.to_ascii_lowercase().as_str()) {
        return true;
    }
    let (mantissa, exponent) = match unsigned.split_once(['e', 'E']) {
        Some((m, e)) => (m, Some(e)),
        None => (unsigned, None),
    };
    let (whole, fraction) = match mantissa.split_once('.') {
        Some((w, f)) => (w, Some(f)),
        None => (mantissa, None),
    };
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let mantissa_ok = match fraction {
        Some(f) => {
            (digits(whole) && (f.is_empty() || digits(f))) || (whole.is_empty() && digits(f))
        }
        None => digits(whole),
    };
    let exponent_ok = exponent.is_none_or(|e| digits(e.strip_prefix(['+', '-']).unwrap_or(e)));
    mantissa_ok && exponent_ok
}

/// The skip description of `validate_vf`'s coordinate test (Parser.pm:476-480)
/// for a VCF record whose POS is not a number: Ensembl VEP's end is the
/// numified POS plus the reference allele's length less one (Parser/VCF.pm:254),
/// `abc` numifying to 0.
pub(crate) fn invalid_vcf_coordinate(pos: &str, ref_allele: &str) -> String {
    let end = perl_numeric(pos) + ref_allele.len() as f64 - 1.0;
    format!(
        "Invalid start '{pos}' or end '{}' coordinate",
        perl_number_string(end)
    )
}

/// The reference sequence over `start..=end` of `chr`, on `strand`, as
/// `_get_ref_allele` reads it from the FASTA-backed slice (Parser.pm:745-753:
/// `sub_Slice($start, $end, $strand)->seq`, reverse-complemented for strand
/// -1); `None` when the FASTA has no such span.
pub(crate) fn reference_sequence(
    fasta: &IndexedFasta,
    chr: &str,
    start: u64,
    end: u64,
    strand: Strand,
) -> Option<Vec<u8>> {
    let sequence = fasta.sequence(chr, start, end)?;
    Some(match strand {
        Strand::Forward => sequence,
        Strand::Reverse => reverse_complement(&sequence),
    })
}

/// The checks Ensembl VEP's parser makes on a record before annotation, applied
/// to the alleles vep-rs parsed from it, and what a record that passes them
/// goes through next.
pub(crate) struct SkipRules<'a> {
    pub transcripts: &'a LazyTranscriptIndexes,
    /// The chromosome synonyms an input name is resolved through before the
    /// checks run (`pipeline::name_source_chromosomes`).
    pub synonyms: &'a ChromosomeSynonyms,
    pub max_sv_size: Option<u64>,
    pub dont_skip: bool,
    pub checks: &'a ReferenceChecks<'a>,
    /// `--minimal`: a validated multi-allelic record is split per allele
    /// (`vcf_parser::split_minimal_record`), after the checks have read it.
    pub minimal: bool,
}

/// The input record as the skip rules read it: `number` Ensembl VEP's count of
/// the records read, `info` the record's INFO column as the VCF reader parsed
/// it (`None` for a format with no INFO column) and `quoted` the bytes a skip
/// warning quotes.
pub(crate) struct RecordLine<'a, 'i> {
    pub number: u64,
    pub info: Option<&'a RecordInfo<'i>>,
    pub quoted: &'a [u8],
}

/// What the checks changed on a record's first allele, to carry to the alleles
/// vep-rs split from the same line.
pub(crate) struct AllelePatch {
    /// The name given to a record whose ID column names it not, when the checks
    /// change the alleles it is built from.
    name: Option<String>,
    upper_cased: bool,
    /// The record was refused before the upper-casing and kept by `--dont_skip`,
    /// so its alleles stand as written.
    written_case: bool,
    /// The reference allele after `--lookup_ref`, when it replaced the record's.
    reference: Option<Vec<u8>>,
    check_ref_failed: bool,
}

impl AllelePatch {
    /// Applies the record's changes to one of its other alleles.
    pub(crate) fn apply(&self, variant: &mut InputVariant) {
        if let Some(name) = &self.name {
            variant.id = Some(name.clone());
        }
        if self.upper_cased {
            upper_case_alleles(variant);
        }
        if self.written_case {
            restore_written_case(variant);
        }
        if let Some(reference) = &self.reference {
            replace_reference(variant, reference);
        }
        variant.check_ref_failed = self.check_ref_failed;
    }

    /// Whether `--lookup_ref` replaced the record's reference, so the overlap-allele
    /// rule runs again over its alleles ([`keep_overlap_allele`]).
    pub(crate) fn replaced_reference(&self) -> bool {
        self.reference.is_some()
    }
}

/// Whether an alternate allele is one Ensembl VEP numbers and writes no row for:
/// a `*` (the VCF marker for a base an upstream deletion covers) or a `<DEL:*>`
/// (VariationFeatureOverlap.pm:483-484). The VCF parser leaves such an allele out
/// of a record's alleles; the parsers of the other formats keep it.
pub(crate) fn rowless_allele(variant: &InputVariant) -> bool {
    matches!(variant.alt_allele(), b"*" | b"<DEL:*>")
}

/// The overlap-allele rule on a record whose alleles may equal its reference or be
/// a `*` or `<DEL:*>`: every record of the parser-fed formats, whose parsers keep
/// such an allele, and a VCF record once `--lookup_ref` replaced its reference.
/// Ensembl VEP builds the overlap alleles at annotation time
/// (VariationFeatureOverlap.pm:459-491), so an ALT equal to the reference writes no
/// row and takes no number (`:481`), and the alleles after it are numbered without
/// it (`:483`); a `*` or `<DEL:*>` takes its number and writes no row (`:484`), so
/// the alleles after it keep theirs. The one ALT of a bi-allelic record, or the
/// last allele left of a record, then makes it the feature with no allele to
/// annotate, written as a non-variant line is. The reference row of
/// `--process_ref_homs` stands whatever the alternates equal: Ensembl VEP builds the
/// reference overlap allele, numbered 0, before it tests any alternate (`:469-474`),
/// so a record whose every alternate equals the looked-up reference keeps that row
/// alone. Returns whether the allele still writes a row; `last` says it is the
/// record's last allele still standing.
pub(crate) fn keep_overlap_allele(variant: &mut InputVariant, last: bool) -> bool {
    if variant.is_structural || variant.no_allele || variant.reference_allele_row {
        return true;
    }
    let rowless = rowless_allele(variant);
    let equal_to_reference = match variant.record_allele_string_multi.as_deref() {
        None => variant.alt_allele() == variant.ref_allele.as_slice(),
        Some(multi) => {
            let reference = String::from_utf8_lossy(&variant.ref_allele);
            let alts: Vec<String> = multi.split('/').skip(1).map(str::to_string).collect();
            let numbers = crate::vcf_parser::perl_allele_numbers(&reference, &alts);
            match numbers.get(variant.allele_index).copied().flatten() {
                Some(number) => {
                    variant.allele_num =
                        (number as usize != variant.allele_index + 1).then_some(number);
                    false
                }
                None => true,
            }
        }
    };
    if !equal_to_reference && !rowless {
        return true;
    }
    if last {
        variant.no_allele = true;
        variant.alt_alleles = vec![b"*".to_vec()];
        return true;
    }
    false
}

impl SkipRules<'_> {
    /// Reports every skip Ensembl VEP reports for the record, applies what the
    /// checks change to `first`, the record's first allele, and returns whether
    /// vep-rs annotates the record and the changes its other alleles take. The
    /// parse-time marks come first: a type outside `%SO_TERMS` (Parser/VCF.pm:477-481),
    /// warned and then, in vep-rs, still annotated, and a deletion without an extent
    /// (`:555-558`), warned and dropped where `validate_svf` drops it
    /// (Parser.pm:821-826). Then `validate_vf` (Parser.pm:457-649) in its order,
    /// each failing check dropping the record unless `--dont_skip` keeps it as it
    /// stands (Parser/VCF.pm:199-201): the `--chr` list, read against the input's
    /// own name and silent (`:470-472`); a start past the end plus one (`:483-488`);
    /// a span over `--max_sv_size` (`:493-500`, warned and in vep-rs still
    /// annotated, unless the span reached the coordinate ceiling, where its SVLEN
    /// or END overflowed and the span is not the record's: that record is
    /// dropped unless `--dont_skip` keeps it); a chromosome the cache does not
    /// carry (`:542-546`); then, for a sequence variant, the allele string
    /// upper-cased (`:554`), an allele string with no `ACGTN-` run (`:556-560`), an
    /// insertion whose start is not the end plus one (`:564-572`), a reference
    /// allele whose length is not the span's (`:585-592`), and the reference allele
    /// against the reference sequence (`--check_ref`, `:596-633`) or replaced by it
    /// (`--lookup_ref`, `:635-646`).
    pub(crate) fn validate(
        &self,
        first: &mut InputVariant,
        line: &RecordLine<'_, '_>,
        report: &mut SkipSink<'_>,
    ) -> (bool, AllelePatch) {
        let mut patch = AllelePatch {
            name: None,
            upper_cased: false,
            written_case: false,
            reference: None,
            check_ref_failed: false,
        };
        let quoted = String::from_utf8_lossy(line.quoted);
        let mut skip = |description: &str| {
            report(Reported::Skip {
                description,
                line_number: line.number,
                quoted: Some(&quoted),
            })
        };
        let mut incomplete = false;
        if first.is_structural {
            let joined_alts = first.uploaded_allele_string.as_deref().unwrap_or("");
            let svtype = first.sv_type.as_deref();
            // `get_SO_term` warns about a type `%SO_TERMS` lacks by its abbreviation
            // (Parser.pm:668-725) whether the record is then carried without
            // consequences (a VCF record, `vep_skip`) or annotated under its word
            // (an Ensembl-format record).
            if first.record_so_term.is_none() {
                let abbrev = sv_type_abbreviation(joined_alts, svtype);
                skip(&format!(
                    "{abbrev} is not a supported structural variant type"
                ));
            }
            // Parser/VCF.pm:477-481: the deletion test runs on the SO term, or on
            // the raw type string when `%SO_TERMS` has no entry for it.
            let so_term = match &first.record_so_term {
                Some(term) => term.clone(),
                None => sv_type_string(joined_alts, svtype),
            };
            incomplete = incomplete_deletion(first, &so_term, line.info);
            if incomplete {
                skip("deletion looks incomplete");
            }
        }
        // The name as the input spelled it, before any rename (Parser.pm:470 runs
        // before `_have_chr`).
        let input_name = first
            .input_chr
            .as_deref()
            .unwrap_or(first.original_chr.as_str());
        if !self.checks.admits(input_name) {
            return self.keep_as_written(first, patch);
        }
        // A start past the end plus one names a span no annotation can read, so
        // the record is dropped under `--dont_skip` too, where Ensembl VEP's run
        // ends on the first transcript the record reaches.
        if first.start > first.end + 1 {
            skip(&format!(
                "start > end+1 : (START={}, END={})",
                first.start, first.end
            ));
            return (false, patch);
        }
        if let (true, Some(max)) = (first.oversize_sv, self.max_sv_size) {
            skip(&format!(
                "variant size ({}) is bigger than --max_sv_size ({max})",
                first.end.saturating_sub(first.start)
            ));
            // A span that saturated at the ceiling runs to the end of every
            // chromosome; the rows it would write belong to no record. `--dont_skip`
            // keeps it as it keeps every refused record.
            return (first.end < MAX_COORD || self.dont_skip, patch);
        }
        if !self.transcripts.contains_key(&first.chr) {
            let chr = &first.original_chr;
            skip(&format!(
                "Chromosome {chr} not found in annotation sources or synonyms; chromosome {chr} does not overlap any features"
            ));
            return self.keep_as_written(first, patch);
        }
        if first.is_structural {
            // `validate_svf` (Parser.pm:821-826) drops a deletion marked incomplete
            // unless `--dont_skip` keeps it; a type `%SO_TERMS` lacks carries the same
            // mark and is still annotated here.
            return (!incomplete || self.dont_skip, patch);
        }
        // Parser.pm:462 names a record whose ID column names it not from its
        // alleles as written, before the upper-casing (:554) and the looked-up
        // reference (:645) change them; the name is kept where those steps would
        // change what the record builds at output.
        let named_by_alleles =
            !self.checks.parser_names_records && !first.id.as_deref().is_some_and(perl_true);
        // The string as written covers every allele of the record, a parser's
        // working copy of a multi-allelic string being upper-cased already; an
        // individual's copy of a record carries its own alleles as written in
        // `written_case_allele_string` alone, its uploaded string being the
        // upper-cased copy (Parser/VCF.pm:639 builds the copy's string from the
        // bits as written, which `:462` then names it by).
        let as_written = first.uploaded_allele_string.as_deref().unwrap_or("");
        let written_case = first.written_case_allele_string.as_deref().unwrap_or("");
        let lower_case = working_allele_string(first)
            .bytes()
            .chain(first.allele_string.bytes())
            .chain(as_written.bytes())
            .chain(written_case.bytes())
            .any(|b| b.is_ascii_lowercase());
        let checks_change_alleles =
            (self.checks.lookup_ref && !self.checks.check_ref) || lower_case;
        if named_by_alleles && checks_change_alleles {
            let name = first.variation_name();
            first.id = Some(name.clone());
            patch.name = Some(name);
        }
        if lower_case {
            upper_case_alleles(first);
            patch.upper_cased = true;
        }
        let allele_string = working_allele_string(first);
        if !allele_string
            .bytes()
            .any(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'-'))
        {
            skip(&format!(
                "Invalid allele string {allele_string} or possible parsing error"
            ));
            return (self.dont_skip, patch);
        }
        if allele_string.starts_with("-/") && first.start != first.end + 1 {
            skip(&format!(
                "Alleles look like an insertion ({allele_string}) but coordinates are not start = end + 1 (START={}, END={}) for variant ({}) ",
                first.start,
                first.end,
                first.variation_name()
            ));
            return (self.dont_skip, patch);
        }
        // Parser.pm:578-580 removes the `-` from every allele: the reference's
        // here, so an insertion's is empty, and the alternates' at the `--lookup_ref`
        // rejoin (`replace_reference`), the one place their stripped form prints.
        let reference: Cow<[u8]> = if first.ref_allele.contains(&b'-') {
            Cow::Owned(
                first
                    .ref_allele
                    .iter()
                    .copied()
                    .filter(|&b| b != b'-')
                    .collect(),
            )
        } else {
            Cow::Borrowed(&first.ref_allele)
        };
        let span = first.end + 1 - first.start;
        if reference
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            && span != reference.len() as u64
        {
            skip(&format!(
                "Length of reference allele ({} length {}) does not match coordinates {}-{}",
                String::from_utf8_lossy(&reference),
                reference.len(),
                first.start,
                first.end
            ));
            return (self.dont_skip, patch);
        }
        let fetch = || {
            self.checks.fasta.and_then(|fasta| {
                reference_sequence(fasta, &first.chr, first.start, first.end, first.strand)
            })
        };
        let span_text = || {
            format!(
                "{}:{}-{}({})",
                first.original_chr,
                first.start,
                first.end,
                first.strand.as_i8()
            )
        };
        if self.checks.check_ref {
            // An insertion has no reference allele to check (Parser.pm:601-603).
            let matched = if reference.is_empty() {
                true
            } else {
                match fetch() {
                    Some(sequence) => sequence.eq_ignore_ascii_case(&reference),
                    None => {
                        skip(&format!("Could not fetch sub-slice from {}", span_text()));
                        false
                    }
                }
            };
            if !matched {
                first.check_ref_failed = true;
                patch.check_ref_failed = true;
                // The expected base never prints: the `my` of Parser.pm:605 shadows
                // the variable `:627-629` reads.
                skip(&format!(
                    "Specified reference allele {} does not match Ensembl reference allele",
                    String::from_utf8_lossy(&reference)
                ));
                return (self.dont_skip, patch);
            }
        } else if self.checks.lookup_ref {
            let looked_up = if reference.is_empty() {
                Some(Vec::new())
            } else {
                fetch()
            };
            let Some(looked_up) = looked_up else {
                skip(&format!(
                    "Could not fetch reference allele for {}",
                    span_text()
                ));
                return (self.dont_skip, patch);
            };
            replace_reference(first, &looked_up);
            patch.reference = Some(looked_up);
        }
        (true, patch)
    }

    /// The verdict on a record the `--chr` list (Parser.pm:470-472) or the
    /// chromosome test (`:542-546`) refused: dropped, or kept by `--dont_skip`
    /// with its alleles as written, since both tests return before the
    /// upper-casing of `:554` and nothing after `validate_vf` changes the case.
    fn keep_as_written(
        &self,
        first: &mut InputVariant,
        mut patch: AllelePatch,
    ) -> (bool, AllelePatch) {
        if self.dont_skip {
            restore_written_case(first);
            patch.written_case = true;
        }
        (self.dont_skip, patch)
    }
}

/// The allele string Ensembl VEP tests after the anchor chop and the minimisation
/// (Parser.pm:554 runs on the post-processed feature): the record's chopped
/// multi-allelic string, else the allele's own.
fn working_allele_string(variant: &InputVariant) -> &str {
    variant
        .record_allele_string_multi
        .as_deref()
        .unwrap_or(&variant.allele_string)
}

/// `tr/[a-z]/[A-Z]/` on the allele string (Parser.pm:554), which every allele of
/// the record shares, so the alleles themselves follow.
fn upper_case_alleles(variant: &mut InputVariant) {
    variant.ref_allele.make_ascii_uppercase();
    for alt in &mut variant.alt_alleles {
        alt.make_ascii_uppercase();
    }
    variant.allele_string.make_ascii_uppercase();
    if let Some(multi) = &variant.record_allele_string_multi {
        if multi.bytes().any(|b| b.is_ascii_lowercase()) {
            variant.record_allele_string_multi = Some(multi.to_ascii_uppercase().into());
        }
    }
    reclassify(variant);
}

/// The alleles in the case the input wrote them, read from the record's
/// written-case string: the reference, the allele whose upper-cased form the
/// variant carries (the reference again on a reference row), the allele's own
/// string and the record's working string. A record whose string has no
/// lower-case base, or whose working string was rebuilt since the parse, is
/// left as it stands. The class reads the alleles' lengths alone and is
/// unchanged.
fn restore_written_case(variant: &mut InputVariant) {
    let working = variant.record_allele_string();
    let Some(written) = variant
        .written_case_allele_string
        .clone()
        .filter(|written| written.eq_ignore_ascii_case(&working))
    else {
        return;
    };
    let mut parts = written.split('/');
    let reference = parts.next().unwrap_or("").to_string();
    if variant.no_allele {
        variant.allele_string = written.to_string();
    } else {
        let alt = if variant.reference_allele_row {
            Some(reference.clone())
        } else {
            let carried = String::from_utf8_lossy(variant.alt_allele());
            parts
                .find(|part| part.eq_ignore_ascii_case(&carried))
                .map(str::to_string)
        };
        if let (Some(alt), Some(carried)) = (alt, variant.alt_alleles.first_mut()) {
            variant.allele_string = format!("{reference}/{alt}");
            *carried = alt.into_bytes();
        }
    }
    variant.ref_allele = reference.into_bytes();
    if variant.record_allele_string_multi.is_some() {
        variant.record_allele_string_multi = Some(written);
    }
}

/// The class of the rewritten alleles. The reference row keeps the record's
/// class: its allele equals the reference, which `classify_variant` would read
/// as a substitution, while Ensembl VEP's reference overlap allele shares the
/// feature overlaps of the record's own class (OutputFactory.pm:509-516), an
/// insertion's included.
fn reclassify(variant: &mut InputVariant) {
    if variant.reference_allele_row {
        return;
    }
    variant.variant_class = classify_variant(&variant.ref_allele, variant.alt_allele());
}

/// `--lookup_ref`'s replacement (Parser.pm:645): the reference allele of the
/// record is the reference sequence, `-` for an empty one, joined to the
/// alternate alleles as `validate_vf` holds them, each with its `-` removed
/// (`:578-580`) and `-` again when nothing is left, so `A-` becomes `A`; every
/// allele string the record prints follows; a record with no allele (a
/// non-variant line's) is the reference alone.
fn replace_reference(variant: &mut InputVariant, reference: &[u8]) {
    let reference: Vec<u8> = if reference.is_empty() {
        b"-".to_vec()
    } else {
        reference.to_vec()
    };
    let reference_text = String::from_utf8_lossy(&reference).into_owned();
    // The reference row's allele is the reference (`VariationFeatureOverlap.pm:469`
    // builds it from the replaced string).
    if variant.reference_allele_row {
        variant.alt_alleles = vec![reference.clone()];
    } else {
        for alt in &mut variant.alt_alleles {
            *alt = without_dashes(alt);
        }
    }
    variant.allele_string = if variant.no_allele {
        reference_text.clone()
    } else {
        format!(
            "{reference_text}/{}",
            String::from_utf8_lossy(variant.alt_allele())
        )
    };
    if let Some(multi) = &variant.record_allele_string_multi {
        // The parts after the reference as Perl's `split` reads them (Parser.pm:575),
        // a trailing empty part dropped: a one-part string (a reference-homozygous
        // copy's) or one ending in `/` (a copy whose genotype names no allele)
        // rejoins as the reference alone.
        let mut parts = vec![reference_text];
        parts.extend(
            perl_split(multi, &['/'])
                .skip(1)
                .map(|alt| String::from_utf8_lossy(&without_dashes(alt.as_bytes())).into_owned()),
        );
        variant.record_allele_string_multi = Some(parts.join("/").into());
    }
    variant.ref_allele = reference;
    reclassify(variant);
}

/// `allele` with its `-` removed (Parser.pm:578-580), `-` when nothing is left
/// (`map {$_ || '-'}`, `:645`).
fn without_dashes(allele: &[u8]) -> Vec<u8> {
    let stripped: Vec<u8> = allele.iter().copied().filter(|&b| b != b'-').collect();
    if stripped.is_empty() {
        b"-".to_vec()
    } else {
        stripped
    }
}

/// Parser/VCF.pm:555: a VCF record whose type reads as a deletion (`/del/i` on the
/// term `so_term`) is incomplete when its start lies past its end or its INFO
/// carries neither SVLEN nor END, read from the column as the reader parsed it
/// (`info`); a line that stops before INFO carries neither key. The test is the
/// VCF parser's alone: a record of another format (`info` is `None`) carries no
/// INFO and Ensembl VEP's other parsers make no such test (Parser/VEP_input.pm:136-148
/// builds a structural feature from a `DEL` word as written).
fn incomplete_deletion(
    variant: &InputVariant,
    so_term: &str,
    info: Option<&RecordInfo<'_>>,
) -> bool {
    let Some(info) = info else {
        return false;
    };
    deletion_looks_incomplete(
        so_term,
        variant.start,
        variant.end,
        info.value_written("SVLEN"),
        info.value_written("END"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::report_to;
    use crate::vcf_parser::{parse_data_line, VcfParse};
    use crate::warnings::{Sink, Target};
    use std::path::Path;

    /// A FASTA of one record `21` whose bases 1 to 20 are `ACGTACGTACGTACGTACGT`,
    /// with its index, in `dir`.
    fn fasta_in(dir: &Path) -> IndexedFasta {
        let path = dir.join("ref.fa");
        std::fs::write(&path, ">21\nACGTACGTACGTACGTACGT\n").unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "21\t20\t4\t20\t21\n").unwrap();
        IndexedFasta::from_path(&path).unwrap()
    }

    /// A transcript set holding chromosome 21 alone.
    fn transcripts(dir: &Path) -> LazyTranscriptIndexes {
        let cache = dir.join("cache");
        std::fs::create_dir_all(cache.join("transcripts").join("21")).unwrap();
        std::fs::write(
            cache.join("transcripts").join("21").join("1-100.json"),
            "[]",
        )
        .unwrap();
        LazyTranscriptIndexes::new(cache.to_str().unwrap(), 5000, 5000).unwrap()
    }

    fn variant(chr: &str, start: u64, end: u64, reference: &str, alt: &str) -> InputVariant {
        let mut v = InputVariant::new(
            chr.to_string(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        v.id = Some("r1".to_string());
        v
    }

    /// Runs the rules over one record and returns the verdict with the warning
    /// lines written for it. The record reports to a sink of its own, so a
    /// warning another test writes to the process-wide sink meanwhile is not
    /// read as this record's.
    fn judge(
        rules: &SkipRules<'_>,
        first: &mut InputVariant,
        dir: &Path,
    ) -> (bool, Vec<String>, AllelePatch) {
        let path = dir.join("w.txt");
        let _ = std::fs::remove_file(&path);
        let sink = Sink::new();
        sink.configure(Target::File(path.clone()), true, false);
        let line = RecordLine {
            number: 7,
            info: None,
            quoted: b"21\t100\tr1\tA\tG\t.\t.\t.",
        };
        let mut report = |reported: Reported<'_>| report_to(&sink, reported);
        let (kept, patch) = rules.validate(first, &line, &mut report);
        let lines = std::fs::read_to_string(&path)
            .unwrap_or_default()
            .lines()
            .map(str::to_string)
            .collect();
        (kept, lines, patch)
    }

    /// The reference checks of a run with the FASTA `fasta` and no `--chr` list.
    fn checks(fasta: &IndexedFasta, check_ref: bool, lookup_ref: bool) -> ReferenceChecks<'_> {
        ReferenceChecks {
            chr: None,
            check_ref,
            lookup_ref,
            fasta: Some(fasta),
            parser_names_records: false,
        }
    }

    static NO_SYNONYMS: std::sync::LazyLock<ChromosomeSynonyms> =
        std::sync::LazyLock::new(ChromosomeSynonyms::new);

    fn rules<'a>(
        transcripts: &'a LazyTranscriptIndexes,
        checks: &'a ReferenceChecks<'a>,
        dont_skip: bool,
    ) -> SkipRules<'a> {
        SkipRules {
            transcripts,
            synonyms: &NO_SYNONYMS,
            max_sv_size: None,
            dont_skip,
            checks,
            minimal: false,
        }
    }

    /// `--chr` admits the names its tokens spell and reads an integer range as
    /// the integers between its bounds, as Ensembl VEP's own help names the shape
    /// (`1-5,10,MT`), without naming each one; a token Ensembl VEP would match
    /// literally stays a name of its own.
    #[test]
    fn chr_range_expands() {
        let list = ChromosomeList::new(&["1-3".to_string(), "21".to_string(), "MT".to_string()]);
        for admitted in ["1", "2", "3", "21", "MT"] {
            assert!(list.admits(admitted), "{admitted}");
        }
        for refused in ["4", "chr21", "X", "01", ""] {
            assert!(!list.admits(refused), "{refused}");
        }
        let wide = ChromosomeList::new(&["1-18446744073709551615".to_string()]);
        assert!(wide.admits("18446744073709551615") && !wide.admits("X"));
        let literal = ChromosomeList::new(&[
            "X-Y".to_string(),
            "5-3".to_string(),
            "+1-3".to_string(),
            String::new(),
        ]);
        assert!(literal.admits("X-Y") && literal.admits("5-3") && literal.admits("+1-3"));
        assert!(!literal.admits("X") && !literal.admits("4") && !literal.admits(""));
        assert!(!literal.admits("2"), "a signed bound is no range");
    }

    /// The list filters on the input's own spelling (Parser.pm:470-472 compares
    /// `$vf->{chr}` before `_have_chr` renames or resolves it): `21` admits `21`,
    /// drops `chr21` and `22` without a warning, and `--dont_skip` keeps every
    /// record, as it keeps any record `validate_vf` refuses.
    #[test]
    fn chr_list_filters_and_chr_prefix_not_normalised() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let admitted = ChromosomeList::new(&["21".to_string()]);
        let checks = ReferenceChecks {
            chr: Some(&admitted),
            check_ref: false,
            lookup_ref: false,
            fasta: None,
            parser_names_records: false,
        };
        let rules = rules(&transcripts, &checks, false);
        for (chr, kept) in [("21", true), ("chr21", false), ("22", false)] {
            let mut v = variant(chr, 5, 5, "A", "G");
            let (verdict, warnings, _) = judge(&rules, &mut v, dir.path());
            assert_eq!(verdict, kept, "{chr}");
            assert!(warnings.is_empty(), "{chr}: {warnings:?}");
        }
        let mut renamed = variant("MT", 5, 5, "A", "G");
        renamed.input_chr = Some("M".into());
        assert!(
            !judge(&rules, &mut renamed, dir.path()).0,
            "M is matched, not MT"
        );
        let keep_all = SkipRules {
            dont_skip: true,
            ..rules
        };
        let mut v = variant("22", 5, 5, "A", "G");
        assert!(judge(&keep_all, &mut v, dir.path()).0);
    }

    /// Parser.pm:483-488: a start past the end plus one is refused with Ensembl
    /// VEP's text, under `--dont_skip` too.
    #[test]
    fn start_gt_end_plus_one() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let rules = rules(&transcripts, &NO_CHECKS, true);
        let mut v = variant("21", 12, 5, "A", "G");
        let (kept, warnings, _) = judge(&rules, &mut v, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): start > end+1 : (START=12, END=5)"]
        );
    }

    /// Parser.pm:564-572: an insertion whose coordinates are not start = end + 1
    /// is refused with Ensembl VEP's text, the record's name quoted and the
    /// trailing blank of the Perl concatenation kept; `--dont_skip` keeps it.
    #[test]
    fn insertion_needs_start_eq_end_plus_one() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let rules = rules(&transcripts, &NO_CHECKS, false);
        let mut v = variant("21", 5, 5, "-", "A");
        let (kept, warnings, _) = judge(&rules, &mut v, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Alleles look like an insertion (-/A) but coordinates are not start = end + 1 (START=5, END=5) for variant (r1) "]
        );
        let mut proper = variant("21", 6, 5, "-", "A");
        assert!(judge(&rules, &mut proper, dir.path()).0);
        let keep = SkipRules {
            dont_skip: true,
            ..rules
        };
        let mut v = variant("21", 5, 5, "-", "A");
        assert!(judge(&keep, &mut v, dir.path()).0);
    }

    /// Parser.pm:585-592: a reference allele of `ACGT` bases whose length is not
    /// the span's is refused with Ensembl VEP's text; one carrying `N` is not
    /// measured.
    #[test]
    fn ref_length_must_equal_span() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let rules = rules(&transcripts, &NO_CHECKS, false);
        let mut v = variant("21", 5, 10, "A", "G");
        let (kept, warnings, _) = judge(&rules, &mut v, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Length of reference allele (A length 1) does not match coordinates 5-10"]
        );
        let mut n = variant("21", 5, 10, "NN", "G");
        assert!(judge(&rules, &mut n, dir.path()).0);
        let mut deletion = variant("21", 5, 7, "ACG", "-");
        assert!(judge(&rules, &mut deletion, dir.path()).0);
    }

    /// Parser.pm:556-560: an allele string with no run of `ACGTN-` is refused with
    /// the upper-cased string in Ensembl VEP's text; `--dont_skip` keeps it.
    #[test]
    fn allele_regex() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let rules = rules(&transcripts, &NO_CHECKS, false);
        let mut v = variant("21", 5, 5, "x", "y");
        let (kept, warnings, _) = judge(&rules, &mut v, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Invalid allele string X/Y or possible parsing error"]
        );
        let keep = SkipRules {
            dont_skip: true,
            ..rules
        };
        let mut v = variant("21", 5, 5, "X", "Y");
        assert!(judge(&keep, &mut v, dir.path()).0);
    }

    /// Parser.pm:554: the alleles are upper-cased before the allele tests, on the
    /// first allele and, through the patch, on the record's others, the shared
    /// multi-allelic string included; the class follows the upper-cased bases.
    #[test]
    fn uppercase() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let rules = rules(&transcripts, &NO_CHECKS, false);
        let mut first = variant("21", 5, 5, "a", "g");
        first.record_allele_string_multi = Some("a/g/t".into());
        first.id = None;
        let (kept, warnings, patch) = judge(&rules, &mut first, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(first.allele_string, "A/G");
        // Parser.pm:462 names the record before :554 upper-cases it.
        assert_eq!(first.variation_name(), "21_5_a/g/t");
        assert_eq!(first.ref_allele, b"A");
        assert_eq!(first.alt_alleles, [b"G".to_vec()]);
        assert_eq!(first.record_allele_string_multi.as_deref(), Some("A/G/T"));
        let mut sibling = variant("21", 5, 5, "a", "t");
        sibling.record_allele_string_multi = Some("a/g/t".into());
        sibling.id = None;
        patch.apply(&mut sibling);
        assert_eq!(sibling.allele_string, "A/T");
        assert_eq!(sibling.record_allele_string_multi.as_deref(), Some("A/G/T"));
        assert_eq!(sibling.variation_name(), "21_5_a/g/t");
        assert!(!sibling.check_ref_failed);
        let mut upper = variant("21", 5, 5, "A", "G");
        upper.id = None;
        assert!(judge(&rules, &mut upper, dir.path()).0);
        assert_eq!(
            upper.id, None,
            "an upper-case record keeps building its name at output"
        );
        // A parser's working copy of the multi-allelic string arrives upper-cased
        // with the chopped string in its written case beside it: the lower-case
        // base reaches the sibling it belongs to, and the name reads the alleles
        // as written (Ensembl VEP 116.2 names `21 30319403 0 tca t,tcg`
        // `21_30319404_ca/-/cg`).
        let mut first = variant("21", 5, 5, "A", "G");
        first.record_allele_string_multi = Some("A/G/T".into());
        first.written_case_allele_string = Some("A/G/t".into());
        first.uploaded_allele_string = Some("A/G/t".into());
        first.id = None;
        let (kept, warnings, patch) = judge(&rules, &mut first, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(first.variation_name(), "21_5_A/G/t");
        assert_eq!(first.record_allele_string_multi.as_deref(), Some("A/G/T"));
        let mut sibling = variant("21", 5, 5, "A", "t");
        sibling.record_allele_string_multi = Some("A/G/T".into());
        sibling.written_case_allele_string = Some("A/G/t".into());
        sibling.uploaded_allele_string = Some("A/G/t".into());
        sibling.id = None;
        patch.apply(&mut sibling);
        assert_eq!(sibling.allele_string, "A/T");
        assert_eq!(sibling.alt_alleles, [b"T".to_vec()]);
        assert_eq!(sibling.variation_name(), "21_5_A/G/t");
        // The bi-allelic record of the VCF parser arrives the same way, its one
        // pair chopped and upper-cased with the written case beside it: Ensembl
        // VEP 116.2 names `21 30319403 0 tca t` `21_30319404_ca/-` and prints
        // that name in every format, its alleles `CA/-`.
        let mut deletion = variant("21", 30319404, 30319405, "CA", "-");
        deletion.written_case_allele_string = Some("ca/-".into());
        deletion.uploaded_allele_string = Some("tca/t".into());
        deletion.id = Some("0".to_string());
        let (kept, warnings, _) = judge(&rules, &mut deletion, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(deletion.id.as_deref(), Some("21_30319404_ca/-"));
        assert_eq!(deletion.variation_name(), "21_30319404_ca/-");
        assert_eq!(deletion.uploaded_variation(), "21_30319404_ca/-");
        assert_eq!(deletion.allele_string, "CA/-");
        assert_eq!(deletion.ref_allele, b"CA");
    }

    /// The `--chr` list (Parser.pm:470-472) and the chromosome test (`:542-546`)
    /// return before the upper-casing of `:554`, so a record `--dont_skip` keeps
    /// after either refused it is annotated with its alleles as written, on the
    /// first allele and through the patch on its siblings, the shared string
    /// included: Ensembl VEP 116.2 prints Allele `g` and JSON `allele_string`
    /// `a/g` for `Z 1000 . a g` under `--dont_skip`, and Allele `t` with
    /// `allele_string` `c/t` for `21 30337048 . c t` under `--chr 22 --dont_skip`.
    /// A record written in upper case has nothing to restore, and a check after
    /// `:554` that refuses a kept record leaves it upper-cased.
    #[test]
    fn dont_skip_keeps_the_written_case_before_the_upper_casing() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let unlisted = rules(&transcripts, &NO_CHECKS, true);
        let mut first = variant("Z", 1000, 1000, "A", "G");
        first.written_case_allele_string = Some("a/g".into());
        first.uploaded_allele_string = Some("a/g".into());
        first.id = None;
        let (kept, warnings, _) = judge(&unlisted, &mut first, dir.path());
        assert!(kept && warnings.len() == 1, "{warnings:?}");
        assert_eq!(first.allele_string, "a/g");
        assert_eq!(first.ref_allele, b"a");
        assert_eq!(first.alt_alleles, [b"g".to_vec()]);
        assert_eq!(first.variation_name(), "Z_1000_a/g");

        let admitted = ChromosomeList::new(&["22".to_string()]);
        let checks = ReferenceChecks {
            chr: Some(&admitted),
            ..NO_CHECKS
        };
        let listed = rules(&transcripts, &checks, true);
        let mut first = variant("21", 5, 5, "C", "T");
        first.record_allele_string_multi = Some("C/T/G".into());
        first.written_case_allele_string = Some("c/t/G".into());
        first.uploaded_allele_string = Some("c/t/G".into());
        let (kept, warnings, patch) = judge(&listed, &mut first, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(first.allele_string, "c/t");
        assert_eq!(first.record_allele_string_multi.as_deref(), Some("c/t/G"));
        let mut sibling = variant("21", 5, 5, "C", "G");
        sibling.allele_index = 1;
        sibling.record_allele_string_multi = Some("C/T/G".into());
        sibling.written_case_allele_string = Some("c/t/G".into());
        sibling.uploaded_allele_string = Some("c/t/G".into());
        patch.apply(&mut sibling);
        assert_eq!(sibling.allele_string, "c/G");
        assert_eq!(sibling.ref_allele, b"c");
        assert_eq!(sibling.alt_alleles, [b"G".to_vec()]);
        assert_eq!(sibling.record_allele_string_multi.as_deref(), Some("c/t/G"));

        let mut upper = variant("21", 5, 5, "C", "T");
        assert!(judge(&listed, &mut upper, dir.path()).0);
        assert_eq!(upper.allele_string, "C/T");
        // `x/y` fails the allele-string test (`:556-560`), after `:554` upper-cased it.
        let mut invalid = variant("21", 5, 5, "X", "Y");
        invalid.written_case_allele_string = Some("x/y".into());
        invalid.uploaded_allele_string = Some("x/y".into());
        let (kept, warnings, _) = judge(&unlisted, &mut invalid, dir.path());
        assert!(kept && warnings.len() == 1, "{warnings:?}");
        assert_eq!(invalid.allele_string, "X/Y");
    }

    /// A parser that names its records itself (Parser/Region.pm:198,
    /// Parser/HGVS.pm:237) gets no name from Parser.pm:462, so a record it left
    /// unnamed keeps no alleles-built name through the upper-casing and the
    /// looked-up reference, on the first allele and through the patch.
    #[test]
    fn parser_named_records_take_no_allele_name() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = ReferenceChecks {
            parser_names_records: true,
            ..checks(&fasta, false, true)
        };
        let rules = rules(&transcripts, &checks, false);
        let mut first = variant("21", 5, 6, "nn", "g");
        first.id = None;
        let (kept, warnings, patch) = judge(&rules, &mut first, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(first.id, None);
        assert_eq!(first.allele_string, "AC/G");
        let mut sibling = variant("21", 5, 6, "nn", "t");
        sibling.id = None;
        patch.apply(&mut sibling);
        assert_eq!(sibling.id, None);
        assert_eq!(sibling.allele_string, "AC/T");
    }

    /// `check_ref` (Parser.pm:596-633): a reference allele that differs from the
    /// FASTA's bases over the span fails the record with Ensembl VEP's text, which
    /// never prints the expected base; the record is dropped, or kept under
    /// `--dont_skip` with `check_ref_failed` on every allele; a span the FASTA
    /// does not hold writes the sub-slice line first.
    #[test]
    fn check_ref_mismatch_skips() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, true, false);
        let rules = rules(&transcripts, &checks, false);
        // Bases 1 to 4 are `ACGT`: `AC` at 1-2 matches, `GG` at 5-6 (`AC`) does not.
        let mut good = variant("21", 1, 2, "AC", "A");
        assert!(judge(&rules, &mut good, dir.path()).0);
        assert!(!good.check_ref_failed);
        let mut bad = variant("21", 5, 6, "GG", "G");
        let (kept, warnings, _) = judge(&rules, &mut bad, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Specified reference allele GG does not match Ensembl reference allele"]
        );
        let keep = SkipRules {
            dont_skip: true,
            ..rules
        };
        let mut kept = variant("21", 5, 6, "GG", "G");
        let (verdict, _, patch) = judge(&keep, &mut kept, dir.path());
        assert!(verdict && kept.check_ref_failed);
        let mut sibling = variant("21", 5, 6, "GG", "T");
        patch.apply(&mut sibling);
        assert!(sibling.check_ref_failed);
        let mut past_end = variant("21", 25, 25, "A", "G");
        let (_, warnings, _) = judge(&rules, &mut past_end, dir.path());
        assert_eq!(
            warnings,
            [
                "WARNING: line 7 skipped (21 100 r1 A G . . .): Could not fetch sub-slice from 21:25-25(1)",
                "WARNING: line 7 skipped (21 100 r1 A G . . .): Specified reference allele A does not match Ensembl reference allele",
            ]
        );
    }

    /// An insertion has no reference allele to check (Parser.pm:601-603), so it
    /// passes `check_ref` whatever the flanking bases.
    #[test]
    fn check_ref_insertion_passes() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, true, false);
        let rules = rules(&transcripts, &checks, false);
        let mut insertion = variant("21", 6, 5, "-", "TTT");
        let (kept, warnings, _) = judge(&rules, &mut insertion, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert!(!insertion.check_ref_failed);
    }

    /// `uc($slice_ref_allele) ne uc($ref_allele)` (Parser.pm:614): the comparison
    /// is case-insensitive, and the record's alleles are upper-cased before it; a
    /// reverse-strand record is checked against the reverse complement
    /// (`sub_Slice` on the record's strand, Parser.pm:751).
    #[test]
    fn check_ref_case_insensitive() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, true, false);
        let rules = rules(&transcripts, &checks, false);
        let mut lower = variant("21", 1, 2, "ac", "a");
        let (kept, warnings, _) = judge(&rules, &mut lower, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(lower.allele_string, "AC/A");
        // Bases 2 to 3 are `CG`, whose reverse complement is `CG`; 1 to 2 are `AC`,
        // reverse complement `GT`.
        let mut reverse = variant("21", 1, 2, "GT", "G");
        reverse.strand = Strand::Reverse;
        assert!(judge(&rules, &mut reverse, dir.path()).0);
        let mut reverse_bad = variant("21", 1, 2, "AC", "A");
        reverse_bad.strand = Strand::Reverse;
        let (kept, warnings, _) = judge(&rules, &mut reverse_bad, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Specified reference allele AC does not match Ensembl reference allele"]
        );
    }

    /// An ALT the looked-up reference now equals writes no row and takes no number
    /// (VariationFeatureOverlap.pm:481): Ensembl VEP 116.2 writes `21 10542415 x T A,G`
    /// over the reference `A` as the one row `G` with ALLELE_NUM 1, `T G,A` the same,
    /// `T A,G,C` as `G` 1 and `C` 2, and `T A` as the feature with no allele (`A/A`,
    /// `?`); an allele the reference does not equal keeps its row.
    #[test]
    fn lookup_ref_drops_an_alt_the_reference_now_equals() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, false, true);
        let rules = rules(&transcripts, &checks, false);
        // Position 1 of the FASTA is `A`.
        let mut first = variant("21", 1, 1, "T", "A");
        first.record_allele_string_multi = Some("T/A/G".into());
        let (kept, _, patch) = judge(&rules, &mut first, dir.path());
        assert!(kept && patch.replaced_reference());
        assert!(!keep_overlap_allele(&mut first, false));
        let mut second = variant("21", 1, 1, "T", "G");
        second.allele_index = 1;
        second.record_allele_string_multi = Some("T/A/G".into());
        patch.apply(&mut second);
        assert!(keep_overlap_allele(&mut second, true));
        assert_eq!(second.allele_string, "A/G");
        assert_eq!(second.allele_num, Some(1));
        assert_eq!(second.uploaded_allele(false).as_deref(), Some("A/A/G"));
        let mut third = variant("21", 1, 1, "T", "C");
        third.allele_index = 2;
        third.record_allele_string_multi = Some("T/A/G/C".into());
        patch.apply(&mut third);
        assert!(keep_overlap_allele(&mut third, true));
        assert_eq!(third.allele_num, Some(2));
        let mut alone = variant("21", 1, 1, "T", "A");
        assert!(judge(&rules, &mut alone, dir.path()).0);
        assert!(keep_overlap_allele(&mut alone, true));
        assert!(alone.no_allele);
        assert_eq!(alone.allele_string, "A/A");
        let mut kept_alt = variant("21", 1, 1, "T", "G");
        assert!(judge(&rules, &mut kept_alt, dir.path()).0);
        assert!(keep_overlap_allele(&mut kept_alt, true));
        assert!(!kept_alt.no_allele);
    }

    /// A `*` or `<DEL:*>` allele of a parser-fed record takes its number and writes
    /// no row (VariationFeatureOverlap.pm:483-484): Ensembl VEP 116.2 writes the
    /// Ensembl-format line `21 30282169 30282169 C/*/A` as the one row `A` with
    /// ALLELE_NUM 2, `C/A/*` as `A` with ALLELE_NUM 1, and `C/*`, `C/<DEL:*>` as
    /// the feature with no allele (no default row, JSON `?`).
    #[test]
    fn star_and_delstar_alleles_write_no_row() {
        let mut star = variant("21", 30282169, 30282169, "C", "*");
        star.record_allele_string_multi = Some("C/*/A".into());
        assert!(rowless_allele(&star));
        assert!(!keep_overlap_allele(&mut star, false));
        let mut after_star = variant("21", 30282169, 30282169, "C", "A");
        after_star.allele_index = 1;
        after_star.record_allele_string_multi = Some("C/*/A".into());
        assert!(!rowless_allele(&after_star));
        assert!(keep_overlap_allele(&mut after_star, true));
        assert!(!after_star.no_allele);
        assert_eq!(after_star.allele_num, None);
        let mut trailing_star = variant("21", 30282169, 30282169, "C", "*");
        trailing_star.allele_index = 1;
        trailing_star.record_allele_string_multi = Some("C/A/*".into());
        assert!(!keep_overlap_allele(&mut trailing_star, false));
        for alt in ["*", "<DEL:*>"] {
            let mut alone = variant("21", 30282169, 30282169, "C", alt);
            assert!(rowless_allele(&alone));
            assert!(keep_overlap_allele(&mut alone, true));
            assert!(alone.no_allele, "{alt}");
            assert_eq!(alone.alt_alleles, [b"*".to_vec()]);
        }
        let mut ordinary = variant("21", 30282169, 30282169, "C", "A");
        assert!(keep_overlap_allele(&mut ordinary, true));
        assert!(!ordinary.no_allele);
    }

    /// The reference row of `--process_ref_homs` stands under `--lookup_ref`
    /// (VariationFeatureOverlap.pm:469-474 builds it before any alternate is
    /// tested against the reference): its allele is the looked-up reference,
    /// numbered 0, and a record whose one ALT that reference now equals keeps the
    /// row alone, the ALT dropped behind it.
    #[test]
    fn lookup_ref_keeps_the_reference_row() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, false, true);
        let rules = rules(&transcripts, &checks, false);
        // Position 1 of the FASTA is `A`, so the one ALT of `T A` equals it.
        let mut alt = variant("21", 1, 1, "T", "A");
        let mut row = alt.reference_allele_row();
        let (kept, warnings, patch) = judge(&rules, &mut row, dir.path());
        assert!(kept && warnings.is_empty() && patch.replaced_reference());
        assert_eq!(row.alt_alleles, [b"A".to_vec()]);
        assert_eq!(row.allele_string, "A/A");
        assert!(keep_overlap_allele(&mut row, false));
        assert!(row.reference_allele_row && !row.no_allele);
        assert_eq!(row.allele_num, Some(0));
        patch.apply(&mut alt);
        assert!(!keep_overlap_allele(&mut alt, false));
        let mut multi = variant("21", 1, 1, "T", "A");
        multi.record_allele_string_multi = Some("T/A/G".into());
        let mut row = multi.reference_allele_row();
        assert!(judge(&rules, &mut row, dir.path()).0);
        assert!(keep_overlap_allele(&mut row, false));
        assert_eq!(row.allele_num, Some(0));
        assert_eq!(row.record_allele_string_multi.as_deref(), Some("A/A/G"));
        // A reference-homozygous copy without `--process_ref_homs` is a feature
        // with no allele carrying the written reference alone; the join of
        // Parser.pm:645 over that one part is the looked-up reference alone
        // (Ensembl VEP 116.2 prints the JSON `allele_string` `T` for a `0|0` copy
        // of `C A` over a `T`), and a trailing empty part, a copy whose genotype
        // names no allele, is no part to Perl's `split` (`:575`).
        for working in ["T", "T/"] {
            let mut hom_ref = variant("21", 1, 1, "T", "*");
            hom_ref.no_allele = true;
            hom_ref.hom_ref = true;
            hom_ref.allele_string = "T".to_string();
            hom_ref.record_allele_string_multi = Some(working.into());
            assert!(judge(&rules, &mut hom_ref, dir.path()).0, "{working}");
            assert_eq!(
                hom_ref.record_allele_string_multi.as_deref(),
                Some("A"),
                "{working}"
            );
            assert_eq!(hom_ref.allele_string, "A", "{working}");
        }
    }

    /// `lookup_ref` (Parser.pm:635-646): the reference allele is replaced by the
    /// FASTA's bases over the span, on the record's strand, every allele string
    /// the record prints following, and the record is kept; the class follows the
    /// new alleles.
    #[test]
    fn lookup_ref_replaces_ref() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, false, true);
        let rules = rules(&transcripts, &checks, false);
        let mut first = variant("21", 5, 6, "GG", "G");
        first.record_allele_string_multi = Some("GG/G/T".into());
        let (kept, warnings, patch) = judge(&rules, &mut first, dir.path());
        assert!(kept && warnings.is_empty(), "{warnings:?}");
        assert_eq!(first.ref_allele, b"AC");
        assert_eq!(first.allele_string, "AC/G");
        assert_eq!(first.record_allele_string_multi.as_deref(), Some("AC/G/T"));
        let mut sibling = variant("21", 5, 6, "GG", "T");
        sibling.record_allele_string_multi = Some("GG/G/T".into());
        patch.apply(&mut sibling);
        assert_eq!(sibling.allele_string, "AC/T");
        assert_eq!(
            sibling.record_allele_string_multi.as_deref(),
            Some("AC/G/T")
        );
        let mut reverse = variant("21", 1, 2, "NN", "A");
        reverse.strand = Strand::Reverse;
        assert!(judge(&rules, &mut reverse, dir.path()).0);
        assert_eq!(reverse.allele_string, "GT/A");
        // The record keeps the name built from the reference as written (:462
        // runs before :645).
        let mut unnamed = variant("21", 5, 6, "GG", "G");
        unnamed.id = None;
        assert!(judge(&rules, &mut unnamed, dir.path()).0);
        assert_eq!(unnamed.allele_string, "AC/G");
        assert_eq!(unnamed.variation_name(), "21_5_GG/G");
        let mut past_end = variant("21", 25, 25, "A", "G");
        let (kept, warnings, _) = judge(&rules, &mut past_end, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Could not fetch reference allele for 21:25-25(1)"]
        );
        // A non-variant record's allele string is its reference alone (Parser/VCF.pm:342),
        // and stays so once looked up.
        let mut non_variant =
            InputVariant::without_allele("21".into(), 5, "GG", std::iter::empty());
        assert!(judge(&rules, &mut non_variant, dir.path()).0);
        assert_eq!(non_variant.allele_string, "AC");
        assert_eq!(non_variant.ref_allele, b"AC");
    }

    /// An insertion's empty reference reads back as `-` (`map {$_ || '-'}`,
    /// Parser.pm:645), so `-/A` stays `-/A`.
    #[test]
    fn lookup_ref_restores_dash() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, false, true);
        let rules = rules(&transcripts, &checks, false);
        let mut insertion = variant("21", 6, 5, "-", "A");
        assert!(judge(&rules, &mut insertion, dir.path()).0);
        assert_eq!(insertion.allele_string, "-/A");
        assert_eq!(insertion.ref_allele, b"-");
    }

    /// `check_ref` wins over `lookup_ref` when both are set (`if ... elsif`,
    /// Parser.pm:596 and 635), as under `--safe`, which skips the pair check.
    #[test]
    fn check_ref_takes_precedence_over_lookup_ref() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let fasta = fasta_in(dir.path());
        let checks = checks(&fasta, true, true);
        let rules = rules(&transcripts, &checks, false);
        let mut bad = variant("21", 5, 6, "GG", "G");
        assert!(!judge(&rules, &mut bad, dir.path()).0);
        assert_eq!(bad.allele_string, "GG/G");
    }

    /// Parser.pm:476-480's text for a VCF POS that is not a number: the end is
    /// the numified POS plus the reference length less one (Parser/VCF.pm:254).
    /// `looks_like_number` reads what Perl's reads as a number, so a negative or
    /// fractional POS, which Ensembl VEP annotates as written, is not this message.
    #[test]
    fn non_numeric_coords() {
        for number in [
            "1", "-1", "+5", "3.5", ".5", "3.", "1e3", "1E-3", " 12 ", "Inf", "-inf", "NaN",
        ] {
            assert!(looks_like_number(number), "{number:?}");
        }
        for text in [
            "abc", "12abc", "1_000", "0x1A", "", ".", "-", "1e", "e3", "1.2.3", "30 425",
        ] {
            assert!(!looks_like_number(text), "{text:?}");
        }
        assert_eq!(
            invalid_vcf_coordinate("abc", "A"),
            "Invalid start 'abc' or end '0' coordinate"
        );
        assert_eq!(
            invalid_vcf_coordinate("12abc", "ACG"),
            "Invalid start '12abc' or end '14' coordinate"
        );
        assert_eq!(
            invalid_vcf_coordinate("1_000", "A"),
            "Invalid start '1_000' or end '1' coordinate"
        );
    }

    /// `validate_vf` transforms a chromosome the cache lacks to the toplevel
    /// coordinate system through the core database (Parser.pm:505-536); offline it
    /// writes the one message (`:542-546`), the branch this port has, so the
    /// database fallback is the documented exclusion.
    #[test]
    fn database_toplevel_transform_not_ported() {
        let dir = tempfile::tempdir().unwrap();
        let transcripts = transcripts(dir.path());
        let rules = rules(&transcripts, &NO_CHECKS, false);
        let mut v = variant("Z", 5, 5, "A", "G");
        let (kept, warnings, _) = judge(&rules, &mut v, dir.path());
        assert!(!kept);
        assert_eq!(
            warnings,
            ["WARNING: line 7 skipped (21 100 r1 A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features"]
        );
        assert!(
            !warnings[0].contains("toplevel"),
            "no transform is attempted offline"
        );
    }

    /// The incomplete-deletion test is the VCF parser's (Parser/VCF.pm:555-558),
    /// read from the INFO column as the reader parsed the record: a `<DEL>` line
    /// with neither SVLEN nor END, or with an END at POS, is incomplete, one with
    /// an END past POS is not, a line that stops before INFO carries neither key,
    /// and the column copied out of a line read through a replacement copy
    /// answers as the parse's view does. A deletion of another format (`info`
    /// None), an Ensembl-format `21 30165851 30165850 DEL + x` whose start is its
    /// end plus one, is never incomplete, since Ensembl VEP 116.2 annotates it
    /// (`feature_truncation` on ENST00000286808) and Parser/VEP_input.pm:136-148
    /// makes no such test.
    #[test]
    fn the_incomplete_deletion_test_reads_a_vcf_line_alone() {
        let parse = VcfParse::default_run(false);
        let mut deletion = variant("21", 30165851, 30165851, "N", "<DEL>");
        deletion.is_structural = true;
        let line = |info: &str| format!("21\t30165850\tx\tN\t<DEL>\t.\t.\t{info}");
        for (so_term, info, incomplete) in [
            ("deletion", "SVTYPE=DEL", true),
            ("deletion", "END=30165851", false),
            ("deletion", "SVLEN=-1", false),
            ("duplication", "SVTYPE=DUP", false),
        ] {
            let text = line(info);
            let split = parse_data_line(&text, &parse);
            assert!(matches!(split.info, Some(RecordInfo::Split(_))));
            let copied = parse_data_line(&text, &parse).with_info_copied(&text);
            assert!(matches!(copied.info, Some(RecordInfo::Copied(_))));
            for record in [&split, &copied] {
                assert_eq!(
                    incomplete_deletion(&deletion, so_term, record.info.as_ref()),
                    incomplete,
                    "{so_term} {info}"
                );
            }
        }
        let short = parse_data_line("21\t30165850\tx\tN\t<DEL>", &parse);
        assert!(incomplete_deletion(
            &deletion,
            "deletion",
            short.info.as_ref()
        ));
        let pair = variant("21", 30165851, 30165850, "N", "<DEL>");
        let end_at_pos_text = line("END=30165850");
        let end_at_pos = parse_data_line(&end_at_pos_text, &parse);
        assert!(incomplete_deletion(
            &pair,
            "deletion",
            end_at_pos.info.as_ref()
        ));
        assert!(!incomplete_deletion(&pair, "deletion", None));
        assert!(!incomplete_deletion(&deletion, "deletion", None));
    }
}
