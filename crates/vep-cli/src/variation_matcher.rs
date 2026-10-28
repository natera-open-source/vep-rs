// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Co-located variant matching: the `--check_existing` stage.
//!
//! Ports Ensembl VEP release/116.2's matching of an input variant against the
//! known variants of its window: `AnnotationSource/Cache/VariationTabix.pm`
//! `_annotate_pm` (the walk, the failed filter, one entry per name),
//! `AnnotationType/Variation.pm` `compare_existing` (span matches, the allele
//! rules) and ensembl-variation release/116's `Utils/Sequence.pm`
//! `get_matched_variant_alleles` and `trim_sequences` (strand-aware allele
//! matching with trimming in both directions). The `Extra` column writer of the
//! default format lives here too.

use std::borrow::Cow;

use vep_core::consequence::FeatureType;
use vep_core::coordinate::Strand;
use vep_core::minimal::rowless_first_feature;
use vep_core::sv_routing::perl_true;
use vep_core::variant::{
    colocated_in_output_order, push_decimal, ColocatedVariant, InputVariant, MatchedAllele,
};
use vep_core::variation::{reverse_comp, CachedVariation};

use crate::variation_source::{ReadFailure, VariationSource};

/// The flags that decide which known variants of a window co-locate with an
/// input variant (`AnnotationType/Variation.pm:146-166, 224-227`); the default is
/// Ensembl VEP's: alleles checked, null alleles kept, `failed` 0.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ColocationRules {
    /// `--no_check_alleles`: a known variant matches on its start and end alone and
    /// prints on every allele's row.
    pub no_check_alleles: bool,
    /// `--exclude_null_alleles`: a known variant without alleles (`COSMIC_MUTATION`,
    /// `HGMD_MUTATION`, `dbSNP_novariation`) never matches; without the flag it
    /// matches on its start and end.
    pub exclude_null_alleles: bool,
    /// `--failed <int>`: a known variant matches only when its `failed` value is at
    /// most this; Ensembl VEP's default is 0.
    pub failed: i64,
}

impl ColocationRules {
    /// The rules a run's configuration asks for.
    pub fn from_config(config: &crate::config::Config) -> Self {
        Self {
            no_check_alleles: config.no_check_alleles,
            exclude_null_alleles: config.exclude_null_alleles,
            failed: config.failed,
        }
    }
}

/// Attaches the known variants of `source` co-located with `variant`:
/// `colocated_variants` holds every match of the input record, in the cache's
/// order, and `existing_variation` the names that print on this allele's rows,
/// in Ensembl VEP's order. A structural variant gets none (`VariationTabix.pm:118`).
///
/// An allele of a `--minimal` split record is matched as its own split feature
/// (InputBuffer.pm:265 splits the record before the sources run), unless the
/// record's first alternate allele writes no row, when the first feature is
/// matched in its place ([`rowless_first_feature`]); `runner::apply_record_rules`
/// then gives every allele of the record the first feature's matches, which are
/// the ones Ensembl VEP's rows print.
///
/// A failure to read a variation file is recorded in `failure` for the batch
/// driver to end the run with (`runner::check_variation_reads`): the file is
/// part of the cache, and a variant annotated without it would read as having
/// no known co-located variant.
pub fn annotate_colocated(
    variant: &mut InputVariant,
    source: &VariationSource,
    rules: &ColocationRules,
    failure: &ReadFailure,
) {
    if variant.is_structural {
        return;
    }
    let first_feature = rowless_first_feature(variant);
    let subject = first_feature.as_ref().unwrap_or(variant);
    let candidates = match source.records_overlapping(&subject.chr, subject.start, subject.end) {
        Ok(candidates) => candidates,
        Err(error) => {
            failure.record(format!("{error:#}"));
            return;
        }
    };
    let colocated = find_colocated_variants(subject, candidates, rules);
    variant.existing_variation = existing_variation_names(&colocated, &variant.display_allele());
    variant.colocated_variants = colocated;
}

/// The names of `colocated` that print on the rows of `allele`, in Ensembl VEP's
/// order: [`colocated_in_output_order`], keeping the records
/// [`ColocatedVariant::applies_to_allele`] admits.
pub(crate) fn existing_variation_names(
    colocated: &[ColocatedVariant],
    allele: &str,
) -> Vec<String> {
    colocated_in_output_order(colocated)
        .into_iter()
        .filter(|cv| cv.applies_to_allele(allele))
        .map(|cv| cv.id.clone())
        .collect()
}

/// The known variants among `candidates`, the records of the variant's window in
/// file order, that Ensembl VEP keeps on the input record
/// (`VariationTabix.pm:232-241`): a record whose `failed` value exceeds the
/// threshold is skipped, a name already kept is not matched again, and
/// `compare_existing` decides the rest. The record's whole allele string is
/// matched, as Ensembl VEP matches the variation feature of a multi-allelic
/// record before splitting its rows by allele; an allele of a `--minimal` split
/// record is matched by its own trimmed pair at its trimmed span, the feature
/// `split_variants` makes of it before the sources run (InputBuffer.pm:265,
/// 658-676).
pub fn find_colocated_variants(
    variant: &InputVariant,
    candidates: Vec<CachedVariation>,
    rules: &ColocationRules,
) -> Vec<ColocatedVariant> {
    let allele_string = if variant.original_start.is_some() {
        variant.allele_string.as_str()
    } else {
        variant
            .record_allele_string_multi
            .as_deref()
            .unwrap_or(&variant.allele_string)
    };
    let input = InputAlleles {
        allele_string,
        start: variant.start,
        end: variant.end,
        strand: variant.strand,
    };
    let mut colocated: Vec<ColocatedVariant> = Vec::new();
    for record in candidates {
        if i64::from(record.failed) > rules.failed {
            continue;
        }
        if colocated.iter().any(|cv| cv.id == record.variation_name) {
            continue;
        }
        if let Some(matched) = compare_existing(&input, &record, rules) {
            colocated.push(cached_to_colocated(record, matched));
        }
    }
    colocated
}

/// The input record as `compare_existing` sees its variation feature.
struct InputAlleles<'a> {
    /// The feature's allele string: every alternate allele after the shared
    /// anchor, so a multi-allelic record is matched whole, or a split feature's
    /// trimmed pair.
    allele_string: &'a str,
    start: u64,
    end: u64,
    strand: Strand,
}

/// Whether a known variant co-locates with the input record, and on which alleles
/// (`AnnotationType/Variation.pm:146-207`). `Some(empty)` is a match on start and
/// end alone, which Ensembl VEP returns without `matched_alleles`: a record
/// without alleles (no `/` in its allele string) unless `--exclude_null_alleles`,
/// or any record under `--no_check_alleles`. Otherwise the alleles decide.
fn compare_existing(
    input: &InputAlleles<'_>,
    record: &CachedVariation,
    rules: &ColocationRules,
) -> Option<Vec<MatchedAllele>> {
    let same_span = record.start == input.start && record.end == input.end;
    if !record.allele_string.contains('/') {
        if rules.exclude_null_alleles {
            return None;
        }
        if same_span {
            return Some(Vec::new());
        }
    }
    if rules.no_check_alleles {
        return if same_span { Some(Vec::new()) } else { None };
    }
    let matched = get_matched_variant_alleles(
        &AlleleSet {
            allele_string: input.allele_string,
            pos: input.start,
            strand: input.strand.as_i8(),
        },
        &AlleleSet {
            allele_string: &record.allele_string,
            pos: record.start,
            strand: record.strand,
        },
    );
    (!matched.is_empty()).then_some(matched)
}

/// A variant as `get_matched_variant_alleles` takes it: its allele string,
/// reference first, its start and its strand.
#[derive(Debug, Clone, Copy)]
pub struct AlleleSet<'a> {
    pub allele_string: &'a str,
    pub pos: u64,
    pub strand: i8,
}

/// The reference and alternate alleles of an allele string as `split('/')`
/// yields them (trailing empty alleles dropped, as Perl's `split` drops them);
/// `None` for an empty string, a string that is `0`, or one beginning with `/`
/// followed by anything (`Utils/Sequence.pm:1151-1157, 1173`).
fn split_allele_string(allele_string: &str) -> Option<(&str, std::str::Split<'_, char>)> {
    if !perl_true(allele_string) || (allele_string.starts_with('/') && allele_string.len() > 1) {
        return None;
    }
    let mut alleles = allele_string.trim_end_matches('/').split('/');
    let reference = alleles.next()?;
    Some((reference, alleles))
}

/// The alternate alleles `a` and `b` share, as `get_matched_variant_alleles`
/// finds them (`Utils/Sequence.pm:1144-1247`, ensembl-variation release/116):
/// `a`'s alleles are reverse-complemented when the strands differ, every
/// reference-alternate pair of each side is trimmed of its shared ends in both
/// orders (left first, then right first) when either allele is longer than one
/// base, and a pair of `b` matches when one of its trimmed forms, tried in order,
/// equals a trimmed form of `a` (the later of `a`'s alleles owning a shared form).
/// Returns one match per matching allele of `b`, in `b`'s order; a position of 0
/// or an unusable allele string matches nothing.
pub fn get_matched_variant_alleles(a: &AlleleSet<'_>, b: &AlleleSet<'_>) -> Vec<MatchedAllele> {
    let (Some((a_ref, a_alts)), Some((b_ref, b_alts))) = (
        split_allele_string(a.allele_string),
        split_allele_string(b.allele_string),
    ) else {
        return Vec::new();
    };
    if a.pos == 0 || b.pos == 0 {
        return Vec::new();
    }
    let flip = a.strand != b.strand;
    let a_ref = reverse_comp_if(flip, a_ref);
    // Each alternate allele of `a` as written, with its reading on `b`'s strand.
    let a_alts: Vec<(&str, Cow<'_, str>)> = a_alts
        .map(|alt| (alt, reverse_comp_if(flip, alt)))
        .collect();

    // Every trimmed form of an input allele with the allele as written; a later
    // allele with the same form replaces an earlier one, as the Perl hash does, so
    // a lookup reads the list from its end.
    let mut minimised: Vec<(TrimmedAlleles<'_>, &str)> = Vec::with_capacity(a_alts.len() * 2);
    for (orig_alt, alt) in &a_alts {
        for &end_first in trim_directions(&a_ref, orig_alt) {
            minimised.push((trim_sequences(&a_ref, alt, a.pos, end_first), *orig_alt));
        }
    }

    let mut matches = Vec::new();
    for (index, orig_alt) in b_alts.enumerate() {
        for &end_first in trim_directions(b_ref, orig_alt) {
            let trimmed = trim_sequences(b_ref, orig_alt, b.pos, end_first);
            let found = minimised
                .iter()
                .rev()
                .find(|(form, _)| *form == trimmed)
                .map(|(_, a_alt)| *a_alt)
                .filter(|a_alt| perl_true(a_alt));
            if let Some(a_alt) = found {
                matches.push(MatchedAllele {
                    input_allele: a_alt.to_string(),
                    // The index is the allele string's last position among `a`'s
                    // alternates, as `$a_indexes{$orig_a_alt} = $i++` leaves it.
                    input_index: a_alts
                        .iter()
                        .rposition(|(alt, _)| *alt == a_alt)
                        .expect("the allele came from this list"),
                    record_allele: orig_alt.to_string(),
                    record_index: index,
                });
                break;
            }
        }
    }
    matches
}

/// The trimming orders tried for an allele pair: left end first, and right end
/// first too when either allele is longer than one base
/// (`Utils/Sequence.pm:1250-1258`).
fn trim_directions(reference: &str, alt: &str) -> &'static [bool] {
    if reference.len() > 1 || alt.len() > 1 {
        &[false, true]
    } else {
        &[false]
    }
}

/// A reference-alternate pair as `trim_sequences` leaves it: the trimmed
/// reference, the trimmed alternate and the new start.
type TrimmedAlleles<'s> = (&'s [u8], &'s [u8], u64);

/// `trim_sequences` for a forward-strand pair (`Utils/Sequence.pm:965-1038`):
/// bytes shared at one end are removed while both alleles are Perl-true
/// (non-empty and not `0`), the left end first and then the right unless
/// `end_first`, the start moving up by one per byte trimmed from the left, as
/// Perl's `substr` moves it over an undecoded string; an emptied allele reads
/// `-`. The returned alleles are slices of the inputs.
fn trim_sequences<'s>(
    reference: &'s str,
    alt: &'s str,
    mut start: u64,
    end_first: bool,
) -> TrimmedAlleles<'s> {
    let mut trim_left = |mut reference: &'s [u8], mut alt: &'s [u8]| {
        while perl_true_bytes(reference) && perl_true_bytes(alt) && reference[0] == alt[0] {
            reference = &reference[1..];
            alt = &alt[1..];
            start += 1;
        }
        (reference, alt)
    };
    let trim_right = |mut reference: &'s [u8], mut alt: &'s [u8]| {
        while perl_true_bytes(reference)
            && perl_true_bytes(alt)
            && reference[reference.len() - 1] == alt[alt.len() - 1]
        {
            reference = &reference[..reference.len() - 1];
            alt = &alt[..alt.len() - 1];
        }
        (reference, alt)
    };
    let (reference, alt) = if end_first {
        let (reference, alt) = trim_right(reference.as_bytes(), alt.as_bytes());
        trim_left(reference, alt)
    } else {
        let (reference, alt) = trim_left(reference.as_bytes(), alt.as_bytes());
        trim_right(reference, alt)
    };
    let dash = |allele: &'s [u8]| if allele.is_empty() { b"-" } else { allele };
    (dash(reference), dash(alt), start)
}

/// Perl's truth of an allele while it is trimmed: the empty string and `0` are
/// false, so trimming stops at either.
fn perl_true_bytes(allele: &[u8]) -> bool {
    !allele.is_empty() && allele != b"0"
}

/// [`reverse_comp`] of `sequence` when `flip`, else `sequence` as written.
fn reverse_comp_if(flip: bool, sequence: &str) -> Cow<'_, str> {
    if flip {
        Cow::Owned(reverse_comp(sequence))
    } else {
        Cow::Borrowed(sequence)
    }
}

/// The `ColocatedVariant` of a matched record: the record's fields, its
/// frequency columns as the cache wrote them, and the matched alleles.
fn cached_to_colocated(
    cv: CachedVariation,
    matched_alleles: Vec<MatchedAllele>,
) -> ColocatedVariant {
    let comma_list = |value: Option<String>| -> Vec<String> {
        value
            .filter(|s| !s.is_empty())
            .map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
            .unwrap_or_default()
    };
    ColocatedVariant {
        chr: cv.chr,
        id: cv.variation_name,
        start: cv.start,
        end: cv.end,
        allele_string: Some(cv.allele_string).filter(|s| !s.is_empty()),
        strand: cv.strand,
        somatic: cv.somatic != 0,
        clin_sig: comma_list(cv.clin_sig),
        phenotype: cv.phenotype_or_disease != 0,
        pubmed: comma_list(cv.pubmed),
        matched_alleles,
        frequency_strings: cv.frequencies,
        clinvar: cv.clinvar,
        var_synonyms: cv.var_synonyms,
    }
}

/// Appends the `Extra` column of one default-format row: `KEY=VALUE` pairs
/// joined by `;`, keys in VEP's order (IMPACT, DISTANCE, STRAND, FLAGS, then the
/// active flag fields, then plugin keys alphabetically), values with `;` escaped
/// as `%3B` and their control bytes percent-encoded, lists joined by `,`, and a
/// lone `-` when no field applies. `tc`
/// is `None` for an intergenic row, which carries IMPACT and no STRAND, and for
/// a summary row, which carries the variant-level keys alone.
///
/// This is the hot path for the default `-o vep` format: it runs once per
/// (variant x transcript consequence). For a sequence variant with no optional
/// field group active, no plugin data and a vocabulary that prints DISTANCE as
/// the SO terms do, it writes the four fixed keys directly, and so does a
/// transcript or intergenic row when the `regulatory` group is the only other
/// one active, since its six fields are empty on such a row; every other row
/// writes each value through [`vep_io::output::fields::write_field_value`], and
/// a row without plugin data that is not structural carries the plan's flag
/// fields alone (a transcript row without the five motif fields), so its key
/// list is not rebuilt.
pub fn write_extra_fields(
    variant: &InputVariant,
    tc: Option<&vep_core::consequence::TranscriptConsequence>,
    plan: &vep_io::output::fields::ExtraFieldsPlan,
    out: &mut Vec<u8>,
) {
    let has_plugin_data =
        !variant.plugin_data.is_empty() || tc.is_some_and(|t| !t.plugin_data.is_empty());
    let summary_row = tc.is_none() && variant.summary_terms.is_some();
    let transcript_row = tc.is_none_or(|t| t.feature_type == FeatureType::Transcript);
    if (!plan.any_optional_group || (plan.transcript_fast_path && transcript_row))
        && !has_plugin_data
        && !variant.is_structural
        && !summary_row
        && plan.options.user_group_is_plain()
    {
        return write_user_group(tc, out);
    }

    let options = &plan.options;
    let mut first = true;
    if has_plugin_data || variant.is_structural {
        for key in vep_io::output::fields::extra_keys_with(&plan.flag_fields, variant, tc, options)
        {
            write_extra_pair(out, key, variant, tc, options, &mut first);
        }
    } else {
        let keys = if transcript_row {
            &plan.transcript_flag_fields
        } else {
            &plan.flag_fields
        };
        for key in keys {
            write_extra_pair(out, key, variant, tc, options, &mut first);
        }
    }
    if first {
        out.push(b'-');
    }
}

/// Appends `;KEY=VALUE` (`KEY=VALUE` for the first pair, clearing `first`) with
/// `;` in the value escaped as `%3B` and its control bytes percent-encoded
/// ([`vep_io::output::fields::percent_encode_control_bytes`]), or nothing when
/// the row has no value and the empty value is not one the format writes
/// ([`vep_io::output::fields::empty_value_is_written`]).
fn write_extra_pair(
    out: &mut Vec<u8>,
    key: &str,
    variant: &InputVariant,
    tc: Option<&vep_core::consequence::TranscriptConsequence>,
    options: &vep_io::output::fields::FieldOptions,
    first: &mut bool,
) {
    let pair_start = out.len();
    if !*first {
        out.push(b';');
    }
    out.extend_from_slice(key.as_bytes());
    out.push(b'=');
    let value_start = out.len();
    vep_io::output::fields::write_field_value(out, key, variant, tc, options, ",");
    if out.len() == value_start
        && !vep_io::output::fields::empty_value_is_written(key, variant, tc, options)
    {
        out.truncate(pair_start);
        return;
    }
    *first = false;
    if out[value_start..].contains(&b';') {
        let value = out.split_off(value_start);
        for &b in &value {
            if b == b';' {
                out.extend_from_slice(b"%3B");
            } else {
                out.push(b);
            }
        }
    }
    vep_io::output::fields::percent_encode_control_bytes(out, value_start);
}

/// The four unconditional keys, written as bytes. IMPACT is always present;
/// DISTANCE only on an up/downstream row; STRAND and FLAGS only against a
/// transcript, FLAGS only when the transcript carries a `cds_` attribute.
fn write_user_group(tc: Option<&vep_core::consequence::TranscriptConsequence>, out: &mut Vec<u8>) {
    let Some(tc) = tc else {
        return out.extend_from_slice(b"IMPACT=MODIFIER");
    };
    out.extend_from_slice(b"IMPACT=");
    out.extend_from_slice(tc.impact.as_str().as_bytes());
    if let Some(dist) = tc.distance {
        out.extend_from_slice(b";DISTANCE=");
        push_decimal(out, dist);
    }
    out.extend_from_slice(b";STRAND=");
    out.extend_from_slice(if tc.strand < 0 { b"-1" } else { b"1" });
    let mut wrote_flags = false;
    let flags_start = out.len();
    for flag in tc.flags.iter().filter(|f| f.starts_with("cds_")) {
        out.extend_from_slice(if wrote_flags { b"," } else { b";FLAGS=" });
        out.extend_from_slice(flag.as_bytes());
        wrote_flags = true;
    }
    vep_io::output::fields::percent_encode_control_bytes(out, flags_start);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cache record as `parse_record` builds it from the release-116 GRCh38
    /// chromosome 21 cache (`chr variation_name failed somatic start end allele_string
    /// strand`): the cited records are quoted from that file.
    fn record(name: &str, start: u64, end: u64, allele_string: &str) -> CachedVariation {
        CachedVariation {
            variation_name: name.to_string(),
            start,
            end,
            allele_string: allele_string.to_string(),
            strand: 1,
            ..Default::default()
        }
    }

    fn somatic(mut record: CachedVariation) -> CachedVariation {
        record.somatic = 1;
        record
    }

    /// The input variants of one VCF data line, as the run parses it.
    fn parsed(line: &str) -> Vec<InputVariant> {
        crate::vcf_parser::parse_vcf_line(line, false).expect("a valid VCF line")
    }

    /// The names that print on the first allele's row under `rules`.
    fn existing(
        line: &str,
        candidates: &[CachedVariation],
        rules: &ColocationRules,
    ) -> Vec<String> {
        let variant = &parsed(line)[0];
        let colocated = find_colocated_variants(variant, candidates.to_vec(), rules);
        existing_variation_names(&colocated, &variant.display_allele())
    }

    fn names(colocated: &[ColocatedVariant]) -> Vec<&str> {
        colocated.iter().map(|cv| cv.id.as_str()).collect()
    }

    const DEFAULT: ColocationRules = ColocationRules {
        no_check_alleles: false,
        exclude_null_alleles: false,
        failed: 0,
    };

    /// `21 14109044 G A` against the record set at 14109044 as the release-116
    /// cache holds it: an HGMD null-allele record, a COSMIC somatic null-allele
    /// record and `rs7278737` (`G/A/T`); Ensembl VEP 116.2 prints the three in the
    /// order dbSNP, HGMD, COSMIC, and with ALT `C` the dbSNP record falls away and
    /// the two null-allele records stay. The COSMIC and HGMD names here are
    /// stand-ins of the real records' prefix classes.
    #[test]
    fn snv_with_different_alt_does_not_match() {
        let candidates = [
            record("CM990001", 14109044, 14109044, "HGMD_MUTATION"),
            somatic(record("COSV9900001", 14109044, 14109044, "COSMIC_MUTATION")),
            record("rs7278737", 14109044, 14109044, "G/A/T"),
        ];
        assert_eq!(
            existing("21\t14109044\t.\tG\tA\t.\t.\t.", &candidates, &DEFAULT),
            ["rs7278737", "CM990001", "COSV9900001"]
        );
        assert_eq!(
            existing("21\t14109044\t.\tG\tC\t.\t.\t.", &candidates, &DEFAULT),
            ["CM990001", "COSV9900001"]
        );
        let variant = &parsed("21\t14109044\t.\tG\tA\t.\t.\t.")[0];
        let colocated = find_colocated_variants(variant, candidates.to_vec(), &DEFAULT);
        assert_eq!(
            colocated[2].matched_alleles,
            [MatchedAllele {
                input_allele: "A".to_string(),
                input_index: 0,
                record_allele: "A".to_string(),
                record_index: 0,
            }]
        );
        assert!(
            colocated[0].matched_alleles.is_empty(),
            "a span match carries no alleles"
        );
    }

    /// The cache keeps `rs2123010281` untrimmed as `GG/G` at 5029945-5029946.
    /// `21 5029945 GG G` minimises to `G/-` at 5029946 and matches it by trimming the
    /// record from the left; `21 5029944 NG N` is `G/-` at 5029945 and matches it by
    /// trimming from the right, the second order tried. Ensembl VEP 116.2 prints
    /// `rs2123010281` for both. With `--no_check_alleles` neither matches: no record
    /// spans exactly 5029945-5029945 or 5029946-5029946, and Ensembl VEP prints `-`.
    #[test]
    fn indel_matches_after_trimming_both_directions() {
        let candidates = [record("rs2123010281", 5029945, 5029946, "GG/G")];
        assert_eq!(
            existing("21\t5029945\t.\tGG\tG\t.\t.\t.", &candidates, &DEFAULT),
            ["rs2123010281"]
        );
        assert_eq!(
            existing("21\t5029944\t.\tNG\tN\t.\t.\t.", &candidates, &DEFAULT),
            ["rs2123010281"]
        );
        let no_check = ColocationRules {
            no_check_alleles: true,
            ..DEFAULT
        };
        assert!(existing("21\t5029945\t.\tGG\tG\t.\t.\t.", &candidates, &no_check).is_empty());
        assert!(existing("21\t5029944\t.\tNG\tN\t.\t.\t.", &candidates, &no_check).is_empty());
    }

    /// `rs71668025` is `TT/T/TTT` at 28382673-28382674. `21 28382673 TT T,TTTT` is
    /// the record `T/-/TTT` at 28382674: the deletion matches the record's `T`, the
    /// three-base insertion matches nothing (the record's `TTT` is a one-base
    /// insertion), so Ensembl VEP 116.2 prints `rs71668025` on the `-` row and `-`
    /// on the `TTT` row.
    #[test]
    fn insertion_length_decides_the_match_of_a_multi_allelic_record() {
        let candidates = [
            record("rs1601149035", 28382673, 28382673, "T/A/G"),
            record("rs71668025", 28382673, 28382674, "TT/T/TTT"),
        ];
        let alleles = parsed("21\t28382673\t.\tTT\tT,TTTT\t.\t.\t.");
        let colocated = find_colocated_variants(&alleles[0], candidates.to_vec(), &DEFAULT);
        assert_eq!(names(&colocated), ["rs71668025"]);
        assert_eq!(
            colocated[0].matched_alleles,
            [MatchedAllele {
                input_allele: "-".to_string(),
                input_index: 0,
                record_allele: "T".to_string(),
                record_index: 0,
            }]
        );
        assert!(colocated[0].applies_to_allele("-"));
        assert!(!colocated[0].applies_to_allele("TTT"));
    }

    /// The records at 38383755 are two COSMIC somatic null alleles then
    /// `rs2516998629` (`T/C`, somatic) in file order; Ensembl VEP 116.2 prints the
    /// dbSNP name first, then the two COSMIC names in file order, for
    /// `21 38383755 T C`. At 14109044 the non-somatic HGMD record and `rs7278737`
    /// precede the somatic COSMIC record, and the dbSNP name precedes the HGMD
    /// name. The COSMIC and HGMD names are stand-ins of the real records' prefix
    /// classes (`OutputFactory.pm:1057-1065` ranks by the first two letters).
    #[test]
    fn existing_variation_sorted_non_somatic_then_prefix_rank() {
        let candidates = [
            somatic(record("COSV9900002", 38383755, 38383755, "COSMIC_MUTATION")),
            somatic(record("COSV9900003", 38383755, 38383755, "COSMIC_MUTATION")),
            somatic(record("rs2516998629", 38383755, 38383755, "T/C")),
        ];
        assert_eq!(
            existing("21\t38383755\t.\tT\tC\t.\t.\t.", &candidates, &DEFAULT),
            ["rs2516998629", "COSV9900002", "COSV9900003"]
        );
        let mixed = [
            record("CM990001", 14109044, 14109044, "HGMD_MUTATION"),
            somatic(record("COSV9900001", 14109044, 14109044, "COSMIC_MUTATION")),
            record("rs7278737", 14109044, 14109044, "G/A/T"),
            record("CD990002", 14109044, 14109044, "HGMD_MUTATION"),
            record("TMP_1", 14109044, 14109044, "G/A"),
        ];
        assert_eq!(
            existing("21\t14109044\t.\tG\tA\t.\t.\t.", &mixed, &DEFAULT),
            ["rs7278737", "CM990001", "CD990002", "TMP_1", "COSV9900001"]
        );
        let variant = &parsed("21\t14109044\t.\tG\tA\t.\t.\t.")[0];
        let colocated = find_colocated_variants(variant, mixed.to_vec(), &DEFAULT);
        assert_eq!(
            names(&colocated),
            ["CM990001", "COSV9900001", "rs7278737", "CD990002", "TMP_1"],
            "the record keeps the cache's order"
        );
    }

    /// No record of the release-116 chromosome 21 cache is on the reverse strand, so
    /// the record is hand-built: `T/C` on strand -1 at 100 is `A/G` on the forward
    /// strand, which `get_matched_variant_alleles` finds by reverse-complementing
    /// the input's alleles (`Utils/Sequence.pm:1182, 1196`); the match names the
    /// alleles as each side spells them. `A/C` does not match it.
    #[test]
    fn reverse_strand_record_matches_complemented_alleles() {
        let mut reverse = record("rsREV", 100, 100, "T/C");
        reverse.strand = -1;
        let candidates = [reverse];
        let variant = &parsed("21\t100\t.\tA\tG\t.\t.\t.")[0];
        let colocated = find_colocated_variants(variant, candidates.to_vec(), &DEFAULT);
        assert_eq!(names(&colocated), ["rsREV"]);
        assert_eq!(
            colocated[0].matched_alleles,
            [MatchedAllele {
                input_allele: "G".to_string(),
                input_index: 0,
                record_allele: "C".to_string(),
                record_index: 0,
            }]
        );
        assert!(existing("21\t100\t.\tA\tC\t.\t.\t.", &candidates, &DEFAULT).is_empty());

        let mut reverse_deletion = record("rsREVDEL", 200, 201, "CA/-");
        reverse_deletion.strand = -1;
        assert_eq!(
            existing("21\t199\t.\tNTG\tN\t.\t.\t.", &[reverse_deletion], &DEFAULT),
            ["rsREVDEL"]
        );
    }

    /// An input on the reverse strand (the Ensembl input format's strand column)
    /// is matched through its reverse complement: `A/G` on strand -1 at 100 is
    /// `T/C` on the forward strand, so it matches a forward `T/C` record and not a
    /// forward `A/G` one, and the match names the input allele as the input spells
    /// it (`Utils/Sequence.pm:1182, 1196`: `a`'s alleles are complemented, `b`'s
    /// are not). A reverse-strand `TT/T` at 200-201 is a forward `AA/A` and matches
    /// the forward `A/-` at 201 by trimming.
    #[test]
    fn reverse_strand_input_matches_a_forward_record() {
        let mut reverse = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        reverse.strand = Strand::Reverse;
        let candidates = [
            record("rsFWD_TC", 100, 100, "T/C"),
            record("rsFWD_AG", 100, 100, "A/G"),
        ];
        let colocated = find_colocated_variants(&reverse, candidates.to_vec(), &DEFAULT);
        assert_eq!(names(&colocated), ["rsFWD_TC"]);
        assert_eq!(
            colocated[0].matched_alleles,
            [MatchedAllele {
                input_allele: "G".to_string(),
                input_index: 0,
                record_allele: "C".to_string(),
                record_index: 0,
            }]
        );
        assert!(colocated[0].applies_to_allele(&reverse.display_allele()));

        let mut deletion = InputVariant::new("21".into(), 200, 201, b"TT".to_vec(), b"T".to_vec());
        deletion.strand = Strand::Reverse;
        let candidates = [record("rsFWD_DEL", 201, 201, "A/-")];
        let colocated = find_colocated_variants(&deletion, candidates.to_vec(), &DEFAULT);
        assert_eq!(names(&colocated), ["rsFWD_DEL"]);
    }

    /// A structural variant is never co-located (`VariationTabix.pm:118` keeps the
    /// sequence variants alone), whatever record starts where it does.
    #[test]
    fn structural_variants_are_not_colocated() {
        let mut variant =
            parsed("21\t5029945\t.\tG\t<DEL>\t.\t.\tSVTYPE=DEL;END=5029950").remove(0);
        assert!(variant.is_structural);
        let dir = tempfile::tempdir().unwrap();
        let source = region_source(
            dir.path(),
            &["rs2123010281\t0\t0\t5029945\t5029946\tGG/G\t1\t\t0\t\t\t\t\t"],
        );
        let failure = ReadFailure::default();
        annotate_colocated(&mut variant, &source, &DEFAULT, &failure);
        assert!(variant.colocated_variants.is_empty());
        assert!(variant.existing_variation.is_empty());
        assert!(failure.take().is_none());

        variant = parsed("21\t5029945\t.\tGG\tG\t.\t.\t.").remove(0);
        annotate_colocated(&mut variant, &source, &DEFAULT, &failure);
        assert_eq!(variant.existing_variation, ["rs2123010281"]);
    }

    /// `21 14109044 G A,C`: `rs7278737` (`G/A/T`) matches the record through `A`
    /// alone, so Ensembl VEP 116.2 prints it on the `A` row and not on the `C` row,
    /// while the two null-allele records (HGMD and COSMIC stand-ins) print on both.
    /// Every allele carries the record's whole match list.
    #[test]
    fn rsid_not_matching_this_allele_is_omitted_from_its_row() {
        let dir = tempfile::tempdir().unwrap();
        let source = region_source(
            dir.path(),
            &[
                "CM990001\t0\t0\t14109044\t\tHGMD_MUTATION\t1\t\t1\t\t\t\t\t",
                "COSV9900001\t0\t1\t14109044\t\tCOSMIC_MUTATION\t1\t\t1\t\t\t\t\t",
                "rs7278737\t0\t0\t14109044\t\tG/A/T\t1\tbenign\t1\tT:benign\tG\t\t\t",
            ],
        );
        let failure = ReadFailure::default();
        let mut alleles = parsed("21\t14109044\t.\tG\tA,C\t.\t.\t.");
        for variant in &mut alleles {
            annotate_colocated(variant, &source, &DEFAULT, &failure);
        }
        assert_eq!(alleles[0].display_allele(), "A");
        assert_eq!(
            alleles[0].existing_variation,
            ["rs7278737", "CM990001", "COSV9900001"]
        );
        assert_eq!(alleles[1].display_allele(), "C");
        assert_eq!(alleles[1].existing_variation, ["CM990001", "COSV9900001"]);
        for variant in &alleles {
            assert_eq!(
                names(&variant.colocated_variants),
                ["CM990001", "COSV9900001", "rs7278737"]
            );
        }
    }

    /// Under `--no_check_alleles` a record matches on its start and end alone and
    /// carries no matched alleles: `21 38383755 T G` gains `rs2516998629` (`T/C`),
    /// `21 40000000 A G` gains `rs775128082` (`T/C`), as Ensembl VEP 116.2 prints;
    /// a record spanning one base more does not match.
    #[test]
    fn no_check_alleles_matches_on_start_and_end_only() {
        let rules = ColocationRules {
            no_check_alleles: true,
            ..DEFAULT
        };
        let candidates = [
            record("rs775128082", 40000000, 40000000, "T/C"),
            record("rsWIDER", 40000000, 40000001, "TA/C"),
        ];
        let variant = &parsed("21\t40000000\t.\tA\tG\t.\t.\t.")[0];
        let colocated = find_colocated_variants(variant, candidates.to_vec(), &rules);
        assert_eq!(names(&colocated), ["rs775128082"]);
        assert!(colocated[0].matched_alleles.is_empty());
        assert!(colocated[0].applies_to_allele("G"));
        assert!(existing("21\t40000000\t.\tA\tG\t.\t.\t.", &candidates, &DEFAULT).is_empty());
    }

    /// A COSMIC stand-in (`COSMIC_MUTATION`, no `/`) at 14109044 matches
    /// `21 14109044 G C` on its span alone, prints on every allele's row, and never
    /// matches under `--exclude_null_alleles`; a null-allele record at another span
    /// never matches, `--no_check_alleles` notwithstanding
    /// (`AnnotationType/Variation.pm:150-157`).
    #[test]
    fn null_allele_record_matches_on_span_unless_excluded() {
        let candidates = [
            somatic(record("COSV9900001", 14109044, 14109044, "COSMIC_MUTATION")),
            somatic(record("COSVWIDE", 14109044, 14109045, "COSMIC_MUTATION")),
        ];
        let line = "21\t14109044\t.\tG\tC\t.\t.\t.";
        assert_eq!(existing(line, &candidates, &DEFAULT), ["COSV9900001"]);
        let exclude = ColocationRules {
            exclude_null_alleles: true,
            ..DEFAULT
        };
        assert!(existing(line, &candidates, &exclude).is_empty());
        let no_check = ColocationRules {
            no_check_alleles: true,
            ..DEFAULT
        };
        assert_eq!(existing(line, &candidates, &no_check), ["COSV9900001"]);
        let both = ColocationRules {
            no_check_alleles: true,
            exclude_null_alleles: true,
            ..DEFAULT
        };
        assert!(existing(line, &candidates, &both).is_empty());
    }

    /// `rs1979841086` at 5230609 is `dbSNP_novariation` with `failed` 1 beside
    /// `rs373402393` (`A/C`): Ensembl VEP 116.2 prints `rs373402393` for
    /// `21 5230609 A C` and `rs1979841086,rs373402393` under `--failed 1`.
    #[test]
    fn failed_threshold_admits_failed_records() {
        let mut failed = record("rs1979841086", 5230609, 5230609, "dbSNP_novariation");
        failed.failed = 1;
        let candidates = [failed, record("rs373402393", 5230609, 5230609, "A/C")];
        let line = "21\t5230609\t.\tA\tC\t.\t.\t.";
        assert_eq!(existing(line, &candidates, &DEFAULT), ["rs373402393"]);
        let admit = ColocationRules {
            failed: 1,
            ..DEFAULT
        };
        assert_eq!(
            existing(line, &candidates, &admit),
            ["rs1979841086", "rs373402393"]
        );
    }

    /// A name already kept is not matched again (`VariationTabix.pm:237`): the first
    /// record of the name stays.
    #[test]
    fn a_name_is_kept_once() {
        let candidates = [
            record("rs1", 100, 100, "A/G"),
            record("rs1", 100, 100, "A/G/T"),
            record("rs2", 100, 100, "A/G"),
        ];
        let variant = &parsed("21\t100\t.\tA\tG\t.\t.\t.")[0];
        let colocated = find_colocated_variants(variant, candidates.to_vec(), &DEFAULT);
        assert_eq!(names(&colocated), ["rs1", "rs2"]);
        assert_eq!(colocated[0].allele_string.as_deref(), Some("A/G"));
    }

    /// A match carries the record's frequency columns as the cache wrote them
    /// and the alleles matched, from which the renderer reads each row's
    /// allele: the record `A/G/T` with `AF T:0.1,G:0.3` matches the `G` and `T`
    /// inputs to their own record alleles; a record matched on its span alone
    /// carries the columns and no matched allele.
    #[test]
    fn matches_carry_the_frequency_columns_and_the_matched_alleles() {
        let mut rec = record("rs1", 100, 100, "A/G/T");
        rec.frequencies
            .insert("AF".to_string(), "T:0.1,G:0.3".to_string());
        let candidates = [rec];
        let alleles = parsed("21\t100\t.\tA\tG,T\t.\t.\t.");
        let for_g = find_colocated_variants(&alleles[0], candidates.to_vec(), &DEFAULT);
        assert_eq!(for_g[0].record_allele_for("G"), Some("G"));
        assert_eq!(for_g[0].record_allele_for("T"), Some("T"));
        let for_t = find_colocated_variants(&alleles[1], candidates.to_vec(), &DEFAULT);
        assert_eq!(for_t[0].matched_alleles, for_g[0].matched_alleles);
        assert_eq!(
            for_t[0].frequency_strings.get("AF").map(String::as_str),
            Some("T:0.1,G:0.3")
        );
        let span_only = find_colocated_variants(
            &parsed("21\t100\t.\tA\tC\t.\t.\t.")[0],
            candidates.to_vec(),
            &ColocationRules {
                no_check_alleles: true,
                ..DEFAULT
            },
        );
        assert!(span_only[0].matched_alleles.is_empty());
        assert_eq!(
            span_only[0].frequency_strings.get("AF").map(String::as_str),
            Some("T:0.1,G:0.3")
        );
    }

    /// `trim_sequences` as `Utils/Sequence.pm:965-1038` trims: `AT/ATT` from the
    /// left is `-/T` at 102, from the right first `-/T` at 101; `GG/G` from the
    /// left is `G/-` at 101, from the right first `G/-` at 100; a one-base pair is
    /// untouched; identical alleles empty to `-/-`.
    #[test]
    fn trim_sequences_follows_perl_in_both_orders() {
        let trimmed = |r: &str, a: &str, start: u64, end_first: bool| {
            let (r, a, pos) = trim_sequences(r, a, start, end_first);
            (
                String::from_utf8_lossy(r).into_owned(),
                String::from_utf8_lossy(a).into_owned(),
                pos,
            )
        };
        let t = |r: &str, a: &str, pos: u64| (r.to_string(), a.to_string(), pos);
        assert_eq!(trimmed("AT", "ATT", 100, false), t("-", "T", 102));
        assert_eq!(trimmed("AT", "ATT", 100, true), t("-", "T", 101));
        assert_eq!(trimmed("GG", "G", 100, false), t("G", "-", 101));
        assert_eq!(trimmed("GG", "G", 100, true), t("G", "-", 100));
        assert_eq!(trimmed("A", "G", 100, false), t("A", "G", 100));
        assert_eq!(trimmed("AC", "AC", 100, false), t("-", "-", 102));
        assert_eq!(trimmed("0", "0", 100, false), t("0", "0", 100));
        assert_eq!(trimmed("A0", "A0", 100, false), t("0", "0", 101));
        // Perl trims the bytes of an undecoded string: a two-byte character shared
        // at the left moves the start by two.
        assert_eq!(trimmed("\u{e9}A", "\u{e9}C", 100, false), t("A", "C", 102));
        assert_eq!(trim_directions("A", "G"), &[false]);
        assert_eq!(trim_directions("AT", "A"), &[false, true]);
    }

    /// `get_matched_variant_alleles` on the shapes `Utils/Sequence.pm:1151-1173`
    /// refuses: a position of 0, an allele string opening with `/`, one without
    /// alternate alleles; and the hash's last-wins rule for a repeated allele.
    #[test]
    fn unusable_allele_sets_match_nothing() {
        let a = AlleleSet {
            allele_string: "A/G",
            pos: 100,
            strand: 1,
        };
        let b = AlleleSet {
            allele_string: "A/G",
            pos: 0,
            strand: 1,
        };
        assert!(get_matched_variant_alleles(&a, &b).is_empty());
        let slash = AlleleSet {
            allele_string: "/G",
            pos: 100,
            strand: 1,
        };
        assert!(get_matched_variant_alleles(&a, &slash).is_empty());
        let null = AlleleSet {
            allele_string: "COSMIC_MUTATION",
            pos: 100,
            strand: 1,
        };
        assert!(get_matched_variant_alleles(&a, &null).is_empty());
        let repeated = AlleleSet {
            allele_string: "A/G/G",
            pos: 100,
            strand: 1,
        };
        let matched = get_matched_variant_alleles(&repeated, &a);
        assert_eq!(matched.len(), 1);
        assert_eq!(
            matched[0].input_index, 1,
            "the later allele owns the shared key"
        );
    }

    /// Every allele of the three committed co-location corpora against the slice
    /// committed beside it: `Existing_variation` as Ensembl VEP 116.2 wrote it for
    /// each `Uploaded_variation` and `Allele` under the flags named here, `-` for
    /// none. The two `colocated` corpora are matched under `--failed 1`; among their
    /// rows the record the window misses (`21_5030482_GGGAGGGGTGG/GG` trims to
    /// 5030484-5030492 and `rs1987113057` starts at 5030482), the failed null allele
    /// `rs1979841086` admitted beside `rs373402393`, the multi-allelic records
    /// whose every allele matches one record, and the four multi-allelic records
    /// matched whole without `--minimal`, each allele finding its own record.
    #[test]
    fn committed_corpus_slices_match_ensembl_vep_existing_variation() {
        const GRCH38: &[(&str, &str, &str)] = &[
            ("21_5029945_GG/G", "-", "rs2123010281"),
            ("21_5030192_A/G", "G", "rs1172282560"),
            ("21_5030287_A/T", "T", "rs1249916603"),
            ("21_5030482_GGGAGGGGTGG/GG", "-", "-"),
            ("21_5030602_A/AA", "A", "rs1370167485"),
            ("21_5030602_A/AA/G", "AA", "rs1370167485"),
            ("21_5030602_A/AA/G", "G", "-"),
            ("21_5030602_A/C", "C", "rs1300916113"),
            ("21_5030602_A/G/AA", "G", "-"),
            ("21_5030602_A/G/AA", "AA", "rs1370167485"),
            ("21_5030625_GCAT/G", "-", "rs1380012774"),
            ("21_5030633_C/T", "T", "rs2123011366"),
            ("21_5030664_A/T", "T", "rs1987116015"),
            ("21_5031926_T/A", "A", "rs1435841083"),
            ("21_5032970_T/TCAGTGAG", "CAGTGAG", "rs1568756826"),
            ("21_5033370_C/CA/CT", "A", "rs1286929552"),
            ("21_5033370_C/CA/CT", "T", "rs1286929552"),
            ("21_5131441_T/C", "C", "rs1555830109"),
            ("21_5131480_A/G", "G", "rs1555830111"),
            ("21_5161864_G/A", "A", "-"),
            ("21_5230609_A/AT/C", "AT", "rs1979841086,rs2123086048"),
            ("21_5230609_A/AT/C", "C", "rs1979841086,rs373402393"),
            ("21_5230609_A/C", "C", "rs1979841086,rs373402393"),
            ("21_5230609_A/C/AT", "C", "rs1979841086,rs373402393"),
            ("21_5230609_A/C/AT", "AT", "rs1979841086,rs2123086048"),
            ("21_6138391_G/C", "C", "rs1555842629"),
            ("21_7759128_T/A/G", "A", "rs1555870136"),
            ("21_7759128_T/A/G", "G", "rs1555870136"),
            ("21_9040338_T/A/C", "A", "rs75275574"),
            ("21_9040338_T/A/C", "C", "rs75275574"),
            ("21_10503115_A/C", "C", "-"),
            ("21_10541121_G/A/T", "A", "rs374325784"),
            ("21_10541121_G/A/T", "T", "rs374325784"),
            ("21_10602054_ATA/A", "-", "rs149363218"),
            ("21_13610458_C/A/G/T", "A", "rs201981080"),
            ("21_13610458_C/A/G/T", "G", "rs201981080"),
            ("21_13610458_C/A/G/T", "T", "rs201981080"),
            ("21_14967531_T/C", "C", "rs139263261"),
            ("21_14981488_A/C/G/T", "C", "rs17274750"),
            ("21_14981488_A/C/G/T", "G", "rs17274750"),
            ("21_14981488_A/C/G/T", "T", "rs17274750"),
            ("21_38403526_T/TTCC", "TCC", "rs2517045861"),
            ("21_44195140_G/A", "A", "rs2838519"),
            ("21_44195140_G/A/C", "A", "rs2838519"),
            ("21_44195140_G/A/C", "C", "rs2838519"),
            ("21_44293150_TGTGTGTG/TGTGTG", "-", "-"),
        ];
        const GRCH37: &[(&str, &str, &str)] = &[
            ("21_9411284_C/A", "A", "rs1980632183"),
            ("21_9411349_T/A", "A", "rs1455295117"),
            ("21_9411379_CT/C", "-", "rs1432311608"),
            ("21_9411405_TCTATC/TC", "-", "-"),
            ("21_9411410_C/A", "A", "-"),
            ("21_9411555_C/T", "T", "rs1980640962"),
            ("21_9412358_T/TG", "G", "rs555324129"),
            ("21_9412370_A/AAC/AT", "AC", "rs1555879625"),
            ("21_9412370_A/AAC/AT", "T", "rs1555879625"),
            ("21_9412503_C/A/G/T", "A", "rs71220884"),
            ("21_9412503_C/A/G/T", "G", "rs71220884"),
            ("21_9412503_C/A/G/T", "T", "rs71220884"),
            ("21_9413584_G/A", "A", "rs80143362"),
            ("21_9492382_C/A/T", "A", "rs1326732194"),
            ("21_9492382_C/A/T", "T", "rs1326732194"),
            ("21_9569350_C/T", "T", "rs538160740"),
            ("21_9569836_T/C", "C", "rs552769896"),
            ("21_9879171_T/A/C", "A", "rs75275574"),
            ("21_9879171_T/A/C", "C", "rs75275574"),
            ("21_10090673_T/A", "A", "-"),
            ("21_10910401_TAT/T", "-", "rs149363218"),
            ("21_15588513_C/A/G/T", "A", "rs753372738"),
            ("21_15588513_C/A/G/T", "G", "rs753372738"),
            ("21_15588513_C/A/G/T", "T", "rs753372738"),
            ("21_16340235_CC/C", "-", "rs1555879360"),
            ("21_16438793_A/G/T", "G", "rs2223043"),
            ("21_16438793_A/G/T", "T", "rs2223043"),
            ("21_43005598_G/A", "A", "rs151008585,rs2146578403"),
        ];
        // GRCh38-colocated-alt under `--no_check_alleles --exclude_null_alleles
        // --failed 1`: a record matches on its span alone, so `21_5030602_A/G` gains
        // both records at 5030602 and `21_5030602_A/AA`, whose minimised span is
        // 5030603-5030602, gains none; the failed null allele `rs1979841086` is
        // excluded while the failed `rs1555857356` (`A/ANA`) is admitted on the span
        // of `21_13492486_A/G`.
        const ALT_SPAN_ONLY: &[(&str, &str, &str)] = &[
            ("21_5030241_CCCCC/CCCC", "-", "-"),
            ("21_5030602_A/AA", "A", "-"),
            ("21_5030602_A/G", "G", "rs1300916113,rs1370167485"),
            ("21_5031987_G/A/T", "A", "rs1190341920"),
            ("21_5031987_G/A/T", "T", "rs1190341920"),
            ("21_5031990_T/C", "C", "rs1601770516"),
            ("21_5063607_A/AA", "A", "rs2123035150"),
            ("21_5227650_CGC/C", "-", "rs1555831393"),
            ("21_5230609_A/C", "C", "rs373402393"),
            ("21_7819485_G/A/C", "A", "rs1983989595"),
            ("21_7819485_G/A/C", "C", "rs1983989595"),
            ("21_8680517_C/A", "A", "rs538160740"),
            ("21_8681003_T/A", "A", "rs552769896"),
            ("21_8701352_A/G", "G", "rs528851961"),
            ("21_8701459_A/T", "T", "rs553814603"),
            ("21_10541121_G/A", "A", "rs374325784"),
            ("21_10541121_G/A/C/T", "A", "rs374325784"),
            ("21_10541121_G/A/C/T", "C", "rs374325784"),
            ("21_10541121_G/A/C/T", "T", "rs374325784"),
            ("21_10541137_G/A", "A", "rs150733676"),
            ("21_10541137_G/A/C/T", "A", "rs150733676"),
            ("21_10541137_G/A/C/T", "C", "rs150733676"),
            ("21_10541137_G/A/C/T", "T", "rs150733676"),
            ("21_12002666_T/A", "A", "-"),
            ("21_13492486_A/ANA", "NA", "-"),
            ("21_13492486_A/G", "G", "rs1555857356"),
            ("21_14109088_AGAGA/AGA", "-", "-"),
            ("21_18281193_CT/C", "-", "rs2074694480"),
            ("21_20257301_G/C", "C", "rs2123213965"),
            ("21_30000028_CC/C", "-", "rs2072291126"),
            ("21_33248813_TCT/TCTC", "C", "-"),
            ("21_41585438_G/A", "A", "rs2146578403"),
        ];
        // The same corpus under `--check_existing --failed 1`: the alleles decide again.
        const ALT_BY_ALLELE: &[(&str, &str, &str)] = &[
            ("21_5030241_CCCCC/CCCC", "-", "-"),
            ("21_5030602_A/AA", "A", "rs1370167485"),
            ("21_5030602_A/G", "G", "-"),
            ("21_5031987_G/A/T", "A", "rs1190341920"),
            ("21_5031987_G/A/T", "T", "rs1190341920"),
            ("21_5031990_T/C", "C", "rs1601770516"),
            ("21_5063607_A/AA", "A", "-"),
            ("21_5227650_CGC/C", "-", "rs1601785900,rs1568769193"),
            ("21_5230609_A/C", "C", "rs1979841086,rs373402393"),
            ("21_7819485_G/A/C", "A", "rs1983989595"),
            ("21_7819485_G/A/C", "C", "-"),
            ("21_8680517_C/A", "A", "-"),
            ("21_8681003_T/A", "A", "-"),
            ("21_8701352_A/G", "G", "rs528851961"),
            ("21_8701459_A/T", "T", "rs553814603"),
            ("21_10541121_G/A", "A", "rs374325784"),
            ("21_10541121_G/A/C/T", "A", "rs374325784"),
            ("21_10541121_G/A/C/T", "C", "rs374325784"),
            ("21_10541121_G/A/C/T", "T", "rs374325784"),
            ("21_10541137_G/A", "A", "rs150733676"),
            ("21_10541137_G/A/C/T", "A", "rs150733676"),
            ("21_10541137_G/A/C/T", "C", "rs150733676"),
            ("21_10541137_G/A/C/T", "T", "rs150733676"),
            ("21_12002666_T/A", "A", "-"),
            ("21_13492486_A/ANA", "NA", "rs1555857356"),
            ("21_13492486_A/G", "G", "-"),
            ("21_14109088_AGAGA/AGA", "-", "-"),
            ("21_18281193_CT/C", "-", "rs2074694480"),
            ("21_20257301_G/C", "C", "-"),
            ("21_30000028_CC/C", "-", "-"),
            ("21_33248813_TCT/TCTC", "C", "-"),
            ("21_41585438_G/A", "A", "rs151008585,rs2146578403"),
        ];
        let admit_failed = ColocationRules {
            failed: 1,
            ..DEFAULT
        };
        let span_only = ColocationRules {
            no_check_alleles: true,
            exclude_null_alleles: true,
            failed: 1,
        };
        for (corpus, rules, expected) in [
            ("GRCh38-colocated", &admit_failed, GRCH38),
            ("GRCh37-colocated", &admit_failed, GRCH37),
            ("GRCh38-colocated-alt", &span_only, ALT_SPAN_ONLY),
            ("GRCh38-colocated-alt", &admit_failed, ALT_BY_ALLELE),
        ] {
            let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../tests/golden/116")
                .join(corpus);
            let source =
                VariationSource::open_json_cache(dir.join("json_cache").to_str().unwrap(), false)
                    .unwrap()
                    .expect("the corpus carries a variation slice");
            let failure = ReadFailure::default();
            let mut got = Vec::new();
            let vcf = std::fs::read_to_string(dir.join("variants.vcf")).unwrap();
            for line in vcf.lines().filter(|l| !l.starts_with('#')) {
                for mut variant in parsed(line) {
                    annotate_colocated(&mut variant, &source, rules, &failure);
                    let names = if variant.existing_variation.is_empty() {
                        "-".to_string()
                    } else {
                        variant.existing_variation.join(",")
                    };
                    got.push((
                        variant.uploaded_variation(),
                        variant.display_allele(),
                        names,
                    ));
                }
            }
            assert!(failure.take().is_none());
            let expected: Vec<(String, String, String)> = expected
                .iter()
                .map(|(name, allele, names)| {
                    (name.to_string(), allele.to_string(), names.to_string())
                })
                .collect();
            assert_eq!(got, expected, "{corpus} under {rules:?}");
        }
    }

    /// Under `--minimal` every row of a split record prints the co-located
    /// variants of the record's first feature, each row keeping the names
    /// matched on its own allele, as Ensembl VEP 116.2 prints on the GRCh38 cache
    /// under `--check_existing --minimal --failed 1` (`split_variants` runs before
    /// the sources, InputBuffer.pm:265, 658-676; the rejoin leaves `existing` on
    /// each feature and the rows read the first feature's, OutputFactory.pm:
    /// 2272-2353, 1046, 1081-1083):
    ///
    /// * `21 5030602 A AA,G`: the first feature `-/A` at 5030603 matches
    ///   rs1370167485 (`A/AA` at 5030602); the A row prints it, the G row `-`.
    /// * `21 5030602 A G,AA`: the first feature `A/G` matches nothing (the
    ///   cache's `A/C` and `A/AA` at 5030602); both rows print `-`, the AA
    ///   allele's own match included.
    /// * `21 5230609 A C,AT`: the first feature `A/C` at 5230609 matches
    ///   rs1979841086 (a failed null allele, on its span) and rs373402393; the
    ///   C row prints both, the T row the span match alone.
    /// * `21 5230609 A AT,C`: the first feature `-/T` at 5230610 matches
    ///   rs2123086048; the T row prints it, the C row `-`.
    #[test]
    fn a_minimal_split_record_prints_its_first_features_matches() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh38-colocated");
        let source =
            VariationSource::open_json_cache(dir.join("json_cache").to_str().unwrap(), false)
                .unwrap()
                .expect("the corpus carries a variation slice");
        let rules = ColocationRules {
            failed: 1,
            ..ColocationRules::default()
        };
        let failure = ReadFailure::default();
        let expected = [
            (
                "21\t5030602\t.\tA\tAA,G\t.\t.\t.",
                vec![("A", "rs1370167485"), ("G", "-")],
            ),
            (
                "21\t5030602\t.\tA\tG,AA\t.\t.\t.",
                vec![("G", "-"), ("A", "-")],
            ),
            (
                "21\t5230609\t.\tA\tC,AT\t.\t.\t.",
                vec![("C", "rs1979841086,rs373402393"), ("T", "rs1979841086")],
            ),
            (
                "21\t5230609\t.\tA\tAT,C\t.\t.\t.",
                vec![("T", "rs2123086048"), ("C", "-")],
            ),
        ];
        for (record, (line, rows)) in expected.into_iter().enumerate() {
            let parsed = crate::vcf_parser::parse_vcf_line(line, true).expect("a valid VCF line");
            let mut alleles = crate::vcf_parser::split_minimal_record(parsed);
            for allele in &mut alleles {
                allele.input_record = record as u64 + 1;
                assert!(allele.original_start.is_some(), "{line}: the record splits");
                annotate_colocated(allele, &source, &rules, &failure);
            }
            crate::runner::apply_record_rules(&mut alleles, None);
            let got: Vec<(String, String)> = alleles
                .iter()
                .map(|allele| {
                    let names = if allele.existing_variation.is_empty() {
                        "-".to_string()
                    } else {
                        allele.existing_variation.join(",")
                    };
                    (allele.display_allele(), names)
                })
                .collect();
            let rows: Vec<(String, String)> = rows
                .iter()
                .map(|(allele, names)| (allele.to_string(), names.to_string()))
                .collect();
            assert_eq!(got, rows, "{line}");
            let ids = |allele: &InputVariant| -> Vec<String> {
                allele
                    .colocated_variants
                    .iter()
                    .map(|cv| cv.id.clone())
                    .collect()
            };
            assert!(
                alleles
                    .windows(2)
                    .all(|pair| ids(&pair[0]) == ids(&pair[1])),
                "{line}: every allele carries the first feature's matches"
            );
        }
        assert!(failure.take().is_none());
    }

    /// An empty input allele (`A//G` splits to `A`, an empty string and `G`) is
    /// stored under its trimmed form but never matched: Perl tests the allele it
    /// finds for truth (`Utils/Sequence.pm:1232`), and the empty string is false, so
    /// a record `A/-`, whose own form is the same `A`, `-`, matches nothing of it.
    #[test]
    fn an_empty_input_allele_matches_nothing() {
        let a = AlleleSet {
            allele_string: "A//G",
            pos: 100,
            strand: 1,
        };
        let b = AlleleSet {
            allele_string: "A/-",
            pos: 100,
            strand: 1,
        };
        assert!(get_matched_variant_alleles(&a, &b).is_empty());
        let g = AlleleSet {
            allele_string: "A/G",
            pos: 100,
            strand: 1,
        };
        assert_eq!(get_matched_variant_alleles(&a, &g).len(), 1);
    }

    /// A region-layout variation source for chromosome 21 over `lines`, each a
    /// record whose tab-separated cells follow the `variation_cols` order below.
    fn region_source(dir: &std::path::Path, lines: &[&str]) -> VariationSource {
        let region_dir = dir.join("variations").join("21");
        std::fs::create_dir_all(&region_dir).unwrap();
        let mut body = Vec::new();
        {
            let mut encoder =
                flate2::write::GzEncoder::new(&mut body, flate2::Compression::default());
            for line in lines {
                std::io::Write::write_all(&mut encoder, line.replace('\t', " ").as_bytes())
                    .unwrap();
                std::io::Write::write_all(&mut encoder, b"\n").unwrap();
            }
            encoder.finish().unwrap();
        }
        std::fs::write(region_dir.join("1-50000000_var.gz"), body).unwrap();
        let info = vep_core::cache_info::CacheInfo {
            variation_cols: [
                "variation_name",
                "failed",
                "somatic",
                "start",
                "end",
                "allele_string",
                "strand",
                "clin_sig",
                "phenotype_or_disease",
                "clin_sig_allele",
                "clin_sig_ref_allele",
                "clinical_impact",
                "pubmed",
                "var_synonyms",
            ]
            .iter()
            .map(|name| name.to_string())
            .collect(),
            ..Default::default()
        };
        VariationSource::open(
            dir,
            &info,
            false,
            std::sync::Arc::new(vep_core::cache_manifest::CacheManifest::absent(dir)),
        )
        .unwrap()
        .expect("a region file is a variation source")
    }

    // The `-o vep` path streams the Extra column straight into the writer, and
    // the end-to-end concordance comparison exercises only the no-flag workload,
    // so these tests pin the bytes of every branch: VEP's key order, the `;`
    // separator, the `-` placeholder, and the `;`-only escaping.

    use vep_core::consequence::{Impact, TranscriptConsequence};
    use vep_io::output::fields::{ExtraFieldsPlan, FieldOptions};

    /// Render `write_extra_fields` into a `String` for byte-level assertions.
    fn extra(
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        options: &FieldOptions,
    ) -> String {
        let mut buf: Vec<u8> = Vec::new();
        let plan = ExtraFieldsPlan::new(options.clone());
        write_extra_fields(variant, tc, &plan, &mut buf);
        String::from_utf8(buf).expect("Extra field is valid UTF-8")
    }

    #[test]
    fn default_row_carries_impact_and_strand() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            impact: Impact::MODERATE,
            strand: -1,
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &FieldOptions::default()),
            "IMPACT=MODERATE;STRAND=-1"
        );
    }

    #[test]
    fn distance_precedes_strand_and_flags_follow_it() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            impact: Impact::MODIFIER,
            strand: 1,
            distance: Some(2500),
            flags: std::sync::Arc::from([
                "cds_start_NF".to_string(),
                "cds_end_NF".to_string(),
                "gencode_basic".to_string(),
            ]),
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &FieldOptions::default()),
            "IMPACT=MODIFIER;DISTANCE=2500;STRAND=1;FLAGS=cds_start_NF,cds_end_NF"
        );
    }

    /// Under `--regulatory` alone a transcript row and the intergenic companion
    /// row print what they print without the flag, through the fast path, while
    /// a motif row prints its own fields and a regulatory row its BIOTYPE, and
    /// `--biotype` beside the flag sends the transcript row down the general
    /// path, where its BIOTYPE prints.
    #[test]
    fn transcript_rows_keep_the_fast_path_under_the_regulatory_group_alone() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let transcript = TranscriptConsequence {
            impact: Impact::MODIFIER,
            strand: 1,
            consequences: std::iter::once(vep_core::consequence::Consequence::UpstreamGeneVariant)
                .collect(),
            distance: Some(120),
            biotype: Some("protein_coding".into()),
            ..Default::default()
        };
        let regulatory = FieldOptions {
            regulatory: true,
            ..Default::default()
        };
        let plan = ExtraFieldsPlan::new(regulatory.clone());
        assert!(plan.any_optional_group && plan.transcript_fast_path);
        assert!(
            !ExtraFieldsPlan::new(FieldOptions {
                regulatory: true,
                biotype: true,
                ..Default::default()
            })
            .transcript_fast_path
        );
        assert_eq!(
            extra(&variant, Some(&transcript), &regulatory),
            extra(&variant, Some(&transcript), &FieldOptions::default())
        );
        assert_eq!(extra(&variant, None, &regulatory), "IMPACT=MODIFIER");
        assert_eq!(
            extra(
                &variant,
                Some(&transcript),
                &FieldOptions {
                    regulatory: true,
                    biotype: true,
                    ..Default::default()
                }
            ),
            "IMPACT=MODIFIER;DISTANCE=120;STRAND=1;BIOTYPE=protein_coding"
        );
        let motif = TranscriptConsequence {
            transcript_id: "ENSM1".into(),
            feature_type: FeatureType::MotifFeature,
            impact: Impact::MODIFIER,
            strand: -1,
            motif: Some(Box::new(vep_core::consequence::MotifFields {
                name: "ENSPFM1".into(),
                pos: Some(3),
                factors: ["A".into(), "B".into()].into(),
            })),
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&motif), &regulatory),
            "IMPACT=MODIFIER;STRAND=-1;MOTIF_NAME=ENSPFM1;MOTIF_POS=3;TRANSCRIPTION_FACTORS=A,B"
        );
        let enhancer = TranscriptConsequence {
            transcript_id: "ENSR1".into(),
            feature_type: FeatureType::RegulatoryFeature,
            impact: Impact::MODIFIER,
            biotype: Some("enhancer".into()),
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&enhancer), &regulatory),
            "IMPACT=MODIFIER;BIOTYPE=enhancer"
        );
    }

    /// A `cds_` flag carrying a tab or a line end prints it percent-encoded on
    /// the fast path of the unconditional group as on the general path, so the
    /// two paths write the same Extra column and the row keeps its columns.
    #[test]
    fn a_flag_with_a_control_byte_is_encoded_on_both_paths() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            impact: Impact::MODIFIER,
            strand: 1,
            flags: std::sync::Arc::from(["cds_\tx".to_string(), "cds_end_NF\n".to_string()]),
            ..Default::default()
        };
        let fast = extra(&variant, Some(&tc), &FieldOptions::default());
        assert_eq!(
            fast,
            "IMPACT=MODIFIER;STRAND=1;FLAGS=cds_%09x,cds_end_NF%0A"
        );
        let general = extra(
            &variant,
            Some(&tc),
            &FieldOptions {
                symbol: true,
                ..Default::default()
            },
        );
        assert_eq!(general, fast);
    }

    #[test]
    fn intergenic_row_carries_impact_alone() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            extra(&variant, None, &FieldOptions::default()),
            "IMPACT=MODIFIER"
        );
    }

    /// The flag groups come out in `@FLAG_FIELDS` order (Constants.pm:60-138);
    /// an all-zero `PHENO` list is dropped (OutputFactory.pm:1136-1138).
    #[test]
    fn optional_groups_follow_vep_order() {
        let mut variant = InputVariant::new("21".into(), 100, 101, b"AC".to_vec(), b"-".to_vec());
        variant.colocated_variants.push(ColocatedVariant {
            id: "rs123".to_string(),
            allele_string: Some("AC/-".to_string()),
            clin_sig: vec!["pathogenic".to_string(), "likely_pathogenic".to_string()],
            somatic: true,
            phenotype: false,
            matched_alleles: vec![MatchedAllele {
                input_allele: "-".to_string(),
                input_index: 0,
                record_allele: "-".to_string(),
                record_index: 0,
            }],
            frequency_strings: [("AF".to_string(), "-:0.003".to_string())].into(),
            ..Default::default()
        });
        let tc = TranscriptConsequence {
            impact: Impact::HIGH,
            strand: 1,
            hgvsc: Some("ENST00000366667.4:c.428C>T".to_string()),
            hgvsp: Some("ENSP00000355627.3:p.Ala143=".to_string()),
            sift: Some("deleterious(0.01)".to_string()),
            ..Default::default()
        };
        let options = FieldOptions {
            variant_class: true,
            hgvs: true,
            sift: true,
            af: true,
            check_existing: true,
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &options),
            "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=deletion;SIFT=deleterious(0.01);HGVSc=ENST00000366667.4:c.428C>T;HGVSp=ENSP00000355627.3:p.Ala143%3D;AF=0.0030;CLIN_SIG=pathogenic,likely_pathogenic;SOMATIC=1"
        );
    }

    #[test]
    fn plugin_keys_follow_flag_fields_alphabetically_and_escape_semicolons_only() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let mut tc = TranscriptConsequence {
            impact: Impact::LOW,
            strand: 1,
            ..Default::default()
        };
        tc.plugin_data
            .insert("ZPlugin".to_string(), "a=b;c".to_string());
        tc.plugin_data
            .insert("APlugin".to_string(), "x".to_string());
        assert_eq!(
            extra(&variant, Some(&tc), &FieldOptions::default()),
            "IMPACT=LOW;STRAND=1;APlugin=x;ZPlugin=a=b%3Bc"
        );
    }

    /// Ensembl VEP sets the identifier keys only for a cache value other than
    /// `-` (OutputFactory.pm:1534-1574 at release/116.2), so the Extra column of
    /// a transcript with no CCDS, RefSeq or UniProt cross-reference carries no
    /// `CCDS=-`, `RefSeq=-`, `SWISSPROT=-`, `TREMBL=-`, `UNIPARC=-` or
    /// `UNIPROT_ISOFORM=-`; the values present print in VEP's key order, lists
    /// joined by `,`.
    #[test]
    fn extra_column_skips_the_dash_sentinel() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            impact: Impact::LOW,
            strand: 1,
            ccds: Some("-".to_string()),
            refseq: Some("-".to_string()),
            swissprot: Some("-".to_string()),
            trembl: Some("L7XE61_HUMAN,L7XCZ9_HUMAN".to_string()),
            uniparc: Some("UPI000002DB1C".to_string()),
            uniprot_isoform: Some("-".to_string()),
            ..Default::default()
        };
        let options = FieldOptions {
            ccds: true,
            uniprot: true,
            xref_refseq: true,
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &options),
            "IMPACT=LOW;STRAND=1;TREMBL=L7XE61_HUMAN,L7XCZ9_HUMAN;UNIPARC=UPI000002DB1C"
        );
    }

    /// A co-located variant whose `clinical_impact` is the defined empty string
    /// (a region file's empty cell) writes `CLINVAR_SOMATIC_CLASSIFICATION=` under
    /// `--clinvar_somatic_classification`, the key Ensembl VEP 116.2 sets to `''`
    /// (OutputFactory.pm:1093) and its default writer prints (VEP_output.pm:117),
    /// and nothing without the flag or with the cell unset.
    #[test]
    fn extra_column_writes_the_empty_clinvar_somatic_classification() {
        use vep_core::variation::ClinVarDetail;

        let mut variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let mut known = ColocatedVariant {
            id: "rs1".into(),
            start: 100,
            end: 100,
            allele_string: Some("A/G".into()),
            strand: 1,
            ..ColocatedVariant::default()
        };
        variant.colocated_variants = vec![known.clone()];
        let options = FieldOptions {
            check_existing: true,
            clinvar_somatic_classification: true,
            ..Default::default()
        };
        assert_eq!(extra(&variant, None, &options), "IMPACT=MODIFIER");
        known.clinvar = Some(Box::new(ClinVarDetail {
            clin_sig_allele: Some(String::new()),
            clin_sig_ref_allele: None,
            clinical_impact: Some(String::new()),
        }));
        variant.colocated_variants = vec![known];
        assert_eq!(
            extra(&variant, None, &options),
            "IMPACT=MODIFIER;CLINVAR_SOMATIC_CLASSIFICATION="
        );
        let without = FieldOptions {
            check_existing: true,
            ..Default::default()
        };
        assert_eq!(extra(&variant, None, &without), "IMPACT=MODIFIER");
    }
}
