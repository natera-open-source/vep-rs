// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's `--minimal` over the alleles vep-rs parsed from one record.
//!
//! Ensembl VEP minimises in two places (ensembl-vep release/116.2). A record with
//! one alternate allele goes through `minimise_alleles` (Parser.pm:964-1019) as
//! it is parsed: `trim_sequences` (ensembl-variation Utils/Sequence.pm:965-1038)
//! trims the bases the two alleles share from the left while both are non-empty,
//! then from the right, the start and end moving in by the bases trimmed on their
//! side and an emptied allele written `-`; the record is marked `minimised`
//! whether or not a base was trimmed (:1011). A record with several alternate
//! alleles is split in the input buffer after it was validated
//! (InputBuffer.pm:638-718): every allele is trimmed against the reference the
//! same way, and when any allele changed (:665-666, :700-705) the record becomes
//! one feature per allele, each annotated on its own trimmed span, every allele
//! after the first numbered by its written position (`_base_allele_number`,
//! :684, read by VariationFeatureOverlap.pm:477), the first carrying the record's
//! chopped allele string and span as its originals (:694-697) and the `minimised`
//! mark. At output the features are merged back (OutputFactory.pm:2272-2354): the
//! record prints its own span and allele string, the alleles print trimmed, and
//! an allele that once trimmed equals the reference, or is a `*` or `<DEL:*>`,
//! writes no row (VariationFeatureOverlap.pm:481, :484). A record none of whose
//! alleles changed is left as it stands, with no mark.

use std::sync::Arc;

use crate::variant::{classify_variant, so_variation_class, InputVariant};

/// What the `--minimal` split made of one allele of a multi-allelic record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MinimalSplit {
    /// No allele of the record changed: the record is annotated as it stands.
    Unchanged,
    /// The record split and this allele is annotated on its trimmed span.
    Trimmed,
    /// The record split and this allele, once trimmed, equals the reference or
    /// is a `*` or `<DEL:*>`, so Ensembl VEP writes no row for it.
    NoAllele,
}

/// `trim_sequences` (Utils/Sequence.pm:965-1038) on one pair with its span:
/// the trimmed reference and allele (`-` when emptied), the span they leave
/// on the forward strand (the start moves by the bases trimmed on the left,
/// the end by those trimmed on the right), and whether a base was trimmed.
pub fn trim_sequences(
    reference: &[u8],
    alt: &[u8],
    start: u64,
    end: u64,
) -> (Vec<u8>, Vec<u8>, u64, u64, bool) {
    // The left loop runs while both alleles have a base (`$ref && $alt`), so the
    // shared prefix stops at the shorter allele; the right loop then runs on what
    // is left of each. A `-` is a byte like any other here, as it is for `eq`.
    let prefix = reference
        .iter()
        .zip(alt)
        .take_while(|(r, a)| r == a)
        .count();
    let (kept_ref, kept_alt) = (&reference[prefix..], &alt[prefix..]);
    let suffix = kept_ref
        .iter()
        .rev()
        .zip(kept_alt.iter().rev())
        .take_while(|(r, a)| r == a)
        .count();
    let dash_if_empty = |allele: &[u8]| {
        if allele.is_empty() {
            b"-".to_vec()
        } else {
            allele.to_vec()
        }
    };
    (
        dash_if_empty(&kept_ref[..kept_ref.len() - suffix]),
        dash_if_empty(&kept_alt[..kept_alt.len() - suffix]),
        start + prefix as u64,
        end.saturating_sub(suffix as u64),
        prefix > 0 || suffix > 0,
    )
}

/// Whether `split_variants` changes a record whose chopped working string is
/// `reference` and `alts` (InputBuffer.pm:665-666, :700): any allele that trims,
/// an allele equal to the reference included.
pub fn split_changes_record<'a>(reference: &str, alts: impl IntoIterator<Item = &'a str>) -> bool {
    alts.into_iter().any(|alt| {
        let (_, _, _, _, changed) = trim_sequences(reference.as_bytes(), alt.as_bytes(), 0, 0);
        changed
    })
}

/// `minimise_alleles` under `--minimal` on a bi-allelic sequence record
/// (Parser.pm:987-1014): the alleles trimmed by [`trim_sequences`], the span
/// moved with them, and the record marked `minimised` whether or not a base was
/// trimmed. The written string stays in `uploaded_allele_string` for
/// `UPLOADED_ALLELE`. A structural record, a record with no allele to annotate
/// and a multi-allelic record (one carrying `record_allele_string_multi`) are
/// left alone: the split handles the last ([`apply_minimal_split`]).
pub fn minimise_bi_allelic(variant: &mut InputVariant) {
    if variant.is_structural
        || variant.no_allele
        || variant.record_allele_string_multi.is_some()
        || variant.alt_alleles.len() != 1
    {
        return;
    }
    let (reference, alt, start, end, changed) = trim_sequences(
        &variant.ref_allele,
        &variant.alt_alleles[0],
        variant.start,
        variant.end,
    );
    if changed {
        set_alleles(variant, reference, alt, start, end);
    }
    variant.minimised = true;
}

/// `split_variants` under `--minimal` on one allele of a multi-allelic record
/// (InputBuffer.pm:638-718), the record's chopped working string in
/// `record_allele_string_multi` and its span in `start` and `end`, `own` the
/// allele's position among that string's alternate alleles (its written position
/// in a VCF record unless a gVCF block before it was left out). When no allele of
/// the record trims, the record is annotated as it stands and carries no
/// `minimised` mark ([`MinimalSplit::Unchanged`]). Otherwise this allele takes
/// its trimmed alleles and span, the record's span as `original_start` and
/// `original_end` for the columns the rejoin restores, the chopped string as its
/// `UPLOADED_ALLELE` (:694), the class `VARIANT_CLASS` prints under the split
/// (the first allele's trimmed pair: OutputFactory.pm:2289 builds the first
/// feature's hash, which caches `class_SO_term`, before the rejoin; a record
/// carrying a class already, a resolved identifier's stored one, keeps it, as
/// `class_SO_term` keeps a set value, VariationFeature.pm:1117), its written
/// position as its ALLELE_NUM (:684) and, after the first allele, the key prefix
/// the rejoin sorts its rows by; an allele that once trimmed equals the reference
/// or is a `*` or `<DEL:*>` (`*C` against `C` leaves `-/*`) is
/// [`MinimalSplit::NoAllele`]: Ensembl VEP builds the overlap alleles from the
/// trimmed string and passes both over (VariationFeatureOverlap.pm:481, :484), the
/// `*` numbered first, so the alleles after it keep their written positions. A
/// structural record, a record with no allele to annotate and a bi-allelic record
/// are left alone.
pub fn apply_minimal_split(variant: &mut InputVariant, own: usize) -> MinimalSplit {
    if variant.is_structural || variant.no_allele {
        return MinimalSplit::Unchanged;
    }
    let Some(working) = variant.record_allele_string_multi.clone() else {
        return MinimalSplit::Unchanged;
    };
    let mut parts = working.split('/');
    let reference = parts.next().unwrap_or("");
    let alts: Vec<&str> = parts.collect();
    if !split_changes_record(reference, alts.iter().copied()) {
        variant.minimised = false;
        return MinimalSplit::Unchanged;
    }
    let Some(own_alt) = alts.get(own) else {
        return MinimalSplit::Unchanged;
    };
    let (start, end) = (variant.start, variant.end);
    let (trimmed_ref, trimmed_alt, new_start, new_end, _) =
        trim_sequences(reference.as_bytes(), own_alt.as_bytes(), start, end);
    if variant.record_so_term.is_none() {
        variant.record_so_term = Some(first_feature_class(&working));
    }
    variant.original_start = Some(start);
    variant.original_end = Some(end);
    variant.minimised = true;
    variant.allele_num = None;
    // The record's string as uploaded is every feature's `original_allele_string`
    // (Parser.pm:982, set before the split and inherited by each copy); the first
    // feature's is then overwritten with the working string (:694), which is what
    // every allele prints after the rejoin. A feature written on its own after a
    // `--check_frequency` drop prints the string it kept.
    if variant.original_allele_string.is_none() {
        variant.original_allele_string = variant
            .uploaded_allele_string
            .as_deref()
            .map(str::to_string);
    }
    variant.uploaded_allele_string = Some(Arc::clone(&working));
    if rowless_trimmed_pair(&trimmed_ref, &trimmed_alt) {
        return MinimalSplit::NoAllele;
    }
    set_split_feature(variant, own, trimmed_ref, trimmed_alt, new_start, new_end);
    MinimalSplit::Trimmed
}

/// Whether the split feature of a pair [`trim_sequences`] left writes no row:
/// the alleles equal (VariationFeatureOverlap.pm:481), or the alternate a `*`
/// or `<DEL:*>` (:484).
fn rowless_trimmed_pair(trimmed_ref: &[u8], trimmed_alt: &[u8]) -> bool {
    trimmed_ref == trimmed_alt || matches!(trimmed_alt, b"*" | b"<DEL:*>")
}

/// The reference allele's row (`--process_ref_homs`) of one allele's split
/// feature, `allele` as [`apply_minimal_split`] left it with `own` the position
/// it was split at. Ensembl VEP builds the reference allele of every feature
/// before it tests the feature's alternate (VariationFeatureOverlap.pm:469-481)
/// and the rejoin carries every split feature's overlaps into the record
/// (OutputFactory.pm:2299-2340), so the record writes one reference row per
/// allele, each the allele's own trimmed reference over its trimmed span, an
/// allele equal to the reference once trimmed included: its feature `-/-` writes
/// the row and nothing else.
pub fn split_reference_row(allele: &InputVariant, own: usize) -> InputVariant {
    split_feature(allele, own)
        .unwrap_or_else(|| allele.clone())
        .reference_allele_row()
}

/// The feature of allele `own` of a `--minimal` split record, built from any
/// allele of the record: the working string's alternate allele `own` trimmed
/// against the reference over the record's span (InputBuffer.pm:658-676), with
/// the key prefix of its place. `None` for an allele the split left alone.
fn split_feature(allele: &InputVariant, own: usize) -> Option<InputVariant> {
    let working = allele.record_allele_string_multi.as_deref()?;
    let (start, end) = (allele.original_start?, allele.original_end?);
    let mut parts = working.split('/');
    let reference = parts.next().unwrap_or("");
    let own_alt = parts.nth(own)?;
    let (trimmed_ref, trimmed_alt, new_start, new_end, _) =
        trim_sequences(reference.as_bytes(), own_alt.as_bytes(), start, end);
    let mut feature = allele.clone();
    set_split_feature(
        &mut feature,
        own,
        trimmed_ref,
        trimmed_alt,
        new_start,
        new_end,
    );
    Some(feature)
}

/// The first feature of a `--minimal` split record whose first alternate
/// allele writes no row, for the lookups the record's rows print. Ensembl VEP
/// annotates every split feature on its own span and its rows print the first
/// feature's co-located variants and overlapping structural variants: the
/// rejoin (OutputFactory.pm:2272-2353) moves the later features' consequences
/// into the first and leaves `existing` and `overlapping_svs` where they are,
/// and the rows are built from the first feature (:342-349). An allele equal
/// to the reference once trimmed, a `*` and a `<DEL:*>` write no row
/// (VariationFeatureOverlap.pm:469-474, :481, :484), so when the first
/// alternate allele is one of those no allele of the batch carries the
/// feature, and each allele looks it up here. `None` for every other allele.
pub fn rowless_first_feature(allele: &InputVariant) -> Option<InputVariant> {
    let working = allele.record_allele_string_multi.as_deref()?;
    allele.original_start?;
    let mut parts = working.split('/');
    let reference = parts.next()?;
    let first = parts.next()?;
    if !rowless_alternate_allele(reference, first) {
        return None;
    }
    split_feature(allele, 0)
}

/// Whether the split feature of `alt`, an alternate allele of a record's working
/// string, writes no row: once the pair is trimmed (InputBuffer.pm:665), an allele
/// equal to `reference` (VariationFeatureOverlap.pm:481) or a `*` or `<DEL:*>`,
/// written so or left by the trim (`*C` against `C`) (:469-474, :484).
pub fn rowless_alternate_allele(reference: &str, alt: &str) -> bool {
    let (trimmed_ref, trimmed_alt, _, _, _) =
        trim_sequences(reference.as_bytes(), alt.as_bytes(), 0, 0);
    rowless_trimmed_pair(&trimmed_ref, &trimmed_alt)
}

/// `variant`, an allele of a `--minimal` split record or the reference row leading
/// it, as the feature of its own that Ensembl VEP writes once a `--check_frequency`
/// drop has reset the buffer and its rejoin mark (BaseCacheVariation.pm:365-369,
/// InputBuffer.pm:535-537 and :610), so that OutputFactory.pm:294 skips the
/// rejoin: the feature prints its own span, string and class, is named by the
/// record's `variation_name` (Parser.pm:462, given before the split and inherited
/// by every copy), and as the first feature keeps the `minimised` mark and the
/// working string as its `UPLOADED_ALLELE` (InputBuffer.pm:694-696), as a later
/// one prints the record's uploaded string and no mark. `own` tells the
/// features of one record apart ([`InputVariant::standalone_feature`]); the
/// first feature is the one without a rejoin key prefix.
pub fn feature_on_its_own(variant: &mut InputVariant, own: u32) {
    variant.id = Some(variant.variation_name());
    let is_first_feature = variant.rejoin_key_prefix.is_none();
    // The class the split gave every allele is the first feature's, so the
    // feature's own class takes its place; a class the record carried before the
    // split, a resolved identifier's, stays on every copy (VariationFeature.pm:1117).
    let split_class = variant
        .record_allele_string_multi
        .as_deref()
        .map(first_feature_class);
    if variant.record_so_term == split_class {
        variant.record_so_term = None;
    }
    variant.standalone_feature = Some(own);
    variant.original_start = None;
    variant.original_end = None;
    variant.record_allele_string_multi = None;
    variant.rejoin_key_prefix = None;
    if !is_first_feature {
        variant.minimised = false;
        if let Some(uploaded) = &variant.original_allele_string {
            variant.uploaded_allele_string = Some(Arc::from(uploaded.as_str()));
        }
    }
}

/// The class `VARIANT_CLASS` prints for a `--minimal` split record: that of the
/// first alternate allele of `working`, the record's chopped allele string, trimmed
/// against the reference (OutputFactory.pm:2289 builds the first feature's hash,
/// which caches `class_SO_term`, before the rejoin).
fn first_feature_class(working: &str) -> String {
    let mut parts = working.split('/');
    let reference = parts.next().unwrap_or("");
    let first = parts.next().unwrap_or("");
    let (trimmed_ref, trimmed_alt, _, _, _) =
        trim_sequences(reference.as_bytes(), first.as_bytes(), 0, 0);
    so_variation_class(&allele_pair(&trimmed_ref, &trimmed_alt)).to_string()
}

/// The `REF_ALLELE` of every row of a split record, `working` the record's
/// chopped allele string: the rows print the record's feature, whose
/// `ref_allele_string` was cached while the first allele's rows were built
/// before the rejoin restored the record's string (OutputFactory.pm:1357 at
/// :2289, VariationFeature.pm:1533-1539), so it is the first allele's trimmed
/// reference, or the record's own when the first allele writes no row
/// ([`rowless_alternate_allele`]: equal to the reference once trimmed, or a `*`
/// or `<DEL:*>`) and so built none to cache it; with `reference_rows`
/// (`--process_ref_homs`, and the row not filtered away) that feature builds the
/// reference allele's row, so its trimmed reference is the one cached.
pub fn rejoined_reference(working: &str, reference_rows: bool) -> String {
    let mut parts = working.split('/');
    let reference = parts.next().unwrap_or("");
    let first = parts.next().unwrap_or("");
    let (trimmed_ref, trimmed_alt, _, _, _) =
        trim_sequences(reference.as_bytes(), first.as_bytes(), 0, 0);
    if rowless_trimmed_pair(&trimmed_ref, &trimmed_alt) && !reference_rows {
        return reference.to_string();
    }
    String::from_utf8_lossy(&trimmed_ref).into_owned()
}

/// `variant` as the split feature of one allele: its trimmed pair over the
/// trimmed span and, after the first allele, the key prefix the rejoin sorts its
/// rows by (OutputFactory.pm:2316, the feature's allele string before the key).
fn set_split_feature(
    variant: &mut InputVariant,
    own: usize,
    trimmed_ref: Vec<u8>,
    trimmed_alt: Vec<u8>,
    start: u64,
    end: u64,
) {
    variant.rejoin_key_prefix =
        (own > 0).then(|| format!("{}_", allele_pair(&trimmed_ref, &trimmed_alt)).into());
    set_alleles(variant, trimmed_ref, trimmed_alt, start, end);
}

/// The allele string of one pair, `reference/alt`.
fn allele_pair(reference: &[u8], alt: &[u8]) -> String {
    format!(
        "{}/{}",
        String::from_utf8_lossy(reference),
        String::from_utf8_lossy(alt)
    )
}

/// Gives `variant` the alleles and span a trim left.
fn set_alleles(variant: &mut InputVariant, reference: Vec<u8>, alt: Vec<u8>, start: u64, end: u64) {
    variant.allele_string = allele_pair(&reference, &alt);
    variant.variant_class = classify_variant(&reference, &alt);
    variant.ref_allele = reference;
    variant.alt_alleles = vec![alt];
    variant.start = start;
    variant.end = end;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `trim_sequences` trims from the left while both alleles have a base, then
    /// from the right, moves the span with the trims and writes an emptied allele
    /// `-`: `GAAG/GGGG` at 14-17 is `AA/GG` at 15-16, `G/GT` at 14 is `-/T` at
    /// 15-14, `GAA/G` at 14-16 is `AA/-` at 15-16 and `A/A` at 15 is `-/-` at 16-15
    /// (every base trimmed, so `changed`); `A/G` and `-/T` are left as they are.
    #[test]
    fn trim_sequences_moves_the_span_with_the_trims() {
        let trim = |r: &str, a: &str, s, e| {
            let (tr, ta, ts, te, changed) = trim_sequences(r.as_bytes(), a.as_bytes(), s, e);
            (
                String::from_utf8(tr).unwrap(),
                String::from_utf8(ta).unwrap(),
                ts,
                te,
                changed,
            )
        };
        assert_eq!(
            trim("GAAG", "GGGG", 14, 17),
            ("AA".into(), "GG".into(), 15, 16, true)
        );
        assert_eq!(
            trim("G", "GT", 14, 14),
            ("-".into(), "T".into(), 15, 14, true)
        );
        assert_eq!(
            trim("GAA", "G", 14, 16),
            ("AA".into(), "-".into(), 15, 16, true)
        );
        assert_eq!(
            trim("A", "A", 15, 15),
            ("-".into(), "-".into(), 16, 15, true)
        );
        assert_eq!(
            trim("A", "G", 15, 15),
            ("A".into(), "G".into(), 15, 15, false)
        );
        assert_eq!(
            trim("-", "T", 15, 14),
            ("-".into(), "T".into(), 15, 14, false)
        );
        assert_eq!(
            trim("A", "*", 15, 15),
            ("A".into(), "*".into(), 15, 15, false)
        );
    }

    /// A record splits when any allele trims, an allele equal to the reference
    /// included, and not when none does: `G/T/C` stays, `G/T/GA`, `A/A/G` and
    /// `GAA/GGA/GGA` split.
    #[test]
    fn a_record_splits_when_any_allele_trims() {
        assert!(!split_changes_record("G", ["T", "C"]));
        assert!(!split_changes_record("G", ["*", "<DEL:*>"]));
        assert!(split_changes_record("G", ["T", "GA"]));
        assert!(split_changes_record("A", ["A", "G"]));
        assert!(split_changes_record("GAA", ["GGA", "GGA"]));
        assert!(!split_changes_record("-", ["A", "*"]));
    }

    /// A bi-allelic record is marked minimised whether or not a base was trimmed:
    /// `GAAG/GGGG` at 14-17 becomes `AA/GG` at 15-16 and `A/G` keeps its alleles;
    /// a structural record and a multi-allelic one are left alone.
    #[test]
    fn a_bi_allelic_record_is_trimmed_and_marked() {
        let mut mnp = InputVariant::new("21".into(), 14, 17, b"GAAG".to_vec(), b"GGGG".to_vec());
        minimise_bi_allelic(&mut mnp);
        assert_eq!(
            (mnp.start, mnp.end, mnp.allele_string.as_str()),
            (15, 16, "AA/GG")
        );
        assert!(mnp.minimised);
        let mut snv = InputVariant::new("21".into(), 15, 15, b"A".to_vec(), b"G".to_vec());
        minimise_bi_allelic(&mut snv);
        assert_eq!(
            (snv.start, snv.end, snv.allele_string.as_str()),
            (15, 15, "A/G")
        );
        assert!(snv.minimised);
        let mut structural =
            InputVariant::new("21".into(), 15, 20, b"N".to_vec(), b"<DEL>".to_vec());
        structural.is_structural = true;
        minimise_bi_allelic(&mut structural);
        assert!(!structural.minimised);
        let mut multi = InputVariant::new("21".into(), 14, 14, b"G".to_vec(), b"GA".to_vec());
        multi.record_allele_string_multi = Some("G/T/GA".into());
        minimise_bi_allelic(&mut multi);
        assert!(!multi.minimised);
        assert_eq!(multi.allele_string, "G/GA");
    }

    /// The split of `G/T/GA` at 14 (Ensembl VEP 116.2 on `21 10542414 x G T,GA`
    /// under `--minimal`): the `T` allele keeps its alleles under the record's
    /// span with no key prefix, the `GA` allele becomes `-/A` at 15-14 with
    /// `original_start` and `original_end` at 14, the prefix `-/A_`, the chopped
    /// string as its `UPLOADED_ALLELE` and the class of the first trimmed pair,
    /// `SNV`; both are minimised and numbered by their written position.
    #[test]
    fn a_split_allele_takes_its_trimmed_span_and_the_record_s_span_as_its_original() {
        let working: Arc<str> = "G/T/GA".into();
        let mut first = InputVariant::new("21".into(), 14, 14, b"G".to_vec(), b"T".to_vec());
        first.record_allele_string_multi = Some(Arc::clone(&working));
        first.uploaded_allele_string = Some("G/T/GA".into());
        assert_eq!(apply_minimal_split(&mut first, 0), MinimalSplit::Trimmed);
        assert_eq!(
            (first.start, first.end, first.allele_string.as_str()),
            (14, 14, "G/T")
        );
        assert_eq!(
            (first.original_start, first.original_end),
            (Some(14), Some(14))
        );
        assert!(first.minimised);
        assert_eq!(first.rejoin_key_prefix, None);
        assert_eq!(first.record_so_term.as_deref(), Some("SNV"));

        let mut second = InputVariant::new("21".into(), 14, 14, b"G".to_vec(), b"GA".to_vec());
        second.allele_index = 1;
        second.allele_num = Some(1);
        second.record_allele_string_multi = Some(Arc::clone(&working));
        assert_eq!(apply_minimal_split(&mut second, 1), MinimalSplit::Trimmed);
        assert_eq!(
            (second.start, second.end, second.allele_string.as_str()),
            (15, 14, "-/A")
        );
        assert_eq!(
            (second.original_start, second.original_end),
            (Some(14), Some(14))
        );
        assert_eq!(second.rejoin_key_prefix.as_deref(), Some("-/A_"));
        assert_eq!(second.uploaded_allele_string.as_deref(), Some("G/T/GA"));
        assert_eq!(second.allele_num, None);
        assert!(second.minimised);
        assert_eq!(second.location(), "21:14");
        assert_eq!(second.display_allele(), "A");
        assert_eq!(second.class_so_term(), "SNV");
    }

    /// `REF_ALLELE` on a split record is the first allele's trimmed reference (`A`
    /// for `AAG/GAG/AGG/-` and for `AG/GG/-`, `-` for `G/GA/T`), or the record's
    /// own reference when the first allele equals it once trimmed (`A` for
    /// `A/A/G`), as Ensembl VEP 116.2 prints under `--minimal --show_ref_allele`
    /// for `21 10542414 AAAG AGAG,AAGG,A`, `21 10542414 AAG AGG,A`,
    /// `21 10542414 G GA,T` and `21 10542415 A A,G`; under `--process_ref_homs`
    /// the same run prints `-` for `A/A/G`, the trimmed reference of the first
    /// allele's feature, whose reference row cached it.
    #[test]
    fn the_rejoined_reference_is_the_first_alleles_unless_it_equals_the_reference() {
        assert_eq!(rejoined_reference("AAG/GAG/AGG/-", false), "A");
        assert_eq!(rejoined_reference("AG/GG/-", false), "A");
        assert_eq!(rejoined_reference("A/A/G", false), "A");
        assert_eq!(rejoined_reference("A/T/AT", false), "A");
        assert_eq!(rejoined_reference("G/GA/T", false), "-");
        assert_eq!(rejoined_reference("A/A/G", true), "-");
        assert_eq!(rejoined_reference("AAG/GAG/AGG/-", true), "A");
    }

    /// Under `--process_ref_homs` every allele of a split record contributes the
    /// reference row of its own feature, as Ensembl VEP 116.2 prints for
    /// `21 10542414 A T,AG`: `A` numbered 0 over 10542414 for the SNV's feature
    /// `A/T`, and `-` numbered 0 over the insertion span 10542415-10542414 for the
    /// feature `-/G`, sorted under that feature's key `-/G_`; and for
    /// `21 10542415 A A,G`, whose first allele equals the reference once trimmed,
    /// `-` over 10542416-10542415 for the feature `-/-`, the row that allele's
    /// feature writes and nothing else.
    #[test]
    fn every_split_allele_contributes_its_own_reference_row() {
        let split = |own: usize| {
            let alt: &[u8] = if own == 0 { b"T" } else { b"AG" };
            let mut allele =
                InputVariant::new("21".into(), 10542414, 10542414, b"A".to_vec(), alt.to_vec());
            allele.record_allele_string_multi = Some("A/T/AG".into());
            let outcome = apply_minimal_split(&mut allele, own);
            (allele, outcome)
        };
        let (snv, outcome) = split(0);
        assert_eq!(outcome, MinimalSplit::Trimmed);
        let row = split_reference_row(&snv, 0);
        assert!(row.reference_allele_row);
        assert!(row.minimised);
        assert_eq!(
            (row.start, row.end, row.allele_num),
            (10542414, 10542414, Some(0))
        );
        assert_eq!(row.allele_string, "A/A");
        assert_eq!(row.rejoin_key_prefix, None);
        let (insertion, outcome) = split(1);
        assert_eq!(outcome, MinimalSplit::Trimmed);
        let row = split_reference_row(&insertion, 1);
        assert_eq!(
            (row.start, row.end, row.allele_num),
            (10542415, 10542414, Some(0))
        );
        assert_eq!(row.allele_string, "-/-");
        assert_eq!(row.rejoin_key_prefix.as_deref(), Some("-/G_"));
        assert_eq!(row.uploaded_allele(true).as_deref(), Some("A/T/AG"));
        assert_eq!(row.location(), "21:10542414");

        let mut equal = InputVariant::new(
            "21".into(),
            10542415,
            10542415,
            b"A".to_vec(),
            b"A".to_vec(),
        );
        equal.record_allele_string_multi = Some("A/A/G".into());
        assert_eq!(apply_minimal_split(&mut equal, 0), MinimalSplit::NoAllele);
        let row = split_reference_row(&equal, 0);
        assert!(row.reference_allele_row);
        assert_eq!(
            (row.start, row.end, row.allele_num),
            (10542416, 10542415, Some(0))
        );
        assert_eq!(row.allele_string, "-/-");
        assert_eq!(row.rejoin_key_prefix, None);
        assert_eq!(row.location(), "21:10542415");
        assert_eq!(row.class_so_term(), "sequence_alteration");
    }

    /// The first feature of a split record is built for an allele of the record
    /// when the first alternate allele writes no row: `A/A/G` at 10542415 gives
    /// `-/-` at 10542416-10542415 and `A/*/G/AT` at 100 gives `A/*` at 100, each
    /// without a key prefix, as Ensembl VEP 116.2 splits `21 10542415 A A,G`
    /// (`-/-` over 10542416-10542415, VARIANT_CLASS `sequence_alteration`) and
    /// trims nothing off `*`. `A/AA/G`, whose first allele writes its rows, and a
    /// record the split left alone give none.
    #[test]
    fn the_first_feature_is_built_when_its_allele_writes_no_row() {
        let mut g = InputVariant::new(
            "21".into(),
            10542415,
            10542415,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        g.allele_index = 1;
        g.record_allele_string_multi = Some("A/A/G".into());
        assert_eq!(apply_minimal_split(&mut g, 1), MinimalSplit::Trimmed);
        let first = rowless_first_feature(&g).unwrap();
        assert_eq!(
            (first.start, first.end, first.allele_string.as_str()),
            (10542416, 10542415, "-/-")
        );
        assert_eq!(first.rejoin_key_prefix, None);
        assert_eq!(first.chr, "21");

        let mut insertion = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"AT".to_vec());
        insertion.allele_index = 2;
        insertion.record_allele_string_multi = Some("A/*/G/AT".into());
        assert_eq!(
            apply_minimal_split(&mut insertion, 2),
            MinimalSplit::Trimmed
        );
        assert_eq!((insertion.start, insertion.end), (101, 100));
        let first = rowless_first_feature(&insertion).unwrap();
        assert_eq!(
            (first.start, first.end, first.allele_string.as_str()),
            (100, 100, "A/*")
        );

        let mut a = InputVariant::new("21".into(), 5030602, 5030602, b"A".to_vec(), b"G".to_vec());
        a.allele_index = 1;
        a.record_allele_string_multi = Some("A/AA/G".into());
        assert_eq!(apply_minimal_split(&mut a, 1), MinimalSplit::Trimmed);
        assert!(rowless_first_feature(&a).is_none());

        let mut unsplit = InputVariant::new("21".into(), 14, 14, b"G".to_vec(), b"C".to_vec());
        unsplit.record_allele_string_multi = Some("G/T/C".into());
        assert_eq!(
            apply_minimal_split(&mut unsplit, 1),
            MinimalSplit::Unchanged
        );
        assert!(rowless_first_feature(&unsplit).is_none());
    }

    /// A `*` the trim leaves writes no row: `21 10542415 x A *A,G` under `--minimal`
    /// splits, as `*A` against `A` trims its shared last base to `-/*` at
    /// 10542415-10542414 (a change, InputBuffer.pm:665-666, the end moving in by
    /// the base trimmed), and Ensembl VEP numbers the feature's allele `*`
    /// and passes it over as it does a `*` written so (VariationFeatureOverlap.pm:484).
    /// So the first allele is [`MinimalSplit::NoAllele`], the `G` allele is `A/G`
    /// at 10542415 numbered by its written position, 2 (`_base_allele_number`,
    /// InputBuffer.pm:684), under the key prefix `A/G_`, with the first trimmed
    /// pair's class, `sequence_alteration`, and the first feature no allele carries
    /// is `-/*`. `REF_ALLELE` is the record's own `A`: the `*` feature built no
    /// row to cache its `-` (OutputFactory.pm:1357 runs per row); under the
    /// reference rows its reference row caches `-`. A `<DEL:*>` the trim leaves
    /// is passed over the same way, and a longer allele with a `*` (`*AA`) is not.
    #[test]
    fn a_star_the_trim_leaves_writes_no_row() {
        let working: Arc<str> = "A/*A/G".into();
        let mut star = InputVariant::new(
            "21".into(),
            10542415,
            10542415,
            b"A".to_vec(),
            b"*A".to_vec(),
        );
        star.record_allele_string_multi = Some(Arc::clone(&working));
        assert_eq!(apply_minimal_split(&mut star, 0), MinimalSplit::NoAllele);
        assert!(star.minimised);
        assert_eq!(
            (star.original_start, star.original_end),
            (Some(10542415), Some(10542415))
        );

        let mut g = InputVariant::new(
            "21".into(),
            10542415,
            10542415,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        g.allele_index = 1;
        g.record_allele_string_multi = Some(Arc::clone(&working));
        assert_eq!(apply_minimal_split(&mut g, 1), MinimalSplit::Trimmed);
        assert_eq!(
            (g.start, g.end, g.allele_string.as_str()),
            (10542415, 10542415, "A/G")
        );
        assert_eq!(g.allele_num, None);
        assert_eq!(g.rejoin_key_prefix.as_deref(), Some("A/G_"));
        assert_eq!(g.class_so_term(), "sequence_alteration");
        assert_eq!(g.uploaded_allele_string.as_deref(), Some("A/*A/G"));
        let first = rowless_first_feature(&g).unwrap();
        assert_eq!(
            (first.start, first.end, first.allele_string.as_str()),
            (10542415, 10542414, "-/*")
        );
        assert_eq!(first.rejoin_key_prefix, None);

        assert!(rowless_alternate_allele("A", "*A"));
        assert!(rowless_alternate_allele("A", "*"));
        assert!(rowless_alternate_allele("C", "<DEL:*>C"));
        assert!(!rowless_alternate_allele("A", "*AA"));
        assert!(!rowless_alternate_allele("A", "G"));
        assert_eq!(rejoined_reference("A/*A/G", false), "A");
        assert_eq!(rejoined_reference("A/*A/G", true), "-");
        assert_eq!(rejoined_reference("A/*/AT", false), "A");
    }

    /// A record none of whose alleles trims (`G/T/C`) is left as it stands with no
    /// mark, and an allele equal to the reference once trimmed (`A` of `A/A/G`)
    /// writes no row while its record's class is that of its trimmed pair, `-/-`,
    /// `sequence_alteration`, as Ensembl VEP 116.2 prints for `21 10542415 x A A,G`;
    /// a record that carries a class already (a resolved identifier's stored
    /// `indel`) keeps it through the split.
    #[test]
    fn an_unchanged_record_keeps_its_alleles_and_a_reference_allele_writes_no_row() {
        let mut snv = InputVariant::new("21".into(), 14, 14, b"G".to_vec(), b"C".to_vec());
        snv.allele_index = 1;
        snv.minimised = true;
        snv.record_allele_string_multi = Some("G/T/C".into());
        assert_eq!(apply_minimal_split(&mut snv, 1), MinimalSplit::Unchanged);
        assert!(!snv.minimised);
        assert_eq!(snv.original_start, None);

        let mut same = InputVariant::new("21".into(), 15, 15, b"A".to_vec(), b"A".to_vec());
        same.record_allele_string_multi = Some("A/A/G".into());
        assert_eq!(apply_minimal_split(&mut same, 0), MinimalSplit::NoAllele);
        assert_eq!(same.record_so_term.as_deref(), Some("sequence_alteration"));
        let mut g = InputVariant::new("21".into(), 15, 15, b"A".to_vec(), b"G".to_vec());
        g.allele_index = 1;
        g.allele_num = Some(1);
        g.record_allele_string_multi = Some("A/A/G".into());
        assert_eq!(apply_minimal_split(&mut g, 1), MinimalSplit::Trimmed);
        assert_eq!(g.allele_num, None);
        assert_eq!(g.rejoin_key_prefix.as_deref(), Some("A/G_"));
        assert_eq!(g.class_so_term(), "sequence_alteration");

        let mut stored = InputVariant::new("21".into(), 14, 15, b"TA".to_vec(), b"T".to_vec());
        stored.record_allele_string_multi = Some("TA/T/TAA".into());
        stored.record_so_term = Some("indel".into());
        assert_eq!(apply_minimal_split(&mut stored, 0), MinimalSplit::Trimmed);
        assert_eq!(
            (stored.start, stored.end, stored.allele_string.as_str()),
            (15, 15, "A/-")
        );
        assert_eq!(stored.class_so_term(), "indel");
    }
}
