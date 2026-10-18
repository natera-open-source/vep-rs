// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Consequence filtering: `--pick`, `--most_severe`, `--summary`, `--per_gene`, `--flag_pick`.
//!
//! Implements VEP's pick algorithm that selects the "best" transcript consequence
//! for each variant using a composite scoring system: canonical status, biotype
//! (protein-coding preferred), and consequence severity rank. Supports both
//! removing non-picked consequences and flag-only modes that annotate the pick
//! without discarding other transcripts.
//!
//! Perl citations name modules of ensembl-vep release/116.2 (`Bio/EnsEMBL/VEP/...`).

use std::collections::HashMap;
use std::sync::Arc;

use smallvec::SmallVec;
use vep_core::consequence::{ConsequenceList, TranscriptConsequence};
use vep_core::variant::InputVariant;

/// Configuration for consequence filtering (`--pick`, `--most_severe`, etc.).
///
/// Each boolean corresponds to a VEP CLI flag that controls how transcript
/// consequences are filtered or flagged after annotation.
#[derive(Debug, Clone, Default)]
pub struct FilterConfig {
    pub pick: bool,
    pub pick_allele: bool,
    pub per_gene: bool,
    pub pick_allele_gene: bool,
    pub most_severe: bool,
    pub summary: bool,
    pub flag_pick: bool,
    pub flag_pick_allele: bool,
    pub flag_pick_allele_gene: bool,
}

impl FilterConfig {
    /// Returns true if any filter is active.
    pub fn any_active(&self) -> bool {
        self.pick
            || self.pick_allele
            || self.per_gene
            || self.pick_allele_gene
            || self.most_severe
            || self.summary
            || self.flag_pick
            || self.flag_pick_allele
            || self.flag_pick_allele_gene
    }
}

/// Apply consequence filtering to a variant based on config.
pub fn apply_filters(variant: &mut InputVariant, config: &FilterConfig) {
    if config.most_severe {
        apply_most_severe(variant);
    } else if config.summary {
        apply_summary(variant);
    } else if config.pick {
        apply_pick(variant);
    } else if config.pick_allele {
        apply_pick_allele(variant);
    } else if config.per_gene {
        apply_per_gene(variant);
    } else if config.pick_allele_gene {
        apply_pick_allele_gene(variant);
    } else if config.flag_pick {
        flag_pick(variant);
    } else if config.flag_pick_allele {
        flag_pick_allele(variant);
    } else if config.flag_pick_allele_gene {
        flag_pick_allele_gene(variant);
    }
}

/// The APPRIS pick tier of an attribute value as the cache spells it
/// (`principal1`, `alternative2`, `principalM`): lower is better, 100 when the
/// value yields no grade (`OutputFactory.pm:757-767` at release/116.2).
///
/// Perl matches `m/([A-Za-z]).+(\d+)/`: the first letter of the value, at
/// least one further character, then a digit. With `.+` greedy, the digit
/// captured is the last one in the value, alone, so `principal1` grades 1 and
/// a two-digit tail would grade by its last digit; a value whose only digit
/// sits right after the letter, such as `P1`, does not match. The grade gains
/// 10 when the captured letter is a lowercase `a` (`substr($type, 0, 1) eq
/// 'a'`), then is kept only when non-zero (`$info->{appris} = $grade if
/// $grade`), so `alternative0` is 10 and `principal0` stays 100. `principalM`
/// and `alternativeM` carry no digit and stay 100.
fn appris_rank(appris: Option<&str>) -> u32 {
    let Some(bytes) = appris.map(str::as_bytes) else {
        return 100;
    };
    let Some(letter) = bytes.iter().position(|b| b.is_ascii_alphabetic()) else {
        return 100;
    };
    let Some(digit) = bytes
        .iter()
        .rposition(|b| b.is_ascii_digit())
        .filter(|&digit| digit >= letter + 2)
    else {
        return 100;
    };
    let mut grade = u32::from(bytes[digit] - b'0');
    if bytes[letter] == b'a' {
        grade += 10;
    }
    if grade == 0 {
        return 100;
    }
    grade
}

/// Compute the pick score for a transcript consequence.
/// Lower score = higher priority (gets picked).
///
/// Mirrors Perl VEP's documented `pick_order`
/// (`Bio/EnsEMBL/VEP/Config.pm:311`):
///   `mane_select, mane_plus_clinical, canonical, appris, tsl, biotype, ccds, rank, length`
///
/// Encoding lifted directly from
/// `Bio/EnsEMBL/VEP/OutputFactory.pm::pick_worst_VariationFeatureOverlapAllele`
/// (`:703-812`; the per-transcript values at `:736-767`), reading the attribute
/// values as the cache spells them (slot 3 through [`appris_rank`]):
///
/// | Slot | Field                | Encoding                                          | Default |
/// |------|----------------------|---------------------------------------------------|---------|
/// | 0    | `mane_select`        | `0` if Some(_), else `1`                          | 1       |
/// | 1    | `mane_plus_clinical` | `0` if Some(_), else `1`                          | 1       |
/// | 2    | `canonical`          | `0` if true, else `1`                             | 1       |
/// | 3    | `appris`             | `principal{N}`→N, `alternative{N}`→N+10, else 100 | 100     |
/// | 4    | `tsl`                | the `tsl(\d+)` digits, `100` when none or zero   | 100     |
/// | 5    | `biotype`            | `0` if `protein_coding`, else `1`                 | 1       |
/// | 6    | `ccds`               | `0` if Some(_) and not `"-"`, else `1`            | 1       |
/// | 7    | `rank`               | min consequence rank across `consequences`        | n/a     |
/// | 8    | `length`             | `0 - protein_or_transcript_length` (longer wins)  | 0       |
///
/// Slot 8 (length) is a `0` placeholder: protein and transcript length are not
/// plumbed onto `TranscriptConsequence`, so two transcripts tied on slots 0-7
/// pick the first iterated where Perl's pick is length-deterministic.
///
/// Slot 4 (TSL) note: Perl matches `m/tsl(\d+)/` and keeps the digits only
/// when non-zero (`$info->{tsl} = $1 if $1`, `:752-753`); the loader applies
/// the same rule, so `tsl1 (assigned to previous version 8)` is 1 and `tslNA`
/// is none, here 100.
///
/// Slot 8 (length) note: Perl computes `0 - length(translateable_seq)` for
/// coding transcripts and `0 - length(transcript)` otherwise. The `i64`
/// inversion lets a 2400-aa protein (`-2400`) win over an 800-aa protein
/// (`-800`).
fn pick_score(tc: &TranscriptConsequence) -> (u32, u32, u32, u32, u32, u32, u32, u32, i64) {
    let mane_select_score = if tc.mane_select.is_some() { 0 } else { 1 };
    let mane_plus_clinical_score = if tc.mane_plus_clinical.is_some() {
        0
    } else {
        1
    };
    let canonical_score = if tc.canonical { 0 } else { 1 };
    let appris_score = appris_rank(tc.appris.as_deref());
    let tsl_score: u32 = match tc.tsl {
        Some(digits) if digits > 0 => u32::from(digits),
        _ => 100,
    };
    let biotype_score = if tc.biotype.as_deref() == Some("protein_coding") {
        0
    } else {
        1
    };
    // CCDS: Perl checks `_ccds && _ccds ne '-'`; the placeholder `"-"` is
    // common in cache exports for transcripts with no CCDS annotation.
    let ccds_score = match tc.ccds.as_deref() {
        Some(s) if !s.is_empty() && s != "-" => 0,
        _ => 1,
    };
    let rank_score = tc
        .consequences
        .iter()
        .map(|c| c.rank())
        .min()
        .unwrap_or(u32::MAX);
    // Length comparator placeholder; see the rustdoc above for slot 8 semantics.
    let length_score = 0i64;

    (
        mane_select_score,
        mane_plus_clinical_score,
        canonical_score,
        appris_score,
        tsl_score,
        biotype_score,
        ccds_score,
        rank_score,
        length_score,
    )
}

/// --most_severe: Remove all transcript_consequences, keep only the
/// most_severe_consequence SO term.
fn apply_most_severe(variant: &mut InputVariant) {
    if variant.most_severe_consequence.is_none() && !variant.transcript_consequences.is_empty() {
        variant.most_severe_consequence = variant
            .transcript_consequences
            .iter()
            .flat_map(|tc| tc.consequences.iter())
            .copied()
            .min_by_key(|c| c.rank());
    }
    variant.transcript_consequences.clear();
}

/// --summary: Remove all transcript_consequences. The output line will have
/// a comma-joined list of all unique consequence terms.
fn apply_summary(variant: &mut InputVariant) {
    let mut unique_terms: ConsequenceList = SmallVec::new();
    for tc in &variant.transcript_consequences {
        for c in &tc.consequences {
            if !unique_terms.contains(c) {
                unique_terms.push(*c);
            }
        }
    }
    unique_terms.sort_by_key(|c| c.rank());

    if !unique_terms.is_empty() {
        variant.most_severe_consequence = Some(unique_terms[0]);
        let summary_tc = TranscriptConsequence {
            consequences: unique_terms,
            ..Default::default()
        };
        variant.transcript_consequences = vec![summary_tc];
    } else {
        variant.transcript_consequences.clear();
    }
}

/// --pick: Keep only the single best-scoring consequence across all transcripts.
fn apply_pick(variant: &mut InputVariant) {
    if variant.transcript_consequences.is_empty() {
        return;
    }
    let best_idx = find_best_index(&variant.transcript_consequences);
    let best = variant.transcript_consequences.swap_remove(best_idx);
    variant.transcript_consequences = vec![best];
}

/// --pick_allele: Keep only the best consequence per allele.
/// Multi-allelic variants are split, so each InputVariant has one allele and
/// this behaves the same as --pick.
fn apply_pick_allele(variant: &mut InputVariant) {
    apply_pick(variant);
}

/// --per_gene: Keep only the best-scoring consequence per gene.
fn apply_per_gene(variant: &mut InputVariant) {
    if variant.transcript_consequences.is_empty() {
        return;
    }
    let mut best_per_gene: HashMap<Arc<str>, usize> = HashMap::new();
    for (i, tc) in variant.transcript_consequences.iter().enumerate() {
        let gene = &tc.gene_id;
        match best_per_gene.get(gene) {
            Some(&prev_idx) => {
                if pick_score(tc) < pick_score(&variant.transcript_consequences[prev_idx]) {
                    best_per_gene.insert(gene.clone(), i);
                }
            }
            None => {
                best_per_gene.insert(gene.clone(), i);
            }
        }
    }
    let mut keep_indices: Vec<usize> = best_per_gene.into_values().collect();
    keep_indices.sort_unstable();
    // Compact in-place by swapping kept elements to the front, avoiding clones.
    for (dest, &src) in keep_indices.iter().enumerate() {
        if dest != src {
            variant.transcript_consequences.swap(dest, src);
        }
    }
    variant.transcript_consequences.truncate(keep_indices.len());
}

/// --pick_allele_gene: Keep only the best per (allele, gene) combination.
/// Multi-allelic variants are split, so this is equivalent to --per_gene.
fn apply_pick_allele_gene(variant: &mut InputVariant) {
    apply_per_gene(variant);
}

/// Append a flag to a TranscriptConsequence's `Arc<[String]>` flags.
fn push_flag(tc: &mut TranscriptConsequence, flag: &str) {
    let mut v: Vec<String> = tc.flags.to_vec();
    v.push(flag.to_string());
    tc.flags = v.into();
}

/// --flag_pick: Don't remove consequences, but set a PICK flag on the best one.
fn flag_pick(variant: &mut InputVariant) {
    if variant.transcript_consequences.is_empty() {
        return;
    }
    let best_idx = find_best_index(&variant.transcript_consequences);
    push_flag(&mut variant.transcript_consequences[best_idx], "PICK");
}

/// --flag_pick_allele: Flag the best consequence per allele.
/// Multi-allelic variants are split, so this is the same as --flag_pick.
fn flag_pick_allele(variant: &mut InputVariant) {
    flag_pick(variant);
}

/// --flag_pick_allele_gene: Flag the best consequence per (allele, gene).
fn flag_pick_allele_gene(variant: &mut InputVariant) {
    if variant.transcript_consequences.is_empty() {
        return;
    }
    let mut best_per_gene: HashMap<Arc<str>, usize> = HashMap::new();
    for (i, tc) in variant.transcript_consequences.iter().enumerate() {
        let gene = &tc.gene_id;
        match best_per_gene.get(gene) {
            Some(&prev_idx) => {
                if pick_score(tc) < pick_score(&variant.transcript_consequences[prev_idx]) {
                    best_per_gene.insert(gene.clone(), i);
                }
            }
            None => {
                best_per_gene.insert(gene.clone(), i);
            }
        }
    }
    for &idx in best_per_gene.values() {
        push_flag(&mut variant.transcript_consequences[idx], "PICK");
    }
}

/// Find the index of the best-scoring transcript consequence.
fn find_best_index(tcs: &[TranscriptConsequence]) -> usize {
    tcs.iter()
        .enumerate()
        .min_by_key(|(_, tc)| pick_score(tc))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::{
        Consequence, ConsequenceList, FeatureType, Impact, TranscriptConsequence,
    };
    use vep_core::variant::InputVariant;

    fn make_tc(
        transcript_id: &str,
        gene_id: &str,
        canonical: bool,
        biotype: &str,
        consequences: ConsequenceList,
    ) -> TranscriptConsequence {
        let impact = consequences
            .first()
            .map(|c| c.impact())
            .unwrap_or(Impact::MODIFIER);
        TranscriptConsequence {
            transcript_id: transcript_id.into(),
            gene_id: gene_id.into(),
            canonical,
            biotype: Some(biotype.into()),
            consequences,
            impact,
            feature_type: FeatureType::Transcript,
            strand: 1,
            ..Default::default()
        }
    }

    fn make_variant_with_tcs(tcs: Vec<TranscriptConsequence>) -> InputVariant {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.transcript_consequences = tcs;
        v
    }

    #[test]
    fn test_pick_canonical_protein_coding_wins() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "lncRNA",
            smallvec![Consequence::MissenseVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc2, tc1]);
        apply_pick(&mut v);
        assert_eq!(v.transcript_consequences.len(), 1);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST00000001");
        assert!(v.transcript_consequences[0].canonical);
    }

    #[test]
    fn test_pick_severity_breaks_tie() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            false,
            "protein_coding",
            smallvec![Consequence::IntronVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc1, tc2]);
        apply_pick(&mut v);
        assert_eq!(v.transcript_consequences.len(), 1);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST00000002");
    }

    #[test]
    fn test_most_severe() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "protein_coding",
            smallvec![Consequence::IntronVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc1, tc2]);
        apply_most_severe(&mut v);
        assert!(v.transcript_consequences.is_empty());
        assert_eq!(
            v.most_severe_consequence,
            Some(Consequence::MissenseVariant)
        );
    }

    #[test]
    fn test_summary() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "protein_coding",
            smallvec![Consequence::IntronVariant, Consequence::MissenseVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc1, tc2]);
        apply_summary(&mut v);
        assert_eq!(v.transcript_consequences.len(), 1);
        let csqs = &v.transcript_consequences[0].consequences;
        assert_eq!(csqs.len(), 2);
        assert!(csqs.contains(&Consequence::MissenseVariant));
        assert!(csqs.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_per_gene() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "protein_coding",
            smallvec![Consequence::IntronVariant],
        );
        let tc3 = make_tc(
            "ENST00000003",
            "ENSG00000002",
            true,
            "protein_coding",
            smallvec![Consequence::SynonymousVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc1, tc2, tc3]);
        apply_per_gene(&mut v);
        assert_eq!(v.transcript_consequences.len(), 2);

        let gene_ids: Vec<&str> = v
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.gene_id)
            .collect();
        assert!(gene_ids.contains(&"ENSG00000001"));
        assert!(gene_ids.contains(&"ENSG00000002"));

        let gene1_tc = v
            .transcript_consequences
            .iter()
            .find(|tc| &*tc.gene_id == "ENSG00000001")
            .unwrap();
        assert_eq!(&*gene1_tc.transcript_id, "ENST00000001");
    }

    #[test]
    fn test_flag_pick() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "lncRNA",
            smallvec![Consequence::IntronVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc1, tc2]);
        flag_pick(&mut v);
        assert_eq!(v.transcript_consequences.len(), 2);
        assert!(v.transcript_consequences[0]
            .flags
            .contains(&"PICK".to_string()));
        assert!(!v.transcript_consequences[1]
            .flags
            .contains(&"PICK".to_string()));
    }

    #[test]
    fn test_flag_pick_allele_gene() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "protein_coding",
            smallvec![Consequence::IntronVariant],
        );
        let tc3 = make_tc(
            "ENST00000003",
            "ENSG00000002",
            true,
            "protein_coding",
            smallvec![Consequence::SynonymousVariant],
        );
        let mut v = make_variant_with_tcs(vec![tc1, tc2, tc3]);
        flag_pick_allele_gene(&mut v);
        assert_eq!(v.transcript_consequences.len(), 3);
        let tc1 = &v.transcript_consequences[0];
        assert!(tc1.flags.contains(&"PICK".to_string()));
        let tc2 = &v.transcript_consequences[1];
        assert!(!tc2.flags.contains(&"PICK".to_string()));
        let tc3 = &v.transcript_consequences[2];
        assert!(tc3.flags.contains(&"PICK".to_string()));
    }

    #[test]
    fn test_filter_config_any_active() {
        let none_active = FilterConfig {
            pick: false,
            pick_allele: false,
            per_gene: false,
            pick_allele_gene: false,
            most_severe: false,
            summary: false,
            flag_pick: false,
            flag_pick_allele: false,
            flag_pick_allele_gene: false,
        };
        assert!(!none_active.any_active());

        let one_active = FilterConfig {
            pick: true,
            ..FilterConfig {
                pick: false,
                pick_allele: false,
                per_gene: false,
                pick_allele_gene: false,
                most_severe: false,
                summary: false,
                flag_pick: false,
                flag_pick_allele: false,
                flag_pick_allele_gene: false,
            }
        };
        assert!(one_active.any_active());
    }

    #[test]
    fn test_apply_filters_dispatches_correctly() {
        let tc1 = make_tc(
            "ENST00000001",
            "ENSG00000001",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
        );
        let tc2 = make_tc(
            "ENST00000002",
            "ENSG00000001",
            false,
            "lncRNA",
            smallvec![Consequence::IntronVariant],
        );

        let mut v = make_variant_with_tcs(vec![tc1.clone(), tc2.clone()]);
        let config = FilterConfig {
            pick: true,
            pick_allele: false,
            per_gene: false,
            pick_allele_gene: false,
            most_severe: false,
            summary: false,
            flag_pick: false,
            flag_pick_allele: false,
            flag_pick_allele_gene: false,
        };
        apply_filters(&mut v, &config);
        assert_eq!(v.transcript_consequences.len(), 1);

        let mut v = make_variant_with_tcs(vec![tc1.clone(), tc2.clone()]);
        let config = FilterConfig {
            pick: true,
            pick_allele: false,
            per_gene: false,
            pick_allele_gene: false,
            most_severe: true,
            summary: false,
            flag_pick: false,
            flag_pick_allele: false,
            flag_pick_allele_gene: false,
        };
        apply_filters(&mut v, &config);
        assert!(v.transcript_consequences.is_empty());
        assert_eq!(
            v.most_severe_consequence,
            Some(Consequence::MissenseVariant)
        );
    }

    #[test]
    fn test_empty_consequences() {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        apply_pick(&mut v);
        assert!(v.transcript_consequences.is_empty());

        apply_per_gene(&mut v);
        assert!(v.transcript_consequences.is_empty());

        flag_pick(&mut v);
        assert!(v.transcript_consequences.is_empty());
    }

    /// Helper: build a TC with all pick-relevant fields set explicitly.
    /// Defaults match Perl's "missing-attribute" defaults so unset fields
    /// score as "worst" in their slot.
    #[allow(clippy::too_many_arguments)]
    fn make_tc_pick(
        transcript_id: &str,
        canonical: bool,
        biotype: &str,
        consequences: ConsequenceList,
        mane_select: Option<&str>,
        mane_plus_clinical: Option<&str>,
        appris: Option<&str>,
        tsl: Option<u8>,
        ccds: Option<&str>,
    ) -> TranscriptConsequence {
        let impact = consequences
            .first()
            .map(|c| c.impact())
            .unwrap_or(Impact::MODIFIER);
        TranscriptConsequence {
            transcript_id: transcript_id.into(),
            gene_id: "ENSG00000001".into(),
            canonical,
            biotype: Some(biotype.into()),
            consequences,
            impact,
            feature_type: FeatureType::Transcript,
            strand: 1,
            mane_select: mane_select.map(String::from),
            mane_plus_clinical: mane_plus_clinical.map(String::from),
            appris: appris.map(String::from),
            tsl,
            ccds: ccds.map(String::from),
            ..Default::default()
        }
    }

    #[test]
    fn test_appris_rank_principal() {
        assert_eq!(appris_rank(Some("principal1")), 1);
        assert_eq!(appris_rank(Some("principal2")), 2);
        assert_eq!(appris_rank(Some("principal9")), 9);
    }

    #[test]
    fn test_appris_rank_alternative() {
        assert_eq!(appris_rank(Some("alternative1")), 11);
        assert_eq!(appris_rank(Some("alternative2")), 12);
        assert_eq!(appris_rank(Some("alternative9")), 19);
    }

    #[test]
    fn test_appris_rank_missing_or_unparseable() {
        assert_eq!(appris_rank(None), 100);
        assert_eq!(appris_rank(Some("")), 100);
        assert_eq!(appris_rank(Some("garbage")), 100);
        // No digit → treated as missing (Perl regex requires (\d+))
        assert_eq!(appris_rank(Some("principal")), 100);
    }

    #[test]
    fn test_appris_rank_short_form() {
        // Perl regex `/([A-Za-z]).+(\d+)/` requires a separator, so the short
        // form "P1" does not match and scores 100.
        assert_eq!(appris_rank(Some("P1")), 100);
        // "p_1" exercises the regex's `.+` (one or more chars between letter
        // and digit). Lowercase 'p' is principal; expected rank = 1.
        assert_eq!(appris_rank(Some("p_1")), 1);
    }

    /// The APPRIS tier follows `m/([A-Za-z]).+(\d+)/`, `$grade += 10` for a
    /// lowercase `a` and `if $grade` (OutputFactory.pm:759-767 at
    /// release/116.2), run on the cache's own spellings: `principalM` and
    /// `alternativeM` carry no digit and keep the default 100; the digit
    /// captured is the value's last, alone; the bump is case-sensitive; a zero
    /// grade is kept only once bumped.
    #[test]
    fn appris_rank_principal_m_is_100() {
        assert_eq!(appris_rank(Some("principalM")), 100);
        assert_eq!(appris_rank(Some("alternativeM")), 100);
        assert_eq!(appris_rank(Some("principal12")), 2);
        assert_eq!(appris_rank(Some("Alternative3")), 3);
        assert_eq!(appris_rank(Some("alternative0")), 10);
        assert_eq!(appris_rank(Some("principal0")), 100);
    }

    /// The TSL tier reads the digits the loader took from the `TSL` attribute
    /// with `m/tsl(\d+)/` and keeps only when non-zero (OutputFactory.pm:750-754
    /// at release/116.2): `tsl1 (assigned to previous version 8)` loads as 1
    /// and scores 1, `tslNA` loads as none and scores the default 100, as does
    /// a zero.
    #[test]
    fn tsl_tier_reads_attribute_digits() {
        let scored = |tsl: Option<u8>| {
            pick_score(&make_tc_pick(
                "ENST00000000001",
                false,
                "protein_coding",
                smallvec![Consequence::MissenseVariant],
                None,
                None,
                None,
                tsl,
                None,
            ))
            .4
        };
        assert_eq!(scored(Some(1)), 1);
        assert_eq!(scored(Some(5)), 5);
        assert_eq!(scored(None), 100);
        assert_eq!(scored(Some(0)), 100);
    }

    #[test]
    fn test_pick_mane_select_beats_canonical() {
        // A non-canonical transcript with a MANE Select tag must win over a
        // canonical non-MANE transcript: Perl pick_order slot 0 precedes slot 2.
        let tc_canonical = make_tc_pick(
            "ENST_canonical",
            true, // canonical
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None, // no MANE_Select
            None, // no MANE_Plus_Clinical
            None,
            None,
            None,
        );
        let tc_mane = make_tc_pick(
            "ENST_mane_select",
            false, // NOT canonical
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            Some("NM_000001.3"), // has MANE_Select
            None,
            None,
            None,
            None,
        );
        let mut v = make_variant_with_tcs(vec![tc_canonical, tc_mane]);
        apply_pick(&mut v);
        assert_eq!(v.transcript_consequences.len(), 1);
        assert_eq!(
            &*v.transcript_consequences[0].transcript_id, "ENST_mane_select",
            "MANE_Select must win over canonical"
        );
    }

    #[test]
    fn test_pick_mane_plus_clinical_beats_canonical() {
        // MANE Plus Clinical, like MANE Select, must beat canonical-only.
        let tc_canonical = make_tc_pick(
            "ENST_canonical",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            None,
            None,
        );
        let tc_mpc = make_tc_pick(
            "ENST_mpc",
            false,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            Some("NM_000002.3"),
            None,
            None,
            None,
        );
        let mut v = make_variant_with_tcs(vec![tc_canonical, tc_mpc]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_mpc");
    }

    #[test]
    fn test_pick_mane_select_beats_mane_plus_clinical() {
        // When both flags are present on different transcripts, MANE_Select
        // (slot 0) wins over MANE_Plus_Clinical (slot 1) per Perl pick_order.
        let tc_mpc = make_tc_pick(
            "ENST_mpc",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            Some("NM_002.1"),
            None,
            None,
            None,
        );
        let tc_mane = make_tc_pick(
            "ENST_mane",
            false,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            Some("NM_001.1"),
            None,
            None,
            None,
            None,
        );
        let mut v = make_variant_with_tcs(vec![tc_mpc, tc_mane]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_mane");
    }

    #[test]
    fn test_pick_appris_p1_beats_p2() {
        // Two canonical, biotype-matched, no-MANE TCs differing only in APPRIS.
        // principal1 (rank=1) must beat principal2 (rank=2).
        let tc_p2 = make_tc_pick(
            "ENST_p2",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            Some("principal2"),
            None,
            None,
        );
        let tc_p1 = make_tc_pick(
            "ENST_p1",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            Some("principal1"),
            None,
            None,
        );
        let mut v = make_variant_with_tcs(vec![tc_p2, tc_p1]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_p1");
    }

    #[test]
    fn test_pick_appris_principal_beats_alternative() {
        // principal{N} ranks 1-9, alternative{N} ranks 11-19, so any principal
        // beats any alternative (including principal9 vs alternative1).
        let tc_alt1 = make_tc_pick(
            "ENST_alt1",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            Some("alternative1"),
            None,
            None,
        );
        let tc_p9 = make_tc_pick(
            "ENST_p9",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            Some("principal9"),
            None,
            None,
        );
        let mut v = make_variant_with_tcs(vec![tc_alt1, tc_p9]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_p9");
    }

    #[test]
    fn test_pick_tsl_1_beats_tsl_5() {
        // Tied through canonical/biotype/appris; TSL 1 must beat TSL 5.
        let tc_t5 = make_tc_pick(
            "ENST_t5",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            Some(5),
            None,
        );
        let tc_t1 = make_tc_pick(
            "ENST_t1",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            Some(1),
            None,
        );
        let mut v = make_variant_with_tcs(vec![tc_t5, tc_t1]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_t1");
    }

    #[test]
    fn test_pick_ccds_present_beats_absent() {
        let tc_no_ccds = make_tc_pick(
            "ENST_no_ccds",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            None,
            None,
        );
        let tc_ccds = make_tc_pick(
            "ENST_ccds",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            None,
            Some("CCDS123.1"),
        );
        let mut v = make_variant_with_tcs(vec![tc_no_ccds, tc_ccds]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_ccds");
    }

    #[test]
    fn test_pick_ccds_dash_treated_as_absent() {
        // Perl checks `_ccds && _ccds ne '-'`: a literal "-" placeholder means
        // "no CCDS" and scores 1 (worst), same as None.
        let tc_dash = make_tc_pick(
            "ENST_dash",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            None,
            Some("-"),
        );
        let tc_real_ccds = make_tc_pick(
            "ENST_real",
            false,
            "protein_coding", // intentionally non-canonical
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            None,
            None,
            Some("CCDS456.2"),
        );
        let mut v = make_variant_with_tcs(vec![tc_dash, tc_real_ccds]);
        apply_pick(&mut v);
        // The canonical slot is ordered ahead of the CCDS slot, and only the
        // dash-CCDS transcript is canonical, so it wins on slot 2 before CCDS is
        // ever compared. Missing CCDS does not cost it the pick.
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_dash");
    }

    #[test]
    fn test_pick_full_priority_chain() {
        // Multi-tier tiebreak: MANE_Select absent on both, MANE_Plus_Clinical
        // absent on both, canonical tied (both true), APPRIS P1 vs P2.
        // P1 must win on slot 3.
        let tc_p2_canonical = make_tc_pick(
            "ENST_p2",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            Some("principal2"),
            Some(1),
            Some("CCDS1.1"),
        );
        let tc_p1_canonical = make_tc_pick(
            "ENST_p1",
            true,
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,
            None,
            Some("principal1"),
            Some(1),
            Some("CCDS1.1"),
        );
        let mut v = make_variant_with_tcs(vec![tc_p2_canonical, tc_p1_canonical]);
        apply_pick(&mut v);
        assert_eq!(&*v.transcript_consequences[0].transcript_id, "ENST_p1");
    }

    #[test]
    fn test_pick_score_slot_independence() {
        // Smoke-test that the 9-tuple comparator orders strictly by slot.
        // A TC with worse slot 0 (mane_select=1) but better slots 1-7 still
        // loses to a TC with mane_select=0 even if everything else is worst.
        let tc_mane_only = make_tc_pick(
            "ENST_mane",
            false,                                 // non-canonical
            "lncRNA",                              // non-protein-coding
            smallvec![Consequence::IntronVariant], // worst rank
            Some("NM_001"),
            None,
            None,
            None,
            None,
        );
        let tc_everything_else = make_tc_pick(
            "ENST_other",
            true, // canonical
            "protein_coding",
            smallvec![Consequence::MissenseVariant],
            None,           // no MANE_Select → score 1
            Some("NM_002"), // has MANE_Plus_Clinical → score 0
            Some("principal1"),
            Some(1),
            Some("CCDS9.9"),
        );
        let mut v = make_variant_with_tcs(vec![tc_everything_else, tc_mane_only]);
        apply_pick(&mut v);
        assert_eq!(
            &*v.transcript_consequences[0].transcript_id, "ENST_mane",
            "MANE_Select (slot 0) must beat all later slots combined"
        );
    }
}
