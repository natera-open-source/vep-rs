// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up outputs of four `<CNV:TR>` records at one position, whose
//! rows share one key: Ensembl VEP names `feature_elongation` on every row, and vep-rs
//! `feature_truncation` on the rows of the alleles that state a loss and neither term on
//! the row of the allele as long as its run, as `cnv_tr_symbolic_loss_named_a_gain`
//! documents. vep-rs writes Ensembl VEP's set on the key as well, so the entry names the
//! records whose rows differ. The nearest keys (another Allele, another term that differs, a
//! second set beside the class's, the terms the other way round) stay unexplained, and a
//! key whose rows agree writes no entry.

mod classify_support;
mod common;

use classify_support::{counted, divergence, key, obj, record, records, row, written, Case};
use common::clean;
use vep_golden_corpus::classify::tandem_repeat_loss::CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN as CLASS;
use vep_golden_corpus::constants::UNEXPLAINED;

const LOCATION: &str = "7:1001-1009";
const FEATURE: &str = "ENST90000000001";

/// The terms of a gain and of a loss as each engine writes them, and as `classify` writes a
/// set: the terms sorted.
const GAIN: &str = "feature_elongation,coding_sequence_variant";
const LOSS: &str = "feature_truncation,coding_sequence_variant";
const GAIN_SET: &str = "coding_sequence_variant,feature_elongation";
const LOSS_SET: &str = "coding_sequence_variant,feature_truncation";

/// The terms vep-rs writes for an allele as long as its run, a set of one term.
const SAME: &str = "coding_sequence_variant";

/// The four records both outputs annotate: a gain, a loss, a record of two alleles whose
/// second states a loss, and an allele as long as the run.
fn held() -> String {
    let line = |name: &str, alleles: &str, units: &str, counts: &str| {
        format!("7\t1000\t{name}\tG\t{alleles}\t.\t.\tSVLEN=9;RUS={units};RUC={counts}")
    };
    records(&[
        record(&line("tr_gain", "<CNV:TR>", "CAG", "4")),
        record(&line("tr_loss", "<CNV:TR>", "CAG", "2")),
        record(&line("tr_multi", "<CNV:TR>,<CNV:TR>", "CAG,CAG", "4,2")),
        record(&line("tr_same", "<CNV:TR>", "CAG", "3")),
    ])
}

/// A default-format output of the four records on one transcript: the rows of `tr_gain`,
/// of `tr_loss`, of the two alleles of `tr_multi` and of `tr_same`, each with its terms.
fn output(allele: &str, terms: [&str; 5]) -> String {
    let names = ["tr_gain", "tr_loss", "tr_multi", "tr_multi", "tr_same"];
    let rows = names.iter().zip(terms);
    rows.map(|(name, terms)| row(name, LOCATION, allele, FEATURE, terms))
        .collect()
}

/// Runs `classify` on two outputs of the four records and returns the manifest's lists, or
/// `None` where the run writes no entry.
fn entry_of(tag: &str, allele: &str, vep: [&str; 5], vep_rs: [&str; 5]) -> Option<String> {
    let case = Case::new(tag);
    let held = held();
    case.manifest(&obj(&[&held]));
    case.outputs(&output(allele, vep), &output(allele, vep_rs));
    let ran = case.classify(&[]);
    if ran == clean(0, &counted((0, "{}"), (0, "{}"), "[]")) {
        assert_eq!(
            case.after(),
            obj(&[&held, &written(&[], &[], ("{}", "{}"))])
        );
        return None;
    }
    for class in [CLASS, UNEXPLAINED] {
        let summary = format!("{{'{class}': 1}}");
        if ran == clean(0, &counted((1, &summary), (0, "{}"), "[]")) {
            let after = case.after();
            // A second run over the same outputs writes the same bytes.
            assert_eq!(case.classify(&[]), ran, "{tag}");
            assert_eq!(case.after(), after, "{tag}");
            return Some(after.replace(&held, "<records>"));
        }
    }
    panic!("{tag}: classify ended {ran:?}");
}

/// The manifest's lists with one entry on the key of `allele`: the records `ordinals`,
/// Ensembl VEP's set, vep-rs's sets and the class.
fn one_entry(
    allele: &str,
    ordinals: &[usize],
    vep_set: &str,
    vep_rs_sets: &[&str],
    class: &str,
) -> Option<String> {
    let expected = written(
        &[divergence(
            &key(LOCATION, allele, (FEATURE, "Transcript"), ordinals),
            vep_set,
            vep_rs_sets,
            class,
        )],
        &[],
        (&format!(r#"{{"{class}": 1}}"#), "{}"),
    );
    Some(obj(&["<records>", &expected]))
}

#[test]
fn the_rows_of_a_loss_on_a_key_the_gains_share_take_the_class() {
    // The loss's own record and the second allele of the record of two: the entry names
    // both records, and vep-rs's sets on the key hold Ensembl VEP's beside the loss's.
    assert_eq!(
        entry_of(
            "takes",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, LOSS, GAIN, LOSS, GAIN]
        ),
        one_entry(
            "tandem_repeat",
            &[1, 2],
            GAIN_SET,
            &[GAIN_SET, LOSS_SET],
            CLASS
        )
    );
    // One record alone.
    assert_eq!(
        entry_of(
            "one_record",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, LOSS, GAIN, GAIN, GAIN]
        ),
        one_entry(
            "tandem_repeat",
            &[1],
            GAIN_SET,
            &[GAIN_SET, LOSS_SET],
            CLASS
        )
    );
}

#[test]
fn the_row_of_an_allele_as_long_as_its_run_on_a_key_the_gains_share_takes_the_class() {
    // The record's one row names neither term: the entry names the record, and vep-rs's
    // sets on the key hold Ensembl VEP's beside the set without the gain term.
    assert_eq!(
        entry_of(
            "as_long",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, GAIN, GAIN, GAIN, SAME]
        ),
        one_entry("tandem_repeat", &[3], GAIN_SET, &[SAME, GAIN_SET], CLASS)
    );
    // Beside the rows of the alleles that state a loss: one entry names every record whose
    // row differs, and both of vep-rs's own sets stand beside Ensembl VEP's.
    assert_eq!(
        entry_of(
            "as_long_beside_losses",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, LOSS, GAIN, LOSS, SAME]
        ),
        one_entry(
            "tandem_repeat",
            &[1, 2, 3],
            GAIN_SET,
            &[SAME, GAIN_SET, LOSS_SET],
            CLASS
        )
    );
}

#[test]
fn the_nearest_keys_stay_unexplained() {
    // Another Allele than a tandem repeat's.
    assert_eq!(
        entry_of(
            "another_allele",
            "duplication",
            [GAIN; 5],
            [GAIN, LOSS, GAIN, GAIN, GAIN]
        ),
        one_entry(
            "duplication",
            &[1],
            GAIN_SET,
            &[GAIN_SET, LOSS_SET],
            UNEXPLAINED
        )
    );
    // Another term that differs beside the two.
    let other = "feature_truncation,3_prime_UTR_variant";
    assert_eq!(
        entry_of(
            "another_term",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, other, GAIN, GAIN, GAIN]
        ),
        one_entry(
            "tandem_repeat",
            &[1],
            GAIN_SET,
            &["3_prime_UTR_variant,feature_truncation", GAIN_SET],
            UNEXPLAINED
        )
    );
    // The loss's set on one record and a set with another term that differs on a second: the
    // entry names both records and takes no class.
    let second = "feature_truncation,intron_variant";
    assert_eq!(
        entry_of(
            "a_second_set",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, LOSS, GAIN, second, GAIN]
        ),
        one_entry(
            "tandem_repeat",
            &[1, 2],
            GAIN_SET,
            &[GAIN_SET, LOSS_SET, second],
            UNEXPLAINED
        )
    );
    // The two terms the other way round: Ensembl VEP's set of the loss is one vep-rs writes
    // nowhere on the key, and the entry is the key's, with every record of it.
    assert_eq!(
        entry_of(
            "the_other_way_round",
            "tandem_repeat",
            [GAIN, LOSS, GAIN, GAIN, GAIN],
            [GAIN; 5]
        ),
        one_entry(
            "tandem_repeat",
            &[0, 1, 2, 3],
            LOSS_SET,
            &[GAIN_SET],
            UNEXPLAINED
        )
    );
}

#[test]
fn the_nearest_keys_of_an_allele_as_long_as_its_run_stay_unexplained() {
    // Another Allele than a tandem repeat's.
    assert_eq!(
        entry_of(
            "as_long_another_allele",
            "duplication",
            [GAIN; 5],
            [GAIN, GAIN, GAIN, GAIN, SAME]
        ),
        one_entry(
            "duplication",
            &[3],
            GAIN_SET,
            &[SAME, GAIN_SET],
            UNEXPLAINED
        )
    );
    // Another term in the place of the two vep-rs leaves out.
    let other = "3_prime_UTR_variant";
    assert_eq!(
        entry_of(
            "as_long_another_term",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, GAIN, GAIN, GAIN, other]
        ),
        one_entry(
            "tandem_repeat",
            &[3],
            GAIN_SET,
            &[other, GAIN_SET],
            UNEXPLAINED
        )
    );
    // The set without the gain term on one record and a set with another term that differs
    // on a second: the entry names both records and takes no class.
    let second = "feature_truncation,intron_variant";
    assert_eq!(
        entry_of(
            "as_long_beside_a_second_set",
            "tandem_repeat",
            [GAIN; 5],
            [GAIN, second, GAIN, GAIN, SAME]
        ),
        one_entry(
            "tandem_repeat",
            &[1, 3],
            GAIN_SET,
            &[SAME, GAIN_SET, second],
            UNEXPLAINED
        )
    );
    // The other way round: Ensembl VEP's set without the gain term is one vep-rs writes
    // nowhere on the key, and the entry is the key's, with every record of it.
    assert_eq!(
        entry_of(
            "as_long_the_other_way_round",
            "tandem_repeat",
            [GAIN, GAIN, GAIN, GAIN, SAME],
            [GAIN; 5]
        ),
        one_entry(
            "tandem_repeat",
            &[0, 1, 2, 3],
            SAME,
            &[GAIN_SET],
            UNEXPLAINED
        )
    );
}

#[test]
fn a_key_whose_rows_agree_writes_no_entry() {
    assert_eq!(
        entry_of("agree", "tandem_repeat", [GAIN; 5], [GAIN; 5]),
        None
    );
    assert_eq!(
        entry_of(
            "agree_on_a_loss",
            "tandem_repeat",
            [GAIN, LOSS, GAIN, LOSS, SAME],
            [GAIN, LOSS, GAIN, LOSS, SAME]
        ),
        None
    );
}
