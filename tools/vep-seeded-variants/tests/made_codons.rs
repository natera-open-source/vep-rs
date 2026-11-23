// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A genetic code that translates a codon unlike the standard code gets up to three
//! records whose substitution makes such a codon, whether or not a coding sequence of the
//! cache holds one: one keeps the residue, one changes it and one stands in the closing
//! stop codon. A code gets each of the three that a substitution can serve under it.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{
    bacterium, category_of, copy_tree, mitochondrion, mitochondrion_under, pair, pilot_cache,
    records_of, run_to_bytes, summary_line, Scratch,
};
use vep_seeded_variants::generate::{Generated, Records};

/// The first and the last base of the made-up cache's one coding sequence: twenty
/// codons in one exon, the whole transcript, under genetic code 4.
const FIRST_BASE: u64 = 101;
const LAST_BASE: u64 = 160;

/// The three records by the rule alone. Code 4 reads one codon unlike the standard code.
/// The coding sequence does not hold it, and three of its codons are one substitution
/// from it: codon 5, counted from 0, by its third base and with the same residue; codon
/// 11 by its second base and with another residue; and codon 19, the stop codon, by its
/// second base. Each entry is the category's name after the code, the codon, the base's
/// place in it, and the base before and after as the transcript reads them.
const MADE: [(&str, u64, u64, u8, u8); 3] = [
    ("made_codon_synonymous", 5, 2, b'G', b'A'),
    ("made_codon_missense", 11, 1, b'T', b'G'),
    ("made_codon_stop", 19, 1, b'A', b'G'),
];

/// What a category says when no substitution of the cache makes such a codon.
const NO_SUBSTITUTION: &str = "no substitution makes a codon this code reads differently";

/// A record as its ID, its position and its two alleles.
type Line = (String, u64, String, String);

/// The three records where the coding sequence is read on the forward strand, or on the
/// reverse strand over the same bases, in the output's order.
fn expected(reverse: bool) -> Vec<Line> {
    let mut lines: Vec<Line> = MADE
        .iter()
        .map(|(kind, codon, offset, before, after)| {
            let into = 3 * codon + offset;
            let (position, before, after) = if reverse {
                (LAST_BASE - into, pair(*before), pair(*after))
            } else {
                (FIRST_BASE + into, *before, *after)
            };
            (
                format!("code4_{kind}_1"),
                position,
                char::from(before).to_string(),
                char::from(after).to_string(),
            )
        })
        .collect();
    lines.sort_by_key(|line| line.1);
    lines
}

/// The records of the three categories, in the output's order.
fn made(generated: &Generated) -> Vec<Line> {
    generated
        .records
        .iter()
        .filter(|record| record.id.contains("made_codon"))
        .map(|record| {
            (
                record.id.clone(),
                record.position,
                record.reference.clone(),
                record.alt.clone(),
            )
        })
        .collect()
}

/// A copy of the made-up cache in `scratch` whose shard is `edit`'s reading of it.
fn edited(scratch: &Scratch, edit: impl Fn(String) -> String) -> PathBuf {
    let copy = scratch.path().join("cache");
    copy_tree(&mitochondrion(), &copy, false);
    let shard = copy.join("transcripts/MT_MADE/1-1000000.json");
    let text = fs::read_to_string(&shard).expect("the shard");
    let changed = edit(text.clone());
    assert_ne!(changed, text, "the edit changes nothing");
    fs::write(&shard, changed).expect("a written file");
    copy
}

/// The shard with `rewrite` applied to the bases of its coding sequence.
fn with_coding(text: String, rewrite: impl Fn(&mut [u8])) -> String {
    let key = r#""translateable_seq": ""#;
    let start = text.find(key).expect("the coding sequence") + key.len();
    let length = text[start..].find('"').expect("its end");
    assert_eq!(length as u64, LAST_BASE - FIRST_BASE + 1);
    let mut bytes = text.into_bytes();
    rewrite(&mut bytes[start..start + length]);
    String::from_utf8(bytes).expect("text")
}

/// Writes codon `from`'s bases over codon `to`.
fn copy_codon(coding: &mut [u8], from: usize, to: usize) {
    coding.copy_within(3 * from..3 * from + 3, 3 * to);
}

#[test]
fn a_coding_sequence_without_such_a_codon_gets_the_three_records_that_make_one() {
    let scratch = Scratch::new("made");
    let (records, summary) = run_to_bytes(&mitochondrion(), &scratch, &["--seed", "1"]);
    let pinned_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/expected");
    let pinned = [
        (pinned_dir.join("made-up-mitochondrion.vcf"), records),
        (
            pinned_dir.join("made-up-mitochondrion.summary.tsv"),
            summary,
        ),
    ];
    let update = std::env::var_os("VEP_SEEDED_VARIANTS_UPDATE").is_some();
    for (file, written) in pinned {
        if update {
            fs::write(&file, &written).expect("a written file");
            continue;
        }
        let expected = fs::read(&file).unwrap_or_default();
        assert!(
            written == expected,
            "the run on the made-up cache under code 4 does not write {}:\n{}",
            file.display(),
            String::from_utf8_lossy(&written)
        );
    }

    // One substitution alone serves each category, so no seed gives another record.
    for seed in 1..=10 {
        let light = records_of(&mitochondrion(), seed, Records::Light);
        assert_eq!(made(&light), expected(false), "seed {seed}");
        // The coding sequence holds no codon the code reads differently.
        assert_eq!(
            summary_line(&light, "code4_changed_codon"),
            (1, 0, "no codon this code reads differently"),
            "seed {seed}"
        );
    }
    // One transcript gives a category one candidate, whatever the mode asks for.
    let wide = records_of(&mitochondrion(), 1, Records::Wide(12));
    assert_eq!(made(&wide), expected(false));
    for (kind, ..) in MADE {
        assert_eq!(
            summary_line(&wide, &format!("code4_{kind}")),
            (12, 1, "few candidates"),
            "{kind}"
        );
    }
}

#[test]
fn a_coding_sequence_on_the_reverse_strand_gets_the_same_three_substitutions() {
    let scratch = Scratch::new("made-reverse");
    let cache = edited(&scratch, |text| {
        assert_eq!(text.matches(r#""strand": 1,"#).count(), 1);
        text.replace(r#""strand": 1,"#, r#""strand": -1,"#)
    });
    for seed in 1..=5 {
        let light = records_of(&cache, seed, Records::Light);
        assert_eq!(made(&light), expected(true), "seed {seed}");
    }
}

#[test]
fn a_category_no_substitution_serves_writes_no_record_and_says_why() {
    // Codon 1 is two substitutions or more from the codon code 4 reads differently.
    // Written over codon 5 it leaves no substitution that keeps a residue, and over
    // codon 11 none that changes one.
    for (replaced, lost) in [(5, "made_codon_synonymous"), (11, "made_codon_missense")] {
        let scratch = Scratch::new("made-none");
        let cache = edited(&scratch, |text| {
            with_coding(text, |coding| copy_codon(coding, 1, replaced))
        });
        let light = records_of(&cache, 1, Records::Light);
        let kept: Vec<Line> = expected(false)
            .into_iter()
            .filter(|line| line.0 != format!("code4_{lost}_1"))
            .collect();
        assert_eq!(made(&light), kept, "{lost}");
        for (kind, ..) in MADE {
            let line = if kind == lost {
                (1, 0, NO_SUBSTITUTION)
            } else {
                (1, 1, "-")
            };
            assert_eq!(
                summary_line(&light, &format!("code4_{kind}")),
                line,
                "{lost}: {kind}"
            );
        }
    }

    // The other stop codon of code 4 is two substitutions from that codon: the closing
    // stop codon is still lost by a substitution, and none makes the codon there.
    let scratch = Scratch::new("made-no-stop");
    let cache = edited(&scratch, |text| {
        with_coding(text, |coding| {
            let last = coding.len() - 1;
            assert_eq!(coding[last], b'A');
            coding[last] = b'G';
        })
    });
    let light = records_of(&cache, 1, Records::Light);
    assert_eq!(made(&light), &expected(false)[..2]);
    assert_eq!(
        summary_line(&light, "code4_made_codon_stop"),
        (1, 0, NO_SUBSTITUTION)
    );
    assert_eq!(summary_line(&light, "code4_stop_lost"), (1, 1, "-"));

    // Codon 11 written over the last codon closes the coding sequence on a sense codon
    // one substitution from that codon: the stop category has no place there.
    let scratch = Scratch::new("made-open-end");
    let cache = edited(&scratch, |text| {
        with_coding(text, |coding| copy_codon(coding, 11, 19))
    });
    let light = records_of(&cache, 1, Records::Light);
    assert_eq!(made(&light), &expected(false)[..2]);
    assert_eq!(
        summary_line(&light, "code4_made_codon_stop"),
        (1, 0, NO_SUBSTITUTION)
    );
    assert_eq!(
        summary_line(&light, "code4_stop_lost"),
        (1, 0, "no complete stop codon")
    );

    // A coding sequence that opens inside a codon is not read: the reason is the one
    // the start-codon category gives.
    let scratch = Scratch::new("made-unread");
    let cache = edited(&scratch, |text| {
        assert_eq!(text.matches(r#""start_phase": 0"#).count(), 1);
        text.replace(r#""start_phase": 0"#, r#""start_phase": 1"#)
    });
    let light = records_of(&cache, 1, Records::Light);
    assert!(made(&light).is_empty());
    for kind in [
        "start_codon",
        "made_codon_synonymous",
        "made_codon_missense",
        "made_codon_stop",
    ] {
        assert_eq!(
            summary_line(&light, &format!("code4_{kind}")),
            (1, 0, "no complete coding sequence"),
            "{kind}"
        );
    }
}

/// The categories of the three that the summary lists, in its order.
fn listed(generated: &Generated) -> Vec<&str> {
    generated
        .summary
        .iter()
        .filter(|line| line.category.contains("made_codon"))
        .map(|line| line.category.as_str())
        .collect()
}

#[test]
fn a_code_has_the_categories_a_substitution_can_serve_under_it() {
    let scratch = Scratch::new("made-codes");
    // Code 12 reads one leucine codon of the standard code as serine. No serine codon
    // and no stop codon is one substitution from it, so the code has one category.
    let twelve = records_of(&mitochondrion_under(&scratch, 12), 1, Records::Light);
    assert_eq!(listed(&twelve), ["code12_made_codon_missense"]);
    assert!(made(&twelve)
        .iter()
        .all(|line| line.0 == "code12_made_codon_missense_1"));

    // Code 23 reads one leucine codon of the standard code as a stop codon, so no
    // substitution that makes it ends on a sense codon. The closing stop codon of the
    // made-up coding sequence becomes it by its second base, and stays a stop codon.
    let other = records_of(&mitochondrion_under(&scratch, 23), 1, Records::Light);
    assert_eq!(listed(&other), ["code23_made_codon_stop"]);
    assert_eq!(
        made(&other),
        [(
            "code23_made_codon_stop_1".to_owned(),
            FIRST_BASE + 3 * 19 + 1,
            "A".to_owned(),
            "T".to_owned()
        )]
    );
}

#[test]
fn a_code_that_translates_every_codon_as_the_standard_code_does_has_no_such_category() {
    // The made-up bacterium names code 11 alone, which differs from the standard code
    // in its start codons; the mouse slice names the standard code and code 2.
    let bacterial = records_of(&bacterium(), 1, Records::Wide(12));
    assert!(made(&bacterial).is_empty());
    assert!(bacterial
        .summary
        .iter()
        .all(|line| !line.category.contains("made_codon")));

    let mouse = records_of(&pilot_cache("GRCm39-mouse"), 1, Records::Light);
    assert_eq!(
        listed(&mouse),
        [
            "code2_made_codon_synonymous",
            "code2_made_codon_missense",
            "code2_made_codon_stop"
        ]
    );
    assert!(made(&mouse)
        .iter()
        .all(|line| category_of(&line.0).starts_with("code2_made_codon_")));
}
