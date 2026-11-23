// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The ninth shape of `hgvsp_stop_codon_protein_comparison` with more than one residue
//! before the inserted stop: the reference writes residues and a stop as inserted between
//! two neighbouring residues and the candidate writes the stop as unchanged. The row does
//! not name the reference residues past the right flank, so the class reads the coding
//! sequence of the row's transcript from the candidate's JSON cache and accepts the pair
//! where that sequence shows the protein kept. Without the cache the pair is a difference.

mod compare_common;
mod support;

use std::ffi::OsString;
use std::path::Path;
use std::path::PathBuf;

use compare_common::default_row;
use compare_common::default_text;
use compare_common::entry;
use support::compact;
use support::member;
use support::report;
use support::run;
use support::scratch;
use vep_formats_compare::error::Error;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;

const STOP_CODON_CLASS: &str = "hgvsp_stop_codon_protein_comparison";
const PROTEIN: &str = "ENSP00000000001.1";
const LOCATION: &str = "21:1009-1010";

/// A made-up transcript whose coding sequence is six codons and a stop codon: under the
/// standard code methionine, alanine, lysine, leucine, serine and valine, the stop at
/// residue 7.
const KEPT: &str = "ENST00000000301";
const KEPT_CODONS: [&str; 7] = ["ATG", "GCT", "AAA", "CTG", "TCT", "GTT", "TAA"];
/// The same six codons with no stop codon after them.
const UNCLOSED: &str = "ENST00000000302";
/// The seven codons under a `codon_table` no table is held for.
const UNHELD_CODE: &str = "ENST00000000303";
/// A transcript under the vertebrate mitochondrial code, whose fourth codon that code reads
/// as methionine and the standard code as isoleucine, and whose fifth it reads as
/// tryptophan and the standard code as a stop.
const MITOCHONDRIAL: &str = "ENST00000000304";
const MITOCHONDRIAL_CODONS: [&str; 7] = ["ATG", "GCT", "AAA", "ATA", "TGA", "GTT", "TAA"];
/// The seven codons and one more codon after them: the stop codon is not the last codon.
const STOP_INSIDE: &str = "ENST00000000305";
/// The seven codons and one base after the stop codon.
const PARTIAL_CODON: &str = "ENST00000000306";

/// Leucine, serine, valine and a stop inserted between lysine 3 and leucine 4 of `KEPT`,
/// against the stop unchanged at residue 7: by the first value the alternate protein is
/// the six residues of the reference protein and a stop at residue 7.
const RESIDUES_AND_STOP: &str = "Lys3_Leu4insLeuSerValTer";
const KEPT_STOP: &str = "Ter7=";

fn text(path: &Path) -> String {
    path.to_str().expect("a path of text").to_string()
}

fn cached(stable_id: &str, codons: &[&str], table: i64) -> String {
    format!(
        r#"{{"stable_id": "{stable_id}", "variation_effect_feature_cache": {{"translateable_seq": "{}", "three_prime_utr": "", "codon_table": {table}}}}}"#,
        codons.concat()
    )
}

/// A candidate JSON cache whose chromosome 21 holds the six made-up transcripts.
fn planted_cache(directory: &Path) -> PathBuf {
    let shard = directory.join("cache/transcripts/21/1-1000000.json");
    std::fs::create_dir_all(shard.parent().expect("a shard has a directory"))
        .expect("the cache directory is made");
    let transcripts = [
        cached(KEPT, &KEPT_CODONS, 1),
        cached(UNCLOSED, &KEPT_CODONS[..6], 1),
        cached(UNHELD_CODE, &KEPT_CODONS, 7),
        cached(MITOCHONDRIAL, &MITOCHONDRIAL_CODONS, 2),
        cached(STOP_INSIDE, &[&KEPT_CODONS[..], &["GTT"]].concat(), 1),
        cached(PARTIAL_CODON, &[&KEPT_CODONS[..], &["G"]].concat(), 1),
    ];
    std::fs::write(&shard, format!("[{}]", transcripts.join(", "))).expect("the shard is written");
    directory.join("cache")
}

/// An HGVSp value in a text format: the protein, `:p.` and the description, its `=` escaped.
fn hgvsp(description: &str, protein: &str) -> String {
    format!("{protein}:p.{description}").replace('=', "%3D")
}

/// What the class says of the differing HGVSp values `ev` and `av` of a row of `feature`
/// whose terms are `terms` and whose FLAGS are `flags`, asked as a comparison asks it
/// beside the candidate cache `cache`, or beside none.
fn class_reads(
    feature: &str,
    terms: &str,
    flags: &str,
    ev: &str,
    av: &str,
    cache: Option<&Path>,
) -> Result<bool, Error> {
    class_reads_at(Some(LOCATION), feature, terms, flags, ev, av, cache)
}

/// `class_reads` on a row at `location`, `None` for an entry that names none, as a VCF or
/// a JSON entry does.
fn class_reads_at(
    location: Option<&str>,
    feature: &str,
    terms: &str,
    flags: &str,
    ev: &str,
    av: &str,
    cache: Option<&Path>,
) -> Result<bool, Error> {
    let held: &[(&str, &str)] = if flags.is_empty() {
        &[]
    } else {
        &[("FLAGS", flags)]
    };
    let e = entry(
        &format!("{LOCATION}|-|{feature}|Transcript"),
        location.is_none().then_some(0),
        location,
        "-",
        feature,
        terms,
        held,
    );
    let side = Side {
        downstream_sequences: cache
            .map(DownstreamSequences::new)
            .transpose()
            .expect("the cache has transcripts"),
        ..Side::default()
    };
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(STOP_CODON_CLASS, &e, "HGVSp", ev, av, None, &context)
}

/// The class on the two descriptions of one protein of `feature`, beside the cache.
fn reads_on(feature: &str, reference: &str, candidate: &str, cache: &Path) -> bool {
    class_reads(
        feature,
        "5_prime_UTR_variant",
        "",
        &hgvsp(reference, PROTEIN),
        &hgvsp(candidate, PROTEIN),
        Some(cache),
    )
    .expect("the class reads the row")
}

#[test]
fn residues_and_a_stop_that_repeat_the_protein_to_its_end_are_read_beside_the_cache() {
    let directory = scratch("inserted-residues-read-beside-the-cache");
    let cache = planted_cache(&directory);
    assert!(reads_on(KEPT, RESIDUES_AND_STOP, KEPT_STOP, &cache));
    // From any residue: two residues before the stop, and all but the first.
    assert!(reads_on(KEPT, "Leu4_Ser5insSerValTer", KEPT_STOP, &cache));
    assert!(reads_on(
        KEPT,
        "Met1_Ala2insAlaLysLeuSerValTer",
        KEPT_STOP,
        &cache
    ));
    // The row's terms are not read, and a Feature that carries a version is found.
    for terms in ["inframe_insertion,stop_gained", "coding_sequence_variant"] {
        let accepted = class_reads(
            &format!("{KEPT}.2"),
            terms,
            "",
            &hgvsp(RESIDUES_AND_STOP, PROTEIN),
            &hgvsp(KEPT_STOP, PROTEIN),
            Some(&cache),
        );
        assert!(accepted.expect("the class reads the row"), "{terms}");
    }
    // An entry that names no Location, as a VCF or a JSON entry, finds its transcript too,
    // and the JSON format's unescaped `=` reads as the escaped one.
    let unplaced = class_reads_at(
        None,
        KEPT,
        "5_prime_UTR_variant",
        "",
        &format!("{PROTEIN}:p.{RESIDUES_AND_STOP}"),
        &format!("{PROTEIN}:p.{KEPT_STOP}"),
        Some(&cache),
    );
    assert!(unplaced.expect("the class reads the row"));
    // Both values inside prediction parentheses read as outside them.
    assert!(reads_on(
        KEPT,
        &format!("({RESIDUES_AND_STOP})"),
        &format!("({KEPT_STOP})"),
        &cache
    ));
}

#[test]
fn without_the_candidate_cache_the_pair_is_a_difference() {
    let accepted = class_reads(
        KEPT,
        "inframe_insertion,stop_retained_variant",
        "",
        &hgvsp(RESIDUES_AND_STOP, PROTEIN),
        &hgvsp(KEPT_STOP, PROTEIN),
        None,
    );
    assert!(!accepted.expect("the class reads the row"));
}

#[test]
fn a_pair_the_coding_sequence_does_not_bear_out_is_a_difference() {
    let directory = scratch("inserted-residues-not-borne-out");
    let cache = planted_cache(&directory);
    let differences = [
        // An inserted residue, a left flank and a right flank that are not the protein's.
        ("Lys3_Leu4insLeuThrValTer", KEPT_STOP),
        ("Arg3_Leu4insLeuSerValTer", KEPT_STOP),
        ("Lys3_Ser4insLeuSerValTer", KEPT_STOP),
        // Residues that end before the protein's last, at either stop number.
        ("Lys3_Leu4insLeuSerTer", KEPT_STOP),
        ("Lys3_Leu4insLeuSerTer", "Ter6="),
        // Residues placed one residue early or late.
        ("Ala2_Lys3insLeuSerValTer", KEPT_STOP),
        ("Leu4_Ser5insLeuSerValTer", KEPT_STOP),
        // The stop at another residue number.
        ("Lys3_Leu4insLeuSerValTer", "Ter8="),
        ("Lys3_Leu4insLeuSerValTer", "Ter6="),
        // Flanks that are not neighbours, and a stop among the inserted residues.
        ("Lys3_Ser5insLeuSerValTer", KEPT_STOP),
        ("Lys3_Leu4insLeuTerValTer", KEPT_STOP),
        // No stop after the inserted residues, and the candidate's another notation.
        ("Lys3_Leu4insLeuSerVal", KEPT_STOP),
        ("Lys3_Leu4insLeuSerValTer", "Ter7ValextTer3"),
        // One value alone inside prediction parentheses.
        ("(Lys3_Leu4insLeuSerValTer)", KEPT_STOP),
        ("Lys3_Leu4insLeuSerValTer", "(Ter7=)"),
        // The sides exchanged.
        (KEPT_STOP, RESIDUES_AND_STOP),
    ];
    for (reference, candidate) in differences {
        assert!(
            !reads_on(KEPT, reference, candidate, &cache),
            "{reference} against {candidate}"
        );
    }
    // Two proteins, and a row whose FLAGS name an incomplete coding sequence.
    let two_proteins = class_reads(
        KEPT,
        "5_prime_UTR_variant",
        "",
        &hgvsp(RESIDUES_AND_STOP, PROTEIN),
        &hgvsp(KEPT_STOP, "ENSP00000000002.1"),
        Some(&cache),
    );
    assert!(!two_proteins.expect("the class reads the row"));
    let incomplete = class_reads(
        KEPT,
        "5_prime_UTR_variant",
        "cds_end_NF",
        &hgvsp(RESIDUES_AND_STOP, PROTEIN),
        &hgvsp(KEPT_STOP, PROTEIN),
        Some(&cache),
    );
    assert!(!incomplete.expect("the class reads the row"));
}

#[test]
fn a_transcript_the_cache_does_not_vouch_for_is_a_difference() {
    let directory = scratch("inserted-residues-no-voucher");
    let cache = planted_cache(&directory);
    // No stop codon after the residues, a stop codon that is not the last codon, a base after
    // the stop codon, a genetic code no table is held for, and a transcript the cache lacks.
    for feature in [
        UNCLOSED,
        STOP_INSIDE,
        PARTIAL_CODON,
        UNHELD_CODE,
        "ENST00000000399",
    ] {
        assert!(
            !reads_on(feature, RESIDUES_AND_STOP, KEPT_STOP, &cache),
            "{feature}"
        );
    }
}

#[test]
fn the_reference_protein_is_read_under_the_transcripts_own_genetic_code() {
    let directory = scratch("inserted-residues-own-genetic-code");
    let cache = planted_cache(&directory);
    assert!(reads_on(
        MITOCHONDRIAL,
        "Lys3_Met4insMetTrpValTer",
        KEPT_STOP,
        &cache
    ));
    // The standard code's reading of the same codons names another residue and stops early.
    assert!(!reads_on(
        MITOCHONDRIAL,
        "Lys3_Ile4insIleTer",
        "Ter5=",
        &cache
    ));
}

#[test]
fn one_residue_before_the_stop_is_read_from_the_row_alone() {
    let directory = scratch("inserted-residue-row-alone");
    let cache = planted_cache(&directory);
    let one_residue = |terms: &str, cache: Option<&Path>| {
        class_reads(
            KEPT,
            terms,
            "",
            &hgvsp("Ser5_Val6insValTer", PROTEIN),
            &hgvsp(KEPT_STOP, PROTEIN),
            cache,
        )
        .expect("the class reads the row")
    };
    assert!(one_residue("inframe_insertion,stop_retained_variant", None));
    assert!(one_residue(
        "inframe_insertion,stop_retained_variant",
        Some(&cache)
    ));
    assert!(!one_residue("5_prime_UTR_variant", Some(&cache)));
}

#[test]
fn a_comparison_reads_the_pair_where_it_is_given_the_candidate_cache() {
    let directory = scratch("inserted-residues-command-line");
    let cache = text(&planted_cache(&directory));
    let row = |description: &str| {
        default_row(
            "made_up_insertion",
            LOCATION,
            "-",
            KEPT,
            "5_prime_UTR_variant",
            &format!("IMPACT=MODIFIER;HGVSp={}", hgvsp(description, PROTEIN)),
        )
    };
    let reference = directory.join("reference.txt");
    let candidate = directory.join("candidate.txt");
    std::fs::write(&reference, default_text(&[row(RESIDUES_AND_STOP)])).unwrap();
    std::fs::write(&candidate, default_text(&[row(KEPT_STOP)])).unwrap();
    let compared = |more: &[&str], name: &str| {
        let path = directory.join(name);
        let mut argv: Vec<OsString> = [
            "--format",
            "default",
            "--expected",
            &text(&reference),
            "--actual",
            &text(&candidate),
            "--masks",
            STOP_CODON_CLASS,
            "--report",
            &text(&path),
        ]
        .iter()
        .map(OsString::from)
        .collect();
        argv.extend(more.iter().map(OsString::from));
        (run(&directory, &argv), path)
    };
    let (ran, path) = compared(&[], "without.json");
    assert_eq!(ran.code, 1, "{}", ran.err());
    assert_eq!(compact(member(&report(&path), &["masks", "applied"])), "{}");
    let (ran, path) = compared(&["--json-cache", &cache], "beside.json");
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        compact(member(&report(&path), &["masks", "applied"])),
        format!(r#"{{"{STOP_CODON_CLASS}": 1}}"#)
    );
}
