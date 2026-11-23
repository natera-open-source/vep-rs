// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Downstream value on a row whose two Consequence sets differ as
//! `start_codon_kept_without_five_prime_utr` takes them, scored by the `score` program end
//! to end on a made-up record: under `--hgvs` Ensembl VEP's plugin reads the CDS pair its
//! protein notation left at the 3' place, vep-rs the record's own pair, and the value is
//! one of `plugin_hgvs_shifted_cds_coordinates` beside either pair of sets.
//!
//! A base is written as its place in the order of a genetic code's codons, `T`, `C`, `A`
//! and `G`.

mod support;

use serde_json::{json, Value};
use std::path::Path;
use support::{
    feature, feature_residuals, message, try_score, PlantedRegistry, Run, Scratch, HEADER,
};
use vep_plugin_pairs::rules::downstream::classes::downstream_class;
use vep_plugin_pairs::rules::downstream::DownstreamSequences;
use vep_plugin_pairs::rules::entry::parse_default;
use vep_plugin_pairs::rules::{BuiltIn, PLUGIN_HGVS_SHIFTED_CDS_COORDINATES};

const BASES: [char; 4] = ['T', 'C', 'A', 'G'];
/// Codons of methionine, cysteine, alanine and asparagine, and a stop codon.
const MET: [usize; 3] = [2, 0, 3];
const CYS: [usize; 3] = [0, 3, 1];
const ALA: [usize; 3] = [3, 1, 0];
const ASN: [usize; 3] = [2, 2, 0];
const STOP: [usize; 3] = [0, 2, 2];

const LOST: &str = "frameshift_variant,start_lost";
const KEPT: &str = "frameshift_variant,start_retained_variant";

fn bases(places: &[usize]) -> String {
    places.iter().map(|place| BASES[*place]).collect()
}

/// A transcript of one exon on the forward strand of chromosome 21 from the genomic base
/// 1001 with no 5' UTR: four codons, a stop codon and six bases of 3' UTR. The second and
/// third coding bases stand again as the fourth and fifth, so a deletion of the first two
/// leaves the start codon and moves two bases 3'.
fn cache(scratch: &Scratch) -> String {
    let transcript = json!({
        "stable_id": "ENST1",
        "strand": 1,
        "start": 1001,
        "end": 1021,
        "cdna_coding_start": 1,
        "cdna_coding_end": 15,
        "exons": [{"start": 1001, "end": 1021}],
        "variation_effect_feature_cache": {
            "translateable_seq": bases(&[MET, CYS, ALA, ASN, STOP].concat()),
            "three_prime_utr": bases(&[1; 6]),
            "codon_table": 1,
            "mapper": {
                "cdna_coding_start": 1,
                "cdna_coding_end": 15,
                "pairs": [
                    {"from_start": 1, "from_end": 21, "to_start": 1001, "to_end": 1021, "ori": 1}
                ]
            },
        },
    });
    scratch.file(
        "cache/transcripts/21/1-1000000.json",
        json!([transcript]).to_string(),
    );
    scratch.path("cache")
}

/// The second and third coding bases deleted, under the identifier `rs7`, with an
/// HGVS_OFFSET of 2: a frameshift on both sides that leaves the start codon, which
/// Ensembl VEP writes `start_lost` on a transcript with no 5' UTR and vep-rs
/// `start_retained_variant`. Ensembl VEP's DownstreamProtein is Downstream.pm on the pair
/// moved by the offset, the second codon on; vep-rs's is its rule on the own pair, the
/// start codon on. The value is set aside under `plugin_hgvs_shifted_cds_coordinates`
/// beside that pair of sets as beside one set. The same two sets on a row whose cDNA span
/// is not its CDS span are no pair of the class of the kept start, and both values stay
/// residuals.
#[test]
fn a_downstream_value_at_the_shifted_pair_is_set_aside_beside_the_terms_of_a_kept_start() {
    let klass = PLUGIN_HGVS_SHIFTED_CDS_COORDINATES;
    let scratch = Scratch::new("kept-start-downstream");
    let registry =
        PlantedRegistry::new(&[("Downstream", &["DownstreamProtein", "ProteinLengthChange"])]);
    let cache = cache(&scratch);
    let row = |terms: &str, cdna: &str, fields: &str| {
        format!(
            "rs7\t21:1002-1003\t-\tENSG1\tENST1\tTranscript\t{terms}\t{cdna}\t2-3\t1\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVS_OFFSET=2{fields}\n"
        )
    };
    // Each engine's value by its own rule, read from the cache: the reference's on the
    // pair two bases 3' of the own one, the candidate's on the own pair.
    let sequences = DownstreamSequences::new(Path::new(&cache)).expect("the cache");
    let unwritten = parse_default(&format!("{HEADER}{}", row(LOST, "2-3", "")))
        .expect("a row of fourteen columns");
    let read = sequences
        .transcript(&unwritten[0])
        .expect("the cache reads")
        .expect("the transcript");
    let value = |pair: (i128, i128), length_change: bool, ensembl: bool| {
        read.downstream_value(pair.0, pair.1, "", length_change, ensembl)
            .expect("a value")
    };
    let (own, moved) = ((2, 3), (4, 5));
    let (theirs, ours) = (value(moved, false, true), value(own, false, false));
    // The start codon the deletion leaves, then the one residue before the new stop.
    assert_eq!((theirs.as_str(), ours.as_str()), ("R", "MR"));
    // Without the flag Ensembl VEP's plugin writes the candidate's value on the own pair.
    assert_eq!(value(own, false, true), ours);
    let change = value(own, true, false);
    assert_eq!(value(moved, true, true), change);
    let fields =
        |protein: &str| format!(";DownstreamProtein={protein};ProteinLengthChange={change}");
    let run_input = scratch.vcf(
        "input.vcf",
        &[["21", "1001", "rs7", &bases(&MET), &bases(&MET[..1])]],
    );
    let table = format!("{klass}\tDownstream\tDownstreamProtein\tmet\t-\n");
    for (name, reference_terms, cdna, taken, set_aside, left) in [
        ("one set", KEPT, "2-3", true, Some(1), 0),
        ("the kept start", LOST, "2-3", true, Some(1), 0),
        ("behind a 5' UTR", LOST, "8-9", false, Some(0), 2),
    ] {
        let reference_row = row(reference_terms, cdna, &fields(&theirs));
        let engine_row = row(KEPT, cdna, &fields(&ours));
        let rows = parse_default(&format!("{HEADER}{reference_row}{engine_row}"))
            .expect("rows of fourteen columns");
        let decided = downstream_class(
            &rows[0],
            &rows[1],
            "DownstreamProtein",
            &theirs,
            &ours,
            &sequences,
            None,
        )
        .expect("no refusal");
        assert_eq!(decided, taken.then_some((klass, None)), "{name}");
        let reference = scratch.output("reference.txt", &[reference_row]);
        let engine = scratch.output("engine.txt", &[engine_row]);
        let residuals = scratch.path(&format!("residuals-{name}"));
        let run = Run::new()
            .table(&table)
            .plugins(&["Downstream"])
            .cache(&cache)
            .residuals(&residuals)
            .run_input(&run_input);
        let report = try_score(
            &scratch,
            &reference,
            &[("engine", &engine)],
            &run,
            &registry,
            &BuiltIn,
        )
        .unwrap_or_else(|stop| panic!("the run was stopped: {}", message(&stop)));
        let masked = &feature(&report, "Downstream")["masked"];
        assert_eq!(
            masked.get(klass).and_then(Value::as_i64),
            set_aside,
            "{name}"
        );
        assert_eq!(
            feature_residuals(&residuals, "Downstream").len(),
            left,
            "{name}"
        );
    }
}
