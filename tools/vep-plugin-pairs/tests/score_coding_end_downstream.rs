// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Downstream value on a row whose two Consequence sets differ as a class of a coding end
//! takes them, scored by the `score` program end to end on a made-up record: the class of a
//! coding end reads the count of reference bases the row's allele stands in for, which a row
//! under an identifier states only through the lines of the run's input.
//!
//! A base is written as its place in the order of a genetic code's codons, `T`, `C`, `A`
//! and `G`.

mod support;

use serde_json::{json, Value};
use std::path::Path;
use support::{
    feature, feature_residuals, message, try_score, PlantedRegistry, Run, Scratch, HEADER,
};
use vep_plugin_pairs::rules::downstream::DownstreamSequences;
use vep_plugin_pairs::rules::entry::parse_default;
use vep_plugin_pairs::rules::{BuiltIn, DOWNSTREAM_START_PAST_OWN_CODON};

const BASES: [char; 4] = ['T', 'C', 'A', 'G'];
/// Codons of methionine, alanine and glycine.
const MET: [usize; 3] = [2, 0, 3];
const ALA: [usize; 3] = [3, 1, 0];
const GLY: [usize; 3] = [3, 3, 0];

fn bases(places: &[usize]) -> String {
    places.iter().map(|place| BASES[*place]).collect()
}

/// A transcript of one exon on the forward strand of chromosome 21 from the genomic base
/// 1001: six bases of 5' UTR, three whole codons whose last names a residue and six bases
/// of 3' UTR, so that no stop codon stands at its coding end.
fn cache(scratch: &Scratch) -> String {
    let transcript = json!({
        "stable_id": "ENST1",
        "strand": 1,
        "start": 1001,
        "end": 1021,
        "cdna_coding_start": 7,
        "cdna_coding_end": 15,
        "exons": [{"start": 1001, "end": 1021}],
        "variation_effect_feature_cache": {
            "translateable_seq": bases(&[MET, ALA, GLY].concat()),
            "three_prime_utr": bases(&[1; 6]),
            "codon_table": 1,
            "mapper": {
                "cdna_coding_start": 7,
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

/// The third base of the second codon and the first base of the last codon replaced by one
/// base, under the identifier `rs7`: a frameshift on both sides, beside which Ensembl VEP
/// reads the last three coding positions as the stop codon and writes `stop_lost`. The
/// span starts on a codon's third base, where the two engines' DownstreamProtein differ as
/// `downstream_start_past_own_codon` takes them. With the lines of the input the row's
/// allele stands in for two reference bases, the class of the stop codon overlap takes the
/// pair of Consequence sets and the value is set aside; without them the count of the
/// reference bases is not known, no class of a coding end takes the pair and both values
/// stay residuals.
#[test]
fn a_downstream_value_beside_a_pair_at_a_coding_end_is_set_aside_where_the_input_states_the_record()
{
    let klass = DOWNSTREAM_START_PAST_OWN_CODON;
    let scratch = Scratch::new("coding-end-downstream");
    let registry =
        PlantedRegistry::new(&[("Downstream", &["DownstreamProtein", "ProteinLengthChange"])]);
    let cache = cache(&scratch);
    let (allele, replaced) = (bases(&[2]), bases(&[ALA[2], GLY[0]]));
    let lower = |places: &[usize]| bases(places).to_lowercase();
    let codons = format!(
        "{}{}{}/{}{}{}",
        lower(&ALA[..2]),
        replaced,
        lower(&GLY[1..]),
        lower(&ALA[..2]),
        allele,
        lower(&GLY[1..])
    );
    let amino_acids = format!("{}{}/{}", 'A', 'G', 'X');
    let row = |terms: &str, fields: &str| {
        format!(
            "rs7\t21:1012-1013\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t12-13\t6-7\t2-3\t{amino_acids}\t{codons}\t-\tIMPACT=HIGH;STRAND=1{fields}\n"
        )
    };
    // Each engine's value by its own rule on the row's CDS pair, read from the cache.
    let sequences = DownstreamSequences::new(Path::new(&cache)).expect("the cache");
    let unwritten = parse_default(&format!("{HEADER}{}", row("frameshift_variant", "")))
        .expect("a row of fourteen columns");
    let read = sequences
        .transcript(&unwritten[0])
        .expect("the cache reads")
        .expect("the transcript");
    let value = |length_change: bool, ensembl: bool| {
        read.downstream_value(6, 7, &allele, length_change, ensembl)
            .expect("a value")
    };
    let (theirs, ours, change) = (value(false, true), value(false, false), value(true, false));
    assert_ne!(theirs, ours);
    let fields =
        |protein: &str| format!(";DownstreamProtein={protein};ProteinLengthChange={change}");
    let reference = scratch.output(
        "reference.txt",
        &[row("frameshift_variant,stop_lost", &fields(&theirs))],
    );
    let engine = scratch.output("engine.txt", &[row("frameshift_variant", &fields(&ours))]);
    let run_input = scratch.vcf("input.vcf", &[["21", "1012", "rs7", &replaced, &allele]]);
    let table = format!("{klass}\tDownstream\tDownstreamProtein\tmet\t-\n");
    for (name, given, set_aside, left) in [
        ("given", Some(&run_input), Some(1), 0),
        ("without", None, Some(0), 2),
    ] {
        let residuals = scratch.path(&format!("residuals-{name}"));
        let mut run = Run::new()
            .table(&table)
            .plugins(&["Downstream"])
            .cache(&cache)
            .residuals(&residuals);
        if let Some(input) = given {
            run = run.run_input(input);
        }
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
