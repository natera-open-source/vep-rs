// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A plugin value that one engine alone writes beside a pair of Consequence sets at a coding
//! end of a transcript, on made-up records. The candidate cache's reading of the transcript
//! decides whether a class of a coding end takes the pair, and each side's terms decide
//! whether its plugin runs: a value is taken where it stands exactly on the side whose
//! terms let the plugin run, beside a pair a class takes.
//!
//! A base is written as its place in the order of a genetic code's codons, `T`, `C`, `A`
//! and `G`, and a codon as the places of its three bases.

use std::collections::BTreeSet;
use vep_plugin_pairs::pytext::{json_loads, Json};
use vep_plugin_pairs::rules::coding_end::coding_end_class;
use vep_plugin_pairs::rules::coding_ends::terms::CodingEndClass;
use vep_plugin_pairs::rules::coding_ends::{StartEnd, StopEnd};
use vep_plugin_pairs::rules::downstream::DownstreamSequences;
use vep_plugin_pairs::rules::entry::{parse_default, Row, ABSENT};
use vep_plugin_pairs::rules::{BuiltIn, Consult, Rules};

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

const FOLLOWS_STOP: &str = "plugin_value_follows_stop_codon_term_pair";
const START_PAST_OWN_CODON: &str = "downstream_start_past_own_codon";

/// The first genomic base of a made-up transcript and the bases of its 5' UTR.
const FIRST: i64 = 1001;
const LEADER: i64 = 6;

const BASES: [char; 4] = ['T', 'C', 'A', 'G'];
type Codon = [usize; 3];
/// Codons of glutamate, methionine, alanine and glycine, and a stop codon of the standard
/// and the vertebrate mitochondrial code.
const GLU: Codon = [3, 2, 2];
const MET: Codon = [2, 0, 3];
const ALA: Codon = [3, 1, 0];
const GLY: Codon = [3, 3, 0];
const STOP: Codon = [0, 2, 2];

/// A first codon that is no start codon, closed by a stop codon: no start codon is shown.
const START_OPEN: &str = "ENST00000000301";
/// Whole codons from a start codon to a stop codon: both codons stand.
const CLOSED: &str = "ENST00000000302";
/// Whole codons whose last names a residue, before a 3' UTR: no stop codon stands.
const OPEN_WHOLE: &str = "ENST00000000303";
/// A coding sequence that ends one base into a stop codon on the transcript's last base,
/// under the vertebrate mitochondrial code: polyadenylation completes the stop codon.
const SHOWN_TAIL: &str = "ENST00000000304";
/// The same end under the standard code: no record shows the stop codon completed.
const TAIL_NOT_SHOWN: &str = "ENST00000000305";
/// One stable id on two records of one sequence, the first shaped as `START_OPEN` and the
/// second, 2,000 bases on, as `CLOSED`.
const TWO_RECORDS: &str = "ENST00000000306";
/// A coding sequence that ends two bases into a codon before a 3' UTR: no stop codon
/// stands, and the last codon lacks one base.
const OPEN_TAIL: &str = "ENST00000000307";

fn bases(codons: &[Codon]) -> String {
    codons.iter().flatten().map(|place| BASES[*place]).collect()
}

/// A cache transcript of one exon on the forward strand from the genomic base `first`:
/// `LEADER` bases of 5' UTR, the coding sequence `coding` and `trailer` bases of 3' UTR,
/// under the genetic code `table`, without a cached peptide.
fn record(stable_id: &str, first: i64, coding: &str, trailer: i64, table: i64) -> Json {
    let last = LEADER + coding.len() as i64 + trailer;
    let (start, end) = (first, first + last - 1);
    let (coding_start, coding_end) = (LEADER + 1, LEADER + coding.len() as i64);
    let utr = BASES[1].to_string().repeat(trailer as usize);
    let text = format!(
        r#"{{"stable_id": "{stable_id}", "strand": 1, "start": {start}, "end": {end}, "cdna_coding_start": {coding_start}, "cdna_coding_end": {coding_end}, "exons": [{{"start": {start}, "end": {end}}}], "variation_effect_feature_cache": {{"translateable_seq": "{coding}", "three_prime_utr": "{utr}", "codon_table": {table}, "mapper": {{"cdna_coding_start": {coding_start}, "cdna_coding_end": {coding_end}, "pairs": [{{"from_start": 1, "from_end": {last}, "to_start": {start}, "to_end": {end}, "ori": 1}}]}}}}}}"#
    );
    json_loads(&text).expect("a transcript object")
}

/// The made-up transcripts, all on chromosome 21.
fn sequences() -> DownstreamSequences {
    let tail = format!("{}{}", bases(&[MET, ALA, GLY]), BASES[STOP[0]]);
    let open_tail = format!("{tail}{}", BASES[STOP[1]]);
    DownstreamSequences::of_transcripts(vec![(
        "21".to_string(),
        vec![
            record(START_OPEN, FIRST, &bases(&[GLU, ALA, GLY, STOP]), 6, 1),
            record(CLOSED, FIRST, &bases(&[MET, ALA, GLY, STOP]), 6, 1),
            record(OPEN_WHOLE, FIRST, &bases(&[MET, ALA, GLY]), 6, 1),
            record(SHOWN_TAIL, FIRST, &tail, 0, 2),
            record(TAIL_NOT_SHOWN, FIRST, &tail, 0, 1),
            record(TWO_RECORDS, FIRST, &bases(&[GLU, ALA, GLY, STOP]), 6, 1),
            record(
                TWO_RECORDS,
                FIRST + 2000,
                &bases(&[MET, ALA, GLY, STOP]),
                6,
                1,
            ),
            record(OPEN_TAIL, FIRST, &open_tail, 6, 1),
        ],
    )])
}

/// A default-format row of a transcript under an identifier: its Location, allele and
/// terms, its five coding cells (`cDNA_position`, `CDS_position`, `Protein_position`,
/// `Amino_acids`, `Codons`) and its Extra text, the allele standing in for the reference
/// bases `reference` where the row is given them.
fn row(
    location: &str,
    allele: &str,
    feature: &str,
    terms: &str,
    cells: &[String; 5],
    extra: &str,
    reference: Option<&str>,
) -> Row {
    let extra = if extra.is_empty() { "-" } else { extra };
    let text = format!(
        "{HEADER}rs1\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{terms}\t{}\t-\t{extra}\n",
        cells.join("\t")
    );
    let row = Row::clone(&parse_default(&text).expect("a row of fourteen columns")[0]);
    *row.reference.borrow_mut() = reference.map(str::to_string);
    row
}

/// An `Amino_acids` cell of one residue over another.
fn residues(reference: char, alternate: char) -> String {
    format!("{reference}/{alternate}")
}

/// A `Codons` cell: the codon with its base at `place` in upper case and its others in
/// lower case, over the same codon with the base `alternate` there, or without a base
/// there for `None`.
fn codons(codon: Codon, place: usize, alternate: Option<usize>) -> String {
    let written = |changed: Option<usize>| -> String {
        (0..3)
            .filter_map(|at| match (at == place, changed) {
                (false, _) => Some(BASES[codon[at]].to_ascii_lowercase()),
                (true, Some(base)) => Some(BASES[base]),
                (true, None) => None,
            })
            .collect()
    };
    format!("{}/{}", written(Some(codon[place])), written(alternate))
}

fn cells(cdna: &str, cds: &str, protein: &str, amino_acids: String, codons: String) -> [String; 5] {
    [
        cdna.to_string(),
        cds.to_string(),
        protein.to_string(),
        amino_acids,
        codons,
    ]
}

/// A substitution of the third base of the first codon on `feature`, a transcript shaped
/// as `START_OPEN` from the genomic base `first`: Ensembl VEP reads the first three coding
/// positions as the start codon and writes `start_lost`; where no start codon is shown
/// vep-rs owes the term of the two residues, `missense_variant`. Each row carries its side's
/// Extra text.
fn start_pair(
    feature: &str,
    first: i64,
    reference_extra: &str,
    candidate_extra: &str,
) -> (Row, Row) {
    let location = format!("21:{}", first + LEADER + 2);
    let cells = cells("9", "3", "1", residues('E', 'D'), codons(GLU, 2, Some(0)));
    let (allele, reference) = (BASES[0].to_string(), BASES[GLU[2]].to_string());
    let side = |terms: &str, extra: &str| {
        row(
            &location,
            &allele,
            feature,
            terms,
            &cells,
            extra,
            Some(&reference),
        )
    };
    (
        side("start_lost", reference_extra),
        side("missense_variant", candidate_extra),
    )
}

/// A substitution of the one tail base on `feature`, a transcript shaped as `SHOWN_TAIL`:
/// Ensembl VEP reads an incomplete terminal codon; where polyadenylation completes the stop
/// codon vep-rs owes `stop_lost` for a base that leaves the codon of a residue.
fn completed_stop_pair(feature: &str, reference_extra: &str, candidate_extra: &str) -> (Row, Row) {
    let (allele, reference) = (BASES[1].to_string(), BASES[STOP[0]].to_string());
    let side = |terms: &str, amino_acids: &str, extra: &str| {
        row(
            "21:1016",
            &allele,
            feature,
            terms,
            &cells("16", "10", "4", amino_acids.to_string(), "-".to_string()),
            extra,
            Some(&reference),
        )
    };
    (
        side(
            "incomplete_terminal_codon_variant,coding_sequence_variant",
            "X",
            reference_extra,
        ),
        side("stop_lost", "-", candidate_extra),
    )
}

/// A deletion of the second base of the last whole codon on `feature`, a transcript shaped
/// as `OPEN_WHOLE`: Ensembl VEP reads the last three coding positions as the stop codon and
/// writes a kept stop codon beside an in-frame term; where no stop codon stands vep-rs owes
/// `frameshift_variant`.
fn kept_stop_pair(feature: &str, reference_extra: &str, candidate_extra: &str) -> (Row, Row) {
    (
        kept_stop_row(
            feature,
            "inframe_deletion,stop_retained_variant",
            reference_extra,
        ),
        kept_stop_row(feature, "frameshift_variant", candidate_extra),
    )
}

/// One side's row of `kept_stop_pair` with the terms `terms`.
fn kept_stop_row(feature: &str, terms: &str, extra: &str) -> Row {
    let cells = cells("14", "8", "3", residues('G', 'X'), codons(GLY, 1, None));
    row("21:1014", "-", feature, terms, &cells, extra, None)
}

/// One question to a class about the field `name` of the reference's row `e` and the
/// candidate's `a`, each side's value the field's own or absent.
fn asked<T>(
    klass: &str,
    e: &Row,
    a: &Row,
    name: &str,
    sequences: Option<&DownstreamSequences>,
    ask: impl FnOnce(&Consult) -> T,
) -> T {
    ask(&Consult {
        klass,
        e,
        name,
        expected: e.fields.get_or(name, ABSENT),
        actual: a.fields.get_or(name, ABSENT),
        a: Some(a),
        record_rows: None,
        candidate_rows: None,
        declared: &BTreeSet::new(),
        assembly: None,
        source: None,
        gwas_pairs: None,
        revel_table: None,
        utrannotator_sweep: None,
        downstream_sequences: sequences,
        splice_windows: None,
        spliceai_table: None,
        exon_spans: None,
    })
}

/// Whether the tool's rules take the field under `klass`: the class is asked only of a
/// field its row of `data/field_classes.tsv` names.
fn accepts(
    klass: &str,
    e: &Row,
    a: &Row,
    name: &str,
    sequences: Option<&DownstreamSequences>,
) -> bool {
    asked(klass, e, a, name, sequences, |consult| {
        BuiltIn.accepts(consult).expect("no refusal")
    })
}

/// Whether `plugin_value_follows_stop_codon_term_pair` takes the field of the two rows.
fn follows(e: &Row, a: &Row, name: &str, sequences: Option<&DownstreamSequences>) -> bool {
    accepts(FOLLOWS_STOP, e, a, name, sequences)
}

/// The class of a coding end that takes the pair beside the cache.
fn class(e: &Row, a: &Row, sequences: &DownstreamSequences) -> Option<CodingEndClass> {
    coding_end_class(e, a, sequences.coding_ends(e).expect("the cache reads"))
}

#[test]
fn a_row_is_answered_with_the_coding_ends_of_the_record_that_holds_its_location() {
    let sequences = sequences();
    let at = |location: &str, feature: &str| {
        let blank = cells("-", "-", "-", "-".to_string(), "-".to_string());
        let row = row(location, "-", feature, "intron_variant", &blank, "", None);
        sequences
            .coding_ends(&row)
            .expect("the cache reads")
            .map(|ends| (ends.stop, ends.start))
    };
    let not_shown = StartEnd::NotShown {
        flagged: false,
        padded: false,
    };
    let open = StopEnd::NoneShown { flagged: false };
    for (feature, ends) in [
        (START_OPEN, (StopEnd::Whole, not_shown)),
        (CLOSED, (StopEnd::Whole, StartEnd::Shown)),
        (OPEN_WHOLE, (open, StartEnd::Shown)),
        (
            SHOWN_TAIL,
            (StopEnd::BegunByTheTail { bases: 1 }, StartEnd::Shown),
        ),
        (TAIL_NOT_SHOWN, (open, StartEnd::Shown)),
    ] {
        assert_eq!(at("21:1009", feature), Some(ends), "{feature}");
        // A Feature with a version the cache's stable id lacks, and a Location whose
        // sequence the cache has no directory for, which reads every directory.
        assert_eq!(
            at("21:1009", &format!("{feature}.4")),
            Some(ends),
            "{feature}"
        );
        assert_eq!(at("7:1009", feature), Some(ends), "{feature}");
        // A Location off the record's span, and a transcript the cache does not hold.
        assert_eq!(at("21:999", feature), None, "{feature}");
        assert_eq!(at("21:1000-1001", feature), Some(ends), "{feature}");
    }
    assert_eq!(at("21:1009", "ENST00000000399"), None);
    // One stable id on two records: the record whose span holds the Location answers,
    // and none where no record or more than one holds it.
    assert_eq!(
        at("21:1009", TWO_RECORDS),
        Some((StopEnd::Whole, not_shown))
    );
    assert_eq!(
        at("21:3009", TWO_RECORDS),
        Some((StopEnd::Whole, StartEnd::Shown))
    );
    assert_eq!(at("21:2000", TWO_RECORDS), None);
    assert_eq!(at("21:1009-3009", TWO_RECORDS), None);
}

#[test]
fn a_value_beside_a_start_codon_pair_is_taken_on_the_side_whose_terms_let_its_plugin_run() {
    let sequences = sequences();
    let scores = "REVEL=0.5;am_class=likely_benign;am_pathogenicity=0.1";
    let (e, a) = start_pair(START_OPEN, FIRST, "", scores);
    assert_eq!(
        class(&e, &a, &sequences),
        Some(CodingEndClass::StartCodonOverlap)
    );
    for name in ["REVEL", "am_class", "am_pathogenicity"] {
        assert!(follows(&e, &a, name, Some(&sequences)), "{name}");
        // A comparison given no cache takes nothing of a class of a coding end.
        assert!(!follows(&e, &a, name, None), "{name}: no cache");
    }
    // The same rows on a transcript whose start codon stands: no class takes the pair.
    for (feature, first) in [(CLOSED, FIRST), (TWO_RECORDS, FIRST + 2000)] {
        let (e, a) = start_pair(feature, first, "", scores);
        assert_eq!(class(&e, &a, &sequences), None, "{feature}");
        assert!(!follows(&e, &a, "REVEL", Some(&sequences)), "{feature}");
    }
    // The record of the same stable id that shows no start codon.
    let (e, a) = start_pair(TWO_RECORDS, FIRST, "", scores);
    assert!(follows(&e, &a, "REVEL", Some(&sequences)));
    // The value on Ensembl VEP's side alone: its terms hold no `missense_variant`.
    let (e, a) = start_pair(START_OPEN, FIRST, scores, "");
    assert!(!follows(&e, &a, "REVEL", Some(&sequences)));
    assert!(!follows(&e, &a, "am_class", Some(&sequences)));
    // ClinPred and dbNSFP run on `start_lost` as on `missense_variant`: a value one side
    // alone writes there does not follow the terms.
    let (e, a) = start_pair(START_OPEN, FIRST, "", "ClinPred=0.9;dbNSFP_SIFT_score=0.2");
    assert!(!follows(&e, &a, "ClinPred", Some(&sequences)));
    assert!(!follows(&e, &a, "dbNSFP_SIFT_score", Some(&sequences)));
    // vep-rs's row keeps Ensembl VEP's term beside the residues' term: no class takes it.
    let (e, a) = start_pair(START_OPEN, FIRST, "", "");
    let kept = row(
        "21:1009",
        &a.allele,
        START_OPEN,
        "start_lost,missense_variant",
        &cells("9", "3", "1", residues('E', 'D'), codons(GLU, 2, Some(0))),
        scores,
        a.reference.borrow().as_deref(),
    );
    assert_eq!(class(&e, &kept, &sequences), None);
    assert!(!follows(&e, &kept, "REVEL", Some(&sequences)));
}

#[test]
fn a_value_beside_a_completed_stop_codon_pair_is_taken_on_the_side_that_writes_the_lost_stop() {
    let sequences = sequences();
    let scores = "ClinPred=0.9;dbNSFP_SIFT_score=0.2;SIFT_score=0.2;REVEL=0.5";
    let (e, a) = completed_stop_pair(SHOWN_TAIL, "", scores);
    assert_eq!(
        class(&e, &a, &sequences),
        Some(CodingEndClass::StopCodonCompleted)
    );
    for name in ["ClinPred", "dbNSFP_SIFT_score"] {
        assert!(follows(&e, &a, name, Some(&sequences)), "{name}");
        assert!(!follows(&e, &a, name, None), "{name}: no cache");
    }
    // REVEL runs on neither side's terms, and a dbNSFP column under its own name is not
    // told from another plugin's field.
    assert!(!follows(&e, &a, "REVEL", Some(&sequences)));
    assert!(!follows(&e, &a, "SIFT_score", Some(&sequences)));
    // The same rows where no record shows the stop codon completed.
    let (e, a) = completed_stop_pair(TAIL_NOT_SHOWN, "", scores);
    assert_eq!(class(&e, &a, &sequences), None);
    assert!(!follows(&e, &a, "ClinPred", Some(&sequences)));
    // The value on Ensembl VEP's side alone: its terms let neither plugin run.
    let (e, a) = completed_stop_pair(SHOWN_TAIL, scores, "");
    assert!(!follows(&e, &a, "ClinPred", Some(&sequences)));
    assert!(!follows(&e, &a, "dbNSFP_SIFT_score", Some(&sequences)));
}

#[test]
fn a_value_beside_a_kept_stop_codon_pair_is_taken_where_no_stop_codon_stands() {
    let sequences = sequences();
    let written = "DownstreamProtein=value;ProteinLengthChange=4;NMD=NMD_escaping_variant;LoF=HC";
    let (e, a) = kept_stop_pair(OPEN_WHOLE, "", written);
    assert_eq!(
        class(&e, &a, &sequences),
        Some(CodingEndClass::StopCodonOverlap)
    );
    for name in ["DownstreamProtein", "ProteinLengthChange", "NMD", "LoF"] {
        assert!(
            accepts(FOLLOWS_STOP, &e, &a, name, Some(&sequences)),
            "{name}"
        );
        assert!(
            !accepts(FOLLOWS_STOP, &a, &e, name, Some(&sequences)),
            "{name}: the sides exchanged"
        );
    }
    // The same rows on a transcript whose stop codon stands: the pair is one of a class of
    // a stop codon that stands, beside the cache as without it.
    let (e, a) = kept_stop_pair(CLOSED, "", written);
    assert_eq!(class(&e, &a, &sequences), None);
    for name in ["DownstreamProtein", "NMD", "LoF"] {
        assert!(
            accepts(FOLLOWS_STOP, &e, &a, name, Some(&sequences)),
            "{name}"
        );
        assert!(
            accepts(FOLLOWS_STOP, &e, &a, name, None),
            "{name}: no cache"
        );
    }
    // Where no stop codon stands such a class takes nothing beside the cache. A row of
    // vep-rs with a gained stop the derivation does not owe has the terms of its pair and
    // is no pair of a class of a coding end: nothing takes it beside the cache, and the
    // two rows alone take it where the comparison is given none.
    let (e, _) = kept_stop_pair(OPEN_WHOLE, "", "");
    let gained = kept_stop_row(OPEN_WHOLE, "stop_gained,frameshift_variant", written);
    assert_eq!(class(&e, &gained, &sequences), None);
    for name in ["DownstreamProtein", "NMD", "LoF"] {
        assert!(
            !accepts(FOLLOWS_STOP, &e, &gained, name, Some(&sequences)),
            "{name}"
        );
        assert!(
            accepts(FOLLOWS_STOP, &e, &gained, name, None),
            "{name}: no cache"
        );
    }
    // A transcript the cache does not hold: beside the cache nothing is taken.
    let (e, a) = kept_stop_pair("ENST00000000399", "", written);
    assert!(!accepts(
        FOLLOWS_STOP,
        &e,
        &a,
        "DownstreamProtein",
        Some(&sequences)
    ));
}

/// A deletion of the last two bases of the last whole codon and of the two bases behind it
/// on `OPEN_TAIL`: Ensembl VEP writes a kept stop codon beside an in-frame term, and vep-rs
/// owes `frameshift_variant`. The reference's alternate codon is the base the deletion
/// leaves and one base of the 3' UTR; the candidate's is `alternate`, and its row carries
/// the Extra text `extra`.
fn deletion_into_the_tail_pair(alternate: &str, extra: &str) -> (Row, Row) {
    let left = BASES[GLY[0]].to_ascii_lowercase();
    let removed: String = [GLY[1], GLY[2], STOP[0], STOP[1]]
        .iter()
        .map(|place| BASES[*place])
        .collect();
    let side = |terms: &str, alternate: &str, extra: &str| {
        let codons = format!("{left}{removed}/{alternate}");
        let cells = cells("14-17", "8-11", "3-4", "GX/X".to_string(), codons);
        row("21:1014-1017", "-", OPEN_TAIL, terms, &cells, extra, None)
    };
    let read_on = format!("{left}{}", BASES[1].to_ascii_lowercase());
    (
        side("inframe_deletion,stop_retained_variant", &read_on, ""),
        side("frameshift_variant", alternate, extra),
    )
}

#[test]
fn a_value_beside_a_kept_stop_codon_pair_is_taken_where_the_alternate_codons_differ_past_the_coding_end(
) {
    let sequences = sequences();
    let written = "DownstreamProtein=value;ProteinLengthChange=-2;NMD=NMD_escaping_variant";
    let left = BASES[GLY[0]].to_ascii_lowercase().to_string();
    // The candidate's alternate codon ends where the coding sequence ends, the pair of
    // `alternate_codon_past_the_coding_end`.
    let (e, a) = deletion_into_the_tail_pair(&left, written);
    assert_ne!(e.fields.get("Codons"), a.fields.get("Codons"));
    assert_eq!(
        class(&e, &a, &sequences),
        Some(CodingEndClass::StopCodonOverlap)
    );
    for name in ["DownstreamProtein", "ProteinLengthChange", "NMD"] {
        assert!(
            accepts(FOLLOWS_STOP, &e, &a, name, Some(&sequences)),
            "{name}"
        );
        assert!(
            !accepts(FOLLOWS_STOP, &a, &e, name, Some(&sequences)),
            "{name}: the sides exchanged"
        );
    }
    // Another alternate codon of the candidate's: a base of the allele, and three bases.
    for alternate in [left.to_ascii_uppercase(), format!("{left}{left}{left}")] {
        let (e, a) = deletion_into_the_tail_pair(&alternate, written);
        assert_eq!(class(&e, &a, &sequences), None, "{alternate}");
        for name in ["DownstreamProtein", "ProteinLengthChange", "NMD"] {
            assert!(
                !accepts(FOLLOWS_STOP, &e, &a, name, Some(&sequences)),
                "{name}: {alternate}"
            );
        }
    }
}

#[test]
fn the_downstream_start_class_takes_a_value_beside_a_start_codon_pair() {
    let sequences = sequences();
    // A deletion of the third base of the first codon: a frameshift on both sides, with
    // `start_lost` on Ensembl VEP's where no start codon is shown.
    let pair = |feature: &str, candidate_terms: &str| {
        let cells = cells("9", "3", "1", residues('E', 'X'), codons(GLU, 2, None));
        let blank = row(
            "21:1009",
            "-",
            feature,
            "frameshift_variant",
            &cells,
            "",
            None,
        );
        let read = sequences
            .transcript(&blank)
            .expect("the cache reads")
            .expect("the transcript");
        let value = |length_change: bool, ensembl: bool| {
            read.downstream_value(3, 3, "", length_change, ensembl)
                .expect("a value")
        };
        let (theirs, ours, change) = (value(false, true), value(false, false), value(true, false));
        assert_ne!(theirs, ours);
        let side = |terms: &str, protein: &str| {
            let extra = format!("DownstreamProtein={protein};ProteinLengthChange={change}");
            row("21:1009", "-", feature, terms, &cells, &extra, None)
        };
        (
            side("frameshift_variant,start_lost", &theirs),
            side(candidate_terms, &ours),
        )
    };
    let (e, a) = pair(START_OPEN, "frameshift_variant");
    assert_eq!(
        class(&e, &a, &sequences),
        Some(CodingEndClass::StartCodonOverlap)
    );
    assert!(accepts(
        START_PAST_OWN_CODON,
        &e,
        &a,
        "DownstreamProtein",
        Some(&sequences)
    ));
    // vep-rs's row keeps a start term, and the same rows where the start codon stands.
    let (e, a) = pair(START_OPEN, "frameshift_variant,start_retained_variant");
    assert!(!accepts(
        START_PAST_OWN_CODON,
        &e,
        &a,
        "DownstreamProtein",
        Some(&sequences)
    ));
    let (e, a) = pair(CLOSED, "frameshift_variant");
    assert!(!accepts(
        START_PAST_OWN_CODON,
        &e,
        &a,
        "DownstreamProtein",
        Some(&sequences)
    ));
}
