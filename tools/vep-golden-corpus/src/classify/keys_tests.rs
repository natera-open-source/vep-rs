// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the key classes: the class file against the rules it names, the set rules
//! by taxonomy class, the mate-side read on a model written by hand, the reference bases a
//! run's records give a row, the rows of a tuple asked pair by pair, the residues and
//! codons a row is put to the rules with, and the rules that are asked.

use super::*;
use crate::constants::UNEXPLAINED;
use crate::mate_context::TranscriptModel;
use crate::rules::row::pair_row_with_codons;

fn run_record(vcf: &str) -> RunRecord {
    RunRecord {
        vcf: vcf.to_string(),
        manifest_index: None,
        reference_warning: None,
    }
}

fn no_records() -> (HashSet<(String, String)>, HashSet<String>) {
    (HashSet::new(), HashSet::new())
}

#[test]
fn the_class_file_names_every_pair_rule_in_the_order_it_is_asked() {
    let listed = |registry: &str| -> Vec<(&str, &str)> {
        pair_classes()
            .iter()
            .filter(|class| class.registry == registry)
            .map(|class| (class.bucket, class.taxonomy_class))
            .collect()
    };
    let set_rules: Vec<&str> = EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| !rule.row_aware())
        .map(|rule| rule.bucket)
        .collect();
    let family = "splice_family_swap";
    assert_eq!(
        listed(SET_RULES),
        [
            ("intronic_overcall_swap", family),
            ("ppt_overcall_swap", family),
            ("frameshift_splice_cascade_swap", family),
            ("splice_lastwrite_swap", family),
            ("start_cooccurrence_swap", "start_cooccurrence_swap"),
            ("covered_splice_region_swap", family),
        ]
    );
    let in_file: Vec<&str> = listed(SET_RULES)
        .iter()
        .map(|(bucket, _)| *bucket)
        .collect();
    assert_eq!(in_file, set_rules);
    let row_rules: Vec<(&str, &str)> = row_pair_rules()
        .iter()
        .map(|rule| (rule.bucket, rule.bucket))
        .collect();
    assert_eq!(listed(ROW_RULES), row_rules);
    assert_eq!(
        row_rules
            .iter()
            .map(|(bucket, _)| *bucket)
            .collect::<Vec<_>>(),
        [
            STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON,
            STOP_CODON_COMPLETED_BY_POLYADENYLATION,
            STOP_CODON_OVERLAP_WITHOUT_STOP_CODON,
            STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON,
            START_CODON_OVERLAP_WITHOUT_START_CODON,
            STOP_CODON_OVERLAP_GENOMIC_SPAN,
            STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
            STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
            STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
            STOP_CODON_OVERLAP_INSERTION_RANGE,
            STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES,
            FRAMESHIFT_INSIDE_A_HELD_STOP_CODON,
            SPLICE_TERM_PAST_THE_RECORD_SPAN,
            "pachysolen_nuclear_code_read_as_standard",
            START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR,
            START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR,
        ]
    );
    assert_eq!(pair_classes().len(), 22);
    for rule in EXCLUSION_REGISTRY.iter() {
        let registry = if rule.row_aware() {
            ROW_RULES
        } else {
            SET_RULES
        };
        assert!(!taxonomy_class(registry, rule.bucket).is_empty());
    }
}

#[test]
fn classify_pair_names_the_taxonomy_class_of_the_first_set_rule() {
    let family = Some("splice_family_swap");
    let start = Some("start_cooccurrence_swap");
    let cases = [
        // Each of the four splice buckets is the one taxonomy class.
        (
            "intron_variant",
            "intron_variant,splice_region_variant",
            "C",
            family,
        ),
        (
            "intron_variant,non_coding_transcript_variant",
            "intron_variant,non_coding_transcript_variant,splice_polypyrimidine_tract_variant",
            "C",
            family,
        ),
        (
            "frameshift_variant",
            "frameshift_variant,splice_region_variant",
            "-",
            family,
        ),
        (
            "missense_variant",
            "missense_variant,splice_region_variant",
            "T",
            family,
        ),
        (
            "splice_donor_variant,coding_sequence_variant",
            "splice_donor_variant,coding_sequence_variant,splice_region_variant,intron_variant",
            "C",
            family,
        ),
        // The start rule reads the allele: a sequence allele keeps `start_retained_variant`,
        // a structural one `start_lost`.
        (
            "start_lost,start_retained_variant",
            "start_retained_variant",
            "T",
            start,
        ),
        ("start_lost,start_retained_variant", "start_lost", "T", None),
        ("start_lost,start_retained_variant", "start_lost", "-", None),
        (
            "start_lost,start_retained_variant",
            "start_lost",
            "<DEL>",
            start,
        ),
        (
            "start_lost,start_retained_variant",
            "start_lost",
            "deletion",
            start,
        ),
        (
            "start_lost,start_retained_variant,5_prime_UTR_variant",
            "start_retained_variant,5_prime_UTR_variant",
            "T",
            start,
        ),
        // The terms the other way round, another difference, a pair of a row-aware rule and
        // equal sets take no class.
        (
            "intron_variant,splice_region_variant",
            "intron_variant",
            "C",
            None,
        ),
        ("missense_variant", "synonymous_variant", "T", None),
        (
            "coding_sequence_variant",
            "stop_retained_variant",
            "TGA",
            None,
        ),
        ("intron_variant", "intron_variant", "C", None),
        ("", "", "", None),
    ];
    for (perl, rust, allele, class) in cases {
        assert_eq!(
            classify_pair(perl, rust, allele),
            class,
            "{perl} | {rust} | {allele}"
        );
    }
}

/// A forward protein-coding transcript on chromosome 21 from 33031935 to 33041244 whose
/// first intron is 33032155-33036102, filed as `ENST1`.
fn models() -> TranscriptModels {
    let mut models = TranscriptModels::default();
    let model = TranscriptModel {
        chrom: "21".to_string(),
        start: 33031935,
        end: 33041244,
        strand: 1,
        biotype: "protein_coding".to_string(),
        exons: vec![(33031935, 33032154), (33036103, 33036199)],
        cds: Some((33031984, 33036182)),
        mirna_mature: Vec::new(),
    };
    models.insert("21", "ENST1", model);
    models
}

#[test]
fn classify_mate_side_pair_names_the_local_read_alone() {
    let models = models();
    let (location, allele) = ("1:14374334", "]21:33034355]A");
    let (local, mate) = ("feature_truncation", "feature_truncation,intron_variant");
    let class = |location, allele, feature, perl, rust| {
        classify_mate_side_pair(location, allele, feature, perl, rust, &models)
    };
    let named = Some(BREAKEND_MATE_LOCAL_READ);
    assert_eq!(named, Some("breakend_mate_local_read"));
    assert_eq!(class(location, allele, "ENST1", local, mate), named);
    assert_eq!(
        class(location, "]chr21:33034355]A", "ENST1", local, mate),
        named
    );
    assert_eq!(class("chr1:14374334", allele, "ENST1", local, mate), named);
    // The terms of a set in any order.
    let unordered = "intron_variant,feature_truncation";
    assert_eq!(class(location, allele, "ENST1", local, unordered), named);
    // No class: the sets the other way round, another set on either side, a transcript the
    // cache does not hold or holds under another chromosome, an allele that names no mate, a
    // ranged Location, a local coordinate that selects the transcript itself, a mate past
    // the distance gate, and a row whose two reads are one set.
    assert_eq!(class(location, allele, "ENST1", mate, local), None);
    assert_eq!(class(location, allele, "ENST1", local, local), None);
    assert_eq!(
        class(location, allele, "ENST1", "intergenic_variant", mate),
        None
    );
    assert_eq!(class(location, allele, "ENST2", local, mate), None);
    assert_eq!(
        class(location, "]22:33034355]A", "ENST1", local, mate),
        None
    );
    assert_eq!(class(location, "deletion", "ENST1", local, mate), None);
    assert_eq!(
        class("1:14374334-14374335", allele, "ENST1", local, mate),
        None
    );
    assert_eq!(class("21:33034000", allele, "ENST1", local, mate), None);
    assert_eq!(
        class(location, "]21:33050000]A", "ENST1", local, mate),
        None
    );
    let same = "downstream_gene_variant";
    assert_eq!(
        class("1:33046244", "]21:33046244]A", "ENST1", same, same),
        None
    );
    assert_eq!(
        classify_mate_side_pair(
            location,
            allele,
            "ENST1",
            local,
            mate,
            &TranscriptModels::default()
        ),
        None
    );
}

#[test]
fn record_reference_reads_a_row_by_position_from_the_records_that_write_it() {
    let replaced = run_record("chr21\t9000\trs1\tAC\tGTTCA\t.\tPASS\t.");
    let inserted = run_record("21\t9000\trs2\tA\tAGTTCA\t.\tPASS\t.");
    let (location, allele) = ("21:9000-9001", "GTTCA");
    let one = std::slice::from_ref;
    let read =
        |records: &[RunRecord], location: &str| record_reference(records, location, allele, "rs9");
    assert_eq!(read(one(&replaced), location).as_deref(), Some("AC"));
    assert_eq!(read(one(&inserted), "chr21:9000-9001").as_deref(), Some(""));
    assert_eq!(read(one(&inserted), "CHR21:9000-9001").as_deref(), Some(""));
    // Two records that read the row differently give nothing under an identifier neither
    // carries, and under the identifier of one give that record's reading; two that read
    // it the same give their reading.
    let both = [replaced.clone(), inserted.clone()];
    assert_eq!(read(&both, location), None);
    assert_eq!(
        record_reference(&both, location, allele, "rs1").as_deref(),
        Some("AC")
    );
    assert_eq!(
        record_reference(&both, location, allele, "rs2").as_deref(),
        Some("")
    );
    assert_eq!(
        read(&[replaced.clone(), replaced.clone()], location).as_deref(),
        Some("AC")
    );
    // Two records of the row's own identifier that read it differently give nothing; a
    // record carries the identifier that its whole ID column spells, and none under `.`.
    let renamed = run_record("21\t9000\trs1\tA\tAGTTCA\t.\tPASS\t.");
    let own = [replaced.clone(), renamed];
    assert_eq!(record_reference(&own, location, allele, "rs1"), None);
    let listed_twice = run_record("21\t9000\trs7;rs8\tAC\tGTTCA\t.\tPASS\t.");
    let beside = [listed_twice, inserted.clone()];
    assert_eq!(
        record_reference(&beside, location, allele, "rs7;rs8").as_deref(),
        Some("AC")
    );
    assert_eq!(record_reference(&beside, location, allele, "rs7"), None);
    let unnamed = run_record("21\t9000\t.\tAC\tGTTCA\t.\tPASS\t.");
    assert_eq!(
        record_reference(&[unnamed, inserted.clone()], location, allele, "."),
        None
    );
    // No record writes the row: another position, another chromosome, a Location without a
    // chromosome, no record at all.
    assert_eq!(read(one(&replaced), "21:9100-9101"), None);
    assert_eq!(read(one(&replaced), "22:9000-9001"), None);
    assert_eq!(read(one(&replaced), "9000-9001"), None);
    assert_eq!(read(&[], location), None);
    // A line of fewer than five columns and a POS vep-rs does not read are no records; a
    // POS it reads in another spelling is one, and so is a line that ends at its ALT.
    let short = run_record("21\t9000\trs1\tAC");
    let unread = run_record("21\t9_000\trs1\tAC\tGTTCA\t.\tPASS\t.");
    assert_eq!(read(&[short, unread], location), None);
    for pos in ["9000.0", "+9000", " 9000", "9e3"] {
        let record = run_record(&format!("21\t{pos}\trs1\tAC\tGTTCA"));
        assert_eq!(read(one(&record), location).as_deref(), Some("AC"), "{pos}");
    }
    // The ALT column is split as Perl splits it: an empty last field names no allele.
    let listed = run_record("21\t9000\trs1\tAC\tT,GTTCA,\t.\tPASS\t.");
    assert_eq!(read(one(&listed), location).as_deref(), Some("AC"));
    let lower = run_record("21\t9000\trs1\tac\tgttca\t.\tPASS\t.");
    assert_eq!(read(one(&lower), location).as_deref(), Some("AC"));
    assert_eq!(
        record_reference(one(&replaced), location, "gttca", "rs9"),
        None
    );
}

#[test]
fn a_row_that_names_no_feature_is_read_by_its_terms_and_the_manifests_own_run() {
    let (location, allele) = ("21:101-105", "-");
    let untested: HashSet<(String, String)> =
        HashSet::from([(location.to_string(), allele.to_string())]);
    let retained: HashSet<String> = HashSet::from([location.to_string()]);
    let nothing = no_records();
    let class = |perl: &str, rust: &str, untested, retained| {
        classify_stop_codon_rows(
            location,
            allele,
            "-",
            perl,
            rust,
            &[],
            &[],
            untested,
            retained,
            &[],
        )
    };
    let shared = "splice_acceptor_variant,intron_variant";
    let with = |term: &str| {
        if term.is_empty() {
            shared.to_string()
        } else {
            format!("{term},{shared}")
        }
    };
    let untested_class = Some(STOP_RETAINED_WITHOUT_STOP_CODON_TEST);
    let removed_class = Some(STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST);
    let insertion_class = Some(STOP_CODON_OVERLAP_INSERTION_RANGE);
    // (the reference's terms, vep-rs's terms, the class on a record documented under an
    // untested class, the class on a record documented under the insertion class)
    let cases = [
        ("stop_retained_variant", "", untested_class, None),
        (
            "stop_retained_variant",
            "coding_sequence_variant",
            untested_class,
            None,
        ),
        (
            "stop_retained_variant",
            "frameshift_variant",
            untested_class,
            None,
        ),
        ("stop_retained_variant", "stop_lost", removed_class, None),
        (
            "stop_retained_variant",
            "stop_lost,frameshift_variant",
            None,
            None,
        ),
        ("stop_retained_variant", "inframe_deletion", None, None),
        (
            "stop_retained_variant,stop_gained",
            "frameshift_variant",
            None,
            None,
        ),
        (
            "coding_sequence_variant",
            "stop_retained_variant",
            None,
            insertion_class,
        ),
        (
            "coding_sequence_variant",
            "inframe_insertion,stop_retained_variant",
            None,
            None,
        ),
        (
            "frameshift_variant",
            "stop_retained_variant",
            None,
            insertion_class,
        ),
        (
            "frameshift_variant",
            "inframe_insertion,stop_retained_variant",
            None,
            insertion_class,
        ),
        (
            "stop_gained,frameshift_variant",
            "inframe_insertion,stop_retained_variant",
            None,
            insertion_class,
        ),
        (
            "frameshift_variant",
            "protein_altering_variant,stop_retained_variant",
            None,
            None,
        ),
        ("", "", None, None),
    ];
    for (perl, rust, on_untested, on_retained) in cases {
        let (perl, rust) = (with(perl), with(rust));
        let label = format!("{perl} | {rust}");
        assert_eq!(
            class(&perl, &rust, &untested, &nothing.1),
            on_untested,
            "{label}"
        );
        assert_eq!(
            class(&perl, &rust, &nothing.0, &retained),
            on_retained,
            "{label}"
        );
        assert_eq!(class(&perl, &rust, &nothing.0, &nothing.1), None, "{label}");
        let both = on_retained.or(on_untested);
        assert_eq!(class(&perl, &rust, &untested, &retained), both, "{label}");
    }
    // The untested classes are the record's, its Location and its Allele; the insertion
    // class is the Location's.
    let other: HashSet<(String, String)> = HashSet::from([(location.to_string(), "A".to_string())]);
    let (perl, rust) = (with("stop_retained_variant"), with(""));
    assert_eq!(class(&perl, &rust, &other, &nothing.1), None);
}

#[test]
fn a_transcript_row_is_asked_of_the_row_aware_rules_and_another_row_of_none() {
    let (untested, retained) = no_records();
    let (location, allele) = ("21:7001-7007", "-");
    let (perl, rust) = ("stop_retained_variant", "coding_sequence_variant");
    let class = |feature_type: &str, reference, candidate| {
        classify_stop_codon_pair(
            location,
            allele,
            feature_type,
            perl,
            rust,
            reference,
            &untested,
            &retained,
            candidate,
            &[],
        )
    };
    let cut = PositionRow::of_positions("?-7057", "-", "rs7000");
    let whole = PositionRow::of_positions("7051-7057", "6974-6980", "rs7000");
    let (cut, whole) = (Some(&cut), Some(&whole));
    let untested_class = Some(STOP_RETAINED_WITHOUT_STOP_CODON_TEST);
    assert_eq!(class("Transcript", cut, cut), untested_class);
    assert_eq!(class("Transcript", cut, whole), None);
    // A pair without the positions of either row is left unclassed: the reference's do
    // not stand in for vep-rs's, which here would name the class the pair of `whole` lacks.
    assert_eq!(class("Transcript", cut, None), None);
    assert_eq!(class("Transcript", None, cut), None);
    assert_eq!(class("Transcript", None, None), None);
    assert_eq!(class("RegulatoryFeature", cut, cut), None);
    assert_eq!(class("MotifFeature", cut, cut), None);
    assert_eq!(class("", cut, cut), None);
}

#[test]
fn the_rows_of_a_tuple_stand_or_fall_together() {
    let (untested, retained) = no_records();
    let (location, allele) = ("21:7001-7007", "-");
    let (perl, rust) = ("stop_retained_variant", "coding_sequence_variant");
    let row = |cdna: &str, cds: &str| PositionRow::of_positions(cdna, cds, "rs7000");
    let (cut, whole) = (row("?-7057", "-"), row("7051-7057", "6974-6980"));
    let class = |reference: &[PositionRow], candidate: &[PositionRow]| {
        classify_stop_codon_rows(
            location,
            allele,
            "Transcript",
            perl,
            rust,
            reference,
            candidate,
            &untested,
            &retained,
            &[],
        )
    };
    let untested_class = Some(STOP_RETAINED_WITHOUT_STOP_CODON_TEST);
    let one = std::slice::from_ref;
    assert_eq!(class(one(&cut), one(&cut)), untested_class);
    assert_eq!(
        class(&[cut.clone(), cut.clone()], &[cut.clone(), cut.clone()]),
        untested_class
    );
    // Rows whose test is made on both sides are the three-base class's.
    let tested_class = Some(STOP_RETAINED_THREE_BASE_STOP_CODON_TEST);
    assert_eq!(class(one(&whole), one(&whole)), tested_class);
    // One pair of rows no class takes, on either side, leaves the tuple unnamed, and so do
    // two pairs of rows that two classes take; so does a side without a row.
    assert_eq!(class(&[cut.clone(), whole.clone()], one(&cut)), None);
    assert_eq!(class(one(&cut), &[cut.clone(), whole.clone()]), None);
    assert_eq!(class(one(&cut), one(&whole)), None);
    let both = [cut.clone(), whole.clone()];
    assert_eq!(class(&both, &both), None);
    assert_eq!(class(&[], one(&cut)), None);
    assert_eq!(class(one(&cut), &[]), None);
    assert_ne!(untested_class, Some(UNEXPLAINED));
}

#[test]
fn a_row_under_an_identifier_is_read_with_the_runs_records() {
    // Two replaced bases inside an exon take the splice term class with the record's
    // alleles, which the first column spells or the run's records give by position; a row
    // no record writes is read without them, and two positions then take no class. Beside
    // an insertion of the same allele at its Location the row is read with the record that
    // carries its identifier, and without the alleles where neither record carries it.
    let (untested, retained) = no_records();
    let (location, allele) = ("21:9000-9001", "GTTCA");
    let perl = "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant";
    let rust = "frameshift_variant,splice_region_variant";
    let rows = |name: &str| vec![PositionRow::of_positions("999-1000", "899-900", name)];
    let class = |name: &str, records: &[RunRecord]| {
        classify_stop_codon_rows(
            location,
            allele,
            "Transcript",
            perl,
            rust,
            &rows(name),
            &rows(name),
            &untested,
            &retained,
            records,
        )
    };
    let record = run_record("21\t9000\trs1\tAC\tGTTCA\t.\tPASS\t.");
    let elsewhere = run_record("21\t109000\trs1\tAC\tGTTCA\t.\tPASS\t.");
    let insertion = run_record("21\t9000\trs2\tA\tAGTTCA\t.\tPASS\t.");
    let named = Some(SPLICE_TERM_PAST_THE_RECORD_SPAN);
    let one = std::slice::from_ref;
    assert_eq!(class("rs1", one(&record)), named);
    assert_eq!(class("21_9000_AC/GTTCA", &[]), named);
    assert_eq!(class("21_9000_AC/GTTCA", one(&insertion)), named);
    assert_eq!(class("rs1", one(&elsewhere)), None);
    assert_eq!(class("rs1", &[]), None);
    assert_eq!(class("rs2", one(&insertion)), None);
    let both = [record.clone(), insertion];
    assert_eq!(class("rs1", &both), named);
    assert_eq!(class("rs2", &both), None);
    assert_eq!(class("rs9", &both), None);
}

#[test]
fn a_row_without_a_feature_takes_the_term_pairs_the_row_rules_take() {
    // The branch for a row that names no feature holds its own term sets, where a
    // Transcript row is asked of the tuple comparator's row rules. Over every way the terms
    // either reads can differ between the two sets, beside one shared term, the row without
    // a feature takes a class exactly where that class's row rule takes a Transcript row
    // that meets the rule's other conditions: a deletion's row whose stop codon test is not
    // made for the two untested classes, an insertion's row with every position for the
    // insertion class. The shared term is not `splice_region_variant`, beside which the
    // insertion rule leaves `coding_sequence_variant` against `stop_retained_variant` to the
    // genomic span class, which no row without a feature takes.
    let read = [
        STOP_RETAINED_TERM,
        CODING_SEQUENCE_TERM,
        FRAMESHIFT_TERM,
        STOP_LOST_TERM,
        "stop_gained",
        "inframe_insertion",
        "inframe_deletion",
        "protein_altering_variant",
        "missense_variant",
    ];
    let nothing = no_records();
    let deletion = ("21:7001-7007", "-");
    let insertion = ("21:7001-7002", "A");
    let untested = HashSet::from([(deletion.0.to_string(), deletion.1.to_string())]);
    let retained = HashSet::from([insertion.0.to_string()]);
    let cut = ("?-7057", "-");
    let whole = ("7051-7052", "6974-6975");
    let untested_classes = [
        STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
        STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
    ];
    let mut taken = [0usize; 3];
    for pair in 0..3usize.pow(read.len() as u32) {
        let (mut perl, mut rust) = (vec!["synonymous_variant"], vec!["synonymous_variant"]);
        let mut side = pair;
        for term in read {
            match side % 3 {
                1 => perl.push(term),
                2 => rust.push(term),
                _ => {}
            }
            side /= 3;
        }
        let (perl, rust) = (perl.join(","), rust.join(","));
        let label = format!("{perl} | {rust}");
        let featureless = |row: (&str, &str), untested, retained| {
            classify_stop_codon_pair(
                row.0,
                row.1,
                "-",
                &perl,
                &rust,
                None,
                untested,
                retained,
                None,
                &[],
            )
        };
        let transcript = |row: (&str, &str), (cdna, cds): (&str, &str)| {
            let carried = PositionRow::of_positions(cdna, cds, "rs7000");
            comparator_stop_codon_class(row.0, row.1, &perl, &rust, &carried, &carried, &[])
        };
        let on_cut_rows =
            transcript(deletion, cut).filter(|class| untested_classes.contains(class));
        assert_eq!(
            featureless(deletion, &untested, &nothing.1),
            on_cut_rows,
            "{label}"
        );
        let on_whole_rows = transcript(insertion, whole)
            .filter(|class| *class == STOP_CODON_OVERLAP_INSERTION_RANGE);
        assert_eq!(
            featureless(insertion, &nothing.0, &retained),
            on_whole_rows,
            "{label}"
        );
        for (count, class) in taken.iter_mut().zip([on_cut_rows, on_whole_rows]) {
            *count += usize::from(class.is_some());
        }
        taken[2] += usize::from(on_cut_rows == Some(STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST));
    }
    // Three pairs of the first untested class and one of the second; five of the insertion
    // class.
    assert_eq!(taken, [4, 5, 1]);
}

#[test]
fn a_transcript_row_is_put_to_the_rules_with_its_residues_and_its_codons() {
    // `stop_codon_settled_before_unknown_bases` holds a pair to the Codons cell of both rows
    // and to their residues: the pair takes the class where the cell settles a stop codon
    // before its unknown bases and vep-rs's residues hold their stop at that codon's place,
    // and takes none without the cell, with a cell that settles no stop codon or with the
    // stop at another place.
    let (location, allele) = ("21:10567700-10567701", "GAATAGNNNN");
    let (perl, rust) = ("frameshift_variant", "frameshift_variant,stop_gained");
    let row = |amino_acids: &str, codons: &str| PositionRow {
        cells: Some(["-", amino_acids, codons].map(str::to_string)),
        ..PositionRow::of_positions("775-776", "423-424", "rs1")
    };
    let class = |reference: &PositionRow, candidate: &PositionRow| {
        comparator_stop_codon_class(location, allele, perl, rust, reference, candidate, &[])
    };
    let settled = "-/GAATAGNNNN";
    assert_eq!(
        class(&row("-", settled), &row("-/E*XX", settled)),
        Some(STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES)
    );
    for (reference, candidate) in [
        (row("-", "-"), row("-/E*XX", "-")),
        (row("-", settled), row("-/E*XX", "-")),
        (row("-", "-/GAAGCTNNNN"), row("-/E*XX", "-/GAAGCTNNNN")),
        (row("-", settled), row("-/EE*X", settled)),
        (row("-", settled), row("-", settled)),
        (row("-/E*XX", settled), row("-/E*XX", settled)),
    ] {
        let label = format!("{reference:?} | {candidate:?}");
        assert_eq!(class(&reference, &candidate), None, "{label}");
    }
    // A row that ends before its Codons column carries none of its cells.
    let short = PositionRow::of_positions("775-776", "423-424", "rs1");
    assert_eq!(class(&short, &short), None);
    assert_eq!(class(&row("-", settled), &short), None);
}

#[test]
fn a_frameshift_inside_a_held_stop_codon_is_named_by_the_residues_and_codons_of_both_rows() {
    // `frameshift_inside_a_held_stop_codon` holds a pair to the residues and the Codons cell
    // of both rows: the pair takes the class where both rows hold the same residues, a stop
    // first and the unknown residue after the `/`, and the same Codons cell, whose halves
    // differ by a number of bases that is no multiple of three. It takes none without the
    // cells, with residues that differ, open with no stop or hold no unknown residue, or
    // with halves of one length or cells that differ.
    let (location, allele) = ("21:12001", "-");
    let (perl, rust) = ("coding_sequence_variant", "frameshift_variant");
    let row = |amino_acids: &str, codons: &str| PositionRow {
        cells: Some(["-", amino_acids, codons].map(str::to_string)),
        ..PositionRow::of_positions("36", "30", "rs1")
    };
    let class = |reference: &PositionRow, candidate: &PositionRow| {
        comparator_stop_codon_class(location, allele, perl, rust, reference, candidate, &[])
    };
    let lost = "tgA/tg";
    assert_eq!(
        class(&row("*/X", lost), &row("*/X", lost)),
        Some(FRAMESHIFT_INSIDE_A_HELD_STOP_CODON)
    );
    for (reference, candidate) in [
        (row("-", "-"), row("-", "-")),
        (row("*/X", "-"), row("*/X", "-")),
        (row("*/X", lost), row("*/LX", lost)),
        (row("L/X", lost), row("L/X", lost)),
        (row("*/L", lost), row("*/L", lost)),
        (row("*/X", "tGa/tAa"), row("*/X", "tGa/tAa")),
        (row("*/X", lost), row("*/X", "tGa/ta")),
    ] {
        let label = format!("{reference:?} | {candidate:?}");
        assert_eq!(class(&reference, &candidate), None, "{label}");
    }
}

#[test]
fn a_lost_start_takes_its_class_by_the_positions_and_the_residues_of_both_rows() {
    // `start_codon_removed_without_five_prime_utr` holds a pair to the two position columns
    // and to the residues of both rows: the pair takes the class where the cDNA and the CDS
    // position open at 1 on both rows and neither row writes an alternate residue, and
    // takes none behind a 5' UTR, after the first codon, with positions that differ between
    // the rows, without residues or beside an alternate residue. A row that names no
    // feature carries neither and takes the class nowhere.
    let (location, allele) = ("21:44817966-44817968", "-");
    let (perl, rust) = ("inframe_deletion", "inframe_deletion,start_lost");
    let row = |cdna: &str, cds: &str, amino_acids: &str| PositionRow {
        cells: Some(["-", amino_acids, "ATG/-"].map(str::to_string)),
        ..PositionRow::of_positions(cdna, cds, "21_44817965_ACAT/A")
    };
    let class = |reference: &PositionRow, candidate: &PositionRow| {
        comparator_stop_codon_class(location, allele, perl, rust, reference, candidate, &[])
    };
    let first_codon = row("1-3", "1-3", "M/-");
    assert_eq!(
        class(&first_codon, &first_codon),
        Some(START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR)
    );
    for (reference, candidate) in [
        (row("101-103", "1-3", "M/-"), row("101-103", "1-3", "M/-")),
        (row("4-6", "4-6", "A/-"), row("4-6", "4-6", "A/-")),
        (row("1-3", "1-3", "M/-"), row("1-6", "1-3", "M/-")),
        (row("1-3", "1-3", "-"), row("1-3", "1-3", "M/-")),
        (row("1-3", "1-3", "M/-"), row("1-3", "1-3", "M/X")),
    ] {
        let label = format!("{reference:?} | {candidate:?}");
        assert_eq!(class(&reference, &candidate), None, "{label}");
    }
    let (untested, retained) = no_records();
    let featureless = classify_stop_codon_pair(
        location,
        allele,
        "-",
        perl,
        rust,
        None,
        &untested,
        &retained,
        None,
        &[],
    );
    assert_eq!(featureless, None);
}

#[test]
fn a_kept_start_takes_its_class_by_the_terms_and_the_positions_of_both_rows() {
    // `start_codon_kept_without_five_prime_utr` holds a pair to the two position columns of
    // both rows and to the two forms of its term difference: the pair takes the class where
    // the cDNA and the CDS position are one span that opens on coding base 1, 2 or 3 on both
    // rows, with the reference's `start_lost` against vep-rs's `start_retained_variant` or
    // with that term added to a set without a start term. It takes none behind a 5' UTR,
    // on a span that opens on the fourth coding base, with positions that differ between
    // the rows, beside another differing term or where vep-rs keeps `start_lost`. A row
    // that names no feature carries no position and takes the class nowhere.
    let (location, allele) = ("21:30001-30010", "-");
    let lost = "frameshift_variant,start_lost";
    let kept = "frameshift_variant,start_retained_variant";
    let row = |cdna: &str, cds: &str| PositionRow::of_positions(cdna, cds, "rs30000");
    let class = |perl: &str, rust: &str, reference: &PositionRow, candidate: &PositionRow| {
        comparator_stop_codon_class(location, allele, perl, rust, reference, candidate, &[])
    };
    let on_the_start_codon = row("3-12", "3-12");
    let taken = Some(START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR);
    for (perl, rust) in [
        (lost, kept),
        ("frameshift_variant", kept),
        (
            "inframe_deletion,start_lost",
            "inframe_deletion,start_retained_variant",
        ),
        (
            "inframe_insertion",
            "inframe_insertion,start_retained_variant",
        ),
    ] {
        let named = class(perl, rust, &on_the_start_codon, &on_the_start_codon);
        assert_eq!(named, taken, "{perl} | {rust}");
    }
    for (reference, candidate) in [
        (row("104-113", "3-12"), row("104-113", "3-12")),
        (row("4-13", "4-13"), row("4-13", "4-13")),
        (row("3-12", "3-12"), row("3-13", "3-12")),
        (row("3-12", "3-12"), row("3-12", "3-13")),
        (row("3-12", "-"), row("3-12", "-")),
    ] {
        let label = format!("{reference:?} | {candidate:?}");
        assert_eq!(class(lost, kept, &reference, &candidate), None, "{label}");
    }
    for (perl, rust) in [
        (lost, "frameshift_variant,start_lost,start_retained_variant"),
        (
            lost,
            "NMD_transcript_variant,frameshift_variant,start_retained_variant",
        ),
        ("start_lost", "start_retained_variant"),
        (kept, lost),
    ] {
        let named = class(perl, rust, &on_the_start_codon, &on_the_start_codon);
        assert_eq!(named, None, "{perl} | {rust}");
    }
    let (untested, retained) = no_records();
    let featureless = classify_stop_codon_pair(
        location,
        allele,
        "-",
        lost,
        kept,
        None,
        &untested,
        &retained,
        None,
        &[],
    );
    assert_eq!(featureless, None);
}

#[test]
fn a_rule_of_a_class_no_corpus_documents_is_asked_of_no_pair() {
    // The rules asked are those of the classes of `CORRECTED_BEHAVIOUR_CLASSES`: the
    // registry's row-aware rules but the rule of `pachysolen_nuclear_code_read_as_standard`,
    // which takes the rows of its own pair and names no class here.
    let outside = "pachysolen_nuclear_code_read_as_standard";
    let asked: Vec<&str> = row_pair_rules()
        .iter()
        .map(|rule| taxonomy_class(ROW_RULES, rule.bucket))
        .filter(|class| CORRECTED_BEHAVIOUR_CLASSES.contains(class))
        .collect();
    assert_eq!(
        asked,
        [
            STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON,
            STOP_CODON_COMPLETED_BY_POLYADENYLATION,
            STOP_CODON_OVERLAP_WITHOUT_STOP_CODON,
            STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON,
            START_CODON_OVERLAP_WITHOUT_START_CODON,
            STOP_CODON_OVERLAP_GENOMIC_SPAN,
            STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
            STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
            STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
            STOP_CODON_OVERLAP_INSERTION_RANGE,
            STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES,
            FRAMESHIFT_INSIDE_A_HELD_STOP_CODON,
            SPLICE_TERM_PAST_THE_RECORD_SPAN,
            START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR,
            START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR,
        ]
    );
    let (location, allele) = ("7:1000", "A");
    let (perl, rust) = ("synonymous_variant", "missense_variant");
    let row = |amino_acids: &str| PositionRow {
        cells: Some(["-", amino_acids, "ctG/ctA"].map(str::to_string)),
        ..PositionRow::of_positions("114", "114", "rs1")
    };
    let (reference, candidate) = (row("L"), row("A/L"));
    let read = |csq: &str, row: &PositionRow| {
        let [protein_position, amino_acids, codons] = row.cells.as_ref().unwrap();
        [pair_row_with_codons(&[
            "rs1",
            location,
            allele,
            "",
            "",
            "Transcript",
            csq,
            &row.cdna_position,
            &row.cds_position,
            protein_position,
            amino_acids,
            codons,
        ])]
    };
    let taking: Vec<&str> = row_pair_rules()
        .iter()
        .filter(|rule| rule.matches_rows(&read(perl, &reference), &read(rust, &candidate)))
        .map(|rule| rule.bucket)
        .collect();
    assert_eq!(taking, [outside]);
    assert_eq!(
        comparator_stop_codon_class(location, allele, perl, rust, &reference, &candidate, &[]),
        None
    );
}
