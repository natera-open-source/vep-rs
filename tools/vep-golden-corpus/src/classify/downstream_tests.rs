// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the Downstream module: the genetic codes against the engine's own list,
//! translation, the two readers of a made-up cache, the span of a record, the allele a
//! shift turns and the coding model.

use super::*;
use crate::TempDir;
use regex::Regex;
use std::path::Path;

/// The transcript the Downstream plugin is tested on in
/// crates/vep-builtins/src/plugins/downstream.rs: eight codons and a stop, and a 3' UTR
/// that reads three codons and a stop from its first base.
const CODING: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";
const UTR: &str = "GGGCCCAAATAGCCC";
/// A second transcript of five codons and a stop.
const STOP_CODING: &str = "ATGTACGGAAAACTGTAA";
/// A third transcript of four codons and a stop under either of the first two tables.
const TABLE_TWO_CODING: &str = "ATGAAACTGACCTAG";

/// Mapper pairs (first and last cDNA base, first and last genomic base, orientation) that
/// lay a transcript of three 5' UTR bases, a coding sequence and the UTR on the genome
/// from base 1001, on one exon of either strand.
const ONE_EXON: [[i64; 5]; 1] = [[1, 45, 1001, 1045, 1]];
const REVERSE_EXON: [[i64; 5]; 1] = [[1, 45, 1001, 1045, -1]];
/// Two exons around a two-base intron after coding base 12.
const TWO_EXONS: [[i64; 5]; 2] = [[1, 15, 1001, 1015, 1], [16, 45, 1018, 1047, 1]];
/// Two pairs that touch on the genome between bases 1015 and 1016 while the model carries
/// coding base 13 between them, and two that touch between bases 1014 and 1015 around
/// coding bases 12 and 13.
const ONE_BASE_SEAM: [[i64; 5]; 2] = [[1, 15, 1001, 1015, 1], [17, 45, 1016, 1044, 1]];
const TWO_BASE_SEAM: [[i64; 5]; 2] = [[1, 14, 1001, 1014, 1], [17, 45, 1015, 1043, 1]];

fn translated(sequence: &str, table: i64) -> String {
    translate_codons(sequence, table).unwrap()
}

/// A cache whose one shard, of chromosome 21, holds the records of the JSON list given.
fn cache_of(tag: &str, records: &str) -> (TempDir, String) {
    let scratch =
        TempDir::create(&std::env::temp_dir(), &format!("golden_downstream_{tag}_")).unwrap();
    let shard = scratch.path.join("transcripts/21");
    std::fs::create_dir_all(&shard).unwrap();
    std::fs::write(shard.join("1-1000000.json"), records).unwrap();
    let cache_dir = scratch.path.to_str().unwrap().to_string();
    (scratch, cache_dir)
}

#[test]
fn codon_table_amino_acids_are_the_codon_modules_genetic_codes() {
    // The table carries, id for id, the amino-acid string of every table `CODES` of
    // crates/vep-core/src/codon.rs holds, and none for an id that list leaves empty (0, 7,
    // 8 and 17 to 20).
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../crates/vep-core/src/codon.rs");
    let source = std::fs::read_to_string(path).unwrap();
    let block = &source[source.find("const CODES:").unwrap()..];
    let block = &block[..block.find("];").unwrap()];
    let entry = Regex::new(r#"(?m)^    (?:None|code\(\s*b"([^"]{64})")"#).unwrap();
    let entries: Vec<&str> = entry
        .captures_iter(block)
        .map(|found| found.get(1).map_or("", |amino_acids| amino_acids.as_str()))
        .collect();
    assert_eq!(entries.len(), 25);
    let listed: Vec<(i64, &str)> = (0i64..)
        .zip(entries)
        .filter(|(_, amino_acids)| !amino_acids.is_empty())
        .collect();
    assert_eq!(listed, CODON_TABLE_AMINO_ACIDS);
    let ids: Vec<i64> = CODON_TABLE_AMINO_ACIDS.iter().map(|(id, _)| *id).collect();
    assert_eq!(
        ids,
        [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24]
    );
    for (id, amino_acids) in CODON_TABLE_AMINO_ACIDS {
        assert_eq!(amino_acids.len(), 64);
        assert_eq!(genetic_code(id), Some(amino_acids));
    }
    assert_eq!(genetic_code(0), None);
    assert_eq!(genetic_code(7), None);
    assert_eq!(genetic_code(25), None);
}

#[test]
fn translate_codons_reads_codons_as_the_codon_table_module_does() {
    // Whole codons are read in order under the named table and a trailing one or two bases
    // are dropped; a stop stays `*`; `U` reads as `T` and lower case as upper; a codon with
    // an ambiguous letter is the one amino acid all its readings code under that table, `B`
    // for aspartate or asparagine, `Z` for glutamate or glutamine and `X` otherwise; three
    // hyphens are a gap, and a codon with any other character that is no nucleotide letter
    // is `X`. Table 2 reads four codons otherwise than the standard table: a stop as
    // tryptophan, arginine twice as a stop and isoleucine as methionine.
    assert_eq!(translated(CODING, 1), "MKTLAVRG*");
    assert_eq!(translated(STOP_CODING, 1), "MYGKL*");
    assert_eq!(translated(TABLE_TWO_CODING, 1), "MKLT*");
    assert_eq!(translated(TABLE_TWO_CODING, 2), "MKLT*");
    assert_eq!(translated("ATGA", 1), "M");
    assert_eq!(translated("ATGAA", 1), "M");
    assert_eq!(translated("AT", 1), "");
    assert_eq!(translated("augCTGtaa", 1), "ML*");
    let each = |codons: &[&str], table: i64| -> Vec<String> {
        codons
            .iter()
            .map(|codon| translated(codon, table))
            .collect()
    };
    assert_eq!(
        each(&["CTN", "TAR", "TRA", "RAY", "SAA", "NAT", "NNN"], 1),
        ["L", "*", "*", "B", "Z", "X", "X"]
    );
    assert_eq!(each(&["---", "A-G", "A.G"], 1), ["-", "X", "X"]);
    let codons = ["TGA", "AGA", "AGG", "ATA", "TGR", "AGR"];
    assert_eq!(each(&codons, 1), ["*", "R", "R", "I", "X", "R"]);
    assert_eq!(each(&codons, 2), ["W", "*", "*", "M", "W", "*"]);
    assert_eq!(translate_codons(CODING, 7), None);
}

#[test]
fn downstream_sequences_reads_the_coding_transcripts_of_the_corpus_cache() {
    // Each transcript carrying a `translateable_seq` comes with that sequence, its
    // `three_prime_utr` (the empty sequence when the cache spells none) and its
    // `codon_table`, a number or the converted cache's text, 1 when the cache names none or
    // 0; a transcript with no coding sequence, or whose table is no number, is left out.
    let coding = format!(r#""translateable_seq": "{CODING}""#);
    let record = |stable_id: &str, members: &str| {
        format!(
            r#"{{"stable_id": "{stable_id}", "variation_effect_feature_cache": {{{members}}}}}"#
        )
    };
    let records = [
        record(
            "ENST1",
            &format!(r#"{coding}, "three_prime_utr": "{UTR}", "codon_table": 1"#),
        ),
        record("ENST2", &format!(r#"{coding}, "three_prime_utr": null"#)),
        record("ENST3", &format!(r#"{coding}, "codon_table": "2""#)),
        record("ENST4", &format!(r#"{coding}, "codon_table": 0"#)),
        record("ENST5", &format!(r#"{coding}, "codon_table": "0""#)),
        record(
            "ENST6",
            &format!(r#"{coding}, "codon_table": "vertebrate""#),
        ),
        record("ENST7", &format!(r#""three_prime_utr": "{UTR}""#)),
        r#"{"stable_id": "ENST8"}"#.to_string(),
    ];
    let (_scratch, cache_dir) = cache_of("sequences", &format!("[{}]", records.join(", ")));
    let held = |utr: &str, codon_table: i64| CodingSequence {
        coding: CODING.to_string(),
        utr: utr.to_string(),
        codon_table,
    };
    let expected: CodingSequences = [
        ("ENST1", held(UTR, 1)),
        ("ENST2", held("", 1)),
        ("ENST3", held("", 2)),
        ("ENST4", held("", 1)),
        ("ENST5", held("", 1)),
    ]
    .into_iter()
    .map(|(stable_id, sequence)| (stable_id.to_string(), sequence))
    .collect();
    assert_eq!(downstream_sequences(&cache_dir).unwrap(), expected);
}

#[test]
fn record_span_is_the_span_the_vcf_parser_gives_a_record() {
    // A record's start is POS and its end the last REF base; where an ALT differs from REF
    // in length and REF and the ALTs open with one base, the start moves one base on
    // (Parser/VCF.pm:295-336), so an insertion's start lies one past its end. An ALT naming
    // `*` is left out of the comparison among several alleles; one ALT opening with another
    // base, or a record whose alleles all have REF's length, keeps POS. `row_record_span`
    // picks the one record span a Location prints, an insertion's two flanks lower first.
    let record = |reference: &str, alternates: &str| {
        format!("21\t100\t.\t{reference}\t{alternates}\t.\t.\t.")
    };
    let span = |reference: &str, alternates: &str| record_span(&record(reference, alternates));
    assert_eq!(span("AT", "A").unwrap(), (101, 101));
    assert_eq!(span("ATTT", "A").unwrap(), (101, 103));
    assert_eq!(span("A", "ATT").unwrap(), (101, 100));
    assert_eq!(span("A", "G").unwrap(), (100, 100));
    assert_eq!(span("AT", "GC").unwrap(), (100, 101));
    assert_eq!(span("AT", "G").unwrap(), (100, 101));
    assert_eq!(span("AT", "A,ATT").unwrap(), (101, 101));
    assert_eq!(span("AT", "A,*").unwrap(), (101, 101));
    assert_eq!(span("AT", "A,G").unwrap(), (100, 101));
    assert_eq!(span("A", "G,T").unwrap(), (100, 100));
    let (insertion, deletion, two_bases) =
        (record("A", "ATT"), record("AT", "A"), record("AT", "GC"));
    let row_span = |location: &str, lines: &[&String]| {
        row_record_span(location, lines.iter().map(|line| line.as_str())).unwrap()
    };
    assert_eq!(row_span("21:100-101", &[&insertion]), Some((101, 100)));
    assert_eq!(row_span("21:100-101", &[&two_bases]), Some((100, 101)));
    assert_eq!(
        row_span("21:101", &[&insertion, &deletion]),
        Some((101, 101))
    );
    assert_eq!(
        row_span("21:101", &[&deletion, &deletion]),
        Some((101, 101))
    );
    assert_eq!(row_span("21:100-101", &[&insertion, &two_bases]), None);
    assert_eq!(row_span("21:102", &[&insertion, &deletion]), None);
    assert_eq!(row_span("21:101", &[]), None);
}

#[test]
fn a_record_line_without_a_span_is_refused() {
    let short = record_span("21\t100\t.\tA").unwrap_err();
    assert_eq!(
        short.line,
        "ERROR: [vep-golden-corpus] the record line '21\\t100\\t.\\tA' has fewer than five columns"
    );
    let unnumbered = record_span("21\tx\trs1\tA\tG").unwrap_err();
    assert_eq!(
        unnumbered.line,
        "ERROR: [vep-golden-corpus] the record 21 x rs1 has a POS that is no integer of 64 bits"
    );
    let last = record_span("21\t9223372036854775807\trs1\tAT\tA").unwrap_err();
    assert_eq!(
        last.line,
        "ERROR: [vep-golden-corpus] the record 21 9223372036854775807 rs1 has a span 64 bits do not hold"
    );
    assert_eq!((short.code, unnumbered.code, last.code), (1, 1, 1));
    // Every line is read, the one whose span the Location does not print as well.
    assert!(row_record_span("21:100", ["21\t100\t.\tA\tG", "21\t200"]).is_err());
}

#[test]
fn hgvs_shifted_allele_turns_the_written_allele_as_the_shift_does() {
    // `shift_feature_seqs` moves the first base of the allele as written to its end once
    // per base of the shift on a forward-strand transcript, and the allele's length less
    // the shift times on a reverse-strand one, not at all where that is not positive
    // (TranscriptVariationAllele.pm:1330-1345); a deletion has no bases to turn.
    let turned = |offsets: &[i128], reverse: bool| -> Vec<String> {
        offsets
            .iter()
            .map(|offset| hgvs_shifted_allele("GTGA", *offset, reverse))
            .collect()
    };
    assert_eq!(
        turned(&[1, 2, 4, 5, -1], false),
        ["TGAG", "GAGT", "GTGA", "TGAG", "TGAG"]
    );
    assert_eq!(
        turned(&[-1, -3, -4, -5, 1], true),
        ["AGTG", "TGAG", "GTGA", "GTGA", "AGTG"]
    );
    assert_eq!(hgvs_shifted_allele("", 3, false), "");
    assert_eq!(hgvs_shifted_allele("", -3, true), "");
    assert_eq!(hgvs_shifted_allele("GTGA", i128::MIN, false), "GTGA");
    assert_eq!(hgvs_shifted_allele("GTGA", i128::MIN, true), "GTGA");
}

#[test]
fn coding_model_places_a_genomic_span_on_the_coding_sequence() {
    // `cds_pair` gives a span of bases the CDS coordinates of its two ends, in the
    // transcript's order on either strand, and an insertion the bases between its flanks'
    // coordinates (Mapper.pm:504-558): a start one past its end between neighbours, the
    // bases the model carries where two pairs touch on the genome, and the place against
    // its one mapped flank at an exon's edge. A span with an end in an intron, a UTR or
    // outside the transcript has no pair, nor has an insertion with a flank in a UTR where
    // the model's bases between its flanks open or close the coding sequence
    // (TranscriptMapper.pm:441-443), and the start exon's phase moves every coordinate by
    // the bases it pads.
    let spanning = |pairs: &[[i64; 5]], coding_start: i64, coding_end: i64| CodingModel {
        pairs: pairs.to_vec(),
        coding_start,
        coding_end,
        peptide_length: Some(8),
    };
    let model = |pairs: &[[i64; 5]]| spanning(pairs, 4, 30);
    let (forward, reverse, exons) = (model(&ONE_EXON), model(&REVERSE_EXON), model(&TWO_EXONS));
    assert_eq!(forward.cdna(1001), Some(1));
    assert_eq!(forward.cdna(1045), Some(45));
    assert_eq!(forward.cdna(1046), None);
    assert_eq!(reverse.cdna(1045), Some(1));
    assert_eq!(reverse.cdna(1001), Some(45));
    assert_eq!(forward.cds_pair(1007, 1007, false, 0), Some((4, 4)));
    assert_eq!(forward.cds_pair(1007, 1009, false, 0), Some((4, 6)));
    assert_eq!(forward.cds_pair(1007, 1009, false, 2), Some((6, 8)));
    assert_eq!(reverse.cds_pair(1037, 1039, true, 0), Some((4, 6)));
    assert_eq!(forward.cds_pair(1007, 1006, false, 0), Some((4, 3)));
    assert_eq!(reverse.cds_pair(1040, 1039, true, 0), Some((4, 3)));
    assert_eq!(forward.cds_pair(1004, 1003, false, 0), None);
    assert_eq!(forward.cds_pair(1004, 1004, false, 0), Some((1, 1)));
    assert_eq!(forward.cds_pair(1030, 1030, false, 0), Some((27, 27)));
    assert_eq!(forward.cds_pair(1031, 1031, false, 0), None);
    assert_eq!(forward.cds_pair(1031, 1030, false, 0), None);
    assert_eq!(forward.cds_pair(1030, 1029, false, 0), Some((27, 26)));
    assert_eq!(forward.cds_pair(1045, 1046, false, 0), None);
    assert_eq!(exons.cds_pair(1015, 1018, false, 0), Some((12, 13)));
    assert_eq!(exons.cds_pair(1015, 1016, false, 0), None);
    assert_eq!(exons.cds_pair(1016, 1015, false, 0), Some((13, 12)));
    assert_eq!(exons.cds_pair(1018, 1017, false, 0), Some((13, 12)));
    assert_eq!(exons.cds_pair(1017, 1016, false, 0), None);
    let one_base = model(&ONE_BASE_SEAM);
    assert_eq!(one_base.cds_pair(1016, 1015, false, 0), Some((13, 13)));
    let two_bases = |coding_start: i64, coding_end: i64| {
        spanning(&TWO_BASE_SEAM, coding_start, coding_end).cds_pair(1015, 1014, false, 0)
    };
    assert_eq!(two_bases(4, 30), Some((12, 13)));
    assert_eq!(two_bases(14, 30), Some((2, 3)));
    assert_eq!(two_bases(15, 30), None);
    assert_eq!(two_bases(4, 17), Some((12, 13)));
    assert_eq!(two_bases(4, 16), None);
}

#[test]
fn downstream_coding_models_reads_the_mapped_coding_transcripts_of_the_corpus_cache() {
    // Each transcript carrying a `cdna_coding_start`, a `cdna_coding_end` and mapper pairs
    // comes, under its chromosome and stable id, with them as numbers, the converted
    // cache's texts included, and the length of its `peptide`, `None` when the cache holds
    // none; a transcript lacking one of the three, or with a pair lacking a key or holding
    // a value that is no number, is left out.
    let pair = r#"{"from_start": 1, "from_end": 45, "to_start": 1001, "to_end": 1045, "ori": 1}"#;
    let text = r#"{"from_start": "1", "from_end": "45", "to_start": "1001", "to_end": "1045", "ori": "1"}"#;
    let signed = pair.replace(r#""ori": 1"#, r#""ori": "+""#);
    let span = r#""cdna_coding_start": 4, "cdna_coding_end": 30"#;
    let mapped = |pairs: &str| format!(r#""mapper": {{"pairs": [{pairs}]}}"#);
    let record = |stable_id: &str, span: &str, cache: &str| {
        format!(
            r#"{{"stable_id": "{stable_id}", {span}, "variation_effect_feature_cache": {{{cache}}}}}"#
        )
    };
    let records = [
        record(
            "ENST1",
            span,
            &format!(r#"{}, "peptide": "MKTLAVRG""#, mapped(pair)),
        ),
        record(
            "ENST2",
            r#""cdna_coding_start": "4", "cdna_coding_end": "30""#,
            &mapped(&format!("{text}, {text}")),
        ),
        record("ENST3", r#""cdna_coding_start": 4"#, &mapped(pair)),
        record("ENST4", span, r#""mapper": {"pair_count": 1}"#),
        record("ENST5", span, &mapped(&signed)),
        record("ENST6", span, &mapped(r#"{"from_start": 1}"#)),
        record(
            "ENST7",
            r#""cdna_coding_start": null, "cdna_coding_end": null"#,
            &mapped(pair),
        ),
        format!(r#"{{"stable_id": "ENST8", {span}}}"#),
    ];
    let (_scratch, cache_dir) = cache_of("models", &format!("[{}]", records.join(", ")));
    let model = |pairs: Vec<[i64; 5]>, peptide_length: Option<usize>| CodingModel {
        pairs,
        coding_start: 4,
        coding_end: 30,
        peptide_length,
    };
    let filed = |stable_id: &str| ("21".to_string(), stable_id.to_string());
    let expected: CodingModels = [
        (filed("ENST1"), model(ONE_EXON.to_vec(), Some(8))),
        (filed("ENST2"), model(vec![ONE_EXON[0], ONE_EXON[0]], None)),
    ]
    .into_iter()
    .collect();
    assert_eq!(downstream_coding_models(&cache_dir).unwrap(), expected);
}

#[test]
fn a_number_of_the_cache_is_read_as_int_reads_it() {
    // A coordinate is a number, a text `int()` reads (blanks around it, a sign, single
    // underscores between its digits), a double cut to its whole part or a boolean; a
    // transcript with any other value is left out. An empty text and an empty object in
    // the pairs' place hold no pair.
    let record = |stable_id: &str, start: &str, pairs: &str| {
        format!(
            r#"{{"stable_id": "{stable_id}", "cdna_coding_start": {start}, "cdna_coding_end": 30, "variation_effect_feature_cache": {{"mapper": {{"pairs": {pairs}}}}}}}"#
        )
    };
    let pair = r#"[{"from_start": " 1 ", "from_end": 45.9, "to_start": "+1_001", "to_end": 1045, "ori": true}]"#;
    let long = "1".repeat(4301);
    let records = [
        record("read", "4", pair),
        record("cut", "-4.7", "[]"),
        record("no_pairs_in_a_text", "4", r#""""#),
        record("no_pairs_in_an_object", "4", "{}"),
        record("fraction_in_a_text", r#""4.5""#, "[]"),
        record("no_number", "NaN", "[]"),
        record("list", "[4]", "[]"),
        record("too_long_for_int", &format!(r#""{long}""#), "[]"),
        record("characters_of_a_text", "4", r#""ab""#),
        record("keys_of_an_object", "4", r#"{"a": 1}"#),
        record("pairs_as_a_number", "4", "7"),
        record("pair_as_a_list", "4", "[[1, 45, 1001, 1045, 1]]"),
        record(
            "second_pair_unread",
            "4",
            r#"[{"from_start": 1, "from_end": 45, "to_start": 1001, "to_end": 1045, "ori": 1}, null]"#,
        ),
    ];
    let (_scratch, cache_dir) = cache_of("numbers", &format!("[{}]", records.join(", ")));
    let model = |pairs: Vec<[i64; 5]>, coding_start: i64| CodingModel {
        pairs,
        coding_start,
        coding_end: 30,
        peptide_length: None,
    };
    let filed = |stable_id: &str| ("21".to_string(), stable_id.to_string());
    let expected: CodingModels = [
        (filed("read"), model(ONE_EXON.to_vec(), 4)),
        (filed("cut"), model(Vec::new(), -4)),
        (filed("no_pairs_in_a_text"), model(Vec::new(), 4)),
        (filed("no_pairs_in_an_object"), model(Vec::new(), 4)),
    ]
    .into_iter()
    .collect();
    assert_eq!(downstream_coding_models(&cache_dir).unwrap(), expected);
}

#[test]
fn the_readers_refuse_a_cache_record_outside_what_they_hold() {
    // What CPython reads as a number the tool does not hold, or ends on a traceback for.
    let refusal = |records: &str, read: fn(&str) -> Option<Failure>| {
        let (_scratch, cache_dir) = cache_of("refused", records);
        let failure = read(&cache_dir).expect("a refusal");
        assert_eq!(failure.code, 1);
        failure.line
    };
    let sequences = |cache_dir: &str| downstream_sequences(cache_dir).err();
    let models = |cache_dir: &str| downstream_coding_models(cache_dir).err();
    assert_eq!(
        refusal(
            r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": [1]}]"#,
            sequences
        ),
        "ERROR: [vep-golden-corpus] the cache record of ENST1 holds a variation_effect_feature_cache that is no object"
    );
    assert_eq!(
        refusal(
            r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": "x"}]"#,
            models
        ),
        "ERROR: [vep-golden-corpus] the cache record of ENST1 holds a variation_effect_feature_cache that is no object"
    );
    let coding = format!(r#""translateable_seq": "{CODING}""#);
    assert_eq!(
        refusal(
            &format!(
                r#"[{{"stable_id": ["ENST1"], "variation_effect_feature_cache": {{{coding}}}}}]"#
            ),
            sequences
        ),
        "ERROR: [vep-golden-corpus] a cache record holds a stable_id that is a list or an object"
    );
    let table = |text: &str| {
        format!(
            r#"[{{"stable_id": "ENST1", "variation_effect_feature_cache": {{{coding}, "codon_table": "{text}"}}}}]"#
        )
    };
    assert_eq!(
        refusal(&table("\u{662}"), sequences),
        "ERROR: [vep-golden-corpus] the cache record of ENST1 holds a codon_table with a digit outside ASCII"
    );
    assert_eq!(
        refusal(&table(&"0".repeat(4301)), sequences),
        "ERROR: [vep-golden-corpus] the cache record of ENST1 holds a codon_table of more than 4,300 digits"
    );
    let start = |value: &str| {
        format!(
            r#"[{{"stable_id": "ENST1", "cdna_coding_start": {value}, "cdna_coding_end": 30, "variation_effect_feature_cache": {{"mapper": {{"pairs": []}}}}}}]"#
        )
    };
    let outside = "ERROR: [vep-golden-corpus] the cache record of ENST1 holds a coordinate that is no number of 64 bits in ASCII digits";
    for value in [
        "9223372036854775808",
        r#""9223372036854775808""#,
        "1e19",
        "Infinity",
        "\"\u{664}\"",
    ] {
        assert_eq!(refusal(&start(value), models), outside, "{value}");
    }
    // A text of the same digits under 64 bits is a number, and the record of a transcript
    // whose id is no text is left out, a list among them.
    let (_scratch, cache_dir) = cache_of("held", &start("\"9223372036854775807\""));
    let held = downstream_coding_models(&cache_dir).unwrap();
    assert_eq!(
        held[&("21".to_string(), "ENST1".to_string())].coding_start,
        i64::MAX
    );
    let unnamed = start("4").replace(r#""ENST1""#, "[1]");
    let (_scratch, cache_dir) = cache_of("unnamed", &unnamed);
    assert!(downstream_coding_models(&cache_dir).unwrap().is_empty());
}
