// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the mate-side derivation: every expected set is derived by hand from the
//! transcript structure quoted in the test, and a cache record is read as `int()` reads its
//! numbers.

use super::*;
use crate::json::loads;
use crate::TempDir;

fn exons(spec: &str) -> Vec<(i128, i128)> {
    spec.split(';')
        .map(|exon| {
            let (start, end) = exon.split_once('-').unwrap();
            (start.parse().unwrap(), end.parse().unwrap())
        })
        .collect()
}

fn model(
    chrom: &str,
    span: (i128, i128),
    strand: i128,
    biotype: &str,
    exon_spec: &str,
    cds: Option<(i128, i128)>,
) -> TranscriptModel {
    TranscriptModel {
        chrom: chrom.to_string(),
        start: span.0,
        end: span.1,
        strand,
        biotype: biotype.to_string(),
        exons: exons(exon_spec),
        cds,
        mirna_mature: Vec::new(),
    }
}

/// ENST00000270142 (chromosome 21, forward strand, protein coding): five exons, CDS
/// 33032083-33040891.
fn app_like() -> TranscriptModel {
    model(
        "21",
        (33031935, 33041244),
        1,
        "protein_coding",
        "33031935-33032154;33036103-33036199;33038762-33038831;33039571-33039688;33040784-33041244",
        Some((33032083, 33040891)),
    )
}

/// ENST00000248935 GSTT1 (chromosome 22, reverse strand, protein coding), CDS
/// 24376423-24384231.
fn gstt1() -> TranscriptModel {
    model(
        "22",
        (24376133, 24384284),
        -1,
        "protein_coding",
        "24376133-24376617;24376822-24376998;24379361-24379511;24381700-24381787;24384120-24384284",
        Some((24376423, 24384231)),
    )
}

/// ENST00000265107 WDR70 (chromosome 5, forward strand, protein coding), CDS
/// 37379368-37752573.
fn wdr70() -> TranscriptModel {
    model(
        "5",
        (37379318, 37753435),
        1,
        "protein_coding",
        concat!(
            "37379318-37379392;37379489-37379554;37381602-37381685;37392000-37392120;",
            "37396375-37396570;37437922-37437981;37443239-37443372;37479834-37479987;",
            "37516514-37516590;37605064-37605238;37697655-37697754;37701058-37701142;",
            "37702949-37703087;37721115-37721215;37722855-37722934;37724934-37725050;",
            "37726883-37727045;37752486-37753435"
        ),
        Some((37379368, 37752573)),
    )
}

/// ENST00000385128 MIR125B2 (chromosome 21, forward strand, miRNA), one exon, mature 17-38
/// and 54-75.
fn mir() -> TranscriptModel {
    TranscriptModel {
        mirna_mature: vec![(17, 38), (54, 75)],
        ..model(
            "21",
            (16590237, 16590325),
            1,
            "miRNA",
            "16590237-16590325",
            None,
        )
    }
}

/// ENST00000687248 PTCHD1-AS (chromosome X, reverse strand, lncRNA), nine exons, no CDS.
fn ptchd1_as() -> TranscriptModel {
    model(
        "X",
        (22227405, 23293174),
        -1,
        "lncRNA",
        concat!(
            "22227405-22227667;22253939-22254065;22258013-22258136;22259755-22259842;",
            "22389663-22389720;22743229-22743338;23064038-23064121;23270051-23270176;",
            "23293014-23293174"
        ),
        None,
    )
}

fn object(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("no JSON object: {other:?}"),
    }
}

/// A cache record of two exons, a coding sequence from cDNA 51 to 202 and one mature miRNA
/// range, with the member `name` replaced by the JSON text `value`.
fn record_with(name: &str, value: &str) -> Object {
    let mut record = object(concat!(
        r#"{"stable_id": "T1", "start": 100, "end": 1000, "strand": 1, "#,
        r#""biotype": "protein_coding", "#,
        r#""exons": [{"start": 100, "end": 200}, {"start": 800, "end": 1000}], "#,
        r#""variation_effect_feature_cache": {"mapper": {"cdna_coding_start": 51, "#,
        r#""cdna_coding_end": 202, "pairs": ["#,
        r#"{"from_start": 1, "from_end": 101, "to_start": 100, "to_end": 200, "ori": 1}, "#,
        r#"{"from_start": 102, "from_end": 302, "to_start": 800, "to_end": 1000, "ori": 1}]}}, "#,
        r#""attributes": [{"code": "miRNA", "value": "3-7"}]}"#
    ));
    let value = loads(value).unwrap_or_else(|error| panic!("{value}: {error}"));
    record.insert(name, value);
    record
}

/// The model of `record_with(name, value)`, filed under chromosome 21.
fn model_with(name: &str, value: &str) -> Result<Option<TranscriptModel>, Failure> {
    model_from_cache_record("21", &record_with(name, value))
}

/// The record as it is written.
fn as_written() -> TranscriptModel {
    TranscriptModel {
        mirna_mature: vec![(3, 7)],
        ..model(
            "21",
            (100, 1000),
            1,
            "protein_coding",
            "100-200;800-1000",
            Some((150, 900)),
        )
    }
}

fn terms(names: &[&'static str]) -> Terms {
    names.iter().copied().collect()
}

#[test]
fn introns_are_the_gaps_between_sorted_exons() {
    let introns = app_like().introns();
    assert_eq!(introns[0], (33032155, 33036102));
    assert_eq!(introns.len(), 4);
    assert_eq!(mir().introns(), []);
}

#[test]
fn genomic_cds_is_mapped_through_the_pairs_in_both_orientations() {
    // Forward: cDNA 149 in a first exon of 220 bases starting at 33031935 is 33032083.
    // Reverse: cDNA 1 of a reverse exon 24384120-24384284 is 24384284, so a coding start at
    // cDNA 54 is 24384231.
    let forward = object(concat!(
        r#"{"variation_effect_feature_cache": {"mapper": {"cdna_coding_start": 149, "#,
        r#""cdna_coding_end": 613, "pairs": ["#,
        r#"{"from_start": 1, "from_end": 220, "to_start": 33031935, "to_end": 33032154, "ori": 1}, "#,
        r#"{"from_start": 221, "from_end": 317, "to_start": 33036103, "to_end": 33036199, "ori": 1}, "#,
        r#"{"from_start": 318, "from_end": 387, "to_start": 33038762, "to_end": 33038831, "ori": 1}, "#,
        r#"{"from_start": 388, "from_end": 505, "to_start": 33039571, "to_end": 33039688, "ori": 1}, "#,
        r#"{"from_start": 506, "from_end": 966, "to_start": 33040784, "to_end": 33041244, "ori": 1}"#,
        "]}}}"
    ));
    assert_eq!(genomic_cds(&forward), Ok(Some((33032083, 33040891))));
    let reverse = object(concat!(
        r#"{"variation_effect_feature_cache": {"mapper": {"cdna_coding_start": 54, "#,
        r#""cdna_coding_end": 195, "pairs": ["#,
        r#"{"from_start": 1, "from_end": 165, "to_start": 24384120, "to_end": 24384284, "ori": -1}, "#,
        r#"{"from_start": 166, "from_end": 253, "to_start": 24381700, "to_end": 24381787, "ori": -1}"#,
        "]}}}"
    ));
    assert_eq!(genomic_cds(&reverse), Ok(Some((24381758, 24384231))));
}

#[test]
fn a_transcript_without_coding_bounds_has_no_cds() {
    let unbounded = object(r#"{"variation_effect_feature_cache": {"mapper": {"pairs": []}}}"#);
    assert_eq!(genomic_cds(&unbounded), Ok(None));
    assert_eq!(genomic_cds(&Object::new()), Ok(None));
}

fn cache(tag: &str, shards: &[(&str, String)]) -> TempDir {
    let scratch = TempDir::create(&std::env::temp_dir(), &format!("golden_mate_{tag}_")).unwrap();
    for (chrom, shard) in shards {
        let directory = scratch.path.join("transcripts").join(chrom);
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(directory.join("1-1000000.json"), shard).unwrap();
    }
    scratch
}

#[test]
fn models_are_keyed_by_chromosome_and_id() {
    // A pseudoautosomal transcript is filed under X and Y with different coordinates; one
    // key a chromosome keeps both.
    let shard = |start: u64| {
        format!(
            concat!(
                r#"[{{"stable_id": "ENST_PAR", "start": {}, "end": {}, "strand": 1, "#,
                r#""biotype": "protein_coding", "exons": [{{"start": {}, "end": {}}}], "#,
                r#""variation_effect_feature_cache": {{"mapper": {{"pairs": [], "pair_count": 0}}}}}}]"#
            ),
            start,
            start + 999,
            start,
            start + 999
        )
    };
    let scratch = cache("keyed", &[("X", shard(200_000)), ("Y", shard(150_000))]);
    let models = load_transcript_models(scratch.path.to_str().unwrap()).unwrap();
    assert_eq!(models.len(), 2);
    assert_eq!(models.get("X", "ENST_PAR").unwrap().start, 200_000);
    assert_eq!(models.get("Y", "ENST_PAR").unwrap().start, 150_000);
    assert_eq!(models.get("21", "ENST_PAR"), None);
    let absent = scratch.path.join("absent");
    let none = load_transcript_models(absent.to_str().unwrap()).unwrap();
    assert!(none.is_empty());
}

#[test]
fn a_storable_shaped_shard_is_read() {
    let shard = concat!(
        r#"{"21": [null, {"stable_id": "ENST_A", "start": 10, "end": 20, "strand": -1, "#,
        r#""biotype": "lncRNA", "exons": [{"start": 10, "end": 20}], "#,
        r#""attributes": [{"code": "miRNA", "value": "3-7"}]}]}"#
    );
    let scratch = cache("storable", &[("21", shard.to_string())]);
    let models = load_transcript_models(scratch.path.to_str().unwrap()).unwrap();
    let model = models.get("21", "ENST_A").unwrap();
    assert_eq!(model.mirna_mature, [(3, 7)]);
    assert_eq!(model.cds, None);
}

#[test]
fn a_later_record_of_a_key_stands_and_a_number_names_no_feature() {
    let shard = concat!(
        r#"[{"stable_id": "T1", "start": 10, "end": 20, "strand": 1}, "#,
        r#"{"stable_id": 5, "start": 30, "end": 40, "strand": 1}, "#,
        r#"{"stable_id": true, "start": 30, "end": 40, "strand": 1}, "#,
        r#"{"stable_id": "T2", "start": null, "end": 40, "strand": 1}, "#,
        r#"{"stable_id": ["T3"], "start": null, "end": 40, "strand": 1}, "#,
        r#"{"stable_id": "T1", "start": 50, "end": 60, "strand": -1}]"#
    );
    let scratch = cache("later", &[("21", shard.to_string())]);
    let models = load_transcript_models(scratch.path.to_str().unwrap()).unwrap();
    assert_eq!(models.len(), 1);
    let model = models.get("21", "T1").unwrap();
    assert_eq!((model.start, model.end, model.strand), (50, 60, -1));
    assert_eq!(model.biotype, "");
    assert_eq!(model.exons, []);

    for (id, name) in [(r#"["T1"]"#, "['T1']"), (r#"{"a": 1}"#, "{'a': 1}")] {
        let shard = format!(r#"[{{"stable_id": {id}, "start": 10, "end": 20, "strand": 1}}]"#);
        let scratch = cache("named", &[("21", shard)]);
        let refusal = load_transcript_models(scratch.path.to_str().unwrap()).unwrap_err();
        let line = format!(
            "ERROR: [vep-golden-corpus] the transcript {name} filed under 21 in the corpus cache \
             is named by a list or an object"
        );
        assert_eq!((refusal.line, refusal.code), (line, 1));
    }
}

#[test]
fn a_number_of_a_record_is_read_as_int_reads_it() {
    assert_eq!(model_with("start", "100"), Ok(Some(as_written())));
    // A text between the blanks `int()` drops, with underscores between its digits, and in
    // the decimal digits of another script; a float cut at its point; `true` as 1.
    let digits = "\u{661}\u{660}\u{660}";
    for text in ["100", " 100\\n", "1_0_0", "+100", digits] {
        let value = format!("\"{text}\"");
        assert_eq!(
            model_with("start", &value),
            Ok(Some(as_written())),
            "{text}"
        );
    }
    assert_eq!(model_with("start", "100.9"), Ok(Some(as_written())));
    let below = TranscriptModel {
        start: -100,
        ..as_written()
    };
    assert_eq!(model_with("start", "-100.9"), Ok(Some(below)));
    let one = TranscriptModel {
        start: 1,
        ..as_written()
    };
    assert_eq!(model_with("start", "true"), Ok(Some(one)));
    let last = TranscriptModel {
        start: i128::from(i64::MAX),
        ..as_written()
    };
    assert_eq!(model_with("start", "9223372036854775807"), Ok(Some(last)));
    // No integer: no model, whatever else the record holds.
    let long = format!("\"{}\"", "1".repeat(4301));
    for value in [
        "null",
        "[100]",
        "{}",
        "NaN",
        "\"\"",
        "\"1__0\"",
        "\"_10\"",
        "\"10_\"",
        "\"1 0\"",
        "\"0x10\"",
        "\"\\u001c100\"",
        long.as_str(),
    ] {
        assert_eq!(model_with("start", value), Ok(None), "{value}");
        assert_eq!(model_with("strand", value), Ok(None), "{value}");
    }
    let mut unnamed = record_with("start", "null");
    unnamed.insert("exons", "100-200");
    assert_eq!(model_from_cache_record("21", &unnamed), Ok(None));
}

#[test]
fn the_members_of_a_record_are_read_as_the_cache_writes_them() {
    let without_exons = TranscriptModel {
        exons: Vec::new(),
        ..as_written()
    };
    for value in ["null", "0", "{}", "[]", "false", "\"\""] {
        assert_eq!(
            model_with("exons", value),
            Ok(Some(without_exons.clone())),
            "{value}"
        );
    }
    let unsorted = r#"[{"start": 800, "end": "1000"}, {"start": 100, "end": 200}, {"start": 5}]"#;
    assert_eq!(model_with("exons", unsorted), Ok(Some(as_written())));

    let without_cds = TranscriptModel {
        cds: None,
        ..as_written()
    };
    let mapper = |members: &str| format!(r#"{{"mapper": {{{members}}}}}"#);
    let pair = r#"{"from_start": 1, "from_end": 302, "to_start": 100, "to_end": 401, "ori": 1}"#;
    for value in [
        "null".to_string(),
        "{}".to_string(),
        r#"{"mapper": null}"#.to_string(),
        mapper(r#""pairs": "no pairs are read without coding bounds""#),
        mapper(&format!(
            r#""cdna_coding_start": 0, "cdna_coding_end": 202, "pairs": [{pair}]"#
        )),
        mapper(&format!(
            r#""cdna_coding_start": -5, "cdna_coding_end": 202, "pairs": [{pair}]"#
        )),
        mapper(&format!(
            r#""cdna_coding_start": 51, "cdna_coding_end": 5000, "pairs": [{pair}]"#
        )),
        mapper(r#""cdna_coding_start": 51, "cdna_coding_end": 202, "pairs": null"#),
        mapper(&format!(
            r#""cdna_coding_start": 51, "cdna_coding_end": 202, "pairs": [{}]"#,
            pair.replace(r#""ori": 1"#, r#""ori": null"#)
        )),
    ] {
        assert_eq!(
            model_with("variation_effect_feature_cache", &value),
            Ok(Some(without_cds.clone())),
            "{value}"
        );
    }
    // The reverse orientation counts down from the pair's end, and the bounds come back low
    // first.
    let reverse = mapper(&format!(
        r#""cdna_coding_start": "51", "cdna_coding_end": 202, "pairs": [{}]"#,
        pair.replace(r#""ori": 1"#, r#""ori": -1"#)
    ));
    let reversed = TranscriptModel {
        cds: Some((200, 351)),
        ..as_written()
    };
    assert_eq!(
        model_with("variation_effect_feature_cache", &reverse),
        Ok(Some(reversed))
    );

    let without_mature = TranscriptModel {
        mirna_mature: Vec::new(),
        ..as_written()
    };
    let attribute = |members: &str| format!("[{{{members}}}]");
    let valued = |value: &str| attribute(&format!("\"code\": \"miRNA\", \"value\": {value}"));
    for (value, ranges) in [
        (valued("\"\\u001c 3-7 \""), 1),
        (valued("\"\u{663}-\u{667}\""), 1),
        (valued("\"3 - 7\""), 0),
        (valued("\"3-7-9\""), 0),
        (valued("\"3-7\u{b2}\""), 0),
        (valued("37"), 0),
        (valued("null"), 0),
        (attribute("\"code\": \"miRNA\""), 0),
        (attribute("\"code\": 5, \"value\": \"3-7\""), 0),
        (attribute("\"code\": \"other\", \"value\": \"3-7\""), 0),
    ] {
        let expected = if ranges == 1 {
            as_written()
        } else {
            without_mature.clone()
        };
        assert_eq!(
            model_with("attributes", &value),
            Ok(Some(expected)),
            "{value}"
        );
    }
    let long = valued(&format!("\"3-{}\"", "7".repeat(4300)));
    let read = model_with("attributes", &long).unwrap().unwrap();
    assert_eq!(read.mirna_mature, [(3, i128::MAX)]);

    for (value, biotype) in [
        ("5", "5"),
        (r#"["a", 1]"#, "['a', 1]"),
        ("null", ""),
        ("0", ""),
    ] {
        let expected = TranscriptModel {
            biotype: biotype.to_string(),
            ..as_written()
        };
        assert_eq!(model_with("biotype", value), Ok(Some(expected)), "{value}");
    }
}

#[test]
fn a_record_whose_read_ends_in_an_exception_is_refused() {
    let refused = |name: &str, value: &str, detail: &str| {
        let line = format!(
            "ERROR: [vep-golden-corpus] the transcript 'T1' filed under 21 in the corpus cache {detail}"
        );
        let refusal = model_with(name, value).unwrap_err();
        assert_eq!((refusal.line, refusal.code), (line, 1), "{name} {value}");
    };
    let infinite = "holds a number that is infinite";
    refused("start", "Infinity", infinite);
    refused("strand", "-Infinity", infinite);
    refused("exons", r#"[{"start": 100, "end": Infinity}]"#, infinite);
    let past = "holds a number that 64 bits do not hold";
    refused("start", "9223372036854775808", past);
    refused("end", "-9223372036854775809", past);
    refused("start", "1e30", past);
    refused("start", r#""1180591620717411303424""#, past);
    refused(
        "start",
        r#""1361129467683753853853498429727072845824""#,
        past,
    );
    let exons = "holds exons that are no list of objects";
    for value in [
        r#""100-200""#,
        r#"{"start": 100, "end": 200}"#,
        "5",
        "true",
        r#"[{"start": 100, "end": 200}, null]"#,
        "[[100, 200]]",
    ] {
        refused("exons", value, exons);
    }
    let attributes = "holds attributes that are no list of objects";
    for value in [
        r#""miRNA""#,
        r#"["miRNA"]"#,
        r#"{"code": "miRNA", "value": "3-7"}"#,
    ] {
        refused("attributes", value, attributes);
    }
    let cache = "holds a variation_effect_feature_cache that is no object";
    for value in ["[1]", r#""x""#, "5"] {
        refused("variation_effect_feature_cache", value, cache);
    }
    let mapper = |members: &str| format!(r#"{{"mapper": {members}}}"#);
    for value in ["[1]", r#""x""#] {
        refused(
            "variation_effect_feature_cache",
            &mapper(value),
            "holds a mapper that is no object",
        );
    }
    let bounds = r#""cdna_coding_start": 51, "cdna_coding_end": 202"#;
    for value in [r#""x""#, r#"{"from_start": 1}"#, "[null]"] {
        refused(
            "variation_effect_feature_cache",
            &mapper(&format!(r#"{{{bounds}, "pairs": {value}}}"#)),
            "holds pairs that are no list of objects",
        );
    }
    refused(
        "variation_effect_feature_cache",
        &mapper(r#"{"cdna_coding_start": Infinity, "cdna_coding_end": 202}"#),
        infinite,
    );
    refused(
        "variation_effect_feature_cache",
        &mapper(&format!(r#"{{{bounds}, "pairs": [{{"ori": Infinity}}]}}"#)),
        infinite,
    );
    let long = format!(
        r#"[{{"code": "miRNA", "value": "3-{}"}}]"#,
        "7".repeat(4301)
    );
    refused(
        "attributes",
        &long,
        "holds a miRNA range of more than 4,300 digits",
    );
}

#[test]
fn cds_utr_by_strand() {
    let coding = terms(&["coding_sequence_variant"]);
    assert_eq!(region_terms(33032100, &app_like()), coding);
    let five = terms(&["5_prime_UTR_variant"]);
    let three = terms(&["3_prime_UTR_variant"]);
    assert_eq!(region_terms(33032000, &app_like()), five);
    assert_eq!(region_terms(33041000, &app_like()), three);
    // Reverse strand: bases beyond the CDS end are the 5' UTR.
    assert_eq!(region_terms(24384250, &gstt1()), five);
    assert_eq!(region_terms(24376300, &gstt1()), three);
}

#[test]
fn intron_terms() {
    // Inside the first intron of the forward transcript (33032155-33036102): the first two
    // and last two bases are the splice sites and carry no intron term; the forward-strand
    // polypyrimidine window is the acceptor side, 16 to 2 bases before the intron end.
    let app = app_like();
    let intron = terms(&["intron_variant"]);
    let tract = terms(&["intron_variant", "splice_polypyrimidine_tract_variant"]);
    assert_eq!(region_terms(33034355, &app), intron);
    assert_eq!(region_terms(33032155, &app), Terms::new());
    assert_eq!(region_terms(33036102, &app), Terms::new());
    assert_eq!(region_terms(33036102 - 5, &app), tract);
    assert_eq!(region_terms(33036102 - 17, &app), intron);
    // Reverse strand: the acceptor is the intron's start side.
    assert_eq!(region_terms(24376618 + 5, &gstt1()), tract);
}

#[test]
fn non_coding_intron_and_exon() {
    let lncrna = ptchd1_as();
    assert_eq!(
        region_terms(22500000, &lncrna),
        terms(&["intron_variant", "non_coding_transcript_variant"])
    );
    assert_eq!(
        region_terms(22227500, &lncrna),
        terms(&["non_coding_transcript_exon_variant"])
    );
}

#[test]
fn nmd_transcript_adds_its_term() {
    let nmd = model(
        "21",
        (100, 1000),
        1,
        "nonsense_mediated_decay",
        "100-200;800-1000",
        Some((150, 900)),
    );
    assert_eq!(
        region_terms(500, &nmd),
        terms(&["intron_variant", "NMD_transcript_variant"])
    );
    assert_eq!(
        region_terms(160, &nmd),
        terms(&["coding_sequence_variant", "NMD_transcript_variant"])
    );
}

#[test]
fn mature_mirna_hit_stands_alone() {
    // cDNA 17-38 of a forward single-exon miRNA starting at 16590237 is 16590253-16590274.
    assert_eq!(
        region_terms(16590260, &mir()),
        terms(&["mature_miRNA_variant"])
    );
    assert_eq!(
        region_terms(16590240, &mir()),
        terms(&["non_coding_transcript_exon_variant"])
    );
}

#[test]
fn inside_with_the_mate_inside_carries_feature_truncation() {
    assert_eq!(
        breakend_context(33034355, &app_like(), true),
        terms(&["feature_truncation", "intron_variant"])
    );
}

#[test]
fn upstream_downstream_are_strand_aware() {
    // 4,108 bases before the start of the reverse-strand transcript: downstream.
    let downstream = terms(&["downstream_gene_variant"]);
    let upstream = terms(&["upstream_gene_variant"]);
    assert_eq!(breakend_context(24372025, &gstt1(), false), downstream);
    assert_eq!(breakend_context(24384284 + 100, &gstt1(), false), upstream);
    assert_eq!(
        breakend_context(33031935 - 100, &app_like(), false),
        upstream
    );
}

#[test]
fn beyond_five_kb_carries_only_the_truncation_term_when_the_mate_is_inside() {
    assert_eq!(
        breakend_context(14374334, &app_like(), true),
        terms(&["feature_truncation"])
    );
    assert_eq!(breakend_context(14374334, &app_like(), false), Terms::new());
}

#[test]
fn a_mature_mirna_hit_takes_no_truncation_term() {
    assert_eq!(
        breakend_context(16590260, &mir(), true),
        terms(&["mature_miRNA_variant"])
    );
}

#[test]
fn cross_chromosome_row_read_at_the_local_coordinate() {
    // `1:14374334 ]21:33034355]A` on the forward transcript: the mate sits in the first
    // intron, the local coordinate on chromosome 1 is nowhere near, so Ensembl VEP writes
    // only `feature_truncation` and the mate supports `feature_truncation,intron_variant`.
    let derived = mate_side_derivation("1:14374334", "]21:33034355]A", &app_like());
    let mate = terms(&["feature_truncation", "intron_variant"]);
    assert_eq!(derived, Some((mate, terms(&["feature_truncation"]))));
}

#[test]
fn the_local_coordinate_landing_in_a_cds_exon_of_the_mate_transcript() {
    // `21:37392051 [5:37717217[A` on WDR70: the mate is in intron 13 (37703088-37721114),
    // the local coordinate lands numerically in exon 4 (37392000-37392120), inside the CDS:
    // Ensembl VEP's read is `feature_truncation,coding_sequence_variant`.
    let derived = mate_side_derivation("21:37392051", "[5:37717217[A", &wdr70());
    let mate = terms(&["feature_truncation", "intron_variant"]);
    let local = terms(&["feature_truncation", "coding_sequence_variant"]);
    assert_eq!(derived, Some((mate, local)));
}

#[test]
fn a_far_local_coordinate_leaves_only_the_truncation_term() {
    let derived = mate_side_derivation("21:5100001", "[5:37717217[T", &wdr70());
    assert_eq!(derived.unwrap().1, terms(&["feature_truncation"]));
}

#[test]
fn a_far_local_coordinate_beside_an_outside_mate_defaults_to_intergenic() {
    // `21:9500001 N[22:24372025[` on GSTT1: the mate is 4,108 bases downstream, the local
    // coordinate fires nothing, and Ensembl VEP writes `intergenic_variant` on the
    // Transcript row.
    let derived = mate_side_derivation("21:9500001", "N[22:24372025[", &gstt1());
    let mate = terms(&["downstream_gene_variant"]);
    assert_eq!(derived, Some((mate, terms(&["intergenic_variant"]))));
}

#[test]
fn equal_numeric_coordinates_give_equal_sets() {
    // `21:24372025 N[22:24372025[`: the two reads coincide, which is why Ensembl VEP's row
    // is right on this record.
    let (mate, local) = mate_side_derivation("21:24372025", "N[22:24372025[", &gstt1()).unwrap();
    assert_eq!(mate, local);
    assert_eq!(mate, terms(&["downstream_gene_variant"]));
}

#[test]
fn a_row_the_local_coordinate_selects_is_not_a_mate_side_row() {
    let derived = mate_side_derivation("21:33034000", "N[21:33034355[", &app_like());
    assert_eq!(derived, None);
}

#[test]
fn a_ranged_location_or_a_non_bracket_allele_is_not_a_mate_side_row() {
    let ranged = mate_side_derivation("21:14604132-41770788", "N[21:41770788[", &app_like());
    assert_eq!(ranged, None);
    assert_eq!(
        mate_side_derivation("1:14374334", "deletion", &app_like()),
        None
    );
}

#[test]
fn a_mate_beyond_five_kb_of_the_transcript_is_not_a_mate_side_row() {
    let derived = mate_side_derivation("1:14374334", "]21:33050000]A", &app_like());
    assert_eq!(derived, None);
}

#[test]
fn a_transcript_on_another_chromosome_than_the_mate_is_not_a_mate_side_row() {
    let derived = mate_side_derivation("1:14374334", "]22:33034355]A", &app_like());
    assert_eq!(derived, None);
}

#[test]
fn helpers() {
    assert_eq!(bracket_mate("N[chr21:41770788["), Some(("21", 41770788)));
    assert_eq!(bracket_mate("]13:123456]T"), Some(("13", 123456)));
    assert_eq!(bracket_mate("N."), None);
    assert_eq!(point_location("chr21:5"), Some(("21", 5)));
    assert_eq!(point_location("21:1-2"), None);
}

#[test]
fn a_position_is_read_in_the_decimal_digits_of_any_script() {
    // The first bracketed chromosome and position of an allele are its mate; a `chr` prefix
    // goes in any case, and with it nothing else.
    assert_eq!(bracket_mate("N[CHR21:5[N]1:7]"), Some(("21", 5)));
    assert_eq!(bracket_mate("N[21:5"), None);
    assert_eq!(
        bracket_mate("N[chromosome21:0012["),
        Some(("omosome21", 12))
    );
    let fullwidth = "N[22:\u{ff11}\u{ff12}[";
    assert_eq!(bracket_mate(fullwidth), Some(("22", 12)));
    let arabic_indic = "\u{663}\u{660}";
    assert_eq!(
        point_location(&format!("21:{arabic_indic}")),
        Some(("21", 30))
    );
    for (digit, value) in [
        ('7', 7),
        ('\u{660}', 0),
        ('\u{669}', 9),
        ('\u{1d7d8}', 0),
        ('\u{1d7e1}', 9),
    ] {
        assert_eq!(decimal_digit(digit), Some(value), "{digit:?}");
    }
    // A superscript two is a digit and no decimal digit: `int()` does not read it.
    for other in ['\u{b2}', 'x', '-', ' '] {
        assert_eq!(decimal_digit(other), None, "{other:?}");
    }
    for location in ["21", "21:", "21:x", "21:5 ", "21:\u{b2}", "21:-5"] {
        assert_eq!(point_location(location), None, "{location}");
    }
    // A position past 128 bits is the largest position there is, and one of more digits
    // than `int()` converts is none.
    let long = "9".repeat(4300);
    assert_eq!(
        point_location(&format!("21:{long}")),
        Some(("21", i128::MAX))
    );
    assert_eq!(point_location(&format!("21:{long}9")), None);
    assert_eq!(
        bracket_mate(&format!("N[21:{long}[")),
        Some(("21", i128::MAX))
    );
    assert_eq!(bracket_mate(&format!("N[21:{long}9[")), None);
    assert_eq!(
        mate_side_derivation(&format!("1:{long}"), "]21:33034355]A", &app_like()),
        Some((
            terms(&["feature_truncation", "intron_variant"]),
            terms(&["feature_truncation"])
        ))
    );
}
