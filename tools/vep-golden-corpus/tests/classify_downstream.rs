// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field` on made-up outputs whose rows the Downstream plugin annotated,
//! and on rows of transcripts under an organelle's codon table: a `DownstreamProtein` or
//! `ProteinLengthChange` that differs takes one of the three Downstream classes where the
//! sequences of the corpus cache give exactly the two values, an `HGVSp` that differs on
//! an organelle transcript takes the organelle class on four peptide paths, and every
//! other difference stays unexplained.

mod classify_support;
mod common;

use classify_support::{
    counted, field_divergence, fields_written, key, nothing_written, obj, record, records,
    transcript, with_member, Case,
};
use common::clean;
use vep_golden_corpus::classify::downstream::{translate_codons, MAPPER_PAIR_KEYS};
use vep_golden_corpus::classify::organelle::{
    organelle_codon_tables, organelle_peptide_path_divergence, peptide_description_shape,
};
use vep_golden_corpus::json::{dumps, Object, Style, Value};

/// The transcript the Downstream plugin is tested on in
/// crates/vep-builtins/src/plugins/downstream.rs: eight codons and a stop, and a 3' UTR
/// that reads three codons and a stop from its first base.
const CODING: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";
const UTR: &str = "GGGCCCAAATAGCCC";
/// A second transcript, of five codons and a stop, whose second codon becomes a stop when
/// its third base goes.
const STOP_CODING: &str = "ATGTACGGAAAACTGTAA";
/// A third transcript, of four codons and a stop, filed under codon table 2 with a 3' UTR
/// of three bases: when the third base of its second codon goes, the read crosses a codon
/// that is tryptophan under that table and a stop under the standard one, and ends at a
/// codon that is a stop under table 2 alone.
const TABLE_TWO_CODING: &str = "ATGAAACTGACCTAG";
const TABLE_TWO_UTR: &str = "ACC";

const START_PAST_OWN_CODON: &str = "downstream_start_past_own_codon";
const UTR_FRAME_PAST_CODING_END: &str = "downstream_utr_frame_past_coding_end";
const HGVS_SHIFTED_CDS_COORDINATES: &str = "plugin_hgvs_shifted_cds_coordinates";
const ORGANELLE: &str = "organelle_peptide_paths_translation_table";
const UNEXPLAINED: &str = "unexplained_residual";

const PROTEIN: &str = "DownstreamProtein";
const LENGTH_CHANGE: &str = "ProteinLengthChange";
/// The fields the Downstream runs compare.
const DOWNSTREAM_RUN: [&str; 3] = [PROTEIN, LENGTH_CHANGE, "HGVS_OFFSET"];

/// Mapper pairs (first and last cDNA base, first and last genomic base, orientation) that
/// lay a transcript of three 5' UTR bases, a coding sequence and the UTR on the genome
/// from base 1001. On one exon coding base N is genomic base 1003 + N; on the reverse
/// strand the first transcript's coding base N is genomic base 1043 - N.
const ONE_EXON: &[[i64; 5]] = &[[1, 45, 1001, 1045, 1]];
const REVERSE_EXON: &[[i64; 5]] = &[[1, 45, 1001, 1045, -1]];
/// Two exons around a two-base intron after coding base 12: coding base 13 is genomic
/// base 1018.
const TWO_EXONS: &[[i64; 5]] = &[[1, 15, 1001, 1015, 1], [16, 45, 1018, 1047, 1]];
/// Two pairs that touch on the genome between bases 1015 and 1016 while the model carries
/// coding base 13 between them, and two that touch between bases 1014 and 1015 around
/// coding bases 12 and 13.
const ONE_BASE_SEAM: &[[i64; 5]] = &[[1, 15, 1001, 1015, 1], [17, 45, 1016, 1044, 1]];
const TWO_BASE_SEAM: &[[i64; 5]] = &[[1, 14, 1001, 1014, 1], [17, 45, 1015, 1043, 1]];
/// The second transcript on one exon: its 18 coding bases end at cDNA base 21.
const STOP_EXON: &[[i64; 5]] = &[[1, 36, 1001, 1036, 1]];

/// One row of a run as both engines write it: the name of its record, its key, its
/// `CDS_position` column, and the Extra column of Ensembl VEP and of vep-rs.
struct Row {
    name: String,
    location: String,
    allele: String,
    feature: String,
    terms: String,
    cds: String,
    extras: [String; 2],
}

/// The class a differing field of a record's row is expected under, with the class of a
/// second cause where the entry names one.
struct Classed {
    name: &'static str,
    field: &'static str,
    class: &'static str,
    second: Option<&'static str>,
}

/// The class of a record's `DownstreamProtein`.
fn protein(name: &'static str, class: &'static str) -> Classed {
    Classed {
        name,
        field: PROTEIN,
        class,
        second: None,
    }
}

/// The class of a record's `ProteinLengthChange`.
fn length_change(name: &'static str, class: &'static str) -> Classed {
    Classed {
        field: LENGTH_CHANGE,
        ..protein(name, class)
    }
}

/// The value of a field of an Extra column, its last where the column gives it twice, `-`
/// for a field it lacks.
fn extra_value<'a>(extra: &'a str, field: &str) -> &'a str {
    let mut items = extra.rsplit(';').filter_map(|item| item.split_once('='));
    items
        .find(|(name, _)| *name == field)
        .map_or("-", |(_, value)| value)
}

/// A row of a default-format output with a `CDS_position`, on a transcript of the gene
/// `ENSG1`.
fn coding_row(row: &Row, extra: &str) -> String {
    let columns = [
        row.name.as_str(),
        &row.location,
        &row.allele,
        "ENSG1",
        &row.feature,
        "Transcript",
        &row.terms,
        "-",
        &row.cds,
        "-",
        "-",
        "-",
        "-",
        extra,
    ];
    columns.join("\t") + "\n"
}

/// What a run is held to: its line of standard output and its manifest.
struct Written {
    line: String,
    manifest: String,
}

/// What `classify` writes for a run over `rows` that compares `fields`: every field whose
/// two values differ is an entry under the class `expected` names for it, the entries in
/// the order of their keys and, on a key, of the fields. `summary` is the count of entries
/// under each class, which the entries must give.
fn written_for(
    held: &[String],
    rows: &[Row],
    fields: &[&str],
    expected: &[Classed],
    summary: &[(&str, usize)],
) -> Written {
    let mut ordered: Vec<(usize, &Row)> = rows.iter().enumerate().collect();
    ordered.sort_by_key(|(_, row)| (&row.location, &row.allele, &row.feature, &row.terms));
    let mut entries = Vec::new();
    let mut counts: Vec<(&str, usize)> = Vec::new();
    for (ordinal, row) in ordered {
        for field in fields {
            let vep = extra_value(&row.extras[0], field);
            let vep_rs = extra_value(&row.extras[1], field);
            if vep == vep_rs {
                continue;
            }
            let classed = expected
                .iter()
                .find(|classed| classed.name == row.name && classed.field == *field)
                .unwrap_or_else(|| panic!("no class is expected of {field} of {}", row.name));
            let entry = field_divergence(
                &key(
                    &row.location,
                    &row.allele,
                    transcript(&row.feature),
                    &[ordinal],
                ),
                &row.terms,
                (field, vep, vep_rs),
                classed.class,
            );
            entries.push(match classed.second {
                Some(second) => with_member(
                    &entry,
                    &format!(r#""second_expected_divergence": "{second}""#),
                ),
                None => entry,
            });
            match counts.iter_mut().find(|(class, _)| *class == classed.class) {
                Some((_, count)) => *count += 1,
                None => counts.push((classed.class, 1)),
            }
        }
    }
    assert_eq!(
        entries.len(),
        expected.len(),
        "a class is expected of a field that does not differ"
    );
    let sorted = |counts: &[(&str, usize)]| {
        let mut counts: Vec<(String, usize)> = counts
            .iter()
            .map(|(class, count)| (class.to_string(), *count))
            .collect();
        counts.sort();
        counts
    };
    assert_eq!(sorted(&counts), sorted(summary));
    let joined = |quote: char| {
        let items: Vec<String> = counts
            .iter()
            .map(|(class, count)| format!("{quote}{class}{quote}: {count}"))
            .collect();
        format!("{{{}}}", items.join(", "))
    };
    let line = counted((0, "{}"), (0, "{}"), "[]");
    Written {
        line: format!(
            "{}; field divergences {} ({})\n",
            line.strip_suffix('\n').unwrap(),
            entries.len(),
            joined('\'')
        ),
        manifest: obj(&[
            &records(held),
            &nothing_written(),
            &fields_written(fields, &entries, &joined('"')),
        ]),
    }
}

/// Writes the manifest of the records `held` and the two outputs of `rows`, and runs
/// `classify` on `fields`.
fn run_fields(case: &Case, held: &[String], rows: &[Row], fields: &[&str]) -> common::Run {
    case.manifest(&obj(&[&records(held)]));
    let output = |engine: usize| -> String {
        rows.iter()
            .map(|row| coding_row(row, &row.extras[engine]))
            .collect()
    };
    case.outputs(&output(0), &output(1));
    let mut args = Vec::new();
    for field in fields {
        args.extend(["--field", field]);
    }
    case.classify(&args)
}

/// Runs `classify` on the Downstream fields over `rows` and holds the run to its exit
/// code, its line of standard output and its manifest.
fn holds(
    case: &Case,
    held: &[String],
    rows: &[Row],
    expected: &[Classed],
    summary: &[(&str, usize)],
) {
    let written = written_for(held, rows, &DOWNSTREAM_RUN, expected, summary);
    assert_eq!(
        run_fields(case, held, rows, &DOWNSTREAM_RUN),
        clean(0, &written.line)
    );
    assert_eq!(case.after(), written.manifest);
}

/// The Extra column of a frameshift row the Downstream plugin annotated: the row's STRAND
/// and the plugin's two fields.
fn extra(strand: i32, protein: &str, change: i32) -> String {
    format!("IMPACT=HIGH;STRAND={strand};DownstreamProtein={protein};ProteinLengthChange={change}")
}

/// The two Extra columns of a frameshift row whose `HGVS_OFFSET` both engines print:
/// Ensembl VEP's with its `DownstreamProtein` and `ProteinLengthChange`, then vep-rs's.
fn shifted(strand: i32, offset: i32, vep: (&str, i32), vep_rs: (&str, i32)) -> [String; 2] {
    [vep, vep_rs]
        .map(|(protein, change)| format!("{};HGVS_OFFSET={offset}", extra(strand, protein, change)))
}

/// A row of a run over the five cached transcripts: its record's name, its transcript,
/// its `CDS_position`, both Extra columns, its allele and its consequence terms.
struct Plain {
    name: &'static str,
    feature: &'static str,
    cds: &'static str,
    extras: [String; 2],
    allele: &'static str,
    terms: &'static str,
}

/// A frameshift row of a deletion.
fn plain(
    name: &'static str,
    feature: &'static str,
    cds: &'static str,
    extras: [String; 2],
) -> Plain {
    Plain {
        name,
        feature,
        cds,
        extras,
        allele: "-",
        terms: "frameshift_variant",
    }
}

/// Runs `classify` over one record per row, a one-base deletion at 21:100, 21:200 and so
/// on, and holds the run to the classes expected. The corpus cache holds the first
/// transcript as `ENST1`, as `ENST4` and, under codon table 7, as `ENST5`, the second
/// transcript as `ENST2` and the third, under codon table 2, as `ENST3`.
fn plain_run_holds(tag: &str, rows: Vec<Plain>, expected: &[Classed], summary: &[(&str, usize)]) {
    let case = Case::new(tag);
    let first = format!(r#""translateable_seq": "{CODING}", "three_prime_utr": "{UTR}""#);
    let cached = |stable_id: &str, members: &str| {
        format!(
            r#"{{"stable_id": "{stable_id}", "variation_effect_feature_cache": {{{members}}}}}"#
        )
    };
    let transcripts = [
        cached("ENST1", &format!(r#"{first}, "codon_table": 1"#)),
        cached(
            "ENST2",
            &format!(r#""translateable_seq": "{STOP_CODING}", "three_prime_utr": "{UTR}""#),
        ),
        cached(
            "ENST3",
            &format!(
                r#""translateable_seq": "{TABLE_TWO_CODING}", "three_prime_utr": "{TABLE_TWO_UTR}", "codon_table": "2""#
            ),
        ),
        cached("ENST4", &first),
        cached("ENST5", &format!(r#"{first}, "codon_table": 7"#)),
    ];
    case.write_in_corpus(
        "json_cache/transcripts/21/1-1000000.json",
        format!("[{}]", transcripts.join(", ")).as_bytes(),
    );
    let position = |ordinal: usize| 100 * (ordinal + 1);
    let held: Vec<String> = (0..rows.len())
        .map(|ordinal| {
            record(&format!(
                "21\t{}\t{}\tAT\tA\t.\tPASS\t.",
                position(ordinal),
                rows[ordinal].name
            ))
        })
        .collect();
    let rows: Vec<Row> = rows
        .into_iter()
        .enumerate()
        .map(|(ordinal, row)| Row {
            name: row.name.to_string(),
            location: format!("21:{}", position(ordinal) + 1),
            allele: row.allele.to_string(),
            feature: row.feature.to_string(),
            terms: row.terms.to_string(),
            cds: row.cds.to_string(),
            extras: row.extras,
        })
        .collect();
    holds(&case, &held, &rows, expected, summary);
}

/// A row of a run whose every row has a transcript of its own, laid on the genome by its
/// mapper pairs: its record's name, the pairs (`None` leaves the transcript out of the
/// cache), the record's POS, REF and ALT, the row's `CDS_position`, both Extra columns,
/// the transcript's coding sequence and the keys its effect cache lacks.
struct Mapped {
    name: &'static str,
    pairs: Option<&'static [[i64; 5]]>,
    record: (i64, &'static str, &'static str),
    cds: &'static str,
    extras: [String; 2],
    coding: &'static str,
    without: &'static [&'static str],
}

/// A row on a transcript of the first coding sequence with every key of its effect cache.
fn mapped(
    name: &'static str,
    pairs: &'static [[i64; 5]],
    record: (i64, &'static str, &'static str),
    cds: &'static str,
    extras: [String; 2],
) -> Mapped {
    Mapped {
        name,
        pairs: Some(pairs),
        record,
        cds,
        extras,
        coding: CODING,
        without: &[],
    }
}

/// Runs `classify` over one record per row and holds the run to the classes expected. REF
/// and ALT share their first base, and the row's Location and Allele are the ones Ensembl
/// VEP prints for such a record. The transcript of row N is `ENST<N>`: three 5' UTR bases,
/// the row's coding sequence with its peptide, and the 3' UTR.
fn mapped_run_holds(tag: &str, rows: Vec<Mapped>, expected: &[Classed], summary: &[(&str, usize)]) {
    let case = Case::new(tag);
    let mut transcripts = Vec::new();
    let mut held = Vec::new();
    let mut written = Vec::new();
    for (ordinal, row) in rows.into_iter().enumerate() {
        let feature = format!("ENST{}", ordinal + 1);
        let (pos, reference, alternate) = row.record;
        if let Some(pairs) = row.pairs {
            let peptide = translate_codons(row.coding, 1).unwrap();
            let pairs: Vec<Value> = pairs
                .iter()
                .map(|pair| {
                    let pair: Object = MAPPER_PAIR_KEYS
                        .iter()
                        .zip(pair)
                        .map(|(key, value)| (*key, value.to_string()))
                        .collect();
                    Value::Object(pair)
                })
                .collect();
            let mapper: Object = [("pairs", Value::List(pairs))].into_iter().collect();
            let mut cache = Object::new();
            cache.insert("translateable_seq", row.coding);
            cache.insert("three_prime_utr", UTR);
            cache.insert("codon_table", "1");
            cache.insert("peptide", peptide.trim_end_matches('*'));
            cache.insert("mapper", mapper);
            for key in row.without {
                cache.remove(key);
            }
            let mut cached = Object::new();
            cached.insert("stable_id", feature.as_str());
            cached.insert("cdna_coding_start", 4i64);
            cached.insert("cdna_coding_end", 3 + row.coding.len());
            cached.insert("variation_effect_feature_cache", cache);
            transcripts.push(Value::Object(cached));
        }
        held.push(record(&format!(
            "21\t{pos}\t{}\t{reference}\t{alternate}\t.\tPASS\t.",
            row.name
        )));
        let ends = [pos + 1, pos + reference.len() as i64 - 1];
        let (low, high) = (ends[0].min(ends[1]), ends[0].max(ends[1]));
        written.push(Row {
            name: row.name.to_string(),
            location: if low == high {
                format!("21:{low}")
            } else {
                format!("21:{low}-{high}")
            },
            allele: match &alternate[1..] {
                "" => "-".to_string(),
                added => added.to_string(),
            },
            feature,
            terms: "frameshift_variant".to_string(),
            cds: row.cds.to_string(),
            extras: row.extras,
        });
    }
    case.write_in_corpus(
        "json_cache/transcripts/21/1-1000000.json",
        dumps(&Value::List(transcripts), &Style::PLAIN).as_bytes(),
    );
    holds(&case, &held, &written, expected, summary);
}

#[test]
fn classify_names_the_organelle_peptide_path_class() {
    // An `HGVSp` difference on a transcript the corpus cache files under a codon table
    // other than the standard one is `organelle_peptide_paths_translation_table` when the
    // row's consequence set and the descriptions place it in one of the four
    // bare-translation paths: any frameshift row (a stop distance, a `Ter` at an immediate
    // stop, a lost-start shape), an extension count on a `stop_lost` row (the deletion of
    // an empty reference window among them), a duplication read as an insertion. The same
    // pairs on a standard-code transcript, a substitution, a `delins` and a plain deletion
    // on the organelle transcript (their residues come through `peptide`, table-aware on
    // both engines) and a pair naming two proteins stay `unexplained_residual`.
    let case = Case::new("organelle");
    let held = [
        record("MT\t3313\tfs\tAT\tA\t.\tPASS\t."),
        record("MT\t3399\tdup\tA\tAATA\t.\tPASS\t."),
        record("MT\t7443\text\tA\tG\t.\tPASS\t."),
        record("MT\t3308\tsub\tT\tC\t.\tPASS\t."),
        record("MT\t3400\tname\tC\tCG\t.\tPASS\t."),
        record("MT\t4260\tstop\tCG\tC\t.\tPASS\t."),
        record("MT\t4470\tlost\tGT\tG\t.\tPASS\t."),
        record("MT\t5000\tdelins\tTGCTA\tTA\t.\tPASS\t."),
        record("21\t100\tnuclear\tAT\tA\t.\tPASS\t."),
    ];
    case.write_in_corpus(
        "json_cache/transcripts/MT/1-1000000.json",
        br#"[{"stable_id": "ENSTMT1", "variation_effect_feature_cache": {"codon_table": 2}}, {"stable_id": "ENSTMT2", "variation_effect_feature_cache": {"codon_table": "2"}}]"#,
    );
    case.write_in_corpus(
        "json_cache/transcripts/21/1-1000000.json",
        br#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"codon_table": 1}}]"#,
    );
    let row = |name: &str,
               location: &str,
               allele: &str,
               feature: &str,
               terms: &str,
               vep: &str,
               vep_rs: &str| {
        Row {
            name: name.to_string(),
            location: location.to_string(),
            allele: allele.to_string(),
            feature: feature.to_string(),
            terms: terms.to_string(),
            cds: "-".to_string(),
            extras: [vep, vep_rs].map(|protein| format!("IMPACT=HIGH;STRAND=1;HGVSp={protein}")),
        }
    };
    let frameshift = "frameshift_variant";
    let rows = [
        row(
            "fs",
            "MT:3314",
            "-",
            "ENSTMT1",
            frameshift,
            "ENSPMT1.1:p.Met3ArgfsTer12",
            "ENSPMT1.1:p.Met3Ter",
        ),
        row(
            "dup",
            "MT:3399-3400",
            "ATA",
            "ENSTMT1",
            "inframe_insertion",
            "ENSPMT1.1:p.Met31_Gln32insMet",
            "ENSPMT1.1:p.Met31dup",
        ),
        row(
            "ext",
            "MT:7443",
            "G",
            "ENSTMT2",
            "stop_lost",
            "ENSPMT2.1:p.Ter514GlyextTer1",
            "ENSPMT2.1:p.Ter514GlyextTer3",
        ),
        row(
            "sub",
            "MT:3308",
            "C",
            "ENSTMT1",
            "missense_variant",
            "ENSPMT1.1:p.Met1Thr",
            "ENSPMT1.1:p.Met1Ile",
        ),
        row(
            "name",
            "MT:3400-3401",
            "G",
            "ENSTMT1",
            frameshift,
            "ENSPMT1.1:p.Gln32ArgfsTer20",
            "ENSPMT1.2:p.Gln32ArgfsTer5",
        ),
        row(
            "stop",
            "MT:4261",
            "-",
            "ENSTMT1",
            frameshift,
            "ENSPMT1.1:p.Asp3Ter",
            "ENSPMT1.1:p.Arg2Ter",
        ),
        row(
            "lost",
            "MT:4471",
            "-",
            "ENSTMT2",
            "frameshift_variant,start_lost",
            "ENSPMT2.1:p.Thr2_?1",
            "ENSPMT2.1:p.Lys3_?1",
        ),
        row(
            "delins",
            "MT:5001-5003",
            "-",
            "ENSTMT1",
            "inframe_deletion",
            "ENSPMT1.1:p.Ala4delinsGly",
            "ENSPMT1.1:p.Ala4delinsSer",
        ),
        row(
            "nuclear",
            "21:101",
            "-",
            "ENST1",
            frameshift,
            "ENSP1.1:p.Leu5ArgfsTer12",
            "ENSP1.1:p.Leu5Ter",
        ),
    ];
    let hgvsp = |name: &'static str, class: &'static str| Classed {
        field: "HGVSp",
        ..protein(name, class)
    };
    let expected = [
        hgvsp("fs", ORGANELLE),
        hgvsp("dup", ORGANELLE),
        hgvsp("ext", ORGANELLE),
        hgvsp("stop", ORGANELLE),
        hgvsp("lost", ORGANELLE),
        hgvsp("sub", UNEXPLAINED),
        hgvsp("name", UNEXPLAINED),
        hgvsp("delins", UNEXPLAINED),
        hgvsp("nuclear", UNEXPLAINED),
    ];
    let written = written_for(
        &held,
        &rows,
        &["HGVSp"],
        &expected,
        &[(ORGANELLE, 5), (UNEXPLAINED, 4)],
    );
    let ran = run_fields(&case, &held, &rows, &["HGVSp"]);
    assert_eq!(ran, clean(0, &written.line));
    assert_eq!(case.after(), written.manifest);
    let cache_dir = format!("{}/json_cache", case.corpus);
    let mut tables: Vec<(String, i64)> = organelle_codon_tables(&cache_dir)
        .unwrap()
        .into_iter()
        .collect();
    tables.sort();
    assert_eq!(
        tables,
        [("ENSTMT1".to_string(), 2), ("ENSTMT2".to_string(), 2)]
    );
    let shapes = [
        "Met31dup",
        "Met31_Gln32insMet",
        "Met31_Gln32dup",
        "Ser513del",
        "Ala4delinsGly",
        "Met3ArgfsTer12",
        "Met1?",
    ]
    .map(peptide_description_shape);
    assert_eq!(
        shapes,
        ["dup", "ins", "dup", "other", "other", "other", "other"]
    );
    assert!(!organelle_peptide_path_divergence(
        "stop_lost",
        "Ter514Gly",
        "Ter514Trp"
    ));
    assert!(organelle_peptide_path_divergence(
        "stop_lost",
        "Ter4delextTer?",
        "Ter4delextTer2"
    ));
    assert!(organelle_peptide_path_divergence(
        "protein_altering_variant",
        "Met31dup",
        "Met31_Gln32insMet"
    ));
    assert!(!organelle_peptide_path_divergence(
        "inframe_deletion",
        "Ser5_Lys6del",
        "Ser5del"
    ));
}

#[test]
fn classify_names_the_downstream_class_of_a_deletion_cut_to_the_last_coding_base() {
    // A deletion with no inserted base from a codon's third base to the last coding base
    // is `downstream_utr_frame_past_coding_end`. Deleting bases 24 to 27 leaves 23 bases:
    // Ensembl VEP's boundary, base 24, lies past them, so it reads the UTR alone, and
    // vep-rs reads the last two bases with the UTR. Deleting base 27 leaves two bases of
    // the stop codon, which take the UTR's first base and make a stop, so vep-rs writes
    // the empty value where Ensembl VEP reads the UTR alone.
    plain_run_holds(
        "coding-end",
        vec![
            plain(
                "to_the_end",
                "ENST1",
                "24-27",
                [extra(1, "GPK", 4), extra(1, "GGPNS", 4)],
            ),
            plain(
                "last_base",
                "ENST1",
                "27",
                [extra(1, "GPK", 0), extra(1, "", 0)],
            ),
        ],
        &[
            protein("to_the_end", UTR_FRAME_PAST_CODING_END),
            protein("last_base", UTR_FRAME_PAST_CODING_END),
        ],
        &[(UTR_FRAME_PAST_CODING_END, 2)],
    );
}

#[test]
fn classify_names_the_downstream_class_of_a_deletion_from_a_codons_third_base() {
    // A deletion whose first base is a codon's third and that leaves coding bases after it
    // is `downstream_start_past_own_codon`: vep-rs starts at the changed codon and Ensembl
    // VEP at the codon after it. Deleting base 6 keeps the residue, the codon taking the
    // next base; deleting bases 6 and 7 changes it; deleting base 6 of the second
    // transcript turns its second codon into a stop, so vep-rs writes the empty value and
    // Ensembl VEP reads past the stop.
    plain_run_holds(
        "third-base",
        vec![
            plain(
                "kept_residue",
                "ENST1",
                "6",
                [extra(1, "LLLSVVKGPNS", 5), extra(1, "KLLLSVVKGPNS", 5)],
            ),
            plain(
                "changed_residue",
                "ENST1",
                "6-7",
                [extra(1, "SCCPWLRAQIA", 5), extra(1, "NSCCPWLRAQIA", 5)],
            ),
            plain(
                "stop",
                "ENST2",
                "6",
                [extra(1, "ENCKGPNS", -4), extra(1, "", -4)],
            ),
        ],
        &[
            protein("kept_residue", START_PAST_OWN_CODON),
            protein("changed_residue", START_PAST_OWN_CODON),
            protein("stop", START_PAST_OWN_CODON),
        ],
        &[(START_PAST_OWN_CODON, 3)],
    );
}

#[test]
fn classify_reads_a_downstream_class_on_the_printed_span_alone() {
    // The two classes of a `DownstreamProtein` are read on the `CDS_position` the row
    // prints, whatever its `HGVS_OFFSET`, since vep-rs computes on that span. Bases 11 and
    // 12 are one base twice: a row printing base 11, a codon's second, carries the pair of
    // a deletion of base 12, a codon's third. No rule gives that pair on base 11, with an
    // offset of 1 or 2 or none, on either strand, so on a transcript whose cache has no
    // mapper pairs to move the span by the row stays `unexplained_residual`.
    let pair = |strand: i32, offset: &str| {
        ["LSVVKGPNS", "LLSVVKGPNS"].map(|value| format!("{}{offset}", extra(strand, value, 5)))
    };
    plain_run_holds(
        "printed-span",
        vec![
            plain("shifted", "ENST1", "11", pair(1, ";HGVS_OFFSET=1")),
            plain(
                "shifted_reverse",
                "ENST4",
                "11",
                pair(-1, ";HGVS_OFFSET=-1"),
            ),
            plain("unshifted", "ENST1", "11", pair(1, "")),
            plain("shifted_further", "ENST1", "11", pair(1, ";HGVS_OFFSET=2")),
        ],
        &[
            protein("shifted", UNEXPLAINED),
            protein("shifted_reverse", UNEXPLAINED),
            protein("unshifted", UNEXPLAINED),
            protein("shifted_further", UNEXPLAINED),
        ],
        &[(UNEXPLAINED, 4)],
    );
}

#[test]
fn classify_reads_a_downstream_allele_on_the_transcripts_strand() {
    // A deletion with inserted bases on a reverse-strand transcript puts the reverse
    // complement of the row's allele into the coding sequence: the two bases of the allele
    // replace bases 6 to 8 with their reverse complement, vep-rs reads from the changed
    // codon and Ensembl VEP from the codon after it. The same row on a forward-strand
    // transcript, where the allele goes in as written, gives another pair and stays
    // `unexplained_residual`.
    let on = |name: &'static str, feature: &'static str, strand: i32| Plain {
        allele: "AC",
        ..plain(
            name,
            feature,
            "6-8",
            [
                extra(strand, "FLLSVVKGPNS", 5),
                extra(strand, "KFLLSVVKGPNS", 5),
            ],
        )
    };
    plain_run_holds(
        "strand",
        vec![on("reverse", "ENST4", -1), on("forward", "ENST1", 1)],
        &[
            protein("reverse", START_PAST_OWN_CODON),
            protein("forward", UNEXPLAINED),
        ],
        &[(START_PAST_OWN_CODON, 1), (UNEXPLAINED, 1)],
    );
}

#[test]
fn classify_leaves_a_downstream_pair_the_sequences_do_not_give_unexplained() {
    // The pair of a deletion of base 6 stays `unexplained_residual` on a row whose span
    // opens on a codon's first or second base, has a gap for an end or runs past the
    // coding sequence; with either value changed or absent; on a row without
    // `frameshift_variant` or with a splice term; and on a transcript the cache lacks. A
    // `ProteinLengthChange` that differs is in neither class, on a row whose
    // `DownstreamProtein` is.
    let pair = || [extra(1, "LLLSVVKGPNS", 5), extra(1, "KLLLSVVKGPNS", 5)];
    let [vep, vep_rs] = pair();
    let absent = || "IMPACT=HIGH;STRAND=1".to_string();
    let unexplained = [
        "first_base",
        "second_base",
        "gap",
        "past_the_end",
        "other_vep_rs_value",
        "other_vep_value",
        "no_vep_rs_value",
        "no_vep_value",
        "inframe",
        "splice",
        "uncached",
    ];
    let mut expected: Vec<Classed> = unexplained
        .into_iter()
        .map(|name| protein(name, UNEXPLAINED))
        .collect();
    expected.extend([
        length_change("no_vep_rs_value", UNEXPLAINED),
        length_change("no_vep_value", UNEXPLAINED),
        protein("length", START_PAST_OWN_CODON),
        length_change("length", UNEXPLAINED),
    ]);
    plain_run_holds(
        "not-given",
        vec![
            plain("first_base", "ENST1", "4", pair()),
            plain("second_base", "ENST1", "5", pair()),
            plain("gap", "ENST1", "6-?", pair()),
            plain(
                "past_the_end",
                "ENST1",
                "24-28",
                [extra(1, "GPK", 4), extra(1, "GGPNS", 4)],
            ),
            plain(
                "other_vep_rs_value",
                "ENST1",
                "6",
                [vep.clone(), extra(1, "KLLLSVVKGPN", 5)],
            ),
            plain(
                "other_vep_value",
                "ENST1",
                "6",
                [extra(1, "LLSVVKGPNS", 5), vep_rs.clone()],
            ),
            plain("no_vep_rs_value", "ENST1", "6", [vep.clone(), absent()]),
            plain("no_vep_value", "ENST1", "6", [absent(), vep_rs]),
            Plain {
                terms: "inframe_deletion",
                ..plain("inframe", "ENST1", "6", pair())
            },
            Plain {
                terms: "frameshift_variant,splice_region_variant",
                ..plain("splice", "ENST1", "6", pair())
            },
            plain("uncached", "ENST9", "6", pair()),
            plain("length", "ENST1", "6", [vep, extra(1, "KLLLSVVKGPNS", 6)]),
        ],
        &expected,
        &[(UNEXPLAINED, 14), (START_PAST_OWN_CODON, 1)],
    );
}

#[test]
fn classify_reads_a_downstream_pair_under_the_transcripts_codon_table() {
    // Both values are read under the transcript's codon table. Deleting base 6 of the
    // third transcript, filed under table 2, leaves a read that crosses a codon table 2
    // reads as tryptophan and ends at one it reads as a stop: the row is
    // `downstream_start_past_own_codon`. Under the standard table the first of the two is
    // a stop, which gives another pair: the row carrying that pair stays
    // `unexplained_residual`. So does the first transcript's pair on a transcript filed
    // under table 7, an id `Bio::Tools::CodonTable` holds no code for.
    plain_run_holds(
        "codon-table",
        vec![
            plain(
                "table_two",
                "ENST3",
                "6",
                [extra(1, "WP", 0), extra(1, "NWP", 0)],
            ),
            plain(
                "standard_table_values",
                "ENST3",
                "6",
                [extra(1, "", -2), extra(1, "N", -2)],
            ),
            plain(
                "table_seven",
                "ENST5",
                "6",
                [extra(1, "LLLSVVKGPNS", 5), extra(1, "KLLLSVVKGPNS", 5)],
            ),
        ],
        &[
            protein("table_two", START_PAST_OWN_CODON),
            protein("standard_table_values", UNEXPLAINED),
            protein("table_seven", UNEXPLAINED),
        ],
        &[(START_PAST_OWN_CODON, 1), (UNEXPLAINED, 2)],
    );
}

#[test]
fn classify_names_the_class_of_a_downstream_value_read_at_the_hgvs_shifted_span() {
    // Under `--hgvs` Ensembl VEP's plugin reads the span the notation's 3' shift leaves,
    // the record's own moved by the row's `HGVS_OFFSET`, and vep-rs the record's own:
    // where that changes Ensembl VEP's value the row is
    // `plugin_hgvs_shifted_cds_coordinates`. Bases 4 to 7 are one base four times. A
    // deletion of base 4 is written at base 7, past the codon's end: vep-rs reads from the
    // second codon and Ensembl VEP from the third. A base inserted before base 4 is
    // written after base 7. Two bases inserted after base 8 are written after base 11, the
    // allele turned once per base of the shift. A deletion of base 11, a codon's second
    // base, is written at base 12, the same codon's third: both rules read one value on
    // base 11 and Ensembl VEP's starts one codon later on base 12.
    mapped_run_holds(
        "shifted",
        vec![
            mapped(
                "deletion",
                ONE_EXON,
                (1006, "GA", "G"),
                "4",
                shifted(1, 3, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5)),
            ),
            mapped(
                "insertion",
                ONE_EXON,
                (1006, "G", "GA"),
                "3-4",
                shifted(1, 4, ("NSCCPWLRAQIA", 6), ("KNSCCPWLRAQIA", 6)),
            ),
            mapped(
                "turned_allele",
                ONE_EXON,
                (1011, "C", "CTC"),
                "8-9",
                shifted(1, 3, ("LLLSVVKGPNS", 6), ("TLLLSVVKGPNS", 6)),
            ),
            mapped(
                "inside_one_codon",
                ONE_EXON,
                (1013, "CT", "C"),
                "11",
                shifted(1, 1, ("LSVVKGPNS", 5), ("LLSVVKGPNS", 5)),
            ),
        ],
        &[
            protein("deletion", HGVS_SHIFTED_CDS_COORDINATES),
            protein("insertion", HGVS_SHIFTED_CDS_COORDINATES),
            protein("turned_allele", HGVS_SHIFTED_CDS_COORDINATES),
            protein("inside_one_codon", HGVS_SHIFTED_CDS_COORDINATES),
        ],
        &[(HGVS_SHIFTED_CDS_COORDINATES, 4)],
    );
}

#[test]
fn classify_reads_a_shifted_downstream_span_on_a_reverse_strand_transcript() {
    // A reverse-strand row prints a negative `HGVS_OFFSET`, the shift along the genome,
    // and its allele turns its length less the shift times before it is read on the
    // transcript's strand. The deletion of coding base 4, genomic base 1039, is written at
    // coding base 7, three genomic bases lower. Four bases inserted between coding bases
    // 11 and 12 are written one base on, turned three times. The same row with the
    // offset's sign turned moves the other way along the transcript, gives another
    // sequence and stays `unexplained_residual`.
    let insertion = |offset: i32| shifted(-1, offset, ("HCCPWLRAQIA", 7), ("LHCCPWLRAQIA", 7));
    mapped_run_holds(
        "reverse-strand",
        vec![
            mapped(
                "deletion",
                REVERSE_EXON,
                (1038, "TT", "T"),
                "4",
                shifted(-1, -3, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5)),
            ),
            mapped(
                "insertion",
                REVERSE_EXON,
                (1031, "A", "AGTGA"),
                "11-12",
                insertion(-1),
            ),
            mapped(
                "other_direction",
                REVERSE_EXON,
                (1031, "A", "AGTGA"),
                "11-12",
                insertion(1),
            ),
        ],
        &[
            protein("deletion", HGVS_SHIFTED_CDS_COORDINATES),
            protein("insertion", HGVS_SHIFTED_CDS_COORDINATES),
            protein("other_direction", UNEXPLAINED),
        ],
        &[(HGVS_SHIFTED_CDS_COORDINATES, 2), (UNEXPLAINED, 1)],
    );
}

#[test]
fn classify_moves_a_shifted_downstream_span_through_the_transcripts_mapper_pairs() {
    // The shifted span is the record's genomic span moved by the offset and placed on the
    // coding sequence by the transcript's mapper pairs, so a shift across an intron moves
    // the CDS pair by fewer bases than the offset. A base inserted between coding bases 11
    // and 12, both that base, is written three genomic bases on, past the last exon base
    // and a two-base intron, before coding base 13. An offset of 2 leaves both flanks in
    // the intron, where the plugin reads the record's own span, and on one exon the offset
    // of 3 moves the pair three coding bases: no rule gives the pair on either row.
    let values = |offset: i32| shifted(1, offset, ("CCPWLRAQIA", 6), ("LCCPWLRAQIA", 6));
    let record = (1014, "T", "TT");
    mapped_run_holds(
        "mapper-pairs",
        vec![
            mapped("across_the_intron", TWO_EXONS, record, "11-12", values(3)),
            mapped("into_the_intron", TWO_EXONS, record, "11-12", values(2)),
            mapped("one_exon", ONE_EXON, record, "11-12", values(3)),
        ],
        &[
            protein("across_the_intron", HGVS_SHIFTED_CDS_COORDINATES),
            protein("into_the_intron", UNEXPLAINED),
            protein("one_exon", UNEXPLAINED),
        ],
        &[(HGVS_SHIFTED_CDS_COORDINATES, 1), (UNEXPLAINED, 2)],
    );
}

#[test]
fn classify_reads_a_downstream_insertion_on_a_seam_of_the_model() {
    // Where two mapper pairs touch on the genome and the model carries bases between them,
    // an insertion at that point replaces those bases (Mapper.pm:504-558,
    // Downstream.pm:108-110). A base inserted between coding bases 11 and 12 is written
    // one base on, at a seam that carries coding base 13: there Ensembl VEP puts the base
    // in that base's place, a sequence of the reference's length, where vep-rs inserts it
    // at the record's own point and reads 14 residues, 6 more than the peptide's 8: both
    // fields are `plugin_hgvs_shifted_cds_coordinates`. A base inserted at a seam that
    // carries coding bases 12 and 13 prints that span in order, a codon's third base
    // first: `downstream_start_past_own_codon` on the row's own span. The same column on a
    // model without the seam is an insertion's two flanks, which no rule gives that pair
    // for.
    let on_the_seam = || [extra(1, "LSVVKGPNS", 5), extra(1, "LLSVVKGPNS", 5)];
    mapped_run_holds(
        "seam",
        vec![
            mapped(
                "onto_the_seam",
                ONE_BASE_SEAM,
                (1014, "T", "TT"),
                "11-12",
                shifted(1, 1, ("SVRG", 0), ("LCCPWLRAQIA", 6)),
            ),
            mapped(
                "on_the_seam",
                TWO_BASE_SEAM,
                (1014, "T", "TA"),
                "12-13",
                on_the_seam(),
            ),
            mapped(
                "no_seam",
                ONE_EXON,
                (1015, "T", "TA"),
                "12-13",
                on_the_seam(),
            ),
        ],
        &[
            protein("onto_the_seam", HGVS_SHIFTED_CDS_COORDINATES),
            length_change("onto_the_seam", HGVS_SHIFTED_CDS_COORDINATES),
            protein("on_the_seam", START_PAST_OWN_CODON),
            protein("no_seam", UNEXPLAINED),
        ],
        &[
            (HGVS_SHIFTED_CDS_COORDINATES, 2),
            (START_PAST_OWN_CODON, 1),
            (UNEXPLAINED, 1),
        ],
    );
}

/// The row of a deletion of base 9 of the second transcript, a codon's third, written at
/// base 12, with Ensembl VEP's `DownstreamProtein` given.
fn third_base_of_the_second_transcript(name: &'static str, vep: &str) -> Mapped {
    Mapped {
        coding: STOP_CODING,
        ..mapped(
            name,
            STOP_EXON,
            (1011, "GA", "G"),
            "9",
            shifted(1, 3, (vep, 5), ("GNCKGPNS", 5)),
        )
    }
}

#[test]
fn classify_keeps_the_downstream_class_of_a_span_whose_shifted_reading_is_its_own() {
    // A row whose own span opens on a codon's third base and whose shifted reading gives
    // Ensembl VEP's own-span value stays `downstream_start_past_own_codon`: the deletion
    // of base 6 is written at base 7 and Ensembl VEP reads from the third codon on either,
    // vep-rs from the second. Where the shift and the start rule both change the value the
    // row is `plugin_hgvs_shifted_cds_coordinates`: on the second transcript bases 9 to 12
    // are one base four times, and the deletion of base 9, a codon's third, is written at
    // base 12; vep-rs reads from the third codon, Ensembl VEP's rule from the fourth on
    // base 9 and from the fifth on base 12. The value that rule gives on base 9 is not
    // what Ensembl VEP writes for the row, so a row carrying it stays
    // `unexplained_residual`.
    mapped_run_holds(
        "own-reading",
        vec![
            mapped(
                "own_value",
                ONE_EXON,
                (1008, "AA", "A"),
                "6",
                shifted(1, 1, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5)),
            ),
            third_base_of_the_second_transcript("both_rules", "CKGPNS"),
            third_base_of_the_second_transcript("start_rule_alone", "NCKGPNS"),
        ],
        &[
            protein("own_value", START_PAST_OWN_CODON),
            Classed {
                second: Some(START_PAST_OWN_CODON),
                ..protein("both_rules", HGVS_SHIFTED_CDS_COORDINATES)
            },
            protein("start_rule_alone", UNEXPLAINED),
        ],
        &[
            (START_PAST_OWN_CODON, 1),
            (HGVS_SHIFTED_CDS_COORDINATES, 1),
            (UNEXPLAINED, 1),
        ],
    );
}

#[test]
fn classify_records_the_second_class_of_a_downstream_value_two_causes_change() {
    // Where the shift and the start rule both change a `DownstreamProtein`, the entry is
    // `plugin_hgvs_shifted_cds_coordinates` and names `downstream_start_past_own_codon`
    // under `second_expected_divergence`, the class the row's own span has. On the second
    // transcript the deletion of base 9, a codon's third, is written at base 12: vep-rs,
    // Ensembl VEP's rule on base 12 and that rule on base 9 give three values. An entry
    // one cause changes carries no second class: the deletion of base 4 written at base 7,
    // where vep-rs's value is Ensembl VEP's on base 4, and the deletion of base 6 written
    // at base 7, where Ensembl VEP's rule gives one value on both. The summary counts each
    // entry once, under its own class.
    let one_codon_late = |offset: i32| shifted(1, offset, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5));
    mapped_run_holds(
        "second-class",
        vec![
            third_base_of_the_second_transcript("both_causes", "CKGPNS"),
            mapped(
                "flag_alone",
                ONE_EXON,
                (1006, "GA", "G"),
                "4",
                one_codon_late(3),
            ),
            mapped(
                "start_rule_alone",
                ONE_EXON,
                (1008, "AA", "A"),
                "6",
                one_codon_late(1),
            ),
        ],
        &[
            Classed {
                second: Some(START_PAST_OWN_CODON),
                ..protein("both_causes", HGVS_SHIFTED_CDS_COORDINATES)
            },
            protein("flag_alone", HGVS_SHIFTED_CDS_COORDINATES),
            protein("start_rule_alone", START_PAST_OWN_CODON),
        ],
        &[(HGVS_SHIFTED_CDS_COORDINATES, 2), (START_PAST_OWN_CODON, 1)],
    );
}

#[test]
fn classify_leaves_a_downstream_pair_the_hgvs_shift_does_not_give_unexplained() {
    // The pair of a deletion of base 4 written at base 7 stays `unexplained_residual` on a
    // row without an `HGVS_OFFSET`; with an offset that ends in the 3' UTR, which does not
    // reach the plugin; with one that leaves the span in its own codon, where the shifted
    // reading is the own one; with either value changed or absent; on a transcript the
    // cache lacks or gives no mapper pairs for. A `ProteinLengthChange` has no class with
    // another value on either side, nor on a transcript whose cache holds no peptide,
    // beside a `DownstreamProtein` that has.
    let pair = [("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5)];
    let deletion = |name: &'static str, extras: [String; 2]| {
        mapped(name, ONE_EXON, (1006, "GA", "G"), "4", extras)
    };
    let seam = |name: &'static str, vep: (&str, i32), vep_rs: (&str, i32)| {
        mapped(
            name,
            ONE_BASE_SEAM,
            (1014, "T", "TT"),
            "11-12",
            shifted(1, 1, vep, vep_rs),
        )
    };
    let absent = || "IMPACT=HIGH;STRAND=1;HGVS_OFFSET=3".to_string();
    let [shifted_vep, shifted_vep_rs] = shifted(1, 3, pair[0], pair[1]);
    let on_the_seam = ("LCCPWLRAQIA", 6);
    let rows = vec![
        deletion(
            "no_offset",
            pair.map(|(protein, change)| extra(1, protein, change)),
        ),
        deletion("into_the_utr", shifted(1, 24, pair[0], pair[1])),
        deletion("own_codon", shifted(1, 1, pair[0], pair[1])),
        deletion("other_vep_value", shifted(1, 3, ("LLLSVVKGPN", 5), pair[1])),
        deletion(
            "other_vep_rs_value",
            shifted(1, 3, pair[0], ("KLLLSVVKGPN", 5)),
        ),
        deletion("no_vep_value", [absent(), shifted_vep_rs.clone()]),
        deletion("no_vep_rs_value", [shifted_vep.clone(), absent()]),
        Mapped {
            pairs: None,
            ..deletion("uncached", [shifted_vep.clone(), shifted_vep_rs.clone()])
        },
        Mapped {
            without: &["mapper"],
            ..deletion("no_mapper", [shifted_vep, shifted_vep_rs])
        },
        seam("other_vep_length", ("SVRG", 1), on_the_seam),
        seam("other_vep_rs_length", ("SVRG", 0), ("LCCPWLRAQIA", 7)),
        Mapped {
            without: &["peptide"],
            ..seam("no_peptide", ("SVRG", 0), on_the_seam)
        },
    ];
    let unexplained = [
        "no_offset",
        "into_the_utr",
        "own_codon",
        "other_vep_value",
        "other_vep_rs_value",
        "no_vep_value",
        "no_vep_rs_value",
        "uncached",
        "no_mapper",
    ];
    let lengths = [
        "no_vep_value",
        "no_vep_rs_value",
        "other_vep_length",
        "other_vep_rs_length",
        "no_peptide",
    ];
    let mut expected: Vec<Classed> = unexplained
        .into_iter()
        .map(|name| protein(name, UNEXPLAINED))
        .collect();
    expected.extend(lengths.map(|name| length_change(name, UNEXPLAINED)));
    expected.extend(
        ["other_vep_length", "other_vep_rs_length", "no_peptide"]
            .map(|name| protein(name, HGVS_SHIFTED_CDS_COORDINATES)),
    );
    mapped_run_holds(
        "shift-not-given",
        rows,
        &expected,
        &[(UNEXPLAINED, 14), (HGVS_SHIFTED_CDS_COORDINATES, 3)],
    );
}
