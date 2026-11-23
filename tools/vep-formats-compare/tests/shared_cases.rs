// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The shared readers and pair rules against stored answers.
//!
//! Each file under `tests/data/shared` is named for a function and holds one case a line:
//! the arguments, a tab, the answer the function must give. A text is
//! its UTF-8 bytes, each as two lowercase hex digits, joined by `_`, and an absent text `-`;
//! a truth value is `1` or `0`, an absent one `-`; a list of texts is its length, a space
//! and the hex of each joined by commas; an entry is eight fields (key, Location, allele,
//! Consequence, its fields as `name=value` pairs, indel, insertion, reference); the parts
//! of a region, a placement or a row are joined by `~`. An answer is `ok <value>`,
//! `err <the error's text>` or `uncaught <the text of an exception nothing catches>`. A
//! whole number is an `i64` here, so an answer that needs a longer one is written as the
//! tool gives it: the error of `integer_of`, or the nearest `i64` of a span or a position;
//! `integer_text_of.tsv` holds the answers of `integer_of` in that file as they are.
//!
//! The ignored test reads the files of the same names, with many more lines, under the
//! directory `VEP_FORMATS_COMPARE_SHAREDDIFF_DIR` names.

use std::fmt::Write as _;
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Once;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Result, Uncaught};
use vep_formats_compare::shared;

const FILES: [&str; 45] = [
    "_ten_digits",
    "dash",
    "allele_is_indel",
    "reference_allele_is_empty",
    "uploaded_alleles",
    "uploaded_position",
    "replaced_reference",
    "blank",
    "integer_of",
    "integer_text_of",
    "consequence_set",
    "feature_kind",
    "indel_row",
    "position_bounds",
    "stop_codon_test_coordinates",
    "stop_codon_test_not_made",
    "stop_codon_overlap_genomic_span_pair",
    "length_changing_row",
    "stop_retained_without_stop_codon_test_pair",
    "stop_retained_three_base_stop_codon_test_pair",
    "stop_codon_removed_without_stop_codon_test_pair",
    "insertion_row",
    "inserted_between_neighbours",
    "stop_codon_overlap_insertion_range_pair",
    "transcript_positions",
    "location_span",
    "longer_allele_inside_one_exon",
    "differing_regions",
    "intron_rows",
    "intron_placements",
    "record_span_placement",
    "splice_term_past_the_record_span_pair",
    "vep_rs_pos",
    "_strip_chr",
    "perl_fields",
    "perl_split",
    "sequence_location",
    "canonical_location",
    "_canonical_coordinate",
    "gvcf_site_rows",
    "is_structural_allele",
    "_is_start_cooccurrence_swap",
    "_is_covered_splice_region_swap",
    "EXCLUDING_RULES",
    "read_record_references",
];

fn unhex(hex: &str) -> Vec<u8> {
    if hex.is_empty() {
        return Vec::new();
    }
    hex.split('_')
        .map(|byte| {
            assert_eq!(byte.len(), 2, "two hex digits a byte in {hex}");
            // `from_str_radix` also takes a sign and a capital digit.
            assert!(
                byte.bytes()
                    .all(|digit| matches!(digit, b'0'..=b'9' | b'a'..=b'f')),
                "lowercase hex digits in {hex}"
            );
            u8::from_str_radix(byte, 16).expect("two hex digits")
        })
        .collect()
}

fn text(hex: &str) -> String {
    String::from_utf8(unhex(hex)).expect("the hex of UTF-8 bytes")
}

fn hex_of(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(3 * bytes.len());
    for (at, byte) in bytes.iter().enumerate() {
        if at > 0 {
            out.push('_');
        }
        let _ = write!(out, "{byte:02x}");
    }
    out
}

fn hex(text: &str) -> String {
    hex_of(text.as_bytes())
}

fn opt_text(field: &str) -> Option<String> {
    (field != "-").then(|| text(field))
}

fn truth(field: &str) -> bool {
    field == "1"
}

fn opt_truth(field: &str) -> Option<bool> {
    (field != "-").then(|| truth(field))
}

fn whole(field: &str) -> i64 {
    field.parse().expect("a whole number")
}

fn counted(field: &str) -> Vec<&str> {
    let (count, joined) = field.split_once(' ').expect("a count and its items");
    if count == "0" {
        Vec::new()
    } else {
        joined.split(',').collect()
    }
}

fn texts(field: &str) -> Vec<String> {
    counted(field).into_iter().map(text).collect()
}

fn entry(fields: &[&str]) -> Entry {
    let mut e = Entry {
        key: text(fields[0]),
        location: opt_text(fields[1]),
        allele: text(fields[2]),
        consequence: text(fields[3]),
        indel: opt_truth(fields[5]),
        insertion: opt_truth(fields[6]),
        reference: opt_text(fields[7]),
        ..Entry::default()
    };
    for pair in counted(fields[4]) {
        let (name, value) = pair.split_once('=').expect("a name and a value");
        e.fields.insert(text(name), text(value));
    }
    e
}

fn bit(value: bool) -> String {
    if value { "1" } else { "0" }.to_string()
}

fn opt_bit(value: Option<bool>) -> String {
    value.map_or("-".to_string(), bit)
}

fn eq_text(value: String) -> String {
    format!("={}", hex(&value))
}

fn opt_eq_text(value: Option<String>) -> String {
    value.map_or("-".to_string(), eq_text)
}

fn opt_whole(value: Option<i64>) -> String {
    value.map_or("-".to_string(), |number| number.to_string())
}

fn items(parts: Vec<String>) -> String {
    let hexes: Vec<String> = parts.iter().map(|part| hex(part)).collect();
    format!("{} {}", parts.len(), hexes.join(","))
}

fn fnv64(rendered: &str) -> String {
    let mut value: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in rendered.bytes() {
        value = (value ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
    }
    hex_of(&value.to_be_bytes())
}

/// Keeps the panic of an uncaught exception off the test's output.
fn quiet_uncaught() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let default = panic::take_hook();
        panic::set_hook(Box::new(move |info| {
            if info.payload().downcast_ref::<Uncaught>().is_none() {
                default(info);
            }
        }));
    });
}

fn outcome<T>(call: impl FnOnce() -> Result<T>, written: impl FnOnce(T) -> String) -> String {
    match panic::catch_unwind(AssertUnwindSafe(call)) {
        Ok(Ok(value)) => format!("ok {}", written(value)),
        Ok(Err(Error::Comparator(what))) => format!("err {}", hex(&what)),
        Err(payload) => match payload.downcast_ref::<Uncaught>() {
            Some(Uncaught(what)) => format!("uncaught {}", hex(what)),
            None => panic::resume_unwind(payload),
        },
    }
}

fn plain<T>(call: impl FnOnce() -> T, written: impl FnOnce(T) -> String) -> String {
    outcome(|| Ok(call()), written)
}

fn pair_answer(rule: fn(&Entry, &Entry) -> Result<bool>, f: &[&str]) -> String {
    let (e, a) = (entry(&f[..8]), entry(&f[8..16]));
    outcome(|| rule(&e, &a), bit)
}

/// A list of placements: its length, the digest of the whole list with the parts of each
/// placement joined by `:`, and its first three placements.
fn placements(places: Vec<(bool, i64, i64)>) -> String {
    let written = |separator: &str| -> Vec<String> {
        places
            .iter()
            .map(|(forward, start, end)| {
                [bit(*forward), start.to_string(), end.to_string()].join(separator)
            })
            .collect()
    };
    let first: Vec<String> = written("~").into_iter().take(3).collect();
    format!(
        "{} {} {}",
        places.len(),
        fnv64(&written(":").join(";")),
        first.join(";")
    )
}

static DIRECTORIES: AtomicUsize = AtomicUsize::new(0);

fn record_references_answer(f: &[&str]) -> String {
    let name = text(f[0]);
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "shared-records-{}-{}",
        std::process::id(),
        DIRECTORIES.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&directory).expect("a directory of the test's own");
    let path = directory.join(name.trim_end_matches('/'));
    if name.ends_with('/') {
        std::fs::create_dir(&path).expect("a directory in the place of the input");
    } else if f[1] != "-" {
        std::fs::write(&path, unhex(f[1])).expect("the input file");
    }
    let (count, joined) = f[2].split_once(' ').expect("a count and its rows");
    let mut rows: Vec<Entry> = joined
        .split(';')
        .filter(|_| count != "0")
        .map(|row| {
            let parts: Vec<&str> = row.split('~').collect();
            Entry {
                key: "k".to_string(),
                location: opt_text(parts[0]),
                allele: text(parts[1]),
                reference: opt_text(parts[2]),
                ..Entry::default()
            }
        })
        .collect();
    let mut listed = rows.clone();
    let answer = |rows: &[Entry]| {
        rows.iter()
            .map(|row| row.reference.as_deref().map_or("-".to_string(), hex))
            .collect::<Vec<_>>()
            .join(";")
    };
    let shown = path.display().to_string();
    let named = |what: String| what.replace(&shown, "{path}");
    let whole = match outcome(
        || shared::read_record_references(&path, &mut rows),
        |()| String::new(),
    ) {
        ok if ok == "ok " => format!("ok {}", answer(&rows)),
        failed => {
            let (kind, what) = failed.split_once(' ').expect("a kind and a text");
            format!("{kind} {}", hex(&named(text(what))))
        }
    };
    let mut borrowed: Vec<&mut Entry> = listed.iter_mut().collect();
    let again = match outcome(
        || shared::read_record_references_of(&path, &mut borrowed),
        |()| String::new(),
    ) {
        ok if ok == "ok " => format!("ok {}", answer(&listed)),
        failed => {
            let (kind, what) = failed.split_once(' ').expect("a kind and a text");
            format!("{kind} {}", hex(&named(text(what))))
        }
    };
    std::fs::remove_dir_all(&directory).expect("the directory is removed");
    assert_eq!(whole, again, "the two forms of the reader");
    whole
}

fn answer(name: &str, f: &[&str]) -> String {
    match name {
        "_ten_digits" => plain(|| shared::_ten_digits(&text(f[0])), bit),
        "dash" => plain(|| shared::dash(&text(f[0])), eq_text),
        "allele_is_indel" => {
            let alternates = texts(f[1]);
            let alternates: Vec<&str> = alternates.iter().map(String::as_str).collect();
            plain(
                || shared::allele_is_indel(&text(f[0]), &alternates, &text(f[2]), truth(f[3])),
                opt_bit,
            )
        }
        "reference_allele_is_empty" => {
            let alternates = texts(f[1]);
            let alternates: Vec<&str> = alternates.iter().map(String::as_str).collect();
            plain(
                || shared::reference_allele_is_empty(&text(f[0]), &alternates, truth(f[2])),
                bit,
            )
        }
        "uploaded_alleles" => plain(
            || shared::uploaded_alleles(&text(f[0])),
            |alleles| match alleles {
                Some((reference, alternates)) => {
                    format!("={} {}", hex(&reference), items(alternates))
                }
                None => "-".to_string(),
            },
        ),
        "uploaded_position" => outcome(|| shared::uploaded_position(&text(f[0])), opt_whole),
        "replaced_reference" => {
            let alternates = texts(f[2]);
            let alternates: Vec<&str> = alternates.iter().map(String::as_str).collect();
            let location = opt_text(f[3]);
            outcome(
                || {
                    shared::replaced_reference(
                        whole(f[0]),
                        &text(f[1]),
                        &alternates,
                        location.as_deref(),
                        &text(f[4]),
                    )
                },
                opt_eq_text,
            )
        }
        "blank" => plain(|| shared::blank(&text(f[0])), bit),
        "integer_of" => outcome(
            || shared::integer_of(&text(f[0]), &text(f[1])),
            |number| number.to_string(),
        ),
        "integer_text_of" => outcome(
            || shared::integer_text_of(&text(f[0]), &text(f[1])),
            |number| number,
        ),
        "consequence_set" => plain(|| shared::consequence_set(&text(f[0])), eq_text),
        "feature_kind" => plain(|| shared::feature_kind(&text(f[0])), eq_text),
        "indel_row" => outcome(|| shared::indel_row(&entry(f)), bit),
        "position_bounds" => plain(
            || shared::position_bounds(&entry(&f[..8]), &text(f[8]), &text(f[9]), &text(f[10])),
            |(start, end)| bit(start) + &bit(end),
        ),
        "stop_codon_test_coordinates" => {
            plain(|| shared::stop_codon_test_coordinates(&entry(f)), bit)
        }
        "stop_codon_test_not_made" => plain(|| shared::stop_codon_test_not_made(&entry(f)), bit),
        "stop_codon_overlap_genomic_span_pair" => {
            pair_answer(shared::stop_codon_overlap_genomic_span_pair, f)
        }
        "length_changing_row" => outcome(|| shared::length_changing_row(&entry(f)), bit),
        "stop_retained_without_stop_codon_test_pair" => {
            pair_answer(shared::stop_retained_without_stop_codon_test_pair, f)
        }
        "stop_retained_three_base_stop_codon_test_pair" => {
            pair_answer(shared::stop_retained_three_base_stop_codon_test_pair, f)
        }
        "stop_codon_removed_without_stop_codon_test_pair" => {
            pair_answer(shared::stop_codon_removed_without_stop_codon_test_pair, f)
        }
        "insertion_row" => plain(|| shared::insertion_row(&entry(f)), bit),
        "inserted_between_neighbours" => {
            plain(|| shared::inserted_between_neighbours(&entry(f)), bit)
        }
        "stop_codon_overlap_insertion_range_pair" => {
            pair_answer(shared::stop_codon_overlap_insertion_range_pair, f)
        }
        "transcript_positions" => plain(
            || shared::transcript_positions(&entry(f)),
            |positions| match positions {
                Some((cdna_start, cdna_end, cds_start, cds_end)) => format!(
                    "={},{},{},{}",
                    hex(&cdna_start),
                    hex(&cdna_end),
                    hex(&cds_start),
                    hex(&cds_end)
                ),
                None => "-".to_string(),
            },
        ),
        "location_span" => outcome(|| shared::location_span(&entry(f)), opt_whole),
        "longer_allele_inside_one_exon" => {
            outcome(|| shared::longer_allele_inside_one_exon(&entry(f)), bit)
        }
        "differing_regions" => plain(
            || shared::differing_regions(&text(f[0]), &text(f[1])),
            |regions| {
                regions
                    .iter()
                    .map(|(first, last)| format!("{first}~{last}"))
                    .collect::<Vec<_>>()
                    .join(",")
            },
        ),
        "intron_rows" => {
            let regions: Vec<(i64, i64)> = f[0]
                .split(',')
                .filter(|region| !region.is_empty())
                .map(|region| {
                    let (first, last) = region.split_once('~').expect("a region");
                    (whole(first), whole(last))
                })
                .collect();
            plain(
                || {
                    shared::intron_rows(
                        &regions,
                        whole(f[1]),
                        whole(f[2]),
                        whole(f[3]),
                        truth(f[4]),
                        truth(f[5]),
                    )
                },
                |(reference, candidate)| {
                    let joined = |terms: std::collections::BTreeSet<String>| {
                        terms.into_iter().collect::<Vec<_>>().join(",")
                    };
                    format!("{}|{}", joined(reference), joined(candidate))
                },
            )
        }
        "intron_placements" => {
            let positions = [text(f[2]), text(f[3]), text(f[4]), text(f[5])];
            outcome(
                || {
                    shared::intron_placements(
                        whole(f[0]),
                        whole(f[1]),
                        (&positions[0], &positions[1], &positions[2], &positions[3]),
                    )
                },
                placements,
            )
        }
        "record_span_placement" => {
            let (e, a) = (entry(&f[..8]), entry(&f[8..16]));
            outcome(
                || shared::record_span_placement(&e, &a, &text(f[16])),
                |place| match place {
                    Some((forward, start, end, wide_exons)) => {
                        format!("{}~{start}~{end}~{}", bit(forward), bit(wide_exons))
                    }
                    None => "-".to_string(),
                },
            )
        }
        "splice_term_past_the_record_span_pair" => {
            pair_answer(shared::splice_term_past_the_record_span_pair, f)
        }
        "vep_rs_pos" => plain(|| shared::vep_rs_pos(&text(f[0])), opt_whole),
        "_strip_chr" => plain(|| shared::_strip_chr(&text(f[0])), eq_text),
        "perl_fields" => plain(|| shared::perl_fields(&texts(f[0])), items),
        "perl_split" => plain(|| shared::perl_split(&text(f[0]), &text(f[1])), items),
        "sequence_location" => plain(
            || shared::sequence_location(&text(f[0]), whole(f[1]), whole(f[2])),
            eq_text,
        ),
        "canonical_location" => plain(|| shared::canonical_location(&text(f[0])), eq_text),
        "_canonical_coordinate" => plain(|| shared::_canonical_coordinate(&text(f[0])), eq_text),
        "gvcf_site_rows" => plain(
            || shared::gvcf_site_rows(whole(f[0]), &text(f[1]), &texts(f[2])),
            |rows| {
                rows.iter()
                    .map(|(allele, start, end)| format!("{}~{start}~{end}", hex(allele)))
                    .collect::<Vec<_>>()
                    .join(";")
            },
        ),
        "is_structural_allele" => plain(|| shared::is_structural_allele(&text(f[0])), bit),
        "_is_start_cooccurrence_swap" => plain(
            || shared::_is_start_cooccurrence_swap(&text(f[0]), &text(f[1]), &text(f[2])),
            bit,
        ),
        "_is_covered_splice_region_swap" => plain(
            || shared::_is_covered_splice_region_swap(&text(f[0]), &text(f[1])),
            bit,
        ),
        "EXCLUDING_RULES" => plain(
            || {
                let (perl, rust, allele) = (text(f[0]), text(f[1]), text(f[2]));
                shared::EXCLUDING_RULES
                    .iter()
                    .map(|rule| {
                        format!(
                            "{}~{}~{}",
                            rule.bucket,
                            rule.taxonomy_class,
                            bit(rule.matches(&perl, &rust, &allele))
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(",")
            },
            |answer| answer,
        ),
        "read_record_references" => record_references_answer(f),
        other => panic!("no function for the file {other}"),
    }
}

/// The cases of every file under `dir`, and the ones answered otherwise than the file says.
fn disagreements(dir: &Path) -> (usize, Vec<String>) {
    quiet_uncaught();
    let (mut cases, mut differing) = (0, Vec::new());
    for name in FILES {
        let path = dir.join(format!("{name}.tsv"));
        let held = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let before = cases;
        for (number, line) in held.lines().enumerate() {
            let mut fields: Vec<&str> = line.split('\t').collect();
            let expected = fields.pop().expect("an answer");
            let given = panic::catch_unwind(AssertUnwindSafe(|| answer(name, &fields)))
                .unwrap_or_else(|_| "the case could not be read or answered".to_string());
            cases += 1;
            if given != expected {
                let cut = |text: &str| text.chars().take(160).collect::<String>();
                differing.push(format!(
                    "{name}.tsv line {}: the file says {}, the tool {}",
                    number + 1,
                    cut(expected),
                    cut(&given)
                ));
            }
        }
        assert!(cases > before, "{name}.tsv holds no case");
    }
    (cases, differing)
}

fn agree(dir: &Path, least: usize) {
    let (cases, differing) = disagreements(dir);
    assert!(
        differing.is_empty(),
        "{} of {cases} cases differ:\n{}",
        differing.len(),
        differing
            .iter()
            .take(40)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(cases >= least, "{cases} cases read, {least} expected");
}

#[test]
fn a_cell_that_is_not_pairs_of_lowercase_hex_digits_is_refused() {
    assert_eq!(unhex(""), [0u8; 0]);
    assert_eq!(unhex("00_7f_ff"), [0x00, 0x7f, 0xff]);
    for damaged in [
        "6", "616", "61_", "_61", "61__62", "61 62", "6g", "6A", "+a",
    ] {
        assert!(panic::catch_unwind(|| unhex(damaged)).is_err(), "{damaged}");
    }
}

#[test]
fn the_committed_cases_are_answered_as_their_files_say() {
    agree(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/shared"),
        5_000,
    );
}

#[test]
#[ignore = "reads the case files under the directory VEP_FORMATS_COMPARE_SHAREDDIFF_DIR names"]
fn every_case_of_the_whole_files_is_answered_as_its_file_says() {
    let dir = PathBuf::from(
        std::env::var_os("VEP_FORMATS_COMPARE_SHAREDDIFF_DIR")
            .expect("VEP_FORMATS_COMPARE_SHAREDDIFF_DIR"),
    );
    agree(&dir, 300_000);
}
