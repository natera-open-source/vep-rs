// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a plugin value at a Location of a record the reference dropped, asked of
//! `dropped_by_the_reference` over the rows the two sides wrote at one Location, the
//! reference's stored output and the input's record lines.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use vep_plugin_pairs::rules::skipped_record::dropped_by_the_reference;
use vep_plugin_pairs::rules::Refusal;

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-skipped-record-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn file(&self, name: &str, content: impl AsRef<[u8]>) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, content).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn row(uploaded: &str, location: &str, allele: &str, feature: &str, extra: &str) -> String {
    format!("{uploaded}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\t{extra}\n")
}

fn gzip(content: impl AsRef<[u8]>) -> Vec<u8> {
    let mut packed = GzEncoder::new(Vec::new(), Compression::default());
    packed.write_all(content.as_ref()).unwrap();
    packed.finish().unwrap()
}

fn none() -> Vec<String> {
    Vec::new()
}

/// Rows as the settling step hands them: without their line ends.
fn handed(rows: &[String]) -> Vec<String> {
    rows.iter()
        .map(|row| row.trim_end_matches('\n').to_string())
        .collect()
}

#[test]
fn the_reference_is_read_for_the_names_at_a_location_the_engine_alone_wrote_a_structural_row_at() {
    let scratch = Scratch::new("names");
    let dropped = vec![row(
        "sv_over",
        "21:7000001-27000000",
        "deletion",
        "ENST1",
        "-",
    )];
    let sequence = vec![row("snv", "21:5", "T", "ENST1", "-")];
    let stored = scratch.file(
        "reference.txt.gz",
        gzip(format!(
            "{HEADER}{}",
            row("sv_else", "21:5", "T", "ENST1", "-")
        )),
    );
    assert_eq!(
        dropped_by_the_reference(&none(), &dropped, &stored, &none()),
        Ok(true)
    );
    assert_eq!(
        dropped_by_the_reference(&none(), &sequence, &stored, &none()),
        Ok(false)
    );
    assert_eq!(
        dropped_by_the_reference(&dropped, &dropped, &stored, &none()),
        Ok(false)
    );

    // Of the three, the first alone opens the output: a path that names no file ends it
    // and leaves the other two as they are.
    let missing = scratch.0.join("no-such-output.txt");
    let refused = dropped_by_the_reference(&none(), &dropped, &missing, &none()).unwrap_err();
    assert!(
        refused.0.starts_with(&format!("{}: ", missing.display())),
        "{}",
        refused.0
    );
    assert_eq!(
        dropped_by_the_reference(&none(), &sequence, &missing, &none()),
        Ok(false)
    );
    assert_eq!(
        dropped_by_the_reference(&dropped, &dropped, &missing, &none()),
        Ok(false)
    );

    // A name the reference carries on a row of a gzip output gives the Location back.
    let stored = scratch.file(
        "reference.txt.gz",
        gzip(format!(
            "{HEADER}{}",
            row("sv_over", "21:9000-9500", "deletion", "ENST1", "-")
        )),
    );
    assert_eq!(
        dropped_by_the_reference(&none(), &dropped, &stored, &none()),
        Ok(false)
    );
}

/// One cell of the grid: its name, the reference's rows and the engine's rows at its
/// Location, the input's record lines that carry the first column of an engine row there,
/// and whether the rows are a dropped record's with those lines and without them.
struct Cell {
    name: &'static str,
    reference: Vec<String>,
    engine: Vec<String>,
    records: Vec<&'static str>,
    with_records: bool,
    without: bool,
}

fn grid() -> Vec<Cell> {
    let cell = |name, reference, engine, records, with_records, without| Cell {
        name,
        reference,
        engine,
        records,
        with_records,
        without,
    };
    vec![
        cell(
            "a structural record over the size limit",
            vec![],
            vec![
                row(
                    "sv_over",
                    "21:7000001-27000000",
                    "deletion",
                    "ENST1",
                    "pHaplo=0.9;pTriplo=0.2",
                ),
                row(
                    "sv_over",
                    "21:7000001-27000000",
                    "deletion",
                    "ENST2",
                    "pHaplo=0.5;Existing_uORFs=2",
                ),
            ],
            vec!["21\t7000000\tsv_over\tN\t<DEL>"],
            true,
            true,
        ),
        cell(
            "a structural record under the limit",
            vec![],
            vec![row(
                "sv_under",
                "21:30000001-30000900",
                "duplication",
                "ENST1",
                "pHaplo=0.7",
            )],
            vec!["21\t30000000\tsv_under\tN\t<DUP>"],
            true,
            true,
        ),
        cell(
            "a sequence allele beside the symbolic one",
            vec![],
            vec![
                row(
                    "sv_mixed",
                    "21:37000001-37000002",
                    "T",
                    "ENST1",
                    "pHaplo=0.12",
                ),
                row(
                    "sv_mixed",
                    "21:37000001-37000002",
                    "deletion",
                    "ENST1",
                    "pHaplo=0.11",
                ),
            ],
            vec!["21\t37000000\tsv_mixed\tN\t<DEL>,NT"],
            true,
            true,
        ),
        cell(
            "every row a sequence allele",
            vec![],
            vec![row("snv_only", "21:31000000", "T", "ENST1", "pHaplo=0.6")],
            vec![],
            false,
            false,
        ),
        cell(
            "the allele tandem_repeat",
            vec![],
            vec![row(
                "tr_1",
                "21:32000001-32000050",
                "tandem_repeat",
                "ENST1",
                "pHaplo=0.4",
            )],
            vec![],
            false,
            false,
        ),
        cell(
            "the reference wrote a row at the Location",
            vec![row(
                "sv_kept",
                "21:33000001-33000500",
                "deletion",
                "ENST1",
                "IMPACT=MODIFIER",
            )],
            vec![row(
                "sv_kept",
                "21:33000001-33000500",
                "deletion",
                "ENST1",
                "pHaplo=0.3",
            )],
            vec![],
            false,
            false,
        ),
        cell(
            "a reference row carries the record's name elsewhere",
            vec![],
            vec![row(
                "sv_moved",
                "21:34000000-34000500",
                "deletion",
                "ENST1",
                "pHaplo=0.8",
            )],
            vec![],
            false,
            false,
        ),
        cell(
            "an engine row at the Location names no record",
            vec![],
            vec![
                row(
                    "sv_named",
                    "21:35000001-35900000",
                    "deletion",
                    "ENST1",
                    "pHaplo=0.2",
                ),
                row(
                    ".",
                    "21:35000001-35900000",
                    "deletion",
                    "ENST2",
                    "pHaplo=0.25",
                ),
            ],
            vec![],
            false,
            false,
        ),
        cell(
            "a field the class does not read, on a dropped record's row all the same",
            vec![],
            vec![row(
                "sv_field",
                "21:36000001-36900000",
                "deletion",
                "ENST1",
                "5UTR_consequence=gained",
            )],
            vec![],
            true,
            true,
        ),
        cell(
            "a gVCF variant site under an identifier",
            vec![],
            vec![
                row("gv_1", "21:38000000", "<NON_REF>", "ENST1", "pHaplo=0.14"),
                row("gv_1", "21:38000000", "T", "ENST1", "pHaplo=0.13"),
            ],
            vec!["21\t38000000\tgv_1\tA\tT,<NON_REF>"],
            false,
            true,
        ),
        cell(
            "a gVCF variant site under a generated name",
            vec![],
            vec![
                row(
                    "21_39000000_A/T/<NON_REF>",
                    "21:39000000",
                    "<NON_REF>",
                    "ENST1",
                    "pHaplo=0.16",
                ),
                row(
                    "21_39000000_A/T/<NON_REF>",
                    "21:39000000",
                    "T",
                    "ENST1",
                    "pHaplo=0.15",
                ),
            ],
            vec![],
            false,
            false,
        ),
    ]
}

#[test]
fn the_grid_of_dropped_and_kept_records_reads_alike_from_a_plain_and_a_gzip_reference() {
    let scratch = Scratch::new("grid");
    let output = format!(
        "{HEADER}{}{}{}{}",
        row("snv_1", "21:100", "T", "ENST9", "pHaplo=0.01"),
        row(
            "sv_ref",
            "21:40000001-40000500",
            "deletion",
            "ENST1",
            "pHaplo=0.1"
        ),
        row(
            "sv_kept",
            "21:33000001-33000500",
            "deletion",
            "ENST1",
            "IMPACT=MODIFIER"
        ),
        row(
            "sv_moved",
            "21:34000001-34000500",
            "deletion",
            "ENST1",
            "IMPACT=MODIFIER"
        ),
    );
    let plain = scratch.file("reference.txt", &output);
    let packed = scratch.file("reference.txt.gz", gzip(&output));
    for stored in [&plain, &packed] {
        for cell in grid() {
            let (reference, engine) = (handed(&cell.reference), handed(&cell.engine));
            let records: Vec<String> = cell.records.iter().map(|line| line.to_string()).collect();
            assert_eq!(
                dropped_by_the_reference(&reference, &engine, stored, &records),
                Ok(cell.with_records),
                "{} with the input's records, {}",
                cell.name,
                stored.display()
            );
            assert_eq!(
                dropped_by_the_reference(&reference, &engine, stored, &none()),
                Ok(cell.without),
                "{} without them, {}",
                cell.name,
                stored.display()
            );
        }
    }
}

#[test]
fn a_stored_reference_is_read_as_text_whatever_its_line_ends_and_its_packing() {
    let scratch = Scratch::new("stored");
    let dropped = vec![row("sv_a", "21:100-200", "deletion", "ENST1", "-")];
    let named = row("sv_a", "21:900-950", "deletion", "ENST1", "-");
    let other = row("rs9", "21:5", "T", "ENST1", "-");
    let asked = |name: &str, content: Vec<u8>, engine: &[String]| {
        dropped_by_the_reference(&none(), engine, &scratch.file(name, content), &none())
    };

    let old_line_ends = format!("{HEADER}{other}{named}{other}").replace('\n', "\r");
    assert_eq!(
        asked("cr.txt", old_line_ends.into_bytes(), &dropped),
        Ok(false)
    );
    let no_last_end = format!("{HEADER}{other}{}", named.trim_end());
    assert_eq!(
        asked("no-end.txt", no_last_end.into_bytes(), &dropped),
        Ok(false)
    );

    // Two members with zero bytes after each, the name split across them.
    let mut members = gzip(format!("{HEADER}sv_"));
    members.extend([0u8; 7]);
    members.extend(gzip(&named[3..]));
    members.extend([0u8; 512]);
    assert_eq!(asked("members.txt.gz", members, &dropped), Ok(false));

    // A row of six columns, a commented row and a row without a Feature name no record.
    let six_columns = "sv_a\t21:900-950\tdeletion\tENSG1\tENST1\tTranscript\n";
    let no_feature = row("sv_a", "21:900-950", "deletion", " ", "-");
    let unnamed = format!("{HEADER}{six_columns}#{named}{no_feature}");
    assert_eq!(
        asked("unnamed.txt", unnamed.into_bytes(), &dropped),
        Ok(true)
    );
    assert_eq!(asked("empty.txt", Vec::new(), &dropped), Ok(true));

    // The name as written between blanks, and a malformed byte read as U+FFFD.
    let blanks = format!(
        "{HEADER}{}",
        row(" \u{1c}sv_a \u{85}", "21:900-950", "deletion", "ENST1", "-")
    );
    assert_eq!(
        asked("blanks.txt", blanks.into_bytes(), &dropped),
        Ok(false)
    );
    let mut malformed = HEADER.as_bytes().to_vec();
    malformed.extend(b"sv\xff\t21:900\tdeletion\tENSG1\tENST1\tTranscript\tintron_variant\n");
    let replaced = vec![row("sv\u{fffd}", "21:100-200", "deletion", "ENST1", "-")];
    assert_eq!(asked("malformed.txt", malformed, &replaced), Ok(false));

    // Two records at one Location stand or fall together, and two Locations are not one.
    let together = vec![
        dropped[0].clone(),
        row("sv_b", "21:100-200", "T", "ENST1", "-"),
    ];
    let names_one = format!("{HEADER}{named}").into_bytes();
    assert_eq!(
        asked("names-one.txt", names_one.clone(), &together),
        Ok(false)
    );
    let apart = vec![
        row("sv_b", "21:100-200", "deletion", "ENST1", "-"),
        row("sv_a", "21:300-400", "duplication", "ENST1", "-"),
    ];
    assert_eq!(asked("names-one.txt", names_one, &apart), Ok(false));
    assert_eq!(
        asked("header.txt", HEADER.as_bytes().to_vec(), &apart),
        Ok(true)
    );
}

#[test]
fn a_row_without_its_columns_and_a_number_past_the_digits_of_an_integer_end_the_reading() {
    let scratch = Scratch::new("ends");
    let stored = scratch.file("reference.txt", HEADER);
    let asked = |reference: &[String], engine: &[String], records: &[String]| {
        dropped_by_the_reference(reference, engine, &stored, records).map_err(|refusal| refusal.0)
    };
    let short = "a row of fewer than three columns states no Location and no allele".to_string();
    assert_eq!(
        asked(&none(), &["sv_c\t21:100".to_string()], &none()),
        Err(short.clone())
    );
    assert_eq!(asked(&["sv_c".to_string()], &none(), &none()), Err(short));
    assert_eq!(
        asked(&none(), &["sv_c\t21:100\tdeletion".to_string()], &none()),
        Ok(true)
    );

    let long = "7".repeat(4301);
    let digits = |what: &str| {
        Err(format!(
            "{what} of 4301 digits is longer than the 4300 digits an integer is read from"
        ))
    };
    let at = |location: &str| vec![row("sv_c", location, "deletion", "ENST1", "-")];
    assert_eq!(
        asked(&none(), &at(&format!("21:5-{long}")), &none()),
        digits("a Location coordinate")
    );
    assert_eq!(
        asked(&none(), &at(&format!("21:5-{}", &long[1..])), &none()),
        Ok(true)
    );
    assert_eq!(
        asked(&none(), &at(&format!("21:x{long}")), &none()),
        Ok(true)
    );
    let mate = vec![row(
        "sv_c",
        "21:100",
        &format!("N[chr5:{long}["),
        "ENST1",
        "-",
    )];
    assert_eq!(
        asked(&none(), &mate, &none()),
        digits("a breakend mate position")
    );
    let record = vec![format!("21\t{long}\t.\tA\tG")];
    assert_eq!(asked(&none(), &at("21:100"), &record), digits("a POS"));
    let blank_before = vec![format!("21\t {long}\t.\tA\tG,<NON_REF>")];
    assert_eq!(asked(&none(), &at("21:100"), &blank_before), Ok(true));

    let unread = dropped_by_the_reference(&none(), &at("21:100"), &scratch.0, &none());
    assert!(
        matches!(unread, Err(Refusal(text)) if text.starts_with(&format!("{}: ", scratch.0.display())))
    );
    let cut = scratch.file(
        "cut.txt.gz",
        &gzip(format!("{HEADER}{}", at("21:5")[0].repeat(400)))[..60],
    );
    let broken = dropped_by_the_reference(&none(), &at("21:100"), &cut, &none());
    assert!(
        matches!(broken, Err(Refusal(text)) if text.starts_with(&format!("{}: cannot read the gzip stream: ", cut.display())))
    );
}
