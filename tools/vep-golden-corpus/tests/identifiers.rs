// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Identifiers under restricted redistribution terms: the two shapes read in a text and the
//! fabricated form written in their place, `select` and `check` on inputs and corpora that
//! carry them, the `identifiers` subcommand on stored files, and the committed corpora. Every
//! identifier here is made up, and is written as its letters and its digits apart.

mod common;

use common::{
    clean, committed_manifests, failed, manifest, records, row, run, scratch, select,
    two_record_dataset, under, vcf_of, write, Run, VCF_HEADER,
};
use flate2::{Compression, GzBuilder};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use vep_golden_corpus::identifiers::stored::{
    corpus_files, read_corpus_text, restricted_identifier_errors,
};
use vep_golden_corpus::identifiers::{
    fabricated_identifiers, is_fabricated, restricted_identifier_count, restricted_identifiers,
    RESTRICTED_IDENTIFIER,
};

/// What a line that refuses a file says after the file and its count.
const RESTRICTED_REFUSAL: &str = "under restricted redistribution terms in another form than \
                                  the fabricated one; the identifiers subcommand replaces them";

/// A made-up identifier: its letters, then its digits.
fn made_up(letters: &str, digits: &str) -> String {
    format!("{letters}{digits}")
}

/// A made-up identifier of the first shape, its digits a counting run.
fn first() -> String {
    made_up("COSV", "12345678")
}

/// A made-up identifier of the second shape, its digits a descending run.
fn second() -> String {
    made_up("CM", "654321")
}

/// The form the tool writes in the place of `first` and of `second`.
fn zeroed_first() -> String {
    made_up("COSV", "00000000")
}

fn zeroed_second() -> String {
    made_up("CM", "000000")
}

/// `template` with `{first}` and `{second}` read as the two made-up identifiers, `{FIRST}`
/// and `{SECOND}` as their fabricated forms, and each `[LETTERS:DIGITS]` as the made-up
/// identifier of those letters and digits.
fn text(template: &str) -> String {
    let named = template
        .replace("{first}", &first())
        .replace("{second}", &second())
        .replace("{FIRST}", &zeroed_first())
        .replace("{SECOND}", &zeroed_second());
    let mut filled = String::new();
    let mut rest = named.as_str();
    while let Some(open) = rest.find('[') {
        let (before, after) = rest.split_at(open);
        filled.push_str(before);
        let inner = after[1..].split_once(']').map(|(inner, _)| inner);
        let pair = inner.and_then(|inner| inner.split_once(':'));
        match pair {
            Some((letters, digits))
                if !letters.is_empty()
                    && letters.bytes().all(|byte| byte.is_ascii_uppercase())
                    && !digits.is_empty()
                    && digits.bytes().all(|byte| byte.is_ascii_digit()) =>
            {
                filled.push_str(&made_up(letters, digits));
                rest = &after[letters.len() + digits.len() + 3..];
            }
            _ => {
                filled.push('[');
                rest = &after[1..];
            }
        }
    }
    filled + rest
}

fn none() -> BTreeSet<String> {
    BTreeSet::new()
}

fn fabricated(text: &str) -> String {
    fabricated_identifiers(text, &none()).into_owned()
}

fn count(text: &str) -> usize {
    restricted_identifier_count(text, &none())
}

fn sha(data: &[u8]) -> String {
    let digest = Sha256::digest(data);
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// A gzip file of one member that holds `data`, under `name` where one is given, with the
/// time stamp `mtime`.
fn gzip(data: &[u8], name: Option<&str>, mtime: u32) -> Vec<u8> {
    let builder = GzBuilder::new().mtime(mtime);
    let builder = match name {
        Some(name) => builder.filename(name),
        None => builder,
    };
    let mut encoder = builder.write(Vec::new(), Compression::best());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn gunzip(stored: &[u8]) -> Vec<u8> {
    let mut data = Vec::new();
    flate2::read::MultiGzDecoder::new(stored)
        .read_to_end(&mut data)
        .unwrap();
    data
}

/// One BGZF block that holds `data`: a gzip member whose extra field carries the subfield
/// `BC` with the block's size.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut deflate = flate2::write::DeflateEncoder::new(Vec::new(), Compression::new(6));
    deflate.write_all(data).unwrap();
    let body = deflate.finish().unwrap();
    let mut checksum = flate2::Crc::new();
    checksum.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend(u16::try_from(body.len() + 25).unwrap().to_le_bytes());
    block.extend(body);
    block.extend(checksum.sum().to_le_bytes());
    block.extend(u32::try_from(data.len()).unwrap().to_le_bytes());
    block
}

fn read(path: &str) -> String {
    fs::read_to_string(path).unwrap()
}

#[test]
fn fabricated_identifiers_replace_both_shapes_and_leave_public_names() {
    // The first shape becomes its letters and a zero for each digit wherever it stands as a
    // word, the second (three of its type letters here) where it is a whole member of a list
    // joined by `&`, `,`, `|` or `;` or a JSON `id`, and from there wherever the same
    // identifier stands as a word in the text or is named to the call.
    let replaced = [
        (
            "vep=G|{first}|{second}&{second}&rs1;AF=0.1",
            "vep=G|{FIRST}|{SECOND}&{SECOND}&rs1;AF=0.1",
        ),
        (
            "21\t5\trs1;{second}\tA\tG\t{first}",
            "21\t5\trs1;{SECOND}\tA\tG\t{FIRST}",
        ),
        (
            "{second},rs2\trs3,{second}\nrs4|{second}\n",
            "{SECOND},rs2\trs3,{SECOND}\nrs4|{SECOND}\n",
        ),
        (
            r#"{"id":"{second}","input":"21\t5\t{first}\t{second}&rs1\tA"}"#,
            r#"{"id":"{SECOND}","input":"21\t5\t{FIRST}\t{SECOND}&rs1\tA"}"#,
        ),
        (
            "DB={second},rs1;see {first}.",
            "DB={SECOND},rs1;see {FIRST}.",
        ),
        ("[COSM:1],[COSM:22]", "[COSM:0],[COSM:00]"),
        ("|[CD:123456]&[CI:123456]|", "|[CD:000000]&[CI:000000]|"),
        (
            "[CD:123456],rs1\t[CI:123456]|rs2",
            "[CD:000000],rs1\t[CI:000000]|rs2",
        ),
        (
            r#"{"id":"[CD:123456]"} {"id": "[CI:123456]"}"#,
            r#"{"id":"[CD:000000]"} {"id": "[CI:000000]"}"#,
        ),
        (r#"{"id": "{second}"}"#, r#"{"id": "{SECOND}"}"#),
        ("rs1&{second}\r\n", "rs1&{SECOND}\r\n"),
        (
            "rs1,{second}\n{second}\t1\n\"{second}\"",
            "rs1,{SECOND}\n{SECOND}\t1\n\"{SECOND}\"",
        ),
    ];
    for (given, expected) in replaced {
        let (given, expected) = (text(given), text(expected));
        assert_eq!(fabricated(&given), expected);
        assert_eq!(
            fabricated(&expected),
            expected,
            "the fabricated form is left as it is"
        );
        assert!(count(&given) > 0 && count(&expected) == 0, "{expected}");
    }
    // The second shape alone in a table cell, opening a field before `;`, with a version
    // after a dot, after a colon, in a run of eight digits or in prose keeps its text.
    let kept = [
        "1\t{second}.1\tchr1",
        "|{second}.2|",
        "CLNDISDB=MedGen:{second},OMIM:1|MedGen:{second}",
        "|{second}78|",
        "|A{second}|",
        "Name={second};rs1",
        "|{second}A|",
        "see {second} here",
        "|{first}X|",
        "|X{first}|",
        "{first}_1",
        "COSV,COSM,CM12345",
        "{second}\t1\n1\t{second}\n",
        r#""synonyms":["{second}","chr1"]"#,
        r#""input":"21\t5\t{second}\tA""#,
    ];
    for (place, given) in kept.into_iter().enumerate() {
        let given = text(given);
        assert_eq!(fabricated(&given), given, "kept case {place}");
        assert_eq!(count(&given), 0, "kept case {place}");
    }
    assert!(is_fabricated(&zeroed_first()) && is_fabricated(&zeroed_second()));
    assert!(!is_fabricated(&first()));
    // An identifier read in a list elsewhere is named to the call, and is then read alone
    // in a cell.
    let named = BTreeSet::from([second()]);
    assert_eq!(
        fabricated_identifiers(&text("{second}\t1\n\"{second}\""), &named),
        text("{SECOND}\t1\n\"{SECOND}\"")
    );
    assert_eq!(
        restricted_identifier_count(&text("{second}\t1\n"), &named),
        1
    );
    for given in ["{second}.1\t1", "MedGen:{second}", "|{second}78|"] {
        let given = text(given);
        assert_eq!(fabricated_identifiers(&given, &named), given);
    }
}

#[test]
fn select_writes_no_restricted_identifier() {
    // A record line whose INFO or ID list carries an identifier of either shape is kept with
    // the fabricated one in its place, in `variants.vcf` and in the manifest alike, and a
    // name after a colon stays. A record named by such an identifier is passed over, by the
    // multi-allelic sample too, when no reference row names it and refused, with the input
    // and a count, when one does: nothing is written.
    let scratch = scratch("identifiers-select");
    let rows = row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + &row("rs2", "21:200", "T", "ENST2", "intron_variant");
    let reference = write(&scratch.path, "ref.txt", &rows);
    let lines = text(concat!(
        "21\t100\trs1\tA\tG\t.\tPASS\tvep=G|{first}|{second}&{second}&rs1;",
        "CLNDISDB=MedGen:{second}\n",
        "21\t200\trs2;{second}\tC\tT\t.\tPASS\t.\n",
        "21\t300\t{first}\tC\tT,G\t.\tPASS\t.\n",
    ));
    let input = write(&scratch.path, "in.vcf", &format!("{VCF_HEADER}{lines}"));
    let out = under(&scratch.path, "corpus");
    let options = ["--covers", "symbol", "--multi-allelic-per-dataset", "1"];
    assert_eq!(select(&reference, &input, &out, &options).code, 0);
    let kept = [
        text(concat!(
            "21\t100\trs1\tA\tG\t.\tPASS\tvep=G|{FIRST}|{SECOND}&{SECOND}&rs1;",
            "CLNDISDB=MedGen:{second}",
        )),
        text("21\t200\trs2;{SECOND}\tC\tT\t.\tPASS\t."),
    ];
    let variants = read(&under(Path::new(&out), "variants.vcf"));
    assert_eq!(variants.lines().skip(2).collect::<Vec<_>>(), kept);
    let written = manifest(&out);
    let lines: Vec<&str> = records(&written).iter().map(vcf_of).collect();
    assert_eq!(lines, kept);
    let unresolved = written.get("unresolved_exemplars").unwrap();
    assert_eq!(unresolved.as_list().map(<[_]>::len), Some(0));
    for name in ["variants.vcf", "manifest.json"] {
        let held = read(&under(Path::new(&out), name));
        assert!(!held.contains(&first()) && count(&held) == 0, "{name}");
    }
    let dataset = format!("s01={reference}={input}");
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &dataset]),
        clean(0, "observed 2, covered by corpus 2, missing 0\n")
    );
    let named = row(&first(), "21:300", "T", "ENST2", "intron_variant");
    let named = write(&scratch.path, "named.txt", &named);
    let refused = under(&scratch.path, "refused");
    assert_eq!(
        select(&named, &input, &refused, &["--covers", "symbol"]),
        failed(&format!(
            "{input} names 1 record the references chose by an identifier under restricted \
             redistribution terms in the ID column; a corpus carries no such record, since its \
             reference rows are named by the identifier"
        ))
    );
    assert!(!Path::new(&refused).exists());
}

#[test]
fn check_refuses_a_restricted_identifier_without_printing_it() {
    // `check` reads every text file of the corpus, a gzip one inflated, and refuses the
    // corpus on an identifier of either shape that is not the fabricated form with one line a
    // file, giving the file and the count and never the identifier, before it reads a
    // reference; a second-shape identifier one file lists is counted alone in another file's
    // cell, and a file whose one identifier has the type letter `D` is read.
    let scratch = scratch("identifiers-check");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]).code,
        0
    );
    let expected = Path::new(&out).join("expected");
    fs::create_dir(&expected).unwrap();
    let stored = expected.join("vcf.vcf.gz");
    let line =
        |a: &str, b: &str, c: &str| format!("21\t100\trs1\tA\tG\t.\tPASS\tCSQ=G|{a}&{b}|{c}\n");
    let fabricated_line = line(&zeroed_second(), "rs1", &zeroed_first());
    fs::write(&stored, gzip(fabricated_line.as_bytes(), None, 0)).unwrap();
    let index = [b"\0".as_slice(), first().as_bytes()].concat();
    fs::write(expected.join("vcf.vcf.gz.csi"), index).unwrap();
    let dataset = format!("s01={reference}={input}");
    // The fabricated form passes, and a file that is not text is not read.
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &dataset]),
        clean(0, "observed 2, covered by corpus 2, missing 0\n")
    );
    let listed = line(&second(), &second(), &first());
    fs::write(&stored, gzip(listed.as_bytes(), None, 0)).unwrap();
    write(Path::new(&out), "notes.txt", &text("x|{first}|y\n"));
    write(
        &expected,
        "tab.txt",
        &text("rs1\t{second}\n1\t[CI:654321]\n"),
    );
    write(&expected, "default.txt", &text("rs1,[CD:123456]\n"));
    // Refused before any reference or input is opened: the identifier the stored output
    // lists is read alone in the tab cell, and one no file lists is not.
    let absent = format!("s01={reference}={}", under(&scratch.path, "absent.vcf"));
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &absent]),
        Run {
            code: 1,
            stdout: String::new(),
            stderr: format!(
                "ERROR: [build_golden_corpus] expected/default.txt holds 1 identifier {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] expected/tab.txt holds 1 identifier {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] expected/vcf.vcf.gz holds 3 identifiers {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] notes.txt holds 1 identifier {RESTRICTED_REFUSAL}\n"
            ),
        }
    );
}

/// The text of a provenance record of three stored files, each under the digest given for
/// it: the two expected outputs uncompressed, the input, and the two outputs as stored.
fn provenance(digests: [&str; 5]) -> String {
    let [json, tab, variants, json_stored, tab_stored] = digests;
    format!(
        "{{\n \"sha256\": {{\n  \"json.jsonl\": \"{json}\",\n  \"tab.txt\": \"{tab}\"\n }},\n \
         \"variants_vcf_sha256\": \"{variants}\",\n \"stored\": {{\n  \
         \"json.jsonl.gz\": \"{json_stored}\",\n  \"tab.txt.gz\": \"{tab_stored}\"\n }}\n}}\n"
    )
}

#[test]
fn identifiers_rewrites_stored_files_and_records_their_digests_again() {
    // `identifiers <dir>` replaces both shapes in every text file under the directory, writes
    // a gzip file back with no name and no time stamp, and records again each digest the
    // directory's `provenance.json` holds of a rewritten file, of its bytes or of its
    // uncompressed bytes; `--check` writes nothing and fails while one is left; a
    // block-compressed file that holds one is refused before anything is written, under any
    // path the command names.
    let scratch = scratch("identifiers-rewrite");
    let corpus = under(&scratch.path, "corpus");
    let expected = Path::new(&corpus).join("expected");
    fs::create_dir_all(&expected).unwrap();
    let variants_text = text(concat!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "21\t100\trs1\tA\tG\t.\tPASS\tvep=G|{second}&rs1\n",
    ));
    let variants = write(Path::new(&corpus), "variants.vcf", &variants_text);
    let output = text("{\"input\":\"21\\t100\\t{first}\\tA\\tG\",\"id\":\"{second}\"}\n");
    let stored = under(&expected, "json.jsonl.gz");
    fs::write(&stored, gzip(output.as_bytes(), Some("json.jsonl"), 1)).unwrap();
    let alone = write(&expected, "default.txt", &text("rs1\t{second}\n"));
    let untouched = under(&expected, "tab.txt.gz");
    let tab_text = text("rs1,{SECOND}\n");
    let tab_bytes = gzip(tab_text.as_bytes(), None, 1);
    fs::write(&untouched, &tab_bytes).unwrap();
    let record = write(
        Path::new(&corpus),
        "provenance.json",
        &provenance([
            &sha(output.as_bytes()),
            &sha(tab_text.as_bytes()),
            &sha(variants_text.as_bytes()),
            &sha(&fs::read(&stored).unwrap()),
            &sha(&tab_bytes),
        ]),
    );

    assert_eq!(
        run(&["identifiers", "--check", &corpus]),
        Run {
            code: 1,
            stdout: format!(
                "{corpus}: 5 files read, 3 hold an identifier in another form than the fabricated one\n"
            ),
            stderr: format!(
                "ERROR: [build_golden_corpus] {alone} holds 1 identifier {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] {stored} holds 2 identifiers {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] {variants} holds 1 identifier {RESTRICTED_REFUSAL}\n"
            ),
        }
    );
    assert_eq!(
        gunzip(&fs::read(&stored).unwrap()),
        output.as_bytes(),
        "--check writes nothing"
    );

    assert_eq!(
        run(&["identifiers", &corpus]),
        clean(
            0,
            &format!("{corpus}: 4 identifiers replaced in 3 files, 3 digests recorded again\n")
        )
    );
    // Read alone in a cell once another file lists it.
    assert_eq!(read(&alone), text("rs1\t{SECOND}\n"));
    let rewritten = output
        .replace(&first(), &zeroed_first())
        .replace(&second(), &zeroed_second());
    let stored_bytes = fs::read(&stored).unwrap();
    assert_eq!(gunzip(&stored_bytes), rewritten.as_bytes());
    assert_eq!(
        stored_bytes[..10],
        [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 2, 3],
        "no name flag and no time stamp, in the header zlib writes at level 9"
    );
    let variants_after = read(&variants);
    assert!(variants_after.ends_with(&text("vep=G|{SECOND}&rs1\n")));
    assert_eq!(
        fs::read(&untouched).unwrap(),
        tab_bytes,
        "a file that holds only the fabricated form is not rewritten"
    );
    assert_eq!(
        read(&record),
        provenance([
            &sha(rewritten.as_bytes()),
            &sha(tab_text.as_bytes()),
            &sha(variants_after.as_bytes()),
            &sha(&stored_bytes),
            &sha(&tab_bytes),
        ])
    );
    assert_eq!(
        run(&["identifiers", "--check", &corpus]),
        clean(
            0,
            &format!(
                "{corpus}: 5 files read, 0 hold an identifier in another form than the fabricated one\n"
            )
        )
    );
    assert_eq!(
        run(&["identifiers", &corpus]),
        clean(
            0,
            &format!("{corpus}: 0 identifiers replaced in 0 files, 0 digests recorded again\n")
        )
    );

    let lone = write(
        &scratch.path,
        "annotated.vcf",
        &text("21\t100\t{first}\tA\tG\n"),
    );
    assert_eq!(run(&["identifiers", &lone]).code, 0);
    assert_eq!(read(&lone), text("21\t100\t{FIRST}\tA\tG\n"));

    let blocked = under(Path::new(&corpus), "all_vars.gz");
    fs::write(&blocked, bgzf_block(text("21\t{first}\t.\n").as_bytes())).unwrap();
    let listing = text("21\t100\trs1\tA\tG\t.\tPASS\tvep=G|{second}\n");
    write(Path::new(&corpus), "variants.vcf", &listing);
    assert_eq!(
        run(&["identifiers", &corpus]),
        failed(&format!(
            "{blocked} holds 1 identifier under restricted redistribution terms and is \
             block-compressed: decompress it, pass the text through this subcommand, then \
             compress and index it again; nothing was written"
        ))
    );
    assert_eq!(
        read(&variants),
        listing,
        "the refusal comes before any file is written"
    );
    let named_first_text = text("21\t100\t{first}\tA\tG\n");
    let named_first = write(&scratch.path, "first.vcf", &named_first_text);
    let absent = under(&scratch.path, "absent");
    for refused in [&corpus, &absent] {
        assert_eq!(run(&["identifiers", &named_first, refused]).code, 1);
        assert_eq!(
            read(&named_first),
            named_first_text,
            "a path named before the refused one is not written either"
        );
    }
    assert_eq!(
        run(&["identifiers", &absent]),
        failed(&format!("{absent} is neither a file nor a directory"))
    );
}

/// Made-up identifiers of the second shape under every other letter pair and digit count it
/// takes, their digits a counting run: each its letters and its digits.
const SECOND_SHAPE_FORMS: [(&str, &str); 25] = [
    ("CS", "123456"),
    ("CR", "123456"),
    ("CX", "123456"),
    ("CG", "123456"),
    ("CP", "123456"),
    ("CP", "12345"),
    ("CM", "1234567"),
    ("CD", "1234567"),
    ("CI", "1234567"),
    ("CS", "1234567"),
    ("CR", "1234567"),
    ("CX", "1234567"),
    ("CG", "1234567"),
    ("BM", "1234567"),
    ("HM", "1234"),
    ("HD", "1234"),
    ("HI", "1234"),
    ("HR", "1234"),
    ("HS", "1234"),
    ("HM", "123456"),
    ("HD", "123456"),
    ("HI", "123456"),
    ("HR", "123456"),
    ("HS", "123456"),
    ("HX", "123456"),
];

/// The header line of an output table in the fixed columns.
const FIXED_COLUMNS_HEADER: &str = concat!(
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\t",
    "cDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\t",
    "Existing_variation\tExtra\n",
);

/// One row of an output table in the fixed columns, with `cell` as its `Existing_variation`
/// cell.
fn fixed_columns_row(cell: &str) -> String {
    format!(
        "rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t5\t5\t2\tA/G\tgCc/gGc\t\
         {cell}\tIMPACT=MODERATE\n"
    )
}

/// The fabricated form of a made-up identifier: its letters, then a zero a digit.
fn zeroed(letters: &str, digits: &str) -> String {
    made_up(letters, &"0".repeat(digits.len()))
}

#[test]
fn fabricated_identifiers_read_every_letter_pair_and_length_of_the_second_shape() {
    // The second shape is read under each letter pair and digit count it takes (`C` and a
    // type letter with six or seven digits, `CP` with five or six, `BM` with seven, `H` and a
    // type letter with four or six, `HX` with six), as a list member and as a JSON `id`,
    // within the bounds of the first pairs: not before a version, after a colon, alone in a
    // cell, in a longer word or in prose.
    for (letters, digits) in SECOND_SHAPE_FORMS {
        let identifier = made_up(letters, digits);
        let zero = zeroed(letters, digits);
        let replaced = [
            "vep=G|rs1&{}&rs2;AF=0.1",
            "|{}|",
            "{},rs1\trs2;{}\n",
            r#"{"id":"{}"} {"id": "{}"}"#,
        ];
        for shape in replaced {
            let (given, expected) = (shape.replace("{}", &identifier), shape.replace("{}", &zero));
            assert_eq!(fabricated(&given), expected);
            assert_eq!(
                fabricated(&expected),
                expected,
                "the fabricated form is left as it is"
            );
            assert!(count(&given) > 0 && count(&expected) == 0, "{expected}");
        }
        assert!(is_fabricated(&zero) && !is_fabricated(&identifier));
        let kept = [
            "1\t{}.1\tchr1",
            "CLNDISDB=MedGen:{},OMIM:1",
            "see {} here",
            "{}\t1\n1\t{}\n",
            r#""synonyms":["{}","chr1"]"#,
            "|A{}|",
            "|{}A|",
        ];
        for shape in kept {
            let given = shape.replace("{}", &identifier);
            assert_eq!(fabricated(&given), given, "{shape} under {zero}");
            assert_eq!(count(&given), 0, "{shape} under {zero}");
        }
    }
    // A pair or a digit count outside those, a sample name among them, keeps its text.
    let outside = [
        ("CN", "123456"),
        ("CE", "123456"),
        ("CM", "12345"),
        ("CM", "12345678"),
        ("CP", "1234"),
        ("CP", "1234567"),
        ("BM", "123456"),
        ("BM", "12345678"),
        ("BD", "1234567"),
        ("HM", "12345"),
        ("HM", "1234567"),
        ("HX", "1234"),
        ("HG", "1234"),
        ("HG", "00096"),
        ("NA", "12878"),
        ("AM", "123456"),
        ("DM", "123456"),
    ];
    for (letters, digits) in outside {
        let token = made_up(letters, digits);
        for shape in ["rs1&{}", "|{}|", "{},{}", r#"{"id":"{}"}"#] {
            let given = shape.replace("{}", &token);
            assert_eq!(fabricated(&given), given);
            assert_eq!(count(&given), 0, "{given}");
        }
    }
}

#[test]
fn a_second_shape_identifier_alone_in_the_existing_variation_cell_is_read() {
    // The `Existing_variation` column of an output table names known variants and never a
    // sequence, so a second-shape identifier that is the whole cell of that column is read,
    // under the fixed columns and under a tab output's own, a line ending in a carriage
    // return too, and from there wherever it stands as a word in the text.
    let (second, zero) = (made_up("CD", "123456"), zeroed("CD", "123456"));
    let (header, row) = (FIXED_COLUMNS_HEADER, fixed_columns_row);
    let described = "## Existing_variation : Identifier(s) of co-located known variants\n";
    let table = format!(
        "{described}{header}{}{}{}",
        row(&second),
        row("-"),
        row("rs1")
    );
    assert_eq!(
        fabricated(&table),
        format!(
            "{described}{header}{}{}{}",
            row(&zero),
            row("-"),
            row("rs1")
        )
    );
    let fabricated_table = format!("{described}{header}{}", row(&zero));
    assert_eq!(fabricated(&fabricated_table), fabricated_table);
    assert_eq!(count(&table), 1);
    assert_eq!(
        restricted_identifiers(&table),
        BTreeSet::from([second.clone()])
    );
    assert_eq!(restricted_identifiers(&fabricated_table), none());

    let own = |cell: &str| {
        format!(
            "#Uploaded_variation\tLocation\tExisting_variation\tSYMBOL\n\
             rs1\t21:100\t{cell}\tGENE1\n{cell}\t21:200\t-\tGENE1\n"
        )
    };
    // Read in the cell, then replaced where it names the second row.
    assert_eq!(fabricated(&own(&second)), own(&zero));
    assert_eq!(count(&own(&second)), 2);
    let last = |cell: &str| format!("#Uploaded_variation\tExisting_variation\r\nrs1\t{cell}\r\n");
    assert_eq!(fabricated(&last(&second)), last(&zero));
    for (letters, digits) in SECOND_SHAPE_FORMS {
        assert_eq!(
            fabricated(&format!("{header}{}", row(&made_up(letters, digits)))),
            format!("{header}{}", row(&zeroed(letters, digits)))
        );
    }

    // The same token keeps its text alone in that column of a table with no header line, in
    // a table whose header does not name the column, in another column of a table that names
    // it, below a `##` line that names it and in the ID column of a variant file.
    let kept = [
        row(&second),
        format!("{described}{}", row(&second)),
        format!("#Uploaded_variation\tLocation\tAllele\n{second}\t21:100\tG\n"),
        format!("{header}{}", row("-").replace("ENSG1", &second)),
        format!(
            "{header}{}#Uploaded_variation\tLocation\n{second}\t21:100\n",
            row("-")
        ),
        format!(
            "{header}{}{}{}",
            row("-"),
            header.replace("Existing_variation", "Known"),
            row(&second)
        ),
        format!("##Uploaded_variation\tExisting_variation\nrs1\t{second}\n"),
        format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             21\t100\t{second}\tA\tG\t.\tPASS\t.\n"
        ),
        format!(
            "{header}{}{}{}",
            row(&format!("{second}.1")),
            row(&made_up("CN", "123456")),
            row(&made_up("HG", "00096"))
        ),
    ];
    for (place, given) in kept.iter().enumerate() {
        assert_eq!(&fabricated(given), given, "kept case {place}");
        assert_eq!(count(given), 0, "kept case {place}");
        assert_eq!(restricted_identifiers(given), none(), "kept case {place}");
    }
}

#[test]
fn identifiers_reads_a_corpus_that_holds_an_identifier_alone_in_a_cell_of_every_file() {
    // A corpus whose stored outputs are tables alone carries a co-located identifier of the
    // second shape alone in the `Existing_variation` cell of every file: `identifiers
    // --check` and the refusal of `check` name each file with its count, and `identifiers`
    // writes the fabricated form into each cell. An identifier of another letter pair or
    // digit count that a stored output lists is refused and replaced the same way, and a
    // file that holds one such identifier and nothing else of either shape is read
    // whichever letters open it.
    let (second, zero) = (made_up("CD", "123456"), zeroed("CD", "123456"));
    let scratch = scratch("identifiers-cells");
    let corpus = under(&scratch.path, "corpus");
    let expected = Path::new(&corpus).join("expected");
    fs::create_dir_all(&expected).unwrap();
    let default = write(
        &expected,
        "default.txt",
        &format!("{FIXED_COLUMNS_HEADER}{}", fixed_columns_row(&second)),
    );
    let tab = under(&expected, "tab.txt.gz");
    let tab_text = |cell: &str| {
        format!(
            "#Uploaded_variation\tLocation\tAllele\tExisting_variation\nrs1\t21:100\tG\t{cell}\n"
        )
    };
    fs::write(&tab, gzip(tab_text(&second).as_bytes(), None, 1)).unwrap();

    assert_eq!(
        run(&["identifiers", "--check", &corpus]),
        Run {
            code: 1,
            stdout: format!(
                "{corpus}: 2 files read, 2 hold an identifier in another form than the fabricated one\n"
            ),
            stderr: format!(
                "ERROR: [build_golden_corpus] {default} holds 1 identifier {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] {tab} holds 1 identifier {RESTRICTED_REFUSAL}\n"
            ),
        }
    );
    assert_eq!(
        restricted_identifier_errors(&corpus_files(&corpus).unwrap(), Some(&corpus)).unwrap(),
        [
            format!("expected/default.txt holds 1 identifier {RESTRICTED_REFUSAL}"),
            format!("expected/tab.txt.gz holds 1 identifier {RESTRICTED_REFUSAL}"),
        ]
    );
    assert_eq!(
        run(&["identifiers", &corpus]),
        clean(
            0,
            &format!("{corpus}: 2 identifiers replaced in 2 files, 0 digests recorded again\n")
        )
    );
    assert_eq!(
        read(&default),
        format!("{FIXED_COLUMNS_HEADER}{}", fixed_columns_row(&zero))
    );
    assert_eq!(gunzip(&fs::read(&tab).unwrap()), tab_text(&zero).as_bytes());
    assert_eq!(run(&["identifiers", "--check", &corpus]).code, 0);

    let listed_text = "21\t100\trs1\tA\tG\t.\tPASS\tCSQ=G|rs1&[CS:123456]|x,G|[CM:1234567]|y\n";
    let listed = write(&expected, "vcf.vcf", &text(listed_text));
    let colocated = concat!(
        r#"{"input":"21\t100\trs1\tA\tG","colocated_variants":"#,
        r#"[{"id":"[HM:1234]"},{"id":"[BM:1234567]"},{"id":"[CP:12345]"}]}"#,
        "\n",
    );
    let stored = write(&expected, "json.jsonl", &text(colocated));
    assert_eq!(
        run(&["identifiers", "--check", &corpus]),
        Run {
            code: 1,
            stdout: format!(
                "{corpus}: 4 files read, 2 hold an identifier in another form than the fabricated one\n"
            ),
            stderr: format!(
                "ERROR: [build_golden_corpus] {stored} holds 3 identifiers {RESTRICTED_REFUSAL}\n\
                 ERROR: [build_golden_corpus] {listed} holds 2 identifiers {RESTRICTED_REFUSAL}\n"
            ),
        }
    );
    assert_eq!(
        run(&["identifiers", &corpus]),
        clean(
            0,
            &format!("{corpus}: 5 identifiers replaced in 2 files, 0 digests recorded again\n")
        )
    );
    let listed_zero = listed_text
        .replace("[CS:123456]", "[CS:000000]")
        .replace("[CM:1234567]", "[CM:0000000]");
    assert_eq!(read(&listed), text(&listed_zero));
    let colocated_zero = colocated
        .replace("[HM:1234]", "[HM:0000]")
        .replace("[BM:1234567]", "[BM:0000000]")
        .replace("[CP:12345]", "[CP:00000]");
    assert_eq!(read(&stored), text(&colocated_zero));
    assert_eq!(run(&["identifiers", "--check", &corpus]).code, 0);

    let lone = under(&scratch.path, "lone.txt");
    for (letters, digits) in SECOND_SHAPE_FORMS {
        fs::write(&lone, format!("rs1,{}\n", made_up(letters, digits))).unwrap();
        assert_eq!(
            restricted_identifier_errors(std::slice::from_ref(&lone), None).unwrap(),
            [format!("{lone} holds 1 identifier {RESTRICTED_REFUSAL}")],
            "{}",
            zeroed(letters, digits)
        );
    }
}

#[test]
fn no_committed_corpus_holds_a_restricted_identifier() {
    // Every text file under `tests/golden`, inflated, carries identifiers of the two
    // restricted shapes in the fabricated form alone, each corpus read as `check` reads it;
    // the corpus that carries that form is the proof the walk reads the files.
    let root = common::repository().join("tests/golden");
    let root = root.to_str().unwrap();
    for manifest in committed_manifests() {
        let corpus = manifest.parent().unwrap().to_str().unwrap();
        let files = corpus_files(corpus).unwrap();
        let refused = restricted_identifier_errors(&files, Some(root)).unwrap();
        assert_eq!(refused, Vec::<String>::new(), "{corpus}");
    }
    let carrying = corpus_files(&format!("{root}/116/GRCh37-hgvs")).unwrap();
    let fabricated: usize = carrying
        .iter()
        .filter_map(|file| read_corpus_text(file).unwrap())
        .map(|read| RESTRICTED_IDENTIFIER.find_iter(&read.text).len())
        .sum();
    assert!(fabricated > 0);
}

#[test]
fn the_files_of_a_directory_are_counted_alike_on_any_count_of_threads() {
    // Forty files, every third one holding an identifier: the lines come in the order of the
    // walk whether one thread reads the files or several do.
    let scratch = scratch("identifiers-threads");
    let corpus = under(&scratch.path, "corpus");
    fs::create_dir(&corpus).unwrap();
    let mut expected = String::new();
    for place in 0..40 {
        let name = format!("f{place:02}.txt");
        let held = if place % 3 == 0 {
            first()
        } else {
            zeroed_first()
        };
        let path = write(Path::new(&corpus), &name, &format!("x|{held}|y\n"));
        if place % 3 == 0 {
            expected += &format!(
                "ERROR: [build_golden_corpus] {path} holds 1 identifier {RESTRICTED_REFUSAL}\n"
            );
        }
    }
    for bound in ["1", "2", "7", "64"] {
        let output = std::process::Command::new(common::BIN)
            .args(["identifiers", "--check", &corpus])
            .env("VEP_GOLDEN_CORPUS_THREADS", bound)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{bound}");
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            expected,
            "{bound}"
        );
        assert_eq!(
            String::from_utf8(output.stdout).unwrap(),
            format!(
                "{corpus}: 40 files read, 14 hold an identifier in another form than the fabricated one\n"
            ),
            "{bound}"
        );
    }
}
