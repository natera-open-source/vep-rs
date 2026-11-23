// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The shift walk, the reference sequence and the sequence-edge class held to stored
//! answers: `tests/data/sequence_edges/cases.jsonl.gz` is a draw of calls, one JSON object
//! a line with the function, its arguments and the answer as `json.dumps` writes it (a
//! tuple a list; `{"raise": class, "text": text}` for an exception that ends the run). Its
//! first line is the fixture the calls read: a
//! reference FASTA, its index and the shards of a candidate cache.
//! `VEP_FORMATS_COMPARE_SEQUENCE_EDGE_CASES` names another draw to read in its place.
//!
//! `tests/data/sequence_edges/construct.jsonl` holds made references and caches, one a
//! line: the files as the hex of their bytes, the directories, the regions asked for, and
//! the stored answer for them, the directory they are planted in written `<dir>`.

use std::io::{Read, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::sequence_edges::{
    edge_transcript, indel_span, perform_shift_walk, turned_allele, EdgeTranscript, EdgeWalks,
    SequenceEdges,
};

/// The lines of `path`, a JSON object each, gunzipped where the file is gzip.
fn lines(path: &Path) -> Vec<Value> {
    let bytes = std::fs::read(path).expect("the cases are read");
    let mut text = String::new();
    if bytes.starts_with(&[0x1f, 0x8b]) {
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .expect("the cases inflate to UTF-8");
    } else {
        text = String::from_utf8(bytes).expect("the cases are UTF-8");
    }
    text.split('\n')
        .filter(|line| !line.is_empty())
        .map(|line| pyjson::loads(line).expect("a case is JSON"))
        .collect()
}

fn arg<'a>(args: &'a Value, name: &str) -> &'a Value {
    args.get(name)
        .unwrap_or_else(|| panic!("no argument {name}"))
}

fn text<'a>(args: &'a Value, name: &str) -> &'a str {
    arg(args, name)
        .as_str()
        .unwrap_or_else(|| panic!("the argument {name} is no text"))
}

fn whole(args: &Value, name: &str) -> i64 {
    match arg(args, name) {
        Value::Int(number) => *number,
        other => panic!("the argument {name} is no whole number: {other:?}"),
    }
}

fn flag(args: &Value, name: &str) -> bool {
    match arg(args, name) {
        Value::Bool(flag) => *flag,
        other => panic!("the argument {name} is no flag: {other:?}"),
    }
}

fn optional_text(value: &Value) -> Option<String> {
    value.as_str().map(str::to_string)
}

fn optional_whole(value: &Value) -> Option<i64> {
    match value {
        Value::Int(number) => Some(*number),
        _ => None,
    }
}

fn optional_flag(value: &Value) -> Option<bool> {
    match value {
        Value::Bool(flag) => Some(*flag),
        _ => None,
    }
}

/// An entry as a case writes it: its members in the order `Entry` declares them.
fn entry(value: &Value) -> Option<Entry> {
    let parts = value.as_list()?;
    let part = |at: usize| parts[at].as_str().expect("a text").to_string();
    Some(Entry {
        key: part(0),
        record: optional_whole(&parts[1]),
        location: optional_text(&parts[2]),
        allele: part(3),
        feature: part(4),
        consequence: part(5),
        fields: parts[6]
            .as_list()
            .expect("the fields")
            .iter()
            .map(|pair| {
                let pair = pair.as_list().expect("a pair");
                (
                    pair[0].as_str().expect("a name").to_string(),
                    pair[1].as_str().expect("a value").to_string(),
                )
            })
            .collect(),
        indel: optional_flag(&parts[7]),
        insertion: optional_flag(&parts[8]),
        reference: optional_text(&parts[9]),
        line: optional_whole(&parts[10]),
        allele_location: optional_text(&parts[11]),
    })
}

fn wide(number: i128) -> Value {
    Value::Int(i64::try_from(number).expect("a case's number fits"))
}

fn optional<T>(value: Option<T>, made: impl FnOnce(T) -> Value) -> Value {
    value.map_or(Value::Null, made)
}

fn pair(first: i64, last: i64) -> Value {
    Value::List(vec![Value::Int(first), Value::Int(last)])
}

/// A transcript as a case writes one: its name, span, strand, exons and coding bounds.
fn transcript_parts(transcript: &EdgeTranscript) -> Vec<Value> {
    vec![
        Value::Str(transcript.name.clone()),
        Value::Int(transcript.start),
        Value::Int(transcript.end),
        Value::Int(transcript.strand),
        Value::List(
            transcript
                .exons
                .iter()
                .map(|(start, end)| pair(*start, *end))
                .collect(),
        ),
        optional(transcript.coding, |(first, last)| pair(first, last)),
    ]
}

fn bodies(bodies: (Option<String>, Option<String>)) -> Value {
    Value::List(vec![
        optional(bodies.0, Value::Str),
        optional(bodies.1, Value::Str),
    ])
}

/// The reference, its index and the cache of the draw's first line, planted under
/// `directory`.
fn planted(fixture: &Value, directory: &Path) -> SequenceEdges {
    let fasta = directory.join("reference.fa");
    std::fs::write(&fasta, text(fixture, "fasta")).expect("the reference is written");
    std::fs::write(directory.join("reference.fa.fai"), text(fixture, "fai"))
        .expect("the index is written");
    let cache = directory.join("cache");
    for shard in arg(fixture, "shards").as_list().expect("the shards") {
        let shard = shard.as_list().expect("a shard");
        let part = |at: usize| shard[at].as_str().expect("a text");
        let region = cache.join("transcripts").join(part(0));
        std::fs::create_dir_all(&region).expect("the region's directory is made");
        let path = region.join(part(1));
        if part(1).ends_with(".gz") {
            let file = std::fs::File::create(&path).expect("the shard is made");
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            encoder
                .write_all(part(2).as_bytes())
                .expect("the shard is written");
            encoder.finish().expect("the shard is closed");
        } else {
            std::fs::write(&path, part(2)).expect("the shard is written");
        }
    }
    SequenceEdges::new(&fasta, &cache).expect("the fixture reads")
}

/// The tool's answer to one case, as the value `json.dumps` is given.
fn answer(function: &str, args: &Value, edges: &SequenceEdges) -> Value {
    let filed = || {
        let key = (
            text(args, "region").to_string(),
            text(args, "stable_id").to_string(),
        );
        edges
            .transcripts
            .get(&key)
            .expect("the cache holds the transcript")
    };
    let row = || {
        (
            entry(arg(args, "e")).expect("an entry"),
            entry(arg(args, "record")),
        )
    };
    let walks = || EdgeWalks {
        own: whole(args, "own"),
        cut: whole(args, "cut"),
    };
    match function {
        "transcripts" => {
            let mut held: Vec<(&(String, String), &EdgeTranscript)> =
                edges.transcripts.iter().collect();
            held.sort_by(|a, b| a.0.cmp(b.0));
            Value::List(
                held.into_iter()
                    .map(|((region, stable_id), transcript)| {
                        let mut parts =
                            vec![Value::Str(region.clone()), Value::Str(stable_id.clone())];
                        parts.extend(transcript_parts(transcript));
                        Value::List(parts)
                    })
                    .collect(),
            )
        }
        "perform_shift_walk" => Value::Int(perform_shift_walk(
            text(args, "head"),
            text(args, "tail"),
            whole(args, "size"),
            text(args, "pre"),
            text(args, "post"),
            flag(args, "forward"),
            flag(args, "whole_flanks"),
        )),
        "turned_allele" => Value::Str(turned_allele(
            text(args, "allele"),
            whole(args, "steps"),
            flag(args, "forward"),
        )),
        "edge_transcript" => optional(edge_transcript(arg(args, "record")), |transcript| {
            Value::List(transcript_parts(&transcript))
        }),
        "notation" => optional(
            filed().notation(
                i128::from(whole(args, "first")),
                i128::from(whole(args, "last")),
                text(args, "change"),
            ),
            Value::Str,
        ),
        "coordinate" => Value::List(vec![
            wide(filed().coordinate(whole(args, "position"))),
            wide(filed().length()),
        ]),
        "bases" => {
            let region = text(args, "region");
            Value::List(vec![
                optional(edges.reference.length(region), Value::Int),
                Value::Str(
                    edges
                        .reference
                        .bases(region, whole(args, "start"), whole(args, "end")),
                ),
            ])
        }
        "walks" => {
            let (e, record) = row();
            optional(edges.walks(&e, record.as_ref()), |walks| {
                pair(walks.own, walks.cut)
            })
        }
        "indel_span" => {
            let (e, record) = row();
            optional(indel_span(&e, record.as_ref()), |(region, start, end)| {
                Value::List(vec![Value::Str(region), Value::Int(start), Value::Int(end)])
            })
        }
        "accepts" => {
            let (e, record) = row();
            let a = entry(arg(args, "a"));
            Value::Bool(edges.accepts(&e, a.as_ref(), record.as_ref(), text(args, "name")))
        }
        "_deletion" => {
            bodies(edges._deletion(filed(), &walks(), whole(args, "start"), whole(args, "end")))
        }
        "_insertion" => bodies(edges._insertion(
            filed(),
            &walks(),
            text(args, "region"),
            whole(args, "start"),
            whole(args, "end"),
            text(args, "allele"),
        )),
        other => panic!("no function {other}"),
    }
}

/// Every case of the draw has its stored answer; a case stored with an exception ends the
/// tool's call with the exception's text.
#[test]
fn every_made_call_has_its_stored_answer() {
    let path = std::env::var_os("VEP_FORMATS_COMPARE_SEQUENCE_EDGE_CASES").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sequence_edges/cases.jsonl.gz"),
        PathBuf::from,
    );
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        if info.payload().downcast_ref::<Uncaught>().is_none() {
            default_hook(info);
        }
    }));
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join("sequence-edge-cases");
    let _ = std::fs::remove_dir_all(&directory);
    std::fs::create_dir_all(&directory).expect("a directory for the fixture");
    let cases = lines(&path);
    let edges = planted(arg(&cases[0], "fixture"), &directory);
    let mut differing = Vec::new();
    let (mut accepted, mut ended) = (0, 0);
    for case in &cases[1..] {
        let function = text(case, "f");
        let args = arg(case, "args");
        let want = arg(case, "want");
        let got = match catch_unwind(AssertUnwindSafe(|| answer(function, args, &edges))) {
            Ok(got) => got,
            Err(payload) => match payload.downcast::<Uncaught>() {
                Ok(uncaught) => Value::Str(uncaught.0),
                Err(payload) => std::panic::resume_unwind(payload),
            },
        };
        let same = match want.get("raise") {
            Some(_) => {
                ended += 1;
                want.get("text") == Some(&got)
            }
            None => *want == got,
        };
        accepted += usize::from(function == "accepts" && got == Value::Bool(true));
        if !same {
            differing.push(format!(
                "{function} {}: stored {}, the tool {}",
                pyjson::dumps_compact(args),
                pyjson::dumps_compact(want),
                pyjson::dumps_compact(&got)
            ));
        }
    }
    assert!(
        differing.is_empty(),
        "{} of {} calls differ:\n{}",
        differing.len(),
        cases.len() - 1,
        differing[..differing.len().min(12)].join("\n")
    );
    // The draw holds calls of every answer.
    assert!(cases.len() >= 500 && accepted >= 40 && ended >= 10);
    std::fs::remove_dir_all(&directory).expect("the fixture is removed");
}

fn bytes_of(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|at| u8::from_str_radix(&hex[2 * at..2 * at + 2], 16).expect("two hex digits"))
        .collect()
}

/// What the tool makes of one made reference and cache, as the value a case holds.
fn constructed(case: &Value, directory: &Path) -> Value {
    for made in arg(case, "dirs").as_list().expect("the directories") {
        std::fs::create_dir_all(directory.join(made.as_str().expect("a path")))
            .expect("a directory is made");
    }
    for file in arg(case, "files").as_list().expect("the files") {
        let file = file.as_list().expect("a file");
        let path = directory.join(file[0].as_str().expect("a path"));
        std::fs::create_dir_all(path.parent().expect("a file has a directory"))
            .expect("the file's directory is made");
        std::fs::write(path, bytes_of(file[1].as_str().expect("the bytes"))).expect("a file");
    }
    let edges = match SequenceEdges::new(
        &directory.join(text(case, "fasta")),
        &directory.join(text(case, "cache")),
    ) {
        Ok(edges) => edges,
        Err(error) => {
            let text = error
                .to_string()
                .replace(&directory.display().to_string(), "<dir>");
            return Value::Dict(
                [("error".to_string(), Value::Str(text))]
                    .into_iter()
                    .collect(),
            );
        }
    };
    let lengths = edges
        .reference
        .regions
        .iter()
        .map(|(name, entry)| Value::List(vec![Value::Str(name.clone()), Value::Int(entry.0)]))
        .collect();
    let bases = arg(case, "bases")
        .as_list()
        .expect("the regions asked for")
        .iter()
        .map(|probe| {
            let probe = probe.as_list().expect("a region and its span");
            let span = |at: usize| optional_whole(&probe[at]).expect("a position");
            let region = probe[0].as_str().expect("a region");
            match catch_unwind(AssertUnwindSafe(|| {
                edges.reference.bases(region, span(1), span(2))
            })) {
                Ok(bases) => Value::Str(bases),
                Err(payload) => match payload.downcast::<Uncaught>() {
                    Ok(uncaught) => Value::Str(format!("raise: {}", uncaught.0)),
                    Err(payload) => std::panic::resume_unwind(payload),
                },
            }
        })
        .collect();
    Value::Dict(
        [
            ("lengths".to_string(), Value::List(lengths)),
            (
                "transcripts".to_string(),
                answer("transcripts", &Value::Null, &edges),
            ),
            ("bases".to_string(), Value::List(bases)),
        ]
        .into_iter()
        .collect(),
    )
}

/// A stored answer with its exception written as the tool's call ends: its text after
/// `raise: `.
fn raised_as_text(want: &Value) -> Value {
    match want {
        Value::Dict(members) => match (members.get("raise"), members.get("text")) {
            (Some(_), Some(Value::Str(text))) => Value::Str(format!("raise: {text}")),
            _ => Value::Dict(
                members
                    .iter()
                    .map(|(name, member)| (name.clone(), raised_as_text(member)))
                    .collect(),
            ),
        },
        Value::List(items) => Value::List(items.iter().map(raised_as_text).collect()),
        other => other.clone(),
    }
}

/// A reference and a cache are read, refused or end the run as their stored answers say:
/// the index's lines and numbers, a reference that is compressed, absent or a directory,
/// the bases of a line the index misstates, the shards of a cache in the order of their
/// names, and the shards and records that are passed over.
#[test]
fn a_reference_and_a_cache_are_read_refused_or_end_the_run_as_stored() {
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        if info.payload().downcast_ref::<Uncaught>().is_none() {
            default_hook(info);
        }
    }));
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sequence_edges/construct.jsonl");
    let cases = lines(&path);
    assert!(cases.len() >= 60);
    let scratch = Path::new(env!("CARGO_TARGET_TMPDIR")).join("sequence-edge-constructions");
    let _ = std::fs::remove_dir_all(&scratch);
    for (at, case) in cases.iter().enumerate() {
        let directory = scratch.join(at.to_string());
        std::fs::create_dir_all(&directory).expect("a directory for the case");
        let got = match catch_unwind(AssertUnwindSafe(|| constructed(case, &directory))) {
            Ok(got) => got,
            Err(payload) => match payload.downcast::<Uncaught>() {
                Ok(uncaught) => Value::Str(format!("raise: {}", uncaught.0)),
                Err(payload) => std::panic::resume_unwind(payload),
            },
        };
        assert_eq!(
            pyjson::dumps_compact(&got),
            pyjson::dumps_compact(&raised_as_text(arg(case, "want"))),
            "{}",
            text(case, "name")
        );
    }
    std::fs::remove_dir_all(&scratch).expect("the cases' files are removed");
}

/// A whole number is held in 64 bits. A position beyond them is outside every region, so
/// its row has no span; a transcript record with such a coordinate makes no transcript; and
/// an index line with such a number is no index line.
#[test]
fn a_number_beyond_64_bits_makes_no_span_no_transcript_and_no_index_line() {
    let beyond = "9223372036854775808";
    let row = Entry {
        key: format!("1:{beyond}|-|ENST00000000001|Transcript"),
        location: Some(format!("1:{beyond}")),
        allele: "-".to_string(),
        ..Entry::default()
    };
    assert_eq!(indel_span(&row, None), None);
    let within = Entry {
        location: Some("1:9223372036854775807".to_string()),
        ..row
    };
    assert_eq!(
        indel_span(&within, None),
        Some(("1".to_string(), i64::MAX, i64::MAX))
    );
    let record = |start: &str| {
        pyjson::loads(&format!(
            r#"{{"stable_id": "ENST00000000001", "start": {start}, "end": 9, "strand": 1}}"#
        ))
        .expect("the record is JSON")
    };
    assert!(edge_transcript(&record("1")).is_some());
    assert!(edge_transcript(&record(beyond)).is_none());
    assert!(edge_transcript(&record(&format!("\"{beyond}\""))).is_none());
    assert!(edge_transcript(&record("1e19")).is_none());

    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join("sequence-edge-wide-numbers");
    let _ = std::fs::remove_dir_all(&directory);
    std::fs::create_dir_all(directory.join("cache/transcripts")).expect("a cache is made");
    let fasta = directory.join("reference.fa");
    std::fs::write(&fasta, ">one\nACGT\n").expect("the reference is written");
    let index = directory.join("reference.fa.fai");
    std::fs::write(&index, format!("one\t{beyond}\t5\t4\t5\n")).expect("the index is written");
    let refused = SequenceEdges::new(&fasta, &directory.join("cache"))
        .expect_err("the index is refused")
        .to_string();
    assert_eq!(
        refused,
        format!(
            "{}: not a FASTA index line: one\t{beyond}\t5\t4\t5",
            index.display()
        )
    );
    std::fs::remove_dir_all(&directory).expect("the files are removed");
}
