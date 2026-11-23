// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp stop-codon predicates held to stored answers:
//! `tests/data/hgvsp_stop/cases.jsonl.gz` is a draw of calls, one JSON object a line with
//! the function, its arguments and the answer as `json.dumps` writes it (a tuple a list;
//! `{"raise": class, "text": text}` for an exception that ends the run).
//! `VEP_FORMATS_COMPARE_STOP_CASES` names another draw to read in its place.

use std::io::Read;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::hgvsp_stop as stop;

/// The cases of `path`, a line each, gunzipped where the file is gzip.
fn cases(path: &Path) -> Vec<Value> {
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

fn text<'a>(args: &'a Value, name: &str) -> &'a str {
    args.get(name)
        .and_then(Value::as_str)
        .unwrap_or_else(|| panic!("no text argument {name}"))
}

fn texts(values: impl IntoIterator<Item = String>) -> Value {
    Value::List(values.into_iter().map(Value::Str).collect())
}

/// The tool's answer to one case, as the value `json.dumps` is given.
fn answer(function: &str, args: &Value) -> Value {
    let row = || {
        (
            text(args, "consequence"),
            text(args, "flags"),
            text(args, "ev"),
            text(args, "av"),
        )
    };
    let pair = || (text(args, "reference"), text(args, "candidate"));
    let terms = || text(args, "consequence");
    match function {
        "stop_codon_protein_comparison" => {
            let (consequence, flags, ev, av) = row();
            Value::Bool(
                stop::stop_codon_protein_comparison(consequence, flags, ev, av)
                    .expect("the class reads the pair"),
            )
        }
        "stop_named_on_incomplete_coding_sequence" => {
            let (consequence, flags, ev, av) = row();
            Value::Bool(stop::stop_named_on_incomplete_coding_sequence(
                consequence,
                flags,
                ev,
                av,
            ))
        }
        "residue_names" => texts(stop::residue_names(text(args, "text"))),
        "stop_extension" => match stop::stop_extension(text(args, "description")) {
            Some((number, first, count)) => texts([number, first, count]),
            None => Value::Null,
        },
        "stop_replaced_by_a_stop" => Value::Bool(stop::stop_replaced_by_a_stop(pair().0, pair().1)),
        "stop_replaced_by_first_residues" => {
            Value::Bool(stop::stop_replaced_by_first_residues(pair().0, pair().1))
        }
        "stop_replaced_through_an_inserted_stop" => Value::Bool(
            stop::stop_replaced_through_an_inserted_stop(pair().0, pair().1),
        ),
        "residues_inserted_before_the_stop" => {
            Value::Bool(stop::residues_inserted_before_the_stop(pair().0, pair().1))
        }
        "lost_stop_written_as_unchanged" => Value::Bool(stop::lost_stop_written_as_unchanged(
            terms(),
            pair().0,
            pair().1,
        )),
        "range_to_the_stop_replaced" => Value::Bool(stop::range_to_the_stop_replaced(
            terms(),
            pair().0,
            pair().1,
        )),
        "unknown_count_where_a_stop_is_met" => Value::Bool(
            stop::unknown_count_where_a_stop_is_met(terms(), pair().0, pair().1),
        ),
        "stop_deleted_with_an_extension_count" => Value::Bool(
            stop::stop_deleted_with_an_extension_count(terms(), pair().0, pair().1),
        ),
        other => panic!("no function {other}"),
    }
}

/// Every case of the draw has its stored answer; a case stored with an exception ends the
/// tool's call with the exception's text.
#[test]
fn every_made_call_has_its_stored_answer() {
    let path = std::env::var_os("VEP_FORMATS_COMPARE_STOP_CASES").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hgvsp_stop/cases.jsonl.gz"),
        PathBuf::from,
    );
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        if info.payload().downcast_ref::<Uncaught>().is_none() {
            default_hook(info);
        }
    }));
    let cases = cases(&path);
    let mut differing = Vec::new();
    let (mut accepted, mut ended) = (0, 0);
    for case in &cases {
        let function = text(case, "f");
        let args = case.get("args").expect("a case has arguments");
        let want = case.get("want").expect("a case has an answer");
        let got = match catch_unwind(AssertUnwindSafe(|| answer(function, args))) {
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
        accepted += usize::from(got == Value::Bool(true));
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
        cases.len(),
        differing[..differing.len().min(20)].join("\n")
    );
    // The draw holds calls of every answer.
    assert!(cases.len() >= 500 && accepted >= 40 && ended >= 3);
}
