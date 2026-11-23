// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The structural-variant functions against stored answers: the predicates of the
//! structural-variant classes, the mate-side derivation and the readers of a Location, an
//! allele and an INFO column.
//!
//! `tests/data/structural/cases.jsonl` holds one case a line, a JSON object: `f` names
//! what is asked, `a` holds the arguments and `x` the stored answer. A whole number of
//! an answer is its decimal text, a float `["float", repr]`, a set a sorted list, a record
//! its CHROM, POS, ID, REF, ALT and INFO columns, a tuple its five members; a run that ends
//! on an input no reading expects is `{"raises": text}` (`*` where the tool words the ending in
//! its own way), the comparator's own error `{"error": text}` and a file that cannot be
//! opened `{"oserror": text}`. The answer of `variant` holds what the record functions
//! say of one record.
//!
//! The ignored test reads a file of the same form, with many more lines and with cases
//! that read planted caches and inputs, that `VEP_FORMATS_COMPARE_STRUCTURAL_CASES` names.

use indexmap::IndexMap;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::io::BufRead;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::pyfmt;
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::structural as s;
use vep_formats_compare::rules::structural::mate_context as m;
use vep_formats_compare::rules::structural::number::{Int, Number};
use vep_formats_compare::rules::structural::outputs as o;
use vep_formats_compare::shared;

fn text(value: &Value) -> &str {
    value.as_str().expect("a text")
}

fn optional(value: &Value) -> Option<&str> {
    value.as_str()
}

fn items(value: &Value) -> &[Value] {
    value.as_list().expect("a list")
}

fn whole(value: &Value) -> i64 {
    match value {
        Value::Int(number) => *number,
        other => panic!("not a whole number: {other:?}"),
    }
}

fn st(text: impl Into<String>) -> Value {
    Value::Str(text.into())
}

fn opt(text: Option<String>) -> Value {
    text.map_or(Value::Null, Value::Str)
}

fn list(values: impl IntoIterator<Item = Value>) -> Value {
    Value::List(values.into_iter().collect())
}

fn dict(members: Vec<(&str, Value)>) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect(),
    )
}

fn number(value: Number) -> Value {
    match value {
        Number::Int(own) => list([st("int"), st(own.to_string())]),
        Number::Float(own) => list([st("float"), st(pyfmt::repr_f64(own))]),
    }
}

fn tuple_of(value: &Value) -> s::AnnotTuple {
    let parts = items(value);
    (
        text(&parts[0]).to_string(),
        text(&parts[1]).to_string(),
        text(&parts[2]).to_string(),
        text(&parts[3]).to_string(),
        text(&parts[4]).to_string(),
    )
}

fn tuples_of(value: &Value) -> BTreeSet<s::AnnotTuple> {
    items(value).iter().map(tuple_of).collect()
}

fn listed(tuples: &BTreeSet<s::AnnotTuple>) -> Value {
    list(tuples.iter().map(|t| {
        list([
            st(t.0.as_str()),
            st(t.1.as_str()),
            st(t.2.as_str()),
            st(t.3.as_str()),
            st(t.4.as_str()),
        ])
    }))
}

fn variant(fields: &Value) -> s::VcfVariant {
    let f = items(fields);
    let (chrom, pos_str, vid, r#ref, alt, info) = (
        text(&f[0]),
        text(&f[1]),
        text(&f[2]),
        text(&f[3]),
        text(&f[4]),
        text(&f[5]),
    );
    let pos = s::perl_pos(pos_str);
    s::VcfVariant {
        chrom: chrom.to_string(),
        pos,
        vid: vid.to_string(),
        r#ref: r#ref.to_string(),
        alt: alt.to_string(),
        variant_type: s::classify_variant(r#ref, alt),
        location: s::normalize_location(&format!("{chrom}:{pos}")),
        info: info.to_string(),
        pos_text: Some(pos_str.to_string()),
    }
}

fn variants(value: &Value) -> Vec<s::VcfVariant> {
    items(value).iter().map(variant).collect()
}

fn span(found: Option<(String, Int, Int)>) -> Value {
    found.map_or(Value::Null, |(chrom, start, end)| {
        list([st(chrom), st(start.to_string()), st(end.to_string())])
    })
}

fn point(found: Option<(String, Int)>) -> Value {
    found.map_or(Value::Null, |(chrom, pos)| {
        list([st(chrom), st(pos.to_string())])
    })
}

fn pair(a: i128, b: i128) -> Value {
    list([st(a.to_string()), st(b.to_string())])
}

fn model_answer(model: Option<&s::TranscriptModel>) -> Value {
    let Some(model) = model else {
        return Value::Null;
    };
    dict(vec![
        ("chrom", st(model.chrom.as_str())),
        ("start", st(model.start.to_string())),
        ("end", st(model.end.to_string())),
        ("strand", st(model.strand.to_string())),
        ("biotype", st(model.biotype.as_str())),
        ("exons", list(model.exons.iter().map(|(a, b)| pair(*a, *b)))),
        ("cds", model.cds.map_or(Value::Null, |(a, b)| pair(a, b))),
        (
            "mirna",
            list(model.mirna_mature.iter().map(|(a, b)| pair(*a, *b))),
        ),
    ])
}

fn id_answer(id: &s::StableId) -> Vec<String> {
    match id {
        s::StableId::Text(own) => vec!["s".to_string(), own.clone()],
        s::StableId::Whole(own) => vec!["w".to_string(), own.clone()],
        s::StableId::Fraction(bits) => {
            vec!["f".to_string(), pyfmt::repr_f64(f64::from_bits(*bits))]
        }
        s::StableId::Nan => vec!["nan".to_string()],
    }
}

fn strings(values: &[String]) -> Value {
    list(values.iter().map(|value| st(value.as_str())))
}

fn authority_answer(cache: &Path) -> Value {
    let (by_chr, spans, models) = s::load_cache_authority(cache);
    let mut chromosomes: Vec<(String, Vec<Vec<String>>)> = by_chr
        .iter()
        .map(|(chrom, ids)| {
            let mut ids: Vec<Vec<String>> = ids.iter().map(id_answer).collect();
            ids.sort();
            (chrom.clone(), ids)
        })
        .collect();
    chromosomes.sort();
    let mut span_rows: Vec<(String, Vec<String>, Value)> = spans
        .iter()
        .map(|((chrom, id), span)| {
            (
                chrom.clone(),
                id_answer(id),
                list([
                    st(chrom.as_str()),
                    strings(&id_answer(id)),
                    st(span.0.as_str()),
                    st(span.1.to_string()),
                    st(span.2.to_string()),
                ]),
            )
        })
        .collect();
    span_rows.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    let mut model_rows: Vec<(String, Vec<String>, Value)> = models
        .iter()
        .map(|((chrom, id), model)| {
            (
                chrom.clone(),
                id_answer(id),
                list([
                    st(chrom.as_str()),
                    strings(&id_answer(id)),
                    model_answer(Some(model)),
                ]),
            )
        })
        .collect();
    model_rows.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    dict(vec![
        (
            "by_chr",
            list(
                chromosomes
                    .into_iter()
                    .map(|(chrom, ids)| list([st(chrom), list(ids.iter().map(|id| strings(id)))])),
            ),
        ),
        ("spans", list(span_rows.into_iter().map(|row| row.2))),
        ("models", list(model_rows.into_iter().map(|row| row.2))),
    ])
}

fn entry(value: &Value) -> Entry {
    let parts = items(value);
    Entry {
        key: text(&parts[0]).to_string(),
        record: match &parts[1] {
            Value::Null => None,
            other => Some(whole(other)),
        },
        location: optional(&parts[2]).map(str::to_string),
        allele: text(&parts[3]).to_string(),
        feature: text(&parts[4]).to_string(),
        consequence: text(&parts[5]).to_string(),
        ..Entry::default()
    }
}

fn side_table(a: &[Value]) -> Result<Value, Error> {
    let classes: Vec<&str> = items(&a[0]).iter().map(text).collect();
    let structural = s::StructuralVariants::new(
        &classes,
        optional(&a[1]).map(Path::new),
        optional(&a[2]).map(Path::new),
    )?;
    let entries: Vec<Entry> = items(&a[3]).iter().map(entry).collect();
    let tuples = structural.tuples(&entries);
    Ok(dict(vec![
        (
            "classes",
            list(structural.classes.iter().map(|class| st(class.as_str()))),
        ),
        (
            "records",
            Value::Int(structural.input_variants.len() as i64),
        ),
        ("spans", Value::Int(structural.spans.len() as i64)),
        (
            "tuples",
            list(tuples.iter().map(|(key, t)| {
                list([
                    st(key.as_str()),
                    list([
                        st(t.0.as_str()),
                        st(t.1.as_str()),
                        st(t.2.as_str()),
                        st(t.3.as_str()),
                        st(t.4.as_str()),
                    ]),
                ])
            })),
        ),
        (
            "locations",
            list(entries.iter().map(|e| opt(structural.location(e)))),
        ),
    ]))
}

fn variant_answer(fields: &Value, max_sv_size: i64) -> Value {
    let v = variant(fields);
    let (chrom, start, end) = s::record_location(&v);
    let (ref_chrom, ref_start) = s::reference_start(&v);
    dict(vec![
        ("pos", st(v.pos.to_string())),
        ("variant_type", st(v.variant_type.as_str())),
        ("location_text", st(v.location.as_str())),
        ("sequence_path", Value::Bool(s::perl_sequence_path(&v))),
        (
            "reference_start",
            list([st(ref_chrom), st(ref_start.to_string())]),
        ),
        ("skip_route", opt(s::perl_skip_route(&v, max_sv_size))),
        ("alleles", list(s::record_alleles(&v).into_iter().map(st))),
        (
            "location",
            list([st(chrom), st(start.to_string()), st(end.to_string())]),
        ),
        (
            "rows",
            list(
                s::vep_rs_record_rows(&v)
                    .into_iter()
                    .map(|(location, allele)| list([st(location), st(allele)])),
            ),
        ),
    ])
}

fn answer(f: &str, a: &[Value]) -> Result<Value, Error> {
    Ok(match f {
        "normalize_location" => st(s::normalize_location(text(&a[0]))),
        "location_span" => span(o::_location_span(text(&a[0]))),
        "non_integer_end" => Value::Bool(o::_non_integer_end(text(&a[0]))),
        "chromosome_of_location" => opt(s::_chromosome_of_location(text(&a[0]))),
        "point_location" => point(m::point_location(text(&a[0]))),
        "breakpoint_position" => {
            s::_breakpoint_position(text(&a[0])).map_or(Value::Null, |pos| st(pos.to_string()))
        }
        "strip_chr_in_bracket" => st(s::_strip_chr_in_bracket(text(&a[0]))),
        "normalize_allele" => st(s::normalize_allele(text(&a[0]), text(&a[1]))),
        "fastvep" => st(s::fastvep_symbolic_allele_to_so_class(text(&a[0]))),
        "own_breakend" => Value::Bool(o::is_own_breakend_allele(text(&a[0]))),
        "breakend" => Value::Bool(o::is_breakend_allele(text(&a[0]))),
        "alt_list" => Value::Bool(o::structural_alt_list(text(&a[0]))),
        "bracket_mate" => point(s::_bracket_mate(text(&a[0]))),
        "structural_allele" => Value::Bool(shared::is_structural_allele(text(&a[0]))),
        "perl_so_abbrev" => opt(s::perl_so_abbrev(text(&a[0]))),
        "perl_record_type" => st(s::perl_record_type(text(&a[0]), optional(&a[1]))),
        "unsupported_alt" => st(s::vep_rs_unsupported_alt_allele(
            text(&a[0]),
            optional(&a[1]),
        )),
        "classify_variant" => st(s::classify_variant(text(&a[0]), text(&a[1]))),
        "gvcf_site" => {
            let alts: Vec<String> = items(&a[0])
                .iter()
                .map(|alt| text(alt).to_string())
                .collect();
            Value::Bool(o::is_gvcf_variant_site(&alts, optional(&a[1])))
        }
        "perl_numeric_value" => number(s::perl_numeric_value(text(&a[0]))),
        "perl_pos" => st(s::perl_pos(text(&a[0])).to_string()),
        "not_an_integer" => Value::Bool(s::perl_number_is_not_an_integer(text(&a[0]))),
        "structural_end" => {
            let pos: Int = text(&a[0]).parse().expect("a position");
            let breakpoint = matches!(a[3], Value::Bool(true));
            s::perl_structural_end(pos, optional(&a[1]), optional(&a[2]), breakpoint)
                .map_or(Value::Null, number)
        }
        "info_fields" => Value::Dict(
            s::_info_fields(text(&a[0]))
                .into_iter()
                .map(|(key, value)| (key, st(value)))
                .collect::<IndexMap<String, Value>>(),
        ),
        "info_written" => opt(o::_info_written(text(&a[0]), text(&a[1]))),
        "variant" => variant_answer(&a[0], whole(&a[1])),
        "skipped" => {
            let (excluded, confirmed, problems, residue) = s::skipped_record_rows(
                &variants(&a[0]),
                &tuples_of(&a[1]),
                &tuples_of(&a[2]),
                whole(&a[3]),
                text(&a[4]),
            );
            list([
                listed(&excluded),
                Value::Int(confirmed),
                strings(&problems),
                listed(&residue),
            ])
        }
        "dot_pairs" => {
            let (rust, perl) = s::filter_dot_breakend_pairs(&tuples_of(&a[0]), &tuples_of(&a[1]));
            list([listed(&rust), listed(&perl)])
        }
        "non_integer" => {
            let (perl, rust) = s::filter_non_integer_end_rows(
                &tuples_of(&a[0]),
                &tuples_of(&a[1]),
                &variants(&a[2]),
                whole(&a[3]),
            );
            list([listed(&perl), listed(&rust)])
        }
        "model" => {
            let tx = pyjson::loads(text(&a[1])).expect("a record");
            model_answer(m::model_from_cache_record(text(&a[0]), &tx).as_ref())
        }
        "mate_side" => {
            let tx = pyjson::loads(text(&a[3])).expect("a record");
            match m::model_from_cache_record(text(&a[2]), &tx) {
                None => st("no model"),
                Some(model) => m::mate_side_derivation(text(&a[0]), text(&a[1]), &model).map_or(
                    Value::Null,
                    |(mate, local)| {
                        list([
                            list(mate.into_iter().map(st)),
                            list(local.into_iter().map(st)),
                        ])
                    },
                ),
            }
        }
        "cache" => authority_answer(Path::new(text(&a[0]))),
        "mate_pairs" => {
            let (_by_chr, _spans, models) = s::load_cache_authority(Path::new(text(&a[2])));
            let (rust, perl, by_kind) = s::filter_breakend_mate_local_read_pairs(
                &tuples_of(&a[0]),
                &tuples_of(&a[1]),
                &models,
            );
            list([
                listed(&rust),
                listed(&perl),
                Value::Dict(
                    by_kind
                        .into_iter()
                        .map(|(kind, pairs)| (kind, Value::Int(pairs)))
                        .collect(),
                ),
            ])
        }
        "own_rows" => {
            let (by_chr, spans, _models) = s::load_cache_authority(Path::new(text(&a[2])));
            listed(&s::filter_bnd_own_rows_lost(
                &tuples_of(&a[0]),
                &tuples_of(&a[1]),
                &by_chr,
                &spans,
            ))
        }
        "parse_vcf" => match s::parse_vcf(Path::new(text(&a[0]))) {
            Err(error) => dict(vec![(
                "oserror",
                st(vep_formats_compare::report::strerror(&error)),
            )]),
            Ok(read) => list(read.iter().map(|v| {
                list([
                    st(v.chrom.as_str()),
                    st(v.pos.to_string()),
                    st(v.vid.as_str()),
                    st(v.r#ref.as_str()),
                    st(v.alt.as_str()),
                    st(v.variant_type.as_str()),
                    st(v.location.as_str()),
                    st(v.info.as_str()),
                    opt(v.pos_text.clone()),
                ])
            })),
        },
        "side_table" => side_table(a)?,
        other => panic!("no function {other}"),
    })
}

fn canonical(value: &Value) -> String {
    pyjson::dumps(
        value,
        &pyjson::DumpOptions {
            indent: None,
            sort_keys: true,
            ensure_ascii: false,
        },
    )
}

fn cut(text: &str) -> String {
    text.chars().take(600).collect()
}

/// The cases of `path` whose answer differs from the stored one, of how many.
fn differing(path: &Path) -> (usize, usize) {
    static HOOK: std::sync::Once = std::sync::Once::new();
    HOOK.call_once(|| {
        let default_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if info.payload().downcast_ref::<Uncaught>().is_none() {
                default_hook(info);
            }
        }));
    });
    let file = std::fs::File::open(path).expect("the cases open");
    let mut compared: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut shown = 0;
    for (number, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line.expect("a line");
        let case = pyjson::loads(&line).expect("a case");
        let f = text(case.get("f").expect("f")).to_string();
        let a = items(case.get("a").expect("a"));
        let expected = canonical(case.get("x").expect("x"));
        let got = match catch_unwind(AssertUnwindSafe(|| answer(&f, a))) {
            Ok(Ok(value)) => value,
            Ok(Err(Error::Comparator(text))) => dict(vec![("error", st(text))]),
            Err(payload) => match payload.downcast_ref::<Uncaught>() {
                Some(Uncaught(what)) => dict(vec![("raises", st(what.as_str()))]),
                None => dict(vec![("panic", st("a panic of the tool's own"))]),
            },
        };
        let mut got = canonical(&got);
        // An ending the tool words in its own way.
        if expected == r#"{"raises": "*"}"# && got.starts_with(r#"{"raises": "#) {
            got = expected.clone();
        }
        let counts = compared.entry(f.clone()).or_default();
        counts.0 += 1;
        if got != expected {
            counts.1 += 1;
            if shown < 20 {
                shown += 1;
                eprintln!(
                    "line {} {f}\n  arguments {}\n  stored     {}\n  the tool   {}",
                    number + 1,
                    cut(&canonical(case.get("a").expect("a"))),
                    cut(&expected),
                    cut(&got)
                );
            }
        }
    }
    let mut cases = 0;
    let mut mismatches = 0;
    for (f, (asked, differ)) in &compared {
        if *differ > 0 {
            eprintln!("{f}: {differ} of {asked} cases differ");
        }
        cases += asked;
        mismatches += differ;
    }
    (mismatches, cases)
}

#[test]
fn the_kept_cases_answer_as_the_script_answers() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/structural/cases.jsonl");
    let (mismatches, cases) = differing(&path);
    assert!(cases >= 3000, "{cases} cases");
    assert_eq!(mismatches, 0, "of {cases} cases");
}

#[test]
#[ignore = "reads the file VEP_FORMATS_COMPARE_STRUCTURAL_CASES names"]
fn the_made_cases_answer_as_the_script_answers() {
    let path = std::env::var_os("VEP_FORMATS_COMPARE_STRUCTURAL_CASES")
        .expect("VEP_FORMATS_COMPARE_STRUCTURAL_CASES names the cases");
    let (mismatches, cases) = differing(Path::new(&path));
    assert_eq!(mismatches, 0, "of {cases} cases");
}
