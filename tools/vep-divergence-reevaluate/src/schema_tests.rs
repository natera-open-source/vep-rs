// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::pyjson::loads;
use crate::testing::{repro, texts, DELETE, IMAGE, SIGHTING_BASIS, UNOBSERVED_BASIS};

fn validate(edits: &[(&str, &str)]) -> Vec<String> {
    Schema::new().validate(&repro(edits))
}

fn has(errors: &[String], line: &str) -> bool {
    errors.iter().any(|error| error == line)
}

/// A TOML basic string of a text.
fn quoted(text: &str) -> String {
    let mut out = String::new();
    crate::pyfmt::push_json_str(&mut out, text);
    out
}

#[test]
fn validate_accepts_the_committed_shape_and_names_each_violation() {
    assert!(IMAGE.contains("@sha256:"));
    assert_eq!(validate(&[]), Vec::<String>::new());
    assert_eq!(
        validate(&[("entry", "'other'")]),
        ["entry is 'other', the file is x.toml"]
    );
    let cases: [(&[(&str, &str)], &str); 33] = [
        (
            &[("image", "'ensemblorg/ensembl-vep:release_116.2'")],
            "image is not pinned by digest",
        ),
        (
            &[("input.text", "'a'"), ("input.path", "'b'")],
            "input needs exactly one of text and path",
        ),
        (
            &[("input.text", DELETE)],
            "input needs exactly one of text and path",
        ),
        (&[("cache.assembly", DELETE)], "cache.assembly is missing"),
        (
            &[("cache.json_cache", DELETE)],
            "cache needs json_cache or gxf",
        ),
        (
            &[("cache.source", "'other'")],
            "cache.source is not ensembl, merged or refseq",
        ),
        (&[("cache.prepare", "['x']")], "unknown key cache.prepare"),
        (&[("colour", "'x'")], "unknown key colour"),
        (
            &[("compare.field", DELETE)],
            "compare.kind field needs compare.field",
        ),
        (
            &[("compare.kind", "'header'"), ("compare.field", DELETE)],
            "compare.kind header needs compare.pattern",
        ),
        (
            &[("compare.kind", "'json'"), ("compare.field", DELETE)],
            "compare.kind json needs compare.path",
        ),
        (
            &[("compare.kind", "'stats'"), ("compare.field", DELETE)],
            "compare.kind stats needs compare.section",
        ),
        (
            &[("compare.where", "{colour = 'x'}")],
            "compare.where has keys colour that compare.kind field does not read",
        ),
        (
            &[
                ("compare.kind", "'rows'"),
                ("compare.field", DELETE),
                ("compare.where", "{id = 'x'}"),
            ],
            "compare.where has keys id that compare.kind rows does not read",
        ),
        (
            &[("compare.kind", "'vcf'"), ("compare.where", "{gene = 'x'}")],
            "compare.where has keys gene that compare.kind vcf does not read",
        ),
        (
            &[("compare.kind", "'exit'"), ("compare.field", DELETE)],
            "compare.where has keys record that compare.kind exit does not read",
        ),
        (
            &[
                ("compare.kind", "'vcf'"),
                ("compare.field", DELETE),
                ("compare.line", "true"),
                ("compare.where", "{id = 'r', feature = 'ENST1'}"),
            ],
            "compare.where has keys feature that compare.kind vcf does not read",
        ),
        (
            &[
                ("expected.ensembl_vep", "'a'"),
                ("expected.vep_rs", "['a']"),
            ],
            "expected.ensembl_vep equals expected.vep_rs: the entry records no divergence",
        ),
        (
            &[("expected.ensembl_vep", "'a'"), ("expected.vep_rs", DELETE)],
            "expected needs ensembl_vep and vep_rs",
        ),
        (&[("run.tool", "'filter'")], "a tool run needs run.args"),
        (&[("prepare", "[1]")], "prepare is not a list of strings"),
        (
            &[("compare.field", DELETE), ("compare.line", "true")],
            "compare.kind field needs compare.field",
        ),
        (&[("image", "5")], "image is missing"),
        (
            &[("run.flags", "'--hgvs'")],
            "run.flags is not a list of strings",
        ),
        (
            &[("run.flags", "['--distance', 500]")],
            "run.flags is not a list of strings",
        ),
        (
            &[("cache.chromosomes", "[21]")],
            "cache.chromosomes is not a list of strings",
        ),
        (&[("prepare", "'true'")], "prepare is not a list of strings"),
        (
            &[
                ("expected.ensembl_vep", "'a'"),
                ("expected.vep_rs", "'b'"),
                ("expected.ensembl_vep_exit", "'2'"),
            ],
            "expected.ensembl_vep_exit is not an integer",
        ),
        (
            &[("cache.json_cache", "5")],
            "cache.json_cache is not a string",
        ),
        (&[("input.text", "5")], "input.text is not a string"),
        (
            &[("compare.where", "'rec'")],
            "compare.where is not a table",
        ),
        (
            &[("run.perl_mode", "'remote'")],
            "run.perl_mode is not offline or database",
        ),
        (
            &[("run.tools", "'samtools'")],
            "run.tools is not a list of strings",
        ),
    ];
    for (edits, line) in cases {
        let errors = validate(edits);
        assert!(has(&errors, line), "{edits:?}: {errors:?} lacks {line:?}");
    }
    let first_starts: [(&[(&str, &str)], &str); 6] = [
        (
            &[("compare.kind", "'flow'")],
            "compare.kind 'flow' is not one of",
        ),
        (&[("run.format", "'csv'")], "run.format 'csv' is not one of"),
        (&[("run.tool", "'grep'")], "run.tool 'grep' is not one of"),
        (
            &[
                ("compare.kind", "'lines'"),
                ("compare.field", DELETE),
                ("compare.where", DELETE),
                ("compare.pattern", "'('"),
            ],
            "compare.pattern is not a regular expression",
        ),
        (
            &[("compare.where", "{pattern = '['}")],
            "compare.where.pattern is not a regular expression",
        ),
        (
            &[("compare.mask", "['ok', '(?P<x']")],
            "compare.mask is not a regular expression",
        ),
    ];
    for (edits, start) in first_starts {
        let errors = validate(edits);
        assert!(
            errors[0].contains(start),
            "{edits:?}: {errors:?} does not open on {start:?}"
        );
    }
    let valid: [&[(&str, &str)]; 5] = [
        &[
            ("cache.json_cache", DELETE),
            ("cache.gxf", "'{corpus}/genes.gtf.gz'"),
        ],
        &[
            ("compare.kind", "'json'"),
            ("compare.field", DELETE),
            ("compare.path", "'x'"),
            ("compare.where", "{input = 'l', id = 'r'}"),
        ],
        &[
            ("run.tool", "'filter'"),
            ("run.args", "['-i', '{input}']"),
            ("cache.assembly", DELETE),
            ("cache.json_cache", DELETE),
        ],
        &[
            ("compare.kind", "'vcf'"),
            ("compare.field", DELETE),
            ("compare.line", "true"),
            ("compare.where", "{id = 'r'}"),
        ],
        &[
            ("expected.ensembl_vep", "'a'"),
            ("expected.vep_rs", "'b'"),
            ("expected.ensembl_vep_exit", "2"),
        ],
    ];
    for edits in valid {
        assert_eq!(validate(edits), Vec::<String>::new(), "{edits:?}");
    }
    assert_eq!(
        validate(&[("prepare", "5")]),
        ["prepare is not a list of strings"]
    );
    assert_eq!(
        validate(&[("input.format", "['vcf']")]),
        ["input.format is not a string"]
    );
    assert_eq!(
        validate(&[("run.tool", "['filter']"), ("compare.kind", "['rows']")]),
        ["run.tool is not a string", "compare.kind is not a string"]
    );
    assert_eq!(validate(&[("cache", "'x'")]), ["cache is not a table"]);
    assert_eq!(validate(&[("compare", "3")]), ["compare is not a table"]);
}

#[test]
fn validate_holds_the_input_class_to_three_values_and_its_basis_to_the_class() {
    let missing = ["input_class is missing", "input_class_basis is missing"];
    assert_eq!(validate(&[("input_class", DELETE)]), missing[..1]);
    assert_eq!(validate(&[("input_class_basis", DELETE)]), missing[1..]);
    assert_eq!(
        validate(&[("input_class", DELETE), ("input_class_basis", DELETE)]),
        missing
    );
    assert_eq!(
        validate(&[("input_class", "''"), ("input_class_basis", "' '")]),
        missing
    );
    assert_eq!(
        validate(&[("input_class", "3"), ("input_class_basis", "['a']")]),
        missing
    );
    for fourth in ["synthetic", "Malformed", "tier 3"] {
        assert_eq!(
            validate(&[("input_class", &quoted(fourth))]),
            [format!(
                "input_class '{fourth}' is not one of observed, unobserved, malformed"
            )]
        );
    }
    let one_line = ["input_class_basis is not one line"];
    for parted in ["\n", "\r", "\u{2028}"] {
        let basis = format!("well-formed under VCF 4.2;{parted}not seen in the measured datasets");
        assert_eq!(
            validate(&[("input_class_basis", &quoted(&basis))]),
            one_line,
            "{parted:?}"
        );
    }
    assert_eq!(
        validate(&[(
            "input_class_basis",
            &quoted(&format!("{UNOBSERVED_BASIS}\n"))
        )]),
        one_line
    );
    let unobserved = ["input_class_basis of an unobserved entry does not read `well-formed under <format or option contract>; not seen in the measured datasets`"];
    for wrong in [
        "allowed by VCF 4.2",
        "well-formed under VCF 4.2",
        "allowed by VCF 4.2; not seen in the measured datasets",
        "well-formed under ; not seen in the measured datasets",
        "well-formed under \t ; not seen in the measured datasets",
    ] {
        assert_eq!(
            validate(&[("input_class_basis", &quoted(wrong))]),
            unobserved,
            "{wrong}"
        );
    }
    let longer = format!("{UNOBSERVED_BASIS}; classed by its most reachable shape, a bracket mate: a mate position that is no number breaks the VCF 4.2 specification, section 5.4");
    assert_eq!(
        validate(&[("input_class_basis", &quoted(&longer))]),
        Vec::<String>::new()
    );
    let observed = ["input_class_basis of an observed entry is neither a tier 1 or tier 2 `Tier and count` cell nor `seen in <the real data its page names>: <the figure its page measures there>`"];
    for right in [
        "1 · 7.3 per million · no flag needed",
        "2 · 156.8 per million · 0 flag sets, 28 corpora",
        SIGHTING_BASIS,
        "seen in a file (region 21:1-100) of a population dataset: 4,158 pairs: two of them on one record",
    ] {
        assert_eq!(
            validate(&[
                ("input_class", "'observed'"),
                ("input_class_basis", &quoted(right))
            ]),
            Vec::<String>::new(),
            "{right}"
        );
    }
    // The words `seen in`, the data name and the figure the page measures there, parted by
    // the first colon a space follows.
    for wrong in [
        "3 · none on the cuts · 9 flag sets, 2 corpora",
        "4,158 pairs in a population dataset",
        "a population dataset: 4,158 pairs",
        "Seen in a population dataset: 4,158 pairs",
        "seen in a population dataset",
        "seen in : 4,158 pairs",
        "seen in : 4,158 pairs: two of them on one record",
        "seen in :4,158 pairs: two of them on one record",
        "seen in  a population dataset: 4,158 pairs",
        "seen in a population dataset:  4,158 pairs",
        "seen in a population dataset: ",
        "seen in a population dataset:4,158 pairs",
        "seen in a population dataset : 4,158 pairs",
        "seen in a population dataset : 4,158 pairs: two of them on one record",
    ] {
        assert_eq!(
            validate(&[
                ("input_class", "'observed'"),
                ("input_class_basis", &quoted(wrong))
            ]),
            observed,
            "{wrong}"
        );
    }
    let sighting = Schema::new()
        .sighting
        .match_start("seen in a: b: c")
        .expect("a sighting");
    assert_eq!(sighting.name("where").as_deref(), Some("a"));
    assert_eq!(sighting.name("what").as_deref(), Some("b: c"));
    let rule = "POS is an Integer (VCF 4.2 specification, section 1.4.1 \"Fixed fields\"): the POS is text that is no number";
    assert_eq!(
        validate(&[
            ("input_class", "'malformed'"),
            ("input_class_basis", &quoted(rule))
        ]),
        Vec::<String>::new()
    );
    // A malformed entry's basis is the rule its input breaks, so the basis of another class
    // left under it is refused.
    for other in [
        UNOBSERVED_BASIS,
        "well-formed under VCF 4.2",
        SIGHTING_BASIS,
        "3 · none on the cuts · 9 flag sets, 2 corpora",
        "2 · 156.8 per million · 0 flag sets, 28 corpora",
    ] {
        assert_eq!(
            validate(&[
                ("input_class", "'malformed'"),
                ("input_class_basis", &quoted(other))
            ]),
            ["input_class_basis of a malformed entry reads as the basis of another class, not as the rule the input breaks"],
            "{other}"
        );
    }
}

/// A record may hold no `implementation` table, since a run of a fixture record needs
/// none; a table it holds names code or says why it names none, and names a test.
#[test]
fn validate_holds_an_implementation_table_to_its_shape() {
    let code = "{ code = ['crates/a/src/b.rs::Type::method'], tests = ['crates/a/tests/c.rs::pins_it'], rules = ['tools/t/src/r.rs::rule'] }";
    assert_eq!(validate(&[("implementation", code)]), Vec::<String>::new());
    let none =
        "{ none = 'The row is the ordinary one.', tests = ['crates/a/tests/c.rs::pins_it'] }";
    assert_eq!(validate(&[("implementation", none)]), Vec::<String>::new());
    assert_eq!(
        validate(&[(
            "implementation",
            "{ none = 'The row is the ordinary one.' }"
        )]),
        ["implementation needs tests"]
    );
    assert_eq!(
        validate(&[(
            "implementation",
            "{ none = 'The row is the ordinary one.', tests = [] }"
        )]),
        ["implementation needs tests"]
    );
    assert_eq!(
        validate(&[(
            "implementation",
            "{ tests = ['crates/a/tests/c.rs::pins_it'] }"
        )]),
        ["implementation needs code or none"]
    );
    assert_eq!(
        validate(&[("implementation", "{}")]),
        [
            "implementation needs code or none",
            "implementation needs tests"
        ]
    );
    assert_eq!(
        validate(&[
            ("implementation", code),
            ("implementation.none", "'The row is the ordinary one.'")
        ]),
        ["implementation holds code and none: none is for an entry with no code of its own"]
    );
    for row in [
        "crates/a/src/b.rs",
        "crates/a/src/b.rs::",
        "::method",
        "/crates/a/src/b.rs::method",
        "crates/../a.rs::method",
        "crates//a.rs::method",
        "crates/a/src/b.rs::Type::",
        "crates/a/src/b.rs::Type<T>::method",
        "crates/a/src/b.rs::a method",
    ] {
        assert_eq!(
            validate(&[
                ("implementation", code),
                ("implementation.rules", &format!("[{}]", quoted(row)))
            ]),
            [format!(
                "implementation.rules names {}, which is not <path>::<item>",
                crate::pyfmt::py_repr_str(row)
            )],
            "{row}"
        );
    }
    assert_eq!(
        validate(&[("implementation", code), ("implementation.none", "' '")]),
        [
            "implementation.none is not one line of text",
            "implementation holds code and none: none is for an entry with no code of its own"
        ]
    );
    assert_eq!(
        validate(&[
            ("implementation", none),
            ("implementation.none", "\"first\\nsecond\"")
        ]),
        ["implementation.none is not one line of text"]
    );
    assert_eq!(
        validate(&[("implementation", "'crates/a/src/b.rs::method'")]),
        ["implementation is not a table"]
    );
    assert_eq!(
        validate(&[
            ("implementation", code),
            ("implementation.code", "'crates/a/src/b.rs::method'"),
            ("implementation.none", "3"),
            ("implementation.files", "['crates/a/src/b.rs']")
        ]),
        [
            "unknown key implementation.files",
            "implementation.none is not a string",
            "implementation.code is not a list of strings"
        ]
    );
}

/// `tests/fixtures/validate_cases.jsonl` holds records edited from the valid one, each
/// with every violation in order: all of a record's lines are held, not one of them.
#[test]
fn every_violation_of_an_edited_record_is_listed_in_order() {
    let table = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validate_cases.jsonl"
    ))
    .expect("the case table");
    let mut cases = 0;
    for line in table.lines() {
        let case = loads(line).expect("a JSON object a line");
        let Some(Value::List(edits)) = case.get("edits") else {
            panic!("no edits: {line}");
        };
        let edits: Vec<(String, String)> = edits
            .iter()
            .map(|edit| {
                let pair = texts(edit);
                (pair[0].clone(), pair[1].clone())
            })
            .collect();
        let borrowed: Vec<(&str, &str)> = edits
            .iter()
            .map(|(path, literal)| (path.as_str(), literal.as_str()))
            .collect();
        assert_eq!(
            validate(&borrowed),
            texts(case.get("errors").expect("errors")),
            "{edits:?}"
        );
        cases += 1;
    }
    assert!(cases > 200, "the table holds {cases} cases");
}

/// A further run is the record with some of its tables replaced, and is held to the
/// schema as the record's own run is, each violation named by the run's place.
#[test]
fn further_runs_are_held_to_the_schema_as_the_records_own_run_is() {
    const EXPECTED: &str = "expected = { ensembl_vep = ['a'], vep_rs = [] }";
    let also = |tables: &[&str]| -> Vec<String> {
        let list: Vec<String> = tables
            .iter()
            .map(|table| format!("{{ {table} }}"))
            .collect();
        validate(&[("also", &format!("[{}]", list.join(", ")))])
    };
    assert_eq!(
        also(&[&format!(
            "run = {{ flags = ['--hgvsg'] }}, notes = 'n', {EXPECTED}"
        )]),
        Vec::<String>::new()
    );
    assert_eq!(also(&[]), Vec::<String>::new());
    for literal in ["'x'", "['x']", "{ run = { flags = [] } }"] {
        assert_eq!(
            validate(&[("also", literal)]),
            ["also is not a list of tables"],
            "{literal}"
        );
    }
    assert_eq!(
        also(&["run = { flags = ['--hgvsg'] }"]),
        ["also[1] needs expected"]
    );
    assert_eq!(
        also(&[EXPECTED, &format!("colour = 'x', notes = 1, {EXPECTED}")]),
        [
            "unknown key also[2].colour",
            "also[2].notes is not a string"
        ]
    );
    // The run's own tables are read in place of the record's.
    assert_eq!(
        also(&[
            EXPECTED,
            "expected = { ensembl_vep = ['a'], vep_rs = ['a'] }"
        ]),
        ["also[2]: expected.ensembl_vep equals expected.vep_rs: the entry records no divergence"]
    );
    assert_eq!(
        also(&[&format!(
            "compare = {{ kind = 'rows', colour = 1 }}, input = {{ text = 'a', path = 'b' }}, {EXPECTED}"
        )]),
        [
            "also[1]: unknown key compare.colour",
            "also[1]: input needs exactly one of text and path"
        ]
    );
    // A record whose own tables fail lists those violations alone.
    assert_eq!(
        validate(&[("colour", "'x'"), ("also", "'x'")]),
        ["unknown key colour"]
    );
}
