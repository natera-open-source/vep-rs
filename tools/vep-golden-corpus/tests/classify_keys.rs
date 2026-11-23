// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up outputs whose keys differ: each row Ensembl VEP alone writes
//! is an entry of `divergences` under its class, each key vep-rs alone writes an entry of
//! `vep_rs_only_tuples`, a with-reference statement written by hand is carried onto the
//! entry of its key, the masks of a run leave out the reference's regulatory rows, and a
//! flag set's run is written into its own entry.

mod classify_support;
mod common;

use classify_support::{
    counted, divergence, key, obj, record, records, row, transcript, vep_rs_only, with_member,
    written, Case,
};
use common::{clean, failed, manifest_file, matrix, object, run};
use vep_golden_corpus::constants::{mask_classes, WITH_REFERENCE_SETS};
use vep_golden_corpus::manifest::with_reference_sets_error;

const SWAP_SET: &str = "start_lost,start_retained_variant";
const BOTH_CLASSES: &str = "{'start_cooccurrence_swap': 1, 'splice_family_swap': 1}";
const BOTH_SUMMARIES: (&str, &str) = (
    r#"{"start_cooccurrence_swap": 1, "splice_family_swap": 1}"#,
    "{}",
);

#[test]
fn classify_marks_documented_and_unexplained_divergences() {
    // Each divergence names the records its key resolves to; the record without an ID
    // resolves through the name of the raw line (`21_400_AT/A`); a key only vep-rs writes
    // for a record Ensembl VEP annotated is an `unexplained_residual`.
    let case = Case::new("marks");
    let held = records(&[
        record("21\t100\trs1\tA\tG\t.\tPASS\t."),
        record("21\t200\trs2\tA\tT\t.\tPASS\t."),
        record("21\t300\trs3\tA\tC\t.\tPASS\t."),
        record("21\t400\t.\tAT\tA\t.\tPASS\t."),
    ]);
    case.manifest(&obj(&[&held]));
    let shared = [
        row("rs1", "21:100", "G", "ENST1", "missense_variant"),
        row("rs2", "21:200", "T", "ENST2", "start_retained_variant"),
        row(
            "rs3",
            "21:300",
            "C",
            "ENST3",
            "intron_variant,splice_region_variant",
        ),
        row("21_400_AT/A", "21:401", "-", "ENST4", "frameshift_variant"),
        row("21_400_AT/A", "21:401", "-", "ENST5", "intron_variant"),
    ];
    case.outputs(
        &[
            row("rs1", "21:100", "G", "ENST1", "missense_variant"),
            row("rs2", "21:200", "T", "ENST2", SWAP_SET),
            row("rs3", "21:300", "C", "ENST3", "intron_variant"),
            row("21_400_AT/A", "21:401", "-", "ENST4", "frameshift_variant"),
        ]
        .concat(),
        &shared.concat(),
    );
    let ran = case.classify(&[]);
    let tuples = (1, "{'unexplained_residual': 1}");
    assert_eq!(ran, clean(0, &counted((2, BOTH_CLASSES), tuples, "[]")));
    let expected = written(
        &[
            divergence(
                &key("21:200", "T", transcript("ENST2"), &[1]),
                SWAP_SET,
                &["start_retained_variant"],
                "start_cooccurrence_swap",
            ),
            divergence(
                &key("21:300", "C", transcript("ENST3"), &[2]),
                "intron_variant",
                &["intron_variant,splice_region_variant"],
                "splice_family_swap",
            ),
        ],
        &[vep_rs_only(
            &key("21:401", "-", transcript("ENST5"), &[3]),
            &["intron_variant"],
            "unexplained_residual",
        )],
        (BOTH_SUMMARIES.0, r#"{"unexplained_residual": 1}"#),
    );
    assert_eq!(case.after(), obj(&[&held, &expected]));
}

/// A two-record corpus with the two outputs `classify` reads: both engines' rows differ on
/// each record, by a pair the set rules name. Returns the case and its `records` member.
fn reference_key_case(tag: &str) -> (Case, String) {
    let case = Case::new(tag);
    let held = records(&[
        record("21\t200\trs2\tA\tT\t.\tPASS\t."),
        record("21\t300\trs3\tA\tC\t.\tPASS\t."),
    ]);
    case.manifest(&obj(&[&held]));
    case.outputs(
        &(row("rs2", "21:200", "T", "ENST2", SWAP_SET)
            + &row("rs3", "21:300", "C", "ENST3", "intron_variant")),
        &(row("rs2", "21:200", "T", "ENST2", "start_retained_variant")
            + &row(
                "rs3",
                "21:300",
                "C",
                "ENST3",
                "intron_variant,splice_region_variant",
            )),
    );
    (case, held)
}

/// The entry of the first record of `reference_key_case`.
fn swap() -> String {
    divergence(
        &key("21:200", "T", transcript("ENST2"), &[0]),
        SWAP_SET,
        &["start_retained_variant"],
        "start_cooccurrence_swap",
    )
}

/// The entry of its second record.
fn splice() -> String {
    divergence(
        &key("21:300", "C", transcript("ENST3"), &[1]),
        "intron_variant",
        &["intron_variant,splice_region_variant"],
        "splice_family_swap",
    )
}

/// The with-reference member of an entry, stating the sets given as JSON.
fn stating(sets: &str) -> String {
    format!(r#""{WITH_REFERENCE_SETS}": {sets}"#)
}

#[test]
fn classify_carries_a_with_reference_key_through_a_regeneration() {
    // `classify` writes no with-reference key of its own, and a key written by hand stays:
    // the entry written for the same Location, Allele, Feature, records and Ensembl VEP set
    // carries it unchanged, as its last member. A flag set's own list carries its key the
    // same way, and a set without a list carries the keys of the manifest's own, the list
    // the golden harness holds such a set to. A run whose outputs give no such entry, or
    // whose regenerated sets the key repeats, stops with the reason and leaves the manifest
    // as it was.
    let (case, held) = reference_key_case("carries");
    let line = counted((2, BOTH_CLASSES), (0, "{}"), "[]");
    assert_eq!(case.classify(&[]), clean(0, &line));
    let first = written(&[swap(), splice()], &[], BOTH_SUMMARIES);
    assert_eq!(case.after(), obj(&[&held, &first]));

    let stated = with_member(&swap(), &stating(&format!(r#"["{SWAP_SET}"]"#)));
    let second = written(&[stated.clone(), splice()], &[], BOTH_SUMMARIES);
    case.manifest(&obj(&[&held, &second]));
    assert_eq!(case.classify(&[]), clean(0, &line));
    assert_eq!(case.after(), obj(&[&held, &second]));
    let carried = case.text();

    // A flag set without a list is held to the manifest's own, so `--flag-set` carries the
    // key that list states into the set's first list; the manifest's own stays.
    let set = |entries: &str| {
        format!(r#""flag_sets": [{{"name": "again", "flags": ["--again"]{entries}}}]"#)
    };
    let listless = obj(&[&held, &second, &set("")]);
    case.manifest(&listless);
    assert_eq!(case.classify(&["--flag-set", "again"]), clean(0, &line));
    let carried_in_set = set(&format!(", {second}"));
    assert_eq!(case.after(), obj(&[&held, &second, &carried_in_set]));
    // A set with a list of its own, an empty one too, carries its own keys alone.
    let empty = set(r#", "divergences": []"#);
    case.manifest(&obj(&[&held, &second, &empty]));
    assert_eq!(case.classify(&["--flag-set", "again"]), clean(0, &line));
    let unstated = set(&format!(", {first}"));
    assert_eq!(case.after(), obj(&[&held, &second, &unstated]));
    let in_set = written(
        &[
            swap(),
            with_member(&splice(), &stating(r#"["intron_variant"]"#)),
        ],
        &[],
        BOTH_SUMMARIES,
    );
    let stated_in_set = set(&format!(", {in_set}"));
    case.manifest(&obj(&[&held, &second, &stated_in_set]));
    assert_eq!(case.classify(&["--flag-set", "again"]), clean(0, &line));
    assert_eq!(case.after(), obj(&[&held, &second, &stated_in_set]));
    case.manifest(&carried);

    // The run's outputs agree on the key: no entry is written to carry the statement on.
    let vep_rs = std::fs::read_to_string(case.path("rs.txt")).unwrap();
    let column = |set: &str| format!("\t{set}\t");
    case.write(
        "agreeing.txt",
        &vep_rs.replace(&column("start_retained_variant"), &column(SWAP_SET)),
    );
    assert_eq!(
        case.classify_against("agreeing.txt", &[]),
        failed(&format!(
            "the entry on 21:200 T ENST2 carries {WITH_REFERENCE_SETS}, and this run writes no entry of its \
             Location, Allele, Feature, records and Ensembl VEP set to carry it on: remove the key, or classify the \
             outputs of a run without a reference sequence"
        ))
    );
    assert_eq!(case.text(), carried);
    // The same run of a set without a list is refused on the key the manifest's own states.
    case.manifest(&listless);
    let unwritten = case.text();
    assert_eq!(
        case.classify_against("agreeing.txt", &["--flag-set", "again"]),
        failed(&format!(
            "the entry on 21:200 T ENST2 carries {WITH_REFERENCE_SETS}, and this run writes no entry of its \
             Location, Allele, Feature, records and Ensembl VEP set to carry it on: remove the key, or classify the \
             outputs of a run without a reference sequence"
        ))
    );
    assert_eq!(case.text(), unwritten);

    // The regenerated sets are the stated ones: the key states no second outcome.
    let repeating = with_member(&swap(), &stating(r#"["start_lost"]"#));
    case.manifest(&obj(&[
        &held,
        &written(&[repeating, splice()], &[], BOTH_SUMMARIES),
    ]));
    let kept = case.text();
    case.write(
        "lost.txt",
        &vep_rs.replace(&column("start_retained_variant"), &column("start_lost")),
    );
    assert_eq!(
        case.classify_against("lost.txt", &[]),
        failed(&format!(
            "the entry on 21:200 T ENST2: {WITH_REFERENCE_SETS} repeats vep_rs_consequence_sets and states no second outcome"
        ))
    );
    assert_eq!(case.text(), kept);
}

#[test]
fn a_malformed_with_reference_key_is_refused_by_every_reader() {
    // The key is a non-empty sorted list of distinct consequence sets, other than the
    // entry's `vep_rs_consequence_sets`, that names Ensembl VEP's set alone or not at all,
    // on a `divergences` entry: any other shape is refused wherever a manifest is read, in
    // its own run's list and in a flag set's alike, before anything is written.
    let (case, held) = reference_key_case("malformed");
    let good = written(&[swap(), splice()], &[], BOTH_SUMMARIES);
    let line = counted((2, BOTH_CLASSES), (0, "{}"), "[]");
    assert_eq!(case.classify(&[]), clean(0, &line));
    assert_eq!(case.after(), obj(&[&held, &good]));
    let own = r#"["start_retained_variant"]"#;
    for value in [
        format!(r#"["{SWAP_SET}"]"#),
        r#"["intron_variant"]"#.to_string(),
        r#"["intron_variant", "start_lost"]"#.to_string(),
    ] {
        let stated = object(&with_member(&swap(), &stating(&value)));
        assert_eq!(with_reference_sets_error(&stated), None, "{value}");
    }
    assert_eq!(with_reference_sets_error(&object(&swap())), None);
    let not_a_list = "is not a non-empty list";
    let not_a_set = "holds an item that is not a sorted, comma-joined set of terms";
    let not_sorted = "is not a sorted list of distinct sets";
    let shapes = [
        (format!(r#""{SWAP_SET}""#), not_a_list),
        ("[]".to_string(), not_a_list),
        (
            format!(r#"{{"vep_rs_consequence_sets": ["{SWAP_SET}"]}}"#),
            not_a_list,
        ),
        (format!(r#"["{SWAP_SET}", null]"#), not_a_set),
        (
            r#"["start_retained_variant,start_lost"]"#.to_string(),
            not_a_set,
        ),
        (r#"[""]"#.to_string(), not_a_set),
        (
            r#"["start_lost", "intron_variant"]"#.to_string(),
            not_sorted,
        ),
        (
            r#"["intron_variant", "intron_variant"]"#.to_string(),
            not_sorted,
        ),
        (
            own.to_string(),
            "repeats vep_rs_consequence_sets and states no second outcome",
        ),
        (
            format!(r#"["intron_variant", "{SWAP_SET}"]"#),
            "names Ensembl VEP's set beside another",
        ),
    ];
    let path = manifest_file(&case.corpus);
    let path = path.to_str().unwrap();
    let matrix = matrix(&case.scratch.path, &["symbol"]);
    for (value, reason) in &shapes {
        let stated = with_member(&swap(), &stating(value));
        assert_eq!(
            with_reference_sets_error(&object(&stated)),
            Some(*reason),
            "{value}"
        );
        let in_the_manifest = written(&[stated.clone(), splice()], &[], BOTH_SUMMARIES);
        let in_a_set = format!(r#""flag_sets": [{{"name": "other", "divergences": [{stated}]}}]"#);
        for (owner, manifest) in [
            ("the manifest's", obj(&[&held, &in_the_manifest])),
            ("flag set other", obj(&[&held, &good, &in_a_set])),
        ] {
            case.manifest(&manifest);
            let refusal = failed(&format!(
                "{path} {owner} divergences entry on 21:200 T ENST2: {WITH_REFERENCE_SETS} {reason}"
            ));
            assert_eq!(case.classify(&[]), refusal, "{owner} {value}");
            assert_eq!(
                run(&["covers", &case.corpus, "--list", "--matrix", &matrix]),
                refusal,
                "{owner} {value}"
            );
            assert_eq!(case.text(), manifest);
        }
    }
    // The key belongs on a `divergences` entry: the other two lists state no Ensembl VEP
    // set.
    let misplaced = format!(
        r#"[{{"location": "21:200", "allele": "T", "feature": "ENST9", {}}}]"#,
        stating(&format!(r#"["{SWAP_SET}"]"#))
    );
    let in_tuples = format!(
        r#""divergences": [{}, {}], "vep_rs_only_tuples": {misplaced}, "divergence_summary": {}, "vep_rs_only_summary": {{}}"#,
        swap(),
        splice(),
        BOTH_SUMMARIES.0
    );
    let in_fields = format!(r#"{good}, "field_divergences": {misplaced}"#);
    for (section, members) in [
        ("vep_rs_only_tuples", in_tuples),
        ("field_divergences", in_fields),
    ] {
        case.manifest(&obj(&[&held, &members]));
        assert_eq!(
            case.classify(&[]),
            failed(&format!(
                "{path} the manifest's {section} entry on 21:200 T ENST9: {WITH_REFERENCE_SETS} belongs on a divergences entry"
            ))
        );
    }
}

#[test]
fn classify_leaves_a_reference_row_without_a_vep_rs_row_unexplained() {
    // A row Ensembl VEP writes on a key vep-rs writes no row for has no documented class,
    // whatever the shape of its allele: a symbolic allele's row, a breakend's and a
    // sequence allele's are each an `unexplained_residual` that lists no vep-rs set.
    let case = Case::new("unexplained");
    let held = records(&[
        record("21\t100\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=500"),
        record("21\t600\tbnd1\tA\tA[21:900[\t.\tPASS\tSVTYPE=BND"),
        record("21\t700\trs1\tA\tG\t.\tPASS\t."),
    ]);
    case.manifest(&obj(&[&held]));
    let on = |feature: &str, terms: [&str; 3]| {
        row("sv1", "21:101-500", "deletion", feature, terms[0])
            + &row("bnd1", "21:600", "A[21:900[", feature, terms[1])
            + &row("rs1", "21:700", "G", feature, terms[2])
    };
    let shared = on(
        "ENST1",
        [
            "transcript_ablation",
            "feature_truncation",
            "intron_variant",
        ],
    );
    let alone = on(
        "ENST2",
        ["intron_variant", "feature_truncation", "intron_variant"],
    );
    case.outputs(&(shared.clone() + &alone), &shared);
    let rows = (3, "{'unexplained_residual': 3}");
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted(rows, (0, "{}"), "[]"))
    );
    let unexplained = |location: &str, allele: &str, ordinal: usize, set: &str| {
        divergence(
            &key(location, allele, transcript("ENST2"), &[ordinal]),
            set,
            &[],
            "unexplained_residual",
        )
    };
    let expected = written(
        &[
            unexplained("21:101-500", "deletion", 0, "intron_variant"),
            unexplained("21:600", "A[21:900[", 1, "feature_truncation"),
            unexplained("21:700", "G", 2, "intron_variant"),
        ],
        &[],
        (r#"{"unexplained_residual": 3}"#, "{}"),
    );
    assert_eq!(case.after(), obj(&[&held, &expected]));
}

#[test]
fn classify_leaves_out_the_masked_regulatory_rows() {
    // Under `masks: ["regulatory_feature_rows"]` the reference's RegulatoryFeature and
    // MotifFeature rows are no keys of Ensembl VEP alone; without the mask they are
    // unexplained residuals, and an unknown class or a `masks` that is not a list of names
    // is refused.
    let regulatory = "rs1\t21:100\tG\t-\tENSR1\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n";
    let motif = "rs1\t21:100\tG\t-\tENSM1\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n";
    let transcript_row = row("rs1", "21:100", "G", "ENST1", "missense_variant");
    let held = records(&[record("21\t100\trs1\tA\tG\t.\tPASS\t.")]);
    let outputs = |case: &Case| {
        case.outputs(
            &(transcript_row.clone() + regulatory + motif),
            &transcript_row,
        );
    };
    let unmasked = Case::named("unmasked", "c0");
    outputs(&unmasked);
    unmasked.manifest(&obj(&[r#""masks": []"#, &held]));
    let rows = (2, "{'unexplained_residual': 2}");
    assert_eq!(
        unmasked.classify(&[]),
        clean(0, &counted(rows, (0, "{}"), "[]"))
    );
    let on = |feature: (&str, &str), set: &str| {
        divergence(
            &key("21:100", "G", feature, &[0]),
            set,
            &[],
            "unexplained_residual",
        )
    };
    let expected = written(
        &[
            on(("ENSM1", "MotifFeature"), "TF_binding_site_variant"),
            on(("ENSR1", "RegulatoryFeature"), "regulatory_region_variant"),
        ],
        &[],
        (r#"{"unexplained_residual": 2}"#, "{}"),
    );
    assert_eq!(unmasked.after(), obj(&[r#""masks": []"#, &held, &expected]));

    let masked = Case::named("masked", "c1");
    outputs(&masked);
    let masks = r#""masks": ["regulatory_feature_rows"]"#;
    masked.manifest(&obj(&[masks, &held]));
    assert_eq!(
        masked.classify(&[]),
        clean(0, &counted((0, "{}"), (0, "{}"), "[]"))
    );
    let nothing = written(&[], &[], ("{}", "{}"));
    assert_eq!(masked.after(), obj(&[masks, &held, &nothing]));

    let unknown = Case::named("unknown", "unknown");
    outputs(&unknown);
    let path = manifest_file(&unknown.corpus);
    let path = path.to_str().unwrap();
    unknown.manifest(r#"{"masks": ["splice_rows"], "records": []}"#);
    assert_eq!(
        unknown.classify(&[]),
        failed(&format!(
            "{path} masks names unknown classes ['splice_rows'] (the classes are {})",
            mask_classes().join(", ")
        ))
    );
    unknown.manifest(r#"{"masks": "regulatory_feature_rows", "records": []}"#);
    assert_eq!(
        unknown.classify(&[]),
        failed(&format!("{path} masks is not a list of names"))
    );
}

#[test]
fn classify_writes_a_flag_sets_divergences_into_its_entry() {
    // `classify --flag-set` documents that set's run: its divergence keys land in the
    // entry, the manifest's own stay as they were, and the default run's `--flag` option
    // is refused.
    let case = Case::new("into-its-entry");
    let held = records(&[record("21\t100\trs1\tA\tG\t.\tPASS\t.")]);
    let own = r#""covers": ["symbol"], "divergences": [{"location": "kept"}]"#;
    let set = |entries: &str| {
        format!(
            r#""flag_sets": [{{"name": "pick", "flags": ["--pick"], "covers": ["pick"]{entries}}}]"#
        )
    };
    case.manifest(&obj(&[&held, own, &set("")]));
    let before = case.text();
    case.outputs(
        &row("rs1", "21:100", "G", "ENST1", SWAP_SET),
        &row("rs1", "21:100", "G", "ENST1", "start_retained_variant"),
    );
    assert_eq!(
        case.classify(&["--flag-set", "pick", "--flag=--pick"]),
        failed(
            "--flag and --fasta-name describe the default run; a flag set's flags are its own entry's"
        )
    );
    assert_eq!(
        case.classify(&["--flag-set", "none"]),
        failed("the manifest has no flag set 'none'")
    );
    assert_eq!(case.text(), before);
    let ran = case.classify(&["--flag-set", "pick"]);
    let rows = (1, "{'start_cooccurrence_swap': 1}");
    assert_eq!(ran, clean(0, &counted(rows, (0, "{}"), "[]")));
    let entries = written(
        &[divergence(
            &key("21:100", "G", transcript("ENST1"), &[0]),
            SWAP_SET,
            &["start_retained_variant"],
            "start_cooccurrence_swap",
        )],
        &[],
        (r#"{"start_cooccurrence_swap": 1}"#, "{}"),
    );
    assert_eq!(
        case.after(),
        obj(&[&held, own, &set(&format!(", {entries}"))])
    );
}
