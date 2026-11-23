// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field` on made-up outputs whose keys agree and whose fields differ:
//! each differing field of a key is an entry of `field_divergences` under the class that
//! takes it, the masks of a run decide what differs, and the run's options travel with
//! the manifest.

mod classify_support;
mod common;

use classify_support::{
    counted, divergence, feature_row, field_divergence, fields_written, key, nothing_written, obj,
    record, records, row_with, transcript, with_member, written, Case,
};
use common::clean;

/// The line of standard output of a run with `--field`.
fn counted_with_fields(rows: (usize, &str), fields: (usize, &str)) -> String {
    let line = counted(rows, (0, "{}"), "[]");
    format!(
        "{}; field divergences {} ({})\n",
        line.strip_suffix('\n').unwrap(),
        fields.0,
        fields.1
    )
}

const UNEXPLAINED: &str = "unexplained_residual";

#[test]
fn classify_names_the_motif_score_and_cell_type_field_classes() {
    // A `MOTIF_SCORE_CHANGE` difference on a reverse-strand motif row whose allele spans two
    // or more bases is `motif_score_reverse_strand_mnv`, an empty Ensembl VEP `CELL_TYPE`
    // beside a vep-rs activity under `--cell_type` is `cell_type_comma_names`, and every
    // other field difference stays an `unexplained_residual`.
    let case = Case::new("motif");
    let held = records(&[
        record("21\t100\trs_rev\tCT\tAG\t.\tPASS\t."),
        record("21\t200\trs_fwd\tCT\tAG\t.\tPASS\t."),
        record("21\t300\trs_snv\tA\tG\t.\tPASS\t."),
        record("21\t400\trs_rf\tA\tG\t.\tPASS\t."),
    ]);
    let comma = r#"{"name": "comma", "flags": ["--regulatory", "--cell_type", "A549_(m,_58_y)"]}"#;
    case.manifest(&obj(&[&held, &format!(r#""flag_sets": [{comma}]"#)]));
    let motif = |strand: &str, score: &str| {
        format!("IMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM1;MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE={score}")
    };
    let binding = "TF_binding_site_variant";
    let region = "regulatory_region_variant";
    let regulatory = ("ENSR1", "RegulatoryFeature");
    let activity = "A549_(m,_58_y):ACTIVE";
    let rows = |scores: [&str; 3], cell_type: &str| {
        feature_row(
            "rs_rev",
            "21:100-101",
            "AG",
            ("ENSM1", "MotifFeature"),
            binding,
            &motif("-1", scores[0]),
        ) + &feature_row(
            "rs_fwd",
            "21:200-201",
            "AG",
            ("ENSM2", "MotifFeature"),
            binding,
            &motif("1", scores[1]),
        ) + &feature_row(
            "rs_snv",
            "21:300",
            "G",
            ("ENSM3", "MotifFeature"),
            binding,
            &motif("-1", scores[2]),
        ) + &feature_row(
            "rs_rf",
            "21:400",
            "G",
            regulatory,
            region,
            &format!("IMPACT=MODIFIER;CELL_TYPE={cell_type}"),
        )
    };
    case.outputs(
        &rows(["-0.092", "-0.016", "-0.089"], ""),
        &rows(["-0.136", "-0.017", "-0.090"], activity),
    );
    let fields = ["--field", "MOTIF_SCORE_CHANGE", "--field", "CELL_TYPE"];
    let compared = ["MOTIF_SCORE_CHANGE", "CELL_TYPE"];
    let score =
        |location: &str, feature: &str, ordinal: usize, values: (&str, &str), class: &str| {
            field_divergence(
                &key(location, "AG", (feature, "MotifFeature"), &[ordinal]),
                binding,
                ("MOTIF_SCORE_CHANGE", values.0, values.1),
                class,
            )
        };
    let scores = [
        score(
            "21:100-101",
            "ENSM1",
            0,
            ("-0.092", "-0.136"),
            "motif_score_reverse_strand_mnv",
        ),
        score("21:200-201", "ENSM2", 1, ("-0.016", "-0.017"), UNEXPLAINED),
        field_divergence(
            &key("21:300", "G", ("ENSM3", "MotifFeature"), &[2]),
            binding,
            ("MOTIF_SCORE_CHANGE", "-0.089", "-0.090"),
            UNEXPLAINED,
        ),
    ];
    let cell_type = |class: &str| {
        field_divergence(
            &key("21:400", "G", regulatory, &[3]),
            region,
            ("CELL_TYPE", "", activity),
            class,
        )
    };
    let entries = |class: &str| {
        let mut entries = scores.to_vec();
        entries.push(cell_type(class));
        entries
    };
    let under_the_flag = fields_written(
        &compared,
        &entries("cell_type_comma_names"),
        r#"{"motif_score_reverse_strand_mnv": 1, "unexplained_residual": 2, "cell_type_comma_names": 1}"#,
    );
    let mut args = fields.to_vec();
    args.extend(["--flag-set", "comma"]);
    assert_eq!(
        case.classify(&args),
        clean(
            0,
            &counted_with_fields(
                (0, "{}"),
                (
                    4,
                    "{'motif_score_reverse_strand_mnv': 1, 'unexplained_residual': 2, 'cell_type_comma_names': 1}"
                )
            )
        )
    );
    let set = with_member(comma, &format!("{}, {under_the_flag}", nothing_written()));
    let sets = format!(r#""flag_sets": [{set}]"#);
    assert_eq!(case.after(), obj(&[&held, &sets]));

    // Without --cell_type among the run's flags the empty reference value is not the
    // comma class.
    let mut args = fields.to_vec();
    args.push("--flag=--regulatory");
    assert_eq!(
        case.classify(&args),
        clean(
            0,
            &counted_with_fields(
                (0, "{}"),
                (
                    4,
                    "{'motif_score_reverse_strand_mnv': 1, 'unexplained_residual': 3}"
                )
            )
        )
    );
    let without_the_flag = fields_written(
        &compared,
        &entries(UNEXPLAINED),
        r#"{"motif_score_reverse_strand_mnv": 1, "unexplained_residual": 3}"#,
    );
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &sets,
            &nothing_written(),
            &without_the_flag,
            r#""flags": ["--regulatory"]"#
        ])
    );
}

#[test]
fn classify_names_the_rejoined_record_plugin_field_class() {
    // A covered plugin's field that differs on a `MINIMISED` row of a multi-allelic record
    // under `--minimal` is `minimal_plugin_sees_the_rejoined_record`; the same difference
    // on a bi-allelic record, on a multi-allelic record the flag left unsplit (no
    // `MINIMISED` on its rows), under a run without `--minimal`, or in a field no covered
    // plugin writes stays an `unexplained_residual`.
    let case = Case::new("rejoined");
    let held = records(&[
        record("21\t100\trs_multi\tC\tCCA,T\t.\tPASS\t."),
        record("21\t200\trs_bi\tC\tT\t.\tPASS\t."),
        record("21\t300\trs_unsplit\tG\tT,C\t.\tPASS\t."),
    ]);
    let minimal = r#"{"name": "minimal", "flags": ["--plugin", "gnomADc,cov.tsv.gz", "--minimal", "--allele_number"]}"#;
    let plain =
        r#"{"name": "plain", "flags": ["--plugin", "gnomADc,cov.tsv.gz", "--allele_number"]}"#;
    case.manifest(&obj(&[
        &held,
        &format!(r#""flag_sets": [{minimal}, {plain}]"#),
    ]));
    let intron = "intron_variant";
    let on = |name: &str, location: &str, allele: &str, extra: &str| {
        feature_row(name, location, allele, transcript("ENST1"), intron, extra)
    };
    case.outputs(
        &(on(
            "rs_multi",
            "21:100",
            "CA",
            "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=1.7938;NOTE=x",
        ) + &on(
            "rs_multi",
            "21:100",
            "T",
            "ALLELE_NUM=2;MINIMISED=1;gnomAD_mean_cov=1.7938;NOTE=x",
        ) + &on(
            "rs_bi",
            "21:200",
            "T",
            "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=2.0000;NOTE=x",
        ) + &on(
            "rs_unsplit",
            "21:300",
            "T",
            "ALLELE_NUM=1;gnomAD_mean_cov=3.0000;NOTE=x",
        )),
        &(on(
            "rs_multi",
            "21:100",
            "CA",
            "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=1.7937;NOTE=y",
        ) + &on(
            "rs_multi",
            "21:100",
            "T",
            "ALLELE_NUM=2;MINIMISED=1;gnomAD_mean_cov=1.7938;NOTE=x",
        ) + &on(
            "rs_bi",
            "21:200",
            "T",
            "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=2.5000;NOTE=x",
        ) + &on(
            "rs_unsplit",
            "21:300",
            "T",
            "ALLELE_NUM=1;gnomAD_mean_cov=3.5000;NOTE=x",
        )),
    );
    let fields = ["--field", "gnomAD_mean_cov", "--field", "NOTE"];
    let compared = ["gnomAD_mean_cov", "NOTE"];
    let differing =
        |location: &str, allele: &str, ordinal: usize, field: (&str, &str, &str), class: &str| {
            field_divergence(
                &key(location, allele, transcript("ENST1"), &[ordinal]),
                intron,
                field,
                class,
            )
        };
    let entries = |split: &str| {
        [
            differing(
                "21:100",
                "CA",
                0,
                ("gnomAD_mean_cov", "1.7938", "1.7937"),
                split,
            ),
            differing("21:100", "CA", 0, ("NOTE", "x", "y"), UNEXPLAINED),
            differing(
                "21:200",
                "T",
                1,
                ("gnomAD_mean_cov", "2.0000", "2.5000"),
                UNEXPLAINED,
            ),
            differing(
                "21:300",
                "T",
                2,
                ("gnomAD_mean_cov", "3.0000", "3.5000"),
                UNEXPLAINED,
            ),
        ]
    };
    let mut args = fields.to_vec();
    args.extend(["--flag-set", "minimal"]);
    assert_eq!(
        case.classify(&args),
        clean(
            0,
            &counted_with_fields(
                (0, "{}"),
                (
                    4,
                    "{'minimal_plugin_sees_the_rejoined_record': 1, 'unexplained_residual': 3}"
                )
            )
        )
    );
    let split = fields_written(
        &compared,
        &entries("minimal_plugin_sees_the_rejoined_record"),
        r#"{"minimal_plugin_sees_the_rejoined_record": 1, "unexplained_residual": 3}"#,
    );
    let minimal_after = with_member(minimal, &format!("{}, {split}", nothing_written()));
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &format!(r#""flag_sets": [{minimal_after}, {plain}]"#)
        ])
    );
    let mut args = fields.to_vec();
    args.extend(["--flag-set", "plain"]);
    assert_eq!(
        case.classify(&args),
        clean(
            0,
            &counted_with_fields((0, "{}"), (4, "{'unexplained_residual': 4}"))
        )
    );
    let unsplit = fields_written(
        &compared,
        &entries(UNEXPLAINED),
        r#"{"unexplained_residual": 4}"#,
    );
    let plain_after = with_member(plain, &format!("{}, {unsplit}", nothing_written()));
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &format!(r#""flag_sets": [{minimal_after}, {plain_after}]"#)
        ])
    );
}

#[test]
fn classify_applies_the_runs_own_masks() {
    // `classify` masks a run as the golden harness does: under `order_only` the `DOMAINS`
    // values of `--field` agree as multisets of their items while an item that differs is
    // still a field divergence, and a flag set's own `masks` (an empty list included)
    // replace the manifest's for `classify --flag-set`, the regulatory class among them.
    let case = Case::new("masks");
    let with_domains = |domains: &str| {
        row_with(
            "rs1",
            "21:100",
            "G",
            "ENST1",
            "missense_variant",
            &format!("IMPACT=MODERATE;DOMAINS={domains}"),
        )
    };
    let regulatory = "rs1\t21:100\tG\t-\tENSR1\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n";
    let listed = "Pfam:PF00059,PANTHER:PTHR14789";
    let reordered = "PANTHER:PTHR14789,Pfam:PF00059";
    let changed = "PANTHER:PTHR14789,Pfam:PF00060";
    case.write("vep.txt", &(with_domains(listed) + regulatory));
    case.write("reordered.txt", &with_domains(reordered));
    case.write("changed.txt", &with_domains(changed));
    let held = records(&[record("21\t100\trs1\tA\tG\t.\tPASS\t.")]);
    let masks = r#""masks": ["order_only"]"#;
    let plain = r#"{"name": "plain", "masks": []}"#;
    let masked = r#"{"name": "regulatory", "masks": ["regulatory_feature_rows"]}"#;
    let inherits = r#"{"name": "inherits"}"#;
    case.manifest(&obj(&[
        masks,
        &held,
        &format!(r#""flag_sets": [{plain}, {masked}, {inherits}]"#),
    ]));
    let one = (1, "{'unexplained_residual': 1}");
    let none = (0, "{}");
    let classify = |vep_rs: &str, flag_set: Option<&str>| {
        let mut args = vec!["--field", "DOMAINS"];
        if let Some(name) = flag_set {
            args.extend(["--flag-set", name]);
        }
        case.classify_against(vep_rs, &args)
    };
    let ended = |rows, fields| clean(0, &counted_with_fields(rows, fields));
    assert_eq!(classify("reordered.txt", None), ended(one, none));
    assert_eq!(classify("changed.txt", None), ended(one, one));
    assert_eq!(
        classify("reordered.txt", Some("inherits")),
        ended(one, none)
    );
    assert_eq!(classify("reordered.txt", Some("plain")), ended(one, one));
    assert_eq!(
        classify("reordered.txt", Some("regulatory")),
        ended(none, one)
    );
    let one_class = r#"{"unexplained_residual": 1}"#;
    let regulatory_row = written(
        &[divergence(
            &key("21:100", "G", ("ENSR1", "RegulatoryFeature"), &[0]),
            "regulatory_region_variant",
            &[],
            UNEXPLAINED,
        )],
        &[],
        (one_class, "{}"),
    );
    let domains = |vep_rs: &str| {
        let entry = field_divergence(
            &key("21:100", "G", transcript("ENST1"), &[0]),
            "missense_variant",
            ("DOMAINS", listed, vep_rs),
            UNEXPLAINED,
        );
        fields_written(&["DOMAINS"], &[entry], one_class)
    };
    let agreeing = fields_written(&["DOMAINS"], &[], "{}");
    let sets = [
        with_member(plain, &format!("{regulatory_row}, {}", domains(reordered))),
        with_member(
            masked,
            &format!("{}, {}", nothing_written(), domains(reordered)),
        ),
        with_member(inherits, &format!("{regulatory_row}, {agreeing}")),
    ];
    assert_eq!(
        case.after(),
        obj(&[
            masks,
            &held,
            &format!(r#""flag_sets": [{}]"#, sets.join(", ")),
            &regulatory_row,
            &domains(changed)
        ])
    );
    case.manifest(&obj(&[&held]));
    assert_eq!(classify("reordered.txt", None), ended(one, one));
    assert_eq!(
        case.after(),
        obj(&[&held, &regulatory_row, &domains(reordered)])
    );
}

#[test]
fn classify_records_field_divergences_on_agreeing_keys() {
    // A named field that differs on a key whose consequence sets agree is a field
    // divergence (not documented here, so `unexplained_residual`); a key whose sets differ
    // is a consequence divergence and its fields are not compared.
    let case = Case::new("agreeing");
    case.manifest(r#"{"records": []}"#);
    let swap = "start_lost,start_retained_variant";
    let missense = row_with(
        "rs1",
        "21:100",
        "G",
        "ENST1",
        "missense_variant",
        "IMPACT=MODERATE;STRAND=1;HGVSc=ENST1.1:c.5C>A;HGVSp=ENSP1.1:p.Ala2Asp",
    );
    let insertion = |protein: &str| {
        row_with(
            "rs2",
            "21:200",
            "-",
            "ENST2",
            "inframe_insertion",
            &format!("IMPACT=MODERATE;STRAND=1;HGVSc=ENST2.1:c.9_10insAGA;HGVSp={protein}"),
        )
    };
    case.outputs(
        &(missense.clone()
            + &insertion("ENSP2.1:p.Gly3_Lys4insArg")
            + &row_with(
                "rs3",
                "21:300",
                "-",
                "ENST3",
                swap,
                "IMPACT=HIGH;STRAND=1;HGVSp=ENSP3.1:p.Met1?",
            )),
        &(missense
            + &insertion("-")
            + &row_with(
                "rs3",
                "21:300",
                "-",
                "ENST3",
                "start_retained_variant",
                "IMPACT=LOW;STRAND=1;HGVSp=ENSP3.1:p.Met1del",
            )),
    );
    let ran = case.classify(&[
        "--field",
        "HGVSc",
        "--field",
        "HGVSp",
        "--flag=--hgvs",
        "--fasta-name",
        "reference.fa.gz",
    ]);
    assert_eq!(
        ran,
        clean(
            0,
            &counted_with_fields(
                (1, "{'start_cooccurrence_swap': 1}"),
                (1, "{'unexplained_residual': 1}")
            )
        )
    );
    let keys = written(
        &[divergence(
            &key("21:300", "-", transcript("ENST3"), &[]),
            swap,
            &["start_retained_variant"],
            "start_cooccurrence_swap",
        )],
        &[],
        (r#"{"start_cooccurrence_swap": 1}"#, "{}"),
    );
    let fields = fields_written(
        &["HGVSc", "HGVSp"],
        &[field_divergence(
            &key("21:200", "-", transcript("ENST2"), &[]),
            "inframe_insertion",
            ("HGVSp", "ENSP2.1:p.Gly3_Lys4insArg", "-"),
            UNEXPLAINED,
        )],
        r#"{"unexplained_residual": 1}"#,
    );
    // The run options the harness replays on both engines travel with the manifest.
    assert_eq!(
        case.after(),
        obj(&[
            r#""records": []"#,
            &keys,
            &fields,
            r#""flags": ["--hgvs"], "fasta": "reference.fa.gz""#
        ])
    );
}
