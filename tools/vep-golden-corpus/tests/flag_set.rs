// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `flag-set` on a copy of the two-record corpus: an entry is written with the keys
//! its options set, a set that exists is rewritten by name and keeps its other keys, the
//! sets are listed one tab-separated line each, and a run the harness could not make, or a
//! rewrite that would discard a hand-written statement, is refused before anything is
//! written.

mod common;

use common::{
    clean, failed, fixture, manifest, manifest_file, manifest_text, matrix, names,
    reference_statement, run, scratch, two_record_corpus, Run,
};
use std::os::unix::fs::PermissionsExt;
use vep_golden_corpus::constants::{mask_classes, WITH_REFERENCE_SETS};
use vep_golden_corpus::json::{dumps, dumps_object, loads, Object, Style, Value};
use vep_golden_corpus::manifest::manifest_shape_error;

/// The manifest's flag sets as `json.dumps` writes the list.
fn sets_text(manifest: &Object) -> String {
    dumps(manifest.get("flag_sets").unwrap(), &Style::PLAIN)
}

/// The entry at `place` of the manifest's flag sets.
fn set_at(manifest: &mut Object, place: usize) -> &mut Object {
    let sets = manifest.get_mut("flag_sets").and_then(Value::as_list_mut);
    sets.unwrap()[place].as_object_mut().unwrap()
}

/// Writes a manifest as `json.dumps(manifest)` writes it, on one line.
fn write_on_one_line(corpus: &str, manifest: &Object) {
    let text = dumps_object(manifest, &Style::PLAIN);
    std::fs::write(manifest_file(corpus), text).unwrap();
}

#[test]
fn flag_set_subcommand_writes_an_entry() {
    let scratch = scratch("flag-set");
    let out = two_record_corpus(&scratch.path);
    let matrix = matrix(&scratch.path, &["symbol", "pick", "per_gene"]);
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "pick",
            "--flag=--pick",
            "--flag=--symbol",
            "--covers",
            "pick",
            "--matrix",
            &matrix
        ]),
        clean(
            0,
            "flag set pick: 2 flags, input variants.vcf, formats default,tab,vcf,json, expected expected_pick, covers 1\n"
        )
    );
    assert_eq!(
        sets_text(&manifest(&out)),
        concat!(
            r#"[{"name": "pick", "flags": ["--pick", "--symbol"], "input": "variants.vcf", "#,
            r#""stdin": false, "unordered": false, "formats": ["default", "tab", "vcf", "json"], "#,
            r#""expected_dir": "expected_pick", "covers": ["pick"]}]"#
        )
    );
    std::fs::write(format!("{out}/variants.ensembl"), "21\t100\t100\tA/G\t+\n").unwrap();
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "ensembl-stdin",
            "--flag=--format",
            "--flag=ensembl",
            "--input",
            "variants.ensembl",
            "--stdin",
            "--format",
            "default,json",
            "--expected-dir",
            "expected_ens",
            "--matrix",
            &matrix
        ])
        .code,
        0
    );
    // A run whose reference orders a record's rows by Ensembl VEP's hash order is marked
    // unordered.
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "pick",
            "--flag=--per_gene",
            "--unordered",
            "--covers",
            "per_gene",
            "--matrix",
            &matrix
        ]),
        clean(
            0,
            "flag set pick: 1 flags, input variants.vcf, unordered, formats default,tab,vcf,json, expected expected_pick, covers 1\n"
        )
    );
    assert_eq!(
        sets_text(&manifest(&out)),
        concat!(
            r#"[{"name": "pick", "flags": ["--per_gene"], "input": "variants.vcf", "#,
            r#""stdin": false, "unordered": true, "formats": ["default", "tab", "vcf", "json"], "#,
            r#""expected_dir": "expected_pick", "covers": ["per_gene"]}, "#,
            r#"{"name": "ensembl-stdin", "flags": ["--format", "ensembl"], "input": "variants.ensembl", "#,
            r#""stdin": true, "unordered": false, "formats": ["default", "json"], "#,
            r#""expected_dir": "expected_ens", "covers": []}]"#
        )
    );
    // `--format=<format>` is the other spelling of the format a stdin run must name.
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "eq",
            "--flag=--format=vcf",
            "--stdin",
            "--matrix",
            &matrix
        ])
        .code,
        0
    );
    // A reference generated against the Ensembl database is marked, and only then.
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "db",
            "--flag=--ga4gh_vrs",
            "--format",
            "default,tab,json",
            "--perl-mode",
            "database",
            "--matrix",
            &matrix
        ]),
        clean(
            0,
            "flag set db: 1 flags, input variants.vcf, reference database, formats default,tab,json, expected expected_db, covers 0\n"
        )
    );
    let mut written = manifest(&out);
    assert_eq!(
        set_at(&mut written, 3).get("perl_mode"),
        Some(&Value::from("database"))
    );
    assert!(!set_at(&mut written, 2).contains_key("perl_mode"));
    assert_eq!(manifest_shape_error(&written), None);
    set_at(&mut written, 3).insert("perl_mode", "cache");
    assert_eq!(
        manifest_shape_error(&written).as_deref(),
        Some("flag set db perl_mode is not one of offline, database, gtf, gff")
    );
    // A reference Ensembl VEP wrote reading its transcripts from a GTF or GFF is marked the
    // same way, and the manifest's own run may carry the mode, which `select` keeps on a
    // regeneration and refuses when it is not one it knows.
    assert_eq!(
        run(&["flag-set", &out, "--name", "gtf", "--flag=--symbol", "--perl-mode", "gtf", "--matrix", &matrix]),
        clean(
            0,
            "flag set gtf: 1 flags, input variants.vcf, reference gtf, formats default,tab,vcf,json, expected expected_gtf, covers 0\n"
        )
    );
    let mut written = manifest(&out);
    assert_eq!(
        set_at(&mut written, 4).get("perl_mode"),
        Some(&Value::from("gtf"))
    );
    written.insert("perl_mode", "gff");
    assert_eq!(manifest_shape_error(&written), None);
    write_on_one_line(&out, &written);
    let dataset = format!(
        "pool={}={}",
        fixture("two_record/ref.txt").to_str().unwrap(),
        fixture("two_record/in.vcf").to_str().unwrap()
    );
    let regenerated = run(&[
        "select",
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--dataset",
        &dataset,
        "--out",
        &out,
        "--k",
        "1",
        "--seed",
        "1",
        "--multi-allelic-per-dataset",
        "0",
        "--covers",
        "symbol",
    ]);
    assert_eq!(
        regenerated,
        clean(
            0,
            "combinations observed 2, covered 2, records 2, unresolved exemplars 0, variant classes 1 (snv 2)\n"
        )
    );
    let mut written = manifest(&out);
    assert_eq!(written.get("perl_mode"), Some(&Value::from("gff")));
    written.insert("perl_mode", "cache");
    assert_eq!(
        manifest_shape_error(&written).as_deref(),
        Some("perl_mode is not one of offline, database, gtf, gff")
    );
    write_on_one_line(&out, &written);
    assert_eq!(
        run(&["flag-set", &out, "--list"]),
        failed(&format!(
            "{out}/manifest.json perl_mode is not one of offline, database, gtf, gff"
        ))
    );
    written.remove("perl_mode");
    write_on_one_line(&out, &written);
    assert_eq!(
        run(&["flag-set", &out, "--list"]),
        clean(
            0,
            concat!(
                "pick\t--per_gene\tvariants.vcf\tfalse\ttrue\tdefault,tab,vcf,json\texpected_pick\tper_gene\n",
                "ensembl-stdin\t--format ensembl\tvariants.ensembl\ttrue\tfalse\tdefault,json\texpected_ens\t\n",
                "eq\t--format=vcf\tvariants.vcf\ttrue\tfalse\tdefault,tab,vcf,json\texpected_eq\t\n",
                "db\t--ga4gh_vrs\tvariants.vcf\tfalse\tfalse\tdefault,tab,json\texpected_db\t\n",
                "gtf\t--symbol\tvariants.vcf\tfalse\tfalse\tdefault,tab,vcf,json\texpected_gtf\t\n",
            )
        )
    );
    let before = manifest_text(&out);
    let not_a_file = format!("input missing.vcf is not a file in {out}");
    let no_row = format!("flag set x covers names ghost, which is not a row of {matrix}");
    let refusals: [(&[&str], &str); 8] = [
        (
            &["--name", "default"],
            "flag set name 'default': letters, digits, _ and -, and not the manifest's own run, default",
        ),
        (
            &["--name", "a/b"],
            "flag set name 'a/b': letters, digits, _ and -, and not the manifest's own run, default",
        ),
        (
            &["--name", "x", "--format", "default,parquet"],
            "--format wants a non-empty subset of default,tab,vcf,json, got 'default,parquet'",
        ),
        // A set comparing no format would certify its covers without a comparison.
        (
            &["--name", "x", "--format", ","],
            "--format wants a non-empty subset of default,tab,vcf,json, got ','",
        ),
        (&["--name", "x", "--input", "missing.vcf"], &not_a_file),
        (&["--name", "x", "--covers", "ghost"], &no_row),
        // Ensembl VEP detects no format on standard input, so a stdin run without --format
        // has no reference.
        (
            &["--name", "x", "--stdin"],
            "a stdin run names its input format: add --flag=--format --flag=<format>",
        ),
        (&[], "flag-set wants --name (or --list)"),
    ];
    for (extra, message) in refusals {
        let mut args = vec!["flag-set", out.as_str()];
        args.extend(extra);
        args.extend(["--matrix", matrix.as_str()]);
        assert_eq!(run(&args), failed(message), "{extra:?}");
    }
    // A path that leaves the corpus directory, or an empty covers name, is refused as well.
    let outside = scratch.path.join("outside.vcf");
    std::fs::write(&outside, "").unwrap();
    let outside = outside.to_str().unwrap();
    let absolute = format!("--input wants a path inside the corpus directory, got '{outside}'");
    let refusals: [(&[&str], &str); 4] = [
        (&["--input", outside], &absolute),
        (
            &["--input", "../outside.vcf"],
            "--input wants a path inside the corpus directory, got '../outside.vcf'",
        ),
        (
            &["--expected-dir", "../expected_x"],
            "--expected-dir wants a path inside the corpus directory, got '../expected_x'",
        ),
        (&["--covers", ""], "--covers wants a non-empty row name"),
    ];
    for (extra, message) in refusals {
        let mut args = vec!["flag-set", out.as_str(), "--name", "x"];
        args.extend(extra);
        args.extend(["--matrix", matrix.as_str()]);
        assert_eq!(run(&args), failed(message), "{extra:?}");
    }
    assert_eq!(manifest_text(&out), before);
}

#[test]
fn flag_set_rewrite_keeps_the_sets_hand_written_keys() {
    let scratch = scratch("flag-set-rewrite");
    let out = two_record_corpus(&scratch.path);
    let matrix = matrix(&scratch.path, &["symbol", "domains", "per_gene"]);
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "domains",
            "--flag=--domains",
            "--covers",
            "domains",
            "--perl-mode",
            "database",
            "--matrix",
            &matrix
        ])
        .code,
        0
    );
    let notice = "found tandem repeats with different references per alternative allele: SVLEN=9,6; only using reference with largest size";
    // The hand-written keys sit before `covers`, as a hand-edited manifest places them.
    let mut written = manifest(&out);
    let entry = set_at(&mut written, 0);
    let covers = entry.remove("covers").unwrap();
    entry.insert("masks", Value::strings(&["order_only"]));
    entry.insert("reference_warnings", Value::strings(&[notice]));
    entry.insert("covers", covers);
    entry.insert(
        "divergences",
        loads(r#"[{"location": "21:100", "record_indices": [0]}]"#).unwrap(),
    );
    entry.insert("divergence_summary", Object::new());
    entry.insert("note", "kept as written");
    write_on_one_line(&out, &written);
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "domains",
            "--flag=--domains",
            "--flag=--per_gene",
            "--unordered",
            "--covers",
            "domains",
            "--covers",
            "per_gene",
            "--matrix",
            &matrix
        ]),
        clean(
            0,
            "flag set domains: 2 flags, input variants.vcf, unordered, formats default,tab,vcf,json, expected expected_domains, covers 2\n"
        )
    );
    let rewritten = manifest(&out);
    assert_eq!(
        sets_text(&rewritten),
        format!(
            concat!(
                r#"[{{"name": "domains", "flags": ["--domains", "--per_gene"], "input": "variants.vcf", "#,
                r#""stdin": false, "unordered": true, "formats": ["default", "tab", "vcf", "json"], "#,
                r#""expected_dir": "expected_domains", "masks": ["order_only"], "#,
                r#""reference_warnings": ["{}"], "covers": ["domains", "per_gene"], "#,
                r#""note": "kept as written"}}]"#
            ),
            notice
        )
    );
    assert_eq!(manifest_shape_error(&rewritten), None);
    // A kept key is still held to its shape at the next read, so a rewrite never writes
    // past a `masks` or `reference_warnings` the harness would refuse.
    let unknown = format!(
        "flag set domains masks names unknown classes ['splice_rows'] (the classes are {})",
        mask_classes().join(", ")
    );
    let broken_keys = [
        ("masks", Value::strings(&["splice_rows"]), unknown.as_str()),
        (
            "masks",
            Value::from("order_only"),
            "flag set domains masks is not a list of names",
        ),
        (
            "reference_warnings",
            Value::strings(&[notice, ""]),
            "flag set domains reference_warnings is not a list of texts",
        ),
        (
            "reference_warnings",
            Value::from(notice),
            "flag set domains reference_warnings is not a list of texts",
        ),
    ];
    for (key, value, reason) in broken_keys {
        let mut broken = manifest(&out);
        set_at(&mut broken, 0).insert(key, value);
        write_on_one_line(&out, &broken);
        assert_eq!(
            run(&[
                "flag-set",
                &out,
                "--name",
                "domains",
                "--flag=--domains",
                "--matrix",
                &matrix
            ]),
            failed(&format!("{out}/manifest.json {reason}"))
        );
        assert_eq!(manifest(&out), broken, "nothing is written");
    }
}

#[test]
fn flag_set_takes_every_format_for_an_empty_list_and_keeps_a_format_once_a_name() {
    let scratch = scratch("flag-set-formats");
    let out = two_record_corpus(&scratch.path);
    assert_eq!(
        run(&["flag-set", &out, "--name", "all", "--format", ""]),
        clean(
            0,
            "flag set all: 0 flags, input variants.vcf, formats default,tab,vcf,json, expected expected_all, covers 0\n"
        )
    );
    assert_eq!(
        run(&["flag-set", &out, "--name", "some", "--format", " json , tab ,json"]),
        clean(
            0,
            "flag set some: 0 flags, input variants.vcf, formats json,tab,json, expected expected_some, covers 0\n"
        )
    );
    let written = manifest(&out);
    let sets = written.get("flag_sets").and_then(Value::as_list).unwrap();
    assert_eq!(
        names(sets[1].as_object().unwrap(), "formats"),
        ["json", "tab", "json"]
    );
    // The entry goes last in a manifest that held no flag sets, after every other key.
    assert_eq!(written.keys().last(), Some("flag_sets"));
}

#[test]
fn flag_set_rewrite_refuses_to_drop_a_with_reference_key() {
    // A rewrite of a set drops that set's divergence keys, so a `divergences` entry of the
    // set that carries the with-reference sets stops it, named, with nothing written;
    // `--drop-reference-statements` names the entry on standard error and rewrites. The
    // manifest's own list and another set's are not the rewrite's: they stay and stop
    // nothing.
    let scratch = scratch("flag-set-statements");
    let out = two_record_corpus(&scratch.path);
    let matrix = matrix(&scratch.path, &["symbol", "pick", "domains"]);
    for name in ["pick", "domains"] {
        let flag = format!("--flag=--{name}");
        let written = run(&[
            "flag-set", &out, "--name", name, &flag, "--covers", name, "--matrix", &matrix,
        ]);
        assert_eq!(written.code, 0, "{name}");
    }
    let mut stated = manifest(&out);
    stated.insert(
        "divergences",
        vec![reference_statement("21:100", "G", "ENST1", 0)],
    );
    set_at(&mut stated, 0).insert(
        "divergences",
        vec![reference_statement("21:200", "T", "ENST2", 1)],
    );
    set_at(&mut stated, 1).insert(
        "divergences",
        vec![reference_statement("21:100", "G", "ENST1", 0)],
    );
    write_on_one_line(&out, &stated);
    let held = manifest_text(&out);
    let rewrite = [
        "flag-set",
        &out,
        "--name",
        "pick",
        "--flag=--pick",
        "--flag=--symbol",
        "--covers",
        "pick",
        "--matrix",
        &matrix,
    ];
    assert_eq!(
        run(&rewrite),
        failed(&format!(
            "a flag-set rewrite would discard {WITH_REFERENCE_SETS}, a key written by hand, on flag set pick divergences entry on 21:200 T ENST2: pass --drop-reference-statements to discard it"
        ))
    );
    assert_eq!(manifest_text(&out), held);
    let rewritten = "flag set pick: 2 flags, input variants.vcf, formats default,tab,vcf,json, expected expected_pick, covers 1\n";
    let mut allowed = rewrite.to_vec();
    allowed.push("--drop-reference-statements");
    // A manifest that cannot be written is found before an entry is named: the run ends on
    // the one line of the write. A user who may write any file skips the check.
    let file = manifest_file(&out);
    let mode = |bits: u32| {
        std::fs::set_permissions(&file, std::fs::Permissions::from_mode(bits)).unwrap();
    };
    mode(0o444);
    if std::fs::OpenOptions::new().write(true).open(&file).is_err() {
        let unwritten = run(&allowed);
        let line = format!(
            "ERROR: [vep-golden-corpus] {} is not writable (",
            file.display()
        );
        assert_eq!((unwritten.code, unwritten.stdout.as_str()), (1, ""));
        assert!(unwritten.stderr.starts_with(&line), "{}", unwritten.stderr);
        assert_eq!(unwritten.stderr.lines().count(), 1, "{}", unwritten.stderr);
        assert_eq!(manifest_text(&out), held);
    }
    mode(0o644);
    assert_eq!(
        run(&allowed),
        Run {
            code: 0,
            stdout: rewritten.to_string(),
            stderr: format!(
                "WARNING: [build_golden_corpus] a flag-set rewrite discards {WITH_REFERENCE_SETS} on flag set pick divergences entry on 21:200 T ENST2\n"
            ),
        }
    );
    let mut after = manifest(&out);
    assert_eq!(
        names(set_at(&mut after, 0), "flags"),
        ["--pick", "--symbol"]
    );
    assert!(!set_at(&mut after, 0).contains_key("divergences"));
    assert_eq!(after.get("divergences"), stated.get("divergences"));
    assert_eq!(set_at(&mut after, 1), set_at(&mut stated, 1));
    // The set holds no statement after the rewrite, and a set written for the first time
    // holds none.
    assert_eq!(run(&rewrite), clean(0, rewritten));
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "both",
            "--flag=--pick",
            "--flag=--domains",
            "--matrix",
            &matrix,
        ]),
        clean(
            0,
            "flag set both: 2 flags, input variants.vcf, formats default,tab,vcf,json, expected expected_both, covers 0\n"
        )
    );
}
