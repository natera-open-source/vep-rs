// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up runs whose records have no row: a record Ensembl VEP warned
//! about is marked skipped and its vep-rs keys take the class of a skipped record, a
//! record that has no row by its nature is kept, any other is refused; the records of a
//! flag set over an input of its own are numbered as that input's lines, and a region
//! line among them names its record.

mod classify_support;
mod common;

use classify_support::{
    counted, divergence, gzip, js, key, nothing_written, obj, record, record_with, records, row,
    row_with, transcript, vep_rs_only, written, Case,
};
use common::{clean, failed, manifest, object, refused, Run};
use vep_golden_corpus::json::{dumps, dumps_object, Style, Value};

const RS1: &str = "21\t100\trs1\tA\tG\t.\tPASS\t.";
const NO_ROWS: (usize, &str) = (0, "{}");

/// The row both outputs hold for the first record.
fn rs1_row() -> String {
    row("rs1", "21:100", "G", "ENST1", "missense_variant")
}

/// The refusal of a record without a row or a warning.
fn rowless(ordinal: usize, label: &str) -> String {
    format!("record {ordinal} ({label}) has no VEP row and no reference_warning")
}

/// A run that ended with exit code 1 on the error lines given.
fn problems(lines: &[String]) -> Run {
    let lines: Vec<String> = lines
        .iter()
        .map(|line| format!("ERROR: [build_golden_corpus] {line}\n"))
        .collect();
    Run {
        code: 1,
        stdout: String::new(),
        stderr: lines.concat(),
    }
}

/// The manifest's `records` as `json.dumps` writes the list.
fn records_text(case: &Case) -> String {
    dumps(
        manifest(&case.corpus).get("records").unwrap(),
        &Style::PLAIN,
    )
}

/// The flag set `name` of the manifest as `json.dumps` writes the entry.
fn set_text(case: &Case, name: &str) -> String {
    let manifest = manifest(&case.corpus);
    let sets = manifest.get("flag_sets").and_then(Value::as_list).unwrap();
    let set = sets
        .iter()
        .find(|set| set.get("name").and_then(Value::as_str) == Some(name));
    dumps(set.unwrap(), &Style::PLAIN)
}

/// A flag set entry as `flag-set` writes one without flags of its own to name.
fn set(name: &str, flags: &str, input: Option<&str>) -> String {
    let input = input.map_or(String::new(), |input| format!(r#""input": "{input}", "#));
    format!(
        r#"{{"name": "{name}", "flags": {flags}, {input}"expected_dir": "expected_{name}", "covers": []}}"#
    )
}

#[test]
fn classify_marks_reference_skipped_records_and_their_tuples() {
    // A record no Ensembl VEP row names carries Ensembl VEP's reason as
    // `reference_warning`; classify writes `reference_rows` 0 on it, classes every key only
    // vep-rs writes for it `reference_skipped_record`, and removes a `reference_rows`
    // another run left on a record that has rows.
    let case = Case::new("skipped");
    let warned = "21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900";
    let warning = r#""reference_warning": "CPX is not a supported structural variant type""#;
    case.manifest(&obj(&[&records(&[
        record_with(RS1, r#""reference_rows": 1"#),
        record_with(warned, warning),
    ])]));
    case.outputs(
        &rs1_row(),
        &(rs1_row()
            + &row("sv1", "21:501-900", "CPX", "ENST2", "feature_truncation")
            + &row(
                "sv1",
                "21:501-900",
                "CPX",
                "ENST3",
                "feature_truncation,intron_variant",
            )),
    );
    let tuples = (2, "{'reference_skipped_record': 2}");
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted(NO_ROWS, tuples, "[1]"))
    );
    // The two markers are the record's last two keys, in this order.
    let held = records(&[
        record(RS1),
        record_with(warned, &format!(r#""reference_rows": 0, {warning}"#)),
    ]);
    let on = |feature: &str, set: &str| {
        vep_rs_only(
            &key("21:501-900", "CPX", transcript(feature), &[1]),
            &[set],
            "reference_skipped_record",
        )
    };
    let expected = written(
        &[],
        &[
            on("ENST2", "feature_truncation"),
            on("ENST3", "feature_truncation,intron_variant"),
        ],
        ("{}", r#"{"reference_skipped_record": 2}"#),
    );
    assert_eq!(case.after(), obj(&[&held, &expected]));
}

#[test]
fn classify_refuses_a_zero_row_record_without_a_warning() {
    // A record no Ensembl VEP row names and no `reference_warning` explains is an error:
    // its name may be unresolved, or the reason was not recorded.
    let case = Case::new("zero-row");
    case.manifest(&obj(&[&records(&[
        record(RS1),
        record("21\t150\trs9\tC\tT\t.\tPASS\t."),
    ])]));
    let before = case.text();
    case.outputs(&rs1_row(), &rs1_row());
    assert_eq!(case.classify(&[]), failed(&rowless(1, "21 150 rs9 C T")));
    assert_eq!(case.text(), before, "a refused run writes nothing");
}

/// Runs the binary with `args` on a standard input that is a pipe holding `rows`, its
/// writing end closed: a read of `/dev/stdin` takes the rows and ends.
#[cfg(unix)]
fn run_on_a_pipe(args: &[&str], rows: &str) -> Run {
    use std::io::Write;
    let (reader, mut writer) = std::io::pipe().expect("a pipe");
    writer.write_all(rows.as_bytes()).expect("room in the pipe");
    drop(writer);
    let output = std::process::Command::new(common::BIN)
        .current_dir(common::ROOT)
        .env_remove("COLUMNS")
        .args(args)
        .stdin(reader)
        .output()
        .expect("the binary runs");
    Run {
        code: output.status.code().expect("an exit code"),
        stdout: String::from_utf8(output.stdout).expect("UTF-8 on standard output"),
        stderr: String::from_utf8(output.stderr).expect("UTF-8 on standard error"),
    }
}

/// The arguments of a `classify` run over `corpus` on the two outputs given.
#[cfg(unix)]
fn classify_arguments<'a>(corpus: &'a str, vep: &'a str, vep_rs: &'a str) -> [&'a str; 7] {
    [
        "classify",
        "--corpus",
        corpus,
        "--vep-default",
        vep,
        "--vep-rs-output",
        vep_rs,
    ]
}

#[cfg(unix)]
#[test]
fn classify_refuses_an_output_that_can_be_read_once() {
    // classify reads each output more than once, and a pipe gives its rows to the first
    // read alone: every later read would find no row and the counts would be those of an
    // empty output. Either output named as a pipe is refused before anything is read.
    let case = Case::new("read-once");
    case.manifest(&obj(&[&records(&[record(RS1)])]));
    let before = case.text();
    case.outputs(&rs1_row(), &rs1_row());
    // The refusal is the tool's own, so it carries the tool's name.
    let refusal = |output: &str| {
        refused(
            1,
            &format!(
                "ERROR: [vep-golden-corpus] {output} can be read once and not again, and \
                 classify reads each output more than once"
            ),
        )
    };
    let (vep, vep_rs) = (case.path("vep.txt"), case.path("rs.txt"));
    let classify = |vep, vep_rs| classify_arguments(&case.corpus, vep, vep_rs);
    // The pipe of standard input comes before the named one: a run that reads it ends, on
    // the counts of its rows, where a run that opens a named pipe waits for a writer.
    for arguments in [
        classify("/dev/stdin", &vep_rs),
        classify(&vep, "/dev/stdin"),
    ] {
        assert_eq!(run_on_a_pipe(&arguments, &rs1_row()), refusal("/dev/stdin"));
        assert_eq!(case.text(), before, "a refused run writes nothing");
    }
    let pipe = case.path("pipe");
    let made = std::process::Command::new("mkfifo")
        .arg(&pipe)
        .status()
        .expect("mkfifo runs");
    assert!(made.success());
    assert_eq!(case.classify_against("pipe", &[]), refusal(&pipe));
    assert_eq!(common::run(&classify(&pipe, &vep_rs)), refusal(&pipe));
    assert_eq!(case.text(), before, "a refused run writes nothing");
}

#[test]
fn classify_keeps_an_unwarned_non_variant_record() {
    // A record whose ALT is `.` carries no row by its nature, so one without a warning is
    // neither refused nor marked skipped; one Ensembl VEP warned about is skipped as any
    // other.
    let case = Case::new("non-variant");
    let non_variant = record("21\t200\tnonvar\tTG\t.\t.\tPASS\t.");
    let warned = "21\t300\tnonvar_bad\tG\t.\t.\tPASS\t.";
    let warning = r#""reference_warning": "Specified reference allele G does not match Ensembl reference allele""#;
    case.manifest(&obj(&[&records(&[
        record(RS1),
        non_variant.clone(),
        record_with(warned, warning),
    ])]));
    case.outputs(&rs1_row(), &rs1_row());
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted(NO_ROWS, NO_ROWS, "[2]"))
    );
    let held = records(&[
        record(RS1),
        non_variant,
        record_with(warned, &format!(r#""reference_rows": 0, {warning}"#)),
    ]);
    assert_eq!(case.after(), obj(&[&held, &nothing_written()]));
}

#[test]
fn classify_refuses_a_warned_record_that_has_rows() {
    let case = Case::new("warned");
    case.manifest(&obj(&[&records(&[record_with(
        RS1,
        r#""reference_warning": "deletion looks incomplete""#,
    )])]));
    let before = case.text();
    case.outputs(&rs1_row(), &rs1_row());
    assert_eq!(
        case.classify(&[]),
        failed(
            "record 0 (21 100 rs1 A G) has VEP rows and reference_warning 'deletion looks incomplete'"
        )
    );
    assert_eq!(case.text(), before);
}

/// One object of a reference's JSON output, as `json.dumps` writes it, on a line.
fn json_line(id: &str, input: &str, most_severe: &str) -> String {
    format!(
        r#"{{"id": {}, "input": {}, "most_severe_consequence": {}}}"#,
        js(id),
        js(input),
        js(most_severe)
    ) + "\n"
}

#[test]
fn classify_keeps_a_record_the_reference_carried_without_a_row() {
    // A record whose JSON object the reference wrote with `most_severe_consequence` `?` has
    // no allele to annotate, so one without a warning is kept, not refused or marked
    // skipped; the object is matched by the record's first five columns, so a twin input
    // carrying sample columns names it too; a record the JSON lacks is refused, and a run
    // whose expected directory has no JSON file falls back to the ALT `.` rule.
    let case = Case::new("carried");
    let held = records(&[
        record(RS1),
        record("21\t200\talt_eq_ref\tC\tC\t.\tPASS\t."),
        record("21\t300\tlone_delstar\tC\t<DEL:*>\t.\tPASS\t."),
        record("21\t400\tstar_only\tA\t*\t.\tPASS\t."),
    ]);
    case.manifest(&obj(&[&held]));
    let before = case.text();
    let objects = json_line("rs1", RS1, "missense_variant")
        + &json_line(
            "alt_eq_ref",
            "21\t200\talt_eq_ref\tC\tC\t.\tPASS\t.\tGT\t0/1",
            "?",
        )
        + &json_line(
            "lone_delstar",
            "21\t300\tlone_delstar\tC\t<DEL:*>\t.\tPASS\t.",
            "?",
        );
    case.write_in_corpus("expected/json.jsonl.gz", &gzip(&objects));
    case.outputs(&rs1_row(), &rs1_row());
    let star = rowless(3, "21 400 star_only A *");
    assert_eq!(case.classify(&[]), failed(&star));
    assert_eq!(case.text(), before);
    let all = objects + &json_line("star_only", "21\t400\tstar_only\tA\t*\t.\tPASS\t.", "?");
    case.write_in_corpus("expected/json.jsonl.gz", &gzip(&all));
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted(NO_ROWS, NO_ROWS, "[]"))
    );
    assert_eq!(case.after(), obj(&[&held, &nothing_written()]));
    std::fs::remove_file(format!("{}/expected/json.jsonl.gz", case.corpus)).unwrap();
    assert_eq!(
        case.classify(&[]),
        problems(&[
            rowless(1, "21 200 alt_eq_ref C C"),
            rowless(2, "21 300 lone_delstar C <DEL:*>"),
            star
        ])
    );
}

#[test]
fn classify_reads_a_flag_sets_own_json_and_refuses_an_unreadable_one() {
    // A flag set's carried records come from its own expected JSON when it compares that
    // format, else from the default run's; a JSON file that is not gzip, or holds a line
    // that is not JSON or not a JSON object, is refused with its path.
    let case = Case::new("own-json");
    let held = records(&[record(RS1), record("21\t200\talt_eq_ref\tC\tC\t.\tPASS\t.")]);
    let pick = r#"["--pick"]"#;
    let shared =
        r#"{"name": "shared", "flags": ["--pick"], "expected_dir": "expected", "covers": []}"#;
    let sets = format!(r#""flag_sets": [{}, {shared}]"#, set("own", pick, None));
    case.manifest(&obj(&[&held, &sets]));
    std::fs::create_dir(format!("{}/expected", case.corpus)).unwrap();
    let carried = json_line("alt_eq_ref", "21\t200\talt_eq_ref\tC\tC\t.\tPASS\t.", "?");
    case.write_in_corpus("expected_own/json.jsonl.gz", &gzip(&carried));
    case.outputs(&rs1_row(), &rs1_row());
    assert_eq!(
        case.classify(&["--flag-set", "own"]),
        clean(0, &counted(NO_ROWS, NO_ROWS, "[]"))
    );
    let own = format!(
        r#"{{"name": "own", "flags": {pick}, "expected_dir": "expected_own", "covers": [], {}}}"#,
        nothing_written()
    );
    let after = obj(&[&held, &format!(r#""flag_sets": [{own}, {shared}]"#)]);
    assert_eq!(case.after(), after);
    let written_once = case.text();
    assert_eq!(
        case.classify(&["--flag-set", "shared"]),
        failed(&rowless(1, "21 200 alt_eq_ref C C"))
    );
    let json = format!("{}/expected/json.jsonl.gz", case.corpus);
    let not_json =
        |reason: &str| failed(&format!("{json} is not a readable JSON output ({reason})"));
    case.write_in_corpus("expected/json.jsonl.gz", b"not gzip");
    assert_eq!(
        case.classify(&["--flag-set", "shared"]),
        not_json("Not a gzipped file (b'no')")
    );
    case.write_in_corpus("expected/json.jsonl.gz", &gzip("{not json\n"));
    assert_eq!(
        case.classify(&[]),
        not_json("Expecting property name enclosed in double quotes: line 1 column 2 (char 1)")
    );
    case.write_in_corpus("expected/json.jsonl.gz", &gzip(&(carried + "[]\n")));
    assert_eq!(case.classify(&[]), not_json("line 2 is not a JSON object"));
    assert_eq!(case.text(), written_once);
}

#[test]
fn classify_keeps_a_record_a_row_filter_left_rowless_on_both_sides() {
    // Under `--coding_only` or `--no_intergenic` a record with no coding (or no intergenic)
    // row has no row on either side, so a flag set or a default run under one keeps such a
    // record without a warning, not refused and not marked skipped, with no JSON to carry
    // it; the same record under no row filter, or with vep-rs rows and no reference row, is
    // refused.
    let case = Case::new("row-filter");
    let held = records(&[record(RS1), record("21\t150\trs9\tC\tT\t.\tPASS\t.")]);
    let coding = set("coding-only", r#"["--coding_only"]"#, None);
    let intergenic = set("no-intergenic", r#"["--no_intergenic"]"#, None);
    let pick = set("pick", r#"["--pick"]"#, None);
    case.manifest(&obj(&[
        &held,
        &format!(r#""flag_sets": [{coding}, {intergenic}, {pick}]"#),
    ]));
    case.outputs(&rs1_row(), &rs1_row());
    let nothing = clean(0, &counted(NO_ROWS, NO_ROWS, "[]"));
    for name in ["coding-only", "no-intergenic"] {
        assert_eq!(case.classify(&["--flag-set", name]), nothing, "{name}");
    }
    assert_eq!(case.classify(&["--flag=--coding_only"]), nothing);
    let done = |entry: &str| with_written(entry);
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &format!(
                r#""flag_sets": [{}, {}, {pick}]"#,
                done(&coding),
                done(&intergenic)
            ),
            &nothing_written(),
            r#""flags": ["--coding_only"]"#
        ])
    );
    let after = case.text();
    let refused = failed(&rowless(1, "21 150 rs9 C T"));
    assert_eq!(case.classify(&["--flag-set", "pick"]), refused);
    case.write(
        "rs.txt",
        &(rs1_row() + &row("rs9", "21:150", "T", "ENST2", "missense_variant")),
    );
    assert_eq!(case.classify(&["--flag-set", "coding-only"]), refused);
    assert_eq!(case.text(), after);
}

#[test]
fn classify_names_a_region_record_by_its_line() {
    // A region line names its record in both outputs (Parser/Region.pm:198 and :213 give
    // the record its line as its name, which OutputFactory.pm:895-901 writes in the first
    // column), so the rows of a flag set over region lines resolve by those lines, with or
    // without `--format region` among its flags, and never by the IDs of `variants.vcf`.
    let case = Case::new("region");
    let held = list_of(&[
        record(RS1),
        record("21\t149\tdel1\tGC\tG\t.\tPASS\t."),
        record("21\t300\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=900"),
    ]);
    let sets = [
        set(
            "region",
            r#"["--format", "region"]"#,
            Some("variants.region.txt"),
        ),
        set("guess", "[]", Some("variants.region.txt")),
    ];
    case.manifest(&format!(
        r#"{{"records": {held}, "flag_sets": [{}]}}"#,
        sets.join(", ")
    ));
    let (snv, deletion, structural) = ("21:100-100:1/G", "21:150-150/-", "21:301-900/DEL");
    case.write_in_corpus(
        "variants.region.txt",
        format!("{snv}\n{deletion}\n{structural}\n").as_bytes(),
    );
    let snv_row = row(snv, "21:100", "G", "ENST1", "missense_variant");
    let deletion_row = row(deletion, "21:150", "-", "ENST1", "frameshift_variant");
    let truncation = |feature: &str| {
        row(
            structural,
            "21:301-900",
            "deletion",
            feature,
            "feature_truncation",
        )
    };
    let both = snv_row.clone() + &deletion_row + &truncation("ENST2");
    case.outputs(&both, &(both.clone() + &truncation("ENST3")));
    let one = (1, "{'unexplained_residual': 1}");
    let entries = written(
        &[],
        &[vep_rs_only(
            &key("21:301-900", "deletion", transcript("ENST3"), &[2]),
            &["feature_truncation"],
            "unexplained_residual",
        )],
        ("{}", r#"{"unexplained_residual": 1}"#),
    );
    for (name, place) in [("region", 0), ("guess", 1)] {
        assert_eq!(
            case.classify(&["--flag-set", name]),
            clean(0, &counted(NO_ROWS, one, "[]")),
            "{name}"
        );
        let after = format!("{}, {entries}}}", sets[place].strip_suffix('}').unwrap());
        assert_eq!(set_text(&case, name), after, "{name}");
    }
    assert_eq!(records_text(&case), held);

    // A region record without a row is refused by its ordinal and its line's position, the
    // strand where the line writes one, and never by a cell that holds an allele.
    let after = case.text();
    case.outputs(&truncation("ENST2"), &truncation("ENST2"));
    assert_eq!(
        case.classify(&["--flag-set", "region"]),
        problems(&[rowless(0, "21:100-100:1"), rowless(1, "21:150-150")])
    );
    // Rows named by the IDs of `variants.vcf` are no rows of a region record.
    let by_id = row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + &row("del1", "21:150", "-", "ENST1", "frameshift_variant")
        + &truncation("ENST2");
    case.outputs(&by_id, &by_id);
    assert_eq!(
        case.classify(&["--flag-set", "region"]),
        problems(&[rowless(0, "21:100-100:1"), rowless(1, "21:150-150")])
    );
    assert_eq!(case.text(), after, "a refused run writes nothing");
}

/// A flag set entry with the members of a run without a difference after its own.
fn with_written(entry: &str) -> String {
    format!(
        "{}, {}}}",
        entry.strip_suffix('}').unwrap(),
        nothing_written()
    )
}

#[test]
fn classify_numbers_a_flag_sets_own_input_records() {
    // A flag set over another `input` runs that file's records, numbered as the golden
    // tests number them (every line that is no header, in order), so its divergences and
    // skipped records resolve by those ordinals; a line the manifest records, a twin of
    // `variants.vcf` in another encoding (CRLF line ends, gzip), keeps its entry, so its
    // `reference_warning` holds and `reference_rows` lands in the manifest; a line that is
    // no VCF record (a space-delimited or an Ensembl-format re-encoding of `variants.vcf`)
    // is the manifest record at its ordinal, and refused when the file's record count is
    // not the manifest's; a line of a VCF input that is neither a twin nor a VCF record is
    // refused with its ordinal; a record the reference skipped in an input the manifest
    // does not record is refused; an unreadable input is refused with its path.
    let case = Case::new("own-input");
    let warned = "21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900";
    let warning = r#""reference_warning": "CPX is not a supported structural variant type""#;
    let held = list_of(&[record(RS1), record_with(warned, warning)]);
    let sets = [
        set("structural", r#"["--coding_only"]"#, Some("structural.vcf")),
        set("crlf", "[]", Some("variants.crlf.vcf")),
        set("missing", "[]", Some("absent.vcf")),
        set(
            "space",
            r#"["--delimiter", " "]"#,
            Some("variants.space.vcf"),
        ),
        set(
            "ensembl",
            r#"["--format", "ensembl"]"#,
            Some("variants.ensembl.txt"),
        ),
        set("gz", "[]", Some("variants.vcf.gz")),
        set(
            "short",
            r#"["--format", "ensembl"]"#,
            Some("one.ensembl.txt"),
        ),
    ];
    case.manifest(&format!(
        r#"{{"records": {held}, "flag_sets": [{}]}}"#,
        sets.join(", ")
    ));
    case.write_in_corpus(
        "structural.vcf",
        concat!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "21\t1000\tsv_a\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=1900\n",
            "21\t3000\tsv_b\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=3900\n",
            "21\t5000\tsv_c\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=5900\n"
        )
        .as_bytes(),
    );
    let both = "coding_sequence_variant,feature_truncation";
    let structural = row("sv_b", "21:3001-3900", "deletion", "ENST2", both);
    case.outputs(
        &structural,
        &(structural.clone()
            + &row(
                "sv_b",
                "21:3001-3900",
                "deletion",
                "ENST3",
                "feature_truncation",
            )),
    );
    let one_tuple = |class: &str| (1, format!("{{'{class}': 1}}"));
    let unexplained = one_tuple("unexplained_residual");
    assert_eq!(
        case.classify(&["--flag-set", "structural"]),
        clean(0, &counted(NO_ROWS, (unexplained.0, &unexplained.1), "[]"))
    );
    let after_structural = format!(
        "{}, {}}}",
        sets[0].strip_suffix('}').unwrap(),
        written(
            &[],
            &[vep_rs_only(
                &key("21:3001-3900", "deletion", transcript("ENST3"), &[1]),
                &["feature_truncation"],
                "unexplained_residual"
            )],
            ("{}", r#"{"unexplained_residual": 1}"#)
        )
    );
    assert_eq!(set_text(&case, "structural"), after_structural);
    assert_eq!(records_text(&case), held);

    // The twin input carries the manifest's records: the skipped one keeps its warning.
    case.write_in_corpus(
        "variants.crlf.vcf",
        concat!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\n",
            "21\t100\trs1\tA\tG\t.\tPASS\t.\r\n",
            "21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900\r\n"
        )
        .as_bytes(),
    );
    case.outputs(
        &rs1_row(),
        &(rs1_row() + &row("sv1", "21:501-900", "CPX", "ENST2", "feature_truncation")),
    );
    let skipped = one_tuple("reference_skipped_record");
    let skipped_line = clean(0, &counted(NO_ROWS, (skipped.0, &skipped.1), "[1]"));
    assert_eq!(case.classify(&["--flag-set", "crlf"]), skipped_line);
    let marked = list_of(&[
        record(RS1),
        record_with(warned, &format!(r#""reference_rows": 0, {warning}"#)),
    ]);
    assert_eq!(records_text(&case), marked);
    let twin_written = written(
        &[],
        &[vep_rs_only(
            &key("21:501-900", "CPX", transcript("ENST2"), &[1]),
            &["feature_truncation"],
            "reference_skipped_record",
        )],
        ("{}", r#"{"reference_skipped_record": 1}"#),
    );
    let after = |entry: &str| format!("{}, {twin_written}}}", entry.strip_suffix('}').unwrap());
    assert_eq!(set_text(&case, "crlf"), after(&sets[1]));

    // A space-delimited or Ensembl-format re-encoding of variants.vcf holds no VCF record
    // and takes the manifest's records in order; the gzipped copy is a twin read through
    // gzip.
    case.write_in_corpus(
        "variants.space.vcf",
        b"#CHROM POS ID REF ALT QUAL FILTER INFO\n21 100 rs1 A G . PASS .\n21 500 sv1 N <CPX> . PASS SVTYPE=CPX;END=900\n",
    );
    case.write_in_corpus(
        "variants.ensembl.txt",
        b"21\t100\t100\tA/G\t1\trs1\n21\t500\t900\tCPX\t1\tsv1\n",
    );
    case.write_in_corpus(
        "variants.vcf.gz",
        &gzip(concat!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "21\t100\trs1\tA\tG\t.\tPASS\t.\n",
            "21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900\n"
        )),
    );
    for (name, place) in [("space", 3), ("ensembl", 4), ("gz", 5)] {
        assert_eq!(case.classify(&["--flag-set", name]), skipped_line, "{name}");
        assert_eq!(set_text(&case, name), after(&sets[place]), "{name}");
    }
    assert_eq!(records_text(&case), marked);
    let corpus = &case.corpus;
    case.write_in_corpus("one.ensembl.txt", b"21\t100\t100\tA/G\t1\trs1\n");
    assert_eq!(
        case.classify(&["--flag-set", "short"]),
        failed(&format!(
            "{corpus}/one.ensembl.txt holds no VCF record and its 1 records are not the manifest's 2 re-encoded"
        ))
    );

    // A standalone VCF of short records (five or seven columns) numbers its own records.
    case.write_in_corpus(
        "own.short.vcf",
        b"#CHROM\tPOS\tID\tREF\tALT\n21\t7000\tsnv_a\tC\tT\n21\t8000\tsnv_b\tG\tA\t.\tPASS\n",
    );
    let own_short = set("own-short", "[]", Some("own.short.vcf"));
    let mut with_set = manifest(corpus);
    let held_sets = with_set.get_mut("flag_sets").and_then(Value::as_list_mut);
    held_sets.unwrap().push(Value::Object(object(&own_short)));
    case.manifest(&dumps_object(&with_set, &Style::PLAIN));
    let snv_a = row("snv_a", "21:7000", "T", "ENST7", "missense_variant");
    case.outputs(
        &(snv_a.clone() + &row("snv_b", "21:8000", "A", "ENST8", "synonymous_variant")),
        &(snv_a + &row("snv_b", "21:8000", "A", "ENST8", "missense_variant")),
    );
    let ran = case.classify(&["--flag-set", "own-short"]);
    let rows = (1, "{'unexplained_residual': 1}");
    assert_eq!(ran, clean(0, &counted(rows, NO_ROWS, "[]")));
    let expected = written(
        &[divergence(
            &key("21:8000", "A", transcript("ENST8"), &[1]),
            "synonymous_variant",
            &["missense_variant"],
            "unexplained_residual",
        )],
        &[],
        (r#"{"unexplained_residual": 1}"#, "{}"),
    );
    assert_eq!(
        set_text(&case, "own-short"),
        format!("{}, {expected}}}", own_short.strip_suffix('}').unwrap())
    );
    assert_eq!(records_text(&case), marked);
    let before = case.text();

    // A line of a VCF input that is neither a twin nor a VCF record is refused with its
    // ordinal.
    case.write_in_corpus(
        "own.short.vcf",
        b"21\t7000\tsnv_a\tC\tT\n21\t.\tsnv_x\tG\tA\n",
    );
    assert_eq!(
        case.classify(&["--flag-set", "own-short"]),
        failed(&format!(
            "{corpus}/own.short.vcf record 1 ('21\\t.\\tsnv_x\\tG\\tA') is no VCF record in a VCF input"
        ))
    );

    // A skipped record of an input the manifest does not record has no warning to carry.
    case.outputs(
        &structural,
        &(structural.clone()
            + &row(
                "sv_c",
                "21:5001-5900",
                "deletion",
                "ENST4",
                "feature_truncation",
            )),
    );
    assert_eq!(
        case.classify(&["--flag-set", "structural"]),
        failed(&rowless(2, "21 5000 sv_c N <DEL>"))
    );
    assert_eq!(
        case.classify(&["--flag-set", "missing"]),
        failed(&format!(
            "{corpus}/absent.vcf is not a readable text file ([Errno 2] No such file or directory: '{corpus}/absent.vcf')"
        ))
    );
    assert_eq!(case.text(), before);
}

/// A JSON list of the items given, as `json.dumps` writes it.
fn list_of(items: &[String]) -> String {
    format!("[{}]", items.join(", "))
}

#[test]
fn classify_refuses_a_record_the_frequency_filter_drops_on_the_reference_side_alone() {
    // Under `--check_frequency` Ensembl VEP drops records by frequency with no warning: a
    // record both engines drop is rowless by the run's nature, a record vep-rs alone keeps
    // is refused as a zero-row record without a warning, whether its FREQS list a
    // frequency of exactly 0 or another, and without the flag the dropped records are
    // refused the same way.
    let case = Case::new("frequency");
    let held = records(&[
        record(RS1),
        record("21\t200\trs2\tA\tT\t.\tPASS\t."),
        record("21\t300\trs3\tA\tG\t.\tPASS\t."),
        record("21\t400\trs4\tA\tC\t.\tPASS\t."),
    ]);
    let include = concat!(
        r#"{"name": "include", "flags": ["--check_frequency", "--freq_pop", "gnomADe_AFR", "#,
        r#""--freq_freq", "0.05", "--freq_gt_lt", "lt", "--freq_filter", "include"]}"#
    );
    case.manifest(&obj(&[&held, &format!(r#""flag_sets": [{include}]"#)]));
    let listed = |allele: &str, frequency: &str| {
        format!("IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:{allele}:{frequency}")
    };
    let intron = "intron_variant";
    let first = row_with("rs1", "21:100", "G", "ENST1", intron, &listed("G", "0.01"));
    // Both engines keep the first record and drop the other three.
    case.outputs(&first, &first);
    assert_eq!(
        case.classify(&["--flag-set", "include"]),
        clean(0, &counted(NO_ROWS, NO_ROWS, "[]"))
    );
    let after = obj(&[
        &held,
        &format!(
            r#""flag_sets": [{}, {}}}]"#,
            include.strip_suffix('}').unwrap(),
            nothing_written()
        ),
    ]);
    assert_eq!(case.after(), after);
    let written_once = case.text();

    // A record vep-rs alone keeps on a frequency of exactly 0 is refused.
    let kept_on_zero = first.clone()
        + &row_with("rs3", "21:300", "G", "ENST3", intron, &listed("G", "0"))
        + &row_with(
            "rs3",
            "21:300",
            "G",
            "ENST4",
            "upstream_gene_variant",
            "IMPACT=MODIFIER;DISTANCE=12;STRAND=1;FREQS=GNOMADE_AFR:G:0",
        );
    case.write("rs.txt", &kept_on_zero);
    assert_eq!(
        case.classify(&["--flag-set", "include"]),
        failed(&rowless(2, "21 300 rs3 A G"))
    );
    // So is one it alone keeps on another frequency, beside the first.
    case.write(
        "rs.txt",
        &(kept_on_zero.clone()
            + &row_with("rs4", "21:400", "C", "ENST5", intron, &listed("C", "0.3"))),
    );
    assert_eq!(
        case.classify(&["--flag-set", "include"]),
        problems(&[rowless(2, "21 300 rs3 A G"), rowless(3, "21 400 rs4 A C")])
    );
    case.write("rs.txt", &first);
    assert_eq!(
        case.classify(&[]),
        problems(&[
            rowless(1, "21 200 rs2 A T"),
            rowless(2, "21 300 rs3 A G"),
            rowless(3, "21 400 rs4 A C")
        ])
    );
    assert_eq!(case.text(), written_once);
}

#[test]
fn classify_reads_flags_and_entries_of_other_shapes_as_they_stand() {
    // The keys of a `flags` object are the run's flags and a `flags` text names none; a
    // `divergences` that is an object holds no entry and is written over, and one that
    // holds a number is refused on a line of the tool's own.
    let case = Case::new("other-shapes");
    let held = records(&[record(RS1), record("21\t200\trs2\tA\tT\t.\t.\t.")]);
    case.outputs(&rs1_row(), &rs1_row());
    let keyed = r#""flags": {"--coding_only": 1}"#;
    case.manifest(&obj(&[&held, keyed]));
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted(NO_ROWS, NO_ROWS, "[]"))
    );
    assert_eq!(case.after(), obj(&[&held, keyed, &nothing_written()]));
    case.manifest(&obj(&[&held, r#""flags": "--coding_only""#]));
    assert_eq!(case.classify(&[]), failed(&rowless(1, "21 200 rs2 A T")));
    case.manifest(&obj(&[&held, r#""flags": 7"#]));
    assert_eq!(
        case.classify(&[]),
        refused(
            1,
            "ERROR: [vep-golden-corpus] the manifest's flags is a number"
        )
    );
    let one = records(&[record(RS1)]);
    case.manifest(&obj(&[&one, r#""divergences": {"a": 1}"#]));
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted(NO_ROWS, NO_ROWS, "[]"))
    );
    assert_eq!(case.after(), obj(&[&one, &nothing_written()]));
    case.manifest(&obj(&[&one, r#""divergences": ["text", 5]"#]));
    assert_eq!(
        case.classify(&[]),
        refused(
            1,
            "ERROR: [vep-golden-corpus] divergences holds what is no entry"
        )
    );
}
