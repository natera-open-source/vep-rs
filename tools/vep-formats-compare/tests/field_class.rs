// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `field-class` command on planted files.
//!
//! `tests/data/field_class/` holds made-up score records as plain VCF text
//! (`snv.synthetic.vcf`, `indel.synthetic.vcf`), the gene symbols of the asked transcripts
//! (`transcripts.json`), questions about the rows of those records (`questions.tsv`) and, one
//! line a question, the answer each must get over the same records (`answers.txt`). The questions
//! are calls the golden corpus builder makes of `spliceai_multi_entry_cell` and questions
//! made from them that the class refuses on each of its grounds the records reach.
//!
//! The questions of `hgvs_shift_sequence_edge_window` are written by their test, beside a
//! made-up reference sequence and a candidate cache of three transcripts.

mod support;

use std::io::Write;
use std::path::{Path, PathBuf};
use support::{run, scratch};

const CLASS: &str = "spliceai_multi_entry_cell";
const EDGE_CLASS: &str = "hgvs_shift_sequence_edge_window";

fn data(name: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/field_class")
        .join(name);
    std::fs::read_to_string(path).expect("the file is read")
}

/// One BGZF block holding `data`.
fn block(data: &[u8]) -> Vec<u8> {
    let mut deflater =
        flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
    deflater.write_all(data).expect("the block is deflated");
    let body = deflater.finish().expect("the block is closed");
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let mut block = vec![
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
    ];
    let size = u16::try_from(body.len() + 25).expect("a block of less than 64 KiB");
    block.extend_from_slice(&size.to_le_bytes());
    block.extend_from_slice(&body);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    block.extend_from_slice(
        &u32::try_from(data.len())
            .expect("a short block")
            .to_le_bytes(),
    );
    block
}

/// Writes the VCF `text` as a BGZF file at `path`, its header lines one block and each
/// record line a block of its own, and beside it a tabix index of chromosome 21 whose one
/// bin holds one chunk over every record.
fn plant_scores(path: &Path, text: &str) {
    let (header, records): (Vec<&str>, Vec<&str>) = text
        .split_inclusive('\n')
        .partition(|line| line.starts_with('#'));
    let header = block(header.concat().as_bytes());
    let body: Vec<u8> = records
        .iter()
        .flat_map(|line| block(line.as_bytes()))
        .collect();
    let mut file = header.clone();
    file.extend_from_slice(&body);
    file.extend_from_slice(&block(b""));
    std::fs::write(path, file).expect("the score file is written");
    let mut index = b"TBI\x01".to_vec();
    // One sequence, a generic layout whose sequence is column 1 and whose position is
    // column 2, `#` before a header line, no line skipped, a name of three bytes.
    for number in [1i32, 0, 1, 2, 2, i32::from(b'#'), 0, 3] {
        index.extend_from_slice(&number.to_le_bytes());
    }
    index.extend_from_slice(b"21\0");
    // One bin, the one that spans the sequence, with one chunk; no linear index.
    index.extend_from_slice(&1i32.to_le_bytes());
    index.extend_from_slice(&0u32.to_le_bytes());
    index.extend_from_slice(&1i32.to_le_bytes());
    index.extend_from_slice(&((header.len() as u64) << 16).to_le_bytes());
    index.extend_from_slice(&(((header.len() + body.len()) as u64) << 16).to_le_bytes());
    index.extend_from_slice(&0i32.to_le_bytes());
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(&index).expect("the index is deflated");
    let mut beside = path.as_os_str().to_owned();
    beside.push(".tbi");
    std::fs::write(beside, gzip.finish().expect("the index is closed"))
        .expect("the index is written");
}

/// The planted files of one test: the two score files, the candidate cache and the
/// questions.
struct Planted {
    directory: PathBuf,
    snv: PathBuf,
    indel: PathBuf,
    cache: PathBuf,
    questions: PathBuf,
    answers: PathBuf,
}

fn planted(name: &str) -> Planted {
    let directory = scratch(name);
    let (snv, indel) = (
        directory.join("snv.synthetic.vcf.gz"),
        directory.join("indel.synthetic.vcf.gz"),
    );
    plant_scores(&snv, &data("snv.synthetic.vcf"));
    plant_scores(&indel, &data("indel.synthetic.vcf"));
    let cache = directory.join("cache");
    let shards = cache.join("transcripts").join("21");
    std::fs::create_dir_all(&shards).expect("the cache is made");
    std::fs::write(shards.join("1-1000000.json"), data("transcripts.json"))
        .expect("the shard is written");
    let questions = directory.join("questions.tsv");
    std::fs::write(&questions, data("questions.tsv")).expect("the questions are written");
    let answers = directory.join("answers.txt");
    Planted {
        directory,
        snv,
        indel,
        cache,
        questions,
        answers,
    }
}

impl Planted {
    /// The command line that asks the planted questions, each path as a text.
    fn arguments(&self) -> Vec<String> {
        let text = |path: &Path| path.to_str().expect("a path of text").to_string();
        vec![
            "field-class".to_string(),
            "--class".to_string(),
            CLASS.to_string(),
            "--spliceai-table".to_string(),
            text(&self.snv),
            "--spliceai-table".to_string(),
            text(&self.indel),
            "--json-cache".to_string(),
            text(&self.cache),
            "--questions".to_string(),
            text(&self.questions),
            "--answers".to_string(),
            text(&self.answers),
        ]
    }

    /// The command line with the value of `option` replaced by `value`.
    fn with(&self, option: &str, value: &Path) -> Vec<String> {
        let mut arguments = self.arguments();
        let at = arguments
            .iter()
            .position(|argument| argument == option)
            .expect("the option is on the command line");
        arguments[at + 1] = value.to_str().expect("a path of text").to_string();
        arguments
    }
}

/// The error line the command ends on.
fn refusal(what: &str) -> String {
    format!("ERROR: [compare_vep_formats] {what}\n")
}

#[test]
fn every_sampled_question_gets_the_answer_stated_for_it() {
    let files = planted("field-class-sample");
    let ran = run(&files.directory, &files.arguments());
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (0, String::new(), String::new())
    );
    let answers = std::fs::read_to_string(&files.answers).expect("the answers are written");
    let stated = data("answers.txt");
    assert_eq!(answers, stated);
    let count = |answer: &str| stated.lines().filter(|line| *line == answer).count();
    assert_eq!(
        (count("1"), count("0"), stated.lines().count()),
        (74, 150, 224)
    );
    assert_eq!(data("questions.tsv").lines().count(), 224);
}

#[test]
fn questions_without_a_line_are_answered_by_an_empty_file() {
    let files = planted("field-class-empty");
    std::fs::write(&files.questions, "").expect("the questions are written");
    let ran = run(&files.directory, &files.arguments());
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (0, String::new(), String::new())
    );
    assert_eq!(
        std::fs::read_to_string(&files.answers).expect("the answers are written"),
        ""
    );
}

#[test]
fn a_table_or_a_cache_that_cannot_be_read_ends_the_run_and_writes_no_answers() {
    let files = planted("field-class-unread");
    let shown = |path: &Path| path.to_str().expect("a path of text").to_string();
    let missing = files.directory.join("none.vcf.gz");
    let seven = files.directory.join("seven.vcf.gz");
    plant_scores(
        &seven,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n21\t100\t.\tA\tC\t.\t.\n",
    );
    let empty = files.directory.join("empty");
    std::fs::create_dir(&empty).expect("the directory is made");
    for (arguments, what) in [
        (
            files.with("--spliceai-table", &missing),
            format!(
                "--spliceai-table: cannot read the tabix index {}.tbi: No such file or directory",
                shown(&missing)
            ),
        ),
        (
            files.with("--spliceai-table", &seven),
            format!(
                "--spliceai-table: {} is not a VCF (a header line of eight columns)",
                shown(&seven)
            ),
        ),
        (
            files.with("--json-cache", &empty),
            format!(
                "--json-cache: {} has no transcripts/ directory",
                shown(&empty)
            ),
        ),
        (
            files.with("--questions", &missing),
            format!("cannot read {}: No such file or directory", shown(&missing)),
        ),
    ] {
        let ran = run(&files.directory, &arguments);
        assert_eq!(
            (ran.code, ran.out(), ran.err()),
            (2, String::new(), refusal(&what))
        );
        assert!(!files.answers.exists());
    }
}

#[test]
fn a_question_of_another_number_of_columns_ends_the_run_and_writes_no_answers() {
    let files = planted("field-class-columns");
    let questions = data("questions.tsv");
    let mut lines: Vec<&str> = questions.lines().collect();
    let short = lines[1].rsplit_once('\t').expect("a line of columns").0;
    lines[1] = short;
    std::fs::write(&files.questions, lines.join("\n") + "\n").expect("the questions are written");
    let ran = run(&files.directory, &files.arguments());
    let what = format!(
        "--questions: {} line 2 has 29 columns, a question of {CLASS} has 30",
        files.questions.to_str().expect("a path of text")
    );
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (2, String::new(), refusal(&what))
    );
    assert!(!files.answers.exists());
}

#[test]
fn answers_that_cannot_be_written_end_the_run() {
    let files = planted("field-class-unwritten");
    let below = files.directory.join("none").join("answers.txt");
    let ran = run(&files.directory, &files.with("--answers", &below));
    let what = format!(
        "cannot write {}: No such file or directory",
        below.to_str().expect("a path of text")
    );
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (2, String::new(), refusal(&what))
    );
}

#[test]
fn the_command_prints_its_help_and_refuses_a_command_line_it_cannot_read() {
    let directory = scratch("field-class-usage");
    let ran = run(&directory, &["field-class", "--help"]);
    assert_eq!((ran.code, ran.err()), (0, String::new()));
    let help = ran.out();
    let usage = "usage: vep-formats-compare field-class [-h]\n";
    assert!(help.starts_with(usage), "{help}");
    for option in [
        "--class {spliceai_multi_entry_cell,hgvs_shift_sequence_edge_window}",
        "--spliceai-table FILE",
        "--reference-fasta FILE",
        "--json-cache DIR",
        "--questions IN",
        "--answers OUT",
    ] {
        assert!(help.contains(option), "{option}");
    }
    let block = help
        .split_once("\n\n")
        .expect("a usage block and the help below it")
        .0;
    let error =
        |message: &str| format!("{block}\nvep-formats-compare field-class: error: {message}\n");
    for (arguments, message) in [
        (
            vec!["field-class"],
            "the following arguments are required: --class, --questions, --answers",
        ),
        (
            vec!["field-class", "--class", CLASS, "--questions", "q.tsv"],
            "the following arguments are required: --spliceai-table, --answers",
        ),
        (
            vec!["field-class", "--class", EDGE_CLASS, "--answers", "a.txt"],
            "the following arguments are required: --reference-fasta, --json-cache, --questions",
        ),
        (
            vec![
                "field-class",
                "--class",
                EDGE_CLASS,
                "--spliceai-table",
                "s.vcf.gz",
                "--reference-fasta",
                "r.fa",
                "--questions",
                "q.tsv",
                "--answers",
                "a.txt",
            ],
            "the following arguments are required: --json-cache",
        ),
        (
            vec!["field-class", "--class", "splice_family_swap"],
            "argument --class: invalid choice: 'splice_family_swap' (choose from 'spliceai_multi_entry_cell', 'hgvs_shift_sequence_edge_window')",
        ),
        (
            vec!["field-class", "--class", CLASS, "--format", "default"],
            "unrecognized arguments: --format default",
        ),
    ] {
        let ran = run(&directory, &arguments);
        assert_eq!((ran.code, ran.out(), ran.err()), (2, String::new(), error(message)));
    }
}

/// The made-up reference and candidate cache of the sequence-edge questions. Region `edge`
/// is forty bases whose bases 10 to 13 are one base repeated between two runs of another;
/// region `inner` is the same bases with 1,000 more after them. The cache holds a
/// forward-strand transcript over bases 5 to 25 of each region and a reverse-strand one over
/// the same bases of `edge`, each of version 1 and without a coding sequence.
struct Edges {
    directory: PathBuf,
    fasta: PathBuf,
    cache: PathBuf,
    questions: PathBuf,
    answers: PathBuf,
}

fn planted_edges(name: &str) -> Edges {
    let directory = scratch(name);
    let forty = "C".repeat(9) + &"A".repeat(4) + &"C".repeat(27);
    let regions = [
        ("edge", forty.clone()),
        ("inner", forty + &"C".repeat(1000)),
    ];
    let (mut text, mut index) = (String::new(), String::new());
    for (region, bases) in &regions {
        text += &format!(">{region}\n");
        let length = bases.len();
        index += &format!(
            "{region}\t{length}\t{}\t{length}\t{}\n",
            text.len(),
            length + 1
        );
        text += &format!("{bases}\n");
    }
    let fasta = directory.join("reference.fa");
    std::fs::write(&fasta, text).expect("the reference is written");
    std::fs::write(directory.join("reference.fa.fai"), index).expect("the index is written");
    let cache = directory.join("cache");
    let record = |name: &str, strand: i64| {
        format!(
            r#"{{"stable_id": "{name}", "version": 1, "start": 5, "end": 25, "strand": {strand}, "biotype": "lncRNA", "exons": [{{"start": 5, "end": 25}}]}}"#
        )
    };
    for (region, records) in [
        ("edge", vec![record("ENSTF", 1), record("ENSTR", -1)]),
        ("inner", vec![record("ENSTI", 1)]),
    ] {
        let shards = cache.join("transcripts").join(region);
        std::fs::create_dir_all(&shards).expect("the cache is made");
        std::fs::write(
            shards.join("1-1000000.json"),
            format!("[{}]", records.join(", ")),
        )
        .expect("the shard is written");
    }
    Edges {
        questions: directory.join("questions.tsv"),
        answers: directory.join("answers.txt"),
        directory,
        fasta,
        cache,
    }
}

impl Edges {
    /// The command line that asks the questions' file, each path as a text.
    fn arguments(&self) -> Vec<String> {
        let text = |path: &Path| path.to_str().expect("a path of text").to_string();
        vec![
            "field-class".to_string(),
            format!("--class={EDGE_CLASS}"),
            "--reference-fasta".to_string(),
            text(&self.fasta),
            "--json-cache".to_string(),
            text(&self.cache),
            "--questions".to_string(),
            text(&self.questions),
            "--answers".to_string(),
            text(&self.answers),
        ]
    }
}

/// One row of a sequence-edge question: where it stands, its Feature_type and strand, each
/// side's notation body and offset, and each side's shift length and cDNA position, an empty
/// text for one the side does not write.
#[derive(Clone, Copy)]
struct EdgeRow<'a> {
    location: &'a str,
    feature: &'a str,
    kind: &'a str,
    strand: &'a str,
    reference: (&'a str, &'a str),
    candidate: (&'a str, &'a str),
    moved: ((&'a str, &'a str), (&'a str, &'a str)),
}

impl EdgeRow<'_> {
    /// The question of the row's field `name`: the twenty-four columns of the class, a
    /// deleted allele, no protein notation and no coding or protein position on either row
    /// and no strand on the candidate's.
    fn question(&self, name: &str) -> String {
        let cell = |text: &str| if text.is_empty() { "-" } else { text }.to_string();
        let notation = |body: &str| match body {
            "" => "-".to_string(),
            body => format!("{}.1:{body}", self.feature),
        };
        let side = |written: (&str, &str), moved: (&str, &str), strand: &str| {
            [
                self.kind.to_string(),
                cell(strand),
                notation(written.0),
                "-".to_string(),
                cell(written.1),
                cell(moved.0),
                cell(moved.1),
                "-".to_string(),
                "-".to_string(),
            ]
        };
        let values = match name {
            "HGVSc" => (notation(self.reference.0), notation(self.candidate.0)),
            "HGVS_OFFSET" => (cell(self.reference.1), cell(self.candidate.1)),
            "SHIFT_LENGTH" => (cell(self.moved.0 .0), cell(self.moved.1 .0)),
            "cDNA_position" => (cell(self.moved.0 .1), cell(self.moved.1 .1)),
            // A start the JSON format keeps, which the candidate's object alone holds.
            "cds_start" => ("-".to_string(), "7".to_string()),
            _ => ("p.one".to_string(), "p.other".to_string()),
        };
        let mut cells = vec![
            self.location.to_string(),
            "-".to_string(),
            self.feature.to_string(),
            name.to_string(),
            values.0,
            values.1,
        ];
        cells.extend(side(self.reference, self.moved.0, self.strand));
        cells.extend(side(self.candidate, self.moved.1, ""));
        cells.join("\t") + "\n"
    }
}

#[test]
fn the_sequence_edge_class_answers_for_the_rows_its_two_walks_explain() {
    // A deleted base at 10 of `edge`, the first of four alike, on the forward-strand
    // transcript: the sequence admits a shift of three bases, and the walk over the cut of
    // a slice clipped at the region's last base, forty bases on, admits none. Ensembl VEP
    // so writes the deletion at its sixth base with no offset and vep-rs at its ninth with
    // an offset of 3. Under `--shift_3prime 1` each engine prints the length of its walk
    // and the base's place where the walk leaves it: none and the sixth, three and the ninth.
    let files = planted_edges("field-class-edges");
    let taken = EdgeRow {
        location: "edge:10",
        feature: "ENSTF",
        kind: "Transcript",
        strand: "1",
        reference: ("n.6del", ""),
        candidate: ("n.9del", "3"),
        moved: (("0", "6"), ("3", "9")),
    };
    let unwritten = EdgeRow {
        reference: ("", ""),
        candidate: ("", ""),
        ..taken
    };
    let asked = [
        // The notation and the offset of the row are the class's.
        (taken.question("HGVSc"), "1"),
        (taken.question("HGVS_OFFSET"), "1"),
        // So are the length of each walk and the position it leaves, beside the notations
        // and on a row that writes none.
        (taken.question("SHIFT_LENGTH"), "1"),
        (taken.question("cDNA_position"), "1"),
        (unwritten.question("SHIFT_LENGTH"), "1"),
        (unwritten.question("cDNA_position"), "1"),
        // A length that is not the candidate's walk is not, nor is a position beside it,
        // nor a position at another base than the walk's on a row that prints no length.
        (
            EdgeRow {
                moved: (("0", "6"), ("2", "9")),
                ..unwritten
            }
            .question("SHIFT_LENGTH"),
            "0",
        ),
        (
            EdgeRow {
                moved: (("0", "6"), ("2", "9")),
                ..unwritten
            }
            .question("cDNA_position"),
            "0",
        ),
        (
            EdgeRow {
                moved: (("", "6"), ("", "8")),
                ..unwritten
            }
            .question("cDNA_position"),
            "0",
        ),
        // A question's rows carry the default format's columns and no start or end of the
        // JSON format: a coding start asked beside them stays a difference, though neither
        // walk gives a transcript without a coding sequence a coding position.
        (taken.question("cds_start"), "0"),
        // A protein notation on a transcript the cache gives no coding sequence stays a
        // difference.
        (taken.question("HGVSp"), "0"),
        // The same values on a row of another feature type, on the other strand's
        // transcript, with another offset on the candidate's side, with the candidate's
        // notation at the unshifted base, and 1,000 bases from the region's last base.
        (
            EdgeRow {
                kind: "RegulatoryFeature",
                ..taken
            }
            .question("HGVSc"),
            "0",
        ),
        (
            EdgeRow {
                feature: "ENSTR",
                strand: "-1",
                ..taken
            }
            .question("HGVSc"),
            "0",
        ),
        (
            EdgeRow {
                candidate: ("n.9del", "2"),
                ..taken
            }
            .question("HGVS_OFFSET"),
            "0",
        ),
        (
            EdgeRow {
                candidate: ("n.6del", "3"),
                ..taken
            }
            .question("HGVS_OFFSET"),
            "0",
        ),
        (
            EdgeRow {
                location: "inner:10",
                feature: "ENSTI",
                ..taken
            }
            .question("HGVSc"),
            "0",
        ),
    ];
    let questions: String = asked
        .iter()
        .map(|(question, _)| question.as_str())
        .collect();
    assert!(questions.lines().all(|line| line.split('\t').count() == 24));
    std::fs::write(&files.questions, questions).expect("the questions are written");
    let ran = run(&files.directory, &files.arguments());
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (0, String::new(), String::new())
    );
    let stated: String = asked
        .iter()
        .map(|(_, answer)| format!("{answer}\n"))
        .collect();
    assert_eq!(
        std::fs::read_to_string(&files.answers).expect("the answers are written"),
        stated
    );
}

#[test]
fn a_reference_or_a_cache_the_sequence_edge_class_cannot_read_ends_the_run() {
    let files = planted_edges("field-class-edges-unread");
    std::fs::write(&files.questions, "").expect("the questions are written");
    let shown = |path: &Path| path.to_str().expect("a path of text").to_string();
    let with = |option: &str, value: &Path| {
        let mut arguments = files.arguments();
        let at = arguments
            .iter()
            .position(|argument| argument == option)
            .expect("the option is on the command line");
        arguments[at + 1] = shown(value);
        arguments
    };
    let missing = files.directory.join("none.fa");
    let packed = files.directory.join("reference.fa.gz");
    let empty = files.directory.join("empty");
    std::fs::create_dir(&empty).expect("the directory is made");
    let short = planted_edges("field-class-edges-columns");
    std::fs::write(&short.questions, "edge:10\t-\tENSTF\tHGVSc\ta\tb\n")
        .expect("the questions are written");
    for (arguments, what) in [
        (
            with("--reference-fasta", &missing),
            format!(
                "--reference-fasta: cannot read {}.fai: No such file or directory",
                shown(&missing)
            ),
        ),
        (
            with("--reference-fasta", &packed),
            format!(
                "--reference-fasta: {} is compressed; an uncompressed FASTA beside its .fai is read",
                shown(&packed)
            ),
        ),
        (
            with("--json-cache", &empty),
            format!(
                "--json-cache: {} has no transcripts/ directory",
                shown(&empty)
            ),
        ),
        (
            short.arguments(),
            format!(
                "--questions: {} line 1 has 6 columns, a question of {EDGE_CLASS} has 24",
                shown(&short.questions)
            ),
        ),
    ] {
        let ran = run(&files.directory, &arguments);
        assert_eq!(
            (ran.code, ran.out(), ran.err()),
            (2, String::new(), refusal(&what))
        );
        assert!(!files.answers.exists() && !short.answers.exists());
    }
    // The command line of the class without a question reads both and answers nothing.
    let ran = run(&files.directory, &files.arguments());
    assert_eq!(
        (ran.code, ran.out(), ran.err()),
        (0, String::new(), String::new())
    );
    assert_eq!(
        std::fs::read_to_string(&files.answers).expect("the answers are written"),
        ""
    );
}

#[test]
fn the_command_is_its_own_as_the_first_argument_alone() {
    let directory = scratch("field-class-first");
    let ran = run(&directory, &["--format", "default", "field-class"]);
    assert_eq!((ran.code, ran.out()), (2, String::new()));
    let error = ran.err();
    assert!(
        error.starts_with("usage: vep-formats-compare [-h]"),
        "{error}"
    );
    assert!(
        error.ends_with("vep-formats-compare: error: unrecognized arguments: field-class\n"),
        "{error}"
    );
}
