// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary on the rows both engines write for the three records of the
//! reproduction `merged_same_dbid_collapse`, with and without the converted cache of that
//! reproduction, on one thread and on the host's default, and holds the rows that stay
//! open and what the reports state of the class.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is held under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 2] = [Some("1"), None];

const OUTPUT_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
    Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\t\
    Existing_variation\tExtra\n";

const DISCORDANT_COLUMNS: &str =
    "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n";

/// The command line of the reference's header for the three records in one buffer.
const ONE_BUFFER: &str = "vep --assembly EXAMPLE1 --cache --cache_version 116 --database 0 \
    --dir_cache /cache --dir_plugins /plugins --force_overwrite --fork --format vcf \
    --input_file [PATH]/records.vcf --merged --no_stats --offline --output_file [PATH]/out.txt \
    --quiet --species two_set_example";

/// The same run with one record a buffer: Ensembl VEP writes a value of 1 as the bare option.
const ONE_RECORD_A_BUFFER: &str = "vep --assembly EXAMPLE1 --buffer_size --cache \
    --cache_version 116 --database 0 --dir_cache /cache --dir_plugins /plugins \
    --force_overwrite --fork --format vcf --input_file [PATH]/records.vcf --merged --no_stats \
    --offline --output_file [PATH]/out.txt --quiet --species two_set_example";

/// The same run in one process with two records a buffer.
const TWO_RECORDS_A_BUFFER: &str = "vep --assembly EXAMPLE1 --buffer_size 2 --cache \
    --cache_version 116 --database 0 --dir_cache /cache --dir_plugins /plugins \
    --force_overwrite --format vcf --input_file [PATH]/records.vcf --merged --no_stats \
    --offline --output_file [PATH]/out.txt --quiet --species two_set_example";

const SHARED_LOCUS_REFSEQ: &str = "shared_locus\tcontig:101500\tG\t900001\tXR_000001.1\t\
    Transcript\tnon_coding_transcript_exon_variant\t1500\t-\t-\t-\t-\t-\t\
    IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq\n";
const SHARED_LOCUS_ENSEMBL: &str = "shared_locus\tcontig:101500\tG\tENSEXG00000000001\t\
    ENSEXT00000000001\tTranscript\tnon_coding_transcript_exon_variant\t500\t-\t-\t-\t-\t-\t\
    IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl\n";
const FIRST_OF_FAR_PAIR: &str = "first_of_far_pair\tcontig:501000\tG\tENSEXG00000000002\t\
    ENSEXT00000000002\tTranscript\tnon_coding_transcript_exon_variant\t1000\t-\t-\t-\t-\t-\t\
    IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl\n";
const SECOND_OF_FAR_PAIR: &str = "second_of_far_pair\tcontig:3201000\tG\t900002\tXR_000002.1\t\
    Transcript\tnon_coding_transcript_exon_variant\t1000\t-\t-\t-\t-\t-\t\
    IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq\n";
const SECOND_OF_FAR_PAIR_INTERGENIC: &str = "second_of_far_pair\tcontig:3201000\tG\t-\t-\t-\t\
    intergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n";

const RECORDS: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    contig\t101500\tshared_locus\tA\tG\t.\t.\t.\n\
    contig\t501000\tfirst_of_far_pair\tA\tG\t.\t.\t.\n\
    contig\t3201000\tsecond_of_far_pair\tA\tG\t.\t.\t.\n";

/// The three rows of the class in `discordant.tsv`, the reference's before vep-rs's: the
/// reference's intergenic row and vep-rs's row of each dropped transcript.
const CLASS_ROWS: &str =
    "records.txt\tmissing_in_rust\tcontig:3201000\tG\t-\t-\tintergenic_variant\n\
    records.txt\textra_in_rust\tcontig:101500\tG\tENSEXT00000000001\tTranscript\t\
    non_coding_transcript_exon_variant\n\
    records.txt\textra_in_rust\tcontig:3201000\tG\tXR_000002.1\tTranscript\t\
    non_coding_transcript_exon_variant\n";

/// The converted cache of the reproduction, which the reproduction's own run reads.
fn converted_cache() -> String {
    let cache = Path::new(env!("CARGO_MANIFEST_DIR")).join(
        "../../scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/json_cache",
    );
    fs::canonicalize(cache)
        .expect("the reproduction's cache")
        .to_string_lossy()
        .into_owned()
}

/// A directory of one test under one thread bound, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, threads: Option<&str>) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-shared-dbid-{}-{tag}-t{}",
            std::process::id(),
            threads.unwrap_or("default")
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(fs::canonicalize(&path).unwrap())
    }

    fn arg(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }

    fn write(&self, name: &str, content: &str) {
        let path = self.0.join(name);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, content).unwrap();
    }

    fn read(&self, name: &str) -> String {
        let path = self.0.join(name);
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
    }

    /// Writes the two outputs and the input: the reference under `command` with
    /// `reference_rows`, and vep-rs's four transcript rows followed by `more_vep_rs_rows`.
    fn outputs(&self, command: &str, reference_rows: &[&str], more_vep_rs_rows: &[&str]) {
        let reference = format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n## VEP command-line: {command}\n\
             {OUTPUT_COLUMNS}{}",
            reference_rows.concat()
        );
        let vep_rs = format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n{OUTPUT_COLUMNS}{SHARED_LOCUS_ENSEMBL}\
             {SHARED_LOCUS_REFSEQ}{FIRST_OF_FAR_PAIR}{SECOND_OF_FAR_PAIR}{}",
            more_vep_rs_rows.concat()
        );
        self.write("reference/records.txt", &reference);
        self.write("vep_rs/records.txt", &vep_rs);
        self.write("records.vcf", RECORDS);
    }

    /// One compare run against the 116.3 reference, with `extra` after the paths.
    fn compare(&self, threads: Option<&str>, extra: &[&str]) -> Output {
        let (reference, vep_rs, report, records) = (
            self.arg("reference"),
            self.arg("vep_rs"),
            self.arg("report"),
            self.arg("records.vcf"),
        );
        let mut command = Command::new(BIN);
        command.current_dir(&self.0).args([
            "--perl-dir",
            &reference,
            "--rust-dir",
            &vep_rs,
            "--report-dir",
            &report,
            "--input-vcf",
            &records,
            "--no-canonical-contigs",
            "--reference-release",
            "116.3",
        ]);
        command.args(extra);
        match threads {
            Some(bound) => command.env(THREADS_ENV, bound),
            None => command.env_remove(THREADS_ENV),
        };
        command.output().unwrap()
    }

    fn open_rows(&self) -> String {
        self.read("report/discordant_open.tsv")
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The class's block of the summary's `excluded_categories`, which the block of the other
/// class of the converted cache follows.
fn class_block(count: u64, records: u64, perl: u64, rust: u64) -> String {
    format!(
        "        \"merged_same_dbid_collapse\": {{\n          \"count\": {count},\n          \
         \"excluded\": {count},\n          \"records\": {records},\n          \"perl\": {perl},\n          \
         \"rust\": {rust}\n        }},\n        \"merged_same_stable_id_collapse\": {{\n"
    )
}

#[test]
fn the_rows_of_the_class_are_set_aside_where_the_converted_cache_is_named() {
    let cache = converted_cache();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("named", threads);
        scratch.outputs(
            ONE_BUFFER,
            &[
                SHARED_LOCUS_REFSEQ,
                FIRST_OF_FAR_PAIR,
                SECOND_OF_FAR_PAIR_INTERGENIC,
            ],
            &[],
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.read("report/discordant.tsv"),
            format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
        );
        assert_eq!(scratch.open_rows(), DISCORDANT_COLUMNS);
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(3, 2, 1, 2)), "{summary}");
        assert!(
            summary.contains("\"excluded_perl_tuples\": 1,"),
            "{summary}"
        );
        assert!(
            summary.contains("\"excluded_rust_tuples\": 2,"),
            "{summary}"
        );
        assert!(
            summary.contains(
                "        \"merged_same_dbid_collapse\": {\n          \"definition\": \"Under \
                 --merged the reference keeps one transcript per dbID"
            ),
            "{summary}"
        );
        assert!(scratch
            .read("report/summary.md")
            .contains("| `merged_same_dbid_collapse` | sets aside | Under --merged"),);
        assert!(
            stdout(&output).contains(
                "adjusted: precision=1.000000 recall=1.000000 f1=1.000000 (excluded 1P + 2R = 3 \
                 intended-divergence tuples)\n"
            ),
            "{}",
            stdout(&output)
        );
    }
}

#[test]
fn without_the_cache_the_rows_stay_open_and_no_report_names_the_class() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("unnamed", threads);
        scratch.outputs(
            ONE_BUFFER,
            &[
                SHARED_LOCUS_REFSEQ,
                FIRST_OF_FAR_PAIR,
                SECOND_OF_FAR_PAIR_INTERGENIC,
            ],
            &[],
        );
        let output = scratch.compare(threads, &[]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
        );
        for report in ["report/summary.json", "report/summary.md"] {
            assert!(!scratch.read(report).contains("merged_same_dbid_collapse"));
        }
        assert!(
            !stdout(&output).contains("adjusted:"),
            "{}",
            stdout(&output)
        );
    }
}

#[test]
fn a_reference_that_read_one_record_a_buffer_dropped_the_near_pair_alone() {
    let cache = converted_cache();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("one-a-buffer", threads);
        scratch.outputs(
            ONE_RECORD_A_BUFFER,
            &[SHARED_LOCUS_REFSEQ, FIRST_OF_FAR_PAIR, SECOND_OF_FAR_PAIR],
            &[],
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(scratch.open_rows(), DISCORDANT_COLUMNS);
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(1, 1, 0, 1)), "{summary}");
    }
}

/// A buffer of its own never loads the other transcript of the far pair, so a reference
/// that wrote the intergenic row there differs for another reason: both rows stay open.
#[test]
fn a_row_the_replay_does_not_account_for_stays_open() {
    let cache = converted_cache();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("unaccounted", threads);
        scratch.outputs(
            ONE_RECORD_A_BUFFER,
            &[
                SHARED_LOCUS_REFSEQ,
                FIRST_OF_FAR_PAIR,
                SECOND_OF_FAR_PAIR_INTERGENIC,
            ],
            &[],
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            scratch.open_rows(),
            format!(
                "{DISCORDANT_COLUMNS}records.txt\tmissing_in_rust\tcontig:3201000\tG\t-\t-\t\
                 intergenic_variant\nrecords.txt\textra_in_rust\tcontig:3201000\tG\tXR_000002.1\t\
                 Transcript\tnon_coding_transcript_exon_variant\n"
            )
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(1, 1, 0, 1)), "{summary}");
    }
}

/// vep-rs writes a row too many under `--merged`: one on a transcript the cache does not
/// hold, and one the reference lacks on the transcript of the far pair its buffer kept.
/// Neither is a row of a dropped transcript, and the intergenic row beside the second is
/// no row of the class either.
#[test]
fn a_transcript_row_too_many_on_a_merged_cache_stays_open() {
    let cache = converted_cache();
    let unknown = "first_of_far_pair\tcontig:501000\tG\tENSEXG00000000009\tENSEXT00000000009\t\
        Transcript\tnon_coding_transcript_exon_variant\t7\t-\t-\t-\t-\t-\t\
        IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl\n";
    let kept_as_intergenic = "first_of_far_pair\tcontig:501000\tG\t-\t-\t-\tintergenic_variant\t\
        -\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n";
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("too-many", threads);
        scratch.outputs(
            ONE_BUFFER,
            &[
                SHARED_LOCUS_REFSEQ,
                kept_as_intergenic,
                SECOND_OF_FAR_PAIR_INTERGENIC,
            ],
            &[unknown],
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            format!(
                "{DISCORDANT_COLUMNS}records.txt\tmissing_in_rust\tcontig:501000\tG\t-\t-\t\
                 intergenic_variant\nrecords.txt\textra_in_rust\tcontig:501000\tG\t\
                 ENSEXT00000000002\tTranscript\tnon_coding_transcript_exon_variant\n\
                 records.txt\textra_in_rust\tcontig:501000\tG\tENSEXT00000000009\tTranscript\t\
                 non_coding_transcript_exon_variant\n"
            )
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(3, 2, 1, 2)), "{summary}");
    }
}

/// The buffer of the first two records drops the second transcript of the far pair and the
/// buffer of the third record keeps it: where the third shares the second's Location, the
/// reference wrote a row of the transcript there and the rows name no record, so vep-rs's row
/// beside the intergenic row, a row it writes too many there and the intergenic row stay
/// open. With the third record one base on, the reference wrote no row of the transcript at
/// the second's Location and the two rows there are set aside.
#[test]
fn a_row_of_a_transcript_the_reference_wrote_at_its_location_stays_open() {
    let cache = converted_cache();
    let transcript_row = |name: &str, position: &str, allele: &str| {
        format!(
            "{name}\tcontig:{position}\t{allele}\t900002\tXR_000002.1\tTranscript\t\
             non_coding_transcript_exon_variant\t1000\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq\n"
        )
    };
    let open_at_the_shared_location = format!(
        "{DISCORDANT_COLUMNS}records.txt\tmissing_in_rust\tcontig:3201000\tG\t-\t-\t\
         intergenic_variant\nrecords.txt\textra_in_rust\tcontig:3201000\tG\tXR_000002.1\t\
         Transcript\tnon_coding_transcript_exon_variant\nrecords.txt\textra_in_rust\t\
         contig:3201000\tT\tXR_000002.1\tTranscript\tnon_coding_transcript_exon_variant\n"
    );
    let open_one_base_on = format!(
        "{DISCORDANT_COLUMNS}records.txt\textra_in_rust\tcontig:3201001\tT\tXR_000002.1\t\
         Transcript\tnon_coding_transcript_exon_variant\n"
    );
    let cases = [
        (
            "3201000",
            open_at_the_shared_location,
            class_block(0, 0, 0, 0),
        ),
        ("3201001", open_one_base_on, class_block(2, 1, 1, 1)),
    ];
    for threads in THREAD_BOUNDS {
        for (position, open_rows, block) in &cases {
            let scratch = Scratch::new(&format!("matched-{position}"), threads);
            let third = transcript_row("third_record", position, "C");
            let too_many = transcript_row("third_record", position, "T");
            scratch.write(
                "reference/records.txt",
                &format!(
                    "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n## VEP command-line: \
                     {TWO_RECORDS_A_BUFFER}\n{OUTPUT_COLUMNS}{FIRST_OF_FAR_PAIR}\
                     {SECOND_OF_FAR_PAIR_INTERGENIC}{third}"
                ),
            );
            scratch.write(
                "vep_rs/records.txt",
                &format!(
                    "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n{OUTPUT_COLUMNS}\
                     {FIRST_OF_FAR_PAIR}{SECOND_OF_FAR_PAIR}{third}{too_many}"
                ),
            );
            scratch.write(
                "records.vcf",
                &format!(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     contig\t501000\tfirst_of_far_pair\tA\tG\t.\t.\t.\n\
                     contig\t3201000\tsecond_of_far_pair\tA\tG\t.\t.\t.\n\
                     contig\t{position}\tthird_record\tA\tC\t.\t.\t.\n"
                ),
            );
            let output = scratch.compare(threads, &["--json-cache", &cache]);
            assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
            assert_eq!(stderr(&output), "");
            assert_eq!(scratch.open_rows(), *open_rows, "{position}");
            let summary = scratch.read("report/summary.json");
            assert!(summary.contains(block), "{position}: {summary}");
        }
    }
}

#[test]
fn a_run_the_replay_does_not_follow_leaves_the_rows_open_and_says_why() {
    let cache = converted_cache();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("not-followed", threads);
        scratch.outputs(
            &format!("{ONE_BUFFER} --transcript_filter stable_id"),
            &[
                SHARED_LOCUS_REFSEQ,
                FIRST_OF_FAR_PAIR,
                SECOND_OF_FAR_PAIR_INTERGENIC,
            ],
            &[],
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            stderr(&output),
            "WARNING: [vep-outputs-compare] records.txt: merged_same_dbid_collapse is not applied: \
             the reference ran under \
             --transcript_filter, which the replay does not follow\n\
             WARNING: [vep-outputs-compare] records.txt: merged_same_stable_id_collapse is not \
             applied: the reference ran under --transcript_filter, which the replay does not \
             follow\n"
        );
        assert_eq!(
            scratch.open_rows(),
            format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(0, 0, 0, 0)), "{summary}");
    }
}

#[test]
fn a_cache_that_is_not_there_is_refused_and_the_older_reference_applies_no_such_class() {
    let scratch = Scratch::new("refused", None);
    scratch.outputs(
        ONE_BUFFER,
        &[
            SHARED_LOCUS_REFSEQ,
            FIRST_OF_FAR_PAIR,
            SECOND_OF_FAR_PAIR_INTERGENIC,
        ],
        &[],
    );
    let absent = scratch.arg("absent");
    let output = scratch.compare(None, &["--json-cache", &absent]);
    assert_eq!(output.status.code(), Some(2));
    assert!(
        stderr(&output).ends_with(&format!(
            ": error: --json-cache {absent} is not a directory\n"
        )),
        "{}",
        stderr(&output)
    );
    let cache = converted_cache();
    let output = scratch.compare(
        None,
        &["--json-cache", &cache, "--reference-release", "115.2"],
    );
    // Under 115.2 the cache gives the run's sequence names alone: its registry holds no
    // such class, the rows stay open and no report names the class.
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    assert_eq!(stderr(&output), "");
    assert_eq!(
        scratch.open_rows(),
        format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
    );
    for report in ["report/summary.json", "report/summary.md"] {
        assert!(!scratch.read(report).contains("merged_same_dbid_collapse"));
    }
}
