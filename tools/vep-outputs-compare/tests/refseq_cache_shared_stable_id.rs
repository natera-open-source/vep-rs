// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary on the rows both engines write under `--refseq` for four records of a
//! made-up converted cache in which four stable ids each name two transcripts, an
//! `ensembl`-sourced one and one of another source, and holds the rows that stay open and
//! what the reports state of the class `merged_same_stable_id_collapse`, on one thread and on
//! the host's default.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is held under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 2] = [Some("1"), None];

const CLASS: &str = "merged_same_stable_id_collapse";

const OUTPUT_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
    Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\t\
    Existing_variation\tExtra\n";

const DISCORDANT_COLUMNS: &str =
    "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n";

/// The command line of the reference's header for the four records in one buffer.
const ONE_BUFFER: &str = "vep --assembly EXAMPLE1 --cache --cache_version 116 --database 0 \
    --dir_cache /cache --dir_plugins /plugins --force_overwrite --input_file [PATH]/records.vcf \
    --no_stats --offline --output_file [PATH]/out.txt --refseq --species one_id_example";

/// The same run with one record a buffer: Ensembl VEP writes a value of 1 as the bare option.
const ONE_RECORD_A_BUFFER: &str = "vep --assembly EXAMPLE1 --buffer_size --cache \
    --cache_version 116 --database 0 --dir_cache /cache --dir_plugins /plugins \
    --force_overwrite --input_file [PATH]/records.vcf --no_stats --offline \
    --output_file [PATH]/out.txt --refseq --species one_id_example";

const INFO: &str = "{\"assembly\": \"EXAMPLE1\", \"cache_version\": 116, \
    \"contig_order\": [\"contig\"], \"source_type\": \"refseq\", \"species\": \"one_id_example\"}\n";

/// One transcript of a region file: its stable id, `dbID`, source, start and end.
type Listed = (&'static str, u32, &'static str, u32, u32);

/// The region files of the cache, each with its transcripts. Under each of the four ids the `ensembl`-sourced transcript stands in one
/// file and the other in another: for the third id and for the first two in the file after
/// it, with the same start; for the fourth five files on, at another locus.
const REGION_FILES: [(&str, &[Listed]); 6] = [
    (
        "2000001-3000000",
        &[("XM_900003.1", 103, "ensembl", 2_998_588, 2_999_361)],
    ),
    (
        "3000001-4000000",
        &[("XM_900003.1", 203, "Gnomon", 2_998_588, 3_017_103)],
    ),
    (
        "4000001-5000000",
        &[
            ("XM_900004.1", 204, "Gnomon", 4_686_346, 4_984_780),
            ("XM_900005.1", 205, "Gnomon", 4_900_000, 4_950_000),
        ],
    ),
    (
        "5000001-6000000",
        &[
            ("XM_900001.1", 101, "ensembl", 5_729_692, 5_820_131),
            ("XM_900002.1", 102, "ensembl", 5_730_751, 5_820_131),
        ],
    ),
    (
        "6000001-7000000",
        &[
            ("XM_900001.1", 201, "Gnomon", 5_729_692, 6_007_324),
            ("XM_900002.1", 202, "Gnomon", 5_730_751, 6_007_324),
        ],
    ),
    (
        "9000001-10000000",
        &[("XM_900004.1", 104, "ensembl", 9_095_991, 9_190_778)],
    ),
];

/// The four records. The third lies past the end of the `ensembl`-sourced transcripts of the
/// first two ids, which the buffer of all four holds, since they lie between its lowest and
/// its highest record.
const RECORDS: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    contig\t2999291\tinside_both_copies\tA\tG\t.\t.\t.\n\
    contig\t4919323\tat_one_locus\tA\tG\t.\t.\t.\n\
    contig\t5996672\tpast_the_kept_copies\tA\tG\t.\t.\t.\n\
    contig\t9100000\tat_the_other_locus\tA\tG\t.\t.\t.\n";

/// One transcript row of an output.
fn transcript_row(name: &str, position: &str, feature: &str, consequence: &str) -> String {
    format!(
        "{name}\tcontig:{position}\tG\t{feature}\t{feature}\tTranscript\t{consequence}\t-\t-\t-\t\
         -\t-\t-\tIMPACT=MODIFIER;STRAND=-1\n"
    )
}

/// The reference's row of the record that its pass left with no transcript.
const PAST_THE_KEPT_COPIES_INTERGENIC: &str = "past_the_kept_copies\tcontig:5996672\tG\t-\t-\t-\t\
    intergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n";

/// The rows the reference writes for the four records in one buffer: the row of the
/// `ensembl`-sourced transcript the third id keeps, the row of the one transcript of the
/// fifth id, the intergenic row, and the row of the `ensembl`-sourced transcript of the
/// fourth id at its own locus.
fn reference_rows_of_one_buffer() -> Vec<String> {
    vec![
        transcript_row(
            "inside_both_copies",
            "2999291",
            "XM_900003.1",
            "non_coding_transcript_exon_variant",
        ),
        transcript_row("at_one_locus", "4919323", "XM_900005.1", "intron_variant"),
        PAST_THE_KEPT_COPIES_INTERGENIC.to_string(),
        transcript_row(
            "at_the_other_locus",
            "9100000",
            "XM_900004.1",
            "non_coding_transcript_exon_variant",
        ),
    ]
}

/// The rows vep-rs writes: a row of every transcript within reach of each record.
fn vep_rs_rows() -> Vec<String> {
    vec![
        transcript_row(
            "inside_both_copies",
            "2999291",
            "XM_900003.1",
            "non_coding_transcript_exon_variant",
        ),
        transcript_row(
            "inside_both_copies",
            "2999291",
            "XM_900003.1",
            "missense_variant",
        ),
        transcript_row("at_one_locus", "4919323", "XM_900004.1", "intron_variant"),
        transcript_row("at_one_locus", "4919323", "XM_900005.1", "intron_variant"),
        transcript_row(
            "past_the_kept_copies",
            "5996672",
            "XM_900001.1",
            "stop_gained",
        ),
        transcript_row(
            "past_the_kept_copies",
            "5996672",
            "XM_900002.1",
            "stop_gained",
        ),
        transcript_row(
            "at_the_other_locus",
            "9100000",
            "XM_900004.1",
            "non_coding_transcript_exon_variant",
        ),
    ]
}

/// The five rows of the class in `discordant.tsv`, the reference's before vep-rs's: the
/// reference's intergenic row and vep-rs's row of each dropped transcript.
const CLASS_ROWS: &str =
    "records.txt\tmissing_in_rust\tcontig:5996672\tG\t-\t-\tintergenic_variant\n\
    records.txt\textra_in_rust\tcontig:2999291\tG\tXM_900003.1\tTranscript\tmissense_variant\n\
    records.txt\textra_in_rust\tcontig:4919323\tG\tXM_900004.1\tTranscript\tintron_variant\n\
    records.txt\textra_in_rust\tcontig:5996672\tG\tXM_900001.1\tTranscript\tstop_gained\n\
    records.txt\textra_in_rust\tcontig:5996672\tG\tXM_900002.1\tTranscript\tstop_gained\n";

/// A directory of one test under one thread bound, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, threads: Option<&str>) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-shared-stable-id-{}-{tag}-t{}",
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

    /// Writes the converted cache, in which `source_of` gives each transcript's source from
    /// the one `REGION_FILES` states, and returns its path.
    fn cache(&self, source_of: impl Fn(&str) -> String) -> String {
        self.write("json_cache/info.json", INFO);
        for (span, transcripts) in REGION_FILES {
            let objects: Vec<String> = transcripts
                .iter()
                .map(|(stable_id, db_id, source, start, end)| {
                    format!(
                        "{{\"stable_id\": \"{stable_id}\", \"dbID\": {db_id}, \
                         \"source\": \"{}\", \"start\": {start}, \"end\": {end}, \"strand\": -1}}",
                        source_of(source)
                    )
                })
                .collect();
            self.write(
                &format!("json_cache/transcripts/contig/{span}.json"),
                &format!("[{}]\n", objects.join(", ")),
            );
        }
        self.arg("json_cache")
    }

    /// Writes the two outputs and the input: the reference under `command`.
    fn outputs(&self, command: &str, reference_rows: &[String], vep_rs_rows: &[String]) {
        let reference = format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n## VEP command-line: {command}\n\
             {OUTPUT_COLUMNS}{}",
            reference_rows.concat()
        );
        let vep_rs = format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n{OUTPUT_COLUMNS}{}",
            vep_rs_rows.concat()
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

/// The sources as `REGION_FILES` states them.
fn as_stated(source: &str) -> String {
    source.to_string()
}

/// The class's block of the summary's `excluded_categories`.
fn class_block(count: u64, records: u64, perl: u64, rust: u64) -> String {
    format!(
        "        \"{CLASS}\": {{\n          \"count\": {count},\n          \
         \"excluded\": {count},\n          \"records\": {records},\n          \"perl\": {perl},\n          \
         \"rust\": {rust}\n        }}\n"
    )
}

/// The three shapes the class takes: a record past the kept transcripts of two ids, whose
/// only rows are vep-rs's beside the reference's intergenic row; a record inside both
/// transcripts of an id, where the reference wrote the kept one's row and vep-rs that row and
/// the dropped one's; and a record at one of two loci of an id, beside a row both engines
/// wrote of another transcript.
#[test]
fn the_rows_of_the_class_are_set_aside_where_the_converted_cache_is_named() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("named", threads);
        let cache = scratch.cache(as_stated);
        scratch.outputs(ONE_BUFFER, &reference_rows_of_one_buffer(), &vep_rs_rows());
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.read("report/discordant.tsv"),
            format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
        );
        assert_eq!(scratch.open_rows(), DISCORDANT_COLUMNS);
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(5, 3, 1, 4)), "{summary}");
        assert!(
            summary.contains("\"excluded_perl_tuples\": 1,"),
            "{summary}"
        );
        assert!(
            summary.contains("\"excluded_rust_tuples\": 4,"),
            "{summary}"
        );
        assert!(
            summary.contains(&format!(
                "        \"{CLASS}\": {{\n          \"definition\": \"Under --refseq or --merged \
                 the reference keeps one transcript per stable identifier"
            )),
            "{summary}"
        );
        assert!(scratch.read("report/summary.md").contains(&format!(
            "| `{CLASS}` | sets aside | Under --refseq or --merged"
        )));
        assert!(
            stdout(&output).contains(
                "adjusted: precision=1.000000 recall=1.000000 f1=1.000000 (excluded 1P + 4R = 5 \
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
        scratch.outputs(ONE_BUFFER, &reference_rows_of_one_buffer(), &vep_rs_rows());
        let output = scratch.compare(threads, &[]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
        );
        for report in ["report/summary.json", "report/summary.md"] {
            assert!(!scratch.read(report).contains(CLASS));
        }
    }
}

/// Where no transcript of an id is `ensembl`-sourced the reference keeps the one of the
/// lowest `dbID`: the transcripts the files call `ensembl` carry the lower numbers, so the
/// same rows are set aside in a cache whose every transcript is of one other source.
#[test]
fn among_transcripts_of_no_ensembl_source_the_lowest_dbid_is_kept() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("lowest", threads);
        let cache = scratch.cache(|_| "BestRefSeq".to_string());
        scratch.outputs(ONE_BUFFER, &reference_rows_of_one_buffer(), &vep_rs_rows());
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(scratch.open_rows(), DISCORDANT_COLUMNS);
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(5, 3, 1, 4)), "{summary}");
    }
}

/// The rows the reference writes with one record a buffer. A buffer's transcripts are those
/// within reach of its lowest and highest record, so the buffer of the third record alone
/// holds no `ensembl`-sourced transcript of its two ids, and the buffer of the second never
/// loads the fourth id's transcript at the other locus: the reference writes those
/// transcripts' rows itself, and drops the third id's transcript as before.
fn reference_rows_of_one_record_a_buffer() -> Vec<String> {
    vec![
        transcript_row(
            "inside_both_copies",
            "2999291",
            "XM_900003.1",
            "non_coding_transcript_exon_variant",
        ),
        transcript_row("at_one_locus", "4919323", "XM_900004.1", "intron_variant"),
        transcript_row("at_one_locus", "4919323", "XM_900005.1", "intron_variant"),
        transcript_row(
            "past_the_kept_copies",
            "5996672",
            "XM_900001.1",
            "stop_gained",
        ),
        transcript_row(
            "past_the_kept_copies",
            "5996672",
            "XM_900002.1",
            "stop_gained",
        ),
        transcript_row(
            "at_the_other_locus",
            "9100000",
            "XM_900004.1",
            "non_coding_transcript_exon_variant",
        ),
    ]
}

#[test]
fn a_reference_that_read_one_record_a_buffer_dropped_the_transcript_beside_its_kept_one_alone() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("one-a-buffer", threads);
        let cache = scratch.cache(as_stated);
        scratch.outputs(
            ONE_RECORD_A_BUFFER,
            &reference_rows_of_one_record_a_buffer(),
            &vep_rs_rows(),
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(scratch.open_rows(), DISCORDANT_COLUMNS);
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(1, 1, 0, 1)), "{summary}");
    }
}

/// Under one record a buffer the replay finds the transcripts of the first, second and
/// fourth id kept at their records, so a reference without their rows there differs for
/// another reason: those rows and the intergenic row stay open.
#[test]
fn a_row_the_replay_does_not_account_for_stays_open() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("unaccounted", threads);
        let cache = scratch.cache(as_stated);
        scratch.outputs(
            ONE_RECORD_A_BUFFER,
            &reference_rows_of_one_buffer(),
            &vep_rs_rows(),
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            scratch.open_rows(),
            format!(
                "{DISCORDANT_COLUMNS}records.txt\tmissing_in_rust\tcontig:5996672\tG\t-\t-\t\
                 intergenic_variant\nrecords.txt\textra_in_rust\tcontig:4919323\tG\tXM_900004.1\t\
                 Transcript\tintron_variant\nrecords.txt\textra_in_rust\tcontig:5996672\tG\t\
                 XM_900001.1\tTranscript\tstop_gained\nrecords.txt\textra_in_rust\tcontig:5996672\t\
                 G\tXM_900002.1\tTranscript\tstop_gained\n"
            )
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(1, 1, 0, 1)), "{summary}");
    }
}

/// One transcript of the third id is dropped within reach of its record, so vep-rs's rows
/// of the id beyond the kept transcript's are one at most: a second one is a row too many,
/// and neither is told from the other, so both stay open.
#[test]
fn more_rows_of_an_id_than_transcripts_dropped_in_reach_stay_open() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("too-many", threads);
        let cache = scratch.cache(as_stated);
        let mut vep_rs = vep_rs_rows();
        vep_rs.push(transcript_row(
            "inside_both_copies",
            "2999291",
            "XM_900003.1",
            "synonymous_variant",
        ));
        scratch.outputs(ONE_BUFFER, &reference_rows_of_one_buffer(), &vep_rs);
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            format!(
                "{DISCORDANT_COLUMNS}records.txt\textra_in_rust\tcontig:2999291\tG\tXM_900003.1\t\
                 Transcript\tmissense_variant\nrecords.txt\textra_in_rust\tcontig:2999291\tG\t\
                 XM_900003.1\tTranscript\tsynonymous_variant\n"
            )
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(4, 2, 1, 3)), "{summary}");
    }
}

/// The reference wrote the kept transcript's row of the third id and vep-rs wrote another
/// row in its place: the row vep-rs lacks shows that its row is no row beside the kept
/// transcript's, and both stay open.
#[test]
fn a_row_of_the_id_that_vep_rs_lacks_keeps_its_row_of_the_id_open() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("lacks", threads);
        let cache = scratch.cache(as_stated);
        let vep_rs: Vec<String> = vep_rs_rows()
            .into_iter()
            .filter(|row| {
                !(row.contains("XM_900003.1") && row.contains("non_coding_transcript_exon_variant"))
            })
            .collect();
        scratch.outputs(ONE_BUFFER, &reference_rows_of_one_buffer(), &vep_rs);
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            format!(
                "{DISCORDANT_COLUMNS}records.txt\tmissing_in_rust\tcontig:2999291\tG\tXM_900003.1\t\
                 Transcript\tnon_coding_transcript_exon_variant\nrecords.txt\textra_in_rust\t\
                 contig:2999291\tG\tXM_900003.1\tTranscript\tmissense_variant\n"
            )
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(4, 2, 1, 3)), "{summary}");
    }
}

#[test]
fn a_run_the_replay_does_not_follow_leaves_the_rows_open_and_says_why() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("not-followed", threads);
        let cache = scratch.cache(as_stated);
        scratch.outputs(
            &format!("{ONE_BUFFER} --transcript_filter stable_id"),
            &reference_rows_of_one_buffer(),
            &vep_rs_rows(),
        );
        let output = scratch.compare(threads, &["--json-cache", &cache]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            stderr(&output),
            format!(
                "WARNING: [vep-outputs-compare] records.txt: {CLASS} is not applied: the \
                 reference ran under --transcript_filter, which the replay does not follow\n"
            )
        );
        assert_eq!(
            scratch.open_rows(),
            format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
        );
        let summary = scratch.read("report/summary.json");
        assert!(summary.contains(&class_block(0, 0, 0, 0)), "{summary}");
    }
}

/// Under 115.2 the cache gives the run's sequence names alone: its registry holds no such
/// class, the rows stay open and no report names the class.
#[test]
fn the_older_reference_applies_no_such_class() {
    let scratch = Scratch::new("older", None);
    let cache = scratch.cache(as_stated);
    scratch.outputs(ONE_BUFFER, &reference_rows_of_one_buffer(), &vep_rs_rows());
    let output = scratch.compare(
        None,
        &["--json-cache", &cache, "--reference-release", "115.2"],
    );
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    assert_eq!(stderr(&output), "");
    assert_eq!(
        scratch.open_rows(),
        format!("{DISCORDANT_COLUMNS}{CLASS_ROWS}")
    );
    for report in ["report/summary.json", "report/summary.md"] {
        assert!(!scratch.read(report).contains(CLASS));
    }
}
