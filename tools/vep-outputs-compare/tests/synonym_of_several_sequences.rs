// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over made-up output pairs beside a made-up Ensembl VEP cache directory
//! whose synonyms file gives one name to two of its sequences, and holds the class
//! `synonym_of_several_sequences_by_hash_order`: the reference's rows of a record under such
//! a name are set aside where vep-rs wrote no row for it, the run names that cache and the
//! reference output's header states a command line that names no `--synonyms`. They stay
//! open without the cache, under the 115.2 reference, where vep-rs wrote a row at the
//! record's Location, on a sequence the cache holds, and where the header names `--synonyms`
//! or states no command line. A cache whose file gives no name to two sequences leaves every
//! report as a run without the class writes it, and a cache whose file cannot be read ends
//! the run before a report file is written. A cache directory that cannot be listed ends the
//! run the same way. An output the class cannot read again, once the discordant rows are
//! written, ends the run on one line that names that output, and the discordance report
//! stands whole under its name.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const CLASS: &str = "synonym_of_several_sequences_by_hash_order";

const OUTPUT_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

/// The header line in which a reference run that names no synonyms file states its command.
const COMMAND_LINE: &str = "## VEP command-line: vep --cache --offline --no_stats\n";

/// One row of an output on a transcript.
fn transcript_row(location: &str, feature: &str, consequence: &str) -> String {
    format!("v\t{location}\tA\tG\t{feature}\tTranscript\t{consequence}\n")
}

/// The row an engine writes for a record no transcript is in reach of.
fn intergenic_row(location: &str) -> String {
    format!("v\t{location}\tA\t-\t-\t-\tintergenic_variant\n")
}

/// A directory of one test, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-shared-synonym-{}-{tag}",
            std::process::id()
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

    /// The reference's output under `command_line` and the column line, and vep-rs's under
    /// the column line.
    fn outputs_under(&self, command_line: &str, reference: &[String], rust: &[String]) {
        self.write(
            "perl/variants.txt",
            &[command_line, OUTPUT_COLUMNS, &reference.concat()].concat(),
        );
        self.write(
            "rust/variants.txt",
            &[OUTPUT_COLUMNS, &rust.concat()].concat(),
        );
    }

    /// [`Self::outputs_under`] the command line of a reference run that names no synonyms
    /// file.
    fn outputs(&self, reference: &[String], rust: &[String]) {
        self.outputs_under(COMMAND_LINE, reference, rust);
    }

    /// An Ensembl VEP cache directory that holds the sequences `I`, `II` and `III` and whose
    /// synonyms file is `synonyms`; none where `synonyms` is `None`.
    fn cache(&self, synonyms: Option<&str>) -> String {
        for sequence in ["I", "II", "III"] {
            fs::create_dir_all(self.0.join("cache").join(sequence)).unwrap();
        }
        self.write("cache/info.txt", "species\tmade_up\n");
        if let Some(synonyms) = synonyms {
            self.write("cache/chr_synonyms.txt", synonyms);
        }
        self.arg("cache")
    }

    /// Compares the two outputs with no contig filter under `extra`.
    fn compare(&self, extra: &[&str]) -> Output {
        let (perl, rust, report) = (self.arg("perl"), self.arg("rust"), self.arg("report"));
        let _ = fs::remove_dir_all(self.0.join("report"));
        Command::new(BIN)
            .current_dir(&self.0)
            .args(["--perl-dir", &perl, "--rust-dir", &rust])
            .args(["--report-dir", &report, "--no-canonical-contigs"])
            .args(extra)
            .output()
            .unwrap()
    }

    /// The data rows left open by the last comparison, each as Location and Feature.
    fn open_rows(&self) -> Vec<(String, String)> {
        let open = self.read("report/discordant_open.tsv");
        let cells = |line: &str| {
            let cells: Vec<&str> = line.split('\t').collect();
            (cells[2].to_string(), cells[4].to_string())
        };
        open.lines().skip(1).map(cells).collect()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn clean(output: &Output) {
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(String::from_utf8_lossy(&output.stderr), "");
}

/// The synonyms of a cache that gives `shared` to the held sequences `I` and `II`, `plain`
/// to `III` alone, and pairs the held sequence `II` with the two others.
const SHARED_BY_TWO: &str = "I\tshared\nII\tshared\nIII\tplain\nI\tI.1\nI\tII\nII\tIII\n";

/// The class's statistics as the summary writes them, one key a line.
fn class_statistics(summary: &str) -> String {
    let marker = format!("\"{CLASS}\": {{");
    let block = &summary[summary.find(&marker).expect("the class's statistics")..];
    block[..block.find('}').unwrap()]
        .split_whitespace()
        .collect()
}

/// Two records under the shared name that vep-rs skipped, one the reference read between
/// genes and one it read on two transcripts, beside a record both engines wrote alike: with
/// the cache named the three reference rows are set aside and the summary states the class
/// last with its counts; without the cache the rows stay open and no report names the
/// class.
#[test]
fn the_reference_rows_of_a_skipped_record_are_set_aside_where_the_cache_is_named() {
    let scratch = Scratch::new("skipped");
    let alike = transcript_row("III:100", "T0", "missense_variant");
    scratch.outputs(
        &[
            alike.clone(),
            intergenic_row("shared:399"),
            transcript_row("shared:800", "T1", "intron_variant"),
            transcript_row("shared:800", "T2", "upstream_gene_variant"),
        ],
        &[alike],
    );
    let cache = scratch.cache(Some(SHARED_BY_TWO));
    let three_open = [
        ("shared:399".to_string(), "-".to_string()),
        ("shared:800".to_string(), "T1".to_string()),
        ("shared:800".to_string(), "T2".to_string()),
    ];

    clean(&scratch.compare(&["--assembly", "MadeUp1"]));
    assert_eq!(scratch.open_rows(), three_open);
    let without = scratch.read("report/summary.json");
    assert!(!without.contains(CLASS), "{without}");
    assert!(!scratch.read("report/summary.md").contains(CLASS));

    clean(&scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]));
    assert_eq!(scratch.open_rows(), []);
    let with = scratch.read("report/summary.json");
    assert_eq!(
        class_statistics(&with),
        format!("\"{CLASS}\":{{\"count\":3,\"excluded\":3,\"records\":2,\"perl\":3,\"rust\":0")
    );
    assert!(with.contains("\"excluded_perl_tuples\": 3"), "{with}");
    assert!(with.contains("\"excluded_rust_tuples\": 0"), "{with}");
    let classes = &with[with.find("\"divergence_classes\"").unwrap()..];
    let last_class = classes.rfind("\"taxonomy_class\": \"").unwrap();
    assert!(
        classes[last_class..].starts_with(&format!("\"taxonomy_class\": \"{CLASS}\"")),
        "the class is stated last: {}",
        &classes[last_class..]
    );
    assert_eq!(scratch.read("report/discordant.tsv").lines().count(), 4);

    // Under the 115.2 reference the class is not applied.
    let older = ["--reference-release", "115.2", "--ensembl-cache", &cache];
    clean(&scratch.compare(&[&["--assembly", "MadeUp1"][..], &older].concat()));
    assert_eq!(scratch.open_rows(), three_open);
    assert!(!scratch.read("report/summary.json").contains(CLASS));

    // A path that is no directory is a usage error.
    let absent = scratch.arg("absent");
    let refused = scratch.compare(&["--ensembl-cache", &absent]);
    assert_eq!(refused.status.code(), Some(2));
    let message = String::from_utf8_lossy(&refused.stderr).into_owned();
    assert!(
        message.contains(&format!("--ensembl-cache {absent} is not a directory")),
        "{message}"
    );
}

/// A record under the shared name that vep-rs wrote a row for is one it did not skip: its
/// rows stay open, whether vep-rs's row is another transcript's or one the reference wrote
/// too. The rows of a name the file gives to one sequence stay open as well, and so do the
/// rows of `II`, a sequence the cache holds that the file pairs with two others: rows
/// vep-rs lost on a held sequence are no rows of the class. The class counts nothing.
#[test]
fn a_record_vep_rs_wrote_a_row_for_stays_open() {
    let scratch = Scratch::new("written");
    let both = transcript_row("shared:500", "T5", "intron_variant");
    scratch.outputs(
        &[
            transcript_row("shared:399", "T1", "intron_variant"),
            both.clone(),
            transcript_row("shared:500", "T6", "intron_variant"),
            transcript_row("plain:10", "T7", "intron_variant"),
            transcript_row("II:10", "T8", "intron_variant"),
        ],
        &[transcript_row("shared:399", "T2", "intron_variant"), both],
    );
    let cache = scratch.cache(Some(SHARED_BY_TWO));
    clean(&scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]));
    let mut open = scratch.open_rows();
    open.sort();
    let named = |location: &str, feature: &str| (location.to_string(), feature.to_string());
    assert_eq!(
        open,
        [
            named("II:10", "T8"),
            named("plain:10", "T7"),
            named("shared:399", "T1"),
            named("shared:399", "T2"),
            named("shared:500", "T6"),
        ]
    );
    assert_eq!(
        class_statistics(&scratch.read("report/summary.json")),
        format!("\"{CLASS}\":{{\"count\":0,\"excluded\":0,\"records\":0,\"perl\":0,\"rust\":0")
    );
}

/// A cache whose synonyms file gives no name to two held sequences, and one without the
/// file, leave the summary as a run that names no cache writes it, the class unnamed.
#[test]
fn a_cache_that_shares_no_name_changes_no_report() {
    let scratch = Scratch::new("unshared");
    scratch.outputs(
        &[
            transcript_row("III:100", "T0", "missense_variant"),
            transcript_row("plain:10", "T7", "intron_variant"),
        ],
        &[transcript_row("III:100", "T0", "missense_variant")],
    );
    let stable = |summary: String| -> String {
        let kept = |line: &&str| !line.contains("\"generated_at_utc\"");
        summary.lines().filter(kept).collect::<Vec<_>>().join("\n")
    };
    clean(&scratch.compare(&["--assembly", "MadeUp1"]));
    let bare = stable(scratch.read("report/summary.json"));
    for synonyms in [
        Some("I\tplain\nII\tI.1\nunheld\tboth\nelsewhere\tboth\n"),
        None,
    ] {
        let _ = fs::remove_dir_all(scratch.0.join("cache"));
        let cache = scratch.cache(synonyms);
        clean(&scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]));
        let named = stable(scratch.read("report/summary.json"));
        assert_eq!(named, bare, "{synonyms:?}");
        assert!(!named.contains(CLASS));
    }
}

/// A cache whose synonyms file cannot be read, a directory of that name here, ends the run
/// with exit code 1 on an error that names the file, and the report directory holds no file:
/// the cache is read before the discordance file is created.
#[test]
fn a_cache_whose_synonyms_cannot_be_read_ends_the_run_before_any_report() {
    let scratch = Scratch::new("unreadable");
    let alike = [transcript_row("III:100", "T0", "missense_variant")];
    scratch.outputs(&alike, &alike);
    let cache = scratch.cache(None);
    fs::create_dir(scratch.0.join("cache/chr_synonyms.txt")).unwrap();
    let stopped = scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]);
    assert_eq!(stopped.status.code(), Some(1));
    let message = String::from_utf8_lossy(&stopped.stderr).into_owned();
    let synonyms = scratch.arg("cache/chr_synonyms.txt");
    let opening = format!("ERROR: [vep-outputs-compare] --ensembl-cache: {synonyms}: ");
    assert!(message.starts_with(&opening), "{message}");
    let report = fs::read_dir(scratch.0.join("report")).unwrap();
    let left: Vec<_> = report.map(|entry| entry.unwrap().file_name()).collect();
    assert!(left.is_empty(), "{left:?}");
}

/// Only a reference run that read the cache's own synonyms file is of the class. Where its
/// header states a command line that names `--synonyms`, or no command line, the rows of a
/// record under the shared name stay open and the class counts nothing.
#[test]
fn a_reference_run_given_a_synonyms_file_keeps_its_rows_open() {
    let scratch = Scratch::new("given-synonyms");
    let cache = scratch.cache(Some(SHARED_BY_TWO));
    let alike = transcript_row("III:100", "T0", "missense_variant");
    let reference = [alike.clone(), intergenic_row("shared:399")];
    let open = [("shared:399".to_string(), "-".to_string())];
    for command_line in [
        "## VEP command-line: vep --cache --offline --synonyms [PATH]/names.txt\n",
        "",
    ] {
        scratch.outputs_under(command_line, &reference, std::slice::from_ref(&alike));
        clean(&scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]));
        assert_eq!(scratch.open_rows(), open, "{command_line:?}");
        assert_eq!(
            class_statistics(&scratch.read("report/summary.json")),
            format!("\"{CLASS}\":{{\"count\":0,\"excluded\":0,\"records\":0,\"perl\":0,\"rust\":0"),
            "{command_line:?}"
        );
    }
    scratch.outputs(&reference, &[alike]);
    clean(&scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]));
    assert_eq!(scratch.open_rows(), []);
}

/// A cache directory whose entries cannot be listed ends the run as a synonyms file that
/// cannot be read does: exit code 1 on an error that names the directory, and the report
/// directory holds no file. A user who may list any directory cannot meet this failure.
#[test]
fn a_cache_directory_that_cannot_be_listed_ends_the_run_before_any_report() {
    use std::os::unix::fs::PermissionsExt;
    let scratch = Scratch::new("unlisted");
    let alike = [transcript_row("III:100", "T0", "missense_variant")];
    scratch.outputs(&alike, &alike);
    let cache = scratch.cache(Some(SHARED_BY_TWO));
    fs::set_permissions(&cache, fs::Permissions::from_mode(0o000)).unwrap();
    let stopped = fs::read_dir(&cache)
        .is_err()
        .then(|| scratch.compare(&["--assembly", "MadeUp1", "--ensembl-cache", &cache]));
    fs::set_permissions(&cache, fs::Permissions::from_mode(0o755)).unwrap();
    let Some(stopped) = stopped else {
        return;
    };
    assert_eq!(stopped.status.code(), Some(1));
    let message = String::from_utf8_lossy(&stopped.stderr).into_owned();
    let opening = format!("ERROR: [vep-outputs-compare] --ensembl-cache: {cache}: ");
    assert!(message.starts_with(&opening), "{message}");
    assert_eq!(message.lines().count(), 1, "{message}");
    let report = fs::read_dir(scratch.0.join("report")).unwrap();
    let left: Vec<_> = report.map(|entry| entry.unwrap().file_name()).collect();
    assert!(left.is_empty(), "{left:?}");
}

/// The environment variable that bounds the worker threads of a run.
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// Holds a run to its end where `gone`, an output of the first of two paired files, is not
/// there when the class reads the outputs again. The first file holds a record under the
/// shared name that vep-rs skipped, and the input VCF of the second is a named pipe, at whose
/// open the run waits: by then it has compared the first file, so `gone` is removed there
/// and the pipe closed. On one thread and on the host's default, the run ends with exit code
/// 1 on one line that names `gone`, and the report directory holds the discordance report
/// alone, whole: the rows of every paired file were written before the pass that reads the
/// outputs again.
fn the_run_ends_on_the_name_of_an_output_gone(gone: &str) {
    for threads in [Some("1"), None] {
        let scratch = Scratch::new(&format!(
            "gone-{}-t{}",
            gone.replace('/', "-"),
            threads.unwrap_or("default")
        ));
        let cache = scratch.cache(Some(SHARED_BY_TWO));
        let alike = transcript_row("III:100", "T0", "missense_variant");
        let skipped = intergenic_row("shared:399");
        let second = transcript_row("III:200", "T9", "intron_variant");
        scratch.write(
            "perl/a.txt",
            &[COMMAND_LINE, OUTPUT_COLUMNS, &alike, &skipped].concat(),
        );
        scratch.write("rust/a.txt", &[OUTPUT_COLUMNS, &alike].concat());
        scratch.write(
            "perl/b.txt",
            &[COMMAND_LINE, OUTPUT_COLUMNS, &second].concat(),
        );
        scratch.write("rust/b.txt", &[OUTPUT_COLUMNS, &second].concat());
        scratch.write("input/a.vcf", "");
        let pipe = scratch.0.join("input/b.vcf");
        let made = Command::new("mkfifo").arg(&pipe).status().unwrap();
        assert!(made.success(), "mkfifo {}", pipe.display());
        let (perl, rust, report) = (
            scratch.arg("perl"),
            scratch.arg("rust"),
            scratch.arg("report"),
        );
        let input = scratch.arg("input");
        let mut command = Command::new(BIN);
        command
            .current_dir(&scratch.0)
            .args(["--perl-dir", &perl, "--rust-dir", &rust])
            .args(["--report-dir", &report, "--no-canonical-contigs"])
            .args(["--assembly", "MadeUp1", "--ensembl-cache", &cache])
            .args(["--input-dir", &input]);
        match threads {
            Some(bound) => command.env(THREADS_ENV, bound),
            None => command.env_remove(THREADS_ENV),
        };
        let run = command
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        let output = scratch.0.join(gone);
        let (removed, was_removed) = std::sync::mpsc::channel();
        let (at_the_pipe, output_to_remove) = (pipe.clone(), output.clone());
        std::thread::spawn(move || {
            // The open for writing returns once the run has opened the pipe for reading.
            let writer = fs::OpenOptions::new()
                .write(true)
                .open(&at_the_pipe)
                .unwrap();
            fs::remove_file(&output_to_remove).unwrap();
            // A later open of the input VCF finds an empty file and waits for no writer.
            fs::remove_file(&at_the_pipe).unwrap();
            fs::write(&at_the_pipe, "").unwrap();
            removed.send(()).unwrap();
            drop(writer);
        });
        let ended = run.wait_with_output().unwrap();
        assert!(
            was_removed.try_recv().is_ok(),
            "{threads:?}: the run ended without opening the input VCF of the second file"
        );
        assert_eq!(ended.status.code(), Some(1), "{threads:?}");
        let missing = fs::File::open(&output).unwrap_err();
        assert_eq!(
            String::from_utf8_lossy(&ended.stderr),
            format!(
                "ERROR: [vep-outputs-compare] {}: {missing}\n",
                output.display()
            ),
            "{threads:?}"
        );
        let report = fs::read_dir(scratch.0.join("report")).unwrap();
        let left: Vec<String> = report
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        assert_eq!(left, ["discordant.tsv"], "{threads:?}");
        assert_eq!(
            scratch.read("report/discordant.tsv"),
            "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n\
             a.txt\tmissing_in_rust\tshared:399\tA\t-\t-\tintergenic_variant\n",
            "{threads:?}"
        );
    }
}

/// The class reads the reference's output again for the command line of its header: one
/// that is gone by then is named in the failure.
#[test]
fn a_reference_output_gone_before_the_class_reads_its_header_is_named_in_the_failure() {
    the_run_ends_on_the_name_of_an_output_gone("perl/a.txt");
}

/// The class reads vep-rs's output again for the rows at the Locations it judges: one that
/// is gone by then is named in the failure.
#[test]
fn a_vep_rs_output_gone_before_the_class_reads_its_rows_is_named_in_the_failure() {
    the_run_ends_on_the_name_of_an_output_gone("rust/a.txt");
}
