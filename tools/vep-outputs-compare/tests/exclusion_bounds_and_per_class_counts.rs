// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over output pairs and saved reports written for each test, on one thread
//! and on the host's default, and holds what it writes and how it ends: the open rows of a
//! compare run beside its discordant rows, the contig filter a `--canonical-contigs` value or
//! a cache's `contig_order` names for an assembly that is not a human one, the breakend mate
//! told from a synonym by the names such a run gives, the per-term totals cache of the
//! per-class mode, and the usage errors of an input path that is not there and of a contig
//! value that names no contig.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::Duration;

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is held under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 2] = [Some("1"), None];

const OUTPUT_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

const DISCORDANT_COLUMNS: &str =
    "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n";

/// One row of an output on a transcript.
fn output_row(location: &str, allele: &str, feature: &str, csq: &str) -> String {
    format!("v\t{location}\t{allele}\tG\t{feature}\tTranscript\t{csq}\n")
}

/// A directory of one test under one thread bound, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, threads: Option<&str>) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-cases-{}-{tag}-t{}",
            std::process::id(),
            threads.unwrap_or("default")
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(fs::canonicalize(&path).unwrap())
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    /// The path of `name` as an argument of the command line.
    fn arg(&self, name: &str) -> String {
        self.path(name).to_string_lossy().into_owned()
    }

    fn write(&self, name: &str, content: &str) {
        let path = self.path(name);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, content).unwrap();
    }

    fn read(&self, name: &str) -> String {
        let path = self.path(name);
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
    }

    /// Runs the binary from this directory, so that a relative path of `args` lies in it.
    fn run(&self, threads: Option<&str>, args: &[&str]) -> Output {
        let mut command = Command::new(BIN);
        command.current_dir(&self.0).args(args);
        match threads {
            Some(bound) => command.env(THREADS_ENV, bound),
            None => command.env_remove(THREADS_ENV),
        };
        command.output().unwrap()
    }

    /// The names of the entries of the directory `name`, sorted.
    fn entries(&self, name: &str) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(self.path(name))
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort_unstable();
        names
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

/// The stdout of a run that ended with exit code 0 and an empty stderr.
fn stdout_of_a_clean_run(output: &Output) -> String {
    assert_eq!(output.status.code(), Some(0), "{}", stderr(output));
    assert_eq!(stderr(output), "");
    stdout(output)
}

/// One compare run under the 116.2 reference of the outputs under `perl` and `rust` of the
/// scratch into its `report`.
fn compare(scratch: &Scratch, threads: Option<&str>, extra: &[&str]) -> Output {
    let (perl, rust, report) = (
        scratch.arg("perl"),
        scratch.arg("rust"),
        scratch.arg("report"),
    );
    let paths = [
        "--perl-dir",
        &perl,
        "--rust-dir",
        &rust,
        "--report-dir",
        &report,
    ];
    let release = ["--reference-release", "116.2"];
    scratch.run(threads, &[&paths[..], &release, extra].concat())
}

/// The lines a compare run under the 116.2 reference prints before its scores: the release
/// and the four files it wrote, the discordant rows under `discordant` and its open rows.
fn wrote_lines(scratch: &Scratch, discordant: &str) -> String {
    let report = scratch.arg("report");
    format!(
        "reference release: 116.2\n\
         wrote summary: {report}/summary.json\n\
         wrote markdown: {report}/summary.md\n\
         wrote discordance: {report}/{discordant}.tsv\n\
         wrote open discordance: {report}/{discordant}_open.tsv\n"
    )
}

/// The first three keys of the `aggregate` block of a summary.
fn aggregate_counts(perl: u64, rust: u64, shared: u64) -> String {
    format!(
        "  \"aggregate\": {{\n    \"perl_tuple_count\": {perl},\n    \"rust_tuple_count\": {rust},\n    \
         \"intersection_count\": {shared},\n"
    )
}

/// The whole number under `key` in the `adjusted` block of a summary.
fn adjusted_count(summary: &str, key: &str) -> usize {
    let block = &summary[summary.find("\"adjusted\": {").expect("the adjusted block")..];
    let marker = format!("\"{key}\": ");
    let number = &block[block.find(&marker).expect("the key") + marker.len()..];
    let end = number.find(|c: char| !c.is_ascii_digit()).unwrap();
    number[..end].parse().unwrap()
}

/// A reference output of four tuples and a vep-rs output of three: one tuple shared, one pair
/// the adjusted F1 sets aside (`splice_region_variant` added), one pair that is counted and
/// stays (the polypyrimidine term added) and one tuple the reference alone wrote.
fn outputs_with_one_excluded_pair(scratch: &Scratch) {
    let shared = output_row("21:1", "A", "ENST0", "intron_variant");
    let perl = [
        OUTPUT_COLUMNS,
        &shared,
        &output_row("21:100", "A", "ENST1", "missense_variant"),
        &output_row("21:300", "T", "ENST3", "intron_variant"),
        &output_row("21:400", "G", "ENST4", "downstream_gene_variant"),
    ];
    let rust = [
        OUTPUT_COLUMNS,
        &shared,
        &output_row(
            "21:100",
            "A",
            "ENST1",
            "missense_variant,splice_region_variant",
        ),
        &output_row(
            "21:300",
            "T",
            "ENST3",
            "intron_variant,splice_polypyrimidine_tract_variant",
        ),
    ];
    scratch.write("perl/output.txt", &perl.concat());
    scratch.write("rust/output.txt", &rust.concat());
    fs::create_dir(scratch.path("report")).unwrap();
}

const DISCORDANT_ROWS: &str = "\
file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set
output.txt\tmissing_in_rust\t21:100\tA\tENST1\tTranscript\tmissense_variant
output.txt\tmissing_in_rust\t21:300\tT\tENST3\tTranscript\tintron_variant
output.txt\tmissing_in_rust\t21:400\tG\tENST4\tTranscript\tdownstream_gene_variant
output.txt\textra_in_rust\t21:100\tA\tENST1\tTranscript\tmissense_variant,splice_region_variant
output.txt\textra_in_rust\t21:300\tT\tENST3\tTranscript\tintron_variant,splice_polypyrimidine_tract_variant
";

const OPEN_ROWS: &str = "\
file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set
output.txt\tmissing_in_rust\t21:300\tT\tENST3\tTranscript\tintron_variant
output.txt\tmissing_in_rust\t21:400\tG\tENST4\tTranscript\tdownstream_gene_variant
output.txt\textra_in_rust\t21:300\tT\tENST3\tTranscript\tintron_variant,splice_polypyrimidine_tract_variant
";

const SCORES_WITH_ONE_EXCLUDED_PAIR: &str = "\
aggregate: precision=0.333333 recall=0.250000 f1=0.285714
adjusted: precision=0.500000 recall=0.333333 f1=0.400000 (excluded 1P + 1R = 2 intended-divergence tuples)
aggregate (context-collapsed): precision=0.333333 recall=0.250000 f1=0.285714
";

#[test]
fn a_compare_run_writes_the_open_rows_beside_the_discordant_rows_and_they_number_the_adjusted_total(
) {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("open-rows", threads);
        outputs_with_one_excluded_pair(&scratch);
        let output = compare(&scratch, threads, &["--assembly", "GRCh37"]);
        assert_eq!(
            stdout_of_a_clean_run(&output),
            wrote_lines(&scratch, "discordant") + SCORES_WITH_ONE_EXCLUDED_PAIR
        );
        let summary = scratch.read("report/summary.json");
        assert!(
            summary.contains(
                "    \"adjusted\": {\n      \"perl_tuple_count\": 3,\n      \"rust_tuple_count\": 2,\n      \
                 \"intersection_count\": 1,\n      \"precision\": 0.5,\n      \
                 \"recall\": 0.3333333333333333,\n      \"f1\": 0.4,\n      \
                 \"excluded_perl_tuples\": 1,\n      \"excluded_rust_tuples\": 1,\n"
            ),
            "{summary}"
        );
        let adjusted = |key| adjusted_count(&summary, key);
        let shared = adjusted("intersection_count");
        let raw = scratch.read("report/discordant.tsv");
        let open = scratch.read("report/discordant_open.tsv");
        assert_eq!(raw, DISCORDANT_ROWS);
        assert_eq!(open, OPEN_ROWS);
        assert_eq!(raw.lines().count() - 1, 5);
        assert_eq!(
            open.lines().count() - 1,
            (adjusted("perl_tuple_count") - shared) + (adjusted("rust_tuple_count") - shared)
        );
    }
}

#[test]
fn the_open_file_takes_its_name_from_the_file_of_the_discordant_rows() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("open-name", threads);
        outputs_with_one_excluded_pair(&scratch);
        let named = ["--assembly", "GRCh37", "--discordant-tsv", "custom.tsv"];
        let output = compare(&scratch, threads, &named);
        assert_eq!(
            stdout_of_a_clean_run(&output),
            wrote_lines(&scratch, "custom") + SCORES_WITH_ONE_EXCLUDED_PAIR
        );
        assert_eq!(
            scratch.entries("report"),
            [
                "custom.tsv",
                "custom_open.tsv",
                "summary.json",
                "summary.md"
            ],
            "no file under the default names"
        );
        assert_eq!(scratch.read("report/custom.tsv"), DISCORDANT_ROWS);
        assert_eq!(scratch.read("report/custom_open.tsv"), OPEN_ROWS);
    }
}

/// The contigs a fruit fly cache lists in its `contig_order`.
const FLY_CONTIGS: [&str; 8] = [
    "2L",
    "2R",
    "3L",
    "3R",
    "4",
    "X",
    "Y",
    "mitochondrion_genome",
];

/// The human set in the order a summary lists it.
const HUMAN_CONTIGS: [&str; 26] = [
    "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "3",
    "4", "5", "6", "7", "8", "9", "M", "MT", "X", "Y",
];

/// The last three keys of a summary: the assembly, whether the contig filter ran and its set.
fn summary_tail(assembly: &str, contigs: Option<&[&str]>) -> String {
    let set = match contigs {
        None => "null".to_string(),
        Some(contigs) => {
            let names: Vec<String> = contigs.iter().map(|c| format!("    \"{c}\"")).collect();
            format!("[\n{}\n  ]", names.join(",\n"))
        }
    };
    format!(
        "  \"assembly\": \"{assembly}\",\n  \"canonical_contigs_filter_active\": {},\n  \
         \"canonical_contigs_set\": {set}\n}}",
        contigs.is_some()
    )
}

/// The warning of a bare flag under `assembly`, which names no set, in a run that names no
/// cache.
fn bare_flag_warning(assembly: &str) -> String {
    format!("WARNING: --canonical-contigs (default ON) requested but --assembly is {assembly}, which names no canonical set. Defensive output-side filter is OFF for this run. The upstream input-prep canonical-contigs filter (vep-prepare-benchmark-vcfs) still applies. Pass --assembly GRCh37 or --assembly GRCh38, name the set with --canonical-contigs <file or list>, or name a converted cache with --json-cache <directory> (the set is its info.json contig_order), to enable the defensive layer, or --no-canonical-contigs to silence this warning.\n")
}

const ALL_SHARED: &str = "\
aggregate: precision=1.000000 recall=1.000000 f1=1.000000
aggregate (context-collapsed): precision=1.000000 recall=1.000000 f1=1.000000
";

const THREE_OF_FOUR_SHARED: &str = "\
aggregate: precision=0.750000 recall=0.750000 f1=0.750000
aggregate (context-collapsed): precision=0.750000 recall=0.750000 f1=0.750000
";

#[test]
fn the_contig_filter_of_an_assembly_that_is_not_human_is_the_set_the_flag_names() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("fly-contigs", threads);
        let shared = [
            output_row("2L:100", "A", "FBtr0", "missense_variant"),
            output_row("X:200", "A", "FBtr1", "intron_variant"),
            output_row(
                "mitochondrion_genome:300",
                "A",
                "FBtr2",
                "synonymous_variant",
            ),
        ]
        .concat();
        let outside_perl = output_row("rDNA:400", "A", "FBtr3", "intron_variant");
        let outside_rust = output_row(
            "Unmapped_Scaffold_8_D1580_D1567:500",
            "A",
            "FBtr4",
            "intron_variant",
        );
        scratch.write(
            "perl/output.txt",
            &[OUTPUT_COLUMNS, &shared, &outside_perl].concat(),
        );
        scratch.write(
            "rust/output.txt",
            &[OUTPUT_COLUMNS, &shared, &outside_rust].concat(),
        );
        fs::create_dir(scratch.path("report")).unwrap();
        let run = |extra: &[&str]| {
            let printed = stdout_of_a_clean_run(&compare(&scratch, threads, extra));
            (scratch.read("report/summary.json"), printed)
        };
        let all_shared = wrote_lines(&scratch, "discordant") + ALL_SHARED;

        // The set as a comma-separated list: the two rows outside it are not keyed, and
        // three tuples a side stay, all shared.
        let list = FLY_CONTIGS.join(",");
        let (summary, printed) = run(&["--assembly", "BDGP6.54", "--canonical-contigs", &list]);
        assert!(
            summary.ends_with(&summary_tail("BDGP6.54", Some(&FLY_CONTIGS))),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(3, 3, 3)), "{summary}");
        assert_eq!(printed, all_shared, "no warning");

        // The same set from a file, one name a line, a comment line ignored.
        let listing = format!("# the cache's karyotype\n{}\n", FLY_CONTIGS.join("\n"));
        scratch.write("contigs.txt", &listing);
        let file = scratch.arg("contigs.txt");
        let (summary, printed) = run(&["--assembly", "BDGP6.54", "--canonical-contigs", &file]);
        assert!(
            summary.ends_with(&summary_tail("BDGP6.54", Some(&FLY_CONTIGS))),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(3, 3, 3)), "{summary}");
        assert_eq!(printed, all_shared);

        // The bare flag: no set follows from the assembly, the filter is off with a warning
        // that names the assembly and the remedy, and the rows outside the set are counted.
        let (summary, printed) = run(&["--assembly", "BDGP6.54"]);
        assert!(
            summary.ends_with(&summary_tail("BDGP6.54", None)),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(4, 4, 3)), "{summary}");
        assert_eq!(
            printed,
            [
                &bare_flag_warning("BDGP6.54"),
                &wrote_lines(&scratch, "discordant"),
                THREE_OF_FOUR_SHARED
            ]
            .concat()
        );

        // The bare flag under a human assembly keeps the human set.
        let (summary, printed) = run(&["--assembly", "GRCh38"]);
        assert!(
            summary.ends_with(&summary_tail("GRCh38", Some(&HUMAN_CONTIGS))),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(1, 1, 1)), "{summary}");
        assert_eq!(printed, all_shared);
    }
}

/// The sequences the made-up cache of these tests lists in its `contig_order`.
const CACHE_CONTIGS: [&str; 3] = ["1", "II", "Mito"];

/// Writes the `info.json` of a converted cache under `cache` of the scratch, its
/// `contig_order` the JSON text `contig_order`, and returns the cache directory as an
/// argument of the command line.
fn cache_listing(scratch: &Scratch, contig_order: &str) -> String {
    scratch.write(
        "cache/info.json",
        &format!(
            "{{\n  \"assembly\": \"MadeUp1\",\n  \"contig_order\": {contig_order},\n  \
             \"species\": \"made_up\"\n}}\n"
        ),
    );
    scratch.arg("cache")
}

#[test]
fn the_contig_filter_of_an_assembly_that_is_not_human_is_the_contig_order_of_the_named_cache() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("cache-contigs", threads);
        let shared = [
            output_row("1:100", "A", "T0", "missense_variant"),
            output_row("II:200", "A", "T1", "intron_variant"),
            output_row("Mito:300", "A", "T2", "synonymous_variant"),
        ]
        .concat();
        let outside_perl = output_row("scaffold_9:400", "A", "T3", "intron_variant");
        let outside_rust = output_row("scaffold_12:500", "A", "T4", "intron_variant");
        scratch.write(
            "perl/output.txt",
            &[OUTPUT_COLUMNS, &shared, &outside_perl].concat(),
        );
        scratch.write(
            "rust/output.txt",
            &[OUTPUT_COLUMNS, &shared, &outside_rust].concat(),
        );
        fs::create_dir(scratch.path("report")).unwrap();
        let cache = cache_listing(&scratch, r#"["1", "II", "Mito"]"#);
        let run = |extra: &[&str]| {
            let printed = stdout_of_a_clean_run(&compare(&scratch, threads, extra));
            (scratch.read("report/summary.json"), printed)
        };
        let all_shared = wrote_lines(&scratch, "discordant") + ALL_SHARED;
        let three_of_four = wrote_lines(&scratch, "discordant") + THREE_OF_FOUR_SHARED;

        // The bare flag beside a cache: the rows of every sequence the cache lists are
        // keyed and the row of an unlisted one on each side is not.
        let (summary, printed) = run(&["--assembly", "MadeUp1", "--json-cache", &cache]);
        assert!(
            summary.ends_with(&summary_tail("MadeUp1", Some(&CACHE_CONTIGS))),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(3, 3, 3)), "{summary}");
        assert_eq!(printed, all_shared, "no warning");

        // Neither a set nor a cache: every row is compared, and the warning says so.
        let (summary, printed) = run(&["--assembly", "MadeUp1"]);
        assert!(
            summary.ends_with(&summary_tail("MadeUp1", None)),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(4, 4, 3)), "{summary}");
        assert_eq!(printed, bare_flag_warning("MadeUp1") + &three_of_four);

        // The filter turned off beside a cache stays off, without a warning.
        let off = ["--no-canonical-contigs", "--json-cache", &cache];
        let (summary, printed) = run(&[&["--assembly", "MadeUp1"][..], &off].concat());
        assert!(
            summary.ends_with(&summary_tail("MadeUp1", None)),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(4, 4, 3)), "{summary}");
        assert_eq!(printed, three_of_four);

        // A set the flag names is the filter's, whatever the cache lists.
        let named = ["--canonical-contigs", "1,II", "--json-cache", &cache];
        let (summary, printed) = run(&[&["--assembly", "MadeUp1"][..], &named].concat());
        assert!(
            summary.ends_with(&summary_tail("MadeUp1", Some(&["1", "II"]))),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(2, 2, 2)), "{summary}");
        assert_eq!(printed, all_shared);

        // A human assembly keeps the human set, whatever the cache lists.
        let (summary, printed) = run(&["--assembly", "GRCh38", "--json-cache", &cache]);
        assert!(
            summary.ends_with(&summary_tail("GRCh38", Some(&HUMAN_CONTIGS))),
            "{summary}"
        );
        assert!(summary.contains(&aggregate_counts(1, 1, 1)), "{summary}");
        assert_eq!(printed, all_shared);
    }
}

#[test]
fn a_named_cache_that_lists_no_sequence_names_ends_the_run() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("cache-without-names", threads);
        let row = output_row("1:100", "A", "T0", "missense_variant");
        scratch.write("perl/output.txt", &[OUTPUT_COLUMNS, &row].concat());
        scratch.write("rust/output.txt", &[OUTPUT_COLUMNS, &row].concat());
        fs::create_dir(scratch.path("report")).unwrap();
        let cases = [
            ("null", "contig_order lists no name"),
            ("[]", "contig_order lists no name"),
            (r#""1,II""#, "contig_order is not a list of names"),
            ("[1, 2]", "contig_order is not a list of names"),
        ];
        for (contig_order, reason) in cases {
            let cache = cache_listing(&scratch, contig_order);
            let refused = compare(
                &scratch,
                threads,
                &["--assembly", "MadeUp1", "--json-cache", &cache],
            );
            assert_eq!(refused.status.code(), Some(1), "{contig_order}");
            assert_eq!(
                stderr(&refused),
                format!("ERROR: [vep-outputs-compare] --json-cache: {cache}/info.json: {reason}\n")
            );
            assert_eq!(stdout(&refused), "");
            assert!(scratch.entries("report").is_empty(), "{contig_order}");

            // Under a human assembly the tool's own names stand and the cache is not read.
            let human = ["--assembly", "GRCh38", "--json-cache", &cache];
            stdout_of_a_clean_run(&compare(&scratch, threads, &human));
            fs::remove_dir_all(scratch.path("report")).unwrap();
            fs::create_dir(scratch.path("report")).unwrap();
        }

        // A file that is not JSON, and a cache directory without the file.
        scratch.write("cache/info.json", "contig_order: 1\n");
        let cache = scratch.arg("cache");
        let named = ["--assembly", "MadeUp1", "--json-cache", &cache];
        let unread = compare(&scratch, threads, &named);
        assert_eq!(unread.status.code(), Some(1));
        let prefix = format!("ERROR: [vep-outputs-compare] --json-cache: {cache}/info.json: ");
        assert!(stderr(&unread).starts_with(&prefix), "{}", stderr(&unread));
        fs::remove_file(scratch.path("cache/info.json")).unwrap();
        let absent = compare(&scratch, threads, &named);
        assert_eq!(absent.status.code(), Some(1));
        assert!(stderr(&absent).starts_with(&prefix), "{}", stderr(&absent));

        // A cache path that is no directory is a usage error.
        let missing = scratch.arg("no-cache");
        let usage = compare(&scratch, threads, &["--json-cache", &missing]);
        assert_eq!(usage.status.code(), Some(2));
        assert_eq!(
            stderr(&usage),
            format!(
                "usage: vep-outputs-compare [options]\nvep-outputs-compare: error: \
                 --json-cache {missing} is not a directory\n"
            )
        );
    }
}

/// The statistics of one class in the `excluded_categories` of a summary, each with its name,
/// in the order the summary writes them.
fn class_statistics(summary: &str, class: &str) -> Vec<(String, u64)> {
    let categories = summary
        .find("\"excluded_categories\"")
        .expect("excluded categories");
    let marker = format!("\"{class}\": {{");
    let block = &summary[categories..];
    let block = &block[block.find(&marker).expect("the class") + marker.len()..];
    let block = &block[..block.find('}').expect("the end of the class")];
    block
        .split(',')
        .map(|member| {
            let (name, value) = member.split_once(':').expect("a statistic");
            (
                name.trim().trim_matches('"').to_string(),
                value.trim().parse().expect("a count"),
            )
        })
        .collect()
}

/// The statistics of a one-sided class as a summary writes them, with the rows that reached
/// a test the run could not make where there are any.
fn one_sided(count: u64, records: u64, not_evaluated: Option<u64>) -> Vec<(String, u64)> {
    let mut statistics = vec![
        ("count", count),
        ("excluded", 0),
        ("records", records),
        ("perl", 0),
        ("rust", 0),
    ];
    statistics.extend(not_evaluated.map(|rows| ("not_evaluated", rows)));
    statistics
        .into_iter()
        .map(|(name, rows)| (name.to_string(), rows))
        .collect()
}

const NOT_EVALUATED_LINE: &str = "- Not evaluated (the run gives no sequence names of its cache): \
     bnd_synonym_mate_rows_lost: 1 row\n";

/// A breakend whose own-end row both engines wrote and whose mate row vep-rs alone wrote is
/// counted under `bnd_synonym_mate_rows_lost` where the mate is written under a name the
/// cache does not file its sequence under. The names are the human ones under a human
/// assembly alone: under any other the run's cache or named set gives them, a mate under the
/// cache's own name is no synonym, and a run that gives none makes no test and says so.
#[test]
fn a_breakend_mate_is_told_from_a_synonym_by_the_sequence_names_the_run_gives() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("breakend-mate-names", threads);
        let cache = cache_listing(&scratch, r#"["1", "II", "Mito"]"#);
        let breakend = |allele: &str, feature: &str| {
            format!("r1\tII:500\t{allele}\tG\t{feature}\tTranscript\tfeature_truncation\n")
        };
        let shared = output_row("1:100", "A", "T0", "missense_variant");
        let own_end = breakend("N.", "T_OWN");
        scratch.write(
            "perl/output.txt",
            &[OUTPUT_COLUMNS, &shared, &own_end].concat(),
        );
        let class = |mate: &str, names: &[&str]| {
            let mate_row = breakend(mate, "T_MATE");
            scratch.write(
                "rust/output.txt",
                &[OUTPUT_COLUMNS, &shared, &own_end, &mate_row].concat(),
            );
            let _ = fs::remove_dir_all(scratch.path("report"));
            fs::create_dir(scratch.path("report")).unwrap();
            let run = compare(
                &scratch,
                threads,
                &[&["--no-canonical-contigs"][..], names].concat(),
            );
            stdout_of_a_clean_run(&run);
            let summary = scratch.read("report/summary.json");
            let stated = scratch
                .read("report/summary.md")
                .contains(NOT_EVALUATED_LINE);
            assert!(
                !scratch.read("report/summary.md").contains("Not evaluated") || stated,
                "the line as the test states it"
            );
            (
                class_statistics(&summary, "bnd_synonym_mate_rows_lost"),
                stated,
            )
        };
        let (own_name, with_prefix, synonym) = ("N[Mito:700[", "N[chrMito:700[", "N[MT:700[");

        // The cache's names: its own name is no synonym, with the prefix the reference
        // removes or without it, and a name the cache lacks is one.
        let cached = ["--assembly", "MadeUp1", "--json-cache", &cache];
        assert_eq!(class(own_name, &cached), (one_sided(0, 0, None), false));
        assert_eq!(class(with_prefix, &cached), (one_sided(0, 0, None), false));
        assert_eq!(class(synonym, &cached), (one_sided(1, 1, None), false));

        // The same names from a set the flag names, the filter keeping every row here.
        let named = ["--assembly", "MadeUp1", "--canonical-contigs", "1,II,Mito"];
        assert_eq!(class(own_name, &named), (one_sided(0, 0, None), false));
        assert_eq!(class(synonym, &named), (one_sided(1, 1, None), false));

        // No names, under another assembly or under none: no test is made, and the summary
        // counts the row that reached it.
        for unnamed in [&["--assembly", "MadeUp1"][..], &[]] {
            assert_eq!(class(own_name, unnamed), (one_sided(0, 0, Some(1)), true));
            assert_eq!(class(synonym, unnamed), (one_sided(0, 0, Some(1)), true));
        }

        // A human assembly keeps the human names, whatever the cache lists: `Mito` is a
        // synonym there and `MT` the cache's name.
        let human = ["--assembly", "GRCh38", "--json-cache", &cache];
        assert_eq!(class(own_name, &human), (one_sided(1, 1, None), false));
        assert_eq!(class(synonym, &human), (one_sided(0, 0, None), false));
    }
}

#[test]
fn a_canonical_contigs_value_that_names_no_contig_ends_the_run_on_its_usage_message() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("names-no-contig", threads);
        scratch.write("nothing.txt", "# nothing\n\n");
        scratch.write("discordant.tsv", DISCORDANT_COLUMNS);
        let listing = scratch.arg("nothing.txt");
        let compare = ["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o"];
        let by_class = [
            "--by-consequence-class",
            "--perl-dir",
            "p",
            "--discordant-in",
            "discordant.tsv",
            "--dataset",
            "s01",
            "--class-csv",
            "class.csv",
            "--sweep-id",
            "sweep-x",
            "--binary-md5",
            "abcd1234",
        ];
        let named_inline = format!("--canonical-contigs={listing}");
        let cases: [(&[&str], &[&str], &str); 4] = [
            (&compare, &["--canonical-contigs", ","], ","),
            (&compare, &["--canonical-contigs=,"], ","),
            (&compare, &[&named_inline], &listing),
            (&by_class, &["--canonical-contigs", ","], ","),
        ];
        for (mode, flag, value) in cases {
            let output = scratch.run(threads, &[mode, flag].concat());
            assert_eq!(output.status.code(), Some(1), "{flag:?}");
            assert_eq!(
                stderr(&output),
                format!("error: --canonical-contigs '{value}' names no contig\n"),
                "the message alone, as a missing mode argument is reported"
            );
            assert_eq!(stdout(&output), "");
        }
        assert!(
            scratch.entries("o").is_empty(),
            "the report directory is made and nothing is written to it"
        );
        assert!(!scratch.path("class.csv").exists());

        // A listing that is not read is a failure of the run, not a usage error.
        fs::write(scratch.path("not-text.txt"), [0xff, 0xfe, b'\n']).unwrap();
        let listing = scratch.arg("not-text.txt");
        let unread = scratch.run(
            threads,
            &[&compare[..], &["--canonical-contigs", &listing]].concat(),
        );
        assert_eq!(unread.status.code(), Some(1));
        let message = stderr(&unread);
        assert!(
            message.starts_with(&format!("ERROR: [vep-outputs-compare] {listing}: ")),
            "{message}"
        );
    }
}

/// A reference output of three tuples, one on a contig outside the human set, and the saved
/// discordant rows of a run against it: one tuple the reference alone wrote.
fn saved_reports(tag: &str, threads: Option<&str>) -> Scratch {
    let scratch = Scratch::new(tag, threads);
    let output = [
        OUTPUT_COLUMNS,
        "v\t21:100\tA\tG\tENST1\tTranscript\tmissense_variant\n",
        "v\t21:200\tT\tG\tENST2\tTranscript\tintron_variant,splice_region_variant\n",
        "v\tKI270728.1:5\tC\tG\tENST3\tTranscript\tintron_variant\n",
    ];
    scratch.write("perl/output.txt", &output.concat());
    let discordant = [
        DISCORDANT_COLUMNS,
        "output.txt\tmissing_in_rust\t21:100\tA\tENST1\tTranscript\tmissense_variant\n",
    ];
    scratch.write("discordant.tsv", &discordant.concat());
    scratch
}

/// One per-class run over the saved reports into `class.csv`, under the human contig filter
/// of GRCh37 unless `extra` says otherwise; returns its stdout.
fn per_class(scratch: &Scratch, threads: Option<&str>, extra: &[&str]) -> String {
    let (perl, discordant, class_csv) = (
        scratch.arg("perl"),
        scratch.arg("discordant.tsv"),
        scratch.arg("class.csv"),
    );
    let saved = [
        "--by-consequence-class",
        "--perl-dir",
        &perl,
        "--discordant-in",
        &discordant,
        "--dataset",
        "s01",
        "--class-csv",
        &class_csv,
        "--sweep-id",
        "sweep-x",
        "--binary-md5",
        "abcd1234",
        "--assembly",
        "GRCh37",
    ];
    stdout_of_a_clean_run(&scratch.run(threads, &[&saved, extra].concat()))
}

/// [`per_class`] with the totals cache of the scratch.
fn per_class_cached(scratch: &Scratch, threads: Option<&str>, extra: &[&str]) -> String {
    let cache = scratch.arg("cache/s01.csv");
    per_class(
        scratch,
        threads,
        &[&["--vep-totals-cache", &cache], extra].concat(),
    )
}

/// What a per-class run prints after its line on the cache, `vep_tuples` the tuples of the
/// reference it counted.
fn per_class_lines(scratch: &Scratch, vep_tuples: u64) -> String {
    format!(
        "wrote per-class counts: {} (3 rows vep-rs/s01)\n\
         terms: 3  VEP tuples: {vep_tuples}  VEP-only: 1  engine-only: 0\n\
         no term reaches 1000 VEP tuples on this dataset, so it contributes no per-dataset floor\n",
        scratch.arg("class.csv")
    )
}

fn cache_miss(scratch: &Scratch, vep_tuples: u64) -> String {
    format!(
        "VEP per-term totals: cache miss, wrote {}\n{}",
        scratch.arg("cache/s01.csv"),
        per_class_lines(scratch, vep_tuples)
    )
}

fn cache_hit(scratch: &Scratch) -> String {
    format!(
        "VEP per-term totals: cache hit {}\n{}",
        scratch.arg("cache/s01.csv"),
        per_class_lines(scratch, 3)
    )
}

/// The cache file of the reference output as it stands, under the contig filter `canonical`
/// and with `intron_variant` on `intron` tuples.
fn cache_file(scratch: &Scratch, canonical: &str, intron: u64) -> String {
    let meta = fs::metadata(scratch.path("perl/output.txt")).unwrap();
    let modified = meta.modified().unwrap();
    let nanos = modified
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    format!(
        "# vep_totals_fingerprint=output.txt:{}:{nanos}|canonical={canonical}\n\
         term,vep_tuples\nintron_variant,{intron}\nmissense_variant,1\nsplice_region_variant,1\n",
        meta.len()
    )
}

/// Sets the modification time of the reference output one second later.
fn touch_reference_output(scratch: &Scratch) {
    let path = scratch.path("perl/output.txt");
    let modified = fs::metadata(&path).unwrap().modified().unwrap();
    let output = fs::File::options().write(true).open(&path).unwrap();
    output
        .set_modified(modified + Duration::from_secs(1))
        .unwrap();
}

const CLASS_CSV: &str = "\
dataset,assembly,engine,sweep_id,binary_md5,term,vep_tuples,vep_only,engine_only,intersection,engine_tuples,precision,recall,f1\r
s01,GRCh37,vep-rs,sweep-x,abcd1234,intron_variant,1,0,0,1,1,1.000000,1.000000,1.000000\r
s01,GRCh37,vep-rs,sweep-x,abcd1234,missense_variant,1,1,0,0,0,0.000000,0.000000,0.000000\r
s01,GRCh37,vep-rs,sweep-x,abcd1234,splice_region_variant,1,0,0,1,1,1.000000,1.000000,1.000000\r
";

#[test]
fn a_second_per_class_run_reads_the_totals_from_the_cache_and_writes_the_same_csv() {
    for threads in THREAD_BOUNDS {
        let scratch = saved_reports("cache-hit", threads);
        let first = per_class_cached(&scratch, threads, &[]);
        assert_eq!(first, cache_miss(&scratch, 3));
        // The tuple on the contig outside the filter is not counted: `intron_variant` is 1.
        assert_eq!(
            scratch.read("cache/s01.csv"),
            cache_file(&scratch, &HUMAN_CONTIGS.join(","), 1)
        );
        let first_csv = scratch.read("class.csv");
        assert_eq!(first_csv, CLASS_CSV);
        let second = per_class_cached(&scratch, threads, &[]);
        assert_eq!(second, cache_hit(&scratch));
        assert_eq!(scratch.read("class.csv"), first_csv);
    }
}

#[test]
fn a_reference_output_modified_since_the_cache_was_written_is_counted_again() {
    for threads in THREAD_BOUNDS {
        let scratch = saved_reports("cache-touched", threads);
        per_class_cached(&scratch, threads, &[]);
        touch_reference_output(&scratch);
        let again = per_class_cached(&scratch, threads, &[]);
        assert_eq!(again, cache_miss(&scratch, 3));
        assert_eq!(
            scratch.read("cache/s01.csv"),
            cache_file(&scratch, &HUMAN_CONTIGS.join(","), 1)
        );
    }
}

#[test]
fn totals_cached_under_one_contig_filter_are_not_read_under_another() {
    for threads in THREAD_BOUNDS {
        let scratch = saved_reports("cache-filter", threads);
        per_class_cached(&scratch, threads, &[]);
        let unfiltered = per_class_cached(&scratch, threads, &["--no-canonical-contigs"]);
        assert_eq!(unfiltered, cache_miss(&scratch, 4));
        assert_eq!(
            scratch.read("cache/s01.csv"),
            cache_file(&scratch, "none", 2),
            "the tuple outside the human set is counted"
        );
    }
}

#[test]
fn a_per_class_run_without_the_cache_option_reads_and_writes_no_cache() {
    for threads in THREAD_BOUNDS {
        let scratch = saved_reports("cache-absent", threads);
        let printed = per_class(&scratch, threads, &[]);
        assert_eq!(printed, per_class_lines(&scratch, 3), "no line on a cache");
        assert!(!scratch.path("cache").exists());
        assert_eq!(scratch.read("class.csv"), CLASS_CSV);
    }
}

/// Holds a compare run given `flag` with `path` to end as a usage error that says the path
/// `is_not`, before a report directory is made.
fn assert_input_path_refused(scratch: &Scratch, flag: &str, path: &str, is_not: &str) {
    let paths = ["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o"];
    let output = scratch.run(None, &[&paths[..], &[flag, path]].concat());
    assert_eq!(output.status.code(), Some(2), "{flag} {path}");
    let message = stderr(&output);
    assert!(
        message.ends_with(&format!(": error: {flag} {path} is not {is_not}\n")),
        "{message}"
    );
    assert_eq!(stdout(&output), "");
    assert!(!scratch.path("o").exists());
}

#[test]
fn an_input_vcf_that_is_not_there_is_a_usage_error() {
    let scratch = Scratch::new("input-vcf-absent", None);
    assert_input_path_refused(
        &scratch,
        "--input-vcf",
        &scratch.arg("absent.vcf"),
        "a file",
    );
}

#[test]
fn an_input_directory_that_is_not_there_is_a_usage_error() {
    let scratch = Scratch::new("input-dir-absent", None);
    assert_input_path_refused(
        &scratch,
        "--input-dir",
        &scratch.arg("absent"),
        "a directory",
    );
}

#[test]
fn an_input_vcf_that_is_a_directory_is_a_usage_error() {
    let scratch = Scratch::new("input-vcf-directory", None);
    assert_input_path_refused(
        &scratch,
        "--input-vcf",
        &scratch.0.to_string_lossy(),
        "a file",
    );
}

#[test]
fn an_input_directory_that_is_a_file_is_a_usage_error() {
    let scratch = Scratch::new("input-dir-file", None);
    scratch.write("absent.vcf", "");
    assert_input_path_refused(
        &scratch,
        "--input-dir",
        &scratch.arg("absent.vcf"),
        "a directory",
    );
}
