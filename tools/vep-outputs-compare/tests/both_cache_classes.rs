// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary on one output pair that holds rows of both classes that read a cache: the
//! three records of the reproduction `merged_same_dbid_collapse`, and two records under a
//! name the synonyms file of a made-up Ensembl VEP cache gives to two of its sequences. A
//! comparison that names the converted cache, the Ensembl VEP cache, both or neither applies
//! each class where its cache is named, with the statistics the class has alone, and the
//! reports state `merged_same_dbid_collapse`, then `merged_same_stable_id_collapse`, the
//! other class of the converted cache, which has no row here, then
//! `synonym_of_several_sequences_by_hash_order`, after every class that reads no cache.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is held under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 2] = [Some("1"), None];

const DBID_CLASS: &str = "merged_same_dbid_collapse";
const STABLE_ID_CLASS: &str = "merged_same_stable_id_collapse";
const SYNONYM_CLASS: &str = "synonym_of_several_sequences_by_hash_order";

const OUTPUT_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
    Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\t\
    Existing_variation\tExtra\n";

/// The command line of the reference's header: the three records of the reproduction in one
/// buffer under `--merged`, and no `--synonyms`.
const COMMAND: &str = "vep --assembly EXAMPLE1 --cache --cache_version 116 --database 0 \
    --dir_cache /cache --dir_plugins /plugins --force_overwrite --fork --format vcf \
    --input_file [PATH]/records.vcf --merged --no_stats --offline --output_file [PATH]/out.txt \
    --quiet --species two_set_example";

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

/// The reference's rows of the two records under the name `shared`, which vep-rs skips: one
/// between genes and one on two transcripts.
const UNDER_THE_SHARED_NAME: &str = "between_genes\tshared:399\tG\t-\t-\t-\tintergenic_variant\t\
    -\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n\
    on_two_transcripts\tshared:800\tG\tENSEXG00000000007\tENSEXT00000000007\tTranscript\t\
    intron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl\n\
    on_two_transcripts\tshared:800\tG\tENSEXG00000000008\tENSEXT00000000008\tTranscript\t\
    upstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=10;STRAND=1;SOURCE=Ensembl\n";

const RECORDS: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    contig\t101500\tshared_locus\tA\tG\t.\t.\t.\n\
    contig\t501000\tfirst_of_far_pair\tA\tG\t.\t.\t.\n\
    contig\t3201000\tsecond_of_far_pair\tA\tG\t.\t.\t.\n\
    shared\t399\tbetween_genes\tA\tG\t.\t.\t.\n\
    shared\t800\ton_two_transcripts\tA\tG\t.\t.\t.\n";

/// A synonyms file that gives the name `shared` to the sequences `contig` and `second`.
const SYNONYMS: &str = "contig\tshared\nsecond\tshared\n";

/// The rows of `merged_same_dbid_collapse` in `discordant.tsv`: the reference's intergenic
/// row and vep-rs's row of each dropped transcript.
const DBID_ROWS: [&str; 3] = [
    "records.txt\tmissing_in_rust\tcontig:3201000\tG\t-\t-\tintergenic_variant",
    "records.txt\textra_in_rust\tcontig:101500\tG\tENSEXT00000000001\tTranscript\t\
     non_coding_transcript_exon_variant",
    "records.txt\textra_in_rust\tcontig:3201000\tG\tXR_000002.1\tTranscript\t\
     non_coding_transcript_exon_variant",
];

/// The rows of `synonym_of_several_sequences_by_hash_order` in `discordant.tsv`: the
/// reference's three rows of the two records under the shared name.
const SYNONYM_ROWS: [&str; 3] = [
    "records.txt\tmissing_in_rust\tshared:399\tG\t-\t-\tintergenic_variant",
    "records.txt\tmissing_in_rust\tshared:800\tG\tENSEXT00000000007\tTranscript\tintron_variant",
    "records.txt\tmissing_in_rust\tshared:800\tG\tENSEXT00000000008\tTranscript\t\
     upstream_gene_variant",
];

/// The statistics of each class over those rows, as [`statistics_of`] reads them.
const DBID_STATISTICS: &str = "\"merged_same_dbid_collapse\":\
    {\"count\":3,\"excluded\":3,\"records\":2,\"perl\":1,\"rust\":2}";
const SYNONYM_STATISTICS: &str = "\"synonym_of_several_sequences_by_hash_order\":\
    {\"count\":3,\"excluded\":3,\"records\":2,\"perl\":3,\"rust\":0}";

/// The statistics of `merged_same_stable_id_collapse` where the comparison names the
/// converted cache: no two transcripts of that cache share a stable id.
const NO_STABLE_ID_ROW: &str = "\"merged_same_stable_id_collapse\":\
    {\"count\":0,\"excluded\":0,\"records\":0,\"perl\":0,\"rust\":0}";

/// The converted cache of the reproduction `merged_same_dbid_collapse`.
fn converted_cache() -> PathBuf {
    let cache = Path::new(env!("CARGO_MANIFEST_DIR")).join(
        "../../scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/json_cache",
    );
    fs::canonicalize(cache).expect("the reproduction's cache")
}

/// Copies the directory `from`, with everything under it, to `to`.
fn copy_directory(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_directory(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// A directory of one test under one thread bound, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, threads: Option<&str>) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-both-caches-{}-{tag}-t{}",
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

    /// Writes the two outputs and the input: the reference's rows of the five records under
    /// `COMMAND`, and vep-rs's four transcript rows of the first three.
    fn outputs(&self) {
        let reference = format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n## VEP command-line: {COMMAND}\n\
             {OUTPUT_COLUMNS}{SHARED_LOCUS_REFSEQ}{FIRST_OF_FAR_PAIR}\
             {SECOND_OF_FAR_PAIR_INTERGENIC}{UNDER_THE_SHARED_NAME}"
        );
        let vep_rs = format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n{OUTPUT_COLUMNS}{SHARED_LOCUS_ENSEMBL}\
             {SHARED_LOCUS_REFSEQ}{FIRST_OF_FAR_PAIR}{SECOND_OF_FAR_PAIR}"
        );
        self.write("reference/records.txt", &reference);
        self.write("vep_rs/records.txt", &vep_rs);
        self.write("records.vcf", RECORDS);
    }

    /// An Ensembl VEP cache directory that holds the sequences `contig` and `second` and
    /// the synonyms file `SYNONYMS`.
    fn ensembl_cache(&self) -> String {
        for sequence in ["contig", "second"] {
            fs::create_dir_all(self.0.join("ensembl_cache").join(sequence)).unwrap();
        }
        self.write("ensembl_cache/chr_synonyms.txt", SYNONYMS);
        self.arg("ensembl_cache")
    }

    /// One compare run against the 116.3 reference into an empty report directory, with
    /// `extra` after the paths.
    fn compare(&self, threads: Option<&str>, extra: &[&str]) -> Output {
        let (reference, vep_rs, report, records) = (
            self.arg("reference"),
            self.arg("vep_rs"),
            self.arg("report"),
            self.arg("records.vcf"),
        );
        let _ = fs::remove_dir_all(self.0.join("report"));
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

    /// The data rows the last comparison left open, sorted.
    fn open_rows(&self) -> Vec<String> {
        let open = self.read("report/discordant_open.tsv");
        let mut rows: Vec<String> = open.lines().skip(1).map(str::to_string).collect();
        rows.sort();
        rows
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

/// `rows`, sorted, as [`Scratch::open_rows`] returns the rows left open.
fn sorted(rows: &[&str]) -> Vec<String> {
    let mut rows: Vec<String> = rows.iter().map(|row| row.to_string()).collect();
    rows.sort();
    rows
}

/// The statistics the summary states for `class`, its white space dropped.
fn statistics_of(summary: &str, class: &str) -> String {
    let opening = format!("\"{class}\": {{");
    let at = summary
        .find(&opening)
        .unwrap_or_else(|| panic!("the summary states no statistics of {class}"));
    let block = &summary[at..];
    let end = block.find('}').expect("the end of the statistics");
    block[..=end].split_whitespace().collect()
}

/// The keys of the summary's object `name`, in the order the summary writes them.
fn keys_of(summary: &str, name: &str) -> Vec<String> {
    let opening = format!("\"{name}\": {{");
    let indent_of = |line: &str| line.len() - line.trim_start().len();
    let mut lines = summary.lines();
    let opened = lines
        .by_ref()
        .find(|line| line.trim_start() == opening)
        .unwrap_or_else(|| panic!("the summary holds no {name}"));
    let indent = indent_of(opened);
    let key_indent = " ".repeat(indent + 2);
    lines
        .take_while(|line| indent_of(line) > indent)
        .filter_map(|line| {
            let key = line.strip_prefix(&key_indent)?.strip_prefix('"')?;
            key.strip_suffix("\": {")
        })
        .map(str::to_string)
        .collect()
}

/// The classes of the summary's table, in the order of its rows.
fn table_classes(markdown: &str) -> Vec<String> {
    markdown
        .lines()
        .filter_map(|line| line.strip_prefix("| `"))
        .filter_map(|row| row.split_once('`'))
        .map(|(class, _)| class.to_string())
        .collect()
}

/// Asserts that the reports state the three classes last, the two of the converted cache
/// before the one of the Ensembl VEP cache: the statistics and the definitions of the
/// summary and the rows of the table.
fn assert_the_cache_classes_are_stated_last(summary: &str, markdown: &str) {
    let last_three = [
        DBID_CLASS.to_string(),
        STABLE_ID_CLASS.to_string(),
        SYNONYM_CLASS.to_string(),
    ];
    let stated = [
        (
            "excluded_categories",
            keys_of(summary, "excluded_categories"),
        ),
        ("divergence_classes", keys_of(summary, "divergence_classes")),
        ("the table of summary.md", table_classes(markdown)),
    ];
    for (place, classes) in stated {
        assert!(classes.ends_with(&last_three), "{place}: {classes:?}");
    }
}

#[test]
fn each_class_is_applied_where_its_cache_is_named_and_all_are_stated_in_registry_order() {
    let converted = converted_cache().to_string_lossy().into_owned();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("named", threads);
        scratch.outputs();
        let ensembl = scratch.ensembl_cache();
        let reports = ["report/summary.json", "report/summary.md"];

        let output = scratch.compare(threads, &[]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            sorted(&[DBID_ROWS, SYNONYM_ROWS].concat()),
            "neither cache"
        );
        for report in reports {
            let text = scratch.read(report);
            assert!(!text.contains(DBID_CLASS), "neither cache: {report}");
            assert!(!text.contains(STABLE_ID_CLASS), "neither cache: {report}");
            assert!(!text.contains(SYNONYM_CLASS), "neither cache: {report}");
        }

        let output = scratch.compare(threads, &["--json-cache", &converted]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            sorted(&SYNONYM_ROWS),
            "the converted cache"
        );
        let summary = scratch.read("report/summary.json");
        assert_eq!(statistics_of(&summary, DBID_CLASS), DBID_STATISTICS);
        assert_eq!(statistics_of(&summary, STABLE_ID_CLASS), NO_STABLE_ID_ROW);
        for report in reports {
            let text = scratch.read(report);
            assert!(
                !text.contains(SYNONYM_CLASS),
                "the converted cache: {report}"
            );
        }

        let output = scratch.compare(threads, &["--ensembl-cache", &ensembl]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), "");
        assert_eq!(
            scratch.open_rows(),
            sorted(&DBID_ROWS),
            "the Ensembl VEP cache"
        );
        let summary = scratch.read("report/summary.json");
        assert_eq!(statistics_of(&summary, SYNONYM_CLASS), SYNONYM_STATISTICS);
        for report in reports {
            let text = scratch.read(report);
            assert!(
                !text.contains(DBID_CLASS),
                "the Ensembl VEP cache: {report}"
            );
            assert!(
                !text.contains(STABLE_ID_CLASS),
                "the Ensembl VEP cache: {report}"
            );
        }

        // The order of the two options on the command line changes no report.
        for caches in [
            ["--json-cache", &converted, "--ensembl-cache", &ensembl],
            ["--ensembl-cache", &ensembl, "--json-cache", &converted],
        ] {
            let output = scratch.compare(threads, &caches);
            assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
            assert_eq!(stderr(&output), "");
            assert_eq!(scratch.open_rows(), sorted(&[]), "{caches:?}");
            let summary = scratch.read("report/summary.json");
            let markdown = scratch.read("report/summary.md");
            assert_eq!(statistics_of(&summary, DBID_CLASS), DBID_STATISTICS);
            assert_eq!(statistics_of(&summary, STABLE_ID_CLASS), NO_STABLE_ID_ROW);
            assert_eq!(statistics_of(&summary, SYNONYM_CLASS), SYNONYM_STATISTICS);
            assert!(
                summary.contains("\"excluded_perl_tuples\": 4,"),
                "{summary}"
            );
            assert!(
                summary.contains("\"excluded_rust_tuples\": 2,"),
                "{summary}"
            );
            assert_the_cache_classes_are_stated_last(&summary, &markdown);
            assert!(
                markdown.contains(
                    "- Divergence buckets: merged_same_dbid_collapse: 3 (3 excluded), \
                     synonym_of_several_sequences_by_hash_order: 3 (3 excluded)\n"
                ),
                "{markdown}"
            );
            assert!(
                stdout(&output).contains(
                    "adjusted: precision=1.000000 recall=1.000000 f1=1.000000 (excluded 4P + 2R = 6 \
                     intended-divergence tuples)\n"
                ),
                "{}",
                stdout(&output)
            );
        }
    }
}

/// A converted cache carries the synonyms file of the cache it is converted from. Where
/// that file gives the records' name to two sequences the converted cache files, the
/// reference's rows do not tell which of the two it read those records on, so the replay
/// does not follow the run: neither class of the converted cache is applied to the file and
/// the run states why for each, while `synonym_of_several_sequences_by_hash_order` sets
/// aside the rows it sets aside without the converted cache.
#[test]
fn a_converted_cache_that_files_both_sequences_of_the_name_is_not_replayed() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("not-replayed", threads);
        scratch.outputs();
        let ensembl = scratch.ensembl_cache();
        copy_directory(&converted_cache(), &scratch.0.join("converted"));
        fs::create_dir_all(scratch.0.join("converted/transcripts/second")).unwrap();
        scratch.write("converted/chr_synonyms.txt", SYNONYMS);
        let converted = scratch.arg("converted");

        let caches = ["--json-cache", &converted, "--ensembl-cache", &ensembl];
        let output = scratch.compare(threads, &caches);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            stderr(&output),
            "WARNING: [vep-outputs-compare] records.txt: merged_same_dbid_collapse is not \
             applied: the cache files two synonyms of the sequence name shared\n\
             WARNING: [vep-outputs-compare] records.txt: merged_same_stable_id_collapse is not \
             applied: the cache files two synonyms of the sequence name shared\n"
        );
        assert_eq!(scratch.open_rows(), sorted(&DBID_ROWS));
        let summary = scratch.read("report/summary.json");
        assert_eq!(
            statistics_of(&summary, DBID_CLASS),
            "\"merged_same_dbid_collapse\":\
             {\"count\":0,\"excluded\":0,\"records\":0,\"perl\":0,\"rust\":0}"
        );
        assert_eq!(statistics_of(&summary, STABLE_ID_CLASS), NO_STABLE_ID_ROW);
        assert_eq!(statistics_of(&summary, SYNONYM_CLASS), SYNONYM_STATISTICS);
        assert_the_cache_classes_are_stated_last(&summary, &scratch.read("report/summary.md"));
    }
}
