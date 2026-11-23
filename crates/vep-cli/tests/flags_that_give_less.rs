// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A flag or a subcommand that gives less than Ensembl VEP's says so. Each
//! test asks for the part vep-rs does not give and reads, on standard error,
//! the one line that names it or the refusal that stops the run; a run that
//! does not ask prints no such line.

mod common;

mod check_svs {
    //! `--check_svs <file>` names the known structural variants of the file it is
    //! given, and Ensembl's file lacks three kinds of mapping its database holds
    //! (`docs/output/known-variants.md`, "What the dump lacks"). A run under the
    //! flag says so in one `INFO:` line on standard error and writes the rows the
    //! committed reference holds; a run without the flag, a `--quiet` run and a
    //! run whose file is refused print no such line.

    use std::path::Path;
    use std::process::Output;

    use crate::common::{corpora, isolated_vep, runs, Corpus};

    /// The line, as a run prints it.
    const NOTICE: &str = "INFO: --check_svs reads names from the file it is given; Ensembl's \
        structural variation file lacks three kinds its database holds: variants on scaffolds and \
        patches, the further locations of a variant mapped to several places, and variants that \
        start at base 1";

    /// The flag set of the corpus that passes the flag with the corpus's slice of
    /// Ensembl's file.
    const FLAG_SET: &str = "check-svs";

    /// The `GRCh38-checksvs` corpus: records on chromosome 21, a slice of
    /// Ensembl's structural-variant file and the rows Ensembl VEP 116.3 wrote for
    /// them with a connection to its database.
    fn corpus() -> Corpus {
        corpora()
            .into_iter()
            .find(|corpus| corpus.release == "116" && corpus.name == "GRCh38-checksvs")
            .expect("tests/golden/116/GRCh38-checksvs")
    }

    /// Runs `vep` in the corpus directory on the corpus input against the corpus
    /// cache, writing `out`, with `flags` and without `--quiet` unless `flags`
    /// carry it; `dir` is the test's own output directory.
    fn run(corpus: &Corpus, dir: &Path, out: &Path, flags: &[String]) -> Output {
        isolated_vep(dir)
            .current_dir(&corpus.dir)
            .arg("-i")
            .arg(corpus.dir.join("variants.vcf"))
            .arg("-o")
            .arg(out)
            .arg("--offline")
            .args(corpus.cache_args())
            .args(["--species", corpus.species.as_str()])
            .args(["--assembly", corpus.assembly.as_str()])
            .args(["--buffer_size", "5000", "--force_overwrite", "--no_stats"])
            .args(flags)
            .output()
            .expect("vep binary runs")
    }

    /// The flags of the corpus's `check-svs` flag set.
    fn flag_set_flags(corpus: &Corpus) -> Vec<String> {
        runs(corpus)
            .into_iter()
            .find(|run| run.name == FLAG_SET)
            .expect("the check-svs flag set")
            .flags
    }

    fn text(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).into_owned()
    }

    /// The lines of an output that are not header lines.
    fn data_lines(output: &str) -> Vec<&str> {
        output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect()
    }

    /// A run under `--check_svs` prints the line once on standard error and
    /// nothing of it on standard output, in the output file or in a warnings
    /// file. Its data lines are the reference's, byte for byte: the rows Ensembl
    /// VEP 116.3 wrote for the same records from its database
    /// (`expected_check-svs/default.txt`).
    #[test]
    fn a_run_names_the_three_kinds_once_and_writes_the_reference_rows() {
        let corpus = corpus();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("with-the-flag.txt");
        let output = run(&corpus, dir.path(), &out, &flag_set_flags(&corpus));
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{stderr}");
        assert_eq!(
            stderr.lines().filter(|line| *line == NOTICE).count(),
            1,
            "{stderr}"
        );
        assert_eq!(stderr.matches("INFO: --check_svs").count(), 1, "{stderr}");
        assert!(output.stdout.is_empty(), "{}", text(&output.stdout));

        let written = std::fs::read_to_string(&out).unwrap();
        assert!(!written.contains("structural variation file"), "{written}");
        let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display()));
        assert!(
            !warnings.unwrap_or_default().contains("--check_svs"),
            "the line is no warning"
        );

        let reference = runs(&corpus)
            .into_iter()
            .find(|run| run.name == FLAG_SET)
            .expect("the check-svs flag set")
            .expected("default.txt");
        let rows = data_lines(&written);
        assert!(rows.iter().any(|row| row.contains("SV=")), "{written}");
        assert_eq!(rows, data_lines(&reference));
    }

    /// `--quiet` withholds the line: no line of standard error names the flag or
    /// Ensembl's file, and the rows are the ones the run writes without it.
    #[test]
    fn quiet_withholds_the_line() {
        let corpus = corpus();
        let dir = tempfile::tempdir().unwrap();
        let flags = flag_set_flags(&corpus);
        let out = dir.path().join("spoken.txt");
        let spoken = run(&corpus, dir.path(), &out, &flags);
        assert_eq!(spoken.status.code(), Some(0), "{}", text(&spoken.stderr));

        let quiet_out = dir.path().join("quiet.txt");
        let mut quiet_flags = flags;
        quiet_flags.push("--quiet".to_string());
        let quiet = run(&corpus, dir.path(), &quiet_out, &quiet_flags);
        assert_eq!(quiet.status.code(), Some(0), "{}", text(&quiet.stderr));
        let stderr = text(&quiet.stderr);
        assert!(!stderr.contains("INFO: --check_svs"), "{stderr}");
        assert!(!stderr.contains("structural variation file"), "{stderr}");
        assert!(quiet.stdout.is_empty(), "{}", text(&quiet.stdout));
        assert_eq!(
            data_lines(&std::fs::read_to_string(&quiet_out).unwrap()),
            data_lines(&std::fs::read_to_string(&out).unwrap())
        );
    }

    /// A run without the flag opens no structural-variant file and prints
    /// nothing about one: no line of its standard error names the flag.
    #[test]
    fn a_run_without_the_flag_prints_no_such_line() {
        let corpus = corpus();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("without-the-flag.txt");
        let output = run(&corpus, dir.path(), &out, &[]);
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{stderr}");
        assert!(
            stderr.contains("Processed "),
            "the run is not quiet: {stderr}"
        );
        assert!(!stderr.contains("check_svs"), "{stderr}");
        assert!(!stderr.contains("structural variation file"), "{stderr}");
        assert!(output.stdout.is_empty(), "{}", text(&output.stdout));
    }

    /// A path that is not a file ends the run with the `check_svs` error line,
    /// exit status 1 and no output; the line about Ensembl's file is printed
    /// once a file is opened, so this run prints none.
    #[test]
    fn a_refused_file_prints_its_error_alone() {
        let corpus = corpus();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("refused.txt");
        let missing = dir.path().join("no-such-file.vcf.gz");
        let flags = ["--check_svs".to_string(), missing.display().to_string()];
        let output = run(&corpus, dir.path(), &out, &flags);
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{stderr}");
        let refusal = format!("ERROR: [check_svs] {} is not a file", missing.display());
        assert_eq!(stderr.lines().last(), Some(refusal.as_str()), "{stderr}");
        assert!(!stderr.contains("INFO: --check_svs"), "{stderr}");
        assert!(!stderr.contains("structural variation file"), "{stderr}");
        assert!(!out.exists());
    }
}

mod debug {
    //! `--debug` through the `vep` binary. Ensembl VEP reads the flag to print each
    //! plugins directory it adds to its search path and the file each plugin module
    //! is loaded from (`Runner.pm:920, 941` at release/116.3). vep-rs prints
    //! neither and says so in one `INFO:` line on standard error, which `--quiet`
    //! withholds; the rows are those of the run without the flag.

    use std::path::PathBuf;
    use std::process::Stdio;

    /// The line a run under `--debug` prints.
    const NOTICE: &str = "INFO: --debug lists no plugin directories or plugin files, the two \
        things Ensembl VEP prints under it";

    /// The GRCh38 golden corpus of release 116: its input records and its cache.
    fn corpus() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
    }

    /// Runs the binary on the corpus's records with the built-in plugin
    /// `TSSDistance` and `flags`, in an output directory of the run's own, and
    /// returns the output's data lines and the lines of standard error.
    fn annotate(flags: &[&str]) -> (Vec<String>, Vec<String>) {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.txt");
        let output = crate::common::isolated_vep(dir.path())
            .arg("-i")
            .arg(corpus().join("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(corpus().join("json_cache"))
            .args(["--assembly", "GRCh38", "--force_overwrite", "--no_stats"])
            .args(["--plugin", "TSSDistance"])
            .args(flags)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(output.status.success(), "{flags:?}: {stderr}");
        let rows = std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect();
        (rows, stderr.lines().map(str::to_string).collect())
    }

    /// The lines of standard error that open as the notice does.
    fn debug_lines(stderr: &[String]) -> Vec<&str> {
        stderr
            .iter()
            .map(String::as_str)
            .filter(|line| line.starts_with("INFO: --debug"))
            .collect()
    }

    /// A run under `--debug` that loads a plugin prints the line once, and its
    /// rows, the plugin's column among them, are those of the run without the
    /// flag.
    #[test]
    fn debug_says_once_that_it_lists_nothing_and_changes_no_row() {
        let (rows, stderr) = annotate(&["--debug"]);
        assert_eq!(debug_lines(&stderr), [NOTICE], "{stderr:#?}");
        let (rows_without, _) = annotate(&[]);
        assert!(
            rows.iter().any(|row| row.contains("TSSDistance=")),
            "the plugin annotates"
        );
        assert_eq!(rows, rows_without);
    }

    /// The same run without the flag prints its status lines and none about
    /// `--debug`.
    #[test]
    fn a_run_without_debug_prints_no_line_about_it() {
        let (_, stderr) = annotate(&[]);
        assert!(!stderr.is_empty(), "the run is not quiet");
        assert!(debug_lines(&stderr).is_empty(), "{stderr:#?}");
    }

    /// `--quiet` withholds the line, as it withholds every notice.
    #[test]
    fn quiet_withholds_the_debug_line() {
        let (_, stderr) = annotate(&["--debug", "--quiet"]);
        assert!(debug_lines(&stderr).is_empty(), "{stderr:#?}");
    }
}

mod recode {
    //! `vep recode` through the binary, on what the cache alone does not give.
    //! Ensembl VEP's `variant_recoder` reads the Ensembl database on every run
    //! (`VariantRecoder.pm:96-109` at release/116.3 presets `database`, `merged`
    //! and `lrg`): the RefSeq transcripts beside the Ensembl ones, the mapping of
    //! a variant between its sequence and an LRG (`Parser.pm:901-944`) and, under
    //! `--var_synonyms`, every synonym source (`OutputFactory.pm:955-959`), the
    //! `Archive dbSNP` one among them, which the dump leaves out of a cache
    //! (`Pipeline/DumpVEP/InitDump.pm:263-266`). `vep recode` says so in one
    //! `INFO:` line on standard error, once in every run that completes; its
    //! output is the one the run writes with the line dropped.

    use std::fs;
    use std::path::PathBuf;
    use std::process::{Output, Stdio};

    use serde_json::Value;

    /// The line of a run whose cache is read as the Ensembl transcript set.
    const ENSEMBL_SET: &str = "INFO: this cache is read as the Ensembl transcript set, and RefSeq \
        notation is written from a merged cache; no variant is mapped to an LRG, which needs the \
        Ensembl database";

    /// The line of such a run that writes `var_synonyms`.
    const ENSEMBL_SET_WITH_SYNONYMS: &str = "INFO: this cache is read as the Ensembl transcript \
        set, and RefSeq notation is written from a merged cache; no variant is mapped to an LRG \
        and no Archive dbSNP synonym is written, which need the Ensembl database";

    /// The line of a run that reads the merged transcript set.
    const MERGED_SET: &str =
        "INFO: no variant is mapped to an LRG, which needs the Ensembl database";

    /// The line of such a run that writes `var_synonyms`.
    const MERGED_SET_WITH_SYNONYMS: &str = "INFO: no variant is mapped to an LRG and no Archive \
        dbSNP synonym is written, which need the Ensembl database";

    /// The line a run over the `GRCh38-merged` corpus writes before its own
    /// where its fields name `id`, as the default fields do: the corpus's
    /// cache, named `json_cache` from its directory, holds no variation data.
    const NO_ID: &str = "INFO: no variation data under json_cache: id, read from a cache's known \
        variants, stays unwritten; a cache from vep-cache-builder --auto or --include-variations, \
        or a conversion of Ensembl VEP's cache, holds it (docs/cache-variation-data.md)";

    /// The line of such a run that names `var_synonyms` as well.
    const NO_ID_AND_NO_SYNONYMS: &str = "INFO: no variation data under json_cache: id and \
        var_synonyms, read from a cache's known variants, stay unwritten; a cache from \
        vep-cache-builder --auto or --include-variations, or a conversion of Ensembl VEP's cache, \
        holds it (docs/cache-variation-data.md)";

    /// The line the annotation run over the `GRCh38-merged` corpus keeps for
    /// the tool, which writes it before its own: the corpus holds 12 insertion
    /// and deletion alleles with a notation on a model edited away from the
    /// genome (`input_left_out_notices.rs` holds the line to the one `vep`
    /// prints over the corpus).
    const MERGED_CORPUS_EDITED_MODELS: &str = "INFO: 12 insertion or deletion alleles have an \
        HGVS notation on a transcript whose sequence is edited away from the genome: the notation \
        is shifted 3' along the genome, and along the transcript's own sequence only on bases the \
        genome lacks";

    /// The cache, FASTA and synonyms flags of the `GRCh38-recoder` corpus, as its
    /// own command names them from the corpus directory.
    const RECODER_RESOURCES: [&str; 6] = [
        "--json_cache",
        "json_cache",
        "--fasta",
        "reference.fa.gz",
        "--synonyms",
        "chr_synonyms.txt",
    ];

    /// The golden corpus `name` of release 116.
    fn corpus(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116")
            .join(name)
    }

    /// Runs `vep recode` with `args` from the directory of the corpus `name`,
    /// under a home of the run's own, with standard error on `error_stream`.
    fn recode(name: &str, args: &[&str], error_stream: Stdio) -> Output {
        let home = tempfile::tempdir().unwrap();
        crate::common::isolated_vep(home.path())
            .arg("recode")
            .args(args)
            .current_dir(corpus(name))
            .stdin(Stdio::null())
            .stderr(error_stream)
            .output()
            .unwrap()
    }

    fn text(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).into_owned()
    }

    /// A run over the nine identifiers of the `GRCh38-recoder` corpus, whose
    /// cache slice is read as the Ensembl transcript set, prints the line once.
    /// Its output is the output of the same command on a standard error that
    /// takes no bytes, where the line is dropped, and a one-identifier run writes
    /// the output `recode_tool.rs::fields_selects_keys_and_pretty_indents` holds.
    #[test]
    fn an_ensembl_set_run_prints_the_line_once_and_changes_no_output_byte() {
        let args = [&["-i", "variants.id.txt"][..], &RECODER_RESOURCES[..]].concat();
        let run = recode("GRCh38-recoder", &args, Stdio::piped());
        assert!(run.status.success(), "{}", text(&run.stderr));
        assert_eq!(text(&run.stderr), format!("{ENSEMBL_SET}\n"));
        let objects: Vec<Value> = serde_json::from_slice(&run.stdout).unwrap();
        assert_eq!(objects.len(), 9, "one object for each input line");
        for stream in crate::common::unwritable_streams() {
            let dropped = recode("GRCh38-recoder", &args, crate::common::unwritable(stream));
            assert!(dropped.status.success(), "{stream}: {}", dropped.status);
            assert_eq!(dropped.stdout, run.stdout, "{stream}");
        }

        let args = [
            &["--id", "rs781261918", "--fields", "hgvsg"][..],
            &RECODER_RESOURCES[..],
        ]
        .concat();
        let run = recode("GRCh38-recoder", &args, Stdio::piped());
        assert_eq!(text(&run.stderr), format!("{ENSEMBL_SET}\n"));
        assert_eq!(
            text(&run.stdout),
            "[{\"T\":{\"hgvsg\":[\"NC_000021.9:g.25702257C>T\"],\"input\":\"rs781261918\"}}]"
        );
    }

    /// Under `--var_synonyms` the line names the `Archive dbSNP` synonyms as
    /// well, once: the run writes the synonyms of the cache's sources and none of
    /// that one.
    #[test]
    fn var_synonyms_adds_the_archive_dbsnp_synonyms_to_the_line() {
        let args = [
            &["-i", "variants.id.txt", "--var_synonyms"][..],
            &RECODER_RESOURCES[..],
        ]
        .concat();
        let run = recode("GRCh38-recoder", &args, Stdio::piped());
        assert!(run.status.success(), "{}", text(&run.stderr));
        assert_eq!(text(&run.stderr), format!("{ENSEMBL_SET_WITH_SYNONYMS}\n"));
        let output = text(&run.stdout);
        assert!(output.contains("\"var_synonyms\":[\"ClinVar: "), "{output}");
        assert!(!output.contains("Archive dbSNP"), "{output}");
    }

    /// A run on the merged cache of the `GRCh38-merged` corpus writes RefSeq
    /// notation, and its line holds no word of the transcript set or of RefSeq:
    /// it names the LRG mapping alone, and the synonyms under `--var_synonyms`.
    /// Before it stand the one line the annotation run keeps for the tool
    /// over this corpus, the count of its notations on edited models, and the
    /// tool's line of the fields that read known variants, of which the
    /// corpus's cache holds none.
    #[test]
    fn a_merged_cache_run_prints_the_line_without_the_transcript_set() {
        let args = [
            "-i",
            "variants.vcf",
            "--json_cache",
            "json_cache",
            "--fasta",
            "reference.fa.gz",
        ];
        let run = recode("GRCh38-merged", &args, Stdio::piped());
        assert!(run.status.success(), "{}", text(&run.stderr));
        assert_eq!(
            text(&run.stderr),
            format!("{MERGED_CORPUS_EDITED_MODELS}\n{NO_ID}\n{MERGED_SET}\n")
        );
        assert!(!text(&run.stderr).contains("RefSeq"));
        assert!(!text(&run.stderr).contains("transcript set"));
        assert!(text(&run.stdout).contains("\"NM_"), "no RefSeq notation");

        let args = [&args[..], &["--var_synonyms"][..]].concat();
        let run = recode("GRCh38-merged", &args, Stdio::piped());
        assert!(run.status.success(), "{}", text(&run.stderr));
        assert_eq!(
            text(&run.stderr),
            format!(
                "{MERGED_CORPUS_EDITED_MODELS}\n{NO_ID_AND_NO_SYNONYMS}\n\
                 {MERGED_SET_WITH_SYNONYMS}\n"
            )
        );
    }

    /// Three runs that end before the output is written print their `ERROR:`
    /// line alone: one the command line refuses before anything is read, `--lrg`,
    /// one the annotation run refuses, an input file that does not exist, and one
    /// that ends after the annotation run, on a warnings file it cannot write.
    #[test]
    fn a_refused_run_prints_its_error_line_alone() {
        let args = [
            &["--id", "rs781261918", "--lrg"][..],
            &RECODER_RESOURCES[..],
        ]
        .concat();
        let run = recode("GRCh38-recoder", &args, Stdio::piped());
        assert_eq!(run.status.code(), Some(1));
        assert_eq!(
            text(&run.stderr),
            "ERROR: [recode] --lrg is not supported: it maps a variant to its LRG transcripts \
             through the Ensembl core database, which vep recode does not read\n"
        );
        assert!(run.stdout.is_empty());

        let args = [&["-i", "no_such_input.txt"][..], &RECODER_RESOURCES[..]].concat();
        let run = recode("GRCh38-recoder", &args, Stdio::piped());
        assert_eq!(run.status.code(), Some(1));
        let error = text(&run.stderr);
        assert!(
            error.starts_with("ERROR: [input] ") && error.matches('\n').count() == 1,
            "{error:?}"
        );
        assert!(run.stdout.is_empty());

        let dir = tempfile::tempdir().unwrap();
        let warnings = dir.path().join("absent").join("warnings.txt");
        let warnings = warnings.to_str().unwrap();
        let args = [
            &["--id", "no_such_identifier", "--warning_file", warnings][..],
            &RECODER_RESOURCES[..],
        ]
        .concat();
        let run = recode("GRCh38-recoder", &args, Stdio::piped());
        assert_eq!(run.status.code(), Some(1));
        assert_eq!(
            text(&run.stderr),
            format!("ERROR: [recode] Could not write to warnings file {warnings}\n")
        );
        assert!(run.stdout.is_empty());
    }

    /// The cache of the `GRCh38-alt-contigs` corpus holds an LRG region, as a
    /// cache Ensembl dumps does (`Pipeline/DumpVEP/CreateDumpJobs.pm:170`). Its 11
    /// records given on that region are annotated there, each keeping the VCF
    /// form and the transcript of the region, and the run prints the line of its
    /// cache's set, which speaks of the mapping and not of the notation.
    #[test]
    fn records_on_an_lrg_region_keep_their_own_notation_under_the_line() {
        let records: Vec<String> =
            fs::read_to_string(corpus("GRCh38-alt-contigs").join("variants.vcf"))
                .unwrap()
                .lines()
                .filter(|line| line.starts_with("LRG_1\t"))
                .map(str::to_string)
                .collect();
        assert_eq!(records.len(), 11);
        let input = records.join("\n");
        let args = [
            "--input_data",
            &input,
            "--format",
            "vcf",
            "--json_cache",
            "json_cache",
            "--fields",
            "transcript_id",
            "--vcf_string",
        ];
        let run = recode("GRCh38-alt-contigs", &args, Stdio::piped());
        assert!(run.status.success(), "{}", text(&run.stderr));
        assert_eq!(text(&run.stderr), format!("{ENSEMBL_SET}\n"));
        assert!(text(&run.stderr).contains("no variant is mapped to an LRG"));
        assert!(!text(&run.stderr).contains("LRG notation"));
        let objects: Vec<Value> = serde_json::from_slice(&run.stdout).unwrap();
        assert_eq!(objects.len(), 11);
        for object in &objects {
            let alleles = object.as_object().unwrap();
            assert!(!alleles.is_empty(), "{object}");
            for fields in alleles.values() {
                for (key, name) in [("vcf_string", "LRG_1-"), ("transcript_id", "LRG_1t")] {
                    let values = fields[key].as_array().unwrap();
                    assert!(!values.is_empty(), "{key}: {fields}");
                    for value in values {
                        assert!(value.as_str().unwrap().starts_with(name), "{key}: {value}");
                    }
                }
            }
        }
    }
}

mod transcript_filter {
    //! `--transcript_filter` through the `vep` binary, in the two places where a
    //! filter reads less than Ensembl VEP's.
    //!
    //! A cache built from a GTF or GFF3 carries neither `_source_cache` nor
    //! `dbID`, which the transcript Ensembl VEP builds from the file holds
    //! (`AnnotationSource/File/BaseGXF.pm:489, 635` at release/116.3). A filter
    //! that names either key over such a cache draws one `INFO:` line for the key
    //! on standard error before the input is read, which `--quiet` withholds, and
    //! the rows are the filter's own. A file read beside a cache gives its
    //! transcripts a `_source_cache` and no `dbID`, so a filter naming `dbID`
    //! draws a line for them over any cache. `vep recode`, whose annotation run
    //! is quiet, writes the lines of its `--transcript_filter` itself, before
    //! its own line, in a run that completes.
    //!
    //! A possessive quantifier in a pattern is Perl's atomic group around the
    //! greedy one, which the `regex` crate lacks: it is refused before the input
    //! is read, under `vep`, `vep filter` and `vep recode` alike.

    use std::path::{Path, PathBuf};
    use std::process::{Output, Stdio};

    /// The corpus whose cache `vep-cache-builder` wrote from a GTF: no transcript
    /// carries a `dbID` or a `_source_cache`.
    const BUILT: &str = "GRCh38-gtf-builder";

    /// The corpus whose cache was converted from the Ensembl cache: every
    /// transcript carries a `dbID` and none a `_source_cache`.
    const CONVERTED: &str = "GRCh38-transcript-filter";

    fn golden(corpus: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116")
            .join(corpus)
    }

    /// The line a filter naming `key` draws over a cache built from a GTF or GFF3.
    fn notice(key: &str) -> String {
        format!(
            "INFO: --transcript_filter names {key}, which no transcript of this cache carries: \
             a cache built from a GTF or GFF3 leaves it undefined, so exists and a comparison on \
             it fail for each of its transcripts"
        )
    }

    /// The line a filter naming `dbID` draws for the transcripts of a file read
    /// beside the cache.
    const FILE_NOTICE: &str =
        "INFO: --transcript_filter names dbID, which no transcript of a file \
        read beside the cache carries: a GTF or GFF3 read under --gtf or --gff leaves it \
        undefined, so exists and a comparison on it fail for each of its transcripts";

    /// The line that refuses `pattern` for its possessive `quantifier`.
    fn refusal(pattern: &str, quantifier: &str) -> String {
        format!(
            "ERROR: [filter] pattern \"{pattern}\" is not supported: \
             the possessive quantifier {quantifier} is not supported"
        )
    }

    /// One run: its exit status, the lines of standard error and the data rows of
    /// the output, `None` when the run wrote no output file.
    struct Run {
        code: Option<i32>,
        stderr: Vec<String>,
        rows: Option<Vec<String>>,
    }

    impl Run {
        /// The lines of standard error that open as the notice does.
        fn filter_lines(&self) -> Vec<&str> {
            self.stderr
                .iter()
                .map(String::as_str)
                .filter(|line| line.starts_with("INFO: --transcript_filter"))
                .collect()
        }

        fn rows(&self) -> &[String] {
            self.rows.as_deref().expect("the run wrote its output")
        }

        /// Whether the run wrote a row and every row is the one a record without
        /// a transcript writes.
        fn is_intergenic_alone(&self) -> bool {
            !self.rows().is_empty()
                && self
                    .rows()
                    .iter()
                    .all(|row| row.split('\t').nth(6) == Some("intergenic_variant"))
        }
    }

    /// Runs the binary against the cache of `corpus` on `input` with `flags`,
    /// writing the tab format under `dir`.
    fn start(dir: &Path, corpus: &str, input: &Path, flags: &[&str]) -> Run {
        let out = dir.join("out.txt");
        let output = crate::common::isolated_vep(dir)
            .arg("-i")
            .arg(input)
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(golden(corpus).join("json_cache"))
            .args(["--assembly", "GRCh38", "--no_stats", "--tab"])
            .args(flags)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let rows = std::fs::read_to_string(&out).ok().map(|text| {
            text.lines()
                .filter(|line| !line.starts_with('#'))
                .map(str::to_string)
                .collect()
        });
        Run {
            code: output.status.code(),
            stderr: String::from_utf8_lossy(&output.stderr)
                .lines()
                .map(str::to_string)
                .collect(),
            rows,
        }
    }

    /// [`start`] over the records of `corpus`, in an output directory of the
    /// run's own; the run must succeed.
    fn annotate(corpus: &str, flags: &[&str]) -> Run {
        let dir = tempfile::tempdir().unwrap();
        let run = start(
            dir.path(),
            corpus,
            &golden(corpus).join("variants.vcf"),
            flags,
        );
        assert_eq!(run.code, Some(0), "{flags:?}: {:#?}", run.stderr);
        run
    }

    /// A filter on `dbID` over the built cache draws the key's line once, as the
    /// field of a unit and as a `#field` value, and writes the rows the same
    /// command writes under `--quiet`: every record intergenic, since `exists`
    /// fails on each transcript. An input that is not there shows the line
    /// standing before the input is read: it precedes the refusal of the input,
    /// and no output is written.
    #[test]
    fn a_filter_on_dbid_over_a_built_cache_says_so_once_before_the_input_is_read() {
        let filter = ["--transcript_filter", "dbID exists"];
        let run = annotate(BUILT, &filter);
        assert_eq!(run.filter_lines(), [notice("dbID")], "{:#?}", run.stderr);
        assert!(run.is_intergenic_alone(), "{:#?}", run.rows());
        let quiet = annotate(BUILT, &["--transcript_filter", "dbID exists", "--quiet"]);
        assert_eq!(run.rows(), quiet.rows());

        let as_value = annotate(BUILT, &["--transcript_filter", "stable_id is #dbID"]);
        assert_eq!(as_value.filter_lines(), [notice("dbID")]);

        let dir = tempfile::tempdir().unwrap();
        let absent = dir.path().join("absent.vcf");
        let refused = start(dir.path(), BUILT, &absent, &filter);
        assert_eq!(refused.code, Some(1), "{:#?}", refused.stderr);
        let position = |opening: &str| {
            refused
                .stderr
                .iter()
                .position(|line| line.starts_with(opening))
                .unwrap_or_else(|| panic!("no {opening} line: {:#?}", refused.stderr))
        };
        assert!(position("INFO: --transcript_filter names dbID") < position("ERROR: "));
        assert!(refused.rows.is_none(), "no output is written");
    }

    /// A filter on `_source_cache` over the built cache draws that key's line
    /// alone; a filter naming both keys draws both, and its rows are the
    /// unfiltered ones, since `nex` and a negated `exists` pass every transcript;
    /// a filter naming neither key draws none.
    #[test]
    fn a_filter_on_source_cache_over_a_built_cache_says_so() {
        let run = annotate(BUILT, &["--transcript_filter", "_source_cache exists"]);
        assert_eq!(
            run.filter_lines(),
            [notice("_source_cache")],
            "{:#?}",
            run.stderr
        );
        assert!(run.is_intergenic_alone(), "{:#?}", run.rows());

        let unfiltered = annotate(BUILT, &[]);
        assert!(unfiltered.filter_lines().is_empty());
        assert!(!unfiltered.is_intergenic_alone());
        let both = annotate(
            BUILT,
            &[
                "--transcript_filter",
                "dbID nex",
                "--transcript_filter",
                "not _source_cache",
            ],
        );
        assert_eq!(
            both.filter_lines(),
            [notice("_source_cache"), notice("dbID")]
        );
        assert_eq!(both.rows(), unfiltered.rows());

        let neither = annotate(BUILT, &["--transcript_filter", "biotype is protein_coding"]);
        assert!(neither.filter_lines().is_empty(), "{:#?}", neither.stderr);
    }

    /// Over a converted cache the same filters draw no line and write the rows
    /// they write under `--quiet`: `dbID exists` the unfiltered rows, since every
    /// transcript carries one, and `_source_cache exists` an intergenic row for
    /// every record, as Ensembl VEP reads a cache that is not a merged one.
    #[test]
    fn the_same_filters_over_a_converted_cache_draw_no_line_and_write_their_rows() {
        let unfiltered = annotate(CONVERTED, &[]);
        assert!(!unfiltered.is_intergenic_alone());
        for (filter, keeps_every_transcript) in
            [("dbID exists", true), ("_source_cache exists", false)]
        {
            let run = annotate(CONVERTED, &["--transcript_filter", filter]);
            assert!(run.filter_lines().is_empty(), "{filter}: {:#?}", run.stderr);
            let quiet = annotate(CONVERTED, &["--transcript_filter", filter, "--quiet"]);
            assert_eq!(run.rows(), quiet.rows(), "{filter}");
            if keeps_every_transcript {
                assert_eq!(run.rows(), unfiltered.rows(), "{filter}");
            } else {
                assert!(run.is_intergenic_alone(), "{filter}");
            }
        }
    }

    /// `--quiet` withholds both lines, as it withholds every notice.
    #[test]
    fn quiet_withholds_the_lines_of_a_filter_on_an_absent_key() {
        let run = annotate(
            BUILT,
            &[
                "--transcript_filter",
                "dbID exists or _source_cache exists",
                "--quiet",
            ],
        );
        assert!(run.filter_lines().is_empty(), "{:#?}", run.stderr);
        assert!(run.is_intergenic_alone());
    }

    /// A GTF read beside the built cache gives its transcripts the file's name
    /// as `_source_cache`, so a filter on that key selects among the sources and
    /// draws no line for it; `dbID`, which the file's transcripts lack as the
    /// cache's do, draws the cache's line and the file's.
    #[test]
    fn a_file_read_beside_the_built_cache_draws_no_source_cache_line() {
        let fasta = golden(BUILT).join("reference.fa.gz");
        let gtf = golden(BUILT).join("chr21_genes.gtf.gz");
        let annotate_beside = |filter: &str| {
            annotate(
                BUILT,
                &[
                    "--fasta",
                    fasta.to_str().unwrap(),
                    "--gtf",
                    gtf.to_str().unwrap(),
                    "--transcript_filter",
                    filter,
                ],
            )
        };
        let of_file = annotate_beside("_source_cache is chr21_genes.gtf.gz");
        assert!(of_file.filter_lines().is_empty(), "{:#?}", of_file.stderr);
        assert!(
            !of_file.is_intergenic_alone(),
            "the file's transcripts pass"
        );

        let both = annotate_beside("_source_cache is chr21_genes.gtf.gz or dbID exists");
        assert_eq!(
            both.filter_lines(),
            [notice("dbID").as_str(), FILE_NOTICE],
            "{:#?}",
            both.stderr
        );
        assert_eq!(both.rows(), of_file.rows());
    }

    /// A GTF read beside a converted cache, whose transcripts carry a `dbID`:
    /// `dbID exists` keeps the cache's transcripts and drops the file's, so the
    /// rows are fewer than the unfiltered ones, and the run draws the file's
    /// line alone, which `--quiet` withholds while the rows stay.
    #[test]
    fn a_file_read_beside_a_converted_cache_draws_the_dbid_line_of_the_file() {
        let fasta = golden(BUILT).join("reference.fa.gz");
        let gtf = golden(BUILT).join("chr21_genes.gtf.gz");
        let input = golden(BUILT).join("variants.vcf");
        let annotate_beside = |flags: &[&str]| {
            let dir = tempfile::tempdir().unwrap();
            let beside = [
                "--fasta",
                fasta.to_str().unwrap(),
                "--gtf",
                gtf.to_str().unwrap(),
            ];
            let run = start(
                dir.path(),
                CONVERTED,
                &input,
                &[&beside[..], flags].concat(),
            );
            assert_eq!(run.code, Some(0), "{flags:?}: {:#?}", run.stderr);
            run
        };
        let unfiltered = annotate_beside(&[]);
        assert!(unfiltered.filter_lines().is_empty());

        let run = annotate_beside(&["--transcript_filter", "dbID exists"]);
        assert_eq!(run.filter_lines(), [FILE_NOTICE], "{:#?}", run.stderr);
        assert!(run.rows().len() < unfiltered.rows().len());

        let quiet = annotate_beside(&["--transcript_filter", "dbID exists", "--quiet"]);
        assert!(quiet.filter_lines().is_empty(), "{:#?}", quiet.stderr);
        assert_eq!(quiet.rows(), run.rows());

        let other_key = annotate_beside(&["--transcript_filter", "_source_cache exists"]);
        assert!(
            other_key.filter_lines().is_empty(),
            "{:#?}",
            other_key.stderr
        );
    }

    /// The line a completed `vep recode` run writes over a cache read as the
    /// Ensembl transcript set, the last line of its standard error.
    const RECODE_NOTICE: &str = "INFO: this cache is read as the Ensembl transcript set, and \
        RefSeq notation is written from a merged cache; no variant is mapped to an LRG, which \
        needs the Ensembl database";

    /// The line a `vep recode` run over the built cache writes before its own:
    /// its default fields name `id` and the cache, named `json_cache` from the
    /// corpus directory, holds no variation data.
    const RECODE_NO_ID: &str = "INFO: no variation data under json_cache: id, read from a cache's \
        known variants, stays unwritten; a cache from vep-cache-builder --auto or \
        --include-variations, or a conversion of Ensembl VEP's cache, holds it \
        (docs/cache-variation-data.md)";

    /// The input, cache and FASTA of the built corpus, named from its directory.
    const BUILT_RESOURCES: [&str; 6] = [
        "-i",
        "variants.vcf",
        "--json_cache",
        "json_cache",
        "--fasta",
        "reference.fa.gz",
    ];

    /// Runs `vep recode` with `args` from the directory of `corpus`, under a home
    /// of the run's own.
    fn recode(corpus: &str, args: &[&str]) -> Output {
        let home = tempfile::tempdir().unwrap();
        crate::common::isolated_vep(home.path())
            .arg("recode")
            .args(args)
            .current_dir(golden(corpus))
            .stdin(Stdio::null())
            .output()
            .unwrap()
    }

    /// [`recode`] over the records of the built corpus under each of `filters`:
    /// the whole standard error and the output of a run that must succeed.
    fn recode_built(filters: &[&str]) -> (String, Vec<u8>) {
        let flags: Vec<&str> = filters
            .iter()
            .flat_map(|filter| ["--transcript_filter", filter])
            .collect();
        let output = recode(BUILT, &[&BUILT_RESOURCES[..], &flags[..]].concat());
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert_eq!(output.status.code(), Some(0), "{filters:?}: {stderr}");
        (stderr, output.stdout)
    }

    /// The whole standard error of a completed `vep recode` run over the built
    /// cache whose filters name `keys`: the line of each key, then the tool's
    /// two, that of a cache without variation data and its own.
    fn recode_lines(keys: &[&str]) -> String {
        let mut lines: String = keys.iter().map(|key| notice(key) + "\n").collect();
        for line in [RECODE_NO_ID, RECODE_NOTICE] {
            lines.push_str(line);
            lines.push('\n');
        }
        lines
    }

    /// `vep recode` runs its annotation quiet, so the line of a filter on a key
    /// the built cache leaves undefined is the tool's to write: once for each key
    /// named, before the tool's two lines, and nothing else. The output is the
    /// one a filter that names neither key and drops every transcript writes;
    /// two filters that pass every transcript draw both lines over the
    /// unfiltered output, and a filter naming neither key draws none.
    #[test]
    fn vep_recode_writes_the_line_of_a_filter_on_an_absent_key_before_its_own() {
        let (stderr, unfiltered) = recode_built(&[]);
        assert_eq!(stderr, recode_lines(&[]));
        let (stderr, dropped) = recode_built(&["stable_id is no_such_transcript"]);
        assert_eq!(stderr, recode_lines(&[]));
        assert_ne!(dropped, unfiltered, "the filter drops the transcripts");

        let (stderr, stdout) = recode_built(&["dbID exists"]);
        assert_eq!(stderr, recode_lines(&["dbID"]));
        assert_eq!(stdout, dropped);

        let (stderr, stdout) = recode_built(&["_source_cache exists"]);
        assert_eq!(stderr, recode_lines(&["_source_cache"]));
        assert_eq!(stdout, dropped);

        let (stderr, _) = recode_built(&["stable_id is #dbID"]);
        assert_eq!(stderr, recode_lines(&["dbID"]));

        let (stderr, stdout) = recode_built(&["dbID nex", "not _source_cache"]);
        assert_eq!(stderr, recode_lines(&["_source_cache", "dbID"]));
        assert_eq!(stdout, unfiltered);
    }

    /// Over the converted cache of the `GRCh38-recoder` corpus, whose
    /// transcripts carry a `dbID`, the same filters draw no line: the tool's own
    /// is all a run writes to standard error. `dbID exists` keeps every
    /// transcript, so its output is the unfiltered one.
    #[test]
    fn vep_recode_over_a_converted_cache_writes_no_line_of_the_filter() {
        let resources = [
            "--id",
            "rs781261918",
            "--json_cache",
            "json_cache",
            "--fasta",
            "reference.fa.gz",
            "--synonyms",
            "chr_synonyms.txt",
        ];
        let unfiltered = recode("GRCh38-recoder", &resources);
        assert_eq!(unfiltered.status.code(), Some(0));
        for (filter, keeps_every_transcript) in
            [("dbID exists", true), ("_source_cache exists", false)]
        {
            let args = [&resources[..], &["--transcript_filter", filter]].concat();
            let run = recode("GRCh38-recoder", &args);
            assert_eq!(run.status.code(), Some(0), "{filter}");
            assert_eq!(
                String::from_utf8_lossy(&run.stderr),
                format!("{RECODE_NOTICE}\n"),
                "{filter}"
            );
            assert_eq!(
                run.stdout == unfiltered.stdout,
                keeps_every_transcript,
                "{filter}"
            );
        }
    }

    /// A `vep recode` run that ends before its output writes its `ERROR:` line
    /// alone, though its filter names a key the built cache leaves undefined: a
    /// run whose input the annotation run refuses, an identifier over a cache
    /// without variation files, and a run that ends after the annotation run, on
    /// a warnings file it cannot write.
    #[test]
    fn a_vep_recode_run_that_ends_before_its_output_writes_no_line_of_the_filter() {
        let filter = ["--transcript_filter", "dbID exists"];
        let resources = &BUILT_RESOURCES[2..];
        let args = [&["--id", "rs781261918"][..], resources, &filter[..]].concat();
        let refused = recode(BUILT, &args);
        assert_eq!(refused.status.code(), Some(1));
        let error = String::from_utf8_lossy(&refused.stderr);
        assert!(
            error.starts_with("ERROR: [input] ") && error.matches('\n').count() == 1,
            "{error:?}"
        );
        assert!(refused.stdout.is_empty());

        // A record on a sequence the cache lacks is skipped with a warning,
        // which the run then cannot write.
        let dir = tempfile::tempdir().unwrap();
        let records = std::fs::read_to_string(golden(BUILT).join("variants.vcf")).unwrap();
        let record = records
            .lines()
            .find(|line| !line.starts_with('#'))
            .expect("the corpus holds a record");
        let (_, columns) = record.split_once('\t').unwrap();
        let input = dir.path().join("off_the_cache.vcf");
        std::fs::write(&input, format!("no_such_sequence\t{columns}\n{record}\n")).unwrap();
        let warnings = dir.path().join("absent").join("warnings.txt");
        let warnings = warnings.to_str().unwrap();
        let args = [
            &["-i", input.to_str().unwrap(), "--format", "vcf"][..],
            resources,
            &filter[..],
            &["--warning_file", warnings][..],
        ]
        .concat();
        let unwritten = recode(BUILT, &args);
        assert_eq!(unwritten.status.code(), Some(1));
        assert_eq!(
            String::from_utf8_lossy(&unwritten.stderr),
            format!("ERROR: [recode] Could not write to warnings file {warnings}\n")
        );
        assert!(unwritten.stdout.is_empty());
    }

    /// A pattern with a possessive quantifier, in each spelling, ends the run
    /// with the one refusal line and exit status 1 before the input is read: the
    /// input named is not there, and no output is written. Under `--quiet` the
    /// line is all the run prints; without it the line follows the banner.
    #[test]
    fn a_possessive_quantifier_is_refused_before_the_input_is_read() {
        for (pattern, quantifier) in [
            ("s++s$", "++"),
            ("ENST0*+", "*+"),
            ("ENST0?+0", "?+"),
            ("ENST0{1,2}+0", "{1,2}+"),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let absent = dir.path().join("absent.vcf");
            let filter = format!("stable_id match {pattern}");
            let line = refusal(pattern, quantifier);
            let quiet = start(
                dir.path(),
                CONVERTED,
                &absent,
                &["--transcript_filter", &filter, "--quiet"],
            );
            assert_eq!(quiet.code, Some(1), "{pattern}: {:#?}", quiet.stderr);
            assert_eq!(quiet.stderr, [line.as_str()], "{pattern}");
            assert!(quiet.rows.is_none(), "{pattern}: no output is written");

            let run = start(
                dir.path(),
                CONVERTED,
                &absent,
                &["--transcript_filter", &filter],
            );
            assert_eq!(run.code, Some(1), "{pattern}: {:#?}", run.stderr);
            let errors: Vec<&String> = run
                .stderr
                .iter()
                .filter(|line| line.starts_with("ERROR: "))
                .collect();
            assert_eq!(errors, [&line], "{pattern}: {:#?}", run.stderr);
            assert_eq!(run.stderr.last(), Some(&line), "{pattern}");
            assert!(run.rows.is_none(), "{pattern}: no output is written");
        }
    }

    /// The same pattern without the marker runs: the greedy and the lazy
    /// quantifier draw no refusal, and `ENST0+` passes every transcript of the
    /// Ensembl cache, so the rows are the unfiltered ones.
    #[test]
    fn a_greedy_or_lazy_quantifier_runs_and_draws_no_refusal() {
        let unfiltered = annotate(CONVERTED, &[]);
        for pattern in ["ENST0+", "ENST0+?", "ENST0{1,2}"] {
            let filter = format!("stable_id match {pattern}");
            let run = annotate(CONVERTED, &["--transcript_filter", &filter]);
            assert!(
                run.stderr.iter().all(|line| !line.contains("possessive")),
                "{pattern}: {:#?}",
                run.stderr
            );
            assert_eq!(run.rows(), unfiltered.rows(), "{pattern}");
        }
    }

    /// `vep filter` compiles its filters through the same reading, so it refuses
    /// the quantifier with the same line before it opens its input.
    #[test]
    fn vep_filter_refuses_a_possessive_quantifier_before_it_opens_its_input() {
        let dir = tempfile::tempdir().unwrap();
        let absent = dir.path().join("absent.vcf");
        let output = crate::common::isolated_vep(dir.path())
            .arg("filter")
            .arg("-i")
            .arg(&absent)
            .args(["--filter", "SYMBOL match s++s$"])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!("{}\n", refusal("s++s$", "++"))
        );
        assert!(output.stdout.is_empty());
    }

    /// `vep recode` hands its `--transcript_filter` to the annotation run, whose
    /// refusal ends it with the same line and nothing on standard output.
    #[test]
    fn vep_recode_refuses_a_possessive_quantifier() {
        let dir = tempfile::tempdir().unwrap();
        let corpus = golden("GRCh38-recoder");
        let output = crate::common::isolated_vep(dir.path())
            .arg("recode")
            .args(["--id", "rs781261918"])
            .args(["--transcript_filter", "stable_id match s++s$"])
            .arg("--json_cache")
            .arg(corpus.join("json_cache"))
            .arg("--fasta")
            .arg(corpus.join("reference.fa.gz"))
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!("{}\n", refusal("s++s$", "++"))
        );
        assert!(output.stdout.is_empty());
    }
}

mod version_short {
    //! `-V` through the `vep` binary. It is the version's short here, where
    //! Getopt::Long's `ignore_case` reads it as `-v` and Ensembl VEP runs verbosely
    //! (`verbose|v`, `Config.pm:159` at release/116.3; `:62` leaves `ignore_case`
    //! on). Alone it prints the version. Beside any other argument the run is
    //! refused with one line that names `-v` and `--version`, before a file is
    //! opened; `-v` in its place annotates, and `--version` spelled out beside the
    //! same arguments prints the version.

    use std::fs;
    use std::path::{Path, PathBuf};
    use std::process::{Output, Stdio};

    /// The refusal, as the binary prints it.
    const REFUSAL: &str = "ERROR: [cli] -V beside other arguments is not supported: it prints \
        the version here, where Ensembl VEP reads it as -v (--verbose) and runs; write -v for a \
        verbose run, or --version for the version\n";

    /// The line the version prints.
    fn version_line() -> String {
        format!("vep {}\n", env!("CARGO_PKG_VERSION"))
    }

    /// The GRCh38 golden corpus of release 116: its input records and its cache.
    fn corpus() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
    }

    fn text(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).into_owned()
    }

    /// Runs the binary inside `dir`, the test's own output directory, on the
    /// corpus's records and cache with `-o <dir>/out.txt` and `tokens` beside them.
    fn run_beside_an_input_and_an_output(dir: &Path, tokens: &[&str]) -> Output {
        crate::common::isolated_vep(dir)
            .current_dir(dir)
            .arg("-i")
            .arg(corpus().join("variants.vcf"))
            .arg("-o")
            .arg(dir.join("out.txt"))
            .args(tokens)
            .arg("--json_cache")
            .arg(corpus().join("json_cache"))
            .args(["--assembly", "GRCh38", "--no_stats"])
            .stdin(Stdio::null())
            .output()
            .unwrap()
    }

    /// What `dir` holds besides the home directory the isolated binary is given.
    fn written(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .filter(|name| name != "home")
            .collect();
        names.sort();
        names
    }

    /// The data lines of an output file.
    fn rows(output_file: &Path) -> Vec<String> {
        fs::read_to_string(output_file)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect()
    }

    /// `-V` as the only argument prints the version line on standard output and
    /// nothing on standard error, exit 0.
    #[test]
    fn alone_prints_the_version_and_nothing_else() {
        let dir = tempfile::tempdir().unwrap();
        let output = crate::common::isolated_vep(dir.path())
            .current_dir(dir.path())
            .arg("-V")
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert_eq!(text(&output.stdout), version_line());
        assert_eq!(text(&output.stderr), "");
        assert_eq!(written(dir.path()), Vec::<String>::new());
    }

    /// `-V` beside an input and an output is refused: the one line on standard
    /// error, nothing on standard output, exit 1, and neither the output file nor
    /// any other file is written.
    #[test]
    fn beside_an_input_and_an_output_is_refused_and_writes_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let output = run_beside_an_input_and_an_output(dir.path(), &["-V"]);
        assert_eq!(text(&output.stderr), REFUSAL);
        assert_eq!(text(&output.stdout), "");
        assert_eq!(output.status.code(), Some(1));
        assert!(!dir.path().join("out.txt").exists());
        assert_eq!(written(dir.path()), Vec::<String>::new());
    }

    /// `-v` in the place of `-V` is the verbose run the refusal names: exit 0, no
    /// refusal, and the output file holds the rows the same command writes without
    /// the flag.
    #[test]
    fn lower_case_short_in_its_place_annotates() {
        let verbose_dir = tempfile::tempdir().unwrap();
        let verbose = run_beside_an_input_and_an_output(verbose_dir.path(), &["-v"]);
        let stderr = text(&verbose.stderr);
        assert_eq!(verbose.status.code(), Some(0), "{stderr}");
        assert!(!stderr.contains("is not supported"), "{stderr}");
        assert_eq!(text(&verbose.stdout), "");

        let plain_dir = tempfile::tempdir().unwrap();
        let plain = run_beside_an_input_and_an_output(plain_dir.path(), &[]);
        assert_eq!(plain.status.code(), Some(0), "{}", text(&plain.stderr));

        let verbose_rows = rows(&verbose_dir.path().join("out.txt"));
        assert!(!verbose_rows.is_empty());
        assert_eq!(verbose_rows, rows(&plain_dir.path().join("out.txt")));
    }

    /// `--version` spelled out beside the same input and output wins over them:
    /// the version line on standard output, nothing on standard error, exit 0, and
    /// no file written.
    #[test]
    fn spelled_out_version_beside_an_input_and_an_output_prints_the_version() {
        let dir = tempfile::tempdir().unwrap();
        let output = run_beside_an_input_and_an_output(dir.path(), &["--version"]);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert_eq!(text(&output.stdout), version_line());
        assert_eq!(text(&output.stderr), "");
        assert_eq!(written(dir.path()), Vec::<String>::new());
    }
}
