// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The proof behind `docs/developers/embedding.md`: the re-exported surface of
//! `vep_cli` writes the bytes the `vep` binary writes. One `Annotator` built
//! from the binary's own command line annotates each record of a corpus, and
//! the `JsonOutputFormatter` renders it as `vep --json` does.

mod common;

use vep_cli::vcf_parser::{parse_vcf_line, VcfParse};
use vep_cli::{parse_vcf_line_with, Annotator, Config, InputVariant, JsonOutputFormatter};
use vep_io::output::fields::FieldOptions;

use common::{corpora, default_run, input_records, output_path, run_vep};

/// The JSON rows the library writes for the GRCh37 corpus equal, byte for byte,
/// the file `vep --json` writes for it from the same command line: the binary's
/// reader parses each data line, so the comparison proves the annotator and the
/// formatter. The binary writes no JSON object for a record whose span exceeds
/// `--max_sv_size`, so the caller leaves those records out, as the page says;
/// their number is the number of `--max_sv_size` warnings the binary wrote, and
/// every other record of the corpus has its row.
#[test]
fn annotator_json_rows_equal_the_vep_binarys_json_output_on_the_corpus() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh37")
        .expect("tests/golden/116/GRCh37");
    let run = default_run(corpus);
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);
    let warnings = std::fs::read_to_string(format!(
        "{}_warnings.txt",
        output_path(&run, "json", out_dir.path()).display()
    ))
    .unwrap();
    let oversize_warnings = warnings
        .lines()
        .filter(|l| l.contains("is bigger than --max_sv_size"))
        .count();

    // The binary's command line without its input and output files: the cache,
    // the species and the corpus flags, which name the assembly.
    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend(["--species".to_string(), "homo_sapiens".to_string()]);
    argv.extend(run.flags.iter().cloned());
    let config = Config::from_argv(argv).unwrap();
    let max_sv_size = config.max_sv_size.expect("the default limit");
    let annotator = Annotator::from_config(&config).unwrap();
    let formatter = JsonOutputFormatter::new(
        Default::default(),
        annotator.plugin_header_fields(),
        Some(corpus.assembly.clone()),
    );

    let records = input_records(corpus);
    let mut library_rows = String::new();
    let mut oversize = 0;
    for line in &records {
        let mut record = parse_vcf_line(line, true).unwrap();
        if record
            .iter()
            .any(|v| v.end.saturating_sub(v.start) > max_sv_size)
        {
            oversize += 1;
            continue;
        }
        annotator.annotate_batch(&mut record);
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let value = formatter.format_record_value(&alleles);
        library_rows.push_str(&serde_json::to_string(&value).unwrap());
        library_rows.push('\n');
    }
    assert_eq!(oversize, oversize_warnings);
    assert!(oversize > 0, "the corpus carries records past the limit");
    assert_eq!(library_rows.lines().count(), records.len() - oversize);

    let differing = library_rows
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(library_rows, binary_rows);
}

/// The same proof over the corpus whose run loads the LoFTEE built-in: the
/// library's `Annotator::from_config` loads the plugin from the command line
/// and attaches the exon and intron context LoFTEE reads to every transcript
/// before the plugin runs, so the `lof` and `lof_info` fields of its rows equal
/// the binary's, byte for byte, with the confidence spread the corpus carries.
#[test]
fn annotator_json_rows_equal_the_vep_binarys_json_output_with_the_loftee_plugin() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-plugin-loftee")
        .expect("tests/golden/116/GRCh38-plugin-loftee");
    let run = default_run(corpus);
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);

    let fasta = corpus
        .dir
        .join(corpus.manifest["fasta"].as_str().expect("a corpus FASTA"));
    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend([
        "--species".to_string(),
        "homo_sapiens".to_string(),
        "--assembly".to_string(),
        corpus.assembly.clone(),
    ]);
    argv.extend(run.flags.iter().cloned());
    argv.extend(["--fasta".to_string(), fasta.display().to_string()]);
    let config = Config::from_argv(argv).unwrap();
    let annotator = Annotator::from_config(&config).unwrap();
    let formatter = JsonOutputFormatter::new(
        Default::default(),
        annotator.plugin_header_fields(),
        Some(corpus.assembly.clone()),
    );

    let records = input_records(corpus);
    let mut library_rows = String::new();
    for line in &records {
        let mut record = parse_vcf_line(line, true).unwrap();
        annotator.annotate_batch(&mut record);
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let value = formatter.format_record_value(&alleles);
        library_rows.push_str(&serde_json::to_string(&value).unwrap());
        library_rows.push('\n');
    }
    assert_eq!(library_rows.lines().count(), records.len());
    for verdict in ["\"lof\":\"HC\"", "\"lof\":\"LC\""] {
        assert!(
            library_rows.contains(verdict),
            "the library wrote no {verdict} field"
        );
    }

    let differing = library_rows
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(library_rows, binary_rows);
}

/// The same proof over a run that reads two gene-model files beside the cache,
/// a `--gtf` and a `--custom` GFF3: the library joins each file's transcripts
/// to the cache's before it indexes a chromosome, the binary as it loads one,
/// and the rows of a record, the cache's and each file's under its `source`,
/// are the same bytes in the same order.
#[test]
fn annotator_json_rows_equal_the_vep_binarys_json_output_beside_gene_model_files() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-gene-model-beside-cache")
        .expect("tests/golden/116/GRCh38-gene-model-beside-cache");
    let file = |name: &str| corpus.dir.join(name).display().to_string();
    let mut run = default_run(corpus);
    run.flags = vec![
        "--gtf".to_string(),
        file("chr21_genes.gtf.gz"),
        "--custom".to_string(),
        format!(
            "file={},format=gff,short_name=novel",
            file("novel_models.gff3.gz")
        ),
    ];
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);

    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend([
        "--species".to_string(),
        "homo_sapiens".to_string(),
        "--assembly".to_string(),
        corpus.assembly.clone(),
    ]);
    argv.extend(run.flags.iter().cloned());
    argv.extend(["--fasta".to_string(), file("reference.fa.gz")]);
    let config = Config::from_argv(argv).unwrap();
    let annotator = Annotator::from_config(&config).unwrap();
    let formatter = JsonOutputFormatter::new(
        {
            let mut options = FieldOptions::default();
            options.custom_sources = annotator.resources().custom_source_infos.clone();
            options
        },
        annotator.plugin_header_fields(),
        Some(corpus.assembly.clone()),
    );

    let records = input_records(corpus);
    let mut library_rows = String::new();
    for line in &records {
        let mut record = parse_vcf_line(line, true).unwrap();
        annotator.annotate_batch(&mut record);
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let value = formatter.format_record_value(&alleles);
        library_rows.push_str(&serde_json::to_string(&value).unwrap());
        library_rows.push('\n');
    }
    assert_eq!(library_rows.lines().count(), records.len());
    for source in ["\"source\":\"chr21_genes.gtf.gz\"", "\"source\":\"novel\""] {
        assert!(
            library_rows.contains(source),
            "the library wrote no {source} row"
        );
    }

    let differing = library_rows
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(library_rows, binary_rows);
}

/// Installs, once for this test binary, a global default subscriber that
/// enables nothing and answers every callsite with `Interest::sometimes`, so
/// the subscriber of the emitting thread decides each event. `tracing` caches
/// a callsite's interest when a thread first reaches it, from that thread's
/// dispatcher while one scoped subscriber is registered: a test building an
/// annotator with no subscriber on its thread would cache the callsite as
/// never enabled, and [`info_messages_of`] on another thread would lose the
/// same line.
fn leave_callsite_interest_to_the_thread() {
    use tracing::level_filters::LevelFilter;
    use tracing::span::{Attributes, Id, Record};
    use tracing::subscriber::Interest;
    use tracing::{Event, Metadata, Subscriber};

    struct Sink;

    impl Subscriber for Sink {
        fn register_callsite(&self, _: &'static Metadata<'static>) -> Interest {
            Interest::sometimes()
        }

        fn max_level_hint(&self) -> Option<LevelFilter> {
            Some(LevelFilter::OFF)
        }

        fn enabled(&self, _: &Metadata<'_>) -> bool {
            false
        }

        fn new_span(&self, _: &Attributes<'_>) -> Id {
            Id::from_u64(1)
        }

        fn record(&self, _: &Id, _: &Record<'_>) {}

        fn record_follows_from(&self, _: &Id, _: &Id) {}

        fn event(&self, _: &Event<'_>) {}

        fn enter(&self, _: &Id) {}

        fn exit(&self, _: &Id) {}
    }

    static SINK: std::sync::Once = std::sync::Once::new();
    SINK.call_once(|| {
        let _ = tracing::subscriber::set_global_default(Sink);
    });
}

/// The messages of the `INFO` events `work` emits on this thread, in order.
fn info_messages_of(work: impl FnOnce()) -> Vec<String> {
    use std::sync::{Arc, Mutex};
    use tracing::field::{Field, Visit};
    use tracing_subscriber::layer::{Context, Layer, SubscriberExt};

    leave_callsite_interest_to_the_thread();

    struct Capture(Arc<Mutex<Vec<String>>>);
    struct Message(String);
    impl Visit for Message {
        fn record_debug(&mut self, field: &Field, value: &dyn std::fmt::Debug) {
            if field.name() == "message" {
                self.0 = format!("{value:?}");
            }
        }
    }
    impl<S: tracing::Subscriber> Layer<S> for Capture {
        fn on_event(&self, event: &tracing::Event<'_>, _ctx: Context<'_, S>) {
            if *event.metadata().level() != tracing::Level::INFO {
                return;
            }
            let mut message = Message(String::new());
            event.record(&mut message);
            self.0.lock().unwrap().push(message.0);
        }
    }
    let lines = Arc::new(Mutex::new(Vec::new()));
    let subscriber = tracing_subscriber::registry().with(Capture(Arc::clone(&lines)));
    tracing::subscriber::with_default(subscriber, work);
    let captured = lines.lock().unwrap().clone();
    captured
}

/// The status lines the library writes as it loads the cache and the two
/// gene-model files of the test above, whole. The files' transcripts join the
/// load before its lines are written, so the line of chromosome 21 counts the
/// cache's 188 transcripts and the 59 of each file together, and it is the
/// line the binary writes for that chromosome under `--verbose`.
#[test]
fn annotator_status_lines_count_gene_model_transcripts_as_the_vep_binarys_do() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-gene-model-beside-cache")
        .expect("tests/golden/116/GRCh38-gene-model-beside-cache");
    let file = |name: &str| corpus.dir.join(name).display().to_string();
    let (gtf, gff, fasta) = (
        file("chr21_genes.gtf.gz"),
        file("novel_models.gff3.gz"),
        file("reference.fa.gz"),
    );
    let mut args = vec!["--offline".to_string()];
    args.extend(corpus.cache_args());
    args.extend(
        [
            "--species",
            "homo_sapiens",
            "--assembly",
            &corpus.assembly,
            "--gtf",
            &gtf,
            "--custom",
            &format!("file={gff},format=gff,short_name=novel"),
            "--fasta",
            &fasta,
        ]
        .map(str::to_string),
    );

    let config = Config::from_argv(std::iter::once("vep".to_string()).chain(args.clone())).unwrap();
    let library_lines = info_messages_of(|| {
        Annotator::from_config(&config).unwrap();
    });
    let chromosome_line = "  chromosome 21: 306 transcripts";
    assert_eq!(
        library_lines,
        [
            format!("Using reference FASTA: {fasta}"),
            format!(
                "Loading transcripts from JSON cache: {}",
                corpus.cache_dir().display()
            ),
            format!("Gene-model source novel: {gff} (gff), 59 transcripts"),
            format!(
                "Skipped the feature types five_prime_utr (44), start_codon (26), \
                 three_prime_utr (25) of {gtf}: not read as a gene, transcript, exon or CDS"
            ),
            format!("Gene-model source chr21_genes.gtf.gz: {gtf} (gtf), 59 transcripts"),
            "Loaded transcripts for 1 chromosome".to_string(),
            chromosome_line.to_string(),
            "Building transcript index".to_string(),
        ]
    );

    let out_dir = tempfile::tempdir().unwrap();
    let output = common::isolated_vep(out_dir.path())
        .args(&args)
        .arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("-o")
        .arg(out_dir.path().join("out.txt"))
        .args(["--verbose", "--fork", "1", "--no_stats"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let binary_lines: Vec<&str> = stderr
        .lines()
        .filter_map(|line| line.split_once(" INFO vep_cli::transcript_index: "))
        .map(|(_, message)| message)
        .collect();
    assert_eq!(binary_lines, [chromosome_line], "{stderr}");
}

/// The corpus whose cache names the alignment file its models were edited
/// against: resolving it switches `--use_transcript_ref` and the `BAM_EDIT`
/// field on, on the annotator's copy of the configuration. The caller's
/// configuration keeps both off and the annotator returns both on. The field
/// options that carry the annotator's two values write the rows of
/// `vep --json` byte for byte, where the options that carry the caller's write
/// no `given_ref` and no `used_ref`.
#[test]
fn annotator_returns_the_flags_a_bam_edited_cache_settles() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-refseq")
        .expect("tests/golden/116/GRCh38-refseq");
    let run = default_run(corpus);
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);

    let fasta = corpus
        .dir
        .join(corpus.manifest["fasta"].as_str().expect("a corpus FASTA"));
    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend([
        "--species".to_string(),
        "homo_sapiens".to_string(),
        "--assembly".to_string(),
        corpus.assembly.clone(),
    ]);
    argv.extend(run.flags.iter().cloned());
    argv.extend(["--fasta".to_string(), fasta.display().to_string()]);
    let config = Config::from_argv(argv).unwrap();
    let max_sv_size = config.max_sv_size.expect("the default limit");
    assert!(!config.use_transcript_ref, "the command line names no flag");
    assert!(!config.bam_edited);

    let annotator = Annotator::from_config(&config).unwrap();
    assert!(annotator.use_transcript_ref());
    assert!(annotator.bam_edited());
    assert_eq!(
        annotator.use_transcript_ref(),
        annotator.resources().effects_config.use_transcript_ref
    );

    let records = input_records(corpus);
    let rows_under = |use_transcript_ref: bool, bam_edited: bool| {
        let formatter = JsonOutputFormatter::new(
            {
                let mut options = FieldOptions::default();
                options.refseq = config.refseq;
                options.use_transcript_ref = use_transcript_ref;
                options.bam_edited = bam_edited;
                options
            },
            annotator.plugin_header_fields(),
            Some(corpus.assembly.clone()),
        );
        let mut rows = String::new();
        for line in &records {
            let mut record = parse_vcf_line(line, true).unwrap();
            if record
                .iter()
                .any(|v| v.end.saturating_sub(v.start) > max_sv_size)
            {
                continue;
            }
            annotator.annotate_batch(&mut record);
            let alleles: Vec<&InputVariant> = record.iter().collect();
            let value = formatter.format_record_value(&alleles);
            rows.push_str(&serde_json::to_string(&value).unwrap());
            rows.push('\n');
        }
        rows
    };
    let settled = rows_under(annotator.use_transcript_ref(), annotator.bam_edited());
    let differing = settled
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(settled, binary_rows);

    let unsettled = rows_under(config.use_transcript_ref, config.bam_edited);
    for key in ["\"given_ref\":", "\"used_ref\":"] {
        assert!(settled.contains(key), "the settled options write no {key}");
        assert!(
            !unsettled.contains(key),
            "the caller's own flags write {key}"
        );
    }
}

/// A program that parses a VCF data line under its own settings builds them
/// with `VcfParse::default_run` and assigns the public fields. Unassigned, the
/// settings parse a line as `parse_vcf_line` does; with `minimal` assigned the
/// parser trims the record and marks it.
#[test]
fn a_callers_own_line_settings_start_from_the_default_run_and_reach_the_parser() {
    let line = "21\t25585733\t.\tCA\tCT\t.\t.\t.";
    let mut settings = VcfParse::default_run(false);
    let as_built = parse_vcf_line_with(line, &settings).unwrap();
    let by_default = parse_vcf_line(line, false).unwrap();
    assert_eq!(format!("{as_built:?}"), format!("{by_default:?}"));
    assert!(as_built.iter().all(|allele| !allele.minimised));

    settings.minimal = true;
    let minimal = parse_vcf_line_with(line, &settings).unwrap();
    assert!(!minimal.is_empty());
    assert!(minimal.iter().all(|allele| allele.minimised));
}
