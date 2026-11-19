// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `--custom` string through the `vep` binary: a malformed one, a remote
//! one, a transcript model and a file that cannot be opened each end the run
//! with one `ERROR: [custom] ...` message, exit status 1, before any output is
//! written. The texts are Ensembl VEP's where it has one
//! (`AnnotationSourceAdaptor.pm:327-460` at release/116.2), less its `ERROR: `
//! prefix. A `num_records` outside the whole numbers from 0 runs, with one
//! warnings-file line naming the value and the reading.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn corpus_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// Runs the binary on a one-record input under `dir` with `custom` as the
/// `--custom` value, writing `<dir>/out.txt`.
fn run_with_custom(dir: &Path, custom: &str, extra: &[&str]) -> Output {
    run_custom_on(dir, "21\t30282169\t.\tC\tT\t.\t.\t.", custom, extra)
}

/// [`run_with_custom`] on the VCF data line `record`.
fn run_custom_on(dir: &Path, record: &str, custom: &str, extra: &[&str]) -> Output {
    let input = dir.join("in.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(dir.join("out.txt"))
        .arg("--json_cache")
        .arg(corpus_cache())
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--quiet",
            "--custom",
            custom,
        ])
        .args(extra)
        .stdin(Stdio::null());
    cmd.output().unwrap()
}

/// Asserts the run exited 1 with `message` as its whole standard error, less
/// the trailing newline (`--quiet` keeps the banner off it), and wrote no
/// output file.
fn assert_refused(dir: &Path, custom: &str, extra: &[&str], message: &str) {
    let output = run_with_custom(dir, custom, extra);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{custom}: {stderr}");
    assert_eq!(stderr.trim_end_matches('\n'), message, "{custom}");
    assert!(
        !dir.join("out.txt").exists(),
        "{custom}: the output file was created"
    );
    assert!(
        !dir.join("out.txt_warnings.txt").exists(),
        "{custom}: the warnings file was created"
    );
}

/// An unknown key, a missing file or format, an invalid type, an unknown
/// format and an unsupported statistic each print Ensembl VEP's text under the
/// `custom` component, the string quoted on the `LINE:` line where Ensembl VEP
/// quotes it (`AnnotationSourceAdaptor.pm:371-377, 431-434`; `File.pm:181,
/// 260`), and `--show_cache_info` refuses the string too, as Ensembl VEP builds
/// its annotation sources before it prints (`Runner.pm:827`).
#[test]
fn malformed_custom_strings_exit_1_with_ensembl_vep_text_before_any_output() {
    let dir = tempfile::tempdir().unwrap();
    assert_refused(
        dir.path(),
        "file=a.vcf.gz,format=vcf,name=X",
        &[],
        "ERROR: [custom] The following options are not supported for custom annotations: name\n\
         LINE: --custom file=a.vcf.gz,format=vcf,name=X",
    );
    assert_refused(
        dir.path(),
        "format=vcf",
        &[],
        "ERROR: [custom] No file was added for custom annotation source.\nLINE: --custom format=vcf",
    );
    assert_refused(
        dir.path(),
        "a.vcf.gz,A",
        &[],
        "ERROR: [custom] No format specified for custom annotation source.\nLINE: --custom a.vcf.gz,A",
    );
    assert_refused(
        dir.path(),
        "a.vcf.gz,A,vcf,near",
        &[],
        "ERROR: [custom] New type \"near\" is not valid",
    );
    assert_refused(
        dir.path(),
        "a.tsv.gz,A,tsv",
        &[],
        "ERROR: [custom] Unknown or unsupported format tsv",
    );
    assert_refused(
        dir.path(),
        "file=a.bed.gz,format=bed,num_records=2.5,summary_stats=median",
        &[],
        "ERROR: [custom] The following summary statistics for custom annotations are not supported: median. \
         Available options are: min, max, mean, count, sum.\n\
         LINE: --custom file=a.bed.gz,format=bed,num_records=2.5,summary_stats=median",
    );
    // The grammar's last refusal (`BaseGXF.pm:109`) leaves no warnings file
    // either.
    assert_refused(
        dir.path(),
        "file=a.gff.gz,format=gff,gff_type=gencode_promoter,type=exact,num_records=-1",
        &[],
        "ERROR: [custom] GXF annotation sources cannot be used as \"exact\" custom annotation type",
    );
    let output = run_with_custom(dir.path(), "format=vcf", &["--show_cache_info"]);
    assert_eq!(output.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .starts_with("ERROR: [custom] No file was added for custom annotation source."),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty(), "the cache info was printed");
}

/// A remote file and a GFF or GTF transcript model are refused with a text of
/// vep-rs's own naming the route: a local copy, or the two-step cache recipe
/// the GTF and GFF3 page (`docs/cache-from-gtf-or-gff3.md`) states.
#[test]
fn remote_files_and_transcript_models_are_refused_with_their_route() {
    let dir = tempfile::tempdir().unwrap();
    assert_refused(
        dir.path(),
        "file=http://example.org/a.vcf.gz,format=vcf",
        &[],
        "ERROR: [custom] remote files are not read; copy http://example.org/a.vcf.gz to a local path",
    );
    for custom in ["file=a.gtf.gz,format=gtf", "a.gff.gz,A,gff"] {
        let output = run_with_custom(dir.path(), custom, &[]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{custom}: {stderr}");
        assert!(
            stderr.starts_with(
                "ERROR: [custom] a GFF or GTF source in transcript mode is not read; build a cache from the file and the genome, then annotate against it: vep-cache-builder --species <species> --assembly <assembly> --release <N> --gtf <file> --genome-fasta <genome.fa> --output-dir <dir> (--gff3 <file> for a GFF3), then vep -i <in> -o <out> --offline --json_cache <dir> --fasta <genome.fa>; the page docs/cache-from-gtf-or-gff3.md states the route"
            ),
            "{custom}: {stderr}"
        );
        assert_eq!(stderr.lines().count(), 1, "{custom}: {stderr}");
    }
    let page = fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/cache-from-gtf-or-gff3.md"),
    )
    .unwrap();
    for word in [
        "vep-cache-builder",
        "--gtf",
        "--gff3",
        "--genome-fasta",
        "--json_cache",
    ] {
        assert!(page.contains(word), "the page names {word}");
    }
}

/// A fractional or negative `num_records` is read as its integer part, at least
/// 0, and the warnings file carries one line naming the value and the reading;
/// Ensembl VEP 116.2 keeps the text and compares the record count with it as a
/// number (`AnnotationSource/File.pm:337`, `:436-439`), so its listing differs
/// there. The run otherwise proceeds: a deletion over three BED records lists
/// two and the `...` record under `2.5`, the `...` record alone under `-1`. A
/// `--show_cache_info` run opens the source as Ensembl VEP builds it before
/// the print (`Runner.pm:827`), so its warnings file carries the line too.
#[test]
fn a_fractional_or_negative_num_records_warns_once_and_lists_its_integer_part() {
    let regions = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-custom/custom/regions.bed.gz");
    let dir = tempfile::tempdir().unwrap();
    let custom = format!(
        "file={},short_name=regions,format=bed,num_records=2.5",
        regions.display()
    );
    let output = run_with_custom(dir.path(), &custom, &["--show_cache_info"]);
    assert_eq!(output.status.code(), Some(0), "{custom}");
    assert!(!output.stdout.is_empty(), "{custom}: no cache info printed");
    assert!(!dir.path().join("out.txt").exists(), "{custom}");
    assert_eq!(
        fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap(),
        format!(
            "WARNING: num_records=2.5 for {} is read as 2: \
             the limit is a whole number of at least 0\n",
            regions.display()
        ),
        "{custom}"
    );
    for (value, count, listed) in [("2.5", "2", "exact_1,overlap_1,..."), ("-1", "0", "...")] {
        let dir = tempfile::tempdir().unwrap();
        let custom = format!(
            "file={},short_name=regions,format=bed,num_records={value}",
            regions.display()
        );
        let output = run_custom_on(dir.path(), "21\t5028998\t.\tACGT\tA\t.\t.\t.", &custom, &[]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{custom}: {stderr}");
        let warnings = fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
        assert_eq!(
            warnings,
            format!(
                "WARNING: num_records={value} for {} is read as {count}: \
                 the limit is a whole number of at least 0\n",
                regions.display()
            ),
            "{custom}"
        );
        let out = fs::read_to_string(dir.path().join("out.txt")).unwrap();
        let rows: Vec<&str> = out.lines().filter(|line| !line.starts_with('#')).collect();
        assert!(!rows.is_empty(), "{custom}: no row\n{out}");
        for row in rows {
            assert!(
                row.contains(&format!("regions={listed}")),
                "{custom}: {row}"
            );
        }
    }
}

/// A source file that is not there, or has no index beside it, ends the run
/// before the input is read and the output is created; a `###CHR###` file
/// matching nothing prints the text of `AnnotationSourceAdaptor.pm:452`, less
/// its `Error: ` prefix. A fractional or negative `num_records` on such a
/// string leaves no warnings file either: its warning is written only once
/// the source opens.
#[test]
fn a_source_file_that_cannot_be_opened_exits_1_before_the_output_is_written() {
    let dir = tempfile::tempdir().unwrap();
    let absent = dir.path().join("absent.vcf.gz");
    let output = run_with_custom(
        dir.path(),
        &format!("file={},format=vcf,num_records=2.5", absent.display()),
        &[],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.starts_with(&format!(
            "ERROR: [custom] cannot open the --custom source {}: ",
            absent.display()
        )),
        "{stderr}"
    );
    assert!(!dir.path().join("out.txt").exists());
    assert!(!dir.path().join("out.txt_warnings.txt").exists());

    let pattern = dir.path().join("###CHR###.bed.gz");
    assert_refused(
        dir.path(),
        &format!(
            "file={},short_name=regions,format=bed,num_records=-1",
            pattern.display()
        ),
        &[],
        &format!(
            "ERROR: [custom] No files with pattern {} were found",
            pattern.display()
        ),
    );
}

/// A bgzipped VCF at `<dir>/<name>` with its tabix index beside it: the
/// `header` lines, the `#CHROM` line, then the one `record`.
fn indexed_vcf(dir: &Path, name: &str, header: &[&str], record: &str) -> PathBuf {
    use std::io::Write;

    use noodles::bgzf;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

    let path = dir.join(name);
    let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
    for line in header {
        writeln!(writer, "{line}").unwrap();
    }
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    let start = writer.virtual_position();
    writeln!(writer, "{record}").unwrap();
    let end = writer.virtual_position();
    writer.finish().unwrap();
    let fields: Vec<&str> = record.split('\t').collect();
    let pos: usize = fields[1].parse().unwrap();
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(HeaderBuilder::vcf().build());
    indexer
        .add_record(
            fields[0],
            Position::try_from(pos).unwrap(),
            Position::try_from(pos + fields[3].len() - 1).unwrap(),
            Chunk::new(start, end),
        )
        .unwrap();
    let mut index_path = path.as_os_str().to_owned();
    index_path.push(".tbi");
    noodles::tabix::fs::write(index_path, &indexer.build()).unwrap();
    path
}

/// Ensembl VEP opens a `--custom` VCF through its input parser, which ends
/// the run on a `##fileformat` line naming a major version other than 4
/// (`VCF4Tabix.pm:55-58`, `BaseVCF4.pm:125`) or a value with no
/// `<digits>.<digits>` (`:129`) before any output is written; a `VCFv4.5` line
/// and a file with no `##fileformat` line are read. Measured on 116.2 with
/// each shape as the first header line of a sites file: exit 2 and no output
/// for `VCFv3.3` and `VCFv4`, the annotated row for the other three.
#[test]
fn a_vcf_source_whose_fileformat_is_not_vcf_4_exits_1_before_the_output_is_written() {
    let dir = tempfile::tempdir().unwrap();
    let info = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">";
    let record = "21\t30282169\tprobe\tC\tT\t.\tPASS\tAF=0.5";
    let custom = |file: &Path| {
        format!(
            "file={},short_name=sites,format=vcf,type=exact,fields=AF",
            file.display()
        )
    };
    for (name, line, text) in [
        (
            "v33.vcf.gz",
            "##fileformat=VCFv3.3",
            "The VCF file format version 3.3 is not compatible with the parser version (VCF v4.2)",
        ),
        (
            "v4.vcf.gz",
            "##fileformat=VCFv4",
            "The script can't read the VCF file format version of 'fileformat'",
        ),
    ] {
        let file = indexed_vcf(dir.path(), name, &[line, info], record);
        assert_refused(
            dir.path(),
            &custom(&file),
            &[],
            &format!(
                "ERROR: [custom] cannot open the --custom source {}: {line}: {text}",
                file.display()
            ),
        );
    }
    for (name, line) in [
        ("v42.vcf.gz", Some("##fileformat=VCFv4.2")),
        ("v45.vcf.gz", Some("##fileformat=VCFv4.5")),
        ("none.vcf.gz", None),
    ] {
        let header: Vec<&str> = line.into_iter().chain([info]).collect();
        let file = indexed_vcf(dir.path(), name, &header, record);
        let output = run_with_custom(dir.path(), &custom(&file), &[]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{name}: {stderr}");
        let text = fs::read_to_string(dir.path().join("out.txt")).unwrap();
        assert!(
            text.lines()
                .any(|row| !row.starts_with('#') && row.contains("sites=probe;sites_AF=0.5")),
            "{name}: {text}"
        );
    }
}
