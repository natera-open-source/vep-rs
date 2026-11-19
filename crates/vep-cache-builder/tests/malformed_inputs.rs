// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Malformed inputs the builder refuses or reads whole rather than failing on:
//! a GFF3 whose exons mix ranked and unranked lines, a seqid that would place a
//! chromosome's directory outside the cache, a FASTA index naming more bases
//! than the file holds, a GFF3 from which no transcript is read, one that
//! holds a sequence section, one with records behind a bare `##FASTA`
//! directive or behind a line that begins with `>`, one whose last line
//! begins with `>`, one with a line longer than a line's limit behind a line
//! that begins with `>`, and one that holds a record inside its sequence
//! section.
//! Each fixture under
//! `tests/fixtures/malformed/` is the shape as a file; a failed build prints
//! one `ERROR: [cache-builder]` line and exits 1, and a refused build writes
//! nothing outside the output directory.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use serde_json::Value;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/malformed")
        .join(name)
}

/// Runs the builder on `gff3` and `peptides` into `out` with `extra` flags.
fn build(gff3: &Path, peptides: &Path, out: &Path, extra: &[(&str, &Path)]) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    cmd.args(["--release", "116", "--assembly", "GRCh37"])
        .arg("--gff3")
        .arg(gff3)
        .arg("--fasta")
        .arg(peptides)
        .arg("--output-dir")
        .arg(out);
    for (flag, path) in extra {
        cmd.arg(flag).arg(path);
    }
    cmd.output().expect("vep-cache-builder runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The names of every entry directly under `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// Exit 1 with `line` as the one `ERROR: [` line on standard error.
fn assert_refused(output: &Output, line: &str) {
    let err = stderr(output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(errors, [line], "{err}");
}

/// The chromosome 21 shard a build of `gff3` into `out` writes, the build
/// succeeding, and what the build logs.
fn shard_and_log(gff3: &Path, out: &Path) -> (Vec<u8>, String) {
    let output = build(gff3, &fixture("mixed_exon_ranks_peptides.fa"), out, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    let shard = std::fs::read(out.join("transcripts/21/1-1000000.json")).unwrap();
    (shard, stderr(&output))
}

/// `model` with `line` between its two halves, and that line's number.
fn with_a_line_mid_file(model: &str, line: &str) -> (String, usize) {
    let lines: Vec<&str> = model.lines().collect();
    let (before, behind) = lines.split_at(lines.len() / 2);
    let text = format!("{}\n{line}\n{}\n", before.join("\n"), behind.join("\n"));
    (text, before.len() + 1)
}

/// Sixty-four exons of one forward-strand transcript, 33 with a `rank`
/// attribute and 31 without, listed in no particular order: the build
/// completes, and the cache carries the exons in genomic order, numbered 1 to
/// 64 in that order, as the transcript order falls back to the start coordinate
/// when any exon lacks a rank and the file's ranks go with it.
#[test]
fn exons_mixing_ranked_and_unranked_lines_build_in_start_order() {
    let out = tempfile::tempdir().unwrap();
    let output = build(
        &fixture("mixed_exon_ranks.gff3"),
        &fixture("mixed_exon_ranks_peptides.fa"),
        out.path(),
        &[],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let shards = entries(&out.path().join("transcripts/21"));
    assert_eq!(shards, ["1-1000000.json"]);
    let text = std::fs::read_to_string(out.path().join("transcripts/21/1-1000000.json")).unwrap();
    let records: Vec<Value> = serde_json::from_str(&text).unwrap();
    assert_eq!(records.len(), 1);
    let exons = records[0]["exons"].as_array().unwrap();
    assert_eq!(exons.len(), 64);
    let starts: Vec<u64> = exons.iter().map(|e| e["start"].as_u64().unwrap()).collect();
    let mut sorted = starts.clone();
    sorted.sort_unstable();
    assert_eq!(starts, sorted);
    let ranks: Vec<u64> = exons.iter().map(|e| e["rank"].as_u64().unwrap()).collect();
    assert_eq!(ranks, (1..=64).collect::<Vec<u64>>());
}

/// A seqid of `../../escaped`, and one that is an absolute path, name a
/// directory outside the output: the build is refused at the GFF3, naming the
/// file, the line and the seqid, and the output directory's parent holds
/// nothing but the output directory, which is empty.
#[test]
fn a_seqid_escaping_the_output_directory_is_refused_and_writes_nothing() {
    let parent = tempfile::tempdir().unwrap();
    let out = parent.path().join("cache").join("deep");
    std::fs::create_dir_all(&out).unwrap();
    let output = build(
        &fixture("seqid_parent_escape.gff3"),
        &fixture("seqid_escape_peptides.fa"),
        &out,
        &[],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line 2: invalid chromosome name '../../escaped'",
            fixture("seqid_parent_escape.gff3").display()
        ),
    );
    assert_eq!(entries(parent.path()), ["cache"]);
    assert_eq!(entries(&parent.path().join("cache")), ["deep"]);
    assert!(entries(&out).is_empty(), "{:?}", entries(&out));

    let target = parent.path().join("elsewhere");
    let absolute = std::fs::read_to_string(fixture("seqid_parent_escape.gff3"))
        .unwrap()
        .replace("../../escaped", target.to_str().unwrap());
    let gff3 = parent.path().join("absolute.gff3");
    std::fs::write(&gff3, absolute).unwrap();
    let output = build(&gff3, &fixture("seqid_escape_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line 2: invalid chromosome name '{}'",
            gff3.display(),
            target.display()
        ),
    );
    assert_eq!(entries(parent.path()), ["absolute.gff3", "cache"]);
    assert!(entries(&out).is_empty(), "{:?}", entries(&out));
}

/// A `.fai` naming a 4e18-base chromosome 21 over a 126-byte FASTA of 120
/// bases, with a GFF3 whose transcript spans 3e18 bases: the genome FASTA is
/// refused as it opens, before any sequence is read, naming the bases the file
/// lacks (the 120 bases after the 4-byte header are held, the rest are not).
#[test]
fn a_fasta_index_naming_bases_beyond_the_file_is_refused() {
    let out = tempfile::tempdir().unwrap();
    let output = build(
        &fixture("index_overstates_contig.gff3"),
        &fixture("index_overstates_contig_peptides.fa"),
        out.path(),
        &[("--genome-fasta", &fixture("index_overstates_contig.fa"))],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to load genome FASTA (needs .fa + .fai index): FASTA index names 3999999999999999880 bases beyond the end of {}",
            fixture("index_overstates_contig.fa").display()
        ),
    );
    assert!(entries(out.path()).is_empty(), "{:?}", entries(out.path()));
}

/// A GFF3 from which no transcript is read, here one whose only line is of a
/// type the model skips, builds no cache: the build is refused before a file
/// or a directory is written, with one `ERROR: [cache-builder]` line naming
/// the file, and the output directory is not created.
#[test]
fn a_model_with_no_transcript_is_refused_and_writes_nothing() {
    let parent = tempfile::tempdir().unwrap();
    let gff3 = parent.path().join("no_transcript.gff3");
    std::fs::write(
        &gff3,
        "##gff-version 3\n21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=region:1\n",
    )
    .unwrap();
    let out = parent.path().join("cache");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] No transcript was read from {}",
            gff3.display()
        ),
    );
    assert_eq!(entries(parent.path()), ["no_transcript.gff3"]);
}

/// A GFF3 with a sequence section, its first header behind a `##FASTA`
/// directive or without one, builds the shard of the same file without the
/// section, and the log names the header's line, where the section opens,
/// so the lines after it go unread with a word. The sequence behind the
/// header opens the section whatever the header's columns: one, several,
/// numbers in the fourth and fifth, or nine that read as a record.
#[test]
fn a_sequence_section_is_named_in_the_log_and_not_read() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let (without, _) = shard_and_log(&fixture("mixed_exon_ranks.gff3"), &dir.path().join("model"));
    for (name, opening, header_line) in [
        ("with_directive", "##FASTA\n>21 placeholder\n", 2),
        ("without_directive", ">21 placeholder\n", 1),
        ("with_tabs", ">21\tplaceholder\tof\tfour columns\n", 1),
        ("with_a_number", ">21\tplaceholder\tof\t100\n", 1),
        ("with_two_numbers", ">21\tplaceholder\tof\t100\t200\n", 1),
        (
            "with_nine_columns",
            ">21\tplaceholder\tof\t100\t200\t.\t+\t.\tnine columns\n",
            1,
        ),
    ] {
        let gff3 = dir.path().join(format!("{name}.gff3"));
        std::fs::write(&gff3, format!("{model}{opening}NNNNNNNN\n")).unwrap();
        let (shard, log) = shard_and_log(&gff3, &dir.path().join(name));
        let named = format!(
            "{} holds a sequence section from line {}: nothing from there on is read",
            gff3.display(),
            model.lines().count() + header_line
        );
        assert!(log.contains(&named), "{log}");
        assert!(shard == without, "{name}: the shard differs");
    }
}

/// A line that begins with `>` with records behind it and no sequence
/// opens no section: one of four, five or eight columns that writes numbers
/// where a start and an end are read is skipped as any short line, one of
/// nine of a skipped type as any such record, and the build reads every
/// record, to the shard of the same file without that line.
#[test]
fn a_line_beginning_with_an_angle_bracket_before_records_opens_no_section() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let (without, _) = shard_and_log(&fixture("mixed_exon_ranks.gff3"), &dir.path().join("model"));
    for (name, line) in [
        ("four_columns", ">21\tplaceholder\tof\t100"),
        ("five_columns", ">21\tplaceholder\tof\t100\t200"),
        ("eight_columns", ">21\tplaceholder\tof\t100\t200\t.\t+\t."),
        (
            "nine_columns",
            ">21\tx\tbiological_region\t100x\t200\t.\t+\t.\tID=region:1",
        ),
    ] {
        let gff3 = dir.path().join(format!("{name}.gff3"));
        std::fs::write(&gff3, with_a_line_mid_file(&model, line).0).unwrap();
        let (shard, log) = shard_and_log(&gff3, &dir.path().join(name));
        assert!(!log.contains("holds a sequence section"), "{log}");
        assert!(shard == without, "{name}: the shard differs");
    }
}

/// A line that begins with `>` opens no section when the line straight
/// behind it is one the index build places, though it is no record, with
/// records behind the two: the `>` line of five columns with numbers or of
/// nine that read as a record, the line behind it of five columns with
/// numbers or of nine whose start only leads with a number. Both are
/// skipped as any such line, and the build reads every record, to the
/// shard of the same file without the two.
#[test]
fn a_placed_line_behind_a_line_beginning_with_an_angle_bracket_opens_no_section() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let (without, _) = shard_and_log(&fixture("mixed_exon_ranks.gff3"), &dir.path().join("model"));
    for (first_name, first) in [
        ("numeric", ">21\tplaceholder\tx\t100\t200"),
        (
            "record",
            ">21\tx\tbiological_region\t100\t200\t.\t+\t.\tID=region:1",
        ),
    ] {
        for (behind_name, behind) in [
            ("short_numeric", "21\tplaceholder\tx\t1\t2"),
            (
                "number_led",
                "21\tx\tbiological_region\t1x\t2\t.\t+\t.\tID=region:9",
            ),
        ] {
            let name = format!("{first_name}_then_{behind_name}");
            let gff3 = dir.path().join(format!("{name}.gff3"));
            let two = format!("{first}\n{behind}");
            std::fs::write(&gff3, with_a_line_mid_file(&model, &two).0).unwrap();
            let (shard, log) = shard_and_log(&gff3, &dir.path().join(&name));
            assert!(!log.contains("holds a sequence section"), "{log}");
            assert!(shard == without, "{name}: the shard differs");
        }
    }
}

/// A line longer than a line's limit straight behind a line that begins
/// with `>` is a sequence's only when no column ends in it, however far
/// past the limit the line runs. One whose first tab stands past the limit,
/// with more of the line behind it than one read of the rest takes, opens
/// no section and is refused for its length. One of a single column that
/// long opens the section, and the record behind it is refused by name at
/// its own line. One whose only tab is the byte read past the limit opens
/// the section too, to the shard of the file without the two lines, and so
/// does one of a single column that ends the file without a line end.
#[test]
fn a_long_line_behind_a_line_beginning_with_an_angle_bracket_is_a_sequences_without_a_column_end() {
    const LIMIT: usize = 16 << 20;
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let end = model.lines().count();
    let filler = "x".repeat(LIMIT);
    let more = "x".repeat(128 << 10);
    let peptides = fixture("mixed_exon_ranks_peptides.fa");

    let gff3 = dir.path().join("tab_past_the_limit.gff3");
    let numeric = ">21\tplaceholder\tx\t100\t200";
    std::fs::write(&gff3, format!("{model}{numeric}\n{filler}x\t{more}\n")).unwrap();
    let output = build(&gff3, &peptides, &dir.path().join("tab"), &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line {} is longer than the {LIMIT}-byte limit for one line",
            gff3.display(),
            end + 2
        ),
    );

    let gff3 = dir.path().join("one_column.gff3");
    let record = "21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=region:1";
    std::fs::write(
        &gff3,
        format!("{model}>21 placeholder\n{filler}x{more}\n{record}\n"),
    )
    .unwrap();
    let output = build(&gff3, &peptides, &dir.path().join("column"), &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line {}: GFF3 annotation after the sequence section, which opens at line {}",
            gff3.display(),
            end + 3,
            end + 1
        ),
    );

    let (without, _) = shard_and_log(&fixture("mixed_exon_ranks.gff3"), &dir.path().join("model"));
    let gff3 = dir.path().join("tab_at_the_limit.gff3");
    std::fs::write(&gff3, format!("{model}>21 placeholder\n{filler}\t\n")).unwrap();
    let (shard, log) = shard_and_log(&gff3, &dir.path().join("limit"));
    let named = format!(
        "{} holds a sequence section from line {}: nothing from there on is read",
        gff3.display(),
        end + 1
    );
    assert!(log.contains(&named), "{log}");
    assert!(shard == without, "the shard differs");

    let gff3 = dir.path().join("no_line_end.gff3");
    std::fs::write(&gff3, format!("{model}>21 placeholder\n{filler}x{more}")).unwrap();
    let (shard, log) = shard_and_log(&gff3, &dir.path().join("end"));
    let named = format!(
        "{} holds a sequence section from line {}: nothing from there on is read",
        gff3.display(),
        end + 1
    );
    assert!(log.contains(&named), "{log}");
    assert!(shard == without, "the shard differs");
}

/// A line that begins with `>`, has fewer than four columns and stands
/// before records is a record cut short, as the index build fails it: the
/// build is refused at that line and writes nothing.
#[test]
fn a_short_line_beginning_with_an_angle_bracket_before_records_is_cut_short() {
    let parent = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let (text, line) = with_a_line_mid_file(&model, ">21 placeholder");
    let gff3 = parent.path().join("short.gff3");
    std::fs::write(&gff3, text).unwrap();
    let out = parent.path().join("cache");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line {line}: GFF3 line ends at column 1, before its start",
            gff3.display()
        ),
    );
    assert!(!out.exists());
}

/// A comment that is not UTF-8 between a line that begins with `>` and a
/// record is read ahead of the `>` line and refused at its own line, as on
/// any line: the build writes nothing.
#[test]
fn a_comment_read_ahead_that_is_not_utf8_is_refused_at_its_line() {
    let parent = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    assert!(!model.contains('@'));
    let (text, line) = with_a_line_mid_file(&model, ">21\tplaceholder\tof\t100\n# caf@");
    let bytes: Vec<u8> = text
        .bytes()
        .map(|byte| if byte == b'@' { 0xe9 } else { byte })
        .collect();
    let gff3 = parent.path().join("comment.gff3");
    std::fs::write(&gff3, bytes).unwrap();
    let out = parent.path().join("cache");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line {}: GFF3 line is not UTF-8",
            gff3.display(),
            line + 1
        ),
    );
    assert!(!out.exists());
}

/// A `##FASTA` directive opens nothing by itself, so a sequence line behind
/// one with no `>` header before it is a record cut short: the build is
/// refused at the sequence line and writes nothing.
#[test]
fn a_sequence_line_behind_a_bare_fasta_directive_is_cut_short() {
    let parent = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let gff3 = parent.path().join("headerless.gff3");
    std::fs::write(&gff3, format!("{model}##FASTA\nNNNNNNNN\n")).unwrap();
    let out = parent.path().join("cache");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line {}: GFF3 line ends at column 1, before its start",
            gff3.display(),
            model.lines().count() + 2
        ),
    );
    assert!(!out.exists());
}

/// A line that begins with `>` on the file's last line has nothing behind
/// it and is read as any other line: one of fewer than four columns is a
/// record cut short, refused, and one the reader skips, four columns with a
/// number or nine of a skipped type, leaves the build the shard of the same
/// file without it, the log naming no section.
#[test]
fn a_line_beginning_with_an_angle_bracket_on_the_last_line_is_a_line_like_any_other() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let gff3 = dir.path().join("short.gff3");
    std::fs::write(&gff3, format!("{model}>21 placeholder\n")).unwrap();
    let out = dir.path().join("short");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: {} line {}: GFF3 line ends at column 1, before its start",
            gff3.display(),
            model.lines().count() + 1
        ),
    );
    assert!(!out.exists());
    let (without, _) = shard_and_log(&fixture("mixed_exon_ranks.gff3"), &dir.path().join("model"));
    for (name, line) in [
        ("four_columns", ">21\tplaceholder\tof\t100"),
        (
            "nine_columns",
            ">21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=region:1",
        ),
    ] {
        let gff3 = dir.path().join(format!("{name}.gff3"));
        std::fs::write(&gff3, format!("{model}{line}\n")).unwrap();
        let (shard, log) = shard_and_log(&gff3, &dir.path().join(name));
        assert!(!log.contains("holds a sequence section"), "{log}");
        assert!(shard == without, "{name}: the shard differs");
    }
}

/// A bare `##FASTA` directive with records behind it is a comment: the
/// build reads every record, to the shard of the same file without the
/// directive, and the log names no section.
#[test]
fn records_behind_a_bare_fasta_directive_are_read() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let lines: Vec<&str> = model.lines().collect();
    let (before, behind) = lines.split_at(lines.len() / 2);
    let gff3 = dir.path().join("directive.gff3");
    std::fs::write(
        &gff3,
        format!("{}\n##FASTA\n{}\n", before.join("\n"), behind.join("\n")),
    )
    .unwrap();
    let shard = |gff3: &Path, name: &str| {
        let out = dir.path().join(name);
        let output = build(gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
        assert!(output.status.success(), "{}", stderr(&output));
        assert!(
            !stderr(&output).contains("holds a sequence section"),
            "{}",
            stderr(&output)
        );
        std::fs::read(out.join("transcripts/21/1-1000000.json")).unwrap()
    };
    assert!(
        shard(&gff3, "with") == shard(&fixture("mixed_exon_ranks.gff3"), "without"),
        "the shards differ"
    );
}

/// A record whose seqid begins with `>` opens no sequence section when a
/// record follows it: placed before the model's lines, it leaves every one
/// of them read, and the log names no section.
#[test]
fn a_record_whose_seqid_begins_with_an_angle_bracket_opens_no_sequence_section() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let gff3 = dir.path().join("record.gff3");
    std::fs::write(
        &gff3,
        format!(">21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=region:1\n{model}"),
    )
    .unwrap();
    let out = dir.path().join("cache");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        !stderr(&output).contains("holds a sequence section"),
        "{}",
        stderr(&output)
    );
    assert_eq!(entries(&out.join("transcripts/21")), ["1-1000000.json"]);
}

/// Inside a sequence section a `>` line is a header whatever its columns,
/// and a comment stays a comment behind any white space: one of each that
/// writes nine columns with numbers in its fourth and fifth refuses
/// nothing, and the build reads the lines before the section.
#[test]
fn a_header_or_a_comment_inside_a_sequence_section_refuses_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let gff3 = dir.path().join("section.gff3");
    std::fs::write(
        &gff3,
        format!(
            "{model}>21 placeholder\nNNNNNNNN\n\
             >22\tplaceholder\tx\t100\t200\t.\t+\t.\tnine columns\nNNNNNNNN\n\
             \u{b}#21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=region:1\n"
        ),
    )
    .unwrap();
    let out = dir.path().join("cache");
    let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(entries(&out.join("transcripts/21")), ["1-1000000.json"]);
}

/// A record inside a sequence section is annotation the section has ended:
/// behind a `>` line and the sequence that makes it a header, with a
/// `##FASTA` directive before the two or without one, the build is refused
/// with one error naming the record's line and the header's line, where the
/// section opens, and writes nothing.
#[test]
fn annotation_inside_a_sequence_section_is_refused_and_writes_nothing() {
    let parent = tempfile::tempdir().unwrap();
    let model = std::fs::read_to_string(fixture("mixed_exon_ranks.gff3")).unwrap();
    let end = model.lines().count();
    for (name, opening, opened, record_line) in [
        ("implied", ">21 placeholder\nNNNNNNNN\n", end + 1, end + 3),
        (
            "stated",
            "##FASTA\n>21 placeholder\nNNNNNNNN\n",
            end + 2,
            end + 4,
        ),
    ] {
        let gff3 = parent.path().join(format!("{name}.gff3"));
        std::fs::write(
            &gff3,
            format!("{model}{opening}21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=region:1\n"),
        )
        .unwrap();
        let out = parent.path().join(name);
        let output = build(&gff3, &fixture("mixed_exon_ranks_peptides.fa"), &out, &[]);
        assert_refused(
            &output,
            &format!(
                "ERROR: [cache-builder] Failed to parse GFF3: {} line {record_line}: GFF3 annotation after the sequence section, which opens at line {opened}",
                gff3.display()
            ),
        );
        assert!(!out.exists());
    }
}
