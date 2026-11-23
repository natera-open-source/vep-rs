// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One line that holds a NUL byte, handed to both readers of a regulatory
//! GFF: the bake (`vep-cache-builder --regulatory-gff`, then a run with
//! `--regulatory` over the cache it writes) and the run-time reader (a run
//! with `--regulatory_gff` over the file and its index). Both hold the line
//! to the start and the end the index build of Ensembl VEP's release reads in
//! it, a column ending at a tab or at a NUL byte there
//! ([`vep_core::regulatory_gff::index_columns`]), so both refuse a line, both
//! pass over it, or both build the same rows from it.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::index::header;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::tabix;
use vep_core::regulatory_gff::{index_columns, unindexed_coordinate};

/// One gene of one single-exon transcript on chromosome 21.
const MODEL: &str = "##gff-version 3\n\
     21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=lncRNA\n\
     21\tx\tlnc_RNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA\n\
     21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T1\n";

/// One record inside each of the file's three lines: the record at the
/// sequence's first base, the line under test and the record behind it.
const RECORDS: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
     21\t40\tfirst\tA\tG\t.\t.\t.\n\
     21\t5045\tchanged\tA\tG\t.\t.\t.\n\
     21\t9040\tlast\tA\tG\t.\t.\t.\n";

/// What both readers make of the line under test.
#[derive(Debug)]
enum Reading {
    /// A feature under this identifier, between the file's other two.
    Built(&'static str),
    /// The file's other two features and none from the line.
    Skipped,
    /// The file is refused. `run` is the text behind the file's name on the
    /// line the run ends on; `bake` is in the bake's error, `FILE` standing
    /// for the file's name.
    Refused { run: String, bake: String },
}

/// A refusal for the start or the end the index build fails: both readers
/// name the column and its text, the run the record and the bake the line.
fn coordinate(column: &str, record: &str, value: &str) -> Reading {
    Reading::Refused {
        run: format!("the {column} of {record} of chromosome 21 is not a whole number: {value}"),
        bake: format!("the {column} of line 2 of FILE is not a whole number: {value}"),
    }
}

/// A refusal for a line of the chromosome with `columns` tab-ended columns:
/// the run names the column it lacks, the bake the count.
fn cut_short(line: &str, columns: usize) -> Reading {
    Reading::Refused {
        run: format!("a record of chromosome 21 has no attributes column: {line}"),
        bake: format!("line 2 of FILE has {columns} columns where a GFF record has 9"),
    }
}

/// Every NUL shape on a line of `feature_type`, `enhancer` (a type a
/// regulatory file yields) or `exon` (one it skips): the shape, the line and
/// the reading both readers give it.
fn shapes(feature_type: &'static str) -> Vec<(&'static str, String, Reading)> {
    let kept = feature_type == "enhancer";
    let sound =
        format!("21\tEnsembl\t{feature_type}\t5000\t5100\t.\t.\t.\tID=CHANGED;color=#faca00");
    let with = |old: &str, new: &str| {
        assert_eq!(sound.matches(old).count(), 1, "{old:?}");
        sound.replacen(old, new, 1)
    };
    let record = format!("{feature_type} record CHANGED");
    // A line of a skipped type is read no further than the start and the
    // end the index build reads in it.
    let built = |id| {
        if kept {
            Reading::Built(id)
        } else {
            Reading::Skipped
        }
    };
    let no_whole_number = |column, value| {
        if kept {
            coordinate(column, &record, value)
        } else {
            Reading::Skipped
        }
    };
    let (split_type, type_tail) = if kept {
        ("enh\0ancer", "ancer")
    } else {
        ("ex\0on", "on")
    };
    let start_is_the_type = || coordinate("start", &record, feature_type);
    let nul_ended = sound.replace('\t', "\0");
    let type_to_start = with(&format!("{feature_type}\t"), &format!("{feature_type}\0"));
    let end_to_score = with("5100\t", "5100\0");
    vec![
        ("no NUL byte", sound.clone(), built("CHANGED")),
        (
            "inside the ID value",
            with("ID=CHANGED", "ID=CHAN\0GED"),
            built("CHAN%00GED"),
        ),
        (
            "ending the ID value",
            with("ID=CHANGED;", "ID=CHANGED\0;"),
            built("CHANGED%00"),
        ),
        (
            "inside another attribute's value",
            with("#faca00", "#fac\0a00"),
            built("CHANGED"),
        ),
        (
            "as the line's last byte",
            format!("{sound}\0"),
            built("CHANGED"),
        ),
        (
            "inside the type",
            with(feature_type, split_type),
            coordinate("start", &format!("{split_type} record CHANGED"), type_tail),
        ),
        (
            "ending the type column",
            with(&format!("{feature_type}\t"), &format!("{feature_type}\0\t")),
            Reading::Skipped,
        ),
        (
            "ending the score",
            with("\t.\t.\t.\tID", "\t.\0\t.\t.\tID"),
            built("CHANGED"),
        ),
        (
            "ending the strand",
            with("\t.\t.\t.\tID", "\t.\t.\0\t.\tID"),
            built("CHANGED"),
        ),
        (
            "ending the phase",
            with("\t.\t.\t.\tID", "\t.\t.\t.\0\tID"),
            built("CHANGED"),
        ),
        (
            "ending the start",
            with("\t5000\t", "\t5000\0\t"),
            no_whole_number("start", "5000\0"),
        ),
        (
            "among the start's digits",
            with("\t5000\t", "\t50\x0000\t"),
            no_whole_number("start", "50\x0000"),
        ),
        (
            "leading the start",
            with("\t5000\t", "\t\x005000\t"),
            coordinate("start", &record, ""),
        ),
        (
            "ending the end",
            with("\t5100\t", "\t5100\0\t"),
            no_whole_number("end", "5100\0"),
        ),
        (
            "among the end's digits",
            with("\t5100\t", "\t51\x0000\t"),
            no_whole_number("end", "51\x0000"),
        ),
        (
            "leading the end",
            with("\t5100\t", "\t\x005100\t"),
            coordinate("end", &record, ""),
        ),
        (
            "ending the sequence name",
            with("21\tEnsembl", "21\0\tEnsembl"),
            start_is_the_type(),
        ),
        (
            "inside the sequence name",
            with("21\tEnsembl", "2\x001\tEnsembl"),
            start_is_the_type(),
        ),
        (
            "ending the source",
            with("Ensembl\t", "Ensembl\0\t"),
            start_is_the_type(),
        ),
        (
            "inside the source",
            with("Ensembl", "Ens\0embl"),
            start_is_the_type(),
        ),
        ("in place of every tab", nul_ended.clone(), Reading::Skipped),
        (
            "in place of every tab, the start no number",
            nul_ended.replacen("5000", "abc", 1),
            Reading::Refused {
                // The line has no third tab-ended column for a type.
                run: "the start of a record of type  of chromosome 21 is not a whole number: abc"
                    .to_string(),
                bake: "line 2 of FILE has 1 columns where a GFF record has 9".to_string(),
            },
        ),
        (
            "in place of the tab between type and start",
            type_to_start.clone(),
            cut_short(&type_to_start, 8),
        ),
        (
            "in place of the tab between end and score",
            end_to_score.clone(),
            cut_short(&end_to_score, 8),
        ),
        (
            "three ending the sequence name, the source a number",
            with("21\tEnsembl", "21\0\0\0\t5050"),
            Reading::Skipped,
        ),
    ]
}

/// Writes `lines` as a bgzip-compressed file under `dir` with its `.tbi`
/// beside it, each line filed as the index build of Ensembl VEP's release
/// files it: by the columns that build reads, and a line whose start or end
/// it fails at the first base of the sequence of the line before it (htslib
/// 1.9 `tbx.c:235-236`).
fn write_indexed_gff(dir: &Path, name: &str, lines: &[&str]) -> PathBuf {
    let path = dir.join(name);
    let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(header::Builder::gff().build());
    let mut sequence = "";
    for line in lines {
        let start = writer.virtual_position();
        writer.write_all(line.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
        let end = writer.virtual_position();
        let columns: Vec<&str> = index_columns(line.as_bytes())
            .map(|column| std::str::from_utf8(column).unwrap_or(""))
            .collect();
        let position = |column: usize| {
            let base = columns.get(column).and_then(|text| text.parse().ok());
            base.and_then(|base: usize| Position::try_from(base).ok())
                .unwrap_or(Position::MIN)
        };
        let (first, second) = (position(3), position(4));
        if sequence.is_empty() || unindexed_coordinate(line.as_bytes()).is_none() {
            sequence = columns[0];
        }
        indexer
            .add_record(
                sequence,
                first.min(second),
                first.max(second),
                Chunk::new(start, end),
            )
            .unwrap();
    }
    writer.finish().unwrap();
    tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();
    path
}

/// Builds the model into `out`, with `regulatory` as the bake flag and its
/// file when given; the error line of a build that fails.
fn build(dir: &Path, out: &Path, regulatory: Option<(&str, &Path)>) -> Result<(), String> {
    let model = dir.join("model.gff3");
    std::fs::write(&model, MODEL).unwrap();
    let peptides = dir.join("peptides.fa");
    std::fs::write(&peptides, ">ENSP1 pep transcript:T1 gene:G1\nMK\n").unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    command
        .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
        .arg(&model)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(out);
    if let Some((flag, file)) = regulatory {
        command.arg(flag).arg(file);
    }
    let output = command.output().expect("vep-cache-builder runs");
    if output.status.success() {
        return Ok(());
    }
    let log = String::from_utf8_lossy(&output.stderr);
    let errors: Vec<&str> = log.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(errors.len(), 1, "{log}");
    Err(errors[0].to_string())
}

/// The regulatory and motif rows a run over `cache` under `flags` writes
/// for the three records, each row whole; the line a run that fails ends on.
fn annotate(dir: &Path, cache: &Path, flags: &[&str]) -> Result<Vec<String>, String> {
    let input = dir.join("input.vcf");
    std::fs::write(&input, RECORDS).unwrap();
    let out = dir.join("out.txt");
    let mut argv = vec![
        "vep".to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        cache.display().to_string(),
        "--no_stats".to_string(),
        "--quiet".to_string(),
        "--force_overwrite".to_string(),
        "-i".to_string(),
        input.display().to_string(),
        "-o".to_string(),
        out.display().to_string(),
    ];
    argv.extend(flags.iter().map(|flag| flag.to_string()));
    let config = vep_cli::Config::from_argv(argv).expect("the command line is accepted");
    vep_cli::runner::Runner::new(config)
        .and_then(|runner| runner.run())
        .map_err(|error| vep_cli::config::CliError::from_run(&error))?;
    let text = std::fs::read_to_string(&out).unwrap();
    Ok(text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| {
            matches!(
                line.split('\t').nth(5),
                Some("RegulatoryFeature" | "MotifFeature")
            )
        })
        .map(str::to_string)
        .collect())
}

/// Each row's record and feature identifier.
fn features(rows: &[String]) -> Vec<(&str, &str)> {
    rows.iter()
        .map(|row| {
            let mut columns = row.split('\t');
            (columns.next().unwrap(), columns.nth(3).unwrap())
        })
        .collect()
}

/// Hands `lines` as one file to both readers and holds both to `expected`:
/// `key` names the file to the run-time reader and `bake_flag` to the bake.
fn read_by_both(
    dir: &Path,
    transcripts: &Path,
    key: &str,
    bake_flag: &str,
    lines: &[&str],
    expected: &Reading,
    shape: &str,
) {
    let file = write_indexed_gff(dir, "features.gff3.gz", lines);
    let name = file.display().to_string();
    let run = annotate(
        dir,
        transcripts,
        &["--regulatory_gff", &format!("{key}={name}")],
    );
    let baked = dir.join("baked");
    if baked.exists() {
        std::fs::remove_dir_all(&baked).unwrap();
    }
    let bake = build(dir, &baked, Some((bake_flag, file.as_path())))
        .and_then(|()| annotate(dir, &baked, &["--regulatory"]));
    match (expected, &run, &bake) {
        (Reading::Refused { run, bake }, Err(run_error), Err(bake_error)) => {
            assert_eq!(
                *run_error,
                format!("ERROR: [cache] {name}: {run}"),
                "{shape}"
            );
            let bake = bake.replace("FILE", &name);
            assert!(bake_error.contains(&bake), "{shape}: {bake_error}");
        }
        (Reading::Built(id), Ok(run_rows), Ok(bake_rows)) => {
            assert_eq!(
                features(run_rows),
                [("first", "FIRST"), ("changed", *id), ("last", "LAST")],
                "{shape}"
            );
            assert_eq!(run_rows, bake_rows, "{shape}");
        }
        (Reading::Skipped, Ok(run_rows), Ok(bake_rows)) => {
            assert_eq!(
                features(run_rows),
                [("first", "FIRST"), ("last", "LAST")],
                "{shape}"
            );
            assert_eq!(run_rows, bake_rows, "{shape}");
        }
        _ => panic!("{shape}: expected {expected:?}, the run gave {run:?} and the bake {bake:?}"),
    }
}

/// Every NUL shape, on a line of a type a regulatory file yields and on one
/// of a type it skips, gets one reading from both readers. The line stands
/// behind a record at its sequence's first base, where the index build of
/// Ensembl VEP's release files a line it cannot parse. A NUL byte behind
/// the end column, or ending the type, moves no column a start or an end is
/// read in: the feature is built with the byte in its identifier, or the
/// line skipped by its type. One in the sequence name or the source, or
/// inside the type, puts the type or its tail where that build reads a
/// start, and one that leads a start or an end leaves it an empty column:
/// both readers refuse the file, naming the column and that text. A start
/// or an end with the byte behind or among its digits is one that build
/// parses and no whole number, refused on a record and unread on a skipped
/// line. A line whose every column a NUL byte ends, and one whose sequence
/// name ends in three before a numeric source, names no sequence of the
/// index where that build parses its start and end, and both pass over it.
/// A motif's strand that ends in a NUL byte is neither `+` nor `-` to both.
#[test]
fn both_readers_of_a_regulatory_gff_agree_on_every_line_that_holds_a_nul_byte() {
    let dir = tempfile::tempdir().unwrap();
    let transcripts = dir.path().join("transcripts-alone");
    build(dir.path(), &transcripts, None).unwrap();
    let first = "21\tEnsembl\tenhancer\t1\t100\t.\t.\t.\tID=FIRST";
    let last = "21\tEnsembl\tenhancer\t9000\t9100\t.\t.\t.\tID=LAST";
    for feature_type in ["enhancer", "exon"] {
        for (shape, line, expected) in shapes(feature_type) {
            assert_eq!(line.contains('\0'), shape != "no NUL byte", "{shape}");
            read_by_both(
                dir.path(),
                &transcripts,
                "file",
                "--regulatory-gff",
                &[first, &line, last],
                &expected,
                &format!("{feature_type}, a NUL byte {shape}"),
            );
        }
    }
    let motif = |start: u32, strand: &str, id: &str| {
        format!(
            "21\tEnsembl\tTF_binding_site\t{start}\t{}\t.\t{strand}\t.\tID={id};\
             binding_matrix_id=M{start};transcription_factor=TF{start}",
            start + 13
        )
    };
    read_by_both(
        dir.path(),
        &transcripts,
        "motifs",
        "--motif-gff",
        &[
            &motif(30, "-", "FIRST"),
            &motif(5040, "-\0", "CHANGED"),
            &motif(9030, "+", "LAST"),
        ],
        &Reading::Built("CHANGED"),
        "a motif, a NUL byte ending the strand",
    );
}
