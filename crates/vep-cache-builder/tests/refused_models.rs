// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript models the builder refuses whole, each with one
//! `ERROR: [cache-builder]` line, exit 1 and no output directory: a transcript
//! whose translation starts or ends in no exon (the CDS line that sets the
//! end starts before the exon it overlaps, ends past it or lies in none), in
//! a GFF3 and in a GTF, on either strand; one whose exon ranks run against
//! its strand; one whose translation would end at an offset outside the exon
//! it is counted in, its exon ranks running against position there or an exon
//! nesting in its last coding one; and a model whose every transcript a build
//! of the cache's content sets aside. Beside them, the models whose CDS lines
//! leave the exons between the translation's two ends: they build.

use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

/// A protein FASTA naming transcript `T1`; reading one makes the build the
/// cache's content.
const PEPTIDES: &str = ">ENSP1 pep transcript:T1 gene:G1\nMK\n";

/// Writes `text` as the model file `name` and builds it under `flag` (`--gff3`
/// or `--gtf`) with the protein FASTA: the build exits 1 with `expected` as
/// its one error line after `ERROR: [cache-builder] `, `{model}` standing for
/// the model's path, and makes no output directory.
fn assert_model_refused(name: &str, flag: &str, text: &str, expected: &str) {
    let parent = tempfile::tempdir().unwrap();
    let model = parent.path().join(name);
    std::fs::write(&model, text).unwrap();
    let peptides = parent.path().join("a_peptides.fa");
    std::fs::write(&peptides, PEPTIDES).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", flag])
        .arg(&model)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(parent.path().join("cache"))
        .output()
        .expect("vep-cache-builder runs");

    let err = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{name}: {err}");
    let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    let expected = format!(
        "ERROR: [cache-builder] {}",
        expected.replace("{model}", &model.display().to_string())
    );
    assert_eq!(errors, [expected], "{name}: {err}");
    let mut left: Vec<String> = std::fs::read_dir(parent.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    left.sort();
    assert_eq!(left, ["a_peptides.fa", name]);
}

/// An exon line's span and, in a GFF3, its `rank`.
type Exon = (u64, u64, Option<u32>);

/// The two exons of every model whose exon lines a test leaves alone.
const TWO_EXONS: [Exon; 2] = [(1000, 1200, None), (1500, 2000, None)];

/// The lines of a GFF3 gene `G1` with one transcript `T1` on `strand`: an
/// exon line per entry of `exons` and a CDS line per span of `cds`.
fn gff3_of(strand: char, exons: &[Exon], cds: &[(u64, u64)]) -> String {
    let mut text = format!(
        "##gff-version 3\n\
         21\tx\tgene\t1000\t2000\t.\t{strand}\t.\tID=gene:G1;biotype=protein_coding\n\
         21\tx\tmRNA\t1000\t2000\t.\t{strand}\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n"
    );
    for (start, end, rank) in exons {
        let rank = rank.map_or(String::new(), |rank| format!(";rank={rank}"));
        text +=
            &format!("21\tx\texon\t{start}\t{end}\t.\t{strand}\t.\tParent=transcript:T1{rank}\n");
    }
    for (start, end) in cds {
        text += &format!(
            "21\tx\tCDS\t{start}\t{end}\t.\t{strand}\t0\tParent=transcript:T1;protein_id=ENSP1\n"
        );
    }
    text
}

/// The same gene, transcript, exons and CDS lines as a GTF, whose exon lines
/// carry no rank.
fn gtf_of(strand: char, exons: &[Exon], cds: &[(u64, u64)]) -> String {
    let ids = "gene_id \"G1\"; transcript_id \"T1\";";
    let mut text = format!(
        "21\tx\tgene\t1000\t2000\t.\t{strand}\t.\tgene_id \"G1\"; gene_biotype \"protein_coding\";\n\
         21\tx\ttranscript\t1000\t2000\t.\t{strand}\t.\t{ids} transcript_biotype \"protein_coding\";\n"
    );
    for (start, end, _) in exons {
        text += &format!("21\tx\texon\t{start}\t{end}\t.\t{strand}\t.\t{ids}\n");
    }
    for (start, end) in cds {
        text +=
            &format!("21\tx\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{ids} protein_id \"ENSP1\";\n");
    }
    text
}

/// The GFF3 model of exons 1000-1200 and 1500-2000 with a CDS line per span
/// of `cds`.
fn gff3(strand: char, cds: &[(u64, u64)]) -> String {
    gff3_of(strand, &TWO_EXONS, cds)
}

/// The same model as a GTF.
fn gtf(strand: char, cds: &[(u64, u64)]) -> String {
    gtf_of(strand, &TWO_EXONS, cds)
}

/// The translation's start and its end each lie inside an exon. A model's
/// last CDS line that sets one of them in the intron refuses the build before
/// a file or a directory is written, as a GFF3 and as a GTF: it starts before
/// the exon it overlaps, ends past it, or is the coding sequence's one line
/// and lies in no exon. Its lowest base is the translation's start on the
/// forward strand and its end on the reverse, its highest base the other.
/// The line names the file, that CDS line, its span, the transcript and the
/// base.
#[test]
fn a_translation_that_starts_or_ends_in_no_exon_is_refused_and_writes_nothing() {
    for (shape, cds, forward, reverse) in [
        (
            "starts_before",
            [(1700, 1800), (1400, 1600)].as_slice(),
            "CDS 1400-1600 of transcript T1 starts the translation at 1400",
            "CDS 1400-1600 of transcript T1 ends the translation at 1400",
        ),
        (
            "ends_past",
            &[(1050, 1100), (1150, 1300)],
            "CDS 1150-1300 of transcript T1 ends the translation at 1300",
            "CDS 1150-1300 of transcript T1 starts the translation at 1300",
        ),
        (
            "in_no_exon",
            &[(1300, 1400)],
            "CDS 1300-1400 of transcript T1 starts the translation at 1300",
            "CDS 1300-1400 of transcript T1 starts the translation at 1400",
        ),
    ] {
        for (strand, name, refusal) in [('+', "forward", forward), ('-', "reverse", reverse)] {
            for (flag, extension, text, lines_before) in [
                ("--gff3", "gff3", gff3(strand, cds), 5),
                ("--gtf", "gtf", gtf(strand, cds), 4),
            ] {
                assert_model_refused(
                    &format!("{shape}_{name}.{extension}"),
                    flag,
                    &text,
                    &format!(
                        "Failed to build cache: {{model}} line {}: {refusal}, in no exon of the \
                         transcript",
                        lines_before + cds.len()
                    ),
                );
            }
        }
    }
}

/// The files a build of `text` writes, each by its path under the output
/// directory: the model is `model.<extension>` and the build reads the protein
/// FASTA.
fn built(flag: &str, extension: &str, text: &str) -> BTreeMap<String, Vec<u8>> {
    let parent = tempfile::tempdir().unwrap();
    let model = parent.path().join(format!("model.{extension}"));
    std::fs::write(&model, text).unwrap();
    let peptides = parent.path().join("peptides.fa");
    std::fs::write(&peptides, PEPTIDES).unwrap();
    let cache = parent.path().join("cache");
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", flag])
        .arg(&model)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(&cache)
        .output()
        .expect("vep-cache-builder runs");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    fn collect(dir: &Path, root: &Path, files: &mut BTreeMap<String, Vec<u8>>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                collect(&path, root, files);
            } else {
                let name = path.strip_prefix(root).unwrap().display().to_string();
                files.insert(name, std::fs::read(&path).unwrap());
            }
        }
    }
    let mut files = BTreeMap::new();
    collect(&cache, &cache, &mut files);
    files
}

/// Between the translation's two ends a CDS line is not held to the exons:
/// one more CDS line in the intron beside the sound lines, and one CDS line
/// across the intron in place of the two it spans, each build the cache the
/// sound model builds, byte for byte, as a GFF3 and as a GTF on either
/// strand.
#[test]
fn cds_lines_off_the_exons_between_the_translations_ends_build_the_sound_cache() {
    let sound = [(1100, 1200), (1500, 1600)];
    let one_more_in_the_intron = [(1100, 1200), (1500, 1600), (1300, 1400)];
    let one_across_the_intron = [(1100, 1600)];
    for strand in ['+', '-'] {
        let from_gff3 = built("--gff3", "gff3", &gff3(strand, &sound));
        let from_gtf = built("--gtf", "gtf", &gtf(strand, &sound));
        for sound in [&from_gff3, &from_gtf] {
            assert!(
                sound.keys().any(|name| name.starts_with("transcripts")),
                "{:?}",
                sound.keys()
            );
        }
        for cds in [
            one_more_in_the_intron.as_slice(),
            one_across_the_intron.as_slice(),
        ] {
            assert!(
                built("--gff3", "gff3", &gff3(strand, cds)) == from_gff3,
                "{strand} {cds:?} as a GFF3"
            );
            assert!(
                built("--gtf", "gtf", &gtf(strand, cds)) == from_gtf,
                "{strand} {cds:?} as a GTF"
            );
        }
    }
}

/// Exon ranks that run against the strand leave every CDS line inside an
/// exon and the coding sequence opening before the exon ranked first: the
/// build is refused, naming the line that opens the coding sequence, the
/// transcript, the whole coding span and that exon.
#[test]
fn exon_ranks_that_run_against_the_strand_are_refused_and_write_nothing() {
    assert_model_refused(
        "ranks.gff3",
        "--gff3",
        "##gff-version 3\n\
         21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
         21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
         21\tx\texon\t1000\t1200\t.\t+\t.\tParent=transcript:T1;rank=2\n\
         21\tx\texon\t1500\t2000\t.\t+\t.\tParent=transcript:T1;rank=1\n\
         21\tx\tCDS\t1100\t1200\t.\t+\t0\tParent=transcript:T1;protein_id=ENSP1\n\
         21\tx\tCDS\t1500\t1600\t.\t+\t0\tParent=transcript:T1;protein_id=ENSP1\n",
        "Failed to build cache: {model} line 6: the coding sequence of transcript T1, \
         1100-1600, starts before 1500-2000, the first exon in transcript order that it \
         overlaps",
    );
}

/// Neither end of the translation is written at an offset outside the exon
/// it is counted in. Exon ranks that run against position at the
/// translation's end, and an exon nested in the last coding exon and coming
/// after it, each leave every CDS line inside an exon and the last exon in
/// transcript order short of the coding sequence's last base: the build is
/// refused before a file or a directory is written, naming the CDS line that
/// sets the end, the offset, that exon and its length. A GTF's exon lines
/// carry no rank and come in position order, so the ranked lines build as a
/// GTF, and the nested exon comes last on the forward strand alone; on the
/// reverse it comes first, where the coding sequence opens outside it.
#[test]
fn a_translation_offset_outside_its_exon_is_refused_and_writes_nothing() {
    let coding = [(1100, 1200), (1400, 1500), (1700, 1800)];
    let end_ranked_second: [Exon; 3] = [
        (1000, 1200, Some(1)),
        (1400, 1500, Some(3)),
        (1700, 2000, Some(2)),
    ];
    assert_model_refused(
        "ranks_forward.gff3",
        "--gff3",
        &gff3_of('+', &end_ranked_second, &coding),
        "Failed to build cache: {model} line 9: CDS 1700-1800 of transcript T1 ends the \
         translation at base 401 of the 101-base exon 1400-1500, the last in transcript order \
         that the coding sequence overlaps",
    );
    let start_ranked_first: [Exon; 3] = [
        (1000, 1200, Some(2)),
        (1400, 1500, Some(3)),
        (1700, 2000, Some(1)),
    ];
    assert_model_refused(
        "ranks_reverse.gff3",
        "--gff3",
        &gff3_of('-', &start_ranked_first, &coding),
        "Failed to build cache: {model} line 7: CDS 1100-1200 of transcript T1 ends the \
         translation at base 401 of the 101-base exon 1400-1500, the last in transcript order \
         that the coding sequence overlaps",
    );
    for strand in ['+', '-'] {
        built("--gtf", "gtf", &gtf_of(strand, &end_ranked_second, &coding));
    }

    let across = [(1100, 1450)];
    let nested: [Exon; 2] = [(1000, 1500, None), (1300, 1350, None)];
    for (flag, extension, text, line) in [
        ("--gff3", "gff3", gff3_of('+', &nested, &across), 6),
        ("--gtf", "gtf", gtf_of('+', &nested, &across), 5),
    ] {
        assert_model_refused(
            &format!("nested_forward.{extension}"),
            flag,
            &text,
            &format!(
                "Failed to build cache: {{model}} line {line}: CDS 1100-1450 of transcript T1 \
                 ends the translation at base 151 of the 51-base exon 1300-1350, the last in \
                 transcript order that the coding sequence overlaps"
            ),
        );
    }
    let nested_and_ranked: [Exon; 2] = [(1000, 1500, Some(1)), (1300, 1350, Some(2))];
    assert_model_refused(
        "nested_reverse.gff3",
        "--gff3",
        &gff3_of('-', &nested_and_ranked, &across),
        "Failed to build cache: {model} line 6: CDS 1100-1450 of transcript T1 ends the \
         translation at base 251 of the 51-base exon 1300-1350, the last in transcript order \
         that the coding sequence overlaps",
    );
    assert_model_refused(
        "nested_reverse.gtf",
        "--gtf",
        &gtf_of('-', &nested, &across),
        "Failed to build cache: {model} line 5: the coding sequence of transcript T1, \
         1100-1450, ends after 1300-1350, the first exon in transcript order that it overlaps",
    );
}

/// Of several transcripts a build refuses, the one named is the same on
/// every run, the first of the gene whose identifier sorts lowest: six genes
/// listed from `G6` down to `G1`, each holding such a transcript, name `T1`
/// at the file's last line in each of eight builds.
#[test]
fn the_same_refused_transcript_is_named_on_every_run() {
    let mut model = String::from("##gff-version 3\n");
    for gene in (1..=6).rev() {
        let at = gene * 10_000;
        model += &format!(
            "21\tx\tgene\t{0}\t{1}\t.\t+\t.\tID=gene:G{gene};biotype=protein_coding\n\
             21\tx\tmRNA\t{0}\t{1}\t.\t+\t.\tID=transcript:T{gene};Parent=gene:G{gene};biotype=protein_coding\n\
             21\tx\texon\t{2}\t{3}\t.\t+\t.\tParent=transcript:T{gene}\n\
             21\tx\tCDS\t{4}\t{5}\t.\t+\t0\tParent=transcript:T{gene};protein_id=ENSP{gene}\n",
            at + 1000,
            at + 2000,
            at + 1100,
            at + 1200,
            at + 1050,
            at + 1150
        );
    }
    for _ in 0..8 {
        assert_model_refused(
            "six_genes.gff3",
            "--gff3",
            &model,
            "Failed to build cache: {model} line 25: CDS 11050-11150 of transcript T1 starts \
             the translation at 11050, in no exon of the transcript",
        );
    }
}

/// A build that reads a protein FASTA leaves out a readthrough transcript and
/// one of biotype `artifact`, so a model holding no other transcript leaves
/// it nothing to write: the build is refused before a file or a directory is
/// written, naming the file and counting the transcripts read and those set
/// aside under each reason. The GTF's one transcript is tagged readthrough;
/// the GFF3 holds an `artifact` transcript and a readthrough one.
#[test]
fn a_model_whose_every_transcript_is_set_aside_is_refused_and_writes_nothing() {
    let ids = "gene_id \"G1\"; transcript_id \"T1\";";
    assert_model_refused(
        "readthrough.gtf",
        "--gtf",
        &format!(
            "21\tx\tgene\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; gene_biotype \"protein_coding\";\n\
             21\tx\ttranscript\t1000\t2000\t.\t+\t.\t{ids} transcript_biotype \"protein_coding\"; tag \"readthrough_transcript\";\n\
             21\tx\texon\t1000\t1200\t.\t+\t.\t{ids}\n"
        ),
        "No transcript was kept from {model}: 1 read, 0 set aside as artifact and 1 as \
         readthrough, as Ensembl VEP's cache holds neither",
    );
    assert_model_refused(
        "artifact_and_readthrough.gff3",
        "--gff3",
        "##gff-version 3\n\
         21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
         21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=artifact\n\
         21\tx\texon\t1000\t1200\t.\t+\t.\tParent=transcript:T1\n\
         21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:T2;Parent=gene:G1;biotype=protein_coding;tag=readthrough_transcript\n\
         21\tx\texon\t1500\t2000\t.\t+\t.\tParent=transcript:T2\n",
        "No transcript was kept from {model}: 2 read, 1 set aside as artifact and 1 as \
         readthrough, as Ensembl VEP's cache holds neither",
    );
}
