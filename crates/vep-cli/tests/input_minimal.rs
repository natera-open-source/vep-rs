// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary runs under `--minimal`: the flag switches on
//! `--allele_number`, every bi-allelic record is trimmed and marked, a
//! multi-allelic record whose alleles trim is split and rejoined as Ensembl VEP
//! 116.2 rejoins it (the record's Location, the trimmed Allele, the written
//! position as ALLELE_NUM, the rows in the rejoin's key order) and the three
//! position columns and the `--coding_only` rows follow the output format.
//! Every expected value is Ensembl VEP 116.2's output on the record under the
//! same flags.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The identifiers corpus: a chromosome 21 cache slice over coding exons of TPTE.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A record whose four ALTs of three shapes split under `--minimal`: chopped to
/// `AAG/GAG/AGG/-` at 21:10542415-10542417, its alleles trim to `A/G` at 10542415
/// (ALLELE_NUM 1, `missense_variant`), `A/G` at 10542416 (2, `synonymous_variant`)
/// and `AAG/-` (3, `inframe_deletion`) on ENST00000427445.
const MULTI_THREE: &str = "21\t10542414\tmulti_three\tAAAG\tAGAG,AAGG,A\t.\t.\t.\n";

/// A bi-allelic substitution trimmed on both sides to `AA/GG` at 10542415-10542416.
const MNP: &str = "21\t10542414\tmnp_two_sided\tAAAG\tAGGG\t.\t.\t.\n";

fn vep(dir: &Path, records: &str, format: &[&str], extra: &[&str]) -> Output {
    vep_on(
        dir,
        "in.vcf",
        &format!("{VCF_HEADER}{records}"),
        format,
        extra,
    )
}

/// A run on the input file `name` holding `input_text` as it stands, `format` the
/// output format flags and `extra` the run's other flags, `--format` among them
/// for a non-VCF input.
fn vep_on(dir: &Path, name: &str, input_text: &str, format: &[&str], extra: &[&str]) -> Output {
    vep_in(&corpus(), dir, name, input_text, format, extra)
}

/// [`vep_on`] against the cache of the corpus at `corpus_dir`.
fn vep_in(
    corpus_dir: &Path,
    dir: &Path,
    name: &str,
    input_text: &str,
    format: &[&str],
    extra: &[&str],
) -> Output {
    let input = dir.join(name);
    std::fs::write(&input, input_text).unwrap();
    let out = dir.join("out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus_dir.join("json_cache"))
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(format)
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

/// The data lines of the run's output.
fn data_lines(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("out"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The rows on one transcript as `(Location, Allele, cDNA_position, CDS_position,
/// Protein_position, Consequence, Extra)`.
fn rows_on(lines: &[String], transcript: &str) -> Vec<Vec<String>> {
    lines
        .iter()
        .filter(|line| line.contains(transcript))
        .map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            [c[1], c[2], c[7], c[8], c[9], c[6], c[13]]
                .iter()
                .map(|s| s.to_string())
                .collect()
        })
        .collect()
}

/// `--minimal` switches on `--allele_number` (the option set at Config.pm:524-529):
/// every row prints `ALLELE_NUM`, and a bi-allelic record `MINIMISED=1` whether a
/// base was trimmed (`AAAG/AGGG` to `AA/GG`) or not (an SNV).
#[test]
fn minimal_implies_allele_number() {
    let dir = tempfile::tempdir().unwrap();
    let snv = "21\t10542415\tsnv\tA\tG\t.\t.\t.\n";
    vep(
        dir.path(),
        &format!("{snv}{MNP}"),
        &[],
        &["--minimal", "--uploaded_allele"],
    );
    let lines = data_lines(dir.path());
    assert!(lines.len() > 2, "{lines:?}");
    for line in &lines {
        let extra = line.split('\t').nth(13).unwrap();
        assert!(extra.contains("ALLELE_NUM=1"), "{line}");
        assert!(extra.contains("MINIMISED=1"), "{line}");
    }
    let mnp = rows_on(&lines, "ENST00000427445");
    let mnp: Vec<&Vec<String>> = mnp
        .iter()
        .filter(|r| r[6].contains("UPLOADED_ALLELE=AAAG/AGGG"))
        .collect();
    assert_eq!(mnp.len(), 1, "{lines:?}");
    assert_eq!(
        &mnp[0][..6],
        [
            "21:10542415-10542416",
            "GG",
            "438-439",
            "86-87",
            "29",
            "missense_variant",
        ]
    );
    // Without the flag the record keeps its alleles and no ALLELE_NUM is printed.
    vep(dir.path(), MNP, &[], &["--uploaded_allele"]);
    let plain = rows_on(&data_lines(dir.path()), "ENST00000427445");
    assert_eq!(plain[0][0], "21:10542414-10542417");
    assert_eq!(plain[0][1], "AGGG");
    assert!(!plain[0][6].contains("ALLELE_NUM"), "{}", plain[0][6]);
}

/// An Ensembl-format record is trimmed in the case it was written: `post_process_vfs`
/// runs `minimise_alleles` on the column as written (Parser.pm:852), whose
/// `trim_sequences` compares bases with `eq` (Sequence.pm:1006, :1019), then
/// `validate_vf` names an ID-less record from the trimmed string (Parser.pm:462) and
/// upper-cases it (:554). So `21 10542414 10542417 AAAG/AGGg +` loses its shared `A`
/// and keeps the `G`/`g` pair, and Ensembl VEP 116.2 writes `21_10542415_AAG/GGg` at
/// `21:10542415-10542417`, Allele `GGG`, positions `438-440`, `86-88` and `29-30`,
/// `missense_variant`, `UPLOADED_ALLELE=AAAG/AGGg` and `MINIMISED=1` on
/// ENST00000427445; the same column in upper case is `AA/GG` at 10542415-10542416.
/// `UPLOADED_ALLELE` prints the column as written for every record here. Ensembl
/// VEP prints it so for a record whose trimmed alleles are of one length or carry a
/// `-`; for an indel that trims to two lengths, `21 10542414 10542417 ATGC/AG +`
/// to `TGC/G`, its second pass over the trimmed record hands `TGC/G` to
/// `minimise_alleles` again and resets the uploaded string to it (Parser.pm:862-881,
/// :1008), printing `TGC/G` where vep-rs prints `ATGC/AG`.
#[test]
fn an_ensembl_format_record_is_trimmed_in_its_written_case() {
    let dir = tempfile::tempdir().unwrap();
    let flags = ["--format", "ensembl", "--minimal", "--uploaded_allele"];
    let annotate = |record: &str| {
        vep_on(dir.path(), "in.txt", record, &[], &flags);
        data_lines(dir.path())
    };
    let lines = annotate("21\t10542414\t10542417\tAAAG/AGGg\t+\n");
    let rows = rows_on(&lines, "ENST00000427445");
    assert_eq!(rows.len(), 1, "{lines:?}");
    assert_eq!(
        &rows[0][..6],
        [
            "21:10542415-10542417",
            "GGG",
            "438-440",
            "86-88",
            "29-30",
            "missense_variant",
        ]
    );
    assert!(
        rows[0][6].contains("UPLOADED_ALLELE=AAAG/AGGg;"),
        "{}",
        rows[0][6]
    );
    assert!(rows[0][6].contains("MINIMISED=1"), "{}", rows[0][6]);
    for line in &lines {
        assert!(line.starts_with("21_10542415_AAG/GGg\t"), "{line}");
    }

    let upper = rows_on(
        &annotate("21\t10542414\t10542417\tAAAG/AGGG\t+\n"),
        "ENST00000427445",
    );
    assert_eq!(&upper[0][..2], ["21:10542415-10542416", "GG"]);

    let indel = rows_on(
        &annotate("21\t10542414\t10542417\tATGC/AG\t+\n"),
        "ENST00000427445",
    );
    assert_eq!(&indel[0][..2], ["21:10542415-10542417", "G"]);
    assert!(
        indel[0][6].contains("UPLOADED_ALLELE=ATGC/AG;"),
        "{}",
        indel[0][6]
    );
}

/// A VCF record is trimmed in its written case under `--minimal` too, and a record
/// whose ID Perl reads as false is named from the trimmed string before the
/// upper-casing (Parser/VCF.pm:343 keeps the string as written, Parser.pm:852, :462,
/// :554): Ensembl VEP 116.2 writes `21 10542414 0 AAAG AGGg` as `21_10542415_AAG/GGg`
/// at `21:10542415-10542417`, Allele `GGG`, `UPLOADED_ALLELE=AAAG/AGGg`.
#[test]
fn a_vcf_record_is_trimmed_in_its_written_case() {
    let dir = tempfile::tempdir().unwrap();
    vep(
        dir.path(),
        "21\t10542414\t0\tAAAG\tAGGg\t.\t.\t.\n",
        &[],
        &["--minimal", "--uploaded_allele"],
    );
    let lines = data_lines(dir.path());
    let rows = rows_on(&lines, "ENST00000427445");
    assert_eq!(rows.len(), 1, "{lines:?}");
    assert_eq!(&rows[0][..2], ["21:10542415-10542417", "GGG"]);
    assert!(
        rows[0][6].contains("UPLOADED_ALLELE=AAAG/AGGg;"),
        "{}",
        rows[0][6]
    );
    for line in &lines {
        assert!(line.starts_with("21_10542415_AAG/GGg\t"), "{line}");
    }
}

/// The rejoin of a split record (OutputFactory.pm:2272-2354): every row prints the
/// record's Location `21:10542415-10542417`, the trimmed Allele and the written
/// position as ALLELE_NUM, the rows sorted by the keys the merged alleles are stored
/// under (`A/G_<id>` and `AAG/-_<id>` before the first allele's bare `<id>`), each
/// `MINIMISED=1` with the chopped string as `UPLOADED_ALLELE`, the first trimmed
/// pair's class, `SNV`, as `VARIANT_CLASS`, and under `--show_ref_allele` the first
/// allele's trimmed reference, `A`, as `REF_ALLELE` (OutputFactory.pm:1357 runs on
/// the first allele before the rejoin restores the record's string and the cached
/// value serves every merged row). The default format prints the record
/// span's positions (`438-440`, `86-88`, `29-30`) on every row; the VCF format
/// prints each allele's own (`439`, `87`, `29` for the synonymous `G`), as the two
/// factories reset the cached positions after and before the rejoin
/// (OutputFactory.pm:294-299, OutputFactory/VCF.pm:305-308).
#[test]
fn a_split_record_rejoins_as_ensembl_vep_rejoins_it() {
    let dir = tempfile::tempdir().unwrap();
    let flags = ["--minimal", "--uploaded_allele", "--variant_class"];
    vep(dir.path(), MULTI_THREE, &[], &flags);
    let rows = rows_on(&data_lines(dir.path()), "ENST00000427445");
    let extra = "UPLOADED_ALLELE=AAG/GAG/AGG/-;IMPACT=";
    let shape: Vec<(String, String, String, String, String, String)> = rows
        .iter()
        .map(|r| {
            let num = r[6]
                .split(';')
                .find_map(|kv| kv.strip_prefix("ALLELE_NUM="))
                .unwrap()
                .to_string();
            (
                r[0].clone(),
                r[1].clone(),
                num,
                [r[2].as_str(), r[3].as_str(), r[4].as_str()].join("|"),
                r[5].clone(),
                r[6].clone(),
            )
        })
        .collect();
    let location = "21:10542415-10542417".to_string();
    let positions = "438-440|86-88|29-30".to_string();
    assert_eq!(
        shape
            .iter()
            .map(|(l, a, n, p, c, _)| (l.as_str(), a.as_str(), n.as_str(), p.as_str(), c.as_str()))
            .collect::<Vec<_>>(),
        [
            (
                location.as_str(),
                "G",
                "2",
                positions.as_str(),
                "synonymous_variant"
            ),
            (
                location.as_str(),
                "-",
                "3",
                positions.as_str(),
                "inframe_deletion"
            ),
            (
                location.as_str(),
                "G",
                "1",
                positions.as_str(),
                "missense_variant"
            ),
        ]
    );
    for (_, _, _, _, _, e) in &shape {
        assert!(e.contains(extra), "{e}");
        assert!(e.contains("VARIANT_CLASS=SNV;MINIMISED=1"), "{e}");
    }

    vep(
        dir.path(),
        MULTI_THREE,
        &[],
        &["--minimal", "--show_ref_allele"],
    );
    let with_reference = rows_on(&data_lines(dir.path()), "ENST00000427445");
    assert_eq!(with_reference.len(), 3);
    for row in &with_reference {
        assert!(row[6].contains("REF_ALLELE=A;"), "{}", row[6]);
    }

    vep(dir.path(), MULTI_THREE, &["--vcf"], &flags);
    let csq: Vec<String> = data_lines(dir.path())
        .pop()
        .unwrap()
        .split('\t')
        .nth(7)
        .unwrap()
        .trim_start_matches("CSQ=")
        .split(',')
        .filter(|entry| entry.contains("ENST00000427445"))
        .map(|entry| {
            let f: Vec<&str> = entry.split('|').collect();
            // Allele, cDNA_position, CDS_position, Protein_position, ALLELE_NUM.
            format!("{}|{}|{}|{}|{}", f[0], f[12], f[13], f[14], f[18])
        })
        .collect();
    assert_eq!(
        csq,
        ["G|439|87|29|2", "-|438-440|86-88|29-30|3", "G|438|86|29|1"]
    );
}

/// A `*` the split's trim leaves writes no row and keeps its number. The ALT `*A`
/// of `21 10542415 star_embedded A *A,G` is exempt from the anchor chop
/// (Parser/VCF.pm:315) and, under `--minimal`, trims against `A` to `-/*`
/// (InputBuffer.pm:665), a change that splits the record; Ensembl VEP numbers
/// the feature's allele `*` and passes it over (VariationFeatureOverlap.pm:484),
/// as it does a `*` written so, and the `G` allele keeps its written position
/// (`_base_allele_number`, InputBuffer.pm:684). So the record writes one row on
/// ENST00000427445: the record's Location, Allele `G`, `ALLELE_NUM=2`,
/// `UPLOADED_ALLELE=A/*A/G`, `MINIMISED=1` and the first feature's class,
/// `sequence_alteration`, as `VARIANT_CLASS`; with the `*A` written second the
/// `G` row is numbered 1 and the class is the first pair's `SNV`. `REF_ALLELE` is
/// the record's own `A`, the `*` feature having built no row to cache its
/// trimmed `-` (OutputFactory.pm:1357 runs per row). Without `--minimal` the
/// record is not split and the allele `*A`, which is not `*`, writes its row,
/// numbered 1.
#[test]
fn a_star_the_split_leaves_writes_no_row_and_keeps_its_number() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t10542415\tstar_embedded\tA\t*A,G\t.\t.\t.\n\
                   21\t10542415\tstar_second\tA\tG,*A\t.\t.\t.\n";
    vep(
        dir.path(),
        records,
        &[],
        &["--minimal", "--uploaded_allele", "--variant_class"],
    );
    let lines = data_lines(dir.path());
    let record = |lines: &[String], name: &str| -> Vec<Vec<String>> {
        let prefix = format!("{name}\t");
        let own: Vec<String> = lines
            .iter()
            .filter(|line| line.starts_with(&prefix))
            .cloned()
            .collect();
        rows_on(&own, "ENST00000427445")
    };
    let embedded = record(&lines, "star_embedded");
    assert_eq!(embedded.len(), 1, "{lines:?}");
    assert_eq!(
        &embedded[0][..6],
        ["21:10542415", "G", "438", "86", "29", "missense_variant"]
    );
    for key in [
        "ALLELE_NUM=2;",
        "UPLOADED_ALLELE=A/*A/G;",
        "VARIANT_CLASS=sequence_alteration;",
        "MINIMISED=1",
    ] {
        assert!(embedded[0][6].contains(key), "{key}: {}", embedded[0][6]);
    }
    let second = record(&lines, "star_second");
    assert_eq!(second.len(), 1, "{lines:?}");
    assert_eq!(&second[0][..2], ["21:10542415", "G"]);
    for key in ["ALLELE_NUM=1;", "VARIANT_CLASS=SNV;", "MINIMISED=1"] {
        assert!(second[0][6].contains(key), "{key}: {}", second[0][6]);
    }
    assert!(
        lines
            .iter()
            .all(|line| line.split('\t').nth(2) != Some("*")),
        "{lines:?}"
    );

    vep(
        dir.path(),
        records,
        &[],
        &["--minimal", "--show_ref_allele"],
    );
    let lines = data_lines(dir.path());
    let rows = rows_on(&lines, "ENST00000427445");
    assert_eq!(rows.len(), 2, "{lines:?}");
    for row in &rows {
        assert!(row[6].contains("REF_ALLELE=A;"), "{}", row[6]);
    }

    vep(
        dir.path(),
        records,
        &[],
        &["--allele_number", "--uploaded_allele"],
    );
    let lines = data_lines(dir.path());
    let embedded = record(&lines, "star_embedded");
    assert_eq!(embedded.len(), 2, "{lines:?}");
    assert_eq!(&embedded[0][..2], ["21:10542415", "*A"]);
    assert!(
        embedded[0][6].contains("ALLELE_NUM=1;"),
        "{}",
        embedded[0][6]
    );
    assert_eq!(&embedded[1][..2], ["21:10542415", "G"]);
    assert!(
        embedded[1][6].contains("ALLELE_NUM=2;"),
        "{}",
        embedded[1][6]
    );
}

/// The attributes corpus: a chromosome 21 cache whose transcripts carry the
/// miRNA structure attributes.
fn attributes_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-attributes")
}

/// A record on the hairpin transcript ENST00000581792 whose alleles split under
/// `--minimal` into the deletion `CT/-` at 8987004-8987005, cDNA 6-7, and the
/// insertion `-/TTA` after 8987005, cDNA 7-8.
const MIRNA_SPLIT: &str = "21\t8987003\tmulti\tACT\tA,ACTTTA\t.\t.\t.\n";

/// `(Allele, cDNA_position, Consequence, ALLELE_NUM, miRNA)` of every row on
/// `transcript`.
fn mirna_rows(lines: &[String], transcript: &str) -> Vec<[String; 5]> {
    rows_on(lines, transcript)
        .iter()
        .map(|r| {
            let extra = |key: &str| {
                r[6].split(';')
                    .find_map(|kv| kv.strip_prefix(key))
                    .unwrap_or_default()
                    .to_string()
            };
            [
                r[1].clone(),
                r[2].clone(),
                r[5].clone(),
                extra("ALLELE_NUM="),
                extra("miRNA="),
            ]
        })
        .collect()
}

/// A split allele's miRNA structure terms are computed over the allele's own
/// span. `mirna_structure_terms` runs once per split allele on the allele it
/// analyses (`calculate_consequences` in
/// `crates/vep-effects/src/consequences.rs`), so the insertion `TTA`, cDNA 7-8,
/// gets the terms `miRNA_loop` and `miRNA_stem` and the deletion, cDNA 6-7,
/// `miRNA_loop` alone, in the default output (`miRNA=miRNA_loop,miRNA_stem`) and
/// the VCF output (`miRNA_loop&miRNA_stem`) alike, while the default format's
/// cDNA_position column prints the rejoined record's span, `6-7`, on both rows.
/// Ensembl VEP 116.2's default factory rejoins the record
/// (OutputFactory.pm:294 `$self->rejoin_variants_in_InputBuffer($buffer) if
/// $buffer->rejoin_required;`) and then clears every allele's cached cDNA
/// coordinates (:297-298 `map {@{$self->reset_shifted_positions($_)}}
/// @{$buffer->buffer};`, through `clear_shifting_variables` at :318, which
/// deletes `cdna_start`, `cdna_end` and `_cdna_coords`,
/// TranscriptVariationAllele.pm:515-517), so its miRNA block
/// (OutputFactory.pm:1630 `my ($cdna_start, $cdna_end) = ($tv->cdna_start,
/// $tv->cdna_end);`) reads the span afresh over the rejoined record, cDNA 6-7,
/// and Ensembl VEP 116.2 prints `miRNA=miRNA_loop` for the `TTA` allele. Its VCF
/// factory clears the coordinates before the rejoin (OutputFactory/VCF.pm:305-308),
/// so the `TTA` CSQ entry carries `miRNA_loop&miRNA_stem` at cDNA `7-8` in both
/// engines.
#[test]
fn a_split_alleles_mirna_terms_come_from_its_own_span_in_the_default_output() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = attributes_corpus();
    let input = format!("{VCF_HEADER}{MIRNA_SPLIT}");
    let flags = ["--mirna", "--minimal"];
    vep_in(&corpus, dir.path(), "in.vcf", &input, &[], &flags);
    let lines = data_lines(dir.path());
    assert_eq!(
        mirna_rows(&lines, "ENST00000581792"),
        [
            [
                "TTA",
                "6-7",
                "non_coding_transcript_exon_variant",
                "2",
                "miRNA_loop,miRNA_stem",
            ],
            [
                "-",
                "6-7",
                "non_coding_transcript_exon_variant",
                "1",
                "miRNA_loop",
            ],
        ]
    );
    for line in &lines {
        assert!(line.contains("\t21:8987004-8987005\t"), "{line}");
    }

    vep_in(&corpus, dir.path(), "in.vcf", &input, &["--vcf"], &flags);
    let csq: Vec<String> = data_lines(dir.path())
        .pop()
        .unwrap()
        .split('\t')
        .nth(7)
        .unwrap()
        .trim_start_matches("CSQ=")
        .split(',')
        .filter(|entry| entry.contains("ENST00000581792"))
        .map(|entry| {
            let f: Vec<&str> = entry.split('|').collect();
            // Allele, cDNA_position, ALLELE_NUM, miRNA.
            format!("{}|{}|{}|{}", f[0], f[12], f[18], f[25])
        })
        .collect();
    assert_eq!(
        csq,
        ["TTA|7-8|2|miRNA_loop&miRNA_stem", "-|6-7|1|miRNA_loop"]
    );
}

/// Three records at the acceptor site of ENST00000427445's exon 6, each splitting
/// under `--minimal` into alleles on both sides of the CDS boundary: an SNV at
/// 10542393, two intronic bases before the exon (`splice_acceptor_variant`), and
/// a deletion of `TACG` over 10542394-10542397, which reaches the CDS
/// (`splice_acceptor_variant,coding_sequence_variant`).
const STRADDLE_A: &str = "21\t10542393\tstraddle_a\tGTACG\tCTACG,G\t.\t.\t.\n";
const STRADDLE_B: &str = "21\t10542393\tstraddle_b\tGTACG\tCTACG,TTACG\t.\t.\t.\n";
const STRADDLE_C: &str = "21\t10542393\tstraddle_c\tGTACG\tG,CTACG\t.\t.\t.\n";

/// `(Uploaded_variation, Allele, Consequence, ALLELE_NUM, REF_ALLELE)` of every
/// row on `transcript`.
fn straddle_rows(lines: &[String], transcript: &str) -> Vec<[String; 5]> {
    lines
        .iter()
        .filter(|line| line.contains(transcript))
        .map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            let extra = |key: &str| {
                c[13]
                    .split(';')
                    .find_map(|kv| kv.strip_prefix(key))
                    .unwrap_or_default()
                    .to_string()
            };
            [
                c[0].to_string(),
                c[2].to_string(),
                c[6].to_string(),
                extra("ALLELE_NUM="),
                extra("REF_ALLELE="),
            ]
        })
        .collect()
}

/// `--coding_only` judges a split record's rows over the span the format rejoins
/// on. The default and tab factories clear the cached CDS coordinates after the
/// rejoin (OutputFactory.pm:297-302), so `affects_cds` maps the record's restored
/// span, and Ensembl VEP 116.2 keeps every row of the seven coding transcripts for
/// both alleles of `straddle_a` and `straddle_b` (14 rows each), the intronic SNV's
/// `splice_acceptor_variant` rows included. The VCF and JSON factories clear them
/// before the rejoin (OutputFactory/VCF.pm:305-308, OutputFactory/JSON.pm:217-220),
/// so each allele is judged on its own span: `straddle_a` keeps the deletion's
/// seven entries alone and `straddle_b` writes no CSQ. `REF_ALLELE` is cached on
/// the first row the first feature builds before the rejoin (OutputFactory.pm:1357
/// at :2289): the record's own `GTACG` where that feature's SNV kept no row
/// (`straddle_a`, `straddle_b`), the first allele's trimmed `TACG` where the
/// deletion came first (`straddle_c`).
#[test]
fn coding_only_judges_a_split_record_over_the_span_the_format_rejoins_on() {
    let dir = tempfile::tempdir().unwrap();
    let flags = [
        "--minimal",
        "--coding_only",
        "--show_ref_allele",
        "--uploaded_allele",
        "--variant_class",
    ];
    vep(
        dir.path(),
        &format!("{STRADDLE_A}{STRADDLE_B}{STRADDLE_C}"),
        &[],
        &flags,
    );
    let lines = data_lines(dir.path());
    let count = |name: &str| lines.iter().filter(|l| l.starts_with(name)).count();
    assert_eq!(
        (
            count("straddle_a"),
            count("straddle_b"),
            count("straddle_c")
        ),
        (14, 14, 14)
    );
    let intronic = ["ENST00000612746", "ENST00000612957"];
    assert!(
        lines
            .iter()
            .all(|l| intronic.iter().all(|t| !l.contains(t))),
        "{lines:?}"
    );
    assert_eq!(
        straddle_rows(&lines, "ENST00000427445"),
        [
            ["straddle_a", "C", "splice_acceptor_variant", "1", "GTACG"],
            [
                "straddle_a",
                "-",
                "splice_acceptor_variant,coding_sequence_variant",
                "2",
                "GTACG",
            ],
            ["straddle_b", "C", "splice_acceptor_variant", "1", "GTACG"],
            ["straddle_b", "T", "splice_acceptor_variant", "2", "GTACG"],
            [
                "straddle_c",
                "-",
                "splice_acceptor_variant,coding_sequence_variant",
                "1",
                "TACG",
            ],
            ["straddle_c", "C", "splice_acceptor_variant", "2", "TACG"],
        ]
    );
    for line in &lines {
        assert!(line.contains("\t21:10542393-10542397\t"), "{line}");
    }

    vep(
        dir.path(),
        &format!("{STRADDLE_A}{STRADDLE_B}"),
        &["--vcf"],
        &flags,
    );
    let vcf = data_lines(dir.path());
    assert_eq!(vcf.len(), 2);
    let info = |line: &String| line.split('\t').nth(7).unwrap().to_string();
    let csq: Vec<String> = info(&vcf[0])
        .trim_start_matches("CSQ=")
        .split(',')
        .map(|entry| {
            let f: Vec<&str> = entry.split('|').collect();
            // Allele, Consequence, ALLELE_NUM, REF_ALLELE.
            format!("{}|{}|{}|{}", f[0], f[1], f[18], f[19])
        })
        .collect();
    assert_eq!(csq.len(), 7, "{csq:?}");
    assert!(
        csq.iter().all(|e| e == "-|splice_acceptor_variant&coding_sequence_variant|2|GTACG"
            || e == "-|splice_acceptor_variant&coding_sequence_variant&NMD_transcript_variant|2|GTACG"),
        "{csq:?}"
    );
    assert_eq!(info(&vcf[1]), ".");
}
