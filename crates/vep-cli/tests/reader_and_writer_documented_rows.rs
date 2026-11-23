// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes on the inputs of three documented
//! differences whose code stands in the readers and writers of `vep-io`: a tab
//! output of one column whose rows read `0`, a region whose allele is the word
//! `BND`, and the JSON rows of a record whose transcripts all read one
//! SpliceVault table row. Each test runs the input of its entry's record under
//! `scripts/divergences/reproductions` over the committed corpus that record
//! names and holds every row of the output.

use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::{json, Value};

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A committed corpus under `tests/golden/116`, by its directory name.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The reference FASTA of the identifiers corpus, which covers the records of
/// the chromosome 21 cache slices.
fn reference_fasta() -> String {
    corpus("GRCh38-identifiers")
        .join("reference.fa.gz")
        .display()
        .to_string()
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// The lines of the output `vep` writes to `<dir>/<name>` for `input` over the
/// JSON cache of the corpus `cache`, with `flags` after the arguments every run
/// of this file takes. The run ends with exit code 0.
fn output_lines(dir: &Path, name: &str, input: &Path, cache: &str, flags: &[&str]) -> Vec<String> {
    let out = dir.join(name);
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus(cache).join("json_cache"))
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name} {flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// The lines of a tab output from its column line on: the `##` lines dropped.
fn from_the_column_line(lines: &[String]) -> Vec<&str> {
    lines
        .iter()
        .map(String::as_str)
        .filter(|line| !line.starts_with("##"))
        .collect()
}

/// A tab row whose whole line is `0` is written, and so is every row after it.
/// Under `--shift_3prime 1 --shift_length` the list `--fields SHIFT_LENGTH`
/// makes the two rows of the substitution the line `0` and the three rows of
/// the insertion the line `10`; the list `Feature,SHIFT_LENGTH` names the
/// transcript of each of the five. Ensembl VEP 116.3 writes the header and no
/// row for the one-column list: its output loop, `while(my $line =
/// $self->next_output_line)` (`Runner.pm:208`), ends at the first line Perl
/// reads as false, and the string `0` is one.
#[test]
fn a_one_column_tab_row_that_reads_zero_is_written_with_the_rows_after_it() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "shifted.vcf",
        &format!(
            "{VCF_HEADER}\
             21\t13615721\tintron_snv_2\tA\tC\t.\t.\t.\n\
             21\t13620283\tcds_tri_ins_fwd_1\tG\tGTCA\t.\t.\t.\n"
        ),
    );
    let fasta = reference_fasta();
    let tab = |name: &str, fields: &str| {
        output_lines(
            dir.path(),
            name,
            &input,
            "GRCh38-identifiers",
            &[
                "--fasta",
                &fasta,
                "--tab",
                "--shift_3prime",
                "1",
                "--shift_length",
                "--fields",
                fields,
            ],
        )
    };
    assert_eq!(
        from_the_column_line(&tab("one-column.txt", "SHIFT_LENGTH")),
        ["#SHIFT_LENGTH", "0", "0", "10", "10", "10"]
    );
    assert_eq!(
        from_the_column_line(&tab("two-columns.txt", "Feature,SHIFT_LENGTH")),
        [
            "#Feature\tSHIFT_LENGTH",
            "ENST00000299443\t0",
            "ENST00000620442\t0",
            "ENST00000299443\t10",
            "ENST00000384745\t10",
            "ENST00000620442\t10",
        ]
    );
}

/// The transcripts of the region-input cache slice that 21:14221000-14222000
/// reaches, in the order of their rows, each with its terms and its `Extra`
/// cell.
const BREAKEND_ROWS: [(&str, &str, &str); 17] = [
    (
        "ENST00000400577",
        "coding_sequence_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00000461088",
        "non_coding_transcript_exon_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00000468643",
        "non_coding_transcript_exon_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00000468788",
        "non_coding_transcript_exon_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00000475864",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=2375;STRAND=1",
    ),
    (
        "ENST00000495055",
        "non_coding_transcript_exon_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00000900829",
        "coding_sequence_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001045222",
        "coding_sequence_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001045223",
        "3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001045224",
        "coding_sequence_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001062921",
        "3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001062922",
        "3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001062923",
        "3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001103284",
        "3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001103285",
        "3_prime_UTR_variant,intron_variant,NMD_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001124786",
        "coding_sequence_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
    (
        "ENST00001131557",
        "coding_sequence_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ),
];

/// A region whose allele is the word `BND` is a `chromosome_breakpoint` over
/// the span as written, one row for each transcript the span reaches, and the
/// region on the next line is annotated after it. Ensembl VEP 116.3 builds the
/// word as a structural variant with no allele string (`Parser/Region.pm` lines
/// 208-216), then asks the record for its breakends, of which a region names
/// none, and dies (`AnnotationSource.pm:251`, exit code 255) with no row
/// written for the record or for the lines after it.
#[test]
fn a_region_whose_allele_is_the_word_bnd_is_a_chromosome_breakpoint_over_its_span() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "breakend.region.txt",
        "21:14221000-14222000/BND\n21:14221107-14221107/A\n",
    );
    let fasta = reference_fasta();
    let lines = output_lines(
        dir.path(),
        "breakend.txt",
        &input,
        "GRCh38-input-region",
        &["--format", "region", "--fasta", &fasta],
    );
    let rows: Vec<&str> = lines
        .iter()
        .map(String::as_str)
        .filter(|line| !line.starts_with('#'))
        .collect();
    let expected: Vec<String> = BREAKEND_ROWS
        .iter()
        .map(|(transcript, terms, extra)| {
            format!(
                "21:14221000-14222000/BND\t21:14221000-14222000\tchromosome_breakpoint\t\
                 ENSG00000185272\t{transcript}\tTranscript\t{terms}\t-\t-\t-\t-\t-\t-\t{extra}"
            )
        })
        .collect();
    let of_the_word: Vec<&str> = rows
        .iter()
        .copied()
        .filter(|row| row.starts_with("21:14221000-14222000/BND\t"))
        .collect();
    assert_eq!(of_the_word, expected);
    assert_eq!(
        &rows[..expected.len()],
        of_the_word,
        "the rows of the word come first"
    );
    let after: Vec<&str> = rows[expected.len()..].to_vec();
    assert_eq!(after.len(), BREAKEND_ROWS.len(), "{after:?}");
    for (row, (transcript, _, _)) in after.iter().zip(BREAKEND_ROWS) {
        assert!(
            row.starts_with(&format!(
                "21:14221107-14221107/A\t21:14221107\tA\tENSG00000185272\t{transcript}\tTranscript\t"
            )),
            "{row}"
        );
    }
}

/// Under `transcript_match=0` every transcript row of a record reads the first
/// SpliceVault table row at the site whose allele matches, and in the JSON
/// format each of them carries the `splicevault` object: the two rows after the
/// record's first hold the object the first holds. Ensembl VEP 116.3 writes the
/// object on the first row alone: the plugin replaces the event list of the
/// table row it cached with the object it builds from it (`SpliceVault.pm`
/// lines 180-191), the next row that reads the cached row meets the object
/// where line 183 expects a list, and the row is written without its plugin
/// fields under the warning
/// `Plugin 'SpliceVault' went wrong: Not an ARRAY reference`. The table is the
/// corpus's made-up one.
#[test]
fn every_transcript_row_of_a_record_carries_its_splicevault_object_in_json() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "donor.vcf",
        &format!("{VCF_HEADER}21\t13618457\tdonor_fwd_two_transcripts\tT\tC\t.\t.\t.\n"),
    );
    let table = corpus("GRCh38-plugin-splicevault-fixture")
        .join("data/SpliceVault_data_GRCh38_synthetic.tsv.gz");
    let lines = output_lines(
        dir.path(),
        "donor.json",
        &input,
        "GRCh38-plugin-splicevault-fixture",
        &[
            "--json",
            "--plugin",
            &format!("SpliceVault,file={},transcript_match=0", table.display()),
        ],
    );
    assert_eq!(lines.len(), 1, "{lines:?}");
    let record: Value = serde_json::from_str(&lines[0]).unwrap();
    let object = json!({
        "SpliceVault_site_sample_count": 44461,
        "SpliceVault_site_max_depth": 191,
        "SpliceVault_SpliceAI_delta": 0.51,
        "SpliceVault_out_of_frame_events": "Frameshift:0/2",
        "SpliceVault_top_events": {
            "Top1": {"site_type": "CA", "transcript_impact": 449,
                     "percent_of_supporting_samples": "0.493%", "frame": "inFrame"},
            "Top2": {"site_type": "CD", "transcript_impact": 44,
                     "percent_of_supporting_samples": "30%", "frame": "inFrame"},
        },
        "SpliceVault_site_pos": "chr21:13618451",
        "SpliceVault_site_type": "Acceptor_loss",
    });
    let donor_row = |transcript: &str| {
        json!({
            "consequence_terms": ["splice_donor_variant"],
            "gene_id": "ENSG00000166351",
            "impact": "HIGH",
            "strand": 1,
            "transcript_id": transcript,
            "variant_allele": "C",
            "splicevault": object,
        })
    };
    assert_eq!(
        record["transcript_consequences"],
        json!([
            donor_row("ENST00000299443"),
            {
                "consequence_terms": ["upstream_gene_variant"],
                "distance": 3120,
                "gene_id": "ENSG00000207476",
                "impact": "MODIFIER",
                "strand": 1,
                "transcript_id": "ENST00000384745",
                "variant_allele": "C",
                "splicevault": object,
            },
            donor_row("ENST00000620442"),
        ])
    );
}
