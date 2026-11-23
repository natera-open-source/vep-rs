// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes on the inputs of the documented
//! differences whose code stands in the built-in plugins of `vep-builtins`.
//! Each test runs the input of its entry's record under
//! `scripts/divergences/reproductions` over the committed corpus that record
//! names, with the record's flags, and holds every row of the output in the
//! columns the entry is about. A file the record prepares is written again in
//! the test's own directory. The record of `gene_score_zero_value` runs
//! LoFtool alone, so the entry's pLI test runs the same input over the pLI
//! corpus. The record of `plugin_cols_names_matched_exactly` runs OpenTargets
//! alone, so the entry's AlphaMissense and TranscriptAnnotator tests run the
//! records its page works for those plugins, each over its plugin's corpus.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The fourteen columns of the default output, in order.
const COLUMNS: [&str; 14] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "Extra",
];

/// A committed corpus under `tests/golden/116`, by its directory name.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// `path` as an argument word.
fn word(path: &Path) -> String {
    path.display().to_string()
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// A VCF of `records` at `<dir>/<name>` under the two header lines.
fn vcf(dir: &Path, name: &str, records: &str) -> PathBuf {
    write(dir, name, &format!("{VCF_HEADER}{records}"))
}

/// A `vep` run of `input` over the JSON cache of the corpus `cache` that
/// writes `out`, with the arguments every run of this file takes; a test adds
/// its record's own.
fn vep(dir: &Path, input: &Path, out: &Path, cache: &str) -> Command {
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--json_cache")
        .arg(corpus(cache).join("json_cache"))
        .args(["--offline", "--force_overwrite", "--no_stats", "--quiet"])
        .stdin(Stdio::null());
    command
}

/// The lines of `out` after `command` has ended with exit code 0.
fn lines_of(mut command: Command, out: &Path) -> Vec<String> {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{command:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// Every data row of a default output in the columns `names`, its cells joined
/// by ` | `. A name outside the fourteen columns is a key of the Extra column,
/// and a row without that key reads `-` there.
fn cells(lines: &[String], names: &[&str]) -> Vec<String> {
    lines
        .iter()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let row: Vec<&str> = line.split('\t').collect();
            assert_eq!(row.len(), COLUMNS.len(), "{line}");
            names
                .iter()
                .map(
                    |name| match COLUMNS.iter().position(|column| column == name) {
                        Some(at) => row[at],
                        None => row[COLUMNS.len() - 1]
                            .split(';')
                            .find_map(|pair| {
                                pair.split_once('=')
                                    .filter(|(key, _)| key == name)
                                    .map(|(_, value)| value)
                            })
                            .unwrap_or("-"),
                    },
                )
                .collect::<Vec<_>>()
                .join(" | ")
        })
        .collect()
}

/// A gene whose LoFtool percentile is written `0.00E+00` scores `0` on every
/// row of that gene. The scores file is the corpus's own with the percentile
/// of `SOD1` rewritten; the record's 47 rows are 26 of that gene, each `0`,
/// and 21 of three other genes the file does not hold. Ensembl VEP 116.3
/// writes no LoFtool on the 26: the plugin returns its score only where Perl
/// reads it as true (`LoFtool.pm:118`), and `0.00E+00` is stored as a `0`
/// (`LoFtool.pm:89`), which is false.
#[test]
fn a_loftool_percentile_of_zero_is_written_on_every_row_of_its_gene() {
    let dir = tempfile::tempdir().unwrap();
    let scores = std::fs::read_to_string(
        corpus("GRCh38-plugin-loftool").join("plugin_data/LoFtool_scores.txt"),
    )
    .unwrap();
    assert_eq!(scores.matches("sod1\t1.26E-01").count(), 1);
    let zeroed = write(
        dir.path(),
        "LoFtool_scores_zero.txt",
        &scores.replace("sod1\t1.26E-01", "sod1\t0.00E+00"),
    );
    let input = vcf(
        dir.path(),
        "sod1.vcf",
        "21\t31667290\t14766\tA\tC\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-loftool");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("LoFtool,{}", word(&zeroed)),
    ]);
    let rows = cells(&lines_of(command, &out), &["Gene", "LoFtool"]);
    let of = |gene: &str, score: &str| {
        rows.iter()
            .filter(|row| **row == format!("{gene} | {score}"))
            .count()
    };
    assert_eq!(rows.len(), 47, "{rows:#?}");
    assert_eq!(of("ENSG00000142168", "0"), 26, "{rows:#?}");
    assert_eq!(of("ENSG00000156304", "-"), 19, "{rows:#?}");
    assert_eq!(of("ENSG00000238390", "-"), 1, "{rows:#?}");
    assert_eq!(of("ENSG00000273271", "-"), 1, "{rows:#?}");
}

/// A gene whose pLI is written `0` scores `0` under `sigfig` on every row of
/// that gene. The input is the record's own, which the pLI corpus holds as
/// well, and the values file is that corpus's with the cell of `Sod1`
/// rewritten: the record's 47 rows are 26 of that gene, each `0`, 19 of a gene
/// whose value rounds to `1` and 2 of genes the file does not hold. Ensembl
/// VEP 116.3 writes no pLI on the 26: under `sigfig` the cell is stored
/// through `%.3g` as a `0` (`pLI.pm:110`, `pLI.pm:145`), and the plugin
/// returns its value only where Perl reads it as true (`pLI.pm:221`). Without
/// `sigfig` the cell is stored as `0.00`, which is true, and both engines
/// write it.
#[test]
fn a_pli_of_zero_is_written_under_sigfig_on_every_row_of_its_gene() {
    let dir = tempfile::tempdir().unwrap();
    let values =
        std::fs::read_to_string(corpus("GRCh38-plugin-pli").join("plugin_data/pLI_values.txt"))
            .unwrap();
    let cell = values
        .lines()
        .find(|line| line.starts_with("Sod1\t"))
        .unwrap();
    let zeroed = write(
        dir.path(),
        "pLI_values_zero.txt",
        &values.replace(cell, "Sod1\t0"),
    );
    let input = vcf(
        dir.path(),
        "sod1.vcf",
        "21\t31667290\t14766\tA\tC\t.\t.\t.\n",
    );
    let rows_under = |name: &str, parameters: &str| {
        let out = dir.path().join(name);
        let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-pli");
        command.args([
            "--assembly",
            "GRCh38",
            "--plugin",
            &format!("pLI,{}{parameters}", word(&zeroed)),
        ]);
        cells(&lines_of(command, &out), &["Gene", "pLI_gene_value"])
    };
    let of = |rows: &[String], gene: &str, value: &str| {
        rows.iter()
            .filter(|row| **row == format!("{gene} | {value}"))
            .count()
    };
    let rounded = rows_under("sigfig.txt", ",sigfig");
    assert_eq!(rounded.len(), 47, "{rounded:#?}");
    assert_eq!(of(&rounded, "ENSG00000142168", "0"), 26, "{rounded:#?}");
    assert_eq!(of(&rounded, "ENSG00000156304", "1"), 19, "{rounded:#?}");
    assert_eq!(of(&rounded, "ENSG00000238390", "-"), 1, "{rounded:#?}");
    assert_eq!(of(&rounded, "ENSG00000273271", "-"), 1, "{rounded:#?}");
    let two_decimals = rows_under("two-decimals.txt", "");
    assert_eq!(two_decimals.len(), 47, "{two_decimals:#?}");
    assert_eq!(
        of(&two_decimals, "ENSG00000142168", "0.00"),
        26,
        "{two_decimals:#?}"
    );
    assert_eq!(
        of(&two_decimals, "ENSG00000156304", "1.00"),
        19,
        "{two_decimals:#?}"
    );
    assert_eq!(
        of(&two_decimals, "ENSG00000238390", "-"),
        1,
        "{two_decimals:#?}"
    );
    assert_eq!(
        of(&two_decimals, "ENSG00000273271", "-"),
        1,
        "{two_decimals:#?}"
    );
}

/// A record whose chromosome is written `Chr21` takes the gnomADc coverage of
/// its position from the corpus's coverage file, whose contig is named
/// `chr21`: each of the record's 14 transcript rows carries the twelve
/// coverage fields, `gnomAD_mean_cov=31.1000` among them, and keeps the
/// input's spelling in its location. Ensembl VEP 116.3 writes none of the
/// twelve on the 14 rows: the plugin strips a lower-case `chr` alone
/// (`gnomADc.pm:158`), and the tabix base class tries a name the file lacks
/// without its prefix, or with `chr` added only where the name has no prefix
/// in any letter case (`BaseVepTabixPlugin.pm:693-708`), so the file is never
/// asked for `chr21`.
#[test]
fn a_gnomadc_record_spelt_with_a_capitalised_chr_prefix_takes_its_coverage() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "cap_chr.vcf",
        "Chr21\t32699959\tcap_chr\tG\tA\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-gnomadc");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!(
            "gnomADc,{}",
            word(&corpus("GRCh38-plugin-gnomadc").join("plugin_data/gnomad.coverage.tsv.gz"))
        ),
    ]);
    let lines = lines_of(command, &out);
    assert_eq!(
        cells(
            &lines,
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "gnomAD_mean_cov",
            ],
        ),
        [
            "cap_chr | Chr21:32699959 | A | ENST00000357345 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000382491 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000429236 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000456084 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000630077 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000674204 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000674308 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000674351 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000902244 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000902245 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00000963072 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00001034287 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00001075151 | missense_variant | 31.1000",
            "cap_chr | Chr21:32699959 | A | ENST00001110443 | missense_variant | 31.1000",
        ]
    );
    assert_eq!(
        cells(
            &lines,
            &[
                "gnomAD_100x_cov",
                "gnomAD_10x_cov",
                "gnomAD_15x_cov",
                "gnomAD_1x_cov",
                "gnomAD_20x_cov",
                "gnomAD_25x_cov",
                "gnomAD_30x_cov",
                "gnomAD_50x_cov",
                "gnomAD_5x_cov",
                "gnomAD_mean_cov",
                "gnomAD_median_approx_cov",
                "gnomAD_total_DP_cov",
            ],
        ),
        ["0.0003 | 1.0000 | 0.9999 | 1.0000 | 0.9996 | 0.9730 | 0.6856 | 0.0030 | 1.0000 | 31.1000 | 31.0000 | 22732808.0000"; 14]
    );
}

/// A copy of the AVI corpus's made-up score table at `<dir>/<name>` in which
/// the one row equal to `row` reads `replacement`, bgzipped and indexed on its
/// first two columns as `tabix -s 1 -b 2 -e 2` indexes it.
fn avi_table_with_row(dir: &Path, name: &str, row: &str, replacement: &str) -> PathBuf {
    use noodles::csi::binning_index::index::header::Builder;
    use std::io::Read;

    let table = corpus("GRCh38-plugin-avi-fixture").join("avi_scores.synthetic.tsv.gz");
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(table).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    assert_eq!(text.lines().filter(|line| *line == row).count(), 1);
    let rows: Vec<&str> = text
        .lines()
        .map(|line| if line == row { replacement } else { line })
        .collect();
    common::indexed_file(
        dir,
        name,
        Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(1))
            .build(),
        (1, 1, 0),
        &rows,
    )
}

/// A score-table row whose alleles are written in lower case scores the record
/// it describes: against a copy of the AVI corpus's made-up table whose `C T`
/// row at `chr21:5030633` reads `c t`, the record `21 5030633 C T` carries
/// `AVI_PHRED=3.465`, `AVI_RAW=1.6133` and a link that names the row's
/// alleles as the table writes them. Ensembl VEP 116.3 writes no `AVI_` field
/// on the row: the plugin hands the row's alleles to
/// `get_matched_variant_alleles` as the line holds them (`AVI.pm:161`,
/// `AVI.pm:139-151`), and the key `c_t_5030633` is compared as text with the
/// record's `C_T_5030633` (`Utils/Sequence.pm:1207`, `Utils/Sequence.pm:1232`).
#[test]
fn an_avi_row_written_in_lower_case_scores_the_record_it_describes() {
    let dir = tempfile::tempdir().unwrap();
    let table = avi_table_with_row(
        dir.path(),
        "avi_scores.lower.synthetic.tsv.gz",
        "chr21\t5030633\tC\tT\t1.6133\t3.465",
        "chr21\t5030633\tc\tt\t1.6133\t3.465",
    );
    let input = vcf(
        dir.path(),
        "substitution.vcf",
        "21\t5030633\t.\tC\tT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-avi-fixture");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("AVI,file={},raw=1,links=1", word(&table)),
    ]);
    let lines = lines_of(command, &out);
    assert_eq!(
        cells(
            &lines,
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "AVI_PHRED",
                "AVI_RAW",
            ],
        ),
        ["21_5030633_C/T | 21:5030633 | T | - | intergenic_variant | 3.465 | 1.6133"]
    );
    let links = cells(&lines, &["AVI_LINK"]);
    let queries: Vec<&str> = links
        .iter()
        .map(|link| link.split_once("?q=").map_or("-", |(_, query)| query))
        .collect();
    assert_eq!(queries, ["chr21:5030633:c>t&m=variant"]);
}

/// Under `--minimal` a record with two alternate alleles splits, and each
/// allele's rows carry the ancestral allele of that allele's own trimmed span.
/// The first allele trims to the deletion of one base and its 60 rows read
/// that base; the second trims to an insertion, whose span is empty, and its
/// 60 rows read `-`. Ensembl VEP 116.3 writes the first allele's base on all
/// 120 rows: the plugin stores its value on the feature it reads
/// (`AncestralAllele.pm:148`), and the rejoin points the later allele's rows
/// at the first allele's feature (`OutputFactory.pm:2323`), whose stored value
/// they print.
#[test]
fn a_later_allele_of_a_minimal_split_record_carries_the_ancestral_allele_of_its_own_span() {
    let dir = tempfile::tempdir().unwrap();
    let ancestor =
        corpus("GRCh38-plugin-ancestral").join("ancestral/homo_sapiens_ancestor_21.slice.fa.gz");
    let input = vcf(
        dir.path(),
        "anc_multi_indel.vcf",
        "21\t42750000\tanc_multi_indel\tTG\tT,TGT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-ancestral");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("AncestralAllele,{}", word(&ancestor)),
        "--minimal",
    ]);
    let lines = lines_of(command, &out);
    let rows = cells(&lines, &["Uploaded_variation", "Location", "Allele", "AA"]);
    let of = |row: &str| rows.iter().filter(|held| *held == row).count();
    assert_eq!(rows.len(), 120, "{rows:#?}");
    assert_eq!(of("anc_multi_indel | 21:42750001 | - | G"), 60, "{rows:#?}");
    assert_eq!(of("anc_multi_indel | 21:42750001 | T | -"), 60, "{rows:#?}");
    let on_the_named_transcript: Vec<String> = cells(&lines, &["Allele", "Feature", "AA"])
        .into_iter()
        .filter(|row| row.contains(" | ENST00000291539 | "))
        .collect();
    assert_eq!(
        on_the_named_transcript,
        ["T | ENST00000291539 | -", "- | ENST00000291539 | G"]
    );
}

/// An Ensembl-format record given on the reverse strand scores the row of the
/// substitution it describes: `A/G` on strand -1 at `21:30282172` is `T>C` on
/// the forward strand, which the ClinPred corpus's made-up table scores
/// `0.163359`, and the record's one row carries `ClinPred=0.163359`. Ensembl
/// VEP 116.3 writes `ClinPred=0.375971`, the score of the table's `T G` row at
/// the position: the plugin reads the allele as the record writes it
/// (`ClinPred.pm:109`) and compares it with a row's `alt` by string equality,
/// with no complement (`ClinPred.pm:117`).
#[test]
fn a_reverse_strand_clinpred_record_scores_the_row_of_its_forward_strand_substitution() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "rev_missense.txt",
        "21\t30282172\t30282172\tA/G\t-\trev_missense\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-clinpred-fixture");
    command.args([
        "--assembly",
        "GRCh38",
        "--format",
        "ensembl",
        "--plugin",
        &format!(
            "ClinPred,file={}",
            word(
                &corpus("GRCh38-plugin-clinpred-fixture")
                    .join("data/ClinPred_hg38_sorted_tabbed_synthetic.tsv.gz")
            )
        ),
    ]);
    assert_eq!(
        cells(
            &lines_of(command, &out),
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "Codons",
                "ClinPred",
            ],
        ),
        ["rev_missense | 21:30282172 | G | ENST00000340345 | missense_variant | tAt/tGt | 0.163359"]
    );
}

/// Under `MAX` a two-base substitution whose two scores in the corpus bigWig
/// are `-3.29` and `-4.96` takes `Conservation=-3.290`, the larger of the two,
/// on each of its three rows. Ensembl VEP 116.3 writes `-4.960` there: the
/// plugin orders the covered scores with a bare `sort`, which compares them as
/// text, and reports the last (`Conservation.pm:242-243`), and `-4.96` sorts
/// after `-3.29`.
#[test]
fn a_conservation_maximum_over_negative_scores_is_the_numeric_maximum() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "all_negative_span.vcf",
        "21\t18383679\tall_negative_span\tTT\tAA\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let scores = corpus("GRCh38-plugin-conservation").join("plugin_data/gerp.bw");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-conservation");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("Conservation,{},MAX", word(&scores)),
    ]);
    assert_eq!(
        cells(
            &lines_of(command, &out),
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "Conservation",
            ],
        ),
        [
            "all_negative_span | 21:18383679-18383680 | AA | ENST00000284885 | missense_variant | -3.290",
            "all_negative_span | 21:18383679-18383680 | AA | ENST00000422787 | missense_variant | -3.290",
            "all_negative_span | 21:18383679-18383680 | AA | ENST00000474775 | 5_prime_UTR_variant | -3.290",
        ]
    );
}

/// Four records over regulatory features and motifs, under `--regulatory`
/// with DosageSensitivity loaded. vep-rs writes the three regulatory rows, the
/// four motif rows and the four intergenic rows without `pHaplo` or `pTriplo`
/// and writes no warning. Ensembl VEP 116.3 writes the same rows and two
/// `Plugin 'DosageSensitivity' went wrong` warnings, one per row kind
/// (`OutputFactory.pm:2253`): the plugin reads the transcript of every row it
/// is run on (`DosageSensitivity.pm:169`), and a regulatory or motif row has
/// none.
#[test]
fn dosage_sensitivity_writes_no_warning_on_regulatory_and_motif_rows() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "regulatory.vcf",
        "21\t5014321\tdel_covers_rf\tGTGGGCCTGTGGGGGGCGCTTACC\tG\t.\t.\t.\n\
         21\t5055683\tlonger_not_repeat\tTGCCCCCAGGGAGG\tTGCCCCCAGGGAGGAAA\t.\t.\t.\n\
         21\t5090720\tfwd_mnv_inside\tGG\tAA\t.\t.\t.\n\
         21\t5090723\trev_mnv_start\tCA\tTT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let scores = corpus("GRCh38-plugin-dosage-sensitivity")
        .join("data/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-regulatory");
    command.args([
        "--assembly",
        "GRCh38",
        "--fasta",
        &word(&corpus("GRCh38-regulatory").join("reference.fa.gz")),
        "--regulatory",
        "--plugin",
        &format!("DosageSensitivity,file={}", word(&scores)),
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Uploaded_variation",
            "Allele",
            "Feature",
            "Feature_type",
            "Consequence",
            "pHaplo",
            "pTriplo",
        ],
    );
    assert_eq!(
        rows,
        [
            "del_covers_rf | - | ENSR21_536SRH | RegulatoryFeature | \
             regulatory_region_ablation,regulatory_region_variant | - | -",
            "del_covers_rf | - | ENSR21_BZTD | RegulatoryFeature | \
             regulatory_region_variant | - | -",
            "del_covers_rf | - | ENSM00000574665 | MotifFeature | \
             TFBS_ablation,TF_binding_site_variant | - | -",
            "del_covers_rf | - | - | - | intergenic_variant | - | -",
            "longer_not_repeat | AAA | ENSR21_D8QT | RegulatoryFeature | \
             regulatory_region_variant | - | -",
            "longer_not_repeat | AAA | - | - | intergenic_variant | - | -",
            "fwd_mnv_inside | AA | ENSM00000574814 | MotifFeature | \
             TF_binding_site_variant | - | -",
            "fwd_mnv_inside | AA | - | - | intergenic_variant | - | -",
            "rev_mnv_start | TT | ENSM00000574814 | MotifFeature | \
             TF_binding_site_variant | - | -",
            "rev_mnv_start | TT | ENSM00000574815 | MotifFeature | \
             TF_binding_site_variant | - | -",
            "rev_mnv_start | TT | - | - | intergenic_variant | - | -",
        ]
    );
    let warnings =
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap_or_default();
    assert_eq!(warnings.lines().collect::<Vec<_>>(), Vec::<&str>::new());
}

/// Under `--minimal` an insertion written beside the deletion of the base
/// before it keeps its own coding coordinates, so the `ProteinLengthChange` on
/// each of its 17 rows is the one the same insertion writes as a record of its
/// own. Ensembl VEP 116.3 writes `0` on the insertion's ten frameshift rows in
/// the default format: the plugin reads the rejoined record's coding
/// coordinates (`Downstream.pm:108`) and writes the inserted base in place of
/// the deleted one (`Downstream.pm:110`).
#[test]
fn a_minimal_split_insertion_changes_the_protein_length_as_the_insertion_alone_does() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "ins_beside_del.vcf",
        "21\t36896882\tins_beside_del\tTG\tTGA,T\t.\t.\t.\n\
         21\t36896883\tins_alone\tG\tGA\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-computed");
    command.args([
        "--assembly",
        "GRCh38",
        "--fasta",
        &word(&corpus("GRCh38-plugin-computed").join("reference.fa.gz")),
        "--minimal",
        "--plugin",
        "Downstream",
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Uploaded_variation",
            "Allele",
            "Feature",
            "Consequence",
            "ProteinLengthChange",
        ],
    );
    let of = |record: &str, allele: &str| {
        let prefix = format!("{record} | {allele} | ");
        rows.iter()
            .filter_map(|row| row.strip_prefix(prefix.as_str()))
            .collect::<Vec<&str>>()
    };
    let insertion = [
        "ENST00000336648 | frameshift_variant | -151",
        "ENST00000399120 | frameshift_variant | -151",
        "ENST00000482273 | non_coding_transcript_exon_variant | -",
        "ENST00000612277 | frameshift_variant | -151",
        "ENST00000674895 | frameshift_variant | -151",
        "ENST00000675057 | frameshift_variant,NMD_transcript_variant | -29",
        "ENST00000675307 | frameshift_variant | -151",
        "ENST00000876941 | frameshift_variant | -151",
        "ENST00000993044 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00000993045 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00000993046 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00001017195 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00001017196 | frameshift_variant | -58",
        "ENST00001017197 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00001083022 | frameshift_variant,NMD_transcript_variant | -29",
        "ENST00001111398 | 3_prime_UTR_variant | -",
        "ENST00001111936 | frameshift_variant | 16",
    ];
    let deletion = [
        "ENST00000336648 | frameshift_variant | -242",
        "ENST00000399120 | frameshift_variant | -242",
        "ENST00000482273 | non_coding_transcript_exon_variant | -",
        "ENST00000612277 | frameshift_variant | -242",
        "ENST00000674895 | frameshift_variant | -242",
        "ENST00000675057 | frameshift_variant,NMD_transcript_variant | -120",
        "ENST00000675307 | frameshift_variant | -242",
        "ENST00000876941 | frameshift_variant | -242",
        "ENST00000993044 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00000993045 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00000993046 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00001017195 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00001017196 | frameshift_variant | -149",
        "ENST00001017197 | 3_prime_UTR_variant,NMD_transcript_variant | -",
        "ENST00001083022 | frameshift_variant,NMD_transcript_variant | -120",
        "ENST00001111398 | 3_prime_UTR_variant | -",
        "ENST00001111936 | frameshift_variant | -36",
    ];
    assert_eq!(rows.len(), 51, "{rows:#?}");
    assert_eq!(of("ins_beside_del", "A"), insertion);
    assert_eq!(of("ins_beside_del", "-"), deletion);
    assert_eq!(of("ins_alone", "A"), insertion);
}

/// The lines of the corpus's gnomAD mitochondrial slice, a bgzip-compressed
/// VCF whose contig is named `chrM`.
fn gnomadmt_slice_lines() -> Vec<String> {
    use std::io::Read;

    let slice = corpus("GRCh38-plugin-gnomadmt").join("data/gnomad.v3.1.chrM.slice.vcf.bgz");
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(slice).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines().map(str::to_string).collect()
}

/// A file record whose `vep` value holds `=` inside it, the `p.Leu45=` of a
/// synonymous change, is written whole under `fields=vep`: each of the input
/// record's 21 rows carries the `vep` item of the slice's record at 3441 as
/// the file holds it, 362 characters in 45 columns. Ensembl VEP 116.3 writes
/// the 147 characters before that sign: the plugin reads each item of the INFO
/// column with `my ($key, $value) = split /=/, $field`
/// (`gnomADMt.pm:327-330`), which Perl splits into at most three parts, and
/// the third is discarded.
#[test]
fn a_gnomadmt_info_value_holding_an_equals_sign_is_written_whole() {
    let dir = tempfile::tempdir().unwrap();
    let slice = gnomadmt_slice_lines();
    let items: Vec<&str> = slice
        .iter()
        .filter(|line| line.starts_with("chrM\t3441\t"))
        .flat_map(|line| line.split('\t').nth(7).unwrap().split(';'))
        .filter_map(|item| item.strip_prefix("vep="))
        .collect();
    let [whole] = items[..] else {
        panic!("{items:#?}")
    };
    let cut = "G|synonymous_variant|LOW|MT-ND1|ENSG00000198888|Transcript|ENST00000361390|protein_coding|1/1||ENST00000361390.2:c.135A>G|ENSP00000354687.2:p.Leu45";
    assert!(whole.starts_with(&format!("{cut}=|")), "{whole}");
    assert_eq!(
        (whole.len(), whole.split('|').count()),
        (362, 45),
        "{whole}"
    );
    let input = vcf(
        dir.path(),
        "nd1_3441_A_G.vcf",
        "MT\t3441\tnd1_3441_A_G\tA\tG\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let file = corpus("GRCh38-plugin-gnomadmt").join("data/gnomad.v3.1.chrM.slice.vcf.bgz");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-gnomadmt");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("gnomADMt,file={},fields=vep", word(&file)),
    ]);
    let expected: Vec<String> = [
        ("ENST00000361390", "synonymous_variant"),
        ("ENST00000361453", "upstream_gene_variant"),
        ("ENST00000361624", "upstream_gene_variant"),
        ("ENST00000361739", "upstream_gene_variant"),
        ("ENST00000361851", "upstream_gene_variant"),
        ("ENST00000386347", "downstream_gene_variant"),
        ("ENST00000387314", "downstream_gene_variant"),
        ("ENST00000387342", "downstream_gene_variant"),
        ("ENST00000387347", "downstream_gene_variant"),
        ("ENST00000387365", "upstream_gene_variant"),
        ("ENST00000387372", "downstream_gene_variant"),
        ("ENST00000387377", "upstream_gene_variant"),
        ("ENST00000387382", "upstream_gene_variant"),
        ("ENST00000387392", "downstream_gene_variant"),
        ("ENST00000387400", "downstream_gene_variant"),
        ("ENST00000387405", "downstream_gene_variant"),
        ("ENST00000387409", "downstream_gene_variant"),
        ("ENST00000387416", "downstream_gene_variant"),
        ("ENST00000387419", "upstream_gene_variant"),
        ("ENST00000387421", "upstream_gene_variant"),
        ("ENST00000389680", "downstream_gene_variant"),
    ]
    .iter()
    .map(|(feature, consequence)| format!("nd1_3441_A_G | G | {feature} | {consequence} | {whole}"))
    .collect();
    assert_eq!(
        cells(
            &lines_of(command, &out),
            &[
                "Uploaded_variation",
                "Allele",
                "Feature",
                "Consequence",
                "gnomAD_vep",
            ],
        ),
        expected
    );
}

/// A data file whose contig is named `MT` annotates a record on `MT`: the
/// file is the corpus's slice with its `##contig` line and the first column of
/// its 373 records rewritten from `chrM`, and each of the record's 10 rows takes
/// the five haplogroup fields of the file's record at position 3. Ensembl VEP
/// 116.3 writes none of the five on those rows: the plugin names a record on
/// `MT` to its file as `M` (`gnomADMt.pm:256-261`), and the base class tries
/// that name and `chrM` alone (`BaseVepTabixPlugin.pm:693-707`).
#[test]
fn a_gnomadmt_file_whose_contig_is_named_mt_annotates_a_record_on_mt() {
    use noodles::csi::binning_index::index::header::Builder;

    let dir = tempfile::tempdir().unwrap();
    let lines: Vec<String> = gnomadmt_slice_lines()
        .iter()
        .map(|line| {
            if let Some(rest) = line.strip_prefix("##contig=<ID=chrM,") {
                format!("##contig=<ID=MT,{rest}")
            } else if let Some(rest) = line.strip_prefix("chrM") {
                format!("MT{rest}")
            } else {
                line.clone()
            }
        })
        .collect();
    let (head, records): (Vec<&str>, Vec<&str>) = lines
        .iter()
        .map(String::as_str)
        .partition(|line| line.starts_with('#'));
    assert!(head.contains(&"##contig=<ID=MT,length=16569,assembly=GRCh38_chrM>"));
    assert_eq!(records.len(), 373);
    assert!(records.iter().all(|record| record.starts_with("MT\t")));
    let renamed = common::indexed_lines(
        dir.path(),
        "gnomad_MT.vcf.bgz",
        Builder::vcf().build(),
        &head,
        &records,
        |fields| {
            let position: usize = fields[1].parse().unwrap();
            (position, position + fields[3].len() - 1)
        },
    );
    let input = vcf(
        dir.path(),
        "MT_3_T_C.vcf",
        "MT\t3\tMT_3_T_C\tT\tC\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-gnomadmt");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("gnomADMt,file={}", word(&renamed)),
    ]);
    let ac_hom = "0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|17|1|0|0|0|0|0";
    let ac_het = "0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0";
    let af_hom = "0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|6.76407e-05|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|5.51948e-03|1.65645e-04|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00";
    let af_het = "0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|3.24675e-04|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00|0.00000e+00";
    let an = "2680|1537|868|603|34|282|91|14784|701|934|3144|2732|663|2977|4724|5672|126|1|1298|366|7|393|3080|6037|1234|819|546|12|89";
    let expected: Vec<String> = [
        ("ENST00000361390", "upstream_gene_variant"),
        ("ENST00000361453", "upstream_gene_variant"),
        ("ENST00000386347", "upstream_gene_variant"),
        ("ENST00000387314", "upstream_gene_variant"),
        ("ENST00000387342", "upstream_gene_variant"),
        ("ENST00000387347", "upstream_gene_variant"),
        ("ENST00000387365", "upstream_gene_variant"),
        ("ENST00000387372", "downstream_gene_variant"),
        ("ENST00000387377", "upstream_gene_variant"),
        ("ENST00000389680", "upstream_gene_variant"),
    ]
    .iter()
    .map(|(feature, consequence)| {
        format!(
            "MT_3_T_C | MT:3 | C | {feature} | {consequence} | {ac_hom} | {ac_het} | {af_hom} | {af_het} | {an}"
        )
    })
    .collect();
    assert_eq!(
        cells(
            &lines_of(command, &out),
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "gnomAD_hap_AC_hom",
                "gnomAD_hap_AC_het",
                "gnomAD_hap_AF_hom",
                "gnomAD_hap_AF_het",
                "gnomAD_hap_AN",
            ],
        ),
        expected
    );
}

/// The lines of the gzip-compressed text file at `path`.
fn loeuf_gzip_lines(path: &Path) -> Vec<String> {
    use std::io::Read;

    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines().map(str::to_string).collect()
}

/// A gnomAD constraint table at `<dir>/<name>`: the `head` line, then `rows`,
/// bgzip-compressed, with a tabix index beside it on the table's coordinate
/// columns, the 114th to the 116th, as `tabix -s 114 -b 115 -e 116` writes it
/// (`LOEUF.pm:66`). The first column of such a table is the gene's name, so
/// each row is indexed under the sequence its 114th column names.
fn loeuf_indexed_table(dir: &Path, name: &str, head: &str, rows: &[String]) -> PathBuf {
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::Builder;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::{bgzf, tabix};
    use std::io::Write;

    const SEQUENCE: usize = 113;
    const START: usize = 114;
    const END: usize = 115;
    let path = dir.join(name);
    let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(
        Builder::default()
            .set_reference_sequence_name_index(SEQUENCE)
            .set_start_position_index(START)
            .set_end_position_index(Some(END))
            .build(),
    );
    writeln!(writer, "{head}").unwrap();
    for row in rows {
        let from = writer.virtual_position();
        writeln!(writer, "{row}").unwrap();
        let cells: Vec<&str> = row.split('\t').collect();
        let position = |column: usize| Position::try_from(cells[column].parse::<usize>().unwrap());
        indexer
            .add_record(
                cells[SEQUENCE],
                position(START).unwrap(),
                position(END).unwrap(),
                Chunk::new(from, writer.virtual_position()),
            )
            .unwrap();
    }
    writer.finish().unwrap();
    let mut index = path.clone().into_os_string();
    index.push(".tbi");
    tabix::fs::write(index, &indexer.build()).unwrap();
    path
}

/// Under `match_by=gene` a gene whose smallest LOEUF cell is written `0` takes
/// `LOEUF=0` on every row of that gene. The table is the corpus's own with the
/// `lof.oe_ci.upper` cell of `ENST00000359726`, the minimum of `APP`, rewritten
/// from `2.5883e-01`; the record's 55 rows are all of that gene and each reads
/// `0`. Ensembl VEP 116.3 writes no LOEUF on the 55: the plugin keeps the
/// minimum as the string the table holds, which its numeric test admits
/// (`LOEUF.pm:192`), and returns it only where Perl reads that string as true
/// (`LOEUF.pm:204`), which `0` is not.
#[test]
fn a_loeuf_gene_minimum_of_zero_is_written_on_every_row_of_its_gene() {
    let dir = tempfile::tempdir().unwrap();
    let lines =
        loeuf_gzip_lines(&corpus("GRCh38-plugin-loeuf").join("data/loeuf_dataset_grch38.tsv.gz"));
    let (head, rows) = lines.split_first().unwrap();
    let mut rewritten = 0;
    let rows: Vec<String> = rows
        .iter()
        .map(|row| {
            let mut cells: Vec<&str> = row.split('\t').collect();
            if cells[2] == "ENST00000359726" {
                assert_eq!(cells[93], "2.5883e-01");
                cells[93] = "0";
                rewritten += 1;
            }
            cells.join("\t")
        })
        .collect();
    assert_eq!(rewritten, 1);
    let table = loeuf_indexed_table(dir.path(), "loeuf_zero.tsv.gz", head, &rows);
    let input = vcf(
        dir.path(),
        "app_exon.vcf",
        "21\t26112000\tapp_exon\tG\tA\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-loeuf");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("LOEUF,file={},match_by=gene", word(&table)),
    ]);
    let lines = lines_of(command, &out);
    let rows = cells(
        &lines,
        &[
            "Uploaded_variation",
            "Allele",
            "Gene",
            "Consequence",
            "LOEUF",
        ],
    );
    let of = |consequence: &str| {
        rows.iter()
            .filter(|row| **row == format!("app_exon | A | ENSG00000142192 | {consequence} | 0"))
            .count()
    };
    assert_eq!(rows.len(), 55, "{rows:#?}");
    assert_eq!(of("synonymous_variant"), 44, "{rows:#?}");
    assert_eq!(of("intron_variant"), 4, "{rows:#?}");
    assert_eq!(of("non_coding_transcript_exon_variant"), 4, "{rows:#?}");
    assert_eq!(
        of("synonymous_variant,NMD_transcript_variant"),
        3,
        "{rows:#?}"
    );
    let minimum: Vec<String> = cells(&lines, &["Feature", "Consequence", "LOEUF"])
        .into_iter()
        .filter(|row| row.starts_with("ENST00000359726 | "))
        .collect();
    assert_eq!(minimum, ["ENST00000359726 | synonymous_variant | 0"]);
}

/// The LoFTEE corpus's GERP bigWig, which carries no zoom level.
fn loftee_corpus_gerp() -> PathBuf {
    corpus("GRCh38-plugin-loftee")
        .join("plugin_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw")
}

/// A `vep` run of `input` over the LoFTEE corpus, its reference FASTA and its
/// ancestral FASTA, with the `LoF` plugin reading GERP from the bigWig `gerp`.
/// `--dir_plugins` and `loftee_path`, where Ensembl VEP finds the plugin's
/// module, name the corpus's plugin data.
fn loftee_run(dir: &Path, input: &Path, out: &Path, gerp: &Path) -> Command {
    let corpus = corpus("GRCh38-plugin-loftee");
    let data = corpus.join("plugin_data");
    let mut command = vep(dir, input, out, "GRCh38-plugin-loftee");
    command
        .args(["--assembly", "GRCh38", "--fasta"])
        .arg(corpus.join("reference.fa.gz"))
        .arg("--dir_plugins")
        .arg(&data)
        .args([
            "--plugin",
            &format!(
                "LoF,loftee_path:{data},human_ancestor_fa:{data}/human_ancestor.fa.gz,gerp_bigwig:{gerp}",
                data = word(&data),
                gerp = word(gerp)
            ),
        ]);
    command
}

/// A stop gain whose base is the ancestral one, alone on its record and beside
/// a sibling allele: `21 32314600 C T` and `21 32314600 C T,A` on
/// `ENST00000339944`, where the corpus's ancestral FASTA reads `T`. vep-rs
/// filters the `T` allele `ANC_ALLELE` on both records. Ensembl VEP 116.3
/// filters it on the one-allele record alone, since `LoF.pm:531` returns
/// before the comparison when the record's allele string, here `C/T/A`, is not
/// three characters long. The two records' other 36 rows carry no `LoF` field.
#[test]
fn loftee_filters_an_ancestral_allele_beside_a_sibling_allele() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "ancestral.vcf",
        "21\t32314600\tone\tC\tT\t.\t.\t.\n21\t32314600\ttwo\tC\tT,A\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let command = loftee_run(dir.path(), &input, &out, &loftee_corpus_gerp());
    let lines = lines_of(command, &out);
    let on_the_transcript: Vec<String> = cells(
        &lines,
        &[
            "Uploaded_variation",
            "Location",
            "Allele",
            "Feature",
            "Consequence",
            "LoF",
            "LoF_filter",
        ],
    )
    .into_iter()
    .filter(|row| row.contains(" | ENST00000339944 | "))
    .collect();
    assert_eq!(
        on_the_transcript,
        [
            "one | 21:32314600 | T | ENST00000339944 | stop_gained | LC | END_TRUNC,ANC_ALLELE",
            "two | 21:32314600 | T | ENST00000339944 | stop_gained | LC | END_TRUNC,ANC_ALLELE",
            "two | 21:32314600 | A | ENST00000339944 | synonymous_variant | - | -",
        ]
    );
    let rows = cells(
        &lines,
        &["Uploaded_variation", "Allele", "LoF", "LoF_filter"],
    );
    let of = |kind: &str| rows.iter().filter(|row| *row == kind).count();
    assert_eq!(rows.len(), 39, "{rows:#?}");
    assert_eq!(of("one | T | LC | END_TRUNC,ANC_ALLELE"), 1, "{rows:#?}");
    assert_eq!(of("two | T | LC | END_TRUNC,ANC_ALLELE"), 1, "{rows:#?}");
    assert_eq!(of("one | T | - | -"), 12, "{rows:#?}");
    assert_eq!(of("two | T | - | -"), 12, "{rows:#?}");
    assert_eq!(of("two | A | - | -"), 13, "{rows:#?}");
}

/// The corpus's GC-donor record with a second alternate allele, `21 46126243
/// C T,A`, on the two transcripts whose intron opens `GC` at 46126242. vep-rs
/// tests each allele's own base: the `T` rows, which make the donor `GT`, read
/// `LoF=LC` with `LoF_filter=GC_TO_GT_DONOR`, and the `A` rows, which make it
/// `GA`, read `LoF=HC` with no filter. Ensembl VEP 116.3 writes the filter on
/// the `A` rows too: `LoF.pm:363` takes the reference and the first alternate
/// of the record's allele list for every allele. The record's other 132 rows
/// carry no `LoF` field.
#[test]
fn loftee_judges_a_gc_donor_by_each_allele_of_its_record() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "donor.vcf",
        "21\t46126243\t497375\tC\tT,A\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let command = loftee_run(dir.path(), &input, &out, &loftee_corpus_gerp());
    let lines = lines_of(command, &out);
    let donor: Vec<String> = cells(
        &lines,
        &[
            "Location",
            "Allele",
            "Feature",
            "Consequence",
            "LoF",
            "LoF_filter",
        ],
    )
    .into_iter()
    .filter(|row| row.contains(" | splice_donor_variant | "))
    .collect();
    assert_eq!(
        donor,
        [
            "21:46126243 | T | ENST00001028357 | splice_donor_variant | LC | GC_TO_GT_DONOR",
            "21:46126243 | A | ENST00001028357 | splice_donor_variant | HC | -",
            "21:46126243 | T | ENST00001028364 | splice_donor_variant | LC | GC_TO_GT_DONOR",
            "21:46126243 | A | ENST00001028364 | splice_donor_variant | HC | -",
        ]
    );
    let rows = cells(
        &lines,
        &["Uploaded_variation", "Allele", "LoF", "LoF_filter"],
    );
    let of = |kind: &str| rows.iter().filter(|row| *row == kind).count();
    assert_eq!(rows.len(), 136, "{rows:#?}");
    assert_eq!(of("497375 | T | LC | GC_TO_GT_DONOR"), 2, "{rows:#?}");
    assert_eq!(of("497375 | A | HC | -"), 2, "{rows:#?}");
    assert_eq!(of("497375 | T | - | -"), 66, "{rows:#?}");
    assert_eq!(of("497375 | A | - | -"), 66, "{rows:#?}");
}

/// One record of a bigWig zoom level: the 0-based start of the span's first
/// scored base, the end of its last, the count of its scored bases and their
/// least value, greatest value, sum and sum of squares.
struct LofteeZoomRecord {
    start: u32,
    end: u32,
    scored: u32,
    least: f64,
    greatest: f64,
    sum: f64,
    squares: f64,
}

/// The corpus's GERP bigWig with zoom levels of 10, 40, 160, 640 and 2,560
/// bases, at `<dir>/gerp-zoomed.bw`. The corpus's bytes stand as they are, so
/// the per-base values are its own; each level's records and their index
/// follow them, ahead of the closing magic number, and the header's zoom
/// slots, level count and buffer size name the levels. A level tiles
/// chromosome 21 in spans of its size from the first scored base and holds one
/// record for each span with a scored base.
fn loftee_gerp_with_zoom_levels(dir: &Path) -> PathBuf {
    use std::io::Write as _;

    let source = loftee_corpus_gerp();
    let mut file = std::fs::read(&source).unwrap();
    // A little-endian bigWig that closes on its magic number, counts no zoom
    // level in bytes 6 to 8 and leaves its ten 24-byte zoom slots, bytes 64
    // to 304, empty.
    let magic = file.split_off(file.len() - 4);
    assert_eq!(magic, 0x888F_FC26_u32.to_le_bytes());
    assert_eq!(file[..4], magic[..]);
    assert!(file[6..8] == [0, 0] && file[64..304].iter().all(|byte| *byte == 0));

    let mut reader = vep_builtins::bigwig::BigWigFile::open(&source).unwrap();
    let chromosomes = reader.chromosome_lengths();
    assert_eq!(chromosomes.len(), 1);
    let length = chromosomes[0].1;
    let values: Vec<(u32, u32, f64)> = reader
        .values("21", 0, length)
        .unwrap()
        .into_iter()
        .map(|(start, end, value)| (start, end, f64::from(value)))
        .collect();
    let first = values[0].0;

    // Bytes 52 to 56 hold the size the largest block inflates to.
    let mut buffer = u32::from_le_bytes(file[52..56].try_into().unwrap());
    for (level, size) in [10_u32, 40, 160, 640, 2560].into_iter().enumerate() {
        let mut records: Vec<LofteeZoomRecord> = Vec::new();
        for &(start, end, value) in &values {
            let mut at = start;
            while at < end {
                let span = first + (at - first) / size * size;
                let to = end.min(span + size);
                let bases = f64::from(to - at);
                match records.last_mut() {
                    Some(record) if record.start >= span => {
                        record.end = to;
                        record.scored += to - at;
                        record.least = record.least.min(value);
                        record.greatest = record.greatest.max(value);
                        record.sum += bases * value;
                        record.squares += bases * value * value;
                    }
                    _ => records.push(LofteeZoomRecord {
                        start: at,
                        end: to,
                        scored: to - at,
                        least: value,
                        greatest: value,
                        sum: bases * value,
                        squares: bases * value * value,
                    }),
                }
                at = to;
            }
        }

        // The level's records, 1,024 to a compressed block, each opening with
        // the id of the file's one chromosome, 0, and one 32-byte index entry
        // for each block: its bounds, offset and size.
        let data = file.len() as u64;
        let mut entries = Vec::new();
        for block in records.chunks(1024) {
            let mut raw = Vec::with_capacity(32 * block.len());
            for record in block {
                for word in [0, record.start, record.end, record.scored] {
                    raw.extend(word.to_le_bytes());
                }
                for value in [record.least, record.greatest, record.sum, record.squares] {
                    raw.extend((value as f32).to_le_bytes());
                }
            }
            buffer = buffer.max(raw.len() as u32);
            let mut encoder =
                flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(&raw).unwrap();
            let packed = encoder.finish().unwrap();
            for word in [0, block[0].start, 0, block[block.len() - 1].end] {
                entries.extend(word.to_le_bytes());
            }
            entries.extend((file.len() as u64).to_le_bytes());
            entries.extend((packed.len() as u64).to_le_bytes());
            file.extend(packed);
        }

        // The level's index: the 48-byte R-tree header (magic number, entries
        // to a node, block count, the bounds of the records, the end of the
        // indexed data, records to a block, a reserved word), then the root, a
        // leaf that holds every entry.
        let index = file.len() as u64;
        let blocks = entries.len() / 32;
        assert!(blocks <= 256);
        file.extend(0x2468_ACE0_u32.to_le_bytes());
        file.extend(256_u32.to_le_bytes());
        file.extend((blocks as u64).to_le_bytes());
        for word in [0, records[0].start, 0, records[records.len() - 1].end] {
            file.extend(word.to_le_bytes());
        }
        file.extend(index.to_le_bytes());
        file.extend(1024_u32.to_le_bytes());
        file.extend(0_u32.to_le_bytes());
        file.extend([1, 0]);
        file.extend((blocks as u16).to_le_bytes());
        file.extend(entries);

        // The level's slot: its span size, a reserved word, the offset of its
        // records and the offset of its index.
        let slot = 64 + 24 * level;
        file[slot..slot + 4].copy_from_slice(&size.to_le_bytes());
        file[slot + 8..slot + 16].copy_from_slice(&data.to_le_bytes());
        file[slot + 16..slot + 24].copy_from_slice(&index.to_le_bytes());
    }
    file[6..8].copy_from_slice(&5_u16.to_le_bytes());
    file[52..56].copy_from_slice(&buffer.to_le_bytes());
    file.extend(magic);
    let path = dir.join("gerp-zoomed.bw");
    std::fs::write(&path, file).unwrap();
    path
}

/// The corpus's stop gain in the single exon of CLDN8, `21 30215764 G T`, run
/// with the corpus's GERP values in a bigWig that carries zoom levels. vep-rs
/// sums the per-base values of the plugin's one interval, the 517 bases
/// `21:30215248-30215764`, and writes `GERP_DIST:-204.02089612931`, the total
/// it writes from the corpus's own file. Ensembl VEP 116.3 totals the interval
/// through one summary bin (`gerp_dist.pl:70-79`), which the bigWig library
/// serves from the 160-base level, and writes `GERP_DIST:-283.796467024088`.
/// The test first reads that total from the level as the library does: each
/// span's sum and count of scored bases enter by the share of the span inside
/// the interval, and the total is the interval's length times their mean.
#[test]
fn loftee_sums_gerp_per_base_on_a_bigwig_with_zoom_levels() {
    let dir = tempfile::tempdir().unwrap();
    let zoomed = loftee_gerp_with_zoom_levels(dir.path());
    let mut reader = vep_builtins::bigwig::BigWigFile::open(&zoomed).unwrap();
    let levels: Vec<u32> = reader.zoom_reduction_levels();
    assert_eq!(levels, [10, 40, 160, 640, 2560]);
    let (from, to) = (30_215_247_u32, 30_215_764_u32);
    let (mut sum, mut scored) = (0.0_f64, 0.0_f64);
    for (start, end, bases_covered, record_sum) in reader.zoom_sums("21", from, to, 160).unwrap() {
        let inside = end.min(to).saturating_sub(start.max(from));
        let share = f64::from(inside) / f64::from(end - start);
        sum += record_sum * share;
        scored += bases_covered as f64 * share;
    }
    assert_eq!(
        format!("{:.12}", f64::from(to - from) * (sum / scored)),
        "-283.796467024088"
    );

    let input = vcf(
        dir.path(),
        "cldn8.vcf",
        "21\t30215764\t64605\tG\tT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let command = loftee_run(dir.path(), &input, &out, &zoomed);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Feature",
            "Consequence",
            "LoF",
            "LoF_filter",
            "LoF_flags",
            "LoF_info",
        ],
    );
    assert_eq!(
        rows,
        [
            "ENST00000399899 | stop_gained | LC | END_TRUNC | SINGLE_EXON | PERCENTILE:0.238938053097345,GERP_DIST:-204.02089612931,BP_DIST:516,DIST_FROM_LAST_EXON:-333,50_BP_RULE:FAIL",
            "ENST00000451410 | upstream_gene_variant | - | - | - | -",
        ]
    );
}

/// The made-up MaxEntScan matrix directory that the `generated` block of the
/// corpus `GRCh38-plugin-maxentscan-fixture` writes, built under `dir` and
/// held to the digests of that corpus's manifest.
fn maxentscan_fixture_matrices(dir: &Path) -> String {
    let fixture = common::corpora()
        .into_iter()
        .find(|held| held.release == "116" && held.name == "GRCh38-plugin-maxentscan-fixture")
        .unwrap();
    let words = common::expand_argv(
        &fixture,
        &["{generated}/fordownload_synthetic".to_string()],
        &dir.join("out.txt"),
        dir,
    );
    words.into_iter().next().unwrap()
}

/// A substitution in a two-base exon lies in the acceptor window of the intron
/// before the exon in transcript order and in the donor window of the intron
/// after it, and its row carries the scores of the first of the two windows,
/// the acceptor window, whose 23 bases `verbose` prints. The matrices are the
/// made-up ones the `generated` block of the corpus
/// `GRCh38-plugin-maxentscan-fixture` writes. Ensembl VEP 116.3 scores the
/// donor window of the intron after the exon on this transcript, 9 bases: the
/// plugin takes the first window met on the list of overlapped introns
/// (`MaxEntScan.pm:270`), whose order is an interval tree's
/// (`BaseTranscriptVariation.pm:813`).
#[test]
fn a_substitution_in_two_maxentscan_windows_is_scored_in_the_acceptor_window_before_its_exon() {
    let dir = tempfile::tempdir().unwrap();
    let matrices = maxentscan_fixture_matrices(dir.path());
    let input = vcf(
        dir.path(),
        "fly_two_base_exon.vcf",
        "2L\t21019776\tfly_two_base_exon\tA\tG\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "BDGP6.54-fly");
    command.args([
        "--assembly",
        "BDGP6.54",
        "--species",
        "drosophila_melanogaster",
        "--fasta",
        &word(&corpus("BDGP6.54-fly").join("reference.fa.gz")),
        "--plugin",
        &format!("MaxEntScan,{matrices},verbose"),
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Uploaded_variation",
            "Location",
            "Allele",
            "Feature",
            "Consequence",
            "MaxEntScan_ref_seq",
            "MaxEntScan_alt_seq",
            "MaxEntScan_ref",
            "MaxEntScan_alt",
            "MaxEntScan_diff",
        ],
    );
    assert_eq!(
        rows,
        [
            "fly_two_base_exon | 2L:21019776 | G | FBtr0290305 | intron_variant | - | - | - | - | -",
            "fly_two_base_exon | 2L:21019776 | G | FBtr0331199 | splice_region_variant,5_prime_UTR_variant | GGTGAGAGAGTGAGGAAGAGATG | GGTGAGAGAGTGAGGAAGAGACG | 4.428 | -4.384 | 8.812",
            "fly_two_base_exon | 2L:21019776 | G | FBtr0347554 | intron_variant | - | - | - | - | -",
        ]
    );
}

/// A substitution nine bases before the first exon of a transcript with two
/// two-base introns reports that exon on its upstream row, whatever records
/// stand around it; here the record inside the first short intron follows it.
/// Ensembl VEP 116.3 writes no `NearestExonJB` on that row in this order: the
/// plugin reads the transcript's overlapped exons (`NearestExonJB.pm:184`),
/// which are stretched by 12 bases only under a flag
/// (`BaseTranscriptVariation.pm:858-865`) that the transcript carries once a
/// record inside it has built its intron index
/// (`BaseTranscriptVariation.pm:1032-1034`).
#[test]
fn a_record_just_outside_a_short_intron_transcript_reports_its_nearest_exon_in_any_order() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "nine_before_first_exon.vcf",
        "21\t6354255\tnine_before_first_exon\tG\tT\t.\t.\t.\n\
         21\t6358575\tin_first_short_intron\tG\tT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-ancestral");
    command.args(["--assembly", "GRCh38", "--plugin", "NearestExonJB"]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Uploaded_variation",
            "Location",
            "Allele",
            "Feature",
            "Consequence",
            "DISTANCE",
            "NearestExonJB",
        ],
    );
    let (bare, valued): (Vec<&str>, Vec<&str>) = rows
        .iter()
        .map(String::as_str)
        .partition(|row| row.ends_with(" | -"));
    let bare_of = |record: &str| bare.iter().filter(|row| row.starts_with(record)).count();
    assert_eq!(rows.len(), 54, "{rows:#?}");
    assert_eq!(
        valued,
        [
            "nine_before_first_exon | 21:6354255 | T | ENST00000624847 | upstream_gene_variant | 9 | ENSE00003757672|9|start|259",
            "in_first_short_intron | 21:6358575 | T | ENST00000624165 | non_coding_transcript_exon_variant | - | ENSE00003760095|20|start|630",
            "in_first_short_intron | 21:6358575 | T | ENST00000624576 | non_coding_transcript_exon_variant | - | ENSE00003755215|20|start|1861",
            "in_first_short_intron | 21:6358575 | T | ENST00000624847 | non_coding_transcript_variant | - | ENSE00003758864|1|end|20",
        ]
    );
    assert_eq!(
        bare_of("nine_before_first_exon | 21:6354255 | T | "),
        44,
        "{rows:#?}"
    );
    assert_eq!(
        bare_of("in_first_short_intron | 21:6358575 | T | "),
        6,
        "{rows:#?}"
    );
}

/// Two alternate alleles of one record that give the same amino acid each take
/// the score of the table's row for their own allele. The record's `C` and `T`
/// are both `M/L` on `ENST00000299443`, where the corpus's table scores `A C`
/// 0.025 and `A T` 0.027; the table names no other transcript at the position,
/// so the rows of `ENST00000620442` carry no score. Ensembl VEP 116.3 writes
/// 0.025 on the `T` row as well: the plugin matches a table row against every
/// alternate allele of the record (`REVEL.pm:165`) and returns the first whose
/// alternate amino acid is the row's (`REVEL.pm:193`).
#[test]
fn a_revel_score_is_read_from_the_table_row_of_its_own_allele() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "two_leucines.vcf",
        "21\t13610388\ttwo_leucines\tA\tC,T\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let table = corpus("GRCh38-plugin-revel").join("plugin_data/revel.tsv.gz");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-revel");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("REVEL,file={}", word(&table)),
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Location",
            "Allele",
            "Feature",
            "Consequence",
            "Amino_acids",
            "REVEL",
        ],
    );
    assert_eq!(
        rows,
        [
            "21:13610388 | C | ENST00000299443 | missense_variant | M/L | 0.025",
            "21:13610388 | T | ENST00000299443 | missense_variant | M/L | 0.027",
            "21:13610388 | C | ENST00000620442 | missense_variant | M/L | -",
            "21:13610388 | T | ENST00000620442 | missense_variant | M/L | -",
        ]
    );
}

/// A record whose alternate alleles all share the reference's first base and
/// mix a deletion with a substitution is read with that base removed, here as
/// `C/-/A` one position on. The substitution's six rows score the `C A` row
/// the corpus's synthetic SNV file holds at that position, and the deletion's
/// six rows carry nothing, the synthetic indel file holding no row for it.
/// Ensembl VEP 116.3 writes no SpliceAI field on any of the twelve rows: seeing
/// the `-`, the plugin rebuilds the whole record with its anchor base
/// (`SpliceAI.pm:259-266`) and compares the substitution as `GC GA` at the
/// anchor's position, which the score row does not equal (`SpliceAI.pm:272`).
#[test]
fn a_spliceai_substitution_beside_a_deletion_in_one_record_scores_its_own_row() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "chopped_multiallelic.vcf",
        "21\t26965499\tchopped_multiallelic\tGC\tG,GA\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let fixture = corpus("GRCh38-plugin-spliceai-fixture");
    let scores = fixture.join("plugin_data");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-spliceai-fixture");
    command.args([
        "--assembly",
        "GRCh38",
        "--fasta",
        &word(&fixture.join("reference.fa.gz")),
        "--plugin",
        &format!(
            "SpliceAI,snv={},indel={}",
            word(&scores.join("spliceai_scores.raw.snv.hg38.synthetic.vcf.gz")),
            word(&scores.join("spliceai_scores.raw.indel.hg38.synthetic.vcf.gz")),
        ),
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &["Location", "Allele", "Feature", "SpliceAI_pred"],
    );
    assert_eq!(
        rows,
        [
            "21:26965500 | - | ENST00000284987 | -",
            "21:26965500 | A | ENST00000284987 | ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17",
            "21:26965500 | - | ENST00000652031 | -",
            "21:26965500 | A | ENST00000652031 | ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17",
            "21:26965500 | - | ENST00000970346 | -",
            "21:26965500 | A | ENST00000970346 | ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17",
            "21:26965500 | - | ENST00001108182 | -",
            "21:26965500 | A | ENST00001108182 | ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17",
            "21:26965500 | - | ENST00001108183 | -",
            "21:26965500 | A | ENST00001108183 | ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17",
            "21:26965500 | - | ENST00001131497 | -",
            "21:26965500 | A | ENST00001131497 | ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17",
        ]
    );
}

/// An Ensembl-format insertion given on the reverse strand scores the row its
/// forward-strand twin scores: `-/G` on strand -1 between `21:21338600` and
/// `21:21338601` is the forward-strand insertion of `C`, the `A AC` row of the
/// SpliceAI corpus's made-up indel table, whose scores each of the record's
/// four rows carries. Ensembl VEP 116.3 writes no `SpliceAI_pred` on the four
/// rows: the plugin reads the forward-strand alleles `to_VCF_record` writes
/// and reverse-complements them again (`SpliceAI.pm:259-271`), so the alleles
/// it compares equal no row of the table.
#[test]
fn a_reverse_strand_spliceai_insertion_scores_the_row_of_its_forward_strand_twin() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "rev_ins.txt",
        "21\t21338601\t21338600\t-/G\t-\trev_ins\n",
    );
    let out = dir.path().join("out.txt");
    let fixture = corpus("GRCh38-plugin-spliceai-fixture");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-spliceai-fixture");
    command.args([
        "--assembly",
        "GRCh38",
        "--format",
        "ensembl",
        "--fasta",
        &word(&fixture.join("reference.fa.gz")),
        "--plugin",
        &format!(
            "SpliceAI,snv={},indel={}",
            word(&fixture.join("plugin_data/spliceai_scores.raw.snv.hg38.synthetic.vcf.gz")),
            word(&fixture.join("plugin_data/spliceai_scores.raw.indel.hg38.synthetic.vcf.gz"))
        ),
    ]);
    assert_eq!(
        cells(
            &lines_of(command, &out),
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "SpliceAI_pred",
            ],
        ),
        [
            "rev_ins | 21:21338600-21338601 | G | ENST00000284894 | intron_variant | NCAM2|0.91|0.48|0.98|0.40|48|-1|11|17",
            "rev_ins | 21:21338600-21338601 | G | ENST00000400546 | intron_variant | NCAM2|0.91|0.48|0.98|0.40|48|-1|11|17",
            "rev_ins | 21:21338600-21338601 | G | ENST00000461281 | downstream_gene_variant | NCAM2|0.91|0.48|0.98|0.40|48|-1|11|17",
            "rev_ins | 21:21338600-21338601 | G | ENST00001131746 | intron_variant | NCAM2|0.91|0.48|0.98|0.40|48|-1|11|17",
        ]
    );
}

/// Each allele of a two-allele record takes the row the corpus's synthetic
/// table holds for that allele and transcript. On `ENST00000299443` the `C`
/// row carries the `T C` row's delta of 0.51 and its site and the `G` row the
/// `T G` row's 0.24 and its site, and each allele's row on `ENST00000620442`
/// its own in the same way; the table holds no row for `ENST00000384745`,
/// whose rows carry no field. Ensembl VEP 116.3 writes the `T C` row's values
/// on the `G` row of `ENST00000299443`: the plugin takes the record's alternate
/// alleles as a list (`SpliceVault.pm:156`) and returns the first table row of
/// the transcript that matches any of them.
#[test]
fn a_splicevault_row_is_written_on_the_allele_it_names() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "sibling.vcf",
        "21\t13618457\tsibling\tT\tC,G\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let table = corpus("GRCh38-plugin-splicevault-fixture")
        .join("data/SpliceVault_data_GRCh38_synthetic.tsv.gz");
    let mut command = vep(
        dir.path(),
        &input,
        &out,
        "GRCh38-plugin-splicevault-fixture",
    );
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("SpliceVault,file={}", word(&table)),
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Allele",
            "Feature",
            "Consequence",
            "SpliceVault_SpliceAI_delta",
            "SpliceVault_site_pos",
        ],
    );
    assert_eq!(
        rows,
        [
            "C | ENST00000299443 | splice_donor_variant | 0.51 | chr21:13618451",
            "G | ENST00000299443 | splice_donor_variant | 0.24 | chr21:13618443",
            "C | ENST00000384745 | upstream_gene_variant | - | -",
            "G | ENST00000384745 | upstream_gene_variant | - | -",
            "C | ENST00000620442 | splice_donor_variant | 0.92 | chr21:13618448",
            "G | ENST00000620442 | splice_donor_variant | 0.24 | chr21:13618437",
        ]
    );
}

/// An allele the table holds no row for takes no field beside a sibling allele
/// that has one. The corpus's table holds one row at the site, for the `T`
/// allele on `ENST00000618007`: its three columns are written on the `T` row
/// of that transcript and on none of the record's other seventeen rows.
/// Ensembl VEP 116.3 writes them on the `A` row of that transcript as well: the
/// plugin takes the record's alternate alleles as a list
/// (`TranscriptAnnotator.pm:215`) and returns the first table row of the
/// transcript that matches any of them (`TranscriptAnnotator.pm:237`).
#[test]
fn a_transcript_annotator_row_is_written_on_the_allele_it_names() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "rs463608.vcf",
        "21\t10521702\trs463608\tC\tT,A\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let table =
        corpus("GRCh38-plugin-transcript-annotator").join("data/transcript_annotator.tsv.gz");
    let mut command = vep(
        dir.path(),
        &input,
        &out,
        "GRCh38-plugin-transcript-annotator",
    );
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("TranscriptAnnotator,file={}", word(&table)),
    ]);
    let rows = cells(
        &lines_of(command, &out),
        &[
            "Allele",
            "Feature",
            "transcript_annotator_SIFT_score",
            "transcript_annotator_SIFT_pred",
            "transcript_annotator_Comment",
        ],
    );
    assert_eq!(
        rows,
        [
            "T | ENST00000427445 | - | - | -",
            "A | ENST00000427445 | - | - | -",
            "T | ENST00000612746 | - | - | -",
            "A | ENST00000612746 | - | - | -",
            "T | ENST00000612957 | - | - | -",
            "A | ENST00000612957 | - | - | -",
            "T | ENST00000618007 | 0.61 | Tolerated | Intronic",
            "A | ENST00000618007 | - | - | -",
            "T | ENST00000622113 | - | - | -",
            "A | ENST00000622113 | - | - | -",
            "T | ENST00001091160 | - | - | -",
            "A | ENST00001091160 | - | - | -",
            "T | ENST00001091161 | - | - | -",
            "A | ENST00001091161 | - | - | -",
            "T | ENST00001091162 | - | - | -",
            "A | ENST00001091162 | - | - | -",
            "T | ENST00001131865 | - | - | -",
            "A | ENST00001131865 | - | - | -",
        ]
    );
}

/// A TranscriptAnnotator table with a record at position 1 of the contig
/// named `1`: the corpus's own table with that record ahead of its rows.
/// vep-rs names the three fields by the table's last header line and writes
/// them on the one row whose transcript the table names. Ensembl VEP 116.3
/// names them `oneone_0.1`, `oneone_Bogus` and `oneone_RecordAtOneOne`, the
/// record's cells: the plugin reads the column names from the last line that
/// `tabix -h` prints for that one-base region (`TranscriptAnnotator.pm:99-105`),
/// which is the record.
#[test]
fn transcript_annotator_names_its_fields_by_the_header_of_a_table_with_a_record_at_1_1() {
    use noodles::csi::binning_index::index::header::Builder;

    let dir = tempfile::tempdir().unwrap();
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut noodles::bgzf::io::Reader::new(
            std::fs::File::open(
                corpus("GRCh38-plugin-transcript-annotator")
                    .join("data/transcript_annotator.tsv.gz"),
            )
            .unwrap(),
        ),
        &mut text,
    )
    .unwrap();
    let (head, mut records): (Vec<&str>, Vec<&str>) =
        text.lines().partition(|line| line.starts_with('#'));
    assert_eq!(head.len(), 2, "{head:#?}");
    records.insert(0, "1\t1\tN\tA\tENST00000000001\t0.1\tBogus\tRecordAtOneOne");
    let table = common::indexed_lines(
        dir.path(),
        "oneone.tsv.gz",
        Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(1))
            .build(),
        &head,
        &records,
        |fields| {
            let position: usize = fields[1].trim().parse().unwrap();
            (position, position)
        },
    );
    let input = vcf(
        dir.path(),
        "rs463608.vcf",
        "21\t10521702\trs463608\tC\tT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(
        dir.path(),
        &input,
        &out,
        "GRCh38-plugin-transcript-annotator",
    );
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("TranscriptAnnotator,file={}", word(&table)),
    ]);
    let lines = lines_of(command, &out);
    let header: Vec<&str> = lines
        .iter()
        .filter(|line| line.starts_with("## oneone_"))
        .map(String::as_str)
        .collect();
    let column = format!("column from {}", word(&table));
    assert_eq!(
        header,
        [
            format!("## oneone_Comment : {column}"),
            format!("## oneone_SIFT_pred : {column}"),
            format!("## oneone_SIFT_score : {column}"),
        ]
    );
    assert_eq!(
        cells(
            &lines,
            &[
                "Uploaded_variation",
                "Allele",
                "Feature",
                "oneone_SIFT_score",
                "oneone_SIFT_pred",
                "oneone_Comment",
            ],
        ),
        [
            "rs463608 | T | ENST00000427445 | - | - | -",
            "rs463608 | T | ENST00000612746 | - | - | -",
            "rs463608 | T | ENST00000612957 | - | - | -",
            "rs463608 | T | ENST00000618007 | 0.61 | Tolerated | Intronic",
            "rs463608 | T | ENST00000622113 | - | - | -",
            "rs463608 | T | ENST00001091160 | - | - | -",
            "rs463608 | T | ENST00001091161 | - | - | -",
            "rs463608 | T | ENST00001091162 | - | - | -",
            "rs463608 | T | ENST00001131865 | - | - | -",
        ]
    );
}

/// OpenTargets under `cols=gwasGeneI.`, a name that holds a regular-expression
/// metacharacter and names no column of the corpus's export. vep-rs compares
/// the name with the column names as text and refuses the plugin: the 17 rows
/// of the record carry no `OpenTargets_gwasGeneI.`, the header has no
/// `## OpenTargets_` line and the warnings file names the column. Ensembl VEP
/// 116.3 writes the field with one `NA` per matched association on every row:
/// it tests the name as the pattern `/^gwasGeneI.$/`, which matches
/// `gwasGeneId` (`OpenTargets.pm:224`).
#[test]
fn open_targets_refuses_a_cols_name_that_matches_a_column_as_a_pattern_alone() {
    let dir = tempfile::tempdir().unwrap();
    let export = corpus("GRCh38-plugin-opentargets").join("open_targets_vep.slice.tsv.gz");
    let input = vcf(
        dir.path(),
        "lead.vcf",
        "21\t14824378\tot_gwas_and_eqtl_leads\tC\tT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-opentargets");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("OpenTargets,file={},cols=gwasGeneI.", word(&export)),
    ]);
    let lines = lines_of(command, &out);
    let header: Vec<&str> = lines
        .iter()
        .filter(|line| line.starts_with("## OpenTargets_"))
        .map(String::as_str)
        .collect();
    assert_eq!(header, Vec::<&str>::new());
    assert_eq!(
        cells(
            &lines,
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "OpenTargets_gwasGeneI.",
            ],
        ),
        [
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000432230 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000454795 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000625278 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000627623 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000630211 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000630354 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000715450 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000715451 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000715452 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000715453 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000715910 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000794242 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000794243 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000794245 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000794246 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000794247 | -",
            "ot_gwas_and_eqtl_leads | 21:14824378 | T | ENST00000794248 | -",
        ]
    );
    let warnings =
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap_or_default();
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [format!(
            "WARNING: Failed to instantiate plugin OpenTargets: columns not found in \
             the header of {}: gwasGeneI.",
            word(&export)
        )]
    );
}

/// AlphaMissense under `cols=am_pathogenicit.`, a name that holds a
/// regular-expression metacharacter and names no column of the corpus's
/// table. vep-rs compares the prefixed name with the column names as text and
/// refuses the plugin: the four missense rows of the record carry no
/// `am_pathogenicit.`, the header has no `## am_` line and the warnings file
/// lists the table's six columns. Ensembl VEP 116.3 writes the field with no
/// value on each of the four rows: it tests the name as the pattern
/// `/^am_pathogenicit.$/`, which matches `am_pathogenicity`
/// (`AlphaMissense.pm:158`).
#[test]
fn alpha_missense_refuses_a_cols_name_that_matches_a_column_as_a_pattern_alone() {
    let dir = tempfile::tempdir().unwrap();
    let table = corpus("GRCh38-plugin-alphamissense").join("plugin_data/alphamissense.tsv.gz");
    let input = vcf(
        dir.path(),
        "two_alts.vcf",
        "21\t13610379\ttwo_alts_13610379\tG\tA,T\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(dir.path(), &input, &out, "GRCh38-plugin-alphamissense");
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("AlphaMissense,file={},cols=am_pathogenicit.", word(&table)),
    ]);
    let lines = lines_of(command, &out);
    let header: Vec<&str> = lines
        .iter()
        .filter(|line| line.starts_with("## am_"))
        .map(String::as_str)
        .collect();
    assert_eq!(header, Vec::<&str>::new());
    assert_eq!(
        cells(
            &lines,
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "Consequence",
                "am_pathogenicit.",
            ],
        ),
        [
            "two_alts_13610379 | 21:13610379 | A | ENST00000299443 | missense_variant | -",
            "two_alts_13610379 | 21:13610379 | T | ENST00000299443 | missense_variant | -",
            "two_alts_13610379 | 21:13610379 | A | ENST00000620442 | missense_variant | -",
            "two_alts_13610379 | 21:13610379 | T | ENST00000620442 | missense_variant | -",
        ]
    );
    let warnings =
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap_or_default();
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [
            "WARNING: Failed to instantiate plugin AlphaMissense: the following \
          columns were not found in file header: am_pathogenicit.; valid columns are: am_genome, \
          am_uniprot_id, am_transcript_id, am_protein_variant, am_pathogenicity, am_class"
        ]
    );
}

/// TranscriptAnnotator under `cols=SIFT_.core`, a name that holds a
/// regular-expression metacharacter and names no column of the corpus's
/// table. vep-rs compares the name with the column names as text, is left
/// with no valid name and refuses the plugin: the nine rows of the record
/// carry no `transcript_annotator_SIFT_.core`, the header has no
/// `## transcript_annotator_` line and the warnings file lists the table's
/// three columns. Ensembl VEP 116.3 writes the field with no value on the row
/// of `ENST00000618007`, the transcript the table names: it tests the name as
/// the pattern `/^SIFT_.core$/`, which matches `SIFT_score`
/// (`TranscriptAnnotator.pm:136`).
#[test]
fn transcript_annotator_refuses_a_cols_name_that_matches_a_column_as_a_pattern_alone() {
    let dir = tempfile::tempdir().unwrap();
    let table =
        corpus("GRCh38-plugin-transcript-annotator").join("data/transcript_annotator.tsv.gz");
    let input = vcf(
        dir.path(),
        "rs463608.vcf",
        "21\t10521702\trs463608\tC\tT\t.\t.\t.\n",
    );
    let out = dir.path().join("out.txt");
    let mut command = vep(
        dir.path(),
        &input,
        &out,
        "GRCh38-plugin-transcript-annotator",
    );
    command.args([
        "--assembly",
        "GRCh38",
        "--plugin",
        &format!("TranscriptAnnotator,file={},cols=SIFT_.core", word(&table)),
    ]);
    let lines = lines_of(command, &out);
    let header: Vec<&str> = lines
        .iter()
        .filter(|line| line.starts_with("## transcript_annotator_"))
        .map(String::as_str)
        .collect();
    assert_eq!(header, Vec::<&str>::new());
    assert_eq!(
        cells(
            &lines,
            &[
                "Uploaded_variation",
                "Location",
                "Allele",
                "Feature",
                "transcript_annotator_SIFT_.core",
            ],
        ),
        [
            "rs463608 | 21:10521702 | T | ENST00000427445 | -",
            "rs463608 | 21:10521702 | T | ENST00000612746 | -",
            "rs463608 | 21:10521702 | T | ENST00000612957 | -",
            "rs463608 | 21:10521702 | T | ENST00000618007 | -",
            "rs463608 | 21:10521702 | T | ENST00000622113 | -",
            "rs463608 | 21:10521702 | T | ENST00001091160 | -",
            "rs463608 | 21:10521702 | T | ENST00001091161 | -",
            "rs463608 | 21:10521702 | T | ENST00001091162 | -",
            "rs463608 | 21:10521702 | T | ENST00001131865 | -",
        ]
    );
    let warnings =
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap_or_default();
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [
            "WARNING: Failed to instantiate plugin TranscriptAnnotator: no \
          valid columns selected. Valid column names: SIFT_score, SIFT_pred, Comment"
        ]
    );
}
