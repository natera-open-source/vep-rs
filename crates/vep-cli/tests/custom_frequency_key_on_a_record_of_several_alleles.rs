// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a frequency column and lists a record of the
//! whole input record, on the rows the binary writes for records of two
//! alternate alleles: lines 9 and 10 of the custom corpus's input, and a twin
//! of line 9 with its alternate alleles in the other order. The corpus's
//! `custom/sites.vcf.gz` holds one record at each site, read as an `overlap`
//! source named `gnomADg` that asks for `AF`, so its value stands under
//! `gnomADg_AF` on every row of the record. The cache of the
//! `GRCh38-everything` corpus holds each record's co-located variant: at line
//! 9 a gnomAD genomes frequency for the allele `C` alone, at line 10 one for
//! each allele.
//!
//! Each row lists the source's value, then the frequency of its own allele.
//! Ensembl VEP 116.3 lists every allele's frequency on every row of the
//! record (the documented difference
//! `custom_key_list_shared_across_rows`): the rows' hashes are
//! shallow copies of the record's (OutputFactory.pm:391), so they share the
//! list the source's record opened (:1003-1014), each row merges its allele's
//! frequency into it (:1245-1248), and the lines are written once every row is
//! built (:266-272). Every other cell is the one Ensembl VEP writes, with the
//! key listed once where Ensembl VEP lists it twice
//! (`duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based lines of the two corpus records in the custom corpus's
/// `variants.vcf`, each with the start of its line.
const RECORD_LINES: [(usize, &str); 2] = [(9, "21\t5380280\t"), (10, "21\t8987421\t")];

/// The twin of the first record: its alternate alleles in the other order.
const TWIN: &str = "21\t5380280\talleles_reversed\tA\tG,C\t.\t.\t.";

/// The rows of the default format under `--af_gnomadg`.
const DEFAULT_ON: &[&str] = &[
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;gnomADg_AF=0.004,7.994e-06;gnomADg_AFR_AF=0;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=1.993e-05;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;gnomADg_AF=0.004;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS",
    "alleles_reversed\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;gnomADg_AF=0.004;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS",
    "alleles_reversed\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;gnomADg_AF=0.004,7.994e-06;gnomADg_AFR_AF=0;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=1.993e-05;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS",
    "rs766964238\t21:8987421\tC\tENSG00000264063\tENST00000577708\tTranscript\tmature_miRNA_variant\t52\t-\t-\t-\t-\trs766964238\tIMPACT=MODIFIER;STRAND=1;gnomADg_AF=0.005,3.593e-05;gnomADg_AFR_AF=2.801e-05;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=8.973e-05;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0;gnomADg=overlap_2;gnomADg_CLASS=overlap;gnomADg_FILTER=PASS",
    "rs766964238\t21:8987421\tT\tENSG00000264063\tENST00000577708\tTranscript\tmature_miRNA_variant\t52\t-\t-\t-\t-\trs766964238\tIMPACT=MODIFIER;STRAND=1;gnomADg_AF=0.005,1.794e-05;gnomADg_AFR_AF=0;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=0;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0.001082;gnomADg=overlap_2;gnomADg_CLASS=overlap;gnomADg_FILTER=PASS",
    "rs766964238\t21:8987421\tC\tENSG00000264462\tENST00000581792\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs766964238\tIMPACT=MODIFIER;DISTANCE=243;STRAND=1;gnomADg_AF=0.005,3.593e-05;gnomADg_AFR_AF=2.801e-05;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=8.973e-05;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0;gnomADg=overlap_2;gnomADg_CLASS=overlap;gnomADg_FILTER=PASS",
    "rs766964238\t21:8987421\tT\tENSG00000264462\tENST00000581792\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs766964238\tIMPACT=MODIFIER;DISTANCE=243;STRAND=1;gnomADg_AF=0.005,1.794e-05;gnomADg_AFR_AF=0;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=0;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0.001082;gnomADg=overlap_2;gnomADg_CLASS=overlap;gnomADg_FILTER=PASS",
];

/// The column line and the rows of the tab format under `--af_gnomadg`.
const TAB_ON: &[&str] = &[
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_REMAINING_AF\tgnomADg_SAS_AF\tCLIN_SIG\tSOMATIC\tPHENO\tgnomADg\tgnomADg_CLASS",
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004,7.994e-06\t0\t0\t0\t0\t0\t0\t0\t1.993e-05\t0\t0\t-\t-\t-\texact_2\texact",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\texact_2\texact",
    "alleles_reversed\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\texact_2\texact",
    "alleles_reversed\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004,7.994e-06\t0\t0\t0\t0\t0\t0\t0\t1.993e-05\t0\t0\t-\t-\t-\texact_2\texact",
    "rs766964238\t21:8987421\tC\tENSG00000264063\tENST00000577708\tTranscript\tmature_miRNA_variant\t52\t-\t-\t-\t-\trs766964238\tMODIFIER\t-\t1\t-\t-\t0.005,3.593e-05\t2.801e-05\t0\t0\t0\t0\t0\t0\t8.973e-05\t0\t0\t-\t-\t-\toverlap_2\toverlap",
    "rs766964238\t21:8987421\tT\tENSG00000264063\tENST00000577708\tTranscript\tmature_miRNA_variant\t52\t-\t-\t-\t-\trs766964238\tMODIFIER\t-\t1\t-\t-\t0.005,1.794e-05\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0.001082\t-\t-\t-\toverlap_2\toverlap",
    "rs766964238\t21:8987421\tC\tENSG00000264462\tENST00000581792\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs766964238\tMODIFIER\t243\t1\t-\t-\t0.005,3.593e-05\t2.801e-05\t0\t0\t0\t0\t0\t0\t8.973e-05\t0\t0\t-\t-\t-\toverlap_2\toverlap",
    "rs766964238\t21:8987421\tT\tENSG00000264462\tENST00000581792\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs766964238\tMODIFIER\t243\t1\t-\t-\t0.005,1.794e-05\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0.001082\t-\t-\t-\toverlap_2\toverlap",
];

/// The `CSQ` format and the records' lines of the VCF format under `--af_gnomadg`.
const VCF_ON: &[&str] = &[
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_REMAINING_AF|gnomADg_SAS_AF|CLIN_SIG|SOMATIC|PHENO|gnomADg|gnomADg_CLASS",
    "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||0.004&7.994e-06|0|0|0|0|0|0|0|1.993e-05|0|0||||exact_2|exact,G|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||0.004||||||||||||||exact_2|exact",
    "21\t5380280\talleles_reversed\tA\tG,C\t.\t.\tCSQ=G|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||0.004||||||||||||||exact_2|exact,C|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||0.004&7.994e-06|0|0|0|0|0|0|0|1.993e-05|0|0||||exact_2|exact",
    "21\t8987421\trs766964238\tG\tC,T\t.\t.\tCSQ=C|mature_miRNA_variant|MODIFIER||ENSG00000264063|Transcript|ENST00000577708|miRNA|1/1||||52|||||rs766964238||1|||||0.005&3.593e-05|2.801e-05|0|0|0|0|0|0|8.973e-05|0|0||||overlap_2|overlap,T|mature_miRNA_variant|MODIFIER||ENSG00000264063|Transcript|ENST00000577708|miRNA|1/1||||52|||||rs766964238||1|||||0.005&1.794e-05|0|0|0|0|0|0|0|0|0|0.001082||||overlap_2|overlap,C|downstream_gene_variant|MODIFIER|MIR3648-2|ENSG00000264462|Transcript|ENST00000581792|miRNA||||||||||rs766964238|243|1||HGNC|HGNC:50843||0.005&3.593e-05|2.801e-05|0|0|0|0|0|0|8.973e-05|0|0||||overlap_2|overlap,T|downstream_gene_variant|MODIFIER|MIR3648-2|ENSG00000264462|Transcript|ENST00000581792|miRNA||||||||||rs766964238|243|1||HGNC|HGNC:50843||0.005&1.794e-05|0|0|0|0|0|0|0|0|0|0.001082||||overlap_2|overlap",
];

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The three records, written under `dir`: the header lines of the custom
/// corpus's input, the first corpus record, its twin and the second corpus
/// record.
fn records_file(dir: &Path) -> PathBuf {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let [first, second] = RECORD_LINES.map(|(line, start)| {
        let record = input.lines().nth(line - 1).expect("the record's line");
        assert!(
            record.starts_with(start),
            "line {line} of the custom corpus's input is the record"
        );
        record
    });
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    for record in [first, TWIN, second] {
        text.push_str(record);
        text.push('\n');
    }
    let path = dir.join("records.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// The output of a run over the three records in `format` (`default`, `tab`,
/// `vcf` or `parquet`) under `--af_gnomadg`, with the custom corpus's site
/// file as an `overlap` source named `gnomADg` that asks for its `AF` and
/// `CLASS` fields.
fn annotate(dir: &Path, format: &str) -> PathBuf {
    let out = dir.join(format!("{format}.out"));
    let source = format!(
        "file={},short_name=gnomADg,format=vcf,type=overlap,fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display()
    );
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(records_file(dir))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-everything").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(["--custom", &source, "--af_gnomadg"]);
    match format {
        "default" => {}
        "parquet" => {
            cmd.args(["--output_format", "parquet"]);
        }
        other => {
            cmd.arg(format!("--{other}"));
        }
    }
    let run = cmd.output().unwrap();
    assert!(
        run.status.success(),
        "vep {format}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// The cells of the tab column `gnomADg_AF`, one per row of [`TAB_ON`].
fn tab_cells() -> Vec<&'static str> {
    let key = TAB_ON[0]
        .split('\t')
        .position(|column| column == "gnomADg_AF")
        .expect("the key's column");
    TAB_ON[1..]
        .iter()
        .map(|row| row.split('\t').nth(key).expect("the key's cell"))
        .collect()
}

/// A row's `gnomADg_AF` holds the source's value and the frequency of the
/// row's own allele: `0.004` alone on the rows of `G`, for which the cache
/// holds no frequency, whichever allele the record writes first, and one
/// frequency on each row of the record whose two alleles each hold one. The
/// default, tab and VCF outputs print the rows whole.
#[test]
fn a_whole_record_custom_key_named_as_a_frequency_column_lists_the_frequency_of_the_rows_own_allele(
) {
    assert_eq!(
        tab_cells(),
        [
            "0.004,7.994e-06",
            "0.004",
            "0.004",
            "0.004,7.994e-06",
            "0.005,3.593e-05",
            "0.005,1.794e-05",
            "0.005,3.593e-05",
            "0.005,1.794e-05",
        ],
        "the cells the expected tab rows hold under the key"
    );
    let tmp = tempfile::tempdir().unwrap();

    let text = fs::read_to_string(annotate(tmp.path(), "default")).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(rows, DEFAULT_ON, "default format");

    let text = fs::read_to_string(annotate(tmp.path(), "tab")).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with("##")).collect();
    assert_eq!(rows, TAB_ON, "tab format");

    let text = fs::read_to_string(annotate(tmp.path(), "vcf")).unwrap();
    let format = text
        .lines()
        .find_map(|l| l.strip_prefix("##INFO=<ID=CSQ,"))
        .and_then(|l| l.split_once("Format: "))
        .map(|(_, format)| format.trim_end_matches("\">"))
        .expect("the CSQ header line");
    let mut lines = vec![format];
    lines.extend(text.lines().filter(|l| !l.starts_with('#')));
    assert_eq!(lines, VCF_ON, "VCF format");
}

/// The Parquet output stores the key as a string that holds the `--tab` cell
/// of each row.
#[cfg(feature = "parquet")]
#[test]
fn the_parquet_cell_of_such_a_key_holds_the_frequency_of_the_rows_own_allele() {
    let columns: Vec<&str> = TAB_ON[0].trim_start_matches('#').split('\t').collect();
    let places = ["Uploaded_variation", "Allele", "Feature", "gnomADg_AF"].map(|name| {
        columns
            .iter()
            .position(|column| *column == name)
            .expect("the column's place")
    });
    let row_and_cell = |row: &str| -> Vec<String> {
        let cells: Vec<&str> = row.split('\t').collect();
        places
            .iter()
            .map(|&place| cells[place].to_string())
            .collect()
    };
    let tmp = tempfile::tempdir().unwrap();
    let out = annotate(tmp.path(), "parquet");
    let mut rows: Vec<Vec<String>> = vep_io::output::parquet::tab_rows(&out)
        .unwrap()
        .iter()
        .map(|row| row_and_cell(row))
        .collect();
    rows.sort_unstable();
    let mut expected: Vec<Vec<String>> = TAB_ON[1..].iter().map(|row| row_and_cell(row)).collect();
    expected.sort_unstable();
    assert_eq!(rows, expected);
}
