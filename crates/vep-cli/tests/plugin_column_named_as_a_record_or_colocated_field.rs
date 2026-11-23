// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A column a plugin declares under the name of a co-located variant field,
//! of a field of the record's own hash or of a `--custom` source's key, on
//! the JSON entries the binary writes for two records of committed corpora,
//! each read over the cache of the `GRCh38-everything` corpus, which holds
//! both records' co-located variants: line 9 of the custom corpus's input,
//! which overlaps no transcript, and line 151 of the `GRCh38-everything`
//! input, which overlaps five. TranscriptAnnotator reads, under `prefix=0`, a
//! table the test writes whose one row names one transcript of the second
//! record: its columns are named `gnomADg_AF`, `CLIN_SIG` and `Comment` in
//! one table, `VARIANT_CLASS` and `Comment` in another, and `SOURCE`,
//! `Comment_FILTER` and `Comment` in a third, which is read beside the custom
//! corpus's site file as a source named `SOURCE` or `Comment`. The expected
//! records are the ones Ensembl VEP 116.3 writes for the same records, tables
//! and flags, their `input` key aside.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// A record of a committed corpus: the corpus, the 1-based line of its
/// `variants.vcf` and the text that line opens with.
struct Record {
    corpus: &'static str,
    line: usize,
    opens_with: &'static str,
}

/// The record that overlaps no transcript.
const INTERGENIC: Record = Record {
    corpus: "GRCh38-custom",
    line: 9,
    opens_with: "21\t5380280\t",
};

/// The record that overlaps five transcripts, one of which the table names.
const ON_TRANSCRIPTS: Record = Record {
    corpus: "GRCh38-everything",
    line: 151,
    opens_with: "21\t44499814\t",
};

/// The transcript of [`ON_TRANSCRIPTS`] the table's row names.
const NAMED_TRANSCRIPT: &str = "ENST00000323084";

/// A table TranscriptAnnotator reads: its file name, the names of the columns
/// the plugin writes, which follow the five it matches a row by, and the
/// cells of its one row under them.
struct Table {
    name: &'static str,
    columns: &'static str,
    cells: &'static str,
}

/// A table whose first two columns are named as a frequency column and as a
/// co-located variant field.
const COLOCATED_COLUMNS: Table = Table {
    name: "colocated_columns.tsv.gz",
    columns: "gnomADg_AF\tCLIN_SIG\tComment",
    cells: "0.61\tfrom_the_table\tNamed",
};

/// A table whose first column is named as a field of the record's own hash.
const RECORD_FIELD_COLUMN: Table = Table {
    name: "record_field_column.tsv.gz",
    columns: "VARIANT_CLASS\tComment",
    cells: "from_the_table\tNamed",
};

/// A table whose columns are named as keys of a `--custom` source named
/// `SOURCE` or `Comment`: the source's own key and the `FILTER` key of its
/// records.
const CUSTOM_KEY_COLUMNS: Table = Table {
    name: "custom_key_columns.tsv.gz",
    columns: "SOURCE\tComment_FILTER\tComment",
    cells: "from_the_table\tfrom_the_table\tNamed",
};

/// Ensembl VEP's record of [`INTERGENIC`] under [`COLOCATED_COLUMNS`] without
/// `--af_gnomadg`.
const INTERGENIC_FLAG_OFF: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// Ensembl VEP's record of [`INTERGENIC`] under [`COLOCATED_COLUMNS`] and
/// `--af_gnomadg`.
const INTERGENIC_FLAG_ON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","colocated_variants":[{"allele_string":"A/C/G","end":5380280,"frequencies":{"C":{"gnomadg":7.994e-06,"gnomadg_afr":0,"gnomadg_ami":0,"gnomadg_amr":0,"gnomadg_asj":0,"gnomadg_eas":0,"gnomadg_fin":0,"gnomadg_mid":0,"gnomadg_nfe":1.993e-05,"gnomadg_remaining":0,"gnomadg_sas":0}},"id":"rs1444031311","seq_region_name":"21","start":5380280,"strand":1}],"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// Ensembl VEP's record of [`ON_TRANSCRIPTS`] under [`COLOCATED_COLUMNS`]
/// without `--af_gnomadg`.
const ON_TRANSCRIPTS_FLAG_OFF: &str = r#"{"allele_string":"A/C","assembly_name":"GRCh38","end":44499814,"id":"rs140818230","most_severe_consequence":"missense_variant","seq_region_name":"21","start":44499814,"strand":1,"transcript_consequences":[{"amino_acids":"L/R","cdna_end":2037,"cdna_start":2037,"cds_end":1979,"cds_start":1979,"clin_sig":"from_the_table","codons":"cTc/cGc","comment":"Named","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","gnomadg_af":0.61,"impact":"MODERATE","protein_end":660,"protein_start":660,"strand":-1,"transcript_id":"ENST00000323084","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000622854","variant_allele":"C"},{"cdna_end":2074,"cdna_start":2074,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000175894","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000642437","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000698336","variant_allele":"C"},{"amino_acids":"L/R","cdna_end":2877,"cdna_start":2877,"cds_end":2105,"cds_start":2105,"codons":"cTc/cGc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","impact":"MODERATE","protein_end":702,"protein_start":702,"strand":-1,"transcript_id":"ENST00000943283","variant_allele":"C"}]}"#;

/// Ensembl VEP's record of [`ON_TRANSCRIPTS`] under [`COLOCATED_COLUMNS`] and
/// `--af_gnomadg`.
const ON_TRANSCRIPTS_FLAG_ON: &str = r#"{"allele_string":"A/C","assembly_name":"GRCh38","colocated_variants":[{"allele_string":"A/C","clin_sig":["uncertain_significance"],"clin_sig_allele":"C:conflicting_classifications_of_pathogenicity;C:uncertain_significance","clin_sig_ref_allele":"A","end":44499814,"frequencies":{"C":{"af":0.0004,"afr":0.0015,"amr":0,"eas":0,"eur":0,"gnomade":3.575e-05,"gnomade_afr":0.001138,"gnomade_amr":4.969e-05,"gnomade_asj":0,"gnomade_eas":0,"gnomade_fin":0,"gnomade_mid":0,"gnomade_nfe":1.825e-06,"gnomade_remaining":0.0001696,"gnomade_sas":0,"gnomadg":0.0004729,"gnomadg_afr":0.001613,"gnomadg_ami":0,"gnomadg_amr":0.000196,"gnomadg_asj":0,"gnomadg_eas":0,"gnomadg_fin":0,"gnomadg_mid":0,"gnomadg_nfe":1.471e-05,"gnomadg_remaining":0.000473,"gnomadg_sas":0,"sas":0}},"id":"rs140818230","phenotype_or_disease":1,"seq_region_name":"21","start":44499814,"strand":1,"var_synonyms":{"ClinVar":["RCV003339673","RCV001545588","VCV001186484"]}}],"end":44499814,"id":"rs140818230","most_severe_consequence":"missense_variant","seq_region_name":"21","start":44499814,"strand":1,"transcript_consequences":[{"amino_acids":"L/R","cdna_end":2037,"cdna_start":2037,"cds_end":1979,"cds_start":1979,"clin_sig":"from_the_table","codons":"cTc/cGc","comment":"Named","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","gnomadg_af":0.61,"impact":"MODERATE","protein_end":660,"protein_start":660,"strand":-1,"transcript_id":"ENST00000323084","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000622854","variant_allele":"C"},{"cdna_end":2074,"cdna_start":2074,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000175894","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000642437","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000698336","variant_allele":"C"},{"amino_acids":"L/R","cdna_end":2877,"cdna_start":2877,"cds_end":2105,"cds_start":2105,"codons":"cTc/cGc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","impact":"MODERATE","protein_end":702,"protein_start":702,"strand":-1,"transcript_id":"ENST00000943283","variant_allele":"C"}]}"#;

/// Ensembl VEP's record of [`ON_TRANSCRIPTS`] under [`RECORD_FIELD_COLUMN`]
/// without `--variant_class`.
const RECORD_FIELD_FLAG_OFF: &str = r#"{"allele_string":"A/C","assembly_name":"GRCh38","end":44499814,"id":"rs140818230","most_severe_consequence":"missense_variant","seq_region_name":"21","start":44499814,"strand":1,"transcript_consequences":[{"amino_acids":"L/R","cdna_end":2037,"cdna_start":2037,"cds_end":1979,"cds_start":1979,"codons":"cTc/cGc","comment":"Named","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","impact":"MODERATE","protein_end":660,"protein_start":660,"strand":-1,"transcript_id":"ENST00000323084","variant_allele":"C","variant_class":"from_the_table"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000622854","variant_allele":"C"},{"cdna_end":2074,"cdna_start":2074,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000175894","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000642437","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000698336","variant_allele":"C"},{"amino_acids":"L/R","cdna_end":2877,"cdna_start":2877,"cds_end":2105,"cds_start":2105,"codons":"cTc/cGc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","impact":"MODERATE","protein_end":702,"protein_start":702,"strand":-1,"transcript_id":"ENST00000943283","variant_allele":"C"}]}"#;

/// Ensembl VEP's record of [`ON_TRANSCRIPTS`] under [`RECORD_FIELD_COLUMN`]
/// and `--variant_class`.
const RECORD_FIELD_FLAG_ON: &str = r#"{"allele_string":"A/C","assembly_name":"GRCh38","end":44499814,"id":"rs140818230","most_severe_consequence":"missense_variant","seq_region_name":"21","start":44499814,"strand":1,"transcript_consequences":[{"amino_acids":"L/R","cdna_end":2037,"cdna_start":2037,"cds_end":1979,"cds_start":1979,"codons":"cTc/cGc","comment":"Named","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","impact":"MODERATE","protein_end":660,"protein_start":660,"strand":-1,"transcript_id":"ENST00000323084","variant_allele":"C","variant_class":"from_the_table"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000622854","variant_allele":"C"},{"cdna_end":2074,"cdna_start":2074,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000175894","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000642437","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4295,"gene_id":"ENSG00000228709","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000698336","variant_allele":"C"},{"amino_acids":"L/R","cdna_end":2877,"cdna_start":2877,"cds_end":2105,"cds_start":2105,"codons":"cTc/cGc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000175894","impact":"MODERATE","protein_end":702,"protein_start":702,"strand":-1,"transcript_id":"ENST00000943283","variant_allele":"C"}],"variant_class":"SNV"}"#;

/// Ensembl VEP's record of [`INTERGENIC`] under [`CUSTOM_KEY_COLUMNS`] beside
/// a `--custom` source named `SOURCE`.
const INTERGENIC_BESIDE_A_SOURCE_NAMED_SOURCE: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"SOURCE":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// Ensembl VEP's record of [`INTERGENIC`] under [`CUSTOM_KEY_COLUMNS`] beside
/// a `--custom` source named `Comment`.
const INTERGENIC_BESIDE_A_SOURCE_NAMED_COMMENT: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"Comment":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The line of `record` in its corpus's input.
fn record_line(record: &Record) -> String {
    let input = fs::read_to_string(corpus(record.corpus).join("variants.vcf")).unwrap();
    let line = input
        .lines()
        .nth(record.line - 1)
        .expect("the record's line");
    assert!(
        line.starts_with(record.opens_with),
        "line {} of the {} input is the record",
        record.line,
        record.corpus
    );
    line.to_string()
}

/// `record` alone, written under `dir`: the header lines of its corpus's
/// input and its own line.
fn record_file(dir: &Path, record: &Record) -> PathBuf {
    let input = fs::read_to_string(corpus(record.corpus).join("variants.vcf")).unwrap();
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    text.push_str(&record_line(record));
    text.push('\n');
    let path = dir.join(format!("{}.vcf", record.corpus));
    fs::write(&path, text).unwrap();
    path
}

/// `table` written under `dir`, bgzip-compressed and indexed on its second
/// column: a header line naming the five columns TranscriptAnnotator matches
/// a row by and then the table's own, and one row at the site and alleles of
/// [`ON_TRANSCRIPTS`] that names [`NAMED_TRANSCRIPT`] and holds the table's
/// cells.
fn table_file(dir: &Path, table: &Table) -> PathBuf {
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;

    let line = record_line(&ON_TRANSCRIPTS);
    let cells: Vec<&str> = line.split('\t').collect();
    let header = format!("#Chrom\tPos\tRef\tAlt\tTranscript\t{}", table.columns);
    let row = format!(
        "{}\t{}\t{}\t{}\t{NAMED_TRANSCRIPT}\t{}",
        cells[0], cells[1], cells[3], cells[4], table.cells
    );
    common::indexed_lines(
        dir,
        table.name,
        HeaderBuilder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .build(),
        &[header.as_str()],
        &[row.as_str()],
        |fields| {
            let position = fields[1].parse().unwrap();
            (position, position)
        },
    )
}

/// The output of a run over `record` alone in `format` (`json` or `tab`)
/// under `flags`, with TranscriptAnnotator reading `table` under `prefix=0`.
fn annotate(dir: &Path, record: &Record, table: &Table, format: &str, flags: &[&str]) -> String {
    let out = dir.join(format!(
        "{}.{}.{format}.{}.out",
        record.corpus,
        table.name,
        flags.len()
    ));
    let plugin = format!(
        "TranscriptAnnotator,file={},prefix=0",
        table_file(dir, table).display()
    );
    let run = common::isolated_vep(dir)
        .arg("-i")
        .arg(record_file(dir, record))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-everything").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(["--plugin", &plugin])
        .arg(format!("--{format}"))
        .args(flags)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "vep {format} {flags:?}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    fs::read_to_string(out).unwrap()
}

/// The JSON record of a run over `record` under `table` and `flags`, its
/// `input` key and the `vep_rs_version` key aside.
fn json_record(dir: &Path, record: &Record, table: &Table, flags: &[&str]) -> serde_json::Value {
    let text = annotate(dir, record, table, "json", flags);
    let mut value: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = value.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    value
}

/// The JSON factory writes a co-located variant's frequencies under
/// `colocated_variants` and sets none on a row's hash
/// (OutputFactory/JSON.pm:278-286, 513-515), and a plugin's keys reach a row
/// only where the plugin ran on it (OutputFactory.pm:2228-2241).
/// TranscriptAnnotator runs on no intergenic row, so neither entry of the
/// record carries a `gnomadg_af`, with `--af_gnomadg` as without it.
#[test]
fn a_frequency_column_a_plugin_declares_stays_off_an_entry_the_plugin_left_unwritten() {
    let tmp = tempfile::tempdir().unwrap();
    for (flags, expected) in [
        (&[][..], INTERGENIC_FLAG_OFF),
        (&["--af_gnomadg"][..], INTERGENIC_FLAG_ON),
    ] {
        let expected: serde_json::Value = serde_json::from_str(expected).unwrap();
        assert_eq!(
            json_record(tmp.path(), &INTERGENIC, &COLOCATED_COLUMNS, flags),
            expected,
            "{flags:?}"
        );
    }
}

/// On a record with transcript rows the entry of the transcript the table
/// names carries the table's cells under `gnomadg_af` and `clin_sig`, and the
/// other four entries carry neither key: under `--af_gnomadg` the co-located
/// variant's frequency and clinical significance stand under
/// `colocated_variants` alone (OutputFactory/JSON.pm:493-515).
#[test]
fn a_colocated_field_a_plugin_declares_is_on_the_entry_the_plugin_wrote_alone() {
    let tmp = tempfile::tempdir().unwrap();
    for (flags, expected) in [
        (&[][..], ON_TRANSCRIPTS_FLAG_OFF),
        (&["--af_gnomadg"][..], ON_TRANSCRIPTS_FLAG_ON),
    ] {
        let expected: serde_json::Value = serde_json::from_str(expected).unwrap();
        assert_eq!(
            json_record(tmp.path(), &ON_TRANSCRIPTS, &COLOCATED_COLUMNS, flags),
            expected,
            "{flags:?}"
        );
    }
}

/// The JSON factory builds each row's hash from an empty one, the record's
/// keys staying on the record (OutputFactory/JSON.pm:261-262, 320), so under
/// `--variant_class` the record carries its class as `variant_class`, the
/// entry of the transcript the table names carries the table's cell under
/// that key, and the other four entries carry no such key.
#[test]
fn a_record_field_a_plugin_declares_is_on_the_entry_the_plugin_wrote_alone() {
    let tmp = tempfile::tempdir().unwrap();
    for (flags, expected) in [
        (&[][..], RECORD_FIELD_FLAG_OFF),
        (&["--variant_class"][..], RECORD_FIELD_FLAG_ON),
    ] {
        let expected: serde_json::Value = serde_json::from_str(expected).unwrap();
        assert_eq!(
            json_record(tmp.path(), &ON_TRANSCRIPTS, &RECORD_FIELD_COLUMN, flags),
            expected,
            "{flags:?}"
        );
    }
}

/// The JSON factory moves the custom annotations to the record before any row
/// hash is built (OutputFactory/JSON.pm:249-258), so the row builder finds no
/// record to push under a source's key (OutputFactory.pm:1426-1435) and an
/// entry carries a key a plugin declares only as the plugin's own. Beside the
/// custom corpus's site file read as a source named `SOURCE` or `Comment`,
/// whose record names the first alternate, neither entry of the record
/// carries a `source`, a `comment` or a `comment_filter`.
#[test]
fn a_custom_key_a_plugin_declares_stays_off_an_entry_the_plugin_left_unwritten() {
    let tmp = tempfile::tempdir().unwrap();
    for (short_name, expected) in [
        ("SOURCE", INTERGENIC_BESIDE_A_SOURCE_NAMED_SOURCE),
        ("Comment", INTERGENIC_BESIDE_A_SOURCE_NAMED_COMMENT),
    ] {
        let source = format!(
            "file={},short_name={short_name},format=vcf,type=exact,fields=AF%CLASS",
            corpus("GRCh38-custom")
                .join("custom/sites.vcf.gz")
                .display()
        );
        let expected: serde_json::Value = serde_json::from_str(expected).unwrap();
        assert_eq!(
            json_record(
                tmp.path(),
                &INTERGENIC,
                &CUSTOM_KEY_COLUMNS,
                &["--custom", &source]
            ),
            expected,
            "{short_name}"
        );
    }
}

/// The tab format reads the same columns from the row's hash, which the base
/// factory fills with the co-located variant's values before the plugins run
/// (OutputFactory.pm:1438-1443, 398): under `--af_gnomadg` a row the plugin
/// left unwritten prints the frequency and the clinical significance in the
/// plugin's columns, as Ensembl VEP's rows of the record do, and the named
/// transcript's row prints the table's cells.
#[test]
fn the_tab_format_prints_the_colocated_value_in_the_plugins_column() {
    let tmp = tempfile::tempdir().unwrap();
    let text = annotate(
        tmp.path(),
        &ON_TRANSCRIPTS,
        &COLOCATED_COLUMNS,
        "tab",
        &["--af_gnomadg"],
    );
    let mut lines = text.lines().filter(|line| !line.starts_with("##"));
    let columns: Vec<&str> = lines
        .next()
        .expect("the column line")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let place = |name: &str| {
        let places: Vec<usize> = (0..columns.len())
            .filter(|&at| columns[at] == name)
            .collect();
        assert_eq!(places.len(), 1, "{name} is listed once: {columns:?}");
        places[0]
    };
    let (feature, frequency, significance) =
        (place("Feature"), place("gnomADg_AF"), place("CLIN_SIG"));
    let cells: Vec<(String, String, String)> = lines
        .map(|row| {
            let row: Vec<&str> = row.split('\t').collect();
            (
                row[feature].to_string(),
                row[frequency].to_string(),
                row[significance].to_string(),
            )
        })
        .collect();
    let unwritten = "conflicting_classifications_of_pathogenicity,uncertain_significance";
    let expected: Vec<(String, String, String)> = [
        (NAMED_TRANSCRIPT, "0.61", "from_the_table"),
        ("ENST00000622854", "0.0004729", unwritten),
        ("ENST00000642437", "0.0004729", unwritten),
        ("ENST00000698336", "0.0004729", unwritten),
        ("ENST00000943283", "0.0004729", unwritten),
    ]
    .iter()
    .map(|(a, b, c)| (a.to_string(), b.to_string(), c.to_string()))
    .collect();
    assert_eq!(cells, expected);
}
