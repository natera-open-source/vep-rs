// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The field classes the dispatch decides in its own branches, on planted rows: through the
//! comparison of two outputs and through the dispatch asked directly.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use vep_formats_compare::classes;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::consts::UTRANNOTATOR_ANNOTATION_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::lists::reordered_items;
use vep_formats_compare::lists::reordered_pair_groups;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::record_alleles::vcf_first_alternate;
use vep_formats_compare::rules::sequence_edges::turned_allele;
use vep_formats_compare::rules::structural::_SO_TERMS;
use vep_formats_compare::rules::utrannotator::UtrannotatorSeedSweep;
use vep_formats_compare::rules::utrannotator::STRUCTURAL_CLASS_WORDS;
use vep_formats_compare::rules::utrannotator::UTRANNOTATOR_CONSEQUENCE_FIELDS;
use vep_formats_compare::rules::utrannotator::UTRANNOTATOR_TERMS;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

const UTR_GAIN: &str = "5_prime_UTR_premature_start_codon_gain_variant";
const UTR_LOSS: &str = "5_prime_UTR_premature_start_codon_loss_variant";
const UTR_STOP_LOSS: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
const UTR_FRAMESHIFT: &str = "5_prime_UTR_uORF_frameshift_variant";
const UTR_STRUCTURAL: &str = "utrannotator_structural_allele_text";
const UTR_PRECEDENCE: &str = "utrannotator_consequence_precedence";
const ITEM_ORDER: &str = "nearest_exon_jb_item_order";
const PAIR_ORDER: &str = "utrannotator_annotation_pair_order";

/// The two exon items of a deletion over three exons of one transcript, which the reference
/// writes in either order over repeated runs of one command.
const NEAREST_EXON_JB_ITEMS: (&str, &str) =
    ("ENSE00000000031|3|start|117", "ENSE00000000041|3|end|207");
const NEAREST_EXON_JB_TERMS: &str = "splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant";
/// A third field, neither order class's own, whose items a test plants in another order.
const DOMAINS_ITEMS: [&str; 2] = ["Pfam:PF1", "PANTHER:PTHR1"];

/// What was applied, class by class.
fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

fn once(klass: &str) -> BTreeMap<String, i64> {
    count_map(&[(klass, 1)])
}

fn nothing() -> BTreeMap<String, i64> {
    BTreeMap::new()
}

/// Whether the column `field` of the report is set aside once under `klass` and under
/// nothing else.
fn masked_once(report: &Value, field: &str, klass: &str) -> bool {
    at(report, &["columns", field])
        .get("masked")
        .is_some_and(|masked| *masked == json(&format!("{{\"{klass}\": 1}}")))
}

fn differing(report: &Value, field: &str) -> i64 {
    number(report, &["columns", field, "differing"])
}

/// Whether the comparison under `klass` is red with nothing set aside and `field` counted
/// once as differing.
fn refused(reference: &str, candidate: &str, klass: &str, field: &str) -> bool {
    let (report, documented) = compare_default(reference, candidate, &[klass]);
    verdict(&report) == "red" && applied(&documented).is_empty() && differing(&report, field) == 1
}

/// The first entry of a default-format output.
fn first_entry(text: &str) -> Entry {
    parse_default(text)
        .expect("the text parses")
        .into_iter()
        .next()
        .expect("an entry")
}

/// What the dispatch reads beside the two values in a test that asks it directly.
#[derive(Default)]
struct Beside<'a> {
    a: Option<&'a Entry>,
    record_rows: Option<&'a [&'a Entry]>,
    assembly: Option<&'a str>,
    columns: &'a [String],
}

fn accepts_beside(
    klass: &str,
    e: &Entry,
    name: &str,
    expected: &str,
    actual: &str,
    beside: &Beside,
) -> bool {
    let side = Side::default();
    let context = RuleContext {
        record_rows: beside.record_rows,
        candidate_rows: None,
        declared: &[],
        assembly: beside.assembly,
        dbnsfp_columns: beside.columns,
        side: &side,
    };
    field_class_accepts(klass, e, name, expected, actual, beside.a, &context)
        .expect("the class answers")
}

/// The dispatch asked of one field with the candidate's entry and nothing else beside it.
fn accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected: &str,
    actual: &str,
    a: Option<&Entry>,
) -> bool {
    accepts_beside(
        klass,
        e,
        name,
        expected,
        actual,
        &Beside {
            a,
            ..Beside::default()
        },
    )
}

/// The header line that declares the CSQ fields `format` of a VCF output.
fn csq_header(format: &str) -> String {
    format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {format}\">\n"
    )
}

/// A VCF output of one line: the CSQ header of `format` and the input line `input` with
/// its INFO column replaced by the CSQ value `csq`.
fn vcf_text(format: &str, input: &str, csq: &str) -> String {
    format!(
        "{}{}CSQ={csq}\n",
        csq_header(format),
        &input[..input.len() - 1]
    )
}

fn vcf_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_vcf(text, inputs).expect("the VCF text parses")
}

fn json_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_json(text, inputs).expect("the JSON text parses")
}

/// The verdict of two entry lists under one class, and what was applied.
fn under_class(
    klass: &str,
    expected: &[Entry],
    actual: &[Entry],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(&[klass]);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), applied(&documented))
}

fn green_once(klass: &str) -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), once(klass))
}

/// A 5' UTR row of the default format with its cDNA_position and the two UTRAnnotator fields.
fn utr_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    cdna: &str,
    consequence: &str,
    annotation: &str,
) -> String {
    let extra = format!(
        "IMPACT=MODIFIER;STRAND=1;5UTR_consequence={consequence};5UTR_annotation={annotation}"
    );
    default_row(
        name,
        location,
        allele,
        feature,
        "5_prime_UTR_variant",
        &extra,
    )
    .replacen(
        "\t-\t-\t-\t-\t-\t-\t",
        &format!("\t{cdna}\t-\t-\t-\t-\t-\t"),
        1,
    )
}

/// The seven pairs of a substitution in the 5' UTR of one transcript that destroys an
/// upstream start codon.
const UTR_ANNOTATION_PAIRS: [(&str, &str); 7] = [
    ("CapDistanceToStart", "407"),
    ("DistanceToCDS", "23"),
    ("DistanceToStop", "24"),
    ("Evidence", "False"),
    ("KozakContext", "ATAATGG"),
    ("KozakStrength", "Strong"),
    ("type", "OutOfFrame_oORF"),
];
/// One order the reference writes the pairs in, and the byte order of their keys.
const UTR_ANNOTATION_HASH_ORDER: [&str; 7] = [
    "DistanceToStop",
    "CapDistanceToStart",
    "KozakContext",
    "KozakStrength",
    "DistanceToCDS",
    "Evidence",
    "type",
];
const UTR_ANNOTATION_KEY_ORDER: [&str; 7] = [
    "CapDistanceToStart",
    "DistanceToCDS",
    "DistanceToStop",
    "Evidence",
    "KozakContext",
    "KozakStrength",
    "type",
];

/// One `5UTR_annotation` group: the pairs of `keys` in that order, `pairs` replacing or
/// adding values.
fn utr_annotation(keys: &[&str], pairs: &[(&str, &str)]) -> String {
    let value = |key: &str| {
        pairs
            .iter()
            .chain(&UTR_ANNOTATION_PAIRS)
            .find(|(held, _)| *held == key)
            .map(|(_, value)| *value)
            .expect("a pair of the key")
    };
    keys.iter()
        .map(|key| format!("{key}={}", value(key)))
        .collect::<Vec<_>>()
        .join(":")
}

fn hash_order() -> String {
    utr_annotation(&UTR_ANNOTATION_HASH_ORDER, &[])
}

fn key_order() -> String {
    utr_annotation(&UTR_ANNOTATION_KEY_ORDER, &[])
}

/// The substitution's default-format row on ENST00000000032 with a `5UTR_annotation` under
/// `consequence` and the DOMAINS items `domains`, or none.
fn utr_annotation_row(annotation: &str, consequence: &str, domains: &[&str]) -> String {
    let row = utr_row(
        "utr_substitution",
        "21:2000984",
        "G",
        "ENST00000000032",
        "410",
        consequence,
        annotation,
    );
    if domains.is_empty() {
        default_text(&[row])
    } else {
        default_text(&[row.replacen(
            "STRAND=1;",
            &format!("STRAND=1;DOMAINS={};", domains.join(",")),
            1,
        )])
    }
}

fn utr_annotation_text(annotation: &str) -> String {
    utr_annotation_row(annotation, UTR_LOSS, &[])
}

/// The deletion's default-format row on ENST00000000031 with the NearestExonJB `items`, or
/// none, and the DOMAINS items `domains`, or none.
fn nearest_exon_jb_row(items: &[&str], domains: &[&str]) -> String {
    let mut extra = "IMPACT=HIGH;STRAND=1".to_string();
    if !domains.is_empty() {
        extra += &format!(";DOMAINS={}", domains.join(","));
    }
    if !items.is_empty() {
        extra += &format!(";NearestExonJB={}", items.join(","));
    }
    default_text(&[default_row(
        "exon_deletion",
        "21:1000011-1001793",
        "-",
        "ENST00000000031",
        NEAREST_EXON_JB_TERMS,
        &extra,
    )])
}

fn nearest_exon_jb_text(items: &[&str]) -> String {
    nearest_exon_jb_row(items, &[])
}

/// The six cells a planted row is made of: its name, Location, Allele, Feature, Consequence
/// and Extra text.
type Row<'a> = [&'a str; 6];

/// One default-format row whose Extra text ends with the NearestExonJB `value`, left out
/// when `value` is empty, and whose cDNA_position is `cdna`.
fn nearest_exon_row_at(row: Row, value: &str, cdna: &str) -> String {
    let extra = if value.is_empty() {
        row[5].to_string()
    } else {
        format!("{};NearestExonJB={value}", row[5])
    };
    let line = default_row(row[0], row[1], row[2], row[3], row[4], &extra);
    let mut columns: Vec<&str> = line.split('\t').collect();
    columns[7] = cdna;
    default_text(&[columns.join("\t")])
}

fn nearest_exon_row(row: Row, value: &str) -> String {
    nearest_exon_row_at(row, value, "-")
}

/// `row` with its Extra text replaced.
fn with_extra<'a>(row: Row<'a>, extra: &'a str) -> Row<'a> {
    [row[0], row[1], row[2], row[3], row[4], extra]
}

const SHIFT_PROBE: Row = [
    "shift_probe",
    "21:37480826",
    "-",
    "ENST00000338785",
    "splice_donor_variant",
    "IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1",
];
const IDLESS_TWO_EXONS: Row = [
    "idless_two_exons",
    "21:10578401-10578530",
    "-",
    "compmerge.2.Heart.chr21",
    "splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,non_coding_transcript_exon_variant,intron_variant",
    "IMPACT=HIGH;STRAND=1;SOURCE=RefSeq",
];
const IDLESS_INTRONIC: Row = [
    "idless_intronic",
    "21:10521800",
    "A",
    "compmerge.2.Heart.chr21",
    "intron_variant,non_coding_transcript_variant",
    "IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq",
];

/// The default-format row of the deletion 21:10578401-10578530 on compmerge.2.Heart.chr21,
/// which begins in one exon and ends in the next.
fn two_exon_row_at(value: &str, cdna: &str) -> String {
    nearest_exon_row_at(IDLESS_TWO_EXONS, value, cdna)
}

fn two_exon_row(value: &str) -> String {
    two_exon_row_at(value, "1266-1297")
}

/// The default-format row of 21 30166352 G>GAG on ENST00000286808 under `--hgvs`.
#[derive(Clone, Copy)]
struct DownstreamRow<'a> {
    strand: &'a str,
    offset: &'a str,
    notation: &'a str,
    consequence: &'a str,
    length_change: &'a str,
}

impl Default for DownstreamRow<'_> {
    fn default() -> Self {
        DownstreamRow {
            strand: "-1",
            offset: "-7",
            notation: "c.271_272dup",
            consequence: "frameshift_variant",
            length_change: "-132",
        }
    }
}

impl DownstreamRow<'_> {
    fn text(&self, protein: &str) -> String {
        let extra = format!(
            "IMPACT=HIGH;STRAND={};HGVSc=ENST00000286808:{};HGVS_OFFSET={};DownstreamProtein={protein};ProteinLengthChange={}",
            self.strand, self.notation, self.offset, self.length_change
        );
        default_text(&[default_row(
            "cldn17_ins_AG",
            "21:30166352-30166353",
            "AG",
            "ENST00000286808",
            self.consequence,
            &extra,
        )])
    }
}

fn downstream_row(protein: &str) -> String {
    DownstreamRow::default().text(protein)
}

/// The default-format row of the record `<contig>` 3460 rs199476118 G A on ENST00000361390,
/// its Extra text ending with `fields` when there are any.
fn mitochondrial_row_at(fields: &str, location: &str) -> String {
    let extra = if fields.is_empty() {
        "IMPACT=MODERATE;STRAND=1".to_string()
    } else {
        format!("IMPACT=MODERATE;STRAND=1;{fields}")
    };
    default_text(&[default_row(
        "rs199476118",
        location,
        "A",
        "ENST00000361390",
        "missense_variant",
        &extra,
    )])
}

fn mitochondrial_row(fields: &str) -> String {
    mitochondrial_row_at(fields, "MT:3460")
}

const MITOCHONDRIAL: &str = "plugin_contig_mitochondrial_spelling";
/// One field of each tabix-backed plugin the class reads but AlphaMissense, with a value of
/// the field's shape and of no file.
const MITOCHONDRIAL_PLUGIN_FIELDS: [&str; 12] = [
    "AVI_PHRED=12.5",
    "BayesDel=0.1",
    "CADD_PHRED=12.5",
    "ClinPred=0.5",
    "ada_score=0.9",
    "LOEUF=0.5",
    "OpenTargets_gwasLocusToGeneScore=0.7",
    "PrimateAI=0.6",
    "REVEL=0.5",
    "SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1",
    "SpliceVault_site_type=donor",
    "gnomAD_mean_cov=31.2",
];
const MITOCHONDRIAL_ALPHAMISSENSE_FIELDS: &str =
    "am_class=likely_pathogenic;am_pathogenicity=0.624";

/// The record's row without and with `fields` compared under the class beside `manifest`, a
/// made-up one that documents nothing.
fn compare_mitochondrial_rows(
    manifest: Option<&str>,
    fields: &str,
    location: &str,
) -> (Value, Documented) {
    let manifest = manifest.map(json);
    let mut documented = Documented::from_manifest(
        manifest.as_ref(),
        &masks(&[MITOCHONDRIAL]),
        None,
        None,
        None,
    )
    .expect("the manifest reads");
    let report = compare_texts(
        &mitochondrial_row_at("", location),
        &mitochondrial_row_at(fields, location),
        &mut documented,
    );
    (report, documented)
}

const DBNSFP_NAMED_M: &str = "dbnsfp_mitochondrial_record_named_m";

/// A made-up manifest naming `assembly` whose own run hands dbNSFP a made-up table and two
/// of its columns.
fn dbnsfp_manifest(assembly: &str) -> String {
    format!(
        "{{\"assembly\": \"{assembly}\", \"flags\": [\"--dir_plugins\", \"{{plugins}}\", \"--plugin\", \"dbNSFP,{{corpus}}/plugin_data/dbNSFP4.9a_grch37_synthetic.gz,SIFT_score,Polyphen2_HDIV_score\"]}}"
    )
}

/// The record's row with the fields `reference` and with the fields `candidate` compared
/// under the dbNSFP class beside `manifest`.
fn compare_dbnsfp_rows(
    manifest: Option<&str>,
    reference: &str,
    candidate: &str,
    location: &str,
) -> (Value, Documented) {
    let manifest = manifest.map(json);
    let mut documented = Documented::from_manifest(
        manifest.as_ref(),
        &masks(&[DBNSFP_NAMED_M]),
        None,
        None,
        None,
    )
    .expect("the manifest reads");
    let report = compare_texts(
        &mitochondrial_row_at(reference, location),
        &mitochondrial_row_at(candidate, location),
        &mut documented,
    );
    (report, documented)
}

/// A value of the field's shape, not a score of any file.
const SPLICEAI_VALUE: &str = "SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1";

/// Default-format rows of the record 21 26965499 GC G,GA, written C/-/A at 21:26965500, on
/// ENST00000284987: (allele, the SpliceAI pair of the Extra text or nothing) each.
fn chopped_rows(name: &str, rows: &[(&str, &str)]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, spliceai)| {
            let extra = if spliceai.is_empty() {
                "IMPACT=MODIFIER;STRAND=1".to_string()
            } else {
                format!("IMPACT=MODIFIER;STRAND=1;{spliceai}")
            };
            default_row(
                name,
                "21:26965500",
                allele,
                "ENST00000284987",
                "intron_variant",
                &extra,
            )
        })
        .collect();
    default_text(&rows)
}

const GNOMAD_VEP_CUT: &str = "G|synonymous_variant|LOW|MT-ND1|ENSG00000198888|Transcript|ENST00000361390|protein_coding|1/1||ENST00000361390.2:c.135A>G|ENSP00000354687.2:p.Leu45";

fn gnomad_vep_whole() -> String {
    format!("{GNOMAD_VEP_CUT}=|135|135|45|L|ctA/ctG|1||1|SNV||HGNC|HGNC:7455|YES")
}

fn gnomad_vep_row(value: &str) -> String {
    let extra = if value.is_empty() {
        "IMPACT=LOW;STRAND=1".to_string()
    } else {
        format!("IMPACT=LOW;STRAND=1;gnomAD_vep={value}")
    };
    default_text(&[default_row(
        "nd1_3441_A_G",
        "MT:3441",
        "G",
        "ENST00000361390",
        "synonymous_variant",
        &extra,
    )])
}

/// The default-format row of a frameshift at 21:36759739 on ENST00001017196: its LoF pairs
/// and a LoF_info carrying the GERP total, the bases walked and the distance from the last
/// exon.
#[derive(Clone, Copy)]
struct GerpRow<'a> {
    walked: &'a str,
    reach: &'a str,
    consequence: &'a str,
}

impl Default for GerpRow<'_> {
    fn default() -> Self {
        GerpRow {
            walked: "118",
            reach: "12",
            consequence: "frameshift_variant",
        }
    }
}

impl GerpRow<'_> {
    fn text(&self, lof: &str, gerp: &str) -> String {
        let info = format!(
            "PERCENTILE:0.95,GERP_DIST:{gerp},BP_DIST:{},DIST_FROM_LAST_EXON:{},50_BP_RULE:FAIL",
            self.walked, self.reach
        );
        default_text(&[default_row(
            "2780984",
            "21:36759739",
            "-",
            "ENST00001017196",
            self.consequence,
            &format!("IMPACT=HIGH;STRAND=1;{lof};LoF_info={info}"),
        )])
    }
}

fn gerp_row(lof: &str, gerp: &str) -> String {
    GerpRow::default().text(lof, gerp)
}

const GERP: &str = "loftee_gerp_zoom_level_summary";
const GERP_REFERENCE: &str = "-69.0908176183701";
const GERP_CANDIDATE: &str = "-51.9470007419586";
const GERP_FILTERED: &str = "LoF=LC;LoF_filter=END_TRUNC";
const GERP_CLEAR: &str = "LoF=HC";

/// A structural record's 5' UTR row under `allele` with the two UTRAnnotator fields and
/// `more` after them in its Extra.
fn structural_utr_row(allele: &str, consequence: &str, annotation: &str, more: &str) -> String {
    utr_row(
        "sv",
        "21:10524589-10524628",
        allele,
        "ENST1",
        "?-208",
        consequence,
        annotation,
    ) + more
}

/// Under each class word a structural row carries as its Allele, and with each of the
/// plugin's five terms, the two fields the reference alone fills are set aside together, once
/// a row, whether the seeded runs are given or not; a row the engines agree on is untouched.
#[test]
fn the_structural_allele_class_accepts_both_fields_the_reference_alone_fills() {
    let spliced = "alt_type=uORF:alt_type_length=39:ref_type=uORF:ref_type_length=51";
    let existing = ";Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6";
    let agreeing = utr_row(
        "rs2",
        "21:400",
        "C",
        "ENST2",
        "40",
        UTR_GAIN,
        "CapDistanceToStart=38:DistanceToCDS=110:KozakStrength=Weak:type=uORF",
    );
    let class_words = [
        "Alu_deletion",
        "Alu_insertion",
        "HERV_deletion",
        "HERV_insertion",
        "LINE1_deletion",
        "LINE1_insertion",
        "SVA_deletion",
        "SVA_insertion",
        "copy_number_variation",
        "deletion",
        "duplication",
        "insertion",
        "inversion",
        "mobile_element_deletion",
        "mobile_element_insertion",
        "tandem_duplication",
        "tandem_repeat",
    ];
    for word in class_words {
        let reference = default_text(&[
            structural_utr_row(word, UTR_FRAMESHIFT, spliced, existing),
            agreeing.clone(),
        ]);
        let candidate = default_text(&[
            structural_utr_row(word, "-", "-", existing),
            agreeing.clone(),
        ]);
        let expected = parse_default(&reference).expect("the reference parses");
        let actual = parse_default(&candidate).expect("the candidate parses");
        for outputs in [0, 1, 2] {
            let mut side = Side::default();
            if outputs > 0 {
                let mut sweep = UtrannotatorSeedSweep::new();
                sweep.add_output(&expected);
                if outputs > 1 {
                    sweep.add_output(&actual);
                }
                side.utrannotator_sweep = Some(sweep);
            }
            let options = CompareOptions::new(&side);
            let mut documented = documented_with(&[UTR_STRUCTURAL]);
            let report = compare_entries(&expected, &actual, &mut documented, &options)
                .expect("the entries compare");
            assert_eq!(verdict(&report), "green", "{word}");
            assert_eq!(applied(&documented), once(UTR_STRUCTURAL), "{word}");
            assert!(
                masked_once(&report, "5UTR_consequence", UTR_STRUCTURAL),
                "{word}"
            );
            assert!(
                masked_once(&report, "5UTR_annotation", UTR_STRUCTURAL),
                "{word}"
            );
            let mut documented = documented_with(&[]);
            let report = compare_entries(&expected, &actual, &mut documented, &options)
                .expect("the entries compare");
            assert_eq!(verdict(&report), "red", "{word}");
        }
    }
    let unfilled = default_text(&[structural_utr_row("duplication", "-", "-", "")]);
    let mut terms: Vec<&str> = UTRANNOTATOR_TERMS.to_vec();
    terms.sort_unstable();
    let mut reference = String::new();
    for term in terms {
        reference = default_text(&[structural_utr_row(
            "duplication",
            term,
            "KozakStrength=Weak:type=uORF",
            "",
        )]);
        let (report, documented) = compare_default(&reference, &unfilled, &[UTR_STRUCTURAL]);
        assert_eq!(verdict(&report), "green", "{term}");
        assert_eq!(applied(&documented), once(UTR_STRUCTURAL), "{term}");
    }
    let words = |words: &[&str]| -> BTreeSet<String> {
        words.iter().map(|word| word.to_string()).collect()
    };
    assert_eq!(words(&class_words), words(STRUCTURAL_CLASS_WORDS));
    // The class words are the terms the structural comparison names a record's class by, the
    // breakpoint's left out.
    let so_terms: BTreeSet<String> = _SO_TERMS
        .iter()
        .map(|(_, term)| term.to_string())
        .filter(|term| term != "chromosome_breakpoint")
        .collect();
    assert_eq!(words(STRUCTURAL_CLASS_WORDS), so_terms);
    let tables = classes::tables();
    let own_fields: Vec<&str> = UTRANNOTATOR_CONSEQUENCE_FIELDS
        .iter()
        .chain(UTRANNOTATOR_ANNOTATION_FIELDS)
        .copied()
        .collect();
    assert_eq!(
        tables.field_classes_of(UTR_STRUCTURAL),
        Some(&own_fields[..])
    );
    // The class follows the plugin's other three, so a row one of them explains is theirs.
    let place = |klass: &str| {
        tables
            .field_classes()
            .iter()
            .position(|(held, _)| *held == klass)
            .expect("a field class")
    };
    assert_eq!(place(UTR_STRUCTURAL), place(UTR_PRECEDENCE) + 1);
    let (e, a) = (first_entry(&reference), first_entry(&unfilled));
    for name in ["5UTR_consequence", "5utr_consequence"] {
        assert!(
            accepts(
                UTR_STRUCTURAL,
                &e,
                name,
                &e.fields["5UTR_consequence"],
                "-",
                Some(&a)
            ),
            "{name}"
        );
    }
    for name in ["5UTR_annotation", "5utr_annotation"] {
        assert!(
            accepts(
                UTR_STRUCTURAL,
                &e,
                name,
                &e.fields["5UTR_annotation"],
                ABSENT,
                Some(&a)
            ),
            "{name}"
        );
    }
    assert!(!accepts(
        UTR_STRUCTURAL,
        &e,
        "5UTR_consequence",
        &e.fields["5UTR_consequence"],
        "-",
        None
    ));
}

/// The item-order class takes a NearestExonJB whose two values carry the same exon items in
/// different orders, counted once per row; without the class the row is a difference.
#[test]
fn the_exon_item_order_class_accepts_the_same_items_in_another_order() {
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    let reference = nearest_exon_jb_text(&[second, first]);
    let candidate = nearest_exon_jb_text(&[first, second]);
    let (report, documented) = compare_default(&reference, &candidate, &[ITEM_ORDER]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(ITEM_ORDER));
    assert!(masked_once(&report, "NearestExonJB", ITEM_ORDER));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    // The same value twice is no reordering.
    let same = format!("{second},{first}");
    assert!(!reordered_items(&same, &same));
    // Three items in another order, one of them written twice on both sides.
    let third = "ENSE00000000042|3|end|90";
    let (report, documented) = compare_default(
        &nearest_exon_jb_text(&[third, first, second, first]),
        &nearest_exon_jb_text(&[first, first, second, third]),
        &[ITEM_ORDER],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(ITEM_ORDER));
}

/// A reordered value in which one item is another item, in its exon, its distance, its side
/// or its length, differs as a multiset and stays a difference.
#[test]
fn the_exon_item_order_class_refuses_an_item_replaced() {
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    let reference = nearest_exon_jb_text(&[second, first]);
    for replaced in [
        "ENSE00000000032|3|start|117",
        "ENSE00000000031|4|start|117",
        "ENSE00000000031|3|end|117",
        "ENSE00000000031|3|start|118",
    ] {
        assert!(
            refused(
                &reference,
                &nearest_exon_jb_text(&[replaced, second]),
                ITEM_ORDER,
                "NearestExonJB"
            ),
            "{replaced}"
        );
        assert!(
            refused(
                &nearest_exon_jb_text(&[second, replaced]),
                &nearest_exon_jb_text(&[first, second]),
                ITEM_ORDER,
                "NearestExonJB"
            ),
            "{replaced}"
        );
    }
}

/// A value that lacks one of the other's items, on either side, and a value one engine does
/// not write at all, differ as a multiset and stay differences.
#[test]
fn the_exon_item_order_class_refuses_an_item_missing() {
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    let third = "ENSE00000000042|3|end|90";
    for (reference, candidate) in [
        (vec![second, first], vec![first]),
        (vec![second], vec![first, second]),
        (vec![third, second, first], vec![first, second]),
        (vec![second, first], vec![]),
    ] {
        assert!(
            refused(
                &nearest_exon_jb_text(&reference),
                &nearest_exon_jb_text(&candidate),
                ITEM_ORDER,
                "NearestExonJB"
            ),
            "{candidate:?}"
        );
    }
    let (report, documented) = compare_default(
        &nearest_exon_jb_text(&[]),
        &nearest_exon_jb_text(&[first, second]),
        &[ITEM_ORDER],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
}

/// A value that writes one of the other's items once more, on either side, or twice in place
/// of the other item, differs as a multiset and stays a difference.
#[test]
fn the_exon_item_order_class_refuses_an_item_duplicated() {
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    let reference = nearest_exon_jb_text(&[second, first]);
    for items in [
        &[first, second, second][..],
        &[second, first, first],
        &[first, first],
        &[second, second],
    ] {
        assert!(
            refused(
                &reference,
                &nearest_exon_jb_text(items),
                ITEM_ORDER,
                "NearestExonJB"
            ),
            "{items:?}"
        );
    }
    let twice = nearest_exon_jb_text(&[second, first, first]);
    assert!(refused(
        &twice,
        &nearest_exon_jb_text(&[first, second]),
        ITEM_ORDER,
        "NearestExonJB"
    ));
    assert!(refused(
        &twice,
        &nearest_exon_jb_text(&[first, second, second]),
        ITEM_ORDER,
        "NearestExonJB"
    ));
}

/// The class reads NearestExonJB alone: the same two texts under any other field's name are
/// refused, and UTRAnnotator's reordered 5UTR_annotation stays a difference under it.
#[test]
fn the_exon_item_order_class_refuses_the_same_reordering_on_another_field_or_plugin() {
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    assert_eq!(
        classes::tables().field_classes_of(ITEM_ORDER),
        Some(&["NearestExonJB"][..])
    );
    let entry = first_entry(&nearest_exon_jb_text(&[second, first]));
    let reference_text = format!("{second},{first}");
    let candidate_text = format!("{first},{second}");
    for name in ["NearestExonJB", "nearestexonjb"] {
        assert!(
            accepts(
                ITEM_ORDER,
                &entry,
                name,
                &reference_text,
                &candidate_text,
                None
            ),
            "{name}"
        );
    }
    for name in [
        "5UTR_annotation",
        "5UTR_consequence",
        "DOMAINS",
        "FREQS",
        "ZYG",
        "TSSDistance",
        "SpliceAI_pred",
        "NearestExonJB_x",
        "Extra",
    ] {
        assert!(
            !accepts(
                ITEM_ORDER,
                &entry,
                name,
                &reference_text,
                &candidate_text,
                None
            ),
            "{name}"
        );
    }
    assert!(refused(
        &utr_annotation_text(&hash_order()),
        &utr_annotation_text(&key_order()),
        ITEM_ORDER,
        "5UTR_annotation"
    ));
}

/// With both order classes named, DOMAINS items planted in another order are counted as a
/// differing DOMAINS and nothing is set aside for them, alone or on a row whose
/// NearestExonJB the class takes.
#[test]
fn a_third_fields_reordering_stays_a_difference_beside_the_item_order_class() {
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    let planted = [DOMAINS_ITEMS[1], DOMAINS_ITEMS[0]];
    let both = [ITEM_ORDER, PAIR_ORDER];
    let (report, documented) = compare_default(
        &nearest_exon_jb_row(&[first, second], &DOMAINS_ITEMS),
        &nearest_exon_jb_row(&[first, second], &planted),
        &both,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
    assert_eq!(differing(&report, "DOMAINS"), 1);
    let (report, documented) = compare_default(
        &nearest_exon_jb_row(&[second, first], &DOMAINS_ITEMS),
        &nearest_exon_jb_row(&[first, second], &planted),
        &both,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), once(ITEM_ORDER));
    assert_eq!(differing(&report, "DOMAINS"), 1);
    assert_eq!(differing(&report, "NearestExonJB"), 0);
}

/// A text format joins a list with `,` and the VCF format with `&`, never both: a value
/// carrying both is split by neither into the items its engine wrote, so the two classes
/// refuse it, whichever separator would make its parts agree.
#[test]
fn the_order_classes_refuse_a_value_carrying_both_list_separators() {
    assert!(reordered_items("b,a", "a,b") && reordered_items("b&a", "a&b"));
    assert!(!reordered_items("a,b&c", "c&a,b") && !reordered_items("b&c,a", "a,b&c"));
    assert!(!reordered_pair_groups("x=1,b=2&a=3", "a=3&x=1,b=2"));
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    let entry = first_entry(&nearest_exon_jb_text(&[first, second]));
    assert!(!accepts(
        ITEM_ORDER,
        &entry,
        "NearestExonJB",
        "a,b&c",
        "c&a,b",
        None
    ));
    assert!(!accepts(
        PAIR_ORDER,
        &entry,
        "5UTR_annotation",
        "x=1,b=2&a=3",
        "a=3&x=1,b=2",
        None
    ));
}

/// The pair-order class takes a 5UTR_annotation whose groups carry the same pairs in
/// different orders, the candidate's pairs in byte order of their keys, counted once per
/// row; without the class the row is a difference.
#[test]
fn the_annotation_pair_order_class_accepts_the_same_pairs_in_another_order() {
    let reference = utr_annotation_text(&hash_order());
    let candidate = utr_annotation_text(&key_order());
    let (report, documented) = compare_default(&reference, &candidate, &[PAIR_ORDER]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(PAIR_ORDER));
    assert!(masked_once(&report, "5UTR_annotation", PAIR_ORDER));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    // Two groups: the groups in another order too, each with its own pairs.
    let second = [("DistanceToCDS", "80"), ("type", "uORF")];
    let two_groups_reference = format!(
        "{},{}",
        utr_annotation(&UTR_ANNOTATION_HASH_ORDER, &second),
        hash_order()
    );
    let two_groups_candidate = format!(
        "{},{}",
        key_order(),
        utr_annotation(&UTR_ANNOTATION_KEY_ORDER, &second)
    );
    let (report, documented) = compare_default(
        &utr_annotation_text(&two_groups_reference),
        &utr_annotation_text(&two_groups_candidate),
        &[PAIR_ORDER],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(PAIR_ORDER));
}

/// A reordered value in which one pair is another pair, in its value or its key, or in which
/// a pair stands in another group, differs as a multiset and stays a difference.
#[test]
fn the_annotation_pair_order_class_refuses_a_pair_replaced() {
    let reference = utr_annotation_text(&hash_order());
    let mut renamed: Vec<&str> = UTR_ANNOTATION_KEY_ORDER
        .iter()
        .map(|key| if *key == "type" { "ref_type" } else { *key })
        .collect();
    renamed.sort_unstable();
    for annotation in [
        utr_annotation(&UTR_ANNOTATION_KEY_ORDER, &[("DistanceToCDS", "24")]),
        utr_annotation(&UTR_ANNOTATION_KEY_ORDER, &[("type", "uORF")]),
        utr_annotation(&renamed, &[("ref_type", "OutOfFrame_oORF")]),
    ] {
        assert!(
            refused(
                &reference,
                &utr_annotation_text(&annotation),
                PAIR_ORDER,
                "5UTR_annotation"
            ),
            "{annotation}"
        );
    }
    let second = [("DistanceToCDS", "80"), ("type", "uORF")];
    let two_groups = format!(
        "{},{}",
        utr_annotation(&UTR_ANNOTATION_HASH_ORDER, &second),
        hash_order()
    );
    let moved = format!(
        "{},{}",
        utr_annotation(&UTR_ANNOTATION_KEY_ORDER, &[("type", "uORF")]),
        utr_annotation(&UTR_ANNOTATION_KEY_ORDER, &[("DistanceToCDS", "80")])
    );
    assert!(refused(
        &utr_annotation_text(&two_groups),
        &utr_annotation_text(&moved),
        PAIR_ORDER,
        "5UTR_annotation"
    ));
}

/// A value that lacks one of the other's pairs or groups, on either side, and a value one
/// engine does not write, differ as a multiset and stay differences.
#[test]
fn the_annotation_pair_order_class_refuses_a_pair_missing() {
    let reference = utr_annotation_text(&hash_order());
    let field = "5UTR_annotation";
    assert!(refused(
        &reference,
        &utr_annotation_text(&utr_annotation(&UTR_ANNOTATION_KEY_ORDER[..6], &[])),
        PAIR_ORDER,
        field
    ));
    assert!(refused(
        &utr_annotation_text(&utr_annotation(&UTR_ANNOTATION_HASH_ORDER[1..], &[])),
        &utr_annotation_text(&key_order()),
        PAIR_ORDER,
        field
    ));
    let second = utr_annotation(
        &UTR_ANNOTATION_KEY_ORDER,
        &[("DistanceToCDS", "80"), ("type", "uORF")],
    );
    assert!(refused(
        &utr_annotation_text(&format!("{},{second}", hash_order())),
        &utr_annotation_text(&key_order()),
        PAIR_ORDER,
        field
    ));
    assert!(refused(
        &reference,
        &utr_annotation_text("-"),
        PAIR_ORDER,
        field
    ));
    let (report, documented) = compare_default(
        &utr_annotation_text("-"),
        &utr_annotation_text(&key_order()),
        &[PAIR_ORDER],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
}

/// A value that writes one of the other's pairs, or its one group, once more, on either
/// side, differs as a multiset and stays a difference.
#[test]
fn the_annotation_pair_order_class_refuses_a_pair_duplicated() {
    let (hash_order, key_order) = (hash_order(), key_order());
    let reference = utr_annotation_text(&hash_order);
    let field = "5UTR_annotation";
    let mut twice: Vec<&str> = UTR_ANNOTATION_KEY_ORDER.to_vec();
    twice.push("type");
    twice.sort_unstable();
    assert!(refused(
        &reference,
        &utr_annotation_text(&utr_annotation(&twice, &[])),
        PAIR_ORDER,
        field
    ));
    let mut hash_twice: Vec<&str> = UTR_ANNOTATION_HASH_ORDER.to_vec();
    hash_twice.push("type");
    assert!(refused(
        &utr_annotation_text(&utr_annotation(&hash_twice, &[])),
        &utr_annotation_text(&key_order),
        PAIR_ORDER,
        field
    ));
    assert!(refused(
        &reference,
        &utr_annotation_text(&format!("{key_order},{key_order}")),
        PAIR_ORDER,
        field
    ));
    assert!(refused(
        &utr_annotation_text(&format!("{hash_order},{hash_order}")),
        &utr_annotation_text(&key_order),
        PAIR_ORDER,
        field
    ));
}

/// The class is the reference's changing order against the candidate's fixed one, the byte
/// order of the keys: a candidate value out of key order, with the same pairs, stays a
/// difference, as does the comparison with its two sides exchanged.
#[test]
fn the_annotation_pair_order_class_refuses_a_candidate_out_of_key_order() {
    let reference = utr_annotation_text(&hash_order());
    let candidate = utr_annotation_text(&key_order());
    let mut reversed = UTR_ANNOTATION_KEY_ORDER;
    reversed.reverse();
    let mut out_of_order = UTR_ANNOTATION_KEY_ORDER;
    out_of_order.rotate_left(1);
    for keys in [reversed, out_of_order] {
        assert!(
            refused(
                &reference,
                &utr_annotation_text(&utr_annotation(&keys, &[])),
                PAIR_ORDER,
                "5UTR_annotation"
            ),
            "{keys:?}"
        );
    }
    assert!(refused(
        &candidate,
        &reference,
        PAIR_ORDER,
        "5UTR_annotation"
    ));
    assert!(!reordered_pair_groups("b=1:a=2", "b=1:a=2"));
    assert!(!reordered_pair_groups("a=2:b=1", "b=1:a=2"));
    // The same value twice, in key order, is no reordering.
    assert!(!reordered_pair_groups("a=2:b=1", "a=2:b=1"));
    // Key order is the keys' own: `k` before `k2`, where the whole pairs sort the other way.
    assert!(reordered_pair_groups("k2=y:k=x", "k=x:k2=y"));
    assert!(!reordered_pair_groups("k=x:k2=y", "k2=y:k=x"));
    // A text that is not key=value pairs is no group of pairs, whatever its order.
    assert!(!reordered_pair_groups("b:a", "a:b"));
    assert!(reordered_pair_groups("b=1:a=2", "a=2:b=1"));
}

/// The class reads 5UTR_annotation alone: the same two texts under any other field's name are
/// refused, NearestExonJB's reordered items stay a difference under it, and a reordered value
/// under another 5UTR_consequence leaves the consequence a difference.
#[test]
fn the_annotation_pair_order_class_refuses_the_same_reordering_on_another_field_or_plugin() {
    let (hash_order, key_order) = (hash_order(), key_order());
    assert_eq!(
        classes::tables().field_classes_of(PAIR_ORDER),
        Some(&["5UTR_annotation", "5utr_annotation"][..])
    );
    let entry = first_entry(&utr_annotation_text(&hash_order));
    assert!(accepts(
        PAIR_ORDER,
        &entry,
        "5UTR_annotation",
        &hash_order,
        &key_order,
        None
    ));
    for name in [
        "5UTR_consequence",
        "NearestExonJB",
        "DOMAINS",
        "FREQS",
        "Existing_uORFs",
        "SpliceVault_site_type",
        "5UTR_annotation_x",
        "Extra",
    ] {
        assert!(
            !accepts(PAIR_ORDER, &entry, name, &hash_order, &key_order, None),
            "{name}"
        );
    }
    let (first, second) = NEAREST_EXON_JB_ITEMS;
    assert!(refused(
        &nearest_exon_jb_text(&[second, first]),
        &nearest_exon_jb_text(&[first, second]),
        PAIR_ORDER,
        "NearestExonJB"
    ));
    let other_effect = utr_annotation_row(&key_order, UTR_STOP_LOSS, &[]);
    let (report, _) = compare_default(
        &utr_annotation_text(&hash_order),
        &other_effect,
        &[PAIR_ORDER],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "5UTR_consequence"), 1);
}

/// With both order classes named, DOMAINS items planted in another order are counted as a
/// differing DOMAINS and nothing is set aside for them, alone or on a row whose
/// 5UTR_annotation the class takes.
#[test]
fn a_third_fields_reordering_stays_a_difference_beside_the_pair_order_class() {
    let (hash_order, key_order) = (hash_order(), key_order());
    let planted = [DOMAINS_ITEMS[1], DOMAINS_ITEMS[0]];
    let both = [ITEM_ORDER, PAIR_ORDER];
    let (report, documented) = compare_default(
        &utr_annotation_row(&key_order, UTR_LOSS, &DOMAINS_ITEMS),
        &utr_annotation_row(&key_order, UTR_LOSS, &planted),
        &both,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
    assert_eq!(differing(&report, "DOMAINS"), 1);
    let (report, documented) = compare_default(
        &utr_annotation_row(&hash_order, UTR_LOSS, &DOMAINS_ITEMS),
        &utr_annotation_row(&key_order, UTR_LOSS, &planted),
        &both,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), once(PAIR_ORDER));
    assert_eq!(differing(&report, "DOMAINS"), 1);
    assert_eq!(differing(&report, "5UTR_annotation"), 0);
}

/// What a class says of one of the plugin's fields on a pair of default-format rows.
fn offset_grid_verdict(klass: &str, reference: &str, candidate: &str, field: &str) -> bool {
    let e = first_entry(&default_text(&[reference]));
    let a = first_entry(&default_text(&[candidate]));
    accepts(
        klass,
        &e,
        field,
        &field_of(Some(&e), &[field]),
        &field_of(Some(&a), &[field]),
        Some(&a),
    )
}

/// The protein gate asks the same HGVSc and HGVS_OFFSET of both rows, or a reference whose
/// insertion is left where its shift stops beside a candidate that describes it further 3':
/// the protein notation the candidate alone writes is the class there too. A candidate
/// HGVS_OFFSET that is the reference's beside another HGVSc stays a difference.
#[test]
fn the_protein_span_gate_reads_a_reference_left_where_its_shift_stops() {
    let klass = "hgvsp_shifted_indel_written_span_gate";
    let allele = "ACG".repeat(333);
    let columns = format!(
        "intron_end_insertion\t21:25001998-25001999\t{allele}\tENSG1\tENST00000000001\tTranscript\tsplice_acceptor_variant\t-\t-\t-\t-\t-\t-"
    );
    let extra = "IMPACT=HIGH;STRAND=1";
    let reference = format!(
        "{columns}\t{extra};HGVSc=ENST00000000001.1:c.251_252ins{};HGVS_OFFSET=2",
        turned_allele(&allele, 2, true)
    );
    let candidate = |offset: &str| {
        format!(
            "{columns}\t{extra};HGVSc=ENST00000000001.1:c.254_1252dup;HGVSp=ENSP00000000001.1:p.Ala84ValfsTer?;HGVS_OFFSET={offset}"
        )
    };
    assert!(offset_grid_verdict(
        klass,
        &reference,
        &candidate("1001"),
        "HGVSp"
    ));
    assert!(!offset_grid_verdict(
        klass,
        &reference,
        &candidate("2"),
        "HGVSp"
    ));
    let with_protein = reference.replace(
        ";HGVS_OFFSET=",
        ";HGVSp=ENSP00000000001.1:p.Ala84ValfsTer?;HGVS_OFFSET=",
    );
    assert!(offset_grid_verdict(
        klass,
        &reference,
        &with_protein,
        "HGVSp"
    ));
}

/// A row on which the candidate writes nothing stays the abutting class's: the reference's
/// pair of a notation the row never prints against no candidate field.
#[test]
fn the_bases_written_outside_class_leaves_the_abutting_class_its_rows() {
    let columns = "cell0\t2:7000000-7000001\tT\tENSG1\tENST00000000071\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-";
    let candidate = format!("{columns}\tIMPACT=MODIFIER;DISTANCE=0;STRAND=1");
    let reference = format!("{candidate};HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0");
    for field in ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"] {
        assert!(
            !offset_grid_verdict(
                "hgvs_added_bases_written_outside_the_transcript",
                &reference,
                &candidate,
                field
            ),
            "{field}"
        );
        assert!(
            offset_grid_verdict(
                "hgvs_intron_offset_abutting_insertion",
                &reference,
                &candidate,
                field
            ),
            "{field}"
        );
    }
}

/// A gnomAD_ field stays a difference when the reference's value is cut anywhere but before
/// the candidate's first `=`, when the candidate's value holds no `=`, when the reference
/// writes the longer value, and when either side leaves the field empty.
#[test]
fn the_second_equals_class_refuses_any_other_pair_of_values() {
    let klass = "gnomadmt_info_value_second_equals";
    let whole = gnomad_vep_whole();
    for (reference, candidate) in [
        (&GNOMAD_VEP_CUT[..GNOMAD_VEP_CUT.len() - 3], whole.as_str()),
        (format!("{GNOMAD_VEP_CUT}=").as_str(), whole.as_str()),
        (GNOMAD_VEP_CUT, format!("{GNOMAD_VEP_CUT}|135").as_str()),
        (whole.as_str(), GNOMAD_VEP_CUT),
        ("", whole.as_str()),
        (GNOMAD_VEP_CUT, ""),
    ] {
        let (report, documented) = compare_default(
            &gnomad_vep_row(reference),
            &gnomad_vep_row(candidate),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{reference} {candidate}");
        assert_eq!(applied(&documented), nothing(), "{reference} {candidate}");
    }
    // Another plugin's field of the same shape is not the class's.
    let reference = default_text(&[default_row(
        "a",
        "MT:3441",
        "G",
        "ENST00000361390",
        "synonymous_variant",
        "IMPACT=LOW;STRAND=1;LoF_info=PERCENTILE:0.4",
    )]);
    assert_eq!(
        verdict_of(
            &reference,
            &reference.replace("PERCENTILE:0.4", "PERCENTILE:0.4=1"),
            &[klass]
        ),
        "red"
    );
    // A `*` of one base beside a two-base REF is an ALT of another length: the record is
    // chopped.
    let alternates = |alleles: &[&str]| -> Vec<String> {
        alleles.iter().map(|allele| allele.to_string()).collect()
    };
    assert_eq!(
        vcf_first_alternate("GC", &alternates(&["GT", "GA", "*"])),
        "T"
    );
    assert_eq!(vcf_first_alternate("C", &alternates(&["T", "*"])), "T");
}

/// A NearestExonJB stays a difference on a row 13 bases out, on a row without a DISTANCE,
/// where the reference alone fills it, and where both fill it differently.
#[test]
fn the_stretch_priming_class_refuses_a_row_farther_out_or_filled_by_the_reference() {
    let klass = "nearest_exon_jb_stretch_priming";
    let row: Row = [
        "nine_before_first_exon",
        "21:6354255",
        "T",
        "ENST00000624847",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=9;STRAND=1",
    ];
    let value = "ENSE00003757672|9|start|259";
    let farther = with_extra(row, "IMPACT=MODIFIER;DISTANCE=13;STRAND=1");
    let inside: Row = [
        row[0],
        row[1],
        row[2],
        row[3],
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    ];
    for (reference, candidate) in [
        (
            nearest_exon_row(farther, ""),
            nearest_exon_row(farther, "ENSE00003757672|13|start|259"),
        ),
        (
            nearest_exon_row(inside, ""),
            nearest_exon_row(inside, value),
        ),
        (nearest_exon_row(row, value), nearest_exon_row(row, "")),
        (
            nearest_exon_row(row, value),
            nearest_exon_row(row, "ENSE00003757672|9|start|260"),
        ),
    ] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
}

/// A NearestExonJB stays a difference on a row whose SHIFT_LENGTH is 0 or not printed, on a
/// row that is not a transcript's, and where both values name the same exons at another
/// distance.
#[test]
fn the_shifted_span_class_refuses_a_row_without_a_shift_or_naming_the_same_exons() {
    let klass = "nearest_exon_jb_shifted_span";
    let value = "ENSE00003584551|0|end|189";
    let unshifted = with_extra(SHIFT_PROBE, "IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=0");
    let unprinted = with_extra(SHIFT_PROBE, "IMPACT=HIGH;STRAND=1");
    let regulatory = default_text(&[default_row_of(
        SHIFT_PROBE[0],
        SHIFT_PROBE[1],
        SHIFT_PROBE[2],
        "ENSR00000141251",
        "regulatory_region_variant",
        "IMPACT=MODIFIER;SHIFT_LENGTH=1",
        "RegulatoryFeature",
        "ENSG1",
    )]);
    for (reference, candidate) in [
        (
            nearest_exon_row(unshifted, ""),
            nearest_exon_row(unshifted, value),
        ),
        (
            nearest_exon_row(unprinted, ""),
            nearest_exon_row(unprinted, value),
        ),
        (
            regulatory.clone(),
            regulatory.replace(
                "SHIFT_LENGTH=1",
                &format!("SHIFT_LENGTH=1;NearestExonJB={value}"),
            ),
        ),
        (
            nearest_exon_row(SHIFT_PROBE, value),
            nearest_exon_row(SHIFT_PROBE, "ENSE00003584551|1|end|189"),
        ),
    ] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
    // Two spellings of an empty value are not a value against none.
    let entry = first_entry(&nearest_exon_row(SHIFT_PROBE, ""));
    assert!(!accepts(klass, &entry, "NearestExonJB", "", "-", None));
    assert!(accepts(klass, &entry, "NearestExonJB", "", value, None));
}

/// On the row of a one-allele record a value the candidate lacks and a value it writes wrong
/// both stay differences under a plugin class: each is counted under its field and the class
/// is applied to nothing.
fn a_planted_value_stays_a_difference(
    klass: &str,
    consequence: &str,
    other_pairs: &str,
    field: &str,
    value: &str,
    wrong_value: &str,
) {
    let row = |pair: &str| {
        let extra = ["IMPACT=MODERATE;STRAND=1", other_pairs, pair]
            .into_iter()
            .filter(|part| !part.is_empty())
            .collect::<Vec<_>>()
            .join(";");
        default_text(&[default_row(
            "21_100_A/G",
            "21:100",
            "G",
            "ENST1",
            consequence,
            &extra,
        )])
    };
    let reference = row(&format!("{field}={value}"));
    for candidate in [row(""), row(&format!("{field}={wrong_value}"))] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
        assert_eq!(differing(&report, field), 1, "{candidate}");
    }
}

#[test]
fn the_alphamissense_sibling_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "alphamissense_sibling_allele_rows",
        "missense_variant",
        "",
        "am_pathogenicity",
        "0.1079",
        "0.9000",
    );
}

#[test]
fn the_revel_sibling_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "revel_sibling_allele_rows",
        "missense_variant",
        "",
        "REVEL",
        "0.025",
        "0.027",
    );
}

#[test]
fn the_revel_other_amino_acid_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "revel_other_reference_amino_acid_rows",
        "missense_variant",
        "",
        "REVEL",
        "0.025",
        "0.027",
    );
}

#[test]
fn the_splicevault_sibling_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "splicevault_sibling_allele_rows",
        "splice_donor_variant",
        "",
        "SpliceVault_site_type",
        "donor",
        "acceptor",
    );
}

#[test]
fn the_second_equals_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "gnomadmt_info_value_second_equals",
        "synonymous_variant",
        "",
        "gnomAD_vep",
        "G|synonymous_variant|p.Leu45",
        "G|synonymous_variant|p.Leu46",
    );
}

#[test]
fn the_ancestral_allele_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "loftee_anc_allele_multiallelic_record",
        "stop_gained",
        "LoF=LC",
        "LoF_filter",
        "END_TRUNC",
        "END_TRUNC,ANC_ALLELE",
    );
}

#[test]
fn the_first_alternate_donor_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "loftee_first_alt_donor_filter",
        "splice_donor_variant",
        "LoF=LC",
        "LoF_filter",
        "GC_TO_GT_DONOR",
        "END_TRUNC",
    );
}

#[test]
fn the_stretch_priming_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "nearest_exon_jb_stretch_priming",
        "upstream_gene_variant",
        "DISTANCE=9",
        "NearestExonJB",
        "ENSE00003757672|9|start|259",
        "ENSE00003757672|9|start|260",
    );
}

#[test]
fn the_shifted_span_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "nearest_exon_jb_shifted_span",
        "splice_donor_variant",
        "",
        "NearestExonJB",
        "ENSE00003584551|0|end|189",
        "ENSE00003675684|0|end|189",
    );
}

#[test]
fn the_id_less_exon_fold_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "nearest_exon_jb_id_less_exon_fold",
        "intron_variant",
        "",
        "NearestExonJB",
        "ENSE00003584551|0|end|189",
        "ENSE00003584551|0|end|188",
    );
}

#[test]
fn the_reverse_strand_rotation_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "downstream_reverse_strand_rotation",
        "frameshift_variant",
        "",
        "DownstreamProtein",
        "SS",
        "SP",
    );
}

#[test]
fn the_chopped_multiallelic_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "spliceai_chopped_multiallelic",
        "intron_variant",
        "",
        "SpliceAI_pred",
        "GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1",
        "GENE1|0.11|0.20|0.30|0.40|-4|-3|2|1",
    );
}

#[test]
fn the_mitochondrial_spelling_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        MITOCHONDRIAL,
        "missense_variant",
        "",
        "ClinPred",
        "0.5",
        "0.6",
    );
}

#[test]
fn the_zero_minimum_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "loeuf_zero_minimum",
        "intron_variant",
        "",
        "LOEUF",
        "0.5",
        "0.6",
    );
}

#[test]
fn the_gerp_summary_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        GERP,
        "stop_gained",
        "LoF=HC",
        "LoF_info",
        "PERCENTILE:0.5,GERP_DIST:-10.5,BP_DIST:60",
        "PERCENTILE:0.5,GERP_DIST:-12.5,BP_DIST:60",
    );
}

#[test]
fn the_haplotype_cross_product_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "gwas_haplotype_cross_product",
        "intron_variant",
        "",
        "GWAS_study",
        "GCST1",
        "GCST2",
    );
}

#[test]
fn the_interaction_whitespace_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "gwas_interaction_first_snp_allele_whitespace",
        "intron_variant",
        "",
        "GWAS_study",
        "GCST1",
        "GCST2",
    );
}

#[test]
fn the_gained_start_window_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "utrannotator_gained_start_window",
        "5_prime_UTR_variant",
        "5UTR_annotation=CapDistanceToStart=64:type=uORF",
        "5UTR_consequence",
        UTR_GAIN,
        UTR_LOSS,
    );
}

#[test]
fn the_consequence_precedence_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        UTR_PRECEDENCE,
        "5_prime_UTR_variant",
        "5UTR_annotation=type=uORF",
        "5UTR_consequence",
        UTR_LOSS,
        UTR_STOP_LOSS,
    );
}

#[test]
fn the_written_span_gate_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "hgvs_shifted_indel_written_span_gate",
        "non_coding_transcript_exon_variant",
        "",
        "HGVS_IntronStartOffset",
        "0",
        "5",
    );
}

#[test]
fn the_clamped_deletion_span_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "hgvs_shifted_deletion_clamped_span",
        "non_coding_transcript_exon_variant",
        "HGVSc=ENST1.1:n.34_43del;HGVS_OFFSET=33",
        "HGVS_IntronStartOffset",
        "2",
        "1",
    );
}

#[test]
fn the_splice_family_swap_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(
        "splice_family_swap",
        "frameshift_variant",
        "",
        "DownstreamProtein",
        "SS",
        "SP",
    );
}

/// A class that reads none of the record's rows, the GWAS catalog's pairs and the seeded runs
/// answers as it does without them.
#[test]
fn a_class_that_reads_nothing_beside_the_values_answers_alike_with_nothing_beside_them() {
    let e = first_entry(&default_text(&[default_row(
        "rs1",
        "21:100",
        "A",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER",
    )]));
    for (value, accepted) in [("0", true), ("0.5", false)] {
        assert_eq!(
            accepts("loeuf_zero_minimum", &e, "LOEUF", ABSENT, value, None),
            accepted
        );
        assert_eq!(
            accepts_beside(
                "loeuf_zero_minimum",
                &e,
                "LOEUF",
                ABSENT,
                value,
                &Beside::default()
            ),
            accepted
        );
    }
}

/// `tabix_allele_case_insensitive_match` takes a field of one of its nine plugins that the
/// candidate alone fills, under the text formats' names and the JSON keys, and leaves a
/// difference where both engines fill the field, where the reference alone fills it and
/// where the field is another plugin's, a `gnomAD_` field of gnomADc among them.
#[test]
fn the_lower_case_table_row_class_takes_a_plugin_field_the_candidate_alone_fills() {
    let e = first_entry(&default_text(&[default_row(
        "rs1",
        "21:100",
        "A",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER",
    )]));
    let klass = "tabix_allele_case_insensitive_match";
    for field in [
        "am_class",
        "alphamissense",
        "AVI_PHRED",
        "avi_raw",
        "CADD_PHRED",
        "cadd_raw",
        "gnomAD_AF_hom",
        "GWAS_risk_allele",
        "gwas_study",
        "OpenTargets_gwasGeneId",
        "opentargets",
        "REVEL",
        "revel",
        "SpliceAI_pred_DS_AG",
        "spliceai",
        "SpliceVault_top_events",
        "splicevault",
    ] {
        assert!(accepts(klass, &e, field, ABSENT, "0.5", None), "{field}");
        assert!(!accepts(klass, &e, field, "0.4", "0.5", None), "{field}");
        assert!(!accepts(klass, &e, field, "0.5", ABSENT, None), "{field}");
    }
    for field in [
        "LOEUF",
        "BayesDel",
        "ClinPred",
        "PrimateAI",
        "Consequence",
        "gnomAD_mean_cov",
        "gnomad_mean_cov",
    ] {
        assert!(!accepts(klass, &e, field, ABSENT, "0.5", None), "{field}");
    }

    // The GWAS plugin writes a cell its catalog row leaves empty as the empty string: the
    // class's on a row the candidate alone annotates, and a difference on a row neither
    // annotates, on a row the reference annotates too and under another plugin's name.
    let annotated = |extra: &str| {
        first_entry(&default_text(&[default_row(
            "rs1",
            "21:100",
            "A",
            "ENST1",
            "intron_variant",
            extra,
        )]))
    };
    let candidate = annotated("IMPACT=MODIFIER;GWAS_odds_ratio=;GWAS_risk_allele=t");
    let reference = annotated("IMPACT=MODIFIER;GWAS_risk_allele=T");
    for field in ["GWAS_odds_ratio", "gwas_odds_ratio"] {
        assert!(
            accepts(klass, &e, field, ABSENT, "", Some(&candidate)),
            "{field}"
        );
        assert!(!accepts(klass, &e, field, ABSENT, "", Some(&e)), "{field}");
        assert!(!accepts(klass, &e, field, ABSENT, "", None), "{field}");
        assert!(
            !accepts(klass, &reference, field, ABSENT, "", Some(&candidate)),
            "{field}"
        );
    }
    assert!(!accepts(
        klass,
        &e,
        "CADD_PHRED",
        ABSENT,
        "",
        Some(&candidate)
    ));
}

/// One row of a deletion over the stop codon of ENST00000635895 as a default-format output,
/// its Extra text `extra`.
fn stop_codon_text(terms: &str, extra: &str) -> String {
    default_text(&[format!(
        "1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{terms}\t7055-7059\t6978-?\t-\t-\t-\t-\t{extra}"
    )])
}

/// A row class of the following cells takes, on a row its predicate admits, the cells read
/// from the two sets of terms: the row is counted once under the class, with its Consequence.
/// A cell the rule does not read stays a difference on that row, the class reads no cell of a
/// row its predicate refuses, and a caller that applies the class by its name reads the same
/// rule.
#[test]
fn a_row_class_of_the_following_cells_reads_them_on_its_own_rows_alone() {
    let klass = "stop_retained_three_base_stop_codon_test";
    assert!(classes::tables().in_term_pair_following_classes(klass));
    let retained = "stop_retained_variant,3_prime_UTR_variant";
    let unknown = "coding_sequence_variant,3_prime_UTR_variant";
    let reference = stop_codon_text(retained, "IMPACT=LOW;STRAND=-1");
    let candidate = stop_codon_text(unknown, "IMPACT=MODIFIER;STRAND=-1");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    assert!(masked_once(&report, "Consequence", klass));
    assert!(masked_once(&report, "Extra", klass));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let (report, documented) = compare_default(
        &reference,
        &candidate.replace("STRAND=-1", "STRAND=1"),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), once(klass));
    assert_eq!(differing(&report, "STRAND"), 1);
    assert!(at(&report, &["columns", "Extra"]).get("masked").is_none());
    // A pair of terms the class's predicate refuses.
    let other = (
        reference.replace("stop_retained_variant", "stop_lost"),
        candidate.replace("stop_retained_variant", "stop_lost"),
    );
    let (report, documented) = compare_default(&other.0, &other.1, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
    assert_eq!(differing(&report, "Consequence"), 1);
    assert_eq!(differing(&report, "Extra"), 1);
    let (e, a) = (first_entry(&reference), first_entry(&candidate));
    assert!(accepts(klass, &e, "IMPACT", "LOW", "MODIFIER", Some(&a)));
    assert!(!accepts(klass, &e, "IMPACT", "LOW", "HIGH", Some(&a)));
    assert!(!accepts(klass, &e, "IMPACT", "LOW", "MODIFIER", None));
    assert!(!accepts(klass, &e, "STRAND", "-1", "1", Some(&a)));
    let (e, a) = (first_entry(&other.0), first_entry(&other.1));
    assert!(!accepts(klass, &e, "IMPACT", "LOW", "MODIFIER", Some(&a)));
}

/// The record MT 3460 G A against a made-up table naming the contig chrM: the reference's
/// plugin finds nothing and the candidate writes the table's score. A record spelt chrM or
/// chrMT, the VCF line's CHROM and the JSON object's seq_region_name are read the same way.
#[test]
fn the_mitochondrial_spelling_class_accepts_a_plugin_field_the_candidate_alone_fills() {
    let klass = MITOCHONDRIAL;
    let bare = mitochondrial_row("");
    let scored = mitochondrial_row("ClinPred=0.5");
    let (report, documented) = compare_default(&bare, &scored, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    assert!(masked_once(&report, "ClinPred", klass));
    assert_eq!(verdict_of(&bare, &scored, &[]), "red");
    let (report, documented) = compare_default(
        &bare,
        &mitochondrial_row(MITOCHONDRIAL_ALPHAMISSENSE_FIELDS),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    for contig in ["chrM", "chrMT"] {
        let location = format!("{contig}:3460");
        assert_eq!(
            verdict_of(
                &mitochondrial_row_at("", &location),
                &mitochondrial_row_at("gnomAD_mean_cov=31.2", &location),
                &[klass]
            ),
            "green",
            "{contig}"
        );
    }
    let inputs = vec!["MT\t3460\trs199476118\tG\tA\t.\t.\t.".to_string()];
    let vcf = vcf_text(
        "Allele|Consequence|Feature|Feature_type|ClinPred",
        &inputs[0],
        "A|missense_variant|ENST00000361390|Transcript|",
    );
    let vcf_scored = vcf.replace("Transcript|\n", "Transcript|0.5\n");
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&vcf, &inputs),
            &vcf_entries(&vcf_scored, &inputs)
        ),
        green_once(klass)
    );
    // A text row names its chromosome itself; a VCF entry takes it from its record's line.
    let entries = vcf_entries(&vcf, &inputs);
    let (record_line, csq_entry) = (&entries[0], &entries[1]);
    assert!(accepts(
        klass,
        &first_entry(&bare),
        "ClinPred",
        ABSENT,
        "0.5",
        None
    ));
    assert!(!accepts(klass, csq_entry, "ClinPred", "", "0.5", None));
    assert!(accepts_beside(
        klass,
        csq_entry,
        "ClinPred",
        "",
        "0.5",
        &Beside {
            record_rows: Some(&[record_line, csq_entry][..]),
            ..Beside::default()
        }
    ));
    let object = |plugin: &str| {
        format!(
            "{{\"input\": \"MT\\t3460\\trs199476118\\tG\\tA\\t.\\t.\\t.\", \"id\": \"rs199476118\", \"seq_region_name\": \"MT\", \"most_severe_consequence\": \"missense_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"A\", \"transcript_id\": \"ENST00000361390\", \"consequence_terms\": [\"missense_variant\"]{plugin}}}]}}\n"
        )
    };
    let alphamissense =
        ", \"alphamissense\": {\"am_class\": \"likely_pathogenic\", \"am_pathogenicity\": 0.624}";
    assert_eq!(
        under_class(
            klass,
            &json_entries(&object(""), &inputs),
            &json_entries(&object(alphamissense), &inputs)
        ),
        green_once(klass)
    );
}

/// A plugin field the candidate alone fills stays a difference on another chromosome and on
/// a record spelt M; on a mitochondrial record a gnomADMt field, a field of a plugin that
/// reads no tabix file, a field the reference alone fills and one both fill stay differences.
#[test]
fn the_mitochondrial_spelling_class_refuses_another_contig_plugin_or_direction() {
    let scored = "ClinPred=0.5";
    for (reference, candidate) in [
        (
            mitochondrial_row_at("", "21:3460"),
            mitochondrial_row_at(scored, "21:3460"),
        ),
        (
            mitochondrial_row_at("", "M:3460"),
            mitochondrial_row_at(scored, "M:3460"),
        ),
        (
            mitochondrial_row(""),
            mitochondrial_row("gnomAD_hap_AF_hom=0.01"),
        ),
        (mitochondrial_row(""), mitochondrial_row("LoFtool=0.5")),
        (mitochondrial_row(scored), mitochondrial_row("")),
        (mitochondrial_row(scored), mitochondrial_row("ClinPred=0.6")),
    ] {
        let (report, documented) = compare_default(&reference, &candidate, &[MITOCHONDRIAL]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
}

/// Beside a manifest naming GRCh38, where the four contig names are one sequence, a field of
/// any tabix-backed plugin that the candidate alone fills on a record spelt MT is set aside;
/// a gnomADMt field, a field of a plugin that reads no tabix file and a record spelt M stay
/// differences.
#[test]
fn the_mitochondrial_spelling_class_takes_every_tabix_backed_plugin_on_grch38() {
    let manifest = Some("{\"assembly\": \"GRCh38\"}");
    for fields in
        std::iter::once(MITOCHONDRIAL_ALPHAMISSENSE_FIELDS).chain(MITOCHONDRIAL_PLUGIN_FIELDS)
    {
        let (report, documented) = compare_mitochondrial_rows(manifest, fields, "MT:3460");
        assert_eq!(documented.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(verdict(&report), "green", "{fields}");
        assert_eq!(applied(&documented), once(MITOCHONDRIAL), "{fields}");
    }
    for (fields, location) in [
        ("gnomAD_hap_AF_hom=0.01", "MT:3460"),
        ("LoFtool=0.5", "MT:3460"),
        ("ClinPred=0.5", "M:3460"),
    ] {
        let (report, documented) = compare_mitochondrial_rows(manifest, fields, location);
        assert_eq!(verdict(&report), "red", "{fields}");
        assert_eq!(applied(&documented), nothing(), "{fields}");
    }
}

/// Beside a manifest naming GRCh37, under every name read as that assembly, a field of any
/// tabix-backed plugin that the candidate alone fills on a record spelt MT, chrM or chrMT is
/// set aside as on GRCh38, in the JSON format as well: the comparison reads the record's
/// spelling and which output fills the field, and neither the name the plugin's file gives
/// the contig nor the bases its rows state. A gnomADMt field, a field of a plugin that reads
/// no tabix file, a record spelt M and a field the reference alone fills stay differences.
#[test]
fn the_mitochondrial_spelling_class_takes_every_tabix_backed_plugin_on_grch37() {
    let klass = MITOCHONDRIAL;
    for assembly in ["GRCh37", "grch37", "hg19", "37"] {
        let manifest = format!("{{\"assembly\": \"{assembly}\"}}");
        let mut refused = Vec::new();
        for fields in
            std::iter::once(MITOCHONDRIAL_ALPHAMISSENSE_FIELDS).chain(MITOCHONDRIAL_PLUGIN_FIELDS)
        {
            for location in ["MT:3460", "chrM:3460", "chrMT:3460"] {
                let (report, documented) =
                    compare_mitochondrial_rows(Some(&manifest), fields, location);
                if verdict(&report) != "green" || applied(&documented) != once(klass) {
                    refused.push(format!("{location} {fields}"));
                }
            }
        }
        assert_eq!(refused, Vec::<String>::new(), "{assembly}");
        for (fields, location) in [
            ("gnomAD_hap_AF_hom=0.01", "MT:3460"),
            ("LoFtool=0.5", "MT:3460"),
            ("ClinPred=0.5", "M:3460"),
        ] {
            let (report, documented) =
                compare_mitochondrial_rows(Some(&manifest), fields, location);
            assert_eq!(verdict(&report), "red", "{assembly} {fields}");
            assert_eq!(applied(&documented), nothing(), "{assembly} {fields}");
        }
    }
    let grch37 = json("{\"assembly\": \"GRCh37\"}");
    for fields in ["ClinPred=0.5", SPLICEAI_VALUE, "gnomAD_mean_cov=31.2"] {
        let mut documented = from_manifest(&grch37, &[klass]);
        let report = compare_texts(
            &mitochondrial_row(fields),
            &mitochondrial_row(""),
            &mut documented,
        );
        assert_eq!(verdict(&report), "red", "{fields}");
        assert_eq!(applied(&documented), nothing(), "{fields}");
    }
    let inputs = vec!["MT\t3460\trs199476118\tG\tA\t.\t.\t.".to_string()];
    let object = |plugin: &str| {
        format!(
            "{{\"input\": \"MT\\t3460\\trs199476118\\tG\\tA\\t.\\t.\\t.\", \"id\": \"rs199476118\", \"seq_region_name\": \"MT\", \"most_severe_consequence\": \"missense_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"A\", \"transcript_id\": \"ENST00000361390\", \"consequence_terms\": [\"missense_variant\"]{plugin}}}]}}\n"
        )
    };
    for plugin in [
        ", \"alphamissense\": {\"am_class\": \"likely_pathogenic\", \"am_pathogenicity\": 0.624}",
        ", \"clinpred\": 0.5",
    ] {
        let mut documented = from_manifest(&grch37, &[klass]);
        let report = compare(
            &json_entries(&object(""), &inputs),
            &json_entries(&object(plugin), &inputs),
            &mut documented,
        );
        assert_eq!(verdict(&report), "green", "{plugin}");
        assert_eq!(applied(&documented), once(klass), "{plugin}");
    }
}

/// A comparison beside no manifest, or beside one that names no assembly, takes the rows the
/// class takes beside one that names an assembly; the dispatch handed an assembly answers
/// as the dispatch handed none.
#[test]
fn the_mitochondrial_spelling_class_takes_every_tabix_backed_plugin_without_a_named_assembly() {
    let klass = MITOCHONDRIAL;
    for manifest in [
        None,
        Some("{\"release\": \"116\"}"),
        Some("{\"assembly\": null}"),
    ] {
        for fields in
            std::iter::once(MITOCHONDRIAL_ALPHAMISSENSE_FIELDS).chain(MITOCHONDRIAL_PLUGIN_FIELDS)
        {
            let (report, documented) = compare_mitochondrial_rows(manifest, fields, "MT:3460");
            assert_eq!(documented.assembly, None);
            assert_eq!(verdict(&report), "green", "{manifest:?} {fields}");
            assert_eq!(applied(&documented), once(klass), "{manifest:?} {fields}");
        }
    }
    let row = first_entry(&mitochondrial_row(""));
    for assembly in [None, Some("GRCh38"), Some("TAIR10"), Some("GRCh37")] {
        let beside = Beside {
            assembly,
            ..Beside::default()
        };
        for (field, value) in [
            ("ClinPred", "0.5"),
            ("SpliceAI_pred", "GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1"),
            ("gnomAD_mean_cov", "31.2"),
            ("am_class", "likely_pathogenic"),
        ] {
            assert!(
                accepts_beside(klass, &row, field, ABSENT, value, &beside),
                "{assembly:?} {field}"
            );
            assert!(
                !accepts_beside(klass, &row, field, value, ABSENT, &beside),
                "{assembly:?} {field}"
            );
        }
        assert!(
            !accepts_beside(klass, &row, "gnomAD_hap_AF_hom", ABSENT, "0.01", &beside),
            "{assembly:?}"
        );
    }
}

/// The command line reads the assembly from `--manifest`: a dbNSFP field the reference alone
/// fills on MT is set aside beside a manifest naming GRCh37 and is a difference beside one
/// naming GRCh38, and the row a ClinPred table fills on MT is set aside beside either.
#[test]
fn the_command_line_reads_the_assembly_of_the_manifest_file() {
    let directory = support::scratch("dispatch-assembly-of-the-manifest");
    let write = |name: &str, text: &str| {
        let path = directory.join(name);
        std::fs::write(&path, text).expect("the file is written");
        path
    };
    write("bare.txt", &mitochondrial_row(""));
    write("clinpred.txt", &mitochondrial_row("ClinPred=0.5"));
    write("dbnsfp.txt", &mitochondrial_row("dbNSFP_SIFT_score=0.01"));
    for (assembly, klass, reference, candidate, status) in [
        ("GRCh37", DBNSFP_NAMED_M, "dbnsfp.txt", "bare.txt", 0),
        ("GRCh38", DBNSFP_NAMED_M, "dbnsfp.txt", "bare.txt", 1),
        ("GRCh37", MITOCHONDRIAL, "bare.txt", "clinpred.txt", 0),
        ("GRCh38", MITOCHONDRIAL, "bare.txt", "clinpred.txt", 0),
    ] {
        let manifest = write(
            &format!("{assembly}.json"),
            &format!("{{\"assembly\": \"{assembly}\"}}"),
        );
        let path = |path: &std::path::Path| path.to_str().expect("a UTF-8 path").to_string();
        let ran = support::run(
            &directory,
            &[
                "--format",
                "default",
                "--expected",
                &path(&directory.join(reference)),
                "--actual",
                &path(&directory.join(candidate)),
                "--manifest",
                &path(&manifest),
                "--masks",
                klass,
                "--report",
                &path(&directory.join("report.json")),
            ],
        );
        assert_eq!(ran.code, status, "{assembly} {klass}: {}", ran.err());
    }
}

/// The record MT 3460 G A against a made-up dbNSFP table naming the contig M whose rows
/// count on another sequence than the assembly's: dbNSFP.pm:313 names the record M to the
/// table and :335-338 match the row at its position, so the reference writes the row's
/// value; the candidate reads the table for the record only where its rows state the
/// assembly's bases and writes nothing. The class takes a column the manifest's run hands
/// the plugin, and a field under the plugin's own prefix beside any manifest naming GRCh37,
/// on a record spelt MT or chrMT.
#[test]
fn the_dbnsfp_named_m_class_accepts_a_dbnsfp_field_the_reference_alone_fills_on_grch37() {
    let klass = DBNSFP_NAMED_M;
    for assembly in ["GRCh37", "grch37", "hg19", "37"] {
        let manifest = dbnsfp_manifest(assembly);
        for location in ["MT:3460", "chrMT:3460"] {
            for field in ["SIFT_score", "Polyphen2_HDIV_score"] {
                let (report, documented) =
                    compare_dbnsfp_rows(Some(&manifest), &format!("{field}=0.01"), "", location);
                let what = format!("{assembly} {location} {field}");
                assert_eq!(verdict(&report), "green", "{what}");
                assert_eq!(applied(&documented), once(klass), "{what}");
                assert!(masked_once(&report, field, klass), "{what}");
            }
        }
    }
    let grch37 = dbnsfp_manifest("GRCh37");
    let bare = mitochondrial_row("");
    let scored = mitochondrial_row("SIFT_score=0.01;Polyphen2_HDIV_score=0.9");
    let mut documented = from_manifest(&json(&grch37), &[klass]);
    let report = compare_texts(&scored, &bare, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    assert_eq!(verdict_of(&scored, &bare, &[]), "red");
    // The prefix is the plugin's own (dbNSFP.pm:159, :257-261), so a field under it needs no
    // word of the manifest.
    let (report, documented) = compare_dbnsfp_rows(
        Some("{\"assembly\": \"GRCh37\"}"),
        "dbNSFP_SIFT_score=0.01",
        "",
        "MT:3460",
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    // The VCF line's CHROM and the JSON object's seq_region_name spell the record, and the
    // JSON format writes the column's name in lower case.
    let inputs = vec!["MT\t3460\trs199476118\tG\tA\t.\t.\t.".to_string()];
    let vcf = vcf_text(
        "Allele|Consequence|Feature|Feature_type|SIFT_score",
        &inputs[0],
        "A|missense_variant|ENST00000361390|Transcript|",
    );
    let vcf_scored = vcf.replace("Transcript|\n", "Transcript|0.01\n");
    let mut documented = from_manifest(&json(&grch37), &[klass]);
    let report = compare(
        &vcf_entries(&vcf_scored, &inputs),
        &vcf_entries(&vcf, &inputs),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    let object = |plugin: &str| {
        format!(
            "{{\"input\": \"MT\\t3460\\trs199476118\\tG\\tA\\t.\\t.\\t.\", \"id\": \"rs199476118\", \"seq_region_name\": \"MT\", \"most_severe_consequence\": \"missense_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"A\", \"transcript_id\": \"ENST00000361390\", \"consequence_terms\": [\"missense_variant\"]{plugin}}}]}}\n"
        )
    };
    let mut documented = from_manifest(&json(&grch37), &[klass]);
    let report = compare(
        &json_entries(&object(", \"sift_score\": \"0.01\""), &inputs),
        &json_entries(&object(""), &inputs),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
}

/// The same field stays a difference beside a manifest naming GRCh38, beside one naming no
/// assembly and beside none; on a record spelt chrM or M, which dbNSFP.pm:313 leaves as
/// written, and on another chromosome; where both outputs fill it or the candidate alone
/// does; and where it is another plugin's or a column the manifest's run does not select.
#[test]
fn the_dbnsfp_named_m_class_refuses_another_assembly_contig_direction_or_plugin() {
    let grch37 = dbnsfp_manifest("GRCh37");
    let grch38 = dbnsfp_manifest("GRCh38");
    let unnamed = "{\"flags\": [\"--plugin\", \"dbNSFP,{corpus}/plugin_data/dbNSFP4.9a_grch37_synthetic.gz,SIFT_score\"]}";
    let scored = "SIFT_score=0.01";
    for (what, manifest, reference, candidate, location) in [
        ("GRCh38", Some(grch38.as_str()), scored, "", "MT:3460"),
        (
            "GRCh38 chrMT",
            Some(grch38.as_str()),
            scored,
            "",
            "chrMT:3460",
        ),
        ("no assembly", Some(unnamed), scored, "", "MT:3460"),
        ("no manifest", None, scored, "", "MT:3460"),
        (
            "no manifest, the prefix",
            None,
            "dbNSFP_SIFT_score=0.01",
            "",
            "MT:3460",
        ),
        ("chrM", Some(grch37.as_str()), scored, "", "chrM:3460"),
        ("M", Some(grch37.as_str()), scored, "", "M:3460"),
        ("21", Some(grch37.as_str()), scored, "", "21:3460"),
        (
            "both fill",
            Some(grch37.as_str()),
            scored,
            "SIFT_score=0.02",
            "MT:3460",
        ),
        (
            "the candidate alone",
            Some(grch37.as_str()),
            "",
            scored,
            "MT:3460",
        ),
        (
            "another plugin",
            Some(grch37.as_str()),
            "ClinPred=0.5",
            "",
            "MT:3460",
        ),
        (
            "a column the run does not select",
            Some(grch37.as_str()),
            "MutationTaster_score=1",
            "",
            "MT:3460",
        ),
    ] {
        let (report, documented) = compare_dbnsfp_rows(manifest, reference, candidate, location);
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&documented), nothing(), "{what}");
    }
}

/// The columns the class reads are those the words of the pair's run hand dbNSFP: the named
/// run's alone, and without a named run the manifest's own beside those of each flag set
/// over the compared input. The table and a word holding `=` are no column, a run that names
/// ALL states none, and the words of another plugin are not read. A run whose output_prefix
/// Perl reads as true writes its columns under the prefix alone (dbNSFP.pm:257-261) and
/// states none.
#[test]
fn the_dbnsfp_named_m_class_reads_the_columns_the_manifests_run_hands_the_plugin() {
    let table = "{corpus}/plugin_data/dbNSFP4.9a_grch37_synthetic.gz";
    let logic = "{corpus}/plugin_data/dbNSFP_replacement_logic";
    let manifest = json(&format!(
        r#"{{"assembly": "GRCh37",
            "flags": ["--dir_plugins", "{{plugins}}", "--plugin", "dbNSFP,{table},SIFT_score"],
            "flag_sets": [
                {{"name": "joined",
                  "flags": ["--plugin=dbNSFP,{table},output_prefix=0,pep_match=0,GERP++_RS"]}},
                {{"name": "prefixed",
                  "flags": ["--plugin", "dbNSFP,{table},output_prefix=1,Polyphen2_HDIV_score"]}},
                {{"name": "every-column",
                  "flags": ["-plugin", "dbNSFP,consequence=ALL,{table},ALL"]}},
                {{"name": "other-input", "input": "other.vcf",
                  "flags": ["--plugin", "dbNSFP,{table},{logic},Interpro_domain,SIFT_score"]}},
                {{"name": "another-plugin",
                  "flags": ["--plugin", "ClinPred,file={{corpus}}/clinpred_synthetic.tsv.gz"]}},
                {{"name": "no-flags"}}
            ]}}"#
    ));
    let read = |run: Option<&str>, input: Option<&str>| {
        Documented::from_manifest(Some(&manifest), &masks(&[DBNSFP_NAMED_M]), input, run, None)
            .expect("the manifest reads")
    };
    let none = Vec::<String>::new();
    assert_eq!(read(Some("default"), None).dbnsfp_columns, ["SIFT_score"]);
    assert_eq!(read(Some("joined"), None).dbnsfp_columns, ["GERP++_RS"]);
    assert_eq!(read(Some("prefixed"), None).dbnsfp_columns, none);
    assert_eq!(read(Some("every-column"), None).dbnsfp_columns, none);
    assert_eq!(
        read(Some("other-input"), Some("other.vcf")).dbnsfp_columns,
        [logic, "Interpro_domain", "SIFT_score"]
    );
    assert_eq!(read(Some("another-plugin"), None).dbnsfp_columns, none);
    assert_eq!(read(Some("no-flags"), None).dbnsfp_columns, none);
    assert_eq!(
        read(None, Some("variants.vcf")).dbnsfp_columns,
        ["SIFT_score", "GERP++_RS"]
    );
    assert_eq!(
        read(None, Some("other.vcf")).dbnsfp_columns,
        [logic, "Interpro_domain", "SIFT_score"]
    );
    assert_eq!(
        read(None, None).dbnsfp_columns,
        ["SIFT_score", "GERP++_RS", logic, "Interpro_domain"]
    );
    assert_eq!(documented_with(&[DBNSFP_NAMED_M]).dbnsfp_columns, none);
    // The last output_prefix of a word decides, the text between its first `=` and the next
    // (dbNSFP.pm:589, :598-599); an empty one and `0` are false.
    for (words, columns) in [
        ("output_prefix=,Made_up", &["Made_up"][..]),
        ("output_prefix=1,output_prefix=0,Made_up", &["Made_up"][..]),
        ("output_prefix==1,Made_up", &["Made_up"][..]),
        ("output_prefix=0,Made_up,output_prefix=00", &[][..]),
        ("output_prefix=1=0,Made_up", &[][..]),
    ] {
        let flags = json(&format!(r#"["--plugin", "dbNSFP,{table},{words}"]"#));
        assert_eq!(
            vep_formats_compare::manifest::dbnsfp_columns(&[&flags]),
            columns,
            "{words}"
        );
    }
    // The pair of the run named `joined` carries that run's column and not the own run's.
    let scored = |fields: &str, run: &str| {
        let mut documented = read(Some(run), None);
        let report = compare_texts(
            &mitochondrial_row(fields),
            &mitochondrial_row(""),
            &mut documented,
        );
        (verdict(&report).to_string(), applied(&documented))
    };
    assert_eq!(
        scored("GERP++_RS=5.1", "joined"),
        green_once(DBNSFP_NAMED_M)
    );
    assert_eq!(
        scored("dbNSFP_GERP++_RS=5.1", "joined"),
        green_once(DBNSFP_NAMED_M)
    );
    assert_eq!(
        scored("dbNSFP_Polyphen2_HDIV_score=0.9", "prefixed"),
        green_once(DBNSFP_NAMED_M)
    );
    assert_eq!(
        scored("Polyphen2_HDIV_score=0.9", "prefixed"),
        ("red".to_string(), nothing())
    );
    assert_eq!(
        scored("GERP++_RS=5.1", "default"),
        ("red".to_string(), nothing())
    );
    assert_eq!(
        scored("SIFT_score=0.01", "default"),
        green_once(DBNSFP_NAMED_M)
    );
    assert_eq!(
        scored("SIFT_score=0.01", "every-column"),
        ("red".to_string(), nothing())
    );
}

/// The dispatch asked directly reads the assembly and the columns it is handed: a column is
/// the class's field as the run's word spells it and in lower case, as the JSON format
/// writes it, and in no other case, which is another field's name; without an assembly, or
/// under another than GRCh37, the class accepts nothing.
#[test]
fn the_dbnsfp_named_m_class_asked_directly_reads_the_assembly_and_the_columns_handed() {
    let klass = DBNSFP_NAMED_M;
    let columns = vec!["SIFT_score".to_string(), "GERP++_RS".to_string()];
    let row = |location: &str| first_entry(&mitochondrial_row_at("", location));
    let asked = |location: &str, name: &str, assembly: Option<&str>, columns: &[String]| {
        let beside = Beside {
            assembly,
            columns,
            ..Beside::default()
        };
        accepts_beside(klass, &row(location), name, "0.01", ABSENT, &beside)
    };
    for assembly in ["GRCh37", "grch37", "hg19", "37", " GRCh37\n"] {
        for location in ["MT:3460", "chrMT:3460"] {
            for name in ["SIFT_score", "sift_score", "GERP++_RS", "dbNSFP_SIFT_score"] {
                assert!(
                    asked(location, name, Some(assembly), &columns),
                    "{assembly} {location} {name}"
                );
            }
        }
        for location in ["chrM:3460", "M:3460", "21:3460", "mt:3460"] {
            assert!(
                !asked(location, "SIFT_score", Some(assembly), &columns),
                "{assembly} {location}"
            );
        }
        assert!(!asked("MT:3460", "SIFT_score", Some(assembly), &[]));
        assert!(!asked("MT:3460", "SIFT_score_x", Some(assembly), &columns));
        for name in ["SIFT_SCORE", "Sift_score", "gerp++_RS"] {
            assert!(
                !asked("MT:3460", name, Some(assembly), &columns),
                "{assembly} {name}"
            );
        }
        assert!(!asked("MT:3460", "ClinPred", Some(assembly), &columns));
        assert!(asked("MT:3460", "dbNSFP_SIFT_score", Some(assembly), &[]));
    }
    for assembly in [None, Some("GRCh38"), Some("TAIR10"), Some("grch37.p13")] {
        for name in ["SIFT_score", "dbNSFP_SIFT_score"] {
            assert!(
                !asked("MT:3460", name, assembly, &columns),
                "{assembly:?} {name}"
            );
        }
    }
    let beside = Beside {
        assembly: Some("GRCh37"),
        columns: &columns,
        ..Beside::default()
    };
    let e = row("MT:3460");
    for (expected, actual) in [("0.01", "0.02"), (ABSENT, "0.01"), ("", ""), ("-", ABSENT)] {
        assert!(
            !accepts_beside(klass, &e, "SIFT_score", expected, actual, &beside),
            "{expected} {actual}"
        );
    }
    for actual in [ABSENT, "", "-"] {
        assert!(accepts_beside(
            klass,
            &e,
            "SIFT_score",
            "0.01",
            actual,
            &beside
        ));
    }
}

/// The record 21 26965499 GC G,GA, whose alleles share their first base: the reference's
/// plugin rebuilds the record and scores neither allele, the candidate scores the
/// substitution's rows. The class reads the record's alleles where the Uploaded_variation,
/// the VCF line or the JSON allele string spells them; under an identifier the default
/// format's rows may be two records', and it leaves them.
#[test]
fn the_chopped_multiallelic_class_accepts_the_substitution_scored_beside_a_deletion() {
    let klass = "spliceai_chopped_multiallelic";
    let spelt = "21_26965500_C/-/A";
    let reference = chopped_rows(spelt, &[("-", ""), ("A", "")]);
    let candidate = chopped_rows(spelt, &[("-", ""), ("A", SPLICEAI_VALUE)]);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    assert!(masked_once(&report, "SpliceAI_pred", klass));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let named = "chopped_multiallelic";
    let (report, documented) = compare_default(
        &chopped_rows(named, &[("-", ""), ("A", "")]),
        &chopped_rows(named, &[("-", ""), ("A", SPLICEAI_VALUE)]),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
    let inputs = vec!["21\t26965499\tchopped_multiallelic\tGC\tG,GA\t.\t.\t.".to_string()];
    let vcf = vcf_text(
        "Allele|Consequence|Feature|Feature_type|SpliceAI_pred_DS_AG",
        &inputs[0],
        "-|intron_variant|ENST00000284987|Transcript|,A|intron_variant|ENST00000284987|Transcript|",
    );
    let vcf_scored = vcf.replace("Transcript|\n", "Transcript|0.10\n");
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&vcf, &inputs),
            &vcf_entries(&vcf_scored, &inputs)
        ),
        green_once(klass)
    );
    let record = |spliceai: &str| {
        format!(
            "{{\"input\": \"21\\t26965499\\tchopped_multiallelic\\tGC\\tG,GA\\t.\\t.\\t.\", \"id\": \"chopped_multiallelic\", \"allele_string\": \"C/-/A\", \"most_severe_consequence\": \"intron_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"-\", \"transcript_id\": \"ENST00000284987\", \"consequence_terms\": [\"intron_variant\"]}}, {{\"variant_allele\": \"A\", \"transcript_id\": \"ENST00000284987\", \"consequence_terms\": [\"intron_variant\"]{spliceai}}}]}}\n"
        )
    };
    assert_eq!(
        under_class(
            klass,
            &json_entries(&record(""), &inputs),
            &json_entries(
                &record(", \"spliceai\": {\"DS_AG\": 0.1, \"SYMBOL\": \"GENE1\"}"),
                &inputs
            )
        ),
        green_once(klass)
    );
}

/// A SpliceAI field the candidate alone fills stays a difference on a record whose rows
/// carry no allele written `-`, on the deletion's own row and on an insertion's row, its
/// alleles spelt or under an identifier; one the reference alone fills or both fill
/// differently stays one beside the deletion too.
#[test]
fn the_chopped_multiallelic_class_refuses_a_record_without_a_deletion_or_another_direction() {
    let klass = "spliceai_chopped_multiallelic";
    let name = "21_26965500_C/-/A";
    let other_value = SPLICEAI_VALUE.replace("0.10", "0.11");
    let two = "21_26965500_C/G/A";
    let one = "21_26965500_C/A";
    let inserted = "21_26965500_C/-/CA";
    let named = "chopped_multiallelic";
    for (reference, candidate) in [
        (
            chopped_rows(two, &[("G", ""), ("A", "")]),
            chopped_rows(two, &[("G", ""), ("A", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows(one, &[("A", "")]),
            chopped_rows(one, &[("A", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows(name, &[("-", ""), ("A", "")]),
            chopped_rows(name, &[("-", SPLICEAI_VALUE), ("A", "")]),
        ),
        (
            chopped_rows(inserted, &[("-", ""), ("CA", "")]),
            chopped_rows(inserted, &[("-", ""), ("CA", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows(named, &[("-", ""), ("CA", "")]),
            chopped_rows(named, &[("-", ""), ("CA", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows(name, &[("-", ""), ("A", SPLICEAI_VALUE)]),
            chopped_rows(name, &[("-", ""), ("A", "")]),
        ),
        (
            chopped_rows(name, &[("-", ""), ("A", SPLICEAI_VALUE)]),
            chopped_rows(name, &[("-", ""), ("A", &other_value)]),
        ),
    ] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
}

/// On a model whose exons carry no stable id: a deletion across an intron prints the first
/// exon's distance and side with the second exon's length in the reference and the first
/// exon's own item in the candidate; under `intronic=1` the reference prints one item for
/// the two flanking exons, either exon's distance and side, where the candidate prints both
/// items. The VCF format joins the parts with `+` and the items with `&`.
#[test]
fn the_id_less_exon_fold_class_accepts_one_item_folded_from_two_id_less_exons() {
    let klass = "nearest_exon_jb_id_less_exon_fold";
    let reference = two_exon_row("|64|start|89");
    let candidate = two_exon_row("|64|start|82");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    assert!(masked_once(&report, "NearestExonJB", klass));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let flanking = nearest_exon_row(IDLESS_INTRONIC, "|106|end|135,|2780|start|109");
    for folded in ["|106|end|109", "|2780|start|109"] {
        let (report, documented) = compare_default(
            &nearest_exon_row(IDLESS_INTRONIC, folded),
            &flanking,
            &[klass],
        );
        assert_eq!(verdict(&report), "green", "{folded}");
        assert_eq!(applied(&documented), once(klass), "{folded}");
    }
    let inputs = vec!["21\t10521800\tidless_intronic\tT\tA\t.\t.\t.".to_string()];
    let vcf = vcf_text(
        "Allele|Consequence|Feature|Feature_type|NearestExonJB",
        &inputs[0],
        "A|intron_variant|compmerge.2.Heart.chr21|Transcript|+106+end+109",
    );
    let vcf_candidate = vcf.replace("+106+end+109", "+106+end+135&+2780+start+109");
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&vcf, &inputs),
            &vcf_entries(&vcf_candidate, &inputs)
        ),
        green_once(klass)
    );
}

/// A NearestExonJB stays a difference where an item names an exon, where the reference
/// prints two items, where its distance or side is none of the candidate's, where either
/// side leaves the field empty, and where one candidate item differs in its length on a row
/// that does not begin in one exon and end in another.
#[test]
fn the_id_less_exon_fold_class_refuses_an_exon_with_an_id_or_another_distance() {
    let klass = "nearest_exon_jb_id_less_exon_fold";
    let one_exon: Row = [
        IDLESS_TWO_EXONS[0],
        IDLESS_TWO_EXONS[1],
        IDLESS_TWO_EXONS[2],
        IDLESS_TWO_EXONS[3],
        "non_coding_transcript_exon_variant",
        IDLESS_TWO_EXONS[5],
    ];
    let (report, documented) = compare_default(
        &nearest_exon_row_at(one_exon, "|64|start|89", "1266-1297"),
        &nearest_exon_row_at(one_exon, "|64|start|82", "1266-1297"),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), nothing());
    for cdna in ["-", "?-1297", "1266-?"] {
        let (report, documented) = compare_default(
            &two_exon_row_at("|64|start|89", cdna),
            &two_exon_row_at("|64|start|82", cdna),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{cdna}");
        assert_eq!(applied(&documented), nothing(), "{cdna}");
    }
    for (reference, candidate) in [
        ("ENSE00001659869|64|start|89", "ENSE00001659869|64|start|82"),
        ("|64|start|89", "ENSE00001659869|64|start|82"),
        ("|106|end|135,|2780|start|109", "|106|end|109"),
        ("|65|start|89", "|64|start|82"),
        ("|64|end|89", "|64|start|82"),
        ("|64|start|89", ""),
        ("", "|64|start|82"),
    ] {
        let (report, documented) =
            compare_default(&two_exon_row(reference), &two_exon_row(candidate), &[klass]);
        assert_eq!(verdict(&report), "red", "{reference} {candidate}");
        assert_eq!(applied(&documented), nothing(), "{reference} {candidate}");
    }
}

/// 21 30166352 G>GAG on a reverse-strand transcript shifts seven bases to `c.271_272dup`:
/// the reference splices the two bases unturned, the candidate turns them one place; a
/// ProteinLengthChange differing on such a row is the class's too, as is a deletion's row
/// and the JSON format's keys.
#[test]
fn the_reverse_strand_rotation_class_accepts_the_protein_of_an_allele_shifted_past_its_length() {
    let klass = "downstream_reverse_strand_rotation";
    let (report, documented) =
        compare_default(&downstream_row("SS"), &downstream_row("SP"), &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    assert!(masked_once(&report, "DownstreamProtein", klass));
    assert_eq!(
        verdict_of(&downstream_row("SS"), &downstream_row("SP"), &[]),
        "red"
    );
    let shorter = DownstreamRow {
        length_change: "-131",
        ..DownstreamRow::default()
    };
    let (report, documented) =
        compare_default(&downstream_row("SS"), &shorter.text("SP"), &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    let deletion = DownstreamRow {
        offset: "-4",
        notation: "c.100_102del",
        ..DownstreamRow::default()
    };
    assert_eq!(
        verdict_of(&deletion.text("SS"), &deletion.text("SP"), &[klass]),
        "green"
    );
    let inputs = vec!["21\t30166352\tcldn17_ins_AG\tG\tGAG\t.\t.\t.".to_string()];
    let record = |protein: &str| {
        format!(
            "{{\"input\": \"21\\t30166352\\tcldn17_ins_AG\\tG\\tGAG\\t.\\t.\\t.\", \"id\": \"cldn17_ins_AG\", \"most_severe_consequence\": \"frameshift_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"AG\", \"transcript_id\": \"ENST00000286808\", \"consequence_terms\": [\"frameshift_variant\"], \"strand\": -1, \"hgvsc\": \"ENST00000286808:c.271_272dup\", \"hgvs_offset\": -7, \"downstreamprotein\": \"{protein}\", \"proteinlengthchange\": -132}}]}}\n"
        )
    };
    assert_eq!(
        under_class(
            klass,
            &json_entries(&record("SS"), &inputs),
            &json_entries(&record("SP"), &inputs)
        ),
        green_once(klass)
    );
}

/// A differing DownstreamProtein stays a difference on a forward-strand row, where the shift
/// does not exceed the allele or is a multiple of its length, without an HGVS_OFFSET or an
/// indel notation, and on a row that is not a frameshift's.
#[test]
fn the_reverse_strand_rotation_class_refuses_a_shift_the_rotation_reads_alike() {
    let klass = "downstream_reverse_strand_rotation";
    let usual = DownstreamRow::default();
    for change in [
        DownstreamRow {
            strand: "1",
            ..usual
        },
        DownstreamRow {
            offset: "-2",
            ..usual
        },
        DownstreamRow {
            offset: "-1",
            ..usual
        },
        DownstreamRow {
            offset: "-6",
            ..usual
        },
        DownstreamRow {
            offset: "-",
            ..usual
        },
        DownstreamRow {
            notation: "c.271A>G",
            ..usual
        },
        DownstreamRow {
            notation: "c.10_5del",
            ..usual
        },
        DownstreamRow {
            consequence: "inframe_insertion",
            ..usual
        },
    ] {
        let (reference, candidate) = (change.text("SS"), change.text("SP"));
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
    // A value one engine lacks stays a difference: the class is two values that differ.
    for (reference, candidate) in [
        (downstream_row("SS"), downstream_row("-")),
        (downstream_row("-"), downstream_row("SP")),
    ] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red");
        assert_eq!(applied(&documented), nothing());
    }
}

/// A frameshift 118 bases walked and 12 from the last exon: the reference's GERP total, read
/// through a zoom level, lies below the cutoff and carries END_TRUNC with LC, the candidate's
/// per-base sum lies above it with HC and no filter. The row is counted once, its LoF_info
/// with its LoF_filter and LoF; two totals on one side of the cutoff differ in LoF_info
/// alone; the lower total may be the candidate's.
#[test]
fn the_gerp_summary_class_accepts_a_total_read_through_a_zoom_level() {
    let klass = GERP;
    let filtered = gerp_row(GERP_FILTERED, GERP_REFERENCE);
    let clear = gerp_row(GERP_CLEAR, GERP_CANDIDATE);
    let (report, documented) = compare_default(&filtered, &clear, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    for field in ["LoF", "LoF_filter", "LoF_info"] {
        assert!(masked_once(&report, field, klass), "{field}");
    }
    assert_eq!(verdict_of(&filtered, &clear, &[]), "red");
    let (report, documented) = compare_default(&clear, &filtered, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    let (report, _) = compare_default(
        &gerp_row(GERP_FILTERED, "-100.5"),
        &gerp_row(GERP_FILTERED, "-90.25"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert!(masked_once(&report, "LoF_info", klass));
    // A stop gain 108 bases from the stop codon and inside the last exon's reach: the
    // reference's total carries the filter, the candidate's does not.
    let stop_gain = GerpRow {
        walked: "108",
        reach: "-16",
        consequence: "stop_gained",
    };
    let (report, documented) = compare_default(
        &stop_gain.text(GERP_FILTERED, "-64.1545345306396"),
        &stop_gain.text(GERP_CLEAR, "-34.7560015916824"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    // An interval of 80 bases, the shortest a zoom level serves, walks 79.
    let shortest = GerpRow {
        walked: "79",
        ..GerpRow::default()
    };
    let (report, documented) = compare_default(
        &shortest.text(GERP_FILTERED, GERP_REFERENCE),
        &shortest.text(GERP_CLEAR, GERP_CANDIDATE),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    let inputs = vec!["21\t36759739\t2780984\tCA\tC\t.\t.\t.".to_string()];
    let line = |lof: &str, lof_filter: &str, gerp: &str| {
        vcf_text(
            "Allele|Consequence|Feature|Feature_type|LoF|LoF_filter|LoF_info",
            &inputs[0],
            &format!(
                "-|frameshift_variant|ENST00001017196|Transcript|{lof}|{lof_filter}|PERCENTILE:0.95&GERP_DIST:{gerp}&BP_DIST:118&DIST_FROM_LAST_EXON:12&50_BP_RULE:FAIL"
            ),
        )
    };
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&line("LC", "END_TRUNC", GERP_REFERENCE), &inputs),
            &vcf_entries(&line("HC", "", GERP_CANDIDATE), &inputs)
        ),
        green_once(klass)
    );
}

/// A GERP total stays a difference on an interval under 80 bases, beside another differing
/// item, on a row that is neither a stop gain nor a frameshift, and where the two totals are
/// one number spelt twice; the filter stays one where the totals do not straddle the cutoff,
/// where it stands beside the higher total, and farther than 50 bases from the last exon.
#[test]
fn the_gerp_summary_class_refuses_a_short_walk_another_item_or_a_filter_off_the_cutoff() {
    let klass = GERP;
    let short = GerpRow {
        walked: "78",
        ..GerpRow::default()
    };
    let longer = GerpRow {
        walked: "119",
        ..GerpRow::default()
    };
    let missense = GerpRow {
        consequence: "missense_variant",
        ..GerpRow::default()
    };
    for (reference, candidate) in [
        (
            short.text(GERP_FILTERED, GERP_REFERENCE),
            short.text(GERP_CLEAR, GERP_CANDIDATE),
        ),
        (
            gerp_row(GERP_FILTERED, GERP_REFERENCE),
            longer.text(GERP_CLEAR, GERP_CANDIDATE),
        ),
        (
            missense.text(GERP_FILTERED, GERP_REFERENCE),
            missense.text(GERP_CLEAR, GERP_CANDIDATE),
        ),
        (
            gerp_row(GERP_FILTERED, "-69.50"),
            gerp_row(GERP_CLEAR, "-69.5"),
        ),
    ] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
    // The total itself is the class's on these rows; the filter beside it is not.
    let farther = GerpRow {
        reach: "51",
        ..GerpRow::default()
    };
    for (reference, candidate) in [
        (
            gerp_row(GERP_FILTERED, "-69.5"),
            gerp_row(GERP_CLEAR, "-60.5"),
        ),
        (
            gerp_row(GERP_CLEAR, GERP_REFERENCE),
            gerp_row(GERP_FILTERED, GERP_CANDIDATE),
        ),
        (
            farther.text(GERP_FILTERED, GERP_REFERENCE),
            farther.text(GERP_CLEAR, GERP_CANDIDATE),
        ),
        (
            gerp_row(GERP_FILTERED, GERP_REFERENCE),
            gerp_row("LoF=LC;LoF_filter=ANC_ALLELE", GERP_CANDIDATE),
        ),
    ] {
        let (report, _) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(differing(&report, "LoF_filter"), 1, "{candidate}");
        assert!(masked_once(&report, "LoF_info", klass), "{candidate}");
    }
}

const WHOLE_EXON_WALK: &str = "loftee_gerp_walk_without_conservation_database";
const WHOLE_EXON_WALK_ITEMS: &str =
    "GERP_DIST:705,BP_DIST:702,PERCENTILE:0.155875299760192,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS";
const WHOLE_EXON_WALK_INPUT: &str = "21\t37442543\tfirst_exon_stop\tC\tT\t.\t.\t.";

/// The default-format row of the stop gain at 21:37442543 on ENST00000290349, whose GERP walk
/// crosses the whole second exon of three: the terms `consequence` and the pairs `lof` after
/// the row's IMPACT and STRAND.
fn whole_exon_walk_row(consequence: &str, lof: &str) -> String {
    let extra = ["IMPACT=HIGH;STRAND=1", lof]
        .into_iter()
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>()
        .join(";");
    default_text(&[default_row(
        "first_exon_stop",
        "21:37442543",
        "T",
        "ENST00000290349",
        consequence,
        &extra,
    )])
}

/// The same row in the tab format, a column for each of the plugin's four fields.
fn whole_exon_walk_tab_row(lof: [&str; 4]) -> Vec<Entry> {
    let header =
        DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tLoF\tLoF_filter\tLoF_flags\tLoF_info");
    let row = default_row(
        "first_exon_stop",
        "21:37442543",
        "T",
        "ENST00000290349",
        "stop_gained",
        "HIGH",
    );
    parse_tab(&format!("{header}\n{row}\t{}\n", lof.join("\t"))).expect("the tab text parses")
}

/// A stop gain in the first of three exons, in a run without a conservation database: the
/// reference's plugin dies in the walk and the row carries none of its fields, the candidate
/// writes the call and its items. The row is counted once with every field the candidate
/// fills, a filter or a flag among them, on a frameshift row and beside another term too,
/// and in the tab, VCF and JSON formats.
#[test]
fn the_whole_exon_walk_class_accepts_the_fields_the_candidate_alone_writes() {
    let klass = WHOLE_EXON_WALK;
    let items = format!("LoF_info={WHOLE_EXON_WALK_ITEMS}");
    let unwritten = whole_exon_walk_row("stop_gained", "");
    let written = whole_exon_walk_row("stop_gained", &format!("LoF=HC;{items}"));
    let (report, documented) = compare_default(&unwritten, &written, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), once(klass));
    for field in ["LoF", "LoF_info", "Extra"] {
        assert!(masked_once(&report, field, klass), "{field}");
    }
    assert_eq!(verdict_of(&unwritten, &written, &[]), "red");
    let three_items = "LoF_info=GERP_DIST:1498,BP_DIST:1496,PERCENTILE:0.240412979351032";
    for (consequence, lof, third) in [
        (
            "stop_gained",
            format!("LoF=LC;LoF_filter=END_TRUNC;{items}"),
            "LoF_filter",
        ),
        (
            "frameshift_variant,splice_region_variant",
            format!("LoF=HC;LoF_flags=NO_EXON_NUMBER;{three_items}"),
            "LoF_flags",
        ),
        ("frameshift_variant", format!("LoF=HC;{items}"), "Extra"),
    ] {
        let reference = whole_exon_walk_row(consequence, "");
        let candidate = whole_exon_walk_row(consequence, &lof);
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green", "{lof}");
        assert_eq!(applied(&documented), once(klass), "{lof}");
        for field in ["LoF", "LoF_info", third] {
            assert!(masked_once(&report, field, klass), "{lof} {field}");
        }
    }
    assert_eq!(
        under_class(
            klass,
            &whole_exon_walk_tab_row(["-", "-", "-", "-"]),
            &whole_exon_walk_tab_row(["LC", "END_TRUNC", "-", WHOLE_EXON_WALK_ITEMS])
        ),
        green_once(klass)
    );
    let inputs = vec![WHOLE_EXON_WALK_INPUT.to_string()];
    let line = |lof: &str| {
        vcf_text(
            "Allele|Consequence|Feature|Feature_type|LoF|LoF_filter|LoF_flags|LoF_info",
            &inputs[0],
            &format!("T|stop_gained|ENST00000290349|Transcript|{lof}"),
        )
    };
    let joined_by_ampersands = WHOLE_EXON_WALK_ITEMS.replace(',', "&");
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&line("|||"), &inputs),
            &vcf_entries(&line(&format!("HC|||{joined_by_ampersands}")), &inputs)
        ),
        green_once(klass)
    );
    let record = |lof: &str| {
        format!(
            "{{\"input\": \"21\\t37442543\\tfirst_exon_stop\\tC\\tT\\t.\\t.\\t.\", \"id\": \"first_exon_stop\", \"most_severe_consequence\": \"stop_gained\", \"transcript_consequences\": [{{\"variant_allele\": \"T\", \"transcript_id\": \"ENST00000290349\", \"consequence_terms\": [\"stop_gained\"]{lof}}}]}}\n"
        )
    };
    let json_items = format!(
        ", \"lof\": \"LC\", \"lof_filter\": \"END_TRUNC\", \"lof_info\": \"{WHOLE_EXON_WALK_ITEMS}\""
    );
    assert_eq!(
        under_class(
            klass,
            &json_entries(&record(""), &inputs),
            &json_entries(&record(&json_items), &inputs)
        ),
        green_once(klass)
    );
}

/// Under `skip_lof_info` no row carries LoF_info: the candidate's row holds the call, with
/// its filter or its flag where one applies, and the class takes it in every format.
#[test]
fn the_whole_exon_walk_class_accepts_a_call_written_without_its_items() {
    let klass = WHOLE_EXON_WALK;
    for (lof, second) in [
        ("LoF=HC", "Extra"),
        ("LoF=LC;LoF_filter=END_TRUNC", "LoF_filter"),
        ("LoF=HC;LoF_flags=NO_EXON_NUMBER", "LoF_flags"),
    ] {
        let reference = whole_exon_walk_row("stop_gained", "");
        let candidate = whole_exon_walk_row("stop_gained", lof);
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green", "{lof}");
        assert_eq!(applied(&documented), once(klass), "{lof}");
        for field in ["LoF", second] {
            assert!(masked_once(&report, field, klass), "{lof} {field}");
        }
        assert_eq!(verdict_of(&reference, &candidate, &[]), "red", "{lof}");
    }
    assert_eq!(
        under_class(
            klass,
            &whole_exon_walk_tab_row(["-", "-", "-", "-"]),
            &whole_exon_walk_tab_row(["HC", "-", "NO_EXON_NUMBER", "-"])
        ),
        green_once(klass)
    );
    let inputs = vec![WHOLE_EXON_WALK_INPUT.to_string()];
    let line = |lof: &str| {
        vcf_text(
            "Allele|Consequence|Feature|Feature_type|LoF|LoF_filter|LoF_flags|LoF_info",
            &inputs[0],
            &format!("T|frameshift_variant&splice_region_variant|ENST00000290349|Transcript|{lof}"),
        )
    };
    assert_eq!(
        under_class(
            klass,
            &vcf_entries(&line("|||"), &inputs),
            &vcf_entries(&line("LC|END_TRUNC||"), &inputs)
        ),
        green_once(klass)
    );
}

/// The fields stay differences on a row where the reference writes one of the four, on a
/// row that is neither a stop gain nor a frameshift, where the candidate writes no call,
/// where its items hold no GERP total, and where the reference alone writes the fields; a
/// field outside the four stays one beside them.
#[test]
fn the_whole_exon_walk_class_refuses_a_row_the_reference_annotated_or_another_term() {
    let klass = WHOLE_EXON_WALK;
    let items = format!("LoF_info={WHOLE_EXON_WALK_ITEMS}");
    let written = format!("LoF=HC;{items}");
    for (consequence, reference, candidate) in [
        (
            "stop_gained",
            "LoF=HC",
            format!("LoF=LC;LoF_filter=END_TRUNC;{items}"),
        ),
        ("stop_gained", "LoF_flags=SINGLE_EXON", written.clone()),
        ("stop_gained", "LoF_filter=END_TRUNC", written.clone()),
        (
            "stop_gained",
            "LoF_info=PERCENTILE:0.155875299760192",
            written.clone(),
        ),
        ("splice_donor_variant", "", written.clone()),
        ("missense_variant", "", written.clone()),
        ("splice_donor_variant", "", "LoF=HC".to_string()),
        ("stop_gained", "", items.clone()),
        ("stop_gained", "", "LoF_filter=END_TRUNC".to_string()),
        (
            "stop_gained",
            "",
            "LoF=HC;LoF_info=BP_DIST:702,PERCENTILE:0.155875299760192".to_string(),
        ),
        ("frameshift_variant", written.as_str(), String::new()),
    ] {
        let reference = whole_exon_walk_row(consequence, reference);
        let candidate = whole_exon_walk_row(consequence, &candidate);
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(applied(&documented), nothing(), "{candidate}");
    }
    // Another plugin's field on a row of the class is no field of the four.
    let (report, documented) = compare_default(
        &whole_exon_walk_row("stop_gained", ""),
        &whole_exon_walk_row(
            "stop_gained",
            &format!("{written};NMD=NMD_escaping_variant"),
        ),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "NMD"), 1);
    assert!(masked_once(&report, "LoF", klass));
    assert!(at(&report, &["columns", "Extra"]).get("masked").is_none());
    assert_eq!(applied(&documented), once(klass));
    // The dispatch reads both rows: a pair in which one row holds neither term, or no
    // candidate row at all, is refused, and so is a value the caller states for the reference.
    let e = first_entry(&whole_exon_walk_row("stop_gained", ""));
    let a = first_entry(&whole_exon_walk_row("stop_gained", &written));
    let other_terms = first_entry(&whole_exon_walk_row("missense_variant", &written));
    let reference_other_terms = first_entry(&whole_exon_walk_row("missense_variant", ""));
    assert!(accepts(klass, &e, "LoF", ABSENT, "HC", Some(&a)));
    assert!(accepts(
        klass,
        &e,
        "lof_info",
        "-",
        WHOLE_EXON_WALK_ITEMS,
        Some(&a)
    ));
    assert!(!accepts(klass, &e, "LoF", ABSENT, "HC", Some(&other_terms)));
    assert!(!accepts(
        klass,
        &reference_other_terms,
        "LoF",
        ABSENT,
        "HC",
        Some(&a)
    ));
    assert!(!accepts(klass, &e, "LoF", ABSENT, "HC", None));
    assert!(!accepts(klass, &e, "LoF", "LC", "HC", Some(&a)));
    assert!(!accepts(klass, &e, "LoF_filter", ABSENT, "-", Some(&a)));
    assert!(!accepts(klass, &e, "IMPACT", ABSENT, "HIGH", Some(&a)));
}

#[test]
fn the_whole_exon_walk_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    a_planted_value_stays_a_difference(WHOLE_EXON_WALK, "stop_gained", "", "LoF", "HC", "LC");
}
