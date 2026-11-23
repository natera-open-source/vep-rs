// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the classes that read the candidate cache's transcripts share: the
//! planted cache, the rows of the Downstream and NMD plugins in the four formats and the
//! readers of a comparison's report.

#![allow(dead_code, unused_imports)]

use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;

use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::counter::Counter;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::Side;

pub use vep_formats_compare::classes::BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES as BLOSUM62_SHARED_SHIFT;
pub use vep_formats_compare::classes::DOWNSTREAM_START_PAST_OWN_CODON as START_PAST_OWN_CODON;
pub use vep_formats_compare::classes::DOWNSTREAM_UTR_FRAME_PAST_CODING_END as UTR_FRAME_PAST_CODING_END;
pub use vep_formats_compare::classes::HGVS_SHIFT_SHARED_ACROSS_ALLELES as SHARED_SHIFT;
pub use vep_formats_compare::classes::PLUGIN_HGVS_SHIFTED_CDS_COORDINATES as SHIFTED_CDS_COORDINATES;

/// The header line of a default-format output.
pub const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

/// A coding sequence of eight codons and a stop codon whose 3' UTR, read from its first
/// base, reaches a stop after three codons; and one of four codons and a stop codon whose
/// second codon, less its third base, spells a stop with the base after it.
pub const TOY_CODING: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";
pub const TOY_UTR: &str = "GGGCCCAAATAGCCC";
pub const SHORT_CODING: &str = "ATGTACAAAGGTTAA";
pub const SHORT_UTR: &str = "GTGACC";
/// A 3' UTR the vertebrate mitochondrial code reads unlike the standard one in both frames
/// the classes read: ATA, methionine there, from its first base; TGA, tryptophan there, and
/// AGA, a stop there, from its second.
pub const TABLE_2_UTR: &str = "GTGACCAAGATAGCC";
pub const TOY: &str = "ENST00000000101";
pub const TOY_NO_UTR: &str = "ENST00000000102";
pub const TOY_TABLE_2: &str = "ENST00000000103";
pub const TOY_UTR_STOP: &str = "ENST00000000104";
pub const SHORT: &str = "ENST00000000105";
pub const TOY_TABLE_7: &str = "ENST00000000106";
/// The transcripts the cache states a model for. Four carry the toy sequences: one exon on
/// the forward strand, coding base c at 21:1000+c; one exon on the reverse strand, coding
/// base c at 21:2043-c; two exons, coding bases 1 to 9 at 21:3001-3009 and the rest from
/// 21:3014; and one exon at 21:4001-4039 whose mapper pairs touch on the genome after coding
/// base 9 while the model carries bases 10 to 12 between them. One holds a run of one base
/// from a codon's third base (coding base c at 21:6000+c) and one a coding sequence of 120
/// bases over three exons, bases 1 to 110 at 21:5001-5110.
pub const TOY_FORWARD: &str = "ENST00000000201";
pub const TOY_REVERSE: &str = "ENST00000000202";
pub const TOY_TWO_EXONS: &str = "ENST00000000203";
pub const TOY_SEAM: &str = "ENST00000000204";
pub const RUN: &str = "ENST00000000205";
pub const THREE_EXONS: &str = "ENST00000000206";
pub const RUN_CODING: &str = "ATGCCAAAAGGTTAA";
pub const DOWNSTREAM_TAB_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tCDS_position\tSTRAND\tDownstreamProtein\tProteinLengthChange";
pub const DOWNSTREAM_CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|CDS_position|STRAND|DownstreamProtein|ProteinLengthChange\">\n";
pub const NMD_VALUE: &str = "NMD_escaping_variant";
pub const ACCEPTOR_TERMS: &str = "splice_acceptor_variant,coding_sequence_variant,intron_variant";

pub fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

pub fn whole(value: i64) -> Value {
    Value::Int(value)
}

pub fn object(members: Vec<(&str, Value)>) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect(),
    )
}

pub fn list(items: Vec<Value>) -> Value {
    Value::List(items)
}

/// `value` as one line of JSON.
pub fn dumped_value(value: &Value) -> String {
    pyjson::dumps(value, &pyjson::DumpOptions::default())
}

/// Writes `text` to `name` under `directory`, with the directories above it.
pub fn written(directory: &Path, name: &str, text: &str) -> PathBuf {
    let path = directory.join(name);
    std::fs::create_dir_all(path.parent().expect("a file has a directory"))
        .expect("the directory is made");
    std::fs::write(&path, text).expect("the file is written");
    path
}

/// A cache transcript without a model: its stable id and what its
/// `variation_effect_feature_cache` holds.
pub fn cached_transcript(stable_id: &str, cached: Vec<(&str, Value)>) -> Value {
    object(vec![
        ("stable_id", text(stable_id)),
        ("variation_effect_feature_cache", object(cached)),
    ])
}

/// What a transcript with a model states beside its mapper pairs.
#[derive(Clone)]
pub struct Modelled<'a> {
    pub coding: &'a str,
    pub utr: &'a str,
    /// The exons, where they are not the genomic runs of the mapper pairs.
    pub exons: &'a [(i64, i64)],
    /// Whether every number is written as a string.
    pub as_text: bool,
    pub peptide: Option<&'a str>,
    /// The `codon_table` the cache names.
    pub table: i64,
}

impl Default for Modelled<'_> {
    fn default() -> Self {
        Modelled {
            coding: TOY_CODING,
            utr: TOY_UTR,
            exons: &[],
            as_text: false,
            peptide: None,
            table: 1,
        }
    }
}

/// A cache transcript with a model: each of `pairs` a run of cDNA coordinates and the
/// genomic bases it maps to, the coding sequence from cDNA base 1.
pub fn modelled(
    stable_id: &str,
    strand: i64,
    pairs: &[(i64, i64, i64, i64)],
    shape: Modelled,
) -> Value {
    let number = |value: i64| {
        if shape.as_text {
            text(&value.to_string())
        } else {
            whole(value)
        }
    };
    let mapper_pairs = pairs
        .iter()
        .map(|(from_start, from_end, to_start, to_end)| {
            object(vec![
                ("from_start", number(*from_start)),
                ("from_end", number(*from_end)),
                ("to_start", number(*to_start)),
                ("to_end", number(*to_end)),
                ("ori", number(strand)),
            ])
        })
        .collect();
    let mut spans: Vec<(i64, i64)> = if shape.exons.is_empty() {
        pairs.iter().map(|pair| (pair.2, pair.3)).collect()
    } else {
        shape.exons.to_vec()
    };
    spans.sort_unstable();
    if strand == -1 {
        spans.reverse();
    }
    let mut cached = vec![
        ("translateable_seq", text(shape.coding)),
        ("three_prime_utr", text(shape.utr)),
        ("codon_table", whole(shape.table)),
        ("mapper", object(vec![("pairs", list(mapper_pairs))])),
    ];
    if let Some(peptide) = shape.peptide {
        cached.push(("peptide", text(peptide)));
    }
    object(vec![
        ("stable_id", text(stable_id)),
        ("strand", number(strand)),
        ("cdna_coding_start", number(1)),
        (
            "cdna_coding_end",
            number(i64::try_from(shape.coding.len()).expect("a length")),
        ),
        (
            "exons",
            list(
                spans
                    .iter()
                    .map(|(start, end)| {
                        object(vec![("start", number(*start)), ("end", number(*end))])
                    })
                    .collect(),
            ),
        ),
        ("variation_effect_feature_cache", object(cached)),
    ])
}

/// A candidate JSON cache of twelve coding transcripts. Chromosome 21 holds the toy
/// transcript (its `codon_table` the string a converted cache writes), a copy without a 3'
/// UTR (the number a native cache writes), a copy under codon table 2 with `TABLE_2_UTR`, a
/// copy whose 3' UTR opens with a stop codon (no `codon_table`) and a copy under table 7,
/// which no engine holds a table for, none with a model, and the six transcripts
/// with one (the reverse-strand one's numbers the strings a converted cache writes, the
/// seam's with its `peptide`); chromosome 22 holds the short transcript (table 0, read as
/// table 1).
pub fn downstream_cache(root: &Path) -> PathBuf {
    let utr_opening_with_a_stop = format!("TGA{}", &TOY_UTR[3..]);
    let three_exons_coding = format!("ATG{}TAA", "GCT".repeat(38));
    let chromosome_21 = list(vec![
        cached_transcript(
            TOY,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("three_prime_utr", text(TOY_UTR)),
                ("codon_table", text("1")),
            ],
        ),
        cached_transcript(
            TOY_NO_UTR,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("codon_table", whole(1)),
            ],
        ),
        cached_transcript(
            TOY_TABLE_2,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("three_prime_utr", text(TABLE_2_UTR)),
                ("codon_table", text("2")),
            ],
        ),
        cached_transcript(
            TOY_UTR_STOP,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("three_prime_utr", text(&utr_opening_with_a_stop)),
            ],
        ),
        cached_transcript(
            TOY_TABLE_7,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("three_prime_utr", text(TOY_UTR)),
                ("codon_table", whole(7)),
            ],
        ),
        modelled(TOY_FORWARD, 1, &[(1, 42, 1001, 1042)], Modelled::default()),
        modelled(
            TOY_REVERSE,
            -1,
            &[(1, 42, 2001, 2042)],
            Modelled {
                as_text: true,
                ..Modelled::default()
            },
        ),
        modelled(
            TOY_TWO_EXONS,
            1,
            &[(1, 9, 3001, 3009), (10, 42, 3014, 3046)],
            Modelled::default(),
        ),
        modelled(
            TOY_SEAM,
            1,
            &[(1, 9, 4001, 4009), (13, 42, 4010, 4039)],
            Modelled {
                exons: &[(4001, 4039)],
                peptide: Some("MKTLAVRG"),
                ..Modelled::default()
            },
        ),
        modelled(
            RUN,
            1,
            &[(1, 21, 6001, 6021)],
            Modelled {
                coding: RUN_CODING,
                utr: SHORT_UTR,
                ..Modelled::default()
            },
        ),
        modelled(
            THREE_EXONS,
            1,
            &[
                (1, 110, 5001, 5110),
                (111, 120, 5201, 5210),
                (121, 130, 5301, 5310),
            ],
            Modelled {
                coding: &three_exons_coding,
                utr: "GGGCCCAAAT",
                ..Modelled::default()
            },
        ),
    ]);
    let chromosome_22 = list(vec![cached_transcript(
        SHORT,
        vec![
            ("translateable_seq", text(SHORT_CODING)),
            ("three_prime_utr", text(SHORT_UTR)),
            ("codon_table", whole(0)),
        ],
    )]);
    written(
        root,
        "transcripts/21/1-1000000.json",
        &dumped_value(&chromosome_21),
    );
    written(
        root,
        "transcripts/22/1-1000000.json",
        &dumped_value(&chromosome_22),
    );
    root.to_path_buf()
}

/// The side tables of a comparison that reads the transcripts of `cache`.
pub fn transcripts_of(cache: &Path) -> Side {
    Side {
        downstream_sequences: Some(
            DownstreamSequences::new(cache).expect("the cache has a transcripts directory"),
        ),
        ..Side::default()
    }
}

/// The transcripts a side table of `transcripts_of` reads.
pub fn sequences_of(side: &Side) -> &DownstreamSequences {
    side.downstream_sequences
        .as_ref()
        .expect("the side tables read a cache")
}

/// One row of a default-format output: a Transcript row of the gene `ENSG1`.
pub fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
) -> String {
    format!(
        "{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{consequence}\t-\t-\t-\t-\t-\t-\t{extra}"
    )
}

/// A default-format output of `rows`.
pub fn default_text<S: AsRef<str>>(rows: &[S]) -> String {
    let mut lines = vec![HEADER];
    lines.extend(rows.iter().map(AsRef::as_ref));
    lines.join("\n") + "\n"
}

/// `row` with the cells from its cDNA_position to its Existing_variation, six `-`, replaced
/// by `cells`.
pub fn with_coding_cells(row: &str, cells: &str) -> String {
    row.replacen("\t-\t-\t-\t-\t-\t-\t", &format!("\t{cells}\t"), 1)
}

/// One row of the Downstream plugin on both engines: its Location, Allele, Feature and
/// CDS_position, the reference's DownstreamProtein and the candidate's, the
/// ProteinLengthChange, and what else the row states.
#[derive(Clone, Debug)]
pub struct DownstreamRow {
    pub location: String,
    pub allele: String,
    pub feature: String,
    pub cds: String,
    pub vep: String,
    pub vep_rs: String,
    pub change: String,
    pub consequence: String,
    pub strand: i64,
    pub offset: Option<i64>,
}

pub fn downstream_row(
    location: &str,
    allele: &str,
    feature: &str,
    cds: &str,
    vep: &str,
    vep_rs: &str,
    change: i64,
) -> DownstreamRow {
    DownstreamRow {
        location: location.to_string(),
        allele: allele.to_string(),
        feature: feature.to_string(),
        cds: cds.to_string(),
        vep: vep.to_string(),
        vep_rs: vep_rs.to_string(),
        change: change.to_string(),
        consequence: "frameshift_variant".to_string(),
        strand: 1,
        offset: None,
    }
}

impl DownstreamRow {
    pub fn at(mut self, location: &str, cds: &str) -> Self {
        self.location = location.to_string();
        self.cds = cds.to_string();
        self
    }

    pub fn allele(mut self, allele: &str) -> Self {
        self.allele = allele.to_string();
        self
    }

    pub fn feature(mut self, feature: &str) -> Self {
        self.feature = feature.to_string();
        self
    }

    pub fn values(mut self, vep: &str, vep_rs: &str, change: i64) -> Self {
        self.vep = vep.to_string();
        self.vep_rs = vep_rs.to_string();
        self.change = change.to_string();
        self
    }

    pub fn vep(mut self, vep: &str) -> Self {
        self.vep = vep.to_string();
        self
    }

    pub fn vep_rs(mut self, vep_rs: &str) -> Self {
        self.vep_rs = vep_rs.to_string();
        self
    }

    pub fn consequence(mut self, consequence: &str) -> Self {
        self.consequence = consequence.to_string();
        self
    }

    pub fn strand(mut self, strand: i64) -> Self {
        self.strand = strand;
        self
    }

    pub fn offset(mut self, offset: i64) -> Self {
        self.offset = Some(offset);
        self
    }

    pub fn no_offset(mut self) -> Self {
        self.offset = None;
        self
    }

    /// The row as one engine writes it, with its DownstreamProtein `protein`.
    pub fn line(&self, protein: &str) -> String {
        downstream_frameshift_row(
            &self.location,
            &self.allele,
            &self.feature,
            &self.cds,
            protein,
            &self.change,
            &self.consequence,
            self.strand,
            self.offset,
        )
    }
}

/// A default-format row of the Downstream plugin: its CDS_position, its STRAND, an
/// HGVS_OFFSET where `offset` gives one, and the plugin's two fields.
#[allow(clippy::too_many_arguments)]
pub fn downstream_frameshift_row(
    location: &str,
    allele: &str,
    feature: &str,
    cds: &str,
    protein: &str,
    change: &str,
    consequence: &str,
    strand: i64,
    offset: Option<i64>,
) -> String {
    let shift = offset.map_or(String::new(), |offset| format!("HGVS_OFFSET={offset};"));
    let extra = format!(
        "IMPACT=HIGH;STRAND={strand};{shift}DownstreamProtein={protein};ProteinLengthChange={change}"
    );
    with_coding_cells(
        &default_row("v", location, allele, feature, consequence, &extra),
        &format!("-\t{cds}\t-\t-\t-\t-"),
    )
}

/// The reference's and the candidate's default-format text of `rows`.
pub fn downstream_texts(rows: &[DownstreamRow]) -> (String, String) {
    let reference: Vec<String> = rows.iter().map(|row| row.line(&row.vep)).collect();
    let candidate: Vec<String> = rows.iter().map(|row| row.line(&row.vep_rs)).collect();
    (default_text(&reference), default_text(&candidate))
}

/// One output's entries in each format a helper writes it in, under the format's name.
pub type Formats = Vec<(&'static str, Vec<Entry>)>;

/// The `INFO` column of `record`, a dot, replaced by the start of a CSQ value.
fn csq_line(record: &str) -> String {
    format!(
        "{}CSQ=",
        record
            .strip_suffix('.')
            .expect("the record's INFO is a dot")
    )
}

/// One forward-strand Downstream row of the input `record` parsed from the tab, VCF and
/// JSON formats; a `protein` of `None` is a row without the value: a tab `-`, an empty CSQ
/// cell, no JSON key.
pub fn downstream_formats(
    record: &str,
    location: &str,
    allele: &str,
    feature: &str,
    cds: &str,
    protein: Option<&str>,
    change: i64,
) -> Formats {
    let tab = format!(
        "{DOWNSTREAM_TAB_HEADER}\nv\t{location}\t{allele}\t{feature}\tTranscript\tframeshift_variant\t{cds}\t1\t{}\t{change}\n",
        protein.unwrap_or("-")
    );
    let vcf = format!(
        "{DOWNSTREAM_CSQ_HEADER}{}{allele}|frameshift_variant|Transcript|{feature}|{cds}|1|{}|{change}\n",
        csq_line(record),
        protein.unwrap_or("")
    );
    let (start, end) = cds.split_once('-').unwrap_or((cds, cds));
    let number = |bound: &str| whole(bound.parse().expect("a CDS bound is a number"));
    let mut entry = vec![
        ("variant_allele", text(allele)),
        ("transcript_id", text(feature)),
        ("consequence_terms", list(vec![text("frameshift_variant")])),
        ("cds_start", number(start)),
        ("cds_end", number(end)),
        ("strand", whole(1)),
        ("proteinlengthchange", whole(change)),
    ];
    if let Some(protein) = protein {
        entry.push(("downstreamprotein", text(protein)));
    }
    let json = dumped_value(&object(vec![
        ("input", text(record)),
        ("id", text("v")),
        ("most_severe_consequence", text("frameshift_variant")),
        ("transcript_consequences", list(vec![object(entry)])),
    ])) + "\n";
    let inputs = [record.to_string()];
    vec![
        ("tab", parse_tab(&tab).expect("the tab row parses")),
        (
            "vcf",
            parse_vcf(&vcf, &inputs).expect("the VCF line parses"),
        ),
        (
            "json",
            parse_json(&json, &inputs).expect("the JSON object parses"),
        ),
    ]
}

/// A plugin field's value as a row states it: a text or a whole number.
#[derive(Clone, Copy, Debug)]
pub enum Cell<'a> {
    Text(&'a str),
    Number(i64),
}

impl Cell<'_> {
    fn written(&self) -> String {
        match self {
            Cell::Text(value) => value.to_string(),
            Cell::Number(value) => value.to_string(),
        }
    }

    fn value(&self) -> Value {
        match self {
            Cell::Text(value) => text(value),
            Cell::Number(value) => whole(*value),
        }
    }
}

/// One forward-strand plugin row of the input `record` parsed from the four formats, with
/// an HGVS_OFFSET where `offset` gives one and the plugin fields `values`, a `None` a field
/// without a value: no Extra pair, a tab `-`, an empty CSQ cell, no JSON key. The JSON entry
/// carries the column's two numbers as cds_start and cds_end, lower first for an insertion
/// too, a `?` bound as no key.
#[allow(clippy::too_many_arguments)]
pub fn plugin_formats(
    record: &str,
    location: &str,
    allele: &str,
    feature: &str,
    cds: &str,
    values: &[(&str, Option<Cell>)],
    offset: Option<i64>,
    consequence: &str,
) -> Formats {
    let mut extra = vec!["IMPACT=HIGH".to_string(), "STRAND=1".to_string()];
    extra.extend(offset.map(|offset| format!("HGVS_OFFSET={offset}")));
    extra.extend(
        values
            .iter()
            .filter_map(|(name, value)| Some(format!("{name}={}", value.as_ref()?.written()))),
    );
    let default = default_text(&[with_coding_cells(
        &default_row(
            "v",
            location,
            allele,
            feature,
            consequence,
            &extra.join(";"),
        ),
        &format!("-\t{cds}\t-\t-\t-\t-"),
    )]);
    let mut columns = vec![
        "Allele",
        "Consequence",
        "Feature_type",
        "Feature",
        "CDS_position",
        "STRAND",
        "HGVS_OFFSET",
    ];
    columns.extend(values.iter().map(|(name, _)| *name));
    let mut cells: Vec<Option<String>> = vec![
        Some(allele.to_string()),
        Some(consequence.to_string()),
        Some("Transcript".to_string()),
        Some(feature.to_string()),
        Some(cds.to_string()),
        Some("1".to_string()),
        offset.map(|offset| offset.to_string()),
    ];
    cells.extend(
        values
            .iter()
            .map(|(_, value)| value.as_ref().map(Cell::written)),
    );
    let tab = format!(
        "#Uploaded_variation\tLocation\t{}\nv\t{location}\t{}\n",
        columns.join("\t"),
        cells
            .iter()
            .map(|cell| cell.as_deref().unwrap_or("-"))
            .collect::<Vec<_>>()
            .join("\t")
    );
    let csq = cells
        .iter()
        .map(|cell| cell.as_deref().unwrap_or("").replace(',', "&"))
        .collect::<Vec<_>>()
        .join("|");
    let vcf = format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {}\">\n{}{csq}\n",
        columns.join("|"),
        csq_line(record)
    );
    let (start, end) = cds.split_once('-').unwrap_or((cds, cds));
    let mut entry = vec![
        ("variant_allele".to_string(), text(allele)),
        ("transcript_id".to_string(), text(feature)),
        (
            "consequence_terms".to_string(),
            list(consequence.split(',').map(text).collect()),
        ),
        ("strand".to_string(), whole(1)),
    ];
    for (key, bound) in [("cds_start", start), ("cds_end", end)] {
        if let Ok(number) = bound.parse() {
            entry.push((key.to_string(), whole(number)));
        }
    }
    if let Some(offset) = offset {
        entry.push(("hgvs_offset".to_string(), whole(offset)));
    }
    for (name, value) in values {
        if let Some(value) = value {
            entry.push((name.to_lowercase(), value.value()));
        }
    }
    let json = dumped_value(&object(vec![
        ("input", text(record)),
        ("id", text("v")),
        (
            "most_severe_consequence",
            text(consequence.split(',').next().unwrap_or(consequence)),
        ),
        (
            "transcript_consequences",
            list(vec![Value::Dict(entry.into_iter().collect())]),
        ),
    ])) + "\n";
    let inputs = [record.to_string()];
    vec![
        (
            "default",
            parse_default(&default).expect("the default row parses"),
        ),
        ("tab", parse_tab(&tab).expect("the tab row parses")),
        (
            "vcf",
            parse_vcf(&vcf, &inputs).expect("the VCF line parses"),
        ),
        (
            "json",
            parse_json(&json, &inputs).expect("the JSON object parses"),
        ),
    ]
}

/// A `Documented` of no manifest whose mask set names `classes`.
pub fn built_in(classes: &[&str]) -> Documented {
    let masks = classes.iter().map(|class| class.to_string()).collect();
    Documented::from_manifest(None, &masks, None, None, None)
        .expect("no manifest documents nothing")
}

/// The report of two entry lists under `documented` beside the side tables of `side`,
/// with the lines of the input where the comparison is handed them.
pub fn compared_under(
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
    side: &Side,
    inputs: Option<&[String]>,
) -> Value {
    let options = CompareOptions {
        input_lines: inputs,
        ..CompareOptions::new(side)
    };
    compare_entries(expected, actual, documented, &options).expect("the entries compare")
}

/// The report of two entry lists under the built-in `classes`, and what was applied.
pub fn compared(
    expected: &[Entry],
    actual: &[Entry],
    side: &Side,
    classes: &[&str],
) -> (Value, Documented) {
    let mut documented = built_in(classes);
    let report = compared_under(expected, actual, &mut documented, side, None);
    (report, documented)
}

/// The report of two default-format texts under the built-in `classes`, and what was
/// applied.
pub fn compare_downstream(
    reference: &str,
    candidate: &str,
    side: &Side,
    classes: &[&str],
) -> (Value, Documented) {
    compared(
        &parse_default(reference).expect("the reference parses"),
        &parse_default(candidate).expect("the candidate parses"),
        side,
        classes,
    )
}

/// What the comparison of two default-format texts under the built-in `classes` says where
/// it ends on a refusal.
pub fn refusal_of(reference: &str, candidate: &str, side: &Side, classes: &[&str]) -> String {
    compare_entries(
        &parse_default(reference).expect("the reference parses"),
        &parse_default(candidate).expect("the candidate parses"),
        &mut built_in(classes),
        &CompareOptions::new(side),
    )
    .expect_err("the comparison is refused")
    .to_string()
}

/// The report of two default-format texts of the input records `inputs` under the built-in
/// `classes`, and what was applied.
pub fn compare_record(
    reference: &str,
    candidate: &str,
    side: &Side,
    classes: &[&str],
    inputs: &[&str],
) -> (Value, Documented) {
    let inputs: Vec<String> = inputs.iter().map(|line| line.to_string()).collect();
    let mut documented = built_in(classes);
    let report = compared_under(
        &parse_default(reference).expect("the reference parses"),
        &parse_default(candidate).expect("the candidate parses"),
        &mut documented,
        side,
        Some(&inputs),
    );
    (report, documented)
}

/// The verdict and the applied counts of each format's comparison under `classes`.
pub fn applied_by_format(
    expected: &Formats,
    actual: &Formats,
    side: &Side,
    classes: &[&str],
) -> Vec<(&'static str, String, BTreeMap<String, i64>)> {
    expected
        .iter()
        .zip(actual)
        .map(|((fmt, entries), (_, others))| {
            let (report, documented) = compared(entries, others, side, classes);
            (
                *fmt,
                verdict_of(&report).to_string(),
                counted(&documented.applied),
            )
        })
        .collect()
}

/// Whether every format of `compared` reads `verdict` with the applied counts `applied`.
pub fn every_format_reads(
    compared: &[(&'static str, String, BTreeMap<String, i64>)],
    verdict: &str,
    applied: &[(&str, i64)],
) -> bool {
    !compared.is_empty()
        && compared
            .iter()
            .all(|(_, read, counts)| read == verdict && *counts == count_map(applied))
}

/// The member of `value` under the names of `path`, one object deep a name.
pub fn member<'a>(value: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(value, |value, name| {
        value
            .get(name)
            .unwrap_or_else(|| panic!("no member {name} of {path:?}"))
    })
}

pub fn verdict_of(report: &Value) -> &str {
    member(report, &["verdict"])
        .as_str()
        .expect("a verdict is text")
}

/// The counts of `counter` other than zero.
pub fn counted(counter: &Counter) -> BTreeMap<String, i64> {
    counter
        .iter()
        .filter(|(_, count)| *count != 0)
        .map(|(name, count)| (name.to_string(), count))
        .collect()
}

pub fn count_map(items: &[(&str, i64)]) -> BTreeMap<String, i64> {
    items
        .iter()
        .map(|(name, count)| (name.to_string(), *count))
        .collect()
}

/// The counts an object of a report holds by name.
pub fn counts_of(value: &Value) -> BTreeMap<String, i64> {
    value
        .as_dict()
        .expect("an object of counts")
        .iter()
        .map(|(name, count)| match count {
            Value::Int(count) => (name.clone(), *count),
            other => panic!("the count of {name} is {other:?}"),
        })
        .collect()
}

/// The rows of `column` set aside by class; none for a column without such rows.
pub fn masked(report: &Value, column: &str) -> BTreeMap<String, i64> {
    member(report, &["columns", column])
        .get("masked")
        .map(counts_of)
        .unwrap_or_default()
}

/// The differing values of `column`.
pub fn differing(report: &Value, column: &str) -> i64 {
    match member(report, &["columns", column, "differing"]) {
        Value::Int(count) => *count,
        other => panic!("the differing count of {column} is {other:?}"),
    }
}

/// The texts of a list of a report.
pub fn texts_of(value: &Value) -> Vec<&str> {
    value
        .as_list()
        .expect("a list")
        .iter()
        .map(|item| item.as_str().expect("a text"))
        .collect()
}

/// The examples a report prints for the rows `class` set aside.
pub fn mask_examples<'a>(report: &'a Value, class: &str) -> Vec<&'a str> {
    texts_of(member(report, &["mask_examples", class]))
}

/// The examples a report prints for the differing values of `column` under the
/// consequence set `terms`; none where it prints none.
pub fn column_examples<'a>(report: &'a Value, column: &str, terms: &str) -> Vec<&'a str> {
    member(report, &["columns", column])
        .get("examples")
        .and_then(|examples| examples.get(terms))
        .map(texts_of)
        .unwrap_or_default()
}

/// A default-format row of the NMD plugin: its CDS_position, an HGVS_OFFSET where `offset`
/// gives one and the plugin's field where `nmd` gives its value.
pub fn nmd_row(
    location: &str,
    allele: &str,
    feature: &str,
    cds: &str,
    nmd: Option<&str>,
    consequence: &str,
    offset: Option<i64>,
) -> String {
    let mut extra = vec!["IMPACT=HIGH".to_string(), "STRAND=1".to_string()];
    extra.extend(offset.map(|offset| format!("HGVS_OFFSET={offset}")));
    extra.extend(nmd.map(|nmd| format!("NMD={nmd}")));
    with_coding_cells(
        &default_row(
            "v",
            location,
            allele,
            feature,
            consequence,
            &extra.join(";"),
        ),
        &format!("-\t{cds}\t-\t-\t-\t-"),
    )
}

/// An empty directory of the test `name`, under the build's own directory for test files.
pub fn scratch_of(name: &str) -> PathBuf {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("downstream-{name}"));
    let _ = std::fs::remove_dir_all(&directory);
    std::fs::create_dir_all(&directory).expect("a directory for the test's files");
    directory
}

/// How one run of the binary ended.
#[derive(Debug)]
pub struct Ran {
    pub code: i32,
    pub stderr: String,
}

/// The binary on `args`.
pub fn run(args: &[String]) -> Ran {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_vep-formats-compare"))
        .args(args)
        .output()
        .expect("the binary runs");
    Ran {
        code: output.status.code().expect("the binary exits"),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// The arguments of a default-format comparison of `reference` and `candidate`, written
/// under `directory`, whose report is `report.json` there.
pub fn default_arguments(directory: &Path, reference: &str, candidate: &str) -> Vec<String> {
    let reference = written(directory, "ref.txt", reference);
    let candidate = written(directory, "cand.txt", candidate);
    [
        "--format",
        "default",
        "--expected",
        path_text(&reference),
        "--actual",
        path_text(&candidate),
        "--report",
        path_text(&directory.join("report.json")),
    ]
    .iter()
    .map(|argument| argument.to_string())
    .collect()
}

/// `arguments` with `more` after them.
pub fn with(arguments: &[String], more: &[&str]) -> Vec<String> {
    let mut all = arguments.to_vec();
    all.extend(more.iter().map(|argument| argument.to_string()));
    all
}

pub fn path_text(path: &Path) -> &str {
    path.to_str().expect("the path is UTF-8")
}

/// The report the last run wrote under `directory`.
pub fn report_of(directory: &Path) -> Value {
    let text = std::fs::read_to_string(directory.join("report.json")).expect("the report is read");
    pyjson::loads(&text).expect("the report is JSON")
}
