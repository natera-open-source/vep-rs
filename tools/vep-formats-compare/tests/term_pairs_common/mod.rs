// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the term-pair classes share: the planted rows of a stop codon and of a
//! record whose span lies beside an intron, in each format.

#![allow(dead_code)]

use crate::compare_common::default_text;
use vep_formats_compare::consts::CONSEQUENCE_RANK;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Result;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::rules::coding_end_pairs::{AtACodingEnd, RowPair};
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

pub type Pair = fn(&Entry, &Entry) -> Result<bool>;

/// Whether `pair` takes the reference entry `e` and the candidate entry `a`.
pub fn holds(pair: Pair, e: &Entry, a: &Entry) -> bool {
    pair(e, a).expect("the pair is read")
}

/// Whether the row rule `pair` takes the two entries in a comparison that was given no
/// cache.
pub fn holds_without_a_cache(pair: RowPair, e: &Entry, a: &Entry) -> bool {
    pair(e, a, &AtACodingEnd::NO_CACHE).expect("the pair is read")
}

/// The one entry of a default-format output of `row`.
pub fn default_entry(row: &str) -> Entry {
    let mut entries = parse_default(&default_text(&[row])).expect("the row parses");
    assert_eq!(entries.len(), 1, "{row}");
    entries.remove(0)
}

/// The value of the entry's field `name`, `ABSENT` where it carries none.
pub fn field<'a>(e: &'a Entry, name: &str) -> &'a str {
    e.fields.get(name).map_or(ABSENT, String::as_str)
}

/// Whether the field class `klass` accepts the field `name` of the reference entry `e`
/// beside the candidate entry `a`, nothing else named.
pub fn class_accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    a: Option<&Entry>,
) -> bool {
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(klass, e, name, expected_value, actual_value, a, &context)
        .expect("the class is asked")
}

pub const STOP_CODON_LOCATION: &str = "19:13208046-13208050";

/// A Transcript row of the record `1993671` at the Location of a deletion over a stop codon.
pub fn stop_codon_row(allele: &str, consequence: &str, cdna: &str, cds: &str) -> String {
    stop_codon_row_at(allele, consequence, cdna, cds, STOP_CODON_LOCATION)
}

pub fn stop_codon_row_at(
    allele: &str,
    consequence: &str,
    cdna: &str,
    cds: &str,
    location: &str,
) -> String {
    stop_codon_row_of(allele, consequence, cdna, cds, location, "Transcript")
}

pub fn stop_codon_row_of(
    allele: &str,
    consequence: &str,
    cdna: &str,
    cds: &str,
    location: &str,
    feature_type: &str,
) -> String {
    format!(
        "1993671\t{location}\t{allele}\tENSG1\tENST00000635895\t{feature_type}\t{consequence}\t{cdna}\t{cds}\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1"
    )
}

/// One row of `stop_codon_row_at` as a default-format output, its Extra text `extra`.
pub fn stop_codon_text(
    terms: &str,
    cdna: &str,
    cds: &str,
    extra: &str,
    allele: &str,
    location: &str,
) -> String {
    default_text(&[stop_codon_row_at(allele, terms, cdna, cds, location)
        .replace("IMPACT=LOW;STRAND=-1", extra)])
}

pub const RETAINED: &str = "stop_retained_variant,3_prime_UTR_variant";
pub const UNKNOWN: &str = "coding_sequence_variant,3_prime_UTR_variant";
pub const SHIFTED: &str = "frameshift_variant,3_prime_UTR_variant";

/// The header of a VCF output whose CSQ entries carry `fields`.
pub fn csq_header(fields: &str) -> String {
    format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {fields}\">\n"
    )
}

/// An input line with its INFO column left open for a CSQ value.
pub fn open_info(input_line: &str) -> &str {
    input_line
        .strip_suffix('.')
        .expect("an input line ends in an empty INFO column")
}

pub fn quoted(text: &str) -> String {
    pyjson::dumps_str(text, true)
}

/// A JSON array of the texts `items`.
pub fn quoted_list(items: &[&str]) -> String {
    let items: Vec<String> = items.iter().map(|item| quoted(item)).collect();
    format!("[{}]", items.join(", "))
}

/// The entries of a VCF output with a feature.
pub fn vcf_feature_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut entries = parse_vcf(text, inputs).expect("the output parses");
    entries.retain(|entry| !entry.feature.is_empty());
    entries
}

/// The entries of a JSON output with a feature.
pub fn json_feature_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut entries = parse_json(text, inputs).expect("the output parses");
    entries.retain(|entry| !entry.feature.is_empty());
    entries
}

/// `JSON.pm:387-388`: the lowest-ranked of `terms`.
pub fn most_severe<'a>(terms: &[&'a str]) -> &'a str {
    let rank = |term: &str| {
        CONSEQUENCE_RANK
            .iter()
            .find(|(ranked, _)| *ranked == term)
            .map_or(1_000_000, |(_, rank)| *rank)
    };
    terms
        .iter()
        .copied()
        .min_by_key(|term| (rank(term), *term))
        .unwrap_or("?")
}

pub const SPLICE_SPAN_CLASS: &str = "splice_term_past_the_record_span";
pub const SPLICE_SPAN_PLUGIN_VALUE_CLASS: &str =
    "plugin_value_follows_splice_term_past_the_record_span";
pub const PLUGIN_VALUE_CLASS: &str = "plugin_value_follows_stop_codon_term_pair";

/// The two term sets of a row of the class: the reference's donor term gone, and the
/// candidate's `splice_region_variant` in its place.
pub const DONOR_PAST_THE_SPAN: &str = "splice_donor_variant,frameshift_variant";
pub const REGION_IN_ITS_PLACE: &str = "frameshift_variant,splice_region_variant";
/// The same beside a change that keeps the reading frame, so that the two sets differ in
/// impact.
pub const DONOR_PAST_THE_SPAN_IN_FRAME: &str = "splice_donor_variant,protein_altering_variant";
pub const REGION_IN_ITS_PLACE_IN_FRAME: &str = "protein_altering_variant,splice_region_variant";

/// The columns of one Transcript row of a made-up record that replaces the three bases of
/// an exon at `location` with `allele`.
#[derive(Debug, Clone, Copy)]
pub struct SpanColumns<'a> {
    pub cdna: &'a str,
    pub cds: &'a str,
    pub allele: &'a str,
    pub location: &'a str,
    pub name: &'a str,
    pub feature_type: &'a str,
    pub extra: &'a str,
}

impl Default for SpanColumns<'_> {
    fn default() -> Self {
        SpanColumns {
            cdna: "998-1000",
            cds: "898-900",
            allele: "ACGTA",
            location: "7:100998-101000",
            name: "rec1",
            feature_type: "Transcript",
            extra: "IMPACT=HIGH;STRAND=1",
        }
    }
}

pub fn splice_span_row(terms: &str, columns: &SpanColumns) -> String {
    let SpanColumns {
        cdna,
        cds,
        allele,
        location,
        name,
        feature_type,
        extra,
    } = columns;
    format!(
        "{name}\t{location}\t{allele}\tENSG1\tENST90000000001\t{feature_type}\t{terms}\t{cdna}\t{cds}\t-\t-\t-\t-\t{extra}"
    )
}

/// A made-up record the class reads with its alleles: the input line's CHROM, POS, REF and
/// ALT, then what the default format writes for its row, the Location, the Allele, the cDNA
/// and CDS positions, and the reference's and the candidate's terms and IMPACT.
#[derive(Debug, Clone)]
pub struct SpliceSpanRecord {
    pub line: [String; 4],
    pub location: String,
    pub allele: String,
    pub cdna: String,
    pub cds: String,
    pub terms: (String, String),
    pub impacts: (String, String),
}

pub const SPLICE_SPAN_CSQ_FIELDS: &str =
    "Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|CDS_position";

impl SpliceSpanRecord {
    pub fn new(
        line: [&str; 4],
        location: &str,
        allele: &str,
        cdna: &str,
        cds: &str,
        terms: (&str, &str),
    ) -> Self {
        SpliceSpanRecord {
            line: line.map(str::to_string),
            location: location.to_string(),
            allele: allele.to_string(),
            cdna: cdna.to_string(),
            cds: cds.to_string(),
            terms: (terms.0.to_string(), terms.1.to_string()),
            impacts: ("HIGH".to_string(), "HIGH".to_string()),
        }
    }

    pub fn with_impacts(mut self, reference: &str, candidate: &str) -> Self {
        self.impacts = (reference.to_string(), candidate.to_string());
        self
    }

    pub fn with_terms(mut self, reference: &str, candidate: &str) -> Self {
        self.terms = (reference.to_string(), candidate.to_string());
        self
    }

    pub fn with_line(mut self, line: [&str; 4]) -> Self {
        self.line = line.map(str::to_string);
        self
    }

    /// The record's line of the input VCF under the identifier `name`.
    pub fn input_line(&self, name: &str) -> String {
        let [chromosome, position, reference, alternate] = &self.line;
        format!("{chromosome}\t{position}\t{name}\t{reference}\t{alternate}\t.\t.\t.")
    }

    /// The first column of the record's rows where the input line carries no identifier.
    pub fn spelt(&self) -> String {
        let [chromosome, position, reference, alternate] = &self.line;
        format!("{chromosome}_{position}_{reference}/{alternate}")
    }

    /// The columns of `splice_span_row` the record gives its row under the first column
    /// `name`.
    pub fn columns<'a>(&'a self, name: &'a str) -> SpanColumns<'a> {
        SpanColumns {
            cdna: &self.cdna,
            cds: &self.cds,
            allele: &self.allele,
            location: &self.location,
            name,
            ..SpanColumns::default()
        }
    }

    /// The reference's and the candidate's default row, under the identifier `name` or
    /// under the first column that spells the record's alleles.
    pub fn rows(&self, name: Option<&str>) -> (String, String) {
        let spelt = self.spelt();
        let name = name.unwrap_or(&spelt);
        let row = |terms: &str, impact: &str| {
            let extra = format!("IMPACT={impact};STRAND=1");
            splice_span_row(
                terms,
                &SpanColumns {
                    extra: &extra,
                    ..self.columns(name)
                },
            )
        };
        (
            row(&self.terms.0, &self.impacts.0),
            row(&self.terms.1, &self.impacts.1),
        )
    }

    pub fn entries(&self, name: Option<&str>) -> (Entry, Entry) {
        let (reference, candidate) = self.rows(name);
        (default_entry(&reference), default_entry(&candidate))
    }

    /// The reference's and the candidate's output of the record's row as the parser of
    /// `fmt` reads them, the default row under the first column that spells the alleles.
    pub fn parsed(&self, fmt: &str) -> (Vec<Entry>, Vec<Entry>) {
        let inputs = [self.input_line(".")];
        let side = |side: usize| match fmt {
            "default" => {
                let rows = self.rows(None);
                let row = if side == 0 { rows.0 } else { rows.1 };
                parse_default(&default_text(&[row])).expect("the row parses")
            }
            "vcf" => parse_vcf(&self.vcf(side), &inputs).expect("the output parses"),
            _ => parse_json(&self.json(side), &inputs).expect("the output parses"),
        };
        (side(0), side(1))
    }

    /// The reference's and the candidate's entry of the record's row in `fmt` (`parsed`).
    pub fn in_format(&self, fmt: &str) -> (Entry, Entry) {
        let (reference, candidate) = self.parsed(fmt);
        let with_feature = |entries: Vec<Entry>| {
            entries
                .into_iter()
                .find(|entry| !entry.feature.is_empty())
                .expect("an entry with a feature")
        };
        (with_feature(reference), with_feature(candidate))
    }

    fn side(&self, side: usize) -> (&str, &str) {
        if side == 0 {
            (&self.terms.0, &self.impacts.0)
        } else {
            (&self.terms.1, &self.impacts.1)
        }
    }

    /// The VCF output of one side, its CSQ entry the row's fields.
    pub fn vcf(&self, side: usize) -> String {
        let (terms, impact) = self.side(side);
        let csq = [
            self.allele.as_str(),
            &terms.replace(',', "&"),
            impact,
            "ENST90000000001",
            "Transcript",
            written(&self.cdna),
            written(&self.cds),
        ]
        .join("|");
        let input_line = self.input_line(".");
        format!(
            "{}{}CSQ={csq}\n",
            csq_header(SPLICE_SPAN_CSQ_FIELDS),
            open_info(&input_line)
        )
    }

    /// The JSON output of one side: the span and the allele string as Ensembl VEP stores
    /// them, without the first base the record's alleles share and an insertion's start one
    /// above its end.
    pub fn json(&self, side: usize) -> String {
        let (terms, impact) = self.side(side);
        let (reference, alternate) = (&self.line[2][1..], &self.line[3][1..]);
        let coordinates = self.location.split_once(':').map_or("", |(_, span)| span);
        let (low, high) = coordinates.split_once('-').unwrap_or((coordinates, ""));
        let high = if high.is_empty() { low } else { high };
        let (start, end) = if reference.is_empty() {
            (high, low)
        } else {
            (low, high)
        };
        let term_list: Vec<&str> = terms.split(',').collect();
        let mut entry = format!(
            "{{\"variant_allele\": {}, \"transcript_id\": \"ENST90000000001\", \"consequence_terms\": {}, \"impact\": {}",
            quoted(&self.allele),
            quoted_list(&term_list),
            quoted(impact)
        );
        for (key, column) in [("cdna", &self.cdna), ("cds", &self.cds)] {
            let (first, last) = column.split_once('-').unwrap_or((column, ""));
            for (bound, text) in [("start", first), ("end", last)] {
                if !text.is_empty() && text.bytes().all(|digit| digit.is_ascii_digit()) {
                    entry.push_str(&format!(", \"{key}_{bound}\": {text}"));
                }
            }
        }
        entry.push('}');
        format!(
            "{{\"input\": {}, \"seq_region_name\": {}, \"start\": {start}, \"end\": {end}, \"allele_string\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
            quoted(&self.input_line(".")),
            quoted(&self.line[0]),
            quoted(&format!(
                "{}/{alternate}",
                if reference.is_empty() { "-" } else { reference }
            )),
            quoted(most_severe(&term_list)),
        )
    }
}

/// A position as a CSQ entry writes it: empty where the row has none.
fn written(position: &str) -> &str {
    if position == "-" {
        ""
    } else {
        position
    }
}

/// Three records whose ALT is their REF and three further bases, and whose rows leave out the
/// first base the two share: each allele opens with its whole reference allele, and the
/// candidate's row is that of the insertion between the span's last base and the base after
/// it. The first span ends one base short of an exon's end: the bases past it reach the donor
/// window of the intron after it, and the insertion stands beside the exon's end. The second
/// ends on the first base of an intron: the insertion stands between the two donor bases,
/// which gives neither the donor region nor the intron term. The third lies in an intron short
/// of its acceptor window, which the bases past it reach.
pub fn splice_span_records() -> Vec<(&'static str, SpliceSpanRecord)> {
    vec![
        (
            "two positions inside an exon",
            SpliceSpanRecord::new(
                ["7", "100998", "ACG", "ACGTAC"],
                "7:100999-101000",
                "CGTAC",
                "999-1000",
                "899-900",
                (DONOR_PAST_THE_SPAN, REGION_IN_ITS_PLACE),
            ),
        ),
        (
            "a span across an exon's end",
            SpliceSpanRecord::new(
                ["7", "200997", "TACG", "TACGTAC"],
                "7:200998-201000",
                "ACGTAC",
                "999-?",
                "899-?",
                (
                    "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                    "splice_donor_variant,coding_sequence_variant",
                ),
            ),
        ),
        (
            "a span in an intron",
            SpliceSpanRecord::new(
                ["7", "300998", "TCG", "TCGTAC"],
                "7:300999-301000",
                "CGTAC",
                "-",
                "-",
                (
                    "splice_acceptor_variant,splice_polypyrimidine_tract_variant,intron_variant",
                    "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
                ),
            )
            .with_impacts("HIGH", "LOW"),
        ),
    ]
}

/// The record of `splice_span_records` named `shape`.
pub fn splice_span_record(shape: &str) -> SpliceSpanRecord {
    splice_span_records()
        .into_iter()
        .find(|(name, _)| *name == shape)
        .map(|(_, record)| record)
        .expect("a record of that shape")
}

/// A pure insertion between two positions.
pub fn inserted_between_two() -> SpliceSpanRecord {
    SpliceSpanRecord::new(
        ["7", "400999", "C", "CGTA"],
        "7:400999-401000",
        "GTA",
        "999-1000",
        "899-900",
        (DONOR_PAST_THE_SPAN, REGION_IN_ITS_PLACE),
    )
}

/// A record that replaces the two bases at the positions of `inserted_between_two` by the
/// same allele: one Location, Allele, positions and pair of term sets, which the record
/// alone tells apart.
pub fn two_replaced() -> SpliceSpanRecord {
    inserted_between_two().with_line(["7", "400998", "ACG", "AGTA"])
}
