// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The GWAS catalog classes.

use crate::argparse::path_text;
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::fields::{entry_kind, field_named, record_label};
use crate::input::{line_slices, read_text};
use crate::pyfmt;
use crate::rules::revel::{decimal_of, decimal_plus};
use crate::rules::spelt_in;
use crate::shared::{blank, gvcf_site_rows};
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::LazyLock;

pub const GWAS_PREFIX: &str = "gwas_";

pub const GWAS_CATALOG_COLUMNS: &[&str] = &[
    "CHR_ID",
    "CHR_POS",
    "STRONGEST SNP-RISK ALLELE",
    "STUDY ACCESSION",
];

pub const GWAS_SNPS_COLUMN: &str = "SNPS";

pub const GWAS_LIST_SPLIT_RE_PATTERN: &str = r"[;,x]";
pub static GWAS_LIST_SPLIT_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(GWAS_LIST_SPLIT_RE_PATTERN).expect("the pattern compiles"));

pub const GWAS_RS_ID_RE_PATTERN: &str = r"rs[0-9]+";
pub static GWAS_RS_ID_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(GWAS_RS_ID_RE_PATTERN).expect("the pattern compiles"));

/// A place and what is paired there: the chromosome without `chr`, the position, the risk
/// allele and the study accession.
pub type GwasPair = (String, String, String, String);

/// Whether `last - first + 1` is `bases`, the two numbers decimal texts of any length
/// without a sign.
fn spans(first: &str, last: &str, bases: usize) -> bool {
    match bases.checked_sub(1) {
        Some(more) => last == decimal_plus(first, more as u64),
        None => first == decimal_plus(last, 1),
    }
}

/// The places where the two engines read a GWAS catalog's rows of two or more SNPs
/// differently, each a (chromosome without `chr`, position, risk allele, study accession):
/// `reference_only` where Ensembl VEP hands a location another SNP's risk allele
/// (GWAS.pm:338-340), `candidate_only` where it keeps the space after an allele written
/// before ` x ` and so never matches the SNP's own (GWAS.pm:343-348).
#[derive(Debug)]
pub struct GwasCatalogPairs {
    pub reference_only: BTreeSet<GwasPair>,
    pub candidate_only: BTreeSet<GwasPair>,
}

impl GwasCatalogPairs {
    pub fn new() -> Self {
        GwasCatalogPairs {
            reference_only: BTreeSet::new(),
            candidate_only: BTreeSet::new(),
        }
    }

    /// The pairs of the associations file `path`, read by its header's column names.
    pub fn from_catalog(path: &Path) -> Result<Self> {
        let text = read_text(path)?;
        let mut lines = line_slices(&text);
        let header: Vec<&str> = lines
            .next()
            .map_or_else(Vec::new, |line| line.split('\t').collect());
        let column = |name: &str| header.iter().position(|held| *held == name);
        let columns: Option<Vec<usize>> = GWAS_CATALOG_COLUMNS
            .iter()
            .map(|name| column(name))
            .collect();
        let Some(columns) = columns else {
            return Err(Error::comparator(format!(
                "{}: not a GWAS catalog associations file (columns {})",
                path_text(path),
                GWAS_CATALOG_COLUMNS.join(", ")
            )));
        };
        let snps = column(GWAS_SNPS_COLUMN);
        let last = columns.iter().copied().max().unwrap_or(0);
        let mut pairs = Self::new();
        for line in lines {
            let cells: Vec<&str> = line.split('\t').collect();
            if cells.len() > last {
                pairs.add_row(
                    cells[columns[0]],
                    cells[columns[1]],
                    cells[columns[2]],
                    cells[columns[3]],
                    snps.and_then(|at| cells.get(at)).copied().unwrap_or(""),
                );
            }
        }
        Ok(pairs)
    }

    /// One catalog row: its lists split as the plugin splits them, the locations
    /// stripped, each risk allele the text after its token's first `-` as written. A row
    /// of several locations whose SNPs, the `rs` ids of `snps` or else its risk alleles,
    /// differ from them in count is one the candidate skips, while the reference writes
    /// every SNP at every location (GWAS.pm:237 compares the chromosomes with the positions
    /// alone), so each of its alleles is the reference's alone at each location. A row
    /// whose chromosomes and positions differ in count, or whose risk alleles differ in
    /// count from as many SNPs as locations, pairs nothing. Neither does a row whose
    /// risk-allele cell holds a `?`: both engines read that cell as empty (GWAS.pm:261)
    /// and write nothing from the row.
    ///
    /// `snps` is the row's `SNPS` cell, empty for a row that has none.
    pub fn add_row(
        &mut self,
        chromosomes: &str,
        positions: &str,
        risk_alleles: &str,
        study: &str,
        snps: &str,
    ) {
        let chromosome_list: Vec<&str> = GWAS_LIST_SPLIT_RE
            .split(chromosomes)
            .map(|chromosome| pyfmt::removeprefix(pyfmt::strip(chromosome), "chr"))
            .collect();
        let position_list: Vec<&str> = GWAS_LIST_SPLIT_RE
            .split(positions)
            .map(pyfmt::strip)
            .collect();
        // A token's allele ends at its second `-`.
        let written: Vec<&str> = GWAS_LIST_SPLIT_RE
            .split(risk_alleles)
            .map(|token| token.split('-').nth(1).unwrap_or(""))
            .collect();
        if position_list.len() < 2
            || chromosome_list.len() != position_list.len()
            || risk_alleles.contains('?')
        {
            return;
        }
        let places = || chromosome_list.iter().zip(&position_list);
        let pair = |chromosome: &str, position: &str, allele: &str| -> GwasPair {
            (
                chromosome.to_string(),
                position.to_string(),
                allele.to_string(),
                study.to_string(),
            )
        };
        let matched_as_written: BTreeSet<&str> = written
            .iter()
            .copied()
            .filter(|allele| spelt_in(allele, b"ACGT"))
            .collect();
        let named = match GWAS_RS_ID_RE.find_iter(snps).count() {
            0 => written.len(),
            ids => ids,
        };
        if named != position_list.len() {
            for (chromosome, position) in places() {
                for allele in &matched_as_written {
                    self.reference_only
                        .insert(pair(chromosome, position, allele));
                }
            }
            return;
        }
        if written.len() != position_list.len() {
            return;
        }
        for ((chromosome, position), allele) in places().zip(&written) {
            let own = Some(pyfmt::strip(allele)).filter(|allele| spelt_in(allele, b"ACGT"));
            for other in &matched_as_written {
                if own != Some(*other) {
                    self.reference_only
                        .insert(pair(chromosome, position, other));
                }
            }
            if let Some(own) = own.filter(|own| !matched_as_written.contains(own)) {
                self.candidate_only.insert(pair(chromosome, position, own));
            }
        }
    }
}

impl Default for GwasCatalogPairs {
    fn default() -> Self {
        Self::new()
    }
}

/// The chromosome, without `chr`, and position at which the plugin matches a row's
/// allele (GWAS.pm:176-195 compare the record's alleles at its start with a risk allele at
/// the catalog's position): in the default and tab formats a one-base Location's position, or
/// the first position of a Location spanning as many bases as the row's sequence allele, a
/// substitution of several bases, which an insertion the parser read from the record's
/// alleles is not; else the `CHROM` of the VCF format's record entry with the start the
/// parser gives the row's allele from the line's `POS`, `REF` and `ALT`
/// (`shared::gvcf_site_rows`; the `POS` itself for an allele it does
/// not place), or the `seq_region_name` and `start` of the JSON format's among
/// `record_rows`. A VCF or JSON row is held to the Location's condition where its record
/// states the span the allele replaces, the parser's start and end of the row's allele or the
/// JSON object's `start` and `end`: a span of another length than one base is read only
/// for a sequence allele of as many bases that is no insertion's or deletion's, so an
/// insertion's row, whose span is empty, states no site. `None` where the row states none.
///
/// A position is a whole number of any length; a text `str.isdigit` takes and `int` does
/// not convert ends the run as the `ValueError` does.
pub fn gwas_site(e: &Entry, record_rows: Option<&[&Entry]>) -> Option<(String, String)> {
    let allele_length = pyfmt::len(&e.allele);
    // A span of several bases is read for a substitution of as many bases alone.
    let substitution = spelt_in(&e.allele, b"ACGTN") && e.indel != Some(true);
    if let Some(location) = e.location.as_deref() {
        let (chromosome, _, span) = pyfmt::rpartition(location, ":");
        let (start, dash_between, end) = pyfmt::partition(span, "-");
        if chromosome.is_empty() || !pyfmt::is_digit(start) {
            return None;
        }
        if !dash_between.is_empty() {
            if !(substitution && pyfmt::is_digit(end)) {
                return None;
            }
            // The end is converted before the start.
            let last = decimal_of(end);
            if !spans(&decimal_of(start), &last, allele_length) {
                return None;
            }
        }
        return Some((
            pyfmt::removeprefix(chromosome, "chr").to_string(),
            start.to_string(),
        ));
    }
    let label = record_label(e);
    let record = record_rows
        .unwrap_or(&[])
        .iter()
        .find(|r| entry_kind(r) == "record" && record_label(r) == label)?;
    let first_of = |names: [&str; 2]| {
        names
            .iter()
            .find_map(|name| record.fields.get(*name))
            .map_or("", String::as_str)
    };
    let chromosome = first_of(["CHROM", "seq_region_name"]);
    let mut position = first_of(["POS", "start"]).to_string();
    // Whether the span the record states is one base, and whether it is as long as the
    // row's allele; `None` where the record states none.
    let mut span: Option<(bool, bool)> = None;
    if record.fields.contains_key("POS") && pyfmt::is_digit(&position) {
        // Parser/VCF.pm:300-336 drop the first base an indel record's alleles share and
        // move its start one base on, so the row's start is not the line's POS there. The
        // rows are placed from position 0 and moved by the line's POS, a number of any
        // length.
        let line_position = decimal_of(&position);
        let reference = record.fields.get("REF").map_or("", String::as_str);
        let alternates: Vec<&str> = record
            .fields
            .get("ALT")
            .map_or("", String::as_str)
            .split(',')
            .collect();
        let allele = pyfmt::upper(&e.allele);
        let placed: Vec<(i64, i64)> = gvcf_site_rows(0, reference, &alternates)
            .into_iter()
            .filter(|(spelt, _, _)| *spelt == allele)
            .map(|(_, start, end)| (start, end))
            .collect();
        let starts: BTreeSet<i64> = placed.iter().map(|(start, _)| *start).collect();
        if let (1, Some(start)) = (starts.len(), starts.first()) {
            position = decimal_plus(&line_position, start.unsigned_abs());
            let lengths: BTreeSet<i64> =
                placed.iter().map(|(start, end)| end - start + 1).collect();
            if let (1, Some(bases)) = (lengths.len(), lengths.first()) {
                span = Some((*bases == 1, *bases == allele_length as i64));
            }
        }
    } else if pyfmt::is_digit(&position) {
        if let Some(end) = record.fields.get("end").filter(|end| pyfmt::is_digit(end)) {
            let last = decimal_of(end);
            let first = decimal_of(&position);
            span = Some((last == first, spans(&first, &last, allele_length)));
        }
    }
    if let Some((one_base, as_long_as_the_allele)) = span {
        if !one_base && !(substitution && as_long_as_the_allele) {
            return None;
        }
    }
    (!chromosome.is_empty() && pyfmt::is_digit(&position))
        .then(|| (pyfmt::removeprefix(chromosome, "chr").to_string(), position))
}

/// Whether `pairs` holds the row's site and allele under one of the studies the
/// `GWAS_study` value `studies` names (joined by `,`; `&` in the VCF format).
pub fn gwas_pair_named(
    e: &Entry,
    studies: &str,
    site: (&str, &str),
    pairs: &BTreeSet<GwasPair>,
) -> bool {
    !blank(studies)
        && studies.split([',', '&']).any(|study| {
            pairs.contains(&(
                site.0.to_string(),
                site.1.to_string(),
                e.allele.clone(),
                study.to_string(),
            ))
        })
}

/// Whether the entry carries a `GWAS_` field with a value: the plugin annotated its row.
pub fn gwas_annotated(e: &Entry) -> bool {
    e.fields
        .iter()
        .any(|(name, value)| field_named(name, &[GWAS_PREFIX]) && !blank(value))
}
