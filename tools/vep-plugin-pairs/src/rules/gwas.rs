// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The GWAS plugin's classes and the catalog they read.
//!
//! Neither output says which catalog row annotated a site, so both classes read the
//! catalog's rows of two or more SNPs (`GwasCatalogPairs`) and accept nothing without them.
//!
//! `gwas_haplotype_cross_product`: GWAS.pm:338-340 hand every location of a haplotype or
//! interaction row to every SNP of the row, so a site carries every SNP's risk allele and
//! study, where the engine pairs SNP i with location i. The reference annotates a superset
//! of sites: the shape is a `GWAS_` field the reference fills and the engine leaves empty,
//! or one both fill differently, on a row at a location of such a catalog row whose allele
//! is another of its SNPs' risk alleles and whose reference `GWAS_study` names its study.
//!
//! `gwas_interaction_first_snp_allele_whitespace`: GWAS.pm:343-348 take each risk allele of
//! a `SNP x SNP` row as the split on `x` leaves it, so every allele but the last keeps the
//! space before the ` x ` and matches no variant, where the engine trims it. The shape is a
//! `GWAS_` field the engine alone fills on a row at the location of such an SNP, the row's
//! allele that SNP's risk allele and the engine's `GWAS_study` naming the catalog row's
//! study.

use super::entry::{blank, field_named, field_of, field_of_any, text_lines, Entry, Row};
use super::shared::{digits_value, read_text, record_entry};
use super::skipped_record::gvcf_site_rows;
use super::{Consult, Refusal};
use crate::pytext;
use regex::Regex;
use std::collections::{BTreeSet, HashSet};
use std::path::Path;
use std::sync::LazyLock;

const GWAS_FIELDS: [&str; 1] = ["gwas_"];
const STUDY_FIELDS: [&str; 1] = ["GWAS_study"];
/// The columns of an associations file the pairs are read from: the chromosomes, the
/// positions, the risk alleles and the study.
const CATALOG_COLUMNS: [&str; 4] = [
    "CHR_ID",
    "CHR_POS",
    "STRONGEST SNP-RISK ALLELE",
    "STUDY ACCESSION",
];
const SNPS_COLUMN: &str = "SNPS";
/// The characters the plugin splits a catalog cell's list at.
const LIST_SEPARATORS: [char; 3] = [';', ',', 'x'];

static RS_ID: LazyLock<Regex> = LazyLock::new(|| Regex::new("rs[0-9]+").expect("regex"));

/// A chromosome without `chr`, a position, a risk allele and a study accession.
pub type Pair = (String, String, String, String);

fn is_bases(text: &str, bases: &[u8]) -> bool {
    !text.is_empty() && text.bytes().all(|b| bases.contains(&b))
}

fn without_chr(chromosome: &str) -> &str {
    chromosome.strip_prefix("chr").unwrap_or(chromosome)
}

/// The places where the two engines read a GWAS catalog's rows of two or more SNPs
/// differently: `reference_only` where Ensembl VEP hands a location another SNP's risk
/// allele (GWAS.pm:338-340), `engine_only` where it keeps the space after an allele written
/// before ` x ` and so never matches the SNP's own (GWAS.pm:343-348).
#[derive(Debug, Default)]
pub struct GwasCatalogPairs {
    reference_only: HashSet<Pair>,
    engine_only: HashSet<Pair>,
}

impl GwasCatalogPairs {
    /// The pairs of the associations file `path`, read by its header's column names. A
    /// file whose first line lacks one of the four columns is refused by name.
    pub fn from_catalog(path: &Path) -> Result<GwasCatalogPairs, Refusal> {
        let text = read_text(path)?;
        let lines = text_lines(&text);
        let header: Vec<&str> = lines
            .first()
            .map(|line| line.split('\t').collect())
            .unwrap_or_default();
        let column = |name: &str| header.iter().position(|cell| *cell == name);
        let named: Option<Vec<usize>> = CATALOG_COLUMNS.iter().map(|name| column(name)).collect();
        let Some(columns) = named else {
            return Err(Refusal(format!(
                "{}: not a GWAS catalog associations file (columns {})",
                path.display(),
                CATALOG_COLUMNS.join(", ")
            )));
        };
        let snps = column(SNPS_COLUMN);
        let last = columns.iter().copied().max().unwrap_or(0);
        let mut pairs = GwasCatalogPairs::default();
        for line in lines.iter().skip(1) {
            let cells: Vec<&str> = line.split('\t').collect();
            if cells.len() > last {
                pairs.add_row(
                    cells[columns[0]],
                    cells[columns[1]],
                    cells[columns[2]],
                    cells[columns[3]],
                    snps.and_then(|at| cells.get(at).copied()).unwrap_or(""),
                );
            }
        }
        Ok(pairs)
    }

    /// One catalog row: its lists split as the plugin splits them, the locations stripped,
    /// each risk allele the text after its token's first `-` as written. A row of several
    /// locations whose SNPs, the `rs` ids of `snps` or else its risk alleles, differ from
    /// them in count is one the engine skips, while the reference writes every SNP at every
    /// location (GWAS.pm:237 compares the chromosomes with the positions alone), so each of
    /// its alleles is the reference's alone at each location. A row whose chromosomes and
    /// positions differ in count, or whose risk alleles differ in count from as many SNPs
    /// as locations, pairs nothing. Neither does a row whose risk-allele cell holds a `?`:
    /// both engines read that cell as empty (GWAS.pm:261) and write nothing from the row.
    pub fn add_row(
        &mut self,
        chromosomes: &str,
        positions: &str,
        risk_alleles: &str,
        study: &str,
        snps: &str,
    ) {
        let chromosome_list: Vec<&str> = chromosomes
            .split(LIST_SEPARATORS)
            .map(|chromosome| without_chr(pytext::strip(chromosome)))
            .collect();
        let position_list: Vec<&str> = positions
            .split(LIST_SEPARATORS)
            .map(pytext::strip)
            .collect();
        let written: Vec<&str> = risk_alleles
            .split(LIST_SEPARATORS)
            .map(|token| token.split('-').nth(1).unwrap_or(""))
            .collect();
        if position_list.len() < 2
            || chromosome_list.len() != position_list.len()
            || risk_alleles.contains('?')
        {
            return;
        }
        let matched_as_written: BTreeSet<&str> = written
            .iter()
            .copied()
            .filter(|allele| is_bases(allele, b"ACGT"))
            .collect();
        let pair = |chromosome: &str, position: &str, allele: &str| {
            (
                chromosome.to_string(),
                position.to_string(),
                allele.to_string(),
                study.to_string(),
            )
        };
        let places = || chromosome_list.iter().zip(&position_list);
        let named = match RS_ID.find_iter(snps).count() {
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
            let own = Some(pytext::strip(allele)).filter(|allele| is_bases(allele, b"ACGT"));
            for other in &matched_as_written {
                if own != Some(*other) {
                    self.reference_only
                        .insert(pair(chromosome, position, other));
                }
            }
            if let Some(own) = own.filter(|own| !matched_as_written.contains(own)) {
                self.engine_only.insert(pair(chromosome, position, own));
            }
        }
    }

    /// The places and alleles the reference alone annotates.
    pub fn reference_only(&self) -> &HashSet<Pair> {
        &self.reference_only
    }

    /// The places and alleles the engine alone annotates.
    pub fn engine_only(&self) -> &HashSet<Pair> {
        &self.engine_only
    }
}

/// Whether a span of `span` bases is one the entry's allele replaces base for base: a
/// sequence allele of as many bases that the record's own alleles do not mark an insertion
/// or a deletion.
fn replaces_its_span(e: &Entry, span: i128) -> bool {
    is_bases(&e.allele, b"ACGTN") && e.indel != Some(true) && span == e.allele.len() as i128
}

/// The chromosome, without `chr`, and position at which the plugin matches a row's allele
/// (GWAS.pm:176-195 compare the record's alleles at its start with a risk allele at the
/// catalog's position): a one-base Location's position, or the first position of a Location
/// spanning as many bases as the row's sequence allele, a substitution of several bases,
/// which an insertion the parser read from the record's alleles is not. An entry without a
/// Location is read from the record entry among `record_rows`: its `CHROM` with the start
/// the parser gives the row's allele from the line's `POS`, `REF` and `ALT` (the `POS`
/// itself for an allele it does not place), or its `seq_region_name` and `start`; such a row
/// is held to the Location's condition where its record states the span the allele
/// replaces, so an insertion's row, whose span is empty, states no site. `None` where the
/// row states none, and for a position `int()` does not read.
pub fn gwas_site(e: &Entry, record_rows: Option<&[Row]>) -> Option<(String, String)> {
    if let Some(location) = e.location.as_deref() {
        let (chromosome, _, span) = pytext::rpartition(location, ":");
        let (start, dash_between, end) = pytext::partition(span, "-");
        if chromosome.is_empty() || !pytext::isdigit(start) {
            return None;
        }
        if !dash_between.is_empty() {
            if !(is_bases(&e.allele, b"ACGTN") && e.indel != Some(true)) {
                return None;
            }
            let bases = (digits_value(end)? - digits_value(start)?).saturating_add(1);
            if !replaces_its_span(e, bases) {
                return None;
            }
        }
        return Some((without_chr(chromosome).to_string(), start.to_string()));
    }
    let record = record_entry(e, record_rows.unwrap_or(&[]))?;
    let chromosome = record.fields.first_of(&["CHROM", "seq_region_name"], "");
    let mut position = record.fields.first_of(&["POS", "start"], "").to_string();
    let mut span: Option<i128> = None;
    if record.fields.contains("POS") && pytext::isdigit(&position) {
        // Parser/VCF.pm:300-336 drop the first base an indel record's alleles share and move
        // its start one base on, so the row's start is not the line's POS there.
        let alternates: Vec<&str> = record.fields.get_or("ALT", "").split(',').collect();
        let allele = e.allele.to_uppercase();
        let placed: Vec<(i128, i128)> = gvcf_site_rows(
            pytext::int(&position)?,
            record.fields.get_or("REF", ""),
            &alternates,
        )
        .into_iter()
        .filter(|(written, _, _)| *written == allele)
        .map(|(_, start, end)| (start, end))
        .collect();
        let starts: BTreeSet<i128> = placed.iter().map(|(start, _)| *start).collect();
        if let (1, Some(start)) = (starts.len(), starts.first()) {
            position = start.to_string();
            let spans: BTreeSet<i128> = placed.iter().map(|(start, end)| end - start + 1).collect();
            span = spans.first().copied().filter(|_| spans.len() == 1);
        }
    } else if pytext::isdigit(&position) && pytext::isdigit(record.fields.get_or("end", "")) {
        let end = pytext::int(record.fields.get_or("end", ""))?;
        span = Some((end - pytext::int(&position)?).saturating_add(1));
    }
    if span.is_some_and(|span| span != 1 && !replaces_its_span(e, span)) {
        return None;
    }
    (!chromosome.is_empty() && pytext::isdigit(&position))
        .then(|| (without_chr(chromosome).to_string(), position))
}

/// Whether `pairs` holds the row's site and allele under one of the studies the
/// `GWAS_study` value `studies` names (joined by `,`, or by `&` in the VCF format).
fn pair_named(e: &Entry, studies: &str, site: (String, String), pairs: &HashSet<Pair>) -> bool {
    if blank(studies) {
        return false;
    }
    let mut pair = (site.0, site.1, e.allele.clone(), String::new());
    studies.split([',', '&']).any(|study| {
        pair.3.clear();
        pair.3.push_str(study);
        pairs.contains(&pair)
    })
}

/// Whether the entry carries a `GWAS_` field with a value: the plugin annotated its row.
fn annotated(e: &Entry) -> bool {
    e.fields
        .iter()
        .any(|(name, value)| field_named(name, &GWAS_FIELDS) && !blank(value))
}

/// The catalog a consult was handed and the site its reference row states under it.
fn catalog_site<'a>(c: &Consult<'a>) -> Option<(&'a GwasCatalogPairs, (String, String))> {
    let pairs = c.gwas_pairs?;
    Some((pairs, gwas_site(c.e, c.record_rows)?))
}

/// `gwas_haplotype_cross_product`: `ev` and `av` are the two values with a blank one as
/// `-`. A field the catalog row leaves empty is a defined empty string on an annotated
/// reference row (`GWAS_beta_coef=`) and no key on a row the engine never annotated: the
/// pair is the class's where the reference row fills another `GWAS_` field, and a key the
/// engine alone writes stays a difference. A value the engine alone writes is its own
/// catalog row's where both rows are annotated and it reports another study than the
/// reference, whose row leaves the field empty.
pub fn haplotype_cross_product(c: &Consult, ev: &str, av: &str) -> bool {
    let Some((pairs, site)) = catalog_site(c) else {
        return false;
    };
    let study = field_of(c.e, &STUDY_FIELDS);
    if !pair_named(c.e, study, site, &pairs.reference_only) {
        return false;
    }
    if !blank(ev) {
        return true;
    }
    if blank(av) {
        return c.expected.is_empty() && blank(c.actual) && annotated(c.e);
    }
    let own_study = field_of_any(c.a.map(|a| &**a), &STUDY_FIELDS);
    annotated(c.e) && !blank(own_study) && own_study != study
}

/// `gwas_interaction_first_snp_allele_whitespace`: `ev` and `av` are the two values with a
/// blank one as `-`. The plugin writes a field the catalog row leaves empty as a defined
/// empty string (`GWAS_odds_ratio=`), which the reference's unannotated row has no key for:
/// the class's beside the engine's filled fields.
pub fn interaction_first_snp_allele_whitespace(c: &Consult, ev: &str, av: &str) -> bool {
    let Some((pairs, site)) = catalog_site(c) else {
        return false;
    };
    if !blank(ev) || annotated(c.e) {
        return false;
    }
    if blank(av) && !(c.actual.is_empty() && c.a.is_some_and(|a| annotated(a))) {
        return false;
    }
    let study = field_of_any(c.a.map(|a| &**a), &STUDY_FIELDS);
    pair_named(c.e, study, site, &pairs.engine_only)
}
