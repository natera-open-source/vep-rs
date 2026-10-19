// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The allele frequency fields of a row: `AF`, the 1000 Genomes continental
//! and gnomAD population columns, `MAX_AF` and `MAX_AF_POPS`, filled as
//! `add_colocated_frequency_data` (OutputFactory.pm:1185-1280 at ensembl-vep
//! release/116.2) fills them from the co-located variants of an input record.
//!
//! A co-located variant stores each population's frequencies as
//! `allele:frequency` pairs. A row reads the frequency of the record allele
//! matched to its own allele; the `AF` column of a biallelic record whose one
//! stored allele is the other one is the remainder `1 - total`. The values of
//! every co-located variant are merged into one list per field, each value
//! once; `MAX_AF` is the largest frequency over the population columns of
//! every frequency group whatever flags are set, the combined populations left
//! out, and `MAX_AF_POPS` lists the populations that reach it.

use std::borrow::Cow;
use std::sync::LazyLock;

use vep_core::sv_routing::{perl_split, perl_true};
use vep_core::variant::{ColocatedVariant, InputVariant};
use vep_core::variation::{perl_number_string, perl_numeric};

/// The gnomAD exome columns the `af_gnomad` and `af_gnomade` groups read
/// (OutputFactory.pm:110-111).
const GNOMADE_COLUMNS: &[&str] = &[
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
];

/// `%FREQUENCY_KEYS` (OutputFactory.pm:107-113): each frequency flag with the
/// cache columns it reads, the groups in the order `sort keys %FREQUENCY_KEYS`
/// visits them and each group's columns in their declared order. A column in
/// two groups is read once per group and merged, so it prints once.
pub const FREQUENCY_KEYS: &[(&str, &[&str])] = &[
    ("af", &["AF"]),
    ("af_1kg", &["AF", "AFR", "AMR", "ASN", "EAS", "EUR", "SAS"]),
    ("af_gnomad", GNOMADE_COLUMNS),
    ("af_gnomade", GNOMADE_COLUMNS),
    (
        "af_gnomadg",
        &[
            "gnomADg",
            "gnomADg_AFR",
            "gnomADg_AMI",
            "gnomADg_AMR",
            "gnomADg_ASJ",
            "gnomADg_EAS",
            "gnomADg_FIN",
            "gnomADg_MID",
            "gnomADg_NFE",
            "gnomADg_REMAINING",
            "gnomADg_SAS",
        ],
    ),
];

/// The combined-population columns `MAX_AF` never takes (OutputFactory.pm:1252).
const COMBINED_POPULATIONS: &[&str] = &["AF", "ExAC", "ExAC_Adj", "gnomADe", "gnomADg"];

/// The population columns `MAX_AF` compares, in the order Perl visits them, each
/// once: Perl visits a column in two groups twice, and the second visit finds
/// the same value already at or below the maximum with its population already
/// listed.
static MAX_AF_POPULATIONS: LazyLock<Vec<&'static str>> = LazyLock::new(|| {
    let mut populations: Vec<&'static str> = Vec::new();
    for column in FREQUENCY_KEYS
        .iter()
        .flat_map(|(_, columns)| columns.iter())
    {
        if !COMBINED_POPULATIONS.contains(column) && !populations.contains(column) {
            populations.push(column);
        }
    }
    populations
});

/// Which frequency groups a run reports and whether it reports the maximum.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FrequencyGroups {
    pub af: bool,
    pub af_1kg: bool,
    pub af_gnomad: bool,
    pub af_gnomade: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
}

impl FrequencyGroups {
    /// Whether the group's fields are reported.
    fn reports(&self, group: &str) -> bool {
        match group {
            "af" => self.af,
            "af_1kg" => self.af_1kg,
            "af_gnomad" => self.af_gnomad,
            "af_gnomade" => self.af_gnomade,
            "af_gnomadg" => self.af_gnomadg,
            _ => false,
        }
    }

    /// Whether a co-located variant matched on none of the row's alleles is
    /// read at all: it is when the `af` group is among the groups scanned, which
    /// are every group under `max_af` and the reported ones otherwise
    /// (OutputFactory.pm:1193-1194, 1202), since its `AF` remainder can still
    /// name the row's allele.
    fn reads_unmatched(&self) -> bool {
        self.max_af || self.af
    }

    /// The groups the JSON factory reports: it switches the 1000 Genomes and
    /// gnomAD groups on in its constructor (OutputFactory/JSON.pm:160-163), so a
    /// `colocated_variants` entry nests those frequencies whatever the flags.
    pub fn for_json(self) -> Self {
        Self {
            af_1kg: true,
            af_gnomad: true,
            af_gnomade: true,
            af_gnomadg: true,
            ..self
        }
    }
}

/// The cache column behind a frequency output field: `AF` for `AF`, the column
/// before the suffix for `<column>_AF` (`AFR_AF` reads `AFR`, `gnomADe_NFE_AF`
/// reads `gnomADe_NFE`), `None` for a name no frequency group produces, which
/// a plugin may own.
pub fn frequency_column(field: &str) -> Option<&'static str> {
    let column = match field.strip_suffix("_AF") {
        None if field == "AF" => "AF",
        Some("AF") | None => return None,
        Some(column) => column,
    };
    FREQUENCY_KEYS
        .iter()
        .flat_map(|(_, columns)| columns.iter())
        .find(|candidate| **candidate == column)
        .copied()
}

/// Whether a reported group reads the column (the `$self->{$group}` test at
/// OutputFactory.pm:1245).
fn column_reported(column: &str, groups: &FrequencyGroups) -> bool {
    FREQUENCY_KEYS
        .iter()
        .any(|(group, columns)| groups.reports(group) && columns.contains(&column))
}

/// One frequency as the loop at OutputFactory.pm:1216-1231 leaves it: a stored
/// pair's value, borrowed from the record's text or formatted `%.4f` in the
/// `AF` column, or the remainder of a biallelic record's `AF`.
#[derive(Debug, Clone, PartialEq)]
pub enum Frequency<'a> {
    Stored(Cow<'a, str>),
    Remainder(f64),
}

impl Frequency<'_> {
    /// The value as Ensembl VEP prints it: the stored text as the cache wrote
    /// it, or the remainder as Perl stringifies a number (`%.15g`).
    pub fn text(&self) -> String {
        match self {
            Frequency::Stored(text) => text.to_string(),
            Frequency::Remainder(value) => perl_number_string(*value),
        }
    }

    /// The value as a Perl numeric context reads it (the `MAX_AF` comparisons,
    /// OutputFactory.pm:1254-1258).
    fn numeric(&self) -> f64 {
        match self {
            Frequency::Stored(text) => perl_numeric(text),
            Frequency::Remainder(value) => *value,
        }
    }

    /// Perl's truth of the value: a stored text other than the empty string
    /// and `0`, a remainder other than zero.
    fn perl_true(&self) -> bool {
        match self {
            Frequency::Stored(text) => perl_true(text),
            Frequency::Remainder(value) => *value != 0.0,
        }
    }
}

/// The frequency one co-located variant contributes to a reported column on the
/// row of `allele`, or `None` when the column is not reported, the record does
/// not carry it, or neither the matched record allele nor the row's allele has
/// a value in it (OutputFactory.pm:1200-1245). A record matched to none of the
/// row's alleles is read only when the `af` group is scanned (:1202).
pub fn reported_frequency<'a>(
    cv: &'a ColocatedVariant,
    column: &str,
    allele: &str,
    groups: &FrequencyGroups,
) -> Option<Frequency<'a>> {
    if !column_reported(column, groups) {
        return None;
    }
    let matched = cv.record_allele_for(allele);
    if matched.is_none() && !groups.reads_unmatched() {
        return None;
    }
    column_frequency(cv, matched, column, allele)
}

/// The frequency of one co-located variant's column on the row of `allele`,
/// whose record allele is `matched`: the value stored for the matched allele,
/// the last pair of a repeated allele winning (OutputFactory.pm:1220), else the
/// row allele's own entry when that is the remainder of a biallelic record's
/// `AF` (:1226-1243).
///
/// Every output field of every row reads its column through here, so the pairs
/// are scanned in place: no map of the record's alleles is built.
fn column_frequency<'a>(
    cv: &'a ColocatedVariant,
    matched: Option<&str>,
    column: &str,
    allele: &str,
) -> Option<Frequency<'a>> {
    let text = cv.frequency_strings.get(column).filter(|t| perl_true(t))?;
    let is_af = column == "AF";
    let pairs = || {
        text.split(',')
            .map(|pair| pair.split_once(':').unwrap_or((pair, "")))
    };
    let mut total = 0.0;
    let mut matched_value: Option<Cow<'a, str>> = None;
    let mut own_value: Option<Cow<'a, str>> = None;
    for (pair_allele, value) in pairs() {
        let value: Cow<'a, str> = if is_af {
            let value = Cow::Owned(format!("{:.4}", perl_numeric(value)));
            total += perl_numeric(&value);
            value
        } else {
            Cow::Borrowed(value)
        };
        if pair_allele == allele {
            own_value = Some(value.clone());
        }
        if matched == Some(pair_allele) {
            matched_value = Some(value);
        }
    }
    let remainder = if is_af {
        af_remainder(cv.allele_string.as_deref().unwrap_or(""), pairs, total)
    } else {
        None
    };
    if let Some(record_allele) = matched {
        if let Some(value) = matched_value {
            return Some(Frequency::Stored(value));
        }
        if let Some((remaining_allele, value)) = remainder {
            if remaining_allele == record_allele {
                return Some(Frequency::Remainder(value));
            }
        }
    }
    let (remaining_allele, value) = remainder?;
    let own = if remaining_allele == allele {
        Frequency::Remainder(value)
    } else {
        Frequency::Stored(own_value?)
    };
    Some(own).filter(Frequency::perl_true)
}

/// The interpolated `AF` of a biallelic record (OutputFactory.pm:1223-1231): the
/// one distinct record allele no pair names (`%remaining`, :1214) with
/// `1 - total`, or `None` when the record has other than two alleles or another
/// count of unnamed ones.
fn af_remainder<'a, P>(
    allele_string: &'a str,
    pairs: impl Fn() -> P,
    total: f64,
) -> Option<(&'a str, f64)>
where
    P: Iterator<Item = (&'a str, &'a str)>,
{
    let record_alleles = || perl_split(allele_string, &['/']);
    if record_alleles().count() != 2 {
        return None;
    }
    let named = |record_allele: &str| pairs().any(|(pair_allele, _)| pair_allele == record_allele);
    let mut unnamed = record_alleles().filter(|record_allele| !named(record_allele));
    let remaining = unnamed.next()?;
    if unnamed.any(|other| other != remaining) {
        return None;
    }
    Some((remaining, 1.0 - total))
}

/// The values of one frequency output field on a row: the reported frequency of
/// each co-located variant in the record's order, each value once
/// (`merge_arrays`, Utils.pm:303-310), as Ensembl VEP prints them.
pub fn frequency_field(
    column: &str,
    variant: &InputVariant,
    groups: &FrequencyGroups,
) -> Vec<String> {
    let allele = variant.display_allele();
    let mut values: Vec<String> = Vec::new();
    for cv in &variant.colocated_variants {
        if let Some(frequency) = reported_frequency(cv, column, &allele, groups) {
            let text = frequency.text();
            if !values.contains(&text) {
                values.push(text);
            }
        }
    }
    values
}

/// `MAX_AF` and `MAX_AF_POPS` of a row (OutputFactory.pm:1197-1198, 1250-1275),
/// or `None` when no co-located variant has a population frequency to compare.
/// Within a variant, a frequency above the running maximum replaces it and
/// names its population alone, an equal one adds its population once; the
/// combined populations never count. Across variants the row's maximum starts
/// at `0`, a larger maximum replaces it and its populations, and an equal one
/// appends its populations as they are. The maximum prints as the cache wrote
/// it, `0` when no population exceeds zero.
pub fn max_af(
    variant: &InputVariant,
    groups: &FrequencyGroups,
) -> Option<(String, Vec<&'static str>)> {
    if !groups.max_af {
        return None;
    }
    let allele = variant.display_allele();
    let mut row_max: Option<(f64, String)> = None;
    let mut row_pops: Vec<&'static str> = Vec::new();
    for cv in &variant.colocated_variants {
        let matched = cv.record_allele_for(&allele);
        if matched.is_none() && !groups.reads_unmatched() {
            continue;
        }
        let mut max: f64 = 0.0;
        let mut max_text = "0".to_string();
        let mut pops: Vec<&'static str> = Vec::new();
        for &column in MAX_AF_POPULATIONS.iter() {
            let Some(frequency) = column_frequency(cv, matched, column, &allele) else {
                continue;
            };
            let value = frequency.numeric();
            if value > max {
                max = value;
                max_text = frequency.text();
                pops = vec![column];
            } else if value == max && !pops.contains(&column) {
                pops.push(column);
            }
        }
        if pops.is_empty() {
            continue;
        }
        match row_max {
            Some((current, _)) if max < current => {}
            Some((current, _)) if max == current => row_pops.extend(pops),
            _ => {
                row_max = Some((max, max_text));
                row_pops = pops;
            }
        }
    }
    row_max.map(|(_, text)| (text, row_pops))
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::MatchedAllele;

    /// A known variant with the allele string and `frequency_strings` given,
    /// matched to the input alleles `matched` (input allele, record allele).
    fn known(
        allele_string: &str,
        matched: &[(&str, &str)],
        columns: &[(&str, &str)],
    ) -> ColocatedVariant {
        ColocatedVariant {
            id: "rs1".to_string(),
            allele_string: Some(allele_string.to_string()),
            strand: 1,
            matched_alleles: matched
                .iter()
                .enumerate()
                .map(|(i, (input, record))| MatchedAllele {
                    input_allele: input.to_string(),
                    input_index: i,
                    record_allele: record.to_string(),
                    record_index: i,
                })
                .collect(),
            frequency_strings: columns
                .iter()
                .map(|(column, value)| (column.to_string(), value.to_string()))
                .collect(),
            ..Default::default()
        }
    }

    /// One allele of an input record carrying the co-located list `known`.
    fn row(reference: &str, alt: &str, known: Vec<ColocatedVariant>) -> InputVariant {
        let mut variant = InputVariant::new(
            "21".into(),
            5_131_441,
            5_131_441,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant.colocated_variants = known;
        variant
    }

    fn reporting(set: impl FnOnce(&mut FrequencyGroups)) -> FrequencyGroups {
        let mut groups = FrequencyGroups::default();
        set(&mut groups);
        groups
    }

    fn field(column: &str, variant: &InputVariant, groups: &FrequencyGroups) -> String {
        frequency_field(column, variant, groups).join(",")
    }

    /// The cache columns of the frequency group `group`.
    fn group_columns(group: &str) -> &'static [&'static str] {
        let (_, columns) = FREQUENCY_KEYS
            .iter()
            .find(|(candidate, _)| *candidate == group)
            .unwrap();
        columns
    }

    /// `AF` prints with four decimals (`sprintf("%.4f")`, OutputFactory.pm:1219),
    /// the other columns as the cache wrote them: Ensembl VEP 116.2 writes
    /// `AF=0.0068;AFR_AF=0.0212` for `21 5131441 T C` of the `GRCh38-colocated`
    /// corpus, whose `rs1555830109` stores `C:0.006789` and `C:0.0212`, and
    /// `gnomADg_AFR_AF=6.875e-05` for `21 5030602 A AA` (`rs1370167485`,
    /// `A:6.875e-05`). A record matched by another allele prints nothing on this
    /// row.
    #[test]
    fn af_prints_four_decimals_for_matched_allele() {
        let groups = reporting(|g| {
            g.af = true;
            g.af_1kg = true;
            g.af_gnomadg = true;
        });
        let rs1555830109 = known(
            "T/C",
            &[("C", "C")],
            &[
                ("AF", "C:0.006789"),
                ("AFR", "C:0.0212"),
                ("gnomADg_AFR", "C:6.875e-05"),
            ],
        );
        let c = row("T", "C", vec![rs1555830109.clone()]);
        assert_eq!(field("AF", &c, &groups), "0.0068");
        assert_eq!(field("AFR", &c, &groups), "0.0212");
        assert_eq!(field("gnomADg_AFR", &c, &groups), "6.875e-05");
        let g = row("T", "G", vec![rs1555830109]);
        assert_eq!(field("AF", &g, &groups), "");
        assert_eq!(field("AFR", &g, &groups), "");
        let one = row(
            "T",
            "C",
            vec![known("T/C", &[("C", "C")], &[("AF", "C:1")])],
        );
        assert_eq!(field("AF", &one, &groups), "1.0000");
    }

    /// A biallelic record whose `AF` stores the other allele yields the remainder
    /// `1 - total` for the unstored one, printed as Perl prints a number
    /// (OutputFactory.pm:1226-1231, 1243); the remainder is read on the row's
    /// own allele, so it reaches a record matched on its span alone, and only
    /// when the `af` group is scanned (:1202). A record with three alleles or a
    /// stored value for both is never interpolated.
    #[test]
    fn af_interpolates_remaining_allele_of_biallelic_record() {
        let groups = reporting(|g| g.af = true);
        let span_only = known("A/G", &[], &[("AF", "A:0.3")]);
        let g = row("A", "G", vec![span_only.clone()]);
        assert_eq!(field("AF", &g, &groups), "0.7");
        let a = row("G", "A", vec![span_only.clone()]);
        assert_eq!(field("AF", &a, &groups), "0.3000");
        assert_eq!(field("AF", &g, &reporting(|g| g.af_1kg = true)), "");
        assert_eq!(
            field(
                "AF",
                &g,
                &reporting(|g| {
                    g.af_1kg = true;
                    g.max_af = true;
                })
            ),
            "0.7"
        );
        let none = known("A/G", &[], &[("AF", "A:1")]);
        assert_eq!(field("AF", &row("A", "G", vec![none]), &groups), "");
        let tri = known("A/G/T", &[], &[("AF", "A:0.3")]);
        assert_eq!(field("AF", &row("A", "G", vec![tri]), &groups), "");
        let both = known("A/G", &[], &[("AF", "A:0.3,G:0.6")]);
        assert_eq!(
            field("AF", &row("A", "G", vec![both]), &groups),
            "",
            "a span-only match with no remainder reports nothing (:1233-1236)"
        );
        let rounded = known("T/C", &[], &[("AF", "C:0.006789")]);
        assert_eq!(
            field("AF", &row("C", "T", vec![rounded]), &groups),
            "0.9932"
        );
    }

    /// The continental columns follow the matched record allele, the last pair of
    /// a repeated allele winning (OutputFactory.pm:1220: GRCh37 `rs538160740`
    /// stores `T:0,T:0.005244` in `gnomADg_AFR`), and the values of every
    /// co-located variant merge in the record's order, each once (:1247,
    /// `merge_arrays`).
    #[test]
    fn continental_frequencies_follow_matched_allele_and_merge_across_variants() {
        let groups = reporting(|g| {
            g.af_1kg = true;
            g.af_gnomadg = true;
        });
        let first = known(
            "A/G/T",
            &[("G", "G"), ("T", "T")],
            &[("EAS", "G:0.3,T:0.1"), ("gnomADg_AFR", "T:0,T:0.005244")],
        );
        let second = known("A/G", &[("G", "G")], &[("EAS", "G:0.3")]);
        let third = known("A/G", &[("G", "G")], &[("EAS", "G:0.25")]);
        let g = row("A", "G", vec![first.clone(), second, third]);
        assert_eq!(field("EAS", &g, &groups), "0.3,0.25");
        let t = row("A", "T", vec![first]);
        assert_eq!(field("EAS", &t, &groups), "0.1");
        assert_eq!(field("gnomADg_AFR", &t, &groups), "0.005244");
        assert_eq!(field("EAS", &t, &reporting(|g| g.af = true)), "");
    }

    /// The gnomAD exome fields carry Ensembl VEP's names, `gnomADe_<POP>_AF`
    /// over the nine populations and `gnomADe_AF` for the combined one
    /// (Constants.pm:120-121), under `--af_gnomade` and its alias `--af_gnomad`
    /// alike (OutputFactory.pm:110-111).
    #[test]
    fn gnomade_columns_named_as_ensembl_vep() {
        let fields = [
            "gnomADe_AF",
            "gnomADe_AFR_AF",
            "gnomADe_AMR_AF",
            "gnomADe_ASJ_AF",
            "gnomADe_EAS_AF",
            "gnomADe_FIN_AF",
            "gnomADe_MID_AF",
            "gnomADe_NFE_AF",
            "gnomADe_REMAINING_AF",
            "gnomADe_SAS_AF",
        ];
        for (field, column) in fields.iter().zip(GNOMADE_COLUMNS) {
            assert_eq!(frequency_column(field), Some(*column), "{field}");
        }
        assert_eq!(frequency_column("AF"), Some("AF"));
        assert_eq!(frequency_column("MAX_AF"), None);
        assert_eq!(frequency_column("AF_AF"), None, "no such field (:1246)");
        assert_eq!(frequency_column("gnomAD_AF"), None, "a plugin's key");
        let record = known("T/A", &[("A", "A")], &[("gnomADe_NFE", "A:4.108e-05")]);
        let a = row("T", "A", vec![record]);
        assert_eq!(
            field("gnomADe_NFE", &a, &reporting(|g| g.af_gnomad = true)),
            "4.108e-05"
        );
        assert_eq!(
            field("gnomADe_NFE", &a, &reporting(|g| g.af_gnomade = true)),
            "4.108e-05"
        );
        assert_eq!(
            field("gnomADe_NFE", &a, &reporting(|g| g.af_gnomadg = true)),
            ""
        );
    }

    /// The gnomAD genome group reads the Amish population (`gnomADg_AMI`,
    /// OutputFactory.pm:112) beside the nine shared with the exomes.
    #[test]
    fn gnomadg_columns_include_ami() {
        let columns = group_columns("af_gnomadg");
        assert!(columns.contains(&"gnomADg_AMI"));
        assert!(!GNOMADE_COLUMNS.contains(&"gnomADe_AMI"));
        assert_eq!(columns.len(), 11);
        let record = known(
            "A/G",
            &[("G", "G")],
            &[("gnomADg_AMI", "G:1"), ("gnomADg", "G:0.9994")],
        );
        let g = row("A", "G", vec![record]);
        assert_eq!(
            field("gnomADg_AMI", &g, &reporting(|g| g.af_gnomadg = true)),
            "1"
        );
        assert_eq!(
            field("gnomADg", &g, &reporting(|g| g.af_gnomadg = true)),
            "0.9994"
        );
    }

    /// Populations tied at the maximum are all listed, in the order
    /// `sort keys %FREQUENCY_KEYS` then each group's declared order
    /// (OutputFactory.pm:1210-1211, 1255-1260): Ensembl VEP 116.2 writes
    /// `MAX_AF=1;MAX_AF_POPS=AFR,AMR,EAS,EUR,SAS,gnomADg_AFR,gnomADg_AMI,gnomADg_ASJ,gnomADg_FIN,gnomADg_MID,gnomADg_REMAINING,gnomADg_SAS`
    /// for the `G` allele of `21 7759128 T A,G` of the `GRCh38-colocated` corpus
    /// (`rs1555870136`), and across variants an equal maximum appends its
    /// populations without deduplication while a larger one replaces them
    /// (:1265-1275).
    #[test]
    fn max_af_ties_list_every_population_in_perl_key_order() {
        let groups = reporting(|g| g.max_af = true);
        let rs1555870136 = known(
            "T/A/G",
            &[("A", "A"), ("G", "G")],
            &[
                ("AF", "G:1"),
                ("AFR", "G:1"),
                ("AMR", "G:1"),
                ("EAS", "G:1"),
                ("EUR", "G:1"),
                ("SAS", "G:1"),
                ("gnomADg", "A:0.0002047,G:0.9994"),
                ("gnomADg_AFR", "A:0,G:1"),
                ("gnomADg_AMI", "A:0,G:1"),
                ("gnomADg_AMR", "A:0.003356,G:0.9967"),
                ("gnomADg_ASJ", "A:0,G:1"),
                ("gnomADg_EAS", "A:0,G:0.9907"),
                ("gnomADg_FIN", "A:0,G:1"),
                ("gnomADg_MID", "A:0,G:1"),
                ("gnomADg_NFE", "A:0,G:0.9996"),
                ("gnomADg_REMAINING", "A:0,G:1"),
                ("gnomADg_SAS", "A:0,G:1"),
            ],
        );
        let g = row("T", "G", vec![rs1555870136.clone()]);
        assert_eq!(
            max_af(&g, &groups),
            Some((
                "1".to_string(),
                vec![
                    "AFR",
                    "AMR",
                    "EAS",
                    "EUR",
                    "SAS",
                    "gnomADg_AFR",
                    "gnomADg_AMI",
                    "gnomADg_ASJ",
                    "gnomADg_FIN",
                    "gnomADg_MID",
                    "gnomADg_REMAINING",
                    "gnomADg_SAS",
                ]
            ))
        );
        let a = row("T", "A", vec![rs1555870136.clone()]);
        assert_eq!(
            max_af(&a, &groups),
            Some(("0.003356".to_string(), vec!["gnomADg_AMR"]))
        );
        let tied = known(
            "T/G",
            &[("G", "G")],
            &[("gnomADe_NFE", "G:1"), ("AFR", "G:1")],
        );
        let g2 = row("T", "G", vec![rs1555870136.clone(), tied.clone()]);
        let (_, pops) = max_af(&g2, &groups).unwrap();
        assert_eq!(&pops[12..], ["AFR", "gnomADe_NFE"]);
        assert_eq!(pops.len(), 14);
        let lower = known("T/A", &[("A", "A")], &[("gnomADe_NFE", "A:0.001")]);
        let a2 = row("T", "A", vec![lower, rs1555870136]);
        assert_eq!(
            max_af(&a2, &groups),
            Some(("0.003356".to_string(), vec!["gnomADg_AMR"]))
        );
        assert_eq!(max_af(&a2, &FrequencyGroups::default()), None);
    }

    /// A variant whose population frequencies are all zero sets `MAX_AF=0` and
    /// lists every zero population (OutputFactory.pm:1197, 1257-1260, 1265-1275):
    /// Ensembl VEP 116.2 writes `MAX_AF=0;MAX_AF_POPS=gnomADg_AFR,gnomADg_AMI,gnomADg_AMR,gnomADg_ASJ,gnomADg_EAS,gnomADg_FIN,gnomADg_MID,gnomADg_NFE,gnomADg_REMAINING,gnomADg_SAS`
    /// for `21 5030192 A G` of the `GRCh38-colocated` corpus (`rs1172282560`); a
    /// record with no population column at all sets neither key.
    #[test]
    fn max_af_zero_lists_zero_populations() {
        let groups = reporting(|g| g.max_af = true);
        let columns: Vec<(&str, &str)> = group_columns("af_gnomadg")
            .iter()
            .map(|column| (*column, "G:0"))
            .collect();
        let rs1172282560 = known("A/G", &[("G", "G")], &columns);
        let g = row("A", "G", vec![rs1172282560]);
        assert_eq!(
            max_af(&g, &groups),
            Some((
                "0".to_string(),
                vec![
                    "gnomADg_AFR",
                    "gnomADg_AMI",
                    "gnomADg_AMR",
                    "gnomADg_ASJ",
                    "gnomADg_EAS",
                    "gnomADg_FIN",
                    "gnomADg_MID",
                    "gnomADg_NFE",
                    "gnomADg_REMAINING",
                    "gnomADg_SAS",
                ]
            ))
        );
        let bare = row(
            "A",
            "G",
            vec![known("A/G", &[("G", "G")], &[("AF", "G:0")])],
        );
        assert_eq!(max_af(&bare, &groups), None);
    }

    /// `MAX_AF` scans every group whether requested or not and skips the
    /// combined populations `AF`, `ExAC`, `ExAC_Adj`, `gnomADe` and `gnomADg`
    /// (OutputFactory.pm:1193, 1252): Ensembl VEP 116.2 writes
    /// `MAX_AF=0.0003339;MAX_AF_POPS=gnomADg_NFE` for `21 5030602 A AA` of the
    /// `GRCh38-colocated` corpus (`rs1370167485`, `gnomADg` `0.0001735`,
    /// `gnomADg_NFE` `0.0003339`) and `MAX_AF=0.0212;MAX_AF_POPS=AFR` for
    /// `21 5131441 T C` (`AF` `0.006789` below `AFR` `0.0212`).
    #[test]
    fn max_af_excludes_combined_population_keys() {
        let groups = reporting(|g| g.max_af = true);
        let rs1370167485 = known(
            "-/A",
            &[("A", "A")],
            &[
                ("gnomADg", "A:0.0001735"),
                ("gnomADg_AFR", "A:6.875e-05"),
                ("gnomADg_NFE", "A:0.0003339"),
                ("gnomADg_AMI", "A:0"),
            ],
        );
        let a = row("-", "A", vec![rs1370167485]);
        assert_eq!(
            max_af(&a, &groups),
            Some(("0.0003339".to_string(), vec!["gnomADg_NFE"]))
        );
        assert_eq!(field("gnomADg_NFE", &a, &groups), "", "no group reported");
        let rs1555830109 = known(
            "T/C",
            &[("C", "C")],
            &[
                ("AF", "C:0.006789"),
                ("AFR", "C:0.0212"),
                ("gnomADe", "C:0.5"),
                ("gnomADg", "C:0.4"),
            ],
        );
        let c = row("T", "C", vec![rs1555830109]);
        assert_eq!(
            max_af(&c, &groups),
            Some(("0.0212".to_string(), vec!["AFR"]))
        );
        let combined_only = row(
            "T",
            "C",
            vec![known(
                "T/C",
                &[("C", "C")],
                &[("AF", "C:0.5"), ("gnomADe", "C:0.5")],
            )],
        );
        assert_eq!(max_af(&combined_only, &groups), None);
    }

    /// Under `--old_maf` the reader rewrites `AFR`, `AMR` and `EUR` to a bare
    /// minor-allele frequency (`BaseCacheVariation.pm:131-137`), which the
    /// pair split reads as an allele with no value (OutputFactory.pm:1217-1219),
    /// so those three fields print empty while `AF`, `EAS` and `SAS`, whose
    /// values keep their alleles, print: Ensembl VEP 116.2 writes
    /// `AF=0.5625;EAS_AF=0.4534;SAS_AF=0.5583` and no `AFR_AF`, `AMR_AF` or
    /// `EUR_AF` for `21 44195140 G A` of the `GRCh38-oldmaf` corpus
    /// (`rs2838554`, `AFR` `0.2852` after the rewrite).
    #[test]
    fn old_maf_blanks_afr_amr_eur_on_a_116_cache() {
        let groups = reporting(|g| g.af_1kg = true);
        let rs2838554 = known(
            "G/A",
            &[("A", "A")],
            &[
                ("AF", "A:0.5625"),
                ("AFR", "0.2852"),
                ("AMR", "0.379"),
                ("EAS", "A:0.4534"),
                ("EUR", "0.3976"),
                ("SAS", "A:0.5583"),
            ],
        );
        let a = row("G", "A", vec![rs2838554]);
        assert_eq!(field("AF", &a, &groups), "0.5625");
        assert_eq!(field("AFR", &a, &groups), "");
        assert_eq!(field("AMR", &a, &groups), "");
        assert_eq!(field("EUR", &a, &groups), "");
        assert_eq!(field("EAS", &a, &groups), "0.4534");
        assert_eq!(field("SAS", &a, &groups), "0.5583");
        assert_eq!(
            max_af(&a, &reporting(|g| g.max_af = true)),
            Some(("0.5583".to_string(), vec!["SAS"]))
        );
    }

    /// The JSON factory's groups: the 1000 Genomes and gnomAD groups on whatever
    /// the flags, `af` and `max_af` as given (OutputFactory/JSON.pm:160-163).
    #[test]
    fn json_groups_force_the_population_groups_on() {
        let json = reporting(|g| g.max_af = true).for_json();
        assert_eq!(
            json,
            FrequencyGroups {
                af: false,
                af_1kg: true,
                af_gnomad: true,
                af_gnomade: true,
                af_gnomadg: true,
                max_af: true,
            }
        );
        assert!(FrequencyGroups::default().for_json().af_1kg);
        assert!(!FrequencyGroups::default().for_json().af);
    }
}
