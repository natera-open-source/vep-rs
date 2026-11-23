// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classes a row or a record is filed under: the variant class of a row, the family
//! classes of `--family-columns`, the genotype classes of the kept samples and the record
//! classes of `--record-classes`.

use crate::constants::{GENOTYPE_CLASSES, RECORD_CLASS_KEYS, REGULATORY_CLASS, SEQUENCE_CLASSES};
use crate::error::{failure, Failure};
use crate::manifest::without_repeats;
use crate::names::{first, is_symbolic_or_breakend, last, length, or_dash, rest, without_last};
use crate::pytext::py_repr_str;
use crate::rows::{parse_extra, Extra, Row};
use regex::Regex;
use std::collections::{BTreeSet, HashSet};
use std::sync::LazyLock;

/// Whether a text is one ASCII digit or more.
fn is_ascii_digits(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit())
}

/// Whether the number `b` is the number `a` plus one, both ASCII digits of any length.
fn is_successor(a: &str, b: &str) -> bool {
    let a = a.trim_start_matches('0').as_bytes();
    let b = b.trim_start_matches('0').as_bytes();
    let (mut i, mut j, mut carry) = (a.len(), b.len(), 1);
    while i > 0 || carry > 0 {
        let mut digit = carry;
        if i > 0 {
            i -= 1;
            digit += a[i] - b'0';
        }
        carry = digit / 10;
        if j == 0 || b[j - 1] - b'0' != digit % 10 {
            return false;
        }
        j -= 1;
    }
    j == 0
}

/// A bracket allele without its mate parts: every stretch from one square bracket to the
/// next is dropped.
fn outside_brackets(allele: &str) -> String {
    let is_bracket = |c: char| c == '[' || c == ']';
    let mut kept = String::new();
    let mut tail = allele;
    while let Some(open) = tail.find(is_bracket) {
        let Some(close) = tail[open + 1..].find(is_bracket) else {
            break;
        };
        kept.push_str(&tail[..open]);
        tail = &tail[open + close + 2..];
    }
    kept.push_str(tail);
    kept
}

/// The variant shape of one row, from its Location and Allele columns: one of
/// `SEQUENCE_CLASSES` for a sequence allele; `breakend`, or `breakend_insertion` when the
/// bracket allele carries more than the anchor base outside the brackets; `single_breakend`;
/// and `symbolic_<word>` for the class word Ensembl VEP prints for a symbolic allele
/// (`symbolic_deletion`, `symbolic_Alu_insertion`). Each is annotated by a path of its own,
/// so each is a stratum of its own. A two-base substitution and a two-base insertion print
/// the same Location and Allele (Ensembl VEP writes an insertion between its two flanking
/// bases), so a span of two adjacent positions is `mnv_or_complex` when `substitution` says
/// the input record holds the allele as a substitution (`holds_allele`), and `insertion`
/// otherwise. A coordinate is a number when it is ASCII digits.
pub fn variant_class(location: &str, allele: &str, substitution: bool) -> String {
    if allele.contains('[') || allele.contains(']') {
        let class = if length(&outside_brackets(allele)) > 1 {
            "breakend_insertion"
        } else {
            "breakend"
        };
        return class.to_string();
    }
    if allele.starts_with('.') || allele.ends_with('.') {
        return "single_breakend".to_string();
    }
    if !allele.is_empty() && !allele.bytes().all(|byte| b"ACGTNacgtn*-".contains(&byte)) {
        return format!("symbolic_{allele}");
    }
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let class = if allele == "-" {
        "deletion"
    } else if let Some((a, b)) = span.split_once('-') {
        let adjacent = is_ascii_digits(a) && is_ascii_digits(b) && is_successor(a, b);
        if adjacent && !substitution {
            "insertion"
        } else {
            "mnv_or_complex"
        }
    } else if length(allele) == 1 {
        "snv"
    } else {
        "mnv_or_complex"
    };
    class.to_string()
}

/// The family class of a `custom_match=<exact>/<overlap>` spec, and its classes.
pub const CUSTOM_MATCH: &str = "custom_match";
pub const CUSTOM_MATCH_CLASSES: [&str; 5] =
    ["exact", "multi_allelic", "overlap", "sv_span", "miss"];

/// The classes of an `<a>/<b>` spec (`agreement_class`).
pub const AGREEMENT_CLASSES: [&str; 3] = ["same", "differ", "-"];

/// Whether a text is one character or more of `[A-Za-z0-9_.-]`.
fn is_column_name(text: &str) -> bool {
    !text.is_empty()
        && text
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"_.-".contains(&byte))
}

/// The two column names of `<a>/<b>`.
fn column_pair(text: &str) -> Option<(&str, &str)> {
    text.split_once('/')
        .filter(|(a, b)| is_column_name(a) && is_column_name(b))
}

/// One spec of `--family-columns`, as the class of a row is read under it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FamilySpec {
    /// `custom_match=<exact>/<overlap>`: how one `--custom` file answered through two sources.
    CustomMatch { exact: String, overlap: String },
    /// `<a>/<b>`: whether two columns agree.
    Agreement { spec: String, a: String, b: String },
    /// `NAME=`: the column's value.
    Value { name: String },
    /// `NAME`: whether the column is present.
    Presence { name: String },
}

impl FamilySpec {
    /// The spec a text states; any text states one, a text of no other form the presence
    /// of a column of its name.
    pub fn read(spec: &str) -> FamilySpec {
        let custom = spec
            .strip_prefix(CUSTOM_MATCH)
            .and_then(|tail| tail.strip_prefix('='))
            .and_then(column_pair);
        if let Some((exact, overlap)) = custom {
            FamilySpec::CustomMatch {
                exact: exact.to_string(),
                overlap: overlap.to_string(),
            }
        } else if let Some((a, b)) = column_pair(spec) {
            FamilySpec::Agreement {
                spec: spec.to_string(),
                a: a.to_string(),
                b: b.to_string(),
            }
        } else if let Some(name) = spec.strip_suffix('=') {
            FamilySpec::Value {
                name: name.to_string(),
            }
        } else {
            FamilySpec::Presence {
                name: spec.to_string(),
            }
        }
    }

    /// The class of a row under this spec: `custom_match=<class>`, `<a>/<b>=<class>`,
    /// `NAME=<value>` (`-` for an absent or empty value), or `NAME+` and `NAME-`.
    pub fn class(&self, extra: &Extra<'_>, location: &str, allele: &str) -> String {
        match self {
            FamilySpec::CustomMatch { exact, overlap } => {
                let class = custom_match_class(extra, exact, overlap, location, allele);
                format!("{CUSTOM_MATCH}={class}")
            }
            FamilySpec::Agreement { spec, a, b } => {
                format!("{spec}={}", agreement_class(extra, a, b))
            }
            FamilySpec::Value { name } => {
                let value = extra.get(name).filter(|value| !value.is_empty());
                format!("{name}={}", value.unwrap_or("-"))
            }
            FamilySpec::Presence { name } => {
                format!("{name}{}", if present(extra, name) { '+' } else { '-' })
            }
        }
    }
}

/// `--family-columns HGVSc,IMPACT=,custom_match=CE/CO,GIVEN_REF/USED_REF` as its specs, the
/// blanks around each dropped: a column name (presence), a name ending in `=` (value),
/// `custom_match=<exact>/<overlap>` (`custom_match_class`) or `<a>/<b>`, two columns classed
/// by their agreement (`agreement_class`). A spec of another form is refused.
pub fn parse_family_columns(arg: Option<&str>) -> Result<Vec<String>, Failure> {
    let specs: Vec<String> = arg
        .unwrap_or("")
        .split(',')
        .map(crate::pyfmt::py_strip)
        .filter(|spec| !spec.is_empty())
        .map(str::to_string)
        .collect();
    for spec in &specs {
        let named = is_column_name(spec.strip_suffix('=').unwrap_or(spec));
        let paired = matches!(
            FamilySpec::read(spec),
            FamilySpec::CustomMatch { .. } | FamilySpec::Agreement { .. }
        );
        if !named && !paired {
            return Err(failure(format!(
                "--family-columns wants NAME, NAME=, custom_match=EXACT/OVERLAP or A/B, got {}",
                py_repr_str(spec)
            )));
        }
    }
    Ok(specs)
}

/// Whether a row carries a value under `name`: a value that is neither empty nor `-`.
fn present(extra: &Extra<'_>, name: &str) -> bool {
    !matches!(extra.get(name), None | Some("" | "-"))
}

/// How two columns of one row relate: `same` when both are present with one value, `differ`
/// when both are present with two (`GIVEN_REF` and `USED_REF` on a row whose transcript
/// sequence differs from the input's reference allele), `-` when either is absent.
pub fn agreement_class(extra: &Extra<'_>, a: &str, b: &str) -> &'static str {
    if !(present(extra, a) && present(extra, b)) {
        "-"
    } else if extra.get(a) == extra.get(b) {
        "same"
    } else {
        "differ"
    }
}

/// The count of items of a comma-joined value.
fn items(value: &str) -> usize {
    value.matches(',').count() + 1
}

/// How one `--custom` file answered for a row through two sources over it: `exact` and
/// `overlap` are the short names of its `type=exact` and `type=overlap` sources, the
/// columns their record names are written under (`<name>_<field>` for a requested field,
/// one value per matched record, comma-joined). `exact`: the exact source matched.
/// `multi_allelic`: it matched and a field of the overlap source lists more values than the
/// source lists records while the exact source lists fewer, so a record carries one value
/// per allele and the exact match kept the matching allele's (`FILTER` lists a record's
/// filters, never alleles, and is passed over). `overlap`: the overlap source alone matched
/// a sequence variant. `sv_span`: it alone matched a structural variant (a symbolic or
/// breakend allele), whose span covers the record. `miss`: neither matched.
pub fn custom_match_class(
    extra: &Extra<'_>,
    exact: &str,
    overlap: &str,
    location: &str,
    allele: &str,
) -> &'static str {
    if present(extra, exact) {
        let records = if present(extra, overlap) {
            extra.get(overlap).map_or(1, items)
        } else {
            1
        };
        let prefix = format!("{overlap}_");
        for (key, value) in extra.iter() {
            let Some(field) = key.strip_prefix(&prefix) else {
                continue;
            };
            if field == "FILTER" {
                continue;
            }
            let overlap_values = items(value);
            let exact_values = items(extra.get(&format!("{exact}_{field}")).unwrap_or(""));
            if overlap_values > records && exact_values < overlap_values {
                return "multi_allelic";
            }
        }
        return "exact";
    }
    if present(extra, overlap) {
        let class = variant_class(location, allele, false);
        return if SEQUENCE_CLASSES.contains(&class.as_str()) {
            "overlap"
        } else {
            "sv_span"
        };
    }
    "miss"
}

/// One class per spec of `--family-columns` for a row with the Extra fields `extra`:
/// `NAME+` or `NAME-` for presence, `NAME=<value>` for a value (`-` when absent),
/// `custom_match=<class>` for a `custom_match=<exact>/<overlap>` spec, `<a>/<b>=<class>` for
/// an agreement spec.
pub fn family_classes<S: AsRef<str>>(
    extra: &Extra<'_>,
    specs: &[S],
    location: &str,
    allele: &str,
) -> Vec<String> {
    specs
        .iter()
        .map(|spec| FamilySpec::read(spec.as_ref()).class(extra, location, allele))
        .collect()
}

/// `family_classes` of a row whose Extra column is `raw_extra`, under specs read once.
pub fn row_family_classes(
    raw_extra: &str,
    specs: &[FamilySpec],
    location: &str,
    allele: &str,
) -> Vec<String> {
    let extra = parse_extra(raw_extra);
    specs
        .iter()
        .map(|spec| spec.class(&extra, location, allele))
        .collect()
}

/// The alleles Ensembl VEP reads as structural: a symbolic or bracket-breakend allele that
/// holds no `*`, a single breakend at the start or at the end of the comma-joined ALT.
/// `\w` is Python's: a letter, a number or `_`; `$` also stands before a last newline.
static STRUCTURAL_ALLELE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"[<\[\]][^*]+[>\]\[]|^\.[\p{L}\p{N}_]+|[\p{L}\p{N}_]+\.\n?$")
        .expect("the structural allele pattern")
});

/// The classes one sample's GT shows at a record, read as Ensembl VEP reads a genotype for
/// `--individual` (ensembl-io BaseVCF4.pm 941-951 and Parser/VCF.pm 627-651, release 116):
/// the bits split on `|` when the GT has one, else on `/`; `.` bits dropped before anything
/// else is read; an allele containing `*` outside the non-reference set. `missing` (no
/// called bit) stands alone; otherwise `partially_missing` when a bit was dropped, `haploid`
/// for a one-bit genotype, `star` for a `*` allele, and one of `multi_allelic` (two or more
/// distinct non-reference alleles), `hom_ref` (none), `hom_alt` (one, every called bit the
/// same) or `het`. A bit that is not an index of the record's alleles in ASCII digits is
/// refused, naming `record`.
pub fn genotype_labels<S: AsRef<str>>(
    gt: &str,
    alts: &[S],
    record: &str,
) -> Result<BTreeSet<&'static str>, Failure> {
    let bits: Vec<&str> = if gt.contains('|') {
        gt.split('|').collect()
    } else {
        gt.split('/').collect()
    };
    let called_bits: Vec<&str> = bits.iter().copied().filter(|bit| *bit != ".").collect();
    if called_bits.is_empty() {
        return Ok(BTreeSet::from(["missing"]));
    }
    let mut called: BTreeSet<usize> = BTreeSet::new();
    for bit in &called_bits {
        // A number past the machine's range is past the alleles too.
        let index = bit.parse::<usize>().ok().filter(|_| is_ascii_digits(bit));
        match index {
            Some(index) if index <= alts.len() => called.insert(index),
            _ => {
                return Err(failure(format!(
                    "record {record} has a GT {} outside its alleles",
                    py_repr_str(gt)
                )))
            }
        };
    }
    let mut labels = BTreeSet::new();
    if called_bits.len() < bits.len() {
        labels.insert("partially_missing");
    }
    if bits.len() == 1 {
        labels.insert("haploid");
    }
    let is_star = |index: &usize| *index > 0 && alts[*index - 1].as_ref().contains('*');
    if called.iter().any(is_star) {
        labels.insert("star");
    }
    let non_ref = called
        .iter()
        .filter(|index| **index != 0 && !is_star(index))
        .count();
    labels.insert(if non_ref > 1 {
        "multi_allelic"
    } else if non_ref == 0 {
        "hom_ref"
    } else if called.len() == 1 {
        "hom_alt"
    } else {
        "het"
    });
    Ok(labels)
}

/// Whether Ensembl VEP reads a record as a structural variant, which it annotates before
/// the per-individual step and so without any individual's genotype (Parser/VCF.pm
/// 236-244): its alleles are not all plain bases and it carries an INFO `SVTYPE` whose value
/// Perl reads as true, or a symbolic, breakend or single-breakend allele (`<*>` is none of
/// these). The value is read as ensembl-io's `get_info` stores it (BaseVCF4.pm 576-589): the
/// piece after the first `=`, so a bare `SVTYPE` flag and `SVTYPE=` hold none, `SVTYPE=0`
/// holds Perl's false `0`, and of a repeated key the last one counts. The allele test runs
/// on the comma-joined ALT string as Ensembl VEP's does, so a single-breakend allele counts
/// only at the string's start or end.
pub fn is_structural_record<S: AsRef<str>>(ref_allele: &str, alts: &[S], info: &str) -> bool {
    let alts: Vec<&str> = alts.iter().map(AsRef::as_ref).collect();
    let plain = |allele: &&str| allele.bytes().all(|byte| b"ACGT".contains(&byte));
    let any_base = !ref_allele.is_empty() || alts.iter().any(|alt| !alt.is_empty());
    if any_base && plain(&ref_allele) && alts.iter().all(plain) {
        return false;
    }
    let mut svtype = "";
    for part in info.split(';') {
        let mut pieces = part.split('=');
        if pieces.next() == Some("SVTYPE") {
            svtype = pieces.next().unwrap_or("");
        }
    }
    if !svtype.is_empty() && svtype != "0" {
        return true;
    }
    STRUCTURAL_ALLELE.is_match(&alts.join(","))
}

/// The genotype classes one record shows among its first `keep` samples, `fields` its
/// columns: every kept sample's `genotype_labels`, plus `all_missing` (no kept sample has a
/// called bit), `no_carrier` (a called bit somewhere but no non-reference allele in any kept
/// sample) and `one_hom_ref_rest_missing` (one kept sample `hom_ref`, every other one
/// `missing`). A structural record shows no class (`is_structural_record`). A record of
/// fewer than `keep` sample columns and one without `GT` in its FORMAT are refused.
pub fn record_genotype_classes<S: AsRef<str>>(
    fields: &[S],
    keep: usize,
) -> Result<BTreeSet<&'static str>, Failure> {
    let fields: Vec<&str> = fields.iter().map(AsRef::as_ref).collect();
    let record = fields[..fields.len().min(5)].join(" ");
    if fields.len() < 9 + keep {
        return Err(failure(format!(
            "record {record} has {} sample columns, {keep} are kept",
            fields.len().saturating_sub(9)
        )));
    }
    let alts: Vec<&str> = fields[4].split(',').collect();
    if is_structural_record(fields[3], &alts, fields[7]) {
        return Ok(BTreeSet::new());
    }
    let Some(gt_index) = fields[8].split(':').position(|key| key == "GT") else {
        return Err(failure(format!("record {record} has no GT in its FORMAT")));
    };
    let mut classes = BTreeSet::new();
    let (mut missing, mut hom_ref) = (0, 0);
    for sample in &fields[9..9 + keep] {
        let gt = sample.split(':').nth(gt_index).unwrap_or(".");
        let labels = genotype_labels(gt, &alts, &record)?;
        missing += usize::from(labels.contains("missing"));
        hom_ref += usize::from(labels.contains("hom_ref"));
        classes.extend(labels);
    }
    if missing == keep {
        classes.insert("all_missing");
    } else if missing + hom_ref == keep {
        classes.insert("no_carrier");
        if hom_ref == 1 {
            classes.insert("one_hom_ref_rest_missing");
        }
    }
    Ok(classes)
}

/// The `--require-genotype-class` names of a command line, each a genotype class, in order
/// and without repeats; another name is refused.
pub fn genotype_class_names(names: &[String]) -> Result<Vec<String>, Failure> {
    let names = without_repeats(names.iter().cloned());
    for name in &names {
        if !GENOTYPE_CLASSES.contains(&name.as_str()) {
            return Err(failure(format!(
                "--require-genotype-class wants one of {}, got {}",
                GENOTYPE_CLASSES.join(", "),
                py_repr_str(name)
            )));
        }
    }
    Ok(names)
}

/// The count of records a record class keeps as candidates for its exemplars.
pub const RECORD_CLASS_RESERVOIR: usize = 16;

/// The `--record-classes` key of a command line, one of `RECORD_CLASS_KEYS`; `None` when
/// the option is not given, and another key is refused.
pub fn record_class_key(name: Option<&str>) -> Result<Option<&str>, Failure> {
    match name {
        Some(name) if !RECORD_CLASS_KEYS.contains(&name) => Err(failure(format!(
            "--record-classes wants one of {}, got {}",
            RECORD_CLASS_KEYS.join(", "),
            py_repr_str(name)
        ))),
        _ => Ok(name),
    }
}

/// The shape of a row's allele, `snv`, `mnv` or `indel`, from the input record's REF and
/// ALT pair the allele came from: a pair of one length is a substitution, of two lengths an
/// indel. A two-base substitution and a two-base insertion print the same Location and
/// Allele (Ensembl VEP writes an insertion between its two flanking bases), so the row alone
/// cannot tell them apart. The pairs are read as Ensembl VEP's parser leaves them (the rule
/// of `record_names`): a bi-allelic indel loses its anchor base and the bases shared at
/// either end, a multi-allelic record with an indel loses the first base of every allele
/// when all share it and nothing more. A row whose allele matches no pair, or a record with
/// no pairs, falls back to `variant_class`, which reads that shape as an insertion.
pub fn allele_shape(row: &Row, alleles: &[(String, String)]) -> &'static str {
    let allele = row.allele.to_uppercase();
    let pairs: Vec<(String, String)> = alleles
        .iter()
        .filter(|(_, alt)| !is_symbolic_or_breakend(alt))
        .map(|(ref_allele, alt)| (ref_allele.to_uppercase(), alt.to_uppercase()))
        .collect();
    let multi_indel = pairs.len() > 1 && pairs.iter().any(|(r, a)| length(r) != length(a));
    let firsts: HashSet<&str> = pairs
        .iter()
        .flat_map(|(r, a)| [r.as_str(), a.as_str()])
        .filter(|allele| !allele.contains('*'))
        .map(first)
        .collect();
    let shared_first = firsts.len() == 1;
    for (ref_allele, alt) in &pairs {
        let mut a = alt.as_str();
        if multi_indel && shared_first {
            // A `*` allele keeps its one character, as `record_names` keeps it.
            if !a.contains('*') {
                a = rest(a);
            }
        } else if !multi_indel && length(ref_allele) != length(a) {
            let mut r = ref_allele.as_str();
            if first(r) == first(a) {
                (r, a) = (rest(r), rest(a));
            }
            while !r.is_empty() && !a.is_empty() && first(r) == first(a) {
                (r, a) = (rest(r), rest(a));
            }
            while !r.is_empty() && !a.is_empty() && last(r) == last(a) {
                (r, a) = (without_last(r), without_last(a));
            }
        }
        if or_dash(a) != allele {
            continue;
        }
        return if length(ref_allele) != length(alt) {
            "indel"
        } else if length(alt) == 1 {
            "snv"
        } else {
            "mnv"
        };
    }
    match variant_class(&row.location, &row.allele, false).as_str() {
        "snv" => "snv",
        "mnv_or_complex" => "mnv",
        _ => "indel",
    }
}

/// The `regulatory_class` classes a record's reference rows show, eight in all: `rf_only`,
/// `mf_only` or `rf_and_mf` for the regulatory feature types among the rows' `Feature_type`;
/// `regulatory_only` when the record has one of them and no Transcript row (its remaining
/// row is the intergenic one Ensembl VEP writes for a variant no transcript overlaps); per
/// MotifFeature row `motif_snv`, `motif_mnv` or `motif_indel` by the shape of its allele
/// (`allele_shape`, over the record's REF and ALT pairs), and `motif_reverse` for a motif on
/// the reverse strand. A record with no regulatory row shows no class.
pub fn regulatory_classes(rows: &[Row], alleles: &[(String, String)]) -> BTreeSet<&'static str> {
    let has = |feature_type: &str| rows.iter().any(|row| row.feature_type == feature_type);
    let (has_rf, has_mf) = (has("RegulatoryFeature"), has("MotifFeature"));
    let mut classes = BTreeSet::new();
    if !(has_rf || has_mf) {
        return classes;
    }
    classes.insert(match (has_rf, has_mf) {
        (true, true) => "rf_and_mf",
        (true, false) => "rf_only",
        _ => "mf_only",
    });
    if !has("Transcript") {
        classes.insert("regulatory_only");
    }
    for row in rows.iter().filter(|row| row.feature_type == "MotifFeature") {
        if SEQUENCE_CLASSES.contains(&row.klass().as_str()) {
            classes.insert(match allele_shape(row, alleles) {
                "snv" => "motif_snv",
                "mnv" => "motif_mnv",
                _ => "motif_indel",
            });
        }
        if row.strand == "-1" {
            classes.insert("motif_reverse");
        }
    }
    classes
}

/// What classes the rows of one input record show under a record-class key, given the REF
/// and ALT pairs of that record.
pub type RecordClassFunction = fn(&[Row], &[(String, String)]) -> BTreeSet<&'static str>;

/// The function of each record-class key, in the order of `RECORD_CLASS_KEYS`.
pub const RECORD_CLASS_FUNCTIONS: [(&str, RecordClassFunction); 1] =
    [(REGULATORY_CLASS, regulatory_classes)];

/// The function of a record-class key, `None` for a text that is no key.
pub fn record_class_function(key: &str) -> Option<RecordClassFunction> {
    RECORD_CLASS_FUNCTIONS
        .iter()
        .find(|(name, _)| *name == key)
        .map(|(_, function)| *function)
}

#[cfg(test)]
#[path = "classes_tests.rs"]
mod tests;
