// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A record's alleles, sibling rows, chromosome, allele lengths and alternates.

use crate::consts::{HGVS_OFFSET_FIELDS, NON_ZERO_INTEGER_RE_PATTERN};
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::{entry_kind, record_label};
use crate::pyfmt;
use crate::rules::downstream::field_of;
use crate::rules::{nearest_i64, whole_number, whole_text};
use crate::shared::{allele_is_indel, dash, indel_row, uploaded_alleles, ABSENT};
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// `NON_ZERO_INTEGER_RE` held to the whole text, as `fullmatch` holds it.
static NON_ZERO_INTEGER_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));

/// The first record entry among `record_rows` under `e`'s own label.
fn record_entry<'a>(e: &Entry, record_rows: &[&'a Entry]) -> Option<&'a Entry> {
    let label = record_label(e);
    record_rows
        .iter()
        .copied()
        .find(|r| entry_kind(r) == "record" && record_label(r) == label)
}

/// The identifier-less `Uploaded_variation` of `e` as its reference and alternate alleles.
fn uploaded(e: &Entry) -> Option<(String, Vec<String>)> {
    uploaded_alleles(
        e.fields
            .get("Uploaded_variation")
            .map_or("", String::as_str),
    )
}

/// The parts of a text between its `separator`s, each owned.
fn parts(text: &str, separator: char) -> Vec<String> {
    text.split(separator).map(str::to_string).collect()
}

/// Whether `allele` is one or more of the bases `ACGTN` and nothing else.
fn sequence(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| b"ACGTN".contains(&b))
}

/// The whole number `offset` writes where it is a non-zero integer (`NON_ZERO_INTEGER_RE`);
/// `None` for any other text.
fn non_zero_integer(offset: &str) -> Option<i64> {
    NON_ZERO_INTEGER_WHOLE_RE
        .is_match(offset)
        .then(|| nearest_i64(whole_number(offset)))
}

/// The alleles the rows of `e`'s record carry, its own among them; an entry of
/// `record_rows` under another record's label (`record_label`) is not read.
pub fn record_alleles(e: &Entry, record_rows: &[&Entry]) -> BTreeSet<String> {
    let label = record_label(e);
    let mut alleles = BTreeSet::from([e.allele.clone()]);
    alleles.extend(
        record_rows
            .iter()
            .filter(|r| entry_kind(r) != "record" && record_label(r) == label)
            .map(|r| r.allele.clone()),
    );
    alleles
}

/// The rows of `e`'s record on its feature under another allele, borrowed from
/// `record_rows` in its order; an entry under another record's label is none of them.
pub fn sibling_rows_of<'a>(e: &Entry, record_rows: &[&'a Entry]) -> Vec<&'a Entry> {
    let (kind, label) = (entry_kind(e), record_label(e));
    record_rows
        .iter()
        .copied()
        .filter(|r| {
            r.allele != e.allele
                && r.feature == e.feature
                && entry_kind(r) == kind
                && record_label(r) == label
        })
        .collect()
}

/// The rows of `e`'s record on its feature under another allele; an entry of
/// `record_rows` under another record's label (`record_label`) is none of them.
pub fn sibling_rows(e: &Entry, record_rows: &[&Entry]) -> Vec<Entry> {
    sibling_rows_of(e, record_rows)
        .into_iter()
        .cloned()
        .collect()
}

/// The chromosome of `e`'s record as its output spells it: the Location's in the default
/// and tab formats, else the `CHROM` of the VCF format's record entry or the
/// `seq_region_name` of the JSON format's among `record_rows`, under `e`'s own label
/// (`record_label`); `None` where none is carried.
pub fn record_chromosome(e: &Entry, record_rows: Option<&[&Entry]>) -> Option<String> {
    if let Some(location) = e.location.as_deref() {
        let chromosome = pyfmt::rpartition(location, ":").0;
        return (!chromosome.is_empty()).then(|| chromosome.to_string());
    }
    let record = record_entry(e, record_rows.unwrap_or(&[]))?;
    ["CHROM", "seq_region_name"]
        .iter()
        .find_map(|field| record.fields.get(*field))
        .cloned()
}

/// `indel_row` read from the record's own alleles where the parser did not read them and
/// `record_rows` carries them: the `REF` and `ALT` of the record entry built from the
/// input line of a default or tab row under an identifier. The Location alone reads an
/// insertion written at a one-base Location, the `CA` of `C/-/CA`, as a substitution.
pub fn record_indel_row(e: &Entry, record_rows: &[&Entry]) -> Result<bool> {
    if e.allele != "-" && e.indel.is_none() {
        let label = record_label(e);
        for r in record_rows {
            if entry_kind(r) != "record" || record_label(r) != label {
                continue;
            }
            let (Some(reference), Some(alternates)) = (r.fields.get("REF"), r.fields.get("ALT"))
            else {
                continue;
            };
            let alternates: Vec<&str> = alternates.split(',').collect();
            if let Some(shape) = allele_is_indel(reference, &alternates, &e.allele, false) {
                return Ok(shape);
            }
        }
    }
    indel_row(e)
}

/// The lengths of the reference allele and of `e`'s allele as Ensembl VEP hands them to
/// a plugin, an empty allele 0: read from the record's own alleles (an identifier-less
/// `Uploaded_variation`, else the record entry among `record_rows`: a VCF line's `REF`
/// and `ALT`, the JSON `allele_string`), else from the row's text where it states them.
/// A deletion's Location spans the deleted bases and a Location of one base or of three or
/// more is the reference allele's; a two-base Location is an insertion's two flanking bases
/// or a two-base reference allele, so `None` there, as for an allele that is no sequence
/// and one that no alternate, or more than one, spells.
///
/// A length past `i64` is the nearest `i64`.
pub fn record_allele_lengths(e: &Entry, record_rows: Option<&[&Entry]>) -> Option<(i64, i64)> {
    if !(sequence(&e.allele) || e.allele == "-") {
        return None;
    }
    let alleles = match uploaded(e) {
        Some((reference, alternates)) => Some(vcf_alleles_as_spelt(&reference, &alternates)),
        None => record_entry(e, record_rows.unwrap_or(&[])).and_then(|record| {
            let fields = &record.fields;
            if let (Some(reference), Some(alternates)) = (fields.get("REF"), fields.get("ALT")) {
                Some(vcf_alleles_as_spelt(reference, &parts(alternates, ',')))
            } else {
                fields.get("allele_string").map(|spelt| {
                    let mut alleles = parts(spelt, '/');
                    let reference = alleles.remove(0);
                    (reference, alleles)
                })
            }
        }),
    };
    let bases = |allele: &str| pyfmt::len(pyfmt::strip_chars(allele, "-")) as i64;
    if let Some((reference, alternates)) = alleles {
        let spelt = alternates.iter().filter(|a| **a == e.allele).count();
        return (spelt == 1).then(|| (bases(&reference), bases(&e.allele)));
    }
    let location = e.location.as_deref()?;
    let coordinates = location.split_once(':').map_or(location, |(_, rest)| rest);
    let (start, _, end) = pyfmt::partition(coordinates, "-");
    if !(pyfmt::is_digit(start) && (pyfmt::is_digit(end) || end.is_empty())) {
        return None;
    }
    let last = whole_number(if end.is_empty() { start } else { end });
    let span = nearest_i64(last.saturating_sub(whole_number(start)).saturating_add(1));
    if e.allele == "-" {
        return Some((span, 0));
    }
    (span != 2).then(|| (span, pyfmt::len(&e.allele) as i64))
}

/// A VCF line's REF and ALTs as its entries spell them: Parser/VCF.pm:300-323 drop the
/// first base of a record that holds an ALT of another length than its REF, a `*` among
/// them, and whose alleles not spelt with `*` all share it, an empty remainder written `-`.
pub fn vcf_alleles_as_spelt(reference: &str, alternates: &[String]) -> (String, Vec<String>) {
    let holds_indel = alternates
        .iter()
        .any(|a| pyfmt::len(a) != pyfmt::len(reference));
    let shared = alternates
        .iter()
        .filter(|a| !a.contains('*'))
        .all(|a| pyfmt::prefix(a, 1) == pyfmt::prefix(reference, 1));
    if !(holds_indel && !reference.is_empty() && shared) {
        return (reference.to_string(), alternates.to_vec());
    }
    let spelt = alternates
        .iter()
        .map(|a| {
            if a.contains('*') {
                a.clone()
            } else {
                dash(pyfmt::slice(a, Some(1), None))
            }
        })
        .collect();
    (dash(pyfmt::slice(reference, Some(1), None)), spelt)
}

/// The first ALT of a VCF line as its entries spell the allele (`vcf_alleles_as_spelt`).
pub fn vcf_first_alternate(reference: &str, alternates: &[String]) -> String {
    let (_, spelt) = vcf_alleles_as_spelt(reference, alternates);
    match spelt.into_iter().next() {
        Some(first) => first,
        None => crate::uncaught!("list index out of range"),
    }
}

/// The first alternate allele of `e`'s record where an entry states the record's alleles
/// in order, spelt as its rows spell the allele: a text row's identifier-less
/// `Uploaded_variation`, which names a VCF record by the line's own REF and ALT, else the
/// `ALT` of the VCF format's record entry or the `allele_string` of the JSON format's
/// among `record_rows`; `None` where none does.
pub fn first_alternate(e: &Entry, record_rows: &[&Entry]) -> Option<String> {
    if let Some((reference, alternates)) = uploaded(e) {
        return Some(vcf_first_alternate(&reference, &alternates));
    }
    let fields = &record_entry(e, record_rows)?.fields;
    if let Some(alternates) = fields.get("ALT") {
        let reference = fields.get("REF").map_or("", String::as_str);
        return Some(vcf_first_alternate(reference, &parts(alternates, ',')));
    }
    let spelt = fields.get("allele_string").map_or("", String::as_str);
    spelt.split('/').nth(1).map(str::to_string)
}

/// The alternate alleles of `e`'s record in the record's own order, spelt as its rows
/// spell the allele, where an entry states them (`first_alternate` names the entries);
/// `None` where none does.
pub fn alternates_in_order(e: &Entry, record_rows: Option<&[&Entry]>) -> Option<Vec<String>> {
    if let Some((reference, alternates)) = uploaded(e) {
        return Some(vcf_alleles_as_spelt(&reference, &alternates).1);
    }
    let fields = &record_entry(e, record_rows.unwrap_or(&[]))?.fields;
    if let Some(alternates) = fields.get("ALT") {
        let reference = fields.get("REF").map_or("", String::as_str);
        return Some(vcf_alleles_as_spelt(reference, &parts(alternates, ',')).1);
    }
    let spelt = fields.get("allele_string").map_or("", String::as_str);
    let alternates: Vec<String> = spelt.split('/').skip(1).map(str::to_string).collect();
    (!alternates.is_empty()).then_some(alternates)
}

/// Whether `allele` is a second or later one of `alternates` that one of them alone spells.
fn later_among(alternates: &[String], allele: &str) -> bool {
    alternates.iter().filter(|a| *a == allele).count() == 1 && alternates[0] != allele
}

/// Whether `e` is the row of a second or later alternate allele of its record: an entry
/// states the record's alternates in order (`alternates_in_order`), one of them alone
/// spells the row's allele, and it is not the first.
pub fn later_alternate(e: &Entry, record_rows: Option<&[&Entry]>) -> bool {
    alternates_in_order(e, record_rows)
        .is_some_and(|alternates| later_among(&alternates, &e.allele))
}

/// The HGVS_OFFSET of the reference's row of the first alternate allele of `e`'s record
/// on `e`'s feature, for the row `e` of a later alternate (`later_alternate`); `None`
/// for any other row, without that one row among `record_rows` and where it carries no
/// non-zero offset.
///
/// An offset past `i64` is the nearest `i64`.
pub fn first_allele_shift(e: &Entry, record_rows: Option<&[&Entry]>) -> Option<i64> {
    let alternates = alternates_in_order(e, record_rows)?;
    if !later_among(&alternates, &e.allele) {
        return None;
    }
    let rows: Vec<&Entry> = sibling_rows_of(e, record_rows.unwrap_or(&[]))
        .into_iter()
        .filter(|r| r.allele == alternates[0])
        .collect();
    let offset = match rows[..] {
        [row] => field_of(Some(row), HGVS_OFFSET_FIELDS),
        _ => ABSENT.to_string(),
    };
    non_zero_integer(&offset)
}

/// The non-zero HGVS_OFFSET of the reference's row, on `e`'s feature, of the last
/// alternate allele before `e`'s own in its record's order whose row carries one, for the
/// row `e` of a later alternate (`later_alternate`): the shift whose CDS pair that
/// allele's protein notation left on the TranscriptVariation the record's alleles share,
/// where each allele's own notation moves it afresh (TranscriptVariationAllele.pm:844-858).
/// `None` for any other row, where no earlier allele's row carries an offset, and where an
/// earlier allele has two rows on the feature among `record_rows`.
///
/// An offset past `i64` is the nearest `i64`.
pub fn earlier_allele_shift(e: &Entry, record_rows: Option<&[&Entry]>) -> Option<i64> {
    let alternates = alternates_in_order(e, record_rows)?;
    if !later_among(&alternates, &e.allele) {
        return None;
    }
    let siblings = sibling_rows_of(e, record_rows.unwrap_or(&[]));
    let own = alternates.iter().position(|a| *a == e.allele)?;
    for allele in alternates[..own].iter().rev() {
        let rows: Vec<&Entry> = siblings
            .iter()
            .copied()
            .filter(|r| r.allele == *allele)
            .collect();
        let offset = match rows[..] {
            [] => ABSENT.to_string(),
            [row] => field_of(Some(row), HGVS_OFFSET_FIELDS),
            _ => return None,
        };
        if let Some(shift) = non_zero_integer(&offset) {
            return Some(shift);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields(pairs: &[(&str, &str)]) -> crate::entry::Fields {
        pairs
            .iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect()
    }

    /// A row of the VCF or JSON format under the record ordinal `record`.
    fn row(record: i64, allele: &str, feature: &str, pairs: &[(&str, &str)]) -> Entry {
        Entry {
            key: format!("{record}|{allele}|{feature}|Transcript"),
            record: Some(record),
            allele: allele.to_string(),
            feature: feature.to_string(),
            fields: fields(pairs),
            ..Entry::default()
        }
    }

    /// The record entry of the VCF or JSON format under the record ordinal `record`.
    fn record(record: i64, pairs: &[(&str, &str)]) -> Entry {
        Entry {
            key: format!("{record}|record"),
            record: Some(record),
            fields: fields(pairs),
            ..Entry::default()
        }
    }

    /// A row of the default or tab format.
    fn text_row(location: &str, allele: &str, uploaded: &str) -> Entry {
        Entry {
            key: format!("{location}|{allele}|ENST1|Transcript"),
            location: Some(location.to_string()),
            allele: allele.to_string(),
            feature: "ENST1".to_string(),
            fields: fields(&[("Uploaded_variation", uploaded)]),
            ..Entry::default()
        }
    }

    fn owned(alleles: &[&str]) -> Vec<String> {
        alleles.iter().map(|a| a.to_string()).collect()
    }

    /// Three alleles of one record on one feature, a row on another feature, a row of
    /// another record and the two record entries.
    fn vcf_rows() -> Vec<Entry> {
        vec![
            row(7, "-", "ENST1", &[("HGVS_OFFSET", "3")]),
            row(7, "TT", "ENST1", &[("HGVS_OFFSET", "-2")]),
            row(7, "G", "ENST1", &[]),
            row(7, "-", "ENST2", &[]),
            row(8, "T", "ENST1", &[]),
            record(
                7,
                &[
                    ("CHROM", "chrM"),
                    ("POS", "5"),
                    ("REF", "AT"),
                    ("ALT", "A,ATT,AG"),
                ],
            ),
            record(
                8,
                &[("CHROM", "21"), ("POS", "9"), ("REF", "C"), ("ALT", "T")],
            ),
        ]
    }

    #[test]
    fn the_alleles_of_a_line_lose_the_first_base_they_all_share() {
        let spelt = |reference: &str, alternates: &[&str]| {
            let (reference, alternates) = vcf_alleles_as_spelt(reference, &owned(alternates));
            (reference, alternates.join(" "))
        };
        assert_eq!(spelt("A", &["T"]), ("A".to_string(), "T".to_string()));
        assert_eq!(
            spelt("AT", &["A", "ATT"]),
            ("T".to_string(), "- TT".to_string())
        );
        // An allele spelt with a star keeps its text and is not asked for the shared base.
        assert_eq!(
            spelt("AT", &["A", "*"]),
            ("T".to_string(), "- *".to_string())
        );
        assert_eq!(
            spelt("A", &["<DEL:*>"]),
            ("-".to_string(), "<DEL:*>".to_string())
        );
        assert_eq!(
            spelt("AT", &["A", "G"]),
            ("AT".to_string(), "A G".to_string())
        );
        assert_eq!(spelt("", &["A"]), (String::new(), "A".to_string()));
        assert_eq!(spelt("A", &[]), ("A".to_string(), String::new()));
        assert_eq!(vcf_first_alternate("AT", &owned(&["A", "ATT"])), "-");
    }

    #[test]
    fn a_record_gives_its_rows_their_alleles_siblings_and_chromosome() {
        let rows = vcf_rows();
        let rows: Vec<&Entry> = rows.iter().collect();
        let e = rows[2];
        let alleles: Vec<String> = record_alleles(e, &rows).into_iter().collect();
        assert_eq!(alleles, ["-", "G", "TT"]);
        let siblings: Vec<(&str, &str)> = sibling_rows_of(e, &rows)
            .iter()
            .map(|r| (r.allele.as_str(), r.feature.as_str()))
            .collect();
        assert_eq!(siblings, [("-", "ENST1"), ("TT", "ENST1")]);
        assert_eq!(sibling_rows(e, &rows).len(), 2);
        assert_eq!(record_chromosome(e, Some(&rows)).as_deref(), Some("chrM"));
        assert_eq!(record_chromosome(e, None), None);
        assert_eq!(first_alternate(e, &rows).as_deref(), Some("-"));
        assert_eq!(
            alternates_in_order(e, Some(&rows)),
            Some(owned(&["-", "TT", "G"]))
        );
        let later: Vec<bool> = rows[..5]
            .iter()
            .map(|r| later_alternate(r, Some(&rows)))
            .collect();
        assert_eq!(later, [false, true, true, false, false]);
        let lengths: Vec<_> = rows[..5]
            .iter()
            .map(|r| record_allele_lengths(r, Some(&rows)))
            .collect();
        assert_eq!(
            lengths,
            [
                Some((1, 0)),
                Some((1, 2)),
                Some((1, 1)),
                Some((1, 0)),
                Some((1, 1))
            ]
        );
        let indel: Vec<bool> = rows[..5]
            .iter()
            .map(|r| record_indel_row(r, &rows).unwrap())
            .collect();
        assert_eq!(indel, [true, true, false, true, false]);
    }

    #[test]
    fn a_later_allele_reads_the_shift_of_the_alleles_before_it() {
        let rows = vcf_rows();
        let rows: Vec<&Entry> = rows.iter().collect();
        let first: Vec<_> = rows[..3]
            .iter()
            .map(|r| first_allele_shift(r, Some(&rows)))
            .collect();
        assert_eq!(first, [None, Some(3), Some(3)]);
        let earlier: Vec<_> = rows[..3]
            .iter()
            .map(|r| earlier_allele_shift(r, Some(&rows)))
            .collect();
        assert_eq!(earlier, [None, Some(3), Some(-2)]);
    }

    #[test]
    fn a_json_record_states_its_alleles_in_its_allele_string() {
        let rows = [
            row(3, "TT", "ENST1", &[]),
            row(3, "-", "ENST1", &[]),
            record(3, &[("seq_region_name", "MT"), ("allele_string", "T/TT/-")]),
        ];
        let rows: Vec<&Entry> = rows.iter().collect();
        assert_eq!(
            record_chromosome(rows[0], Some(&rows)).as_deref(),
            Some("MT")
        );
        assert_eq!(first_alternate(rows[1], &rows).as_deref(), Some("TT"));
        assert_eq!(
            alternates_in_order(rows[0], Some(&rows)),
            Some(owned(&["TT", "-"]))
        );
        assert_eq!(record_allele_lengths(rows[0], Some(&rows)), Some((1, 2)));
        assert_eq!(record_allele_lengths(rows[1], Some(&rows)), Some((1, 0)));
        // An allele string of the reference alone states no alternate.
        let bare = record(4, &[("allele_string", "T")]);
        let e = row(4, "T", "ENST1", &[]);
        assert_eq!(first_alternate(&e, &[&bare]), None);
        assert_eq!(alternates_in_order(&e, Some(&[&bare])), None);
    }

    #[test]
    fn a_text_row_states_its_alleles_in_a_generated_name_or_its_location() {
        let inserted = text_row("21:100-101", "CA", "21_100_C/-/CA");
        let deleted = text_row("21:100-101", "-", "21_100_C/-/CA");
        assert_eq!(first_alternate(&inserted, &[]).as_deref(), Some("-"));
        assert_eq!(
            alternates_in_order(&inserted, None),
            Some(owned(&["-", "CA"]))
        );
        assert_eq!(record_allele_lengths(&inserted, None), Some((1, 2)));
        assert_eq!(record_allele_lengths(&deleted, None), Some((1, 0)));
        assert_eq!(record_chromosome(&inserted, None).as_deref(), Some("21"));
        assert_eq!(record_chromosome(&text_row("100", "A", "rs1"), None), None);
        let lengths = |location: &str, allele: &str| {
            record_allele_lengths(&text_row(location, allele, "rs1"), None)
        };
        assert_eq!(lengths("21:100", "A"), Some((1, 1)));
        // Two bases are an insertion's flanks or a two-base reference allele.
        assert_eq!(lengths("21:100-101", "A"), None);
        assert_eq!(lengths("21:100-102", "A"), Some((3, 1)));
        assert_eq!(lengths("21:100-102", "-"), Some((3, 0)));
        assert_eq!(lengths("21:x", "A"), None);
        assert_eq!(lengths("21:100-", "ACG"), Some((1, 3)));
        assert_eq!(lengths("21:5-3", "A"), Some((-1, 1)));
        assert_eq!(lengths("21:100", "<DEL>"), None);
    }

    #[test]
    fn a_row_under_an_identifier_reads_its_shape_from_the_record_entry() {
        let named = text_row("21:100", "CA", "rs9");
        let record = Entry {
            key: "21:100|record".to_string(),
            location: Some("21:100".to_string()),
            fields: fields(&[("Uploaded_variation", "rs9"), ("REF", "C"), ("ALT", "-,CA")]),
            ..Entry::default()
        };
        assert!(!record_indel_row(&named, &[]).unwrap());
        assert!(record_indel_row(&named, &[&record]).unwrap());
        assert_eq!(
            record_allele_lengths(&named, Some(&[&record])),
            Some((1, 2))
        );
    }
}
