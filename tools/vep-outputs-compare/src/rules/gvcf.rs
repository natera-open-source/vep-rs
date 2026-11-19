// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! gVCF variant sites: a `<NON_REF>` block beside a real allele, which Ensembl VEP routes to
//! the structural path and drops, recognised from the input VCF's shape or from the generated
//! name of an ID-less record in vep-rs's output.

use super::breakend::vep_rs_pos;
use crate::input::{contig_of, for_each_record, open_vcf, push_canonical_location};
use crate::pyfmt::py_strip;
use regex::Regex;
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::LazyLock;

/// `CHROM_POS_REF/ALT...`, the name VEP gives an ID-less record; the chromosome may carry
/// underscores, the position is the first `_digits_` group after it.
static GENERATED_NAME: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^.+?_[0-9]+_(.+)$").expect("regex"));

/// Parser/VCF.pm:241, the ALT forms that route a line to the structural path.
static STRUCTURAL_ALT_LIST: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$").expect("regex"));

const NON_REF_BLOCK: &str = "<NON_REF>";

fn is_non_ref_block(alt: &str) -> bool {
    alt.to_uppercase() == NON_REF_BLOCK
}

/// Whether an ALT names a variant: not the `<NON_REF>` block, `*` or `.`.
fn names_an_allele(alt: &str) -> bool {
    !is_non_ref_block(alt) && alt != "*" && alt != "."
}

/// Whether Perl reads an INFO value as true: defined, not empty and not `0`.
fn perl_true(value: Option<&str>) -> bool {
    !matches!(value, None | Some("") | Some("0"))
}

/// Drops the trailing empty fields, as Perl's list-context `split` drops them.
fn drop_trailing_empty(fields: &mut Vec<&str>) {
    while fields.last() == Some(&"") {
        fields.pop();
    }
}

/// Perl's list-context `split`: the trailing empty fields dropped, the middle ones kept, so
/// `G,,T,,` is `G`, the empty field, `T`, and the empty string has no field.
pub fn perl_split(text: &str, sep: char) -> Vec<&str> {
    let mut fields: Vec<&str> = text.split(sep).collect();
    drop_trailing_empty(&mut fields);
    fields
}

/// The value of `key=value` in an INFO column, the last occurrence winning as in the hash
/// Ensembl VEP's parser builds (BaseVCF4.pm:581-583), the text up to a second `=` kept; None
/// for an absent key or a bare flag.
fn info_written<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    let mut value = None;
    for item in info.split(';') {
        match item.split_once('=') {
            Some((name, text)) if name == key => {
                value = Some(text.split_once('=').map_or(text, |(v, _)| v));
            }
            None if item == key => value = None,
            _ => {}
        }
    }
    value
}

/// Whether a record is a gVCF variant site vep-rs reads on its sequence path: `<NON_REF>`
/// beside alleles that name a variant, no true INFO/SVTYPE, and no remaining allele written
/// in a structural form (Parser/VCF.pm:237-245).
pub fn is_gvcf_variant_site(alts: &[&str], svtype: Option<&str>) -> bool {
    let mut alts = alts.to_vec();
    drop_trailing_empty(&mut alts);
    let real: Vec<&str> = alts
        .iter()
        .copied()
        .filter(|a| names_an_allele(a))
        .collect();
    if real.is_empty() || !alts.iter().any(|a| is_non_ref_block(a)) || perl_true(svtype) {
        return false;
    }
    !STRUCTURAL_ALT_LIST.is_match(&real.join(","))
}

fn same_letter(a: char, b: char) -> bool {
    a.to_uppercase().eq(b.to_uppercase())
}

/// A bi-allelic indel minimised as vep-rs minimises it: the shared prefix trimmed while both
/// alleles have bases, then the shared suffix, an emptied allele reading `-`; returns the two
/// alleles and the bases trimmed in front.
fn minimised(reference: &[char], alt: &[char]) -> (Vec<char>, Vec<char>, usize) {
    let mut shared = 0;
    while shared < reference.len().min(alt.len()) && same_letter(reference[shared], alt[shared]) {
        shared += 1;
    }
    let (mut reference, mut alt) = (&reference[shared..], &alt[shared..]);
    while let (Some(r), Some(a)) = (reference.last(), alt.last()) {
        if !same_letter(*r, *a) {
            break;
        }
        reference = &reference[..reference.len() - 1];
        alt = &alt[..alt.len() - 1];
    }
    let or_dash = |s: &[char]| if s.is_empty() { vec!['-'] } else { s.to_vec() };
    (or_dash(reference), or_dash(alt), shared)
}

/// The (allele, start, end) vep-rs writes for each real allele of a gVCF variant site, the
/// sequence path of its parser (Parser/VCF.pm:296-328; the minimisation of Parser.pm:978-984).
pub fn gvcf_site_rows(pos: i128, reference: &str, alts: &[&str]) -> Vec<(String, i128, i128)> {
    let mut alts = alts.to_vec();
    drop_trailing_empty(&mut alts);
    let named: Vec<&str> = alts
        .iter()
        .copied()
        .filter(|a| !is_non_ref_block(a))
        .collect();
    let real: Vec<&str> = named
        .iter()
        .copied()
        .filter(|a| names_an_allele(a))
        .collect();
    if real.is_empty() {
        return Vec::new();
    }
    let reference: Vec<char> = reference.chars().collect();
    let is_indel = named.iter().any(|a| a.chars().count() != reference.len());
    let compared: Vec<&str> = if named.len() == 1 {
        named.clone()
    } else {
        named.iter().copied().filter(|a| !a.contains('*')).collect()
    };
    let anchored = !reference.is_empty()
        && compared.iter().all(|a| {
            a.chars()
                .next()
                .is_some_and(|c| same_letter(c, reference[0]))
        });
    let prefix = usize::from(is_indel && anchored);
    let trimmed_ref: Vec<char> = if reference.len() > prefix {
        reference[prefix..].to_vec()
    } else {
        vec!['-']
    };
    let mut rows = Vec::new();
    for alt in real {
        let alt_chars: Vec<char> = alt.chars().collect();
        let mut allele: Vec<char> = if alt.contains('*') {
            alt_chars
        } else if alt_chars.len() > prefix {
            alt_chars[prefix..].to_vec()
        } else {
            vec!['-']
        };
        let mut start = pos + prefix as i128;
        let mut ref_part = trimmed_ref.clone();
        if named.len() == 1
            && ref_part != ['-']
            && allele != ['-']
            && ref_part.len() != allele.len()
        {
            let (r, a, shared) = minimised(&ref_part, &allele);
            ref_part = r;
            allele = a;
            start += shared as i128;
        }
        let end = if ref_part == ['-'] {
            start - 1
        } else {
            start + ref_part.len() as i128 - 1
        };
        let allele: String = allele.into_iter().collect();
        rows.push((allele.to_uppercase(), start, end));
    }
    rows
}

/// The Location both engines print for a span: `chrom:start`, or `chrom:lo-hi` with an
/// insertion's `end = start - 1` written first.
pub fn sequence_location(chrom: &str, start: i128, end: i128) -> String {
    let (lo, hi) = (start.min(end), start.max(end));
    if lo == hi {
        format!("{chrom}:{lo}")
    } else {
        format!("{chrom}:{lo}-{hi}")
    }
}

/// The Location of one row of vep-rs's output when its record is a gVCF variant site: an
/// ID-less record is named `CHROM_POS_REF/ALT...` from its own line, so a name whose ALT list
/// carries `<NON_REF>` beside another allele is the site. None for any other row and for a
/// Location outside `canonical_contigs`; `location` is scratch space.
pub fn gvcf_site_location(
    record: &[u8],
    canonical_contigs: Option<&HashSet<String>>,
    location: &mut Vec<u8>,
) -> Option<String> {
    // A row without `<` names no `<NON_REF>` block.
    if record.first() == Some(&b'#') || memchr::memchr(b'<', record).is_none() {
        return None;
    }
    let text = String::from_utf8_lossy(record);
    let fields: Vec<&str> = text.split('\t').collect();
    if fields.len() < 7 {
        return None;
    }
    let name = GENERATED_NAME.captures(fields[0])?;
    let parts = perl_split(name.get(1).expect("alleles group").as_str(), '/');
    if !is_gvcf_variant_site(parts.get(1..).unwrap_or(&[]), None) {
        return None;
    }
    location.clear();
    push_canonical_location(py_strip(fields[1]), location);
    if let Some(contigs) = canonical_contigs {
        let contig = std::str::from_utf8(contig_of(location)).expect("utf-8");
        if !contigs.contains(contig) {
            return None;
        }
    }
    Some(String::from_utf8(location.clone()).expect("utf-8"))
}

/// The Locations of the rows in vep-rs's output whose record is a gVCF variant site. A
/// missing file yields the empty set.
pub fn gvcf_site_locations(
    path: &Path,
    canonical_contigs: Option<&HashSet<String>>,
) -> io::Result<HashSet<String>> {
    let mut sites = HashSet::new();
    let Ok(file) = File::open(path) else {
        return Ok(sites);
    };
    let mut location = Vec::new();
    for_each_record(file, |record| {
        if let Some(site) = gvcf_site_location(record, canonical_contigs, &mut location) {
            sites.insert(site);
        }
        Ok(())
    })?;
    Ok(sites)
}

/// The Locations vep-rs writes for every gVCF variant site of an input VCF (plain or gzipped),
/// recognised by the record's shape whatever its ID, at the POS as vep-rs reads it; a record
/// whose POS vep-rs refuses has no rows.
pub fn gvcf_site_locations_from_vcf(
    path: &Path,
    canonical_contigs: Option<&HashSet<String>>,
) -> io::Result<HashSet<String>> {
    let mut sites = HashSet::new();
    let mut location = Vec::new();
    for_each_record(open_vcf(path)?, |record| {
        // A record without `<` lists no `<NON_REF>` block.
        if record.first() == Some(&b'#') || memchr::memchr(b'<', record).is_none() {
            return Ok(());
        }
        let text = String::from_utf8_lossy(record);
        let fields: Vec<&str> = text.split('\t').collect();
        if fields.len() < 5 {
            return Ok(());
        }
        let (chrom, reference, alt) = (fields[0], fields[3], fields[4]);
        let alts = perl_split(alt, ',');
        let svtype = fields.get(7).and_then(|info| info_written(info, "SVTYPE"));
        let Some(pos) = vep_rs_pos(fields[1]) else {
            return Ok(());
        };
        if !is_gvcf_variant_site(&alts, svtype) {
            return Ok(());
        }
        if canonical_contigs.is_some_and(|contigs| !contigs.contains(chrom)) {
            return Ok(());
        }
        for (_allele, start, end) in gvcf_site_rows(pos, reference, &alts) {
            location.clear();
            push_canonical_location(&sequence_location(chrom, start, end), &mut location);
            sites.insert(String::from_utf8(location.clone()).expect("utf-8"));
        }
        Ok(())
    })?;
    Ok(sites)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_generated_name_holds_its_position_in_the_digits_0_to_9() {
        assert!(GENERATED_NAME.is_match("21_100_A/G,<NON_REF>"));
        assert!(!GENERATED_NAME.is_match("21_١٠٠_A/G,<NON_REF>"));
    }
}
