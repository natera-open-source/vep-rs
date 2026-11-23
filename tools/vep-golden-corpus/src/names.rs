// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The names a record can carry in the first column of an Ensembl VEP output: the first
//! element of its ID, the names Ensembl VEP writes for a record without one, and the index
//! from a name to the records of a run it names.

use crate::error::Failure;
use crate::manifest::without_repeats;
use crate::pytext::{py_int, py_repr_str};
use std::collections::{HashMap, HashSet};

/// The strata of a multi-allelic record, in the order a manifest lists them.
pub const MULTI_ALLELIC_STRATA: [&str; 3] = [
    "multi_allelic_snv",
    "multi_allelic_indel",
    "multi_allelic_indel_trimmable",
];

/// The length of a text in characters: `len(s)`.
pub fn length(text: &str) -> usize {
    text.chars().count()
}

/// The first character of a text as a text, empty for an empty one: `s[:1]`.
pub fn first(text: &str) -> &str {
    &text[..text.chars().next().map_or(0, char::len_utf8)]
}

/// A text without its first character: `s[1:]`.
pub fn rest(text: &str) -> &str {
    &text[first(text).len()..]
}

/// The last character of a text as a text, empty for an empty one: `s[-1:]`.
pub fn last(text: &str) -> &str {
    &text[text.len() - text.chars().next_back().map_or(0, char::len_utf8)..]
}

/// A text without its last character: `s[:-1]`.
pub fn without_last(text: &str) -> &str {
    &text[..text.len() - last(text).len()]
}

/// `text or "-"`.
pub fn or_dash(text: &str) -> &str {
    if text.is_empty() {
        "-"
    } else {
        text
    }
}

/// Whether an ALT allele is symbolic or a bracket breakend: it starts with `<` or holds a
/// square bracket.
pub fn is_symbolic_or_breakend(alt: &str) -> bool {
    alt.starts_with('<') || alt.contains('[') || alt.contains(']')
}

/// Whether the alleles that hold no `*` share one first character, an empty allele's being
/// the empty text.
pub fn share_first<'a>(alleles: impl IntoIterator<Item = &'a str>) -> bool {
    let firsts: HashSet<&str> = alleles
        .into_iter()
        .filter(|allele| !allele.contains('*'))
        .map(first)
        .collect();
    firsts.len() == 1
}

/// The names Ensembl VEP writes in Uploaded_variation for a VCF record without an ID, in
/// the forms of releases 116 and 115, each once.
///
/// Release 116 names the record from the raw line, `CHR_POS_REF/ALT1/ALT2`, and that name
/// comes first. Release 115's names follow: the alleles are the raw REF and ALTs joined by
/// `/` and the start is the variant's after parsing. It is POS when no allele differs in
/// length from REF. For a bi-allelic indel the first shared base is chopped and the start
/// then moves over the shared leading bases. For a multi-allelic indel the start is POS+1
/// when every allele without `*` shares its first base. A record with a symbolic or
/// breakend allele is named `CHR_START_<ALT...>` with start POS and POS+1, both given.
pub fn vep_auto_names<S: AsRef<str>>(
    chrom: &str,
    pos: i64,
    ref_allele: &str,
    alts: &[S],
) -> Vec<String> {
    let alts: Vec<&str> = alts.iter().map(AsRef::as_ref).collect();
    let pos = i128::from(pos);
    let mut alleles = vec![ref_allele];
    alleles.extend(&alts);
    let raw = alleles.join("/");
    let mut names = vec![format!("{chrom}_{pos}_{raw}")];
    if alts.iter().any(|alt| is_symbolic_or_breakend(alt)) {
        let joined = alts.join("/");
        names.push(format!("{chrom}_{pos}_{joined}"));
        names.push(format!("{chrom}_{}_{joined}", pos + 1));
    } else if alts.len() == 1 {
        let (mut r, mut a, mut start) = (ref_allele, alts[0], pos);
        if length(r) != length(a) {
            if first(r) == first(a) {
                (r, a, start) = (or_dash(rest(r)), or_dash(rest(a)), start + 1);
            }
            while !r.is_empty() && !a.is_empty() && first(r) == first(a) {
                (r, a, start) = (rest(r), rest(a), start + 1);
            }
            names.push(format!("{chrom}_{start}_{raw}"));
        }
    } else if alts.iter().any(|alt| length(alt) != length(ref_allele))
        && share_first(alleles.iter().copied())
    {
        names.push(format!("{chrom}_{}_{raw}", pos + 1));
    }
    without_repeats(names)
}

/// The stratum of a multi-allelic record: `multi_allelic_snv` when no ALT differs in length
/// from REF; else `multi_allelic_indel_trimmable` when, after the first base every allele
/// shares is chopped, REF and some ALT of another length still share their first or their
/// last base, the shape on which a per-allele minimiser and Ensembl VEP's untouched
/// multi-allelic record disagree on positions; else `multi_allelic_indel`. An empty allele
/// among alleles that do not share their first base is refused.
pub fn multi_allelic_stratum<S: AsRef<str>>(
    ref_allele: &str,
    alts: &[S],
) -> Result<&'static str, Failure> {
    let alts: Vec<&str> = alts.iter().map(AsRef::as_ref).collect();
    if !alts.iter().any(|alt| length(alt) != length(ref_allele)) {
        return Ok(MULTI_ALLELIC_STRATA[0]);
    }
    let mut alleles = vec![ref_allele];
    alleles.extend(&alts);
    let firsts: HashSet<&str> = alleles.iter().map(|allele| first(allele)).collect();
    if firsts.len() == 1 {
        for allele in &mut alleles {
            *allele = or_dash(rest(allele));
        }
    }
    let r = alleles[0];
    for a in &alleles[1..] {
        if r == "-" || *a == "-" || length(r) == length(a) {
            continue;
        }
        if r.is_empty() || a.is_empty() {
            return Err(Failure::own(format!(
                "a record with REF {} and ALT {} has an empty allele",
                py_repr_str(ref_allele),
                py_repr_str(&alts.join(","))
            )));
        }
        if first(r) == first(a) || last(r) == last(a) {
            return Ok(MULTI_ALLELIC_STRATA[2]);
        }
    }
    Ok(MULTI_ALLELIC_STRATA[1])
}

/// The POS of a record as `int()` reads it within 64 bits; another text is refused,
/// naming the record by its first three columns.
pub fn record_pos(chrom: &str, pos: &str, id: &str) -> Result<i64, Failure> {
    py_int(pos).ok_or_else(|| {
        Failure::own(format!(
            "the record {chrom} {pos} {id} has a POS that is no integer of 64 bits"
        ))
    })
}

/// `record_names` of a record given as its first five columns.
pub fn names_of_record(
    chrom: &str,
    pos: &str,
    id: &str,
    ref_allele: &str,
    alt: &str,
) -> Result<Vec<String>, Failure> {
    let pos = record_pos(chrom, pos, id)?;
    let first_id = id.split(';').next().unwrap_or("");
    if ![".", "", "0"].contains(&first_id) {
        return Ok(vec![first_id.to_string()]);
    }
    let alts: Vec<&str> = alt.split(',').collect();
    let mut names = vep_auto_names(chrom, pos, ref_allele, &alts);
    if first_id != "." && !alts.iter().any(|alt| is_symbolic_or_breakend(alt)) {
        let mut start = i128::from(pos);
        let mut r = ref_allele;
        let mut parsed: Vec<&str> = alts.clone();
        if alts.len() == 1 && length(r) != length(alts[0]) {
            let mut a = alts[0];
            if first(r) == first(a) {
                (r, a, start) = (or_dash(rest(r)), or_dash(rest(a)), start + 1);
            }
            if r != "-" && a != "-" {
                while !r.is_empty() && !a.is_empty() && first(r) == first(a) {
                    (r, a, start) = (rest(r), rest(a), start + 1);
                }
                while !r.is_empty() && !a.is_empty() && last(r) == last(a) {
                    (r, a) = (without_last(r), without_last(a));
                }
                (r, a) = (or_dash(r), or_dash(a));
            }
            parsed = vec![a];
        } else if alts.len() > 1
            && alts.iter().any(|alt| length(alt) != length(r))
            && share_first(std::iter::once(r).chain(alts.iter().copied()))
        {
            r = or_dash(rest(r));
            for alt in &mut parsed {
                if !alt.contains('*') {
                    *alt = or_dash(rest(alt));
                }
            }
            start += 1;
        }
        names.push(format!("{chrom}_{start}_{r}/{}", parsed.join("/")));
    }
    Ok(without_repeats(names))
}

/// The Uploaded_variation names a VCF record line can carry: the first `;`-separated
/// element of its ID, else, when that element is `.`, empty or `0`, the names Ensembl VEP
/// writes for a record without an ID (`vep_auto_names`; Parser/VCF.pm:347 takes the first
/// element, Parser.pm:462 replaces a name Perl reads as false, OutputFactory.pm:900 a `.`).
/// An empty or `0` element also offers the name Parser.pm:462 builds, `CHR_START_ALLELES`
/// from the alleles the parser left (the anchor chop of Parser/VCF.pm:300-336, then
/// `trim_sequences` on a bi-allelic indel) in the case written, which the output keeps when
/// upper-casing changed the string (OutputFactory.pm:900). A line of fewer than five
/// columns and a POS that is no integer are refused.
pub fn record_names(vcf_line: &str) -> Result<Vec<String>, Failure> {
    let mut columns = vcf_line.split('\t');
    match (
        columns.next(),
        columns.next(),
        columns.next(),
        columns.next(),
        columns.next(),
    ) {
        (Some(chrom), Some(pos), Some(id), Some(ref_allele), Some(alt)) => {
            names_of_record(chrom, pos, id, ref_allele, alt)
        }
        _ => Err(Failure::own(format!(
            "the record line {} has fewer than five columns",
            py_repr_str(vcf_line)
        ))),
    }
}

/// Uploaded_variation (the ID, else a name Ensembl VEP writes) to the ordinals of the
/// record lines it names, in order.
pub fn record_name_index<'a>(
    vcf_lines: impl IntoIterator<Item = &'a str>,
) -> Result<HashMap<String, Vec<usize>>, Failure> {
    let mut index: HashMap<String, Vec<usize>> = HashMap::new();
    for (ordinal, line) in vcf_lines.into_iter().enumerate() {
        for name in record_names(line)? {
            index.entry(name).or_default().push(ordinal);
        }
    }
    Ok(index)
}

#[cfg(test)]
#[path = "names_tests.rs"]
mod tests;
