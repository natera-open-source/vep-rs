// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SPDI notation (`SPDI`) for one allele of a record, Ensembl VEP's
//! `spdi_genomic`, and the GA4GH Variation Representation Specification allele
//! (`GA4GH_VRS`) Ensembl VEP derives from it.
//!
//! Perl citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`).

use std::borrow::Cow;

use crate::codon::reverse_complement;

/// The longest allele a tandem notation expands to; a count asking for more
/// leaves the notation as written, so the allele gets no SPDI.
const MAX_EXPANDED_ALLELE_BYTES: usize = 64 << 20;

/// `expand` (ensembl core Utils/Sequence.pm:116-119 at release/116,
/// `s/(\w*)\((\w+)\)(\d+)/$1.$2 x $3/eg`): every `(unit)count` group of a
/// tandem-notation allele written out, `(AC)3` as `ACACAC`, the text around a
/// group kept; an allele without `(`, and a group whose unit or count is not
/// word characters then digits, come back as they are.
fn expand_tandem(allele: &[u8]) -> Cow<'_, [u8]> {
    if !allele.contains(&b'(') {
        return Cow::Borrowed(allele);
    }
    let word = |b: &u8| b.is_ascii_alphanumeric() || *b == b'_';
    let mut out = Vec::with_capacity(allele.len());
    let mut rest = allele;
    while let Some(open) = rest.iter().position(|b| *b == b'(') {
        out.extend_from_slice(&rest[..open]);
        let after_open = &rest[open + 1..];
        let unit_len = after_open.iter().take_while(|b| word(b)).count();
        let count_start = unit_len + 1;
        let count_len = after_open.get(count_start..).map_or(0, |tail| {
            tail.iter().take_while(|b| b.is_ascii_digit()).count()
        });
        let is_group = unit_len > 0 && after_open.get(unit_len) == Some(&b')') && count_len > 0;
        let expanded = if is_group {
            let unit = &after_open[..unit_len];
            let digits = &after_open[count_start..count_start + count_len];
            repeat_unit(unit, digits, out.len())
        } else {
            None
        };
        match expanded {
            Some(bases) => {
                out.extend_from_slice(&bases);
                rest = &after_open[count_start + count_len..];
            }
            None => {
                out.push(b'(');
                rest = after_open;
            }
        }
    }
    out.extend_from_slice(rest);
    Cow::Owned(out)
}

/// `unit` written out `digits` times; `None` when the count does not fit a
/// `usize` or the allele, `written` bytes long so far, would pass
/// [`MAX_EXPANDED_ALLELE_BYTES`].
fn repeat_unit(unit: &[u8], digits: &[u8], written: usize) -> Option<Vec<u8>> {
    let count: usize = std::str::from_utf8(digits).ok()?.parse().ok()?;
    let total = unit.len().checked_mul(count)?;
    if written.checked_add(total)? > MAX_EXPANDED_ALLELE_BYTES {
        return None;
    }
    Some(unit.repeat(count))
}

/// `spdi_genomic` (VariationFeature.pm 2129) for one alternate allele of a
/// record: `reference_name`, the 0-based position, the deleted bases and the
/// inserted bases joined by `:`. `start..=end` is the record's span and
/// `reference` and `alt` its alleles as the record writes them (`-` for no
/// bases), which `reverse` reads on the reverse strand. Each allele has its
/// tandem notation expanded first (2156, 2170, [`expand_tandem`]), so
/// `(AC)3` is the six bases `ACACAC`. A single-base span is a substitution or a
/// one-base deletion (2185-2191); a reference of `-` is an insertion positioned
/// after `end` (2193-2197); a deletion of more than 20 bases writes their count
/// (2201-2203); any other pair writes both alleles as they are (2207-2210).
/// Alleles on the reverse strand are complemented when the deleted part is all
/// bases (2214-2217): the count of a long deletion is not, and a reference
/// carrying a character outside `ACGTN` leaves both the deleted and the
/// inserted part as written. `None` where Ensembl VEP writes no value: a span
/// before the first base (2144), an expanded alternate allele with a character
/// outside `ACGTN-` (2173), and an expanded reference allele ending in one,
/// where Ensembl VEP stops the run (2159); a reference with such a character
/// before its last is written as it is, as Ensembl VEP writes it.
pub fn spdi_genomic(
    reference_name: &str,
    start: u64,
    end: u64,
    reference: &[u8],
    alt: &[u8],
    reverse: bool,
) -> Option<String> {
    if start < 1 || end < 1 {
        return None;
    }
    let reference = &*expand_tandem(reference);
    let alt = &*expand_tandem(alt);
    let base = |b: &u8| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N');
    let base_or_dash = |b: &u8| *b == b'-' || base(b);
    if !reference.last().is_none_or(base_or_dash) || !alt.iter().all(base_or_dash) {
        return None;
    }
    let mut position = start - 1;
    let (mut deleted, mut inserted): (Vec<u8>, Vec<u8>) = if start == end {
        let inserted = if alt == b"-" {
            Vec::new()
        } else {
            alt.to_vec()
        };
        (reference.to_vec(), inserted)
    } else if reference == b"-" {
        position = end;
        (Vec::new(), alt.to_vec())
    } else if alt == b"-" {
        let deleted = if reference.len() <= 20 {
            reference.to_vec()
        } else {
            reference.len().to_string().into_bytes()
        };
        (deleted, Vec::new())
    } else {
        (reference.to_vec(), alt.to_vec())
    };
    // Only a deleted part that is all bases flips (2214): a count stays, and so
    // does a reference carrying a character outside ACGTN, with its insertion.
    if reverse && deleted.iter().all(base) {
        deleted = reverse_complement(&deleted);
        inserted = reverse_complement(&inserted);
    }
    Some(format!(
        "{reference_name}:{position}:{}:{}",
        String::from_utf8_lossy(&deleted),
        String::from_utf8_lossy(&inserted)
    ))
}

/// The notations of one other alternate allele of a multi-allelic record, keyed
/// by the allele as the record's allele string writes it. Ensembl VEP keeps a
/// record's SPDI and VRS alleles by allele and reads them by the row's `Allele`
/// (ensembl-vep OutputFactory.pm 1396-1400), so a row whose `Allele` a genomic
/// shift rotated onto a sibling allele prints the sibling's notations.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SiblingSpdi {
    /// The sibling allele as the record's allele string writes it.
    pub allele: String,
    /// Its SPDI, on the record's span.
    pub spdi: String,
    /// Its GA4GH VRS allele, when the run computes them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ga4gh_vrs: Option<VrsAllele>,
}

/// The GA4GH VRS `Allele` Ensembl VEP builds from an SPDI
/// (`ga4gh_vrs_from_spdi`, Utils/Sequence.pm 1430): a `SequenceLocation` on
/// `refseq:<sequence_id>` over the `SimpleInterval` `start..end` and a
/// `SequenceState` holding the inserted `sequence`; no digest is computed.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct VrsAllele {
    /// The SPDI's sequence name, an `NC_` accession.
    pub sequence_id: String,
    /// The SPDI's 0-based position.
    pub start: u64,
    /// `start` plus the deleted length (1453).
    pub end: u64,
    /// The inserted bases, empty for a deletion.
    pub sequence: String,
}

/// The allele of an SPDI whose sequence name starts with `NC` (1434) and whose
/// four fields parse (1437): the deleted part is a count or the deleted bases,
/// whose length is the count (1448-1452). `None` for any other string, which is
/// what every offline SPDI is, named by its chromosome, and for an end past
/// `u64::MAX`.
pub fn ga4gh_vrs_from_spdi(spdi: &str) -> Option<VrsAllele> {
    if !spdi.starts_with("NC") {
        return None;
    }
    let mut parts = spdi.split(':');
    let (sequence_id, position, deleted, inserted) =
        (parts.next()?, parts.next()?, parts.next()?, parts.next()?);
    if parts.next().is_some()
        || !position.bytes().all(|b| b.is_ascii_digit())
        || !inserted
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
    {
        return None;
    }
    // An empty position fails the parse, as Perl's `\d+` fails to match it.
    let start: u64 = position.parse().ok()?;
    let deleted_is_count = !deleted.is_empty() && deleted.bytes().all(|b| b.is_ascii_digit());
    let deleted_length: u64 = if deleted_is_count {
        deleted.parse().ok()?
    } else {
        deleted.len() as u64
    };
    Some(VrsAllele {
        sequence_id: sequence_id.to_string(),
        start,
        end: start.checked_add(deleted_length)?,
        sequence: inserted.to_string(),
    })
}

impl VrsAllele {
    /// The object as the text output formats carry it, JSON-encoded on one line
    /// in the order the Perl literal names its keys (1455-1470). Ensembl VEP
    /// encodes the hash with its keys in hash order, so a reader compares the
    /// parsed object, not the text. `start` is a JSON string and `end` a number:
    /// Perl's JSON encoder writes a scalar last used as text as a string, and
    /// the position is the SPDI's text while the end is a sum (1446, 1453).
    pub fn encoded(&self) -> String {
        format!(
            concat!(
                r#"{{"location":{{"type":"SequenceLocation","interval":"#,
                r#"{{"type":"SimpleInterval","start":"{start}","end":{end}}},"#,
                r#""sequence_id":"refseq:{id}"}},"#,
                r#""state":{{"type":"SequenceState","sequence":"{seq}"}},"type":"Allele"}}"#,
            ),
            start = self.start,
            end = self.end,
            id = json_text(&self.sequence_id),
            seq = json_text(&self.sequence),
        )
    }
}

/// `s` as the body of a JSON string: the quote, the backslash and control
/// characters escaped.
fn json_text(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ensembl VEP 116.2 under `--spdi` on these records (the chromosome named
    /// `21` offline): a substitution at a 0-based position (`21 10538885 A C`),
    /// an insertion positioned after the anchor base (`21 10578406 T TA`, the
    /// record's `-/A` at 10578407-10578406), a one-base deletion
    /// (`21 20000325 AT A`, `T/-` at 20000326), a 25-base deletion written as
    /// its count (`21 30000000 GCCCACTTTTACTTACTCACGTCTTT G`, the 25 bases behind
    /// the anchor at 30000001-30000025), and a
    /// multi-allelic pair with no shared anchor written as the record's alleles
    /// (`21 14144642 AT CT,A`, `AT/CT/A`). Alleles given on the reverse strand
    /// are complemented (`21 10538885 10538885 T/G -` prints `A:C`,
    /// `21 20003758 20003757 -/A -` prints `::T`), a long deletion's count is
    /// not. A `*` or symbolic allele gets no value; a reference whose last
    /// character is outside `ACGTN-` gets none either, while `21 20006769 RA CA`
    /// prints `21:20006768:RA:CA` as Ensembl VEP does, on either strand.
    #[test]
    fn spdi_snv_ins_del_long_del_and_reverse_strand_match_vep() {
        let spdi = |start, end, reference: &str, alt: &str, reverse| {
            spdi_genomic(
                "21",
                start,
                end,
                reference.as_bytes(),
                alt.as_bytes(),
                reverse,
            )
        };
        assert_eq!(
            spdi(10538885, 10538885, "A", "C", false).as_deref(),
            Some("21:10538884:A:C")
        );
        assert_eq!(
            spdi(10578407, 10578406, "-", "A", false).as_deref(),
            Some("21:10578406::A")
        );
        assert_eq!(
            spdi(10578407, 10578406, "-", "AA", false).as_deref(),
            Some("21:10578406::AA")
        );
        assert_eq!(
            spdi(20000326, 20000326, "T", "-", false).as_deref(),
            Some("21:20000325:T:")
        );
        assert_eq!(
            spdi(30000001, 30000025, "CCCACTTTTACTTACTCACGTCTTT", "-", false).as_deref(),
            Some("21:30000000:25:")
        );
        assert_eq!(
            spdi(14144642, 14144643, "AT", "CT", false).as_deref(),
            Some("21:14144641:AT:CT")
        );
        assert_eq!(
            spdi(14144642, 14144643, "AT", "A", false).as_deref(),
            Some("21:14144641:AT:A")
        );
        assert_eq!(
            spdi(10538885, 10538885, "T", "G", true).as_deref(),
            Some("21:10538884:A:C")
        );
        assert_eq!(
            spdi(20003758, 20003757, "-", "A", true).as_deref(),
            Some("21:20003757::T")
        );
        assert_eq!(
            spdi(20000326, 20000326, "A", "-", true).as_deref(),
            Some("21:20000325:T:")
        );
        assert_eq!(
            spdi(30000001, 30000025, "AAAGACGTGAGTAAGTAAAAGTGGG", "-", true).as_deref(),
            Some("21:30000000:25:")
        );
        // Twenty deleted bases are still written out; the count starts at 21.
        assert_eq!(
            spdi(101, 120, "ACGTACGTACGTACGTACGT", "-", false).as_deref(),
            Some("21:100:ACGTACGTACGTACGTACGT:")
        );
        assert_eq!(
            spdi(101, 121, "ACGTACGTACGTACGTACGTA", "-", false).as_deref(),
            Some("21:100:21:")
        );
        assert_eq!(spdi(100, 100, "A", "*", false), None);
        assert_eq!(spdi(100, 100, "A", "<DEL>", false), None);
        assert_eq!(spdi(20006769, 20006770, "AR", "CA", false), None);
        assert_eq!(
            spdi(20006769, 20006770, "RA", "CA", false).as_deref(),
            Some("21:20006768:RA:CA")
        );
        assert_eq!(
            spdi(20006769, 20006770, "RA", "CA", true).as_deref(),
            Some("21:20006768:RA:CA")
        );
        assert_eq!(spdi(0, 0, "A", "C", false), None);
        assert_eq!(
            spdi_genomic("NC_000021.9", 10538885, 10538885, b"A", b"C", false).as_deref(),
            Some("NC_000021.9:10538884:A:C")
        );
    }

    /// Ensembl VEP 116.2 under `--spdi` on the Ensembl-format record
    /// `21 10538885 10538890 (AC)3/- +` prints `21:10538884:ACACAC:` on every
    /// row: `spdi_genomic` expands a tandem notation (`expand`, Utils/Sequence.pm
    /// 116-119) before it tests the characters (VariationFeature.pm 2156, 2170).
    /// The alternate allele expands the same way, the text around a group is
    /// kept, a reverse-strand allele is complemented after the expansion, and a
    /// parenthesis that opens no `(unit)count` group is left as written: it
    /// fails the character test where that test reaches it (the alternate
    /// allele, the reference's last character) and is printed where it does
    /// not, as Ensembl VEP prints it.
    #[test]
    fn spdi_expands_a_tandem_notation_before_the_character_test() {
        let spdi = |start, end, reference: &str, alt: &str, reverse| {
            spdi_genomic(
                "21",
                start,
                end,
                reference.as_bytes(),
                alt.as_bytes(),
                reverse,
            )
        };
        assert_eq!(
            spdi(10538885, 10538890, "(AC)3", "-", false).as_deref(),
            Some("21:10538884:ACACAC:")
        );
        assert_eq!(
            spdi(10538885, 10538890, "(AC)3", "-", true).as_deref(),
            Some("21:10538884:GTGTGT:")
        );
        assert_eq!(
            spdi(10538885, 10538885, "A", "(AC)3", false).as_deref(),
            Some("21:10538884:A:ACACAC")
        );
        assert_eq!(
            spdi(10538885, 10538892, "A(AC)3G", "T(GA)2", false).as_deref(),
            Some("21:10538884:AACACACG:TGAGA")
        );
        assert_eq!(
            spdi(10538885, 10538894, "(AC)3(GT)2", "-", false).as_deref(),
            Some("21:10538884:ACACACGTGT:")
        );
        assert_eq!(spdi(10538885, 10538886, "(AC)", "-", false), None);
        assert_eq!(
            spdi(10538885, 10538886, "(AC", "-", false).as_deref(),
            Some("21:10538884:(AC:")
        );
        assert_eq!(spdi(10538885, 10538886, "AC(", "-", false), None);
        assert_eq!(spdi(10538885, 10538885, "A", "(AC", false), None);
        assert_eq!(
            spdi(10538885, 10538885, "A", "(AC)99999999999", false),
            None
        );
        assert_eq!(
            &*expand_tandem(b"(AC)3"),
            b"ACACAC",
            "the expansion alone, as `expand` rewrites the string"
        );
        assert_eq!(&*expand_tandem(b"AC(G)0T"), b"ACT");
        assert_eq!(&*expand_tandem(b"ACGT"), b"ACGT");
    }

    /// Ensembl VEP 116.2 in database mode under `--spdi --ga4gh_vrs` for
    /// `21 10538885 A C` (SPDI `NC_000021.9:10538884:A:C`) writes the object
    /// `{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval",
    /// "start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},
    /// "state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}` with
    /// its keys in hash order; a deletion's `sequence` is empty and its `end`
    /// is the start plus the deleted length, a count as much as bases. An SPDI
    /// named by a chromosome, one with fewer than four fields or a position
    /// that is not all digits, and one whose end passes `u64::MAX` yield
    /// nothing.
    #[test]
    fn ga4gh_vrs_object_matches_vep_shape() {
        let snv = ga4gh_vrs_from_spdi("NC_000021.9:10538884:A:C").unwrap();
        assert_eq!(
            snv,
            VrsAllele {
                sequence_id: "NC_000021.9".to_string(),
                start: 10538884,
                end: 10538885,
                sequence: "C".to_string(),
            }
        );
        assert_eq!(
            snv.encoded(),
            concat!(
                r#"{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval","#,
                r#""start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},"#,
                r#""state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}"#
            )
        );
        let deletion = ga4gh_vrs_from_spdi("NC_000021.9:20000325:T:").unwrap();
        assert_eq!((deletion.start, deletion.end), (20000325, 20000326));
        assert_eq!(deletion.sequence, "");
        assert_eq!(
            deletion.encoded(),
            concat!(
                r#"{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval","#,
                r#""start":"20000325","end":20000326},"sequence_id":"refseq:NC_000021.9"},"#,
                r#""state":{"type":"SequenceState","sequence":""},"type":"Allele"}"#
            )
        );
        let counted = ga4gh_vrs_from_spdi("NC_000021.9:30000000:25:").unwrap();
        assert_eq!((counted.start, counted.end), (30000000, 30000025));
        let insertion = ga4gh_vrs_from_spdi("NC_000021.9:10578406::AA").unwrap();
        assert_eq!((insertion.start, insertion.end), (10578406, 10578406));
        assert_eq!(insertion.sequence, "AA");
        assert_eq!(ga4gh_vrs_from_spdi("21:10538884:A:C"), None);
        assert_eq!(ga4gh_vrs_from_spdi("NC_000021.9:x:A:C"), None);
        assert_eq!(ga4gh_vrs_from_spdi("NC_000021.9::A:C"), None);
        assert_eq!(ga4gh_vrs_from_spdi("NC_000021.9:+5:A:C"), None);
        assert_eq!(
            ga4gh_vrs_from_spdi("NC_000021.9:18446744073709551615:1:"),
            None
        );
        assert_eq!(ga4gh_vrs_from_spdi("NC_000021.9:1:A"), None);
        assert_eq!(ga4gh_vrs_from_spdi(""), None);
    }
}
