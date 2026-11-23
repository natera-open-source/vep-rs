// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The allele trimming of Ensembl VEP's variant matcher.
//!
//! `get_matched_variant_alleles` (ensembl-variation `Utils/Sequence.pm:1144-1247`)
//! decides whether two records name one allele change: it trims the bases a
//! reference and an alternate allele share, in the orders [`trim_directions`]
//! gives, and compares what is left with the position the trim reaches.
//! [`trim_sequences`] is that trim on a forward-strand pair, returning what is
//! left of each allele and how many bases left each end, so a caller moves a
//! start or an end of its own kind by the count.
//!
//! Ensembl VEP compares bases as written (`eq`, `Utils/Sequence.pm:978`, `:991`,
//! `:1006`, `:1019`). A caller that reads a data file's alleles in either letter
//! case passes [`CaseRule::EitherCase`]; the alleles returned are slices of the
//! caller's own, in the case the caller gave them.

/// How two bases are compared while a pair is trimmed and when two trimmed
/// alleles are matched.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaseRule {
    /// Byte for byte, as Ensembl VEP's `eq` compares them.
    AsWritten,
    /// With ASCII letters equal in either case.
    EitherCase,
}

impl CaseRule {
    /// Whether `a` and `b` are one base under this rule.
    pub fn same_base(self, a: u8, b: u8) -> bool {
        match self {
            Self::AsWritten => a == b,
            Self::EitherCase => a.eq_ignore_ascii_case(&b),
        }
    }

    /// Whether `a` and `b` are one allele under this rule.
    pub fn same_allele(self, a: &[u8], b: &[u8]) -> bool {
        match self {
            Self::AsWritten => a == b,
            Self::EitherCase => a.eq_ignore_ascii_case(b),
        }
    }
}

/// An allele pair as [`trim_sequences`] leaves it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Trimmed<'a> {
    /// What is left of the reference allele, `-` when nothing is.
    pub reference: &'a [u8],
    /// What is left of the alternate allele, `-` when nothing is.
    pub alt: &'a [u8],
    /// The bases trimmed from the left, each a `$start++` on the forward
    /// strand (`Utils/Sequence.pm:998`, `:1013`).
    pub left: usize,
    /// The bases trimmed from the right, each an `$end--` on the forward
    /// strand (`Utils/Sequence.pm:985`, `:1026`).
    pub right: usize,
}

/// An empty allele as Ensembl VEP writes it under `$empty_to_dash`
/// (`Utils/Sequence.pm:1032-1035`): `-`, and any other allele as it is.
pub fn dashed(allele: &[u8]) -> &[u8] {
    if allele.is_empty() {
        b"-"
    } else {
        allele
    }
}

/// `_get_trim_directions` (`Utils/Sequence.pm:1250-1258`): the values of
/// `$end_first` a pair is trimmed under, the left-first order alone when
/// neither allele is longer than one byte and the right-first order after it
/// otherwise.
pub fn trim_directions(reference: &[u8], alt: &[u8]) -> &'static [bool] {
    if reference.len() > 1 || alt.len() > 1 {
        &[false, true]
    } else {
        &[false]
    }
}

/// `trim_sequences` (`Utils/Sequence.pm:965-1038`) with `$empty_to_dash` set:
/// the bytes the pair shares leave the left end and then the right, or under
/// `end_first` the right end and then the left, each loop running while both
/// alleles are Perl-true, neither empty nor the one-character zero.
///
/// The bytes are those of the undecoded string, as Perl's `substr` reads them,
/// so a shared character of two bytes counts two.
pub fn trim_sequences<'a>(
    reference: &'a [u8],
    alt: &'a [u8],
    end_first: bool,
    case: CaseRule,
) -> Trimmed<'a> {
    let mut trimmed = Trimmed {
        reference,
        alt,
        left: 0,
        right: 0,
    };
    if end_first {
        trimmed.trim_right(case);
        trimmed.trim_left(case);
    } else {
        trimmed.trim_left(case);
        trimmed.trim_right(case);
    }
    trimmed.reference = dashed(trimmed.reference);
    trimmed.alt = dashed(trimmed.alt);
    trimmed
}

impl Trimmed<'_> {
    fn trim_left(&mut self, case: CaseRule) {
        while perl_true(self.reference)
            && perl_true(self.alt)
            && case.same_base(self.reference[0], self.alt[0])
        {
            self.reference = &self.reference[1..];
            self.alt = &self.alt[1..];
            self.left += 1;
        }
    }

    fn trim_right(&mut self, case: CaseRule) {
        while perl_true(self.reference)
            && perl_true(self.alt)
            && case.same_base(
                self.reference[self.reference.len() - 1],
                self.alt[self.alt.len() - 1],
            )
        {
            self.reference = &self.reference[..self.reference.len() - 1];
            self.alt = &self.alt[..self.alt.len() - 1];
            self.right += 1;
        }
    }
}

/// Perl's truth of an allele while it is trimmed (`$ref && $alt`): the empty
/// string and `0` are false.
fn perl_true(allele: &[u8]) -> bool {
    !allele.is_empty() && allele != b"0"
}

/// The pair without its first base, as the `parse_data` of a plugin whose
/// table writes an insertion or a deletion with its anchor base strips it
/// before the match (AlphaMissense.pm:293-301, gnomADMt.pm:306-314,
/// CADD.pm:279-289, OpenTargets.pm:401-411): when the alleles differ in length
/// and share their first byte under `case`, that byte goes from both, and what
/// is left of an allele is `-` when it is Perl-false, empty or the
/// one-character zero (`substr($ref, 1) || "-"`). `None` for any other pair,
/// which is handed on as it is. The caller steps its own position by one.
pub fn strip_shared_first_base<'a>(
    reference: &'a [u8],
    alt: &'a [u8],
    case: CaseRule,
) -> Option<(&'a [u8], &'a [u8])> {
    let (&first_reference, &first_alt) = (reference.first()?, alt.first()?);
    let rest = |allele: &'a [u8]| -> &'a [u8] {
        if perl_true(&allele[1..]) {
            &allele[1..]
        } else {
            b"-"
        }
    };
    (reference.len() != alt.len() && case.same_base(first_reference, first_alt))
        .then(|| (rest(reference), rest(alt)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `(reference, alt, left, right)` of a trim under the reference's rule.
    fn as_written(reference: &str, alt: &str, end_first: bool) -> (String, String, usize, usize) {
        under(CaseRule::AsWritten, reference, alt, end_first)
    }

    fn under(
        case: CaseRule,
        reference: &str,
        alt: &str,
        end_first: bool,
    ) -> (String, String, usize, usize) {
        let trimmed = trim_sequences(reference.as_bytes(), alt.as_bytes(), end_first, case);
        (
            String::from_utf8_lossy(trimmed.reference).into_owned(),
            String::from_utf8_lossy(trimmed.alt).into_owned(),
            trimmed.left,
            trimmed.right,
        )
    }

    fn pair(
        reference: &str,
        alt: &str,
        left: usize,
        right: usize,
    ) -> (String, String, usize, usize) {
        (reference.to_string(), alt.to_string(), left, right)
    }

    /// `trim_sequences` as `Utils/Sequence.pm:965-1038` trims in each order: an
    /// insertion after a shared base, a deletion inside a repeat, a
    /// substitution between two shared bases, a pair with nothing shared and
    /// identical alleles, which empty to `-/-`.
    #[test]
    fn both_orders_trim_both_ends_and_count_each() {
        assert_eq!(as_written("AT", "ATT", false), pair("-", "T", 2, 0));
        assert_eq!(as_written("AT", "ATT", true), pair("-", "T", 1, 1));
        assert_eq!(as_written("GG", "G", false), pair("G", "-", 1, 0));
        assert_eq!(as_written("GG", "G", true), pair("G", "-", 0, 1));
        assert_eq!(as_written("A", "G", false), pair("A", "G", 0, 0));
        assert_eq!(as_written("AC", "AC", false), pair("-", "-", 2, 0));
        assert_eq!(as_written("AC", "AC", true), pair("-", "-", 0, 2));

        assert_eq!(as_written("ACGT", "AGT", false), pair("C", "-", 1, 2));
        assert_eq!(as_written("ACGT", "AGT", true), pair("C", "-", 1, 2));
        assert_eq!(as_written("ATA", "A", false), pair("TA", "-", 1, 0));
        assert_eq!(as_written("ATA", "A", true), pair("AT", "-", 0, 1));

        assert_eq!(as_written("ATG", "AAG", false), pair("T", "A", 1, 1));
        assert_eq!(as_written("ATG", "AAG", true), pair("T", "A", 1, 1));
        assert_eq!(as_written("AT", "A", false), pair("T", "-", 1, 0));
        assert_eq!(as_written("AT", "A", true), pair("T", "-", 1, 0));
        assert_eq!(as_written("TA", "TATA", false), pair("-", "TA", 2, 0));
        assert_eq!(as_written("TA", "TATA", true), pair("-", "TA", 0, 2));

        assert_eq!(as_written("G", "GT", false), pair("-", "T", 1, 0));
        assert_eq!(as_written("GT", "G", false), pair("T", "-", 1, 0));
        assert_eq!(as_written("C", "CC", false), pair("-", "C", 1, 0));
        assert_eq!(as_written("C", "CC", true), pair("-", "C", 0, 1));
        assert_eq!(as_written("AAA", "A", false), pair("AA", "-", 1, 0));
        assert_eq!(as_written("AAA", "A", true), pair("AA", "-", 0, 1));
        assert_eq!(as_written("ACA", "ACACA", false), pair("-", "CA", 3, 0));
        assert_eq!(as_written("ACA", "ACACA", true), pair("-", "AC", 0, 3));
    }

    /// Each loop runs while both alleles are Perl-true (`$ref && $alt`), so a
    /// pair of the one-character zero is not trimmed and a trim stops at a
    /// zero it leaves alone; an allele already `-` shares nothing.
    #[test]
    fn a_trim_stops_at_a_perl_false_allele() {
        assert_eq!(as_written("0", "0", false), pair("0", "0", 0, 0));
        assert_eq!(as_written("0", "0", true), pair("0", "0", 0, 0));
        assert_eq!(as_written("A0", "A0", false), pair("0", "0", 1, 0));
        assert_eq!(as_written("0A", "0A", true), pair("0", "0", 0, 1));
        assert_eq!(as_written("0", "0T", false), pair("0", "0T", 0, 0));
        assert_eq!(as_written("-", "T", false), pair("-", "T", 0, 0));
        assert_eq!(as_written("", "T", false), pair("-", "T", 0, 0));
    }

    /// Ensembl VEP trims the bytes of an undecoded string, so a two-byte
    /// character shared at the left counts two.
    #[test]
    fn a_shared_two_byte_character_counts_two() {
        assert_eq!(
            as_written("\u{e9}A", "\u{e9}C", false),
            pair("A", "C", 2, 0)
        );
    }

    #[test]
    fn a_one_base_pair_is_trimmed_in_one_order() {
        assert_eq!(trim_directions(b"A", b"G"), &[false]);
        assert_eq!(trim_directions(b"-", b"G"), &[false]);
        assert_eq!(trim_directions(b"AT", b"A"), &[false, true]);
        assert_eq!(trim_directions(b"A", b"AT"), &[false, true]);
    }

    /// As written a base in lower case is another base than its capital, so a
    /// pair written in two cases shares nothing; in either case it shares what
    /// it shares as letters, and what is left keeps the case it was given in.
    #[test]
    fn the_case_rule_decides_which_bases_are_shared() {
        assert_eq!(as_written("AT", "atT", false), pair("A", "at", 0, 1));
        assert_eq!(
            under(CaseRule::EitherCase, "AT", "atT", false),
            pair("-", "T", 2, 0)
        );
        assert_eq!(
            under(CaseRule::EitherCase, "acgt", "AGT", true),
            pair("c", "-", 1, 2)
        );
        assert!(CaseRule::EitherCase.same_allele(b"acgt", b"ACGT"));
        assert!(!CaseRule::AsWritten.same_allele(b"acgt", b"ACGT"));
        assert!(CaseRule::AsWritten.same_allele(b"ACGT", b"ACGT"));
        assert!(!CaseRule::EitherCase.same_allele(b"ACG", b"ACGT"));
        assert!(CaseRule::EitherCase.same_base(b'a', b'A'));
        assert!(!CaseRule::AsWritten.same_base(b'a', b'A'));
    }

    /// The first-base strip: a pair of two lengths that shares its first base
    /// loses it, an allele left empty or as the one-character zero reading `-`;
    /// a pair of one length, a pair with another first base and a pair with an
    /// empty allele are left alone.
    #[test]
    fn a_pair_of_two_lengths_loses_its_shared_first_base() {
        let strip = |reference: &'static str, alt: &'static str, case| {
            strip_shared_first_base(reference.as_bytes(), alt.as_bytes(), case)
        };
        assert_eq!(
            strip("AT", "A", CaseRule::AsWritten),
            Some((&b"T"[..], &b"-"[..]))
        );
        assert_eq!(
            strip("A", "ATG", CaseRule::AsWritten),
            Some((&b"-"[..], &b"TG"[..]))
        );
        assert_eq!(strip("AT", "AG", CaseRule::AsWritten), None);
        assert_eq!(strip("AT", "G", CaseRule::AsWritten), None);
        assert_eq!(strip("", "G", CaseRule::AsWritten), None);
        assert_eq!(strip("AT", "", CaseRule::AsWritten), None);
        assert_eq!(strip("at", "A", CaseRule::AsWritten), None);
        assert_eq!(
            strip("at", "A", CaseRule::EitherCase),
            Some((&b"t"[..], &b"-"[..]))
        );
        assert_eq!(
            strip("A0", "A", CaseRule::AsWritten),
            Some((&b"-"[..], &b"-"[..]))
        );
        assert_eq!(
            strip("A00", "A", CaseRule::AsWritten),
            Some((&b"00"[..], &b"-"[..]))
        );
    }
}
