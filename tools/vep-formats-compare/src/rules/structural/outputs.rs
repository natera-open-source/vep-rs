// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The readers of a Location, an allele and an INFO column that the structural-variant
//! classes share with the tuple comparison (`tools/vep-outputs-compare`).

use super::number::{whole_number, Int};
use crate::pyfmt;
use crate::shared::{_is_non_ref_block, _names_an_allele, perl_fields};
use std::sync::LazyLock;

/// `_MAX_SV_SIZE`: the largest span Ensembl VEP annotates by default.
pub const _MAX_SV_SIZE: Int = 10_000_000;

static LOCATION_SPAN: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"^([^:]+):(-?[0-9]+)(?:-(-?[0-9]+))?\z").expect("the pattern compiles")
});

static INTEGER_COORDINATES: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"^-?[0-9]+(?:--?[0-9]+)?\z").expect("the pattern compiles")
});

/// The bracket mate of an allele as `tools/vep-outputs-compare` reads it: the mate position
/// in the digits `0` to `9` alone.
static ASCII_BRACKET_MATE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]").expect("the pattern compiles")
});

/// `_OWN_BREAKEND_RE`; a `$` of `re` also stands before a last line end.
static OWN_BREAKEND: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"^(?:\.[ACGTN]+|[ACGTN]+\.|chromosome_breakpoint)\n?\z")
        .expect("the pattern compiles")
});

/// `chrom:start-end` or `chrom:pos` as (chrom, start, end), a coordinate the reference
/// numified below zero included (`21:-4-295`); None when unparseable.
pub fn _location_span(loc: &str) -> Option<(String, Int, Int)> {
    let found = LOCATION_SPAN.captures(loc)?;
    let start = whole_number(&found[2]);
    let end = found.get(3).map_or(start, |end| whole_number(end.as_str()));
    Some((found[1].to_string(), start, end))
}

/// Whether a coordinate of the Location is not an integer (`21:30425651-30425949.5`,
/// `21:30425651-NaN`, `21:30425949.0`); a coordinate the reference numified below zero
/// (`21:-4-295`) is an integer.
pub fn _non_integer_end(location: &str) -> bool {
    match location.split_once(':') {
        Some((_chrom, coords)) => !INTEGER_COORDINATES.is_match(coords),
        None => false,
    }
}

/// Whether the allele column names a breakend record's own end rather than a mate.
pub fn is_own_breakend_allele(allele: &str) -> bool {
    OWN_BREAKEND.is_match(allele)
}

/// Whether the allele column names a breakend's own end or its mate.
pub fn is_breakend_allele(allele: &str) -> bool {
    is_own_breakend_allele(allele)
        || ASCII_BRACKET_MATE.captures(allele).is_some_and(|mate| {
            // `bracket_mate` converts the position it then leaves unread.
            whole_number(&mate[2]);
            true
        })
}

/// The Locations at which vep-rs can write the point `c:(P+1)` of a dot-notation breakend
/// without SVTYPE whose reference row spans `span` (chromosome, low, high): the span's
/// start and the base after it.
pub fn dot_breakend_points(span: &(String, Int, Int)) -> [String; 2] {
    let (chrom, start, _end) = span;
    [
        format!("{chrom}:{start}"),
        format!("{chrom}:{}", start.saturating_add(1)),
    ]
}

/// Whether Perl reads an INFO value as true: defined, not empty and not `0`.
pub fn perl_true(value: Option<&str>) -> bool {
    !matches!(value, None | Some("") | Some("0"))
}

/// The value of `key=value` in an INFO column, Perl's `defined $info->{key}`: `None` for
/// an absent key or a bare flag, the last occurrence winning, and the text up to a second
/// `=` as the value (`END2=30425700=5` is `30425700`, `SVTYPE==DEL` the empty value).
pub fn _info_written(info: &str, key: &str) -> Option<String> {
    let mut value = None;
    for item in info.split(';') {
        match item.split_once('=') {
            Some((name, text)) if name == key => {
                value = Some(text.split('=').next().unwrap_or(text).to_string());
            }
            None if item == key => value = None,
            _ => {}
        }
    }
    value
}

/// A character of `\w` of `re` on text: a letter, a digit of any script, a numeric
/// character or the underscore.
fn is_word(c: char) -> bool {
    c == '_'
        || pyfmt::char_is_alpha(c)
        || pyfmt::char_is_decimal(c)
        || pyfmt::char_is_digit(c)
        || pyfmt::char_is_numeric(c)
}

/// `_STRUCTURAL_ALT_LIST.search(text)`, the ALT forms that route a line to the structural
/// path (Parser/VCF.pm:241): `[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$`.
pub fn structural_alt_list(text: &str) -> bool {
    let chars: Vec<char> = text.chars().collect();
    // An opening bracket, one character or more that is no `*`, a closing bracket.
    let mut run_end = chars.len();
    for at in (0..chars.len()).rev() {
        if chars[at] == '*' {
            run_end = at;
            continue;
        }
        if matches!(chars[at], '<' | '[' | ']')
            && chars[(at + 2).min(run_end)..run_end]
                .iter()
                .any(|c| matches!(c, '>' | ']' | '['))
        {
            return true;
        }
    }
    if chars.first() == Some(&'.') && chars.get(1).is_some_and(|c| is_word(*c)) {
        return true;
    }
    // `$` also stands before a last line end.
    let body = match chars.split_last() {
        Some((&'\n', body)) if body.last() == Some(&'.') => body,
        _ => &chars[..],
    };
    match body {
        [.., before, '.'] => is_word(*before),
        _ => false,
    }
}

/// Whether a record is a gVCF variant site vep-rs reads on its sequence path: `<NON_REF>`
/// beside alleles that name a variant, no true INFO/SVTYPE, and no remaining allele written
/// in a structural form (Parser/VCF.pm:237-245). A trailing empty field of the ALT list is
/// no field (`perl_fields`).
pub fn is_gvcf_variant_site(alts: &[String], svtype: Option<&str>) -> bool {
    let alts = perl_fields(alts);
    let real: Vec<&str> = alts
        .iter()
        .map(String::as_str)
        .filter(|alt| _names_an_allele(alt))
        .collect();
    if real.is_empty() || !alts.iter().any(|alt| _is_non_ref_block(alt)) || perl_true(svtype) {
        return false;
    }
    !structural_alt_list(&real.join(","))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_location_reads_as_its_span() {
        assert_eq!(_location_span("21:100"), Some(("21".to_string(), 100, 100)));
        assert_eq!(
            _location_span("21:-4-295"),
            Some(("21".to_string(), -4, 295))
        );
        assert_eq!(
            _location_span("chr1:5--6"),
            Some(("chr1".to_string(), 5, -6))
        );
        assert_eq!(_location_span("a\nb:007"), Some(("a\nb".to_string(), 7, 7)));
        for unread in [
            "21", "21:", ":5", "21:5-", "21:5.0", "21:5\n", "21:٣", "21:1:2",
        ] {
            assert_eq!(_location_span(unread), None, "{unread:?}");
        }
    }

    #[test]
    fn a_coordinate_that_is_no_integer_is_told_from_a_negative_one() {
        for integer in ["21:5", "21:-4-295", "21:5--6", "21"] {
            assert!(!_non_integer_end(integer), "{integer}");
        }
        for other in [
            "21:5.0",
            "21:30425651-NaN",
            "21:",
            "21:5-",
            "21:5\n",
            "21:1:2",
        ] {
            assert!(_non_integer_end(other), "{other:?}");
        }
    }

    #[test]
    fn the_structural_alt_forms_route_a_line() {
        for structural in [
            "<DEL>", "N[22:5[", "]22:5]N", ".N", "N.", "N.\n", "G,<DEL>", "<a*>,<b>", "._", "².",
        ] {
            assert!(structural_alt_list(structural), "{structural:?}");
        }
        for sequence in [
            "G", "<>", "<*>", "<*", ".", "..", "-.", ". N", "G,T", "<a*>", "N.\n\n", "<",
        ] {
            assert!(!structural_alt_list(sequence), "{sequence:?}");
        }
    }

    #[test]
    fn an_info_value_is_read_as_the_reference_reads_it() {
        let info = "SVTYPE=DEL;END2=30425700=5;FLAG;SVLEN=;SVTYPE=BND";
        assert_eq!(_info_written(info, "SVTYPE").as_deref(), Some("BND"));
        assert_eq!(_info_written(info, "END2").as_deref(), Some("30425700"));
        assert_eq!(_info_written(info, "SVLEN").as_deref(), Some(""));
        assert_eq!(_info_written(info, "FLAG"), None);
        assert_eq!(_info_written("END=5;END", "END"), None);
        assert_eq!(_info_written("SVTYPE==DEL", "SVTYPE").as_deref(), Some(""));
        assert!(perl_true(Some("DEL")));
        assert!(!perl_true(Some("0")) && !perl_true(Some("")) && !perl_true(None));
    }

    #[test]
    fn a_breakend_allele_names_an_own_end_or_a_mate() {
        for own in ["N.", ".N", "chromosome_breakpoint", ".ACGTN", "N.\n"] {
            assert!(is_own_breakend_allele(own), "{own:?}");
        }
        for other in ["n.", "N", ".", "deletion", "N[22:5[", "N.\n\n"] {
            assert!(!is_own_breakend_allele(other), "{other:?}");
        }
        assert!(is_breakend_allele("N[22:5["));
        assert!(is_breakend_allele("]chr22:5]N"));
        assert!(!is_breakend_allele("N[22:٣["));
        assert!(!is_breakend_allele("N[22["));
    }
}
