// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What is compared of a row: the Location, the Allele and the consequence set as both engines are known to spell them.

use std::borrow::Cow;

use crate::input::{push_canonical_location, KeyPrep};
use crate::locus::{strip_chr, CNV_TR_RUST_ALLELE};
use crate::pyfmt::py_strip;
use crate::text::INT_MAX_STR_DIGITS;

/// The `--scored-engine` value that turns the symbolic-allele mapping on, compared without
/// regard to case or surrounding blanks; every other engine name leaves the Allele column
/// as written.
pub const FASTVEP_ENGINE: &str = "fastvep";

/// Whether `engine` names fastVEP, the one engine whose Allele column is mapped: through
/// v0.5.0 fastVEP writes raw symbolic tokens there.
pub fn is_fastvep(engine: &str) -> bool {
    py_strip(engine).eq_ignore_ascii_case(FASTVEP_ENGINE)
}

/// Appends the Location as it is compared: each coordinate written as its number, a
/// leading `chr` in any case dropped from the chromosome, and `M` named `MT`. A text
/// without `:` is appended as written.
pub fn push_normalized_location(loc: &str, out: &mut Vec<u8>) {
    let Some((chrom, _)) = loc.split_once(':') else {
        out.extend_from_slice(loc.as_bytes());
        return;
    };
    let bare = strip_chr(chrom);
    let bare = if bare.eq_ignore_ascii_case("m") {
        "MT"
    } else {
        bare
    };
    out.extend_from_slice(bare.as_bytes());
    push_canonical_location(&loc[chrom.len()..], out);
}

/// The Location as it is compared; see [`push_normalized_location`].
pub fn normalize_location(loc: &str) -> String {
    let mut out = Vec::with_capacity(loc.len());
    push_normalized_location(loc, &mut out);
    String::from_utf8(out).expect("a normalised Location is cut at ASCII bytes of a str")
}

/// The digits of the coordinate Python's integer conversion refuses while it normalises
/// `loc`: its first coordinate written in more than [`INT_MAX_STR_DIGITS`] digits.
pub fn location_integer_digits(loc: &str) -> Option<usize> {
    if loc.len() <= INT_MAX_STR_DIGITS {
        return None;
    }
    let (_, coords) = loc.split_once(':')?;
    coords
        .split('-')
        .map(|coordinate| coordinate.strip_prefix('+').unwrap_or(coordinate))
        .find(|digits| {
            digits.len() > INT_MAX_STR_DIGITS && digits.bytes().all(|b| b.is_ascii_digit())
        })
        .map(str::len)
}

/// The allele without the `chr` that follows a bracket of breakend notation: Ensembl VEP
/// and vep-rs both write the mate's chromosome there as the record spells it, and the
/// prefix is dropped so that an allele compares alike under either spelling.
pub fn strip_chr_in_bracket(allele: &str) -> Cow<'_, str> {
    let bytes = allele.as_bytes();
    let mut out: Option<String> = None;
    let mut copied = 0;
    let mut i = 0;
    while i < bytes.len() {
        if matches!(bytes[i], b'[' | b']') && bytes[i + 1..].starts_with(b"chr") {
            out.get_or_insert_with(|| String::with_capacity(allele.len()))
                .push_str(&allele[copied..=i]);
            copied = i + 4;
            i = copied;
        } else {
            i += 1;
        }
    }
    match out {
        Some(mut out) => {
            out.push_str(&allele[copied..]);
            Cow::Owned(out)
        }
        None => Cow::Borrowed(allele),
    }
}

/// The name Ensembl VEP writes for a mobile-element subtype, by the upper-cased third
/// segment of `<INS:ME:...>` or `<DEL:ME:...>`; any other subtype is the generic
/// `mobile_element`.
fn mobile_element_allele(subtype: Option<&str>, deletion: bool) -> &'static str {
    match (subtype, deletion) {
        (Some("ALU"), false) => "Alu_insertion",
        (Some("ALU"), true) => "Alu_deletion",
        (Some("L1" | "LINE1"), false) => "LINE1_insertion",
        (Some("L1" | "LINE1"), true) => "LINE1_deletion",
        (Some("SVA"), false) => "SVA_insertion",
        (Some("SVA"), true) => "SVA_deletion",
        (Some("HERV"), false) => "HERV_insertion",
        (Some("HERV"), true) => "HERV_deletion",
        (_, false) => "mobile_element_insertion",
        (_, true) => "mobile_element_deletion",
    }
}

/// The digits of a `<CNn>` or `<CN=n>` token, given the token's upper-cased text between
/// the angle brackets; `None` when what follows `CN` and an optional `=` is not a run of
/// the digits `0` to `9`.
fn copy_number_digits(inner: &str) -> Option<&str> {
    let rest = inner.strip_prefix("CN")?;
    let digits = rest.strip_prefix('=').unwrap_or(rest);
    (!digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())).then_some(digits)
}

/// The Allele Ensembl VEP and vep-rs write for a raw symbolic ALT token, which is what
/// fastVEP through v0.5.0 writes in that column: the Sequence Ontology class of the token,
/// read without regard to case. From v0.6.0 fastVEP writes the class name itself.
///
/// | token | Allele |
/// | --- | --- |
/// | `<INS:ME>`, `<INS:ME:other>` | `mobile_element_insertion` |
/// | `<INS:ME:ALU>`, `:L1`, `:LINE1`, `:SVA`, `:HERV` | `Alu_insertion`, `LINE1_insertion`, `LINE1_insertion`, `SVA_insertion`, `HERV_insertion` |
/// | `<DEL:ME[:x]>` | as above, with `_deletion` |
/// | `<DEL...>` | `deletion` |
/// | `<INS...>` | `insertion` |
/// | `<DUP:TANDEM...>` | `tandem_duplication` |
/// | `<DUP...>` | `duplication` |
/// | `<INV...>` | `inversion` |
/// | `<CNV:TR...>` | `tandem_repeat` |
/// | `<CN0>`, `<CN=0>` | `deletion` |
/// | `<CN2>`, `<CN=2>` | `duplication` |
/// | `<CNn>` for any other n, `<CNV...>`, `<CN...>` | `copy_number_variation` |
/// | `<CPX[:x]>` | `CPX` |
///
/// Returned as written: `<NON_REF>` and `<*>`, which both engines write literally;
/// `<BND>`, whose Allele Ensembl VEP builds from INFO fields no row carries; any other
/// `<...>` token; and every text that is not a `<...>` token, the class names themselves
/// included, so that mapping a mapped value changes nothing.
pub fn fastvep_symbolic_allele_to_so_class(allele: &str) -> &str {
    if allele.len() < 3 || !(allele.starts_with('<') && allele.ends_with('>')) {
        return allele;
    }
    let upper = allele.to_uppercase();
    let inner = &upper[1..upper.len() - 1];
    if matches!(upper.as_str(), "<*>" | "<NON_REF>" | "<BND>") {
        return allele;
    }
    // Mobile elements before the generic deletion and insertion: `<DEL:ME>` starts with
    // `<DEL`.
    if upper.starts_with("<INS:ME") || upper.starts_with("<DEL:ME") {
        let mut parts = inner.split(':');
        let deletion = parts.next() == Some("DEL");
        return mobile_element_allele(parts.nth(1), deletion);
    }
    if upper.starts_with("<DEL") {
        return "deletion";
    }
    if upper.starts_with("<INS") {
        return "insertion";
    }
    if upper.starts_with("<DUP:TANDEM") {
        return "tandem_duplication";
    }
    if upper.starts_with("<DUP") {
        return "duplication";
    }
    if upper.starts_with("<INV") {
        return "inversion";
    }
    if upper.starts_with("<CNV:TR") {
        return CNV_TR_RUST_ALLELE;
    }
    if upper.starts_with("<CN") && !upper.starts_with("<CNV") {
        if let Some(digits) = copy_number_digits(inner) {
            return match digits.trim_start_matches('0') {
                "" => "deletion",
                "2" => "duplication",
                _ => "copy_number_variation",
            };
        }
    }
    if upper.starts_with("<CN") {
        return "copy_number_variation";
    }
    if upper.starts_with("<CPX") {
        // A type with no Sequence Ontology term prints Ensembl VEP's type string, the ALT
        // as written or INFO/SVTYPE (Parser/VCF.pm:469-481); the bare abbreviation is the
        // spelling the record's own alleles list.
        return "CPX";
    }
    allele
}

/// The allele as it is compared when its row was written by fastVEP (`fastvep`) or by any
/// other engine: the `chr` after a breakend bracket dropped for every engine, the symbolic
/// token mapped to its class for fastVEP alone.
pub fn normalize_allele_of(allele: &str, fastvep: bool) -> Cow<'_, str> {
    let stripped = strip_chr_in_bracket(allele);
    if !fastvep {
        return stripped;
    }
    match stripped {
        Cow::Borrowed(text) => Cow::Borrowed(fastvep_symbolic_allele_to_so_class(text)),
        Cow::Owned(text) => {
            let mapped = fastvep_symbolic_allele_to_so_class(&text);
            if std::ptr::eq(mapped, text.as_str()) {
                Cow::Owned(text)
            } else {
                Cow::Owned(mapped.to_string())
            }
        }
    }
}

/// The allele as it is compared; `engine` names the engine that wrote the row. See
/// [`normalize_allele_of`].
pub fn normalize_allele(allele: &str, engine: &str) -> String {
    normalize_allele_of(allele, is_fastvep(engine)).into_owned()
}

/// The digits of the copy number Python's integer conversion refuses while it maps a
/// fastVEP `allele`: a copy-number token of more than [`INT_MAX_STR_DIGITS`] digits.
pub fn fastvep_allele_integer_digits(allele: &str) -> Option<usize> {
    if allele.len() <= INT_MAX_STR_DIGITS {
        return None;
    }
    let stripped = strip_chr_in_bracket(allele);
    if !(stripped.starts_with('<') && stripped.ends_with('>')) {
        return None;
    }
    let upper = stripped.to_uppercase();
    let digits = copy_number_digits(&upper[1..upper.len() - 1])?;
    (digits.len() > INT_MAX_STR_DIGITS).then_some(digits.len())
}

/// The consequence set as it is compared: its terms stripped of blanks, each once, in
/// code-point order, joined by commas.
pub fn normalize_consequence_set(raw: &str) -> String {
    let mut out = Vec::with_capacity(raw.len());
    KeyPrep::default().normalized_set(raw, &mut out);
    String::from_utf8(out).expect("the terms of a str")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The table the mapping implements, token to Allele; the keys are the upper-case
    /// spellings.
    const MAPPING: [(&str, &str); 30] = [
        ("<DEL>", "deletion"),
        ("<DEL:ME>", "mobile_element_deletion"),
        ("<DEL:ME:ALU>", "Alu_deletion"),
        ("<DEL:ME:L1>", "LINE1_deletion"),
        ("<DEL:ME:LINE1>", "LINE1_deletion"),
        ("<DEL:ME:SVA>", "SVA_deletion"),
        ("<DEL:ME:HERV>", "HERV_deletion"),
        ("<INS>", "insertion"),
        ("<INS:ME>", "mobile_element_insertion"),
        ("<INS:ME:ALU>", "Alu_insertion"),
        ("<INS:ME:L1>", "LINE1_insertion"),
        ("<INS:ME:LINE1>", "LINE1_insertion"),
        ("<INS:ME:SVA>", "SVA_insertion"),
        ("<INS:ME:HERV>", "HERV_insertion"),
        ("<INS:ME:LINE>", "mobile_element_insertion"),
        ("<DUP>", "duplication"),
        ("<DUP:TANDEM>", "tandem_duplication"),
        ("<DUP:ISP>", "duplication"),
        ("<INV>", "inversion"),
        ("<CNV>", "copy_number_variation"),
        ("<CNV:TR>", "tandem_repeat"),
        ("<CN0>", "deletion"),
        ("<CN=0>", "deletion"),
        ("<CN1>", "copy_number_variation"),
        ("<CN2>", "duplication"),
        ("<CN=2>", "duplication"),
        ("<CN3>", "copy_number_variation"),
        ("<CN9>", "copy_number_variation"),
        ("<CPX>", "CPX"),
        ("<CPX:DEL>", "CPX"),
    ];

    /// Tokens the mapping leaves as written.
    const IDENTITY: [&str; 9] = [
        "<NON_REF>",
        "<*>",
        "<BND>",
        "<CTX>",
        "<STR>",
        "A",
        "ATAT",
        "-",
        "deletion",
    ];

    fn assert_maps(token: &str, expected: &str) {
        assert_eq!(fastvep_symbolic_allele_to_so_class(token), expected);
        assert_eq!(normalize_allele(token, "fastvep"), expected);
    }

    #[test]
    fn test_chr_prefix_stripped() {
        assert_eq!(normalize_location("chr21:100"), "21:100");
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_chrM_becomes_MT() {
        assert_eq!(normalize_location("chrM:100"), "MT:100");
    }

    #[test]
    fn test_bnd_bracket_chr_stripped() {
        assert_eq!(normalize_allele("N[chr21:1234[", "vep-rs"), "N[21:1234[");
    }

    #[test]
    fn test_consequence_set_sorted_and_deduped() {
        assert_eq!(
            normalize_consequence_set("intron_variant,stop_gained,intron_variant"),
            "intron_variant,stop_gained"
        );
    }

    #[test]
    fn test_consequence_order_does_not_affect_identity() {
        let tuple = |csq: &str| {
            (
                normalize_location("21:1"),
                normalize_allele("A", "vep-rs"),
                format!("E{}", 1),
                "Transcript",
                normalize_consequence_set(csq),
            )
        };
        assert_eq!(tuple("b_term,a_term"), tuple("a_term,b_term"));
    }

    #[test]
    fn a_location_is_spelled_one_way_whatever_the_engine_wrote() {
        let cases = [
            ("CHRm:5", "MT:5"),
            ("m:5", "MT:5"),
            ("M", "M"),
            ("chr:5", ":5"),
            ("chr", "chr"),
            ("chrchr1:5", "chr1:5"),
            ("x:y:z", "x:y:z"),
            (":5", ":5"),
            ("chr1:+005-0010", "1:5-10"),
            ("1:5-", "1:5-"),
            ("1:-5", "1:-5"),
            ("1:--5", "1:--5"),
            ("Chr1:1e3", "1:1e3"),
            ("chrMT:1", "MT:1"),
            ("chrm:1:2", "MT:1:2"),
            ("Khr1:5", "Khr1:5"),
            ("1:\u{ff15}", "1:\u{ff15}"),
            ("\u{e9}:5", "\u{e9}:5"),
            ("21:030425949", "21:30425949"),
            ("21:000", "21:0"),
        ];
        for (written, compared) in cases {
            assert_eq!(normalize_location(written), compared, "{written}");
        }
    }

    #[test]
    fn a_coordinate_past_the_integer_limit_is_named_by_its_digits() {
        let zeros = "0".repeat(INT_MAX_STR_DIGITS + 1);
        assert_eq!(location_integer_digits(&format!("1:5-{zeros}")), Some(4301));
        assert_eq!(location_integer_digits(&format!("1:+{zeros}")), Some(4301));
        assert_eq!(location_integer_digits(&format!("1:{}", &zeros[1..])), None);
        assert_eq!(
            location_integer_digits(&format!("1:+{}", &zeros[1..])),
            None
        );
        assert_eq!(location_integer_digits(&zeros), None);
        assert_eq!(location_integer_digits(&format!("1:{zeros}x")), None);
        assert_eq!(
            fastvep_allele_integer_digits(&format!("<cn={zeros}>")),
            Some(4301)
        );
        assert_eq!(
            fastvep_allele_integer_digits(&format!("<CN{}>", &zeros[1..])),
            None
        );
        assert_eq!(
            fastvep_allele_integer_digits(&format!("<DEL{zeros}>")),
            None
        );
        assert_eq!(fastvep_allele_integer_digits(&zeros), None);
    }

    #[test]
    fn the_chr_after_a_bracket_is_dropped_wherever_it_stands() {
        let cases = [
            ("[chrchr1:5[", "[chr1:5["),
            ("]chr1:5]N", "]1:5]N"),
            ("[CHR1:5[", "[CHR1:5["),
            ("[chr[chr]chr", "[[]"),
            ("chr[chr", "chr["),
            ("N[chr21:5[N]chr3:4]", "N[21:5[N]3:4]"),
            ("N[21:5[", "N[21:5["),
            ("[ch", "[ch"),
            ("", ""),
        ];
        for (written, compared) in cases {
            assert_eq!(strip_chr_in_bracket(written), compared, "{written}");
        }
        assert!(matches!(strip_chr_in_bracket("N[21:5["), Cow::Borrowed(_)));
    }

    #[test]
    fn a_token_is_read_by_the_first_rule_it_meets() {
        let cases = [
            ("<INS:MEX>", "mobile_element_insertion"),
            ("<INS:MEI:ALU>", "Alu_insertion"),
            ("<ins:me:alu:x>", "Alu_insertion"),
            ("<DEL:ME::ALU>", "mobile_element_deletion"),
            ("<CN002>", "duplication"),
            ("<CN=>", "copy_number_variation"),
            ("<CN==2>", "copy_number_variation"),
            ("<CN\u{662}>", "copy_number_variation"),
            ("<CN2 >", "copy_number_variation"),
            ("<>", "<>"),
            ("<a", "<a"),
            ("<\u{df}>", "<\u{df}>"),
            ("<de\u{17f}>", "<de\u{17f}>"),
            ("<CNX2>", "copy_number_variation"),
            ("<CN=02>", "duplication"),
            ("<CN+2>", "copy_number_variation"),
            ("<DUP:TANDEMX>", "tandem_duplication"),
            ("<CPXX>", "CPX"),
            ("<non_ref>", "<non_ref>"),
            ("<bnd>", "<bnd>"),
            ("< DEL>", "< DEL>"),
            ("<in\u{17f}>", "insertion"),
        ];
        for (token, expected) in cases {
            assert_eq!(
                fastvep_symbolic_allele_to_so_class(token),
                expected,
                "{token}"
            );
        }
    }

    #[test]
    fn test_the_perl_reference_is_never_mapped() {
        for engine in ["vep-rs", "perl", "VEP-RS"] {
            for (token, _) in MAPPING {
                assert_eq!(normalize_allele(token, engine), token, "{engine} {token}");
            }
        }
    }

    #[test]
    fn test_maps_del() {
        assert_maps("<DEL>", "deletion");
    }

    #[test]
    fn test_maps_ins() {
        assert_maps("<INS>", "insertion");
    }

    #[test]
    fn test_maps_dup() {
        assert_maps("<DUP>", "duplication");
    }

    #[test]
    fn test_maps_dup_tandem_before_generic_dup() {
        assert_maps("<DUP:TANDEM>", "tandem_duplication");
    }

    #[test]
    fn test_maps_dup_subtype_other_than_tandem_to_duplication() {
        assert_maps("<DUP:ISP>", "duplication");
    }

    #[test]
    fn test_maps_inv() {
        assert_maps("<INV>", "inversion");
    }

    #[test]
    fn test_maps_cnv() {
        assert_maps("<CNV>", "copy_number_variation");
    }

    #[test]
    fn test_maps_cnv_tr_before_generic_cnv() {
        assert_maps("<CNV:TR>", CNV_TR_RUST_ALLELE);
        assert_eq!(CNV_TR_RUST_ALLELE, "tandem_repeat");
    }

    #[test]
    fn test_maps_cn0_to_deletion() {
        assert_maps("<CN0>", "deletion");
        assert_maps("<CN=0>", "deletion");
    }

    #[test]
    fn test_maps_cn2_to_duplication() {
        assert_maps("<CN2>", "duplication");
        assert_maps("<CN=2>", "duplication");
    }

    #[test]
    fn test_maps_other_copy_numbers_to_the_generic_class() {
        for n in [1, 3, 4, 5, 6, 7, 8, 9, 12] {
            assert_maps(&format!("<CN{n}>"), "copy_number_variation");
            assert_maps(&format!("<CN={n}>"), "copy_number_variation");
        }
    }

    #[test]
    fn test_unparseable_copy_number_falls_through_to_the_generic_class() {
        assert_maps("<CNX>", "copy_number_variation");
        assert_maps("<CN>", "copy_number_variation");
        assert_maps("<CN=X>", "copy_number_variation");
    }

    #[test]
    fn test_maps_generic_mobile_element_insertion() {
        assert_maps("<INS:ME>", "mobile_element_insertion");
    }

    #[test]
    fn test_maps_alu_insertion() {
        assert_maps("<INS:ME:ALU>", "Alu_insertion");
    }

    #[test]
    fn test_maps_line1_insertion_from_both_spellings() {
        assert_maps("<INS:ME:L1>", "LINE1_insertion");
        assert_maps("<INS:ME:LINE1>", "LINE1_insertion");
    }

    #[test]
    fn test_maps_sva_insertion() {
        assert_maps("<INS:ME:SVA>", "SVA_insertion");
    }

    #[test]
    fn test_maps_herv_insertion() {
        assert_maps("<INS:ME:HERV>", "HERV_insertion");
    }

    #[test]
    fn test_unknown_mobile_element_subtype_is_the_generic_class() {
        assert_maps("<INS:ME:LINE>", "mobile_element_insertion");
        assert_maps("<INS:ME:UNKNOWN>", "mobile_element_insertion");
    }

    #[test]
    fn test_maps_mobile_element_deletion_before_generic_del() {
        assert_maps("<DEL:ME>", "mobile_element_deletion");
        assert_maps("<DEL:ME:ALU>", "Alu_deletion");
        assert_maps("<DEL:ME:SVA>", "SVA_deletion");
    }

    #[test]
    fn test_maps_cpx_to_the_bracket_stripped_abbreviation() {
        assert_maps("<CPX>", "CPX");
        assert_maps("<CPX:DEL>", "CPX");
    }

    #[test]
    fn test_the_whole_mapping_table_is_pinned() {
        let got: Vec<(&str, &str)> = MAPPING
            .iter()
            .map(|(token, _)| (*token, fastvep_symbolic_allele_to_so_class(token)))
            .collect();
        assert_eq!(got, MAPPING);
    }

    #[test]
    fn test_mapped_values_are_exactly_the_perl_so_terms() {
        let perl_so_terms: std::collections::BTreeSet<&str> = [
            "insertion",
            "mobile_element_insertion",
            "Alu_insertion",
            "HERV_insertion",
            "LINE1_insertion",
            "SVA_insertion",
            "deletion",
            "mobile_element_deletion",
            "Alu_deletion",
            "HERV_deletion",
            "LINE1_deletion",
            "SVA_deletion",
            "tandem_repeat",
            "tandem_duplication",
            "duplication",
            "copy_number_variation",
            "inversion",
            "CPX",
        ]
        .into();
        let mapped: std::collections::BTreeSet<&str> =
            MAPPING.iter().map(|(_, value)| *value).collect();
        assert_eq!(mapped, perl_so_terms);
    }

    #[test]
    fn test_identity_tokens_are_unchanged() {
        for token in IDENTITY {
            assert_eq!(fastvep_symbolic_allele_to_so_class(token), token);
        }
    }

    #[test]
    fn test_non_ref_and_star_are_literal_in_both_engines() {
        assert_maps("<NON_REF>", "<NON_REF>");
        assert_maps("<*>", "<*>");
    }

    #[test]
    fn test_symbolic_bnd_is_not_mapped() {
        assert_maps("<BND>", "<BND>");
    }

    #[test]
    fn test_an_unknown_symbolic_token_is_left_alone_not_guessed() {
        assert_maps("<CTX>", "<CTX>");
        assert_maps("<STR>", "<STR>");
        assert_maps("<FOO:BAR>", "<FOO:BAR>");
    }

    #[test]
    fn test_mapping_is_idempotent_on_already_normalised_values() {
        let values: std::collections::BTreeSet<&str> =
            MAPPING.iter().map(|(_, value)| *value).collect();
        assert_eq!(values.len(), 18);
        for value in values {
            assert_eq!(fastvep_symbolic_allele_to_so_class(value), value);
        }
    }

    #[test]
    fn test_token_matching_is_case_insensitive() {
        assert_maps("<del>", "deletion");
        assert_maps("<Dup:Tandem>", "tandem_duplication");
        assert_maps("<ins:me:alu>", "Alu_insertion");
        assert_maps("<cnv:tr>", "tandem_repeat");
    }

    #[test]
    fn test_engine_name_matching_is_case_insensitive() {
        for engine in ["fastvep", "fastVEP", "FASTVEP", " fastvep "] {
            assert_eq!(normalize_allele("<DEL>", engine), "deletion", "{engine}");
        }
        assert!(is_fastvep("\u{a0}fastvep\u{3000}"));
        assert!(is_fastvep("fastvep\u{1f}"));
        assert!(!is_fastvep("fast vep"));
        assert!(!is_fastvep("fa\u{17f}tvep"));
    }

    #[test]
    fn test_breakend_bracket_normalisation_still_applies_under_fastvep() {
        assert_eq!(normalize_allele("N[chr21:1234[", "fastvep"), "N[21:1234[");
        assert_eq!(normalize_allele("N[chr21:1234[", "vep-rs"), "N[21:1234[");
    }

    #[test]
    fn test_literal_bases_are_untouched_under_fastvep() {
        for allele in ["A", "ACGT", "-", "ATATATATATATATAT"] {
            assert_eq!(normalize_allele(allele, "fastvep"), allele);
        }
    }

    #[test]
    fn a_consequence_set_keeps_each_stripped_term_once_in_order() {
        let cases = [
            ("b,a,,b", "a,b"),
            (" a , b ", "a,b"),
            ("", ""),
            (",", ""),
            ("a\u{a0},\u{a0}a", "a"),
            ("a\u{1c},b", "a,b"),
            ("\u{e9},e,Z,a", "Z,a,e,\u{e9}"),
        ];
        for (raw, expected) in cases {
            assert_eq!(normalize_consequence_set(raw), expected, "{raw:?}");
        }
    }
}
