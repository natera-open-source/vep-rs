// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How Ensembl VEP routes a VCF record to its structural-variant path and names the
//! record's type, shared by every VCF reader.
//!
//! `Parser/VCF.pm` lines 237-245 route a line to `create_StructuralVariationFeatures`
//! when REF and the ALTs together are not a run of `ACGT` and INFO/SVTYPE is true or an
//! ALT is written in a structural form. Lines 469-481 read the type from the ALTs joined
//! by `/`, INFO/SVTYPE standing in for a string outside the VCF 4.4 forms, through
//! `get_SO_term` (`Parser.pm` lines 663-727), whose abbreviation is looked up in
//! `%SO_TERMS` (ensembl-variation `Utils/Config.pm`); a type with no entry is dropped
//! as unsupported.

use crate::variant::VariantClass;

/// Whether Perl reads an INFO value as true: anything but an empty string and `0`, the
/// test `$info->{SVTYPE}` passes at Parser/VCF.pm:240 and :473, so `SVTYPE=`, `SVTYPE=0`
/// and a bare flag name no type.
pub fn perl_true(value: &str) -> bool {
    !value.is_empty() && value != "0"
}

/// The fields of `text` as Perl's list-context `split` returns them (ensembl-io
/// BaseVCF4.pm:489 `[split(',',$alt_allele)]`, Parser.pm:669 `split(/[\/,]/, $type)`):
/// a trailing empty field is no field, a middle one is kept, so `<DEL>,` is one ALT and
/// `G,,T` three.
pub fn perl_split<'a>(text: &'a str, separators: &'a [char]) -> impl Iterator<Item = &'a str> + 'a {
    let kept = text.trim_end_matches(separators);
    (!kept.is_empty())
        .then(|| kept.split(separators))
        .into_iter()
        .flatten()
}

/// The ALT column's fields with their written positions, read as Perl reads the list.
pub fn written_alts(alt_field: &str) -> Vec<(usize, &str)> {
    perl_split(alt_field, &[',']).enumerate().collect()
}

/// Keeps the ALTs a record names its variants by. A gVCF variant site lists `<NON_REF>`
/// beside its real alleles; the block is the absence of a variant, so it names no
/// allele of such a record, which is read as its other ALTs alone, each keeping its
/// written position. A list of `<NON_REF>` alone is a reference block and is kept
/// whole. Ensembl VEP drops both shapes as unsupported types (`G/NON_REF`, `NON_REF`).
/// Under a true INFO/SVTYPE a record carrying the block is structural
/// (Parser/VCF.pm:237-245) and every ALT is an allele of it (`:469`), the block
/// included: Ensembl VEP 116.2 writes two `deletion` rows per transcript for
/// `G,<NON_REF>` with `SVTYPE=DEL`.
pub fn variant_alts(alts: &mut Vec<(usize, &str)>, svtype: Option<&str>) {
    let is_block = |alt: &str| alt.eq_ignore_ascii_case("<NON_REF>");
    let blocks = alts.iter().filter(|(_, alt)| is_block(alt)).count();
    if !svtype.is_some_and(perl_true) && blocks > 0 && blocks < alts.len() {
        alts.retain(|(_, alt)| !is_block(alt));
    }
}

/// Parser/VCF.pm:237-245: the line is structural when REF and the ALTs together are
/// not a run of `ACGT` and INFO/SVTYPE is true or an ALT is written in a structural
/// form, so `N G` with `SVTYPE=DEL` is a deletion, `<*>` with `SVTYPE=DEL` a deletion
/// over its END, and `N G` alone a sequence variant.
pub fn routes_to_structural(
    ref_allele: &str,
    alts: &[(usize, &str)],
    svtype: Option<&str>,
) -> bool {
    let plain_acgt = ref_allele
        .bytes()
        .chain(alts.iter().flat_map(|(_, alt)| alt.bytes()))
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'));
    !plain_acgt
        && (svtype.is_some_and(perl_true) || is_perl_structural_alt_list(&joined(alts, ",")))
}

/// The ALTs joined by `sep`, in their written order.
pub fn joined(alts: &[(usize, &str)], sep: &str) -> String {
    alts.iter()
        .map(|(_, alt)| *alt)
        .collect::<Vec<_>>()
        .join(sep)
}

/// The class of every allele a structural record writes a row for, with the ALT's
/// written position: under a type `%SO_TERMS` names, every ALT that names an allele
/// takes the record's class (one StructuralVariationFeature carries the record,
/// Parser/VCF.pm:560-574); under a type the table lacks, each symbolic ALT keeps its own
/// class under the `vep_skip` mark, INFO/SVTYPE naming the class of a mixed record whose
/// joined ALTs open with a sequence allele (`:472-475`). `valid_alts` excludes `*` and `.`.
pub fn per_alt_classes<'a>(
    valid_alts: &[(usize, &'a str)],
    record_class: Option<VariantClass>,
    svtype: Option<&str>,
) -> Vec<(usize, &'a str, VariantClass)> {
    if let Some(class) = record_class {
        return valid_alts.iter().map(|&(i, alt)| (i, alt, class)).collect();
    }
    let symbolic: Vec<(usize, &str)> = valid_alts
        .iter()
        .copied()
        .filter(|(_, alt)| is_symbolic_alt(alt))
        .collect();
    let mixed = !symbolic.is_empty() && symbolic.len() < valid_alts.len();
    symbolic
        .into_iter()
        .map(|(i, alt)| {
            let class = match svtype {
                Some(svtype) if mixed => {
                    classify_symbolic_alt(&format!("<{svtype}>"), Some(svtype))
                }
                _ => classify_symbolic_alt(alt, svtype),
            };
            (i, alt, class)
        })
        .collect()
}

/// Parser/VCF.pm:241, `[<\[\]][^\*]+[>\]\[]|^\.\w+|\w+\.$` over the ALTs joined by
/// `,`: an angle or bracket token with at least one character and no `*` inside it
/// anywhere in the list, or a dot-bounded breakend as the list's first or last ALT
/// (`G,.N,T` is not routed; `.N,G` is). `<*>` and `*` are not structural forms.
fn is_perl_structural_alt_list(alts: &str) -> bool {
    // The token pattern cannot cross a `*`, so it matches when a `*`-free stretch holds
    // an opening bracket and, at least one character later, a closing one.
    let token = alts.split('*').any(|stretch| {
        let open = stretch.find(['<', '[', ']']);
        let close = stretch.rfind(['>', ']', '[']);
        matches!((open, close), (Some(open), Some(close)) if close >= open + 2)
    });
    let word = |c: char| c.is_ascii_alphanumeric() || c == '_';
    let leading_dot = alts
        .strip_prefix('.')
        .is_some_and(|rest| rest.starts_with(word));
    let trailing_dot = alts
        .strip_suffix('.')
        .is_some_and(|rest| rest.ends_with(word));
    token || leading_dot || trailing_dot
}

fn is_paired_breakend_alt(alt: &str) -> bool {
    alt.contains('[') || alt.contains(']')
}

fn is_native_single_breakend_alt(alt: &str) -> bool {
    alt.len() > 1 && (alt.ends_with('.') || alt.starts_with('.')) && !is_paired_breakend_alt(alt)
}

/// Whether an ALT allele is written as a symbolic or breakend allele.
pub fn is_symbolic_alt(alt: &str) -> bool {
    if alt.is_empty() {
        return false;
    }
    // Perl VEP's SV regex `[<\[\]][^\*]+[>\]\[]` excludes `<*>`: without INFO/SVTYPE
    // it is a non-SV allele at POS.
    if alt.starts_with('<') && alt.ends_with('>') && alt != "<*>" {
        return true;
    }
    is_paired_breakend_alt(alt) || is_native_single_breakend_alt(alt)
}

/// The class of one symbolic ALT allele, with INFO/SVTYPE standing in for an ALT that
/// names no type, for a record whose joined type has no Sequence Ontology term.
pub fn classify_symbolic_alt(alt: &str, svtype: Option<&str>) -> VariantClass {
    let upper = alt.to_uppercase();

    // Mobile elements must be checked before generic DEL/INS: <DEL:ME> starts with <DEL.
    if upper.starts_with("<INS:ME") {
        return VariantClass::MobileElementInsertion;
    }
    if upper.starts_with("<DEL:ME") {
        return VariantClass::MobileElementDeletion;
    }
    if upper.starts_with("<DEL") {
        return VariantClass::StructuralDeletion;
    }
    if upper.starts_with("<INS") {
        return VariantClass::StructuralInsertion;
    }
    if upper.starts_with("<DUP:TANDEM") {
        return VariantClass::TandemDuplication;
    }
    if upper.starts_with("<DUP") {
        return VariantClass::Duplication;
    }
    if upper.starts_with("<INV") {
        return VariantClass::Inversion;
    }
    if upper.starts_with("<CNV:TR") {
        return VariantClass::TandemRepeat;
    }
    // Perl VEP maps CN0 -> DEL ("deletion"), CN2 -> DUP ("duplication")
    if upper.starts_with("<CN") && !upper.starts_with("<CNV") {
        let inner = &upper[3..upper.len().saturating_sub(1)]; // strip "<CN" and ">"
        let num_str = inner.strip_prefix('=').unwrap_or(inner);
        if let Ok(cn) = num_str.parse::<u32>() {
            return match cn {
                0 => VariantClass::StructuralDeletion,
                2 => VariantClass::Duplication,
                _ => VariantClass::CopyNumberVariation,
            };
        }
    }
    if upper.starts_with("<CNV") || upper.starts_with("<CN") {
        return VariantClass::CopyNumberVariation;
    }
    if upper.starts_with("<CPX") {
        return VariantClass::ComplexStructural;
    }
    if upper.starts_with("<BND") {
        return VariantClass::Translocation;
    }
    if is_paired_breakend_alt(alt) || is_native_single_breakend_alt(alt) {
        return VariantClass::Translocation;
    }
    if let Some(svt) = svtype {
        match svt.to_uppercase().as_str() {
            "DEL" => return VariantClass::StructuralDeletion,
            "INS" => return VariantClass::StructuralInsertion,
            "DUP" => return VariantClass::Duplication,
            "INV" => return VariantClass::Inversion,
            "BND" => return VariantClass::Translocation,
            "CNV" => return VariantClass::CopyNumberVariation,
            _ => {}
        }
    }
    // A `<NON_REF>` with a supported SVTYPE took that class above, as in VEP, where
    // INFO/SVTYPE replaces an ALT outside the VCF 4.4 forms before the SO lookup
    // (Parser/VCF.pm `create_StructuralVariationFeatures`).
    if upper == "<NON_REF>" {
        return VariantClass::ReferenceBlock;
    }
    VariantClass::ComplexStructural
}

/// The structural variant type abbreviation VEP derives for a record: the port of
/// `Parser::get_SO_term` (ensembl-vep) applied to every ALT joined by `/`, with
/// INFO/SVTYPE replacing a joined string that opens outside the VCF 4.4 forms
/// (`Parser/VCF.pm`). Identical parts collapse to one (`<INS>/<INS>` reads `INS`),
/// and a string naming both DEL and DUP, `<DEL>/<DUP:TANDEM>` included, is `CNV`.
pub fn sv_type_abbreviation(joined_alts: &str, svtype: Option<&str>) -> String {
    let svtype = svtype.filter(|v| perl_true(v));
    // Parser/VCF.pm:473, `/^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)/` without `/i`: a
    // lowercase `<del>` is outside the VCF 4.4 forms, so SVTYPE names its type.
    let follows_vcf44 = {
        let t = joined_alts.strip_prefix('<').unwrap_or(joined_alts);
        ["DEL", "INS", "DUP", "INV", "CNV"]
            .iter()
            .any(|p| t.starts_with(p))
            || t.strip_prefix("CN")
                .map(|c| c.strip_prefix('=').unwrap_or(c))
                .is_some_and(|c| c.starts_with(|ch: char| ch.is_ascii_digit()))
    };
    let mut kind = match svtype {
        Some(st) if !follows_vcf44 => st.to_string(),
        _ => joined_alts.to_string(),
    };
    // Repeated parts collapse to the first as written (Parser.pm:668-681): each part is
    // trimmed, loses one leading `<` and one trailing `>` (`s/^\<|\>$//g`) and is
    // compared in upper case, so `<<DEL>>/<DEL>` stays a list; a trailing empty part is
    // no part (`DEL/` is `DEL`).
    if kind.contains(['/', ',']) {
        let parts: Vec<&str> = perl_split(&kind, &['/', ',']).collect();
        let normalized = |part: &str| {
            let part = part.trim();
            let part = part.strip_prefix('<').unwrap_or(part);
            let part = part.strip_suffix('>').unwrap_or(part);
            part.to_ascii_uppercase()
        };
        if let Some(first) = parts.first().map(|p| normalized(p)) {
            if parts.iter().all(|p| normalized(p) == first) {
                kind = parts[0].to_string();
            }
        }
    }
    // The pattern branches of `get_SO_term` (Parser.pm:685-709) match without case
    // and name an upper-case abbreviation; the two fall-through branches (:712-717)
    // keep the type as written, so `del` or `bnd` finds no `%SO_TERMS` key.
    let upper = kind.to_ascii_uppercase();
    let mobile_element = ["INS:ME", "DEL:ME"]
        .into_iter()
        .filter_map(|needle| upper.find(needle))
        .min();
    if let Some(pos) = mobile_element {
        let base = &upper[pos..pos + 3];
        // `(ME):?([A-Z0-9]+)?` (Parser.pm:685): one optional colon, then the element's
        // alphanumerics, so a second colon (`INS:ME::ALU`) leaves the element empty and
        // the type `INS_ME`.
        let rest = &upper[pos + 6..];
        let rest = rest.strip_prefix(':').unwrap_or(rest);
        let element: String = rest
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric())
            .collect();
        let element = if element == "L1" {
            "LINE1".to_string()
        } else {
            element
        };
        let subtype = if ["ALU", "HERV", "LINE1", "SVA"].contains(&element.as_str()) {
            element
        } else {
            "ME".to_string()
        };
        format!("{base}_{subtype}")
    } else if upper.contains("DEL") && upper.contains("DUP") {
        "CNV".to_string()
    } else if upper.contains("DUP:TANDEM") {
        "TDUP".to_string()
    } else if upper.contains("CNV:TR") {
        "TREP".to_string()
    } else if names_a_copy_number(&upper) {
        // Perl's `/^<?CN=?0>?$/` and `/^<?CN=?2>?$/` (Parser.pm:705-706) match the
        // whole string with at most one bracket on each side, in upper case only, so a
        // list of copy numbers, a `<cn0>` or a `<<CN0>>` is `CNV`.
        let bare = kind.strip_prefix('<').unwrap_or(&kind);
        let bare = bare.strip_suffix('>').unwrap_or(bare);
        let copies = bare
            .strip_prefix("CN")
            .map(|c| c.strip_prefix('=').unwrap_or(c));
        match copies {
            Some("0") => "DEL".to_string(),
            Some("2") => "DUP".to_string(),
            _ => "CNV".to_string(),
        }
    } else if upper.contains("CNV") {
        "CNV".to_string()
    } else if kind.contains(['[', ']']) || kind.starts_with('.') || kind.ends_with('.') {
        "BND".to_string()
    } else if kind.starts_with('<') || kind.ends_with('>') {
        // `s/\<|\>//g` then `s/\:.+//g` (Parser.pm:714-715): every bracket goes, and
        // the text from the first colon on when at least one character follows it, so
        // `<DEL:X>` is `DEL` and `<DEL:>` keeps its colon and names no type.
        let mut stripped: String = kind.chars().filter(|c| *c != '<' && *c != '>').collect();
        if let Some(colon) = stripped.find(':') {
            if colon + 1 < stripped.len() {
                stripped.truncate(colon);
            }
        }
        stripped
    } else {
        kind
    }
}

/// `/CN=?[0-9]/i` (Parser.pm:702): `CN`, an optional `=` and a digit, anywhere in
/// the string, so `<CN=X>` names no copy number and `<CNX>/<CN3>` does.
fn names_a_copy_number(upper: &str) -> bool {
    upper.match_indices("CN").any(|(i, _)| {
        let rest = &upper[i + 2..];
        let rest = rest.strip_prefix('=').unwrap_or(rest);
        rest.starts_with(|c: char| c.is_ascii_digit())
    })
}

/// VEP's `%SO_TERMS` (ensembl-variation `Utils/Config.pm`): the SO term of a
/// record whose type abbreviation is one of its keys, with the class that
/// annotates the term, `None` for an unsupported type. Ensembl VEP 116 drops a
/// record of unknown type before annotation (`Parser.pm` `validate_svf`) and
/// writes nothing for it in any format; vep-rs parses and annotates it, carrying
/// its VCF line without consequences and writing an `input`-only JSON object
/// (`docs/intended-divergences.md`).
pub fn sv_so_term_and_class(abbrev: &str) -> Option<(&'static str, VariantClass)> {
    use VariantClass::*;
    Some(match abbrev {
        "INS" => ("insertion", StructuralInsertion),
        "INS_ME" => ("mobile_element_insertion", MobileElementInsertion),
        "INS_ALU" => ("Alu_insertion", MobileElementInsertion),
        "INS_HERV" => ("HERV_insertion", MobileElementInsertion),
        "INS_LINE1" => ("LINE1_insertion", MobileElementInsertion),
        "INS_SVA" => ("SVA_insertion", MobileElementInsertion),
        "DEL" => ("deletion", StructuralDeletion),
        "DEL_ME" => ("mobile_element_deletion", MobileElementDeletion),
        "DEL_ALU" => ("Alu_deletion", MobileElementDeletion),
        "DEL_HERV" => ("HERV_deletion", MobileElementDeletion),
        "DEL_LINE1" => ("LINE1_deletion", MobileElementDeletion),
        "DEL_SVA" => ("SVA_deletion", MobileElementDeletion),
        "TREP" => ("tandem_repeat", TandemRepeat),
        "TDUP" => ("tandem_duplication", TandemDuplication),
        "DUP" => ("duplication", Duplication),
        "CNV" => ("copy_number_variation", CopyNumberVariation),
        "INV" => ("inversion", Inversion),
        "BND" => ("chromosome_breakpoint", Translocation),
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn term(alts: &str, svtype: Option<&str>) -> Option<&'static str> {
        sv_so_term_and_class(&sv_type_abbreviation(alts, svtype)).map(|(term, _)| term)
    }

    /// `/CN=?[0-9]/i` (Parser.pm:702) needs a digit, so `<CN=X>` reaches the
    /// fall-through branch and names the type `CN=X`, which `%SO_TERMS` lacks: Ensembl
    /// VEP 116.2 skips `<CN=X>` as `CN=X is not a supported structural variant type`,
    /// writes `copy_number_variation` for `<CN=5>` and for `<CNX>,<CN3>`, and reads
    /// `SVTYPE=DEL` for `<CN=X>` (outside the VCF 4.4 forms) as `deletion`.
    #[test]
    fn a_copy_number_needs_a_digit() {
        assert_eq!(sv_type_abbreviation("<CN=X>", None), "CN=X");
        assert_eq!(term("<CN=X>", None), None);
        assert_eq!(term("<CN=X>", Some("DEL")), Some("deletion"));
        assert_eq!(term("<CN=5>", None), Some("copy_number_variation"));
        assert_eq!(term("<CNX>/<CN3>", None), Some("copy_number_variation"));
        assert_eq!(term("<CN0>", None), Some("deletion"));
        assert_eq!(term("<CN=2>", None), Some("duplication"));
    }

    /// The repeated-part normalisation strips one `<` and one `>` from each part
    /// (`s/^\<|\>$//g`, Parser.pm:674), so `<<DEL>>/<DEL>` is a list the fall-through
    /// branch reads as `DEL/DEL`, while the fall-through strips every bracket of a single
    /// `<<DEL>>`; `/^<?CN=?0>?$/` allows one bracket a side. Ensembl VEP 116.2 writes
    /// `deletion` for `<<DEL>>`, skips `<<DEL>>,<DEL>` as `DEL/DEL is not a supported
    /// structural variant type`, and writes `copy_number_variation` for `<<CN0>>`.
    #[test]
    fn a_doubled_bracket_is_one_bracket_too_many_for_the_list_and_the_copy_number() {
        assert_eq!(term("<<DEL>>", None), Some("deletion"));
        assert_eq!(sv_type_abbreviation("<<DEL>>/<DEL>", None), "DEL/DEL");
        assert_eq!(term("<<DEL>>/<DEL>", None), None);
        assert_eq!(term("<<CN0>>", None), Some("copy_number_variation"));
        assert_eq!(term("<DEL>/<DEL>", None), Some("deletion"));
        assert_eq!(term("<CN0>/<CN0>", None), Some("deletion"));
    }

    /// `/(INS|DEL):(ME):?([A-Z0-9]+)?/i` (Parser.pm:685) takes one optional colon
    /// before the element, so a second colon leaves the element empty: Ensembl VEP
    /// 116.2 writes `mobile_element_insertion` for `<INS:ME::ALU>` and for `<INS:ME:>`,
    /// and `Alu_insertion` for `<INS:ME:ALU>`.
    #[test]
    fn a_mobile_element_takes_one_colon_before_its_name() {
        assert_eq!(sv_type_abbreviation("<INS:ME::ALU>", Some("INS")), "INS_ME");
        assert_eq!(
            term("<INS:ME::ALU>", Some("INS")),
            Some("mobile_element_insertion")
        );
        assert_eq!(
            term("<INS:ME:>", Some("INS")),
            Some("mobile_element_insertion")
        );
        assert_eq!(term("<INS:ME:ALU>", Some("INS")), Some("Alu_insertion"));
        assert_eq!(term("<INS:MEALU>", None), Some("Alu_insertion"));
        assert_eq!(term("<DEL:ME:L1>", None), Some("LINE1_deletion"));
        assert_eq!(
            term("<INS:ME:ALUY>", None),
            Some("mobile_element_insertion")
        );
    }

    /// `s/\:.+//g` (Parser.pm:715) needs a character after the colon: Ensembl VEP 116.2
    /// writes `deletion` for `<DEL:X>` and skips `<DEL:>` as `DEL: is not a supported
    /// structural variant type`.
    #[test]
    fn a_trailing_colon_stays_in_the_type() {
        assert_eq!(sv_type_abbreviation("<DEL:>", None), "DEL:");
        assert_eq!(term("<DEL:>", None), None);
        assert_eq!(term("<DEL:X>", None), Some("deletion"));
    }

    /// `<NON_REF>` beside another allele is dropped; a list of blocks is kept whole, as
    /// is a list with no block and any list under a true SVTYPE.
    #[test]
    fn a_non_ref_block_beside_an_allele_names_no_variant() {
        let kept = |alt_field: &str, svtype: Option<&str>| {
            let mut alts = written_alts(alt_field);
            variant_alts(&mut alts, svtype);
            alts.iter()
                .map(|(i, alt)| format!("{i}:{alt}"))
                .collect::<Vec<_>>()
        };
        assert_eq!(kept("G,<NON_REF>", None), ["0:G"]);
        assert_eq!(kept("<NON_REF>,G", None), ["1:G"]);
        assert_eq!(kept("<non_ref>,G,T", Some("")), ["1:G", "2:T"]);
        assert_eq!(kept("<NON_REF>", None), ["0:<NON_REF>"]);
        assert_eq!(
            kept("<NON_REF>,<NON_REF>", None),
            ["0:<NON_REF>", "1:<NON_REF>"]
        );
        assert_eq!(kept("G,T", None), ["0:G", "1:T"]);
        assert_eq!(kept("G,<NON_REF>", Some("DEL")), ["0:G", "1:<NON_REF>"]);
        assert_eq!(term("<NON_REF>/<NON_REF>", None), None);
        assert_eq!(sv_type_abbreviation("<NON_REF>/<NON_REF>", None), "NON_REF");
    }

    /// Perl's list-context `split` drops trailing empty fields and keeps middle ones
    /// (`split(",", "<DEL>,")` is one field, `"G,,T,,"` three): Ensembl VEP 116.2
    /// annotates `<DEL>,` with `SVTYPE=DEL` as a deletion, writes the one allele `G` for
    /// `A G,`, and reads `SVTYPE=DEL/` as `DEL` (Parser.pm:669).
    #[test]
    fn a_trailing_empty_field_is_no_field() {
        assert_eq!(written_alts("<DEL>,"), [(0, "<DEL>")]);
        assert_eq!(written_alts("G,,T,,"), [(0, "G"), (1, ""), (2, "T")]);
        assert_eq!(written_alts(","), []);
        assert_eq!(written_alts(""), []);
        assert_eq!(perl_split("DEL/", &['/', ',']).collect::<Vec<_>>(), ["DEL"]);
        assert_eq!(sv_type_abbreviation("DEL/", None), "DEL");
        assert_eq!(term("<DEL>/", Some("DEL")), Some("deletion"));
        assert_eq!(sv_type_abbreviation("/", None), "/");
    }

    /// Perl's string truth: the empty string and `0` are false, every other string is
    /// true, `0.0`, `00` and `false` included.
    #[test]
    fn perl_truth_is_string_truth() {
        assert!(!perl_true(""));
        assert!(!perl_true("0"));
        assert!(perl_true("1"));
        assert!(perl_true("0.0"));
        assert!(perl_true("00"));
        assert!(perl_true("false"));
        assert!(perl_true(" 0"));
    }

    /// Parser/VCF.pm:237-245 over REF, the ALTs and a true SVTYPE.
    #[test]
    fn routing_follows_the_reference_and_alleles_and_a_true_svtype() {
        let routes = |reference: &str, alt_field: &str, svtype: Option<&str>| {
            routes_to_structural(reference, &written_alts(alt_field), svtype)
        };
        assert!(routes("N", "G", Some("DEL")));
        assert!(!routes("A", "G", Some("DEL")));
        assert!(!routes("N", "G", None));
        assert!(!routes("N", "G", Some("")));
        assert!(!routes("N", "G", Some("0")));
        assert!(routes("N", "<*>", Some("DEL")));
        assert!(!routes("N", "<*>", None));
        assert!(routes("A", "<DEL>", None));
        assert!(routes("A", "G,<NON_REF>", None));
        assert!(routes("A", ".N,G", None));
        assert!(!routes("A", "G,.N,T", None));
        assert!(routes("A", "g", Some("DEL")));
        // The dot forms bind to the joined list's ends (`^\.\w+|\w+\.$` over
        // `join(",", @$alts)`), a token never crosses a `*`.
        assert!(routes("A", "G,N.", None));
        assert!(!routes("A", "<*>,G", None));
        assert!(routes("A", "G,<DEL>", None) && routes("A", "N[21:100[,G", None));
    }

    /// One StructuralVariationFeature carries a supported record (Parser/VCF.pm:560-574),
    /// so every allele takes its class; an unsupported list keeps each symbolic ALT's own
    /// class, INFO/SVTYPE naming it for a mixed record (`:472-475`).
    #[test]
    fn every_allele_of_a_supported_record_takes_its_class() {
        let alts = written_alts("G,<DEL>");
        assert_eq!(
            per_alt_classes(&alts, Some(VariantClass::StructuralDeletion), Some("DEL")),
            [
                (0, "G", VariantClass::StructuralDeletion),
                (1, "<DEL>", VariantClass::StructuralDeletion)
            ]
        );
        assert_eq!(
            per_alt_classes(&written_alts("<INS>,<DEL>"), None, None),
            [
                (0, "<INS>", VariantClass::StructuralInsertion),
                (1, "<DEL>", VariantClass::StructuralDeletion)
            ]
        );
        assert_eq!(
            per_alt_classes(&written_alts("G,<INV>"), None, Some("FOO")),
            [(1, "<INV>", VariantClass::ComplexStructural)]
        );
        assert_eq!(per_alt_classes(&written_alts("G"), None, Some("FOO")), []);
    }
}
