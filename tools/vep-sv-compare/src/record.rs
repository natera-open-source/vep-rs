// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What both engines make of an input record: the route Ensembl VEP's parser takes it, the
//! Location and the alleles each prints.
//!
//! The reference's rules are those of release 116: the type string of Parser/VCF.pm:468-481,
//! the abbreviations of Parser.pm 116.2:663-728 (`get_SO_term`), the routing of a record to
//! the structural path at Parser/VCF.pm:237-245 and the coordinates of ensembl-io BaseVCF4.pm
//! `get_start` and `get_end`: the base before the event is the start, SVLEN sets the end
//! before END does, and a breakend's END is discarded. A coordinate is held in 128 bits and
//! a sum past that word stays at its bound.

use std::collections::{BTreeSet, HashMap};
use std::sync::LazyLock;

use regex::Regex;

use crate::locus::{location_span, strip_chr};
use crate::normalize::{normalize_allele, normalize_location};
use crate::perl::{
    info_fields, info_value, info_written, is_perl_infinity, is_perl_nan, is_perl_number,
    perl_fields, perl_numeric_value, perl_true, PerlNumber,
};
use crate::pyfmt::py_strip;
use crate::rules::gvcf::{gvcf_site_rows, perl_split, sequence_location};
use crate::strings::{Id, Strings, Tuple, TupleSet};
use crate::vcf::VcfVariant;

/// The Sequence Ontology term of each structural type abbreviation (ensembl-variation
/// Utils/Config.pm:1417-1438 at release 116).
pub const SO_TERMS: [(&str, &str); 18] = [
    ("INS", "insertion"),
    ("INS_ME", "mobile_element_insertion"),
    ("INS_ALU", "Alu_insertion"),
    ("INS_HERV", "HERV_insertion"),
    ("INS_LINE1", "LINE1_insertion"),
    ("INS_SVA", "SVA_insertion"),
    ("DEL", "deletion"),
    ("DEL_ME", "mobile_element_deletion"),
    ("DEL_ALU", "Alu_deletion"),
    ("DEL_HERV", "HERV_deletion"),
    ("DEL_LINE1", "LINE1_deletion"),
    ("DEL_SVA", "SVA_deletion"),
    ("TREP", "tandem_repeat"),
    ("TDUP", "tandem_duplication"),
    ("DUP", "duplication"),
    ("CNV", "copy_number_variation"),
    ("INV", "inversion"),
    ("BND", "chromosome_breakpoint"),
];

/// The Alleles either engine writes for a breakend record's own end.
const OWN_END_ALLELES: [&str; 3] = ["N.", ".N", "chromosome_breakpoint"];

/// The class vep-rs gives a symbolic ALT by its own prefix, read without regard to case
/// (`classify_symbolic_alt` in crates/vep-core/src/sv_routing.rs), first match first; a
/// `<CN...>` with a copy number the table knows (`<CN0>`, `<CN2>`) is a supported type and
/// never reaches this list.
const VEP_RS_ALT_PREFIX_TERMS: [(&str, &str); 11] = [
    ("<INS:ME", "mobile_element_insertion"),
    ("<DEL:ME", "mobile_element_deletion"),
    ("<DEL", "deletion"),
    ("<INS", "insertion"),
    ("<DUP:TANDEM", "tandem_duplication"),
    ("<DUP", "duplication"),
    ("<INV", "inversion"),
    ("<CNV:TR", "tandem_repeat"),
    ("<CNV", "copy_number_variation"),
    ("<CN", "copy_number_variation"),
    ("<BND", "chromosome_breakpoint"),
];

const VEP_RS_SVTYPE_TERMS: [(&str, &str); 6] = [
    ("DEL", "deletion"),
    ("INS", "insertion"),
    ("DUP", "duplication"),
    ("INV", "inversion"),
    ("BND", "chromosome_breakpoint"),
    ("CNV", "copy_number_variation"),
];

fn regex(pattern: &str) -> Regex {
    Regex::new(pattern).expect("regex")
}

/// The ALT forms VCF 4.4 names a type by (Parser/VCF.pm:471).
static VCF44_SYMBOLIC: LazyLock<Regex> =
    LazyLock::new(|| regex(r"^<?(?:DEL|INS|DUP|INV|CNV|CN=?[0-9]+)"));

/// Parser/VCF.pm:241, the ALT forms that route a line to the structural path. A word
/// character is a letter, a number of any kind or `_`; the dot at the end may stand before
/// a final line feed.
static STRUCTURAL_ALT_LIST: LazyLock<Regex> =
    LazyLock::new(|| regex(r"[<\[\]][^*]+[>\]\[]|^\.[\p{L}\p{N}_]+|[\p{L}\p{N}_]+\.\n?$"));

/// A bracket-notation mate, its position in the decimal digits of any script.
static BRACKET_MATE: LazyLock<Regex> = LazyLock::new(|| regex(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]"));

/// A dotted capital and a dotless small `i` stand for `I` beside the letters the pattern
/// names, as the case-insensitive match this reading is defined by takes them.
static MOBILE_ELEMENT: LazyLock<Regex> =
    LazyLock::new(|| regex(r"(?i)([I\x{130}\x{131}]NS|DEL):(ME):?([A-Z0-9\x{130}\x{131}]+)?"));
static HAS_DEL: LazyLock<Regex> = LazyLock::new(|| regex(r"(?i)DEL"));
static HAS_DUP: LazyLock<Regex> = LazyLock::new(|| regex(r"(?i)DUP"));
static HAS_TANDEM_DUP: LazyLock<Regex> = LazyLock::new(|| regex(r"(?i)DUP:TANDEM"));
static HAS_TANDEM_REPEAT: LazyLock<Regex> = LazyLock::new(|| regex(r"(?i)CNV:TR"));
static HAS_COPY_NUMBER: LazyLock<Regex> = LazyLock::new(|| regex(r"(?i)CN=?[0-9]"));
static COPY_NUMBER_ZERO: LazyLock<Regex> = LazyLock::new(|| regex(r"^<?CN=?0>?\n?$"));
static COPY_NUMBER_TWO: LazyLock<Regex> = LazyLock::new(|| regex(r"^<?CN=?2>?\n?$"));
static HAS_CNV: LazyLock<Regex> = LazyLock::new(|| regex(r"(?i)CNV"));
static BREAKEND_FORM: LazyLock<Regex> = LazyLock::new(|| regex(r"[\[\]]|^\.|\.\n?$"));
static BRACKETED: LazyLock<Regex> = LazyLock::new(|| regex(r"^<|>\n?$"));
static FROM_COLON: LazyLock<Regex> = LazyLock::new(|| regex(r":.+"));
static DOT_BREAKEND: LazyLock<Regex> = LazyLock::new(|| regex(r"^\.?[A-Za-z]+\.?$"));

/// The term of an abbreviation of [`SO_TERMS`].
pub fn so_term(abbrev: &str) -> Option<&'static str> {
    SO_TERMS
        .iter()
        .find(|(held, _)| *held == abbrev)
        .map(|(_, term)| *term)
}

/// The chromosome of a record as its Location names it.
pub fn record_chromosome(v: &VcfVariant) -> &str {
    v.location
        .split_once(':')
        .map_or(v.location.as_str(), |(chrom, _)| chrom)
}

/// The type string Ensembl VEP keeps for a record (Parser/VCF.pm:469-475): INFO/SVTYPE when
/// Perl reads it as true and the ALT string does not open with a VCF 4.4 form, else the ALT
/// string (`/`-joined when a list).
pub fn perl_record_type<'a>(sv_type: &'a str, svtype: Option<&'a str>) -> &'a str {
    match svtype {
        Some(named) if perl_true(svtype) && !VCF44_SYMBOLIC.is_match(sv_type) => named,
        _ => sv_type,
    }
}

/// The Sequence Ontology abbreviation `get_SO_term` derives for a type string, or `None`
/// when the reference has no term for it (the type is then `vep_skip`).
pub fn perl_so_abbrev(sv_type: &str) -> Option<&'static str> {
    let mut sv_type = sv_type;
    if sv_type.contains(['/', ',']) {
        let parts = perl_fields(sv_type.split(['/', ',']).collect());
        let normalized = |part: &str| {
            let part = py_strip(part);
            let part = part.strip_prefix('<').unwrap_or(part);
            part.strip_suffix('>').unwrap_or(part).to_uppercase()
        };
        if let Some(first) = parts.first() {
            let named = normalized(first);
            if parts.iter().all(|part| normalized(part) == named) {
                sv_type = first;
            }
        }
    }
    let abbrev = if let Some(me) = MOBILE_ELEMENT.captures(sv_type) {
        let mut subtype = me[2].to_uppercase();
        if let Some(element) = me.get(3) {
            let element = match element.as_str().to_uppercase() {
                short if short == "L1" => String::from("LINE1"),
                other => other,
            };
            if matches!(element.as_str(), "ALU" | "HERV" | "LINE1" | "SVA") {
                subtype = element;
            }
        }
        format!("{}_{subtype}", me[1].to_uppercase())
    } else if HAS_DEL.is_match(sv_type) && HAS_DUP.is_match(sv_type) {
        String::from("CNV")
    } else if HAS_TANDEM_DUP.is_match(sv_type) {
        String::from("TDUP")
    } else if HAS_TANDEM_REPEAT.is_match(sv_type) {
        String::from("TREP")
    } else if HAS_COPY_NUMBER.is_match(sv_type) {
        if COPY_NUMBER_TWO.is_match(sv_type) {
            String::from("DUP")
        } else if COPY_NUMBER_ZERO.is_match(sv_type) {
            String::from("DEL")
        } else {
            String::from("CNV")
        }
    } else if HAS_CNV.is_match(sv_type) {
        String::from("CNV")
    } else if BREAKEND_FORM.is_match(sv_type) {
        String::from("BND")
    } else if BRACKETED.is_match(sv_type) {
        FROM_COLON
            .replace_all(&sv_type.replace(['<', '>'], ""), "")
            .into_owned()
    } else {
        sv_type.to_string()
    };
    SO_TERMS
        .iter()
        .find(|(held, _)| *held == abbrev)
        .map(|(held, _)| *held)
}

/// The end `get_end` derives for a symbolic record at `pos` (BaseVCF4.pm:297-317) with the
/// breakpoint rule of Parser/VCF.pm:487-489: `pos + |SVLEN|` from the first value as Perl
/// numifies it (`SVLEN=` and `SVLEN=.` are 0, the pair before the start); else END as a
/// number, `None` when it is not one (`validate_vf` skips the record); else the start. A
/// breakpoint reads END only when Perl reads it as false (`0` or empty), which keeps the
/// record on the skip route.
pub fn perl_structural_end(
    pos: i128,
    svlen: Option<&str>,
    end: Option<&str>,
    breakpoint: bool,
) -> Option<PerlNumber> {
    let start = PerlNumber::Integer(pos.saturating_add(1));
    if let Some(svlen) = svlen {
        let first = svlen.split(',').next().unwrap_or("");
        return Some(perl_numeric_value(first).abs() + pos);
    }
    if breakpoint {
        return match end {
            None => Some(start),
            Some(end) if perl_true(Some(end)) => Some(start),
            Some("0") => Some(PerlNumber::Integer(0)),
            Some(_) => None,
        };
    }
    let Some(end) = end else {
        return Some(start);
    };
    is_perl_number(py_strip(end)).then(|| perl_numeric_value(end))
}

/// Whether Parser/VCF.pm:237-245 routes the record to the sequence path: REF and every ALT
/// together a run of `ACGT`, or no true INFO/SVTYPE and no ALT written in a structural
/// form. A sequence record's Location starts at POS; a structural record's at POS + 1.
pub fn perl_sequence_path(v: &VcfVariant) -> bool {
    let alts = perl_split(&v.alt, ',');
    let svtype = info_written(&v.info, "SVTYPE");
    let mut letters = v
        .reference
        .bytes()
        .chain(alts.iter().flat_map(|alt| alt.bytes()))
        .peekable();
    let sequence = letters.peek().is_some()
        && letters.all(|letter| matches!(letter, b'A' | b'C' | b'G' | b'T'));
    sequence || !(perl_true(svtype) || STRUCTURAL_ALT_LIST.is_match(&alts.join(",")))
}

/// Whether an ALT is the `<NON_REF>` block of a gVCF, in any case.
fn is_non_ref_block(alt: &str) -> bool {
    alt.to_uppercase() == "<NON_REF>"
}

/// Whether an ALT list is a gVCF variant site vep-rs reads on its sequence path: `<NON_REF>`
/// beside alleles that name a variant (not the block, `*` or `.`), no true INFO/SVTYPE, and
/// no remaining allele written in a structural form (Parser/VCF.pm:237-245). The structural
/// forms are those of [`perl_sequence_path`], word characters included.
fn is_gvcf_variant_site(alts: &[&str], svtype: Option<&str>) -> bool {
    let real: Vec<&str> = alts
        .iter()
        .copied()
        .filter(|alt| !is_non_ref_block(alt) && *alt != "*" && *alt != ".")
        .collect();
    if real.is_empty() || !alts.iter().any(|alt| is_non_ref_block(alt)) || perl_true(svtype) {
        return false;
    }
    !STRUCTURAL_ALT_LIST.is_match(&real.join(","))
}

/// The `(chromosome, start)` the reference prints for a record: POS on the sequence path,
/// POS + 1 on the structural path (BaseVCF4.pm:263).
pub fn reference_start(v: &VcfVariant) -> (String, i128) {
    let start = if perl_sequence_path(v) {
        v.pos
    } else {
        v.pos.saturating_add(1)
    };
    (record_chromosome(v).to_string(), start)
}

/// The route by which Ensembl VEP 116.2 drops a VCF record before its input buffer, or
/// `None` when it keeps the record: `unsupported_type` (Parser/VCF.pm:477-481), or
/// `gvcf_variant_site` when the unsupported type lists `<NON_REF>` beside another ALT,
/// `incomplete_deletion` (:555-558), `invalid_span` (Parser.pm:475-488: a start or end that
/// is not a number, or a start past the end by more than one), `oversize` (Parser.pm:491-501
/// over the span `get_end` derives: SVLEN first, else END, else the start; a breakend's END
/// is discarded). An SVLEN or END spelled `Inf` is `oversize` under a cap (Perl's
/// `abs("Inf")`); one spelled `NaN` keeps the record, every comparison on NaN being false.
/// `max_sv_size` is the cap the reference ran with; -1 is no cap, as for the flag. A record
/// the parser routes to the sequence path is kept.
pub fn perl_skip_route(v: &VcfVariant, max_sv_size: i128) -> Option<&'static str> {
    if perl_sequence_path(v) {
        return None;
    }
    let alts = perl_split(&v.alt, ',');
    let svtype = info_written(&v.info, "SVTYPE");
    let Some(abbrev) = perl_so_abbrev(perl_record_type(&alts.join("/"), svtype)) else {
        return Some(if is_gvcf_variant_site(&alts, svtype) {
            "gvcf_variant_site"
        } else {
            "unsupported_type"
        });
    };
    let start = v.pos.saturating_add(1);
    let svlen = info_written(&v.info, "SVLEN");
    let end_field = info_written(&v.info, "END");
    let is_breakend = abbrev == "BND";
    let end = perl_structural_end(v.pos, svlen, end_field, is_breakend);
    // `abs("NaN")` and `abs("Inf")` are the values the spellings name. Every comparison on
    // NaN is false (Parser/VCF.pm:555, Parser.pm:483, 491-494), so the record is kept and
    // printed with the text; Inf is a span over any cap, and under no cap the reference run
    // dies on it.
    let span_text = match svlen {
        Some(svlen) => svlen.split(',').next().unwrap_or(""),
        None if is_breakend => "",
        None => end_field.unwrap_or(""),
    };
    if is_perl_nan(span_text) {
        return None;
    }
    if is_perl_infinity(span_text) {
        return (max_sv_size != -1).then_some("oversize");
    }
    // Parser/VCF.pm:555-558 compares `$start > $end` before `validate_vf` runs, numifying an
    // END that is not a number to 0, so a deletion with `END=abc` is incomplete first.
    let compared_end = end.unwrap_or_else(|| perl_numeric_value(end_field.unwrap_or("")));
    let deletion = so_term(abbrev).is_some_and(|term| term.contains("del"));
    if deletion && (compared_end < start || (svlen.is_none() && end_field.is_none())) {
        return Some("incomplete_deletion");
    }
    let Some(end) = end else {
        return Some("invalid_span");
    };
    if end + 1 < start {
        return Some("invalid_span");
    }
    if max_sv_size != -1 && end - start > max_sv_size {
        return Some("oversize");
    }
    None
}

/// The Allele vep-rs writes for a symbolic ALT whose type Ensembl VEP has no term for
/// (`classify_symbolic_alt` with `display_allele`): the class term of the ALT's own prefix
/// read without regard to case (`<del>`, `<DEL:>` and `<DEL >` are `deletion`, `<CN=X>` and
/// `<CN>` `copy_number_variation`, `<Bnd>` `chromosome_breakpoint`), else the class
/// INFO/SVTYPE names, else the type string Ensembl VEP keeps for the record
/// (Parser/VCF.pm:469-481): INFO/SVTYPE as written when it is set and the ALT does not open
/// with a VCF 4.4 form, else the ALT as written (`<<DEL>>`, `<:DEL>`, `<CPX:X>`).
pub fn vep_rs_unsupported_alt_allele<'a>(alt: &'a str, svtype: Option<&'a str>) -> &'a str {
    let upper = alt.to_uppercase();
    if let Some((_, term)) = VEP_RS_ALT_PREFIX_TERMS
        .iter()
        .find(|(prefix, _)| upper.starts_with(prefix))
    {
        return term;
    }
    if let Some(named) = svtype {
        let named = named.to_uppercase();
        if let Some((_, term)) = VEP_RS_SVTYPE_TERMS.iter().find(|(held, _)| *held == named) {
            return term;
        }
    }
    perl_record_type(alt, svtype)
}

/// Every Allele string either engine writes for a VCF record, normalised as the comparison
/// normalises them: the Sequence Ontology class of each symbolic ALT and of the whole ALT
/// list (`copy_number_variation` for a DEL beside a DUP), the raw abbreviation of a type
/// with no term and its bracketed spelling (`CPX`, `<CPX>`, `<NON_REF>`) and the string
/// vep-rs prints for it ([`vep_rs_unsupported_alt_allele`]), each bracket breakend string
/// with `chr` stripped inside the brackets, the own-end strings of a breakend record (`N.`,
/// `.N`, `chromosome_breakpoint`, the REF base with the dot), the mate string a `<BND>`
/// record's CHR2 with END2, POS2 or END gives, and a sequence ALT as written.
pub fn record_alleles(v: &VcfVariant) -> BTreeSet<String> {
    let mut alleles: BTreeSet<String> = BTreeSet::new();
    let info = info_fields(&v.info);
    let alts = perl_split(&v.alt, ',');
    let svtype = info_value(&info, "SVTYPE");
    let joined = alts.join("/");
    let record_type = perl_record_type(&joined, svtype);
    let whole = perl_so_abbrev(record_type);
    let own_end = |alleles: &mut BTreeSet<String>| {
        alleles.extend(OWN_END_ALLELES.map(String::from));
        alleles.insert(format!("{}.", v.reference));
        alleles.insert(format!(".{}", v.reference));
    };
    // Under a type the table lacks, vep-rs keeps each symbolic ALT's own class, INFO/SVTYPE
    // naming the class of a mixed list (a symbolic ALT beside a sequence allele).
    let mixed = alts.iter().any(|alt| alt.starts_with('<'))
        && alts.iter().any(|alt| {
            !alt.starts_with('<') && !BRACKET_MATE.is_match(alt) && *alt != "*" && *alt != "."
        });
    for alt in &alts {
        if alt.starts_with('<') && alt.ends_with('>') {
            match perl_so_abbrev(alt).and_then(so_term) {
                Some(term) => {
                    alleles.insert(term.to_string());
                }
                None => {
                    let raw = FROM_COLON.replace_all(&alt[1..alt.len() - 1], "");
                    alleles.insert(format!("<{raw}>"));
                    alleles.insert(raw.into_owned());
                    alleles.insert(alt.to_string());
                }
            }
            if whole.is_none() {
                let bracketed;
                let own = match svtype {
                    Some(named) if mixed && perl_true(svtype) => {
                        bracketed = format!("<{named}>");
                        bracketed.as_str()
                    }
                    _ => alt,
                };
                let written = vep_rs_unsupported_alt_allele(own, svtype);
                alleles.insert(written.to_string());
                if written == own || Some(written) == svtype {
                    // An ALT with no class term prints the record's type string.
                    alleles.insert(record_type.to_string());
                }
            }
        } else if BRACKET_MATE.is_match(alt) {
            alleles.insert(normalize_allele(alt, "vep-rs"));
            own_end(&mut alleles);
        } else if DOT_BREAKEND.is_match(alt) && (alt.starts_with('.') || alt.ends_with('.')) {
            alleles.insert(alt.to_string());
            alleles.insert(String::from("chromosome_breakpoint"));
        } else {
            alleles.insert(alt.to_string());
        }
    }
    if let Some(term) = whole.and_then(so_term) {
        alleles.insert(term.to_string());
    }
    if alleles.contains("") {
        // The tab format writes an empty Allele as `-`.
        alleles.insert(String::from("-"));
    }
    if whole == Some("BND") || svtype == Some("BND") {
        own_end(&mut alleles);
        if let Some(chr2) = info_value(&info, "CHR2").filter(|chr2| !chr2.is_empty()) {
            // The mate string the parser builds from INFO (Parser/VCF.pm:488-497): the mate
            // position is END2 when present, else the record's END, which a breakend record
            // from a caller that writes no END2 carries as its mate position.
            let chr2 = strip_chr(chr2);
            for key in ["END2", "POS2", "END"] {
                if let Some(pos2) = info_value(&info, key).filter(|pos2| !pos2.is_empty()) {
                    alleles.insert(format!("N[{chr2}:{pos2}["));
                }
            }
        }
    }
    alleles
}

/// The `(chromosome, start, end)` both engines write for a structural record: the base
/// after POS, and the end `get_end` derives (SVLEN first, else END, else the start; a
/// breakend's END discarded). A breakend is a record whose type string
/// ([`perl_record_type`]) has the term `chromosome_breakpoint`: a `<BND>` with no SVTYPE is
/// one, a bracket ALT under `SVTYPE=DEL` is not.
pub fn record_location(v: &VcfVariant) -> (String, i128, i128) {
    let start = v.pos.saturating_add(1);
    let svlen = info_written(&v.info, "SVLEN");
    let end_field = info_written(&v.info, "END");
    let joined = perl_split(&v.alt, ',').join("/");
    let record_type = perl_record_type(&joined, info_written(&v.info, "SVTYPE"));
    let is_breakend = perl_so_abbrev(record_type) == Some("BND");
    let end = perl_structural_end(v.pos, svlen, end_field, is_breakend);
    (
        record_chromosome(v).to_string(),
        start,
        end.map_or(start, PerlNumber::int),
    )
}

/// The normalised Location of a span: `chrom:start` for a point, else `chrom:start-end`.
pub fn span_location(chrom: &str, start: i128, end: i128) -> String {
    normalize_location(&if end == start {
        format!("{chrom}:{start}")
    } else {
        format!("{chrom}:{start}-{end}")
    })
}

/// The normalised Location both engines print for a span, its low end first: `chrom:start`
/// for a point, else `chrom:lo-hi`.
pub fn printed_location(chrom: &str, start: i128, end: i128) -> String {
    normalize_location(&sequence_location(chrom, start, end))
}

/// Whether a reference row's span lies on `chrom` within one base of `start`..`end`: the
/// reference wrote the record (a `<CNV:TR>` expanded to a literal allele sits at the run's
/// 3' end, one base past it).
pub fn reference_row_beside(
    spans: &[(&str, i128, i128)],
    chrom: &str,
    start: i128,
    end: i128,
) -> bool {
    spans.iter().any(|span| {
        span.0 == chrom && span.1 <= end.saturating_add(1) && span.2 >= start.saturating_sub(1)
    })
}

/// The `(Location, Allele)` pairs vep-rs prints for a record: on the structural path the
/// span [`record_location`] derives, printed low end first (`21:30425650-30425651` for a
/// deletion whose end is before its start), under every allele of [`record_alleles`]; on
/// the sequence path each real allele at the place the parser gives it.
pub fn vep_rs_record_rows(v: &VcfVariant) -> BTreeSet<(String, String)> {
    let chrom = record_chromosome(v);
    if perl_sequence_path(v) {
        return gvcf_site_rows(v.pos, &v.reference, &perl_split(&v.alt, ','))
            .into_iter()
            .map(|(allele, start, end)| (printed_location(chrom, start, end), allele))
            .collect();
    }
    let (_, start, end) = record_location(v);
    let location = printed_location(chrom, start, end);
    record_alleles(v)
        .into_iter()
        .map(|allele| (location.clone(), allele))
        .collect()
}

/// The Locations of a tuple set that name a span, by the chromosome and the start they
/// name, each with its end.
pub struct LocatedSet<'a> {
    set: &'a TupleSet,
    by_chrom: HashMap<&'a str, EndsByStart>,
}

/// The Locations of one chromosome by their start, each with its end.
type EndsByStart = HashMap<i128, Vec<(Id, i128)>>;

impl<'a> LocatedSet<'a> {
    pub fn new(set: &'a TupleSet, strings: &'a Strings) -> LocatedSet<'a> {
        let mut by_chrom: HashMap<&'a str, EndsByStart> = HashMap::new();
        for run in set.as_slice().chunk_by(|a, b| a.location() == b.location()) {
            let location = run[0].location();
            if let Some((chrom, start, end)) = location_span(strings.text(location)) {
                by_chrom
                    .entry(chrom)
                    .or_default()
                    .entry(start)
                    .or_default()
                    .push((location, end));
            }
        }
        LocatedSet { set, by_chrom }
    }

    /// The Locations starting at `start` on `chrom`, each with its end.
    pub fn locations_at(&self, chrom: &str, start: i128) -> &[(Id, i128)] {
        self.by_chrom
            .get(chrom)
            .and_then(|starts| starts.get(&start))
            .map_or(&[], Vec::as_slice)
    }

    /// Whether a Location of the set starts at `start` on `chrom`.
    pub fn starts_at(&self, chrom: &str, start: i128) -> bool {
        !self.locations_at(chrom, start).is_empty()
    }

    /// The tuples of the set at one Location.
    pub fn tuples_at(&self, location: Id) -> &'a [Tuple] {
        self.set.with_prefix(&[location])
    }

    /// The tuples of the set at one Location under one Allele.
    pub fn tuples_under(&self, location: Id, allele: Id) -> &'a [Tuple] {
        self.set.with_prefix(&[location, allele])
    }

    /// Every span of the set as `(chromosome, start, end)`, in no fixed order.
    pub fn spans(&self) -> Vec<(&'a str, i128, i128)> {
        let mut spans = Vec::new();
        for (chrom, starts) in &self.by_chrom {
            for (start, ends) in starts {
                spans.extend(ends.iter().map(|(_, end)| (*chrom, *start, *end)));
            }
        }
        spans
    }
}

/// The tuples of a set at the Location of this text; none when no tuple of the run carries
/// the text.
pub fn tuples_at<'a>(set: &'a TupleSet, strings: &Strings, location: &str) -> &'a [Tuple] {
    match strings.get(location) {
        Some(location) => set.with_prefix(&[location]),
        None => &[],
    }
}

/// The tuples of a set at the Location and under the Allele of these texts; none when no
/// tuple of the run carries one of them.
pub fn tuples_under<'a>(
    set: &'a TupleSet,
    strings: &Strings,
    location: &str,
    allele: &str,
) -> &'a [Tuple] {
    match (strings.get(location), strings.get(allele)) {
        (Some(location), Some(allele)) => set.with_prefix(&[location, allele]),
        _ => &[],
    }
}

/// The reference's tuples at the start it prints for a record (POS + 1 on the structural
/// path, POS on the sequence path) on the record's chromosome, at Locations the scored side
/// wrote nothing for; the end is not matched, since the reference's arithmetic on a
/// numified POS (a negative start) is not reproduced.
pub fn reference_rows_starting_at(
    v: &VcfVariant,
    perl: &LocatedSet<'_>,
    rust: &TupleSet,
) -> Vec<Tuple> {
    let (chrom, start) = reference_start(v);
    perl.locations_at(&chrom, start)
        .iter()
        .filter(|(location, _)| rust.with_prefix(&[*location]).is_empty())
        .flat_map(|(location, _)| perl.tuples_at(*location).iter().copied())
        .collect()
}

#[cfg(test)]
pub(crate) mod testing {
    use super::*;
    use crate::strings::OwnedSides;

    /// An input record as the VCF reader builds it.
    pub(crate) fn rec_on(
        chrom: &str,
        pos: i128,
        reference: &str,
        alt: &str,
        info: &str,
        vid: &str,
    ) -> VcfVariant {
        VcfVariant {
            chrom: chrom.to_string(),
            pos,
            vid: vid.to_string(),
            reference: reference.to_string(),
            alt: alt.to_string(),
            variant_type: crate::vcf::classify_variant(reference, alt).to_string(),
            location: format!("{chrom}:{pos}"),
            info: info.to_string(),
            pos_text: None,
        }
    }

    /// A record of chromosome 21 named `vid`.
    pub(crate) fn rec_id(
        pos: i128,
        reference: &str,
        alt: &str,
        info: &str,
        vid: &str,
    ) -> VcfVariant {
        rec_on("21", pos, reference, alt, info, vid)
    }

    /// A record of chromosome 21 named `r`.
    pub(crate) fn rec(pos: i128, reference: &str, alt: &str, info: &str) -> VcfVariant {
        rec_id(pos, reference, alt, info, "r")
    }

    /// One annotation tuple as an output reader builds it: the Location and the Allele
    /// normalised, the consequence terms sorted and each once.
    pub(crate) fn tuple_on(
        loc: &str,
        allele: &str,
        feature: &str,
        feature_type: &str,
        csq: &str,
    ) -> [String; 5] {
        let terms: BTreeSet<&str> = csq
            .split(',')
            .map(py_strip)
            .filter(|term| !term.is_empty())
            .collect();
        [
            normalize_location(loc),
            normalize_allele(allele, "vep-rs"),
            feature.to_string(),
            feature_type.to_string(),
            terms.into_iter().collect::<Vec<_>>().join(","),
        ]
    }

    /// A tuple on a transcript.
    pub(crate) fn t(loc: &str, allele: &str, feature: &str, csq: &str) -> [String; 5] {
        tuple_on(loc, allele, feature, "Transcript", csq)
    }

    /// The row a record with no overlapping feature gets.
    pub(crate) fn intergenic(loc: &str, allele: &str) -> [String; 5] {
        tuple_on(loc, allele, "-", "-", "intergenic_variant")
    }

    /// The records of a VCF written from its eight-column data lines, as the VCF reader
    /// reads them.
    pub(crate) fn vcf(records: &[[&str; 8]]) -> Vec<VcfVariant> {
        let dir = crate::TempDir::new("vep_sv_compare_records_").expect("a scratch directory");
        let path = dir.path.join("in.vcf");
        let mut text =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for record in records {
            text.push_str(&record.join("\t"));
            text.push('\n');
        }
        std::fs::write(&path, text).expect("the VCF is written");
        crate::vcf::parse_vcf(&path).expect("the VCF is read")
    }

    /// A Location as the comparison holds it.
    pub(crate) fn normalised(loc: &str) -> String {
        normalize_location(loc)
    }

    /// The two sides of a comparison from the tuples of each.
    pub(crate) fn owned(perl: &[[String; 5]], rust: &[[String; 5]]) -> OwnedSides {
        OwnedSides::from_rows(&rows(perl), &rows(rust))
    }

    /// A set of tuples the sides hold, for an expected value.
    pub(crate) fn set_of(sides: &OwnedSides, tuples: &[[String; 5]]) -> TupleSet {
        sides.set(&rows(tuples))
    }

    /// Rows of owned texts as the rows of borrowed texts the tuple sets are built from.
    pub(crate) fn rows(owned: &[[String; 5]]) -> Vec<[&str; 5]> {
        owned
            .iter()
            .map(|row| [0, 1, 2, 3, 4].map(|column| row[column].as_str()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::testing::{intergenic, rec, rows, t};
    use super::*;
    use crate::strings::OwnedSides;

    const SATURATED: i128 = (1 << 63) - 1;

    fn alleles(v: &VcfVariant) -> Vec<String> {
        record_alleles(v).into_iter().collect()
    }

    fn sorted(texts: &[&str]) -> Vec<String> {
        let set: BTreeSet<String> = texts.iter().map(|text| text.to_string()).collect();
        set.into_iter().collect()
    }

    fn location(v: &VcfVariant) -> (String, i128, i128) {
        record_location(v)
    }

    fn at(start: i128, end: i128) -> (String, i128, i128) {
        (String::from("21"), start, end)
    }

    fn route(pos: i128, alt: &str, info: &str) -> Option<&'static str> {
        perl_skip_route(&rec(pos, "N", alt, info), crate::masks::MAX_SV_SIZE)
    }

    #[test]
    fn test_perl_numification_and_the_structural_end() {
        let cases: [(&str, f64); 8] = [
            ("300", 300.0),
            ("-300", -300.0),
            ("300abc", 300.0),
            (" 3e2", 300.0),
            ("abc", 0.0),
            ("", 0.0),
            (".", 0.0),
            ("300.7", 300.7),
        ];
        for (text, value) in cases {
            assert_eq!(
                perl_numeric_value(text),
                PerlNumber::Float(value),
                "{text:?}"
            );
        }
        assert_eq!(perl_numeric_value("1e400"), SATURATED);
        assert_eq!(perl_numeric_value("-1e400"), -SATURATED);
        assert_eq!(perl_numeric_value("1.7e308"), PerlNumber::Float(1.7e308));
        let end = |svlen, end, breakpoint| perl_structural_end(1000, svlen, end, breakpoint);
        let saturated = end(Some("1e400"), None, false).expect("an end");
        assert!(matches!(saturated, PerlNumber::Integer(_)));
        assert_eq!(saturated, 1000 + SATURATED);
        assert_eq!(end(None, Some("1e400"), false).expect("an end"), SATURATED);
        assert_eq!(
            end(Some("-300"), Some("2000"), false).expect("an end"),
            1300
        );
        assert_eq!(end(Some("abc"), Some("2000"), false).expect("an end"), 1000);
        assert_eq!(end(None, Some("2000"), false).expect("an end"), 2000);
        assert_eq!(
            end(None, Some("2000.5"), false),
            Some(PerlNumber::Float(2000.5))
        );
        assert_eq!(end(None, Some("abc"), false), None);
        assert_eq!(end(None, Some(""), false), None);
        assert_eq!(end(None, Some("2000"), true).expect("an end"), 1001);
        assert_eq!(end(None, Some("0"), true).expect("an end"), 0);
        assert_eq!(end(None, Some(""), true), None);
        assert_eq!(end(Some("100"), Some("0"), true).expect("an end"), 1100);
    }

    #[test]
    fn an_end_is_a_double_from_a_text_and_an_integer_from_the_start() {
        let end = |svlen, end, breakpoint| perl_structural_end(1000, svlen, end, breakpoint);
        assert!(matches!(
            end(None, None, false),
            Some(PerlNumber::Integer(1001))
        ));
        assert!(matches!(
            end(None, None, true),
            Some(PerlNumber::Integer(1001))
        ));
        assert!(matches!(
            end(Some("300"), None, false),
            Some(PerlNumber::Float(_))
        ));
        assert!(matches!(
            end(Some(""), None, true),
            Some(PerlNumber::Float(_))
        ));
        assert_eq!(end(Some("5,9000"), None, false).expect("an end"), 1005);
        assert_eq!(end(Some(",9000"), None, false).expect("an end"), 1000);
        assert_eq!(end(None, Some(" 2000 "), false).expect("an end"), 2000);
        assert_eq!(end(None, Some(" 2000"), false).expect("an end"), 2000);
        // The blanks of another script pass the test of the stripped text, and the text as
        // written then numifies to 0.
        assert_eq!(end(None, Some("\u{a0}2000"), false).expect("an end"), 0);
        assert_eq!(end(None, Some("00"), true).expect("an end"), 1001);
    }

    #[test]
    fn test_symbolic_types_with_and_without_a_term() {
        assert_eq!(
            alleles(&rec(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL")),
            sorted(&["deletion"])
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<CPX>", "END=2000;SVTYPE=CPX")),
            sorted(&["CPX", "<CPX>"])
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<NON_REF>", "END=1100")),
            sorted(&["NON_REF", "<NON_REF>"])
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<DEL>,<DUP>", "END=2000")),
            sorted(&["deletion", "duplication", "copy_number_variation"])
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<INS>,<DEL>", "END=2000")),
            sorted(&["insertion", "deletion"])
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<DEL>,T", "END=2000")),
            sorted(&["deletion", "T"])
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<INS:ME:ALU>", "SVLEN=300")),
            sorted(&["Alu_insertion"])
        );
    }

    #[test]
    fn test_breakend_records() {
        assert_eq!(
            alleles(&rec(
                18457512,
                "N",
                "<BND>",
                "SVTYPE=BND;CHR2=21;POS2=32484691;END2=32484692;END=18457513"
            )),
            sorted(&[
                "chromosome_breakpoint",
                "N.",
                ".N",
                "N[21:32484692[",
                "N[21:32484691[",
                "N[21:18457513["
            ])
        );
        assert_eq!(
            alleles(&rec(100, "A", "A[chr22:5000[", "SVTYPE=BND")),
            sorted(&[
                "A[22:5000[",
                "N.",
                ".N",
                "A.",
                ".A",
                "chromosome_breakpoint"
            ])
        );
        assert_eq!(
            alleles(&rec(100, "T", ".T", "SVTYPE=BND")),
            sorted(&[".T", "chromosome_breakpoint", "N.", ".N", "T."])
        );
    }

    #[test]
    fn test_an_svtype_of_zero_is_false_as_perl_reads_it() {
        let record = rec(1000, "N", "<BND>,A", "SVTYPE=0");
        let alleles = record_alleles(&record);
        assert!(!alleles.contains("<0>"));
        assert!(!alleles.contains("<BND>/A"));
        assert!(!alleles.contains("0"));
        assert_eq!(route(1000, "<BND>,A", "SVTYPE=0"), Some("unsupported_type"));
    }

    #[test]
    fn test_start_is_the_base_after_pos_and_end_follows_get_end() {
        assert_eq!(
            location(&rec(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL")),
            at(1001, 2000)
        );
        assert_eq!(
            location(&rec(1000, "N", "<DEL>", "END=2000;SVLEN=-500;SVTYPE=DEL")),
            at(1001, 1500),
            "SVLEN before END"
        );
        assert_eq!(
            location(&rec(30165850, "N", "<DEL>", "SVTYPE=DEL")),
            at(30165851, 30165851)
        );
        assert_eq!(
            location(&rec(
                18457512,
                "N",
                "<BND>",
                "END=18457513;SVLEN=14027179;CHR2=21;END2=32484692"
            )),
            at(18457513, 32484691)
        );
        assert_eq!(
            location(&rec(
                18457512,
                "N",
                "<BND>",
                "END=18457513;CHR2=21;END2=32484692"
            )),
            at(18457513, 18457513),
            "a breakend's END is discarded"
        );
    }

    #[test]
    fn test_a_symbolic_bnd_without_svtype_discards_its_end() {
        assert_eq!(
            location(&rec(1000, "N", "<BND>", "END=5000")),
            at(1001, 1001)
        );
        assert_eq!(
            location(&rec(1000, "N", "<BND>", "END=5000;SVTYPE=0")),
            at(1001, 1001),
            "a false SVTYPE leaves the ALT"
        );
        assert_eq!(
            location(&rec(1000, "N", "N[22:5000[", "END=5000;SVTYPE=DEL")),
            at(1001, 5000)
        );
        assert_eq!(
            location(&rec(1000, "N", "<DEL>", "END=5000;SVTYPE=BND")),
            at(1001, 5000)
        );
        assert_eq!(
            location(&rec(1000, "N", "<DEL>", "END=5000;SVTYPE=BND")),
            location(&rec(1000, "N", "<DEL>", "END=5000;SVTYPE=DEL"))
        );
    }

    #[test]
    fn a_saturated_or_fractional_reading_gives_an_integer_end() {
        assert_eq!(
            location(&rec(1000, "N", "<DUP>", "SVTYPE=DUP;SVLEN=1e400")),
            at(1001, 1000 + SATURATED)
        );
        assert_eq!(
            location(&rec(1000, "N", "<DUP>", "SVTYPE=DUP;END=1e400")),
            at(1001, SATURATED)
        );
        assert_eq!(
            location(&rec(1000, "N", "<DUP>", "SVTYPE=DUP;END=-1e400")),
            at(1001, -SATURATED)
        );
        assert_eq!(
            location(&rec(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=300.7")),
            at(30425651, 30425950)
        );
    }

    #[test]
    fn test_the_four_skip_routes_and_the_kept_shapes() {
        assert_eq!(
            route(
                18457512,
                "<BND>",
                "END=18457513;SVTYPE=BND;SVLEN=14027179;CHR2=21;END2=32484692"
            ),
            Some("oversize")
        );
        assert_eq!(
            route(1000, "<DEL>", "END=12000000;SVTYPE=DEL"),
            Some("oversize")
        );
        assert_eq!(
            route(1000, "<DEL>", "END=10001001;SVTYPE=DEL"),
            None,
            "a span of exactly 10,000,000 is not over the cap"
        );
        assert_eq!(
            route(1000, "<DUP>", "SVTYPE=DUP;SVLEN=1e400"),
            Some("oversize")
        );
        assert_eq!(
            route(1000, "<DUP>", "SVTYPE=DUP;END=1e400"),
            Some("oversize")
        );
        assert_eq!(
            route(45716815, "<CPX>", "END=45730272;SVTYPE=CPX"),
            Some("unsupported_type")
        );
        assert_eq!(
            route(70000, "<NON_REF>", "END=70100"),
            Some("unsupported_type")
        );
        assert_eq!(route(1000, "<CTX>", "SVTYPE=CTX"), Some("unsupported_type"));
        assert_eq!(
            route(30165800, "<INS>,<DEL>", "END=30166000"),
            Some("unsupported_type")
        );
        assert_eq!(
            route(30165820, "<INV>,<DUP>", "END=30166020"),
            Some("unsupported_type")
        );
        assert_eq!(
            route(1000, "<DEL>,T", "END=2000"),
            Some("unsupported_type"),
            "a mixed ALT list has no term"
        );
        assert_eq!(
            route(30165850, "<DEL>", "SVTYPE=DEL"),
            Some("incomplete_deletion")
        );
        assert_eq!(
            route(2000, "<DEL>", "END=1500;SVTYPE=DEL"),
            Some("incomplete_deletion")
        );
        for info in [
            "SVTYPE=DEL;SVLEN=abc",
            "SVTYPE=DEL;SVLEN=",
            "SVTYPE=DEL;SVLEN=.;END=2000",
            "SVTYPE=DEL;SVLEN=0;END=2000",
        ] {
            assert_eq!(
                route(1000, "<DEL>", info),
                Some("incomplete_deletion"),
                "{info}"
            );
        }
        for info in [
            "SVTYPE=DUP;SVLEN=abc",
            "SVTYPE=DUP;SVLEN=",
            "SVTYPE=DUP;SVLEN=.",
            "SVTYPE=DUP;END=1000",
        ] {
            assert_eq!(route(1000, "<DUP>", info), None, "{info}");
        }
        for (alt, info) in [
            ("<DUP>", "SVTYPE=DUP;END=abc"),
            ("<DUP>", "SVTYPE=DUP;END=0"),
            ("<DUP>", "SVTYPE=DUP;END=900"),
            ("<BND>", "SVTYPE=BND;END=0"),
            ("<BND>", "SVTYPE=BND;END="),
            ("<DUP>", "SVTYPE=DUP;END="),
        ] {
            assert_eq!(route(1000, alt, info), Some("invalid_span"), "{alt} {info}");
        }
        assert_eq!(
            route(1000, "<DEL>", "SVTYPE=DEL;END=abc"),
            Some("incomplete_deletion")
        );
        assert_eq!(
            route(1000, "<DUP>", "SVTYPE=DUP;END=1500.5"),
            None,
            "a fractional END is annotated"
        );
        assert_eq!(route(1000, "G,<NON_REF>", "."), Some("gvcf_variant_site"));
        assert_eq!(
            perl_skip_route(
                &rec(1000, "A", "G,<NON_REF>", "."),
                crate::masks::MAX_SV_SIZE
            ),
            Some("gvcf_variant_site")
        );
        assert_eq!(
            route(1000, "G", "SVTYPE=DEL;END=2000"),
            None,
            "an N reference base routes the record to the structural path"
        );
        assert_eq!(
            route(1000, "<del>", "SVTYPE=DUP;END=2000"),
            None,
            "a lowercase ALT is outside the VCF 4.4 forms, so SVTYPE names the type"
        );
        assert_eq!(route(1000, "<del>", "."), Some("unsupported_type"));
        for (alt, info) in [
            ("<DEL>,<DUP>", "END=2000"),
            ("<DEL>", "END=2000;SVTYPE=DEL"),
            ("<CN0>", "END=2000"),
            ("<CN0>,<CN2>", "END=2000"),
            ("<INS:ME:ALU>", "SVLEN=300"),
            ("<DUP:TANDEM>", "END=5000"),
            ("<CNV:TR>", "END=1050;SVLEN=50;RUS=AC;RUC=10"),
            ("N[22:50000000[", "SVTYPE=BND"),
            ("N[21:50000000[", "SVTYPE=BND;END=50000000"),
            ("<INS>,<INS>", "SVLEN=300"),
        ] {
            assert_eq!(route(1000, alt, info), None, "{alt}");
        }
        assert_eq!(
            perl_skip_route(&rec(1000, "A", "T", "AF=0.1"), crate::masks::MAX_SV_SIZE),
            None,
            "a sequence variant takes the sequence path"
        );
        assert_eq!(
            perl_skip_route(&rec(1000, "AT", "A", "."), crate::masks::MAX_SV_SIZE),
            None
        );
    }

    #[test]
    fn test_so_abbreviations_follow_get_so_term() {
        for (sv_type, abbrev) in [
            ("<INS:ME:L1>", "INS_LINE1"),
            ("<DEL:ME:SVA>", "DEL_SVA"),
            ("<INS:ME>", "INS_ME"),
            ("<DUP:TANDEM>", "TDUP"),
            ("<CNV:TR>", "TREP"),
            ("<CN3>", "CNV"),
            ("<CN0>", "DEL"),
            ("<CN2>", "DUP"),
            ("<CNV>", "CNV"),
            ("]13:123456]T", "BND"),
            (".N", "BND"),
            ("<INV>", "INV"),
            ("DUP", "DUP"),
            ("<INS>/<INS>", "INS"),
            ("<DEL>/<DUP:TANDEM>", "CNV"),
            ("<DUP:INT>", "DUP"),
        ] {
            assert_eq!(perl_so_abbrev(sv_type), Some(abbrev), "{sv_type}");
        }
        for sv_type in ["<CPX>", "<NON_REF>", "CTX", "<INS>/<DEL>", "<DEL>/T"] {
            assert_eq!(perl_so_abbrev(sv_type), None, "{sv_type}");
        }
    }

    #[test]
    fn an_abbreviation_reads_letters_as_the_reference_does() {
        for (sv_type, abbrev) in [
            ("<ins:me:alu>", Some("INS_ALU")),
            ("<del:me:Herv>", Some("DEL_HERV")),
            ("<INS:ME:FOO>", Some("INS_ME")),
            ("<INS:MEALU>", Some("INS_ALU")),
            ("<\u{131}ns:me:l1>", Some("INS_LINE1")),
            ("<in\u{17f}:me:\u{17f}va>", Some("INS_SVA")),
            ("<\u{130}NS:ME>", None),
            ("<INS:ME:\u{212a}>", Some("INS_ME")),
            ("<cn2>", Some("CNV")),
            ("<CN=0>", Some("DEL")),
            ("CN2", Some("DUP")),
            ("<CN2>\n", Some("DUP")),
            ("<CN20>", Some("CNV")),
            ("<cnv:tr>", Some("TREP")),
            ("N.\n", Some("BND")),
            ("<DEL:FOO:BAR>", Some("DEL")),
            ("<DEL:>", None),
            ("DEL>", Some("DEL")),
            ("<D<E>L>", Some("DEL")),
            ("<DEL> / DEL", None),
            ("\u{a0}<DEL>\u{a0}/<del>", None),
            ("<DUP>/<dup>", Some("DUP")),
            ("/", None),
            ("", None),
        ] {
            assert_eq!(perl_so_abbrev(sv_type), abbrev, "{sv_type:?}");
        }
        assert_eq!(perl_record_type("<DEL>", Some("BND")), "<DEL>");
        assert_eq!(perl_record_type("<del>", Some("BND")), "BND");
        assert_eq!(perl_record_type("CN=12x", Some("BND")), "CN=12x");
        assert_eq!(perl_record_type("<FOO>", Some("0")), "<FOO>");
        assert_eq!(perl_record_type("<FOO>", Some("")), "<FOO>");
        assert_eq!(perl_record_type("<FOO>", None), "<FOO>");
    }

    #[test]
    fn test_the_class_term_vep_rs_prints_for_an_unsupported_alt() {
        for (alt, term) in [
            ("<del>", "deletion"),
            ("<Del>", "deletion"),
            ("<dEL>", "deletion"),
            ("<DEL:>", "deletion"),
            ("<DEL >", "deletion"),
            ("<DEL\u{a0}>", "deletion"),
            ("<DEL\u{200b}>", "deletion"),
            ("<Ins>", "insertion"),
            ("<INS:>", "insertion"),
            ("<Bnd>", "chromosome_breakpoint"),
            ("<CN=X>", "copy_number_variation"),
            ("<CNX>", "copy_number_variation"),
            ("<CN>", "copy_number_variation"),
            ("<Dup>", "duplication"),
            ("<Inv>", "inversion"),
            ("<<DEL>>", "<<DEL>>"),
            ("<:DEL>", "<:DEL>"),
            ("<CPX>", "<CPX>"),
            ("<CPX:X>", "<CPX:X>"),
            ("<D\u{c9}L>", "<D\u{c9}L>"),
            ("<\u{ff24}\u{ff25}\u{ff2c}>", "<\u{ff24}\u{ff25}\u{ff2c}>"),
            ("< DEL >", "< DEL >"),
        ] {
            assert_eq!(vep_rs_unsupported_alt_allele(alt, None), term, "{alt}");
        }
        assert_eq!(
            vep_rs_unsupported_alt_allele("<FOO>", Some("DEL")),
            "deletion",
            "INFO/SVTYPE names the class of a type the ALT does not"
        );
        assert_eq!(
            vep_rs_unsupported_alt_allele("<FOO>", Some("foo")),
            "foo",
            "INFO/SVTYPE is the type string when it names no class"
        );
        assert_eq!(vep_rs_unsupported_alt_allele("<CPX>", Some("CPX")), "CPX");
        assert_eq!(
            vep_rs_unsupported_alt_allele("<FOO>", Some("0")),
            "<FOO>",
            "a false SVTYPE names nothing, as Perl reads `0`"
        );
        assert_eq!(
            vep_rs_unsupported_alt_allele("<DEL:X>", Some("CPX")),
            "deletion",
            "a VCF 4.4 form keeps its own class"
        );
        assert_eq!(
            alleles(&rec(30425650, "N", "<Del>", "END=30425949")),
            sorted(&["Del", "<Del>", "deletion"])
        );
        assert_eq!(
            alleles(&rec(30425650, "N", "<CN=X>", "END=30425949")),
            sorted(&["CN=X", "<CN=X>", "copy_number_variation"])
        );
        assert_eq!(
            alleles(&rec(30425650, "N", "<:DEL>", "END=30425949")),
            sorted(&["", "<>", "<:DEL>", "-"]),
            "the tab format writes the empty Allele as `-`"
        );
        assert_eq!(
            alleles(&rec(30425650, "N", "<DEL>,<FOO>", "END=30425949")),
            sorted(&["deletion", "FOO", "<FOO>", "<DEL>/<FOO>"]),
            "the typeless ALT's rows print the list as written, the DEL's its class"
        );
        assert_eq!(
            alleles(&rec(
                30425650,
                "N",
                "<FOO>,<BAR>",
                "SVTYPE=baz;END=30425949"
            )),
            sorted(&["FOO", "<FOO>", "BAR", "<BAR>", "baz"]),
            "INFO/SVTYPE names a list that opens outside the VCF 4.4 forms"
        );
        assert_eq!(
            alleles(&rec(30425650, "N", "<<DEL>>,<DEL>", "END=30425949")),
            sorted(&["deletion", "<<DEL>>", "<<DEL>>/<DEL>"])
        );
        assert_eq!(
            alleles(&rec(30425650, "N", "<DEL>,T", "SVTYPE=INS;END=30425949")),
            sorted(&["deletion", "insertion", "T"]),
            "a mixed list under a type the table lacks takes INFO/SVTYPE's class in vep-rs"
        );
        assert_eq!(
            alleles(&rec(1000, "N", "<CPX>", "END=2000;SVTYPE=CPX")),
            sorted(&["CPX", "<CPX>"]),
            "the raw spelling stays"
        );
    }

    #[test]
    fn a_record_is_on_the_sequence_path_by_its_letters_or_its_forms() {
        let path = |reference: &str, alt: &str, info: &str| {
            perl_sequence_path(&rec(1000, reference, alt, info))
        };
        assert!(path("A", "T", "SVTYPE=DEL"));
        assert!(path("A", "T,", "SVTYPE=DEL"));
        assert!(!path("A", "t", "SVTYPE=DEL"));
        assert!(!path("N", "T", "SVTYPE=DEL"));
        assert!(path("N", "T", "SVTYPE=0"));
        assert!(path("N", "T", "SVTYPE"));
        assert!(path("N", "*", "."));
        assert!(path("N", "<*>", "."));
        assert!(!path("N", "<DEL>", "."));
        assert!(!path("N", "N[5:10[", "."));
        assert!(!path("N", ".N", "."));
        assert!(!path("N", "N.", "."));
        assert!(!path("N", "N.\n", "."));
        assert!(path("N", ".", "."));
        assert!(path("", "", "."));
        assert!(!path("N", "G,N.", "."));
        assert!(path("N", "N.,G", "."));
        assert!(!path("N", ".\u{b2}", "."));
        assert!(path("N", ".\u{301}", "."));
        let site = |alt: &str, info: &str| route(1000, alt, info);
        assert_eq!(site("G,<non_ref>", "."), Some("gvcf_variant_site"));
        assert_eq!(site("<NON_REF>,.\u{b2}", "."), Some("unsupported_type"));
        assert_eq!(site("<NON_REF>,.\u{301}", "."), Some("gvcf_variant_site"));
        assert_eq!(site("<NON_REF>,*", "."), Some("unsupported_type"));
        assert_eq!(site("<NON_REF>,.", "."), None);
        assert_eq!(
            reference_start(&rec(1000, "A", "T", ".")),
            (String::from("21"), 1000)
        );
        assert_eq!(
            reference_start(&rec(1000, "N", "<DEL>", ".")),
            (String::from("21"), 1001)
        );
    }

    #[test]
    fn a_location_is_normalised_and_a_span_is_printed_low_end_first() {
        assert_eq!(span_location("21", 5, 5), "21:5");
        assert_eq!(span_location("chr21", 5, 9), "21:5-9");
        assert_eq!(span_location("21", 9, 5), "21:9-5");
        assert_eq!(span_location("chrM", 5, 5), "MT:5");
        assert_eq!(span_location("m", 5, 5), "MT:5");
        assert_eq!(span_location("CHR5", -4, 295), "5:-4-295");
        assert_eq!(printed_location("21", 9, 5), "21:5-9");
        assert_eq!(
            normalize_location("21:30425651-030425949"),
            "21:30425651-30425949"
        );
        assert_eq!(
            normalize_location("21:30425651-+30425949"),
            "21:30425651-30425949"
        );
        assert_eq!(
            normalize_location("21:30425651-30425949.5"),
            "21:30425651-30425949.5"
        );
        assert_eq!(normalize_location("chr21"), "chr21");
        assert_eq!(
            normalize_allele("N[chr22:5[,]chr1:7]N", "vep-rs"),
            "N[22:5[,]1:7]N"
        );
        assert_eq!(normalize_allele("N[Chr22:5[", "vep-rs"), "N[Chr22:5[");
        let spans = [("21", 100, 200), ("22", 150, 150)];
        assert!(reference_row_beside(&spans, "21", 201, 300));
        assert!(!reference_row_beside(&spans, "21", 202, 300));
        assert!(reference_row_beside(&spans, "21", 50, 99));
        assert!(!reference_row_beside(&spans, "21", 50, 98));
        assert!(!reference_row_beside(&spans, "22", 100, 148));
        assert!(!reference_row_beside(&spans, "X", 100, 200));
    }

    #[test]
    fn the_rows_vep_rs_prints_for_a_record() {
        let rows = |v: &VcfVariant| -> Vec<(String, String)> {
            vep_rs_record_rows(v).into_iter().collect()
        };
        let pair = |location: &str, allele: &str| (location.to_string(), allele.to_string());
        assert_eq!(
            rows(&rec(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=0")),
            vec![pair("21:30425650-30425651", "deletion")]
        );
        assert_eq!(
            rows(&rec(30425650, "N", "<DEL>", "SVTYPE=DEL;END=30425949")),
            vec![pair("21:30425651-30425949", "deletion")]
        );
        assert_eq!(
            rows(&rec(30425650, "ATG", "AG,<NON_REF>", ".")),
            vec![
                pair("21:30425651", "<NON_REF>"),
                pair("21:30425651", "AG"),
                pair("21:30425651", "AG/<NON_REF>"),
                pair("21:30425651", "NON_REF"),
            ]
        );
        assert_eq!(
            rows(&rec(30425650, "ATG", "AG", ".")),
            vec![pair("21:30425651", "-")]
        );
        assert_eq!(
            rows(&testing::rec_on("chrM", 100, "A", "G", ".", "m")),
            vec![pair("MT:100", "G")]
        );
    }

    #[test]
    fn the_reference_rows_at_the_start_it_prints_for_a_record() {
        let numified = [
            t(
                "21:1-30425949",
                "deletion",
                "ENST00000262354",
                "transcript_ablation",
            ),
            t(
                "21:31-30425949",
                "deletion",
                "ENST00000262354",
                "transcript_ablation",
            ),
            intergenic("21:-4-295", "deletion"),
            intergenic("21:-5", "G"),
        ];
        let shared = [intergenic("21:1-300", "deletion")];
        let perl: Vec<[String; 5]> = numified.iter().chain(&shared).cloned().collect();
        let owned = OwnedSides::from_rows(&rows(&perl), &rows(&shared));
        let located = LocatedSet::new(&owned.perl, &owned.strings);
        let at = |v: &VcfVariant| -> TupleSet {
            reference_rows_starting_at(v, &located, &owned.rust)
                .into_iter()
                .collect()
        };
        let deletion = "SVTYPE=DEL;END=30425949";
        assert_eq!(
            at(&rec(0, "N", "<DEL>", deletion)),
            owned.set(&rows(&numified[..1]))
        );
        assert_eq!(
            at(&rec(30, "N", "<DEL>", deletion)),
            owned.set(&rows(&numified[1..2]))
        );
        assert_eq!(
            at(&rec(-5, "N", "<DEL>", "SVTYPE=DEL;SVLEN=-300")),
            owned.set(&rows(&numified[2..3]))
        );
        assert_eq!(
            at(&rec(-5, "A", "G", ".")),
            owned.set(&rows(&numified[3..]))
        );
        assert!(at(&rec(7, "A", "G", ".")).is_empty());
        assert!(at(&testing::rec_on("22", 0, "N", "<DEL>", deletion, "r")).is_empty());
        assert!(located.starts_at("21", 1));
        assert!(!located.starts_at("21", 2));
        assert!(!located.starts_at("22", 1));
        assert_eq!(located.locations_at("21", 1).len(), 2);
        let mut spans = located.spans();
        spans.sort_unstable();
        assert_eq!(
            spans,
            vec![
                ("21", -5, -5),
                ("21", -4, 295),
                ("21", 1, 300),
                ("21", 1, 30425949),
                ("21", 31, 30425949),
            ]
        );
    }
}
