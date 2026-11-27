// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The structural-variant classes.
//!
//! These classes are decided with the predicates of the structural-variant comparison
//! (`tools/vep-sv-compare`); this module holds those predicates and what they call:
//! `mate_context` derives what a breakend's mate-side row supports from the cache,
//! `outputs` holds the readers of a Location, an allele and an INFO column, and
//! `number` the whole numbers and floats their arithmetic keeps apart.

pub mod mate_context;
pub mod number;
pub mod outputs;
pub mod zcat;

use crate::argparse::path_text;
use crate::cache;
use crate::classes;
use crate::classes::REFERENCE_SKIPPED_RECORD;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::Error;
use crate::error::Result;
use crate::input;
use crate::pyfmt;
use crate::pyjson::Value;
use crate::report::strerror;
use crate::rules::whole_text;
use crate::shared::canonical_location;
use crate::shared::consequence_set;
use crate::shared::dash;
use crate::shared::gvcf_site_rows;
use crate::shared::is_structural_allele;
use crate::shared::perl_fields;
use crate::shared::perl_split;
use crate::stream::lines::lines;
use crate::stream::source::{Buffers, Source, BLOCK};
use indexmap::IndexMap;
use number::minus;
use number::plus;
use number::whole_number;
use number::Int;
use number::Number;
use outputs::_info_written;
use outputs::_location_span;
use outputs::_non_integer_end;
use outputs::dot_breakend_points;
use outputs::is_breakend_allele;
use outputs::is_gvcf_variant_site;
use outputs::is_own_breakend_allele;
use outputs::perl_true;
use outputs::structural_alt_list;
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::sync::LazyLock;

pub use mate_context::iter_cache_transcripts;
pub use mate_context::read_cache_transcripts;

/// A transcript as the cache files it, with the genomic CDS bounds resolved.
#[derive(Debug)]
pub struct TranscriptModel {
    pub chrom: String,
    pub start: Int,
    pub end: Int,
    pub strand: Int,
    pub biotype: String,
    pub exons: Vec<(Int, Int)>,
    pub cds: Option<(Int, Int)>,
    pub mirna_mature: Vec<(Int, Int)>,
}

/// Membership only.
pub const _ACGT: &[&str] = &["A", "C", "G", "T"];

/// The `--scored-engine` value that turns the symbolic-allele mapping on. Compared
/// case-insensitively; every other engine name (the default `vep-rs`, `perl`) leaves
/// the Allele column as written.
pub const FASTVEP_ENGINE: &str = "fastvep";

/// Mobile-element subtypes Ensembl VEP names in the Allele column, keyed by the
/// upper-cased third segment of `<INS:ME:...>` / `<DEL:ME:...>`. Mirrors
/// `sv_type_abbreviation` in crates/vep-core/src/sv_routing.rs, which mirrors the
/// `@mobile_elements` list and the `L1 -> LINE1` alias in Perl's
/// `Parser.pm::get_SO_term`. Any other subtype falls back to the generic
/// `mobile_element_<operation>`, as both engines do.
pub const _ME_SUBTYPE_NAMES: &[(&str, &str)] = &[
    ("ALU", "Alu"),
    ("L1", "LINE1"),
    ("LINE1", "LINE1"),
    ("SVA", "SVA"),
    ("HERV", "HERV"),
];

/// Comparison tuple: (Location, Allele, Feature, Feature_type, Consequence)
pub type AnnotTuple = (String, String, String, String, String);

/// Transcript span type: (chromosome, start, end) as the cache files it.
pub type TranscriptSpan = (String, Int, Int);

/// Ensembl VEP's breakend-to-feature distance gate (ensembl-variation
/// StructuralVariationOverlap.pm:144-160, MAX_DISTANCE_FROM_TRANSCRIPT in
/// Utils/VariationEffect.pm:60).
pub const _BREAKEND_FEATURE_DISTANCE: i64 = 5000;

pub const _BRACKET_MATE_RE_PATTERN: &str = r"[\[\]]([^:\[\]]+):(\d+)[\[\]]";
pub static _BRACKET_MATE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(_BRACKET_MATE_RE_PATTERN).expect("the pattern compiles"));

/// vep-rs writes this literal string in the Allele column for <CNV:TR> records.
pub const _CNV_TR_RUST_ALLELE: &str = "tandem_repeat";

/// The reference's own skip routes, reproduced for the cross-check below. `_SO_TERMS` is
/// ensembl-variation Utils/Config.pm:1417-1438 at release 116; the abbreviation rules are
/// Parser.pm 116.2:663-728 (`get_SO_term`), the type derivation Parser/VCF.pm:468-481, the
/// routing of a record to the structural path Parser/VCF.pm:237-245, the coordinates
/// ensembl-io BaseVCF4.pm `get_start` / `get_end` (the base before the event is the
/// start; SVLEN sets the end before END does; a breakend's END is discarded).
pub const _SO_TERMS: &[(&str, &str)] = &[
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

pub const _VCF44_SYMBOLIC_PATTERN: &str = r"^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)";
pub static _VCF44_SYMBOLIC: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(_VCF44_SYMBOLIC_PATTERN).expect("the pattern compiles"));

pub const _PERL_NUMBER_PATTERN: &str =
    r"^[ \t\n\r\x0b\x0c]*[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?";
pub static _PERL_NUMBER: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(_PERL_NUMBER_PATTERN).expect("the pattern compiles"));

/// An infinite reading saturates to the machine word, as vep-rs's reading does; a span that
/// large is over every --max_sv_size, as Perl's `Inf` is.
pub const _SATURATED_MAGNITUDE: i64 = 9223372036854775807;

/// The spellings `looks_like_number` admits beside a decimal: Perl's NaN and its infinity, with
/// an optional sign and blanks around.
pub const _PERL_NAN_PATTERN: &str = r"(?i)[\s\x1c-\x1f]*[+-]?nan[\s\x1c-\x1f]*";
pub static _PERL_NAN: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(_PERL_NAN_PATTERN).expect("the pattern compiles"));

/// Under `re.IGNORECASE` an `i` of the pattern also takes the dotted capital and the dotless
/// small letter.
pub const _PERL_INFINITY_PATTERN: &str =
    r"(?i)[\s\x1c-\x1f]*[+-]?[iİı]nf(?:[iİı]n[iİı]ty)?[\s\x1c-\x1f]*";
pub static _PERL_INFINITY: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(_PERL_INFINITY_PATTERN).expect("the pattern compiles"));

pub const _OWN_END_ALLELES: &[&str] = &[".N", "N.", "chromosome_breakpoint"];

/// The class vep-rs gives a symbolic ALT by its own prefix, read without regard to case
/// (`classify_symbolic_alt` in crates/vep-core/src/sv_routing.rs), first match first; a
/// `<CN...>` with a copy number the table knows (`<CN0>`, `<CN2>`) is a supported type and
/// never reaches this list.
pub const _VEP_RS_ALT_PREFIX_TERMS: &[(&str, &str)] = &[
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

pub const _VEP_RS_SVTYPE_TERMS: &[(&str, &str)] = &[
    ("DEL", "deletion"),
    ("INS", "insertion"),
    ("DUP", "duplication"),
    ("INV", "inversion"),
    ("BND", "chromosome_breakpoint"),
    ("CNV", "copy_number_variation"),
];

static PERL_NUMBER_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(_PERL_NUMBER_PATTERN));
static PERL_NAN_WHOLE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(_PERL_NAN_PATTERN));
static PERL_INFINITY_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(_PERL_INFINITY_PATTERN));

/// The patterns `perl_so_abbrev` reads a type string with. Under `re.IGNORECASE` an `i`
/// and a class of the capital letters also take the dotted capital and the dotless small
/// letter; a `$` of `re` also stands before a last line end.
static MOBILE_ELEMENT: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"(?i)([iİı]NS|DEL):(ME):?([A-Z0-9İı]+)?").expect("the pattern compiles")
});
static NAMES_DEL: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"(?i)DEL").expect("the pattern compiles"));
static NAMES_DUP: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"(?i)DUP").expect("the pattern compiles"));
static NAMES_TANDEM_DUP: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"(?i)DUP:TANDEM").expect("the pattern compiles"));
static NAMES_TANDEM_REPEAT: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"(?i)CNV:TR").expect("the pattern compiles"));
static NAMES_COPY_NUMBER: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"(?i)CN=?[0-9]").expect("the pattern compiles"));
static COPY_NUMBER_ZERO: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"^<?CN=?0>?\n?\z").expect("the pattern compiles"));
static COPY_NUMBER_TWO: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"^<?CN=?2>?\n?\z").expect("the pattern compiles"));
static NAMES_CNV: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"(?i)CNV").expect("the pattern compiles"));
static BREAKEND_FORM: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"[\[\]]|^\.|\.\n?\z").expect("the pattern compiles"));
static BRACKETED_FORM: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"^<|>\n?\z").expect("the pattern compiles"));

/// The term of `_SO_TERMS` for an abbreviation.
fn so_term(abbrev: &str) -> Option<&'static str> {
    _SO_TERMS
        .iter()
        .find(|(name, _)| *name == abbrev)
        .map(|(_, term)| *term)
}

/// `re.sub(r":.+", "", text)`: every `:` with a character after it on its line goes, with
/// the rest of the line.
fn without_subtype(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find(':') {
        let after = &rest[at + 1..];
        if after.is_empty() || after.starts_with('\n') {
            out.push_str(&rest[..=at]);
            rest = after;
        } else {
            out.push_str(&rest[..at]);
            rest = after.find('\n').map_or("", |end| &after[end..]);
        }
    }
    out.push_str(rest);
    out
}

/// `chrom[3:] if chrom.lower().startswith("chr") else chrom`.
fn without_chr(chrom: &str) -> &str {
    match chrom.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &chrom[3..],
        _ => chrom,
    }
}

/// The Location both engines print for a span: `chrom:start`, or `chrom:lo-hi` with an
/// insertion's `end = start - 1` written first (`shared::sequence_location`, at a position
/// of any size).
fn sequence_location(chrom: &str, start: Int, end: Int) -> String {
    let (lo, hi) = (start.min(end), start.max(end));
    if lo == hi {
        format!("{chrom}:{lo}")
    } else {
        format!("{chrom}:{lo}-{hi}")
    }
}

/// `shared::gvcf_site_rows` at a position of any size.
fn site_rows(v: &VcfVariant) -> Vec<(String, Int, Int)> {
    gvcf_site_rows(0, &v.r#ref, &perl_split(&v.alt, ","))
        .into_iter()
        .map(|(allele, start, end)| {
            (
                allele,
                v.pos.saturating_add(Int::from(start)),
                v.pos.saturating_add(Int::from(end)),
            )
        })
        .collect()
}

/// The chromosome of a record's Location.
fn record_chromosome(v: &VcfVariant) -> &str {
    v.location.split(':').next().unwrap_or(&v.location)
}

/// Strip chr prefix from Location for comparison (vep-rs and Ensembl VEP both print a
/// `chr` prefix the input record carries),
/// and write each coordinate as its number (`canonical_location`: Ensembl VEP prints an
/// INFO/END as given, `030425949` or `+30425949`).
pub fn normalize_location(loc: &str) -> String {
    let loc = canonical_location(loc);
    let Some((chrom, rest)) = loc.split_once(':') else {
        return loc;
    };
    // Strip chr prefix (case-insensitive)
    let chrom = without_chr(chrom);
    // Ensembl convention: M -> MT
    let chrom = if chrom.eq_ignore_ascii_case("M") {
        "MT"
    } else {
        chrom
    };
    format!("{chrom}:{rest}")
}

/// Strip chr prefix inside BND bracket notation for allele comparison.
///
/// Ensembl VEP and vep-rs both keep a `chr` prefix the record writes on a BND allele's
/// mate chromosome (e.g., N[chr21:1234[). Normalize every allele to the chr-stripped
/// form (e.g., N[21:1234[) so alleles match in concordance tuples under either spelling.
pub fn _strip_chr_in_bracket(allele: &str) -> String {
    let mut out = String::with_capacity(allele.len());
    let mut rest = allele;
    while let Some(at) = rest.find(['[', ']']) {
        out.push_str(&rest[..=at]);
        rest = &rest[at + 1..];
        rest = rest.strip_prefix("chr").unwrap_or(rest);
    }
    out.push_str(rest);
    out
}

pub fn _is_fastvep(engine: &str) -> bool {
    pyfmt::lower(pyfmt::strip(engine)) == FASTVEP_ENGINE
}

/// Map a raw symbolic VCF ALT token to the Allele string Ensembl VEP writes for it.
///
/// Ensembl VEP writes an SV's Allele column as its Sequence Ontology class term
/// (`OutputFactory.pm`: `$vfoa->{breakend}->{string} || $svf->class_SO_term`), and
/// vep-rs reproduces that through `VariantClass::so_term` plus `display_allele` in
/// crates/vep-core/src/variant.rs. fastVEP through v0.5.0 instead writes the raw VCF
/// token, so without this mapping every symbolic-allele tuple of such an output is
/// a guaranteed miss whatever its consequence set: `<DEL>` can never equal
/// `deletion`. From v0.6.0 fastVEP writes the class term itself, which the mapping
/// returns as written.
///
/// The rules below are vep-rs's own typing (`vep_core::sv_routing`: `sv_so_term_and_class`
/// over `sv_type_abbreviation` for a type with a Sequence Ontology term,
/// `classify_symbolic_alt` for one without) followed by its display rule, so the
/// mapping is what the scored vep-rs arm writes, the `<CPX>` row excepted:
///
///     token (case-insensitive)   vep-rs VariantClass       Allele written
///     -------------------------  ------------------------  --------------------------
///     <INS:ME>                   MobileElementInsertion    mobile_element_insertion
///     <INS:ME:ALU>               MobileElementInsertion    Alu_insertion
///     <INS:ME:L1>, <INS:ME:LINE1>  MobileElementInsertion  LINE1_insertion
///     <INS:ME:SVA>               MobileElementInsertion    SVA_insertion
///     <INS:ME:HERV>              MobileElementInsertion    HERV_insertion
///     <INS:ME:other>             MobileElementInsertion    mobile_element_insertion
///     <DEL:ME[:x]>               MobileElementDeletion     as above, with _deletion
///     <DEL...>                   StructuralDeletion        deletion
///     <INS...>                   StructuralInsertion       insertion
///     <DUP:TANDEM...>            TandemDuplication         tandem_duplication
///     <DUP...>                   Duplication               duplication
///     <INV...>                   Inversion                 inversion
///     <CNV:TR...>                TandemRepeat              tandem_repeat
///     <CN0> / <CN=0>             StructuralDeletion        deletion
///     <CN2> / <CN=2>             Duplication               duplication
///     <CNn> (any other n)        CopyNumberVariation       copy_number_variation
///     <CNV...>, <CN...>          CopyNumberVariation       copy_number_variation
///     <CPX[:x]>                  ComplexStructural         CPX (<CPX[:x]> as written
///                                                          without INFO/SVTYPE; both
///                                                          spellings are set aside)
///
/// Left UNCHANGED, deliberately:
///
/// * `<NON_REF>` and `<*>`: both engines write the literal token (Ensembl VEP has no SO
///   term for either; vep-rs `display_allele` returns the raw allele), so identity
///   is already the match.
/// * `<BND>`: Ensembl VEP synthesises bracket-notation and single-breakend alleles from
///   `INFO/CHR2`, `END2`/`END` and `SVLEN`, none of which a per-row output
///   token carries, so no token-level mapping can reproduce them.
/// * Any other `<...>` token: vep-rs classifies it by `INFO/SVTYPE` when present
///   and only falls back to the stripped token, and the comparator cannot see
///   INFO. A miss stays a miss rather than inventing an Ensembl VEP string.
/// * Everything that is not a `<...>` token: literal bases, breakend notation,
///   `-`, and the SO class strings themselves are returned as-is, so applying the
///   mapping to an already-normalised value is a no-op.
///
/// Two Ensembl VEP behaviours are out of reach of a per-row mapping and score as
/// divergences on the fastVEP arm: a multi-allelic `<CN0>,<CN2>` record, which
/// Ensembl VEP coalesces into ONE `copy_number_variation` feature, and `<CNV:TR>` run
/// with `--fasta`, where Ensembl VEP materialises the literal repeat sequence (the
/// class `filter_cnv_tr_expansion_divergences` masks for vep-rs).
pub fn fastvep_symbolic_allele_to_so_class(allele: &str) -> String {
    if pyfmt::len(allele) < 3 || !(allele.starts_with('<') && allele.ends_with('>')) {
        return allele.to_string();
    }
    let upper = pyfmt::upper(allele);
    let inner = &upper[1..upper.len() - 1];

    if matches!(upper.as_str(), "<*>" | "<NON_REF>" | "<BND>") {
        return allele.to_string();
    }

    // Mobile elements before generic DEL/INS: `<DEL:ME>` starts with `<DEL`.
    if upper.starts_with("<INS:ME") || upper.starts_with("<DEL:ME") {
        let parts: Vec<&str> = inner.split(':').collect();
        let operation = if parts[0] == "DEL" {
            "deletion"
        } else {
            "insertion"
        };
        let subtype = parts.get(2).and_then(|subtype| {
            _ME_SUBTYPE_NAMES
                .iter()
                .find(|(name, _)| name == subtype)
                .map(|(_, written)| *written)
        });
        return match subtype {
            Some(written) => format!("{written}_{operation}"),
            None => format!("mobile_element_{operation}"),
        };
    }
    for (prefix, term) in [
        ("<DEL", "deletion"),
        ("<INS", "insertion"),
        ("<DUP:TANDEM", "tandem_duplication"),
        ("<DUP", "duplication"),
        ("<INV", "inversion"),
        ("<CNV:TR", "tandem_repeat"),
    ] {
        if upper.starts_with(prefix) {
            return term.to_string();
        }
    }
    // `<CNn>` / `<CN=n>`: Ensembl VEP maps CN0 -> deletion and CN2 -> duplication, every
    // other copy number to the generic class. An unparseable number falls through
    // to the generic rule below, as vep-rs's `parse::<u32>()` failure does.
    if upper.starts_with("<CN") && !upper.starts_with("<CNV") {
        let num_str = &inner[2..];
        let num_str = num_str.strip_prefix('=').unwrap_or(num_str);
        if !num_str.is_empty() && num_str.bytes().all(|digit| digit.is_ascii_digit()) {
            return match whole_number(num_str) {
                0 => "deletion",
                2 => "duplication",
                _ => "copy_number_variation",
            }
            .to_string();
        }
    }
    if upper.starts_with("<CNV") || upper.starts_with("<CN") {
        return "copy_number_variation".to_string();
    }
    if upper.starts_with("<CPX") {
        // A type with no Sequence Ontology term prints Ensembl VEP's type string: the ALT
        // as written (`<CPX>`) or INFO/SVTYPE (`CPX`, Parser/VCF.pm:469-481). Ensembl VEP
        // skips such a record and `record_alleles` lists both spellings to set its rows
        // aside, so the fastVEP token maps to the bare abbreviation.
        return "CPX".to_string();
    }
    allele.to_string()
}

/// Normalize allele string for concordance comparison.
///
/// `engine` names the engine that wrote the row. Breakend bracket normalisation
/// applies to every engine; the symbolic-token mapping applies ONLY when `engine`
/// is fastVEP (`FASTVEP_ENGINE`), because vep-rs and Ensembl VEP already write the SO
/// class strings and mapping them again would be a no-op at best and a masking
/// hazard at worst.
///
/// `engine` is `vep-rs` where the caller names none.
pub fn normalize_allele(allele: &str, engine: &str) -> String {
    let allele = _strip_chr_in_bracket(allele);
    if _is_fastvep(engine) {
        return fastvep_symbolic_allele_to_so_class(&allele);
    }
    allele
}

/// Classify a VCF variant by REF/ALT into one of the defined types.
pub fn classify_variant(r#ref: &str, alt: &str) -> String {
    // Multi-allelic (more than one ALT field, a trailing empty one being none) -- check
    // first since other rules assume a single ALT allele.
    let alts = perl_split(alt, ",");
    if alts.len() > 1 {
        return "Multi_allelic".to_string();
    }
    let alt = alts.first().map_or(alt, String::as_str);

    // Spanning deletion
    if alt == "*" {
        return "Spanning".to_string();
    }

    // Reference-only
    if alt == "." {
        return "RefOnly".to_string();
    }

    // NON_REF / <*>
    if alt.starts_with("<NON_REF") || alt == "<*>" {
        return "NON_REF".to_string();
    }

    let alt_upper = pyfmt::upper(alt);
    // Symbolic alleles (angle-bracket notation)
    if alt.starts_with('<') {
        let upper = &alt_upper;
        let named = if upper.starts_with("<DEL") {
            "Symbolic_DEL"
        } else if upper.starts_with("<INS") {
            // <INS:ME:*> handled below, but plain <INS> first
            if upper.contains("ME") {
                "Mobile_Element"
            } else {
                "Symbolic_INS"
            }
        } else if upper.starts_with("<DUP") {
            "Symbolic_DUP"
        } else if upper.starts_with("<INV") {
            "Symbolic_INV"
        } else if upper.starts_with("<CNV") || upper.starts_with("<CN") {
            "Symbolic_CNV"
        } else if upper.starts_with("<CPX") {
            "Complex_SV"
        } else if upper.starts_with("<BND") {
            // A symbolic breakend (`<BND>` with INFO/CHR2 and END2) is a breakend record
            // like its bracket-notation form.
            "BND"
        } else {
            "Other"
        };
        return named.to_string();
    }

    // Breakend notation
    if alt.contains('[') || alt.contains(']') {
        return "BND".to_string();
    }
    // Single breakend: ALT is e.g. "A." or ".A"
    if let Some(bases) = alt.strip_suffix('.') {
        if pyfmt::is_alpha(bases) {
            return "BND".to_string();
        }
    }
    if let Some(bases) = alt.strip_prefix('.') {
        if pyfmt::is_alpha(bases) {
            return "BND".to_string();
        }
    }

    let ref_upper = pyfmt::upper(r#ref);
    let is_acgt = |text: &str| text.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'));

    if is_acgt(&ref_upper) && is_acgt(&alt_upper) {
        let (lr, la) = (ref_upper.len(), alt_upper.len());
        let named = if lr == 1 && la == 1 {
            "SNV"
        } else if lr == la && lr > 1 {
            "MNP"
        } else if lr > la {
            "Small_DEL"
        } else if lr < la {
            "Small_INS"
        } else {
            // Same length > 1 already covered by MNP
            "Complex"
        };
        return named.to_string();
    }

    "Other".to_string()
}

#[derive(Debug)]
pub struct VcfVariant {
    pub chrom: String,
    pub pos: Int,
    pub vid: String,
    pub r#ref: String,
    pub alt: String,
    pub variant_type: String,
    pub location: String,
    pub info: String,
    pub pos_text: Option<String>,
}

/// Whether `path.suffix` is `.gz`: the name ends in `.gz` after a character that is no
/// leading dot.
fn names_gzip(path: &Path) -> bool {
    let raw = path.as_os_str().as_bytes();
    let name = raw.rsplit(|byte| *byte == b'/').next().unwrap_or(raw);
    name.strip_suffix(b".gz")
        .is_some_and(|stem| stem.iter().any(|byte| *byte != b'.'))
}

/// The text of an input VCF as `parse_vcf` reads its lines: a file named `.gz` as `zcat`
/// writes it (`zcat::text`), a file that cannot be opened having no line; any other
/// opened as text.
///
/// The text is decoded as a stream, so a character that the end of the text cuts short
/// is refused at position 0, the place it has in the last piece decoded.
fn vcf_text(path: &Path) -> std::io::Result<String> {
    let raw = path.as_os_str().as_bytes();
    if raw.contains(&0) {
        crate::uncaught!("embedded null byte");
    }
    let bytes = if names_gzip(path) {
        // `zcat` takes a name that opens with `-` for its options and writes nothing.
        let file = (!raw.starts_with(b"-"))
            .then(|| std::fs::read(path).ok())
            .flatten();
        file.map(|file| zcat::text(&file)).unwrap_or_default()
    } else {
        std::fs::read(path)?
    };
    match String::from_utf8(bytes) {
        Ok(text) => Ok(input::universal_newlines(text)),
        Err(error) => {
            let refused = error.utf8_error();
            let piece = match refused.error_len() {
                Some(_) => 0,
                None => refused.valid_up_to(),
            };
            crate::uncaught!(
                "{}",
                pyfmt::decode_utf8(&error.as_bytes()[piece..])
                    .err()
                    .unwrap_or_default()
            )
        }
    }
}

/// Parse a VCF file (plain or bgzipped) and return classified variants. POS is read as
/// Ensembl VEP numifies it (`perl_pos`), so a record whose POS is not an integer is
/// measured, never raised on.
///
/// The error is the `OSError` of a file that cannot be opened or read.
pub fn parse_vcf(path: &Path) -> std::io::Result<Vec<VcfVariant>> {
    let text = vcf_text(path)?;
    Ok(text
        .split('\n')
        .filter(|line| !line.starts_with('#'))
        .filter_map(vcf_variant)
        .collect())
}

/// The record of one data line of the input; `None` for a line of fewer than five columns.
fn vcf_variant(line: &str) -> Option<VcfVariant> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 5 {
        return None;
    }
    let (chrom, pos_str, vid, r#ref, alt) = (parts[0], parts[1], parts[2], parts[3], parts[4]);
    let pos = perl_pos(pos_str);
    Some(VcfVariant {
        chrom: chrom.to_string(),
        pos,
        vid: vid.to_string(),
        r#ref: r#ref.to_string(),
        alt: alt.to_string(),
        variant_type: classify_variant(r#ref, alt),
        location: normalize_location(&format!("{chrom}:{pos}")),
        info: parts.get(7).copied().unwrap_or_default().to_string(),
        pos_text: Some(pos_str.to_string()),
    })
}

/// Each data line of the input with its ordinal, counted from 0 as the stream counts them
/// (`stream::Plan::force_records`): the lines the stream's own blocks cut, every line that
/// does not open with `#`, read once without holding the file.
fn each_input_record(input: &Path, mut each: impl FnMut(u64, &str)) -> crate::error::Result<()> {
    let source = Source::open(input)?;
    let buffers = Buffers::new(BLOCK);
    let mut blocks = source.blocks(&buffers)?;
    let mut ordinal = 0u64;
    while let Some(block) = blocks.next().map_err(|detail| source.failed(&detail))? {
        let bytes = block.bytes();
        for line in lines(bytes) {
            let line = &bytes[line];
            if line.first() != Some(&b'#') {
                each(ordinal, &String::from_utf8_lossy(line));
                ordinal += 1;
            }
        }
        buffers.give(block.into_buffer());
    }
    Ok(())
}

/// The chromosome of a record and the lowest and highest coordinate its rows can name.
fn record_reach(v: &VcfVariant) -> (String, Int, Int) {
    let (_chrom, start, end) = record_location(v);
    let (mut lo, mut hi) = (start.min(end).min(v.pos), start.max(end).max(v.pos));
    for (_allele, start, end) in site_rows(v) {
        (lo, hi) = (lo.min(start.min(end)), hi.max(start.max(end)));
    }
    (record_chromosome(v).to_string(), lo, hi)
}

/// The ordinals of the input's data lines whose groups the classes that read across records
/// must see: a record the reference may skip (`perl_skip_route`) under
/// `reference_skipped_record`, a breakend (a bracket, `<BND>`, a single breakend, or
/// `SVTYPE=BND`, whose rows carry bracket mates) under `breakend_mate_local_read` or
/// `bnd_own_chromosome_rows_lost`, and every record whose span lies within two bases of one
/// of those, whose rows those classes read beside it (the rivals at a start of
/// `skipped_record_rows`, the tuples at a Location of the two breakend classes). The rest
/// of a file's groups are touched by none of them. The input is read twice, the second time
/// only where it holds such a record.
pub fn records_read_across(
    input: &Path,
    skipped: bool,
    breakends: bool,
    max_sv_size: i64,
) -> crate::error::Result<Vec<u64>> {
    let mut seeds: IndexMap<String, Vec<(Int, Int)>> = IndexMap::new();
    each_input_record(input, |_, line| {
        let Some(v) = vcf_variant(line) else {
            return;
        };
        let breakend = v.variant_type == "BND"
            || v.alt.contains(['[', ']'])
            || v.info.split(';').any(|item| item == "SVTYPE=BND");
        if (skipped && perl_skip_route(&v, max_sv_size).is_some()) || (breakends && breakend) {
            let (chrom, lo, hi) = record_reach(&v);
            seeds
                .entry(chrom)
                .or_default()
                .push((minus(lo, 2), plus(hi, 2)));
        }
    })?;
    if seeds.is_empty() {
        return Ok(Vec::new());
    }
    // The seeds' spans merged per chromosome, so that one search finds a record's.
    for merged in seeds.values_mut() {
        merged.sort_unstable();
        let mut joined: Vec<(Int, Int)> = Vec::new();
        for (lo, hi) in merged.drain(..) {
            match joined.last_mut() {
                Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
                _ => joined.push((lo, hi)),
            }
        }
        *merged = joined;
    }
    let mut read = Vec::new();
    each_input_record(input, |ordinal, line| {
        let Some(v) = vcf_variant(line) else {
            return;
        };
        let (chrom, lo, hi) = record_reach(&v);
        let beside = seeds.get(&chrom).is_some_and(|merged| {
            let at = merged.partition_point(|span| span.1 < lo);
            merged.get(at).is_some_and(|span| span.0 <= hi)
        });
        if beside {
            read.push(ordinal);
        }
    })?;
    Ok(read)
}

/// A transcript's `stable_id` as a key of the cache's views: a text, or a number by its
/// value, as the interpreter holds a `stable_id` that is no text (`1`, `1.0` and `true` are
/// one key, every NaN another).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StableId {
    Text(String),
    /// A whole number, as its decimal text.
    Whole(String),
    /// A float that is no whole number, as its bits.
    Fraction(u64),
    Nan,
}

impl StableId {
    /// The key of a `stable_id` that reads as true. A list or an object is no key and ends
    /// the run as the interpreter's `TypeError`.
    fn of(stable_id: &Value) -> StableId {
        match stable_id {
            Value::Str(text) => StableId::Text(text.clone()),
            Value::Bool(flag) => StableId::Whole(i64::from(*flag).to_string()),
            Value::Int(whole) => StableId::Whole(whole.to_string()),
            Value::BigInt(whole) => StableId::Whole(whole.clone()),
            Value::Float(float) if float.is_nan() => StableId::Nan,
            Value::Float(float) if float.is_finite() && float.fract() == 0.0 => {
                StableId::Whole(pyfmt::int_of_f64(*float).unwrap_or_default())
            }
            Value::Float(float) => StableId::Fraction(float.to_bits()),
            Value::List(_) => {
                crate::uncaught!("cannot use 'list' as a set element (unhashable type: 'list')")
            }
            Value::Dict(_) => {
                crate::uncaught!("cannot use 'dict' as a set element (unhashable type: 'dict')")
            }
            Value::Null => crate::uncaught!("'stable_id'"),
        }
    }

    /// The key of a transcript named in an output.
    pub fn text(feature: &str) -> StableId {
        StableId::Text(feature.to_string())
    }
}

/// `(int(tx["start"]), int(tx["end"]))`, `None` where either is absent or no number; the
/// end is read after a start that is one.
fn span_of(tx: &Value) -> Option<(Int, Int)> {
    let start = mate_context::_int(tx.get("start"))?;
    let end = mate_context::_int(tx.get("end"))?;
    Some((start, end))
}

/// The three views `load_cache_authority` gives of a cache: the transcript ids of each
/// chromosome, and each transcript's span and structure by chromosome and id.
pub type CacheAuthority = (
    IndexMap<String, BTreeSet<StableId>>,
    IndexMap<(String, StableId), TranscriptSpan>,
    IndexMap<(String, StableId), TranscriptModel>,
);

/// Read the vep-rs JSON cache once and return three views of it.
///
/// The first maps each chromosome to the transcript ids the cache holds for it; the
/// second maps `(chromosome, transcript id)` to the transcript's `(chromosome, start,
/// end)`; the third maps the same key to the transcript's structure (`TranscriptModel`:
/// exons, strand, biotype, genomic CDS bounds, mature miRNA ranges). The chromosome is
/// part of the key because the pseudoautosomal transcripts are filed under both X and Y
/// with their own coordinates on each. Both engines are
/// run against this one cache, so a transcript it files under chr21 is a chr21
/// transcript and its structure is the structure both engines annotated against,
/// whatever either engine's output says.
///
/// Returns three empty mappings when the directory is absent or unreadable, and the
/// callers treat that as "cannot judge" rather than "nothing to exclude".
pub fn load_cache_authority(cache_dir: &Path) -> CacheAuthority {
    let mut by_chr: IndexMap<String, BTreeSet<StableId>> = IndexMap::new();
    let mut spans = IndexMap::new();
    let mut models = IndexMap::new();
    let read = read_cache_transcripts(cache_dir, |chrom, tx| {
        let id = StableId::of(tx.get("stable_id").unwrap_or(&Value::Null));
        let span = span_of(&tx);
        let model = mate_context::model_from_cache_record(chrom, &tx);
        Some((chrom.to_string(), id, span, model))
    });
    for (chrom, id, span, model) in read {
        by_chr.entry(chrom.clone()).or_default().insert(id.clone());
        if let Some((start, end)) = span {
            spans.insert((chrom.clone(), id.clone()), (chrom.clone(), start, end));
        }
        if let Some(model) = model {
            models.insert((chrom, id), model);
        }
    }
    (by_chr, spans, models)
}

/// The chromosome a VEP Location names, normalised to the cache's convention.
pub fn _chromosome_of_location(loc: &str) -> Option<String> {
    let (chrom, _) = loc.split_once(':')?;
    Some(without_chr(chrom).to_string())
}

/// The mate `(chromosome, position)` a bracket-notation breakend allele names.
pub fn _bracket_mate(allele: &str) -> Option<(String, Int)> {
    let found = _BRACKET_MATE_RE.captures(allele)?;
    Some((without_chr(&found[1]).to_string(), whole_number(&found[2])))
}

/// The tuples one side wrote and the other did not, in the tuples' order.
fn one_sided<'a>(
    own: &'a BTreeSet<AnnotTuple>,
    other: &'a BTreeSet<AnnotTuple>,
) -> impl Iterator<Item = &'a AnnotTuple> {
    own.iter().filter(move |t| !other.contains(*t))
}

/// A point breakend's mate-side Transcript row, read by Ensembl VEP at the local coordinate.
///
/// THE SHAPE. Ensembl VEP builds a Transcript row for every transcript within 5 kb of a
/// breakend's mate on the mate chromosome (StructuralVariationOverlap.pm:86-102) and then
/// evaluates every positional predicate on it with the LOCAL variation feature
/// (BaseVariationFeatureOverlapAllele.pm:257,273); only `feature_truncation` reads the
/// mate (Utils/VariationEffect.pm:358). The row therefore carries the region term the
/// local coordinate lands on when applied numerically to the mate transcript, or
/// `intergenic_variant` on a Transcript row when nothing fires. Moving the local
/// position of one record while its mate stays put moves the row's terms; on a
/// chromosome-21 coordinate against a chromosome-5 transcript Ensembl VEP computes a CDS
/// position and a protein position. vep-rs evaluates the mate row at the mate
/// coordinate, the position the row describes.
///
/// THE PREDICATE, on a pair (one tuple per side at a Location, Allele, Feature,
/// Feature_type): the row is a mate-side row of a point breakend
/// (`mate_side_derivation`), Ensembl VEP's set equals the local-coordinate derivation, vep-rs's
/// set equals the mate-coordinate derivation, and the two derivations differ. Both
/// members are excluded. A vep-rs set that differs from the mate-coordinate derivation
/// in any term stays charged, as does an Ensembl VEP set the local read does not explain.
///
/// Returns (excluded_rust, excluded_perl, pairs_by_kind) with the pairs counted under
/// `cross_chromosome` and `same_chromosome`. With no `transcript_models` (no cache)
/// it excludes nothing.
pub fn filter_breakend_mate_local_read_pairs(
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
    transcript_models: &IndexMap<(String, StableId), TranscriptModel>,
) -> (
    BTreeSet<AnnotTuple>,
    BTreeSet<AnnotTuple>,
    IndexMap<String, i64>,
) {
    let mut by_kind: IndexMap<String, i64> = IndexMap::from([
        ("cross_chromosome".to_string(), 0),
        ("same_chromosome".to_string(), 0),
    ]);
    let mut excluded_rust = BTreeSet::new();
    let mut excluded_perl = BTreeSet::new();
    if transcript_models.is_empty() {
        return (excluded_rust, excluded_perl, by_kind);
    }
    // A row by its Location, Allele, Feature and Feature_type.
    type Key<'a> = (&'a str, &'a str, &'a str, &'a str);
    fn key_of(t: &AnnotTuple) -> Key<'_> {
        (&t.0, &t.1, &t.2, &t.3)
    }
    let mut perl_by_key: IndexMap<Key<'_>, Vec<&AnnotTuple>> = IndexMap::new();
    let mut rust_by_key: HashMap<Key<'_>, Vec<&AnnotTuple>> = HashMap::new();
    for t in one_sided(perl_tuples, rust_tuples) {
        perl_by_key.entry(key_of(t)).or_default().push(t);
    }
    for t in one_sided(rust_tuples, perl_tuples) {
        rust_by_key.entry(key_of(t)).or_default().push(t);
    }
    for (key, perl_side) in &perl_by_key {
        let (Some([rt]), [pt]) = (rust_by_key.get(key).map(Vec::as_slice), &perl_side[..]) else {
            continue;
        };
        let (loc, allele, feature, _ftype) = *key;
        let Some(mate) = _bracket_mate(allele) else {
            continue;
        };
        let Some(model) = transcript_models.get(&(mate.0.clone(), StableId::text(feature))) else {
            continue;
        };
        let Some((mate_set, local_set)) = mate_context::mate_side_derivation(loc, allele, model)
        else {
            continue;
        };
        if mate_set == local_set {
            continue;
        }
        let terms = |t: &'_ AnnotTuple| -> BTreeSet<String> {
            t.4.split(',').map(str::to_string).collect()
        };
        let same = |written: BTreeSet<String>, derived: &BTreeSet<&'static str>| {
            written.len() == derived.len()
                && written.iter().all(|term| derived.contains(term.as_str()))
        };
        if !same(terms(pt), &local_set) || !same(terms(rt), &mate_set) {
            continue;
        }
        excluded_perl.insert((*pt).clone());
        excluded_rust.insert((*rt).clone());
        let local_chrom = _chromosome_of_location(loc);
        let kind = if local_chrom.as_deref() == Some(mate.0.as_str()) {
            "same_chromosome"
        } else {
            "cross_chromosome"
        };
        *by_kind.entry(kind.to_string()).or_default() += 1;
    }
    (excluded_rust, excluded_perl, by_kind)
}

pub fn _tuples_by_location(
    tuples: &BTreeSet<AnnotTuple>,
) -> IndexMap<String, BTreeSet<AnnotTuple>> {
    let mut by_loc: IndexMap<String, BTreeSet<AnnotTuple>> = IndexMap::new();
    for t in tuples {
        by_loc.entry(t.0.clone()).or_default().insert(t.clone());
    }
    by_loc
}

/// `_tuples_by_location` over tuples that stay where they are.
fn by_location<'a>(
    tuples: impl Iterator<Item = &'a AnnotTuple>,
) -> IndexMap<&'a str, Vec<&'a AnnotTuple>> {
    let mut by_loc: IndexMap<&str, Vec<&AnnotTuple>> = IndexMap::new();
    for t in tuples {
        by_loc.entry(t.0.as_str()).or_default().push(t);
    }
    by_loc
}

/// The type string Ensembl VEP keeps for a record (Parser/VCF.pm:469-475): INFO/SVTYPE
/// when Perl reads it as true and the ALT string does not open with a VCF 4.4 form, else
/// the ALT string (`/`-joined when a list).
pub fn perl_record_type(sv_type: &str, svtype: Option<&str>) -> String {
    match svtype {
        Some(named) if perl_true(svtype) && !_VCF44_SYMBOLIC.is_match(sv_type) => named,
        _ => sv_type,
    }
    .to_string()
}

/// The Sequence Ontology abbreviation `get_SO_term` derives for a type string, or None
/// when the reference has no term for it (the type is then `vep_skip`).
pub fn perl_so_abbrev(sv_type: &str) -> Option<String> {
    let mut sv_type = sv_type;
    let parts: Vec<String>;
    if sv_type.contains(['/', ',']) {
        let fields: Vec<String> = sv_type.split(['/', ',']).map(str::to_string).collect();
        parts = perl_fields(&fields);
        let normalized: Vec<String> = parts
            .iter()
            .map(|part| {
                let part = pyfmt::strip(part);
                let part = part.strip_prefix('<').unwrap_or(part);
                pyfmt::upper(part.strip_suffix('>').unwrap_or(part))
            })
            .collect();
        if let Some(first) = normalized.first() {
            if normalized.iter().all(|n| n == first) {
                sv_type = &parts[0];
            }
        }
    }
    let abbrev = if let Some(me) = MOBILE_ELEMENT.captures(sv_type) {
        let abbrev = pyfmt::upper(&me[1]);
        let mut subtype = pyfmt::upper(&me[2]);
        // An element the pattern takes as empty is none.
        if let Some(element) = me.get(3).filter(|element| !element.is_empty()) {
            let element = pyfmt::upper(element.as_str());
            let element = if element == "L1" {
                "LINE1".to_string()
            } else {
                element
            };
            if matches!(element.as_str(), "ALU" | "HERV" | "LINE1" | "SVA") {
                subtype = element;
            }
        }
        format!("{abbrev}_{subtype}")
    } else if NAMES_DEL.is_match(sv_type) && NAMES_DUP.is_match(sv_type) {
        "CNV".to_string()
    } else if NAMES_TANDEM_DUP.is_match(sv_type) {
        "TDUP".to_string()
    } else if NAMES_TANDEM_REPEAT.is_match(sv_type) {
        "TREP".to_string()
    } else if NAMES_COPY_NUMBER.is_match(sv_type) {
        if COPY_NUMBER_TWO.is_match(sv_type) {
            "DUP"
        } else if COPY_NUMBER_ZERO.is_match(sv_type) {
            "DEL"
        } else {
            "CNV"
        }
        .to_string()
    } else if NAMES_CNV.is_match(sv_type) {
        "CNV".to_string()
    } else if BREAKEND_FORM.is_match(sv_type) {
        "BND".to_string()
    } else if BRACKETED_FORM.is_match(sv_type) {
        without_subtype(&sv_type.replace(['<', '>'], ""))
    } else {
        sv_type.to_string()
    };
    so_term(&abbrev).map(|_| abbrev)
}

/// The INFO column as the hash Ensembl VEP's parser builds (BaseVCF4.pm:581-583): exact
/// keys, the last occurrence of a repeated key, `KEY=` as the empty value, a bare `KEY` as
/// the empty value too (Perl's undef; callers that test `defined` use `_info_written`), and
/// the text up to a second `=` as the value, as the two-element list assignment
/// `my ($key,$value) = split('=',$info)` keeps it (`END2=30425700=5` is `30425700`).
pub fn _info_fields(info: &str) -> Fields {
    let mut out = Fields::new();
    for item in info.split(';') {
        if !item.is_empty() && item != "." {
            let (key, value) = item.split_once('=').unwrap_or((item, ""));
            let value = value.split('=').next().unwrap_or(value);
            out.insert(key.to_string(), value.to_string());
        }
    }
    out
}

/// Perl's numeric reading of a string, as `abs($svlen)` and `max(@svlen)` apply it:
/// leading whitespace, an optional sign, digits with an optional fraction and exponent;
/// anything else, `.` and the empty string included, is 0 (`300abc` is 300). A reading
/// beyond a double (`1e400`, Perl's `Inf`) is saturated to `_SATURATED_MAGNITUDE` with its
/// sign, so it stays a number the coordinate arithmetic can hold. `NaN`, `Inf` and
/// `Infinity`, Perl's NaN and infinity, read as 0 here, the length vep-rs gives the
/// spellings, so a Location derived from them is the one vep-rs prints; `perl_skip_route`
/// reads the spellings themselves.
///
/// The saturated reading is a whole number and every other one a float.
pub fn perl_numeric_value(text: &str) -> Number {
    let Some(found) = _PERL_NUMBER.find(text) else {
        return Number::Float(0.0);
    };
    let value = pyfmt::float_of(found.as_str()).unwrap_or(0.0);
    if value.is_infinite() {
        let magnitude = Int::from(_SATURATED_MAGNITUDE);
        return Number::Int(if value > 0.0 { magnitude } else { -magnitude });
    }
    Number::Float(value)
}

/// The POS column as Ensembl VEP numifies it (BaseVCF4.pm:263 `$start++` on the structural
/// path): a decimal integer exactly (`4611686018427387904`), any other spelling by
/// `perl_numeric_value` with its fraction dropped (`3e7` is 30000000, `30425949.0` and
/// ` 30425949` are 30425949, `.`, `0x1D04A1D` and `３０４２５９４９` are 0, `-5` is -5).
pub fn perl_pos(text: &str) -> Int {
    let stripped = pyfmt::strip(text);
    let digits = stripped.strip_prefix(['+', '-']).unwrap_or(stripped);
    if !digits.is_empty() && digits.bytes().all(|digit| digit.is_ascii_digit()) {
        return whole_number(stripped);
    }
    perl_numeric_value(text).int()
}

/// Whether Perl reads the text as a number that is not an integer: a fraction (`9.5`,
/// `300.7`), a value past a double (`1e400`), `NaN` or an infinity (`Inf`, `Infinity`), with an
/// optional sign and blanks around; `9`, `9.0`, `3e2` are integers, and `abc`, `.`, the empty
/// string and a fullwidth digit are 0.
pub fn perl_number_is_not_an_integer(text: &str) -> bool {
    let stripped = pyfmt::strip(text);
    if PERL_NAN_WHOLE.is_match(stripped) || PERL_INFINITY_WHOLE.is_match(stripped) {
        return true;
    }
    if !PERL_NUMBER_WHOLE.is_match(stripped) {
        return false;
    }
    pyfmt::float_of(stripped).is_some_and(|value| value.is_infinite() || value.fract() != 0.0)
}

/// The end `get_end` derives for a symbolic record at `pos` (BaseVCF4.pm:297-317) with
/// the breakpoint rule of Parser/VCF.pm:487-489: `pos + |SVLEN|` from the first value as
/// Perl numifies it (`SVLEN=` and `SVLEN=.` are 0, the pair before the start); else END as
/// a number, `None` when it is not one (`validate_vf` skips the record); else the start.
/// A breakpoint reads END only when Perl reads it as false (`0` or empty), which keeps
/// the record on the skip route.
pub fn perl_structural_end(
    pos: Int,
    svlen: Option<&str>,
    end: Option<&str>,
    breakpoint: bool,
) -> Option<Number> {
    let start = pos.saturating_add(1);
    if let Some(svlen) = svlen {
        let first = svlen.split(',').next().unwrap_or(svlen);
        return Some(perl_numeric_value(first).abs().added_to(pos));
    }
    if breakpoint {
        if end.is_none() || perl_true(end) {
            return Some(Number::Int(start));
        }
        return (end == Some("0")).then_some(Number::Int(0));
    }
    let Some(end) = end else {
        return Some(Number::Int(start));
    };
    if !PERL_NUMBER_WHOLE.is_match(pyfmt::strip(end)) {
        return None;
    }
    Some(perl_numeric_value(end))
}

/// Whether Parser/VCF.pm:237-245 routes the record to the sequence path: REF and every ALT
/// together a run of `ACGT`, or no true INFO/SVTYPE and no ALT written in a structural
/// form. A sequence record's Location starts at POS; a structural record's at POS + 1.
pub fn perl_sequence_path(v: &VcfVariant) -> bool {
    let alts = perl_split(&v.alt, ",");
    let svtype = _info_written(&v.info, "SVTYPE");
    let is_acgt = |text: &str| {
        text.bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    };
    let bases = !(v.r#ref.is_empty() && alts.iter().all(String::is_empty))
        && is_acgt(&v.r#ref)
        && alts.iter().all(|alt| is_acgt(alt));
    bases || !(perl_true(svtype.as_deref()) || structural_alt_list(&alts.join(",")))
}

/// The (chromosome, start) the reference prints for a record: POS on the sequence path,
/// POS + 1 on the structural path (BaseVCF4.pm:263).
pub fn reference_start(v: &VcfVariant) -> (String, Int) {
    let start = if perl_sequence_path(v) {
        v.pos
    } else {
        v.pos.saturating_add(1)
    };
    (record_chromosome(v).to_string(), start)
}

/// The route by which Ensembl VEP 116.2 drops a VCF record before its input buffer, or
/// None when it keeps the record: `unsupported_type` (Parser/VCF.pm:477-481), or
/// `gvcf_variant_site` when the unsupported type lists `<NON_REF>` beside another ALT,
/// `incomplete_deletion` (:555-558), `invalid_span` (Parser.pm:475-488: a start or
/// end that is not a number, or a start past the end by more than one), `oversize`
/// (Parser.pm:491-501 over the span `get_end` derives: SVLEN first, else END, else the
/// start; a breakend's END is discarded). An SVLEN or END spelled `Inf` is `oversize`
/// under a cap (Perl's `abs("Inf")`); one spelled `NaN` keeps the record, every comparison
/// on NaN being false. `max_sv_size` is the cap the reference ran with; -1 is no cap, as
/// for the flag. A record the parser routes to the sequence path is kept.
///
/// `max_sv_size` is `_MAX_SV_SIZE` where the caller names none.
pub fn perl_skip_route(v: &VcfVariant, max_sv_size: i64) -> Option<String> {
    if perl_sequence_path(v) {
        return None;
    }
    let alts = perl_split(&v.alt, ",");
    let svtype = _info_written(&v.info, "SVTYPE");
    let Some(abbrev) = perl_so_abbrev(&perl_record_type(&alts.join("/"), svtype.as_deref())) else {
        if is_gvcf_variant_site(&alts, svtype.as_deref()) {
            return Some("gvcf_variant_site".to_string());
        }
        return Some("unsupported_type".to_string());
    };
    let breakend = abbrev == "BND";
    let start = v.pos.saturating_add(1);
    let svlen = _info_written(&v.info, "SVLEN");
    let end_field = _info_written(&v.info, "END");
    let end = perl_structural_end(v.pos, svlen.as_deref(), end_field.as_deref(), breakend);
    // `abs("NaN")` and `abs("Inf")` are the values the spellings name. Every comparison on NaN
    // is false (Parser/VCF.pm:555, Parser.pm:483, 491-494), so the record is kept and printed
    // with the text; Inf is a span over any cap, and under no cap the reference run dies on it.
    let span_text = match &svlen {
        Some(svlen) => svlen.split(',').next().unwrap_or(svlen),
        None if breakend => "",
        None => end_field.as_deref().unwrap_or(""),
    };
    if PERL_NAN_WHOLE.is_match(span_text) {
        return None;
    }
    if PERL_INFINITY_WHOLE.is_match(span_text) {
        return (max_sv_size != -1).then(|| "oversize".to_string());
    }
    // Parser/VCF.pm:555-558 compares `$start > $end` before `validate_vf` runs, numifying
    // an END that is not a number to 0, so a deletion with `END=abc` is incomplete first.
    let compared_end =
        end.unwrap_or_else(|| perl_numeric_value(end_field.as_deref().unwrap_or("")));
    if so_term(&abbrev).is_some_and(|term| term.contains("del"))
        && (compared_end.is_below(start) || (svlen.is_none() && end_field.is_none()))
    {
        return Some("incomplete_deletion".to_string());
    }
    let Some(end) = end else {
        return Some("invalid_span".to_string());
    };
    if end.plus_one().is_below(start) {
        return Some("invalid_span".to_string());
    }
    if max_sv_size != -1 && end.minus(start).is_above(Int::from(max_sv_size)) {
        return Some("oversize".to_string());
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
pub fn vep_rs_unsupported_alt_allele(alt: &str, svtype: Option<&str>) -> String {
    let upper = pyfmt::upper(alt);
    for (prefix, term) in _VEP_RS_ALT_PREFIX_TERMS {
        if upper.starts_with(prefix) {
            return term.to_string();
        }
    }
    if let Some(named) = svtype {
        let named = pyfmt::upper(named);
        if let Some((_, term)) = _VEP_RS_SVTYPE_TERMS.iter().find(|(name, _)| *name == named) {
            return term.to_string();
        }
    }
    perl_record_type(alt, svtype)
}

/// `re.fullmatch(r"\.?[A-Za-z]+\.?", alt)` on an ALT with a dot at either end.
fn is_dot_breakend(alt: &str) -> bool {
    if !(alt.starts_with('.') || alt.ends_with('.')) {
        return false;
    }
    let bases = alt.strip_prefix('.').unwrap_or(alt);
    let bases = bases.strip_suffix('.').unwrap_or(bases);
    !bases.is_empty() && bases.bytes().all(|base| base.is_ascii_alphabetic())
}

/// Every Allele string either engine writes for a VCF record, normalised as the
/// comparator normalises them: the Sequence Ontology class of each symbolic ALT and of the
/// whole ALT list (`copy_number_variation` for a DEL beside a DUP), the raw abbreviation of
/// a type with no term and its bracketed spelling (`CPX`, `<CPX>`, `<NON_REF>`) and the
/// string vep-rs prints for it (`vep_rs_unsupported_alt_allele`: `deletion` for `<del>`,
/// else INFO/SVTYPE or the ALT as written), each bracket breakend string with `chr`
/// stripped inside the brackets, the own-end strings of a breakend record (`N.`, `.N`,
/// `chromosome_breakpoint`, the REF base with the dot), the mate string a `<BND>` record's
/// CHR2/END2 give, and a sequence ALT as written.
pub fn record_alleles(v: &VcfVariant) -> BTreeSet<String> {
    let mut alleles: BTreeSet<String> = BTreeSet::new();
    let info = _info_fields(&v.info);
    let alts = perl_split(&v.alt, ",");
    let svtype = info.get("SVTYPE").map(String::as_str);
    let record_type = perl_record_type(&alts.join("/"), svtype);
    let whole = perl_so_abbrev(&record_type);
    let own_ends = |alleles: &mut BTreeSet<String>| {
        alleles.extend(_OWN_END_ALLELES.iter().map(|allele| allele.to_string()));
        alleles.insert(format!("{}.", v.r#ref));
        alleles.insert(format!(".{}", v.r#ref));
    };
    // Under a type the table lacks, vep-rs keeps each symbolic ALT's own class, INFO/SVTYPE
    // naming the class of a mixed list (a symbolic ALT beside a sequence allele).
    let mixed = alts.iter().any(|alt| alt.starts_with('<'))
        && alts.iter().any(|alt| {
            !alt.starts_with('<')
                && !_BRACKET_MATE_RE.is_match(alt)
                && !matches!(&alt[..], "*" | ".")
        });
    for alt in &alts {
        if alt.starts_with('<') && alt.ends_with('>') {
            match perl_so_abbrev(alt).and_then(|abbrev| so_term(&abbrev)) {
                Some(term) => {
                    alleles.insert(term.to_string());
                }
                None => {
                    let raw = without_subtype(alt.get(1..alt.len() - 1).unwrap_or(""));
                    alleles.insert(format!("<{raw}>"));
                    alleles.insert(raw);
                    alleles.insert(alt.clone());
                }
            }
            if whole.is_none() {
                let own = match svtype {
                    Some(named) if mixed && perl_true(svtype) => format!("<{named}>"),
                    _ => alt.clone(),
                };
                let written = vep_rs_unsupported_alt_allele(&own, svtype);
                if written == own || Some(written.as_str()) == svtype {
                    // An ALT with no class term prints the record's type string.
                    alleles.insert(record_type.clone());
                }
                alleles.insert(written);
            }
        } else if _BRACKET_MATE_RE.is_match(alt) {
            alleles.insert(normalize_allele(alt, "vep-rs"));
            own_ends(&mut alleles);
        } else if is_dot_breakend(alt) {
            alleles.insert(alt.clone());
            alleles.insert("chromosome_breakpoint".to_string());
        } else {
            alleles.insert(alt.clone());
        }
    }
    if let Some(term) = whole.as_deref().and_then(so_term) {
        alleles.insert(term.to_string());
    }
    if alleles.contains("") {
        // The tab format writes an empty Allele as `-`.
        alleles.insert("-".to_string());
    }
    if whole.as_deref() == Some("BND") || svtype == Some("BND") {
        own_ends(&mut alleles);
        if let Some(chr2) = info.get("CHR2").filter(|chr2| !chr2.is_empty()) {
            // The mate string the parser builds from INFO (Parser/VCF.pm:488-497): the mate
            // position is END2 when present, else the record's END, which a breakend record
            // from a caller that writes no END2 carries as its mate position.
            let chr2 = without_chr(chr2);
            for key in ["END2", "POS2", "END"] {
                if let Some(pos2) = info.get(key).filter(|pos2| !pos2.is_empty()) {
                    alleles.insert(format!("N[{chr2}:{pos2}["));
                }
            }
        }
    }
    alleles
}

/// The (chromosome, start, end) both engines write for a structural record: the base
/// after POS, and the end `get_end` derives (SVLEN first, else END, else the start; a
/// breakend's END discarded). A breakend is a record whose type string (`perl_record_type`)
/// has the term `chromosome_breakpoint`: a `<BND>` with no SVTYPE is one, a bracket ALT
/// under `SVTYPE=DEL` is not.
pub fn record_location(v: &VcfVariant) -> (String, Int, Int) {
    let chrom = record_chromosome(v).to_string();
    let start = v.pos.saturating_add(1);
    let svlen = _info_written(&v.info, "SVLEN");
    let end_field = _info_written(&v.info, "END");
    let record_type = perl_record_type(
        &perl_split(&v.alt, ",").join("/"),
        _info_written(&v.info, "SVTYPE").as_deref(),
    );
    let is_bnd = perl_so_abbrev(&record_type).as_deref() == Some("BND");
    let end = perl_structural_end(v.pos, svlen.as_deref(), end_field.as_deref(), is_bnd);
    (chrom, start, end.map_or(start, Number::int))
}

/// The Location both engines write for the rows of a VCF line on the structural path,
/// with the chromosome as the line spells it: the span `record_location` derives from the
/// line's POS and its INFO/SVLEN or INFO/END, low end first. `None` for a line on the
/// sequence path (`perl_sequence_path`), whose rows stand where each allele does.
pub fn structural_line_location(
    chrom: &str,
    pos: &str,
    reference: &str,
    alt: &str,
    info: &str,
) -> Option<String> {
    let line = VcfVariant {
        chrom: chrom.to_string(),
        pos: perl_pos(pos),
        vid: String::new(),
        r#ref: reference.to_string(),
        alt: alt.to_string(),
        variant_type: classify_variant(reference, alt),
        location: format!("{chrom}:{pos}"),
        info: info.to_string(),
        pos_text: Some(pos.to_string()),
    };
    if perl_sequence_path(&line) {
        return None;
    }
    let (_, start, end) = record_location(&line);
    Some(sequence_location(chrom, start, end))
}

/// Whether a reference row's span lies on `chrom` within one base of `start`..`end`:
/// the reference wrote the record (a `<CNV:TR>` expanded to a literal allele sits at the
/// run's 3' end, one base past it).
pub fn _reference_row_beside(
    spans: &[(String, Int, Int)],
    chrom: &str,
    start: Int,
    end: Int,
) -> bool {
    spans
        .iter()
        .any(|s| s.0 == chrom && s.1 <= end.saturating_add(1) && s.2 >= start.saturating_sub(1))
}

/// The (Location, Allele) pairs vep-rs prints for a record: on the structural path the
/// span `record_location` derives, printed low end first (`21:30425650-30425651` for a
/// deletion whose end is before its start), under every allele of `record_alleles`; on the
/// sequence path each real allele at the place the parser gives it (`gvcf_site_rows`).
pub fn vep_rs_record_rows(v: &VcfVariant) -> BTreeSet<(String, String)> {
    let chrom = record_chromosome(v);
    if perl_sequence_path(v) {
        return site_rows(v)
            .into_iter()
            .map(|(allele, start, end)| {
                (
                    normalize_location(&sequence_location(chrom, start, end)),
                    allele,
                )
            })
            .collect();
    }
    let (_chrom, start, end) = record_location(v);
    let location = normalize_location(&sequence_location(chrom, start, end));
    record_alleles(v)
        .into_iter()
        .map(|allele| (location.clone(), allele))
        .collect()
}

/// A record at its start: its place in the input, its skip route, the alleles it produces
/// and the far end of its span.
struct Placed {
    record: usize,
    route: Option<String>,
    alleles: BTreeSet<String>,
    end: Int,
}

/// The tuples of `tuples` by the (chromosome, start) of their Location, each with the
/// Location's end.
fn located(tuples: &BTreeSet<AnnotTuple>) -> HashMap<(String, Int), Vec<(&AnnotTuple, Int)>> {
    let mut out: HashMap<(String, Int), Vec<(&AnnotTuple, Int)>> = HashMap::new();
    for t in tuples {
        if let Some((chrom, start, end)) = _location_span(&t.0) {
            out.entry((chrom, start)).or_default().push((t, end));
        }
    }
    out
}

/// The tuples of a start's bucket under one of `alleles`, at `end` alone where the start
/// is contested.
fn rows_of<'a>(
    bucket: Option<&Vec<(&'a AnnotTuple, Int)>>,
    alleles: &BTreeSet<String>,
    end: Int,
    contested: bool,
) -> Vec<&'a AnnotTuple> {
    bucket
        .into_iter()
        .flatten()
        .filter(|(t, t_end)| alleles.contains(&t.1) && (!contested || *t_end == end))
        .map(|(t, _)| *t)
        .collect()
}

/// The records the 116 reference dropped as `vep_skip`, their vep-rs rows, every
/// disagreement between the records and the outputs, and the skip-route residue of a
/// scored engine other than vep-rs.
///
/// THE SHAPE. 116.2 writes at least one row for every record it keeps and no row in any
/// format for a record it marks `vep_skip`; vep-rs annotates all of them. A record is
/// named by the input: `perl_skip_route` over its ALT and INFO. It is CONFIRMED skipped
/// when the reference wrote no row under the record's own Location start and an allele the
/// record produces (`record_alleles`); a kept record sharing the start, or the whole
/// Location string, with a skipped one therefore never hides it and is never mistaken for
/// it. When two records at one start produce the same allele, the end of the Location
/// tells them apart. A record on the sequence path is placed at POS under its sequence
/// alleles (`gvcf_site_rows`), whatever INFO it carries: an SNV with SVTYPE and END is not
/// at the structural span those would give.
///
/// SET ASIDE: every vep-rs tuple at a confirmed record's start under one of its alleles
/// (the same end rule applies). The start is the one vep-rs prints, the lower end of the
/// span: an incomplete deletion (`END=POS`, `SVLEN=0`, `SVLEN=abc`) is the pair
/// `POS-(POS+1)`, keyed at POS. `<CNV:TR>` records are outside the class on both sides:
/// the reference keeps them and writes them at another Location
/// (`filter_cnv_tr_expansion_divergences`).
///
/// PROBLEMS, one line each, make the rule's blind spots loud: a skip-route record whose
/// own rows the reference wrote (a reference run under --dont_skip or another
/// --max_sv_size, declared with --reference-max-sv-size); a vep-rs structural row at a
/// Location the reference wrote nothing for that no skip-route record explains and no
/// other class pairs with a reference row elsewhere (the dot-breakend pair, the record
/// whose reference Location carries a non-integer coordinate), so a record dropped on a
/// route 116.2 shares with 115.2 (`start > end+1`) is reported; two records at one start
/// whose alleles and Locations coincide, one skipped and one kept.
///
/// That residue is a problem for vep-rs alone, whose rows the reference's drop rules
/// account for. Under any other `scored_engine` the engine's structural rows at such
/// Locations are its own disagreements with the reference, returned as the fourth member
/// for the report to count; they stay in both figures and stop nothing.
///
/// The SNP/indel comparator opens the input only to recognise gVCF variant sites and
/// judges this class by the output-only rule (`ReferenceRecordScope`): a Location with
/// vep-rs structural rows and no reference row.
///
/// The gVCF variant site (`gvcf_variant_site`) is its own class,
/// `GVCF_VARIANT_SITE_DROPPED_BUCKET`: its record is confirmed, so it explains its Location
/// in the residue check, and is neither counted nor set aside here.
///
/// Returns (excluded_rust, records, problems, scored_engine_residue).
///
/// `max_sv_size` is `_MAX_SV_SIZE` and `scored_engine` is `vep-rs` where the caller names none.
pub fn skipped_record_rows(
    input_variants: &[VcfVariant],
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
    max_sv_size: i64,
    scored_engine: &str,
) -> (BTreeSet<AnnotTuple>, i64, Vec<String>, BTreeSet<AnnotTuple>) {
    let perl_at = located(perl_tuples);
    let rust_at = located(rust_tuples);
    let mut records_at: IndexMap<(String, Int), Vec<Placed>> = IndexMap::new();
    for (record, v) in input_variants.iter().enumerate() {
        let chrom = record_chromosome(v).to_string();
        let (alleles, lo, hi) = if perl_sequence_path(v) {
            let placed = site_rows(v);
            let bounds = placed
                .iter()
                .map(|(_allele, start, end)| (*start.min(end), *start.max(end)))
                .reduce(|(lo, hi), (start, end)| (lo.min(start), hi.max(end)));
            let Some((lo, hi)) = bounds else {
                continue;
            };
            let alleles = placed.into_iter().map(|(allele, _, _)| allele).collect();
            (alleles, lo, hi)
        } else {
            let (_chrom, start, end) = record_location(v);
            (record_alleles(v), start.min(end), start.max(end))
        };
        records_at.entry((chrom, lo)).or_default().push(Placed {
            record,
            route: perl_skip_route(v, max_sv_size),
            alleles,
            end: hi,
        });
    }

    let mut excluded: BTreeSet<AnnotTuple> = BTreeSet::new();
    let mut problems: Vec<String> = Vec::new();
    let mut confirmed_starts: BTreeSet<(String, Int)> = BTreeSet::new();
    let mut confirmed = 0;
    for (key, records) in &records_at {
        for placed in records {
            let Some(route) = &placed.route else {
                continue;
            };
            let v = &input_variants[placed.record];
            // A kept record at this start producing one of this record's alleles: the end
            // decides whose rows are whose, and equal ends cannot be told apart.
            let rivals: Vec<&Placed> = records
                .iter()
                .filter(|rival| {
                    rival.record != placed.record
                        && rival.route.is_none()
                        && !rival.alleles.is_disjoint(&placed.alleles)
                })
                .collect();
            let contested = !rivals.is_empty();
            if let Some(twin) = rivals.iter().find(|rival| rival.end == placed.end) {
                let twin = &input_variants[twin.record];
                problems.push(format!(
                    "{} {}>{} and {} {}>{} share {}:{}-{} and an allele; \
                     the skipped one cannot be told from the kept one",
                    v.vid, v.r#ref, v.alt, twin.vid, twin.r#ref, twin.alt, key.0, key.1, placed.end
                ));
                continue;
            }
            let kept = rows_of(perl_at.get(key), &placed.alleles, placed.end, contested);
            if !kept.is_empty() {
                let written: BTreeSet<&str> = kept.iter().map(|t| t.1.as_str()).collect();
                problems.push(format!(
                    "{} {}>{} at {}:{} is on the {route} route, yet the reference wrote \
                     {} row(s) under its own allele(s) ({}): a reference run \
                     under --dont_skip or another --max_sv_size (pass --reference-max-sv-size)",
                    v.vid,
                    v.r#ref,
                    v.alt,
                    key.0,
                    key.1,
                    kept.len(),
                    written.into_iter().collect::<Vec<_>>().join(", ")
                ));
                continue;
            }
            confirmed_starts.insert(key.clone());
            if route == "gvcf_variant_site" {
                continue;
            }
            confirmed += 1;
            excluded.extend(
                rows_of(rust_at.get(key), &placed.alleles, placed.end, contested)
                    .into_iter()
                    .cloned(),
            );
        }
    }
    // The output-only rule's residue: a Location with vep-rs structural rows and no
    // reference row that no confirmed record explains and no pair class reaches.
    let mut paired = filter_dot_breakend_pairs(perl_tuples, rust_tuples).0;
    paired.extend(
        filter_non_integer_end_rows(perl_tuples, rust_tuples, input_variants, max_sv_size).1,
    );
    let mut explained = confirmed_starts;
    explained.extend(
        paired
            .iter()
            .filter_map(|t| _location_span(&t.0))
            .map(|(chrom, start, _end)| (chrom, start)),
    );
    let perl_locations: BTreeSet<&str> = perl_tuples.iter().map(|t| t.0.as_str()).collect();
    let mut residue: BTreeSet<AnnotTuple> = BTreeSet::new();
    for t in rust_tuples {
        let span = _location_span(&t.0);
        if perl_locations.contains(t.0.as_str())
            || !is_structural_allele(&t.1)
            || t.1 == _CNV_TR_RUST_ALLELE
        {
            continue;
        }
        let Some((chrom, start, _end)) = span else {
            continue;
        };
        if explained.contains(&(chrom, start)) {
            continue;
        }
        if scored_engine == "vep-rs" {
            problems.push(format!(
                "{} {}: the reference wrote nothing at this Location and no input record on \
                 a skip route starts there",
                t.0, t.1
            ));
        } else {
            residue.insert(t.clone());
        }
    }
    problems.sort();
    problems.dedup();
    (excluded, confirmed, problems, residue)
}

/// The mask over `skipped_record_rows` without its problems (`apply_masks` calls
/// `skipped_record_rows` itself, to stop on them): the vep-rs tuples of every record the
/// reference dropped as `vep_skip` or refused in `validate_vf`, the gVCF variant site
/// excepted, and the number of such records.
///
/// `max_sv_size` is `_MAX_SV_SIZE` where the caller names none.
pub fn filter_reference_skipped_records(
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
    input_variants: &[VcfVariant],
    max_sv_size: i64,
) -> (BTreeSet<AnnotTuple>, i64) {
    let (excluded, records, _problems, _residue) = skipped_record_rows(
        input_variants,
        perl_tuples,
        rust_tuples,
        max_sv_size,
        "vep-rs",
    );
    (excluded, records)
}

/// A dot-notation breakend without SVTYPE: the reference's two-coordinate tuple and
/// vep-rs's point under the same breakend allele on the same feature, set aside as a pair
/// (`DOT_BREAKEND_PAIR_BUCKET`). The point `c:(P+1)` is the span's start or the base after
/// it: the end of `c:P-(P+1)` for a one-base REF, the start of `c:(P+1)-(P+len(REF)-1)` for
/// a longer REF beside an ALT of another length, the second base of `c:P-(P+len(REF)-1)`
/// beside an ALT of the REF's length (`dot_breakend_points`). The allele is
/// `chromosome_breakpoint` for a bare `N.`, and the `N.` and `N[21:30425700[` strings of one
/// whose INFO names a mate through CHR2 and END2 (Parser/VCF.pm:491-497). Returns
/// (excluded_rust, excluded_perl).
///
/// Where two tuples of the reference can take one point, the first in the tuples' order
/// takes it.
pub fn filter_dot_breakend_pairs(
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
) -> (BTreeSet<AnnotTuple>, BTreeSet<AnnotTuple>) {
    // A point by its Location, Allele, Feature and Feature_type.
    type Point<'a> = (String, &'a str, &'a str, &'a str);
    let mut rust_points: HashMap<Point<'_>, Vec<&AnnotTuple>> = HashMap::new();
    for t in one_sided(rust_tuples, perl_tuples) {
        let span = _location_span(&t.0);
        if is_breakend_allele(&t.1) && span.is_some_and(|(_, start, end)| start == end) {
            rust_points
                .entry((t.0.clone(), &t.1, &t.2, &t.3))
                .or_default()
                .push(t);
        }
    }
    let mut excluded_rust: BTreeSet<AnnotTuple> = BTreeSet::new();
    let mut excluded_perl: BTreeSet<AnnotTuple> = BTreeSet::new();
    for pt in one_sided(perl_tuples, rust_tuples) {
        let span = _location_span(&pt.0);
        let breakend = is_breakend_allele(&pt.1);
        let Some(span) = span.filter(|(_, start, end)| breakend && start < end) else {
            continue;
        };
        for point in dot_breakend_points(&span) {
            let partners =
                rust_points.remove(&(point, pt.1.as_str(), pt.2.as_str(), pt.3.as_str()));
            if let Some(partners) = partners.filter(|partners| !partners.is_empty()) {
                excluded_perl.insert(pt.clone());
                excluded_rust.extend(partners.into_iter().cloned());
                break;
            }
        }
    }
    (excluded_rust, excluded_perl)
}

/// The reference's tuples whose Location carries a coordinate that is not an integer
/// (`21:30425651-30425949.5`, `21:30425651-NaN`, `21:30425949.0`), an INFO/END, SVLEN or
/// POS taken as written, and vep-rs's tuples of the same input records, annotated at the
/// numified value (`SVLEN=300.7` ends at 30425950, `SVLEN=NaN` is 0 and the pair
/// `21:30425650-30425651`, POS `30425949.0` on `A G` is `21:30425949`)
/// (`NON_INTEGER_END_BUCKET`). A record is a reference row's when the row's start, read as
/// its number, is the start the reference prints for the record (POS + 1 on the structural
/// path, POS on the sequence path); the record's vep-rs tuples are the one-sided ones at
/// the Locations and under the alleles vep-rs prints for it (`vep_rs_record_rows`). With a
/// FASTA the reference expands a `<CNV:TR>` whose SVLEN or END is a number that is not an
/// integer (`SVLEN=9.5`, `SVLEN=NaN`) to a literal insertion whose coordinates fail
/// `start = end + 1` (Parser.pm:564-571) and writes nothing for the record; vep-rs's rows at
/// the numified run (`21:30425588-30425596`, the pair `21:30425587-30425588` for NaN) are
/// set aside when the record is on no skip route under `max_sv_size` and the reference
/// wrote nothing from the base before the span to the base after it.
/// Returns (excluded_perl, excluded_rust).
///
/// `max_sv_size` is `_MAX_SV_SIZE` where the caller names none.
pub fn filter_non_integer_end_rows(
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
    input_variants: &[VcfVariant],
    max_sv_size: i64,
) -> (BTreeSet<AnnotTuple>, BTreeSet<AnnotTuple>) {
    let excluded_perl: BTreeSet<AnnotTuple> = perl_tuples
        .iter()
        .filter(|t| _non_integer_end(&t.0))
        .cloned()
        .collect();
    let mut starts: BTreeSet<(String, Int)> = BTreeSet::new();
    for t in &excluded_perl {
        let (chrom, coords) = t.0.split_once(':').unwrap_or((&t.0, ""));
        let first = coords.split('-').next().unwrap_or(coords);
        starts.insert((chrom.to_string(), perl_numeric_value(first).int()));
    }
    let mut excluded_rust: BTreeSet<AnnotTuple> = BTreeSet::new();
    let rust_only: Vec<&AnnotTuple> = one_sided(rust_tuples, perl_tuples).collect();
    let mut set_aside = |v: &VcfVariant| {
        let rows = vep_rs_record_rows(v);
        excluded_rust.extend(
            rust_only
                .iter()
                .filter(|t| {
                    rows.iter()
                        .any(|(location, allele)| *location == t.0 && *allele == t.1)
                })
                .map(|t| (*t).clone()),
        );
    };
    if !starts.is_empty() {
        for v in input_variants {
            if starts.contains(&reference_start(v)) {
                set_aside(v);
            }
        }
    }
    let perl_spans: Vec<(String, Int, Int)> = perl_tuples
        .iter()
        .filter_map(|t| _location_span(&t.0))
        .collect();
    for v in input_variants {
        let alts = perl_split(&pyfmt::upper(&v.alt), ",");
        if alts.is_empty()
            || alts.iter().any(|alt| alt != "<CNV:TR>")
            || perl_skip_route(v, max_sv_size).is_some()
        {
            continue;
        }
        let svlen = _info_written(&v.info, "SVLEN");
        let end_field = _info_written(&v.info, "END");
        let span_text = match &svlen {
            Some(svlen) => Some(svlen.split(',').next().unwrap_or(svlen)),
            None => end_field.as_deref(),
        };
        if !span_text.is_some_and(perl_number_is_not_an_integer) {
            continue;
        }
        let (chrom, start, end) = record_location(v);
        if _reference_row_beside(&perl_spans, &chrom, start.min(end), start.max(end)) {
            continue;
        }
        set_aside(v);
    }
    (excluded_perl, excluded_rust)
}

pub fn _breakpoint_position(loc: &str) -> Option<Int> {
    _location_span(loc).map(|(_chrom, start, _end)| start)
}

/// A bracket breakend whose own-end rows the reference lost to its mate chromosome.
///
/// THE SHAPE. A record's slice is named by the first transcript that reaches it
/// (AnnotationType/Transcript.pm 116.2:122, `$vf->{slice} ||= $tr->{slice}`). A buffer
/// holds the records of one chromosome (InputBuffer.pm 116.2:191-207, 223-234), its
/// regions load in the order its records first ask for them, each record's own region
/// before its mates' and a region listed once at its first request (AnnotationSource.pm
/// 116.2:235-258), and the transcripts are visited in that order. When an earlier record
/// of the buffer asked for a region holding a transcript within reach of the record's
/// mate before the record's own region was asked for, that mate-chromosome transcript
/// names the slice, every own-chromosome transcript then fails the slice test (:125-143:
/// neither the record's slice nor any breakend's is on the transcript's chromosome) and
/// the reference writes nothing under the own-end allele while keeping the mate rows.
/// The result depends on the order of the input, which no annotation should; the record
/// alone, with its mate, or under --buffer_size 1 keeps the rows, and --fork does not
/// change the result.
///
/// THE PREDICATE, on a Location whose bracket allele names a mate on another
/// chromosome: the reference wrote a Transcript tuple under the bracket allele and no
/// tuple under any own-end allele; then every vep-rs-only Transcript tuple under an
/// own-end allele whose transcript the cache files on the Location's chromosome within
/// `_BREAKEND_FEATURE_DISTANCE` of the breakpoint, the reference's own admission rule
/// (StructuralVariationOverlap.pm 116:144-160), is excluded. A transcript farther away
/// than that, or on another chromosome, stays charged. With no cache it excludes
/// nothing.
pub fn filter_bnd_own_rows_lost(
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
    transcripts_by_chr: &IndexMap<String, BTreeSet<StableId>>,
    transcript_spans: &IndexMap<(String, StableId), TranscriptSpan>,
) -> BTreeSet<AnnotTuple> {
    let mut excluded: BTreeSet<AnnotTuple> = BTreeSet::new();
    if transcripts_by_chr.is_empty() || transcript_spans.is_empty() {
        return excluded;
    }
    let perl_by_loc = by_location(perl_tuples.iter());
    let no_tuples = Vec::new();
    for (loc, rust_here) in &by_location(one_sided(rust_tuples, perl_tuples)) {
        let own = _chromosome_of_location(loc);
        let pos = _breakpoint_position(loc);
        let (Some(own), Some(pos)) = (own, pos) else {
            continue;
        };
        let perl_here = perl_by_loc.get(loc).unwrap_or(&no_tuples);
        let mut other_chromosome = false;
        for t in perl_here.iter().chain(rust_here) {
            if _bracket_mate(&t.1).is_some_and(|(mate, _)| mate != own) {
                other_chromosome = true;
            }
        }
        if !other_chromosome {
            continue;
        }
        if perl_here.iter().any(|t| is_own_breakend_allele(&t.1)) {
            continue;
        }
        if !perl_here
            .iter()
            .any(|t| _bracket_mate(&t.1).is_some() && t.2 != "-")
        {
            continue;
        }
        let local = transcripts_by_chr.get(&own);
        for t in rust_here {
            if !is_own_breakend_allele(&t.1) || t.2 == "-" {
                continue;
            }
            let id = StableId::text(&t.2);
            if !local.is_some_and(|local| local.contains(&id)) {
                continue;
            }
            let Some(span) = transcript_spans.get(&(own.clone(), id)) else {
                continue;
            };
            let reach = Int::from(_BREAKEND_FEATURE_DISTANCE);
            if minus(span.1, reach) <= pos && pos <= plus(span.2, reach) {
                excluded.insert((*t).clone());
            }
        }
    }
    excluded
}

/// A bracket breakend's mate rows lost to a mate chromosome spelled in a way the reference
/// does not resolve (`tools/vep-sv-compare`'s `filter_bnd_synonym_mate_rows`).
///
/// THE SHAPE. The reference strips a lower-case `chr` from the mate's name and keeps any other
/// spelling (StructuralVariationOverlap.pm 116:86-103), so `[CHR17:7675000[` or `Chr21` names
/// no chromosome of the cache; `_close_to_feature`'s name test (Utils/VariationEffect.pm
/// 116:111-124) then refuses every transcript of the mate, and the reference writes no mate
/// row. Where the record's own end reaches no transcript either, it writes no CSQ at all. The
/// candidate resolves the name through the synonym table and writes the mate rows.
///
/// THE PREDICATE: a candidate-only Transcript tuple under a bracket allele the reference wrote
/// no tuple under at the Location, whose mate name as the reference reads it is no chromosome
/// of the cache. Where the reference wrote a tuple at the Location, the transcript is on the
/// chromosome the candidate resolved the mate to or, for a name the cache files nowhere, not
/// on the Location's own chromosome. Where it wrote none, the transcript is on the mate's
/// chromosome within `_BREAKEND_FEATURE_DISTANCE` of the mate (StructuralVariationOverlap.pm
/// 116:144-160). With no cache it excludes nothing.
pub fn filter_bnd_synonym_mate_rows(
    perl_tuples: &BTreeSet<AnnotTuple>,
    rust_tuples: &BTreeSet<AnnotTuple>,
    transcripts_by_chr: &IndexMap<String, BTreeSet<StableId>>,
    transcript_spans: &IndexMap<(String, StableId), TranscriptSpan>,
) -> BTreeSet<AnnotTuple> {
    let mut excluded: BTreeSet<AnnotTuple> = BTreeSet::new();
    if transcripts_by_chr.is_empty() {
        return excluded;
    }
    let perl_by_loc = by_location(perl_tuples.iter());
    let no_tuples = Vec::new();
    for (loc, rust_here) in &by_location(one_sided(rust_tuples, perl_tuples)) {
        let Some(own) = _chromosome_of_location(loc) else {
            continue;
        };
        let perl_here = perl_by_loc.get(loc).unwrap_or(&no_tuples);
        for t in rust_here {
            let Some((mate, mate_pos)) = _bracket_mate(&t.1) else {
                continue;
            };
            let resolved = _BRACKET_MATE_RE
                .captures(&t.1)
                .map(|found| {
                    found[1]
                        .strip_prefix("chr")
                        .unwrap_or(&found[1])
                        .to_string()
                })
                .is_some_and(|name| transcripts_by_chr.contains_key(&name));
            if resolved || t.2 == "-" || perl_here.iter().any(|p| p.1 == t.1) {
                continue;
            }
            let id = StableId::text(&t.2);
            let on = |chrom: &str| {
                transcripts_by_chr
                    .get(chrom)
                    .is_some_and(|ids| ids.contains(&id))
            };
            let lost = if perl_here.is_empty() {
                let reach = Int::from(_BREAKEND_FEATURE_DISTANCE);
                transcript_spans
                    .get(&(mate.clone(), id.clone()))
                    .is_some_and(|span| {
                        minus(span.1, reach) <= mate_pos && mate_pos <= plus(span.2, reach)
                    })
            } else if transcripts_by_chr.contains_key(&mate) {
                on(&mate)
            } else {
                !on(&own)
            };
            if lost {
                excluded.insert((*t).clone());
            }
        }
    }
    excluded
}

/// The structural-variant comparator's tuple: Location, Allele, Feature, Feature_type and the
/// Consequence set.
pub type StructuralTuple = (String, String, String, String, String);

/// The structural-variant comparator (`tools/vep-sv-compare`): its
/// predicates decide the classes of `STRUCTURAL_VARIANT_CLASSES` here, so the rows set aside
/// are the rows its adjusted figures set aside.
///
/// The predicates are this module's own functions, so there is nothing to import.
pub fn structural_variant_comparator() {}

/// What the built-in structural-variant classes read beside the entries: the input records
/// (the comparator's `VcfVariant`, whose ALT and INFO the skip routes read), the candidate
/// cache's three views (the transcript ids per chromosome, the spans and the models), and each
/// consequence entry's tuple as the structural-variant comparator keys a row. A text entry's
/// Location is its own; a VCF or JSON entry's is the one vep-rs prints for its input record
/// under the entry's allele (`vep_rs_record_rows`), so an entry of an output line no input
/// record matches has no tuple and is outside every class. `classes` is the subset of
/// `STRUCTURAL_VARIANT_CLASSES` requested; the skipped-record class needs the input VCF, the
/// two breakend classes the candidate JSON cache.
#[derive(Debug)]
pub struct StructuralVariants {
    pub classes: BTreeSet<String>,
    pub csc: (),
    pub input_variants: Vec<VcfVariant>,
    pub by_chr: IndexMap<String, BTreeSet<StableId>>,
    pub spans: IndexMap<(String, StableId), TranscriptSpan>,
    pub models: IndexMap<(String, StableId), TranscriptModel>,
    /// A record's Locations by allele and its one Location, read at the first entry that
    /// asks for them.
    pub _record_locations: RefCell<IndexMap<i64, (Fields, Option<String>)>>,
}

impl StructuralVariants {
    pub fn new(classes: &[&str], input_path: Option<&Path>, cache: Option<&Path>) -> Result<Self> {
        let tables = classes::tables();
        let requested: BTreeSet<String> = tables
            .structural_variant_classes()
            .iter()
            .filter(|class| classes.contains(class))
            .map(|class| class.to_string())
            .collect();
        if requested.contains(REFERENCE_SKIPPED_RECORD) && input_path.is_none() {
            return Err(Error::comparator(format!(
                "--masks {REFERENCE_SKIPPED_RECORD} needs --input, the input VCF"
            )));
        }
        let cache_reading: Vec<&str> = tables
            .cache_reading_classes()
            .iter()
            .copied()
            .filter(|class| requested.contains(*class))
            .collect();
        if !cache_reading.is_empty() && cache.is_none() {
            return Err(Error::comparator(format!(
                "--masks {} needs --json-cache, the candidate JSON cache",
                cache_reading.join(",")
            )));
        }
        structural_variant_comparator();
        let mut input_variants = Vec::new();
        if let Some(input_path) = input_path {
            input_variants = parse_vcf(input_path).map_err(|error| {
                Error::comparator(format!(
                    "cannot read {}: {}",
                    path_text(input_path),
                    strerror(&error)
                ))
            })?;
        }
        let mut by_chr = IndexMap::new();
        let mut spans = IndexMap::new();
        let mut models = IndexMap::new();
        if let Some(cache) = cache {
            let transcripts = cache::joined(cache, "transcripts");
            if !std::fs::metadata(&transcripts).is_ok_and(|found| found.is_dir()) {
                return Err(Error::comparator(format!(
                    "--json-cache: {} has no transcripts/ directory",
                    path_text(cache)
                )));
            }
            (by_chr, spans, models) = load_cache_authority(cache);
        }
        Ok(StructuralVariants {
            classes: requested,
            csc: (),
            input_variants,
            by_chr,
            spans,
            models,
            _record_locations: RefCell::new(IndexMap::new()),
        })
    }

    /// The entry's Location as the structural-variant comparator spells it: the text
    /// formats' own, normalised; a VCF or JSON entry's from its input record, the Location of
    /// the record's row under the entry's allele, else the record's one Location, else none.
    pub fn location(&self, e: &Entry) -> Option<String> {
        if let Some(location) = &e.location {
            return Some(normalize_location(location));
        }
        let record = e.record?;
        let records = self.input_variants.len() as i64;
        if record >= records {
            return None;
        }
        // A record counted from the end names one, as an index below zero does.
        let index = if record < 0 { record + records } else { record };
        if index < 0 {
            crate::uncaught!("list index out of range");
        }
        let mut read = self._record_locations.borrow_mut();
        let (by_allele, only) = read.entry(record).or_insert_with(|| {
            let rows = vep_rs_record_rows(&self.input_variants[index as usize]);
            let mut by_allele = Fields::new();
            let mut locations: BTreeSet<&str> = BTreeSet::new();
            for (location, allele) in &rows {
                by_allele.insert(allele.clone(), location.clone());
                locations.insert(location);
            }
            let only = match locations.len() {
                1 => locations.first().map(|location| location.to_string()),
                _ => None,
            };
            (by_allele, only)
        });
        by_allele
            .get(&normalize_allele(&dash(&e.allele), "vep-rs"))
            .or(only.as_ref())
            .cloned()
    }

    /// Each consequence entry's tuple by its key; a record entry has none.
    pub fn tuples(&self, entries: &[Entry]) -> IndexMap<String, StructuralTuple> {
        let mut out: IndexMap<String, StructuralTuple> = IndexMap::new();
        for e in entries {
            if e.key.ends_with("|record") {
                continue;
            }
            let Some(location) = self.location(e) else {
                continue;
            };
            let last = e.key.rsplit('|').next().unwrap_or(&e.key);
            let feature_type = last.split('#').next().unwrap_or(last);
            let allele = normalize_allele(&dash(&e.allele), "vep-rs");
            out.insert(
                e.key.clone(),
                (
                    location,
                    allele,
                    dash(&e.feature),
                    feature_type.to_string(),
                    consequence_set(&e.consequence),
                ),
            );
        }
        out
    }
}
