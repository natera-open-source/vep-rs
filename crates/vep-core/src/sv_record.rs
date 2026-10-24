// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The alleles a structural VCF record writes rows for, built as Ensembl VEP builds the
//! one StructuralVariationFeature that carries the record (`Parser/VCF.pm`
//! `create_StructuralVariationFeatures`, lines 466-575) and its breakends
//! (ensembl-variation `StructuralVariationFeature.pm` `_parse_breakends`, lines
//! 1262-1323), shared by every VCF reader: the record's type and class, each allele's
//! span, a breakpoint record's own breakend and mates, a tandem repeat's expanded
//! allele string, and the allele numbering Ensembl VEP prints.

use std::fmt;
use std::sync::Arc;

use crate::coordinate::MAX_COORD;
use crate::sv_routing::{
    joined, per_alt_classes, perl_split, perl_true, sv_so_term_and_class, sv_type_abbreviation,
    sv_type_string, variant_alts,
};
use crate::sv_span::{
    deletion_looks_incomplete, perl_numeric_f64, structural_end, structural_start, svlen_values,
    SpanClass, SpanError,
};
use crate::variant::{InputVariant, VariantClass};
use crate::vcf_info::InfoFields;

/// Ensembl VEP's default `--max_sv_size` (`Config.pm:315`), in bases.
pub const DEFAULT_MAX_SV_SIZE: u64 = 10_000_000;

/// Why a structural record writes no allele.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StructuralRecordError {
    /// Ensembl VEP skips the record on its span (`validate_vf`).
    Span(SpanError),
    /// The record's type is in no `%SO_TERMS` table and no ALT is symbolic, so there is
    /// no allele to write even under the `vep_skip` mark; the text is the type as
    /// `get_SO_term` derives it.
    UnsupportedType(String),
}

impl fmt::Display for StructuralRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructuralRecordError::Span(error) => error.fmt(f),
            StructuralRecordError::UnsupportedType(abbrev) => {
                write!(f, "{abbrev} is not a supported structural variant type")
            }
        }
    }
}

impl std::error::Error for StructuralRecordError {}

impl From<SpanError> for StructuralRecordError {
    fn from(error: SpanError) -> Self {
        StructuralRecordError::Span(error)
    }
}

/// The variants of a structural record at VCF position `pos` on `chr`: `written` are
/// the record's ALTs with their written positions (`sv_routing::written_alts`), `info`
/// the record's INFO column and `id` the ID column when it names the record. One
/// variant per ALT that names an allele under the record's class after the gVCF rule
/// (`sv_routing::variant_alts`), or under its own class with the `vep_skip` mark for a
/// type `%SO_TERMS` lacks; a breakpoint record with parsed mates writes its own
/// breakend and one variant per mate instead, at its first ALT. The record's allele
/// strings keep every written ALT. The caller sets what the record's line adds (`id`,
/// `raw_input`, `input_line_name`).
pub fn structural_variants(
    chr: &str,
    pos: u64,
    ref_allele: &str,
    written: &[(usize, &str)],
    info: &InfoFields<'_>,
    id: Option<&str>,
) -> Result<Vec<InputVariant>, StructuralRecordError> {
    let svtype = info.value("SVTYPE").filter(|v| perl_true(v));
    let mut alts = written.to_vec();
    variant_alts(&mut alts, svtype);
    let mut variants = Vec::with_capacity(alts.len());
    // One structural variant record has one class and one SO term, read from the
    // joined ALTs with INFO/SVTYPE standing in for a string outside the VCF 4.4 forms
    // (Parser.pm `get_SO_term`, Parser/VCF.pm:472-475), so `N.` with `SVTYPE=DEL`
    // is a deletion; the term labels the row of every ALT, each ALT keeping its own
    // bytes (`per_alt_classes`).
    let joined_alts = joined(&alts, "/");
    let record_term = sv_so_term_and_class(&sv_type_abbreviation(&joined_alts, svtype));
    let record_class = record_term.map(|(_, class)| class);
    // The deletion test of Parser/VCF.pm:555 runs on the term, or on the raw type
    // string when `%SO_TERMS` has no entry for it (`:477-481`).
    let deletion_test_term = record_term
        .map(|(term, _)| term.to_string())
        .unwrap_or_else(|| sv_type_string(&joined_alts, svtype));
    // VEP's allele string for a structural variant record is every ALT of the
    // record joined by `/`, `*`, a `<NON_REF>` block and the sequence alleles of a
    // mixed record included (Parser/VCF.pm:469); a breakend record's carries the mate
    // INFO/CHR2 names and, as the record's string, the reference base in front
    // (`:491-499`).
    let breakpoint = (record_class == Some(VariantClass::Translocation))
        .then(|| breakend_record(ref_allele, written, info));
    let uploaded: Arc<str> = breakpoint
        .as_ref()
        .map_or_else(|| joined(written, "/"), |record| record.alts.clone())
        .into();
    let record_string: Option<Arc<str>> = breakpoint.as_ref().map(|record| {
        if record.reference_prefixed {
            format!("{ref_allele}/{}", record.alts).into()
        } else {
            record.alts.as_str().into()
        }
    });
    // Every ALT is an allele of the record, `*` and `.` included: Ensembl VEP builds
    // the record's allele string from every ALT (Parser/VCF.pm:469) and
    // StructuralVariationOverlap.pm:64-77 makes one allele per `/`-part of it with no
    // filter (VariationFeatureOverlap.pm:483-484 filters `*` on the sequence path
    // alone), so `G/*` and `G/.` under a true INFO/SVTYPE each write two rows. A
    // record whose only ALT is `.` never reaches here (`vcf_parser.rs` reads it as a
    // non-variant record first).
    let classes = per_alt_classes(&alts, record_class, svtype);
    let svlen_text = info.value_written("SVLEN");
    let end_text = info.value_written("END");
    let ci_pos = parse_ci_field(info.value("CIPOS"));
    let ci_end = parse_ci_field(info.value("CIEND"));

    for (k, &(allele_idx, alt, sv_class)) in classes.iter().enumerate() {
        // SVLEN=0 or END=POS leaves the pair `start = POS + 1, end = POS` Perl carries
        // for every class: the event between the two bases, read with the predicates
        // of a sequence insertion.
        let sv_start = structural_start(pos, sv_class);
        if sv_start > MAX_COORD {
            return Err(SpanError::StartPastCeiling(pos).into());
        }
        let sv_end_val = structural_end(pos, svlen_text, end_text, SpanClass::of(sv_class))?;

        let mut variant = InputVariant::new(
            chr.to_string(),
            sv_start,
            sv_end_val,
            ref_allele.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant.variant_class = sv_class;
        variant.is_structural = true;
        // Ensembl VEP's `vep_skip` mark (Parser/VCF.pm:575): a type `%SO_TERMS` lacks
        // (`:477-481`) or a deletion without an extent (`:555-558`).
        variant.vep_skip = record_term.is_none()
            || deletion_looks_incomplete(
                &deletion_test_term,
                sv_start,
                sv_end_val,
                svlen_text,
                end_text,
            );
        variant.record_so_term = record_term.map(|(term, _)| term.to_string());
        variant.uploaded_allele_string = Some(uploaded.clone());
        variant.record_allele_string_multi = record_string.clone();
        variant.sv_end = Some(sv_end_val);
        variant.sv_type = svtype.map(str::to_string);
        variant.sv_len = svlen_text.and_then(|text| svlen_values(text).first().copied());
        variant.ci_pos = ci_pos;
        variant.ci_end = ci_end;
        variant.allele_index = allele_idx;
        // ALLELE_NUM numbers the parts of Ensembl VEP's allele string from 1
        // (StructuralVariationOverlap.pm:67-76), so a breakpoint record's reference
        // base shifts its ALTs by one (`N/<BND>` prints 2 for the `<BND>`); a record
        // with parsed breakends numbers its own breakend 0 and the mates from 1
        // (`:89-101`).
        variant.allele_num = breakpoint
            .as_ref()
            .map(|record| allele_idx as u32 + 1 + u32::from(record.reference_prefixed));
        // A translocation allele of a record with no parsed breakend is a single
        // breakend: the native `N.` and `.N` forms, and a symbolic `<BND>` without a
        // mate, which Perl holds as `N/<BND>` with no breakend string
        // (Parser/VCF.pm:499, `StructuralVariationFeature.pm` `_parse_breakends`).
        variant.is_single_breakend = sv_class == VariantClass::Translocation;
        match breakpoint
            .as_ref()
            .filter(|record| !record.mates.is_empty())
        {
            // The breakends come from the record's whole allele string, so they are
            // built once, at its first ALT: the own breakend ahead of the mates, as
            // `StructuralVariationOverlap::new` annotates `($vf, @$breakends)`. An ALT
            // that is not a bracket (`.N` beside `N[21:x[`) lends the own breakend its
            // base and gets no allele of its own.
            Some(record) => {
                if k == 0 {
                    mark_paired_breakend_record(&mut variant, chr, info, id);
                    variant.allele_num = Some(0);
                    variants.push(own_breakend(&variant, ref_allele, &record.own));
                    for (number, mate) in (1..).zip(&record.mates) {
                        let mut variant = variant.clone();
                        variant.allele_num = Some(number);
                        variant.alt_alleles = vec![mate.text.as_bytes().to_vec()];
                        variant.allele_string = format!("{ref_allele}/{}", mate.text);
                        // The mate keeps its written name for the cache lookup to resolve:
                        // Ensembl VEP fetches a mate's region through `get_source_chr_name`
                        // (AnnotationSource.pm:251-254, BaseVEP.pm:689-738: the name, its
                        // synonyms, the `chr` rule), which has no `M` to `MT` rename
                        // (Parser.pm:794-796 renames the record's own name alone).
                        variant.mate_chr = Some(mate.chromosome.clone());
                        variant.mate_pos = Some(mate.position);
                        variant.is_single_breakend = false;
                        variant.allele_index = mate.allele_index.unwrap_or(allele_idx);
                        variants.push(variant);
                    }
                }
            }
            None => variants.push(variant),
        }
    }

    if variants.is_empty() {
        return Err(StructuralRecordError::UnsupportedType(
            sv_type_abbreviation(&joined_alts, svtype),
        ));
    }
    Ok(variants)
}

/// The record's own breakend: a copy of `variant` whose only allele is `single_alt`
/// (`t.` or `.t`) over `ref_str`, flagged single, with no mate coordinates.
fn own_breakend(variant: &InputVariant, ref_str: &str, single_alt: &str) -> InputVariant {
    let mut own = variant.clone();
    own.alt_alleles = vec![single_alt.as_bytes().to_vec()];
    own.allele_string = format!("{ref_str}/{single_alt}");
    own.mate_chr = None;
    own.mate_pos = None;
    own.is_single_breakend = true;
    own
}

/// One bracket part of a breakpoint record's allele string: the text Ensembl VEP
/// prints as the Allele, the mate chromosome as written, the mate position, and the
/// ALT the part came from (`None` for the part INFO/CHR2 appends).
struct BreakendMate {
    text: String,
    chromosome: String,
    position: u64,
    allele_index: Option<usize>,
}

/// The allele string Ensembl VEP holds for a breakpoint record and what
/// `_parse_breakends` reads from it (Parser/VCF.pm:484-499,
/// StructuralVariationFeature.pm:1262-1323).
struct BreakendRecord {
    /// The ALT side of the string: the ALTs joined by `/`, and with INFO/CHR2 written
    /// and END2 true (else END written), `N[CHR2:pos[` in place of a lone `<BND>` or
    /// appended to anything else.
    alts: String,
    /// Whether the record's string is the reference base joined to `alts`
    /// (Parser/VCF.pm:499, `$alt = $ref . "/$alt" unless $alt =~ /^\.|\.$/`): false
    /// when the ALT side opens or closes with `.`, a native single breakend kept bare.
    reference_prefixed: bool,
    /// The record's own breakend, `t.`, or `.t` when any mate is left-placed; `t` is
    /// the base of the last part the bracket match refuses, the reference base
    /// included, when that part is a bare or dot-bounded ACGTN run, and `N` otherwise.
    own: String,
    /// One per part the breakend regex matches, in order.
    mates: Vec<BreakendMate>,
}

fn breakend_record(
    ref_allele: &str,
    alts: &[(usize, &str)],
    info: &InfoFields<'_>,
) -> BreakendRecord {
    let mut parts: Vec<(Option<usize>, String)> = alts
        .iter()
        .map(|&(idx, alt)| (Some(idx), alt.to_string()))
        .collect();
    // Parser/VCF.pm:491-497: `defined $info->{CHR2}` holds for `CHR2=` (Perl's split
    // keeps the empty field) and fails for a bare `CHR2`; the position is END2 when
    // Perl reads it as true (not empty, not `0`), else END as written, empty or `0`
    // included.
    if let Some(chr2) = info.value_written("CHR2") {
        let mate_pos = info
            .value("END2")
            .filter(|end2| perl_true(end2))
            .or_else(|| info.value_written("END"));
        if let Some(mate_pos) = mate_pos {
            let (allele_index, lone_bnd) = match parts.as_slice() {
                [(idx, alt)] if is_bare_bnd_alt(alt) => (*idx, true),
                _ => (None, false),
            };
            if lone_bnd {
                parts.clear();
            }
            parts.push((allele_index, format!("N[{chr2}:{mate_pos}[")));
        }
    }
    let alts = parts
        .iter()
        .map(|(_, part)| part.as_str())
        .collect::<Vec<_>>()
        .join("/");
    let reference_prefixed = !alts.starts_with('.') && !alts.ends_with('.');
    let reference_part = reference_prefixed.then(|| (None, ref_allele.to_string()));
    let mut base: Option<&str> = None;
    let mut any_left = false;
    let mut mates = Vec::new();
    for (allele_index, part) in reference_part.iter().chain(parts.iter()) {
        match bracket_breakend(part) {
            Some((chromosome, position)) => {
                any_left |= part.starts_with(['[', ']']);
                mates.push(BreakendMate {
                    text: part.clone(),
                    chromosome: chromosome.to_string(),
                    position,
                    allele_index: *allele_index,
                });
            }
            None => base = single_breakend_base(part),
        }
    }
    let base = base.unwrap_or("N");
    let own = if any_left {
        format!(".{base}")
    } else {
        format!("{base}.")
    };
    BreakendRecord {
        alts,
        reference_prefixed,
        own,
        mates,
    }
}

/// Perl's `/^<?BND>?$/i` (Parser/VCF.pm:495).
fn is_bare_bnd_alt(alt: &str) -> bool {
    let inner = alt.strip_prefix('<').unwrap_or(alt);
    inner
        .strip_suffix('>')
        .unwrap_or(inner)
        .eq_ignore_ascii_case("BND")
}

/// `_parse_breakends`' bracket match (StructuralVariationFeature.pm:1276,
/// `[A-Za-z]*\s*[\[\]]\s*(\S+)\s*:\s*([0-9]+)\s*[\[\]]`): a bracket, a chromosome name
/// of at least one non-space character, `:`, a run of digits, a bracket, with whitespace
/// allowed around the name, the colon and the digits (`N[ 21 : 30425700 [` names the
/// mate `21:30425700`); the chromosome is the longest name that leaves a digit run
/// before the second bracket. A position that is not digits alone (`-7`, `+5`, `abc`), an
/// empty chromosome (`N[:100[`) or a space inside one names no breakend. The position
/// saturates at [`MAX_COORD`]: the text stays as written in the Allele.
fn bracket_breakend(part: &str) -> Option<(&str, u64)> {
    let space = |c: char| c.is_ascii_whitespace();
    let open = part.find(['[', ']'])?;
    let inner = &part[open + 1..];
    let close = inner.find(['[', ']'])?;
    let inner = inner[..close].trim_matches(space);
    let digits_start = inner.trim_end_matches(|c: char| c.is_ascii_digit()).len();
    let digits = &inner[digits_start..];
    let chromosome = inner[..digits_start]
        .trim_end_matches(space)
        .strip_suffix(':')?
        .trim_end_matches(space);
    if digits.is_empty() || chromosome.is_empty() || chromosome.contains(space) {
        return None;
    }
    let position = digits
        .parse::<u64>()
        .map_or(MAX_COORD, |p| p.min(MAX_COORD));
    Some((chromosome, position))
}

/// `_parse_breakends`' reading of a part its bracket match refuses
/// (StructuralVariationFeature.pm:1281, `^[\.]?([ACGTN]+)[\.]?$`): the ACGTN run of a
/// bare base or a dot-bounded single breakend, and `None` for anything else, which
/// Perl assigns as undef.
fn single_breakend_base(part: &str) -> Option<&str> {
    let run = part.strip_prefix('.').unwrap_or(part);
    let run = run.strip_suffix('.').unwrap_or(run);
    (!run.is_empty()
        && run
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')))
    .then_some(run)
}

/// Names the mate of a breakend record that has one: INFO/MATEID, else PARID, else
/// the record's ID or `chr:start` (gnomAD `<BND>` records carry no MATEID). A set
/// `mate_id` is how `display_allele()` and the breakend predicates tell the record's
/// own breakend, written `N.` or `.N` (`StructuralVariationFeature.pm`
/// `_parse_breakends` names it from the bracket or synthesised mate allele alone),
/// from a record whose only ALT is a single breakend, whose Allele is the bare class
/// term `chromosome_breakpoint` (`OutputFactory.pm:2029`,
/// `$vfoa->{breakend}->{string} || $svf->class_SO_term`). Ensembl VEP never reads
/// MATEID, so one on a record without a bracket mate names nothing and is not kept.
fn mark_paired_breakend_record(
    variant: &mut InputVariant,
    chr: &str,
    info: &InfoFields<'_>,
    id: Option<&str>,
) {
    let mate_id = info
        .value("MATEID")
        .or_else(|| info.value("PARID"))
        .or(id)
        .map_or_else(|| format!("{}:{}", chr, variant.start), str::to_string);
    variant.mate_id = Some(mate_id);
}

/// The allele string `_expand_tandem_repeat_allele_string` (`Parser/VCF.pm:382-431`)
/// builds for a `<CNV:TR>` record once the reference sequence over its run is to hand:
/// `reference` followed by one explicit repeat sequence per INFO/RN entry, joined by
/// `/`. INFO/RN gives each allele its number of repeat sequences (one each for the
/// record's `alleles` tandem-repeat ALTs when RN is absent or false, `:394-397`; a
/// `.` entry adds no allele, `:410`; the count is a Perl range, truncated), consumed in
/// order from the RUS, RUC and RB lists (`:412-414`): a `.` unit is the one base `N`
/// (`:413`), a unit past the RUS list is Perl's undef, no bases; each sequence is its
/// unit repeated its RUC entry times whenever the record carries RUC, else its RB entry
/// divided by the unit's length (`:414`), the count numified and truncated as `x`
/// truncates (`:419`), a missing count no repeats. `None` keeps the symbolic reading:
/// without RUS or without RUC and RB (`:385`), when any sequence's unit length times
/// its count exceeds `max_sv_size` (`:403`, `:417-418`, `:428`; `None`, the flag's
/// `-1`, oversizes every sequence and `0` falls to 5000), and where Perl's RB division
/// by a unit of no bases ends its run. The string is Perl's: an allele of no sequences
/// is an empty part (`RN=0,1` gives `REF//<repeats>`), which the reader of the string
/// takes as a `-` allele (VariationFeatureOverlap.pm:466).
pub fn tandem_repeat_allele_string(
    info: &InfoFields<'_>,
    alleles: usize,
    max_sv_size: Option<u64>,
    reference: &str,
) -> Option<String> {
    let value = |key: &str| info.value(key).filter(|v| perl_true(v));
    // The lists as Perl's `split` reads the values (`:395-400`): a trailing empty
    // field is no field, so `RUC=,` is a true value whose list is empty, and the
    // record then counts by RB (`:401`, `:414`).
    let list = |key: &str| -> Vec<&str> {
        value(key)
            .map(|v| perl_split(v, &[',']).collect())
            .unwrap_or_default()
    };
    if value("RUS").is_none() || (value("RUC").is_none() && value("RB").is_none()) {
        return None;
    }
    let rus = list("RUS");
    let ruc = list("RUC");
    let rb = list("RB");
    // RN defaults to one sequence per allele only when the value is false (`:394`),
    // so `RN=,`, a true value whose list is empty, builds no allele.
    let counts: Vec<Option<f64>> = match value("RN") {
        None => vec![Some(1.0); alleles],
        Some(rn) => perl_split(rn, &[','])
            .map(|entry| (entry != ".").then(|| perl_numeric_f64(entry).trunc()))
            .collect(),
    };
    let cap = match max_sv_size {
        None => -1.0,
        Some(0) => 5000.0,
        Some(n) => n as f64,
    };
    let count_list = if ruc.is_empty() { &rb } else { &ruc };
    let (mut next_unit, mut next_count) = (0usize, 0usize);
    let mut repeats: Vec<String> = Vec::new();
    for records in counts.into_iter().flatten() {
        let mut repeat = String::new();
        let mut sequences = 0.0;
        // Once both lists are shifted out every further sequence is undef times undef,
        // no bases, so the walk stops there where Perl's range would run on.
        while sequences < records && (next_unit < rus.len() || next_count < count_list.len()) {
            sequences += 1.0;
            let unit: &str = match rus.get(next_unit).copied() {
                Some(".") => "N",
                Some(unit) => unit,
                None => "",
            };
            next_unit += 1;
            let entry = count_list
                .get(next_count)
                .map_or(0.0, |v| perl_numeric_f64(v));
            next_count += 1;
            let count = if ruc.is_empty() {
                if unit.is_empty() {
                    return None;
                }
                entry / unit.len() as f64
            } else {
                entry
            };
            let bases = unit.len() as f64 * count;
            // A product past the address space is Perl's out-of-memory end.
            if bases > cap || bases > isize::MAX as f64 {
                return None;
            }
            repeat.push_str(&unit.repeat(count.trunc().max(0.0) as usize));
        }
        repeats.push(repeat);
    }
    Some(format!("{reference}/{}", repeats.join("/")))
}

/// Parse a VCF confidence interval INFO field (CIPOS or CIEND).
///
/// Format: two comma-separated integers, e.g., "-100,50".
/// Returns `Some((lo, hi))` on success, `None` if missing or malformed.
pub fn parse_ci_field(value: Option<&str>) -> Option<(i64, i64)> {
    let s = value?;
    let (a, b) = s.split_once(',')?;
    Some((a.parse().ok()?, b.parse().ok()?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sv_routing::written_alts;

    /// `_parse_breakends`' match allows whitespace around the brackets, the name and the
    /// colon (`\s*`): Ensembl VEP 116.2 names the mate `21:30425700` for `N[ 21 : 30425700 [`
    /// and for the `N[ 21:30425700[` that `CHR2= 21` builds, printing the text as written;
    /// a space inside the name, an empty name or a signed position names no mate.
    #[test]
    fn a_bracket_mate_is_read_through_the_whitespace_perl_allows() {
        assert_eq!(
            bracket_breakend("N[ 21 : 30425700 ["),
            Some(("21", 30425700))
        );
        assert_eq!(bracket_breakend("N[ 21:30425700["), Some(("21", 30425700)));
        assert_eq!(bracket_breakend("]21:100]A"), Some(("21", 100)));
        assert_eq!(bracket_breakend("N[HLA-A:01:100["), Some(("HLA-A:01", 100)));
        for part in [
            "N[2 1:100[",
            "N[:100[",
            "N[21:+5[",
            "N[21:-7[",
            "N[21:abc[",
            "N[21:[",
            "<BND>",
        ] {
            assert_eq!(bracket_breakend(part), None, "{part}");
        }
    }

    /// BaseVCF4.pm:582's two-element list assignment keeps the text between the first
    /// and second `=`, so `END2=30425700=5` names the mate at 30425700 (Ensembl VEP 116.2
    /// writes `N[21:30425700[` beside the own `N.`).
    #[test]
    fn an_info_value_ends_at_a_second_equals_sign() {
        let record = breakend_record(
            "N",
            &written_alts("<BND>"),
            &InfoFields::new("SVTYPE=BND;CHR2=21;END2=30425700=5"),
        );
        assert_eq!(record.alts, "N[21:30425700[");
        assert_eq!(record.mates.len(), 1);
        assert_eq!(record.mates[0].position, 30425700);
        assert!(record.reference_prefixed);
        let spaced = breakend_record(
            "N",
            &written_alts("<BND>"),
            &InfoFields::new("SVTYPE=BND;CHR2= 21;END2=30425700"),
        );
        assert_eq!(spaced.alts, "N[ 21:30425700[");
        assert_eq!(
            (
                spaced.mates[0].chromosome.as_str(),
                spaced.mates[0].position
            ),
            ("21", 30425700)
        );
        let dot = breakend_record("N", &written_alts(".N"), &InfoFields::new("SVTYPE=BND"));
        assert!(!dot.reference_prefixed);
    }

    /// A mate's chromosome is carried as written, `chr` prefix and `M` included, for
    /// the cache lookup to resolve (`get_source_chr_name`, AnnotationSource.pm:251-254);
    /// the record's own name is normalised on construction.
    #[test]
    fn a_breakend_mate_keeps_its_written_chromosome_name() {
        let bracket = structural_variants(
            "chr21",
            30425600,
            "N",
            &written_alts("N[chrM:100["),
            &InfoFields::new("SVTYPE=BND"),
            None,
        )
        .unwrap();
        assert_eq!(bracket.len(), 2);
        assert_eq!(bracket[0].chr, "21");
        assert!(bracket[0].mate_chr.is_none());
        assert_eq!(bracket[1].mate_chr.as_deref(), Some("chrM"));
        assert_eq!(bracket[1].mate_pos, Some(100));
        let chr2 = structural_variants(
            "21",
            30425600,
            "N",
            &written_alts("<BND>"),
            &InfoFields::new("SVTYPE=BND;CHR2=chrX;END2=50000000"),
            None,
        )
        .unwrap();
        assert_eq!(chr2[1].mate_chr.as_deref(), Some("chrX"));
        assert_eq!(chr2[1].alt_allele(), b"N[chrX:50000000[");
    }

    /// `_expand_tandem_repeat_allele_string` reads its counts as Perl numifies them and
    /// repeats with `x`: `"AGC" x "2abc"` is six bases, `(1..1.9)` walks one sequence, a
    /// shifted-out RUC entry (`RUC=1,`) is undef and repeats nothing, a unit past the RUS
    /// list is undef and adds no bases; a unit of no bases (`RUS=,AC`) is an empty part
    /// under RUC (`R/`, which Perl's `split` reads as no allele) and a division by zero
    /// under RB, which keeps the symbolic reading; a false RN (`RN=0`) is an absent one; `RN=0,1` gives its first
    /// allele no sequence (Ensembl VEP 116.2 reads `21 30425587 x N <CNV:TR>,<CNV:TR>`
    /// with `SVLEN=9;RN=0,1;RUS=AGC;RUC=4` and a FASTA as `AGCAGCAGC/-`, the run deleted)
    /// and its second one sequence of four units; a count the list lacks is no repeats.
    /// A RUC of empty fields alone (`RUC=,`) is a true value whose list Perl's `split`
    /// leaves empty, so the record counts by RB, and by nothing without one.
    #[test]
    fn repeat_counts_are_read_as_perl_numifies_them() {
        // The parts as Perl's `split` reads the string: a trailing empty part is none.
        let alleles = |info: &str, alleles: usize| {
            tandem_repeat_allele_string(&InfoFields::new(info), alleles, Some(10_000_000), "R").map(
                |s| {
                    crate::sv_routing::perl_split(&s, &['/'])
                        .skip(1)
                        .map(str::to_string)
                        .collect::<Vec<_>>()
                },
            )
        };
        let one = |info: &str| alleles(info, 1);
        assert_eq!(one("RUS=AGC;RUC=2abc"), Some(vec!["AGCAGC".into()]));
        assert_eq!(one("RUS=AGC;RUC=1.9"), Some(vec!["AGC".into()]));
        assert_eq!(
            one("RN=1.9;RUS=AGC,AGC;RUC=2,2"),
            Some(vec!["AGCAGC".into()])
        );
        assert_eq!(one("RN=2;RUS=AGC,AGC;RUC=1,"), Some(vec!["AGC".into()]));
        assert_eq!(one("RN=2;RUS=AGC;RUC=1,2"), Some(vec!["AGC".into()]));
        assert_eq!(one("RUS=,AC;RUC=5"), Some(vec![]), "`R/`: no allele");
        assert_eq!(one("RUS=,AC;RB=12"), None);
        assert_eq!(one("RUS=AGC;RUC=0"), None);
        assert_eq!(one("RUS=AGC;RUC=abc"), Some(vec![]), "`R/`: no allele");
        assert_eq!(
            one("RN=0;RUS=AGC;RUC=4"),
            Some(vec!["AGCAGCAGCAGC".into()]),
            "a false RN is an absent one"
        );
        assert_eq!(
            alleles("RN=0,1;RUS=AGC;RUC=4", 2),
            Some(vec!["".into(), "AGCAGCAGCAGC".into()])
        );
        assert_eq!(
            alleles("RN=1,0;RUS=AGC;RUC=4", 2),
            Some(vec!["AGCAGCAGCAGC".into()]),
            "Perl's split drops the trailing empty part"
        );
        assert_eq!(
            alleles("RUS=AGC,AGC;RUC=4", 2),
            Some(vec!["AGCAGCAGCAGC".into()]),
            "the second allele's count is shifted out: an empty trailing part"
        );
        assert_eq!(
            alleles("RN=1,1;RUS=AGC,CGA;RB=12", 2),
            Some(vec!["AGCAGCAGCAGC".into()]),
            "a shifted-out RB entry is undef, no repeats"
        );
        assert_eq!(one("RN=abc;RUS=AGC;RUC=4"), Some(vec![]));
        assert_eq!(one("RUS=AGC;RB=10.5"), Some(vec!["AGCAGCAGC".into()]));
        assert_eq!(
            one("RUS=AGC;RUC=,;RB=12"),
            Some(vec!["AGCAGCAGCAGC".into()]),
            "an empty RUC list counts by RB"
        );
        assert_eq!(one("RUS=AGC;RUC=,"), Some(vec![]), "`R/`: no allele");
        assert_eq!(one("RUS=AGC"), None);
        assert_eq!(
            one("RN=18446744073709551615;RUS=AGC;RUC=4"),
            Some(vec!["AGCAGCAGCAGC".into()]),
            "a count past the lists consumes nothing more"
        );
    }

    /// The allele string Ensembl VEP 116.2 builds over the run `TCATCATCA` of
    /// `21 13620283 x G <CNV:TR>` with a FASTA: `RUC=4` and `RB=12` give four units,
    /// `RUC=4.9` four (`x` truncates), `RN=2;RUS=TCA,TC;RUC=3,1` the two sequences
    /// joined, `RUS=.` the base `N`, `RUS=tca` the unit as written; two alleles give
    /// two parts, a `.` RN entry none and a zero one an empty part; `RN=0`, a false
    /// value, defaults to one sequence per allele where `RN=,`, a true value whose
    /// list is empty, builds no allele (`REF/`); the string is `None` without RUS or
    /// a count list, under `--max_sv_size -1`, under a limit a sequence exceeds
    /// (`RUC=4` at 10, `RUC=2` not) and where RB divides by an empty unit.
    #[test]
    fn tandem_repeat_allele_string_follows_the_perl_expansion() {
        let run = "TCATCATCA";
        let expand = |info: &str, alleles: usize, max: Option<u64>| {
            tandem_repeat_allele_string(&InfoFields::new(info), alleles, max, run)
        };
        let default = Some(DEFAULT_MAX_SV_SIZE);
        let four = Some("TCATCATCA/TCATCATCATCA".to_string());
        assert_eq!(expand("SVLEN=9;RN=1;RUS=TCA;RUC=4", 1, default), four);
        assert_eq!(expand("SVLEN=9;RUS=TCA;RB=12", 1, default), four);
        assert_eq!(expand("SVLEN=9;RN=1;RUS=TCA;RUC=4.9", 1, default), four);
        assert_eq!(
            expand("SVLEN=9;RUS=TCA;RUC=5", 1, default),
            Some("TCATCATCA/TCATCATCATCATCA".into())
        );
        assert_eq!(
            expand("SVLEN=9;RN=2;RUS=TCA,TC;RUC=3,1", 1, default),
            Some("TCATCATCA/TCATCATCATC".into())
        );
        assert_eq!(
            expand("SVLEN=9;RN=1;RUS=.;RUC=4", 1, default),
            Some("TCATCATCA/NNNN".into())
        );
        assert_eq!(
            expand("SVLEN=9;RN=1;RUS=tca;RUC=4", 1, default),
            Some("TCATCATCA/tcatcatcatca".into())
        );
        assert_eq!(
            expand("SVLEN=9,9;RN=1,1;RUS=TCA,TCA;RUC=4,2", 2, default),
            Some("TCATCATCA/TCATCATCATCA/TCATCA".into())
        );
        assert_eq!(
            expand("SVLEN=9;RUS=TCA,TCA;RUC=4,2", 2, default),
            Some("TCATCATCA/TCATCATCATCA/TCATCA".into())
        );
        assert_eq!(expand("RN=.,1;RUS=TCA;RUC=4", 2, default), four);
        assert_eq!(
            expand("RN=0,1;RUS=TCA;RUC=4", 2, default),
            Some("TCATCATCA//TCATCATCATCA".into())
        );
        assert_eq!(expand("SVLEN=9;RN=0;RUS=TCA;RUC=4", 1, default), four);
        assert_eq!(
            expand("SVLEN=9;RN=,;RUS=TCA;RUC=4", 1, default),
            Some("TCATCATCA/".into())
        );
        assert_eq!(expand("SVLEN=9", 1, default), None);
        assert_eq!(expand("SVLEN=9;RUS=TCA", 1, default), None);
        assert_eq!(expand("SVLEN=9;RUC=4", 1, default), None);
        assert_eq!(expand("SVLEN=9;RN=1;RUS=TCA;RUC=4", 1, None), None);
        assert_eq!(expand("SVLEN=9;RN=1;RUS=TCA;RUC=4", 1, Some(10)), None);
        assert_eq!(
            expand("SVLEN=9;RN=1;RUS=TCA;RUC=2", 1, Some(10)),
            Some("TCATCATCA/TCATCA".into())
        );
        assert_eq!(expand("SVLEN=9;RN=2;RUS=TCA;RB=9,3", 1, default), None);
        assert_eq!(
            expand("SVLEN=9;RN=1;RUS=TCA;RUC=2", 1, Some(0)),
            Some("TCATCATCA/TCATCA".into())
        );
    }

    /// `([0-9]+)` matches ASCII digits alone: a mate written with fullwidth digits
    /// (`21 30425650 x N <BND> . . SVTYPE=BND;CHR2=21;END2=\u{ff13}\u{ff10}`) names no
    /// breakend, as in Ensembl VEP 116.2, whose rows are the bare `chromosome_breakpoint`;
    /// the multi-byte digit before the closing bracket is no boundary to slice at.
    #[test]
    fn a_mate_position_in_fullwidth_digits_names_no_breakend() {
        assert_eq!(bracket_breakend("N[21:\u{ff13}\u{ff10}["), None);
        assert_eq!(bracket_breakend("N[21:30\u{ff10}["), None);
        assert_eq!(bracket_breakend("N[\u{ff12}1:30["), Some(("\u{ff12}1", 30)));
        let record = breakend_record(
            "N",
            &written_alts("<BND>"),
            &InfoFields::new("SVTYPE=BND;CHR2=21;END2=\u{ff13}\u{ff10}"),
        );
        assert!(record.mates.is_empty());
        assert_eq!(record.alts, "N[21:\u{ff13}\u{ff10}[");
    }

    /// A trailing empty ALT field is no field: `<DEL>,` under `SVTYPE=DEL` is one deletion
    /// allele, not an unsupported `DEL/`.
    #[test]
    fn a_trailing_empty_alt_names_no_allele() {
        let variants = structural_variants(
            "21",
            30425650,
            "N",
            &written_alts("<DEL>,"),
            &InfoFields::new("SVTYPE=DEL;END=30425949"),
            None,
        )
        .unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].record_so_term.as_deref(), Some("deletion"));
        assert!(!variants[0].vep_skip);
        assert_eq!(variants[0].uploaded_allele_string.as_deref(), Some("<DEL>"));
    }

    /// A `*` beside another ALT under a true INFO/SVTYPE is an allele of the structural
    /// record (StructuralVariationOverlap.pm:64-77 numbers every `/`-part of `G/*`):
    /// Ensembl VEP 116.2 writes two `deletion` rows per transcript for
    /// `21 30289183 x C G,* . . SVTYPE=DEL;END=30289190`, `ALLELE_NUM=1` and
    /// `ALLELE_NUM=2`, both at 21:30289184-30289190 with `frameshift_variant,stop_lost,
    /// feature_truncation` and the JSON `"allele_string":"G/*"`; `C *,G` the same, and
    /// so are `C G,.` and `C .,G` (two `deletion` rows, `ALLELE_NUM=1` and `2`, JSON
    /// `"allele_string":"G/."`), the `.` a `/`-part like any other.
    #[test]
    fn a_star_beside_an_alt_is_an_allele_of_the_structural_record() {
        let record = |alts: &str| {
            structural_variants(
                "21",
                30289183,
                "C",
                &written_alts(alts),
                &InfoFields::new("SVTYPE=DEL;END=30289190"),
                Some("x"),
            )
            .unwrap()
        };
        for (alts, part_index, part_text) in [
            ("G,*", 1, "*"),
            ("*,G", 0, "*"),
            ("G,.", 1, "."),
            (".,G", 0, "."),
        ] {
            let variants = record(alts);
            assert_eq!(variants.len(), 2, "{alts}");
            let part = &variants[part_index];
            assert_eq!(part.alt_allele(), part_text.as_bytes(), "{alts}");
            assert_eq!(part.allele_index, part_index, "{alts}");
            assert_eq!(part.display_allele(), "deletion", "{alts}");
            for variant in &variants {
                assert!(variant.is_structural && !variant.vep_skip, "{alts}");
                assert_eq!((variant.start, variant.end), (30289184, 30289190), "{alts}");
                assert_eq!(
                    variant.record_so_term.as_deref(),
                    Some("deletion"),
                    "{alts}"
                );
                assert_eq!(
                    variant.uploaded_allele_string.as_deref(),
                    Some(alts.replace(',', "/").as_str()),
                    "{alts}"
                );
            }
        }
    }

    /// Ensembl VEP 116.2 writes `21 30425650 x N N[21:30425700[,<NON_REF> . . .` as the
    /// breakpoint `N/N[21:30425700[/<NON_REF>` (Parser/VCF.pm:469 joins every ALT), its
    /// rows the own `N.` and the mate: the block part is in the string and names no row.
    /// A gVCF site `<DEL>,<NON_REF>` without SVTYPE is its `<DEL>` allele, a deletion.
    #[test]
    fn a_block_beside_a_bracket_alt_stays_in_the_record_string() {
        let variants = structural_variants(
            "21",
            30425650,
            "N",
            &written_alts("N[21:30425700[,<NON_REF>"),
            &InfoFields::new("."),
            None,
        )
        .unwrap();
        let alleles: Vec<&str> = variants
            .iter()
            .map(|v| std::str::from_utf8(v.alt_allele()).unwrap())
            .collect();
        assert_eq!(alleles, ["N.", "N[21:30425700["]);
        assert_eq!(
            variants[0].record_allele_string_multi.as_deref(),
            Some("N/N[21:30425700[/<NON_REF>")
        );
        assert_eq!(variants[0].mate_id.as_deref(), Some("21:30425651"));
        let site = structural_variants(
            "21",
            30425650,
            "N",
            &written_alts("<DEL>,<NON_REF>"),
            &InfoFields::new("END=30425949"),
            Some("site"),
        )
        .unwrap();
        assert_eq!(site.len(), 1);
        assert_eq!(site[0].record_so_term.as_deref(), Some("deletion"));
        assert_eq!(
            site[0].uploaded_allele_string.as_deref(),
            Some("<DEL>/<NON_REF>")
        );
    }

    /// A structural record at POS 2^62 - 1 has no first base below the ceiling; Ensembl
    /// VEP writes `21:4611686018427387904-4611686018427388203 deletion - intergenic_variant`
    /// for `21 4611686018427387903 x N <DEL> . . SVTYPE=DEL;SVLEN=-300`, a span vep-rs
    /// cannot hold, so the line fails instead of printing an end below its start; one
    /// base lower the record is read.
    #[test]
    fn a_first_base_past_the_ceiling_fails_the_record() {
        let at = |pos: u64, alt: &str, info: &str| {
            structural_variants(
                "21",
                pos,
                "N",
                &written_alts(alt),
                &InfoFields::new(info),
                None,
            )
        };
        for (alt, info) in [("<DEL>", "SVTYPE=DEL;SVLEN=-300"), ("<BND>", "SVTYPE=BND")] {
            let error = at(MAX_COORD, alt, info).unwrap_err();
            assert_eq!(
                error.to_string(),
                format!("Invalid POS: '{MAX_COORD}'"),
                "{alt}"
            );
            assert_eq!(
                at(MAX_COORD - 1, alt, info).unwrap()[0].start,
                MAX_COORD,
                "{alt}"
            );
        }
    }
}
