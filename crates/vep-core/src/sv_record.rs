// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The alleles a structural VCF record writes rows for, built as Ensembl VEP builds the
//! one StructuralVariationFeature that carries the record (`Parser/VCF.pm`
//! `create_StructuralVariationFeatures`, lines 466-575) and its breakends
//! (ensembl-variation `StructuralVariationFeature.pm` `_parse_breakends`, lines
//! 1262-1323), shared by every VCF reader: the record's type and class, each allele's
//! span, a breakpoint record's own breakend and mates, a tandem repeat's expanded
//! length, and the allele numbering Ensembl VEP prints.

use std::fmt;
use std::sync::Arc;

use crate::coordinate::{normalize_chromosome, MAX_COORD};
use crate::sv_routing::{
    joined, per_alt_classes, perl_true, sv_so_term_and_class, sv_type_abbreviation, variant_alts,
};
use crate::sv_span::{
    perl_numeric_f64, structural_end, structural_start, svlen_values, SpanClass, SpanError,
};
use crate::variant::{InputVariant, VariantClass};
use crate::vcf_info::{info_value, info_value_written};

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
/// the record's ALTs with their written positions (`sv_routing::written_alts`), `id`
/// the ID column when it names the record, and `max_sv_size` the run's `--max_sv_size`
/// (`None` for `-1`), which decides whether a `<CNV:TR>` allele's length is expanded or
/// kept symbolic. One variant per ALT that names an allele under the record's class
/// after the gVCF rule (`sv_routing::variant_alts`), or under its own class with the
/// `vep_skip` mark for a type `%SO_TERMS` lacks; a breakpoint record with parsed mates
/// writes its own breakend and one variant per mate instead, at its first ALT. The
/// record's allele strings keep every written ALT. The caller sets what the record's
/// line adds (`id`, `raw_input`, `input_line_name`).
pub fn structural_variants(
    chr: &str,
    pos: u64,
    ref_allele: &str,
    written: &[(usize, &str)],
    info_field: &str,
    max_sv_size: Option<u64>,
    id: Option<&str>,
) -> Result<Vec<InputVariant>, StructuralRecordError> {
    let svtype = info_value(info_field, "SVTYPE").filter(|v| perl_true(v));
    let mut alts = written.to_vec();
    variant_alts(&mut alts, svtype);
    // `*` and `.` are written as no allele here; both stay in the record's allele
    // string, which Ensembl VEP builds from every ALT (Parser/VCF.pm:469). Ensembl VEP
    // writes no row for `*` and reads a `.` beside another ALT as an allele of its own.
    let valid_alts: Vec<(usize, &str)> = alts
        .iter()
        .copied()
        .filter(|(_, alt)| *alt != "*" && *alt != ".")
        .collect();
    let mut variants = Vec::with_capacity(valid_alts.len());
    // One structural variant record has one class and one SO term, read from the
    // joined ALTs with INFO/SVTYPE standing in for a string outside the VCF 4.4 forms
    // (Parser.pm `get_SO_term`, Parser/VCF.pm:472-475), so `N.` with `SVTYPE=DEL`
    // is a deletion; the term labels the row of every ALT, each ALT keeping its own
    // bytes (`per_alt_classes`).
    let joined_alts = joined(&alts, "/");
    let record_term = sv_so_term_and_class(&sv_type_abbreviation(&joined_alts, svtype));
    let record_class = record_term.map(|(_, class)| class);
    // VEP's allele string for a structural variant record is every ALT of the
    // record joined by `/`, `*`, a `<NON_REF>` block and the sequence alleles of a
    // mixed record included (Parser/VCF.pm:469); a breakend record's carries the mate
    // INFO/CHR2 names and, as the record's string, the reference base in front
    // (`:491-499`).
    let breakpoint = (record_class == Some(VariantClass::Translocation))
        .then(|| breakend_record(ref_allele, written, info_field));
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
    let classes = per_alt_classes(&valid_alts, record_class, svtype);
    let tandem_repeat_alt_count = |classes: &[(usize, &str, VariantClass)]| {
        classes
            .iter()
            .filter(|(_, a, _)| a.to_ascii_uppercase().starts_with("<CNV:TR"))
            .count()
    };
    let record_tandem_repeat_alts = tandem_repeat_alt_count(&classes);
    let svlen_text = info_value_written(info_field, "SVLEN");
    let end_text = info_value_written(info_field, "END");
    let ci_pos = parse_ci_field(info_value(info_field, "CIPOS"));
    let ci_end = parse_ci_field(info_value(info_field, "CIEND"));

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
        variant.vep_skip = record_term.is_none();
        variant.record_so_term = record_term.map(|(term, _)| term.to_string());
        variant.uploaded_allele_string = Some(uploaded.clone());
        variant.record_allele_string_multi = record_string.clone();
        variant.sv_end = Some(sv_end_val);
        variant.sv_type = svtype.map(str::to_string);
        variant.sv_len = svlen_text.and_then(|text| svlen_values(text).first().copied());
        if sv_class == VariantClass::TandemRepeat {
            variant.tr_alt_bases = tandem_repeat_alt_bases(
                info_field,
                tandem_repeat_alt_count(&classes[..k]),
                record_tandem_repeat_alts,
                max_sv_size,
            );
        }
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
                    mark_paired_breakend_record(&mut variant, chr, info_field, id);
                    variant.allele_num = Some(0);
                    variants.push(own_breakend(&variant, ref_allele, &record.own));
                    for (number, mate) in (1..).zip(&record.mates) {
                        let mut variant = variant.clone();
                        variant.allele_num = Some(number);
                        variant.alt_alleles = vec![mate.text.as_bytes().to_vec()];
                        variant.allele_string = format!("{ref_allele}/{}", mate.text);
                        variant.mate_chr = Some(normalize_chromosome(&mate.chromosome));
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

fn breakend_record(ref_allele: &str, alts: &[(usize, &str)], info: &str) -> BreakendRecord {
    let mut parts: Vec<(Option<usize>, String)> = alts
        .iter()
        .map(|&(idx, alt)| (Some(idx), alt.to_string()))
        .collect();
    // Parser/VCF.pm:491-497: `defined $info->{CHR2}` holds for `CHR2=` (Perl's split
    // keeps the empty field) and fails for a bare `CHR2`; the position is END2 when
    // Perl reads it as true (not empty, not `0`), else END as written, empty or `0`
    // included.
    if let Some(chr2) = info_value_written(info, "CHR2") {
        let mate_pos = info_value(info, "END2")
            .filter(|end2| perl_true(end2))
            .or_else(|| info_value_written(info, "END"));
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
    info: &str,
    id: Option<&str>,
) {
    let mate_id = info_value(info, "MATEID")
        .or_else(|| info_value(info, "PARID"))
        .or(id)
        .map_or_else(|| format!("{}:{}", chr, variant.start), str::to_string);
    variant.mate_id = Some(mate_id);
}

/// The alternate allele's length in bases for the `ordinal`th `<CNV:TR>` allele of a
/// record carrying `alleles` of them: the length of the allele
/// `_expand_tandem_repeat_allele_string` (`Parser/VCF.pm:385-428`) builds. INFO/RN gives
/// each tandem-repeat allele its number of repeat sequences, a Perl range truncating it
/// (`RN=1.9` walks one; one each when RN is absent or false, `:394-397`; an allele past
/// the list's end has none, `:408`), consumed in order from the RUS, RUC and RB lists
/// (`:412-414`); a `.` unit is the one base `N` (`:413`); each sequence's count is its
/// RUC entry whenever the record carries RUC, else its RB entry divided by the unit's
/// length (`:414`), read as Perl numifies the text (`2abc` is 2) and truncated the way
/// `x` truncates (`:419`); a count the RUC or RB list lacks is Perl's undef, no repeats.
/// RUL is never read. An allele of no bases (`RN=0,1`'s first, `RUS=,AC;RUC=5`) is
/// `Some(0)`, the loss of the whole run, the direction of Perl's empty allele. `None`
/// keeps the symbolic reading: without RUS or without RUC and RB (`:385`), for an allele
/// with `.` in RN, with more sequences than the RUS list has units or with no count entry
/// for any of its sequences (`RUS=AGC,AGC;RUC=4` states one allele's length), for every
/// allele of a record whose RB divides by a unit of no bases (Perl dies), when the RN
/// counts through the allele sum past `usize::MAX` or its length past `u64::MAX`, and for
/// every allele of the record when any repeat sequence of the record is longer than
/// `max_sv_size` (`:403`, `:417-418`, `:428`): under `None`, the flag's `-1`, no sequence
/// is ever short enough.
pub fn tandem_repeat_alt_bases(
    info: &str,
    ordinal: usize,
    alleles: usize,
    max_sv_size: Option<u64>,
) -> Option<u64> {
    // `$info->{RUS} && ($info->{RUC} || $info->{RB})` and `$info->{RN} ? ... : (1) x N`
    // (`:385`, `:394-400`): a value Perl reads as false is no list.
    let list = |key: &str| -> Vec<&str> {
        info_value(info, key)
            .filter(|v| perl_true(v))
            .map(|v| v.split(',').collect())
            .unwrap_or_default()
    };
    let rus = list("RUS");
    let ruc = list("RUC");
    let rb = list("RB");
    if rus.is_empty() || (ruc.is_empty() && rb.is_empty()) {
        return None;
    }
    let rn = list("RN");
    // `for ( 1 .. $records )` (`:411`): the range truncates the count and runs no
    // iteration below one; `.` skips the allele (`:410`).
    let sequences = |allele: usize| -> Option<usize> {
        if rn.is_empty() {
            return Some(1);
        }
        let entry = *rn.get(allele)?;
        (entry != ".").then(|| perl_numeric_f64(entry).trunc().max(0.0) as usize)
    };
    // The `j`th repeat sequence's unit length and count, `(length($seq), $num)`: a count
    // list shifted past its end yields undef, a count of none.
    let sequence = |j: usize| -> (f64, f64) {
        let unit_len = match rus[j] {
            "." => 1.0,
            unit => unit.len() as f64,
        };
        let count = if ruc.is_empty() {
            rb.get(j).map_or(0.0, |v| perl_numeric_f64(v)) / unit_len
        } else {
            ruc.get(j).map_or(0.0, |v| perl_numeric_f64(v))
        };
        (unit_len, count)
    };
    // `$self->param('max_sv_size') || 5000` (`:403`): `-1` passes `||` unchanged, so
    // every sequence is oversized under it, and `0` falls to 5000.
    let cap = match max_sv_size {
        None => return None,
        Some(0) => 5000.0,
        Some(n) => n as f64,
    };
    let record_alleles = if rn.is_empty() { alleles } else { rn.len() };
    let mut total = 0usize;
    for allele in 0..record_alleles {
        total = total.checked_add(sequences(allele).unwrap_or(0))?;
    }
    // Under RB, `shift(@RB) / length($seq)` (`:414`) on a unit of no bases, written or
    // past the RUS list, divides by zero and ends the Perl run for every allele of the
    // record. `length($seq) * $num > $max_sv_size` (`:417`) tests the count before `x`
    // truncates it.
    if ruc.is_empty() && total > rus.len() {
        return None;
    }
    for j in 0..total.min(rus.len()) {
        let (unit_len, count) = sequence(j);
        if (ruc.is_empty() && unit_len == 0.0) || unit_len * count > cap {
            return None;
        }
    }
    let mut first = 0usize;
    for earlier in 0..ordinal {
        first = first.checked_add(sequences(earlier)?)?;
    }
    let last = first.checked_add(sequences(ordinal)?)?;
    if last > rus.len() {
        return None;
    }
    // An allele none of whose sequences has a count entry states no length: Perl's
    // `shift` past the list yields undef, read as no repeats (`:414`, `:419`), the
    // empty part that `:428` joins and `split` drops.
    let counts = if ruc.is_empty() { &rb } else { &ruc };
    if first < last && counts.len() <= first {
        return None;
    }
    let mut bases = 0u64;
    for j in first..last {
        let (unit_len, count) = sequence(j);
        bases = bases.checked_add((count.trunc().max(0.0) * unit_len) as u64)?;
    }
    Some(bases)
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
            "SVTYPE=BND;CHR2=21;END2=30425700=5",
        );
        assert_eq!(record.alts, "N[21:30425700[");
        assert_eq!(record.mates.len(), 1);
        assert_eq!(record.mates[0].position, 30425700);
        assert!(record.reference_prefixed);
        let spaced = breakend_record(
            "N",
            &written_alts("<BND>"),
            "SVTYPE=BND;CHR2= 21;END2=30425700",
        );
        assert_eq!(spaced.alts, "N[ 21:30425700[");
        assert_eq!(
            (
                spaced.mates[0].chromosome.as_str(),
                spaced.mates[0].position
            ),
            ("21", 30425700)
        );
        let dot = breakend_record("N", &written_alts(".N"), "SVTYPE=BND");
        assert!(!dot.reference_prefixed);
    }

    /// `_expand_tandem_repeat_allele_string` reads its counts as Perl numifies them and
    /// repeats with `x`: `"AGC" x "2abc"` is six bases, `(1..1.9)` walks one sequence, a
    /// shifted-out RUC entry (`RUC=1,`) is undef and repeats nothing; an allele with more
    /// sequences than the RUS list has units states no length; a unit of no bases
    /// (`RUS=,AC`) is the empty allele under RUC, no bases and so the loss of the run, and
    /// a division by zero under RB, which keeps the symbolic reading; `RN=0,1` gives its
    /// first allele no sequence (Ensembl VEP 116.2 reads `21 30425587 x N <CNV:TR>,<CNV:TR>`
    /// with `SVLEN=9;RN=0,1;RUS=AGC;RUC=4` and a FASTA as `AGCAGCAGC/-`, the run deleted)
    /// and its second one sequence of four units.
    #[test]
    fn repeat_counts_are_read_as_perl_numifies_them() {
        let bases = |info: &str| tandem_repeat_alt_bases(info, 0, 1, Some(10_000_000));
        assert_eq!(bases("RUS=AGC;RUC=2abc"), Some(6));
        assert_eq!(bases("RUS=AGC;RUC=1.9"), Some(3));
        assert_eq!(bases("RN=1.9;RUS=AGC,AGC;RUC=2,2"), Some(6));
        assert_eq!(bases("RN=2;RUS=AGC,AGC;RUC=1,"), Some(3));
        assert_eq!(bases("RN=2;RUS=AGC;RUC=1,2"), None);
        assert_eq!(bases("RUS=,AC;RUC=5"), Some(0));
        assert_eq!(bases("RUS=,AC;RB=12"), None);
        assert_eq!(bases("RUS=AGC;RUC=0"), None);
        assert_eq!(bases("RUS=AGC;RUC=abc"), Some(0));
        assert_eq!(
            bases("RN=0;RUS=AGC;RUC=4"),
            Some(12),
            "a false RN is an absent one"
        );
        let two = |info: &str, ordinal: usize| {
            tandem_repeat_alt_bases(info, ordinal, 2, Some(10_000_000))
        };
        assert_eq!(two("RN=0,1;RUS=AGC;RUC=4", 0), Some(0));
        assert_eq!(two("RN=0,1;RUS=AGC;RUC=4", 1), Some(12));
        assert_eq!(two("RN=1,0;RUS=AGC;RUC=4", 1), Some(0));
        assert_eq!(two("RUS=AGC,AGC;RUC=4", 0), Some(12));
        assert_eq!(two("RUS=AGC,AGC;RUC=4", 1), None);
        assert_eq!(two("RN=1,1;RUS=AGC,CGA;RB=12", 1), None);
        assert_eq!(bases("RN=abc;RUS=AGC;RUC=4"), Some(0));
        assert_eq!(bases("RUS=AGC;RB=10.5"), Some(9));
        assert_eq!(bases("RUS=AGC"), None);
        assert_eq!(bases("RUS=AGC;RUC=4"), Some(12));
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
            "SVTYPE=BND;CHR2=21;END2=\u{ff13}\u{ff10}",
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
            "SVTYPE=DEL;END=30425949",
            Some(DEFAULT_MAX_SV_SIZE),
            None,
        )
        .unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].record_so_term.as_deref(), Some("deletion"));
        assert!(!variants[0].vep_skip);
        assert_eq!(variants[0].uploaded_allele_string.as_deref(), Some("<DEL>"));
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
            ".",
            Some(DEFAULT_MAX_SV_SIZE),
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
            "END=30425949",
            Some(DEFAULT_MAX_SV_SIZE),
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
                info,
                Some(DEFAULT_MAX_SV_SIZE),
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
