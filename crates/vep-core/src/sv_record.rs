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
    deletion_looks_incomplete, end_coordinate, perl_numeric_f64, structural_end, structural_start,
    svlen_values, SpanClass, SpanError,
};
use crate::variant::{InputVariant, VariantClass, WrittenEnd};
use crate::vcf_info::InfoFields;

/// Ensembl VEP's default `--max_sv_size` (`Config.pm:315`), in bases.
pub const DEFAULT_MAX_SV_SIZE: u64 = 10_000_000;

/// Why a structural record writes no allele.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
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
    // The length each allele states for itself is read where every written ALT of the
    // record is `<CNV:TR>`, a block or a `*` beside one included, so that the INFO lists
    // follow the ALTs one to one.
    let (stated_bases, as_long_as_the_run) = if written.iter().all(|&(_, alt)| alt == "<CNV:TR>")
        && classes
            .iter()
            .all(|&(_, _, class)| class == VariantClass::TandemRepeat)
    {
        (
            tandem_repeat_stated_bases(info, classes.len()),
            tandem_repeat_alleles_as_long_as_the_run(info, classes.len(), pos),
        )
    } else {
        (Vec::new(), Vec::new())
    };
    let svlen_text = info.value_written("SVLEN");
    let end_text = info.value_written("END");
    let ci_pos = parse_ci_field(info.value("CIPOS"));
    let ci_end = parse_ci_field(info.value("CIEND"));
    let imprecise = info.field("IMPRECISE").is_some();

    for (k, &(allele_idx, alt, sv_class)) in classes.iter().enumerate() {
        // SVLEN=0 or END=POS leaves the pair `start = POS + 1, end = POS` Ensembl VEP carries
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
        // An allele that states fewer bases than the run from `sv_start` to
        // `sv_end_val` is a loss of repeat units, whatever the reference sequence.
        variant.states_a_loss = stated_bases
            .get(k)
            .copied()
            .flatten()
            .is_some_and(|bases| sv_end_val >= sv_start && bases < sv_end_val - sv_start + 1);
        // An allele the record states to be exactly as long as that run leaves the
        // run's length as it is.
        variant.states_the_run_length = as_long_as_the_run.get(k).copied().unwrap_or(false);
        variant.record_so_term = record_term.map(|(term, _)| term.to_string());
        variant.uploaded_allele_string = Some(uploaded.clone());
        variant.record_allele_string_multi = record_string.clone();
        variant.sv_end = Some(sv_end_val);
        variant.sv_type = svtype.map(str::to_string);
        variant.sv_len = svlen_text.and_then(|text| svlen_values(text).first().copied());
        variant.ci_pos = ci_pos;
        variant.ci_end = ci_end;
        variant.imprecise = imprecise;
        variant.written_end = WrittenEnd::of(end_text, pos, sv_end_val);
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
        // mate, which Ensembl VEP holds as `N/<BND>` with no breakend string
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
                    variants.push(own_breakend(&variant, ref_allele, record));
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

/// The record's own breakend: a copy of `variant` whose only allele is the
/// record's own (`t.` or `.t`) over `ref_str`, flagged single, with no mate
/// coordinates of its own and the record's mates listed for the feature match.
fn own_breakend(variant: &InputVariant, ref_str: &str, record: &BreakendRecord) -> InputVariant {
    let mut own = variant.clone();
    own.alt_alleles = vec![record.own.as_bytes().to_vec()];
    own.allele_string = format!("{ref_str}/{}", record.own);
    own.mate_chr = None;
    own.mate_pos = None;
    own.is_single_breakend = true;
    own.breakend_mates = record
        .mates
        .iter()
        .map(|mate| (mate.chromosome.clone(), mate.position))
        .collect();
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

/// Ensembl VEP's `/^<?BND>?$/i` (Parser/VCF.pm:495).
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
/// Ensembl VEP assigns as undef.
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
/// term `chromosome_breakpoint` (`OutputFactory.pm:2038`,
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
/// `-1`, oversizes every sequence and `0` falls to 5000), and where Ensembl VEP's RB division
/// by a unit of no bases ends its run. The string is Ensembl VEP's: an allele of no sequences
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
        // no bases, so the walk stops there where Ensembl VEP's range would run on.
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
            // A product past the address space is Ensembl VEP's out-of-memory end.
            if bases > cap || bases > isize::MAX as f64 {
                return None;
            }
            repeat.push_str(&unit.repeat(count.trunc().max(0.0) as usize));
        }
        repeats.push(repeat);
    }
    Some(format!("{reference}/{}", repeats.join("/")))
}

/// The bases each alternate allele of a `<CNV:TR>` record states for itself, one
/// entry per tandem-repeat ALT in the record's order: the sum, over the repeat
/// sequences INFO/RN gives the allele (one where RN is absent), of the unit's length
/// times its RUC entry, or of its RB entry where the record carries no RUC, the
/// lists consumed in order as `tandem_repeat_allele_string` consumes them. `None`
/// where the record does not state the allele's length in whole units of a written
/// unit: no RUS, no RUC and no RB, a CIRUC or CIRB interval around the counts, an RN
/// that is not one whole number per allele, an allele of no sequence, a `.` or
/// empty unit, a count that is not a whole number, an RB that is no whole multiple
/// of its unit, or a list that ends before the allele's sequences do.
pub fn tandem_repeat_stated_bases(info: &InfoFields<'_>, alleles: usize) -> Vec<Option<u64>> {
    stated_alleles(info, alleles)
        .into_iter()
        .map(|allele| allele.map(|allele| allele.bases))
        .collect()
}

/// Whether each alternate allele of the `<CNV:TR>` record at VCF position `pos` is
/// exactly as long as the reference run the record spans, one entry per
/// tandem-repeat ALT in the record's order. The run is the bases from `pos + 1` to
/// the end SVLEN or END gives ([`structural_end`]), and an allele is as long as it
/// where all of these hold:
///
/// - the record states the allele's length ([`tandem_repeat_stated_bases`]) and
///   that length is the run's;
/// - the run is longer than one base. VCF 4.4 (section 5.7) writes a repeat the
///   reference lacks with an SVLEN of 1 at the base before the inserted sequence,
///   and a record with neither SVLEN nor END spans one base, so a span of one base
///   names no reference run of one base;
/// - every SVLEN value is the same number, and an END written beside SVLEN names
///   the run's end. SVLEN holds one reference length for each allele and Ensembl
///   VEP spans every allele of the record with the largest
///   (Parser/VCF.pm:505-511), so under unequal values the run is not each allele's
///   own reference;
/// - no other field states another length for the allele: each RUC entry of its
///   repeat sequences counts the sequence's units, each RB entry is the
///   sequence's bases, each RUL entry the length of the sequence's unit, the
///   record carries no RUB, which gives each unit a length of its own, and the
///   allele's entry of INFO/CN, its length over its reference's, is the number 1.
///
/// RUC, RB, RUL, RUB and CN count as written wherever the INFO column carries the
/// key at all, each list split on its commas with an entry at the place of its
/// sequence in RUS (CN: of its allele). A value of `0` or `.`, an empty value, a
/// bare flag and a list that ends before the allele's entries each leave the
/// allele unmarked.
pub fn tandem_repeat_alleles_as_long_as_the_run(
    info: &InfoFields<'_>,
    alleles: usize,
    pos: u64,
) -> Vec<bool> {
    let svlen = info.value_written("SVLEN");
    let end_text = info.value_written("END");
    let start = structural_start(pos, VariantClass::TandemRepeat);
    let end = structural_end(pos, svlen, end_text, SpanClass::Tandem).ok();
    let run = end.filter(|end| *end >= start).map(|end| end - start + 1);
    let longer_than_one_base = run.is_some_and(|run| run > 1);
    let one_reference_length = svlen.is_none_or(|text| {
        let values: Vec<f64> = perl_split(text, &[',']).map(perl_numeric_f64).collect();
        values.windows(2).all(|pair| pair[0] == pair[1])
    });
    let end_names_the_run = end_text.is_none_or(|text| end_coordinate(text) == end);
    let entries = |key: &str| -> Option<Vec<&str>> {
        info.value(key)
            .map(|text| perl_split(text, &[',']).collect())
    };
    let ruc = entries("RUC");
    let rb = entries("RB");
    let rul = entries("RUL");
    let no_rub = info.field("RUB").is_none();
    let copy_numbers = entries("CN");
    stated_alleles(info, alleles)
        .iter()
        .enumerate()
        .map(|(allele, stated)| {
            stated.as_ref().is_some_and(|stated| {
                let ruc_agrees = ruc.as_ref().is_none_or(|ruc| {
                    stated.sequences.iter().all(|sequence| {
                        ruc.get(sequence.place)
                            .copied()
                            .and_then(whole_number)
                            .and_then(|count| sequence.unit_bases.checked_mul(count))
                            == Some(sequence.bases)
                    })
                });
                let rb_agrees = rb.as_ref().is_none_or(|rb| {
                    stated.sequences.iter().all(|sequence| {
                        rb.get(sequence.place).copied().and_then(whole_number)
                            == Some(sequence.bases)
                    })
                });
                let rul_agrees = rul.as_ref().is_none_or(|rul| {
                    stated.sequences.iter().all(|sequence| {
                        rul.get(sequence.place).copied().and_then(whole_number)
                            == Some(sequence.unit_bases)
                    })
                });
                let one_copy = copy_numbers.as_ref().is_none_or(|copy_numbers| {
                    copy_numbers
                        .get(allele)
                        .and_then(|entry| entry.parse::<f64>().ok())
                        == Some(1.0)
                });
                Some(stated.bases) == run
                    && longer_than_one_base
                    && one_reference_length
                    && end_names_the_run
                    && ruc_agrees
                    && rb_agrees
                    && rul_agrees
                    && no_rub
                    && one_copy
            })
        })
        .collect()
}

/// A repeat sequence of an allele whose length its record states: its place in the
/// record's RUS, RUC, RB and RUL lists, the length of its unit and its own length,
/// both in bases.
struct StatedSequence {
    place: usize,
    unit_bases: u64,
    bases: u64,
}

/// An alternate allele of a `<CNV:TR>` record whose length the record states: its
/// bases and its repeat sequences in the record's order.
struct StatedAllele {
    bases: u64,
    sequences: Vec<StatedSequence>,
}

/// The number `text` writes in digits alone, where 64 bits hold it.
fn whole_number(text: &str) -> Option<u64> {
    (!text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit()))
        .then(|| text.parse().ok())
        .flatten()
}

/// Each alternate allele of a `<CNV:TR>` record as the record states it, by the
/// rules of [`tandem_repeat_stated_bases`]: `None` where the record does not state
/// the allele's length.
fn stated_alleles(info: &InfoFields<'_>, alleles: usize) -> Vec<Option<StatedAllele>> {
    let value = |key: &str| info.value(key).filter(|v| perl_true(v));
    let list = |key: &str| -> Vec<&str> {
        value(key)
            .map(|v| perl_split(v, &[',']).collect())
            .unwrap_or_default()
    };
    let unstated = || (0..alleles).map(|_| None).collect();
    if value("RUS").is_none()
        || (value("RUC").is_none() && value("RB").is_none())
        || value("CIRUC").is_some()
        || value("CIRB").is_some()
    {
        return unstated();
    }
    let sequences: Vec<u64> = match value("RN") {
        None => vec![1; alleles],
        Some(rn) => match perl_split(rn, &[','])
            .map(whole_number)
            .collect::<Option<Vec<u64>>>()
        {
            Some(counts) if counts.len() == alleles => counts,
            _ => return unstated(),
        },
    };
    let rus = list("RUS");
    let ruc = list("RUC");
    let rb = list("RB");
    let by_count = !ruc.is_empty();
    let counts = if by_count { &ruc } else { &rb };
    let mut next = 0usize;
    sequences
        .into_iter()
        .map(|own| {
            // An allele that asks for more sequences than the lists still hold states
            // no length and leaves none for the alleles after it.
            let Some(own) = usize::try_from(own)
                .ok()
                .filter(|own| *own <= rus.len().saturating_sub(next))
            else {
                next = rus.len();
                return None;
            };
            let mut bases = 0u64;
            let mut stated = true;
            let mut own_sequences = Vec::with_capacity(own);
            for (place, &unit) in rus.iter().enumerate().skip(next).take(own) {
                let sequence = match counts.get(place).copied().and_then(whole_number) {
                    Some(entry) if !unit.is_empty() && unit != "." => {
                        let unit_bases = unit.len() as u64;
                        let sequence_bases = if by_count {
                            unit_bases.checked_mul(entry)
                        } else {
                            (entry % unit_bases == 0).then_some(entry)
                        };
                        sequence_bases.map(|bases| StatedSequence {
                            place,
                            unit_bases,
                            bases,
                        })
                    }
                    _ => None,
                };
                match sequence
                    .as_ref()
                    .and_then(|sequence| bases.checked_add(sequence.bases))
                {
                    Some(sum) => bases = sum,
                    None => stated = false,
                }
                own_sequences.extend(sequence);
            }
            next += own;
            (stated && own > 0).then_some(StatedAllele {
                bases,
                sequences: own_sequences,
            })
        })
        .collect()
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
    /// colon (`\s*`): Ensembl VEP 116.3 names the mate `21:30425700` for `N[ 21 : 30425700 [`
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
    /// and second `=`, so `END2=30425700=5` names the mate at 30425700 (Ensembl VEP 116.3
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
        assert_eq!(bracket[0].breakend_mates, [("chrM".to_string(), 100)]);
        assert_eq!(bracket[1].mate_chr.as_deref(), Some("chrM"));
        assert_eq!(bracket[1].mate_pos, Some(100));
        assert!(bracket[1].breakend_mates.is_empty());
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
    /// under RB, which keeps the symbolic reading; a false RN (`RN=0`) is an absent one; `RN=0,1`
    /// gives its first allele no sequence (Ensembl VEP 116.3 reads
    /// `21 30425587 x N <CNV:TR>,<CNV:TR>` with `SVLEN=9;RN=0,1;RUS=AGC;RUC=4` and a FASTA as
    /// `AGCAGCAGC/-`, the run deleted) and its second one sequence of four units; a count the list
    /// lacks is no repeats.
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

    /// The allele string Ensembl VEP 116.3 builds over the run `TCATCATCA` of
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
    /// breakend, as in Ensembl VEP 116.3, whose rows are the bare `chromosome_breakpoint`;
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
    /// Ensembl VEP 116.3 writes two `deletion` rows per transcript for
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

    /// Ensembl VEP 116.3 writes `21 30425650 x N N[21:30425700[,<NON_REF> . . .` as the
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

    /// The length an allele states for itself: the unit's length times its RUC, or its
    /// RB where the record carries no RUC, summed over the sequences RN gives it.
    #[test]
    fn a_tandem_repeat_allele_states_its_bases_in_whole_units_of_a_written_unit() {
        let stated = |info: &str, alleles: usize| {
            tandem_repeat_stated_bases(&InfoFields::new(info), alleles)
        };
        assert_eq!(stated("SVLEN=9;RUS=XYZ;RUC=2", 1), [Some(6)]);
        assert_eq!(
            stated("SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=4,2", 2),
            [Some(12), Some(6)]
        );
        assert_eq!(stated("SVLEN=9;RUS=XYZ;RB=12", 1), [Some(12)]);
        assert_eq!(stated("SVLEN=9;RN=2;RUS=XYZ,XY;RUC=3,1", 1), [Some(11)]);
        // An allele of no sequence states nothing, and the next takes the first unit.
        assert_eq!(stated("SVLEN=9;RN=0,1;RUS=XYZ;RUC=4", 2), [None, Some(12)]);
        // A list that ends before an allele's sequences leaves that allele unstated.
        assert_eq!(
            stated("SVLEN=9;RN=1,1;RUS=XYZ;RUC=4,2", 2),
            [Some(12), None]
        );
        // No unit, no count, a `.` unit, a fractional count, an RB that is no whole
        // number of units, an RN that is not one whole number per allele and a
        // confidence interval around the count state nothing.
        for info in [
            "SVLEN=9",
            "SVLEN=9;RUS=XYZ",
            "SVLEN=9;RUC=2",
            "SVLEN=9;RN=1;RUS=.;RUC=4",
            "SVLEN=9;RN=1;RUS=XYZ;RUC=4.9",
            "SVLEN=9;RUS=XYZ;RB=10",
            "SVLEN=9;RN=.;RUS=XYZ;RUC=2",
            "SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=2,2",
            "SVLEN=9;RN=1;RUS=XYZ;RUC=2;CIRUC=-1,2",
            "SVLEN=9;RN=1;RUS=XYZ;RB=6;CIRB=-3,6",
        ] {
            assert_eq!(stated(info, 1), [None], "{info}");
        }
        // A stated length holds up to 64 bits: a count, a product or a sum of sequences
        // past them states nothing, and neither does an RN past the lists.
        assert_eq!(
            stated("SVLEN=9;RUS=XYZ;RUC=6148914691236517205", 1),
            [Some(u64::MAX)]
        );
        for info in [
            "SVLEN=9;RUS=XYZ;RUC=18446744073709551616",
            "SVLEN=9;RUS=XYZ;RUC=6148914691236517206",
            "SVLEN=9;RN=2;RUS=XYZ,XYZ;RUC=6148914691236517205,1",
            "SVLEN=9;RN=18446744073709551615;RUS=XYZ;RUC=4",
        ] {
            assert_eq!(stated(info, 1), [None], "{info}");
        }
    }

    /// An allele is marked a loss where its stated bases are fewer than the run the
    /// largest SVLEN gives; an allele as long as the run, a longer one and one whose
    /// length the record does not state are not, nor is a record of another class.
    #[test]
    fn an_allele_that_states_fewer_bases_than_its_run_is_marked_a_loss() {
        let marks = |info: &str, alts: &str| -> Vec<bool> {
            structural_variants(
                "21",
                100,
                "N",
                &written_alts(alts),
                &InfoFields::new(info),
                None,
            )
            .unwrap()
            .iter()
            .map(|variant| variant.states_a_loss)
            .collect()
        };
        assert_eq!(marks("SVLEN=9;RN=1;RUS=XYZ;RUC=2", "<CNV:TR>"), [true]);
        assert_eq!(marks("SVLEN=9;RN=1;RUS=XYZ;RUC=3", "<CNV:TR>"), [false]);
        assert_eq!(marks("SVLEN=9;RN=1;RUS=XYZ;RUC=4", "<CNV:TR>"), [false]);
        assert_eq!(marks("END=109;RN=1;RUS=XYZ;RUC=2", "<CNV:TR>"), [true]);
        assert_eq!(
            marks("SVLEN=9,6;RN=1,1;RUS=XYZ,XYZ;RUC=4,1", "<CNV:TR>,<CNV:TR>"),
            [false, true]
        );
        assert_eq!(marks("SVLEN=9;RN=1;RUS=.;RUC=4", "<CNV:TR>"), [false]);
        assert_eq!(marks("SVLEN=9", "<CNV:TR>"), [false]);
        // A block written beside the allele leaves the record unmarked.
        assert!(marks("SVLEN=9;RN=1;RUS=XYZ;RUC=2", "<CNV:TR>,<NON_REF>")
            .iter()
            .all(|marked| !marked));
        assert_eq!(marks("SVTYPE=DEL;END=109", "<DEL>"), [false]);
    }

    /// The marks of each allele of the record `21 100 . N <alts>` with INFO `info`: as
    /// long as its run, then a loss.
    fn length_marks(info: &str, alts: &str) -> Vec<(bool, bool)> {
        structural_variants(
            "21",
            100,
            "N",
            &written_alts(alts),
            &InfoFields::new(info),
            None,
        )
        .unwrap()
        .iter()
        .map(|variant| (variant.states_the_run_length, variant.states_a_loss))
        .collect()
    }

    const SAME: (bool, bool) = (true, false);
    const LOSS: (bool, bool) = (false, true);
    const NEITHER: (bool, bool) = (false, false);

    /// An allele is marked as long as its run where its stated bases equal the run
    /// SVLEN or END gives, by RUC or by RB; an allele one base longer or shorter, one
    /// whose length the record does not state, one on the pair of a record with no
    /// run and one of another class are not.
    #[test]
    fn an_allele_that_states_as_many_bases_as_its_run_is_marked_and_no_other() {
        for info in [
            "SVLEN=9;RN=1;RUS=XYZ;RUC=3",
            "END=109;RN=1;RUS=XYZ;RUC=3",
            "SVLEN=9;RUS=XYZ;RUC=3",
            "SVLEN=9;RUS=XYZ;RB=9",
            "SVLEN=9;RUS=X;RUC=9",
            // Two sequences of one allele, three bases and six.
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUC=1,3",
            // An RN of `0` is false to Perl and reads as no RN (Parser/VCF.pm:394-397).
            "SVLEN=9;RN=0;RUS=XYZ;RUC=3",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [SAME], "{info}");
        }
        // One base more than the run and one base fewer, by the count and by the run.
        assert_eq!(length_marks("SVLEN=9;RUS=X;RUC=10", "<CNV:TR>"), [NEITHER]);
        assert_eq!(length_marks("SVLEN=9;RUS=X;RUC=8", "<CNV:TR>"), [LOSS]);
        assert_eq!(length_marks("SVLEN=8;RUS=XYZ;RUC=3", "<CNV:TR>"), [NEITHER]);
        assert_eq!(length_marks("SVLEN=10;RUS=XYZ;RUC=3", "<CNV:TR>"), [LOSS]);
        // No unit, a confidence interval around the count, a count of four and a half
        // two-base units, a `.` unit counted nine times and an RB of nine bases in
        // two-base units state no length.
        for info in [
            "SVLEN=9",
            "SVLEN=9;RUC=3",
            "SVLEN=9;RN=1;RUS=XYZ;RUC=3;CIRUC=-1,1",
            "SVLEN=9;RUS=XYZ;RB=9;CIRB=-3,3",
            "SVLEN=9;RN=1;RUS=XY;RUC=4.5",
            "SVLEN=9;RN=1;RUS=.;RUC=9",
            "SVLEN=9;RUS=XY;RB=9",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
        // Each allele of a record is read against the record's one run, and the first
        // allele of `RN=0,1`, one of no sequence, states no length.
        for (info, expected) in [
            ("SVLEN=9,9;RN=1,1;RUS=XYZ,XYZ;RUC=3,2", [SAME, LOSS]),
            ("SVLEN=9,9;RN=1,1;RUS=XYZ,XYZ;RUC=2,3", [LOSS, SAME]),
            ("SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=4,3", [NEITHER, SAME]),
            ("SVLEN=9;RN=0,1;RUS=XYZ;RUC=3", [NEITHER, SAME]),
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>,<CNV:TR>"), expected, "{info}");
        }
        // An allele of no base on the pair `start = POS + 1, end = POS` has no run.
        for info in [
            "SVLEN=0;RN=2;RUS=XYZ,XYZ;RUC=0,0",
            "END=100;RN=2;RUS=XYZ,XYZ;RUC=0,0",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
        // A block written beside the allele leaves the record unmarked.
        assert!(
            length_marks("SVLEN=9;RN=1;RUS=XYZ;RUC=3", "<CNV:TR>,<NON_REF>")
                .iter()
                .all(|marked| *marked == NEITHER)
        );
        assert_eq!(
            length_marks("SVTYPE=DUP;SVLEN=9;RUS=XYZ;RUC=3", "<DUP>"),
            [NEITHER]
        );
    }

    /// SVLEN holds one reference length for each allele and the run is the largest
    /// (Parser/VCF.pm:505-511), so under unequal values no allele is marked as long as
    /// its run, whichever allele states the run's nine bases; the loss mark stands.
    #[test]
    fn an_allele_of_a_record_with_unequal_svlen_values_is_not_marked_as_long_as_its_run() {
        for (info, expected) in [
            ("SVLEN=9,6;RN=1,1;RUS=XYZ,XYZ;RUC=2,3", [LOSS, NEITHER]),
            ("SVLEN=9,6;RN=1,1;RUS=XYZ,XYZ;RUC=3,2", [NEITHER, LOSS]),
            ("SVLEN=6,9;RN=1,1;RUS=XYZ,XYZ;RUC=4,3", [NEITHER, NEITHER]),
            // Values a fraction apart are two numbers, and `9.0` is the number 9.
            ("SVLEN=9,9.7;RN=1,1;RUS=XYZ,XYZ;RUC=3,3", [NEITHER, NEITHER]),
            ("SVLEN=9.7,9;RN=1,1;RUS=XYZ,XYZ;RUC=3,3", [NEITHER, NEITHER]),
            ("SVLEN=9,9.0;RN=1,1;RUS=XYZ,XYZ;RUC=3,3", [SAME, SAME]),
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>,<CNV:TR>"), expected, "{info}");
        }
        assert_eq!(
            length_marks("SVLEN=9,6;RUS=XYZ;RUC=3", "<CNV:TR>"),
            [NEITHER]
        );
    }

    /// A span of one base names no one-base reference run: VCF 4.4 writes a repeat
    /// the reference lacks with an SVLEN of 1. A run of two bases is marked.
    #[test]
    fn an_allele_of_one_base_on_a_run_of_one_base_is_not_marked_as_long_as_its_run() {
        assert_eq!(length_marks("SVLEN=1;RUS=X;RUC=1", "<CNV:TR>"), [NEITHER]);
        assert_eq!(length_marks("END=101;RUS=X;RUC=1", "<CNV:TR>"), [NEITHER]);
        assert_eq!(length_marks("SVLEN=2;RUS=X;RUC=2", "<CNV:TR>"), [SAME]);
    }

    /// A record with neither SVLEN nor END spans the one base after POS, which is no
    /// stated run.
    #[test]
    fn an_allele_of_a_record_with_neither_svlen_nor_end_is_not_marked_as_long_as_its_run() {
        assert_eq!(length_marks("RUS=X;RUC=1", "<CNV:TR>"), [NEITHER]);
        assert_eq!(length_marks("RN=1;RUS=X;RB=1", "<CNV:TR>"), [NEITHER]);
    }

    /// An END written beside SVLEN names the end of the run SVLEN gives, or the
    /// record states two runs and no allele is marked as long as either.
    #[test]
    fn an_allele_of_a_record_whose_end_disagrees_with_its_svlen_is_not_marked() {
        assert_eq!(
            length_marks("END=109;SVLEN=9;RUS=XYZ;RUC=3", "<CNV:TR>"),
            [SAME]
        );
        for info in [
            "END=112;SVLEN=9;RUS=XYZ;RUC=3",
            "END=109;SVLEN=12;RUS=XYZ;RUC=4",
            "END=abc;SVLEN=9;RUS=XYZ;RUC=3",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
    }

    /// An RB written beside RUC states each sequence's bases a second time: the
    /// allele is marked where every entry is its unit's length times its count, and
    /// not where one differs, is `.` or `0`, or the list ends early.
    #[test]
    fn an_allele_whose_rb_disagrees_with_its_ruc_is_not_marked_as_long_as_its_run() {
        for info in [
            "SVLEN=9;RUS=XYZ;RUC=3;RB=9",
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUC=1,3;RB=3,6",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [SAME], "{info}");
        }
        for info in [
            "SVLEN=9;RUS=XYZ;RUC=3;RB=12",
            "SVLEN=9;RUS=XYZ;RUC=3;RB=.",
            "SVLEN=9;RUS=XYZ;RUC=3;RB=0",
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUC=1,3;RB=3,5",
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUC=1,3;RB=3",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
        assert_eq!(
            length_marks(
                "SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=3,3;RB=9,6",
                "<CNV:TR>,<CNV:TR>"
            ),
            [SAME, NEITHER]
        );
        // Ensembl VEP counts by RB where it reads RUC as false or as an empty list
        // (Parser/VCF.pm:399-401, 414). The record still writes that RUC, a count its
        // RB does not agree with.
        assert_eq!(length_marks("SVLEN=9;RUS=XYZ;RB=9", "<CNV:TR>"), [SAME]);
        for ruc in ["RUC=0", "RUC=", "RUC", "RUC=,"] {
            let info = format!("SVLEN=9;RUS=XYZ;{ruc};RB=9");
            assert_eq!(length_marks(&info, "<CNV:TR>"), [NEITHER], "{info}");
        }
    }

    /// An RUL states each unit's length a second time: the allele is marked where
    /// every entry is the length of its RUS unit, and not where one differs, is `.`
    /// or `0`, or the list ends early.
    #[test]
    fn an_allele_whose_rul_disagrees_with_its_rus_is_not_marked_as_long_as_its_run() {
        for info in [
            "SVLEN=9;RUS=XYZ;RUL=3;RUC=3",
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUL=3,2;RUC=1,3",
            "SVLEN=9;RUS=XYZ;RUL=3;RB=9",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [SAME], "{info}");
        }
        for info in [
            "SVLEN=9;RUS=XYZ;RUL=4;RUC=3",
            "SVLEN=9;RUS=XYZ;RUL=.;RUC=3",
            "SVLEN=9;RUS=XYZ;RUL=0;RUC=3",
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUL=3,3;RUC=1,3",
            "SVLEN=9;RN=2;RUS=XYZ,XY;RUL=3;RUC=1,3",
            "SVLEN=9;RUS=XYZ;RUL=1;RB=9",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
    }

    /// An RUB gives each repeat unit a length of its own, so a record that carries
    /// one states no allele's length by its RUS and RUC alone.
    #[test]
    fn an_allele_of_a_record_with_an_rub_is_not_marked_as_long_as_its_run() {
        for info in [
            "SVLEN=9;RUS=XYZ;RUC=3;RUB=3,3,3",
            "SVLEN=9;RUS=XYZ;RUC=3;RUB=2,3,4",
            "SVLEN=9;RUS=XYZ;RUC=3;RUB=0",
            "SVLEN=9;RUS=XYZ;RUC=3;RUB",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
    }

    /// INFO/CN is an allele's length over its reference's (VCF 4.4, section 5.7), 1
    /// for an allele as long as its run: an allele is marked where its CN entry is the
    /// number 1 and not where it is another number, `.`, `0` or absent from the list.
    #[test]
    fn an_allele_whose_copy_number_is_not_one_is_not_marked_as_long_as_its_run() {
        for info in ["SVLEN=9;RUS=XYZ;RUC=3;CN=1", "SVLEN=9;RUS=XYZ;RUC=3;CN=1.0"] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [SAME], "{info}");
        }
        for info in [
            "SVLEN=9;RUS=XYZ;RUC=3;CN=3",
            "SVLEN=9;RUS=XYZ;RUC=3;CN=0.9666",
            "SVLEN=9;RUS=XYZ;RUC=3;CN=.",
            "SVLEN=9;RUS=XYZ;RUC=3;CN=0",
            "SVLEN=9;RUS=XYZ;RUC=3;CN=1x",
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>"), [NEITHER], "{info}");
        }
        for (info, expected) in [
            ("SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=3,3;CN=1,1", [SAME, SAME]),
            ("SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=3,3;CN=1,2", [SAME, NEITHER]),
            ("SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=3,3;CN=1", [SAME, NEITHER]),
            (
                "SVLEN=9;RN=1,1;RUS=XYZ,XYZ;RUC=3,2;CN=1,0.6667",
                [SAME, LOSS],
            ),
        ] {
            assert_eq!(length_marks(info, "<CNV:TR>,<CNV:TR>"), expected, "{info}");
        }
    }

    /// RB beside RUC, RUL, RUB and CN count as written wherever the INFO column
    /// carries the key: a value of `0`, which Perl reads as false, a `.`, an empty
    /// value and a bare flag each leave an allele unmarked that the record marks
    /// without the field.
    #[test]
    fn a_field_written_as_zero_a_dot_an_empty_value_or_a_bare_flag_leaves_the_allele_unmarked() {
        assert_eq!(length_marks("SVLEN=9;RUS=XYZ;RUC=3", "<CNV:TR>"), [SAME]);
        let mut marked = Vec::new();
        for field in ["RB", "RUL", "RUB", "CN"] {
            for spelling in ["=0", "=.", "=", ""] {
                let info = format!("SVLEN=9;RUS=XYZ;RUC=3;{field}{spelling}");
                if length_marks(&info, "<CNV:TR>") != [NEITHER] {
                    marked.push(info);
                }
            }
        }
        assert_eq!(marked, Vec::<String>::new());
    }
}
