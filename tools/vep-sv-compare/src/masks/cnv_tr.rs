// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tandem-repeat classes: the literal and the symbolic reading of one repeat, the size
//! of that class, the expansion the reference writes nothing for and the allele it drops.

use std::collections::{BTreeSet, HashMap, HashSet};

use crate::cache::CacheAuthority;
use crate::locus::{location_span, strip_chr, CNV_TR_RUST_ALLELE};
use crate::mate_context::TranscriptModel;
use crate::perl::{info_written, perl_numeric_value, perl_true, PerlNumber};
use crate::record::{record_location, span_location};
use crate::rules::gvcf::perl_split;
use crate::strings::{Id, Sides, Strings, Tuple, TupleSet};
use crate::vcf::VcfVariant;

/// The ALT of a tandem-repeat record.
const CNV_TR_ALT: &str = "<CNV:TR>";

/// The term the symbolic reading adds for a gain of copies. Ensembl VEP without a FASTA
/// reads every tandem repeat as a gain (Utils/VariationEffect.pm:246-256), so it adds this
/// term for a record that loses copies or keeps them too; vep-rs adds it unless the record
/// states the allele as long as its run or shorter.
const GAIN_DIRECTION_TERM: &str = "feature_elongation";

/// The term a symbolic reading adds for a loss of copies where the scored run takes the
/// direction from the record's unit count.
const LOSS_DIRECTION_TERM: &str = "feature_truncation";

/// The coding effect the literal reading names for an insertion of whole repeat units
/// inside the CDS, where the symbolic reading names the region.
const GAIN_LITERAL_CODING_TERMS: [&str; 2] = ["inframe_insertion", "frameshift_variant"];

/// The coding effect the literal reading names for a deletion of whole repeat units inside
/// the CDS.
const LOSS_LITERAL_CODING_TERMS: [&str; 2] = ["inframe_deletion", "frameshift_variant"];

/// The coding effects of a literal allele that is not a unit indel at the run's end: a
/// multi-allelic record, a unit that is not the reference's repeat, an `N` unit or an SVLEN
/// that disagrees with the unit count substitutes bases across the run and reads as a
/// coding change of any kind.
const SUBSTITUTION_LITERAL_CODING_TERMS: [&str; 12] = [
    "inframe_insertion",
    "inframe_deletion",
    "frameshift_variant",
    "protein_altering_variant",
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
    "start_lost",
    "start_retained_variant",
    "coding_sequence_variant",
];

/// The coding terms the symbolic reading names inside the CDS: the region, the whole-unit
/// indel of a run whose length the scored engine knows, and the stop codon a deletion of
/// the run's tail removes.
const SYMBOLIC_CODING_TERMS: [&str; 5] = [
    "coding_sequence_variant",
    "inframe_deletion",
    "inframe_insertion",
    "frameshift_variant",
    "stop_lost",
];

/// The terms two symbolic rows of one run may differ in. Without a FASTA Ensembl VEP never
/// expands a `<CNV:TR>`: every record is a symbolic `tandem_repeat` over the run and a
/// copy-number gain, so a loss reads `feature_elongation` with the region terms
/// (Parser/VCF.pm:516-527, Utils/VariationEffect.pm:246-256).
const UNEXPANDED_TERMS: [&str; 7] = [
    "feature_elongation",
    "feature_truncation",
    "coding_sequence_variant",
    "inframe_deletion",
    "frameshift_variant",
    "inframe_insertion",
    "stop_lost",
];

/// The splice term a symbolic reading can carry. Ensembl VEP evaluates the polypyrimidine
/// tract for a structural allele, over its whole span, and every other splice term for a
/// sequence allele alone (Utils/Constants.pm:985 beside :688, :707, :930, :948 and :966;
/// StructuralVariationOverlapAllele.pm:210-214).
const TRACT_TERM: &str = "splice_polypyrimidine_tract_variant";

/// The tract's two ends, in bases from the first base of an intron at its acceptor end: the
/// third base and the seventeenth (BaseTranscriptVariationAllele.pm:155-157, 159-161).
const TRACT_NEAR: i128 = 2;
const TRACT_FAR: i128 = 16;

/// The most bases an intron's end lies past its start where Ensembl VEP reads the intron as
/// a frameshift intron (BaseTranscriptVariation.pm:1032-1035).
const FRAMESHIFT_INTRON_REACH: i128 = 12;

/// Which reading of a tandem repeat a reference tuple is, beside the scored side's symbolic
/// tuple over the run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairKind {
    /// An insertion of whole units between the run's last base and the next.
    Gain,
    /// A deletion of whole units ending at the run's last base.
    Loss,
    /// Any other literal allele at the run.
    Literal,
    /// The reference's own symbolic row, written without a FASTA.
    Unexpanded,
}

/// Whether an Allele is a literal reading's: at least one character, none of them a
/// lowercase letter `a` to `z`, an underscore, a dot, an angle bracket or a square bracket.
/// Ensembl VEP writes a literal allele in upper case (Parser.pm:554), so a unit of digits,
/// of ambiguity codes or of asterisks is literal where a Sequence Ontology class term or a
/// breakend string is not.
fn is_literal_allele(allele: &str) -> bool {
    !allele.is_empty()
        && !allele
            .chars()
            .any(|c| matches!(c, 'a'..='z' | '_' | '.' | '<' | '>' | '[' | ']'))
}

/// A Location read as its chromosome, start and end.
type Span<'a> = (&'a str, i128, i128);

/// Whether a reference tuple, given by its span and its Allele, is another reading of the
/// tandem repeat a scored tuple reads symbolically over the span `run`, and which; `None`
/// when it is no reading of it, a span on another chromosome among them.
///
/// Ensembl VEP expands the repeat units against the reference run and trims the common
/// prefix, so its tuple sits at the run's 3' end: a gain is an insertion of whole units
/// between the run's last base and the next (Location `end-(end+1)`, Allele the inserted
/// bases), a loss a deletion of whole units ending at the run's last base (Location inside
/// the run, Allele `-`). A literal reading is any other literal allele whose span starts at
/// or inside the run, or at the base after it, and ends at or after the run's start; the
/// trim moves a start forward only, so an allele whose span begins before the run is
/// another record's. An unexpanded reading is the reference's own `tandem_repeat` tuple
/// over the run.
fn cnv_tr_pair_kind(perl: Span<'_>, perl_allele: &str, run: Span<'_>) -> Option<PairKind> {
    let (perl_chrom, p_start, p_end) = perl;
    let (run_chrom, run_start, run_end) = run;
    if perl_chrom != run_chrom {
        return None;
    }
    if perl_allele == CNV_TR_RUST_ALLELE {
        return ((p_start, p_end) == (run_start, run_end)).then_some(PairKind::Unexpanded);
    }
    if perl_allele == "-" {
        return (p_end == run_end && p_start >= run_start && p_start <= p_end)
            .then_some(PairKind::Loss);
    }
    if !is_literal_allele(perl_allele) {
        return None;
    }
    // `None` where the base after the run is past the coordinates read.
    let after_run = run_end.checked_add(1);
    if p_start == run_end && Some(p_end) == after_run {
        return Some(PairKind::Gain);
    }
    (run_start <= p_start && after_run.is_none_or(|after| p_start <= after) && p_end >= run_start)
        .then_some(PairKind::Literal)
}

/// The terms of a consequence column, an empty field being none.
fn terms(csq: &str) -> BTreeSet<&str> {
    csq.split(',').filter(|term| !term.is_empty()).collect()
}

/// Whether a paired literal and symbolic reading of one tandem repeat differ only in
/// representation.
///
/// The two term sets must agree once the symbolic reading's direction term is set aside
/// (`feature_elongation` on a gain, one of the two on a loss, either on a literal
/// reading), with one correspondence inside the CDS: where both sides hold terms the other
/// lacks, the literal side's must be coding effects of its kind and the symbolic side's
/// among [`SYMBOLIC_CODING_TERMS`]. An unexpanded pair may differ in the terms of
/// [`UNEXPANDED_TERMS`] alone. Any other difference, a start or splice term one side
/// carries and the other lacks, is a consequence disagreement, the tract term of the
/// symbolic side excepted where `tract_of_the_run` says that the run alone reaches the
/// tract ([`tract_is_of_the_run_alone`]).
///
/// On a loss the term is `feature_truncation` where `loss_named` says that the scored side
/// names a loss over the run on the feature, and `feature_elongation` where it names none:
/// a scored run without a FASTA reads a loss it does not take from the record's counts as
/// a gain, while one that names a loss names a gain for another allele of the run alone.
fn is_cnv_tr_representation_pair(
    perl_csq: &str,
    rust_csq: &str,
    kind: PairKind,
    loss_named: bool,
    tract_of_the_run: bool,
) -> bool {
    let (perl_terms, rust_terms) = (terms(perl_csq), terms(rust_csq));
    let (direction, literal_coding): (&[&str], &[&str]) = match kind {
        PairKind::Unexpanded => {
            return perl_terms
                .symmetric_difference(&rust_terms)
                .all(|term| UNEXPANDED_TERMS.contains(term));
        }
        PairKind::Gain => (&[GAIN_DIRECTION_TERM], &GAIN_LITERAL_CODING_TERMS),
        PairKind::Loss if loss_named => (&[LOSS_DIRECTION_TERM], &LOSS_LITERAL_CODING_TERMS),
        PairKind::Loss => (&[GAIN_DIRECTION_TERM], &LOSS_LITERAL_CODING_TERMS),
        PairKind::Literal => (
            &[GAIN_DIRECTION_TERM, LOSS_DIRECTION_TERM],
            &SUBSTITUTION_LITERAL_CODING_TERMS,
        ),
    };
    let perl_extra: Vec<&str> = perl_terms.difference(&rust_terms).copied().collect();
    let rust_extra: Vec<&str> = rust_terms
        .difference(&perl_terms)
        .copied()
        .filter(|term| !direction.contains(term))
        .filter(|term| !(tract_of_the_run && *term == TRACT_TERM))
        .collect();
    if !rust_extra.is_empty() && !perl_extra.is_empty() {
        return rust_extra
            .iter()
            .all(|term| SYMBOLIC_CODING_TERMS.contains(term))
            && perl_extra.iter().all(|term| literal_coding.contains(term));
    }
    rust_extra.is_empty() && perl_extra.is_empty()
}

/// Whether the span `low` to `high` overlaps the polypyrimidine tract of an intron of the
/// transcript: the fifteen bases that begin two bases inside the intron's acceptor end, its
/// higher end on the forward strand and its lower end on the reverse
/// (BaseTranscriptVariationAllele.pm:155-157, 159-161; Utils/VariationEffect.pm:607-614).
///
/// `None` where the structure does not settle it: a transcript with an intron of thirteen
/// bases or fewer, whose tract Ensembl VEP passes over for a span that overlaps the intron
/// (BaseTranscriptVariationAllele.pm:138-141, 171-174), or a strand that is neither.
fn span_overlaps_a_tract(model: &TranscriptModel, low: i128, high: i128) -> Option<bool> {
    let forward = match model.strand {
        1 => true,
        -1 => false,
        _ => return None,
    };
    let introns = model.introns();
    if introns
        .iter()
        .any(|&(start, end)| end.saturating_sub(start) <= FRAMESHIFT_INTRON_REACH)
    {
        return None;
    }
    Some(introns.iter().any(|&(start, end)| {
        let (tract_low, tract_high) = if forward {
            (
                end.saturating_sub(TRACT_FAR),
                end.saturating_sub(TRACT_NEAR),
            )
        } else {
            (
                start.saturating_add(TRACT_NEAR),
                start.saturating_add(TRACT_FAR),
            )
        };
        high >= tract_low && low <= tract_high
    }))
}

/// Whether the cache shows the tract term to be the symbolic reading's alone on `feature`:
/// the run overlaps a tract of the transcript and no exon of it, and the literal reading's
/// span overlaps no tract ([`span_overlaps_a_tract`]). Ensembl VEP tests the symbolic row on
/// the whole run (StructuralVariationOverlapAllele.pm:210-214), and the bases a literal
/// allele changes lie inside its Location. Nothing is shown without the transcript's
/// structure in the cache.
///
/// Ensembl VEP evaluates the term for a span that overlaps the transcript and no exon of it
/// (Utils/Constants.pm:977-980; BaseVariationFeatureOverlapAllele.pm:317-333, 465-477,
/// 488-517), so the symbolic reading of a run that overlaps an exon carries no tract term.
fn tract_is_of_the_run_alone(
    cache: &CacheAuthority,
    feature: &str,
    run: Span<'_>,
    literal: Span<'_>,
) -> bool {
    let Some(model) = cache.model(strip_chr(run.0), feature) else {
        return false;
    };
    let (run_low, run_high) = (run.1.min(run.2), run.1.max(run.2));
    let reached = |&(start, end): &(i128, i128)| run_high >= start && run_low <= end;
    let evaluated = reached(&(model.start, model.end)) && !model.exons.iter().any(reached);
    let overlaps =
        |(_, start, end): Span<'_>| span_overlaps_a_tract(model, start.min(end), start.max(end));
    evaluated && overlaps(run) == Some(true) && overlaps(literal) == Some(false)
}

/// A reference tuple beside the scored side's symbolic tuple of the same tandem repeat on
/// the same feature.
struct Pair {
    perl: Tuple,
    /// Whether the reference alone holds its tuple.
    perl_only: bool,
    rust: Tuple,
    kind: PairKind,
    /// Whether one of the scored side's symbolic tuples over the run on the feature
    /// carries `feature_truncation`.
    loss_named: bool,
}

/// The scored side's symbolic tuples with their runs, under their chromosome, Feature and
/// Feature type.
type SymbolicByKey<'a> = HashMap<(&'a str, Id, Id), Vec<(Tuple, Span<'a>)>>;

/// Every pair of one tandem repeat on one feature: a scored-only tuple with the symbolic
/// allele and a reference tuple that is a reading of it by [`cnv_tr_pair_kind`], the two on
/// one chromosome with one Feature and Feature type.
///
/// The reference tuples are its one-sided ones and, with `with_shared`, its `tandem_repeat`
/// tuples that both sides hold: without a FASTA the reference writes one symbolic row per
/// allele of a multi-ALT record, identical rows being one tuple, so the tuple the scored
/// side matched for one allele is the unexpanded reading of the record's other alleles.
fn cnv_tr_pairs(sides: &Sides<'_>, with_shared: bool) -> Vec<Pair> {
    let strings = sides.strings;
    let mut pairs = Vec::new();
    let Some(symbolic) = strings.get(CNV_TR_RUST_ALLELE) else {
        return pairs;
    };
    let mut rust_by_key: SymbolicByKey<'_> = HashMap::new();
    let mut rust_features: HashSet<(Id, Id)> = HashSet::new();
    for rust in sides.rust_only.iter().filter(|t| t.allele() == symbolic) {
        let Some(run) = location_span(strings.text(rust.location())) else {
            continue;
        };
        rust_by_key
            .entry((run.0, rust.feature(), rust.feature_type()))
            .or_default()
            .push((rust, run));
        rust_features.insert((rust.feature(), rust.feature_type()));
    }
    if rust_by_key.is_empty() {
        return pairs;
    }
    let shared: Vec<Tuple> = if with_shared {
        sides
            .perl
            .iter()
            .filter(|t| t.allele() == symbolic && sides.rust.contains(*t))
            .collect()
    } else {
        Vec::new()
    };
    let one_sided = sides.perl_only.iter().map(|perl| (perl, true));
    for (perl, perl_only) in one_sided.chain(shared.into_iter().map(|perl| (perl, false))) {
        if !rust_features.contains(&(perl.feature(), perl.feature_type())) {
            continue;
        }
        let Some(span) = location_span(strings.text(perl.location())) else {
            continue;
        };
        let Some(readings) = rust_by_key.get(&(span.0, perl.feature(), perl.feature_type())) else {
            continue;
        };
        let allele = strings.text(perl.allele());
        for &(rust, run) in readings {
            if let Some(kind) = cnv_tr_pair_kind(span, allele, run) {
                let loss_named = readings.iter().any(|&(other, other_run)| {
                    other_run == run
                        && terms(strings.text(other.consequence())).contains(LOSS_DIRECTION_TERM)
                });
                pairs.push(Pair {
                    perl,
                    perl_only,
                    rust,
                    kind,
                    loss_named,
                });
            }
        }
    }
    pairs
}

/// The pairs of a literal and a symbolic reading that differ only in representation:
/// `(excluded_rust, excluded_perl)`.
///
/// A pair is a reference tuple and a scored tuple of one tandem repeat on one feature
/// (`cnv_tr_pairs`); both are excluded when `is_cnv_tr_representation_pair` holds for
/// their consequence sets. A reference tuple both sides hold (the one symbolic row of a
/// multi-ALT record) pairs the record's other scored tuples and stays on its side.
///
/// The cache is read for one term: a `splice_polypyrimidine_tract_variant` the scored tuple
/// alone carries is set aside where the transcript's structure shows the run alone to reach
/// the tract (`tract_is_of_the_run_alone`), and stays charged without a cache.
pub fn filter_cnv_tr_expansion_divergences(
    sides: &Sides<'_>,
    cache: &CacheAuthority,
) -> (TupleSet, TupleSet) {
    let strings = sides.strings;
    let mut excluded_rust = Vec::new();
    let mut excluded_perl = Vec::new();
    for pair in cnv_tr_pairs(sides, true) {
        let perl_csq = strings.text(pair.perl.consequence());
        let rust_csq = strings.text(pair.rust.consequence());
        let tract_of_the_run = terms(rust_csq).contains(TRACT_TERM)
            && !terms(perl_csq).contains(TRACT_TERM)
            && match (
                location_span(strings.text(pair.rust.location())),
                location_span(strings.text(pair.perl.location())),
            ) {
                (Some(run), Some(literal)) => tract_is_of_the_run_alone(
                    cache,
                    strings.text(pair.rust.feature()),
                    run,
                    literal,
                ),
                _ => false,
            };
        if is_cnv_tr_representation_pair(
            perl_csq,
            rust_csq,
            pair.kind,
            pair.loss_named,
            tract_of_the_run,
        ) {
            excluded_rust.push(pair.rust);
            if pair.perl_only {
                excluded_perl.push(pair.perl);
            }
        }
    }
    (
        TupleSet::from_unsorted(excluded_rust),
        TupleSet::from_unsorted(excluded_perl),
    )
}

/// The size of the literal-versus-symbolic class on each side before the mask:
/// `(rust_total, perl_total)`.
///
/// The scored side's count is its one-sided tuples whose Allele is `tandem_repeat`, with or
/// without a partner. A reference tuple has no allele of its own to be told by (a loss is
/// `-`), so the reference's count is its one-sided tuples that pair with one of those
/// (`cnv_tr_pairs`). These are the two populations the mask draws its pairs from, so each
/// excluded count is at most its total. Nothing is excluded here.
pub fn count_cnv_tr_swap_population(sides: &Sides<'_>) -> (u64, u64) {
    let rust_total = sides.strings.get(CNV_TR_RUST_ALLELE).map_or(0, |symbolic| {
        sides
            .rust_only
            .iter()
            .filter(|t| t.allele() == symbolic)
            .count()
    });
    let paired: TupleSet = cnv_tr_pairs(sides, false)
        .iter()
        .map(|pair| pair.perl)
        .collect();
    (rust_total as u64, paired.len() as u64)
}

/// The number of ALTs of a record whose every ALT is `<CNV:TR>`, in any case since the
/// reference reads the column upper-cased, a trailing empty ALT being none; `None` for any
/// other record.
fn cnv_tr_alts(alt: &str) -> Option<usize> {
    let alts = perl_split(alt, ',');
    // No character outside ASCII upper-cases into the characters of `<CNV:TR>`, so the
    // comparison in ASCII is the comparison of the upper-cased column.
    (!alts.is_empty() && alts.iter().all(|alt| alt.eq_ignore_ascii_case(CNV_TR_ALT)))
        .then_some(alts.len())
}

/// The length of a repeat unit in bases: its characters, `.` being one `N`
/// (Parser/VCF.pm:413).
fn unit_length(unit: &str) -> i128 {
    if unit == "." {
        1
    } else {
        unit.chars().count() as i128
    }
}

/// The copies one repeat sequence makes of its unit (Parser/VCF.pm:414): RUC as written, or
/// RB over the unit length, a quotient in floating point. `None` for RB over a unit of no
/// bases, where the reference run ends on a division by zero.
fn unit_copies(unit_len: i128, count: &str, by_ruc: bool) -> Option<PerlNumber> {
    let count = perl_numeric_value(count);
    if by_ruc {
        Some(count)
    } else {
        count.div_int(unit_len)
    }
}

/// The bases one repeat sequence adds to a `<CNV:TR>` allele, `$seq x $num`
/// (Parser/VCF.pm:414, 419): the unit length times its copies ([`unit_copies`]), the copies
/// truncated toward zero as Perl's `x` operator truncates them and never below zero. `None`
/// as for [`unit_copies`].
fn repeat_bases(unit_len: i128, count: &str, by_ruc: bool) -> Option<i128> {
    let copies = unit_copies(unit_len, count, by_ruc)?;
    Some(unit_len.saturating_mul(copies.int().max(0)))
}

/// Whether `_expand_tandem_repeat_allele_string` (Parser/VCF.pm:382-431) builds an allele
/// string the reference writes nothing for, on a `<CNV:TR>` record of one or several ALTs.
///
/// Every part of the string must be absent (an RN of `.`, `:410`), empty or equal to the
/// reference run. A part is the bases of one allele's repeat sequences, each unit length
/// times its count (`repeat_bases`): a count of `0.5` repeats nothing, `3.5` three times,
/// and a count the RUC or RB list has run out of repeats nothing. The run is SVLEN, its
/// largest value (`:505-511`). The split drops a trailing empty part, while an empty part
/// before another is kept as the loss of the run (Parser.pm:881) and is not this shape.
///
/// The repeat sequences take RUS and RUC (or RB) in order across the alleles (`shift`,
/// `:412-414`). An RN Perl reads as false is one repeat sequence per allele, as an absent
/// one is (`:394-397`), and a negative RN none (`:411`). A RUC Perl reads as false falls
/// through to RB (`:385`, `:414`). A unit of no bases repeats nothing under RUC (`:419`);
/// under RB it divides by zero (`:414`), which ends the reference run, so such a record is
/// not judged, as is one whose repeat sequences outrun the RUS list.
pub fn cnv_tr_expansion_is_empty(v: &VcfVariant) -> bool {
    let Some(alts) = cnv_tr_alts(&v.alt) else {
        return false;
    };
    let rus = info_written(&v.info, "RUS");
    let ruc = info_written(&v.info, "RUC");
    let rb = info_written(&v.info, "RB");
    let by_ruc = perl_true(ruc);
    let (Some(rus), Some(counts)) = (rus, if by_ruc { ruc } else { rb }) else {
        return false;
    };
    if !perl_true(Some(rus)) || !perl_true(Some(counts)) {
        return false;
    }
    let rn = info_written(&v.info, "RN");
    let rn_parts = match rn {
        Some(rn) if perl_true(Some(rn)) => perl_split(rn, ','),
        _ => vec!["1"; alts],
    };
    let mut units = rus.split(',');
    let mut counts = counts.split(',');
    let run = info_written(&v.info, "SVLEN").and_then(|svlen| {
        svlen
            .split(',')
            .map(perl_numeric_value)
            .reduce(|largest, value| if value > largest { value } else { largest })
    });
    let mut parts: Vec<i128> = Vec::new();
    for records in rn_parts {
        if records == "." {
            continue;
        }
        let mut bases: i128 = 0;
        for _ in 0..perl_numeric_value(records).int().max(0) {
            let Some(unit) = units.next() else {
                return false;
            };
            let Some(count) = counts.next() else {
                continue;
            };
            let Some(added) = repeat_bases(unit_length(unit), count, by_ruc) else {
                return false;
            };
            bases = bases.saturating_add(added);
        }
        parts.push(bases);
    }
    while parts.last() == Some(&0) {
        parts.pop();
    }
    parts
        .iter()
        .all(|bases| run.is_some_and(|run| run == *bases))
}

/// The spans a tuple set's Locations name, by chromosome: whether one of them lies beside
/// a run is one search, whatever their number.
struct SpanIndex<'a> {
    /// A chromosome's spans in the order of their starts, each as its start and the
    /// largest end among the spans up to it.
    by_chrom: HashMap<&'a str, Vec<(i128, i128)>>,
}

impl<'a> SpanIndex<'a> {
    /// The index over the Locations that name a span and hold a tuple `keep` admits.
    fn new(strings: &'a Strings, tuples: &TupleSet, keep: impl Fn(Tuple) -> bool) -> SpanIndex<'a> {
        let mut by_chrom: HashMap<&'a str, Vec<(i128, i128)>> = HashMap::new();
        for here in tuples
            .as_slice()
            .chunk_by(|a, b| a.location() == b.location())
        {
            if !here.iter().any(|tuple| keep(*tuple)) {
                continue;
            }
            if let Some((chrom, start, end)) = location_span(strings.text(here[0].location())) {
                by_chrom.entry(chrom).or_default().push((start, end));
            }
        }
        for spans in by_chrom.values_mut() {
            spans.sort_unstable();
            let mut reach = i128::MIN;
            for span in spans.iter_mut() {
                reach = reach.max(span.1);
                span.1 = reach;
            }
        }
        SpanIndex { by_chrom }
    }

    /// Whether a span of the index lies on `chrom` within one base of `start`..`end`, as
    /// [`reference_row_beside`](crate::record::reference_row_beside) reads a list of spans.
    fn beside(&self, chrom: &str, start: i128, end: i128) -> bool {
        let Some(spans) = self.by_chrom.get(chrom) else {
            return false;
        };
        let (from, to) = (start.saturating_sub(1), end.saturating_add(1));
        let starting_in_reach = spans.partition_point(|span| span.0 <= to);
        starting_in_reach > 0 && spans[starting_in_reach - 1].1 >= from
    }
}

/// The scored side's `tandem_repeat` tuples at the run of a record whose expansion is empty
/// or unchanged, where the reference wrote nothing beside the run.
///
/// With a FASTA the reference drops a `<CNV:TR>` record whose expansion is empty or
/// unchanged ([`cnv_tr_expansion_is_empty`]) without a warning; without one it writes the
/// symbolic rows at the run. The scored side's tuples at the record's Location under the
/// symbolic allele are excluded when no reference tuple's span lies on the record's
/// chromosome from the base before the run to the base after it.
pub fn filter_cnv_tr_empty_expansion(sides: &Sides<'_>, records: &[VcfVariant]) -> TupleSet {
    let strings = sides.strings;
    let Some(symbolic) = strings.get(CNV_TR_RUST_ALLELE) else {
        return TupleSet::new();
    };
    let mut perl_spans: Option<SpanIndex<'_>> = None;
    let mut excluded = Vec::new();
    for v in records {
        if !cnv_tr_expansion_is_empty(v) {
            continue;
        }
        let (chrom, start, end) = record_location(v);
        let Some(location) = strings.get(&span_location(&chrom, start, end)) else {
            continue;
        };
        let rows = sides.rust.with_prefix(&[location, symbolic]);
        if rows.is_empty() {
            continue;
        }
        let perl_spans =
            perl_spans.get_or_insert_with(|| SpanIndex::new(strings, sides.perl, |_| true));
        if !perl_spans.beside(&chrom, start, end) {
            excluded.extend_from_slice(rows);
        }
    }
    TupleSet::from_unsorted(excluded)
}

/// The repeat fields of a two-ALT `<CNV:TR>,<CNV:TR>` record that both readers of the
/// dropped allele judge.
#[derive(Debug, PartialEq, Eq)]
struct TwoAlleleRn<'a> {
    rus: &'a str,
    ruc: Option<&'a str>,
    rb: Option<&'a str>,
    /// The two RN counts as Perl numifies them for `for (1 .. $records)`.
    counts: [i128; 2],
}

impl<'a> TwoAlleleRn<'a> {
    /// Whether the counts are RUC's: a RUC Perl reads as false falls through to RB.
    fn by_ruc(&self) -> bool {
        perl_true(self.ruc)
    }

    /// The counts of the repeat sequences as written, RUC's or RB's.
    fn counts_written(&self) -> impl Iterator<Item = &'a str> {
        let written = if self.by_ruc() { self.ruc } else { self.rb };
        written.unwrap_or("").split(',')
    }

    /// The first `sequences` repeat units.
    fn units(&self, sequences: i128) -> Vec<&'a str> {
        self.rus.split(',').take(saturated_len(sequences)).collect()
    }
}

/// A count of list items as a length: a count past every list is a length past every list.
fn saturated_len(count: i128) -> usize {
    usize::try_from(count.max(0)).unwrap_or(usize::MAX)
}

/// The repeat fields of a two-ALT `<CNV:TR>,<CNV:TR>` record with a RUS, with a RUC or an
/// RB, and with a two-part RN neither part of which is `.`, each of the three true as Perl
/// reads it; `None` for a record of any other shape.
fn two_allele_cnv_tr_rn(v: &VcfVariant) -> Option<TwoAlleleRn<'_>> {
    if cnv_tr_alts(&v.alt) != Some(2) {
        return None;
    }
    let rus = info_written(&v.info, "RUS");
    let ruc = info_written(&v.info, "RUC");
    let rb = info_written(&v.info, "RB");
    let rn = info_written(&v.info, "RN");
    if !perl_true(rus) || !(perl_true(ruc) || perl_true(rb)) || !perl_true(rn) {
        return None;
    }
    let counts_rn = perl_split(rn?, ',');
    if counts_rn.len() != 2 || counts_rn.contains(&".") {
        return None;
    }
    Some(TwoAlleleRn {
        rus: rus?,
        ruc,
        rb,
        counts: [counts_rn[0], counts_rn[1]].map(|count| perl_numeric_value(count).int()),
    })
}

/// The direction of a record's copy change.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CopyDirection {
    Gain,
    Loss,
    /// No change: the allele the parser builds is as long as the run.
    Same,
}

/// Whether a consequence column holds `term`.
fn names(csq: &str, term: &str) -> bool {
    csq.split(',').any(|held| held == term)
}

/// The scored side's `tandem_repeat` tuples among `rows`, the scored-only tuples at one
/// run, that are the rows of an allele as long as the run beside an allele of no repeat
/// sequence: a tuple without a direction term whose terms with `feature_elongation` are
/// those of another of the tuples on its feature. vep-rs writes the allele of no repeat
/// sequence, which states no length, as a gain, and an allele its record states as long as
/// the run with the same terms less `feature_elongation`.
fn rows_without_the_gain_term(strings: &Strings, rows: &[Tuple], symbolic: Id) -> Vec<Tuple> {
    let of_the_run: Vec<(Tuple, BTreeSet<&str>)> = rows
        .iter()
        .filter(|tuple| tuple.allele() == symbolic)
        .map(|tuple| (*tuple, terms(strings.text(tuple.consequence()))))
        .collect();
    let with_the_gain_term = |tuple: Tuple, own: &BTreeSet<&str>| {
        of_the_run.iter().any(|(other, theirs)| {
            (other.feature(), other.feature_type()) == (tuple.feature(), tuple.feature_type())
                && theirs.len() == own.len() + 1
                && theirs.contains(GAIN_DIRECTION_TERM)
                && own.is_subset(theirs)
        })
    };
    of_the_run
        .iter()
        .filter(|(tuple, own)| {
            !own.contains(GAIN_DIRECTION_TERM)
                && !own.contains(LOSS_DIRECTION_TERM)
                && with_the_gain_term(*tuple, own)
        })
        .map(|(tuple, _)| *tuple)
        .collect()
}

/// The scored side's `tandem_repeat` tuples among `rows`, the scored-only tuples at the run
/// of a two-allele record, that are the rows of its dropped allele.
///
/// A dropped gain or loss carries its direction term. A dropped allele as long as the run
/// carries none where the scored side takes the allele's length from the record
/// (`rows_without_the_gain_term`). Where it wrote no such tuple at the run it read both
/// alleles as a gain, as it reads an allele whose record states no length, and the
/// `feature_elongation` tuples are the dropped allele's as well.
fn dropped_symbolic_rows(
    strings: &Strings,
    rows: &[Tuple],
    symbolic: Id,
    direction: CopyDirection,
) -> Vec<Tuple> {
    let carrying = |term: &str| -> Vec<Tuple> {
        rows.iter()
            .filter(|tuple| tuple.allele() == symbolic)
            .filter(|tuple| names(strings.text(tuple.consequence()), term))
            .copied()
            .collect()
    };
    match direction {
        CopyDirection::Gain => carrying(GAIN_DIRECTION_TERM),
        CopyDirection::Loss => carrying(LOSS_DIRECTION_TERM),
        CopyDirection::Same => {
            let told = rows_without_the_gain_term(strings, rows, symbolic);
            if told.is_empty() {
                carrying(GAIN_DIRECTION_TERM)
            } else {
                told
            }
        }
    }
}

/// [`cnv_tr_rn_zero_dropped_direction`] for a record of the shape, `run` being the bases of
/// its reference run.
fn dropped_direction(shape: &TwoAlleleRn<'_>, run: i128) -> Option<CopyDirection> {
    let empty = match shape.counts {
        [first, second] if first <= 0 && second > 0 => 0,
        [first, second] if first > 0 && second <= 0 => 1,
        _ => return None,
    };
    let sequences = shape.counts[1 - empty];
    let by_ruc = shape.by_ruc();
    let units: Vec<i128> = shape
        .units(sequences)
        .into_iter()
        .map(unit_length)
        .collect();
    let counts: Vec<&str> = shape.counts_written().take(units.len()).collect();
    if (units.len() as i128) < sequences
        || counts.len() < units.len()
        || (units.contains(&0) && !by_ruc)
    {
        return None;
    }
    let mut bases: i128 = 0;
    for (unit_len, count) in units.iter().zip(&counts) {
        bases = bases.saturating_add(repeat_bases(*unit_len, count, by_ruc)?);
    }
    if bases < run {
        return None;
    }
    if empty == 0 {
        return Some(if bases == run {
            CopyDirection::Same
        } else {
            CopyDirection::Gain
        });
    }
    (bases > run).then_some(CopyDirection::Loss)
}

/// The bases of a record's reference run, from its start to its end.
fn run_length(start: i128, end: i128) -> i128 {
    end.saturating_sub(start).saturating_add(1)
}

/// The direction of the allele Ensembl VEP drops from a two-ALT `<CNV:TR>,<CNV:TR>` record
/// of which one allele's RN is a count `for (1 .. $records)` (Parser/VCF.pm:411) runs no
/// iteration for, zero or negative; an RN of `.` is skipped at `:410` and leaves no part.
///
/// That allele's expansion is an empty part of the allele string (`:428`), which
/// Parser.pm:978 does not read as multi-allelic: when the empty part leads, `:881` keeps the
/// first minimised allele alone and the real allele is dropped; when it trails, the split
/// drops the empty part and the allele of no repeat sequence is dropped. The real allele's
/// bases are its repeat sequences' units times their counts, read from the start of RUS and
/// RUC (or RB), the empty allele having taken none.
///
/// `None` when the record is not that shape, or where this rule tells no dropped allele's
/// rows apart: beside a real allele that is a loss, whose direction term a scored run that
/// reads no counts shares with the empty allele's rows. An unchanged real allele is the
/// dropped allele when the empty part leads, and has no direction ([`CopyDirection::Same`]),
/// while with the empty part trailing the reference writes nothing for the record
/// ([`cnv_tr_expansion_is_empty`]).
pub fn cnv_tr_rn_zero_dropped_direction(v: &VcfVariant) -> Option<CopyDirection> {
    let shape = two_allele_cnv_tr_rn(v)?;
    let (_, start, end) = record_location(v);
    dropped_direction(&shape, run_length(start, end))
}

/// [`cnv_tr_rn_zero_dropped_literal`] for a record of the shape.
fn dropped_literal(shape: &TwoAlleleRn<'_>, max_sv_size: i128) -> Option<String> {
    let [empty, sequences] = shape.counts;
    if empty > 0 || sequences <= 0 {
        return None;
    }
    let by_ruc = shape.by_ruc();
    let units: Vec<&str> = shape
        .units(sequences)
        .into_iter()
        .map(|unit| if unit == "." { "N" } else { unit })
        .collect();
    let counts: Vec<&str> = shape.counts_written().take(units.len()).collect();
    if (units.len() as i128) < sequences
        || counts.len() < units.len()
        || (units.contains(&"") && !by_ruc)
    {
        return None;
    }
    let cap = if max_sv_size == 0 { 5000 } else { max_sv_size };
    let mut literal = String::new();
    for (unit, count) in units.iter().zip(&counts) {
        let unit_len = unit.chars().count() as i128;
        let copies = unit_copies(unit_len, count, by_ruc)?;
        if copies * unit_len > cap {
            return None;
        }
        if !unit.is_empty() {
            literal.push_str(&unit.repeat(saturated_len(copies.int())));
        }
    }
    Some(literal.to_uppercase())
}

/// The real allele behind a leading RN of no repeat sequence, as the parser expands it: the
/// allele's repeat sequences, each unit repeated by its count (`$seq x $num`,
/// Parser/VCF.pm:414, 419) with `.` read as `N` (`:413`), joined and upper-cased
/// (Parser.pm:554).
///
/// `None` for a record of another shape, a trailing empty RN included, and for a record the
/// parser leaves symbolic: one with a repeat sequence of more bases, its unit length times
/// its count as written, than `max_sv_size` (`:403`, 5000 when the cap reads false, so -1
/// expands nothing; `:417-418`, `:428`).
pub fn cnv_tr_rn_zero_dropped_literal(v: &VcfVariant, max_sv_size: i128) -> Option<String> {
    dropped_literal(&two_allele_cnv_tr_rn(v)?, max_sv_size)
}

/// The scored side's tuples of the allele the reference dropped from a two-allele record of
/// which one allele repeats nothing.
///
/// The mask applies where the reference expanded the record: a reference tuple at or beside
/// the run under an allele other than `tandem_repeat`, which only a run with a FASTA
/// writes. A scored run without a FASTA writes the dropped allele as the `tandem_repeat`
/// tuples carrying its direction term at the run ([`cnv_tr_rn_zero_dropped_direction`]), and
/// a dropped allele as long as the run as the tuples that lack the `feature_elongation` of
/// the empty allele's tuples (`dropped_symbolic_rows`). A scored run with one writes the
/// real allele dropped behind a leading empty RN as its literal reading over the run
/// ([`cnv_tr_rn_zero_dropped_literal`]), untrimmed because the scored engine reads the
/// record's allele string as multi-allelic. The scored-only tuples of either kind at the
/// record's Location are excluded. `max_sv_size` is the cap the reference ran with, which
/// bounds the repeat sequences it expands.
pub fn filter_cnv_tr_rn_zero_allele_dropped(
    sides: &Sides<'_>,
    records: &[VcfVariant],
    max_sv_size: i128,
) -> TupleSet {
    let strings = sides.strings;
    if sides.rust_only.is_empty() {
        return TupleSet::new();
    }
    let symbolic = strings.get(CNV_TR_RUST_ALLELE);
    let mut expanded: Option<SpanIndex<'_>> = None;
    let mut excluded = Vec::new();
    for v in records {
        let Some(shape) = two_allele_cnv_tr_rn(v) else {
            continue;
        };
        let (chrom, start, end) = record_location(v);
        let Some(location) = strings.get(&span_location(&chrom, start, end)) else {
            continue;
        };
        let rows = sides.rust_only.with_prefix(&[location]);
        if rows.is_empty() {
            continue;
        }
        let expanded = expanded.get_or_insert_with(|| {
            SpanIndex::new(strings, sides.perl, |tuple| {
                Some(tuple.allele()) != symbolic
            })
        });
        if !expanded.beside(&chrom, start, end) {
            continue;
        }
        if let (Some(direction), Some(symbolic)) =
            (dropped_direction(&shape, run_length(start, end)), symbolic)
        {
            excluded.extend(dropped_symbolic_rows(strings, rows, symbolic, direction));
        }
        let literal =
            dropped_literal(&shape, max_sv_size).and_then(|literal| strings.get(&literal));
        excluded.extend(
            rows.iter()
                .copied()
                .filter(|tuple| literal == Some(tuple.allele())),
        );
    }
    TupleSet::from_unsorted(excluded)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::cache::CacheAuthority;
    use crate::masks::{
        apply_masks, MaskContext, MaskResult, MASKS_115_2, MASKS_116_2, MAX_SV_SIZE,
    };
    use crate::strings::OwnedSides;

    /// One tuple as five texts, written as the comparator normalises them.
    type Row = [String; 5];

    /// A Transcript tuple.
    fn t(loc: &str, allele: &str, feature: &str, csq: &str) -> Row {
        [loc, allele, feature, "Transcript", csq].map(String::from)
    }

    fn borrowed(rows: &[Row]) -> Vec<[&str; 5]> {
        rows.iter()
            .map(|row| row.each_ref().map(String::as_str))
            .collect()
    }

    fn owned_sides(perl: &[Row], rust: &[Row]) -> OwnedSides {
        OwnedSides::from_rows(&borrowed(perl), &borrowed(rust))
    }

    fn set(owned: &OwnedSides, rows: &[Row]) -> TupleSet {
        owned.set(&borrowed(rows))
    }

    /// The kind of reading of a reference Location and Allele beside a scored Location; a
    /// Location that is no span is no reading.
    fn pair_kind(perl_loc: &str, perl_allele: &str, rust_loc: &str) -> Option<PairKind> {
        cnv_tr_pair_kind(
            location_span(perl_loc)?,
            perl_allele,
            location_span(rust_loc)?,
        )
    }

    /// The mask and the class size on one pair of sides, with no cache.
    fn masked(perl: &[Row], rust: &[Row]) -> (OwnedSides, TupleSet, TupleSet, (u64, u64)) {
        masked_beside(perl, rust, &CacheAuthority::default())
    }

    /// The mask and the class size on one pair of sides, beside a cache.
    fn masked_beside(
        perl: &[Row],
        rust: &[Row],
        cache: &CacheAuthority,
    ) -> (OwnedSides, TupleSet, TupleSet, (u64, u64)) {
        let owned = owned_sides(perl, rust);
        let ((ex_rust, ex_perl), totals) = {
            let sides = owned.sides();
            (
                filter_cnv_tr_expansion_divergences(&sides, cache),
                count_cnv_tr_swap_population(&sides),
            )
        };
        (owned, ex_rust, ex_perl, totals)
    }

    /// A record of chromosome 21 with the reference base `N`, as the input reader builds it.
    fn rec(pos: i128, alt: &str, info: &str) -> VcfVariant {
        let variant_type = if alt.contains(',') {
            "Multi_allelic"
        } else {
            "Symbolic_CNV"
        };
        VcfVariant {
            chrom: "21".to_string(),
            pos,
            vid: "r".to_string(),
            reference: "N".to_string(),
            alt: alt.to_string(),
            variant_type: variant_type.to_string(),
            location: format!("21:{pos}"),
            info: info.to_string(),
            pos_text: None,
        }
    }

    /// The one-sided tuples every 116.2 mask leaves in place, each side, with the masks'
    /// result: `(perl_open, rust_open, masks)`.
    fn adjusted(
        owned: &OwnedSides,
        records: &[VcfVariant],
        max_sv_size: i128,
    ) -> (TupleSet, TupleSet, MaskResult) {
        let sides = owned.sides();
        let cache = CacheAuthority::default();
        let context = MaskContext {
            reference_release: "116.2",
            cache: &cache,
            scored_engine: "vep-rs",
            input_variants: Some(records),
            max_sv_size,
            perl_paths: &[],
            rust_paths: &[],
            input_by_output: &BTreeMap::new(),
            coding_ends: None,
        };
        let masks = apply_masks(&sides, &context).expect("the 116.2 masks run");
        (
            sides.perl_only.difference(&masks.all_perl()),
            sides.rust_only.difference(&masks.all_rust()),
            masks,
        )
    }

    mod cnv_tr_swap_population {
        use super::*;

        const TR: &str = "tandem_repeat";

        /// The scored side's symbolic tuple over an eighteen-base run ending at `end`.
        fn run(end: i128, feature: &str, rust_csq: &str) -> Row {
            t(&format!("21:{}-{end}", end - 17), TR, feature, rust_csq)
        }

        /// The reference's insertion of whole units after the run.
        fn gain(end: i128, feature: &str, perl_csq: &str) -> Row {
            t(&format!("21:{end}-{}", end + 1), "CACA", feature, perl_csq)
        }

        /// The reference's deletion of four bases ending at the run's last base.
        fn loss(end: i128, feature: &str, perl_csq: &str) -> Row {
            t(&format!("21:{}-{end}", end - 3), "-", feature, perl_csq)
        }

        #[test]
        fn test_pair_kinds() {
            let run = "21:1001-1018";
            let kind = |loc, allele| pair_kind(loc, allele, run);
            assert_eq!(kind("21:1018-1019", "CACA"), Some(PairKind::Gain));
            assert_eq!(kind("21:1015-1018", "-"), Some(PairKind::Loss));
            assert_eq!(kind("21:1010-1011", "CACA"), Some(PairKind::Literal));
            assert_eq!(
                kind("21:1001-1018", "CACACACACACACACACACA"),
                Some(PairKind::Literal)
            );
            assert_eq!(kind("21:1019-1024", "CAC"), Some(PairKind::Literal));
            assert_eq!(kind("21:1001-1018", "NNNN"), Some(PairKind::Literal));
            assert_eq!(kind("21:1001-1018", TR), Some(PairKind::Unexpanded));
            assert_eq!(kind("21:1001-1017", TR), None);
            assert_eq!(kind("21:1015-1017", "-"), None);
            assert_eq!(kind("21:1001-1018", "deletion"), None);
            assert_eq!(kind("21:1020-1021", "CACA"), None);
            assert_eq!(kind("22:1018-1019", "CACA"), None);
        }

        #[test]
        fn test_a_reading_never_starts_before_the_run() {
            let run = "21:27268731-27268748";
            let whole_run = "CACACACACACACACACACA";
            assert_eq!(pair_kind("21:27254789-27276862", "CPX", run), None);
            assert_eq!(pair_kind("21:27268730-27268748", whole_run, run), None);
            assert_eq!(
                pair_kind("21:27268731-27268748", whole_run, run),
                Some(PairKind::Literal)
            );
            let perl = vec![t("21:27254789-27276862", "CPX", "ENST_A", "intron_variant")];
            let rust = vec![t(run, TR, "ENST_A", "feature_elongation,intron_variant")];
            let (_, ex_rust, ex_perl, totals) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
            assert_eq!(totals, (1, 0));
        }

        #[test]
        fn test_representation_pairs_are_masked_and_consequence_disagreements_are_not() {
            use PairKind::{Gain, Literal, Loss, Unexpanded};
            for (perl_csq, rust_csq, kind, representation) in [
                ("intron_variant", "intron_variant", Gain, true),
                (
                    "3_prime_UTR_variant",
                    "3_prime_UTR_variant,feature_elongation",
                    Gain,
                    true,
                ),
                (
                    "3_prime_UTR_variant",
                    "3_prime_UTR_variant,feature_truncation",
                    Loss,
                    true,
                ),
                (
                    "inframe_insertion",
                    "coding_sequence_variant,feature_elongation",
                    Gain,
                    true,
                ),
                (
                    "frameshift_variant",
                    "coding_sequence_variant,feature_truncation",
                    Loss,
                    true,
                ),
                (
                    "intron_variant",
                    "intron_variant,feature_truncation",
                    Gain,
                    false,
                ),
                (
                    "coding_sequence_variant,start_lost,start_retained_variant",
                    "coding_sequence_variant,feature_truncation",
                    Loss,
                    false,
                ),
                (
                    "splice_region_variant,intron_variant",
                    "intron_variant,feature_elongation",
                    Gain,
                    false,
                ),
                ("inframe_insertion", "feature_elongation", Gain, false),
                (
                    "protein_altering_variant",
                    "coding_sequence_variant,feature_elongation",
                    Literal,
                    true,
                ),
                (
                    "missense_variant",
                    "coding_sequence_variant,feature_truncation",
                    Literal,
                    true,
                ),
                (
                    "5_prime_UTR_variant,coding_sequence_variant",
                    "5_prime_UTR_variant,coding_sequence_variant",
                    Literal,
                    true,
                ),
                (
                    "splice_region_variant,missense_variant",
                    "coding_sequence_variant,feature_elongation",
                    Literal,
                    false,
                ),
                (
                    "feature_elongation,coding_sequence_variant",
                    "feature_truncation,inframe_deletion",
                    Unexpanded,
                    true,
                ),
                (
                    "feature_elongation,coding_sequence_variant",
                    "feature_elongation,coding_sequence_variant",
                    Unexpanded,
                    true,
                ),
                (
                    "feature_elongation,coding_sequence_variant",
                    "feature_truncation,intron_variant",
                    Unexpanded,
                    false,
                ),
            ] {
                let loss_named = terms(rust_csq).contains(LOSS_DIRECTION_TERM);
                assert_eq!(
                    is_cnv_tr_representation_pair(perl_csq, rust_csq, kind, loss_named, false),
                    representation,
                    "{perl_csq} against {rust_csq} as {kind:?}"
                );
            }
        }

        /// A scored run without a FASTA that reads no counts reads a record that loses units
        /// as a gain, so its symbolic tuple carries `feature_elongation` beside the
        /// reference's deletion of the units: Ensembl VEP 116.3 writes the record `tr_loss`
        /// of the corpus
        /// `tests/golden/116/GRCh38-vcf-tandem-repeats` as `inframe_deletion` with a FASTA
        /// and as `feature_elongation,coding_sequence_variant` without one. Where the scored
        /// side names a loss over the run, its `feature_elongation` tuple is another
        /// allele's and stays a consequence disagreement beside the deletion.
        #[test]
        fn a_loss_pair_sets_aside_the_direction_term_the_scored_side_reads_the_loss_as() {
            use PairKind::{Gain, Loss};
            let gain_read = "feature_elongation,coding_sequence_variant";
            for (perl_csq, rust_csq, kind, loss_named, representation) in [
                ("inframe_deletion", gain_read, Loss, false, true),
                ("frameshift_variant", gain_read, Loss, false, true),
                (
                    "3_prime_UTR_variant",
                    "3_prime_UTR_variant,feature_elongation",
                    Loss,
                    false,
                    true,
                ),
                (
                    "3_prime_UTR_variant,NMD_transcript_variant",
                    "3_prime_UTR_variant,NMD_transcript_variant,feature_elongation",
                    Loss,
                    false,
                    true,
                ),
                (
                    "non_coding_transcript_exon_variant",
                    "feature_elongation,non_coding_transcript_exon_variant",
                    Loss,
                    false,
                    true,
                ),
                (
                    "inframe_deletion",
                    "feature_truncation,coding_sequence_variant",
                    Loss,
                    true,
                    true,
                ),
                ("inframe_deletion", gain_read, Loss, true, false),
                ("inframe_deletion", "feature_elongation", Loss, false, false),
                ("inframe_insertion", gain_read, Loss, false, false),
                (
                    "coding_sequence_variant,start_lost,start_retained_variant",
                    gain_read,
                    Loss,
                    false,
                    false,
                ),
                (
                    "intron_variant",
                    "intron_variant,splice_polypyrimidine_tract_variant",
                    Loss,
                    false,
                    false,
                ),
                (
                    "inframe_insertion",
                    "feature_truncation,coding_sequence_variant",
                    Gain,
                    true,
                    false,
                ),
            ] {
                assert_eq!(
                    is_cnv_tr_representation_pair(perl_csq, rust_csq, kind, loss_named, false),
                    representation,
                    "{perl_csq} against {rust_csq} as {kind:?}, a loss named: {loss_named}"
                );
            }
        }

        /// A transcript of chromosome 21 with one intron, `201-300`, between two exons.
        fn one_intron(strand: i128) -> TranscriptModel {
            TranscriptModel {
                chrom: "21".to_string(),
                start: 100,
                end: 400,
                strand,
                biotype: "protein_coding".to_string(),
                exons: vec![(100, 200), (301, 400)],
                cds: None,
                mirna_mature: Vec::new(),
            }
        }

        /// A cache that files `model` under chromosome 21 as `ENST_A`.
        fn cache_of(model: TranscriptModel) -> CacheAuthority {
            let mut cache = CacheAuthority::default();
            cache
                .models
                .entry("21".to_string())
                .or_default()
                .insert("ENST_A".to_string(), model);
            cache
        }

        /// The tract is the intron's third to seventeenth base from its acceptor end
        /// (BaseTranscriptVariationAllele.pm:155-157, 159-161): `284-298` of the intron
        /// `201-300` on the forward strand and `203-217` on the reverse.
        #[test]
        fn a_span_overlaps_the_tract_at_the_acceptor_end_of_an_intron() {
            for (strand, low, high, overlaps) in [
                (1, 270, 284, true),
                (1, 270, 283, false),
                (1, 298, 310, true),
                (1, 299, 310, false),
                (1, 203, 217, false),
                (-1, 217, 230, true),
                (-1, 218, 230, false),
                (-1, 190, 203, true),
                (-1, 190, 202, false),
                (-1, 284, 298, false),
            ] {
                assert_eq!(
                    span_overlaps_a_tract(&one_intron(strand), low, high),
                    Some(overlaps),
                    "{low}-{high} on strand {strand}"
                );
            }
        }

        /// Ensembl VEP passes over the tract of an intron of thirteen bases or fewer for a
        /// span that overlaps it (BaseTranscriptVariationAllele.pm:138-141), so a transcript
        /// with such an intron settles nothing, as a strand that is neither settles nothing.
        #[test]
        fn a_short_intron_or_an_unknown_strand_settles_no_tract() {
            let mut short = one_intron(-1);
            short.exons = vec![(100, 200), (214, 400), (501, 600)];
            assert_eq!(span_overlaps_a_tract(&short, 403, 410), None);
            let mut longer = one_intron(-1);
            longer.exons = vec![(100, 200), (215, 400), (501, 600)];
            assert_eq!(span_overlaps_a_tract(&longer, 403, 410), Some(true));
            assert_eq!(span_overlaps_a_tract(&one_intron(0), 203, 217), None);
            let mut one_exon = one_intron(-1);
            one_exon.exons = vec![(100, 400)];
            assert_eq!(span_overlaps_a_tract(&one_exon, 203, 217), Some(false));
        }

        /// The tract term is the one splice term Ensembl VEP evaluates for a structural
        /// allele (Utils/Constants.pm:985), so on the symbolic side it is the run's span;
        /// the predicate sets it aside only where the caller has shown the run alone to
        /// reach the tract.
        #[test]
        fn the_tract_term_of_the_symbolic_side_is_set_aside_where_the_run_alone_reaches_it() {
            use PairKind::{Gain, Literal, Loss, Unexpanded};
            let tract = "intron_variant,splice_polypyrimidine_tract_variant";
            let exon_and_tract =
                "coding_sequence_variant,feature_elongation,splice_polypyrimidine_tract_variant";
            for (perl_csq, rust_csq, kind, tract_of_the_run, representation) in [
                ("intron_variant", tract, Loss, true, true),
                ("intron_variant", tract, Gain, true, true),
                ("intron_variant", tract, Literal, true, true),
                ("intron_variant", tract, Loss, false, false),
                ("intron_variant", tract, Gain, false, false),
                ("intron_variant", tract, Literal, false, false),
                ("intron_variant", tract, Unexpanded, true, false),
                (tract, "intron_variant", Loss, true, false),
                (
                    "intron_variant",
                    "intron_variant,splice_region_variant",
                    Loss,
                    true,
                    false,
                ),
                (
                    "intron_variant",
                    "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
                    Loss,
                    true,
                    false,
                ),
                ("inframe_deletion", exon_and_tract, Loss, true, true),
                ("inframe_deletion", exon_and_tract, Loss, false, false),
                (tract, tract, Loss, false, true),
            ] {
                assert_eq!(
                    is_cnv_tr_representation_pair(
                        perl_csq,
                        rust_csq,
                        kind,
                        false,
                        tract_of_the_run
                    ),
                    representation,
                    "{perl_csq} against {rust_csq} as {kind:?}, the tract shown to be the \
                     run's alone: {tract_of_the_run}"
                );
            }
        }

        /// The filter shows the tract to be the run's alone from the cache: the run
        /// `210-219` overlaps the reverse-strand tract `203-217`, and the deletion of its
        /// last four bases, `216-219`, overlaps it too, while the deletion of its last two
        /// does not.
        #[test]
        fn the_filter_reads_the_cache_for_a_tract_the_run_alone_reaches() {
            let symbolic = t(
                "21:210-219",
                TR,
                "ENST_A",
                "intron_variant,splice_polypyrimidine_tract_variant",
            );
            let past_the_tract = t("21:218-219", "-", "ENST_A", "intron_variant");
            let in_the_tract = t("21:216-219", "-", "ENST_A", "intron_variant");
            let after_the_run = t("21:219-220", "CACA", "ENST_A", "intron_variant");
            let reverse = cache_of(one_intron(-1));

            for literal in [&past_the_tract, &after_the_run] {
                let (perl, rust) = (vec![literal.clone()], vec![symbolic.clone()]);
                let (owned, ex_rust, ex_perl, totals) = masked_beside(&perl, &rust, &reverse);
                assert_eq!(ex_rust, set(&owned, &rust), "{literal:?}");
                assert_eq!(ex_perl, set(&owned, &perl), "{literal:?}");
                assert_eq!(totals, (1, 1));
            }

            let mut short = one_intron(-1);
            short.exons = vec![(50, 80), (90, 200), (301, 400)];
            let mut other_feature = CacheAuthority::default();
            other_feature
                .models
                .entry("21".to_string())
                .or_default()
                .insert("ENST_B".to_string(), one_intron(-1));
            let charged: [(&Row, CacheAuthority); 5] = [
                (&in_the_tract, cache_of(one_intron(-1))),
                (&past_the_tract, CacheAuthority::default()),
                (&past_the_tract, cache_of(one_intron(1))),
                (&past_the_tract, cache_of(short)),
                (&past_the_tract, other_feature),
            ];
            for (literal, cache) in &charged {
                let (perl, rust) = (vec![(*literal).clone()], vec![symbolic.clone()]);
                let (_, ex_rust, ex_perl, totals) = masked_beside(&perl, &rust, cache);
                assert!(ex_rust.is_empty() && ex_perl.is_empty(), "{literal:?}");
                assert_eq!(totals, (1, 1));
            }
        }

        /// Ensembl VEP evaluates the tract term for a span that overlaps the transcript and
        /// no exon of it (Utils/Constants.pm:977-980,
        /// BaseVariationFeatureOverlapAllele.pm:317-333). The term therefore stays charged
        /// on a run that touches an exon, covers the transcript or lies past its end.
        #[test]
        fn the_tract_term_stays_charged_where_the_run_overlaps_an_exon() {
            let reverse = cache_of(one_intron(-1));
            // The tract of the intron `385-398` is `387-401`, one base past the transcript.
            let mut short_last_exon = one_intron(-1);
            short_last_exon.exons = vec![(100, 384), (399, 400)];
            let short_last_exon = cache_of(short_last_exon);
            for (run, deleted, cache, set_aside) in [
                ("21:201-225", "21:220-225", &reverse, true),
                ("21:200-225", "21:220-225", &reverse, false),
                ("21:50-450", "21:440-450", &reverse, false),
                ("21:401-410", "21:405-410", &short_last_exon, false),
            ] {
                let perl = vec![t(deleted, "-", "ENST_A", "intron_variant")];
                let rust = vec![t(
                    run,
                    TR,
                    "ENST_A",
                    "intron_variant,splice_polypyrimidine_tract_variant",
                )];
                let (_, ex_rust, ex_perl, _) = masked_beside(&perl, &rust, cache);
                assert_eq!(
                    (ex_rust.len(), ex_perl.len()),
                    if set_aside { (1, 1) } else { (0, 0) },
                    "the run {run}"
                );
            }
        }

        /// The filter reads whether the scored side names a loss from its symbolic tuples
        /// over the same run on the same feature, the pair's own among them.
        #[test]
        fn the_filter_pairs_a_gain_tuple_with_a_deletion_where_no_scored_tuple_names_the_loss() {
            let deletion = vec![loss(1018, "ENST_A", "inframe_deletion")];
            let as_a_gain = run(1018, "ENST_A", "feature_elongation,coding_sequence_variant");
            let as_a_loss = run(1018, "ENST_A", "feature_truncation,inframe_deletion");

            let rust = vec![as_a_gain.clone()];
            let (owned, ex_rust, ex_perl, totals) = masked(&deletion, &rust);
            assert_eq!(ex_rust, set(&owned, &rust));
            assert_eq!(ex_perl, set(&owned, &deletion));
            assert_eq!(totals, (1, 1));

            let rust = vec![as_a_gain.clone(), as_a_loss.clone()];
            let (owned, ex_rust, ex_perl, totals) = masked(&deletion, &rust);
            assert_eq!(ex_rust, set(&owned, &[as_a_loss]));
            assert_eq!(ex_perl, set(&owned, &deletion));
            assert_eq!(totals, (2, 1));

            let elsewhere = [
                run(1018, "ENST_B", "feature_truncation,inframe_deletion"),
                run(2018, "ENST_A", "feature_truncation,inframe_deletion"),
            ];
            let mut rust = vec![as_a_gain.clone()];
            rust.extend(elsewhere);
            let (owned, ex_rust, ex_perl, totals) = masked(&deletion, &rust);
            assert_eq!(ex_rust, set(&owned, &[as_a_gain]));
            assert_eq!(ex_perl, set(&owned, &deletion));
            assert_eq!(totals, (3, 1));
        }

        #[test]
        fn test_total_counts_every_class_tuple_not_only_the_masked_pairs() {
            let truncated = "coding_sequence_variant,feature_truncation";
            let masked_gain = (
                gain(2000, "ENST_A", "intron_variant"),
                run(2000, "ENST_A", "intron_variant"),
            );
            let masked_loss = (
                loss(3000, "ENST_B", "inframe_deletion"),
                run(3000, "ENST_B", truncated),
            );
            let charged = (
                loss(
                    4000,
                    "ENST_C",
                    "coding_sequence_variant,start_lost,start_retained_variant",
                ),
                run(4000, "ENST_C", truncated),
            );
            let perl = vec![masked_gain.0, masked_loss.0, charged.0.clone()];
            let rust = vec![masked_gain.1, masked_loss.1, charged.1.clone()];
            let (owned, ex_rust, ex_perl, totals) = masked(&perl, &rust);
            assert_eq!(totals, (3, 3));
            assert_eq!((ex_rust.len(), ex_perl.len()), (2, 2));
            assert!(!set(&owned, &[charged.0]).is_subset(&ex_perl));
            assert!(!set(&owned, &[charged.1]).is_subset(&ex_rust));
        }

        #[test]
        fn test_total_is_zero_when_the_corpus_has_no_tandem_repeats() {
            let perl = vec![t(
                "21:1000-2000",
                "deletion",
                "ENST_A",
                "feature_truncation",
            )];
            let rust = vec![t(
                "21:1000-2000",
                "deletion",
                "ENST_B",
                "feature_truncation",
            )];
            let (_, _, _, totals) = masked(&perl, &rust);
            assert_eq!(totals, (0, 0));
        }

        #[test]
        fn test_a_concordant_tandem_repeat_tuple_is_not_a_swap() {
            let shared = vec![t("21:1000-2000", TR, "ENST_A", "feature_elongation")];
            let (_, _, _, totals) = masked(&shared, &shared);
            assert_eq!(totals, (0, 0));
        }

        #[test]
        fn test_a_symbolic_tuple_with_no_literal_partner_counts_on_its_side_only() {
            let perl = vec![t("21:1001-1004", "-", "ENST_A", "intron_variant")];
            let rust = vec![run(2000, "ENST_A", "intron_variant")];
            let (_, ex_rust, ex_perl, totals) = masked(&perl, &rust);
            assert_eq!(totals, (1, 0));
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_excluded_is_a_subset_of_total_on_each_side() {
            let mut perl = Vec::new();
            let mut rust = Vec::new();
            for i in 2..7 {
                let feature = format!("ENST_{i}");
                perl.push(gain(1000 * i, &feature, "intron_variant"));
                rust.push(run(1000 * i, &feature, "intron_variant"));
            }
            perl.push(loss(9000, "ENST_X", "intron_variant,start_lost"));
            rust.push(run(9000, "ENST_X", "intron_variant"));
            let (_, ex_rust, ex_perl, (tot_rust, tot_perl)) = masked(&perl, &rust);
            assert!(ex_rust.len() as u64 <= tot_rust);
            assert!(ex_perl.len() as u64 <= tot_perl);
            assert_eq!(
                (tot_rust, tot_perl, ex_rust.len(), ex_perl.len()),
                (6, 6, 5, 5)
            );
        }

        #[test]
        fn test_the_counter_is_not_a_filter() {
            assert!(MASKS_115_2
                .iter()
                .chain(&MASKS_116_2)
                .all(|mask| mask.starts_with("filter_") && !mask.contains("swap_population")));
        }

        #[test]
        fn a_reading_is_told_by_its_allele_and_its_span() {
            let run = "21:1001-1018";
            let kind = |loc, allele| pair_kind(loc, allele, run);
            for (loc, allele, expected) in [
                ("21:1018-1019", "-", None),
                ("21:1018-1018", "-", Some(PairKind::Loss)),
                ("21:1019-1018", "-", None),
                ("21:1000-1018", "-", None),
                ("21:1019-1019", "CA", Some(PairKind::Literal)),
                ("21:1020-1020", "CA", None),
                ("21:1005-1000", "CA", None),
                ("21:1005-1001", "CA", Some(PairKind::Literal)),
                ("21:1018-1019", TR, None),
                ("chr21:1018-1019", "CA", None),
                ("21:1018", "CA", Some(PairKind::Literal)),
                ("21:1018", "-", Some(PairKind::Loss)),
                ("21:1018-1019", "", None),
                ("21:1001-1018", "ca", None),
            ] {
                assert_eq!(kind(loc, allele), expected, "{loc} {allele:?}");
            }
            let last = i128::MAX;
            let at_the_last_coordinate = format!("21:{}-{last}", last - 8);
            assert_eq!(
                pair_kind(&format!("21:{last}"), "CA", &at_the_last_coordinate),
                Some(PairKind::Literal)
            );
        }

        #[test]
        fn a_literal_allele_has_no_lowercase_letter_underscore_dot_or_bracket() {
            for (allele, literal) in [
                ("", false),
                ("-", true),
                ("ACGT", true),
                ("a", false),
                ("A_B", false),
                ("A.B", false),
                ("<DEL>", false),
                ("[x", false),
                ("x]", false),
                ("A\nC", true),
                ("\u{e9}", true),
                ("\u{c9}", true),
                ("\u{ff5a}", true),
                ("1", true),
                ("*", true),
                (" ", true),
                ("A-C", true),
            ] {
                assert_eq!(is_literal_allele(allele), literal, "{allele:?}");
            }
        }

        #[test]
        fn a_term_set_is_read_without_its_empty_fields() {
            use PairKind::{Gain, Literal, Loss, Unexpanded};
            let both_directions = "a,feature_elongation,feature_truncation";
            for (perl_csq, rust_csq, kind, representation) in [
                ("", "", Gain, true),
                ("a,,b", "b,a", Loss, true),
                ("a", both_directions, Literal, true),
                ("a", both_directions, Gain, false),
                (
                    "stop_lost",
                    "stop_lost,coding_sequence_variant",
                    Gain,
                    false,
                ),
                (
                    "inframe_insertion,x",
                    "coding_sequence_variant,x",
                    Gain,
                    true,
                ),
                ("inframe_deletion", "coding_sequence_variant", Gain, false),
                (
                    "inframe_deletion",
                    "stop_lost,feature_truncation",
                    Loss,
                    true,
                ),
                (" a", "a", Gain, false),
                ("feature_truncation", "feature_elongation", Unexpanded, true),
                ("intron_variant", "", Unexpanded, false),
            ] {
                let loss_named = terms(rust_csq).contains(LOSS_DIRECTION_TERM);
                assert_eq!(
                    is_cnv_tr_representation_pair(perl_csq, rust_csq, kind, loss_named, false),
                    representation,
                    "{perl_csq:?} against {rust_csq:?} as {kind:?}"
                );
            }
        }
    }

    mod tandem_repeat_reading {
        use super::*;

        const RUN_A: &str = "21:30425588-30425596";
        const RUN_B: &str = "21:30425756-30425763";
        const RUN_C: &str = "21:30430410-30430417";
        const TR: &str = "tandem_repeat";
        const UPSTREAM: &str = "upstream_gene_variant";
        const GAIN: &str = "coding_sequence_variant,feature_elongation";
        const LOSS: &str = "feature_truncation,inframe_deletion";

        /// The scored side's two rows over the first run.
        fn rust_a(csq: &str) -> Vec<Row> {
            vec![
                t(RUN_A, TR, "ENST00000334068", UPSTREAM),
                t(RUN_A, TR, "ENST00000390690", csq),
            ]
        }

        /// A literal reading over the whole run beside the symbolic one: the kind is
        /// literal and the four rows pair.
        fn literal_pairs(
            record: VcfVariant,
            perl_allele: &str,
            perl_csq: &str,
            rust_csq: &str,
            run: &str,
        ) {
            let perl = vec![
                t(run, perl_allele, "ENST00000334068", UPSTREAM),
                t(run, perl_allele, "ENST00000390690", perl_csq),
            ];
            let rust = vec![
                t(run, TR, "ENST00000334068", UPSTREAM),
                t(run, TR, "ENST00000390690", rust_csq),
            ];
            assert_eq!(pair_kind(run, perl_allele, run), Some(PairKind::Literal));
            let owned = owned_sides(&perl, &rust);
            let (perl_open, rust_open, masks) = adjusted(&owned, &[record], MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_perl, set(&owned, &perl));
            assert_eq!(masks.cnvtr_rust, set(&owned, &rust));
        }

        /// Whether every one-sided tuple is set aside.
        fn nothing_open(perl: &[Row], rust: &[Row], record: VcfVariant) -> bool {
            let owned = owned_sides(perl, rust);
            let (perl_open, rust_open, _) = adjusted(&owned, &[record], MAX_SV_SIZE);
            perl_open.is_empty() && rust_open.is_empty()
        }

        #[test]
        fn test_a_digit_unit_is_a_literal_reading() {
            literal_pairs(
                rec(
                    30425587,
                    "<CNV:TR>",
                    "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=123;RUC=4",
                ),
                "123123123123",
                "coding_sequence_variant",
                GAIN,
                RUN_A,
            );
        }

        #[test]
        fn test_an_iupac_unit_is_a_literal_reading() {
            literal_pairs(
                rec(
                    30425587,
                    "<CNV:TR>",
                    "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=R;RUC=10",
                ),
                "RRRRRRRRRR",
                "frameshift_variant",
                GAIN,
                RUN_A,
            );
        }

        #[test]
        fn test_a_star_unit_is_a_literal_reading() {
            literal_pairs(
                rec(
                    30425587,
                    "<CNV:TR>",
                    "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=*;RUC=9",
                ),
                "*********",
                "coding_sequence_variant",
                GAIN,
                RUN_A,
            );
        }

        #[test]
        fn test_an_x_unit_is_a_literal_reading() {
            literal_pairs(
                rec(
                    30425587,
                    "<CNV:TR>",
                    "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=X;RUC=10",
                ),
                "XXXXXXXXXX",
                "frameshift_variant",
                GAIN,
                RUN_A,
            );
        }

        #[test]
        fn test_a_loss_read_as_a_protein_altering_substitution_pairs() {
            let record = rec(
                30425755,
                "<CNV:TR>",
                "SVTYPE=CNV;END=30425763;SVLEN=8;RN=3;RUS=AG,AG,C;RUC=1,1,1",
            );
            let trimmed = "21:30425760-30425763";
            let perl = vec![
                t(trimmed, "C", "ENST00000334068", UPSTREAM),
                t(trimmed, "C", "ENST00000390690", "protein_altering_variant"),
            ];
            let rust = vec![
                t(RUN_B, TR, "ENST00000334068", UPSTREAM),
                t(RUN_B, TR, "ENST00000390690", LOSS),
            ];
            assert_eq!(pair_kind(trimmed, "C", RUN_B), Some(PairKind::Literal));
            assert!(nothing_open(&perl, &rust, record));
        }

        #[test]
        fn test_a_loss_read_as_a_stop_gained_substitution_pairs() {
            let record = rec(
                30430409,
                "<CNV:TR>",
                "SVTYPE=CNV;END=30430417;SVLEN=8;RN=3;RUS=CT,AG,C;RUC=1,1,1",
            );
            let trimmed = "21:30430412-30430417";
            let perl = vec![
                t(trimmed, "AGC", "ENST00000334068", "stop_gained"),
                t(trimmed, "AGC", "ENST00000390690", UPSTREAM),
            ];
            let rust = vec![
                t(RUN_C, TR, "ENST00000334068", LOSS),
                t(RUN_C, TR, "ENST00000390690", UPSTREAM),
            ];
            assert!(nothing_open(&perl, &rust, record));
        }

        #[test]
        fn test_a_dot_unit_under_rb_pairs_as_n_bases() {
            literal_pairs(
                rec(
                    30425755,
                    "<CNV:TR>",
                    "SVTYPE=CNV;END=30425763;SVLEN=8;RUS=.;RB=5",
                ),
                "NNNNN",
                "coding_sequence_variant",
                LOSS,
                RUN_B,
            );
        }

        #[test]
        fn test_a_loss_over_the_stop_codon_pairs_with_and_without_a_fasta() {
            let record = rec(30430409, "<CNV:TR>", "SVTYPE=CNV;SVLEN=4000;RUS=CT;RUC=5");
            let run = "21:30430410-30434409";
            let rust = vec![
                t(
                    run,
                    TR,
                    "ENST00000334068",
                    "3_prime_UTR_variant,coding_sequence_variant,feature_truncation,stop_lost",
                ),
                t(run, TR, "ENST00000390690", UPSTREAM),
                t(run, TR, "ENST00000418755", "downstream_gene_variant"),
            ];
            let trimmed = "21:30430418-30434408";
            let fasta = vec![
                t(
                    trimmed,
                    "C",
                    "ENST00000334068",
                    "3_prime_UTR_variant,stop_retained_variant",
                ),
                t(trimmed, "C", "ENST00000390690", UPSTREAM),
                t(trimmed, "C", "ENST00000418755", "downstream_gene_variant"),
            ];
            assert!(nothing_open(&fasta, &rust, record.clone()));
            let symbolic = vec![
                t(
                    run,
                    TR,
                    "ENST00000334068",
                    "3_prime_UTR_variant,coding_sequence_variant",
                ),
                t(run, TR, "ENST00000390690", UPSTREAM),
                t(run, TR, "ENST00000418755", "downstream_gene_variant"),
            ];
            assert!(nothing_open(&symbolic, &rust, record));
        }

        /// The tuples the two engines write for that record: Ensembl VEP's literal ones
        /// from a run with a FASTA beside the symbolic ones vep-rs writes without one, whose
        /// row on ENST00000334068 carries the region term and no stop term.
        #[test]
        fn test_the_literal_stop_retained_tuple_pairs_with_the_symbolic_tuple_of_a_run_without_a_fasta(
        ) {
            let record = rec(30430409, "<CNV:TR>", "SVTYPE=CNV;SVLEN=4000;RUS=CT;RUC=5");
            let run = "21:30430410-30434409";
            let trimmed = "21:30430418-30434408";
            let reference = vec![
                t(
                    trimmed,
                    "C",
                    "ENST00000334068",
                    "3_prime_UTR_variant,stop_retained_variant",
                ),
                t(trimmed, "C", "ENST00000390690", UPSTREAM),
                t(trimmed, "C", "ENST00000418755", "downstream_gene_variant"),
            ];
            let scored = vec![
                t(
                    run,
                    TR,
                    "ENST00000334068",
                    "3_prime_UTR_variant,coding_sequence_variant",
                ),
                t(run, TR, "ENST00000390690", UPSTREAM),
                t(run, TR, "ENST00000418755", "downstream_gene_variant"),
            ];
            assert_eq!(pair_kind(trimmed, "C", run), Some(PairKind::Literal));
            let owned = owned_sides(&reference, &scored);
            let (reference_open, scored_open, masks) = adjusted(&owned, &[record], MAX_SV_SIZE);
            assert!(reference_open.is_empty() && scored_open.is_empty());
            assert_eq!(masks.cnvtr_perl, set(&owned, &reference));
            assert_eq!(masks.cnvtr_rust, set(&owned, &scored));
        }

        #[test]
        fn test_the_collapsed_symbolic_row_of_a_multi_alt_record_pairs_every_vep_rs_allele() {
            let record = rec(
                30425587,
                "<CNV:TR>,<CNV:TR>,<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,1,1;RUS=AGC,AGC,AGC;RUC=4,2,3",
            );
            let perl = rust_a(GAIN);
            // The second allele states a loss, and the third as many bases as the run.
            let loss = t(RUN_A, TR, "ENST00000390690", LOSS);
            let same = t(RUN_A, TR, "ENST00000390690", "coding_sequence_variant");
            let mut rust = perl.clone();
            rust.extend([loss.clone(), same.clone()]);
            let (owned, ex_rust, ex_perl, totals) = masked(&perl, &rust);
            let the_other_alleles = set(&owned, &[loss, same]);
            assert_eq!(ex_rust, the_other_alleles);
            assert!(ex_perl.is_empty());
            assert_eq!(totals, (2, 0));
            let (perl_open, rust_open, masks) = adjusted(&owned, &[record], MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_rust, the_other_alleles);
            assert!(masks.cnvtr_perl.is_empty());
            let two = rec(
                30425587,
                "<CNV:TR>,<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC,AGC;RUC=4,2",
            );
            assert!(nothing_open(&perl, &rust, two));
        }

        /// A multi-ALT record read with a FASTA is one literal tuple for each allele that
        /// differs from the run, written whole over the run, and one symbolic tuple without
        /// one. The pairing starts from a scored tuple whose Allele is `tandem_repeat`: the
        /// reference's literal tuples pair with the scored side's symbolic one, and with the
        /// sides exchanged, the reference symbolic and the scored side literal, nothing
        /// pairs and every tuple stays open.
        #[test]
        fn test_a_multi_alt_record_pairs_only_where_the_scored_side_reads_it_symbolically() {
            let record = rec(
                30425587,
                "<CNV:TR>,<CNV:TR>,<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,1,1;RUS=AGC,AGC,AGC;RUC=4,2,3",
            );
            let symbolic = rust_a(GAIN);
            let literal = vec![
                t(RUN_A, "AGCAGCAGCAGC", "ENST00000334068", UPSTREAM),
                t(RUN_A, "AGCAGC", "ENST00000334068", UPSTREAM),
                t(
                    RUN_A,
                    "AGCAGCAGCAGC",
                    "ENST00000390690",
                    "protein_altering_variant",
                ),
                t(RUN_A, "AGCAGC", "ENST00000390690", "inframe_deletion"),
            ];

            let owned = owned_sides(&literal, &symbolic);
            let (perl_open, rust_open, masks) =
                adjusted(&owned, std::slice::from_ref(&record), MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_perl, set(&owned, &literal));
            assert_eq!(masks.cnvtr_rust, set(&owned, &symbolic));

            let (owned, ex_rust, ex_perl, totals) = masked(&symbolic, &literal);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
            assert_eq!(totals, (0, 0));
            let (perl_open, rust_open, masks) =
                adjusted(&owned, std::slice::from_ref(&record), MAX_SV_SIZE);
            assert_eq!(perl_open, set(&owned, &symbolic));
            assert_eq!(rust_open, set(&owned, &literal));
            assert!(masks.cnvtr_perl.is_empty() && masks.cnvtr_rust.is_empty());
        }

        #[test]
        fn test_an_unchanged_expansion_over_three_repeat_sequences_is_empty() {
            let record = rec(
                30425587,
                "<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1",
            );
            assert!(cnv_tr_expansion_is_empty(&record));
            assert!(nothing_open(&[], &rust_a(GAIN), record));
        }

        #[test]
        fn test_an_unchanged_expansion_whose_second_sequence_repeats_nothing_is_empty() {
            let record = rec(
                30425587,
                "<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RN=2;RUS=AGC,AG;RUC=3,0",
            );
            assert!(cnv_tr_expansion_is_empty(&record));
            assert!(nothing_open(&[], &rust_a(GAIN), record));
        }

        #[test]
        fn test_a_unit_of_no_bases_under_ruc_is_empty() {
            let record = rec(
                30425587,
                "<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AGC;RUC=4,4",
            );
            assert!(cnv_tr_expansion_is_empty(&record));
            assert!(nothing_open(&[], &rust_a(LOSS), record));
        }

        #[test]
        fn test_a_fullwidth_count_repeats_nothing() {
            let record = rec(
                30425587,
                "<CNV:TR>",
                "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=\u{ff14}",
            );
            assert!(perl_numeric_value("\u{ff14}") == 0);
            assert!(cnv_tr_expansion_is_empty(&record));
            assert!(nothing_open(&[], &rust_a(LOSS), record));
        }

        #[test]
        fn test_an_svlen_numified_to_thirty_pairs_as_a_literal_reading() {
            let record = rec(
                30425587,
                "<CNV:TR>",
                "SVTYPE=CNV;SVLEN=30_425_949;RUS=AGC;RUC=4",
            );
            assert_eq!(
                record_location(&record),
                ("21".to_string(), 30425588, 30425617)
            );
            let trimmed = "21:30425597-30425617";
            let perl = vec![
                t(trimmed, "AGC", "ENST00000334068", UPSTREAM),
                t(
                    trimmed,
                    "AGC",
                    "ENST00000390690",
                    "protein_altering_variant",
                ),
            ];
            let run = "21:30425588-30425617";
            let rust = vec![
                t(run, TR, "ENST00000334068", UPSTREAM),
                t(run, TR, "ENST00000390690", LOSS),
            ];
            assert!(nothing_open(&perl, &rust, record));
        }

        #[test]
        fn test_a_two_alt_record_whose_parts_are_the_run_and_nothing_is_empty() {
            let two = "<CNV:TR>,<CNV:TR>";
            let record = rec(
                30425755,
                two,
                "SVTYPE=CNV;END=30425763;SVLEN=8;RN=2,0;RUS=AG,AG;RUC=2,2",
            );
            let rust = vec![
                t(RUN_B, TR, "ENST00000334068", UPSTREAM),
                t(RUN_B, TR, "ENST00000390690", GAIN),
                t(
                    RUN_B,
                    TR,
                    "ENST00000390690",
                    "feature_truncation,frameshift_variant",
                ),
            ];
            assert!(cnv_tr_expansion_is_empty(&record));
            assert!(cnv_tr_expansion_is_empty(&rec(
                30430409,
                two,
                "SVTYPE=CNV;END=30430417;SVLEN=8;RN=2,0;RUS=CT,CT;RUC=2,2"
            )));
            assert!(cnv_tr_expansion_is_empty(&rec(
                30425587,
                two,
                "SVTYPE=CNV;END=30425596;SVLEN=9;RN=.,1;RUS=AGC;RUC=3"
            )));
            assert!(!cnv_tr_expansion_is_empty(&rec(
                30425755,
                two,
                "SVTYPE=CNV;END=30425763;SVLEN=8;RN=0,2;RUS=AG,AG;RUC=2,2"
            )));
            assert!(!cnv_tr_expansion_is_empty(&rec(
                30425587,
                two,
                "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC,AGC;RUC=4"
            )));
            let owned = owned_sides(&[], &rust);
            let (perl_open, rust_open, masks) = adjusted(&owned, &[record], MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_empty_rust, set(&owned, &rust));
        }
    }

    /// The empty-expansion mask on rows of texts.
    fn empty_expansion_rows(
        perl: &[Row],
        rust: &[Row],
        records: &[VcfVariant],
    ) -> (OwnedSides, TupleSet) {
        let owned = owned_sides(perl, rust);
        let excluded = filter_cnv_tr_empty_expansion(&owned.sides(), records);
        (owned, excluded)
    }

    /// The dropped-allele mask on rows of texts.
    fn dropped_allele_rows(
        perl: &[Row],
        rust: &[Row],
        records: &[VcfVariant],
        max_sv_size: i128,
    ) -> (OwnedSides, TupleSet) {
        let owned = owned_sides(perl, rust);
        let excluded = filter_cnv_tr_rn_zero_allele_dropped(&owned.sides(), records, max_sv_size);
        (owned, excluded)
    }

    /// The rows of both lists.
    fn with(rows: &[Row], more: &[Row]) -> Vec<Row> {
        rows.iter().chain(more).cloned().collect()
    }

    mod fasta_named_slice_intergenic {
        use super::*;

        fn single(info: &str) -> VcfVariant {
            rec(30425587, "<CNV:TR>", info)
        }

        #[test]
        fn test_an_empty_tandem_repeat_expansion_is_set_aside_only_when_the_reference_wrote_nothing(
        ) {
            let run = "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;";
            let on_run = |rest: &str| single(&format!("{run}{rest}"));
            let empties = [
                on_run("RUC=3"),
                on_run("RN=.;RUC=4"),
                on_run("RUC=abc"),
                on_run("RUC=-2"),
            ];
            let gain = on_run("RUC=4");
            for v in &empties {
                assert!(cnv_tr_expansion_is_empty(v), "{}", v.info);
            }
            assert!(!cnv_tr_expansion_is_empty(&gain));
            assert!(!cnv_tr_expansion_is_empty(&rec(
                30425587,
                "<CNV:TR>,<CNV:TR>",
                &format!("{run}RUC=3,3")
            )));
            assert!(!cnv_tr_expansion_is_empty(&on_run("RUC=0;RB=12")));
            assert!(cnv_tr_expansion_is_empty(&on_run("RUC=0;RB=9")));
            let symbolic = vec![t(
                "21:30425588-30425596",
                "tandem_repeat",
                "ENST_1",
                "coding_sequence_variant,feature_elongation",
            )];
            let (owned, excluded) = empty_expansion_rows(&[], &symbolic, &empties[..1]);
            assert_eq!(excluded, set(&owned, &symbolic));
            let (_, excluded) = empty_expansion_rows(&symbolic, &symbolic, &empties[..1]);
            assert!(excluded.is_empty());
            let literal_reference = vec![t(
                "21:30425596-30425597",
                "AGC",
                "ENST_1",
                "inframe_insertion",
            )];
            let (_, excluded) = empty_expansion_rows(&literal_reference, &symbolic, &[gain]);
            assert!(excluded.is_empty());
        }

        #[test]
        fn test_an_expansion_count_is_truncated_as_perls_repetition_truncates_it() {
            let run = "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;";
            let on_run = |rest: &str| format!("{run}{rest}");
            let other = |rest: &str| format!("SVTYPE=CNV;END=30425596;{rest}");
            for (info, empty) in [
                (on_run("RUC=3.5"), true),
                (on_run("RUC=3.9"), true),
                (on_run("RUC=0.5"), true),
                (on_run("RUC=4.5"), false),
                (on_run("RUC=2.5"), false),
                (other("SVLEN=9;RUS=.;RUC=9"), true),
                (other("SVLEN=9;RUS=.;RUC=10"), false),
                (other("SVLEN=9;RN=2;RUS=AGC,T;RUC=3,2"), false),
                (other("SVLEN=9;RN=2;RUS=AGC,T;RUC=0,0.9"), true),
                (other("SVLEN=9;RN=2;RUS=AGC;RUC=3,2"), false),
                (other("SVLEN=9;RUS=,AC;RB=12"), false),
                (other("SVLEN=9;RUS=,AC;RUC=0;RB=12"), false),
                (other("SVLEN=9;RUS=AGC;RUC=0;RB=10.5"), true),
                (other("SVLEN=9,12;RUS=AGC;RUC=4"), true),
                (other("SVLEN=9;RN=0;RUS=AGC;RUC=4"), false),
                (other("SVLEN=9;RN=;RUS=AGC;RUC=3"), true),
                (other("SVLEN=9;RUS=0;RUC=3"), false),
                (other("SVLEN=9;RN=-1;RUS=AGC,T;RUC=3,2"), true),
            ] {
                assert_eq!(cnv_tr_expansion_is_empty(&single(&info)), empty, "{info}");
            }
        }

        #[test]
        fn test_a_repeat_unit_of_no_bases_is_empty_under_ruc_and_unjudged_under_rb() {
            let run = "SVTYPE=CNV;END=30425596;SVLEN=9;";
            for (rest, empty) in [
                ("RUS=,AGC;RUC=5", true),
                ("RUS=,AC;RB=12", false),
                ("RUS=,AGC;RN=2;RUC=5,2", false),
                ("RUS=,AGC;RN=2;RUC=5,0", true),
                ("RUS=,AC;RUC=0;RB=12", false),
            ] {
                let info = format!("{run}{rest}");
                assert_eq!(cnv_tr_expansion_is_empty(&single(&info)), empty, "{info}");
            }
        }

        #[test]
        fn an_expansion_is_read_from_the_fields_as_the_parser_reads_them() {
            let run = "SVTYPE=CNV;END=30425596;SVLEN=9;";
            for (alt, info, empty) in [
                ("<cnv:tr>", format!("{run}RUS=AGC;RUC=3"), true),
                ("<CNV:TR>,", format!("{run}RUS=AGC;RUC=3"), true),
                ("<CNV:TR>,,", format!("{run}RUS=AGC;RUC=3"), true),
                (",<CNV:TR>", format!("{run}RUS=AGC;RUC=3"), false),
                ("", format!("{run}RUS=AGC;RUC=3"), false),
                ("<CNV:TR>", format!("{run}RUS=AGC;RUC=3;RN=,"), true),
                ("<CNV:TR>", format!("{run}RUS=AGC;RUC=3;RN=1e30"), false),
                ("<CNV:TR>", format!("{run}RUS=AGC;RUC=1e300"), false),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=1e30;RUS=A;RUC=1e30".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=1e400;RUS=A;RUC=1e400".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=3e400;RUS=AGC;RUC=1e400".to_string(),
                    false,
                ),
                ("<CNV:TR>", "SVTYPE=CNV;RUS=AGC;RUC=3".to_string(), false),
                ("<CNV:TR>", "SVTYPE=CNV;RUS=AGC;RUC=0.5".to_string(), true),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=,9;RUS=AGC;RUC=3".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=-9,-12;RUS=AGC;RUC=3".to_string(),
                    false,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=\u{e9};RUC=9".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RUC=3;RUC=4".to_string(),
                    false,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RUC=3=4".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RB=9;RUC".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RB=9;RUC=".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RB=-9".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC,;RUC=3,7".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RUC= 3".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RUC=+3".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RUC=3abc".to_string(),
                    true,
                ),
                (
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=9;RUS=AGC;RUC=0x3".to_string(),
                    true,
                ),
            ] {
                assert_eq!(
                    cnv_tr_expansion_is_empty(&rec(30425587, alt, &info)),
                    empty,
                    "{alt:?} {info}"
                );
            }
        }

        #[test]
        fn a_count_past_the_word_is_held_at_its_bound() {
            let past_the_word = rec(
                30425587,
                "<CNV:TR>",
                "SVTYPE=CNV;SVLEN=1e300;RUS=A;RUC=1e300",
            );
            assert!(!cnv_tr_expansion_is_empty(&past_the_word));
        }

        /// The spans of the reference tuples `keep` admits, one a tuple.
        fn spans_of<'a>(
            owned: &'a OwnedSides,
            keep: &dyn Fn(Tuple) -> bool,
        ) -> Vec<(&'a str, i128, i128)> {
            owned
                .perl
                .iter()
                .filter(|tuple| keep(*tuple))
                .filter_map(|tuple| location_span(owned.strings.text(tuple.location())))
                .collect()
        }

        #[test]
        fn the_span_index_reads_a_run_as_the_list_of_spans_does() {
            use crate::record::reference_row_beside;
            let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
            let mut next = |below: u64| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (state >> 33) % below
            };
            let mut rows = Vec::new();
            for _ in 0..300 {
                let chrom = ["21", "22", "X"][next(3) as usize];
                let start = next(2000) as i128;
                let end = start + next(40) as i128 - 5;
                let allele = ["tandem_repeat", "-", "CA"][next(3) as usize];
                rows.push(t(&format!("{chrom}:{start}-{end}"), allele, "ENST_A", "x"));
            }
            rows.push(t("nowhere", "CA", "ENST_A", "x"));
            let owned = owned_sides(&rows, &[]);
            let symbolic = owned.strings.get("tandem_repeat");
            let literal_only = |tuple: Tuple| Some(tuple.allele()) != symbolic;
            let all = SpanIndex::new(&owned.strings, &owned.perl, |_| true);
            let literal = SpanIndex::new(&owned.strings, &owned.perl, literal_only);
            let all_spans = spans_of(&owned, &|_| true);
            let literal_spans = spans_of(&owned, &literal_only);
            for _ in 0..4000 {
                let chrom = ["21", "22", "X", "7"][next(4) as usize];
                let start = next(2100) as i128 - 50;
                let end = start + next(30) as i128 - 3;
                assert_eq!(
                    all.beside(chrom, start, end),
                    reference_row_beside(&all_spans, chrom, start, end),
                    "{chrom}:{start}-{end}"
                );
                assert_eq!(
                    literal.beside(chrom, start, end),
                    reference_row_beside(&literal_spans, chrom, start, end),
                    "{chrom}:{start}-{end} beside a literal allele"
                );
            }
            assert!(!all.beside("21", i128::MAX, i128::MAX));
            assert!(all.beside("21", i128::MIN, i128::MAX));
        }
    }

    mod rn_zero_allele_dropped {
        use super::*;

        const RUN: &str = "21:30425588-30425596";
        const INFO: &str = "SVTYPE=CNV;END=30425596;SVLEN=9;";
        const TR: &str = "tandem_repeat";
        const UPSTREAM: &str = "upstream_gene_variant";
        const GAIN: &str = "coding_sequence_variant,feature_elongation";
        const LOSS: &str = "feature_truncation,inframe_deletion";
        const ALTERING: &str = "protein_altering_variant";
        const FIRST: &str = "ENST00000334068";
        const SECOND: &str = "ENST00000390690";

        /// A two-ALT record over the run.
        fn two(info: &str) -> VcfVariant {
            rec(30425587, "<CNV:TR>,<CNV:TR>", &format!("{INFO}{info}"))
        }

        /// The scored side's rows without a FASTA: the upstream row both alleles share, the
        /// loss and the gain.
        fn rust() -> Vec<Row> {
            vec![
                t(RUN, TR, FIRST, UPSTREAM),
                t(RUN, TR, SECOND, LOSS),
                t(RUN, TR, SECOND, GAIN),
            ]
        }

        /// The reference's rows for the allele that repeats nothing.
        fn deletion_rows() -> Vec<Row> {
            vec![
                t(RUN, "-", FIRST, UPSTREAM),
                t(RUN, "-", SECOND, "inframe_deletion"),
            ]
        }

        /// The reference's rows for a gain of one unit.
        fn insertion_rows() -> Vec<Row> {
            let after_run = "21:30425596-30425597";
            vec![
                t(after_run, "AGC", FIRST, UPSTREAM),
                t(after_run, "AGC", SECOND, "inframe_insertion"),
            ]
        }

        /// The dropped-allele mask and what every 116.2 mask leaves open, on one record.
        fn set_aside(perl: &[Row], rust: &[Row], record: &VcfVariant, expected: &[Row]) {
            let owned = owned_sides(perl, rust);
            let records = std::slice::from_ref(record);
            let expected = set(&owned, expected);
            assert_eq!(
                filter_cnv_tr_rn_zero_allele_dropped(&owned.sides(), records, MAX_SV_SIZE),
                expected
            );
            let (perl_open, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.rn_zero_rust, expected);
        }

        #[test]
        fn test_a_leading_rn_zero_drops_the_real_allele() {
            let record = two("RN=0,1;RUS=AGC;RUC=4");
            let gain = t(RUN, TR, SECOND, GAIN);
            assert_eq!(
                cnv_tr_rn_zero_dropped_direction(&record),
                Some(CopyDirection::Gain)
            );
            set_aside(&deletion_rows(), &rust(), &record, &[gain]);
        }

        /// The real allele is the run unchanged. A scored side that takes the allele's length
        /// from the record writes it without a direction term beside the gain's row of the
        /// allele of no repeat sequence, and that row is the dropped allele's. A scored side
        /// that takes no length writes one gain's row for both alleles, which is set aside.
        #[test]
        fn test_a_leading_rn_zero_beside_an_unchanged_real_allele_drops_the_allele_of_no_direction()
        {
            let record = two("RN=0,1;RUS=AGC;RUC=3");
            let gain = t(RUN, TR, SECOND, GAIN);
            let same = t(RUN, TR, SECOND, "coding_sequence_variant");
            assert!(!cnv_tr_expansion_is_empty(&record));
            assert_eq!(
                cnv_tr_rn_zero_dropped_direction(&record),
                Some(CopyDirection::Same)
            );
            let told = vec![t(RUN, TR, FIRST, UPSTREAM), gain.clone(), same.clone()];
            set_aside(&deletion_rows(), &told, &record, &[same]);
            set_aside(
                &deletion_rows(),
                &rust(),
                &record,
                std::slice::from_ref(&gain),
            );
            // A row without a direction term whose other terms are not the gain's row's is
            // no row of the dropped allele: the gain's row is set aside and the row stays.
            let other = t(RUN, TR, SECOND, "intron_variant");
            let rows = vec![gain.clone(), other.clone()];
            let owned = owned_sides(&deletion_rows(), &rows);
            let records = std::slice::from_ref(&record);
            let (_, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
            assert_eq!(masks.rn_zero_rust, set(&owned, &[gain]));
            assert_eq!(rust_open, set(&owned, &[other]));
        }

        #[test]
        fn test_a_trailing_rn_zero_drops_the_rn_zero_allele() {
            let loss = t(RUN, TR, SECOND, LOSS);
            for info in ["RN=1,0;RUS=AGC;RUC=4", "RN=2,0;RUS=AGC,AGC;RUC=2,2"] {
                let record = two(info);
                assert_eq!(
                    cnv_tr_rn_zero_dropped_direction(&record),
                    Some(CopyDirection::Loss),
                    "{info}"
                );
                let owned = owned_sides(&insertion_rows(), &rust());
                let (perl_open, rust_open, masks) = adjusted(&owned, &[record], MAX_SV_SIZE);
                assert!(perl_open.is_empty() && rust_open.is_empty(), "{info}");
                assert_eq!(
                    masks.rn_zero_rust,
                    set(&owned, std::slice::from_ref(&loss)),
                    "{info}"
                );
            }
        }

        #[test]
        fn test_a_vep_rs_run_with_a_fasta_writes_the_dropped_allele_literally() {
            let record = two("RN=0,1;RUS=AGC;RUC=4");
            let perl = deletion_rows();
            let whole = "AGCAGCAGCAGC";
            let literal = vec![
                t(RUN, whole, FIRST, UPSTREAM),
                t(RUN, whole, SECOND, ALTERING),
            ];
            let rust = with(&perl, &literal);
            set_aside(&perl, &rust, &record, &literal);
            let other = t(RUN, "AGCAGCAGCAGCAGC", SECOND, ALTERING);
            let elsewhere = t("21:30425596-30425597", whole, SECOND, "inframe_insertion");
            let (owned, excluded) = dropped_allele_rows(
                &perl,
                &with(&rust, &[other, elsewhere]),
                &[record],
                MAX_SV_SIZE,
            );
            assert_eq!(excluded, set(&owned, &literal));
            let unchanged = two("RN=0,1;RUS=AGC;RUC=3");
            let (_, excluded) =
                dropped_allele_rows(&perl, &perl, std::slice::from_ref(&unchanged), MAX_SV_SIZE);
            assert!(excluded.is_empty());
            let run_itself = t(RUN, "AGCAGCAGC", SECOND, ALTERING);
            let (owned, excluded) = dropped_allele_rows(
                &perl,
                &with(&perl, std::slice::from_ref(&run_itself)),
                &[unchanged],
                MAX_SV_SIZE,
            );
            assert_eq!(excluded, set(&owned, &[run_itself]));
        }

        #[test]
        fn test_the_dropped_allele_is_expanded_as_the_parser_expands_it() {
            let perl = deletion_rows();
            for (info, literal) in [
                ("RN=0,1;RUS=AGC;RB=12", "AGCAGCAGCAGC"),
                ("RN=0,1;RUS=AGC;RUC=4.9", "AGCAGCAGCAGC"),
                ("RN=0,2;RUS=AGC,AG;RUC=3,2", "AGCAGCAGCAGAG"),
                ("RN=0,1;RUS=.;RUC=12", "NNNNNNNNNNNN"),
                ("RN=0,1;RUS=agc;RUC=4", "AGCAGCAGCAGC"),
                ("RN=0,1;RUS=AGC;RUC=2", "AGCAGC"),
            ] {
                let row = t(RUN, literal, SECOND, ALTERING);
                let decoy = t(RUN, &format!("{literal}AGC"), SECOND, ALTERING);
                let (owned, excluded) = dropped_allele_rows(
                    &perl,
                    &with(&perl, &[row.clone(), decoy]),
                    &[two(info)],
                    MAX_SV_SIZE,
                );
                assert_eq!(excluded, set(&owned, &[row]), "{info}");
            }
            let gain_rows = insertion_rows();
            let loss = t(RUN, TR, SECOND, LOSS);
            let row = t(RUN, "AGCAGCAGCAGC", SECOND, ALTERING);
            for info in ["RN=1,0;RUS=AGC;RUC=4", "RN=2,0;RUS=AGC,AGC;RUC=2,2"] {
                let (owned, excluded) = dropped_allele_rows(
                    &gain_rows,
                    &with(&gain_rows, &[loss.clone(), row.clone()]),
                    &[two(info)],
                    MAX_SV_SIZE,
                );
                assert_eq!(excluded, set(&owned, std::slice::from_ref(&loss)), "{info}");
            }
            for info in [
                "RN=1,1;RUS=AGC,AGC;RUC=4,2",
                "RN=0,0;RUS=AGC;RUC=4",
                "RN=.,1;RUS=AGC;RUC=4",
            ] {
                let (_, excluded) = dropped_allele_rows(
                    &perl,
                    &with(&perl, std::slice::from_ref(&row)),
                    &[two(info)],
                    MAX_SV_SIZE,
                );
                assert!(excluded.is_empty(), "{info}");
            }
        }

        #[test]
        fn test_shapes_outside_the_class_are_not_judged() {
            for info in [
                "RN=0,1;RUS=AGC;RUC=2",
                "RN=1,0;RUS=AGC;RUC=3",
                "RN=.,1;RUS=AGC;RUC=4",
                "RN=0,0;RUS=AGC;RUC=4",
                "RN=1,1;RUS=AGC,AGC;RUC=4,2",
                "RUS=AGC;RUC=4",
            ] {
                assert_eq!(cnv_tr_rn_zero_dropped_direction(&two(info)), None, "{info}");
            }
            let trailing = two("RN=1,0;RUS=AGC;RUC=3");
            assert!(cnv_tr_expansion_is_empty(&trailing));
            let owned = owned_sides(&[], &rust());
            let (perl_open, rust_open, masks) = adjusted(&owned, &[trailing], MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_empty_rust, set(&owned, &rust()));
            assert!(masks.rn_zero_rust.is_empty());
            let single = rec(30425587, "<CNV:TR>", &format!("{INFO}RN=0;RUS=AGC;RUC=4"));
            assert_eq!(cnv_tr_rn_zero_dropped_direction(&single), None);
            let record = two("RN=0,1;RUS=AGC;RUC=4");
            let symbolic = vec![t(RUN, TR, SECOND, GAIN)];
            let (_, excluded) = dropped_allele_rows(&symbolic, &rust(), &[record], MAX_SV_SIZE);
            assert!(excluded.is_empty());
        }

        #[test]
        fn test_dropped_literal_and_direction_judge_the_same_shapes() {
            use CopyDirection::{Gain, Loss, Same};
            let two_alts = "<CNV:TR>,<CNV:TR>";
            for (alt, info) in [
                ("<CNV:TR>", "RN=0;RUS=AGC;RUC=4"),
                (two_alts, "RUS=AGC;RUC=4"),
                (two_alts, "RN=.,1;RUS=AGC;RUC=4"),
                (two_alts, "RN=0,1;RUC=4"),
                (two_alts, "RN=0,1;RUS=AGC"),
                ("<CNV:TR>,<CNV:TR>,<CNV:TR>", "RN=0,1,1;RUS=AGC;RUC=4"),
            ] {
                let record = rec(30425587, alt, &format!("{INFO}{info}"));
                assert_eq!(two_allele_cnv_tr_rn(&record), None, "{alt} {info}");
                assert_eq!(
                    cnv_tr_rn_zero_dropped_direction(&record),
                    None,
                    "{alt} {info}"
                );
                assert_eq!(
                    cnv_tr_rn_zero_dropped_literal(&record, MAX_SV_SIZE),
                    None,
                    "{alt} {info}"
                );
            }
            let agc = |units: usize| Some("AGC".repeat(units));
            for (info, direction, literal) in [
                ("RN=0,1;RUS=AGC;RUC=4", Some(Gain), agc(4)),
                ("RN=0,1;RUS=AGC;RUC=3", Some(Same), agc(3)),
                ("RN=0,1;RUS=AGC;RUC=2", None, agc(2)),
                ("RN=0,1;RUS=AGC;RB=12", Some(Gain), agc(4)),
                ("RN=1,0;RUS=AGC;RUC=4", Some(Loss), None),
                ("RN=2,0;RUS=AGC,AGC;RUC=2,2", Some(Loss), None),
                ("RN=1,0;RUS=AGC;RUC=3", None, None),
                ("RN=0,0;RUS=AGC;RUC=4", None, None),
                ("RN=1,1;RUS=AGC,AGC;RUC=4,2", None, None),
                ("RN=0,1;RUS=AGC;RUC=1e400", Some(Gain), None),
                ("RN=0,1;RUS=AGC;RUC=1e18", Some(Gain), None),
                ("RN=0,1;RUS=AGC;RB=10000002", Some(Gain), None),
                ("RN=0,1;RUS=AGC;RUC=3333333.9", Some(Gain), None),
                ("RN=0,1;RUS=AGC;RB=9999999", Some(Gain), agc(3_333_333)),
                (
                    "RN=0,2;RUS=AGC,AGC;RUC=2000000,2000000",
                    Some(Gain),
                    agc(4_000_000),
                ),
            ] {
                let record = two(info);
                assert!(two_allele_cnv_tr_rn(&record).is_some(), "{info}");
                assert_eq!(
                    cnv_tr_rn_zero_dropped_direction(&record),
                    direction,
                    "{info}"
                );
                assert!(
                    cnv_tr_rn_zero_dropped_literal(&record, MAX_SV_SIZE) == literal,
                    "the literal of {info}"
                );
            }
            let gain = two("RN=0,1;RUS=AGC;RUC=2000");
            assert!(cnv_tr_rn_zero_dropped_literal(&gain, 6000) == agc(2000));
            assert_eq!(cnv_tr_rn_zero_dropped_literal(&gain, 5999), None);
            assert_eq!(cnv_tr_rn_zero_dropped_literal(&gain, -1), None);
            assert_eq!(
                cnv_tr_rn_zero_dropped_literal(&gain, 0),
                None,
                "a cap of 0 reads as 5000, under the 6,000 bases"
            );
            let perl = deletion_rows();
            let row = t(RUN, &"AGC".repeat(2000), SECOND, ALTERING);
            let rust = with(&perl, std::slice::from_ref(&row));
            let records = [gain];
            let (owned, excluded) = dropped_allele_rows(&perl, &rust, &records, 6000);
            let only_the_row = set(&owned, &[row]);
            assert_eq!(excluded, only_the_row);
            let (_, excluded) = dropped_allele_rows(&perl, &rust, &records, 5999);
            assert!(excluded.is_empty());
            assert_eq!(
                adjusted(&owned, &records, 6000).2.rn_zero_rust,
                only_the_row
            );
            assert!(adjusted(&owned, &records, 5999).2.rn_zero_rust.is_empty());
            assert_eq!(
                two_allele_cnv_tr_rn(&two("RN=0,2.9;RUS=AGC,AG;RB=6,4")),
                Some(TwoAlleleRn {
                    rus: "AGC,AG",
                    ruc: None,
                    rb: Some("6,4"),
                    counts: [0, 2],
                })
            );
        }

        #[test]
        fn the_two_readers_at_the_edges_of_a_count() {
            use CopyDirection::{Gain, Same};
            let agc = |units: usize| Some("AGC".repeat(units));
            for (info, cap, direction, literal) in [
                (
                    "RN=0,1;RUS=AGC;RUC=0.5",
                    MAX_SV_SIZE,
                    None,
                    Some(String::new()),
                ),
                (
                    "RN=0,1;RUS=AGC;RUC=-5",
                    MAX_SV_SIZE,
                    None,
                    Some(String::new()),
                ),
                ("RN=0,1;RUS=AGC;RUC=-5", -1, None, Some(String::new())),
                ("RN=0,2;RUS=,AGC;RUC=1e30,2", -1, None, None),
                ("RN=0,2;RUS=,AGC;RB=5,6", MAX_SV_SIZE, None, None),
                ("RN=0,2;RUS=,AGC;RUC=1e400,2", MAX_SV_SIZE, None, agc(2)),
                ("RN=-3,1;RUS=AGC;RUC=4", MAX_SV_SIZE, Some(Gain), agc(4)),
                ("RN=0,1e30;RUS=AGC;RUC=4", MAX_SV_SIZE, None, None),
                (
                    "RN=0,1;RUS=\u{df};RUC=2",
                    MAX_SV_SIZE,
                    None,
                    Some("SSSS".to_string()),
                ),
                ("RN=0,1;RUS=AGC;RB=10", MAX_SV_SIZE, Some(Same), agc(3)),
                ("RN=0,1;RUS=AGC;RB=1e400", MAX_SV_SIZE, Some(Gain), None),
                ("RN=1e400,0;RUS=AGC;RUC=4", MAX_SV_SIZE, None, None),
                (
                    "RN=0,1;RUS=AGC;RUC=3333333.4",
                    MAX_SV_SIZE,
                    Some(Gain),
                    None,
                ),
                (
                    "RN=0,1;RUS=AGC;RUC=3333333.34",
                    MAX_SV_SIZE,
                    Some(Gain),
                    None,
                ),
                (
                    "RN=0,1;RUS=AGC;RB=10000000.9",
                    MAX_SV_SIZE,
                    Some(Gain),
                    None,
                ),
                ("RN=0,1;RUS=AGC;RB=10000001", MAX_SV_SIZE, Some(Gain), None),
                ("RN=0,1;RUS=AGC;RUC=4", 0, Some(Gain), agc(4)),
                ("RN=0,1;RUS=AGC;RUC=1667", 0, Some(Gain), None),
                ("RN=0,1;RUS=AGC;RUC=1666.9", 0, Some(Gain), None),
            ] {
                let record = two(info);
                assert_eq!(
                    cnv_tr_rn_zero_dropped_direction(&record),
                    direction,
                    "{info}"
                );
                assert!(
                    cnv_tr_rn_zero_dropped_literal(&record, cap) == literal,
                    "the literal of {info} under {cap}"
                );
            }
        }

        #[test]
        fn a_unit_of_no_bases_repeats_to_nothing_whatever_its_count() {
            let record = two("RN=0,2;RUS=,AGC;RUC=1e30,2");
            assert!(
                cnv_tr_rn_zero_dropped_literal(&record, MAX_SV_SIZE) == Some("AGCAGC".to_string())
            );
        }
    }

    /// The scored side's row of a `<CNV:TR>` allele whose record states it as long as its
    /// run, read without a FASTA: the run's Location, the Allele `tandem_repeat`, the region
    /// terms and neither direction term.
    mod allele_as_long_as_its_run {
        use super::*;

        const RUN: &str = "21:30425588-30425596";
        const INFO: &str = "SVTYPE=CNV;END=30425596;SVLEN=9;";
        const TR: &str = "tandem_repeat";
        const UPSTREAM: &str = "upstream_gene_variant";
        const GAIN: &str = "coding_sequence_variant,feature_elongation";
        const SAME: &str = "coding_sequence_variant";
        const FIRST: &str = "ENST00000334068";
        const SECOND: &str = "ENST00000390690";

        /// A record of one allele over the run.
        fn one(info: &str) -> VcfVariant {
            rec(30425587, "<CNV:TR>", &format!("{INFO}{info}"))
        }

        /// A record of two alleles over the run.
        fn two(info: &str) -> VcfVariant {
            rec(30425587, "<CNV:TR>,<CNV:TR>", &format!("{INFO}{info}"))
        }

        /// The row on a transcript the run lies upstream of, the same for every allele.
        fn upstream() -> Row {
            t(RUN, TR, FIRST, UPSTREAM)
        }

        /// The row of an allele as long as the run inside a coding sequence.
        fn same() -> Row {
            t(RUN, TR, SECOND, SAME)
        }

        /// The row of a gain there, which Ensembl VEP writes for every allele it reads
        /// symbolically (Utils/VariationEffect.pm:246-256).
        fn gain() -> Row {
            t(RUN, TR, SECOND, GAIN)
        }

        /// Beside the reference's own symbolic row over the run the two rows are an
        /// unexpanded pair that differs in `feature_elongation` alone, and both are set
        /// aside. A gain's row both sides hold pairs the row as well and stays on its side.
        #[test]
        fn the_row_pairs_with_the_reference_symbolic_row_that_names_the_gain() {
            assert_eq!(pair_kind(RUN, TR, RUN), Some(PairKind::Unexpanded));
            assert!(is_cnv_tr_representation_pair(
                GAIN,
                SAME,
                PairKind::Unexpanded,
                false,
                false
            ));
            let record = one("RUS=AGC;RUC=3");
            let perl = vec![upstream(), gain()];
            let owned = owned_sides(&perl, &[upstream(), same()]);
            let records = std::slice::from_ref(&record);
            let (perl_open, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_perl, set(&owned, &[gain()]));
            assert_eq!(masks.cnvtr_rust, set(&owned, &[same()]));
            assert!(masks.cnvtr_empty_rust.is_empty());

            let owned = owned_sides(&perl, &[upstream(), gain(), same()]);
            let (perl_open, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert!(masks.cnvtr_perl.is_empty());
            assert_eq!(masks.cnvtr_rust, set(&owned, &[same()]));

            // A region term that differs is a consequence disagreement.
            let other = t(RUN, TR, SECOND, "intron_variant");
            let owned = owned_sides(&perl, &[upstream(), other.clone()]);
            let (perl_open, rust_open, _) = adjusted(&owned, &[], MAX_SV_SIZE);
            assert_eq!(perl_open, set(&owned, &[gain()]));
            assert_eq!(rust_open, set(&owned, &[other]));
        }

        /// With a FASTA the reference writes nothing for a record whose expansion is the run
        /// itself (VariationFeatureOverlap.pm:481): the row is set aside with the record's
        /// other symbolic rows, and stays where a reference row lies beside the run.
        #[test]
        fn the_row_is_set_aside_where_the_reference_wrote_nothing_for_the_unchanged_expansion() {
            let rust = vec![upstream(), same()];
            for record in [
                one("RUS=AGC;RUC=3"),
                one("RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1"),
                one("RUS=AGC;RB=9"),
            ] {
                assert!(cnv_tr_expansion_is_empty(&record), "{}", record.info);
                let owned = owned_sides(&[], &rust);
                let records = std::slice::from_ref(&record);
                let (perl_open, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
                assert!(
                    perl_open.is_empty() && rust_open.is_empty(),
                    "{}",
                    record.info
                );
                assert_eq!(
                    masks.cnvtr_empty_rust,
                    set(&owned, &rust),
                    "{}",
                    record.info
                );
            }
            // A trailing allele of no repeat sequence states no length and names the gain.
            let trailing = two("RN=1,0;RUS=AGC;RUC=3");
            assert!(cnv_tr_expansion_is_empty(&trailing));
            let both = vec![upstream(), same(), gain()];
            let owned = owned_sides(&[], &both);
            let (perl_open, rust_open, masks) = adjusted(&owned, &[trailing], MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.cnvtr_empty_rust, set(&owned, &both));
            assert!(masks.rn_zero_rust.is_empty());
            // An intron's row of another record at the base after the run: the reference
            // wrote a row beside the run, and the mask leaves the record's rows.
            let beside = vec![t("21:30425597", "A", SECOND, "intron_variant")];
            let record = one("RUS=AGC;RUC=3");
            let (_, excluded) = empty_expansion_rows(&beside, &rust, &[record]);
            assert!(excluded.is_empty());
        }

        /// A two-allele record whose first allele has no repeat sequence and whose real
        /// allele is as long as the run. With a FASTA the reference keeps the deletion of the
        /// run and drops the real allele (Parser/VCF.pm:411, 428; Parser.pm:881, 978). The
        /// scored side writes the first allele, which states no length, as a gain and the
        /// real allele without a direction term: the row without the term is the dropped
        /// allele's, and the gain's row pairs with the reference's deletion.
        #[test]
        fn the_row_is_the_dropped_allele_behind_a_leading_allele_of_no_repeat_sequence() {
            let record = two("RN=0,1;RUS=AGC;RUC=3");
            assert!(!cnv_tr_expansion_is_empty(&record));
            assert_eq!(
                cnv_tr_rn_zero_dropped_direction(&record),
                Some(CopyDirection::Same)
            );
            let perl = vec![
                t(RUN, "-", FIRST, UPSTREAM),
                t(RUN, "-", SECOND, "inframe_deletion"),
            ];
            let rust = vec![upstream(), gain(), same()];
            let owned = owned_sides(&perl, &rust);
            let records = std::slice::from_ref(&record);
            let (perl_open, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
            assert!(perl_open.is_empty() && rust_open.is_empty());
            assert_eq!(masks.rn_zero_rust, set(&owned, &[same()]));
            assert_eq!(masks.cnvtr_perl, set(&owned, &perl));
            assert_eq!(masks.cnvtr_rust, set(&owned, &rust));
        }

        /// Beside a literal reading of equal length: a record whose units are not the run's
        /// is written with a FASTA as its whole allele over the run, and the row pairs with
        /// it as the symbolic reading of the same repeat.
        #[test]
        fn the_row_pairs_with_a_literal_reading_as_long_as_the_run() {
            let record = one("RUS=123;RUC=3");
            let literal = "123123123";
            assert_eq!(pair_kind(RUN, literal, RUN), Some(PairKind::Literal));
            for perl_csq in [
                "protein_altering_variant",
                "missense_variant",
                "synonymous_variant",
                "coding_sequence_variant",
            ] {
                let perl = vec![
                    t(RUN, literal, FIRST, UPSTREAM),
                    t(RUN, literal, SECOND, perl_csq),
                ];
                let rust = vec![upstream(), same()];
                let owned = owned_sides(&perl, &rust);
                let records = std::slice::from_ref(&record);
                let (perl_open, rust_open, masks) = adjusted(&owned, records, MAX_SV_SIZE);
                assert!(perl_open.is_empty() && rust_open.is_empty(), "{perl_csq}");
                assert_eq!(masks.cnvtr_perl, set(&owned, &perl), "{perl_csq}");
                assert_eq!(masks.cnvtr_rust, set(&owned, &rust), "{perl_csq}");
                assert!(masks.cnvtr_empty_rust.is_empty(), "{perl_csq}");
            }
            // A splice term the literal reading alone carries is a consequence disagreement.
            assert!(!is_cnv_tr_representation_pair(
                "splice_region_variant,missense_variant",
                SAME,
                PairKind::Literal,
                false,
                false
            ));
        }
    }
}
