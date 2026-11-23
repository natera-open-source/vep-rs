// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The terms a point breakend's mate-side Transcript row supports, derived from the
//! transcript models of a JSON cache: read at the mate position, and read at the local
//! position as Ensembl VEP reads it.
//!
//! A mate-side row is a Transcript row of a bracket-notation breakend whose Feature lies
//! within [`MAX_DISTANCE_FROM_TRANSCRIPT`] of the mate coordinate on the mate chromosome
//! while the local coordinate does not itself select that transcript (another chromosome,
//! or farther than that distance on the same one). Ensembl VEP builds such a row for every
//! transcript near the mate (ensembl-variation StructuralVariationOverlap.pm:86-102,
//! `_close_to_feature`) and then evaluates every positional predicate on it with the local
//! variation feature (BaseVariationFeatureOverlapAllele.pm:257,273 hand the local feature
//! to each predicate; `upstream` and `downstream` at Utils/VariationEffect.pm:444-457 and
//! `_bvfo_preds` at BaseVariationFeatureOverlapAllele.pm:454 read its start and end); only
//! `feature_truncation` (Utils/VariationEffect.pm:358) reads the mate breakend.
//!
//! [`breakend_context`] gives the term set a point breakend at a position supports on a
//! transcript, with `feature_truncation` decided by whether the mate lies inside the
//! transcript. Called with the mate position it is the set the mate breakend supports;
//! called with the local position it is Ensembl VEP's read, whose empty result becomes the
//! default `intergenic_variant` on the row. [`mate_side_derivation`] applies both to a row.
//!
//! The transcript structure comes from the JSON cache both engines annotated against:
//! exons, strand, biotype, the genomic CDS bounds mapped from the mapper's
//! `cdna_coding_start` and `cdna_coding_end` through its cDNA-to-genome pairs, and the
//! mature miRNA ranges of the `miRNA` attribute (cDNA coordinates).
//!
//! What the derivation does not model, so that a row reaching one of these edges takes no
//! class from it (a pair is named only when both engines' sets equal the derivations):
//!
//! - a breakend inside a start or stop codon: Ensembl VEP fires `start_lost` and
//!   `stop_lost` on span overlap (Utils/VariationEffect.pm:886-896, the stop analogue) and
//!   `start_retained_variant` through `_ins_del_start_altered` returning 0 on a structural
//!   allele (:1037);
//! - the frameshift-intron exon lookup: a transcript with an intron of 13 bases or fewer
//!   stretches its exon lookup by 12 bases (BaseTranscriptVariation.pm:911-914); the exon
//!   and intron tests here are exact;
//! - a multi-exon miRNA: the mature range is a cDNA interval that is not read off a
//!   multi-exon span, so the non-coding exon term stands;
//! - a same-chromosome breakend whose local coordinate is also within the gate of the
//!   transcript: the row is a local-side row and [`mate_side_derivation`] gives nothing.
//!
//! `_close_to_feature` expands the feature slice by 5,000 on each side and tests overlap
//! inclusively, as [`TranscriptModel::near`] does: a breakend exactly 5,000 bases outside
//! is inside the gate on both readings.
//!
//! A number of a cache record is read as `int()` reads the value and is held within 64
//! bits; a position written in a Location or in an allele is read at any length up to the
//! digits `int()` converts, one past 128 bits as the largest position there is.

use crate::classify::cache::iter_cache_transcripts;
use crate::error::Failure;
use crate::json::{Object, Value, MAX_INT_DIGITS};
use crate::pyfmt::{is_py_space, py_strip};
use crate::pytext::{py_repr, py_str};
use crate::rules::breakend::strip_chr;
use regex::Regex;
use std::collections::{BTreeSet, HashMap};
use std::sync::LazyLock;

/// Ensembl VEP's breakend-to-feature distance (Utils/VariationEffect.pm:60, used by
/// StructuralVariationOverlap.pm:157), which is also its upstream and downstream window.
pub const MAX_DISTANCE_FROM_TRANSCRIPT: i128 = 5000;

/// The two bases at each end of an intron: a position inside them is no `intron_variant`.
const SPLICE_SITE: i128 = 2;

/// The polypyrimidine tract: the acceptor-side window 3 to 17 bases into the intron
/// (Utils/VariationEffect.pm `_intron_overlap`, the tract branch).
const PPT_INNER: i128 = 2;
const PPT_OUTER: i128 = 16;

/// The class of a mate-side row Ensembl VEP reads at the local coordinate.
pub const BREAKEND_MATE_LOCAL_READ: &str = "breakend_mate_local_read";

/// A set of consequence terms.
pub type Terms = BTreeSet<&'static str>;

/// A transcript as the cache files it, with the genomic CDS bounds resolved.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptModel {
    /// The chromosome the cache files the transcript under.
    pub chrom: String,
    pub start: i128,
    pub end: i128,
    pub strand: i128,
    pub biotype: String,
    /// The exons as `(start, end)`, ascending.
    pub exons: Vec<(i128, i128)>,
    /// The genomic bounds of the coding sequence, low first; `None` without one.
    pub cds: Option<(i128, i128)>,
    /// The mature miRNA ranges as `(first, last)` in cDNA coordinates.
    pub mirna_mature: Vec<(i128, i128)>,
}

impl TranscriptModel {
    /// The gaps between neighbouring exons, as `(start, end)`; exons that abut or overlap
    /// leave none.
    pub fn introns(&self) -> Vec<(i128, i128)> {
        self.exons
            .windows(2)
            .filter_map(|pair| {
                let (start, end) = (pair[0].1.saturating_add(1), pair[1].0.saturating_sub(1));
                (end >= start).then_some((start, end))
            })
            .collect()
    }

    /// Whether `pos` is inside the transcript or within [`MAX_DISTANCE_FROM_TRANSCRIPT`] of
    /// either end.
    pub fn near(&self, pos: i128) -> bool {
        self.start.saturating_sub(MAX_DISTANCE_FROM_TRANSCRIPT) <= pos
            && pos <= self.end.saturating_add(MAX_DISTANCE_FROM_TRANSCRIPT)
    }

    pub fn inside(&self, pos: i128) -> bool {
        self.start <= pos && pos <= self.end
    }
}

/// The transcripts of a cache by chromosome and stable id. The chromosome is part of the
/// key because a pseudoautosomal transcript is filed under both X and Y with its own
/// coordinates on each.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct TranscriptModels {
    by_chromosome: HashMap<String, HashMap<String, TranscriptModel>>,
}

impl TranscriptModels {
    /// The model filed under `chrom` with the stable id `stable_id`.
    pub fn get(&self, chrom: &str, stable_id: &str) -> Option<&TranscriptModel> {
        self.by_chromosome.get(chrom)?.get(stable_id)
    }

    /// Files `model` under `chrom` and `stable_id`, in the place of one filed there before.
    pub fn insert(&mut self, chrom: &str, stable_id: &str, model: TranscriptModel) {
        self.by_chromosome
            .entry(chrom.to_string())
            .or_default()
            .insert(stable_id.to_string(), model);
    }

    pub fn len(&self) -> usize {
        self.by_chromosome.values().map(HashMap::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

static DECIMAL_DIGIT: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\d$").expect("regex"));

fn is_decimal_digit(c: char) -> bool {
    let mut buffer = [0u8; 4];
    DECIMAL_DIGIT.is_match(c.encode_utf8(&mut buffer))
}

/// The value of a decimal digit of any script, the characters `int()` reads. The decimal
/// digits of a script are ten consecutive characters from its zero, and the runs of
/// neighbouring scripts follow each other in whole tens, so a digit is its distance from
/// the start of its run, modulo ten.
pub fn decimal_digit(c: char) -> Option<u32> {
    if c.is_ascii() {
        return c.to_digit(10);
    }
    if !is_decimal_digit(c) {
        return None;
    }
    let mut zero = c as u32;
    while char::from_u32(zero - 1).is_some_and(is_decimal_digit) {
        zero -= 1;
    }
    Some((c as u32 - zero) % 10)
}

/// Whether `text` is one decimal digit or more, of any script.
fn is_decimal_run(text: &str) -> bool {
    !text.is_empty() && text.chars().all(|c| decimal_digit(c).is_some())
}

/// `int()` of a run of decimal digits: its value, one past 128 bits read as the largest
/// value there is; `None` for a run of more digits than `int()` converts and for a text that
/// is no run.
fn decimal_value(digits: &str) -> Option<i128> {
    if !is_decimal_run(digits) || digits.chars().count() > MAX_INT_DIGITS {
        return None;
    }
    let mut value: i128 = 0;
    for digit in digits.chars().filter_map(decimal_digit) {
        value = value.saturating_mul(10).saturating_add(i128::from(digit));
    }
    Some(value)
}

/// What `int()` makes of a text.
enum IntText {
    Value(i128),
    /// An integer 128 bits do not hold.
    Long,
    /// No integer: `int()` raises `ValueError`.
    Refused,
}

/// `int(text)`: the blanks `int()` drops around the number (those of `str.isspace()` but
/// U+001C to U+001F), an optional sign, then decimal digits of any script with single
/// underscores between them, 4,300 digits at most.
fn int_of_text(text: &str) -> IntText {
    let body = text.trim_matches(|c| is_py_space(c) && !('\u{1c}'..='\u{1f}').contains(&c));
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
    {
        return IntText::Refused;
    }
    let mut magnitude = Some(0i128);
    let mut count = 0usize;
    for c in digits.chars().filter(|c| *c != '_') {
        let Some(digit) = decimal_digit(c) else {
            return IntText::Refused;
        };
        count += 1;
        magnitude = magnitude
            .and_then(|value| value.checked_mul(10))
            .and_then(|value| value.checked_add(i128::from(digit)));
    }
    if count > MAX_INT_DIGITS {
        return IntText::Refused;
    }
    match magnitude {
        Some(value) => IntText::Value(if negative { -value } else { value }),
        None => IntText::Long,
    }
}

const INFINITE: &str = "holds a number that is infinite";
const PAST_64_BITS: &str = "holds a number that 64 bits do not hold";

/// A number of a cache record as `int()` reads the value: an integer, a float cut at its
/// point, `true` and `false` as 1 and 0, a text as `int()` reads one; `None` for a value
/// `int()` refuses with a `TypeError` or a `ValueError`. The reason is returned for a number
/// the tool does not hold: an infinite one, on which `int()` ends in an `OverflowError`, and
/// one past 64 bits.
fn cache_int(value: Option<&Value>) -> Result<Option<i128>, &'static str> {
    let number = match value {
        Some(Value::Bool(flag)) => i128::from(*flag),
        Some(Value::Int(digits)) => digits.text().parse().map_err(|_| PAST_64_BITS)?,
        Some(Value::Float(number)) if number.is_nan() => return Ok(None),
        Some(Value::Float(number)) if number.is_infinite() => return Err(INFINITE),
        Some(Value::Float(number)) => number.trunc() as i128,
        Some(Value::Str(text)) => match int_of_text(text) {
            IntText::Value(number) => number,
            IntText::Long => return Err(PAST_64_BITS),
            IntText::Refused => return Ok(None),
        },
        _ => return Ok(None),
    };
    match i64::try_from(number) {
        Ok(_) => Ok(Some(number)),
        Err(_) => Err(PAST_64_BITS),
    }
}

/// The objects of the list `holder` keeps under `key`, none for a value Python reads as
/// false. The reason is returned for any other value and for a list with an item that is
/// no object, on which the read of an item's member ends in an exception.
fn objects_under<'a>(holder: &'a Object, key: &str) -> Result<Vec<&'a Object>, String> {
    let refused = || format!("holds {key} that are no list of objects");
    match holder.get(key).filter(|held| held.truthy()) {
        None => Ok(Vec::new()),
        Some(Value::List(items)) => items
            .iter()
            .map(|item| item.as_object().ok_or_else(refused))
            .collect(),
        Some(_) => Err(refused()),
    }
}

/// The object `holder` keeps under `key`, `None` for a value Python reads as false; the
/// reason is returned for any other value.
fn object_under<'a>(holder: &'a Object, key: &str) -> Result<Option<&'a Object>, String> {
    match holder.get(key).filter(|held| held.truthy()) {
        None => Ok(None),
        Some(Value::Object(object)) => Ok(Some(object)),
        Some(_) => Err(format!("holds a {key} that is no object")),
    }
}

/// The genomic bounds of the coding sequence, low first: the mapper's two cDNA coding
/// bounds, each carried through the pair that holds it (the last such pair, should several
/// hold it). A pair maps the cDNA interval `from_start` to `from_end` onto the genomic
/// interval `to_start` to `to_end` in its orientation `ori`: a cDNA offset inside the pair
/// lands at `to_start + offset` on the forward strand and at `to_end - offset` on the
/// reverse. `None` when a bound is missing, is not above zero or lies in no pair.
pub fn genomic_cds(tx: &Object) -> Result<Option<(i128, i128)>, String> {
    let Some(cache) = object_under(tx, "variation_effect_feature_cache")? else {
        return Ok(None);
    };
    let Some(mapper) = object_under(cache, "mapper")? else {
        return Ok(None);
    };
    let coding_start = cache_int(mapper.get("cdna_coding_start"))?;
    let coding_end = cache_int(mapper.get("cdna_coding_end"))?;
    let (Some(cs), Some(ce)) = (coding_start, coding_end) else {
        return Ok(None);
    };
    if cs <= 0 || ce <= 0 {
        return Ok(None);
    }
    let (mut gs, mut ge) = (None, None);
    for pair in objects_under(mapper, "pairs")? {
        let mut numbers = [None; 5];
        for (number, key) in
            numbers
                .iter_mut()
                .zip(["from_start", "from_end", "to_start", "to_end", "ori"])
        {
            *number = cache_int(pair.get(key))?;
        }
        let [Some(fs), Some(fe), Some(ts), Some(te), Some(ori)] = numbers else {
            continue;
        };
        let genomic = |cdna: i128| {
            if ori >= 0 {
                ts + (cdna - fs)
            } else {
                te - (cdna - fs)
            }
        };
        if fs <= cs && cs <= fe {
            gs = Some(genomic(cs));
        }
        if fs <= ce && ce <= fe {
            ge = Some(genomic(ce));
        }
    }
    let (Some(gs), Some(ge)) = (gs, ge) else {
        return Ok(None);
    };
    Ok(Some(if gs <= ge { (gs, ge) } else { (ge, gs) }))
}

/// The mature miRNA ranges of a record: the value of each attribute whose code is `miRNA`,
/// written `first-last` in decimal digits with blanks around it dropped.
fn mirna_mature(tx: &Object) -> Result<Vec<(i128, i128)>, String> {
    let mut ranges = Vec::new();
    for attribute in objects_under(tx, "attributes")? {
        if attribute.get("code").and_then(Value::as_str) != Some("miRNA") {
            continue;
        }
        let value = attribute.get("value").map(py_str).unwrap_or_default();
        let Some((first, last)) = py_strip(&value).split_once('-') else {
            continue;
        };
        if !(is_decimal_run(first) && is_decimal_run(last)) {
            continue;
        }
        match (decimal_value(first), decimal_value(last)) {
            (Some(first), Some(last)) => ranges.push((first, last)),
            _ => return Err("holds a miRNA range of more than 4,300 digits".to_string()),
        }
    }
    Ok(ranges)
}

/// The failure of a cache record the tool does not read.
fn refused(chrom: &str, tx: &Object, detail: &str) -> Failure {
    let name = tx
        .get("stable_id")
        .map_or_else(|| "without a stable_id".to_string(), py_repr);
    Failure::own(format!(
        "the transcript {name} filed under {chrom} in the corpus cache {detail}"
    ))
}

/// The model of one transcript record of a cache shard filed under `chrom`; `None` for a
/// record that gives no integer start, end or strand. A record is refused in the tool's
/// own words where its read ends in a CPython exception: a number that is infinite; `exons`,
/// `attributes` or the mapper's `pairs` that are no list of objects; a
/// `variation_effect_feature_cache` or a `mapper` that is no object; a miRNA range of more
/// than 4,300 digits. A number past 64 bits is refused as well.
pub fn model_from_cache_record(
    chrom: &str,
    tx: &Object,
) -> Result<Option<TranscriptModel>, Failure> {
    let read = || -> Result<Option<TranscriptModel>, String> {
        let start = cache_int(tx.get("start"))?;
        let end = cache_int(tx.get("end"))?;
        let strand = cache_int(tx.get("strand"))?;
        let (Some(start), Some(end), Some(strand)) = (start, end, strand) else {
            return Ok(None);
        };
        let mut exons = Vec::new();
        for exon in objects_under(tx, "exons")? {
            let (first, last) = (cache_int(exon.get("start"))?, cache_int(exon.get("end"))?);
            if let (Some(first), Some(last)) = (first, last) {
                exons.push((first, last));
            }
        }
        exons.sort_unstable();
        let biotype = tx.get("biotype").filter(|held| held.truthy());
        Ok(Some(TranscriptModel {
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            biotype: biotype.map(py_str).unwrap_or_default(),
            exons,
            cds: genomic_cds(tx)?,
            mirna_mature: mirna_mature(tx)?,
        }))
    };
    read().map_err(|detail| refused(chrom, tx, &detail))
}

/// Every transcript of the cache at `cache_dir` as a model, by the chromosome it is filed
/// under and its stable id; a later record of a key stands in the place of an earlier one.
/// None for a cache without a `transcripts` directory. A record that has a model and whose
/// stable id is a list or an object is refused in the tool's own words; one whose stable id
/// is a number is the model of no Feature.
pub fn load_transcript_models(cache_dir: &str) -> Result<TranscriptModels, Failure> {
    let mut models = TranscriptModels::default();
    for item in iter_cache_transcripts(cache_dir) {
        let (chrom, tx) = item?;
        let Some(model) = model_from_cache_record(&chrom, &tx)? else {
            continue;
        };
        match tx.get("stable_id") {
            Some(Value::Str(stable_id)) => models.insert(&chrom, stable_id, model),
            Some(Value::List(_) | Value::Object(_)) => {
                return Err(refused(&chrom, &tx, "is named by a list or an object"))
            }
            _ => {}
        }
    }
    Ok(models)
}

static BRACKET_MATE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]").expect("regex"));

/// The mate `(chromosome, position)` a bracket-notation breakend allele names, the
/// chromosome without a leading `chr` in any case; `None` for any other allele and for a
/// position of more digits than `int()` converts.
pub fn bracket_mate(allele: &str) -> Option<(&str, i128)> {
    let found = BRACKET_MATE.captures(allele)?;
    let chrom = found.get(1).expect("chromosome group").as_str();
    Some((strip_chr(chrom), decimal_value(&found[2])?))
}

/// `(chromosome, position)` of a point Location, the chromosome without a leading `chr` in
/// any case; `None` for a ranged Location, for any other text and for a position `int()`
/// does not convert.
pub fn point_location(loc: &str) -> Option<(&str, i128)> {
    let (chrom, coords) = loc.split_once(':')?;
    Some((strip_chr(chrom), decimal_value(coords)?))
}

/// Whether `pos` lies in a mature range of a single-exon miRNA; a multi-exon miRNA has no
/// cDNA position to read off its span.
fn mature_mirna_hit(pos: i128, model: &TranscriptModel) -> bool {
    if model.biotype != "miRNA" || model.exons.len() != 1 {
        return false;
    }
    let cdna = if model.strand == 1 {
        pos.saturating_sub(model.start).saturating_add(1)
    } else {
        model.end.saturating_sub(pos).saturating_add(1)
    };
    model
        .mirna_mature
        .iter()
        .any(|&(first, last)| first <= cdna && cdna <= last)
}

/// The terms a point breakend at `pos`, inside the transcript, carries besides
/// `feature_truncation`: the exon terms (CDS, 5' or 3' UTR by strand, the mature-miRNA term
/// alone, or the non-coding exon term), or the intron terms (`intron_variant` strictly
/// inside the splice sites, the polypyrimidine tract window, and the non-coding transcript
/// term for a transcript without a CDS that is not NMD), plus `NMD_transcript_variant` on a
/// nonsense-mediated-decay transcript.
pub fn region_terms(pos: i128, model: &TranscriptModel) -> Terms {
    let mut terms = Terms::new();
    let nmd = model.biotype == "nonsense_mediated_decay";
    let within = |&(start, end): &(i128, i128)| start <= pos && pos <= end;
    if model.exons.iter().any(within) {
        match model.cds {
            Some((cs, ce)) => {
                terms.insert(if cs <= pos && pos <= ce {
                    "coding_sequence_variant"
                } else if (model.strand == 1 && pos < cs) || (model.strand == -1 && pos > ce) {
                    "5_prime_UTR_variant"
                } else {
                    "3_prime_UTR_variant"
                });
            }
            None => {
                if mature_mirna_hit(pos, model) {
                    return Terms::from(["mature_miRNA_variant"]);
                }
                terms.insert("non_coding_transcript_exon_variant");
            }
        }
    } else if let Some((start, end)) = model.introns().into_iter().find(within) {
        if start.saturating_add(SPLICE_SITE) <= pos && pos <= end.saturating_sub(SPLICE_SITE) {
            terms.insert("intron_variant");
        }
        let in_tract = match model.strand {
            1 => end.saturating_sub(PPT_OUTER) <= pos && pos <= end.saturating_sub(PPT_INNER),
            -1 => start.saturating_add(PPT_INNER) <= pos && pos <= start.saturating_add(PPT_OUTER),
            _ => false,
        };
        if in_tract {
            terms.insert("splice_polypyrimidine_tract_variant");
        }
        if model.cds.is_none() && !nmd {
            terms.insert("non_coding_transcript_variant");
        }
    }
    if nmd {
        terms.insert("NMD_transcript_variant");
    }
    terms
}

/// The term set a point breakend evaluated at `pos` supports on the transcript.
///
/// `mate_inside` says whether the mate breakend lies inside the transcript, which is what
/// decides `feature_truncation` (Utils/VariationEffect.pm:358 reads the breakend, not the
/// evaluated position). Inside the transcript the region terms apply; within
/// [`MAX_DISTANCE_FROM_TRANSCRIPT`] outside it the strand-aware upstream or downstream
/// term; beyond that nothing but the truncation term. A mature-miRNA hit stands alone.
pub fn breakend_context(pos: i128, model: &TranscriptModel, mate_inside: bool) -> Terms {
    if model.inside(pos) {
        let mut terms = region_terms(pos, model);
        if mate_inside && terms != Terms::from(["mature_miRNA_variant"]) {
            terms.insert("feature_truncation");
        }
        return terms;
    }
    let mut terms = Terms::new();
    let forward = model.strand == 1;
    if model.start.saturating_sub(MAX_DISTANCE_FROM_TRANSCRIPT) <= pos && pos < model.start {
        terms.insert(if forward {
            "upstream_gene_variant"
        } else {
            "downstream_gene_variant"
        });
    } else if model.end < pos && pos <= model.end.saturating_add(MAX_DISTANCE_FROM_TRANSCRIPT) {
        terms.insert(if forward {
            "downstream_gene_variant"
        } else {
            "upstream_gene_variant"
        });
    }
    if mate_inside {
        terms.insert("feature_truncation");
    }
    terms
}

/// `(mate_set, local_set)` for a mate-side row, or `None` when the row is not one.
///
/// The row is a mate-side row when the Location is a point, the Allele is bracket notation,
/// the transcript is on the mate chromosome within the distance gate of the mate, and the
/// local coordinate does not itself select the transcript. `local_set` is Ensembl VEP's
/// read with its empty result replaced by `intergenic_variant`, the default a Transcript
/// row carries when no predicate fires.
pub fn mate_side_derivation(
    loc: &str,
    allele: &str,
    model: &TranscriptModel,
) -> Option<(Terms, Terms)> {
    let (local_chrom, local_pos) = point_location(loc)?;
    let (mate_chrom, mate_pos) = bracket_mate(allele)?;
    if model.chrom != mate_chrom || !model.near(mate_pos) {
        return None;
    }
    if local_chrom == model.chrom && model.near(local_pos) {
        return None;
    }
    let mate_inside = model.inside(mate_pos);
    let mate_set = breakend_context(mate_pos, model, mate_inside);
    let mut local_set = breakend_context(local_pos, model, mate_inside);
    if local_set.is_empty() {
        local_set.insert("intergenic_variant");
    }
    Some((mate_set, local_set))
}

#[cfg(test)]
#[path = "mate_context_tests.rs"]
mod tests;
