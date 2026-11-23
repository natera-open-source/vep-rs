// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a point breakend's mate-side Transcript row supports, derived from the candidate
//! JSON cache as `tools/vep-sv-compare` derives it (`mate_context`): the transcripts
//! of the cache, each one's model, and the two term sets of a mate-side row.

use super::number::{minus, nearest, plus, whole_number, Int};
use super::TranscriptModel;
use crate::cache;
use crate::fields::{member, python_type};
use crate::input;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;
use std::panic;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::LazyLock;

/// Ensembl VEP's breakend-to-feature distance (Utils/VariationEffect.pm:60, used by
/// StructuralVariationOverlap.pm:157) and its up/downstream window, which is the same value.
pub const MAX_DISTANCE_FROM_TRANSCRIPT: Int = 5000;

/// The intron positions Ensembl VEP's structural path reads as `intron_variant`: strictly
/// inside the two-base splice sites. The polypyrimidine tract is the acceptor-side window 3
/// to 17 bases into the intron (Utils/VariationEffect.pm `_intron_overlap`, the tract branch).
const _SPLICE_SITE: Int = 2;
const _PPT_INNER: Int = 2;
const _PPT_OUTER: Int = 16;

static BRACKET_MATE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]").expect("the pattern compiles")
});

static MIRNA_RANGE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(r"^[\s\x1c-\x1f]*(\d+)-(\d+)[\s\x1c-\x1f]*\z").expect("the pattern compiles")
});

impl TranscriptModel {
    /// The gaps between neighbouring exons, in the exons' order.
    pub fn introns(&self) -> Vec<(Int, Int)> {
        self.exons
            .iter()
            .zip(self.exons.iter().skip(1))
            .map(|((_s1, e1), (s2, _e2))| (plus(*e1, 1), minus(*s2, 1)))
            .filter(|(start, end)| end >= start)
            .collect()
    }

    pub fn near(&self, pos: Int) -> bool {
        minus(self.start, MAX_DISTANCE_FROM_TRANSCRIPT) <= pos
            && pos <= plus(self.end, MAX_DISTANCE_FROM_TRANSCRIPT)
    }

    pub fn inside(&self, pos: Int) -> bool {
        self.start <= pos && pos <= self.end
    }
}

/// `int(v)`, `None` where it raises `TypeError` or `ValueError`. An infinite float ends
/// the run as the interpreter's `OverflowError`.
pub(super) fn _int(v: Option<&Value>) -> Option<Int> {
    match v? {
        Value::Null | Value::List(_) | Value::Dict(_) => None,
        Value::Bool(flag) => Some(Int::from(*flag)),
        Value::Int(whole) => Some(Int::from(*whole)),
        Value::BigInt(whole) => Some(nearest(whole)),
        Value::Float(float) if float.is_nan() => None,
        Value::Float(float) if float.is_infinite() => {
            crate::uncaught!("cannot convert float infinity to integer")
        }
        Value::Float(float) => Some(float.trunc() as Int),
        Value::Str(text) => pyfmt::int_checked(text).ok().map(|whole| nearest(&whole)),
    }
}

/// The items of `for item in value or []` where the loop's body calls `item.get`: the
/// items of a list, none for a value that is absent or reads as false. A number cannot be
/// walked, and a text or an object is walked as texts, which have no `get`.
fn records(value: Option<&Value>) -> &[Value] {
    match value {
        Some(Value::List(items)) => items,
        Some(found @ (Value::Bool(_) | Value::Int(_) | Value::BigInt(_) | Value::Float(_)))
            if pyjson::truthy(found) =>
        {
            crate::uncaught!("'{}' object is not iterable", python_type(found))
        }
        Some(found @ (Value::Str(_) | Value::Dict(_))) if pyjson::truthy(found) => {
            crate::uncaught!("'str' object has no attribute 'get'")
        }
        _ => &[],
    }
}

/// Map the mapper's cDNA coding bounds to genomic coordinates through its pairs.
///
/// Each pair maps a cDNA interval (`from_start`..`from_end`) onto a genomic interval
/// (`to_start`..`to_end`) in the pair's orientation; a cDNA offset inside the pair lands
/// at `to_start + offset` on the forward strand and `to_end - offset` on the reverse.
/// The two mapped bounds are returned in ascending genomic order.
fn _genomic_cds(tx: &Value) -> Option<(Int, Int)> {
    let empty = Value::Dict(IndexMap::new());
    let cache = match tx.get("variation_effect_feature_cache") {
        Some(found) if pyjson::truthy(found) => found,
        _ => &empty,
    };
    let mapper = match member(cache, "mapper") {
        Some(found) if pyjson::truthy(found) => found,
        _ => &empty,
    };
    let cs = _int(member(mapper, "cdna_coding_start"));
    let ce = _int(member(mapper, "cdna_coding_end"));
    let (cs, ce) = match (cs, ce) {
        (Some(cs), Some(ce)) if cs > 0 && ce > 0 => (cs, ce),
        _ => return None,
    };
    let mut gs = None;
    let mut ge = None;
    for pair in records(member(mapper, "pairs")) {
        let fs = _int(member(pair, "from_start"));
        let fe = _int(member(pair, "from_end"));
        let ts = _int(member(pair, "to_start"));
        let te = _int(member(pair, "to_end"));
        let ori = _int(member(pair, "ori"));
        let (Some(fs), Some(fe), Some(ts), Some(te), Some(ori)) = (fs, fe, ts, te, ori) else {
            continue;
        };
        let mapped = |cdna: Int| {
            if ori >= 0 {
                plus(ts, minus(cdna, fs))
            } else {
                minus(te, minus(cdna, fs))
            }
        };
        if fs <= cs && cs <= fe {
            gs = Some(mapped(cs));
        }
        if fs <= ce && ce <= fe {
            ge = Some(mapped(ce));
        }
    }
    let (gs, ge) = (gs?, ge?);
    Some(if gs <= ge { (gs, ge) } else { (ge, gs) })
}

fn _mirna_mature(tx: &Value) -> Vec<(Int, Int)> {
    let mut out = Vec::new();
    for attr in records(tx.get("attributes")) {
        if member(attr, "code").and_then(Value::as_str) != Some("miRNA") {
            continue;
        }
        let value = member(attr, "value")
            .map(pyjson::py_str)
            .unwrap_or_default();
        if let Some(range) = MIRNA_RANGE.captures(&value) {
            out.push((whole_number(&range[1]), whole_number(&range[2])));
        }
    }
    out
}

/// Build a `TranscriptModel` from one transcript record of a cache shard.
pub fn model_from_cache_record(chrom: &str, tx: &Value) -> Option<TranscriptModel> {
    let start = _int(tx.get("start"));
    let end = _int(tx.get("end"));
    let strand = _int(tx.get("strand"));
    let (start, end, strand) = (start?, end?, strand?);
    let mut exons = Vec::new();
    for e in records(tx.get("exons")) {
        let s = _int(member(e, "start"));
        let t = _int(member(e, "end"));
        if let (Some(s), Some(t)) = (s, t) {
            exons.push((s, t));
        }
    }
    exons.sort_unstable();
    exons.shrink_to_fit();
    let biotype = match tx.get("biotype") {
        Some(found) if pyjson::truthy(found) => pyjson::py_str(found),
        _ => String::new(),
    };
    Some(TranscriptModel {
        chrom: chrom.to_string(),
        start,
        end,
        strand,
        biotype,
        exons,
        cds: _genomic_cds(tx),
        mirna_mature: _mirna_mature(tx),
    })
}

/// A name as the code points of the text the interpreter holds it as, a byte that is not
/// UTF-8 standing as the lone surrogate U+DC00 plus its value. Two names sort as these
/// sequences do.
fn name_order(name: &OsStr) -> Vec<u32> {
    let mut points = Vec::new();
    for chunk in name.as_bytes().utf8_chunks() {
        points.extend(chunk.valid().chars().map(u32::from));
        points.extend(chunk.invalid().iter().map(|byte| 0xDC00 + u32::from(*byte)));
    }
    points
}

/// The shards of the cache in the order `iter_cache_transcripts` reads them, each with
/// the name of its directory: the directories of `transcripts/` by name, and in each its
/// `*.json` and `*.json.gz` together by name.
fn shards_in_order(root: &Path) -> Vec<(String, PathBuf)> {
    let mut by_directory: IndexMap<PathBuf, Vec<PathBuf>> = IndexMap::new();
    for shard in cache::shards(root, None).unwrap_or_default() {
        let Some(directory) = shard.parent() else {
            continue;
        };
        by_directory
            .entry(directory.to_path_buf())
            .or_default()
            .push(shard);
    }
    let last_name = |path: &Path| name_order(path.file_name().unwrap_or_default());
    let mut directories: Vec<(PathBuf, Vec<PathBuf>)> = by_directory.into_iter().collect();
    directories.sort_by_cached_key(|(directory, _)| last_name(directory));
    let mut ordered = Vec::new();
    for (directory, mut shards) in directories {
        let chromosome = directory
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        shards.sort_by_cached_key(|shard| last_name(shard));
        ordered.extend(shards.into_iter().map(|shard| (chromosome.clone(), shard)));
    }
    ordered
}

/// The value a shard holds, `None` for a shard `iter_cache_transcripts` skips: one that
/// cannot be read, is not gzip where its name says so, is cut short, is not UTF-8 or is not
/// JSON. A deflate stream zlib refuses, and a text nested deeper than the interpreter reads,
/// end the run.
fn shard_value(shard: &Path) -> Option<Value> {
    let mut bytes = std::fs::read(shard).ok()?;
    if shard.as_os_str().as_bytes().ends_with(b".gz") {
        bytes = match input::gunzip(&bytes) {
            Ok(inflated) => inflated,
            Err(refusal) if refusal.starts_with("Error ") => crate::uncaught!("{refusal}"),
            Err(_) => return None,
        };
    }
    let text = String::from_utf8(bytes).ok()?;
    match pyjson::loads(&text) {
        Ok(parsed) => Some(parsed),
        Err(error) if error.is_value_error() => None,
        Err(error) => crate::uncaught!("{}: {error}", crate::argparse::path_text(shard)),
    }
}

/// The transcript records of one shard: the items of a bare array, or of every array a
/// range-keyed map holds, each one an object with a `stable_id` that reads as true.
fn shard_records(shard: &Path) -> Vec<Value> {
    let Some(parsed) = shard_value(shard) else {
        return Vec::new();
    };
    let items = match parsed.into_dict() {
        Ok(ranges) => {
            let mut items = Vec::new();
            for (_, range) in ranges {
                if !pyjson::truthy(&range) {
                    continue;
                }
                match range.into_list() {
                    Ok(listed) => items.extend(listed),
                    // A text or an object is walked as texts, none of them a record.
                    Err(Value::Str(_) | Value::Dict(_)) => {}
                    Err(other) => {
                        crate::uncaught!("'{}' object is not iterable", python_type(&other))
                    }
                }
            }
            items
        }
        Err(parsed) => parsed.into_list().unwrap_or_default(),
    };
    items
        .into_iter()
        .filter(|tx| tx.get("stable_id").is_some_and(pyjson::truthy))
        .collect()
}

/// `read(chromosome, record)` of every transcript of a candidate JSON cache, in the order
/// of its shards (`shards_in_order`) and of each shard's records, without the ones
/// `read` returns `None` for.
///
/// A shard is filed as `{"chr": [tx, null, ...]}` or as a list; a shard may be plain or
/// gzipped (`.json.gz`); an unreadable shard is skipped. The
/// shards are read on the workers of `cache::threads`, and a run that ends in a shard or in
/// `read` ends on the first such shard in the shards' order.
pub fn read_cache_transcripts<T, F>(cache_dir: &Path, read: F) -> Vec<T>
where
    T: Send,
    F: Fn(&str, Value) -> Option<T> + Sync,
{
    let root = cache::joined(cache_dir, "transcripts");
    if !std::fs::metadata(&root).is_ok_and(|found| found.is_dir()) {
        return Vec::new();
    }
    let shards = shards_in_order(&root);
    let one = |index: usize| -> std::thread::Result<Vec<T>> {
        let (chromosome, shard) = &shards[index];
        panic::catch_unwind(panic::AssertUnwindSafe(|| {
            shard_records(shard)
                .into_iter()
                .filter_map(|tx| read(chromosome, tx))
                .collect()
        }))
    };
    let workers = cache::threads().clamp(1, shards.len().max(1));
    let mut ended: Vec<Option<std::thread::Result<Vec<T>>>> = Vec::new();
    ended.resize_with(shards.len(), || None);
    if workers == 1 {
        for (index, slot) in ended.iter_mut().enumerate() {
            let result = one(index);
            let failed = result.is_err();
            *slot = Some(result);
            if failed {
                break;
            }
        }
    } else {
        let next = AtomicUsize::new(0);
        let first_failure = AtomicUsize::new(usize::MAX);
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|_| {
                    scope.spawn(|| {
                        let mut done = Vec::new();
                        loop {
                            let index = next.fetch_add(1, Ordering::Relaxed);
                            if index >= shards.len()
                                || index > first_failure.load(Ordering::Relaxed)
                            {
                                break;
                            }
                            let result = one(index);
                            if result.is_err() {
                                first_failure.fetch_min(index, Ordering::Relaxed);
                            }
                            done.push((index, result));
                        }
                        done
                    })
                })
                .collect();
            for handle in handles {
                match handle.join() {
                    Ok(done) => {
                        for (index, result) in done {
                            ended[index] = Some(result);
                        }
                    }
                    Err(payload) => panic::resume_unwind(payload),
                }
            }
        });
    }
    let mut read_all = Vec::new();
    for result in ended {
        match result {
            Some(Ok(records)) => read_all.extend(records),
            Some(Err(payload)) => panic::resume_unwind(payload),
            // A shard is passed over only after an earlier one ended the run.
            None => break,
        }
    }
    read_all
}

/// `(chromosome, transcript record)` for every transcript of a candidate JSON cache, in
/// the order `read_cache_transcripts` reads them.
pub fn iter_cache_transcripts(cache_dir: &Path) -> Vec<(String, Value)> {
    read_cache_transcripts(cache_dir, |chromosome, tx| {
        Some((chromosome.to_string(), tx))
    })
}

/// The mate `(chromosome, position)` a bracket-notation breakend allele names.
pub fn bracket_mate(allele: &str) -> Option<(String, Int)> {
    let found = BRACKET_MATE.captures(allele)?;
    Some((
        crate::shared::_strip_chr(&found[1]),
        whole_number(&found[2]),
    ))
}

/// `(chromosome, position)` of a point Location; None for a ranged one.
pub fn point_location(loc: &str) -> Option<(String, Int)> {
    let (chrom, coords) = loc.split_once(':')?;
    if coords.contains('-') || !pyfmt::is_digit(coords) {
        return None;
    }
    Some((crate::shared::_strip_chr(chrom), whole_number(coords)))
}

/// Where a position of a single-exon miRNA with mature ranges lies.
#[derive(Debug, PartialEq, Eq)]
enum MatureMirna {
    Hit,
    Miss,
    MultiExon,
}

/// `Hit` or `Miss` for a single-exon miRNA with mature ranges; `MultiExon` when the cDNA
/// position cannot be read off the span; None for any other transcript.
fn _mature_mirna_hit(pos: Int, model: &TranscriptModel) -> Option<MatureMirna> {
    if model.biotype != "miRNA" || model.mirna_mature.is_empty() {
        return None;
    }
    if model.exons.len() != 1 {
        return Some(MatureMirna::MultiExon);
    }
    let cdna = if model.strand == 1 {
        plus(minus(pos, model.start), 1)
    } else {
        plus(minus(model.end, pos), 1)
    };
    let hit = model
        .mirna_mature
        .iter()
        .any(|(lo, hi)| *lo <= cdna && cdna <= *hi);
    Some(if hit {
        MatureMirna::Hit
    } else {
        MatureMirna::Miss
    })
}

/// Terms a point breakend at `pos`, inside the transcript, carries besides
/// `feature_truncation`: the exon terms (CDS, 5' or 3' UTR by strand, the mature-miRNA
/// term alone, or the non-coding exon term), or the intron terms (`intron_variant`
/// strictly inside the splice sites, the polypyrimidine tract window, and the
/// non-coding transcript term for a transcript without a CDS that is not NMD), plus
/// `NMD_transcript_variant` on a nonsense-mediated-decay transcript.
pub fn region_terms(pos: Int, model: &TranscriptModel) -> BTreeSet<&'static str> {
    let mut terms = BTreeSet::new();
    let exon = model.exons.iter().any(|(s, e)| *s <= pos && pos <= *e);
    if exon {
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
                if _mature_mirna_hit(pos, model) == Some(MatureMirna::Hit) {
                    return BTreeSet::from(["mature_miRNA_variant"]);
                }
                terms.insert("non_coding_transcript_exon_variant");
            }
        }
    } else if let Some((s, e)) = model
        .introns()
        .into_iter()
        .find(|(s, e)| *s <= pos && pos <= *e)
    {
        if plus(s, _SPLICE_SITE) <= pos && pos <= minus(e, _SPLICE_SITE) {
            terms.insert("intron_variant");
        }
        if (model.strand == 1 && minus(e, _PPT_OUTER) <= pos && pos <= minus(e, _PPT_INNER))
            || (model.strand == -1 && plus(s, _PPT_INNER) <= pos && pos <= plus(s, _PPT_OUTER))
        {
            terms.insert("splice_polypyrimidine_tract_variant");
        }
        if model.cds.is_none() && model.biotype != "nonsense_mediated_decay" {
            terms.insert("non_coding_transcript_variant");
        }
    }
    if model.biotype == "nonsense_mediated_decay" {
        terms.insert("NMD_transcript_variant");
    }
    terms
}

/// The term set a point breakend evaluated at `pos` supports on the transcript.
///
/// `mate_inside` says whether the MATE breakend lies inside the transcript, which is
/// what decides `feature_truncation` (Utils/VariationEffect.pm:358 reads the breakend,
/// not the evaluated position). Inside the transcript the region terms apply; within
/// `MAX_DISTANCE_FROM_TRANSCRIPT` outside it the strand-aware upstream or downstream
/// term; beyond that nothing but the truncation term. A mature-miRNA hit stands alone.
pub fn breakend_context(
    pos: Int,
    model: &TranscriptModel,
    mate_inside: bool,
) -> BTreeSet<&'static str> {
    let (start, end) = (model.start, model.end);
    if model.inside(pos) {
        let mut terms = region_terms(pos, model);
        if mate_inside && terms != BTreeSet::from(["mature_miRNA_variant"]) {
            terms.insert("feature_truncation");
        }
        return terms;
    }
    let mut terms = BTreeSet::new();
    if minus(start, MAX_DISTANCE_FROM_TRANSCRIPT) <= pos && pos < start {
        terms.insert(if model.strand == 1 {
            "upstream_gene_variant"
        } else {
            "downstream_gene_variant"
        });
    } else if end < pos && pos <= plus(end, MAX_DISTANCE_FROM_TRANSCRIPT) {
        terms.insert(if model.strand == 1 {
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

/// `(mate_set, local_set)` for a mate-side row, or None when the row is not one.
///
/// The row is a mate-side row when the Location is a point, the Allele is bracket
/// notation, the transcript is on the mate chromosome within the distance gate of the
/// mate, and the local coordinate does not itself select the transcript. `local_set`
/// is Ensembl VEP's read with its empty result replaced by `intergenic_variant`, the
/// default a Transcript row carries when no predicate fires.
pub fn mate_side_derivation(
    loc: &str,
    allele: &str,
    model: &TranscriptModel,
) -> Option<(BTreeSet<&'static str>, BTreeSet<&'static str>)> {
    let local = point_location(loc);
    let mate = bracket_mate(allele);
    let ((lchr, lpos), (mchr, mpos)) = (local?, mate?);
    if model.chrom != mchr || !model.near(mpos) {
        return None;
    }
    if lchr == model.chrom && model.near(lpos) {
        return None;
    }
    let mate_inside = model.inside(mpos);
    let mate_set = breakend_context(mpos, model, mate_inside);
    let mut local_set = breakend_context(lpos, model, mate_inside);
    if local_set.is_empty() {
        local_set.insert("intergenic_variant");
    }
    Some((mate_set, local_set))
}
