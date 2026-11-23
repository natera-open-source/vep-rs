// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The settling of one engine's remainder and the counts it gives each source.

use super::api::{Counted, Extracted, NamedLines, PluginReads, PROJECTIONS};
use super::class_table::Classes;
use super::context::RowContext;
use super::error::{Scored, ScorerError};
use super::merge::merge_spans;
use super::projection::{AlleleScore, FeatureScore, Verdicts};
use super::ranges::{Ranges, RANGE_BYTES};
use super::residuals::{ResidualLines, Residuals};
use super::sources::Sources;
use super::tally::Tally;
use crate::ordered::in_order;
use crate::rules::{Reads, Rules};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// The counts of one source, engine and projection under the report's count key names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceEntry {
    /// The reference's tuples.
    pub perl: i64,
    /// The engine's tuples.
    pub rust: i64,
    pub intersection: i64,
    /// `adj_perl`, `adj_rust` and `adj_intersection`; `None` where a class that is not
    /// applied would set a tuple aside.
    pub adjusted: Option<[i64; 3]>,
    /// Every applied class of the source in table order, with the keys under which it set
    /// a tuple aside.
    pub masked: Vec<(String, i64)>,
    /// The classes not applied that would set a tuple aside, sorted, with their keys.
    pub not_applied: Vec<(String, i64)>,
    pub text_only_differences: i64,
    pub residual_rows: i64,
    pub shared_keys: i64,
}

fn count(value: u64) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

/// How the settling step of one engine is spread. The tallies and the residual files are
/// the same whatever the spread.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Spread {
    /// The threads that settle ranges of Locations at a time.
    pub threads: usize,
    /// The bytes of one sorted file a range holds, a Location's lines kept whole.
    pub range_bytes: u64,
}

impl Spread {
    /// Over `threads` threads, in ranges of the program's own size.
    pub fn over(threads: usize) -> Spread {
        Spread {
            threads,
            range_bytes: RANGE_BYTES,
        }
    }
}

/// The steps of one range, in the order the steps of a walk over every Location at once
/// would come: every range's feature projection, the sorting of the verdict lines, every
/// range's allele projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Step {
    Feature,
    Sort,
    Allele,
}

/// The refusals met so far, by which a range knows the steps it can leave out: none that
/// could raise the run's refusal, the first in the order of `Step` and then of the ranges.
struct Refused {
    /// The least range whose feature projection refused.
    feature: AtomicUsize,
    sort: AtomicBool,
    /// The least range whose allele projection refused.
    allele: AtomicUsize,
}

impl Refused {
    fn none() -> Refused {
        Refused {
            feature: AtomicUsize::new(usize::MAX),
            sort: AtomicBool::new(false),
            allele: AtomicUsize::new(usize::MAX),
        }
    }

    fn note(&self, step: Step, range: usize) {
        match step {
            Step::Feature => {
                self.feature.fetch_min(range, Ordering::Relaxed);
            }
            Step::Sort => self.sort.store(true, Ordering::Relaxed),
            Step::Allele => {
                self.allele.fetch_min(range, Ordering::Relaxed);
            }
        }
    }

    /// Whether a refusal met comes before any the feature projection of `range` could
    /// raise.
    fn before_feature(&self, range: usize) -> bool {
        self.feature.load(Ordering::Relaxed) < range
    }

    /// Whether a refusal met comes before any the allele projection of `range` could
    /// raise.
    fn before_allele(&self, range: usize) -> bool {
        self.feature.load(Ordering::Relaxed) != usize::MAX
            || self.sort.load(Ordering::Relaxed)
            || self.allele.load(Ordering::Relaxed) < range
    }
}

/// What one range leaves: per projection its tally and, where the run keeps them, its
/// residual lines.
struct Settled {
    tallies: [Tally; 2],
    residuals: [Option<ResidualLines>; 2],
}

enum Outcome {
    Settled(Box<Settled>),
    Refused(Step, ScorerError),
    /// Left out after a refusal that comes first.
    Skipped,
}

/// A file removed when the value is dropped.
struct Removed<'p>(&'p Path);

impl Drop for Removed<'_> {
    fn drop(&mut self) {
        let _ = fs::remove_file(self.0);
    }
}

/// One engine against the reference, as every range reads it.
struct Engine<'a> {
    reference: &'a Extracted,
    engine: &'a Extracted,
    verdicts: &'a Path,
    sources: &'a Sources<'a>,
    classes: &'a Classes,
    rules: &'a dyn Rules,
    reads: &'a Reads,
    plugin_reads: &'a PluginReads,
    run_input: Option<&'a dyn NamedLines>,
    /// Whether the ranges keep their residual lines.
    residuals: bool,
    ranges: Ranges,
    refused: Refused,
}

/// The places of the two sides' files in the order `Ranges` is given them.
const FEATURE: [usize; 2] = [0, 1];
const ALLELE: [usize; 2] = [2, 3];
const ROWS: [usize; 2] = [4, 5];

impl Engine<'_> {
    /// The file of the sorted verdict lines of one range.
    fn verdicts_of(&self, range: usize) -> PathBuf {
        let mut name = self.verdicts.as_os_str().to_os_string();
        name.push(format!(".{range}"));
        PathBuf::from(name)
    }

    /// Settles one range: its feature projection, whose verdict lines are sorted into a
    /// file of the range's own, then its allele projection, which reads them.
    fn settle(&self, range: usize) -> Outcome {
        match self.settled(range) {
            Ok(Some(settled)) => Outcome::Settled(Box::new(settled)),
            Ok(None) => Outcome::Skipped,
            Err((step, refusal)) => {
                self.refused.note(step, range);
                Outcome::Refused(step, refusal)
            }
        }
    }

    fn settled(&self, range: usize) -> Result<Option<Settled>, (Step, ScorerError)> {
        if self.refused.before_feature(range) {
            return Ok(None);
        }
        let at = |step: Step| move |refusal: ScorerError| (step, refusal);
        let spans = self.ranges.spans(range).map_err(at(Step::Feature))?;
        let context = RowContext::new(
            self.reference,
            self.engine,
            self.sources,
            self.classes,
            self.rules,
            self.reads,
            self.plugin_reads,
            self.run_input,
        )
        .within(spans[ROWS[0]], spans[ROWS[1]]);
        let verdicts = self.verdicts_of(range);
        let mut feature = FeatureScore::new(
            context,
            Verdicts::new(&verdicts),
            self.sources,
            self.classes,
            self.residuals,
        );
        merge_spans(
            (&self.reference.feature, spans[FEATURE[0]]),
            (&self.engine.feature, spans[FEATURE[1]]),
            &mut feature,
        )
        .map_err(at(Step::Feature))?;
        let (feature_tally, feature_lines, unsorted) = feature.done();
        let _removed = Removed(&verdicts);
        unsorted.finish().map_err(at(Step::Sort))?;
        if self.refused.before_allele(range) {
            return Ok(None);
        }
        let mut allele = AlleleScore::new(&verdicts, self.sources, self.classes, self.residuals)
            .map_err(at(Step::Allele))?;
        merge_spans(
            (&self.reference.allele, spans[ALLELE[0]]),
            (&self.engine.allele, spans[ALLELE[1]]),
            &mut allele,
        )
        .map_err(at(Step::Allele))?;
        let (allele_tally, allele_lines) = allele.done();
        Ok(Some(Settled {
            tallies: [feature_tally, allele_tally],
            residuals: [feature_lines, allele_lines],
        }))
    }
}

/// What the ranges leave, joined in their order: the tallies summed, the residual lines
/// appended to the engine's files, and the first refusal in the order of `Step` and then of
/// the ranges.
struct Joined<'r> {
    tallies: [Tally; 2],
    residuals: Option<&'r mut Residuals>,
    refused: Option<(Step, ScorerError)>,
}

impl Joined<'_> {
    fn refuse(&mut self, step: Step, refusal: ScorerError) {
        if self.refused.as_ref().is_none_or(|(held, _)| step < *held) {
            self.refused = Some((step, refusal));
        }
    }

    /// Takes the outcome of the next range.
    fn take(&mut self, outcome: Outcome) {
        let settled = match outcome {
            Outcome::Settled(settled) if self.refused.is_none() => *settled,
            Outcome::Settled(_) | Outcome::Skipped => return,
            Outcome::Refused(step, refusal) => return self.refuse(step, refusal),
        };
        let steps = [Step::Feature, Step::Allele];
        for (place, (tally, lines)) in settled
            .tallies
            .into_iter()
            .zip(settled.residuals)
            .enumerate()
        {
            self.tallies[place].add(tally);
            let (Some(files), Some(lines)) = (self.residuals.as_deref_mut(), lines) else {
                continue;
            };
            for (source, lines) in lines.by_source() {
                if let Err(refusal) = files.append(PROJECTIONS[place], source, lines) {
                    return self.refuse(steps[place], refusal);
                }
            }
        }
    }

    fn done(self) -> Scored<[Tally; 2]> {
        match self.refused {
            Some((_, refusal)) => Err(refusal),
            None => Ok(self.tallies),
        }
    }
}

/// The tallies of one engine against the reference, the feature projection's and the
/// allele projection's. The tuple files of `reference` and `engine` are the remainder the
/// passes left; every line of them goes to the settling step, and what the passes
/// `counted` is added to what that step tallies. The sorted files are cut into ranges of
/// whole Locations, which `spread.threads` threads settle, each range with readers and
/// rows of its own; what the ranges leave is joined in their order, so the tallies and
/// the residual files are the same whatever the spread, and a run that a range refuses
/// ends in the refusal a walk over every Location in order meets first. A range's verdict
/// lines are written beside `verdicts`, in a file of the range's own that its allele
/// projection reads and that is removed with the range. The shards of a JSON cache that a
/// class asks `reads` for are read on as many threads. With `residuals`, the tuples left
/// unexplained are written to the engine's residual files.
#[allow(clippy::too_many_arguments)]
pub fn score_engine(
    reference: &Extracted,
    engine: &Extracted,
    counted: &Counted,
    verdicts: &Path,
    sources: &Sources,
    classes: &Classes,
    rules: &dyn Rules,
    reads: &Reads,
    plugin_reads: &PluginReads,
    residuals: Option<&mut Residuals>,
    run_input: Option<&dyn NamedLines>,
    spread: Spread,
) -> Scored<[Tally; 2]> {
    let sorted = [
        &reference.feature,
        &engine.feature,
        &reference.allele,
        &engine.allele,
        &reference.rows,
        &engine.rows,
    ];
    let engine = Engine {
        reference,
        engine,
        verdicts,
        sources,
        classes,
        rules,
        reads,
        plugin_reads,
        run_input,
        residuals: residuals.is_some(),
        ranges: Ranges::new(&sorted.map(PathBuf::as_path), spread.range_bytes)?,
        refused: Refused::none(),
    };
    let mut joined = Joined {
        tallies: [Tally::default(), Tally::default()],
        residuals,
        refused: None,
    };
    reads.read_on(spread.threads);
    in_order(
        engine.ranges.count(),
        spread.threads,
        &|range| engine.settle(range),
        &mut |outcome| {
            joined.take(outcome);
            true
        },
    );
    let [mut feature_tally, mut allele_tally] = joined.done()?;
    feature_tally.take(&counted.feature)?;
    allele_tally.take(&counted.allele)?;
    Ok([feature_tally, allele_tally])
}

/// The tally summed per source, the sources in their names' order. `masked` lists every
/// applied class of the source, and an entry in which a class that is not applied would set
/// a tuple aside carries no adjusted count.
pub fn source_entries(
    tally: &Tally,
    sources: &Sources,
    classes: &Classes,
) -> Scored<Vec<(String, SourceEntry)>> {
    let mut sums: BTreeMap<String, [u64; 4]> = BTreeMap::new();
    let held =
        |counts: &BTreeMap<Vec<u8>, u64>, name: &[u8]| counts.get(name).copied().unwrap_or(0);
    for name in tally.fields() {
        let counts = sums.entry(sources.named(&name)?).or_default();
        counts[0] += held(&tally.paired, &name);
        counts[1] += held(&tally.ref_left, &name);
        counts[2] += held(&tally.eng_left, &name);
        counts[3] += held(&tally.text_only, &name);
    }
    let mut entries = Vec::new();
    for (source, [paired, ref_left, eng_left, text_only]) in sums {
        let aside =
            |counts: &BTreeMap<String, u64>| count(counts.get(&source).copied().unwrap_or(0));
        let (perl, rust) = (count(paired + ref_left), count(paired + eng_left));
        let adj_perl = perl - aside(&tally.ref_aside);
        let adj_rust = rust - aside(&tally.eng_aside);
        let not_applied: Vec<(String, i64)> = tally
            .not_applied
            .get(&source)
            .into_iter()
            .flatten()
            .map(|(class, keys)| (class.clone(), count(*keys)))
            .collect();
        let masked = tally.masked.get(&source);
        let entry = SourceEntry {
            perl,
            rust,
            intersection: count(paired),
            adjusted: not_applied
                .is_empty()
                .then_some([adj_perl, adj_rust, count(paired)]),
            masked: classes
                .of(&source)
                .into_iter()
                .map(|class| {
                    let keys = masked.and_then(|masked| masked.get(&class)).copied();
                    (class, count(keys.unwrap_or(0)))
                })
                .collect(),
            not_applied,
            text_only_differences: count(text_only),
            residual_rows: adj_perl + adj_rust - 2 * count(paired),
            shared_keys: aside(&tally.shared),
        };
        entries.push((source, entry));
    }
    Ok(entries)
}
