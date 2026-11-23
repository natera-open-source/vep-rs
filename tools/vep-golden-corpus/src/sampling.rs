// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The draws of `select`: the Mersenne Twister as `random.Random(seed)` seeds and reads it,
//! `randrange`, `shuffle` and the reservoir step that keeps a bounded sample of a stream;
//! the pass over the reference outputs that samples their rows, and the choice of the
//! exemplars among the sampled rows.

pub mod random;

pub use random::Random;

use crate::cache_cost::CacheCost;
use crate::constants::{FAMILY_RESERVOIR, RESERVOIR};
use crate::error::Failure;
use crate::rows::blocks::{lines, read_blocks, RowReader};
use crate::rows::{substitution_alts, Dataset, Row};
use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;

/// One reservoir step under a capacity of any sign: `bucket` stays a uniform sample of at
/// most `capacity` of the `offered` items seen so far for its key. A bucket below its
/// capacity takes the item; a full one draws a place among the items offered and takes the
/// item there when the place is inside the bucket. `item` is asked for only when it is kept.
pub fn reservoir_step_with<T>(
    rng: &mut Random,
    bucket: &mut Vec<T>,
    offered: usize,
    capacity: i64,
    item: impl FnOnce() -> T,
) {
    if i64::try_from(bucket.len()).is_ok_and(|held| held < capacity) {
        bucket.push(item());
        return;
    }
    let place = rng.randrange(offered as u64);
    if i64::try_from(place).is_ok_and(|place| place < capacity) {
        bucket[place as usize] = item();
    }
}

/// One reservoir step (`reservoir_step`): `bucket` stays a uniform sample of at most
/// `capacity` of the `offered` items seen so far for its key.
pub fn reservoir_step<T>(
    rng: &mut Random,
    bucket: &mut Vec<T>,
    item: T,
    offered: usize,
    capacity: usize,
) {
    let capacity = i64::try_from(capacity).unwrap_or(i64::MAX);
    reservoir_step_with(rng, bucket, offered, capacity, || item);
}

/// A map that keeps its keys in the order they were first set, as a `dict` does.
#[derive(Debug, Clone)]
pub struct OrderedMap<K, V> {
    entries: Vec<(K, V)>,
    places: HashMap<K, usize>,
}

impl<K, V> Default for OrderedMap<K, V> {
    fn default() -> Self {
        OrderedMap {
            entries: Vec::new(),
            places: HashMap::new(),
        }
    }
}

impl<K: Hash + Eq + Clone, V> OrderedMap<K, V> {
    pub fn new() -> Self {
        OrderedMap::default()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The place of `key` among the keys, `None` for a key the map lacks.
    pub fn place<Q: Hash + Eq + ?Sized>(&self, key: &Q) -> Option<usize>
    where
        K: Borrow<Q>,
    {
        self.places.get(key).copied()
    }

    pub fn contains_key<Q: Hash + Eq + ?Sized>(&self, key: &Q) -> bool
    where
        K: Borrow<Q>,
    {
        self.places.contains_key(key)
    }

    /// `d.get(key)`.
    pub fn get<Q: Hash + Eq + ?Sized>(&self, key: &Q) -> Option<&V>
    where
        K: Borrow<Q>,
    {
        self.place(key).map(|place| &self.entries[place].1)
    }

    /// The place of `key`, set to `value()` after the keys the map holds when it lacks the
    /// key: `d.setdefault(key, value)`.
    pub fn place_or_insert_with<Q>(&mut self, key: &Q, value: impl FnOnce() -> V) -> usize
    where
        K: Borrow<Q>,
        Q: Hash + Eq + ?Sized + ToOwned<Owned = K>,
    {
        if let Some(place) = self.place(key) {
            return place;
        }
        self.places.insert(key.to_owned(), self.entries.len());
        self.entries.push((key.to_owned(), value()));
        self.entries.len() - 1
    }

    /// `d.setdefault(key, V())`, open to change.
    pub fn entry<Q>(&mut self, key: &Q) -> &mut V
    where
        K: Borrow<Q>,
        Q: Hash + Eq + ?Sized + ToOwned<Owned = K>,
        V: Default,
    {
        let place = self.place_or_insert_with(key, V::default);
        &mut self.entries[place].1
    }

    /// `d[key] = value`: a key the map holds keeps its place.
    pub fn insert(&mut self, key: K, value: V) {
        match self.places.get(&key) {
            Some(place) => self.entries[*place].1 = value,
            None => {
                self.places.insert(key.clone(), self.entries.len());
                self.entries.push((key, value));
            }
        }
    }

    /// The key and the value at a place.
    pub fn at(&self, place: usize) -> (&K, &V) {
        let (key, value) = &self.entries[place];
        (key, value)
    }

    /// The value at a place, open to change.
    pub fn at_mut(&mut self, place: usize) -> &mut V {
        &mut self.entries[place].1
    }

    /// `d.items()`.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.entries.iter().map(|(key, value)| (key, value))
    }

    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.entries.iter().map(|(key, _)| key)
    }

    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.entries.iter().map(|(_, value)| value)
    }

    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut V> {
        self.entries.iter_mut().map(|(_, value)| value)
    }

    /// Keeps the entries at the places `keep` answers true for, in their order.
    pub fn retain_places(&mut self, keep: impl Fn(usize) -> bool) {
        let entries = std::mem::take(&mut self.entries);
        self.places.clear();
        for (place, (key, value)) in entries.into_iter().enumerate() {
            if keep(place) {
                self.places.insert(key.clone(), self.entries.len());
                self.entries.push((key, value));
            }
        }
    }
}

/// A bounded uniform sample of the rows of one key, and how many rows the key showed.
#[derive(Debug, Clone, Default)]
pub struct Sample {
    pub seen: usize,
    pub rows: Vec<Row>,
}

/// The rows and the consequence sets one dataset showed.
#[derive(Debug, Clone, Default)]
pub struct DatasetCount {
    pub rows: usize,
    /// The places of its consequence sets in the pass's `reservoir`.
    pub combinations: HashSet<usize>,
}

/// What one pass over the references yields (`reservoir_pass`).
#[derive(Debug, Clone, Default)]
pub struct ReservoirPass {
    /// Per consequence set, in the order the sets were first seen: the count of its rows
    /// and a sample of at most `RESERVOIR` of them.
    pub reservoir: OrderedMap<String, Sample>,
    /// Per dataset id, the rows and the consequence sets of its reference.
    pub per_dataset: HashMap<String, DatasetCount>,
    /// Per family stratum, with family specs: the count of its rows and a sample of at
    /// most `FAMILY_RESERVOIR` of them, so that a stratum too rare to reach a consequence
    /// set's sample has candidates of its own.
    pub families: OrderedMap<String, Sample>,
    /// Every row of the named records, by dataset id and Uploaded_variation, which the
    /// samples' bounds never touch.
    pub named_rows: OrderedMap<(String, String), Vec<Row>>,
}

/// What a worker reads of one row: where its line stands in the block, the block's own
/// numbers for its consequence set and its family stratum, and whether its name is wanted.
struct RowFact {
    start: usize,
    length: usize,
    set: u32,
    stratum: u32,
    named: bool,
}

/// What a worker reads of one block: its rows, and the distinct consequence sets and
/// family strata in the order the block first shows them.
#[derive(Default)]
struct BlockFacts {
    rows: Vec<RowFact>,
    sets: Vec<String>,
    strata: Vec<String>,
}

/// The number of `name` among `names`, added behind them when it is not there.
fn number_of(numbers: &mut HashMap<String, u32>, names: &mut Vec<String>, name: &str) -> u32 {
    if let Some(number) = numbers.get(name) {
        return *number;
    }
    let number = names.len() as u32;
    numbers.insert(name.to_string(), number);
    names.push(name.to_string());
    number
}

fn block_facts(
    block: &str,
    reader: &RowReader<'_>,
    with_strata: bool,
    wanted: Option<&HashSet<String>>,
) -> BlockFacts {
    let mut facts = BlockFacts::default();
    let mut set_numbers: HashMap<String, u32> = HashMap::new();
    let mut stratum_numbers: HashMap<String, u32> = HashMap::new();
    for line in lines(block) {
        let Some(view) = reader.view(line) else {
            continue;
        };
        // Rows of one record follow each other and often share their consequence set.
        let set = match facts.rows.last() {
            Some(last) if facts.sets[last.set as usize] == view.consequence_set => last.set,
            _ => number_of(&mut set_numbers, &mut facts.sets, &view.consequence_set),
        };
        let stratum = if with_strata {
            let stratum = view.family_stratum();
            number_of(&mut stratum_numbers, &mut facts.strata, &stratum)
        } else {
            0
        };
        facts.rows.push(RowFact {
            start: line.as_ptr() as usize - block.as_ptr() as usize,
            length: line.len(),
            set,
            stratum,
            named: wanted.is_some_and(|names| names.contains(view.uploaded_variation)),
        });
    }
    facts
}

/// One streaming pass over the references (`reservoir_pass`): per consequence set a bounded
/// uniform sample of rows plus counts, per family stratum (when `family_specs` are given) a
/// smaller sample of its own, and every row of the records `named` (per dataset id, their
/// Uploaded_variation names). The draws are those of `random.Random(seed)` over the rows in
/// file order; the rows are read on at most `threads` threads.
pub fn reservoir_pass(
    datasets: &[Dataset],
    seed: i64,
    family_specs: &[String],
    named: &HashMap<String, HashSet<String>>,
    threads: usize,
) -> Result<ReservoirPass, Failure> {
    let mut rng = Random::new(seed);
    let mut pass = ReservoirPass::default();
    let with_strata = !family_specs.is_empty();
    let (reservoir, families) = (RESERVOIR as i64, FAMILY_RESERVOIR as i64);
    for (sid, reference, input) in datasets {
        let wanted = named.get(sid).filter(|names| !names.is_empty());
        let substitutions = substitution_alts(input)?;
        let reader = RowReader::new(sid, family_specs, Some(&substitutions));
        let read = |block: &str| block_facts(block, &reader, with_strata, wanted);
        let take = |block: &str, facts: BlockFacts| {
            let sets: Vec<usize> = facts
                .sets
                .iter()
                .map(|set| {
                    pass.reservoir
                        .place_or_insert_with(set.as_str(), Sample::default)
                })
                .collect();
            let strata: Vec<usize> = facts
                .strata
                .iter()
                .map(|stratum| {
                    pass.families
                        .place_or_insert_with(stratum.as_str(), Sample::default)
                })
                .collect();
            let count = pass.per_dataset.entry(sid.clone()).or_default();
            count.rows += facts.rows.len();
            count.combinations.extend(&sets);
            for fact in &facts.rows {
                let row = || {
                    let line = &block[fact.start..fact.start + fact.length];
                    reader.view(line).expect("a row").to_row()
                };
                let sample = pass.reservoir.at_mut(sets[fact.set as usize]);
                sample.seen += 1;
                reservoir_step_with(&mut rng, &mut sample.rows, sample.seen, reservoir, row);
                if with_strata {
                    let sample = pass.families.at_mut(strata[fact.stratum as usize]);
                    sample.seen += 1;
                    reservoir_step_with(&mut rng, &mut sample.rows, sample.seen, families, row);
                }
                if fact.named {
                    let row = row();
                    let key = (sid.clone(), row.uploaded_variation.clone());
                    pass.named_rows.entry(&key).push(row);
                }
            }
            Ok(())
        };
        read_blocks(reference, threads, read, take)?;
    }
    Ok(pass)
}

/// The key candidates are sorted by: whether the row is wider than `max_span`, the bytes it
/// adds to the fixture cache, and its preference, negated.
pub type ExemplarOrder = (u8, usize, i64);

/// The sort key of a candidate row (`exemplar_order`): no wider than `max_span` (when it is
/// not 0) first, then the fewest cache bytes under `cost`, then the highest preference.
pub fn exemplar_order(
    row: &Row,
    max_span: i64,
    cost: Option<&CacheCost>,
) -> Result<ExemplarOrder, Failure> {
    let wide = u8::from(max_span != 0 && row.span() > max_span);
    let bytes = match cost {
        Some(cost) => cost.row_cost(row)?,
        None => 0,
    };
    Ok((wide, bytes, -row.preference()))
}

/// The rows chosen as exemplars, by consequence set, in the order the sets were first given
/// rows.
pub type Chosen = OrderedMap<String, Vec<Row>>;

/// Up to `k` rows a combination (`choose_exemplars`), distinct (class, strand, family
/// classes) first, in `exemplar_order`. The candidates of a combination are shuffled with
/// `rng` before they are sorted, so that rows the order does not tell apart fall in random
/// order. A combination that holds any of `focus_terms` takes up to `focus_k` rows in place
/// of `k`, when `focus_k` is given. A count of zero or less takes every candidate.
pub fn choose_exemplars(
    reservoir: &OrderedMap<String, Sample>,
    k: i64,
    rng: &mut Random,
    max_span: i64,
    focus_terms: &[String],
    focus_k: Option<i64>,
    cost: Option<&CacheCost>,
) -> Result<Chosen, Failure> {
    let focus: HashSet<&str> = focus_terms.iter().map(String::as_str).collect();
    let mut combos: Vec<&String> = reservoir.keys().collect();
    combos.sort();
    let mut chosen = Chosen::new();
    for combo in combos {
        let focused = combo.split(',').any(|term| focus.contains(term));
        let want = match focus_k {
            Some(focus_k) if focused => focus_k,
            _ => k,
        };
        let full = |picked: &Vec<&Row>| i64::try_from(picked.len()).is_ok_and(|n| n == want);
        let sample = reservoir.get(combo.as_str()).expect("a key of the map");
        let mut candidates: Vec<&Row> = sample.rows.iter().collect();
        rng.shuffle(&mut candidates);
        let mut ordered = Vec::with_capacity(candidates.len());
        for row in candidates {
            ordered.push((exemplar_order(row, max_span, cost)?, row));
        }
        ordered.sort_by_key(|(order, _)| *order);
        let mut picked: Vec<&Row> = Vec::new();
        let mut strata_seen: HashSet<(String, &str, &[String])> = HashSet::new();
        for (_, row) in &ordered {
            if !strata_seen.insert((row.klass(), row.strand.as_str(), row.family.as_slice())) {
                continue;
            }
            picked.push(row);
            if full(&picked) {
                break;
            }
        }
        for (_, row) in &ordered {
            if full(&picked) {
                break;
            }
            if !picked.contains(row) {
                picked.push(row);
            }
        }
        chosen.insert(combo.clone(), picked.into_iter().cloned().collect());
    }
    Ok(chosen)
}

/// Every variant class the references show keeps an exemplar (`ensure_class_exemplars`): a
/// class no chosen row has takes, from the whole reservoir, the row that sorts first in
/// `exemplar_order`, added to the combination it was sampled for.
pub fn ensure_class_exemplars(
    chosen: &mut Chosen,
    reservoir: &OrderedMap<String, Sample>,
    max_span: i64,
    cost: Option<&CacheCost>,
) -> Result<(), Failure> {
    let covered: HashSet<String> = chosen.values().flatten().map(Row::klass).collect();
    let mut best: OrderedMap<String, (ExemplarOrder, &String, &Row)> = OrderedMap::new();
    for (combo, sample) in reservoir.iter() {
        for row in &sample.rows {
            let klass = row.klass();
            if covered.contains(&klass) {
                continue;
            }
            let order = exemplar_order(row, max_span, cost)?;
            if best.get(&klass).is_none_or(|held| order < held.0) {
                best.insert(klass, (order, combo, row));
            }
        }
    }
    for (_, combo, row) in best.values() {
        chosen.entry(combo.as_str()).push((*row).clone());
    }
    Ok(())
}

/// Every family stratum the references show keeps an exemplar (`ensure_family_exemplars`):
/// a stratum no chosen row has takes, from its own sample, the row that sorts first in
/// `exemplar_order`, added to the combination it was sampled for.
pub fn ensure_family_exemplars(
    chosen: &mut Chosen,
    families: &OrderedMap<String, Sample>,
    max_span: i64,
    cost: Option<&CacheCost>,
) -> Result<(), Failure> {
    let covered: HashSet<String> = chosen.values().flatten().map(Row::family_stratum).collect();
    let mut strata: Vec<(&String, &Sample)> = families.iter().collect();
    strata.sort_by_key(|(stratum, _)| *stratum);
    for (stratum, sample) in strata {
        if covered.contains(stratum) {
            continue;
        }
        let mut first: Option<(ExemplarOrder, &Row)> = None;
        for row in &sample.rows {
            let order = exemplar_order(row, max_span, cost)?;
            if first.is_none_or(|(held, _)| order < held) {
                first = Some((order, row));
            }
        }
        if let Some((_, row)) = first {
            chosen.entry(row.consequence_set.as_str()).push(row.clone());
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "sampling_tests.rs"]
mod tests;
