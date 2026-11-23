// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript comparison: the records of a shard under any of its containers, what is
//! kept of a record, the fields two records of one stable id differ in, and the counts of
//! one shard name both sides are asked for.

use crate::input;
use crate::json::{py_ne, truthy, Reader, Res, Shape, Value};
use crate::paths::Pair;
use std::collections::HashMap;
use std::path::Path;

/// The members of a record compared as they stand.
pub const COMPARE_FIELDS: [&str; 9] = [
    "start",
    "end",
    "strand",
    "biotype",
    "source",
    "gene_stable_id",
    "gene_symbol",
    "hgnc_id",
    "canonical",
];

/// The members of a record's feature cache compared as they stand.
pub const VEFC_FIELDS: [&str; 7] = [
    "introns",
    "sorted_exons",
    "five_prime_utr",
    "three_prime_utr",
    "codon_table",
    "protein_features",
    "seq_edits",
];

const VEFC_INTRONS: usize = 0;
const VEFC_SORTED_EXONS: usize = 1;

/// The members counted where the feature cache lacks them, in the order of `top`.
const COUNTED: [&str; 3] = ["sorted_exons", "introns", "mapper"];

/// The names a mismatch is counted under, in the order a record's comparison meets them:
/// the nine members, the three counts, then for each member of the feature cache its
/// `_count` name and its plain name, of which a record takes one.
pub const MISMATCH_NAMES: usize = COMPARE_FIELDS.len() + COUNTED.len() + 2 * VEFC_FIELDS.len();

const EXON_COUNT: usize = COMPARE_FIELDS.len();
const INTRON_COUNT: usize = EXON_COUNT + 1;
const MAPPER_PAIR_COUNT: usize = EXON_COUNT + 2;
const VEFC_FIRST: usize = EXON_COUNT + 3;

/// The name mismatch `rank` is counted under.
pub fn mismatch_name(rank: usize) -> String {
    match rank {
        EXON_COUNT => "exon_count".to_string(),
        INTRON_COUNT => "intron_count".to_string(),
        MAPPER_PAIR_COUNT => "mapper_pair_count".to_string(),
        _ if rank < EXON_COUNT => COMPARE_FIELDS[rank].to_string(),
        _ => {
            let field = VEFC_FIELDS[(rank - VEFC_FIRST) / 2];
            if (rank - VEFC_FIRST).is_multiple_of(2) {
                format!("vefc_{field}_count")
            } else {
                format!("vefc_{field}")
            }
        }
    }
}

/// What is kept of a record's feature cache: its seven compared members, `None` for one
/// it lacks, and the shape of its `mapper`.
struct FeatureCache {
    fields: [Value; 7],
    mapper: Shape,
}

impl FeatureCache {
    fn empty() -> FeatureCache {
        FeatureCache {
            fields: [const { Value::Null }; 7],
            mapper: Shape::ABSENT,
        }
    }
}

/// A member that should hold a feature cache: an object, or any other value.
enum CacheMember {
    Object {
        cache: Box<FeatureCache>,
        truthy: bool,
    },
    Other {
        truthy: bool,
    },
}

/// What is kept of a record that carries an id.
pub struct Transcript {
    pub id: Box<[u8]>,
    fields: [Value; 9],
    /// The shapes of the record's own `sorted_exons`, `introns` and `mapper`.
    top: [Shape; 3],
    /// The feature cache the comparison reads, `None` where the member it would read is
    /// no object.
    cache: Option<Box<FeatureCache>>,
}

/// One object in a record's place.
enum Record {
    /// It carries no id and is passed over.
    Nameless,
    /// Its id is a value Python counts true that is no string.
    UnreadableId,
    Named(Box<Transcript>),
}

/// A member that may hold an id: a string, or any other value with its truth.
enum IdMember {
    Text(Vec<u8>),
    Other { truthy: bool },
}

fn id_member(reader: &mut Reader) -> Res<IdMember> {
    if reader.value_start()? == b'"' {
        let mut text = Vec::new();
        reader.string(Some(&mut text))?;
        Ok(IdMember::Text(text))
    } else {
        Ok(IdMember::Other {
            truthy: reader.shape()?.truthy,
        })
    }
}

fn cache_member(reader: &mut Reader) -> Res<CacheMember> {
    if reader.value_start()? != b'{' {
        return Ok(CacheMember::Other {
            truthy: reader.shape()?.truthy,
        });
    }
    let mut cache = Box::new(FeatureCache::empty());
    let mut any = false;
    if reader.enter(b'}')? {
        let mut name = Vec::new();
        loop {
            any = true;
            name.clear();
            reader.key(&mut name)?;
            if name == b"mapper" {
                cache.mapper = reader.shape()?;
            } else if let Some(field) = VEFC_FIELDS.iter().position(|f| f.as_bytes() == name) {
                cache.fields[field] = reader.value()?;
            } else {
                reader.skip_value()?;
            }
            if !reader.another(b'}')? {
                break;
            }
        }
    }
    Ok(CacheMember::Object { cache, truthy: any })
}

/// Reads the object the reader stands at as a record.
fn record(reader: &mut Reader) -> Res<Record> {
    let mut stable_id = None;
    let mut plain_id = None;
    let mut fields = [const { Value::Null }; 9];
    let mut top = [Shape::ABSENT; 3];
    let mut long_cache = None;
    let mut short_cache = None;
    if reader.enter(b'}')? {
        let mut name = Vec::new();
        loop {
            name.clear();
            reader.key(&mut name)?;
            if name == b"stable_id" {
                stable_id = Some(id_member(reader)?);
            } else if name == b"id" {
                plain_id = Some(id_member(reader)?);
            } else if name == b"variation_effect_feature_cache" {
                long_cache = Some(cache_member(reader)?);
            } else if name == b"vefc" {
                short_cache = Some(cache_member(reader)?);
            } else if let Some(field) = COMPARE_FIELDS.iter().position(|f| f.as_bytes() == name) {
                fields[field] = reader.value()?;
            } else if let Some(counted) = COUNTED.iter().position(|f| f.as_bytes() == name) {
                top[counted] = reader.shape()?;
            } else {
                reader.skip_value()?;
            }
            if !reader.another(b'}')? {
                break;
            }
        }
    }
    // `record.get("stable_id") or record.get("id", "")`, kept where Python counts it true.
    let id = match stable_id {
        Some(IdMember::Text(text)) if !text.is_empty() => text,
        Some(IdMember::Other { truthy: true }) => return Ok(Record::UnreadableId),
        _ => match plain_id {
            Some(IdMember::Text(text)) if !text.is_empty() => text,
            Some(IdMember::Other { truthy: true }) => return Ok(Record::UnreadableId),
            _ => return Ok(Record::Nameless),
        },
    };
    // `record.get("variation_effect_feature_cache") or record.get("vefc", {})`.
    let chosen = match long_cache {
        Some(CacheMember::Object {
            cache,
            truthy: true,
        }) => Some(cache),
        Some(CacheMember::Other { truthy: true }) => None,
        _ => match short_cache {
            None => Some(Box::new(FeatureCache::empty())),
            Some(CacheMember::Object { cache, .. }) => Some(cache),
            Some(CacheMember::Other { .. }) => None,
        },
    };
    Ok(Record::Named(Box::new(Transcript {
        id: id.into(),
        fields,
        top,
        cache: chosen,
    })))
}

/// The records among the items of the array the reader stands at: its objects.
fn records_of_array(reader: &mut Reader, records: &mut Vec<Record>) -> Res<()> {
    if reader.enter(b']')? {
        loop {
            if reader.value_start()? == b'{' {
                records.push(record(reader)?);
            } else {
                reader.skip_value()?;
            }
            if !reader.another(b']')? {
                break;
            }
        }
    }
    Ok(())
}

/// The members of an object in a `dict`'s order, each with the records its value holds.
#[derive(Default)]
struct Members {
    order: Vec<(Vec<u8>, Vec<Record>)>,
}

impl Members {
    /// A repeated name keeps its first place and takes its last value.
    fn set(&mut self, name: &[u8], records: Vec<Record>) {
        match self.order.iter_mut().find(|(held, _)| held == name) {
            Some((_, held)) => *held = records,
            None => self.order.push((name.to_vec(), records)),
        }
    }

    fn into_records(self) -> Vec<Record> {
        self.order
            .into_iter()
            .flat_map(|(_, records)| records)
            .collect()
    }
}

/// The records a member of an object holds: the objects of an array, or the member's own
/// object as one record.
fn records_of_member(reader: &mut Reader) -> Res<Vec<Record>> {
    let mut records = Vec::new();
    match reader.value_start()? {
        b'[' => records_of_array(reader, &mut records)?,
        b'{' => records.push(record(reader)?),
        _ => reader.skip_value()?,
    }
    Ok(records)
}

/// The records of the object the reader stands at, member by member; with `wrapper`, a
/// member of that name replaces the object: its array's objects, or its object's members.
fn records_of_object(reader: &mut Reader, wrapper: Option<&[u8]>) -> Res<Vec<Record>> {
    let mut members = Members::default();
    let mut wrapped = None;
    if reader.enter(b'}')? {
        let mut name = Vec::new();
        loop {
            name.clear();
            reader.key(&mut name)?;
            if wrapper == Some(&name[..]) {
                let mut records = Vec::new();
                match reader.value_start()? {
                    b'[' => records_of_array(reader, &mut records)?,
                    b'{' => records = records_of_object(reader, None)?,
                    _ => reader.skip_value()?,
                }
                wrapped = Some(records);
            } else {
                let records = records_of_member(reader)?;
                members.set(&name, records);
            }
            if !reader.another(b'}')? {
                break;
            }
        }
    }
    Ok(wrapped.unwrap_or_else(|| members.into_records()))
}

/// The transcripts of a shard's text, whatever container holds them: a bare array, an
/// object that wraps the array under `transcripts`, an object of arrays keyed by
/// chromosome, or an object of records keyed by id. An item that is no object is passed
/// over, and so is a record without an id.
pub fn read_shard(text: &[u8]) -> Res<Vec<Transcript>> {
    let mut reader = Reader::new(text);
    reader.begin()?;
    let mut records = Vec::new();
    match reader.value_start()? {
        b'[' => records_of_array(&mut reader, &mut records)?,
        b'{' => records = records_of_object(&mut reader, Some(b"transcripts"))?,
        _ => reader.skip_value()?,
    }
    reader.finish()?;
    let mut transcripts = Vec::with_capacity(records.len());
    for record in records {
        match record {
            Record::Nameless => {}
            Record::UnreadableId => {
                return Err(
                    "a record's id is neither a string nor a value Python counts false".to_string(),
                )
            }
            Record::Named(transcript) => transcripts.push(*transcript),
        }
    }
    Ok(transcripts)
}

/// `len(x) if isinstance(x, list) else 0` of `cache.get(name) or record.get(name)`.
fn list_count(in_cache: &Value, top: &Shape) -> usize {
    if truthy(in_cache) {
        match in_cache {
            Value::List(items) => items.len(),
            _ => 0,
        }
    } else {
        top.list_len()
    }
}

/// `len()` of a list or a dict, `None` for any other value.
fn container_len(value: &Value) -> Option<usize> {
    match value {
        Value::List(items) => Some(items.len()),
        Value::Object(object) => Some(object.len()),
        _ => None,
    }
}

/// The mismatches of two records of one id, a bit a name of [`mismatch_name`]. Fails where
/// either record's feature cache is no object.
pub fn compare(baseline: &Transcript, builder: &Transcript) -> Res<u32> {
    let mut mismatches = 0u32;
    for field in 0..COMPARE_FIELDS.len() {
        if py_ne(&baseline.fields[field], &builder.fields[field]) {
            mismatches |= 1 << field;
        }
    }
    let (Some(b), Some(c)) = (&baseline.cache, &builder.cache) else {
        return Err("the feature cache of a record is no object".to_string());
    };
    let exons = |cache: &FeatureCache, record: &Transcript| {
        list_count(&cache.fields[VEFC_SORTED_EXONS], &record.top[0])
    };
    if exons(b, baseline) != exons(c, builder) {
        mismatches |= 1 << EXON_COUNT;
    }
    let introns = |cache: &FeatureCache, record: &Transcript| {
        list_count(&cache.fields[VEFC_INTRONS], &record.top[1])
    };
    if introns(b, baseline) != introns(c, builder) {
        mismatches |= 1 << INTRON_COUNT;
    }
    let mapper = |cache: &FeatureCache, record: &Transcript| {
        if cache.mapper.truthy {
            cache.mapper.container_len()
        } else {
            record.top[2].container_len()
        }
    };
    if mapper(b, baseline) != mapper(c, builder) {
        mismatches |= 1 << MAPPER_PAIR_COUNT;
    }
    for field in 0..VEFC_FIELDS.len() {
        let (bval, cval) = (&b.fields[field], &c.fields[field]);
        if py_ne(bval, cval) {
            // Two lists or dicts of different lengths are counted under the `_count` name.
            let lengths_differ = match (container_len(bval), container_len(cval)) {
                (Some(x), Some(y)) => x != y,
                _ => false,
            };
            mismatches |= 1 << (VEFC_FIRST + 2 * field + usize::from(!lengths_differ));
        }
    }
    Ok(mismatches)
}

/// The comparison of one shard name: the counts of the ids each side holds under it and of
/// the mismatches of the ids both hold.
pub struct Region {
    pub baseline_count: usize,
    pub builder_count: usize,
    pub shared: usize,
    pub only_baseline: usize,
    pub only_builder: usize,
    /// The first ten ids only one side holds, in sorted order.
    pub first_only_baseline: Vec<Box<[u8]>>,
    pub first_only_builder: Vec<Box<[u8]>>,
    pub matched: usize,
    pub differed: usize,
    /// How many shared ids differ under each name of [`mismatch_name`].
    pub mismatches: [usize; MISMATCH_NAMES],
    /// The place, among the shared ids in sorted order, of the first id that differs under
    /// each name.
    pub first_mismatch: [Option<usize>; MISMATCH_NAMES],
}

impl Region {
    pub fn has_differences(&self) -> bool {
        self.only_baseline > 0 || self.only_builder > 0 || self.differed > 0
    }
}

fn named(path: &Path, detail: impl std::fmt::Display) -> String {
    format!("{}: {detail}", path.display())
}

fn load(path: Option<&Path>) -> Res<Vec<Transcript>> {
    let Some(path) = path else {
        return Ok(Vec::new());
    };
    match input::read_shard(path).map_err(|e| named(path, e))? {
        None => Ok(Vec::new()),
        Some(text) => read_shard(&text).map_err(|e| named(path, e)),
    }
}

/// The transcripts by id; an id a shard holds twice is the later record.
fn by_id(transcripts: &[Transcript]) -> HashMap<&[u8], &Transcript> {
    transcripts
        .iter()
        .map(|transcript| (&*transcript.id, transcript))
        .collect()
}

/// How many ids `only` keeps of those one side alone holds.
const ONLY_KEPT: usize = 10;

/// Compares the shards of the two sides under one name.
pub fn compare_region(pair: &Pair) -> Res<Region> {
    let baseline = load(pair.baseline.as_deref())?;
    let builder = load(pair.builder.as_deref())?;
    let (b, c) = (by_id(&baseline), by_id(&builder));
    let mut ids: Vec<&[u8]> = b.keys().copied().collect();
    ids.sort_unstable();
    let mut region = Region {
        baseline_count: b.len(),
        builder_count: c.len(),
        shared: 0,
        only_baseline: 0,
        only_builder: 0,
        first_only_baseline: Vec::new(),
        first_only_builder: Vec::new(),
        matched: 0,
        differed: 0,
        mismatches: [0; MISMATCH_NAMES],
        first_mismatch: [None; MISMATCH_NAMES],
    };
    for id in ids {
        let Some(other) = c.get(id) else {
            region.only_baseline += 1;
            if region.first_only_baseline.len() < ONLY_KEPT {
                region.first_only_baseline.push(id.into());
            }
            continue;
        };
        let mismatches = compare(b[id], other).map_err(|e| {
            let path = pair
                .baseline
                .as_deref()
                .expect("a shared id has a baseline shard");
            named(
                path,
                format!("transcript {}: {e}", String::from_utf8_lossy(id)),
            )
        })?;
        if mismatches == 0 {
            region.matched += 1;
        } else {
            region.differed += 1;
            for name in 0..MISMATCH_NAMES {
                if mismatches & (1 << name) != 0 {
                    region.mismatches[name] += 1;
                    region.first_mismatch[name].get_or_insert(region.shared);
                }
            }
        }
        region.shared += 1;
    }
    region.only_builder = c.len() - region.shared;
    if region.only_builder > 0 {
        let mut only: Vec<&[u8]> = c.keys().copied().filter(|id| !b.contains_key(id)).collect();
        only.sort_unstable();
        only.truncate(ONLY_KEPT);
        region.first_only_builder = only.into_iter().map(Box::from).collect();
    }
    Ok(region)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(text: &str) -> Vec<String> {
        read_shard(text.as_bytes())
            .unwrap_or_else(|e| panic!("{text}: {e}"))
            .iter()
            .map(|t| String::from_utf8_lossy(&t.id).into_owned())
            .collect()
    }

    #[test]
    fn records_are_found_under_every_container() {
        let cases = [
            (
                r#"[{"stable_id": "A"}, null, 5, "x", [], {"stable_id": "B"}]"#,
                vec!["A", "B"],
            ),
            (
                r#"{"transcripts": [{"stable_id": "A"}, null, {"id": "B"}]}"#,
                vec!["A", "B"],
            ),
            (
                r#"{"21": [{"stable_id": "A"}, 7], "22": [{"stable_id": "B"}], "x": null}"#,
                vec!["A", "B"],
            ),
            (
                r#"{"A": {"stable_id": "A"}, "B": {"id": "B"}, "n": 1, "s": "t"}"#,
                vec!["A", "B"],
            ),
            (
                r#"{"transcripts": {"21": [{"stable_id": "A"}], "B": {"stable_id": "B"}, "n": 3}}"#,
                vec!["A", "B"],
            ),
            // The wrapper replaces the object, wherever it stands and whatever it holds.
            (
                r#"{"other": [{"stable_id": "X"}], "transcripts": [{"stable_id": "A"}]}"#,
                vec!["A"],
            ),
            (
                r#"{"other": [{"stable_id": "X"}], "transcripts": null}"#,
                vec![],
            ),
            (
                r#"{"other": [{"stable_id": "X"}], "transcripts": 5}"#,
                vec![],
            ),
            (
                r#"{"transcripts": [{"stable_id": "A"}], "transcripts": [{"stable_id": "B"}]}"#,
                vec!["B"],
            ),
            // The wrapper is read at the top alone.
            (
                r#"{"transcripts": {"transcripts": [{"stable_id": "A"}], "k": {"id": "B"}}}"#,
                vec!["A", "B"],
            ),
            // A repeated member keeps its first place and its last value.
            (
                r#"{"1": [{"id": "A"}], "2": [{"id": "B"}], "1": [{"id": "C"}]}"#,
                vec!["C", "B"],
            ),
            ("[]", vec![]),
            ("{}", vec![]),
            ("null", vec![]),
            ("5", vec![]),
            (r#""text""#, vec![]),
        ];
        for (text, expected) in cases {
            assert_eq!(ids(text), expected, "{text}");
        }
    }

    #[test]
    fn a_record_is_named_by_its_stable_id_or_its_id() {
        let cases = [
            (r#"[{"stable_id": "S", "id": "I"}]"#, vec!["S"]),
            (r#"[{"id": "I", "stable_id": "S"}]"#, vec!["S"]),
            (r#"[{"id": "I"}]"#, vec!["I"]),
            (r#"[{"stable_id": "", "id": "I"}]"#, vec!["I"]),
            (r#"[{"stable_id": null, "id": "I"}]"#, vec!["I"]),
            (r#"[{"stable_id": 0, "id": "I"}]"#, vec!["I"]),
            (r#"[{"stable_id": [], "id": "I"}]"#, vec!["I"]),
            (r#"[{"stable_id": "A", "stable_id": "B"}]"#, vec!["B"]),
            (
                r#"[{"stable_id": ""}, {"id": ""}, {"id": null}, {"id": false}, {}]"#,
                vec![],
            ),
            (r#"[{"stable_id": "", "id": 0.0}]"#, vec![]),
        ];
        for (text, expected) in cases {
            assert_eq!(ids(text), expected, "{text}");
        }
        for text in [
            r#"[{"stable_id": 5}]"#,
            r#"[{"stable_id": true}]"#,
            r#"[{"stable_id": ["A"]}]"#,
            r#"[{"stable_id": {"a": 1}}]"#,
            r#"[{"id": 1.5}]"#,
            r#"[{"stable_id": "", "id": [0]}]"#,
        ] {
            assert!(read_shard(text.as_bytes()).is_err(), "{text}");
        }
        // A record the wrapper sets aside is never looked at.
        assert_eq!(
            ids(r#"{"other": [{"stable_id": 5}], "transcripts": [{"stable_id": "A"}]}"#),
            ["A"]
        );
    }

    fn one(text: &str) -> Transcript {
        let mut all =
            read_shard(format!("[{text}]").as_bytes()).unwrap_or_else(|e| panic!("{text}: {e}"));
        assert_eq!(all.len(), 1, "{text}");
        all.pop().unwrap()
    }

    fn mismatches(baseline: &str, builder: &str) -> Vec<String> {
        let found = compare(&one(baseline), &one(builder)).unwrap();
        (0..MISMATCH_NAMES)
            .filter(|name| found & (1 << name) != 0)
            .map(mismatch_name)
            .collect()
    }

    #[test]
    fn the_names_stand_in_the_order_a_comparison_meets_them() {
        let names: Vec<String> = (0..MISMATCH_NAMES).map(mismatch_name).collect();
        assert_eq!(
            names,
            [
                "start",
                "end",
                "strand",
                "biotype",
                "source",
                "gene_stable_id",
                "gene_symbol",
                "hgnc_id",
                "canonical",
                "exon_count",
                "intron_count",
                "mapper_pair_count",
                "vefc_introns_count",
                "vefc_introns",
                "vefc_sorted_exons_count",
                "vefc_sorted_exons",
                "vefc_five_prime_utr_count",
                "vefc_five_prime_utr",
                "vefc_three_prime_utr_count",
                "vefc_three_prime_utr",
                "vefc_codon_table_count",
                "vefc_codon_table",
                "vefc_protein_features_count",
                "vefc_protein_features",
                "vefc_seq_edits_count",
                "vefc_seq_edits",
            ]
        );
    }

    #[test]
    fn members_are_compared_as_python_compares_them() {
        let same = r#"{"id": "T", "start": 1, "end": 2.0, "strand": true, "canonical": null}"#;
        assert!(mismatches(same, r#"{"id": "T", "start": 1.0, "end": 2, "strand": 1}"#).is_empty());
        assert_eq!(
            mismatches(
                r#"{"id": "T", "start": 1, "biotype": "a", "hgnc_id": null, "canonical": 0}"#,
                r#"{"id": "T", "start": 2, "biotype": "b", "hgnc_id": "", "canonical": null}"#,
            ),
            ["start", "biotype", "hgnc_id", "canonical"]
        );
        // NaN differs from NaN as a member, and equals it inside a list.
        assert_eq!(
            mismatches(
                r#"{"id": "T", "source": NaN}"#,
                r#"{"id": "T", "source": NaN}"#
            ),
            ["source"]
        );
        assert!(mismatches(
            r#"{"id": "T", "source": [NaN]}"#,
            r#"{"id": "T", "source": [NaN]}"#
        )
        .is_empty());
    }

    #[test]
    fn counts_come_from_the_feature_cache_or_the_record() {
        let in_cache = r#"{"id": "T", "vefc": {"sorted_exons": [1, 2], "introns": [1], "mapper": {"a": 1, "b": 2}}}"#;
        let on_record = r#"{"id": "T", "sorted_exons": [3, 4], "introns": [5], "mapper": [6, 7]}"#;
        // The same counts on both sides; the cache members themselves differ.
        assert_eq!(
            mismatches(in_cache, on_record),
            ["vefc_introns", "vefc_sorted_exons"]
        );
        // An empty member of the cache gives way to the record's.
        let emptied = r#"{"id": "T", "sorted_exons": [3, 4], "introns": [5], "mapper": [6, 7],
            "vefc": {"sorted_exons": [], "introns": [], "mapper": {}}}"#;
        assert_eq!(
            mismatches(in_cache, emptied),
            ["vefc_introns_count", "vefc_sorted_exons_count"]
        );
        // A member that is no list counts zero; a mapper counts as a list or a dict.
        assert_eq!(
            mismatches(
                in_cache,
                r#"{"id": "T", "vefc": {"sorted_exons": {"a": 1, "b": 2}, "introns": "x", "mapper": "ab"}}"#,
            ),
            [
                "exon_count",
                "intron_count",
                "mapper_pair_count",
                "vefc_introns",
                "vefc_sorted_exons"
            ]
        );
        // A dict that names a key twice counts it once.
        assert!(mismatches(
            in_cache,
            r#"{"id": "T", "vefc": {"sorted_exons": [1, 2], "introns": [1], "mapper": {"a": 1, "b": 2, "a": 3}}}"#,
        )
        .is_empty());
    }

    #[test]
    fn the_feature_cache_is_the_long_member_where_python_counts_it_true() {
        let long = r#"{"id": "T", "variation_effect_feature_cache": {"codon_table": 1}, "vefc": {"codon_table": 2}}"#;
        let short =
            r#"{"id": "T", "variation_effect_feature_cache": {}, "vefc": {"codon_table": 2}}"#;
        let neither = r#"{"id": "T"}"#;
        assert_eq!(mismatches(long, short), ["vefc_codon_table"]);
        assert!(mismatches(short, r#"{"id": "T", "vefc": {"codon_table": 2.0}}"#).is_empty());
        assert!(mismatches(neither, r#"{"id": "T", "vefc": {}}"#).is_empty());
        assert!(mismatches(
            neither,
            r#"{"id": "T", "variation_effect_feature_cache": null, "vefc": {"other": 1}}"#
        )
        .is_empty());
        // A member that is no object ends the comparison where it is the one read.
        for text in [
            r#"{"id": "T", "vefc": null}"#,
            r#"{"id": "T", "vefc": []}"#,
            r#"{"id": "T", "vefc": "x"}"#,
            r#"{"id": "T", "variation_effect_feature_cache": [1]}"#,
            r#"{"id": "T", "variation_effect_feature_cache": 5, "vefc": {}}"#,
        ] {
            assert!(compare(&one(text), &one(neither)).is_err(), "{text}");
            assert!(compare(&one(neither), &one(text)).is_err(), "{text}");
        }
        assert!(mismatches(
            neither,
            r#"{"id": "T", "variation_effect_feature_cache": {"x": 1}, "vefc": null}"#
        )
        .is_empty());
    }

    #[test]
    fn cache_members_take_the_count_name_where_two_containers_differ_in_length() {
        let cache = |members: &str| format!(r#"{{"id": "T", "vefc": {{{members}}}}}"#);
        assert_eq!(
            mismatches(
                &cache(r#""seq_edits": [1]"#),
                &cache(r#""seq_edits": [1, 2]"#)
            ),
            ["vefc_seq_edits_count"]
        );
        assert_eq!(
            mismatches(&cache(r#""seq_edits": [1]"#), &cache(r#""seq_edits": [2]"#)),
            ["vefc_seq_edits"]
        );
        assert_eq!(
            mismatches(
                &cache(r#""seq_edits": [1]"#),
                &cache(r#""seq_edits": {"a": 1, "b": 2}"#)
            ),
            ["vefc_seq_edits_count"]
        );
        assert_eq!(
            mismatches(
                &cache(r#""seq_edits": [1]"#),
                &cache(r#""seq_edits": {"a": 1}"#)
            ),
            ["vefc_seq_edits"]
        );
        assert_eq!(
            mismatches(
                &cache(r#""codon_table": 1"#),
                &cache(r#""codon_table": [1]"#)
            ),
            ["vefc_codon_table"]
        );
        assert_eq!(
            mismatches(&cache(r#""three_prime_utr": "x""#), &cache("")),
            ["vefc_three_prime_utr"]
        );
        assert!(mismatches(&cache(r#""protein_features": null"#), &cache("")).is_empty());
    }
}
