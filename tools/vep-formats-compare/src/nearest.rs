// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The nearest transcription start sites of a cache.

use crate::argparse::path_text;
use crate::cache;
use crate::error::{Error, Result};
use crate::fields::{extra_accepted, python_type};
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::shared::ABSENT;
use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::Path;

/// The values of a `NEAREST` field as every format joins them (`,` or `&`).
pub fn nearest_list(value: &str) -> Vec<String> {
    if value == "-" || value == ABSENT {
        return Vec::new();
    }
    value
        .split([',', '&'])
        .filter(|item| !item.is_empty())
        .map(str::to_string)
        .collect()
}

/// The protein-coding transcription start sites of a candidate JSON cache, by transcript
/// stable id, gene stable id and gene symbol: a transcript's site is its start on the forward
/// strand and its end on the reverse (`AnnotationType/Transcript.pm:749`), a gene's or
/// symbol's the sites of its transcripts.
#[derive(Debug)]
pub struct NearestSites {
    pub sites: IndexMap<String, BTreeSet<i64>>,
    /// The ids and symbols the cache writes as numbers: keys that no `NEAREST` value names,
    /// counted in `len` and held nowhere.
    pub numbered: usize,
    /// The sites of each sequence, under the name of the directory its shards stand in.
    pub by_sequence: IndexMap<String, BTreeSet<i64>>,
}

/// What a transcript files a site under.
enum Name {
    Text(String),
    Number(String),
}

impl NearestSites {
    pub fn new(sites: &IndexMap<String, BTreeSet<i64>>) -> Self {
        NearestSites {
            sites: sites.clone(),
            numbered: 0,
            by_sequence: IndexMap::new(),
        }
    }

    /// The sites `named` file with every one of them on the sequence `sequence`.
    pub fn on_one_sequence(sequence: &str, named: &IndexMap<String, BTreeSet<i64>>) -> Self {
        let mut sites = NearestSites::new(named);
        sites.by_sequence.insert(
            sequence.to_string(),
            named.values().flatten().copied().collect(),
        );
        sites
    }

    /// The count of the names that file a site, the numbered ones among them.
    pub fn len(&self) -> usize {
        self.sites.len() + self.numbered
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn from_cache(cache: &Path) -> Result<Self> {
        let transcripts = cache::joined(cache, "transcripts");
        if !transcripts.is_dir() {
            return Err(Error::comparator(format!(
                "--nearest-sites: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        let shards = cache::shards_by_suffix(&transcripts);
        let read = cache::map_shards(&shards, cache::threads(), |shard, transcripts| {
            let sequence = shard
                .parent()
                .and_then(Path::file_name)
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_default();
            let mut found = Vec::new();
            let mut placed = Vec::new();
            for transcript in &transcripts {
                if transcript.get("biotype").and_then(Value::as_str) != Some("protein_coding") {
                    continue;
                }
                let forward = transcript
                    .get("strand")
                    .is_none_or(|strand| pyjson::py_str(strand) == "1");
                let site =
                    Self::coordinate(shard, transcript, if forward { "start" } else { "end" })?;
                placed.push(site);
                for key in ["stable_id", "gene_stable_id", "gene_symbol"] {
                    if let Some(value) = transcript.get(key).filter(|value| pyjson::truthy(value)) {
                        found.push((name(value), site));
                    }
                }
            }
            Ok((sequence, placed, found))
        })?;
        let mut sites: IndexMap<String, BTreeSet<i64>> = IndexMap::new();
        let mut numbers: BTreeSet<String> = BTreeSet::new();
        let mut by_sequence: IndexMap<String, BTreeSet<i64>> = IndexMap::new();
        for (sequence, placed, found) in read {
            if !placed.is_empty() {
                by_sequence.entry(sequence).or_default().extend(placed);
            }
            for (name, site) in found {
                match name {
                    Name::Text(text) => {
                        sites.entry(text).or_default().insert(site);
                    }
                    Name::Number(text) => {
                        numbers.insert(text);
                    }
                }
            }
        }
        Ok(NearestSites {
            sites,
            numbered: numbers.len(),
            by_sequence,
        })
    }

    /// The transcript objects of a parsed shard: a bare array, or a range-keyed map of
    /// arrays; a `null` in either is a placeholder and is dropped.
    pub fn shard_transcripts(shard: &Path, parsed: &Value) -> Result<Vec<Value>> {
        cache::transcripts_of(shard, parsed.clone())
    }

    /// The `key` coordinate of a transcript, an integer or a string of digits.
    pub fn coordinate(shard: &Path, transcript: &Value, key: &str) -> Result<i64> {
        let coordinate = match transcript.get(key) {
            Some(Value::Int(number)) if *number >= 0 => Some(*number),
            Some(Value::Str(digits)) if pyfmt::is_decimal(digits) => {
                match pyfmt::int_checked(digits) {
                    Ok(number) => number.parse::<i64>().ok(),
                    Err(error) => crate::uncaught!("{error}"),
                }
            }
            _ => None,
        };
        coordinate.ok_or_else(|| {
            Error::comparator(format!(
                "{}: transcript {} has no integer {key}",
                path_text(shard),
                transcript
                    .get("stable_id")
                    .map_or_else(|| "None".to_string(), pyjson::py_repr)
            ))
        })
    }

    /// Whether two `NEAREST` values name the same winning sites: every expected value is
    /// among the actual values and every actual value starts at a site of an expected value;
    /// a value the cache lacks starts nowhere, so it never agrees.
    pub fn agree(&self, expected: &str, actual: &str) -> bool {
        let expected_values = nearest_list(expected);
        let actual_values = nearest_list(actual);
        if expected_values.is_empty() || actual_values.is_empty() {
            return expected == actual;
        }
        if expected_values
            .iter()
            .any(|value| !actual_values.contains(value))
        {
            return false;
        }
        let expected_sites: BTreeSet<i64> = expected_values
            .iter()
            .filter_map(|value| self.sites.get(value))
            .flatten()
            .copied()
            .collect();
        actual_values.iter().all(|value| {
            self.sites
                .get(value)
                .is_some_and(|sites| !sites.is_disjoint(&expected_sites))
        })
    }
}

impl NearestSites {
    /// The sites of the sequence a Location names, under that name or under it with a
    /// leading `chr` taken off or put on.
    fn sequence_sites(&self, sequence: &str) -> Option<&BTreeSet<i64>> {
        self.by_sequence
            .get(sequence)
            .or_else(|| match sequence.strip_prefix("chr") {
                Some(bare) => self.by_sequence.get(bare),
                None => self.by_sequence.get(&format!("chr{sequence}")),
            })
    }

    /// Whether a `NEAREST` pair has the shape of `nearest_start_site_one_base_above` on a
    /// row at `location`: `site_above_start_span` over the Location's sequence and span.
    pub fn site_above_start_pair(&self, location: &str, expected: &str, actual: &str) -> bool {
        location_span(location).is_some_and(|(sequence, low, high)| {
            self.site_above_start_span(sequence, low, high, expected, actual)
        })
    }

    /// Whether a `NEAREST` pair has the shape of `nearest_start_site_one_base_above` on a
    /// record of `sequence` from `low` to `high`: a site stands on the span's first base or
    /// one base above it, the two searches of `TranscriptTree.pm:233-236` keep one set of
    /// sites, those two searches and a third that begins on the first base keep another,
    /// the reference's values name the first set and the candidate's the second. A set is
    /// named when every value starts at one of its sites and every one of its sites starts
    /// a value; the empty set is named by no value.
    pub fn site_above_start_span(
        &self,
        sequence: &str,
        low: i64,
        high: i64,
        expected: &str,
        actual: &str,
    ) -> bool {
        let Some(sites) = self.sequence_sites(sequence) else {
            return false;
        };
        let by_the_start =
            sites.contains(&low) || low.checked_add(1).is_some_and(|site| sites.contains(&site));
        if !by_the_start {
            return false;
        }
        let passed_over = reached(sites, low, high, &[2]);
        let counted = reached(sites, low, high, &[0, 2]);
        passed_over != counted
            && self.names(&nearest_list(expected), &passed_over)
            && self.names(&nearest_list(actual), &counted)
    }

    /// Whether `values` name the sites `reached` and no other.
    fn names(&self, values: &[String], reached: &BTreeSet<i64>) -> bool {
        let starts_at = |value: &String, site: &i64| {
            self.sites
                .get(value)
                .is_some_and(|sites| sites.contains(site))
        };
        values
            .iter()
            .all(|value| reached.iter().any(|site| starts_at(value, site)))
            && reached
                .iter()
                .all(|site| values.iter().any(|value| starts_at(value, site)))
    }
}

/// A position written in ASCII digits alone.
fn position(text: &str) -> Option<i64> {
    (!text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit()))
        .then(|| text.parse::<i64>().ok())
        .flatten()
}

/// The sequence and the first and last position of a Location, `21:5000` or `21:4990-5010`;
/// `None` for a Location with no position.
fn location_span(location: &str) -> Option<(&str, i64, i64)> {
    let (sequence, span) = location.rsplit_once(':')?;
    let (first, last) = span.split_once('-').unwrap_or((span, span));
    let (first, last) = (position(first)?, position(last)?);
    Some((sequence, first.min(last), first.max(last)))
}

/// The sequence and span of a JSON record's entry, from its `seq_region_name`, `start` and
/// `end`, an insertion's start standing one above its end; `None` for an entry without them.
pub fn record_span(
    sequence: Option<&str>,
    start: Option<&str>,
    end: Option<&str>,
) -> Option<(String, i64, i64)> {
    let (first, last) = (position(start?)?, position(end?)?);
    Some((sequence?.to_string(), first.min(last), first.max(last)))
}

/// The sites the searches of `TranscriptTree.pm:233-236` keep for the span `low` to `high`:
/// the highest site at or below `high` and, for each offset of `above`, the lowest site at
/// or above `low` plus that offset, each measured to its nearer end of the span, the
/// nearest kept and all on a tie (`:242-253`). `above` is `[2]` for the reference, whose
/// tree passes over the interval holding the start and the one ending at it, and `[0, 2]`
/// for the candidate, which takes the passed-over site as a further candidate.
fn reached(sites: &BTreeSet<i64>, low: i64, high: i64, above: &[i64]) -> BTreeSet<i64> {
    let down = sites.range(..=high).next_back().copied();
    let ups = above.iter().filter_map(|offset| {
        low.checked_add(*offset)
            .and_then(|from| sites.range(from..).next().copied())
    });
    let found: Vec<i64> = down.into_iter().chain(ups).collect();
    let distance = |site: i64| low.abs_diff(site).min(high.abs_diff(site));
    let nearest = found.iter().copied().map(distance).min();
    found
        .into_iter()
        .filter(|site| Some(distance(*site)) == nearest)
        .collect()
}

/// The key a truthy `stable_id`, `gene_stable_id` or `gene_symbol` files a site under: a
/// whole float, `true` and the integer they equal are one key, and every `NaN` of a JSON
/// text is the one object, so one key.
fn name(value: &Value) -> Name {
    match value {
        Value::Str(text) => Name::Text(text.clone()),
        Value::Int(number) => Name::Number(number.to_string()),
        Value::BigInt(digits) => Name::Number(digits.clone()),
        Value::Bool(_) => Name::Number("1".to_string()),
        Value::Float(number) if number.is_finite() && number.fract() == 0.0 => {
            match pyfmt::int_of_f64(*number) {
                Ok(whole) => Name::Number(whole),
                Err(_) => Name::Number(pyfmt::repr_f64(*number)),
            }
        }
        Value::Float(number) => Name::Number(pyfmt::repr_f64(*number)),
        Value::Null | Value::List(_) | Value::Dict(_) => crate::uncaught!(
            "cannot use '{0}' as a dict key (unhashable type: '{0}')",
            python_type(value)
        ),
    }
}

/// `extra_accepted` for an agreeing `NEAREST` value.
pub fn extra_nearest_accepted(expected: &str, actual: &str, sites: &NearestSites) -> bool {
    extra_accepted(expected, actual, &mut |key, value, other| {
        key == "NEAREST" && sites.agree(value, other)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sites(named: &[(&str, &[i64])]) -> NearestSites {
        let sites: IndexMap<String, BTreeSet<i64>> = named
            .iter()
            .map(|(name, sites)| (name.to_string(), sites.iter().copied().collect()))
            .collect();
        NearestSites::new(&sites)
    }

    #[test]
    fn nearest_values_agree_on_the_winning_sites_not_on_the_ids() {
        let sites = sites(&[
            ("A", &[1000]),
            ("B", &[1000]),
            ("R", &[5000]),
            ("F", &[9000]),
            ("G1", &[1000]),
            ("G2", &[5000]),
            ("G3", &[9000]),
            ("SYM1", &[1000]),
            ("SYM3", &[9000]),
        ]);
        assert!(sites.agree("A", "A,B"));
        assert!(sites.agree("B", "A,B"));
        assert!(sites.agree("A", "A&B"));
        assert!(sites.agree("R,F", "R,F"));
        assert!(sites.agree("G1", "G1"));
        assert!(!sites.agree("A", "B"));
        assert!(!sites.agree("A", "A,R"));
        assert!(!sites.agree("R", "R,F"));
        assert!(!sites.agree("A", "A,X"));
        assert!(!sites.agree("X", "X"));
        assert!(sites.agree("-", "-") && sites.agree(ABSENT, ABSENT));
        assert!(!sites.agree("A", "-") && !sites.agree("-", "A"));
        assert!(!sites.agree("-", ABSENT));
        assert_eq!(nearest_list("A&B"), ["A", "B"]);
        assert!(nearest_list("-").is_empty());
        assert!(extra_nearest_accepted(
            "IMPACT=HIGH;NEAREST=A",
            "IMPACT=HIGH;NEAREST=A,B",
            &sites
        ));
        assert!(!extra_nearest_accepted(
            "IMPACT=HIGH;NEAREST=A",
            "NEAREST=A,B;IMPACT=HIGH",
            &sites
        ));
        assert!(!extra_nearest_accepted(
            "IMPACT=HIGH;NEAREST=A",
            "IMPACT=LOW;NEAREST=A,B",
            &sites
        ));
        assert!(!extra_nearest_accepted(
            "IMPACT=HIGH",
            "IMPACT=HIGH",
            &sites
        ));
    }

    /// The table of the class's tests: one transcript a site but for the two at 1000, a gene
    /// and a symbol a transcript.
    fn one_sequence() -> NearestSites {
        let named: IndexMap<String, BTreeSet<i64>> = [
            ("A", 1000),
            ("B", 1000),
            ("GENE_AB", 1000),
            ("R", 5000),
            ("GENE_R", 5000),
            ("F", 9000),
            ("GENE_F", 9000),
        ]
        .iter()
        .map(|(name, site)| (name.to_string(), BTreeSet::from([*site])))
        .collect();
        NearestSites::on_one_sequence("21", &named)
    }

    #[test]
    fn a_site_one_base_above_the_start_is_the_pair_the_class_takes() {
        let sites = one_sequence();
        // 4999: the reference's searches reach 1000 (3999 below) and 9000 (4001 above) and
        // keep 1000; the candidate's reach 5000, one base above.
        assert!(sites.site_above_start_pair("21:4999", "A", "R"));
        assert!(sites.site_above_start_pair("21:4999", "B", "R"));
        assert!(sites.site_above_start_pair("21:4999", "GENE_AB", "GENE_R"));
        assert!(sites.site_above_start_pair("chr21:4999", "A&B", "R"));
        // Another sequence holds no site, and another position no site one base above.
        assert!(!sites.site_above_start_pair("22:4999", "A", "R"));
        assert!(!sites.site_above_start_pair("21:4998", "A", "R"));
        assert!(!sites.site_above_start_pair("21:5000", "A", "R"));
        // Each value names its own search's sites and no other.
        assert!(!sites.site_above_start_pair("21:4999", "F", "R"));
        assert!(!sites.site_above_start_pair("21:4999", "A", "F"));
        assert!(!sites.site_above_start_pair("21:4999", "A", "A"));
        assert!(!sites.site_above_start_pair("21:4999", "R", "R"));
        assert!(!sites.site_above_start_pair("21:4999", "A", "R,F"));
        assert!(!sites.site_above_start_pair("21:4999", "A,X", "R"));
        assert!(!sites.site_above_start_pair("21:4999", "-", "R"));
        assert!(!sites.site_above_start_pair("21:4999", "A", "-"));
        // A Location with no position names no span.
        for location in ["21", "21:", "21:x", "21:4999-", "21:-4999", "4999"] {
            assert!(
                !sites.site_above_start_pair(location, "A", "R"),
                "{location}"
            );
        }
        // A table read with no sequence takes nothing.
        let unplaced = NearestSites::new(&sites.sites);
        assert!(!unplaced.site_above_start_pair("21:4999", "A", "R"));
    }

    #[test]
    fn the_two_searches_differ_only_where_the_site_above_the_start_decides() {
        let named = |sites: &[(&str, i64)]| -> NearestSites {
            let named: IndexMap<String, BTreeSet<i64>> = sites
                .iter()
                .map(|(name, site)| (name.to_string(), BTreeSet::from([*site])))
                .collect();
            NearestSites::on_one_sequence("21", &named)
        };
        // A site one base below ties with the site one base above: the reference names the
        // lower site alone, the candidate both.
        let beside = named(&[("BELOW", 4998), ("ABOVE", 5000), ("FAR", 9000)]);
        assert!(beside.site_above_start_pair("21:4999", "BELOW", "BELOW,ABOVE"));
        assert!(!beside.site_above_start_pair("21:4999", "BELOW", "ABOVE"));
        // A site on the record is at distance zero in both searches: no pair is taken.
        let on = named(&[("ON", 4999), ("ABOVE", 5000)]);
        assert!(!on.site_above_start_pair("21:4999", "ON", "ABOVE"));
        assert!(!on.site_above_start_pair("21:4999", "ON", "ON,ABOVE"));
        // A lone site: the reference's searches reach nothing and it writes no value.
        let lone = named(&[("ABOVE", 5000)]);
        assert!(lone.site_above_start_pair("21:4999", "-", "ABOVE"));
        assert!(lone.site_above_start_pair("21:4999", "", "ABOVE"));
        assert!(!lone.site_above_start_pair("21:4999", "ABOVE", "ABOVE"));
        // A span: the site at start + 1 is one base from the start, a higher site inside
        // the span six; the reference's searches both reach the higher one.
        let inside = named(&[("BY_START", 4991), ("HIGHER", 4996)]);
        assert!(inside.site_above_start_pair("21:4990-5010", "HIGHER", "BY_START"));
        // A span that ends on the site reaches it from below in both searches, and so does
        // an insertion written between the record and the site.
        let table = one_sequence();
        assert!(!table.site_above_start_pair("21:4999-5000", "R", "R"));
        assert!(!table.site_above_start_pair("21:4999-5000", "A", "R"));
        assert_eq!(location_span("21:4990-5010"), Some(("21", 4990, 5010)));
        assert_eq!(location_span("21:5010-4990"), Some(("21", 4990, 5010)));
        assert_eq!(location_span("a:b:7"), Some(("a:b", 7, 7)));
        assert_eq!(
            record_span(Some("21"), Some("5000"), Some("4999")),
            Some(("21".to_string(), 4999, 5000))
        );
        assert_eq!(record_span(None, Some("5000"), Some("4999")), None);
        assert_eq!(record_span(Some("21"), Some("5000"), None), None);
        assert_eq!(record_span(Some("21"), Some("x"), Some("4999")), None);
    }

    #[test]
    fn a_site_on_the_first_base_of_a_longer_span_is_the_pair_the_class_takes() {
        let named = |sites: &[(&str, i64)]| -> NearestSites {
            let named: IndexMap<String, BTreeSet<i64>> = sites
                .iter()
                .map(|(name, site)| (name.to_string(), BTreeSet::from([*site])))
                .collect();
            NearestSites::on_one_sequence("21", &named)
        };
        // Sites on both positions of a two-position span, an insertion's too: the
        // reference's searches reach the upper site and the site at 9000 and keep the upper;
        // the candidate's reach both at distance zero.
        let adjacent = named(&[("LOW", 5000), ("HIGH", 5001), ("FAR", 9000)]);
        assert!(adjacent.site_above_start_pair("21:5000-5001", "HIGH", "HIGH,LOW"));
        assert!(adjacent.site_above_start_pair("21:5000-5001", "HIGH", "HIGH&LOW"));
        assert!(!adjacent.site_above_start_pair("21:5000-5001", "HIGH", "LOW"));
        assert!(!adjacent.site_above_start_pair("21:5000-5001", "LOW", "HIGH,LOW"));
        assert!(!adjacent.site_above_start_pair("21:5000-5001", "HIGH", "HIGH"));
        assert!(!adjacent.site_above_start_pair("21:5000-5001", "HIGH,LOW", "HIGH,LOW"));
        // A site on the first base and a site five bases below the end: the reference's
        // searches both reach the higher one, the candidate's keep the one at distance zero.
        let inside = named(&[("FIRST", 5000), ("INSIDE", 5030), ("FAR", 9000)]);
        assert!(inside.site_above_start_pair("21:5000-5035", "INSIDE", "FIRST"));
        assert!(!inside.site_above_start_pair("21:5000-5035", "FIRST", "INSIDE"));
        assert!(!inside.site_above_start_pair("21:5000-5035", "INSIDE", "INSIDE,FIRST"));
        // With the higher site on the last base the two ends tie in the candidate's search.
        assert!(inside.site_above_start_pair("21:5000-5030", "INSIDE", "INSIDE,FIRST"));
        // The only site a span covers, on its first base, is reached from below the end in
        // both searches, and so is the site a one-base record stands on.
        let lone = named(&[("FIRST", 5000), ("FAR", 9000)]);
        assert!(!lone.site_above_start_pair("21:5000-5003", "FIRST", "FIRST"));
        assert!(!lone.site_above_start_pair("21:5000-5003", "FAR", "FIRST"));
        assert!(!lone.site_above_start_pair("21:5000", "FAR", "FIRST"));
        // A span that starts one base above the lower site covers the higher site alone, and
        // both searches return the covered site.
        assert!(!inside.site_above_start_pair("21:5001-5033", "INSIDE", "FIRST"));
        // Sites on the first base and one base above it, and a higher site inside: the
        // reference's searches reach the higher one alone.
        let both = named(&[("FIRST", 5000), ("NEXT", 5001), ("INSIDE", 5030)]);
        assert!(both.site_above_start_pair("21:5000-5035", "INSIDE", "FIRST"));
        assert!(!both.site_above_start_pair("21:5000-5035", "INSIDE", "NEXT"));
        // A site on the second base alone inside a span, and a site one base above its end:
        // the reference's search above the start returns the outer site, tied with the
        // covered one at one base, and the candidate keeps that search's answer beside its
        // own, so the two name the same sites and no pair is taken.
        let outer = named(&[("SECOND", 5001), ("OUTER", 5011)]);
        assert!(!outer.site_above_start_pair("21:5000-5010", "SECOND,OUTER", "SECOND"));
        assert!(!outer.site_above_start_pair("21:5000-5010", "SECOND,OUTER", "SECOND,OUTER"));
        assert!(!outer.site_above_start_pair("21:5000-5010", "SECOND", "SECOND,OUTER"));
        // A higher site inside such a span, three bases below its end: the reference
        // names it, the candidate the site on the second base, one base from the start.
        let inside_too = named(&[("SECOND", 5001), ("INSIDE", 5007), ("OUTER", 5011)]);
        assert!(inside_too.site_above_start_pair("21:5000-5010", "INSIDE", "SECOND"));
    }

    #[test]
    fn the_last_position_has_no_base_above_it() {
        let named: IndexMap<String, BTreeSet<i64>> = [("BELOW", 1000), ("LAST", i64::MAX)]
            .iter()
            .map(|(name, site)| (name.to_string(), BTreeSet::from([*site])))
            .collect();
        let sites = NearestSites::on_one_sequence("21", &named);
        // A record on the last position has no base above its first base.
        let last = format!("21:{}", i64::MAX);
        assert!(!sites.site_above_start_pair(&last, "BELOW", "LAST"));
        // One base below it the reference's search above the start has no position to
        // start from, so it reaches the lower site alone.
        let below = format!("21:{}", i64::MAX - 1);
        assert!(sites.site_above_start_pair(&below, "BELOW", "LAST"));
        assert!(!sites.site_above_start_pair(&below, "LAST", "LAST"));
    }

    #[test]
    fn sites_of_a_cache_drop_nulls_and_refuse_a_coordinate_that_is_not_an_integer() {
        use crate::cache::testing::{one_shard_cache, scratch};
        let root = scratch("nearest-sites");
        let plus = r#"{"biotype": "protein_coding", "strand": 1, "start": 10, "end": 20, "stable_id": "P", "gene_stable_id": "G", "gene_symbol": "S"}"#;
        let minus = r#"{"biotype": "protein_coding", "strand": -1, "start": "30", "end": "40", "stable_id": "M", "gene_stable_id": "G"}"#;
        let other = r#"{"biotype": "lncRNA", "strand": 1, "stable_id": "L"}"#;
        let expected = sites(&[("P", &[10]), ("M", &[40]), ("G", &[10, 40]), ("S", &[10])]).sites;
        let array = one_shard_cache(
            &root.join("array"),
            &format!("[null, {plus}, null, {minus}, {other}]"),
        );
        assert_eq!(NearestSites::from_cache(&array).unwrap().sites, expected);
        let wrapped = one_shard_cache(
            &root.join("map"),
            &format!(
                r#"{{"1-1000000": [null, {plus}], "1000001-2000000": [{minus}, null, {other}]}}"#
            ),
        );
        let read = NearestSites::from_cache(&wrapped).unwrap();
        assert_eq!(read.sites, expected);
        assert_eq!(read.len(), 4);
        let sequences: Vec<(&String, &BTreeSet<i64>)> = read.by_sequence.iter().collect();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].1, &BTreeSet::from([10, 40]));
        let without_an_integer_site = [
            (
                r#"{"biotype": "protein_coding", "strand": 1, "end": 20, "stable_id": "P"}"#,
                "start",
            ),
            (
                r#"{"biotype": "protein_coding", "strand": -1, "start": 10, "stable_id": "P"}"#,
                "end",
            ),
            (
                r#"{"biotype": "protein_coding", "strand": 1, "start": "abc", "end": 20, "stable_id": "P"}"#,
                "start",
            ),
            (
                r#"{"biotype": "protein_coding", "strand": 1, "start": 12.5, "end": 20, "stable_id": "P"}"#,
                "start",
            ),
            (
                r#"{"biotype": "protein_coding", "strand": 1, "start": true, "end": 20, "stable_id": "P"}"#,
                "start",
            ),
            (
                r#"{"biotype": "protein_coding", "strand": 1, "start": null, "end": 20, "stable_id": "P"}"#,
                "start",
            ),
        ];
        for (n, (transcript, key)) in without_an_integer_site.iter().enumerate() {
            let cache = one_shard_cache(&root.join(format!("bad{n}")), &format!("[{transcript}]"));
            let refused = NearestSites::from_cache(&cache).unwrap_err().to_string();
            assert!(
                refused.ends_with(&format!(
                    "1-1000000.json: transcript 'P' has no integer {key}"
                )),
                "{refused}"
            );
        }
        let refused = NearestSites::from_cache(&root.join("absent"))
            .unwrap_err()
            .to_string();
        assert!(
            refused.starts_with("--nearest-sites: ")
                && refused.ends_with("absent has no transcripts/ directory")
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_coordinate_is_a_whole_number_or_a_text_of_decimal_digits() {
        let shard = Path::new("c/transcripts/1/1-5.json");
        let coordinate = |text: &str| {
            NearestSites::coordinate(shard, &pyjson::loads(text).unwrap(), "start")
                .map_err(|refused| refused.to_string())
        };
        assert_eq!(coordinate(r#"{"start": 10}"#), Ok(10));
        assert_eq!(coordinate(r#"{"start": 0}"#), Ok(0));
        assert_eq!(coordinate(r#"{"start": "007"}"#), Ok(7));
        assert_eq!(coordinate("{\"start\": \"\u{663}\u{660}\"}"), Ok(30));
        for (text, id) in [
            (r#"{"start": -5, "stable_id": "P"}"#, "'P'"),
            (r#"{"start": "+5", "stable_id": "P"}"#, "'P'"),
            (r#"{"start": "", "stable_id": "P"}"#, "'P'"),
            (r#"{"start": 1.0}"#, "None"),
            (r#"{"start": [1], "stable_id": null}"#, "None"),
            (
                r#"{"end": 1, "stable_id": ["a", 1.5, true]}"#,
                "['a', 1.5, True]",
            ),
            (
                r#"{"start": 99999999999999999999, "stable_id": "it's"}"#,
                "\"it's\"",
            ),
        ] {
            assert_eq!(
                coordinate(text),
                Err(format!(
                    "c/transcripts/1/1-5.json: transcript {id} has no integer start"
                )),
                "{text}"
            );
        }
    }
}
