// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The coding ends of the transcripts of a candidate JSON cache: for each, whether its
//! coding sequence ends inside a codon, and the stop codon and the start codon vep-rs reads
//! at its two ends.

use crate::argparse::path_text;
use crate::cache;
use crate::coding_ends::{coding_ends_with_three_prime_utr, CodingEnds, MapperFacts, RecordFacts};
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::genetic_codes::CodonTable;
use crate::pyjson::Value;
use crate::rules::downstream::codon_table_id;
use crate::rules::intron_ends::VERSION_SUFFIX_RE;
use crate::rules::transcript_model::{cache_integer, TranscriptModel};
use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::hash_map;
use std::collections::HashMap;
use std::path::Path;
use std::path::PathBuf;

/// What the cache's records of one stable id show.
#[derive(Debug, Clone)]
struct Kept {
    /// A record of the id holds a coding sequence of a length three does not divide.
    partial_last_codon: bool,
    /// The coding ends of each record of the id, a reading that two records share once:
    /// a record that two shards repeat is one, a copy of the transcript on another span or
    /// another transcript under the same id is a second.
    records: Vec<CodingEnds>,
}

/// The lowest and the highest genomic base of a Location, `None` for a text that names
/// no position after its sequence name.
pub(crate) fn location_span(location: &str) -> Option<(i64, i64)> {
    let (_, position) = location.rsplit_once(':')?;
    let (first, last) = match position.split_once('-') {
        Some((first, last)) => (first.parse::<i64>().ok()?, last.parse::<i64>().ok()?),
        None => {
            let base = position.parse::<i64>().ok()?;
            (base, base)
        }
    };
    Some((first.min(last), first.max(last)))
}

/// The coding ends of the transcripts a candidate JSON cache holds, by stable id
/// (`CodingEnds`), and with them the transcripts whose coding sequence is no whole number
/// of codons: a transcript's `variation_effect_feature_cache.translateable_seq` is text
/// whose length three does not divide. Such a coding sequence ends inside a codon and no
/// stop codon closes it, whatever its last three bases spell; no row of the transcript need
/// say so, since Ensembl VEP writes `cds_end_NF` from the transcript's attributes and
/// `incomplete_terminal_codon_variant` only for an allele written in the partial codon. A
/// chromosome's shards are read on the first lookup of an entry whose Location names its
/// directory, every shard on the first lookup of an entry without one (a VCF or JSON
/// entry).
#[derive(Debug)]
pub struct PartialCodonModels {
    root: PathBuf,
    read: RefCell<IndexMap<Option<String>, HashMap<String, Kept>>>,
}

impl PartialCodonModels {
    pub fn new(cache: &Path) -> Result<Self> {
        let root = cache::joined(cache, "transcripts");
        if !root.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        Ok(PartialCodonModels {
            root,
            read: RefCell::new(IndexMap::new()),
        })
    }

    fn models(&self, chromosome: Option<&str>) -> Result<HashMap<String, Kept>> {
        let shards = cache::shards(&self.root, chromosome)?;
        let read = cache::map_shards(&shards, cache::threads(), |_, transcripts| {
            Ok(transcripts
                .iter()
                .filter_map(coding_ends_of)
                .collect::<Vec<_>>())
        })?;
        let mut models: HashMap<String, Kept> = HashMap::new();
        for (stable_id, ends) in read.into_iter().flatten() {
            match models.entry(stable_id) {
                hash_map::Entry::Vacant(place) => {
                    place.insert(Kept {
                        partial_last_codon: ends.tail_length != 0,
                        records: vec![ends],
                    });
                }
                hash_map::Entry::Occupied(mut place) => {
                    let kept = place.get_mut();
                    kept.partial_last_codon |= ends.tail_length != 0;
                    if !kept.records.contains(&ends) {
                        kept.records.push(ends);
                    }
                }
            }
        }
        Ok(models)
    }

    /// What the cache shows of the transcript `feature`, under the Feature as the row
    /// writes it and under the Feature without a version suffix (`--transcript_version`),
    /// each `None` for an id the cache does not hold, handed to `answer`. The transcripts
    /// are those of the sequence `location` names, or of every sequence without a Location
    /// and for a Location that names no directory of the cache.
    fn kept<T>(
        &self,
        feature: &str,
        location: Option<&str>,
        answer: impl FnOnce(Option<&Kept>, Option<&Kept>) -> T,
    ) -> Result<T> {
        let mut chromosome = location
            .filter(|location| !location.is_empty())
            .map(|location| location.split(':').next().unwrap_or("").to_string());
        if chromosome
            .as_deref()
            .is_some_and(|name| !cache::joined(&self.root, name).is_dir())
        {
            chromosome = None;
        }
        if !self.read.borrow().contains_key(&chromosome) {
            let models = self.models(chromosome.as_deref())?;
            self.read.borrow_mut().insert(chromosome.clone(), models);
        }
        let read = self.read.borrow();
        let models = &read[&chromosome];
        Ok(answer(
            models.get(feature),
            models.get(VERSION_SUFFIX_RE.replace_all(feature, "").as_ref()),
        ))
    }

    /// Whether the entry's transcript is one whose coding sequence ends inside a codon; a
    /// Feature carrying a version the cache's stable id lacks (`--transcript_version`) is
    /// looked up without it. A transcript the cache does not hold is none.
    pub fn holds(&self, e: &Entry) -> Result<bool> {
        self.kept(&e.feature, e.location.as_deref(), |written, unversioned| {
            [written, unversioned]
                .into_iter()
                .flatten()
                .any(|kept| kept.partial_last_codon)
        })
    }

    /// The coding ends of the entry's transcript: of the record of its Feature, on the
    /// sequence its Location names, whose span shares a base with the Location, an
    /// insertion's being the two bases it stands between. A VCF or JSON entry is looked up
    /// at the Location of its allele, and a Feature carrying a version the cache's stable
    /// id lacks without it. `None` for a transcript the cache does not hold, and where no
    /// record of the stable id holds the Location or more than one does: one id stands on
    /// several records of a sequence, and a row names one of them by where it lies. An
    /// entry without a position is answered where the id names one record.
    pub fn coding_ends(&self, e: &Entry) -> Result<Option<CodingEnds>> {
        let location = e
            .location
            .as_deref()
            .filter(|location| !location.is_empty())
            .or(e.allele_location.as_deref());
        let span = location.and_then(location_span);
        self.kept(&e.feature, location, |written, unversioned| {
            let records = &written.or(unversioned)?.records;
            let mut holding = records.iter().filter(|ends| match (span, ends.span) {
                (None, _) => true,
                (Some((first, last)), Some((start, end))) => first <= end && start <= last,
                (Some(_), None) => false,
            });
            match (holding.next(), holding.next()) {
                (Some(ends), None) => Some(ends.clone()),
                _ => None,
            }
        })
    }
}

/// A member that is text, `Some(None)` for one the record lacks or holds as `null`; `None`
/// for any other value, which vep-rs refuses.
fn text(member: Option<&Value>) -> Option<Option<&str>> {
    match member {
        None | Some(Value::Null) => Some(None),
        Some(Value::Str(text)) => Some(Some(text)),
        Some(_) => None,
    }
}

/// The facts of a record whose coding sequence is `coding`, read into its coding ends;
/// `None` for a record that states one of them in a form vep-rs reads otherwise or
/// refuses: an attribute whose `code` or `value` is no text, a `ccds` or a `peptide` that is
/// neither text nor `null`, a `codon_table` no genetic code is held for, a `start_phase`
/// or a coding bound of the mapper that is no small whole number, a mapper beside which
/// the record states no `TranscriptModel`, or `introns` that are no list of two whole
/// numbers each. The introns are the cached ones, or where the record caches none the
/// stretches between its exons.
fn read_coding_ends(
    transcript: &Value,
    cached: &IndexMap<String, Value>,
    coding: &str,
) -> Option<CodingEnds> {
    let attribute_codes: Vec<&str> = match transcript.get("attributes") {
        None => Vec::new(),
        Some(attributes) => attributes
            .as_list()?
            .iter()
            .map(|attribute| {
                attribute.get("value")?.as_str()?;
                attribute.get("code")?.as_str()
            })
            .collect::<Option<_>>()?,
    };
    let CodonTable::Held(codon_table) =
        codon_table_id(cached.get("codon_table").unwrap_or(&Value::Null))
    else {
        return None;
    };
    let span = record_span(transcript);
    let mut introns: Vec<(i64, i64)> = match cached.get("introns") {
        None | Some(Value::Null) => Vec::new(),
        Some(introns) => introns
            .as_list()?
            .iter()
            .map(|intron| {
                Some((
                    intron.get("start").and_then(cache_integer)?,
                    intron.get("end").and_then(cache_integer)?,
                ))
            })
            .collect::<Option<_>>()?,
    };
    let model;
    let mapper = match cached.get("mapper") {
        None | Some(Value::Null) => None,
        Some(mapper) => {
            mapper.as_dict()?;
            let bound = |name: &str| match mapper.get(name) {
                None | Some(Value::Null) => Some(0),
                Some(bound) => cache_integer(bound).filter(|bound| *bound >= 0),
            };
            let start_phase = match mapper.get("start_phase") {
                None => 0,
                Some(phase) => cache_integer(phase).filter(|phase| i8::try_from(*phase).is_ok())?,
            };
            model = TranscriptModel::from_cache(transcript, start_phase.max(0))?;
            if introns.is_empty() {
                let mut exons = model.exons.clone();
                exons.sort_unstable();
                introns = exons
                    .windows(2)
                    .filter_map(|exons| {
                        Some((exons[0].1.checked_add(1)?, exons[1].0.checked_sub(1)?))
                    })
                    .collect();
            }
            Some(MapperFacts {
                start_phase,
                coding_start: bound("cdna_coding_start")?,
                coding_end: bound("cdna_coding_end")?,
                pairs: &model.pairs,
            })
        }
    };
    let three_prime_utr = text(cached.get("three_prime_utr"))?;
    Some(coding_ends_with_three_prime_utr(
        &RecordFacts {
            attribute_codes: &attribute_codes,
            consensus_id: text(transcript.get("ccds"))?,
            strand: transcript
                .get("strand")
                .and_then(cache_integer)
                .filter(|strand| matches!(strand, 1 | -1)),
            codon_table: Some(codon_table),
            coding: Some(coding),
            peptide: text(cached.get("peptide"))?,
            mapper,
            span,
            introns: &introns,
        },
        three_prime_utr,
    ))
}

/// The record's `start` and `end`, `None` where it states either in no form an integer
/// is read from.
fn record_span(transcript: &Value) -> Option<(i64, i64)> {
    Some((
        transcript.get("start").and_then(cache_integer)?,
        transcript.get("end").and_then(cache_integer)?,
    ))
}

/// The stable id of one cache transcript object and its coding ends, read as vep-rs loads
/// the record; `None` for a transcript without a stable id that is text. The mapper pairs
/// are the transcript's `TranscriptModel`'s. A transcript without a coding sequence that is
/// text is not read, nor is one whose record `read_coding_ends` refuses.
pub fn coding_ends_of(transcript: &Value) -> Option<(String, CodingEnds)> {
    let stable_id = transcript.get("stable_id")?.as_str()?.to_string();
    let cached = transcript
        .get("variation_effect_feature_cache")
        .and_then(Value::as_dict);
    let coding = cached
        .and_then(|cached| cached.get("translateable_seq"))
        .and_then(Value::as_str);
    let unread = |tail_length: u8| CodingEnds {
        span: record_span(transcript),
        ..CodingEnds::unread(tail_length)
    };
    let ends = match (cached, coding) {
        (Some(cached), Some(coding)) => read_coding_ends(transcript, cached, coding)
            .unwrap_or_else(|| unread((coding.chars().count() % 3) as u8)),
        _ => unread(0),
    };
    Some((stable_id, ends))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::testing::{plant, scratch};
    use crate::classes::HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE;
    use crate::rules::{field_class_accepts, RuleContext, Side};

    /// A reference row of a frameshift on `feature` at chromosome 3, its FLAGS empty and
    /// its terms without `incomplete_terminal_codon_variant`.
    fn frameshift_row(location: Option<&str>, feature: &str) -> Entry {
        Entry {
            key: "k".to_string(),
            location: location.map(str::to_string),
            feature: feature.to_string(),
            consequence: "frameshift_variant".to_string(),
            ..Entry::default()
        }
    }

    /// A cache of one chromosome: `ENST1` ends two bases into a codon, `ENST2` one base,
    /// `ENST3` is whole codons, `ENST4` carries no coding sequence and `ENST5` a coding
    /// sequence that is no text.
    fn planted(name: &str) -> PathBuf {
        let root = scratch(name);
        plant(
            &root,
            "transcripts/3/1-5.json",
            r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"translateable_seq": "NNNNNNNN"}},
                {"stable_id": "ENST2", "variation_effect_feature_cache": {"translateable_seq": "NNNN"}},
                {"stable_id": "ENST3", "variation_effect_feature_cache": {"translateable_seq": "NNNNNN"}},
                {"stable_id": "ENST4", "variation_effect_feature_cache": {}},
                {"stable_id": "ENST5", "variation_effect_feature_cache": {"translateable_seq": 8}},
                {"stable_id": "ENST6"},
                {"variation_effect_feature_cache": {"translateable_seq": "NNNN"}},
                null]"#,
        );
        root
    }

    #[test]
    fn a_coding_sequence_that_ends_inside_a_codon_is_read_from_the_cache() {
        let root = planted("partial-codon-models");
        plant(
            &root,
            "transcripts/21/1-5.json",
            r#"{"1-5": [{"stable_id": "ENST7", "variation_effect_feature_cache": {"translateable_seq": "NNNNN"}}]}"#,
        );
        let models = PartialCodonModels::new(&root).unwrap();
        let holds = |location: Option<&str>, feature: &str| {
            models
                .holds(&frameshift_row(location, feature))
                .map_err(|refused| refused.to_string())
        };
        for feature in ["ENST1", "ENST2", "ENST1.4"] {
            assert_eq!(holds(Some("3:7-8"), feature), Ok(true), "{feature}");
        }
        for feature in [
            "ENST3", "ENST4", "ENST5", "ENST6", "ENST7", "ENST9", "ENST1.",
        ] {
            assert_eq!(holds(Some("3:7-8"), feature), Ok(false), "{feature}");
        }
        assert_eq!(holds(Some("21:7"), "ENST7"), Ok(true));
        // A row without a Location, or one whose Location names no directory of the cache,
        // reads every chromosome.
        for location in [None, Some(""), Some("nowhere:1")] {
            assert_eq!(holds(location, "ENST1"), Ok(true), "{location:?}");
            assert_eq!(holds(location, "ENST7"), Ok(true), "{location:?}");
            assert_eq!(holds(location, "ENST3"), Ok(false), "{location:?}");
        }
        let refused = PartialCodonModels::new(&root.join("absent"))
            .unwrap_err()
            .to_string();
        assert!(
            refused.starts_with("--json-cache: ")
                && refused.ends_with("absent has no transcripts/ directory"),
            "{refused}"
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_stop_named_on_a_transcript_the_cache_shows_incomplete_reads_under_the_class() {
        let root = planted("partial-codon-class");
        let cached = Side {
            partial_codon_models: Some(PartialCodonModels::new(&root).unwrap()),
            ..Side::default()
        };
        let bare = Side::default();
        let accepts = |side: &Side, feature: &str, reference: &str, candidate: &str| {
            let row = frameshift_row(Some("3:7-8"), feature);
            let context = RuleContext {
                record_rows: None,
                candidate_rows: None,
                declared: &[],
                assembly: None,
                dbnsfp_columns: &[],
                side,
            };
            field_class_accepts(
                HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE,
                &row,
                "HGVSp",
                reference,
                candidate,
                None,
                &context,
            )
            .unwrap()
        };
        // The four ways a description names the stop alone, on a row that shows nothing.
        for description in [
            "ENSP1.2:p.Ter612%3D",
            "ENSP1.2:p.Ter612=",
            "ENSP1.2:p.Ter612CysextTer1",
            "ENSP1.2:p.Ter107ArgextTer?",
            "ENSP1.2:p.Ter68del",
            "ENSP1.2:p.Ter68Gly",
            "ENSP1.2:p.(Ter612=)",
        ] {
            assert!(accepts(&cached, "ENST1", description, "-"), "{description}");
            assert!(accepts(&cached, "ENST2", description, "-"), "{description}");
            // No cache, a transcript of whole codons and one the cache lacks.
            assert!(!accepts(&bare, "ENST1", description, "-"), "{description}");
            assert!(
                !accepts(&cached, "ENST3", description, "-"),
                "{description}"
            );
            assert!(
                !accepts(&cached, "ENST9", description, "-"),
                "{description}"
            );
        }
        // Another description, a candidate that writes a notation and the sides exchanged
        // stay differences on a transcript the cache shows incomplete.
        assert!(!accepts(
            &cached,
            "ENST1",
            "ENSP1.2:p.Lys174ArgfsTer12",
            "-"
        ));
        assert!(!accepts(
            &cached,
            "ENST1",
            "ENSP1.2:p.Ter612=",
            "ENSP1.2:p.Ter612="
        ));
        assert!(!accepts(&cached, "ENST1", "-", "ENSP1.2:p.Ter612="));
        assert!(!accepts(
            &cached,
            "ENST1",
            "ENSP1.2:p.Ter612=",
            "ENSP1.2:p.Lys612Ter"
        ));
        std::fs::remove_dir_all(root).unwrap();
    }

    /// A start codon, a codon of alanine and a stop codon of the standard code.
    const CLOSED: [&str; 3] = ["ATG", "GCT", "TAA"];

    /// A cache record of one exon on the forward strand whose cDNA is its coding sequence
    /// `coding`, its first base the genomic base 101 and its span that exon, with a cached
    /// peptide that opens with methionine. `top` and `cached` are further members of the record and of its
    /// `variation_effect_feature_cache`, each opening with a comma; a member named there
    /// replaces the one written here. With `quoted` every number is the text a converted
    /// cache writes.
    fn record(stable_id: &str, coding: &str, top: &str, cached: &str, quoted: bool) -> String {
        let number = |value: usize| match quoted {
            true => format!("\"{value}\""),
            false => value.to_string(),
        };
        let (one, last, first_base, last_base) = (
            number(1),
            number(coding.len()),
            number(101),
            number(100 + coding.len()),
        );
        format!(
            r#"{{"stable_id": "{stable_id}", "strand": {one}, "cdna_coding_start": {one},
                "start": {first_base}, "end": {last_base},
                "cdna_coding_end": {last}, "exons": [{{"start": {first_base}, "end": {last_base}}}]
                {top}, "variation_effect_feature_cache": {{"translateable_seq": "{coding}",
                "peptide": "MA", "mapper": {{"start_phase": {zero}, "cdna_coding_start": {one},
                "cdna_coding_end": {last}, "pairs": [{{"from_start": {one}, "from_end": {last},
                "to_start": {first_base}, "to_end": {last_base}, "ori": {one}}}]}}{cached}}}}}"#,
            zero = number(0)
        )
    }

    /// The models of a cache whose chromosome 3 holds `records` in one shard.
    fn models_of(name: &str, records: &[String]) -> (PathBuf, PartialCodonModels) {
        let root = scratch(name);
        plant(
            &root,
            "transcripts/3/1-5.json",
            &format!("[{}]", records.join(",")),
        );
        let models = PartialCodonModels::new(&root).unwrap();
        (root, models)
    }

    fn ends_of(models: &PartialCodonModels, feature: &str) -> Option<CodingEnds> {
        models
            .coding_ends(&frameshift_row(Some("3:102-103"), feature))
            .unwrap()
    }

    /// The ends of a record `record` makes that is not read, whose coding sequence has
    /// `coding_length` bases.
    fn unread(coding_length: i64) -> CodingEnds {
        CodingEnds {
            span: Some((101, 100 + coding_length)),
            ..CodingEnds::unread((coding_length % 3) as u8)
        }
    }

    #[test]
    fn a_flag_is_read_from_an_attribute_code_whatever_the_value_beside_it() {
        use crate::coding_ends::{StartEnd, StopEnd};
        let closed = CLOSED.concat();
        let attributes = |items: &str| format!(r#", "attributes": [{items}]"#);
        let flag = |code: &str, value: &str| format!(r#"{{"code": "{code}", "value": "{value}"}}"#);
        let made = |stable_id: &str, items: &str| {
            record(stable_id, &closed, &attributes(items), "", false)
        };
        let (root, models) = models_of(
            "coding-ends-flags",
            &[
                made("END1", &flag("cds_end_NF", "1")),
                made("END0", &flag("cds_end_NF", "0")),
                made("START1", &flag("cds_start_NF", "1")),
                made("START0", &flag("cds_start_NF", "0")),
                made(
                    "BOTH",
                    &[
                        flag("TSL", "tsl1"),
                        flag("cds_start_NF", "0"),
                        flag("cds_end_NF", ""),
                    ]
                    .join(","),
                ),
                made("OTHER", &flag("TSL", "tsl1")),
                made("EMPTY", ""),
                record("BARE", &closed, "", "", false),
                // An attribute in a form the loader refuses leaves the record unread.
                made("NO_VALUE", r#"{"code": "cds_end_NF"}"#),
                made("NO_CODE", r#"{"value": "1"}"#),
                record("NO_LIST", &closed, r#", "attributes": {}"#, "", false),
            ],
        );
        let flagged_end = StopEnd::NoneShown { flagged: true };
        let flagged_start = StartEnd::NotShown {
            flagged: true,
            padded: false,
        };
        for (feature, stop, start) in [
            ("END1", flagged_end, StartEnd::Shown),
            ("END0", flagged_end, StartEnd::Shown),
            ("START1", StopEnd::Whole, flagged_start),
            ("START0", StopEnd::Whole, flagged_start),
            ("BOTH", flagged_end, flagged_start),
            ("OTHER", StopEnd::Whole, StartEnd::Shown),
            ("EMPTY", StopEnd::Whole, StartEnd::Shown),
            ("BARE", StopEnd::Whole, StartEnd::Shown),
            ("NO_VALUE", StopEnd::Unread, StartEnd::Unread),
            ("NO_CODE", StopEnd::Unread, StartEnd::Unread),
            ("NO_LIST", StopEnd::Unread, StartEnd::Unread),
        ] {
            let ends = ends_of(&models, feature).expect(feature);
            assert_eq!((ends.stop, ends.start), (stop, start), "{feature}");
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_consensus_id_is_text_that_is_neither_empty_nor_a_dash() {
        use crate::coding_ends::{StartEnd, StopEnd};
        // A coding sequence that ends one base into a stop codon on the transcript's last
        // base, with a cached peptide that opens with alanine on a codon of alanine.
        let begun = ["GCT", "GCT", "T"].concat();
        let made = |stable_id: &str, top: &str| {
            record(stable_id, &begun, top, r#", "peptide": "AA""#, false)
        };
        let (root, models) = models_of(
            "coding-ends-consensus",
            &[
                made("ID", r#", "ccds": "CCDS1.1""#),
                made("EMPTY", r#", "ccds": """#),
                made("DASH", r#", "ccds": "-""#),
                made("NULL", r#", "ccds": null"#),
                made("ABSENT", ""),
                made("NUMBER", r#", "ccds": 5"#),
            ],
        );
        let not_shown = StartEnd::NotShown {
            flagged: false,
            padded: false,
        };
        for (feature, shown) in [
            ("ID", true),
            ("EMPTY", false),
            ("DASH", false),
            ("NULL", false),
            ("ABSENT", false),
        ] {
            let ends = ends_of(&models, feature).expect(feature);
            let stop = match shown {
                true => StopEnd::BegunByTheTail { bases: 1 },
                false => StopEnd::NoneShown { flagged: false },
            };
            assert_eq!(ends.stop, stop, "{feature}");
            let start = if shown { StartEnd::Shown } else { not_shown };
            assert_eq!(ends.start, start, "{feature}");
            assert_eq!(ends.tail_bases, [Some(b'T'), None], "{feature}");
            assert_eq!(ends.tail, [Some(107), None], "{feature}");
        }
        assert_eq!(ends_of(&models, "NUMBER"), Some(unread(7)));
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_record_is_read_alike_as_numbers_and_as_the_text_of_a_converted_cache() {
        use crate::coding_ends::{StartEnd, StopEnd};
        // A tail that begins a stop codon of the vertebrate mitochondrial code alone.
        let begun = ["ATG", "GCT", "AG"].concat();
        let (root, models) = models_of(
            "coding-ends-forms",
            &[
                record("NUMBERS", &begun, "", r#", "codon_table": 2"#, false),
                record("TEXT", &begun, "", r#", "codon_table": "2""#, true),
                record("STANDARD", &begun, "", "", false),
                record("ZERO", &begun, "", r#", "codon_table": 0"#, false),
                record("NO_TABLE_HELD", &begun, "", r#", "codon_table": 7"#, false),
                record("NO_NUMBER", &begun, "", r#", "codon_table": "two""#, false),
                record(
                    "NO_MAPPER",
                    &begun,
                    "",
                    r#", "codon_table": 2, "mapper": null"#,
                    false,
                ),
                record("NO_EXONS", &begun, r#", "exons": []"#, "", false),
                record("NO_PEPTIDE", &begun, "", r#", "peptide": null"#, false),
                record("PEPTIDE_NO_TEXT", &begun, "", r#", "peptide": 7"#, false),
                record(
                    "PHASE_NO_NUMBER",
                    &begun,
                    "",
                    r#", "mapper": {"start_phase": "one", "pairs": []}"#,
                    false,
                ),
            ],
        );
        let none = StopEnd::NoneShown { flagged: false };
        let shown = StopEnd::BegunByTheTail { bases: 2 };
        let numbers = ends_of(&models, "NUMBERS").unwrap();
        assert_eq!((numbers.stop, numbers.start), (shown, StartEnd::Shown));
        assert_eq!(numbers.strand, Some(1));
        assert_eq!(numbers.first_codon, [Some(101), Some(102), Some(103)]);
        assert_eq!(numbers.last_codon, [Some(106), Some(107), Some(108)]);
        assert_eq!(numbers.coding_region_end, Some((106, 108)));
        assert_eq!(numbers.tail, [Some(107), Some(108)]);
        assert_eq!(ends_of(&models, "TEXT"), Some(numbers.clone()));
        for feature in ["STANDARD", "ZERO"] {
            assert_eq!(ends_of(&models, feature).unwrap().stop, none, "{feature}");
        }
        // Without a mapper the transcript is not known to end with its coding sequence,
        // and no position is known.
        let bare = ends_of(&models, "NO_MAPPER").unwrap();
        assert_eq!((bare.stop, bare.start), (none, StartEnd::Shown));
        assert_eq!((bare.first_codon, bare.tail), ([None; 3], [None; 2]));
        assert_eq!(bare.tail_bases, [Some(b'A'), Some(b'G')]);
        // Without a cached peptide the first codon's own residue opens the translation.
        assert_eq!(
            ends_of(&models, "NO_PEPTIDE").unwrap().start,
            StartEnd::Shown
        );
        for feature in [
            "NO_TABLE_HELD",
            "NO_NUMBER",
            "NO_EXONS",
            "PEPTIDE_NO_TEXT",
            "PHASE_NO_NUMBER",
        ] {
            assert_eq!(ends_of(&models, feature), Some(unread(8)), "{feature}");
            // A record that is not read still ends inside a codon.
            assert!(models
                .holds(&frameshift_row(Some("3:7-8"), feature))
                .unwrap());
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_row_is_answered_by_its_sequence_and_its_feature() {
        use crate::coding_ends::StopEnd;
        let closed = CLOSED.concat();
        let open = ["ATG", "GCT", "GCT"].concat();
        let partial = ["ATG", "GCT", "G"].concat();
        let plain = |stable_id: &str, coding: &str| record(stable_id, coding, "", "", false);
        let (root, models) = models_of(
            "coding-ends-rows",
            &[
                plain("ENST1", &closed),
                plain("ENST2", &open),
                plain("NM_1.2", &closed),
                // Two records of one id on one span name no one transcript.
                plain("ENST3", &closed),
                plain("ENST3", &partial),
                // Two that read alike are one.
                plain("ENST4", &open),
                plain("ENST4", &open),
                r#"{"stable_id": "ENST5"}"#.to_string(),
            ],
        );
        // A second shard of the sequence repeats a record, and another sequence holds
        // another transcript under the same id.
        plant(
            &root,
            "transcripts/3/6-9.json",
            &format!("[{}]", plain("ENST1", &closed)),
        );
        plant(
            &root,
            "transcripts/21/1-5.json",
            &format!("[{}, {}]", plain("ENST1", &open), plain("ENST6", &closed)),
        );
        let stop = |location: Option<&str>, feature: &str| {
            models
                .coding_ends(&frameshift_row(location, feature))
                .unwrap()
                .map(|ends| ends.stop)
        };
        let none = StopEnd::NoneShown { flagged: false };
        for (location, feature, expected) in [
            (Some("3:104-105"), "ENST1", Some(StopEnd::Whole)),
            (Some("3:104-105"), "ENST1.4", Some(StopEnd::Whole)),
            (Some("3:104-105"), "ENST2", Some(none)),
            (Some("3:104-105"), "NM_1.2", Some(StopEnd::Whole)),
            (Some("3:104-105"), "NM_1", None),
            (Some("3:104-105"), "ENST3", None),
            (Some("3:104-105"), "ENST4", Some(none)),
            (Some("3:104-105"), "ENST6", None),
            (Some("3:104-105"), "ENST9", None),
            (Some("21:104"), "ENST1", Some(none)),
            (Some("21:104"), "ENST6", Some(StopEnd::Whole)),
            (Some("21:104"), "ENST2", None),
            // A Location is held where it shares a base with the record's span, 101 to
            // 109: a span that runs on past either end, and an insertion one of whose
            // two bases is the record's first or last.
            (Some("3:90-101"), "ENST1", Some(StopEnd::Whole)),
            (Some("3:109-140"), "ENST1", Some(StopEnd::Whole)),
            (Some("3:50-500"), "ENST1", Some(StopEnd::Whole)),
            (Some("3:100-101"), "ENST1", Some(StopEnd::Whole)),
            (Some("3:109-110"), "ENST1", Some(StopEnd::Whole)),
            (Some("3:7-8"), "ENST1", None),
            (Some("3:100"), "ENST1", None),
            (Some("3:110-111"), "ENST1", None),
            // A record that states no span holds no Location.
            (Some("3:104-105"), "ENST5", None),
            // A row that names no sequence of the cache reads every sequence, where one
            // id stands on two records that both hold the position.
            (Some("nowhere:104"), "ENST1", None),
            (Some("nowhere:104"), "ENST2", Some(none)),
            (Some("nowhere:7"), "ENST2", None),
            // A row without a position is answered where its id names one record.
            (None, "ENST1", None),
            (Some(""), "ENST1", None),
            (None, "ENST2", Some(none)),
            (None, "ENST5", Some(StopEnd::Unread)),
            (None, "ENST6", Some(StopEnd::Whole)),
            (None, "ENST6.1", Some(StopEnd::Whole)),
        ] {
            assert_eq!(stop(location, feature), expected, "{location:?} {feature}");
        }
        // A VCF or JSON entry is looked up at the Location of its allele.
        let entry = Entry {
            allele_location: Some("21:103-104".to_string()),
            ..frameshift_row(None, "ENST1")
        };
        assert_eq!(
            models.coding_ends(&entry).unwrap().map(|ends| ends.stop),
            Some(none)
        );
        // The question of the incomplete coding sequence is answered from the same records:
        // one record of the id that ends inside a codon is enough.
        let holds = |feature: &str| {
            models
                .holds(&frameshift_row(Some("3:7-8"), feature))
                .unwrap()
        };
        assert!(holds("ENST3"));
        for feature in ["ENST1", "ENST2", "ENST4", "ENST5", "ENST9"] {
            assert!(!holds(feature), "{feature}");
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn one_id_on_two_records_of_a_sequence_is_answered_by_the_record_that_holds_the_location() {
        use crate::coding_ends::StopEnd;
        // The id's first record is closed by a stop codon on the bases 101 to 109, and
        // its second, whose last codon names a residue, lies 400 bases further on.
        let closed = record("ENST1", &CLOSED.concat(), "", "", false);
        let further_on = record("ENST1", &["ATG", "GCT", "GCT"].concat(), "", "", false)
            .replace("101", "501")
            .replace("109", "509");
        let (root, models) = models_of("coding-ends-two-records", &[closed, further_on]);
        let stop = |location: &str| {
            models
                .coding_ends(&frameshift_row(Some(location), "ENST1"))
                .unwrap()
                .map(|ends| ends.stop)
        };
        assert_eq!(stop("3:104"), Some(StopEnd::Whole));
        assert_eq!(stop("3:109-110"), Some(StopEnd::Whole));
        assert_eq!(stop("3:504"), Some(StopEnd::NoneShown { flagged: false }));
        assert_eq!(
            stop("3:500-501"),
            Some(StopEnd::NoneShown { flagged: false })
        );
        // A Location between the two records is held by neither, and one over both by
        // both.
        assert_eq!(stop("3:300"), None);
        assert_eq!(stop("3:105-505"), None);
        std::fs::remove_dir_all(root).unwrap();
    }
}
