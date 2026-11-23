// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A second reading of a transcript shard, for what the cache reader of
//! `vep-seeded-variants` does not give: the record's `dbID`, whether the exons, in the
//! order the cache lists them, run with position along the strand, and the first residue
//! of the cached peptide.

use std::fmt;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use serde::de::{self, Deserializer, IgnoredAny, MapAccess, SeqAccess, Visitor};
use serde::Deserialize;

/// What a shard states of one transcript beyond the cache reader's `Transcript`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extra {
    pub stable_id: String,
    /// The transcript's start as the cache reader gives it: 1 for a start stored
    /// below 1.
    pub start: u64,
    pub end: u64,
    /// The record's `dbID`, as the cache reader reads it: `None` for a record that
    /// carries none, a null or an empty string.
    pub db_id: Option<i64>,
    /// An exon, in the cache's order, starts on the wrong side of the exon listed
    /// before it: before it on the forward strand, after it on the reverse.
    pub exons_against_position: bool,
    /// The first residue of the cached peptide, `None` where the cache carries none.
    pub peptide_first: Option<u8>,
}

/// A whole number as a JSON number or as a string that spells one; any other value
/// reads as none.
struct Whole(Option<i64>);

impl<'de> Deserialize<'de> for Whole {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Whole, D::Error> {
        struct Shape;

        impl<'de> Visitor<'de> for Shape {
            type Value = Whole;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("any JSON value")
            }

            fn visit_i64<E: de::Error>(self, value: i64) -> Result<Whole, E> {
                Ok(Whole(Some(value)))
            }

            fn visit_u64<E: de::Error>(self, value: u64) -> Result<Whole, E> {
                Ok(Whole(i64::try_from(value).ok()))
            }

            fn visit_f64<E: de::Error>(self, value: f64) -> Result<Whole, E> {
                let whole = value.fract() == 0.0 && value.abs() < 9.0e15;
                Ok(Whole(whole.then_some(value as i64)))
            }

            fn visit_str<E: de::Error>(self, value: &str) -> Result<Whole, E> {
                Ok(Whole(value.trim().parse().ok()))
            }

            fn visit_bool<E: de::Error>(self, _: bool) -> Result<Whole, E> {
                Ok(Whole(None))
            }

            fn visit_unit<E: de::Error>(self) -> Result<Whole, E> {
                Ok(Whole(None))
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Whole, A::Error> {
                while seq.next_element::<IgnoredAny>()?.is_some() {}
                Ok(Whole(None))
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Whole, A::Error> {
                while map.next_entry::<IgnoredAny, IgnoredAny>()?.is_some() {}
                Ok(Whole(None))
            }
        }

        deserializer.deserialize_any(Shape)
    }
}

/// The first byte of a JSON string, read without keeping the string; any other value
/// reads as none.
struct FirstByte(Option<u8>);

impl<'de> Deserialize<'de> for FirstByte {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<FirstByte, D::Error> {
        struct Shape;

        impl<'de> Visitor<'de> for Shape {
            type Value = FirstByte;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("any JSON value")
            }

            fn visit_str<E: de::Error>(self, value: &str) -> Result<FirstByte, E> {
                Ok(FirstByte(value.bytes().next()))
            }

            fn visit_i64<E: de::Error>(self, _: i64) -> Result<FirstByte, E> {
                Ok(FirstByte(None))
            }

            fn visit_u64<E: de::Error>(self, _: u64) -> Result<FirstByte, E> {
                Ok(FirstByte(None))
            }

            fn visit_f64<E: de::Error>(self, _: f64) -> Result<FirstByte, E> {
                Ok(FirstByte(None))
            }

            fn visit_bool<E: de::Error>(self, _: bool) -> Result<FirstByte, E> {
                Ok(FirstByte(None))
            }

            fn visit_unit<E: de::Error>(self) -> Result<FirstByte, E> {
                Ok(FirstByte(None))
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<FirstByte, A::Error> {
                while seq.next_element::<IgnoredAny>()?.is_some() {}
                Ok(FirstByte(None))
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<FirstByte, A::Error> {
                while map.next_entry::<IgnoredAny, IgnoredAny>()?.is_some() {}
                Ok(FirstByte(None))
            }
        }

        deserializer.deserialize_any(Shape)
    }
}

#[derive(Deserialize)]
struct RawExon {
    #[serde(default)]
    start: Option<Whole>,
}

#[derive(Deserialize)]
struct RawFeatureCache {
    #[serde(default)]
    peptide: Option<FirstByte>,
}

/// A transcript record, down to the keys this reading keeps.
#[derive(Deserialize)]
struct RawRecord {
    #[serde(default)]
    stable_id: Option<String>,
    #[serde(default, rename = "dbID")]
    db_id: Option<Whole>,
    #[serde(default)]
    start: Option<Whole>,
    #[serde(default)]
    end: Option<Whole>,
    #[serde(default)]
    strand: Option<Whole>,
    #[serde(default)]
    exons: Option<Vec<Option<RawExon>>>,
    #[serde(default, rename = "variation_effect_feature_cache")]
    feature_cache: Option<RawFeatureCache>,
}

impl RawRecord {
    /// The record's facts, `None` for a record without an identifier, a start or an
    /// end above zero, which the cache reader refuses or never gives.
    fn into_extra(self) -> Option<Extra> {
        let stable_id = self.stable_id?;
        let stored_start = self.start.and_then(|start| start.0)?;
        let end = u64::try_from(self.end.and_then(|end| end.0)?).ok()?;
        let reverse = self.strand.and_then(|strand| strand.0).unwrap_or(1) < 0;
        let starts: Vec<i64> = self
            .exons
            .unwrap_or_default()
            .into_iter()
            .flatten()
            .filter_map(|exon| exon.start.and_then(|start| start.0))
            .collect();
        let exons_against_position = starts.windows(2).any(|pair| {
            if reverse {
                pair[1] > pair[0]
            } else {
                pair[1] < pair[0]
            }
        });
        Some(Extra {
            stable_id,
            start: u64::try_from(stored_start).unwrap_or(1).max(1),
            end,
            db_id: self.db_id.and_then(|db_id| db_id.0),
            exons_against_position,
            peptide_first: self
                .feature_cache
                .and_then(|cache| cache.peptide)
                .and_then(|first| first.0),
        })
    }
}

/// A shard's records in either shape the converter writes: an array, or a map from a
/// region name to an array; a null entry of an array is no record.
struct RawShard(Vec<Extra>);

impl<'de> Deserialize<'de> for RawShard {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<RawShard, D::Error> {
        struct Shape;

        impl<'de> Visitor<'de> for Shape {
            type Value = RawShard;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("an array of transcripts or a map of such arrays")
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<RawShard, A::Error> {
                let mut extras = Vec::new();
                while let Some(record) = seq.next_element::<Option<RawRecord>>()? {
                    extras.extend(record.and_then(RawRecord::into_extra));
                }
                Ok(RawShard(extras))
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<RawShard, A::Error> {
                let mut extras = Vec::new();
                while let Some((_, list)) =
                    map.next_entry::<IgnoredAny, Vec<Option<RawRecord>>>()?
                {
                    extras.extend(list.into_iter().flatten().filter_map(RawRecord::into_extra));
                }
                Ok(RawShard(extras))
            }
        }

        deserializer.deserialize_any(Shape)
    }
}

/// Reads the shard at `path`, plain or gzip by its first two bytes, one record in
/// memory at a time. The error is the file's path and what is wrong with it.
pub fn read_shard(path: &Path) -> Result<Vec<Extra>, String> {
    let failed = |what: String| format!("{}: {what}", path.display());
    let mut file = File::open(path).map_err(|error| failed(error.to_string()))?;
    let mut magic = [0u8; 2];
    let mut filled = 0;
    while filled < magic.len() {
        match file.read(&mut magic[filled..]) {
            Ok(0) => break,
            Ok(read) => filled += read,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(failed(error.to_string())),
        }
    }
    let whole = std::io::Cursor::new(magic[..filled].to_vec()).chain(file);
    let parsed: Result<RawShard, serde_json::Error> = if magic == [0x1F, 0x8B] {
        serde_json::from_reader(BufReader::new(MultiGzDecoder::new(whole)))
    } else {
        serde_json::from_reader(BufReader::new(whole))
    };
    parsed
        .map(|RawShard(extras)| extras)
        .map_err(|error| failed(error.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn extras(json: &str) -> Vec<Extra> {
        let RawShard(extras) = serde_json::from_str(json).expect("a shard");
        extras
    }

    #[test]
    fn exons_listed_from_the_five_prime_end_run_with_position_on_either_strand() {
        let read = extras(
            r#"[
              {"stable_id":"F","start":100,"end":400,"strand":1,
               "exons":[{"start":100,"end":150},{"start":300,"end":400}]},
              {"stable_id":"R","start":"100","end":"400","strand":-1,
               "exons":[{"start":300,"end":400},{"start":100,"end":150}]}
            ]"#,
        );
        assert_eq!(read.len(), 2);
        assert!(!read[0].exons_against_position);
        assert!(!read[1].exons_against_position);
        assert_eq!((read[1].start, read[1].end), (100, 400));
    }

    #[test]
    fn an_exon_listed_on_the_wrong_side_of_the_one_before_it_runs_against_position() {
        let read = extras(
            r#"[
              {"stable_id":"F","start":100,"end":400,"strand":1,
               "exons":[{"start":300,"end":400},{"start":100,"end":150}]},
              {"stable_id":"R","start":100,"end":400,"strand":-1,
               "exons":[{"start":100,"end":150},{"start":300,"end":400}]},
              {"stable_id":"THIRD","start":100,"end":900,"strand":-1,
               "exons":[{"start":800,"end":900},{"start":100,"end":150},{"start":300,"end":400}]}
            ]"#,
        );
        assert!(read.iter().all(|extra| extra.exons_against_position));
    }

    #[test]
    fn the_peptide_is_read_down_to_its_first_residue() {
        let read = extras(
            r#"{"2:1-1000000":[
              {"stable_id":"A","start":1,"end":9,"variation_effect_feature_cache":{"peptide":"MK"}},
              null,
              {"stable_id":"B","start":1,"end":9,"variation_effect_feature_cache":{"peptide":""}},
              {"stable_id":"C","start":1,"end":9,"variation_effect_feature_cache":{"peptide":null}},
              {"stable_id":"D","start":1,"end":9}
            ]}"#,
        );
        let firsts: Vec<Option<u8>> = read.iter().map(|extra| extra.peptide_first).collect();
        assert_eq!(firsts, [Some(b'M'), None, None, None]);
    }

    #[test]
    fn the_database_id_is_read_as_a_number_or_a_string_that_spells_one() {
        let read = extras(
            r#"[
              {"stable_id":"A","dbID":41,"start":1,"end":9},
              {"stable_id":"A","dbID":"42","start":1,"end":9},
              {"stable_id":"A","dbID":" 43 ","start":1,"end":9},
              {"stable_id":"A","dbID":-44,"start":1,"end":9},
              {"stable_id":"A","dbID":45.0,"start":1,"end":9},
              {"stable_id":"A","dbID":null,"start":1,"end":9},
              {"stable_id":"A","dbID":"","start":1,"end":9},
              {"stable_id":"A","start":1,"end":9}
            ]"#,
        );
        let ids: Vec<Option<i64>> = read.iter().map(|extra| extra.db_id).collect();
        assert_eq!(
            ids,
            [
                Some(41),
                Some(42),
                Some(43),
                Some(-44),
                Some(45),
                None,
                None,
                None
            ]
        );
    }

    #[test]
    fn a_start_stored_below_one_reads_as_one_and_a_record_without_a_place_is_passed_over() {
        let read = extras(
            r#"[
              {"stable_id":"O","start":-20,"end":300,"strand":1},
              {"stable_id":"X","end":300},
              {"start":1,"end":300},
              {"stable_id":"Y","start":1,"end":-3}
            ]"#,
        );
        assert_eq!(read.len(), 1);
        assert_eq!(
            (read[0].stable_id.as_str(), read[0].start, read[0].end),
            ("O", 1, 300)
        );
    }
}
