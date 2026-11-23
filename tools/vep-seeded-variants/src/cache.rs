// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader of a converted cache, the JSON layout `tools/vep-storable-to-json`
//! writes: `info.json`, `chr_synonyms.txt`, the optional `sequence_regions.tsv` and the
//! shards `transcripts/<sequence>/<start>-<end>.json`, plain or gzip. A transcript is read
//! down to its place, its exons and its coding sequence, one shard in memory at a time
//! beside the stable id, start, end and `dbID` of each of the sequence's transcripts already
//! given.

use std::collections::{BTreeSet, HashSet};
use std::fmt;
use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use serde::de::{self, DeserializeOwned, Deserializer, IgnoredAny, MapAccess, SeqAccess, Visitor};
use serde::Deserialize;
use vep_core::codon::complement_base;

/// The name of the side file that holds each sequence's length and circular flag.
pub const SEQUENCE_REGIONS_FILE: &str = "sequence_regions.tsv";

/// The header line of the side file.
pub const SEQUENCE_REGIONS_HEADER: &str = "name\tlength\tcircular";

/// A file of the cache that cannot be read as the layout states it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CacheError {
    /// The file or directory.
    pub path: PathBuf,
    /// What is wrong with it.
    pub what: String,
}

impl CacheError {
    fn new(path: &Path, what: impl Into<String>) -> CacheError {
        CacheError {
            path: path.to_path_buf(),
            what: what.into(),
        }
    }
}

impl fmt::Display for CacheError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.path.display(), self.what)
    }
}

impl std::error::Error for CacheError {}

/// What `info.json` says of the cache.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct Info {
    /// The species name, as the cache directory spells it.
    #[serde(default)]
    pub species: Option<String>,
    /// The assembly name.
    #[serde(default)]
    pub assembly: Option<String>,
    /// The converted sequences in the converter's order; empty in a cache whose
    /// `info.json` lists none.
    #[serde(default)]
    pub contig_order: Vec<String>,
    /// Every other key of the file.
    #[serde(flatten)]
    pub other: serde_json::Map<String, serde_json::Value>,
}

/// One row of the side file: a sequence, its length in bases and whether it is circular,
/// `None` where the Ensembl cache states neither.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRegion {
    pub name: String,
    pub length: u64,
    pub circular: Option<bool>,
}

/// The strand a transcript is read on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

/// A run of bases on a sequence, both ends included, counted from 1.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Span {
    pub start: u64,
    pub end: u64,
}

impl Span {
    /// How many bases the run holds.
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }

    /// Whether `position` lies in the run.
    pub fn contains(&self, position: u64) -> bool {
        self.start <= position && position <= self.end
    }
}

/// The coding model of a transcript, on the transcript's own strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Coding {
    /// The coding sequence's first base on the spliced transcript, counted from 1.
    pub cdna_start: u64,
    /// Its last base there.
    pub cdna_end: u64,
    /// The mapper's start phase: above zero, the cached coding sequence opens with that
    /// many padding bases.
    pub start_phase: i8,
    /// The cached coding sequence (`translateable_seq`), upper-cased.
    pub sequence: Vec<u8>,
    /// The 5' untranslated region's bases, where the cache carries them.
    pub five_prime_utr: Option<Vec<u8>>,
    /// The 3' untranslated region's bases, where the cache carries them.
    pub three_prime_utr: Option<Vec<u8>>,
}

/// Where a transcript stored across the origin of a circular sequence begins, among the
/// sequence's last bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BeforeOrigin {
    /// How many of the transcript's bases lie before the origin: 1 less its stored start.
    pub bases: u64,
    /// The transcript's first base as its mapper pairs give it: the lowest pair start
    /// above the transcript's end, `None` where no pair starts there. The pairs keep the
    /// sequence's own positions, the pair across the origin with its start above its end
    /// (`TranscriptMapper.pm:160-161` of Ensembl's core API).
    pub pair_start: Option<u64>,
}

impl BeforeOrigin {
    /// The position of the transcript's first base on a sequence of `length` bases, where
    /// the side file gives the length: the lower of what the mapper pairs and the length
    /// say, `None` where neither says.
    pub fn first_base(&self, length: Option<u64>) -> Option<u64> {
        let by_length = length.and_then(|length| length.checked_add(1)?.checked_sub(self.bases));
        match (self.pair_start, by_length) {
            (Some(by_pairs), Some(by_length)) => Some(by_pairs.min(by_length)),
            (by_pairs, by_length) => by_pairs.or(by_length),
        }
    }
}

/// A transcript, down to the fields a record's place and alleles are read from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub stable_id: String,
    pub gene_stable_id: Option<String>,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub biotype: String,
    /// The exons from the transcript's 5' end to its 3' end; none for a transcript one of
    /// whose exons the record places outside the transcript's own span.
    pub exons: Vec<Span>,
    /// The genetic code the record names (`codon_table`), `None` where it names none.
    pub codon_table: Option<u8>,
    /// The coding model, `None` for a transcript without a coding sequence.
    pub coding: Option<Coding>,
    /// The `cds_start_NF` attribute is present.
    pub cds_start_nf: bool,
    /// The `cds_end_NF` attribute is present.
    pub cds_end_nf: bool,
    /// The transcript carries a sequence edit, an `_rna_edit` attribute or a BAM edit
    /// status, so its bases are not the genome's at every position.
    pub edited: bool,
    /// The record's start is below 1: the transcript lies across the origin of a circular
    /// sequence, stored with the sequence's length taken from its start. `start` is then
    /// 1 and no exon is kept.
    pub across_origin: bool,
    /// Where such a transcript begins; `None` for every other transcript.
    pub before_origin: Option<BeforeOrigin>,
}

impl Transcript {
    /// The genetic code the transcript is translated under: the one it names, and the
    /// standard code where it names none or names 0.
    pub fn code(&self) -> u8 {
        match self.codon_table {
            None | Some(0) => 1,
            Some(table) => table,
        }
    }

    /// The length of the spliced transcript.
    pub fn cdna_length(&self) -> u64 {
        self.exons.iter().map(Span::length).sum()
    }

    /// Each exon's first and last base on the spliced transcript, counted from 1, in the
    /// exons' order.
    pub fn exon_cdna_spans(&self) -> Vec<Span> {
        let mut before = 0;
        self.exons
            .iter()
            .map(|exon| {
                let span = Span {
                    start: before + 1,
                    end: before + exon.length(),
                };
                before = span.end;
                span
            })
            .collect()
    }

    /// The sequence position of base `cdna` of the spliced transcript, counted from 1.
    pub fn genomic_position(&self, cdna: u64) -> Option<u64> {
        let mut before = 0;
        for exon in &self.exons {
            let length = exon.length();
            if cdna > before && cdna <= before + length {
                let offset = cdna - before - 1;
                return Some(match self.strand {
                    Strand::Forward => exon.start + offset,
                    Strand::Reverse => exon.end - offset,
                });
            }
            before += length;
        }
        None
    }

    /// The base of the spliced transcript, counted from 1, that sequence position
    /// `genomic` holds; `None` outside every exon.
    pub fn cdna_position(&self, genomic: u64) -> Option<u64> {
        let mut before = 0;
        for exon in &self.exons {
            if exon.contains(genomic) {
                let offset = match self.strand {
                    Strand::Forward => genomic - exon.start,
                    Strand::Reverse => exon.end - genomic,
                };
                return Some(before + offset + 1);
            }
            before += exon.length();
        }
        None
    }

    /// The introns from the 5' end to the 3' end: the bases between two exons that
    /// follow each other.
    pub fn introns(&self) -> Vec<Span> {
        self.exons
            .windows(2)
            .filter_map(|pair| {
                let (low, high) = match self.strand {
                    Strand::Forward => (pair[0], pair[1]),
                    Strand::Reverse => (pair[1], pair[0]),
                };
                (low.end + 1 < high.start).then(|| Span {
                    start: low.end + 1,
                    end: high.start - 1,
                })
            })
            .collect()
    }

    /// Whether the cache gives this transcript's bases as the genome's and its codons
    /// from the first base of a complete coding sequence: a coding model without a start
    /// phase above zero, without either incomplete-coding-sequence flag and without an
    /// edit, whose coding sequence is as long as its span on the spliced transcript and
    /// holds the four bases alone, and whose untranslated regions, where carried, fill
    /// the rest of the exons.
    pub fn reads_codons(&self) -> bool {
        let Some(coding) = &self.coding else {
            return false;
        };
        if self.edited || self.cds_start_nf || self.cds_end_nf || coding.start_phase > 0 {
            return false;
        }
        let length = self.cdna_length();
        if coding.cdna_start < 1 || coding.cdna_end < coding.cdna_start || coding.cdna_end > length
        {
            return false;
        }
        if coding.sequence.len() as u64 != coding.cdna_end - coding.cdna_start + 1
            || !coding.sequence.iter().all(|base| b"ACGT".contains(base))
        {
            return false;
        }
        let fills = |utr: &Option<Vec<u8>>, bases: u64| {
            utr.as_ref().is_none_or(|utr| utr.len() as u64 == bases)
        };
        fills(&coding.five_prime_utr, coding.cdna_start - 1)
            && fills(&coding.three_prime_utr, length - coding.cdna_end)
            && self.exons_are_apart()
    }

    /// Whether the exons follow each other along the strand without touching.
    fn exons_are_apart(&self) -> bool {
        !self.exons.is_empty()
            && self.exons.windows(2).all(|pair| match self.strand {
                Strand::Forward => pair[0].end < pair[1].start,
                Strand::Reverse => pair[1].end < pair[0].start,
            })
    }

    /// The base of the spliced transcript at `cdna`, on the transcript's strand, where
    /// [`Transcript::reads_codons`] holds and the cache carries that stretch: the coding
    /// sequence always, an untranslated region where its string is present.
    pub fn cdna_base(&self, cdna: u64) -> Option<u8> {
        if !self.reads_codons() {
            return None;
        }
        let coding = self.coding.as_ref()?;
        let base = if cdna < coding.cdna_start {
            *coding
                .five_prime_utr
                .as_ref()?
                .get(usize::try_from(cdna.checked_sub(1)?).ok()?)?
        } else if cdna <= coding.cdna_end {
            *coding
                .sequence
                .get(usize::try_from(cdna - coding.cdna_start).ok()?)?
        } else {
            *coding
                .three_prime_utr
                .as_ref()?
                .get(usize::try_from(cdna - coding.cdna_end - 1).ok()?)?
        };
        b"ACGT".contains(&base).then_some(base)
    }

    /// The forward-strand base at sequence position `genomic`, under the conditions of
    /// [`Transcript::cdna_base`].
    pub fn forward_base(&self, genomic: u64) -> Option<u8> {
        let base = self.cdna_base(self.cdna_position(genomic)?)?;
        Some(match self.strand {
            Strand::Forward => base,
            Strand::Reverse => complement_base(base),
        })
    }
}

/// A JSON number, a string that spells one, or null.
struct Number;

impl<'de> Visitor<'de> for Number {
    type Value = Option<i64>;

    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a whole number, a string that spells one, or null")
    }

    fn visit_i64<E: de::Error>(self, value: i64) -> Result<Self::Value, E> {
        Ok(Some(value))
    }

    fn visit_u64<E: de::Error>(self, value: u64) -> Result<Self::Value, E> {
        i64::try_from(value)
            .map(Some)
            .map_err(|_| E::invalid_value(de::Unexpected::Unsigned(value), &self))
    }

    fn visit_f64<E: de::Error>(self, value: f64) -> Result<Self::Value, E> {
        if value.fract() == 0.0 && value.abs() < 9.0e15 {
            Ok(Some(value as i64))
        } else {
            Err(E::invalid_value(de::Unexpected::Float(value), &self))
        }
    }

    fn visit_str<E: de::Error>(self, value: &str) -> Result<Self::Value, E> {
        let text = value.trim();
        if text.is_empty() {
            return Ok(None);
        }
        text.parse::<i64>()
            .map(Some)
            .map_err(|_| E::invalid_value(de::Unexpected::Str(value), &self))
    }

    fn visit_none<E: de::Error>(self) -> Result<Self::Value, E> {
        Ok(None)
    }

    fn visit_unit<E: de::Error>(self) -> Result<Self::Value, E> {
        Ok(None)
    }

    fn visit_some<D: Deserializer<'de>>(self, deserializer: D) -> Result<Self::Value, D::Error> {
        deserializer.deserialize_any(Number)
    }
}

fn number<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Option<i64>, D::Error> {
    deserializer.deserialize_any(Number)
}

#[derive(Deserialize)]
struct RawSpan {
    #[serde(default, deserialize_with = "number")]
    start: Option<i64>,
    #[serde(default, deserialize_with = "number")]
    end: Option<i64>,
}

#[derive(Deserialize)]
struct RawAttribute {
    #[serde(default)]
    code: Option<String>,
}

/// A mapper pair, down to the sequence position its run of bases starts at.
#[derive(Deserialize)]
struct RawPair {
    #[serde(default, deserialize_with = "number")]
    to_start: Option<i64>,
}

#[derive(Deserialize)]
struct RawMapper {
    #[serde(default, deserialize_with = "number")]
    start_phase: Option<i64>,
    #[serde(default, deserialize_with = "number")]
    cdna_coding_start: Option<i64>,
    #[serde(default, deserialize_with = "number")]
    cdna_coding_end: Option<i64>,
    #[serde(default)]
    pairs: Option<Vec<Option<RawPair>>>,
}

#[derive(Deserialize)]
struct RawFeatureCache {
    #[serde(default, deserialize_with = "number")]
    codon_table: Option<i64>,
    #[serde(default)]
    five_prime_utr: Option<String>,
    #[serde(default)]
    three_prime_utr: Option<String>,
    #[serde(default)]
    translateable_seq: Option<String>,
    #[serde(default)]
    seq_edits: Option<Vec<IgnoredAny>>,
    #[serde(default)]
    mapper: Option<RawMapper>,
}

/// A transcript record of a shard: the keys this reader keeps, every other key skipped.
#[derive(Deserialize)]
struct RawTranscript {
    stable_id: String,
    #[serde(default, rename = "dbID", deserialize_with = "number")]
    db_id: Option<i64>,
    #[serde(default)]
    gene_stable_id: Option<String>,
    #[serde(default, deserialize_with = "number")]
    start: Option<i64>,
    #[serde(default, deserialize_with = "number")]
    end: Option<i64>,
    #[serde(default, deserialize_with = "number")]
    strand: Option<i64>,
    #[serde(default)]
    biotype: Option<String>,
    #[serde(default)]
    exons: Option<Vec<RawSpan>>,
    #[serde(default, deserialize_with = "number")]
    cdna_coding_start: Option<i64>,
    #[serde(default, deserialize_with = "number")]
    cdna_coding_end: Option<i64>,
    #[serde(default)]
    attributes: Option<Vec<RawAttribute>>,
    #[serde(default)]
    bam_edit_status: Option<IgnoredAny>,
    #[serde(default, rename = "variation_effect_feature_cache")]
    feature_cache: Option<RawFeatureCache>,
}

/// A shard's records in either shape: an array, or a map from a region name to an array;
/// a null entry of an array is no record.
struct RawShard(Vec<RawTranscript>);

impl<'de> Deserialize<'de> for RawShard {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<RawShard, D::Error> {
        struct Shape;

        impl<'de> Visitor<'de> for Shape {
            type Value = RawShard;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("an array of transcripts or a map of such arrays")
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<RawShard, A::Error> {
                let mut records = Vec::new();
                while let Some(record) = seq.next_element::<Option<RawTranscript>>()? {
                    records.extend(record);
                }
                Ok(RawShard(records))
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<RawShard, A::Error> {
                let mut records = Vec::new();
                while let Some((_, list)) =
                    map.next_entry::<IgnoredAny, Vec<Option<RawTranscript>>>()?
                {
                    records.extend(list.into_iter().flatten());
                }
                Ok(RawShard(records))
            }
        }

        deserializer.deserialize_any(Shape)
    }
}

fn upper(text: String) -> Vec<u8> {
    let mut bytes = text.into_bytes();
    bytes.make_ascii_uppercase();
    bytes
}

fn positive(value: Option<i64>) -> Option<u64> {
    value.and_then(|value| u64::try_from(value).ok())
}

impl RawTranscript {
    fn into_transcript(self) -> Result<Transcript, String> {
        let id = self.stable_id;
        let (Some(stored_start), Some(end)) = (self.start, positive(self.end)) else {
            return Err(format!("transcript {id} has no start or no end"));
        };
        if end < 1 || i64::try_from(end).is_ok_and(|end| stored_start > end) {
            return Err(format!("transcript {id} runs from {stored_start} to {end}"));
        }
        let across_origin = stored_start < 1;
        let start = u64::try_from(stored_start).unwrap_or(1).max(1);
        let strand = if self.strand.unwrap_or(1) >= 0 {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        let mut exons = Vec::new();
        let mut placed = !across_origin;
        for exon in self.exons.unwrap_or_default() {
            match (positive(exon.start), positive(exon.end)) {
                (Some(first), Some(last)) if start <= first && first <= last && last <= end => {
                    exons.push(Span {
                        start: first,
                        end: last,
                    });
                }
                _ => placed = false,
            }
        }
        if !placed {
            exons.clear();
        }
        exons.sort();
        if strand == Strand::Reverse {
            exons.reverse();
        }
        let attributes = self.attributes.unwrap_or_default();
        let has = |code: &str| {
            attributes
                .iter()
                .any(|attribute| attribute.code.as_deref() == Some(code))
        };
        let mut edited = has("_rna_edit") || self.bam_edit_status.is_some();
        let mut codon_table = None;
        let mut coding = None;
        let before_origin = across_origin.then(|| BeforeOrigin {
            bases: stored_start.unsigned_abs().saturating_add(1),
            pair_start: self
                .feature_cache
                .as_ref()
                .and_then(|cache| cache.mapper.as_ref())
                .and_then(|mapper| mapper.pairs.as_ref())
                .into_iter()
                .flatten()
                .flatten()
                .filter_map(|pair| positive(pair.to_start))
                .filter(|pair_start| *pair_start > end)
                .min(),
        });
        if let Some(cache) = self.feature_cache {
            if let Some(table) = cache.codon_table {
                codon_table = Some(
                    u8::try_from(table)
                        .map_err(|_| format!("transcript {id} names the genetic code {table}"))?,
                );
            }
            edited |= cache.seq_edits.is_some_and(|edits| !edits.is_empty());
            let mapper = cache.mapper.as_ref();
            let cdna_start = positive(self.cdna_coding_start)
                .or_else(|| positive(mapper.and_then(|mapper| mapper.cdna_coding_start)));
            let cdna_end = positive(self.cdna_coding_end)
                .or_else(|| positive(mapper.and_then(|mapper| mapper.cdna_coding_end)));
            let sequence = cache.translateable_seq.filter(|text| !text.is_empty());
            if let (Some(sequence), Some(cdna_start), Some(cdna_end)) =
                (sequence, cdna_start, cdna_end)
            {
                let start_phase = mapper
                    .and_then(|mapper| mapper.start_phase)
                    .map_or(-1, |phase| phase.clamp(-1, 2) as i8);
                coding = Some(Coding {
                    cdna_start,
                    cdna_end,
                    start_phase,
                    sequence: upper(sequence),
                    five_prime_utr: cache.five_prime_utr.map(upper),
                    three_prime_utr: cache.three_prime_utr.map(upper),
                });
            }
        }
        Ok(Transcript {
            stable_id: id,
            gene_stable_id: self.gene_stable_id,
            start,
            end,
            strand,
            biotype: self.biotype.unwrap_or_default(),
            exons,
            codon_table,
            coding,
            cds_start_nf: has("cds_start_NF"),
            cds_end_nf: has("cds_end_NF"),
            edited,
            across_origin,
            before_origin,
        })
    }
}

/// A shard is `<start>-<end>.json`, plain or with a `.gz` ending.
fn is_shard_name(name: &str) -> bool {
    name.ends_with(".json") || name.ends_with(".json.gz")
}

/// The order shards are read in: by the start their name states, a name that states none
/// last, equal starts by name.
fn shard_order(path: &Path) -> (u64, String) {
    let name = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let start = name
        .split('-')
        .next()
        .and_then(|start| start.parse().ok())
        .unwrap_or(u64::MAX);
    (start, name)
}

/// Reads a JSON file, plain or gzip by its first two bytes, without holding its text.
fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T, CacheError> {
    let mut file = File::open(path).map_err(|error| CacheError::new(path, error.to_string()))?;
    let mut magic = [0u8; 2];
    let read = file
        .read(&mut magic)
        .map_err(|error| CacheError::new(path, error.to_string()))?;
    let head = std::io::Cursor::new(magic[..read].to_vec());
    let whole = head.chain(file);
    let parsed = if magic == [0x1F, 0x8B] {
        serde_json::from_reader(BufReader::new(MultiGzDecoder::new(whole)))
    } else {
        serde_json::from_reader(BufReader::new(whole))
    };
    parsed.map_err(|error| CacheError::new(path, error.to_string()))
}

/// Reads the side file: the header line, then one row a sequence.
pub fn read_sequence_regions(path: &Path) -> Result<Vec<SequenceRegion>, CacheError> {
    let text =
        fs::read_to_string(path).map_err(|error| CacheError::new(path, error.to_string()))?;
    let mut lines = text.lines();
    if lines.next() != Some(SEQUENCE_REGIONS_HEADER) {
        return Err(CacheError::new(
            path,
            "the first line is not the header name, length, circular",
        ));
    }
    let mut regions: Vec<SequenceRegion> = Vec::new();
    for (index, line) in lines.enumerate() {
        let row = index + 2;
        let fields: Vec<&str> = line.split('\t').collect();
        let [name, length, circular] = fields.as_slice() else {
            return Err(CacheError::new(
                path,
                format!("line {row} does not have three fields"),
            ));
        };
        let length: u64 = match length.parse() {
            Ok(length) if length >= 1 => length,
            _ => {
                return Err(CacheError::new(
                    path,
                    format!("line {row}: the length {length} is not a whole number above zero"),
                ))
            }
        };
        let circular = match *circular {
            "1" => Some(true),
            "0" => Some(false),
            "-" => None,
            other => {
                return Err(CacheError::new(
                    path,
                    format!("line {row}: the circular flag {other} is not 1, 0 or -"),
                ))
            }
        };
        if name.is_empty() || regions.iter().any(|region| region.name == *name) {
            return Err(CacheError::new(
                path,
                format!("line {row}: the name is empty or stands twice"),
            ));
        }
        regions.push(SequenceRegion {
            name: (*name).to_owned(),
            length,
            circular,
        });
    }
    Ok(regions)
}

/// A converted cache directory.
#[derive(Debug, Clone)]
pub struct Cache {
    dir: PathBuf,
    info: Info,
    synonyms: Vec<(String, String)>,
    regions: Option<Vec<SequenceRegion>>,
    sequences: Vec<String>,
}

impl Cache {
    /// Opens the cache at `dir`, with the side file beside `info.json` where it exists.
    pub fn open(dir: &Path) -> Result<Cache, CacheError> {
        let side_file = dir.join(SEQUENCE_REGIONS_FILE);
        let side_file = side_file.is_file().then_some(side_file);
        Cache::open_with_regions(dir, side_file.as_deref())
    }

    /// Opens the cache at `dir` with the side file at `regions`, or with none.
    pub fn open_with_regions(dir: &Path, regions: Option<&Path>) -> Result<Cache, CacheError> {
        let info: Info = read_json(&dir.join("info.json"))?;
        let transcripts = dir.join("transcripts");
        if !transcripts.is_dir() {
            return Err(CacheError::new(&transcripts, "no transcripts directory"));
        }
        let mut directories: BTreeSet<String> = BTreeSet::new();
        let entries = fs::read_dir(&transcripts)
            .map_err(|error| CacheError::new(&transcripts, error.to_string()))?;
        for entry in entries {
            let entry = entry.map_err(|error| CacheError::new(&transcripts, error.to_string()))?;
            let name = entry.file_name().to_string_lossy().into_owned();
            if !name.starts_with('.') && entry.path().is_dir() {
                directories.insert(name);
            }
        }
        let mut sequences = Vec::with_capacity(directories.len());
        for name in &info.contig_order {
            if directories.remove(name) {
                sequences.push(name.clone());
            }
        }
        sequences.extend(directories);
        let synonyms_file = dir.join("chr_synonyms.txt");
        let synonyms = if synonyms_file.is_file() {
            let text = fs::read_to_string(&synonyms_file)
                .map_err(|error| CacheError::new(&synonyms_file, error.to_string()))?;
            synonym_pairs(&text)
        } else {
            Vec::new()
        };
        let regions = regions.map(read_sequence_regions).transpose()?;
        Ok(Cache {
            dir: dir.to_path_buf(),
            info,
            synonyms,
            regions,
            sequences,
        })
    }

    /// The directory the cache was opened at.
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    /// What `info.json` says.
    pub fn info(&self) -> &Info {
        &self.info
    }

    /// The sequences in reading order: the names of `contig_order` that have a directory
    /// under `transcripts`, each once, then every directory it does not list, in bytewise
    /// name order. Ensembl VEP and vep-rs know a sequence by its directory, so a name
    /// without one is no sequence of the cache.
    pub fn sequences(&self) -> &[String] {
        &self.sequences
    }

    /// The pairs of `chr_synonyms.txt` in file order: a line's first name with each
    /// later name of the line. Empty for a cache without the file.
    pub fn synonym_pairs(&self) -> &[(String, String)] {
        &self.synonyms
    }

    /// The rows of the side file, `None` for a cache opened without one.
    pub fn sequence_regions(&self) -> Option<&[SequenceRegion]> {
        self.regions.as_deref()
    }

    /// The side file's row for `sequence`.
    pub fn sequence_region(&self, sequence: &str) -> Option<&SequenceRegion> {
        self.regions
            .as_ref()?
            .iter()
            .find(|region| region.name == sequence)
    }

    /// The shards of `sequence` in reading order, whatever order the directory lists its
    /// files in; none for a name without a directory.
    pub fn shards(&self, sequence: &str) -> Result<Vec<PathBuf>, CacheError> {
        let dir = self.dir.join("transcripts").join(sequence);
        if !dir.is_dir() {
            return Ok(Vec::new());
        }
        let mut shards = Vec::new();
        let entries =
            fs::read_dir(&dir).map_err(|error| CacheError::new(&dir, error.to_string()))?;
        for entry in entries {
            let entry = entry.map_err(|error| CacheError::new(&dir, error.to_string()))?;
            if is_shard_name(&entry.file_name().to_string_lossy()) {
                shards.push(entry.path());
            }
        }
        shards.sort_by_key(|path| shard_order(path));
        Ok(shards)
    }

    /// The transcripts of `sequence`: shard after shard in reading order, a shard's
    /// transcripts by start, then stable id, then end, then `dbID`. A transcript is
    /// given once. Ensembl VEP writes into a region file every transcript of each gene
    /// that overlaps the region (`AnnotationSource/Database/Transcript.pm:237` and
    /// `:250`), so a shard repeats transcripts of earlier shards, wherever they end: a
    /// record under the stable id, start, end and `dbID` of a transcript already given
    /// is that transcript again. Ensembl VEP tells the transcripts it loads apart by
    /// `dbID` (`AnnotationType/Transcript.pm:283-297`): a record under the stable id of
    /// a transcript already given and another `dbID`, at the same place or at another,
    /// is another transcript. A record without a `dbID` is told apart by its stable
    /// id, start and end.
    pub fn transcripts(&self, sequence: &str) -> Result<Transcripts, CacheError> {
        Ok(Transcripts {
            shards: self.shards(sequence)?.into_iter(),
            current: Vec::new().into_iter(),
            given: HashSet::new(),
            failed: false,
        })
    }
}

/// The pairs a synonyms file states: each later name of a line with the line's first.
fn synonym_pairs(text: &str) -> Vec<(String, String)> {
    let mut pairs = Vec::new();
    for line in text.lines() {
        let mut names = line
            .split(|c: char| c.is_ascii_whitespace() || c == '\u{0B}')
            .filter(|name| !name.is_empty());
        let Some(first) = names.next() else {
            continue;
        };
        // A line that opens with white space has an empty first name in Ensembl VEP's
        // reading and links nothing to a sequence.
        if line.starts_with(|c: char| c.is_ascii_whitespace() || c == '\u{0B}') {
            continue;
        }
        for name in names {
            pairs.push((first.to_owned(), name.to_owned()));
        }
    }
    pairs
}

/// The iterator of [`Cache::transcripts`].
#[derive(Debug)]
pub struct Transcripts {
    shards: std::vec::IntoIter<PathBuf>,
    current: std::vec::IntoIter<Transcript>,
    /// The stable id, start, end and `dbID` of each transcript already given.
    given: HashSet<(String, u64, u64, Option<i64>)>,
    failed: bool,
}

impl Transcripts {
    fn load(&mut self, path: &Path) -> Result<(), CacheError> {
        let RawShard(records) = read_json(path)?;
        let mut transcripts = Vec::with_capacity(records.len());
        for record in records {
            let db_id = record.db_id;
            let transcript = record
                .into_transcript()
                .map_err(|what| CacheError::new(path, what))?;
            transcripts.push((transcript, db_id));
        }
        transcripts.sort_by(|(a, a_id), (b, b_id)| {
            (a.start, &a.stable_id, a.end, a_id).cmp(&(b.start, &b.stable_id, b.end, b_id))
        });
        transcripts.retain(|(transcript, db_id)| {
            self.given.insert((
                transcript.stable_id.clone(),
                transcript.start,
                transcript.end,
                *db_id,
            ))
        });
        let transcripts: Vec<Transcript> = transcripts
            .into_iter()
            .map(|(transcript, _)| transcript)
            .collect();
        self.current = transcripts.into_iter();
        Ok(())
    }
}

impl Iterator for Transcripts {
    type Item = Result<Transcript, CacheError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.failed {
                return None;
            }
            if let Some(transcript) = self.current.next() {
                return Some(Ok(transcript));
            }
            let shard = self.shards.next()?;
            if let Err(error) = self.load(&shard) {
                self.failed = true;
                return Some(Err(error));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn transcript(json: &str) -> Transcript {
        serde_json::from_str::<RawTranscript>(json)
            .expect("a transcript record")
            .into_transcript()
            .expect("a transcript")
    }

    const FORWARD: &str = r#"{"stable_id":"T1","start":"101","end":230,"strand":1,
        "biotype":"protein_coding","cdna_coding_start":4,"cdna_coding_end":"15",
        "exons":[{"start":101,"end":110},{"start":221,"end":230}],
        "variation_effect_feature_cache":{"codon_table":"11","five_prime_utr":"ggg",
          "translateable_seq":"ATGAAACCCTAA","three_prime_utr":"TTTTT",
          "mapper":{"start_phase":-1}}}"#;

    const REVERSE: &str = r#"{"stable_id":"T2","start":101,"end":230,"strand":-1,
        "biotype":"protein_coding","cdna_coding_start":1,"cdna_coding_end":12,
        "exons":[{"start":221,"end":230},{"start":101,"end":110}],
        "variation_effect_feature_cache":{"translateable_seq":"ATGAAACCCTAA",
          "three_prime_utr":"TTTTTTTT","mapper":{"start_phase":0}}}"#;

    #[test]
    fn a_record_is_read_down_to_its_place_and_coding_model() {
        let forward = transcript(FORWARD);
        assert_eq!((forward.start, forward.end), (101, 230));
        assert_eq!(forward.strand, Strand::Forward);
        assert_eq!(forward.code(), 11);
        assert_eq!(forward.cdna_length(), 20);
        let coding = forward.coding.as_ref().expect("a coding model");
        assert_eq!((coding.cdna_start, coding.cdna_end), (4, 15));
        assert_eq!(coding.five_prime_utr.as_deref(), Some(&b"GGG"[..]));
        assert!(forward.reads_codons());
        assert_eq!(
            forward.introns(),
            [Span {
                start: 111,
                end: 220
            }]
        );
        let reverse = transcript(REVERSE);
        assert_eq!(reverse.code(), 1);
        assert_eq!(reverse.exons[0].start, 221);
        assert_eq!(reverse.introns(), forward.introns());
    }

    #[test]
    fn positions_map_between_the_spliced_transcript_and_the_sequence() {
        let forward = transcript(FORWARD);
        assert_eq!(forward.genomic_position(1), Some(101));
        assert_eq!(forward.genomic_position(10), Some(110));
        assert_eq!(forward.genomic_position(11), Some(221));
        assert_eq!(forward.genomic_position(21), None);
        assert_eq!(forward.cdna_position(221), Some(11));
        assert_eq!(forward.cdna_position(150), None);
        let reverse = transcript(REVERSE);
        assert_eq!(reverse.genomic_position(1), Some(230));
        assert_eq!(reverse.genomic_position(11), Some(110));
        assert_eq!(reverse.cdna_position(101), Some(20));
        for cdna in 1..=20 {
            let genomic = reverse.genomic_position(cdna).expect("an exon base");
            assert_eq!(reverse.cdna_position(genomic), Some(cdna));
        }
    }

    #[test]
    fn a_base_is_given_on_the_forward_strand() {
        let forward = transcript(FORWARD);
        assert_eq!(forward.forward_base(101), Some(b'G'));
        assert_eq!(forward.forward_base(104), Some(b'A'));
        assert_eq!(forward.forward_base(226), Some(b'T'));
        assert_eq!(forward.forward_base(150), None);
        let reverse = transcript(REVERSE);
        assert_eq!(reverse.cdna_base(1), Some(b'A'));
        assert_eq!(reverse.forward_base(230), Some(b'T'));
        assert_eq!(reverse.forward_base(101), Some(b'A'));
    }

    #[test]
    fn a_transcript_whose_bases_are_not_the_genomes_reads_no_codon() {
        let padded = FORWARD.replace(r#""start_phase":-1"#, r#""start_phase":2"#);
        assert!(!transcript(&padded).reads_codons());
        let short = FORWARD.replace("ATGAAACCCTAA", "ATGAAACCC");
        assert!(!transcript(&short).reads_codons());
        let flagged = FORWARD.replace(
            r#""biotype""#,
            r#""attributes":[{"code":"cds_start_NF","value":"1"}],"biotype""#,
        );
        assert!(!transcript(&flagged).reads_codons());
        let edited = FORWARD.replace(
            r#""codon_table""#,
            r#""seq_edits":[{"start":2,"end":2,"alt_seq":"U"}],"codon_table""#,
        );
        assert!(transcript(&edited).edited);
        assert_eq!(transcript(&edited).forward_base(104), None);
        let noncoding = r#"{"stable_id":"T3","start":5,"end":9,"strand":1,
            "exons":[{"start":5,"end":9}],"variation_effect_feature_cache":{}}"#;
        assert!(transcript(noncoding).coding.is_none());
    }

    #[test]
    fn a_transcript_stored_across_the_origin_keeps_its_end_and_no_exon() {
        let across = FORWARD
            .replace(r#""start":"101""#, r#""start":-199"#)
            .replace(r#"{"start":101,"end":110}"#, r#"{"start":-199,"end":110}"#);
        let across = transcript(&across);
        assert!(across.across_origin);
        assert_eq!((across.start, across.end), (1, 230));
        assert!(across.exons.is_empty());
        assert!(!across.reads_codons());
        assert!(!transcript(FORWARD).across_origin);
        let outside = FORWARD.replace(r#"{"start":221,"end":230}"#, r#"{"start":221,"end":260}"#);
        assert!(transcript(&outside).exons.is_empty());
    }

    #[test]
    fn a_transcript_stored_across_the_origin_says_where_it_begins() {
        // A sequence of 1,000 bases: the transcript begins at base 801, 200 bases before
        // the origin, so its stored start is 801 - 1,000 = -199.
        let stored = FORWARD.replace(r#""start":"101""#, r#""start":-199"#);
        let pairs = |pairs: &str| {
            let record = stored.replace(
                r#""mapper":{"start_phase":-1}"#,
                &format!(r#""mapper":{{"start_phase":-1,"pairs":{pairs}}}"#),
            );
            assert_ne!(record, stored);
            transcript(&record).before_origin.expect("where it begins")
        };
        // The pair across the origin starts above its end; a pair that starts at or
        // below the transcript's end, 230, lies among the first bases.
        let across = pairs(
            r#"[{"from_start":1,"from_end":310,"to_start":801,"to_end":110},null,
                {"from_start":311,"from_end":320,"to_start":"221","to_end":230}]"#,
        );
        assert_eq!((across.bases, across.pair_start), (200, Some(801)));
        assert_eq!(across.first_base(None), Some(801));
        assert_eq!(across.first_base(Some(1000)), Some(801));
        // An exon wholly before the origin: the lowest start above the end is the first.
        let two = pairs(
            r#"[{"to_start":950,"to_end":990},{"to_start":801,"to_end":900},
                {"to_start":230,"to_end":230},{"to_start":-199,"to_end":110}]"#,
        );
        assert_eq!(two.pair_start, Some(801));
        // No pair above the end: the side file's length says where it begins, or nothing.
        for none in [
            r#"[{"to_start":-199,"to_end":230}]"#,
            r#"[{"to_end":230}]"#,
            "[]",
            "null",
        ] {
            let unplaced = pairs(none);
            assert_eq!((unplaced.bases, unplaced.pair_start), (200, None), "{none}");
            assert_eq!(unplaced.first_base(Some(1000)), Some(801));
            assert_eq!(unplaced.first_base(None), None);
        }
        let absent = transcript(&stored).before_origin.expect("where it begins");
        assert_eq!((absent.bases, absent.pair_start), (200, None));
        // The two readings differ: the lower one stands. A length below the bases
        // stored before the origin places nothing.
        assert_eq!(across.first_base(Some(900)), Some(701));
        assert_eq!(across.first_base(Some(2000)), Some(801));
        assert_eq!(unplaced_before(i64::MIN).first_base(Some(1000)), None);
        assert_eq!(unplaced_before(0).bases, 1);
        assert_eq!(unplaced_before(0).first_base(Some(1000)), Some(1000));
        assert!(transcript(FORWARD).before_origin.is_none());
    }

    /// Where a transcript with this stored start and no mapper pair begins.
    fn unplaced_before(stored_start: i64) -> BeforeOrigin {
        let record = FORWARD.replace(r#""start":"101""#, &format!(r#""start":{stored_start}"#));
        transcript(&record).before_origin.expect("where it begins")
    }

    #[test]
    fn both_shard_shapes_are_read_and_a_null_entry_is_no_record() {
        let array = format!("[{FORWARD},null,{REVERSE}]");
        let RawShard(records) = serde_json::from_str(&array).expect("an array shard");
        assert_eq!(records.len(), 2);
        let map = format!(r#"{{"1":[{FORWARD},null],"2":[{REVERSE}]}}"#);
        let RawShard(records) = serde_json::from_str(&map).expect("a map shard");
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn shard_names_state_their_order() {
        assert!(shard_order(Path::new("2-3.json")) < shard_order(Path::new("10-20.json")));
        assert!(shard_order(Path::new("10-20.json")) < shard_order(Path::new("all.json")));
        assert!(is_shard_name("1-2.json.gz") && !is_shard_name("1-2.json.tmp"));
    }

    #[test]
    fn synonym_pairs_keep_file_order_and_skip_one_name_lines() {
        let pairs = synonym_pairs("1\tchr1\tCM1\nalone\n\nNC_2 2\n");
        let expected: Vec<(String, String)> = [("1", "chr1"), ("1", "CM1"), ("NC_2", "2")]
            .iter()
            .map(|(a, b)| ((*a).to_owned(), (*b).to_owned()))
            .collect();
        assert_eq!(pairs, expected);
    }
}
