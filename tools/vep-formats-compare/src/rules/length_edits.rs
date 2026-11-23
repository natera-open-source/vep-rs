// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The edits of a candidate JSON cache's transcripts that change the length of the cDNA, and
//! the HGVSc position class that reads them.

use crate::argparse::path_text;
use crate::cache;
use crate::classes::USE_TRANSCRIPT_REF_CACHE_CDNA;
use crate::consts::CDNA_END_FIELD;
use crate::consts::CDNA_START_FIELDS;
use crate::consts::HGVS_OFFSET_FIELDS;
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::pyfmt;
use crate::pyjson::Value;
use crate::rules::intron_ends::VERSION_SUFFIX_RE;
use crate::rules::term_pairs::cds_bounds;
use crate::rules::term_pairs::transcript_row;
use crate::rules::transcript_reference::failed_edit_reference_cell;
use crate::rules::transcript_reference::failed_edit_row;
use crate::shared::blank;
use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::HashMap;
use std::path::Path;
use std::path::PathBuf;
use std::sync::LazyLock;

/// A transcript notation in coding numbering whose coordinates carry no intron offset: the
/// sequence name, one coordinate or two joined by `_`, and the text after them.
pub const CODING_NOTATION_RE_PATTERN: &str =
    r"\A([^:]+):c\.(\*?-?\d+)(?:_(\*?-?\d+))?((?:[^0-9+_-].*)?)\z";
pub static CODING_NOTATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(CODING_NOTATION_RE_PATTERN).expect("the pattern compiles"));

/// The text after the coordinate of a substitution of one base.
pub const ONE_BASE_SUBSTITUTION_RE_PATTERN: &str = r"\A[ACGTN]>[ACGTN]\z";
pub static ONE_BASE_SUBSTITUTION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(ONE_BASE_SUBSTITUTION_RE_PATTERN).expect("the pattern compiles")
});

/// One `_rna_edit` that changes the length of the cDNA: its first position, the positions the
/// cached numbering gains past it (negative where it loses them), what
/// `get_misalignment_offset` adds for it (TranscriptVariationAllele.pm:1614, 1617) and the
/// cached positions from its first that are the edit's own: the positions an edit adds, the
/// text of an edit that writes its span shorter, none for an edit that takes its span out,
/// whose first position is the base after the ones taken out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LengthEdit {
    pub start: i64,
    pub gained: i64,
    pub summed: i64,
    pub own: i64,
}

/// Where an edit stands against a span of cached cDNA positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Standing {
    Before,
    Inside,
    After,
}

impl LengthEdit {
    /// The edit against the cached positions `first` to `last`: before them where its own
    /// positions end before `first`, inside where the span holds one of its own positions or
    /// runs across the bases it takes out, after otherwise.
    fn standing(&self, first: i64, last: i64) -> Standing {
        if first >= self.start + self.own {
            Standing::Before
        } else if last >= self.start {
            Standing::Inside
        } else {
            Standing::After
        }
    }
}

/// What the class reads of one cached transcript: the cDNA numbers of its first and last
/// coding positions, on the cached numbering, and its edits that change the cDNA's length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LengthEditedModel {
    pub coding_start: i64,
    pub coding_end: i64,
    pub edits: Vec<LengthEdit>,
}

impl LengthEditedModel {
    /// What the edits before a row whose cDNA positions run from `first` to `last` come to:
    /// the positions the cached numbering gains over the count along the genome's exons, and
    /// the offset Ensembl VEP sums for the row, which counts an edit that starts at or before
    /// the row's cDNA start (TranscriptVariationAllele.pm:1608-1610). `None` where no edit
    /// stands before the row, or the row holds an edit's own position or runs across the
    /// bases an edit takes out.
    pub fn before(&self, first: i64, last: i64) -> Option<(i64, i64)> {
        if self
            .edits
            .iter()
            .any(|edit| edit.standing(first, last) == Standing::Inside)
        {
            return None;
        }
        let before = self
            .edits
            .iter()
            .filter(|edit| edit.standing(first, last) == Standing::Before)
            .collect::<Vec<_>>();
        (!before.is_empty()).then(|| {
            (
                before.iter().map(|edit| edit.gained).sum(),
                before.iter().map(|edit| edit.summed).sum(),
            )
        })
    }

    /// The first and last coding positions on the count along the genome's exons: the cached
    /// ones less what the cached numbering gains before each. `None` where either is an
    /// edit's own position, which has no count.
    pub fn counted_coding(&self) -> Option<(i64, i64)> {
        let counted = |cached: i64| {
            if self
                .edits
                .iter()
                .any(|edit| edit.standing(cached, cached) == Standing::Inside)
            {
                return None;
            }
            let gained: i64 = self
                .edits
                .iter()
                .filter(|edit| edit.standing(cached, cached) == Standing::Before)
                .map(|edit| edit.gained)
                .sum();
            cached.checked_sub(gained)
        };
        Some((counted(self.coding_start)?, counted(self.coding_end)?))
    }
}

/// The cDNA number of a coding-numbering coordinate between the coding positions `coding`:
/// `*5` counts from the last coding position, `-5` back from the first and `5` from it
/// (`_get_cDNA_position`, TranscriptVariationAllele.pm:2846-2868, read backwards).
fn cdna_number(coding: (i64, i64), coordinate: &str) -> Option<i64> {
    if let Some(past) = coordinate.strip_prefix('*') {
        let past = past.parse::<i64>().ok().filter(|past| *past > 0)?;
        return coding.1.checked_add(past);
    }
    let number = coordinate.parse::<i64>().ok()?;
    match number {
        0 => None,
        before if before < 0 => coding.0.checked_add(before),
        inside => coding.0.checked_add(inside - 1),
    }
}

/// The coding-numbering coordinate of cDNA number `cdna` between the coding positions
/// `coding` (`_get_cDNA_position`, TranscriptVariationAllele.pm:2846-2868).
fn coordinate(coding: (i64, i64), cdna: i64) -> String {
    if cdna > coding.1 {
        format!("*{}", cdna - coding.1)
    } else {
        (cdna + i64::from(cdna >= coding.0) - coding.0).to_string()
    }
}

/// The transcripts of a candidate JSON cache that carry an `_rna_edit` that changes the
/// cDNA's length, by stable id. A chromosome's shards are read on the first lookup of an
/// entry whose Location names its directory, every shard on the first lookup of an entry
/// without one (a VCF or JSON entry).
#[derive(Debug)]
pub struct LengthEdits {
    root: PathBuf,
    read: RefCell<IndexMap<Option<String>, HashMap<String, LengthEditedModel>>>,
}

impl LengthEdits {
    pub fn new(cache: &Path) -> Result<Self> {
        let root = cache::joined(cache, "transcripts");
        if !root.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        Ok(LengthEdits {
            root,
            read: RefCell::new(IndexMap::new()),
        })
    }

    fn models(&self, chromosome: Option<&str>) -> Result<HashMap<String, LengthEditedModel>> {
        let shards = cache::shards(&self.root, chromosome)?;
        let read = cache::map_shards(&shards, cache::threads(), |_, transcripts| {
            Ok(transcripts
                .iter()
                .filter_map(length_edited_model)
                .collect::<Vec<_>>())
        })?;
        Ok(read.into_iter().flatten().collect())
    }

    /// The entry's transcript where the cache holds it with such an edit; a Feature carrying
    /// a version the cache's stable id lacks (`--transcript_version`) is looked up without it.
    pub fn model(&self, e: &Entry) -> Result<Option<LengthEditedModel>> {
        let mut chromosome = e
            .location
            .as_deref()
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
        Ok(models
            .get(&e.feature)
            .or_else(|| models.get(VERSION_SUFFIX_RE.replace_all(&e.feature, "").as_ref()))
            .cloned())
    }
}

/// A whole number a cache record holds under `name`.
fn whole(transcript: &Value, name: &str) -> Option<i64> {
    match transcript.get(name)? {
        Value::Int(number) => Some(*number),
        _ => None,
    }
}

/// How an `_rna_edit` value is read. The value is split on blanks as
/// `get_misalignment_offset` splits it (TranscriptVariationAllele.pm:1605), an empty text
/// after the last blank being no field. A value of three fields writes its span as its text:
/// one whose span is as long as its text changes no length (:1606), and any other gains the
/// text's length less the span and is summed as the text's length (:1617). A value of two
/// fields takes its span out, which is what the numbering loses and what is summed (:1614).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EditValue {
    Changes(LengthEdit),
    KeepsLength,
    Unread,
}

fn edit_value(value: &str) -> EditValue {
    let mut fields = value.split(' ').collect::<Vec<_>>();
    while fields.last().is_some_and(|field| field.is_empty()) {
        fields.pop();
    }
    let position = |field: &str| field.parse::<i64>().ok();
    let (start, span, written) = match *fields.as_slice() {
        [start, end] => match (position(start), position(end)) {
            (Some(start), Some(end)) => (start, end - start + 1, None),
            _ => return EditValue::Unread,
        },
        [start, end, text] => match (position(start), position(end)) {
            (Some(start), Some(end)) => (
                start,
                end - start + 1,
                i64::try_from(text.chars().count()).ok(),
            ),
            _ => return EditValue::Unread,
        },
        _ => return EditValue::Unread,
    };
    if span < 0 {
        return EditValue::Unread;
    }
    match written {
        Some(written) if written == span => EditValue::KeepsLength,
        Some(written) => EditValue::Changes(LengthEdit {
            start,
            gained: written - span,
            summed: written,
            own: if written > span {
                written - span
            } else {
                written
            },
        }),
        None if span == 0 => EditValue::Unread,
        None => EditValue::Changes(LengthEdit {
            start,
            gained: -span,
            summed: -span,
            own: 0,
        }),
    }
}

/// The stable id and the reading of a transcript that carries such an edit and states its
/// coding positions; `None` for any other record, and for one with an edit value that is
/// not read. An attribute whose description names a mismatch (`op=X`) is no edit of the sum
/// (TranscriptVariationAllele.pm:1603).
fn length_edited_model(transcript: &Value) -> Option<(String, LengthEditedModel)> {
    let stable_id = transcript.get("stable_id")?.as_str()?.to_string();
    let values = transcript
        .get("attributes")?
        .as_list()?
        .iter()
        .filter(|attribute| {
            attribute
                .get("code")
                .and_then(Value::as_str)
                .is_some_and(|code| code.starts_with("_rna_edit"))
                && !attribute
                    .get("description")
                    .and_then(Value::as_str)
                    .is_some_and(|description| description.contains("op=X"))
        })
        .map(|attribute| {
            attribute
                .get("value")
                .and_then(Value::as_str)
                .map_or(EditValue::Unread, edit_value)
        })
        .collect::<Vec<_>>();
    if values.contains(&EditValue::Unread) {
        return None;
    }
    let edits = values
        .into_iter()
        .filter_map(|value| match value {
            EditValue::Changes(edit) => Some(edit),
            _ => None,
        })
        .collect::<Vec<_>>();
    if edits.is_empty() {
        return None;
    }
    let model = LengthEditedModel {
        coding_start: whole(transcript, "cdna_coding_start")?,
        coding_end: whole(transcript, "cdna_coding_end")?,
        edits,
    };
    Some((stable_id, model))
}

/// The first and last cDNA positions a row states, both plain numbers. The text formats
/// write them in one `cDNA_position` (`23-43`, `23` for one base, `23-43/1500` under
/// --total_length); the JSON format writes `cdna_start` and `cdna_end` apart.
fn cdna_positions(e: &Entry) -> Option<(i64, i64)> {
    let (start, end) = match e.fields.get(CDNA_START_FIELDS[0]) {
        Some(position) => {
            let written = pyfmt::partition(position, "/").0;
            let (start, separator, end) = pyfmt::partition(written, "-");
            (start, if separator.is_empty() { start } else { end })
        }
        None => {
            let start = e.fields.get(CDNA_START_FIELDS[1])?.as_str();
            (
                start,
                e.fields.get(CDNA_END_FIELD).map_or(start, String::as_str),
            )
        }
    };
    let (start, end) = (start.parse::<i64>().ok()?, end.parse::<i64>().ok()?);
    Some((start.min(end), start.max(end)))
}

/// How far a row's notation can stand from the row's own positions: the length of the 3'
/// shift the row prints as `HGVS_OFFSET` and the length of its allele, which a duplication
/// counts back from the point its bases are added at. `None` for an offset that is no number.
fn notation_reach(a: &Entry) -> Option<i64> {
    let shifted = match HGVS_OFFSET_FIELDS
        .iter()
        .find_map(|name| a.fields.get(*name))
        .filter(|value| !blank(value))
    {
        Some(value) => value.parse::<i64>().ok()?.checked_abs()?,
        None => 0,
    };
    shifted.checked_add(i64::try_from(pyfmt::len(&a.allele)).ok()?)
}

/// Whether two HGVSc values of the reference entry `e` and the candidate entry `a` differ as
/// `refseq_failed_edit_model_hgvsc_offset` reads them. Both rows carry `BAM_EDIT` `FAILED`
/// and one pair of cDNA positions; the row's transcript carries, in the candidate's cache, an
/// edit that changes the cDNA's length before the row's first position; and both values name
/// one sequence in coding numbering with as many coordinates, none with an intron offset.
///
/// The candidate counts each coordinate along the genome's exons, between the coding
/// positions counted the same way, and each stands at the row: a substitution of one base
/// at the count of the row's one cDNA position, any other notation within the row's 3' shift
/// and its allele's length of the row's positions. The reference writes a substitution of
/// one base on a row with one CDS position at that CDS position, the cached model's number
/// (TranscriptVariationAllele.pm:1519-1522). It writes any other notation, and a
/// substitution of one base on a row with no CDS position, at the candidate's count plus the
/// offset it sums for the row (:1524-1525), between the cached coding positions. What
/// follows the coordinates is the same text, or differs as `use_transcript_ref_cache_cdna`
/// reads the base changed where `declared` names that class.
pub fn failed_edit_offset_cell(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
    edits: &LengthEdits,
    declared: &[&str],
) -> Result<bool> {
    if blank(expected_value) || blank(actual_value) || !transcript_row(e) {
        return Ok(false);
    }
    if !failed_edit_row(e) || !failed_edit_row(a) {
        return Ok(false);
    }
    let (Some(positions), Some(candidate_positions)) = (cdna_positions(e), cdna_positions(a))
    else {
        return Ok(false);
    };
    if positions != candidate_positions {
        return Ok(false);
    }
    let Some(model) = edits.model(e)? else {
        return Ok(false);
    };
    let (Some((gained, summed)), Some(counted)) = (
        model.before(positions.0, positions.1),
        model.counted_coding(),
    ) else {
        return Ok(false);
    };
    let cached = (model.coding_start, model.coding_end);
    let (Some(expected), Some(actual)) = (
        CODING_NOTATION_RE.captures(expected_value),
        CODING_NOTATION_RE.captures(actual_value),
    ) else {
        return Ok(false);
    };
    let part = |notation: &regex::Captures, group: usize| {
        notation.get(group).map(|found| found.as_str().to_string())
    };
    let after = |notation: &regex::Captures| part(notation, 4).unwrap_or_default();
    if part(&expected, 1) != part(&actual, 1) {
        return Ok(false);
    }
    let one_base = |notation: &regex::Captures| {
        part(notation, 3).is_none() && ONE_BASE_SUBSTITUTION_RE.is_match(&after(notation))
    };
    let numbered = if positions.0 == positions.1 && one_base(&expected) && one_base(&actual) {
        // A substitution of one base stands at its row: the candidate's at the count of the
        // row's cDNA position, the reference's at the row's CDS position where it has one.
        let count = positions.0 - gained;
        let reference = match (cds_bounds(e), cds_bounds(a)) {
            (None, None) => coordinate(cached, count + summed),
            (Some((start, end)), Some(candidate)) if start == end && candidate == (start, end) => {
                start.to_string()
            }
            _ => return Ok(false),
        };
        part(&actual, 2) == Some(coordinate(counted, count))
            && part(&expected, 2) == Some(reference)
    } else {
        let Some(reach) = notation_reach(a) else {
            return Ok(false);
        };
        let at_the_row = (positions.0 - gained - reach)..=(positions.1 - gained + reach);
        [2, 3].into_iter().all(
            |group| match (part(&expected, group), part(&actual, group)) {
                (None, None) => true,
                (Some(expected), Some(actual)) => {
                    cdna_number(counted, &actual).is_some_and(|count| {
                        at_the_row.contains(&count)
                            && coordinate(counted, count) == actual
                            && coordinate(cached, count + summed) == expected
                    })
                }
                _ => false,
            },
        )
    };
    if !numbered {
        return Ok(false);
    }
    if after(&expected) == after(&actual) {
        return Ok(true);
    }
    // The base changed is another difference of such a row, read by its own class on the
    // reference's notation written at the candidate's coordinates.
    let coordinates = actual_value
        .strip_suffix(after(&actual).as_str())
        .unwrap_or(actual_value);
    let renumbered = format!("{coordinates}{}", after(&expected));
    Ok(declared.contains(&USE_TRANSCRIPT_REF_CACHE_CDNA)
        && failed_edit_reference_cell(e, a, "HGVSc", &renumbered, actual_value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::testing::{plant, scratch};
    use crate::classes::REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET;
    use crate::rules::{field_class_accepts, RuleContext, Side};

    /// A cache of one chromosome whose models state their coding sequence at cDNA 51 to 301.
    /// `NM_1.1` writes cDNA 100 as two bases; `NM_2.1` writes cDNA 100 to 102 as one base;
    /// `NM_3.1` adds two bases before cDNA 101 and takes cDNA 120 and 121 out, beside a
    /// substitution, a mismatch attribute and an attribute of another code; `NM_5.1` takes
    /// cDNA 100 out by a value that ends in a blank. `NM_4.1` states no coding positions and
    /// `NM_6.1` carries a value of one field.
    fn planted(name: &str) -> PathBuf {
        let root = scratch(name);
        plant(
            &root,
            "transcripts/7/1-1000000.json",
            r#"[{"stable_id": "NM_1.1", "cdna_coding_start": 51, "cdna_coding_end": 301,
                 "attributes": [{"code": "_rna_edit", "value": "100 100 NN"}]},
                {"stable_id": "NM_2.1", "cdna_coding_start": 51, "cdna_coding_end": 301,
                 "attributes": [{"code": "_rna_edit", "value": "100 102 N"}]},
                {"stable_id": "NM_3.1", "cdna_coding_start": 51, "cdna_coding_end": 301,
                 "attributes": [{"code": "_rna_edit", "value": "101 100 NN"},
                                {"code": "_rna_edit", "value": "120 121"},
                                {"code": "_rna_edit", "value": "130 131 NN"},
                                {"code": "_rna_edit", "value": "140 140 NN", "description": "op=X"},
                                {"code": "other", "value": "150 150 NN"}]},
                {"stable_id": "NM_4.1", "attributes": [{"code": "_rna_edit", "value": "100 100 NN"}]},
                {"stable_id": "NM_5.1", "cdna_coding_start": 51, "cdna_coding_end": 301,
                 "attributes": [{"code": "_rna_edit", "value": "100 100 "}]},
                {"stable_id": "NM_6.1", "cdna_coding_start": 51, "cdna_coding_end": 301,
                 "attributes": [{"code": "_rna_edit", "value": "100 100 NN"},
                                {"code": "_rna_edit", "value": "100"}]},
                null]"#,
        );
        root
    }

    /// A row of `feature` on chromosome 7 carrying `BAM_EDIT` `status`, `cdna` as its
    /// `cDNA_position` and, where given, `cds` as its `CDS_position`.
    fn row(feature: &str, status: &str, cdna: &str, cds: Option<&str>) -> Entry {
        let mut entry = Entry {
            key: "7:500|A|NM|Transcript".to_string(),
            location: Some("7:500".to_string()),
            allele: "A".to_string(),
            feature: feature.to_string(),
            consequence: "3_prime_UTR_variant".to_string(),
            ..Entry::default()
        };
        entry
            .fields
            .insert("BAM_EDIT".to_string(), status.to_string());
        entry
            .fields
            .insert("cDNA_position".to_string(), cdna.to_string());
        if let Some(cds) = cds {
            entry
                .fields
                .insert("CDS_position".to_string(), cds.to_string());
        }
        entry
    }

    fn edit(start: i64, gained: i64, summed: i64, own: i64) -> LengthEdit {
        LengthEdit {
            start,
            gained,
            summed,
            own,
        }
    }

    #[test]
    fn the_edits_that_change_the_length_are_read_from_the_cache() {
        let root = planted("length-edits");
        let edits = LengthEdits::new(&root).unwrap();
        let model = |feature: &str, location: Option<&str>| {
            let mut entry = row(feature, "FAILED", "400", None);
            entry.location = location.map(str::to_string);
            edits.model(&entry).unwrap()
        };
        let lengthened = LengthEditedModel {
            coding_start: 51,
            coding_end: 301,
            edits: vec![edit(100, 1, 2, 1)],
        };
        assert_eq!(model("NM_1.1", Some("7:500")), Some(lengthened.clone()));
        assert_eq!(model("NM_1.1", None), Some(lengthened.clone()));
        assert_eq!(model("NM_1.1", Some("nowhere:1")), Some(lengthened.clone()));
        let edits_of = |feature: &str| model(feature, Some("7:500")).map(|model| model.edits);
        assert_eq!(edits_of("NM_2.1"), Some(vec![edit(100, -2, 1, 1)]));
        assert_eq!(
            edits_of("NM_3.1"),
            Some(vec![edit(101, 2, 2, 2), edit(120, -2, -2, 0)])
        );
        assert_eq!(edits_of("NM_5.1"), Some(vec![edit(100, -1, -1, 0)]));
        for feature in ["NM_4.1", "NM_6.1", "NM_9.1"] {
            assert_eq!(model(feature, Some("7:500")), None, "{feature}");
        }
        // The coding positions on the count along the genome's exons, and a row's sums.
        assert_eq!(lengthened.counted_coding(), Some((51, 300)));
        assert_eq!(lengthened.before(400, 400), Some((1, 2)));
        assert_eq!(lengthened.before(101, 105), Some((1, 2)));
        for (first, last) in [(90, 90), (99, 100), (100, 100), (100, 104)] {
            assert_eq!(lengthened.before(first, last), None, "{first}-{last}");
        }
        let at_the_coding_end = LengthEditedModel {
            coding_end: 100,
            ..lengthened
        };
        assert_eq!(at_the_coding_end.counted_coding(), None);
        // An edit that takes a position out has no position of its own: its first position
        // is the base after the one taken out, a row there is past it, a row that runs
        // across the base taken out is not, and a coding end there has its count.
        let taken_out = LengthEditedModel {
            coding_start: 51,
            coding_end: 100,
            edits: vec![edit(100, -1, -1, 0)],
        };
        assert_eq!(taken_out.before(100, 100), Some((-1, -1)));
        assert_eq!(taken_out.before(100, 104), Some((-1, -1)));
        assert_eq!(taken_out.before(99, 100), None);
        assert_eq!(taken_out.before(99, 99), None);
        assert_eq!(taken_out.counted_coding(), Some((51, 101)));
        let refused = LengthEdits::new(&root.join("absent"))
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
    fn a_position_past_a_length_edit_reads_under_the_class_and_no_other_pair_does() {
        let root = planted("length-edits-class");
        let cached = Side {
            length_edits: Some(LengthEdits::new(&root).unwrap()),
            ..Side::default()
        };
        let bare = Side::default();
        let accepts = |side: &Side,
                       declared: &[&str],
                       rows: (&Entry, &Entry),
                       reference: &str,
                       candidate: &str| {
            let context = RuleContext {
                record_rows: None,
                candidate_rows: None,
                declared,
                assembly: None,
                dbnsfp_columns: &[],
                side,
            };
            field_class_accepts(
                REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET,
                rows.0,
                "HGVSc",
                reference,
                candidate,
                Some(rows.1),
                &context,
            )
            .unwrap()
        };
        // `NM_1.1` gains one position at its edit and Ensembl VEP sums two. In the 3' UTR
        // the candidate's count from the stop codon is the cached one and the reference
        // stands one above; a notation that is not a substitution of one base stands two
        // above in the coding sequence, the last coding bases carried over the stop codon.
        let in_the_utr = row("NM_1.1", "FAILED", "400", None);
        let coding = row("NM_1.1", "FAILED", "250", Some("200"));
        let coding_insertion = row("NM_1.1", "FAILED", "250-251", Some("200-201"));
        let last_coding = row("NM_1.1", "FAILED", "301", Some("251"));
        let last_coding_insertion = row("NM_1.1", "FAILED", "300-301", Some("250-251"));
        let utr_insertion = row("NM_1.1", "FAILED", "400-401", None);
        for (rows, reference, candidate) in [
            (&in_the_utr, "NM_1.1:c.*100del", "NM_1.1:c.*99del"),
            (&in_the_utr, "NM_1.1:c.*100T>A", "NM_1.1:c.*99T>A"),
            (
                &utr_insertion,
                "NM_1.1:c.*100_*101insC",
                "NM_1.1:c.*99_*100insC",
            ),
            (&coding, "NM_1.1:c.201del", "NM_1.1:c.199del"),
            (&last_coding, "NM_1.1:c.*1del", "NM_1.1:c.250del"),
            (
                &coding_insertion,
                "NM_1.1:c.201_202insA",
                "NM_1.1:c.199_200insA",
            ),
            (
                &last_coding_insertion,
                "NM_1.1:c.251_*1insA",
                "NM_1.1:c.249_250insA",
            ),
            // A substitution of one base on a row with one CDS position stands at that CDS
            // position, one above the candidate's count.
            (&coding, "NM_1.1:c.200A>C", "NM_1.1:c.199A>C"),
        ] {
            let pair = (rows, rows);
            assert!(
                accepts(&cached, &[], pair, reference, candidate),
                "{reference}"
            );
            assert!(
                !accepts(&bare, &[], pair, reference, candidate),
                "{reference}"
            );
            assert!(
                !accepts(&cached, &[], pair, candidate, reference),
                "{reference}"
            );
        }
        // Another distance, the offset's distance on a substitution of one base that has a
        // CDS position, the CDS position's distance on a deletion, another text after the
        // coordinates, another sequence, an intron offset, `n.` numbering, one coordinate
        // against two, no position and a candidate written past its own coding end stay
        // differences. So does a pair at the class's distance that stands off the row: a
        // notation is taken at its row's positions alone.
        for (rows, reference, candidate) in [
            (&coding, "NM_1.1:c.202del", "NM_1.1:c.199del"),
            (&coding, "NM_1.1:c.199del", "NM_1.1:c.199del"),
            (&coding, "NM_1.1:c.201A>C", "NM_1.1:c.199A>C"),
            (&coding, "NM_1.1:c.200A>C", "NM_1.1:c.198A>C"),
            (&coding, "NM_1.1:c.200del", "NM_1.1:c.199del"),
            (&coding, "NM_1.1:c.201dup", "NM_1.1:c.199del"),
            (&coding, "NM_9.1:c.201del", "NM_1.1:c.199del"),
            (&coding, "NM_1.1:c.201+2del", "NM_1.1:c.199+2del"),
            (&coding, "NM_1.1:n.252del", "NM_1.1:n.250del"),
            (&coding, "NM_1.1:c.201_202del", "NM_1.1:c.199del"),
            (&coding, "NM_1.1:c.2del", "NM_1.1:c.0del"),
            (&last_coding, "NM_1.1:c.*2del", "NM_1.1:c.251del"),
            (&in_the_utr, "NM_1.1:c.*101del", "NM_1.1:c.*99del"),
            (&in_the_utr, "NM_1.1:c.*6del", "NM_1.1:c.*5del"),
            (&in_the_utr, "NM_1.1:c.*6T>A", "NM_1.1:c.*5T>A"),
            (&coding, "NM_1.1:c.211del", "NM_1.1:c.209del"),
        ] {
            assert!(
                !accepts(&cached, &[], (rows, rows), reference, candidate),
                "{reference}"
            );
        }
        // A notation the 3' shift moved stands off its row by the offset the row prints, and
        // a duplication counts its allele back from there.
        let mut shifted = row("NM_1.1", "FAILED", "250", Some("200"));
        shifted
            .fields
            .insert("HGVS_OFFSET".to_string(), "9".to_string());
        assert!(accepts(
            &cached,
            &[],
            (&shifted, &shifted),
            "NM_1.1:c.210del",
            "NM_1.1:c.208del"
        ));
        assert!(!accepts(
            &cached,
            &[],
            (&shifted, &shifted),
            "NM_1.1:c.213del",
            "NM_1.1:c.211del"
        ));
        shifted
            .fields
            .insert("HGVS_OFFSET".to_string(), "soon".to_string());
        assert!(!accepts(
            &cached,
            &[],
            (&shifted, &shifted),
            "NM_1.1:c.201del",
            "NM_1.1:c.199del"
        ));
        // A row before the edit, two rows that hold the edit's own position, a row that is
        // not FAILED on both sides, two rows at other positions and a row with an open side.
        let before = row("NM_1.1", "FAILED", "90", Some("40"));
        let across = row("NM_1.1", "FAILED", "99-100", Some("49-50"));
        let at = row("NM_1.1", "FAILED", "100", Some("50"));
        let corrected = row("NM_1.1", "OK", "250", Some("200"));
        let elsewhere = row("NM_1.1", "FAILED", "251", Some("200"));
        let open = row("NM_1.1", "FAILED", "?-250", Some("200"));
        for rows in [
            (&before, &before),
            (&across, &across),
            (&at, &at),
            (&corrected, &corrected),
            (&coding, &corrected),
            (&coding, &elsewhere),
            (&open, &open),
        ] {
            assert!(!accepts(
                &cached,
                &[],
                rows,
                "NM_1.1:c.201del",
                "NM_1.1:c.199del"
            ));
        }
        // An edit that writes three positions as one base loses two positions and is summed
        // as one: the count stands two above the cached number, the reference three above
        // the count.
        let shortened = row("NM_2.1", "FAILED", "400", None);
        let pair = (&shortened, &shortened);
        assert!(accepts(
            &cached,
            &[],
            pair,
            "NM_2.1:c.*102del",
            "NM_2.1:c.*99del"
        ));
        assert!(!accepts(
            &cached,
            &[],
            pair,
            "NM_2.1:c.*100del",
            "NM_2.1:c.*99del"
        ));
        // An edit that takes a position out is summed as the numbering loses it: the
        // reference stands one below the count in the coding sequence, and at the count in
        // the 3' UTR, where the two values are equal. The row at the edit's first position,
        // the base after the one taken out, reads as any row past it; a row that runs
        // across the base taken out stays a difference.
        let taken_out = row("NM_5.1", "FAILED", "250", Some("200"));
        let first_past = row("NM_5.1", "FAILED", "100", Some("50"));
        let across_the_base = row("NM_5.1", "FAILED", "99-100", Some("49-50"));
        for (rows, reference, candidate) in [
            (&taken_out, "NM_5.1:c.200del", "NM_5.1:c.201del"),
            (&taken_out, "NM_5.1:c.200A>C", "NM_5.1:c.201A>C"),
            (&first_past, "NM_5.1:c.50del", "NM_5.1:c.51del"),
            (&first_past, "NM_5.1:c.50A>C", "NM_5.1:c.51A>C"),
        ] {
            let pair = (rows, rows);
            assert!(
                accepts(&cached, &[], pair, reference, candidate),
                "{reference}"
            );
            assert!(
                !accepts(&cached, &[], pair, candidate, reference),
                "{reference}"
            );
        }
        assert!(!accepts(
            &cached,
            &[],
            (&across_the_base, &across_the_base),
            "NM_5.1:c.49_50del",
            "NM_5.1:c.49_51del"
        ));
        // Between two edits a row takes the first alone, and past both their sum.
        let between = row("NM_3.1", "FAILED", "110", Some("60"));
        assert!(accepts(
            &cached,
            &[],
            (&between, &between),
            "NM_3.1:c.60del",
            "NM_3.1:c.58del"
        ));
        let past_both = row("NM_3.1", "FAILED", "250", Some("200"));
        assert!(!accepts(
            &cached,
            &[],
            (&past_both, &past_both),
            "NM_3.1:c.202del",
            "NM_3.1:c.200del"
        ));
        // The JSON format states the positions apart.
        let mut json = row("NM_1.1", "FAILED", "250", None);
        json.fields.shift_remove("cDNA_position");
        for (name, value) in [
            ("cdna_start", "250"),
            ("cdna_end", "250"),
            ("cds_start", "200"),
            ("cds_end", "200"),
        ] {
            json.fields.insert(name.to_string(), value.to_string());
        }
        for (reference, candidate) in [
            ("NM_1.1:c.201del", "NM_1.1:c.199del"),
            ("NM_1.1:c.200A>C", "NM_1.1:c.199A>C"),
        ] {
            assert!(
                accepts(&cached, &[], (&json, &json), reference, candidate),
                "{reference}"
            );
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_substitution_whose_base_differs_too_is_read_only_beside_the_class_of_the_base() {
        let root = planted("length-edits-base");
        let cached = Side {
            length_edits: Some(LengthEdits::new(&root).unwrap()),
            ..Side::default()
        };
        // The reference reads the base one position on (`USED_REF` `G`), the candidate the
        // genome's (`T`), on the forward strand.
        let with = |used: &str| {
            let mut entry = row("NM_1.1", "FAILED", "400", None);
            for (name, value) in [("GIVEN_REF", "T"), ("USED_REF", used), ("STRAND", "1")] {
                entry.fields.insert(name.to_string(), value.to_string());
            }
            entry
        };
        let (reference, candidate) = (with("G"), with("T"));
        let accepts = |declared: &[&str], expected: &str, actual: &str| {
            let context = RuleContext {
                record_rows: None,
                candidate_rows: None,
                declared,
                assembly: None,
                dbnsfp_columns: &[],
                side: &cached,
            };
            field_class_accepts(
                REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET,
                &reference,
                "HGVSc",
                expected,
                actual,
                Some(&candidate),
                &context,
            )
            .unwrap()
        };
        let both = [
            REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET,
            USE_TRANSCRIPT_REF_CACHE_CDNA,
        ];
        assert!(accepts(&both, "NM_1.1:c.*100G>A", "NM_1.1:c.*99T>A"));
        // The base is the other class's, and without it the cell stays a difference.
        assert!(!accepts(
            &[REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET],
            "NM_1.1:c.*100G>A",
            "NM_1.1:c.*99T>A"
        ));
        // A base that is not the reference's `USED_REF`, and another distance, stay
        // differences beside both classes.
        assert!(!accepts(&both, "NM_1.1:c.*100C>A", "NM_1.1:c.*99T>A"));
        assert!(!accepts(&both, "NM_1.1:c.*101G>A", "NM_1.1:c.*99T>A"));
        // The same base on both sides needs the position class alone.
        assert!(accepts(
            &[REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET],
            "NM_1.1:c.*100T>A",
            "NM_1.1:c.*99T>A"
        ));
        std::fs::remove_dir_all(root).unwrap();
    }
}
