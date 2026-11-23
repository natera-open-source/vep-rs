// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A regulatory shard as JSON: its shape, the names of its features' activity maps,
//! its features and the binding matrices they name.
//!
//! A shard holds chromosome => feature class => list of features. Every number
//! written passes through Perl's `+ 0` and every string through `"..."`, so a field is
//! a number or a string by the conversion's rule and not by how the cache stored it.

use std::collections::{BTreeSet, HashMap};

use crate::json::{write_plain, write_string, Object};
use crate::number::Num;
use crate::perl::{Die, OwnedStr, Perl, Sv};
use crate::storable::{Image, NodeId};
use crate::transcripts::{fetch_string, key_characters, live_entries};

/// One code per activity value of the funcgen `regulatory_activity` table.
fn activity_code(value: &[u8]) -> Option<u8> {
    match value {
        b"ACTIVE" => Some(b'A'),
        b"INACTIVE" => Some(b'I'),
        b"POISED" => Some(b'P'),
        b"REPRESSED" => Some(b'R'),
        b"NA" => Some(b'N'),
        _ => None,
    }
}

/// The names of the info.json `cell_types` list, as Perl characters in UTF-8 in sorted
/// order, with the place of each.
pub struct CellTypes {
    pub names: Vec<Vec<u8>>,
    index: HashMap<Vec<u8>, usize>,
}

impl CellTypes {
    pub fn new(names: BTreeSet<Vec<u8>>) -> Self {
        let names: Vec<Vec<u8>> = names.into_iter().collect();
        let index = names
            .iter()
            .enumerate()
            .map(|(place, name)| (name.clone(), place))
            .collect();
        Self { names, index }
    }
}

/// One binding matrix a motif named: its stable id as Perl characters in UTF-8 and
/// its canonical JSON.
pub struct MatrixRecord {
    pub id: Vec<u8>,
    pub json: Vec<u8>,
}

/// One converted regulatory shard. The matrices are every motif's, in feature order,
/// up to the feature that failed when one did: the caller replays them against the
/// matrices of earlier shards before it looks at `error`, since a matrix is recorded
/// when its motif is met.
pub struct RegulatoryShard {
    pub json: Vec<u8>,
    pub matrices: Vec<MatrixRecord>,
    pub error: Option<Die>,
}

fn text(parts: &[&[u8]]) -> Die {
    Die(parts.concat())
}

/// The shape every regulatory shard must have: a hash of hashes of arrays.
pub fn check_shape(image: &Image) -> Result<(), Die> {
    let perl = Perl::new(image);
    let root = image.root();
    if image.entries(root).is_none() {
        return Err(Die::new("invalid regulatory Storable image"));
    }
    for chromosome in live_entries(image, root) {
        let Some(classes) = perl.hash_target(Sv::Node(chromosome.value)) else {
            return Err(Die::new("invalid regulatory chromosome entry"));
        };
        for class in live_entries(image, classes) {
            if perl.array_target(Sv::Node(class.value)).is_none() {
                return Err(Die::new("invalid regulatory feature list"));
            }
        }
    }
    Ok(())
}

/// The keys of every regulatory feature's activity map in a shard whose shape is
/// checked, as Perl characters in UTF-8.
pub fn cell_type_names(image: &Image) -> BTreeSet<Vec<u8>> {
    let perl = Perl::new(image);
    let mut names = BTreeSet::new();
    for chromosome in live_entries(image, image.root()) {
        let Some(classes) = perl.hash_target(Sv::Node(chromosome.value)) else {
            continue;
        };
        let Some(features) = perl.array_target(perl.fetch(classes, b"RegulatoryFeature")) else {
            continue;
        };
        for feature in image.items(features).unwrap_or(&[]) {
            let Some(feature) = perl.hash_target(perl.element(*feature)) else {
                continue;
            };
            let Some(map) = perl.hash_target(perl.fetch(feature, b"cell_types")) else {
                continue;
            };
            // A key stored twice is one name, and a deleted key of a restricted hash none.
            for entry in image.entries(map).unwrap_or(&[]) {
                if !entry.flags.placeholder() {
                    names.insert(key_characters(image, entry).into_owned());
                }
            }
        }
    }
    names
}

/// Converts a shard whose shape is checked.
pub fn convert(image: &Image, cell_types: &CellTypes) -> RegulatoryShard {
    let mut shard = RegulatoryShard {
        json: vec![b'['],
        matrices: Vec::new(),
        error: None,
    };
    let mut writer = Writer {
        perl: Perl::new(image),
        cell_types,
        encode_error: None,
        written: 0,
    };
    match writer.shard(&mut shard) {
        Ok(()) => {
            shard.json.push(b']');
            // A string JSON::XS refuses fails the shard when the whole list is encoded,
            // after every feature is built.
            shard.error = writer.encode_error;
        }
        Err(error) => shard.error = Some(error),
    }
    shard
}

struct Writer<'a> {
    perl: Perl<'a>,
    cell_types: &'a CellTypes,
    encode_error: Option<Die>,
    written: u64,
}

impl Writer<'_> {
    fn string(&mut self, string: &OwnedStr<'_>, out: &mut Vec<u8>) {
        if let Err(error) = write_string(string.as_perl(), out) {
            self.encode_error.get_or_insert(Die(error.0));
        }
    }

    /// `$hash->{key} + 0`.
    fn number(&self, hash: NodeId, key: &[u8]) -> Result<Num, Die> {
        Ok(self
            .perl
            .number(self.perl.fetch(hash, key))?
            .plus(Num::Int(0)))
    }

    fn shard(&mut self, shard: &mut RegulatoryShard) -> Result<(), Die> {
        let image = self.perl.image;
        let mut chromosomes = live_entries(image, image.root());
        chromosomes.sort_by_key(|entry| key_characters(image, entry).into_owned());
        for chromosome in chromosomes {
            let Some(classes) = self.perl.hash_target(Sv::Node(chromosome.value)) else {
                continue;
            };
            // RegulatoryFeature before MotifFeature: the class names in reverse order.
            let mut classes = live_entries(image, classes);
            classes
                .sort_by_key(|entry| std::cmp::Reverse(key_characters(image, entry).into_owned()));
            for class in classes {
                let kind = OwnedStr {
                    bytes: image.key(class).into(),
                    utf8: class.flags.utf8(),
                };
                let features = self
                    .perl
                    .array_target(Sv::Node(class.value))
                    .and_then(|list| image.items(list))
                    .unwrap_or(&[]);
                for feature in features {
                    self.feature(self.perl.element(*feature), &kind, shard)?;
                }
            }
        }
        Ok(())
    }

    /// One feature: a regulatory feature with its activity string, or a motif with its
    /// sequence and the id of its matrix.
    fn feature(
        &mut self,
        feature: Sv,
        kind: &OwnedStr<'_>,
        shard: &mut RegulatoryShard,
    ) -> Result<(), Die> {
        let perl = self.perl;
        let kind_text = kind.characters();
        let Some(hash) = perl.hash_target(feature) else {
            return Err(text(&[b"invalid ", &kind_text, b" object"]));
        };
        for key in [&b"stable_id"[..], b"start", b"end", b"dbID"] {
            if !perl.defined(perl.fetch(hash, key)) {
                return Err(text(&[&kind_text, b" lacks ", key]));
            }
        }
        let stable_id = perl.string(perl.fetch(hash, b"stable_id"))?;
        let id_text = stable_id.characters();
        let start = self.number(hash, b"start")?;
        let end = self.number(hash, b"end")?;
        let dbid = self.number(hash, b"dbID")?;
        if kind_text == b"RegulatoryFeature" {
            let mut biotype = perl.fetch(hash, b"feature_type");
            if let Some(object) = perl.hash_target(biotype) {
                biotype = perl.fetch(object, b"so_name");
            }
            let biotype = if perl.defined(biotype) && !perl.ref_is_true(biotype) {
                Some(perl.string(biotype)?)
            } else {
                None
            };
            let activity = self.activity(hash, &id_text)?;
            if self.written > 0 {
                shard.json.push(b',');
            }
            let mut object = Object::open(&mut shard.json);
            write_plain(&activity, object.key("activity"));
            if let Some(biotype) = &biotype {
                self.string(biotype, object.key("biotype"));
            }
            dbid.write_json(object.key("dbid"));
            end.write_json(object.key("end"));
            self.string(kind, object.key("feature_type"));
            self.string(&stable_id, object.key("stable_id"));
            start.write_json(object.key("start"));
            object.close();
        } else if kind_text == b"MotifFeature" {
            if !perl.defined(perl.fetch(hash, b"strand")) {
                return Err(text(&[b"MotifFeature ", &id_text, b" lacks strand"]));
            }
            let strand = self.number(hash, b"strand")?;
            if !perl.defined(perl.fetch(hash, b"score")) {
                return Err(text(&[b"MotifFeature ", &id_text, b" lacks score"]));
            }
            let score = self.number(hash, b"score")?;
            let sequence =
                match perl.hash_target(perl.fetch(hash, b"_variation_effect_feature_cache")) {
                    Some(cache) => perl.fetch(cache, b"seq"),
                    None => Sv::Missing,
                };
            if !perl.defined(sequence) || perl.ref_is_true(sequence) {
                return Err(text(&[b"MotifFeature ", &id_text, b" lacks its sequence"]));
            }
            let sequence = perl.string(sequence)?;
            let span = end.minus(start).plus(Num::Int(1));
            let bases = sequence.as_perl().char_len();
            if Num::Int(bases as i128).compare(span) != Some(std::cmp::Ordering::Equal) {
                return Err(text(&[
                    b"MotifFeature ",
                    &id_text,
                    b" has a sequence of ",
                    bases.to_string().as_bytes(),
                    b" bases on a span of ",
                    &span.to_perl_string(),
                ]));
            }
            let (matrix_id, matrix_json) =
                self.binding_matrix(perl.fetch(hash, b"binding_matrix"), &id_text)?;
            shard.matrices.push(MatrixRecord {
                id: matrix_id.characters(),
                json: matrix_json,
            });
            if self.written > 0 {
                shard.json.push(b',');
            }
            let mut object = Object::open(&mut shard.json);
            dbid.write_json(object.key("dbid"));
            end.write_json(object.key("end"));
            self.string(kind, object.key("feature_type"));
            self.string(&matrix_id, object.key("matrix"));
            score.write_json(object.key("score"));
            self.string(&sequence, object.key("seq"));
            self.string(&stable_id, object.key("stable_id"));
            start.write_json(object.key("start"));
            strand.write_json(object.key("strand"));
            object.close();
        } else {
            return Err(text(&[b"unknown regulatory feature type ", &kind_text]));
        }
        self.written += 1;
        Ok(())
    }

    /// The activity string: one code per name of the cell-type list, `-` where the
    /// feature's map lacks the name.
    fn activity(&self, feature: NodeId, stable_id: &[u8]) -> Result<Vec<u8>, Die> {
        let perl = self.perl;
        let image = perl.image;
        let mut codes = vec![b'-'; self.cell_types.names.len()];
        // The names whose value has no code, by place in the list; the first in
        // list order is the one reported.
        let mut refused: Vec<(usize, Die)> = Vec::new();
        if let Some(map) = perl.hash_target(perl.fetch(feature, b"cell_types")) {
            for entry in image.entries(map).unwrap_or(&[]) {
                if entry.flags.placeholder() {
                    continue;
                }
                let name = key_characters(image, entry);
                let Some(&place) = self.cell_types.index.get(name.as_ref()) else {
                    continue;
                };
                // A later entry under a key stored twice replaces the earlier one.
                refused.retain(|(earlier, _)| *earlier != place);
                codes[place] = b'-';
                let value = Sv::Node(entry.value);
                if !perl.defined(value) {
                    continue;
                }
                let value = perl.string(value)?.characters();
                match activity_code(&value) {
                    Some(code) => codes[place] = code,
                    None => refused.push((
                        place,
                        text(&[
                            b"RegulatoryFeature ",
                            stable_id,
                            b" has activity '",
                            &value,
                            b"' in ",
                            &name,
                            b"; expected one of ACTIVE, INACTIVE, NA, POISED, REPRESSED",
                        ]),
                    )),
                }
            }
        }
        match refused.into_iter().min_by_key(|(place, _)| *place) {
            Some((_, error)) => Err(error),
            None => Ok(codes),
        }
    }

    /// A motif's binding matrix: its stable id and its entry's canonical JSON.
    fn binding_matrix(
        &self,
        matrix: Sv,
        motif_id: &[u8],
    ) -> Result<(OwnedStr<'static>, Vec<u8>), Die> {
        let perl = self.perl;
        let image = perl.image;
        if !perl.defined(matrix) {
            return Err(text(&[
                b"MotifFeature ",
                motif_id,
                b" lacks binding_matrix",
            ]));
        }
        let Some(hash) = perl.hash_target(matrix) else {
            return Err(Die::new("invalid binding matrix"));
        };
        let stable_id = perl.fetch(hash, b"stable_id");
        if !perl.defined(stable_id) {
            return Err(Die::new("binding matrix lacks stable_id"));
        }
        let id = perl.string(stable_id)?;
        let id_text = id.characters();
        let optional_string = |key: &[u8]| -> Result<Option<OwnedStr<'_>>, Die> {
            let sv = perl.fetch(hash, key);
            if perl.defined(sv) {
                Ok(Some(perl.string(sv)?))
            } else {
                Ok(None)
            }
        };
        let name = optional_string(b"name")?;
        let source = optional_string(b"source")?;
        let threshold = if perl.defined(perl.fetch(hash, b"threshold")) {
            Some(self.number(hash, b"threshold")?)
        } else {
            None
        };
        if let Some(unit) = optional_string(b"unit")? {
            let unit = unit.characters();
            if unit != b"Frequencies" {
                return Err(text(&[
                    b"binding matrix ",
                    &id_text,
                    b" has unit '",
                    &unit,
                    b"'; only Frequencies converts",
                ]));
            }
        }
        let mut factors = Vec::new();
        let complexes = perl.fetch(hash, b"associated_transcription_factor_complexes");
        if perl.defined(complexes) {
            let Some(list) = perl.array_target(complexes) else {
                return Err(Die::new("invalid transcription factor complexes"));
            };
            for complex in image.items(list).unwrap_or(&[]) {
                let Some(complex) = perl.hash_target(perl.element(*complex)) else {
                    continue;
                };
                let display_name = perl.fetch(complex, b"display_name");
                if perl.defined(display_name) {
                    factors.push(perl.string(display_name)?);
                }
            }
        }
        let elements = perl.fetch(hash, b"elements");
        let length = perl.fetch(hash, b"length");
        let mut rows: Vec<[Num; 4]> = Vec::new();
        let length_value;
        if perl.defined(elements) {
            let Some(elements) = perl.hash_target(elements) else {
                return Err(Die::new("invalid binding matrix elements"));
            };
            let positions = live_entries(image, elements).len();
            let (declared, declared_text) = if perl.defined(length) {
                (perl.number(length)?, perl.string(length)?.characters())
            } else {
                (
                    Num::Int(positions as i128),
                    positions.to_string().into_bytes(),
                )
            };
            if declared.compare(Num::Int(positions as i128)) != Some(std::cmp::Ordering::Equal) {
                return Err(text(&[
                    b"binding matrix ",
                    &id_text,
                    b" declares length ",
                    &declared_text,
                    b" but has ",
                    positions.to_string().as_bytes(),
                    b" positions",
                ]));
            }
            for position in 1..=positions {
                let position_text = position.to_string().into_bytes();
                let row = fetch_string(image, elements, &OwnedStr::ascii(position_text.clone()));
                let Some(row) = perl.hash_target(row) else {
                    return Err(text(&[
                        b"binding matrix ",
                        &id_text,
                        b" lacks position ",
                        &position_text,
                    ]));
                };
                let mut counts = [Num::Int(0); 4];
                for (slot, base) in [&b"A"[..], b"C", b"G", b"T"].iter().enumerate() {
                    if !perl.defined(perl.fetch(row, base)) {
                        return Err(text(&[
                            b"binding matrix ",
                            &id_text,
                            b" lacks ",
                            base,
                            b" at position ",
                            &position_text,
                        ]));
                    }
                    counts[slot] = self.number(row, base)?;
                }
                rows.push(counts);
            }
            length_value = declared.plus(Num::Int(0));
        } else {
            if !perl.defined(length) {
                return Err(text(&[
                    b"binding matrix ",
                    &id_text,
                    b" has neither elements nor a length",
                ]));
            }
            length_value = perl.number(length)?.plus(Num::Int(0));
        }

        let mut json = Vec::new();
        let mut refused = None;
        let mut string = |value: &OwnedStr<'_>, out: &mut Vec<u8>| {
            if let Err(error) = write_string(value.as_perl(), out) {
                refused.get_or_insert(Die(error.0));
            }
        };
        let mut object = Object::open(&mut json);
        let out = object.key("frequencies");
        out.push(b'[');
        for (index, row) in rows.iter().enumerate() {
            if index > 0 {
                out.push(b',');
            }
            out.push(b'[');
            for (slot, count) in row.iter().enumerate() {
                if slot > 0 {
                    out.push(b',');
                }
                count.write_json(out);
            }
            out.push(b']');
        }
        out.push(b']');
        length_value.write_json(object.key("length"));
        if let Some(name) = &name {
            string(name, object.key("name"));
        }
        if let Some(source) = &source {
            string(source, object.key("source"));
        }
        if let Some(threshold) = threshold {
            threshold.write_json(object.key("threshold"));
        }
        let out = object.key("transcription_factors");
        out.push(b'[');
        for (index, factor) in factors.iter().enumerate() {
            if index > 0 {
                out.push(b',');
            }
            string(factor, out);
        }
        out.push(b']');
        object.close();
        // A matrix is encoded when it is recorded, while its feature is built.
        if let Some(error) = refused {
            return Err(error);
        }
        let id = OwnedStr {
            bytes: std::borrow::Cow::Owned(id.bytes.into_owned()),
            utf8: id.utf8,
        };
        Ok((id, json))
    }
}
