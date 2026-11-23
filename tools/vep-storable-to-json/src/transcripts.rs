// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A transcript shard as JSON: every field of a transcript under its own rule, written as
//! `JSON->new->utf8->canonical->allow_blessed->convert_blessed` writes the result.
//!
//! The conversion is defined as the build of every transcript of a shard and then the
//! encoding of all, so a value that cannot be walked (a reference of the wrong kind)
//! fails the shard before a value JSON::XS cannot write does. The writer here does both
//! in one pass and keeps that precedence: it remembers the first value it could not
//! write and still walks the rest of the shard for a value it cannot walk.

use std::cmp::Ordering;

use crate::json::{write_base64, write_plain, write_string, Object, PerlStr};
use crate::number::{write_g15, write_int, Num};
use crate::perl::{Die, OwnedStr, Perl, Sv};
use crate::regions::{self, Noted};
use crate::storable::{HashEntry, Image, NodeId, Value};

/// JSON::XS's default nesting limit.
const MAX_DEPTH: usize = 512;

/// One converted transcript shard.
pub struct TranscriptShard {
    /// The shard file's bytes.
    pub json: Vec<u8>,
    /// Transcripts written.
    pub transcripts: u64,
    /// Whether a SIFT matrix was written.
    pub sift: bool,
    /// Whether a PolyPhen matrix was written.
    pub polyphen: bool,
    /// The first value JSON::XS refuses. It is met when the shard is encoded, after
    /// the shard's transcripts are counted.
    pub encode_error: Option<Die>,
}

/// The walk over a transcript shard: its conversion, and what the transcripts walked
/// state of their chromosome's sequence.
pub struct Walk {
    /// The converted shard, or the die that ended the walk.
    pub shard: Result<TranscriptShard, Die>,
    /// Read from each transcript before it is converted, so a walk that dies keeps
    /// what the transcripts before the die state.
    pub sequence: Noted,
}

/// The entries of a hash as Perl's `keys` and `values` see them: placeholders skipped
/// and a key stored twice counted once, with its later value. The order is the
/// image's; Perl's own order is not reproducible.
pub fn live_entries(image: &Image, hash: NodeId) -> Vec<&HashEntry> {
    let entries = image.entries(hash).unwrap_or(&[]);
    // The index of the last entry under each key; a small hash is searched in place.
    let last: Option<std::collections::HashMap<std::borrow::Cow<'_, [u8]>, usize>> =
        (entries.len() > 8).then(|| {
            entries
                .iter()
                .enumerate()
                .filter(|(_, entry)| !entry.flags.placeholder())
                .map(|(index, entry)| (key_characters(image, entry), index))
                .collect()
        });
    let mut live: Vec<&HashEntry> = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        if entry.flags.placeholder() {
            continue;
        }
        let key = key_characters(image, entry);
        let repeated = match &last {
            Some(last) => last.get(&key) != Some(&index),
            None => entries[index + 1..]
                .iter()
                .any(|later| !later.flags.placeholder() && key_characters(image, later) == key),
        };
        if !repeated {
            live.push(entry);
        }
    }
    live
}

/// An entry's key as Perl characters in UTF-8.
pub fn key_characters<'a>(image: &'a Image, entry: &HashEntry) -> std::borrow::Cow<'a, [u8]> {
    PerlStr {
        bytes: image.key(entry),
        utf8: entry.flags.utf8(),
    }
    .characters()
}

/// A fetch by a key that is a Perl string of any characters.
pub fn fetch_string(image: &Image, hash: NodeId, key: &OwnedStr<'_>) -> Sv {
    let wanted = key.characters();
    image
        .entries(hash)
        .unwrap_or(&[])
        .iter()
        .rev()
        .find(|entry| {
            !entry.flags.placeholder() && key_characters(image, entry).as_ref() == wanted.as_slice()
        })
        .map_or(Sv::Missing, |entry| Sv::Node(entry.value))
}

struct Writer<'a> {
    perl: Perl<'a>,
    /// The first value JSON::XS would refuse, kept while the walk goes on.
    encode_error: Option<Die>,
    sift: bool,
    polyphen: bool,
    sequence: Noted,
}

/// Converts the root of a transcript shard: a hash of
/// arrays of transcripts (one array a chromosome), or an array of transcripts.
pub fn convert(image: &Image) -> Result<TranscriptShard, Die> {
    walk(image).shard
}

/// Converts the root of a transcript shard as [`convert`] does, and reads the
/// chromosome's sequence from the slice of each transcript on the way.
pub fn walk(image: &Image) -> Walk {
    let mut writer = Writer {
        perl: Perl::new(image),
        encode_error: None,
        sift: false,
        polyphen: false,
        sequence: Noted::Nothing,
    };
    let shard = shard(&mut writer, image);
    Walk {
        shard,
        sequence: writer.sequence,
    }
}

fn shard(writer: &mut Writer<'_>, image: &Image) -> Result<TranscriptShard, Die> {
    let perl = writer.perl;
    let root = image.root();
    let class = image.class(root);
    // `ref($data)` is false for a root blessed into the class `0`.
    let root_ref_is_true = class != Some(b"0");
    let is_hash = matches!(image.value(root), Value::Hash { .. });
    let is_array = matches!(image.value(root), Value::Array(_));
    let mut json = vec![b'['];
    let mut transcripts = 0u64;
    let mut write_list = |writer: &mut Writer<'_>, items: &[NodeId], json: &mut Vec<u8>| {
        for item in items {
            let transcript = perl.element(*item);
            if perl.defined(transcript) && perl.ref_is_true(transcript) {
                if transcripts > 0 {
                    json.push(b',');
                }
                writer.transcript(transcript, json)?;
                transcripts += 1;
            }
        }
        Ok::<(), Die>(())
    };
    if root_ref_is_true && (is_hash || class.is_some()) {
        if !is_hash {
            // `keys %{$data}` on a blessed value that is no hash.
            return Err(Die::new("Not a HASH reference"));
        }
        for entry in live_entries(image, root) {
            let list = Sv::Node(entry.value);
            if !perl.ref_is_true(list) {
                continue;
            }
            let Some(array) = perl.array_target(list) else {
                continue;
            };
            write_list(writer, image.items(array).unwrap_or(&[]), &mut json)?;
        }
    } else if root_ref_is_true && is_array {
        write_list(writer, image.items(root).unwrap_or(&[]), &mut json)?;
    }
    json.push(b']');
    Ok(TranscriptShard {
        json,
        transcripts,
        sift: writer.sift,
        polyphen: writer.polyphen,
        encode_error: writer.encode_error.take(),
    })
}

impl Writer<'_> {
    fn refuse(&mut self, error: Die) {
        self.encode_error.get_or_insert(error);
    }

    /// The hash `$sv->{...}` reads from, for a scalar already known to be a
    /// reference.
    fn hash(&self, sv: Sv) -> Result<NodeId, Die> {
        self.perl
            .hash_of(sv)?
            .ok_or_else(|| Die::new("Not a HASH reference"))
    }

    /// A scalar copied as it is, written as JSON::XS writes it: a string
    /// scalar as a string, an integer in decimal, a double with `%.15g`, a reference as
    /// the structure it points at.
    fn scalar(&mut self, sv: Sv, out: &mut Vec<u8>) {
        self.scalar_at(sv, out, 0);
    }

    fn scalar_at(&mut self, sv: Sv, out: &mut Vec<u8>, depth: usize) {
        let image = self.perl.image;
        match self.perl.value(sv) {
            None | Some(Value::Undef | Value::ImmortalUndef | Value::ElementUndef) => {
                out.extend_from_slice(b"null");
            }
            Some(Value::Integer(value)) => write_int(i128::from(value), out),
            Some(Value::Double(value)) => write_g15(value, out),
            Some(Value::String { bytes, utf8 }) => {
                let string = PerlStr {
                    bytes: image.bytes(bytes),
                    utf8,
                };
                if let Err(error) = write_string(string, out) {
                    self.refuse(Die(error.0));
                }
            }
            // The immortals are strings to JSON::XS, which looks at the string flag first.
            Some(Value::Yes | Value::Bool(true)) => write_plain(b"1", out),
            Some(Value::No | Value::Bool(false)) => write_plain(b"", out),
            Some(Value::Ref { target, .. }) => self.referent(target, out, depth + 1),
            Some(Value::Array(_) | Value::Hash { .. } | Value::Regexp(_)) => {
                self.refuse(Die::new(
                    "encountered a container where a scalar is stored, which JSON cannot handle",
                ));
            }
        }
    }

    /// What a copied reference points at, as JSON::XS encodes it under
    /// `allow_blessed` and `convert_blessed`: a blessed value as `null` (no class the
    /// conversion loads has a `TO_JSON`), a hash with its keys sorted, an array with
    /// `null` for an element that does not exist, a reference to the scalar 1 or 0 as
    /// `true` or `false`.
    fn referent(&mut self, target: NodeId, out: &mut Vec<u8>, depth: usize) {
        let image = self.perl.image;
        if depth > MAX_DEPTH {
            self.refuse(Die::new(
                "json text or perl structure exceeds maximum nesting level (max_depth set too low?)",
            ));
            return;
        }
        if image.class(target).is_some() {
            out.extend_from_slice(b"null");
            return;
        }
        match image.value(target) {
            Value::Hash { .. } => {
                let mut members: Vec<(std::borrow::Cow<'_, [u8]>, NodeId)> =
                    live_entries(image, target)
                        .into_iter()
                        .map(|entry| (key_characters(image, entry), entry.value))
                        .collect();
                members.sort_by(|left, right| left.0.cmp(&right.0));
                out.push(b'{');
                for (index, (key, value)) in members.iter().enumerate() {
                    if index > 0 {
                        out.push(b',');
                    }
                    let key = PerlStr {
                        bytes: key,
                        utf8: true,
                    };
                    if let Err(error) = write_string(key, out) {
                        self.refuse(Die(error.0));
                    }
                    out.push(b':');
                    self.scalar_at(Sv::Node(*value), out, depth);
                }
                out.push(b'}');
            }
            Value::Array(_) => {
                out.push(b'[');
                for (index, item) in image.items(target).unwrap_or(&[]).iter().enumerate() {
                    if index > 0 {
                        out.push(b',');
                    }
                    self.scalar_at(self.perl.element(*item), out, depth);
                }
                out.push(b']');
            }
            _ => {
                let text = self
                    .perl
                    .string(Sv::Node(target))
                    .map(|text| text.bytes.into_owned());
                match text.as_deref() {
                    Ok(b"1") => out.extend_from_slice(b"true"),
                    Ok(b"0") => out.extend_from_slice(b"false"),
                    _ => self.refuse(Die::new(
                        "cannot encode reference to scalar unless the scalar is 0 or 1",
                    )),
                }
            }
        }
    }

    /// `$out{name} = $hash->{key} if defined $hash->{key}`.
    fn copy(&mut self, object: &mut Object<'_>, name: &'static str, hash: NodeId, key: &[u8]) {
        let sv = self.perl.fetch(hash, key);
        if self.perl.defined(sv) {
            self.scalar(sv, object.key(name));
        }
    }

    /// `[map { step($_) } grep { defined $_ && ref $_ } @{$list}]`.
    fn list(
        &mut self,
        list: Sv,
        out: &mut Vec<u8>,
        step: impl Fn(&mut Self, Sv, &mut Vec<u8>) -> Result<(), Die>,
    ) -> Result<u64, Die> {
        let mut written = 0;
        out.push(b'[');
        for item in self.perl.items_of(list)? {
            if self.perl.defined(item) && self.perl.ref_is_true(item) {
                if written > 0 {
                    out.push(b',');
                }
                step(self, item, out)?;
                written += 1;
            }
        }
        out.push(b']');
        Ok(written)
    }

    fn transcript(&mut self, transcript: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let perl = self.perl;
        let hash = self.hash(transcript)?;
        // The slice is read before the transcript is converted: a transcript that ends
        // the walk has stated its sequence by then.
        if self.sequence == Noted::Nothing {
            self.sequence = match regions::of_transcript(perl, hash) {
                Ok(Some(region)) => Noted::Region(region),
                Ok(None) => Noted::Nothing,
                Err(die) => Noted::Dies(die),
            };
        }
        // Exons, translation and attributes are read in this order, before the
        // feature cache; each is held until its place among the sorted keys.
        let exons = {
            let own = perl.fetch(hash, b"exons");
            if perl.is_true(own) {
                own
            } else {
                perl.fetch(hash, b"_trans_exon_array")
            }
        };
        let mut exons_json = Vec::new();
        let has_exons = perl.ref_is(exons, b"ARRAY");
        if has_exons {
            self.list(exons, &mut exons_json, Self::exon)?;
        }
        let translation = perl.fetch(hash, b"translation");
        let mut translation_json = Vec::new();
        let has_translation = perl.ref_is_true(translation);
        if has_translation {
            self.translation(translation, &mut translation_json)?;
        }
        let attributes = perl.fetch(hash, b"attributes");
        let mut attributes_json = Vec::new();
        let has_attributes = perl.ref_is(attributes, b"ARRAY");
        if has_attributes {
            self.list(attributes, &mut attributes_json, Self::attribute)?;
        }
        let cache = {
            let own = perl.fetch(hash, b"variation_effect_feature_cache");
            if perl.is_true(own) {
                own
            } else {
                perl.fetch(hash, b"_variation_effect_feature_cache")
            }
        };

        let mut object = Object::open(out);
        if has_attributes {
            object.key("attributes").extend_from_slice(&attributes_json);
        }
        self.copy(&mut object, "bam_edit_status", hash, b"_bam_edit_status");
        self.biotype(&mut object, hash)?;
        self.copy(&mut object, "ccds", hash, b"_ccds");
        self.copy(&mut object, "cdna_coding_end", hash, b"cdna_coding_end");
        self.copy(&mut object, "cdna_coding_start", hash, b"cdna_coding_start");
        self.copy(&mut object, "dbID", hash, b"dbID");
        self.copy(&mut object, "end", hash, b"end");
        if has_exons {
            object.key("exons").extend_from_slice(&exons_json);
        }
        self.copy(&mut object, "gene_hgnc_id", hash, b"_gene_hgnc_id");
        self.copy(&mut object, "gene_phenotype", hash, b"_gene_phenotype");
        self.copy(&mut object, "gene_stable_id", hash, b"_gene_stable_id");
        self.copy(&mut object, "gene_symbol", hash, b"_gene_symbol");
        self.copy(&mut object, "gene_symbol_id", hash, b"_gene_symbol_id");
        self.copy(
            &mut object,
            "gene_symbol_source",
            hash,
            b"_gene_symbol_source",
        );
        self.copy(&mut object, "gene_version", hash, b"_gene_version");
        self.copy(&mut object, "is_canonical", hash, b"is_canonical");
        self.copy(&mut object, "protein", hash, b"_protein");
        self.copy(&mut object, "refseq", hash, b"_refseq");
        self.copy(&mut object, "source", hash, b"source");
        self.copy(&mut object, "source_cache", hash, b"_source_cache");
        self.copy(&mut object, "stable_id", hash, b"stable_id");
        self.copy(&mut object, "start", hash, b"start");
        self.copy(&mut object, "strand", hash, b"strand");
        self.copy(&mut object, "swissprot", hash, b"_swissprot");
        if has_translation {
            object
                .key("translation")
                .extend_from_slice(&translation_json);
        }
        self.copy(&mut object, "trembl", hash, b"_trembl");
        self.copy(&mut object, "uniparc", hash, b"_uniparc");
        self.copy(&mut object, "uniprot_isoform", hash, b"_uniprot_isoform");
        if perl.ref_is_true(cache) {
            let strand = perl.fetch(hash, b"strand");
            self.feature_cache(cache, strand, object.key("variation_effect_feature_cache"))?;
        }
        self.copy(&mut object, "version", hash, b"version");
        object.close();
        Ok(())
    }

    /// The `biotype` member. A biotype held as a reference to a hash, the shape of a
    /// `Bio::EnsEMBL::Biotype` object, is written as the string of its `name`, the value
    /// `Transcript::biotype` returns for it (ensembl core Transcript.pm:3398-3408), and
    /// is left out when the name is undef or a reference. Any other biotype is copied.
    fn biotype(&mut self, object: &mut Object<'_>, transcript: NodeId) -> Result<(), Die> {
        let perl = self.perl;
        let Some(held) = perl.hash_target(perl.fetch(transcript, b"biotype")) else {
            self.copy(object, "biotype", transcript, b"biotype");
            return Ok(());
        };
        let name = perl.fetch(held, b"name");
        if perl.defined(name) && !perl.ref_is_true(name) {
            let name = perl.string(name)?;
            if let Err(error) = write_string(name.as_perl(), object.key("biotype")) {
                self.refuse(Die(error.0));
            }
        }
        Ok(())
    }

    fn exon(&mut self, exon: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let hash = self.hash(exon)?;
        let mut object = Object::open(out);
        self.copy(&mut object, "end", hash, b"end");
        self.copy(&mut object, "end_phase", hash, b"end_phase");
        self.copy(&mut object, "phase", hash, b"phase");
        self.copy(&mut object, "rank", hash, b"rank");
        self.copy(&mut object, "stable_id", hash, b"stable_id");
        self.copy(&mut object, "start", hash, b"start");
        object.close();
        Ok(())
    }

    fn attribute(&mut self, attribute: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let hash = self.hash(attribute)?;
        let mut object = Object::open(out);
        self.copy(&mut object, "code", hash, b"code");
        self.copy(&mut object, "name", hash, b"name");
        self.copy(&mut object, "value", hash, b"value");
        object.close();
        Ok(())
    }

    fn intron(&mut self, intron: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let hash = self.hash(intron)?;
        let mut object = Object::open(out);
        self.copy(&mut object, "end", hash, b"end");
        self.copy(&mut object, "start", hash, b"start");
        object.close();
        Ok(())
    }

    fn mapper_pair(&mut self, pair: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let hash = self.hash(pair)?;
        let mut object = Object::open(out);
        self.copy(&mut object, "from_end", hash, b"from_end");
        self.copy(&mut object, "from_start", hash, b"from_start");
        self.copy(&mut object, "ori", hash, b"ori");
        self.copy(&mut object, "to_end", hash, b"to_end");
        self.copy(&mut object, "to_start", hash, b"to_start");
        object.close();
        Ok(())
    }

    /// A translation: its own fields, and `seq` as the string
    /// itself or, for a sequence object, its `seq` or its `primary_seq`'s. The
    /// object's `seq` method is not called: the conversion loads no class that has one.
    fn translation(&mut self, translation: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let perl = self.perl;
        let hash = self.hash(translation)?;
        let sequence = perl.fetch(hash, b"seq");
        // `None`: no `seq` key is written; `Some(Sv::Missing)`: it is written as null.
        let mut written: Option<Sv> = None;
        if perl.defined(sequence) {
            if perl.ref_is_true(sequence) {
                let object = self.hash(sequence)?;
                let own = perl.fetch(object, b"seq");
                written = Some(if perl.defined(own) {
                    own
                } else {
                    // `$seq->{primary_seq}{seq}`: an absent `primary_seq` reads as undef.
                    perl.get(perl.fetch(object, b"primary_seq"), b"seq")?
                });
            } else {
                written = Some(sequence);
            }
        }
        let mut object = Object::open(out);
        self.copy(&mut object, "dbID", hash, b"dbID");
        self.copy(&mut object, "end", hash, b"end");
        if let Some(sequence) = written {
            self.scalar(sequence, object.key("seq"));
        }
        self.copy(&mut object, "stable_id", hash, b"stable_id");
        self.copy(&mut object, "start", hash, b"start");
        self.copy(&mut object, "version", hash, b"version");
        object.close();
        Ok(())
    }

    /// A sequence value: the scalar itself, or the `seq` of a sequence object or of
    /// its `primary_seq`; `None` when there is no string to take.
    fn sequence_string(&self, value: Sv) -> Result<Option<Sv>, Die> {
        let perl = self.perl;
        if !perl.defined(value) {
            return Ok(None);
        }
        if !perl.ref_is_true(value) {
            return Ok(Some(value));
        }
        let hash_like = |sv: Sv| {
            perl.ref_is(sv, b"HASH")
                || perl
                    .target(sv)
                    .is_some_and(|target| perl.image.class(target).is_some())
                    && perl.hash_target(sv).is_some()
        };
        if hash_like(value) {
            let own = perl.get(value, b"seq")?;
            if perl.defined(own) && !perl.ref_is_true(own) {
                return Ok(Some(own));
            }
            let primary = perl.get(value, b"primary_seq")?;
            if perl.ref_is_true(primary) && hash_like(primary) {
                let inner = perl.get(primary, b"seq")?;
                if perl.defined(inner) && !perl.ref_is_true(inner) {
                    return Ok(Some(inner));
                }
            }
        }
        Ok(None)
    }

    /// The variation effect feature cache of a transcript.
    fn feature_cache(&mut self, cache: Sv, strand: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let perl = self.perl;
        let outer = self.hash(cache)?;
        let image = perl.image;
        let exists = |key: &[u8]| image.fetch(outer, key).is_some();
        let mut cache = cache;
        if exists(b"_variation_effect_feature_cache") || exists(b"variation_effect_feature_cache") {
            let underscored = perl.fetch(outer, b"_variation_effect_feature_cache");
            let plain = perl.fetch(outer, b"variation_effect_feature_cache");
            if perl.is_true(underscored) {
                cache = underscored;
            } else if perl.is_true(plain) {
                cache = plain;
            }
        }
        let hash = self.hash(cache)?;

        let five_prime = {
            let value = perl.fetch(hash, b"five_prime_utr");
            self.sequence_string(value)?
        };
        let three_prime = {
            let value = perl.fetch(hash, b"three_prime_utr");
            self.sequence_string(value)?
        };
        // The parts are built in the order seq_edits, introns, sorted_exons, mapper,
        // protein_features, prediction matrices; the two that sort after the large
        // ones are held.
        let mut seq_edits_json = Vec::new();
        let seq_edits = perl.fetch(hash, b"seq_edits");
        if perl.ref_is(seq_edits, b"ARRAY") {
            self.list(seq_edits, &mut seq_edits_json, Self::seq_edit)?;
        } else {
            seq_edits_json.extend_from_slice(b"[]");
        }

        let mut object = Object::open(out);
        self.copy(&mut object, "codon_table", hash, b"codon_table");
        if let Some(sequence) = five_prime {
            self.scalar(sequence, object.key("five_prime_utr"));
        }
        let introns = perl.fetch(hash, b"introns");
        if perl.ref_is(introns, b"ARRAY") {
            self.list(introns, object.key("introns"), Self::intron)?;
        }
        let sorted_exons = perl.fetch(hash, b"sorted_exons");
        let mut sorted_exons_json = Vec::new();
        let has_sorted_exons = perl.ref_is(sorted_exons, b"ARRAY");
        if has_sorted_exons {
            self.list(sorted_exons, &mut sorted_exons_json, Self::exon)?;
        }
        let mapper = perl.fetch(hash, b"mapper");
        if perl.ref_is_true(mapper) {
            self.mapper(mapper, sorted_exons, strand, object.key("mapper"))?;
        }
        self.copy(&mut object, "peptide", hash, b"peptide");
        let protein_features = perl.fetch(hash, b"protein_features");
        if perl.ref_is(protein_features, b"ARRAY") {
            self.list(
                protein_features,
                object.key("protein_features"),
                Self::protein_feature,
            )?;
        }
        let predictions = perl.fetch(hash, b"protein_function_predictions");
        if perl.ref_is(predictions, b"HASH")
            || (perl.ref_is_true(predictions) && perl.hash_target(predictions).is_some())
        {
            let predictions = self.hash(predictions)?;
            // Read in this order: sift, polyphen_humvar, polyphen_humdiv.
            let mut held: [Option<Vec<u8>>; 3] = [None, None, None];
            for (slot, key) in [&b"sift"[..], b"polyphen_humvar", b"polyphen_humdiv"]
                .iter()
                .enumerate()
            {
                let matrix = perl.fetch(predictions, key);
                if perl.ref_is_true(matrix) {
                    let mut json = Vec::new();
                    self.prediction_matrix(matrix, slot == 0, &mut json)?;
                    held[slot] = Some(json);
                }
            }
            if held.iter().any(Option::is_some) {
                let mut inner = Object::open(object.key("protein_function_predictions"));
                let [sift, humvar, humdiv] = held;
                if let Some(json) = humdiv {
                    inner.key("polyphen_humdiv").extend_from_slice(&json);
                }
                if let Some(json) = humvar {
                    inner.key("polyphen_humvar").extend_from_slice(&json);
                }
                if let Some(json) = sift {
                    inner.key("sift").extend_from_slice(&json);
                }
                inner.close();
            }
        }
        object.key("seq_edits").extend_from_slice(&seq_edits_json);
        if has_sorted_exons {
            object
                .key("sorted_exons")
                .extend_from_slice(&sorted_exons_json);
        }
        self.copy(&mut object, "spliced_seq", hash, b"spliced_seq");
        if let Some(sequence) = three_prime {
            self.scalar(sequence, object.key("three_prime_utr"));
        }
        self.copy(&mut object, "translateable_seq", hash, b"translateable_seq");
        object.close();
        Ok(())
    }

    /// A sequence edit: `start` and `end` as numbers, `alt_seq` as a string, and the
    /// length difference `SeqEdit::length_diff` computes.
    fn seq_edit(&mut self, edit: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let perl = self.perl;
        let hash = self.hash(edit)?;
        let number = |key: &[u8]| -> Result<Option<Num>, Die> {
            let sv = perl.fetch(hash, key);
            if !perl.defined(sv) {
                return Ok(None);
            }
            Ok(Some(perl.number(sv)?.plus(Num::Int(0))))
        };
        let start = number(b"start")?;
        let end = number(b"end")?;
        let alt_seq = perl.fetch(hash, b"alt_seq");
        let alt_seq = if perl.defined(alt_seq) {
            perl.string(alt_seq)?
        } else {
            OwnedStr::ascii(Vec::new())
        };
        let length_diff = match (start, end) {
            (Some(start), Some(end)) => Num::Int(alt_seq.as_perl().char_len() as i128)
                .minus(end.minus(start).plus(Num::Int(1))),
            _ => Num::Int(0),
        };
        let mut object = Object::open(out);
        if let Err(error) = write_string(alt_seq.as_perl(), object.key("alt_seq")) {
            self.refuse(Die(error.0));
        }
        self.copy(&mut object, "code", hash, b"code");
        if let Some(end) = end {
            end.write_json(object.key("end"));
        }
        length_diff.write_json(object.key("length_diff"));
        self.copy(&mut object, "name", hash, b"name");
        if let Some(start) = start {
            start.write_json(object.key("start"));
        }
        object.close();
        Ok(())
    }

    /// A protein feature: its fields and its analysis as the display
    /// label of an analysis object, or the scalar itself.
    fn protein_feature(&mut self, feature: Sv, out: &mut Vec<u8>) -> Result<(), Die> {
        let perl = self.perl;
        let hash = self.hash(feature)?;
        let analysis = perl.fetch(hash, b"analysis");
        let label = match perl.hash_target(analysis) {
            Some(object) => {
                let label = perl.fetch(object, b"_display_label");
                perl.defined(label).then_some(label)
            }
            None => (perl.defined(analysis) && !perl.ref_is_true(analysis)).then_some(analysis),
        };
        let mut object = Object::open(out);
        if let Some(label) = label {
            self.scalar(label, object.key("analysis"));
        }
        self.copy(&mut object, "end", hash, b"end");
        self.copy(&mut object, "hseqname", hash, b"hseqname");
        self.copy(&mut object, "start", hash, b"start");
        object.close();
        Ok(())
    }

    /// A prediction matrix: its fields and its blob in base64, gzipped
    /// first when it does not carry the gzip magic.
    fn prediction_matrix(&mut self, matrix: Sv, sift: bool, out: &mut Vec<u8>) -> Result<(), Die> {
        let perl = self.perl;
        let hash = self.hash(matrix)?;
        let blob = perl.fetch(hash, b"matrix");
        let mut encoded = None;
        if perl.defined(blob) {
            let blob = perl.string(blob)?;
            let bytes = downgrade(&blob)?;
            let mut text = Vec::new();
            if bytes.starts_with(&[0x1f, 0x8b]) {
                write_base64(&bytes, &mut text);
            } else {
                write_base64(&crate::gzip::mem_gzip(&bytes), &mut text);
            }
            encoded = Some(text);
            if sift {
                self.sift = true;
            } else {
                self.polyphen = true;
            }
        }
        let mut object = Object::open(out);
        self.copy(&mut object, "analysis", hash, b"analysis");
        self.copy(&mut object, "peptide_length", hash, b"peptide_length");
        if let Some(text) = encoded {
            write_plain(&text, object.key("predictions_data"));
        }
        self.copy(&mut object, "sub_analysis", hash, b"sub_analysis");
        self.copy(&mut object, "translation_md5", hash, b"translation_md5");
        object.close();
        Ok(())
    }

    /// A mapper: its own fields and its pairs from the first of three
    /// sources that has any: the mapper's `pairs` list, the pairs of its
    /// `exon_coord_mapper` object, the exons.
    fn mapper(
        &mut self,
        mapper: Sv,
        sorted_exons: Sv,
        strand: Sv,
        out: &mut Vec<u8>,
    ) -> Result<(), Die> {
        let perl = self.perl;
        let hash = self.hash(mapper)?;
        let object_pairs = self.mapper_object_pairs(hash)?;
        let own_pairs = perl.fetch(hash, b"pairs");
        let mut pairs_json = Vec::new();
        let pair_count;
        if perl.ref_is(own_pairs, b"ARRAY") && !perl.items_of(own_pairs)?.is_empty() {
            pair_count = self.list(own_pairs, &mut pairs_json, Self::mapper_pair)?;
        } else if !object_pairs.is_empty() {
            pair_count = object_pairs.len() as u64;
            pairs_json.push(b'[');
            for (index, pair) in object_pairs.iter().enumerate() {
                if index > 0 {
                    pairs_json.push(b',');
                }
                let mut object = Object::open(&mut pairs_json);
                pair.from_end.write_json(object.key("from_end"));
                pair.from_start.write_json(object.key("from_start"));
                pair.ori.write_json(object.key("ori"));
                pair.to_end.write_json(object.key("to_end"));
                pair.to_start.write_json(object.key("to_start"));
                object.close();
            }
            pairs_json.push(b']');
        } else if perl.ref_is(sorted_exons, b"ARRAY") && !perl.items_of(sorted_exons)?.is_empty() {
            pair_count = self.exon_pairs(sorted_exons, strand, &mut pairs_json)?;
        } else {
            pair_count = 0;
            pairs_json.extend_from_slice(b"[]");
        }
        let mut object = Object::open(out);
        self.copy(&mut object, "cdna_coding_end", hash, b"cdna_coding_end");
        self.copy(&mut object, "cdna_coding_start", hash, b"cdna_coding_start");
        write_int(i128::from(pair_count), object.key("pair_count"));
        object.key("pairs").extend_from_slice(&pairs_json);
        self.copy(&mut object, "start_phase", hash, b"start_phase");
        object.close();
        Ok(())
    }

    /// The pairs under `_pair_<from>` of the mapper's
    /// `exon_coord_mapper`, as numbers, in cDNA order.
    fn mapper_object_pairs(&self, mapper: NodeId) -> Result<Vec<ObjectPair>, Die> {
        let perl = self.perl;
        let image = perl.image;
        let Some(coord_mapper) = perl.hash_target(perl.fetch(mapper, b"exon_coord_mapper")) else {
            return Ok(Vec::new());
        };
        let from = perl.fetch(coord_mapper, b"from");
        if !perl.defined(from) {
            return Ok(Vec::new());
        }
        // A `from` that is a reference names no list: its string is an address.
        let Ok(from) = perl.string(from) else {
            return Ok(Vec::new());
        };
        let mut key = b"_pair_".to_vec();
        key.extend_from_slice(&from.bytes);
        let key = OwnedStr {
            bytes: key.into(),
            utf8: from.utf8,
        };
        let Some(lists) = perl.hash_target(fetch_string(image, coord_mapper, &key)) else {
            return Ok(Vec::new());
        };
        let mut pairs = Vec::new();
        for entry in live_entries(image, lists) {
            let Some(list) = perl.array_target(Sv::Node(entry.value)) else {
                continue;
            };
            for item in image.items(list).unwrap_or(&[]) {
                let Some(pair) = perl.hash_target(perl.element(*item)) else {
                    continue;
                };
                let (Some(from), Some(to)) = (
                    perl.hash_target(perl.fetch(pair, b"from")),
                    perl.hash_target(perl.fetch(pair, b"to")),
                ) else {
                    continue;
                };
                let number = |hash: NodeId, key: &[u8]| -> Result<Num, Die> {
                    Ok(perl.number(perl.fetch(hash, key))?.plus(Num::Int(0)))
                };
                let ori = perl.fetch(pair, b"ori");
                pairs.push(ObjectPair {
                    from_start: number(from, b"start")?,
                    from_end: number(from, b"end")?,
                    to_start: number(to, b"start")?,
                    to_end: number(to, b"end")?,
                    ori: if perl.defined(ori) {
                        perl.number(ori)?.plus(Num::Int(0))
                    } else {
                        Num::Int(1)
                    },
                });
            }
        }
        // Perl's sort is a stable merge sort.
        pairs.sort_by(|left, right| {
            left.from_start
                .compare(right.from_start)
                .unwrap_or(Ordering::Equal)
        });
        Ok(pairs)
    }

    /// The pairs of a mapper object that carries none: one an exon in transcript
    /// order, cDNA positions counted from 1, the genomic ends copied as the exon holds
    /// them.
    fn exon_pairs(&mut self, sorted_exons: Sv, strand: Sv, out: &mut Vec<u8>) -> Result<u64, Die> {
        let perl = self.perl;
        let negative = |perl: &Perl<'_>| -> Result<bool, Die> {
            Ok(perl.is_true(strand)
                && perl.number(strand)?.compare(Num::Int(0)) == Some(Ordering::Less))
        };
        let mut exons = perl.items_of(sorted_exons)?;
        if exons.len() > 1 {
            // The sort block reads `$a->{start}` and `$b->{start}` as numbers.
            let reverse = negative(&perl)?;
            let mut keyed = Vec::with_capacity(exons.len());
            for exon in &exons {
                keyed.push((perl.number(perl.get(*exon, b"start")?)?, *exon));
            }
            keyed.sort_by(|left, right| {
                let order = left.0.compare(right.0).unwrap_or(Ordering::Equal);
                if reverse {
                    order.reverse()
                } else {
                    order
                }
            });
            exons = keyed.into_iter().map(|(_, exon)| exon).collect();
        }
        let mut cdna_position = Num::Int(1);
        let mut written = 0;
        out.push(b'[');
        for exon in exons {
            if !perl.ref_is_true(exon) {
                continue;
            }
            let start = perl.get(exon, b"start")?;
            if !perl.is_true(start) {
                continue;
            }
            let end = perl.get(exon, b"end")?;
            if !perl.is_true(end) {
                continue;
            }
            let length = perl
                .number(end)?
                .minus(perl.number(start)?)
                .plus(Num::Int(1));
            let ori = if negative(&perl)? { -1 } else { 1 };
            if written > 0 {
                out.push(b',');
            }
            let mut object = Object::open(out);
            cdna_position
                .plus(length)
                .minus(Num::Int(1))
                .write_json(object.key("from_end"));
            cdna_position.write_json(object.key("from_start"));
            write_int(ori, object.key("ori"));
            self.scalar(end, object.key("to_end"));
            self.scalar(start, object.key("to_start"));
            object.close();
            cdna_position = cdna_position.plus(length);
            written += 1;
        }
        out.push(b']');
        Ok(written)
    }
}

struct ObjectPair {
    from_start: Num,
    from_end: Num,
    to_start: Num,
    to_end: Num,
    ori: Num,
}

/// The bytes of a blob as `MIME::Base64` and `Compress::Zlib` take them: a string
/// under the UTF-8 flag is brought down to one byte a character, and one with a
/// character past 255 is refused.
fn downgrade<'a>(blob: &'a OwnedStr<'_>) -> Result<std::borrow::Cow<'a, [u8]>, Die> {
    if !blob.utf8 || blob.bytes.is_ascii() {
        return Ok(std::borrow::Cow::Borrowed(&blob.bytes));
    }
    let mut bytes = Vec::with_capacity(blob.bytes.len());
    let mut at = 0;
    while at < blob.bytes.len() {
        let byte = blob.bytes[at];
        if byte < 0x80 {
            bytes.push(byte);
            at += 1;
        } else if (byte == 0xc2 || byte == 0xc3) && at + 1 < blob.bytes.len() {
            bytes.push(((byte & 0x03) << 6) | (blob.bytes[at + 1] & 0x3f));
            at += 2;
        } else {
            return Err(Die::new("Wide character in subroutine entry"));
        }
    }
    Ok(std::borrow::Cow::Owned(bytes))
}
