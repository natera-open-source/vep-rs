// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The audit: a seeded sample of transcripts of the Storable cache, the same
//! transcripts in the JSON cache, and every field path of each classed against the
//! converter's rules.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::OsString;
use std::io::Write;
use std::os::unix::ffi::OsStringExt;
use std::path::PathBuf;
use std::sync::Arc;

use crate::fields::{json_fields, json_get, storable_fields, Fields};
use crate::find::{find, name_ends};
use crate::jsonpp::{decode, string_text, write_string, Value as Json};
use crate::number::numify;
use crate::perl::{Die, Perl, Sv};
use crate::perlsort::{shuffled, srand_seed};
use crate::pool::run_ordered;
use crate::shard::{fd_retrieve, gzip_dc_output};
use crate::shards::Shard;
use crate::storable::{NodeId, Value};

/// The exit code of a `die` with no error number behind it.
pub const DIE: i32 = 255;

/// What the command line asks for.
pub struct Request {
    pub storable_root: Vec<u8>,
    pub json_root: Vec<u8>,
    /// The third argument's text: how many transcripts to sample.
    pub sample: Option<Vec<u8>>,
    /// The fourth argument's text: the seed of the files' order.
    pub seed: Option<Vec<u8>>,
    pub threads: usize,
}

fn os_path(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

/// A path the converter writes under another name: the first key of a Storable path
/// and the key it has in the JSON cache.
const RENAME: [(&str, &str); 18] = [
    ("_gene_stable_id", "gene_stable_id"),
    ("_gene_symbol", "gene_symbol"),
    ("_gene_symbol_source", "gene_symbol_source"),
    ("_gene_hgnc_id", "gene_hgnc_id"),
    ("_gene_phenotype", "gene_phenotype"),
    ("_gene_version", "gene_version"),
    ("_ccds", "ccds"),
    ("_refseq", "refseq"),
    ("_protein", "protein"),
    ("_swissprot", "swissprot"),
    ("_trembl", "trembl"),
    ("_uniparc", "uniparc"),
    ("_uniprot_isoform", "uniprot_isoform"),
    ("_source_cache", "source_cache"),
    ("_gene_symbol_id", "gene_symbol_id"),
    ("_bam_edit_status", "bam_edit_status"),
    ("_trans_exon_array", "exons"),
    (
        "_variation_effect_feature_cache",
        "variation_effect_feature_cache",
    ),
];

const CACHE: &str = "_variation_effect_feature_cache";
const JSON_CACHE: &str = "variation_effect_feature_cache";

/// `\w+`: one or more ASCII letters, digits or underscores, the whole of `text`.
fn is_word(text: &[u8]) -> bool {
    !text.is_empty()
        && text
            .iter()
            .all(|byte| byte.is_ascii_alphanumeric() || *byte == b'_')
}

/// Perl's `$`: the end of the text, or the place before one final line feed.
fn at_end(rest: &[u8]) -> bool {
    rest.is_empty() || rest == b"\n"
}

/// `^<head>(?:\.|$)`: the path is `head` or lies below it.
fn is_or_below(path: &[u8], head: &str) -> bool {
    path.strip_prefix(head.as_bytes())
        .is_some_and(|rest| rest.starts_with(b".") || at_end(rest))
}

/// `^<whole>$`.
fn is(path: &[u8], whole: &str) -> bool {
    path.strip_prefix(whole.as_bytes()).is_some_and(at_end)
}

/// `^<head>\.\w+\.<tail>$`: one key of word characters between `head` and `tail`.
fn has_one_key_between(path: &[u8], head: &str, tail: &str) -> bool {
    let Some(rest) = path
        .strip_prefix(head.as_bytes())
        .and_then(|rest| rest.strip_prefix(b"."))
    else {
        return false;
    };
    let rest = rest.strip_suffix(b"\n").unwrap_or(rest);
    let Some(key) = rest
        .strip_suffix(tail.as_bytes())
        .and_then(|key| key.strip_suffix(b"."))
    else {
        return false;
    };
    is_word(key)
}

/// The converter's rule that accounts for a Storable path the JSON cache lacks: the
/// first that fits.
fn dropped_rule(path: &[u8]) -> Option<&'static str> {
    let cache = |tail: &str| format!("{CACHE}.{tail}");
    if is_or_below(path, "slice") {
        return Some("slice: the chromosome is the shard directory and no slice field is read");
    }
    if is_or_below(path, "_gene") {
        return Some("gene object: its stable_id and version are carried by the _gene_* renames; start, end and strand are read by nothing at run time");
    }
    if is(path, "_vep_lazy_loaded") {
        return Some("runtime flag set by the Perl cache loader");
    }
    if is_or_below(path, "translation.start_exon") || is_or_below(path, "translation.end_exon") {
        return Some("translation exons: derivable from exons and the translation start and end");
    }
    if is(path, "translation.transcript") {
        return Some("back-reference to the enclosing transcript");
    }
    if is_or_below(path, &cache("mapper.exon_coord_mapper")) {
        return Some("exon coordinate mapper: replaced by mapper.pairs");
    }
    if is(path, "_trans_exon_array[].strand") || is(path, &cache("sorted_exons[].strand")) {
        return Some("exon strand: the transcript strand");
    }
    if is(path, &cache("introns[].strand")) {
        return Some("intron strand: the transcript strand");
    }
    if is_or_below(path, &cache("introns[].prev")) || is_or_below(path, &cache("introns[].next")) {
        return Some("intron exon references: the flanking exons are the exons array");
    }
    if is_or_below(path, &cache("introns[].slice")) {
        return Some(
            "intron slice: the chromosome is the shard directory and no slice field is read",
        );
    }
    if is(path, &cache("seq_edits[].description")) {
        return Some("sequence edit description: read by nothing at run time");
    }
    if has_one_key_between(
        path,
        &cache("protein_function_predictions"),
        "matrix_compressed",
    ) {
        return Some(
            "matrix_compressed: constant, every exported predictions_data is the gzip blob",
        );
    }
    if path.starts_with(cache("five_prime_utr.").as_bytes())
        || path.starts_with(cache("three_prime_utr.").as_bytes())
    {
        return Some("UTR object flattened to its sequence string");
    }
    None
}

/// The converter's rule that produces a JSON path no Storable path maps to.
fn derived_rule(path: &[u8]) -> Option<&'static str> {
    let cache = |tail: &str| format!("{JSON_CACHE}.{tail}");
    let pairs = cache("mapper.pairs");
    if path
        .strip_prefix(pairs.as_bytes())
        .is_some_and(|rest| rest.starts_with(b"[]") || rest.starts_with(b".") || at_end(rest))
    {
        return Some("pair list taken from mapper.pairs or computed from sorted_exons");
    }
    if is(path, &cache("mapper.pair_count")) {
        return Some("length of mapper.pairs");
    }
    if is(path, &cache("seq_edits")) {
        return Some("seq_edits: written as [] when the transcript has no translation");
    }
    if is(path, &cache("seq_edits[].length_diff")) {
        return Some("seq_edits length_diff: computed from start, end and alt_seq as SeqEdit::length_diff does");
    }
    None
}

/// The JSON path a Storable path is written under: its first key renamed, a
/// prediction's `matrix` as `predictions_data`. The path is split at its dots and
/// joined again, which drops dots at its end.
fn mapped_path(path: &[u8]) -> Vec<u8> {
    let mut parts: Vec<Vec<u8>> = path
        .split(|byte| *byte == b'.')
        .map(<[u8]>::to_vec)
        .collect();
    while parts.last().is_some_and(Vec::is_empty) {
        parts.pop();
    }
    if let Some(first) = parts.first_mut() {
        // The first key without a final `[]`. Perl's `$` also stands before one final
        // line feed, which a renamed key loses; a key with a line feed inside it is
        // not looked up.
        let body = first.strip_suffix(b"\n").unwrap_or(first);
        if !body.contains(&b'\n') {
            let (head, suffix) = match body.strip_suffix(b"[]") {
                Some(head) => (head, &b"[]"[..]),
                None => (body, &b""[..]),
            };
            if let Some((_, to)) = RENAME.iter().find(|(from, _)| from.as_bytes() == head) {
                *first = [to.as_bytes(), suffix].concat();
            }
        }
    }
    let mut mapped = parts.join(&b'.');
    let newline = mapped.ends_with(b"\n");
    let body_len = mapped.len() - usize::from(newline);
    if has_one_key_between(
        &mapped[..body_len],
        &format!("{JSON_CACHE}.protein_function_predictions"),
        "matrix",
    ) {
        mapped.truncate(body_len - "matrix".len());
        mapped.extend_from_slice(b"predictions_data");
        if newline {
            mapped.push(b'\n');
        }
    }
    mapped
}

/// `my ($parent) = ($k =~ /^(.*)\.[^.]+$/)`: the path above the last key.
fn parent(path: &[u8]) -> Option<&[u8]> {
    let dot = path.iter().rposition(|byte| *byte == b'.')?;
    (dot + 1 < path.len() && !path[..dot].contains(&b'\n')).then(|| &path[..dot])
}

/// Counts by name, written in the order of the names.
type Counts = BTreeMap<Vec<u8>, u64>;

fn count(counts: &mut Counts, name: Vec<u8>) {
    *counts.entry(name).or_insert(0) += 1;
}

fn write_counts(counts: &Counts, out: &mut Vec<u8>) {
    out.push(b'{');
    for (index, (name, number)) in counts.iter().enumerate() {
        if index > 0 {
            out.push(b',');
        }
        write_string(name, out);
        out.extend_from_slice(format!(":{number}").as_bytes());
    }
    out.push(b'}');
}

/// One sampled transcript of the Storable cache.
struct Sampled {
    /// The directory the shard stands in, or `unknown`.
    contig: Vec<u8>,
    /// The stable id, where it is a true value.
    id: Option<SampledId>,
    fields: Fields,
}

/// The stable id of a sampled transcript.
#[derive(Clone)]
enum SampledId {
    /// The id's characters, and whether the cache stores it as a number: `JSON::PP`
    /// writes such an id without quotes.
    Text { characters: Vec<u8>, number: bool },
    /// A reference: Perl reads it as its address, which names no transcript.
    Address,
}

/// One Storable file of the sample's order: its transcripts, where it holds any.
struct Chunk {
    shard: Arc<Shard>,
    transcripts: Vec<NodeId>,
}

/// The directory of a shard: `m{/(\w+)/[^/]+\.gz$}` on its name. A name without a
/// directory of word characters between two slashes before it has none.
fn contig_of(file: &[u8]) -> Vec<u8> {
    let found = (|| {
        let slash = file.iter().rposition(|byte| *byte == b'/')?;
        let name = &file[slash + 1..];
        let name = name.strip_suffix(b"\n").unwrap_or(name);
        if name.len() < 4 || !name.ends_with(b".gz") {
            return None;
        }
        let directory = &file[..slash];
        let from = directory.iter().rposition(|byte| *byte == b'/')? + 1;
        is_word(&directory[from..]).then(|| directory[from..].to_vec())
    })();
    found.unwrap_or_else(|| b"unknown".to_vec())
}

/// The transcripts of one Storable file as the sample takes them: every item of a list
/// that is a reference to a hash with a `stable_id` key, and a hash that has that key
/// itself. `Ok(None)` is a file that is no shard with an unblessed hash at its root.
fn read_chunk(file: &[u8]) -> Result<Option<Chunk>, Die> {
    // A file that cannot be read gives gunzip nothing to write.
    let compressed = std::fs::read(os_path(file)).unwrap_or_default();
    let Ok(Some(image)) = fd_retrieve(gzip_dc_output(&compressed)) else {
        return Ok(None);
    };
    let root = image.root();
    let plain_hash = matches!(image.value(root), Value::Hash { .. })
        && matches!(image.class(root), None | Some(b"HASH"));
    if !plain_hash {
        return Ok(None);
    }
    let shard = Arc::new(Shard::new(image));
    let image = &shard.image;
    let perl = Perl::new(image);
    let has_id = |sv: Sv| -> Result<Option<NodeId>, Die> {
        Ok(perl
            .hash_of(sv)?
            .filter(|hash| image.fetch(*hash, b"stable_id").is_some()))
    };
    let mut transcripts = Vec::new();
    for (_, value) in crate::dumper::hash_entries_in_stored_order(image, root) {
        let entry = shard.read(Sv::Node(value));
        if perl.ref_is(entry, b"ARRAY") {
            for item in perl.items_of(entry)? {
                let item = shard.read(item);
                // `ref($_) && exists($_->{stable_id})`.
                if perl.ref_is_true(item) {
                    transcripts.extend(has_id(item)?);
                }
            }
        } else if perl.ref_is(entry, b"HASH") {
            transcripts.extend(has_id(entry)?);
        }
    }
    Ok(Some(Chunk { shard, transcripts }))
}

/// What one JSON file says: the id of every transcript it holds, and the fields of the
/// sampled ones, in the file's order.
#[derive(Default)]
struct JsonFile {
    ids: Vec<IdKey>,
    sampled: Vec<(Vec<u8>, Fields)>,
}

/// The key a transcript is indexed under.
enum IdKey {
    Text(Vec<u8>),
    /// A reference: Perl indexes it under its address, which no other transcript shares.
    Address,
}

/// The key of a transcript's `stable_id` where the value is true: its text.
fn id_key(text: &[u8], value: &Json) -> Option<IdKey> {
    match value {
        Json::Null | Json::Bool(false) => None,
        Json::Bool(true) => Some(IdKey::Text(b"1".to_vec())),
        Json::Str(from, to) => {
            let id = string_text(&text[*from..*to]);
            (!id.is_empty() && id != b"0").then_some(IdKey::Text(id))
        }
        Json::Number(from, to) => {
            let raw = &text[*from..*to];
            // An integer of twenty characters or more stays the text it is.
            let integer = !raw.iter().any(|byte| matches!(byte, b'.' | b'e' | b'E'));
            if integer && raw.len() > 19 {
                return Some(IdKey::Text(raw.to_vec()));
            }
            // A number that is zero, with whatever sign and exponent, is a false value.
            let number = numify(raw);
            (number.as_f64() != 0.0).then(|| IdKey::Text(number.to_perl_string()))
        }
        Json::List(_) | Json::Dict(_) => Some(IdKey::Address),
    }
}

fn scan_json(file: &[u8], wanted: &HashSet<Vec<u8>>) -> JsonFile {
    let mut scanned = JsonFile::default();
    let Ok(text) = std::fs::read(os_path(file)) else {
        return scanned;
    };
    let Ok(data) = decode(&text) else {
        return scanned;
    };
    let mut list = |items: &[Json]| {
        for item in items {
            let Json::Dict(entries) = item else {
                continue;
            };
            let Some(key) = json_get(&text, entries, b"stable_id").and_then(|id| id_key(&text, id))
            else {
                continue;
            };
            if let IdKey::Text(id) = &key {
                if wanted.contains(id) {
                    scanned
                        .sampled
                        .push((id.clone(), json_fields(&text, entries)));
                }
            }
            scanned.ids.push(key);
        }
    };
    match &data {
        Json::List(items) => list(items),
        Json::Dict(entries) => {
            // The lists of an object, each key once with its last value.
            let mut last: HashMap<Vec<u8>, usize> = HashMap::new();
            let mut order = Vec::new();
            for (index, ((from, to), _)) in entries.iter().enumerate() {
                let key = string_text(&text[*from..*to]);
                if last.insert(key.clone(), index).is_none() {
                    order.push(key);
                }
            }
            for key in order {
                if let Json::List(items) = &entries[last[&key]].1 {
                    list(items);
                }
            }
        }
        _ => {}
    }
    scanned
}

/// Runs the audit. `stdout` takes the manifest and `stderr` the lines of the walk and
/// of a run that died. Returns the exit code.
pub fn run(request: &Request, stdout: &mut dyn Write, stderr: &mut dyn Write) -> i32 {
    let seed_text = request.seed.clone().unwrap_or_else(|| b"1".to_vec());
    let sample_text = request.sample.clone().unwrap_or_else(|| b"1000".to_vec());
    let mut warn = |lines: Vec<Vec<u8>>| {
        for mut line in lines {
            line.push(b'\n');
            stderr.write_all(&line).ok();
        }
    };

    let storable = find(&request.storable_root, |name| {
        name_ends(name, b".gz") && !name_ends(name, b"_reg.gz") && !name_ends(name, b"_var.gz")
    });
    warn(storable.warnings);
    if storable.names.is_empty() {
        let mut line = b"No .gz Storable chunks found under ".to_vec();
        line.extend_from_slice(&request.storable_root);
        warn(vec![line]);
        return DIE;
    }
    let files: Vec<&Vec<u8>> = shuffled(storable.names.len(), srand_seed(&seed_text))
        .into_iter()
        .map(|place| &storable.names[place])
        .collect();

    // The quota is filled in the order of the files; the files are read ahead of it on
    // the worker threads.
    let mut remaining = numify(&sample_text).as_f64();
    let mut sampled: Vec<Sampled> = Vec::new();
    let died: Option<Die> = run_ordered(
        &files,
        request.threads,
        |file| read_chunk(file),
        |pool| {
            for (index, file) in files.iter().enumerate() {
                if remaining <= 0.0 {
                    break;
                }
                let chunk = match pool.take(index) {
                    Ok(Some(chunk)) => chunk,
                    Ok(None) => continue,
                    Err(die) => return Some(die),
                };
                let perl = Perl::new(&chunk.shard.image);
                for transcript in chunk.transcripts {
                    if remaining <= 0.0 {
                        break;
                    }
                    remaining -= 1.0;
                    // `my $sid = $st->{stable_id} or next`.
                    let id = chunk.shard.read(perl.fetch(transcript, b"stable_id"));
                    let id = if !perl.is_true(id) {
                        None
                    } else if perl.target(id).is_some() {
                        Some(SampledId::Address)
                    } else {
                        let number = matches!(
                            perl.value(id),
                            Some(Value::Integer(_)) | Some(Value::Double(_))
                        );
                        match perl.string(id) {
                            Ok(text) => Some(SampledId::Text {
                                characters: text.characters(),
                                number,
                            }),
                            Err(die) => return Some(die),
                        }
                    };
                    let fields = match storable_fields(&chunk.shard, transcript) {
                        Ok(fields) => fields,
                        Err(_) => {
                            return Some(Die::new(
                                "a transcript holds an array that holds itself: Perl walks it without end",
                            ))
                        }
                    };
                    sampled.push(Sampled {
                        contig: contig_of(file),
                        id,
                        fields,
                    });
                }
            }
            None
        },
    );
    if let Some(die) = died {
        warn(vec![die.0]);
        return DIE;
    }

    let wanted: HashSet<Vec<u8>> = sampled
        .iter()
        .filter_map(|entry| match &entry.id {
            Some(SampledId::Text { characters, .. }) => Some(characters.clone()),
            _ => None,
        })
        .collect();
    let json = find(&request.json_root, |name| name_ends(name, b".json"));
    warn(json.warnings);
    let mut json_files = json.names;
    json_files.sort();
    let mut indexed: HashSet<Vec<u8>> = HashSet::new();
    let mut addresses = 0usize;
    let mut json_by_id: HashMap<Vec<u8>, Fields> = HashMap::new();
    run_ordered(
        &json_files,
        request.threads,
        |file| scan_json(file, &wanted),
        |pool| {
            for index in 0..json_files.len() {
                let scanned = pool.take(index);
                for key in scanned.ids {
                    match key {
                        IdKey::Text(id) => {
                            indexed.insert(id);
                        }
                        IdKey::Address => addresses += 1,
                    }
                }
                // A later transcript of an id replaces an earlier one.
                json_by_id.extend(scanned.sampled);
            }
        },
    );

    let (mut identical, mut renamed, mut converted) = (Counts::new(), Counts::new(), Counts::new());
    let (mut dropped_undef, mut dropped, mut derived) =
        (Counts::new(), Counts::new(), Counts::new());
    let mut gaps = Counts::new();
    // Up to three ids a gap, each with whether it is written as a number.
    let mut gap_examples: BTreeMap<Vec<u8>, Vec<(Vec<u8>, bool)>> = BTreeMap::new();
    let (mut matched, mut unmatched_canonical, mut unmatched_other) = (0u64, 0u64, 0u64);
    for entry in &sampled {
        let Some(id) = &entry.id else {
            continue;
        };
        let found = match id {
            SampledId::Text { characters, number } => json_by_id
                .get(characters)
                .map(|fields| (fields, characters, *number)),
            SampledId::Address => None,
        };
        let Some((json_fields, id, id_is_number)) = found else {
            let digits = !entry.contig.is_empty() && entry.contig.iter().all(u8::is_ascii_digit);
            if digits || matches!(entry.contig.as_slice(), b"X" | b"Y" | b"MT") {
                unmatched_canonical += 1;
            } else {
                unmatched_other += 1;
            }
            continue;
        };
        matched += 1;
        let mut note_gap = |gap: &[u8]| {
            let examples = gap_examples.entry(gap.to_vec()).or_default();
            if examples.len() < 3 {
                examples.push((id.clone(), id_is_number));
            }
        };
        // Shortest paths first, so that the fields of a container the converter wrote
        // as one string are passed over with it.
        let mut paths: Vec<&Vec<u8>> = entry.fields.keys().collect();
        paths.sort_by_key(|path| (String::from_utf8_lossy(path).chars().count(), *path));
        let mut accounted: HashSet<Vec<u8>> = HashSet::new();
        let mut flattened: HashSet<&[u8]> = HashSet::new();
        for path in paths {
            if parent(path).is_some_and(|parent| flattened.contains(parent)) {
                continue;
            }
            let kind = &entry.fields[path];
            let mapped = mapped_path(path);
            if let Some(json_kind) = json_fields.get(&mapped) {
                if kind == json_kind {
                    count(
                        if mapped == *path {
                            &mut identical
                        } else {
                            &mut renamed
                        },
                        path.clone(),
                    );
                } else {
                    count(
                        &mut converted,
                        [path.as_slice(), b":", kind, b"->", json_kind].concat(),
                    );
                    let container =
                        |name: &[u8]| json_kind.windows(name.len()).any(|part| part == name);
                    if !container(b"HASH") && !container(b"ARRAY") {
                        flattened.insert(path.as_slice());
                    }
                }
                accounted.insert(mapped);
                continue;
            }
            if kind == b"undef" {
                count(&mut dropped_undef, path.clone());
            } else if let Some(rule) = dropped_rule(path) {
                count(
                    &mut dropped,
                    [path.as_slice(), b": ", rule.as_bytes()].concat(),
                );
            } else {
                let gap = [b"missing_in_json:", path.as_slice()].concat();
                note_gap(&gap);
                count(&mut gaps, gap);
            }
        }
        for path in json_fields.keys() {
            if accounted.contains(path) {
                continue;
            }
            if let Some(rule) = derived_rule(path) {
                count(
                    &mut derived,
                    [path.as_slice(), b": ", rule.as_bytes()].concat(),
                );
            } else {
                let gap = [b"only_in_json:", path.as_slice()].concat();
                note_gap(&gap);
                count(&mut gaps, gap);
            }
        }
    }

    // The manifest: one object, its keys and the keys of every object in it in the
    // order of their names.
    let mut out = Vec::new();
    let key = |out: &mut Vec<u8>, name: &str, first: bool| {
        out.extend_from_slice(if first { b"{" } else { b"," });
        write_string(name.as_bytes(), out);
        out.push(b':');
    };
    key(&mut out, "conversion_rules", true);
    write_counts(&converted, &mut out);
    key(&mut out, "derived_on_json_side", false);
    write_counts(&derived, &mut out);
    key(&mut out, "dropped_by_converter", false);
    write_counts(&dropped, &mut out);
    key(&mut out, "dropped_undef_counts", false);
    write_counts(&dropped_undef, &mut out);
    key(&mut out, "identical_field_counts", false);
    write_counts(&identical, &mut out);
    key(&mut out, "json_transcripts_indexed", false);
    out.extend_from_slice((indexed.len() + addresses).to_string().as_bytes());
    key(&mut out, "matched_in_json", false);
    out.extend_from_slice(matched.to_string().as_bytes());
    key(&mut out, "parity_gap_examples", false);
    out.push(b'{');
    for (index, (gap, ids)) in gap_examples.iter().enumerate() {
        if index > 0 {
            out.push(b',');
        }
        write_string(gap, &mut out);
        out.extend_from_slice(b":[");
        for (place, (id, number)) in ids.iter().enumerate() {
            if place > 0 {
                out.push(b',');
            }
            if *number {
                out.extend_from_slice(id);
            } else {
                write_string(id, &mut out);
            }
        }
        out.push(b']');
    }
    out.push(b'}');
    key(&mut out, "parity_gaps", false);
    write_counts(&gaps, &mut out);
    key(&mut out, "renamed_field_counts", false);
    write_counts(&renamed, &mut out);
    key(&mut out, "sample_seed", false);
    match &request.seed {
        // The seed given is a string; the seed of a run without one is the number 1.
        Some(seed) => write_string(seed, &mut out),
        None => out.push(b'1'),
    }
    key(&mut out, "sampled_transcripts", false);
    out.extend_from_slice(sampled.len().to_string().as_bytes());
    key(&mut out, "storable_files_scanned", false);
    out.extend_from_slice(storable.names.len().to_string().as_bytes());
    key(&mut out, "unmatched_in_json_by_contig_class", false);
    out.extend_from_slice(
        format!("{{\"canonical_contig\":{unmatched_canonical},\"non_canonical_contig\":{unmatched_other}}}")
            .as_bytes(),
    );
    out.push(b'}');
    if stdout
        .write_all(&out)
        .and_then(|()| stdout.flush())
        .is_err()
    {
        return DIE;
    }
    0
}

#[cfg(test)]
#[path = "audit_tests.rs"]
mod tests;
