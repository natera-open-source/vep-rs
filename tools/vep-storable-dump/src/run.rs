// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run: the ids asked for, every transcript shard of the canonical chromosomes in
//! the order 1 to 22, X, Y, MT, and one block a wanted transcript the first time its id
//! is met.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::io::Write;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

use crate::dumper::{dump, hash_entries, Indent, Key, Overlay, VNode, Val};
use crate::json::PerlStr;
use crate::perl::{Die, Perl, Sv};
use crate::pool::run_ordered;
use crate::shards::{
    join, read_shard, transcript_shards, transcripts, Shard, ShardError, CANONICAL,
};
use crate::storable::NodeId;
use crate::text::{strip_version, Text};

/// The words that open the summary line.
const TAG: &str = "[dump_storable_transcripts]";

/// The exit code of a `die` with no error number behind it.
pub const DIE: i32 = 255;

/// What the command line asks for.
pub struct Request {
    pub cache_dir: Vec<u8>,
    pub ids: Vec<u8>,
    pub out: Vec<u8>,
    pub threads: usize,
}

fn path_of(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

/// The text of an operating system error without Rust's ` (os error <n>)`.
fn os_text(error: &std::io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(at) if text.ends_with(')') => text[..at].to_string(),
        _ => text,
    }
}

/// The exit code of a Perl `die`: the error number when there is one, else 255.
fn die_code(error: &std::io::Error) -> i32 {
    match error.raw_os_error() {
        Some(code) if code & 0xff != 0 => code & 0xff,
        _ => DIE,
    }
}

/// The ids of the list: every line that starts with `ENS`, without its line feed and
/// without a final `.<digits>`. The list is read as bytes, one character a byte; an id
/// is held as its characters in UTF-8, the form ids met in the cache are compared in.
pub fn wanted_ids(list: &[u8]) -> HashSet<Vec<u8>> {
    let mut ids = HashSet::new();
    for line in list.split_inclusive(|byte| *byte == b'\n') {
        let mut id = line.strip_suffix(b"\n").unwrap_or(line).to_vec();
        if !id.starts_with(b"ENS") {
            continue;
        }
        strip_version(&mut id);
        ids.insert(PerlStr::bytes(&id).characters().into_owned());
    }
    ids
}

/// A transcript of a shard that the list names: its id and its block.
pub struct Block {
    /// The id's characters: what the list of ids and the ids met are compared by.
    pub id: Vec<u8>,
    pub bytes: Vec<u8>,
}

/// The transcripts of one shard that the list names, and what ended the shard early if
/// anything did.
#[derive(Default)]
pub struct ShardDump {
    pub blocks: Vec<Block>,
    pub failure: Option<Die>,
}

/// `ref $x eq 'HASH'`: a reference to an unblessed hash, or to a value blessed into a
/// class of that name.
fn plain_hash(shard: &Shard, sv: Sv) -> Option<NodeId> {
    let perl = Perl::new(&shard.image);
    if !perl.ref_is(shard.read(sv), b"HASH") {
        return None;
    }
    perl.target(sv)
}

/// The scalar a `%s` or a string operator makes of `sv`.
fn string_of<'a>(perl: &Perl<'a>, sv: Sv) -> Result<(Cow<'a, [u8]>, bool), Die> {
    let text = perl.string(sv)?;
    Ok((text.bytes, text.utf8))
}

/// The summary that stands in the place of one prediction matrix:
/// `matrix:<length> md5:<translation_md5>`, or `none` for a value that is no reference.
fn matrix_summary<'i>(shard: &'i Shard, value: NodeId) -> Result<Val<'i>, Die> {
    let perl = Perl::new(&shard.image);
    let matrix = shard.read(Sv::Node(value));
    if !perl.ref_is_true(matrix) {
        return Ok(Val::Str {
            bytes: Cow::Borrowed(b"none"),
            utf8: false,
        });
    }
    let held = shard.read(perl.get(matrix, b"matrix")?);
    let mut text = Text::bytes(&b"matrix:"[..]);
    if perl.defined(held) {
        // `length`: the characters of the scalar's text.
        let (bytes, utf8) = string_of(&perl, held)?;
        let length = PerlStr {
            bytes: &bytes,
            utf8,
        }
        .char_len();
        text.push(length.to_string().as_bytes(), false);
    } else {
        text.push(b"undef", false);
    }
    text.push(b" md5:", false);
    let md5 = shard.read(perl.get(matrix, b"translation_md5")?);
    let (bytes, utf8) = string_of(&perl, md5)?;
    text.push(&bytes, utf8);
    Ok(Val::Str {
        bytes: Cow::Owned(text.bytes),
        utf8: text.utf8,
    })
}

/// The copy of a transcript that is dumped: an unblessed hash of the transcript's own
/// values, with a copy of `_variation_effect_feature_cache` in which every prediction
/// matrix is its summary.
fn dumped_copy<'i>(
    shard: &'i Shard,
    transcript: NodeId,
    overlay: &mut Overlay<'i>,
) -> Result<Val<'i>, Die> {
    let image = &shard.image;
    let mut top: Vec<(Key<'i>, Val<'i>)> = Vec::new();
    for (key, value) in hash_entries(image, transcript) {
        let is_cache = key.characters().as_ref() == b"_variation_effect_feature_cache";
        let cache = is_cache
            .then(|| plain_hash(shard, Sv::Node(value)))
            .flatten();
        let Some(cache) = cache else {
            top.push((key, Val::CopyOf(value)));
            continue;
        };
        let mut copy: Vec<(Key<'i>, Val<'i>)> = Vec::new();
        for (inner_key, inner) in hash_entries(image, cache) {
            let is_predictions = inner_key.characters().as_ref() == b"protein_function_predictions";
            let predictions = is_predictions
                .then(|| plain_hash(shard, Sv::Node(inner)))
                .flatten();
            let Some(predictions) = predictions else {
                copy.push((inner_key, Val::CopyOf(inner)));
                continue;
            };
            let mut summaries: Vec<(Key<'i>, Val<'i>)> = Vec::new();
            for (name, matrix) in hash_entries(image, predictions) {
                summaries.push((name, matrix_summary(shard, matrix)?));
            }
            copy.push((inner_key, Val::Ref(overlay.push(VNode::Hash(summaries)))));
        }
        top.push((key, Val::Ref(overlay.push(VNode::Hash(copy)))));
    }
    Ok(Val::Ref(overlay.push(VNode::Hash(top))))
}

/// The blocks of the wanted transcripts of one shard, in the shard's order; an id met
/// twice in the shard is dumped once.
pub fn dump_shard(shard: Option<&Shard>, wanted: &HashSet<Vec<u8>>) -> ShardDump {
    let mut dumped = ShardDump::default();
    dumped.failure = dump_into(shard, wanted, &mut dumped.blocks).err();
    dumped
}

fn dump_into(
    shard: Option<&Shard>,
    wanted: &HashSet<Vec<u8>>,
    blocks: &mut Vec<Block>,
) -> Result<(), Die> {
    let list = transcripts(shard)?;
    let Some(shard) = shard else {
        return Ok(());
    };
    let perl = Perl::new(&shard.image);
    let mut met: HashSet<Vec<u8>> = HashSet::new();
    for transcript in list {
        // `$t->{stable_id} // ''`.
        let id = shard.read(perl.get(transcript, b"stable_id")?);
        if perl.target(id).is_some() {
            // A reference reads as its address in memory, which names no transcript.
            continue;
        }
        let id = perl.string(id)?;
        let mut id = Text {
            bytes: id.bytes.into_owned(),
            utf8: id.utf8,
        };
        strip_version(&mut id.bytes);
        let characters = id.characters();
        if !wanted.contains(&characters) || !met.insert(characters.clone()) {
            continue;
        }
        let Some(node) = perl.hash_of(transcript)? else {
            continue;
        };
        let mut overlay = Overlay::default();
        let copy = dumped_copy(shard, node, &mut overlay)?;
        let mut head = Text::bytes(&b"### "[..]);
        head.push(&id.bytes, id.utf8);
        head.push(b"\n", false);
        let mut bytes = head.printed();
        dump(shard, &overlay, &copy, Indent::Fixed, &mut bytes)
            .map_err(|error| Die::new(error.to_string()))?;
        bytes.push(b'\n');
        blocks.push(Block {
            id: characters,
            bytes,
        });
    }
    Ok(())
}

/// The shards of the run in its order: the chromosomes 1 to 22, X, Y, MT, a
/// chromosome's shards by name.
fn tasks(cache_dir: &Path) -> Vec<PathBuf> {
    CANONICAL
        .iter()
        .flat_map(|contig| transcript_shards(&join(cache_dir, contig.as_bytes())))
        .collect()
}

fn die_line(stderr: &mut dyn Write, mut line: Vec<u8>) {
    line.push(b'\n');
    stderr.write_all(&line).ok();
}

/// Runs the dump. Returns the exit code; `stderr` takes the summary line and the line of
/// a run that died.
pub fn run(request: &Request, stderr: &mut dyn Write) -> i32 {
    let list = match std::fs::read(path_of(&request.ids)) {
        Ok(list) => list,
        Err(error) => {
            let mut line = b"cannot read ".to_vec();
            line.extend_from_slice(&request.ids);
            line.extend_from_slice(format!(": {}", os_text(&error)).as_bytes());
            die_line(stderr, line);
            return die_code(&error);
        }
    };
    let wanted = wanted_ids(&list);
    let cannot_write = |error: &std::io::Error| {
        let mut line = b"cannot write ".to_vec();
        line.extend_from_slice(&request.out);
        line.extend_from_slice(format!(": {}", os_text(error)).as_bytes());
        line
    };
    let file = match std::fs::File::create(path_of(&request.out)) {
        Ok(file) => file,
        Err(error) => {
            die_line(stderr, cannot_write(&error));
            return die_code(&error);
        }
    };
    let mut out = std::io::BufWriter::with_capacity(1 << 16, file);
    let tasks = tasks(&path_of(&request.cache_dir));
    let work = |path: &PathBuf| -> Result<ShardDump, ShardError> {
        let shard = read_shard(path)?;
        Ok(dump_shard(shard.as_ref(), &wanted))
    };
    let mut seen: HashMap<Vec<u8>, ()> = HashMap::new();
    let died: Option<Vec<u8>> = run_ordered(&tasks, request.threads, work, |pool| {
        for (index, path) in tasks.iter().enumerate() {
            let dumped = match pool.take(index) {
                Ok(dumped) => dumped,
                Err(ShardError::Gunzip(_)) => {
                    let mut line = b"gunzip ".to_vec();
                    line.extend_from_slice(path.as_os_str().as_bytes());
                    line.extend_from_slice(b" failed");
                    return Some(line);
                }
                Err(ShardError::Retrieve(die)) => return Some(die.0),
            };
            for block in dumped.blocks {
                if seen.insert(block.id, ()).is_some() {
                    continue;
                }
                if let Err(error) = out.write_all(&block.bytes) {
                    return Some(cannot_write(&error));
                }
            }
            if let Some(die) = dumped.failure {
                return Some(die.0);
            }
        }
        None
    });
    let flushed = out.flush();
    if let Some(line) = died {
        die_line(stderr, line);
        return DIE;
    }
    if let Err(error) = flushed {
        die_line(stderr, cannot_write(&error));
        return DIE;
    }
    let mut line = format!("{TAG} ").into_bytes();
    line.extend_from_slice(&request.cache_dir);
    line.extend_from_slice(
        format!(
            ": {} of {} requested transcripts dumped\n",
            seen.len(),
            wanted.len()
        )
        .as_bytes(),
    );
    stderr.write_all(&line).ok();
    0
}

#[cfg(test)]
mod tests {
    use super::wanted_ids;

    #[test]
    fn the_list_names_lines_that_start_with_ens_without_their_version() {
        let ids = wanted_ids(
            b"transcript_id\nENST00000000001\nENST00000000002.7\n ENST00000000003\nENSMUST00000000004.1\r\n\nENST00000000005\tchr\nENST00000000001.2\nens\nENS",
        );
        let mut ids: Vec<String> = ids
            .into_iter()
            .map(|id| String::from_utf8(id).expect("ASCII"))
            .collect();
        ids.sort();
        assert_eq!(
            ids,
            [
                "ENS",
                "ENSMUST00000000004.1\r",
                "ENST00000000001",
                "ENST00000000002",
                "ENST00000000005\tchr",
            ]
        );
    }

    #[test]
    fn a_byte_above_127_of_the_list_is_one_character() {
        let ids = wanted_ids(b"ENST\xe9\n");
        assert!(ids.contains("ENST\u{e9}".as_bytes()));
    }
}
