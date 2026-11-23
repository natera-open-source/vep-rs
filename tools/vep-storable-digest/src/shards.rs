// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript shards of an Ensembl VEP Storable cache: the files of a chromosome
//! directory, one file as `IO::Uncompress::Gunzip`'s `gunzip` and Storable's
//! `fd_retrieve` give it, and the list of transcripts a shard holds.

use std::collections::HashSet;
use std::ffi::OsString;
use std::io::Read;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

use flate2::read::DeflateDecoder;

use crate::perl::{Die, Perl, Sv};
use crate::shard::fd_retrieve;
use crate::storable::{Image, NodeId, Value};

/// The chromosomes whose shards are read: 1 to 22, X, Y and MT.
pub const CANONICAL: [&str; 25] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// A retrieved shard: its image, and what Perl's reference counts leave of it.
pub struct Shard {
    pub image: Image,
    /// The weak references whose target nothing else holds: Perl frees the target when
    /// the image is read, and the reference reads as undef.
    cleared: HashSet<NodeId>,
}

impl Shard {
    pub fn new(image: Image) -> Self {
        let cleared = cleared_weak_references(&image);
        Shard { image, cleared }
    }

    /// Whether `id` is a weak reference that reads as undef.
    pub fn cleared(&self, id: NodeId) -> bool {
        !self.cleared.is_empty() && self.cleared.contains(&id)
    }

    /// The scalar a program reads at `sv`: undef where it is such a weak reference.
    pub fn read(&self, sv: Sv) -> Sv {
        match sv {
            Sv::Node(id) if self.cleared(id) => Sv::Missing,
            sv => sv,
        }
    }
}

/// The weak references of an image whose target is freed once the image is read.
/// Perl counts the holders of a value: a container that holds it, a reference to it
/// that is not weak, and the reference `retrieve` returns for the root. A value with
/// no holder is freed, and each value it held has one holder fewer.
fn cleared_weak_references(image: &Image) -> HashSet<NodeId> {
    let mut weak = Vec::new();
    for id in image.node_ids() {
        if let Value::Ref { weak: true, .. } = image.value(id) {
            weak.push(id);
        }
    }
    if weak.is_empty() {
        return HashSet::new();
    }
    let index = |id: NodeId| id.tag() as usize;
    let mut holders = vec![0u32; image.node_count()];
    holders[index(image.root())] += 1;
    let held = |id: NodeId| -> Vec<NodeId> {
        match image.value(id) {
            Value::Array(_) => image.items(id).unwrap_or(&[]).to_vec(),
            Value::Hash { .. } => image
                .entries(id)
                .unwrap_or(&[])
                .iter()
                .map(|entry| entry.value)
                .collect(),
            Value::Ref {
                target,
                weak: false,
                ..
            } => vec![target],
            _ => Vec::new(),
        }
    };
    for id in image.node_ids() {
        for value in held(id) {
            holders[index(value)] += 1;
        }
    }
    let mut freed = vec![false; image.node_count()];
    let mut unheld: Vec<NodeId> = image
        .node_ids()
        .filter(|id| holders[index(*id)] == 0)
        .collect();
    while let Some(id) = unheld.pop() {
        if std::mem::replace(&mut freed[index(id)], true) {
            continue;
        }
        for value in held(id) {
            holders[index(value)] -= 1;
            if holders[index(value)] == 0 {
                unheld.push(value);
            }
        }
    }
    weak.into_iter()
        .filter(|id| match image.value(*id) {
            Value::Ref { target, .. } => freed[index(target)],
            _ => false,
        })
        .collect()
}

/// Why a shard could not be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShardError {
    /// `gunzip` returned false; the text says why.
    Gunzip(String),
    /// Storable died with this text.
    Retrieve(Die),
}

/// The files `sort glob("$dir/*.gz")` gives whose name is `<digits>-<digits>.gz`: the
/// transcript shards of one chromosome directory, in byte order. A directory that is
/// absent or cannot be listed has none.
pub fn transcript_shards(directory: &Path) -> Vec<PathBuf> {
    let Ok(listing) = std::fs::read_dir(directory) else {
        return Vec::new();
    };
    let mut names: Vec<OsString> = listing
        .filter_map(Result::ok)
        .map(|entry| entry.file_name())
        .filter(|name| is_transcript_shard(name.as_bytes()))
        .collect();
    names.sort();
    names
        .into_iter()
        .map(|name| join(directory, name.as_bytes()))
        .collect()
}

/// `"$directory/$name"`: Perl joins the two with one slash and keeps what it is given.
pub fn join(directory: &Path, name: &[u8]) -> PathBuf {
    let mut path = directory.as_os_str().as_bytes().to_vec();
    path.push(b'/');
    path.extend_from_slice(name);
    PathBuf::from(OsString::from_vec(path))
}

fn is_transcript_shard(name: &[u8]) -> bool {
    let Some(stem) = name.strip_suffix(b".gz") else {
        return false;
    };
    let Some(dash) = stem.iter().position(|byte| *byte == b'-') else {
        return false;
    };
    let digits = |part: &[u8]| !part.is_empty() && part.iter().all(u8::is_ascii_digit);
    digits(&stem[..dash]) && digits(&stem[dash + 1..])
}

/// What `IO::Uncompress::Gunzip::gunzip($file => \$buffer)` leaves in the buffer: the
/// first member of a gzip file, without a check of its trailer; the file's own bytes
/// when it does not open with the gzip magic. The error is why `gunzip` returns false,
/// in this module's words.
pub fn gunzip(path: &Path) -> Result<Vec<u8>, String> {
    let bytes = std::fs::read(path).map_err(|error| {
        let text = error.to_string();
        let text = match text.rfind(" (os error ") {
            Some(at) if text.ends_with(')') => &text[..at],
            _ => &text,
        };
        format!("cannot read the file: {text}")
    })?;
    if bytes.len() < 2 || bytes[..2] != [0x1f, 0x8b] {
        return Ok(bytes);
    }
    let body = gzip_body(&bytes)?;
    let mut out = Vec::new();
    // The last four bytes of a member hold its length modulo 2^32: right for a file of
    // one member under 4 GiB, and a first size for the buffer.
    if bytes.len() >= 18 {
        let tail: [u8; 4] = bytes[bytes.len() - 4..].try_into().expect("four bytes");
        let hint = u32::from_le_bytes(tail) as usize;
        if hint as u64 <= (bytes.len() as u64).saturating_mul(1032) {
            out.reserve_exact(hint);
        }
    }
    DeflateDecoder::new(body)
        .read_to_end(&mut out)
        .map_err(|_| "the deflate stream is damaged or short".to_string())?;
    Ok(out)
}

/// The deflate stream of a gzip member: what follows its header.
fn gzip_body(bytes: &[u8]) -> Result<&[u8], String> {
    const FHCRC: u8 = 0x02;
    const FEXTRA: u8 = 0x04;
    const FNAME: u8 = 0x08;
    const FCOMMENT: u8 = 0x10;
    const RESERVED: u8 = 0xe0;
    let short = || "the gzip header is short".to_string();
    if bytes.len() < 10 {
        return Err(short());
    }
    if bytes[2] != 8 {
        return Err(format!(
            "the gzip header names compression method {}",
            bytes[2]
        ));
    }
    let flags = bytes[3];
    if flags & RESERVED != 0 {
        return Err("the gzip header sets reserved flag bits".to_string());
    }
    let mut at = 10;
    if flags & FEXTRA != 0 {
        let length = bytes.get(at..at + 2).ok_or_else(short)?;
        at += 2 + usize::from(u16::from_le_bytes([length[0], length[1]]));
    }
    for flag in [FNAME, FCOMMENT] {
        if flags & flag != 0 {
            let rest = bytes.get(at..).ok_or_else(short)?;
            let end = rest.iter().position(|byte| *byte == 0).ok_or_else(short)?;
            at += end + 1;
        }
    }
    if flags & FHCRC != 0 {
        at += 2;
    }
    bytes.get(at..).ok_or_else(short)
}

/// `gunzip` on the file, then `fd_retrieve` on its output. `Ok(None)` is an image that
/// ends early: `fd_retrieve` returns undef there.
pub fn read_shard(path: &Path) -> Result<Option<Shard>, ShardError> {
    let bytes = gunzip(path).map_err(ShardError::Gunzip)?;
    Ok(fd_retrieve(bytes)
        .map_err(ShardError::Retrieve)?
        .map(Shard::new))
}

/// The transcripts of a retrieved shard, as
/// `ref($data) eq 'HASH' ? (map { @{ $_ || [] } } values %$data) : @{ $data || [] }`
/// lists them: the items of every list of an unblessed hash, in the order the image
/// stores its lists, else the items of the root array. An image that ended early holds
/// none.
pub fn transcripts(shard: Option<&Shard>) -> Result<Vec<Sv>, Die> {
    let Some(shard) = shard else {
        return Ok(Vec::new());
    };
    let image = &shard.image;
    let perl = Perl::new(image);
    let root = image.root();
    let plain_hash = matches!(image.value(root), Value::Hash { .. })
        && matches!(image.class(root), None | Some(b"HASH"));
    if !plain_hash {
        return match image.items(root) {
            Some(items) => Ok(items
                .iter()
                .map(|item| shard.read(perl.element(*item)))
                .collect()),
            None => Err(Die::new("Not an ARRAY reference")),
        };
    }
    let mut all = Vec::new();
    for (_, value) in crate::dumper::hash_entries_in_stored_order(image, root) {
        let list = shard.read(Sv::Node(value));
        // `$_ || []`: a false value is an empty list.
        if perl.is_true(list) {
            all.extend(
                perl.items_of(list)?
                    .into_iter()
                    .map(|item| shard.read(item)),
            );
        }
    }
    Ok(all)
}
