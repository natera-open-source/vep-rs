// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run: every transcript shard of the canonical chromosomes, in the order the
//! chromosomes sort as strings and the shards sort by name, one line a transcript the
//! first time its id is met.

use std::collections::HashMap;
use std::ffi::OsString;
use std::io::Write;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

use crate::digest::{digest_shard, ContentRule, ShardDigest};
use crate::pool::run_ordered;
use crate::shards::{join, read_shard, transcript_shards, ShardError, CANONICAL};
use crate::text::Text;

/// The words that open the summary line.
const TAG: &str = "[compare_storable_cache_versions]";

/// What the command line asks for.
pub struct Request {
    pub cache_dir: Vec<u8>,
    pub out: Vec<u8>,
    pub threads: usize,
    pub rule: ContentRule,
}

struct Task {
    contig: &'static str,
    path: PathBuf,
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

/// The exit code of a `die` with no error number behind it.
pub const DIE: i32 = 255;

fn path_of(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

/// The shards of the run in its order: the chromosomes as strings sort, a
/// chromosome's shards by name.
fn tasks(cache_dir: &Path) -> Vec<Task> {
    let mut contigs = CANONICAL;
    contigs.sort_unstable();
    let mut all = Vec::new();
    for contig in contigs {
        let directory = join(cache_dir, contig.as_bytes());
        if !directory.is_dir() {
            continue;
        }
        all.extend(
            transcript_shards(&directory)
                .into_iter()
                .map(|path| Task { contig, path }),
        );
    }
    all
}

/// Runs the digest. Returns the exit code; `stderr` takes the warnings, the summary
/// line and the line of a run that died.
pub fn run(request: &Request, stderr: &mut dyn Write) -> i32 {
    let out_path = path_of(&request.out);
    let file = match std::fs::File::create(&out_path) {
        Ok(file) => file,
        Err(error) => {
            let mut line = b"cannot write ".to_vec();
            line.extend_from_slice(&request.out);
            line.extend_from_slice(format!(": {}\n", os_text(&error)).as_bytes());
            stderr.write_all(&line).ok();
            return die_code(&error);
        }
    };
    let mut out = std::io::BufWriter::with_capacity(1 << 16, file);
    let cache_dir = path_of(&request.cache_dir);
    let tasks = tasks(&cache_dir);
    let rule = request.rule;
    let work = |task: &Task| -> Result<ShardDigest, ShardError> {
        let shard = read_shard(&task.path)?;
        Ok(digest_shard(shard.as_ref(), task.contig.as_bytes(), rule))
    };

    let mut seen: HashMap<Vec<u8>, [u8; 32]> = HashMap::new();
    let (mut files, mut records, mut duplicates, mut conflicts) = (0u64, 0u64, 0u64, 0u64);
    let died: Option<Vec<u8>> = run_ordered(&tasks, request.threads, work, |pool| {
        for (index, task) in tasks.iter().enumerate() {
            let digest = match pool.take(index) {
                Ok(digest) => digest,
                Err(ShardError::Gunzip(reason)) => {
                    let mut line = b"gunzip ".to_vec();
                    line.extend_from_slice(task.path.as_os_str().as_bytes());
                    line.extend_from_slice(format!(" failed: {reason}").as_bytes());
                    return Some(line);
                }
                Err(ShardError::Retrieve(die)) => return Some(die.0),
            };
            files += 1;
            for record in digest.records {
                records += 1;
                match seen.get(&record.id.characters()) {
                    Some(first) => {
                        duplicates += 1;
                        if *first != record.content {
                            conflicts += 1;
                            let mut line = Text::bytes(&b"WARN: "[..]);
                            line.push(&record.id.bytes, record.id.utf8);
                            line.push(b" differs between two shards of ", false);
                            line.push(&request.cache_dir, false);
                            line.push(b"\n", false);
                            stderr.write_all(&line.printed()).ok();
                        }
                    }
                    None => {
                        seen.insert(record.id.characters(), record.content);
                        if let Err(error) = out.write_all(&record.line) {
                            return Some(write_failure(&request.out, &error));
                        }
                    }
                }
            }
            if let Some(die) = digest.failure {
                return Some(die.0);
            }
        }
        None
    });
    let flushed = out.flush();
    if let Some(mut line) = died {
        line.push(b'\n');
        stderr.write_all(&line).ok();
        return DIE;
    }
    if let Err(error) = flushed {
        let mut line = write_failure(&request.out, &error);
        line.push(b'\n');
        stderr.write_all(&line).ok();
        return DIE;
    }
    let mut line = format!("{TAG} ").into_bytes();
    line.extend_from_slice(&request.cache_dir);
    line.extend_from_slice(
        format!(
            ": {files} shard files, {records} transcript records, {} distinct transcripts, {duplicates} shard-boundary duplicates, {conflicts} conflicting duplicates\n",
            seen.len()
        )
        .as_bytes(),
    );
    stderr.write_all(&line).ok();
    i32::from(conflicts > 0)
}

fn write_failure(out: &[u8], error: &std::io::Error) -> Vec<u8> {
    let mut line = b"cannot write ".to_vec();
    line.extend_from_slice(out);
    line.extend_from_slice(format!(": {}", os_text(error)).as_bytes());
    line
}
