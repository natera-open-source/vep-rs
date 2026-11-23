// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reading one shard: `gzip -dc` on the file, Storable's `fd_retrieve` on its output,
//! then the check of gzip's exit.
//!
//! A damaged file is heard of in that order. `fd_retrieve` dies first when what gzip
//! wrote is no Storable image; when the image is whole or merely short, gzip's exit
//! status is reported if gzip failed; and a short image from a
//! gzip that succeeded is `fd_retrieve returned no data`. The outcomes here follow the
//! same order, with gzip's part done in this process.

use std::io::Read;
use std::path::Path;

use flate2::bufread::GzDecoder;

use crate::perl::Die;
use crate::storable::{ErrorKind, Header, Image, StorableError, Value};

/// gzip writes its output a window at a time. On data that is not a deflate stream it
/// exits without writing the window it was filling; on a file that ends early it
/// writes that window first.
const GZIP_WINDOW: usize = 32 * 1024;
/// What a pipe holds: `gzip -dc` finishes when no more than this follows the image,
/// and is ended by a signal when its reader closes the pipe on more.
const PIPE_CAPACITY: usize = 64 * 1024;
/// The binary minor of the Storable whose messages these are (3.23, binary 2.11) and
/// the last opcode it knows.
const STORABLE_MINOR: u8 = 11;
const STORABLE_LAST_OPCODE: u8 = 33;

enum GzipEnd {
    /// Exit 0.
    Clean,
    /// Exit 1: a stream that is damaged, short or no gzip at all.
    Error,
    /// Exit 2: a whole member followed by bytes that are neither a member nor zeros.
    Warning,
}

/// What `gzip -dc` writes for the file's bytes, and how it exits.
fn gzip_dc(compressed: &[u8]) -> (Vec<u8>, GzipEnd) {
    let mut out = Vec::new();
    if compressed.len() >= 18 {
        // The last four bytes of a member hold its length modulo 2^32.
        let tail: [u8; 4] = compressed[compressed.len() - 4..]
            .try_into()
            .expect("four bytes");
        let hint = u32::from_le_bytes(tail) as usize;
        if hint as u64 <= (compressed.len() as u64).saturating_mul(1032) {
            out.reserve_exact(hint);
        }
    }
    let mut rest = compressed;
    let mut members = 0;
    loop {
        if rest.is_empty() {
            return (
                out,
                if members == 0 {
                    GzipEnd::Error
                } else {
                    GzipEnd::Clean
                },
            );
        }
        if members > 0 && !rest.starts_with(&[0x1f, 0x8b]) {
            // Zero padding after a member is ignored; anything else is a warning.
            let end = if rest.iter().all(|byte| *byte == 0) {
                GzipEnd::Clean
            } else {
                GzipEnd::Warning
            };
            return (out, end);
        }
        let member_start = out.len();
        let mut decoder = GzDecoder::new(rest);
        match decoder.read_to_end(&mut out) {
            Ok(_) => {
                rest = decoder.into_inner();
                members += 1;
            }
            Err(error) => {
                // A short file and a member whose checksum or length does not match
                // leave everything decoded behind; damaged data leaves whole windows.
                let all_written = error.kind() == std::io::ErrorKind::UnexpectedEof
                    || error.to_string().contains("checksum");
                if !all_written {
                    let whole_windows = (out.len() - member_start) / GZIP_WINDOW * GZIP_WINDOW;
                    out.truncate(member_start + whole_windows);
                }
                return (out, GzipEnd::Error);
            }
        }
    }
}

/// The bytes `gzip -dc` writes for a file's bytes, whatever its exit status: what a
/// reader of its output gets when it does not ask how gzip ended.
pub fn gzip_dc_output(compressed: &[u8]) -> Vec<u8> {
    gzip_dc(compressed).0
}

/// The text Storable dies with for an image the reader refuses, where it has one of
/// its own; else the reader's.
fn storable_text(error: &StorableError, available: usize) -> Die {
    match &error.kind {
        // Storable reads eleven bytes before it tells another kind of file from a
        // file that ends inside the magic.
        ErrorKind::NotStorable if available >= 11 => Die::new("File is not a perl storable"),
        ErrorKind::NotStorable => Die::new("Magic number checking on storable file failed"),
        ErrorKind::UnsupportedVersion { major, minor } if *major > 2 => Die::new(format!(
            "Storable binary image v{major}.{minor} more recent than I am (v2.{STORABLE_MINOR})"
        )),
        _ => Die::new(error.to_string()),
    }
}

/// The checks Storable 3.23 on a little-endian perl with 8-byte integers makes of a
/// native-order header: an image of another layout is refused.
fn native_layout(header: &Header) -> Result<(), Die> {
    let Some(native) = &header.native else {
        return Ok(());
    };
    let refuse = |text: &str| Err(Die::new(text));
    if native.byte_order_string != "12345678" {
        return refuse("Byte order is not compatible");
    }
    if native.int_size != 4 {
        return refuse("Integer size is not compatible");
    }
    if native.long_size != 8 {
        return refuse("Long integer size is not compatible");
    }
    if native.pointer_size != 8 {
        return refuse("Pointer size is not compatible");
    }
    if native.double_size.is_some_and(|size| size != 8) {
        return refuse("Double size is not compatible");
    }
    Ok(())
}

/// Storable's `fd_retrieve` on a stream that holds `bytes`. `Ok(None)` is an image that
/// ends early, where `fd_retrieve` returns undef; `Err` is the text Storable dies with.
pub fn fd_retrieve(bytes: Vec<u8>) -> Result<Option<Image>, Die> {
    let available = bytes.len();
    // The header first: Storable refuses a layout it was not built for before it
    // reads a value.
    match Header::parse(&bytes) {
        Ok(header) => native_layout(&header)?,
        Err(error) if error.kind == ErrorKind::Truncated => {
            return Err(Die::new("Magic number checking on storable file failed"));
        }
        Err(error) => return Err(storable_text(&error, available)),
    }
    let image = match Image::read(bytes) {
        Ok(image) => image,
        Err(error) if error.kind == ErrorKind::Truncated => return Ok(None),
        Err(error) => return Err(storable_text(&error, available)),
    };
    // A boolean opcode is past the last one that Storable knows.
    if image.header().minor > STORABLE_MINOR {
        let boolean = image.node_ids().find_map(|id| match image.value(id) {
            Value::Bool(true) => Some(34),
            Value::Bool(false) => Some(35),
            _ => None,
        });
        if let Some(opcode) = boolean {
            let header = image.header();
            return Err(Die::new(format!(
                "Storable binary image v{}.{} contains data of type {opcode}. This Storable is v2.{STORABLE_MINOR} and can only handle data types up to {STORABLE_LAST_OPCODE}",
                header.major, header.minor
            )));
        }
    }
    Ok(Some(image))
}

/// Reads a shard file. `Ok(None)` is an image that ends early from a gzip that
/// succeeded: `fd_retrieve` returns undef there, and what that means depends on the
/// kind of shard. `Err` is the text Perl dies with.
pub fn load_shard(path: &Path) -> Result<Option<Image>, Die> {
    // A file gzip cannot read gives `fd_retrieve` an empty stream.
    let compressed = std::fs::read(path).unwrap_or_default();
    let (bytes, gzip_end) = gzip_dc(&compressed);
    drop(compressed);
    let gzip_failure = match gzip_end {
        GzipEnd::Clean => None,
        GzipEnd::Error => Some("gzip -dc failed (wait status 256)"),
        GzipEnd::Warning => Some("gzip -dc failed (wait status 512)"),
    };
    // What Storable dies with comes before gzip's exit status is looked at.
    let image = fd_retrieve(bytes)?;
    if let Some(text) = gzip_failure {
        return Err(Die::new(text));
    }
    let Some(image) = image else {
        return Ok(None);
    };
    if image.trailing_len() > PIPE_CAPACITY {
        // gzip is still writing when the pipe is closed.
        return Err(Die::new("gzip -dc failed (wait status 13)"));
    }
    Ok(Some(image))
}
