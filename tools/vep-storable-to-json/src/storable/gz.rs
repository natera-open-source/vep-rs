// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Opening an image from a file: plain, or gzip-compressed as an Ensembl VEP cache
//! stores its shards.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use super::error::{ErrorKind, StorableError};
use super::image::Image;
use super::reader;

fn io_error(error: &std::io::Error) -> StorableError {
    StorableError::new(ErrorKind::Io(error.to_string()), 0)
}

impl Image {
    /// Decodes an image held in memory, header first. Bytes after the image stay in the
    /// buffer and are counted by [`Image::trailing_len`].
    pub fn read(bytes: Vec<u8>) -> Result<Image, StorableError> {
        reader::decode(bytes)
    }

    /// Reads and decodes a file `store` or `nstore` wrote.
    pub fn read_file(path: &Path) -> Result<Image, StorableError> {
        Image::read(std::fs::read(path).map_err(|error| io_error(&error))?)
    }

    /// Reads, decompresses and decodes a gzip-compressed image, every member of the
    /// gzip file in turn as `gzip -dc` writes them.
    pub fn read_gz_file(path: &Path) -> Result<Image, StorableError> {
        Image::read(gunzip_file(path)?)
    }
}

/// The decompressed bytes of a gzip file, every member in turn. A damaged or short
/// stream is [`ErrorKind::Gzip`]; a file that cannot be read is [`ErrorKind::Io`].
pub fn gunzip_file(path: &Path) -> Result<Vec<u8>, StorableError> {
    let mut file = File::open(path).map_err(|error| io_error(&error))?;
    // The last four bytes of a gzip member hold its decompressed length modulo 2^32:
    // a first size for the buffer, right for a file of one member under 4 GiB.
    let mut capacity = 0usize;
    let compressed_len = file.metadata().map_err(|error| io_error(&error))?.len();
    if compressed_len >= 18 {
        let mut tail = [0u8; 4];
        file.seek(SeekFrom::End(-4))
            .and_then(|_| file.read_exact(&mut tail))
            .and_then(|()| file.seek(SeekFrom::Start(0)))
            .map_err(|error| io_error(&error))?;
        capacity = u32::from_le_bytes(tail) as usize;
        // A length that no deflate stream of this size can reach is a damaged trailer.
        if capacity as u64 > compressed_len.saturating_mul(1032) {
            capacity = 0;
        }
    }
    let mut bytes = Vec::with_capacity(capacity);
    MultiGzDecoder::new(file)
        .read_to_end(&mut bytes)
        .map_err(|error| StorableError::new(ErrorKind::Gzip(error.to_string()), 0))?;
    Ok(bytes)
}
