// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader's one error type.

use std::fmt;

/// Why an image could not be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ErrorKind {
    /// The file could not be opened or read; the text is the operating system's.
    Io(String),
    /// The gzip stream is damaged or ends early; the text is the decoder's.
    Gzip(String),
    /// The bytes do not start with the `pst0` magic.
    NotStorable,
    /// A binary major other than 2: the formats before Storable 0.7 and any later major.
    UnsupportedVersion { major: u8, minor: u8 },
    /// A native-order header whose byte-order string is none of `1234`, `4321`,
    /// `12345678` and `87654321`.
    UnsupportedByteOrder(String),
    /// The image ends inside the header or inside a value.
    Truncated,
    /// An opcode the reader does not decode: a tied value, a hook, a code reference, a
    /// large array or hash, or a number the format does not define.
    UnsupportedOpcode { opcode: u8, name: &'static str },
    /// `SX_OBJECT` names a tag no earlier value carries.
    UnknownTag(u64),
    /// `SX_IX_BLESS` names a class index no earlier `SX_BLESS` introduced.
    UnknownClass(u64),
    /// A length the reader cannot hold or the format cannot mean: a negative key length,
    /// a string of 4 GiB or more, more than `u32::MAX` values, a double that is not 8 bytes.
    BadLength(&'static str),
    /// A flagged hash whose key is itself a stored scalar (`SHV_K_ISSV`).
    ScalarHashKey,
}

/// An [`ErrorKind`] with the byte offset in the decompressed image where it was met: the
/// offset of the opcode for a value, of the field for a header, 0 for a file that could
/// not be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StorableError {
    pub kind: ErrorKind,
    pub offset: u64,
}

impl StorableError {
    pub(super) fn new(kind: ErrorKind, offset: usize) -> Self {
        Self {
            kind,
            offset: offset as u64,
        }
    }
}

impl fmt::Display for StorableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let at = self.offset;
        match &self.kind {
            ErrorKind::Io(text) => write!(f, "{text}"),
            ErrorKind::Gzip(text) => write!(f, "gzip stream: {text}"),
            ErrorKind::NotStorable => write!(f, "File is not a perl storable"),
            ErrorKind::UnsupportedVersion { major, minor } => write!(
                f,
                "Storable binary image v{major}.{minor} is not a version 2 image"
            ),
            ErrorKind::UnsupportedByteOrder(order) => write!(
                f,
                "Storable image has byte order '{order}' at offset {at}; expected 1234, 4321, 12345678 or 87654321"
            ),
            ErrorKind::Truncated => write!(f, "Storable image ends early at offset {at}"),
            ErrorKind::UnsupportedOpcode { opcode, name } => write!(
                f,
                "unsupported Storable opcode {name} ({opcode}) at offset {at}"
            ),
            ErrorKind::UnknownTag(tag) => write!(
                f,
                "Object #{tag} should have been retrieved already (SX_OBJECT at offset {at})"
            ),
            ErrorKind::UnknownClass(index) => write!(
                f,
                "Class name #{index} should have been seen already (SX_IX_BLESS at offset {at})"
            ),
            ErrorKind::BadLength(what) => write!(f, "{what} at offset {at}"),
            ErrorKind::ScalarHashKey => write!(
                f,
                "unsupported Storable hash key held as a scalar (SHV_K_ISSV) at offset {at}"
            ),
        }
    }
}

impl std::error::Error for StorableError {}
