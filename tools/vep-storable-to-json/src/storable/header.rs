// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The file header `store` and `nstore` write before the first value: the `pst0` magic,
//! one byte holding the binary major and the network-order bit, the binary minor, and
//! for a native-order image the writer's byte-order string and type sizes.

use super::error::{ErrorKind, StorableError};

pub(super) const MAGIC: &[u8; 4] = b"pst0";

/// The byte order of the multi-byte fields that follow the image's own order: lengths,
/// native integers and doubles. Object tags and network integers are big-endian in
/// every image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteOrder {
    Little,
    Big,
}

/// What a native-order header records of the machine that wrote the image. A
/// network-order header records none of it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeLayout {
    /// The writer's byte-order string as stored (`12345678` on a little-endian perl
    /// with 8-byte integers). Its length is the size of a native integer.
    pub byte_order_string: String,
    pub int_size: u8,
    pub long_size: u8,
    pub pointer_size: u8,
    /// Absent before binary minor 2.
    pub double_size: Option<u8>,
}

/// The decoded file header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// Binary major; the reader accepts 2.
    pub major: u8,
    /// Binary minor. A minor past the opcodes the reader knows is accepted, and an
    /// opcode it does not know then fails by number.
    pub minor: u8,
    /// True for an image written by `nstore`.
    pub network_order: bool,
    /// The order of lengths, native integers and doubles: big-endian in network order.
    pub byte_order: ByteOrder,
    /// Bytes of an `SX_INTEGER` value: the length of the byte-order string, 8 in network
    /// order, where the opcode is not written.
    pub integer_size: u8,
    /// Present for a native-order image.
    pub native: Option<NativeLayout>,
    /// Offset of the first value.
    pub body_offset: usize,
}

impl Header {
    /// Reads the header alone from the first bytes of an image.
    pub fn parse(data: &[u8]) -> Result<Header, StorableError> {
        if data.len() < MAGIC.len() || &data[..MAGIC.len()] != MAGIC {
            // A file shorter than the magic that is a prefix of it ends early; anything
            // else is another kind of file.
            let kind = if data.len() < MAGIC.len() && MAGIC.starts_with(data) {
                ErrorKind::Truncated
            } else {
                ErrorKind::NotStorable
            };
            return Err(StorableError::new(kind, 0));
        }
        let mut pos = MAGIC.len();
        let next = |pos: &mut usize| -> Result<u8, StorableError> {
            let byte = *data
                .get(*pos)
                .ok_or_else(|| StorableError::new(ErrorKind::Truncated, *pos))?;
            *pos += 1;
            Ok(byte)
        };
        let version_byte = next(&mut pos)?;
        let network_order = version_byte & 1 == 1;
        let major = version_byte >> 1;
        // Only a major above 1 is followed by a minor byte.
        let minor = if major > 1 { next(&mut pos)? } else { 0 };
        if major != 2 {
            return Err(StorableError::new(
                ErrorKind::UnsupportedVersion { major, minor },
                MAGIC.len(),
            ));
        }
        if network_order {
            return Ok(Header {
                major,
                minor,
                network_order,
                byte_order: ByteOrder::Big,
                integer_size: 8,
                native: None,
                body_offset: pos,
            });
        }
        let order_offset = pos;
        let order_len = next(&mut pos)? as usize;
        let order_bytes = data
            .get(pos..pos + order_len)
            .ok_or_else(|| StorableError::new(ErrorKind::Truncated, pos))?;
        pos += order_len;
        let byte_order = match order_bytes {
            b"1234" | b"12345678" => ByteOrder::Little,
            b"4321" | b"87654321" => ByteOrder::Big,
            other => {
                return Err(StorableError::new(
                    ErrorKind::UnsupportedByteOrder(String::from_utf8_lossy(other).into_owned()),
                    order_offset,
                ))
            }
        };
        let int_size = next(&mut pos)?;
        let long_size = next(&mut pos)?;
        let pointer_size = next(&mut pos)?;
        let double_size = if minor >= 2 {
            Some(next(&mut pos)?)
        } else {
            None
        };
        Ok(Header {
            major,
            minor,
            network_order,
            byte_order,
            integer_size: order_len as u8,
            native: Some(NativeLayout {
                byte_order_string: String::from_utf8_lossy(order_bytes).into_owned(),
                int_size,
                long_size,
                pointer_size,
                double_size,
            }),
            body_offset: pos,
        })
    }
}
