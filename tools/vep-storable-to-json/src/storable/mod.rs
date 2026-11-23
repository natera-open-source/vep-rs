// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A reader of Perl Storable images: the files `Storable::store` and `Storable::nstore`
//! write, binary major 2, with no Perl at run time.
//!
//! [`Image::read`], [`Image::read_file`] and [`Image::read_gz_file`] decode a whole
//! image into an [`Image`]: one node a stored value, numbered by its Storable tag, with
//! what Perl sees of it after `retrieve` kept apart: whether a scalar was stored as an
//! integer, a double or a string and with which UTF-8 flag, the class a value is blessed
//! into, a reference as a node of its own that names its target, and a value reached
//! twice as one node. Arrays keep the difference between an element that does not exist
//! and one that is undef; hashes keep their stored order, their key flags and the
//! restricted flag.
//!
//! Decoded: scalars (`SX_BYTE`, `SX_INTEGER`, `SX_NETINT`, `SX_DOUBLE`, the four string
//! opcodes, `SX_UNDEF`, the immortals `SX_SV_UNDEF`, `SX_SV_YES`, `SX_SV_NO`, the
//! booleans), references (`SX_REF`, `SX_WEAKREF`, `SX_OVERLOAD`, `SX_WEAKOVERLOAD`),
//! `SX_ARRAY` with `SX_SVUNDEF_ELEM`, `SX_HASH`, `SX_FLAG_HASH`, `SX_BLESS`,
//! `SX_IX_BLESS`, `SX_OBJECT`, version strings, `SX_REGEXP`, and the string and tag
//! forms of `SX_LOBJECT`. Refused with the opcode's name, its number and its offset:
//! tied values, `SX_HOOK`, `SX_CODE`, the array and hash forms of `SX_LOBJECT`, and any
//! number the format does not define. Images of a binary major other than 2 are refused.
//!
//! The modules of this directory use nothing else of the package, so another tool
//! compiles them with `#[path = ".../src/storable/mod.rs"] mod storable;` and a
//! dependency on `flate2`.

mod error;
mod gz;
mod header;
mod image;
mod reader;

pub use error::{ErrorKind, StorableError};
pub use gz::gunzip_file;
pub use header::{ByteOrder, Header, NativeLayout};
pub use image::{HashEntry, Image, KeyFlags, NodeId, Run, Span, Value};
pub use reader::op;
