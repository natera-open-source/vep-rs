// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One compared entry and its fields.

use indexmap::IndexMap;

/// Field values by name in insertion order; a name set twice keeps its first place and
/// takes the last value.
pub type Fields = IndexMap<String, String>;

/// One row of the default or tab format, one CSQ entry or record line of the VCF format,
/// or one consequence object or record object of the JSON format.
///
/// `key` is the entry's identity within its format: the Location (default, tab) or the
/// input record ordinal (VCF, JSON), the allele, the feature and the feature type, plus a
/// `#n` counter on exact duplicates.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Entry {
    pub key: String,
    pub record: Option<i64>,
    pub location: Option<String>,
    pub allele: String,
    pub feature: String,
    pub consequence: String,
    pub fields: Fields,
    /// Whether the allele is an insertion's or a deletion's, read from the record's own
    /// alleles; `None` where they say nothing.
    pub indel: Option<bool>,
    /// Whether the record's alleles make every sequence allele an insertion; `None` where
    /// the entry's record spells no alleles.
    pub insertion: Option<bool>,
    /// The reference bases the allele stands in for at the entry's Location; empty for a
    /// pure insertion's entry, `None` where no record was read.
    pub reference: Option<String>,
    /// The ordinal of the input line a VCF entry's line is written from; `None` for any
    /// other entry.
    pub line: Option<i64>,
    /// The Location the default format writes for a VCF or JSON entry's allele; `None` for
    /// any other entry.
    pub allele_location: Option<String>,
}
