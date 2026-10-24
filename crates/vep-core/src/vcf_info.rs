// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The INFO column read as the hash Ensembl VEP's parser builds from it (ensembl-io
//! `BaseVCF4.pm` lines 581-583: over the `;`-split fields, `my ($key,$value) =
//! split('=',$info); $info_data{$key} = $value;`), shared by every VCF reader. The last
//! occurrence of a key wins; a bare key sets undef; `KEY=` sets the empty string; the
//! two-element list assignment keeps the text between the first `=` and a second one
//! (`END2=100=5` is `100`). Keys match exactly, case and spaces included.
//!
//! A reader that looks up several keys of one record holds its column as an
//! [`InfoFields`], split into fields once by the first lookup and read by every
//! later one; a single lookup ([`info_field`]) scans the column for its key and
//! allocates nothing. Both read the hash the same way.

use std::cell::OnceCell;

/// The most fields a column is split into: a column of more is read by scans of
/// its text, so the split's memory is bounded whatever a line holds.
const SPLIT_FIELD_LIMIT: usize = 1 << 16;

/// One `;`-separated field as Perl's hash holds it: the key and `Some(value)`
/// for `key=value`, `None` for a bare key.
type Field<'a> = (&'a str, Option<&'a str>);

/// One record's INFO column, split into its `Field`s on its first lookup and
/// held for every lookup after it; a record nothing looks up is never split. The
/// slices borrow the column text. A column of more than `SPLIT_FIELD_LIMIT`
/// fields holds no split (`None`) and every lookup on it scans the text.
#[derive(Debug)]
pub struct InfoFields<'a> {
    text: &'a str,
    fields: OnceCell<Option<Vec<Field<'a>>>>,
}

impl<'a> InfoFields<'a> {
    /// The view over a record's INFO column, nothing split yet.
    pub fn new(text: &'a str) -> Self {
        InfoFields {
            text,
            fields: OnceCell::new(),
        }
    }

    /// The fields in written order, each as [`field_entry`] reads it; none for `.`
    /// or an empty column; `None` for a column over [`SPLIT_FIELD_LIMIT`] fields.
    fn fields(&self) -> Option<&[Field<'a>]> {
        self.fields
            .get_or_init(|| {
                if self.text == "." || self.text.is_empty() {
                    return Some(Vec::new());
                }
                let count = self.text.matches(';').count() + 1;
                if count > SPLIT_FIELD_LIMIT {
                    return None;
                }
                let mut fields = Vec::with_capacity(count);
                fields.extend(self.text.split(';').map(field_entry));
                Some(fields)
            })
            .as_deref()
    }

    /// The last occurrence of `key`: `Some(Some(value))` for `key=value`, `Some(None)`
    /// for a bare flag, `None` when the column is `.`, empty or lacks the key.
    pub fn field(&self, key: &str) -> Option<Option<&'a str>> {
        match self.fields() {
            Some(fields) => fields
                .iter()
                .rev()
                .find(|(k, _)| *k == key)
                .map(|(_, value)| *value),
            None => info_field(self.text, key),
        }
    }

    /// The value of `key`, `""` for a bare flag, `None` when absent.
    pub fn value(&self, key: &str) -> Option<&'a str> {
        self.field(key).map(|value| value.unwrap_or(""))
    }

    /// The value of `key` when it is written `key=value`, the empty value included:
    /// Perl's `defined $info->{KEY}`, which a bare `KEY` flag fails because its hash
    /// value is undef.
    pub fn value_written(&self, key: &str) -> Option<&'a str> {
        self.field(key).flatten()
    }

    /// Whether a lookup has split the column into held fields.
    #[cfg(test)]
    fn is_split(&self) -> bool {
        self.fields.get().is_some_and(Option::is_some)
    }
}

/// One `;`-separated field as the hash entry Perl's two-element list assignment
/// makes of it: `(key, Some(value))` for `key=value`, the value cut at a second
/// `=`, and `(key, None)` for a bare key.
fn field_entry(field: &str) -> Field<'_> {
    let Some((key, rest)) = field.split_once('=') else {
        return (field, None);
    };
    let value = rest.split_once('=').map_or(rest, |(value, _)| value);
    (key, Some(value))
}

/// The last occurrence of `key`, as [`InfoFields::field`] reads it, by one scan of
/// the column from its end.
pub fn info_field<'a>(info: &'a str, key: &str) -> Option<Option<&'a str>> {
    if info == "." || info.is_empty() {
        return None;
    }
    info.rsplit(';')
        .map(field_entry)
        .find(|(k, _)| *k == key)
        .map(|(_, value)| value)
}

/// The value of `key`, `""` for a bare flag, `None` when absent.
pub fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info_field(info, key).map(|value| value.unwrap_or(""))
}

/// The value of `key` when it is written `key=value`, the empty value included: Perl's
/// `defined $info->{KEY}`, which a bare `KEY` flag fails because its hash value is undef.
pub fn info_value_written<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info_field(info, key).flatten()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The field both readers find, asserting the scan and the split view agree.
    fn field<'a>(info: &'a str, key: &str) -> Option<Option<&'a str>> {
        let scanned = info_field(info, key);
        assert_eq!(
            InfoFields::new(info).field(key),
            scanned,
            "{info:?} {key:?}"
        );
        scanned
    }

    /// `my ($key, $value) = split('=', $info)` (BaseVCF4.pm:582): `END=` is the empty
    /// string, `END2=100=5` is `100`, a bare `CHR2` is undef, and `$info_data{$key} =
    /// $value` keeps the last occurrence, a bare key written last leaving undef.
    #[test]
    fn the_info_hash_is_read_as_perl_builds_it() {
        assert_eq!(field("END=", "END"), Some(Some("")));
        assert_eq!(field("END2=100=5", "END2"), Some(Some("100")));
        assert_eq!(field("CHR2", "CHR2"), Some(None));
        assert_eq!(field("SVTYPE=DEL", "END"), None);
        assert_eq!(field(".", "END"), None);
        assert_eq!(field("", "END"), None);
        assert_eq!(field("END=100;END=200", "END"), Some(Some("200")));
        assert_eq!(field("END=100;END", "END"), Some(None));
        assert_eq!(field("END;END=100", "END"), Some(Some("100")));
        assert_eq!(info_value("CHR2", "CHR2"), Some(""));
        assert_eq!(InfoFields::new("CHR2").value("CHR2"), Some(""));
        assert_eq!(info_value_written("CHR2", "CHR2"), None);
        assert_eq!(InfoFields::new("CHR2").value_written("CHR2"), None);
        assert_eq!(info_value_written("CHR2=", "CHR2"), Some(""));
        assert_eq!(InfoFields::new("CHR2=").value_written("CHR2"), Some(""));
        assert_eq!(info_value("END2=100=5;X=1", "END2"), Some("100"));
        assert_eq!(InfoFields::new("END2=100=5;X=1").value("END2"), Some("100"));
    }

    /// Keys match as written: Ensembl VEP's hash has no entry `END` for `end=100` or
    /// ` END=100`, so a deletion carrying only those is incomplete (Parser/VCF.pm:555-558).
    #[test]
    fn a_key_matches_exactly() {
        for info in [
            "svtype=DEL;end=100",
            "SVTYPE=DEL; END=100",
            "SVTYPE=DEL;END =100",
        ] {
            assert_eq!(field(info, "END"), None, "{info}");
        }
        assert_eq!(field("SVTYPE=DEL; END=100", " END"), Some(Some("100")));
    }

    /// Every field is an entry, the empty field between two `;` (key `""`, undef)
    /// and a field opening with `=` (key `""`, the rest) included, as in Perl's hash;
    /// the lone `.` of an absent column is read as no entry (Perl's hash holds it as
    /// the key `.`, which nothing reads), while the `.` of `A=1;.` is a bare key.
    #[test]
    fn every_field_is_an_entry_and_the_absent_column_none() {
        assert_eq!(field("A=1;;=5;.", "A"), Some(Some("1")));
        assert_eq!(field("A=1;;=5;.", ""), Some(Some("5")));
        assert_eq!(field("A=1;;=5;.", "."), Some(None));
        assert_eq!(field(".", "."), None);
        assert_eq!(field("", ""), None);
    }

    /// The scan and the split view read every key of a column alike, the keys the
    /// column lacks included, over columns of repeated, bare, empty and `=`-bearing
    /// fields.
    #[test]
    fn the_scan_and_the_split_view_agree_on_every_key() {
        let columns = [
            "SVTYPE=DEL;END=200;SVLEN=-100;END=300;CIPOS=-5,5;IMPRECISE",
            "A=1;;=5;.;A;B==;C=x=y=z;A=",
            "END=;END;END=1",
            ";",
            "=",
            "key=value",
            ".",
            "",
        ];
        for column in columns {
            let fields = InfoFields::new(column);
            let keys = column
                .split(';')
                .map(|f| f.split_once('=').map_or(f, |(k, _)| k))
                .chain(["SVTYPE", "END", "MISSING", "", ".", "=", "A"]);
            for key in keys {
                assert_eq!(
                    fields.field(key),
                    info_field(column, key),
                    "{column:?} {key:?}"
                );
            }
        }
    }

    /// The column is split by the first lookup alone: a view nothing reads holds no
    /// fields, the first lookup splits it, and a second lookup on another key and a
    /// repeat of the first read the same split.
    #[test]
    fn the_column_is_split_once_on_the_first_lookup() {
        let fields = InfoFields::new("SVTYPE=DEL;END=200;SVLEN=-100;END=300");
        assert!(!fields.is_split());
        assert_eq!(fields.value("SVTYPE"), Some("DEL"));
        assert!(fields.is_split());
        let split = fields.fields().unwrap().as_ptr();
        assert_eq!(fields.value_written("END"), Some("300"));
        assert_eq!(fields.value("SVLEN"), Some("-100"));
        assert_eq!(fields.value("SVTYPE"), Some("DEL"));
        assert_eq!(fields.value("CIPOS"), None);
        assert_eq!(fields.fields().unwrap().as_ptr(), split);
        assert_eq!(fields.fields().unwrap().len(), 4);
    }

    /// A column of more than `SPLIT_FIELD_LIMIT` fields is never split: its lookups
    /// scan the text and read every key as the split would, the last occurrence of
    /// a repeated key and a key the column lacks included; one field fewer is split.
    #[test]
    fn a_column_over_the_field_limit_is_read_by_scans() {
        let fields = |count: usize| {
            (0..count)
                .map(|i| format!("K{i}={i}"))
                .collect::<Vec<_>>()
                .join(";")
        };
        let over = format!("{};END=1;END=2;X", fields(SPLIT_FIELD_LIMIT - 2));
        let view = InfoFields::new(&over);
        assert_eq!(view.field("END"), Some(Some("2")));
        assert!(!view.is_split());
        assert_eq!(view.field("K0"), Some(Some("0")));
        assert_eq!(view.field("X"), Some(None));
        assert_eq!(view.field("MISSING"), None);
        assert_eq!(view.field("K0"), info_field(&over, "K0"));
        let at_limit = fields(SPLIT_FIELD_LIMIT);
        let view = InfoFields::new(&at_limit);
        assert_eq!(view.field("K0"), Some(Some("0")));
        assert!(view.is_split());
    }
}
