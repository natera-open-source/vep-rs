// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The INFO column read as the hash Ensembl VEP's parser builds from it (ensembl-io
//! `BaseVCF4.pm` lines 581-583: over the `;`-split fields, `my ($key,$value) =
//! split('=',$info); $info_data{$key} = $value;`), shared by every VCF reader. The last
//! occurrence of a key wins; a bare key sets undef; `KEY=` sets the empty string; the
//! two-element list assignment keeps the text between the first `=` and a second one
//! (`END2=100=5` is `100`). Keys match exactly, case and spaces included.

/// The last occurrence of `key`: `Some(Some(value))` for `key=value`, `Some(None)` for
/// a bare flag, `None` when the column is `.`, empty or lacks the key.
pub fn info_field<'a>(info: &'a str, key: &str) -> Option<Option<&'a str>> {
    if info == "." || info.is_empty() {
        return None;
    }
    info.rsplit(';')
        .find_map(|field| match field.split_once('=') {
            Some((k, v)) if k == key => Some(Some(v.split_once('=').map_or(v, |(v, _)| v))),
            None if field == key => Some(None),
            _ => None,
        })
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

    /// `my ($key, $value) = split('=', $info)` (BaseVCF4.pm:582): `END=` is the empty
    /// string, `END2=100=5` is `100`, a bare `CHR2` is undef, and `$info_data{$key} =
    /// $value` keeps the last occurrence, a bare key written last leaving undef.
    #[test]
    fn the_info_hash_is_read_as_perl_builds_it() {
        assert_eq!(info_field("END=", "END"), Some(Some("")));
        assert_eq!(info_field("END2=100=5", "END2"), Some(Some("100")));
        assert_eq!(info_field("CHR2", "CHR2"), Some(None));
        assert_eq!(info_field("SVTYPE=DEL", "END"), None);
        assert_eq!(info_field(".", "END"), None);
        assert_eq!(info_field("", "END"), None);
        assert_eq!(info_field("END=100;END=200", "END"), Some(Some("200")));
        assert_eq!(info_field("END=100;END", "END"), Some(None));
        assert_eq!(info_field("END;END=100", "END"), Some(Some("100")));
        assert_eq!(info_value("CHR2", "CHR2"), Some(""));
        assert_eq!(info_value_written("CHR2", "CHR2"), None);
        assert_eq!(info_value_written("CHR2=", "CHR2"), Some(""));
        assert_eq!(info_value("END2=100=5;X=1", "END2"), Some("100"));
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
            assert_eq!(info_field(info, "END"), None, "{info}");
        }
        assert_eq!(info_field("SVTYPE=DEL; END=100", " END"), Some(Some("100")));
    }
}
