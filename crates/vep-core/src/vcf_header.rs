// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `##fileformat` test ensembl-io's VCF parser makes on every `##` header
//! line it reads (`BaseVCF4.pm:93-131`), shared by every vep-rs reader of a VCF
//! that Ensembl VEP opens through that parser: its input, and through the
//! parser's tabix form (`VCF4Tabix.pm:55-58`) a `--custom` source, a plugin's
//! data file and the cache dumper's frequency file.

use std::fmt;

/// Why Ensembl VEP ends a run on a VCF header line, in its words, or `None`
/// for a line it reads on. Its ensembl-io parser reads every `##key=value`
/// header line (BaseVCF4.pm:93-94: `##`, blanks, a word, `=`, a non-empty value)
/// and tests the `fileformat` key's value (`:116-130`): the first
/// `<digits>.<digits>` in it is the version, and a major number other than the
/// parser's own 4, compared as numbers, ends the run (`:125`), as does a value
/// with no such pair (`:129`); the minor number is not read (`:126`), so `VCFv4.5`
/// passes, and a line with no value, another key or `FILEFORMAT` carries no
/// version. The check runs on every such line of the file, and in Ensembl VEP
/// before the output opens.
pub fn fileformat_refusal(line: &[u8]) -> Option<String> {
    let rest = line.strip_prefix(b"##")?;
    let key_start = rest
        .iter()
        .position(|&b| !is_ascii_space(b))
        .unwrap_or(rest.len());
    let value = rest[key_start..].strip_prefix(b"fileformat=")?;
    if value.is_empty() {
        return None;
    }
    let Some((major, minor)) = first_version(value) else {
        return Some(
            "The script can't read the VCF file format version of 'fileformat'".to_string(),
        );
    };
    let significant = major.iter().position(|&b| b != b'0').unwrap_or(major.len());
    if &major[significant..] == b"4" {
        return None;
    }
    Some(incompatible_version(
        String::from_utf8_lossy(major),
        String::from_utf8_lossy(minor),
    ))
}

/// Ensembl VEP's words for a `##fileformat` version whose major number is not
/// the parser's 4 (`BaseVCF4.pm:125`), the version printed as the line gave it.
pub fn incompatible_version(major: impl fmt::Display, minor: impl fmt::Display) -> String {
    format!(
        "The VCF file format version {major}.{minor} is not compatible with the parser version (VCF v4.2)"
    )
}

/// `char::is_whitespace` over the ASCII range, Perl's `\s` on the bytes of a line.
fn is_ascii_space(b: u8) -> bool {
    matches!(b, b'\t' | b'\n' | 0x0B | 0x0C | b'\r' | b' ')
}

/// The first `<digits>.<digits>` of `value`, Ensembl VEP's `/(\d+)\.(\d+)/`: the
/// leftmost run of digits a `.` and a digit follow, with the digits after the `.`.
fn first_version(value: &[u8]) -> Option<(&[u8], &[u8])> {
    let digits_end = |from: usize| {
        value[from..]
            .iter()
            .position(|b| !b.is_ascii_digit())
            .map_or(value.len(), |n| from + n)
    };
    let mut i = 0;
    while i < value.len() {
        if !value[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let major_end = digits_end(i);
        if value.get(major_end) == Some(&b'.')
            && value.get(major_end + 1).is_some_and(u8::is_ascii_digit)
        {
            let minor_start = major_end + 1;
            let minor_end = digits_end(minor_start);
            return Some((&value[i..major_end], &value[minor_start..minor_end]));
        }
        i = major_end;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    const INCOMPATIBLE_3_3: &str =
        "The VCF file format version 3.3 is not compatible with the parser version (VCF v4.2)";
    const UNREADABLE: &str = "The script can't read the VCF file format version of 'fileformat'";

    /// Ensembl VEP 116.3, measured on each line as the first header line of a
    /// `--custom` VCF source: `VCFv3.3`, `VCFv5.0` and `##  fileformat=VCFv3.3`
    /// end the run with the version's text, `VCFv4` with the unreadable-version
    /// text; `VCFv4.2`, `VCFv4.5`, `VCFv04.2`, an empty value and `FILEFORMAT`
    /// are read on, as is a file with no such line.
    #[test]
    fn fileformat_refusal_follows_the_parser_on_every_measured_shape() {
        assert_eq!(
            fileformat_refusal(b"##fileformat=VCFv3.3").as_deref(),
            Some(INCOMPATIBLE_3_3)
        );
        assert_eq!(
            fileformat_refusal(b"##fileformat=VCFv5.0").as_deref(),
            Some("The VCF file format version 5.0 is not compatible with the parser version (VCF v4.2)")
        );
        assert_eq!(
            fileformat_refusal(b"##  fileformat=VCFv3.3").as_deref(),
            Some(INCOMPATIBLE_3_3)
        );
        assert_eq!(
            fileformat_refusal(b"##fileformat=VCFv4").as_deref(),
            Some(UNREADABLE)
        );
        for read_on in [
            "##fileformat=VCFv4.2",
            "##fileformat=VCFv4.5",
            "##fileformat=VCFv04.2",
            "##fileformat=",
            "##FILEFORMAT=VCFv3.3",
            "##source=VCFv3.3",
            "#CHROM\tPOS",
            "21\t5\t.\tA\tG\t.\t.\t.",
        ] {
            assert_eq!(fileformat_refusal(read_on.as_bytes()), None, "{read_on}");
        }
    }

    /// The version is the leftmost `<digits>.<digits>`, whatever precedes it
    /// (`/(\d+)\.(\d+)/` is unanchored): `v1.0 of VCFv4.2` reads 1.0, a bare
    /// `4.` or `.2` reads nothing.
    #[test]
    fn first_version_is_the_leftmost_dotted_number_pair() {
        assert_eq!(
            first_version(b"version 1.0 of VCFv4.2"),
            Some((&b"1"[..], &b"0"[..]))
        );
        assert_eq!(first_version(b"VCFv4.2"), Some((&b"4"[..], &b"2"[..])));
        assert_eq!(first_version(b"VCFv10.1"), Some((&b"10"[..], &b"1"[..])));
        assert_eq!(first_version(b"VCFv4."), None);
        assert_eq!(first_version(b"VCFv.2"), None);
        assert_eq!(first_version(b"VCF"), None);
    }
}
