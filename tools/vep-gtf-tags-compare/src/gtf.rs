// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The GTF side: the kept attribute codes each `transcript` line's `tag` values imply, by
//! unversioned transcript id.

use crate::input::for_each_record;
use crate::pyfmt::py_strip;
use std::collections::HashMap;
use std::io::{self, Read};

/// The four attribute codes compared, in report order. Bit `i` of a code set stands for
/// `CODES[i]`.
pub const CODES: [&str; 4] = [
    "gencode_basic",
    "gencode_primary",
    "cds_start_NF",
    "cds_end_NF",
];

/// Unversioned transcript id to the set of kept codes present on it.
pub type CodeSets = HashMap<Box<[u8]>, u8>;

/// The code bit a GTF `tag` value stands for: `basic` and `gencode_basic` are one code.
fn tag_bit(value: &str) -> Option<u8> {
    match value {
        "basic" | "gencode_basic" => Some(1),
        "gencode_primary" => Some(2),
        "cds_start_NF" => Some(4),
        "cds_end_NF" => Some(8),
        _ => None,
    }
}

/// The id and code set of one record, when it is a `transcript` line that names an id. A
/// line is a transcript line when it does not start with `#`, has nine tab-separated
/// columns or more and its third is `transcript`. In the ninth column each `;`-separated
/// attribute is stripped and split at its first space; the value is stripped, then freed of
/// the double quotes at its ends. The last `transcript_id` names the line, cut at its first
/// `.`, and an id may be empty.
fn transcript_line(line: &str) -> Option<(&str, u8)> {
    if line.starts_with('#') {
        return None;
    }
    let mut columns = line.split('\t');
    if columns.nth(2) != Some("transcript") {
        return None;
    }
    let attributes = columns.nth(5)?;
    let mut id = None;
    let mut codes = 0u8;
    for attribute in attributes.split(';') {
        let Some((key, value)) = py_strip(attribute).split_once(' ') else {
            continue;
        };
        let value = py_strip(value).trim_matches('"');
        if key == "transcript_id" {
            id = value.split('.').next();
        } else if key == "tag" {
            codes |= tag_bit(value).unwrap_or(0);
        }
    }
    id.map(|id| (id, codes))
}

/// Reads a GTF as text: UTF-8, refused where it is not, with universal newlines. The lines
/// of one transcript id add up to one code set.
pub fn read_gtf(reader: impl Read) -> io::Result<CodeSets> {
    let mut sets = CodeSets::new();
    for_each_record(reader, |record| {
        let line = std::str::from_utf8(record)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if let Some((id, codes)) = transcript_line(line) {
            match sets.get_mut(id.as_bytes()) {
                Some(set) => *set |= codes,
                None => {
                    sets.insert(id.as_bytes().into(), codes);
                }
            }
        }
        Ok(())
    })?;
    Ok(sets)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(attributes: &str) -> String {
        format!("21\tsrc\ttranscript\t1\t9\t.\t+\t.\t{attributes}")
    }

    #[test]
    fn a_transcript_line_yields_its_unversioned_id_and_kept_codes() {
        let text = line(
            r#"gene_id "G1"; transcript_id "T1.4"; tag "basic"; tag "Ensembl_canonical"; tag "cds_end_NF";"#,
        );
        assert_eq!(transcript_line(&text), Some(("T1", 1 | 8)));
        let text = line(r#"transcript_id "T2"; tag "gencode_basic"; tag "gencode_primary""#);
        assert_eq!(transcript_line(&text), Some(("T2", 1 | 2)));
        let text = line(r#"transcript_id "T3"; tag "cds_start_NF"; tag "mRNA_start_NF";"#);
        assert_eq!(transcript_line(&text), Some(("T3", 4)));
        assert_eq!(
            transcript_line(&line(r#"transcript_id "T4";"#)),
            Some(("T4", 0))
        );
    }

    #[test]
    fn other_lines_yield_nothing() {
        let commented = format!("#{}", line(r#"transcript_id "T1";"#));
        assert_eq!(transcript_line(&commented), None);
        let exon = line(r#"transcript_id "T1"; tag "basic";"#).replace("transcript\t", "exon\t");
        assert_eq!(transcript_line(&exon), None);
        let padded = line(r#"transcript_id "T1";"#).replace("transcript\t", "transcript \t");
        assert_eq!(transcript_line(&padded), None);
        assert_eq!(transcript_line("21\tsrc\ttranscript\t1\t9\t.\t+\t."), None);
        assert_eq!(
            transcript_line(&line(r#"gene_id "G1"; tag "basic";"#)),
            None
        );
        assert_eq!(transcript_line(""), None);
    }

    #[test]
    fn attributes_are_read_with_pythons_strip_and_split() {
        // A tenth column is ignored, the last transcript_id wins, and a value keeps what
        // lies inside its quotes.
        let text = format!(
            "{}\textra",
            line(
                "transcript_id \"A\";\u{a0}transcript_id  \"B.2.3\" ;tag basic;tag\u{1f};\
                 tag \" gencode_primary \"; tag \"\"cds_end_NF\"\";TAG \"cds_start_NF\""
            )
        );
        assert_eq!(transcript_line(&text), Some(("B", 1 | 8)));
        assert_eq!(
            transcript_line(&line(r#"transcript_id "";"#)),
            Some(("", 0))
        );
        assert_eq!(
            transcript_line(&line(r#"transcript_id ".7"; tag "basic""#)),
            Some(("", 1))
        );
        assert_eq!(
            transcript_line(&line("transcript_id a b.c d")),
            Some(("a b", 0))
        );
    }

    #[test]
    fn lines_of_one_id_add_up_and_text_must_be_utf8() {
        let text = [
            "#!genome-build TEST\n",
            &line(r#"transcript_id "T1.1"; tag "basic";"#),
            "\r\n",
            &line(r#"transcript_id "T1.2"; tag "cds_start_NF";"#),
            "\r",
            &line(r#"transcript_id "T2";"#),
        ]
        .concat();
        let sets = read_gtf(text.as_bytes()).unwrap();
        assert_eq!(sets.len(), 2);
        assert_eq!(sets[&b"T1"[..]], 1 | 4);
        assert_eq!(sets[&b"T2"[..]], 0);
        assert!(read_gtf(&b"21\tsrc\tgene\t\xff\n"[..]).is_err());
    }
}
