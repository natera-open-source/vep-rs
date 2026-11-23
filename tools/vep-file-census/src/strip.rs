// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines of a VCF with the `chr` prefix taken off its contig names, as this `awk`
//! program writes them:
//!
//! ```text
//! awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'
//! ```
//!
//! A `##contig=` line loses the `chr` of its first `ID=chr`. Any other line that starts
//! with `#` is written as it is. A record whose first column starts with `chr` loses that
//! prefix and is written with one tab between its columns, the columns being what blanks
//! (spaces and tabs) separate: `awk` builds the record again from its fields once one of
//! them changes. A record whose first column does not start with `chr` is written as it
//! is. Every line written ends in a newline.

use memchr::{memchr, memchr2, memmem};

fn is_blank(byte: u8) -> bool {
    byte == b' ' || byte == b'\t'
}

fn record(line: &[u8], out: &mut Vec<u8>) {
    let start = line.iter().take_while(|byte| is_blank(**byte)).count();
    if !line[start..].starts_with(b"chr") {
        out.extend_from_slice(line);
        return;
    }
    let mut rest = &line[start + 3..];
    let mut first = true;
    loop {
        let end = memchr2(b' ', b'\t', rest).unwrap_or(rest.len());
        if !first {
            out.push(b'\t');
        }
        first = false;
        out.extend_from_slice(&rest[..end]);
        rest = &rest[end..];
        let blanks = rest.iter().take_while(|byte| is_blank(**byte)).count();
        rest = &rest[blanks..];
        if rest.is_empty() {
            return;
        }
    }
}

/// The lines of `data`, a run of whole lines, as the program writes them.
pub fn strip(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() + 1);
    let mut rest = data;
    while !rest.is_empty() {
        let (line, after) = match memchr(b'\n', rest) {
            Some(end) => (&rest[..end], &rest[end + 1..]),
            None => (rest, &rest[rest.len()..]),
        };
        rest = after;
        if line.starts_with(b"##contig=") {
            match memmem::find(line, b"ID=chr") {
                Some(at) => {
                    out.extend_from_slice(&line[..at + 3]);
                    out.extend_from_slice(&line[at + 6..]);
                }
                None => out.extend_from_slice(line),
            }
        } else if line.first() == Some(&b'#') {
            out.extend_from_slice(line);
        } else {
            record(line, &mut out);
        }
        out.push(b'\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stripped(text: &str) -> String {
        String::from_utf8(strip(text.as_bytes())).unwrap()
    }

    #[test]
    fn a_contig_line_loses_the_prefix_of_its_first_id_alone() {
        assert_eq!(
            stripped("##contig=<ID=chr1,length=5,assembly=ID=chr1>\n##contig=<ID=1>\n"),
            "##contig=<ID=1,length=5,assembly=ID=chr1>\n##contig=<ID=1>\n"
        );
    }

    #[test]
    fn other_header_lines_are_written_as_they_are() {
        let text = "##INFO=<ID=chrX  y>\n#CHROM\tPOS  ID\n#chr1\t5\n";
        assert_eq!(stripped(text), text);
    }

    #[test]
    fn a_record_on_a_prefixed_contig_is_built_again_from_its_fields() {
        assert_eq!(
            stripped("chr1\t5\t.\tA B\t\tend \n"),
            "1\t5\t.\tA\tB\tend\n"
        );
        assert_eq!(stripped("  chrX 7\n"), "X\t7\n");
        assert_eq!(stripped("chr\t7\n"), "\t7\n");
        assert_eq!(stripped("chr1"), "1\n");
    }

    #[test]
    fn a_record_on_any_other_contig_keeps_its_bytes() {
        let text = "1\t5\t.\tA B\t\tend \n  \n\n21 chr1\n";
        assert_eq!(stripped(text), text);
        assert_eq!(stripped("X\t7"), "X\t7\n");
    }
}
