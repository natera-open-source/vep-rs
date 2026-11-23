// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rewrite of one run of whole lines.
//!
//! A line ends at `\n`, at `\r\n` or at a lone `\r`, as universal newlines read it. A line
//! that starts with `#` is written as it stands. Every other line is a row: a row of six
//! cells or more whose Gene and Feature cells are `-` and whose Feature_type cell is
//! `Transcript` has that cell written `-`, and no other byte of a row changes. Every line
//! ending is written `\n`; a last row without an ending gains one, a last `#` line does not.

use memchr::{memchr2, memchr_iter};
use std::ops::AddAssign;

/// The rows of a run and how many of them were rewritten.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    /// Lines that do not start with `#`, an empty line among them.
    pub rows: u64,
    /// Rows whose Feature_type cell was rewritten.
    pub rewritten: u64,
}

impl AddAssign for Counts {
    fn add_assign(&mut self, other: Counts) {
        self.rows += other.rows;
        self.rewritten += other.rewritten;
    }
}

/// The Gene, Feature and Feature_type cells of an intergenic row as fastVEP through v0.4.0
/// writes them.
const INTERGENIC_CELLS: &[u8] = b"-\t-\tTranscript";
/// The offset of the Feature_type cell in them, and its length.
const FEATURE_TYPE_AT: usize = 4;
const FEATURE_TYPE_LEN: usize = 10;

/// The offset in `line` of the Feature_type cell to rewrite, where the line holds one. The
/// cells are the fourth, fifth and sixth, so they start after the line's third tab, and the
/// sixth ends at a tab or at the end of the line.
fn feature_type_to_rewrite(line: &[u8]) -> Option<usize> {
    let third_tab = memchr_iter(b'\t', line).nth(2)?;
    let cells = &line[third_tab + 1..];
    match cells.strip_prefix(INTERGENIC_CELLS)?.first() {
        None | Some(b'\t') => Some(third_tab + 1 + FEATURE_TYPE_AT),
        Some(_) => None,
    }
}

/// Writes the rewrite of `run` into `out`, which is emptied first, and counts its rows.
///
/// `run` starts at the start of a line. A `\r` that is its last byte ends a line by itself,
/// so a run never stops between the `\r` and the `\n` of one line ending.
pub fn rewrite(run: &[u8], out: &mut Vec<u8>) -> Counts {
    out.clear();
    let mut counts = Counts::default();
    // `run[..copied]` is in `out` as rewritten; what follows it is copied in one piece when
    // the next byte to change is met.
    let mut copied = 0;
    let mut start = 0;
    while start < run.len() {
        let (end, next) = match memchr2(b'\n', b'\r', &run[start..]) {
            None => (run.len(), run.len()),
            Some(found) => {
                let end = start + found;
                let crlf = run[end] == b'\r' && run.get(end + 1) == Some(&b'\n');
                (end, end + 1 + usize::from(crlf))
            }
        };
        let line = &run[start..end];
        let row = line.first() != Some(&b'#');
        if row {
            counts.rows += 1;
            if let Some(cell) = feature_type_to_rewrite(line) {
                out.extend_from_slice(&run[copied..start + cell]);
                out.push(b'-');
                copied = start + cell + FEATURE_TYPE_LEN;
                counts.rewritten += 1;
            }
        }
        let ended = next > end;
        if !(ended && run[end] == b'\n') {
            out.extend_from_slice(&run[copied..end]);
            if ended || row {
                out.push(b'\n');
            }
            copied = next;
        }
        start = next;
    }
    out.extend_from_slice(&run[copied..]);
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rewritten(run: &str) -> (String, u64, u64) {
        let mut out = b"left over from the run before".to_vec();
        let counts = rewrite(run.as_bytes(), &mut out);
        (
            String::from_utf8(out).unwrap(),
            counts.rows,
            counts.rewritten,
        )
    }

    fn held(run: &str, expected: &str, rows: u64, changed: u64) {
        assert_eq!(
            rewritten(run),
            (expected.to_string(), rows, changed),
            "{run:?}"
        );
    }

    #[test]
    fn an_intergenic_row_has_its_feature_type_written_as_a_dash() {
        held(
            "v1\t21:100\tG\t-\t-\tTranscript\tintergenic_variant\t-\tIMPACT=MODIFIER\n",
            "v1\t21:100\tG\t-\t-\t-\tintergenic_variant\t-\tIMPACT=MODIFIER\n",
            1,
            1,
        );
    }

    #[test]
    fn a_row_of_six_cells_is_rewritten_at_its_last_cell() {
        held("a\tb\tc\t-\t-\tTranscript\n", "a\tb\tc\t-\t-\t-\n", 1, 1);
        held("\t\t\t-\t-\tTranscript\n", "\t\t\t-\t-\t-\n", 1, 1);
    }

    #[test]
    fn a_row_with_a_gene_a_feature_or_another_feature_type_is_kept() {
        for row in [
            "v1\t21:100\tG\tENSG1\tENST1\tTranscript\tmissense_variant\n",
            "v1\t21:100\tG\tENSG1\t-\tTranscript\tmissense_variant\n",
            "v1\t21:100\tG\t-\tENST1\tTranscript\tmissense_variant\n",
            "v1\t21:100\tG\t-\t-\t-\tintergenic_variant\n",
            "v1\t21:100\tG\t-\t-\tRegulatoryFeature\tregulatory_region_variant\n",
            "v1\t21:100\tG\t-\t-\tTranscripts\tintergenic_variant\n",
            "v1\t21:100\tG\t-\t-\ttranscript\tintergenic_variant\n",
            "v1\t21:100\tG\t-\t-\t Transcript\tintergenic_variant\n",
            "v1\t21:100\tG\t-\t-\tTranscript \tintergenic_variant\n",
            "v1\t21:100\tG\t--\t-\tTranscript\tintergenic_variant\n",
        ] {
            held(row, row, 1, 0);
        }
    }

    #[test]
    fn the_three_cells_count_only_as_the_fourth_fifth_and_sixth() {
        for row in [
            "v1\t21:100\t-\t-\tTranscript\tintergenic_variant\n",
            "v1\t21:100\tG\tx\t-\t-\tTranscript\tintergenic_variant\n",
            "v1\t21:100\tG\tENSG1\tENST1\tTranscript\t-\t-\tTranscript\n",
        ] {
            held(row, row, 1, 0);
        }
    }

    #[test]
    fn a_row_under_six_cells_is_kept_and_counted() {
        held("a\tb\tc\t-\t-\n", "a\tb\tc\t-\t-\n", 1, 0);
        held("a\tb\tc\t-\n", "a\tb\tc\t-\n", 1, 0);
        held("Transcript\n", "Transcript\n", 1, 0);
        held("\n", "\n", 1, 0);
        held("\n\n\n", "\n\n\n", 3, 0);
    }

    #[test]
    fn a_line_that_starts_with_a_hash_is_kept_and_not_counted() {
        held(
            "## fastVEP output\n#a\tb\tc\t-\t-\tTranscript\n",
            "## fastVEP output\n#a\tb\tc\t-\t-\tTranscript\n",
            0,
            0,
        );
        held(
            "a\tb\tc\t-\t-\tTranscript\n#late\nd\te\tf\t-\t-\tTranscript\tx\n",
            "a\tb\tc\t-\t-\t-\n#late\nd\te\tf\t-\t-\t-\tx\n",
            2,
            2,
        );
    }

    #[test]
    fn a_hash_after_the_first_byte_starts_no_comment() {
        held(
            " #a\tb\tc\t-\t-\tTranscript\n",
            " #a\tb\tc\t-\t-\t-\n",
            1,
            1,
        );
        held(
            "\u{feff}#a\tb\tc\t-\t-\tTranscript\n",
            "\u{feff}#a\tb\tc\t-\t-\t-\n",
            1,
            1,
        );
    }

    #[test]
    fn every_line_ending_is_written_as_a_line_feed() {
        held("#h\r\na\tb\r\nc\td\n", "#h\na\tb\nc\td\n", 2, 0);
        held("#h\ra\tb\rc\td\r", "#h\na\tb\nc\td\n", 2, 0);
        held("a\r\r\nb\n\rc\r\n\r\n", "a\n\nb\n\nc\n\n", 6, 0);
        held(
            "a\tb\tc\t-\t-\tTranscript\r\nd\te\tf\t-\t-\tTranscript\rg\n",
            "a\tb\tc\t-\t-\t-\nd\te\tf\t-\t-\t-\ng\n",
            3,
            2,
        );
    }

    #[test]
    fn a_carriage_return_splits_a_line_in_two() {
        held("#h\ra\tb\n", "#h\na\tb\n", 1, 0);
        held(
            "a\tb\tc\t-\t-\tTrans\rcript\n",
            "a\tb\tc\t-\t-\tTrans\ncript\n",
            2,
            0,
        );
        held("x\r#y\n", "x\n#y\n", 1, 0);
    }

    #[test]
    fn a_last_row_without_an_ending_gains_one() {
        held("a\tb", "a\tb\n", 1, 0);
        held("a\tb\tc\t-\t-\tTranscript", "a\tb\tc\t-\t-\t-\n", 1, 1);
        held(
            "a\n\tb\tc\t-\t-\tTranscript\tx",
            "a\n\tb\tc\t-\t-\t-\tx\n",
            2,
            1,
        );
        held("a\nb\tc", "a\nb\tc\n", 2, 0);
    }

    #[test]
    fn a_last_hash_line_without_an_ending_gains_none() {
        held("#h", "#h", 0, 0);
        held("a\n#h", "a\n#h", 1, 0);
        held("a\r\n#h", "a\n#h", 1, 0);
    }

    #[test]
    fn an_empty_run_writes_nothing() {
        held("", "", 0, 0);
    }

    #[test]
    fn bytes_outside_ascii_pass_through() {
        held(
            "v\u{e9}\t21:100\t\u{2192}\t-\t-\tTranscript\t\u{1f9ec}\n",
            "v\u{e9}\t21:100\t\u{2192}\t-\t-\t-\t\u{1f9ec}\n",
            1,
            1,
        );
        let mut out = Vec::new();
        let counts = rewrite(b"a\xff\tb\tc\t-\t-\tTranscript\n", &mut out);
        assert_eq!(out, b"a\xff\tb\tc\t-\t-\t-\n");
        assert_eq!(
            counts,
            Counts {
                rows: 1,
                rewritten: 1
            }
        );
    }

    #[test]
    fn counts_add_up() {
        let mut counts = Counts {
            rows: 3,
            rewritten: 1,
        };
        counts += Counts {
            rows: 4,
            rewritten: 2,
        };
        assert_eq!(
            counts,
            Counts {
                rows: 7,
                rewritten: 3
            }
        );
    }
}
