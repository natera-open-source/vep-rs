// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `gxf_reverse_strand_utr_sequence`: the count after `fsTer` or `extTer` of an
//! HGVSp on a reverse-strand row of a transcript the reference read from a GTF or GFF file.
//!
//! The reference's row states where its transcript comes from. BaseGXF.pm:635 stamps every
//! transcript built from a GTF or GFF file with the file source's short name, which under
//! `--gtf` and `--gff` is the file's last path component (Config.pm:531-543,
//! File.pm:154-158), and OutputFactory.pm:1592 writes the stamp as the row's `SOURCE`. A
//! cache transcript carries a stamp in a merged cache alone, `Ensembl` or `RefSeq`
//! (MergeVEP.pm:103 and :110), and its 3' UTR is stored in the cache, so
//! BaseTranscriptVariation.pm:1101-1113 never reads it again.

use crate::compare::{REFSEQ_SOURCE, SOURCE_FIELDS};
use crate::consts::STRAND_FIELDS;
use crate::entry::Entry;
use crate::rules::{first_field, TER_COUNT_RE};
use crate::shared::blank;

/// The `SOURCE` of a row whose transcript is of the Ensembl set of a merged cache
/// (MergeVEP.pm:103); `REFSEQ_SOURCE` is the other set's (MergeVEP.pm:110).
pub const ENSEMBL_SOURCE: &str = "Ensembl";

/// Whether the row's `SOURCE` names a GTF or GFF file as its transcript's source: the cell
/// is filled and holds neither name a merged cache writes there. A file whose own name is
/// one of those two cannot be told from a merged cache by the row.
pub fn names_a_file_source(e: &Entry) -> bool {
    let source = first_field(e, SOURCE_FIELDS);
    !blank(source) && source != ENSEMBL_SOURCE && source != REFSEQ_SOURCE
}

/// The text of an HGVSp before the count that follows its `fsTer` or `extTer`; `None` for a
/// value that ends in no such count.
fn before_the_stop_count(hgvsp: &str) -> Option<&str> {
    Some(TER_COUNT_RE.captures(hgvsp)?.get(1)?.as_str())
}

/// Whether the two HGVSp values of the reference's row `e` are the class's pair: both are
/// filled and equal up to the count after `fsTer` or `extTer`, the row's `STRAND` is -1 and
/// its `SOURCE` names a file. A row of a cache transcript, a forward-strand row and a pair
/// that differs before the count stay differences.
pub fn file_transcript_stop_count_pair(e: &Entry, expected: &str, actual: &str) -> bool {
    if blank(expected) || blank(actual) {
        return false;
    }
    let reverse_strand = STRAND_FIELDS
        .iter()
        .any(|name| e.fields.get(*name).is_some_and(|held| held == "-1"));
    if !reverse_strand || !names_a_file_source(e) {
        return false;
    }
    match (
        before_the_stop_count(expected),
        before_the_stop_count(actual),
    ) {
        (Some(reference), Some(candidate)) => reference == candidate,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(fields: &[(&str, &str)]) -> Entry {
        Entry {
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            ..Entry::default()
        }
    }

    const REFERENCE: &str = "ENSP1:p.Gly326ArgfsTer?";
    const CANDIDATE: &str = "ENSP1:p.Gly326ArgfsTer31";

    #[test]
    fn a_source_cell_names_a_file_unless_it_is_empty_or_a_merged_caches_set() {
        for (cell, file) in [
            ("genes.gtf.gz", true),
            ("GCF_000001405.40.gff.gz", true),
            ("ensembl", true),
            ("Ensembl", false),
            ("RefSeq", false),
            ("-", false),
            ("", false),
        ] {
            assert_eq!(
                names_a_file_source(&row(&[("SOURCE", cell)])),
                file,
                "{cell}"
            );
            assert_eq!(
                names_a_file_source(&row(&[("source", cell)])),
                file,
                "{cell}"
            );
        }
        assert!(!names_a_file_source(&row(&[("STRAND", "-1")])));
    }

    #[test]
    fn the_pair_is_taken_on_a_reverse_strand_row_of_a_file_transcript_alone() {
        let file_row = row(&[("STRAND", "-1"), ("SOURCE", "genes.gtf.gz")]);
        assert!(file_transcript_stop_count_pair(
            &file_row, REFERENCE, CANDIDATE
        ));
        assert!(file_transcript_stop_count_pair(
            &row(&[("strand", "-1"), ("source", "genes.gtf.gz")]),
            "ENSP1:p.Ter452ArgextTer71",
            "ENSP1:p.Ter452ArgextTer9"
        ));
        for cache_row in [
            row(&[("STRAND", "-1")]),
            row(&[("STRAND", "-1"), ("SOURCE", "-")]),
            row(&[("STRAND", "-1"), ("SOURCE", "Ensembl")]),
            row(&[("STRAND", "-1"), ("SOURCE", "RefSeq")]),
        ] {
            assert!(
                !file_transcript_stop_count_pair(&cache_row, REFERENCE, CANDIDATE),
                "{cache_row:?}"
            );
        }
        let forward = row(&[("STRAND", "1"), ("SOURCE", "genes.gtf.gz")]);
        assert!(!file_transcript_stop_count_pair(
            &forward, REFERENCE, CANDIDATE
        ));
        let unstated = row(&[("SOURCE", "genes.gtf.gz")]);
        assert!(!file_transcript_stop_count_pair(
            &unstated, REFERENCE, CANDIDATE
        ));
    }

    #[test]
    fn a_pair_that_differs_before_the_count_or_lacks_a_value_is_not_taken() {
        let file_row = row(&[("STRAND", "-1"), ("SOURCE", "genes.gtf.gz")]);
        for (expected, actual) in [
            (REFERENCE, "ENSP1:p.Gly326LysfsTer31"),
            (REFERENCE, "ENSP1:p.Gly326Arg"),
            (REFERENCE, "-"),
            ("", CANDIDATE),
            ("ENSP1:p.Ter452ArgextTer71", "ENSP1:p.Ter452ArgfsTer71"),
        ] {
            assert!(
                !file_transcript_stop_count_pair(&file_row, expected, actual),
                "{expected} {actual}"
            );
        }
    }
}
