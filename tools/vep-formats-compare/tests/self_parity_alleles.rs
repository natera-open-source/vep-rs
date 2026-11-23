// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reverse complement of an inserted allele in made-up rows: an insertion that the
//! reference describes where its shift stops, on either strand, and the `dup` of such a row
//! counted against the exon spans of a made cache.

mod clipped_added_bases_common;
mod compare_common;
mod support;

use std::path::PathBuf;

use clipped_added_bases_common::{plant_shard, transcript_json, MadeTranscript};
use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::{
    tables, HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS as SHIFT_STOP,
};
use vep_formats_compare::consts::{
    HGVSC_FIELDS, HGVSP_FIELDS, HGVS_INTRON_OFFSET_FIELDS, HGVS_OFFSET_FIELDS,
};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::{self, DumpOptions, Value};
use vep_formats_compare::rules::clipped::shift_bound;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::sequence_edges::turned_allele;
use vep_formats_compare::rules::transcript_edges::TranscriptEdges;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};
use vep_formats_compare::self_parity::reverse_complement;

const SHIFT_STOP_FEATURE: &str = "ENST00000000062";
const SHIFT_STOP_REFERENCES: [&str; 3] = ["ins", "dup", "delins"];
const SHIFT_STOP_CANDIDATES: [&str; 7] = [
    "further dup",
    "further ins",
    "ins not turned",
    "shorter dup",
    "del",
    "none",
    "fields without a notation",
];

/// An allele of 999 bases, which `perform_shift` walks at most twice
/// (TranscriptVariationAllele.pm:303-306).
fn shift_stop_allele() -> String {
    "ACG".repeat(333)
}

/// One cell of the made-up grid of `hgvs_insertion_left_where_the_shift_stops`: the
/// reference's notation where its walk stops (an `ins` in an intron, a `dup` of the
/// allele's count of exonic bases, a `delins`); its HGVS_OFFSET (the bound for the allele's
/// length, one more, none); the candidate's row (a `dup` of as many bases or an `ins` of the
/// allele's bases turned, three bases further 3'; that `ins` with its bases not turned; a
/// `dup` of one base fewer; a `del`; nothing; the plugin's two fields without a notation);
/// the candidate's HGVS_OFFSET (the bases its notation stands from the written point, the
/// bound, the other strand's sign, none); the strand; and the plugin's fields (each row's
/// the intron parts of its own notation, the candidate's others, the reference's others, or
/// written on neither row).
#[derive(Clone, Copy, Debug)]
struct Cell {
    reference: &'static str,
    stopped: Option<&'static str>,
    candidate: &'static str,
    moved: Option<&'static str>,
    strand: &'static str,
    plugin: Option<&'static str>,
}

fn cell(
    reference: &'static str,
    stopped: Option<&'static str>,
    candidate: &'static str,
    moved: Option<&'static str>,
    strand: &'static str,
    plugin: Option<&'static str>,
) -> Cell {
    Cell {
        reference,
        stopped,
        candidate,
        moved,
        strand,
        plugin,
    }
}

fn shift_stop_grid() -> Vec<Cell> {
    let mut cells = Vec::new();
    for reference in SHIFT_STOP_REFERENCES {
        for stopped in [Some("the bound"), Some("one more"), None] {
            for candidate in SHIFT_STOP_CANDIDATES {
                for moved in [Some("further"), Some("the bound"), Some("other sign"), None] {
                    for strand in ["1", "-1"] {
                        for plugin in [
                            Some("own"),
                            Some("candidate other"),
                            Some("reference other"),
                            None,
                        ] {
                            cells.push(cell(reference, stopped, candidate, moved, strand, plugin));
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the reference writes an `ins` or
/// a `dup` beside an HGVS_OFFSET that is the bound, and its two fields, where written, are
/// its notation's; against a candidate `dup` of as many bases or `ins` of the allele's
/// bases turned, beside a larger HGVS_OFFSET of the same sign and its own notation's fields,
/// every field that differs; against a candidate that writes nothing, every field that
/// differs, the HGVS_OFFSET where the candidate writes none or the reference's own.
fn shift_stop_set_aside(cell: &Cell, field: &str) -> bool {
    if cell.reference == "delins"
        || cell.stopped != Some("the bound")
        || cell.plugin == Some("reference other")
    {
        return false;
    }
    if cell.candidate == "none"
        || (cell.candidate == "fields without a notation" && cell.plugin.is_none())
    {
        return cell.moved.is_none() || (cell.moved == Some("the bound") && field != "HGVS_OFFSET");
    }
    if !["further dup", "further ins"].contains(&cell.candidate) {
        return false;
    }
    cell.moved == Some("further") && cell.plugin != Some("candidate other")
}

/// The reference's and the candidate's default-format rows of a cell: an insertion of the
/// allele between two bases, described on a transcript without a coding sequence.
fn shift_stop_rows(cell: &Cell, number: usize) -> (String, String) {
    let allele = shift_stop_allele();
    let reverse = cell.strand == "-1";
    let sign = if reverse { "-" } else { "" };
    let oriented = if reverse {
        reverse_complement(&allele)
    } else {
        allele.clone()
    };
    let position = 7_000_000 + 10 * number;
    let columns = format!(
        "cell{number}\t2:{position}-{}\t{allele}\tENSG1\t{SHIFT_STOP_FEATURE}\tTranscript\tnon_coding_transcript_variant\t-\t-\t-\t-\t-\t-",
        position + 1
    );
    let extra = format!("IMPACT=MODIFIER;STRAND={}", cell.strand);
    let fields = |pair: Option<(&str, &str)>| match (cell.plugin, pair) {
        (Some(_), Some((start, end))) => {
            format!(";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}")
        }
        _ => String::new(),
    };
    let (notation, pair) = match cell.reference {
        "ins" => (
            format!("n.2012+1_2012+2ins{}", turned_allele(&oriented, 2, true)),
            ("1", "2"),
        ),
        "dup" => ("n.1014_2012dup".to_string(), ("0", "0")),
        _ => (format!("n.2012_2013delins{oriented}"), ("0", "0")),
    };
    let stopped_text = match cell.stopped {
        Some("the bound") => format!(";HGVS_OFFSET={sign}2"),
        Some(_) => format!(";HGVS_OFFSET={sign}3"),
        None => String::new(),
    };
    let reference_pair = if cell.plugin == Some("reference other") {
        ("7", "7")
    } else {
        pair
    };
    let reference_text = format!(
        "{extra};HGVSc={SHIFT_STOP_FEATURE}.1:{notation}{stopped_text}{}",
        fields(Some(reference_pair))
    );
    let (written, mut own): (Option<String>, Option<(&str, &str)>) = match cell.candidate {
        "further dup" => (Some("n.1017_2015dup".to_string()), Some(("0", "0"))),
        "further ins" => (
            Some(format!(
                "n.2012+4_2012+5ins{}",
                turned_allele(&oriented, 5, true)
            )),
            Some(("4", "5")),
        ),
        "ins not turned" => (
            Some(format!(
                "n.2012+4_2012+5ins{}",
                turned_allele(&oriented, 4, true)
            )),
            Some(("4", "5")),
        ),
        "shorter dup" => (Some("n.1018_2015dup".to_string()), Some(("0", "0"))),
        "del" => (Some("n.1017_2015del".to_string()), Some(("0", "0"))),
        "none" => (None, None),
        _ => (None, Some(("0", "0"))),
    };
    let moved_text = match cell.moved {
        Some("further") => format!(";HGVS_OFFSET={sign}5"),
        Some("the bound") => format!(";HGVS_OFFSET={sign}2"),
        Some(_) => format!(";HGVS_OFFSET={}5", if reverse { "" } else { "-" }),
        None => String::new(),
    };
    if written.is_some() && cell.plugin == Some("candidate other") {
        own = Some(("7", "7"));
    }
    let named = written.map_or_else(String::new, |written| {
        format!(";HGVSc={SHIFT_STOP_FEATURE}.1:{written}")
    });
    let candidate_text = format!("{extra}{named}{moved_text}{}", fields(own));
    (
        format!("{columns}\t{reference_text}"),
        format!("{columns}\t{candidate_text}"),
    )
}

/// The value a default-format row's Extra text holds for `field`.
fn extra_value<'a>(row: &'a str, field: &str) -> Option<&'a str> {
    row.rsplit('\t')
        .next()
        .unwrap_or("")
        .split(';')
        .filter_map(|item| item.split_once('='))
        .filter(|(name, _)| *name == field)
        .map(|(_, value)| value)
        .next_back()
}

/// The cell's two rows with each field of the class that differs between them.
fn shift_stop_asks(cell: &Cell, number: usize) -> Vec<(String, String, &'static str)> {
    let (reference, candidate) = shift_stop_rows(cell, number);
    ["HGVSc", "HGVS_OFFSET"]
        .into_iter()
        .chain(HGVS_INTRON_OFFSET_FIELDS.iter().copied())
        .filter(|field| extra_value(&reference, field) != extra_value(&candidate, field))
        .map(|field| (reference.clone(), candidate.clone(), field))
        .collect()
}

/// What the class says of one of the plugin's fields on a pair of default-format rows.
fn offset_grid_verdict(klass: &str, reference: &str, candidate: &str, field: &str) -> bool {
    let first = |row: &str| {
        parse_default(&default_text(&[row]))
            .expect("the row parses")
            .remove(0)
    };
    let (e, a) = (first(reference), first(candidate));
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(
        klass,
        &e,
        field,
        &field_of(Some(&e), &[field]),
        &field_of(Some(&a), &[field]),
        Some(&a),
        &context,
    )
    .expect("the class reads the rows")
}

/// The rows of an insertion of 1,001 bases, which a walk towards 5' never starts on.
struct Unwalked {
    strand: &'static str,
    stopped: Option<&'static str>,
    moved: Option<&'static str>,
    notation: Option<&'static str>,
    allele: String,
    location: &'static str,
    uploaded: &'static str,
}

fn unwalked_allele() -> String {
    "ACGTA".repeat(200) + "C"
}

fn unwalked(
    strand: &'static str,
    stopped: Option<&'static str>,
    moved: Option<&'static str>,
) -> Unwalked {
    Unwalked {
        strand,
        stopped,
        moved,
        notation: Some("n.1017_2017dup"),
        allele: unwalked_allele(),
        location: "2:7000000-7000001",
        uploaded: "unwalked",
    }
}

impl Unwalked {
    fn rows(&self) -> (String, String) {
        let columns = format!(
            "{}\t{}\t{}\tENSG1\t{SHIFT_STOP_FEATURE}\tTranscript\tnon_coding_transcript_variant\t-\t-\t-\t-\t-\t-",
            self.uploaded, self.location, self.allele
        );
        let extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        let offset = |value: Option<&str>| {
            value.map_or_else(String::new, |value| format!(";HGVS_OFFSET={value}"))
        };
        let reference_text = format!(
            "{extra};HGVSc={SHIFT_STOP_FEATURE}.1:n.1014_2014dup{}",
            offset(self.stopped)
        );
        let named = self.notation.map_or_else(String::new, |notation| {
            format!(";HGVSc={SHIFT_STOP_FEATURE}.1:{notation}")
        });
        let candidate_text = format!("{extra}{named}{}", offset(self.moved));
        (
            format!("{columns}\t{reference_text}"),
            format!("{columns}\t{candidate_text}"),
        )
    }

    fn accepted(&self, field: &str) -> bool {
        let (reference, candidate) = self.rows();
        offset_grid_verdict(SHIFT_STOP, &reference, &candidate, field)
    }
}

/// Over the made-up grid the class sets aside the HGVSc, the HGVS_OFFSET and the
/// plugin's two fields exactly where the rule's statement does. Ensembl VEP reads 1,000
/// bases of flank for the shift (TranscriptVariationAllele.pm:427-437) and walks an
/// insertion at most 1001 less its length times (:303-306), so an insertion whose run of
/// placements goes further is described where the walk stops, beside an HGVS_OFFSET that is
/// that bound. The candidate describes it at the end of the run, or writes nothing where
/// that placement lies past the transcript's last base. Going 5' an allele of 1,001 bases
/// is not walked at all (the limit of :303 is 0, which :304 leaves), so on a reverse-strand
/// row the reference of such an allele carries no HGVS_OFFSET and the class reads it as a
/// walk stopped at once. A reference HGVS_OFFSET of another size, a reference `delins`, a
/// candidate HGVS_OFFSET no larger than the bound or of the other sign, an `ins` whose bases
/// are not the allele's turned by the steps walked, a `dup` of another count of bases, a
/// `del` and fields that are not a notation's own stay differences.
#[test]
fn hgvs_insertion_left_where_the_shift_stops_takes_the_row_at_the_bound_and_nothing_wider() {
    let mut wrong = Vec::new();
    for (number, cell) in shift_stop_grid().iter().enumerate() {
        for (reference, candidate, field) in shift_stop_asks(cell, number) {
            let accepted = offset_grid_verdict(SHIFT_STOP, &reference, &candidate, field);
            if accepted != shift_stop_set_aside(cell, field) {
                wrong.push((*cell, field));
            }
        }
    }
    assert!(wrong.is_empty(), "{wrong:?}");
    let fields: Vec<&str> = [
        HGVSC_FIELDS,
        HGVS_OFFSET_FIELDS,
        HGVS_INTRON_OFFSET_FIELDS,
        HGVSP_FIELDS,
    ]
    .concat();
    assert_eq!(
        tables().field_classes_of(SHIFT_STOP),
        Some(fields.as_slice())
    );
    let bounds = |lengths: &[i64], towards_5prime: bool| -> Vec<i64> {
        lengths
            .iter()
            .map(|length| shift_bound(*length, towards_5prime))
            .collect()
    };
    assert_eq!(
        bounds(&[1, 999, 1000, 1001, 5000], false),
        [1000, 2, 1, 1000, 1000]
    );
    assert_eq!(
        bounds(&[1, 999, 1000, 1001, 1002, 5000], true),
        [1000, 2, 1, 0, 1000, 1000]
    );
    // The walk that never starts: 1,001 bases on a reverse-strand row, the reference without an
    // HGVS_OFFSET. The candidate's notation further 3', or none, is the class's there alone.
    for field in ["HGVSc", "HGVS_OFFSET"] {
        assert!(unwalked("-1", None, Some("-3")).accepted(field));
    }
    // The same duplication one copy further: the HGVS_OFFSET alone differs.
    let one_copy_further = Unwalked {
        notation: Some("n.1014_2014dup"),
        ..unwalked("-1", None, Some("-1001"))
    };
    assert!(one_copy_further.accepted("HGVS_OFFSET"));
    // The most 3' placement past the last base: no notation.
    let no_notation = || Unwalked {
        notation: None,
        ..unwalked("-1", None, None)
    };
    assert!(no_notation().accepted("HGVSc"));
    // On a forward-strand row the same allele is walked 1,000 steps.
    assert!(unwalked("1", Some("1000"), Some("1003")).accepted("HGVSc"));
    // A record that spells an insertion's alleles is read the same way.
    let spelt = Unwalked {
        uploaded: "2_7000000_A/ACG",
        ..unwalked("-1", None, Some("-3"))
    };
    assert!(spelt.accepted("HGVSc"));
    // An allele that keeps bases of its span is never walked and carries no HGVS_OFFSET either:
    // a Location of one, three or five positions, or a record whose reference keeps a base.
    let kept = [
        (Some("2:7000000"), None),
        (Some("2:7000000-7000002"), None),
        (Some("2:7000000-7000004"), None),
        (None, Some("2_7000000_AC/ACGGT")),
    ];
    for (location, uploaded) in kept {
        for (moved, notation) in [(Some("-3"), Some("n.1017_2017dup")), (None, None)] {
            let base = unwalked("-1", None, moved);
            let rows = Unwalked {
                notation,
                location: location.unwrap_or(base.location),
                uploaded: uploaded.unwrap_or(base.uploaded),
                ..base
            };
            assert!(!rows.accepted("HGVSc"), "{location:?} {uploaded:?}");
        }
    }
    // A reference notation of fewer bases than the allele holds is no insertion's own: a dup of
    // 999 exonic bases, an ins of 999 bases.
    let turned = reverse_complement(&unwalked_allele());
    let (reference, candidate) = no_notation().rows();
    for notation in [
        "n.1016_2014dup".to_string(),
        format!("n.1014_1015ins{}", &turned[..999]),
    ] {
        let reference = reference.replace("n.1014_2014dup", &notation);
        assert!(!offset_grid_verdict(
            SHIFT_STOP, &reference, &candidate, "HGVSc"
        ));
    }
    let whole = reference.replace("n.1014_2014dup", &format!("n.1014_1015ins{turned}"));
    assert!(offset_grid_verdict(SHIFT_STOP, &whole, &candidate, "HGVSc"));
    let allele = unwalked_allele();
    let refused = [
        unwalked("1", None, Some("3")),
        unwalked("-1", None, Some("3")),
        unwalked("-1", None, None),
        unwalked("-1", Some("-1000"), Some("-1003")),
        Unwalked {
            notation: Some("n.1017_2016dup"),
            allele: allele[..allele.len() - 1].to_string(),
            ..unwalked("-1", None, Some("-3"))
        },
        Unwalked {
            notation: Some("n.1017_2018dup"),
            allele: format!("{allele}A"),
            ..unwalked("-1", None, Some("-3"))
        },
        Unwalked {
            notation: Some("n.1017_2017del"),
            ..unwalked("-1", None, Some("-3"))
        },
    ];
    for rows in &refused {
        assert!(!rows.accepted("HGVSc"));
    }
    // A deletion's row is not read.
    let walked = cell(
        "dup",
        Some("the bound"),
        "further dup",
        Some("further"),
        "1",
        Some("own"),
    );
    let (reference, candidate) = shift_stop_rows(&walked, 0);
    assert!(offset_grid_verdict(
        SHIFT_STOP, &reference, &candidate, "HGVSc"
    ));
    let deleted = |row: &str| row.replace(&format!("\t{}\t", shift_stop_allele()), "\t-\t");
    assert!(!offset_grid_verdict(
        SHIFT_STOP,
        &deleted(&reference),
        &deleted(&candidate),
        "HGVSc"
    ));
}

/// Through the comparison itself the class sets aside the differing fields of a row at
/// the bound, counted once, with the plugin's fields and without them, and the row stays a
/// difference without the class.
#[test]
fn hgvs_insertion_left_where_the_shift_stops_counts_a_row_once_through_the_comparison() {
    let cells = [
        cell(
            "ins",
            Some("the bound"),
            "further ins",
            Some("further"),
            "1",
            Some("own"),
        ),
        cell(
            "dup",
            Some("the bound"),
            "further dup",
            Some("further"),
            "-1",
            None,
        ),
        cell("ins", Some("the bound"), "none", None, "1", Some("own")),
        cell(
            "ins",
            Some("the bound"),
            "none",
            Some("the bound"),
            "-1",
            None,
        ),
    ];
    for cell in &cells {
        let (reference, candidate) = shift_stop_rows(cell, 0);
        let (reference, candidate) = (default_text(&[reference]), default_text(&[candidate]));
        let (report, documented) = compare_default(&reference, &candidate, &[SHIFT_STOP]);
        assert_eq!(verdict(&report), "green");
        assert_eq!(counts(&documented.applied), count_map(&[(SHIFT_STOP, 1)]));
        assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    }
}

const COPY_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];
const COPY_TERM: &str = "non_coding_transcript_variant";

/// A reverse-strand transcript whose exon at the point 7000000 | 7000001 has an intron of
/// two bases 5' of it and one of 100 bases 3' of it.
const AROUND_INTRONS: &str = "ENST00000000062";
/// A forward-strand transcript with an intron of 100 bases 3' of the point's exon.
const FORWARD_INTRON: &str = "ENST00000000063";
/// A reverse-strand transcript of one exon.
const ONE_EXON: &str = "ENST00000000064";
/// A reverse-strand transcript whose intron 5' of the point's exon holds the base 1,001
/// bases 5' of the point; the base 1,001 bases 3' of it lies past the last exon.
const INTRON_START: &str = "ENST00000000065";
/// A transcript the made cache does not hold.
const NOT_IN_THE_CACHE: &str = "ENST00000000066";
/// A reverse-strand transcript whose intron 5' of the point's exon is 104 bases long and
/// holds the base 1,001 bases 5' of the point as its fifth base.
const INTRON_AFTER_AN_EXON: &str = "ENST00000000067";

const COPY_TRANSCRIPTS: [MadeTranscript; 5] = [
    (
        INTRON_AFTER_AN_EXON,
        -1,
        &[(6999800, 7000901), (7001006, 7003000)],
    ),
    (
        AROUND_INTRONS,
        -1,
        &[(6990000, 6999699), (6999800, 7000500), (7000503, 7003000)],
    ),
    (FORWARD_INTRON, 1, &[(6990000, 7001000), (7001101, 7005000)]),
    (ONE_EXON, -1, &[(6990000, 7005000)]),
    (INTRON_START, -1, &[(6999800, 7000500), (7002001, 7003000)]),
];

/// The side tables of a comparison given a made cache of `COPY_TRANSCRIPTS`.
fn copy_side(name: &str) -> (Side, PathBuf) {
    let cache = scratch(name);
    let transcripts: Vec<String> = COPY_TRANSCRIPTS
        .iter()
        .map(|(stable_id, strand, exons)| transcript_json(stable_id, *strand, exons))
        .collect();
    plant_shard(
        &cache,
        "2",
        "1-10000000.json",
        &format!("[{}]", transcripts.join(", ")),
    );
    let side = Side {
        transcript_edges: Some(TranscriptEdges::new(&cache, None).expect("the cache reads")),
        ..Side::default()
    };
    (side, cache)
}

/// One made row of the insertion of 1,001 bases between 2:7000000 and 2:7000001, written
/// under an identifier: its transcript, its STRAND, the coordinates of its `dup` on that
/// transcript or none for no HGVSc, and its HGVS_OFFSET or none.
#[derive(Debug, Clone)]
struct CopyRow {
    feature: &'static str,
    strand: i64,
    coordinates: Option<String>,
    offset: Option<i64>,
}

fn copy_row(
    feature: &'static str,
    strand: i64,
    coordinates: Option<&str>,
    offset: Option<i64>,
) -> CopyRow {
    CopyRow {
        feature,
        strand,
        coordinates: coordinates.map(str::to_string),
        offset,
    }
}

impl CopyRow {
    /// The row's HGVSc: its coordinates in the `n.` numbering, or as written where they
    /// carry the `c.` numbering.
    fn hgvsc(&self) -> Option<String> {
        self.coordinates.as_ref().map(|coordinates| {
            let numbering = if coordinates.starts_with("c.") {
                ""
            } else {
                "n."
            };
            format!("{}.1:{numbering}{coordinates}dup", self.feature)
        })
    }

    /// The columns of the row in the text formats up to the Extra text.
    fn text_columns(&self) -> String {
        format!(
            "ins1\t2:7000000-7000001\t{}\tENSG1\t{}\tTranscript\t{COPY_TERM}\t-\t-\t-\t-\t-\t-",
            unwalked_allele(),
            self.feature
        )
    }

    /// A default-format output of the row alone.
    fn default_output(&self) -> String {
        let mut extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        if let Some(hgvsc) = self.hgvsc() {
            extra.push_str(&format!(";HGVSc={hgvsc}"));
        }
        if let Some(offset) = self.offset {
            extra.push_str(&format!(";HGVS_OFFSET={offset}"));
        }
        default_text(&[format!("{}\t{extra}", self.text_columns())])
    }

    /// The row's entry as the format `fmt` writes and the tool parses it.
    fn entry(&self, fmt: &str) -> Entry {
        let (allele, feature, hgvsc) = (unwalked_allele(), self.feature, self.hgvsc());
        let record = format!("2\t7000000\tins1\tG\tG{allele}\t.\t.\t.");
        let columns = self.text_columns();
        let offset = |absent: &str| {
            self.offset
                .map_or(absent.to_string(), |offset| offset.to_string())
        };
        let entries = match fmt {
            "default" => parse_default(&self.default_output()),
            "tab" => {
                let header =
                    DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSc\tHGVS_OFFSET");
                let row = format!(
                    "{columns}\tMODIFIER\t{}\t{}\t{}",
                    self.strand,
                    hgvsc.as_deref().unwrap_or("-"),
                    offset("-")
                );
                parse_tab(&format!("{header}\n{row}\n"))
            }
            "vcf" => {
                let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|STRAND|HGVSc|HGVS_OFFSET\">\n";
                let line = format!(
                    "{}CSQ={allele}|{COPY_TERM}|{feature}|Transcript|{}|{}|{}\n",
                    &record[..record.len() - 1],
                    self.strand,
                    hgvsc.as_deref().unwrap_or(""),
                    offset("")
                );
                parse_vcf(&format!("{header}{line}"), &[record])
            }
            _ => {
                let text = |value: &str| Value::Str(value.to_string());
                let mut item = IndexMap::from([
                    ("variant_allele".to_string(), text(&allele)),
                    ("transcript_id".to_string(), text(feature)),
                    (
                        "consequence_terms".to_string(),
                        Value::List(vec![text(COPY_TERM)]),
                    ),
                    ("strand".to_string(), Value::Int(self.strand)),
                ]);
                if let Some(hgvsc) = &hgvsc {
                    item.insert("hgvsc".to_string(), text(hgvsc));
                }
                if let Some(offset) = self.offset {
                    item.insert("hgvs_offset".to_string(), Value::Int(offset));
                }
                let object = Value::Dict(IndexMap::from([
                    ("input".to_string(), text(&record)),
                    ("id".to_string(), text("ins1")),
                    ("seq_region_name".to_string(), text("2")),
                    ("end".to_string(), Value::Int(7000000)),
                    ("strand".to_string(), Value::Int(1)),
                    ("start".to_string(), Value::Int(7000001)),
                    ("allele_string".to_string(), text(&format!("-/{allele}"))),
                    ("most_severe_consequence".to_string(), text(COPY_TERM)),
                    (
                        "transcript_consequences".to_string(),
                        Value::List(vec![Value::Dict(item)]),
                    ),
                ]));
                let line = pyjson::dumps(&object, &DumpOptions::default()) + "\n";
                parse_json(&line, &[record])
            }
        };
        entries
            .expect("the row parses")
            .into_iter()
            .find(|entry| entry.feature == feature)
            .expect("the row's entry")
    }

    /// The row's entry in each format of `COPY_FORMATS`.
    fn entries(&self) -> Vec<Entry> {
        COPY_FORMATS.iter().map(|fmt| self.entry(fmt)).collect()
    }
}

/// What the class says of each differing field of HGVSc and HGVS_OFFSET on a made pair of
/// rows, the same pair read in the default, tab, VCF and JSON formats, beside `side`.
fn copy_verdicts(reference: &[Entry], candidate: &[Entry], side: &Side) -> Vec<Vec<bool>> {
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side,
    };
    reference
        .iter()
        .zip(candidate)
        .map(|(e, a)| {
            ["HGVSc", "HGVS_OFFSET"]
                .into_iter()
                .filter(|field| field_of(Some(e), &[field]) != field_of(Some(a), &[field]))
                .map(|field| {
                    field_class_accepts(
                        SHIFT_STOP,
                        e,
                        field,
                        &field_of(Some(e), &[field]),
                        &field_of(Some(a), &[field]),
                        Some(a),
                        &context,
                    )
                    .expect("the class reads the rows")
                })
                .collect()
        })
        .collect()
}

fn in_each_format(expected: &[bool]) -> Vec<Vec<bool>> {
    vec![expected.to_vec(); COPY_FORMATS.len()]
}

/// The rows a labelled list holds, each put to the class beside `side`, against what the
/// list says of them in every format.
fn assert_copies(rows: &[(&str, CopyRow, CopyRow, &[bool])], side: &Side) {
    for (label, reference, candidate, expected) in rows {
        assert_eq!(
            copy_verdicts(&reference.entries(), &candidate.entries(), side),
            in_each_format(expected),
            "{label}"
        );
    }
}

/// A reverse-strand row of the walk that never starts whose `dup` on `feature` starts at
/// the coordinate 4000 and counts `count` bases, and the candidate's row without a notation.
fn counted(feature: &'static str, count: i64) -> (CopyRow, CopyRow) {
    let coordinates = if count == 1 {
        "4000".to_string()
    } else {
        format!("4000_{}", 3999 + count)
    };
    (
        copy_row(feature, -1, Some(&coordinates), None),
        copy_row(feature, -1, None, None),
    )
}

/// Each end of a `dup` is mapped to a coordinate alone
/// (TranscriptVariationAllele.pm:1524-1525), so the two coordinates are those of the first
/// and last of as many consecutive bases of the chromosome as the allele holds, the run
/// that follows the notation's point or the run that precedes it. Given the candidate's
/// cache, the class places the two runs by the row's Location and the HGVS_OFFSET beside
/// the notation and
/// takes a `dup` whose coordinates are a run's: 1,001 bases of one exon, also from the 5'
/// UTR into the coding sequence, where `c.-1` and `c.1` are neighbours; 999 around an
/// intron of two bases, the count an allele that keeps the two bases of its span adds; 901
/// around an intron of 100 bases; and a `dup` that starts in an intron, 501 bases before an
/// exon's first base. An intron base numbered from the stop codon's last base is written
/// without the `+` of its part (:2854-2859), so `c.*5` is also the fifth base of the intron
/// after an exon's last base, and a `dup` from it is a run's. The candidate's `dup` further
/// 3' is counted the same way on a reverse-strand row and on a forward-strand row.
#[test]
fn the_shift_stop_class_counts_a_copy_less_the_introns_between_its_coordinates() {
    let (side, cache) = copy_side("shift-stop-copy-counted");
    let nothing = |feature: &'static str| copy_row(feature, -1, None, None);
    let rows: [(&str, CopyRow, CopyRow, &[bool]); 8] = [
        (
            "a start in the intron after the stop codon's last base",
            copy_row(INTRON_AFTER_AN_EXON, -1, Some("c.*5_*901"), None),
            nothing(INTRON_AFTER_AN_EXON),
            &[true],
        ),
        (
            "one exon",
            copy_row(ONE_EXON, -1, Some("4000_5000"), None),
            nothing(ONE_EXON),
            &[true],
        ),
        (
            "one exon, from the 5' UTR into the coding sequence",
            copy_row(ONE_EXON, -1, Some("c.-500_501"), None),
            nothing(ONE_EXON),
            &[true],
        ),
        (
            "two fewer around an intron of two bases",
            copy_row(AROUND_INTRONS, -1, Some("2000_2998"), None),
            nothing(AROUND_INTRONS),
            &[true],
        ),
        (
            "the run on the other side of the point, around an intron of 100 bases",
            copy_row(AROUND_INTRONS, -1, Some("2999_3899"), None),
            nothing(AROUND_INTRONS),
            &[true],
        ),
        (
            "a candidate around an intron, reverse strand",
            copy_row(AROUND_INTRONS, -1, Some("2999_3899"), None),
            copy_row(AROUND_INTRONS, -1, Some("3002_3902"), Some(-3)),
            &[true, true],
        ),
        (
            "a candidate around an intron, forward strand",
            copy_row(FORWARD_INTRON, 1, Some("10001_11001"), Some(1000)),
            copy_row(FORWARD_INTRON, 1, Some("10801_11701"), Some(1800)),
            &[true, true],
        ),
        (
            "a start in an intron",
            copy_row(INTRON_START, -1, Some("1001-501_1500"), None),
            nothing(INTRON_START),
            &[true],
        ),
    ];
    assert_copies(&rows, &side);
    std::fs::remove_dir_all(cache).expect("the made cache is removed");
}

/// Given the candidate's cache, a `dup` whose coordinates are no run's stays a difference.
/// On a transcript of one exon that is every count but the allele's own 1,001: each count
/// from 2 to 998, 999 and 1,000. Around an intron it is a count one off the run's, and the
/// allele's own count where each run holds an intron. So do a candidate `dup` one base off
/// its run's count, a start whose intron part or whose exon coordinate is not the run's,
/// written with its `+` or `-` or as `c.*5`, a coordinate that numbers no base, a `dup`
/// with one coordinate in the 3' UTR and
/// the other before it, whose distance the notation does not carry, a row whose transcript
/// the cache does not hold, and a row whose STRAND is not the cache's.
#[test]
fn the_shift_stop_class_refuses_a_count_no_intron_explains() {
    let (side, cache) = copy_side("shift-stop-copy-refused");
    for count in 2..=1000 {
        let (reference, candidate) = counted(ONE_EXON, count);
        assert_eq!(
            copy_verdicts(&reference.entries(), &candidate.entries(), &side),
            in_each_format(&[false]),
            "{count}"
        );
    }
    let nothing = |feature: &'static str| copy_row(feature, -1, None, None);
    let reference =
        |feature: &'static str, coordinates: &str| copy_row(feature, -1, Some(coordinates), None);
    let mut rows: Vec<(&str, CopyRow, CopyRow, &[bool])> = Vec::new();
    for coordinates in [
        "2000_2997",
        "2000_2999",
        "2000_3000",
        "2999_3898",
        "2999_3900",
    ] {
        rows.push((
            "around an intron",
            reference(AROUND_INTRONS, coordinates),
            nothing(AROUND_INTRONS),
            &[false],
        ));
    }
    for coordinates in ["1001-500_1500", "1001-501_1499", "1001-501_1501"] {
        rows.push((
            "a start in an intron",
            reference(INTRON_START, coordinates),
            nothing(INTRON_START),
            &[false],
        ));
    }
    for coordinates in ["c.*6_*901", "c.*5_*900", "c.*5_*902", "c.*4_*900"] {
        rows.push((
            "a start in the intron after the stop codon's last base",
            reference(INTRON_AFTER_AN_EXON, coordinates),
            nothing(INTRON_AFTER_AN_EXON),
            &[false],
        ));
    }
    for coordinates in ["c.-500_500", "c.-500_502", "c.900_*100"] {
        rows.push((
            "one exon, a count across the start codon one off, and the 3' UTR after a coding base",
            reference(ONE_EXON, coordinates),
            nothing(ONE_EXON),
            &[false],
        ));
    }
    rows.push((
        "a coordinate before the stop codon's last base that carries a second signed part",
        reference(ONE_EXON, "c.*-5+3_*1000"),
        nothing(ONE_EXON),
        &[false],
    ));
    for coordinates in ["3002_3901", "3002_3903"] {
        rows.push((
            "a candidate one off its run, reverse strand",
            reference(AROUND_INTRONS, "2999_3899"),
            copy_row(AROUND_INTRONS, -1, Some(coordinates), Some(-3)),
            &[false, false],
        ));
    }
    for coordinates in ["10801_11700", "10801_11702"] {
        rows.push((
            "a candidate one off its run, forward strand",
            copy_row(FORWARD_INTRON, 1, Some("10001_11001"), Some(1000)),
            copy_row(FORWARD_INTRON, 1, Some(coordinates), Some(1800)),
            &[false, false],
        ));
    }
    rows.push((
        "a transcript the cache does not hold",
        reference(NOT_IN_THE_CACHE, "4000_5000"),
        nothing(NOT_IN_THE_CACHE),
        &[false],
    ));
    rows.push((
        "a STRAND that is not the cache's",
        reference(FORWARD_INTRON, "4000_5000"),
        nothing(FORWARD_INTRON),
        &[false],
    ));
    assert_copies(&rows, &side);
    std::fs::remove_dir_all(cache).expect("the made cache is removed");
}

/// A Location whose positions 64 bits do not hold names no base of a cache's exon spans, so
/// given the cache a `dup` the class places by its Location stays a difference there: the
/// candidate's `dup` of another count than the allele's on a walked row, and the reference's
/// `dup` on the row of the walk that never starts.
#[test]
fn the_shift_stop_class_places_no_copy_by_a_location_past_64_bits() {
    let (side, cache) = copy_side("shift-stop-copy-long-location");
    let far = "2:170141183460469231731687303715884105726-170141183460469231731687303715884105727";
    let entry = |row: String| {
        parse_default(&default_text(&[row]))
            .expect("the row parses")
            .remove(0)
    };
    let walked = cell(
        "dup",
        Some("the bound"),
        "shorter dup",
        Some("further"),
        "1",
        None,
    );
    let (reference, candidate) = shift_stop_rows(&walked, 0);
    let moved = |row: String| row.replace("2:7000000-7000001", far);
    let never_started = Unwalked {
        notation: None,
        location: far,
        ..unwalked("-1", None, None)
    };
    let (unwalked_reference, unwalked_candidate) = never_started.rows();
    let rows = [
        (moved(reference), moved(candidate), vec![false, false]),
        (unwalked_reference, unwalked_candidate, vec![false]),
    ];
    for (reference, candidate, expected) in rows {
        assert!(reference.contains(far) && candidate.contains(far));
        assert_eq!(
            copy_verdicts(&[entry(reference)], &[entry(candidate)], &side),
            [expected]
        );
    }
    std::fs::remove_dir_all(cache).expect("the made cache is removed");
}

/// The reading without the cache: the class has no exon spans to place a run by and reads a
/// `dup` of the walk that never starts by its count alone. Two exonic coordinates of one
/// region are an insertion's own from 2 bases up to the allele's 1,001, whatever intron
/// would explain the count, but for 999, the count an allele that keeps the two bases of
/// its span adds; a single coordinate and 1,002 stay differences. A `dup` with a coordinate
/// in an intron, or with its coordinates in two regions, is not counted, and a candidate
/// `dup` of two exonic coordinates names as many bases as the allele holds.
#[test]
fn without_the_cache_the_shift_stop_class_reads_a_copy_by_its_count_alone() {
    let side = Side::default();
    for count in 1..=1002 {
        let (reference, candidate) = counted(ONE_EXON, count);
        let taken = (2..=1001).contains(&count) && count != 999;
        assert_eq!(
            copy_verdicts(&reference.entries(), &candidate.entries(), &side),
            in_each_format(&[taken]),
            "{count}"
        );
    }
    let nothing = |feature: &'static str| copy_row(feature, -1, None, None);
    let mut rows: Vec<(&str, CopyRow, CopyRow, &[bool])> = vec![
        (
            "two fewer around an intron of two bases",
            copy_row(AROUND_INTRONS, -1, Some("2000_2998"), None),
            nothing(AROUND_INTRONS),
            &[false],
        ),
        (
            "a candidate around an intron, reverse strand",
            copy_row(AROUND_INTRONS, -1, Some("2999_3899"), None),
            copy_row(AROUND_INTRONS, -1, Some("3002_3902"), Some(-3)),
            &[false, false],
        ),
        (
            "a candidate around an intron, forward strand",
            copy_row(FORWARD_INTRON, 1, Some("10001_11001"), Some(1000)),
            copy_row(FORWARD_INTRON, 1, Some("10801_11701"), Some(1800)),
            &[false, false],
        ),
        (
            "a transcript no cache is asked for",
            copy_row(NOT_IN_THE_CACHE, -1, Some("4000_5000"), None),
            nothing(NOT_IN_THE_CACHE),
            &[true],
        ),
    ];
    for coordinates in ["1001-501_1500", "1001-500_1500", "1001-501_1499"] {
        rows.push((
            "a start in an intron",
            copy_row(INTRON_START, -1, Some(coordinates), None),
            nothing(INTRON_START),
            &[true],
        ));
    }
    for coordinates in ["c.-500_501", "c.-500_500", "c.900_*100"] {
        rows.push((
            "coordinates in two regions",
            copy_row(ONE_EXON, -1, Some(coordinates), None),
            nothing(ONE_EXON),
            &[true],
        ));
    }
    assert_copies(&rows, &side);
}

/// Through the command line the class counts by the cache `--json-cache` names. The row
/// whose `dup` counts two fewer around an intron of two bases is set aside with the cache
/// and stays a difference without it; a `dup` of 500 bases on a transcript of one exon is
/// set aside without the cache and stays a difference with it.
#[test]
fn the_shift_stop_class_reads_the_cache_the_command_line_names() {
    let (_, cache) = copy_side("shift-stop-copy-command-line");
    let ended = |feature: &'static str, coordinates: &str, named: bool| {
        let path = |file: &str| cache.join(file).display().to_string();
        let reference = copy_row(feature, -1, Some(coordinates), None).default_output();
        let candidate = copy_row(feature, -1, None, None).default_output();
        std::fs::write(path("reference.txt"), reference).expect("the reference is written");
        std::fs::write(path("candidate.txt"), candidate).expect("the candidate is written");
        let _ = std::fs::remove_file(path("report.json"));
        let mut args: Vec<String> = [
            "--format",
            "default",
            "--expected",
            &path("reference.txt"),
            "--actual",
            &path("candidate.txt"),
            "--masks",
            SHIFT_STOP,
            "--report",
            &path("report.json"),
        ]
        .map(str::to_string)
        .to_vec();
        if named {
            args.extend(["--json-cache".to_string(), cache.display().to_string()]);
        }
        let ran = support::run(&cache, &args);
        let written = support::report(&cache.join("report.json"));
        let applied: Vec<(String, i64)> = support::member(&written, &["masks", "applied"])
            .as_dict()
            .expect("the applied classes")
            .iter()
            .map(|(klass, count)| match count {
                Value::Int(count) => (klass.clone(), *count),
                other => panic!("{klass} counts {other:?}"),
            })
            .collect();
        (ran.code, applied, ran.err())
    };
    let set_aside = (0, vec![(SHIFT_STOP.to_string(), 1)], String::new());
    let difference = (1, Vec::new(), String::new());
    assert_eq!(ended(AROUND_INTRONS, "2000_2998", true), set_aside);
    assert_eq!(ended(AROUND_INTRONS, "2000_2998", false), difference);
    assert_eq!(ended(ONE_EXON, "4000_4499", true), difference);
    assert_eq!(ended(ONE_EXON, "4000_4499", false), set_aside);
    std::fs::remove_dir_all(cache).expect("the made cache is removed");
}
