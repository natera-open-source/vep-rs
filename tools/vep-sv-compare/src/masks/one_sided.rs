// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The one-sided classes of the 116.2 reference that are read off a record's own text: the
//! dot-notation breakend, the missing ALT, the lone star allele, the coordinate that is not
//! an integer and the POS that is not a number.

use std::collections::{HashMap, HashSet};

use crate::locus::location_span;
use crate::perl::{
    info_written, perl_number_is_not_an_integer, perl_numeric_value, perl_pos_is_a_number,
    perl_true,
};
use crate::record::{
    perl_skip_route, record_location, reference_row_beside, reference_rows_starting_at,
    reference_start, span_location, tuples_at, tuples_under, vep_rs_record_rows, LocatedSet,
};
use crate::rules::breakend::{is_breakend_allele, non_integer_end};
use crate::rules::gvcf::perl_split;
use crate::strings::{Id, Sides, Tuple, TupleSet};
use crate::vcf::VcfVariant;

/// A dot-notation breakend without SVTYPE: the reference's two-coordinate tuple and the
/// scored side's point under the same breakend allele on the same feature, set aside as a
/// pair: `(excluded_rust, excluded_perl)`. The point `c:(P+1)` is the span's start or the
/// base after it: the end of `c:P-(P+1)` for a one-base REF, the start of
/// `c:(P+1)-(P+len(REF)-1)` for a longer REF beside an ALT of another length, the second
/// base of `c:P-(P+len(REF)-1)` beside an ALT of the REF's length. The allele is
/// `chromosome_breakpoint` for a bare `N.`, and the `N.` and `N[21:30425700[` strings of
/// one whose INFO names a mate through CHR2 and END2 (Parser/VCF.pm:491-497).
///
/// The points of one feature under one allele go to the first reference tuple that reaches
/// them, the reference tuples taken in the order of the set.
pub fn filter_dot_breakend_pairs(sides: &Sides<'_>) -> (TupleSet, TupleSet) {
    let strings = sides.strings;
    let mut breakend: HashMap<Id, bool> = HashMap::new();
    let mut taken: HashSet<[Id; 4]> = HashSet::new();
    let mut excluded_rust: Vec<Tuple> = Vec::new();
    let mut excluded_perl: Vec<Tuple> = Vec::new();
    for run in sides
        .perl_only
        .as_slice()
        .chunk_by(|a, b| a.location() == b.location())
    {
        let Some((chrom, start, end)) = location_span(strings.text(run[0].location())) else {
            continue;
        };
        if start >= end {
            continue;
        }
        let points =
            [start, start.saturating_add(1)].map(|point| strings.get(&format!("{chrom}:{point}")));
        for &reference in run {
            let is_breakend = *breakend
                .entry(reference.allele())
                .or_insert_with(|| is_breakend_allele(strings.text(reference.allele())));
            if !is_breakend {
                continue;
            }
            for point in points.into_iter().flatten() {
                let key = [
                    point,
                    reference.allele(),
                    reference.feature(),
                    reference.feature_type(),
                ];
                let partners = sides.rust_only.with_prefix(&key);
                if !partners.is_empty() && taken.insert(key) {
                    excluded_perl.push(reference);
                    excluded_rust.extend_from_slice(partners);
                    break;
                }
            }
        }
    }
    (
        TupleSet::from_unsorted(excluded_rust),
        TupleSet::from_unsorted(excluded_perl),
    )
}

/// The reference's tuples at the Location of a record whose ALT column is `.` and whose
/// true INFO/SVTYPE routes it to the structural path.
pub fn filter_missing_alt_annotated(sides: &Sides<'_>, records: &[VcfVariant]) -> TupleSet {
    let mut excluded: Vec<Tuple> = Vec::new();
    for v in records {
        if v.alt != "." || !perl_true(info_written(&v.info, "SVTYPE")) {
            continue;
        }
        let (chrom, start, end) = record_location(v);
        let location = span_location(&chrom, start, end);
        excluded.extend_from_slice(tuples_at(sides.perl, sides.strings, &location));
    }
    TupleSet::from_unsorted(excluded)
}

/// The reference's tuples at the Location of a record whose every ALT is `*` and whose true
/// INFO/SVTYPE routes it to the structural path, where the type names the allele
/// (`SVTYPE=DEL` a deletion over END, `SVTYPE=BND` a `chromosome_breakpoint` at POS + 1;
/// Parser/VCF.pm:237-245, 469-475); `*` names no allele of the record to vep-rs, and
/// without SVTYPE the reference writes nothing.
pub fn filter_star_alone_svtype_annotated(sides: &Sides<'_>, records: &[VcfVariant]) -> TupleSet {
    let mut excluded: Vec<Tuple> = Vec::new();
    for v in records {
        if !v.alt.starts_with('*') {
            continue;
        }
        let alts = perl_split(&v.alt, ',');
        if alts.iter().any(|alt| *alt != "*") || !perl_true(info_written(&v.info, "SVTYPE")) {
            continue;
        }
        let (chrom, start, end) = record_location(v);
        let location = span_location(&chrom, start, end);
        excluded.extend_from_slice(tuples_at(sides.perl, sides.strings, &location));
    }
    TupleSet::from_unsorted(excluded)
}

/// The scored side's one-sided tuples at the Locations and under the alleles vep-rs prints
/// for a record.
fn push_record_rows(v: &VcfVariant, sides: &Sides<'_>, excluded: &mut Vec<Tuple>) {
    for (location, allele) in vep_rs_record_rows(v) {
        excluded.extend_from_slice(tuples_under(
            &sides.rust_only,
            sides.strings,
            &location,
            &allele,
        ));
    }
}

/// The reference's tuples whose Location carries a coordinate that is not an integer
/// (`21:30425651-30425949.5`, `21:30425651-NaN`, `21:30425949.0`), an INFO/END, SVLEN or
/// POS taken as written, and vep-rs's tuples of the same input records, annotated at the
/// numified value: `(excluded_perl, excluded_rust)`. A record is a reference row's when the
/// row's start, read as its number, is the start the reference prints for the record (POS +
/// 1 on the structural path, POS on the sequence path); the record's vep-rs tuples are the
/// one-sided ones at the Locations and under the alleles vep-rs prints for it. With a FASTA
/// the reference expands a `<CNV:TR>` whose SVLEN or END is a number that is not an integer
/// to a literal insertion whose coordinates fail `start = end + 1` (Parser.pm:564-571) and
/// writes nothing for the record; vep-rs's rows at the numified run are set aside when the
/// record is on no skip route under `max_sv_size` and the reference wrote nothing from the
/// base before the span to the base after it.
pub fn filter_non_integer_end_rows(
    sides: &Sides<'_>,
    records: &[VcfVariant],
    max_sv_size: i128,
) -> (TupleSet, TupleSet) {
    let strings = sides.strings;
    let mut excluded_perl: Vec<Tuple> = Vec::new();
    let mut starts: HashSet<(&str, i128)> = HashSet::new();
    for run in sides
        .perl
        .as_slice()
        .chunk_by(|a, b| a.location() == b.location())
    {
        let location = strings.text(run[0].location());
        if !non_integer_end(location) {
            continue;
        }
        excluded_perl.extend_from_slice(run);
        let (chrom, coords) = location.split_once(':').unwrap_or((location, ""));
        let first = coords.split('-').next().unwrap_or("");
        starts.insert((chrom, perl_numeric_value(first).int()));
    }
    let mut excluded_rust: Vec<Tuple> = Vec::new();
    if !starts.is_empty() {
        for v in records {
            let (chrom, start) = reference_start(v);
            if starts.contains(&(chrom.as_str(), start)) {
                push_record_rows(v, sides, &mut excluded_rust);
            }
        }
    }
    let mut perl_spans: Option<Vec<(&str, i128, i128)>> = None;
    for v in records {
        if !v.alt.starts_with('<') {
            continue;
        }
        let upper = v.alt.to_uppercase();
        let alts = perl_split(&upper, ',');
        if alts.is_empty()
            || alts.iter().any(|alt| *alt != "<CNV:TR>")
            || perl_skip_route(v, max_sv_size).is_some()
        {
            continue;
        }
        let span_text = match info_written(&v.info, "SVLEN") {
            Some(svlen) => svlen.split(',').next(),
            None => info_written(&v.info, "END"),
        };
        if !span_text.is_some_and(perl_number_is_not_an_integer) {
            continue;
        }
        let (chrom, start, end) = record_location(v);
        let spans = perl_spans.get_or_insert_with(|| LocatedSet::new(sides.perl, strings).spans());
        if reference_row_beside(spans, &chrom, start.min(end), start.max(end)) {
            continue;
        }
        push_record_rows(v, sides, &mut excluded_rust);
    }
    (
        TupleSet::from_unsorted(excluded_perl),
        TupleSet::from_unsorted(excluded_rust),
    )
}

/// The reference's tuples of a record whose POS is not a non-negative number (`.`, empty,
/// `0x1D04A1D`, fullwidth digits, `30_425_949`, `-5`) at Locations the scored side wrote
/// nothing for: the reference numifies the text (`$start++`, BaseVCF4.pm:263: `.` to 0 then
/// 1, `30_425_949` to 30 then 31, `-5` to -4) and annotates the record from there, where
/// vep-rs refuses the line. The tuples are those at the start the reference prints for the
/// record; on the sequence path only a negative POS passes `looks_like_number`
/// (Parser.pm:475) and is annotated at POS.
pub fn filter_non_numeric_pos_numified(sides: &Sides<'_>, records: &[VcfVariant]) -> TupleSet {
    let mut located: Option<LocatedSet<'_>> = None;
    let mut excluded: Vec<Tuple> = Vec::new();
    for v in records {
        if perl_pos_is_a_number(v.pos_text.as_deref()) {
            continue;
        }
        let perl = located.get_or_insert_with(|| LocatedSet::new(sides.perl, sides.strings));
        excluded.extend(reference_rows_starting_at(v, perl, sides.rust));
    }
    TupleSet::from_unsorted(excluded)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masks::skipped::testing::{adjusted, problems};
    use crate::masks::MAX_SV_SIZE;
    use crate::record::testing::{intergenic, normalised, owned, rec, rec_id, set_of, t, vcf};

    const NO_PROBLEM: [&str; 0] = [];
    const RUN_A: &str = "21:30425588-30425596";
    const RUN_B: &str = "21:30425756-30425763";
    const RUN_C: &str = "21:30430410-30430417";

    fn none() -> TupleSet {
        TupleSet::new()
    }

    /// The tuples of `rows` moved to another Location.
    fn moved(rows: &[[String; 5]], location: &str) -> Vec<[String; 5]> {
        rows.iter()
            .map(|row| t(location, &row[1], &row[2], &row[4]))
            .collect()
    }

    fn rust_a(csq: &str) -> [[String; 5]; 2] {
        [
            t(
                RUN_A,
                "tandem_repeat",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(RUN_A, "tandem_repeat", "ENST00000390690", csq),
        ]
    }

    #[test]
    fn test_a_dot_breakend_pair_is_set_aside_on_both_sides() {
        let perl = [
            t(
                "21:30425650-30425651",
                "chromosome_breakpoint",
                "ENST_1",
                "coding_sequence_variant",
            ),
            t(
                "21:30425650-30425651",
                "chromosome_breakpoint",
                "ENST_2",
                "upstream_gene_variant",
            ),
        ];
        let rust = [
            t(
                "21:30425651",
                "chromosome_breakpoint",
                "ENST_1",
                "coding_sequence_variant",
            ),
            t(
                "21:30425651",
                "N[21:30425700[",
                "ENST_2",
                "upstream_gene_variant",
            ),
        ];
        let sets = owned(&perl, &rust);
        let (excluded_rust, excluded_perl) = filter_dot_breakend_pairs(&sets.sides());
        assert_eq!(excluded_rust, set_of(&sets, &rust[..1]));
        assert_eq!(excluded_perl, set_of(&sets, &perl[..1]));
    }

    #[test]
    fn test_a_missing_alt_the_reference_annotated_is_set_aside_from_the_input() {
        let records = [
            rec(30425650, "N", ".", "SVTYPE=BND"),
            rec(40000000, "N", ".", "SVTYPE=DEL;END=40000100"),
            rec(50000000, "N", ".", "."),
        ];
        let perl = [
            t(
                "21:30425651",
                "chromosome_breakpoint",
                "ENST_1",
                "coding_sequence_variant",
            ),
            t(
                "21:40000001-40000100",
                "deletion",
                "ENST_2",
                "intron_variant",
            ),
            t("21:50000001", "deletion", "ENST_3", "intron_variant"),
        ];
        let sets = owned(&perl, &[]);
        assert_eq!(
            filter_missing_alt_annotated(&sets.sides(), &records),
            set_of(&sets, &perl[..2])
        );
    }

    #[test]
    fn test_a_non_integer_end_row_is_set_aside_and_a_text_end_is_read_as_its_number() {
        let perl = [
            t(
                "21:30425651-30425949.5",
                "deletion",
                "ENST_1",
                "intron_variant",
            ),
            t(
                "21:30425651-3.0425949e7",
                "deletion",
                "ENST_2",
                "intron_variant",
            ),
            t(
                "21:30425651-30425949",
                "deletion",
                "ENST_3",
                "intron_variant",
            ),
        ];
        let sets = owned(&perl, &[]);
        assert_eq!(
            filter_non_integer_end_rows(&sets.sides(), &[], MAX_SV_SIZE),
            (set_of(&sets, &perl[..2]), none())
        );
        for text in [
            "21:30425651-030425949",
            "21:30425651-+30425949",
            "chr21:30425651-30425949",
        ] {
            assert_eq!(normalised(text), "21:30425651-30425949", "{text}");
        }
        assert_eq!(
            normalised("21:30425651-30425949.5"),
            "21:30425651-30425949.5"
        );
    }

    #[test]
    fn test_the_reference_rows_of_a_numified_pos_are_set_aside_from_the_input() {
        let deletion = "SVTYPE=DEL;END=30425949";
        let short = "SVTYPE=DEL;SVLEN=-300";
        let records = vcf(&[
            ["21", ".", "dot", "N", "<DEL>", ".", ".", deletion],
            [
                "21",
                "30_425_949",
                "underscore",
                "N",
                "<DEL>",
                ".",
                ".",
                deletion,
            ],
            ["21", "-5", "neg", "N", "<DEL>", ".", ".", short],
            ["21", "-2147483648", "negbig", "N", "<DEL>", ".", ".", short],
            ["21", "-5", "neg_snv", "A", "G", ".", ".", "."],
            ["21", "0", "zero", "N", "<DEL>", ".", ".", short],
        ]);
        let numified = [
            t(
                "21:1-30425949",
                "deletion",
                "ENST00000262354",
                "transcript_ablation",
            ),
            t(
                "21:1-30425949",
                "deletion",
                "ENST00000282964",
                "transcript_ablation",
            ),
            t(
                "21:31-30425949",
                "deletion",
                "ENST00000262354",
                "transcript_ablation",
            ),
            intergenic("21:-4-295", "deletion"),
            intergenic("21:-2147483647--2147483348", "deletion"),
            intergenic("21:-5", "G"),
        ];
        let shared = [intergenic("21:1-300", "deletion")];
        let sets = owned(&[&numified[..], &shared[..]].concat(), &shared);
        assert_eq!(
            filter_non_numeric_pos_numified(&sets.sides(), &records),
            set_of(&sets, &numified)
        );
        assert_eq!(
            filter_non_integer_end_rows(&sets.sides(), &records, MAX_SV_SIZE),
            (none(), none())
        );
        let (perl_open, rust_open, masks) = adjusted(&sets, &records, -1);
        assert_eq!((perl_open, rust_open), (none(), none()));
        assert_eq!(masks.numified_pos_perl, set_of(&sets, &numified));
        assert_eq!(masks.skipped_problems, NO_PROBLEM);
    }

    #[test]
    fn test_the_reference_rows_of_a_star_only_alt_under_a_true_svtype_are_set_aside() {
        let records = [
            rec_id(30425650, "N", "*", "SVTYPE=DEL;END=30425949", "star"),
            rec_id(30425650, "A", "*", "SVTYPE=DEL;END=30425949", "star_a"),
            rec_id(30425650, "N", "*", "SVTYPE=BND", "star_bnd"),
            rec_id(30425650, "N", "*,*", "SVTYPE=DEL;END=30425949", "two_stars"),
            rec_id(30425650, "N", "*", ".", "bare"),
            rec_id(11000000, "N", "<DEL>", "SVTYPE=DEL;END=11000299", "kept"),
        ];
        let set_aside = [
            t(
                "21:30425651-30425949",
                "deletion",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(
                "21:30425651-30425949",
                "deletion",
                "ENST00000390690",
                "feature_truncation,coding_sequence_variant,5_prime_UTR_variant",
            ),
            t(
                "21:30425651",
                "chromosome_breakpoint",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(
                "21:30425651",
                "chromosome_breakpoint",
                "ENST00000390690",
                "coding_sequence_variant",
            ),
        ];
        let kept = [intergenic("21:11000001-11000299", "deletion")];
        let sets = owned(&[&set_aside[..], &kept[..]].concat(), &kept);
        assert_eq!(
            filter_star_alone_svtype_annotated(&sets.sides(), &records),
            set_of(&sets, &set_aside)
        );
        let deletion = owned(&set_aside[..2], &[]);
        assert_eq!(
            filter_star_alone_svtype_annotated(&deletion.sides(), &records[4..]),
            none()
        );
        let (perl_open, rust_open, masks) = adjusted(&sets, &records, MAX_SV_SIZE);
        assert_eq!((perl_open, rust_open), (none(), none()));
        assert_eq!(masks.star_alone_perl, set_of(&sets, &set_aside));
    }

    #[test]
    fn test_a_dot_breakend_with_an_info_mate_pairs_under_its_own_and_mate_alleles() {
        let mut perl = Vec::new();
        for allele in ["N.", "N[21:30425700["] {
            for (feature, csq) in [
                ("ENST00000334068", "upstream_gene_variant"),
                (
                    "ENST00000390690",
                    "feature_truncation,coding_sequence_variant",
                ),
            ] {
                perl.push(t("21:30425650-30425651", allele, feature, csq));
            }
        }
        let rust = moved(&perl, "21:30425651");
        let sets = owned(&perl, &rust);
        assert_eq!(
            filter_dot_breakend_pairs(&sets.sides()),
            (set_of(&sets, &rust), set_of(&sets, &perl))
        );
        let record = [rec(30425650, "N", "N.", "CHR2=21;END2=30425700")];
        assert_eq!(problems(&sets, &record), NO_PROBLEM);
        let (perl_open, rust_open, masks) = adjusted(&sets, &record, MAX_SV_SIZE);
        assert_eq!((perl_open, rust_open), (none(), none()));
        assert_eq!(
            (masks.dot_pair_perl.len(), masks.dot_pair_rust.len()),
            (4, 4)
        );
        let bare_perl = [
            t(
                "21:30425650-30425651",
                "chromosome_breakpoint",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(
                "21:30425650-30425651",
                "chromosome_breakpoint",
                "ENST00000390690",
                "coding_sequence_variant",
            ),
        ];
        let bare_rust = moved(&bare_perl, "21:30425651");
        let bare = owned(&bare_perl, &bare_rust);
        assert_eq!(
            problems(&bare, &[rec(30425650, "N", "N.", ".")]),
            NO_PROBLEM
        );
        assert_eq!(
            filter_dot_breakend_pairs(&bare.sides()),
            (set_of(&bare, &bare_rust), set_of(&bare, &bare_perl))
        );
        let sequence = owned(
            &[t(
                "21:30425650-30425651",
                "ACG",
                "ENST00000390690",
                "inframe_insertion",
            )],
            &[t(
                "21:30425651",
                "ACG",
                "ENST00000390690",
                "inframe_insertion",
            )],
        );
        assert_eq!(
            filter_dot_breakend_pairs(&sequence.sides()),
            (none(), none()),
            "a sequence allele is not a breakend"
        );
    }

    #[test]
    fn test_a_dot_breakend_pairs_at_the_spans_start_or_the_base_after_it() {
        let rust = [
            t(
                "21:30425651",
                "chromosome_breakpoint",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(
                "21:30425651",
                "chromosome_breakpoint",
                "ENST00000390690",
                "coding_sequence_variant",
            ),
        ];
        for (reference, alt, span) in [
            ("NAC", ".N", "21:30425651-30425652"),
            ("NACG", ".N", "21:30425651-30425653"),
            ("NACGT", "N.", "21:30425651-30425654"),
            ("NAC", "NA.", "21:30425650-30425652"),
            ("NAC", ".AC", "21:30425650-30425652"),
            ("NACG", "NAC.", "21:30425650-30425653"),
            ("N", "N.", "21:30425650-30425651"),
            ("NA", "N.", "21:30425650-30425651"),
        ] {
            let record = [rec(30425650, reference, alt, ".")];
            let perl = moved(&rust, span);
            let sets = owned(&perl, &rust);
            assert_eq!(
                filter_dot_breakend_pairs(&sets.sides()),
                (set_of(&sets, &rust), set_of(&sets, &perl)),
                "{reference} {alt}"
            );
            assert_eq!(problems(&sets, &record), NO_PROBLEM, "{reference} {alt}");
            let (perl_open, rust_open, masks) = adjusted(&sets, &record, MAX_SV_SIZE);
            assert_eq!(
                (perl_open, rust_open),
                (none(), none()),
                "{reference} {alt}"
            );
            assert_eq!(
                (masks.dot_pair_perl.len(), masks.dot_pair_rust.len()),
                (2, 2)
            );
        }
        let wide_end = owned(&moved(&rust, "21:30425649-30425651"), &rust);
        assert_eq!(
            filter_dot_breakend_pairs(&wide_end.sides()),
            (none(), none())
        );
    }

    #[test]
    fn test_a_fractional_svlen_expanded_under_a_fasta_leaves_no_reference_row() {
        for text in ["9.5", "300.7", "-8.5", "NaN", "-Inf", "Infinity", "1e400"] {
            assert!(perl_number_is_not_an_integer(text), "{text}");
        }
        for text in ["9", "9.0", "3e2", "+30425949", "abc", ".", "", "\u{ff14}"] {
            assert!(!perl_number_is_not_an_integer(text), "{text}");
        }
        let elongation = "coding_sequence_variant,feature_elongation";
        let record = [rec(
            30425587,
            "N",
            "<CNV:TR>",
            "SVTYPE=CNV;END=30425596;SVLEN=9.5;RUS=AGC;RUC=4",
        )];
        let rust = rust_a(elongation);
        let sets = owned(&[], &rust);
        assert_eq!(
            filter_non_integer_end_rows(&sets.sides(), &record, MAX_SV_SIZE),
            (none(), set_of(&sets, &rust))
        );
        let (perl_open, rust_open, masks) = adjusted(&sets, &record, MAX_SV_SIZE);
        assert_eq!((perl_open, rust_open), (none(), none()));
        assert_eq!(masks.non_integer_end_rust, set_of(&sets, &rust));
        for (pos, info, run) in [
            (
                30425755,
                "SVTYPE=CNV;END=30425763;SVLEN=8.5;RUS=AG;RUC=5",
                RUN_B,
            ),
            (
                30430409,
                "SVTYPE=CNV;END=30430417;SVLEN=8.5;RUS=CT;RUC=5",
                RUN_C,
            ),
        ] {
            let rows = [t(run, "tandem_repeat", "ENST00000390690", elongation)];
            let sets = owned(&[], &rows);
            assert_eq!(
                filter_non_integer_end_rows(
                    &sets.sides(),
                    &[rec(pos, "N", "<CNV:TR>", info)],
                    MAX_SV_SIZE
                ),
                (none(), set_of(&sets, &rows)),
                "{info}"
            );
        }
        let symbolic = [t(
            "21:30425588-30425596.5",
            "tandem_repeat",
            "ENST00000390690",
            "feature_elongation,coding_sequence_variant",
        )];
        let sets = owned(&symbolic, &rust);
        assert_eq!(
            filter_non_integer_end_rows(&sets.sides(), &record, MAX_SV_SIZE),
            (set_of(&sets, &symbolic), set_of(&sets, &rust))
        );
        let nan = [rec(
            30425587,
            "N",
            "<CNV:TR>",
            "SVTYPE=CNV;SVLEN=NaN;RUS=AGC;RUC=4",
        )];
        let pair = [
            t(
                "21:30425587-30425588",
                "tandem_repeat",
                "ENST00000390690",
                elongation,
            ),
            t(
                "21:30425587-30425588",
                "tandem_repeat",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
        ];
        let sets = owned(&[], &pair);
        assert_eq!(perl_skip_route(&nan[0], MAX_SV_SIZE), None);
        assert_eq!(
            filter_non_integer_end_rows(&sets.sides(), &nan, MAX_SV_SIZE),
            (none(), set_of(&sets, &pair))
        );
        let (perl_open, rust_open, _) = adjusted(&sets, &nan, MAX_SV_SIZE);
        assert_eq!((perl_open, rust_open), (none(), none()));
        let deletion = [rec(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=300.7")];
        let deletion_rows = [t(
            "21:30425651-30425950",
            "deletion",
            "ENST00000390690",
            "feature_truncation,coding_sequence_variant,5_prime_UTR_variant",
        )];
        assert_eq!(
            filter_non_integer_end_rows(
                &owned(&[], &deletion_rows).sides(),
                &deletion,
                MAX_SV_SIZE
            ),
            (none(), none())
        );
    }

    #[test]
    fn a_contested_point_goes_to_the_first_reference_tuple_of_the_set() {
        let point = [t(
            "21:101",
            "chromosome_breakpoint",
            "ENST_1",
            "intron_variant",
        )];
        let narrow = t(
            "21:100-101",
            "chromosome_breakpoint",
            "ENST_1",
            "intron_variant",
        );
        let wide = t(
            "21:100-105",
            "chromosome_breakpoint",
            "ENST_1",
            "intron_variant",
        );
        let later = t(
            "21:101-105",
            "chromosome_breakpoint",
            "ENST_1",
            "intron_variant",
        );
        let sets = owned(&[narrow.clone(), wide.clone()], &point);
        assert_eq!(
            filter_dot_breakend_pairs(&sets.sides()),
            (
                set_of(&sets, &point),
                set_of(&sets, std::slice::from_ref(&narrow))
            )
        );
        let sets = owned(&[wide.clone(), later.clone()], &point);
        assert_eq!(
            filter_dot_breakend_pairs(&sets.sides()),
            (
                set_of(&sets, &point),
                set_of(&sets, std::slice::from_ref(&wide))
            )
        );
        let both = [
            point[0].clone(),
            t(
                "21:102",
                "chromosome_breakpoint",
                "ENST_1",
                "intron_variant",
            ),
        ];
        let sets = owned(&[wide.clone(), later.clone()], &both);
        assert_eq!(
            filter_dot_breakend_pairs(&sets.sides()),
            (set_of(&sets, &both), set_of(&sets, &[wide, later]))
        );
    }
}
