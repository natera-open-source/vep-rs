// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The records the 116.2 reference dropped before its input buffer: their rows on the
//! scored side, the cross-check of the records against the outputs, and the gVCF variant
//! site.

use std::collections::{BTreeSet, HashMap, HashSet};

use crate::locus::{location_span, CNV_TR_RUST_ALLELE};
use crate::masks::one_sided::{filter_dot_breakend_pairs, filter_non_integer_end_rows};
use crate::record::{
    perl_sequence_path, perl_skip_route, printed_location, record_alleles, record_chromosome,
    record_location, tuples_under, LocatedSet,
};
use crate::rules::gvcf::{gvcf_site_rows, perl_split};
use crate::rules::pair::is_structural_allele;
use crate::strings::{Id, Sides, Tuple, TupleSet};
use crate::threads::crew_map;
use crate::vcf::VcfVariant;

/// The engine whose rows the reference's drop rules account for.
const VEP_RS: &str = "vep-rs";

/// One input record at the start vep-rs prints for it.
struct Placed<'a> {
    record: &'a VcfVariant,
    route: Option<&'static str>,
    /// The upper end of the span vep-rs prints.
    end: i128,
}

/// The alleles under which vep-rs prints a record, and the low and high end of the span it
/// prints them at; `None` for a sequence record that names no allele.
fn placed_span(v: &VcfVariant) -> Option<(i128, i128)> {
    if perl_sequence_path(v) {
        let placed = gvcf_site_rows(v.pos, &v.reference, &perl_split(&v.alt, ','));
        let low = placed.iter().map(|(_, start, end)| *start.min(end)).min()?;
        let high = placed.iter().map(|(_, start, end)| *start.max(end)).max()?;
        Some((low, high))
    } else {
        let (_, start, end) = record_location(v);
        Some((start.min(end), start.max(end)))
    }
}

/// The alleles under which vep-rs prints a record: its sequence alleles on the sequence
/// path, every allele either engine writes on the structural path.
fn placed_alleles(v: &VcfVariant) -> BTreeSet<String> {
    if perl_sequence_path(v) {
        gvcf_site_rows(v.pos, &v.reference, &perl_split(&v.alt, ','))
            .into_iter()
            .map(|(allele, _, _)| allele)
            .collect()
    } else {
        record_alleles(v)
    }
}

/// The tuples of a side at a start under a record's alleles; when a kept record contests
/// the alleles, only those whose Location ends where the record's does.
fn rows_of(
    side: &LocatedSet<'_>,
    chrom: &str,
    start: i128,
    alleles: &[Id],
    end: i128,
    contested: bool,
) -> Vec<Tuple> {
    let mut rows = Vec::new();
    for &(location, location_end) in side.locations_at(chrom, start) {
        if contested && location_end != end {
            continue;
        }
        for &allele in alleles {
            rows.extend_from_slice(side.tuples_under(location, allele));
        }
    }
    rows
}

/// The records the 116 reference dropped as `vep_skip`, their vep-rs rows, every
/// disagreement between the records and the outputs, and the skip-route residue of a scored
/// engine other than vep-rs: `(excluded_rust, records, problems, scored_engine_residue)`.
///
/// 116.2 writes at least one row for every record it keeps and no row in any format for a
/// record it marks `vep_skip`; vep-rs annotates all of them. A record is named by the
/// input: its skip route over its ALT and INFO. It is confirmed skipped when the reference
/// wrote no row under the record's own Location start and an allele the record produces; a
/// kept record sharing the start, or the whole Location string, with a skipped one
/// therefore never hides it and is never mistaken for it. When two records at one start
/// produce the same allele, the end of the Location tells them apart. A record on the
/// sequence path is placed at POS under its sequence alleles, whatever INFO it carries.
///
/// Set aside: every vep-rs tuple at a confirmed record's start under one of its alleles
/// (the same end rule applies). The start is the one vep-rs prints, the lower end of the
/// span: an incomplete deletion (`END=POS`, `SVLEN=0`, `SVLEN=abc`) is the pair
/// `POS-(POS+1)`, keyed at POS. The gVCF variant site is its own class: its record is
/// confirmed, so it explains its Location in the residue check, and is neither counted nor
/// set aside here.
///
/// The problems, sorted and each once: a skip-route record whose own rows the reference
/// wrote; a vep-rs structural row at a Location the reference wrote nothing for that no
/// skip-route record explains and no other class pairs with a reference row elsewhere (the
/// dot-breakend pair, the record whose reference Location carries a non-integer
/// coordinate); two records at one start whose alleles and Locations coincide, one skipped
/// and one kept. That residue is a problem for vep-rs alone: under any other
/// `scored_engine` the engine's structural rows at such Locations are returned as the
/// fourth member instead.
pub fn skipped_record_rows(
    sides: &Sides<'_>,
    records: &[VcfVariant],
    max_sv_size: i128,
    scored_engine: &str,
) -> (TupleSet, u64, Vec<String>, TupleSet) {
    let strings = sides.strings;
    let perl_at = LocatedSet::new(sides.perl, strings);
    let rust_at = LocatedSet::new(sides.rust, strings);
    // A record's span and route depend on that record alone.
    let placed = crew_map(records, |v| {
        placed_span(v).map(|span| (span, perl_skip_route(v, max_sv_size)))
    });
    let mut records_at: HashMap<(&str, i128), Vec<Placed<'_>>> = HashMap::new();
    for (v, placed) in records.iter().zip(placed) {
        let Some(((low, high), route)) = placed else {
            continue;
        };
        records_at
            .entry((record_chromosome(v), low))
            .or_default()
            .push(Placed {
                record: v,
                route,
                end: high,
            });
    }

    let mut excluded: Vec<Tuple> = Vec::new();
    let mut problems: Vec<String> = Vec::new();
    let mut explained: HashSet<(&str, i128)> = HashSet::new();
    let mut confirmed: u64 = 0;
    for (&(chrom, start), group) in &records_at {
        if group.iter().all(|placed| placed.route.is_none()) {
            continue;
        }
        let alleles: Vec<BTreeSet<String>> = group
            .iter()
            .map(|placed| placed_alleles(placed.record))
            .collect();
        for (index, placed) in group.iter().enumerate() {
            let Some(route) = placed.route else {
                continue;
            };
            let v = placed.record;
            let end = placed.end;
            // A kept record at this start producing one of this record's alleles: the end
            // decides whose rows are whose, and equal ends cannot be told apart.
            let mut rivals = group.iter().enumerate().filter(|(other, rival)| {
                *other != index
                    && rival.route.is_none()
                    && !alleles[*other].is_disjoint(&alleles[index])
            });
            let contested = rivals.clone().next().is_some();
            if let Some((_, twin)) = rivals.find(|(_, rival)| rival.end == end) {
                let twin = twin.record;
                problems.push(format!(
                    "{} {}>{} and {} {}>{} share {chrom}:{start}-{end} and an allele; \
                     the skipped one cannot be told from the kept one",
                    v.vid, v.reference, v.alt, twin.vid, twin.reference, twin.alt
                ));
                continue;
            }
            let own: Vec<Id> = alleles[index]
                .iter()
                .filter_map(|allele| strings.get(allele))
                .collect();
            let kept = rows_of(&perl_at, chrom, start, &own, end, contested);
            if !kept.is_empty() {
                let written: BTreeSet<&str> =
                    kept.iter().map(|t| strings.text(t.allele())).collect();
                problems.push(format!(
                    "{} {}>{} at {chrom}:{start} is on the {route} route, yet the reference wrote \
                     {} row(s) under its own allele(s) ({}): a reference run \
                     under --dont_skip or another --max_sv_size (pass --reference-max-sv-size)",
                    v.vid,
                    v.reference,
                    v.alt,
                    kept.len(),
                    written.into_iter().collect::<Vec<_>>().join(", ")
                ));
                continue;
            }
            explained.insert((chrom, start));
            if route == "gvcf_variant_site" {
                continue;
            }
            confirmed += 1;
            excluded.extend(rows_of(&rust_at, chrom, start, &own, end, contested));
        }
    }

    // The residue: a Location with scored structural rows and no reference row that no
    // confirmed record explains and no pair class reaches.
    let (dot_pair_rust, _) = filter_dot_breakend_pairs(sides);
    let (_, non_integer_rust) = filter_non_integer_end_rows(sides, records, max_sv_size);
    for paired in dot_pair_rust.iter().chain(non_integer_rust.iter()) {
        if let Some((chrom, start, _)) = location_span(strings.text(paired.location())) {
            explained.insert((chrom, start));
        }
    }
    let mut residue: Vec<Tuple> = Vec::new();
    for run in sides
        .rust
        .as_slice()
        .chunk_by(|a, b| a.location() == b.location())
    {
        let location = run[0].location();
        if !sides.perl.with_prefix(&[location]).is_empty() {
            continue;
        }
        let text = strings.text(location);
        let Some((chrom, start, _)) = location_span(text) else {
            continue;
        };
        if explained.contains(&(chrom, start)) {
            continue;
        }
        for &scored in run {
            let allele = strings.text(scored.allele());
            if !is_structural_allele(allele) || allele == CNV_TR_RUST_ALLELE {
                continue;
            }
            if scored_engine == VEP_RS {
                problems.push(format!(
                    "{text} {allele}: the reference wrote nothing at this Location and no \
                     input record on a skip route starts there"
                ));
            } else {
                residue.push(scored);
            }
        }
    }
    problems.sort_unstable();
    problems.dedup();
    (
        TupleSet::from_unsorted(excluded),
        confirmed,
        problems,
        TupleSet::from_unsorted(residue),
    )
}

/// The mask over [`skipped_record_rows`] without its problems: the vep-rs tuples of every
/// record the reference dropped as `vep_skip` or refused in `validate_vf`, the gVCF variant
/// site excepted, and the number of such records.
pub fn filter_reference_skipped_records(
    sides: &Sides<'_>,
    records: &[VcfVariant],
    max_sv_size: i128,
) -> (TupleSet, u64) {
    let (excluded, confirmed, _, _) = skipped_record_rows(sides, records, max_sv_size, VEP_RS);
    (excluded, confirmed)
}

/// The problems of [`skipped_record_rows`]; empty when the records and the outputs agree.
pub fn check_skip_routes(
    sides: &Sides<'_>,
    records: &[VcfVariant],
    max_sv_size: i128,
    scored_engine: &str,
) -> Vec<String> {
    skipped_record_rows(sides, records, max_sv_size, scored_engine).2
}

/// The vep-rs tuples of every gVCF variant site (`G,<NON_REF>`) the reference dropped as an
/// unsupported type, and the number of such records. vep-rs reads the site as its sequence
/// alleles, placed and minimised as its parser places them (the anchor base trimmed from an
/// indel, a lone allele minimised from both ends), so the tuples set aside are those at each
/// allele's own Location under its own allele string; a site the reference wrote any row
/// for, from POS to the base after REF, is not set aside.
pub fn filter_gvcf_variant_sites(
    sides: &Sides<'_>,
    records: &[VcfVariant],
    max_sv_size: i128,
) -> (TupleSet, u64) {
    // A site lists the `<NON_REF>` block; whether a record is one depends on that record
    // alone.
    let on_the_route = crew_map(records, |v| {
        v.alt.contains('<') && perl_skip_route(v, max_sv_size) == Some("gvcf_variant_site")
    });
    let mut perl_at: Option<LocatedSet<'_>> = None;
    let mut excluded: Vec<Tuple> = Vec::new();
    let mut sites: u64 = 0;
    for (v, on_the_route) in records.iter().zip(on_the_route) {
        if !on_the_route {
            continue;
        }
        let chrom = record_chromosome(v);
        let perl = perl_at.get_or_insert_with(|| LocatedSet::new(sides.perl, sides.strings));
        let bases = v.reference.chars().count() as i128;
        if (0..=bases).any(|offset| perl.starts_at(chrom, v.pos.saturating_add(offset))) {
            continue;
        }
        sites += 1;
        for (allele, start, end) in gvcf_site_rows(v.pos, &v.reference, &perl_split(&v.alt, ',')) {
            let location = printed_location(chrom, start, end);
            excluded.extend_from_slice(tuples_under(sides.rust, sides.strings, &location, &allele));
        }
    }
    (TupleSet::from_unsorted(excluded), sites)
}

#[cfg(test)]
pub(crate) mod testing {
    use std::collections::BTreeMap;

    use super::*;
    use crate::cache::CacheAuthority;
    use crate::masks::{apply_masks, MaskContext, MaskResult, MAX_SV_SIZE};
    use crate::strings::OwnedSides;

    /// The mask over the skipped records under the default cap.
    pub(crate) fn skipped(owned: &OwnedSides, records: &[VcfVariant]) -> (TupleSet, u64) {
        filter_reference_skipped_records(&owned.sides(), records, MAX_SV_SIZE)
    }

    /// The problems of the cross-check under the default cap, vep-rs the scored engine.
    pub(crate) fn problems(owned: &OwnedSides, records: &[VcfVariant]) -> Vec<String> {
        check_skip_routes(&owned.sides(), records, MAX_SV_SIZE, VEP_RS)
    }

    /// The one-sided tuples every 116.2 mask leaves in place, each side, with the mask result.
    pub(crate) fn adjusted(
        owned: &OwnedSides,
        records: &[VcfVariant],
        max_sv_size: i128,
    ) -> (TupleSet, TupleSet, MaskResult) {
        let sides = owned.sides();
        let cache = CacheAuthority::default();
        let context = MaskContext {
            reference_release: "116.2",
            cache: &cache,
            scored_engine: VEP_RS,
            input_variants: Some(records),
            max_sv_size,
            perl_paths: &[],
            rust_paths: &[],
            input_by_output: &BTreeMap::new(),
            coding_ends: None,
        };
        let masks = apply_masks(&sides, &context).expect("the 116.2 masks run");
        (
            sides.perl_only.difference(&masks.all_perl()),
            sides.rust_only.difference(&masks.all_rust()),
            masks,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::testing::{adjusted, problems, skipped};
    use super::*;
    use crate::masks::MAX_SV_SIZE;
    use crate::perl::info_fields;
    use crate::record::testing::{intergenic, owned, rec, rec_id, set_of, t, vcf};
    use crate::record::{perl_so_abbrev, vep_rs_record_rows};
    use crate::vcf::classify_variant;

    const SPAN: &str = "21:30425651-30425949";
    const NO_PROBLEM: [&str; 0] = [];

    fn none() -> TupleSet {
        TupleSet::new()
    }

    fn route(v: &VcfVariant) -> Option<&'static str> {
        perl_skip_route(v, MAX_SV_SIZE)
    }

    #[test]
    fn test_every_vep_rs_row_of_a_skipped_record_is_excluded() {
        let records = [
            rec_id(50000, "N", "<CPX>", "END=60000;SVTYPE=CPX", "cpx"),
            rec_id(30000000, "N", "<DEL>", "SVTYPE=DEL", "del_no_end"),
            rec_id(70000, "N", "<NON_REF>", "END=70100", "gvcf"),
            rec_id(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", "kept"),
        ];
        let rust = [
            t("21:50001-60000", "CPX", "ENST_A", "intron_variant"),
            t("21:50001-60000", "CPX", "ENST_B", "coding_sequence_variant"),
            t("21:30000001", "deletion", "ENST_C", "feature_truncation"),
            intergenic("21:70001-70100", "<NON_REF>"),
            t("21:1001-2000", "deletion", "ENST_D", "feature_truncation"),
        ];
        let sets = owned(&rust[4..], &rust);
        let (excluded, records_n) = skipped(&sets, &records);
        assert_eq!(excluded, set_of(&sets, &rust[..4]));
        assert_eq!(records_n, 3);
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
    }

    #[test]
    fn test_a_kept_record_at_the_same_start_neither_hides_nor_joins_a_skipped_one() {
        let records = [
            rec_id(27330442, "N", "<NON_REF>", "END=27330631", "nonref"),
            rec_id(
                27330442,
                "N",
                "<DEL>",
                "END=27330700;SVTYPE=DEL",
                "multi_span",
            ),
            rec_id(16450356, "N", "<CPX>", "END=16456471;SVTYPE=CPX", "cpx"),
            rec_id(16450356, "N", "<DEL>", "END=16456471;SVTYPE=DEL", "del"),
        ];
        let perl = [
            t(
                "21:27330443-27330700",
                "deletion",
                "ENST_M",
                "intron_variant",
            ),
            t(
                "21:16450357-16456471",
                "deletion",
                "ENST_D",
                "feature_truncation",
            ),
        ];
        let only_rust = [
            intergenic("21:27330443-27330631", "<NON_REF>"),
            t(
                "21:16450357-16456471",
                "CPX",
                "ENST_D",
                "coding_sequence_variant",
            ),
            t("21:16450357-16456471", "CPX", "ENST_E", "intron_variant"),
        ];
        let sets = owned(&perl, &[&perl[..], &only_rust[..]].concat());
        let (excluded, n) = skipped(&sets, &records);
        assert_eq!(excluded, set_of(&sets, &only_rust));
        assert_eq!(n, 2);
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
    }

    #[test]
    fn test_a_skipped_breakend_record_is_matched_by_its_own_end_and_mate_strings() {
        let records = [
            rec_id(
                18457512,
                "N",
                "<BND>",
                "END=18457513;SVTYPE=BND;SVLEN=14027179;CHR2=21;POS2=32484691;END2=32484692",
                "giant",
            ),
            rec_id(
                16316682,
                "N",
                "<BND>",
                "CHR2=chr21;END=16316682;SVLEN=15787275;SVTYPE=BND",
                "giant_no_end2",
            ),
            rec_id(18457512, "N", "N[22:5000[", "SVTYPE=BND", "kept_bnd"),
        ];
        let perl = [t("21:18457513", "N[22:5000[", "ENST_K", "intron_variant")];
        let only_rust = [
            t("21:18457513-32484691", "N.", "ENST_1", "intron_variant"),
            t(
                "21:18457513-32484691",
                "N[21:32484692[",
                "ENST_2",
                "downstream_gene_variant",
            ),
            t("21:16316683-32103957", "N.", "ENST_3", "intron_variant"),
            t(
                "21:16316683-32103957",
                "N[21:16316682[",
                "ENST_4",
                "intron_variant",
            ),
        ];
        let sets = owned(&perl, &[&perl[..], &only_rust[..]].concat());
        let (excluded, n) = skipped(&sets, &records);
        assert_eq!(excluded, set_of(&sets, &only_rust));
        assert_eq!(n, 2);
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
    }

    #[test]
    fn test_a_mixed_alt_list_is_set_aside_whole_including_its_sequence_allele_rows() {
        let records = [rec_id(30165800, "N", "<DEL>,T", "END=30166000", "mixed")];
        let rust = [
            t(
                "21:30165801-30166000",
                "deletion",
                "ENST_1",
                "feature_truncation",
            ),
            t("21:30165801-30166000", "T", "ENST_1", "missense_variant"),
        ];
        let sets = owned(&[], &rust);
        assert_eq!(skipped(&sets, &records), (set_of(&sets, &rust), 1));
    }

    #[test]
    fn test_a_record_the_reference_kept_is_never_touched() {
        let records = [
            rec_id(30165700, "N", "<DEL>,<DUP>", "END=30165900", "cnv"),
            rec_id(
                18457512,
                "N",
                "<BND>",
                "SVTYPE=BND;CHR2=21;END2=32484692",
                "bnd_in_cap",
            ),
        ];
        let perl = [
            t(
                "21:30165701-30165900",
                "copy_number_variation",
                "ENST_1",
                "3_prime_UTR_variant",
            ),
            t("21:18457513", "N.", "ENST_2", "intron_variant"),
        ];
        let rust = [
            t(
                "21:30165701-30165900",
                "deletion",
                "ENST_1",
                "feature_truncation,3_prime_UTR_variant",
            ),
            t(
                "21:30165701-30165900",
                "duplication",
                "ENST_1",
                "feature_elongation,3_prime_UTR_variant",
            ),
            t("21:18457513", "N.", "ENST_2", "intron_variant"),
            t("21:18457513", "N.", "ENST_3", "transcript_ablation"),
        ];
        let sets = owned(&perl, &rust);
        assert_eq!(skipped(&sets, &records), (none(), 0));
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
    }

    #[test]
    fn test_a_tandem_repeat_record_is_outside_the_class() {
        let records = [rec_id(
            899,
            "N",
            "<CNV:TR>",
            "END=950;SVLEN=51;RUS=AC;RUC=10",
            "tr",
        )];
        let rust = [t(
            "21:900-950",
            "tandem_repeat",
            "ENST_T",
            "feature_elongation",
        )];
        let perl = [t("21:950-951", "ACG", "ENST_T", "inframe_insertion")];
        let sets = owned(&perl, &rust);
        assert_eq!(skipped(&sets, &records), (none(), 0));
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
    }

    #[test]
    fn test_an_info_value_is_read_to_a_second_equals_sign() {
        assert_eq!(
            info_fields("SVTYPE==DEL;END=30425949"),
            vec![("SVTYPE", ""), ("END", "30425949")]
        );
        let fields = info_fields("SVTYPE=BND;CHR2=21;END2=30425700=5");
        assert_eq!(crate::perl::info_value(&fields, "END2"), Some("30425700"));
        assert_eq!(
            info_fields("END2=30425700=5=7;CHR2"),
            vec![("END2", "30425700"), ("CHR2", "")]
        );
        let eq_eq = rec_id(30425650, "N", "G", "SVTYPE==DEL;END=30425949", "eq_eq");
        let seq = rec_id(30425650, "A", "G", "SVTYPE=DEL;END=30425949", "seq");
        assert!(perl_sequence_path(&eq_eq));
        assert_eq!(route(&eq_eq), None);
        assert_eq!(
            vep_rs_record_rows(&eq_eq).into_iter().collect::<Vec<_>>(),
            vec![(String::from("21:30425650"), String::from("G"))]
        );
        let rows = [
            t(
                "21:30425650",
                "G",
                "ENST00000390690",
                "coding_sequence_variant",
            ),
            t(
                "21:30425650",
                "G",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t("21:30425650", "G", "ENST00000390690", "missense_variant"),
        ];
        let sets = owned(&rows, &rows);
        let records = [eq_eq, seq];
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
        assert_eq!(skipped(&sets, &records), (none(), 0));
        let bnd = [rec_id(
            30425650,
            "N",
            "<BND>",
            "SVTYPE=BND;CHR2=21;END2=30425700=5",
            "bnd",
        )];
        assert_eq!(
            record_alleles(&bnd[0]).into_iter().collect::<Vec<_>>(),
            [".N", "N.", "N[21:30425700[", "chromosome_breakpoint"]
        );
        let mut mate_rows = Vec::new();
        for allele in ["N.", "N[21:30425700["] {
            for (feature, csq) in [
                ("ENST00000334068", "upstream_gene_variant"),
                (
                    "ENST00000390690",
                    "feature_truncation,coding_sequence_variant",
                ),
            ] {
                mate_rows.push(t("21:30425651", allele, feature, csq));
            }
        }
        let sets = owned(&mate_rows, &mate_rows);
        assert_eq!(problems(&sets, &bnd), NO_PROBLEM);
        let (perl_open, rust_open, _) = adjusted(&sets, &bnd, MAX_SV_SIZE);
        assert_eq!((perl_open, rust_open), (none(), none()));
    }

    #[test]
    fn test_an_infinite_svlen_or_end_is_a_saturated_integer_end() {
        let saturated: i128 = (1 << 63) - 1;
        let location = |info: &str| record_location(&rec(1000, "N", "<DUP>", info));
        let at = |end: i128| (String::from("21"), 1001, end);
        assert_eq!(location("SVTYPE=DUP;SVLEN=1e400"), at(1000 + saturated));
        assert_eq!(location("SVTYPE=DUP;END=1e400"), at(saturated));
        assert_eq!(location("SVTYPE=DUP;END=-1e400"), at(-saturated));
        // The end is an integer by its type; each is exact where a double would round.
        assert_eq!(location("SVTYPE=DUP;SVLEN=1e400").2 % 1000, 807);
        assert_eq!(location("SVTYPE=DUP;END=1e400").2 % 1000, 807);
        let oversize = [rec(1000, "N", "<DUP>", "SVTYPE=DUP;SVLEN=1e400")];
        let rust = [t(
            "21:1001-9223372036854775807",
            "duplication",
            "ENST_1",
            "intron_variant",
        )];
        let sets = owned(&[], &rust);
        assert_eq!(skipped(&sets, &oversize), (set_of(&sets, &rust), 1));
        assert_eq!(problems(&sets, &oversize), NO_PROBLEM);
    }

    #[test]
    fn test_a_trailing_empty_alt_field_is_no_field() {
        let r = |alt: &str, info: &str| rec(1000, "N", alt, info);
        assert_eq!(route(&r("<DEL>,", "SVTYPE=DEL;END=1100")), None);
        assert_eq!(route(&r("N.", "SVTYPE=DEL/;END=1100")), None);
        assert_eq!(
            route(&r("<DEL>,,", "END=1100")),
            None,
            "two trailing empty fields are none"
        );
        assert_eq!(
            route(&r("<DEL>,,T", "END=1100")),
            Some("unsupported_type"),
            "a middle empty field is kept"
        );
        assert_eq!(perl_so_abbrev("<DEL>/"), Some("DEL"));
        assert_eq!(perl_so_abbrev("DEL/"), Some("DEL"));
        assert_eq!(perl_so_abbrev("<DEL>//<DUP>"), Some("CNV"));
        assert_eq!(classify_variant("N", "<DEL>,"), "Symbolic_DEL");
        assert_eq!(classify_variant("N", "<DEL>,<DUP>"), "Multi_allelic");
        assert_eq!(
            record_alleles(&r("<DEL>,", "END=2000;SVTYPE=DEL"))
                .into_iter()
                .collect::<Vec<_>>(),
            ["deletion"]
        );
        assert_eq!(
            record_location(&r("N.,", "END=2000")),
            (String::from("21"), 1001, 1001)
        );
        assert_eq!(
            record_location(&r("<DEL>,", "END=2000")),
            (String::from("21"), 1001, 2000)
        );
        let rust = [t("21:30425651", "-", "ENST_1", "frameshift_variant")];
        let sets = owned(&[], &rust);
        let site = |alt: &str| {
            filter_gvcf_variant_sites(
                &sets.sides(),
                &[rec(30425650, "ATG", alt, ".")],
                MAX_SV_SIZE,
            )
        };
        assert_eq!(site("AG,<NON_REF>,"), (set_of(&sets, &rust), 1));
        assert_eq!(site("AG,"), (none(), 0));
    }

    #[test]
    fn test_an_absent_location_no_skip_route_record_explains_is_a_problem() {
        let kept = [rec_id(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", "kept")];
        let rust = [
            t("21:1001-2000", "deletion", "ENST_D", "feature_truncation"),
            t("21:9999", "deletion", "ENST_X", "intron_variant"),
        ];
        let found = problems(&owned(&[], &rust), &kept);
        assert_eq!(found.len(), 2);
        assert!(found
            .iter()
            .all(|problem| problem.contains("no input record on a skip route starts there")));
        assert_eq!(
            found[0],
            "21:1001-2000 deletion: the reference wrote nothing at this Location and no input \
             record on a skip route starts there"
        );
    }

    #[test]
    fn test_a_skip_route_record_with_reference_rows_under_its_own_allele_is_a_problem() {
        let cpx = [rec_id(
            45716815,
            "N",
            "<CPX>",
            "END=45730272;SVTYPE=CPX",
            "cpx",
        )];
        let perl = [t(
            "21:45716816-45730272",
            "CPX",
            "ENST_C1",
            "intron_variant",
        )];
        let found = problems(&owned(&perl, &[]), &cpx);
        assert_eq!(
            found,
            [
                "cpx N><CPX> at 21:45716816 is on the unsupported_type route, yet the reference \
              wrote 1 row(s) under its own allele(s) (CPX): a reference run under --dont_skip \
              or another --max_sv_size (pass --reference-max-sv-size)"
            ]
        );
    }

    #[test]
    fn test_reference_rows_of_a_kept_record_at_the_same_start_are_not_the_skipped_records() {
        let records = [
            rec_id(16450356, "N", "<CPX>", "END=16456471;SVTYPE=CPX", "cpx"),
            rec_id(16450356, "N", "<DEL>", "END=16456471;SVTYPE=DEL", "del"),
        ];
        let perl = [t(
            "21:16450357-16456471",
            "deletion",
            "ENST_D",
            "feature_truncation",
        )];
        assert_eq!(problems(&owned(&perl, &perl), &records), NO_PROBLEM);
    }

    #[test]
    fn test_twins_sharing_start_end_and_allele_are_a_problem() {
        let mixed = rec_id(1000, "N", "<DEL>,T", "END=2000", "mixed");
        let kept_shorter = rec_id(
            1000,
            "N",
            "<DEL>",
            "END=2000;SVLEN=-500;SVTYPE=DEL",
            "kept_shorter",
        );
        let perl = [t(
            "21:1001-1500",
            "deletion",
            "ENST_D",
            "feature_truncation",
        )];
        let only_rust = [
            t("21:1001-2000", "deletion", "ENST_E", "transcript_ablation"),
            t("21:1001-2000", "T", "ENST_E", "missense_variant"),
        ];
        let sets = owned(&perl, &[&perl[..], &only_rust[..]].concat());
        let records = [mixed.clone(), kept_shorter];
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
        assert_eq!(skipped(&sets, &records), (set_of(&sets, &only_rust), 1));
        let twin = rec_id(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", "twin");
        let twin_perl = [t(
            "21:1001-2000",
            "deletion",
            "ENST_D",
            "feature_truncation",
        )];
        let twin_rust = [
            twin_perl[0].clone(),
            t("21:1001-2000", "T", "ENST_D", "missense_variant"),
        ];
        assert_eq!(
            problems(&owned(&twin_perl, &twin_rust), &[mixed, twin]),
            [
                "mixed N><DEL>,T and twin N><DEL> share 21:1001-2000 and an allele; the skipped \
              one cannot be told from the kept one"
            ]
        );
    }

    #[test]
    fn test_a_kept_sequence_record_carrying_svtype_and_end_is_not_the_twin_of_a_skipped_record() {
        let kept = rec_id(30425650, "A", "G", "SVTYPE=DEL;END=30425949", "seq");
        let skipped_record = rec_id(30425650, "N", "G", "SVTYPE=DEL ;END=30425949", "space");
        let control = rec_id(30425650, "N", "<DEL>", "SVTYPE=DEL;END=30425949", "del");
        let rows = [
            t("21:30425650", "G", "ENST00000390690", "missense_variant"),
            t(
                "21:30425650",
                "G",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(SPAN, "deletion", "ENST00000334068", "upstream_gene_variant"),
            t(
                SPAN,
                "deletion",
                "ENST00000390690",
                "feature_truncation,coding_sequence_variant,5_prime_UTR_variant",
            ),
        ];
        assert!(perl_sequence_path(&kept));
        assert!(!perl_sequence_path(&skipped_record));
        assert_eq!(route(&skipped_record), Some("unsupported_type"));
        let records = [kept, skipped_record, control];
        let sets = owned(&rows, &rows);
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
        assert_eq!(skipped(&sets, &records), (none(), 1));
        let (perl_open, rust_open, masks) = adjusted(&sets, &records, MAX_SV_SIZE);
        assert_eq!((perl_open, rust_open), (none(), none()));
        assert_eq!(masks.skipped_records, 1);
        assert_eq!(masks.skipped_problems, NO_PROBLEM);
        let snv = rec_id(1001, "A", "T", ".", "snv");
        let mixed = rec_id(1000, "N", "<DEL>,T", "END=2000", "mixed");
        let shared = [t("21:1001", "T", "ENST_D", "missense_variant")];
        let mixed_rows = [
            t("21:1001-2000", "deletion", "ENST_E", "transcript_ablation"),
            t("21:1001-2000", "T", "ENST_E", "missense_variant"),
        ];
        let sets = owned(&shared, &[&shared[..], &mixed_rows[..]].concat());
        let records = [snv.clone(), mixed];
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
        assert_eq!(skipped(&sets, &records), (set_of(&sets, &mixed_rows), 1));
        let twin = rec_id(1000, "N", "<DEL>,T", "END=1001", "twin");
        let twin_rust = [
            shared[0].clone(),
            t("21:1001", "deletion", "ENST_D", "feature_truncation"),
        ];
        assert_eq!(
            problems(&owned(&shared, &twin_rust), &[snv, twin]),
            [
                "twin N><DEL>,T and snv A>T share 21:1001-1001 and an allele; the skipped one \
              cannot be told from the kept one"
            ]
        );
    }

    #[test]
    fn test_an_infinite_span_is_oversize_and_its_pair_rows_are_set_aside() {
        let pair = |location: &str, allele: &str| {
            [
                t(
                    location,
                    allele,
                    "ENST00000390690",
                    "coding_sequence_variant,feature_elongation",
                ),
                t(location, allele, "ENST00000334068", "upstream_gene_variant"),
            ]
        };
        let cases = [
            (
                rec_id(30425650, "N", "<INS>", "SVTYPE=INS;SVLEN=Inf", "ins"),
                pair("21:30425650-30425651", "insertion"),
            ),
            (
                rec_id(
                    30425650,
                    "N",
                    "<DUP:TANDEM>",
                    "SVTYPE=DUP;SVLEN=Inf",
                    "tdup",
                ),
                pair("21:30425650-30425651", "tandem_duplication"),
            ),
            (
                rec_id(
                    30425587,
                    "N",
                    "<CNV:TR>",
                    "SVTYPE=CNV;SVLEN=Inf;RUS=AGC;RUC=4",
                    "tr",
                ),
                pair("21:30425587-30425588", "tandem_repeat"),
            ),
        ];
        for (record, rust) in cases {
            assert_eq!(route(&record), Some("oversize"), "{}", record.alt);
            assert_eq!(
                perl_skip_route(&record, -1),
                None,
                "without a cap the reference run dies on the span"
            );
            let sets = owned(&[], &rust);
            let records = [record];
            assert_eq!(skipped(&sets, &records), (set_of(&sets, &rust), 1));
            assert_eq!(problems(&sets, &records), NO_PROBLEM);
        }
        let deletion = |info: &str| route(&rec(30425650, "N", "<DEL>", info));
        assert_eq!(deletion("SVTYPE=DEL;SVLEN=NaN"), None);
        assert_eq!(deletion("SVTYPE=DEL;END=NaN"), None);
        assert_eq!(deletion("SVTYPE=DEL;END=Inf"), Some("oversize"));
        assert_eq!(deletion("SVTYPE=DEL;SVLEN=-Infinity"), Some("oversize"));
        assert_eq!(
            route(&rec(30425650, "N", "<BND>", "SVTYPE=BND;END=Inf")),
            None
        );
    }

    #[test]
    fn test_agreement_yields_no_problem() {
        let records = [
            rec_id(45716815, "N", "<CPX>", "END=45730272;SVTYPE=CPX", "cpx"),
            rec_id(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", "kept"),
        ];
        let perl = [t(
            "21:1001-2000",
            "deletion",
            "ENST_D",
            "feature_truncation",
        )];
        let rust = [
            perl[0].clone(),
            t("21:45716816-45730272", "CPX", "ENST_C", "intron_variant"),
        ];
        assert_eq!(problems(&owned(&perl, &rust), &records), NO_PROBLEM);
    }

    #[test]
    fn test_the_rows_of_an_unsupported_alt_are_set_aside_under_the_class_term() {
        let upstream = |location: &str, allele: &str| {
            t(location, allele, "ENST00000334068", "upstream_gene_variant")
        };
        let coding =
            |location: &str, allele: &str, csq: &str| t(location, allele, "ENST00000390690", csq);
        let truncation = "feature_truncation,coding_sequence_variant,5_prime_UTR_variant";
        let utr = "coding_sequence_variant,5_prime_UTR_variant";
        let cases = [
            (
                rec_id(30425650, "N", "<:DEL>", "END=30425949", "colon"),
                vec![upstream(SPAN, "<:DEL>"), coding(SPAN, "<:DEL>", utr)],
            ),
            (
                rec_id(30425650, "N", "<Del>", "END=30425949", "lower"),
                vec![
                    upstream(SPAN, "deletion"),
                    coding(SPAN, "deletion", truncation),
                ],
            ),
            (
                rec_id(30425650, "N", "<DEL >", "SVTYPE=DEL;END=30425949", "blank"),
                vec![
                    upstream(SPAN, "deletion"),
                    coding(SPAN, "deletion", truncation),
                ],
            ),
            (
                rec_id(30425650, "N", "<CN=X>", "END=30425949", "cn"),
                vec![
                    upstream(SPAN, "copy_number_variation"),
                    coding(SPAN, "copy_number_variation", utr),
                ],
            ),
            (
                rec_id(30425650, "N", "<Bnd>", ".", "bnd"),
                vec![
                    upstream("21:30425651", "chromosome_breakpoint"),
                    coding(
                        "21:30425651",
                        "chromosome_breakpoint",
                        "coding_sequence_variant",
                    ),
                ],
            ),
            (
                rec_id(30425650, "N", "<<DEL>>,<DEL>", "END=30425949", "doubled"),
                vec![
                    upstream(SPAN, "<<DEL>>/<DEL>"),
                    upstream(SPAN, "deletion"),
                    coding(SPAN, "<<DEL>>/<DEL>", utr),
                    coding(SPAN, "deletion", truncation),
                ],
            ),
        ];
        for (record, rust) in cases {
            assert_eq!(route(&record), Some("unsupported_type"), "{}", record.alt);
            let sets = owned(&[], &rust);
            let (perl_open, rust_open, masks) = adjusted(&sets, &[record], MAX_SV_SIZE);
            assert_eq!((perl_open, rust_open), (none(), none()));
            assert_eq!(masks.skipped_rust, set_of(&sets, &rust));
            assert_eq!(masks.skipped_records, 1);
            assert_eq!(masks.skipped_problems, NO_PROBLEM);
        }
    }

    #[test]
    fn test_an_incomplete_deletions_pair_rows_are_set_aside() {
        let pair = [
            t(
                "21:30425650-30425651",
                "deletion",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(
                "21:30425650-30425651",
                "deletion",
                "ENST00000390690",
                "coding_sequence_variant",
            ),
        ];
        let sets = owned(&[], &pair);
        for info in [
            "SVTYPE=DEL;END=30425650",
            "SVTYPE=DEL;SVLEN=0",
            "SVTYPE=DEL;SVLEN=abc",
            "SVTYPE=DEL;SVLEN=.",
            "SVTYPE=DEL;SVLEN=",
            "SVTYPE=DEL;SVLEN=e5",
            "SVTYPE=DEL;SVLEN=0x1D04A1D",
            "SVTYPE=DEL;SVLEN=.;END=30425949",
            "SVTYPE=DEL;SVLEN=;END=30425949",
            "SVTYPE=DEL;SVLEN=abc;END=30425949",
        ] {
            let record = [rec(30425650, "N", "<DEL>", info)];
            assert_eq!(route(&record[0]), Some("incomplete_deletion"), "{info}");
            assert_eq!(skipped(&sets, &record), (set_of(&sets, &pair), 1), "{info}");
            assert_eq!(problems(&sets, &record), NO_PROBLEM, "{info}");
        }
        let plus = vcf(&[[
            "21",
            "+30425949",
            "plus",
            "N",
            "<DEL>",
            ".",
            ".",
            "SVTYPE=DEL;END=30425949",
        ]]);
        let plus_pair = [
            t(
                "21:30425949-30425950",
                "deletion",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            t(
                "21:30425949-30425950",
                "deletion",
                "ENST00000390690",
                "upstream_gene_variant",
            ),
        ];
        let sets = owned(&[], &plus_pair);
        assert_eq!(skipped(&sets, &plus), (set_of(&sets, &plus_pair), 1));
        assert_eq!(problems(&sets, &plus), NO_PROBLEM);
    }

    #[test]
    fn test_a_gvcf_variant_site_is_its_own_class_apart_from_the_skipped_records() {
        let records = [
            rec(30425650, "A", "G,<NON_REF>", "."),
            rec(50000, "N", "<CPX>", "END=60000;SVTYPE=CPX"),
        ];
        let rust = [
            t("21:30425650", "G", "ENST_1", "missense_variant"),
            t("21:30425650", "G", "ENST_2", "upstream_gene_variant"),
            t("21:50001-60000", "CPX", "ENST_3", "intron_variant"),
        ];
        let sets = owned(&[], &rust);
        let (excluded, sites) = filter_gvcf_variant_sites(&sets.sides(), &records, MAX_SV_SIZE);
        assert_eq!((sites, excluded), (1, set_of(&sets, &rust[..2])));
        let (skipped_rows, records_n) = skipped(&sets, &records);
        assert_eq!((records_n, skipped_rows), (1, set_of(&sets, &rust[2..])));
        assert_eq!(problems(&sets, &records), NO_PROBLEM);
        let (_, _, masks) = adjusted(&sets, &records, MAX_SV_SIZE);
        assert_eq!(
            (masks.skipped_records, masks.gvcf_records),
            (1, 1),
            "the mask run names the site once"
        );
    }

    #[test]
    fn test_a_gvcf_site_allele_is_upper_cased_as_the_parser_writes_it() {
        let rust = [
            t("21:30425650", "G", "ENST00000390690", "missense_variant"),
            t(
                "21:30425650",
                "G",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
        ];
        let sets = owned(&[], &rust);
        for (reference, alt) in [("A", "g,<NON_REF>"), ("a", "g,<non_ref>")] {
            let record = [rec(30425650, reference, alt, ".")];
            assert_eq!(route(&record[0]), Some("gvcf_variant_site"), "{alt}");
            assert_eq!(
                filter_gvcf_variant_sites(&sets.sides(), &record, MAX_SV_SIZE),
                (set_of(&sets, &rust), 1),
                "{alt}"
            );
            let (perl_open, rust_open, _) = adjusted(&sets, &record, MAX_SV_SIZE);
            assert_eq!((perl_open, rust_open), (none(), none()), "{alt}");
        }
    }

    #[test]
    fn test_a_gvcf_site_is_set_aside_at_the_alleles_the_parser_writes() {
        let records = [
            rec(30425650, "ATG", "AG,<NON_REF>", "."),
            rec(40000000, "A", "AT,ATT,<NON_REF>", "."),
            rec(50000000, "N", "<NON_REF>,<NON_REF>", "END=50000100"),
        ];
        let rust = [
            t("21:30425651", "-", "ENST_1", "frameshift_variant"),
            t("21:40000000-40000001", "T", "ENST_2", "frameshift_variant"),
            t("21:40000000-40000001", "TT", "ENST_2", "inframe_insertion"),
            t("21:30425651", "G", "ENST_1", "missense_variant"),
            intergenic("21:50000001-50000100", "<NON_REF>"),
        ];
        let sets = owned(&[], &rust);
        let (excluded, sites) = filter_gvcf_variant_sites(&sets.sides(), &records, MAX_SV_SIZE);
        assert_eq!(sites, 2);
        assert_eq!(excluded, set_of(&sets, &rust[..3]));
        assert_eq!(route(&records[2]), Some("unsupported_type"));
        let r = |alt: &str, info: &str| route(&rec(1000, "N", alt, info));
        assert_eq!(r("*,<NON_REF>", "."), Some("unsupported_type"));
        assert_eq!(r("<DEL>,<NON_REF>", "END=1100"), Some("unsupported_type"));
        // Under a true SVTYPE the block is a deletion allele in both engines: the record is
        // kept.
        assert_eq!(r("G,<NON_REF>", "SVTYPE=DEL;END=1100"), None);
        assert_eq!(
            r("G,<NON_REF>", "SVTYPE=FOO;END=1100"),
            Some("unsupported_type")
        );
        assert_eq!(
            r("G,<NON_REF>", "SVTYPE=;END=1100"),
            Some("gvcf_variant_site")
        );
        let (skipped_rows, n) = skipped(&sets, &records);
        assert_eq!((n, skipped_rows), (1, set_of(&sets, &rust[4..])));
    }

    #[test]
    fn a_scored_engine_other_than_vep_rs_owns_its_unexplained_rows() {
        let kept = [rec_id(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", "kept")];
        let rust = [
            t("21:1001-2000", "deletion", "ENST_D", "feature_truncation"),
            t("21:9999", "deletion", "ENST_X", "intron_variant"),
            t("21:9999", "tandem_repeat", "ENST_X", "intron_variant"),
            t("21:9999", "G", "ENST_X", "intron_variant"),
        ];
        let sets = owned(&[], &rust);
        let (excluded, confirmed, found, residue) =
            skipped_record_rows(&sets.sides(), &kept, MAX_SV_SIZE, "fastvep");
        assert_eq!((excluded, confirmed), (none(), 0));
        assert_eq!(found, NO_PROBLEM);
        assert_eq!(residue, set_of(&sets, &rust[..2]));
        assert_eq!(
            check_skip_routes(&sets.sides(), &kept, MAX_SV_SIZE, "fastvep"),
            NO_PROBLEM
        );
        let (_, _, found, residue) = skipped_record_rows(&sets.sides(), &kept, MAX_SV_SIZE, VEP_RS);
        assert_eq!(found.len(), 2);
        assert_eq!(residue, none());
    }
}
