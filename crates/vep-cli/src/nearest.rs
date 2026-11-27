// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--nearest`: the protein-coding transcription start sites nearest a variant
//! and the identifiers reported for them, as Ensembl VEP's `TranscriptTree`
//! finds them (`TranscriptTree.pm:222-254`, `AnnotationType/Transcript.pm:700-715`
//! at release/116.3).
//!
//! Ensembl VEP keeps every protein-coding transcript's start site in an interval
//! tree, one point per transcript, and asks the tree for one site on each side of
//! the variant: `fetch_nearest_down(end)` and `fetch_nearest_up(start)`. The
//! distance of a site is the smallest difference between it and either end of the
//! variant (`_get_dist`, `:303-317`); every candidate at the minimum is reported,
//! the one found below the end first. Here the sites are the chromosome's sorted
//! [`TssPoint`] table and the two searches are binary searches.
//!
//! A third candidate joins Ensembl VEP's two: the lowest site at or above the
//! start, as the search below the end takes the highest at or below it. Ensembl
//! VEP's search above the start passes over a site on the start and a site one
//! base above it. Its tree holds a site at `tss` as the half-open interval
//! `[tss - 1, tss)` (`:165`), the search is given the 1-based start (`:235`),
//! and the library returns only an interval that begins above its query, so it
//! returns the lowest site at or above `start + 2`. A one-base variant one base
//! below a start site is then given a farther site, or none where neither
//! search returns one (`:239`), and a longer variant whose first base is a
//! start site another site it covers. With the third candidate
//! the passed-over site is named where it is as near, and no site Ensembl VEP
//! names at the smallest distance is lost.
//!
//! Where several transcripts start at the winning coordinate, Ensembl VEP reports
//! the one its tree happens to hold first, which depends on the order its cache
//! directory was listed in; this module reports every transcript at the
//! coordinate, in stable-id order. `gene` and `symbol` de-duplicate their values
//! as `get_nearest` does, so the two agree whenever the transcripts share a gene.

use std::str::FromStr;

use vep_core::minimal::rowless_first_feature;
use vep_core::variant::InputVariant;

use crate::runner::structural_variant_bounds;
use crate::transcript_index::{ChromosomeIndexes, TssPoint};

/// The identifier `--nearest` reports for a transcription start site
/// (`Config.pm:602`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum NearestType {
    /// The transcript stable id.
    Transcript,
    /// The gene stable id.
    Gene,
    /// The gene symbol, else the gene stable id.
    Symbol,
}

impl NearestType {
    fn value(self, point: &TssPoint) -> &str {
        match self {
            Self::Transcript => &point.transcript,
            Self::Gene => &point.gene,
            Self::Symbol => &point.symbol,
        }
    }
}

impl FromStr for NearestType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "transcript" => Ok(Self::Transcript),
            "gene" => Ok(Self::Gene),
            "symbol" => Ok(Self::Symbol),
            other => Err(format!(
                "\"{other}\" is not a valid value for --nearest\nValid values: transcript, gene, symbol"
            )),
        }
    }
}

/// The sites nearest the span `[start, end]` (either order) among `points`,
/// which are sorted by coordinate: every transcript at each winning coordinate,
/// in stable-id order, the candidate found below the end before the ones found
/// above the start (`TranscriptTree.pm:233-236` lists its two in that order,
/// and all can lie inside a span).
///
/// The candidate below is the highest site at or below the span's end, the
/// site `fetch_nearest_down(end)` returns. Ensembl VEP's
/// `fetch_nearest_up(start)` returns the lowest site at or above `start + 2`,
/// since a site at `start + 1` is the interval `[start, start + 1)` of its
/// tree, which holds the query, and a site at `start` is the interval that
/// ends at it. That site is a candidate here too, after the lowest site at or
/// above `start`, which that search passes over. Each candidate's distance is
/// the smaller of its differences to the two ends; all at the smallest are
/// kept.
pub(crate) fn nearest_points(points: &[TssPoint], start: u64, end: u64) -> Vec<&TssPoint> {
    let (low, high) = (start.min(end), start.max(end));
    let down = points[..points.partition_point(|p| p.tss <= high)].last();
    let lowest_from = |from: u64| points.get(points.partition_point(|p| p.tss < from));
    let up = lowest_from(low);
    let above = low.checked_add(2).and_then(lowest_from);
    let mut candidates: Vec<u64> = Vec::with_capacity(3);
    for site in down.into_iter().chain(up).chain(above) {
        if !candidates.contains(&site.tss) {
            candidates.push(site.tss);
        }
    }
    let distance = |tss: u64| low.abs_diff(tss).min(high.abs_diff(tss));
    let Some(best) = candidates.iter().map(|&tss| distance(tss)).min() else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for coordinate in candidates.into_iter().filter(|&tss| distance(tss) == best) {
        let from = points.partition_point(|p| p.tss < coordinate);
        let to = points.partition_point(|p| p.tss <= coordinate);
        let mut at_coordinate: Vec<&TssPoint> = points[from..to].iter().collect();
        at_coordinate.sort_by(|a, b| a.transcript.cmp(&b.transcript));
        out.extend(at_coordinate);
    }
    out
}

/// The `NEAREST` values of a variant: the identifiers of [`nearest_points`] of
/// its span, each once in order of first appearance (`get_nearest`,
/// `AnnotationType/Transcript.pm:707-712`). A structural variant's span runs to
/// its `END`. Empty on a chromosome with no protein-coding transcription start
/// site, which Ensembl VEP writes as an empty list.
pub(crate) fn nearest_values(
    points: &[TssPoint],
    variant: &InputVariant,
    kind: NearestType,
) -> Vec<String> {
    let (start, end) = structural_variant_bounds(variant);
    let mut values: Vec<String> = Vec::new();
    for point in nearest_points(points, start, end) {
        let value = kind.value(point);
        if !value.is_empty() && !values.iter().any(|v| v == value) {
            values.push(value.to_string());
        }
    }
    values
}

/// Sets [`InputVariant::nearest`] from the chromosome's transcription start
/// sites when `--nearest` is on; a chromosome the cache lacks yields the empty
/// list, as Ensembl VEP's empty tree does. An allele of a `--minimal` split
/// record is measured on its own trimmed span, the span of the feature Ensembl
/// VEP measures (AnnotationType/Transcript.pm:170-174), unless the record's
/// first alternate allele writes no row, when the first feature's span is
/// measured ([`rowless_first_feature`]); `runner::apply_record_rules` then gives
/// every allele of the record the first feature's list, the one Ensembl VEP's
/// rows print (OutputFactory.pm:2295-2305, :1017).
pub(crate) fn annotate(
    variant: &mut InputVariant,
    transcripts: &impl ChromosomeIndexes,
    kind: Option<NearestType>,
) {
    let Some(kind) = kind else {
        return;
    };
    let points = transcripts
        .get(&variant.chr)
        .map(|index| index.tss_points())
        .unwrap_or(&[]);
    let first_feature = rowless_first_feature(variant);
    let subject = first_feature.as_ref().unwrap_or(variant);
    variant.nearest = Some(nearest_values(points, subject, kind));
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;

    fn point(tss: u64, transcript: &str, gene: &str, symbol: &str) -> TssPoint {
        TssPoint {
            tss,
            transcript: Arc::from(transcript),
            gene: Arc::from(gene),
            symbol: Arc::from(symbol),
        }
    }

    /// A sorted table: two transcripts of one forward-strand gene, both
    /// starting at 1000, a reverse-strand gene whose start site is its end at
    /// 5000, and a lone site at 9000.
    fn table() -> Vec<TssPoint> {
        vec![
            point(1000, "ENST_B", "ENSG_1", "GENE1"),
            point(1000, "ENST_A", "ENSG_1", "GENE1"),
            point(5000, "ENST_R", "ENSG_2", "GENE2"),
            point(9000, "ENST_F", "ENSG_3", "ENSG_3"),
        ]
    }

    fn ids<'a>(points: &[&'a TssPoint]) -> Vec<&'a str> {
        points.iter().map(|p| &*p.transcript).collect()
    }

    fn snv(pos: u64) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec())
    }

    /// Ensembl VEP's `_get_dist` minimum over the site below the end and the
    /// site above the start (`TranscriptTree.pm:233-253`), reported as the
    /// transcript, the gene and the symbol (`Transcript.pm:747-751`: the symbol
    /// falls back to the gene id).
    #[test]
    fn nearest_symbol_gene_transcript_on_a_hand_built_cache() {
        let table = table();
        let v = snv(6000);
        assert_eq!(
            nearest_values(&table, &v, NearestType::Transcript),
            ["ENST_R"]
        );
        assert_eq!(nearest_values(&table, &v, NearestType::Gene), ["ENSG_2"]);
        assert_eq!(nearest_values(&table, &v, NearestType::Symbol), ["GENE2"]);
        // 8000 is 1000 from 9000 and 3000 from 5000: the site above wins.
        let v = snv(8000);
        assert_eq!(
            nearest_values(&table, &v, NearestType::Symbol),
            ["ENSG_3"],
            "a gene without a symbol reports its id"
        );
        assert_eq!(
            nearest_values(&table, &v, NearestType::Transcript),
            ["ENST_F"]
        );
    }

    /// Both candidates at the minimum are kept, the one below the end first
    /// (`TranscriptTree.pm:233-236, 250-253`); inside a span the two can come
    /// out in descending coordinate order.
    #[test]
    fn ties_at_equal_distance_are_all_returned() {
        let table = table();
        // 7000 is 2000 from 5000 and 2000 from 9000.
        assert_eq!(
            ids(&nearest_points(&table, 7000, 7000)),
            ["ENST_R", "ENST_F"]
        );
        assert_eq!(
            nearest_values(&table, &snv(7000), NearestType::Symbol),
            ["GENE2", "ENSG_3"]
        );
        // Two sites inside one span, each 5 from its nearer end: the highest site
        // at or below the end (5005) before the lowest above the start (4995).
        let inside = vec![
            point(4995, "ENST_LOW", "G", "G"),
            point(5005, "ENST_HIGH", "G", "G"),
        ];
        assert_eq!(
            ids(&nearest_points(&inside, 4990, 5010)),
            ["ENST_HIGH", "ENST_LOW"]
        );
    }

    /// A site the span covers is found by `down` when it is the highest at or
    /// below the end and by `up` when it is the lowest at or above the start;
    /// its distance is the smaller gap to the two ends, zero on an end.
    #[test]
    fn a_site_the_span_covers_is_reported() {
        let table = table();
        let deletion = InputVariant::new("21".into(), 4990, 5010, b"A".repeat(21), b"A".to_vec());
        assert_eq!(
            nearest_values(&table, &deletion, NearestType::Transcript),
            ["ENST_R"]
        );
        // The site at the span's own position.
        assert_eq!(ids(&nearest_points(&table, 5000, 5000)), ["ENST_R"]);
        // An insertion is given as end < start.
        assert_eq!(ids(&nearest_points(&table, 5001, 5000)), ["ENST_R"]);
        // A site on the end of a span is at distance zero and beats a site one
        // base below the start.
        let ends = vec![
            point(4999, "ENST_BELOW", "G", "G"),
            point(5010, "ENST_END", "G", "G"),
        ];
        assert_eq!(ids(&nearest_points(&ends, 5000, 5010)), ["ENST_END"]);
    }

    /// A site one base above the start is the `up` candidate, at distance one
    /// from a one-base variant. Ensembl VEP's tree holds a site at `t` as
    /// `[t - 1, t)` and is asked `fetch_nearest_up(s)` with the 1-based start
    /// (`TranscriptTree.pm:235`), which passes over the interval holding `s`,
    /// so it reports the site below instead, 3999 bases away in this table.
    #[test]
    fn a_site_one_base_above_the_start_is_the_up_candidate() {
        let table = table();
        // 4999: the site at 5000 is one base away and the site at 1000 3999.
        assert_eq!(ids(&nearest_points(&table, 4999, 4999)), ["ENST_R"]);
        // 4998: the site at 5000 wins at distance 2.
        assert_eq!(ids(&nearest_points(&table, 4998, 4998)), ["ENST_R"]);
        // A two-base span ending at 5000 reaches the site from below.
        assert_eq!(ids(&nearest_points(&table, 4999, 5000)), ["ENST_R"]);
        // Sites one base below and one base above a variant tie.
        let beside = vec![
            point(4998, "ENST_BELOW", "G", "G"),
            point(5000, "ENST_ABOVE", "G", "G"),
            point(9000, "ENST_FAR", "G", "G"),
        ];
        assert_eq!(
            ids(&nearest_points(&beside, 4999, 4999)),
            ["ENST_BELOW", "ENST_ABOVE"]
        );
        // A site on the variant is at distance zero and beats the site above it.
        let on = vec![
            point(4999, "ENST_ON", "G", "G"),
            point(5000, "ENST_ABOVE", "G", "G"),
        ];
        assert_eq!(ids(&nearest_points(&on, 4999, 4999)), ["ENST_ON"]);
        // Inside a span the site at start + 1 is one base from the start and
        // beats a higher site six bases from it.
        let inside = vec![
            point(4991, "ENST_BY_START", "G", "G"),
            point(4996, "ENST_HIGHER", "G", "G"),
        ];
        assert_eq!(ids(&nearest_points(&inside, 4990, 5010)), ["ENST_BY_START"]);
    }

    /// A site on the first base of a span of two or more positions is the `up`
    /// candidate, at distance zero. Ensembl VEP's `fetch_nearest_up(s)` returns
    /// the lowest site at or above `s + 2` (`TranscriptTree.pm:235`), and its
    /// search below the end reaches the site on the first base only where the
    /// span covers no higher site, so it names a covered site farther from
    /// both ends, or one of two sites that both lie at distance zero.
    #[test]
    fn a_site_on_the_first_base_is_the_up_candidate() {
        // Sites on both positions of a two-position span, which an insertion
        // between them has too: both at distance zero, the one found below
        // the end first.
        let adjacent = vec![
            point(5000, "ENST_LOW", "G", "G"),
            point(5001, "ENST_HIGH", "G", "G"),
            point(9000, "ENST_FAR", "G", "G"),
        ];
        assert_eq!(
            ids(&nearest_points(&adjacent, 5000, 5001)),
            ["ENST_HIGH", "ENST_LOW"]
        );
        assert_eq!(
            ids(&nearest_points(&adjacent, 5001, 5000)),
            ["ENST_HIGH", "ENST_LOW"]
        );
        // A site on the first base is at distance zero and beats a covered
        // site five bases below the end.
        let inside = vec![
            point(5000, "ENST_FIRST", "G", "G"),
            point(5030, "ENST_INSIDE", "G", "G"),
            point(9000, "ENST_FAR", "G", "G"),
        ];
        assert_eq!(ids(&nearest_points(&inside, 5000, 5035)), ["ENST_FIRST"]);
        // With a covered site on the last base too, the two ends tie.
        assert_eq!(
            ids(&nearest_points(&inside, 5000, 5030)),
            ["ENST_INSIDE", "ENST_FIRST"]
        );
        // The only site a span covers, on its first base, is found below the
        // end as well.
        let lone = vec![
            point(5000, "ENST_FIRST", "G", "G"),
            point(9000, "ENST_FAR", "G", "G"),
        ];
        assert_eq!(ids(&nearest_points(&lone, 5000, 5003)), ["ENST_FIRST"]);
        // A span that starts one base above the lower site covers the higher
        // site alone, and every search returns the covered site, as Ensembl
        // VEP's do: a site below the span is no candidate beside it.
        assert_eq!(ids(&nearest_points(&inside, 5001, 5033)), ["ENST_INSIDE"]);
    }

    /// Ensembl VEP's own candidate above the start stays a candidate. A span
    /// whose one covered site is on its second base, with a site one base
    /// above its end: `fetch_nearest_up(s)` passes over the covered site and
    /// returns the outer one (`TranscriptTree.pm:235`), both lie one base from
    /// an end, and both are named, the covered site first.
    #[test]
    fn ensembl_veps_candidate_above_the_start_is_kept() {
        let outer = vec![
            point(5001, "ENST_SECOND", "G", "G"),
            point(5011, "ENST_OUTER", "G", "G"),
        ];
        assert_eq!(
            ids(&nearest_points(&outer, 5000, 5010)),
            ["ENST_SECOND", "ENST_OUTER"]
        );
        // One base longer, the outer site is covered, on the last base.
        assert_eq!(ids(&nearest_points(&outer, 5000, 5011)), ["ENST_OUTER"]);
        // One base shorter, the outer site lies two bases above the end.
        assert_eq!(ids(&nearest_points(&outer, 5000, 5009)), ["ENST_SECOND"]);
        // A site on the first base is nearer than both.
        let first = vec![
            point(5000, "ENST_FIRST", "G", "G"),
            point(5001, "ENST_SECOND", "G", "G"),
            point(5011, "ENST_OUTER", "G", "G"),
        ];
        assert_eq!(ids(&nearest_points(&first, 5000, 5010)), ["ENST_FIRST"]);
        // The last coordinate has no base two above it.
        let last = vec![point(u64::MAX, "ENST_LAST", "G", "G")];
        assert_eq!(
            ids(&nearest_points(&last, u64::MAX - 1, u64::MAX)),
            ["ENST_LAST"]
        );
    }

    /// The table is built from protein-coding transcripts alone, so a chromosome
    /// whose transcripts are all non-coding has no site and reports the empty
    /// list; a lone site is reported from any distance.
    #[test]
    fn only_protein_coding_tss_enter_the_array() {
        use vep_core::coordinate::Strand;

        use crate::runner::make_runner_test_transcript;
        use crate::transcript_index::TranscriptIndex;

        let transcript = |id: &str, biotype: &str, strand: Strand| {
            let mut t = make_runner_test_transcript(id, "ENSG_X", "21", 100);
            t.end = 900;
            t.biotype = Arc::from(biotype);
            t.strand = strand;
            t
        };
        let index = TranscriptIndex::new(
            vec![
                transcript("ENST_LNC", "lncRNA", Strand::Forward),
                transcript("ENST_PC", "protein_coding", Strand::Reverse),
                transcript("ENST_PSEUDO", "processed_pseudogene", Strand::Forward),
            ],
            5000,
            5000,
        );
        let points = index.tss_points();
        assert_eq!(ids(&points.iter().collect::<Vec<_>>()), ["ENST_PC"]);
        assert_eq!(points[0].tss, 900, "a reverse-strand site is the end");
        assert_eq!(
            nearest_values(points, &snv(50_000_000), NearestType::Transcript),
            ["ENST_PC"]
        );
        let none: Vec<TssPoint> = Vec::new();
        assert!(nearest_values(&none, &snv(1), NearestType::Symbol).is_empty());
    }

    /// Every transcript at the winning coordinate is reported, in stable-id
    /// order, and the gene and symbol values collapse to one per gene.
    #[test]
    fn transcripts_sharing_a_tss_are_all_reported_in_stable_id_order() {
        let table = table();
        assert_eq!(
            ids(&nearest_points(&table, 1500, 1500)),
            ["ENST_A", "ENST_B"]
        );
        let v = snv(1500);
        assert_eq!(
            nearest_values(&table, &v, NearestType::Transcript),
            ["ENST_A", "ENST_B"]
        );
        assert_eq!(nearest_values(&table, &v, NearestType::Gene), ["ENSG_1"]);
        assert_eq!(nearest_values(&table, &v, NearestType::Symbol), ["GENE1"]);
    }

    /// Transcripts of different genes starting at one site: every gene and
    /// every symbol at the site is reported, in the order of the transcripts'
    /// stable ids, where Ensembl VEP reports the gene of the one transcript its
    /// tree holds first.
    #[test]
    fn genes_sharing_a_site_are_all_reported() {
        let shared = vec![
            point(2000, "ENST_Y", "ENSG_Y", "GENEY"),
            point(2000, "ENST_X1", "ENSG_X", "GENEX"),
            point(2000, "ENST_X2", "ENSG_X", "GENEX"),
        ];
        let v = snv(2500);
        assert_eq!(
            nearest_values(&shared, &v, NearestType::Transcript),
            ["ENST_X1", "ENST_X2", "ENST_Y"]
        );
        assert_eq!(
            nearest_values(&shared, &v, NearestType::Gene),
            ["ENSG_X", "ENSG_Y"]
        );
        assert_eq!(
            nearest_values(&shared, &v, NearestType::Symbol),
            ["GENEX", "GENEY"]
        );
    }

    /// A structural variant is measured over its whole span.
    #[test]
    fn a_structural_variant_spans_to_its_end() {
        let table = table();
        let mut sv = InputVariant::new("21".into(), 2000, 2000, b"N".to_vec(), b"<DEL>".to_vec());
        sv.sv_end = Some(5500);
        sv.is_structural = true;
        assert_eq!(nearest_values(&table, &sv, NearestType::Gene), ["ENSG_2"]);
    }

    /// The deletion `TC>T` trims to `C/-` at 1021, five bases above the site at
    /// 1016 and six below the site at 1027, so it reports the site at 1016
    /// alone. As the first feature of `TC>T,AC` it is measured there; as the
    /// second feature of `TC>TC,T`, whose first feature `-/-` at 1022-1021
    /// writes no row, it is measured at that feature's span, whose end at 1022
    /// is five bases below the site at 1027, so the two sites tie.
    #[test]
    fn a_rowless_first_feature_lends_its_span_to_the_record() {
        use std::collections::HashMap;

        use crate::runner::make_runner_test_transcript;
        use crate::transcript_index::TranscriptIndex;
        use crate::vcf_parser::{parse_vcf_line, split_minimal_record};

        let coding = |id: &str, start: u64| {
            let mut t = make_runner_test_transcript(id, "ENSG_X", "21", start);
            t.biotype = Arc::from("protein_coding");
            t
        };
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(
                vec![coding("ENST_BELOW", 1016), coding("ENST_ABOVE", 1027)],
                5000,
                5000,
            ),
        );
        let transcripts = crate::transcript_index::LazyTranscriptIndexes::from_indexes(by_chr);
        let kind = Some(NearestType::Transcript);

        let line = "21\t1020\t.\tTC\tT,AC\t.\t.\t.";
        let mut alleles = split_minimal_record(parse_vcf_line(line, true).unwrap());
        assert_eq!(alleles[0].allele_string, "C/-");
        annotate(&mut alleles[0], &transcripts, kind);
        assert_eq!(alleles[0].nearest, Some(vec!["ENST_BELOW".to_string()]));

        let line = "21\t1020\t.\tTC\tTC,T\t.\t.\t.";
        let alleles = split_minimal_record(parse_vcf_line(line, true).unwrap());
        let [mut deletion] = <[InputVariant; 1]>::try_from(alleles).unwrap();
        assert_eq!(
            (
                deletion.start,
                deletion.end,
                deletion.allele_string.as_str()
            ),
            (1021, 1021, "C/-")
        );
        annotate(&mut deletion, &transcripts, kind);
        assert_eq!(
            deletion.nearest,
            Some(vec!["ENST_BELOW".to_string(), "ENST_ABOVE".to_string()])
        );
    }

    #[test]
    fn nearest_type_parses_its_three_values() {
        assert_eq!("transcript".parse(), Ok(NearestType::Transcript));
        assert_eq!("gene".parse(), Ok(NearestType::Gene));
        assert_eq!("symbol".parse(), Ok(NearestType::Symbol));
        assert_eq!(
            "Symbol".parse::<NearestType>().unwrap_err(),
            "\"Symbol\" is not a valid value for --nearest\nValid values: transcript, gene, symbol"
        );
    }
}
