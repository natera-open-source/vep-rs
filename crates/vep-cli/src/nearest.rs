// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--nearest`: the protein-coding transcription start sites nearest a variant
//! and the identifiers reported for them, as Ensembl VEP's `TranscriptTree`
//! finds them (`TranscriptTree.pm:222-254`, `AnnotationType/Transcript.pm:700-715`
//! at release/116.2).
//!
//! Ensembl VEP keeps every protein-coding transcript's start site in an interval
//! tree, one point per transcript, and asks the tree for one site on each side of
//! the variant: `fetch_nearest_down(end)` and `fetch_nearest_up(start)`. The
//! distance of a site is the smallest difference between it and either end of the
//! variant (`_get_dist`, `:303-317`); every candidate at the minimum is reported,
//! the one found below the end first. Here the sites are the chromosome's sorted
//! [`TssPoint`] table and the two searches are binary searches.
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
use crate::transcript_index::{LazyTranscriptIndexes, TssPoint};

/// The identifier `--nearest` reports for a transcription start site
/// (`Config.pm:602`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
/// in stable-id order, the candidate found below the end before the one found
/// above the start (`TranscriptTree.pm:233-236` lists them in that order, and
/// both can lie inside a span).
///
/// The tree holds a site at `tss` as the half-open interval `[tss - 1, tss)`.
/// `fetch_nearest_down(end)` returns the interval with the highest upper bound
/// at or below `end`, so the highest site at or below the span's end;
/// `fetch_nearest_up(start)` returns the interval with the lowest lower bound
/// above `start`, so the lowest site at or above `start + 2`: a site at
/// `start + 1` is the interval `[start, start + 1)`, which contains the query
/// and is never returned. Each candidate's distance is the smaller of its
/// differences to the two ends; both are kept when they tie.
pub fn nearest_points(points: &[TssPoint], start: u64, end: u64) -> Vec<&TssPoint> {
    let (low, high) = (start.min(end), start.max(end));
    let down = points[..points.partition_point(|p| p.tss <= high)].last();
    let up = points.get(points.partition_point(|p| p.tss < low.saturating_add(2)));
    let mut candidates: Vec<u64> = Vec::with_capacity(2);
    for site in down.into_iter().chain(up) {
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
pub fn nearest_values(
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
/// rows print (OutputFactory.pm:2286-2296, :1017).
pub(crate) fn annotate(
    variant: &mut InputVariant,
    transcripts: &LazyTranscriptIndexes,
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
        // at or below the end (5005) before the lowest at or above start + 2 (4995).
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
    /// below the end and by `up` when it is the lowest at or above `start + 2`;
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

    /// The tree holds a site at `t` as `[t - 1, t)`, and `fetch_nearest_up(s)`
    /// skips the interval containing `s`: a site at `s + 1` is never the `up`
    /// candidate, so a variant one base below a start site reports the site
    /// below it instead.
    #[test]
    fn a_site_one_base_above_the_start_is_not_the_up_candidate() {
        let table = table();
        // 4999: the site at 5000 is [4999, 5000), which contains 4999; the up
        // candidate is 9000 and the down candidate 1000 (distance 3999) wins.
        assert_eq!(
            ids(&nearest_points(&table, 4999, 4999)),
            ["ENST_A", "ENST_B"]
        );
        // 4998: the site at 5000 is above 4998 + 1 and wins at distance 2.
        assert_eq!(ids(&nearest_points(&table, 4998, 4998)), ["ENST_R"]);
        // A two-base span ending at 5000 reaches the site from below.
        assert_eq!(ids(&nearest_points(&table, 4999, 5000)), ["ENST_R"]);
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

    /// The deletion `TC>T` trims to `C/-` at 1021, a span that ends one base
    /// below the site at 1022 and so reaches the site at 996 alone. As the
    /// first feature of `TC>T,AC` it is measured there; as the second feature of
    /// `TC>TC,T`, whose first feature `-/-` at 1022-1021 writes no row, it is
    /// measured at that feature's span, which the site at 1022 ends.
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
                vec![coding("ENST_BELOW", 996), coding("ENST_AT", 1022)],
                5000,
                5000,
            ),
        );
        let transcripts = LazyTranscriptIndexes::from_indexes(by_chr);
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
        assert_eq!(deletion.nearest, Some(vec!["ENST_AT".to_string()]));
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
