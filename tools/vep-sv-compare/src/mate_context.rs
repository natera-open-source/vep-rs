// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The terms a point breakend's mate-side Transcript row supports, derived from the cache
//! alone: read at the mate position, and read at the local position.
//!
//! A mate-side row is a Transcript row of a bracket-notation breakend whose Feature lies
//! within [`MAX_DISTANCE_FROM_TRANSCRIPT`] of the mate coordinate on the mate chromosome
//! while the local coordinate does not itself select that transcript (another chromosome,
//! or farther than that distance on the same one). Ensembl VEP builds such a row for every
//! transcript near the mate (ensembl-variation StructuralVariationOverlap.pm:86-102,
//! `_close_to_feature`) and then evaluates every positional predicate on it with the local
//! variation feature (BaseVariationFeatureOverlapAllele.pm:257,273 hand the local feature
//! to each predicate; `upstream` and `downstream` at Utils/VariationEffect.pm:444-457 and
//! `_bvfo_preds` at BaseVariationFeatureOverlapAllele.pm:454 read its start and end); only
//! `feature_truncation` (Utils/VariationEffect.pm:358) reads the mate breakend.
//!
//! [`breakend_context`] gives the term set a point breakend at a position supports on a
//! transcript, with `feature_truncation` decided by whether the mate lies inside the
//! transcript. Called with the mate position it is the set the mate breakend supports;
//! called with the local position it is Ensembl VEP's read, whose empty result becomes the
//! default `intergenic_variant` on the row. [`mate_side_derivation`] applies both to a row.
//!
//! The transcript structure comes from the JSON cache both engines annotated against:
//! exons, strand, biotype, the genomic CDS bounds mapped from the mapper's
//! `cdna_coding_start` and `cdna_coding_end` through its cDNA-to-genome pairs, and the
//! mature miRNA ranges of the `miRNA` attribute (cDNA coordinates).
//!
//! What the derivation does not model, so that a record reaching one of these edges stays
//! charged (a pair is set aside only when both engines' sets equal the derivations):
//!
//! - a breakend inside a start or stop codon: Ensembl VEP fires `start_lost` and
//!   `stop_lost` on span overlap (Utils/VariationEffect.pm:886-896, the stop analogue) and
//!   `start_retained_variant` through `_ins_del_start_altered` returning 0 on a structural
//!   allele (:1037);
//! - the frameshift-intron exon lookup: a transcript with an intron of 13 bases or fewer
//!   stretches its exon lookup by 12 bases (BaseTranscriptVariation.pm:911-914); the exon
//!   and intron tests here are exact;
//! - a multi-exon miRNA: the mature range is a cDNA interval that is not read off a
//!   multi-exon span, so the non-coding exon term stands;
//! - a same-chromosome breakend whose local coordinate is also within the gate of the
//!   transcript: the row is a local-side row and [`mate_side_derivation`] gives nothing.
//!
//! `_close_to_feature` expands the feature slice by 5,000 on each side and tests overlap
//! inclusively, as [`TranscriptModel::near`] does: a breakend exactly 5,000 bases outside
//! is inside the gate on both readings.
//!
//! A position written in a Location is read at any length, one past 128 bits as the
//! largest position there is; sums stop at the same bound.

use std::collections::BTreeSet;
use std::sync::LazyLock;

use regex::Regex;

use crate::failure::Failure;
pub use crate::locus::bracket_mate;
use crate::locus::{strip_chr, BREAKEND_FEATURE_DISTANCE};
use crate::pyfmt::py_strip;
use crate::text::{int_conversion_failure, INT_MAX_STR_DIGITS};

/// Ensembl VEP's breakend-to-feature distance (Utils/VariationEffect.pm:60, used by
/// StructuralVariationOverlap.pm:157), which is also its upstream and downstream window.
pub const MAX_DISTANCE_FROM_TRANSCRIPT: i128 = BREAKEND_FEATURE_DISTANCE;

/// The two bases at each end of an intron: a position inside them is no `intron_variant`.
const SPLICE_SITE: i128 = 2;

/// The polypyrimidine tract: the acceptor-side window 3 to 17 bases into the intron
/// (Utils/VariationEffect.pm `_intron_overlap`, the tract branch).
const PPT_INNER: i128 = 2;
const PPT_OUTER: i128 = 16;

/// A set of consequence terms.
pub type Terms = BTreeSet<&'static str>;

/// The structure of one cached transcript: exons, strand, biotype, genomic CDS bounds and
/// mature miRNA ranges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptModel {
    /// The chromosome the cache files the transcript under.
    pub chrom: String,
    pub start: i128,
    pub end: i128,
    pub strand: i128,
    pub biotype: String,
    /// The exons as `(start, end)`, ascending.
    pub exons: Vec<(i128, i128)>,
    /// The genomic bounds of the coding sequence, low first; `None` without one.
    pub cds: Option<(i128, i128)>,
    /// The mature miRNA ranges as `(first, last)` in cDNA coordinates.
    pub mirna_mature: Vec<(i128, i128)>,
}

impl TranscriptModel {
    fn intron_spans(&self) -> impl Iterator<Item = (i128, i128)> + '_ {
        self.exons.windows(2).filter_map(|pair| {
            let (start, end) = (pair[0].1.saturating_add(1), pair[1].0.saturating_sub(1));
            (end >= start).then_some((start, end))
        })
    }

    /// The gaps between neighbouring exons, as `(start, end)`; exons that abut or overlap
    /// leave none.
    pub fn introns(&self) -> Vec<(i128, i128)> {
        self.intron_spans().collect()
    }

    /// Whether `pos` is inside the transcript or within [`MAX_DISTANCE_FROM_TRANSCRIPT`] of
    /// either end.
    pub fn near(&self, pos: i128) -> bool {
        self.start.saturating_sub(MAX_DISTANCE_FROM_TRANSCRIPT) <= pos
            && pos <= self.end.saturating_add(MAX_DISTANCE_FROM_TRANSCRIPT)
    }

    pub fn inside(&self, pos: i128) -> bool {
        self.start <= pos && pos <= self.end
    }
}

/// One cDNA-to-genome pair of a transcript's mapper: the cDNA interval `from_start` to
/// `from_end` lies on the genomic interval `to_start` to `to_end`, in the orientation
/// `ori` (below zero for the reverse strand).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MapperPair {
    pub from_start: i128,
    pub from_end: i128,
    pub to_start: i128,
    pub to_end: i128,
    pub ori: i128,
}

/// One transcript record of a cache shard, reduced to what the comparison reads. A number
/// is `None` where the record gives no integer for it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CacheRecord {
    /// `stable_id` when it is a text; `None` for an id of another kind.
    pub stable_id: Option<String>,
    /// `stable_id` when it is a number other than zero or `true`: the text that stands for
    /// the number in a count of distinct ids. `None` for a text id.
    pub number_id: Option<String>,
    pub start: Option<i128>,
    pub end: Option<i128>,
    pub strand: Option<i128>,
    /// `biotype` when it is a text, else empty.
    pub biotype: String,
    /// The exons that give both ends as integers, in file order.
    pub exons: Vec<(i128, i128)>,
    /// The mapper's `cdna_coding_start` and `cdna_coding_end`.
    pub cdna_coding_start: Option<i128>,
    pub cdna_coding_end: Option<i128>,
    /// The mapper's pairs that give all five numbers as integers, in file order.
    pub pairs: Vec<MapperPair>,
    /// The text values of the attributes whose code is `miRNA`, in file order.
    pub mirna: Vec<String>,
}

/// The genomic bounds of the coding sequence, low first: the mapper's two cDNA coding
/// bounds, each carried through the pair that holds it (the last such pair, should several
/// hold it). A cDNA offset inside a pair lands at `to_start + offset` on the forward strand
/// and at `to_end - offset` on the reverse. `None` when a bound is missing, is not above
/// zero or lies in no pair.
pub fn genomic_cds(tx: &CacheRecord) -> Option<(i128, i128)> {
    let (cs, ce) = (tx.cdna_coding_start?, tx.cdna_coding_end?);
    if cs <= 0 || ce <= 0 {
        return None;
    }
    let (mut gs, mut ge) = (None, None);
    for pair in &tx.pairs {
        let genomic = |cdna: i128| {
            let offset = cdna.saturating_sub(pair.from_start);
            if pair.ori >= 0 {
                pair.to_start.saturating_add(offset)
            } else {
                pair.to_end.saturating_sub(offset)
            }
        };
        if pair.from_start <= cs && cs <= pair.from_end {
            gs = Some(genomic(cs));
        }
        if pair.from_start <= ce && ce <= pair.from_end {
            ge = Some(genomic(ce));
        }
    }
    let (gs, ge) = (gs?, ge?);
    Some(if gs <= ge { (gs, ge) } else { (ge, gs) })
}

/// The two runs of decimal digits of a mature range written `first-last`, blanks around it
/// dropped; `None` for any other text.
fn mature_range_digits(value: &str) -> Option<(&str, &str)> {
    let (first, last) = py_strip(value).split_once('-')?;
    let digits = |run: &str| !run.is_empty() && run.chars().all(|c| decimal_digit(c).is_some());
    (digits(first) && digits(last)).then_some((first, last))
}

/// The failure Python's `int` ends in on a mature range one of whose runs has more digits
/// than it converts, the first such run deciding.
pub(crate) fn mature_range_failure(value: &str) -> Option<Failure> {
    let (first, last) = mature_range_digits(value)?;
    [first, last]
        .into_iter()
        .map(|run| run.chars().count())
        .find(|digits| *digits > INT_MAX_STR_DIGITS)
        .map(int_conversion_failure)
}

/// The mature miRNA ranges of a record: each `miRNA` attribute value written `first-last`
/// in decimal digits, blanks around it dropped.
fn mirna_mature(tx: &CacheRecord) -> Vec<(i128, i128)> {
    tx.mirna
        .iter()
        .filter_map(|value| {
            let (first, last) = mature_range_digits(value)?;
            Some((decimal_value(first)?, decimal_value(last)?))
        })
        .collect()
}

/// The structure of a cache record filed under `chrom`; `None` for a record that gives no
/// integer start, end or strand.
pub fn model_from_cache_record(chrom: &str, tx: &CacheRecord) -> Option<TranscriptModel> {
    let (start, end, strand) = (tx.start?, tx.end?, tx.strand?);
    let mut exons = tx.exons.clone();
    exons.sort_unstable();
    Some(TranscriptModel {
        chrom: chrom.to_string(),
        start,
        end,
        strand,
        biotype: tx.biotype.clone(),
        exons,
        cds: genomic_cds(tx),
        mirna_mature: mirna_mature(tx),
    })
}

/// The value of a decimal digit of any script, the characters Python's `int` reads. The
/// decimal digits of a script are ten consecutive characters from its zero, and the runs of
/// neighbouring scripts follow each other in whole tens, so a digit is its distance from
/// the start of its run, modulo ten.
pub(crate) fn decimal_digit(c: char) -> Option<u32> {
    static DIGIT: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\d$").expect("regex"));
    let is_digit = |c: char| {
        let mut buffer = [0u8; 4];
        DIGIT.is_match(c.encode_utf8(&mut buffer))
    };
    if c.is_ascii() {
        return c.to_digit(10);
    }
    if !is_digit(c) {
        return None;
    }
    let mut zero = c as u32;
    while char::from_u32(zero - 1).is_some_and(is_digit) {
        zero -= 1;
    }
    Some((c as u32 - zero) % 10)
}

/// The value of a text of decimal digits of any script, one past 128 bits read as the
/// largest value; `None` for an empty text or one holding any other character.
pub(crate) fn decimal_value(digits: &str) -> Option<i128> {
    if digits.is_empty() {
        return None;
    }
    let mut value: i128 = 0;
    for c in digits.chars() {
        let digit = decimal_digit(c)?;
        value = value.saturating_mul(10).saturating_add(i128::from(digit));
    }
    Some(value)
}

/// `(chromosome, position)` of a point Location, the chromosome without a leading `chr`;
/// `None` for a ranged Location or any other text.
pub fn point_location(loc: &str) -> Option<(&str, i128)> {
    let (chrom, coords) = loc.split_once(':')?;
    Some((strip_chr(chrom), decimal_value(coords)?))
}

/// Whether `pos` lies in a mature range of a single-exon miRNA; a multi-exon miRNA has no
/// cDNA position to read off its span.
fn mature_mirna_hit(pos: i128, model: &TranscriptModel) -> bool {
    if model.biotype != "miRNA" || model.exons.len() != 1 {
        return false;
    }
    let cdna = if model.strand == 1 {
        pos.saturating_sub(model.start).saturating_add(1)
    } else {
        model.end.saturating_sub(pos).saturating_add(1)
    };
    model
        .mirna_mature
        .iter()
        .any(|&(first, last)| first <= cdna && cdna <= last)
}

/// The terms a point breakend at `pos`, inside the transcript, carries besides
/// `feature_truncation`: the exon terms (CDS, 5' or 3' UTR by strand, the mature-miRNA term
/// alone, or the non-coding exon term), or the intron terms (`intron_variant` strictly
/// inside the splice sites, the polypyrimidine tract window, and the non-coding transcript
/// term for a transcript without a CDS that is not NMD), plus `NMD_transcript_variant` on a
/// nonsense-mediated-decay transcript.
pub fn region_terms(pos: i128, model: &TranscriptModel) -> Terms {
    let mut terms = Terms::new();
    let nmd = model.biotype == "nonsense_mediated_decay";
    if model.exons.iter().any(|&(s, e)| s <= pos && pos <= e) {
        match model.cds {
            Some((cs, ce)) => {
                terms.insert(if cs <= pos && pos <= ce {
                    "coding_sequence_variant"
                } else if (model.strand == 1 && pos < cs) || (model.strand == -1 && pos > ce) {
                    "5_prime_UTR_variant"
                } else {
                    "3_prime_UTR_variant"
                });
            }
            None => {
                if mature_mirna_hit(pos, model) {
                    return Terms::from(["mature_miRNA_variant"]);
                }
                terms.insert("non_coding_transcript_exon_variant");
            }
        }
    } else if let Some((s, e)) = model.intron_spans().find(|&(s, e)| s <= pos && pos <= e) {
        if s.saturating_add(SPLICE_SITE) <= pos && pos <= e.saturating_sub(SPLICE_SITE) {
            terms.insert("intron_variant");
        }
        let in_tract = match model.strand {
            1 => e.saturating_sub(PPT_OUTER) <= pos && pos <= e.saturating_sub(PPT_INNER),
            -1 => s.saturating_add(PPT_INNER) <= pos && pos <= s.saturating_add(PPT_OUTER),
            _ => false,
        };
        if in_tract {
            terms.insert("splice_polypyrimidine_tract_variant");
        }
        if model.cds.is_none() && !nmd {
            terms.insert("non_coding_transcript_variant");
        }
    }
    if nmd {
        terms.insert("NMD_transcript_variant");
    }
    terms
}

/// The term set a point breakend evaluated at `pos` supports on the transcript.
///
/// `mate_inside` says whether the mate breakend lies inside the transcript, which is what
/// decides `feature_truncation` (Utils/VariationEffect.pm:358 reads the breakend, not the
/// evaluated position). Inside the transcript the region terms apply; within
/// [`MAX_DISTANCE_FROM_TRANSCRIPT`] outside it the strand-aware upstream or downstream
/// term; beyond that nothing but the truncation term. A mature-miRNA hit stands alone.
pub fn breakend_context(pos: i128, model: &TranscriptModel, mate_inside: bool) -> Terms {
    if model.inside(pos) {
        let mut terms = region_terms(pos, model);
        if mate_inside && terms != Terms::from(["mature_miRNA_variant"]) {
            terms.insert("feature_truncation");
        }
        return terms;
    }
    let mut terms = Terms::new();
    let forward = model.strand == 1;
    if model.start.saturating_sub(MAX_DISTANCE_FROM_TRANSCRIPT) <= pos && pos < model.start {
        terms.insert(if forward {
            "upstream_gene_variant"
        } else {
            "downstream_gene_variant"
        });
    } else if model.end < pos && pos <= model.end.saturating_add(MAX_DISTANCE_FROM_TRANSCRIPT) {
        terms.insert(if forward {
            "downstream_gene_variant"
        } else {
            "upstream_gene_variant"
        });
    }
    if mate_inside {
        terms.insert("feature_truncation");
    }
    terms
}

/// `(mate_set, local_set)` for a mate-side row, or `None` when the row is not one.
///
/// The row is a mate-side row when the Location is a point, the Allele is bracket notation,
/// the transcript is on the mate chromosome within the distance gate of the mate, and the
/// local coordinate does not itself select the transcript. `local_set` is Ensembl VEP's
/// read with its empty result replaced by `intergenic_variant`, the default a Transcript
/// row carries when no predicate fires.
pub fn mate_side_derivation(
    loc: &str,
    allele: &str,
    model: &TranscriptModel,
) -> Option<(Terms, Terms)> {
    let (local_chrom, local_pos) = point_location(loc)?;
    let (mate_chrom, mate_pos) = bracket_mate(allele)?;
    if model.chrom != mate_chrom || !model.near(mate_pos) {
        return None;
    }
    if local_chrom == model.chrom && model.near(local_pos) {
        return None;
    }
    let mate_inside = model.inside(mate_pos);
    let mate_set = breakend_context(mate_pos, model, mate_inside);
    let mut local_set = breakend_context(local_pos, model, mate_inside);
    if local_set.is_empty() {
        local_set.insert("intergenic_variant");
    }
    Some((mate_set, local_set))
}

#[cfg(test)]
mod tests {
    //! Every expected set is derived by hand from the transcript structure quoted in the
    //! test.

    use super::*;

    fn exons(spec: &str) -> Vec<(i128, i128)> {
        spec.split(';')
            .map(|exon| {
                let (start, end) = exon.split_once('-').expect("start-end");
                (start.parse().expect("start"), end.parse().expect("end"))
            })
            .collect()
    }

    fn model(
        chrom: &str,
        span: (i128, i128),
        strand: i128,
        biotype: &str,
        exon_spec: &str,
        cds: Option<(i128, i128)>,
        mirna_mature: &[(i128, i128)],
    ) -> TranscriptModel {
        TranscriptModel {
            chrom: chrom.to_string(),
            start: span.0,
            end: span.1,
            strand,
            biotype: biotype.to_string(),
            exons: exons(exon_spec),
            cds,
            mirna_mature: mirna_mature.to_vec(),
        }
    }

    /// ENST00000270142 (chromosome 21, forward strand, protein coding): five exons, CDS
    /// 33032083-33040891.
    fn app_like() -> TranscriptModel {
        model(
            "21",
            (33031935, 33041244),
            1,
            "protein_coding",
            "33031935-33032154;33036103-33036199;33038762-33038831;33039571-33039688;\
             33040784-33041244",
            Some((33032083, 33040891)),
            &[],
        )
    }

    /// ENST00000248935 GSTT1 (chromosome 22, reverse strand, protein coding), CDS
    /// 24376423-24384231.
    fn gstt1() -> TranscriptModel {
        model(
            "22",
            (24376133, 24384284),
            -1,
            "protein_coding",
            "24376133-24376617;24376822-24376998;24379361-24379511;24381700-24381787;\
             24384120-24384284",
            Some((24376423, 24384231)),
            &[],
        )
    }

    /// ENST00000265107 WDR70 (chromosome 5, forward strand, protein coding), CDS
    /// 37379368-37752573.
    fn wdr70() -> TranscriptModel {
        model(
            "5",
            (37379318, 37753435),
            1,
            "protein_coding",
            "37379318-37379392;37379489-37379554;37381602-37381685;37392000-37392120;\
             37396375-37396570;37437922-37437981;37443239-37443372;37479834-37479987;\
             37516514-37516590;37605064-37605238;37697655-37697754;37701058-37701142;\
             37702949-37703087;37721115-37721215;37722855-37722934;37724934-37725050;\
             37726883-37727045;37752486-37753435",
            Some((37379368, 37752573)),
            &[],
        )
    }

    /// ENST00000385128 MIR125B2 (chromosome 21, forward strand, miRNA), one exon, mature
    /// 17-38 and 54-75.
    fn mir() -> TranscriptModel {
        model(
            "21",
            (16590237, 16590325),
            1,
            "miRNA",
            "16590237-16590325",
            None,
            &[(17, 38), (54, 75)],
        )
    }

    /// ENST00000687248 PTCHD1-AS (chromosome X, reverse strand, lncRNA), nine exons, no
    /// CDS.
    fn ptchd1_as() -> TranscriptModel {
        model(
            "X",
            (22227405, 23293174),
            -1,
            "lncRNA",
            "22227405-22227667;22253939-22254065;22258013-22258136;22259755-22259842;\
             22389663-22389720;22743229-22743338;23064038-23064121;23270051-23270176;\
             23293014-23293174",
            None,
            &[],
        )
    }

    fn terms(names: &[&'static str]) -> Terms {
        names.iter().copied().collect()
    }

    /// The one record of a shard written as a JSON object.
    fn record(json: &str) -> CacheRecord {
        let shard = crate::cache::parse_shard(&format!("[{json}]"));
        let mut records = shard.expect("a record Python reads").expect("a shard");
        assert_eq!(records.len(), 1);
        records.remove(0)
    }

    mod model_tests {
        use super::*;

        #[test]
        fn test_introns_are_the_gaps_between_sorted_exons() {
            assert_eq!(app_like().introns()[0], (33032155, 33036102));
            assert_eq!(app_like().introns().len(), 4);
            assert_eq!(mir().introns(), []);
        }

        /// Forward: cDNA 149 in a first exon of 220 bases starting at 33031935 is 33032083.
        /// Reverse: cDNA 1 of a reverse exon 24384120-24384284 is 24384284, so a coding
        /// start at cDNA 54 is 24384231.
        #[test]
        fn test_genomic_cds_is_mapped_through_the_pairs_in_both_orientations() {
            let fwd = record(
                r#"{"stable_id": "fwd", "variation_effect_feature_cache": {"mapper": {
                "cdna_coding_start": 149, "cdna_coding_end": 613,
                "pairs": [
                 {"from_start": 1, "from_end": 220, "to_start": 33031935, "to_end": 33032154, "ori": 1},
                 {"from_start": 221, "from_end": 317, "to_start": 33036103, "to_end": 33036199, "ori": 1},
                 {"from_start": 318, "from_end": 387, "to_start": 33038762, "to_end": 33038831, "ori": 1},
                 {"from_start": 388, "from_end": 505, "to_start": 33039571, "to_end": 33039688, "ori": 1},
                 {"from_start": 506, "from_end": 966, "to_start": 33040784, "to_end": 33041244, "ori": 1}]}}}"#,
            );
            assert_eq!(genomic_cds(&fwd), Some((33032083, 33040891)));
            let rev = record(
                r#"{"stable_id": "rev", "variation_effect_feature_cache": {"mapper": {
                "cdna_coding_start": 54, "cdna_coding_end": 195,
                "pairs": [
                 {"from_start": 1, "from_end": 165, "to_start": 24384120, "to_end": 24384284, "ori": -1},
                 {"from_start": 166, "from_end": 253, "to_start": 24381700, "to_end": 24381787, "ori": -1}]}}}"#,
            );
            assert_eq!(genomic_cds(&rev), Some((24381758, 24384231)));
        }

        #[test]
        fn test_a_transcript_without_coding_bounds_has_no_cds() {
            let no_bounds = record(
                r#"{"stable_id": "t", "variation_effect_feature_cache": {"mapper": {"pairs": []}}}"#,
            );
            assert_eq!(genomic_cds(&no_bounds), None);
            assert_eq!(genomic_cds(&record(r#"{"stable_id": "t"}"#)), None);
        }
    }

    mod region_term_tests {
        use super::*;

        #[test]
        fn test_cds_utr_by_strand() {
            let (app, gstt1) = (app_like(), gstt1());
            assert_eq!(
                region_terms(33032100, &app),
                terms(&["coding_sequence_variant"])
            );
            assert_eq!(
                region_terms(33032000, &app),
                terms(&["5_prime_UTR_variant"])
            );
            assert_eq!(
                region_terms(33041000, &app),
                terms(&["3_prime_UTR_variant"])
            );
            // Reverse strand: bases beyond the CDS end are the 5' UTR.
            assert_eq!(
                region_terms(24384250, &gstt1),
                terms(&["5_prime_UTR_variant"])
            );
            assert_eq!(
                region_terms(24376300, &gstt1),
                terms(&["3_prime_UTR_variant"])
            );
        }

        /// Inside intron 1 of the forward transcript (33032155-33036102): the first two and
        /// last two bases are the splice sites and carry no intron term; the forward-strand
        /// polypyrimidine window is the acceptor side, 16 to 2 bases before the intron end.
        #[test]
        fn test_intron_terms() {
            let app = app_like();
            assert_eq!(region_terms(33034355, &app), terms(&["intron_variant"]));
            assert_eq!(region_terms(33032155, &app), terms(&[]));
            assert_eq!(region_terms(33036102, &app), terms(&[]));
            assert_eq!(
                region_terms(33036102 - 5, &app),
                terms(&["intron_variant", "splice_polypyrimidine_tract_variant"])
            );
            assert_eq!(
                region_terms(33036102 - 17, &app),
                terms(&["intron_variant"])
            );
            // Reverse strand: the acceptor is the intron start side.
            assert_eq!(
                region_terms(24376618 + 5, &gstt1()),
                terms(&["intron_variant", "splice_polypyrimidine_tract_variant"])
            );
        }

        #[test]
        fn test_non_coding_intron_and_exon() {
            let lnc = ptchd1_as();
            assert_eq!(
                region_terms(22500000, &lnc),
                terms(&["intron_variant", "non_coding_transcript_variant"])
            );
            assert_eq!(
                region_terms(22227500, &lnc),
                terms(&["non_coding_transcript_exon_variant"])
            );
        }

        #[test]
        fn test_nmd_transcript_adds_its_term() {
            let nmd = model(
                "21",
                (100, 1000),
                1,
                "nonsense_mediated_decay",
                "100-200;800-1000",
                Some((150, 900)),
                &[],
            );
            assert_eq!(
                region_terms(500, &nmd),
                terms(&["intron_variant", "NMD_transcript_variant"])
            );
            assert_eq!(
                region_terms(160, &nmd),
                terms(&["coding_sequence_variant", "NMD_transcript_variant"])
            );
        }

        /// cDNA 17-38 of a forward single-exon miRNA starting at 16590237 is
        /// 16590253-16590274.
        #[test]
        fn test_mature_mirna_hit_stands_alone() {
            let mir = mir();
            assert_eq!(
                region_terms(16590260, &mir),
                terms(&["mature_miRNA_variant"])
            );
            assert_eq!(
                region_terms(16590240, &mir),
                terms(&["non_coding_transcript_exon_variant"])
            );
        }
    }

    mod breakend_context_tests {
        use super::*;

        #[test]
        fn test_inside_with_the_mate_inside_carries_feature_truncation() {
            assert_eq!(
                breakend_context(33034355, &app_like(), true),
                terms(&["feature_truncation", "intron_variant"])
            );
        }

        #[test]
        fn test_upstream_downstream_are_strand_aware() {
            // 4,108 bases before the reverse-strand transcript's start: downstream.
            assert_eq!(
                breakend_context(24372025, &gstt1(), false),
                terms(&["downstream_gene_variant"])
            );
            assert_eq!(
                breakend_context(24384284 + 100, &gstt1(), false),
                terms(&["upstream_gene_variant"])
            );
            assert_eq!(
                breakend_context(33031935 - 100, &app_like(), false),
                terms(&["upstream_gene_variant"])
            );
        }

        #[test]
        fn test_beyond_five_kb_carries_only_the_truncation_term_when_the_mate_is_inside() {
            assert_eq!(
                breakend_context(14374334, &app_like(), true),
                terms(&["feature_truncation"])
            );
            assert_eq!(breakend_context(14374334, &app_like(), false), terms(&[]));
        }

        #[test]
        fn test_a_mature_mirna_hit_takes_no_truncation_term() {
            assert_eq!(
                breakend_context(16590260, &mir(), true),
                terms(&["mature_miRNA_variant"])
            );
        }
    }

    mod mate_side_derivation_tests {
        use super::*;

        /// `1:14374334 ]21:33034355]A` on ENST00000270142: the mate sits in intron 1, the
        /// local coordinate on chromosome 1 is nowhere near, so Ensembl VEP writes only
        /// `feature_truncation` and the mate supports `feature_truncation,intron_variant`.
        #[test]
        fn test_cross_chromosome_row_perl_read_at_the_local_coordinate() {
            let got = mate_side_derivation("1:14374334", "]21:33034355]A", &app_like());
            assert_eq!(
                got,
                Some((
                    terms(&["feature_truncation", "intron_variant"]),
                    terms(&["feature_truncation"])
                ))
            );
        }

        /// `21:37392051 [5:37717217[A` on WDR70: the mate is in intron 13
        /// (37703088-37721114), the local coordinate lands numerically in exon 4
        /// (37392000-37392120), inside the CDS: Ensembl VEP's read is
        /// `feature_truncation,coding_sequence_variant`.
        #[test]
        fn test_the_local_coordinate_landing_in_a_cds_exon_of_the_mate_transcript() {
            let got = mate_side_derivation("21:37392051", "[5:37717217[A", &wdr70());
            assert_eq!(
                got,
                Some((
                    terms(&["feature_truncation", "intron_variant"]),
                    terms(&["feature_truncation", "coding_sequence_variant"])
                ))
            );
        }

        #[test]
        fn test_a_far_local_coordinate_leaves_only_the_truncation_term() {
            let got = mate_side_derivation("21:5100001", "[5:37717217[T", &wdr70());
            assert_eq!(
                got.expect("a mate-side row").1,
                terms(&["feature_truncation"])
            );
        }

        /// `21:9500001 N[22:24372025[` on GSTT1: the mate is 4,108 bases downstream, the
        /// local coordinate fires nothing, and Ensembl VEP writes `intergenic_variant` on
        /// the Transcript row.
        #[test]
        fn test_a_far_local_coordinate_beside_an_outside_mate_defaults_to_intergenic() {
            let got = mate_side_derivation("21:9500001", "N[22:24372025[", &gstt1());
            assert_eq!(
                got,
                Some((
                    terms(&["downstream_gene_variant"]),
                    terms(&["intergenic_variant"])
                ))
            );
        }

        /// `21:24372025 N[22:24372025[`: the two reads coincide, which is why Ensembl VEP's
        /// row is right on this record.
        #[test]
        fn test_equal_numeric_coordinates_give_equal_sets() {
            let (mate_set, local_set) =
                mate_side_derivation("21:24372025", "N[22:24372025[", &gstt1())
                    .expect("a mate-side row");
            assert_eq!(mate_set, local_set);
            assert_eq!(mate_set, terms(&["downstream_gene_variant"]));
        }

        #[test]
        fn test_a_row_the_local_coordinate_selects_is_not_a_mate_side_row() {
            assert_eq!(
                mate_side_derivation("21:33034000", "N[21:33034355[", &app_like()),
                None
            );
        }

        #[test]
        fn test_a_ranged_location_or_a_non_bracket_allele_is_not_a_mate_side_row() {
            assert_eq!(
                mate_side_derivation("21:14604132-41770788", "N[21:41770788[", &app_like()),
                None
            );
            assert_eq!(
                mate_side_derivation("1:14374334", "deletion", &app_like()),
                None
            );
        }

        #[test]
        fn test_a_mate_beyond_five_kb_of_the_transcript_is_not_a_mate_side_row() {
            assert_eq!(
                mate_side_derivation("1:14374334", "]21:33050000]A", &app_like()),
                None
            );
        }

        #[test]
        fn test_a_transcript_on_another_chromosome_than_the_mate_is_not_a_mate_side_row() {
            assert_eq!(
                mate_side_derivation("1:14374334", "]22:33034355]A", &app_like()),
                None
            );
        }

        #[test]
        fn test_helpers() {
            assert_eq!(bracket_mate("N[chr21:41770788["), Some(("21", 41770788)));
            assert_eq!(bracket_mate("]13:123456]T"), Some(("13", 123456)));
            assert_eq!(bracket_mate("N."), None);
            assert_eq!(point_location("chr21:5"), Some(("21", 5)));
            assert_eq!(point_location("21:1-2"), None);
        }
    }

    #[test]
    fn a_position_is_read_in_the_digits_of_any_script_and_at_any_length() {
        assert_eq!(point_location("21:\u{0663}\u{0660}"), Some(("21", 30)));
        assert_eq!(point_location("21:\u{ff11}2"), Some(("21", 12)));
        assert_eq!(point_location("21:"), None);
        assert_eq!(point_location("21:5:6"), None);
        assert_eq!(point_location("21:\u{b2}"), None);
        assert_eq!(point_location("21"), None);
        let long = format!("21:{}", "9".repeat(60));
        assert_eq!(point_location(&long), Some(("21", i128::MAX)));
        let far = mate_side_derivation(&long, "]21:33034355]A", &app_like());
        assert_eq!(
            far.expect("a mate-side row").1,
            terms(&["feature_truncation"])
        );
    }

    #[test]
    fn a_mature_range_is_two_runs_of_digits_around_one_hyphen() {
        let ranges = |values: &[&str]| {
            mirna_mature(&CacheRecord {
                mirna: values.iter().map(|value| value.to_string()).collect(),
                ..CacheRecord::default()
            })
        };
        assert_eq!(
            ranges(&["3-7", " 17-38\n", "\u{1c}54-75\u{a0}"]),
            [(3, 7), (17, 38), (54, 75)]
        );
        assert_eq!(ranges(&["3", "3-", "-7", "3-7-9", "3 -7", "a-7", ""]), []);
        assert_eq!(ranges(&["\u{0663}-\u{ff17}"]), [(3, 7)]);
        let long = "1".repeat(INT_MAX_STR_DIGITS + 1);
        let limit = |digits: usize| Some(int_conversion_failure(digits));
        assert_eq!(mature_range_failure(&format!("{long}-5")), limit(4301));
        assert_eq!(mature_range_failure(&format!(" 5-{long}7\n")), limit(4302));
        assert_eq!(
            mature_range_failure(&format!("{long}7-{long}")),
            limit(4302)
        );
        assert_eq!(mature_range_failure(&format!("{}-5", &long[1..])), None);
        assert_eq!(mature_range_failure(&format!("{long}-x")), None);
        assert_eq!(mature_range_failure(&long), None);
        let other_script = "\u{0663}".repeat(INT_MAX_STR_DIGITS + 1);
        assert_eq!(
            mature_range_failure(&format!("{other_script}-5")),
            limit(4301)
        );
    }

    #[test]
    fn a_record_without_an_integer_start_end_or_strand_has_no_model() {
        let whole = CacheRecord {
            start: Some(10),
            end: Some(20),
            strand: Some(-1),
            exons: vec![(15, 20), (10, 12)],
            ..CacheRecord::default()
        };
        let model = model_from_cache_record("21", &whole).expect("a model");
        assert_eq!(model.exons, [(10, 12), (15, 20)]);
        assert_eq!((model.chrom.as_str(), model.cds), ("21", None));
        let lacking = [
            CacheRecord {
                start: None,
                ..whole.clone()
            },
            CacheRecord {
                end: None,
                ..whole.clone()
            },
            CacheRecord {
                strand: None,
                ..whole.clone()
            },
        ];
        for record in &lacking {
            assert_eq!(model_from_cache_record("21", record), None);
        }
    }
}
