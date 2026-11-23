// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `RegulatoryFeature` and `MotifFeature` rows a variant gains from the
//! features it overlaps (`AnnotationType/RegFeat.pm:114-150`), with the six
//! regulatory terms of `Utils/Constants.pm:1271-1372` and the motif fields of
//! `OutputFactory.pm:1909-1952`.

use std::sync::Arc;

use smallvec::SmallVec;
use vep_core::consequence::{Activity, CellTypeActivity, Consequence, TranscriptConsequence};
use vep_core::variant::{InputVariant, VariantClass};
use vep_effects::regulatory::{motif_fields, Allele, Motif};

use super::source::{Feature, FeatureKind, RegulatorySource};
use crate::config::CliError;

/// The epigenomes `--cell_type` named, resolved against the cache's cell-type
/// list: each name with its position in the list, in the order given.
#[derive(Debug)]
pub struct CellTypes {
    requested: Vec<(Arc<str>, usize)>,
}

impl CellTypes {
    /// Reads the `--cell_type` value against the list of the regulatory source
    /// the run opened ([`crate::config::split_cell_types`]), as Ensembl VEP
    /// checks the names when it builds its regulatory source
    /// (`Cache/RegFeat.pm:113-117`): `None` without the flag, and without a
    /// source, a cache with no regulatory data, where Ensembl VEP builds no
    /// source and checks nothing. Under `--regulatory_gff` the cache source is
    /// not built either (`CacheDir.pm:156-176`), so the GFF source's own
    /// warning stands for the flag, no name is checked, and the list is empty:
    /// a GFF feature carries no activity map, so `get_cell_types` keeps none of
    /// the requested names (OutputFactory.pm:1968-1978) and every regulatory
    /// and motif row carries `CELL_TYPE` as the empty list.
    pub(crate) fn resolve(
        value: Option<&str>,
        source: Option<&RegulatorySource>,
    ) -> Result<Option<Arc<Self>>, CliError> {
        let (Some(value), Some(source)) = (value, source) else {
            return Ok(None);
        };
        let requested = if source.reads_gff() {
            Vec::new()
        } else {
            let available = source.cell_types();
            crate::config::split_cell_types(value, available)?
                .into_iter()
                .map(|index| (Arc::from(available[index].as_str()), index))
                .collect()
        };
        Ok(Some(Arc::new(Self { requested })))
    }

    /// The `CELL_TYPE` list of a regulatory feature whose activity string is
    /// `activity`: each requested epigenome with the feature's activity in it,
    /// an epigenome the feature's map lacks (`-`) skipped, as `get_cell_types`
    /// keeps the requested names the map has a value for (OutputFactory.pm:1974).
    fn project(&self, activity: &str) -> Box<[CellTypeActivity]> {
        let codes = activity.as_bytes();
        self.requested
            .iter()
            .filter_map(|(name, index)| {
                let activity = match codes.get(*index) {
                    Some(b'A') => Activity::Active,
                    Some(b'I') => Activity::Inactive,
                    Some(b'P') => Activity::Poised,
                    Some(b'R') => Activity::Repressed,
                    Some(b'N') => Activity::Na,
                    _ => return None,
                };
                Some(CellTypeActivity {
                    name: Arc::clone(name),
                    activity,
                })
            })
            .collect()
    }
}

/// The record a query builds rows for, with the `--cell_type` list its
/// regulatory and motif rows carry (`None` without the flag).
#[derive(Clone, Copy)]
struct RowSubject<'a> {
    variant: &'a InputVariant,
    cell_types: Option<&'a CellTypes>,
}

impl RegulatorySource {
    /// Annotates a variant's regulatory and motif overlaps after its transcript
    /// rows, then, where the record gained its companion intergenic row, its
    /// `--shift_3prime` shift, which the transcript path computes for an
    /// intergenic verdict alone. Both annotation paths call this one method.
    pub(crate) fn annotate_variant(
        &self,
        variant: &mut InputVariant,
        effects: &vep_effects::EffectsConfig,
        cell_types: Option<&CellTypes>,
    ) {
        self.annotate_rows(variant, cell_types);
        if variant.intergenic_overlap {
            variant.intergenic_shift = crate::runner::intergenic_shift(variant, effects);
        }
    }

    /// [`Self::annotate_rows`] of a run without `--cell_type`.
    #[cfg(test)]
    pub(crate) fn annotate(&self, variant: &mut InputVariant) {
        self.annotate_rows(variant, None);
    }

    /// Adds the variant's regulatory and motif rows; under `--cell_type` each
    /// carries its `CELL_TYPE` list, a regulatory row the projection of the
    /// requested epigenomes onto its feature's activity and a motif row the
    /// empty list (OutputFactory.pm:1888, 1936, 2107-2108).
    pub(crate) fn annotate_rows(&self, variant: &mut InputVariant, cell_types: Option<&CellTypes>) {
        if variant.no_allele || variant.variant_class == VariantClass::ReferenceBlock {
            return;
        }
        let mut rows = Vec::new();
        let subject = RowSubject {
            variant,
            cell_types,
        };
        match breakend_mate(variant) {
            // The bracket allele of a paired breakend is built for a feature only
            // when its mate lies within MAX_DISTANCE_FROM_TRANSCRIPT of the feature
            // on the feature's chromosome (StructuralVariationOverlap.pm:90-102,
            // `_close_to_feature` `:144-160`), so a mate on another chromosome has
            // no row against the features its own span overlaps; and each
            // breakend position is a query of its own (InputBuffer.pm:335-342),
            // matched on its chromosome.
            Some((mate_chr, mate_pos)) => {
                if mate_chr == variant.chr {
                    self.rows_for_span(
                        subject,
                        &variant.chr,
                        (variant.start, variant.end),
                        &|feature| feature.close_to(mate_pos),
                        RowTerms::Overlap,
                        &mut rows,
                    );
                }
                self.rows_for_span(
                    subject,
                    mate_chr,
                    (mate_pos, mate_pos),
                    &|_| true,
                    RowTerms::Overlap,
                    &mut rows,
                );
            }
            None => {
                self.rows_for_span(
                    subject,
                    &variant.chr,
                    (variant.start, variant.end),
                    &|_| true,
                    RowTerms::Overlap,
                    &mut rows,
                );
                // The record's own breakend is added to every feature one of its
                // mates lies in when its own position is close to the feature on
                // the same chromosome (`StructuralVariationOverlap.pm:90-102` walks
                // `($vf, @$breakends)` for each feature a breakend found). Its own
                // position lies outside such a feature (an overlapped one was found
                // above), so no overlap predicate holds and the row carries the
                // default term (`BaseVariationFeatureOverlapAllele.pm:285`,
                // `intergenic_variant`).
                for (mate_chr, mate_pos) in &variant.breakend_mates {
                    if *mate_chr == variant.chr {
                        self.rows_for_span(
                            subject,
                            mate_chr,
                            (*mate_pos, *mate_pos),
                            &|feature| feature.close_to(variant.start),
                            RowTerms::Default,
                            &mut rows,
                        );
                    }
                }
            }
        }
        if rows.is_empty() {
            return;
        }
        if variant.transcript_consequences.is_empty() {
            variant.intergenic_overlap = true;
        }
        // The transcript path left the record's verdict at the least rank of its
        // rows; the added rows fold into it, so a record that gained no row, the
        // common case, is not rescanned.
        let worst_added = rows
            .iter()
            .flat_map(|r| &r.consequences)
            .copied()
            .min_by_key(|c| c.rank());
        variant.transcript_consequences.extend(rows);
        if let Some(worst_added) = worst_added {
            let outranks_current = variant
                .most_severe_consequence
                .is_none_or(|current| worst_added.rank() < current.rank());
            if outranks_current {
                variant.most_severe_consequence = Some(worst_added);
            }
        }
    }

    /// One row per feature of `chr` overlapping the span that `wanted`
    /// accepts, appended to `rows`; a feature an earlier query left in `rows`,
    /// known by its type, stable id and span, is not repeated (one query yields
    /// each feature once, so the first query of a variant scans nothing).
    /// Ensembl VEP writes a breakend's regulatory row once per breakend
    /// position lying in the feature (`InputBuffer.pm:323-342` lists the
    /// variant once per overlapping position and `RegFeat.pm:122-135` builds
    /// one overlap per listing); vep-rs writes the feature once.
    /// Two distinct features of one file sharing a stable id, two records
    /// under one `ID`, give a sequence variant one regulatory row, the record
    /// later in the file: Ensembl VEP keys the variant's regulatory overlaps on
    /// the stable id and the later feature's write wins
    /// (`VariationFeature.pm:913-919`, the features arriving in file order), so
    /// the feature of the greater ordinal keeps the row whichever the
    /// descending-start walk yields first (a promoter widened under
    /// `--extended_promoters` sorts by its widened start, not its line). A
    /// structural variant keeps both rows (`AnnotationType/RegFeat.pm:128-136`
    /// lists its overlaps), as does a motif pair
    /// (`MotifFeatureVariation.pm:70-76` keys on the database id), and features
    /// of two files sharing an identifier are two rows, since an identifier is
    /// unique within its file alone.
    /// The terms are [`vep_effects::regulatory::regulatory_terms`] for a row
    /// whose allele overlaps the feature: all six are tier 2
    /// (Utils/Constants.pm:1285-1369), so the tier loop
    /// (BaseVariationFeatureOverlapAllele.pm:267-279) writes the ablation or
    /// amplification term and the overlap term together; a row the allele does
    /// not overlap carries the default term alone.
    fn rows_for_span(
        &self,
        subject: RowSubject<'_>,
        chr: &str,
        (start, end): (u64, u64),
        wanted: &dyn Fn(&Feature) -> bool,
        terms: RowTerms,
        rows: &mut Vec<TranscriptConsequence>,
    ) {
        let Some(index) = self.index(chr) else {
            return;
        };
        let RowSubject {
            variant,
            cell_types,
        } = subject;
        let earlier = rows.len();
        // The regulatory features this query wrote a row for, each with its
        // row's index: a second feature of the same file and stable id is the
        // same feature to a sequence variant, and the record later in the file
        // keeps the row.
        let mut written: SmallVec<[(&Feature, usize); 4]> = SmallVec::new();
        for (feature, ranks) in index.overlapping(start, end) {
            let feature_type = feature.feature_type();
            if !wanted(feature)
                || rows[..earlier].iter().any(|row| {
                    row.feature_type == feature_type
                        && row.transcript_id() == &*feature.stable_id
                        && (row.feature_start, row.feature_end) == (feature.start, feature.end)
                })
            {
                continue;
            }
            let mut replaces = None;
            if matches!(feature.kind, FeatureKind::Regulatory { .. }) && !variant.is_structural {
                match written.iter_mut().find(|(kept, _)| {
                    kept.file == feature.file && kept.stable_id == feature.stable_id
                }) {
                    Some((kept, at)) if kept.ordinal < feature.ordinal => {
                        *kept = feature;
                        replaces = Some(*at);
                    }
                    Some(_) => continue,
                    None => written.push((feature, rows.len())),
                }
            }
            let consequences = match terms {
                RowTerms::Overlap => vep_effects::regulatory::regulatory_terms(
                    variant,
                    start,
                    end,
                    feature_type,
                    feature.start,
                    feature.end,
                ),
                RowTerms::Default => smallvec::smallvec![Consequence::IntergenicVariant],
            };
            let mut row = TranscriptConsequence::default();
            row.feature_id = Some(Arc::clone(&feature.stable_id));
            row.feature_type = feature_type;
            // StructuralVariationOverlapAllele_to_output_hash does not call the
            // regulatory/motif sequence-allele factories, so a structural row has
            // no BIOTYPE and no motif fields.
            row.feature_biotype = match &feature.kind {
                FeatureKind::Regulatory { biotype, .. } if !variant.is_structural => {
                    biotype.clone()
                }
                _ => None,
            };
            row.strand = feature.strand;
            row.feature_start = feature.start;
            row.feature_end = feature.end;
            row.impact = consequences[0].impact();
            row.consequences = consequences;
            // A structural variant's rows come in the order the features were
            // met (`StructuralVariationFeature.pm:495-509` reads the lists
            // `AnnotationType/RegFeat.pm:135` fills as pushed); a sequence
            // variant's motif rows by the motif's database id.
            row.id_rank = if variant.is_structural {
                ranks.met
            } else {
                ranks.motif
            };
            row.motif_db_id = match feature.kind {
                FeatureKind::Motif { .. } => feature.dbid,
                FeatureKind::Regulatory { .. } => None,
            };
            row.cell_type = cell_types.map(|requested| match &feature.kind {
                FeatureKind::Regulatory { activity, .. } => requested.project(activity),
                FeatureKind::Motif { .. } => Box::default(),
            });
            if let FeatureKind::Motif { seq, matrix, .. } = &feature.kind {
                if !variant.is_structural {
                    let motif = Motif::new(
                        feature.start,
                        feature.end,
                        feature.strand,
                        seq.as_deref(),
                        matrix,
                    );
                    // The alleles on the forward strand of the genome: a record
                    // given on the reverse strand is complemented onto it, then
                    // `feature_seq` onto a reverse motif, together Ensembl VEP's
                    // one complement when the record's strand differs from the
                    // motif's (VariationFeatureOverlapAllele.pm:253-256).
                    let (reference, alternate) = variant.forward_strand_alleles();
                    let allele = Allele::new(start, end, &reference, &alternate);
                    row.motif = Some(Box::new(motif_fields(&motif, &allele)));
                }
            }
            match replaces {
                Some(at) => rows[at] = row,
                None => rows.push(row),
            }
        }
    }
}

/// How a row's terms are derived: from the allele's overlap of the feature, or
/// the default term of a row built for a feature the allele lies near without
/// overlapping it.
#[derive(Clone, Copy)]
enum RowTerms {
    Overlap,
    Default,
}

/// The distance within which a breakend is close to a feature
/// (`MAX_DISTANCE_FROM_TRANSCRIPT`, Utils/VariationEffect.pm:60, read by
/// StructuralVariationOverlap.pm:157 for every feature kind).
const BREAKEND_FEATURE_DISTANCE: u64 = 5000;

impl Feature {
    /// Whether a breakend at `pos` on the feature's chromosome lies within
    /// [`BREAKEND_FEATURE_DISTANCE`] of it (`_close_to_feature`,
    /// StructuralVariationOverlap.pm:144-160).
    fn close_to(&self, pos: u64) -> bool {
        self.start.saturating_sub(BREAKEND_FEATURE_DISTANCE) <= pos
            && pos <= self.end + BREAKEND_FEATURE_DISTANCE
    }
}

/// The mate of a paired breakend's bracket allele: its chromosome as written and
/// its position. `None` for the record's own breakend, a single breakend and
/// every other variant.
fn breakend_mate(variant: &InputVariant) -> Option<(&str, u64)> {
    if variant.variant_class != VariantClass::Translocation || variant.is_single_breakend {
        return None;
    }
    Some((variant.mate_chr.as_deref()?, variant.mate_pos?))
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{resolved_cache, source, variant, write_fixture_cache};
    use super::*;
    use vep_core::consequence::FeatureType;
    use vep_core::coordinate::Strand;

    /// Ensembl VEP reads a motif allele through `feature_seq`, complemented
    /// when the record's strand differs from the motif's
    /// (`VariationFeatureOverlapAllele.pm:253-256`), so a record given on the
    /// reverse strand scores the genomic change it names: its motif fields are
    /// its forward twin's, on a forward motif and on a reverse one.
    #[test]
    fn reverse_strand_record_scores_its_forward_strand_alleles() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let fields = |v: &InputVariant, motif_id: &str| {
            let row = v
                .transcript_consequences
                .iter()
                .find(|r| r.transcript_id() == motif_id)
                .unwrap();
            let motif = row.motif.as_deref().unwrap();
            (motif.pos, motif.high_inf_pos, motif.score_change)
        };
        // 106 is the second base of the forward motif (`C` of `ACGT`); 116 the
        // second forward-strand base of the reverse motif (`G`, the `C` at
        // position 3 of its own-strand `TGCA`), changed to the base its matrix
        // favours there.
        for (pos, motif_id, (forward_ref, forward_alt), (reverse_ref, reverse_alt)) in [
            (106, "ENSM_FORWARD", (b"C", b"G"), (b"G", b"C")),
            (116, "ENSM_REVERSE", (b"G", b"C"), (b"C", b"G")),
        ] {
            let mut forward = variant(pos, pos, forward_ref, forward_alt);
            source.annotate(&mut forward);
            let expected = fields(&forward, motif_id);
            assert!(expected.2.is_some_and(|delta| delta != 0.0), "{expected:?}");
            let mut reverse = variant(pos, pos, reverse_ref, reverse_alt);
            reverse.strand = Strand::Reverse;
            source.annotate(&mut reverse);
            assert_eq!(fields(&reverse, motif_id), expected, "{motif_id}");
        }
    }

    #[test]
    fn insertion_boundaries_and_reverse_positions() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        for (start, end, expected) in [(105, 104, false), (106, 105, true), (109, 108, false)] {
            let mut v = variant(start, end, b"-", b"A");
            source.annotate(&mut v);
            assert_eq!(
                v.transcript_consequences
                    .iter()
                    .any(|r| r.transcript_id() == "ENSM_FORWARD"),
                expected
            );
        }
        let mut reverse = variant(115, 116, b"AC", b"-");
        source.annotate(&mut reverse);
        let row = reverse
            .transcript_consequences
            .iter()
            .find(|r| r.transcript_id() == "ENSM_REVERSE")
            .unwrap();
        assert_eq!(row.motif.as_deref().unwrap().pos, Some(4)); // variant start, not end
        let mut before = variant(114, 116, b"ACC", b"-");
        source.annotate(&mut before);
        let row = before
            .transcript_consequences
            .iter()
            .find(|r| r.transcript_id() == "ENSM_REVERSE")
            .unwrap();
        assert_eq!(row.motif.as_deref().unwrap().pos, None);
    }

    #[test]
    fn complete_deletion_and_duplication_keep_generic_overlap_term() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        for (class, term) in [
            (
                VariantClass::StructuralDeletion,
                Consequence::RegulatoryRegionAblation,
            ),
            (
                VariantClass::Duplication,
                Consequence::RegulatoryRegionAmplification,
            ),
        ] {
            let mut v = variant(100, 120, b"A", b"<DEL>");
            v.is_structural = true;
            v.variant_class = class;
            source.annotate(&mut v);
            let row = v
                .transcript_consequences
                .iter()
                .find(|r| r.feature_type == FeatureType::RegulatoryFeature)
                .unwrap();
            assert_eq!(
                &row.consequences[..],
                &[term, Consequence::RegulatoryRegionVariant]
            );
            assert!(row.biotype().is_none());
            for row in &v.transcript_consequences {
                assert!(row.motif.is_none(), "structural rows omit motif extras");
            }
        }
    }

    /// A sequence variant whose alternate allele repeats its reference allele
    /// is Ensembl VEP's `tandem_repeat` (`Utils/VariationEffect.pm:282-298`),
    /// a copy-number gain, so covering a feature with it is an amplification
    /// (`feature_amplification`, `:330-335`) beside the overlap term; a longer
    /// allele that is not a repeat is the overlap term alone.
    #[test]
    fn repeated_reference_allele_covering_a_feature_amplifies_it() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let reference = b"ACGTACGTACGTACGTACGTA";
        let repeated = [&reference[..], &reference[..]].concat();
        let mut v = variant(100, 120, reference, &repeated);
        source.annotate(&mut v);
        // The writer orders the rows; the annotator's own order is not part of
        // the contract, so compare them sorted by id.
        let mut terms: Vec<_> = v
            .transcript_consequences
            .iter()
            .map(|r| (r.transcript_id(), r.consequences.to_vec()))
            .collect();
        terms.sort();
        assert_eq!(
            terms,
            [
                (
                    "ENSM_FORWARD",
                    vec![
                        Consequence::TfbsAmplification,
                        Consequence::TfBindingSiteVariant
                    ]
                ),
                (
                    "ENSM_REVERSE",
                    vec![
                        Consequence::TfbsAmplification,
                        Consequence::TfBindingSiteVariant
                    ]
                ),
                (
                    "ENSR_FIXTURE",
                    vec![
                        Consequence::RegulatoryRegionAmplification,
                        Consequence::RegulatoryRegionVariant
                    ]
                ),
            ]
        );
        let mut longer = variant(100, 120, reference, b"ACGTACGTACGTACGTACGTAC");
        source.annotate(&mut longer);
        assert!(longer
            .transcript_consequences
            .iter()
            .all(|r| r.consequences.len() == 1));
        let gain = |reference: &[u8], alternate: &[u8]| {
            vep_effects::sv::perl_gain_class(&variant(100, 101, reference, alternate))
        };
        assert!(gain(b"AC", b"ACACAC"));
        assert!(!gain(b"AC", b"ACAC-"));
        assert!(!gain(b"-", b"ACAC"));
        assert!(!gain(b"AC", b"AC"));
        assert!(!gain(b"AC", b"ACA"));
    }

    /// Ensembl VEP stores a variant's motif overlaps under the motif's database
    /// id and writes them in the string order of those ids
    /// (`VariationFeature.pm:937`, `:750`): a motif with the id `1000000`
    /// precedes one with `999999`, whatever their stable ids sort as. Each
    /// motif row carries its id's rank among the chromosome's motifs, which
    /// the writers order motif rows by, and the id itself, which keys the row
    /// behind the rejoin prefix of a `--minimal` split allele
    /// (`OutputFactory.pm:2325`): a prefix opening with a letter sorts after
    /// the first allele's digit-led keys, one opening with `-` before them. A
    /// regulatory row carries neither.
    #[test]
    fn motif_rows_follow_ensembl_vep_dbid_order() {
        let dir = tempfile::tempdir().unwrap();
        let motif = |stable_id: &str, dbid: u64| {
            format!(
                r#"[{{"stable_id":"{stable_id}","dbid":{dbid},"start":200,"end":203,"strand":1,"feature_type":"MotifFeature","score":1.0,"seq":"ACGT","matrix":"ENSPFM_FORWARD"}}]"#
            )
        };
        let early = motif("ENSM00000999999", 999_999);
        let late = motif("ENSM00001000000", 1_000_000);
        let same = motif("ENSM_SAME_ID", 1_000_000);
        write_fixture_cache(
            dir.path(),
            &[
                ("2000-2100.json", &early),
                ("3000-3100.json", &late),
                ("4000-4100.json", &same),
            ],
        );
        let source = RegulatorySource::open(&resolved_cache(dir.path()), true)
            .unwrap()
            .unwrap();
        source.prewarm("21").unwrap();
        let mut v = variant(201, 201, b"C", b"G");
        source.annotate(&mut v);
        let ranks: Vec<(&str, u32)> = v
            .transcript_consequences
            .iter()
            .map(|row| (row.transcript_id(), row.id_rank))
            .collect();
        assert!(ranks.contains(&("ENSM00000999999", 4)), "{ranks:?}");
        assert!(ranks.contains(&("ENSM00001000000", 1)), "{ranks:?}");
        assert!(ranks.contains(&("ENSM_SAME_ID", 1)), "{ranks:?}");
        let db_ids: Vec<(&str, Option<u64>)> = v
            .transcript_consequences
            .iter()
            .map(|row| (row.transcript_id(), row.motif_db_id))
            .collect();
        assert!(
            db_ids.contains(&("ENSM00000999999", Some(999_999))),
            "{db_ids:?}"
        );
        assert!(
            db_ids.contains(&("ENSM00001000000", Some(1_000_000))),
            "{db_ids:?}"
        );
        assert!(
            db_ids.contains(&("ENSM_SAME_ID", Some(1_000_000))),
            "{db_ids:?}"
        );
        // Two motifs under one id cannot come from a database, so their order
        // between themselves is not pinned; both precede the lower rank.
        let written: Vec<&str> = vep_io::output::fields::record_rows(&[&v])
            .iter()
            .filter_map(|(_, row)| row.map(|row| row.transcript_id()))
            .collect();
        assert_eq!(written.len(), 3);
        assert_eq!(written[2], "ENSM00000999999", "{written:?}");
        assert!(written[..2].contains(&"ENSM00001000000"), "{written:?}");
        for (prefix, expected_alleles) in
            [("A/G_", [0, 0, 0, 1, 1, 1]), ("-/G_", [1, 1, 1, 0, 0, 0])]
        {
            let mut second = variant(201, 201, b"A", b"G");
            second.allele_index = 1;
            second.rejoin_key_prefix = Some(prefix.into());
            source.annotate(&mut second);
            let written_alleles: Vec<usize> = vep_io::output::fields::record_rows(&[&v, &second])
                .iter()
                .filter(|(_, row)| row.is_some())
                .map(|(allele, _)| allele.allele_index)
                .collect();
            assert_eq!(written_alleles, expected_alleles, "{prefix}");
        }
        let mut fixture = variant(116, 116, b"G", b"A");
        source.annotate(&mut fixture);
        let regulatory = fixture
            .transcript_consequences
            .iter()
            .find(|row| row.feature_type == FeatureType::RegulatoryFeature)
            .unwrap();
        assert_eq!(regulatory.id_rank, 0);
        assert_eq!(regulatory.motif_db_id, None);
    }

    /// Under `--cell_type` a regulatory row's `CELL_TYPE` list is the requested
    /// epigenomes projected onto the feature's activity string, in the order
    /// requested, each with its activity decoded (A ACTIVE, I INACTIVE, P
    /// POISED, R REPRESSED, N NA) and a name the feature's map lacks (`-`)
    /// skipped, as `get_cell_types` keeps the names the map has a value for
    /// (OutputFactory.pm:1968-1978); a motif row's list is empty, since a motif
    /// carries no map (`Database/RegFeat.pm:212-221`); without the flag no row
    /// carries one. The names come from the cache's list through the greedy
    /// split, so `A549_(m,_58_y)` is one epigenome.
    #[test]
    fn cell_type_projects_the_requested_epigenomes_onto_the_activity() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let requested = CellTypes::resolve(
            Some(
                "prostate_gland_(m,_37_y),astrocyte,psoas_muscle_(m,_3_y),A549_(m,_58_y),\
                 placenta_(85-113_dpf),neuronal_stem_cell_(m),astrocyte",
            ),
            Some(&source),
        )
        .unwrap()
        .unwrap();
        let mut v = variant(116, 116, b"G", b"A");
        source.annotate_rows(&mut v, Some(&requested));
        let mut rows: Vec<(&str, Vec<String>)> = v
            .transcript_consequences
            .iter()
            .map(|r| {
                (
                    r.transcript_id(),
                    r.cell_type
                        .as_deref()
                        .unwrap()
                        .iter()
                        .map(ToString::to_string)
                        .collect(),
                )
            })
            .collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        assert_eq!(
            rows,
            [
                ("ENSM_REVERSE", vec![]),
                (
                    "ENSR_FIXTURE",
                    [
                        "prostate_gland_(m,_37_y):NA",
                        "astrocyte:INACTIVE",
                        "A549_(m,_58_y):ACTIVE",
                        "placenta_(85-113_dpf):REPRESSED",
                        "neuronal_stem_cell_(m):POISED",
                        "astrocyte:INACTIVE",
                    ]
                    .map(str::to_string)
                    .to_vec()
                ),
            ]
        );
        // A structural variant's regulatory and motif rows carry the list too
        // (OutputFactory.pm:2107-2108).
        let mut sv = variant(100, 120, b"A", b"<DEL>");
        sv.is_structural = true;
        sv.variant_class = VariantClass::StructuralDeletion;
        source.annotate_rows(&mut sv, Some(&requested));
        for row in &sv.transcript_consequences {
            let list = row.cell_type.as_deref().unwrap();
            assert_eq!(
                list.is_empty(),
                row.feature_type == FeatureType::MotifFeature,
                "{}",
                row.transcript_id()
            );
        }
        let mut plain = variant(116, 116, b"G", b"A");
        source.annotate(&mut plain);
        assert!(plain
            .transcript_consequences
            .iter()
            .all(|r| r.cell_type.is_none()));

        // The check runs against the source's list: an unknown name is Ensembl
        // VEP's refusal, and without a source (a cache with no regulatory data)
        // nothing is checked, as Ensembl VEP builds no regulatory source there.
        let error = CellTypes::resolve(Some("astrocyte,HUVEC"), Some(&source))
            .unwrap_err()
            .to_string();
        assert!(
            error.starts_with(
                "ERROR: [cli] Cell type \"HUVEC\" unavailable; available cell types are:\n\
                 A549_(m,_58_y),A673_(f,_15_y),"
            ),
            "{error}"
        );
        assert!(CellTypes::resolve(Some("HUVEC"), None).unwrap().is_none());
        assert!(CellTypes::resolve(None, Some(&source)).unwrap().is_none());
    }

    /// A regulatory row carries the feature's SO term as its biotype and strand
    /// 0; a motif row carries the motif's strand, its matrix id as the motif
    /// name and the matrix's factors, and no biotype.
    #[test]
    fn rows_carry_the_feature_kind_fields() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let mut v = variant(116, 116, b"G", b"A");
        source.annotate(&mut v);
        let mut rows: Vec<_> = v
            .transcript_consequences
            .iter()
            .map(|r| {
                (
                    r.transcript_id(),
                    r.feature_type,
                    r.strand,
                    r.biotype(),
                    r.motif
                        .as_deref()
                        .map(|m| (m.name.as_ref(), m.pos, m.factors.len())),
                )
            })
            .collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        assert_eq!(
            rows,
            [
                (
                    "ENSM_REVERSE",
                    FeatureType::MotifFeature,
                    -1,
                    None,
                    Some(("ENSPFM_REVERSE", Some(3), 0))
                ),
                (
                    "ENSR_FIXTURE",
                    FeatureType::RegulatoryFeature,
                    0,
                    Some("enhancer"),
                    None
                ),
            ]
        );
    }

    /// Each breakend of a paired record is a query of its own on its chromosome
    /// (InputBuffer.pm:335-342 tests every breakend position; `_close_to_feature`,
    /// StructuralVariationOverlap.pm:144-160, matches its chromosome): the bracket
    /// allele whose mate lies in a feature gets that feature's row with the overlap
    /// term alone, a breakpoint being neither a deletion nor a gain; its own
    /// position's features are kept only when the mate is within 5000 bases of
    /// them on the same chromosome, and never for a mate on another chromosome; the
    /// record's own breakend keeps its own position's rows and gains the rows of
    /// the features a same-chromosome mate lies in when they are within 5000
    /// bases of its own position.
    #[test]
    fn breakend_positions_overlap_regulatory_features() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let bracket = |pos: u64, mate_chr: &str, mate_pos: u64| {
            let mut v = variant(
                pos,
                pos,
                b"N",
                format!("N[{mate_chr}:{mate_pos}[").as_bytes(),
            );
            v.is_structural = true;
            v.variant_class = VariantClass::Translocation;
            v.sv_end = Some(pos);
            v.mate_id = Some("mate".into());
            v.mate_chr = Some(mate_chr.into());
            v.mate_pos = Some(mate_pos);
            v
        };
        let ids = |v: &InputVariant| -> Vec<String> {
            let mut ids: Vec<String> = v
                .transcript_consequences
                .iter()
                .map(|r| format!("{}:{}", r.transcript_id(), r.consequences[0].so_term()))
                .collect();
            ids.sort();
            ids
        };

        // Own position outside every feature, mate inside the regulatory feature and
        // the forward motif.
        let mut mate_inside = bracket(50_000, "21", 106);
        source.annotate(&mut mate_inside);
        assert_eq!(
            ids(&mate_inside),
            [
                "ENSM_FORWARD:TF_binding_site_variant",
                "ENSR_FIXTURE:regulatory_region_variant"
            ]
        );

        // Own position inside the regulatory feature, mate far away on the same
        // chromosome: the bracket allele is not close to the feature.
        let mut mate_far = bracket(110, "21", 50_000);
        source.annotate(&mut mate_far);
        assert!(mate_far.transcript_consequences.is_empty());

        // Own position inside, mate within 5000 bases: the own position's rows stay,
        // the feature counted once when both positions lie in it.
        let mut mate_near = bracket(110, "21", 118);
        source.annotate(&mut mate_near);
        assert_eq!(
            ids(&mate_near),
            [
                "ENSM_REVERSE:TF_binding_site_variant",
                "ENSR_FIXTURE:regulatory_region_variant"
            ]
        );

        // A mate on another chromosome: no row for the own position, and the mate
        // chromosome's features where that chromosome is loaded.
        let mut other_chromosome = bracket(110, "22", 110);
        source.annotate(&mut other_chromosome);
        assert!(other_chromosome.transcript_consequences.is_empty());

        // The record's own breakend queries its own position, and the position of
        // each mate on its chromosome for the features within 5000 bases of it.
        let own = |pos: u64, mates: &[(&str, u64)]| {
            let mut v = variant(pos, pos, b"N", b"N.");
            v.is_structural = true;
            v.variant_class = VariantClass::Translocation;
            v.is_single_breakend = true;
            v.sv_end = Some(pos);
            v.breakend_mates = mates
                .iter()
                .map(|(chr, pos)| (chr.to_string(), *pos))
                .collect();
            v
        };
        let mut inside = own(110, &[("21", 50_000)]);
        source.annotate(&mut inside);
        assert_eq!(ids(&inside), ["ENSR_FIXTURE:regulatory_region_variant"]);
        let mut near_mate_inside = own(2_000, &[("21", 106)]);
        source.annotate(&mut near_mate_inside);
        assert_eq!(
            ids(&near_mate_inside),
            [
                "ENSM_FORWARD:intergenic_variant",
                "ENSR_FIXTURE:intergenic_variant"
            ]
        );
        assert_eq!(
            near_mate_inside.most_severe_consequence,
            Some(Consequence::IntergenicVariant)
        );
        let mut far_mate_inside = own(50_000, &[("21", 106)]);
        source.annotate(&mut far_mate_inside);
        assert!(far_mate_inside.transcript_consequences.is_empty());
        let mut other_chromosome_mate = own(2_000, &[("22", 106)]);
        source.annotate(&mut other_chromosome_mate);
        assert!(other_chromosome_mate.transcript_consequences.is_empty());
    }

    /// A mate written under another spelling of the record's chromosome (`chr21`
    /// against a cache that names it `21`) is the same chromosome once the run
    /// has named the record's chromosomes (`pipeline::name_source_chromosomes`;
    /// Ensembl VEP resolves a breakend's name through `get_source_chr_name` as
    /// the features at its position are fetched, AnnotationSource.pm:251-254,
    /// and tests the own allele's closeness on the record's resolved name,
    /// StructuralVariationOverlap.pm:144-152), so the record's own breakend
    /// gains the rows of the features that mate lies in under either spelling.
    #[test]
    fn a_mate_spelled_otherwise_than_the_cache_reaches_the_own_breakend() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let transcripts = crate::transcript_index::LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let synonyms = vep_core::synonyms::ChromosomeSynonyms::new();
        let own_rows = |line: &str| -> Vec<String> {
            let mut alleles = crate::vcf_parser::parse_vcf_line(line, false).unwrap();
            let own = &mut alleles[0];
            crate::pipeline::name_source_chromosomes(
                own,
                &synonyms,
                &transcripts,
                crate::validate::NO_CHECKS.source_chromosomes,
                None,
            );
            assert!(own.is_single_breakend && own.mate_chr.is_none());
            source.annotate(own);
            let mut ids: Vec<String> = own
                .transcript_consequences
                .iter()
                .map(|r| r.transcript_id().to_string())
                .collect();
            ids.sort();
            ids
        };
        let expected = ["ENSM_FORWARD", "ENSR_FIXTURE"];
        assert_eq!(
            own_rows("21\t2000\tbare\tN\tN[21:106[\t.\t.\tSVTYPE=BND"),
            expected
        );
        assert_eq!(
            own_rows("chr21\t2000\tprefixed\tN\tN[chr21:106[\t.\t.\tSVTYPE=BND"),
            expected
        );
    }
}
