// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SPDI notation (`SPDI`) for each allele of a record, Ensembl VEP's
//! `spdi_genomic` on the variation feature as the parser leaves it, and the
//! GA4GH VRS allele (`GA4GH_VRS`) built from it.
//!
//! Perl citations name modules of ensembl-vep release/116.2
//! (`Bio/EnsEMBL/VEP/...`) unless another repository is named.

use std::collections::HashMap;

use vep_core::consequence::FeatureType;
use vep_core::coordinate::Strand;
use vep_core::spdi::{ga4gh_vrs_from_spdi, spdi_genomic, SiblingSpdi};
use vep_core::variant::InputVariant;

use crate::consequences::EffectsConfig;

/// What `--spdi` and `--ga4gh_vrs` compute.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SpdiOptions {
    /// Build `InputVariant::ga4gh_vrs` (`--ga4gh_vrs`) from the SPDI named by
    /// the RefSeq accession each cache chromosome name maps to; a chromosome
    /// absent from the map gets no object, as Ensembl VEP writes none for an
    /// SPDI that names no `NC_` sequence (OutputFactory.pm 1403). `None` builds
    /// no object.
    pub vrs_accessions: Option<HashMap<String, String>>,
}

/// Sets `variant.spdi`, and `variant.ga4gh_vrs` when asked, for a sequence
/// variant, and `variant.sibling_spdi` for the other alleles of a multi-allelic
/// record; a structural variant never carries them (Ensembl VEP builds its rows
/// outside `VariationFeatureOverlapAllele_to_output_hash`). Runs after the
/// transcript consequences are known, because they decide the chromosome's
/// name (`reference_name`).
pub fn annotate(variant: &mut InputVariant, config: &EffectsConfig) {
    let Some(options) = &config.spdi else {
        return;
    };
    // Ensembl VEP builds a notation per allele after the first of the feature's
    // string and looks the row's allele up among them (OutputFactory.pm:1396-1400;
    // VariationFeature.pm:2153 drops the first), so the reference row of
    // `--process_ref_homs` carries none, except on a reference-homozygous
    // individual's copy, whose string is the reference twice (Parser/VCF.pm:646)
    // and whose one notation is the reference against itself.
    if variant.is_structural || (variant.reference_allele_row && !variant.hom_ref) {
        return;
    }
    let spdi = |name: &str, alt: &[u8]| {
        spdi_genomic(
            name,
            variant.start,
            variant.end,
            &variant.ref_allele,
            alt,
            variant.strand == Strand::Reverse,
        )
    };
    let chromosome_name = reference_name(variant, config);
    // With a database Ensembl VEP names the SPDI by the slice's RefSeq synonym
    // (ensembl-variation VariationFeature.pm 2150) and the object follows from
    // that name; offline the accession comes from the synonyms table instead.
    let accession = options
        .vrs_accessions
        .as_ref()
        .and_then(|accessions| accessions.get(&variant.chr));
    let vrs = |alt: &[u8]| {
        accession
            .and_then(|accession| spdi(accession, alt))
            .and_then(|named| ga4gh_vrs_from_spdi(&named))
    };
    let own_allele = variant.alt_allele();
    let own_spdi = spdi(&chromosome_name, own_allele);
    let own_vrs = vrs(own_allele);
    // Ensembl VEP writes every alternate allele of the record's allele string
    // into one table read by the row's Allele (OutputFactory.pm 1396-1400); the
    // siblings share the record's span and reference allele.
    let record_alt_alleles = variant
        .record_allele_string_multi
        .iter()
        .flat_map(|record| record.split('/').skip(1));
    let siblings: Vec<SiblingSpdi> = record_alt_alleles
        .filter(|allele| allele.as_bytes() != own_allele)
        .filter_map(|allele| {
            Some(SiblingSpdi {
                allele: allele.to_string(),
                spdi: spdi(&chromosome_name, allele.as_bytes())?,
                ga4gh_vrs: vrs(allele.as_bytes()),
            })
        })
        .collect();
    variant.spdi = own_spdi;
    variant.ga4gh_vrs = own_vrs;
    variant.sibling_spdi = siblings;
}

/// The chromosome name an SPDI carries offline, `seq_region_name` of the slice
/// the variation feature ends up on (ensembl-variation VariationFeature.pm
/// 2150-2151): a variant that overlaps a transcript takes the transcript's
/// slice, the cache's name for the chromosome (AnnotationType/Transcript.pm
/// 122); one that overlaps none is given a slice built from the reference FASTA
/// under the FASTA's own name for the input chromosome (InputBuffer.pm 518,
/// BaseVEP.pm 409-437). Without a FASTA Ensembl VEP has no slice for such a
/// variant and stops; the cache's name is the one value left to print. A
/// regulatory or motif row is not a transcript, so an allele with such rows
/// alone takes the FASTA's name as an intergenic one does.
fn reference_name(variant: &InputVariant, config: &EffectsConfig) -> String {
    let overlaps_a_transcript = variant
        .transcript_consequences
        .iter()
        .any(|tc| tc.feature_type == FeatureType::Transcript);
    if !overlaps_a_transcript {
        if let Some(name) = config
            .reference_fasta
            .as_deref()
            .and_then(|fasta| fasta.record_name(&variant.original_chr))
        {
            return name.to_string();
        }
    }
    variant.chr.clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use vep_core::consequence::TranscriptConsequence;
    use vep_core::spdi::VrsAllele;

    /// A reference whose one record is named `chr21`, 120 bases long.
    fn chr21_fasta() -> (tempfile::TempDir, Arc<vep_fasta::IndexedFasta>) {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let seq = "ACGT".repeat(30);
        std::fs::write(&fa, format!(">chr21\n{seq}\n")).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "chr21\t120\t7\t120\t121\n").unwrap();
        (
            dir,
            Arc::new(vep_fasta::IndexedFasta::from_path(&fa).unwrap()),
        )
    }

    fn variant(chr: &str, start: u64, end: u64, reference: &str, alt: &str) -> InputVariant {
        InputVariant::new(
            chr.to_string(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn with_transcript_row(mut variant: InputVariant) -> InputVariant {
        variant
            .transcript_consequences
            .push(TranscriptConsequence::default());
        variant
    }

    fn spdi_only() -> EffectsConfig {
        EffectsConfig {
            spdi: Some(SpdiOptions::default()),
            ..Default::default()
        }
    }

    fn spdi_of(variant: &mut InputVariant, config: &EffectsConfig) -> Option<String> {
        annotate(variant, config);
        variant.spdi.clone()
    }

    /// Ensembl VEP 116.2 under `--spdi --offline` names the chromosome as the
    /// cache does on every transcript row (`21` for an input `chr21`), and as
    /// the FASTA does on an intergenic row (`chr21` when the FASTA's record is
    /// `chr21`, whatever the input wrote); an intergenic row without a FASTA,
    /// where Ensembl VEP stops, takes the cache's name.
    #[test]
    fn chromosome_is_named_by_the_cache_on_transcript_rows_and_by_the_fasta_elsewhere() {
        let (_dir, fasta) = chr21_fasta();
        let with_fasta = EffectsConfig {
            reference_fasta: Some(Arc::clone(&fasta)),
            ..spdi_only()
        };
        let mut genic = with_transcript_row(variant("chr21", 50, 50, "A", "C"));
        assert_eq!(
            spdi_of(&mut genic, &with_fasta).as_deref(),
            Some("21:49:A:C")
        );
        let mut intergenic = variant("21", 50, 50, "A", "C");
        assert_eq!(
            spdi_of(&mut intergenic, &with_fasta).as_deref(),
            Some("chr21:49:A:C")
        );
        let mut regulatory_only = variant("21", 50, 50, "A", "C");
        regulatory_only
            .transcript_consequences
            .push(TranscriptConsequence {
                feature_type: FeatureType::RegulatoryFeature,
                ..TranscriptConsequence::default()
            });
        assert_eq!(
            spdi_of(&mut regulatory_only, &with_fasta).as_deref(),
            Some("chr21:49:A:C")
        );
        let mut intergenic = variant("chr21", 50, 50, "A", "C");
        assert_eq!(
            spdi_of(&mut intergenic, &spdi_only()).as_deref(),
            Some("21:49:A:C")
        );
        let mut off = variant("21", 50, 50, "A", "C");
        assert_eq!(spdi_of(&mut off, &EffectsConfig::default()), None);
    }

    /// The alleles are the record's as the parser leaves them: a bi-allelic
    /// indel's minimised pair at its own span (`21 10561150 AT A` is `T/-` at
    /// 10561151), a multi-allelic record's anchor-chopped pair at the shared
    /// span (`21 14144642 AT CT,A` keeps `AT/CT` and `AT/A` at 14144642-14144643,
    /// as Ensembl VEP prints `21:14144641:AT:A`); an allele given on the reverse
    /// strand is complemented; a structural variant and a `*` allele carry none.
    #[test]
    fn alleles_follow_the_record_as_parsed() {
        let config = spdi_only();
        let mut deletion = variant("21", 10561151, 10561151, "T", "-");
        assert_eq!(
            spdi_of(&mut deletion, &config).as_deref(),
            Some("21:10561150:T:")
        );
        let mut shared = variant("21", 14144642, 14144643, "AT", "A");
        assert_eq!(
            spdi_of(&mut shared, &config).as_deref(),
            Some("21:14144641:AT:A")
        );
        let mut reverse = variant("21", 20003758, 20003757, "-", "A");
        reverse.strand = Strand::Reverse;
        assert_eq!(
            spdi_of(&mut reverse, &config).as_deref(),
            Some("21:20003757::T")
        );
        let mut structural = variant("21", 101, 200, "N", "<DEL>");
        structural.is_structural = true;
        assert_eq!(spdi_of(&mut structural, &config), None);
        assert_eq!(structural.ga4gh_vrs, None);
        let mut spanning = variant("21", 100, 100, "A", "*");
        assert_eq!(spdi_of(&mut spanning, &config), None);
    }

    /// Ensembl VEP 116.2 with a database under `--ga4gh_vrs` writes, for
    /// `21 10538885 A C`, the allele of `NC_000021.9:10538884:A:C`; here the
    /// accession is the one the synonyms map gives the cache's chromosome name,
    /// and a chromosome with none gets no object. The SPDI column keeps the
    /// offline name.
    #[test]
    fn ga4gh_vrs_is_built_from_the_spdi_named_by_the_accession() {
        let config = EffectsConfig {
            spdi: Some(SpdiOptions {
                vrs_accessions: Some(HashMap::from([(
                    "21".to_string(),
                    "NC_000021.9".to_string(),
                )])),
            }),
            ..Default::default()
        };
        let mut mapped = with_transcript_row(variant("chr21", 10538885, 10538885, "A", "C"));
        annotate(&mut mapped, &config);
        assert_eq!(mapped.spdi.as_deref(), Some("21:10538884:A:C"));
        assert_eq!(
            mapped.ga4gh_vrs,
            Some(VrsAllele {
                sequence_id: "NC_000021.9".to_string(),
                start: 10538884,
                end: 10538885,
                sequence: "C".to_string(),
            })
        );
        let mut unmapped = variant("X", 100, 100, "A", "C");
        annotate(&mut unmapped, &config);
        assert_eq!(unmapped.spdi.as_deref(), Some("X:99:A:C"));
        assert_eq!(unmapped.ga4gh_vrs, None);
        let mut plain = variant("21", 10538885, 10538885, "A", "C");
        annotate(&mut plain, &spdi_only());
        assert_eq!(plain.ga4gh_vrs, None);
    }

    /// Ensembl VEP keys the record's notations by every allele of its allele
    /// string (`-/TCA/CAT/*` for `21 13620283 G GTCA,GCAT,*`) and reads them by
    /// the row's Allele, so the `TCA` allele carries the `CAT` allele's SPDI and
    /// VRS allele for a row a shift rotated onto it, on the record's own span;
    /// the `*` allele has no notation, and a bi-allelic record has no siblings.
    #[test]
    fn siblings_carry_the_other_alleles_notations_on_the_records_span() {
        let config = EffectsConfig {
            spdi: Some(SpdiOptions {
                vrs_accessions: Some(HashMap::from([(
                    "21".to_string(),
                    "NC_000021.9".to_string(),
                )])),
            }),
            ..Default::default()
        };
        let mut tri_allelic = variant("21", 13620284, 13620283, "-", "TCA");
        tri_allelic.record_allele_string_multi = Some(Arc::from("-/TCA/CAT/*"));
        annotate(&mut tri_allelic, &config);
        assert_eq!(tri_allelic.spdi.as_deref(), Some("21:13620283::TCA"));
        assert_eq!(
            tri_allelic.sibling_spdi,
            vec![SiblingSpdi {
                allele: "CAT".to_string(),
                spdi: "21:13620283::CAT".to_string(),
                ga4gh_vrs: Some(VrsAllele {
                    sequence_id: "NC_000021.9".to_string(),
                    start: 13620283,
                    end: 13620283,
                    sequence: "CAT".to_string(),
                }),
            }]
        );
        let mut bi_allelic = variant("21", 13620284, 13620283, "-", "TCA");
        annotate(&mut bi_allelic, &config);
        assert_eq!(bi_allelic.spdi.as_deref(), Some("21:13620283::TCA"));
        assert!(bi_allelic.sibling_spdi.is_empty());
    }

    /// Ensembl VEP keys the notations by the alleles after the first of the
    /// feature's string (VariationFeature.pm:2153-2167), so the reference row of
    /// `--process_ref_homs` finds none on a record or an individual's copy whose
    /// string names an alternate, and on a reference-homozygous individual's copy,
    /// whose string is the reference twice (Parser/VCF.pm:646), finds the
    /// reference described against itself, in each class's own shape (:2185-2220):
    /// `21:99:A:A` for a base, `21:100:AC:AC` for a deletion's bases, and
    /// `21:30284295::-` for an insertion's empty reference.
    #[test]
    fn reference_row_carries_a_notation_only_on_a_reference_homozygous_copy() {
        let config = spdi_only();
        let mut plain = variant("21", 100, 100, "A", "G").reference_allele_row();
        assert_eq!(spdi_of(&mut plain, &config), None);
        let mut het = variant("21", 100, 100, "A", "G");
        het.individual = Some("S1".into());
        het.genotype = vec!["A".into(), "G".into()];
        let mut het_reference = het.reference_allele_row();
        assert_eq!(spdi_of(&mut het_reference, &config), None);
        assert!(het_reference.sibling_spdi.is_empty());

        let hom_ref_copy = |start: u64, end: u64, reference: &str| {
            let mut copy = variant("21", start, end, reference, "*");
            copy.record_allele_string_multi = Some(Arc::from(format!("{reference}/{reference}")));
            copy.individual = Some("S2".into());
            copy.genotype = vec![reference.into(), reference.into()];
            copy.hom_ref = true;
            copy.reference_allele_row()
        };
        let mut snv = hom_ref_copy(100, 100, "A");
        assert_eq!(spdi_of(&mut snv, &config).as_deref(), Some("21:99:A:A"));
        assert!(snv.sibling_spdi.is_empty());
        let mut deletion = hom_ref_copy(101, 102, "AC");
        assert_eq!(
            spdi_of(&mut deletion, &config).as_deref(),
            Some("21:100:AC:AC")
        );
        let mut insertion = hom_ref_copy(30284296, 30284295, "-");
        assert_eq!(
            spdi_of(&mut insertion, &config).as_deref(),
            Some("21:30284295::-")
        );
    }
}
