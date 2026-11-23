// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` source over a tabix-indexed BED, read as
//! `AnnotationSource/File/BED.pm` reads one through ensembl-io's `Parser/Bed.pm`:
//! a record spans `chromStart + 1..=chromEnd`, its fourth column names it and
//! its fifth is its score (`Bed.pm:159-257`); every type matches through the
//! base predicates (`File.pm:534-575`).

use std::sync::Arc;

use indexmap::IndexMap;
use vep_core::custom_annotation::{CustomRecord, CustomValue};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::InputVariant;

use super::grammar::CustomSpec;
use super::source::{
    base_overlap, CustomError, CustomSource, FormatRules, RecordGroup, TabixSource, Verdict,
};
use crate::tabix::TabixRecord;
use crate::traits::PrefetchData;

/// A `--custom` source of `format=bed`.
pub struct BedSource {
    base: TabixSource,
}

impl BedSource {
    pub fn open(
        spec: CustomSpec,
        synonyms: Arc<ChromosomeSynonyms>,
        assembly: Option<vep_core::assembly::Assembly>,
    ) -> Result<Self, CustomError> {
        Ok(BedSource {
            base: TabixSource::open(spec, synonyms, assembly)?,
        })
    }
}

impl CustomSource for BedSource {
    fn spec(&self) -> &CustomSpec {
        self.base.spec()
    }

    fn chromosomes(&self) -> Vec<String> {
        self.base.chromosomes()
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<PrefetchData, CustomError> {
        self.base.prefetch(variants)
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), CustomError> {
        self.base.annotate(variants, data, self)
    }
}

impl FormatRules for BedSource {
    /// The base `_record_overlaps_VF` on `chromStart + 1..=chromEnd`
    /// (`File.pm:523, 534-575`), the walk ending at a record starting past the
    /// feature's end plus one (`:333-335`, `get_start` being `chromStart + 1`,
    /// `Bed.pm:159-162`), and the base `_create_records` (`:460-497`): the
    /// record is named by its fourth column, or by its coordinates under
    /// `coords` or when it has no fourth column, the chromosome without a
    /// leading `chr` (`Bed.pm:118-122`); its score is the fifth column, none
    /// when absent or `.` (`:249-257`). `chromStart` and `chromEnd` are the
    /// numbers Perl reads from them ([`TabixSource::start_of`],
    /// [`TabixSource::end_of`]).
    fn match_record(&self, record: &TabixRecord, group: &RecordGroup<'_>, scores: bool) -> Verdict {
        let (start, end) = (self.base.start_of(record), self.base.end_of(record));
        if start.is_past(group.end as i64 + 1) {
            return Verdict::PastEnd;
        }
        let Some(percentage) = base_overlap(
            self.base.spec(),
            start,
            end,
            group.start as i64,
            group.end as i64,
        ) else {
            return Verdict::NoMatch;
        };
        let name = match record.columns.get(3) {
            Some(name) if !self.base.spec().report_coords => name.clone(),
            _ => format!("{}:{}-{}", seqname(&record.chr), start.whole(), end.whole()),
        };
        let score = scores.then(|| match record.columns.get(4) {
            Some(score) if score != "." => CustomValue::Text(score.clone()),
            _ => CustomValue::Undefined,
        });
        Verdict::Match(
            vec![CustomRecord::new(name, None, IndexMap::new(), score)],
            Some(percentage),
        )
    }
}

/// `get_seqname` (`Bed.pm:118-122`): the chromosome without a leading `chr`
/// in any case.
fn seqname(chr: &str) -> &str {
    match chr.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &chr[3..],
        _ => chr,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::custom::grammar::parse_custom;
    use crate::custom::source::CustomSources;
    use crate::custom::CustomRunOptions;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use vep_core::custom_annotation::{CustomAnnotations, CustomFieldValue, CustomStats};

    fn open(argument: &str) -> CustomSources {
        CustomSources::open(
            vec![parse_custom(argument).unwrap()],
            Arc::default(),
            CustomRunOptions::default(),
        )
        .unwrap()
    }

    fn variant(start: u64, end: u64, record: u64) -> InputVariant {
        let mut variant = InputVariant::new("21".into(), start, end, b"A".to_vec(), b"G".to_vec());
        variant.input_record = record;
        variant
    }

    fn only(variant: &InputVariant) -> &CustomAnnotations {
        assert_eq!(variant.custom_annotations.len(), 1, "{variant:?}");
        &variant.custom_annotations[0]
    }

    /// `Bed.pm:159-162, 197-200`: the record `21 99 110 region_1 7` spans
    /// 100-110, so a variant at 99 misses it and one at 100 matches it; the
    /// record is named by its fourth column, or `21:100-110` under `coords`
    /// and for a three-column record, a `chr21` file losing the prefix in that
    /// name alone (`:118-122`); under `summary_stats` the score is the fifth
    /// column as written, none for a record without one (`:249-257`), and a
    /// `chr`-named file answers an input named `21` (`BaseVEP.pm:728-729`).
    #[test]
    fn bed_records_are_zero_based_with_name_and_score() {
        let dir = tempfile::tempdir().unwrap();
        let bed = write_tabix_fixture(
            dir.path(),
            "regions.bed.gz",
            &["21\t99\t110\tregion_1\t7", "21\t200\t210"],
            IndexSpec::standard(),
        );
        let sources = open(&format!(
            "file={},format=bed,summary_stats=min%count",
            bed.display()
        ));
        let mut batch = vec![
            variant(99, 99, 1),
            variant(100, 100, 2),
            variant(205, 205, 3),
        ];
        sources.run_batch(&mut batch).unwrap();
        // A statistic is computed over a record that matched nothing too
        // (`File.pm:344-352`).
        assert!(only(&batch[0]).records.is_empty());
        assert_eq!(
            only(&batch[0]).stats,
            Some(CustomStats::Values(
                [("count".to_string(), CustomValue::Number(0.0))]
                    .into_iter()
                    .collect()
            ))
        );
        let region = &only(&batch[1]).records[0];
        assert_eq!(region.name, "region_1");
        assert_eq!(region.score, Some(CustomValue::Text("7".into())));
        assert_eq!(
            only(&batch[1]).stats,
            Some(CustomStats::Values(
                [
                    ("min".to_string(), CustomValue::Text("7".into())),
                    ("count".to_string(), CustomValue::Number(1.0)),
                ]
                .into_iter()
                .collect()
            ))
        );
        let bare = &only(&batch[2]).records[0];
        assert_eq!(bare.name, "21:201-210");
        assert_eq!(bare.score, Some(CustomValue::Undefined));
        assert_eq!(
            only(&batch[2]).stats,
            Some(CustomStats::Values(
                [("count".to_string(), CustomValue::Number(1.0))]
                    .into_iter()
                    .collect()
            ))
        );

        let coords = open(&format!("file={},format=bed,coords=1", bed.display()));
        let mut batch = vec![variant(100, 100, 1)];
        coords.run_batch(&mut batch).unwrap();
        let region = &only(&batch[0]).records[0];
        assert_eq!(region.name, "21:100-110");
        assert_eq!(region.score, None);

        let prefixed = write_tabix_fixture(
            dir.path(),
            "regions_chr.bed.gz",
            &["chr21\t99\t110"],
            IndexSpec::standard(),
        );
        let sources = open(&format!("file={},format=bed", prefixed.display()));
        let mut batch = vec![variant(100, 100, 1)];
        sources.run_batch(&mut batch).unwrap();
        assert_eq!(only(&batch[0]).records[0].name, "21:100-110");
    }

    /// A record whose `chromStart` is a minus sign and digits (`-5`) is listed
    /// as Ensembl VEP 116.3 lists it: for every record from the first base,
    /// where `tabix` indexes it (htslib 1.9 `tbx.c:99`), to its end and for
    /// none past it, under its fourth column or, without one or under
    /// `coords`, under its coordinates, the start `chromStart + 1` on the
    /// number written (`Bed.pm:159-162`, `File.pm:489-495`). The percentage is
    /// the share of the span from that start (`File.pm:552-563`): `-5` to 10
    /// is 15 bases, where 0 to 10 is 10.
    #[test]
    fn a_record_whose_chrom_start_a_minus_sign_leads_is_listed_from_the_first_base() {
        let dir = tempfile::tempdir().unwrap();
        let listed = |lines: &[&str], options: &str| -> Vec<Vec<(String, Option<String>)>> {
            let bed =
                write_tabix_fixture(dir.path(), "regions.bed.gz", lines, IndexSpec::standard());
            let sources = open(&format!("file={},format=bed{options}", bed.display()));
            let mut batch = vec![
                variant(1, 1, 1),
                variant(10, 10, 2),
                variant(11, 11, 3),
                variant(100, 100, 4),
            ];
            sources.run_batch(&mut batch).unwrap();
            batch
                .iter()
                .map(|variant| {
                    variant
                        .custom_annotations
                        .iter()
                        .flat_map(|annotation| &annotation.records)
                        .map(|record| {
                            let percentage = match record.fields.get("PC") {
                                Some(CustomFieldValue::Text(percentage)) => {
                                    Some(percentage.clone())
                                }
                                _ => None,
                            };
                            (record.name.clone(), percentage)
                        })
                        .collect()
                })
                .collect()
        };
        let named =
            |name: &str, percentage: &str| vec![(name.to_string(), Some(percentage.to_string()))];
        let lines = ["21\t-5\t10\tbelow", "21\t99\t110\tregion_1"];
        assert_eq!(
            listed(&lines, ",fields=PC"),
            [
                named("below", "6.667"),
                named("below", "6.667"),
                vec![],
                named("region_1", "9.091"),
            ]
        );
        assert_eq!(
            listed(&["21\t0\t10\tbelow", "21\t99\t110\tregion_1"], ",fields=PC")[0],
            named("below", "10.000")
        );
        assert_eq!(
            listed(&lines, ",coords=1,fields=PC")[..2],
            [named("21:-4-10", "6.667"), named("21:-4-10", "6.667")]
        );
        assert_eq!(
            listed(&["21\t-5\t10", "21\t99\t110"], ",fields=PC")[1],
            named("21:-4-10", "6.667")
        );
    }

    /// On GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence
    /// and `chrM` and `M` another one. A file that names the contig `M` is
    /// read there for a record written `chrM` (the `chr` step of
    /// `BaseVEP.pm:728-729`) and for none written `MT` or `chrMT`, which
    /// `BaseVEP.pm:689-738` do not reach from either; the synonym `chrM` that
    /// a cache's table gives `MT` reaches a file naming `chrM`, as the same
    /// lines reach it. On GRCh38, and in a run without an assembly, the file
    /// naming `M` is read for a record written `MT` as well.
    #[test]
    fn on_grch37_a_mitochondrial_record_reads_a_file_of_the_other_pair_through_a_synonym_alone() {
        use vep_core::assembly::Assembly;
        let dir = tempfile::tempdir().unwrap();
        let named = |contig: &str| {
            write_tabix_fixture(
                dir.path(),
                &format!("regions_{contig}_synthetic.bed.gz"),
                &[format!("{contig}\t99\t110\tregion_1\t7")],
                IndexSpec::standard(),
            )
        };
        let reads = |file: &std::path::Path,
                     assembly: Option<Assembly>,
                     synonyms: &str,
                     record: &str|
         -> usize {
            let sources = CustomSources::open(
                vec![parse_custom(&format!("file={},format=bed", file.display())).unwrap()],
                Arc::new(ChromosomeSynonyms::parse(synonyms)),
                CustomRunOptions {
                    assembly,
                    ..CustomRunOptions::default()
                },
            )
            .unwrap();
            let mut batch = vec![InputVariant::new(
                record.into(),
                100,
                100,
                b"A".to_vec(),
                b"G".to_vec(),
            )];
            sources.run_batch(&mut batch).unwrap();
            batch[0]
                .custom_annotations
                .iter()
                .map(|annotations| annotations.records.len())
                .sum()
        };
        let grch37 = Some(Assembly::Grch37);
        let named_m = named("M");
        for (record, read) in [("MT", 0), ("chrMT", 0), ("chrM", 1), ("M", 1)] {
            assert_eq!(reads(&named_m, grch37, "", record), read, "{record}");
        }
        for assembly in [Some(Assembly::Grch38), None] {
            assert_eq!(reads(&named_m, assembly, "", "MT"), 1, "{assembly:?}");
        }
        let named_chr_m = named("chrM");
        assert_eq!(reads(&named_chr_m, grch37, "", "MT"), 0);
        assert_eq!(reads(&named_chr_m, grch37, "MT\tchrM\n", "MT"), 1);
    }

    /// A record written bare `M` is the one the resolver renames: a file that
    /// names the contig `MT`, and neither `M` nor `chrM`, is asked for `MT`
    /// (`Parser.pm:794-796`), with no synonym listed, on GRCh37 as on any
    /// assembly. A record written `chrM` reaches no such file there.
    #[test]
    fn on_grch37_a_record_written_bare_m_reads_a_file_naming_mt_through_the_rename() {
        use vep_core::assembly::Assembly;
        let dir = tempfile::tempdir().unwrap();
        let file = write_tabix_fixture(
            dir.path(),
            "regions_MT_synthetic.bed.gz",
            &["MT\t99\t110\tregion_1\t7"],
            IndexSpec::standard(),
        );
        let reads = |record: &str| -> usize {
            let sources = CustomSources::open(
                vec![parse_custom(&format!("file={},format=bed", file.display())).unwrap()],
                Arc::default(),
                CustomRunOptions {
                    assembly: Some(Assembly::Grch37),
                    ..CustomRunOptions::default()
                },
            )
            .unwrap();
            let mut batch = vec![InputVariant::new(
                record.into(),
                100,
                100,
                b"A".to_vec(),
                b"G".to_vec(),
            )];
            sources.run_batch(&mut batch).unwrap();
            batch[0]
                .custom_annotations
                .iter()
                .map(|annotations| annotations.records.len())
                .sum()
        };
        assert_eq!(reads("M"), 1);
        assert_eq!(reads("chrM"), 0);
        assert_eq!(reads("MT"), 1);
    }
}
