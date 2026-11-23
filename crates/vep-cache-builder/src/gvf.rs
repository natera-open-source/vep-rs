// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Known variants read from a GVF file of Ensembl's variation export, each
//! feature as the record the variation database holds of it.
//!
//! The export writes a variation feature as one line of nine tab-separated
//! columns (`scripts/export/release/dump_gvf.pl` of ensembl-variation): the
//! sequence, the source, the class, the start, the end, a score, the strand,
//! a phase and the attributes, `key=value` pairs joined by `;`. This module
//! reads back what the export wrote from the database's own values:
//!
//! - the name, the text after the source in `Dbxref` (`add_source`:
//!   `<source>_<version>:<name>`, or `<source>:<name>` for a source without a
//!   version);
//! - the span: the export collapses a start past the end onto the end
//!   (`add_coords`), which only an insertion has, so an `insertion` feature
//!   starts one base after its end;
//! - the allele string (`add_alleles`): `Reference_seq`, unless it is the `.`
//!   the export writes for an allele string of one allele or none, then the
//!   alleles of `Variant_seq`, a `.` among them standing for the next string
//!   of `invalid_variant_strings` (an allele that is no sequence, such as a
//!   class name), all joined by `/`. A feature whose `comment` is
//!   `HGMD_MUTATION` has that allele string.
//!
//! Two things the export does not keep are read as written: a reference
//! allele of more than 50 bases is the `~` the export writes for it, and a
//! tandem-repeat allele is the sequence the export expanded it to.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;

use crate::core_dumps::BoundedLines;
use crate::line::MAX_LINE_BYTES;
use crate::variation::VariationRecord;

/// The allele string the variation database holds for an HGMD-PUBLIC variant,
/// which the export names in a feature's `comment`.
const HGMD_MUTATION: &str = "HGMD_MUTATION";

/// One feature of a GVF file as a variation record.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GvfFeature {
    /// The sequence the feature lies on.
    pub(crate) sequence: String,
    pub(crate) record: VariationRecord,
    /// The feature's start column, the position its alleles are matched to a
    /// frequency file's at: `None` for a feature without a reference allele
    /// or without an alternate one, and for one whose reference allele the
    /// export did not write.
    pub(crate) allele_position: Option<u64>,
}

/// Reads the features of the GVF at `path`, gzip when its name ends `.gz`,
/// and hands each to `feature` in file order; lines opening with `#` and
/// empty lines are passed over. A line that is no feature (fewer than nine
/// columns, a start or an end that is no number, no `Dbxref` naming the
/// variant) is an error naming the file and the line. Returns the number of
/// features read.
pub(crate) fn read_features(
    path: &Path,
    mut feature: impl FnMut(GvfFeature) -> Result<()>,
) -> Result<u64> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };
    let mut features = 0u64;
    let lines = BoundedLines::new(reader, MAX_LINE_BYTES, path.display().to_string());
    for (index, line) in lines.enumerate() {
        let line = line?;
        let line = String::from_utf8_lossy(&line);
        let line = line.strip_suffix('\r').unwrap_or(&line);
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let read = feature_of(line)
            .map_err(|reason| anyhow!("{} line {}: {reason}", path.display(), index + 1))?;
        features += 1;
        feature(read)?;
    }
    Ok(features)
}

/// The record of one feature line, or what the line lacks.
fn feature_of(line: &str) -> std::result::Result<GvfFeature, String> {
    let columns: Vec<&str> = line.split('\t').collect();
    let [sequence, _source, class, start, end, _score, strand, _phase, attributes] =
        columns.as_slice()
    else {
        return Err(format!(
            "a feature has 9 columns, this line {}",
            columns.len()
        ));
    };
    let number = |column: &str, value: &str| -> std::result::Result<u64, String> {
        value
            .parse()
            .map_err(|_| format!("the {column} {value:?} is not a number"))
    };
    let written_start = number("start", start)?;
    let end = number("end", end)?;

    let attributes: HashMap<&str, &str> = attributes
        .split(';')
        .filter_map(|pair| pair.split_once('='))
        .collect();
    let name = attributes
        .get("Dbxref")
        .and_then(|dbxref| dbxref.split_once(':'))
        .map(|(_, name)| name)
        .filter(|name| !name.is_empty())
        .ok_or("no Dbxref attribute names the variant")?;

    let reference = attributes
        .get("Reference_seq")
        .copied()
        .filter(|reference| *reference != ".");
    let mut invalid = attributes
        .get("invalid_variant_strings")
        .into_iter()
        .flat_map(|strings| strings.split(','));
    let variants: Vec<&str> = attributes
        .get("Variant_seq")
        .filter(|variants| !variants.is_empty())
        .into_iter()
        .flat_map(|variants| variants.split(','))
        .map(|variant| match variant {
            "." => invalid.next().unwrap_or(variant),
            _ => variant,
        })
        .collect();
    let hgmd = attributes.get("comment") == Some(&HGMD_MUTATION);
    let allele_string = if hgmd {
        HGMD_MUTATION.to_string()
    } else {
        let alleles: Vec<&str> = reference
            .into_iter()
            .chain(variants.iter().copied())
            .collect();
        alleles.join("/")
    };
    let matched_by_allele =
        !hgmd && reference.is_some_and(|reference| reference != "~") && !variants.is_empty();

    Ok(GvfFeature {
        sequence: sequence.to_string(),
        record: VariationRecord {
            variation_name: name.to_string(),
            start: if *class == "insertion" {
                end + 1
            } else {
                written_start
            },
            end,
            strand: if *strand == "-" { -1 } else { 1 },
            allele_string,
            ..Default::default()
        },
        allele_position: matched_by_allele.then_some(written_start),
    })
}

/// The failed variants of a genome, from the GVF file a release publishes of
/// them (`<species>_failed.gvf.gz`): every variation feature of a variant
/// the variation database marks as failed, which Ensembl VEP's cache holds
/// with `failed` 1 and the variation VCFs leave out unless the variant is
/// kept on display.
#[derive(Debug, Default)]
pub(crate) struct FailedVariants {
    /// The names of the file's variants, whatever sequence they lie on.
    names: HashSet<String>,
    /// The features of each sequence kept, by start then name, one a variant
    /// and place.
    by_sequence: HashMap<String, Vec<GvfFeature>>,
    /// The features the file holds.
    pub(crate) features: u64,
    /// The features on a sequence that is not kept.
    pub(crate) elsewhere: u64,
}

impl FailedVariants {
    /// Reads the file at `path`, keeping every variant's name and the
    /// features on the sequences `chromosomes` names, each as a record with
    /// `failed` 1. A feature written twice is kept once.
    pub(crate) fn load(path: &Path, chromosomes: &[String]) -> Result<Self> {
        let mut failed = Self::default();
        failed.features = read_features(path, |mut feature| {
            failed.names.insert(feature.record.variation_name.clone());
            if !chromosomes.contains(&feature.sequence) {
                failed.elsewhere += 1;
                return Ok(());
            }
            feature.record.failed = 1;
            failed
                .by_sequence
                .entry(feature.sequence.clone())
                .or_default()
                .push(feature);
            Ok(())
        })?;
        for features in failed.by_sequence.values_mut() {
            let place = |feature: &GvfFeature| {
                let record = &feature.record;
                (
                    record.start,
                    record.variation_name.clone(),
                    record.end,
                    record.allele_string.clone(),
                )
            };
            features.sort_by_key(place);
            features.dedup_by_key(|feature| place(feature));
        }
        Ok(failed)
    }

    /// Whether the file names the variant `name`.
    pub(crate) fn names(&self, name: &str) -> bool {
        self.names.contains(name)
    }

    /// The features kept on `sequence`, by start then name.
    pub(crate) fn on(&self, sequence: &str) -> &[GvfFeature] {
        self.by_sequence.get(sequence).map_or(&[], Vec::as_slice)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A feature line of `class` at `start`-`end` on sequence `7` with
    /// `attributes`.
    fn line(class: &str, start: u64, end: u64, attributes: &str) -> String {
        format!("7\tmadeup\t{class}\t{start}\t{end}\t.\t+\t.\t{attributes}")
    }

    fn feature(class: &str, start: u64, end: u64, attributes: &str) -> GvfFeature {
        feature_of(&line(class, start, end, attributes)).unwrap()
    }

    /// The name, span, strand and allele string of a feature's record, and
    /// the position its alleles are matched at.
    fn read(feature: &GvfFeature) -> (&str, u64, u64, i8, &str, Option<u64>) {
        let record = &feature.record;
        (
            record.variation_name.as_str(),
            record.start,
            record.end,
            record.strand,
            record.allele_string.as_str(),
            feature.allele_position,
        )
    }

    /// A feature is the database's record of it. A substitution and a
    /// deletion keep the span written; an insertion, whose start the export
    /// collapsed onto its end, starts one base after its end. The allele
    /// string is the reference allele and the variant alleles, the commas of
    /// `Variant_seq` read `/`. A name holding a colon is the text after the
    /// first one, and a source without a version names it the same way. The
    /// attributes are read in any order, and a feature on the reverse strand
    /// keeps that strand.
    #[test]
    fn a_feature_reads_as_the_database_record_it_was_written_from() {
        let snv = feature(
            "SNV",
            500,
            500,
            "ID=1;Variant_seq=C,T;Dbxref=madeup_7:mv500;Reference_seq=A",
        );
        assert_eq!(read(&snv), ("mv500", 500, 500, 1, "A/C/T", Some(500)));
        assert_eq!(snv.sequence, "7");
        assert_eq!((snv.record.failed, snv.record.somatic), (0, 0));

        let deletion = feature(
            "deletion",
            600,
            602,
            "Dbxref=madeup:mv:600;Reference_seq=CAT;Variant_seq=-;ID=2",
        );
        assert_eq!(read(&deletion), ("mv:600", 600, 602, 1, "CAT/-", Some(600)));

        let insertion = feature(
            "insertion",
            700,
            700,
            "ID=3;Variant_seq=GG;Dbxref=madeup_7:mv700;Reference_seq=-",
        );
        assert_eq!(read(&insertion), ("mv700", 701, 700, 1, "-/GG", Some(700)));

        let reverse = feature_of(
            "7\tmadeup\tSNV\t800\t800\t.\t-\t.\tID=4;Variant_seq=C;Dbxref=madeup_7:mv800;Reference_seq=A",
        )
        .unwrap();
        assert_eq!(read(&reverse), ("mv800", 800, 800, -1, "A/C", Some(800)));
    }

    /// An allele string of one allele, or of none, has no reference allele:
    /// the export writes `.` for it, and the record's allele string is the
    /// variant alleles alone, empty where `Variant_seq` is. An allele that is
    /// no sequence is written `.` in `Variant_seq` and named, in order, in
    /// `invalid_variant_strings`. A reference allele past 50 bases is the `~`
    /// the export writes. An HGMD-PUBLIC feature has the allele string its
    /// comment names. Of these only the feature that has a reference allele
    /// written as bases has alleles a frequency file's are matched to.
    #[test]
    fn alleles_the_export_wrote_another_way_are_read_back() {
        let one = feature(
            "sequence_alteration",
            10,
            10,
            "ID=1;Variant_seq=T;Dbxref=madeup:mv10;Reference_seq=.",
        );
        assert_eq!(read(&one), ("mv10", 10, 10, 1, "T", None));

        let none = feature(
            "sequence_alteration",
            20,
            20,
            "ID=2;Variant_seq=;Dbxref=madeup:mv20;Reference_seq=.",
        );
        assert_eq!(read(&none), ("mv20", 20, 20, 1, "", None));

        let class_name = feature(
            "sequence_alteration",
            30,
            34,
            "ID=3;Variant_seq=.;invalid_variant_strings=madeup_variant;Dbxref=madeup_7:mv30;Reference_seq=.",
        );
        assert_eq!(
            read(&class_name),
            ("mv30", 30, 34, 1, "madeup_variant", None)
        );

        let mixed = feature(
            "sequence_alteration",
            40,
            40,
            "ID=4;Variant_seq=.,G,.;invalid_variant_strings=(LONG)1,(LONG)2;Dbxref=madeup_7:mv40;Reference_seq=A",
        );
        assert_eq!(
            read(&mixed),
            ("mv40", 40, 40, 1, "A/(LONG)1/G/(LONG)2", Some(40))
        );

        let long = feature(
            "deletion",
            50,
            110,
            "ID=5;Variant_seq=-;Dbxref=madeup_7:mv50;Reference_seq=~",
        );
        assert_eq!(read(&long), ("mv50", 50, 110, 1, "~/-", None));

        let hgmd = feature(
            "insertion",
            60,
            60,
            "ID=6;Variant_seq=.;comment=HGMD_MUTATION;Dbxref=madeup_7:mv60;Reference_seq=-",
        );
        assert_eq!(read(&hgmd), ("mv60", 61, 60, 1, "HGMD_MUTATION", None));
    }

    /// A line that is no feature refuses the read, naming the file, the line
    /// and what the line lacks; header lines and empty lines are passed
    /// over.
    #[test]
    fn a_line_that_is_no_feature_is_refused_by_file_and_line() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("made_up_failed.gvf");
        let good = line(
            "SNV",
            5,
            5,
            "ID=1;Variant_seq=C;Dbxref=madeup_7:mv5;Reference_seq=A",
        );
        for (bad, reason) in [
            (
                "7\tmadeup\tSNV\t5\t5\t.\t+\t.".to_string(),
                "a feature has 9 columns, this line 8",
            ),
            (
                line("SNV", 5, 5, "ID=1;Variant_seq=C;Reference_seq=A"),
                "no Dbxref attribute names the variant",
            ),
            (
                line(
                    "SNV",
                    5,
                    5,
                    "ID=1;Variant_seq=C;Dbxref=madeup_7;Reference_seq=A",
                ),
                "no Dbxref attribute names the variant",
            ),
            (
                good.replace("\t5\t5\t", "\tfive\t5\t"),
                "the start \"five\" is not a number",
            ),
            (
                good.replace("\t5\t5\t", "\t5\t\t"),
                "the end \"\" is not a number",
            ),
        ] {
            std::fs::write(&path, format!("##gff-version 3\n\n{good}\r\n{bad}\n")).unwrap();
            let mut read = Vec::new();
            let error = read_features(&path, |feature| {
                read.push(feature.record.variation_name);
                Ok(())
            })
            .unwrap_err();
            assert_eq!(
                format!("{error:#}"),
                format!("{} line 4: {reason}", path.display())
            );
            assert_eq!(read, ["mv5"]);
        }
    }

    /// The failed variants of a file: every name, whatever its sequence, and
    /// the features of the sequences asked for, each flagged failed, in
    /// start order then name order, a feature written twice kept once and a
    /// variant with two places kept at both. A gzip file reads as a plain
    /// one.
    #[test]
    fn failed_variants_keep_every_name_and_the_features_of_the_sequences_asked_for() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("made_up_failed.gvf.gz");
        let on = |sequence: &str, start: u64, name: &str| {
            format!(
                "{sequence}\tmadeup\tSNV\t{start}\t{start}\t.\t+\t.\tID=1;Variant_seq=C;Dbxref=madeup_7:{name};Reference_seq=A;ensembl_failure_reason=a made-up reason\n"
            )
        };
        let text = [
            "##gff-version 3\n##gvf-version 1.07\n".to_string(),
            on("7", 900, "mv900"),
            on("7", 300, "mvb"),
            on("7", 300, "mva"),
            on("7", 900, "mv900"),
            on("7", 100, "mv900"),
            on("8", 50, "mv50"),
            on("9", 60, "mv60"),
        ]
        .concat();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, text.as_bytes()).unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();

        let failed = FailedVariants::load(&path, &["7".to_string(), "9".to_string()]).unwrap();
        assert_eq!((failed.features, failed.elsewhere), (7, 1));
        for name in ["mv900", "mva", "mvb", "mv50", "mv60"] {
            assert!(failed.names(name), "{name}");
        }
        assert!(!failed.names("mv"));
        let placed = |sequence: &str| -> Vec<(u64, &str, u8)> {
            failed
                .on(sequence)
                .iter()
                .map(|feature| {
                    let record = &feature.record;
                    (record.start, record.variation_name.as_str(), record.failed)
                })
                .collect()
        };
        assert_eq!(
            placed("7"),
            [
                (100, "mv900", 1),
                (300, "mva", 1),
                (300, "mvb", 1),
                (900, "mv900", 1)
            ]
        );
        assert_eq!(placed("9"), [(60, "mv60", 1)]);
        assert!(placed("8").is_empty());
    }
}
