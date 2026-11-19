// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` source over a tabix-indexed GFF in `gff_type=gencode_promoter`
//! mode, read as `AnnotationSource/File/GFF.pm` reads one through ensembl-io's
//! `Parser/GFF3Tabix.pm`: the one GFF mode that overlays its records on the
//! input instead of standing in for the transcript cache (`GFF.pm:154-164`).
//!
//! A record spans its fourth to its fifth column and matches through the base
//! predicates (`File.pm:534-575`). A matching `promoter_window` record
//! (`GFF.pm:119-121`) becomes one output record named by its coordinates
//! (`File.pm:483-497`: the parser has no name accessor) with three fields:
//! `type`, `feature_id` (the `ID` attribute) and `associated_gene` (the
//! `gene_id` attribute, else `.`) (`GFF.pm:181-204`), the two attribute values
//! with every `%XX` escape decoded over the record's bytes, as the other GFF3
//! readers decode them, where Ensembl VEP prints them as written (the raw hash
//! of `GFF3.pm:190-198`). A matching record of any other type yields nothing
//! and raises Ensembl VEP's warning (`:191-195`).

use std::borrow::Cow;
use std::fmt::Write as _;
use std::sync::Arc;

use indexmap::IndexMap;
use tracing::warn;
use vep_core::custom_annotation::{CustomFieldValue, CustomRecord};
use vep_core::regulatory_gff::{decode_attribute, decode_type};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::InputVariant;

use super::grammar::{truthy, CustomSpec};
use super::source::{
    base_overlap, CustomError, CustomSource, FormatRules, RecordGroup, TabixSource, Verdict,
};
use crate::tabix::TabixRecord;
use crate::traits::PrefetchData;

/// The feature types a `gencode_promoter` source annotates
/// (`%INCLUDE_FEATURE_TYPES_GENCODE_PROMOTER`, `GFF.pm:119-121`).
const PROMOTER_FEATURE_TYPES: [&str; 1] = ["promoter_window"];

/// The 0-based column of a GFF record's type.
const TYPE_COLUMN: usize = 2;

/// The 0-based column of a GFF record's attributes.
const ATTRIBUTES_COLUMN: usize = 8;

/// A `--custom` source of `format=gff,gff_type=gencode_promoter`.
pub struct GffOverlaySource {
    base: TabixSource,
}

impl GffOverlaySource {
    /// The reader builds a record's type cell from the line's bytes through
    /// [`decode_type`], the `get_type` of `BaseGXF.pm:288` (`GFF3.pm:161-164`):
    /// the escapes decoded over the bytes, then read as text, so an escape
    /// completing a raw first byte of a character forms the character as
    /// Ensembl VEP's bytes do. The attributes cell is the line's bytes as
    /// [`attributes_text`] spells them, split as `get_attributes` splits them
    /// (`GFF3.pm:190-198`) and decoded over those bytes where a value is read
    /// ([`FormatRules::match_record`]).
    pub fn open(spec: CustomSpec, synonyms: Arc<ChromosomeSynonyms>) -> Result<Self, CustomError> {
        Ok(GffOverlaySource {
            base: TabixSource::open(spec, synonyms)?.with_reader(|reader| {
                reader
                    .decoding_column_from_bytes(TYPE_COLUMN, decode_type)
                    .decoding_column_from_bytes(ATTRIBUTES_COLUMN, attributes_text)
            }),
        })
    }
}

/// The attributes column as text that keeps the column's bytes for
/// [`decode_attribute`]: a UTF-8 column borrowed as it is; in a column holding
/// bytes that are not UTF-8, each such byte spelled as its `%XX` escape, which
/// the decoder reads back to the byte. A value's escapes then decode over the
/// bytes the file holds, as the `--regulatory_gff` reader decodes them: an
/// escape beside a raw byte completes the character, and a raw byte that is no
/// UTF-8 once decoded reads as one U+FFFD. The escapes written in the column
/// stay written, so the split at `;` and `=` and the truth test read the value
/// as Ensembl VEP's hash holds it (`GFF3.pm:190-198`).
pub(crate) fn attributes_text(column: &[u8]) -> Cow<'_, str> {
    if let Ok(text) = std::str::from_utf8(column) {
        return Cow::Borrowed(text);
    }
    let mut text = String::with_capacity(column.len() + 2);
    for chunk in column.utf8_chunks() {
        text.push_str(chunk.valid());
        for byte in chunk.invalid() {
            let _ = write!(text, "%{byte:02X}");
        }
    }
    Cow::Owned(text)
}

impl CustomSource for GffOverlaySource {
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

impl FormatRules for GffOverlaySource {
    /// The base `_record_overlaps_VF` on the record's fourth to fifth column
    /// (`File.pm:523, 534-575`), the walk ending at a record starting past the
    /// feature's end plus one (`:333-335`), then `GFF.pm:181-204`: a record of
    /// a type outside [`PROMOTER_FEATURE_TYPES`] warns and yields nothing; a
    /// promoter window is named `<seqname>:<start>-<end>`, its seqname without
    /// a leading `chr` (`File.pm:489-495`, `GXF.pm:88-94`), its start and end
    /// the whole numbers of the fourth and fifth columns as Perl reads them
    /// ([`TabixSource::start_of`], [`TabixSource::end_of`]), and carries
    /// `type`, `feature_id` (`ID`, undefined when the attribute is not there)
    /// and `associated_gene` (`gene_id` when true in Perl's sense on the value
    /// as written, else `.`), the two values then decoded over the record's
    /// bytes ([`attributes_text`], [`decode_attribute`]), which Ensembl VEP
    /// leaves as written (`GFF.pm:198-201` over the raw hash of
    /// `GFF3.pm:190-198`). No score is read, so the statistics count the
    /// records alone.
    fn match_record(
        &self,
        record: &TabixRecord,
        group: &RecordGroup<'_>,
        _scores: bool,
    ) -> Verdict {
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
        let column = |index: usize| record.columns.get(index).map_or("", String::as_str);
        let feature_type = column(TYPE_COLUMN);
        if !PROMOTER_FEATURE_TYPES.contains(&feature_type) {
            warn!(
                "Ignoring '{feature_type}' feature_type from {} GFF/GTF file. \
                 This feature_type is not supported in VEP.",
                self.base.spec().file
            );
            return Verdict::NoMatch;
        }
        let attributes = column(ATTRIBUTES_COLUMN);
        let text = |value: &str| CustomFieldValue::Text(value.to_string());
        let decoded =
            |value: &str| CustomFieldValue::Text(decode_attribute(value.as_bytes()).into_owned());
        let feature_id = attribute(attributes, "ID").map_or(CustomFieldValue::Undefined, decoded);
        let associated_gene = match attribute(attributes, "gene_id") {
            Some(gene_id) if truthy(gene_id) => decoded(gene_id),
            _ => text("."),
        };
        let fields: IndexMap<String, CustomFieldValue> = [
            ("type", text(feature_type)),
            ("feature_id", feature_id),
            ("associated_gene", associated_gene),
        ]
        .into_iter()
        .map(|(field, value)| (field.to_string(), value))
        .collect();
        Verdict::Match(
            vec![CustomRecord {
                name: format!("{}:{}-{}", seqname(&record.chr), start.whole(), end.whole()),
                allele: None,
                fields,
                score: None,
            }],
            Some(percentage),
        )
    }
}

/// `get_seqname` (`GXF.pm:88-94`): the sequence name without a leading `chr`.
fn seqname(chr: &str) -> &str {
    chr.strip_prefix("chr").unwrap_or(chr)
}

/// The value of `key` among the record's attributes as ensembl-io reads the
/// ninth column (`GFF3.pm:190-198`): the pairs split at `;`, a pair's key the
/// text before its first `=` and its value the text between its first and
/// second `=`, a later pair with the key replacing an earlier one. `None` for a
/// key that is not there or whose last pair has no `=`.
fn attribute<'a>(attributes: &'a str, key: &str) -> Option<&'a str> {
    let mut value = None;
    for pair in attributes.split(';') {
        let mut pieces = pair.splitn(3, '=');
        if pieces.next() == Some(key) {
            value = pieces.next();
        }
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::custom::grammar::parse_custom;
    use crate::custom::source::CustomSources;
    use crate::custom::CustomRunOptions;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec, Warnings};
    use vep_core::custom_annotation::{CustomAnnotations, CustomStats, CustomValue};

    /// The index layout of `tabix -p gff`: sequence in column 1, span in
    /// columns 4 and 5.
    fn gff_index() -> IndexSpec {
        IndexSpec::standard().with_span_cols(1, 4, 5)
    }

    fn open(argument: &str) -> CustomSources {
        CustomSources::open(
            vec![parse_custom(argument).unwrap()],
            Arc::default(),
            CustomRunOptions::default(),
        )
        .unwrap()
    }

    fn variant(chr: &str, start: u64, end: u64, record: u64) -> InputVariant {
        let mut variant = InputVariant::new(chr.into(), start, end, b"A".to_vec(), b"G".to_vec());
        variant.input_record = record;
        variant
    }

    fn only(variant: &InputVariant) -> &CustomAnnotations {
        assert_eq!(variant.custom_annotations.len(), 1, "{variant:?}");
        &variant.custom_annotations[0]
    }

    fn text_fields(record: &CustomRecord) -> Vec<(&str, &str)> {
        record
            .fields
            .iter()
            .map(|(field, value)| match value {
                CustomFieldValue::Text(text) => (field.as_str(), text.as_str()),
                CustomFieldValue::List(_) => panic!("{field} is a list"),
                CustomFieldValue::Undefined => panic!("{field} is undefined"),
            })
            .collect()
    }

    /// `GFF.pm:181-204` on the record `21 . promoter_window 90 110 . + .
    /// ID=PW1;gene_id=ENSG00000141510` over a variant at 100: one record under
    /// `GENCODE_PROMOTER`, named `21:90-110` as `File.pm:489-495` names a record
    /// whose parser has no `get_name` (Ensembl VEP's own test expects
    /// `21:25584733-25595733` for its promoter fixture), with `type`
    /// `promoter_window`, `feature_id` `PW1` and `associated_gene`
    /// `ENSG00000141510`; a window without `gene_id` gives `.` (`:200`) and one
    /// without `ID` an undefined `feature_id` (`:198`). A `chr21` seqname loses
    /// its prefix in the name alone (`GXF.pm:92`). A variant at 89, one base before
    /// the window, matches nothing (`File.pm:569`). Under `summary_stats` the
    /// records count but add no score (`GFF.pm:186`: no `score` is set).
    #[test]
    fn promoter_windows_annotate_as_overlap_records_with_three_fields() {
        let dir = tempfile::tempdir().unwrap();
        let gff = write_tabix_fixture(
            dir.path(),
            "promoters.gff.gz",
            &[
                "##gff-version 3",
                "21\t.\tpromoter_window\t90\t110\t.\t+\t.\tID=PW1;gene_id=ENSG00000141510",
                "21\t.\tpromoter_window\t200\t210\t.\t-\t.\tID=PW2",
                "21\t.\tpromoter_window\t300\t310\t.\t+\t.\tgene_id=ENSG00000000003;Name=PW3",
            ],
            gff_index(),
        );
        let sources = open(&format!(
            "file={},format=gff,gff_type=gencode_promoter",
            gff.display()
        ));
        let mut batch = vec![
            variant("21", 89, 89, 1),
            variant("21", 100, 100, 2),
            variant("21", 205, 205, 3),
            variant("21", 305, 305, 4),
        ];
        sources.run_batch(&mut batch).unwrap();

        assert!(batch[0].custom_annotations.is_empty());
        let annotation = only(&batch[1]);
        assert_eq!(annotation.short_name.as_ref(), "GENCODE_PROMOTER");
        assert_eq!(annotation.records.len(), 1);
        let window = &annotation.records[0];
        assert_eq!(window.name, "21:90-110");
        assert_eq!(window.allele, None);
        assert_eq!(window.score, None);
        assert_eq!(
            text_fields(window),
            [
                ("type", "promoter_window"),
                ("feature_id", "PW1"),
                ("associated_gene", "ENSG00000141510"),
            ]
        );
        assert_eq!(
            text_fields(&only(&batch[2]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", "PW2"),
                ("associated_gene", "."),
            ]
        );
        let unnamed = &only(&batch[3]).records[0];
        assert_eq!(
            unnamed.fields["type"],
            CustomFieldValue::Text("promoter_window".into())
        );
        assert_eq!(unnamed.fields["feature_id"], CustomFieldValue::Undefined);
        assert_eq!(
            unnamed.fields["associated_gene"],
            CustomFieldValue::Text("ENSG00000000003".into())
        );

        let prefixed = write_tabix_fixture(
            dir.path(),
            "promoters_chr.gff.gz",
            &["chr21\t.\tpromoter_window\t90\t110\t.\t+\t.\tID=PW1;gene_id=ENSG00000141510"],
            gff_index(),
        );
        let sources = open(&format!(
            "file={},format=gff,gff_type=gencode_promoter,summary_stats=count%max",
            prefixed.display()
        ));
        let mut batch = vec![variant("21", 100, 100, 1)];
        sources.run_batch(&mut batch).unwrap();
        let annotation = only(&batch[0]);
        assert_eq!(annotation.records[0].name, "21:90-110");
        assert_eq!(
            annotation.stats,
            Some(CustomStats::Values(
                [("count".to_string(), CustomValue::Number(1.0))]
                    .into_iter()
                    .collect()
            ))
        );
    }

    /// The names each record of `batch` lists, in the batch's order.
    fn names(batch: &[InputVariant]) -> Vec<Vec<&str>> {
        batch
            .iter()
            .map(|variant| {
                variant
                    .custom_annotations
                    .iter()
                    .flat_map(|annotation| &annotation.records)
                    .map(|record| record.name.as_str())
                    .collect()
            })
            .collect()
    }

    /// A window whose start a minus sign leads is listed as Ensembl VEP 116.2
    /// lists it, on a file's first sequence and on a later one: for every
    /// record from the first base, where `tabix` indexes the window (htslib
    /// 1.9 `tbx.c:99`), to the window's end and for none past it, under the
    /// number Perl reads from the column as written (`GXF.pm:163-166`,
    /// `File.pm:489-495`): the decimal digits after the sign, after any blank,
    /// where no fraction or exponent follows them. Two such windows of one
    /// sequence are listed in the file's order. The match reads the same
    /// number (`File.pm:523`), so a window from `-5` to 10 is 16 bases to
    /// `overlap_cutoff` and eight bases from a record at 3 to `distance`. A
    /// window from -9223372036854775807 is matched on a span that saturates
    /// and listed as Ensembl VEP lists that one; below -2^31 Ensembl VEP's
    /// overlap test reads the low 32 bits of a start alone.
    #[test]
    fn a_window_whose_start_a_minus_sign_leads_is_listed_under_the_number_written() {
        let dir = tempfile::tempdir().unwrap();
        let fixture = |lines: &[String]| {
            let gff = write_tabix_fixture(dir.path(), "promoters.gff.gz", lines, gff_index());
            format!(
                "file={},format=gff,gff_type=gencode_promoter",
                gff.display()
            )
        };
        let window = |chr: &str, start: &str, end: u64, id: &str| {
            format!("{chr}\t.\tpromoter_window\t{start}\t{end}\t.\t+\t.\tID={id}")
        };
        for (written, read) in [
            ("-5", "-5"),
            (" -5", "-5"),
            ("-5x", "-5"),
            ("-05", "-5"),
            ("-010", "-10"),
            ("-0", "0"),
            ("-0x1A", "0"),
            ("-9223372036854775807", "-9223372036854775807"),
        ] {
            let sources = open(&fixture(&[
                window("20", written, 500, "PW20"),
                window("20", "1000", 1010, "SOUND20"),
                window("21", written, 400, "PW21"),
                window("21", "1000", 1010, "SOUND21"),
            ]));
            let mut batch = vec![
                variant("20", 1, 1, 1),
                variant("20", 250, 250, 2),
                variant("20", 500, 500, 3),
                variant("20", 501, 501, 4),
                variant("20", 1005, 1005, 5),
                variant("21", 1, 1, 6),
                variant("21", 400, 400, 7),
                variant("21", 401, 401, 8),
                variant("21", 1005, 1005, 9),
            ];
            sources.run_batch(&mut batch).unwrap();
            let (first, later) = (format!("20:{read}-500"), format!("21:{read}-400"));
            assert_eq!(
                names(&batch),
                [
                    vec![first.as_str()],
                    vec![first.as_str()],
                    vec![first.as_str()],
                    vec![],
                    vec!["20:1000-1010"],
                    vec![later.as_str()],
                    vec![later.as_str()],
                    vec![],
                    vec!["21:1000-1010"],
                ],
                "{written:?}"
            );
            assert_eq!(
                text_fields(&only(&batch[0]).records[0]),
                [
                    ("type", "promoter_window"),
                    ("feature_id", "PW20"),
                    ("associated_gene", "."),
                ]
            );
        }

        let sources = open(&fixture(&[
            window("20", "90", 110, "PW20"),
            window("21", "-5", 400, "PWA"),
            window("21", "-3", 1010, "PWB"),
        ]));
        let mut batch = vec![
            variant("20", 100, 100, 1),
            variant("21", 1, 1, 2),
            variant("21", 401, 401, 3),
        ];
        sources.run_batch(&mut batch).unwrap();
        assert_eq!(
            names(&batch),
            [
                vec!["20:90-110"],
                vec!["21:-5-400", "21:-3-1010"],
                vec!["21:-3-1010"],
            ]
        );

        let short = fixture(&[
            window("21", "-5", 10, "SHORT"),
            window("21", "20", 29, "SOUND"),
        ]);
        for (options, listed) in [
            ("", [vec!["21:-5-10"], vec!["21:-5-10"], vec!["21:20-29"]]),
            (",overlap_cutoff=8", [vec![], vec![], vec!["21:20-29"]]),
            (
                ",overlap_cutoff=6.25",
                [vec!["21:-5-10"], vec!["21:-5-10"], vec!["21:20-29"]],
            ),
            (",distance=7", [vec![], vec![], vec!["21:20-29"]]),
            (",distance=8", [vec!["21:-5-10"], vec![], vec!["21:20-29"]]),
        ] {
            let sources = open(&format!("{short}{options}"));
            let mut batch = vec![
                variant("21", 3, 3, 1),
                variant("21", 4, 4, 2),
                variant("21", 25, 25, 3),
            ];
            sources.run_batch(&mut batch).unwrap();
            assert_eq!(names(&batch), listed, "{options:?}");
        }

        let far = fixture(&[window("21", "-9223372036854775807", 10, "FAR")]);
        let sources = open(&format!("{far},distance=8"));
        let mut batch = vec![variant("21", 3, 3, 1)];
        sources.run_batch(&mut batch).unwrap();
        assert_eq!(names(&batch), [Vec::<&str>::new()]);
    }

    /// `GFF.pm:190-195`: a matching record of a type outside the promoter list
    /// yields no record and warns `Ignoring '<type>' feature_type from <file>
    /// GFF/GTF file. This feature_type is not supported in VEP.`, the file as
    /// the string gave it; the input record keeps the promoter windows it
    /// matches and the warning counts nothing against `num_records`
    /// (`File.pm:402, 339`). A record that does not match is passed over
    /// without a warning (`:394-396`), and the type is read percent-decoded
    /// (`GFF3.pm:161-164`).
    #[test]
    fn other_feature_types_warn_and_are_skipped() {
        let dir = tempfile::tempdir().unwrap();
        let gff = write_tabix_fixture(
            dir.path(),
            "promoters.gff.gz",
            &[
                "21\t.\tenhancer\t80\t105\t.\t+\t.\tID=EH1;gene_id=ENSG00000141510",
                "21\t.\tpromoter_window\t90\t110\t.\t+\t.\tID=PW1;gene_id=ENSG00000141510",
                "21\t.\tenhancer\t200\t210\t.\t+\t.\tID=EH2",
                "21\t.\tpromoter%5Fwindow\t300\t310\t.\t+\t.\tID=PW3",
            ],
            gff_index(),
        );
        let argument = format!(
            "file={},format=gff,gff_type=gencode_promoter,num_records=1",
            gff.display()
        );
        let sources = open(&argument);
        let mut batch = vec![
            variant("21", 100, 100, 1),
            variant("21", 205, 205, 2),
            variant("21", 305, 305, 3),
        ];
        sources.run_batch(&mut batch).unwrap();
        let names = |variant: &InputVariant| -> Vec<String> {
            only(variant)
                .records
                .iter()
                .map(|record| record.name.clone())
                .collect()
        };
        assert_eq!(names(&batch[0]), ["21:90-110"]);
        assert!(batch[1].custom_annotations.is_empty());
        assert_eq!(names(&batch[2]), ["21:300-310"]);
        assert_eq!(
            text_fields(&only(&batch[2]).records[0])[0],
            ("type", "promoter_window")
        );

        let source =
            GffOverlaySource::open(parse_custom(&argument).unwrap(), Arc::default()).unwrap();
        let record = |feature_type: &str, start: u64, end: u64| TabixRecord {
            chr: "21".to_string(),
            start,
            end,
            columns: vec![
                String::new(),
                ".".to_string(),
                feature_type.to_string(),
                start.to_string(),
                end.to_string(),
                ".".to_string(),
                "+".to_string(),
                ".".to_string(),
                "ID=EH1".to_string(),
            ],
            header: None,
        };
        let variants = [variant("21", 100, 100, 1)];
        let group = RecordGroup::of(&variants).unwrap();
        let warnings = Arc::new(Warnings::default());
        let verdicts = tracing::subscriber::with_default(Arc::clone(&warnings), || {
            [
                source.match_record(&record("enhancer", 80, 95), &group, false),
                source.match_record(&record("enhancer", 80, 105), &group, false),
                source.match_record(&record("enhancer", 102, 120), &group, false),
                source.match_record(&record("promoter_window", 90, 110), &group, false),
            ]
            .map(|verdict| match verdict {
                Verdict::PastEnd => "past end".to_string(),
                Verdict::NoMatch => "no match".to_string(),
                Verdict::Match(records, _) => records[0].name.clone(),
            })
        });
        assert_eq!(verdicts, ["no match", "no match", "past end", "21:90-110"]);
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [format!(
                "Ignoring 'enhancer' feature_type from {} GFF/GTF file. \
                 This feature_type is not supported in VEP.",
                gff.display()
            )]
        );
    }

    /// `GFF3.pm:161-164` through `BaseGXF.pm:288`: the type is decoded over
    /// the record's bytes, so an escape completing a raw first byte of a
    /// character forms the character, as Ensembl VEP's bytes do, and a byte
    /// that is not UTF-8, raw or escaped, is U+FFFD once in the warning. An
    /// attribute value is decoded over the record's bytes the same way, as
    /// `--regulatory_gff` decodes one: `%80` beside a raw first byte and `%CA`
    /// before a raw continuation byte each complete the character, a raw byte
    /// that completes none reads as one U+FFFD beside the character an escape
    /// spells, and two raw bytes completed by an escape form one character
    /// where the cell's text would hold two U+FFFD. Ensembl VEP writes the
    /// value's bytes (`GFF3.pm:190-198`). A value written in UTF-8 decodes as
    /// the cell's text does.
    #[test]
    fn a_type_column_is_decoded_over_the_records_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let lines: &[&[u8]] = &[
            b"21\t.\tcaf%FF\t90\t110\t.\t+\t.\tID=A",
            b"21\t.\tcaf\xFF\t90\t110\t.\t+\t.\tID=B",
            b"21\t.\tcaf%CA%80\t90\t110\t.\t+\t.\tID=C",
            b"21\t.\tmix\xCA%80\t90\t110\t.\t+\t.\tID=D",
            b"21\t.\tmix%CA\x80\t90\t110\t.\t+\t.\tID=E",
            b"21\t.\tpromoter%5Fwindow\t90\t110\t.\t+\t.\tID=caf\xCA%80;gene_id=caf%CA\x80",
            b"21\t.\tpromoter_window\t200\t210\t.\t+\t.\tID=caf\xFF%41;gene_id=\xE3\x80%80",
            b"21\t.\tpromoter_window\t300\t310\t.\t+\t.\tID=caf%C3%A9;gene_id=\xFF",
        ];
        let gff = write_tabix_fixture(dir.path(), "promoters.gff.gz", lines, gff_index());
        let sources = open(&format!(
            "file={},format=gff,gff_type=gencode_promoter",
            gff.display()
        ));
        let mut batch = vec![
            variant("21", 100, 100, 1),
            variant("21", 205, 205, 2),
            variant("21", 305, 305, 3),
        ];
        // The batch is matched on the rayon pool, so the warnings are collected
        // on a pool of the one thread the collector is installed on.
        let warnings = Arc::new(Warnings::default());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        pool.install(|| {
            tracing::subscriber::with_default(Arc::clone(&warnings), || {
                sources.run_batch(&mut batch)
            })
        })
        .unwrap();

        let window = &only(&batch[0]).records[0];
        assert_eq!(window.name, "21:90-110");
        assert_eq!(
            text_fields(window),
            [
                ("type", "promoter_window"),
                ("feature_id", "caf\u{280}"),
                ("associated_gene", "caf\u{280}"),
            ]
        );
        assert_eq!(
            text_fields(&only(&batch[1]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", "caf\u{FFFD}A"),
                ("associated_gene", "\u{3000}"),
            ]
        );
        assert_eq!(
            text_fields(&only(&batch[2]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", "caf\u{E9}"),
                ("associated_gene", "\u{FFFD}"),
            ]
        );
        let ignored = |feature_type: &str| {
            format!(
                "Ignoring '{feature_type}' feature_type from {} GFF/GTF file. \
                 This feature_type is not supported in VEP.",
                gff.display()
            )
        };
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                ignored("caf\u{FFFD}"),
                ignored("caf\u{FFFD}"),
                ignored("caf\u{280}"),
                ignored("mix\u{280}"),
                ignored("mix\u{280}"),
            ]
        );
    }

    /// `GFF3.pm:190-198` through Perl's two-variable `split`: the value is the
    /// text between a pair's first and second `=`, empty for `ID=`, absent for
    /// a pair without `=`; a later pair replaces an earlier one with its key.
    #[test]
    fn attributes_split_as_perl_splits_them() {
        assert_eq!(attribute("ID=PW1;gene_id=ENSG1", "gene_id"), Some("ENSG1"));
        assert_eq!(attribute("ID=a=b;x=1", "ID"), Some("a"));
        assert_eq!(attribute("ID==x", "ID"), Some(""));
        assert_eq!(attribute("ID=", "ID"), Some(""));
        assert_eq!(attribute("ID", "ID"), None);
        assert_eq!(attribute("ID=A;ID=B", "ID"), Some("B"));
        assert_eq!(attribute("ID=A;ID", "ID"), None);
        assert_eq!(attribute("Name=PW1", "ID"), None);
        assert_eq!(attribute("", "ID"), None);
        assert_eq!(seqname("chr21"), "21");
        assert_eq!(seqname("CHR21"), "CHR21");
        assert_eq!(seqname("21"), "21");
    }

    /// The attributes cell keeps the column's bytes for the decoder: a UTF-8
    /// column is borrowed as written, escapes included; each byte that is not
    /// UTF-8 is spelled as its escape, so `decode_attribute` over a value of
    /// the cell gives the bytes the file holds, read as text.
    #[test]
    fn the_attributes_cell_spells_a_byte_that_is_not_utf8_as_its_escape() {
        assert!(matches!(
            attributes_text(b"ID=PW1;gene_id=ENSG%3B1"),
            Cow::Borrowed("ID=PW1;gene_id=ENSG%3B1")
        ));
        assert_eq!(
            attributes_text(b"ID=caf\xCA%80;gene_id=caf%CA\x80"),
            "ID=caf%CA%80;gene_id=caf%CA%80"
        );
        assert_eq!(
            attributes_text(b"ID=\xEF\xBC%81;x=\xFF"),
            "ID=%EF%BC%81;x=%FF"
        );
        assert_eq!(attributes_text(b"ID=caf\xC3\xA9"), "ID=caf\u{E9}");
        let cell = attributes_text(b"ID=\xE3\x80%80;gene_id=caf\xFF%41");
        assert_eq!(
            decode_attribute(attribute(&cell, "ID").unwrap().as_bytes()),
            "\u{3000}"
        );
        assert_eq!(
            decode_attribute(attribute(&cell, "gene_id").unwrap().as_bytes()),
            "caf\u{FFFD}A"
        );
    }

    /// Ensembl VEP 116.2 stores the `ID` attribute as the window's `feature_id`
    /// (`File/GFF.pm:198`), so a record with no `ID` pair carries an undefined
    /// value, where `ID=` carries the empty string its two-variable `split`
    /// keeps (`GFF3.pm:190-198`): the text formats print both as the empty
    /// entry, the JSON format `null` and `""`.
    #[test]
    fn a_window_without_an_id_attribute_leaves_feature_id_undefined() {
        let dir = tempfile::tempdir().unwrap();
        let gff = write_tabix_fixture(
            dir.path(),
            "promoters.gff.gz",
            &[
                "21\t.\tpromoter_window\t90\t110\t.\t+\t.\tgene_id=ENSG00000000003;Name=PW3",
                "21\t.\tpromoter_window\t200\t210\t.\t+\t.\tID=;gene_id=ENSG00000000004",
            ],
            gff_index(),
        );
        let sources = open(&format!(
            "file={},format=gff,gff_type=gencode_promoter",
            gff.display()
        ));
        let mut batch = vec![variant("21", 100, 100, 1), variant("21", 205, 205, 2)];
        sources.run_batch(&mut batch).unwrap();
        let without_id = &only(&batch[0]).records[0];
        assert_eq!(without_id.fields["feature_id"], CustomFieldValue::Undefined);
        let empty_id = &only(&batch[1]).records[0];
        assert_eq!(
            empty_id.fields["feature_id"],
            CustomFieldValue::Text(String::new())
        );
    }

    /// Every `%XX` escape is decoded in `ID` and `gene_id`
    /// ([`decode_attribute`]), the six GFF3 reserves and any other, in either
    /// letter case, so `ID=PROM%2C1;gene_id=ENSG%3B1` carries `feature_id`
    /// `PROM,1` and `associated_gene` `ENSG;1`, and `ID=PROM%2c1%3E;gene_id=
    /// ENSG%3b%30` carries `PROM,1>` and `ENSG;0`, where Ensembl VEP 116.2
    /// prints each as written, the raw hash of `GFF3.pm:190-198`
    /// (`GFF.pm:198-201`). The split comes first, so `%3B` and `%3D` never cut
    /// a pair; the escapes decode in one pass, so `%2520%30` is `%200`; a
    /// `gene_id` is tested for Perl's truth as written, so `%30` is set where
    /// `0` is `.`.
    #[test]
    fn reserved_escapes_in_the_identifiers_are_decoded() {
        let dir = tempfile::tempdir().unwrap();
        let gff = write_tabix_fixture(
            dir.path(),
            "promoters.gff.gz",
            &[
                "21\t.\tpromoter_window\t90\t110\t.\t+\t.\tID=PROM%2C1;gene_id=ENSG%3B1",
                "21\t.\tpromoter_window\t200\t210\t.\t+\t.\tID=a%3Db%3Bc%3D1;gene_id=%2520%30",
                "21\t.\tpromoter_window\t300\t310\t.\t+\t.\tID=%20;gene_id=0",
                "21\t.\tpromoter_window\t400\t410\t.\t+\t.\tID=PROM%2c1%3E;gene_id=ENSG%3b%30",
            ],
            gff_index(),
        );
        let sources = open(&format!(
            "file={},format=gff,gff_type=gencode_promoter",
            gff.display()
        ));
        let mut batch = vec![
            variant("21", 100, 100, 1),
            variant("21", 205, 205, 2),
            variant("21", 305, 305, 3),
            variant("21", 405, 405, 4),
        ];
        sources.run_batch(&mut batch).unwrap();
        assert_eq!(
            text_fields(&only(&batch[0]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", "PROM,1"),
                ("associated_gene", "ENSG;1"),
            ]
        );
        assert_eq!(
            text_fields(&only(&batch[1]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", "a=b;c=1"),
                ("associated_gene", "%200"),
            ]
        );
        assert_eq!(
            text_fields(&only(&batch[2]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", " "),
                ("associated_gene", "."),
            ]
        );
        assert_eq!(
            text_fields(&only(&batch[3]).records[0]),
            [
                ("type", "promoter_window"),
                ("feature_id", "PROM,1>"),
                ("associated_gene", "ENSG;0"),
            ]
        );
    }
}
