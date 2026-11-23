// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The known variants of a sequence no variation VCF holds a record on, read
//! from the `variation_feature` table of the variation database.
//!
//! A species other than `homo_sapiens` publishes one genome-wide variation
//! VCF, and that file can hold no record on a sequence the database holds
//! known variants on (a mitochondrion, in release 116). `variation_feature`
//! is the table Ensembl VEP's dumper reads a record from
//! (`AnnotationSource/Database/Variation.pm:167-179`), so a row of it gives
//! every column a VCF line gives, in the database's own form:
//!
//! - the name, the start, the end and the strand as the row holds them;
//! - the allele string, each run of white space in it an underscore (`:188`);
//! - `somatic`;
//! - the clinical significance, each space of it an underscore (`:175`);
//! - `phenotype_or_disease`, 1 where the row's `evidence_attribs` hold the
//!   identifier of the `attrib` row whose value is `Phenotype_or_Disease`
//!   (`:176, 391-409`).
//!
//! The row's sequence is the name of its `seq_region_id` in the variation
//! database's `seq_region` table. Only the rows on display are read, the rows
//! a variation VCF holds: the export leaves the others out, and a failed
//! variant among them comes from the failed-variants GVF, as on any other
//! sequence.
//!
//! The table is read once as a stream; the rows of the sequences asked for
//! are sorted on disk by sequence and start, so neither the table nor the
//! records of a sequence are held in memory.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use tracing::info;

use super::{AllVarsWriter, ClinVarAttribs, GenomeWideInputs, ThousandGenomes, VariationRecord};
use crate::count::counted;
use crate::dump_rows::{read_rows, unsigned, DumpRow};
use crate::gvf::GvfFeature;
use crate::sorted_runs::{Merge, RunWriter};
use crate::variant_cells::VariantCells;

/// The bytes of rows buffered before a sorted run is written.
const RUN_BUFFER_BYTES: usize = 64 << 20;

/// The `attrib` value of the evidence the phenotype-or-disease flag is read
/// from.
const PHENOTYPE_EVIDENCE: &str = "Phenotype_or_Disease";

/// The dumps of the tables the known variants of a sequence are read from
/// where no variation VCF holds a record on it.
#[derive(Debug, Clone)]
pub struct FeatureTables {
    /// The dump of `variation_feature`.
    pub variation_feature: PathBuf,
    /// The dump of the variation database's `seq_region`, which names the
    /// sequence of each row.
    pub seq_region: PathBuf,
    /// The dump of `attrib`, which names the evidence the
    /// phenotype-or-disease flag is read from; without it no record carries
    /// the flag.
    pub attrib: Option<PathBuf>,
    /// The directory the kept rows are sorted under, created when it does
    /// not exist; the files written there are removed before the call
    /// returns.
    pub scratch: PathBuf,
}

/// Writes the variation file of each of `sequences` from the rows `tables`
/// hold on it, through an [`AllVarsWriter`] with `columns`: one record a row
/// on display, in start order and, at one start, in name order, a row
/// written twice kept once. Each record takes its ClinVar columns from
/// `clinvar`, its 1000 Genomes columns from the frequency file of `inputs`
/// and its `pubmed`, `var_synonyms` and `failed` cells from `cells` by its
/// name, as a record read from a VCF does. The features of the
/// failed-variants GVF on the sequence join its records: one whose name a
/// row carries at its start flags that record, and any other is a record of
/// its own.
///
/// A sequence with no row and no such feature gets no file. Returns each
/// sequence with the number of records written; on an error the file under
/// way is removed.
pub fn write_sequences_from_tables(
    tables: &FeatureTables,
    sequences: &[String],
    inputs: &GenomeWideInputs,
    clinvar: Option<&ClinVarAttribs>,
    cells: Option<&VariantCells>,
    columns: &[String],
    output_dir: &Path,
) -> Result<Vec<(String, usize)>> {
    let from = Sources {
        inputs,
        clinvar,
        cells,
    };
    write_sequences(
        tables,
        sequences,
        &from,
        columns,
        output_dir,
        RUN_BUFFER_BYTES,
    )
}

/// What a record takes beside its row.
struct Sources<'a> {
    inputs: &'a GenomeWideInputs,
    clinvar: Option<&'a ClinVarAttribs>,
    cells: Option<&'a VariantCells>,
}

fn write_sequences(
    tables: &FeatureTables,
    sequences: &[String],
    from: &Sources<'_>,
    columns: &[String],
    output_dir: &Path,
    run_buffer: usize,
) -> Result<Vec<(String, usize)>> {
    let runs = read_rows_of(tables, sequences, run_buffer)?;
    let mut rows = runs.merge()?;
    let mut written = Vec::with_capacity(sequences.len());
    for (index, sequence) in sequences.iter().enumerate() {
        let failed = from
            .cells
            .map_or(&[][..], |cells| cells.failed_features(sequence));
        let holds_a_row = rows.head().is_some_and(|((of, _), _)| of == index as u64);
        if !holds_a_row && failed.is_empty() {
            written.push((sequence.clone(), 0));
            continue;
        }
        let mut writer = AllVarsWriter::create(output_dir, sequence, columns)?;
        let records = write_sequence(&mut rows, index as u64, sequence, failed, from, &mut writer);
        match records {
            Ok(()) => written.push((sequence.clone(), writer.finish()?)),
            Err(error) => {
                writer.discard();
                return Err(error);
            }
        }
    }
    Ok(written)
}

/// Writes the records of the sequence at `index` of the sorted `rows`, with
/// the features of `failed` that no row carries, to `writer`.
fn write_sequence(
    rows: &mut Merge,
    index: u64,
    sequence: &str,
    failed: &[GvfFeature],
    from: &Sources<'_>,
    writer: &mut AllVarsWriter,
) -> Result<()> {
    let frequencies = match &from.inputs.frequencies {
        Some(path) => Some(ThousandGenomes::load(path, sequence)?),
        None => None,
    };
    let mut write = |mut record: VariationRecord| -> Result<()> {
        if let Some(clinvar) = from.clinvar {
            clinvar.apply(&mut record, sequence);
        }
        if let Some(cells) = from.cells {
            cells.apply(&mut record)?;
        }
        writer.write(&record)
    };
    let mut failed = failed.iter().peekable();
    let mut added = 0u64;
    loop {
        let start = match rows.head() {
            Some(((of, start), _)) if of == index => Some(start),
            _ => None,
        };
        // The records of one start: its rows, then the features no row of
        // that start carries, in name order.
        let mut at_start: Vec<VariationRecord> = Vec::new();
        while let Some(((of, at), row)) = rows.head() {
            if of != index || Some(at) != start {
                break;
            }
            let mut record = decoded(at, row)?;
            // The dumper joins frequencies onto non-somatic variants only
            // (`Dumper/Variation.pm:280`), matching a record's own alleles
            // at its own position.
            if let (Some(frequencies), 0) = (&frequencies, record.somatic) {
                let alleles = record.allele_string.clone();
                if let Some((reference, alternates)) = alleles.split_once('/') {
                    let alternates: Vec<&str> = alternates.split('/').collect();
                    let position = record.start.min(record.end);
                    frequencies.apply(&mut record, position, reference, &alternates);
                }
            }
            at_start.push(record);
            rows.advance()?;
        }
        let place = |record: &VariationRecord| {
            (
                record.variation_name.clone(),
                record.end,
                record.allele_string.clone(),
            )
        };
        at_start.sort_by_key(place);
        at_start.dedup_by_key(|record| place(record));
        let rows_at_start = at_start.len();
        while let Some(feature) = failed.peek() {
            let before = start.is_none_or(|start| feature.record.start < start);
            if !before && Some(feature.record.start) != start {
                break;
            }
            let carried = !before
                && at_start[..rows_at_start]
                    .iter()
                    .any(|row| row.variation_name == feature.record.variation_name);
            if before {
                added += 1;
                write(super::feature_record(feature, frequencies.as_ref()))?;
            } else if !carried {
                added += 1;
                at_start.push(super::feature_record(feature, frequencies.as_ref()));
            }
            failed.next();
        }
        at_start.sort_by(|a, b| a.variation_name.cmp(&b.variation_name));
        for record in at_start {
            write(record)?;
        }
        if start.is_none() {
            break;
        }
    }
    if let Some(cells) = from.cells {
        cells.count_added(added);
    }
    Ok(())
}

/// Reads the rows of `tables` on `sequences` that are on display and spills
/// each as a record, sorted by the sequence's place in `sequences` and the
/// row's start.
fn read_rows_of(
    tables: &FeatureTables,
    sequences: &[String],
    run_buffer: usize,
) -> Result<crate::sorted_runs::Runs> {
    std::fs::create_dir_all(&tables.scratch)
        .with_context(|| format!("Failed to create {}", tables.scratch.display()))?;

    // Every identifier a name has: the table lists a sequence once for each
    // coordinate system that holds it.
    let mut wanted: HashMap<u32, u64> = HashMap::new();
    let seq_region = &tables.seq_region;
    read_rows(seq_region, 3, &[0, 1], |row| {
        let id = number(row, 0, seq_region, "seq_region_id")?;
        let name = row.value(1).unwrap_or_default();
        if let Some(index) = sequences.iter().position(|s| s.as_bytes() == name) {
            wanted.insert(id as u32, index as u64);
        }
        Ok(())
    })?;

    // `SELECT attrib_id FROM attrib WHERE value = 'Phenotype_or_Disease'`
    // (`Database/Variation.pm:395-397`), the first row of that value.
    let mut phenotype_evidence: Option<Vec<u8>> = None;
    if let Some(attrib) = &tables.attrib {
        read_rows(attrib, 3, &[0, 2], |row| {
            let named = row
                .value(1)
                .is_some_and(|value| value.eq_ignore_ascii_case(PHENOTYPE_EVIDENCE.as_bytes()));
            if named && phenotype_evidence.is_none() {
                phenotype_evidence = row.value(0).map(<[u8]>::to_vec);
            }
            Ok(())
        })?;
    }

    let table = &tables.variation_feature;
    info!(
        "Reading the known variants of {} from {}",
        counted(sequences.len(), "sequence"),
        table.display()
    );
    let mut runs = RunWriter::new(&tables.scratch, "features", run_buffer);
    let mut payload = Vec::new();
    let read = read_rows(table, 23, &[1, 2, 3, 4, 6, 8, 15, 20, 21, 22], |row| {
        let region = number(row, 0, table, "seq_region_id")?;
        let Some(index) = wanted.get(&(region as u32)) else {
            return Ok(());
        };
        let is = |nth: usize, value: &str| row.value(nth) == Some(value.as_bytes());
        if !is(9, "1") {
            return Ok(());
        }
        let start = number(row, 1, table, "seq_region_start")?;
        let text = |nth: usize| String::from_utf8_lossy(row.value(nth).unwrap_or_default());
        let evidence = row.value(7).unwrap_or_default();
        let record = VariationRecord {
            variation_name: text(5).into_owned(),
            somatic: u8::from(is(6, "1")),
            start,
            end: number(row, 2, table, "seq_region_end")?,
            strand: if is(3, "-1") { -1 } else { 1 },
            allele_string: white_space_as_underscores(&text(4)),
            clin_sig: row
                .value(8)
                .filter(|terms| !terms.is_empty())
                .map(|terms| String::from_utf8_lossy(terms).replace(' ', "_")),
            phenotype_or_disease: u8::from(
                phenotype_evidence
                    .as_deref()
                    .is_some_and(|id| evidence.split(|byte| *byte == b',').any(|held| held == id)),
            ),
            ..Default::default()
        };
        encode(&record, &mut payload);
        runs.push((*index, start), &payload)
    })?;
    let runs = runs.finish()?;
    info!(
        "Known variants read: {} of {} on those sequences and on display",
        runs.records,
        counted(read, "row")
    );
    Ok(runs)
}

/// `text` with each run of white space an underscore, as the dumper's
/// `s/\s+/_/g` leaves an allele string (`Database/Variation.pm:188`).
fn white_space_as_underscores(text: &str) -> String {
    let mut written = String::with_capacity(text.len());
    let mut in_a_run = false;
    for character in text.chars() {
        let space = matches!(character, ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c');
        if !space {
            written.push(character);
        } else if !in_a_run {
            written.push('_');
        }
        in_a_run = space;
    }
    written
}

/// The number in the `nth` kept column of `row`, or an error naming the
/// file, the row and the column.
fn number(row: &DumpRow<'_>, nth: usize, path: &Path, column: &str) -> Result<u64> {
    row.value(nth)
        .and_then(unsigned)
        .map(u64::from)
        .with_context(|| {
            format!(
                "{} row {} holds no number in {column}: {}",
                path.display(),
                row.number,
                match row.value(nth) {
                    Some(value) => format!("{:?}", String::from_utf8_lossy(value)),
                    None => "NULL".to_string(),
                }
            )
        })
}

/// `record` as a sorted row carries it beside its start: the end, the
/// strand, the two flags and whether it has a clinical significance, then
/// the name, the allele string and the clinical significance, each behind
/// its length.
fn encode(record: &VariationRecord, payload: &mut Vec<u8>) {
    payload.clear();
    payload.extend_from_slice(&record.end.to_le_bytes());
    payload.push(record.strand as u8);
    payload.push(record.somatic);
    payload.push(record.phenotype_or_disease);
    payload.push(u8::from(record.clin_sig.is_some()));
    for text in [
        record.variation_name.as_str(),
        record.allele_string.as_str(),
        record.clin_sig.as_deref().unwrap_or_default(),
    ] {
        payload.extend_from_slice(&(text.len() as u32).to_le_bytes());
        payload.extend_from_slice(text.as_bytes());
    }
}

/// The record [`encode`] wrote, at `start`.
fn decoded(start: u64, payload: &[u8]) -> Result<VariationRecord> {
    let cut = || anyhow::anyhow!("a sorted row of variation_feature is cut short");
    let (fixed, mut texts) = payload.split_at_checked(12).ok_or_else(cut)?;
    let mut text = || -> Result<String> {
        let (length, rest) = texts.split_at_checked(4).ok_or_else(cut)?;
        let length = u32::from_le_bytes(length.try_into().expect("4 bytes")) as usize;
        let (text, rest) = rest.split_at_checked(length).ok_or_else(cut)?;
        texts = rest;
        Ok(String::from_utf8_lossy(text).into_owned())
    };
    let variation_name = text()?;
    let allele_string = text()?;
    let clin_sig = text()?;
    Ok(VariationRecord {
        variation_name,
        somatic: fixed[9],
        start,
        end: u64::from_le_bytes(fixed[..8].try_into().expect("8 bytes")),
        strand: fixed[8] as i8,
        allele_string,
        clin_sig: (fixed[11] != 0).then_some(clin_sig),
        phenotype_or_disease: fixed[10],
        ..Default::default()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant_cells::VariantCellInputs;
    use crate::variation::variation_columns;
    use noodles::bgzf;
    use std::io::Read;

    /// A made-up row of `variation_feature`: the values a record is read
    /// from, among the schema's 23 columns.
    #[derive(Clone, Copy)]
    struct Row<'a> {
        region: u32,
        start: &'a str,
        end: u64,
        strand: &'a str,
        alleles: &'a str,
        name: &'a str,
        somatic: &'a str,
        evidence: &'a str,
        significance: &'a str,
        display: &'a str,
    }

    /// A substitution on display on region 11, with no evidence and no
    /// clinical significance.
    const SNV: Row<'static> = Row {
        region: 11,
        start: "500",
        end: 500,
        strand: "1",
        alleles: "A/G",
        name: "mv500",
        somatic: "0",
        evidence: "\\N",
        significance: "\\N",
        display: "1",
    };

    fn table(rows: &[Row<'_>]) -> String {
        rows.iter()
            .enumerate()
            .map(|(id, row)| {
                let Row {
                    region,
                    start,
                    end,
                    strand,
                    alleles,
                    name,
                    somatic,
                    evidence,
                    significance,
                    display,
                } = row;
                format!(
                    "{id}\t{region}\t{start}\t{end}\t{strand}\t{id}\t{alleles}\t\\N\t{name}\t1\t\\N\t\
                     1\tintergenic_variant\t\t2\t{somatic}\t\\N\t\\N\t\\N\t\\N\t{evidence}\t\
                     {significance}\t{display}\n"
                )
            })
            .collect()
    }

    /// The made-up tables in a directory: `seq_region` names `seqA` under
    /// two identifiers, 11 and 31, `seqB` under 12 and `seqC` under 13, and
    /// `attrib` gives the phenotype evidence the identifier 418.
    fn tables(dir: &Path, rows: &[Row<'_>]) -> FeatureTables {
        let write = |name: &str, text: &str| {
            let path = dir.join(name);
            std::fs::write(&path, text).unwrap();
            path
        };
        FeatureTables {
            variation_feature: write("variation_feature.txt", &table(rows)),
            seq_region: write(
                "seq_region.txt",
                "11\tseqA\t2\n12\tseqB\t2\n31\tseqA\t3\n13\tseqC\t2\n",
            ),
            attrib: Some(write(
                "attrib.txt",
                "367\t497\tMultiple_observations\n418\t497\tPhenotype_or_Disease\n",
            )),
            scratch: dir.join("scratch"),
        }
    }

    /// The lines of the variation file written for `sequence` under `out`,
    /// each without its chromosome; `None` where no file was written.
    fn lines(out: &Path, sequence: &str) -> Option<Vec<String>> {
        let path = out.join("variations").join(sequence).join("all_vars.gz");
        let file = std::fs::File::open(path).ok()?;
        let mut text = String::new();
        bgzf::io::Reader::new(file)
            .read_to_string(&mut text)
            .unwrap();
        Some(
            text.lines()
                .map(|line| {
                    line.strip_prefix(&format!("{sequence}\t"))
                        .unwrap()
                        .to_string()
                })
                .collect(),
        )
    }

    fn write(
        tables: &FeatureTables,
        sequences: &[&str],
        from: &Sources<'_>,
        out: &Path,
        run_buffer: usize,
    ) -> Vec<(String, usize)> {
        let sequences: Vec<String> = sequences.iter().map(|s| s.to_string()).collect();
        let columns = variation_columns(from.inputs.frequencies.is_some());
        write_sequences(tables, &sequences, from, &columns, out, run_buffer).unwrap()
    }

    const NO_INPUTS: GenomeWideInputs = GenomeWideInputs {
        frequencies: None,
        somatic: None,
        germline: Vec::new(),
    };

    const ROWS_ALONE: Sources<'static> = Sources {
        inputs: &NO_INPUTS,
        clinvar: None,
        cells: None,
    };

    /// A row on display reads as the record the dumper reads from it: the
    /// name, the span, the strand and `somatic` as the row holds them, an
    /// insertion with its start past its end; the allele string with each
    /// run of white space an underscore, empty where the row holds none;
    /// (a tab among it written behind a backslash, as a dump writes one);
    /// the clinical significance with each space an underscore; and the
    /// phenotype-or-disease flag where the row's evidence holds the
    /// identifier `attrib` gives that evidence, as one of its items and not
    /// inside another. A row that is not on display, a row whose display is
    /// NULL and a row on another sequence are passed over. The rows of both
    /// identifiers of a sequence are its rows. The records stand in start
    /// order whatever the table's order, those of one start in name order,
    /// and a row written twice is one record, another row between the two
    /// or not.
    #[test]
    fn a_row_on_display_reads_as_the_record_the_dumper_reads_from_it() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            Row {
                start: "700",
                end: 700,
                name: "mv700",
                alleles: "G/A",
                evidence: "4180,367",
                ..SNV
            },
            Row {
                evidence: "367,418",
                significance: "likely pathogenic,pathogenic",
                ..SNV
            },
            Row {
                region: 31,
                start: "300",
                end: 299,
                name: "mv300",
                alleles: "-/TT",
                evidence: "367",
                ..SNV
            },
            Row {
                start: "400",
                end: 402,
                strand: "-1",
                name: "mv400",
                alleles: "CAT/-",
                somatic: "1",
                ..SNV
            },
            Row {
                start: "450",
                end: 459,
                name: "mv450",
                alleles: "(LARGE \\\t DELETION)/A",
                ..SNV
            },
            Row {
                start: "460",
                end: 460,
                name: "mv460",
                alleles: "\\N",
                ..SNV
            },
            Row {
                start: "470",
                end: 470,
                name: "mv470",
                display: "0",
                ..SNV
            },
            Row {
                start: "480",
                end: 480,
                name: "mv480",
                display: "\\N",
                ..SNV
            },
            Row {
                region: 12,
                start: "100",
                end: 100,
                name: "mv100",
                ..SNV
            },
            Row {
                start: "600",
                end: 600,
                name: "mva",
                alleles: "C/T",
                ..SNV
            },
            Row {
                start: "600",
                end: 600,
                name: "mvb",
                alleles: "C/T",
                ..SNV
            },
            Row {
                start: "600",
                end: 600,
                name: "mva",
                alleles: "C/T",
                ..SNV
            },
        ];
        let out = dir.path().join("cache");
        let written = write(
            &tables(dir.path(), &rows),
            &["seqA"],
            &ROWS_ALONE,
            &out,
            1 << 20,
        );
        assert_eq!(written, [("seqA".to_string(), 8)]);
        assert_eq!(
            lines(&out, "seqA").unwrap(),
            [
                "mv300\t.\t.\t300\t299\t-/TT\t.\t.\t.\t.\t.\t.\t.\t.",
                "mv400\t.\t1\t400\t402\tCAT/-\t-1\t.\t.\t.\t.\t.\t.\t.",
                "mv450\t.\t.\t450\t459\t(LARGE_DELETION)/A\t.\t.\t.\t.\t.\t.\t.\t.",
                "mv460\t.\t.\t460\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                "mv500\t.\t.\t500\t.\tA/G\t.\tlikely_pathogenic,pathogenic\t1\t.\t.\t.\t.\t.",
                "mva\t.\t.\t600\t.\tC/T\t.\t.\t.\t.\t.\t.\t.\t.",
                "mvb\t.\t.\t600\t.\tC/T\t.\t.\t.\t.\t.\t.\t.\t.",
                "mv700\t.\t.\t700\t.\tG/A\t.\t.\t.\t.\t.\t.\t.\t.",
            ]
        );
        assert!(out.join("variations/seqA/all_vars.gz.csi").is_file());
        // The sorted rows are removed with the call.
        assert_eq!(
            std::fs::read_dir(dir.path().join("scratch"))
                .unwrap()
                .count(),
            0
        );
    }

    /// Each sequence asked for is written from its own rows, in whatever
    /// order the sequences are asked for, and one with no row gets no file.
    /// A run buffer that holds two rows spills the rows as several runs,
    /// which merge into the files one run gives. Without `attrib` no record
    /// carries the phenotype-or-disease flag.
    #[test]
    fn each_sequence_is_written_from_its_own_rows_however_they_are_spilled() {
        let dir = tempfile::tempdir().unwrap();
        let mut rows = Vec::new();
        let names: Vec<String> = (0..12).map(|n| format!("mv{n:02}")).collect();
        let starts: Vec<String> = (0..12).map(|n| (900 - 50 * n).to_string()).collect();
        for n in 0..12 {
            rows.push(Row {
                region: [11, 12, 31][n % 3],
                start: &starts[n],
                end: 900 - 50 * n as u64,
                name: &names[n],
                evidence: "418",
                ..SNV
            });
        }
        let tables = tables(dir.path(), &rows);
        let asked = ["seqC", "seqB", "seqA"];
        let whole = dir.path().join("whole");
        assert_eq!(
            write(&tables, &asked, &ROWS_ALONE, &whole, 1 << 20),
            [
                ("seqC".to_string(), 0),
                ("seqB".to_string(), 4),
                ("seqA".to_string(), 8)
            ]
        );
        assert!(lines(&whole, "seqC").is_none());
        assert!(!whole.join("variations/seqC").exists());
        let of_b = lines(&whole, "seqB").unwrap();
        assert_eq!(
            of_b,
            [
                "mv10\t.\t.\t400\t.\tA/G\t.\t.\t1\t.\t.\t.\t.\t.",
                "mv07\t.\t.\t550\t.\tA/G\t.\t.\t1\t.\t.\t.\t.\t.",
                "mv04\t.\t.\t700\t.\tA/G\t.\t.\t1\t.\t.\t.\t.\t.",
                "mv01\t.\t.\t850\t.\tA/G\t.\t.\t1\t.\t.\t.\t.\t.",
            ]
        );

        let runs = read_rows_of(&tables, &asked.map(String::from), 120).unwrap();
        assert_eq!((runs.count(), runs.records), (6, 12));
        drop(runs);
        let spilled = dir.path().join("spilled");
        write(&tables, &asked, &ROWS_ALONE, &spilled, 120);
        for sequence in ["seqA", "seqB"] {
            assert_eq!(
                lines(&spilled, sequence),
                lines(&whole, sequence),
                "{sequence}"
            );
        }

        let no_attrib = FeatureTables {
            attrib: None,
            ..tables
        };
        let unflagged = dir.path().join("unflagged");
        write(&no_attrib, &["seqB"], &ROWS_ALONE, &unflagged, 1 << 20);
        assert_eq!(
            lines(&unflagged, "seqB").unwrap(),
            of_b.iter()
                .map(|line| line.replace("\t1\t", "\t.\t"))
                .collect::<Vec<_>>()
        );
    }

    /// The records take what a VCF record takes beside its line. The
    /// features of the failed-variants GVF on the sequence join them: a
    /// feature a row carries by name at its start flags that record; a
    /// feature before the first row, between two rows, at a row's start
    /// under another name and after the last row is a record of its own, in
    /// start order and name order; a sequence that has features and no row
    /// is written from its features. A record takes the 1000 Genomes
    /// frequencies of its alleles, matched at its own position.
    #[test]
    fn records_of_the_table_take_their_cells_and_the_failed_features_join_them() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            Row {
                start: "200",
                end: 200,
                name: "mv200",
                ..SNV
            },
            SNV,
        ];
        let feature = |sequence: &str, at: u64, name: &str| {
            format!(
                "{sequence}\tmadeup\tSNV\t{at}\t{at}\t.\t+\t.\t\
                 ID=1;Dbxref=madeup_1:{name};Reference_seq=C;Variant_seq=T\n"
            )
        };
        let gvf = dir.path().join("made_up_failed.gvf");
        std::fs::write(
            &gvf,
            [
                feature("seqA", 100, "mv100"),
                feature("seqA", 200, "mv200"),
                feature("seqA", 200, "mv199"),
                feature("seqA", 350, "mv350"),
                feature("seqA", 500, "mv900"),
                feature("seqA", 900, "mv900"),
                feature("seqC", 40, "mv40"),
            ]
            .concat(),
        )
        .unwrap();
        let frequencies = dir.path().join("made_up_frequencies.vcf");
        std::fs::write(
            &frequencies,
            "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             seqA\t500\tf\tA\tG\t.\t.\tAFR=0.5;AMR=0.25;EAS=0.125;EUR=0.75;SAS=0.0625\n",
        )
        .unwrap();
        let cells = VariantCells::read(
            &VariantCellInputs {
                failed_gvf: Some(gvf),
                scratch: dir.path().join("cells"),
                ..Default::default()
            },
            &["seqA".to_string(), "seqC".to_string()],
        )
        .unwrap();
        let inputs = GenomeWideInputs {
            frequencies: Some(frequencies),
            ..Default::default()
        };
        let from = Sources {
            inputs: &inputs,
            clinvar: None,
            cells: Some(&cells),
        };
        let out = dir.path().join("cache");
        let written = write(
            &tables(dir.path(), &rows),
            &["seqA", "seqC"],
            &from,
            &out,
            1 << 20,
        );
        assert_eq!(written, [("seqA".to_string(), 7), ("seqC".to_string(), 1)]);
        let read = |sequence: &str| -> Vec<(String, String, String, String, String)> {
            lines(&out, sequence)
                .unwrap()
                .iter()
                .map(|line| {
                    let cells: Vec<&str> = line.split('\t').collect();
                    // name, failed, start, allele string, AFR.
                    (
                        cells[0].to_string(),
                        cells[1].to_string(),
                        cells[3].to_string(),
                        cells[5].to_string(),
                        cells[15].to_string(),
                    )
                })
                .collect()
        };
        let record = |name: &str, failed: &str, start: &str, alleles: &str, afr: &str| {
            (
                name.to_string(),
                failed.to_string(),
                start.to_string(),
                alleles.to_string(),
                afr.to_string(),
            )
        };
        assert_eq!(
            read("seqA"),
            [
                record("mv100", "1", "100", "C/T", "."),
                record("mv199", "1", "200", "C/T", "."),
                record("mv200", "1", "200", "A/G", "."),
                record("mv350", "1", "350", "C/T", "."),
                record("mv500", ".", "500", "A/G", "G:0.5"),
                record("mv900", "1", "500", "C/T", "."),
                record("mv900", "1", "900", "C/T", "."),
            ]
        );
        assert_eq!(read("seqC"), [record("mv40", "1", "40", "C/T", ".")]);
        let facts = cells.facts();
        assert_eq!(
            (facts.records_flagged_failed, facts.records_added_from_gvf),
            (7, 6)
        );
    }

    /// A row whose start is no number refuses the write, naming the file,
    /// the row, the column and the value, and leaves no variation file; so
    /// does a row short of the schema's columns.
    #[test]
    fn a_row_that_does_not_read_refuses_the_write() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            SNV,
            Row {
                start: "soon",
                ..SNV
            },
        ];
        let tables = tables(dir.path(), &rows);
        let out = dir.path().join("cache");
        let columns = variation_columns(false);
        let sequences = ["seqA".to_string()];
        let refused = |tables: &FeatureTables| {
            let error = write_sequences(tables, &sequences, &ROWS_ALONE, &columns, &out, 1 << 20)
                .unwrap_err();
            assert!(!out.join("variations").exists());
            format!("{error:#}")
        };
        assert_eq!(
            refused(&tables),
            format!(
                "{} row 2 holds no number in seq_region_start: \"soon\"",
                tables.variation_feature.display()
            )
        );
        std::fs::write(&tables.variation_feature, "1\t11\t500\t500\t1\n").unwrap();
        assert_eq!(
            refused(&tables),
            format!(
                "{} row 1 has 5 columns, the schema has 23",
                tables.variation_feature.display()
            )
        );
    }

    /// An allele string keeps what is not white space and takes one
    /// underscore for each run of it, at its ends too.
    #[test]
    fn each_run_of_white_space_is_one_underscore() {
        assert_eq!(white_space_as_underscores("A/G"), "A/G");
        assert_eq!(white_space_as_underscores(" A \t\n G  "), "_A_G_");
        assert_eq!(white_space_as_underscores(""), "");
    }
}
