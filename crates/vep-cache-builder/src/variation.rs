// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Build the variation files of a native cache from Ensembl's variation VCFs.
//!
//! Each chromosome's records are read from Ensembl's per-chromosome germline
//! VCF, joined with the records of the genome-wide somatic, phenotype- and
//! clinically-associated VCFs that the per-chromosome file lacks, given their
//! 1000 Genomes frequencies from the frequency VCF and, with the variation
//! database dumps, their ClinVar allele and somatic columns; then written as
//! the file Ensembl VEP's dumper writes (`Pipeline/DumpVEP/Dumper/Variation.pm:
//! 98-135, 211-270`): `variations/<chr>/all_vars.gz`, tab-delimited bgzip with
//! the chromosome first and `.` for an empty cell, sorted by start and indexed
//! by a `.csi` on the start column alone, its columns the `variation_cols` of
//! `info.json`.

use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::csi::binning_index::index::header::format::CoordinateSystem;
use noodles::csi::binning_index::index::header::Format;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
use noodles::csi::binning_index::index::Header;
use noodles::csi::binning_index::Indexer;
use noodles::vcf;
use tracing::{debug, info, warn};
use vep_core::cache_info::FIXED_VARIATION_COLUMNS;
use vep_core::coordinate::is_cache_path_component;

use crate::core_dumps::read_table;

/// The allele string of a COSMIC variant, whose alleles the variation database
/// does not hold (`VariationFeature.pm:2355`).
const COSMIC_MUTATION: &str = "COSMIC_MUTATION";

/// The allele string of an HGMD-PUBLIC variant, likewise without alleles.
const HGMD_MUTATION: &str = "HGMD_MUTATION";

/// The 1000 Genomes phase 3 populations of the frequency VCF, in the order the
/// dumper lists their columns (`DumpVEP_conf.pm:150`), each with the number of
/// samples of the phase 3 panel: 2,504 in all.
const THOUSAND_GENOMES_POPULATIONS: [(&str, u64); 5] = [
    ("AFR", 661),
    ("AMR", 347),
    ("EAS", 504),
    ("EUR", 503),
    ("SAS", 489),
];

/// The frequency columns the 1000 Genomes file yields: the whole panel's
/// frequency under the name the dumper gives its unnamed population
/// (`Dumper/Variation.pm:191`), then the five continental populations.
pub const THOUSAND_GENOMES_COLUMNS: [&str; 6] = ["AF", "AFR", "AMR", "EAS", "EUR", "SAS"];

/// The chromosome and start columns of an `all_vars.gz` record, 0-based.
const CHR_COLUMN: usize = 0;
const START_COLUMN: usize = 4;

/// A variation record in the form the variation database holds it, the form
/// Ensembl VEP's dumper writes.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct VariationRecord {
    pub variation_name: String,
    pub failed: u8,
    pub somatic: u8,
    pub start: u64,
    /// `start - 1` for an insertion.
    pub end: u64,
    pub strand: i8,
    /// `/`-separated alleles, or a bare class name for a variant without
    /// alleles.
    pub allele_string: String,
    /// Clinical significance terms joined by `,`.
    pub clin_sig: Option<String>,
    pub phenotype_or_disease: u8,
    /// ClinVar significance per risk allele, `allele:term` pairs joined by `;`.
    pub clin_sig_allele: Option<String>,
    /// The ClinVar reference allele behind `clin_sig_allele`.
    pub clin_sig_ref_allele: Option<String>,
    /// ClinVar somatic classification, `phenotype(classification)` entries joined by `;`.
    pub clinical_impact: Option<String>,
    /// Population frequencies by column name, `allele:frequency` entries joined
    /// by `,`.
    pub frequencies: HashMap<String, String>,
}

/// The columns of the variation files: `chr`, the fixed columns and, when the
/// frequency VCF was read, the 1000 Genomes columns.
pub fn variation_columns(with_frequencies: bool) -> Vec<String> {
    let mut columns: Vec<String> = FIXED_VARIATION_COLUMNS
        .iter()
        .map(|c| c.to_string())
        .collect();
    if with_frequencies {
        columns.extend(THOUSAND_GENOMES_COLUMNS.iter().map(|c| c.to_string()));
    }
    columns
}

/// One germline ClinVar phenotype feature: the attributes
/// `PhenotypeFeatureAdaptor::get_clinsig_alleles_by_location` returns per feature.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct GermlineClinVar {
    risk_allele: Option<String>,
    clinvar_clin_sig: Option<String>,
    clinvar_ref_allele: Option<String>,
}

/// One somatic ClinVar phenotype feature: the attributes
/// `PhenotypeFeatureAdaptor::get_somatic_clin_impact_by_location` returns per feature.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct SomaticClinVar {
    phenotype: String,
    somatic_clin_sig: Option<String>,
    oncogenic_clin_sig: Option<String>,
}

/// The ClinVar phenotype features of a variation database.
///
/// Ensembl VEP's dumper joins `clin_sig_allele`, `clin_sig_ref_allele` and
/// `clinical_impact` onto each variation feature from the `phenotype_feature`
/// and `phenotype_feature_attrib` tables (`Database/Variation.pm:161-162,
/// 196-221`), looking them up by the feature's `chromosome:start-end` span; the
/// Ensembl variation VCFs carry none of the three, so this type reads the same
/// tables from the database's text dumps. The tables hold the database's span,
/// while a record read from the VCF export may hold the trimmed one, so the
/// germline features, which the dumper keeps only under the variant's own name,
/// are keyed by `chromosome:variation_name`, and the somatic features, which it
/// keeps whatever their name, by span, a record's database span being read from
/// its own somatic features when it has any.
#[derive(Debug, Default)]
pub struct ClinVarAttribs {
    germline: HashMap<String, Vec<GermlineClinVar>>,
    somatic: HashMap<String, Vec<SomaticClinVar>>,
    /// By `chromosome:variation_name`, the spans of the variant's own somatic
    /// features.
    somatic_spans: HashMap<String, Vec<String>>,
}

impl ClinVarAttribs {
    /// Read the ClinVar phenotype features from the variation database dumps in `dir`:
    /// `source`, `seq_region`, `attrib_type`, `phenotype`, `phenotype_feature` and
    /// `phenotype_feature_attrib`, each as `.txt.gz` or `.txt`, on the sequence regions
    /// named by `chromosomes` (every region when the list is empty).
    ///
    /// Keeps the features of the `ClinVar` source whatever their object type, as the
    /// adaptor's two location queries do: a germline feature (`DNA_type` `Germline` or
    /// NULL) when it carries a `risk_allele` or a `clinvar_ref_allele` attribute
    /// (`PhenotypeFeatureAdaptor.pm:1253-1254`), a somatic feature (`DNA_type`
    /// `Somatic`) with its phenotype description (`:1303, 1320`). Errors when a table is
    /// missing; the rows are read with the core dump reader, MySQL's escapes decoded.
    pub fn load(dir: &Path, chromosomes: &[String]) -> Result<Self> {
        let mut clinvar_source_id = None;
        for row in read_table(dir, "source", 8)? {
            let row = row?;
            if row.field(1) == "ClinVar" {
                clinvar_source_id = Some(row.field(0).to_string());
            }
        }
        let Some(clinvar_source_id) = clinvar_source_id else {
            info!("No ClinVar source in {}", dir.display());
            return Ok(Self::default());
        };

        let mut seq_region_names: HashMap<String, String> = HashMap::new();
        for row in read_table(dir, "seq_region", 3)? {
            let row = row?;
            if chromosomes.is_empty() || chromosomes.iter().any(|c| c == row.field(1)) {
                seq_region_names.insert(row.field(0).to_string(), row.field(1).to_string());
            }
        }

        let mut attrib_codes: HashMap<String, String> = HashMap::new();
        for row in read_table(dir, "attrib_type", 4)? {
            let row = row?;
            attrib_codes.insert(row.field(0).to_string(), row.field(1).to_string());
        }

        let mut phenotypes: HashMap<String, String> = HashMap::new();
        for row in read_table(dir, "phenotype", 5)? {
            let row = row?;
            phenotypes.insert(
                row.field(0).to_string(),
                row.optional(3).unwrap_or("").to_string(),
            );
        }

        struct Feature {
            span_key: String,
            name_key: String,
            phenotype_id: String,
            somatic: bool,
            risk_allele: Option<String>,
            clinvar_clin_sig: Option<String>,
            clinvar_ref_allele: Option<String>,
            somatic_clin_sig: Option<String>,
            oncogenic_clin_sig: Option<String>,
        }
        let mut features: HashMap<String, Feature> = HashMap::new();
        for row in read_table(dir, "phenotype_feature", 12)? {
            let row = row?;
            if row.field(2) != clinvar_source_id {
                continue;
            }
            let Some(chr) = seq_region_names.get(row.field(7)) else {
                continue;
            };
            features.insert(
                row.field(0).to_string(),
                Feature {
                    span_key: format!("{chr}:{}-{}", row.field(8), row.field(9)),
                    name_key: format!("{chr}:{}", row.field(5)),
                    phenotype_id: row.field(1).to_string(),
                    somatic: row.optional(11) == Some("Somatic"),
                    risk_allele: None,
                    clinvar_clin_sig: None,
                    clinvar_ref_allele: None,
                    somatic_clin_sig: None,
                    oncogenic_clin_sig: None,
                },
            );
        }

        // A code set twice on a feature keeps the last row, as the adaptor's
        // per-feature hash does (PhenotypeFeatureAdaptor.pm:1273-1276 assigns each
        // `code=value` in turn).
        for row in read_table(dir, "phenotype_feature_attrib", 3)? {
            let row = row?;
            let (Some(code), Some(f)) = (
                attrib_codes.get(row.field(1)),
                features.get_mut(row.field(0)),
            ) else {
                continue;
            };
            let value = Some(row.field(2).to_string());
            match code.as_str() {
                "risk_allele" => f.risk_allele = value,
                "clinvar_clin_sig" => f.clinvar_clin_sig = value,
                "clinvar_ref_allele" => f.clinvar_ref_allele = value,
                "somatic_clin_sig" => f.somatic_clin_sig = value,
                "oncogenic_clin_sig" => f.oncogenic_clin_sig = value,
                _ => {}
            }
        }

        let mut ordered: Vec<(u64, Feature)> = features
            .into_iter()
            .map(|(id, f)| (id.parse::<u64>().unwrap_or(u64::MAX), f))
            .collect();
        ordered.sort_by_key(|(id, _)| *id);

        let mut out = Self::default();
        for (_, f) in ordered {
            if f.somatic {
                let feature = SomaticClinVar {
                    phenotype: phenotypes.get(&f.phenotype_id).cloned().unwrap_or_default(),
                    somatic_clin_sig: f.somatic_clin_sig,
                    oncogenic_clin_sig: f.oncogenic_clin_sig,
                };
                let spans = out.somatic_spans.entry(f.name_key).or_default();
                if !spans.contains(&f.span_key) {
                    spans.push(f.span_key.clone());
                }
                out.somatic.entry(f.span_key).or_default().push(feature);
            } else if f.risk_allele.is_some() || f.clinvar_ref_allele.is_some() {
                out.germline
                    .entry(f.name_key)
                    .or_default()
                    .push(GermlineClinVar {
                        risk_allele: f.risk_allele,
                        clinvar_clin_sig: f.clinvar_clin_sig,
                        clinvar_ref_allele: f.clinvar_ref_allele,
                    });
            }
        }
        info!(
            "Read ClinVar features of {} germline variants and {} somatic spans from {}",
            out.germline.len(),
            out.somatic.len(),
            dir.display()
        );
        Ok(out)
    }

    /// Whether no feature was read.
    pub fn is_empty(&self) -> bool {
        self.germline.is_empty() && self.somatic.is_empty()
    }

    /// Fill `clin_sig_allele`, `clin_sig_ref_allele` and `clinical_impact` on a record
    /// at `chr` as `Database/Variation.pm:196-221` does.
    ///
    /// A germline feature counts when its variation name is the record's: each
    /// `risk_allele:clinvar_clin_sig` pair (spaces in the term read `_`) once, joined
    /// by `;`, and its `clinvar_ref_allele` as the reference allele. Every somatic
    /// feature on the record's database span counts, whatever its name, formatted by
    /// `_format_clinical_impact` (`:252-275`): the one span the record's own somatic
    /// features lie on when they lie on one, the tables holding the database's span
    /// where the record's may be the trimmed one, else the record's own. The
    /// germline pairs keep their first-seen order where Perl's hash keys have none.
    fn apply(&self, var: &mut VariationRecord, chr: &str) {
        let name_key = format!("{chr}:{}", var.variation_name);

        if let Some(features) = self.germline.get(&name_key) {
            let mut pairs: Vec<String> = Vec::new();
            let mut ref_allele = None;
            for f in features {
                let Some(clin_sig) = &f.clinvar_clin_sig else {
                    continue;
                };
                let pair = format!(
                    "{}:{}",
                    f.risk_allele.as_deref().unwrap_or(""),
                    clin_sig.replace(' ', "_")
                );
                if !pairs.contains(&pair) {
                    pairs.push(pair);
                }
                if let Some(r) = f.clinvar_ref_allele.as_deref().filter(|r| !r.is_empty()) {
                    ref_allele = Some(r.to_string());
                }
            }
            if !pairs.is_empty() {
                var.clin_sig_allele = Some(pairs.join(";"));
            }
            var.clin_sig_ref_allele = ref_allele;
        }

        let span_key = format!("{chr}:{}-{}", var.start, var.end);
        let database_span = match self.somatic_spans.get(&name_key).map(Vec::as_slice) {
            Some([span]) => span,
            _ => &span_key,
        };
        if let Some(features) = self.somatic.get(database_span) {
            let impact = format_clinical_impact(features);
            if !impact.is_empty() {
                var.clinical_impact = Some(impact);
            }
        }
    }
}

/// `Database/Variation.pm::_format_clinical_impact` (`:252-275`): for each somatic
/// feature, `phenotype(somatic_clin_sig)` with the classification's commas read `;`,
/// then `phenotype(oncogenicity:oncogenic_clin_sig)`, every entry joined by `;`.
fn format_clinical_impact(features: &[SomaticClinVar]) -> String {
    let mut entries = Vec::new();
    for f in features {
        if let Some(sig) = f.somatic_clin_sig.as_deref().filter(|s| !s.is_empty()) {
            entries.push(format!("{}({})", f.phenotype, sig.replace(',', ";")));
        }
        if let Some(onc) = f.oncogenic_clin_sig.as_deref().filter(|s| !s.is_empty()) {
            entries.push(format!("{}(oncogenicity:{onc})", f.phenotype));
        }
    }
    entries.join(";")
}

/// The genome-wide variation VCFs read beside a chromosome's own file.
#[derive(Debug, Default)]
pub struct GenomeWideInputs {
    /// The 1000 Genomes frequency VCF.
    pub frequencies: Option<PathBuf>,
    /// The somatic VCF: every record of it is a somatic variant.
    pub somatic: Option<PathBuf>,
    /// VCFs whose records are germline variants the per-chromosome file may
    /// lack (the phenotype- and clinically-associated files).
    pub germline: Vec<PathBuf>,
}

/// Streams the variation records of chromosome `chr` to `emit` in start order:
/// the per-chromosome VCF's, plus the records of the genome-wide files absent
/// from it (a record is the same one when its name, position and alleles are),
/// the somatic file's flagged `somatic`, every non-somatic record given its 1000
/// Genomes frequencies and, with `clinvar`, its ClinVar columns. The
/// per-chromosome file is read in its position order and reordered through a
/// buffer, since a record's start lies within two bases of its position;
/// records sharing a start come out in name order ([`StartOrder`]). Returns the
/// number of records emitted.
pub fn chromosome_records(
    chr: &str,
    vcf: &Path,
    inputs: &GenomeWideInputs,
    clinvar: Option<&ClinVarAttribs>,
    mut emit: impl FnMut(VariationRecord) -> Result<()>,
) -> Result<usize> {
    let frequencies = match &inputs.frequencies {
        Some(path) => Some(ThousandGenomes::load(path, chr)?),
        None => None,
    };

    let mut held: HashMap<RecordIdentity, VariationRecord> = HashMap::new();
    // The somatic file first, so a variant it shares with a germline file
    // keeps its flag.
    let genome_wide = inputs
        .somatic
        .iter()
        .map(|path| (path, true))
        .chain(inputs.germline.iter().map(|path| (path, false)));
    for (path, somatic) in genome_wide {
        parse_chromosome(
            path,
            Some(chr),
            somatic,
            frequencies.as_ref(),
            |fields, record| {
                held.entry(fields.identity()).or_insert(record);
                Ok(())
            },
        )?;
    }

    let mut emitted = 0usize;
    let mut merge = StartOrder::new(held, |mut record: VariationRecord| {
        if let Some(clinvar) = clinvar {
            clinvar.apply(&mut record, chr);
        }
        emitted += 1;
        emit(record)
    });
    parse_chromosome(vcf, None, false, frequencies.as_ref(), |fields, record| {
        merge.push(fields, record)
    })?;
    merge.finish()?;
    info!("chr {}: {} variation records", chr, emitted);
    Ok(emitted)
}

/// The start-ordered merge of a chromosome's records: the per-chromosome
/// file's, each held back until the file's position passes its start, and the
/// genome-wide files' records the per-chromosome file lacks. A record's start
/// never precedes its position, so once the file's position passes a start
/// nothing before that start is still to come. Records sharing a start come out
/// in name order, the order the dumper's sort on the start column leaves them in
/// (`Dumper/Variation.pm:120`: `sort -k5,5n`, whose equal keys fall back to the
/// whole line, which the name follows the chromosome on).
struct StartOrder<F: FnMut(VariationRecord) -> Result<()>> {
    /// The per-chromosome file's records held back until the file's position
    /// passes their start, keyed by start, name and the order they were read in.
    buffered: BTreeMap<(u64, String, usize), VariationRecord>,
    sequence: usize,
    /// The genome-wide records by start and name, descending, so the next one is
    /// last.
    added: Vec<(u64, String, RecordIdentity)>,
    /// The genome-wide records; one leaves, written in the copy's place, when
    /// the per-chromosome file yields its copy.
    held: HashMap<RecordIdentity, VariationRecord>,
    emit: F,
}

impl<F: FnMut(VariationRecord) -> Result<()>> StartOrder<F> {
    fn new(held: HashMap<RecordIdentity, VariationRecord>, emit: F) -> Self {
        let mut added: Vec<(u64, String, RecordIdentity)> = held
            .iter()
            .map(|(identity, record)| {
                (
                    record.start,
                    record.variation_name.clone(),
                    identity.clone(),
                )
            })
            .collect();
        added.sort_by(|a, b| b.cmp(a));
        Self {
            buffered: BTreeMap::new(),
            sequence: 0,
            added,
            held,
            emit,
        }
    }

    /// Takes the next per-chromosome record, read at `fields.pos`: emits every
    /// record starting before that position, then holds this one back. A
    /// genome-wide copy of the record takes its place, so a variant the somatic
    /// file holds keeps that file's flag and goes without frequencies.
    fn push(&mut self, fields: &VcfFields, record: VariationRecord) -> Result<()> {
        let record = if self.held.is_empty() {
            record
        } else {
            self.held.remove(&fields.identity()).unwrap_or(record)
        };
        self.flush(Some(fields.pos))?;
        self.sequence += 1;
        let key = (record.start, record.variation_name.clone(), self.sequence);
        self.buffered.insert(key, record);
        Ok(())
    }

    /// Emits everything still held.
    fn finish(mut self) -> Result<()> {
        self.flush(None)
    }

    /// Emits, by start then name, the held records whose start precedes `before`
    /// (every held record when `None`).
    fn flush(&mut self, before: Option<u64>) -> Result<()> {
        loop {
            let next_buffered = self
                .buffered
                .keys()
                .next()
                .map(|(start, name, _)| (*start, name.as_str()));
            let next_added = self
                .added
                .last()
                .map(|(start, name, _)| (*start, name.as_str()));
            let (take_buffered, start) = match (next_buffered, next_added) {
                (None, None) => return Ok(()),
                (Some(buffered), None) => (true, buffered.0),
                (None, Some(added)) => (false, added.0),
                (Some(buffered), Some(added)) if buffered <= added => (true, buffered.0),
                (Some(_), Some(added)) => (false, added.0),
            };
            if before.is_some_and(|position| start >= position) {
                return Ok(());
            }
            if take_buffered {
                let (_, record) = self.buffered.pop_first().expect("peeked");
                (self.emit)(record)?;
            } else {
                let (_, _, identity) = self.added.pop().expect("peeked");
                // A genome-wide record the per-chromosome file also holds
                // left the table when its copy was read.
                if let Some(record) = self.held.remove(&identity) {
                    (self.emit)(record)?;
                }
            }
        }
    }
}

/// The `variations/<chr>/all_vars.gz` of one chromosome and its `.csi`, written
/// record by record in start order.
pub struct AllVarsWriter {
    path: PathBuf,
    index_path: PathBuf,
    chr: String,
    columns: Vec<String>,
    writer: bgzf::io::Writer<std::fs::File>,
    indexer: Indexer<BinnedIndex>,
    line: String,
    last_start: u64,
    count: usize,
}

impl AllVarsWriter {
    /// Creates `<output_dir>/variations/<chr>/all_vars.gz` for the records of
    /// `chr`, with a cell per column of `columns`.
    pub fn create(output_dir: &Path, chr: &str, columns: &[String]) -> Result<Self> {
        if !is_cache_path_component(chr) {
            bail!("invalid chromosome name '{chr}'");
        }
        let dir = output_dir.join("variations").join(chr);
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("Failed to create {}", dir.display()))?;
        let path = dir.join("all_vars.gz");
        let file = std::fs::File::create(&path)
            .with_context(|| format!("Failed to create {}", path.display()))?;
        Ok(Self {
            index_path: dir.join("all_vars.gz.csi"),
            path,
            chr: chr.to_string(),
            columns: columns.to_vec(),
            writer: bgzf::io::Writer::new(file),
            // The index parameters of the Ensembl VEP cache's own `.csi` files.
            indexer: Indexer::<BinnedIndex>::new(14, 8),
            line: String::new(),
            last_start: 0,
            count: 0,
        })
    }

    /// Writes one record as `Dumper/Variation.pm:262-264` prints it: the
    /// chromosome, then a tab-separated cell per column with `.` for an empty
    /// one; and indexes it by its start. A record starting before the previous
    /// one is an error: the index is built on the start column in file order.
    pub fn write(&mut self, record: &VariationRecord) -> Result<()> {
        if record.start < self.last_start {
            anyhow::bail!(
                "{} starts at {}, before the record written before it ({}): the variation \
                 records of chr {} are not in position order",
                record.variation_name,
                record.start,
                self.last_start,
                self.chr
            );
        }
        self.last_start = record.start;
        self.line.clear();
        self.line.push_str(&self.chr);
        for column in self.columns.iter().skip(1) {
            self.line.push('\t');
            let cell = record.cell(column);
            self.line
                .push_str(if cell.is_empty() { "." } else { &cell });
        }
        self.line.push('\n');
        let begin = self.writer.virtual_position();
        self.writer
            .write_all(self.line.as_bytes())
            .with_context(|| format!("Failed to write {}", self.path.display()))?;
        let end = self.writer.virtual_position();
        let start = Position::try_from(record.start as usize)
            .with_context(|| format!("{} has start {}", record.variation_name, record.start))?;
        self.indexer
            .add_record(Some((0, start, start, true)), Chunk::new(begin, end))
            .with_context(|| format!("Failed to index {}", self.path.display()))?;
        self.count += 1;
        Ok(())
    }

    /// Finishes the bgzip file and writes its `.csi` (`:131`: `tabix -C -s 1 -b 5
    /// -e 5`, the start column alone, so a record is found by its start and never
    /// by its span); returns the number of records written.
    pub fn finish(mut self) -> Result<usize> {
        self.writer
            .try_finish()
            .with_context(|| format!("Failed to finish {}", self.path.display()))?;
        let header = Header::builder()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_name_index(CHR_COLUMN)
            .set_start_position_index(START_COLUMN)
            .set_end_position_index(Some(START_COLUMN))
            .set_reference_sequence_names([self.chr.as_bytes().into()].into_iter().collect())
            .build();
        let index = self.indexer.set_header(header).build(1);
        noodles::csi::fs::write(&self.index_path, &index)
            .with_context(|| format!("Failed to write {}", self.index_path.display()))?;
        info!(
            "chr {}: wrote {} ({} records)",
            self.chr,
            self.path.display(),
            self.count
        );
        Ok(self.count)
    }

    /// Removes the file and any index beside it, leaving the chromosome without
    /// a variation file rather than with a partial one.
    pub fn discard(self) {
        let _ = std::fs::remove_file(&self.path);
        let _ = std::fs::remove_file(&self.index_path);
    }
}

/// Writes chromosome `chr`'s variation file from its inputs: the records of
/// [`chromosome_records`] through an [`AllVarsWriter`] with `columns`. Returns
/// the number of records written; on an error the partial file and any index
/// beside it are removed.
pub fn write_chromosome(
    chr: &str,
    vcf: &Path,
    inputs: &GenomeWideInputs,
    clinvar: Option<&ClinVarAttribs>,
    columns: &[String],
    output_dir: &Path,
) -> Result<usize> {
    let mut writer = AllVarsWriter::create(output_dir, chr, columns)?;
    match chromosome_records(chr, vcf, inputs, clinvar, |record| writer.write(&record)) {
        Ok(_) => writer.finish(),
        Err(error) => {
            writer.discard();
            Err(error)
        }
    }
}

impl VariationRecord {
    /// The text of one cell, empty when the dumper leaves it empty
    /// (`Dumper/Variation.pm:233-250`: `failed`, `somatic` and
    /// `phenotype_or_disease` when 0, `end` when it equals `start`, `strand`
    /// when 1).
    fn cell(&self, column: &str) -> String {
        let flag = |value: u8| {
            if value == 0 {
                String::new()
            } else {
                value.to_string()
            }
        };
        match column {
            "variation_name" => self.variation_name.clone(),
            "failed" => flag(self.failed),
            "somatic" => flag(self.somatic),
            "start" => self.start.to_string(),
            "end" => {
                if self.end == self.start {
                    String::new()
                } else {
                    self.end.to_string()
                }
            }
            "allele_string" => self.allele_string.clone(),
            "strand" => {
                if self.strand == 1 {
                    String::new()
                } else {
                    self.strand.to_string()
                }
            }
            "clin_sig" => self.clin_sig.clone().unwrap_or_default(),
            "phenotype_or_disease" => flag(self.phenotype_or_disease),
            "clin_sig_allele" => self.clin_sig_allele.clone().unwrap_or_default(),
            "clin_sig_ref_allele" => self.clin_sig_ref_allele.clone().unwrap_or_default(),
            "clinical_impact" => self.clinical_impact.clone().unwrap_or_default(),
            "pubmed" | "var_synonyms" => String::new(),
            population => self
                .frequencies
                .get(population)
                .cloned()
                .unwrap_or_default(),
        }
    }
}

/// What makes two VCF records the same variant across Ensembl's files, which
/// export every variant identically wherever it appears.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct RecordIdentity {
    name: String,
    pos: u64,
    reference: String,
    alternates: String,
}

/// The fields of a VCF record this module reads.
struct VcfFields {
    chr: String,
    pos: u64,
    name: String,
    reference: String,
    alternates: String,
    filter: String,
    info: String,
}

impl VcfFields {
    fn read(record: &vcf::Record) -> Result<Self> {
        Ok(Self {
            chr: record.reference_sequence_name().to_string(),
            pos: record
                .variant_start()
                .context("Missing POS")?
                .context("Invalid POS")?
                .get() as u64,
            name: record.ids().as_ref().to_string(),
            reference: record.reference_bases().to_string(),
            alternates: record.alternate_bases().as_ref().to_string(),
            filter: record.filters().as_ref().to_string(),
            info: record.info().as_ref().to_string(),
        })
    }

    fn identity(&self) -> RecordIdentity {
        RecordIdentity {
            name: self.name.clone(),
            pos: self.pos,
            reference: self.reference.clone(),
            alternates: self.alternates.clone(),
        }
    }

    /// The INFO fields, each `key` or `key=value`.
    fn info_fields(&self) -> impl Iterator<Item = &str> {
        self.info
            .split(';')
            .filter(|field| !field.is_empty() && *field != ".")
    }

    fn has_flag(&self, key: &str) -> bool {
        self.info_fields().any(|field| field == key)
    }

    fn info_value(&self, key: &str) -> Option<&str> {
        self.info_fields()
            .find_map(|field| field.split_once('=').filter(|(k, _)| *k == key))
            .map(|(_, value)| value)
    }
}

/// Reads the records of `path` on chromosome `chr` (every record when `None`)
/// and hands each, with its VCF fields, to `f` as a variation record: through
/// the file's `.csi` or `.tbi` when one lies beside it, else by reading the file
/// through and keeping the chromosome's records. Records whose FILTER is not
/// `PASS`, `.` or empty are skipped, as are those the database form has no
/// shape for.
fn parse_chromosome(
    path: &Path,
    chr: Option<&str>,
    somatic: bool,
    frequencies: Option<&ThousandGenomes>,
    mut f: impl FnMut(&VcfFields, VariationRecord) -> Result<()>,
) -> Result<()> {
    let mut skipped = 0u64;
    let mut parsed = 0u64;
    let malformed = for_each_record(path, chr, |fields| {
        match parse_vcf_fields(&fields, somatic, frequencies) {
            Some(record) => {
                parsed += 1;
                f(&fields, record)
            }
            None => {
                skipped += 1;
                Ok(())
            }
        }
    })?;
    info!(
        "Parsed {} variations from {} ({} skipped)",
        parsed,
        path.display(),
        skipped + malformed
    );
    Ok(())
}

/// Calls `f` with each record of `path` on `chr`, or with every record when
/// `chr` is `None`; returns the number of records skipped as malformed.
fn for_each_record(
    path: &Path,
    chr: Option<&str>,
    mut f: impl FnMut(VcfFields) -> Result<()>,
) -> Result<u64> {
    if let Some(chr) = chr {
        if !has_index(path) {
            warn!(
                "No index beside {}: reading the file through for chromosome {}",
                path.display(),
                chr
            );
        } else if let Some(malformed) = query_chromosome(path, chr, &mut f)? {
            return Ok(malformed);
        } else {
            warn!(
                "The index beside {} names no sequence: reading the file through for \
                 chromosome {}",
                path.display(),
                chr
            );
        }
    }
    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open VCF: {}", path.display()))?;
    let reader: Box<dyn BufRead> = if is_compressed(path) {
        Box::new(std::io::BufReader::new(bgzf::io::Reader::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let mut vcf_reader = vcf::io::Reader::new(reader);
    vcf_reader
        .read_header()
        .with_context(|| format!("Failed to read VCF header: {}", path.display()))?;
    read_records(
        path,
        |record| vcf_reader.read_record(record),
        |fields| {
            if chr.is_some_and(|chr| fields.chr != chr) {
                return Ok(());
            }
            f(fields)
        },
    )
}

/// Calls `f` with each record of chromosome `chr` of the indexed VCF `path` and
/// returns the number skipped as malformed, or `None` without a call when the
/// index names no sequence at all (an index without the tabix header, which the
/// caller then reads around); a chromosome the index names no records for yields
/// none.
fn query_chromosome(
    path: &Path,
    chr: &str,
    f: &mut impl FnMut(VcfFields) -> Result<()>,
) -> Result<Option<u64>> {
    let mut reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("Failed to open indexed VCF: {}", path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read VCF header: {}", path.display()))?;
    let Some(index_header) = reader.index().header() else {
        return Ok(None);
    };
    if !index_header
        .reference_sequence_names()
        .contains(chr.as_bytes())
    {
        return Ok(Some(0));
    }
    let region = Region::new(chr, ..);
    let mut query = reader
        .query(&header, &region)
        .with_context(|| format!("Failed to query {} in {}", chr, path.display()))?;
    read_records(path, |record| query.read_record(record), f).map(Some)
}

/// Calls `f` with the fields of each record `read_next` yields until it yields
/// none; returns the number of records skipped because their fields do not read.
fn read_records(
    path: &Path,
    mut read_next: impl FnMut(&mut vcf::Record) -> std::io::Result<usize>,
    mut f: impl FnMut(VcfFields) -> Result<()>,
) -> Result<u64> {
    let mut record = vcf::Record::default();
    let mut malformed = 0u64;
    loop {
        let read = read_next(&mut record)
            .with_context(|| format!("Failed to read a record of {}", path.display()))?;
        if read == 0 {
            return Ok(malformed);
        }
        match VcfFields::read(&record) {
            Ok(fields) => f(fields)?,
            Err(e) => {
                malformed += 1;
                debug!("Skipping malformed VCF record: {}", e);
            }
        }
    }
}

fn is_compressed(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
}

/// Whether a `.csi` or `.tbi` index lies beside the compressed file `path`.
fn has_index(path: &Path) -> bool {
    is_compressed(path)
        && ["csi", "tbi"].iter().any(|ext| {
            let mut name = path.as_os_str().to_owned();
            name.push(".");
            name.push(ext);
            Path::new(&name).is_file()
        })
}

/// Parse one VCF record into its database form, or `None` when it is skipped:
/// a filtered record, or one without alleles whose source names no class.
fn parse_vcf_fields(
    fields: &VcfFields,
    somatic: bool,
    frequencies: Option<&ThousandGenomes>,
) -> Option<VariationRecord> {
    if !fields.filter.is_empty() && fields.filter != "." && fields.filter != "PASS" {
        return None;
    }
    let variation_name = if fields.name == "." {
        String::new()
    } else {
        fields.name.split(';').next().unwrap_or("").to_string()
    };
    let alts: Vec<&str> = fields.alternates.split(',').collect();
    let source = Source::of(fields);
    let (start, end, allele_string) = database_form(
        fields.pos,
        &fields.reference,
        &alts,
        fields.info_value("TSA"),
        source,
    )?;

    let clin_sig_terms: Vec<&str> = fields
        .info_fields()
        .filter_map(|field| field.strip_prefix("CLIN_"))
        .collect();
    let mut record = VariationRecord {
        variation_name,
        failed: 0,
        somatic: u8::from(somatic),
        start,
        end,
        strand: 1,
        allele_string,
        clin_sig: (!clin_sig_terms.is_empty()).then(|| clin_sig_terms.join(",")),
        phenotype_or_disease: u8::from(fields.has_flag("E_Phenotype_or_Disease")),
        ..Default::default()
    };
    // The dumper joins frequencies onto non-somatic variants only
    // (`Dumper/Variation.pm:280`).
    if let Some(frequencies) = frequencies.filter(|_| !somatic) {
        frequencies.apply(&mut record, fields.pos, &fields.reference, &alts);
    }
    Some(record)
}

/// The source of a variant without alleles, read from the source flag the
/// Ensembl VCFs set on every record (`COSMIC_102`, `HGMD-PUBLIC_20204`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Source {
    Cosmic,
    Hgmd,
    Other,
}

impl Source {
    fn of(fields: &VcfFields) -> Self {
        let mut source = Self::Other;
        for field in fields.info_fields() {
            if field.contains('=') {
                continue;
            }
            if field.starts_with("COSMIC") {
                source = Self::Cosmic;
            } else if field.starts_with("HGMD") {
                source = Self::Hgmd;
            }
        }
        source
    }
}

/// The `start`, `end` and allele string the variation database holds for a VCF
/// record, or `None` for a record without alleles whose source names no class.
///
/// The export writes a variant's database alleles into a VCF line
/// (`VariationFeature.pm:2340-2420`, `to_VCF_record`): an indel's alleles are
/// trimmed of their shared trailing bases down to one base, and the base before
/// the variant is prepended to every allele, with `POS` moved onto it, only when
/// the alleles do not all share a first base, which an empty allele never does.
/// A record with alleles therefore takes the form Ensembl VEP's VCF parser gives
/// it ([`ensembl_alleles`]: the shared first base chopped) when the database
/// holds an empty allele (`TSA=deletion`, `TSA=insertion`, the mixed
/// `sequence_alteration`), except that a record marked `TSA=insertion` lies one
/// base further right: the database holds such a variant with its start one past
/// its end, the GVF export collapses the start onto the end, and the VCF export
/// then anchors the record on the base before that, so the base the VCF names is
/// two before the database start. A record marked `TSA=indel` has no empty allele
/// in the database, so it keeps the export's alleles at `POS` ([`indel_form`])
/// unless the chop leaves alleles that neither include an empty one nor share a
/// first base, the shape the export's anchor leaves behind (a database record
/// spelling one shared base before alleles that then differ, `AC/AGT`, would be
/// read as anchored; the sources' minimal spellings carry no such base).
/// An HGMD-PUBLIC variant keeps its coordinates but its allele string reads
/// `HGMD_MUTATION`, the string the database holds for every variant of that
/// source.
///
/// A record whose ALT is `.` or `<.>` is a COSMIC or HGMD-PUBLIC variant whose
/// alleles the database does not hold: its allele string is the source's class
/// name and its span is read back from the shape the export gives it. `.` on a
/// reference base is a substitution or one-base deletion of that base, or, under
/// `TSA=insertion`, an insertion after it; `<.>` after an anchor base spans the
/// bases that follow the anchor, or, under `TSA=insertion`, an insertion two
/// bases on, as for a record with alleles.
fn database_form(
    pos: u64,
    reference: &str,
    alts: &[&str],
    tsa: Option<&str>,
    source: Source,
) -> Option<(u64, u64, String)> {
    let class_name = match source {
        Source::Cosmic => Some(COSMIC_MUTATION),
        Source::Hgmd => Some(HGMD_MUTATION),
        Source::Other => None,
    };
    let marked_insertion = tsa == Some("insertion");
    let reference_end = pos + reference.len() as u64 - 1;
    match alts {
        ["."] | [""] => {
            let name = class_name?;
            if marked_insertion {
                Some((pos + 1, pos, name.to_string()))
            } else {
                Some((pos, reference_end, name.to_string()))
            }
        }
        ["<.>"] => {
            let name = class_name?;
            if marked_insertion {
                Some((pos + 2, pos + 1, name.to_string()))
            } else {
                Some((pos + 1, reference_end, name.to_string()))
            }
        }
        _ => {
            let alts: Vec<String> = alts.iter().map(|alt| alt.to_string()).collect();
            let (mut start, mut end, allele_string) = if tsa == Some("indel") {
                indel_form(pos, reference, &alts)
            } else {
                ensembl_alleles(pos, reference, &alts)
            };
            if marked_insertion && allele_string.starts_with("-/") {
                start += 1;
                end += 1;
            }
            let allele_string = match source {
                Source::Hgmd => HGMD_MUTATION.to_string(),
                _ => allele_string,
            };
            Some((start, end, allele_string))
        }
    }
}

/// The `start`, `end` and allele string Ensembl VEP's VCF parser gives a record
/// (`Parser/VCF.pm:254, 297-334`).
///
/// `end` is `POS + len(REF) - 1`. The record is an indel when any ALT's length
/// differs from REF's (`:297`); then, when REF and every ALT share their first base
/// (with one ALT, `:328`; with several, `:306-308`, an ALT carrying `*` left out of the
/// comparison), that base is chopped from each allele but a `*` one, `start` becomes
/// `POS + 1` and an emptied allele reads `-` (`:309-316, 331-334`), so an insertion
/// ends one base before it starts. Same-length alleles are never chopped, nor is a
/// record whose REF opens with a non-ASCII character: the chop is a one-byte slice, a
/// character boundary only in ASCII.
fn ensembl_alleles(pos: u64, ref_allele: &str, alts: &[String]) -> (u64, u64, String) {
    let end = pos + ref_allele.len() as u64 - 1;
    let is_indel = alts.iter().any(|alt| alt.len() != ref_allele.len());
    let first_base = ref_allele.as_bytes().first().filter(|b| b.is_ascii());
    let shares_first_base = |alt: &str| alt.as_bytes().first() == first_base;
    let carries_star = |allele: &str| allele.contains('*');
    let anchored = is_indel
        && first_base.is_some()
        && match alts {
            [single] => shares_first_base(single),
            _ => alts
                .iter()
                .filter(|alt| !carries_star(alt))
                .all(|alt| shares_first_base(alt)),
        };

    let chop = |allele: &str| -> String {
        if !anchored || carries_star(allele) {
            return allele.to_string();
        }
        match &allele[1..] {
            "" => "-".to_string(),
            rest => rest.to_string(),
        }
    };

    let start = if anchored { pos + 1 } else { pos };
    let mut allele_string = chop(ref_allele);
    for alt in alts {
        allele_string.push('/');
        allele_string.push_str(&chop(alt));
    }
    (start, end, allele_string)
}

/// The `start`, `end` and allele string of a record marked `TSA=indel`: the
/// export's alleles at `POS`, with `end` at `POS + len(REF) - 1`, unless the
/// chopped form ([`ensembl_alleles`]) leaves alleles that neither include an empty
/// one nor all share a first base, which is then the database's form.
///
/// The database holds an indel with every allele non-empty, and the export
/// prepends an anchor base only when the alleles do not all share a first base
/// (`VariationFeature.pm:2400-2408`). A chop that empties an allele or leaves a
/// shared first base therefore removed a base the database holds (`A/AA` stays
/// `A/AA`; `AGA/AGAGA` stays), while a chop leaving alleles that differ at once
/// (`CC/CTTCT` to `C/TTCT`) removed the anchor the export added. The record's
/// start and alleles are then the database's, trailing bases the export trimmed
/// aside, so a left-aligned input co-locates with it exactly where Ensembl VEP
/// co-locates it, which the chopped form, one base to the right, never does.
fn indel_form(pos: u64, ref_allele: &str, alts: &[String]) -> (u64, u64, String) {
    let (start, end, chopped) = ensembl_alleles(pos, ref_allele, alts);
    let anchor_chopped = start > pos;
    let alleles: Vec<&str> = chopped
        .split('/')
        .filter(|allele| !allele.contains('*'))
        .collect();
    let none_empty = alleles.iter().all(|allele| *allele != "-");
    let share_first_base = alleles
        .windows(2)
        .all(|pair| pair[0].as_bytes().first() == pair[1].as_bytes().first());
    if anchor_chopped && none_empty && !share_first_base {
        (start, end, chopped)
    } else {
        let mut allele_string = ref_allele.to_string();
        for alt in alts {
            allele_string.push('/');
            allele_string.push_str(alt);
        }
        (pos, end, allele_string)
    }
}

/// The 1000 Genomes frequencies of one chromosome, keyed by minimised allele.
///
/// The frequency VCF gives each ALT its frequency in the five continental
/// populations and no frequency in the whole panel. The dumper's join
/// (`Dumper/Variation.pm:272-420`) matches a VCF allele to a variant's allele
/// through `get_matched_variant_alleles` (`Utils/Sequence.pm:1144-1245`): each
/// REF/ALT pair is trimmed to its minimal form (`trim_sequences`, `:965-1031`),
/// from the left and, when an allele exceeds one base, from the right as well,
/// and two alleles match when a minimal form and its position agree.
struct ThousandGenomes {
    /// Minimised `ref_alt_pos` key to the frequency record (its index in `alleles`)
    /// and the ALT the key came from.
    by_allele: HashMap<String, (u32, u8)>,
    /// Per frequency record in file order, per ALT, the five populations'
    /// frequencies in [`THOUSAND_GENOMES_POPULATIONS`] order, NaN for a `.`.
    alleles: Vec<Vec<[f64; 5]>>,
}

impl ThousandGenomes {
    /// Reads the frequency VCF's records on chromosome `chr`.
    fn load(path: &Path, chr: &str) -> Result<Self> {
        let mut by_allele = HashMap::new();
        let mut alleles: Vec<Vec<[f64; 5]>> = Vec::new();
        let malformed = for_each_record(path, Some(chr), |fields| {
            let values: Vec<Vec<&str>> = THOUSAND_GENOMES_POPULATIONS
                .iter()
                .map(|(population, _)| {
                    fields
                        .info_value(population)
                        .map(|value| value.split(',').collect())
                        .unwrap_or_default()
                })
                .collect();
            let record_index = u32::try_from(alleles.len()).context("frequency record count")?;
            let mut record_alleles = Vec::new();
            for (index, alt) in fields.alternates.split(',').enumerate() {
                let mut entry = [f64::NAN; 5];
                for (slot, population_values) in entry.iter_mut().zip(&values) {
                    if let Some(value) = population_values.get(index) {
                        *slot = value.parse().unwrap_or(f64::NAN);
                    }
                }
                record_alleles.push(entry);
                let alt_index = u8::try_from(index).context("ALT count")?;
                for key in minimised_keys(fields.pos, &fields.reference, alt) {
                    by_allele.entry(key).or_insert((record_index, alt_index));
                }
            }
            alleles.push(record_alleles);
            Ok(())
        })?;
        info!(
            "chr {}: {} frequency records read from {} ({} skipped)",
            chr,
            alleles.len(),
            path.display(),
            malformed
        );
        Ok(Self { by_allele, alleles })
    }

    /// Fills the record's 1000 Genomes columns from the frequencies of its
    /// alleles: each population's text is `allele:frequency` per matched allele,
    /// in the order of the frequency records matched and of their ALTs, as the
    /// dumper joins them (`:349-399`), with the frequency printed as it prints it
    /// (`:378`: `0`, else `%.4g`). An allele whose frequency is 0 in every
    /// population gets no entry: the frequency VCF lists a frequency for every
    /// allele of a variant, the 1000 Genomes release only for the alleles its call
    /// set holds, and 0 everywhere is the one sign of an allele it does not list.
    /// `AF`, the whole panel's frequency, is the allele's count over the panel's
    /// 5,008 chromosomes: a population's count is its frequency times its
    /// chromosomes, rounded to a whole number.
    fn apply(&self, record: &mut VariationRecord, pos: u64, reference: &str, alts: &[&str]) {
        let cache_alleles: Vec<&str> = record.allele_string.split('/').skip(1).collect();
        let mut matched: Vec<((u32, u8), &str)> = alts
            .iter()
            .zip(cache_alleles)
            .filter_map(|(alt, cache_allele)| {
                minimised_keys(pos, reference, alt)
                    .into_iter()
                    .find_map(|key| self.by_allele.get(&key))
                    .map(|&entry| (entry, cache_allele))
            })
            .collect();
        matched.sort_by_key(|(entry, _)| *entry);
        matched.dedup_by_key(|(entry, _)| *entry);

        let mut columns: Vec<Vec<String>> = vec![Vec::new(); THOUSAND_GENOMES_COLUMNS.len()];
        for ((record_index, alt_index), cache_allele) in matched {
            let frequencies = &self.alleles[record_index as usize][alt_index as usize];
            if frequencies.iter().all(|f| f.is_nan() || *f == 0.0) {
                continue;
            }
            let mut count = 0.0f64;
            let mut chromosomes = 0u64;
            for ((_, samples), frequency) in THOUSAND_GENOMES_POPULATIONS.iter().zip(frequencies) {
                if frequency.is_nan() {
                    continue;
                }
                let population_chromosomes = 2 * samples;
                count += (frequency * population_chromosomes as f64).round();
                chromosomes += population_chromosomes;
            }
            for (column, frequency) in columns.iter_mut().skip(1).zip(frequencies) {
                if !frequency.is_nan() {
                    column.push(format!("{cache_allele}:{}", perl_g4(*frequency)));
                }
            }
            if chromosomes > 0 {
                columns[0].push(format!(
                    "{cache_allele}:{}",
                    perl_g4(count / chromosomes as f64)
                ));
            }
        }
        for (column, entries) in THOUSAND_GENOMES_COLUMNS.iter().zip(columns) {
            if !entries.is_empty() {
                record
                    .frequencies
                    .insert(column.to_string(), entries.join(","));
            }
        }
    }
}

/// The `ref_alt_pos` keys `get_matched_variant_alleles` builds for one REF/ALT
/// pair (`Utils/Sequence.pm:1198-1210`): the pair trimmed from the left then the
/// right and, when either allele exceeds one base, from the right then the left,
/// an emptied allele reading `-`.
fn minimised_keys(pos: u64, reference: &str, alt: &str) -> Vec<String> {
    let mut keys = vec![trim_key(pos, reference, alt, false)];
    if reference.len() > 1 || alt.len() > 1 {
        let right_first = trim_key(pos, reference, alt, true);
        if right_first != keys[0] {
            keys.push(right_first);
        }
    }
    keys
}

/// `trim_sequences` (`Utils/Sequence.pm:965-1031`) on the forward strand, as the
/// key `ref_alt_pos`: shared bases are removed from the left, moving the position
/// right by one per base, and from the right, in the order `end_first` says, while
/// both alleles still have a base.
fn trim_key(pos: u64, reference: &str, alt: &str, end_first: bool) -> String {
    let (mut r, mut a, mut pos) = (reference.as_bytes(), alt.as_bytes(), pos);
    let trim_left = |r: &mut &[u8], a: &mut &[u8], pos: &mut u64| {
        while let (Some((r_first, r_rest)), Some((a_first, a_rest))) =
            (r.split_first(), a.split_first())
        {
            if r_first != a_first {
                break;
            }
            *r = r_rest;
            *a = a_rest;
            *pos += 1;
        }
    };
    let trim_right = |r: &mut &[u8], a: &mut &[u8]| {
        while let (Some((r_last, r_rest)), Some((a_last, a_rest))) =
            (r.split_last(), a.split_last())
        {
            if r_last != a_last {
                break;
            }
            *r = r_rest;
            *a = a_rest;
        }
    };
    if end_first {
        trim_right(&mut r, &mut a);
        trim_left(&mut r, &mut a, &mut pos);
    } else {
        trim_left(&mut r, &mut a, &mut pos);
        trim_right(&mut r, &mut a);
    }
    let dash_if_empty = |s: &[u8]| {
        if s.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(s).into_owned()
        }
    };
    format!("{}_{}_{}", dash_if_empty(r), dash_if_empty(a), pos)
}

/// A number as Perl's `sprintf('%.4g')` prints it: four significant digits,
/// trailing zeros dropped, exponent form below 0.0001 and from 10,000 up.
fn perl_g4(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let scientific = format!("{value:.3e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("exponent form");
    let exponent: i32 = exponent.parse().expect("exponent digits");
    if !(-4..4).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{sign}{:02}",
            trim_fraction_zeros(mantissa),
            exponent.abs()
        )
    } else {
        let decimals = (3 - exponent).max(0) as usize;
        trim_fraction_zeros(&format!("{value:.decimals$}"))
    }
}

fn trim_fraction_zeros(text: &str) -> String {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        text.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Read;

    /// INFO lines of an Ensembl variation VCF header: the source flags, the
    /// `E_*` evidence flags, the `CLIN_*` flags, `TSA` and the 1000 Genomes
    /// populations.
    const ENSEMBL_HEADER: &str = "##fileformat=VCFv4.1\n\
        ##INFO=<ID=dbSNP_156,Number=0,Type=Flag,Description=\"dbSNP\">\n\
        ##INFO=<ID=ClinVar_202509,Number=0,Type=Flag,Description=\"ClinVar\">\n\
        ##INFO=<ID=COSMIC_102,Number=0,Type=Flag,Description=\"COSMIC\">\n\
        ##INFO=<ID=HGMD-PUBLIC_20204,Number=0,Type=Flag,Description=\"HGMD-PUBLIC\">\n\
        ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
        ##INFO=<ID=E_Freq,Number=0,Type=Flag,Description=\"Frequency\">\n\
        ##INFO=<ID=E_Phenotype_or_Disease,Number=0,Type=Flag,Description=\"Phenotype_or_Disease\">\n\
        ##INFO=<ID=CLIN_pathogenic,Number=0,Type=Flag,Description=\"pathogenic\">\n\
        ##INFO=<ID=CLIN_likely_pathogenic,Number=0,Type=Flag,Description=\"likely pathogenic\">\n\
        ##INFO=<ID=CLIN_uncertain_significance,Number=0,Type=Flag,Description=\"uncertain significance\">\n\
        ##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n\
        ##INFO=<ID=AFR,Number=A,Type=Float,Description=\"AFR\">\n\
        ##INFO=<ID=AMR,Number=A,Type=Float,Description=\"AMR\">\n\
        ##INFO=<ID=EAS,Number=A,Type=Float,Description=\"EAS\">\n\
        ##INFO=<ID=EUR,Number=A,Type=Float,Description=\"EUR\">\n\
        ##INFO=<ID=SAS,Number=A,Type=Float,Description=\"SAS\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    /// Writes a plain-text VCF of `records` under `name` in `dir`.
    fn write_vcf(dir: &Path, name: &str, records: &str) -> PathBuf {
        let path = dir.join(name);
        fs::write(&path, format!("{ENSEMBL_HEADER}{records}")).unwrap();
        path
    }

    /// The records [`chromosome_records`] emits, in order.
    fn collect(
        vcf: &Path,
        inputs: &GenomeWideInputs,
        clinvar: Option<&ClinVarAttribs>,
    ) -> Vec<VariationRecord> {
        let mut records = Vec::new();
        let emitted = chromosome_records("21", vcf, inputs, clinvar, |record| {
            records.push(record);
            Ok(())
        })
        .unwrap();
        assert_eq!(emitted, records.len());
        records
    }

    fn records_of(records: &str) -> Vec<VariationRecord> {
        let dir = tempfile::tempdir().unwrap();
        let vcf = write_vcf(dir.path(), "homo_sapiens-chr21.vcf", records);
        collect(&vcf, &GenomeWideInputs::default(), None)
    }

    fn span(record: &VariationRecord) -> (u64, u64, &str) {
        (record.start, record.end, record.allele_string.as_str())
    }

    /// A deletion the database holds with an empty allele is exported behind an
    /// anchor base, which Ensembl VEP's VCF parser chops, moving `start` one past
    /// `POS`, leaving `end` at `POS + len(REF) - 1` and spelling the emptied allele
    /// `-` (`Parser/VCF.pm:254, 309-316, 331-334`): the export's `TA/T` at 5030218 is
    /// `A/-` at 5030219-5030219, where the release-116 cache holds rs1416749598. A
    /// record marked `TSA=insertion` lies one base further right than the chop puts
    /// it, Ensembl's export anchoring it two bases before the database start: `G/GG`
    /// at 5031747 is `-/G` at 5031749-5031748, the cache's rs1987130899, and a
    /// two-allele insertion `A/AT,AG` at 100 is `-/T/G` at 102-101. A record marked
    /// `indel` has no empty allele in the database, and the export adds no anchor to
    /// alleles that share a first base (`VariationFeature.pm:2400-2408`), so its
    /// alleles stay as exported at `POS`: `A/AA` at 5030602 is the cache's `A/AA` at
    /// 5030602-5030602 (rs1370167485), and `AAA/A,AA,AAAA` at 5030913 is
    /// `AAA/A/AA/AAAA` at 5030913-5030915, the cache's ten-base run `AAAAAAAAAA/
    /// AAAAAAAA/AAAAAAAAA/AAAAAAAAAAA` at 5030913-5030922 (rs1232130144) trimmed of
    /// the trailing bases the export drops; both minimise, right end first, to the
    /// cache record's key for a left-aligned input (`-/A` at 5030602; `AA/-` at
    /// 5030913). A `G/GT,A` record whose ALTs disagree on the first base keeps every
    /// allele at `POS`, and a same-length `AC/GT` is never chopped.
    #[test]
    fn indel_records_minimised_to_ensembl_coordinates() {
        let records = records_of(
            "21\t100\tins2\tA\tAT,AG\t.\t.\tdbSNP_156;TSA=insertion\n\
             21\t300\tmixed1\tG\tGT,A\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t400\tmnv1\tAC\tGT\t.\t.\tdbSNP_156;TSA=substitution\n\
             21\t500\tsnv1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t5030218\trs1416749598\tTA\tT\t.\t.\tdbSNP_156;TSA=deletion\n\
             21\t5030602\trs1370167485\tA\tAA\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t5030913\trs1232130144\tAAA\tA,AA,AAAA\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t5031747\trs1987130899\tG\tGG\t.\t.\tdbSNP_156;TSA=insertion\n",
        );
        let by_name: HashMap<&str, &VariationRecord> = records
            .iter()
            .map(|r| (r.variation_name.as_str(), r))
            .collect();
        assert_eq!(span(by_name["rs1416749598"]), (5030219, 5030219, "A/-"));
        assert_eq!(span(by_name["rs1370167485"]), (5030602, 5030602, "A/AA"));
        assert_eq!(
            span(by_name["rs1232130144"]),
            (5030913, 5030915, "AAA/A/AA/AAAA")
        );
        assert_eq!(span(by_name["mixed1"]), (300, 300, "G/GT/A"));
        assert_eq!(span(by_name["mnv1"]), (400, 401, "AC/GT"));
        assert_eq!(span(by_name["rs1987130899"]), (5031749, 5031748, "-/G"));
        assert_eq!(span(by_name["ins2"]), (102, 101, "-/T/G"));
        assert_eq!(span(by_name["snv1"]), (500, 500, "A/G"));
        let starts: Vec<u64> = records.iter().map(|r| r.start).collect();
        assert_eq!(
            starts,
            [102, 300, 400, 500, 5030219, 5030602, 5030913, 5031749]
        );
    }

    /// The export prepends an anchor base to an indel's alleles only when they do
    /// not all share a first base (`VariationFeature.pm:2400-2408`), so a `TSA=indel`
    /// record is chopped only when the chop leaves alleles that neither include an
    /// empty one nor share a first base: the export's `CC/CTTCT` at 5227450 is the
    /// cache's `C/TTCT` at 5227451 (rs1554602202) and `ATG/AC` at 5230561 its `TG/C`
    /// at 5230562-5230563 (rs1555831546), while `A/AGA`, whose chop empties an
    /// allele, and `ACT/ACCG`, whose chop leaves `CT/CCG` sharing `C`, keep the
    /// export's alleles at `POS` (the cache's `AGA/AGAGA` at 15575270-15575272,
    /// rs1361934799, which the export right-trims to `A/AGA`). A `*` allele takes no
    /// part in the comparison, as it takes none in the export's.
    #[test]
    fn indel_records_keep_the_exported_alleles_unless_the_export_anchored_them() {
        let records = records_of(
            "21\t5227450\trs1554602202\tCC\tCTTCT\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t5230561\trs1555831546\tATG\tAC\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t6000000\tstar1\tATG\tAC,*\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t7000000\tshared1\tACT\tACCG\t.\t.\tdbSNP_156;TSA=indel\n\
             21\t15575270\trs1361934799\tA\tAGA\t.\t.\tdbSNP_156;TSA=indel\n",
        );
        let by_name: HashMap<&str, &VariationRecord> = records
            .iter()
            .map(|r| (r.variation_name.as_str(), r))
            .collect();
        assert_eq!(span(by_name["rs1554602202"]), (5227451, 5227451, "C/TTCT"));
        assert_eq!(span(by_name["rs1555831546"]), (5230562, 5230563, "TG/C"));
        assert_eq!(span(by_name["star1"]), (6000001, 6000002, "TG/C/*"));
        assert_eq!(span(by_name["shared1"]), (7000000, 7000002, "ACT/ACCG"));
        assert_eq!(span(by_name["rs1361934799"]), (15575270, 15575270, "A/AGA"));
    }

    /// An allele whose first character is not ASCII never anchors the chop: the
    /// `ÄT/Ä` record keeps both alleles whole with `start` at `POS` rather than
    /// panicking on a byte slice through the character.
    #[test]
    fn non_ascii_allele_is_left_unchopped() {
        let records = records_of("21\t100\trs1\tÄT\tÄ\t.\t.\tdbSNP_156\n");
        assert_eq!(records[0].allele_string, "ÄT/Ä");
        assert_eq!(records[0].start, 100);
    }

    /// The somatic file's records are flagged `somatic` and kept with the span and
    /// allele string the variation database holds for a variant without alleles,
    /// `COSMIC_MUTATION`: on the reference base for a `.` substitution (`C/.` at
    /// 1000 is 1000-1000), after it for a `.` insertion (`T/.` with `TSA=insertion`
    /// at 2000 is 2001-2000), over the bases after the anchor for a `<.>` deletion
    /// (`TTCT/<.>` at 3000 is 3001-3003); a somatic record with alleles keeps them,
    /// an insertion two bases on (`A/ATCC` at 5000 is `-/TCC` at 5002-5001). An
    /// HGMD-PUBLIC record reads `HGMD_MUTATION` with `somatic` 0, two bases on for a
    /// `<.>` insertion (`A/<.>` at 4000 is 4002-4001) and on the deleted base for a
    /// deletion with alleles (`AT/A` at 1500 is 1501-1501); a `.` record of no such
    /// source is skipped. `phenotype_or_disease` and `clin_sig` read the record's
    /// own flags. The identifiers are made up: the shapes are what Ensembl's export
    /// writes for the two sources, whose identifiers the tests do not reproduce.
    #[test]
    fn somatic_vcf_records_are_flagged_and_null_alleles_kept() {
        let dir = tempfile::tempdir().unwrap();
        let germline = write_vcf(
            dir.path(),
            "homo_sapiens-chr21.vcf",
            "21\t100\trs1\tA\t.\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t1500\tCD990002\tAT\tA\t.\t.\tHGMD-PUBLIC_20204;TSA=deletion;E_Phenotype_or_Disease;AA=TA\n\
             21\t4000\tCI990003\tA\t<.>\t.\t.\tHGMD-PUBLIC_20204;TSA=insertion;E_Phenotype_or_Disease\n",
        );
        let somatic = write_vcf(
            dir.path(),
            "homo_sapiens_somatic.vcf",
            "1\t100\tCOSV9900001\tT\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease;AA=T\n\
             21\t1000\tCOSV9900002\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n\
             21\t2000\tCOSV9900003\tT\t.\t.\t.\tCOSMIC_102;TSA=insertion;E_Phenotype_or_Disease\n\
             21\t3000\tCOSV9900004\tTTCT\t<.>\t.\t.\tCOSMIC_102;TSA=deletion;E_Phenotype_or_Disease\n\
             21\t5000\trs2517045861\tA\tATCC\t.\t.\tClinVar_202509;TSA=insertion;E_Phenotype_or_Disease;CLIN_uncertain_significance\n",
        );
        let inputs = GenomeWideInputs {
            somatic: Some(somatic),
            ..Default::default()
        };
        let records = collect(&germline, &inputs, None);
        let summary: Vec<String> = records
            .iter()
            .map(|r| {
                format!(
                    "{} {}-{} {} somatic={} pheno={} clin_sig={}",
                    r.variation_name,
                    r.start,
                    r.end,
                    r.allele_string,
                    r.somatic,
                    r.phenotype_or_disease,
                    r.clin_sig.as_deref().unwrap_or("-")
                )
            })
            .collect();
        assert_eq!(
            summary,
            [
                "COSV9900002 1000-1000 COSMIC_MUTATION somatic=1 pheno=1 clin_sig=-",
                "CD990002 1501-1501 HGMD_MUTATION somatic=0 pheno=1 clin_sig=-",
                "COSV9900003 2001-2000 COSMIC_MUTATION somatic=1 pheno=1 clin_sig=-",
                "COSV9900004 3001-3003 COSMIC_MUTATION somatic=1 pheno=1 clin_sig=-",
                "CI990003 4002-4001 HGMD_MUTATION somatic=0 pheno=1 clin_sig=-",
                "rs2517045861 5002-5001 -/TCC somatic=1 pheno=1 clin_sig=uncertain_significance",
            ]
        );
    }

    /// `clin_sig` is the record's `CLIN_<term>` flags in their order, joined by `,`,
    /// the string `Database/Variation.pm:177` reads from the database (its spaces
    /// already `_`); `phenotype_or_disease` is 1 when the record carries
    /// `E_Phenotype_or_Disease` (`:178`); a record with neither leaves both unset.
    #[test]
    fn clin_sig_and_phenotype_read_the_records_flags() {
        let records = records_of(
            "21\t100\trs1\tC\tA,T\t.\t.\tdbSNP_156;TSA=SNV;E_Freq;E_Phenotype_or_Disease;CLIN_uncertain_significance;CLIN_likely_pathogenic;AA=C\n\
             21\t200\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;E_Freq;CLIN_pathogenic\n\
             21\t300\trs3\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;E_Freq\n",
        );
        assert_eq!(
            records[0].clin_sig.as_deref(),
            Some("uncertain_significance,likely_pathogenic")
        );
        assert_eq!(records[0].phenotype_or_disease, 1);
        assert_eq!(records[1].clin_sig.as_deref(), Some("pathogenic"));
        assert_eq!(records[1].phenotype_or_disease, 0);
        assert_eq!(records[2].clin_sig, None);
        assert_eq!(records[2].phenotype_or_disease, 0);
    }

    /// Records whose FILTER is neither `PASS`, `.` nor empty are skipped; the first
    /// of several IDs names the record and a `.` ID leaves the name empty.
    #[test]
    fn filtered_records_are_skipped_and_the_first_id_names_a_record() {
        let records = records_of(
            "21\t100\trs1;rs2\tA\tG\t.\tPASS\tdbSNP_156\n\
             21\t200\trs3\tA\tG\t.\tLowQual\tdbSNP_156\n\
             21\t300\t.\tA\tG\t.\t.\tdbSNP_156\n",
        );
        let names: Vec<&str> = records.iter().map(|r| r.variation_name.as_str()).collect();
        assert_eq!(names, ["rs1", ""]);
    }

    /// The 1000 Genomes file gives rs2838519's A allele `AFR=0.7148 AMR=0.379
    /// EAS=0.4534 EUR=0.6024 SAS=0.5583`; the release-116 cache holds the record
    /// as `A:0.5625 A:0.7148 A:0.379 A:0.4534 A:0.6024 A:0.5583` under `AF, AFR,
    /// AMR, EAS, EUR, SAS`. The whole-panel frequency is the recovered counts
    /// (945 + 263 + 457 + 606 + 546 = 2,817 of 5,008 chromosomes) printed `%.4g`;
    /// the C allele, absent from the frequency file, gets no entry. rs142566111's
    /// A allele, a singleton (`AFR=0.0008`, one of 1,322 chromosomes), reads
    /// `0.0001997`, the cache's figure for one allele in 5,008, with `A:0` in the
    /// populations that lack it; its G allele, 0 in every population, gets no
    /// entry, as the cache has none. A variant the file does not carry has no
    /// frequency column, and a somatic variant is given none
    /// (`Dumper/Variation.pm:280`).
    #[test]
    fn thousand_genomes_file_populates_af_and_continental_columns() {
        let dir = tempfile::tempdir().unwrap();
        let germline = write_vcf(
            dir.path(),
            "homo_sapiens-chr21.vcf",
            "21\t44195140\trs2838519\tG\tA,C\t.\t.\tdbSNP_156;TSA=SNV;E_Freq\n\
             21\t44229058\trs142566111\tT\tA,G\t.\t.\tdbSNP_156;TSA=SNV;E_Freq\n\
             21\t44300000\trs9\tC\tT\t.\t.\tdbSNP_156;TSA=SNV\n",
        );
        let frequencies = write_vcf(
            dir.path(),
            "1000GENOMES-phase_3.vcf",
            "21\t44195140\trs2838519\tG\tA\t.\t.\tdbSNP_156;TSA=SNV;E_Freq;E_1000G;AFR=0.7148;AMR=0.379;EAS=0.4534;EUR=0.6024;SAS=0.5583\n\
             21\t44229058\trs142566111\tT\tA,G\t.\t.\tdbSNP_156;TSA=SNV;E_Freq;E_1000G;AFR=0.0008,0;AMR=0,0;EAS=0,0;EUR=0,0;SAS=0,0\n\
             21\t44400000\tCOSV1\tA\tT\t.\t.\tCOSMIC_102;TSA=SNV;AFR=0.5;AMR=0.5;EAS=0.5;EUR=0.5;SAS=0.5\n",
        );
        let somatic = write_vcf(
            dir.path(),
            "homo_sapiens_somatic.vcf",
            "21\t44400000\tCOSV1\tA\tT\t.\t.\tCOSMIC_102;TSA=SNV\n",
        );
        let inputs = GenomeWideInputs {
            frequencies: Some(frequencies),
            somatic: Some(somatic),
            ..Default::default()
        };
        let records = collect(&germline, &inputs, None);
        let frequencies = |name: &str| -> BTreeMap<&str, &str> {
            let record = records.iter().find(|r| r.variation_name == name).unwrap();
            record
                .frequencies
                .iter()
                .map(|(column, text)| (column.as_str(), text.as_str()))
                .collect()
        };
        assert_eq!(
            frequencies("rs2838519"),
            BTreeMap::from([
                ("AF", "A:0.5625"),
                ("AFR", "A:0.7148"),
                ("AMR", "A:0.379"),
                ("EAS", "A:0.4534"),
                ("EUR", "A:0.6024"),
                ("SAS", "A:0.5583"),
            ])
        );
        assert_eq!(
            frequencies("rs142566111"),
            BTreeMap::from([
                ("AF", "A:0.0001997"),
                ("AFR", "A:0.0008"),
                ("AMR", "A:0"),
                ("EAS", "A:0"),
                ("EUR", "A:0"),
                ("SAS", "A:0"),
            ])
        );
        assert!(frequencies("rs9").is_empty());
        assert!(frequencies("COSV1").is_empty());
    }

    /// A frequency VCF allele matches a variant's allele through the minimised
    /// REF/ALT pair, so a record the frequency file spells with another anchor or
    /// a subset of the alleles still joins: `G/GT` at 300 (the frequency file) and
    /// `G/GT,A` at 300 (the variant, unchopped because `A` shares no first base)
    /// both minimise to `-_T_301`. A pair an allele of more than one base belongs
    /// to is keyed both ways: `AAA/A` at 5030913 trimmed from the left is `AA/-`
    /// at 5030914, from the right `AA/-` at 5030913; `GT/G` at 200 has one key,
    /// nothing being shared on the right; `ACA/ACACA` at 100 is `-/CA` at 103 from
    /// the left and `-/AC` at 100 from the right.
    #[test]
    fn frequency_alleles_match_through_minimised_pairs() {
        assert_eq!(minimised_keys(300, "G", "GT"), ["-_T_301"]);
        assert_eq!(minimised_keys(300, "G", "A"), ["G_A_300"]);
        assert_eq!(minimised_keys(200, "GT", "G"), ["T_-_201"]);
        assert_eq!(minimised_keys(10178, "C", "CC"), ["-_C_10179", "-_C_10178"]);
        assert_eq!(
            minimised_keys(5030913, "AAA", "A"),
            ["AA_-_5030914", "AA_-_5030913"]
        );
        assert_eq!(
            minimised_keys(100, "ACA", "ACACA"),
            ["-_CA_103", "-_AC_100"]
        );

        let dir = tempfile::tempdir().unwrap();
        let germline = write_vcf(
            dir.path(),
            "homo_sapiens-chr21.vcf",
            "21\t300\tmixed1\tG\tGT,A\t.\t.\tdbSNP_156;TSA=indel\n",
        );
        let frequencies = write_vcf(
            dir.path(),
            "1000GENOMES-phase_3.vcf",
            "21\t300\tmixed1\tG\tGT\t.\t.\tdbSNP_156;TSA=insertion;AFR=0.1;AMR=0.2;EAS=0.3;EUR=0.4;SAS=0.5\n",
        );
        let inputs = GenomeWideInputs {
            frequencies: Some(frequencies),
            ..Default::default()
        };
        let records = collect(&germline, &inputs, None);
        assert_eq!(records[0].allele_string, "G/GT/A");
        assert_eq!(records[0].frequencies["AFR"], "GT:0.1");
        assert_eq!(records[0].frequencies["SAS"], "GT:0.5");
        // 132 + 139 + 302 + 402 + 489 of 5,008.
        assert_eq!(records[0].frequencies["AF"], "GT:0.2923");
    }

    /// `%.4g` as Perl prints it: four significant digits, trailing zeros gone,
    /// exponent form below 0.0001 and from 10,000.
    #[test]
    fn perl_g4_prints_four_significant_digits() {
        assert_eq!(perl_g4(0.0), "0");
        assert_eq!(perl_g4(0.5625), "0.5625");
        assert_eq!(perl_g4(0.379), "0.379");
        assert_eq!(perl_g4(1.0 / 5008.0), "0.0001997");
        assert_eq!(perl_g4(0.0000121), "1.21e-05");
        assert_eq!(perl_g4(1.0), "1");
        assert_eq!(perl_g4(12345.0), "1.234e+04");
    }

    /// A record of the phenotype- or clinically-associated file that the
    /// per-chromosome file already holds (same name, position and alleles) is
    /// written once; one it lacks is added, sorted into place by start; one the
    /// somatic file also holds is written once, flagged `somatic`, whether or
    /// not the per-chromosome file holds a copy too.
    #[test]
    fn genome_wide_germline_records_are_added_once() {
        let dir = tempfile::tempdir().unwrap();
        let germline = write_vcf(
            dir.path(),
            "homo_sapiens-chr21.vcf",
            "21\t100\trs1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;E_Phenotype_or_Disease;CLIN_pathogenic\n\
             21\t250\tCOSV9\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n\
             21\t300\trs3\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
        );
        let phenotype = write_vcf(
            dir.path(),
            "homo_sapiens_phenotype_associated.vcf",
            "21\t100\trs1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;E_Phenotype_or_Disease;CLIN_pathogenic\n\
             21\t200\trs2516155147\tG\tA\t.\t.\tClinVar_202509;TSA=SNV;E_Phenotype_or_Disease;CLIN_uncertain_significance\n\
             21\t250\tCOSV9\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n\
             22\t200\trs22\tG\tA\t.\t.\tClinVar_202509;TSA=SNV;E_Phenotype_or_Disease\n",
        );
        let clinical = write_vcf(
            dir.path(),
            "homo_sapiens_clinically_associated.vcf",
            "21\t200\trs2516155147\tG\tA\t.\t.\tClinVar_202509;TSA=SNV;E_Phenotype_or_Disease;CLIN_uncertain_significance\n",
        );
        let somatic = write_vcf(
            dir.path(),
            "homo_sapiens_somatic.vcf",
            "21\t250\tCOSV9\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n",
        );
        let inputs = GenomeWideInputs {
            germline: vec![phenotype, clinical],
            somatic: Some(somatic),
            ..Default::default()
        };
        let records = collect(&germline, &inputs, None);
        let names: Vec<(&str, u64, u8)> = records
            .iter()
            .map(|r| (r.variation_name.as_str(), r.start, r.somatic))
            .collect();
        assert_eq!(
            names,
            [
                ("rs1", 100, 0),
                ("rs2516155147", 200, 0),
                ("COSV9", 250, 1),
                ("rs3", 300, 0)
            ]
        );
        assert_eq!(
            records[1].clin_sig.as_deref(),
            Some("uncertain_significance")
        );
    }

    /// Write each table's tab-separated rows to `<table>.txt.gz` under `dir`,
    /// gzipped as the Ensembl FTP publishes the dumps, with `\N` for NULL.
    fn write_dump_dir(dir: &Path, tables: &[(&str, &str)]) {
        for (table, rows) in tables {
            let path = dir.join(format!("{table}.txt.gz"));
            let file = fs::File::create(&path).unwrap();
            let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
            gz.write_all(rows.as_bytes()).unwrap();
            gz.finish().unwrap();
        }
    }

    /// Variation dump rows for one ClinVar span: 21:100-100 carries two germline
    /// features of rs1 (risk alleles T and C, significances `pathogenic` and
    /// `uncertain significance`, reference allele A), one germline feature of another
    /// variation on the same span, one germline feature whose risk allele is set twice
    /// (A, then C), and four somatic features (a somatic classification with a comma,
    /// an oncogenicity, a classification on a phenotype whose description carries
    /// MySQL-escaped backslashes, and an oncogenicity on a `Gene`-type feature of the
    /// span, which the adaptor's location query joins like any other); on the
    /// database span 21:200-204 of rs5 (held untrimmed as `AGAGA/AGA`), one germline
    /// feature (risk allele `GA`, `benign`, reference allele `AGAGA`) and one somatic
    /// feature; and on 21:201-202, the span rs5's export trims to, one somatic
    /// feature of rs7.
    const DUMP_TABLES: &[(&str, &str)] = &[
        ("source", "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n32\tClinVar\t202509\t\\N\t\\N\t\\N\tmixed\tphenotype_feature\n"),
        ("seq_region", "131537\t21\t4\n131538\t22\t4\n"),
        ("attrib_type", "480\tclinvar_clin_sig\tClinVar clinical significance\t\\N\n483\trisk_allele\tRisk allele\t\\N\n577\tsomatic_clin_sig\tClinVar somatic classification\t\\N\n578\toncogenic_clin_sig\tClinVar somatic classification of oncogenicity\t\\N\n583\tclinvar_ref_allele\tClinVar reference allele\t\\N\n"),
        ("phenotype", "7\t\\N\t\\N\tBreast cancer\t\\N\n8\t\\N\t\\N\tLung carcinoma\t\\N\n9\t\\N\t\\N\tNot provided\t\\N\n10\t\\N\t\\N\tCarcinoma\\\\of\\\\lung\t\\N\n"),
        ("phenotype_feature", "\
1001\t7\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tGermline\n\
1002\t9\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\t\\N\n\
1003\t7\t32\t\\N\tVariation\trs9\t1\t131537\t100\t100\t1\tGermline\n\
1004\t7\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tGermline\n\
1005\t8\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tSomatic\n\
1006\t7\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tSomatic\n\
1007\t7\t1\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tGermline\n\
1008\t7\t32\t\\N\tGene\tENSG1\t1\t131537\t100\t100\t1\tGermline\n\
1009\t10\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tSomatic\n\
1010\t9\t32\t\\N\tGene\tENSG1\t1\t131537\t100\t100\t1\tSomatic\n\
1011\t7\t32\t\\N\tVariation\trs5\t1\t131537\t200\t204\t1\tGermline\n\
1012\t8\t32\t\\N\tVariation\trs5\t1\t131537\t200\t204\t1\tSomatic\n\
1013\t9\t32\t\\N\tVariation\trs7\t1\t131537\t201\t202\t1\tSomatic\n"),
        ("phenotype_feature_attrib", "\
1001\t480\tpathogenic\n1001\t483\tT\n1001\t583\tA\n\
1002\t480\tuncertain significance\n1002\t483\tC\n\
1003\t480\tbenign\n1003\t483\tG\n\
1004\t480\tnot provided\n\
1005\t577\tTier I - Strong,Tier II - Potential\n\
1006\t578\tOncogenic\n\
1007\t480\tbenign\n1007\t483\tG\n\
1009\t577\tTier II - Potential\n1004\t483\tA\n1004\t483\tC\n\
1010\t578\tLikely oncogenic\n\
1011\t480\tbenign\n1011\t483\tGA\n1011\t583\tAGAGA\n\
1012\t577\tTier I - Strong\n\
1013\t577\tTier III - Unknown\n"),
    ];

    /// `clin_sig_allele` joins each `risk_allele:clinvar_clin_sig` pair of the germline
    /// ClinVar features that carry the record's own name on its span, spaces in the
    /// term as `_`, and `clin_sig_ref_allele` is their `clinvar_ref_allele`
    /// (`Database/Variation.pm:196-214`; the pair grammar is what Ensembl VEP reports
    /// for rs113488022, `T:pathogenic;...;C:uncertain_significance`). A feature of
    /// another variation on the same span and a dbSNP-source feature are left out; a
    /// feature whose risk allele is set twice keeps the last row
    /// (`PhenotypeFeatureAdaptor.pm:1273-1276`), so rs1's `not provided` pairs with C.
    /// `clinical_impact` is `_format_clinical_impact` (`:252-275`) over the somatic
    /// features: `phenotype(classification)` with `,` read `;`, then
    /// `phenotype(oncogenicity:term)`, joined by `;`, the dump's `\\` read as one
    /// backslash. A record on another span gets none. rs5, exported trimmed (`AGA/A`
    /// at 200, `GA/-` at 201-202) where the tables hold the database's 200-204,
    /// still takes the features of that span, found through its own somatic feature,
    /// not those of rs7, whose database span is the 201-202 the export gives rs5 and
    /// whose own record, the same `AGA/A` at 200, takes them.
    #[test]
    fn clinical_impact_from_clinvar_dumps() {
        let dir = tempfile::tempdir().unwrap();
        write_dump_dir(dir.path(), DUMP_TABLES);
        let clinvar = ClinVarAttribs::load(dir.path(), &[]).unwrap();
        assert!(!clinvar.is_empty());

        let vcf = write_vcf(
            dir.path(),
            "homo_sapiens-chr21.vcf",
            "21\t100\trs1\tA\tC,T\t.\t.\tdbSNP_156;TSA=SNV;E_Phenotype_or_Disease;CLIN_pathogenic\n\
             21\t150\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t200\trs5\tAGA\tA\t.\t.\tdbSNP_156;TSA=deletion;E_Phenotype_or_Disease;CLIN_benign\n\
             21\t200\trs7\tAGA\tA\t.\t.\tdbSNP_156;TSA=deletion;E_Phenotype_or_Disease\n",
        );
        let records = collect(&vcf, &GenomeWideInputs::default(), Some(&clinvar));
        assert_eq!(records.len(), 4);

        let r = &records[0];
        assert_eq!(
            r.clin_sig_allele,
            Some("T:pathogenic;C:uncertain_significance;C:not_provided".to_string())
        );
        assert_eq!(r.clin_sig_ref_allele, Some("A".to_string()));
        assert_eq!(
            r.clinical_impact,
            Some("Lung carcinoma(Tier I - Strong;Tier II - Potential);Breast cancer(oncogenicity:Oncogenic);Carcinoma\\of\\lung(Tier II - Potential);Not provided(oncogenicity:Likely oncogenic)".to_string())
        );

        let other = &records[1];
        assert_eq!(other.clin_sig_allele, None);
        assert_eq!(other.clin_sig_ref_allele, None);
        assert_eq!(other.clinical_impact, None);

        let trimmed = &records[2];
        assert_eq!(trimmed.variation_name, "rs5");
        assert_eq!(span(trimmed), (201, 202, "GA/-"));
        assert_eq!(trimmed.clin_sig_allele.as_deref(), Some("GA:benign"));
        assert_eq!(trimmed.clin_sig_ref_allele.as_deref(), Some("AGAGA"));
        assert_eq!(
            trimmed.clinical_impact.as_deref(),
            Some("Lung carcinoma(Tier I - Strong)")
        );

        let same_export = &records[3];
        assert_eq!(same_export.variation_name, "rs7");
        assert_eq!(span(same_export), (201, 202, "GA/-"));
        assert_eq!(same_export.clin_sig_allele, None);
        assert_eq!(
            same_export.clinical_impact.as_deref(),
            Some("Not provided(Tier III - Unknown)")
        );
    }

    /// The dump reader accepts the uncompressed `.txt` spelling, keeps nothing from a
    /// chromosome outside the list it is given, and reports an empty table set when
    /// the source table names no ClinVar.
    #[test]
    fn clinvar_dump_reader_accepts_plain_text_and_absent_clinvar() {
        let dir = tempfile::tempdir().unwrap();
        for (table, rows) in DUMP_TABLES {
            fs::write(dir.path().join(format!("{table}.txt")), rows).unwrap();
        }
        let clinvar = ClinVarAttribs::load(dir.path(), &[]).unwrap();
        assert!(!clinvar.is_empty());
        let other_chromosome = ClinVarAttribs::load(dir.path(), &["22".to_string()]).unwrap();
        assert!(other_chromosome.is_empty());

        let no_clinvar = tempfile::tempdir().unwrap();
        fs::write(
            no_clinvar.path().join("source.txt"),
            "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n",
        )
        .unwrap();
        assert!(ClinVarAttribs::load(no_clinvar.path(), &[])
            .unwrap()
            .is_empty());
    }

    /// A dump directory without the tables is an error naming the missing table.
    #[test]
    fn clinvar_dump_reader_names_the_missing_table() {
        let dir = tempfile::tempdir().unwrap();
        let err = ClinVarAttribs::load(dir.path(), &[]).unwrap_err();
        assert!(err.to_string().contains("source.txt"), "{err}");
    }

    /// A chromosome that cannot name a directory of the cache (`..`, a name with
    /// a separator) is refused before the variations directory is created, so a
    /// chromosome read from a file name never writes outside `variations/`.
    #[test]
    fn a_chromosome_that_is_not_a_path_component_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let columns = variation_columns(false);
        for chr in ["..", ".", "21/22", "", "21\\22"] {
            let error = match AllVarsWriter::create(dir.path(), chr, &columns) {
                Ok(_) => panic!("{chr:?} was accepted"),
                Err(error) => error,
            };
            assert_eq!(
                error.to_string(),
                format!("invalid chromosome name '{chr}'"),
                "{chr:?}"
            );
        }
        assert!(!dir.path().join("variations").exists());
    }

    /// The file is written as the dumper prints a record (`Dumper/Variation.pm:
    /// 233-264`): the chromosome first, one cell per column, `.` for an empty one,
    /// `failed`, `somatic` and `phenotype_or_disease` only when 1, `end` only when
    /// it differs from `start`, `strand` only when -1, and a frequency column's
    /// text as joined. A record starting before the one written before it is
    /// refused, the index being built in file order.
    #[test]
    fn all_vars_lines_follow_the_dumpers_cells() {
        let dir = tempfile::tempdir().unwrap();
        let somatic = VariationRecord {
            variation_name: "COSV1".to_string(),
            somatic: 1,
            start: 100,
            end: 100,
            strand: 1,
            allele_string: COSMIC_MUTATION.to_string(),
            phenotype_or_disease: 1,
            ..Default::default()
        };
        let snv = VariationRecord {
            variation_name: "rs1".to_string(),
            start: 150,
            end: 150,
            strand: 1,
            allele_string: "A/G".to_string(),
            clin_sig: Some("pathogenic".to_string()),
            clin_sig_allele: Some("G:pathogenic".to_string()),
            clin_sig_ref_allele: Some("A".to_string()),
            clinical_impact: Some("Breast cancer(Tier I - Strong)".to_string()),
            ..Default::default()
        };
        let mut insertion = VariationRecord {
            variation_name: "rs2".to_string(),
            start: 201,
            end: 200,
            strand: 1,
            allele_string: "-/G".to_string(),
            ..Default::default()
        };
        insertion
            .frequencies
            .insert("AF".to_string(), "G:0.1".to_string());
        let columns = variation_columns(true);
        assert_eq!(columns.len(), 15 + 6);
        let mut writer = AllVarsWriter::create(dir.path(), "21", &columns).unwrap();
        for record in [&somatic, &snv, &insertion] {
            writer.write(record).unwrap();
        }
        let error = writer.write(&snv).unwrap_err();
        assert!(
            error.to_string().contains("not in position order"),
            "{error}"
        );
        assert_eq!(writer.finish().unwrap(), 3);

        let path = dir.path().join("variations/21/all_vars.gz");
        let mut text = String::new();
        bgzf::io::Reader::new(fs::File::open(&path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(
            text,
            "21\tCOSV1\t.\t1\t100\t.\tCOSMIC_MUTATION\t.\t.\t1\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n\
             21\trs1\t.\t.\t150\t.\tA/G\t.\tpathogenic\t.\tG:pathogenic\tA\tBreast cancer(Tier I - Strong)\t.\t.\t.\t.\t.\t.\t.\t.\n\
             21\trs2\t.\t.\t201\t200\t-/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.1\t.\t.\t.\t.\t.\n"
        );
        assert!(dir.path().join("variations/21/all_vars.gz.csi").is_file());
    }

    /// A per-chromosome file out of position order fails the write with the
    /// out-of-order error and leaves no variation file or index behind.
    #[test]
    fn failed_write_leaves_no_variation_file() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = write_vcf(
            dir.path(),
            "homo_sapiens-chr21.vcf",
            "21\t500\trs5\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t600\trs6\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t100\trs1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
        );
        let out = dir.path().join("cache");
        let error = write_chromosome(
            "21",
            &vcf,
            &GenomeWideInputs::default(),
            None,
            &variation_columns(false),
            &out,
        )
        .unwrap_err();
        assert!(
            error.to_string().contains("not in position order"),
            "{error}"
        );
        let written: Vec<_> = fs::read_dir(out.join("variations/21"))
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect();
        assert!(written.is_empty(), "{written:?}");
    }

    /// The per-chromosome file's records come out in start order although the
    /// file is in position order: an insertion at 100 (start 102) follows the
    /// substitutions at 100 and 101, and a deletion at 102 (start 103) the
    /// insertion. Records sharing a start come out in name order, the dumper's
    /// (`sort -k5,5n` falls back to the whole line): at 102, `ins` before
    /// `rs10` before `rs9`, whatever the file's order.
    #[test]
    fn records_are_emitted_in_start_order() {
        let records = records_of(
            "21\t100\tins\tA\tAT\t.\t.\tdbSNP_156;TSA=insertion\n\
             21\t100\tsnv1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t101\tsnv2\tC\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t102\tdel\tCA\tC\t.\t.\tdbSNP_156;TSA=deletion\n\
             21\t102\trs9\tC\tT\t.\t.\tdbSNP_156;TSA=SNV\n\
             21\t102\trs10\tC\tA\t.\t.\tdbSNP_156;TSA=SNV\n",
        );
        let order: Vec<(&str, u64)> = records
            .iter()
            .map(|r| (r.variation_name.as_str(), r.start))
            .collect();
        assert_eq!(
            order,
            [
                ("snv1", 100),
                ("snv2", 101),
                ("ins", 102),
                ("rs10", 102),
                ("rs9", 102),
                ("del", 103)
            ]
        );
    }

    /// A `.csi` without the tabix header names no sequence, so a chromosome
    /// cannot be queried through it: the file is read through instead and only
    /// the chromosome's records are kept.
    #[test]
    fn index_without_a_header_falls_back_to_reading_through() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("homo_sapiens_somatic.vcf.gz");
        let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        writer.write_all(ENSEMBL_HEADER.as_bytes()).unwrap();
        let mut indexer = Indexer::<BinnedIndex>::new(14, 5);
        let records = [
            (0, 100, "1\t100\tCOSV1\tA\t.\t.\t.\tCOSMIC_102;TSA=SNV"),
            (1, 200, "21\t200\tCOSV2\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV"),
        ];
        for (sequence, pos, line) in records {
            let begin = writer.virtual_position();
            writeln!(writer, "{line}").unwrap();
            let end = writer.virtual_position();
            let pos = Position::try_from(pos).unwrap();
            indexer
                .add_record(Some((sequence, pos, pos, true)), Chunk::new(begin, end))
                .unwrap();
        }
        writer.try_finish().unwrap();
        noodles::csi::fs::write(format!("{}.csi", path.display()), &indexer.build(2)).unwrap();
        assert!(has_index(&path));

        let mut names = Vec::new();
        for_each_record(&path, Some("21"), |fields| {
            names.push(fields.name);
            Ok(())
        })
        .unwrap();
        assert_eq!(names, ["COSV2"]);
    }
}
