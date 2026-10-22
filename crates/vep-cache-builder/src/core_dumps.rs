// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Ensembl core and variation table dumps: every transcript fact the
//! GFF3 does not carry.
//!
//! Ensembl publishes each MySQL database as one tab-separated file per table
//! (`<table>.txt.gz`, or `<table>.txt` where a release leaves it
//! uncompressed), with `\N` for NULL and the MySQL escapes `\t`, `\n`, `\r`,
//! `\0`, `\b`, `\Z` and `\\` inside a value. The column order is the
//! database schema's (`sql/table.sql` of the core API and of the variation
//! API).
//!
//! From the core database: `gene` (display xref, canonical transcript),
//! `transcript` (database identifier to stable identifier),
//! `transcript_attrib` with `attrib_type` (the attributes Ensembl VEP's
//! dumper keeps), `translation` and `translation_attrib` (the ENSP, its
//! version and the sequence edits), `protein_feature` with
//! `analysis_description` (the domains and their display labels), `xref`,
//! `object_xref` and `external_db` (CCDS, RefSeq, UniProt, UniParc, HGNC and
//! the symbol source), `seq_region` with `seq_region_attrib` (the codon
//! table) and `meta` (the source versions). From the variation database:
//! `phenotype_feature` and `phenotype` (the gene phenotype count).
//!
//! Each fact is read the way `AnnotationSource::Database::Transcript` reads
//! it through the core API before the dumper writes the cache:
//! `prefetch_gene_ids`, `prefetch_transcript_ids`,
//! `prefetch_translation_ids` and `prefetch_translation_data`.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use tracing::{info, warn};

/// The longest dump line read, 16 MiB. A row's widest values are an xref
/// description and a protein feature's alignment string, each well under a
/// megabyte, so the cap refuses a file with no line end in it long before the
/// row buffer grows to its size.
const MAX_LINE_BYTES: usize = 16 << 20;

/// The lines of a text file, each without its `\n`, read one at a time with
/// the length of a line capped: a line past `cap` bytes is an error naming
/// `display` and the line number, as is a read failure, and either ends the
/// iteration.
pub(crate) struct BoundedLines<R> {
    reader: R,
    cap: usize,
    display: String,
    line_number: u64,
    done: bool,
}

impl<R: BufRead> BoundedLines<R> {
    pub(crate) fn new(reader: R, cap: usize, display: String) -> Self {
        Self {
            reader,
            cap,
            display,
            line_number: 0,
            done: false,
        }
    }
}

impl<R: BufRead> Iterator for BoundedLines<R> {
    type Item = Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        // The cap plus the line end: a line of exactly the cap is read whole.
        let limit = self.cap as u64 + 1;
        let mut line = Vec::new();
        self.line_number += 1;
        let item = match (&mut self.reader).take(limit).read_until(b'\n', &mut line) {
            Ok(0) => None,
            Ok(read) if read as u64 == limit && line.last() != Some(&b'\n') => {
                Some(Err(anyhow::anyhow!(
                    "{} line {} is longer than the {}-byte limit for one line",
                    self.display,
                    self.line_number,
                    self.cap
                )))
            }
            Ok(_) => {
                if line.last() == Some(&b'\n') {
                    line.pop();
                }
                Some(Ok(line))
            }
            Err(e) => Some(Err(
                anyhow::Error::from(e).context(format!("Failed to read {}", self.display))
            )),
        };
        self.done = !matches!(item, Some(Ok(_)));
        item
    }
}

/// The core tables `load` reads, plus `seq_region_synonym`, which the
/// chromosome synonyms file is written from (`info::write_chr_synonyms`).
pub const CORE_TABLES: [&str; 15] = [
    "external_db",
    "attrib_type",
    "seq_region",
    "seq_region_attrib",
    "seq_region_synonym",
    "meta",
    "gene",
    "transcript",
    "xref",
    "object_xref",
    "translation",
    "translation_attrib",
    "transcript_attrib",
    "analysis_description",
    "protein_feature",
];

/// The variation tables read for the gene phenotype count.
pub const VARIATION_TABLES: [&str; 2] = ["phenotype", "phenotype_feature"];

/// The transcript attribute codes the Ensembl VEP dumper keeps in its cache
/// (`Pipeline/DumpVEP/Dumper/Core.pm` `clean_transcript`).
pub const KEPT_ATTRIBUTE_CODES: [&str; 22] = [
    "gencode_basic",
    "gencode_primary",
    "miRNA",
    "ncRNA",
    "cds_start_NF",
    "cds_end_NF",
    "TSL",
    "appris",
    "rseq_mrna_match",
    "rseq_mrna_nonmatch",
    "rseq_5p_mismatch",
    "rseq_cds_mismatch",
    "rseq_3p_mismatch",
    "rseq_nctran_mismatch",
    "rseq_no_comparison",
    "rseq_ens_match_wt",
    "rseq_ens_match_cds",
    "rseq_ens_no_match",
    "enst_refseq_compare",
    "_rna_edit",
    "MANE_Select",
    "MANE_Plus_Clinical",
];

/// The translation attribute codes `Translation::get_all_SeqEdits` reads, in
/// its order.
const SEQ_EDIT_CODES: [&str; 4] = [
    "initial_met",
    "_selenocysteine",
    "amino_acid_sub",
    "_stop_codon_rt",
];

/// A transcript attribute as the cache stores it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attribute {
    pub code: String,
    pub value: String,
}

/// One protein feature of a translation, the four fields the dumper keeps.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinFeature {
    /// Protein coordinates (`seq_start`, `seq_end`).
    pub start: u64,
    pub end: u64,
    /// `hit_name`, the feature's accession in its source database.
    pub hseqname: String,
    /// The analysis display label (`analysis_description.display_label`).
    pub display_label: Option<String>,
}

/// A translation sequence edit as `Bio::EnsEMBL::SeqEdit` reads the
/// attribute: the value is `<start> <end> <alt_seq>`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqEdit {
    pub code: String,
    pub value: String,
}

/// What the core database knows about a transcript's translation.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TranslationFacts {
    pub stable_id: String,
    pub version: Option<u32>,
    /// `Uniprot/SWISSPROT` xref labels, in object xref row order.
    pub swissprot: Vec<String>,
    /// `Uniprot/SPTREMBL` xref labels.
    pub trembl: Vec<String>,
    /// `UniParc` xref labels.
    pub uniparc: Vec<String>,
    /// `Uniprot_isoform` xref labels.
    pub uniprot_isoform: Vec<String>,
    pub protein_features: Vec<ProteinFeature>,
    pub seq_edits: Vec<SeqEdit>,
}

/// What the core database knows about a transcript beyond its GFF3 line.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TranscriptFacts {
    /// The kept attributes, in table order, `cds_start_NF` before `cds_end_NF`.
    pub attributes: Vec<Attribute>,
    /// The first `CCDS` xref label.
    pub ccds: Option<String>,
    /// Every `RefSeq_mRNA` xref label, in object xref row order.
    pub refseq: Vec<String>,
    /// Whether the gene's `canonical_transcript_id` is this transcript.
    pub is_canonical: bool,
    /// Whether the transcript carries a `readthrough_tra` attribute; the
    /// Ensembl VEP dumper drops such a transcript.
    pub readthrough: bool,
    pub translation: Option<TranslationFacts>,
}

/// What the core and variation databases know about a gene.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GeneFacts {
    /// The display xref's label, or the `RefSeq_gene_name` xref's label when
    /// the gene has no display xref.
    pub symbol: Option<String>,
    /// The display xref's database name.
    pub symbol_source: Option<String>,
    /// The display xref's primary accession when its database is HGNC.
    pub hgnc_id: Option<String>,
    /// Significant phenotype features of type `Gene` on this gene, or `None`
    /// when no variation tables were read.
    pub phenotype_count: Option<u32>,
}

/// The facts of one Ensembl release, keyed by stable identifier.
#[derive(Debug, Default)]
pub struct CoreDumps {
    transcripts: HashMap<String, TranscriptFacts>,
    genes: HashMap<String, GeneFacts>,
    /// Sequence region name to codon table, for regions carrying the
    /// `codon_table` attribute.
    codon_tables: HashMap<String, u8>,
    /// `meta_key` to `meta_value`.
    meta: HashMap<String, String>,
}

impl CoreDumps {
    /// Read the tables under `dir`; the variation tables are optional and
    /// their absence leaves every gene's phenotype count `None`.
    pub fn load(dir: &Path) -> Result<CoreDumps> {
        let mut dumps = CoreDumps::default();

        let external_db = id_to_name(dir, "external_db", 10)?;
        let attrib_type = id_to_name(dir, "attrib_type", 4)?;

        let seq_region = id_to_name(dir, "seq_region", 4)?;
        for row in read_table(dir, "seq_region_attrib", 3)? {
            let row = row?;
            if attrib_type.get(row.field(1)).map(String::as_str) != Some("codon_table") {
                continue;
            }
            let Some(name) = seq_region.get(row.field(0)) else {
                continue;
            };
            let codon_table = u8::try_from(row.number(2)?)
                .with_context(|| format!("codon table {} of {name}", row.field(2)))?;
            dumps.codon_tables.insert(name.clone(), codon_table);
        }
        for row in read_table(dir, "meta", 4)? {
            let row = row?;
            dumps
                .meta
                .insert(row.field(2).to_string(), row.field(3).to_string());
        }

        struct GeneRow {
            stable_id: String,
            display_xref_id: Option<String>,
            canonical_transcript_id: Option<String>,
        }
        let mut gene_rows: HashMap<String, GeneRow> = HashMap::new();
        for row in read_table(dir, "gene", 16)? {
            let row = row?;
            gene_rows.insert(
                row.field(0).to_string(),
                GeneRow {
                    stable_id: row.field(12).to_string(),
                    display_xref_id: row.optional(7).map(str::to_string),
                    canonical_transcript_id: row.optional(11).map(str::to_string),
                },
            );
        }

        struct TranscriptRow {
            stable_id: String,
            gene_id: String,
            canonical_translation_id: Option<String>,
        }
        let mut transcript_rows: HashMap<String, TranscriptRow> = HashMap::new();
        for row in read_table(dir, "transcript", 17)? {
            let row = row?;
            transcript_rows.insert(
                row.field(0).to_string(),
                TranscriptRow {
                    stable_id: row.field(13).to_string(),
                    gene_id: row.field(1).to_string(),
                    canonical_translation_id: row.optional(12).map(str::to_string),
                },
            );
        }

        // Xrefs: the databases the prefetches read, plus every gene's
        // display xref.
        const WANTED_DBS: [&str; 7] = [
            "CCDS",
            "RefSeq_mRNA",
            "RefSeq_gene_name",
            "Uniprot/SWISSPROT",
            "Uniprot/SPTREMBL",
            "UniParc",
            "Uniprot_isoform",
        ];
        let display_xref_ids: HashSet<&str> = gene_rows
            .values()
            .filter_map(|g| g.display_xref_id.as_deref())
            .collect();
        struct Xref {
            db_name: String,
            primary_acc: String,
            display_label: String,
        }
        let mut xrefs: HashMap<String, Xref> = HashMap::new();
        for row in read_table(dir, "xref", 8)? {
            let row = row?;
            let db_name = external_db.get(row.field(1)).cloned().unwrap_or_default();
            if WANTED_DBS.contains(&db_name.as_str()) || display_xref_ids.contains(row.field(0)) {
                xrefs.insert(
                    row.field(0).to_string(),
                    Xref {
                        db_name,
                        primary_acc: row.field(2).to_string(),
                        display_label: row.field(3).to_string(),
                    },
                );
            }
        }

        // Object xrefs per object, in the dump's line order, each xref once:
        // the order the core API returns a feature's DBEntries in. Its query
        // has no ORDER BY, so the rows come as stored, which is neither the
        // object_xref nor the xref identifier order.
        let mut transcript_xrefs: HashMap<String, Vec<String>> = HashMap::new();
        let mut translation_xrefs: HashMap<String, Vec<String>> = HashMap::new();
        let mut gene_xrefs: HashMap<String, Vec<String>> = HashMap::new();
        for row in read_table(dir, "object_xref", 6)? {
            let row = row?;
            let xref_id = row.field(3);
            if !xrefs.contains_key(xref_id) {
                continue;
            }
            let target = match row.field(2) {
                "Transcript" => &mut transcript_xrefs,
                "Translation" => &mut translation_xrefs,
                "Gene" => &mut gene_xrefs,
                _ => continue,
            };
            target
                .entry(row.field(1).to_string())
                .or_default()
                .push(xref_id.to_string());
        }
        let labels = |xref_ids: Option<&Vec<String>>, db_name: &str| -> Vec<String> {
            let mut seen = HashSet::new();
            xref_ids
                .into_iter()
                .flatten()
                .filter(|id| seen.insert(id.as_str()))
                .filter_map(|id| xrefs.get(id))
                .filter(|x| x.db_name == db_name)
                .map(|x| x.display_label.clone())
                .collect()
        };

        struct TranslationRow {
            transcript_id: String,
            stable_id: String,
            version: Option<u32>,
        }
        let mut translation_rows: HashMap<String, TranslationRow> = HashMap::new();
        for row in read_table(dir, "translation", 10)? {
            let row = row?;
            let version = if row.optional(7).is_some() {
                Some(u32::try_from(row.number(7)?).with_context(|| {
                    format!("translation version {} of {}", row.field(7), row.field(6))
                })?)
            } else {
                None
            };
            translation_rows.insert(
                row.field(0).to_string(),
                TranslationRow {
                    transcript_id: row.field(1).to_string(),
                    stable_id: row.field(6).to_string(),
                    version,
                },
            );
        }
        let mut seq_edits: HashMap<String, Vec<SeqEdit>> = HashMap::new();
        for row in read_table(dir, "translation_attrib", 3)? {
            let row = row?;
            let Some(code) = attrib_type.get(row.field(1)) else {
                continue;
            };
            if SEQ_EDIT_CODES.contains(&code.as_str()) {
                seq_edits
                    .entry(row.field(0).to_string())
                    .or_default()
                    .push(SeqEdit {
                        code: code.clone(),
                        value: row.field(2).to_string(),
                    });
            }
        }

        // Transcript attributes, the dumper's kept codes only, and the
        // readthrough marker it drops a transcript for.
        let mut attributes: HashMap<String, Vec<Attribute>> = HashMap::new();
        let mut readthrough: HashSet<String> = HashSet::new();
        for row in read_table(dir, "transcript_attrib", 3)? {
            let row = row?;
            let Some(code) = attrib_type.get(row.field(1)) else {
                continue;
            };
            if code == "readthrough_tra" {
                readthrough.insert(row.field(0).to_string());
            }
            if KEPT_ATTRIBUTE_CODES.contains(&code.as_str()) {
                attributes
                    .entry(row.field(0).to_string())
                    .or_default()
                    .push(Attribute {
                        code: code.clone(),
                        value: row.field(2).to_string(),
                    });
            }
        }
        for attrs in attributes.values_mut() {
            order_cds_flags(attrs);
        }

        let display_labels: HashMap<String, Option<String>> =
            read_table(dir, "analysis_description", 5)?
                .map(|r| {
                    let r = r?;
                    Ok((r.field(0).to_string(), r.optional(2).map(str::to_string)))
                })
                .collect::<Result<_>>()?;
        let mut protein_features: HashMap<String, Vec<(u64, ProteinFeature)>> = HashMap::new();
        for row in read_table(dir, "protein_feature", 15)? {
            let row = row?;
            let feature = ProteinFeature {
                start: row.number(2)?,
                end: row.number(3)?,
                hseqname: row.field(6).to_string(),
                display_label: display_labels.get(row.field(7)).cloned().flatten(),
            };
            protein_features
                .entry(row.field(1).to_string())
                .or_default()
                .push((row.number(0)?, feature));
        }

        // A transcript's translation is the one its canonical_translation_id
        // names (TranslationAdaptor::fetch_by_Transcript joins on it); any
        // other translation row of the transcript is not its translation.
        let mut translations: HashMap<String, TranslationFacts> = HashMap::new();
        for (translation_id, tl) in &translation_rows {
            let canonical = transcript_rows
                .get(&tl.transcript_id)
                .and_then(|tr| tr.canonical_translation_id.as_deref());
            if canonical != Some(translation_id.as_str()) {
                continue;
            }
            let mut features = protein_features.remove(translation_id).unwrap_or_default();
            features.sort_by_key(|(id, _)| *id);
            let mut edits: Vec<SeqEdit> = Vec::new();
            if let Some(found) = seq_edits.get(translation_id) {
                for code in SEQ_EDIT_CODES {
                    edits.extend(found.iter().filter(|e| e.code == code).cloned());
                }
            }
            let xref_ids = translation_xrefs.get(translation_id);
            translations.insert(
                tl.transcript_id.clone(),
                TranslationFacts {
                    stable_id: tl.stable_id.clone(),
                    version: tl.version,
                    swissprot: labels(xref_ids, "Uniprot/SWISSPROT"),
                    trembl: labels(xref_ids, "Uniprot/SPTREMBL"),
                    uniparc: labels(xref_ids, "UniParc"),
                    uniprot_isoform: labels(xref_ids, "Uniprot_isoform"),
                    protein_features: features.into_iter().map(|(_, f)| f).collect(),
                    seq_edits: edits,
                },
            );
        }

        // Genes: symbol, source, HGNC id (`prefetch_gene_ids`).
        for (gene_id, g) in &gene_rows {
            let mut facts = GeneFacts::default();
            if let Some(xref) = g.display_xref_id.as_ref().and_then(|id| xrefs.get(id)) {
                facts.symbol = Some(xref.display_label.clone());
                facts.symbol_source = Some(xref.db_name.clone());
                if xref.db_name == "HGNC" {
                    facts.hgnc_id = Some(xref.primary_acc.clone());
                }
            } else {
                facts.symbol = labels(gene_xrefs.get(gene_id), "RefSeq_gene_name")
                    .into_iter()
                    .next();
            }
            dumps.genes.insert(g.stable_id.clone(), facts);
        }

        // Transcripts: attributes, xrefs, canonical flag and translation by
        // stable id (`prefetch_transcript_ids`).
        let canonical: HashSet<&str> = gene_rows
            .values()
            .filter_map(|g| g.canonical_transcript_id.as_deref())
            .collect();
        for (transcript_id, tr) in &transcript_rows {
            if !gene_rows.contains_key(&tr.gene_id) {
                continue;
            }
            let xref_ids = transcript_xrefs.get(transcript_id);
            dumps.transcripts.insert(
                tr.stable_id.clone(),
                TranscriptFacts {
                    attributes: attributes.remove(transcript_id).unwrap_or_default(),
                    ccds: labels(xref_ids, "CCDS").into_iter().next(),
                    refseq: labels(xref_ids, "RefSeq_mRNA"),
                    is_canonical: canonical.contains(transcript_id.as_str()),
                    readthrough: readthrough.contains(transcript_id),
                    translation: translations.remove(transcript_id),
                },
            );
        }

        info!(
            "Core dumps: {} transcripts, {} genes, {} codon-table regions",
            dumps.transcripts.len(),
            dumps.genes.len(),
            dumps.codon_tables.len()
        );

        dumps.load_phenotypes(dir)?;
        Ok(dumps)
    }

    /// The gene phenotype count, as `PhenotypeFeatureAdaptor::fetch_all_by_Gene`
    /// counts it: significant `phenotype_feature` rows of type `Gene` whose
    /// `object_id` is the gene's stable id and whose phenotype has a class.
    fn load_phenotypes(&mut self, dir: &Path) -> Result<()> {
        if VARIATION_TABLES
            .iter()
            .any(|t| table_path(dir, t).is_none())
        {
            warn!(
                "No phenotype tables under {}: gene_phenotype is not written",
                dir.display()
            );
            return Ok(());
        }
        let mut classed: HashSet<String> = HashSet::new();
        for row in read_table(dir, "phenotype", 5)? {
            let row = row?;
            if row.optional(4).is_some() {
                classed.insert(row.field(0).to_string());
            }
        }
        let mut counts: HashMap<String, u32> = HashMap::new();
        for row in read_table(dir, "phenotype_feature", 12)? {
            let row = row?;
            if row.field(4) == "Gene"
                && row.field(6) == "1"
                && row.optional(1).is_some_and(|p| classed.contains(p))
            {
                *counts.entry(row.field(5).to_string()).or_default() += 1;
            }
        }
        for (stable_id, facts) in self.genes.iter_mut() {
            facts.phenotype_count = Some(counts.get(stable_id).copied().unwrap_or(0));
        }
        info!(
            "Phenotype features: {} genes with a significant phenotype",
            counts.len()
        );
        Ok(())
    }

    pub fn transcript(&self, stable_id: &str) -> Option<&TranscriptFacts> {
        self.transcripts.get(stable_id)
    }

    pub fn gene(&self, stable_id: &str) -> Option<&GeneFacts> {
        self.genes.get(stable_id)
    }

    /// The codon table of a sequence region, when it carries the attribute
    /// (table 2 on `MT`).
    pub fn codon_table(&self, seq_region_name: &str) -> Option<u8> {
        self.codon_tables.get(seq_region_name).copied()
    }

    /// The core source versions Ensembl VEP writes to its cache info
    /// (`Database/Transcript.pm` `info`): `assembly.name`, `gencode.version`,
    /// `genebuild.initial_release_date` and `genebuild.version`, each under
    /// the key before its first dot, a later key overwriting an earlier one.
    pub fn source_versions(&self) -> BTreeMap<String, String> {
        let mut out = BTreeMap::new();
        for key in [
            "assembly.name",
            "gencode.version",
            "genebuild.initial_release_date",
            "genebuild.version",
        ] {
            if let Some(value) = self.meta.get(key) {
                let short = key.split('.').next().unwrap_or(key);
                out.insert(short.to_string(), value.clone());
            }
        }
        out
    }
}

/// Puts `cds_start_NF` before `cds_end_NF` when a transcript carries both,
/// the other attributes staying where the dump lists them. FLAGS prints the
/// two in list order (OutputFactory.pm:1477); the Ensembl cache holds the order
/// its database query returned (Dumper/Core.pm:203-207), start first on most
/// transcripts and end first on a few, and the dumps do not carry that query
/// order, so one fixed order is written.
fn order_cds_flags(attributes: &mut [Attribute]) {
    let position = |code: &str| attributes.iter().position(|a| a.code == code);
    if let (Some(start), Some(end)) = (position("cds_start_NF"), position("cds_end_NF")) {
        if end < start {
            attributes.swap(start, end);
        }
    }
}

/// A table's first two columns as a map, its identifier to its name:
/// `external_db`, `attrib_type` and `seq_region` all start that way.
pub fn id_to_name(
    dir: &Path,
    table: &'static str,
    columns: usize,
) -> Result<HashMap<String, String>> {
    read_table(dir, table, columns)?
        .map(|row| {
            let row = row?;
            Ok((row.field(0).to_string(), row.field(1).to_string()))
        })
        .collect()
}

/// `<dir>/<table>.txt.gz`, else `<dir>/<table>.txt`, when one exists.
pub fn table_path(dir: &Path, table: &str) -> Option<PathBuf> {
    [format!("{table}.txt.gz"), format!("{table}.txt")]
        .into_iter()
        .map(|name| dir.join(name))
        .find(|p| p.is_file())
}

/// One dump row: its unescaped fields, `\N` kept as the NULL marker.
pub struct Row {
    table: &'static str,
    fields: Vec<String>,
}

impl Row {
    /// The field at `index`, or an empty string past the row's end.
    pub fn field(&self, index: usize) -> &str {
        self.fields.get(index).map(String::as_str).unwrap_or("")
    }

    /// The field at `index`, `None` for NULL (`\N`) or a missing column.
    pub fn optional(&self, index: usize) -> Option<&str> {
        match self.fields.get(index) {
            Some(v) if v != "\\N" => Some(v.as_str()),
            _ => None,
        }
    }

    /// The field at `index` as an unsigned integer; a value that is not one
    /// is an error naming the table and column.
    pub fn number(&self, index: usize) -> Result<u64> {
        let value = self.field(index);
        value.parse().with_context(|| {
            format!(
                "{} column {index} holds {value:?}, not a number",
                self.table
            )
        })
    }
}

/// The rows of a table dump, each line split on tabs with MySQL's escapes
/// decoded (a literal `\N` field stays as written, see [`Row::optional`]);
/// a row with fewer than `columns` fields, the schema's column count, is an
/// error naming the table and line.
pub fn read_table(
    dir: &Path,
    table: &'static str,
    columns: usize,
) -> Result<impl Iterator<Item = Result<Row>>> {
    let Some(path) = table_path(dir, table) else {
        bail!(
            "Table dump {}.txt.gz (or {}.txt) not found under {}",
            table,
            table,
            dir.display()
        );
    };
    let file =
        std::fs::File::open(&path).with_context(|| format!("Cannot open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::with_capacity(1 << 20, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };
    let lines = BoundedLines::new(reader, MAX_LINE_BYTES, path.display().to_string());
    Ok(lines.enumerate().map(move |(index, line)| {
        let fields = split_escaped(&String::from_utf8_lossy(&line?));
        if fields.len() < columns {
            bail!(
                "{table} line {} has {} columns, the schema has {columns}",
                index + 1,
                fields.len()
            );
        }
        Ok(Row { table, fields })
    }))
}

/// Split one dump line on its unescaped tabs and decode the MySQL escapes
/// in each field; the NULL marker `\N` is left as the two characters.
fn split_escaped(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars();
    while let Some(c) = chars.next() {
        match c {
            '\t' => fields.push(std::mem::take(&mut field)),
            '\\' => match chars.next() {
                Some('t') => field.push('\t'),
                Some('n') => field.push('\n'),
                Some('r') => field.push('\r'),
                Some('0') => field.push('\0'),
                Some('b') => field.push('\u{8}'),
                Some('Z') => field.push('\u{1a}'),
                Some('\\') => field.push('\\'),
                Some('N') => field.push_str("\\N"),
                Some(other) => field.push(other),
                None => field.push('\\'),
            },
            _ => field.push(c),
        }
    }
    fields.push(field);
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Writes `rows` as `<table>.txt` under `dir`.
    fn write_table(dir: &Path, table: &str, rows: &[&str]) {
        let mut text = rows.join("\n");
        text.push('\n');
        std::fs::write(dir.join(format!("{table}.txt")), text).unwrap();
    }

    /// A core dump slice: chromosome 21's TPTE (ENSG00000274391, display
    /// xref HGNC:12023), its transcript ENST00000622113 with translation
    /// ENSP00000482040, a lncRNA gene with no display xref, a gene whose
    /// display xref is RFAM, one whose display xref is EntrezGene, and MT.
    /// Rows follow the release 116 GRCh38 core schema column order.
    fn fixture(dir: &Path) {
        write_table(
            dir,
            "external_db",
            &[
                "1100\tHGNC\t\\N\tKNOWNXREF\t5\tHGNC Symbol\tMISC\t\\N\t\\N\t\\N",
                "1300\tEntrezGene\t\\N\tKNOWNXREF\t5\tNCBI gene (formerly Entrezgene)\tMISC\t\\N\t\\N\t\\N",
                "4200\tRFAM\t\\N\tXREF\t5\tRFAM\tMISC\t\\N\t\\N\t\\N",
                "3800\tCCDS\t\\N\tKNOWNXREF\t5\tCCDS\tMISC\t\\N\t\\N\t\\N",
                "1801\tRefSeq_mRNA\t\\N\tKNOWNXREF\t5\tRefSeq mRNA\tMISC\t\\N\t\\N\t\\N",
                "2200\tUniprot/SWISSPROT\t\\N\tKNOWNXREF\t5\tUniProtKB/Swiss-Prot\tMISC\t\\N\t\\N\t\\N",
                "2000\tUniprot/SPTREMBL\t\\N\tKNOWNXREF\t5\tUniProtKB/TrEMBL\tMISC\t\\N\t\\N\t\\N",
                "2300\tUniParc\t\\N\tKNOWNXREF\t5\tUniParc\tMISC\t\\N\t\\N\t\\N",
                "2250\tUniprot_isoform\t\\N\tKNOWNXREF\t5\tUniProtKB isoform\tMISC\t\\N\t\\N\t\\N",
                "1810\tRefSeq_gene_name\t\\N\tXREF\t5\tRefSeq gene name\tMISC\t\\N\t\\N\t\\N",
            ],
        );
        write_table(
            dir,
            "attrib_type",
            &[
                "1\tcodon_table\tCodon Table\tAlternate codon table",
                "2\tgencode_basic\tGENCODE basic annotation\tGENCODE basic transcript",
                "3\tTSL\tTranscript Support Level\t",
                "4\tappris\tAPPRIS\t",
                "5\tMANE_Select\tMANE Select\t",
                "6\tcds_start_NF\tCDS start not found\t",
                "7\tcds_end_NF\tCDS end not found\t",
                "8\tmiRNA\tMicro RNA\t",
                "9\tncRNA\tStructure\t",
                "10\tremark\tRemark\t",
                "11\tinitial_met\tInitial methionine\t",
                "12\t_selenocysteine\tSelenocysteine\t",
                "13\tgencode_primary\tGENCODE primary annotation\t",
                "14\tMANE_Plus_Clinical\tMANE Plus Clinical\t",
                "15\tMolecularWeight\tMolecular weight\t",
                "16\treadthrough_tra\treadthrough transcript\t",
            ],
        );
        write_table(
            dir,
            "seq_region",
            &["131539\t21\t4\t46709983", "132907\tMT\t4\t16569"],
        );
        write_table(dir, "seq_region_attrib", &["132907\t1\t2"]);
        write_table(
            dir,
            "meta",
            &[
                "1\t\\N\tschema_version\t116",
                "21\t1\tassembly.name\tGRCh38.p14",
                "300\t1\tgencode.version\tGENCODE 50",
                "16\t1\tgenebuild.last_geneset_update\t2025-11",
                "6410\t1\tgenebuild.version\tGENCODE50",
            ],
        );
        // gene: gene_id biotype analysis_id seq_region_id start end strand
        // display_xref_id source description is_current canonical_transcript_id
        // stable_id version created_date modified_date
        write_table(
            dir,
            "gene",
            &[
                "100\tprotein_coding\t1\t131539\t10521553\t10606140\t1\t9001\tensembl_havana\ttransmembrane phosphatase with tensin homology\t1\t200\tENSG00000274391\t6\t\\N\t\\N",
                "101\tlncRNA\t1\t131539\t5232668\t5243903\t-1\t\\N\thavana\tnovel transcript\t1\t201\tENSG00000279669\t2\t\\N\t\\N",
                "102\tmisc_RNA\t1\t131539\t8212572\t8212679\t-1\t9002\tensembl\t\\N\t1\t202\tENSG00000223254\t1\t\\N\t\\N",
                "103\tunprocessed_pseudogene\t1\t131539\t9000000\t9001000\t1\t9003\thavana\t\\N\t1\t203\tENSG00000290541\t1\t\\N\t\\N",
                "104\tMt_tRNA\t1\t132907\t577\t647\t1\t9004\tinsdc\tmitochondrially encoded tRNA-Phe\t1\t204\tENSG00000210049\t1\t\\N\t\\N",
                "105\tlncRNA\t1\t131539\t6000000\t6001000\t1\t\\N\thavana\t\\N\t1\t205\tENSG00000999999\t1\t\\N\t\\N",
            ],
        );
        // transcript: transcript_id gene_id analysis_id seq_region_id start end
        // strand display_xref_id source biotype description is_current
        // canonical_translation_id stable_id version created_date modified_date
        write_table(
            dir,
            "transcript",
            &[
                "200\t100\t1\t131539\t10521553\t10606140\t1\t\\N\tensembl_havana\tprotein_coding\t\\N\t1\t300\tENST00000622113\t4\t\\N\t\\N",
                "210\t100\t1\t131539\t10521553\t10600000\t1\t\\N\thavana\tprotein_coding\t\\N\t1\t\\N\tENST00000000210\t1\t\\N\t\\N",
                "201\t101\t1\t131539\t5232668\t5243903\t-1\t\\N\thavana_tagene\tlncRNA\t\\N\t1\t\\N\tENST00000623753\t2\t\\N\t\\N",
                "202\t102\t1\t131539\t8212572\t8212679\t-1\t\\N\tensembl\tmisc_RNA\t\\N\t1\t\\N\tENST00000411223\t1\t\\N\t\\N",
                "203\t103\t1\t131539\t9000000\t9001000\t1\t\\N\thavana\tunprocessed_pseudogene\t\\N\t1\t\\N\tENST00000000203\t1\t\\N\t\\N",
                "204\t104\t1\t132907\t577\t647\t1\t\\N\tinsdc\tMt_tRNA\t\\N\t1\t\\N\tENST00000387314\t1\t\\N\t\\N",
                "205\t105\t1\t131539\t6000000\t6001000\t1\t\\N\thavana\tlncRNA\t\\N\t1\t\\N\tENST00000000205\t1\t\\N\t\\N",
            ],
        );
        // xref: xref_id external_db_id dbprimary_acc display_label version
        // description info_type info_text
        write_table(
            dir,
            "xref",
            &[
                "9001\t1100\tHGNC:12023\tTPTE\t0\t\\N\tDIRECT\t",
                "9002\t4200\tRF00019\tY_RNA\t0\t\\N\tDIRECT\t",
                "9003\t1300\t100132340\tCYP4F29P\t0\t\\N\tDEPENDENT\t",
                "9004\t1100\tHGNC:7481\tMT-TF\t0\t\\N\tDIRECT\t",
                "9010\t3800\tCCDS74772\tCCDS74772.1\t0\t\\N\tDIRECT\t",
                "9011\t1801\tNM_199261\tNM_199261.4\t0\t\\N\tDIRECT\t",
                "9012\t1801\tNM_001290224\tNM_001290224.2\t0\t\\N\tDIRECT\t",
                "9020\t2200\tP56180\tP56180.197\t0\t\\N\tDIRECT\t",
                "9021\t2000\tA0A0G2JMT5\tA0A0G2JMT5.47\t0\t\\N\tDIRECT\t",
                "9022\t2300\tUPI000050C1B9\tUPI000050C1B9\t0\t\\N\tCHECKSUM\t",
                "9023\t2250\tP56180-1\tP56180-1\t0\t\\N\tDIRECT\t",
                "9024\t2250\tP56180-2\tP56180-2\t0\t\\N\tDIRECT\t",
                "9030\t1810\tLOC105372760\tLOC105372760\t0\t\\N\tDIRECT\t",
            ],
        );
        // object_xref: object_xref_id ensembl_id ensembl_object_type xref_id
        // linkage_annotation analysis_id
        write_table(
            dir,
            "object_xref",
            &[
                "1\t200\tTranscript\t9010\t\\N\t1",
                "3\t200\tTranscript\t9012\t\\N\t1",
                "2\t200\tTranscript\t9011\t\\N\t1",
                "4\t300\tTranslation\t9020\t\\N\t1",
                "5\t300\tTranslation\t9021\t\\N\t1",
                "6\t300\tTranslation\t9022\t\\N\t1",
                "8\t300\tTranslation\t9024\t\\N\t1",
                "7\t300\tTranslation\t9023\t\\N\t1",
                "9\t300\tTranslation\t9023\t\\N\t2",
                "10\t105\tGene\t9030\t\\N\t1",
            ],
        );
        // translation: translation_id transcript_id seq_start start_exon_id
        // seq_end end_exon_id stable_id version created_date modified_date
        write_table(
            dir,
            "translation",
            &["300\t200\t44\t5001\t54\t5009\tENSP00000482040\t1\t\\N\t\\N"],
        );
        write_table(
            dir,
            "translation_attrib",
            &["300\t12\t523 523 U", "300\t15\t51204.3", "300\t11\t1 1 M"],
        );
        write_table(
            dir,
            "transcript_attrib",
            &[
                "200\t2\tGENCODE basic",
                "200\t3\ttsl1 (assigned to previous version 3)",
                "200\t10\tLEAP_5",
                "200\t4\tprincipal2",
                "200\t5\tNM_199261.4",
                "200\t13\tGENCODE primary",
                "210\t6\t1",
                "210\t7\t0",
                "202\t9\t1:104\\t.11(3.(9.4(.(3.4)3.).)5.",
                "202\t8\t43-63",
                "204\t3\ttslNA",
                "205\t16\treadthrough",
            ],
        );
        write_table(
            dir,
            "analysis_description",
            &[
                "8633\tProtein features from Pfam\tPfam\t1\t\\N",
                "8650\tProtein features from Gene3D\tGene3D\t1\t\\N",
            ],
        );
        // protein_feature: protein_feature_id translation_id seq_start seq_end
        // hit_start hit_end hit_name analysis_id ...
        write_table(
            dir,
            "protein_feature",
            &[
                "70002\t300\t200\t350\t1\t150\tG3DSA:3.90.190.10\t8650\t\\N\t\\N\t\\N\t\t\t\t\\N",
                "70001\t300\t187\t334\t1\t148\tPF00102\t8633\t\\N\t\\N\t\\N\t\t\t\t\\N",
                "70003\t300\t10\t20\t1\t11\tunknown.1\t9999\t\\N\t\\N\t\\N\t\t\t\t\\N",
            ],
        );
    }

    /// The variation slice: phenotype 59727 has a class, 60000 has none;
    /// TPTE has four significant gene phenotype features and one that is
    /// not significant, one on a phenotype without a class, and one
    /// variation feature; MT-TF has one.
    fn phenotype_fixture(dir: &Path) {
        write_table(
            dir,
            "phenotype",
            &[
                "59727\t\\N\tHYPOPARATHYROIDISM\thypoparathyroidism\t663",
                "60000\t\\N\t\\N\tunclassified\t\\N",
            ],
        );
        // phenotype_feature: phenotype_feature_id phenotype_id source_id
        // study_id type object_id is_significant seq_region_id start end strand DNA_type
        write_table(
            dir,
            "phenotype_feature",
            &[
                "1\t59727\t36\t\\N\tGene\tENSG00000274391\t1\t131539\t10521553\t10606140\t1\t\\N",
                "2\t59727\t36\t\\N\tGene\tENSG00000274391\t1\t131539\t10521553\t10606140\t1\t\\N",
                "3\t59727\t36\t\\N\tGene\tENSG00000274391\t1\t131539\t10521553\t10606140\t1\t\\N",
                "4\t59727\t36\t\\N\tGene\tENSG00000274391\t1\t131539\t10521553\t10606140\t1\t\\N",
                "5\t59727\t36\t\\N\tGene\tENSG00000274391\t0\t131539\t10521553\t10606140\t1\t\\N",
                "6\t60000\t36\t\\N\tGene\tENSG00000274391\t1\t131539\t10521553\t10606140\t1\t\\N",
                "7\t59727\t36\t\\N\tVariation\trs1234\t1\t131539\t10521553\t10521553\t1\t\\N",
                "8\t59727\t36\t\\N\tGene\tENSG00000210049\t1\t132907\t577\t647\t1\t\\N",
            ],
        );
    }

    fn load_fixture(with_phenotypes: bool) -> (tempfile::TempDir, CoreDumps) {
        let dir = tempfile::tempdir().unwrap();
        fixture(dir.path());
        if with_phenotypes {
            phenotype_fixture(dir.path());
        }
        let dumps = CoreDumps::load(dir.path()).unwrap();
        (dir, dumps)
    }

    /// The HGNC id is the display xref's primary accession when the xref's
    /// database is HGNC (`Database/Transcript.pm` `prefetch_gene_ids`:
    /// `_gene_hgnc_id = $xref->primary_id if $xref->dbname eq 'HGNC'`);
    /// a gene whose display xref is another database has none, and a gene
    /// without a display xref takes its symbol from a `RefSeq_gene_name`
    /// xref and has no source.
    #[test]
    fn hgnc_id_from_display_xref() {
        let (_dir, dumps) = load_fixture(false);
        let tpte = dumps.gene("ENSG00000274391").unwrap();
        assert_eq!(tpte.symbol.as_deref(), Some("TPTE"));
        assert_eq!(tpte.symbol_source.as_deref(), Some("HGNC"));
        assert_eq!(tpte.hgnc_id.as_deref(), Some("HGNC:12023"));
        assert_eq!(tpte.phenotype_count, None);

        let lnc = dumps.gene("ENSG00000279669").unwrap();
        assert_eq!(lnc.symbol, None);
        assert_eq!(lnc.symbol_source, None);
        assert_eq!(lnc.hgnc_id, None);

        let refseq_named = dumps.gene("ENSG00000999999").unwrap();
        assert_eq!(refseq_named.symbol.as_deref(), Some("LOC105372760"));
        assert_eq!(refseq_named.symbol_source, None);
        assert_eq!(refseq_named.hgnc_id, None);
    }

    /// The symbol source is the display xref's `external_db.db_name`
    /// (`prefetch_gene_ids`: `_gene_symbol_source = $xref->dbname`), so an
    /// RFAM display xref gives `RFAM` and an EntrezGene one `EntrezGene`,
    /// neither with an HGNC id.
    #[test]
    fn symbol_source_from_xref_external_db() {
        let (_dir, dumps) = load_fixture(false);
        let y_rna = dumps.gene("ENSG00000223254").unwrap();
        assert_eq!(y_rna.symbol.as_deref(), Some("Y_RNA"));
        assert_eq!(y_rna.symbol_source.as_deref(), Some("RFAM"));
        assert_eq!(y_rna.hgnc_id, None);
    }

    /// Ensembl VEP's golden output lists `RFAM` and `EntrezGene` beside
    /// `HGNC` as symbol sources; both come from the display xref's database.
    #[test]
    fn symbol_source_rfam_and_entrezgene() {
        let (_dir, dumps) = load_fixture(false);
        let pseudo = dumps.gene("ENSG00000290541").unwrap();
        assert_eq!(pseudo.symbol.as_deref(), Some("CYP4F29P"));
        assert_eq!(pseudo.symbol_source.as_deref(), Some("EntrezGene"));
        assert_eq!(pseudo.hgnc_id, None);
        assert_eq!(
            dumps
                .gene("ENSG00000223254")
                .unwrap()
                .symbol_source
                .as_deref(),
            Some("RFAM")
        );
    }

    /// `_gene_phenotype` is the count `PhenotypeFeatureAdaptor::fetch_all_by_Gene`
    /// returns (`Database/Transcript.pm` `get_features_by_regions_uncached`):
    /// rows of type `Gene` on the gene's stable id, `is_significant = 1`,
    /// joined to a phenotype whose `class_attrib_id` is set. TPTE counts
    /// 4 (not the insignificant row, the unclassed phenotype or the
    /// variation row); MT-TF counts 1; a gene with no row counts 0.
    #[test]
    fn gene_phenotype_count_from_phenotype_feature() {
        let (_dir, dumps) = load_fixture(true);
        assert_eq!(
            dumps.gene("ENSG00000274391").unwrap().phenotype_count,
            Some(4)
        );
        assert_eq!(
            dumps.gene("ENSG00000210049").unwrap().phenotype_count,
            Some(1)
        );
        assert_eq!(
            dumps.gene("ENSG00000279669").unwrap().phenotype_count,
            Some(0)
        );
    }

    /// A `readthrough_tra` attribute marks the transcript the Ensembl VEP
    /// dumper drops (`Database/Transcript.pm` `get_features_by_regions_uncached`);
    /// the GRCh37 release-87 GFF3 carries no tag for it, so the dumps are its
    /// source there.
    #[test]
    fn readthrough_marker_from_transcript_attrib() {
        let (_dir, dumps) = load_fixture(false);
        assert!(dumps.transcript("ENST00000000205").unwrap().readthrough);
        assert!(!dumps.transcript("ENST00000622113").unwrap().readthrough);
    }

    /// Without the variation tables no count is known: `None`, not 0.
    #[test]
    fn gene_phenotype_unknown_without_variation_tables() {
        let (_dir, dumps) = load_fixture(false);
        assert_eq!(dumps.gene("ENSG00000274391").unwrap().phenotype_count, None);
    }

    /// `_uniparc` is the `UniParc` xrefs of the translation
    /// (`prefetch_translation_ids`, `$tl->get_all_DBEntries`), by display
    /// label.
    #[test]
    fn uniparc_from_translation_xref() {
        let (_dir, dumps) = load_fixture(false);
        let tl = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .translation
            .as_ref()
            .unwrap();
        assert_eq!(tl.stable_id, "ENSP00000482040");
        assert_eq!(tl.version, Some(1));
        assert_eq!(tl.uniparc, vec!["UPI000050C1B9"]);
        // The isoform xrefs come in object xref row order, each xref once
        // (xref 9023 is linked twice).
        assert_eq!(tl.uniprot_isoform, vec!["P56180-2", "P56180-1"]);
    }

    /// Only the translation `transcript.canonical_translation_id` names is
    /// the transcript's (`TranslationAdaptor::fetch_by_Transcript` joins
    /// `translation_id` to it); another translation row of the transcript, or
    /// a translation row of a transcript whose column is NULL, is not.
    #[test]
    fn translation_follows_canonical_translation_id() {
        let dir = tempfile::tempdir().unwrap();
        fixture(dir.path());
        let mut rows = std::fs::read_to_string(dir.path().join("translation.txt")).unwrap();
        rows.push_str("301\t200\t44\t5001\t54\t5009\tENSP00000000301\t1\t\\N\t\\N\n");
        rows.push_str("302\t210\t44\t5001\t54\t5009\tENSP00000000302\t1\t\\N\t\\N\n");
        std::fs::write(dir.path().join("translation.txt"), rows).unwrap();
        let dumps = CoreDumps::load(dir.path()).unwrap();
        let tl = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .translation
            .as_ref();
        assert_eq!(tl.map(|t| t.stable_id.as_str()), Some("ENSP00000482040"));
        assert!(dumps
            .transcript("ENST00000000210")
            .unwrap()
            .translation
            .is_none());
    }

    /// A NULL `translation.version` loads as no version.
    #[test]
    fn null_translation_version_is_none() {
        let dir = tempfile::tempdir().unwrap();
        fixture(dir.path());
        let rows = std::fs::read_to_string(dir.path().join("translation.txt")).unwrap();
        let rows = rows.replace("ENSP00000482040\t1\t", "ENSP00000482040\t\\N\t");
        assert!(rows.contains("ENSP00000482040\t\\N\t"));
        std::fs::write(dir.path().join("translation.txt"), rows).unwrap();
        let dumps = CoreDumps::load(dir.path()).unwrap();
        let tl = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .translation
            .as_ref()
            .unwrap();
        assert_eq!(tl.version, None);
    }

    /// `_swissprot` and `_trembl` are the `Uniprot/SWISSPROT` and
    /// `Uniprot/SPTREMBL` xrefs of the translation, by display label, which
    /// carries the UniProt entry version (`P56180.197`).
    #[test]
    fn swissprot_trembl_from_translation_xrefs() {
        let (_dir, dumps) = load_fixture(false);
        let tl = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .translation
            .as_ref()
            .unwrap();
        assert_eq!(tl.swissprot, vec!["P56180.197"]);
        assert_eq!(tl.trembl, vec!["A0A0G2JMT5.47"]);
        assert!(dumps
            .transcript("ENST00000623753")
            .unwrap()
            .translation
            .is_none());
    }

    /// `_ccds` is the first `CCDS` xref's display label and `_refseq` every
    /// `RefSeq_mRNA` label (`prefetch_transcript_ids`), in the object xref
    /// table's row order (xref 9012's row precedes 9011's), the order the
    /// core API returns a feature's DBEntries in.
    #[test]
    fn ccds_and_refseq_from_transcript_xrefs() {
        let (_dir, dumps) = load_fixture(false);
        let tr = dumps.transcript("ENST00000622113").unwrap();
        assert_eq!(tr.ccds.as_deref(), Some("CCDS74772.1"));
        assert_eq!(tr.refseq, vec!["NM_001290224.2", "NM_199261.4"]);
        assert!(tr.is_canonical);
        assert!(!tr.readthrough);
        let other = dumps.transcript("ENST00000000210").unwrap();
        assert_eq!(other.ccds, None);
        assert!(other.refseq.is_empty());
        assert!(!other.is_canonical);
    }

    /// `cds_start_NF` and `cds_end_NF` are transcript attributes kept with
    /// their values (`Pipeline/DumpVEP/Dumper/Core.pm` `clean_transcript`), a
    /// `0` included: Ensembl VEP's FLAGS column tests the code's presence.
    #[test]
    fn cds_nf_from_transcript_attrib() {
        let (_dir, dumps) = load_fixture(false);
        let attrs = &dumps.transcript("ENST00000000210").unwrap().attributes;
        assert_eq!(
            attrs,
            &vec![
                Attribute {
                    code: "cds_start_NF".into(),
                    value: "1".into()
                },
                Attribute {
                    code: "cds_end_NF".into(),
                    value: "0".into()
                },
            ]
        );
    }

    /// A transcript whose dump rows list `cds_end_NF` before `cds_start_NF`
    /// carries them start first, the builder's fixed order where the Ensembl
    /// cache's varies by transcript; the attributes around them keep their
    /// places.
    #[test]
    fn cds_start_nf_precedes_cds_end_nf() {
        let dir = tempfile::tempdir().unwrap();
        fixture(dir.path());
        let mut rows = std::fs::read_to_string(dir.path().join("transcript_attrib.txt")).unwrap();
        rows.push_str("201\t2\tGENCODE basic\n201\t7\t1\n201\t3\ttsl5\n201\t6\t1\n");
        std::fs::write(dir.path().join("transcript_attrib.txt"), rows).unwrap();
        let dumps = CoreDumps::load(dir.path()).unwrap();
        let codes: Vec<&str> = dumps
            .transcript("ENST00000623753")
            .unwrap()
            .attributes
            .iter()
            .map(|a| a.code.as_str())
            .collect();
        assert_eq!(
            codes,
            ["gencode_basic", "cds_start_NF", "TSL", "cds_end_NF"]
        );
    }

    /// The kept attributes of a coding transcript, in table order, without
    /// the `remark` the dumper drops: the GENCODE labels, the TSL value with
    /// its `assigned to previous version` suffix, APPRIS and the MANE Select
    /// accession.
    #[test]
    fn kept_attributes_in_table_order() {
        let (_dir, dumps) = load_fixture(false);
        let codes: Vec<(&str, &str)> = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .attributes
            .iter()
            .map(|a| (a.code.as_str(), a.value.as_str()))
            .collect();
        assert_eq!(
            codes,
            vec![
                ("gencode_basic", "GENCODE basic"),
                ("TSL", "tsl1 (assigned to previous version 3)"),
                ("appris", "principal2"),
                ("MANE_Select", "NM_199261.4"),
                ("gencode_primary", "GENCODE primary"),
            ]
        );
        assert_eq!(
            dumps.transcript("ENST00000387314").unwrap().attributes,
            vec![Attribute {
                code: "TSL".into(),
                value: "tslNA".into()
            }]
        );
    }

    /// The `miRNA` attribute (`<start>-<end>` of the mature arm in cDNA
    /// coordinates) is a transcript attribute the dumper keeps.
    #[test]
    fn mirna_attribute_from_transcript_attrib() {
        let (_dir, dumps) = load_fixture(false);
        let attrs = &dumps.transcript("ENST00000411223").unwrap().attributes;
        assert!(attrs.contains(&Attribute {
            code: "miRNA".into(),
            value: "43-63".into()
        }));
    }

    /// The `ncRNA` attribute (`<start>:<end> <structure>`) is kept, and the
    /// MySQL escape `\t` in its value is decoded to the tab the database
    /// holds.
    #[test]
    fn ncrna_structure_from_transcript_attrib() {
        let (_dir, dumps) = load_fixture(false);
        let attrs = &dumps.transcript("ENST00000411223").unwrap().attributes;
        assert_eq!(
            attrs[0],
            Attribute {
                code: "ncRNA".into(),
                value: "1:104\t.11(3.(9.4(.(3.4)3.).)5.".into()
            }
        );
    }

    /// Protein features carry `seq_start`, `seq_end`, `hit_name` and the
    /// analysis display label, in feature id order; an analysis without a
    /// description has no label (`Database/Transcript.pm`
    /// `prefetch_translation_data` keeps `_display_label` alone).
    #[test]
    fn protein_features_with_display_label() {
        let (_dir, dumps) = load_fixture(false);
        let tl = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .translation
            .as_ref()
            .unwrap();
        assert_eq!(
            tl.protein_features,
            vec![
                ProteinFeature {
                    start: 187,
                    end: 334,
                    hseqname: "PF00102".into(),
                    display_label: Some("Pfam".into())
                },
                ProteinFeature {
                    start: 200,
                    end: 350,
                    hseqname: "G3DSA:3.90.190.10".into(),
                    display_label: Some("Gene3D".into())
                },
                ProteinFeature {
                    start: 10,
                    end: 20,
                    hseqname: "unknown.1".into(),
                    display_label: None
                },
            ]
        );
    }

    /// Sequence edits are the translation attributes `initial_met`,
    /// `_selenocysteine`, `amino_acid_sub` and `_stop_codon_rt`, in that code
    /// order (`Translation::get_all_SeqEdits`); the molecular weight is not
    /// one.
    #[test]
    fn seq_edits_in_get_all_seqedits_order() {
        let (_dir, dumps) = load_fixture(false);
        let tl = dumps
            .transcript("ENST00000622113")
            .unwrap()
            .translation
            .as_ref()
            .unwrap();
        assert_eq!(
            tl.seq_edits,
            vec![
                SeqEdit {
                    code: "initial_met".into(),
                    value: "1 1 M".into()
                },
                SeqEdit {
                    code: "_selenocysteine".into(),
                    value: "523 523 U".into()
                },
            ]
        );
    }

    /// The codon table is the `codon_table` sequence region attribute: 2 on
    /// MT, none on 21 (`Database/Transcript.pm` `prefetch_transcript_data`,
    /// default 1).
    #[test]
    fn codon_table_from_seq_region_attrib() {
        let (_dir, dumps) = load_fixture(false);
        assert_eq!(dumps.codon_table("MT"), Some(2));
        assert_eq!(dumps.codon_table("21"), None);
    }

    /// The source versions `Database/Transcript.pm` `info` writes:
    /// `genebuild.version` overwrites `genebuild.initial_release_date` under
    /// `genebuild`, and a key the table lacks is absent.
    #[test]
    fn source_versions_from_core_meta() {
        let (_dir, dumps) = load_fixture(false);
        let versions = dumps.source_versions();
        assert_eq!(versions["assembly"], "GRCh38.p14");
        assert_eq!(versions["gencode"], "GENCODE 50");
        assert_eq!(versions["genebuild"], "GENCODE50");
        assert_eq!(
            versions.keys().map(String::as_str).collect::<Vec<_>>(),
            ["assembly", "gencode", "genebuild"]
        );
    }

    /// A dump whose third line runs past the cap stops the read at that line,
    /// naming the file and the line; the rows before it are read, and a cap
    /// read at a small limit shows a line of exactly the cap passes.
    #[test]
    fn dump_line_past_the_cap_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let mut text = b"1\ta\n2\tb\n3\t".to_vec();
        text.extend(std::iter::repeat_n(b'x', MAX_LINE_BYTES));
        text.extend_from_slice(b"\n4\td\n");
        let path = dir.path().join("meta.txt");
        std::fs::write(&path, text).unwrap();
        let mut rows = read_table(dir.path(), "meta", 2).unwrap();
        assert_eq!(rows.next().unwrap().unwrap().field(1), "a");
        assert_eq!(rows.next().unwrap().unwrap().field(1), "b");
        let Some(Err(error)) = rows.next() else {
            panic!("line 3 is refused");
        };
        assert_eq!(
            error.to_string(),
            format!(
                "{} line 3 is longer than the {MAX_LINE_BYTES}-byte limit for one line",
                path.display()
            )
        );
        assert!(rows.next().is_none(), "the refused line ends the read");
        let lines: Vec<Vec<u8>> = BoundedLines::new(&b"ab\ncd"[..], 2, "x".into())
            .map(Result::unwrap)
            .collect();
        assert_eq!(lines, [b"ab".to_vec(), b"cd".to_vec()]);
        let refused = BoundedLines::new(&b"abc\n"[..], 2, "x".into())
            .next()
            .unwrap()
            .unwrap_err();
        assert_eq!(
            refused.to_string(),
            "x line 1 is longer than the 2-byte limit for one line"
        );
    }

    #[test]
    fn mysql_escapes_are_decoded() {
        assert_eq!(
            split_escaped("a\\tb\tc\\\\d\t\\N\te\\nf"),
            vec!["a\tb", "c\\d", "\\N", "e\nf"]
        );
        let row = Row {
            table: "t",
            fields: split_escaped("1\t\\N\tx"),
        };
        assert_eq!(row.optional(1), None);
        assert_eq!(row.optional(2), Some("x"));
        assert_eq!(row.optional(5), None);
        assert_eq!(row.field(5), "");
    }

    /// A numeric column that is not a number stops the build naming the
    /// table and the column.
    #[test]
    fn malformed_number_names_the_table_and_column() {
        let dir = tempfile::tempdir().unwrap();
        fixture(dir.path());
        let mut rows = std::fs::read_to_string(dir.path().join("protein_feature.txt")).unwrap();
        rows.push_str("70004\t300\tabc\t20\t1\t11\tx\t8633\t\\N\t\\N\t\\N\t\t\t\t\\N\n");
        std::fs::write(dir.path().join("protein_feature.txt"), rows).unwrap();
        let err = format!("{:#}", CoreDumps::load(dir.path()).unwrap_err());
        assert!(
            err.contains("protein_feature column 2 holds \"abc\", not a number"),
            "{err}"
        );
    }

    #[test]
    fn short_row_names_the_table_and_line() {
        let dir = tempfile::tempdir().unwrap();
        fixture(dir.path());
        let mut rows = std::fs::read_to_string(dir.path().join("gene.txt")).unwrap();
        rows.push_str("106\tprotein_coding\t1\n");
        std::fs::write(dir.path().join("gene.txt"), rows).unwrap();
        let err = format!("{:#}", CoreDumps::load(dir.path()).unwrap_err());
        assert!(
            err.contains("gene line 7 has 3 columns, the schema has 16"),
            "{err}"
        );
    }

    #[test]
    fn missing_core_table_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let err = CoreDumps::load(dir.path()).unwrap_err().to_string();
        assert!(err.contains("external_db.txt.gz"), "{err}");
    }
}
