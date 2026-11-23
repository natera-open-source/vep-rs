// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse GFF3 files into gene -> transcript -> exon/CDS hierarchy.
//!
//! A file is read as Ensembl VEP reads a GFF3 or GTF given to it as an
//! annotation source (`AnnotationSource/File/BaseGXF.pm`, `GFF.pm`,
//! `GTF.pm`): only the feature types of `GFF3_FEATURE_TYPES` are taken; a
//! record is named by `ID`, `Name`, `id` or `name` and placed under each of
//! its `Parent` records; a record without a parent is a gene when its type
//! names one (`gene`, `pseudogene`, `ncRNA_gene`, a gene segment,
//! `processed_transcript`, `RNA`) and otherwise a transcript of its own,
//! and one whose parent list opens with `0` and names nothing set is
//! dropped; a transcript's biotype is its `biotype`, `transcript_type` or
//! `transcript_biotype` attribute, else the RefSeq rules of
//! `Records::biotype`; its gene identifier is its `gene_id`, its `Dbxref`
//! `GeneID`, or the gene line's; its symbol is its `gene_name` or `gene`, or
//! the gene line's `name` or `Name`. An attribute whose value is empty or `0`
//! is absent, as Ensembl VEP's truth tests read it (`is_set`), so each of
//! those readers falls to its next alternative; the tests read a GFF3 value as
//! written, and its escapes are decoded where the value becomes the cache's
//! text (`text`). The same tree is built from a GTF by
//! [`crate::gtf::parse_gtf_model`].
//!
//! An Ensembl GFF3 carries the shape the Ensembl core API writes
//! (`Feature.pm`, `Gene.pm`, `Transcript.pm`, `Exon.pm` and `CDS.pm`
//! `summary_as_hash`): a gene is a line whose `ID` is `gene:<id>` and a
//! transcript a line whose `ID` is `transcript:<id>`, the Sequence Ontology
//! term in the type column following the biotype (`gene`, `ncRNA_gene`,
//! `pseudogene`; `mRNA`, `lnc_RNA`, `pseudogenic_transcript`,
//! `unconfirmed_transcript` and the rest), and the GRCh37 release-87 file
//! uses `pseudogene` for both a gene and a transcript line. A build of the
//! Ensembl cache's content ([`Reading::EnsemblCache`]) reads such a line
//! whatever its type, since the cache carries the TEC and scaRNA transcripts
//! whose terms the list lacks; a build of what Ensembl VEP's `--gff` source
//! annotates ([`Reading::GffSource`]) reads it by its type like any other.
//! Either way the line keeps the feature `version`, the transcript `tag`
//! list, `ccdsid` and `transcript_support_level`, the exon `ensembl_phase`
//! and `ensembl_end_phase`, and the CDS `protein_id` and translation
//! `version`. A transcript of biotype `artifact` or tagged readthrough is a
//! transcript of the tree like any other, as it is of Ensembl VEP's `--gff`
//! source (`AnnotationSource/File/BaseGXF.pm` `_create_transcript`); whether
//! the cache keeps it is decided where the cache is built, by the build's
//! inputs ([`crate::builder::build_cache`]).

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use tracing::{debug, info, warn};
use vep_core::coordinate::is_cache_path_component;
use vep_core::regulatory_gff::{
    decode_attribute, decode_type, undecodable_attribute_error, uri_unescape,
};
use vep_core::sv_routing::perl_split;
use vep_core::variation::perl_numeric;

use crate::count::{counted, one_or};
pub(crate) use crate::line::without_byte_order_marks;
use crate::line::{
    check_indexed_coordinates, index_build_places, index_columns_may_differ, is_a_comment,
    is_blank_or_a_comment, read_line, read_line_of_any_length, record_columns, Line,
    MAX_LINE_BYTES,
};

/// The GFF3 feature types Ensembl VEP reads from a file given as an
/// annotation source (`GFF.pm` `%INCLUDE_FEATURE_TYPES`; `ncRNA_gene` joins
/// the list at release 116.3). A line of any other type is skipped and
/// counted (`BaseGXF.pm` `_get_records_by_coords`), unless the file is read
/// for the Ensembl cache's content and the line's `ID` carries the Ensembl
/// `gene:` or `transcript:` prefix ([`Reading`]).
const GFF3_FEATURE_TYPES: [&str; 36] = [
    "aberrant_processed_transcript",
    "CDS",
    "C_gene_segment",
    "D_gene_segment",
    "exon",
    "gene",
    "J_gene_segment",
    "lnc_RNA",
    "lincRNA",
    "lincRNA_gene",
    "miRNA",
    "miRNA_gene",
    "mRNA",
    "mt_gene",
    "ncRNA",
    "ncRNA_gene",
    "nc_primary_transcript",
    "NMD_transcript_variant",
    "primary_transcript",
    "processed_pseudogene",
    "processed_transcript",
    "protein_coding_gene",
    "pseudogene",
    "pseudogenic_transcript",
    "RNA",
    "rRNA",
    "rRNA_gene",
    "snoRNA",
    "snoRNA_gene",
    "snRNA",
    "snRNA_gene",
    "supercontig",
    "transcript",
    "tRNA",
    "VD_gene_segment",
    "V_gene_segment",
];

/// A parsed gene record with its child transcripts.
#[derive(Debug)]
#[allow(dead_code)]
pub struct GeneRecord {
    pub gene_id: String,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    /// The gene's `Name`, its display xref's label.
    pub name: Option<String>,
    /// The gene line's `version`.
    pub version: Option<u32>,
    pub transcripts: Vec<TranscriptRecord>,
}

/// A parsed transcript record with its child exons and CDS regions.
#[derive(Debug)]
pub struct TranscriptRecord {
    pub transcript_id: String,
    /// The gene stable identifier the record names (`BaseGXF.pm`
    /// `_add_identifiers`): the transcript line's `gene_id`, else the
    /// `GeneID` of its `Dbxref`, else the gene line's `gene_id`, else the
    /// gene line's identifier; empty for a transcript with no gene line and
    /// none of those.
    pub gene_id: String,
    /// The gene symbol the record names: the transcript line's `gene_name`
    /// or `gene`, else the gene line's `name` or `Name`.
    pub symbol: Option<String>,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    pub source: String,
    /// `version`, else `transcript_version`.
    pub version: Option<u32>,
    /// `Ensembl_canonical` in `tag`.
    pub is_canonical: bool,
    /// `MANE_Select` in `tag`: the membership only. The RefSeq accession is
    /// the `MANE_Select` transcript attribute, which the GFF3 does not carry.
    pub mane_select: bool,
    /// `MANE_Plus_Clinical` in `tag`, the membership only.
    pub mane_plus_clinical: bool,
    /// `gencode_basic` in `tag`, or `basic`, the GRCh37 release-87 spelling,
    /// on an Ensembl gene or transcript line.
    pub gencode_basic: bool,
    /// `gencode_primary` in `tag`.
    pub gencode_primary: bool,
    /// `transcript_support_level`, else `tsl`: the `TSL` attribute's value
    /// with its `tsl` prefix removed (`Transcript.pm` `tsl`), so `1`, `NA` or
    /// `5 (assigned to previous version 6)`.
    pub tsl: Option<String>,
    /// `ccdsid`, the first CCDS xref's display label on an Ensembl GFF3 line
    /// (`Transcript.pm` `summary_as_hash`), else `ccds_id`, the xref's primary
    /// identifier on a GTF line (ensembl-io `GTFSerializer.pm`), each as the
    /// line carries it: the GFF3 writes `CCDS42939.1`, the GTF `CCDS42939`.
    pub ccds: Option<String>,
    /// `readthrough_transcript` or `readthrough_tra` in `tag`: the
    /// `readthrough_tra` transcript attribute, as the GRCh38 file and the
    /// core API spell it.
    pub readthrough: bool,
    /// `cds_start_NF` in `tag`: the model states that the coding sequence
    /// is incomplete at its start, so its first codon need be no start
    /// codon.
    pub cds_start_nf: bool,
    pub exons: Vec<ExonRecord>,
    pub cds_regions: Vec<CdsRecord>,
}

/// A parsed exon record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ExonRecord {
    pub exon_id: Option<String>,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub rank: Option<u32>,
    /// `ensembl_phase`, the exon's own phase (`Exon.pm` `summary_as_hash`):
    /// -1 outside the coding sequence.
    pub ensembl_phase: Option<i8>,
    /// `ensembl_end_phase`.
    pub ensembl_end_phase: Option<i8>,
}

/// A parsed CDS record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CdsRecord {
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub phase: i8,
    /// `protein_id`, the translation's stable identifier (`CDS.pm`
    /// `summary_as_hash`).
    pub protein_id: Option<String>,
    /// `version`, the translation's version; the GRCh37 release-87 file
    /// carries none. 1 when absent, as Ensembl VEP's `--gff` and `--gtf`
    /// sources give it (`BaseGXF.pm` `_add_translation`), except on a line
    /// read for the cache's content ([`Reading::EnsemblCache`]), which on a
    /// GTF takes its `protein_version` and otherwise keeps none, so the build
    /// takes the core table dumps' version when they hold the translation,
    /// else the protein FASTA header's.
    pub protein_version: Option<u32>,
    /// The record's line in the model's file, counted from 1 with the header
    /// lines: the place a refusal of the transcript's coding sequence names.
    pub line: usize,
}

/// The transcript model a file holds, which decides the biotype rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Model {
    /// A GFF3: `GFF.pm` `_record_get_biotype`.
    Gff3,
    /// A GTF: `GTF.pm` `_record_get_biotype`.
    Gtf,
}

/// Which lines of a GFF3 make the tree, whether a parentless record is a
/// gene by its type alone or by an Ensembl `gene:` identifier too
/// (`is_gene_line`), and where a translation's version comes from when its
/// CDS lines carry none, on a GTF as on a GFF3.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Reading {
    /// The lines of the types in `GFF3_FEATURE_TYPES`, as Ensembl VEP's
    /// `--gff` source reads a file (`BaseGXF.pm` `_get_records_by_coords`,
    /// `GFF.pm` `%INCLUDE_FEATURE_TYPES`): a line of another type is skipped
    /// whatever its `ID`, so an Ensembl `unconfirmed_transcript` (TEC) or
    /// `scRNA` (scaRNA) line makes no transcript and its exons are orphans,
    /// and a parentless record is a gene by its type alone (`BaseGXF.pm`
    /// `_record_is_gene`). A translation whose CDS lines carry no `version`
    /// has the version 1 (`BaseGXF.pm` `_add_translation`), as the `--gff`
    /// and `--gtf` sources give it.
    #[default]
    GffSource,
    /// Those lines and every line whose `ID` carries the Ensembl `gene:` or
    /// `transcript:` prefix, whatever its type: the Ensembl cache carries the
    /// TEC and scaRNA transcripts the type list leaves out. A parentless
    /// `gene:` line is a gene whatever its type: Ensembl writes the prefix on
    /// gene lines alone. A translation keeps no version when its CDS lines
    /// carry none (a GTF line's `protein_version` counts,
    /// [`crate::gtf::parse_gtf_model`]); the build takes the core table dumps'
    /// version when they hold the translation, else the protein FASTA
    /// header's.
    EnsemblCache,
}

/// One line of the transcript model that is not an exon or a CDS: a gene, a
/// transcript, or a feature the tree places under one of them.
#[derive(Debug)]
pub(crate) struct Feature {
    /// The record's identifier: `ID`, `Name`, `id` or `name` on a GFF3
    /// line, `<type>_id` on a GTF line, else the MD5 of the line
    /// (`GFF.pm` and `GTF.pm` `_record_get_id`).
    pub id: String,
    /// The set identifiers of its parent records (`parent_ids`): the
    /// `Parent` list of a GFF3 line; `gene_id` on a GTF transcript line.
    pub parents: Vec<String>,
    /// Whether the parent list opens with a non-empty field
    /// (`names_a_parent`), which makes the record a child, never a gene or
    /// a transcript of its own, whether or not a parent is set; a child whose
    /// `parents` is empty is dropped. A record whose list opens with an empty
    /// field is a child only when it is no gene and one of its `parents` is
    /// a gene line of the file ([`Records::assemble`]).
    pub child: bool,
    pub feature_type: String,
    /// The record's version, read under `feature_version_keys` as the line
    /// is, so a version that is not a whole number, or one the cache cannot
    /// hold, is refused at its line whether or not the tree keeps the record.
    pub version: Option<u32>,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub source: String,
    pub attrs: HashMap<String, String>,
    /// The `tag` values: the comma-separated GFF3 list, or the repeated GTF
    /// attribute.
    pub tags: Vec<String>,
    /// An Ensembl GFF3 gene or transcript line, whose `ID` carries the
    /// `gene:` or `transcript:` prefix.
    pub ensembl: bool,
}

/// The records of one file, collected line by line and assembled into the
/// gene tree once the file is read, since a child may precede its parent.
#[derive(Debug, Default)]
pub(crate) struct Records {
    /// How the file is read.
    reading: Reading,
    features: Vec<Feature>,
    /// Exons by parent record identifier; an exon with several parents is
    /// listed under each.
    exons: HashMap<String, Vec<ExonRecord>>,
    /// CDS segments by parent record identifier.
    cds: HashMap<String, Vec<CdsRecord>>,
    /// The lines skipped, counted by feature type.
    skipped: BTreeMap<String, usize>,
    /// The identifiers of the skipped lines that carry one, so a child
    /// naming such a parent is told its parent was skipped by type.
    skipped_ids: HashSet<String>,
}

/// What becomes of a record naming a parent in no line read, stated by the
/// two warnings that list such parents: a record whose list opens with such
/// a parent is placed as its other parents place it (`Parent=x,g1` is
/// `g1`'s transcript; an exon's `Parent=x,t1` is `t1`'s exon; a record
/// under a transcript alone makes no transcript, as with no `x`), and
/// dropped with its exons and CDS when it has none (`Parent=x`). A list
/// opening with an empty field is no such record ([`Feature::child`]).
const ORPHAN_RULE: &str = "a record whose parent list opens with one of them is placed as its other parents place it, else skipped with its children";

/// The transcripts the tree drops, reported once after the file is read.
#[derive(Default)]
struct Dropped {
    /// Transcripts whose biotype no rule names.
    no_biotype: Vec<String>,
    /// `protein_coding` transcripts without an exon line.
    no_exons: Vec<String>,
}

impl Records {
    /// An empty collection read under `reading`.
    pub(crate) fn new(reading: Reading) -> Self {
        Records {
            reading,
            ..Records::default()
        }
    }

    /// The feature types of the lines skipped so far, in name order.
    pub(crate) fn skipped_types(&self) -> Vec<String> {
        self.skipped.keys().cloned().collect()
    }

    /// Counts a line of a feature type the model does not read; `id` is the
    /// identifier its attributes name (`named_id`), empty when they name
    /// none.
    pub(crate) fn skip(&mut self, feature_type: &str, id: &str) {
        *self.skipped.entry(feature_type.to_string()).or_default() += 1;
        if !id.is_empty() {
            self.skipped_ids.insert(id.to_string());
        }
    }

    pub(crate) fn add_feature(&mut self, feature: Feature) {
        self.features.push(feature);
    }

    /// Lists `exon` under each of its parents.
    pub(crate) fn add_exon(&mut self, parents: &[String], exon: ExonRecord) {
        for parent in parents {
            self.exons
                .entry(parent.clone())
                .or_default()
                .push(exon.clone());
        }
    }

    /// Lists `cds` under each of its parents.
    pub(crate) fn add_cds(&mut self, parents: &[String], cds: CdsRecord) {
        for parent in parents {
            self.cds
                .entry(parent.clone())
                .or_default()
                .push(cds.clone());
        }
    }

    /// Builds the gene tree (`BaseGXF.pm` `_get_parent_child_structure` and
    /// `_create_transcripts`): an exon or CDS is placed under each parent
    /// its list names and a transcript under the first, in file order, of
    /// the gene lines its list names, a later record of a repeated
    /// identifier standing for the earlier; a parent in no line, or in a
    /// line of a type not read, is named in a warning, and a record whose
    /// list opens with such a parent is placed as its other parents place
    /// it and otherwise dropped with its children (`ORPHAN_RULE`); one
    /// whose parent list opens with `0` and names nothing set is dropped
    /// without a warning; each record without a parent is a gene whose
    /// children are its transcripts, or, when its type names no gene, a
    /// transcript with no gene line. A record whose parent list opens with
    /// an empty field ([`Feature::child`] false, `parents` set) is a child
    /// when it is no gene and one of its parents is a gene line, and a
    /// record without a parent otherwise.
    pub(crate) fn assemble(self, model: Model, path: &Path) -> Result<HashMap<String, GeneRecord>> {
        let Records {
            reading,
            features,
            mut exons,
            mut cds,
            skipped,
            skipped_ids,
        } = self;
        if !skipped.is_empty() {
            let skipped_types: Vec<String> = skipped
                .iter()
                .map(|(feature_type, count)| format!("{feature_type} ({count})"))
                .collect();
            info!(
                "Skipped the feature types {} of {}: not read as a gene, transcript, exon or CDS",
                skipped_types.join(", "),
                path.display()
            );
        }

        // The last record of a repeated identifier is the one the tree holds.
        let last: HashMap<&str, usize> = features
            .iter()
            .enumerate()
            .map(|(i, f)| (f.id.as_str(), i))
            .collect();
        let mut children: HashMap<&str, Vec<usize>> = HashMap::new();
        let mut missing: BTreeSet<&str> = BTreeSet::new();
        for (i, feature) in features.iter().enumerate() {
            if last[feature.id.as_str()] != i {
                continue;
            }
            for parent in &feature.parents {
                if last.contains_key(parent.as_str()) {
                    children.entry(parent.as_str()).or_default().push(i);
                } else {
                    missing.insert(parent.as_str());
                }
            }
        }
        for parent in exons.keys().chain(cds.keys()) {
            if !last.contains_key(parent.as_str()) {
                missing.insert(parent.as_str());
            }
        }
        let (skipped_parents, absent): (Vec<&str>, Vec<&str>) = missing
            .into_iter()
            .partition(|parent| skipped_ids.contains(*parent));
        if !absent.is_empty() {
            warn!(
                "{} {} no line of {}; {ORPHAN_RULE}: {}",
                counted(absent.len(), "parent identifier"),
                one_or(absent.len(), "names", "name"),
                path.display(),
                listed(absent.iter().copied())
            );
        }
        if !skipped_parents.is_empty() {
            warn!(
                "{} {} a line of {} of a type not read; {ORPHAN_RULE}: {}",
                counted(skipped_parents.len(), "parent identifier"),
                one_or(skipped_parents.len(), "names", "name"),
                path.display(),
                listed(skipped_parents.iter().copied())
            );
        }

        let mut genes: HashMap<String, GeneRecord> = HashMap::new();
        let mut dropped = Dropped::default();
        // A record whose parent list opens with an empty field stays at the top
        // level for Ensembl VEP (`BaseGXF.pm` `_get_parent_child_structure`,
        // `$parent_ids->[0] ne ''`) even when a later field is linked, so
        // `_create_transcripts` meets a transcript under a gene twice and its
        // gene record is that gene or none by hash order; here it is the gene's
        // child. A record under a transcript, or a gene under a gene, is a
        // record of its own in both: `_create_transcripts` reads a top-level
        // gene's children off before `_create_transcript` meets the gene's copy
        // under its parent gene, so that copy makes no transcript.
        let is_own_record = |feature: &Feature| {
            if feature.child {
                return false;
            }
            let has_gene_parent = feature.parents.iter().any(|parent| {
                last.get(parent.as_str())
                    .is_some_and(|&index| is_gene(reading, &features[index]))
            });
            is_gene(reading, feature) || !has_gene_parent
        };
        for (i, feature) in features.iter().enumerate() {
            if feature.child && feature.parents.is_empty() {
                debug!("Record {} names no set parent", feature.id);
            }
            if !is_own_record(feature) || last[feature.id.as_str()] != i {
                continue;
            }
            if is_gene(reading, feature) {
                let mut gene = GeneRecord {
                    gene_id: strip_gxf_prefix(&feature.id).to_string(),
                    chr: feature.chr.clone(),
                    start: feature.start,
                    end: feature.end,
                    strand: feature.strand,
                    biotype: attribute_text(model, &feature.attrs, &["biotype", "gene_biotype"])
                        .unwrap_or_default(),
                    name: attribute_text(model, &feature.attrs, &["name", "Name"]),
                    version: feature.version,
                    transcripts: Vec::new(),
                };
                for &child_index in children.get(feature.id.as_str()).into_iter().flatten() {
                    let child = &features[child_index];
                    if is_own_record(child) {
                        continue;
                    }
                    if let Some(tr) = transcript_record(
                        model,
                        reading,
                        child,
                        Some(feature),
                        &mut exons,
                        &mut cds,
                        &mut dropped,
                    ) {
                        gene.transcripts.push(tr);
                    }
                }
                // A transcript line naming this gene without a `Parent` was
                // placed first under a record of its own; the gene line's
                // record takes it over.
                if let Some(placed) = genes.remove(&gene.gene_id) {
                    gene.transcripts.extend(placed.transcripts);
                }
                genes.insert(gene.gene_id.clone(), gene);
            } else if let Some(tr) = transcript_record(
                model,
                reading,
                feature,
                None,
                &mut exons,
                &mut cds,
                &mut dropped,
            ) {
                // A transcript with no gene line hangs on a gene record of its
                // own coordinates, shared by the transcripts naming one gene.
                let key = if tr.gene_id.is_empty() {
                    format!("transcript:{}", feature.id)
                } else {
                    tr.gene_id.clone()
                };
                genes
                    .entry(key)
                    .or_insert_with(|| GeneRecord {
                        gene_id: tr.gene_id.clone(),
                        chr: tr.chr.clone(),
                        start: tr.start,
                        end: tr.end,
                        strand: tr.strand,
                        biotype: String::new(),
                        name: None,
                        version: None,
                        transcripts: Vec::new(),
                    })
                    .transcripts
                    .push(tr);
            }
        }

        if !dropped.no_biotype.is_empty() {
            warn!(
                "Unable to determine the biotype of {}, skipped: {}",
                counted(dropped.no_biotype.len(), "transcript"),
                listed(dropped.no_biotype.iter().map(String::as_str))
            );
        }
        if !dropped.no_exons.is_empty() {
            warn!(
                "No exons found for {}, skipped: {}",
                counted(dropped.no_exons.len(), "protein_coding transcript"),
                listed(dropped.no_exons.iter().map(String::as_str))
            );
        }
        Ok(genes)
    }

    /// The biotype of a transcript record, or `None` when no rule names
    /// one. A GFF3 record's is its `biotype`, `transcript_type` or
    /// `transcript_biotype`, else by its type (`GFF.pm` `_record_get_biotype`):
    /// `mRNA` is `protein_coding`, `ncRNA` takes its `ncrna_class`, a
    /// `<X>_gene_segment` is `IG_<X>_gene`, a transcript of a gene described
    /// `microRNA ...` is `miRNA`, and otherwise `gbkey` is the biotype. A GTF
    /// record's is its `transcript_biotype`, `transcript_type` or
    /// `biotype`, else its source column (`GTF.pm` `_record_get_biotype`).
    fn biotype(model: Model, tr: &Feature, gene: Option<&Feature>) -> Option<String> {
        match model {
            Model::Gtf => attribute_text(
                model,
                &tr.attrs,
                &["transcript_biotype", "transcript_type", "biotype"],
            )
            .or_else(|| Some(tr.source.clone()).filter(|s| is_set(s))),
            Model::Gff3 => {
                if let Some(biotype) = attribute_text(
                    model,
                    &tr.attrs,
                    &["biotype", "transcript_type", "transcript_biotype"],
                ) {
                    return Some(biotype);
                }
                let feature_type = tr.feature_type.to_ascii_lowercase();
                if feature_type == "mrna" {
                    return Some("protein_coding".to_string());
                }
                if feature_type == "ncrna" {
                    return attribute_text(model, &tr.attrs, &["ncrna_class"]);
                }
                if let Some(segment) = feature_type
                    .strip_suffix("_gene_segment")
                    .filter(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_lowercase()))
                {
                    return Some(format!("IG_{}_gene", segment.to_ascii_uppercase()));
                }
                if gene.is_some_and(|g| {
                    g.attrs
                        .get("description")
                        .is_some_and(|d| d.starts_with("microRNA"))
                }) {
                    return Some("miRNA".to_string());
                }
                attribute_text(model, &tr.attrs, &["gbkey"])
            }
        }
    }
}

/// The prefix of an Ensembl GFF3 gene line's `ID`.
const GENE_PREFIX: &str = "gene:";

/// The prefix of an Ensembl GFF3 transcript line's `ID`.
const TRANSCRIPT_PREFIX: &str = "transcript:";

/// Whether `id` opens with `prefix` in any letter case, the reading
/// `BaseGXF.pm` gives the two prefixes (`s/^(gene|transcript)://i`, lines 465
/// and 657), so every test of a prefix agrees with the strip.
fn has_prefix(id: &str, prefix: &str) -> bool {
    id.get(..prefix.len())
        .is_some_and(|head| head.eq_ignore_ascii_case(prefix))
}

/// Whether `id` carries the Ensembl `gene:` or `transcript:` prefix, in any
/// letter case ([`has_prefix`]).
fn has_gxf_prefix(id: &str) -> bool {
    has_prefix(id, GENE_PREFIX) || has_prefix(id, TRANSCRIPT_PREFIX)
}

/// Whether a GFF3 line's `ID` carries the Ensembl prefix ([`has_gxf_prefix`]),
/// read off the decoded identifier, the one the record keeps ([`record_id`],
/// [`is_gene_line`]).
fn has_ensembl_id(attrs: &HashMap<String, String>) -> bool {
    attrs
        .get("ID")
        .is_some_and(|id| has_gxf_prefix(&text(Model::Gff3, id)))
}

/// Whether a record without a parent is a gene: its type contains `gene`, or
/// is `processed_transcript` or `RNA`, the three tests of `BaseGXF.pm`
/// `_record_is_gene` (lines 788-796), read in the letter case Ensembl VEP reads
/// them. Under [`Reading::GffSource`] the type decides alone, as it does for
/// Ensembl VEP's `--gff` source, so a parentless `ID=gene:x` line of type
/// `mRNA` is a transcript. Under [`Reading::EnsemblCache`] an Ensembl `gene:`
/// line is a gene whatever its type, the convention of Ensembl's dumper,
/// which writes the prefix on gene-typed lines alone.
fn is_gene(reading: Reading, feature: &Feature) -> bool {
    is_gene_line(reading, &feature.id, &feature.feature_type)
}

/// [`is_gene`] on a line's identifier and type, before its record is built.
fn is_gene_line(reading: Reading, id: &str, feature_type: &str) -> bool {
    (reading == Reading::EnsemblCache && has_prefix(id, GENE_PREFIX))
        || feature_type.contains("gene")
        || feature_type == "processed_transcript"
        || feature_type == "RNA"
}

/// The keys a line's version is read from ([`Feature::version`]): `version`
/// on a line the tree makes a gene of its own, which is one that names no
/// parent in its first field and is a gene by [`is_gene_line`]; `version`,
/// else `transcript_version`, on every other.
pub(crate) fn feature_version_keys(
    reading: Reading,
    child: bool,
    id: &str,
    feature_type: &str,
) -> &'static [&'static str] {
    if !child && is_gene_line(reading, id, feature_type) {
        &["version"]
    } else {
        &["version", "transcript_version"]
    }
}

/// The transcript a record makes with its exons and CDS segments
/// (`BaseGXF.pm` `_create_transcript`), or `None` when it has neither, when
/// no rule names its biotype, or when it is `protein_coding` without an exon.
fn transcript_record(
    model: Model,
    reading: Reading,
    record: &Feature,
    gene: Option<&Feature>,
    exons_by_parent: &mut HashMap<String, Vec<ExonRecord>>,
    cds_by_parent: &mut HashMap<String, Vec<CdsRecord>>,
    dropped: &mut Dropped,
) -> Option<TranscriptRecord> {
    let mut exons = exons_by_parent.remove(&record.id).unwrap_or_default();
    let mut cds_regions = cds_by_parent.remove(&record.id).unwrap_or_default();
    if exons.is_empty() && cds_regions.is_empty() {
        debug!("Record {} has no exon or CDS line", record.id);
        return None;
    }
    // The prefix comes off the `transcript_id` attribute as it does off the
    // record's identifier (`BaseGXF.pm` `_create_transcript`, the
    // substitution on the `||` result).
    let named = attribute_text(model, &record.attrs, &["transcript_id"])
        .unwrap_or_else(|| record.id.clone());
    let transcript_id = strip_gxf_prefix(&named).to_string();
    let Some(biotype) = Records::biotype(model, record, gene) else {
        dropped.no_biotype.push(transcript_id);
        return None;
    };
    if biotype == "protein_coding" && exons.is_empty() {
        dropped.no_exons.push(transcript_id);
        return None;
    }

    let gene_id = attribute_text(model, &record.attrs, &["gene_id"])
        .or_else(|| dbxref_gene_id(model, &record.attrs))
        .or_else(|| gene.and_then(|g| attribute_text(model, &g.attrs, &["gene_id"])))
        .or_else(|| gene.map(|g| strip_gxf_prefix(&g.id).to_string()))
        .unwrap_or_default();
    let symbol = attribute_text(model, &record.attrs, &["gene_name", "gene"])
        .or_else(|| gene.and_then(|g| attribute_text(model, &g.attrs, &["name", "Name"])));

    // A rank compared with a start is no order at all, so the two keys
    // never share one comparator.
    if exons.iter().all(|e| e.rank.is_some()) {
        exons.sort_by_key(|e| e.rank);
    } else {
        if exons.iter().any(|e| e.rank.is_some()) {
            warn!("Transcript {transcript_id}: exons with and without a rank are ordered by start");
        }
        exons.sort_by_key(|e| e.start);
    }
    cds_regions.sort_by_key(|c| c.start);
    // Ensembl VEP's `--gff` and `--gtf` sources give a translation the
    // `version` of its CDS line, else 1 (`BaseGXF.pm` `_add_translation`); a
    // build of the Ensembl cache's content keeps an absent version absent
    // and takes the core table dumps' version when they hold the translation,
    // else the protein FASTA header's (`builder::translation_version`).
    if reading == Reading::GffSource {
        for segment in &mut cds_regions {
            segment.protein_version = segment.protein_version.or(Some(1));
        }
    }

    let tags = &record.tags;
    let tagged = |tag: &str| tags.iter().any(|t| t == tag);
    Some(TranscriptRecord {
        transcript_id,
        gene_id,
        symbol,
        chr: record.chr.clone(),
        start: record.start,
        end: record.end,
        strand: record.strand,
        biotype,
        source: record.source.clone(),
        version: record.version,
        is_canonical: tagged("Ensembl_canonical"),
        mane_select: tagged("MANE_Select"),
        mane_plus_clinical: tagged("MANE_Plus_Clinical"),
        gencode_basic: tagged("gencode_basic") || (record.ensembl && tagged("basic")),
        gencode_primary: tagged("gencode_primary"),
        tsl: attribute_text(model, &record.attrs, &["transcript_support_level", "tsl"]),
        ccds: attribute_text(model, &record.attrs, &["ccdsid", "ccds_id"]),
        readthrough: tagged("readthrough_transcript") || tagged("readthrough_tra"),
        cds_start_nf: tagged("cds_start_NF"),
        exons,
        cds_regions,
    })
}

/// Whether an attribute value is set, as Ensembl VEP's truth tests read it
/// (`if $x`, `||`): not empty and not the string `0`, so `ID=0;Name=foo`
/// names the record `foo` (`GFF.pm` `_record_get_id`).
pub(crate) fn is_set(value: &str) -> bool {
    !value.is_empty() && value != "0"
}

/// The first set value among `keys` ([`is_set`]), as a `||` chain over the
/// attribute hash reads them, as the hash holds it: a GFF3 value as written.
pub(crate) fn attribute(attrs: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
    first_set(attrs, keys).map(str::to_string)
}

/// [`attribute`] as the cache carries it ([`text`]): the first set value
/// among `keys`, a GFF3 value decoded after the truth test on its written
/// form.
pub(crate) fn attribute_text(
    model: Model,
    attrs: &HashMap<String, String>,
    keys: &[&str],
) -> Option<String> {
    first_set_text(model, attrs, keys).map(Cow::into_owned)
}

/// [`attribute_text`] without its copy: the value as [`text`] gives it,
/// borrowed from the hash where that reading borrows it.
fn first_set_text<'a>(
    model: Model,
    attrs: &'a HashMap<String, String>,
    keys: &[&str],
) -> Option<Cow<'a, str>> {
    first_set(attrs, keys).map(|value| text(model, value))
}

/// The first set value among `keys`, borrowed from the hash.
fn first_set<'a>(attrs: &'a HashMap<String, String>, keys: &[&str]) -> Option<&'a str> {
    keys.iter()
        .find_map(|k| attrs.get(*k).map(String::as_str).filter(|v| is_set(v)))
}

/// The `tag` list of a GFF3 line as Ensembl VEP reads it (`BaseGXF.pm`
/// lines 683-685): a set value ([`is_set`], line 683's truth test on the
/// value as written, so `0` and the empty value list nothing where `%30`
/// lists the tag `0`) split on the `,` written, so an escaped comma stays
/// inside its tag, each member then decoded ([`text`]) and compared whole
/// (lines 690-691 and 700), so a member padded with a blank names no tag.
fn tag_list(attrs: &HashMap<String, String>) -> Vec<String> {
    first_set(attrs, &["tag"])
        .map(|list| {
            perl_split(list, &[','])
                .map(|tag| text(Model::Gff3, tag).into_owned())
                .collect()
        })
        .unwrap_or_default()
}

/// The `GeneID` of a `dbxref` or `Dbxref` list (`GeneID:7157,Genbank:...`),
/// when the first `GeneID` pair's value is set (`BaseGXF.pm`
/// `_add_identifiers`, the `last` after the first pair, then `if(!...)`), the
/// list split and tested as the hash holds it and the identifier then
/// decoded ([`text`]).
fn dbxref_gene_id(model: Model, attrs: &HashMap<String, String>) -> Option<String> {
    first_set(attrs, &["dbxref", "Dbxref"])?
        .split(',')
        .find_map(|pair| pair.strip_prefix("GeneID:"))
        .map(|id| id.split(':').next().unwrap_or(id))
        .filter(|id| is_set(id))
        .map(|id| text(model, id).into_owned())
}

/// Up to ten of `ids`, comma-separated, then the count left out.
fn listed<'a>(ids: impl Iterator<Item = &'a str>) -> String {
    let ids: Vec<&str> = ids.collect();
    let mut text = ids.iter().take(10).copied().collect::<Vec<_>>().join(", ");
    if ids.len() > 10 {
        text.push_str(&format!(" and {} more", ids.len() - 10));
    }
    text
}

/// Parse a GFF3 file (plain, gzipped or bgzipped) into a gene map, from the
/// lines `reading` admits.
///
/// Returns `HashMap<gene_id, GeneRecord>`. Fails when the file cannot be
/// opened or read, and on a line it refuses, each error opening with the file
/// and the line (`located`): a line longer than `MAX_LINE_BYTES`, a line
/// that is not UTF-8, and, on a line whose type `reading` admits, a `start` or
/// `end` column that is not a number, naming the column and the value, a
/// start of 0 or an end before its start, a seqid that cannot name a
/// chromosome's cache directory ([`is_cache_path_component`]), naming the
/// value, a source column or an attribute value whose escapes decode to bytes
/// that are not UTF-8, naming the column or the key and the value, and a gene,
/// transcript or CDS `version` that is not a whole number or is above the
/// largest the cache holds (`parse_version`), naming the text or the limit,
/// and a CDS phase Ensembl VEP's source does not take (`parse_phase`),
/// naming the text.
/// A line of a type it skips is read no further than its type, its `ID` and
/// whether the index build parses its `start` and its `end`, and a line of
/// fewer than nine columns no further than those two
/// (`check_indexed_coordinates`): a column no number leads and an end
/// below zero are refused there as on a kept line, since no htslib's index
/// build parses either (`tbx.c` `tbx_parse1`). That build ends a column at
/// a tab or at a NUL byte and a reader at a tab (`crate::line`): a line
/// whose columns NUL bytes end is one the build places and one of fewer
/// than nine columns to the reader, skipped as such, and a line of a kept
/// type in which a NUL byte stands is held to the build's reading of its
/// start and its end (`index_columns_may_differ`). Nothing else of the
/// line refuses: the index build of Ensembl VEP's release (htslib 1.9)
/// parses a zero start, an end before its start and any other column a
/// number leads, takes any text as a seqid, and Ensembl VEP skips the line
/// by its type. A byte-order mark at the start of a line is taken off
/// before the line is read (`without_byte_order_marks`). A blank line is
/// skipped, and nothing is read from the first line that begins with `>`
/// and has a sequence behind it to the file's end: the next line that is
/// neither blank nor a comment exists and is one the index build could not
/// place, of fewer than four columns or with a start or an end that build
/// fails (`sequence_follows`). That is the sequence section the GFF3
/// specification allows a file, and the log names the line it opens at. A
/// `>` line with a line the index build places behind it, or with nothing,
/// is read as any other line: one that reads as a record is a record whose
/// seqid begins with `>`, which the index build files as a sequence of its
/// own while Ensembl VEP reads the records after it. So no section opens
/// in a tab-delimited file whose every data line the index build places.
/// A `##FASTA` directive is a comment and opens nothing: the index build
/// skips it as any comment, and the records behind one are read. A line
/// inside the section that reads as a record, and is neither a comment nor
/// led by `>`, is refused, naming its line and the line the section opens
/// at (`refuse_annotation_in_section`).
/// Any other line that is no comment and
/// writes fewer than four columns, as the index build ends them, is a
/// record without its start, which that build fails: it is refused, naming
/// the column it ends at. A line
/// that opens with a tab is read at the columns the file writes, its seqid
/// empty, as the index build reads it: a stray tab before a record puts its
/// type where a start is read, and the line refuses, the record's type kept
/// or skipped, unless a number leads the type's name. White space and a tab
/// before a record refuse it for that type alike
/// (`crate::line::opens_with_a_blank_column`); a line led by spaces alone,
/// no tab before its seqid, is read trimmed. A record's sequence is named as
/// `names` names its seqid ([`SequenceNames`]).
pub fn parse_gff3(
    path: &Path,
    reading: Reading,
    names: &SequenceNames,
) -> Result<HashMap<String, GeneRecord>> {
    Ok(parse_gff3_listing_skipped(path, reading, names)?.0)
}

/// [`parse_gff3`], with the feature types of the lines it skipped, in name
/// order: the types Ensembl VEP's `--gff` source warns of as it passes over
/// each such line (`BaseGXF.pm:239-243`).
pub(crate) fn parse_gff3_listing_skipped(
    path: &Path,
    reading: Reading,
    names: &SequenceNames,
) -> Result<(HashMap<String, GeneRecord>, Vec<String>)> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut records = Records::new(reading);
    // The line read ahead of the loop, behind a line that begins with `>`,
    // with its number, and the count of lines read off the file.
    let mut ahead = None;
    let mut read = 0;
    let mut section_opens_at = None;

    loop {
        let (line_number, line) = match ahead.take() {
            Some(numbered) => numbered,
            None => match read_line(&mut reader, path)? {
                Some(line) => {
                    read += 1;
                    (read, line)
                }
                None => break,
            },
        };
        if !line.whole {
            bail!(
                "{} line {line_number} is longer than the {MAX_LINE_BYTES}-byte limit for one line",
                path.display()
            );
        }
        let written = String::from_utf8(line.written)
            .map_err(|_| anyhow!("GFF3 line is not UTF-8"))
            .with_context(|| located(path, line_number))?;
        // A line that begins with `>` opens the sequence section when a
        // sequence follows it; with a line the index build places behind
        // it, or with nothing, it is a line like any other.
        if written.starts_with('>') && sequence_follows(&mut reader, path, &mut read, &mut ahead)? {
            info!(
                "{} holds a sequence section from line {line_number}: nothing from there on is read",
                path.display()
            );
            section_opens_at = Some(line_number);
            break;
        }
        if is_blank_or_a_comment(written.as_bytes()) {
            continue;
        }

        // A record's start is the fourth column as the index build ends
        // them, which fails a line that ends before it. A line that opens
        // with a tab is read at the columns the file writes, where the
        // index build reads them: its seqid is empty.
        let (line, fields) = record_columns(&written)
            .map_err(|columns| anyhow!("GFF3 line ends at column {columns}, before its start"))
            .with_context(|| located(path, line_number))?;
        let coordinate = |column: &str, value: &str| {
            parse_coordinate(column, value).with_context(|| located(path, line_number))
        };
        // A line the index build places in fewer than nine columns as a
        // reader ends them, one whose columns NUL bytes end among them, is
        // no record.
        if fields.len() < 9 {
            check_indexed_coordinates(&written, coordinate)?;
            continue;
        }

        // The include test runs on the decoded type, as Ensembl VEP's does, so
        // a type whose escapes spell bytes that are not UTF-8 matches no listed
        // type and its line is skipped. It runs before any other column is
        // read: Ensembl VEP's source reads a record's raw start, for its loop's
        // bound, and its type, and nothing else of one it skips
        // (`BaseGXF.pm:217-245`), so a skipped line's seqid and values refuse
        // nothing, nor does a coordinate the index build parses.
        let feature_type = decode_type(fields[2].as_bytes());
        let feature_type = feature_type.as_ref();
        let listed_type = GFF3_FEATURE_TYPES.contains(&feature_type);
        // A line of an unlisted type is kept by the Ensembl prefix of its
        // `ID` under the cache's reading alone, so its hash is built there,
        // to read the prefix. Under the `--gff` reading it is skipped
        // whatever its `ID`, and names itself off the column.
        let unlisted_hash =
            (!listed_type && reading == Reading::EnsemblCache).then(|| parse_attributes(fields[8]));
        let admitted = listed_type || unlisted_hash.as_ref().is_some_and(has_ensembl_id);
        if !admitted {
            check_indexed_coordinates(&written, coordinate)?;
            let id = match &unlisted_hash {
                Some(attrs) => named_id(attrs),
                None => skipped_id(fields[8]),
            };
            records.skip(feature_type, id.as_deref().unwrap_or(""));
            continue;
        }
        // White space and a tab before a record are a column to the index
        // build, which then reads the record's type as its start, and a NUL
        // byte ends a column to that build alone: the line is held to that
        // reading, as its twin of a skipped type is.
        if index_columns_may_differ(&written) {
            check_indexed_coordinates(&written, coordinate)?;
        }
        let chr = names.name(fields[0]);
        if !is_cache_path_component(chr) {
            return Err(anyhow!("invalid chromosome name '{}'", fields[0]))
                .with_context(|| located(path, line_number));
        }
        let start = coordinate("start", fields[3])?;
        let end = coordinate("end", fields[4])?;
        // Region binning subtracts one from a start and spans run from start to
        // end, so a start of 0 or an end before its start has no place to go.
        if start == 0 || end < start {
            return Err(anyhow!(
                "GFF3 coordinates {start}-{end} are not 1-based and ordered"
            ))
            .with_context(|| located(path, line_number));
        }
        let strand: i8 = match fields[6] {
            "-" => -1,
            _ => 1,
        };
        let attrs = unlisted_hash.unwrap_or_else(|| parse_attributes(fields[8]));
        let decoded = |key: &str| attrs.get(key).map(|value| text(Model::Gff3, value));
        let ensembl = has_ensembl_id(&attrs);
        check_attributes(fields[8], &attrs).with_context(|| located(path, line_number))?;
        let parents = parent_ids(Model::Gff3, &attrs, &["Parent", "parent"]);

        match feature_type {
            "exon" => records.add_exon(
                &parents,
                ExonRecord {
                    exon_id: decoded("Name").map(Cow::into_owned),
                    start,
                    end,
                    strand,
                    rank: decoded("rank").and_then(|r| r.parse::<u32>().ok()),
                    ensembl_phase: decoded("ensembl_phase").and_then(|p| p.parse::<i8>().ok()),
                    ensembl_end_phase: decoded("ensembl_end_phase")
                        .and_then(|p| p.parse::<i8>().ok()),
                },
            ),
            "CDS" => records.add_cds(
                &parents,
                CdsRecord {
                    start,
                    end,
                    strand,
                    phase: parse_phase(Model::Gff3, fields[7])
                        .with_context(|| located(path, line_number))?,
                    protein_id: attribute_text(Model::Gff3, &attrs, &["protein_id"]),
                    protein_version: parse_version(Model::Gff3, &attrs, &["version"])
                        .with_context(|| located(path, line_number))?,
                    line: line_number,
                },
            ),
            _ => {
                let tags = tag_list(&attrs);
                let id = record_id(&attrs, line);
                let child = names_a_parent(attribute(&attrs, &["Parent", "parent"]).as_deref());
                let version = parse_version(
                    Model::Gff3,
                    &attrs,
                    feature_version_keys(reading, child, &id, feature_type),
                )
                .with_context(|| located(path, line_number))?;
                records.add_feature(Feature {
                    id,
                    parents,
                    child,
                    feature_type: feature_type.to_string(),
                    version,
                    chr: chr.to_string(),
                    start,
                    end,
                    strand,
                    // Read on the lines the tree keeps, as `get_source` is read
                    // once a record passes the include test (`BaseGXF.pm:224,287`).
                    source: decoded_column("source", fields[1])
                        .with_context(|| located(path, line_number))?,
                    attrs,
                    tags,
                    ensembl,
                });
            }
        }
    }
    if let Some(opened) = section_opens_at {
        refuse_annotation_in_section(&mut reader, path, opened, read)?;
    }

    let skipped = records.skipped_types();
    Ok((records.assemble(Model::Gff3, path)?, skipped))
}

/// Whether a sequence follows a line that begins with `>`, behind which
/// `reader` stands, so that the line opens the sequence section: the next
/// line that is neither blank nor a comment exists and is one the index
/// build could not place ([`index_build_places`]). A line it places there,
/// or the file's end, leaves the `>` line a line like any other: the index
/// build places a `>` line whose start and end it parses and indexes the
/// lines behind it, which Ensembl VEP reads. A line longer than
/// `MAX_LINE_BYTES` is held to that many bytes, too few to tell what the
/// index build does with it: it is a sequence's when no tab stands in the
/// bytes held or in the rest behind the one byte [`read_line`] reads past
/// them ([`a_tab_stands_in_the_rest`]), and otherwise counts as placed, for
/// the loop to refuse for its length. `read` counts the lines read to tell.
/// The first of them that [`parse_gff3`] would not skip is left in `ahead`,
/// empty at the call, with its number, for the loop to read in its turn:
/// the line that decides or, before it, a blank line or a comment the loop
/// refuses, one longer than `MAX_LINE_BYTES` or not UTF-8.
fn sequence_follows(
    reader: &mut impl BufRead,
    path: &Path,
    read: &mut usize,
    ahead: &mut Option<(usize, Line)>,
) -> Result<bool> {
    debug_assert!(ahead.is_none());
    while let Some(line) = read_line(reader, path)? {
        *read += 1;
        let placed = if line.whole {
            index_build_places(&line.written)
        } else {
            let in_the_rest = a_tab_stands_in_the_rest(reader, path)?;
            in_the_rest || line.written.contains(&b'\t')
        };
        let passed_over = is_blank_or_a_comment(&line.written);
        let skipped = passed_over && line.whole && std::str::from_utf8(&line.written).is_ok();
        if !skipped && ahead.is_none() {
            *ahead = Some((*read, line));
        }
        if !passed_over {
            return Ok(!placed);
        }
    }
    Ok(false)
}

/// Passes over the rest of a line longer than `MAX_LINE_BYTES`, behind
/// which `reader` stands, and tells whether a tab stands in it.
fn a_tab_stands_in_the_rest(reader: &mut impl BufRead, path: &Path) -> Result<bool> {
    const PART_BYTES: u64 = 64 << 10;
    let mut a_tab_stands = false;
    let mut part = Vec::new();
    loop {
        part.clear();
        let read = (&mut *reader)
            .take(PART_BYTES)
            .read_until(b'\n', &mut part)
            .with_context(|| format!("Failed to read {}", path.display()))?;
        a_tab_stands = a_tab_stands || part.contains(&b'\t');
        if read == 0 || part.last() == Some(&b'\n') {
            return Ok(a_tab_stands);
        }
    }
}

/// Refuses annotation inside the sequence section that opens at line
/// `opened`, behind which `reader` stands with `read` lines of the file
/// read: the first line there that reads as a record
/// ([`reads_as_a_record`]) and is neither a comment ([`is_a_comment`]) nor
/// led by `>` is an error naming the file, the line and the line the
/// section opens at, since the section ends the annotation and the record
/// would go unread without a word. A `>` line there is a header whatever
/// its columns. A sequence may fill one line of any length
/// ([`read_line_of_any_length`]).
fn refuse_annotation_in_section(
    reader: &mut impl BufRead,
    path: &Path,
    opened: usize,
    read: usize,
) -> Result<()> {
    let mut line_number = read;
    while let Some(line) = read_line_of_any_length(reader, path)? {
        line_number += 1;
        let line = line.written;
        if reads_as_a_record(&line) && !line.starts_with(b">") && !is_a_comment(&line) {
            return Err(anyhow!(
                "GFF3 annotation after the sequence section, which opens at line {opened}"
            ))
            .with_context(|| located(path, line_number));
        }
    }
    Ok(())
}

/// Whether `line`, as the file writes it, reads as a record: nine or more
/// columns, the fourth and the fifth whole numbers. It is the test of a
/// record on a line inside the sequence section
/// ([`refuse_annotation_in_section`]).
fn reads_as_a_record(line: &[u8]) -> bool {
    let number = |column: Option<&[u8]>| {
        column
            .and_then(|text| std::str::from_utf8(text).ok())
            .is_some_and(|text| text.parse::<u64>().is_ok())
    };
    let mut columns = line.split(|byte| *byte == b'\t');
    number(columns.nth(3)) && number(columns.next()) && columns.nth(3).is_some()
}

/// How a build names the sequences its annotation files name.
///
/// A leading `chr` is the UCSC spelling of a name where a digit or a capital
/// letter follows it (`chr1`, `chrX`, `chrUn_KI270302v1`, `chrI`, `chr2L`) and
/// is dropped there, and wherever the genome FASTA holds the shortened name
/// and not the written one. Any other name that begins with the three letters
/// (`chromosome_1`, `chrom3`) is a name of its own and is kept as written.
///
/// A mitochondrion written `M` or `chrM` is named `MT` where the build reads
/// no genome FASTA, where the FASTA reader reaches a record under `MT` (a
/// record named `M` or `chrM`, or `MT` or `chrMT` with those two letters in
/// either case, as `Mt`), and where the FASTA names no organelle genome at
/// all (a part of a genome, such as one chromosome). Under a FASTA that
/// reaches no record under `MT` and names an organelle genome of its own
/// (`Mito`, `mitochondrion_genome`, `Pt`:
/// `crate::builder::is_organelle_name`) the name is kept as written, since
/// `MT` would name a sequence the genome does not hold.
///
/// On GRCh37 `MT` is the assembly's mitochondrial sequence, of
/// [`vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH`] bases, and `chrM` and `M` also
/// name the one of the UCSC hg19 naming, of another length, on which a
/// file's models of that name can count. There the record the FASTA reader
/// reaches under `MT` counts only where it has the assembly's length; over
/// any other the name is kept as written, as under a FASTA that names an
/// organelle genome of its own ([`Self::keeps_apart_from_mt`]).
///
/// The default reads no FASTA: the names of a build without one, and the
/// reading `--codon-table` gives its regions.
#[derive(Debug, Default)]
pub struct SequenceNames {
    fasta: Option<FastaNames>,
    /// Whether the naming is GRCh37's over a FASTA whose reader reaches a
    /// record of another length than the assembly's under `MT`, or reaches
    /// none and names an organelle genome.
    apart_on_grch37: bool,
}

/// What [`SequenceNames`] reads of a genome FASTA: its records' names, and
/// whether it names an organelle genome while its reader reaches no record
/// under `MT`, or on GRCh37 its reader reaches one of another length than
/// the assembly's.
#[derive(Debug)]
struct FastaNames {
    names: HashSet<String>,
    names_its_own_organelle: bool,
}

impl SequenceNames {
    /// The names of a build or a run that reads `fasta` as its genome, on
    /// `assembly` where it names one.
    pub fn of_fasta(fasta: &vep_fasta::IndexedFasta, assembly: Option<&str>) -> Self {
        let names: HashSet<String> = fasta.record_names().map(str::to_string).collect();
        let on_grch37 = assembly == Some("GRCh37");
        let reached = fasta.sequence_length("MT");
        // A record the run's synonym table gives for `MT` can be named like
        // no organelle genome, so its length alone keeps the names apart.
        let of_another_length = on_grch37
            && reached.is_some_and(|length| length != vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH);
        let names_its_own_organelle = of_another_length
            || (reached.is_none()
                && names
                    .iter()
                    .any(|name| crate::builder::is_organelle_name(name)));
        Self {
            fasta: Some(FastaNames {
                names,
                names_its_own_organelle,
            }),
            apart_on_grch37: on_grch37 && names_its_own_organelle,
        }
    }

    /// Whether `named`, a name `Self::name` gives, is a mitochondrion this
    /// naming keeps apart from GRCh37's `MT`: an `M` or a `chrM` kept as
    /// written because the FASTA gives no record of the assembly's length
    /// for `MT`.
    pub fn keeps_apart_from_mt(&self, named: &str) -> bool {
        self.apart_on_grch37 && matches!(named, "M" | "chrM")
    }

    /// The cache's name for the sequence an annotation file writes as
    /// `written`.
    pub(crate) fn name<'a>(&self, written: &'a str) -> &'a str {
        let shortened = match written.strip_prefix("chr") {
            None => written,
            Some(rest) if is_ucsc_designation(rest) || self.holds_the_shortened(written, rest) => {
                rest
            }
            Some(_) => return written,
        };
        if shortened != "M" {
            return shortened;
        }
        match &self.fasta {
            Some(fasta) if fasta.names_its_own_organelle => written,
            _ => "MT",
        }
    }

    /// Whether the genome FASTA holds a record named `rest` and none named
    /// `written`.
    fn holds_the_shortened(&self, written: &str, rest: &str) -> bool {
        self.fasta
            .as_ref()
            .is_some_and(|fasta| fasta.names.contains(rest) && !fasta.names.contains(written))
    }
}

/// Whether `rest`, what follows a leading `chr`, is the chromosome's
/// designation of a UCSC name: it opens with a digit or a capital letter.
fn is_ucsc_designation(rest: &str) -> bool {
    rest.as_bytes()
        .first()
        .is_some_and(|byte| byte.is_ascii_digit() || byte.is_ascii_uppercase())
}

/// A column's text with its escapes decoded ([`uri_unescape`], the
/// `decode_string` of `GFF3.pm:153`); escapes spelling bytes that are not
/// UTF-8 are an error naming the column and the value, as an attribute
/// value's are ([`check_attributes`]).
fn decoded_column(column: &str, value: &str) -> Result<String> {
    match uri_unescape(value) {
        Cow::Borrowed(_) => Ok(value.to_string()),
        Cow::Owned(bytes) => String::from_utf8(bytes)
            .map_err(|_| anyhow!("GFF3 {column} {value:?} decodes to text that is not UTF-8")),
    }
}

/// The `start` or `end` column as a number; a value that is not one is an
/// error naming the column and the value.
fn parse_coordinate(column: &str, value: &str) -> Result<u64> {
    value
        .parse()
        .map_err(|_| anyhow!("GFF3 {column} {value:?} is not a number"))
}

/// The longest leading part of a phase column read as a number: text past it
/// has no part in the comparison with 1 and 2.
const PHASE_HEAD_CHARS: usize = 64;

/// The phase column of a CDS line as a number, -1 when it states none. The
/// column is text Ensembl VEP's source takes for an exon's phase: text that
/// holds a `.` states none (ensembl-io `GXF.pm:264`), text that reads as the
/// number 1 or 2 is exchanged for the other by a comparison of numbers
/// (`BaseGXF.pm:815-820`), so it is that number here whatever follows the
/// digit, and any other text is `-1` or `0`, the values the exon takes
/// beside them (core `Exon.pm:244-250`, whose `throw` ends Ensembl VEP's run
/// with `Bad value (...) for exon phase`). Text the source does not take is
/// an error naming the value.
pub(crate) fn parse_phase(model: Model, value: &str) -> Result<i8> {
    let head = value
        .char_indices()
        .nth(PHASE_HEAD_CHARS)
        .map_or(value, |(at, _)| &value[..at]);
    let number = perl_numeric(head);
    if value.contains('.') {
        return Ok(-1);
    }
    if number == 1.0 || number == 2.0 || matches!(value, "0" | "-1") {
        return Ok(number as i8);
    }
    let format = match model {
        Model::Gff3 => "GFF3",
        Model::Gtf => "GTF",
    };
    bail!("{format} phase {value:?} is not a phase")
}

/// Holds the `start` and `end` of a line a reader skips, for its type or for
/// having fewer than nine columns, to the one rule they keep there: `parse`,
/// the reader's own reading of a column, runs on the first of the two the
/// index build fails ([`vep_core::regulatory_gff::unindexed_coordinate`]),
/// so its error refuses the line as it would a kept one. That build fails a
/// column no number leads and an end that reads below zero, in every htslib,
/// and parses any other text a number leads (`100x`, `-5`, `0x1A`, an empty
/// start before a numeric end), which stays unread here. `written` is the
/// line as the file writes it, whose fourth and fifth columns the index
/// build reads: a blank one that ends the line is a column, and a line that
/// opens with a tab keeps each column in its place. A kept line that opens
/// with white space and a tab is held to the rule too
/// ([`crate::line::opens_with_a_blank_column`]): the index build reads its
/// type where a reader reads its start.
pub(crate) fn check_skipped_coordinates(
    written: &str,
    parse: impl Fn(&str, &str) -> Result<u64>,
) -> Result<()> {
    // A carriage return that ends the line is no part of its last column.
    let written = written.strip_suffix('\r').unwrap_or(written);
    if let Some((column, value)) =
        vep_core::regulatory_gff::unindexed_coordinate(written.as_bytes())
    {
        parse(column, &String::from_utf8_lossy(value))?;
    }
    Ok(())
}

/// An identifier without its `gene:` or `transcript:` prefix, in any letter
/// case (`BaseGXF.pm` `s/^(gene|transcript)://i`).
fn strip_gxf_prefix(id: &str) -> &str {
    for prefix in [GENE_PREFIX, TRANSCRIPT_PREFIX] {
        if has_prefix(id, prefix) {
            return &id[prefix.len()..];
        }
    }
    id
}

/// The first set value among `keys` ([`first_set_text`]) as a version number
/// ([`version_number`]); no set value is no version.
pub(crate) fn parse_version(
    model: Model,
    attrs: &HashMap<String, String>,
    keys: &[&str],
) -> Result<Option<u32>> {
    first_set_text(model, attrs, keys)
        .map(|version| version_number(&version))
        .transpose()
}

/// A version as a file writes it: ASCII digits alone, read as the number
/// (`01` is 1). Any other text is an error naming it, and digits that do not
/// fit `u32`, the cache's version field, an error naming them and the limit:
/// a version read as none would lose what the file states, where Ensembl VEP
/// prints the text as written.
pub(crate) fn version_number(text: &str) -> Result<u32> {
    if !is_whole_number(text) {
        bail!("version {text:?} is not a whole number");
    }
    text.parse::<u32>().map_err(|_| {
        anyhow!(
            "version {text} is above {}, the largest version the cache holds",
            u32::MAX
        )
    })
}

/// Whether `text` is a whole number as a file writes one: ASCII digits
/// alone, at least one.
pub(crate) fn is_whole_number(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit())
}

/// The place of a refused line, `<file> line <n>`, for the context of its
/// error; `line` is counted from 1 with the header lines.
pub(crate) fn located(path: &Path, line: usize) -> String {
    format!("{} line {line}", path.display())
}

/// The keys a GFF3 line's identifier is read from, in the order they are
/// tried (`GFF.pm` `_record_get_id`).
const ID_KEYS: [&str; 4] = ["ID", "Name", "id", "name"];

/// The identifier a GFF3 line's attributes name (`GFF.pm` `_record_get_id`):
/// the first of `ID`, `Name`, `id` and `name` set, when that value as written
/// holds a word character (`/\w+/`), decoded ([`text`]).
fn named_id(attrs: &HashMap<String, String>) -> Option<String> {
    identifier(first_set(attrs, &ID_KEYS))
}

/// [`named_id`] read off the attribute column of a line whose hash is not
/// built: the value the hash would hold under each of [`ID_KEYS`], which is
/// its last pair's, none when that pair has no `=` ([`parse_attributes`]),
/// then the first of them set.
fn skipped_id(column: &str) -> Option<String> {
    let mut held: [Option<&str>; ID_KEYS.len()] = [None; ID_KEYS.len()];
    for (key, value) in attribute_pairs(column) {
        if let Some(slot) = ID_KEYS.iter().position(|id_key| *id_key == key) {
            held[slot] = value;
        }
    }
    identifier(held.into_iter().flatten().find(|value| is_set(value)))
}

/// The identifier the first set value of [`ID_KEYS`] names: `value` decoded
/// ([`text`]), when as written it holds a word character.
fn identifier(value: Option<&str>) -> Option<String> {
    value
        .filter(|id| id.chars().any(|c| c.is_alphanumeric() || c == '_'))
        .map(|id| text(Model::Gff3, id).into_owned())
}

/// A GFF3 record's identifier: [`named_id`], or the MD5 of the line when the
/// attributes name none.
fn record_id(attrs: &HashMap<String, String>, line: &str) -> String {
    named_id(attrs).unwrap_or_else(|| format!("{:x}", md5::compute(line)))
}

/// The records an exon, CDS or transcript line belongs to: the comma-separated
/// list under the first of `keys` set (`Parent` or `parent` on a GFF3 line,
/// `GFF.pm` `_record_get_parent_ids`; `transcript_id` or `gene_id` on a GTF
/// line, `GTF.pm` `%PARENTS`), split and tested as the hash holds it, without
/// the identifiers that are not set (`BaseGXF.pm`
/// `_get_parent_child_structure`, `next unless ($parent_id)`), each then
/// decoded ([`text`]).
pub(crate) fn parent_ids(
    model: Model,
    attrs: &HashMap<String, String>,
    keys: &[&str],
) -> Vec<String> {
    first_set(attrs, keys)
        .map(|list| {
            list.split(',')
                .filter(|parent| is_set(parent))
                .map(|parent| text(model, parent).into_owned())
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

/// Whether a parent list makes its record a child by Ensembl VEP's prune of
/// the top level (`BaseGXF.pm` `_get_parent_child_structure`,
/// `scalar(@$parent_ids) && $parent_ids->[0] ne ''`): the list's first
/// field is non-empty, set or not. `list` is the value Ensembl VEP splits on
/// `,`: the `Parent || parent` result of a GFF3 line (`GFF.pm`
/// `_record_get_parent_ids`), the `gene_id` of a GTF transcript line as
/// written (`GTF.pm` `_record_get_parent_ids`). An unset first field,
/// `0`, is still non-empty, while `next unless ($parent_id)` links the
/// record to nothing, so a GTF transcript line of `gene_id "0"` is dropped.
pub(crate) fn names_a_parent(list: Option<&str>) -> bool {
    let first_field = list.and_then(|list| list.split(',').next());
    first_field.is_some_and(|field| !field.is_empty())
}

/// A GFF3 attribute column as ensembl-io's `Parser/GFF3.pm` `get_attributes`
/// splits it (lines 190-198): the column split on `;`, each pair on `=` into
/// two variables (line 194, a `split` Perl runs with a limit of three), the
/// key the text before the first `=` and the value the text between the first
/// and the second, so `Name=a=b` is `a`; a later pair with a key replaces an
/// earlier one (line 195), a pair without `=` with the undefined value, so it
/// sets nothing and unsets an earlier pair with its key. Nothing is trimmed:
/// a key opening with a space is that key. The values are kept as written,
/// the hash every test of Ensembl VEP's `--gff` source reads (its truth
/// tests, its word-character test, its splits on `,`), and are decoded where
/// a value becomes the cache's text ([`text`]), once [`check_attributes`] has
/// admitted them.
pub(crate) fn parse_attributes(attrs_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for (key, value) in attribute_pairs(attrs_str) {
        match value {
            Some(value) => {
                map.insert(key.to_string(), value.to_string());
            }
            None => {
                map.remove(key);
            }
        }
    }
    map
}

/// The pairs of a GFF3 attribute column in written order, each its key and
/// its value as [`parse_attributes`] cuts them, `None` for a pair without `=`.
fn attribute_pairs(column: &str) -> impl Iterator<Item = (&str, Option<&str>)> {
    // Perl's `split` in Ensembl VEP drops the trailing empty fields, so a closing
    // `;` is no pair.
    column.trim_end_matches(';').split(';').map(|pair| {
        let mut fields = pair.splitn(3, '=');
        (fields.next().unwrap_or_default(), fields.next())
    })
}

/// Refuses a value of the hash built from `column` whose escapes spell bytes
/// that are not UTF-8, with the text both readers share
/// ([`undecodable_attribute_error`]); of several such values, the one whose
/// pair comes first in the column as written, among the pairs whose value the
/// hash holds (the last pair of each key), as the GTF reader names one. The
/// hash is tested first, so a line with nothing to refuse, every line of the
/// files Ensembl distributes, is read once; the column is walked only to name
/// the offender. Run on the lines a reader keeps, the lines whose hash Ensembl
/// VEP builds (`BaseGXF.pm:224,289`), so a skipped line's values refuse
/// nothing.
pub(crate) fn check_attributes(column: &str, attrs: &HashMap<String, String>) -> Result<()> {
    if attrs.values().all(|value| decodes_to_text(value)) {
        return Ok(());
    }
    let kept = |key: &str, value: &str| attrs.get(key).is_some_and(|held| held == value);
    let refused = attribute_pairs(column)
        .filter_map(|(key, value)| value.map(|value| (key, value)))
        .find(|(key, value)| kept(key, value) && !decodes_to_text(value));
    match refused {
        Some((key, value)) => bail!(undecodable_attribute_error("GFF3", key, value)),
        None => Ok(()),
    }
}

/// Whether a value's escapes spell UTF-8 once decoded ([`uri_unescape`]); a
/// value without an escape does.
fn decodes_to_text(value: &str) -> bool {
    match uri_unescape(value) {
        Cow::Owned(bytes) => std::str::from_utf8(&bytes).is_ok(),
        Cow::Borrowed(_) => true,
    }
}

/// A value of the attribute hash as the cache carries it. A GFF3 value is
/// decoded ([`decode_attribute`]) here, once every test Ensembl VEP's `--gff`
/// source runs on the hash has run on the value as written (`GFF3.pm:190-198`
/// stores it so): the truth tests ([`is_set`]), the word-character test of
/// [`named_id`] and the splits on `,` of [`parent_ids`], [`tag_list`] and
/// [`dbxref_gene_id`], so `%30` is a set value reading `0` and `a%2Cb` one
/// name. A GTF value is as the hash holds it, decoded when the hash was built
/// (`GTF.pm:126`), so the tests there ran on the decoded value, as Ensembl
/// VEP's do.
pub(crate) fn text(model: Model, value: &str) -> Cow<'_, str> {
    match model {
        Model::Gff3 => decode_attribute(value.as_bytes()),
        Model::Gtf => Cow::Borrowed(value),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parses `lines` (GFF3 body lines, no header) from a temporary file
    /// under `reading`.
    fn try_parse_as(lines: &str, reading: Reading) -> Result<HashMap<String, GeneRecord>> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gff3");
        std::fs::write(&path, lines).unwrap();
        parse_gff3(&path, reading, &SequenceNames::default())
    }

    /// Parses `lines` as Ensembl VEP's `--gff` source reads them.
    fn try_parse(lines: &str) -> Result<HashMap<String, GeneRecord>> {
        try_parse_as(lines, Reading::GffSource)
    }

    fn parse(lines: &str) -> HashMap<String, GeneRecord> {
        try_parse(lines).unwrap()
    }

    fn transcript<'a>(genes: &'a HashMap<String, GeneRecord>, id: &str) -> &'a TranscriptRecord {
        genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .find(|t| t.transcript_id == id)
            .unwrap_or_else(|| panic!("transcript {id} not parsed"))
    }

    // Lines of Homo_sapiens.GRCh38.116.gff3.gz, chromosome 21 (the transcript `Name`
    // attribute, which is not read, left out).
    const TPTE_GENE: &str = "21\tensembl_havana\tgene\t10521553\t10606140\t.\t+\t.\tID=gene:ENSG00000274391;Name=TPTE;biotype=protein_coding;description=transmembrane phosphatase with tensin homology [Source:HGNC Symbol%3BAcc:HGNC:12023];gene_id=ENSG00000274391;logic_name=ensembl_havana_gene_homo_sapiens;version=6\n";
    const TPTE_205: &str = "21\tensembl_havana\tmRNA\t10521553\t10606140\t.\t+\t.\tID=transcript:ENST00000622113;Parent=gene:ENSG00000274391;biotype=protein_coding;ccdsid=CCDS74772.1;tag=gencode_basic,gencode_primary;transcript_id=ENST00000622113;transcript_support_level=1;version=4\n";
    const TPTE_205_EXON: &str = "21\tensembl_havana\texon\t10521553\t10521694\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00003729301;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003729301;rank=1;version=1\n";
    const TPTE_205_CDS: &str = "21\tensembl_havana\tCDS\t10538724\t10538734\t.\t+\t0\tID=CDS:ENSP00000482040;Parent=transcript:ENST00000622113;protein_id=ENSP00000482040;version=1\n\
21\tensembl_havana\tCDS\t10541112\t10541165\t.\t+\t1\tID=CDS:ENSP00000482040;Parent=transcript:ENST00000622113;protein_id=ENSP00000482040;version=1\n";

    // Lines of Homo_sapiens.GRCh37.87.gff3.gz, chromosome 21: a `pseudogene`
    // gene line and a `pseudogene` transcript line, the transcript tagged `basic`.
    const GRCH37_PSEUDOGENE: &str = "21\thavana\tpseudogene\t10872119\t10872429\t.\t+\t.\tID=gene:ENSG00000235476;Name=AF254982.2;biotype=pseudogene;gene_id=ENSG00000235476;logic_name=havana;version=1\n\
21\thavana\tpseudogene\t10872119\t10872429\t.\t+\t.\tID=transcript:ENST00000411878;Parent=gene:ENSG00000235476;Name=AF254982.2-001;biotype=unprocessed_pseudogene;havana_transcript=OTTHUMT00000157411;havana_version=1;tag=basic;transcript_id=ENST00000411878;version=1\n\
21\thavana\texon\t10872119\t10872429\t.\t+\t.\tParent=transcript:ENST00000411878;Name=ENSE00001618434;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001618434;rank=1;version=1\n";

    /// Lines in the shape of the RefSeq GRCh38 annotation (`GCF_000001405`):
    /// TP53's gene, mRNA, exon and CDS lines, an `ncRNA` with its
    /// `ncrna_class`, a `V_gene_segment`, a `miRNA` whose gene is described
    /// `microRNA 548X`, a `transcript` carrying only `gbkey`, and a record
    /// of no biotype at all.
    const REFSEQ: &str = "NC_000017.11\tBestRefSeq\tgene\t7668402\t7687550\t.\t-\t.\tID=gene-TP53;Dbxref=GeneID:7157,HGNC:HGNC:11998;Name=TP53;description=tumor protein p53;gbkey=Gene;gene=TP53;gene_biotype=protein_coding\n\
NC_000017.11\tBestRefSeq\tmRNA\t7668402\t7687550\t.\t-\t.\tID=rna-NM_000546.6;Parent=gene-TP53;Dbxref=GeneID:7157,Genbank:NM_000546.6,HGNC:HGNC:11998;Name=NM_000546.6;gbkey=mRNA;gene=TP53;product=tumor protein p53%2C transcript variant 1;tag=MANE Select;transcript_id=NM_000546.6\n\
NC_000017.11\tBestRefSeq\texon\t7687377\t7687550\t.\t-\t.\tID=exon-NM_000546.6-1;Parent=rna-NM_000546.6;Dbxref=GeneID:7157,Genbank:NM_000546.6;gbkey=mRNA;gene=TP53;transcript_id=NM_000546.6\n\
NC_000017.11\tBestRefSeq\texon\t7676521\t7676622\t.\t-\t.\tID=exon-NM_000546.6-2;Parent=rna-NM_000546.6;Dbxref=GeneID:7157,Genbank:NM_000546.6;gbkey=mRNA;gene=TP53;transcript_id=NM_000546.6\n\
NC_000017.11\tBestRefSeq\tCDS\t7676521\t7676594\t.\t-\t0\tID=cds-NP_000537.3;Parent=rna-NM_000546.6;Dbxref=GeneID:7157,Genbank:NP_000537.3;Name=NP_000537.3;gbkey=CDS;gene=TP53;product=cellular tumor antigen p53 isoform a;protein_id=NP_000537.3\n\
NC_000017.11\tBestRefSeq\tgene\t7700000\t7701000\t.\t+\t.\tID=gene-LINC1;Dbxref=GeneID:1;Name=LINC1;gbkey=Gene;gene=LINC1;gene_biotype=lncRNA\n\
NC_000017.11\tBestRefSeq\tncRNA\t7700000\t7701000\t.\t+\t.\tID=rna-NR_1;Parent=gene-LINC1;Dbxref=GeneID:1,Genbank:NR_1;Name=NR_1;gbkey=ncRNA;gene=LINC1;ncrna_class=lncRNA;transcript_id=NR_1\n\
NC_000017.11\tBestRefSeq\texon\t7700000\t7701000\t.\t+\t.\tID=exon-NR_1-1;Parent=rna-NR_1;gbkey=ncRNA;gene=LINC1;transcript_id=NR_1\n\
NC_000017.11\tBestRefSeq\tgene\t7800000\t7800500\t.\t+\t.\tID=gene-IGHV1;Dbxref=GeneID:2;Name=IGHV1;gbkey=Gene;gene=IGHV1;gene_biotype=V_segment\n\
NC_000017.11\tBestRefSeq\tV_gene_segment\t7800000\t7800500\t.\t+\t.\tID=id-IGHV1;Parent=gene-IGHV1;Dbxref=GeneID:2;gbkey=V_segment;gene=IGHV1\n\
NC_000017.11\tBestRefSeq\texon\t7800000\t7800500\t.\t+\t.\tID=exon-IGHV1-1;Parent=id-IGHV1;gbkey=V_segment;gene=IGHV1\n\
NC_000017.11\tBestRefSeq\tgene\t7900000\t7900090\t.\t+\t.\tID=gene-MIR548X;Dbxref=GeneID:3;Name=MIR548X;description=microRNA 548X;gbkey=Gene;gene=MIR548X;gene_biotype=miRNA\n\
NC_000017.11\tBestRefSeq\tprimary_transcript\t7900000\t7900090\t.\t+\t.\tID=rna-NR_3;Parent=gene-MIR548X;Dbxref=GeneID:3,Genbank:NR_3;Name=NR_3;gene=MIR548X;product=microRNA 548X;transcript_id=NR_3\n\
NC_000017.11\tBestRefSeq\texon\t7900000\t7900090\t.\t+\t.\tID=exon-NR_3-1;Parent=rna-NR_3;gene=MIR548X;transcript_id=NR_3\n\
NC_000017.11\tBestRefSeq\tmiRNA\t7900010\t7900031\t.\t+\t.\tID=rna-MIR548X-5p;Parent=rna-NR_3;Dbxref=GeneID:3,miRBase:MIMAT0022000;Name=MIR548X-5p;gbkey=ncRNA;gene=MIR548X;product=hsa-miR-548x-5p\n\
NC_000017.11\tBestRefSeq\tgene\t8000000\t8001000\t.\t+\t.\tID=gene-MISC1;Dbxref=GeneID:4;Name=MISC1;gbkey=Gene;gene=MISC1;gene_biotype=misc_RNA\n\
NC_000017.11\tBestRefSeq\ttranscript\t8000000\t8001000\t.\t+\t.\tID=rna-NR_4;Parent=gene-MISC1;Dbxref=GeneID:4,Genbank:NR_4;Name=NR_4;gbkey=misc_RNA;gene=MISC1;transcript_id=NR_4\n\
NC_000017.11\tBestRefSeq\texon\t8000000\t8001000\t.\t+\t.\tID=exon-NR_4-1;Parent=rna-NR_4;gbkey=misc_RNA;gene=MISC1;transcript_id=NR_4\n\
NC_000017.11\tBestRefSeq\tgene\t8100000\t8101000\t.\t+\t.\tID=gene-NOTYPE;Dbxref=GeneID:5;Name=NOTYPE;gbkey=Gene;gene=NOTYPE\n\
NC_000017.11\tBestRefSeq\ttranscript\t8100000\t8101000\t.\t+\t.\tID=rna-NR_5;Parent=gene-NOTYPE;Dbxref=GeneID:5;Name=NR_5;gene=NOTYPE;transcript_id=NR_5\n\
NC_000017.11\tBestRefSeq\texon\t8100000\t8101000\t.\t+\t.\tID=exon-NR_5-1;Parent=rna-NR_5;gene=NOTYPE;transcript_id=NR_5\n\
NC_000017.11\tRefSeq\tregion\t1\t83257441\t.\t+\t.\tID=NC_000017.11:1..83257441;Dbxref=taxon:9606;Name=17;chromosome=17;gbkey=Src;genome=chromosome;mol_type=genomic DNA\n\
NC_000017.11\tRefSeq\tcDNA_match\t7668402\t7687550\t.\t-\t.\tID=aln1;Target=NM_000546.6 1 2512 +;gap_count=0\n";

    /// A genome FASTA of one 60-base record per name, with its index.
    fn fasta_holding(names: &[&str]) -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genome.fa");
        let mut text = String::new();
        let mut index = String::new();
        for name in names {
            text.push_str(&format!(">{name}\n"));
            index.push_str(&format!("{name}\t60\t{}\t60\t61\n", text.len()));
            text.push_str(&"N".repeat(60));
            text.push('\n');
        }
        std::fs::write(&path, text).unwrap();
        std::fs::write(dir.path().join("genome.fa.fai"), index).unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&path).unwrap();
        (dir, fasta)
    }

    /// Without a genome FASTA a UCSC name loses its `chr` (a digit or a
    /// capital letter follows it), `chrM` and a bare `M` are `MT`, and a name
    /// that merely begins with the three letters, or with `Chr`, is kept as
    /// written.
    #[test]
    fn a_ucsc_prefix_is_dropped_and_a_name_that_begins_with_chr_is_kept() {
        let names = SequenceNames::default();
        for (written, named) in [
            ("chr1", "1"),
            ("chrX", "X"),
            ("chrM", "MT"),
            ("chrMT", "MT"),
            ("M", "MT"),
            ("1", "1"),
            ("MT", "MT"),
            ("chrUn_KI270302v1", "Un_KI270302v1"),
            ("chrI", "I"),
            ("chr2L", "2L"),
            ("chromosome_1", "chromosome_1"),
            ("chrom3", "chrom3"),
            ("chr_1", "chr_1"),
            ("chr", "chr"),
            ("Chr1", "Chr1"),
            ("ChrM", "ChrM"),
            ("Mito", "Mito"),
        ] {
            assert_eq!(names.name(written), named, "{written}");
        }
    }

    /// Under a genome FASTA the name is the one the FASTA holds. A genome
    /// whose mitochondrion is `Mito` keeps a `chrM` or an `M` as written,
    /// where `MT` would name a sequence it does not hold, and its `chrI` is
    /// `I`; a FASTA that spells its records `chr1` and `chrM`, or `1` and
    /// `MT`, gives `1` and `MT`, the names a build without one gives, and so
    /// does a FASTA that names no organelle genome at all (one chromosome of
    /// a genome); a sequence named `chromosome_1` keeps its name; and a name
    /// the three letters lead in lower case is shortened only where the
    /// FASTA holds the shortened name and not the written one.
    #[test]
    fn a_sequence_is_named_as_the_genome_fasta_holds_it() {
        let named = |records: &[&str], cases: &[(&str, &str)]| {
            let (_dir, fasta) = fasta_holding(records);
            let names = SequenceNames::of_fasta(&fasta, None);
            for (written, named) in cases {
                assert_eq!(names.name(written), *named, "{written} under {records:?}");
            }
        };
        named(
            &["I", "II", "Mito"],
            &[
                ("chrI", "I"),
                ("II", "II"),
                ("chrM", "chrM"),
                ("M", "M"),
                ("Mito", "Mito"),
                ("chrMito", "Mito"),
                ("MT", "MT"),
            ],
        );
        named(
            &["2L", "mitochondrion_genome"],
            &[("chr2L", "2L"), ("chrM", "chrM"), ("M", "M")],
        );
        named(
            &["chr1", "chrM"],
            &[("chr1", "1"), ("chrM", "MT"), ("M", "MT"), ("1", "1")],
        );
        named(
            &["1", "MT"],
            &[("chr1", "1"), ("chrM", "MT"), ("chrMT", "MT"), ("M", "MT")],
        );
        named(&["21"], &[("chr21", "21"), ("chrM", "MT"), ("M", "MT")]);
        named(
            &["chromosome_1", "chromosome_2"],
            &[("chromosome_1", "chromosome_1"), ("chrM", "MT")],
        );
        named(
            &["chromosome_1", "Pt"],
            &[("chromosome_1", "chromosome_1"), ("chrM", "chrM")],
        );
        named(&["scaffold_1", "x"], &[("chrx", "x")]);
        named(&["scaffold_1", "chrx"], &[("chrx", "chrx")]);
        named(&["x", "chrx"], &[("chrx", "chrx")]);
    }

    /// A genome FASTA of one record `name` of `length` unknown bases, with
    /// its index.
    fn fasta_of_length(name: &str, length: u64) -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genome.fa");
        let header = format!(">{name}\n");
        let mut text = header.clone();
        let mut left = length;
        while left > 0 {
            let line = left.min(60);
            text.push_str(&"N".repeat(line as usize));
            text.push('\n');
            left -= line;
        }
        std::fs::write(&path, text).unwrap();
        std::fs::write(
            dir.path().join("genome.fa.fai"),
            format!("{name}\t{length}\t{}\t60\t61\n", header.len()),
        )
        .unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&path).unwrap();
        (dir, fasta)
    }

    /// On GRCh37 a `chrM` or an `M` is `MT` only over a FASTA record of the
    /// assembly's length, whatever name holds it, or over a FASTA that gives
    /// no record for `MT` and names no organelle genome (the test after
    /// this one). It is kept as written over a record of another length:
    /// also where the run's synonym table gives that record for `MT`, as
    /// the table of Ensembl's GRCh37 cache links `MT` and `chrM`, and where
    /// the record is named `MT` itself. A name kept so is one the reader
    /// keeps apart from `MT`. Another assembly, and a naming that states
    /// none, reads either record as `MT`.
    #[test]
    fn on_grch37_a_chrm_is_mt_only_over_a_record_of_the_assemblys_length() {
        let assembly_length = vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH;
        let named = |fasta: &vep_fasta::IndexedFasta, assembly: Option<&str>| {
            let names = SequenceNames::of_fasta(fasta, assembly);
            ["chrM", "M", "MT", "chrMT"].map(|written| {
                let name = names.name(written).to_string();
                let apart = names.keeps_apart_from_mt(&name);
                (name, apart)
            })
        };
        let all_mt = ["MT", "MT", "MT", "MT"].map(|name| (name.to_string(), false));
        let kept = [("chrM", true), ("M", true), ("MT", false), ("MT", false)]
            .map(|(name, apart)| (name.to_string(), apart));

        for held in ["chrM", "M", "MT", "chrMT"] {
            let (_dir, fasta) = fasta_of_length(held, assembly_length);
            assert_eq!(named(&fasta, Some("GRCh37")), all_mt, "{held}");
        }
        for held in ["chrM", "M", "MT"] {
            let (_dir, fasta) = fasta_of_length(held, assembly_length + 2);
            assert_eq!(named(&fasta, Some("GRCh37")), kept, "{held}");
            assert_eq!(named(&fasta, Some("GRCh38")), all_mt, "{held}");
            assert_eq!(named(&fasta, None), all_mt, "{held}");
        }

        // The run's handle: the cache's synonyms, then the assembly.
        let (_dir, mut fasta) = fasta_of_length("chrM", assembly_length + 2);
        fasta.add_synonyms(&vep_core::synonyms::ChromosomeSynonyms::parse("MT\tchrM\n"));
        fasta.set_assembly("GRCh37");
        assert_eq!(named(&fasta, Some("GRCh37")), kept);

        // Both sequences in one FASTA: `MT` is read from the record of the
        // assembly's length.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genome.fa");
        let mut text = String::new();
        let mut index = String::new();
        for (name, length) in [("chrM", assembly_length + 2), ("chrMT", assembly_length)] {
            text.push_str(&format!(">{name}\n"));
            index.push_str(&format!("{name}\t{length}\t{}\t60\t61\n", text.len()));
            let mut left = length;
            while left > 0 {
                let line = left.min(60);
                text.push_str(&"N".repeat(line as usize));
                text.push('\n');
                left -= line;
            }
        }
        std::fs::write(&path, text).unwrap();
        std::fs::write(dir.path().join("genome.fa.fai"), index).unwrap();
        let mut both = vep_fasta::IndexedFasta::from_path(&path).unwrap();
        both.set_assembly("GRCh37");
        assert_eq!(named(&both, Some("GRCh37")), all_mt);
    }

    /// On GRCh37, over a FASTA that gives no record for `MT`, a `chrM` or an
    /// `M` is `MT` where the FASTA names no organelle genome, as one that
    /// holds `21` alone, and is kept as written, apart from `MT`, where it
    /// names an organelle genome of its own, as one that holds `Mito`. The
    /// record's length decides neither.
    #[test]
    fn on_grch37_a_fasta_with_no_record_for_mt_gives_mt_only_where_it_names_no_organelle_genome() {
        let assembly_length = vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH;
        let named = |held: &str, length: u64| {
            let (_dir, fasta) = fasta_of_length(held, length);
            assert_eq!(fasta.sequence_length("MT"), None, "{held}");
            let names = SequenceNames::of_fasta(&fasta, Some("GRCh37"));
            ["chrM", "M", "MT", "chrMT"].map(|written| {
                let name = names.name(written).to_string();
                let apart = names.keeps_apart_from_mt(&name);
                (name, apart)
            })
        };
        let all_mt = ["MT", "MT", "MT", "MT"].map(|name| (name.to_string(), false));
        let kept = [("chrM", true), ("M", true), ("MT", false), ("MT", false)]
            .map(|(name, apart)| (name.to_string(), apart));
        for length in [assembly_length, assembly_length + 2] {
            assert_eq!(named("21", length), all_mt, "21 of {length} bases");
            assert_eq!(named("Mito", length), kept, "Mito of {length} bases");
        }
    }

    /// A record the run's synonym table gives for `MT` counts by its length
    /// whatever its name: on GRCh37 a `chrM` or an `M` is kept as written,
    /// apart from `MT`, over a record of another length that is named like
    /// no organelle genome, and is `MT` over one of the assembly's length.
    /// Another assembly reads either record as `MT`.
    #[test]
    fn on_grch37_a_record_the_synonym_table_gives_for_mt_counts_by_its_length() {
        let assembly_length = vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH;
        let under_the_table = |length: u64, assembly: &str| {
            let (dir, mut fasta) = fasta_of_length("made_up_accession", length);
            fasta.add_synonyms(&vep_core::synonyms::ChromosomeSynonyms::parse(
                "MT\tmade_up_accession\n",
            ));
            fasta.set_assembly(assembly);
            assert_eq!(fasta.sequence_length("MT"), Some(length));
            (dir, SequenceNames::of_fasta(&fasta, Some(assembly)))
        };
        let (_dir, names) = under_the_table(assembly_length + 2, "GRCh37");
        for written in ["chrM", "M"] {
            assert_eq!(names.name(written), written);
            assert!(names.keeps_apart_from_mt(written), "{written}");
        }
        assert_eq!(names.name("MT"), "MT");
        assert!(!names.keeps_apart_from_mt("MT"));
        for (length, assembly) in [(assembly_length, "GRCh37"), (assembly_length + 2, "GRCh38")] {
            let (_dir, names) = under_the_table(length, assembly);
            for written in ["chrM", "M"] {
                assert_eq!(names.name(written), "MT", "{written} on {assembly}");
                assert!(
                    !names.keeps_apart_from_mt(written),
                    "{written} on {assembly}"
                );
            }
        }
    }

    /// The model reader names each record's sequence by that rule: of two
    /// genes written on `chromosome_1` and `chrM`, the first keeps its name
    /// without a FASTA and under one that holds it, and the second is on `MT`
    /// without a FASTA and under one that holds a mitochondrion of that
    /// family, and on `chrM` under one whose mitochondrion is `Mito`.
    #[test]
    fn the_model_reader_names_a_sequence_as_the_fasta_holds_it() {
        let lines = format!(
            "##gff-version 3\n{}{}",
            TPTE_GENE.replacen("21", "chromosome_1", 1),
            TPTE_GENE
                .replacen("21", "chrM", 1)
                .replace("ENSG00000274391", "ENSG00000000001")
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gff3");
        std::fs::write(&path, lines).unwrap();
        let read = |names: &SequenceNames| {
            let genes = parse_gff3(&path, Reading::GffSource, names).unwrap();
            (
                genes["ENSG00000274391"].chr.clone(),
                genes["ENSG00000000001"].chr.clone(),
            )
        };
        let on = |first: &str, second: &str| (first.to_string(), second.to_string());
        assert_eq!(read(&SequenceNames::default()), on("chromosome_1", "MT"));
        let (_dir, own_names) = fasta_holding(&["chromosome_1", "Mito"]);
        assert_eq!(
            read(&SequenceNames::of_fasta(&own_names, None)),
            on("chromosome_1", "chrM")
        );
        let (_dir, as_written) = fasta_holding(&["chromosome_1", "chrM"]);
        assert_eq!(
            read(&SequenceNames::of_fasta(&as_written, None)),
            on("chromosome_1", "MT")
        );
    }

    /// The `gene:` and `transcript:` prefixes come off in any letter case,
    /// as `BaseGXF.pm` strips them (`s/^(gene|transcript)://i`), and any
    /// other identifier is kept whole.
    #[test]
    fn gxf_prefix_is_stripped_in_any_case() {
        assert_eq!(strip_gxf_prefix("gene:ENSG00000123"), "ENSG00000123");
        assert_eq!(strip_gxf_prefix("transcript:ENST00000456"), "ENST00000456");
        assert_eq!(strip_gxf_prefix("Transcript:ENST00000456"), "ENST00000456");
        assert_eq!(strip_gxf_prefix("ENSG00000123"), "ENSG00000123");
        assert_eq!(strip_gxf_prefix("rna-NM_000546.6"), "rna-NM_000546.6");
    }

    /// The Ensembl `gene:` and `transcript:` prefixes are read in any letter
    /// case wherever they are tested, the case `BaseGXF.pm` strips them in
    /// (`s/^(gene|transcript)://i`, lines 465 and 657): under
    /// [`Reading::EnsemblCache`] a `Gene:G1` line of a type naming no gene is a
    /// gene, a `TRANSCRIPT:T2` line of a type outside the list is admitted with
    /// its `tag=basic` read as `gencode_basic`, and both identifiers lose their
    /// prefix; a word that merely opens with a prefix's letters is no prefix.
    /// Under the `--gff` reading both lines are skipped by their types, as
    /// Ensembl VEP skips them.
    #[test]
    fn the_ensembl_prefixes_are_read_in_any_letter_case() {
        let cache = Reading::EnsemblCache;
        assert!(is_gene_line(cache, "Gene:G1", "transcript"));
        assert!(is_gene_line(cache, "GENE:G1", "biological_region"));
        assert!(!is_gene_line(cache, "genes:G1", "transcript"));
        assert!(has_gxf_prefix("TRANSCRIPT:T1") && has_gxf_prefix("Gene:G1"));
        assert!(!has_gxf_prefix("transcript_id:T1") && !has_gxf_prefix("gen"));
        let lines = "21\tx\tbiological_region\t1000\t2000\t.\t+\t.\tID=Gene:G1;biotype=protein_coding\n\
21\tx\tunconfirmed_transcript\t1000\t2000\t.\t+\t.\tID=TRANSCRIPT:T2;Parent=Gene:G1;biotype=TEC;tag=basic\n\
21\tx\texon\t1000\t2000\t.\t+\t.\tParent=TRANSCRIPT:T2;Name=X2\n";
        let genes = try_parse_as(lines, Reading::EnsemblCache).unwrap();
        let tec = transcript(&genes, "T2");
        assert_eq!(
            (
                tec.gene_id.as_str(),
                tec.biotype.as_str(),
                tec.gencode_basic
            ),
            ("G1", "TEC", true)
        );
        assert_eq!(genes["G1"].transcripts.len(), 1);
        assert!(parse(lines).is_empty());
    }

    /// A record without a parent is a gene by its type alone under the `--gff`
    /// reading, the three tests of `BaseGXF.pm` `_record_is_gene` (lines
    /// 788-796), so a parentless `ID=gene:x` line of type `mRNA` is a
    /// transcript, as it is to Ensembl VEP's `--gff` source, and its `version`
    /// is read under the transcript keys. Under [`Reading::EnsemblCache`] the
    /// `gene:` prefix makes the line a gene whatever its type, Ensembl's
    /// dumper convention, so the same line is a gene with no transcript and
    /// `version` is read under the gene key alone. The type tests read the
    /// type as written: `Gene` and `rna` name no gene.
    #[test]
    fn a_parentless_gene_prefixed_mrna_is_a_transcript_under_the_gff_reading() {
        let lines = "21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=gene:x;biotype=protein_coding;transcript_version=7\n\
21\tx\texon\t1000\t2000\t.\t+\t.\tParent=gene:x;Name=X1\n";
        let genes = parse(lines);
        let tr = transcript(&genes, "x");
        assert_eq!(
            (tr.biotype.as_str(), tr.version),
            ("protein_coding", Some(7))
        );
        assert_eq!(tr.exons.len(), 1);
        assert!(!genes.contains_key("x"), "no gene record of its own");

        let genes = try_parse_as(lines, Reading::EnsemblCache).unwrap();
        let gene = &genes["x"];
        assert_eq!((gene.transcripts.len(), gene.version), (0, None));

        assert!(!is_gene_line(Reading::GffSource, "gene:x", "mRNA"));
        assert!(is_gene_line(Reading::GffSource, "x", "ncRNA_gene"));
        assert!(!is_gene_line(Reading::GffSource, "x", "Gene"));
        assert!(!is_gene_line(Reading::GffSource, "x", "rna"));
        assert_eq!(
            feature_version_keys(Reading::GffSource, false, "gene:x", "mRNA"),
            ["version", "transcript_version"]
        );
        assert_eq!(
            feature_version_keys(Reading::EnsemblCache, false, "gene:x", "mRNA"),
            ["version"]
        );
    }

    /// A gene, its transcript and one exon on sequence `chr`, three lines.
    fn three_lines(chr: &str, tag: &str) -> String {
        format!(
            "{chr}\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G{tag};biotype=protein_coding\n\
{chr}\tx\ttranscript\t1000\t2000\t.\t+\t.\tID=transcript:T{tag};Parent=gene:G{tag};biotype=protein_coding\n\
{chr}\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T{tag};Name=X{tag}\n"
        )
    }

    /// The model of two sequences with one more line after each one's own
    /// three: `on_first` is line 4 of the file and `on_later` line 8.
    fn two_sequences(on_first: &str, on_later: &str) -> String {
        format!(
            "{}{on_first}{}{on_later}",
            three_lines("20", "0"),
            three_lines("21", "1")
        )
    }

    /// A line of type `feature_type` on sequence `chr` under transcript
    /// `T<tag>` with the start and end as written.
    fn line_of(chr: &str, feature_type: &str, start: &str, end: &str, tag: &str) -> String {
        format!("{chr}\tx\t{feature_type}\t{start}\t{end}\t.\t+\t.\tParent=transcript:T{tag}\n")
    }

    const SOUND_UTR_ON_FIRST: &str =
        "20\tx\tfive_prime_UTR\t1000\t1100\t.\t+\t.\tParent=transcript:T0\n";
    const SOUND_UTR_ON_LATER: &str =
        "21\tx\tfive_prime_UTR\t1000\t1100\t.\t+\t.\tParent=transcript:T1\n";

    /// A start or end no number leads refuses the parse on a line of a type
    /// the reading skips, on a file's first sequence and on a later one,
    /// naming the file, the line, the column and the value: no htslib's index
    /// build parses such a column (`tbx.c` `tbx_parse1`), and the 1.9 build of
    /// Ensembl VEP's release passes it only at a file's first data line, under
    /// a sequence id it never set, or, for an end, by losing the later lines
    /// of its block to every query. An empty column is read on into the next,
    /// so an empty start refuses before an end that is no number and an empty
    /// end before the score `.`. An end below zero fails the index build as
    /// one no number leads does, and refuses with it.
    #[test]
    fn a_coordinate_the_index_build_fails_is_refused_on_a_skipped_line_of_any_sequence() {
        for (start, end, message) in [
            ("abc", "1100", "GFF3 start \"abc\" is not a number"),
            ("1000", "xyz", "GFF3 end \"xyz\" is not a number"),
            ("", "xyz", "GFF3 start \"\" is not a number"),
            ("1000", "", "GFF3 end \"\" is not a number"),
            ("1000", "-", "GFF3 end \"-\" is not a number"),
            ("1000", "-9", "GFF3 end \"-9\" is not a number"),
        ] {
            let on_first = line_of("20", "five_prime_UTR", start, end, "0");
            let err = try_parse(&two_sequences(&on_first, SOUND_UTR_ON_LATER)).unwrap_err();
            assert_eq!(
                format!("{err:#}")
                    .rsplit_once('/')
                    .map(|(_, tail)| tail.to_string()),
                Some(format!("t.gff3 line 4: {message}")),
                "{start:?} {end:?} on the first sequence"
            );
            let on_later = line_of("21", "five_prime_UTR", start, end, "1");
            let err = try_parse(&two_sequences(SOUND_UTR_ON_FIRST, &on_later)).unwrap_err();
            assert_eq!(
                format!("{err:#}")
                    .rsplit_once('/')
                    .map(|(_, tail)| tail.to_string()),
                Some(format!("t.gff3 line 8: {message}")),
                "{start:?} {end:?} on a later sequence"
            );
        }
        // An empty end is read on into its score, and fails below zero.
        let scored = "20\tx\tfive_prime_UTR\t1000\t\t-3\t+\t.\tParent=transcript:T0\n";
        let err = try_parse(&two_sequences(scored, SOUND_UTR_ON_LATER)).unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gff3 line 4: GFF3 end \"\" is not a number"),
            "{err:#}"
        );
        // The columns are the line's as written: one that opens with a tab
        // holds its type where the index build reads a start.
        let shifted = format!("\t{SOUND_UTR_ON_FIRST}");
        let err = try_parse(&two_sequences(&shifted, SOUND_UTR_ON_LATER)).unwrap_err();
        assert!(
            format!("{err:#}")
                .ends_with("t.gff3 line 4: GFF3 start \"five_prime_UTR\" is not a number"),
            "{err:#}"
        );
    }

    /// Nothing else of a skipped line is read, as Ensembl VEP's source reads a
    /// record's raw start, for its loop's bound, and its type, and nothing
    /// else of one it skips (`BaseGXF.pm:217-245`): a `five_prime_UTR` line
    /// whose start is 0, whose end precedes its start, whose seqid could name
    /// no cache directory, or whose start or end a number leads without being
    /// one (a sign, a blank, trailing text, a hexadecimal prefix, an empty or
    /// blank start before its end), the end not below zero, refuses nothing
    /// on either sequence, the release's `tabix` indexing each. The same
    /// columns on an exon line, which the tree keeps, refuse the parse.
    #[test]
    fn a_skipped_lines_seqid_and_number_led_coordinates_are_not_read() {
        for (chr, start, end) in [
            ("20", "0", "1100"),
            ("20", "1100", "1000"),
            ("20/x", "1000", "1100"),
            ("20", "1000x", "1100"),
            ("20", "-5", "1100"),
            ("20", " 1000", "1100"),
            ("20", "0x3E8", "1100"),
            ("20", "+1000", "1100"),
            ("20", "001000", "001100"),
            ("20", "1000", "1100x"),
            ("20", "", "1100"),
            ("20", " ", "1100"),
            ("20", "1000", "-0"),
        ] {
            for (first, later) in [
                (
                    line_of(chr, "five_prime_UTR", start, end, "0"),
                    SOUND_UTR_ON_LATER.to_string(),
                ),
                (
                    SOUND_UTR_ON_FIRST.to_string(),
                    line_of(&chr.replace("20", "21"), "five_prime_UTR", start, end, "1"),
                ),
            ] {
                let genes = try_parse(&two_sequences(&first, &later))
                    .unwrap_or_else(|err| panic!("{chr:?} {start:?} {end:?}: {err:#}"));
                assert_eq!(transcript(&genes, "T0").exons.len(), 1);
                assert_eq!(transcript(&genes, "T1").exons.len(), 1);
            }
        }
        // An empty end is read on into a score that is a number.
        let scored = "20\tx\tfive_prime_UTR\t1000\t\t7\t+\t.\tParent=transcript:T0\n";
        let genes = try_parse(&two_sequences(scored, SOUND_UTR_ON_LATER)).unwrap();
        assert_eq!(transcript(&genes, "T0").exons.len(), 1);
        for (chr, start, end, message) in [
            ("21", "abc", "2000", "GFF3 start \"abc\" is not a number"),
            (
                "21",
                "1000x",
                "2000",
                "GFF3 start \"1000x\" is not a number",
            ),
            ("21", "-5", "2000", "GFF3 start \"-5\" is not a number"),
            ("21", "", "2000", "GFF3 start \"\" is not a number"),
            (
                "21",
                "0",
                "2000",
                "GFF3 coordinates 0-2000 are not 1-based and ordered",
            ),
            (
                "21",
                "2000",
                "1000",
                "GFF3 coordinates 2000-1000 are not 1-based and ordered",
            ),
            ("21/x", "1000", "2000", "invalid chromosome name '21/x'"),
        ] {
            let exon = line_of(chr, "exon", start, end, "1");
            let err = try_parse(&two_sequences(SOUND_UTR_ON_FIRST, &exon)).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gff3 line 8: {message}")),
                "{err:#}"
            );
        }
    }

    /// A line of four to eight columns is skipped whatever its type and is
    /// read no further than its start and its end: one the index build fails,
    /// a blank one ending the line included, refuses the parse on either
    /// sequence, as on a line of nine columns, and any other spelling
    /// refuses nothing. The columns are the line's as written, so a line of
    /// a skipped type whose seqid is empty is read at its own start and end.
    #[test]
    fn a_short_line_is_refused_for_a_coordinate_the_index_build_fails_alone() {
        for (columns, message) in [
            (
                "\tx\texon\tabc\t1100\t.\t+\t.\n",
                "GFF3 start \"abc\" is not a number",
            ),
            (
                "\tx\tfive_prime_UTR\t1000\txyz\n",
                "GFF3 end \"xyz\" is not a number",
            ),
            (
                "\tx\tfive_prime_UTR\t-\n",
                "GFF3 start \"-\" is not a number",
            ),
            ("\tx\tfive_prime_UTR\t\n", "GFF3 start \"\" is not a number"),
            (
                "\tx\tfive_prime_UTR\t1000\t \t\n",
                "GFF3 end \" \" is not a number",
            ),
            (
                "\tx\tfive_prime_UTR\t1000\t\r\n",
                "GFF3 end \"\" is not a number",
            ),
            (
                "\tx\tfive_prime_UTR\t1000\t-9\n",
                "GFF3 end \"-9\" is not a number",
            ),
        ] {
            let on_first = format!("20{columns}");
            let err = try_parse(&two_sequences(&on_first, SOUND_UTR_ON_LATER)).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gff3 line 4: {message}")),
                "{err:#}"
            );
            let on_later = format!("21{columns}");
            let err = try_parse(&two_sequences(SOUND_UTR_ON_FIRST, &on_later)).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gff3 line 8: {message}")),
                "{err:#}"
            );
        }
        for columns in [
            "\tx\texon\t0\t1100\t.\t+\t.\n",
            "\tx\texon\t1100\t1000\t.\t+\t.\n",
            "\tx\texon\t1000x\t1100\n",
            "\tx\texon\t1000\n",
        ] {
            for (first, later) in [
                (format!("20{columns}"), SOUND_UTR_ON_LATER.to_string()),
                (SOUND_UTR_ON_FIRST.to_string(), format!("21{columns}")),
            ] {
                let genes = try_parse(&two_sequences(&first, &later))
                    .unwrap_or_else(|err| panic!("{columns:?}: {err:#}"));
                assert_eq!(transcript(&genes, "T0").exons.len(), 1);
                assert_eq!(transcript(&genes, "T1").exons.len(), 1);
            }
        }
        // A line of a skipped type whose seqid is empty is read at its columns
        // as written.
        let unnamed = line_of("", "five_prime_UTR", "1000", "1100", "0");
        assert!(try_parse(&two_sequences(&unnamed, SOUND_UTR_ON_LATER)).is_ok());
        let unnamed = line_of("", "five_prime_UTR", "abc", "1100", "0");
        let err = try_parse(&two_sequences(&unnamed, SOUND_UTR_ON_LATER)).unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gff3 line 4: GFF3 start \"abc\" is not a number"),
            "{err:#}"
        );
    }

    /// A line that is no comment and ends before its fourth column is a
    /// record without its start: the index build fails it, and the parse
    /// refuses it on either sequence, whatever its type, naming the column it
    /// ends at. A blank line, empty or of white space alone, is skipped.
    #[test]
    fn a_line_that_ends_before_its_start_column_is_refused() {
        for (columns, message) in [
            (
                "\tx\texon\n",
                "GFF3 line ends at column 3, before its start",
            ),
            ("\tx\t \n", "GFF3 line ends at column 3, before its start"),
            ("\tx\n", "GFF3 line ends at column 2, before its start"),
            ("\n", "GFF3 line ends at column 1, before its start"),
        ] {
            let on_first = format!("20{columns}");
            let err = try_parse(&two_sequences(&on_first, SOUND_UTR_ON_LATER)).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gff3 line 4: {message}")),
                "{err:#}"
            );
            let on_later = format!("21{columns}");
            let err = try_parse(&two_sequences(SOUND_UTR_ON_FIRST, &on_later)).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gff3 line 8: {message}")),
                "{err:#}"
            );
        }
        for blank in ["\n", "  \n", "\t\t\t\n", "\r\n"] {
            let genes = try_parse(&two_sequences(blank, SOUND_UTR_ON_LATER))
                .unwrap_or_else(|err| panic!("{blank:?}: {err:#}"));
            assert_eq!(transcript(&genes, "T0").exons.len(), 1);
            assert_eq!(transcript(&genes, "T1").exons.len(), 1);
        }
    }

    /// A byte-order mark at the start of a line changes nothing: a model
    /// with a mark before each of its lines, a comment, a blank line and
    /// every record among them, one mark behind white space, parses to the
    /// genes of the same lines without the marks, and a marked record cut
    /// short refuses as it does unmarked, on a later line as on the first.
    #[test]
    fn a_byte_order_mark_at_the_start_of_a_line_changes_nothing() {
        let model = format!(
            "##gff-version 3\n\n  # x\n{}",
            two_sequences(SOUND_UTR_ON_FIRST, SOUND_UTR_ON_LATER)
        );
        let marked: String = model.lines().map(|l| format!("\u{feff}{l}\n")).collect();
        let marked = marked.replacen("\u{feff}  # x", "  \u{feff}# x", 1);
        assert_eq!(marked.matches('\u{feff}').count(), 11);
        let (plain, marked) = (parse(&model), parse(&marked));
        assert_eq!(marked.len(), 2);
        for gene in ["G0", "G1"] {
            assert_eq!(format!("{:?}", marked[gene]), format!("{:?}", plain[gene]));
        }
        for (text, message) in [
            (
                format!("\u{feff}20\tx\n{model}"),
                "t.gff3 line 1: GFF3 line ends at column 2, before its start",
            ),
            (
                format!("{model}\u{feff}20\tx\n"),
                "t.gff3 line 12: GFF3 line ends at column 2, before its start",
            ),
        ] {
            let err = try_parse(&text).unwrap_err();
            assert!(format!("{err:#}").ends_with(message), "{err:#}");
        }
    }

    /// The marks that open a line are taken off, alone or behind white
    /// space, however many stand there; a mark further on is the line's
    /// text, and a line without one is returned as it is.
    #[test]
    fn the_marks_that_open_a_line_are_taken_off() {
        for (line, unmarked) in [
            ("\u{feff}21\tx", "21\tx"),
            ("\u{feff}\u{feff}##x", "##x"),
            ("  \t\u{feff}21", "  \t21"),
            ("\u{feff}", ""),
            ("21\t\u{feff}x", "21\t\u{feff}x"),
            ("x\u{feff}", "x\u{feff}"),
            ("21\tx", "21\tx"),
            ("", ""),
        ] {
            assert_eq!(
                without_byte_order_marks(line.as_bytes().to_vec()),
                unmarked.as_bytes(),
                "{line:?}"
            );
        }
    }

    /// Nothing is read from the first line that begins with `>` and has a
    /// sequence behind it to the file's end, with a `##FASTA` directive
    /// before it or without one: the sequence section the GFF3 specification
    /// gives a file (version 1.26, the `##FASTA` directive and the directive
    /// a `>` line implies). What stands behind the header decides, a blank
    /// line and a comment between the two or not, so the header is of one
    /// column or of several, numbers in its fourth and fifth or not, nine
    /// columns that read as a record among them. A line there that would
    /// refuse as a record refuses nothing; a record before the section is
    /// read, and refuses when it ends before its start.
    #[test]
    fn a_sequence_section_is_not_read() {
        let body = "NNNNNNNN\n20\tx\n20\tx\texon\tabc\t2000\t.\t+\t.\tParent=transcript:T0\n";
        for header in [
            ">20 placeholder\n",
            ">20\tplaceholder\tof\tfour columns\n",
            ">20\tplaceholder\tof\tfive\tcolumns\n",
            ">20\tplaceholder\tof\t100\t-9\n",
            ">20\tplaceholder\tof\t100\n",
            ">20\tplaceholder\tof\t100\t200\n",
            ">20\tplaceholder\tof\t100\t200\t.\t+\t.\n",
            ">20\tplaceholder\tof\t100\t200\t.\t+\t.\tnine columns\n",
            ">20 placeholder\n\n# x\n",
        ] {
            for directive in ["##FASTA\n", ""] {
                let section = format!("{directive}{header}{body}");
                let text = format!("{}{section}", three_lines("20", "0"));
                let genes = try_parse(&text).unwrap_or_else(|err| panic!("{section:?}: {err:#}"));
                assert_eq!(transcript(&genes, "T0").exons.len(), 1);
                let text = format!("{}20\tx\n{section}", three_lines("20", "0"));
                let err = try_parse(&text).unwrap_err();
                assert!(
                    format!("{err:#}")
                        .ends_with("t.gff3 line 4: GFF3 line ends at column 2, before its start"),
                    "{err:#}"
                );
            }
        }
    }

    /// A `##FASTA` directive is a comment and opens no sequence section: the
    /// records behind a bare one are read as behind any comment, between
    /// two sequences as on the file's last line.
    #[test]
    fn a_fasta_directive_alone_opens_no_sequence_section() {
        let genes = parse(&two_sequences("##FASTA\n", SOUND_UTR_ON_LATER));
        assert_eq!(genes.len(), 2);
        assert_eq!(transcript(&genes, "T1").exons.len(), 1);
        let genes = parse(&two_sequences(SOUND_UTR_ON_FIRST, "##FASTA\n"));
        assert_eq!(genes.len(), 2);
    }

    /// A line that begins with `>` with a record behind it, or with nothing,
    /// opens no sequence section and is read as any other line, a blank
    /// line and a comment between the two or not. One that reads as a
    /// record is a record whose seqid begins with `>`, read by its type,
    /// kept or skipped: the index build files it as a sequence of its own
    /// and Ensembl VEP reads the records after it, so the parse reads them
    /// too. One of four to eight columns is skipped where the index build
    /// parses its start and its end, and refused where it fails one; one of
    /// fewer than four is a record cut short, refused by name. A record or
    /// a comment longer than a line's limit behind such a line refuses as
    /// on any line, and a record behind lines passed over refuses at its
    /// own line.
    #[test]
    fn a_line_beginning_with_an_angle_bracket_before_a_record_is_a_line_like_any_other() {
        let utr = line_of(">20", "five_prime_UTR", "1000", "1100", "0");
        let genes = try_parse(&two_sequences(&utr, SOUND_UTR_ON_LATER)).unwrap();
        assert_eq!(genes.len(), 2);
        assert_eq!(transcript(&genes, "T1").exons.len(), 1);
        let exon = line_of(">20", "exon", "1500", "1600", "0");
        let genes = try_parse(&two_sequences(&exon, SOUND_UTR_ON_LATER)).unwrap();
        assert_eq!(genes.len(), 2);
        assert_eq!(transcript(&genes, "T0").exons.len(), 2);
        for skipped in [
            ">20\tx\texon\t1000\n",
            ">20\tx\texon\t1000\t1100\n",
            ">20\tx\texon\t1000\t1100\t.\t+\t.\n",
            ">20\tx\tfive_prime_UTR\t100x\t1100\t.\t+\t.\tParent=transcript:T0\n",
            ">20\tx\texon\t1000\n\n# x\n",
        ] {
            let genes = try_parse(&two_sequences(skipped, SOUND_UTR_ON_LATER))
                .unwrap_or_else(|err| panic!("{skipped:?}: {err:#}"));
            assert_eq!(genes.len(), 2, "{skipped:?}");
            assert_eq!(transcript(&genes, "T1").exons.len(), 1);
            let at_the_end = format!("{}{skipped}", three_lines("20", "0"));
            assert_eq!(parse(&at_the_end).len(), 1, "{skipped:?}");
        }
        for (refused, message) in [
            (
                ">20 placeholder\n",
                "GFF3 line ends at column 1, before its start",
            ),
            (
                ">20\tx\texon\n",
                "GFF3 line ends at column 3, before its start",
            ),
            (">20\tx\texon\tabc\n", "GFF3 start \"abc\" is not a number"),
        ] {
            for text in [
                two_sequences(refused, SOUND_UTR_ON_LATER),
                format!("{}{refused}", three_lines("20", "0")),
            ] {
                let err = try_parse(&text).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gff3 line 4: {message}")),
                    "{err:#}"
                );
            }
        }
        let padding: String = std::iter::repeat_n('x', MAX_LINE_BYTES).collect();
        let long = line_of("21", "exon", "1500", "1600", &format!("1;Note={padding}"));
        let text = format!("{}>20\tx\texon\t1000\n{long}", three_lines("20", "0"));
        let err = try_parse(&text).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "t.gff3 line 5 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
        let comment = format!("#{padding}x");
        let text = format!(
            "{}>20\tx\texon\t1000\n{comment}\n{SOUND_UTR_ON_LATER}",
            three_lines("20", "0")
        );
        let err = try_parse(&text).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "t.gff3 line 5 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
        let text = format!(
            "{}>20\tx\texon\t1000\n\n# x\n{}",
            three_lines("20", "0"),
            line_of("21", "exon", "0", "5", "9")
        );
        let err = try_parse(&text).unwrap_err();
        assert!(
            format!("{err:#}")
                .ends_with("t.gff3 line 7: GFF3 coordinates 0-5 are not 1-based and ordered"),
            "{err:#}"
        );
    }

    /// The line behind a line that begins with `>` decides by what the index
    /// build does with it. One that build places opens no section, a record
    /// by the test of a line inside a section or not: four, five or eight
    /// columns with numbers, and nine whose start only leads with a number,
    /// reads below zero or is empty before a numeric end. Both lines are
    /// then read as any other and the records behind them too, the `>` line
    /// being of five columns or one that reads as a record; a `>` line of
    /// fewer than four columns is refused there, a record cut short. A line
    /// of four or more columns whose start or end the index build fails is
    /// a sequence's: the section opens, the line refuses nothing, and a
    /// record behind it is refused by name. A line longer than a line's
    /// limit in which a tab stands counts as placed, its fourth column or
    /// its ninth opening past the limit, a comment among them: it is
    /// refused for its length, as on any line, and the rest of it is passed
    /// over.
    #[test]
    fn the_line_behind_an_angle_bracket_line_decides_by_what_the_index_build_places() {
        let placed = [
            "21\tplaceholder\tx\t1\n",
            "21\tplaceholder\tx\t1\t2\n",
            "21\tplaceholder\tx\t1\t2\t.\t+\t.\n",
            "21\tx\tbiological_region\t1x\t2\t.\t+\t.\tID=region:9\n",
            "21\tx\tbiological_region\t-5\t2\t.\t+\t.\tID=region:9\n",
            "21\tx\tbiological_region\t\t2\t.\t+\t.\tID=region:9\n",
        ];
        for line in placed {
            for first in [
                ">20\tplaceholder\tx\t100\t200\n",
                ">20\tx\tbiological_region\t100\t200\t.\t+\t.\tID=region:1\n",
            ] {
                let both = format!("{first}{line}");
                let genes = try_parse(&two_sequences(&both, SOUND_UTR_ON_LATER))
                    .unwrap_or_else(|err| panic!("{both:?}: {err:#}"));
                assert_eq!(genes.len(), 2, "{both:?}");
                assert_eq!(transcript(&genes, "T1").exons.len(), 1);
                let at_the_end = format!("{}{both}", three_lines("20", "0"));
                assert_eq!(parse(&at_the_end).len(), 1, "{both:?}");
            }
            let short = format!(">20 placeholder\n{line}");
            let err = try_parse(&two_sequences(&short, SOUND_UTR_ON_LATER)).unwrap_err();
            assert!(
                format!("{err:#}")
                    .ends_with("t.gff3 line 4: GFF3 line ends at column 1, before its start"),
                "{short:?}: {err:#}"
            );
        }
        for failed in [
            "21\tplaceholder\tx\tabc\n",
            "21\tplaceholder\tx\t1\tabc\n",
            "21\tplaceholder\tx\t1\t-2\n",
            "21\tx\tgene\tabc\t2\t.\t+\t.\tID=gene:G9\n",
        ] {
            let text = format!(
                "{}>20\tplaceholder\tx\t100\t200\n{failed}",
                three_lines("20", "0")
            );
            assert_eq!(parse(&text).len(), 1, "{failed:?}");
            let err = try_parse(&format!("{text}{SOUND_UTR_ON_FIRST}")).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(
                    "t.gff3 line 6: GFF3 annotation after the sequence section, which opens at line 4"
                ),
                "{failed:?}: {err:#}"
            );
        }
        let filler: String = std::iter::repeat_n('x', MAX_LINE_BYTES).collect();
        for long in [
            format!("21\tx\tgene\t1\t2\t.\t+\t{filler}\tID=gene:G9"),
            format!("{filler}x\tx\tgene\t1\t2"),
            format!("{filler}x\tx"),
            format!("#\tx{filler}"),
        ] {
            let text = format!(
                "{}>20\tplaceholder\tx\t100\t200\n{long}\n",
                three_lines("20", "0")
            );
            let err = try_parse(&text).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!(
                    "t.gff3 line 5 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
                )),
                "{err:#}"
            );
        }
    }

    /// No sequence section opens in a tab-delimited file whose every data
    /// line the index build places, so the file is read line by line as one
    /// without a `>` line is. The cases are drawn from seventeen lines:
    /// eight that build places, each written on a sequence and on one whose
    /// name begins with `>`, and a comment. The eight are of four, five and
    /// eight columns with numbers and of nine: a kept record, a skipped
    /// one, and skipped ones whose start only leads with a number, reads
    /// below zero or is empty before a numeric end. A run is one of the
    /// seventeen, or two of them straight behind one another or with a
    /// comment between: 595 runs. Each stands between the records of two
    /// sequences and, again, behind them at the file's end: 1,190 files. In
    /// each, no `>` line has a sequence behind it, and the parse reads what
    /// it reads of the file's kept records alone.
    #[test]
    fn no_section_opens_where_the_index_build_places_every_data_line() {
        const KEPT: &str = "\tx\texon\t1500\t1600\t.\t+\t.\tParent=transcript:T0";
        let placed = [
            "\tplaceholder\tx\t100",
            "\tplaceholder\tx\t100\t200",
            "\tplaceholder\tx\t100\t200\t.\t+\t.",
            "\tx\tbiological_region\t100\t200\t.\t+\t.\tID=region:1",
            "\tx\tbiological_region\t100x\t200\t.\t+\t.\tID=region:1",
            "\tx\tbiological_region\t-5\t200\t.\t+\t.\tID=region:1",
            "\tx\tbiological_region\t\t200\t.\t+\t.\tID=region:1",
            KEPT,
        ];
        const COMMENT: &str = "# x\n";
        let mut lines = vec![COMMENT.to_string()];
        for seqid in ["20", ">20"] {
            lines.extend(placed.iter().map(|columns| format!("{seqid}{columns}\n")));
        }
        let mut runs = lines.clone();
        for one in &lines {
            for other in &lines {
                runs.push(format!("{one}{other}"));
                runs.push(format!("{one}{COMMENT}{other}"));
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gff3");
        let parsed = |text: &str| {
            std::fs::write(&path, text).unwrap();
            let genes = parse_gff3(&path, Reading::GffSource, &SequenceNames::default())
                .unwrap_or_else(|err| panic!("{text:?}: {err:#}"));
            let mut genes: Vec<_> = genes.into_iter().collect();
            genes.sort_by(|(one, _), (other, _)| one.cmp(other));
            format!("{genes:?}")
        };
        let mut parsed_of_kept: HashMap<String, String> = HashMap::new();
        let (first, later) = (three_lines("20", "0"), three_lines("21", "1"));
        let mut files = 0;
        for run in &runs {
            let kept: String = run
                .split_inclusive('\n')
                .filter(|line| line.trim_end().ends_with(KEPT))
                .collect();
            for (text, kept) in [
                (
                    format!("{first}{run}{later}"),
                    format!("{first}{kept}{later}"),
                ),
                (
                    format!("{first}{later}{run}"),
                    format!("{first}{later}{kept}"),
                ),
            ] {
                let mut behind = text.as_str();
                while let Some((line, rest)) = behind.split_once('\n') {
                    behind = rest;
                    if !line.starts_with('>') {
                        continue;
                    }
                    let (mut read, mut ahead) = (0, None);
                    let mut reader = std::io::Cursor::new(behind);
                    let follows = sequence_follows(&mut reader, &path, &mut read, &mut ahead);
                    assert!(!follows.unwrap(), "{text:?}: behind {line:?}");
                }
                let expected = parsed_of_kept
                    .entry(kept)
                    .or_insert_with_key(|kept| parsed(kept));
                assert_eq!(parsed(&text), *expected, "{text:?}");
                files += 1;
            }
        }
        assert_eq!(files, 1_190);
    }

    /// A line inside a sequence section that reads as a record, nine or more
    /// columns of which the fourth and the fifth are whole numbers, refuses
    /// the parse, naming its line and the line the section opens at, the
    /// `>` line behind a `##FASTA` directive as one without it, and whatever
    /// the record's type: the section ends the annotation, and the record
    /// would go unread without a word. A comment there, behind a mark or
    /// any white space, a `>` header whatever its columns, a line of fewer
    /// columns, a line whose start is no number and a sequence longer than
    /// a line's limit, straight behind the header or further on, refuse
    /// nothing, and a line behind that sequence is still counted. A line
    /// longer than the limit is held to the limit: a ninth column that
    /// opens past it makes no record of the line.
    #[test]
    fn annotation_inside_the_sequence_section_is_refused() {
        for (opening, opened, record_line) in [
            (">20 placeholder\nNNNNNNNN\n", 4, 6),
            ("##FASTA\n>20 placeholder\nNNNNNNNN\n", 5, 7),
            (
                ">20\tplaceholder\tof\tfour columns\n\u{feff}# x\n\nNNNNNNNN\n",
                4,
                8,
            ),
        ] {
            let err = try_parse(&two_sequences(opening, SOUND_UTR_ON_LATER)).unwrap_err();
            let refusal = format!(
                "t.gff3 line {record_line}: GFF3 annotation after the sequence section, which opens at line {opened}"
            );
            assert!(format!("{err:#}").ends_with(&refusal), "{err:#}");
        }
        let skipped_type = format!(
            "{}>20\nNNNNNNNN\n{SOUND_UTR_ON_FIRST}",
            three_lines("20", "0")
        );
        let err = try_parse(&skipped_type).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 6: GFF3 annotation after the sequence section, which opens at line 4"
            ),
            "{err:#}"
        );
        let unread = "#21\tx\tgene\t1\t2\t.\t+\t.\tID=gene:G9\n\
\u{feff} #21\tx\tgene\t1\t2\t.\t+\t.\tID=gene:G9\n\
\u{b}#21\tx\tgene\t1\t2\t.\t+\t.\tID=gene:G9\n\
\u{a0}#21\tx\tgene\t1\t2\t.\t+\t.\tID=gene:G9\n\
>21\tx\tgene\t1\t2\t.\t+\t.\tID=gene:G9\n\
21\tx\tgene\t1\t2\t.\t+\t.\n\
21\tx\tgene\tabc\t2\t.\t+\t.\tID=gene:G9\n\
21\tx\tgene\t1\t-2\t.\t+\t.\tID=gene:G9\n";
        let sequence: String = std::iter::repeat_n('N', MAX_LINE_BYTES + 1).collect();
        let text = format!(
            "{}>20\nNNNNNNNN\n{unread}{sequence}\n",
            three_lines("20", "0")
        );
        assert_eq!(parse(&text).len(), 1);
        let err = try_parse(&format!("{text}{SOUND_UTR_ON_FIRST}")).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 15: GFF3 annotation after the sequence section, which opens at line 4"
            ),
            "{err:#}"
        );
        let text = format!("{}>20\n{sequence}\n", three_lines("20", "0"));
        assert_eq!(parse(&text).len(), 1);
        let err = try_parse(&format!("{text}{SOUND_UTR_ON_FIRST}")).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 6: GFF3 annotation after the sequence section, which opens at line 4"
            ),
            "{err:#}"
        );
        let held = |eighth_column: usize| {
            let filler: String = std::iter::repeat_n('x', eighth_column).collect();
            format!(
                "{}>20\nNNNNNNNN\n21\tx\tgene\t1\t2\t.\t+\t{filler}\tID=gene:G9\n",
                three_lines("20", "0")
            )
        };
        assert_eq!(parse(&held(MAX_LINE_BYTES - 18)).len(), 1);
        let err = try_parse(&held(MAX_LINE_BYTES - 19)).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 6: GFF3 annotation after the sequence section, which opens at line 4"
            ),
            "{err:#}"
        );
    }

    /// A line reads as a record when it has nine or more columns and its
    /// fourth and fifth are whole numbers, whatever opens it; a line of
    /// fewer columns does not, numbers in its fourth and fifth or not, nor
    /// one whose start or end is no whole number.
    #[test]
    fn a_line_reads_as_a_record_with_nine_columns_and_a_numeric_start_and_end() {
        for record in [
            ">20\tx\texon\t1000\t1100\t.\t+\t.\tParent=transcript:T0",
            "20\tx\texon\t1000\t1100\t.\t+\t.\tParent=transcript:T0\textra",
            "20\tx\texon\t0\t0\t.\t+\t.\t",
        ] {
            assert!(reads_as_a_record(record.as_bytes()), "{record:?}");
        }
        for other in [
            ">20 placeholder",
            ">20\tx\texon\t1000",
            ">20\tx\texon\t1000\t1100",
            ">20\tx\texon\t1000\t1100\t.\t+\t.",
            ">20\tx\texon\tabc\t1100\t.\t+\t.\tParent=transcript:T0",
            ">20\tx\texon\t1000\t-9\t.\t+\t.\tParent=transcript:T0",
            ">20\tx\texon\t1000x\t1100\t.\t+\t.\tParent=transcript:T0",
            ">20\tx\texon\t\t1100\t.\t+\t.\tParent=transcript:T0",
        ] {
            assert!(!reads_as_a_record(other.as_bytes()), "{other:?}");
        }
    }

    /// The index build places a line of four or more columns whose start a
    /// number leads and whose end, where the line has one, a number leads
    /// that does not read below zero, whatever opens the line and whatever
    /// bytes stand in its other columns; every line that reads as a record
    /// is one. It places no line of fewer than four columns, and none whose
    /// start or end it fails.
    #[test]
    fn the_index_build_places_a_line_whose_start_and_end_it_parses() {
        for placed in [
            &b">20\tx\texon\t1000"[..],
            b">20\tx\texon\t1000\t1100",
            b"20\tx\texon\t1000x\t1100\t.\t+\t.",
            b"20\tx\texon\t-5\t0x1A\t.\t+\t.\tParent=transcript:T0",
            b"20\tx\texon\t\t1100\t.\t+\t.\tParent=transcript:T0",
            b"20\tx\texon\t1000\t\t5",
            b"\t\t\t7",
            b"caf\xe9\tx\texon\t1000\t1100\t.\t+\t.\tParent=transcript:T0",
            b"20\tx\texon\t1000\t1100\t.\t+\t.\tParent=transcript:T0\r",
        ] {
            assert!(index_build_places(placed), "{placed:?}");
        }
        for failed in [
            &b""[..],
            b">20 placeholder",
            b"20\tx\texon",
            b"20\tx\texon\t",
            b"20\tx\texon\tabc\t1100",
            b"20\tx\texon\t1000\tabc",
            b"20\tx\texon\t1000\t-9\t.\t+\t.\tParent=transcript:T0",
            b"20\tx\texon\t\xe91000\t1100",
        ] {
            assert!(!index_build_places(failed), "{failed:?}");
        }
    }

    /// A line that opens with a tab is read at the columns the file writes,
    /// where the index build reads them: a stray tab before a record puts its
    /// type in the start column, so the line refuses, the record's type kept
    /// or skipped, on either sequence, and a record whose seqid alone is
    /// empty refuses on a kept type for a seqid that names no chromosome. A
    /// line led by spaces alone, no tab before its seqid, is read trimmed.
    #[test]
    fn a_line_that_opens_with_a_tab_is_read_at_its_columns_as_written() {
        let exon = line_of("20", "exon", "1500", "1600", "0");
        let err = try_parse(&two_sequences(&format!("\t{exon}"), SOUND_UTR_ON_LATER)).unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gff3 line 4: GFF3 start \"exon\" is not a number"),
            "{err:#}"
        );
        let later = format!("\t{}", line_of("21", "exon", "1500", "1600", "1"));
        let err = try_parse(&two_sequences(SOUND_UTR_ON_FIRST, &later)).unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gff3 line 8: GFF3 start \"exon\" is not a number"),
            "{err:#}"
        );
        for kept in [
            "\t20\tgene\texon\t1500\t1600\t.\t+\t.\tParent=transcript:T0\n".to_string(),
            line_of("", "exon", "1500", "1600", "0"),
        ] {
            let err = try_parse(&two_sequences(&kept, SOUND_UTR_ON_LATER)).unwrap_err();
            assert!(
                format!("{err:#}").ends_with("t.gff3 line 4: invalid chromosome name ''"),
                "{err:#}"
            );
        }
        let genes = try_parse(&two_sequences(&format!("  {exon}"), SOUND_UTR_ON_LATER)).unwrap();
        assert_eq!(transcript(&genes, "T0").exons.len(), 2);
    }

    /// White space and then a tab before a record are a column to the index
    /// build, which reads the record's type where its start is and fails the
    /// line: the parse refuses it on a kept type as on a skipped one, on
    /// either sequence, naming the type as the start.
    #[test]
    fn a_kept_line_led_by_white_space_and_a_tab_is_refused_as_its_skipped_twin_is() {
        for lead in [" \t", "   \t"] {
            for feature_type in ["exon", "five_prime_UTR"] {
                let message = format!("GFF3 start \"{feature_type}\" is not a number");
                let on_first = line_of("20", feature_type, "1500", "1600", "0");
                let text = two_sequences(&format!("{lead}{on_first}"), SOUND_UTR_ON_LATER);
                let err = try_parse(&text).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gff3 line 4: {message}")),
                    "{err:#}"
                );
                let on_later = line_of("21", feature_type, "1500", "1600", "1");
                let text = two_sequences(SOUND_UTR_ON_FIRST, &format!("{lead}{on_later}"));
                let err = try_parse(&text).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gff3 line 8: {message}")),
                    "{err:#}"
                );
            }
        }
    }

    /// A line led by spaces alone, no tab before its seqid, is read trimmed:
    /// a gene line's seqid is the one behind the spaces, a kept record joins
    /// its transcript and a skipped one is skipped.
    #[test]
    fn a_line_led_by_spaces_alone_is_read_trimmed() {
        let model = format!(
            "  {}",
            two_sequences(SOUND_UTR_ON_FIRST, SOUND_UTR_ON_LATER)
        );
        assert_eq!(parse(&model)["G0"].chr, "20");
        let exon = format!("  {}", line_of("20", "exon", "1500", "1600", "0"));
        let genes = parse(&two_sequences(&exon, SOUND_UTR_ON_LATER));
        assert_eq!(transcript(&genes, "T0").exons.len(), 2);
        let utr = format!("  {SOUND_UTR_ON_FIRST}");
        let genes = parse(&two_sequences(&utr, SOUND_UTR_ON_LATER));
        assert_eq!(transcript(&genes, "T0").exons.len(), 1);
    }

    /// A line the Ensembl-cache reading keeps by its `gene:` or `transcript:`
    /// identifier is a kept line and holds the whole rule, whatever its type:
    /// a start of 0, an end before its start and a start a number leads
    /// without being one refuse the parse there, where the `--gff` reading
    /// skips the line by its type and reads none of the three. A start no
    /// number leads refuses under both readings.
    #[test]
    fn a_line_kept_by_its_ensembl_prefix_holds_the_whole_coordinate_rule() {
        let tec = |start: &str, end: &str| {
            format!(
                "{}21\tx\tunconfirmed_transcript\t{start}\t{end}\t.\t+\t.\tID=transcript:T2;Parent=gene:G1;biotype=TEC\n",
                three_lines("21", "1")
            )
        };
        for (start, end, message) in [
            (
                "0",
                "2000",
                "GFF3 coordinates 0-2000 are not 1-based and ordered",
            ),
            (
                "2000",
                "1000",
                "GFF3 coordinates 2000-1000 are not 1-based and ordered",
            ),
            ("1000x", "2000", "GFF3 start \"1000x\" is not a number"),
        ] {
            let err = try_parse_as(&tec(start, end), Reading::EnsemblCache).unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gff3 line 4: {message}")),
                "{err:#}"
            );
            let genes = try_parse(&tec(start, end))
                .unwrap_or_else(|err| panic!("{start:?} {end:?}: {err:#}"));
            assert_eq!(genes["G1"].transcripts.len(), 1);
        }
        for reading in [Reading::EnsemblCache, Reading::GffSource] {
            let err = try_parse_as(&tec("abc", "2000"), reading).unwrap_err();
            assert!(
                format!("{err:#}").ends_with("t.gff3 line 4: GFF3 start \"abc\" is not a number"),
                "{err:#}"
            );
        }
        assert!(try_parse_as(&tec("1000", "2000"), Reading::EnsemblCache).is_ok());
    }

    /// The identifier read off a skipped line's column is the one its hash
    /// names, on each of these shapes: the keys tried in their own order
    /// whatever the column's, a later pair replacing an earlier one, a pair
    /// without `=` unsetting its key, a `0` and an empty value that are not
    /// set, a value set as written that reads `0` once decoded, a set value
    /// without a word character, which names nothing whatever a later key
    /// holds, an escaped value, a value cut at its second `=`, a key that
    /// opens with a blank, and no identifier at all.
    #[test]
    fn a_skipped_lines_identifier_is_the_one_its_hash_names() {
        for (column, id) in [
            ("ID=gene:G1;Name=X", Some("gene:G1")),
            ("Name=N1;ID=I1", Some("I1")),
            ("id=i1;Name=N1", Some("N1")),
            ("ID=0;Name=foo", Some("foo")),
            ("ID=%30;Name=n", Some("0")),
            ("ID=A;ID=B", Some("B")),
            ("ID=A;ID", None),
            ("ID=A;ID;Name=n", Some("n")),
            ("ID=-;Name=n", None),
            ("ID=%2C;Name=w%5F1", Some(",")),
            ("ID=;Name=w%5F1", Some("w_1")),
            ("id=low;name=lower", Some("low")),
            ("ID==;Name=x=y", Some("x")),
            ("biotype=protein_coding", None),
            ("", None),
            ("ID=;Name=0;id=0;name=last;", Some("last")),
            (" ID=spaced;name=n", Some("n")),
        ] {
            assert_eq!(skipped_id(column).as_deref(), id, "{column:?}");
            assert_eq!(
                skipped_id(column),
                named_id(&parse_attributes(column)),
                "{column:?}"
            );
        }
    }

    #[test]
    fn test_parse_attributes() {
        let attrs = parse_attributes("ID=gene:ENSG00000123;Name=BRCA2;biotype=protein_coding");
        assert_eq!(attrs.get("ID").unwrap(), "gene:ENSG00000123");
        assert_eq!(attrs.get("Name").unwrap(), "BRCA2");
        assert_eq!(attrs.get("biotype").unwrap(), "protein_coding");
    }

    /// The hash holds a value as written and [`text`] decodes it: every
    /// `%XX` escape, in either letter case, in one pass, so `%2520` is `%20`
    /// and `%3E`, an escape of a character GFF3 does not reserve, is `>`; a
    /// GTF value, decoded when its hash was built, is left as it stands. A
    /// value whose escapes spell bytes that are not UTF-8 is refused by
    /// [`check_attributes`], naming the key and the value, as the GTF reader
    /// names one; a line holding two such values names the one written first,
    /// `Name` before `ID` when the line writes them so, on every run, whatever
    /// order the hash iterates in, and a pair that a later pair with its key
    /// replaced is not the hash's value, so it names nothing.
    #[test]
    fn test_parse_attributes_url_decode() {
        let attrs = parse_attributes("description=Gene%20name%2C%20variant");
        assert_eq!(
            attrs.get("description").unwrap(),
            "Gene%20name%2C%20variant"
        );
        assert_eq!(
            text(Model::Gff3, attrs.get("description").unwrap()),
            "Gene name, variant"
        );
        assert_eq!(text(Model::Gff3, "A%2cB%3e%2520%30"), "A,B>%200");
        assert_eq!(text(Model::Gff3, "a%2Gb%2"), "a%2Gb%2");
        assert_eq!(text(Model::Gtf, "A%2cB%2520"), "A%2cB%2520");
        let check = |column: &str| check_attributes(column, &parse_attributes(column));
        assert_eq!(
            check("ID=g1;Name=caf%FF").unwrap_err().to_string(),
            "GFF3 attribute Name \"caf%FF\" decodes to text that is not UTF-8"
        );
        check("ID=g1;Name=caf%C3%A9;x=%2G").unwrap();
        // Each parse builds a hash with an iteration order of its own, so a
        // pass over several holds the key named to the one written first
        // whatever the hash's order is.
        for _ in 0..32 {
            assert_eq!(
                check("Name=caf%FF;ID=g%FF").unwrap_err().to_string(),
                "GFF3 attribute Name \"caf%FF\" decodes to text that is not UTF-8"
            );
            assert_eq!(
                check("Name=caf%FF;ID=g%FF;Name=ok")
                    .unwrap_err()
                    .to_string(),
                "GFF3 attribute ID \"g%FF\" decodes to text that is not UTF-8"
            );
        }
        check("Name=caf%FF;Name=ok;Name").unwrap();
        check("Name=caf%FF;Name").unwrap();
    }

    /// The pairs split as `GFF3.pm:194` splits them: a value runs from the
    /// first `=` to the second (`ID=a=b` is `a`, `Name=x=y=z` is `x`,
    /// `ID==x` is empty), a pair without `=` sets nothing and unsets the key
    /// an earlier pair set, a later pair with the key wins, nothing is
    /// trimmed (` biotype` is the key written), a closing `;` is no pair (so
    /// the empty key keeps its value) and an escape is decoded after the cut,
    /// so `%3D` stays an `=` of the value.
    #[test]
    fn attribute_values_are_cut_at_a_second_equals_sign_as_the_parser_cuts_them() {
        let attrs = parse_attributes(
            "ID=a=b;Name=x=y=z;tag=;flag; biotype=p;ID==x;gene_id=A;gene_id=B;desc=k%3Dv=w;logic=l;logic;=e%2C;",
        );
        let get = |key: &str| attrs.get(key).map(|value| text(Model::Gff3, value));
        assert_eq!(get("ID").as_deref(), Some(""));
        assert_eq!(get("Name").as_deref(), Some("x"));
        assert_eq!(get("tag").as_deref(), Some(""));
        assert_eq!(get("flag"), None);
        assert_eq!(get("biotype"), None);
        assert_eq!(get(" biotype").as_deref(), Some("p"));
        assert_eq!(get("gene_id").as_deref(), Some("B"));
        assert_eq!(get("desc").as_deref(), Some("k=v"));
        assert_eq!(get("logic"), None);
        assert_eq!(get("").as_deref(), Some("e,"));
        assert_eq!(parse_attributes("ID=a=b;Name=n").get("ID").unwrap(), "a");
        assert_eq!(parse_attributes("").len(), 0);
    }

    /// Through the model: a transcript whose line reads `biotype=lncRNA=x`
    /// has the biotype `lncRNA`, the field before the second `=`, as Ensembl
    /// VEP's `--gff` source reads it; a gene `Name=A%2CB` is the symbol
    /// `A,B`, which that source prints as written.
    #[test]
    fn a_value_with_a_second_equals_sign_is_cut_before_it_reaches_the_transcript() {
        let genes = parse(
            "21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;Name=A%2CB;biotype=protein_coding\n\
21\tx\ttranscript\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA=x\n\
21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T1;Name=X1\n",
        );
        let tr = transcript(&genes, "T1");
        assert_eq!(tr.biotype, "lncRNA");
        assert_eq!(tr.symbol.as_deref(), Some("A,B"));
    }

    /// Through the model, every escape of a read value is decoded, in either
    /// letter case and whether or not GFF3 reserves the character, and each
    /// test Ensembl VEP's `--gff` source runs on its hash is taken on the
    /// value as written: a gene `Name=A%2cB%3E` is the symbol `A,B>`; a
    /// transcript's `Parent=gene%3AG1` names the gene line `ID=gene%3AG1`, one
    /// parent whose escaped `:` is decoded on both sides; its `tag=%30` is set
    /// (`%30` is true in Perl) and reads `0`, so no membership follows; its
    /// `gene_id=%30` is set and is the cache's gene identifier `0` where a
    /// `gene_id=0` would fall to the gene line's. The Ensembl `gene:` and
    /// `transcript:` prefixes are read off the decoded identifier, so under
    /// [`Reading::EnsemblCache`] an `unconfirmed_transcript` line
    /// `ID=transcript%3AT2` is admitted and its `tag=basic` is `gencode_basic`.
    /// A line whose attribute value decodes to bytes that are not UTF-8 refuses
    /// the parse, naming the file and the line.
    #[test]
    fn every_escape_of_a_read_value_is_decoded_after_the_tests_on_the_written_value() {
        let lines = "21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene%3AG1;Name=A%2cB%3E;biotype=protein_coding\n\
21\tx\ttranscript\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene%3AG1;biotype=lncRNA;tag=%30;gene_id=%30\n\
21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T1;Name=X%2d1\n\
21\tx\tunconfirmed_transcript\t1000\t2000\t.\t+\t.\tID=transcript%3AT2;Parent=gene%3AG1;biotype=TEC;tag=basic\n\
21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript%3AT2;Name=X2\n";
        let genes = parse(lines);
        let tr = transcript(&genes, "T1");
        assert_eq!(tr.symbol.as_deref(), Some("A,B>"));
        assert_eq!(tr.gene_id, "0");
        assert!(!tr.gencode_basic && !tr.is_canonical);
        assert_eq!(tr.exons[0].exon_id.as_deref(), Some("X-1"));
        let gene = &genes["G1"];
        assert_eq!(gene.name.as_deref(), Some("A,B>"));
        assert_eq!(gene.transcripts.len(), 1, "the unlisted type is skipped");
        let genes = try_parse_as(lines, Reading::EnsemblCache).unwrap();
        let tec = transcript(&genes, "T2");
        assert_eq!((tec.biotype.as_str(), tec.gencode_basic), ("TEC", true));

        let err = try_parse(
            "21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;Name=caf%FF;biotype=protein_coding\n",
        )
        .unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 1: GFF3 attribute Name \"caf%FF\" decodes to text that is not UTF-8"
            ),
            "{err:#}"
        );
    }

    /// The type and source columns are read decoded, as `get_type` and
    /// `get_source` read them for the `--gff` source (`GFF3.pm:153,161-164`
    /// through `BaseGXF.pm:224,287-288`): a `g%65ne` line is a gene, a
    /// `transcri%70t` line its transcript and `havana%5ftagene` its source,
    /// where the written text matches no listed type. A type whose escapes
    /// spell bytes that are not UTF-8 matches no listed type and its line is
    /// skipped, as the line is on Ensembl VEP. A skipped line's source is never
    /// read and its attribute values are never checked, so the `region` line's
    /// `caf%FF` refuses nothing, where a gene or transcript line whose source
    /// does so refuses the parse, naming the file and the line.
    #[test]
    fn type_and_source_columns_are_read_decoded() {
        let genes = parse(
            "21\thavana%5ftagene\tg%65ne\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
21\thavana%5ftagene\ttranscri%70t\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
21\thavana%5ftagene\t%65xon\t1000\t2000\t.\t+\t.\tParent=transcript:T1;Name=X1\n\
21\tx\tgene%FF\t1000\t2000\t.\t+\t.\tID=skipped\n\
21\tcaf%FF\tregion\t1\t5000\t.\t.\t.\tID=21;Note=caf%FF\n",
        );
        assert_eq!(
            genes.len(),
            1,
            "the line of a type that is no UTF-8 is skipped"
        );
        let tr = transcript(&genes, "T1");
        assert_eq!(
            (tr.biotype.as_str(), tr.source.as_str()),
            ("protein_coding", "havana_tagene")
        );
        assert_eq!(tr.exons.len(), 1);

        let err = try_parse("21\tcaf%FF\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1\n").unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 1: GFF3 source \"caf%FF\" decodes to text that is not UTF-8"
            ),
            "{err:#}"
        );
    }

    /// `ccdsid` is the first CCDS xref's display label (`Transcript.pm`
    /// `summary_as_hash` through `ccds`), the same entry Ensembl VEP stores
    /// as `_ccds` (`Database/Transcript.pm` `prefetch_transcript_ids`), so
    /// `CCDS74772.1` on ENST00000622113 is the cache's `ccds`.
    #[test]
    fn ccdsid_becomes_ccds() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}"));
        let tr = transcript(&genes, "ENST00000622113");
        assert_eq!(tr.ccds.as_deref(), Some("CCDS74772.1"));
        let no_ccds = TPTE_205.replace("ccdsid=CCDS74772.1;", "");
        let genes = parse(&format!("{TPTE_GENE}{no_ccds}{TPTE_205_EXON}"));
        assert_eq!(transcript(&genes, "ENST00000622113").ccds, None);
    }

    /// The gene line's `version` (`Feature.pm` `summary_as_hash`) is the
    /// `_gene_version` Ensembl VEP appends to `Gene` under `--gene_version`:
    /// `6` for ENSG00000274391.
    #[test]
    fn gene_version_from_gff3() {
        let genes = parse(TPTE_GENE);
        let gene = &genes["ENSG00000274391"];
        assert_eq!(gene.version, Some(6));
        assert_eq!(gene.name.as_deref(), Some("TPTE"));
        assert_eq!(gene.biotype, "protein_coding");
        let unversioned = TPTE_GENE.replace(";version=6", "");
        assert_eq!(parse(&unversioned)["ENSG00000274391"].version, None);
    }

    /// A `start` or `end` column that is not a number fails the whole parse
    /// with an error naming the file, the line, the column and the value, so
    /// a malformed line is never binned at position 0.
    #[test]
    fn non_numeric_coordinate_fails_the_parse() {
        for (lines, message) in [
            (
                TPTE_GENE.replace("\t10521553\t", "\tabc\t"),
                "t.gff3 line 1: GFF3 start \"abc\" is not a number",
            ),
            (
                TPTE_GENE.replace("\t10606140\t", "\t1060614x\t"),
                "t.gff3 line 1: GFF3 end \"1060614x\" is not a number",
            ),
        ] {
            let err = try_parse(&lines).unwrap_err();
            assert!(format!("{err:#}").ends_with(message), "{err:#}");
        }
    }

    /// A start of 0 or an end before its start fails the parse, naming the
    /// file and the line: the region of a start is binned from `start - 1`,
    /// which a 0 has no room for, and a span whose end precedes its start has
    /// no length.
    #[test]
    fn a_start_of_zero_or_an_end_before_its_start_fails_the_parse() {
        for (lines, message) in [
            (
                TPTE_GENE.replace("\t10521553\t", "\t0\t"),
                "t.gff3 line 1: GFF3 coordinates 0-10606140 are not 1-based and ordered",
            ),
            (
                TPTE_GENE.replace("\t10606140\t", "\t10521552\t"),
                "t.gff3 line 1: GFF3 coordinates 10521553-10521552 are not 1-based and ordered",
            ),
        ] {
            let err = try_parse(&lines).unwrap_err();
            assert!(format!("{err:#}").ends_with(message), "{err:#}");
        }
        let genes = parse(&TPTE_GENE.replace("\t10606140\t", "\t10521553\t"));
        assert_eq!(genes["ENSG00000274391"].end, 10521553);
    }

    /// A seqid that would place the chromosome's directory outside the cache (a
    /// parent-directory walk, a rooted path) fails the parse, naming the file,
    /// the line, counted from the file's first line with the header included,
    /// and the seqid as written. A seqid of the three letters `chr` alone is
    /// a name, kept as written.
    #[test]
    fn seqid_that_is_not_one_path_component_fails_the_parse() {
        for (seqid, message) in [
            (
                "../../escaped",
                "t.gff3 line 2: invalid chromosome name '../../escaped'",
            ),
            (
                "/tmp/escaped",
                "t.gff3 line 2: invalid chromosome name '/tmp/escaped'",
            ),
            ("..", "t.gff3 line 2: invalid chromosome name '..'"),
        ] {
            let lines = format!("##gff-version 3\n{}", TPTE_GENE.replacen("21", seqid, 1));
            let err = try_parse(&lines).unwrap_err();
            assert!(format!("{err:#}").ends_with(message), "{err:#}");
        }
        let genes = parse(&format!(
            "##gff-version 3\n{}",
            TPTE_GENE.replacen("21", "HSCHR6_MHC_COX", 1)
        ));
        assert_eq!(genes["ENSG00000274391"].chr, "HSCHR6_MHC_COX");
        let genes = parse(&format!(
            "##gff-version 3\n{}",
            TPTE_GENE.replacen("21", "chr", 1)
        ));
        assert_eq!(genes["ENSG00000274391"].chr, "chr");
    }

    /// A line that is not UTF-8 (a Latin-1 `é`, 0xE9, in the gene's `Name`)
    /// fails the parse, naming the file and the line, as every refusal of a
    /// line does.
    #[test]
    fn a_line_that_is_not_utf8_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gff3");
        let (before_name, after_name) = TPTE_GENE.split_once("TPTE").unwrap();
        let mut bytes = b"##gff-version 3\n".to_vec();
        bytes.extend_from_slice(before_name.as_bytes());
        bytes.extend_from_slice(b"TPT\xe9");
        bytes.extend_from_slice(after_name.as_bytes());
        std::fs::write(&path, bytes).unwrap();
        let err = parse_gff3(&path, Reading::GffSource, &SequenceNames::default()).unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gff3 line 2: GFF3 line is not UTF-8"),
            "{err:#}"
        );
    }

    /// Exons with a rank sort by it; once any exon of the transcript lacks one,
    /// every exon sorts by its start, as a rank compared with a start orders
    /// nothing, and the file is parsed whole.
    #[test]
    fn exons_mixing_ranked_and_unranked_sort_by_start() {
        let exon = |start: u64, rank: Option<u32>| {
            let rank = rank.map(|r| format!(";rank={r}")).unwrap_or_default();
            format!(
                "21\tensembl_havana\texon\t{start}\t{}\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE{start}{rank}\n",
                start + 99
            )
        };
        let exon_starts = |genes: &HashMap<String, GeneRecord>| -> Vec<u64> {
            transcript(genes, "ENST00000622113")
                .exons
                .iter()
                .map(|e| e.start)
                .collect()
        };
        // Sixty-four exons, every other one ranked, listed in a scrambled order.
        let mut lines = format!("{TPTE_GENE}{TPTE_205}");
        for i in (0..64u64).map(|i| (i * 37) % 64) {
            let rank = (i % 2 == 0).then_some(64 - i as u32);
            lines.push_str(&exon(10_521_600 + i * 1_000, rank));
        }
        let expected: Vec<u64> = (0..64).map(|i| 10_521_600 + i * 1_000).collect();
        assert_eq!(exon_starts(&parse(&lines)), expected);

        let ranked = format!(
            "{TPTE_GENE}{TPTE_205}{}{}",
            exon(10_540_000, Some(2)),
            exon(10_600_000, Some(1))
        );
        assert_eq!(exon_starts(&parse(&ranked)), [10_600_000, 10_540_000]);
    }

    /// A bgzipped GFF3 is a run of gzip members, the shape Ensembl VEP's
    /// `--gff` takes: every member is read, so the gene tree is the plain
    /// file's, not the first block's alone.
    #[test]
    fn bgzipped_gff3_is_read_past_its_first_member() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gff3.gz");
        let mut bytes = Vec::new();
        for member in [
            "##gff-version 3\n".to_string(),
            format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}"),
        ] {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(member.as_bytes()).unwrap();
            bytes.extend(encoder.finish().unwrap());
        }
        std::fs::write(&path, bytes).unwrap();
        let genes = parse_gff3(&path, Reading::GffSource, &SequenceNames::default()).unwrap();
        assert_eq!(genes["ENSG00000274391"].transcripts.len(), 1);
        assert_eq!(transcript(&genes, "ENST00000622113").exons.len(), 1);
    }

    /// A line past the cap (the gene line with an attribute value of the cap's
    /// length) stops the parse naming the file and its line, 2.
    #[test]
    fn line_past_the_cap_fails_the_parse() {
        let padding: String = std::iter::repeat_n('x', MAX_LINE_BYTES).collect();
        let lines = format!(
            "{}{}",
            TPTE_205,
            TPTE_GENE.replace("Name=TPTE;", &format!("Name=TPTE;Note={padding};"))
        );
        let err = try_parse(&lines).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "t.gff3 line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
    }

    /// The GENCODE basic tag is `gencode_basic` on the GRCh38 file (the
    /// attribute code, `Transcript.pm` `summary_as_hash`) and `basic` on the
    /// GRCh37 release-87 file; either marks the transcript.
    #[test]
    fn gencode_basic_tag_both_spellings() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}"));
        assert!(transcript(&genes, "ENST00000622113").gencode_basic);

        let genes = parse(GRCH37_PSEUDOGENE);
        let tr = transcript(&genes, "ENST00000411878");
        assert!(tr.gencode_basic);
        assert!(!tr.gencode_primary);

        let untagged = TPTE_205.replace("tag=gencode_basic,gencode_primary;", "");
        let genes = parse(&format!("{TPTE_GENE}{untagged}{TPTE_205_EXON}"));
        assert!(!transcript(&genes, "ENST00000622113").gencode_basic);
    }

    /// `gencode_primary` in `tag` is the GENCODE primary attribute's code
    /// (`Transcript.pm` `summary_as_hash`).
    #[test]
    fn gencode_primary_tag_from_gff3() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}"));
        assert!(transcript(&genes, "ENST00000622113").gencode_primary);
        let basic_only =
            TPTE_205.replace("tag=gencode_basic,gencode_primary;", "tag=gencode_basic;");
        let genes = parse(&format!("{TPTE_GENE}{basic_only}{TPTE_205_EXON}"));
        let tr = transcript(&genes, "ENST00000622113");
        assert!(tr.gencode_basic);
        assert!(!tr.gencode_primary);
    }

    /// `cds_start_NF` in a transcript line's `tag` list is carried on the
    /// record; a line without it carries none.
    #[test]
    fn the_start_incomplete_tag_is_carried_from_the_transcript_line() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}"));
        assert!(!transcript(&genes, "ENST00000622113").cds_start_nf);
        let flagged = TPTE_205.replace(
            "tag=gencode_basic,gencode_primary;",
            "tag=gencode_basic,cds_start_NF;",
        );
        assert_ne!(flagged, TPTE_205);
        let genes = parse(&format!("{TPTE_GENE}{flagged}{TPTE_205_EXON}"));
        let tr = transcript(&genes, "ENST00000622113");
        assert!(tr.cds_start_nf && tr.gencode_basic);
    }

    /// A `tag` list is split on `,` alone and each member compared whole
    /// (`BaseGXF.pm` lines 683-685, then 690-691 `/^$attr$/` and 700
    /// `/^Ensembl_canonical$/`), so a member padded with a blank,
    /// ` Ensembl_canonical` or `MANE_Select `, names no tag, while the
    /// member beside it reads as in the unpadded list.
    #[test]
    fn a_tag_list_member_padded_with_a_blank_names_no_tag() {
        let memberships = |list: &str| {
            let line = TPTE_205.replace(
                "tag=gencode_basic,gencode_primary;",
                &format!("tag={list};"),
            );
            let genes = parse(&format!("{TPTE_GENE}{line}{TPTE_205_EXON}"));
            let tr = transcript(&genes, "ENST00000622113");
            (tr.gencode_basic, tr.is_canonical, tr.mane_select)
        };
        assert_eq!(
            memberships("gencode_basic, Ensembl_canonical,MANE_Select "),
            (true, false, false),
            "a padded member is not the tag"
        );
        assert_eq!(
            memberships("gencode_basic,Ensembl_canonical,MANE_Select"),
            (true, true, true)
        );
    }

    /// A `tag` value is split only when it is true to Perl (`BaseGXF.pm` line
    /// 683 `if($tr_record->{attributes}->{tag})`), a test on the value as
    /// written, so `tag=0` and `tag=` list no tag where the list
    /// `0,gencode_basic` opens with the member `0` and `tag=%30`, the set
    /// value `%30`, lists the one tag that reads `0` once decoded. The split is
    /// on the commas written, so `gencode%5Fbasic,a%2Cb` lists `gencode_basic`
    /// and `a,b`.
    #[test]
    fn a_tag_value_that_is_not_true_to_perl_lists_no_tag() {
        let tags = |column: &str| tag_list(&parse_attributes(column));
        let none: [&str; 0] = [];
        assert_eq!(tags("ID=transcript:t1;tag=0"), none, "`0` is not true");
        assert_eq!(
            tags("ID=transcript:t1;tag="),
            none,
            "the empty value lists nothing"
        );
        assert_eq!(tags("ID=transcript:t1;tag=%30"), ["0"]);
        assert_eq!(
            tags("ID=transcript:t1;tag=gencode%5Fbasic,a%2Cb"),
            ["gencode_basic", "a,b"]
        );
        assert_eq!(
            tags("ID=transcript:t1;tag=0,gencode_basic"),
            ["0", "gencode_basic"]
        );
        assert_eq!(
            tags("ID=transcript:t1;tag=gencode_basic,Ensembl_canonical"),
            ["gencode_basic", "Ensembl_canonical"]
        );
    }

    /// `transcript_support_level` is the `TSL` attribute value without its
    /// `tsl` prefix (`Transcript.pm` `tsl`): `1` here, and the `assigned to
    /// previous version` suffix is kept verbatim, as it is in the attribute
    /// Ensembl VEP prints the digit from (`OutputFactory.pm` `m/tsl(\d+)/`).
    #[test]
    fn tsl_from_transcript_support_level() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}"));
        assert_eq!(
            transcript(&genes, "ENST00000622113").tsl.as_deref(),
            Some("1")
        );

        let col6a2 = "21\tensembl_havana\tgene\t46096722\t46132851\t.\t+\t.\tID=gene:ENSG00000142173;Name=COL6A2;biotype=protein_coding;description=collagen type VI alpha 2 chain [Source:HGNC Symbol%3BAcc:HGNC:2212];gene_id=ENSG00000142173;logic_name=ensembl_havana_gene_homo_sapiens;version=19\n\
21\thavana\tmRNA\t46098112\t46129801\t.\t+\t.\tID=transcript:ENST00000397763;Parent=gene:ENSG00000142173;Name=COL6A2-202;biotype=protein_coding;ccdsid=CCDS13729.1;tag=gencode_basic,gencode_primary,MANE_Plus_Clinical;transcript_id=ENST00000397763;transcript_support_level=5 (assigned to previous version 5);version=6\n\
21\thavana\texon\t46098112\t46098300\t.\t+\t.\tParent=transcript:ENST00000397763;Name=ENSE00001512888;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001512888;rank=1;version=1\n";
        let genes = parse(col6a2);
        let tr = transcript(&genes, "ENST00000397763");
        assert_eq!(
            tr.tsl.as_deref(),
            Some("5 (assigned to previous version 5)")
        );
        assert!(tr.mane_plus_clinical);
        assert!(!tr.mane_select);
        assert!(!tr.is_canonical);

        let no_tsl = TPTE_205.replace("transcript_support_level=1;", "");
        let genes = parse(&format!("{TPTE_GENE}{no_tsl}{TPTE_205_EXON}"));
        assert_eq!(transcript(&genes, "ENST00000622113").tsl, None);
    }

    /// A CDS line carries `protein_id`, the translation's stable identifier,
    /// and `version`, the translation's version (`CDS.pm` `summary_as_hash`):
    /// ENSP00000482040 version 1 for ENST00000622113, the identifier Ensembl
    /// VEP stores as `_protein` (`Database/Transcript.pm`
    /// `prefetch_translation_ids`). The GRCh37 release-87 file has no
    /// `version` on its CDS lines: Ensembl VEP's `--gff` source gives such a
    /// translation the version 1 (`BaseGXF.pm` `_add_translation`, `version
    /// || 1`), and a build of the cache's content takes the core table dumps'
    /// version when they hold the translation, else the protein FASTA
    /// header's.
    #[test]
    fn protein_id_and_version_from_cds_line() {
        let genes = parse(&format!(
            "{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}{TPTE_205_CDS}"
        ));
        let tr = transcript(&genes, "ENST00000622113");
        assert_eq!(tr.cds_regions.len(), 2);
        for cds in &tr.cds_regions {
            assert_eq!(cds.protein_id.as_deref(), Some("ENSP00000482040"));
            assert_eq!(cds.protein_version, Some(1));
        }
        assert_eq!(tr.cds_regions[0].phase, 0);
        assert_eq!(tr.cds_regions[1].phase, 1);

        let grch37_cds = TPTE_205_CDS.replace(";version=1", "");
        let lines = format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}{grch37_cds}");
        let genes = parse(&lines);
        let cds = &transcript(&genes, "ENST00000622113").cds_regions[0];
        assert_eq!(cds.protein_id.as_deref(), Some("ENSP00000482040"));
        assert_eq!(cds.protein_version, Some(1));
        let genes = try_parse_as(&lines, Reading::EnsemblCache).unwrap();
        let cds = &transcript(&genes, "ENST00000622113").cds_regions[0];
        assert_eq!(cds.protein_id.as_deref(), Some("ENSP00000482040"));
        assert_eq!(cds.protein_version, None);
    }

    /// A phase column that reads as the number 1 or 2 is that number, blanks
    /// before the digit or text behind it or not, as the comparison of
    /// numbers reads it (`BaseGXF.pm:815-820`); text holding a `.` states
    /// none; `0` and `-1` are themselves; any other text is refused, naming
    /// the format and the value.
    #[test]
    fn a_phase_that_reads_as_one_or_two_is_that_number() {
        for (value, phase) in [
            ("0", 0),
            ("1", 1),
            ("2", 2),
            ("-1", -1),
            (".", -1),
            ("0.5", -1),
            ("1.0", -1),
            ("1x", 1),
            ("2x", 2),
            (" 1", 1),
            ("2 ", 2),
            ("+1", 1),
            ("01", 1),
            ("1e0", 1),
            ("1\0", 1),
        ] {
            assert_eq!(parse_phase(Model::Gff3, value).unwrap(), phase, "{value:?}");
        }
        for value in ["", "x", "7", "3", "00", "0\0", "-1x", "0x1"] {
            assert_eq!(
                parse_phase(Model::Gtf, value).unwrap_err().to_string(),
                format!("GTF phase {value:?} is not a phase")
            );
        }
    }

    /// A `version` whose digits do not fit the cache's 32-bit field fails the
    /// parse, naming the file, the line and the limit, on a gene line, a
    /// transcript line (its `transcript_version` too) and a CDS line alike,
    /// as the line is read, so a transcript line the tree would drop for
    /// having no exon or CDS line is refused as well: a version read as none
    /// would lose what the file states. 4294967295 is held, and a gene line's
    /// `transcript_version`, which no gene reads, builds.
    #[test]
    fn a_version_above_the_cache_limit_fails_the_parse() {
        let above = ";version=4294967296";
        let limit = "version 4294967296 is above 4294967295, the largest version the cache holds";
        let gene = TPTE_GENE.replace(";version=6", above);
        let transcript_line = TPTE_205.replace(";version=4", above);
        let fallback = TPTE_205.replace(";version=4", ";transcript_version=4294967296");
        let cds = TPTE_205_CDS.replacen(";version=1", above, 1);
        for (lines, line) in [
            (format!("{gene}{TPTE_205}{TPTE_205_EXON}{TPTE_205_CDS}"), 1),
            (
                format!("{TPTE_GENE}{transcript_line}{TPTE_205_EXON}{TPTE_205_CDS}"),
                2,
            ),
            (
                format!("{TPTE_GENE}{fallback}{TPTE_205_EXON}{TPTE_205_CDS}"),
                2,
            ),
            (format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}{cds}"), 4),
            (format!("{TPTE_GENE}{transcript_line}"), 2),
        ] {
            let err = try_parse(&lines).unwrap_err();
            let message = format!("{err:#}");
            assert!(
                message.ends_with(&format!("t.gff3 line {line}: {limit}")),
                "{message}"
            );
        }
        let held = TPTE_GENE.replace(";version=6", ";version=4294967295");
        let genes = parse(&format!("{held}{TPTE_205}{TPTE_205_EXON}{TPTE_205_CDS}"));
        assert_eq!(genes["ENSG00000274391"].version, Some(4294967295));
        let unread = TPTE_GENE.replace(";version=6", ";transcript_version=4294967296");
        let genes = parse(&format!("{unread}{TPTE_205}{TPTE_205_EXON}{TPTE_205_CDS}"));
        assert_eq!(genes["ENSG00000274391"].version, None);
    }

    /// A `version` whose text is not a whole number fails the parse, naming
    /// the file, the line and the text, on a gene line, a transcript line
    /// (its `transcript_version` too) and a CDS line alike, as the line is
    /// read, where Ensembl VEP prints the text after the identifier: a
    /// version read as none would lose what the file states. A text whose
    /// digits would overflow but carries a non-digit is no whole number, not
    /// a number above the limit. A whole number written with leading zeros
    /// is that number, so `01` is 1.
    #[test]
    fn a_version_that_is_not_a_whole_number_fails_the_parse() {
        let gene = TPTE_GENE.replace(";version=6", ";version=v2");
        let transcript_line = TPTE_205.replace(";version=4", ";version=v2");
        let fallback = TPTE_205.replace(";version=4", ";transcript_version=v2");
        let cds = TPTE_205_CDS.replacen(";version=1", ";version=4294967296.0", 1);
        for (lines, line, text) in [
            (
                format!("{gene}{TPTE_205}{TPTE_205_EXON}{TPTE_205_CDS}"),
                1,
                "v2",
            ),
            (
                format!("{TPTE_GENE}{transcript_line}{TPTE_205_EXON}{TPTE_205_CDS}"),
                2,
                "v2",
            ),
            (
                format!("{TPTE_GENE}{fallback}{TPTE_205_EXON}{TPTE_205_CDS}"),
                2,
                "v2",
            ),
            (
                format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}{cds}"),
                4,
                "4294967296.0",
            ),
        ] {
            let err = try_parse(&lines).unwrap_err();
            let message = format!("{err:#}");
            assert!(
                message.ends_with(&format!(
                    "t.gff3 line {line}: version {text:?} is not a whole number"
                )),
                "{message}"
            );
        }
        let zeros = TPTE_205.replace(";version=4", ";version=01");
        let genes = parse(&format!("{TPTE_GENE}{zeros}{TPTE_205_EXON}{TPTE_205_CDS}"));
        assert_eq!(transcript(&genes, "ENST00000622113").version, Some(1));
    }

    /// A GFF3 `version` is parsed from its decoded value ([`text`]), on a
    /// gene line, a transcript line (its `transcript_version` too) and a CDS
    /// line alike: `%31` reads 1, `%30`, set where the written `0` is not,
    /// reads 0, and escapes that spell no whole number fail the parse naming
    /// the decoded text. A GTF value is parsed as the hash holds it, so an
    /// escape left in one is no digit.
    #[test]
    fn an_escaped_version_is_parsed_from_its_decoded_value() {
        let gene = TPTE_GENE.replace(";version=6", ";version=%31");
        let transcript_line = TPTE_205.replace(";version=4", ";version=%31%32");
        let fallback = TPTE_205.replace(";version=4", ";transcript_version=%33");
        let cds = TPTE_205_CDS.replace(";version=1", ";version=%37");
        let genes = parse(&format!("{gene}{transcript_line}{TPTE_205_EXON}{cds}"));
        assert_eq!(genes["ENSG00000274391"].version, Some(1));
        let tr = transcript(&genes, "ENST00000622113");
        assert_eq!(tr.version, Some(12));
        assert_eq!(tr.cds_regions.len(), 2);
        for segment in &tr.cds_regions {
            assert_eq!(segment.protein_version, Some(7));
        }
        let genes = parse(&format!(
            "{TPTE_GENE}{fallback}{TPTE_205_EXON}{TPTE_205_CDS}"
        ));
        assert_eq!(transcript(&genes, "ENST00000622113").version, Some(3));

        let zero = TPTE_GENE.replace(";version=6", ";version=%30");
        assert_eq!(parse(&zero)["ENSG00000274391"].version, Some(0));
        let unset = TPTE_GENE.replace(";version=6", ";version=0");
        assert_eq!(parse(&unset)["ENSG00000274391"].version, None);

        let lettered = TPTE_GENE.replace(";version=6", ";version=%76%32");
        let message = format!("{:#}", try_parse(&lettered).unwrap_err());
        assert!(
            message.ends_with("t.gff3 line 1: version \"v2\" is not a whole number"),
            "{message}"
        );

        let attrs = parse_attributes("version=%31");
        assert_eq!(
            parse_version(Model::Gff3, &attrs, &["version"]).unwrap(),
            Some(1)
        );
        assert_eq!(
            parse_version(Model::Gtf, &attrs, &["version"])
                .unwrap_err()
                .to_string(),
            "version \"%31\" is not a whole number"
        );
    }

    /// A transcript of biotype `artifact` and a readthrough transcript are
    /// transcripts of the tree, with the biotype and the tags the line
    /// carries, as Ensembl VEP's `--gff` source builds them (`BaseGXF.pm`
    /// `_create_transcript` drops a record without children, without a
    /// biotype, or `protein_coding` without an exon, and nothing else). The
    /// GRCh38 release-116 file spells the readthrough tag
    /// `readthrough_transcript` and the core API's attribute code is
    /// `readthrough_tra`; either marks the record `readthrough`.
    #[test]
    fn readthrough_and_artifact_transcripts_are_kept() {
        let lines = "21\tensembl_havana\tgene\t5022531\t5036771\t.\t+\t.\tID=gene:ENSG00000277117;biotype=artifact;description=novel protein%2C similar to inducible T-cell co-stimulator ligand ICOSLG;gene_id=ENSG00000277117;logic_name=ensembl_havana_gene_homo_sapiens;version=6\n\
21\thavana\ttranscript\t5022531\t5036771\t.\t+\t.\tID=transcript:ENST00000623960;Parent=gene:ENSG00000277117;biotype=artifact;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000623960;transcript_support_level=1 (assigned to previous version 3);version=4\n\
21\thavana\texon\t5022531\t5022693\t.\t+\t.\tParent=transcript:ENST00000623960;Name=ENSE00003755585;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003755585;rank=1;version=1\n\
21\thavana\tncRNA_gene\t10482738\t10605716\t.\t+\t.\tID=gene:ENSG00000273840;biotype=lncRNA;description=novel transcript;gene_id=ENSG00000273840;logic_name=havana_homo_sapiens;version=2\n\
21\thavana\tlnc_RNA\t10482738\t10605716\t.\t+\t.\tID=transcript:ENST00000612267;Parent=gene:ENSG00000273840;biotype=lncRNA;tag=gencode_basic,gencode_primary,Ensembl_canonical,readthrough_transcript;transcript_id=ENST00000612267;transcript_support_level=5;version=1\n\
21\thavana\texon\t10482738\t10482900\t.\t+\t.\tParent=transcript:ENST00000612267;Name=ENSE00003755586;rank=1\n\
21\thavana\tlnc_RNA\t10482738\t10500000\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000273840;biotype=lncRNA;tag=gencode_basic,readthrough_tra;transcript_id=ENST00000000001;version=1\n\
21\thavana\texon\t10482738\t10482900\t.\t+\t.\tParent=transcript:ENST00000000001;Name=ENSE00003755586;rank=1\n\
21\thavana\tlnc_RNA\t10482738\t10490000\t.\t+\t.\tID=transcript:ENST00000000002;Parent=gene:ENSG00000273840;biotype=lncRNA;tag=gencode_basic;transcript_id=ENST00000000002;version=1\n\
21\thavana\texon\t10482738\t10482900\t.\t+\t.\tParent=transcript:ENST00000000002;Name=ENSE00003755586;rank=1\n";
        let genes = parse(lines);
        let artifact = transcript(&genes, "ENST00000623960");
        assert_eq!(artifact.biotype, "artifact");
        assert!(artifact.is_canonical && artifact.gencode_basic && artifact.gencode_primary);
        assert!(!artifact.readthrough);
        let kept: Vec<(&str, bool)> = genes["ENSG00000273840"]
            .transcripts
            .iter()
            .map(|t| (t.transcript_id.as_str(), t.readthrough))
            .collect();
        assert_eq!(
            kept,
            vec![
                ("ENST00000612267", true),
                ("ENST00000000001", true),
                ("ENST00000000002", false),
            ]
        );
        let readthrough = transcript(&genes, "ENST00000612267");
        assert_eq!(readthrough.biotype, "lncRNA");
        assert!(readthrough.is_canonical && readthrough.gencode_primary);
        assert_eq!(readthrough.tsl.as_deref(), Some("5"));
    }

    /// Read for the Ensembl cache's content, a gene is any line whose `ID`
    /// is `gene:<id>`, so a lncRNA gene's `ncRNA_gene` line and a
    /// `pseudogene` gene line are genes; a transcript is any
    /// `transcript:<id>` line, so a `pseudogenic_transcript` line is a
    /// transcript, as is the GRCh37 release-87 file's `pseudogene` transcript
    /// line.
    #[test]
    fn gene_and_transcript_lines_by_id_prefix() {
        let lines = "21\thavana\tncRNA_gene\t5232668\t5243903\t.\t-\t.\tID=gene:ENSG00000279669;biotype=lncRNA;description=novel transcript;gene_id=ENSG00000279669;logic_name=havana_homo_sapiens;version=2\n\
21\thavana_tagene\tlnc_RNA\t5232668\t5243903\t.\t-\t.\tID=transcript:ENST00000623753;Parent=gene:ENSG00000279669;biotype=lncRNA;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000623753;transcript_support_level=3 (assigned to previous version 1);version=2\n\
21\thavana_tagene\texon\t5243800\t5243903\t.\t-\t.\tParent=transcript:ENST00000623753;Name=ENSE00003755587;rank=1\n\
21\tensembl_havana\tpseudogene\t9068414\t9073609\t.\t-\t.\tID=gene:ENSG00000188681;Name=TEKT4P2;biotype=transcribed_unprocessed_pseudogene;description=tektin 4 pseudogene 2 [Source:HGNC Symbol%3BAcc:HGNC:40046];gene_id=ENSG00000188681;logic_name=ensembl_havana_gene_homo_sapiens;version=14\n\
21\thavana\tpseudogenic_transcript\t9068414\t9073609\t.\t-\t.\tID=transcript:ENST00000559466;Parent=gene:ENSG00000188681;Name=TEKT4P2-201;biotype=transcribed_unprocessed_pseudogene;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000559466;transcript_support_level=NA;version=1\n\
21\thavana\texon\t9073500\t9073609\t.\t-\t.\tParent=transcript:ENST00000559466;Name=ENSE00003755588;rank=1\n\
21\t.\tbiological_region\t5000000\t5000100\t.\t+\t.\tlogic_name=cpg\n";
        let genes = try_parse_as(lines, Reading::EnsemblCache).unwrap();
        assert_eq!(genes.len(), 2);
        let lnc = transcript(&genes, "ENST00000623753");
        assert_eq!(lnc.biotype, "lncRNA");
        assert!(lnc.is_canonical);
        assert_eq!(lnc.source, "havana_tagene");
        assert_eq!(
            lnc.tsl.as_deref(),
            Some("3 (assigned to previous version 1)")
        );
        assert_eq!(genes["ENSG00000279669"].version, Some(2));
        assert_eq!(genes["ENSG00000279669"].name, None);
        let pseudo = transcript(&genes, "ENST00000559466");
        assert_eq!(pseudo.gene_id, "ENSG00000188681");
        assert_eq!(pseudo.symbol.as_deref(), Some("TEKT4P2"));
        assert_eq!(pseudo.tsl.as_deref(), Some("NA"));
        assert_eq!(genes["ENSG00000188681"].name.as_deref(), Some("TEKT4P2"));

        let genes = try_parse_as(GRCH37_PSEUDOGENE, Reading::EnsemblCache).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(
            transcript(&genes, "ENST00000411878").biotype,
            "unprocessed_pseudogene"
        );
    }

    /// An Ensembl line of a type outside `GFF.pm`'s list is read by its type
    /// when the file is read as the `--gff` source reads it (`BaseGXF.pm`
    /// `_get_records_by_coords` keeps the listed types alone), its
    /// `transcript:` prefix notwithstanding: the TEC `unconfirmed_transcript`
    /// and the scaRNA `scRNA` lines are skipped and counted, their exons name
    /// a parent whose line is skipped, and the `gene` and `ncRNA_gene` lines
    /// are genes with no transcript. Read for the Ensembl cache's content,
    /// the same lines are the two transcripts the cache carries.
    #[test]
    fn an_unlisted_type_with_an_ensembl_id_is_skipped_by_the_source_reading() {
        let lines = "21\thavana\tgene\t18917213\t18917429\t.\t+\t.\tID=gene:ENSG00000279156;biotype=TEC;description=TEC;gene_id=ENSG00000279156;logic_name=havana_homo_sapiens;version=1\n\
21\thavana\tunconfirmed_transcript\t18917213\t18917429\t.\t+\t.\tID=transcript:ENST00000623771;Parent=gene:ENSG00000279156;biotype=TEC;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000623771;transcript_support_level=NA;version=1\n\
21\thavana\texon\t18917213\t18917429\t.\t+\t.\tParent=transcript:ENST00000623771;Name=ENSE00003758460;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003758460;rank=1;version=1\n\
1\tensembl\tncRNA_gene\t9082696\t9082879\t.\t+\t.\tID=gene:ENSG00000252404;biotype=scaRNA;gene_id=ENSG00000252404;logic_name=ncrna_homo_sapiens;version=1\n\
1\tensembl\tscRNA\t9082696\t9082879\t.\t+\t.\tID=transcript:ENST00000516595;Parent=gene:ENSG00000252404;biotype=scaRNA;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000516595;transcript_support_level=NA;version=1\n\
1\tensembl\texon\t9082696\t9082879\t.\t+\t.\tParent=transcript:ENST00000516595;Name=ENSE00002088872;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002088872;rank=1;version=1\n";
        let genes = parse(lines);
        assert_eq!(genes.len(), 2);
        assert!(genes["ENSG00000279156"].transcripts.is_empty());
        assert!(genes["ENSG00000252404"].transcripts.is_empty());

        let genes = try_parse_as(lines, Reading::EnsemblCache).unwrap();
        assert_eq!(genes.len(), 2);
        let tec = transcript(&genes, "ENST00000623771");
        assert_eq!(tec.biotype, "TEC");
        assert!(tec.is_canonical);
        assert_eq!(tec.exons.len(), 1);
        assert_eq!(transcript(&genes, "ENST00000516595").biotype, "scaRNA");
    }

    /// An exon line's `ensembl_phase` and `ensembl_end_phase` are the exon's
    /// own phases (`Exon.pm` `summary_as_hash`), the values the Ensembl
    /// transcript cache stores on each exon: -1 on a non-coding exon, and on
    /// ENST00000622113 exon 4 starts in the 5' UTR (-1) and ends in frame 2.
    #[test]
    fn ensembl_phases_read_from_exon_line() {
        let lines = format!(
            "{TPTE_GENE}{TPTE_205}\
21\tensembl_havana\texon\t10521553\t10521694\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00003729301;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003729301;rank=1;version=1\n\
21\tensembl_havana\texon\t10538681\t10538734\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00003735534;constitutive=0;ensembl_end_phase=2;ensembl_phase=-1;exon_id=ENSE00003735534;rank=4;version=1\n\
21\tensembl_havana\texon\t10541112\t10541165\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00002295075;constitutive=0;ensembl_end_phase=2;ensembl_phase=2;exon_id=ENSE00002295075;rank=5;version=2\n"
        );
        let genes = parse(&lines);
        let exons = &transcript(&genes, "ENST00000622113").exons;
        let phases: Vec<(Option<i8>, Option<i8>)> = exons
            .iter()
            .map(|e| (e.ensembl_phase, e.ensembl_end_phase))
            .collect();
        assert_eq!(
            phases,
            vec![
                (Some(-1), Some(-1)),
                (Some(-1), Some(2)),
                (Some(2), Some(2))
            ]
        );
        assert_eq!(exons[1].exon_id.as_deref(), Some("ENSE00003735534"));
        assert_eq!(exons[1].rank, Some(4));
    }

    /// A transcript line without an exon or CDS line is no transcript
    /// (`BaseGXF.pm` `_create_transcript`, `return unless
    /// $tr_record->{_children}`), whatever its shape.
    #[test]
    fn a_transcript_without_exon_or_cds_lines_is_dropped() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}"));
        assert!(genes["ENSG00000274391"].transcripts.is_empty());
    }

    /// A GFF3 that is not Ensembl's follows `GFF.pm` `_record_get_biotype`
    /// (lines 305-340): an `mRNA` is `protein_coding`; an `ncRNA` takes its
    /// `ncrna_class` (`lncRNA`); a `V_gene_segment` is `IG_V_gene`, the
    /// segment letter upper-cased into the Ensembl biotype; a transcript of
    /// a gene whose `description` opens `microRNA` is `miRNA`; a `transcript`
    /// with none of those takes its `gbkey` (`misc_RNA`); and a record with
    /// none of them has no biotype and is skipped, as Ensembl VEP skips it
    /// with `Unable to determine biotype`. A `miRNA` line under a transcript
    /// is a child Ensembl VEP passes over (`_create_transcript`, `mirna`), so
    /// it makes no transcript of its own, and a `region` or `cDNA_match` line
    /// is not read at all.
    #[test]
    fn refseq_gff_biotypes_follow_gff_pm_rules() {
        let genes = parse(REFSEQ);
        let mut biotypes: Vec<(&str, &str)> = genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .map(|t| (t.transcript_id.as_str(), t.biotype.as_str()))
            .collect();
        biotypes.sort();
        assert_eq!(
            biotypes,
            [
                ("NM_000546.6", "protein_coding"),
                ("NR_1", "lncRNA"),
                ("NR_3", "miRNA"),
                ("NR_4", "misc_RNA"),
                ("id-IGHV1", "IG_V_gene"),
            ]
        );
        assert!(genes["gene-NOTYPE"].transcripts.is_empty());
        assert!(!genes.contains_key("NC_000017.11:1..83257441"));
        assert!(!genes.contains_key("aln1"));
    }

    /// The gene identifier a RefSeq transcript prints is the `GeneID` of its
    /// `Dbxref` list, `7157` for TP53 (`BaseGXF.pm` `_add_identifiers`, lines
    /// 638-660: the transcript line's `gene_id`, else `GeneID`, else the gene
    /// line's `gene_id`, else the gene line's identifier); its symbol is the
    /// transcript line's `gene`, `TP53`, else the gene line's `Name` (lines
    /// 663-671). The CDS line's `protein_id` names the protein and, with no
    /// `version` on the line, its version is 1 (`_add_translation`, line 730).
    #[test]
    fn gene_id_falls_back_to_dbxref_geneid_and_symbol_to_the_gene_name() {
        let genes = parse(REFSEQ);
        let tp53 = transcript(&genes, "NM_000546.6");
        assert_eq!(tp53.gene_id, "7157");
        assert_eq!(tp53.symbol.as_deref(), Some("TP53"));
        assert_eq!(tp53.version, None);
        assert_eq!(tp53.exons.len(), 2);
        assert_eq!(tp53.cds_regions.len(), 1);
        assert_eq!(
            tp53.cds_regions[0].protein_id.as_deref(),
            Some("NP_000537.3")
        );
        assert_eq!(tp53.cds_regions[0].protein_version, Some(1));
        assert_eq!(genes["gene-TP53"].name.as_deref(), Some("TP53"));
        assert!(!tp53.mane_select, "`MANE Select` is not the Ensembl tag");

        // Without `gene` on the transcript line the symbol is the gene line's
        // `Name`; without `GeneID` the gene identifier is the gene line's.
        let plain = REFSEQ
            .replace(";gene=TP53;product=tumor", ";product=tumor")
            .replace(
                "Dbxref=GeneID:7157,Genbank:NM_000546.6,HGNC:HGNC:11998;",
                "",
            );
        let genes = parse(&plain);
        let tp53 = transcript(&genes, "NM_000546.6");
        assert_eq!(tp53.gene_id, "gene-TP53");
        assert_eq!(tp53.symbol.as_deref(), Some("TP53"));
    }

    /// Records are placed by their own identifiers, whatever prefix they
    /// carry: `rna1` is a transcript of `g1` because its `Parent` names
    /// `g1`, and an exon whose `Parent` list is `rna1,rna2` joins both
    /// transcripts (`GFF.pm` `_record_get_parent_ids`, lines 249-252). A
    /// record whose parent is in no line is dropped, as Ensembl VEP drops it
    /// with its `Parent entries ... not found` warning, and a record with no
    /// parent whose type names no gene is a transcript with no gene line,
    /// its gene identifier its own `gene_id`.
    #[test]
    fn unprefixed_ids_and_parent_lists_build_the_tree() {
        let lines = "1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1;Name=G1;gene_id=G1\n\
1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=rna1;Parent=g1\n\
1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=rna2;Parent=g1\n\
1\tsrc\texon\t100\t200\t.\t+\t.\tID=e1;Parent=rna1,rna2\n\
1\tsrc\texon\t300\t500\t.\t+\t.\tID=e2;Parent=rna2\n\
1\tsrc\texon\t700\t900\t.\t+\t.\tID=e3;Parent=rna1\n\
1\tsrc\tmRNA\t2000\t2500\t.\t+\t.\tID=orphan;Parent=nowhere\n\
1\tsrc\texon\t2000\t2500\t.\t+\t.\tID=e4;Parent=orphan\n\
1\tsrc\tmRNA\t3000\t3500\t.\t-\t.\tID=alone;gene_id=GENELESS;biotype=lncRNA\n\
1\tsrc\texon\t3000\t3500\t.\t-\t.\tParent=alone\n";
        let genes = parse(lines);
        let rna1 = transcript(&genes, "rna1");
        assert_eq!(rna1.gene_id, "G1");
        assert_eq!(rna1.symbol.as_deref(), Some("G1"));
        assert_eq!(
            rna1.exons.iter().map(|e| e.start).collect::<Vec<_>>(),
            [100, 700]
        );
        let rna2 = transcript(&genes, "rna2");
        assert_eq!(
            rna2.exons.iter().map(|e| e.start).collect::<Vec<_>>(),
            [100, 300]
        );
        assert_eq!(genes["g1"].transcripts.len(), 2);
        assert!(genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .all(|t| t.transcript_id != "orphan"));
        let alone = transcript(&genes, "alone");
        assert_eq!(alone.gene_id, "GENELESS");
        assert_eq!(alone.symbol, None);
        assert_eq!(genes["GENELESS"].strand, -1);
        assert_eq!(genes["GENELESS"].version, None);
    }

    /// A repeated identifier names its last record (`BaseGXF.pm`
    /// `_get_parent_child_structure`, `$top_level{$record_id} = $record` over
    /// the records in file order): the second `g1` line is the gene, with its
    /// own coordinates, and `rna1` and its exon hang on it once.
    #[test]
    fn a_repeated_identifier_is_its_last_record() {
        let lines = "1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1;Name=FIRST\n\
1\tsrc\tgene\t2000\t2900\t.\t+\t.\tID=g1;Name=SECOND\n\
1\tsrc\tmRNA\t2000\t2900\t.\t+\t.\tID=rna1;Parent=g1\n\
1\tsrc\texon\t2000\t2900\t.\t+\t.\tParent=rna1\n";
        let genes = parse(lines);
        assert_eq!(genes.len(), 1);
        let gene = &genes["g1"];
        assert_eq!((gene.start, gene.end), (2000, 2900));
        assert_eq!(gene.name.as_deref(), Some("SECOND"));
        assert_eq!(gene.transcripts.len(), 1);
        assert_eq!(gene.transcripts[0].symbol.as_deref(), Some("SECOND"));
    }

    /// A transcript line with no `Parent` whose `gene_id` names a gene line is
    /// a transcript of its own (`BaseGXF.pm` `_create_transcripts`, a record
    /// without a parent) and prints that gene identifier, whichever line comes
    /// first in the file; the gene line's record carries it beside the gene's
    /// own transcripts.
    #[test]
    fn a_parentless_transcript_naming_a_gene_line_joins_its_record() {
        let gene = "1\tsrc\tgene\t100\t900\t.\t+\t.\tID=G1;Name=G1\n";
        let child = "1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=rna1;Parent=G1\n\
1\tsrc\texon\t100\t900\t.\t+\t.\tParent=rna1\n";
        let loose = "1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=rna0;gene_id=G1;biotype=lncRNA\n\
1\tsrc\texon\t100\t500\t.\t+\t.\tParent=rna0\n";
        for lines in [
            format!("{gene}{child}{loose}"),
            format!("{loose}{gene}{child}"),
        ] {
            let genes = parse(&lines);
            assert_eq!(genes.len(), 1, "{lines}");
            let mut ids: Vec<&str> = genes["G1"]
                .transcripts
                .iter()
                .map(|t| t.transcript_id.as_str())
                .collect();
            ids.sort_unstable();
            assert_eq!(ids, ["rna0", "rna1"], "{lines}");
            assert_eq!(genes["G1"].name.as_deref(), Some("G1"), "{lines}");
            assert_eq!(transcript(&genes, "rna0").gene_id, "G1");
        }
    }

    /// An attribute whose value is `0` is absent, as Ensembl VEP's truth
    /// tests read it, so each reader falls to its next alternative: the gene
    /// line `ID=0;Name=g1` is named `g1` (`GFF.pm` line 278, `ID || Name ||
    /// id || name`) and its `version=0` is no version (`BaseGXF.pm` line 639
    /// and `OutputFactory.pm` line 1478); the transcript's `Parent=0,g1`
    /// names `g1` alone (`BaseGXF.pm` line 389, `next unless ($parent_id)`),
    /// its gene identifier passes `gene_id=0`, `GeneID:0` and the gene line's
    /// `gene_id=0` to reach the gene line's identifier (lines 638-659), its
    /// symbol passes `gene_name=0` and `gene=0` to the gene line's `Name`
    /// (lines 662-670), `biotype=0` leaves an `mRNA` `protein_coding`
    /// (`GFF.pm` lines 311-320), `version=0` yields to `transcript_version`
    /// (`BaseGXF.pm` line 481), `transcript_support_level=0` to `tsl` (line
    /// 674), the CDS line's `protein_id=0` is no protein (line 723) and its
    /// `version=0` is version 1 (line 727). A `Parent=0` alone is no parent
    /// (`GFF.pm` lines 249 and 252, the `||` then the truth test before the
    /// split), so the record is a transcript of its own, and
    /// `transcript_id=0` leaves the record its `ID` (`BaseGXF.pm` line 464).
    /// A `Parent=0,0` is a parent list with no set identifier: the record
    /// is linked to nothing (`BaseGXF.pm` line 389) and, opening with a
    /// field, is no record of its own (lines 430-435), so it is dropped;
    /// `Parent=,0` opens with an empty field, so the prune leaves the record
    /// a transcript of its own (line 432, `$parent_ids->[0] ne ''`).
    #[test]
    fn a_zero_value_is_an_absent_attribute() {
        let lines = "1\tsrc\tgene\t100\t900\t.\t+\t.\tID=0;Name=g1;gene_id=0;version=0\n\
1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=rna1;Parent=0,g1;gene_id=0;Dbxref=GeneID:0,HGNC:HGNC:1;gene_name=0;gene=0;biotype=0;version=0;transcript_version=7;transcript_support_level=0;tsl=3\n\
1\tsrc\texon\t100\t900\t.\t+\t.\tParent=rna1\n\
1\tsrc\tCDS\t200\t400\t.\t+\t0\tParent=rna1;protein_id=0;version=0\n\
1\tsrc\tmRNA\t2000\t2900\t.\t-\t.\tID=rna2;Parent=0;transcript_id=0;gene_id=G2;biotype=lncRNA\n\
1\tsrc\texon\t2000\t2900\t.\t-\t.\tParent=rna2\n\
1\tsrc\tmRNA\t4000\t4900\t.\t+\t.\tID=rna3;Parent=0,0;biotype=lncRNA\n\
1\tsrc\texon\t4000\t4900\t.\t+\t.\tParent=rna3\n\
1\tsrc\tmRNA\t6000\t6900\t.\t+\t.\tID=rna4;Parent=,0;biotype=lncRNA\n\
1\tsrc\texon\t6000\t6900\t.\t+\t.\tParent=rna4\n";
        let genes = parse(lines);
        assert_eq!(genes.len(), 3, "{:?}", genes.keys());
        let mut ids: Vec<&str> = genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .map(|t| t.transcript_id.as_str())
            .collect();
        ids.sort_unstable();
        assert_eq!(ids, ["rna1", "rna2", "rna4"]);
        let gene = &genes["g1"];
        assert_eq!(gene.name.as_deref(), Some("g1"));
        assert_eq!(gene.version, None);
        assert_eq!(gene.transcripts.len(), 1);
        assert_eq!(transcript(&genes, "rna4").gene_id, "");
        let tr = transcript(&genes, "rna1");
        assert_eq!(tr.gene_id, "g1");
        assert_eq!(tr.symbol.as_deref(), Some("g1"));
        assert_eq!(tr.biotype, "protein_coding");
        assert_eq!(tr.version, Some(7));
        assert_eq!(tr.tsl.as_deref(), Some("3"));
        assert_eq!(tr.cds_regions.len(), 1);
        assert_eq!(tr.cds_regions[0].protein_id, None);
        assert_eq!(tr.cds_regions[0].protein_version, Some(1));
        let alone = transcript(&genes, "rna2");
        assert_eq!(alone.gene_id, "G2");
        assert_eq!(genes["G2"].transcripts.len(), 1);
    }

    /// A parent list opening with an empty field leaves its record at the top
    /// level for Ensembl VEP (`BaseGXF.pm` line 432, `$parent_ids->[0] ne
    /// ''`) while line 391 links it under each set parent in a line, so a
    /// transcript naming a gene that way is met twice there, its gene record
    /// the gene or none by the order the top level is read; here `rna5`
    /// (`Parent=,g1`) is `g1`'s transcript once. The record is a transcript
    /// of its own when its parent is in no line (`rna6`, `Parent=,x`, the
    /// `x` warned about) or is a transcript (`mir1`, `Parent=,rna1`, which
    /// `_create_transcript` passes over as a `mirna` child, lines 513-514),
    /// and a gene of its own when its type names one (`g2`, `Parent=,g1`,
    /// whose child `rna7` is its transcript): `_create_transcripts` reads a
    /// top-level gene's children off, line 341, before `_create_transcript`
    /// meets the copy under `g1` and returns at line 462, so `g2`'s own exon
    /// line is no transcript of `g1`.
    #[test]
    fn a_parent_list_opening_with_an_empty_field() {
        let lines = "1\tsrc\tgene\t100\t9000\t.\t+\t.\tID=g1\n\
1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=rna1;Parent=g1;biotype=lncRNA\n\
1\tsrc\texon\t100\t900\t.\t+\t.\tParent=rna1\n\
1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=rna5;Parent=,g1;biotype=lncRNA\n\
1\tsrc\texon\t100\t500\t.\t+\t.\tParent=rna5\n\
1\tsrc\tmRNA\t2000\t2900\t.\t+\t.\tID=rna6;Parent=,x;biotype=lncRNA\n\
1\tsrc\texon\t2000\t2900\t.\t+\t.\tParent=rna6\n\
1\tsrc\tmiRNA\t300\t320\t.\t+\t.\tID=mir1;Parent=,rna1;biotype=miRNA\n\
1\tsrc\texon\t300\t320\t.\t+\t.\tParent=mir1\n\
1\tsrc\tgene\t5000\t5900\t.\t+\t.\tID=g2;Parent=,g1;biotype=lncRNA\n\
1\tsrc\texon\t5000\t5900\t.\t+\t.\tParent=g2\n\
1\tsrc\tlnc_RNA\t5000\t5900\t.\t+\t.\tID=rna7;Parent=g2;biotype=lncRNA\n\
1\tsrc\texon\t5000\t5900\t.\t+\t.\tParent=rna7\n";
        let genes = parse(lines);
        let mut placed: Vec<(&str, &str)> = genes
            .iter()
            .flat_map(|(key, g)| {
                g.transcripts
                    .iter()
                    .map(move |t| (key.as_str(), t.transcript_id.as_str()))
            })
            .collect();
        placed.sort_unstable();
        assert_eq!(
            placed,
            [
                ("g1", "rna1"),
                ("g1", "rna5"),
                ("g2", "rna7"),
                ("transcript:mir1", "mir1"),
                ("transcript:rna6", "rna6"),
            ]
        );
        assert_eq!(transcript(&genes, "rna5").gene_id, "g1");
        assert_eq!(transcript(&genes, "rna7").gene_id, "g2");
    }

    /// The transcript identifier is `transcript_id`, else the record's
    /// identifier, with a `gene:` or `transcript:` prefix removed from
    /// either (`BaseGXF.pm` lines 464-465, `$id =~ s/^(gene|transcript)://i`
    /// on the `||` result).
    #[test]
    fn transcript_id_attribute_loses_its_prefix() {
        let lines = "1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1\n\
1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=rna1;Parent=g1;transcript_id=transcript:T1\n\
1\tsrc\texon\t100\t900\t.\t+\t.\tParent=rna1\n";
        let genes = parse(lines);
        assert_eq!(genes["g1"].transcripts[0].transcript_id, "T1");
    }

    /// A kept line that would draw two refusals draws the one read first, its
    /// start and end columns before its attribute column and its attribute
    /// column before its version, under either reading: an `exon` line whose
    /// `Name` decodes to bytes that are not UTF-8 is refused for a start or
    /// an end that is no number and for an end before its start, and a `CDS`
    /// line and a gene line whose `version` is no whole number are refused
    /// for such an attribute.
    #[test]
    fn a_kept_line_draws_the_refusal_of_the_column_read_first() {
        let exon = |span: &str| {
            let exon = TPTE_205_EXON.replace("Name=ENSE00003729301", "Name=caf%FF");
            exon.replace("\t10521553\t10521694\t", span)
        };
        let cds = TPTE_205_CDS.replacen(
            "protein_id=ENSP00000482040;version=1",
            "protein_id=caf%FF;version=v2",
            1,
        );
        let gene = TPTE_GENE
            .replace("Name=TPTE", "Name=caf%FF")
            .replace(";version=6", ";version=v2");
        for (lines, refusal) in [
            (
                format!("{TPTE_GENE}{TPTE_205}{}", exon("\tabc\t10521694\t")),
                "line 3: GFF3 start \"abc\" is not a number",
            ),
            (
                format!("{TPTE_GENE}{TPTE_205}{}", exon("\t10521553\txyz\t")),
                "line 3: GFF3 end \"xyz\" is not a number",
            ),
            (
                format!("{TPTE_GENE}{TPTE_205}{}", exon("\t10521694\t10521553\t")),
                "line 3: GFF3 coordinates 10521694-10521553 are not 1-based and ordered",
            ),
            (
                format!("{TPTE_GENE}{TPTE_205}{TPTE_205_EXON}{cds}"),
                "line 4: GFF3 attribute protein_id \"caf%FF\" decodes to text that is not UTF-8",
            ),
            (
                format!("{gene}{TPTE_205}{TPTE_205_EXON}"),
                "line 1: GFF3 attribute Name \"caf%FF\" decodes to text that is not UTF-8",
            ),
        ] {
            for reading in [Reading::GffSource, Reading::EnsemblCache] {
                let err = try_parse_as(&lines, reading).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gff3 {refusal}")),
                    "{reading:?}: {err:#}"
                );
            }
        }
    }

    /// An attribute value whose escapes spell bytes that are not UTF-8
    /// refuses the parse on a line of an unlisted type only where the reading
    /// keeps the line: a `region` line's `Note` refuses nothing under either
    /// reading; the `Name` of an `unconfirmed_transcript` line whose `ID` is
    /// `transcript:T2` refuses the parse under [`Reading::EnsemblCache`],
    /// naming the file, the line, the key and the value, and refuses nothing
    /// under the `--gff` reading.
    #[test]
    fn an_unlisted_lines_undecodable_value_refuses_where_the_line_is_kept() {
        let region = format!(
            "{}21\tx\tregion\t1\t5000\t.\t.\t.\tID=21;Note=caf%FF\n",
            three_lines("21", "1")
        );
        for reading in [Reading::GffSource, Reading::EnsemblCache] {
            let genes =
                try_parse_as(&region, reading).unwrap_or_else(|err| panic!("{reading:?}: {err:#}"));
            assert_eq!(transcript(&genes, "T1").exons.len(), 1, "{reading:?}");
        }
        let tec = format!(
            "{}21\tx\tunconfirmed_transcript\t1000\t2000\t.\t+\t.\tID=transcript:T2;Parent=gene:G1;Name=caf%FF\n",
            three_lines("21", "1")
        );
        let err = try_parse_as(&tec, Reading::EnsemblCache).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "t.gff3 line 4: GFF3 attribute Name \"caf%FF\" decodes to text that is not UTF-8"
            ),
            "{err:#}"
        );
        let genes = try_parse(&tec).unwrap();
        assert_eq!(genes["G1"].transcripts.len(), 1);
    }

    /// A gene line that would draw a refusal for its `version` and one for
    /// its source column draws the version's, the one read first, under
    /// either reading: a gene line whose `version` is no whole number and
    /// whose source decodes to bytes that are not UTF-8 is refused for the
    /// version.
    #[test]
    fn a_gene_line_is_refused_for_its_version_before_its_source() {
        let gene = TPTE_GENE
            .replace("\tensembl_havana\t", "\tcaf%FF\t")
            .replace(";version=6", ";version=v2");
        for reading in [Reading::GffSource, Reading::EnsemblCache] {
            let err = try_parse_as(&gene, reading).unwrap_err();
            assert!(
                format!("{err:#}").ends_with("t.gff3 line 1: version \"v2\" is not a whole number"),
                "{reading:?}: {err:#}"
            );
        }
    }
}
