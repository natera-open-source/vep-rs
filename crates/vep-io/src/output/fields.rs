// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The VEP output field table shared by every text format.
//!
//! VEP renders one hash of named fields per (variant, allele, feature) and then
//! lays that hash out per format: the default format keeps thirteen fixed
//! columns and folds every other key into `Extra` as `KEY=VALUE` pairs, `--tab`
//! prints every key as its own column, VCF joins the CSQ subfields with `|`, and
//! JSON lowercases the keys. Which keys exist is decided by the command-line
//! flags in a fixed order (`@FLAG_FIELDS` in `Bio/EnsEMBL/VEP/Constants.pm`,
//! ensembl-vep release/115), and the `user` group (IMPACT, DISTANCE, STRAND,
//! FLAGS) is always on. This module is that table: the column names, their
//! descriptions, the flag-to-field order, the value of any field for a row, and
//! the order in which one input record's rows come out.

use std::collections::BTreeSet;

use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::variant::{InputVariant, VariantClass};

/// The thirteen fixed columns of the default and tab formats, in order.
pub const DEFAULT_OUTPUT_COLS: [&str; 13] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
];

/// Which optional field groups are switched on. Field order across groups is
/// fixed by [`flag_fields`]; a group's fields appear only when its flag is set.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FieldOptions {
    pub allele_number: bool,
    pub show_ref_allele: bool,
    pub uploaded_allele: bool,
    pub include_pick: bool,
    pub variant_class: bool,
    pub minimal: bool,
    pub symbol: bool,
    pub biotype: bool,
    pub canonical: bool,
    pub mane_select: bool,
    pub mane: bool,
    pub tsl: bool,
    pub gencode_primary: bool,
    pub appris: bool,
    pub ccds: bool,
    pub protein: bool,
    pub uniprot: bool,
    pub xref_refseq: bool,
    pub gene_phenotype: bool,
    pub sift: bool,
    pub polyphen: bool,
    pub numbers: bool,
    pub domains: bool,
    pub mirna: bool,
    pub hgvs: bool,
    pub hgvsg: bool,
    pub af: bool,
    pub af_1kg: bool,
    pub af_gnomade: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
    pub check_existing: bool,
    pub pubmed: bool,
    pub dont_skip: bool,
    pub overlaps: bool,
    pub regulatory: bool,
    /// VEP's `--no_escape` (implied by `--json`): HGVSp keeps its `=` instead of
    /// the `%3D` the text formats print.
    pub no_escape: bool,
}

/// The flag-driven field groups in VEP's order (`@FLAG_FIELDS`). The `user`
/// group is unconditional and comes first among the groups that can fire here.
const FLAG_FIELD_GROUPS: &[(&str, &[&str])] = &[
    ("allele_number", &["ALLELE_NUM"]),
    ("show_ref_allele", &["REF_ALLELE"]),
    ("uploaded_allele", &["UPLOADED_ALLELE"]),
    ("user", &["IMPACT", "DISTANCE", "STRAND", "FLAGS"]),
    ("include_pick", &["PICK"]),
    ("variant_class", &["VARIANT_CLASS"]),
    ("minimal", &["MINIMISED"]),
    ("symbol", &["SYMBOL", "SYMBOL_SOURCE", "HGNC_ID"]),
    ("biotype", &["BIOTYPE"]),
    ("canonical", &["CANONICAL"]),
    ("mane_select", &["MANE", "MANE_SELECT"]),
    ("mane", &["MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL"]),
    ("tsl", &["TSL"]),
    ("gencode_primary", &["GENCODE_PRIMARY"]),
    ("appris", &["APPRIS"]),
    ("ccds", &["CCDS"]),
    ("protein", &["ENSP"]),
    (
        "uniprot",
        &["SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"],
    ),
    ("xref_refseq", &["RefSeq"]),
    ("gene_phenotype", &["GENE_PHENO"]),
    ("sift", &["SIFT"]),
    ("polyphen", &["PolyPhen"]),
    ("numbers", &["EXON", "INTRON"]),
    ("domains", &["DOMAINS"]),
    ("mirna", &["miRNA"]),
    ("hgvs", &["HGVSc", "HGVSp", "HGVS_OFFSET"]),
    ("hgvsg", &["HGVSg"]),
    ("af", &["AF"]),
    (
        "af_1kg",
        &["AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"],
    ),
    (
        "af_gnomade",
        &[
            "gnomADe_AF",
            "gnomADe_AFR_AF",
            "gnomADe_AMR_AF",
            "gnomADe_ASJ_AF",
            "gnomADe_EAS_AF",
            "gnomADe_FIN_AF",
            "gnomADe_MID_AF",
            "gnomADe_NFE_AF",
            "gnomADe_REMAINING_AF",
            "gnomADe_SAS_AF",
        ],
    ),
    (
        "af_gnomadg",
        &[
            "gnomADg_AF",
            "gnomADg_AFR_AF",
            "gnomADg_AMI_AF",
            "gnomADg_AMR_AF",
            "gnomADg_ASJ_AF",
            "gnomADg_EAS_AF",
            "gnomADg_FIN_AF",
            "gnomADg_MID_AF",
            "gnomADg_NFE_AF",
            "gnomADg_REMAINING_AF",
            "gnomADg_SAS_AF",
        ],
    ),
    ("max_af", &["MAX_AF", "MAX_AF_POPS"]),
    ("check_existing", &["CLIN_SIG", "SOMATIC", "PHENO"]),
    ("pubmed", &["PUBMED"]),
    ("dont_skip", &["CHECK_REF"]),
    ("overlaps", &["OverlapBP", "OverlapPC"]),
    (
        "regulatory",
        &[
            "MOTIF_NAME",
            "MOTIF_POS",
            "HIGH_INF_POS",
            "MOTIF_SCORE_CHANGE",
            "TRANSCRIPTION_FACTORS",
        ],
    ),
];

impl FieldOptions {
    fn group_active(&self, group: &str) -> bool {
        match group {
            "user" => true,
            "allele_number" => self.allele_number,
            "show_ref_allele" => self.show_ref_allele,
            "uploaded_allele" => self.uploaded_allele,
            "include_pick" => self.include_pick,
            "variant_class" => self.variant_class,
            "minimal" => self.minimal,
            "symbol" => self.symbol,
            "biotype" => self.biotype,
            "canonical" => self.canonical,
            "mane_select" => self.mane_select,
            "mane" => self.mane,
            "tsl" => self.tsl,
            "gencode_primary" => self.gencode_primary,
            "appris" => self.appris,
            "ccds" => self.ccds,
            "protein" => self.protein,
            "uniprot" => self.uniprot,
            "xref_refseq" => self.xref_refseq,
            "gene_phenotype" => self.gene_phenotype,
            "sift" => self.sift,
            "polyphen" => self.polyphen,
            "numbers" => self.numbers,
            "domains" => self.domains,
            "mirna" => self.mirna,
            "hgvs" => self.hgvs,
            "hgvsg" => self.hgvsg,
            "af" => self.af,
            "af_1kg" => self.af_1kg,
            "af_gnomade" => self.af_gnomade,
            "af_gnomadg" => self.af_gnomadg,
            "max_af" => self.max_af,
            "check_existing" => self.check_existing,
            "pubmed" => self.pubmed,
            "dont_skip" => self.dont_skip,
            "overlaps" => self.overlaps,
            "regulatory" => self.regulatory,
            _ => false,
        }
    }

    /// True when any group beyond the unconditional `user` group is on, i.e.
    /// when a row can carry a field that is not IMPACT, DISTANCE, STRAND or FLAGS.
    pub fn any_optional_group(&self) -> bool {
        FLAG_FIELD_GROUPS
            .iter()
            .any(|(g, _)| *g != "user" && self.group_active(g))
    }
}

/// The row-invariant part of a configuration's `Extra` column, computed once
/// per run so the per-row writer neither walks the flag table nor rebuilds the
/// flag-field list.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtraFieldsPlan {
    pub options: FieldOptions,
    /// [`flag_fields`] of `options`.
    pub flag_fields: Vec<&'static str>,
    /// [`FieldOptions::any_optional_group`] of `options`.
    pub any_optional_group: bool,
}

impl ExtraFieldsPlan {
    pub fn new(options: FieldOptions) -> Self {
        let flag_fields = flag_fields(&options);
        let any_optional_group = options.any_optional_group();
        Self {
            options,
            flag_fields,
            any_optional_group,
        }
    }
}

/// The optional fields a configuration emits, in VEP's order, each name once.
pub fn flag_fields(options: &FieldOptions) -> Vec<&'static str> {
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    let mut out = Vec::new();
    for (group, fields) in FLAG_FIELD_GROUPS {
        if !options.group_active(group) {
            continue;
        }
        for f in *fields {
            if seen.insert(f) {
                out.push(*f);
            }
        }
    }
    out
}

/// VEP's description of a field, for the `## Column descriptions:` and
/// `## Extra column keys:` header blocks.
pub fn field_description(field: &str) -> Option<&'static str> {
    FIELD_DESCRIPTIONS
        .iter()
        .find(|(k, _)| *k == field)
        .map(|(_, d)| *d)
}

const FIELD_DESCRIPTIONS: &[(&str, &str)] = &[
    ("Uploaded_variation", "Identifier of uploaded variant"),
    (
        "Location",
        "Location of variant in standard coordinate format (chr:start or chr:start-end)",
    ),
    ("Allele", "The variant allele used to calculate the consequence"),
    ("Gene", "Stable ID of affected gene"),
    ("Feature", "Stable ID of feature"),
    (
        "Feature_type",
        "Type of feature - Transcript, RegulatoryFeature or MotifFeature",
    ),
    ("Consequence", "Consequence type"),
    ("cDNA_position", "Relative position of base pair in cDNA sequence"),
    ("CDS_position", "Relative position of base pair in coding sequence"),
    ("Protein_position", "Relative position of amino acid in protein"),
    ("Amino_acids", "Reference and variant amino acids"),
    ("Codons", "Reference and variant codon sequence"),
    ("Existing_variation", "Identifier(s) of co-located known variants"),
    ("IMPACT", "Subjective impact classification of consequence type"),
    ("CANONICAL", "Indicates if transcript is canonical for this gene"),
    (
        "MANE",
        "MANE (Matched Annotation from NCBI and EMBL-EBI) set(s) the transcript belongs to",
    ),
    (
        "MANE_SELECT",
        "MANE Select (Matched Annotation from NCBI and EMBL-EBI) Transcript",
    ),
    (
        "MANE_PLUS_CLINICAL",
        "MANE Plus Clinical (Matched Annotation from NCBI and EMBL-EBI) Transcript",
    ),
    ("TSL", "Transcript support level"),
    (
        "APPRIS",
        "Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods",
    ),
    ("CCDS", "Indicates if transcript is a CCDS transcript"),
    ("SYMBOL", "Gene symbol (e.g. HGNC)"),
    ("SYMBOL_SOURCE", "Source of gene symbol"),
    ("SOURCE", "Source of transcript"),
    ("HGNC_ID", "Stable identifer of HGNC gene symbol"),
    ("ENSP", "Protein identifer"),
    ("FLAGS", "Transcript quality flags"),
    ("SWISSPROT", "UniProtKB/Swiss-Prot accession"),
    ("TREMBL", "UniProtKB/TrEMBL accession"),
    ("UNIPARC", "UniParc accession"),
    ("UNIPROT_ISOFORM", "Direct mappings to UniProtKB isoforms"),
    ("miRNA", "SO terms of overlapped miRNA secondary structure feature(s)"),
    ("HGVSc", "HGVS coding sequence name"),
    ("HGVSp", "HGVS protein sequence name"),
    ("HGVSg", "HGVS genomic sequence name"),
    ("SIFT", "SIFT prediction and/or score"),
    ("PolyPhen", "PolyPhen prediction and/or score"),
    ("EXON", "Exon number(s) / total"),
    ("INTRON", "Intron number(s) / total"),
    (
        "DOMAINS",
        "The source and identifer of any overlapping protein domains",
    ),
    (
        "MOTIF_NAME",
        "The stable identifier of a transcription factor binding profile (TFBP) aligned at this position",
    ),
    (
        "MOTIF_POS",
        "The relative position of the variation in the aligned TFBP",
    ),
    (
        "HIGH_INF_POS",
        "A flag indicating if the variant falls in a high information position of the TFBP",
    ),
    (
        "MOTIF_SCORE_CHANGE",
        "The difference in motif score of the reference and variant sequences for the TFBP",
    ),
    (
        "TRANSCRIPTION_FACTORS",
        "List of transcription factors which bind to the transcription factor binding profile",
    ),
    (
        "AF",
        "Frequency of existing variant in 1000 Genomes combined population",
    ),
    (
        "AFR_AF",
        "Frequency of existing variant in 1000 Genomes combined African population",
    ),
    (
        "AMR_AF",
        "Frequency of existing variant in 1000 Genomes combined American population",
    ),
    (
        "EAS_AF",
        "Frequency of existing variant in 1000 Genomes combined East Asian population",
    ),
    (
        "EUR_AF",
        "Frequency of existing variant in 1000 Genomes combined European population",
    ),
    (
        "SAS_AF",
        "Frequency of existing variant in 1000 Genomes combined South Asian population",
    ),
    (
        "gnomADe_AF",
        "Frequency of existing variant in gnomAD exomes combined population",
    ),
    (
        "gnomADe_AFR_AF",
        "Frequency of existing variant in gnomAD exomes African/American population",
    ),
    (
        "gnomADe_AMR_AF",
        "Frequency of existing variant in gnomAD exomes American population",
    ),
    (
        "gnomADe_ASJ_AF",
        "Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population",
    ),
    (
        "gnomADe_EAS_AF",
        "Frequency of existing variant in gnomAD exomes East Asian population",
    ),
    (
        "gnomADe_FIN_AF",
        "Frequency of existing variant in gnomAD exomes Finnish population",
    ),
    (
        "gnomADe_MID_AF",
        "Frequency of existing variant in gnomAD exomes Mid-eastern population",
    ),
    (
        "gnomADe_NFE_AF",
        "Frequency of existing variant in gnomAD exomes Non-Finnish European population",
    ),
    (
        "gnomADe_REMAINING_AF",
        "Frequency of existing variant in gnomAD exomes remaining combined populations",
    ),
    (
        "gnomADe_SAS_AF",
        "Frequency of existing variant in gnomAD exomes South Asian population",
    ),
    (
        "gnomADg_AF",
        "Frequency of existing variant in gnomAD genomes combined population",
    ),
    (
        "gnomADg_AFR_AF",
        "Frequency of existing variant in gnomAD genomes African/American population",
    ),
    (
        "gnomADg_AMI_AF",
        "Frequency of existing variant in gnomAD genomes Amish population",
    ),
    (
        "gnomADg_AMR_AF",
        "Frequency of existing variant in gnomAD genomes American population",
    ),
    (
        "gnomADg_ASJ_AF",
        "Frequency of existing variant in gnomAD genomes Ashkenazi Jewish population",
    ),
    (
        "gnomADg_EAS_AF",
        "Frequency of existing variant in gnomAD genomes East Asian population",
    ),
    (
        "gnomADg_FIN_AF",
        "Frequency of existing variant in gnomAD genomes Finnish population",
    ),
    (
        "gnomADg_MID_AF",
        "Frequency of existing variant in gnomAD genomes Mid-eastern population",
    ),
    (
        "gnomADg_NFE_AF",
        "Frequency of existing variant in gnomAD genomes Non-Finnish European population",
    ),
    (
        "gnomADg_REMAINING_AF",
        "Frequency of existing variant in gnomAD genomes remaining combined populations",
    ),
    (
        "gnomADg_SAS_AF",
        "Frequency of existing variant in gnomAD genomes South Asian population",
    ),
    (
        "MAX_AF",
        "Maximum observed allele frequency in 1000 Genomes, ESP and ExAC/gnomAD",
    ),
    (
        "MAX_AF_POPS",
        "Populations in which maximum allele frequency was observed",
    ),
    ("DISTANCE", "Shortest distance from variant to transcript"),
    ("CLIN_SIG", "ClinVar clinical significance of the dbSNP variant"),
    ("BIOTYPE", "Biotype of transcript or regulatory feature"),
    ("PUBMED", "Pubmed ID(s) of publications that cite existing variant"),
    (
        "ALLELE_NUM",
        "Allele number from input; 0 is reference, 1 is first alternate etc",
    ),
    ("REF_ALLELE", "Reference allele"),
    ("UPLOADED_ALLELE", "The variant allele as it was uploaded"),
    ("STRAND", "Strand of the feature (1/-1)"),
    (
        "PICK",
        "Indicates if this consequence has been picked as the most severe",
    ),
    ("SOMATIC", "Somatic status of existing variant"),
    ("VARIANT_CLASS", "SO variant class"),
    (
        "PHENO",
        "Indicates if existing variant(s) is associated with a phenotype, disease or trait; multiple values correspond to multiple variants",
    ),
    (
        "GENE_PHENO",
        "Indicates if gene is associated with a phenotype, disease or trait",
    ),
    (
        "MINIMISED",
        "Alleles in this variant have been converted to minimal representation before consequence calculation",
    ),
    (
        "HGVS_OFFSET",
        "Indicates by how many bases the HGVS notations for this variant have been shifted",
    ),
    (
        "OverlapBP",
        "Number of base pairs overlapping with the corresponding structural variation feature",
    ),
    (
        "OverlapPC",
        "Percentage of corresponding structural variation feature overlapped by the given input",
    ),
    (
        "CHECK_REF",
        "Reports variants where the input reference does not match the expected reference",
    ),
    ("GENCODE_PRIMARY", "Indicates if transcript is a GENCODE Primary transcript"),
];

/// One output row: an allele of the record and the transcript consequence it is
/// reported against, or `None` for the record's intergenic row.
pub type Row<'a> = (&'a InputVariant, Option<&'a TranscriptConsequence>);

/// The rows one input record produces, in VEP's order: transcript-major (by
/// transcript id) and allele-minor (ALT order) across every allele split from
/// the record, then one intergenic row per allele that overlapped nothing.
///
/// An allele whose chromosome had no transcripts at all (no `most_severe_consequence`)
/// produces no row, as VEP writes nothing for it.
pub fn record_rows<'a>(record: &[&'a InputVariant]) -> Vec<Row<'a>> {
    let mut rows = Vec::new();
    record_rows_into(record, &mut rows);
    rows
}

/// [`record_rows`] into a caller-owned vector, cleared first, so a writer that
/// renders many records reuses one allocation.
pub fn record_rows_into<'a>(record: &[&'a InputVariant], rows: &mut Vec<Row<'a>>) {
    rows.clear();
    for variant in record {
        for tc in &variant.transcript_consequences {
            rows.push((variant, Some(tc)));
        }
    }
    // Stable, so rows with one transcript and allele keep their consequence order.
    rows.sort_by(|a, b| {
        let (ta, tb) = (a.1.expect("feature row"), b.1.expect("feature row"));
        ta.transcript_id
            .as_ref()
            .cmp(tb.transcript_id.as_ref())
            .then(a.0.allele_index.cmp(&b.0.allele_index))
    });
    for variant in record {
        if variant.transcript_consequences.is_empty() && variant.most_severe_consequence.is_some() {
            rows.push((variant, None));
        }
    }
}

/// Consequence terms of a row joined by `,`.
pub fn consequence_string(tc: &TranscriptConsequence) -> String {
    let mut out = String::new();
    for (i, c) in tc.consequences.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(c.so_term());
    }
    out
}

fn first_frequency(variant: &InputVariant, pop: &str) -> Option<f64> {
    for cv in &variant.colocated_variants {
        if let Some(&f) = cv.frequencies.get(pop) {
            return Some(f);
        }
    }
    None
}

/// Frequencies print as VEP prints them (`%.4g`): four significant digits,
/// trailing zeros trimmed, scientific notation outside [1e-4, 1e4).
pub fn format_freq(f: f64) -> String {
    if f == 0.0 || !f.is_finite() {
        return "0".to_string();
    }
    let abs = f.abs();
    let mut out = if !(0.0001..10_000.0).contains(&abs) {
        format!("{f:.3e}")
    } else {
        let digits_before_decimal = abs.log10().floor() as i32 + 1;
        let decimals = (4 - digits_before_decimal).max(0) as usize;
        format!("{f:.decimals$}")
    };
    if let Some((mantissa, exponent)) = out.split_once('e') {
        let trimmed = mantissa.trim_end_matches('0').trim_end_matches('.');
        out = format!("{trimmed}e{exponent}");
    } else if out.contains('.') {
        out = out.trim_end_matches('0').trim_end_matches('.').to_string();
    }
    out
}

fn joined_unique(values: impl IntoIterator<Item = String>, sep: &str) -> String {
    let mut seen: Vec<String> = Vec::new();
    for v in values {
        if !v.is_empty() && v != "." && !seen.contains(&v) {
            seen.push(v);
        }
    }
    seen.join(sep)
}

/// The value of one named field for one row, or an empty string when the field
/// is absent on this row. `list_sep` joins list-valued fields (`,` in the default
/// and tab formats, `&` inside a VCF CSQ entry). Field names not in the table are
/// looked up in the row's plugin data.
pub fn field_value(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
    list_sep: &str,
) -> String {
    match field {
        "Uploaded_variation" => variant.uploaded_variation(),
        "Location" => variant.location(),
        "Allele" => variant.display_allele(),
        "Consequence" => tc
            .map(consequence_string)
            .or_else(|| {
                variant
                    .most_severe_consequence
                    .map(|c| c.so_term().to_string())
            })
            .unwrap_or_else(|| Consequence::IntergenicVariant.so_term().to_string()),
        "IMPACT" => tc
            .map(|x| x.impact.to_string())
            .or_else(|| {
                variant
                    .most_severe_consequence
                    .map(|c| c.impact().to_string())
            })
            .unwrap_or_else(|| Consequence::IntergenicVariant.impact().to_string()),
        "SYMBOL" => tc
            .and_then(|x| x.gene_symbol.as_deref())
            .unwrap_or("")
            .to_string(),
        "Gene" => tc.map(|x| x.gene_id.as_ref()).unwrap_or("").to_string(),
        "Feature_type" => tc.map(|x| x.feature_type.to_string()).unwrap_or_default(),
        "Feature" => tc
            .map(|x| x.transcript_id.as_ref())
            .unwrap_or("")
            .to_string(),
        "BIOTYPE" => tc
            .and_then(|x| x.biotype.as_deref())
            .unwrap_or("")
            .to_string(),
        "EXON" => tc.and_then(|x| x.exon.clone()).unwrap_or_default(),
        "INTRON" => tc.and_then(|x| x.intron.clone()).unwrap_or_default(),
        "HGVSc" if !options.hgvs => String::new(),
        "HGVSp" if !options.hgvs => String::new(),
        "HGVSc" => tc.and_then(|x| x.hgvsc.clone()).unwrap_or_default(),
        "HGVSp" => {
            let v = tc.and_then(|x| x.hgvsp.clone()).unwrap_or_default();
            if options.no_escape {
                v
            } else {
                v.replace('=', "%3D")
            }
        }
        "cDNA_position" => tc.and_then(|x| x.cdna_position.clone()).unwrap_or_default(),
        "CDS_position" => tc.and_then(|x| x.cds_position.clone()).unwrap_or_default(),
        "Protein_position" => tc
            .and_then(|x| x.protein_position.clone())
            .unwrap_or_default(),
        "Amino_acids" => tc.and_then(|x| x.amino_acids.clone()).unwrap_or_default(),
        "Codons" => tc.and_then(|x| x.codons.clone()).unwrap_or_default(),
        "Existing_variation" => variant.existing_variation.join(list_sep),
        "ALLELE_NUM" => (variant.allele_index + 1).to_string(),
        "REF_ALLELE" => String::from_utf8_lossy(&variant.ref_allele).to_string(),
        "UPLOADED_ALLELE" => variant
            .uploaded_allele_string
            .clone()
            .unwrap_or_else(|| variant.allele_string.clone()),
        "DISTANCE" => tc
            .and_then(|x| x.distance.map(|d| d.to_string()))
            .unwrap_or_default(),
        "OverlapBP" => sv_overlap(variant, tc)
            .map(|(bp, _)| bp.to_string())
            .unwrap_or_default(),
        "OverlapPC" => sv_overlap(variant, tc)
            .map(|(_, pc)| format!("{pc:.2}"))
            .unwrap_or_default(),
        "STRAND" => tc.map(|x| x.strand.to_string()).unwrap_or_default(),
        "FLAGS" => tc
            .map(|x| {
                x.flags
                    .iter()
                    .filter(|f| f.starts_with("cds_"))
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(list_sep)
            })
            .unwrap_or_default(),
        "PICK" => {
            if options.include_pick
                && tc
                    .map(|x| x.flags.iter().any(|f| f == "PICK"))
                    .unwrap_or(false)
            {
                "1".to_string()
            } else {
                String::new()
            }
        }
        "VARIANT_CLASS" => variant.variant_class.so_term().to_string(),
        "MINIMISED" => {
            if variant.minimised {
                "1".to_string()
            } else {
                String::new()
            }
        }
        "SYMBOL_SOURCE" => tc
            .and_then(|x| x.gene_symbol_source.as_deref())
            .unwrap_or("")
            .to_string(),
        "HGNC_ID" => tc
            .and_then(|x| x.hgnc_id.as_deref())
            .unwrap_or("")
            .to_string(),
        "CANONICAL" => {
            if tc.map(|x| x.canonical).unwrap_or(false) {
                "YES".to_string()
            } else {
                String::new()
            }
        }
        "MANE" => {
            let mut sets = Vec::new();
            if tc.and_then(|x| x.mane_select.as_ref()).is_some() {
                sets.push("MANE_Select");
            }
            if tc.and_then(|x| x.mane_plus_clinical.as_ref()).is_some() {
                sets.push("MANE_Plus_Clinical");
            }
            sets.join(list_sep)
        }
        "MANE_SELECT" => tc.and_then(|x| x.mane_select.clone()).unwrap_or_default(),
        "MANE_PLUS_CLINICAL" => tc
            .and_then(|x| x.mane_plus_clinical.clone())
            .unwrap_or_default(),
        "TSL" => tc
            .and_then(|x| x.tsl.map(|t| t.to_string()))
            .unwrap_or_default(),
        "APPRIS" => tc.and_then(|x| x.appris.clone()).unwrap_or_default(),
        "CCDS" => tc.and_then(|x| x.ccds.clone()).unwrap_or_default(),
        "ENSP" => tc
            .and_then(|x| x.protein_id.as_deref())
            .unwrap_or("")
            .to_string(),
        "SWISSPROT" => tc.and_then(|x| x.swissprot.clone()).unwrap_or_default(),
        "TREMBL" => tc.and_then(|x| x.trembl.clone()).unwrap_or_default(),
        "UNIPARC" | "UNIPROT_ISOFORM" | "GENCODE_PRIMARY" | "miRNA" | "HGVSg" | "CHECK_REF" => {
            String::new()
        }
        "RefSeq" => tc.and_then(|x| x.refseq.clone()).unwrap_or_default(),
        "GENE_PHENO" => {
            if options.gene_phenotype && variant.colocated_variants.iter().any(|cv| cv.phenotype) {
                "1".to_string()
            } else {
                String::new()
            }
        }
        "DOMAINS" => tc
            .map(|x| {
                x.domains
                    .iter()
                    .map(|(src, id)| format!("{src}:{id}"))
                    .collect::<Vec<_>>()
                    .join(list_sep)
            })
            .unwrap_or_default(),
        "HGVS_OFFSET" => tc
            .and_then(|x| x.hgvs_offset.map(|o| o.to_string()))
            .unwrap_or_default(),
        "CLIN_SIG" => joined_unique(
            variant
                .colocated_variants
                .iter()
                .flat_map(|cv| cv.clin_sig.iter().cloned()),
            list_sep,
        ),
        "SOMATIC" => {
            if variant.colocated_variants.is_empty() {
                String::new()
            } else if variant.colocated_variants.iter().any(|cv| cv.somatic) {
                "1".to_string()
            } else {
                "0".to_string()
            }
        }
        "PHENO" => {
            if variant.colocated_variants.is_empty() {
                String::new()
            } else if variant.colocated_variants.iter().any(|cv| cv.phenotype) {
                "1".to_string()
            } else {
                "0".to_string()
            }
        }
        "PUBMED" => joined_unique(
            variant
                .colocated_variants
                .iter()
                .flat_map(|cv| cv.pubmed.iter().cloned()),
            list_sep,
        ),
        "SIFT" => tc.and_then(|x| x.sift.clone()).unwrap_or_default(),
        "PolyPhen" => tc.and_then(|x| x.polyphen.clone()).unwrap_or_default(),
        "MOTIF_NAME"
        | "MOTIF_POS"
        | "HIGH_INF_POS"
        | "MOTIF_SCORE_CHANGE"
        | "TRANSCRIPTION_FACTORS" => String::new(),
        "AF" => first_frequency(variant, "AF")
            .map(format_freq)
            .unwrap_or_default(),
        f if f.ends_with("_AF") => {
            // `AFR_AF` -> population `AFR`; `gnomADe_NFE_AF` -> `gnomADe_NFE`; `gnomADe_AF` -> `gnomADe`.
            let pop = &f[..f.len() - 3];
            first_frequency(variant, pop)
                .map(format_freq)
                .unwrap_or_default()
        }
        "MAX_AF" => {
            let mut max_af: Option<f64> = None;
            for cv in &variant.colocated_variants {
                for value in cv.frequencies.values() {
                    if max_af.is_none_or(|m| *value > m) {
                        max_af = Some(*value);
                    }
                }
            }
            max_af.map(format_freq).unwrap_or_default()
        }
        "MAX_AF_POPS" => {
            let mut max_af: Option<f64> = None;
            let mut max_pop = String::new();
            for cv in &variant.colocated_variants {
                for (pop, value) in &cv.frequencies {
                    if max_af.is_none_or(|m| *value > m) {
                        max_af = Some(*value);
                        max_pop = pop.clone();
                    }
                }
            }
            max_pop
        }
        _ => tc
            .and_then(|x| x.plugin_data.get(field))
            .or_else(|| variant.plugin_data.get(field))
            .cloned()
            .unwrap_or_default(),
    }
}

/// `OverlapBP` and `OverlapPC` of a structural variant row: the base pairs of
/// the variant span inside the feature and the percentage of the feature they
/// cover. VEP computes them for every structural variant except breakpoints
/// (ensembl-vep `OutputFactory.pm`) and prints them only when the overlap is positive.
fn sv_overlap(variant: &InputVariant, tc: Option<&TranscriptConsequence>) -> Option<(u64, f64)> {
    if !variant.is_structural || variant.variant_class == VariantClass::Translocation {
        return None;
    }
    tc?.feature_overlap(
        variant.start.min(variant.end),
        variant.start.max(variant.end),
    )
}

/// The keys of a row's `Extra` column in VEP's order: the active flag fields
/// first, then every other key the row carries (the structural variant overlap
/// fields and plugin keys), alphabetically.
pub fn extra_keys<'a>(
    variant: &'a InputVariant,
    tc: Option<&'a TranscriptConsequence>,
    options: &FieldOptions,
) -> Vec<&'a str> {
    extra_keys_with(&flag_fields(options), variant, tc, options)
}

/// [`extra_keys`] given the configuration's [`flag_fields`] already computed.
pub fn extra_keys_with<'a>(
    flagged: &[&'static str],
    variant: &'a InputVariant,
    tc: Option<&'a TranscriptConsequence>,
    options: &FieldOptions,
) -> Vec<&'a str> {
    let mut keys: Vec<&str> = flagged.to_vec();
    let mut extra: BTreeSet<&str> = BTreeSet::new();
    if !options.overlaps && sv_overlap(variant, tc).is_some() {
        extra.insert("OverlapBP");
        extra.insert("OverlapPC");
    }
    for k in variant.plugin_data.keys() {
        if !flagged.contains(&k.as_str()) {
            extra.insert(k.as_str());
        }
    }
    if let Some(tc) = tc {
        for k in tc.plugin_data.keys() {
            if !flagged.contains(&k.as_str()) {
                extra.insert(k.as_str());
            }
        }
    }
    keys.extend(extra);
    keys
}

/// The header lines of the default format, in VEP's order, ending with the
/// column line. `--no_headers` callers write none of them.
pub fn default_format_header(
    options: &FieldOptions,
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
    command_line: &str,
) -> Vec<String> {
    let mut lines = common_header_prefix(cache_dir, cache_info);
    lines.push("## Column descriptions:".to_string());
    for col in DEFAULT_OUTPUT_COLS {
        lines.push(format!(
            "## {} : {}",
            col,
            field_description(col).unwrap_or("?")
        ));
    }
    lines.push("## Extra column keys:".to_string());
    for key in flag_fields(options) {
        lines.push(format!(
            "## {} : {}",
            key,
            field_description(key).unwrap_or("?")
        ));
    }
    lines.push(format!("## VEP command-line: {}", command_line));
    lines.push(format!("#{}\tExtra", DEFAULT_OUTPUT_COLS.join("\t")));
    lines
}

/// The header lines of the tab format: every column, flag fields included, gets
/// a description line, and the column line names them all.
pub fn tab_header(
    options: &FieldOptions,
    plugin_fields: &[String],
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
    command_line: &str,
) -> Vec<String> {
    let cols = tab_columns(options, plugin_fields);
    let mut lines = common_header_prefix(cache_dir, cache_info);
    lines.push("## Column descriptions:".to_string());
    for col in &cols {
        lines.push(format!(
            "## {} : {}",
            col,
            field_description(col).unwrap_or("?")
        ));
    }
    lines.push(format!("## VEP command-line: {}", command_line));
    lines.push(format!("#{}", cols.join("\t")));
    lines
}

/// The tab format's columns: the thirteen defaults, the flag fields, then plugin fields.
pub fn tab_columns(options: &FieldOptions, plugin_fields: &[String]) -> Vec<String> {
    let mut cols: Vec<String> = DEFAULT_OUTPUT_COLS.iter().map(|s| s.to_string()).collect();
    cols.extend(flag_fields(options).into_iter().map(String::from));
    for p in plugin_fields {
        if !cols.iter().any(|c| c == p) {
            cols.push(p.clone());
        }
    }
    cols
}

/// The `##VEP=` meta line of a VCF output: version, API version, time, the
/// cache directory and every cache source version as `key="value"`, sorted.
pub fn vcf_vep_meta_line(
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
) -> String {
    let mut line = format!(
        "##VEP=\"v{}.{}\" API=\"v{}\" time=\"{}\"",
        vep_core::VEP_VERSION,
        vep_core::VEP_SUB_VERSION,
        vep_core::VEP_VERSION,
        timestamp()
    );
    if !cache_dir.is_empty() {
        line.push_str(&format!(" cache=\"{cache_dir}\""));
    }
    if let Some(info) = cache_info {
        for (key, value) in &info.source_versions {
            line.push_str(&format!(" {key}=\"{value}\""));
        }
    }
    line
}

fn common_header_prefix(
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
) -> Vec<String> {
    let mut lines = vec![
        format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v{}.{}",
            vep_core::VEP_VERSION,
            vep_core::VEP_SUB_VERSION
        ),
        format!("## Output produced at {}", timestamp()),
        format!("## Using cache in {}", cache_dir),
        format!(
            "## Using API version {}, DB version ?",
            vep_core::VEP_VERSION
        ),
    ];
    if let Some(info) = cache_info {
        for (key, value) in &info.source_versions {
            lines.push(format!("## {} version {}", key, value));
        }
    }
    lines
}

/// The current UTC time in the form VEP's `Output produced at` line uses.
fn timestamp() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Civil date from days since the epoch (proleptic Gregorian), UTC.
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let (h, m, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    format!("{y:04}-{mo:02}-{d:02} {h:02}:{m:02}:{s:02}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn user_group_fields_come_first_and_once() {
        let opts = FieldOptions {
            symbol: true,
            af_1kg: true,
            af: true,
            ..Default::default()
        };
        let f = flag_fields(&opts);
        assert_eq!(&f[..4], &["IMPACT", "DISTANCE", "STRAND", "FLAGS"]);
        assert_eq!(&f[4..7], &["SYMBOL", "SYMBOL_SOURCE", "HGNC_ID"]);
        assert_eq!(f.iter().filter(|x| **x == "AF").count(), 1);
    }

    #[test]
    fn default_options_emit_only_the_user_group() {
        assert_eq!(
            flag_fields(&FieldOptions::default()),
            vec!["IMPACT", "DISTANCE", "STRAND", "FLAGS"]
        );
        assert!(!FieldOptions::default().any_optional_group());
    }

    #[test]
    fn record_rows_are_transcript_major_allele_minor_intergenic_last() {
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.allele_index = 0;
        let mut b = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        b.allele_index = 1;
        let tc = |id: &str| TranscriptConsequence {
            transcript_id: std::sync::Arc::from(id),
            ..Default::default()
        };
        a.transcript_consequences = vec![tc("ENST2"), tc("ENST1")];
        b.transcript_consequences = vec![tc("ENST1"), tc("ENST2")];
        let mut c = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"C".to_vec());
        c.allele_index = 2;
        c.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let rows = record_rows(&[&a, &b, &c]);
        let order: Vec<(String, usize)> = rows
            .iter()
            .map(|(v, t)| {
                (
                    t.map(|t| t.transcript_id.to_string()).unwrap_or_default(),
                    v.allele_index,
                )
            })
            .collect();
        assert_eq!(
            order,
            vec![
                ("ENST1".to_string(), 0),
                ("ENST1".to_string(), 1),
                ("ENST2".to_string(), 0),
                ("ENST2".to_string(), 1),
                (String::new(), 2),
            ]
        );
    }

    /// The reusable form clears what an earlier record left behind and orders
    /// rows exactly as `record_rows` does.
    #[test]
    fn record_rows_into_reuses_the_vector_and_matches_record_rows() {
        let tc = |id: &str| TranscriptConsequence {
            transcript_id: std::sync::Arc::from(id),
            ..Default::default()
        };
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.transcript_consequences = vec![tc("ENST9"), tc("ENST1")];
        let mut b = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        b.allele_index = 1;
        b.transcript_consequences = vec![tc("ENST1")];
        let mut c = InputVariant::new("21".into(), 300, 300, b"C".to_vec(), b"T".to_vec());
        c.most_severe_consequence = Some(Consequence::IntergenicVariant);

        let mut rows = Vec::new();
        record_rows_into(&[&c], &mut rows);
        assert_eq!(rows.len(), 1);
        record_rows_into(&[&a, &b], &mut rows);
        let key = |rows: &[Row<'_>]| -> Vec<(String, usize)> {
            rows.iter()
                .map(|(v, t)| {
                    (
                        t.map(|t| t.transcript_id.to_string()).unwrap_or_default(),
                        v.allele_index,
                    )
                })
                .collect()
        };
        assert_eq!(key(&rows), key(&record_rows(&[&a, &b])));
        assert_eq!(
            key(&rows),
            vec![
                ("ENST1".to_string(), 0),
                ("ENST1".to_string(), 1),
                ("ENST9".to_string(), 0),
            ]
        );
    }

    #[test]
    fn extra_fields_plan_carries_the_options_flag_fields() {
        let options = FieldOptions {
            symbol: true,
            hgvs: true,
            ..Default::default()
        };
        let plan = ExtraFieldsPlan::new(options.clone());
        assert_eq!(plan.flag_fields, flag_fields(&options));
        assert!(plan.any_optional_group);
        assert!(!ExtraFieldsPlan::new(FieldOptions::default()).any_optional_group);
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            extra_keys_with(&plan.flag_fields, &v, None, &options),
            extra_keys(&v, None, &options)
        );
    }

    #[test]
    fn flags_keep_cds_codes_only() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            flags: std::sync::Arc::from(["cds_start_NF".to_string(), "gencode_basic".to_string()]),
            ..Default::default()
        };
        assert_eq!(
            field_value("FLAGS", &v, Some(&tc), &FieldOptions::default(), ","),
            "cds_start_NF"
        );
    }

    #[test]
    fn header_ends_with_the_column_line() {
        let lines = default_format_header(&FieldOptions::default(), "/c", None, "vep -i x");
        assert!(lines[0].starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR v"));
        assert!(lines.contains(&"## Column descriptions:".to_string()));
        assert!(lines.contains(&"## Extra column keys:".to_string()));
        assert!(lines.contains(&"## FLAGS : Transcript quality flags".to_string()));
        assert_eq!(
            lines.last().unwrap(),
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra"
        );
    }
}
