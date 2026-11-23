// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells one of the thirteen columns of the default
//! format, on the rows the binary writes for two records of the custom
//! corpus's input, for each of whose sites the corpus's `custom/sites.vcf.gz`
//! holds a record on one allele: line 9, whose two rows are intergenic, and
//! line 64, whose eight rows are on transcripts. Ensembl VEP keys a row's hash
//! by the whole string, so the column and the source's key are one cell. The
//! expected rows are the ones Ensembl VEP 116.3 writes for the same record,
//! source and flags; where it lists the key a second time, vep-rs lists it
//! once (the documented difference `duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// A record of the custom corpus's `variants.vcf` and the committed cache a
/// run over it reads.
struct Record {
    /// The 1-based line of the record in the input.
    line: usize,
    /// The chromosome and position the line opens with.
    opens_with: &'static str,
    /// The corpus whose `json_cache` holds the transcripts around the record.
    cache: &'static str,
}

/// Two alternate alleles, each on one intergenic row; the site file's record
/// names the first.
const INTERGENIC_RECORD: Record = Record {
    line: 9,
    opens_with: "21\t5380280\t",
    cache: "GRCh38-everything",
};

/// One alternate allele on eight transcripts of one gene: four rows inside
/// the coding sequence, three inside the 3' UTR and one upstream of the
/// transcript. The site file's record names the allele.
const TRANSCRIPT_RECORD: Record = Record {
    line: 64,
    opens_with: "21\t43683285\t",
    cache: "GRCh38-custom",
};

/// What Ensembl VEP 116.3 writes for a record under one source and one flag
/// set: the rows of the default format, the column line and the rows of the
/// tab format, and the `CSQ` format and the record's line of the VCF. A list
/// is empty where the tests read no output of that format. A record line of
/// [`TRANSCRIPT_RECORD`] is stated from its `CSQ` key on, as
/// `vcf_record_line` reads it.
struct Reference {
    short_name: &'static str,
    /// The key of the source that spells a column.
    key: &'static str,
    /// The source's `type`: an `exact` record names the allele it matched, an
    /// `overlap` record is the whole input record's.
    source_type: &'static str,
    flags: &'static [&'static str],
    default: &'static [&'static str],
    tab: &'static [&'static str],
    vcf: &'static [&'static str],
}

/// `Gene` under an `exact` source of that name, whose record names one allele.
const GENE_EXACT: Reference = Reference {
    short_name: "Gene",
    key: "Gene",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Gene_AF=0.004;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tGene\tGene_AF\tGene_CLASS",
        "rs1444031311\t21:5380280\tC\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Gene|Gene_AF|Gene_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||exact_2||||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Feature` under an `exact` source of that name, whose record names one allele.
const FEATURE_EXACT: Reference = Reference {
    short_name: "Feature",
    key: "Feature",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\texact_2\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Feature_AF=0.004;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tFeature\tFeature_AF\tFeature_CLASS",
        "rs1444031311\t21:5380280\tC\t-\texact_2\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Feature|Feature_AF|Feature_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||exact_2||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Feature_type` under an `exact` source of that name, whose record names one allele.
const FEATURE_TYPE_EXACT: Reference = Reference {
    short_name: "Feature_type",
    key: "Feature_type",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\texact_2\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Feature_type_AF=0.004;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tFeature_type\tFeature_type_AF\tFeature_type_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\texact_2\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Feature_type|Feature_type_AF|Feature_type_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||exact_2|||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `cDNA_position` under an `exact` source of that name, whose record names one allele.
const CDNA_POSITION_EXACT: Reference = Reference {
    short_name: "cDNA_position",
    key: "cDNA_position",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\texact_2\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;cDNA_position_AF=0.004;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tcDNA_position\tcDNA_position_AF\tcDNA_position_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\texact_2\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|cDNA_position|cDNA_position_AF|cDNA_position_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||exact_2||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `CDS_position` under an `exact` source of that name, whose record names one allele.
const CDS_POSITION_EXACT: Reference = Reference {
    short_name: "CDS_position",
    key: "CDS_position",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\texact_2\t-\t-\t-\t-\tIMPACT=MODIFIER;CDS_position_AF=0.004;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCDS_position\tCDS_position_AF\tCDS_position_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\texact_2\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|CDS_position|CDS_position_AF|CDS_position_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||exact_2|||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Protein_position` under an `exact` source of that name, whose record names one allele.
const PROTEIN_POSITION_EXACT: Reference = Reference {
    short_name: "Protein_position",
    key: "Protein_position",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\texact_2\t-\t-\t-\tIMPACT=MODIFIER;Protein_position_AF=0.004;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tProtein_position\tProtein_position_AF\tProtein_position_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\texact_2\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Protein_position|Protein_position_AF|Protein_position_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||exact_2||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Amino_acids` under an `exact` source of that name, whose record names one allele.
const AMINO_ACIDS_EXACT: Reference = Reference {
    short_name: "Amino_acids",
    key: "Amino_acids",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\texact_2\t-\t-\tIMPACT=MODIFIER;Amino_acids_AF=0.004;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tAmino_acids\tAmino_acids_AF\tAmino_acids_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\texact_2\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Amino_acids|Amino_acids_AF|Amino_acids_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||exact_2|||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Codons` under an `exact` source of that name, whose record names one allele.
const CODONS_EXACT: Reference = Reference {
    short_name: "Codons",
    key: "Codons",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\texact_2\t-\tIMPACT=MODIFIER;Codons_AF=0.004;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCodons\tCodons_AF\tCodons_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\texact_2\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Codons|Codons_AF|Codons_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||||exact_2||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Consequence` under an `exact` source of that name: the record's name after the row's term.
const CONSEQUENCE_EXACT: Reference = Reference {
    short_name: "Consequence",
    key: "Consequence",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant,exact_2\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Consequence_AF=0.004;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tConsequence\tConsequence_AF\tConsequence_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant,exact_2\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\tintergenic_variant,exact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\tintergenic_variant\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Consequence|Consequence_AF|Consequence_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant&exact_2|MODIFIER||||||||||||||||||||||intergenic_variant&exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||intergenic_variant||",
    ],
};

/// `Existing_variation` under an `exact` source of that name, without `--check_existing`.
const EXISTING_VARIATION_EXACT: Reference = Reference {
    short_name: "Existing_variation",
    key: "Existing_variation",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\texact_2\tIMPACT=MODIFIER;Existing_variation_AF=0.004;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tExisting_variation\tExisting_variation_AF\tExisting_variation_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\texact_2\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Existing_variation|Existing_variation_AF|Existing_variation_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||exact_2|||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `Existing_variation` under `--check_existing`: the record's name before the co-located identifier.
const EXISTING_VARIATION_EXACT_CHECKED: Reference = Reference {
    short_name: "Existing_variation",
    key: "Existing_variation",
    source_type: "exact",
    flags: &["--check_existing"],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\texact_2,rs1444031311\tIMPACT=MODIFIER;Existing_variation_AF=0.004;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCLIN_SIG\tSOMATIC\tPHENO\tExisting_variation\tExisting_variation_AF\tExisting_variation_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\texact_2,rs1444031311\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\texact_2,rs1444031311\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\trs1444031311\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|CLIN_SIG|SOMATIC|PHENO|Existing_variation|Existing_variation_AF|Existing_variation_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||exact_2&rs1444031311||||||||||exact_2&rs1444031311|0.004|exact,G|intergenic_variant|MODIFIER|||||||||||||||rs1444031311||||||||||rs1444031311||",
    ],
};

/// `Gene` under an `overlap` source of that name, whose record is the whole input record's.
const GENE_OVERLAP: Reference = Reference {
    short_name: "Gene",
    key: "Gene",
    source_type: "overlap",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Gene_AF=0.004;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Gene_AF=0.004;Gene_CLASS=exact;Gene_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tGene\tGene_AF\tGene_CLASS",
        "rs1444031311\t21:5380280\tC\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Gene|Gene_AF|Gene_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||exact_2||||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||exact_2||||||||||||||||||||exact_2|0.004|exact",
    ],
};

/// `Gene` under an `overlap` source on the `--most_severe` row, which is the record's hash.
const GENE_OVERLAP_MOST_SEVERE: Reference = Reference {
    short_name: "Gene",
    key: "Gene",
    source_type: "overlap",
    flags: &["--most_severe"],
    default: &[
        "rs1444031311\t21:5380280\t-\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tGene_AF=0.004;Gene_CLASS=exact;Gene_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tGene\tGene_AF\tGene_CLASS",
        "rs1444031311\t21:5380280\t-\texact_2\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\texact_2\t0.004\texact",
    ],
    vcf: &[],
};

/// `Consequence` under an `overlap` source: the row's term list is assigned over the whole record's values.
const CONSEQUENCE_OVERLAP: Reference = Reference {
    short_name: "Consequence",
    key: "Consequence",
    source_type: "overlap",
    flags: &[],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Consequence_AF=0.004;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;Consequence_AF=0.004;Consequence_CLASS=exact;Consequence_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tConsequence\tConsequence_AF\tConsequence_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\tintergenic_variant\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\tintergenic_variant\t0.004\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Consequence|Consequence_AF|Consequence_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||||||||||||intergenic_variant|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||intergenic_variant|0.004|exact",
    ],
};

/// `Existing_variation` under an `overlap` source, without `--check_existing`.
const EXISTING_VARIATION_OVERLAP: Reference = Reference {
    short_name: "Existing_variation",
    key: "Existing_variation",
    source_type: "overlap",
    flags: &[],
    default: &[],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tExisting_variation\tExisting_variation_AF\tExisting_variation_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\texact_2\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\texact_2\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
    ],
    vcf: &[],
};

/// `Allele` under an `overlap` source: the row's allele is assigned over the whole record's values.
const ALLELE_OVERLAP: Reference = Reference {
    short_name: "Allele",
    key: "Allele",
    source_type: "overlap",
    flags: &[],
    default: &[],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tAllele\tAllele_AF\tAllele_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\tC\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\tG\t0.004\texact",
    ],
    vcf: &[],
};

/// `Allele` under an `overlap` source on the `--most_severe` row, which sets no allele.
const ALLELE_OVERLAP_MOST_SEVERE: Reference = Reference {
    short_name: "Allele",
    key: "Allele",
    source_type: "overlap",
    flags: &["--most_severe"],
    default: &[
        "rs1444031311\t21:5380280\texact_2\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tAllele_AF=0.004;Allele_CLASS=exact;Allele_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tAllele\tAllele_AF\tAllele_CLASS",
        "rs1444031311\t21:5380280\texact_2\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\texact_2\t0.004\texact",
    ],
    vcf: &[],
};

/// The record Ensembl VEP 116.3 writes in JSON under an `exact` source named `Gene`, its `input` key aside.
const GENE_EXACT_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"Gene":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The record Ensembl VEP 116.3 writes in JSON under an `exact` source named `Consequence`, its `input` key aside.
const CONSEQUENCE_EXACT_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"Consequence":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The record Ensembl VEP 116.3 writes in JSON under an `exact` source named `Feature`, its `input` key aside.
const FEATURE_EXACT_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"Feature":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The record Ensembl VEP 116.3 writes in JSON under an `exact` source named `Existing_variation` and `--check_existing`, its `input` key aside.
const EXISTING_VARIATION_EXACT_CHECKED_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","colocated_variants":[{"allele_string":"A/C/G","end":5380280,"frequencies":{"C":{"gnomadg":7.994e-06,"gnomadg_afr":0,"gnomadg_ami":0,"gnomadg_amr":0,"gnomadg_asj":0,"gnomadg_eas":0,"gnomadg_fin":0,"gnomadg_mid":0,"gnomadg_nfe":1.993e-05,"gnomadg_remaining":0,"gnomadg_sas":0}},"id":"rs1444031311","seq_region_name":"21","start":5380280,"strand":1}],"custom_annotations":{"Existing_variation":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The record Ensembl VEP 116.3 writes in JSON under an `overlap` source named `Gene`, its `input` key aside.
const GENE_OVERLAP_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"Gene":[{"fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The record Ensembl VEP 116.3 writes in JSON under an `overlap` source named `Gene` and `--most_severe`, its `input` key aside.
const GENE_OVERLAP_MOST_SEVERE_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"Gene":[{"fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"]}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// `Gene` under an `exact` source of that name, on the eight transcript rows of line 64.
const GENE_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Gene",
    key: "Gene",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Gene_AF=0.034;Gene_CLASS=exact;Gene_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tGene\tGene_AF\tGene_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENSG00000160208\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Gene|Gene_AF|Gene_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||||||||2847|1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||||||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||||||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||||||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENSG00000160208|0.034|exact",
    ],
};

/// `Feature` under an `exact` source of that name, on the eight transcript rows of line 64.
const FEATURE_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Feature",
    key: "Feature",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_AF=0.034;Feature_CLASS=exact;Feature_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tFeature\tFeature_AF\tFeature_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENST00000340648\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\tENST00000470886\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tENST00000989746\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tENST00001043852\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tENST00001043853\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENST00001062422\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENST00001126721\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tENST00001131407\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Feature|Feature_AF|Feature_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENST00000340648|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||||||||2847|1||HGNC|HGNC:23818||ENST00000470886|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||||||1||HGNC|HGNC:23818||ENST00000989746|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||||||1||HGNC|HGNC:23818||ENST00001043852|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||||||1||HGNC|HGNC:23818||ENST00001043853|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENST00001062422|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENST00001126721|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||ENST00001131407|0.034|exact",
    ],
};

/// `Feature_type` under an `exact` source of that name, on the eight transcript rows of line 64.
const FEATURE_TYPE_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Feature_type",
    key: "Feature_type",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Feature_type_AF=0.034;Feature_type_CLASS=exact;Feature_type_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tFeature_type\tFeature_type_AF\tFeature_type_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tTranscript\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tTranscript\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Feature_type|Feature_type_AF|Feature_type_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||Transcript|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||||||||2847|1||HGNC|HGNC:23818||Transcript|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||||||1||HGNC|HGNC:23818||Transcript|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||||||1||HGNC|HGNC:23818||Transcript|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||||||1||HGNC|HGNC:23818||Transcript|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||Transcript|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||Transcript|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||Transcript|0.034|exact",
    ],
};

/// `cDNA_position` under an `exact` source of that name, on the eight transcript rows of line 64.
const CDNA_POSITION_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "cDNA_position",
    key: "cDNA_position",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\texact_12\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;cDNA_position_AF=0.034;cDNA_position_CLASS=exact;cDNA_position_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tcDNA_position\tcDNA_position_AF\tcDNA_position_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t908\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\texact_12\t-\t-\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\t963\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\t868\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\t909\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t882\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t919\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t919\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|cDNA_position|cDNA_position_AF|cDNA_position_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||908|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||exact_12||||||2847|1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||||||1||HGNC|HGNC:23818||963|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||||||1||HGNC|HGNC:23818||868|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||||||1||HGNC|HGNC:23818||909|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||882|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||919|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||919|0.034|exact",
    ],
};

/// `CDS_position` under an `exact` source of that name, on the eight transcript rows of line 64.
const CDS_POSITION_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "CDS_position",
    key: "CDS_position",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\texact_12\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\texact_12\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\texact_12\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\texact_12\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;CDS_position_AF=0.034;CDS_position_CLASS=exact;CDS_position_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCDS_position\tCDS_position_AF\tCDS_position_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t803\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\texact_12\t-\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\texact_12\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\texact_12\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\texact_12\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t803\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t803\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t803\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|CDS_position|CDS_position_AF|CDS_position_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||803|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron||||||exact_12|||||2847|1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|exact_12||||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|exact_12||||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|exact_12||||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||803|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||803|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||803|0.034|exact",
    ],
};

/// `Protein_position` under an `exact` source of that name, on the eight transcript rows of line 64.
const PROTEIN_POSITION_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Protein_position",
    key: "Protein_position",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\texact_12\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\texact_12\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\texact_12\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\texact_12\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Protein_position_AF=0.034;Protein_position_CLASS=exact;Protein_position_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tProtein_position\tProtein_position_AF\tProtein_position_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t268\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\texact_12\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\texact_12\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\texact_12\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\texact_12\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t268\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t268\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\t268\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Protein_position|Protein_position_AF|Protein_position_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||268|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||||exact_12||||2847|1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963||exact_12|||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868||exact_12|||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909||exact_12|||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||268|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||268|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||268|0.034|exact",
    ],
};

/// `Amino_acids` under an `exact` source of that name, on the eight transcript rows of line 64.
const AMINO_ACIDS_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Amino_acids",
    key: "Amino_acids",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\texact_12\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\texact_12\t-\t-\tIMPACT=MODIFIER;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\texact_12\t-\t-\tIMPACT=MODIFIER;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\texact_12\t-\t-\tIMPACT=MODIFIER;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Amino_acids_AF=0.034;Amino_acids_CLASS=exact;Amino_acids_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tAmino_acids\tAmino_acids_AF\tAmino_acids_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tG/D\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\texact_12\t-\t-\tMODIFIER\t2847\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\texact_12\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\texact_12\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\texact_12\t-\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tG/D\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tG/D\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tG/D\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Amino_acids|Amino_acids_AF|Amino_acids_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||G/D|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron||||||||exact_12|||2847|1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||exact_12||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||exact_12||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||exact_12||||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||G/D|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||G/D|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||G/D|0.034|exact",
    ],
};

/// `Codons` under an `exact` source of that name, on the eight transcript rows of line 64.
const CODONS_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Codons",
    key: "Codons",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\texact_12\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\texact_12\t-\tIMPACT=MODIFIER;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\texact_12\t-\tIMPACT=MODIFIER;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\texact_12\t-\tIMPACT=MODIFIER;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Codons_AF=0.034;Codons_CLASS=exact;Codons_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCodons\tCodons_AF\tCodons_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tgGc/gAc\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\texact_12\t-\tMODIFIER\t2847\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\texact_12\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\texact_12\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\texact_12\t-\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tgGc/gAc\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tgGc/gAc\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tgGc/gAc\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Codons|Codons_AF|Codons_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||gGc/gAc|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||||||exact_12||2847|1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963||||exact_12|||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868||||exact_12|||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909||||exact_12|||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||gGc/gAc|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||gGc/gAc|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||gGc/gAc|0.034|exact",
    ],
};

/// `Consequence` under an `exact` source of that name, on the eight transcript rows of line 64.
const CONSEQUENCE_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Consequence",
    key: "Consequence",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant,exact_12\t908\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant,exact_12\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t963\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t868\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t909\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant,exact_12\t882\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant,exact_12\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant,exact_12\t919\t803\t268\tG/D\tgGc/gAc\t-\tIMPACT=MODERATE;STRAND=1;Consequence_AF=0.034;Consequence_CLASS=exact;Consequence_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tConsequence\tConsequence_AF\tConsequence_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant,exact_12\t908\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tmissense_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant,exact_12\t-\t-\t-\t-\t-\t-\tMODIFIER\t2847\t1\t-\t-\tupstream_gene_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t963\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t868\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t909\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\t3_prime_UTR_variant,NMD_transcript_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant,exact_12\t882\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tmissense_variant,NMD_transcript_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant,exact_12\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tmissense_variant,exact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant,exact_12\t919\t803\t268\tG/D\tgGc/gAc\t-\tMODERATE\t-\t1\t-\t-\tmissense_variant,exact_12\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Consequence|Consequence_AF|Consequence_CLASS",
        "CSQ=A|missense_variant&exact_12|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||missense_variant&exact_12|0.034|exact,A|upstream_gene_variant&exact_12|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron|||||||||||2847|1||HGNC|HGNC:23818||upstream_gene_variant&exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant&exact_12|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||||||1||HGNC|HGNC:23818||3_prime_UTR_variant&NMD_transcript_variant&exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant&exact_12|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||||||1||HGNC|HGNC:23818||3_prime_UTR_variant&NMD_transcript_variant&exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant&exact_12|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||||||1||HGNC|HGNC:23818||3_prime_UTR_variant&NMD_transcript_variant&exact_12|0.034|exact,A|missense_variant&NMD_transcript_variant&exact_12|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||missense_variant&NMD_transcript_variant&exact_12|0.034|exact,A|missense_variant&exact_12|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||missense_variant&exact_12|0.034|exact,A|missense_variant&exact_12|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|||1||HGNC|HGNC:23818||missense_variant&exact_12|0.034|exact",
    ],
};

/// `Existing_variation` under an `exact` source of that name, on the eight transcript rows of line 64.
const EXISTING_VARIATION_EXACT_ON_TRANSCRIPTS: Reference = Reference {
    short_name: "Existing_variation",
    key: "Existing_variation",
    source_type: "exact",
    flags: &[],
    default: &[
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\texact_12\tIMPACT=MODERATE;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\texact_12\tIMPACT=MODIFIER;DISTANCE=2847;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\texact_12\tIMPACT=MODIFIER;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\texact_12\tIMPACT=MODIFIER;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\texact_12\tIMPACT=MODIFIER;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\texact_12\tIMPACT=MODERATE;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\texact_12\tIMPACT=MODERATE;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\texact_12\tIMPACT=MODERATE;STRAND=1;Existing_variation_AF=0.034;Existing_variation_CLASS=exact;Existing_variation_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tExisting_variation\tExisting_variation_AF\tExisting_variation_CLASS",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000340648\tTranscript\tmissense_variant\t908\t803\t268\tG/D\tgGc/gAc\texact_12\tMODERATE\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000470886\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\texact_12\tMODIFIER\t2847\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00000989746\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t963\t-\t-\t-\t-\texact_12\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043852\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t868\t-\t-\t-\t-\texact_12\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001043853\tTranscript\t3_prime_UTR_variant,NMD_transcript_variant\t909\t-\t-\t-\t-\texact_12\tMODIFIER\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001062422\tTranscript\tmissense_variant,NMD_transcript_variant\t882\t803\t268\tG/D\tgGc/gAc\texact_12\tMODERATE\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001126721\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\texact_12\tMODERATE\t-\t1\t-\t-\texact_12\t0.034\texact",
        "rs755653178\t21:43683285\tA\tENSG00000160208\tENST00001131407\tTranscript\tmissense_variant\t919\t803\t268\tG/D\tgGc/gAc\texact_12\tMODERATE\t-\t1\t-\t-\texact_12\t0.034\texact",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|Existing_variation|Existing_variation_AF|Existing_variation_CLASS",
        "CSQ=A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00000340648|protein_coding|9/16||||908|803|268|G/D|gGc/gAc|exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|upstream_gene_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000470886|retained_intron||||||||||exact_12|2847|1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00000989746|nonsense_mediated_decay|9/13||||963|||||exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043852|nonsense_mediated_decay|8/15||||868|||||exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|3_prime_UTR_variant&NMD_transcript_variant|MODIFIER|RRP1B|ENSG00000160208|Transcript|ENST00001043853|nonsense_mediated_decay|9/13||||909|||||exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant&NMD_transcript_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001062422|nonsense_mediated_decay|9/16||||882|803|268|G/D|gGc/gAc|exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001126721|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact,A|missense_variant|MODERATE|RRP1B|ENSG00000160208|Transcript|ENST00001131407|protein_coding|9/16||||919|803|268|G/D|gGc/gAc|exact_12||1||HGNC|HGNC:23818||exact_12|0.034|exact",
    ],
};

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Gene`, its `input` key aside.
const GENE_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Gene":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Feature`, its `input` key aside.
const FEATURE_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Feature":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Feature_type`, its `input` key aside.
const FEATURE_TYPE_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Feature_type":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `cDNA_position`, its `input` key aside.
const CDNA_POSITION_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"cDNA_position":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `CDS_position`, its `input` key aside.
const CDS_POSITION_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"CDS_position":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Protein_position`, its `input` key aside.
const PROTEIN_POSITION_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Protein_position":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Amino_acids`, its `input` key aside.
const AMINO_ACIDS_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Amino_acids":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Codons`, its `input` key aside.
const CODONS_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Codons":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Consequence`, its `input` key aside.
const CONSEQUENCE_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Consequence":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The record Ensembl VEP 116.3 writes in JSON for line 64 under an `exact` source named `Existing_variation`, its `input` key aside.
const EXISTING_VARIATION_EXACT_ON_TRANSCRIPTS_JSON: &str = r#"{"allele_string":"G/A","assembly_name":"GRCh38","custom_annotations":{"Existing_variation":[{"allele":"A","fields":{"AF":0.034,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_12"}]},"end":43683285,"id":"rs755653178","most_severe_consequence":"missense_variant","seq_region_name":"21","start":43683285,"strand":1,"transcript_consequences":[{"amino_acids":"G/D","cdna_end":908,"cdna_start":908,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00000340648","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":2847,"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000470886","variant_allele":"A"},{"cdna_end":963,"cdna_start":963,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000989746","variant_allele":"A"},{"cdna_end":868,"cdna_start":868,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043852","variant_allele":"A"},{"cdna_end":909,"cdna_start":909,"consequence_terms":["3_prime_UTR_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001043853","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":882,"cdna_start":882,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant","NMD_transcript_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001062422","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001126721","variant_allele":"A"},{"amino_acids":"G/D","cdna_end":919,"cdna_start":919,"cds_end":803,"cds_start":803,"codons":"gGc/gAc","consequence_terms":["missense_variant"],"gene_id":"ENSG00000160208","impact":"MODERATE","protein_end":268,"protein_start":268,"strand":1,"transcript_id":"ENST00001131407","variant_allele":"A"}]}"#;

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The custom corpus's input and the line of `record` in it.
fn input_and_line(record: &Record) -> (String, String) {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let line = input
        .lines()
        .nth(record.line - 1)
        .expect("the record's line")
        .to_string();
    assert!(
        line.starts_with(record.opens_with),
        "line {} of the custom corpus's input is the record",
        record.line
    );
    (input, line)
}

/// The record alone, written under `dir`: the header lines of the custom
/// corpus's input and the record's line.
fn record_file(dir: &Path, record: &Record) -> PathBuf {
    let (input, line) = input_and_line(record);
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    text.push_str(&line);
    text.push('\n');
    let path = dir.join(format!("record.{}.vcf", record.line));
    fs::write(&path, text).unwrap();
    path
}

/// The output of a run over `record` in `format` (`default`, `tab`, `vcf`,
/// `json` or `parquet`) under the flags `flags`, with the custom corpus's site
/// file as a source of the type `source_type` under the name `short_name`
/// that asks for its `AF` and `CLASS` fields.
fn annotate(
    dir: &Path,
    record: &Record,
    format: &str,
    short_name: &str,
    source_type: &str,
    flags: &[&str],
) -> PathBuf {
    let out = dir.join(format!(
        "{}.{short_name}.{source_type}.{}.{format}.out",
        record.line,
        flags.len()
    ));
    let source = format!(
        "file={},short_name={short_name},format=vcf,type={source_type},fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display(),
    );
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir, record))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus(record.cache).join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(["--custom", &source])
        .args(flags);
    match format {
        "default" => {}
        "parquet" => {
            cmd.args(["--output_format", "parquet"]);
        }
        other => {
            cmd.arg(format!("--{other}"));
        }
    }
    let run = cmd.output().unwrap();
    assert!(
        run.status.success(),
        "vep {format} {short_name} {source_type} {flags:?}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// `line` cut at `separator` without the cell at `place`, counted from 0.
fn without_cell(line: &str, separator: char, place: usize) -> String {
    let cells: Vec<&str> = line
        .split(separator)
        .enumerate()
        .filter(|(i, _)| *i != place)
        .map(|(_, cell)| cell)
        .collect();
    cells.join(&separator.to_string())
}

/// The place, counted from 0, at which `names` lists `key` a second time.
fn second_listing(names: &str, separator: char, key: &str) -> Option<usize> {
    names
        .trim_start_matches('#')
        .split(separator)
        .enumerate()
        .filter(|(_, name)| *name == key)
        .map(|(i, _)| i)
        .nth(1)
}

/// The reference's tab lines as vep-rs writes them: the key's second column,
/// which Ensembl VEP lists after the flag columns, left out of the column line
/// and of each row.
fn tab_listed_once(reference: &Reference) -> Vec<String> {
    let place = second_listing(reference.tab[0], '\t', reference.key)
        .expect("Ensembl VEP lists the key twice");
    reference
        .tab
        .iter()
        .map(|line| without_cell(line, '\t', place))
        .collect()
}

/// The reference's `CSQ` format and record line as vep-rs writes them: the
/// key's second field left out of the format and of each allele's entry.
fn vcf_listed_once(reference: &Reference) -> Vec<String> {
    let (format, record) = (reference.vcf[0], reference.vcf[1]);
    let place =
        second_listing(format, '|', reference.key).expect("Ensembl VEP lists the key twice");
    let (columns, entries) = record.split_once("CSQ=").expect("the CSQ key");
    let entries: Vec<String> = entries
        .split(',')
        .map(|entry| without_cell(entry, '|', place))
        .collect();
    vec![
        without_cell(format, '|', place),
        format!("{columns}CSQ={}", entries.join(",")),
    ]
}

/// The whole line of `record` in a VCF output, from the line a reference
/// states. Ensembl VEP writes the `CSQ` key in place of an INFO column that
/// holds `.` and after a semicolon behind one that holds values, so a stated
/// line that opens with the key follows the input line.
fn vcf_record_line(record: &Record, stated: &str) -> String {
    if stated.starts_with("CSQ=") {
        format!("{};{stated}", input_and_line(record).1)
    } else {
        stated.to_string()
    }
}

/// Holds the binary's default, tab and VCF outputs over `record` for the
/// reference's source and flags to the reference's rows, in each format the
/// reference states.
fn assert_rows_are_the_references(dir: &Path, record: &Record, reference: &Reference) {
    let case = format!(
        "line {}, {} as an {} source, {:?}",
        record.line, reference.key, reference.source_type, reference.flags
    );
    let run = |format: &str| {
        let out = annotate(
            dir,
            record,
            format,
            reference.short_name,
            reference.source_type,
            reference.flags,
        );
        fs::read_to_string(out).unwrap()
    };
    if !reference.default.is_empty() {
        let text = run("default");
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(rows, reference.default, "default format, {case}");
    }
    if !reference.tab.is_empty() {
        let text = run("tab");
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with("##")).collect();
        assert_eq!(rows, tab_listed_once(reference), "tab format, {case}");
    }
    if !reference.vcf.is_empty() {
        let text = run("vcf");
        let format = text
            .lines()
            .find_map(|l| l.strip_prefix("##INFO=<ID=CSQ,"))
            .and_then(|l| l.split_once("Format: "))
            .map(|(_, format)| format.trim_end_matches("\">"))
            .expect("the CSQ header line");
        let line = text
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("the record's line");
        let listed_once = vcf_listed_once(reference);
        assert_eq!(
            vec![format, line],
            vec![
                listed_once[0].as_str(),
                vcf_record_line(record, &listed_once[1]).as_str()
            ],
            "VCF format, {case}"
        );
    }
}

/// Holds the record the binary writes in JSON over `record` for the
/// reference's source and flags to `json`, the record Ensembl VEP writes: as
/// parsed values, without the `input` key and vep-rs's own `vep_rs_version`.
fn assert_json_record_is_the_references(
    dir: &Path,
    record: &Record,
    reference: &Reference,
    json: &str,
) {
    let expected: serde_json::Value = serde_json::from_str(json).unwrap();
    let out = annotate(
        dir,
        record,
        "json",
        reference.short_name,
        reference.source_type,
        reference.flags,
    );
    let text = fs::read_to_string(out).unwrap();
    let mut written: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = written.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    assert_eq!(
        written, expected,
        "line {}, {} as an {} source, {:?}",
        record.line, reference.key, reference.source_type, reference.flags
    );
}

/// A source's values are pushed under its key (OutputFactory.pm:1300-1304)
/// and a feature assigns its own columns afterwards (`Gene` at :1477,
/// `Feature` at :1473, `Feature_type` at :1472, the positions at :1735, :1746
/// and :1750, `Amino_acids` and `Codons` at :1742-1743), so a row that names
/// no feature keeps the source's values in the column. An `exact` record
/// prints on the row of the allele it names; a whole record's values stand in
/// the record's hash (:1004-1014) and so on every row and on the
/// `--most_severe` row. The default format lists a column's name in no `Extra`
/// item (VEP_output.pm:104-109).
#[test]
fn a_custom_key_named_as_a_column_a_feature_sets_prints_the_sources_value_on_a_row_without_one() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [
        &GENE_EXACT,
        &FEATURE_EXACT,
        &FEATURE_TYPE_EXACT,
        &CDNA_POSITION_EXACT,
        &CDS_POSITION_EXACT,
        &PROTEIN_POSITION_EXACT,
        &AMINO_ACIDS_EXACT,
        &CODONS_EXACT,
        &GENE_OVERLAP,
        &GENE_OVERLAP_MOST_SEVERE,
    ] {
        assert_rows_are_the_references(tmp.path(), &INTERGENIC_RECORD, reference);
    }
}

/// A row's term list is assigned (OutputFactory.pm:1341) and the records of
/// the row's allele are pushed onto it afterwards (:1426-1435), so an `exact`
/// record's name follows the row's terms; the assignment drops the values of a
/// whole record, which stood in the record's hash before it.
#[test]
fn a_custom_key_named_consequence_lists_the_allele_records_after_the_rows_terms() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [&CONSEQUENCE_EXACT, &CONSEQUENCE_OVERLAP] {
        assert_rows_are_the_references(tmp.path(), &INTERGENIC_RECORD, reference);
    }
}

/// The co-located identifiers are pushed under `Existing_variation` after the
/// source's values (OutputFactory.pm:1086, called at :1438), so the source's
/// values come first under `--check_existing` and stand alone without it.
#[test]
fn a_custom_key_named_existing_variation_lists_the_sources_values_before_the_identifiers() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [
        &EXISTING_VARIATION_EXACT,
        &EXISTING_VARIATION_EXACT_CHECKED,
        &EXISTING_VARIATION_OVERLAP,
    ] {
        assert_rows_are_the_references(tmp.path(), &INTERGENIC_RECORD, reference);
    }
}

/// A row's allele is assigned over the values of a whole record
/// (OutputFactory.pm:1347); the `--most_severe` row is the record's hash,
/// which sets no allele (:426-458), so its `Allele` column holds the source's
/// values.
#[test]
fn a_custom_key_named_allele_prints_the_rows_allele_and_the_sources_value_on_a_summary_row() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [&ALLELE_OVERLAP, &ALLELE_OVERLAP_MOST_SEVERE] {
        assert_rows_are_the_references(tmp.path(), &INTERGENIC_RECORD, reference);
    }
}

/// The JSON factory moves the custom annotations to the record before any row
/// hash is built (OutputFactory/JSON.pm:249-258), so the JSON output prints
/// the source's record under `custom_annotations` and each entry's keys hold
/// the column's own value alone.
#[test]
fn a_custom_key_named_as_a_default_column_stays_out_of_the_json_entries() {
    let tmp = tempfile::tempdir().unwrap();
    for (reference, json) in [
        (&GENE_EXACT, GENE_EXACT_JSON),
        (&CONSEQUENCE_EXACT, CONSEQUENCE_EXACT_JSON),
        (&FEATURE_EXACT, FEATURE_EXACT_JSON),
        (
            &EXISTING_VARIATION_EXACT_CHECKED,
            EXISTING_VARIATION_EXACT_CHECKED_JSON,
        ),
        (&GENE_OVERLAP, GENE_OVERLAP_JSON),
        (&GENE_OVERLAP_MOST_SEVERE, GENE_OVERLAP_MOST_SEVERE_JSON),
    ] {
        assert_json_record_is_the_references(tmp.path(), &INTERGENIC_RECORD, reference, json);
    }
}

/// A transcript assigns `Feature_type`, `Feature` and `Gene` on each of its
/// rows (OutputFactory.pm:1472-1473, :1477), `cDNA_position` on a row inside
/// an exon (:1735), and `Amino_acids`, `Codons`, `CDS_position` and
/// `Protein_position` on a row inside the coding sequence (:1742-1743, :1746,
/// :1750), each after the source's values are pushed under the key
/// (:1300-1304, called at :1426-1435). A row therefore holds the feature's
/// own value in a column its feature sets and the source's value in a column
/// it does not: on line 64 the upstream row in the five position, residue and
/// codon columns, and a 3' UTR row in the four of them that a coding row
/// alone sets.
#[test]
fn a_custom_key_named_as_a_column_a_feature_sets_prints_the_features_value_on_a_row_that_sets_it() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [
        &GENE_EXACT_ON_TRANSCRIPTS,
        &FEATURE_EXACT_ON_TRANSCRIPTS,
        &FEATURE_TYPE_EXACT_ON_TRANSCRIPTS,
        &CDNA_POSITION_EXACT_ON_TRANSCRIPTS,
        &CDS_POSITION_EXACT_ON_TRANSCRIPTS,
        &PROTEIN_POSITION_EXACT_ON_TRANSCRIPTS,
        &AMINO_ACIDS_EXACT_ON_TRANSCRIPTS,
        &CODONS_EXACT_ON_TRANSCRIPTS,
    ] {
        assert_rows_are_the_references(tmp.path(), &TRANSCRIPT_RECORD, reference);
    }
}

/// The term list of a transcript's row is assigned (OutputFactory.pm:1341)
/// before the records of the row's allele are pushed under the key
/// (:1426-1435), as on an intergenic row, so each of the eight rows lists the
/// record's name after its own terms.
#[test]
fn a_custom_key_named_consequence_lists_the_allele_records_after_each_transcript_rows_terms() {
    let tmp = tempfile::tempdir().unwrap();
    assert_rows_are_the_references(
        tmp.path(),
        &TRANSCRIPT_RECORD,
        &CONSEQUENCE_EXACT_ON_TRANSCRIPTS,
    );
}

/// No feature assigns `Existing_variation`, and a record without co-located
/// variants pushes no identifier under it (OutputFactory.pm:1046, :1086), so
/// each row of a transcript holds the source's value.
#[test]
fn a_custom_key_named_existing_variation_prints_the_sources_value_on_each_transcript_row() {
    let tmp = tempfile::tempdir().unwrap();
    assert_rows_are_the_references(
        tmp.path(),
        &TRANSCRIPT_RECORD,
        &EXISTING_VARIATION_EXACT_ON_TRANSCRIPTS,
    );
}

/// The JSON factory takes the custom annotations off the record before a
/// transcript's row hash is built (OutputFactory/JSON.pm:249-258), so each
/// entry holds the keys its transcript sets and no value of the source: on a
/// row whose feature sets the column, on a row whose feature does not, and in
/// `consequence_terms`, which lists the row's terms alone.
#[test]
fn a_custom_key_named_as_a_default_column_stays_out_of_the_json_entries_of_a_transcript() {
    let tmp = tempfile::tempdir().unwrap();
    for (reference, json) in [
        (&GENE_EXACT_ON_TRANSCRIPTS, GENE_EXACT_ON_TRANSCRIPTS_JSON),
        (
            &FEATURE_EXACT_ON_TRANSCRIPTS,
            FEATURE_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &FEATURE_TYPE_EXACT_ON_TRANSCRIPTS,
            FEATURE_TYPE_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &CDNA_POSITION_EXACT_ON_TRANSCRIPTS,
            CDNA_POSITION_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &CDS_POSITION_EXACT_ON_TRANSCRIPTS,
            CDS_POSITION_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &PROTEIN_POSITION_EXACT_ON_TRANSCRIPTS,
            PROTEIN_POSITION_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &AMINO_ACIDS_EXACT_ON_TRANSCRIPTS,
            AMINO_ACIDS_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &CODONS_EXACT_ON_TRANSCRIPTS,
            CODONS_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &CONSEQUENCE_EXACT_ON_TRANSCRIPTS,
            CONSEQUENCE_EXACT_ON_TRANSCRIPTS_JSON,
        ),
        (
            &EXISTING_VARIATION_EXACT_ON_TRANSCRIPTS,
            EXISTING_VARIATION_EXACT_ON_TRANSCRIPTS_JSON,
        ),
    ] {
        assert_json_record_is_the_references(tmp.path(), &TRANSCRIPT_RECORD, reference, json);
    }
}

/// The helpers that drop Ensembl VEP's second listing take one cell out of the
/// column line and of each row, and one field out of the `CSQ` format and of
/// each allele's entry.
#[test]
fn the_second_listing_is_dropped_from_the_column_line_and_from_each_row() {
    let tab = tab_listed_once(&GENE_EXACT);
    assert_eq!(GENE_EXACT.tab[0].split('\t').count(), 21);
    assert_eq!(tab[0].split('\t').filter(|c| *c == "Gene").count(), 1);
    assert_eq!(tab[0].split('\t').nth(3), Some("Gene"));
    assert_eq!(tab[1].split('\t').nth(3), Some("exact_2"));
    assert_eq!(tab[2].split('\t').nth(3), Some("-"));
    for line in &tab {
        assert_eq!(line.split('\t').count(), 20);
    }

    let vcf = vcf_listed_once(&CONSEQUENCE_EXACT);
    assert_eq!(vcf[0].split('|').filter(|f| *f == "Consequence").count(), 1);
    assert_eq!(vcf[0].split('|').count(), 26);
    let entries = vcf[1].split_once("CSQ=").unwrap().1;
    for (entry, terms) in entries
        .split(',')
        .zip(["intergenic_variant&exact_2", "intergenic_variant"])
    {
        assert_eq!(entry.split('|').count(), 26);
        assert_eq!(entry.split('|').nth(1), Some(terms));
    }
}

/// A record line stated from its `CSQ` key on follows the eight columns of the
/// input line after a semicolon, and a line stated whole stands as stated.
#[test]
fn a_record_line_stated_from_its_csq_key_follows_the_input_line() {
    let line = vcf_record_line(&TRANSCRIPT_RECORD, "CSQ=A|term");
    let columns: Vec<&str> = line.split('\t').collect();
    assert_eq!(columns.len(), 8);
    assert_eq!(columns[..2], ["21", "43683285"]);
    assert!(columns[7].ends_with(";CSQ=A|term"));
    assert_eq!(columns[7].matches("CSQ=").count(), 1);

    let whole = GENE_EXACT.vcf[1];
    assert_eq!(vcf_record_line(&INTERGENIC_RECORD, whole), whole);
}

/// The Parquet output stores the column as the `--tab` cell: the cell of each
/// row under the key's column is the reference's.
#[cfg(feature = "parquet")]
#[test]
fn a_custom_key_named_as_a_default_column_is_stored_in_parquet_as_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [
        &GENE_EXACT,
        &CONSEQUENCE_EXACT,
        &EXISTING_VARIATION_EXACT_CHECKED,
        &GENE_OVERLAP,
    ] {
        let tab = tab_listed_once(reference);
        let key = tab[0]
            .trim_start_matches('#')
            .split('\t')
            .position(|column| column == reference.key)
            .expect("the key's column");
        let cell = |row: &str| {
            row.split('\t')
                .nth(key)
                .expect("the key's cell")
                .to_string()
        };
        let mut expected: Vec<String> = tab[1..].iter().map(|row| cell(row)).collect();
        expected.sort_unstable();
        let out = annotate(
            tmp.path(),
            &INTERGENIC_RECORD,
            "parquet",
            reference.short_name,
            reference.source_type,
            reference.flags,
        );
        let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
            .unwrap()
            .iter()
            .map(|row| cell(row))
            .collect();
        cells.sort_unstable();
        assert_eq!(
            cells, expected,
            "{} as an {} source, {:?}",
            reference.key, reference.source_type, reference.flags
        );
    }
}
