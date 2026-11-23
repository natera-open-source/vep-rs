// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` source whose short name is `ENSP`, in a run under
//! `--protein_version` and not under `--protein`, on the rows the binary
//! writes for one record: line 42 of the custom corpus's input, a
//! substitution with rows against sixteen transcripts, thirteen of which have
//! a translation, with one record of the corpus's `custom/sites.vcf.gz` on
//! its allele. The key holds the source's values on every row: no flag of the
//! run asks for a protein identifier, so the version flag has none to add a
//! version to.
//!
//! Ensembl VEP 116.3 appends the translation's version to whatever its row
//! hash holds under `ENSP` (OutputFactory.pm:1569-1574), which without
//! `--protein` is the list the source pushed (OutputFactory.pm:1300-1304), and
//! the append turns the list into the text of its memory address: the
//! documented difference `custom_key_named_ensp_under_protein_version`. Each
//! expected row is the row Ensembl VEP 116.3 writes for the same record and
//! flags, with the source's record name in the cells where its row holds that
//! address; every other cell is Ensembl VEP's.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 42;

/// The source's short name, the key of its record names.
const KEY: &str = "ENSP";

/// The name of the record the source's file holds at the input record's
/// position, on the input record's allele.
const SOURCE_RECORD: &str = "multi_7";

/// The corpus's site file as a VCF source named [`KEY`] of the type `kind`
/// (`exact` or `overlap`), the custom corpus's `custom` directory written
/// `{custom}`.
fn source(kind: &str) -> String {
    format!("file={{custom}}/sites.vcf.gz,short_name={KEY},format=vcf,type={kind},fields=AF%CLASS")
}

/// The rows of the default format under an `exact` source and
/// `--protein_version`.
const EXACT_DEFAULT: &[&str] = &[
    "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1037;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2623;STRAND=1;FLAGS=cds_start_NF;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2873;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2586;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=11;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=10;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;ENSP=multi_7;ENSP_AF=0.021;ENSP_CLASS=multi_allelic;ENSP_FILTER=PASS",
];

/// The column line and the rows of the tab format under an `exact` source
/// and `--protein_version`.
const EXACT_TAB: &[&str] = &[
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tENSP\tENSP_AF\tENSP_CLASS",
    "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t1037\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2623\t1\tcds_start_NF\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2873\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2586\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t11\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t10\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
];

/// The `CSQ` format and the record's line of the VCF under an `exact` source
/// and `--protein_version`.
const EXACT_VCF: &[&str] = &[
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|ENSP|ENSP_AF|ENSP_CLASS",
    "21\t33549625\trs2145814066\tC\tT\t.\t.\tCSQ=T|upstream_gene_variant|MODIFIER|MIR6501|ENSG00000284448|Transcript|ENST00000290239|miRNA|||||||||||1037|1||HGNC|HGNC:50033||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000300278|protein_coding|3/7||||1273|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000356577|protein_coding|3/12||||449|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000381679|protein_coding|3/5||||1273|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000381692|protein_coding||2/10||||||||||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|upstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000436227|protein_coding|||||||||||2623|1|cds_start_NF|HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|stop_gained&NMD_transcript_variant|HIGH|SON|ENSG00000159140|Transcript|ENST00000455528|nonsense_mediated_decay|3/13||||423|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000475072|retained_intron|||||||||||2873|1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000492229|protein_coding_CDS_not_defined|||||||||||2586|1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|intron_variant&NMD_transcript_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000695559|nonsense_mediated_decay||2/12||||||||||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000704334|protein_coding||2/10||||||||||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036731|protein_coding||2/5||||||||||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036732|protein_coding||2/3||||||||||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129604|protein_coding|||||||||||11|1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129605|protein_coding|||||||||||10|1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00001131645|protein_coding|3/12||||1273|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183||multi_7|0.021|multi_allelic",
];

/// The column line and the rows of the tab format under an `overlap` source
/// and `--protein_version`: the record's name stands on the record and every
/// row holds it.
const OVERLAP_TAB: &[&str] = &[
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tENSP\tENSP_AF\tENSP_CLASS",
    "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t1037\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2623\t1\tcds_start_NF\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2873\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2586\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t11\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t10\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021,0.521\tmulti_allelic",
];

/// The column line and the row of the tab format under an `exact` source,
/// `--protein_version` and `--pick`.
const PICK_TAB: &[&str] = &[
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tENSP\tENSP_AF\tENSP_CLASS",
    "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tmulti_7\t0.021\tmulti_allelic",
];

/// The key's cell on each of the sixteen rows under an `exact` source,
/// `--protein` and `--protein_version`, as Ensembl VEP 116.3 writes it: the
/// identifier with its version on a row whose transcript has a translation,
/// the source's record name on the other three.
const WITH_PROTEIN_CELLS: &[&str] = &[
    "multi_7",
    "ENSP00000300278.2",
    "ENSP00000348984.4",
    "ENSP00000371095.4",
    "ENSP00000371111.2",
    "ENSP00000400985.1",
    "ENSP00000399783.1",
    "multi_7",
    "multi_7",
    "ENSP00000512016.1",
    "ENSP00000515848.1",
    "ENSP00000706548.1",
    "ENSP00000706549.1",
    "ENSP00000799409.1",
    "ENSP00000799410.1",
    "ENSP00000808470.1",
];

/// The JSON record Ensembl VEP 116.3 writes under an `exact` source and
/// `--protein_version`, its keys sorted, without its copy of the input line.
const EXACT_JSON: &str = r#"{"allele_string":"C/T","assembly_name":"GRCh38","custom_annotations":{"ENSP":[{"allele":"T","fields":{"AF":0.021,"CLASS":"multi_allelic","FILTER":["PASS"]},"name":"multi_7"}]},"end":33549625,"id":"rs2145814066","most_severe_consequence":"stop_gained","seq_region_name":"21","start":33549625,"strand":1,"transcript_consequences":[{"consequence_terms":["upstream_gene_variant"],"distance":1037,"gene_id":"ENSG00000284448","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000290239","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":1273,"cdna_start":1273,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000300278","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":449,"cdna_start":449,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000356577","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":1273,"cdna_start":1273,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000381679","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000381692","variant_allele":"T"},{"consequence_terms":["upstream_gene_variant"],"distance":2623,"flags":["cds_start_NF"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000436227","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":423,"cdna_start":423,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained","NMD_transcript_variant"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000455528","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":2873,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000475072","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":2586,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000492229","variant_allele":"T"},{"consequence_terms":["intron_variant","NMD_transcript_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000695559","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000704334","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001036731","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001036732","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":11,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001129604","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":10,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001129605","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":1273,"cdna_start":1273,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00001131645","variant_allele":"T"}]}"#;

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The record alone, written under `dir`: the header lines of the custom
/// corpus's input and its line [`RECORD_LINE`].
fn record_file(dir: &Path) -> PathBuf {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let record = input
        .lines()
        .nth(RECORD_LINE - 1)
        .expect("the record's line");
    assert!(
        record.starts_with("21\t33549625\t"),
        "line {RECORD_LINE} of the custom corpus's input is the record"
    );
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    text.push_str(record);
    text.push('\n');
    let path = dir.join("record.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// The output of a run over the record in `format` (`default`, `tab`, `vcf`,
/// `json` or `parquet`) under a source of the type `kind` and `flags`, read
/// over the custom corpus's own cache.
fn annotate(dir: &Path, format: &str, kind: &str, flags: &[&str]) -> PathBuf {
    let out = dir.join(format!("{kind}.{format}.{}.out", flags.join("")));
    let custom = corpus("GRCh38-custom").join("custom");
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-custom").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(flags)
        .args([
            "--custom",
            &source(kind).replace("{custom}", &custom.display().to_string()),
        ]);
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
        "vep {format} {kind} {flags:?}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// The rows of a default-format output.
fn default_rows(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

/// The column line and the rows of a tab output.
fn tab_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with("##")).collect()
}

/// The `CSQ` format and the record's line of a VCF output.
fn vcf_lines(text: &str) -> Vec<&str> {
    let format = text
        .lines()
        .find_map(|l| l.strip_prefix("##INFO=<ID=CSQ,"))
        .and_then(|l| l.split_once("Format: "))
        .map(|(_, format)| format.trim_end_matches("\">"))
        .expect("the CSQ header line");
    let record = text
        .lines()
        .find(|l| !l.starts_with('#'))
        .expect("the record's line");
    vec![format, record]
}

/// The key's cell on each row of the tab lines `lines`, whose first line is
/// the column line.
fn key_cells<'a>(lines: &[&'a str]) -> Vec<&'a str> {
    let place = lines[0]
        .trim_start_matches('#')
        .split('\t')
        .position(|name| name == KEY)
        .expect("the key's column");
    lines[1..]
        .iter()
        .map(|row| row.split('\t').nth(place).expect("the key's cell"))
        .collect()
}

/// Under `--protein_version` alone the key holds the source's record name on
/// each of the sixteen rows, in the default, tab and VCF formats; every other
/// cell is the one Ensembl VEP 116.3 writes.
#[test]
fn a_custom_key_named_ensp_keeps_the_sources_values_under_the_protein_version_flag() {
    let tmp = tempfile::tempdir().unwrap();
    let flags = ["--protein_version"];

    let text = fs::read_to_string(annotate(tmp.path(), "default", "exact", &flags)).unwrap();
    assert_eq!(default_rows(&text), EXACT_DEFAULT, "default format");

    let text = fs::read_to_string(annotate(tmp.path(), "tab", "exact", &flags)).unwrap();
    assert_eq!(tab_lines(&text), EXACT_TAB, "tab format");
    assert_eq!(key_cells(EXACT_TAB), [SOURCE_RECORD; 16]);

    let text = fs::read_to_string(annotate(tmp.path(), "vcf", "exact", &flags)).unwrap();
    assert_eq!(vcf_lines(&text), EXACT_VCF, "VCF format");
}

/// A source of the type `overlap` puts its record on the input record, and
/// every row holds the record's name under the key.
#[test]
fn a_whole_record_source_named_ensp_keeps_its_values_on_every_row() {
    let tmp = tempfile::tempdir().unwrap();
    let out = annotate(tmp.path(), "tab", "overlap", &["--protein_version"]);
    let text = fs::read_to_string(out).unwrap();
    assert_eq!(tab_lines(&text), OVERLAP_TAB);
    assert_eq!(key_cells(OVERLAP_TAB), [SOURCE_RECORD; 16]);
}

/// The one row `--pick` keeps is a row of a transcript that has a
/// translation, and its cell holds the source's record name.
#[test]
fn the_row_pick_keeps_holds_the_sources_value_under_the_key() {
    let tmp = tempfile::tempdir().unwrap();
    let flags = ["--protein_version", "--pick"];
    let text = fs::read_to_string(annotate(tmp.path(), "tab", "exact", &flags)).unwrap();
    assert_eq!(tab_lines(&text), PICK_TAB);
    assert_eq!(key_cells(PICK_TAB), [SOURCE_RECORD]);
}

/// With `--protein` the field sets its value on a row whose transcript has a
/// translation, and the version flag adds the version to it: the identifier
/// stands over the source's values there, in both engines, and vep-rs's
/// column line names the key once.
#[test]
fn the_protein_flag_writes_the_identifier_with_its_version_over_the_sources_values() {
    let tmp = tempfile::tempdir().unwrap();
    let flags = ["--protein", "--protein_version"];
    let text = fs::read_to_string(annotate(tmp.path(), "tab", "exact", &flags)).unwrap();
    let lines = tab_lines(&text);
    let listed = lines[0].split('\t').filter(|name| *name == KEY).count();
    assert_eq!(listed, 1, "the key heads one column");
    assert_eq!(key_cells(&lines), WITH_PROTEIN_CELLS);

    let named: Vec<&&str> = WITH_PROTEIN_CELLS
        .iter()
        .filter(|cell| **cell == SOURCE_RECORD)
        .collect();
    assert_eq!(named.len(), 3, "three transcripts have no translation");
}

/// A JSON entry carries no key of the source and, without `--protein`, no
/// `protein_id`: the record is the one Ensembl VEP 116.3 writes.
#[test]
fn a_custom_key_named_ensp_stays_out_of_the_json_entries() {
    let tmp = tempfile::tempdir().unwrap();
    let expected: serde_json::Value = serde_json::from_str(EXACT_JSON).unwrap();
    let out = annotate(tmp.path(), "json", "exact", &["--protein_version"]);
    let text = fs::read_to_string(out).unwrap();
    let mut record: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = record.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    assert_eq!(record, expected);

    let entries = expected["transcript_consequences"]
        .as_array()
        .expect("the transcript entries");
    assert_eq!(entries.len(), 16);
    assert!(entries
        .iter()
        .all(|entry| entry.get("protein_id").is_none()));
}

/// The Parquet output stores the tab cell: the source's record name on each
/// of the sixteen rows.
#[cfg(feature = "parquet")]
#[test]
fn the_parquet_cell_of_a_custom_key_named_ensp_is_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    let place = EXACT_TAB[0]
        .trim_start_matches('#')
        .split('\t')
        .position(|name| name == KEY)
        .expect("the key's column");
    let out = annotate(tmp.path(), "parquet", "exact", &["--protein_version"]);
    let rows = vep_io::output::parquet::tab_rows(&out).unwrap();
    let cells: Vec<&str> = rows
        .iter()
        .map(|row| row.split('\t').nth(place).expect("the key's cell"))
        .collect();
    assert_eq!(cells, [SOURCE_RECORD; 16]);
}
