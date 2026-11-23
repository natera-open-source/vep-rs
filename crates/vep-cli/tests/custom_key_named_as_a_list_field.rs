// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a flag field Ensembl VEP holds as a list, on
//! the rows the binary writes for four records. The first is line 42 of the
//! custom corpus's input, a substitution with rows against sixteen
//! transcripts, one of them in the MANE Select set, with one co-located
//! variant that carries synonyms and one record of the corpus's
//! `custom/sites.vcf.gz` on its allele. The second is line 111 of the
//! zygosity corpus's input, which carries eight genotype columns, and the
//! third line 188 of the statistics corpus's input, whose co-located variant
//! holds a frequency of 0 in nine populations and none above it; the fourth
//! is line 9 of the custom corpus's input, whose co-located variant holds a
//! frequency above 0 for one allele. Each of the last three is read beside a
//! made-up source of one record on its allele. Ensembl VEP pushes a source's
//! values onto its key (OutputFactory.pm:1300-1304) and pushes the items of
//! such a field onto the same key, so the cell holds both. The expected rows
//! are the ones Ensembl VEP 116.3 writes for the same record and flags, with
//! three documented differences: where it lists the key a second time,
//! vep-rs lists it once (`duplicate_output_column_collapsed`); where it lists
//! an allele's record once for each row, vep-rs lists it once
//! (`custom_exact_source_allele_record_accumulation`, whose reproduction
//! record's further run one test makes); and where it lists the samples'
//! zygosities in a hash's order, vep-rs lists them in the header's
//! (`individual_zyg_list_order`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// A record of a committed corpus's input and the cache its rows are read
/// over.
struct Record {
    /// The corpus whose `variants.vcf` holds the record.
    corpus: &'static str,
    /// The 1-based line of the record in that file.
    line: usize,
    /// The record's chromosome and position, as its line opens.
    opens: &'static str,
    /// The corpus whose `json_cache` the run reads.
    cache: &'static str,
}

/// The record with a transcript in the MANE Select set and a co-located
/// variant that carries synonyms.
const WITH_SYNONYMS: Record = Record {
    corpus: "GRCh38-custom",
    line: 42,
    opens: "21\t33549625\t",
    cache: "GRCh38-everything",
};

/// The record with eight genotype columns.
const WITH_GENOTYPES: Record = Record {
    corpus: "GRCh38-individual-zyg",
    line: 111,
    opens: "21\t42951004\t",
    cache: "GRCh38-individual-zyg",
};

/// The record whose co-located variant holds a frequency of 0 in nine
/// populations and none above it.
const WITH_A_ZERO_MAXIMUM: Record = Record {
    corpus: "GRCh38-statistics",
    line: 188,
    opens: "21\t44285966\t",
    cache: "GRCh38-statistics",
};

/// The record whose co-located variant holds a frequency above 0 for the
/// first of its two alternate alleles, in one population.
const WITH_A_POSITIVE_MAXIMUM: Record = Record {
    corpus: "GRCh38-custom",
    line: 9,
    opens: "21\t5380280\t",
    cache: "GRCh38-everything",
};

/// What Ensembl VEP 116.3 writes for the record under one set of sources and
/// flags: the rows of the default format, the column line and the rows of the
/// tab format, and the `CSQ` format and the record's line of the VCF. A
/// format the reference was not run in is empty.
struct Reference {
    /// The record the run annotates.
    record: &'static Record,
    /// The `--custom` values of the run, the custom corpus's `custom`
    /// directory written `{custom}` and the directory of the made-up sources
    /// `{made_up}`.
    sources: &'static [&'static str],
    /// The key of the sources that spells a flag field.
    key: &'static str,
    flags: &'static [&'static str],
    default: &'static [&'static str],
    tab: &'static [&'static str],
    vcf: &'static [&'static str],
}

/// `MANE` under `--mane_select` and an `exact` source of that name: the record's name,
/// then the set, on the row of the transcript that is in the set.
const MANE_SELECT: Reference = Reference {
    record: &WITH_SYNONYMS,
    sources: &["file={custom}/sites.vcf.gz,short_name=MANE,format=vcf,type=exact,fields=AF%CLASS"],
    key: "MANE",
    flags: &["--mane_select"],
    default: &[
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1037;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;MANE=multi_7,MANE_Select;MANE_SELECT=NM_138927.4;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2623;STRAND=1;FLAGS=cds_start_NF;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2873;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2586;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=11;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=10;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tIMPACT=HIGH;STRAND=1;MANE=multi_7;MANE_AF=0.021;MANE_CLASS=multi_allelic;MANE_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tMANE\tMANE_SELECT\tSOURCE\tMANE\tMANE_AF\tMANE_CLASS",
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t1037\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\tmulti_7,MANE_Select\tNM_138927.4\t-\tmulti_7,MANE_Select\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2623\t1\tcds_start_NF\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2873\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2586\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t11\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t10\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\tmulti_7\t-\t-\tmulti_7\t0.021\tmulti_allelic",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|MANE|MANE_SELECT|SOURCE|MANE|MANE_AF|MANE_CLASS",
        "21\t33549625\trs2145814066\tC\tT\t.\t.\tCSQ=T|upstream_gene_variant|MODIFIER|MIR6501|ENSG00000284448|Transcript|ENST00000290239|miRNA|||||||||||1037|1||HGNC|HGNC:50033|multi_7|||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000300278|protein_coding|3/7||||1273|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000356577|protein_coding|3/12||||449|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183|multi_7&MANE_Select|NM_138927.4||multi_7&MANE_Select|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000381679|protein_coding|3/5||||1273|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000381692|protein_coding||2/10||||||||||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|upstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000436227|protein_coding|||||||||||2623|1|cds_start_NF|HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|stop_gained&NMD_transcript_variant|HIGH|SON|ENSG00000159140|Transcript|ENST00000455528|nonsense_mediated_decay|3/13||||423|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000475072|retained_intron|||||||||||2873|1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000492229|protein_coding_CDS_not_defined|||||||||||2586|1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|intron_variant&NMD_transcript_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000695559|nonsense_mediated_decay||2/12||||||||||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000704334|protein_coding||2/10||||||||||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036731|protein_coding||2/5||||||||||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036732|protein_coding||2/3||||||||||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129604|protein_coding|||||||||||11|1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129605|protein_coding|||||||||||10|1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00001131645|protein_coding|3/12||||1273|394|132|Q/*|Cag/Tag|||1||HGNC|HGNC:11183|multi_7|||multi_7|0.021|multi_allelic",
    ],
};

/// `VAR_SYNONYMS` under `--check_existing --var_synonyms` and an `exact` source of that
/// name: the record's name and the co-located variant's synonyms, joined by two hyphens.
const VAR_SYNONYMS_ON: Reference = Reference {
    record: &WITH_SYNONYMS,
    sources: &["file={custom}/sites.vcf.gz,short_name=VAR_SYNONYMS,format=vcf,type=exact,fields=AF%CLASS"],
    key: "VAR_SYNONYMS",
    flags: &["--check_existing", "--var_synonyms"],
    default: &[
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=1037;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=2623;STRAND=1;FLAGS=cds_start_NF;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=2873;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=2586;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=11;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=10;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--ClinVar::RCV002273474,VCV001699617;VAR_SYNONYMS_AF=0.021;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCLIN_SIG\tSOMATIC\tPHENO\tVAR_SYNONYMS\tVAR_SYNONYMS\tVAR_SYNONYMS_AF\tVAR_SYNONYMS_CLASS",
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t1037\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t2623\t1\tcds_start_NF\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t2873\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t2586\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t11\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t10\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--ClinVar::RCV002273474,VCV001699617\tmulti_7--ClinVar::RCV002273474,VCV001699617\t0.021\tmulti_allelic",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|CLIN_SIG|SOMATIC|PHENO|VAR_SYNONYMS|VAR_SYNONYMS|VAR_SYNONYMS_AF|VAR_SYNONYMS_CLASS",
        "21\t33549625\trs2145814066\tC\tT\t.\t.\tCSQ=T|upstream_gene_variant|MODIFIER|MIR6501|ENSG00000284448|Transcript|ENST00000290239|miRNA||||||||||rs2145814066|1037|1||HGNC|HGNC:50033||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000300278|protein_coding|3/7||||1273|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000356577|protein_coding|3/12||||449|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000381679|protein_coding|3/5||||1273|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000381692|protein_coding||2/10||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|upstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000436227|protein_coding||||||||||rs2145814066|2623|1|cds_start_NF|HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|stop_gained&NMD_transcript_variant|HIGH|SON|ENSG00000159140|Transcript|ENST00000455528|nonsense_mediated_decay|3/13||||423|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000475072|retained_intron||||||||||rs2145814066|2873|1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000492229|protein_coding_CDS_not_defined||||||||||rs2145814066|2586|1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|intron_variant&NMD_transcript_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000695559|nonsense_mediated_decay||2/12||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000704334|protein_coding||2/10||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036731|protein_coding||2/5||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036732|protein_coding||2/3||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129604|protein_coding||||||||||rs2145814066|11|1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129605|protein_coding||||||||||rs2145814066|10|1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00001131645|protein_coding|3/12||||1273|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--ClinVar::RCV002273474&VCV001699617|multi_7--ClinVar::RCV002273474&VCV001699617|0.021|multi_allelic",
    ],
};

/// `VAR_SYNONYMS` under `--check_existing`, without `--var_synonyms`, and two `overlap`
/// sources of that name: their two record names, joined by two hyphens.
const VAR_SYNONYMS_OF_TWO_SOURCES: Reference = Reference {
    record: &WITH_SYNONYMS,
    sources: &["file={custom}/sites.vcf.gz,short_name=VAR_SYNONYMS,format=vcf,type=overlap,fields=AF%CLASS", "file={custom}/regions.bed.gz,short_name=VAR_SYNONYMS,format=bed,type=overlap"],
    key: "VAR_SYNONYMS",
    flags: &["--check_existing"],
    default: &[
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=1037;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=2623;STRAND=1;FLAGS=cds_start_NF;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=2873;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=2586;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=11;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tIMPACT=MODIFIER;DISTANCE=10;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tIMPACT=HIGH;STRAND=1;CLIN_SIG=pathogenic;PHENO=1;VAR_SYNONYMS=multi_7--multi_7;VAR_SYNONYMS_AF=0.021,0.521;VAR_SYNONYMS_CLASS=multi_allelic;VAR_SYNONYMS_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tCLIN_SIG\tSOMATIC\tPHENO\tVAR_SYNONYMS\tVAR_SYNONYMS_AF\tVAR_SYNONYMS_CLASS",
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t1037\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t2623\t1\tcds_start_NF\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t2873\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t2586\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t11\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\trs2145814066\tMODIFIER\t10\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\trs2145814066\tHIGH\t-\t1\t-\t-\tpathogenic\t-\t1\tmulti_7--multi_7\t0.021,0.521\tmulti_allelic",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|CLIN_SIG|SOMATIC|PHENO|VAR_SYNONYMS|VAR_SYNONYMS_AF|VAR_SYNONYMS_CLASS",
        "21\t33549625\trs2145814066\tC\tT\t.\t.\tCSQ=T|upstream_gene_variant|MODIFIER|MIR6501|ENSG00000284448|Transcript|ENST00000290239|miRNA||||||||||rs2145814066|1037|1||HGNC|HGNC:50033||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000300278|protein_coding|3/7||||1273|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000356577|protein_coding|3/12||||449|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00000381679|protein_coding|3/5||||1273|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000381692|protein_coding||2/10||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|upstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000436227|protein_coding||||||||||rs2145814066|2623|1|cds_start_NF|HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|stop_gained&NMD_transcript_variant|HIGH|SON|ENSG00000159140|Transcript|ENST00000455528|nonsense_mediated_decay|3/13||||423|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000475072|retained_intron||||||||||rs2145814066|2873|1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000492229|protein_coding_CDS_not_defined||||||||||rs2145814066|2586|1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|intron_variant&NMD_transcript_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000695559|nonsense_mediated_decay||2/12||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00000704334|protein_coding||2/10||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036731|protein_coding||2/5||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|intron_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001036732|protein_coding||2/3||||||||rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129604|protein_coding||||||||||rs2145814066|11|1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|downstream_gene_variant|MODIFIER|SON|ENSG00000159140|Transcript|ENST00001129605|protein_coding||||||||||rs2145814066|10|1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic,T|stop_gained|HIGH|SON|ENSG00000159140|Transcript|ENST00001131645|protein_coding|3/12||||1273|394|132|Q/*|Cag/Tag|rs2145814066||1||HGNC|HGNC:11183||pathogenic||1|multi_7--multi_7|0.021&0.521|multi_allelic",
    ],
};

/// `NEAREST` under `--nearest symbol` and an `exact` source of that name, as Ensembl VEP
/// writes it: the nearest gene's symbol, then the record's name once for each of the
/// sixteen rows of its allele, on every row.
const NEAREST_SYMBOL: Reference = Reference {
    record: &WITH_SYNONYMS,
    sources: &["file={custom}/sites.vcf.gz,short_name=NEAREST,format=vcf,type=exact,fields=AF%CLASS"],
    key: "NEAREST",
    flags: &["--nearest", "symbol"],
    default: &[],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tNEAREST\tNEAREST\tNEAREST_AF\tNEAREST_CLASS",
        "rs2145814066\t21:33549625\tT\tENSG00000284448\tENST00000290239\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t1037\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000300278\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000356577\tTranscript\tstop_gained\t449\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381679\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000381692\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000436227\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2623\t1\tcds_start_NF\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000455528\tTranscript\tstop_gained,NMD_transcript_variant\t423\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000475072\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2873\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000492229\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t2586\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000695559\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00000704334\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036731\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001036732\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129604\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t11\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001129605\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t10\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
        "rs2145814066\t21:33549625\tT\tENSG00000159140\tENST00001131645\tTranscript\tstop_gained\t1273\t394\t132\tQ/*\tCag/Tag\t-\tHIGH\t-\t1\t-\t-\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\tSON,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7,multi_7\t0.021\tmulti_allelic",
    ],
    vcf: &[],
};

/// The record Ensembl VEP 116.3 writes in JSON under the source named `MANE` and
/// `--mane_select`, its `input` key aside.
const MANE_SELECT_JSON: &str = r#"{"allele_string":"C/T","assembly_name":"GRCh38","custom_annotations":{"MANE":[{"allele":"T","fields":{"AF":0.021,"CLASS":"multi_allelic","FILTER":["PASS"]},"name":"multi_7"}]},"end":33549625,"id":"rs2145814066","most_severe_consequence":"stop_gained","seq_region_name":"21","start":33549625,"strand":1,"transcript_consequences":[{"consequence_terms":["upstream_gene_variant"],"distance":1037,"gene_id":"ENSG00000284448","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000290239","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":1273,"cdna_start":1273,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000300278","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":449,"cdna_start":449,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","mane":["MANE_Select"],"mane_select":"NM_138927.4","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000356577","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":1273,"cdna_start":1273,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000381679","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000381692","variant_allele":"T"},{"consequence_terms":["upstream_gene_variant"],"distance":2623,"flags":["cds_start_NF"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000436227","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":423,"cdna_start":423,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained","NMD_transcript_variant"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00000455528","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":2873,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000475072","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":2586,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000492229","variant_allele":"T"},{"consequence_terms":["intron_variant","NMD_transcript_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000695559","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000704334","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001036731","variant_allele":"T"},{"consequence_terms":["intron_variant"],"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001036732","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":11,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001129604","variant_allele":"T"},{"consequence_terms":["downstream_gene_variant"],"distance":10,"gene_id":"ENSG00000159140","impact":"MODIFIER","strand":1,"transcript_id":"ENST00001129605","variant_allele":"T"},{"amino_acids":"Q/*","cdna_end":1273,"cdna_start":1273,"cds_end":394,"cds_start":394,"codons":"Cag/Tag","consequence_terms":["stop_gained"],"gene_id":"ENSG00000159140","impact":"HIGH","protein_end":132,"protein_start":132,"strand":1,"transcript_id":"ENST00001131645","variant_allele":"T"}]}"#;

/// `ZYG` under `--individual_zyg all` and an `exact` source of that name, as Ensembl VEP
/// writes it: the eight samples' zygosities in the order of one run's hash, then the
/// record's name once for each of the two rows of its allele, on both rows.
const ZYG_OF_EVERY_SAMPLE: Reference = Reference {
    record: &WITH_GENOTYPES,
    sources: &["file={made_up}/genotypes_site.vcf.gz,short_name=ZYG,format=vcf,type=exact,fields=AF"],
    key: "ZYG",
    flags: &["--individual_zyg", "all"],
    default: &[
        "21:42951004:T:C\t21:42951004\tC\tENSG00000264580\tENST00000579137\tTranscript\tmature_miRNA_variant\t11\t-\t-\t-\t-\t-\tZYG=S7:HET,S4:HOM,S2:HOM,S3:HOM,S8:HOM,S6:HOM,S5:HET,S1:HET,made_up_1,made_up_1;IMPACT=MODIFIER;STRAND=-1;ZYG_AF=0.25;ZYG_FILTER=PASS",
        "21:42951004:T:C\t21:42951004\tC\tENSG00000296063\tENST00000735996\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tZYG=S7:HET,S4:HOM,S2:HOM,S3:HOM,S8:HOM,S6:HOM,S5:HET,S1:HET,made_up_1,made_up_1;IMPACT=MODIFIER;DISTANCE=4349;STRAND=-1;ZYG_AF=0.25;ZYG_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tZYG\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tZYG\tZYG_AF",
        "21:42951004:T:C\t21:42951004\tC\tENSG00000264580\tENST00000579137\tTranscript\tmature_miRNA_variant\t11\t-\t-\t-\t-\t-\tS7:HET,S4:HOM,S2:HOM,S3:HOM,S8:HOM,S6:HOM,S5:HET,S1:HET,made_up_1,made_up_1\tMODIFIER\t-\t-1\t-\t-\tS7:HET,S4:HOM,S2:HOM,S3:HOM,S8:HOM,S6:HOM,S5:HET,S1:HET,made_up_1,made_up_1\t0.25",
        "21:42951004:T:C\t21:42951004\tC\tENSG00000296063\tENST00000735996\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tS7:HET,S4:HOM,S2:HOM,S3:HOM,S8:HOM,S6:HOM,S5:HET,S1:HET,made_up_1,made_up_1\tMODIFIER\t4349\t-1\t-\t-\tS7:HET,S4:HOM,S2:HOM,S3:HOM,S8:HOM,S6:HOM,S5:HET,S1:HET,made_up_1,made_up_1\t0.25",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ZYG|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|ZYG|ZYG_AF",
        "21\t42951004\t21:42951004:T:C\tT\tC\t.\t.\tAC=5588;AF=0.87258;CM=67.1501;AN=6404;AN_EAS=1170;AN_AMR=980;AN_EUR=1266;AN_AFR=1786;AN_SAS=1202;AN_EUR_unrel=1006;AN_EAS_unrel=1008;AN_AMR_unrel=694;AN_SAS_unrel=978;AN_AFR_unrel=1322;AF_EAS=0.876068;AF_AMR=0.780612;AF_EUR=0.865719;AF_AFR=0.93729;AF_SAS=0.855241;AF_EUR_unrel=0.868787;MAF_EUR_unrel=0.131213;AF_EAS_unrel=0.873016;MAF_EAS_unrel=0.126984;AF_AMR_unrel=0.783862;MAF_AMR_unrel=0.216138;AF_SAS_unrel=0.858896;MAF_SAS_unrel=0.141104;AF_AFR_unrel=0.939486;MAF_AFR_unrel=0.0605144;AC_EAS=1025;AC_AMR=765;AC_EUR=1096;AC_AFR=1674;AC_SAS=1028;AC_EUR_unrel=874;AC_EAS_unrel=880;AC_AMR_unrel=544;AC_SAS_unrel=840;AC_AFR_unrel=1242;AC_Het_EAS=125;AC_Het_AMR=159;AC_Het_EUR=146;AC_Het_AFR=108;AC_Het_SAS=144;AC_Het_EUR_unrel=122;AC_Het_EAS_unrel=112;AC_Het_AMR_unrel=112;AC_Het_SAS_unrel=112;AC_Het_AFR_unrel=76;AC_Het=682;AC_Hom_EAS=900;AC_Hom_AMR=606;AC_Hom_EUR=950;AC_Hom_AFR=1566;AC_Hom_SAS=884;AC_Hom_EUR_unrel=752;AC_Hom_EAS_unrel=768;AC_Hom_AMR_unrel=432;AC_Hom_SAS_unrel=728;AC_Hom_AFR_unrel=1166;AC_Hom=4906;HWE_EAS=0.70239;ExcHet_EAS=0.733052;HWE_AMR=0.236842;ExcHet_AMR=0.90531;HWE_EUR=0.863829;ExcHet_EUR=0.661913;HWE_AFR=0.570694;ExcHet_AFR=0.306469;HWE_SAS=0.410951;ExcHet_SAS=0.836983;HWE=0.0208391;ExcHet=0.99202;CSQ=C|mature_miRNA_variant|MODIFIER|MIR5692B|ENSG00000264580|Transcript|ENST00000579137|miRNA|1/1||||11||||||S7:HET&S4:HOM&S2:HOM&S3:HOM&S8:HOM&S6:HOM&S5:HET&S1:HET&made_up_1&made_up_1||-1||HGNC|HGNC:43535||S7:HET&S4:HOM&S2:HOM&S3:HOM&S8:HOM&S6:HOM&S5:HET&S1:HET&made_up_1&made_up_1|0.25,C|downstream_gene_variant|MODIFIER||ENSG00000296063|Transcript|ENST00000735996|lncRNA|||||||||||S7:HET&S4:HOM&S2:HOM&S3:HOM&S8:HOM&S6:HOM&S5:HET&S1:HET&made_up_1&made_up_1|4349|-1|||||S7:HET&S4:HOM&S2:HOM&S3:HOM&S8:HOM&S6:HOM&S5:HET&S1:HET&made_up_1&made_up_1|0.25\tGT\t1|0\t1|1\t1|1\t1|1\t1|0\t1|1\t1|0\t1|1",
    ],
};

/// The record Ensembl VEP 116.3 writes in JSON under the source named `ZYG` and
/// `--individual_zyg all`, its `input` key aside.
const ZYG_OF_EVERY_SAMPLE_JSON: &str = r#"{"allele_string":"T/C","assembly_name":"GRCh38","custom_annotations":{"ZYG":[{"allele":"C","fields":{"AF":0.25,"FILTER":["PASS"]},"name":"made_up_1"}]},"end":42951004,"id":"21:42951004:T:C","most_severe_consequence":"mature_miRNA_variant","seq_region_name":"21","start":42951004,"strand":1,"transcript_consequences":[{"cdna_end":11,"cdna_start":11,"consequence_terms":["mature_miRNA_variant"],"gene_id":"ENSG00000264580","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000579137","variant_allele":"C"},{"consequence_terms":["downstream_gene_variant"],"distance":4349,"gene_id":"ENSG00000296063","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000735996","variant_allele":"C"}],"zyg":["S7:HET","S4:HOM","S2:HOM","S3:HOM","S8:HOM","S6:HOM","S5:HET","S1:HET"]}"#;

/// `ZYG` under `--individual_zyg S3` and the same source, as Ensembl VEP writes it: the one
/// sample's zygosity, then the record's name twice.
const ZYG_OF_ONE_SAMPLE: Reference = Reference {
    record: &WITH_GENOTYPES,
    sources: &["file={made_up}/genotypes_site.vcf.gz,short_name=ZYG,format=vcf,type=exact,fields=AF"],
    key: "ZYG",
    flags: &["--individual_zyg", "S3"],
    default: &[],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tZYG\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tZYG\tZYG_AF",
        "21:42951004:T:C\t21:42951004\tC\tENSG00000264580\tENST00000579137\tTranscript\tmature_miRNA_variant\t11\t-\t-\t-\t-\t-\tS3:HOM,made_up_1,made_up_1\tMODIFIER\t-\t-1\t-\t-\tS3:HOM,made_up_1,made_up_1\t0.25",
        "21:42951004:T:C\t21:42951004\tC\tENSG00000296063\tENST00000735996\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tS3:HOM,made_up_1,made_up_1\tMODIFIER\t4349\t-1\t-\t-\tS3:HOM,made_up_1,made_up_1\t0.25",
    ],
    vcf: &[],
};

/// `MAX_AF_POPS` under `--max_af` and an `exact` source of that name, on a record whose
/// co-located variant holds a frequency of 0 in nine populations and none above it: the
/// record's name, then the nine populations.
const MAX_AF_POPS_OF_A_ZERO_MAXIMUM: Reference = Reference {
    record: &WITH_A_ZERO_MAXIMUM,
    sources: &["file={made_up}/zero_maximum_site.vcf.gz,short_name=MAX_AF_POPS,format=vcf,type=exact,fields=AF"],
    key: "MAX_AF_POPS",
    flags: &["--max_af"],
    default: &[
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000291582\tTranscript\t5_prime_UTR_variant\t91\t-\t-\t-\t-\trs2146374912\tIMPACT=MODIFIER;STRAND=1;MAX_AF=0;MAX_AF_POPS=made_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS;CLIN_SIG=benign;PHENO=1;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000337909\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2146374912\tIMPACT=MODIFIER;DISTANCE=4789;STRAND=1;MAX_AF=0;MAX_AF_POPS=made_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS;CLIN_SIG=benign;PHENO=1;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000397994\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2146374912\tIMPACT=MODIFIER;DISTANCE=4789;STRAND=1;MAX_AF=0;MAX_AF_POPS=made_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS;CLIN_SIG=benign;PHENO=1;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000527919\tTranscript\tnon_coding_transcript_exon_variant\t121\t-\t-\t-\t-\trs2146374912\tIMPACT=MODIFIER;STRAND=1;MAX_AF=0;MAX_AF_POPS=made_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS;CLIN_SIG=benign;PHENO=1;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000530812\tTranscript\tnon_coding_transcript_exon_variant\t129\t-\t-\t-\t-\trs2146374912\tIMPACT=MODIFIER;STRAND=1;MAX_AF=0;MAX_AF_POPS=made_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS;CLIN_SIG=benign;PHENO=1;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000966178\tTranscript\t5_prime_UTR_variant\t91\t-\t-\t-\t-\trs2146374912\tIMPACT=MODIFIER;STRAND=1;MAX_AF=0;MAX_AF_POPS=made_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS;CLIN_SIG=benign;PHENO=1;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tMAX_AF_POPS\tMAX_AF_POPS_AF",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000291582\tTranscript\t5_prime_UTR_variant\t91\t-\t-\t-\t-\trs2146374912\tMODIFIER\t-\t1\t-\t-\t0\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\tbenign\t-\t1\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\t0.5",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000337909\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2146374912\tMODIFIER\t4789\t1\t-\t-\t0\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\tbenign\t-\t1\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\t0.5",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000397994\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\trs2146374912\tMODIFIER\t4789\t1\t-\t-\t0\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\tbenign\t-\t1\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\t0.5",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000527919\tTranscript\tnon_coding_transcript_exon_variant\t121\t-\t-\t-\t-\trs2146374912\tMODIFIER\t-\t1\t-\t-\t0\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\tbenign\t-\t1\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\t0.5",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000530812\tTranscript\tnon_coding_transcript_exon_variant\t129\t-\t-\t-\t-\trs2146374912\tMODIFIER\t-\t1\t-\t-\t0\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\tbenign\t-\t1\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\t0.5",
        "1292092\t21:44285966\tA\tENSG00000160224\tENST00000966178\tTranscript\t5_prime_UTR_variant\t91\t-\t-\t-\t-\trs2146374912\tMODIFIER\t-\t1\t-\t-\t0\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\tbenign\t-\t1\tmade_up_2,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS\t0.5",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|MAX_AF_POPS|MAX_AF_POPS_AF",
        "21\t44285966\t1292092\tC\tA\t.\t.\tCSQ=A|5_prime_UTR_variant|MODIFIER|AIRE|ENSG00000160224|Transcript|ENST00000291582|protein_coding|1/14||||91|||||rs2146374912||1||HGNC|HGNC:360||0|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|benign||1|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|0.5,A|upstream_gene_variant|MODIFIER|AIRE|ENSG00000160224|Transcript|ENST00000337909|retained_intron||||||||||rs2146374912|4789|1||HGNC|HGNC:360||0|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|benign||1|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|0.5,A|upstream_gene_variant|MODIFIER|AIRE|ENSG00000160224|Transcript|ENST00000397994|retained_intron||||||||||rs2146374912|4789|1||HGNC|HGNC:360||0|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|benign||1|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|0.5,A|non_coding_transcript_exon_variant|MODIFIER|AIRE|ENSG00000160224|Transcript|ENST00000527919|retained_intron|1/14||||121|||||rs2146374912||1||HGNC|HGNC:360||0|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|benign||1|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|0.5,A|non_coding_transcript_exon_variant|MODIFIER|AIRE|ENSG00000160224|Transcript|ENST00000530812|retained_intron|1/12||||129|||||rs2146374912||1||HGNC|HGNC:360||0|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|benign||1|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|0.5,A|5_prime_UTR_variant|MODIFIER|AIRE|ENSG00000160224|Transcript|ENST00000966178|protein_coding|1/14||||91|||||rs2146374912||1||HGNC|HGNC:360||0|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|benign||1|made_up_2&gnomADe_AFR&gnomADe_AMR&gnomADe_ASJ&gnomADe_EAS&gnomADe_FIN&gnomADe_MID&gnomADe_NFE&gnomADe_REMAINING&gnomADe_SAS|0.5",
    ],
};

/// The record Ensembl VEP 116.3 writes in JSON under the source named `MAX_AF_POPS` and
/// `--max_af`, its `input` key aside.
const MAX_AF_POPS_OF_A_ZERO_MAXIMUM_JSON: &str = r#"{"allele_string":"C/A","assembly_name":"GRCh38","colocated_variants":[{"allele_string":"C/A","clin_sig":["benign"],"clin_sig_allele":"A:benign","clin_sig_ref_allele":"C","end":44285966,"frequencies":{"A":{"gnomade":0,"gnomade_afr":0,"gnomade_amr":0,"gnomade_asj":0,"gnomade_eas":0,"gnomade_fin":0,"gnomade_mid":0,"gnomade_nfe":0,"gnomade_remaining":0,"gnomade_sas":0}},"id":"rs2146374912","phenotype_or_disease":1,"seq_region_name":"21","start":44285966,"strand":1,"var_synonyms":{"ClinVar":["RCV001714268","VCV001292092"]}}],"custom_annotations":{"MAX_AF_POPS":[{"allele":"A","fields":{"AF":0.5,"FILTER":["PASS"]},"name":"made_up_2"}]},"end":44285966,"id":"1292092","most_severe_consequence":"5_prime_UTR_variant","seq_region_name":"21","start":44285966,"strand":1,"transcript_consequences":[{"cdna_end":91,"cdna_start":91,"consequence_terms":["5_prime_UTR_variant"],"gene_id":"ENSG00000160224","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000291582","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":4789,"gene_id":"ENSG00000160224","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000337909","variant_allele":"A"},{"consequence_terms":["upstream_gene_variant"],"distance":4789,"gene_id":"ENSG00000160224","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000397994","variant_allele":"A"},{"cdna_end":121,"cdna_start":121,"consequence_terms":["non_coding_transcript_exon_variant"],"gene_id":"ENSG00000160224","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000527919","variant_allele":"A"},{"cdna_end":129,"cdna_start":129,"consequence_terms":["non_coding_transcript_exon_variant"],"gene_id":"ENSG00000160224","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000530812","variant_allele":"A"},{"cdna_end":91,"cdna_start":91,"consequence_terms":["5_prime_UTR_variant"],"gene_id":"ENSG00000160224","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000966178","variant_allele":"A"}]}"#;

/// `MAX_AF_POPS` under `--max_af` and an `exact` source of that name, on a record whose
/// co-located variant holds a frequency above 0 for the first of its two alternate alleles
/// in one population: that population alone on the allele's row, without the record's name.
const MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM: Reference = Reference {
    record: &WITH_A_POSITIVE_MAXIMUM,
    sources: &["file={made_up}/positive_maximum_site.vcf.gz,short_name=MAX_AF_POPS,format=vcf,type=exact,fields=AF"],
    key: "MAX_AF_POPS",
    flags: &["--max_af"],
    default: &[
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;MAX_AF=1.993e-05;MAX_AF_POPS=gnomADg_NFE;MAX_AF_POPS_AF=0.5;MAX_AF_POPS_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER",
    ],
    tab: &[
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tMAX_AF_POPS\tMAX_AF_POPS_AF",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t1.993e-05\tgnomADg_NFE\t-\t-\t-\tgnomADg_NFE\t0.5",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: &[
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|MAX_AF_POPS|MAX_AF_POPS_AF",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||1.993e-05|gnomADg_NFE||||gnomADg_NFE|0.5,G|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||||||||",
    ],
};

/// The record Ensembl VEP 116.3 writes in JSON under the same source and flag, its
/// `input` key aside.
const MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","colocated_variants":[{"allele_string":"A/C/G","end":5380280,"frequencies":{"C":{"gnomadg":7.994e-06,"gnomadg_afr":0,"gnomadg_ami":0,"gnomadg_amr":0,"gnomadg_asj":0,"gnomadg_eas":0,"gnomadg_fin":0,"gnomadg_mid":0,"gnomadg_nfe":1.993e-05,"gnomadg_remaining":0,"gnomadg_sas":0}},"id":"rs1444031311","seq_region_name":"21","start":5380280,"strand":1}],"custom_annotations":{"MAX_AF_POPS":[{"allele":"C","fields":{"AF":0.5,"FILTER":["PASS"]},"name":"made_up_3"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The record alone, written under `dir`: the header lines of its corpus's
/// input and the record's own line.
fn record_file(dir: &Path, record: &Record) -> PathBuf {
    let input = fs::read_to_string(corpus(record.corpus).join("variants.vcf")).unwrap();
    let line = input
        .lines()
        .nth(record.line - 1)
        .expect("the record's line");
    assert!(
        line.starts_with(record.opens),
        "line {} of the input of {} is the record",
        record.line,
        record.corpus
    );
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    text.push_str(line);
    text.push('\n');
    let path = dir.join(format!("record.{}.vcf", record.line));
    fs::write(&path, text).unwrap();
    path
}

/// The three made-up sources, written under `dir` with their indexes: one
/// record on the allele of the record with genotype columns, one on the
/// allele of the record with a maximum of 0 and one on the first alternate
/// allele of the record with a maximum above 0.
fn made_up_sources(dir: &Path) {
    let frequency = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"A made-up allele frequency\">";
    common::indexed_vcf(
        dir,
        "genotypes_site.vcf.gz",
        &["##fileformat=VCFv4.2", frequency, "##contig=<ID=21>"],
        "21\t42951004\tmade_up_1\tT\tC\t.\tPASS\tAF=0.25",
    );
    common::indexed_vcf(
        dir,
        "zero_maximum_site.vcf.gz",
        &[
            "##fileformat=VCFv4.2",
            frequency,
            "##INFO=<ID=POPS,Number=1,Type=String,Description=\"A made-up population name\">",
            "##contig=<ID=21>",
        ],
        "21\t44285966\tmade_up_2\tC\tA\t.\tPASS\tAF=0.5;POPS=made_up_population",
    );
    common::indexed_vcf(
        dir,
        "positive_maximum_site.vcf.gz",
        &[
            "##fileformat=VCFv4.2",
            frequency,
            "##INFO=<ID=POPS,Number=1,Type=String,Description=\"A made-up population name\">",
            "##contig=<ID=21>",
        ],
        "21\t5380280\tmade_up_3\tA\tC\t.\tPASS\tAF=0.5;POPS=made_up_population",
    );
}

/// The output of a run over the reference's record in `format` (`default`,
/// `tab`, `vcf`, `json` or `parquet`) under the reference's sources and
/// flags, read over the cache the record names.
fn annotate(dir: &Path, format: &str, reference: &Reference) -> PathBuf {
    let out = dir.join(format!(
        "{}.{format}.{}.{}.out",
        reference.key,
        reference.sources.len(),
        reference.flags.join("")
    ));
    let custom = corpus("GRCh38-custom").join("custom");
    made_up_sources(dir);
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir, reference.record))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus(reference.record.cache).join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(reference.flags);
    for source in reference.sources {
        cmd.args([
            "--custom",
            &source
                .replace("{custom}", &custom.display().to_string())
                .replace("{made_up}", &dir.display().to_string()),
        ]);
    }
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
        "vep {format} {} {:?}:\n{}",
        reference.key,
        reference.flags,
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

/// The places, counted from 0, at which `names` lists `key`.
fn listings(names: &str, separator: char, key: &str) -> Vec<usize> {
    names
        .trim_start_matches('#')
        .split(separator)
        .enumerate()
        .filter(|(_, name)| *name == key)
        .map(|(i, _)| i)
        .collect()
}

/// The reference's tab lines as vep-rs writes them: the key's second column,
/// where Ensembl VEP lists one, left out of the column line and of each row.
fn tab_listed_once(reference: &Reference) -> Vec<String> {
    match listings(reference.tab[0], '\t', reference.key).get(1) {
        Some(&place) => reference
            .tab
            .iter()
            .map(|line| without_cell(line, '\t', place))
            .collect(),
        None => reference.tab.iter().map(|line| line.to_string()).collect(),
    }
}

/// The reference's `CSQ` format and record line as vep-rs writes them: the
/// key's second field, where Ensembl VEP lists one, left out of the format and
/// of each row's entry.
fn vcf_listed_once(reference: &Reference) -> Vec<String> {
    let (format, record) = (reference.vcf[0], reference.vcf[1]);
    let Some(&place) = listings(format, '|', reference.key).get(1) else {
        return vec![format.to_string(), record.to_string()];
    };
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

/// Holds the default, tab and VCF outputs of a run under the reference's
/// sources and flags to the reference's rows, its second listing of the key
/// dropped.
fn assert_text_formats(dir: &Path, reference: &Reference) {
    let case = format!("{} {:?}", reference.key, reference.flags);

    let text = fs::read_to_string(annotate(dir, "default", reference)).unwrap();
    assert_eq!(
        default_rows(&text),
        reference.default,
        "default format, {case}"
    );

    let text = fs::read_to_string(annotate(dir, "tab", reference)).unwrap();
    assert_eq!(
        tab_lines(&text),
        tab_listed_once(reference),
        "tab format, {case}"
    );

    let text = fs::read_to_string(annotate(dir, "vcf", reference)).unwrap();
    assert_eq!(
        vcf_lines(&text),
        vcf_listed_once(reference),
        "VCF format, {case}"
    );
}

/// The key's cell on each row of the tab lines `lines`, the column line first.
fn key_cells<'a>(lines: &'a [String], key: &str) -> Vec<&'a str> {
    let place = listings(&lines[0], '\t', key)[0];
    lines[1..]
        .iter()
        .map(|line| line.split('\t').nth(place).expect("the key's cell"))
        .collect()
}

/// The column line and the rows of the tab output the binary writes under the
/// reference's sources and flags.
fn written_tab_lines(dir: &Path, reference: &Reference) -> Vec<String> {
    let text = fs::read_to_string(annotate(dir, "tab", reference)).unwrap();
    tab_lines(&text)
        .iter()
        .map(|line| line.to_string())
        .collect()
}

/// The record of the JSON output the binary writes under the reference's
/// sources and flags, without its `input` key and without the key that states
/// the writer's version.
fn json_record(dir: &Path, reference: &Reference) -> serde_json::Value {
    let text = fs::read_to_string(annotate(dir, "json", reference)).unwrap();
    let mut record: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = record.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    record
}

/// Ensembl VEP pushes the sets a transcript is in onto `MANE`
/// (OutputFactory.pm:1601, :1609) after the source's values were pushed under
/// the same key (:1300-1304, the row's allele records at :1426-1435), so the
/// row of the transcript in the MANE Select set holds the record's name and
/// then `MANE_Select`, and every other row the record's name alone.
#[test]
fn a_custom_key_named_mane_lists_the_sources_values_before_the_rows_sets() {
    let tmp = tempfile::tempdir().unwrap();
    assert_text_formats(tmp.path(), &MANE_SELECT);

    let tab = written_tab_lines(tmp.path(), &MANE_SELECT);
    let cells = key_cells(&tab, "MANE");
    assert_eq!(cells.len(), 16);
    let both: Vec<&&str> = cells
        .iter()
        .filter(|cell| **cell == "multi_7,MANE_Select")
        .collect();
    assert_eq!(both.len(), 1, "one transcript of the sixteen is in the set");
    assert_eq!(cells.iter().filter(|cell| **cell == "multi_7").count(), 15);
}

/// Ensembl VEP pushes a co-located variant's synonyms onto `VAR_SYNONYMS`
/// after the source's values (OutputFactory.pm:1090) and joins the whole list
/// by two hyphens wherever the record has a co-located variant, with
/// `--var_synonyms` off too (:1141). Under the flag the cell is the record's
/// name, two hyphens and the synonyms; without it, the names of two sources'
/// records joined by two hyphens in every format.
#[test]
fn a_custom_key_named_var_synonyms_is_joined_by_two_hyphens_beside_a_co_located_variant() {
    let tmp = tempfile::tempdir().unwrap();
    assert_text_formats(tmp.path(), &VAR_SYNONYMS_ON);
    assert_text_formats(tmp.path(), &VAR_SYNONYMS_OF_TWO_SOURCES);

    let tab = written_tab_lines(tmp.path(), &VAR_SYNONYMS_ON);
    let cells = key_cells(&tab, "VAR_SYNONYMS");
    assert_eq!(cells.len(), 16);
    for cell in cells {
        assert_eq!(cell, "multi_7--ClinVar::RCV002273474,VCV001699617");
    }
    let tab = written_tab_lines(tmp.path(), &VAR_SYNONYMS_OF_TWO_SOURCES);
    let cells = key_cells(&tab, "VAR_SYNONYMS");
    assert_eq!(cells.len(), 16);
    for cell in cells {
        assert_eq!(cell, "multi_7--multi_7");
    }
}

/// The cell vep-rs writes where Ensembl VEP writes `cell` under a key that
/// spells a list the input record's hash holds: each item once, the samples'
/// zygosities in the order `samples` names them. Ensembl VEP lists the
/// zygosities in the order of a hash's keys (OutputFactory.pm:983, the
/// documented difference `individual_zyg_list_order`) and each row pushes its
/// allele's record onto the one list every row of the record shares (:391,
/// :1426-1435), which the rows print once all are built (:266-271): the
/// documented difference `custom_exact_source_allele_record_accumulation`.
fn each_item_once_in_sample_order(cell: &str, separator: &str, samples: &[&str]) -> String {
    let mut items: Vec<&str> = Vec::new();
    for item in cell.split(separator) {
        if !items.contains(&item) {
            items.push(item);
        }
    }
    let sample_of = |item: &str| {
        item.split_once(':')
            .and_then(|(name, _)| samples.iter().position(|sample| *sample == name))
    };
    let mut zygosities: Vec<&str> = items
        .iter()
        .copied()
        .filter(|item| sample_of(item).is_some())
        .collect();
    zygosities.sort_by_key(|item| sample_of(item));
    let records = items
        .iter()
        .copied()
        .filter(|item| sample_of(item).is_none());
    let ordered: Vec<&str> = zygosities.iter().copied().chain(records).collect();
    ordered.join(separator)
}

/// The sample names of the record's corpus input, in the order of its column
/// line.
fn samples(record: &Record) -> Vec<String> {
    let input = fs::read_to_string(corpus(record.corpus).join("variants.vcf")).unwrap();
    input
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("the column line")
        .split('\t')
        .skip(9)
        .map(str::to_string)
        .collect()
}

/// `lines`, the column line first, with the cell under `key` on each row
/// rewritten by `rewrite`.
fn with_key_cells(lines: &[String], key: &str, rewrite: impl Fn(&str) -> String) -> Vec<String> {
    let place = listings(&lines[0], '\t', key)[0];
    let mut rewritten = vec![lines[0].clone()];
    for line in &lines[1..] {
        let cells: Vec<String> = line
            .split('\t')
            .enumerate()
            .map(|(i, cell)| {
                if i == place {
                    rewrite(cell)
                } else {
                    cell.to_string()
                }
            })
            .collect();
        rewritten.push(cells.join("\t"));
    }
    rewritten
}

/// The record's hash sets the `--individual_zyg` list under `ZYG`
/// (OutputFactory.pm:981-988) before any record of a source is pushed, and
/// each row pushes its allele's record onto that list (:1426-1435). Ensembl
/// VEP's rows share the one list (:391), so each of the record's two rows
/// prints the eight zygosities and then the record's name twice. vep-rs
/// writes the zygosities in the header's order and then the record's name
/// once on each row of its allele, and writes every other cell as Ensembl VEP
/// does. The JSON output carries the list on the input record and the
/// source's record under `custom_annotations` in both engines.
#[test]
fn a_custom_key_named_zyg_lists_the_alleles_record_once_after_the_samples_zygosities() {
    let tmp = tempfile::tempdir().unwrap();
    let names = samples(&WITH_GENOTYPES);
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    assert_eq!(names.len(), 8);

    for (reference, cell) in [
        (
            &ZYG_OF_EVERY_SAMPLE,
            "S1:HET,S2:HOM,S3:HOM,S4:HOM,S5:HET,S6:HOM,S7:HET,S8:HOM,made_up_1",
        ),
        (&ZYG_OF_ONE_SAMPLE, "S3:HOM,made_up_1"),
    ] {
        let expected = with_key_cells(&tab_listed_once(reference), "ZYG", |cell| {
            each_item_once_in_sample_order(cell, ",", &names)
        });
        let written = written_tab_lines(tmp.path(), reference);
        assert_eq!(written, expected, "tab format, {:?}", reference.flags);
        assert_eq!(
            key_cells(&written, "ZYG"),
            [cell, cell],
            "{:?}",
            reference.flags
        );
    }

    let expected: Vec<String> = ZYG_OF_EVERY_SAMPLE
        .default
        .iter()
        .map(|row| {
            let (columns, extra) = row.rsplit_once('\t').expect("the Extra column");
            let extra: Vec<String> = extra
                .split(';')
                .map(|item| match item.strip_prefix("ZYG=") {
                    Some(cell) => {
                        format!("ZYG={}", each_item_once_in_sample_order(cell, ",", &names))
                    }
                    None => item.to_string(),
                })
                .collect();
            format!("{columns}\t{}", extra.join(";"))
        })
        .collect();
    let text = fs::read_to_string(annotate(tmp.path(), "default", &ZYG_OF_EVERY_SAMPLE)).unwrap();
    assert_eq!(default_rows(&text), expected, "default format");

    let listed_once = vcf_listed_once(&ZYG_OF_EVERY_SAMPLE);
    let place = listings(&listed_once[0], '|', "ZYG")[0];
    let (columns, entries) = listed_once[1].split_once("CSQ=").expect("the CSQ key");
    let entries: Vec<String> = entries
        .split(',')
        .map(|entry| {
            let cells: Vec<String> = entry
                .split('|')
                .enumerate()
                .map(|(i, cell)| {
                    if i == place {
                        each_item_once_in_sample_order(cell, "&", &names)
                    } else {
                        cell.to_string()
                    }
                })
                .collect();
            cells.join("|")
        })
        .collect();
    let expected = [
        listed_once[0].clone(),
        format!("{columns}CSQ={}", entries.join(",")),
    ];
    let text = fs::read_to_string(annotate(tmp.path(), "vcf", &ZYG_OF_EVERY_SAMPLE)).unwrap();
    assert_eq!(vcf_lines(&text), expected, "VCF format");

    let mut expected: serde_json::Value = serde_json::from_str(ZYG_OF_EVERY_SAMPLE_JSON).unwrap();
    let listed: Vec<&str> = expected["zyg"]
        .as_array()
        .expect("the record's zygosities")
        .iter()
        .map(|item| item.as_str().expect("a zygosity"))
        .collect();
    let in_header_order = each_item_once_in_sample_order(&listed.join(","), ",", &names);
    expected["zyg"] = in_header_order.split(',').collect();
    let record = json_record(tmp.path(), &ZYG_OF_EVERY_SAMPLE);
    assert_eq!(record, expected);
    let entries = record["transcript_consequences"]
        .as_array()
        .expect("the transcript entries");
    assert_eq!(entries.len(), 2);
    for entry in entries {
        assert!(entry.get("zyg").is_none(), "{entry}");
    }
}

/// The further run of the reproduction record of
/// `custom_exact_source_allele_record_accumulation`: one substitution with one
/// heterozygous sample, under `--individual_zyg all` and an `exact` source
/// named `ZYG` that holds the record `seqrec` on its allele. Ensembl VEP's
/// recorded cell is the zygosity and then `seqrec` once for each row of the
/// allele (OutputFactory.pm:391, :1426-1435); vep-rs writes the zygosity and
/// `seqrec` once on each row.
#[test]
fn the_reproduction_run_under_a_key_named_zyg_lists_the_alleles_record_once_on_each_row() {
    let record = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../scripts/divergences/reproductions")
            .join("custom_exact_source_allele_record_accumulation.toml"),
    )
    .unwrap();
    let input = "21\t10542415\t.\tA\tG\t.\t.\t.\tGT\t0|1";
    let source = "file={stage}/allele.vcf.gz,short_name=ZYG,format=vcf,type=exact,fields=TAG";
    for stated in [
        input.to_string(),
        "samples = [\"S1\"]".to_string(),
        format!("flags = [\"--individual_zyg\", \"all\", \"--custom\", \"{source}\"]"),
        r"INFO\\n21\\t10542415\\tseqrec\\tA\\tG\\t.\\tPASS\\tTAG=seq\\n' | bgzip >".to_string(),
        "json_cache = \"tests/golden/116/GRCh38-custom/json_cache\"".to_string(),
        "vep_rs = [\"S1:HET,seqrec\"]".to_string(),
    ] {
        assert!(record.contains(&stated), "the record states {stated:?}");
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    common::indexed_vcf(
        dir,
        "allele.vcf.gz",
        &[
            "##fileformat=VCFv4.2",
            "##INFO=<ID=TAG,Number=1,Type=String,Description=\"Record tag\">",
            "##contig=<ID=21>",
        ],
        "21\t10542415\tseqrec\tA\tG\t.\tPASS\tTAG=seq",
    );
    let variants = dir.join("input.vcf");
    fs::write(
        &variants,
        format!(
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n{input}\n"
        ),
    )
    .unwrap();
    let out = dir.join("out.txt");
    let run = common::isolated_vep(dir)
        .arg("-i")
        .arg(&variants)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-custom").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(["--individual_zyg", "all", "--custom"])
        .arg(source.replace("{stage}", &dir.display().to_string()))
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let text = fs::read_to_string(&out).unwrap();
    let rows = default_rows(&text);
    assert_eq!(rows.len(), 9);
    for row in &rows {
        let extra = row.rsplit('\t').next().expect("the Extra column");
        let cells: Vec<&str> = extra
            .split(';')
            .filter_map(|item| item.strip_prefix("ZYG="))
            .collect();
        assert_eq!(cells, ["S1:HET,seqrec"], "{row}");
    }
    let repeated = format!("ensembl_vep = [\"S1:HET{}\"]", ",seqrec".repeat(rows.len()));
    assert!(record.contains(&repeated), "the record states {repeated:?}");
}

/// A row pushes the populations of its maximum onto `MAX_AF_POPS` where the
/// maximum equals the one the row holds, which starts at 0
/// (OutputFactory.pm:1267, :1274), after the row's allele record was pushed
/// under the same key (:1426-1435, called before :1441-1443). At a maximum of
/// 0 the cell therefore holds the record's name and then the populations, in
/// both engines and in every text format; the JSON output keeps the
/// frequencies under `colocated_variants` and the source's record under
/// `custom_annotations`.
#[test]
fn a_custom_key_named_max_af_pops_lists_the_sources_values_before_the_populations_of_a_zero_maximum(
) {
    let tmp = tempfile::tempdir().unwrap();
    assert_text_formats(tmp.path(), &MAX_AF_POPS_OF_A_ZERO_MAXIMUM);

    let populations = "gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,\
                       gnomADe_MID,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS";
    let written = written_tab_lines(tmp.path(), &MAX_AF_POPS_OF_A_ZERO_MAXIMUM);
    let cells = key_cells(&written, "MAX_AF_POPS");
    assert_eq!(cells.len(), 6);
    for cell in cells {
        assert_eq!(cell, format!("made_up_2,{populations}"));
    }
    for cell in key_cells(&written, "MAX_AF") {
        assert_eq!(cell, "0");
    }

    let expected: serde_json::Value =
        serde_json::from_str(MAX_AF_POPS_OF_A_ZERO_MAXIMUM_JSON).unwrap();
    assert_eq!(
        json_record(tmp.path(), &MAX_AF_POPS_OF_A_ZERO_MAXIMUM),
        expected
    );
}

/// A maximum above the one the row holds assigns `MAX_AF_POPS` anew
/// (OutputFactory.pm:1269-1272), so the record's name, pushed under the key
/// before (:1426-1435), is gone: the row of the allele with a frequency above
/// 0 holds the maximum's population alone, in both engines and in every text
/// format, and the row of the other allele, which has neither a record of the
/// source nor a frequency, holds nothing under the key.
#[test]
fn a_custom_key_named_max_af_pops_holds_the_populations_alone_beside_a_maximum_above_zero() {
    let tmp = tempfile::tempdir().unwrap();
    assert_text_formats(tmp.path(), &MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM);

    let written = written_tab_lines(tmp.path(), &MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM);
    assert_eq!(key_cells(&written, "MAX_AF_POPS"), ["gnomADg_NFE", "-"]);
    assert_eq!(key_cells(&written, "MAX_AF"), ["1.993e-05", "-"]);
    assert_eq!(key_cells(&written, "MAX_AF_POPS_AF"), ["0.5", "-"]);

    let expected: serde_json::Value =
        serde_json::from_str(MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM_JSON).unwrap();
    assert_eq!(
        json_record(tmp.path(), &MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM),
        expected
    );
}

/// Ensembl VEP assigns the nearest names over the whole record's values in
/// the record's hash (OutputFactory.pm:1017), and each row pushes its
/// allele's record onto that one list (:1426-1435), which every row's hash
/// shares (:391): each of its sixteen rows lists the symbol and then the
/// record's name sixteen times. vep-rs lists the symbol and then the record
/// once on each row of its allele, and writes every other cell as Ensembl VEP
/// does.
#[test]
fn a_custom_key_named_nearest_lists_the_alleles_record_once_after_the_nearest_names() {
    let tmp = tempfile::tempdir().unwrap();
    let reference = tab_listed_once(&NEAREST_SYMBOL);
    let place = listings(&reference[0], '\t', "NEAREST")[0];
    let repeated = format!("SON{}", ",multi_7".repeat(16));
    let cells = key_cells(&reference, "NEAREST");
    assert_eq!(cells.len(), 16);
    for cell in cells {
        assert_eq!(cell, repeated);
    }

    let expected: Vec<String> = reference
        .iter()
        .map(|line| {
            let cells: Vec<&str> = line
                .split('\t')
                .enumerate()
                .map(|(i, cell)| {
                    if i == place && cell == repeated {
                        "SON,multi_7"
                    } else {
                        cell
                    }
                })
                .collect();
            cells.join("\t")
        })
        .collect();
    let text = fs::read_to_string(annotate(tmp.path(), "tab", &NEAREST_SYMBOL)).unwrap();
    assert_eq!(tab_lines(&text), expected);
}

/// The JSON factory moves the custom annotations to the record before any row
/// hash is built (OutputFactory/JSON.pm:249-258), so the JSON output prints
/// the source's record under `custom_annotations` and each entry's `mane`
/// list holds the sets alone.
#[test]
fn a_custom_key_named_as_a_list_field_stays_out_of_the_json_entries() {
    let tmp = tempfile::tempdir().unwrap();
    let expected: serde_json::Value = serde_json::from_str(MANE_SELECT_JSON).unwrap();
    let record = json_record(tmp.path(), &MANE_SELECT);
    assert_eq!(record, expected);

    let entries = record["transcript_consequences"]
        .as_array()
        .expect("the transcript entries");
    let listed: Vec<&serde_json::Value> = entries
        .iter()
        .filter_map(|entry| entry.get("mane"))
        .collect();
    assert_eq!(listed, [&serde_json::json!(["MANE_Select"])]);
}

/// The helpers that drop Ensembl VEP's second listing take one cell out of
/// each line of a reference that lists the key twice and leave one that lists
/// it once as it is.
#[test]
fn the_second_listing_is_dropped_from_the_column_line_and_from_each_row() {
    assert_eq!(
        tab_listed_once(&VAR_SYNONYMS_OF_TWO_SOURCES),
        VAR_SYNONYMS_OF_TWO_SOURCES.tab
    );
    assert_eq!(
        vcf_listed_once(&VAR_SYNONYMS_OF_TWO_SOURCES),
        VAR_SYNONYMS_OF_TWO_SOURCES.vcf
    );

    for reference in [&MANE_SELECT, &VAR_SYNONYMS_ON, &NEAREST_SYMBOL] {
        let key = reference.key;
        let twice = reference.tab[0].split('\t').count();
        assert_eq!(listings(reference.tab[0], '\t', key).len(), 2, "{key}");
        let tab = tab_listed_once(reference);
        assert_eq!(listings(&tab[0], '\t', key).len(), 1, "{key}");
        assert_eq!(tab.len(), 17, "{key}");
        for line in &tab {
            assert_eq!(line.split('\t').count(), twice - 1, "{key}");
        }
    }

    assert_eq!(
        each_item_once_in_sample_order("B:HOM,A:HET,name,name", ",", &["A", "B"]),
        "A:HET,B:HOM,name"
    );
    assert_eq!(
        each_item_once_in_sample_order("B:HOM&name&name&name", "&", &["A", "B"]),
        "B:HOM&name"
    );
    let lines = ["#one\tKEY\tthree".to_string(), "1\tb,a,a\t3".to_string()];
    assert_eq!(
        with_key_cells(&lines, "KEY", |cell| cell.replace(",a,a", ",a")),
        ["#one\tKEY\tthree", "1\tb,a\t3"]
    );

    for reference in [&MANE_SELECT, &VAR_SYNONYMS_ON] {
        let key = reference.key;
        let twice = reference.vcf[0].split('|').count();
        assert_eq!(listings(reference.vcf[0], '|', key).len(), 2, "{key}");
        let vcf = vcf_listed_once(reference);
        assert_eq!(listings(&vcf[0], '|', key).len(), 1, "{key}");
        let entries = vcf[1].split_once("CSQ=").unwrap().1;
        assert_eq!(entries.split(',').count(), 16, "{key}");
        for entry in entries.split(',') {
            assert_eq!(entry.split('|').count(), twice - 1, "{key}");
        }
    }
}

/// The Parquet output stores the key as a string that holds the `--tab` cell
/// of each row.
#[cfg(feature = "parquet")]
#[test]
fn a_custom_key_named_as_a_list_field_is_stored_in_parquet_as_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    let names = samples(&WITH_GENOTYPES);
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    let zygosities = with_key_cells(&tab_listed_once(&ZYG_OF_EVERY_SAMPLE), "ZYG", |cell| {
        each_item_once_in_sample_order(cell, ",", &names)
    });
    for (reference, tab) in [
        (&MANE_SELECT, tab_listed_once(&MANE_SELECT)),
        (&VAR_SYNONYMS_ON, tab_listed_once(&VAR_SYNONYMS_ON)),
        (
            &VAR_SYNONYMS_OF_TWO_SOURCES,
            tab_listed_once(&VAR_SYNONYMS_OF_TWO_SOURCES),
        ),
        (
            &MAX_AF_POPS_OF_A_ZERO_MAXIMUM,
            tab_listed_once(&MAX_AF_POPS_OF_A_ZERO_MAXIMUM),
        ),
        (
            &MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM,
            tab_listed_once(&MAX_AF_POPS_OF_A_POSITIVE_MAXIMUM),
        ),
        (&ZYG_OF_EVERY_SAMPLE, zygosities),
    ] {
        let place = listings(&tab[0], '\t', reference.key)[0];
        let mut expected: Vec<&str> = key_cells(&tab, reference.key);
        expected.sort_unstable();
        let out = annotate(tmp.path(), "parquet", reference);
        let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
            .unwrap()
            .iter()
            .map(|row| {
                row.split('\t')
                    .nth(place)
                    .expect("the key's cell")
                    .to_string()
            })
            .collect();
        cells.sort_unstable();
        assert_eq!(cells, expected, "{} {:?}", reference.key, reference.flags);
    }
}
