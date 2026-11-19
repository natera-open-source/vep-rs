#!/usr/bin/env python3
"""Write a synthetic plugin data file for dbNSFP or dbscSNV at a corpus's own sites.

dbNSFP is free for academic use and needs a paid agreement for commercial use, and dbscSNV
publishes no licence text. Neither data file is redistributed here, so their plugins are verified
on made-up fixtures: a file with the layout's
header and column order, prepared, sorted and indexed as Ensembl VEP's plugin documents for
the assembly (dbNSFP.pm:62-74, dbscSNV.pm:49-59), whose values are fabricated from a seed.
The rows sit at the corpus's own
single-nucleotide sites, and the columns a plugin matches on beyond position and allele (the
transcript ids, the amino acid change) are read from a VEP-format output of the corpus, so a
row the plugin must pick carries the transcript and the amino acids of a real consequence row
and sits beside rows it must not pick. A regeneration from the same inputs and seed writes the
same bytes.

Layouts (`--layout`), each prepared for `--assembly GRCh37` or `GRCh38`:

    dbnsfp      a made-up dbNSFP table under the 458 column names this generator writes
                (`#chr pos(1-based) ref alt aaref aaalt rs_dbSNP hg19_chr hg19_pos(1-based)
                hg18_chr hg18_pos(1-based) aapos genename Ensembl_geneid Ensembl_transcriptid
                ...`); one row per alternate base
                at every single-nucleotide site: for the input's own alternate allele one row
                per distinct amino acid change the output shows, listing the transcripts of
                that change, and for the two other bases one row with a fabricated change; the
                per-transcript columns carry one `;`-joined value per listed transcript. The
                GRCh38 file carries the corpus coordinates in `#chr`/`pos(1-based)`, is sorted
                on them and indexed `tabix -s 1 -b 2 -e 2`; the GRCh37 file carries them in
                `hg19_chr`/`hg19_pos(1-based)`, is sorted on those and indexed
                `tabix -s 8 -b 9 -e 9`; the other assembly's columns carry a shifted position.
                One row is planted on chromosome 1 inside the region Ensembl VEP's plugin reads
                to detect the `;`-joined columns, with two transcripts in every one of them.
    dbscsnv     a made-up dbscSNV table under the 18 column names this generator writes
                (`chr pos ref alt hg38_chr hg38_pos RefSeq? Ensembl? ... Ensembl_gene ...
                ada_score rf_score`), the header without a
                leading `#` (tabix is told it starts with `c`, dbscSNV.pm:61-62) and every line
                the generator writes ending in CRLF; one row per
                alternate base at every single-nucleotide site. The GRCh37 file carries the
                corpus coordinates in `chr`/`pos` and is indexed `tabix -s 1 -b 2 -e 2 -c c`;
                the GRCh38 file carries them in `hg38_chr`/`hg38_pos`, is sorted on those and
                indexed `tabix -s 5 -b 6 -e 6 -c c`.
    cadd_snv, cadd_indel, cadd_sv, spliceai_snv, spliceai_indel
                registered names whose writers are not implemented here; asking for one is an
                error naming it.

A site is every record of the VCF whose REF and ALT are single bases (one per ALT of a
multi-allelic record); a record of any other shape writes nothing. A site named by
`--miss chr:pos` gets no row. The generator spells chromosomes with the `chr` prefix dropped
and `MT` written `M` for dbnsfp, and with the `chr` prefix dropped for
dbscsnv. A row's generator is seeded from the seed, the layout and the row's chromosome,
position, reference and alternate allele (and, for a dbnsfp row, its amino acid change)
alone, so a row does not move when another row is added.

The VEP-format output (`--reference`, a default- or tab-format file, plain or gzip) supplies,
per site and alternate allele, the transcript rows: `Feature`, `Gene`, `Consequence`,
`Protein_position`, `Amino_acids` and, where the output carries it, `SYMBOL`. A site with no
row in the output gets rows listing no transcript.

bgzip and tabix write the file and its index; the file's spec (layout, assembly, seed, the
preparation, the input digests, the row counts, the planted row, the shift) is written beside
it as `<file>.spec.json`, and `--provenance <corpus>/provenance.json` records the fixture under
`plugin_fixtures.<file name>` there, keeping every other key. An existing output is overwritten
only under --force.

Usage:
    make_plugin_fixture.py --layout dbnsfp --assembly GRCh38 --seed 38 \\
        --vcf tests/golden/116/<corpus>/variants.vcf \\
        --reference tests/golden/116/<corpus>/expected/default.txt.gz \\
        --out tests/golden/116/<corpus>/plugin_data/dbNSFP4.9a_grch38_synthetic.gz \\
        [--provenance tests/golden/116/<corpus>/provenance.json] [--miss 21:19268555]... [--force]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

GENERATOR = "scripts/fixtures/make_plugin_fixture.py"
TOOLS = ("bgzip", "tabix")
BASES = "ACGT"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
ASSEMBLIES = ("GRCh37", "GRCh38")
SEQUENCE = re.compile(r"^[ACGTacgt]+$")
SITE = re.compile(r"[^:\s]+:[0-9]+")
# The position written in the columns of the assembly the file is not prepared for: the
# corpus position plus this many bases, so a reader of the wrong column matches nothing.
OTHER_ASSEMBLY_SHIFT = 100_000

# The region Ensembl VEP's dbNSFP plugin reads at start-up to find the `;`-joined columns
# (dbNSFP.pm:162), and the row planted inside it on both assemblies.
DETECTION_REGION = ("1", 1_008_170, 1_082_927)
PLANTED_ROW = {"chr": "1", "pos": 1_045_000, "ref": "A", "alt": "G"}
PLANTED_TRANSCRIPTS = ("ENST00000900001", "ENST00000900002")
PLANTED_GENES = ("ENSG00000900001", "ENSG00000900002")

# The header of the generator's made-up dbNSFP table: 458 column names, in order.
DBNSFP_COLUMNS = (
    "#chr pos(1-based) ref alt aaref aaalt rs_dbSNP hg19_chr hg19_pos(1-based) hg18_chr "
    "hg18_pos(1-based) aapos genename Ensembl_geneid Ensembl_transcriptid Ensembl_proteinid "
    "Uniprot_acc Uniprot_entry HGVSc_ANNOVAR HGVSp_ANNOVAR HGVSc_snpEff HGVSp_snpEff HGVSc_VEP "
    "HGVSp_VEP APPRIS GENCODE_basic TSL VEP_canonical cds_strand refcodon codonpos codon_degeneracy "
    "Ancestral_allele AltaiNeandertal Denisova VindijiaNeandertal ChagyrskayaNeandertal SIFT_score "
    "SIFT_converted_rankscore SIFT_pred SIFT4G_score SIFT4G_converted_rankscore SIFT4G_pred "
    "Polyphen2_HDIV_score Polyphen2_HDIV_rankscore Polyphen2_HDIV_pred Polyphen2_HVAR_score "
    "Polyphen2_HVAR_rankscore Polyphen2_HVAR_pred LRT_score LRT_converted_rankscore LRT_pred "
    "LRT_Omega MutationTaster_score MutationTaster_converted_rankscore MutationTaster_pred "
    "MutationTaster_model MutationTaster_AAE MutationAssessor_score MutationAssessor_rankscore "
    "MutationAssessor_pred FATHMM_score FATHMM_converted_rankscore FATHMM_pred PROVEAN_score "
    "PROVEAN_converted_rankscore PROVEAN_pred VEST4_score VEST4_rankscore MetaSVM_score "
    "MetaSVM_rankscore MetaSVM_pred MetaLR_score MetaLR_rankscore MetaLR_pred Reliability_index "
    "MetaRNN_score MetaRNN_rankscore MetaRNN_pred M-CAP_score M-CAP_rankscore M-CAP_pred "
    "REVEL_score REVEL_rankscore MutPred_score MutPred_rankscore MutPred_protID MutPred_AAchange "
    "MutPred_Top5features MVP_score MVP_rankscore gMVP_score gMVP_rankscore MPC_score MPC_rankscore "
    "PrimateAI_score PrimateAI_rankscore PrimateAI_pred DEOGEN2_score DEOGEN2_rankscore "
    "DEOGEN2_pred BayesDel_addAF_score BayesDel_addAF_rankscore BayesDel_addAF_pred "
    "BayesDel_noAF_score BayesDel_noAF_rankscore BayesDel_noAF_pred ClinPred_score "
    "ClinPred_rankscore ClinPred_pred LIST-S2_score LIST-S2_rankscore LIST-S2_pred VARITY_R_score "
    "VARITY_R_rankscore VARITY_ER_score VARITY_ER_rankscore VARITY_R_LOO_score "
    "VARITY_R_LOO_rankscore VARITY_ER_LOO_score VARITY_ER_LOO_rankscore ESM1b_score ESM1b_rankscore "
    "ESM1b_pred EVE_score EVE_rankscore EVE_Class10_pred EVE_Class20_pred EVE_Class25_pred "
    "EVE_Class30_pred EVE_Class40_pred EVE_Class50_pred EVE_Class60_pred EVE_Class70_pred "
    "EVE_Class75_pred EVE_Class80_pred EVE_Class90_pred AlphaMissense_score AlphaMissense_rankscore "
    "AlphaMissense_pred PHACTboost_score PHACTboost_rankscore MutFormer_score MutFormer_rankscore "
    "MutScore_score MutScore_rankscore Aloft_Fraction_transcripts_affected Aloft_prob_Tolerant "
    "Aloft_prob_Recessive Aloft_prob_Dominant Aloft_pred Aloft_Confidence CADD_raw "
    "CADD_raw_rankscore CADD_phred CADD_raw_hg19 CADD_raw_rankscore_hg19 CADD_phred_hg19 DANN_score "
    "DANN_rankscore fathmm-MKL_coding_score fathmm-MKL_coding_rankscore fathmm-MKL_coding_pred "
    "fathmm-MKL_coding_group fathmm-XF_coding_score fathmm-XF_coding_rankscore "
    "fathmm-XF_coding_pred Eigen-raw_coding Eigen-raw_coding_rankscore Eigen-phred_coding "
    "Eigen-PC-raw_coding Eigen-PC-raw_coding_rankscore Eigen-PC-phred_coding GenoCanyon_score "
    "GenoCanyon_rankscore integrated_fitCons_score integrated_fitCons_rankscore "
    "integrated_confidence_value GM12878_fitCons_score GM12878_fitCons_rankscore "
    "GM12878_confidence_value H1-hESC_fitCons_score H1-hESC_fitCons_rankscore "
    "H1-hESC_confidence_value HUVEC_fitCons_score HUVEC_fitCons_rankscore HUVEC_confidence_value "
    "LINSIGHT LINSIGHT_rankscore GERP++_NR GERP++_RS GERP++_RS_rankscore GERP_91_mammals "
    "GERP_91_mammals_rankscore phyloP100way_vertebrate phyloP100way_vertebrate_rankscore "
    "phyloP470way_mammalian phyloP470way_mammalian_rankscore phyloP17way_primate "
    "phyloP17way_primate_rankscore phastCons100way_vertebrate phastCons100way_vertebrate_rankscore "
    "phastCons470way_mammalian phastCons470way_mammalian_rankscore phastCons17way_primate "
    "phastCons17way_primate_rankscore SiPhy_29way_pi SiPhy_29way_logOdds "
    "SiPhy_29way_logOdds_rankscore bStatistic bStatistic_converted_rankscore 1000Gp3_AC 1000Gp3_AF "
    "1000Gp3_AFR_AC 1000Gp3_AFR_AF 1000Gp3_EUR_AC 1000Gp3_EUR_AF 1000Gp3_AMR_AC 1000Gp3_AMR_AF "
    "1000Gp3_EAS_AC 1000Gp3_EAS_AF 1000Gp3_SAS_AC 1000Gp3_SAS_AF TWINSUK_AC TWINSUK_AF ALSPAC_AC "
    "ALSPAC_AF UK10K_AC UK10K_AF ESP6500_AA_AC ESP6500_AA_AF ESP6500_EA_AC ESP6500_EA_AF ExAC_AC "
    "ExAC_AF ExAC_Adj_AC ExAC_Adj_AF ExAC_AFR_AC ExAC_AFR_AF ExAC_AMR_AC ExAC_AMR_AF ExAC_EAS_AC "
    "ExAC_EAS_AF ExAC_FIN_AC ExAC_FIN_AF ExAC_NFE_AC ExAC_NFE_AF ExAC_SAS_AC ExAC_SAS_AF "
    "ExAC_nonTCGA_AC ExAC_nonTCGA_AF ExAC_nonTCGA_Adj_AC ExAC_nonTCGA_Adj_AF ExAC_nonTCGA_AFR_AC "
    "ExAC_nonTCGA_AFR_AF ExAC_nonTCGA_AMR_AC ExAC_nonTCGA_AMR_AF ExAC_nonTCGA_EAS_AC "
    "ExAC_nonTCGA_EAS_AF ExAC_nonTCGA_FIN_AC ExAC_nonTCGA_FIN_AF ExAC_nonTCGA_NFE_AC "
    "ExAC_nonTCGA_NFE_AF ExAC_nonTCGA_SAS_AC ExAC_nonTCGA_SAS_AF ExAC_nonpsych_AC ExAC_nonpsych_AF "
    "ExAC_nonpsych_Adj_AC ExAC_nonpsych_Adj_AF ExAC_nonpsych_AFR_AC ExAC_nonpsych_AFR_AF "
    "ExAC_nonpsych_AMR_AC ExAC_nonpsych_AMR_AF ExAC_nonpsych_EAS_AC ExAC_nonpsych_EAS_AF "
    "ExAC_nonpsych_FIN_AC ExAC_nonpsych_FIN_AF ExAC_nonpsych_NFE_AC ExAC_nonpsych_NFE_AF "
    "ExAC_nonpsych_SAS_AC ExAC_nonpsych_SAS_AF gnomAD_exomes_flag gnomAD_exomes_AC gnomAD_exomes_AN "
    "gnomAD_exomes_AF gnomAD_exomes_nhomalt gnomAD_exomes_POPMAX_AC gnomAD_exomes_POPMAX_AN "
    "gnomAD_exomes_POPMAX_AF gnomAD_exomes_POPMAX_nhomalt gnomAD_exomes_AFR_AC gnomAD_exomes_AFR_AN "
    "gnomAD_exomes_AFR_AF gnomAD_exomes_AFR_nhomalt gnomAD_exomes_AMR_AC gnomAD_exomes_AMR_AN "
    "gnomAD_exomes_AMR_AF gnomAD_exomes_AMR_nhomalt gnomAD_exomes_ASJ_AC gnomAD_exomes_ASJ_AN "
    "gnomAD_exomes_ASJ_AF gnomAD_exomes_ASJ_nhomalt gnomAD_exomes_EAS_AC gnomAD_exomes_EAS_AN "
    "gnomAD_exomes_EAS_AF gnomAD_exomes_EAS_nhomalt gnomAD_exomes_FIN_AC gnomAD_exomes_FIN_AN "
    "gnomAD_exomes_FIN_AF gnomAD_exomes_FIN_nhomalt gnomAD_exomes_MID_AC gnomAD_exomes_MID_AN "
    "gnomAD_exomes_MID_AF gnomAD_exomes_MID_nhomalt gnomAD_exomes_NFE_AC gnomAD_exomes_NFE_AN "
    "gnomAD_exomes_NFE_AF gnomAD_exomes_NFE_nhomalt gnomAD_exomes_SAS_AC gnomAD_exomes_SAS_AN "
    "gnomAD_exomes_SAS_AF gnomAD_exomes_SAS_nhomalt gnomAD_exomes_non_ukb_AC "
    "gnomAD_exomes_non_ukb_AN gnomAD_exomes_non_ukb_AF gnomAD_exomes_non_ukb_nhomalt "
    "gnomAD_exomes_non_ukb_AFR_AC gnomAD_exomes_non_ukb_AFR_AN gnomAD_exomes_non_ukb_AFR_AF "
    "gnomAD_exomes_non_ukb_AFR_nhomalt gnomAD_exomes_non_ukb_AMR_AC gnomAD_exomes_non_ukb_AMR_AN "
    "gnomAD_exomes_non_ukb_AMR_AF gnomAD_exomes_non_ukb_AMR_nhomalt gnomAD_exomes_non_ukb_ASJ_AC "
    "gnomAD_exomes_non_ukb_ASJ_AN gnomAD_exomes_non_ukb_ASJ_AF gnomAD_exomes_non_ukb_ASJ_nhomalt "
    "gnomAD_exomes_non_ukb_EAS_AC gnomAD_exomes_non_ukb_EAS_AN gnomAD_exomes_non_ukb_EAS_AF "
    "gnomAD_exomes_non_ukb_EAS_nhomalt gnomAD_exomes_non_ukb_FIN_AC gnomAD_exomes_non_ukb_FIN_AN "
    "gnomAD_exomes_non_ukb_FIN_AF gnomAD_exomes_non_ukb_FIN_nhomalt gnomAD_exomes_non_ukb_MID_AC "
    "gnomAD_exomes_non_ukb_MID_AN gnomAD_exomes_non_ukb_MID_AF gnomAD_exomes_non_ukb_MID_nhomalt "
    "gnomAD_exomes_non_ukb_NFE_AC gnomAD_exomes_non_ukb_NFE_AN gnomAD_exomes_non_ukb_NFE_AF "
    "gnomAD_exomes_non_ukb_NFE_nhomalt gnomAD_exomes_non_ukb_SAS_AC gnomAD_exomes_non_ukb_SAS_AN "
    "gnomAD_exomes_non_ukb_SAS_AF gnomAD_exomes_non_ukb_SAS_nhomalt gnomAD_genomes_flag "
    "gnomAD_genomes_AC gnomAD_genomes_AN gnomAD_genomes_AF gnomAD_genomes_nhomalt "
    "gnomAD_genomes_POPMAX_AC gnomAD_genomes_POPMAX_AN gnomAD_genomes_POPMAX_AF "
    "gnomAD_genomes_POPMAX_nhomalt gnomAD_genomes_AFR_AC gnomAD_genomes_AFR_AN "
    "gnomAD_genomes_AFR_AF gnomAD_genomes_AFR_nhomalt gnomAD_genomes_AMI_AC gnomAD_genomes_AMI_AN "
    "gnomAD_genomes_AMI_AF gnomAD_genomes_AMI_nhomalt gnomAD_genomes_AMR_AC gnomAD_genomes_AMR_AN "
    "gnomAD_genomes_AMR_AF gnomAD_genomes_AMR_nhomalt gnomAD_genomes_ASJ_AC gnomAD_genomes_ASJ_AN "
    "gnomAD_genomes_ASJ_AF gnomAD_genomes_ASJ_nhomalt gnomAD_genomes_EAS_AC gnomAD_genomes_EAS_AN "
    "gnomAD_genomes_EAS_AF gnomAD_genomes_EAS_nhomalt gnomAD_genomes_FIN_AC gnomAD_genomes_FIN_AN "
    "gnomAD_genomes_FIN_AF gnomAD_genomes_FIN_nhomalt gnomAD_genomes_MID_AC gnomAD_genomes_MID_AN "
    "gnomAD_genomes_MID_AF gnomAD_genomes_MID_nhomalt gnomAD_genomes_NFE_AC gnomAD_genomes_NFE_AN "
    "gnomAD_genomes_NFE_AF gnomAD_genomes_NFE_nhomalt gnomAD_genomes_SAS_AC gnomAD_genomes_SAS_AN "
    "gnomAD_genomes_SAS_AF gnomAD_genomes_SAS_nhomalt ALFA_European_AC ALFA_European_AN "
    "ALFA_European_AF ALFA_African_Others_AC ALFA_African_Others_AN ALFA_African_Others_AF "
    "ALFA_East_Asian_AC ALFA_East_Asian_AN ALFA_East_Asian_AF ALFA_African_American_AC "
    "ALFA_African_American_AN ALFA_African_American_AF ALFA_Latin_American_1_AC "
    "ALFA_Latin_American_1_AN ALFA_Latin_American_1_AF ALFA_Latin_American_2_AC "
    "ALFA_Latin_American_2_AN ALFA_Latin_American_2_AF ALFA_Other_Asian_AC ALFA_Other_Asian_AN "
    "ALFA_Other_Asian_AF ALFA_South_Asian_AC ALFA_South_Asian_AN ALFA_South_Asian_AF ALFA_Other_AC "
    "ALFA_Other_AN ALFA_Other_AF ALFA_African_AC ALFA_African_AN ALFA_African_AF ALFA_Asian_AC "
    "ALFA_Asian_AN ALFA_Asian_AF ALFA_Total_AC ALFA_Total_AN ALFA_Total_AF clinvar_id "
    "clinvar_clnsig clinvar_trait clinvar_review clinvar_hgvs clinvar_var_source clinvar_MedGen_id "
    "clinvar_OMIM_id clinvar_Orphanet_id Interpro_domain GTEx_V8_eQTL_gene GTEx_V8_eQTL_tissue "
    "GTEx_V8_sQTL_gene GTEx_V8_sQTL_tissue eQTLGen_snp_id eQTLGen_gene_id eQTLGen_gene_symbol "
    "eQTLGen_cis_or_trans Geuvadis_eQTL_target_gene"
).split(" ")

# The header of the generator's made-up dbscSNV table.
DBSCSNV_COLUMNS = (
    "chr pos ref alt hg38_chr hg38_pos RefSeq? Ensembl? RefSeq_region RefSeq_gene "
    "RefSeq_functional_consequence RefSeq_id_c.change_p.change Ensembl_region Ensembl_gene "
    "Ensembl_functional_consequence Ensembl_id_c.change_p.change ada_score rf_score"
).split(" ")

# The dbNSFP columns this writer fills with one value per listed transcript, `;`-joined
# (the plugin splits them and Ensembl_transcriptid on `;`, dbNSFP.pm:353 and 376);
# every other column holds `.`, or one of the single values of SINGLE_FILLERS.
TRANSCRIPT_FILLERS = (
    "aapos", "genename", "Ensembl_geneid", "Ensembl_transcriptid", "Ensembl_proteinid",
    "Uniprot_acc", "HGVSc_VEP", "HGVSp_VEP", "APPRIS", "GENCODE_basic", "TSL", "VEP_canonical",
    "SIFT_score", "SIFT_pred", "SIFT4G_score", "SIFT4G_pred", "Polyphen2_HDIV_score",
    "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "FATHMM_score",
    "FATHMM_pred", "PROVEAN_score", "PROVEAN_pred", "VEST4_score", "REVEL_score", "MVP_score",
)
SPLICEAI_AND_CADD = ("cadd_snv", "cadd_indel", "cadd_sv", "spliceai_snv", "spliceai_indel")
LAYOUTS = ("dbnsfp", "dbscsnv", *SPLICEAI_AND_CADD)


def fail(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [make_plugin_fixture] {message}")


@dataclass(frozen=True)
class Preparation:
    """How a layout's file is prepared for one assembly: the 0-based columns that carry the
    corpus coordinates, the columns that carry the other assembly's, the tabix arguments and
    the line ending."""

    header: tuple[str, ...]
    chr_col: int
    pos_col: int
    other_chr_col: int
    other_pos_col: int
    tabix: tuple[str, ...]
    line_end: str

    @property
    def index_text(self) -> str:
        return "tabix " + " ".join(self.tabix)


def preparation(layout: str, assembly: str) -> Preparation:
    if layout == "dbnsfp":
        header = tuple(DBNSFP_COLUMNS)
        if assembly == "GRCh38":
            return Preparation(header, 0, 1, 7, 8, ("-s", "1", "-b", "2", "-e", "2"), "\n")
        return Preparation(header, 7, 8, 0, 1, ("-s", "8", "-b", "9", "-e", "9"), "\n")
    if layout == "dbscsnv":
        header = tuple(DBSCSNV_COLUMNS)
        if assembly == "GRCh37":
            return Preparation(header, 0, 1, 4, 5, ("-s", "1", "-b", "2", "-e", "2", "-c", "c"), "\r\n")
        return Preparation(header, 4, 5, 0, 1, ("-s", "5", "-b", "6", "-e", "6", "-c", "c"), "\r\n")
    raise fail(f"layout {layout} is registered but has no writer here")


@dataclass(frozen=True)
class Site:
    """One single-nucleotide site of the corpus: the chromosome as the input spells it, the
    position, the reference base and the input's own alternate bases."""

    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]


def read_sites(vcf: Path) -> list[Site]:
    sites: dict[tuple[str, int, str], list[str]] = {}
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise fail(f"{vcf}: a record has fewer than five columns: {line!r}")
        chrom, pos, ref = fields[0], int(fields[1]), fields[3].upper()
        if len(ref) != 1 or not SEQUENCE.match(ref):
            continue
        alts = sites.setdefault((chrom, pos, ref), [])
        for alt in fields[4].split(","):
            alt = alt.upper()
            if len(alt) == 1 and SEQUENCE.match(alt) and alt != ref and alt not in alts:
                alts.append(alt)
    return [Site(c, p, r, tuple(alts)) for (c, p, r), alts in sites.items() if alts]


@dataclass(frozen=True)
class TranscriptRow:
    """One transcript row of the VEP-format output at a site and allele."""

    feature: str
    gene: str
    consequence: str
    protein_position: str
    amino_acids: str
    symbol: str


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def read_reference(path: Path) -> dict[tuple[str, int, str], list[TranscriptRow]]:
    """The transcript rows of a default- or tab-format VEP output keyed by (chromosome,
    position, allele) of a single-base `Location`; `SYMBOL` from its column or the `Extra`
    field when the output carries it, else empty."""
    rows: dict[tuple[str, int, str], list[TranscriptRow]] = {}
    columns: dict[str, int] | None = None
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##") or not line:
                continue
            if line.startswith("#"):
                columns = {name: i for i, name in enumerate(line[1:].split("\t"))}
                continue
            if columns is None:
                raise fail(f"{path}: no #Uploaded_variation header line before the first row")
            fields = line.split("\t")

            def cell(name: str) -> str:
                index = columns.get(name)
                return fields[index] if index is not None and index < len(fields) else ""

            if cell("Feature_type") != "Transcript":
                continue
            location, allele = cell("Location"), cell("Allele")
            if "-" in location.rpartition(":")[2] or allele not in BASES:
                continue
            chrom, _, pos = location.rpartition(":")
            symbol = cell("SYMBOL")
            if not symbol:
                extra = dict(pair.split("=", 1) for pair in cell("Extra").split(";") if "=" in pair)
                symbol = extra.get("SYMBOL", "")
            rows.setdefault((chrom, int(pos), allele), []).append(
                TranscriptRow(cell("Feature"), cell("Gene"), cell("Consequence"), cell("Protein_position"), cell("Amino_acids"), symbol)
            )
    if columns is None:
        raise fail(f"{path}: no #Uploaded_variation header line")
    return rows


def spell(chrom: str, layout: str) -> str:
    """The chromosome as the generator spells it for the layout."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if layout == "dbnsfp" and name.upper() == "MT":
        return "M"
    return name


def row_rng(seed: int, layout: str, chrom: str, pos: int, ref: str, alt: str, tag: str = "") -> random.Random:
    """One generator per row, seeded from the row's identity, so a value depends on nothing
    but the seed and the row."""
    return random.Random(f"{seed}:{layout}:{chrom}:{pos}:{ref}:{alt}:{tag}")


def dbnsfp_amino_acid(peptide: str) -> str:
    """A peptide letter as the plugin reads it from `aaref` and `aaalt`: the stop codon is
    `X` (dbNSFP.pm:342)."""
    return "X" if peptide == "*" else peptide


@dataclass
class DbnsfpRow:
    """One dbNSFP row before it is laid out: the alleles, the amino acid change, the
    transcripts whose per-transcript columns the row lists, in order, and whether the change
    is one the output shows for the input's own allele."""

    chrom: str
    pos: int
    ref: str
    alt: str
    aaref: str
    aaalt: str
    transcripts: list[TranscriptRow] = field(default_factory=list)
    from_output: bool = False


def dbnsfp_rows_for_site(site: Site, reference: dict, layout: str, seed: int) -> list[DbnsfpRow]:
    """The dbNSFP rows of one site: for the input's own alternate bases one row per distinct
    amino acid change among the output's transcript rows (`X/Y`, the reference peptide and
    the allele's), listing those transcripts; an alternate base with no such row, the input's
    or another, gets one row with a fabricated change listing every transcript the output
    shows at the site."""
    chrom = spell(site.chrom, layout)
    every: list[TranscriptRow] = []
    seen: set[str] = set()
    for alt in BASES:
        for row in reference.get((site.chrom, site.pos, alt), []):
            if row.feature not in seen:
                seen.add(row.feature)
                every.append(row)
    every.sort(key=lambda row: row.feature)
    rows: list[DbnsfpRow] = []
    for alt in BASES:
        if alt == site.ref:
            continue
        changes: dict[tuple[str, str], list[TranscriptRow]] = {}
        if alt in site.alts:
            for row in reference.get((site.chrom, site.pos, alt), []):
                parts = row.amino_acids.split("/")
                if len(parts) == 2 and parts[0] != parts[1] and all(len(p) == 1 for p in parts):
                    changes.setdefault((parts[0], parts[1]), []).append(row)
        if changes:
            for (aaref, aaalt), transcripts in sorted(changes.items()):
                transcripts.sort(key=lambda row: row.feature)
                rows.append(DbnsfpRow(chrom, site.pos, site.ref, alt, dbnsfp_amino_acid(aaref), dbnsfp_amino_acid(aaalt), transcripts, True))
        else:
            rng = row_rng(seed, layout, chrom, site.pos, site.ref, alt, "aa")
            aaref = rng.choice(AMINO_ACIDS)
            aaalt = rng.choice(AMINO_ACIDS.replace(aaref, ""))
            rows.append(DbnsfpRow(chrom, site.pos, site.ref, alt, aaref, aaalt, list(every)))
    return rows


def planted_dbnsfp_row() -> DbnsfpRow:
    transcripts = [
        TranscriptRow(feature, gene, "missense_variant", str(100 + i), "A/G", f"GENE{i + 1}")
        for i, (feature, gene) in enumerate(zip(PLANTED_TRANSCRIPTS, PLANTED_GENES))
    ]
    return DbnsfpRow(PLANTED_ROW["chr"], PLANTED_ROW["pos"], PLANTED_ROW["ref"], PLANTED_ROW["alt"], "A", "G", transcripts)


def transcript_value(column: str, i: int, transcript: TranscriptRow, row: DbnsfpRow, rng: random.Random) -> str:
    """The fabricated value of a per-transcript column for the `i`th listed transcript, drawn
    from the row's generator in column order."""
    if column == "aapos":
        return transcript.protein_position if transcript.protein_position.isdigit() else str(rng.randrange(1, 900))
    if column == "genename":
        return transcript.symbol or transcript.gene or "."
    if column == "Ensembl_geneid":
        return transcript.gene or "."
    if column == "Ensembl_transcriptid":
        return transcript.feature
    if column == "Ensembl_proteinid":
        return f"ENSP{rng.randrange(100000000, 999999999):011d}"
    if column == "Uniprot_acc":
        return f"{rng.choice('OPQ')}{rng.randrange(10000, 99999)}"
    if column == "HGVSc_VEP":
        return f"c.{rng.randrange(1, 3000)}{row.ref}>{row.alt}"
    if column == "HGVSp_VEP":
        return f"p.{row.aaref}{rng.randrange(1, 900)}{row.aaalt}"
    if column == "APPRIS":
        return rng.choice(("principal1", "alternative1", "."))
    if column == "GENCODE_basic":
        return rng.choice("YN")
    if column == "TSL":
        return rng.choice(("1", "2", "3", "5", "NA"))
    if column == "VEP_canonical":
        return "YES" if i == 0 else "NO"
    if column in ("SIFT_score", "SIFT4G_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "VEST4_score", "REVEL_score", "MVP_score"):
        return f"{rng.uniform(0, 1):.3f}"
    if column in ("SIFT_pred", "SIFT4G_pred", "FATHMM_pred"):
        return rng.choice("DT")
    if column in ("Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred"):
        return rng.choice("DPB")
    if column == "FATHMM_score":
        return f"{rng.uniform(-10, 5):.2f}"
    if column == "PROVEAN_score":
        return f"{rng.uniform(-12, 2):.2f}"
    if column == "PROVEAN_pred":
        return rng.choice("DN")
    raise fail(f"no filler for the per-transcript column {column}")


# The dbNSFP columns that carry one fabricated value for the whole row; every column outside
# TRANSCRIPT_FILLERS and this tuple holds `.`.
SINGLE_FILLERS = (
    "rs_dbSNP", "SIFT_converted_rankscore", "SIFT4G_converted_rankscore", "Polyphen2_HDIV_rankscore",
    "Polyphen2_HVAR_rankscore", "MutationTaster_score", "MutationTaster_pred", "MetaSVM_pred",
    "REVEL_rankscore", "CADD_raw", "CADD_raw_rankscore", "CADD_phred", "GERP++_NR", "GERP++_RS",
    "GERP++_RS_rankscore", "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore",
    "clinvar_clnsig", "clinvar_trait", "Interpro_domain",
)
CLINVAR_TRAITS = ("Cardiomyopathy", "not_provided", "Hereditary_cancer-predisposing_syndrome", "Inborn_genetic_diseases", "not_specified")
INTERPRO_DOMAINS = ("Protein kinase domain", "Zinc finger C2H2-type", "Immunoglobulin-like domain", "P-loop NTPase", "EGF-like domain", "WD40 repeat")


def single_value(column: str, rng: random.Random) -> str:
    """The fabricated value of a whole-row column, drawn from the row's generator in column
    order; `MutationTaster_score` and `MutationTaster_pred` are `;`-joined lists that follow no
    transcript, `clinvar_trait` a `|`-joined list and `Interpro_domain` groups joined by `;` of
    names joined by `|`."""
    if column == "rs_dbSNP":
        number = rng.randrange(1000000, 999999999)
        return f"rs{number}" if rng.random() < 0.6 else "."
    if column.endswith("rankscore"):
        return f"{rng.uniform(0, 1):.5f}"
    if column == "MutationTaster_score":
        return ";".join(f"{rng.uniform(0, 1):.6f}".rstrip("0").rstrip(".") or "0" for _ in range(rng.randrange(1, 4)))
    if column == "MutationTaster_pred":
        return ";".join(rng.choice("ADNP") for _ in range(rng.randrange(1, 4)))
    if column == "MetaSVM_pred":
        return rng.choice("DT")
    if column == "CADD_raw":
        return f"{rng.uniform(-3, 10):.6f}"
    if column == "CADD_phred":
        return f"{rng.uniform(0, 45):.1f}"
    if column == "GERP++_NR":
        return f"{rng.uniform(0, 10):.2f}"
    if column == "GERP++_RS":
        return f"{rng.uniform(-12.3, 6.17):.2f}"
    if column == "phyloP100way_vertebrate":
        return f"{rng.uniform(-5, 10):.3f}"
    if column == "clinvar_clnsig":
        return rng.choice(("Pathogenic", "Likely_pathogenic", "Benign", "Uncertain_significance", ".", "."))
    if column == "clinvar_trait":
        count = rng.randrange(0, 3)
        return "|".join(rng.choice(CLINVAR_TRAITS) for _ in range(count)) or "."
    if column == "Interpro_domain":
        if rng.random() < 0.4:
            return "."
        groups = []
        for _ in range(rng.randrange(1, 3)):
            names = [f"{rng.choice(INTERPRO_DOMAINS)} ({rng.randrange(1, 4)})" for _ in range(rng.randrange(1, 3))]
            groups.append("|".join(names))
        return ";".join(groups)
    raise fail(f"no filler for the column {column}")


def layout_dbnsfp_row(row: DbnsfpRow, prep: Preparation, seed: int, layout: str, planted: bool) -> list[str]:
    """The 458 cells of one row: the corpus coordinates in the prepared assembly's columns, the
    shifted position in the other assembly's, the alleles and the amino acid change, then the
    fillers; a planted row fills every per-transcript column so the detection run sees `;` in
    each of them."""
    rng = row_rng(seed, layout, row.chrom, row.pos, row.ref, row.alt, f"{row.aaref}{row.aaalt}")
    cells = ["."] * len(prep.header)
    other_pos = row.pos if planted else row.pos + OTHER_ASSEMBLY_SHIFT
    cells[prep.chr_col] = row.chrom
    cells[prep.pos_col] = str(row.pos)
    cells[prep.other_chr_col] = row.chrom
    cells[prep.other_pos_col] = str(other_pos)
    index = {name: i for i, name in enumerate(prep.header)}
    cells[index["ref"]] = row.ref
    cells[index["alt"]] = row.alt
    cells[index["aaref"]] = row.aaref
    cells[index["aaalt"]] = row.aaalt
    for column in TRANSCRIPT_FILLERS:
        values = [transcript_value(column, i, transcript, row, rng) for i, transcript in enumerate(row.transcripts)]
        cells[index[column]] = ";".join(values) if values else "."
    for column in SINGLE_FILLERS:
        cells[index[column]] = single_value(column, rng)
    if planted:
        planted_lists = {
            "MutationTaster_score": "1;0.999912",
            "MutationTaster_pred": "D;D",
            "Interpro_domain": "Protein kinase domain (2)|P-loop NTPase (1);EGF-like domain (1)",
            "clinvar_trait": "Cardiomyopathy|not_provided",
        }
        for column, value in planted_lists.items():
            cells[index[column]] = value
    return cells


def dbnsfp_lines(sites: list[Site], reference: dict, prep: Preparation, seed: int, layout: str, misses: set[str]) -> tuple[list[str], dict]:
    """The sorted data lines of a dbNSFP file and the row counts: `own_allele` (rows for an
    input allele with a real amino acid change), `fabricated` (the other rows), `planted` and
    `miss` (sites left without a row)."""
    counts = {"own_allele": 0, "fabricated": 0, "planted": 1, "miss": 0}
    keyed: dict[tuple, list[str]] = {}
    for site in sites:
        if f"{site.chrom}:{site.pos}" in misses:
            counts["miss"] += 1
            continue
        for row in dbnsfp_rows_for_site(site, reference, layout, seed):
            counts["own_allele" if row.from_output else "fabricated"] += 1
            key = (row.chrom, row.pos, row.ref, row.alt, row.aaref, row.aaalt)
            keyed[key] = layout_dbnsfp_row(row, prep, seed, layout, planted=False)
    planted = planted_dbnsfp_row()
    planted_key = (planted.chrom, planted.pos, planted.ref, planted.alt, planted.aaref, planted.aaalt)
    keyed[planted_key] = layout_dbnsfp_row(planted, prep, seed, layout, planted=True)
    return sorted_lines(keyed, prep), counts


def sorted_lines(keyed: dict[tuple, list[str]], prep: Preparation) -> list[str]:
    """The rows as the generator sorts them: by the prepared assembly's
    chromosome column as text, then its position column as a number, then the row's key."""
    ordered = sorted(keyed.items(), key=lambda item: (item[1][prep.chr_col], int(item[1][prep.pos_col]), item[0][2:]))
    return ["\t".join(cells) for _, cells in ordered]


def dbscsnv_row(site: Site, alt: str, prep: Preparation, seed: int, layout: str, reference: dict) -> list[str]:
    chrom = spell(site.chrom, layout)
    rng = row_rng(seed, layout, chrom, site.pos, site.ref, alt)
    genes: list[str] = []
    for allele in BASES:
        for row in reference.get((site.chrom, site.pos, allele), []):
            if row.gene and row.gene not in genes:
                genes.append(row.gene)
    genes.sort()
    cells = ["."] * len(prep.header)
    cells[prep.chr_col] = chrom
    cells[prep.pos_col] = str(site.pos)
    cells[prep.other_chr_col] = chrom
    cells[prep.other_pos_col] = str(site.pos + OTHER_ASSEMBLY_SHIFT)
    index = {name: i for i, name in enumerate(prep.header)}
    cells[index["ref"]] = site.ref
    cells[index["alt"]] = alt
    cells[index["RefSeq?"]] = rng.choice("YN")
    cells[index["Ensembl?"]] = "Y" if genes else "N"
    region = rng.choice(("exon", "intron", "exon", "intron"))
    cells[index["RefSeq_region"]] = region
    cells[index["RefSeq_gene"]] = f"GENE{rng.randrange(1, 999)}"
    cells[index["RefSeq_functional_consequence"]] = rng.choice(("synonymous", "missense", "intron"))
    cells[index["RefSeq_id_c.change_p.change"]] = f"NM_{rng.randrange(100000, 999999)}:c.{rng.randrange(1, 3000)}{site.ref}>{alt}"
    cells[index["Ensembl_region"]] = region
    cells[index["Ensembl_gene"]] = ";".join(genes) if genes else "."
    cells[index["Ensembl_functional_consequence"]] = rng.choice(("synonymous", "missense", "intron"))
    cells[index["Ensembl_id_c.change_p.change"]] = f"ENST{rng.randrange(100000000, 999999999):011d}:c.{rng.randrange(1, 3000)}{site.ref}>{alt}"
    scores = [f"{rng.uniform(0, 1):.6f}", f"{rng.uniform(0, 1):.6f}"]
    if rng.random() < 0.15:
        scores[rng.randrange(2)] = "."
    cells[index["ada_score"]], cells[index["rf_score"]] = scores
    return cells


def dbscsnv_lines(sites: list[Site], reference: dict, prep: Preparation, seed: int, layout: str, misses: set[str]) -> tuple[list[str], dict]:
    """The sorted data lines of a dbscSNV file and the row counts: `snv` (one per alternate
    base per site) and `miss`."""
    counts = {"snv": 0, "miss": 0}
    keyed: dict[tuple, list[str]] = {}
    for site in sites:
        if f"{site.chrom}:{site.pos}" in misses:
            counts["miss"] += 1
            continue
        for alt in BASES:
            if alt == site.ref:
                continue
            cells = dbscsnv_row(site, alt, prep, seed, layout, reference)
            keyed[(cells[prep.chr_col], site.pos, site.ref, alt)] = cells
            counts["snv"] += 1
    return sorted_lines(keyed, prep), counts


def tool_version(tool: str) -> str:
    """The first line `<tool> --version` prints."""
    if shutil.which(tool) is None:
        raise fail(f"{tool} is not on PATH; bgzip and tabix write the file and its index")
    out = subprocess.run([tool, "--version"], capture_output=True, text=True, check=False)
    lines = (out.stdout or out.stderr).strip().splitlines()
    if out.returncode != 0 or not lines:
        raise fail(f"{tool} --version exited {out.returncode} and reported no version")
    return lines[0]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fixture(lines: list[str], out: Path, prep: Preparation) -> None:
    """Writes the header and `lines` through bgzip to `out`, each line ended as the layout
    ends it, and the layout's tabix index to `out.tbi`."""
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as scratch_name:
        scratch = Path(scratch_name)
        plain = scratch / "fixture.tsv"
        body = "".join(f"{line}{prep.line_end}" for line in ["\t".join(prep.header), *lines])
        plain.write_bytes(body.encode("utf-8"))
        compressed = scratch / out.name
        with compressed.open("wb") as handle:
            bgzip_run = subprocess.run(["bgzip", "-c", str(plain)], stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
        if bgzip_run.returncode != 0:
            raise fail(f"bgzip -c exited {bgzip_run.returncode}: {bgzip_run.stderr.strip()}")
        tabix_run = subprocess.run(["tabix", "-f", *prep.tabix, str(compressed)], capture_output=True, text=True, check=False)
        if tabix_run.returncode != 0:
            raise fail(f"tabix exited {tabix_run.returncode}: {tabix_run.stderr.strip()}")
        shutil.move(str(compressed), str(out))
        shutil.move(f"{compressed}.tbi", f"{out}.tbi")


def spec_path(out: Path) -> Path:
    return out.with_name(f"{out.name}.spec.json")


def shown_path(path: Path, corpus: Path | None) -> str:
    """A path as the spec records it: relative to the corpus when it lies inside one, else
    its file name alone, so the record names no directory of the machine it was written on."""
    resolved = path.resolve()
    if corpus is not None and resolved.is_relative_to(corpus.resolve()):
        return str(resolved.relative_to(corpus.resolve()))
    return path.name


def write_spec(out: Path, args: argparse.Namespace, prep: Preparation, counts: dict, versions: dict[str, str]) -> dict:
    corpus = args.provenance.parent if args.provenance else None
    spec = {
        "generator": GENERATOR,
        "layout": args.layout,
        "assembly": args.assembly,
        "seed": args.seed,
        "preparation": {
            "index": prep.index_text,
            "coordinate_columns": [prep.header[prep.chr_col].lstrip("#"), prep.header[prep.pos_col]],
            "other_assembly_columns": [prep.header[prep.other_chr_col].lstrip("#"), prep.header[prep.other_pos_col]],
            "other_assembly_shift": OTHER_ASSEMBLY_SHIFT,
            "line_end": "CRLF" if prep.line_end == "\r\n" else "LF",
            "columns": len(prep.header),
        },
        "inputs": {
            "vcf": {shown_path(args.vcf, corpus): sha256(args.vcf)},
            "reference": {shown_path(args.reference, corpus): sha256(args.reference)},
        },
        "misses": sorted(args.miss or []),
        "rows": counts,
        "files": {out.name: sha256(out), f"{out.name}.tbi": sha256(Path(f"{out}.tbi"))},
        "tools": versions,
    }
    if args.layout == "dbnsfp":
        spec["planted_row"] = {**PLANTED_ROW, "transcripts": list(PLANTED_TRANSCRIPTS), "region": f"{DETECTION_REGION[0]}:{DETECTION_REGION[1]}-{DETECTION_REGION[2]}"}
    spec_path(out).write_text(json.dumps(spec, indent=1) + "\n", encoding="utf-8")
    return spec


def update_provenance(provenance: Path, out: Path, spec: dict) -> None:
    """Rewrites `plugin_fixtures.<file name>` in a corpus's provenance.json from the spec, with
    the file paths relative to the corpus, and keeps every other key."""
    record = json.loads(provenance.read_text(encoding="utf-8")) if provenance.is_file() else {}
    if not isinstance(record, dict):
        raise fail(f"{provenance}: the record is not a JSON object")
    fixtures = record.setdefault("plugin_fixtures", {})
    if not isinstance(fixtures, dict):
        raise fail(f"{provenance}: plugin_fixtures is not a JSON object")
    corpus = provenance.resolve().parent
    rel = out.resolve().relative_to(corpus)
    entry = {k: v for k, v in spec.items() if k not in ("files", "inputs")}
    entry["spec"] = str(spec_path(out).resolve().relative_to(corpus))
    entry["files"] = {
        str(rel): spec["files"][out.name],
        f"{rel}.tbi": spec["files"][f"{out.name}.tbi"],
        entry["spec"]: sha256(spec_path(out)),
    }
    fixtures[out.name] = entry
    tmp = provenance.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    tmp.replace(provenance)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layout", required=True, choices=LAYOUTS, help="the file layout written")
    ap.add_argument("--assembly", required=True, choices=ASSEMBLIES, help="the assembly the file is prepared for")
    ap.add_argument("--seed", required=True, type=int, help="the seed every fabricated value derives from")
    ap.add_argument("--vcf", required=True, type=Path, help="the corpus VCF whose single-nucleotide sites get rows")
    ap.add_argument("--reference", required=True, type=Path, help="a default- or tab-format VEP output of the VCF (plain or .gz)")
    ap.add_argument("--out", required=True, type=Path, help="the file written (its .tbi and .spec.json sit beside it)")
    ap.add_argument("--provenance", type=Path, help="a corpus provenance.json to record the fixture in")
    ap.add_argument("--miss", action="append", metavar="CHR:POS", help="a site (the input's spelling) left without a row (repeatable)")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output")
    args = ap.parse_args(argv)

    if args.layout in SPLICEAI_AND_CADD:
        raise fail(f"layout {args.layout} is registered but has no writer here")
    for path in (args.vcf, args.reference):
        if not path.is_file():
            raise fail(f"input file not found: {path}")
    for miss in args.miss or []:
        if not SITE.fullmatch(miss):
            raise fail(f"--miss names a site as chr:pos, got {miss!r}")
    if args.out.exists() and not args.force:
        raise fail(f"{args.out} exists; pass --force to overwrite it")
    if args.provenance and not args.out.resolve().is_relative_to(args.provenance.resolve().parent):
        raise fail(f"--out {args.out} lies outside the corpus of --provenance {args.provenance}")
    versions = {tool: tool_version(tool) for tool in TOOLS}

    prep = preparation(args.layout, args.assembly)
    sites = read_sites(args.vcf)
    reference = read_reference(args.reference)
    misses = set(args.miss or [])
    if args.layout == "dbnsfp":
        lines, counts = dbnsfp_lines(sites, reference, prep, args.seed, args.layout, misses)
    else:
        lines, counts = dbscsnv_lines(sites, reference, prep, args.seed, args.layout, misses)
    if not any(v for k, v in counts.items() if k not in ("miss", "planted")):
        raise fail(f"{args.vcf}: no single-nucleotide site to write a row for (or every site missed)")
    write_fixture(lines, args.out, prep)
    spec = write_spec(args.out, args, prep, counts, versions)
    if args.provenance:
        update_provenance(args.provenance, args.out, spec)
    summary = ", ".join(f"{k} {v}" for k, v in counts.items())
    print(f"wrote {args.out} and its .tbi and .spec.json: {len(lines)} rows ({summary}), seed {args.seed}, {prep.index_text}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
