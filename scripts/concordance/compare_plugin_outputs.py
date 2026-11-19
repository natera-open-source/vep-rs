#!/usr/bin/env python3
"""Compare plugin output fields between Perl VEP and Rust VEP.

Supports batch mode (multiple plugins), per-directory processing (all matching
file pairs), and JSON output for integration with the concordance harness.

Usage:
    # Single plugin, single file pair
    python3 compare_plugin_outputs.py \
        --perl-output perl_out.txt --rust-output rust_out.txt \
        --plugins CADD

    # Multiple plugins, directory mode
    python3 compare_plugin_outputs.py \
        --perl-dir tmp/concordance/perl --rust-dir tmp/concordance/rust \
        --plugins CADD,REVEL,gnomADc \
        --json-output tmp/concordance/reports/plugin_concordance.json

    # All known plugins with custom tolerance
    python3 compare_plugin_outputs.py \
        --perl-dir perl/ --rust-dir rust/ \
        --plugins all --tolerance 1e-4
"""

import argparse
import functools
import gzip
import json
import re
import sys
from pathlib import Path
from urllib.parse import unquote

# Plugin field registry

# The columns of a plugin's data file that carry the name of one of Ensembl VEP's own
# Extra keys, each with the flags that turn that key on. No entry of PLUGIN_FIELDS
# names them: the name is the column's in a run that loads the plugin and sets none
# of the flags, and both the column's and the key's in a run that sets one.
CORE_NAMED_COLUMNS = {
    "dbNSFP": {
        "APPRIS": ["--appris", "--everything"],
        "TSL": ["--tsl", "--everything"],
    },
}

# The columns of the golden corpora's dbNSFP fixtures as the corpora's references
# declare them, without dbNSFP's CORE_NAMED_COLUMNS.
_DBNSFP_COLUMNS = """
1000Gp3_AC 1000Gp3_AF 1000Gp3_AFR_AC 1000Gp3_AFR_AF 1000Gp3_AMR_AC 1000Gp3_AMR_AF
1000Gp3_EAS_AC 1000Gp3_EAS_AF 1000Gp3_EUR_AC 1000Gp3_EUR_AF 1000Gp3_SAS_AC
1000Gp3_SAS_AF ALFA_African_AC ALFA_African_AF ALFA_African_AN
ALFA_African_American_AC ALFA_African_American_AF ALFA_African_American_AN
ALFA_African_Others_AC ALFA_African_Others_AF ALFA_African_Others_AN ALFA_Asian_AC
ALFA_Asian_AF ALFA_Asian_AN ALFA_East_Asian_AC ALFA_East_Asian_AF ALFA_East_Asian_AN
ALFA_European_AC ALFA_European_AF ALFA_European_AN ALFA_Latin_American_1_AC
ALFA_Latin_American_1_AF ALFA_Latin_American_1_AN ALFA_Latin_American_2_AC
ALFA_Latin_American_2_AF ALFA_Latin_American_2_AN ALFA_Other_AC ALFA_Other_AF
ALFA_Other_AN ALFA_Other_Asian_AC ALFA_Other_Asian_AF ALFA_Other_Asian_AN
ALFA_South_Asian_AC ALFA_South_Asian_AF ALFA_South_Asian_AN ALFA_Total_AC
ALFA_Total_AF ALFA_Total_AN ALSPAC_AC ALSPAC_AF Aloft_Confidence
Aloft_Fraction_transcripts_affected Aloft_pred Aloft_prob_Dominant
Aloft_prob_Recessive Aloft_prob_Tolerant AlphaMissense_pred AlphaMissense_rankscore
AlphaMissense_score AltaiNeandertal Ancestral_allele BayesDel_addAF_pred
BayesDel_addAF_rankscore BayesDel_addAF_score BayesDel_noAF_pred
BayesDel_noAF_rankscore BayesDel_noAF_score CADD_phred CADD_phred_hg19 CADD_raw
CADD_raw_hg19 CADD_raw_rankscore CADD_raw_rankscore_hg19 ChagyrskayaNeandertal
ClinPred_pred ClinPred_rankscore ClinPred_score DANN_rankscore DANN_score
DEOGEN2_pred DEOGEN2_rankscore DEOGEN2_score Denisova ESM1b_pred ESM1b_rankscore
ESM1b_score ESP6500_AA_AC ESP6500_AA_AF ESP6500_EA_AC ESP6500_EA_AF EVE_Class10_pred
EVE_Class20_pred EVE_Class25_pred EVE_Class30_pred EVE_Class40_pred EVE_Class50_pred
EVE_Class60_pred EVE_Class70_pred EVE_Class75_pred EVE_Class80_pred EVE_Class90_pred
EVE_rankscore EVE_score Eigen-PC-phred_coding Eigen-PC-raw_coding
Eigen-PC-raw_coding_rankscore Eigen-phred_coding Eigen-raw_coding
Eigen-raw_coding_rankscore Ensembl_geneid Ensembl_proteinid Ensembl_transcriptid
ExAC_AC ExAC_AF ExAC_AFR_AC ExAC_AFR_AF ExAC_AMR_AC ExAC_AMR_AF ExAC_Adj_AC
ExAC_Adj_AF ExAC_EAS_AC ExAC_EAS_AF ExAC_FIN_AC ExAC_FIN_AF ExAC_NFE_AC ExAC_NFE_AF
ExAC_SAS_AC ExAC_SAS_AF ExAC_nonTCGA_AC ExAC_nonTCGA_AF ExAC_nonTCGA_AFR_AC
ExAC_nonTCGA_AFR_AF ExAC_nonTCGA_AMR_AC ExAC_nonTCGA_AMR_AF ExAC_nonTCGA_Adj_AC
ExAC_nonTCGA_Adj_AF ExAC_nonTCGA_EAS_AC ExAC_nonTCGA_EAS_AF ExAC_nonTCGA_FIN_AC
ExAC_nonTCGA_FIN_AF ExAC_nonTCGA_NFE_AC ExAC_nonTCGA_NFE_AF ExAC_nonTCGA_SAS_AC
ExAC_nonTCGA_SAS_AF ExAC_nonpsych_AC ExAC_nonpsych_AF ExAC_nonpsych_AFR_AC
ExAC_nonpsych_AFR_AF ExAC_nonpsych_AMR_AC ExAC_nonpsych_AMR_AF ExAC_nonpsych_Adj_AC
ExAC_nonpsych_Adj_AF ExAC_nonpsych_EAS_AC ExAC_nonpsych_EAS_AF ExAC_nonpsych_FIN_AC
ExAC_nonpsych_FIN_AF ExAC_nonpsych_NFE_AC ExAC_nonpsych_NFE_AF ExAC_nonpsych_SAS_AC
ExAC_nonpsych_SAS_AF FATHMM_converted_rankscore FATHMM_pred FATHMM_score
GENCODE_basic GERP++_NR GERP++_RS GERP++_RS_rankscore GERP_91_mammals
GERP_91_mammals_rankscore GM12878_confidence_value GM12878_fitCons_rankscore
GM12878_fitCons_score GTEx_V8_eQTL_gene GTEx_V8_eQTL_tissue GTEx_V8_sQTL_gene
GTEx_V8_sQTL_tissue GenoCanyon_rankscore GenoCanyon_score Geuvadis_eQTL_target_gene
H1-hESC_confidence_value H1-hESC_fitCons_rankscore H1-hESC_fitCons_score
HGVSc_ANNOVAR HGVSc_VEP HGVSc_snpEff HGVSp_ANNOVAR HGVSp_VEP HGVSp_snpEff
HUVEC_confidence_value HUVEC_fitCons_rankscore HUVEC_fitCons_score Interpro_domain
LINSIGHT LINSIGHT_rankscore LIST-S2_pred LIST-S2_rankscore LIST-S2_score LRT_Omega
LRT_converted_rankscore LRT_pred LRT_score M-CAP_pred M-CAP_rankscore M-CAP_score
MPC_rankscore MPC_score MVP_rankscore MVP_score MetaLR_pred MetaLR_rankscore
MetaLR_score MetaRNN_pred MetaRNN_rankscore MetaRNN_score MetaSVM_pred
MetaSVM_rankscore MetaSVM_score MutFormer_rankscore MutFormer_score MutPred_AAchange
MutPred_Top5features MutPred_protID MutPred_rankscore MutPred_score
MutScore_rankscore MutScore_score MutationAssessor_pred MutationAssessor_rankscore
MutationAssessor_score MutationTaster_AAE MutationTaster_converted_rankscore
MutationTaster_model MutationTaster_pred MutationTaster_score PHACTboost_rankscore
PHACTboost_score PROVEAN_converted_rankscore PROVEAN_pred PROVEAN_score
Polyphen2_HDIV_pred Polyphen2_HDIV_rankscore Polyphen2_HDIV_score
Polyphen2_HVAR_pred Polyphen2_HVAR_rankscore Polyphen2_HVAR_score PrimateAI_pred
PrimateAI_rankscore PrimateAI_score REVEL_rankscore REVEL_score Reliability_index
SIFT4G_converted_rankscore SIFT4G_pred SIFT4G_score SIFT_converted_rankscore
SIFT_pred SIFT_score SiPhy_29way_logOdds SiPhy_29way_logOdds_rankscore
SiPhy_29way_pi TWINSUK_AC TWINSUK_AF UK10K_AC UK10K_AF Uniprot_acc Uniprot_entry
VARITY_ER_LOO_rankscore VARITY_ER_LOO_score VARITY_ER_rankscore VARITY_ER_score
VARITY_R_LOO_rankscore VARITY_R_LOO_score VARITY_R_rankscore VARITY_R_score
VEP_canonical VEST4_rankscore VEST4_score VindijiaNeandertal aaalt aapos aaref alt
bStatistic bStatistic_converted_rankscore cds_strand chr clinvar_MedGen_id
clinvar_OMIM_id clinvar_Orphanet_id clinvar_clnsig clinvar_hgvs clinvar_id
clinvar_review clinvar_trait clinvar_var_source codon_degeneracy codonpos
eQTLGen_cis_or_trans eQTLGen_gene_id eQTLGen_gene_symbol eQTLGen_snp_id
fathmm-MKL_coding_group fathmm-MKL_coding_pred fathmm-MKL_coding_rankscore
fathmm-MKL_coding_score fathmm-XF_coding_pred fathmm-XF_coding_rankscore
fathmm-XF_coding_score gMVP_rankscore gMVP_score genename gnomAD_exomes_AC
gnomAD_exomes_AF gnomAD_exomes_AFR_AC gnomAD_exomes_AFR_AF gnomAD_exomes_AFR_AN
gnomAD_exomes_AFR_nhomalt gnomAD_exomes_AMR_AC gnomAD_exomes_AMR_AF
gnomAD_exomes_AMR_AN gnomAD_exomes_AMR_nhomalt gnomAD_exomes_AN gnomAD_exomes_ASJ_AC
gnomAD_exomes_ASJ_AF gnomAD_exomes_ASJ_AN gnomAD_exomes_ASJ_nhomalt
gnomAD_exomes_EAS_AC gnomAD_exomes_EAS_AF gnomAD_exomes_EAS_AN
gnomAD_exomes_EAS_nhomalt gnomAD_exomes_FIN_AC gnomAD_exomes_FIN_AF
gnomAD_exomes_FIN_AN gnomAD_exomes_FIN_nhomalt gnomAD_exomes_MID_AC
gnomAD_exomes_MID_AF gnomAD_exomes_MID_AN gnomAD_exomes_MID_nhomalt
gnomAD_exomes_NFE_AC gnomAD_exomes_NFE_AF gnomAD_exomes_NFE_AN
gnomAD_exomes_NFE_nhomalt gnomAD_exomes_POPMAX_AC gnomAD_exomes_POPMAX_AF
gnomAD_exomes_POPMAX_AN gnomAD_exomes_POPMAX_nhomalt gnomAD_exomes_SAS_AC
gnomAD_exomes_SAS_AF gnomAD_exomes_SAS_AN gnomAD_exomes_SAS_nhomalt
gnomAD_exomes_flag gnomAD_exomes_nhomalt gnomAD_exomes_non_ukb_AC
gnomAD_exomes_non_ukb_AF gnomAD_exomes_non_ukb_AFR_AC gnomAD_exomes_non_ukb_AFR_AF
gnomAD_exomes_non_ukb_AFR_AN gnomAD_exomes_non_ukb_AFR_nhomalt
gnomAD_exomes_non_ukb_AMR_AC gnomAD_exomes_non_ukb_AMR_AF
gnomAD_exomes_non_ukb_AMR_AN gnomAD_exomes_non_ukb_AMR_nhomalt
gnomAD_exomes_non_ukb_AN gnomAD_exomes_non_ukb_ASJ_AC gnomAD_exomes_non_ukb_ASJ_AF
gnomAD_exomes_non_ukb_ASJ_AN gnomAD_exomes_non_ukb_ASJ_nhomalt
gnomAD_exomes_non_ukb_EAS_AC gnomAD_exomes_non_ukb_EAS_AF
gnomAD_exomes_non_ukb_EAS_AN gnomAD_exomes_non_ukb_EAS_nhomalt
gnomAD_exomes_non_ukb_FIN_AC gnomAD_exomes_non_ukb_FIN_AF
gnomAD_exomes_non_ukb_FIN_AN gnomAD_exomes_non_ukb_FIN_nhomalt
gnomAD_exomes_non_ukb_MID_AC gnomAD_exomes_non_ukb_MID_AF
gnomAD_exomes_non_ukb_MID_AN gnomAD_exomes_non_ukb_MID_nhomalt
gnomAD_exomes_non_ukb_NFE_AC gnomAD_exomes_non_ukb_NFE_AF
gnomAD_exomes_non_ukb_NFE_AN gnomAD_exomes_non_ukb_NFE_nhomalt
gnomAD_exomes_non_ukb_SAS_AC gnomAD_exomes_non_ukb_SAS_AF
gnomAD_exomes_non_ukb_SAS_AN gnomAD_exomes_non_ukb_SAS_nhomalt
gnomAD_exomes_non_ukb_nhomalt gnomAD_genomes_AC gnomAD_genomes_AF
gnomAD_genomes_AFR_AC gnomAD_genomes_AFR_AF gnomAD_genomes_AFR_AN
gnomAD_genomes_AFR_nhomalt gnomAD_genomes_AMI_AC gnomAD_genomes_AMI_AF
gnomAD_genomes_AMI_AN gnomAD_genomes_AMI_nhomalt gnomAD_genomes_AMR_AC
gnomAD_genomes_AMR_AF gnomAD_genomes_AMR_AN gnomAD_genomes_AMR_nhomalt
gnomAD_genomes_AN gnomAD_genomes_ASJ_AC gnomAD_genomes_ASJ_AF gnomAD_genomes_ASJ_AN
gnomAD_genomes_ASJ_nhomalt gnomAD_genomes_EAS_AC gnomAD_genomes_EAS_AF
gnomAD_genomes_EAS_AN gnomAD_genomes_EAS_nhomalt gnomAD_genomes_FIN_AC
gnomAD_genomes_FIN_AF gnomAD_genomes_FIN_AN gnomAD_genomes_FIN_nhomalt
gnomAD_genomes_MID_AC gnomAD_genomes_MID_AF gnomAD_genomes_MID_AN
gnomAD_genomes_MID_nhomalt gnomAD_genomes_NFE_AC gnomAD_genomes_NFE_AF
gnomAD_genomes_NFE_AN gnomAD_genomes_NFE_nhomalt gnomAD_genomes_POPMAX_AC
gnomAD_genomes_POPMAX_AF gnomAD_genomes_POPMAX_AN gnomAD_genomes_POPMAX_nhomalt
gnomAD_genomes_SAS_AC gnomAD_genomes_SAS_AF gnomAD_genomes_SAS_AN
gnomAD_genomes_SAS_nhomalt gnomAD_genomes_flag gnomAD_genomes_nhomalt hg18_chr
hg18_pos(1-based) hg19_chr hg19_pos(1-based) integrated_confidence_value
integrated_fitCons_rankscore integrated_fitCons_score phastCons100way_vertebrate
phastCons100way_vertebrate_rankscore phastCons17way_primate
phastCons17way_primate_rankscore phastCons470way_mammalian
phastCons470way_mammalian_rankscore phyloP100way_vertebrate
phyloP100way_vertebrate_rankscore phyloP17way_primate phyloP17way_primate_rankscore
phyloP470way_mammalian phyloP470way_mammalian_rankscore pos(1-based) ref refcodon
rs_dbSNP
""".split()  # noqa: SIM905

# An entry is a field name as Ensembl VEP writes it, or a pattern: `<N>` stands for a
# run of digits, any other `<word>` for a name the run chooses.
PLUGIN_FIELDS = {
    "HGVS": ["HGVSc", "HGVSp"],
    "CADD": ["CADD_PHRED", "CADD_RAW"],
    "REVEL": ["REVEL"],
    # One joined field in the default and tab formats, nine fields in the VCF format
    # and under split_output, the flag under cutoff; a field absent from both
    # outputs is not compared.
    "SpliceAI": [
        "SpliceAI_pred",
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
        "SpliceAI_pred_DP_AG",
        "SpliceAI_pred_DP_AL",
        "SpliceAI_pred_DP_DG",
        "SpliceAI_pred_DP_DL",
        "SpliceAI_pred_SYMBOL",
        "SpliceAI_cutoff",
    ],
    # The coverage file's columns choose the keys: a median or an approximate median,
    # a total depth where the file has one, and one field per depth threshold.
    "gnomADc": [
        "gnomAD_mean_cov",
        "gnomAD_median_cov",
        "gnomAD_median_approx_cov",
        "gnomAD_total_DP_cov",
        "gnomAD_<N>x_cov",
    ],
    "AlphaMissense": [
        "am_class",
        "am_pathogenicity",
        "am_genome",
        "am_uniprot_id",
        "am_transcript_id",
        "am_protein_variant",
    ],
    "AVI": ["AVI_PHRED", "AVI_RAW", "AVI_LINK"],
    "GWAS": [
        "GWAS_accessions",
        "GWAS_associated_gene",
        "GWAS_beta_coef",
        "GWAS_odds_ratio",
        "GWAS_p_value",
        "GWAS_pmid",
        "GWAS_risk_allele",
        "GWAS_study",
    ],
    "Blosum62": ["BLOSUM62"],
    "Downstream": ["DownstreamProtein", "ProteinLengthChange"],
    "HGVSIntronOffset": ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"],
    "NearestExonJB": ["NearestExonJB"],
    "TSSDistance": ["TSSDistance"],
    # The first five are the default columns; `cols=` selects among all of them.
    "OpenTargets": [
        "OpenTargets_gwasLocusToGeneScore",
        "OpenTargets_gwasGeneId",
        "OpenTargets_gwasDiseases",
        "OpenTargets_qtlGeneId",
        "OpenTargets_qtlBiosampleName",
        "OpenTargets_beta",
        "OpenTargets_finemappingMethod",
        "OpenTargets_is95CredibleSet",
        "OpenTargets_is99CredibleSet",
        "OpenTargets_isLead",
        "OpenTargets_pValueExponent",
        "OpenTargets_pValueMantissa",
        "OpenTargets_posteriorProbability",
        "OpenTargets_studyId",
        "OpenTargets_studyLocusId",
        "OpenTargets_studyType",
        "OpenTargets_variantId",
    ],
    # dbNSFP writes the selected columns under their own names; a run with
    # output_prefix=1 prefixes them dbNSFP_. NoSuchColumn is a name the fixtures lack,
    # which Ensembl VEP writes as the run gives it.
    "dbNSFP": ["dbNSFP_<column>", *_DBNSFP_COLUMNS, "NoSuchColumn"],
    "dbscSNV": ["ada_score", "rf_score"],
    "LoFtool": ["LoFtool"],
    # pLI writes one of the two: the gene field by default, the transcript
    # field in its transcript mode.
    "pLI": ["pLI_gene_value", "pLI_transcript_value"],
    # LoFTEE emits the same field names as the Perl LoF plugin, so vep-rs
    # (--plugin LoFTEE) and Perl (--plugin LoF) compare directly. Compare the
    # LoF (HC/LC) verdict and LoF_filter per-field; LoF_flags/LoF_info are
    # informational, and GERP_DIST in LoF_info differs under a GERP bigWig
    # with zoom levels (use --fields LoF to isolate).
    "LoFTEE": ["LoF", "LoF_filter", "LoF_flags", "LoF_info"],
    # The three groups of MaxEntScan.pm: the canonical-window scores, the
    # sliding-window scores (SWA) and the nearest canonical sites (NCSS), with the
    # sequences, positions and contexts `verbose` adds.
    "MaxEntScan": [
        "MaxEntScan_ref",
        "MaxEntScan_alt",
        "MaxEntScan_diff",
        "MaxEntScan_ref_seq",
        "MaxEntScan_alt_seq",
        "MES-SWA_donor_ref",
        "MES-SWA_donor_alt",
        "MES-SWA_donor_ref_comp",
        "MES-SWA_donor_diff",
        "MES-SWA_donor_ref_seq",
        "MES-SWA_donor_ref_frame",
        "MES-SWA_donor_ref_context",
        "MES-SWA_donor_alt_seq",
        "MES-SWA_donor_alt_frame",
        "MES-SWA_donor_alt_context",
        "MES-SWA_donor_ref_comp_seq",
        "MES-SWA_acceptor_ref",
        "MES-SWA_acceptor_alt",
        "MES-SWA_acceptor_ref_comp",
        "MES-SWA_acceptor_diff",
        "MES-SWA_acceptor_ref_seq",
        "MES-SWA_acceptor_ref_frame",
        "MES-SWA_acceptor_ref_context",
        "MES-SWA_acceptor_alt_seq",
        "MES-SWA_acceptor_alt_frame",
        "MES-SWA_acceptor_alt_context",
        "MES-SWA_acceptor_ref_comp_seq",
        "MES-NCSS_upstream_acceptor",
        "MES-NCSS_upstream_donor",
        "MES-NCSS_downstream_acceptor",
        "MES-NCSS_downstream_donor",
        "MES-NCSS_upstream_acceptor_seq",
        "MES-NCSS_upstream_donor_seq",
        "MES-NCSS_downstream_acceptor_seq",
        "MES-NCSS_downstream_donor_seq",
    ],
    "AncestralAllele": ["AA"],
    "Conservation": ["Conservation"],
    "BayesDel": ["BayesDel"],
    "ClinPred": ["ClinPred"],
    "DosageSensitivity": ["pHaplo", "pTriplo"],
    # The file's INFO fields: five haplogroup fields by default, any of the others
    # under `fields=`. Each is named, because gnomADc's fields and dbNSFP's gnomAD
    # columns share the prefix.
    "gnomADMt": [
        "gnomAD_hap_AC_het",
        "gnomAD_hap_AC_hom",
        "gnomAD_hap_AF_het",
        "gnomAD_hap_AF_hom",
        "gnomAD_hap_AN",
        "gnomAD_AC_het",
        "gnomAD_AC_hom",
        "gnomAD_AF_het",
        "gnomAD_AF_hom",
        "gnomAD_AN",
        "gnomAD_age_hist_het_bin_freq",
        "gnomAD_age_hist_het_n_larger",
        "gnomAD_age_hist_het_n_smaller",
        "gnomAD_age_hist_hom_bin_freq",
        "gnomAD_age_hist_hom_n_larger",
        "gnomAD_age_hist_hom_n_smaller",
        "gnomAD_base_qual_hist",
        "gnomAD_common_low_heteroplasmy",
        "gnomAD_contamination_hist",
        "gnomAD_dp_hist_all_bin_freq",
        "gnomAD_dp_hist_all_n_larger",
        "gnomAD_dp_hist_alt_bin_freq",
        "gnomAD_dp_hist_alt_n_larger",
        "gnomAD_dp_mean",
        "gnomAD_excluded_AC",
        "gnomAD_faf_hapmax_hom",
        "gnomAD_filters",
        "gnomAD_hap_defining_variant",
        "gnomAD_hap_faf_hom",
        "gnomAD_hap_hl_hist",
        "gnomAD_hapmax_AF_het",
        "gnomAD_hapmax_AF_hom",
        "gnomAD_heteroplasmy_below_min_het_threshold_hist",
        "gnomAD_hl_hist",
        "gnomAD_max_hl",
        "gnomAD_mitotip_score",
        "gnomAD_mitotip_trna_prediction",
        "gnomAD_mq_mean",
        "gnomAD_pon_ml_probability_of_pathogenicity",
        "gnomAD_pon_mt_trna_prediction",
        "gnomAD_pop_AC_het",
        "gnomAD_pop_AC_hom",
        "gnomAD_pop_AF_het",
        "gnomAD_pop_AF_hom",
        "gnomAD_pop_AN",
        "gnomAD_pop_hl_hist",
        "gnomAD_position_hist",
        "gnomAD_strand_bias_hist",
        "gnomAD_tlod_mean",
        "gnomAD_variant_collapsed",
        "gnomAD_vep",
        "gnomAD_weak_evidence_hist",
    ],
    "LOEUF": ["LOEUF"],
    "NMD": ["NMD"],
    "PrimateAI": ["PrimateAI"],
    "SpliceVault": [
        "SpliceVault_SpliceAI_delta",
        "SpliceVault_out_of_frame_events",
        "SpliceVault_site_max_depth",
        "SpliceVault_site_pos",
        "SpliceVault_site_sample_count",
        "SpliceVault_site_type",
        "SpliceVault_top_events",
    ],
    # One field per selected column of the run's own table, under the
    # transcript_annotator_ prefix by default; `prefix=0` writes the bare column
    # names, here those of the golden corpus's table, two of which dbNSFP writes too.
    "TranscriptAnnotator": [
        "transcript_annotator_<column>",
        "Comment",
        "SIFT_pred",
        "SIFT_score",
    ],
    "UTRAnnotator": [
        "5UTR_annotation",
        "5UTR_consequence",
        "Existing_InFrame_oORFs",
        "Existing_OutOfFrame_oORFs",
        "Existing_uORFs",
    ],
}

# `HGVS` is in PLUGIN_FIELDS because this script doubles as the HGVS comparator,
# but it is not a plugin: `--plugins all` leaves it out, because a pseudo-plugin
# evaluated against plugin-run outputs can only report zero.
NON_PLUGIN_FIELD_SETS = {"HGVS"}
ALL_PLUGIN_NAMES = [k for k in PLUGIN_FIELDS if k not in NON_PLUGIN_FIELD_SETS]

# The name Ensembl VEP's command line loads a plugin under, where it is not the
# registry's.
PLUGIN_ALIASES = {"LoF": "LoFTEE"}

_PLACEHOLDER = re.compile(r"<(\w+)>")


def field_matcher(entries):
    """Return a predicate telling whether a field name matches a registry entry."""
    names = set()
    patterns = []
    for entry in entries:
        parts = _PLACEHOLDER.split(entry)
        if len(parts) == 1:
            names.add(entry)
            continue
        # The split alternates literal text with placeholder names.
        regex = ""
        for index, part in enumerate(parts):
            if index % 2 == 0:
                regex += re.escape(part)
            elif part == "N":
                regex += r"\d+"
            else:
                regex += r".+"
        patterns.append(regex)

    # The predicate runs once per Extra key of every row: the set's own method
    # serves an entry list without patterns, one alternation the others.
    if not patterns:
        return names.__contains__
    pattern = re.compile("|".join(patterns))

    def claims(field):
        return field in names or pattern.fullmatch(field) is not None

    return claims


# The lookup's matchers cover a plugin's CORE_NAMED_COLUMNS too. The comparison reads
# PLUGIN_FIELDS alone, where such a name would count a core key as the plugin's.
_FIELD_MATCHERS = {
    name: field_matcher([*PLUGIN_FIELDS[name], *CORE_NAMED_COLUMNS.get(name, ())])
    for name in ALL_PLUGIN_NAMES
}

_COLLECTIONS = (list, tuple, set, frozenset)


_ALIASES_OF = {
    name: [alias for alias, target in PLUGIN_ALIASES.items() if target == name]
    for name in PLUGIN_ALIASES.values()
}


@functools.cache
def _plugins_of_field(field):
    return tuple(name for name, claims in _FIELD_MATCHERS.items() if claims(field))


def plugin_of_field(field, plugins, flags=None):
    """Return the one plugin of a run that writes the Extra key `field`, or None.

    `plugins` names the plugins the run loaded, as the registry names them or as
    Ensembl VEP's command line does (`PLUGIN_ALIASES`); the answer is the registry's
    name. A name two plugins write, such as a dbNSFP column that a
    TranscriptAnnotator table under `prefix=0` also carries, belongs to the one of
    them the run loaded: with both loaded the answer is None, as it is for a key
    none of the run's plugins writes. A key of the HGVS field set belongs to no
    plugin.

    A column of CORE_NAMED_COLUMNS answers its plugin when the run loaded it.
    `flags` holds the words of the run's command line: a run that loaded the plugin
    and set a flag that turns the core key of that name on is refused with
    ValueError, since the name is then both the column's and the key's. Without
    `flags` no run can be refused.

    `plugins`, and `flags` when given, are a list, tuple or set: anything else
    raises TypeError, since a string would be searched by substring and an iterator
    only once.

    The answer comes from the registry as imported, so the `--fields` override of a
    command-line run does not reach it.
    """
    if not isinstance(plugins, _COLLECTIONS):
        raise TypeError("plugins is a list, tuple or set of the run's plugin names")
    if flags is not None and not isinstance(flags, _COLLECTIONS):
        raise TypeError("flags is a list, tuple or set of the run's command-line words")
    candidates = _plugins_of_field(field)
    if not candidates:
        return None
    writers = [
        name
        for name in candidates
        if name in plugins
        or any(alias in plugins for alias in _ALIASES_OF.get(name, ()))
    ]
    if len(writers) != 1:
        return None
    writer = writers[0]
    if flags is not None:
        for flag in CORE_NAMED_COLUMNS.get(writer, {}).get(field, ()):
            if flag in flags:
                raise ValueError(
                    f"the run loads {writer} and sets {flag}: {field} is both the "
                    "plugin's column and Ensembl VEP's own key"
                )
    return writer


def parse_extra_field(extra_str):
    """Parse VEP Extra column (key=value;key=value) into a dict."""
    if extra_str == "-" or not extra_str:
        return {}
    result = {}
    for pair in extra_str.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            result[key] = value
    return result


def is_vcf_path(filepath):
    lower = str(filepath).lower()
    return lower.endswith(".vcf") or lower.endswith(".vcf.gz")


def open_text_auto(path):
    """Open plain or gzip text file in UTF-8 with replacement mode."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def decode_csq_value(value):
    """Decode VCF CSQ percent-encoding (VEP-compatible)."""
    return unquote(value)


def parse_csq_header_fields(filepath):
    """Parse CSQ field names from a VCF header."""
    with open_text_auto(filepath) as f:
        for line in f:
            if not line.startswith("##INFO=<ID=CSQ"):
                continue
            marker = "Format: "
            idx = line.find(marker)
            if idx == -1:
                return []
            tail = line[idx + len(marker) :]
            # Header line is quoted: ... Format: A|B|C">
            fields = tail.split('"')[0]
            return fields.split("|")
    return []


def extract_info_value(info_str, key):
    """Extract key=value from VCF INFO; returns None if missing."""
    prefix = f"{key}="
    for field in info_str.split(";"):
        if field.startswith(prefix):
            return field[len(prefix) :]
    return None


def parse_vep_output(filepath, claims):
    """Parse a VEP text output file and extract the fields `claims` accepts.

    Returns a dict mapping (location, allele, feature) -> {field: value}.
    """
    records = {}
    with open_text_auto(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 14:
                continue

            location = parts[1]
            allele = parts[2]
            feature = parts[4]
            extra = parse_extra_field(parts[13])

            key = (location, allele, feature)
            # Interned, so the rows share one string per field name.
            plugin_data = {
                sys.intern(field): value
                for field, value in extra.items()
                if claims(field)
            }

            if plugin_data:
                records[key] = plugin_data

    return records


def parse_vep_vcf_output(filepath, claims):
    """Parse VCF output and extract the CSQ fields `claims` accepts.

    Returns a dict mapping (location, allele, feature) -> {field: value}.
    """
    records = {}
    csq_fields = parse_csq_header_fields(filepath)
    if not csq_fields:
        return records

    with open_text_auto(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = parts[1]
            info = parts[7]
            csq_value = extract_info_value(info, "CSQ")
            if not csq_value:
                continue

            entries = csq_value.split(",")
            for entry in entries:
                raw_values = entry.split("|")
                # Preserve empty trailing columns by padding to header width.
                if len(raw_values) < len(csq_fields):
                    raw_values += [""] * (len(csq_fields) - len(raw_values))
                elif len(raw_values) > len(csq_fields):
                    raw_values = raw_values[: len(csq_fields)]

                decoded_values = [decode_csq_value(v) for v in raw_values]
                csq_map = dict(zip(csq_fields, decoded_values))

                location = f"{chrom}:{pos}"
                allele = csq_map.get("Allele", "")
                feature = csq_map.get("Feature", "")
                key = (location, allele, feature)

                plugin_data = {
                    field: value
                    for field, value in csq_map.items()
                    if value not in ("", ".") and claims(field)
                }

                if plugin_data:
                    records[key] = plugin_data

    return records


def compare_records(perl_records, rust_records, tolerance):
    """Compare plugin fields between Perl and Rust outputs for a single plugin."""
    all_keys = set(perl_records.keys()) | set(rust_records.keys())

    stats = {
        "total_keys": len(all_keys),
        "match": 0,
        "mismatch": 0,
        "perl_only": 0,
        "rust_only": 0,
        "both_empty": 0,
    }
    mismatches = []

    for key in sorted(all_keys):
        perl_data = perl_records.get(key, {})
        rust_data = rust_records.get(key, {})

        if not perl_data and not rust_data:
            stats["both_empty"] += 1
            continue

        if perl_data and not rust_data:
            stats["perl_only"] += 1
            mismatches.append(
                {
                    "key": list(key),
                    "type": "perl_only",
                    "perl": dict(perl_data),
                    "rust": {},
                }
            )
            continue

        if rust_data and not perl_data:
            stats["rust_only"] += 1
            mismatches.append(
                {
                    "key": list(key),
                    "type": "rust_only",
                    "perl": {},
                    "rust": dict(rust_data),
                }
            )
            continue

        # Both have data -- compare field values
        # Use union of fields present in either side
        check_fields = set(perl_data.keys()) | set(rust_data.keys())

        field_match = True
        for field in check_fields:
            perl_val = perl_data.get(field, "")
            rust_val = rust_data.get(field, "")
            if perl_val == rust_val:
                continue
            # Numeric tolerance
            try:
                if abs(float(perl_val) - float(rust_val)) < tolerance:
                    continue
            except (ValueError, TypeError):
                pass
            field_match = False
            break

        if field_match:
            stats["match"] += 1
        else:
            stats["mismatch"] += 1
            mismatches.append(
                {
                    "key": list(key),
                    "type": "mismatch",
                    "perl": dict(perl_data),
                    "rust": dict(rust_data),
                }
            )

    annotated = (
        stats["match"] + stats["mismatch"] + stats["perl_only"] + stats["rust_only"]
    )
    # A zero-annotation run is NOT 100% concordant: it is unmeasured. Reporting
    # 100.0 here would let a completely dead plugin read as a pass: a contig-name
    # mismatch, an unparsed header, or a wrong position column all produce zero
    # records and no error. `None` forces callers to treat it as "no data".
    stats["concordance_pct"] = (
        (stats["match"] / annotated * 100) if annotated > 0 else None
    )
    stats["annotated_keys"] = annotated
    # Per-ENGINE annotated counts. The union count above cannot distinguish "both
    # engines annotated the same 12 records" from "only one engine annotated
    # anything", and a one-sided run is a failure even when every shared key
    # agrees.
    stats["perl_annotated"] = stats["match"] + stats["mismatch"] + stats["perl_only"]
    stats["rust_annotated"] = stats["match"] + stats["mismatch"] + stats["rust_only"]

    return stats, mismatches


def process_file_pair(perl_path, rust_path, plugins, tolerance, max_mismatches):
    """Process a single Perl/Rust output file pair for all requested plugins.

    `max_mismatches` caps the mismatches kept per plugin under `mismatches_sample`;
    `total_mismatches` and the statistics count every one.
    """
    results = {}

    # Each compared record needs one owner, so a name two requested plugins write
    # counts for the first of them in registry order, where plugin_of_field answers
    # None; a plugin the registry lacks has no entry.
    fields_of = {name: [] for name in plugins}
    matchers = {}
    for name, entries in PLUGIN_FIELDS.items():
        if name in fields_of:
            fields_of[name] = [
                entry
                for entry in entries
                if not any(claims(entry) for claims in matchers.values())
            ]
            matchers[name] = field_matcher(entries)

    compared_fields = set()
    for plugin_name in plugins:
        registered_fields = fields_of[plugin_name]
        claims = field_matcher(registered_fields)

        if is_vcf_path(perl_path):
            perl_records = parse_vep_vcf_output(perl_path, claims)
        else:
            perl_records = parse_vep_output(perl_path, claims)

        if is_vcf_path(rust_path):
            rust_records = parse_vep_vcf_output(rust_path, claims)
        else:
            rust_records = parse_vep_output(rust_path, claims)

        stats, mismatches = compare_records(perl_records, rust_records, tolerance)

        # The fields either output carries, or the entries looked for when neither
        # carries one.
        found = {
            field
            for records in (perl_records, rust_records)
            for data in records.values()
            for field in data
        }
        compared = [entry for entry in registered_fields if entry in found]
        compared += sorted(found.difference(registered_fields))
        compared_fields |= found

        results[plugin_name] = {
            "fields": compared or registered_fields,
            "shared_fields": [],
            "perl_annotated_records": len(perl_records),
            "rust_annotated_records": len(rust_records),
            "stats": stats,
            "mismatches_sample": mismatches[:max_mismatches],
            "total_mismatches": len(mismatches),
        }

    # A compared field that is an entry of two requested plugins is marked on both,
    # the one that owns its records and the one that does not.
    contested = {
        field
        for field in compared_fields
        if sum(writes(field) for writes in matchers.values()) > 1
    }
    for plugin_name, result in results.items():
        if plugin_name in matchers:
            result["shared_fields"] = sorted(filter(matchers[plugin_name], contested))

    return results


def canonical_stem(path_obj):
    """Normalize benchmark file names for perl/rust pair matching."""
    name = path_obj.name
    for suffix in [".vep.vcf.gz", ".vep.vcf", ".vcf.gz", ".vcf", ".txt", ".gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.startswith("perl_"):
        name = name[len("perl_") :]
    elif name.startswith("rust_"):
        name = name[len("rust_") :]
    return name


def candidate_output_files(root):
    root_path = Path(root)
    files = []
    for pattern in ("*.txt", "*.vep.vcf", "*.vep.vcf.gz", "*.vcf", "*.vcf.gz"):
        files.extend(root_path.glob(pattern))
    return files


def find_matching_pairs(perl_dir, rust_dir):
    """Find matching file pairs in perl and rust directories."""
    perl_files = {canonical_stem(f): f for f in candidate_output_files(perl_dir)}
    rust_files = {canonical_stem(f): f for f in candidate_output_files(rust_dir)}
    common = sorted(set(perl_files.keys()) & set(rust_files.keys()))
    return [(str(perl_files[name]), str(rust_files[name]), name) for name in common]


def print_report(all_results, plugins, min_annotated=1):
    """Print human-readable summary to stdout.

    `min_annotated` is the per-engine floor of annotated records below which a
    plugin FAILS. It defaults to 1 because the interesting failure is zero: a
    contig-name mismatch, an unparsed header, or a wrong position column all make
    a plugin annotate nothing while exiting 0.
    """
    any_failure = False

    for plugin_name in plugins:
        print(f"\n{'=' * 60}")
        print(f"  Plugin: {plugin_name}")
        print(f"{'=' * 60}")

        # Aggregate across files
        total_match = 0
        total_mismatch = 0
        total_perl_only = 0
        total_rust_only = 0
        total_annotated = 0
        total_perl_annotated = 0
        total_rust_annotated = 0

        for file_name, file_results in all_results.items():
            if plugin_name not in file_results:
                continue
            pr = file_results[plugin_name]
            s = pr["stats"]
            total_match += s["match"]
            total_mismatch += s["mismatch"]
            total_perl_only += s["perl_only"]
            total_rust_only += s["rust_only"]
            total_annotated += s["annotated_keys"]
            total_perl_annotated += s.get("perl_annotated", 0)
            total_rust_annotated += s.get("rust_annotated", 0)

            print(f"\n  File: {file_name}")
            print(
                f"    Fields: {', '.join(pr['fields'][:5])}{'...' if len(pr['fields']) > 5 else ''}"
            )
            if pr["shared_fields"]:
                print(
                    "    Also entries of another requested plugin: "
                    + ", ".join(pr["shared_fields"])
                )
            print(f"    Perl records: {pr['perl_annotated_records']}")
            print(f"    Rust records: {pr['rust_annotated_records']}")
            print(
                f"    Match: {s['match']}  Mismatch: {s['mismatch']}  "
                f"Perl-only: {s['perl_only']}  Rust-only: {s['rust_only']}"
            )
            if s["concordance_pct"] is None:
                print("    Concordance: n/a (no records annotated by either engine)")
            else:
                print(f"    Concordance: {s['concordance_pct']:.2f}%")

            if pr["mismatches_sample"]:
                print("    First mismatches:")
                for m in pr["mismatches_sample"][:3]:
                    loc, allele, feat = m["key"]
                    print(f"      {loc} {allele} {feat}: {m['type']}")
                    if m["perl"]:
                        print(f"        Perl: {m['perl']}")
                    if m["rust"]:
                        print(f"        Rust: {m['rust']}")

        if total_annotated > 0:
            agg = f"{total_match / total_annotated * 100:.2f}%"
        else:
            agg = "n/a"

        print(f"\n  AGGREGATE: {total_match}/{total_annotated} match ({agg})")
        print(
            f"  ANNOTATED: perl={total_perl_annotated}  rust={total_rust_annotated}"
            f"  (floor {min_annotated})"
        )

        # Reasons a plugin FAILS. Each is a shape that would otherwise pass silently:
        #   - either engine annotated fewer than `min_annotated` records: the run
        #     measured nothing, which is not agreement;
        #   - rust_only > 0: vep-rs over-annotated where Perl did not.
        reasons = []
        if total_perl_annotated < min_annotated:
            reasons.append(
                f"Perl annotated {total_perl_annotated} records (< {min_annotated})"
            )
        if total_rust_annotated < min_annotated:
            reasons.append(
                f"vep-rs annotated {total_rust_annotated} records (< {min_annotated})"
            )
        if total_mismatch > 0:
            reasons.append(f"{total_mismatch} value mismatches")
        if total_perl_only > 0:
            reasons.append(f"{total_perl_only} Perl-only records")
        if total_rust_only > 0:
            reasons.append(f"{total_rust_only} vep-rs-only records")

        if reasons:
            print(f"  STATUS: FAIL ({'; '.join(reasons)})")
            any_failure = True
        else:
            print("  STATUS: PASS")

    return any_failure


def main():
    parser = argparse.ArgumentParser(
        description="Compare VEP plugin outputs between Perl and Rust"
    )
    # Input modes: single file pair or directory
    parser.add_argument("--perl-output", help="Single Perl VEP output file")
    parser.add_argument("--rust-output", help="Single Rust VEP output file")
    parser.add_argument("--perl-dir", help="Directory of Perl VEP output files")
    parser.add_argument("--rust-dir", help="Directory of Rust VEP output files")

    # Plugin selection
    parser.add_argument(
        "--plugins",
        required=True,
        help="Comma-separated plugin names, or 'all' for all known plugins",
    )
    parser.add_argument(
        "--fields",
        default=None,
        help="Override: comma-separated field names (for single-plugin manual use)",
    )

    # Output
    parser.add_argument(
        "--json-output", default=None, help="Write JSON results to this file"
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=10,
        help=(
            "Max mismatches to include in JSON output per plugin per file "
            "(default: 10; 0 keeps none). The printed report shows at most three "
            "of them, so a value under 3 shortens it as well."
        ),
    )

    # Tolerance
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Numeric tolerance for floating-point comparison (default: 1e-6)",
    )
    parser.add_argument(
        "--min-annotated",
        type=int,
        default=1,
        help=(
            "Minimum annotated records EACH engine must produce, or the plugin "
            "FAILS (default: 1). A plugin that annotates nothing has not been "
            "measured -- a contig-name mismatch, an unparsed header, or a wrong "
            "position column all yield zero records without an error."
        ),
    )

    args = parser.parse_args()
    if args.max_mismatches < 0:
        # As a slice bound a negative value drops mismatches from the end.
        parser.error(
            "argument --max-mismatches: must be 0 or greater, "
            f"got {args.max_mismatches}"
        )

    # Resolve plugins
    if args.plugins.lower() == "all":
        plugins = ALL_PLUGIN_NAMES
    else:
        plugins = [p.strip() for p in args.plugins.split(",")]
        unknown = [p for p in plugins if p not in PLUGIN_FIELDS]
        if unknown:
            # A typo would otherwise reach the comparison without an entry, a field
            # list that accepts nothing, so the run would compare nothing under the
            # misspelt name.
            print(
                f"ERROR: unknown plugin(s): {', '.join(unknown)}. "
                f"Known: {', '.join(sorted(PLUGIN_FIELDS))}",
                file=sys.stderr,
            )
            return 2

    # Override fields for single-plugin manual use
    if args.fields and len(plugins) == 1:
        PLUGIN_FIELDS[plugins[0]] = [f.strip() for f in args.fields.split(",")]

    # Resolve file pairs
    if args.perl_dir and args.rust_dir:
        pairs = find_matching_pairs(args.perl_dir, args.rust_dir)
        if not pairs:
            print(
                f"ERROR: no matching file pairs found in {args.perl_dir} and {args.rust_dir}",
                file=sys.stderr,
            )
            sys.exit(1)
    elif args.perl_output and args.rust_output:
        name = Path(args.perl_output).stem
        pairs = [(args.perl_output, args.rust_output, name)]
    else:
        print(
            "ERROR: provide either --perl-dir/--rust-dir or --perl-output/--rust-output",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Plugins: {', '.join(plugins)}")
    print(f"File pairs: {len(pairs)}")
    print(f"Tolerance: {args.tolerance}")

    # Process all file pairs
    all_results = {}
    for perl_path, rust_path, file_name in pairs:
        all_results[file_name] = process_file_pair(
            perl_path, rust_path, plugins, args.tolerance, args.max_mismatches
        )

    # Print human-readable report
    any_failure = print_report(all_results, plugins, args.min_annotated)

    # Write JSON output
    if args.json_output:
        json_data = {
            "plugins": plugins,
            "tolerance": args.tolerance,
            "min_annotated": args.min_annotated,
            "file_count": len(pairs),
            "results": all_results,
        }

        # Add per-plugin aggregate stats
        aggregates = {}
        for plugin_name in plugins:
            agg = {
                "match": 0,
                "mismatch": 0,
                "perl_only": 0,
                "rust_only": 0,
                "annotated": 0,
                "perl_annotated": 0,
                "rust_annotated": 0,
            }
            for file_results in all_results.values():
                if plugin_name in file_results:
                    s = file_results[plugin_name]["stats"]
                    agg["match"] += s["match"]
                    agg["mismatch"] += s["mismatch"]
                    agg["perl_only"] += s["perl_only"]
                    agg["rust_only"] += s["rust_only"]
                    agg["annotated"] += s["annotated_keys"]
                    agg["perl_annotated"] += s.get("perl_annotated", 0)
                    agg["rust_annotated"] += s.get("rust_annotated", 0)
            # None, not 100.0, when nothing was annotated: unmeasured is not perfect.
            agg["concordance_pct"] = (
                agg["match"] / agg["annotated"] * 100 if agg["annotated"] > 0 else None
            )
            # Mirrors print_report's predicate: both engines must clear the floor,
            # and rust_only counts against the run.
            agg["pass"] = (
                agg["mismatch"] == 0
                and agg["perl_only"] == 0
                and agg["rust_only"] == 0
                and agg["perl_annotated"] >= args.min_annotated
                and agg["rust_annotated"] >= args.min_annotated
            )
            aggregates[plugin_name] = agg

        json_data["aggregates"] = aggregates

        Path(args.json_output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_output).write_text(
            json.dumps(json_data, indent=2), encoding="utf-8"
        )
        print(f"\nJSON results written to: {args.json_output}")

    # Exit code
    if any_failure:
        print("\nOVERALL: FAIL (one or more plugins have discordances)")
        sys.exit(1)
    else:
        print("\nOVERALL: PASS")
        sys.exit(0)


if __name__ == "__main__":
    # main() returns a status for its early-exit paths (e.g. an unknown plugin
    # name returns 2). Called bare, that status is discarded: a typo'd --plugins
    # prints ERROR yet exits 0, and a script gating on exit status reads a
    # misconfigured run as a pass. The success/failure paths call sys.exit()
    # directly, so this only ever propagates an early return.
    sys.exit(main() or 0)
