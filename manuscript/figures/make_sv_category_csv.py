"""Build manuscript/data/sv_concordance_by_category.csv.

One row per assembly for each of the 13 synthetic SV VCF files and the three
real-world SV sources: the vep-rs raw and adjusted F1 against the per-VCF VEP 115.2
reference output at cache version 115 (the GRCh37 set is generated against the
full-genome Ensembl GRCh37 cache) with the tuple counts each divides
(``perl_tuples``, ``vep_rs_tuples``, ``intersection`` and their ``adj_`` forms), and
the same for fastVEP (``fastvep_raw_f1``, ``fastvep_adj_f1`` with ``fastvep_tuples``,
``fastvep_intersection`` and the ``adj_`` forms; its adjusted VEP-side denominator is
``adj_perl_tuples_fastvep``, because the VEP-side mask under ``--scored-engine fastvep``
removes a different set). Backs supplementary Table S5; no figure reads it.

The values come from scripts/validation/compare_sv_concordance.py as run by
run_clone_measurement.sh over the VCFs in tests/sv_validation/, which writes one
report per assembly at s07/report/concordance_report.json (GRCh37) and
s08/report/concordance_report.json (GRCh38), one report set per engine. Table S4's
corpus figure is each report's ``overall`` and ``adjusted`` F1, taken once from the 16
files' pooled tuple counts; F1 is a ratio, so no mean of the rows here reproduces it,
and the reports' ``*_union`` counts, which deduplicate a tuple key seen in several
files, are not F1 denominators. Every F1 here is 2I / (P + R) of the counts beside it.

Set VEP_SV_REPORT_DIR (vep-rs) and VEP_SV_FASTVEP_REPORT_DIR (fastVEP) to run
directories with that layout to rebuild from live reports; a missing report raises
rather than emitting a partial CSV. Unset, BAKED_ROWS is emitted: the pinned sweep's
values, written from a live rebuild by ``--bake`` and never edited by hand. Both paths
emit the same rows in the same order.
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
from pathlib import Path

_sv_report_env = os.environ.get("VEP_SV_REPORT_DIR") or None
_fv_report_env = os.environ.get("VEP_SV_FASTVEP_REPORT_DIR") or None
USE_LIVE_REPORTS = _sv_report_env is not None
SV_REPORT_DIR = Path(_sv_report_env) if _sv_report_env else None
FASTVEP_REPORT_DIR = Path(_fv_report_env) if _fv_report_env else None

# Report suite id -> assembly; each report covers all 16 files of its assembly.
SV_SUITES = {"s07": "GRCh37", "s08": "GRCh38"}

# The reports name the gnomAD real-world source by its release (the GRCh37 set
# is gnomAD-SV v2.1, the GRCh38 set v4.1); the CSV and Table S5's row labels use
# one assembly-agnostic key so the table has a single gnomAD row.
REPORT_FILE_ALIASES = {
    "gnomad_sv_v2.1_chr21": "gnomad_sv_chr21",
    "gnomad_sv_v4.1_chr21": "gnomad_sv_chr21",
}

SYNTHETIC_FILES = [
    "01_snv",
    "02_mnp_complex",
    "03_small_indels",
    "04_large_explicit",
    "05_symbolic_del_ins",
    "06_symbolic_dup_inv",
    "07_cnv_repeat",
    "08_mobile_elements",
    "09_breakends",
    "10_special_alleles",
    "11_multi_allelic",
    "12_complex_imprecise",
    "13_vcf45_features",
]

FIELDS = [
    "file",
    "suite",
    "assembly",
    "n_variants",
    "raw_f1",
    "adj_f1",
    "perl_tuples",
    "vep_rs_tuples",
    "intersection",
    "adj_perl_tuples",
    "adj_vep_rs_tuples",
    "adj_intersection",
    "fastvep_raw_f1",
    "fastvep_adj_f1",
    "fastvep_tuples",
    "fastvep_intersection",
    "adj_perl_tuples_fastvep",
    "adj_fastvep_tuples",
    "adj_fastvep_intersection",
]

# Per-file values from the vep-rs binary and fastVEP release the sweeps
# VEP_RS_SV_SWEEP_ID and FASTVEP_SV_F1_SWEEP_ID in generate_figures.py name, written by
# `--bake` from the live reports. The gnomAD rows sit far below the rest on raw F1
# because the genome-wide reference lets VEP's cross-chromosome annotation fire on both
# assemblies, which raw F1 charges to vep-rs and the mask removes; the residual under
# adjusted F1 on those two files is transcript selection on breakends spanning more than
# 10 Mb, which vep-rs annotates against every overlapping transcript and VEP against the
# transcripts its batch loaded. 07_cnv_repeat's two columns differ by the <CNV:TR>
# representation pairs the symmetric mask removes from both sides.
# BAKED-BEGIN
BAKED_ROWS: list[dict[str, str]] = [
    {"file": "01_snv", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "1000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "7488", "vep_rs_tuples": "7488", "intersection": "7488", "adj_perl_tuples": "7488", "adj_vep_rs_tuples": "7488", "adj_intersection": "7488", "fastvep_raw_f1": "0.840025", "fastvep_adj_f1": "0.840025", "fastvep_tuples": "5464", "fastvep_intersection": "5440", "adj_perl_tuples_fastvep": "7488", "adj_fastvep_tuples": "5464", "adj_fastvep_intersection": "5440"},
    {"file": "01_snv", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "1000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "15755", "vep_rs_tuples": "15755", "intersection": "15755", "adj_perl_tuples": "15755", "adj_vep_rs_tuples": "15755", "adj_intersection": "15755", "fastvep_raw_f1": "1.000000", "fastvep_adj_f1": "1.000000", "fastvep_tuples": "15755", "fastvep_intersection": "15755", "adj_perl_tuples_fastvep": "15755", "adj_fastvep_tuples": "15755", "adj_fastvep_intersection": "15755"},
    {"file": "02_mnp_complex", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "2000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "19664", "vep_rs_tuples": "19664", "intersection": "19664", "adj_perl_tuples": "19664", "adj_vep_rs_tuples": "19664", "adj_intersection": "19664", "fastvep_raw_f1": "0.755493", "fastvep_adj_f1": "0.755493", "fastvep_tuples": "16016", "fastvep_intersection": "13478", "adj_perl_tuples_fastvep": "19664", "adj_fastvep_tuples": "16016", "adj_fastvep_intersection": "13478"},
    {"file": "02_mnp_complex", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "2000", "raw_f1": "0.999822", "adj_f1": "0.999822", "perl_tuples": "56171", "vep_rs_tuples": "56171", "intersection": "56161", "adj_perl_tuples": "56171", "adj_vep_rs_tuples": "56171", "adj_intersection": "56161", "fastvep_raw_f1": "0.843425", "fastvep_adj_f1": "0.843425", "fastvep_tuples": "56171", "fastvep_intersection": "47376", "adj_perl_tuples_fastvep": "56171", "adj_fastvep_tuples": "56171", "adj_fastvep_intersection": "47376"},
    {"file": "03_small_indels", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "2000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "19763", "vep_rs_tuples": "19763", "intersection": "19763", "adj_perl_tuples": "19763", "adj_vep_rs_tuples": "19763", "adj_intersection": "19763", "fastvep_raw_f1": "0.436462", "fastvep_adj_f1": "0.436462", "fastvep_tuples": "15594", "fastvep_intersection": "7716", "adj_perl_tuples_fastvep": "19763", "adj_fastvep_tuples": "15594", "adj_fastvep_intersection": "7716"},
    {"file": "03_small_indels", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "2000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "54896", "vep_rs_tuples": "54896", "intersection": "54896", "adj_perl_tuples": "54896", "adj_vep_rs_tuples": "54896", "adj_intersection": "54896", "fastvep_raw_f1": "0.492331", "fastvep_adj_f1": "0.492331", "fastvep_tuples": "54896", "fastvep_intersection": "27027", "adj_perl_tuples_fastvep": "54896", "adj_fastvep_tuples": "54896", "adj_fastvep_intersection": "27027"},
    {"file": "04_large_explicit", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "400", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "3985", "vep_rs_tuples": "3985", "intersection": "3985", "adj_perl_tuples": "3985", "adj_vep_rs_tuples": "3985", "adj_intersection": "3985", "fastvep_raw_f1": "0.293398", "fastvep_adj_f1": "0.293398", "fastvep_tuples": "3043", "fastvep_intersection": "1031", "adj_perl_tuples_fastvep": "3985", "adj_fastvep_tuples": "3043", "adj_fastvep_intersection": "1031"},
    {"file": "04_large_explicit", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "400", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "10790", "vep_rs_tuples": "10790", "intersection": "10790", "adj_perl_tuples": "10790", "adj_vep_rs_tuples": "10790", "adj_intersection": "10790", "fastvep_raw_f1": "0.356997", "fastvep_adj_f1": "0.356997", "fastvep_tuples": "10790", "fastvep_intersection": "3852", "adj_perl_tuples_fastvep": "10790", "adj_fastvep_tuples": "10790", "adj_fastvep_intersection": "3852"},
    {"file": "05_symbolic_del_ins", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "2000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "73007", "vep_rs_tuples": "73007", "intersection": "73007", "adj_perl_tuples": "73007", "adj_vep_rs_tuples": "73007", "adj_intersection": "73007", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "32046", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "73007", "adj_fastvep_tuples": "32046", "adj_fastvep_intersection": "0"},
    {"file": "05_symbolic_del_ins", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "2000", "raw_f1": "0.999991", "adj_f1": "0.999991", "perl_tuples": "117575", "vep_rs_tuples": "117575", "intersection": "117574", "adj_perl_tuples": "117575", "adj_vep_rs_tuples": "117575", "adj_intersection": "117574", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "96548", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "117575", "adj_fastvep_tuples": "96548", "adj_fastvep_intersection": "0"},
    {"file": "06_symbolic_dup_inv", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "2250", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "67036", "vep_rs_tuples": "67036", "intersection": "67036", "adj_perl_tuples": "67036", "adj_vep_rs_tuples": "67036", "adj_intersection": "67036", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "38002", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "67036", "adj_fastvep_tuples": "38002", "adj_fastvep_intersection": "0"},
    {"file": "06_symbolic_dup_inv", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "2250", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "106265", "vep_rs_tuples": "106265", "intersection": "106265", "adj_perl_tuples": "106265", "adj_vep_rs_tuples": "106265", "adj_intersection": "106265", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "106014", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "106265", "adj_fastvep_tuples": "106014", "adj_fastvep_intersection": "0"},
    {"file": "07_cnv_repeat", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "900", "raw_f1": "0.936322", "adj_f1": "1.000000", "perl_tuples": "21483", "vep_rs_tuples": "21483", "intersection": "20115", "adj_perl_tuples": "20115", "adj_vep_rs_tuples": "20115", "adj_intersection": "20115", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "12613", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "21396", "adj_fastvep_tuples": "12526", "adj_fastvep_intersection": "0"},
    {"file": "07_cnv_repeat", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "900", "raw_f1": "0.923333", "adj_f1": "0.999969", "perl_tuples": "34500", "vep_rs_tuples": "34500", "intersection": "31855", "adj_perl_tuples": "31856", "adj_vep_rs_tuples": "31856", "adj_intersection": "31855", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "34500", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "34335", "adj_fastvep_tuples": "34335", "adj_fastvep_intersection": "0"},
    {"file": "08_mobile_elements", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "800", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "7607", "vep_rs_tuples": "7607", "intersection": "7607", "adj_perl_tuples": "7607", "adj_vep_rs_tuples": "7607", "adj_intersection": "7607", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "6240", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "7607", "adj_fastvep_tuples": "6240", "adj_fastvep_intersection": "0"},
    {"file": "08_mobile_elements", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "800", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "22042", "vep_rs_tuples": "22042", "intersection": "22042", "adj_perl_tuples": "22042", "adj_vep_rs_tuples": "22042", "adj_intersection": "22042", "fastvep_raw_f1": "0.001497", "fastvep_adj_f1": "0.001497", "fastvep_tuples": "22042", "fastvep_intersection": "33", "adj_perl_tuples_fastvep": "22042", "adj_fastvep_tuples": "22042", "adj_fastvep_intersection": "33"},
    {"file": "09_breakends", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "1400", "raw_f1": "0.997981", "adj_f1": "0.998969", "perl_tuples": "24266", "vep_rs_tuples": "24268", "intersection": "24218", "adj_perl_tuples": "24218", "adj_vep_rs_tuples": "24268", "adj_intersection": "24218", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "10279", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "24218", "adj_fastvep_tuples": "10279", "adj_fastvep_intersection": "0"},
    {"file": "09_breakends", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "1400", "raw_f1": "0.998294", "adj_f1": "0.998992", "perl_tuples": "71483", "vep_rs_tuples": "71511", "intersection": "71375", "adj_perl_tuples": "71383", "adj_vep_rs_tuples": "71511", "adj_intersection": "71375", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "37143", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "71383", "adj_fastvep_tuples": "37143", "adj_fastvep_intersection": "0"},
    {"file": "10_special_alleles", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "1000", "raw_f1": "0.955571", "adj_f1": "0.955571", "perl_tuples": "2865", "vep_rs_tuples": "3077", "intersection": "2839", "adj_perl_tuples": "2865", "adj_vep_rs_tuples": "3077", "adj_intersection": "2839", "fastvep_raw_f1": "0.142330", "fastvep_adj_f1": "0.142330", "fastvep_tuples": "7969", "fastvep_intersection": "771", "adj_perl_tuples_fastvep": "2865", "adj_fastvep_tuples": "7969", "adj_fastvep_intersection": "771"},
    {"file": "10_special_alleles", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "1000", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "7998", "vep_rs_tuples": "7998", "intersection": "7998", "adj_perl_tuples": "7998", "adj_vep_rs_tuples": "7998", "adj_intersection": "7998", "fastvep_raw_f1": "0.143547", "fastvep_adj_f1": "0.143547", "fastvep_tuples": "27628", "fastvep_intersection": "2557", "adj_perl_tuples_fastvep": "7998", "adj_fastvep_tuples": "27628", "adj_fastvep_intersection": "2557"},
    {"file": "11_multi_allelic", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "700", "raw_f1": "0.999533", "adj_f1": "0.999533", "perl_tuples": "8565", "vep_rs_tuples": "8565", "intersection": "8561", "adj_perl_tuples": "8565", "adj_vep_rs_tuples": "8565", "adj_intersection": "8561", "fastvep_raw_f1": "0.153738", "fastvep_adj_f1": "0.153738", "fastvep_tuples": "12822", "fastvep_intersection": "1644", "adj_perl_tuples_fastvep": "8565", "adj_fastvep_tuples": "12822", "adj_fastvep_intersection": "1644"},
    {"file": "11_multi_allelic", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "700", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "21619", "vep_rs_tuples": "21619", "intersection": "21619", "adj_perl_tuples": "21619", "adj_vep_rs_tuples": "21619", "adj_intersection": "21619", "fastvep_raw_f1": "0.170714", "fastvep_adj_f1": "0.170714", "fastvep_tuples": "43238", "fastvep_intersection": "5536", "adj_perl_tuples_fastvep": "21619", "adj_fastvep_tuples": "43238", "adj_fastvep_intersection": "5536"},
    {"file": "12_complex_imprecise", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "300", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "6548", "vep_rs_tuples": "6548", "intersection": "6548", "adj_perl_tuples": "6548", "adj_vep_rs_tuples": "6548", "adj_intersection": "6548", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "4157", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "6548", "adj_fastvep_tuples": "4157", "adj_fastvep_intersection": "0"},
    {"file": "12_complex_imprecise", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "300", "raw_f1": "0.999008", "adj_f1": "0.999008", "perl_tuples": "13099", "vep_rs_tuples": "13099", "intersection": "13086", "adj_perl_tuples": "13099", "adj_vep_rs_tuples": "13099", "adj_intersection": "13086", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "13099", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "13099", "adj_fastvep_tuples": "13099", "adj_fastvep_intersection": "0"},
    {"file": "13_vcf45_features", "suite": "synthetic_grch37", "assembly": "GRCh37", "n_variants": "650", "raw_f1": "0.999014", "adj_f1": "0.999014", "perl_tuples": "8615", "vep_rs_tuples": "8632", "intersection": "8615", "adj_perl_tuples": "8615", "adj_vep_rs_tuples": "8632", "adj_intersection": "8615", "fastvep_raw_f1": "0.156590", "fastvep_adj_f1": "0.156590", "fastvep_tuples": "10671", "fastvep_intersection": "1510", "adj_perl_tuples_fastvep": "8615", "adj_fastvep_tuples": "10671", "adj_fastvep_intersection": "1510"},
    {"file": "13_vcf45_features", "suite": "synthetic_grch38", "assembly": "GRCh38", "n_variants": "650", "raw_f1": "0.999843", "adj_f1": "0.999843", "perl_tuples": "22274", "vep_rs_tuples": "22281", "intersection": "22274", "adj_perl_tuples": "22274", "adj_vep_rs_tuples": "22281", "adj_intersection": "22274", "fastvep_raw_f1": "0.153373", "fastvep_adj_f1": "0.153373", "fastvep_tuples": "39445", "fastvep_intersection": "4733", "adj_perl_tuples_fastvep": "22274", "adj_fastvep_tuples": "39445", "adj_fastvep_intersection": "4733"},
    {"file": "1kg_sv_chr21", "suite": "real_world_grch37", "assembly": "GRCh37", "n_variants": "877", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "3057", "vep_rs_tuples": "3057", "intersection": "3057", "adj_perl_tuples": "3057", "adj_vep_rs_tuples": "3057", "adj_intersection": "3057", "fastvep_raw_f1": "0.078308", "fastvep_adj_f1": "0.078308", "fastvep_tuples": "2332", "fastvep_intersection": "211", "adj_perl_tuples_fastvep": "3057", "adj_fastvep_tuples": "2332", "adj_fastvep_intersection": "211"},
    {"file": "1kg_sv_chr21", "suite": "real_world_grch38", "assembly": "GRCh38", "n_variants": "2386", "raw_f1": "0.999762", "adj_f1": "0.999762", "perl_tuples": "25183", "vep_rs_tuples": "25183", "intersection": "25177", "adj_perl_tuples": "25183", "adj_vep_rs_tuples": "25183", "adj_intersection": "25177", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "30154", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "25183", "adj_fastvep_tuples": "30154", "adj_fastvep_intersection": "0"},
    {"file": "clinvar_sv_chr21", "suite": "real_world_grch37", "assembly": "GRCh37", "n_variants": "2794", "raw_f1": "0.999952", "adj_f1": "0.999952", "perl_tuples": "20913", "vep_rs_tuples": "20913", "intersection": "20912", "adj_perl_tuples": "20913", "adj_vep_rs_tuples": "20913", "adj_intersection": "20912", "fastvep_raw_f1": "0.507189", "fastvep_adj_f1": "0.507189", "fastvep_tuples": "16016", "fastvep_intersection": "9365", "adj_perl_tuples_fastvep": "20913", "adj_fastvep_tuples": "16016", "adj_fastvep_intersection": "9365"},
    {"file": "clinvar_sv_chr21", "suite": "real_world_grch38", "assembly": "GRCh38", "n_variants": "2796", "raw_f1": "1.000000", "adj_f1": "1.000000", "perl_tuples": "56361", "vep_rs_tuples": "56361", "intersection": "56361", "adj_perl_tuples": "56361", "adj_vep_rs_tuples": "56361", "adj_intersection": "56361", "fastvep_raw_f1": "0.586112", "fastvep_adj_f1": "0.586112", "fastvep_tuples": "56928", "fastvep_intersection": "33200", "adj_perl_tuples_fastvep": "56361", "adj_fastvep_tuples": "56928", "adj_fastvep_intersection": "33200"},
    {"file": "gnomad_sv_chr21", "suite": "real_world_grch37", "assembly": "GRCh37", "n_variants": "4921", "raw_f1": "0.832168", "adj_f1": "0.990660", "perl_tuples": "35002", "vep_rs_tuples": "44530", "intersection": "33092", "adj_perl_tuples": "33212", "adj_vep_rs_tuples": "33596", "adj_intersection": "33092", "fastvep_raw_f1": "0.000000", "fastvep_adj_f1": "0.000000", "fastvep_tuples": "12798", "fastvep_intersection": "0", "adj_perl_tuples_fastvep": "33212", "adj_fastvep_tuples": "12798", "adj_fastvep_intersection": "0"},
    {"file": "gnomad_sv_chr21", "suite": "real_world_grch38", "assembly": "GRCh38", "n_variants": "32630", "raw_f1": "0.813283", "adj_f1": "0.995552", "perl_tuples": "481026", "vep_rs_tuples": "644840", "intersection": "457824", "adj_perl_tuples": "458645", "adj_vep_rs_tuples": "461094", "adj_intersection": "457824", "fastvep_raw_f1": "0.000005", "fastvep_adj_f1": "0.000005", "fastvep_tuples": "338058", "fastvep_intersection": "2", "adj_perl_tuples_fastvep": "458645", "adj_fastvep_tuples": "338058", "adj_fastvep_intersection": "2"},
]
# BAKED-END


def _load_live_report(base: Path, suite: str, engine: str) -> dict:
    """Read one assembly's per-VCF SV concordance report for one engine.

    Raises if the report is absent: a live rebuild that quietly skipped an
    assembly would emit a CSV missing that assembly's rows.
    """
    path = base / suite / "report" / "concordance_report.json"
    if not path.is_file():
        raise SystemExit(
            f"the {engine} report directory is set but {path} is missing; point it at a run "
            "directory containing s07/report/concordance_report.json and "
            "s08/report/concordance_report.json"
        )
    with path.open() as fh:
        return json.load(fh)


def _f6(x: float) -> str:
    return f"{x:.6f}"


def _live_rows() -> list[dict[str, str]]:
    """Build CSV rows from the live per-assembly SV concordance reports of both engines.

    Rows are emitted in the same order as the baked path: each synthetic file's GRCh37
    then GRCh38 row, in file order, then the real-world sources by source and assembly,
    so a live rebuild and a baked rebuild of the same values are byte-identical.
    """
    if FASTVEP_REPORT_DIR is None:
        raise SystemExit(
            "VEP_SV_REPORT_DIR is set but VEP_SV_FASTVEP_REPORT_DIR is not; the CSV carries both "
            "engines' per-file values, so a rebuild needs both report sets"
        )
    rows: list[dict[str, str]] = []
    for suite, assembly in SV_SUITES.items():
        vr = _load_live_report(SV_REPORT_DIR, suite, "vep-rs")
        fr = _load_live_report(FASTVEP_REPORT_DIR, suite, "fastVEP")
        if set(vr["per_file"]) != set(fr["per_file"]):
            raise SystemExit(
                f"{suite}: the vep-rs and fastVEP reports cover different files: "
                f"{sorted(set(vr['per_file']) ^ set(fr['per_file']))}"
            )
        for report_name, raw in vr["per_file"].items():
            adj = vr["per_file_adjusted"][report_name]
            fraw = fr["per_file"][report_name]
            fadj = fr["per_file_adjusted"][report_name]
            if fraw["perl_tuples"] != raw["perl_tuples"] or fraw["total_input"] != raw["total_input"]:
                raise SystemExit(
                    f"{suite}/{report_name}: the two engines were scored against different VEP "
                    f"reference rows ({raw['perl_tuples']} vs {fraw['perl_tuples']} tuples)"
                )
            name = REPORT_FILE_ALIASES.get(report_name, report_name)
            kind = "synthetic" if name in SYNTHETIC_FILES else "real_world"
            rows.append(
                {
                    "file": name,
                    "suite": f"{kind}_{assembly.lower()}",
                    "assembly": assembly,
                    "n_variants": str(raw["total_input"]),
                    "raw_f1": _f6(raw["f1"]),
                    "adj_f1": _f6(adj["f1"]),
                    "perl_tuples": str(raw["perl_tuples"]),
                    "vep_rs_tuples": str(raw["rust_tuples"]),
                    "intersection": str(raw["intersection"]),
                    "adj_perl_tuples": str(adj["perl_tuples"]),
                    "adj_vep_rs_tuples": str(adj["rust_tuples"]),
                    "adj_intersection": str(adj["intersection"]),
                    "fastvep_raw_f1": _f6(fraw["f1"]),
                    "fastvep_adj_f1": _f6(fadj["f1"]),
                    "fastvep_tuples": str(fraw["rust_tuples"]),
                    "fastvep_intersection": str(fraw["intersection"]),
                    "adj_perl_tuples_fastvep": str(fadj["perl_tuples"]),
                    "adj_fastvep_tuples": str(fadj["rust_tuples"]),
                    "adj_fastvep_intersection": str(fadj["intersection"]),
                }
            )
    real_order = sorted({r["file"] for r in rows if r["file"] not in SYNTHETIC_FILES})
    order = {name: i for i, name in enumerate(SYNTHETIC_FILES + real_order)}
    rows.sort(key=lambda r: (order[r["file"]], r["assembly"]))
    return rows


def _check_arithmetic(rows: list[dict[str, str]]) -> None:
    """Every F1 must be 2I / (P + R) of the counts on its row, to six decimals."""
    for r in rows:
        for f1, p, rr, i in (
            ("raw_f1", "perl_tuples", "vep_rs_tuples", "intersection"),
            ("adj_f1", "adj_perl_tuples", "adj_vep_rs_tuples", "adj_intersection"),
            ("fastvep_raw_f1", "perl_tuples", "fastvep_tuples", "fastvep_intersection"),
            ("fastvep_adj_f1", "adj_perl_tuples_fastvep", "adj_fastvep_tuples", "adj_fastvep_intersection"),
        ):
            denom = int(r[p]) + int(r[rr])
            want = _f6(2 * int(r[i]) / denom) if denom else _f6(0.0)
            if want != r[f1]:
                raise SystemExit(
                    f"{r['file']} {r['assembly']}: {f1} {r[f1]} is not 2*{r[i]}/({r[p]}+{r[rr]}) = {want}"
                )


def _bake(rows: list[dict[str, str]]) -> None:
    """Rewrite BAKED_ROWS in this file from the live rows."""
    src_path = Path(__file__).resolve()
    src = src_path.read_text(encoding="utf-8")
    body = "BAKED_ROWS: list[dict[str, str]] = [\n"
    for r in rows:
        body += "    {" + ", ".join(f'"{k}": "{r[k]}"' for k in FIELDS) + "},\n"
    body += "]\n"
    new = re.sub(
        r"# BAKED-BEGIN\n.*?# BAKED-END\n",
        "# BAKED-BEGIN\n" + body + "# BAKED-END\n",
        src,
        flags=re.DOTALL,
    )
    if new == src:
        print("baked rows unchanged")
    else:
        src_path.write_text(new, encoding="utf-8")
        print(f"baked {len(rows)} rows into {src_path.name}")


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "sv_concordance_by_category.csv"

    if USE_LIVE_REPORTS:
        rows = _live_rows()
        if "--bake" in sys.argv[1:]:
            _bake(rows)
    else:
        if "--bake" in sys.argv[1:]:
            raise SystemExit("--bake needs VEP_SV_REPORT_DIR and VEP_SV_FASTVEP_REPORT_DIR set to live reports")
        if not BAKED_ROWS:
            raise SystemExit("BAKED_ROWS is empty; run once with live reports and --bake")
        rows = [dict(r) for r in BAKED_ROWS]
    if len(rows) != 32:
        raise SystemExit(f"{len(rows)} rows, expected 32 (16 files x 2 assemblies)")
    _check_arithmetic(rows)

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()
