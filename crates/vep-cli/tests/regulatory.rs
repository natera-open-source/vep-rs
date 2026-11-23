// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use serde_json::Value;
use std::fs;
use std::path::Path;
mod common;

fn cache(dir: &Path) -> std::path::PathBuf {
    common::regulatory_cache(dir)
}

fn annotate(dir: &Path, cache: &Path, format: &str, regulatory: bool) -> String {
    annotate_with_args(dir, cache, format, regulatory, &[])
}

fn annotate_with_args(
    dir: &Path,
    cache: &Path,
    format: &str,
    regulatory: bool,
    extra_args: &[&str],
) -> String {
    let output = dir.join(format!("out-{format}-{regulatory}"));
    let mut command = common::isolated_vep(dir);
    command
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--force_overwrite",
            "--format",
            "vcf",
            "--output_format",
            if format == "default" { "vep" } else { format },
        ])
        .args(if extra_args.contains(&"--fasta") {
            &[][..]
        } else {
            &["--no_fasta"][..]
        })
        .arg("--input_file")
        .arg(dir.join("in.vcf"))
        .arg("--output_file")
        .arg(&output);
    if regulatory {
        command.arg("--regulatory");
    }
    command.args(extra_args);
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!String::from_utf8_lossy(&result.stderr).contains("no regulatory data"));
    fs::read_to_string(output).unwrap()
}

#[test]
fn regulatory_rows_and_extras_reach_every_text_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(dir.path().join("in.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n21\t115\treverse\tA\tG\t.\t.\t.\n").unwrap();
    let json = annotate(dir.path(), &cache, "json", true);
    let records: Vec<Value> = json
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    let first = &records[0];
    assert_eq!(
        first["regulatory_feature_consequences"][0]["regulatory_feature_id"],
        "ENSR_FIXTURE"
    );
    assert_eq!(
        first["regulatory_feature_consequences"][0]["biotype"],
        "enhancer"
    );
    assert!(first["regulatory_feature_consequences"][0]
        .get("strand")
        .is_none());
    let motif = &first["motif_feature_consequences"][0];
    assert_eq!(motif["motif_feature_id"], "ENSM_FORWARD");
    assert_eq!(motif["motif_name"], "ENSPFM_FORWARD");
    assert_eq!(motif["motif_pos"], 2);
    assert_eq!(
        motif["transcription_factors"],
        serde_json::json!(["TF_A", "TF_B"])
    );
    // The fixture matrices hold one count per position, so a base the motif's
    // sequence holds scores the best site and any other base a quarter less;
    // no position of a matrix that flat is informative.
    assert_eq!(motif["motif_score_change"], serde_json::json!(-0.25));
    assert_eq!(motif["high_inf_pos"], "N");
    assert_eq!(first["most_severe_consequence"], "TF_binding_site_variant");
    assert_eq!(records[1]["motif_feature_consequences"][0]["motif_pos"], 4);
    assert_eq!(
        records[1]["motif_feature_consequences"][0]["transcription_factors"],
        serde_json::json!([])
    );
    assert_eq!(
        records[1]["motif_feature_consequences"][0]["motif_score_change"],
        serde_json::json!(0)
    );
    assert_eq!(
        records[1]["motif_feature_consequences"][0]["high_inf_pos"],
        "N"
    );
    for format in ["default", "tab", "vcf"] {
        let text = annotate(dir.path(), &cache, format, true);
        let rows = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            rows.contains("regulatory_region_variant"),
            "{format}: {rows}"
        );
        assert!(rows.contains("TF_binding_site_variant"), "{format}: {rows}");
        assert!(rows.contains("enhancer"), "{format}: {rows}");
        assert!(rows.contains("ENSPFM_FORWARD"), "{format}: {rows}");
        assert!(rows.contains("TF_A"), "{format}: {rows}");
    }
    let disabled = annotate(dir.path(), &cache, "json", false);
    for line in disabled.lines() {
        let record: Value = serde_json::from_str(line).unwrap();
        assert!(record.get("regulatory_feature_consequences").is_none());
        assert!(record.get("motif_feature_consequences").is_none());
        assert_eq!(record["most_severe_consequence"], "intergenic_variant");
    }
    // Metadata is authoritative even when a regulatory tree exists.
    fs::write(cache.join("info.json"), r#"{"assembly":"GRCh38"}"#).unwrap();
    fs::write(cache.join("regulatory/21/1-1000000.json"), "[").unwrap();
    let mut command = common::isolated_vep(dir.path());
    // An output of STDOUT writes its warnings file beside the working directory.
    let result = command
        .current_dir(dir.path())
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--regulatory",
            "--no_fasta",
            "--no_stats",
            "--json",
            "--input_file",
        ])
        .arg(dir.path().join("in.vcf"))
        .args(["--output_file", "STDOUT"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("this cache holds no regulatory data"));
    assert!(!String::from_utf8_lossy(&result.stdout).contains("motif_feature_consequences"));
}

/// Two records of the gnomAD v4.1 chromosome 21 dataset against the slice of
/// the converted release 116 cache (`fixtures/regulatory/chr21-slice.json`):
/// `rs868224642` (`21:46699255 G>C`) inside the enhancer `ENSR21_C6M55` and
/// its four `NRF1` motifs on both strands, and `21:5090724 A>T` inside two
/// motifs of different matrices and no regulatory feature. Every row is held
/// whole, in the default, tab, VCF and JSON formats, to the row Ensembl VEP
/// 116.3 writes for the same record on the same cache, `HIGH_INF_POS` and
/// `MOTIF_SCORE_CHANGE` included (`Y`, `-0.094`, `Y`, `-0.094`, `N`, `-0.070`,
/// `N`, `-0.082` on the NRF1 rows; `Y`, `-0.089` and `N`, `0.039` on the two
/// others).
#[test]
fn regulatory_rows_are_byte_exact_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache_of_real_features(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t5090724\t.\tA\tT\t.\t.\t.\n\
         21\t46699255\trs868224642\tG\tC\t.\t.\t.\n",
    )
    .unwrap();
    let rows = |format: &str| -> Vec<String> {
        annotate(dir.path(), &cache, format, true)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_owned)
            .collect()
    };
    let nrf1 = |id: &str, strand: &str, pos: &str, informative: &str, delta: &str| {
        format!(
            "rs868224642\t21:46699255\tC\t-\t{id}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM0427;MOTIF_POS={pos};HIGH_INF_POS={informative};MOTIF_SCORE_CHANGE={delta};TRANSCRIPTION_FACTORS=NRF1"
        )
    };
    assert_eq!(
        rows("default"),
        [
            "21_5090724_A/T\t21:5090724\tT\t-\tENSM00000574814\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM0290;MOTIF_POS=7;HIGH_INF_POS=Y;MOTIF_SCORE_CHANGE=-0.089;TRANSCRIPTION_FACTORS=HNF4A,NR2F1,RXRA".to_string(),
            "21_5090724_A/T\t21:5090724\tT\t-\tENSM00000574815\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;MOTIF_NAME=ENSPFM0484;MOTIF_POS=15;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=0.039;TRANSCRIPTION_FACTORS=RUNX3".to_string(),
            "21_5090724_A/T\t21:5090724\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
            "rs868224642\t21:46699255\tC\t-\tENSR21_C6M55\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
            nrf1("ENSM00000588223", "-1", "5", "Y", "-0.094"),
            nrf1("ENSM00000588224", "1", "8", "Y", "-0.094"),
            nrf1("ENSM00000588225", "-1", "11", "N", "-0.070"),
            nrf1("ENSM00000588226", "1", "2", "N", "-0.082"),
            "rs868224642\t21:46699255\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        ]
    );
    let tab_motif = |uploaded: &str,
                     location: &str,
                     allele: &str,
                     id: &str,
                     strand: &str,
                     matrix: &str,
                     pos: &str,
                     informative: &str,
                     delta: &str,
                     factors: &str| {
        format!(
            "{uploaded}\t{location}\t{allele}\t-\t{id}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t{strand}\t-\t-\t{matrix}\t{pos}\t{informative}\t{delta}\t{factors}"
        )
    };
    assert_eq!(
        rows("tab"),
        [
            tab_motif("21_5090724_A/T", "21:5090724", "T", "ENSM00000574814", "1", "ENSPFM0290", "7", "Y", "-0.089", "HNF4A,NR2F1,RXRA"),
            tab_motif("21_5090724_A/T", "21:5090724", "T", "ENSM00000574815", "-1", "ENSPFM0484", "15", "N", "0.039", "RUNX3"),
            "21_5090724_A/T\t21:5090724\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-".to_string(),
            "rs868224642\t21:46699255\tC\t-\tENSR21_C6M55\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\tenhancer\t-\t-\t-\t-\t-".to_string(),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588223", "-1", "ENSPFM0427", "5", "Y", "-0.094", "NRF1"),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588224", "1", "ENSPFM0427", "8", "Y", "-0.094", "NRF1"),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588225", "-1", "ENSPFM0427", "11", "N", "-0.070", "NRF1"),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588226", "1", "ENSPFM0427", "2", "N", "-0.082", "NRF1"),
            "rs868224642\t21:46699255\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-".to_string(),
        ]
    );
    let vcf: Vec<String> = rows("vcf")
        .iter()
        .map(|line| line.split('\t').nth(7).unwrap().to_string())
        .collect();
    let csq_motif = |allele: &str,
                     id: &str,
                     strand: &str,
                     matrix: &str,
                     pos: &str,
                     informative: &str,
                     delta: &str,
                     factors: &str| {
        format!("{allele}|TF_binding_site_variant|MODIFIER|||MotifFeature|{id}|||||||||||||{strand}||||{matrix}|{pos}|{informative}|{delta}|{factors}")
    };
    assert_eq!(
        vcf,
        [
            format!(
                "CSQ={},{},T|intergenic_variant|MODIFIER|||||||||||||||||||||||||",
                csq_motif("T", "ENSM00000574814", "1", "ENSPFM0290", "7", "Y", "-0.089", "HNF4A&NR2F1&RXRA"),
                csq_motif("T", "ENSM00000574815", "-1", "ENSPFM0484", "15", "N", "0.039", "RUNX3")
            ),
            format!(
                "CSQ=C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR21_C6M55|enhancer||||||||||||||||||||,{},{},{},{},C|intergenic_variant|MODIFIER|||||||||||||||||||||||||",
                csq_motif("C", "ENSM00000588223", "-1", "ENSPFM0427", "5", "Y", "-0.094", "NRF1"),
                csq_motif("C", "ENSM00000588224", "1", "ENSPFM0427", "8", "Y", "-0.094", "NRF1"),
                csq_motif("C", "ENSM00000588225", "-1", "ENSPFM0427", "11", "N", "-0.070", "NRF1"),
                csq_motif("C", "ENSM00000588226", "1", "ENSPFM0427", "2", "N", "-0.082", "NRF1")
            ),
        ]
    );
    let json: Vec<Value> = rows("json")
        .iter()
        .map(|line| {
            let mut record: Value = serde_json::from_str(line).unwrap();
            // The input line and the version key are vep-rs's own, outside the comparison.
            record.as_object_mut().unwrap().remove("input");
            record.as_object_mut().unwrap().remove("vep_rs_version");
            record
        })
        .collect();
    let motif = |id: &str,
                 matrix: &str,
                 pos: u32,
                 informative: &str,
                 delta: f64,
                 strand: i32,
                 factors: &[&str],
                 allele: &str| {
        serde_json::json!({
            "consequence_terms": ["TF_binding_site_variant"], "high_inf_pos": informative,
            "impact": "MODIFIER", "motif_feature_id": id, "motif_name": matrix, "motif_pos": pos,
            "motif_score_change": delta, "strand": strand, "transcription_factors": factors,
            "variant_allele": allele
        })
    };
    assert_eq!(
        json,
        [
            serde_json::json!({
                "allele_string": "A/T", "assembly_name": "GRCh38", "end": 5090724, "id": ".",
                "intergenic_consequences": [{"consequence_terms": ["intergenic_variant"], "impact": "MODIFIER", "variant_allele": "T"}],
                "most_severe_consequence": "TF_binding_site_variant",
                "motif_feature_consequences": [
                    motif("ENSM00000574814", "ENSPFM0290", 7, "Y", -0.089, 1, &["HNF4A", "NR2F1", "RXRA"], "T"),
                    motif("ENSM00000574815", "ENSPFM0484", 15, "N", 0.039, -1, &["RUNX3"], "T"),
                ],
                "seq_region_name": "21", "start": 5090724, "strand": 1
            }),
            serde_json::json!({
                "allele_string": "G/C", "assembly_name": "GRCh38", "end": 46699255, "id": "rs868224642",
                "intergenic_consequences": [{"consequence_terms": ["intergenic_variant"], "impact": "MODIFIER", "variant_allele": "C"}],
                "most_severe_consequence": "TF_binding_site_variant",
                "motif_feature_consequences": [
                    motif("ENSM00000588223", "ENSPFM0427", 5, "Y", -0.094, -1, &["NRF1"], "C"),
                    motif("ENSM00000588224", "ENSPFM0427", 8, "Y", -0.094, 1, &["NRF1"], "C"),
                    motif("ENSM00000588225", "ENSPFM0427", 11, "N", -0.070, -1, &["NRF1"], "C"),
                    motif("ENSM00000588226", "ENSPFM0427", 2, "N", -0.082, 1, &["NRF1"], "C"),
                ],
                "regulatory_feature_consequences": [{"biotype": "enhancer", "consequence_terms": ["regulatory_region_variant"], "impact": "MODIFIER", "regulatory_feature_id": "ENSR21_C6M55", "variant_allele": "C"}],
                "seq_region_name": "21", "start": 46699255, "strand": 1
            }),
        ]
    );
}

/// Multi-base alleles over the motifs of `fixtures/regulatory/chr21-slice.json`,
/// beside the rows Ensembl VEP 116.3 writes for the same records on the same
/// cache. On a forward-strand motif the two engines agree: `21:5090720 GG>AA`
/// scores `-0.108` on ENSM00000574814, `21:5090723 CA>TT` `-0.111` on it, and
/// `21:46699250 CG>TA` `-0.149` on ENSM00000588224. On a reverse-strand motif
/// vep-rs splices the reverse-complemented allele at the positions it covers,
/// while Ensembl VEP splices it one position too far along the motif
/// (`MotifFeatureVariationAllele.pm:245` starts at `motif_start`, the allele's
/// last motif position): `21:5090730 CT>AG` scores `-0.136` here on
/// ENSM00000574815 where Ensembl VEP writes `-0.092`, `21:46699250 CG>TA`
/// `-0.143` on ENSM00000588223 where it writes `-0.126`, and `21:5090723
/// CA>TT`, which ends the motif, `0.044` where Ensembl VEP's splice lengthens
/// the sequence and writes no score. MOTIF_POS and HIGH_INF_POS agree on
/// every row.
#[test]
fn reverse_strand_mnv_scores_where_the_allele_lies() {
    let dir = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache_of_real_features(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t5090720\tfwd_mnv_inside\tGG\tAA\t.\t.\t.\n\
         21\t5090723\trev_mnv_start\tCA\tTT\t.\t.\t.\n\
         21\t5090730\trev_mnv_inside\tCT\tAG\t.\t.\t.\n\
         21\t46699250\tnrf1_rev_mnv\tCG\tTA\t.\t.\t.\n",
    )
    .unwrap();
    let motif_rows: Vec<String> = annotate(dir.path(), &cache, "default", true)
        .lines()
        .filter(|line| line.contains("\tMotifFeature\t"))
        .map(str::to_owned)
        .collect();
    let row = |id: &str, location: &str, allele: &str, feature: &str, extra: &str| {
        format!(
            "{id}\t{location}\t{allele}\t-\t{feature}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;{extra}"
        )
    };
    assert_eq!(
        motif_rows,
        [
            row("fwd_mnv_inside", "21:5090720-5090721", "AA", "ENSM00000574814", "STRAND=1;MOTIF_NAME=ENSPFM0290;MOTIF_POS=3;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.108;TRANSCRIPTION_FACTORS=HNF4A,NR2F1,RXRA"),
            row("rev_mnv_start", "21:5090723-5090724", "TT", "ENSM00000574814", "STRAND=1;MOTIF_NAME=ENSPFM0290;MOTIF_POS=6;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.111;TRANSCRIPTION_FACTORS=HNF4A,NR2F1,RXRA"),
            row("rev_mnv_start", "21:5090723-5090724", "TT", "ENSM00000574815", "STRAND=-1;MOTIF_NAME=ENSPFM0484;MOTIF_POS=16;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=0.044;TRANSCRIPTION_FACTORS=RUNX3"),
            row("rev_mnv_inside", "21:5090730-5090731", "AG", "ENSM00000574814", "STRAND=1;MOTIF_NAME=ENSPFM0290;MOTIF_POS=13;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.016;TRANSCRIPTION_FACTORS=HNF4A,NR2F1,RXRA"),
            row("rev_mnv_inside", "21:5090730-5090731", "AG", "ENSM00000574815", "STRAND=-1;MOTIF_NAME=ENSPFM0484;MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.136;TRANSCRIPTION_FACTORS=RUNX3"),
            row("nrf1_rev_mnv", "21:46699250-46699251", "TA", "ENSM00000588223", "STRAND=-1;MOTIF_NAME=ENSPFM0427;MOTIF_POS=10;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.143;TRANSCRIPTION_FACTORS=NRF1"),
            row("nrf1_rev_mnv", "21:46699250-46699251", "TA", "ENSM00000588224", "STRAND=1;MOTIF_NAME=ENSPFM0427;MOTIF_POS=3;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.149;TRANSCRIPTION_FACTORS=NRF1"),
        ]
    );
}

/// The motif rows of a run over the real-feature cache on `input`, written as
/// `text` and read under `extra_args`, the input format among them.
fn motif_rows_of(
    dir: &Path,
    cache: &Path,
    input: &str,
    text: &str,
    extra_args: &[&str],
) -> Vec<String> {
    fs::write(dir.join(input), text).unwrap();
    let output = dir.join(format!("{input}.out"));
    let mut command = common::isolated_vep(dir);
    command
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--force_overwrite",
            "--regulatory",
        ])
        .arg("--input_file")
        .arg(dir.join(input))
        .arg("--output_file")
        .arg(&output)
        .args(extra_args);
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    fs::read_to_string(output)
        .unwrap()
        .lines()
        .filter(|line| line.contains("\tMotifFeature\t"))
        .map(str::to_owned)
        .collect()
}

/// The two motif rows Ensembl VEP writes for the change A>T at 21:5090724,
/// under the record name and `Allele` cell given: position 7 of the forward
/// motif ENSM00000574814 and position 15 of the reverse motif ENSM00000574815.
fn motif_rows_at_5090724(records: &[(&str, &str)]) -> Vec<String> {
    let row = |id: &str, allele: &str, feature: &str, extra: &str| {
        format!(
            "{id}\t21:5090724\t{allele}\t-\t{feature}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;{extra}"
        )
    };
    records
        .iter()
        .flat_map(|(id, allele)| {
            [
                row(id, allele, "ENSM00000574814", "STRAND=1;MOTIF_NAME=ENSPFM0290;MOTIF_POS=7;HIGH_INF_POS=Y;MOTIF_SCORE_CHANGE=-0.089;TRANSCRIPTION_FACTORS=HNF4A,NR2F1,RXRA"),
                row(id, allele, "ENSM00000574815", "STRAND=-1;MOTIF_NAME=ENSPFM0484;MOTIF_POS=15;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=0.039;TRANSCRIPTION_FACTORS=RUNX3"),
            ]
        })
        .collect()
}

/// An Ensembl-format record given on the reverse strand (`-` or `-1`) names
/// the genomic change its forward twin names, and its motif rows score that
/// change: Ensembl VEP's `feature_seq` complements the allele when the record's
/// strand differs from the motif's (`VariationFeatureOverlapAllele.pm:253-256`),
/// so `T/A -` at a forward `A` scores A>T on the forward motif and on the
/// reverse one, with the `Allele` cell as written. The rows are Ensembl VEP's
/// on this input.
#[test]
fn reverse_strand_ensembl_record_scores_the_genomic_change() {
    let dir = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache_of_real_features(dir.path());
    let rows = motif_rows_of(
        dir.path(),
        &cache,
        "in.txt",
        "21\t5090724\t5090724\tA/T\t+\tfwd_strand_input\n\
         21\t5090724\t5090724\tT/A\t-\trev_strand_input\n\
         21\t5090724\t5090724\tT/A\t-1\trev_strand_minus_one\n",
        &["--format", "ensembl", "--no_fasta"],
    );
    assert_eq!(
        rows,
        motif_rows_at_5090724(&[
            ("fwd_strand_input", "T"),
            ("rev_strand_input", "A"),
            ("rev_strand_minus_one", "A"),
        ])
    );
}

/// A region on strand `-1` scores as its forward twins do, the reference base
/// read from the FASTA; the rows are Ensembl VEP's on this input.
#[test]
fn reverse_strand_region_scores_the_genomic_change() {
    let dir = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache_of_real_features(dir.path());
    let fasta = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-regulatory/reference.fa.gz");
    let rows = motif_rows_of(
        dir.path(),
        &cache,
        "in.region.txt",
        "21:5090724-5090724:1/T\n21:5090724-5090724:-1/A\n21:5090724-5090724/T\n",
        &["--format", "region", "--fasta", fasta.to_str().unwrap()],
    );
    assert_eq!(
        rows,
        motif_rows_at_5090724(&[
            ("21:5090724-5090724:1/T", "T"),
            ("21:5090724-5090724:-1/A", "A"),
            ("21:5090724-5090724/T", "T"),
        ])
    );
}

/// An HGVS description on a reverse-strand transcript names its alleles on
/// that strand, as Ensembl VEP's parser builds the variant, so `c.81G>A` on
/// ENST00000303775 is the genomic change 21:28885265 C>T and its motif rows
/// are the rows Ensembl VEP writes for the VCF record `21 28885265 C T`: the
/// forward motif ENSM00000579187 and the reverse motif ENSM00000579188, both
/// ENSPFM0365, position 9, score change 0.002; the `Allele` cell is `A` as
/// written. The run reads the regulatory corpus's cache for the transcript.
#[test]
fn reverse_strand_transcript_hgvs_scores_the_genomic_change() {
    let dir = tempfile::tempdir().unwrap();
    let corpus =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-regulatory");
    let fasta = corpus.join("reference.fa.gz");
    let rows = motif_rows_of(
        dir.path(),
        &corpus.join("json_cache"),
        "in.hgvs.txt",
        "ENST00000303775.10:c.81G>A\n",
        &["--format", "hgvs", "--fasta", fasta.to_str().unwrap()],
    );
    let row = |feature: &str, strand: &str| {
        format!(
            "ENST00000303775.10:c.81G>A\t21:28885265\tA\t-\t{feature}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM0365;MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=0.002;TRANSCRIPTION_FACTORS=MAX"
        )
    };
    assert_eq!(
        rows,
        [row("ENSM00000579187", "1"), row("ENSM00000579188", "-1")]
    );
}

/// The Parquet output holds the regulatory, motif and intergenic rows the
/// `--tab` output holds, read back row by row through the cell contract and
/// compared as a multiset (the file's row order is the sort order). A motif
/// row's factor list comes back as the tab cell, and a factorless motif's
/// empty factor cell, a written value the tab format leaves empty rather than
/// `-`, is stored as the empty string and comes back empty.
#[cfg(feature = "parquet")]
#[test]
fn regulatory_rows_round_trip_through_parquet() {
    use vep_io::output::parquet::{canonical_tab_text, tab_column_types, tab_rows};
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t106\tforward\tC\tG\t.\t.\t.\n\
         21\t116\treverse\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let tab = annotate(dir.path(), &cache, "tab", true);
    let columns: Vec<&str> = tab
        .lines()
        .find(|line| line.starts_with('#') && !line.starts_with("##"))
        .expect("the tab column line")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let parquet = dir.path().join("out.parquet");
    let result = common::isolated_vep(dir.path())
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--no_fasta",
            "--force_overwrite",
            "--format",
            "vcf",
            "--regulatory",
            "--output_format",
            "parquet",
        ])
        .arg("--input_file")
        .arg(dir.path().join("in.vcf"))
        .arg("--output_file")
        .arg(&parquet)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "parquet: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    // Each tab cell in the form of the type the part stores its column under.
    let stored = tab_column_types(&parquet).unwrap();
    assert_eq!(
        stored
            .iter()
            .map(|(name, _)| name.as_str())
            .collect::<Vec<_>>(),
        columns
    );
    let mut expected: Vec<String> = tab
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            line.split('\t')
                .zip(&stored)
                .map(|(cell, (_, column_type))| canonical_tab_text(*column_type, cell))
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    expected.sort_unstable();
    assert!(
        expected.iter().any(|row| row.contains("\tTF_A,TF_B")),
        "the forward motif row carries its factors:\n{tab}"
    );
    assert!(
        expected
            .iter()
            .any(|row| row.contains("ENSM_REVERSE") && row.ends_with('\t')),
        "the factorless motif row ends in the empty factor cell:\n{tab}"
    );
    let mut actual = tab_rows(&parquet).unwrap();
    actual.sort_unstable();
    assert_eq!(actual, expected, "parquet round trip");
}

/// Every structural class and the sequence-variant edges over the hand-built
/// fixture (enhancer 100-120, motifs 105-108 and 115-118), as Ensembl VEP
/// 116.3 writes them: `<DEL>` and `<CN0>` covering a feature are its ablation
/// (`deletion`, `copy_number_loss`: Utils/VariationEffect.pm:192-216), `<DUP>`,
/// `<CN2>`, `<DUP:TANDEM>` and `<CNV:TR>` its amplification (`copy_number_gain`,
/// :246-256), `<CNV>` and `<INV>` the overlap term alone (a VCF `<CNV>` has no
/// supporting variants), each with `OverlapBP` and `OverlapPC`
/// (OutputFactory.pm:2096-2104) and none of `BIOTYPE`, `STRAND` and the motif
/// fields (:2078-2111 writes no such key); `--overlaps` adds the two as columns
/// to the tab and VCF outputs, `-` and empty on the intergenic row. A sequence
/// deletion covering the features is their ablation with the motif fields and
/// no `MOTIF_POS` (its start lies before the matrix); a multi-allelic SNV writes
/// each allele's rows; an insertion overlaps a feature only with both flanks
/// inside it (InputBuffer.pm:324): 100-101 and 104-105 overlap the enhancer and
/// not the motif at 105-108, 99-100 and 120-121 overlap nothing.
#[test]
fn sv_and_indel_classes_over_a_regulatory_feature() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let symbolic = [
        "<DEL>",
        "<DUP>",
        "<CNV>",
        "<INV>",
        "<CN0>",
        "<CN2>",
        "<CNV:TR>",
        "<DUP:TANDEM>",
    ];
    let mut input = header.to_string();
    for alt in symbolic {
        input.push_str(&format!("21\t94\t{alt}\tN\t{alt}\t.\t.\tEND=125\n"));
    }
    input.push_str("21\t104\tpartial\tN\t<DEL>\t.\t.\tEND=112\n");
    input.push_str("21\t99\tcover_del\tAAAAAAAAAAAAAAAAAAAAAA\tA\t.\t.\t.\n");
    input.push_str("21\t106\tmulti\tC\tG,T\t.\t.\t.\n");
    for (id, pos) in [
        ("ins_in", 100),
        ("ins_before", 99),
        ("ins_after", 120),
        ("ins_motif", 104),
    ] {
        input.push_str(&format!("21\t{pos}\t{id}\tA\tAT\t.\t.\t.\n"));
    }
    fs::write(dir.path().join("in.vcf"), input).unwrap();

    // (id, class term, regulatory terms, motif terms, motif impact)
    let classes = [
        (
            "<DEL>",
            "deletion",
            "regulatory_region_ablation,regulatory_region_variant",
            "TFBS_ablation,TF_binding_site_variant",
            "MODERATE",
        ),
        (
            "<DUP>",
            "duplication",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<CNV>",
            "copy_number_variation",
            "regulatory_region_variant",
            "TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<INV>",
            "inversion",
            "regulatory_region_variant",
            "TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<CN0>",
            "deletion",
            "regulatory_region_ablation,regulatory_region_variant",
            "TFBS_ablation,TF_binding_site_variant",
            "MODERATE",
        ),
        (
            "<CN2>",
            "duplication",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<CNV:TR>",
            "tandem_repeat",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<DUP:TANDEM>",
            "tandem_duplication",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
    ];
    let mut expected = Vec::new();
    for (id, class, regulatory, motif, motif_impact) in classes {
        let row = |feature: &str, kind: &str, terms: &str, extra: &str| {
            format!("{id}\t21:95-125\t{class}\t-\t{feature}\t{kind}\t{terms}\t-\t-\t-\t-\t-\t-\t{extra}")
        };
        expected.push(row(
            "ENSR_FIXTURE",
            "RegulatoryFeature",
            regulatory,
            "IMPACT=MODIFIER;OverlapBP=21;OverlapPC=100.00",
        ));
        for feature in ["ENSM_FORWARD", "ENSM_REVERSE"] {
            expected.push(row(
                feature,
                "MotifFeature",
                motif,
                &format!("IMPACT={motif_impact};OverlapBP=4;OverlapPC=100.00"),
            ));
        }
        expected.push(row("-", "-", "intergenic_variant", "IMPACT=MODIFIER"));
    }
    expected.extend([
        "partial\t21:105-112\tdeletion\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;OverlapBP=8;OverlapPC=38.10".to_string(),
        "partial\t21:105-112\tdeletion\t-\tENSM_FORWARD\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;OverlapBP=4;OverlapPC=100.00".to_string(),
        "partial\t21:105-112\tdeletion\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "cover_del\t21:100-120\t-\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_ablation,regulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "cover_del\t21:100-120\t-\t-\tENSM_FORWARD\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;STRAND=1;MOTIF_NAME=ENSPFM_FORWARD;HIGH_INF_POS=N;TRANSCRIPTION_FACTORS=TF_A,TF_B".to_string(),
        "cover_del\t21:100-120\t-\t-\tENSM_REVERSE\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;STRAND=-1;MOTIF_NAME=ENSPFM_REVERSE;HIGH_INF_POS=N;TRANSCRIPTION_FACTORS=".to_string(),
        "cover_del\t21:100-120\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "multi\t21:106\tG\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "multi\t21:106\tT\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "multi\t21:106\tG\t-\tENSM_FORWARD\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM_FORWARD;MOTIF_POS=2;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.250;TRANSCRIPTION_FACTORS=TF_A,TF_B".to_string(),
        "multi\t21:106\tT\t-\tENSM_FORWARD\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM_FORWARD;MOTIF_POS=2;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.250;TRANSCRIPTION_FACTORS=TF_A,TF_B".to_string(),
        "multi\t21:106\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "multi\t21:106\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_in\t21:100-101\tT\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "ins_in\t21:100-101\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_before\t21:99-100\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_after\t21:120-121\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_motif\t21:104-105\tT\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "ins_motif\t21:104-105\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
    ]);
    let rows: Vec<String> = annotate(dir.path(), &cache, "default", true)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    assert_eq!(rows, expected);

    // --overlaps: the two values become columns of the tab and VCF outputs, after
    // FLAGS and before BIOTYPE, `-` and empty on the intergenic row.
    let tab = annotate_with_args(dir.path(), &cache, "tab", true, &["--overlaps"]);
    let header = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap();
    assert_eq!(
        header,
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tOverlapBP\tOverlapPC\tBIOTYPE\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tTRANSCRIPTION_FACTORS"
    );
    let partial: Vec<&str> = tab.lines().filter(|l| l.starts_with("partial\t")).collect();
    assert_eq!(
        partial,
        [
            "partial\t21:105-112\tdeletion\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t8\t38.10\t-\t-\t-\t-\t-\t-",
            "partial\t21:105-112\tdeletion\t-\tENSM_FORWARD\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tMODERATE\t-\t-\t-\t4\t100.00\t-\t-\t-\t-\t-\t-",
            "partial\t21:105-112\tdeletion\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
        ]
    );
    let vcf = annotate_with_args(dir.path(), &cache, "vcf", true, &["--overlaps"]);
    let format = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .unwrap();
    assert_eq!(
        format,
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|OverlapBP|OverlapPC|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS\">"
    );
    let csq = |id: &str| -> Vec<String> {
        vcf.lines()
            .find(|l| l.split('\t').nth(2) == Some(id))
            .and_then(|l| l.split('\t').nth(7))
            .and_then(|info| info.split(';').find_map(|item| item.strip_prefix("CSQ=")))
            .unwrap()
            .split(',')
            .map(str::to_owned)
            .collect()
    };
    assert_eq!(
        csq("partial"),
        [
            "deletion|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR_FIXTURE|||||||||||||||||8|38.10|||||",
            "deletion|TFBS_ablation&TF_binding_site_variant|MODERATE|||MotifFeature|ENSM_FORWARD|||||||||||||||||4|100.00|||||",
            "deletion|intergenic_variant|MODIFIER|||||||||||||||||||||||||||",
        ]
    );
    assert_eq!(
        csq("cover_del"),
        [
            "-|regulatory_region_ablation&regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR_FIXTURE|enhancer||||||||||||||||||||||",
            "-|TFBS_ablation&TF_binding_site_variant|MODERATE|||MotifFeature|ENSM_FORWARD|||||||||||||1||||||ENSPFM_FORWARD||N||TF_A&TF_B",
            "-|TFBS_ablation&TF_binding_site_variant|MODERATE|||MotifFeature|ENSM_REVERSE|||||||||||||-1||||||ENSPFM_REVERSE||N||",
            "-|intergenic_variant|MODIFIER|||||||||||||||||||||||||||",
        ]
    );
}

/// The pick family through the command line, on a record inside two
/// transcripts of one gene and the fixture's features shifted onto it, and on a
/// record inside the same features shifted between the genes: `--pick` keeps
/// the one worst allele of the record (the intron row; the motif row where no
/// transcript overlaps, the intergenic row then dropped); `--per_gene` and
/// `--pick_allele_gene` keep one transcript row per gene and every regulatory,
/// motif and intergenic row (OutputFactory.pm:837-845, non-transcript alleles
/// pass through), written before the gene picks as the per-gene list returns
/// them (`:838-856`); `--flag_pick` and `--flag_pick_allele` flag the one pick;
/// `--flag_pick_allele_gene` flags the per-gene pick and every pass-through row,
/// the intergenic row included (:623-627); `--most_severe` and `--summary` read
/// the regulatory terms with the transcript terms.
#[test]
fn pick_family_flags_through_the_cli_keep_regulatory_rows() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    fs::create_dir_all(cache.join("transcripts/21")).unwrap();
    fs::create_dir_all(cache.join("regulatory/21")).unwrap();
    // Two protein-coding transcripts of ENSG00000166351 at 13609777-13645823 and
    // three lncRNA transcripts of another gene at 13843132-13897298.
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join(
            "../../tests/golden/116/GRCh38/json_cache/transcripts/21/13000001-14000000.json.gz",
        ),
        cache.join("transcripts/21/13000001-14000000.json.gz"),
    )
    .unwrap();
    let fixture: Vec<Value> =
        serde_json::from_slice(include_bytes!("fixtures/regulatory/1-1000000.json")).unwrap();
    let shifted: Vec<Value> = [(13_619_894u64, "_IN_GENE"), (13_699_894, "_BETWEEN")]
        .iter()
        .flat_map(|(offset, suffix)| {
            fixture.iter().map(move |feature| {
                let mut feature = feature.clone();
                for key in ["start", "end"] {
                    let value = feature[key].as_u64().unwrap() + offset;
                    feature[key] = Value::from(value);
                }
                let id = format!("{}{suffix}", feature["stable_id"].as_str().unwrap());
                feature["stable_id"] = Value::from(id);
                feature
            })
        })
        .collect();
    fs::write(
        cache.join("regulatory/21/13000001-14000000.json"),
        serde_json::to_vec(&shifted).unwrap(),
    )
    .unwrap();
    common::write_regulatory_matrices(&cache);
    fs::write(cache.join("info.json"), common::regulatory_info_json()).unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t13620000\tin_gene\tC\tG\t.\t.\t.\n\
         21\t13700000\tbetween\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    // (record, Feature, Feature_type, Consequence, PICK)
    let rows = |flags: &[&str]| -> Vec<(String, String, String, String, String)> {
        annotate_with_args(dir.path(), &cache, "default", true, flags)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                let pick = cols[13]
                    .split(';')
                    .find_map(|item| item.strip_prefix("PICK="))
                    .unwrap_or("-");
                (
                    cols[0].to_string(),
                    cols[4].to_string(),
                    cols[5].to_string(),
                    cols[6].to_string(),
                    pick.to_string(),
                )
            })
            .collect()
    };
    let row = |record: &str, feature: &str, kind: &str, terms: &str, pick: &str| {
        (
            record.to_string(),
            feature.to_string(),
            kind.to_string(),
            terms.to_string(),
            pick.to_string(),
        )
    };
    let in_gene_t1 = |pick| {
        row(
            "in_gene",
            "ENST00000299443",
            "Transcript",
            "intron_variant",
            pick,
        )
    };
    let in_gene_t2 = |pick| {
        row(
            "in_gene",
            "ENST00000620442",
            "Transcript",
            "intron_variant",
            pick,
        )
    };
    let in_gene_rf = |pick| {
        row(
            "in_gene",
            "ENSR_FIXTURE_IN_GENE",
            "RegulatoryFeature",
            "regulatory_region_variant",
            pick,
        )
    };
    let in_gene_mf = |pick| {
        row(
            "in_gene",
            "ENSM_FORWARD_IN_GENE",
            "MotifFeature",
            "TF_binding_site_variant",
            pick,
        )
    };
    let between_rf = |pick| {
        row(
            "between",
            "ENSR_FIXTURE_BETWEEN",
            "RegulatoryFeature",
            "regulatory_region_variant",
            pick,
        )
    };
    let between_mf = |pick| {
        row(
            "between",
            "ENSM_FORWARD_BETWEEN",
            "MotifFeature",
            "TF_binding_site_variant",
            pick,
        )
    };
    let between_ig = |pick| row("between", "-", "-", "intergenic_variant", pick);

    assert_eq!(
        rows(&[]),
        [
            in_gene_t1("-"),
            in_gene_t2("-"),
            in_gene_rf("-"),
            in_gene_mf("-"),
            between_rf("-"),
            between_mf("-"),
            between_ig("-")
        ]
    );
    assert_eq!(rows(&["--pick"]), [in_gene_t1("-"), between_mf("-")]);
    for flag in ["--per_gene", "--pick_allele_gene"] {
        assert_eq!(
            rows(&[flag]),
            [
                in_gene_rf("-"),
                in_gene_mf("-"),
                in_gene_t1("-"),
                between_rf("-"),
                between_mf("-"),
                between_ig("-")
            ],
            "{flag}"
        );
    }
    for flag in ["--flag_pick", "--flag_pick_allele"] {
        assert_eq!(
            rows(&[flag]),
            [
                in_gene_t1("1"),
                in_gene_t2("-"),
                in_gene_rf("-"),
                in_gene_mf("-"),
                between_rf("-"),
                between_mf("1"),
                between_ig("-")
            ],
            "{flag}"
        );
    }
    assert_eq!(
        rows(&["--flag_pick_allele_gene"]),
        [
            in_gene_t1("1"),
            in_gene_t2("-"),
            in_gene_rf("1"),
            in_gene_mf("1"),
            between_rf("1"),
            between_mf("1"),
            between_ig("1")
        ]
    );
    assert_eq!(
        rows(&["--most_severe"]),
        [
            row("in_gene", "-", "-", "intron_variant", "-"),
            row("between", "-", "-", "TF_binding_site_variant", "-")
        ]
    );
    assert_eq!(
        rows(&["--summary"]),
        [
            row(
                "in_gene",
                "-",
                "-",
                "intron_variant,TF_binding_site_variant,regulatory_region_variant",
                "-"
            ),
            row(
                "between",
                "-",
                "-",
                "TF_binding_site_variant,regulatory_region_variant,intergenic_variant",
                "-"
            ),
        ]
    );
}

/// Ensembl VEP matches a record to the regulatory features its span as read
/// overlaps (InputBuffer.pm:324, with `:328-332` reading an unshifted span no
/// code sets, so the one span tested is the input's), and `--shift_3prime` moves
/// no regulatory row: a deletion of the forward motif's last two bases, which a
/// 3' shift carries two bases past the motif where the reference repeats them,
/// keeps its regulatory and motif rows at the input position, while under
/// `--shift_genomic 1` the intergenic row prints the shifted location, as
/// `IntergenicVariationAllele_to_output_hash` shifts it
/// (OutputFactory.pm:2001-2007) and the regulatory factories do not.
#[test]
fn unshifted_span_is_tested_under_shift_3prime() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    // The fixture's forward motif spans 105..=108 inside the regulatory feature
    // 100..=120. Bases 107..=110 read `AAAA`, so the deletion of 107..=108 shifts
    // two bases to the right, past the motif.
    let mut sequence = vec![b'C'; 200];
    sequence[106..110].copy_from_slice(b"AAAA");
    let fasta = dir.path().join("21.fa");
    let mut contents = b">21\n".to_vec();
    contents.extend_from_slice(&sequence);
    contents.push(b'\n');
    fs::write(&fasta, contents).unwrap();
    fs::write(dir.path().join("21.fa.fai"), "21\t200\t4\t200\t201\n").unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t106\tedge\tCAA\tC\t.\t.\t.\n",
    )
    .unwrap();
    let rows = |extra: &[&str]| -> Vec<(String, String, String)> {
        let mut args = vec!["--fasta", fasta.to_str().unwrap()];
        args.extend_from_slice(extra);
        annotate_with_args(dir.path(), &cache, "default", true, &args)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                (
                    cols[1].to_string(),
                    cols[2].to_string(),
                    cols[4].to_string(),
                )
            })
            .collect()
    };
    let row = |location: &str, allele: &str, feature: &str| {
        (
            location.to_string(),
            allele.to_string(),
            feature.to_string(),
        )
    };
    let plain = rows(&[]);
    assert_eq!(
        plain,
        [
            row("21:107-108", "-", "ENSR_FIXTURE"),
            row("21:107-108", "-", "ENSM_FORWARD"),
            row("21:107-108", "-", "-"),
        ]
    );
    assert_eq!(rows(&["--shift_3prime", "1"]), plain);
    assert_eq!(
        rows(&["--shift_3prime", "1", "--shift_genomic", "1"]),
        [
            row("21:107-108", "-", "ENSR_FIXTURE"),
            row("21:107-108", "-", "ENSM_FORWARD"),
            row("21:109-110", "-", "-"),
        ]
    );
}

/// A regulatory shard the run reads must be a JSON array of regulatory and
/// motif features: a truncated shard, a feature type other than the two (a
/// transcript among them), an element without the feature fields and a motif
/// without its binding matrix each fail the run with exit 1 and the shard's
/// path in the error, writing no data row; an empty array is a chromosome
/// with no features, so its records write their intergenic rows; a
/// gzip-compressed shard loads as a plain one does.
#[test]
fn malformed_regulatory_shards_fail_the_run_by_name() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tr\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let run = |shard: &[u8], name: &str| -> (std::process::ExitStatus, String, Vec<String>) {
        let cache = dir.path().join("cache");
        let _ = fs::remove_dir_all(&cache);
        fs::create_dir_all(cache.join("transcripts/21")).unwrap();
        fs::create_dir_all(cache.join("regulatory/21")).unwrap();
        fs::write(cache.join("transcripts/21/1-1000000.json"), "[]").unwrap();
        fs::write(cache.join("regulatory/21").join(name), shard).unwrap();
        common::write_regulatory_matrices(&cache);
        fs::write(cache.join("info.json"), common::regulatory_info_json()).unwrap();
        let out = dir.path().join("out.txt");
        let _ = fs::remove_file(&out);
        let result = common::isolated_vep(dir.path())
            .args(["--json_cache", cache.to_str().unwrap()])
            .args([
                "--assembly",
                "GRCh38",
                "--no_fasta",
                "--no_stats",
                "--force_overwrite",
            ])
            .arg("--regulatory")
            .arg("--input_file")
            .arg(dir.path().join("in.vcf"))
            .arg("--output_file")
            .arg(&out)
            .output()
            .unwrap();
        let rows = fs::read_to_string(&out)
            .unwrap_or_default()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_owned)
            .collect();
        (
            result.status,
            String::from_utf8_lossy(&result.stderr).into_owned(),
            rows,
        )
    };
    for (shard, phrase) in [
        (&b"["[..], "EOF while parsing"),
        (
            br#"[{"stable_id":"X","start":100,"end":120,"strand":0,"feature_type":"Enhancer"}]"#,
            "unknown variant `Enhancer`, expected `RegulatoryFeature` or `MotifFeature`",
        ),
        (
            br#"[{"stable_id":"T1","start":100,"end":120,"strand":1,"feature_type":"Transcript"}]"#,
            "unknown variant `Transcript`",
        ),
        (br#"[{"foo":1}]"#, "missing field `feature_type`"),
        (
            br#"[{"foo":1,"feature_type":"RegulatoryFeature"}]"#,
            "missing field `stable_id`",
        ),
        (
            br#"[{"stable_id":"ENSM_X","start":100,"end":120,"strand":1,"feature_type":"MotifFeature"}]"#,
            "missing field `matrix`",
        ),
        (
            br#"[{"stable_id":"ENSM_X","start":100,"end":120,"strand":1,"feature_type":"MotifFeature","matrix":"ENSPFM_ABSENT"}]"#,
            "motif feature ENSM_X names the binding matrix ENSPFM_ABSENT, which binding_matrices.json lacks",
        ),
        (
            br#"[{"stable_id":"ENSR_X","start":100,"end":120,"feature_type":"RegulatoryFeature","activity":"A"}]"#,
            "regulatory feature ENSR_X carries 1 activity codes where the cache lists 112 cell types",
        ),
    ] {
        let (status, stderr, rows) = run(shard, "1-1000000.json");
        assert_eq!(status.code(), Some(1), "{stderr}");
        assert!(
            stderr.contains("ERROR: [cache] parsing")
                && stderr.contains("regulatory/21/1-1000000.json")
                && stderr.contains(phrase),
            "{stderr}"
        );
        assert!(rows.is_empty(), "{rows:?}");
    }
    let (status, stderr, rows) = run(b"[]", "1-1000000.json");
    assert!(status.success(), "{stderr}");
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].contains("\t-\t-\t-\tintergenic_variant\t"),
        "{}",
        rows[0]
    );
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(
        &mut encoder,
        include_bytes!("fixtures/regulatory/1-1000000.json"),
    )
    .unwrap();
    let (status, stderr, rows) = run(&encoder.finish().unwrap(), "1-1000000.json.gz");
    assert!(status.success(), "{stderr}");
    assert_eq!(rows.len(), 3, "{rows:?}");
    assert!(
        rows[0].contains("\tENSR_FIXTURE\tRegulatoryFeature\t"),
        "{}",
        rows[0]
    );
}

#[test]
fn regulatory_and_motif_rows_use_dash_for_absent_gene() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    for format in ["default", "tab"] {
        let text = annotate(dir.path(), &cache, format, true);
        let features: Vec<_> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|columns| matches!(columns[5], "RegulatoryFeature" | "MotifFeature"))
            .map(|columns| (columns[4], columns[3]))
            .collect();
        assert_eq!(
            features,
            [("ENSR_FIXTURE", "-"), ("ENSM_FORWARD", "-")],
            "{format}: {text}"
        );
    }
}

/// `--cell_type` prints `CELL_TYPE` on every regulatory and motif row in every
/// text format: a regulatory row's list is the requested epigenomes projected
/// onto the enhancer's activity in the order given (a name the feature's map
/// lacks, `psoas_muscle_(m,_3_y)`, skipped; a name holding a comma one
/// epigenome), a motif row's list is empty and still a value (`CELL_TYPE=` in
/// the default Extra, an empty tab cell, an empty CSQ subfield, `[]` in JSON),
/// and the header describes the column. The flag switches `--regulatory` on by
/// itself; a name the cache lacks ends the run with Ensembl VEP's error (exit
/// 1); on a cache without regulatory data the names are not checked, as Ensembl
/// VEP builds no regulatory source there and checks none.
#[test]
fn cell_type_projection_reaches_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(dir.path().join("in.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n").unwrap();
    let names = "GM12878_(f),astrocyte,neuronal_stem_cell_(m),psoas_muscle_(m,_3_y),\
                 placenta_(85-113_dpf),prostate_gland_(m,_37_y)";
    let list = [
        "GM12878_(f):ACTIVE",
        "astrocyte:INACTIVE",
        "neuronal_stem_cell_(m):POISED",
        "placenta_(85-113_dpf):REPRESSED",
        "prostate_gland_(m,_37_y):NA",
    ];
    let rows = |text: &str, column: usize| -> Vec<(String, String)> {
        text.lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|columns| matches!(columns[5], "RegulatoryFeature" | "MotifFeature"))
            .map(|columns| (columns[4].to_string(), columns[column].to_string()))
            .collect()
    };

    let text = annotate_with_args(
        dir.path(),
        &cache,
        "default",
        false,
        &["--cell_type", names],
    );
    let extras = rows(&text, 13);
    assert_eq!(extras.len(), 2, "{text}");
    assert_eq!(extras[0].0, "ENSR_FIXTURE");
    assert!(
        extras[0]
            .1
            .ends_with(&format!(";BIOTYPE=enhancer;CELL_TYPE={}", list.join(","))),
        "{text}"
    );
    assert_eq!(extras[1].0, "ENSM_FORWARD");
    assert!(
        extras[1]
            .1
            .ends_with(";TRANSCRIPTION_FACTORS=TF_A,TF_B;CELL_TYPE="),
        "{text}"
    );
    assert!(
        text.contains(
            "## CELL_TYPE : List of cell types and classifications for regulatory feature"
        ),
        "{text}"
    );

    let text = annotate_with_args(dir.path(), &cache, "tab", true, &["--cell_type", names]);
    let header = text.lines().find(|line| line.starts_with("#U")).unwrap();
    assert_eq!(header.split('\t').nth(23), Some("CELL_TYPE"), "{header}");
    assert_eq!(
        rows(&text, 23),
        [
            ("ENSR_FIXTURE".to_string(), list.join(",")),
            ("ENSM_FORWARD".to_string(), String::new()),
        ],
        "{text}"
    );
    let without = annotate_with_args(dir.path(), &cache, "tab", true, &[]);
    let header = without.lines().find(|line| line.starts_with("#U")).unwrap();
    assert_eq!(header.split('\t').count(), 23, "{header}");
    assert!(!without.contains("CELL_TYPE"), "{without}");

    let text = annotate_with_args(dir.path(), &cache, "vcf", true, &["--cell_type", names]);
    let format_line = text
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ"))
        .unwrap();
    assert!(
        format_line.ends_with("|TRANSCRIPTION_FACTORS|CELL_TYPE\">"),
        "{format_line}"
    );
    let csq = text
        .lines()
        .find(|line| !line.starts_with('#'))
        .and_then(|line| line.split('\t').nth(7))
        .and_then(|info| info.strip_prefix("CSQ="))
        .unwrap();
    let entries: Vec<Vec<&str>> = csq.split(',').map(|e| e.split('|').collect()).collect();
    let by_feature = |id: &str| entries.iter().find(|e| e.contains(&id)).unwrap().clone();
    let regulatory = by_feature("ENSR_FIXTURE");
    // The VCF writer joins the list with `&` and, as Ensembl VEP's does
    // (`OutputFactory/VCF.pm:401`), writes a value's own comma as `&` too.
    assert_eq!(
        *regulatory.last().unwrap(),
        list.join("&").replace(',', "&"),
        "{csq}"
    );
    let motif = by_feature("ENSM_FORWARD");
    assert_eq!(*motif.last().unwrap(), "", "{csq}");

    let json = annotate_with_args(dir.path(), &cache, "json", true, &["--cell_type", names]);
    let record: Value = serde_json::from_str(json.lines().next().unwrap()).unwrap();
    assert_eq!(
        record["regulatory_feature_consequences"][0]["cell_type"],
        serde_json::json!(list)
    );
    assert_eq!(
        record["motif_feature_consequences"][0]["cell_type"],
        serde_json::json!([])
    );

    let mut command = common::isolated_vep(dir.path());
    let result = command
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--no_fasta",
            "--force_overwrite",
            "--cell_type",
            "astrocyte,HUVEC",
            "--input_file",
        ])
        .arg(dir.path().join("in.vcf"))
        .arg("--output_file")
        .arg(dir.path().join("refused.txt"))
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains(
            "ERROR: [cli] Cell type \"HUVEC\" unavailable; available cell types are:\n\
             A549_(m,_58_y),A673_(f,_15_y),"
        ),
        "{stderr}"
    );

    // No regulatory data: the warning, no check, exit 0.
    fs::write(cache.join("info.json"), r#"{"assembly":"GRCh38"}"#).unwrap();
    let mut command = common::isolated_vep(dir.path());
    let result = command
        .current_dir(dir.path())
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--cell_type",
            "HUVEC",
            "--no_fasta",
            "--no_stats",
            "--tab",
            "--input_file",
        ])
        .arg(dir.path().join("in.vcf"))
        .args(["--output_file", "STDOUT"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("this cache holds no regulatory data"));
    let stdout = String::from_utf8_lossy(&result.stdout);
    let header = stdout.lines().find(|line| line.starts_with("#U")).unwrap();
    assert_eq!(header.split('\t').nth(23), Some("CELL_TYPE"), "{header}");
}

fn assert_regulatory_only_rows(extra_args: &[&str], regulatory: bool, intergenic: bool) {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t110\tregulatory_only\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    for format in ["default", "tab"] {
        let text = annotate_with_args(dir.path(), &cache, format, true, extra_args);
        let rows: Vec<_> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect();
        let mut expected = Vec::new();
        if regulatory {
            expected.push((
                "ENSR_FIXTURE",
                "RegulatoryFeature",
                "regulatory_region_variant",
            ));
        }
        if intergenic {
            expected.push(("-", "-", "intergenic_variant"));
        }
        assert_eq!(
            rows.iter()
                .map(|row| (row[4], row[5], row[6]))
                .collect::<Vec<_>>(),
            expected,
            "{format} {extra_args:?}: {text}"
        );
        for row in &rows {
            assert_eq!(row[3], "-");
            if row[5] == "RegulatoryFeature" {
                if format == "default" {
                    assert!(row[13].split(';').any(|item| item == "BIOTYPE=enhancer"));
                } else {
                    let header = text
                        .lines()
                        .find(|line| line.starts_with("#Uploaded_variation"))
                        .unwrap();
                    let biotype = header
                        .split('\t')
                        .position(|column| column == "BIOTYPE")
                        .unwrap();
                    assert_eq!(row[biotype], "enhancer");
                }
            }
        }
    }
    let text = annotate_with_args(dir.path(), &cache, "json", true, extra_args);
    let records: Vec<Value> = text
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert!(record.get("transcript_consequences").is_none());
    assert!(record.get("motif_feature_consequences").is_none());
    assert_eq!(
        record["most_severe_consequence"],
        if regulatory {
            "regulatory_region_variant"
        } else {
            "?"
        }
    );
    if regulatory {
        assert_eq!(
            record["regulatory_feature_consequences"],
            serde_json::json!([{
                "biotype": "enhancer",
                "consequence_terms": ["regulatory_region_variant"],
                "impact": "MODIFIER",
                "regulatory_feature_id": "ENSR_FIXTURE",
                "variant_allele": "G"
            }])
        );
    } else {
        assert!(record.get("regulatory_feature_consequences").is_none());
    }
    if intergenic {
        assert_eq!(
            record["intergenic_consequences"],
            serde_json::json!([{
                "consequence_terms": ["intergenic_variant"],
                "impact": "MODIFIER",
                "variant_allele": "G"
            }])
        );
    } else {
        assert!(record.get("intergenic_consequences").is_none());
    }
    let text = annotate_with_args(dir.path(), &cache, "vcf", true, extra_args);
    let rows: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 1);
    let info = rows[0].split('\t').nth(7).unwrap();
    if regulatory {
        let consequences: Vec<_> = info.strip_prefix("CSQ=").unwrap().split(',').collect();
        assert_eq!(consequences.len(), 1 + usize::from(intergenic));
        let regulatory_row: Vec<_> = consequences[0].split('|').collect();
        assert_eq!(
            &regulatory_row[..8],
            &[
                "G",
                "regulatory_region_variant",
                "MODIFIER",
                "",
                "",
                "RegulatoryFeature",
                "ENSR_FIXTURE",
                "enhancer"
            ]
        );
        if intergenic {
            let intergenic_row: Vec<_> = consequences[1].split('|').collect();
            assert_eq!(
                &intergenic_row[..8],
                &["G", "intergenic_variant", "MODIFIER", "", "", "", "", ""]
            );
        }
    } else {
        assert_eq!(info, ".");
    }
}

#[test]
fn regulatory_only_retains_intergenic_row() {
    assert_regulatory_only_rows(&[], true, true);
}

#[test]
fn regulatory_only_pick_retains_regulatory_row() {
    assert_regulatory_only_rows(&["--pick"], true, false);
}

#[test]
fn regulatory_only_coding_only_removes_all_rows() {
    assert_regulatory_only_rows(&["--coding_only"], false, false);
}

#[test]
fn regulatory_only_canonical_retains_intergenic_and_regulatory_rows() {
    assert_regulatory_only_rows(&["--canonical"], true, true);
}

#[test]
fn regulatory_only_biotype_retains_intergenic_and_regulatory_rows() {
    assert_regulatory_only_rows(&["--biotype"], true, true);
}

#[test]
fn regulatory_only_no_intergenic_removes_all_rows() {
    assert_regulatory_only_rows(&["--no_intergenic"], false, false);
}

/// A corrupt regulatory shard fails the embedding annotator's construction, as a
/// corrupt transcript shard does, instead of annotating every variant of the
/// chromosome without its regulatory rows and reporting the failure as an event.
#[test]
fn embedding_annotator_refuses_a_corrupt_regulatory_shard_at_construction() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    let config = vep_cli::Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--regulatory",
        "--no_fasta",
    ])
    .unwrap();
    let error = match vep_cli::Annotator::from_config(&config) {
        Err(error) => format!("{error:#}"),
        Ok(_) => panic!("a cache with a corrupt regulatory shard built an annotator"),
    };
    assert!(
        error.contains("parsing") && error.contains("regulatory/22/1-1000000.json"),
        "{error}"
    );
}

#[test]
fn embedding_batches_load_regulatory_features_and_apply_summary_filters() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    // The embedding annotator loads every chromosome when it is built.
    fs::remove_dir_all(cache.join("regulatory/22")).unwrap();
    let config = vep_cli::Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--regulatory",
        "--no_fasta",
    ])
    .unwrap();
    let annotator = vep_cli::Annotator::from_config(&config).unwrap();
    let mut variants = vec![
        vep_cli::InputVariant::new("21".into(), 106, 106, b"C".to_vec(), b"G".to_vec()),
        vep_cli::InputVariant::new("21".into(), 115, 115, b"A".to_vec(), b"G".to_vec()),
    ];
    annotator.annotate_batch(&mut variants);
    assert_eq!(
        variants[0].most_severe_consequence,
        Some(vep_cli::Consequence::TfBindingSiteVariant)
    );
    assert!(variants[0].intergenic_overlap);
    assert!(variants[1]
        .transcript_consequences
        .iter()
        .any(|r| r.motif.as_deref().is_some_and(|m| m.pos == Some(4))));
    let mut summary = variants.clone();
    let mut summary_alone = vep_cli::FilterConfig::default();
    summary_alone.summary = true;
    vep_cli::apply_to_batch(&mut summary, &summary_alone);
    assert_eq!(
        summary[0].summary_terms.as_deref(),
        Some(
            &[
                vep_cli::Consequence::TfBindingSiteVariant,
                vep_cli::Consequence::RegulatoryRegionVariant,
                vep_cli::Consequence::IntergenicVariant,
            ][..]
        )
    );
    let mut coding = variants.clone();
    let mut coding_only = vep_cli::FilterConfig::default();
    coding_only.row_filters.coding_only = true;
    vep_cli::apply_to_batch(&mut coding, &coding_only);
    assert!(coding
        .iter()
        .all(|v| v.transcript_consequences.is_empty() && v.most_severe_consequence.is_none()));
    let mut no_intergenic = vep_cli::FilterConfig::default();
    no_intergenic.row_filters.no_intergenic = true;
    vep_cli::apply_to_batch(&mut variants, &no_intergenic);
    assert!(variants
        .iter()
        .all(|v| v.transcript_consequences.is_empty() && v.most_severe_consequence.is_none()));
}

/// The three GFF fixture files under their keys, as a `--regulatory_gff` value.
fn regulatory_gff_value() -> String {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regulatory-gff");
    format!(
        "file={},emars={},motifs={}",
        fixtures.join("regulatory_features.gff3.gz").display(),
        fixtures.join("EMARs.gff.gz").display(),
        fixtures.join("motif_features.gff3.gz").display()
    )
}

/// A run over `in.vcf` under `dir` with `cache`, in the tab format with the
/// extra arguments, returning the data rows, the warnings file's text and the
/// run's standard error; the run must succeed.
fn tab_rows(dir: &Path, cache: &Path, extra_args: &[&str]) -> (Vec<String>, String, String) {
    let output = dir.join("out.tab");
    let _ = fs::remove_file(&output);
    let result = common::isolated_vep(dir)
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--no_fasta",
            "--force_overwrite",
            "--tab",
        ])
        .arg("--input_file")
        .arg(dir.join("in.vcf"))
        .arg("--output_file")
        .arg(&output)
        .args(extra_args)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "{stderr}");
    let rows = fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    let warnings = fs::read_to_string(dir.join("out.tab_warnings.txt")).unwrap_or_default();
    (rows, warnings, stderr)
}

/// `--regulatory_gff` replaces the cache's regulatory data through the binary:
/// the warning of `CacheDir.pm:161-164` is written once, the cache's enhancer
/// at 100-120 yields no row, and the GFF records yield theirs: an enhancer and
/// the EMAR of the same span as two RegulatoryFeature rows with their types as
/// BIOTYPE, two motifs as MotifFeature rows carrying MOTIF_NAME, MOTIF_POS and
/// TRANSCRIPTION_FACTORS and no HIGH_INF_POS or MOTIF_SCORE_CHANGE, and a base
/// in a promoter's extension a row only under `--extended_promoters`, whose
/// absence of effect without `--regulatory_gff` is byte-identical output. A
/// structural deletion over the enhancer, its EMAR and four motifs lists its
/// rows as the sources' features were met, the `file=` record before the
/// `emars=` record of the same span and the motifs in file order
/// (`AnnotationType/RegFeat.pm:135`, `StructuralVariationFeature.pm:495-509`),
/// each with its overlap columns.
#[test]
fn regulatory_gff_replaces_the_cache_through_the_cli() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t106\tcache_enhancer\tC\tG\t.\t.\t.\n\
         21\t5010240\tsv_del_over_gff\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=5011600\n\
         21\t5010250\tmotif_pair\tA\tC\t.\t.\t.\n\
         21\t5011400\tenhancer_and_emar\tA\tC\t.\t.\t.\n\
         21\t5553160\tpromoter_extension\tA\tC\t.\t.\t.\n",
    )
    .unwrap();
    let gff = regulatory_gff_value();
    let (rows, warnings, _) = tab_rows(dir.path(), &cache, &["--regulatory_gff", &gff]);
    // The fixture cache carries no checksum manifest, which is warned about too.
    let gff_warnings: Vec<&str> = warnings
        .lines()
        .filter(|line| line.contains("regulatory"))
        .collect();
    assert_eq!(
        gff_warnings,
        ["WARNING: Using regulatory features from --regulatory_gff; ignoring the regulatory data in this cache"]
    );
    let feature_rows = |id: &str| -> Vec<&String> {
        rows.iter()
            .filter(|row| row.starts_with(&format!("{id}\t")))
            .collect()
    };
    let cache_rows = feature_rows("cache_enhancer");
    assert_eq!(cache_rows.len(), 1);
    assert!(
        cache_rows[0].contains("\tintergenic_variant\t"),
        "{}",
        cache_rows[0]
    );
    let motif_rows = feature_rows("motif_pair");
    let motif_features: Vec<&str> = motif_rows
        .iter()
        .filter(|row| row.contains("\tMotifFeature\t"))
        .map(|row| row.split('\t').nth(4).unwrap())
        .collect();
    assert_eq!(motif_features, ["ENSM00000574652", "ENSM00000574653"]);
    for row in &motif_rows {
        if row.contains("\tMotifFeature\t") {
            let fields: Vec<&str> = row.split('\t').collect();
            assert!(
                fields.contains(&"ENSPFM0085")
                    && fields.contains(&"EBF1")
                    && fields.contains(&"TF_binding_site_variant"),
                "{row}"
            );
        }
    }
    let header = fs::read_to_string(dir.path().join("out.tab")).unwrap();
    let columns = header
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap()
        .to_string();
    // The columns follow the flags, as Ensembl VEP's do under `--regulatory`,
    // which `--regulatory_gff` switches on; the two scored fields hold no value.
    let column_index = |name: &str| columns.split('\t').position(|c| c == name).unwrap();
    for column in [
        "MOTIF_NAME",
        "MOTIF_POS",
        "HIGH_INF_POS",
        "MOTIF_SCORE_CHANGE",
        "TRANSCRIPTION_FACTORS",
        "BIOTYPE",
    ] {
        assert!(columns.contains(column), "{columns}");
    }
    assert!(!columns.contains("CELL_TYPE"), "{columns}");
    for row in motif_rows
        .iter()
        .filter(|row| row.contains("\tMotifFeature\t"))
    {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[column_index("HIGH_INF_POS")], "-", "{row}");
        assert_eq!(fields[column_index("MOTIF_SCORE_CHANGE")], "-", "{row}");
        assert_eq!(fields[column_index("MOTIF_NAME")], "ENSPFM0085", "{row}");
    }
    let regulatory: Vec<(&str, &str)> = feature_rows("enhancer_and_emar")
        .iter()
        .filter(|row| row.contains("\tRegulatoryFeature\t"))
        .map(|row| {
            let fields: Vec<&str> = row.split('\t').collect();
            (fields[4], fields[column_index("BIOTYPE")])
        })
        .collect();
    assert_eq!(
        regulatory,
        [("ENSR00001793149", "EMAR"), ("ENSR21_85ZF5", "enhancer")]
    );
    // The structural deletion's rows in the default format: feature, type, terms
    // and the Extra column, Ensembl VEP 116.3's own rows and order over these
    // files (the `file=` record before the `emars=` record of the same span,
    // then the motifs in file order).
    let default = annotate_with_args(
        dir.path(),
        &cache,
        "default",
        false,
        &["--regulatory_gff", &gff],
    );
    let structural: Vec<(&str, &str, &str, &str)> = default
        .lines()
        .filter(|line| {
            line.starts_with("sv_del_over_gff\t") && !line.contains("\tintergenic_variant\t")
        })
        .map(|row| {
            let fields: Vec<&str> = row.split('\t').collect();
            (fields[4], fields[5], fields[6], fields[13])
        })
        .collect();
    let regulatory_terms = "regulatory_region_ablation,regulatory_region_variant";
    let motif_terms = "TFBS_ablation,TF_binding_site_variant";
    assert_eq!(
        structural,
        [
            (
                "ENSR21_85ZF5",
                "RegulatoryFeature",
                regulatory_terms,
                "IMPACT=MODIFIER;OverlapBP=231;OverlapPC=100.00"
            ),
            (
                "ENSR00001793149",
                "RegulatoryFeature",
                regulatory_terms,
                "IMPACT=MODIFIER;OverlapBP=231;OverlapPC=100.00"
            ),
            (
                "ENSM00000574652",
                "MotifFeature",
                motif_terms,
                "IMPACT=MODERATE;OverlapBP=14;OverlapPC=100.00"
            ),
            (
                "ENSM00000574653",
                "MotifFeature",
                motif_terms,
                "IMPACT=MODERATE;OverlapBP=14;OverlapPC=100.00"
            ),
            (
                "ENSM00000574654",
                "MotifFeature",
                motif_terms,
                "IMPACT=MODERATE;OverlapBP=14;OverlapPC=100.00"
            ),
            (
                "ENSM00000574655",
                "MotifFeature",
                motif_terms,
                "IMPACT=MODERATE;OverlapBP=14;OverlapPC=100.00"
            ),
        ]
    );
    let structural_tab: Vec<&str> = feature_rows("sv_del_over_gff")
        .iter()
        .filter(|row| !row.contains("\tintergenic_variant\t"))
        .map(|row| row.split('\t').nth(4).unwrap())
        .collect();
    assert_eq!(
        structural_tab,
        structural.iter().map(|row| row.0).collect::<Vec<_>>(),
        "the tab rows come in the same order"
    );
    // The EMAR starting at the promoter's extended start covers the base; the
    // promoter itself does not without the flag.
    let extension_rows = feature_rows("promoter_extension");
    assert!(
        extension_rows
            .iter()
            .any(|row| row.contains("\tENSR00001793220\t")),
        "{extension_rows:?}"
    );
    assert!(
        extension_rows
            .iter()
            .all(|row| !row.contains("\tENSR21_B34QG\t")),
        "{extension_rows:?}"
    );

    let (extended, _, _) = tab_rows(
        dir.path(),
        &cache,
        &["--regulatory_gff", &gff, "--extended_promoters"],
    );
    // The EMAR row stays beside the promoter's.
    let promoter: Vec<&String> = extended
        .iter()
        .filter(|row| row.starts_with("promoter_extension\t") && row.contains("\tENSR21_B34QG\t"))
        .collect();
    assert_eq!(promoter.len(), 1, "{extended:?}");
    assert!(
        promoter[0].contains("\tRegulatoryFeature\t") && promoter[0].contains("\tpromoter"),
        "{}",
        promoter[0]
    );

    // `--cell_type` with the GFF files: no name is checked against the cache
    // (HUVEC is in no list), the source's warning is written once, and every
    // regulatory and motif row carries CELL_TYPE as the empty list, which
    // Ensembl VEP sets on each such row under the flag (OutputFactory.pm:1888,
    // 1936): an empty tab cell where every other row prints `-` (Tab.pm:128),
    // `CELL_TYPE=` in the default format's Extra column and `[]` in JSON.
    let is_feature_row = |line: &str| {
        matches!(
            line.split('\t').nth(5),
            Some("RegulatoryFeature" | "MotifFeature")
        )
    };
    let (with_cell_type, warnings, _) = tab_rows(
        dir.path(),
        &cache,
        &["--regulatory_gff", &gff, "--cell_type", "HUVEC,HUVEC"],
    );
    assert_eq!(
        warnings
            .lines()
            .filter(|line| line.contains("--cell_type"))
            .collect::<Vec<_>>(),
        ["WARNING: --cell_type is ignored for regulatory features from --regulatory_gff; a GFF carries no epigenome activity data"]
    );
    let header = fs::read_to_string(dir.path().join("out.tab")).unwrap();
    let header = header
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap()
        .to_string();
    let cell_type = header.split('\t').position(|c| c == "CELL_TYPE").unwrap();
    assert!(
        with_cell_type.iter().all(|row| {
            let expected = if is_feature_row(row) { "" } else { "-" };
            row.split('\t').nth(cell_type) == Some(expected)
        }),
        "{with_cell_type:?}"
    );
    assert_eq!(
        with_cell_type.len(),
        rows.len(),
        "the same rows as without the flag"
    );
    let default = annotate_with_args(
        dir.path(),
        &cache,
        "default",
        false,
        &["--regulatory_gff", &gff, "--cell_type", "HUVEC"],
    );
    let feature_rows = rows.iter().filter(|row| is_feature_row(row)).count();
    assert!(feature_rows > 0);
    let extras: Vec<&str> = default
        .lines()
        .filter(|line| !line.starts_with('#') && is_feature_row(line))
        .map(|line| line.split('\t').nth(13).unwrap())
        .collect();
    assert_eq!(extras.len(), feature_rows, "{default}");
    assert!(
        extras
            .iter()
            .all(|extra| extra.split(';').any(|field| field == "CELL_TYPE=")),
        "{default}"
    );
    let json = annotate_with_args(
        dir.path(),
        &cache,
        "json",
        false,
        &["--regulatory_gff", &gff, "--cell_type", "HUVEC"],
    );
    let lists: Vec<Value> = json
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .flat_map(|record| {
            [
                "regulatory_feature_consequences",
                "motif_feature_consequences",
            ]
            .into_iter()
            .filter_map(|key| record[key].as_array().cloned())
            .flatten()
            .map(|row| row["cell_type"].clone())
            .collect::<Vec<_>>()
        })
        .collect();
    assert_eq!(lists.len(), feature_rows, "{json}");
    assert!(
        lists.iter().all(|list| *list == serde_json::json!([])),
        "{json}"
    );

    // Without --regulatory_gff the flag reaches nothing: the output bytes are
    // the same, the cache's regulatory data read the same way.
    let (with_flag, _, _) = tab_rows(
        dir.path(),
        &cache,
        &["--regulatory", "--extended_promoters"],
    );
    let with_flag_text = fs::read_to_string(dir.path().join("out.tab")).unwrap();
    let (without_flag, _, _) = tab_rows(dir.path(), &cache, &["--regulatory"]);
    let without_flag_text = fs::read_to_string(dir.path().join("out.tab")).unwrap();
    assert_eq!(with_flag, without_flag);
    // The header's time and command line are the run's own; every other line
    // and every byte of the rows are the same.
    fn steady(text: &str) -> Vec<&str> {
        text.lines()
            .filter(|line| {
                !line.starts_with("## Output produced at ")
                    && !line.starts_with("## VEP command-line: ")
            })
            .collect()
    }
    assert_eq!(steady(&with_flag_text), steady(&without_flag_text));
    assert!(with_flag.iter().any(|row| row.contains("\tENSR_FIXTURE\t")));
}

/// A `--regulatory_gff` value the grammar refuses, a file that does not exist
/// and a file without its tabix index each exit 1 with the contract's line
/// under the `cli` component.
#[test]
fn regulatory_gff_refusals_exit_with_their_message() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tr\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let unindexed = dir.path().join("features.gff3.gz");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regulatory-gff/EMARs.gff.gz"),
        &unindexed,
    )
    .unwrap();
    let cache = cache(dir.path());
    let missing = dir.path().join("missing.gff3.gz");
    for (value, message) in [
        (
            "file=a.gff3.gz,motifs".to_string(),
            "ERROR: [cli] Failed to parse --regulatory_gff parameter \"motifs\"; expected key=value".to_string(),
        ),
        (
            "matrices=m.json".to_string(),
            "ERROR: [cli] Unsupported --regulatory_gff key(s): matrices\nSupported keys: emars, file, motifs".to_string(),
        ),
        (
            "file=0".to_string(),
            "ERROR: [cli] --regulatory_gff requires at least one of: emars, file, motifs".to_string(),
        ),
        (
            format!("emars={}", missing.display()),
            format!("ERROR: [cli] --regulatory_gff emars={} not found", missing.display()),
        ),
        (
            unindexed.display().to_string(),
            format!("ERROR: [cli] --regulatory_gff file={} has no tabix index", unindexed.display()),
        ),
    ] {
        let output = dir.path().join("out.txt");
        let result = common::isolated_vep(dir.path())
            .args(["--json_cache", cache.to_str().unwrap(), "--assembly", "GRCh38", "--no_fasta", "--no_stats", "--force_overwrite"])
            .args(["--regulatory_gff", &value])
            .arg("--input_file")
            .arg(dir.path().join("in.vcf"))
            .arg("--output_file")
            .arg(&output)
            .arg("--quiet")
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1), "{value}");
        assert_eq!(String::from_utf8_lossy(&result.stderr).trim_end(), message, "{value}");
    }
}

/// A refused regulatory source leaves no file behind. Ensembl VEP 116.3 builds
/// every annotation source in `init` (Runner.pm:128): a `--regulatory_gff` file
/// that is not there throws at AnnotationSourceAdaptor.pm:266, a `--cell_type`
/// the cache lacks at AnnotationType/RegFeat.pm:87-93 (from Cache/RegFeat.pm:117),
/// each before the statistics report opens (Runner.pm:144-148) and before `run`
/// opens the output and prints the headers (Runner.pm:202-205); a report path
/// that cannot be written fails there (BaseRunner.pm:280), still before the
/// output. Each refusal exits 1 with no output file and no report; the same run
/// with its sources present writes the header lines, then its rows, and the
/// report.
#[test]
fn refused_regulatory_source_leaves_no_output_file_or_report() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tr\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let unindexed = dir.path().join("features.gff3.gz");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regulatory-gff/EMARs.gff.gz"),
        &unindexed,
    )
    .unwrap();
    let cache = cache(dir.path());
    let missing = dir.path().join("missing.gff3.gz");
    let output = dir.path().join("out.txt");
    let report = dir.path().join("out.txt_summary.html");
    let run = |args: &[&str]| {
        common::isolated_vep(dir.path())
            .args([
                "--json_cache",
                cache.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--no_fasta",
                "--force_overwrite",
                "--quiet",
            ])
            .args(args)
            .arg("--input_file")
            .arg(dir.path().join("in.vcf"))
            .arg("--output_file")
            .arg(&output)
            .output()
            .unwrap()
    };
    let missing_value = format!("file={}", missing.display());
    let unwritable_report = dir.path().join("no-such-dir/report.html");
    for (args, message) in [
        (
            ["--stats_file", unwritable_report.to_str().unwrap()],
            format!(
                "ERROR: [output] Could not write to stats file {}",
                unwritable_report.display()
            ),
        ),
        (
            ["--regulatory_gff", missing_value.as_str()],
            format!(
                "ERROR: [cli] --regulatory_gff file={} not found",
                missing.display()
            ),
        ),
        (
            ["--regulatory_gff", unindexed.to_str().unwrap()],
            format!(
                "ERROR: [cli] --regulatory_gff file={} has no tabix index",
                unindexed.display()
            ),
        ),
        (
            ["--cell_type", "bogus"],
            "ERROR: [cli] Cell type \"bogus\" unavailable; available cell types are:".to_string(),
        ),
    ] {
        let result = run(&args);
        assert_eq!(result.status.code(), Some(1), "{args:?}");
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(stderr.starts_with(&message), "{args:?}: {stderr}");
        assert!(!output.exists(), "{args:?} left {}", output.display());
        assert!(!report.exists(), "{args:?} left {}", report.display());
    }

    let result = run(&["--regulatory"]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let first_row = lines
        .iter()
        .position(|line| !line.starts_with('#'))
        .unwrap();
    assert!(
        text.starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR"),
        "{text}"
    );
    assert_eq!(
        lines[first_row - 1].split('\t').next(),
        Some("#Uploaded_variation"),
        "{text}"
    );
    assert!(lines[first_row].contains("\tRegulatoryFeature\t"), "{text}");
    assert!(fs::metadata(&report).unwrap().len() > 0);
}

/// A `--regulatory_gff` record the reader refuses leaves no output or report,
/// whichever chromosome holds it: every record of the files meets the reader's
/// refusals before the statistics report and the output are created. The file
/// holds a sound enhancer on chromosome 21, which the input names, and on
/// chromosome 22, which it does not, a record without its attributes column,
/// which the index admits. The run exits 1 with the reader's refusal on
/// standard error, no output file and no report; the same run without that
/// record writes the enhancer's row and the report.
#[test]
fn a_refused_regulatory_gff_record_leaves_no_output_file_or_report() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tr\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let cache = cache(dir.path());
    let output = dir.path().join("out.txt");
    let report = dir.path().join("out.txt_summary.html");
    let run = |gff: &Path| {
        common::isolated_vep(dir.path())
            .args([
                "--json_cache",
                cache.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--no_fasta",
                "--force_overwrite",
                "--quiet",
            ])
            .arg("--regulatory_gff")
            .arg(format!("file={}", gff.display()))
            .arg("--input_file")
            .arg(dir.path().join("in.vcf"))
            .arg("--output_file")
            .arg(&output)
            .output()
            .unwrap()
    };
    let sound = "21\tEnsembl\tenhancer\t100\t120\t.\t.\t.\tID=ENSR_SOUND";
    let refused = "22\tEnsembl\tenhancer\t100\t120\t.\t.\t.";

    let with_refused = common::indexed_gff(dir.path(), "refused.gff3.gz", &[sound, refused]);
    let result = run(&with_refused);
    assert_eq!(result.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&result.stderr).trim_end(),
        format!(
            "ERROR: [cache] {}: a record of chromosome 22 has no attributes column: {refused}",
            with_refused.display()
        )
    );
    assert!(!output.exists(), "the refusal left {}", output.display());
    assert!(!report.exists(), "the refusal left {}", report.display());

    let sound_alone = common::indexed_gff(dir.path(), "sound.gff3.gz", &[sound]);
    let result = run(&sound_alone);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = fs::read_to_string(&output).unwrap();
    let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert!(
        rows.iter()
            .any(|row| row.contains("\tENSR_SOUND\tRegulatoryFeature\t")),
        "{text}"
    );
    assert!(fs::metadata(&report).unwrap().len() > 0);
}

/// A type no regulatory file yields draws its warning when its chromosome is
/// first read, and the check every record meets before the output is created
/// draws none: a run whose file holds such a line on chromosome 22, which the
/// input does not name, leaves the warnings file of the run without the line.
/// The same file under an input record on chromosome 22 draws the warning.
#[test]
fn a_regulatory_gff_type_on_an_unread_chromosome_draws_no_warning() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    let enhancer = "21\tEnsembl\tenhancer\t100\t120\t.\t.\t.\tID=ENSR_SOUND";
    let exon = "22\tEnsembl\texon\t100\t120\t.\t.\t.\tID=EXON";
    let gff = dir.path().join("features.gff3.gz");
    let warnings_of = |records: &[&str], chromosome: &str| -> String {
        fs::write(
            dir.path().join("in.vcf"),
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 {chromosome}\t106\tr\tC\tG\t.\t.\t.\n"
            ),
        )
        .unwrap();
        common::indexed_gff(dir.path(), "features.gff3.gz", records);
        let _ = fs::remove_file(dir.path().join("out.tab_warnings.txt"));
        let value = format!("file={}", gff.display());
        let (_, warnings, _) = tab_rows(dir.path(), &cache, &["--regulatory_gff", &value]);
        warnings
    };
    let type_warning = format!(
        "WARNING: Ignoring 'exon' feature type from {}; not a supported regulatory feature type",
        gff.display()
    );

    let with_line = warnings_of(&[enhancer, exon], "21");
    let without_line = warnings_of(&[enhancer], "21");
    assert_eq!(with_line, without_line);
    assert!(!with_line.contains(&type_warning), "{with_line}");

    let read = warnings_of(&[enhancer, exon], "22");
    assert_eq!(
        read.lines().filter(|line| *line == type_warning).count(),
        1,
        "{read}"
    );
}
