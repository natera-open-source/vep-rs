// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! End-to-end integration tests for the VEP CLI binary.
//!
//! These tests invoke the compiled `vep` binary as a subprocess and verify
//! correct behaviour for help/version flags, output formats, and annotation.

use std::process::Command;

fn vep_binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vep"))
}

fn test_input_vcf() -> String {
    format!(
        "{}/../../tests/sv_validation/grch37/01_snv.vcf.gz",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// A JSON transcript cache named by `VEP_RS_TEST_CACHE`; `None` when the
/// variable is unset or the directory is absent, and the caller skips.
fn json_cache_path() -> Option<String> {
    let path = std::env::var("VEP_RS_TEST_CACHE").ok()?;
    if std::path::Path::new(&path).is_dir() {
        Some(path)
    } else {
        None
    }
}

#[test]
fn test_help_flag() {
    let output = vep_binary().arg("--help").output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("Predict functional effects of genomic variants"),
        "Expected help text to contain 'Predict functional effects of genomic variants', got:\n{}",
        stdout
    );
}

#[test]
fn test_version_flag() {
    let output = vep_binary().arg("--version").output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("vep"),
        "Expected version output to contain 'vep', got:\n{}",
        stdout
    );
}

#[test]
fn test_end_to_end_with_cache() {
    let json_cache = match json_cache_path() {
        Some(p) => p,
        None => {
            eprintln!("Skipping: VEP_RS_TEST_CACHE not set or directory absent");
            return;
        }
    };

    let input = test_input_vcf();
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("output.txt");

    let result = vep_binary()
        .args(["-i", &input])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--json_cache", &json_cache])
        .args(["--offline", "--force", "--quiet"])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "VEP failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let output = std::fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = output.lines().collect();

    assert!(lines.len() > 5, "Output too short: {} lines", lines.len());

    assert!(
        lines[0].starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR"),
        "Expected header line, got: {}",
        lines[0]
    );

    let data_lines: Vec<&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .copied()
        .collect();
    assert!(!data_lines.is_empty(), "Expected data lines in output");

    let has_coding = data_lines
        .iter()
        .any(|l| l.contains("missense_variant") || l.contains("synonymous_variant"));
    assert!(has_coding, "Expected coding variants in output");

    for line in &data_lines {
        if line.contains("missense_variant") {
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(
                fields.len() >= 13,
                "Expected at least 13 tab-separated fields"
            );
            assert_ne!(fields[10], "-", "Missense should have amino acid change");
            assert_ne!(fields[11], "-", "Missense should have codon change");
            break;
        }
    }
}

#[test]
fn test_json_output_format() {
    let json_cache = match json_cache_path() {
        Some(p) => p,
        None => {
            eprintln!("Skipping: VEP_RS_TEST_CACHE not set or directory absent");
            return;
        }
    };

    let input = test_input_vcf();
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("output.json");

    let result = vep_binary()
        .args(["-i", &input])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--json_cache", &json_cache])
        .args(["--offline", "--force", "--quiet", "--json"])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "VEP JSON output failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let output = std::fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = output.lines().filter(|l| !l.is_empty()).collect();
    assert!(!lines.is_empty(), "Expected JSON output lines");

    for (i, line) in lines.iter().enumerate() {
        let parsed: Result<serde_json::Value, _> = serde_json::from_str(line);
        assert!(
            parsed.is_ok(),
            "Line {} is not valid JSON: {}\nParse error: {}",
            i + 1,
            line,
            parsed.unwrap_err()
        );
    }
}

#[test]
fn test_vcf_output_format() {
    let json_cache = match json_cache_path() {
        Some(p) => p,
        None => {
            eprintln!("Skipping: VEP_RS_TEST_CACHE not set or directory absent");
            return;
        }
    };

    let input = test_input_vcf();
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("output.vcf");

    let result = vep_binary()
        .args(["-i", &input])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--json_cache", &json_cache])
        .args(["--offline", "--force", "--quiet", "--vcf"])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "VEP VCF output failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let output = std::fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = output.lines().collect();

    let has_vcf_header = lines.iter().any(|l| l.starts_with("##fileformat=VCF"));
    assert!(has_vcf_header, "Expected VCF fileformat header in output");

    let has_csq_header = lines.iter().any(|l| l.contains("ID=CSQ"));
    assert!(has_csq_header, "Expected CSQ INFO header in VCF output");

    let data_lines: Vec<&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .copied()
        .collect();
    assert!(!data_lines.is_empty(), "Expected data lines in VCF output");

    let has_csq_data = data_lines.iter().any(|l| l.contains("CSQ="));
    assert!(has_csq_data, "Expected CSQ annotation in VCF data lines");
}

#[test]
fn test_pick_reduces_output() {
    let json_cache = match json_cache_path() {
        Some(p) => p,
        None => {
            eprintln!("Skipping: VEP_RS_TEST_CACHE not set or directory absent");
            return;
        }
    };

    let input = test_input_vcf();

    let output_dir_no_pick = tempfile::tempdir().unwrap();
    let output_no_pick = output_dir_no_pick.path().join("no_pick.txt");

    let result_no_pick = vep_binary()
        .args(["-i", &input])
        .args(["-o", output_no_pick.to_str().unwrap()])
        .args(["--json_cache", &json_cache])
        .args(["--offline", "--force", "--quiet"])
        .output()
        .unwrap();

    assert!(
        result_no_pick.status.success(),
        "VEP (no pick) failed: {}",
        String::from_utf8_lossy(&result_no_pick.stderr)
    );

    let output_dir_pick = tempfile::tempdir().unwrap();
    let output_pick = output_dir_pick.path().join("pick.txt");

    let result_pick = vep_binary()
        .args(["-i", &input])
        .args(["-o", output_pick.to_str().unwrap()])
        .args(["--json_cache", &json_cache])
        .args(["--offline", "--force", "--quiet", "--pick"])
        .output()
        .unwrap();

    assert!(
        result_pick.status.success(),
        "VEP (--pick) failed: {}",
        String::from_utf8_lossy(&result_pick.stderr)
    );

    let no_pick_output = std::fs::read_to_string(&output_no_pick).unwrap();
    let pick_output = std::fs::read_to_string(&output_pick).unwrap();

    let no_pick_data: Vec<&str> = no_pick_output
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    let pick_data: Vec<&str> = pick_output
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();

    assert!(
        pick_data.len() <= no_pick_data.len(),
        "--pick produced {} lines but no-pick produced {} lines",
        pick_data.len(),
        no_pick_data.len()
    );

    assert!(!pick_data.is_empty(), "--pick produced no output lines");
}

#[test]
fn test_missing_input_file() {
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("output.txt");

    let result = vep_binary()
        .args(["-i", "/nonexistent/file.vcf"])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--offline", "--force", "--quiet"])
        .output()
        .unwrap();

    assert!(
        !result.status.success(),
        "Expected VEP to fail with nonexistent input file"
    );
}

#[test]
fn test_output_file_exists_without_force() {
    let json_cache = match json_cache_path() {
        Some(p) => p,
        None => {
            eprintln!("Skipping: VEP_RS_TEST_CACHE not set or directory absent");
            return;
        }
    };

    let input = test_input_vcf();
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("existing.txt");

    std::fs::write(&output_path, "existing content").unwrap();

    let result = vep_binary()
        .args(["-i", &input])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--json_cache", &json_cache])
        .args(["--offline", "--quiet"])
        .output()
        .unwrap();

    assert!(
        !result.status.success(),
        "Expected VEP to fail when output file exists without --force"
    );
}

/// The three flags vep-rs declares but refuses exit with a usage error naming the
/// feature they belong to, before any input is read.
#[test]
fn refused_flags_exit_with_their_message() {
    for (flag, feature) in [
        (
            vec!["--regulatory_gff", "reg.gff3.gz"],
            "regulatory annotation",
        ),
        (vec!["--extended_promoters"], "regulatory annotation"),
        (vec!["--custom_suppress_filter"], "custom annotation"),
    ] {
        let result = vep_binary()
            .args(&flag)
            .args(["-i", "/nonexistent/file.vcf", "--quiet"])
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(2), "{flag:?}");
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(
            stderr.contains(&format!("{} is not supported", flag[0])) && stderr.contains(feature),
            "{flag:?}: {stderr}"
        );
    }
}

/// Annotates one VCF record in the given format against the release 116 golden
/// corpus's own pruned cache, which unlike `VEP_RS_TEST_CACHE` is always
/// present, and returns the output's data lines.
fn annotate_record(record: &str, format_flags: &[&str], extra: &[&str]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    let output = dir.path().join("out");
    let cache = format!(
        "{}/../../tests/golden/116/GRCh38-release116/json_cache",
        env!("CARGO_MANIFEST_DIR")
    );
    let result = vep_binary()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args(["--json_cache", &cache])
        .args(["--offline", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(format_flags)
        .args(extra)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{format_flags:?} {extra:?}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The Location, Allele and Consequence columns of the default-format rows that
/// contain `text`, a transcript ID or the record's own ID.
fn rows_on(rows: &[String], text: &str) -> Vec<(String, String, String)> {
    rows.iter()
        .filter(|row| row.contains(text))
        .map(|row| {
            let c: Vec<&str> = row.split('\t').collect();
            (c[1].to_string(), c[2].to_string(), c[6].to_string())
        })
        .collect()
}

/// A [`rows_on`] expectation.
fn row(location: &str, allele: &str, consequence: &str) -> (String, String, String) {
    (
        location.to_string(),
        allele.to_string(),
        consequence.to_string(),
    )
}

/// A structural variant wider than the default `--max_sv_size` (`del_10_6_mb`
/// of `tests/golden/116/GRCh38-release116`) is annotated, not dropped: the
/// default output lists its transcript rows, the VCF output carries its line
/// without `CSQ`, and the JSON output omits it. Ensembl VEP 116.2 writes nothing
/// for the record in any format (`Parser.pm` `validate_vf` returns 0 at the
/// size test); this is the documented divergence of `--max_sv_size`. With `-1`
/// the limit is lifted and every format annotates the record, as Ensembl VEP
/// 116.2 does under `-1`.
#[test]
fn oversize_structural_variant_is_annotated_not_dropped() {
    let record = "21\t17750000\tdel_10_6_mb\tN\t<DEL>\t.\t.\tEND=28350000;SVTYPE=DEL";

    let default_rows = annotate_record(record, &[], &[]);
    assert!(
        default_rows.len() > 1 && default_rows.iter().all(|r| r.starts_with("del_10_6_mb\t")),
        "default output rows: {}",
        default_rows.len()
    );
    let vcf_lines = annotate_record(record, &["--vcf"], &[]);
    assert_eq!(vcf_lines.len(), 1, "{vcf_lines:?}");
    assert!(
        vcf_lines[0].starts_with("21\t17750000\tdel_10_6_mb\tN\t<DEL>\t")
            && !vcf_lines[0].contains("CSQ="),
        "{}",
        vcf_lines[0]
    );
    let json_lines = annotate_record(record, &["--json"], &[]);
    assert!(json_lines.is_empty(), "{json_lines:?}");

    let lifted = ["--max_sv_size", "-1"];
    let vcf_lifted = annotate_record(record, &["--vcf"], &lifted);
    assert!(
        vcf_lifted.len() == 1 && vcf_lifted[0].contains("CSQ=deletion|"),
        "{vcf_lifted:?}"
    );
    let json_lifted = annotate_record(record, &["--json"], &lifted);
    assert!(
        json_lifted.len() == 1 && json_lifted[0].contains("\"transcript_consequences\""),
        "{json_lifted:?}"
    );
}

/// A record whose only ALT is a single breakend prints `chromosome_breakpoint` and
/// takes no `feature_truncation` whether or not INFO carries MATEID: Ensembl VEP
/// 116.2 never reads MATEID, `_parse_breakends` finds no bracket in `.N`
/// (`StructuralVariationFeature.pm:1279-1283`), and the row's Allele falls to
/// `class_SO_term` (`OutputFactory.pm:2029`).
#[test]
fn single_breakend_with_mateid_prints_the_class_term() {
    let without = annotate_record("21\t30425650\tsb\tN\t.N\t.\t.\tSVTYPE=BND", &[], &[]);
    let with = annotate_record(
        "21\t30425650\tsb\tN\t.N\t.\t.\tSVTYPE=BND;MATEID=bnd_x",
        &[],
        &[],
    );
    assert_eq!(with, without);
    assert!(!without.is_empty());
    for row in &without {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns[2], "chromosome_breakpoint", "{row}");
        assert!(!columns[6].contains("feature_truncation"), "{row}");
    }
}

/// A breakend record listing two bracket ALTs writes three rows, the record's own
/// breakend `N.` and one per mate: Ensembl VEP 116.2 annotates `($vf, @$breakends)`
/// (`StructuralVariationOverlap.pm:90`), one `$vf` entry per record however many
/// brackets the ALT column lists. Its JSON `allele_string` is the whole record, REF
/// and both ALTs joined (`Parser/VCF.pm:469`, `:499`; `OutputFactory/JSON.pm:241`).
#[test]
fn two_bracket_breakend_alts_write_three_rows() {
    let record = "21\t11000000\tmb2\tN\tN[21:11500000[,N]21:11600000]\t.\t.\tSVTYPE=BND";
    let rows = annotate_record(record, &[], &[]);
    let alleles: Vec<&str> = rows.iter().map(|r| r.split('\t').nth(2).unwrap()).collect();
    assert_eq!(
        alleles,
        ["N.", "N[21:11500000[", "N]21:11600000]"],
        "{rows:?}"
    );
    let json = annotate_record(record, &["--json"], &[]);
    assert_eq!(json.len(), 1, "{json:?}");
    assert!(
        json[0].contains("\"allele_string\":\"N/N[21:11500000[/N]21:11600000]\""),
        "{}",
        json[0]
    );
}

/// A bare `<BND>`, with neither CHR2 nor a span, prints `chromosome_breakpoint`:
/// Ensembl VEP 116.2 keeps the ALT as written without CHR2 (the rewrite at
/// `Parser/VCF.pm:491-498` is skipped), `_parse_breakends` finds no bracket, and the
/// Allele falls to `class_SO_term` (`OutputFactory.pm:2029`), with or without
/// INFO/MATEID.
#[test]
fn bare_symbolic_bnd_prints_the_class_term() {
    for info in ["SVTYPE=BND", "SVTYPE=BND;MATEID=bnd_x"] {
        let rows = annotate_record(
            &format!("21\t30425650\tbare\tN\t<BND>\t.\t.\t{info}"),
            &[],
            &[],
        );
        assert!(!rows.is_empty(), "{info}");
        for row in &rows {
            assert_eq!(
                row.split('\t').nth(2),
                Some("chromosome_breakpoint"),
                "{info}: {row}"
            );
        }
    }
}

/// `--max_sv_size` governs a `<CNV:TR>` record's expansion as it does Ensembl VEP
/// 116.2's (`Parser/VCF.pm:403`, `:417-418`, `:428`): under the default a 200-base
/// allele against a 300-base run is a 100-base contraction, `frameshift_variant`; under
/// `-1`, and under a limit the 200-base allele exceeds, the record keeps the symbolic
/// reading, a gain over the run, `feature_elongation`.
#[test]
fn max_sv_size_governs_tandem_repeat_expansion() {
    let record = "21\t30425600\ttr\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425899;SVLEN=300;RUS=AT;RUC=100;RB=300";
    let consequences = |extra: &[&str]| -> String {
        annotate_record(record, &[], extra)
            .iter()
            .find(|row| row.contains("ENST00000390690"))
            .map(|row| row.split('\t').nth(6).unwrap().to_string())
            .expect("the run overlaps ENST00000390690")
    };
    assert_eq!(consequences(&[]), "frameshift_variant,feature_truncation");
    for limit in ["-1", "100"] {
        assert_eq!(
            consequences(&["--max_sv_size", limit]),
            "feature_elongation,coding_sequence_variant",
            "--max_sv_size {limit}"
        );
    }
}

/// The rows Ensembl VEP 116.2 writes on ENST00000390690 (strand -1, 30425380-30425926) and
/// ENST00000334068 (4.5 kb upstream on +) for the record shapes the parser routes and
/// classes as Perl does: an `N` reference base with `SVTYPE=DEL;END` is a deletion over
/// END and alone a `coding_sequence_variant` SNV (no peptide from an ambiguous reference);
/// `<*>` with `SVTYPE=DEL;END` is the same deletion; a lowercase `<ins:me:alu>` takes
/// `SVTYPE=INS` as an `insertion` and `<del>` with `SVTYPE=DUP` is a `duplication`; a
/// gVCF site `N G,<NON_REF>` is its `G` allele; `<BND>,<DEL>` with `SVTYPE=BND` puts both
/// alleles at the point; a `<CNV:TR>` gain crossing the transcript's edge reads the
/// region terms over its run; and the `<INS>` pair's DISTANCE is measured from its nearer
/// flank (4569, not 4570).
#[test]
fn routing_class_and_distance_follow_ensembl_vep() {
    let rows = |record: &str| annotate_record(record, &[], &[]);
    let deletion = row(
        "21:30425651-30425949",
        "deletion",
        "feature_truncation,coding_sequence_variant,5_prime_UTR_variant",
    );
    for record in [
        "21\t30425650\tnref\tN\tG\t.\t.\tSVTYPE=DEL;END=30425949",
        "21\t30425650\tstar\tN\t<*>\t.\t.\tSVTYPE=DEL;END=30425949",
        "21\t30425650\tdotalt\tN\t.N\t.\t.\tSVTYPE=DEL;END=30425949",
    ] {
        assert_eq!(
            rows_on(&rows(record), "ENST00000390690"),
            std::slice::from_ref(&deletion),
            "{record}"
        );
    }
    assert_eq!(
        rows_on(
            &rows("21\t30425650\tnref_snv\tN\tG\t.\t.\t."),
            "ENST00000390690"
        ),
        [row("21:30425650", "G", "coding_sequence_variant")]
    );
    assert_eq!(
        rows_on(
            &rows("21\t30425650\tgvcf\tN\tG,<NON_REF>\t.\t.\t."),
            "ENST00000390690"
        ),
        [row("21:30425650", "G", "coding_sequence_variant")]
    );
    assert_eq!(
        rows_on(
            &rows("21\t30425650\tlower_ins\tN\t<ins:me:alu>\t.\t.\tSVTYPE=INS"),
            "ENST00000390690"
        ),
        [row(
            "21:30425651",
            "insertion",
            "feature_elongation,coding_sequence_variant"
        )]
    );
    assert_eq!(
        rows_on(
            &rows("21\t30425650\tlower_del\tN\t<del>\t.\t.\tSVTYPE=DUP;END=30425949"),
            "ENST00000390690"
        ),
        [row(
            "21:30425651-30425949",
            "duplication",
            "coding_sequence_variant,5_prime_UTR_variant"
        )]
    );
    let breakpoint = row(
        "21:30425651",
        "chromosome_breakpoint",
        "coding_sequence_variant",
    );
    assert_eq!(
        rows_on(
            &rows("21\t30425650\tbnd_del\tN\t<BND>,<DEL>\t.\t.\tSVTYPE=BND;END=30425949"),
            "ENST00000390690"
        ),
        [breakpoint.clone(), breakpoint]
    );
    assert_eq!(
        rows_on(
            &rows("21\t30425650\ttr_gain\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425949;SVLEN=300;RUS=AT;RUC=200"),
            "ENST00000390690"
        ),
        [row("21:30425651-30425950", "tandem_repeat", "coding_sequence_variant,5_prime_UTR_variant")]
    );
    let pair = rows("21\t30425650\tins_pair\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=0");
    let upstream = pair
        .iter()
        .find(|r| r.contains("ENST00000334068"))
        .expect("the upstream row");
    assert!(upstream.contains("DISTANCE=4569"), "{upstream}");
}

/// A coordinate past the parser's ceiling (2^62 - 1): an END of `i64::MAX` or
/// `u64::MAX` fails its record, which is skipped as Ensembl VEP skips an oversize or
/// invalid one, and a span from SVLEN saturates at the ceiling, so SVLEN at either
/// extreme annotates a span to it in every output format without a panic (the
/// effects layer adds one to the end in `i64` on a partial transcript overlap).
#[test]
fn coordinates_past_the_ceiling_are_annotated_without_a_panic() {
    for end in ["END=9223372036854775807", "END=18446744073709551615"] {
        let record = format!("21\t30425600\tdel\tN\t<DEL>\t.\t.\tSVTYPE=DEL;{end}");
        for format in [&[][..], &["--json"], &["--vcf"], &["--tab"]] {
            let rows = annotate_record(&record, format, &["--max_sv_size", "-1"]);
            assert!(rows.is_empty(), "{end} {format:?}: {rows:?}");
        }
    }
    for svlen in ["SVLEN=-9223372036854775808", "SVLEN=9223372036854775807"] {
        let record = format!("21\t30425600\tdel\tN\t<DEL>\t.\t.\tSVTYPE=DEL;{svlen}");
        for format in [&[][..], &["--json"], &["--vcf"], &["--tab"]] {
            let rows = annotate_record(&record, format, &["--max_sv_size", "-1"]);
            assert!(!rows.is_empty(), "{svlen} {format:?}");
        }
        let location = annotate_record(&record, &[], &["--max_sv_size", "-1"])
            .iter()
            .find(|row| row.contains("ENST00000390690"))
            .map(|row| row.split('\t').nth(1).unwrap().to_string())
            .expect("ENST00000390690 overlaps POS + 1");
        assert_eq!(location, "21:30425601-4611686018427387903", "{svlen}");
    }
}

/// The Allele Ensembl VEP 116.2 writes for the type shapes at the edges of
/// `get_SO_term`'s regexes: `mobile_element_insertion` for `<INS:ME::ALU>` (one colon
/// before the element), `copy_number_variation` for `<<CN0>>` (one bracket a side in the
/// copy-number test) and `deletion` for `<<DEL>>` and `<DEL:X>`; and the reference block
/// `<NON_REF>,<NON_REF>`, which Ensembl VEP skips as `NON_REF`, is one intergenic row per
/// allele, never an error.
#[test]
fn the_type_and_the_block_rows_follow_ensembl_vep() {
    let allele_of = |record: &str| -> Vec<String> {
        annotate_record(record, &[], &[])
            .iter()
            .filter(|row| row.contains("ENST00000390690"))
            .map(|row| row.split('\t').nth(2).unwrap().to_string())
            .collect()
    };
    for (record, allele) in [
        (
            "21\t30425650\tme\tN\t<INS:ME::ALU>\t.\t.\tSVTYPE=INS",
            "mobile_element_insertion",
        ),
        (
            "21\t30425650\tcn0\tN\t<<CN0>>\t.\t.\tEND=30425949",
            "copy_number_variation",
        ),
        (
            "21\t30425650\tdel\tN\t<<DEL>>\t.\t.\tEND=30425949",
            "deletion",
        ),
        (
            "21\t30425650\tdelx\tN\t<DEL:X>\t.\t.\tEND=30425949",
            "deletion",
        ),
    ] {
        assert_eq!(allele_of(record), [allele], "{record}");
    }
    let blocks = annotate_record(
        "21\t30425650\tblocks\tN\t<NON_REF>,<NON_REF>\t.\t.\tEND=30425700",
        &[],
        &[],
    );
    let block = row("21:30425651-30425700", "<NON_REF>", "intergenic_variant");
    assert_eq!(rows_on(&blocks, "blocks\t"), [block.clone(), block]);
}

/// The rows Ensembl VEP 116.2 writes on ENST00000390690 for the shapes Perl's `split` and
/// breakend regex admit: `A G,` is the one allele `G` (`missense_variant`); `<DEL>,` with
/// `SVTYPE=DEL` and `N.` with `SVTYPE=DEL/` are the deletion over END; `N[ 21 : 30425700 [`
/// is the own breakend `N.` beside the mate `N[ 21 : 30425700 [` as written, both
/// `feature_truncation,coding_sequence_variant`, as are `<BND>` with `CHR2= 21;END2=30425700`
/// (mate `N[ 21:30425700[`) and with `END2=30425700=5` (mate `N[21:30425700[`). The
/// gVCF site `<NON_REF>,G`, which Ensembl VEP drops, keeps `G`'s written ALLELE_NUM 2.
#[test]
fn trailing_fields_and_breakend_whitespace_follow_ensembl_vep() {
    let on = |record: &str, extra: &[&str]| -> Vec<(String, String, String, String)> {
        annotate_record(record, &[], extra)
            .iter()
            .filter(|row| row.contains("ENST00000390690"))
            .map(|row| {
                let c: Vec<&str> = row.split('\t').collect();
                (
                    c[1].to_string(),
                    c[2].to_string(),
                    c[6].to_string(),
                    c[13].to_string(),
                )
            })
            .collect()
    };
    let rows = on("21\t30425650\ttrailing\tA\tG,\t.\t.\t.", &[]);
    assert_eq!(rows.len(), 1);
    assert_eq!(
        (rows[0].1.as_str(), rows[0].2.as_str()),
        ("G", "missense_variant")
    );
    let deletion = row(
        "21:30425651-30425949",
        "deletion",
        "feature_truncation,coding_sequence_variant,5_prime_UTR_variant",
    );
    for record in [
        "21\t30425650\tdel_trailing\tN\t<DEL>,\t.\t.\tSVTYPE=DEL;END=30425949",
        "21\t30425650\tsvtype_slash\tN\tN.\t.\t.\tSVTYPE=DEL/;END=30425949",
    ] {
        let rows = on(record, &[]);
        assert_eq!(rows.len(), 1, "{record}");
        assert_eq!(
            (rows[0].0.clone(), rows[0].1.clone(), rows[0].2.clone()),
            deletion,
            "{record}"
        );
    }
    for (record, mate) in [
        (
            "21\t30425650\tws_bracket\tN\tN[ 21 : 30425700 [\t.\t.\tSVTYPE=BND",
            "N[ 21 : 30425700 [",
        ),
        (
            "21\t30425650\tws_chr2\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2= 21;END2=30425700",
            "N[ 21:30425700[",
        ),
        (
            "21\t30425650\tend2_eq\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700=5",
            "N[21:30425700[",
        ),
    ] {
        let rows = on(record, &["--allele_number"]);
        let alleles: Vec<(&str, &str)> = rows
            .iter()
            .map(|(_, allele, csq, extra)| {
                assert_eq!(
                    csq, "feature_truncation,coding_sequence_variant",
                    "{record}"
                );
                let num = extra
                    .split(';')
                    .find_map(|kv| kv.strip_prefix("ALLELE_NUM="))
                    .unwrap_or("");
                (allele.as_str(), num)
            })
            .collect();
        assert_eq!(alleles, [("N.", "0"), (mate, "1")], "{record}");
    }
    let site = on(
        "21\t30425650\tnonref_first\tA\t<NON_REF>,G\t.\t.\t.",
        &["--allele_number"],
    );
    assert_eq!(site.len(), 1);
    assert_eq!(site[0].1, "G");
    assert!(site[0].3.contains("ALLELE_NUM=2"), "{}", site[0].3);
}

/// The rows Ensembl VEP 116.2 writes for the coordinate and allele spellings it numifies
/// or upper-cases: a `<CNV:TR>` contraction's DISTANCE is the run's, 4484 for
/// `21 30430409 x N <CNV:TR> . . SVTYPE=CNV;END=30430417;SVLEN=8;RUS=CT;RUC=3` on
/// ENST00000390690 (30430410 - 30425926) and 75 for `21 30431100 x N <CNV:TR> . .
/// SVTYPE=CNV;END=30431108;SVLEN=8;RUS=AG;RUC=2` on ENST00000334068, not the lost
/// tail's 4490 and 79; a blank closing the data line stays in the last INFO value, so
/// `END2=30425700 ` names the mate `N[21:30425700 [`; a lowercase `g` is the Allele `G`
/// and `agc a` the JSON allele string `GC/-`; `3e7` is POS 30000000 (the deletion
/// `21:30000001-30425949`) and `30425949.0` is 30425949 (`upstream_gene_variant
/// DISTANCE=23` on ENST00000390690); a structural record at POS 2^62 - 1 fails its line;
/// and the JSON allele string keeps `*` at its place (`N/*/G`) and a `<NON_REF>` beside a
/// bracket ALT (`N/N[21:30425700[/<NON_REF>`).
#[test]
fn numified_and_upper_cased_spellings_follow_ensembl_vep() {
    let field = |record: &str, feature: &str, column: usize, extra: &[&str]| -> String {
        annotate_record(record, &[], extra)
            .iter()
            .find(|row| row.contains(feature))
            .map(|row| row.split('\t').nth(column).unwrap().to_string())
            .unwrap_or_else(|| panic!("{record}: no row on {feature}"))
    };
    let distance = |record: &str, feature: &str| -> String {
        field(record, feature, 13, &[])
            .split(';')
            .find_map(|kv| kv.strip_prefix("DISTANCE="))
            .unwrap_or_else(|| panic!("{record}: no DISTANCE on {feature}"))
            .to_string()
    };
    assert_eq!(
        distance(
            "21\t30430409\tloss\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30430417;SVLEN=8;RUS=CT;RUC=3",
            "ENST00000390690"
        ),
        "4484"
    );
    assert_eq!(
        distance(
            "21\t30431100\tloss\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30431108;SVLEN=8;RUS=AG;RUC=2",
            "ENST00000334068"
        ),
        "75"
    );
    let blank = "21\t30425650\tblank\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700 ";
    let alleles: Vec<String> = annotate_record(blank, &[], &[])
        .iter()
        .filter(|row| row.contains("ENST00000390690"))
        .map(|row| row.split('\t').nth(2).unwrap().to_string())
        .collect();
    assert_eq!(alleles, ["N.", "N[21:30425700 ["]);
    let info = annotate_record(blank, &["--vcf"], &[])
        .pop()
        .map(|line| line.split('\t').nth(7).unwrap().to_string())
        .unwrap();
    assert!(
        info.starts_with("SVTYPE=BND;CHR2=21;END2=30425700 ;CSQ="),
        "{info}"
    );
    assert_eq!(
        field(
            "21\t30425650\tlower\tA\tg\t.\t.\t.",
            "ENST00000390690",
            2,
            &[]
        ),
        "G"
    );
    let json_allele_string = |record: &str| -> String {
        let line = annotate_record(record, &["--json"], &[]).pop().unwrap();
        serde_json::from_str::<serde_json::Value>(&line).unwrap()["allele_string"]
            .as_str()
            .unwrap()
            .to_string()
    };
    assert_eq!(
        json_allele_string("21\t30425650\tlower\tagc\ta\t.\t.\t."),
        "GC/-"
    );
    assert_eq!(
        field(
            "21\t3e7\texp\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425949",
            "ENST00000390690",
            1,
            &["--max_sv_size", "-1"]
        ),
        "21:30000001-30425949"
    );
    let float = "21\t30425949.0\tfloat\tA\tG\t.\t.\t.";
    assert_eq!(field(float, "ENST00000390690", 1, &[]), "21:30425949");
    assert_eq!(distance(float, "ENST00000390690"), "23");
    assert!(annotate_record(
        "21\t4611686018427387903\tceiling\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-300",
        &[],
        &[]
    )
    .is_empty());
    assert_eq!(
        json_allele_string("21\t30425650\tstar\tN\t*,G\t.\t.\t."),
        "N/*/G"
    );
    assert_eq!(
        json_allele_string("21\t30425650\tblock\tN\tN[21:30425700[,<NON_REF>\t.\t.\t."),
        "N/N[21:30425700[/<NON_REF>"
    );
}

/// A record whose every ALT is `*` names no allele. Ensembl VEP 116.2 writes
/// `21 30425650 x N * . . .` with no row in the default format, carries the VCF line
/// without CSQ and writes a JSON object with `"allele_string":"N/*"`, start and end
/// 30425650 and `"most_severe_consequence":"?"`; with `SVTYPE=DEL;END=30425949` it
/// annotates a deletion from INFO alone, which vep-rs does not.
#[test]
fn a_star_only_record_is_carried_without_rows() {
    for record in [
        "21\t30425650\tstar\tN\t*\t.\t.\t.",
        "21\t30425650\tstar\tN\t*\t.\t.\tSVTYPE=DEL;END=30425949",
    ] {
        assert!(annotate_record(record, &[], &[]).is_empty(), "{record}");
        assert_eq!(
            annotate_record(record, &["--vcf"], &[]),
            [record],
            "{record}"
        );
        let line = annotate_record(record, &["--json"], &[]).pop().unwrap();
        let object: serde_json::Value = serde_json::from_str(&line).unwrap();
        assert_eq!(object["allele_string"], "N/*", "{record}");
        assert_eq!(object["start"], 30425650, "{record}");
        assert_eq!(object["end"], 30425650, "{record}");
        assert_eq!(object["most_severe_consequence"], "?", "{record}");
        assert_eq!(object["id"], "star", "{record}");
        assert!(object.get("transcript_consequences").is_none(), "{record}");
        assert!(object.get("intergenic_consequences").is_none(), "{record}");
    }
}

/// With a FASTA, Ensembl VEP 116.2 reads `21 30425587 x N <CNV:TR>,<CNV:TR> . .
/// SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4` as `AGCAGCAGC//AGCAGCAGCAGC`,
/// keeps the first minimised allele alone and writes `21:30425588-30425596 -
/// ENST00000390690 inframe_deletion ALLELE_NUM=1`: the RN=0 allele as the run's loss, the
/// real allele dropped. Both alleles are written here, the RN=0 allele as the loss of
/// the run and the four-unit allele as the gain.
#[test]
fn an_rn_zero_allele_is_the_loss_of_the_run_beside_its_sibling() {
    let rows: Vec<(String, String)> = annotate_record(
        "21\t30425587\trn0\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4",
        &[],
        &["--allele_number"],
    )
    .iter()
    .filter(|row| row.contains("ENST00000390690"))
    .map(|row| {
        let c: Vec<&str> = row.split('\t').collect();
        let num = c[13]
            .split(';')
            .find_map(|kv| kv.strip_prefix("ALLELE_NUM="))
            .unwrap()
            .to_string();
        (num, c[6].to_string())
    })
    .collect();
    assert_eq!(
        rows,
        [
            (
                "1".to_string(),
                "feature_truncation,inframe_deletion".to_string()
            ),
            (
                "2".to_string(),
                "feature_elongation,coding_sequence_variant".to_string()
            ),
        ]
    );
}

/// The rows Ensembl VEP 116.2 writes on ENST00000390690 for the shapes the parser reads
/// from INFO as Perl does: a mate named by CHR2 is `N[21:30425700[` beside the own `A.`
/// of a record whose REF is `A`, and keeps `END2=00030425700` as written; `SVLEN=0` or
/// `SVLEN=abc` is the pair `21:30425650-30425651`, `feature_elongation` for a
/// duplication and the region term alone for an inversion or a mate-less breakpoint;
/// `.N` beside `N[21:30425700[` lends the own breakend its base and writes no row of its
/// own; a repeated SVLEN reads its last value.
#[test]
fn the_pair_and_the_mate_text_follow_ensembl_vep() {
    let on_390690 = |record: &str| rows_on(&annotate_record(record, &[], &[]), "ENST00000390690");
    let truncating = "feature_truncation,coding_sequence_variant";
    assert_eq!(
        on_390690("21\t30425650\tref_a\tA\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700"),
        [
            row("21:30425651", "A.", truncating),
            row("21:30425651", "N[21:30425700[", truncating)
        ]
    );
    assert_eq!(
        on_390690("21\t30425650\tpadded\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=00030425700"),
        [
            row("21:30425651", "N.", truncating),
            row("21:30425651", "N[21:00030425700[", truncating)
        ]
    );
    assert_eq!(
        on_390690("21\t30425650\tmixed\tN\t.N,N[21:30425700[\t.\t.\tSVTYPE=BND"),
        [
            row("21:30425651", "N.", truncating),
            row("21:30425651", "N[21:30425700[", truncating)
        ]
    );
    for (alt, svlen, allele, csq) in [
        (
            "<DUP>",
            "0",
            "duplication",
            "feature_elongation,coding_sequence_variant",
        ),
        (
            "<DUP>",
            "abc",
            "duplication",
            "feature_elongation,coding_sequence_variant",
        ),
        (
            "<INS:ME:ALU>",
            "0",
            "Alu_insertion",
            "feature_elongation,coding_sequence_variant",
        ),
        ("<INV>", "0", "inversion", "coding_sequence_variant"),
        (
            "<CNV>",
            "abc",
            "copy_number_variation",
            "coding_sequence_variant",
        ),
        (
            "<BND>",
            "0",
            "chromosome_breakpoint",
            "coding_sequence_variant",
        ),
    ] {
        assert_eq!(
            on_390690(&format!(
                "21\t30425650\tpair\tN\t{alt}\t.\t.\tSVLEN={svlen}"
            )),
            [row("21:30425650-30425651", allele, csq)],
            "{alt} SVLEN={svlen}"
        );
    }
    assert_eq!(
        on_390690("21\t30425650\trepeated\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-100;SVLEN=-300"),
        [row(
            "21:30425651-30425950",
            "deletion",
            "feature_truncation,coding_sequence_variant,5_prime_UTR_variant"
        )]
    );
}

/// The rows of a breakend record follow what it names, as Ensembl VEP 116.2 writes
/// them (`Parser/VCF.pm:484-499`): a `<BND>` with SVLEN and no CHR2 is one allele
/// over its span, Allele `chromosome_breakpoint`, `coding_sequence_variant,
/// 5_prime_UTR_variant` on ENST00000390690 and no `feature_truncation`; an INFO/END
/// never ends a breakpoint's span (`:487-489`), so `END=30425750` is the point
/// `21:30425651` and `CHR2=2;END=40000000` the point with its mate on chromosome 2
/// (two local rows, not every transcript to 40 Mb); POS2 and an END2 that is not a
/// coordinate name no mate (`:493`), so those records print `chromosome_breakpoint`
/// alone.
#[test]
fn breakend_rows_follow_what_the_record_names() {
    let rows_for = |info: &str| {
        annotate_record(
            &format!("21\t30425650\tbnd\tN\t<BND>\t.\t.\t{info}"),
            &[],
            &[],
        )
    };
    let on_390690 = |rows: &[String]| rows_on(rows, "ENST00000390690");
    let point = row(
        "21:30425651",
        "chromosome_breakpoint",
        "coding_sequence_variant",
    );
    let span = rows_for("SVTYPE=BND;SVLEN=300");
    assert_eq!(
        on_390690(&span),
        [row(
            "21:30425651-30425950",
            "chromosome_breakpoint",
            "coding_sequence_variant,5_prime_UTR_variant"
        )]
    );
    let end_only = rows_for("SVTYPE=BND;END=30425750");
    assert_eq!(on_390690(&end_only), std::slice::from_ref(&point));
    let far_mate = rows_for("SVTYPE=BND;CHR2=2;END=40000000");
    assert_eq!(far_mate.len(), 2, "{far_mate:?}");
    assert_eq!(
        on_390690(&far_mate),
        [row(
            "21:30425651",
            "N.",
            "feature_truncation,coding_sequence_variant"
        )]
    );
    for info in [
        "SVTYPE=BND;CHR2=2;POS2=40000000",
        "SVTYPE=BND;CHR2=21;END2=-7",
    ] {
        assert_eq!(
            on_390690(&rows_for(info)),
            std::slice::from_ref(&point),
            "{info}"
        );
    }
}

/// INFO/SVTYPE names the class of an ALT outside the VCF 4.4 forms
/// (`Parser/VCF.pm:472-475`): `N.` with `SVTYPE=DEL` and END is a deletion over its
/// span, `deletion` with `frameshift_variant,feature_truncation` on ENST00000390690,
/// as Ensembl VEP 116.2 writes it. A comma-separated SVLEN reads its first value
/// (BaseVCF4.pm:304) and a tandem repeat its largest (`Parser/VCF.pm:505-511`), the
/// Locations Ensembl VEP 116.2 writes for `SVLEN=-300,400` and `SVLEN=300,400`.
#[test]
fn svtype_and_a_comma_separated_svlen_are_read_as_ensembl_vep_reads_them() {
    let location_and_terms =
        |record: &str| rows_on(&annotate_record(record, &[], &[]), "ENST00000390690");
    assert_eq!(
        location_and_terms("21\t30425600\tgnomad_del\tN\tN.\t.\t.\tSVTYPE=DEL;END=30425899"),
        [row(
            "21:30425601-30425899",
            "deletion",
            "frameshift_variant,feature_truncation"
        )]
    );
    let cnv =
        location_and_terms("21\t30425600\tcnv\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;SVLEN=-300,400");
    assert_eq!(cnv.len(), 2);
    assert!(
        cnv.iter()
            .all(|(location, _, _)| location == "21:30425601-30425900"),
        "{cnv:?}"
    );
    let repeat = location_and_terms(
        "21\t30425600\ttr\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;SVLEN=300,400;RUS=AT,AT;RUC=100,50",
    );
    assert_eq!(repeat.len(), 2);
    assert!(
        repeat
            .iter()
            .all(|(location, _, _)| location == "21:30425601-30426000"),
        "{repeat:?}"
    );
}

/// A record the parser refuses is skipped with a warning the user sees, as Ensembl
/// VEP's `skipped_variant_msg` reports one, and the rest of the file is annotated.
#[test]
fn a_refused_record_is_skipped_with_a_visible_warning() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30425600\tbad_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=-5\n\
         21\t-1\tbad_pos\tA\tG\t.\t.\t.\n\
         21\t30425600\tkept\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425899\n",
    )
    .unwrap();
    let output = dir.path().join("out");
    let cache = format!(
        "{}/../../tests/golden/116/GRCh38-release116/json_cache",
        env!("CARGO_MANIFEST_DIR")
    );
    let result = vep_binary()
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .args(["--json_cache", &cache, "--offline", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats"])
        .output()
        .unwrap();
    assert!(result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("WARN") && stderr.contains("Invalid END: '-5'"),
        "{stderr}"
    );
    assert!(stderr.contains("Invalid POS: '-1'"), "{stderr}");
    let rows: Vec<String> = std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    assert!(
        rows.iter().all(|r| r.starts_with("kept\t")) && rows.len() == 2,
        "{rows:?}"
    );
}

/// The `--max_sv_size` row of `docs/cli-reference.md` states the documented
/// exclusion: Ensembl VEP 116 drops an oversize record and vep-rs annotates it.
#[test]
fn max_sv_size_exclusion_is_documented() {
    let page = std::fs::read_to_string(format!(
        "{}/../../docs/cli-reference.md",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let row = page
        .lines()
        .find(|l| l.starts_with("| `--max_sv_size`"))
        .expect("docs/cli-reference.md has a --max_sv_size row");
    for phrase in [
        "Ensembl VEP 116 drops an oversize record before annotation and writes nothing for it in any format",
        "vep-rs annotates it",
        "intended-divergences.md",
    ] {
        assert!(row.contains(phrase), "the --max_sv_size row lacks {phrase:?}: {row}");
    }
}

#[test]
fn test_dont_skip_fails_on_invalid_vcf_line() {
    let output_dir = tempfile::tempdir().unwrap();
    let input_path = output_dir.path().join("invalid.vcf");
    let output_path = output_dir.path().join("output.txt");

    let vcf = "##fileformat=VCFv4.3\n\
               #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
               21\t100\t.\tA\n";
    std::fs::write(&input_path, vcf).unwrap();

    let result = vep_binary()
        .args(["-i", input_path.to_str().unwrap()])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--dont_skip", "--quiet"])
        .output()
        .unwrap();

    assert!(
        !result.status.success(),
        "Expected VEP to fail with --dont_skip on invalid VCF line"
    );
}

#[test]
fn test_allow_non_variant_passthrough_vcf() {
    let output_dir = tempfile::tempdir().unwrap();
    let input_path = output_dir.path().join("non_variant.vcf");
    let output_path = output_dir.path().join("output.vcf");

    let data_line = "21\t100\t.\tA\t.\t.\tPASS\t.";
    let vcf = format!(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
        data_line
    );
    std::fs::write(&input_path, vcf).unwrap();

    let result = vep_binary()
        .args(["-i", input_path.to_str().unwrap()])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--vcf", "--allow_non_variant", "--force", "--quiet"])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "VEP failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let output = std::fs::read_to_string(&output_path).unwrap();
    assert!(
        output.lines().any(|l| l == data_line),
        "Expected non-variant line to be passed through"
    );
}

#[test]
fn test_unknown_plugin_warns_but_does_not_fail() {
    // An unknown plugin name should warn but not crash the process.
    let dir = tempfile::tempdir().unwrap();
    let output_path = dir.path().join("output.txt");

    let result = vep_binary()
        .args(["-i", &test_input_vcf()])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--force", "--quiet", "--offline"])
        .args(["--plugin", "NonExistentPlugin,param1"])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "VEP should not fail for unknown plugins: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn test_builtin_plugin_names_are_recognized() {
    // A builtin plugin name resolves without an unknown-plugin warning; with no
    // data-file params it fails at init with a plugin error instead.
    let dir = tempfile::tempdir().unwrap();
    let output_path = dir.path().join("output.txt");

    let result = vep_binary()
        .args(["-i", &test_input_vcf()])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--force", "--quiet", "--offline"])
        .args(["--plugin", "CADD"])
        .output()
        .unwrap();

    let stderr = String::from_utf8_lossy(&result.stderr);
    // CADD requires file params, so it should fail with an init error
    // containing "CADD" (not "not found as a built-in").
    assert!(
        !stderr.contains("not found as a built-in"),
        "CADD should be recognized as a builtin plugin"
    );
}

fn duckdb_available() -> bool {
    std::process::Command::new("duckdb")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

#[test]
fn parquet_round_trip_nested() {
    // VCF in, TSV intermediate, DuckDB finalize, Parquet directory, schema probe
    // via `duckdb -c`. Skipped when `duckdb` is not on PATH.
    if !duckdb_available() {
        eprintln!("Skipping parquet_round_trip_nested: duckdb not on PATH");
        return;
    }
    let json_cache = match json_cache_path() {
        Some(p) => p,
        None => {
            eprintln!(
                "Skipping parquet_round_trip_nested: VEP_RS_TEST_CACHE not set or directory absent"
            );
            return;
        }
    };

    let vcf = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        21\t25000100\trs1\tA\tG\t.\tPASS\t.\n\
        21\t25000200\t.\tAAA\tA\t.\tPASS\t.\n\
        21\t25000300\t.\tT\tTGGG\t.\tPASS\t.\n";
    let dir = tempfile::tempdir().unwrap();
    let input_path = dir.path().join("in.vcf");
    std::fs::write(&input_path, vcf).unwrap();
    let output_path = dir.path().join("out.parquet");

    // `--no_headers` is passed deliberately: the runner special-cases parquet so
    // the TSV intermediate always carries a header, or DuckDB's
    // `read_csv(header=true)` would consume the first data row as column names.
    let result = vep_binary()
        .args(["-i", input_path.to_str().unwrap()])
        .args(["-o", output_path.to_str().unwrap()])
        .args(["--output_format", "parquet"])
        .args(["--parquet_shape", "nested"])
        .args(["--json_cache", &json_cache])
        .args([
            "--offline",
            "--force",
            "--no_stats",
            "--quiet",
            "--no_headers",
            "--fork",
            "1",
        ])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "Parquet run failed: stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );

    assert!(
        output_path.is_dir(),
        "Parquet output path should be a directory: {:?}",
        output_path
    );

    // Probe schema via duckdb: Gene must be a LIST<VARCHAR> in nested shape.
    let probe_sql = format!(
        "SELECT typeof(Gene) FROM read_parquet('{}/**/*.parquet') LIMIT 1",
        output_path.display()
    );
    let probe = std::process::Command::new("duckdb")
        .args(["-csv", "-c", &probe_sql])
        .output()
        .unwrap();
    assert!(
        probe.status.success(),
        "duckdb schema probe failed: {}",
        String::from_utf8_lossy(&probe.stderr)
    );
    let probe_stdout = String::from_utf8_lossy(&probe.stdout);
    assert!(
        probe_stdout.contains("[]") || probe_stdout.contains("LIST"),
        "Expected LIST type for nested Gene column, got: {}",
        probe_stdout
    );

    // Run the adapter against the Parquet dir; require the rsID for
    // the first variant to surface as Uploaded_variation (rsID branch
    // of variant.rs::uploaded_variation).
    let adapter = format!(
        "{}/../../scripts/adapters/parquet_to_vep_tab.py",
        env!("CARGO_MANIFEST_DIR")
    );
    let adapter_out = dir.path().join("adapter.tsv");
    let adapter_result = std::process::Command::new("python3")
        .arg(&adapter)
        .args(["--input", output_path.to_str().unwrap()])
        .args(["--output", adapter_out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        adapter_result.status.success(),
        "Adapter failed: stderr={}",
        String::from_utf8_lossy(&adapter_result.stderr)
    );
    let adapter_text = std::fs::read_to_string(&adapter_out).unwrap();
    assert!(
        adapter_text.contains("rs1"),
        "Adapter output should contain rsID 'rs1' as Uploaded_variation, got:\n{}",
        adapter_text
    );
}

/// A region-format deletion, `chr:start-end:strand/-`, is annotated against the
/// reference FASTA's bases over its span, as Ensembl VEP 116.2 annotates it: each
/// of the two below writes the rows the GRCh37-hgvs golden corpus expects for the
/// same deletion written as a VCF record (`829531`, `rs149363218`), consequences
/// and HGVS alike, under the name the region line gives it (`21_27062184_-/-`).
#[test]
fn region_deletion_is_annotated_against_the_reference_bases() {
    use std::io::Read;

    let corpus = format!(
        "{}/../../tests/golden/116/GRCh37-hgvs",
        env!("CARGO_MANIFEST_DIR")
    );
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("reference.fa");
    let mut bases = Vec::new();
    flate2::read::GzDecoder::new(std::fs::File::open(format!("{corpus}/reference.fa.gz")).unwrap())
        .read_to_end(&mut bases)
        .unwrap();
    std::fs::write(&fasta, bases).unwrap();
    std::fs::copy(
        format!("{corpus}/reference.fa.fai"),
        dir.path().join("reference.fa.fai"),
    )
    .unwrap();
    let input = dir.path().join("in.txt");
    std::fs::write(
        &input,
        "21:27062184-27062184:1/-\n21:10910401-10910402:1/-\n",
    )
    .unwrap();
    let output = dir.path().join("out");
    let result = vep_binary()
        .args(["-i", input.to_str().unwrap(), "--format", "region"])
        .args(["-o", output.to_str().unwrap()])
        .args(["--json_cache", &format!("{corpus}/json_cache")])
        .args(["--offline", "--assembly", "GRCh37", "--hgvs"])
        .args(["--fasta", fasta.to_str().unwrap()])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let mut expected_text = String::new();
    flate2::read::GzDecoder::new(
        std::fs::File::open(format!("{corpus}/expected/default.txt.gz")).unwrap(),
    )
    .read_to_string(&mut expected_text)
    .unwrap();
    let actual_text = std::fs::read_to_string(&output).unwrap();
    for (vcf_id, region_name) in [
        ("829531", "21_27062184_-/-"),
        ("rs149363218", "21_10910401_-/-"),
    ] {
        let rows_named = |text: &str, name: &str| -> Vec<String> {
            text.lines()
                .filter(|l| l.split('\t').next() == Some(name))
                .map(|l| l.split_once('\t').unwrap().1.to_string())
                .collect()
        };
        let expected = rows_named(&expected_text, vcf_id);
        assert!(!expected.is_empty(), "{vcf_id} in the golden output");
        assert_eq!(
            rows_named(&actual_text, region_name),
            expected,
            "{region_name}"
        );
    }
}
