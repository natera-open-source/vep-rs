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

/// The release 116 golden corpus's own pruned cache, so a test can annotate
/// chromosome 21 records without an environment variable.
fn release116_corpus_cache() -> String {
    format!(
        "{}/../../tests/golden/116/GRCh38-release116/json_cache",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// Annotates one VCF record against the release 116 corpus cache in the given
/// format and returns the output's data lines.
fn annotate_record(record: &str, format_flags: &[&str], extra: &[&str]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    let output = dir.path().join("out");
    let result = vep_binary()
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .args(["--json_cache", &release116_corpus_cache()])
        .args([
            "--offline",
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
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

/// A structural variant wider than `--max_sv_size` (the 10,600,000-base
/// `del_10_6_mb` of `tests/golden/116/GRCh38-release116` under the default
/// limit of 10,000,000) is annotated, not dropped: the default output lists its
/// transcript rows, the VCF output carries its line without `CSQ`, and the JSON
/// output omits it. Ensembl VEP 116.2 writes nothing for the record in any
/// format (`Parser.pm` `validate_vf` returns 0 at the size test); this is the
/// documented divergence of `--max_sv_size`. With `-1` the limit is lifted and
/// every format annotates the record, as Ensembl VEP 116.2 does under `-1`
/// (the corpus runs with it: 823 reference rows for the record).
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
