// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The flags that read data Ensembl ships for human alone, run against a cache
//! of another species through the `vep` binary: refused where Ensembl VEP
//! refuses them (SIFT, PolyPhen, `--af_gnomad`), switched off where it switches
//! them off (`--everything`, `--safe`), and otherwise as the flag table states: a column flag
//! declares its columns and leaves every cell empty, a transcript-set flag drops every transcript,
//! a frequency filter whose population is accepted keeps or drops every record and `--regulatory`
//! writes no regulatory row, as Ensembl VEP writes a flag whose data the cache lacks. The cache is
//! the synthetic one under `fixtures/species/json_cache` (its README describes it): a species that
//! is not `homo_sapiens`, no SIFT or PolyPhen, variation columns without a frequency column, and
//! transcripts carrying none of the GENCODE, MANE, APPRIS, TSL or CCDS attributes. Each test cites
//! the Ensembl VEP source (release/116.2) for the behaviour it proves.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The synthetic cache of another species.
fn species_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/species/json_cache")
}

/// Five records on the cache's one chromosome: upstream of both transcripts, a
/// missense change inside the coding transcript (co-located with the cache's
/// one known variant, `rs1`), an intronic base, a 3' UTR base and a base no
/// transcript lies within 5,000 bases of.
const RECORDS: &str = "\
1\t400\tup\tG\tC\t.\t.\t.
1\t1014\tmis\tC\tT\t.\t.\t.
1\t1100\tintron\tA\tG\t.\t.\t.
1\t1250\tutr3\tG\tA\t.\t.\t.
1\t50000\tfar\tG\tA\t.\t.\t.
";

/// The default-format rows the five records annotate to without a filter: the
/// coding transcript (`SYNT0000001`) and the non-coding one (`SYNT0000002`)
/// give two rows to each of the first four, the fifth is intergenic.
const ROWS_WITHOUT_A_FILTER: usize = 9;

fn write_input(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    fs::write(
        &path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORDS}"),
    )
    .unwrap();
    path
}

/// Runs the binary on `in.vcf` under `dir` against the species cache with
/// `extra`, writing `<dir>/<label>.out` and `<dir>/<label>_warnings.txt`;
/// nothing is asserted about the exit status. The notices (`INFO:` lines) are
/// printed, so `--quiet` is not passed.
fn run(dir: &Path, label: &str, extra: &[&str]) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.out")))
        .arg("--warning_file")
        .arg(dir.join(format!("{label}_warnings.txt")))
        .arg("--json_cache")
        .arg(species_cache())
        .args([
            "--species",
            "synthetic_species",
            "--assembly",
            "SYN1",
            "--offline",
            "--no_stats",
            "--force_overwrite",
        ])
        .args(extra)
        .stdin(Stdio::null());
    cmd.output().unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// [`run`] asserting success; the output text, the standard error and the
/// warnings file's text (empty when the run wrote none).
fn annotate(dir: &Path, label: &str, extra: &[&str]) -> (String, String, String) {
    let output = run(dir, label, extra);
    let err = stderr(&output);
    assert!(output.status.success(), "{extra:?}: {err}");
    let text = fs::read_to_string(dir.join(format!("{label}.out"))).unwrap();
    let warnings =
        fs::read_to_string(dir.join(format!("{label}_warnings.txt"))).unwrap_or_default();
    (text, err, warnings)
}

/// Asserts the run exits 1 with `message` as its one `ERROR: [cache]` line and
/// writes no output file.
fn assert_cache_refusal(dir: &Path, label: &str, extra: &[&str], message: &str) {
    let output = run(dir, label, extra);
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(1), "{extra:?}: {err}");
    let error_lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(error_lines, [format!("ERROR: [cache] {message}")], "{err}");
    assert!(
        !dir.join(format!("{label}.out")).exists(),
        "{extra:?}: no output for a refused run"
    );
}

fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// The `#Uploaded_variation ...` column names of a `--tab` output.
fn tab_columns(text: &str) -> Vec<&str> {
    text.lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .expect("a --tab output carries its column header")
        .trim_start_matches('#')
        .split('\t')
        .collect()
}

/// The `## <field> : <description>` header lines name `field`.
fn describes_field(text: &str, field: &str) -> bool {
    let prefix = format!("## {field} : ");
    text.lines().any(|line| line.starts_with(&prefix))
}

/// The `KEY=` names of a default-format row's Extra column.
fn extra_keys(row: &str) -> Vec<&str> {
    row.split('\t')
        .nth(13)
        .unwrap_or("")
        .split(';')
        .filter_map(|pair| pair.split('=').next())
        .collect()
}

/// Asserts a `--tab` output names every `field` as a column, describes each in
/// its header, and leaves its cell empty (`-`) on every row.
fn assert_columns_present_and_empty(text: &str, fields: &[&str], what: &str) {
    let columns = tab_columns(text);
    let rows = data_lines(text);
    assert_eq!(
        rows.len(),
        ROWS_WITHOUT_A_FILTER,
        "{what}: every row:\n{text}"
    );
    for field in fields {
        assert!(
            describes_field(text, field),
            "{what}: the header describes {field}:\n{text}"
        );
        let at = columns
            .iter()
            .position(|column| column == field)
            .unwrap_or_else(|| panic!("{what}: the column header names {field}:\n{text}"));
        for row in &rows {
            let cell = row.split('\t').nth(at).unwrap_or("");
            assert_eq!(cell, "-", "{what}: {field} is empty on {row}");
        }
    }
}

/// `--sift` on a cache whose metadata names no SIFT version is refused with
/// Ensembl VEP's text, `ERROR: SIFT not available`
/// (`AnnotationSource/Cache/Transcript.pm:148-176`: `check_sift_polyphen`
/// throws unless `everything` or `safe` is set), as `--polyphen` is for
/// PolyPhen, before the input is read and under `--show_cache_info` too, which
/// builds the annotation sources before it prints (BaseRunner.pm:415). Under
/// `--everything` or `--safe` the tool is switched off instead with
/// `INFO: disabling SIFT` (`:161-166`), so the run annotates with no SIFT or
/// PolyPhen column and no `## SIFT` or `## PolyPhen` header line.
#[test]
fn sift_refuses_on_a_cache_without_sift_and_everything_disables_it() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    assert_cache_refusal(dir.path(), "sift", &["--sift", "b"], "SIFT not available");
    assert_cache_refusal(
        dir.path(),
        "polyphen",
        &["--polyphen", "b"],
        "PolyPhen not available",
    );
    assert_cache_refusal(
        dir.path(),
        "info",
        &["--sift", "b", "--show_cache_info"],
        "SIFT not available",
    );

    let (text, err, _) = annotate(dir.path(), "everything", &["--everything"]);
    for notice in ["INFO: disabling SIFT", "INFO: disabling PolyPhen"] {
        assert!(
            err.lines().any(|line| line == notice),
            "--everything prints {notice}:\n{err}"
        );
    }
    assert!(!err.contains("ERROR:"), "{err}");
    assert_eq!(data_lines(&text).len(), ROWS_WITHOUT_A_FILTER, "{text}");
    for field in ["SIFT", "PolyPhen"] {
        assert!(
            !describes_field(&text, field),
            "no {field} header line once the tool is disabled:\n{text}"
        );
        assert!(
            data_lines(&text)
                .iter()
                .all(|row| !extra_keys(row).contains(&field)),
            "no {field} key in any row once the tool is disabled:\n{text}"
        );
    }

    let (text, err, _) = annotate(
        dir.path(),
        "safe",
        &["--sift", "b", "--polyphen", "b", "--safe"],
    );
    for notice in ["INFO: disabling SIFT", "INFO: disabling PolyPhen"] {
        assert!(
            err.lines().any(|line| line == notice),
            "--safe prints {notice}:\n{err}"
        );
    }
    assert!(
        !describes_field(&text, "SIFT") && !describes_field(&text, "PolyPhen"),
        "{text}"
    );
}

/// `--af_gnomad` on a cache whose declared variation columns lack `gnomADe` is
/// refused with Ensembl VEP's text, `ERROR: gnomAD data is not available in
/// this cache` (CacheDir.pm:188-197), while `--af_gnomade` and `--af_gnomadg`
/// run: their columns are declared in the header (Constants.pm:121-122 lists
/// the fields of each flag whatever the cache holds) and every cell stays
/// empty, since `add_colocated_frequency_data` writes a frequency key only when
/// the co-located variant carries it (OutputFactory.pm:1208). The co-located
/// variant itself is reported, so the empty cells are the frequencies' absence
/// and not the variant's.
#[test]
fn af_gnomad_refuses_and_af_gnomade_prints_empty_cells() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    assert_cache_refusal(
        dir.path(),
        "af_gnomad",
        &["--af_gnomad"],
        "gnomAD data is not available in this cache",
    );

    let gnomade = [
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
    ];
    let (text, _, _) = annotate(dir.path(), "af_gnomade", &["--af_gnomade", "--tab"]);
    assert_columns_present_and_empty(&text, &gnomade, "--af_gnomade");
    let colocated: Vec<&str> = data_lines(&text)
        .into_iter()
        .filter(|row| row.starts_with("mis\t"))
        .map(|row| row.split('\t').nth(12).unwrap())
        .collect();
    assert_eq!(
        colocated,
        ["rs1", "rs1"],
        "the known variant is reported:\n{text}"
    );

    let gnomadg = [
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
    ];
    let (text, _, _) = annotate(dir.path(), "af_gnomadg", &["--af_gnomadg", "--tab"]);
    assert_columns_present_and_empty(&text, &gnomadg, "--af_gnomadg");
}

/// Every frequency and attribute flag reading data the cache lacks declares
/// its fields in the header (Constants.pm:85-128, `@FLAG_FIELDS`) and writes
/// nothing into them: the frequency keys are skipped on a co-located variant
/// without them (OutputFactory.pm:1208), `MAX_AF` needs a frequency to compare
/// (`:1252-1259`), and each transcript column reads an attribute or xref the
/// transcript does not carry. In the default format the Extra column of every
/// row lacks the keys; under `--tab` every column is present and `-`.
#[test]
fn frequency_and_attribute_columns_print_headers_and_empty_cells() {
    let flags = [
        "--af",
        "--af_1kg",
        "--af_gnomade",
        "--af_gnomadg",
        "--max_af",
        "--pubmed",
        "--var_synonyms",
        "--mane",
        "--mane_select",
        "--appris",
        "--tsl",
        "--ccds",
        "--gene_phenotype",
        "--uniprot",
        "--domains",
    ];
    let fields = [
        "AF",
        "AFR_AF",
        "AMR_AF",
        "EAS_AF",
        "EUR_AF",
        "SAS_AF",
        "gnomADe_AF",
        "gnomADe_NFE_AF",
        "gnomADg_AF",
        "gnomADg_NFE_AF",
        "MAX_AF",
        "MAX_AF_POPS",
        "PUBMED",
        "VAR_SYNONYMS",
        "MANE",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "APPRIS",
        "TSL",
        "CCDS",
        "GENE_PHENO",
        "SWISSPROT",
        "TREMBL",
        "UNIPARC",
        "UNIPROT_ISOFORM",
        "DOMAINS",
    ];
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());

    let (text, err, _) = annotate(dir.path(), "default", &flags);
    assert!(!err.contains("ERROR:"), "{err}");
    let rows = data_lines(&text);
    assert_eq!(rows.len(), ROWS_WITHOUT_A_FILTER, "{text}");
    for field in fields {
        assert!(
            describes_field(&text, field),
            "the header describes {field}:\n{text}"
        );
        for row in &rows {
            assert!(
                !extra_keys(row).contains(&field),
                "{field} is absent from the Extra column of {row}"
            );
        }
    }
    let colocated = rows
        .iter()
        .filter(|row| row.starts_with("mis\t"))
        .filter(|row| row.split('\t').nth(12) == Some("rs1"))
        .count();
    assert_eq!(
        colocated, 2,
        "the known variant is reported without frequencies:\n{text}"
    );

    let mut tab = flags.to_vec();
    tab.push("--tab");
    let (text, _, _) = annotate(dir.path(), "tab", &tab);
    assert_columns_present_and_empty(&text, &fields, "--tab");
}

/// `--gencode_basic` on a cache none of whose transcripts carries the
/// `gencode_basic` attribute drops every transcript as the cache loads, as
/// Ensembl VEP's `filter_transcript` drops each one
/// (`AnnotationType/Transcript.pm:222-225`), so every record annotates as
/// `intergenic_variant`; the run says so once, through the warnings file
/// (Ensembl VEP annotates silently). `--gencode_primary` does the same for its
/// attribute, after the assembly warning Ensembl VEP prints for it
/// (Config.pm:721-724), and `--flag_gencode_primary` keeps every transcript and
/// writes an empty `GENCODE_PRIMARY` column.
#[test]
fn gencode_basic_on_a_cache_without_the_attribute_warns_once_and_keeps_every_variant_intergenic() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let records = RECORDS.lines().count();
    let filtered_warning = |attribute: &str| {
        format!(
            "WARNING: no transcript of the cache carries the {attribute} attribute; \
             every transcript is filtered"
        )
    };

    let (text, err, warnings) = annotate(dir.path(), "basic", &["--gencode_basic"]);
    let rows = data_lines(&text);
    assert_eq!(
        rows.len(),
        records,
        "one intergenic row per record:\n{text}"
    );
    for row in &rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns[4], "-", "Feature of {row}");
        assert_eq!(columns[6], "intergenic_variant", "Consequence of {row}");
    }
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [filtered_warning("gencode_basic").as_str()],
        "one warning, in the warnings file"
    );
    assert_eq!(
        err.lines()
            .filter(|line| line.starts_with("WARNING:"))
            .count(),
        1,
        "the warning is echoed once on standard error:\n{err}"
    );

    let (text, _, warnings) = annotate(dir.path(), "primary", &["--gencode_primary"]);
    assert!(
        data_lines(&text)
            .iter()
            .all(|row| row.split('\t').nth(6) == Some("intergenic_variant")),
        "{text}"
    );
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [
            "WARNING: --gencode_primary option is currently only available for human on the GRCh38 assembly",
            filtered_warning("gencode_primary").as_str(),
        ],
        "the assembly warning, then the one filter warning"
    );

    let (text, _, warnings) = annotate(dir.path(), "flagged", &["--flag_gencode_primary", "--tab"]);
    assert_columns_present_and_empty(&text, &["GENCODE_PRIMARY"], "--flag_gencode_primary");
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        ["WARNING: --flag_gencode_primary option is currently only available for human on the GRCh38 assembly"],
        "the column flag filters nothing and warns about the assembly alone"
    );
}

/// `--check_frequency` on a cache whose variation columns hold no frequency
/// finds no frequency data for any record, so under `--freq_filter exclude`
/// every record is kept (BaseCacheVariation.pm:340-344: "no frequency data so
/// we have to include it") and under `include` every record is dropped
/// (`:321-322`); `1KG_ALL` is the one population the check accepts on such a
/// cache (`:285-288`: the others must be declared columns), so `--filter_common`
/// filters nothing too. The `FREQS` column is declared and empty.
#[test]
fn check_frequency_without_frequencies_filters_nothing() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let filter = [
        "--check_frequency",
        "--freq_pop",
        "1KG_ALL",
        "--freq_freq",
        "0.01",
        "--freq_gt_lt",
        "gt",
        "--freq_filter",
        "exclude",
        "--tab",
    ];
    let (text, err, _) = annotate(dir.path(), "exclude", &filter);
    assert!(!err.contains("ERROR:"), "{err}");
    assert_columns_present_and_empty(&text, &["FREQS"], "--freq_filter exclude");
    let kept: Vec<&str> = data_lines(&text)
        .iter()
        .map(|row| row.split('\t').next().unwrap())
        .collect();
    assert_eq!(
        kept,
        ["up", "up", "mis", "mis", "intron", "intron", "utr3", "utr3", "far"],
        "every record is kept:\n{text}"
    );

    let (text, _, _) = annotate(dir.path(), "common", &["--filter_common", "--tab"]);
    assert_eq!(data_lines(&text).len(), ROWS_WITHOUT_A_FILTER, "{text}");

    let include: Vec<&str> = filter
        .iter()
        .map(|flag| if *flag == "exclude" { "include" } else { flag })
        .collect();
    let (text, _, _) = annotate(dir.path(), "include", &include);
    assert!(
        data_lines(&text).is_empty(),
        "every record lacks frequency data, so none passes an include filter:\n{text}"
    );
    assert!(
        tab_columns(&text).contains(&"FREQS"),
        "the header is still written:\n{text}"
    );

    assert_cache_refusal(
        dir.path(),
        "population",
        &["--check_frequency", "--freq_pop", "gnomADe_NFE"],
        "Invalid population gnomADe_NFE specified",
    );
}
