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
//! the Ensembl VEP source (release/116.3) for the behaviour it proves.
//!
//! The release caches differ in what they hold beside their transcripts, so the flags that read
//! variants, SIFT and regulatory data are also run on five copies of that cache, each cut or
//! extended to one combination the releases ship ([`SHAPES`]): transcripts alone, variants,
//! variants with SIFT, regulatory data, and variants with SIFT and regulatory data. No copy holds
//! PolyPhen, which Ensembl ships for human alone.

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

/// Runs the binary on `in.vcf` under `dir` against `cache`, the species cache
/// or a copy of it, with `extra`, writing `<dir>/<label>.out` and
/// `<dir>/<label>_warnings.txt`; nothing is asserted about the exit status.
/// The notices (`INFO:` lines) are printed, so `--quiet` is not passed.
fn run_on(dir: &Path, cache: &Path, label: &str, extra: &[&str]) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.out")))
        .arg("--warning_file")
        .arg(dir.join(format!("{label}_warnings.txt")))
        .arg("--json_cache")
        .arg(cache)
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

/// [`run_on`] the species cache asserting success; the output text, the
/// standard error and the warnings file's text (empty when the run wrote none).
fn annotate(dir: &Path, label: &str, extra: &[&str]) -> (String, String, String) {
    annotate_on(dir, &species_cache(), label, extra)
}

/// [`annotate`] against `cache`.
fn annotate_on(dir: &Path, cache: &Path, label: &str, extra: &[&str]) -> (String, String, String) {
    let output = run_on(dir, cache, label, extra);
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
    assert_cache_refusal_on(dir, &species_cache(), label, extra, message);
}

/// [`assert_cache_refusal`] against `cache`.
fn assert_cache_refusal_on(dir: &Path, cache: &Path, label: &str, extra: &[&str], message: &str) {
    let output = run_on(dir, cache, label, extra);
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

/// What a copy of the species cache holds beside its transcripts.
#[derive(Clone, Copy, Debug)]
struct Shape {
    /// The directory the copy is written to.
    name: &'static str,
    /// The fixture's variation file and its declared columns, none a frequency.
    variants: bool,
    /// A SIFT version in the metadata and a SIFT matrix on the coding transcript.
    sift: bool,
    /// A regulatory tree holding one promoter over the coding transcript's start.
    regulatory: bool,
}

const TRANSCRIPTS_ALONE: Shape = Shape {
    name: "transcripts_alone",
    variants: false,
    sift: false,
    regulatory: false,
};
const VARIANTS: Shape = Shape {
    name: "variants",
    variants: true,
    sift: false,
    regulatory: false,
};
const VARIANTS_AND_SIFT: Shape = Shape {
    name: "variants_and_sift",
    variants: true,
    sift: true,
    regulatory: false,
};
const REGULATORY: Shape = Shape {
    name: "regulatory",
    variants: false,
    sift: false,
    regulatory: true,
};
const VARIANTS_SIFT_AND_REGULATORY: Shape = Shape {
    name: "variants_sift_and_regulatory",
    variants: true,
    sift: true,
    regulatory: true,
};

/// The combinations of variant, SIFT and regulatory data the Ensembl and
/// Ensembl Genomes release caches ship.
const SHAPES: [Shape; 5] = [
    TRANSCRIPTS_ALONE,
    VARIANTS,
    VARIANTS_AND_SIFT,
    REGULATORY,
    VARIANTS_SIFT_AND_REGULATORY,
];

/// The SIFT cell the shaped caches hold: the change the `mis` record makes in
/// the coding transcript.
const SIFT_CELL: &str = "deleterious(0.01)";

/// The promoter of the shaped caches' regulatory tree, which the `mis` record
/// alone lies in.
const PROMOTER: &str = "SYNR0000001";

/// The warning of a run that asks for regulatory rows from a cache without
/// regulatory data.
const NO_REGULATORY_DATA: &str = "WARNING: --regulatory: this cache holds no regulatory data; \
     no RegulatoryFeature or MotifFeature row is written";

/// The warning of a run that asks for known variants from `cache`, which holds
/// no variation data, under the flags `asked` lists as the line names them.
fn no_variation_data(cache: &Path, asked: &str) -> String {
    format!(
        "WARNING: [cache] no variation data under {}: Existing_variation and the frequency \
         columns stay empty under {asked}; a cache from vep-cache-builder --auto or \
         --include-variations, or a conversion of Ensembl VEP's cache, holds it \
         (docs/cache-variation-data.md)",
        cache.display()
    )
}

/// Copies the directory tree at `from` to `to`.
fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A SIFT matrix in the cache's form, base64 of a gzip stream over `VEP` and
/// twenty little-endian 16-bit cells per residue, each the category in its top
/// two bits over the score times 1000 (`ProteinFunctionPredictionMatrix.pm`,
/// `pack 'v'`), a cell of all ones holding no prediction: every cell is empty
/// but residue 2 to valine, which reads [`SIFT_CELL`].
fn sift_matrix() -> String {
    use base64::Engine;
    use std::io::Write as _;

    const RESIDUES: usize = 32;
    // A residue's cells follow the amino acids' letters in alphabetical order.
    const VALINE: usize = 17;
    let mut raw = b"VEP".to_vec();
    raw.extend(std::iter::repeat_n(0xFFu8, 2 * 20 * RESIDUES));
    let cell = 3 + 2 * (20 + VALINE);
    raw[cell..cell + 2].copy_from_slice(&((1u16 << 14) | 10).to_le_bytes());
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&raw).unwrap();
    base64::engine::general_purpose::STANDARD.encode(gz.finish().unwrap())
}

/// A copy of the species cache at `<dir>/<shape.name>` holding what `shape`
/// names, its `info.json` declaring each kind of data only with its files, as
/// the converter writes it: the variation columns with the variation file, the
/// SIFT version and marker with the matrix, and the regulatory marker, schema
/// and cell types with the tree. The checksum manifest is written over the
/// result, so a run reads the copy with no manifest warning.
fn shaped_cache(dir: &Path, shape: Shape) -> PathBuf {
    let cache = dir.join(shape.name);
    copy_tree(&species_cache(), &cache);
    let read_json = |path: &Path| -> serde_json::Value {
        serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap()
    };
    let info_path = cache.join("info.json");
    let mut info = read_json(&info_path);
    if !shape.variants {
        fs::remove_dir_all(cache.join("variations")).unwrap();
        info.as_object_mut().unwrap().remove("variation_cols");
    }
    if shape.sift {
        info["sift"] = "b".into();
        info["source_versions"]["sift"] = "synthetic_sift_1".into();
        let shard = cache.join("transcripts/1/1-1000000.json");
        let mut transcripts = read_json(&shard);
        let coding = transcripts
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|transcript| transcript["stable_id"] == "SYNT0000001")
            .expect("the coding transcript");
        coding["variation_effect_feature_cache"]["protein_function_predictions"] = serde_json::json!({
            "sift": {"analysis": "sift", "predictions_data": sift_matrix()},
        });
        fs::write(&shard, transcripts.to_string()).unwrap();
    }
    if shape.regulatory {
        info["regulatory"] = true.into();
        info["regulatory_schema"] = vep_core::cache_info::REGULATORY_SCHEMA.into();
        info["cell_types"] = serde_json::json!(["synthetic_cell"]);
        fs::create_dir_all(cache.join("regulatory/1")).unwrap();
        fs::write(cache.join("regulatory/binding_matrices.json"), "{}").unwrap();
        let promoter = serde_json::json!([{
            "feature_type": "RegulatoryFeature",
            "stable_id": PROMOTER,
            "biotype": "promoter",
            "dbid": 1,
            "start": 990,
            "end": 1030,
            "activity": "A",
        }]);
        fs::write(
            cache.join("regulatory/1/1-1000000.json"),
            promoter.to_string(),
        )
        .unwrap();
    }
    fs::write(&info_path, info.to_string()).unwrap();
    vep_core::cache_manifest::CacheManifest::write(&cache).unwrap();
    cache
}

/// The values of `key` in the Extra column of a default-format output's rows,
/// in row order.
fn extra_values<'a>(text: &'a str, key: &str) -> Vec<&'a str> {
    data_lines(text)
        .into_iter()
        .filter_map(|row| row.split('\t').nth(13))
        .flat_map(|extra| extra.split(';'))
        .filter_map(|pair| pair.strip_prefix(key)?.strip_prefix('='))
        .collect()
}

/// The Feature column of the rows whose Feature_type is `feature_type`.
fn features_of_type<'a>(text: &'a str, feature_type: &str) -> Vec<&'a str> {
    data_lines(text)
        .into_iter()
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .filter(|columns| columns.get(5) == Some(&feature_type))
        .map(|columns| columns[4])
        .collect()
}

/// The lines of `err` that switch a prediction tool off.
fn disabling_notices(err: &str) -> Vec<&str> {
    err.lines()
        .filter(|line| line.starts_with("INFO: disabling "))
        .collect()
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
/// every record is kept (BaseCacheVariation.pm:341-345: "no frequency data so
/// we have to include it") and under `include` every record is dropped
/// (`:322-323`); `1KG_ALL` is the one population the check accepts on such a
/// cache (`:286-289`: the others must be declared columns), so `--filter_common`
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

/// `--regulatory` on a cache whose metadata declares no regulatory data builds
/// no regulatory source, as Ensembl VEP builds one only where the cache's
/// metadata carries the marker (CacheDir.pm:168-176, `if
/// $self->param('regulatory') and $info->{regulatory}`): no RegulatoryFeature
/// or MotifFeature row is written and every transcript row is the row of a run
/// without the flag, while the header still describes the flag's five motif
/// fields (Constants.pm:136). The run says so once for its five records,
/// through the warnings file and standard error, where Ensembl VEP annotates
/// silently. On the copy that holds a regulatory tree the same records get the
/// promoter's row and no warning, so the absence is the cache's.
#[test]
fn regulatory_on_a_cache_without_regulatory_data_warns_once_and_writes_no_row() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let bare = shaped_cache(dir.path(), TRANSCRIPTS_ALONE);

    let (plain, _, plain_warnings) = annotate_on(dir.path(), &bare, "plain", &[]);
    assert_eq!(plain_warnings, "", "no warning without the flag");
    let (text, err, warnings) = annotate_on(dir.path(), &bare, "regulatory", &["--regulatory"]);
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [NO_REGULATORY_DATA],
        "one warning, in the warnings file"
    );
    assert_eq!(
        err.lines()
            .filter(|line| *line == NO_REGULATORY_DATA)
            .count(),
        1,
        "the warning is echoed once on standard error:\n{err}"
    );
    assert!(!err.contains("ERROR:"), "{err}");
    for feature_type in ["RegulatoryFeature", "MotifFeature"] {
        assert!(
            features_of_type(&text, feature_type).is_empty(),
            "no {feature_type} row:\n{text}"
        );
    }
    assert_eq!(
        data_lines(&text),
        data_lines(&plain),
        "every transcript row is unchanged"
    );
    assert_eq!(data_lines(&text).len(), ROWS_WITHOUT_A_FILTER, "{text}");
    for field in [
        "MOTIF_NAME",
        "MOTIF_POS",
        "HIGH_INF_POS",
        "MOTIF_SCORE_CHANGE",
        "TRANSCRIPTION_FACTORS",
    ] {
        assert!(
            describes_field(&text, field),
            "the header describes {field}:\n{text}"
        );
    }

    let with_tree = shaped_cache(dir.path(), REGULATORY);
    let (text, _, warnings) = annotate_on(dir.path(), &with_tree, "tree", &["--regulatory"]);
    assert_eq!(warnings, "", "no warning on a cache holding the data");
    assert_eq!(
        features_of_type(&text, "RegulatoryFeature"),
        [PROMOTER],
        "{text}"
    );
    let transcript_rows: Vec<&str> = data_lines(&text)
        .into_iter()
        .filter(|row| !row.contains(PROMOTER))
        .collect();
    assert_eq!(
        transcript_rows,
        data_lines(&plain),
        "the regulatory row is added to the transcript rows"
    );
}

/// `check_sift_polyphen` reads each tool's version from the cache's metadata
/// (`AnnotationSource/Cache/Transcript.pm:148-176`), so what a flag does
/// follows what the cache holds. `--polyphen b` is refused on every shape, and
/// `--sift b` on each one without SIFT, with no input file to read: the refusal
/// comes before the input is opened. On a cache with SIFT `--sift b` writes the
/// prediction, and `--sift b --polyphen b` is still refused for PolyPhen.
/// Under `--everything` or `--safe` each tool the cache lacks is switched off
/// with its `INFO: disabling` line (`:161-166`) and the other tool's column is
/// written.
#[test]
fn sift_and_polyphen_are_refused_or_switched_off_by_what_the_cache_holds() {
    let dir = tempfile::tempdir().unwrap();
    let caches: Vec<(Shape, PathBuf)> = SHAPES
        .iter()
        .map(|shape| (*shape, shaped_cache(dir.path(), *shape)))
        .collect();
    let label = |shape: Shape, run: &str| format!("{}_{run}", shape.name);

    for (shape, cache) in &caches {
        assert_cache_refusal_on(
            dir.path(),
            cache,
            &label(*shape, "polyphen"),
            &["--polyphen", "b"],
            "PolyPhen not available",
        );
        let (flags, message): (&[&str], &str) = if shape.sift {
            (
                &["--sift", "b", "--polyphen", "b"],
                "PolyPhen not available",
            )
        } else {
            (&["--sift", "b"], "SIFT not available")
        };
        assert_cache_refusal_on(dir.path(), cache, &label(*shape, "sift"), flags, message);
    }

    write_input(dir.path());
    for (shape, cache) in &caches {
        let sift_cells: &[&str] = if shape.sift { &[SIFT_CELL] } else { &[] };
        if shape.sift {
            let (text, err, _) =
                annotate_on(dir.path(), cache, &label(*shape, "sift"), &["--sift", "b"]);
            assert!(disabling_notices(&err).is_empty(), "{shape:?}: {err}");
            assert!(describes_field(&text, "SIFT"), "{shape:?}:\n{text}");
            assert_eq!(
                extra_values(&text, "SIFT"),
                sift_cells,
                "{shape:?}:\n{text}"
            );
        }
        let mut notices = Vec::new();
        if !shape.sift {
            notices.push("INFO: disabling SIFT");
        }
        notices.push("INFO: disabling PolyPhen");
        for (run, flags) in [
            ("everything", &["--everything"][..]),
            ("safe", &["--sift", "b", "--polyphen", "b", "--safe"][..]),
        ] {
            let (text, err, _) = annotate_on(dir.path(), cache, &label(*shape, run), flags);
            assert_eq!(
                disabling_notices(&err),
                notices,
                "{shape:?} {flags:?}: {err}"
            );
            assert!(!err.contains("ERROR:"), "{shape:?} {flags:?}: {err}");
            assert_eq!(
                describes_field(&text, "SIFT"),
                shape.sift,
                "{shape:?} {flags:?}: a SIFT header line only where the cache holds SIFT:\n{text}"
            );
            assert_eq!(
                extra_values(&text, "SIFT"),
                sift_cells,
                "{shape:?} {flags:?}:\n{text}"
            );
            assert!(
                !describes_field(&text, "PolyPhen") && extra_values(&text, "PolyPhen").is_empty(),
                "{shape:?} {flags:?}: no PolyPhen column:\n{text}"
            );
        }
    }
}

/// The known-variant flags follow the cache's variation columns. Ensembl VEP
/// attaches a variation source only to a cache that declares them
/// (CacheDir.pm:179, `if($self->param('check_existing') &&
/// $info->{variation_cols})`), so on the cache of transcripts alone
/// `--check_existing` declares its columns (Constants.pm:127) and leaves every
/// cell empty, with one warning that Ensembl VEP does not print, and
/// `--af_gnomad` is accepted, its check sitting inside that branch
/// (CacheDir.pm:188-197); on the cache whose columns hold no `gnomADe` the flag is
/// refused (`:193-196`) and `--check_existing` reports the known variant. On
/// both, `--af`, `--max_af`, `--pubmed` and `--var_synonyms` declare their
/// columns (Constants.pm:118, 123, 128, 133) and leave them empty, no co-located
/// variant carrying a frequency, a citation or a synonym.
#[test]
fn known_variant_flags_follow_the_variation_columns() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let bare = shaped_cache(dir.path(), TRANSCRIPTS_ALONE);
    let with_variants = shaped_cache(dir.path(), VARIANTS);
    let existing_of = |text: &str| -> Vec<String> {
        let at = tab_columns(text)
            .iter()
            .position(|column| *column == "Existing_variation")
            .expect("the Existing_variation column");
        data_lines(text)
            .iter()
            .map(|row| row.split('\t').nth(at).unwrap().to_string())
            .collect()
    };
    let assert_warned_once = |err: &str, warnings: &str, what: &str| {
        let missing = no_variation_data(&bare, what);
        assert_eq!(
            warnings.lines().collect::<Vec<_>>(),
            [missing.as_str()],
            "{what}: one warning, in the warnings file"
        );
        assert_eq!(
            err.lines().filter(|line| *line == missing).count(),
            1,
            "{what}: the warning is echoed once on standard error:\n{err}"
        );
    };

    let known = ["Existing_variation", "CLIN_SIG", "SOMATIC", "PHENO"];
    let (text, err, warnings) = annotate_on(
        dir.path(),
        &bare,
        "bare_existing",
        &["--check_existing", "--tab"],
    );
    assert_columns_present_and_empty(&text, &known, "--check_existing without variation data");
    assert_warned_once(&err, &warnings, "--check_existing");

    let (text, _, warnings) = annotate_on(
        dir.path(),
        &with_variants,
        "variants_existing",
        &["--check_existing", "--tab"],
    );
    assert_eq!(warnings, "", "no warning on a cache holding variants");
    assert_columns_present_and_empty(&text, &known[1..], "--check_existing");
    assert_eq!(
        existing_of(&text),
        ["-", "-", "rs1", "rs1", "-", "-", "-", "-", "-"],
        "the known variant is reported on its record's rows:\n{text}"
    );

    assert_cache_refusal_on(
        dir.path(),
        &with_variants,
        "variants_gnomad",
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
    let (text, err, warnings) =
        annotate_on(dir.path(), &bare, "bare_gnomad", &["--af_gnomad", "--tab"]);
    assert_columns_present_and_empty(&text, &gnomade, "--af_gnomad without variation columns");
    assert_warned_once(&err, &warnings, "--af_gnomad");

    let flags = ["--af", "--max_af", "--pubmed", "--var_synonyms", "--tab"];
    let fields = ["AF", "MAX_AF", "MAX_AF_POPS", "PUBMED", "VAR_SYNONYMS"];
    let (text, err, warnings) = annotate_on(dir.path(), &bare, "bare_columns", &flags);
    assert_columns_present_and_empty(&text, &fields, "the column flags without variation data");
    assert_columns_present_and_empty(&text, &known, "the column flags without variation data");
    assert_warned_once(
        &err,
        &warnings,
        "--af, --max_af, --pubmed and --var_synonyms",
    );
    let (text, _, warnings) = annotate_on(dir.path(), &with_variants, "variants_columns", &flags);
    assert_eq!(warnings, "", "no warning on a cache holding variants");
    assert_columns_present_and_empty(&text, &fields, "the column flags");
    assert_eq!(
        existing_of(&text),
        ["-", "-", "rs1", "rs1", "-", "-", "-", "-", "-"],
        "{text}"
    );
}

/// `--everything` runs on every shape, and what it says names only what the
/// shape lacks: an `INFO: disabling` line for each prediction tool the cache
/// has no version of (`AnnotationSource/Cache/Transcript.pm:161-166`) and one
/// warning each for absent variation and regulatory data, where Ensembl VEP
/// builds neither source and says nothing (CacheDir.pm:176, 179). The header
/// describes SIFT and names its version only where the cache holds it, never
/// PolyPhen, and the rows carry the SIFT prediction, the known variant and the
/// promoter's row exactly where the shape holds each.
#[test]
fn everything_runs_on_every_shape() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    for shape in SHAPES {
        let cache = shaped_cache(dir.path(), shape);
        let (text, err, warnings) = annotate_on(dir.path(), &cache, shape.name, &["--everything"]);
        assert!(!err.contains("ERROR:"), "{shape:?}: {err}");

        let mut notices = Vec::new();
        if !shape.sift {
            notices.push("INFO: disabling SIFT");
        }
        notices.push("INFO: disabling PolyPhen");
        assert_eq!(disabling_notices(&err), notices, "{shape:?}: {err}");

        let mut expected_warnings = Vec::new();
        if !shape.variants {
            expected_warnings.push(no_variation_data(&cache, "--everything"));
        }
        if !shape.regulatory {
            expected_warnings.push(NO_REGULATORY_DATA.to_string());
        }
        assert_eq!(
            warnings.lines().collect::<Vec<_>>(),
            expected_warnings,
            "{shape:?}: one warning for each kind of data the cache lacks"
        );

        assert_eq!(
            describes_field(&text, "SIFT"),
            shape.sift,
            "{shape:?}:\n{text}"
        );
        assert_eq!(
            text.lines()
                .any(|line| line == "## sift version synthetic_sift_1"),
            shape.sift,
            "{shape:?}: the SIFT version line:\n{text}"
        );
        assert!(!describes_field(&text, "PolyPhen"), "{shape:?}:\n{text}");
        assert!(
            extra_values(&text, "PolyPhen").is_empty(),
            "{shape:?}:\n{text}"
        );

        let sift_cells: &[&str] = if shape.sift { &[SIFT_CELL] } else { &[] };
        assert_eq!(
            extra_values(&text, "SIFT"),
            sift_cells,
            "{shape:?}:\n{text}"
        );
        let promoters: &[&str] = if shape.regulatory { &[PROMOTER] } else { &[] };
        assert_eq!(
            features_of_type(&text, "RegulatoryFeature"),
            promoters,
            "{shape:?}:\n{text}"
        );
        assert_eq!(
            features_of_type(&text, "Transcript").len(),
            ROWS_WITHOUT_A_FILTER - 1,
            "{shape:?}: every transcript row:\n{text}"
        );
        let known: Vec<&str> = data_lines(&text)
            .iter()
            .filter(|row| row.starts_with("mis\t"))
            .map(|row| row.split('\t').nth(12).unwrap())
            .collect();
        let expected_known = if shape.variants { "rs1" } else { "-" };
        assert!(
            !known.is_empty() && known.iter().all(|cell| *cell == expected_known),
            "{shape:?}: Existing_variation of the co-located record is {expected_known}:\n{text}"
        );
    }
}
