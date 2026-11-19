// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The list-valued columns of a co-located variant in every output format,
//! held to what Ensembl VEP 116.2 writes for the same stored values.
//!
//! Each test stores one shape in one column of a known variant at `21:500` and
//! annotates `21 500 A G` against it. The expected values are the ones Ensembl
//! VEP 116.2 wrote for that stored value; every identifier is invented.
//!
//! The tests of a `pubmed` or `var_synonyms` value the cache cut at 1,024
//! bytes hold vep-rs to the stored list without its cut last item, where
//! Ensembl VEP 116.2 prints that item as an id.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use serde_json::{json, Value};

/// The columns of the test cache's variation file, in file order.
const COLUMNS: &[&str] = &[
    "variation_name",
    "failed",
    "somatic",
    "start",
    "end",
    "allele_string",
    "strand",
    "clin_sig",
    "phenotype_or_disease",
    "clin_sig_allele",
    "clin_sig_ref_allele",
    "clinical_impact",
    "pubmed",
    "var_synonyms",
    "AF",
    "AFR",
    "gnomADe",
    "gnomADe_AFR",
];

/// The flags every run passes: the ones writing the co-located columns the
/// tests read, and `--failed 1`, which keeps a known variant the cache flags
/// as failed.
const FLAGS: &[&str] = &[
    "--check_existing",
    "--af",
    "--af_1kg",
    "--af_gnomade",
    "--max_af",
    "--pubmed",
    "--var_synonyms",
    "--clinvar_somatic_classification",
    "--failed",
    "1",
];

/// The cells one known variant sets, by column name.
type Known<'a> = &'a [(&'a str, &'a str)];

/// A cache at `dir/cache` whose one variation file,
/// `variations/21/1-1000000_var.gz`, holds one known variant per entry of
/// `known` at `21:500`, named `rs9920000001` onwards, with the allele string
/// `A/G` unless the entry sets another; every column an entry does not set is
/// unset. An empty transcript shard lets the cache resolve, so the input's one
/// row is intergenic.
fn cache(dir: &Path, known: &[Known]) -> PathBuf {
    let cache = dir.join("cache");
    let shard = cache.join("transcripts").join("21").join("1-1000000.json");
    std::fs::create_dir_all(shard.parent().unwrap()).unwrap();
    std::fs::write(&shard, b"[]").unwrap();
    let region_file = cache.join("variations").join("21").join("1-1000000_var.gz");
    std::fs::create_dir_all(region_file.parent().unwrap()).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(&region_file).unwrap(),
        flate2::Compression::default(),
    );
    for (serial, columns) in known.iter().enumerate() {
        let name = format!("rs99200000{:02}", serial + 1);
        let mut cells: BTreeMap<&str, &str> = BTreeMap::from([
            ("variation_name", name.as_str()),
            ("start", "500"),
            ("allele_string", "A/G"),
        ]);
        for (column, value) in columns.iter() {
            assert!(COLUMNS.contains(column), "{column} is a cache column");
            assert!(!value.contains(' '), "a region file's cells hold no space");
            cells.insert(column, value);
        }
        let line: Vec<&str> = COLUMNS
            .iter()
            .map(|column| cells.get(column).copied().unwrap_or("."))
            .collect();
        writeln!(encoder, "{}", line.join(" ")).unwrap();
    }
    encoder.finish().unwrap();
    std::fs::write(
        cache.join("info.json"),
        json!({
            "assembly": "GRCh38",
            "source_versions": {"dbSNP": "156"},
            "variation_cols": COLUMNS,
        })
        .to_string(),
    )
    .unwrap();
    cache
}

/// Runs `vep` over `21 500 A G` against `cache` in `format` with [`FLAGS`] and
/// `extra`, and returns the output path; the run exits 0.
fn annotate(dir: &Path, cache: &Path, format: &str, extra: &[&str]) -> PathBuf {
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t500\t.\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.join(format!("out.{format}"));
    let output = Command::new(env!("CARGO_BIN_EXE_vep"))
        .stdin(Stdio::null())
        .env("HOME", dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args(["--output_format", format])
        .args(["--no_stats", "--force_overwrite", "--quiet"])
        .args(FLAGS)
        .args(extra)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "vep --output_format {format} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    out
}

/// The column names of a tab-separated output, read from its `#` column line.
fn column_names(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .find(|l| l.starts_with('#') && !l.starts_with("##"))
        .expect("a column line")
        .trim_start_matches('#')
        .split('\t')
}

/// The cells of the first data row of a tab-separated output by the names of
/// its `#` column line.
fn first_row(text: &str) -> BTreeMap<String, String> {
    let names = column_names(text);
    let row = text
        .lines()
        .find(|l| !l.starts_with('#'))
        .expect("a data row");
    names
        .zip(row.split('\t'))
        .map(|(name, cell)| (name.to_string(), cell.to_string()))
        .collect()
}

/// The elements of the list column `column` on the one row of the Parquet
/// output at `out`, `None` for a null cell.
#[cfg(feature = "parquet")]
fn parquet_list(out: &Path, column: &str) -> Option<Vec<String>> {
    use parquet::file::reader::FileReader;
    use parquet::record::Field;
    let parts = vep_io::output::parquet::part_files(out).unwrap();
    assert_eq!(parts.len(), 1, "one Parquet part");
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        std::fs::File::open(&parts[0]).unwrap(),
    )
    .unwrap();
    let mut rows = reader.get_row_iter(None).unwrap();
    let row = rows.next().expect("one Parquet row").unwrap();
    assert!(rows.next().is_none(), "one Parquet row");
    let (_, field) = row
        .get_column_iter()
        .find(|(name, _)| name.as_str() == column)
        .unwrap_or_else(|| panic!("the Parquet output has a {column} column"));
    match field {
        Field::Null => None,
        Field::ListInternal(list) => Some(
            list.elements()
                .iter()
                .map(|element| match element {
                    Field::Str(text) => text.clone(),
                    other => panic!("{column} holds a non-string element {other:?}"),
                })
                .collect(),
        ),
        other => panic!("{column} is not a list: {other:?}"),
    }
}

/// What each format writes for the one input record.
struct Written {
    /// The default format's `Extra` pairs.
    default: BTreeMap<String, String>,
    /// The `--tab` cells by column name.
    tab: BTreeMap<String, String>,
    /// The VCF `CSQ` subfields by name.
    vcf: BTreeMap<String, String>,
    /// The JSON record's `colocated_variants`, `Null` when it has none.
    json: Value,
    /// The Parquet row as `--tab` cells by column name.
    #[cfg(feature = "parquet")]
    parquet: BTreeMap<String, String>,
    /// The Parquet row's `PUBMED` and `CLIN_SIG` lists.
    #[cfg(feature = "parquet")]
    parquet_lists: BTreeMap<&'static str, Option<Vec<String>>>,
}

impl Written {
    /// Asserts every row format writes `value` for `field`, `None` being no
    /// value: an absent `Extra` key, a `-` cell, an empty subfield, a null
    /// Parquet cell. The VCF format joins a list with `&` where the others
    /// join it with `,`.
    fn assert_field(&self, field: &str, value: Option<&str>) {
        let vcf = value.map(|value| value.replace(',', "&"));
        let cell = |cells: &BTreeMap<String, String>, unset: &str| {
            cells.get(field).filter(|v| v.as_str() != unset).cloned()
        };
        assert_eq!(
            [
                self.default.get(field).cloned(),
                cell(&self.tab, "-"),
                cell(&self.vcf, ""),
            ],
            [value.map(str::to_string), value.map(str::to_string), vcf],
            "{field} in the default, tab and VCF formats"
        );
        // A typed Parquet cell prints as its type's canonical text: a DOUBLE
        // as its shortest decimal, where the tab format keeps `0.2500`.
        #[cfg(feature = "parquet")]
        assert_eq!(
            cell(&self.parquet, "-"),
            value.map(|value| {
                use vep_io::output::parquet::{canonical_tab_text, column_type};
                canonical_tab_text(column_type(field), value)
            }),
            "{field} in the Parquet format"
        );
    }

    /// Asserts the Parquet row's list column `field` holds `elements`.
    fn assert_parquet_list(&self, field: &str, elements: Option<&[&str]>) {
        #[cfg(feature = "parquet")]
        assert_eq!(
            self.parquet_lists[field],
            elements.map(|elements| elements.iter().map(|e| e.to_string()).collect::<Vec<_>>()),
            "the elements of the Parquet {field} list"
        );
        #[cfg(not(feature = "parquet"))]
        let _ = (field, elements);
    }
}

/// Annotates the input against a cache holding `known`, in every format, with
/// [`FLAGS`] and `extra`.
fn written(known: &[Known], extra: &[&str]) -> Written {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let cache = cache(dir, known);
    let read =
        |format: &str| std::fs::read_to_string(annotate(dir, &cache, format, extra)).unwrap();

    let default = first_row(&read("vep"))["Extra"]
        .split(';')
        .filter_map(|pair| pair.split_once('='))
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect();

    let tab_text = read("tab");
    let tab = first_row(&tab_text);

    let vcf_text = read("vcf");
    let names = vcf_text
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .and_then(|l| l.split('"').next())
        .expect("the CSQ header");
    let entry = vcf_text
        .lines()
        .find(|l| !l.starts_with('#'))
        .and_then(|l| l.split('\t').nth(7))
        .and_then(|info| info.split(';').find_map(|kv| kv.strip_prefix("CSQ=")))
        .expect("a CSQ entry");
    assert!(!entry.contains(','), "one CSQ entry: {entry}");
    let vcf = names
        .split('|')
        .zip(entry.split('|'))
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect();

    let record: Value = serde_json::from_str(read("json").lines().next().unwrap()).unwrap();
    let json = record
        .get("colocated_variants")
        .cloned()
        .unwrap_or(Value::Null);

    #[cfg(feature = "parquet")]
    let (parquet, parquet_lists) = {
        let out = annotate(dir, &cache, "parquet", extra);
        let rows = vep_io::output::parquet::tab_rows(&out).unwrap();
        assert_eq!(rows.len(), 1, "one Parquet row");
        let cells = column_names(&tab_text)
            .zip(rows[0].split('\t'))
            .map(|(name, cell)| (name.to_string(), cell.to_string()))
            .collect();
        let lists = ["PUBMED", "CLIN_SIG"]
            .into_iter()
            .map(|column| (column, parquet_list(&out, column)))
            .collect();
        (cells, lists)
    };

    Written {
        default,
        tab,
        vcf,
        json,
        #[cfg(feature = "parquet")]
        parquet,
        #[cfg(feature = "parquet")]
        parquet_lists,
    }
}

/// A `pubmed` value ending in `,` lists its ids and no empty one after them:
/// Perl's `split(',', $ex->{pubmed})` drops every trailing empty field
/// (OutputFactory.pm:1126 for the row formats, OutputFactory/JSON.pm:469 for
/// JSON). Ensembl VEP 116.2 writes `PUBMED=99000001,99000002` and
/// `"pubmed":[99000001,99000002]` for the stored `99000001,99000002,`.
#[test]
fn pubmed_ending_in_a_comma_lists_no_empty_id() {
    let out = written(&[&[("pubmed", "99000001,99000002,")]], &[]);
    out.assert_field("PUBMED", Some("99000001,99000002"));
    out.assert_parquet_list("PUBMED", Some(&["99000001", "99000002"]));
    assert_eq!(out.json[0]["pubmed"], json!([99000001, 99000002]));
}

/// A `clin_sig` value ending in `,` lists its terms and no empty one after
/// them (OutputFactory.pm:1125, OutputFactory/JSON.pm:469). Ensembl VEP 116.2
/// writes `CLIN_SIG=benign,likely_benign` and
/// `"clin_sig":["benign","likely_benign"]` for the stored
/// `benign,likely_benign,`.
#[test]
fn clin_sig_ending_in_a_comma_lists_no_empty_term() {
    let out = written(&[&[("clin_sig", "benign,likely_benign,")]], &[]);
    out.assert_field("CLIN_SIG", Some("benign,likely_benign"));
    out.assert_parquet_list("CLIN_SIG", Some(&["benign", "likely_benign"]));
    assert_eq!(out.json[0]["clin_sig"], json!(["benign", "likely_benign"]));
}

/// The lists of two co-located variants join with nothing between them when
/// the first one's stored value ends in `,` (OutputFactory.pm:1125-1126 push
/// each variant's fields onto one list). Ensembl VEP 116.2 writes
/// `PUBMED=99000001,99000002,99000003` and `CLIN_SIG=benign,pathogenic`.
#[test]
fn lists_of_two_variants_join_with_no_empty_element() {
    let out = written(
        &[
            &[("pubmed", "99000001,99000002,"), ("clin_sig", "benign,")],
            &[("pubmed", "99000003"), ("clin_sig", "pathogenic")],
        ],
        &[],
    );
    out.assert_field("PUBMED", Some("99000001,99000002,99000003"));
    out.assert_field("CLIN_SIG", Some("benign,pathogenic"));
    assert_eq!(out.json[0]["pubmed"], json!([99000001, 99000002]));
    assert_eq!(out.json[1]["pubmed"], json!([99000003]));
    assert_eq!(out.json[0]["clin_sig"], json!(["benign"]));
    assert_eq!(out.json[1]["clin_sig"], json!(["pathogenic"]));
}

/// A leading or a middle empty element is kept, as Perl's `split` keeps it.
/// Ensembl VEP 116.2 writes `PUBMED=,99000001` with `"pubmed":["",99000001]`
/// for the stored `,99000001`, and `CLIN_SIG=benign,,pathogenic` with
/// `"clin_sig":["benign","","pathogenic"]` for the stored
/// `benign,,pathogenic`.
#[test]
fn leading_and_middle_empty_elements_are_kept() {
    let out = written(
        &[&[("pubmed", ",99000001"), ("clin_sig", "benign,,pathogenic")]],
        &[],
    );
    out.assert_field("PUBMED", Some(",99000001"));
    out.assert_field("CLIN_SIG", Some("benign,,pathogenic"));
    out.assert_parquet_list("PUBMED", Some(&["", "99000001"]));
    assert_eq!(out.json[0]["pubmed"], json!(["", 99000001]));
    assert_eq!(out.json[0]["clin_sig"], json!(["benign", "", "pathogenic"]));
}

/// A value of separators alone is a list of no element: the row formats write
/// no field (the pushed list is empty and OutputFactory.pm:1136-1138 deletes
/// it), and JSON writes the empty list `split` returns (JSON.pm:469). Ensembl
/// VEP 116.2 writes no `PUBMED` or `CLIN_SIG` and `"pubmed":[]`,
/// `"clin_sig":[]` for the stored `,,` and `,`.
#[test]
fn separators_alone_are_a_list_of_no_element() {
    let out = written(&[&[("pubmed", ",,"), ("clin_sig", ",")]], &[]);
    out.assert_field("PUBMED", None);
    out.assert_field("CLIN_SIG", None);
    out.assert_parquet_list("PUBMED", None);
    assert_eq!(out.json[0]["pubmed"], json!([]));
    assert_eq!(out.json[0]["clin_sig"], json!([]));
}

/// A value Perl reads as false is no list in any format: `0` is neither pushed
/// (`if $ex->{pubmed}`, OutputFactory.pm:1125-1126) nor kept by the JSON
/// factory (JSON.pm:433-435). Ensembl VEP 116.2 writes no `PUBMED`, no
/// `CLIN_SIG` and neither JSON key for the stored `0`.
#[test]
fn a_value_of_zero_is_no_list() {
    let out = written(&[&[("pubmed", "0"), ("clin_sig", "0")]], &[]);
    out.assert_field("PUBMED", None);
    out.assert_field("CLIN_SIG", None);
    assert_eq!(out.json[0].get("pubmed"), None);
    assert_eq!(out.json[0].get("clin_sig"), None);
}

/// A list none of whose elements Perl reads as true prints no field in the
/// row formats (OutputFactory.pm:1136-1138 deletes it), while JSON, which
/// tests the whole stored value (JSON.pm:433-435), lists its elements.
/// Ensembl VEP 116.2 writes no `PUBMED` or `CLIN_SIG` with `"pubmed":[0,0]`
/// and `"clin_sig":["","",0]` for the stored `0,0` and `,,0`, and
/// `PUBMED=0,99000001` once one element is true.
#[test]
fn a_list_of_false_elements_prints_no_row_field() {
    let out = written(&[&[("pubmed", "0,0"), ("clin_sig", ",,0")]], &[]);
    out.assert_field("PUBMED", None);
    out.assert_field("CLIN_SIG", None);
    assert_eq!(out.json[0]["pubmed"], json!([0, 0]));
    assert_eq!(out.json[0]["clin_sig"], json!(["", "", 0]));

    let out = written(&[&[("pubmed", "0,99000001")]], &[]);
    out.assert_field("PUBMED", Some("0,99000001"));
    assert_eq!(out.json[0]["pubmed"], json!([0, 99000001]));
}

/// A frequency pair is read as two scalars, `my ($a, $f) = split(':', $pair)`
/// (OutputFactory.pm:1218): the value is the text between the pair's first and
/// second `:`, so a pair ending in `:` or carrying a third field prints its
/// frequency alone. Ensembl VEP 116.2 writes `AFR_AF=0.25`,
/// `gnomADe_AFR_AF=0.25` and `MAX_AF=0.25` for the stored `G:0.25:` and
/// `G:0.25:9`.
#[test]
fn frequency_value_ends_at_the_pairs_second_colon() {
    let out = written(&[&[("AFR", "G:0.25:")]], &[]);
    out.assert_field("AFR_AF", Some("0.25"));
    out.assert_field("MAX_AF", Some("0.25"));
    out.assert_field("MAX_AF_POPS", Some("AFR"));
    assert_eq!(out.json[0]["frequencies"], json!({"G": {"afr": 0.25}}));

    let out = written(&[&[("gnomADe_AFR", "G:0.25:9")]], &[]);
    out.assert_field("gnomADe_AFR_AF", Some("0.25"));
    out.assert_field("MAX_AF", Some("0.25"));
    assert_eq!(
        out.json[0]["frequencies"],
        json!({"G": {"gnomade_afr": 0.25}})
    );
}

/// A frequency column ending in `,` holds no empty pair after its last one
/// (`split(',', $ex->{$key})`, OutputFactory.pm:1217). Ensembl VEP 116.2
/// writes `AF=0.2500`, `AFR_AF=0.25` and `gnomADe_AF=0.125` for the stored
/// `G:0.25,`, `G:0.25,` and `G:0.125,`.
#[test]
fn frequency_column_ending_in_a_comma_prints_its_pairs() {
    let out = written(
        &[&[
            ("AF", "G:0.25,"),
            ("AFR", "G:0.25,"),
            ("gnomADe", "G:0.125,"),
        ]],
        &[],
    );
    out.assert_field("AF", Some("0.2500"));
    out.assert_field("AFR_AF", Some("0.25"));
    out.assert_field("gnomADe_AF", Some("0.125"));
    assert_eq!(
        out.json[0]["frequencies"],
        json!({"G": {"af": 0.25, "afr": 0.25, "gnomade": 0.125}})
    );
}

/// An empty pair after the last `,` names no allele. A known variant whose
/// allele string is `/G` has an empty first allele; with the stored `AF` of
/// `G:0.3,` that allele is the one no pair names, so the record's `AF` is
/// interpolated for it and the row of `G` prints its own (OutputFactory.pm:
/// 1214-1243). Ensembl VEP 116.2 writes `AF=0.3000` under
/// `--no_check_alleles`, which prints a known variant matched on its span.
#[test]
fn frequency_column_ending_in_a_comma_names_no_empty_allele() {
    let out = written(
        &[&[("allele_string", "/G"), ("AF", "G:0.3,")]],
        &["--no_check_alleles"],
    );
    out.assert_field("AF", Some("0.3000"));
    assert_eq!(out.json[0]["frequencies"], json!({"G": {"af": 0.3}}));
}

/// The frequency filter reads a column the same way (`split(',', $tmp)` then
/// `my ($a, $f) = split(':', $a_f)`, BaseCacheVariation.pm:234-235): a column
/// ending in `,` holds no empty pair, and a pair ending in `:` or carrying a
/// third field gives its frequency alone. Ensembl VEP 116.2 keeps the record
/// and writes `FREQS=1KG_AFR:G:0.25` for each of the stored `G:0.25,`,
/// `G:0.25:`, `G:0.25:9` and `G:0.25:,` when the filter includes a 1000
/// Genomes African frequency above 0.2.
#[test]
fn frequency_filter_reads_a_column_ending_in_a_separator() {
    let filter = [
        "--check_frequency",
        "--freq_pop",
        "1KG_AFR",
        "--freq_freq",
        "0.2",
        "--freq_gt_lt",
        "gt",
        "--freq_filter",
        "include",
    ];
    for stored in ["G:0.25,", "G:0.25:", "G:0.25:9", "G:0.25:,"] {
        let out = written(&[&[("AFR", stored)]], &filter);
        out.assert_field("FREQS", Some("1KG_AFR:G:0.25"));
    }
}

/// `var_synonyms` is written whole by the row formats (OutputFactory.pm:1090)
/// and split by the JSON factory on `--`, `::` and `,` (JSON.pm:452-456), each
/// split dropping its trailing empty fields. Ensembl VEP 116.2 writes
/// `VAR_SYNONYMS=ClinVar::RCV999000001,VCV999000002--` and
/// `"var_synonyms":{"ClinVar":["RCV999000001","VCV999000002"]}` for a stored
/// value ending in `--`, and one id for a group ending in `,`.
#[test]
fn var_synonyms_ending_in_a_separator_nest_no_empty_entry() {
    let out = written(
        &[&[("var_synonyms", "ClinVar::RCV999000001,VCV999000002--")]],
        &[],
    );
    assert_eq!(
        out.default["VAR_SYNONYMS"],
        "ClinVar::RCV999000001,VCV999000002--"
    );
    assert_eq!(
        out.tab["VAR_SYNONYMS"],
        "ClinVar::RCV999000001,VCV999000002--"
    );
    assert_eq!(
        out.vcf["VAR_SYNONYMS"],
        "ClinVar::RCV999000001&VCV999000002--"
    );
    assert_eq!(
        out.json[0]["var_synonyms"],
        json!({"ClinVar": ["RCV999000001", "VCV999000002"]})
    );

    let out = written(&[&[("var_synonyms", "ClinVar::RCV999000001,")]], &[]);
    assert_eq!(out.default["VAR_SYNONYMS"], "ClinVar::RCV999000001,");
    assert_eq!(
        out.json[0]["var_synonyms"],
        json!({"ClinVar": ["RCV999000001"]})
    );
}

/// `clin_sig_allele` is split on `;` for the row formats' `CLIN_SIG`
/// (OutputFactory.pm:1106) and written whole in JSON, which lists no such
/// field among those it splits (JSON.pm:123-128). Ensembl VEP 116.2 writes
/// `CLIN_SIG=benign` and `"clin_sig_allele":"G:benign;"` for the stored
/// `G:benign;`.
#[test]
fn clin_sig_allele_ending_in_a_semicolon_names_no_empty_entry() {
    let out = written(
        &[&[
            ("clin_sig", "benign"),
            ("clin_sig_allele", "G:benign;"),
            ("clin_sig_ref_allele", "A"),
        ]],
        &[],
    );
    out.assert_field("CLIN_SIG", Some("benign"));
    assert_eq!(out.json[0]["clin_sig_allele"], json!("G:benign;"));
    assert_eq!(out.json[0]["clin_sig"], json!(["benign"]));
}

/// `clinical_impact` is written whole in every format (OutputFactory.pm:1094;
/// JSON.pm:444-447 renames it before the list loop reaches it), the default
/// and VCF formats escaping its `;`. Ensembl VEP 116.2 writes
/// `CLINVAR_SOMATIC_CLASSIFICATION=Neoplasm(Tier_I_-_Strong)%3B` and
/// `"clinvar_somatic_classification":"Neoplasm(Tier_I_-_Strong);"` for the
/// stored `Neoplasm(Tier_I_-_Strong);`.
#[test]
fn clinical_impact_ending_in_a_semicolon_is_written_whole() {
    let out = written(&[&[("clinical_impact", "Neoplasm(Tier_I_-_Strong);")]], &[]);
    assert_eq!(
        out.default["CLINVAR_SOMATIC_CLASSIFICATION"],
        "Neoplasm(Tier_I_-_Strong)%3B"
    );
    assert_eq!(
        out.vcf["CLINVAR_SOMATIC_CLASSIFICATION"],
        "Neoplasm(Tier_I_-_Strong)%3B"
    );
    assert_eq!(
        out.json[0]["clinvar_somatic_classification"],
        json!("Neoplasm(Tier_I_-_Strong);")
    );
}

/// An allele string ending in `/` names no empty allele (`split('/', ...)`,
/// OutputFactory.pm:1191), so the record still matches the input's `G`, and
/// JSON writes the string whole. Ensembl VEP 116.2 writes `AF=0.2500` and
/// `"allele_string":"A/G/"` for the stored `A/G/` with `AF` `G:0.25`.
#[test]
fn allele_string_ending_in_a_slash_names_no_empty_allele() {
    let out = written(&[&[("allele_string", "A/G/"), ("AF", "G:0.25")]], &[]);
    out.assert_field("AF", Some("0.2500"));
    assert_eq!(out.json[0]["allele_string"], json!("A/G/"));
    assert_eq!(out.json[0]["frequencies"], json!({"G": {"af": 0.25}}));
}

/// The length in bytes at which the Ensembl cache dump ends a stored `pubmed`
/// or `var_synonyms` value.
const CAP: usize = 1024;

/// `count` invented eight-digit PubMed ids, `99000001` onwards, joined by `,`.
fn pubmed_ids(count: usize) -> String {
    let ids: Vec<String> = (1..=count).map(|n| (99_000_000 + n).to_string()).collect();
    ids.join(",")
}

/// A `var_synonyms` value of exactly [`CAP`] bytes ending in `tail`, with the
/// text before `tail`: a `UniProt` group whose one invented id is as long as
/// the length needs, then a `ClinVar` group of 70 invented ids.
fn synonyms_ending(tail: &str) -> (String, String) {
    let ids: Vec<String> = (1..=70)
        .map(|n| format!("RCV{}", 999_000_000 + n))
        .collect();
    let clinvar = format!("--ClinVar::{}", ids.join(","));
    let filler = CAP - tail.len() - "UniProt::VAR_".len() - clinvar.len();
    let body = format!("UniProt::VAR_{}{clinvar}", "9".repeat(filler));
    let stored = format!("{body}{tail}");
    assert_eq!(stored.len(), CAP, "a stored value at the cap");
    (stored, body)
}

/// The ids JSON lists under the source `source` of the first known variant.
fn json_ids<'a>(out: &'a Written, source: &str) -> &'a Vec<Value> {
    out.json[0]["var_synonyms"][source]
        .as_array()
        .unwrap_or_else(|| panic!("{source} ids in {}", out.json[0]["var_synonyms"]))
}

/// A `pubmed` value the cache cut inside an id is written without that id.
/// The stored value is 113 ids and the first seven digits of another, 1,024
/// bytes in all; Ensembl VEP 116.2 prints the seven digits as a 114th id
/// (`split(',', $ex->{pubmed})`, OutputFactory.pm:1126 and
/// OutputFactory/JSON.pm:469), and vep-rs writes the 113.
#[test]
fn pubmed_cut_inside_an_id_ends_before_it() {
    let whole = pubmed_ids(113);
    let stored = format!("{whole},9900011");
    assert_eq!(stored.len(), CAP);
    let out = written(&[&[("pubmed", stored.as_str())]], &[]);
    out.assert_field("PUBMED", Some(whole.as_str()));
    let ids: Vec<&str> = whole.split(',').collect();
    out.assert_parquet_list("PUBMED", Some(ids.as_slice()));
    let listed = out.json[0]["pubmed"].as_array().expect("a pubmed list");
    assert_eq!(listed.len(), 113, "the ids JSON lists");
    assert_eq!(listed[112], json!(99000113));
}

/// A `pubmed` value of 1,024 bytes whose last item is as long as a whole id
/// loses it too: nothing in the stored value tells such an item from the
/// first eight digits of a longer id.
#[test]
fn pubmed_at_the_cap_loses_a_last_id_of_whole_length() {
    let kept = format!("{},9900113", pubmed_ids(112));
    let stored = format!("{kept},99000114");
    assert_eq!(stored.len(), CAP);
    let out = written(&[&[("pubmed", stored.as_str())]], &[]);
    out.assert_field("PUBMED", Some(kept.as_str()));
    let listed = out.json[0]["pubmed"].as_array().expect("a pubmed list");
    assert_eq!(listed.len(), 113, "the ids JSON lists");
    assert_eq!(listed[112], json!(9900113));
}

/// A `pubmed` value of 1,024 bytes with no `,` is one cut item: no format
/// writes the field.
#[test]
fn pubmed_of_one_cut_item_writes_no_field() {
    let out = written(&[&[("pubmed", "9".repeat(CAP).as_str())]], &[]);
    out.assert_field("PUBMED", None);
    out.assert_parquet_list("PUBMED", None);
    assert_eq!(out.json[0].get("pubmed"), None);
}

/// The lists of two known variants join without the id the cache cut from
/// the first one's `pubmed` and `var_synonyms` values, and each JSON entry
/// lists its own variant's whole items.
#[test]
fn lists_of_two_variants_join_without_the_first_ones_cut_item() {
    let ids = pubmed_ids(113);
    let pubmed = format!("{ids},9900011");
    let (synonyms, body) = synonyms_ending(",RCV999");
    let out = written(
        &[
            &[
                ("pubmed", pubmed.as_str()),
                ("var_synonyms", synonyms.as_str()),
            ],
            &[
                ("pubmed", "99000201"),
                ("var_synonyms", "UniProt::VAR_999002"),
            ],
        ],
        &[],
    );
    out.assert_field("PUBMED", Some(format!("{ids},99000201").as_str()));
    out.assert_field(
        "VAR_SYNONYMS",
        Some(format!("{body}--UniProt::VAR_999002").as_str()),
    );
    assert_eq!(out.json[0]["pubmed"].as_array().map(Vec::len), Some(113));
    assert_eq!(out.json[1]["pubmed"], json!([99000201]));
    assert_eq!(json_ids(&out, "ClinVar").len(), 70);
    assert_eq!(
        out.json[1]["var_synonyms"],
        json!({"UniProt": ["VAR_999002"]})
    );
}

/// Control: a `pubmed` value of 1,024 bytes ending in `,` ends after a
/// whole id and is written as stored, the split dropping the empty field
/// after the last `,`.
#[test]
fn pubmed_at_the_cap_ending_in_a_comma_is_written_whole() {
    let whole = format!("{},990011", pubmed_ids(113));
    let stored = format!("{whole},");
    assert_eq!(stored.len(), CAP);
    let out = written(&[&[("pubmed", stored.as_str())]], &[]);
    out.assert_field("PUBMED", Some(whole.as_str()));
    assert_eq!(out.json[0]["pubmed"].as_array().map(Vec::len), Some(114));
}

/// Control: a `pubmed` or `var_synonyms` value one byte under the cap, or one
/// byte over it, is written as stored, its last item included.
#[test]
fn a_value_beside_the_cap_is_written_whole() {
    let under = format!("{},990011", pubmed_ids(113));
    assert_eq!(under.len(), CAP - 1);
    let (at_the_cap, _) = synonyms_ending(",RCV999");
    let synonyms_under = at_the_cap.replacen('9', "", 1);
    assert_eq!(synonyms_under.len(), CAP - 1);
    let out = written(
        &[&[
            ("pubmed", under.as_str()),
            ("var_synonyms", synonyms_under.as_str()),
        ]],
        &[],
    );
    out.assert_field("PUBMED", Some(under.as_str()));
    out.assert_field("VAR_SYNONYMS", Some(synonyms_under.as_str()));
    assert_eq!(json_ids(&out, "ClinVar").len(), 71);

    let over = format!("{},99000114", pubmed_ids(113));
    assert_eq!(over.len(), CAP + 1);
    let synonyms_over = format!("{at_the_cap}0");
    let out = written(
        &[&[
            ("pubmed", over.as_str()),
            ("var_synonyms", synonyms_over.as_str()),
        ]],
        &[],
    );
    out.assert_field("PUBMED", Some(over.as_str()));
    out.assert_field("VAR_SYNONYMS", Some(synonyms_over.as_str()));
}

/// Control: a short list is written as stored.
#[test]
fn a_short_list_is_written_whole() {
    let out = written(
        &[&[
            ("pubmed", "99000001,99000002,99000003"),
            ("var_synonyms", "ClinVar::RCV999000001--UniProt::VAR_999001"),
        ]],
        &[],
    );
    out.assert_field("PUBMED", Some("99000001,99000002,99000003"));
    out.assert_field(
        "VAR_SYNONYMS",
        Some("ClinVar::RCV999000001--UniProt::VAR_999001"),
    );
    assert_eq!(
        out.json[0]["var_synonyms"],
        json!({"ClinVar": ["RCV999000001"], "UniProt": ["VAR_999001"]})
    );
}

/// A `var_synonyms` value the cache cut inside an id is written without that
/// id and the `,` before it. Ensembl VEP 116.2 prints the stored value whole
/// in the row formats (OutputFactory.pm:1090) and lists the fragment as an id
/// in JSON (OutputFactory/JSON.pm:452-456).
#[test]
fn synonyms_cut_inside_an_id_end_before_it() {
    let (stored, body) = synonyms_ending(",RCV999");
    let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
    out.assert_field("VAR_SYNONYMS", Some(body.as_str()));
    let ids = json_ids(&out, "ClinVar");
    assert_eq!(ids.len(), 70, "the ClinVar ids JSON lists");
    assert_eq!(ids[69], json!("RCV999000070"));
}

/// A `var_synonyms` value of 1,024 bytes whose last id is as long as a whole
/// one loses it too.
#[test]
fn synonyms_at_the_cap_lose_a_last_id_of_whole_length() {
    let (stored, body) = synonyms_ending(",RCV999000999");
    let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
    out.assert_field("VAR_SYNONYMS", Some(body.as_str()));
    assert_eq!(json_ids(&out, "ClinVar").len(), 70);
}

/// What each format writes for a `var_synonyms` value of 1,024 bytes cut
/// inside the only id of its last group, `--OMIM::6059`, with the text
/// before that group. A source whose only id is the cut one goes with its
/// id, so no format writes a source name without an id: Ensembl VEP 116.2
/// prints none from a whole list, and a reader of the field expects the pair.
fn only_id_of_its_source_cut() -> (Written, String) {
    let (stored, body) = synonyms_ending("--OMIM::6059");
    let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
    (out, body)
}

/// The default format writes the groups before the cut one and no bare
/// source name after them.
#[test]
fn a_source_whose_only_id_is_cut_goes_in_the_default_format() {
    let (out, body) = only_id_of_its_source_cut();
    assert_eq!(out.default.get("VAR_SYNONYMS"), Some(&body));
}

/// The `--tab` cell holds the groups before the cut one.
#[test]
fn a_source_whose_only_id_is_cut_goes_in_the_tab_format() {
    let (out, body) = only_id_of_its_source_cut();
    assert_eq!(out.tab.get("VAR_SYNONYMS"), Some(&body));
}

/// The VCF `CSQ` subfield holds the groups before the cut one, its ids
/// joined by `&`.
#[test]
fn a_source_whose_only_id_is_cut_goes_in_the_vcf_format() {
    let (out, body) = only_id_of_its_source_cut();
    assert_eq!(out.vcf.get("VAR_SYNONYMS"), Some(&body.replace(',', "&")));
}

/// The JSON object has no key for the cut group's source, and so no empty
/// array under it; the sources before it keep their ids.
#[test]
fn a_source_whose_only_id_is_cut_goes_in_json() {
    let (out, _) = only_id_of_its_source_cut();
    let sources = out.json[0]["var_synonyms"]
        .as_object()
        .expect("a var_synonyms object");
    assert_eq!(sources.get("OMIM"), None, "no source without an id");
    assert_eq!(sources.len(), 2, "the two sources before the cut one");
    assert_eq!(json_ids(&out, "ClinVar").len(), 70);
}

/// The Parquet cell follows the `--tab` text: the groups before the cut one.
#[cfg(feature = "parquet")]
#[test]
fn a_source_whose_only_id_is_cut_goes_in_parquet() {
    use vep_io::output::parquet::{canonical_tab_text, column_type};
    let (out, body) = only_id_of_its_source_cut();
    assert_eq!(
        out.parquet.get("VAR_SYNONYMS"),
        Some(&canonical_tab_text(column_type("VAR_SYNONYMS"), &body))
    );
}

/// A last group with no `::` is a source name the cut may have ended, in its
/// middle or after its last letter: it goes with the `--` before it.
#[test]
fn synonyms_cut_inside_a_source_name_end_before_it() {
    for tail in ["--OMI", "--OMIM"] {
        let (stored, body) = synonyms_ending(tail);
        let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
        out.assert_field("VAR_SYNONYMS", Some(body.as_str()));
        assert_eq!(out.json[0]["var_synonyms"].get("OMIM"), None, "{tail}");
        assert_eq!(out.json[0]["var_synonyms"].get("OMI"), None, "{tail}");
    }
}

/// A value ending on one character of `--` or of `::` does not end in a
/// separator: its last item goes, an id with the `-` after it or a source
/// name with the `:` after it.
#[test]
fn synonyms_cut_on_half_a_separator_lose_their_last_item() {
    let (stored, body) = synonyms_ending(",RCV999000999-");
    let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
    out.assert_field("VAR_SYNONYMS", Some(body.as_str()));
    assert_eq!(json_ids(&out, "ClinVar").len(), 70);

    let (stored, body) = synonyms_ending("--OMIM:");
    let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
    out.assert_field("VAR_SYNONYMS", Some(body.as_str()));
    assert_eq!(out.json[0]["var_synonyms"].get("OMIM:"), None);
}

/// A `var_synonyms` value of 1,024 bytes that holds no whole id writes no
/// field: one item with no separator, and one group whose only id is cut.
#[test]
fn synonyms_of_one_cut_item_write_no_field() {
    let no_separator = format!("ClinVar{}", "9".repeat(CAP - "ClinVar".len()));
    let only_id = format!("ClinVar::RCV{}", "9".repeat(CAP - "ClinVar::RCV".len()));
    for stored in [no_separator, only_id] {
        assert_eq!(stored.len(), CAP);
        let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
        out.assert_field("VAR_SYNONYMS", None);
        assert_eq!(out.json[0].get("var_synonyms"), None);
    }
}

/// Control: a `var_synonyms` value of 1,024 bytes ending in `,`, `::` or `--`
/// ends after a separator and is written as stored: whole in the row
/// formats, and in JSON with no empty id and no id under a source name that
/// `::` ends.
#[test]
fn synonyms_at_the_cap_ending_in_a_separator_are_written_whole() {
    for (tail, omim) in [(",", None), ("--OMIM::", Some(json!([]))), ("--", None)] {
        let (stored, _) = synonyms_ending(tail);
        let out = written(&[&[("var_synonyms", stored.as_str())]], &[]);
        out.assert_field("VAR_SYNONYMS", Some(stored.as_str()));
        assert_eq!(json_ids(&out, "ClinVar").len(), 70, "{tail}");
        assert_eq!(
            out.json[0]["var_synonyms"].get("OMIM"),
            omim.as_ref(),
            "{tail}"
        );
    }
}
