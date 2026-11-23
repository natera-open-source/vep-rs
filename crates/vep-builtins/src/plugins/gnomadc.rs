// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! gnomADc plugin: gnomAD coverage statistics, as Ensembl VEP's `gnomADc.pm`
//! (VEP_plugins release/116) writes them.
//!
//! The coverage file is read by position alone, with no allele matching, and
//! the result is stored at the variant level. The output keys are
//! `gnomAD_<key>_cov`, chosen by the header's column count (`gnomADc.pm:112,
//! 136-148, 187-198`): a fourteen-column file, the GRCh38 v3 and v4 layout
//! once its `locus` column is split, prints `mean`, `median_approx`, `total_DP`
//! and the nine depth fractions; any other count prints `mean`, `median` and
//! the nine fractions, the GRCh37 v2.1 layout. The values are read by column
//! position from the third column on, not by the header's names (`:185, 202`).
//!
//! A variant spanning several positions (a deletion) averages every record of
//! its span, each value printed with four decimals (`:170-175`), so a
//! single-position variant prints its record's values rounded to four
//! decimals, `1.1177e-01` as `0.1118`. An insertion, whose end lies before
//! its start, reads the record at its start (`:161`). A cell is read as Perl's
//! addition reads a string: its leading number, zero when it opens with none.
//!
//! **Usage:** `--plugin gnomADc,/path/to/gnomad.coverage.tsv.gz`

use std::path::Path;

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::perl_format;
use crate::tabix::{table_contig, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;

const NAME: &str = "gnomADc";

/// The prefix of every output key (`gnomADc.pm:99`).
const PREFIX: &str = "gnomAD";

/// The value columns of a fourteen-column file, in column order (`:191`).
const FOURTEEN_COLUMN_KEYS: &[&str] = &[
    "mean",
    "median_approx",
    "total_DP",
    "1x",
    "5x",
    "10x",
    "15x",
    "20x",
    "25x",
    "30x",
    "50x",
    "100x",
];

/// The value columns of a file with any other column count (`:197`).
const OTHER_COLUMN_KEYS: &[&str] = &[
    "mean", "median", "1x", "5x", "10x", "15x", "20x", "25x", "30x", "50x", "100x",
];

/// The number of value columns in the longer key list, the size of a variant's
/// sums.
const MAX_VALUE_COLUMNS: usize = if FOURTEEN_COLUMN_KEYS.len() > OTHER_COLUMN_KEYS.len() {
    FOURTEEN_COLUMN_KEYS.len()
} else {
    OTHER_COLUMN_KEYS.len()
};

/// The header's description of a key: the depth fractions, a depth followed by
/// `x` (`:146-148`), else the depth statistics (`:137-144`).
fn description(key: &str) -> String {
    let is_fraction = key
        .strip_suffix('x')
        .is_some_and(|depth| !depth.is_empty() && depth.bytes().all(|b| b.is_ascii_digit()));
    if is_fraction {
        format!("Fraction of samples at {key} coverage")
    } else {
        format!("{key} coverage")
    }
}

/// The positions a variant reads: its own span, an insertion's end, which lies
/// before its start, raised to it (`:161`). An allele of a `--minimal` split
/// record reads its own trimmed span (the documented class
/// `minimal_plugin_sees_the_rejoined_record`).
fn coverage_span(variant: &InputVariant) -> (u64, u64) {
    (variant.start, variant.end.max(variant.start))
}

/// A cell as Perl's addition reads it: its longest leading number (`12abc` is
/// 12), or zero for text opening with none. A well-formed cell parses whole at
/// the first attempt.
fn numeric(cell: &str) -> f64 {
    let text = cell.trim();
    (0..=text.len())
        .rev()
        .filter(|&end| text.is_char_boundary(end))
        .find_map(|end| text[..end].parse::<f64>().ok())
        .unwrap_or(0.0)
}

/// gnomAD coverage statistics plugin.
pub struct GnomadcPlugin {
    /// Annotation store, initialized after `init()` is called.
    store: Option<Box<dyn AnnotationStore>>,
    /// The output key of each value column, from the third column on.
    keys: Vec<String>,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for GnomadcPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl GnomadcPlugin {
    /// Create an uninitialized gnomADc plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            keys: Vec::new(),
            assembly: None,
        }
    }

    /// The output keys of a file whose header has `columns` whitespace-separated
    /// names (`:112, 187-198`).
    fn keys_for(columns: usize) -> Vec<String> {
        let names = if columns == 14 {
            FOURTEEN_COLUMN_KEYS
        } else {
            OTHER_COLUMN_KEYS
        };
        names
            .iter()
            .map(|name| format!("{PREFIX}_{name}_cov"))
            .collect()
    }

    /// Writes each value column's mean over `records`, in key order, into the
    /// variant's fields with four decimals (`:167-175`). A record short of a
    /// column adds nothing to that column's sum and still counts toward its
    /// divisor, as Perl's undefined value does. The sums live on the stack and
    /// the map is sized once, so a variant costs the key and value strings
    /// alone beyond the lookup.
    fn write_means(&self, records: &[&TabixRecord], variant: &mut InputVariant) {
        let mut sums = [0.0f64; MAX_VALUE_COLUMNS];
        let sums = &mut sums[..self.keys.len()];
        for record in records {
            for (sum, cell) in sums.iter_mut().zip(record.columns.iter().skip(2)) {
                *sum += numeric(cell);
            }
        }
        let count = records.len() as f64;
        variant.plugin_data.reserve(self.keys.len());
        for (key, sum) in self.keys.iter().zip(sums.iter()) {
            variant.set_plugin_value(
                key.clone(),
                perl_format::fixed(sum / count, 4),
                false,
                self.feature_types(),
            );
        }
    }
}

impl crate::traits::sealed::Sealed for GnomadcPlugin {}

impl BuiltinPlugin for GnomadcPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`: every row,
    /// the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.keys
            .iter()
            .map(|key| {
                let name = &key[PREFIX.len() + 1..key.len() - "_cov".len()];
                (key.clone(), description(name))
            })
            .collect()
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        if let Some(store) = self.store.as_mut() {
            store.set_reference_fasta(fasta);
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &[], params)?;
        let file = params
            .positional()
            .first()
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: the coverage file is the first parameter (--plugin {NAME},<file>)"
                ))
            })?;
        let path = Path::new(file);
        if !path.is_file() {
            return Err(PluginError::Init(format!(
                "{NAME}: coverage file {file} not found"
            )));
        }

        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        };
        // A row keeps the whitespace before its line end, which Ensembl VEP's
        // `split /\t/` leaves in the last cell (`:185`), so the table reads as a
        // store built from it does; the addition of `:171` reads the number.
        let store = open_best_store_with(
            path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;

        // The count is of whitespace-separated names, as Perl's `split` counts
        // them (`:108`).
        let columns = store
            .header()
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: coverage file {file} has no header line naming its columns"
                ))
            })?
            .iter()
            .flat_map(|name| name.split_ascii_whitespace())
            .count();
        self.keys = Self::keys_for(columns);

        debug!(file, columns, keys = ?self.keys, "gnomADc plugin initialized");

        self.store = Some(store);
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        if variants.is_empty() {
            return Ok(None);
        }

        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|variant| self.annotates_variant(variant))
            .map(|variant| {
                let (start, end) = coverage_span(variant);
                (table_contig(variant, self.assembly), start, end)
            })
            .collect();

        let batch_result = store.query_batch(&regions)?;
        Ok(Some(PrefetchData::new(batch_result)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch_result: BatchQueryResult = data.downcast()?;

        if batch_result.is_empty() {
            return Ok(());
        }

        for variant in variants
            .iter_mut()
            .filter(|variant| self.annotates_variant(variant))
        {
            let (start, end) = coverage_span(variant);
            let records = batch_result.lookup(table_contig(variant, self.assembly), start, end);
            if records.is_empty() {
                continue;
            }
            self.write_means(&records, variant);
        }

        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        if let Some(data) = self.prefetch(variants)? {
            self.annotate(variants, data)?;
        }
        Ok(())
    }

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Variant-level plugin: annotation happens in `run_batch`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, HashMap};

    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

    fn record(start: u64, cells: &[&str]) -> TabixRecord {
        let mut columns = vec!["21".to_string(), start.to_string()];
        columns.extend(cells.iter().map(|cell| cell.to_string()));
        TabixRecord {
            chr: "21".to_string(),
            start,
            end: start,
            columns,
            header: None,
        }
    }

    fn batch(records: Vec<TabixRecord>) -> BatchQueryResult {
        let mut by_start: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
        for record in records {
            by_start.entry(record.start).or_default().push(record);
        }
        BatchQueryResult::from_records(HashMap::from([("21".to_string(), by_start)]))
    }

    fn variant(start: u64, end: u64) -> InputVariant {
        let mut variant = InputVariant::default();
        variant.chr = "21".to_string();
        variant.original_chr = "21".to_string();
        variant.start = start;
        variant.end = end;
        variant
    }

    fn plugin(columns: usize) -> GnomadcPlugin {
        GnomadcPlugin {
            store: None,
            keys: GnomadcPlugin::keys_for(columns),
            assembly: None,
        }
    }

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = GnomadcPlugin::new();
        assert_eq!(plugin.name(), "gnomADc");
        assert!(plugin.store.is_none());
        assert!(plugin.keys.is_empty());
    }

    #[test]
    fn init_requires_the_file_as_the_first_parameter() {
        let mut plugin = GnomadcPlugin::new();
        let error = plugin.init(&[]).unwrap_err().to_string();
        assert!(error.contains("gnomADc: the coverage file is the first parameter"));
        let error = plugin
            .init(&["file=coverage.tsv.gz".to_string()])
            .unwrap_err()
            .to_string();
        assert!(error.contains("gnomADc: unknown parameter 'file'"));
        let error = plugin
            .init(&["/no/such/coverage.tsv.gz".to_string()])
            .unwrap_err()
            .to_string();
        assert!(error.contains("gnomADc: coverage file /no/such/coverage.tsv.gz not found"));
    }

    /// `gnomADc.pm:136-140, 146-148, 187-192`.
    #[test]
    fn header_keys_follow_ensembl_names_for_14_columns() {
        let headers = plugin(14).header_info();
        let keys: Vec<&str> = headers.iter().map(|(key, _)| key.as_str()).collect();
        assert_eq!(
            keys,
            [
                "gnomAD_mean_cov",
                "gnomAD_median_approx_cov",
                "gnomAD_total_DP_cov",
                "gnomAD_1x_cov",
                "gnomAD_5x_cov",
                "gnomAD_10x_cov",
                "gnomAD_15x_cov",
                "gnomAD_20x_cov",
                "gnomAD_25x_cov",
                "gnomAD_30x_cov",
                "gnomAD_50x_cov",
                "gnomAD_100x_cov",
            ]
        );
        assert_eq!(headers[1].1, "median_approx coverage");
        assert_eq!(headers[2].1, "total_DP coverage");
        assert_eq!(headers[3].1, "Fraction of samples at 1x coverage");
    }

    /// `gnomADc.pm:141-148, 193-198`: thirteen columns, and any count but
    /// fourteen, print `median` and no `total_DP`.
    #[test]
    fn header_keys_follow_ensembl_names_for_13_columns() {
        let headers = plugin(13).header_info();
        let keys: Vec<&str> = headers.iter().map(|(key, _)| key.as_str()).collect();
        assert_eq!(
            keys,
            [
                "gnomAD_mean_cov",
                "gnomAD_median_cov",
                "gnomAD_1x_cov",
                "gnomAD_5x_cov",
                "gnomAD_10x_cov",
                "gnomAD_15x_cov",
                "gnomAD_20x_cov",
                "gnomAD_25x_cov",
                "gnomAD_30x_cov",
                "gnomAD_50x_cov",
                "gnomAD_100x_cov",
            ]
        );
        assert_eq!(headers[0].1, "mean coverage");
        assert_eq!(headers[1].1, "median coverage");
        assert_eq!(headers[10].1, "Fraction of samples at 100x coverage");
        assert_eq!(plugin(15).header_info().len(), 11);
    }

    /// `gnomADc.pm:170-175`: a two-position deletion averages the two records.
    #[test]
    fn multi_position_variant_reports_four_decimal_mean() {
        let plugin = plugin(13);
        let mut variants = vec![variant(100, 101)];
        let data = PrefetchData::new(batch(vec![
            record(100, &["30.1", "30", "1.0000", "0.9", "0.5"]),
            record(101, &["31.3", "32", "1.0000", "0.8", "0.5"]),
        ]));
        plugin.annotate(&mut variants, data).unwrap();
        let fields = &variants[0].plugin_data;
        assert_eq!(fields["gnomAD_mean_cov"], "30.7000");
        assert_eq!(fields["gnomAD_median_cov"], "31.0000");
        assert_eq!(fields["gnomAD_1x_cov"], "1.0000");
        assert_eq!(fields["gnomAD_5x_cov"], "0.8500");
        assert_eq!(fields["gnomAD_10x_cov"], "0.5000");
        // A column both records lack sums to nothing over the two records.
        assert_eq!(fields["gnomAD_100x_cov"], "0.0000");
        assert_eq!(fields.len(), 11);
    }

    /// `gnomADc.pm:161`: an insertion (`start` one past `end`) reads the record
    /// at its start alone, not the two records its reversed span would cover.
    #[test]
    fn insertion_clamps_end_to_start() {
        assert_eq!(coverage_span(&variant(101, 100)), (101, 101));
        assert_eq!(coverage_span(&variant(100, 102)), (100, 102));
        let plugin = plugin(14);
        let mut variants = vec![variant(101, 100)];
        let data = PrefetchData::new(batch(vec![
            record(100, &["10.0", "10", "100"]),
            record(101, &["20.0", "20", "200"]),
        ]));
        plugin.annotate(&mut variants, data).unwrap();
        let fields = &variants[0].plugin_data;
        assert_eq!(fields["gnomAD_mean_cov"], "20.0000");
        assert_eq!(fields["gnomAD_median_approx_cov"], "20.0000");
        assert_eq!(fields["gnomAD_total_DP_cov"], "200.0000");
    }

    /// An allele of a `--minimal` split record averages its own trimmed
    /// positions, the ones its row's alleles describe, as the same allele does
    /// when it is the record's only alternate; the record's span, the one its
    /// `Location` prints, is not read. An insertion allele reads the base after
    /// its own insertion point.
    #[test]
    fn a_split_allele_reads_its_own_span() {
        let mut split = variant(101, 102);
        split.original_start = Some(100);
        split.original_end = Some(102);
        assert_eq!(coverage_span(&split), (101, 102));
        let mut insertion = variant(101, 100);
        insertion.original_start = Some(100);
        insertion.original_end = Some(100);
        assert_eq!(coverage_span(&insertion), (101, 101));
        let plugin = plugin(14);
        let mut variants = vec![split, insertion, variant(101, 102)];
        let data = PrefetchData::new(batch(vec![
            record(100, &["8.2345e+00", "6", "6018955"]),
            record(101, &["8.2555e+00", "7", "6034327"]),
            record(102, &["8.5104e+00", "7", "6220647"]),
        ]));
        plugin.annotate(&mut variants, data).unwrap();
        let split = &variants[0].plugin_data;
        assert_eq!(split["gnomAD_mean_cov"], "8.3830");
        assert_eq!(split["gnomAD_median_approx_cov"], "7.0000");
        assert_eq!(split["gnomAD_total_DP_cov"], "6127487.0000");
        let insertion = &variants[1].plugin_data;
        assert_eq!(insertion["gnomAD_mean_cov"], "8.2555");
        assert_eq!(insertion["gnomAD_total_DP_cov"], "6034327.0000");
        let own = &variants[2].plugin_data;
        assert_eq!(own["gnomAD_mean_cov"], "8.3830");
        assert_eq!(own["gnomAD_median_approx_cov"], "7.0000");
        assert_eq!(own["gnomAD_total_DP_cov"], "6127487.0000");
    }

    /// `gnomADc.pm:175`: every value, an integer depth included, prints with
    /// four decimals, and the file's exponent forms are read as numbers.
    #[test]
    fn integer_total_dp_prints_four_decimals() {
        let plugin = plugin(14);
        let mut variants = vec![variant(5032011, 5032011)];
        let data = PrefetchData::new(batch(vec![record(
            5032011,
            &[
                "1.0757e+01",
                "8",
                "12345",
                "8.3068e-01",
                "6.9381e-01",
                "4.0174e-01",
                "2.8219e-01",
                "2.2035e-01",
                "1.5021e-01",
                "6.0152e-02",
                "5.3355e-05",
                "8.2085e-06",
            ],
        )]));
        plugin.annotate(&mut variants, data).unwrap();
        let fields = &variants[0].plugin_data;
        assert_eq!(fields["gnomAD_total_DP_cov"], "12345.0000");
        assert_eq!(fields["gnomAD_mean_cov"], "10.7570");
        assert_eq!(fields["gnomAD_median_approx_cov"], "8.0000");
        assert_eq!(fields["gnomAD_1x_cov"], "0.8307");
        assert_eq!(fields["gnomAD_50x_cov"], "0.0001");
        assert_eq!(fields["gnomAD_100x_cov"], "0.0000");
        assert_eq!(numeric("NA"), 0.0);
        assert_eq!(numeric(" 7 "), 7.0);
        assert_eq!(numeric("12abc"), 12.0);
        assert_eq!(numeric("1e5x"), 100000.0);
        assert_eq!(numeric("1.5.3"), 1.5);
        assert_eq!(numeric("e5"), 0.0);
    }

    #[test]
    fn run_batch_skips_when_uninitialized() {
        let plugin = GnomadcPlugin::new();
        let mut variants = vec![];
        assert!(plugin.run_batch(&mut variants).is_ok());
    }

    /// Ensembl VEP 116.3 never hands gnomADc a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). A `<DEL>` over the covered
    /// positions gets no coverage field, while the sequence variant beside it
    /// is averaged.
    #[test]
    fn a_structural_variant_gets_no_coverage_field() {
        let plugin = plugin(13);
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = variant(100, 101);
        deletion.ref_allele = b"N".to_vec();
        deletion.alt_alleles = vec![b"<DEL>".to_vec()];
        deletion.is_structural = true;
        let mut variants = vec![deletion, variant(100, 101)];
        let data = PrefetchData::new(batch(vec![
            record(100, &["30.1", "30", "1.0000", "0.9", "0.5"]),
            record(101, &["31.3", "32", "1.0000", "0.8", "0.5"]),
        ]));
        plugin.annotate(&mut variants, data).unwrap();
        assert!(variants[0].plugin_data.is_empty());
        assert_eq!(variants[1].plugin_data["gnomAD_mean_cov"], "30.7000");
    }

    /// `gnomADc.pm:185` splits the line the tabix iterator hands over with its
    /// line end alone removed, so a blank after the last cell stays in that
    /// cell, a cell of one space is a cell, and a line ending in a tab ends in
    /// an empty cell (`split` drops it, and the key it would fill adds nothing
    /// at `:171`); the addition of `:171` reads the number before a blank. The
    /// table and a store built from it hand the plugin the same cells, and the
    /// means are the ones Ensembl VEP prints for these rows.
    #[test]
    fn a_trailing_blank_row_reads_alike_through_the_table_and_a_store() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};

        const HEADER: &str = "#chr\tpos\tmean\tmedian_approx\ttotal_DP\tover_1\tover_5\tover_10\
                              \tover_15\tover_20\tover_25\tover_30\tover_50\tover_100";
        const ROWS: [&str; 8] = [
            "21\t13826645\t2.1859e+00\t0\t1597793\t4.0428e-01\t2.1310e-01\t3.3116e-02\t6.1304e-03\
             \t1.6198e-03\t7.4971e-04\t2.6951e-04\t1.3681e-06\t1.0000e-01 ",
            "21\t13826646\t2.1890e+00\t0\t1600055\t4.0428e-01\t2.1345e-01\t3.3284e-02\t6.2371e-03\
             \t1.6485e-03\t7.6066e-04\t2.7772e-04\t1.3681e-06\t2.0000e-01\t",
            "21\t13826647\t2.1938e+00\t0\t1603542\t4.0429e-01\t2.1383e-01\t3.3648e-02\t6.3944e-03\
             \t1.6896e-03\t7.7981e-04\t2.8456e-04\t1.3681e-06\t3.0000e-01",
            "21\t14011314\t3.1100e+01\t31\t22700000\t1.0000e+00\t9.9999e-01\t9.9990e-01\t9.9900e-01\
             \t9.9000e-01\t9.0000e-01\t5.0000e-01\t1.0000e-01\t1.0000e-02 ",
            "21\t14011316\t3.2200e+01\t32\t23500000\t1.0000e+00\t9.9998e-01\t9.9980e-01\t9.9800e-01\
             \t9.8000e-01\t8.0000e-01\t4.0000e-01\t2.0000e-01\t2.0000e-02\t",
            "21\t14011318\t3.3300e+01\t33\t24300000\t1.0000e+00\t9.9997e-01\t9.9970e-01\t9.9700e-01\
             \t9.7000e-01\t7.0000e-01\t3.0000e-01\t3.0000e-01\t ",
            "21\t14011320\t3.4400e+01\t34\t25100000\t1.0000e+00\t9.9996e-01\t9.9960e-01\t9.9600e-01\
             \t9.6000e-01\t6.0000e-01\t2.0000e-01\t4.0000e-01\t4.0000e+00  ",
            "21\t14011322\t3.5500e+01\t35\t25900000\t1.0000e+00\t9.9995e-01\t9.9950e-01\t9.9500e-01\
             \t9.5000e-01\t5.0000e-01\t1.0000e-01\t5.0000e-01\t5.0000e-02",
        ];
        let mut lines = vec![HEADER];
        lines.extend(ROWS);

        let table_dir = tempfile::tempdir().unwrap();
        let table = write_tabix_fixture(
            table_dir.path(),
            "coverage.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let store_dir = tempfile::tempdir().unwrap();
        let stored_table = write_tabix_fixture(
            store_dir.path(),
            "coverage.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let names: Vec<String> = HEADER
            .trim_start_matches('#')
            .split('\t')
            .map(String::from)
            .collect();
        let records: Vec<TabixRecord> = ROWS
            .iter()
            .map(|row| {
                let cells: Vec<&str> = row.split('\t').collect();
                record(cells[1].parse().unwrap(), &cells[2..])
            })
            .collect();
        BinaryStoreWriter::new(store_dir.path().join("coverage.tsv.vpd"))
            .write(
                &names,
                &[("21".to_string(), records)],
                &TabixConfig::default(),
                &ContentIdentity::of_file(&stored_table).unwrap(),
            )
            .unwrap();

        let mut through_table = GnomadcPlugin::new();
        through_table.init(&[table.display().to_string()]).unwrap();
        let mut through_store = GnomadcPlugin::new();
        through_store
            .init(&[stored_table.display().to_string()])
            .unwrap();

        // The tabix reader moves the chromosome out of a record's cells and a
        // store keeps it in them, so the first cell tells which path answered:
        // a store refused at `init` falls back to the table.
        for row in ROWS {
            let cells: Vec<&str> = row.split('\t').collect();
            let position: u64 = cells[1].parse().unwrap();
            for (source, plugin, chromosome_cell) in [
                ("table", &through_table, ""),
                ("store", &through_store, "21"),
            ] {
                let records = plugin
                    .store
                    .as_ref()
                    .unwrap()
                    .query("21", position, position)
                    .unwrap();
                assert_eq!(records.len(), 1, "{source} at {position}");
                assert_eq!(
                    records[0].columns[0], chromosome_cell,
                    "{source} at {position}"
                );
                assert_eq!(
                    records[0].columns[1..],
                    cells[1..],
                    "{source} at {position}"
                );
            }
        }

        let mut through_table_variants = vec![
            variant(13826645, 13826647),
            variant(14011314, 14011314),
            variant(14011316, 14011316),
            variant(14011318, 14011318),
            variant(14011320, 14011320),
            variant(14011322, 14011322),
        ];
        let mut through_store_variants = through_table_variants.clone();
        through_table
            .run_batch(&mut through_table_variants)
            .unwrap();
        through_store
            .run_batch(&mut through_store_variants)
            .unwrap();
        let hundred_x: Vec<&str> = through_table_variants
            .iter()
            .map(|variant| variant.plugin_data["gnomAD_100x_cov"].as_str())
            .collect();
        assert_eq!(
            hundred_x,
            ["0.2000", "0.0100", "0.0200", "0.0000", "4.0000", "0.0500"]
        );
        let deletion = &through_table_variants[0].plugin_data;
        assert_eq!(deletion["gnomAD_mean_cov"], "2.1896");
        assert_eq!(deletion["gnomAD_total_DP_cov"], "1600463.3333");
        assert_eq!(deletion["gnomAD_50x_cov"], "0.0000");
        assert_eq!(deletion.len(), 12);
        for (through_table, through_store) in
            through_table_variants.iter().zip(&through_store_variants)
        {
            assert_eq!(through_table.plugin_data, through_store.plugin_data);
        }
    }

    /// The gate of OutputFactory.pm:2229 is applied before the coverage file is
    /// read: a structural variant's span, one record per base in the file, is
    /// left out of the prefetch query, so a `<DEL>` over eleven covered
    /// positions beside a substitution at one of them fetches that one record
    /// alone, and the batch run writes on the substitution only.
    #[test]
    fn a_structural_variant_span_is_not_read_in_prefetch() {
        let dir = tempfile::tempdir().unwrap();
        let mut lines = vec![
            "#chrom\tpos\tmean\tmedian\tover_1\tover_5\tover_10\tover_15\tover_20\tover_25\tover_30\tover_50\tover_100"
                .to_string(),
        ];
        for pos in 100..=110 {
            lines.push(format!("21\t{pos}\t30\t30\t1\t1\t1\t1\t1\t1\t1\t0.5\t0"));
        }
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let file = write_tabix_fixture(dir.path(), "coverage.tsv.gz", &refs, IndexSpec::standard());
        let mut plugin = GnomadcPlugin::new();
        plugin.init(&[file.display().to_string()]).unwrap();

        let mut deletion = variant(100, 110);
        deletion.ref_allele = b"N".to_vec();
        deletion.alt_alleles = vec![b"<DEL>".to_vec()];
        deletion.is_structural = true;
        let mut variants = vec![deletion, variant(105, 105)];

        let fetched: BatchQueryResult = plugin
            .prefetch(&variants)
            .unwrap()
            .expect("an open store prefetches")
            .downcast()
            .unwrap();
        assert_eq!(fetched.lookup("21", 100, 110).len(), 1);

        plugin.run_batch(&mut variants).unwrap();
        assert!(variants[0].plugin_data.is_empty());
        assert_eq!(variants[1].plugin_data["gnomAD_mean_cov"], "30.0000");
    }
}
