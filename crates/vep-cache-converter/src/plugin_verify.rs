// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Query concordance against the original tabix store, including interval
//! boundaries, absence, record order and multiplicity in different batches.

use std::collections::BTreeSet;
use std::path::Path;

use anyhow::{Context, Result};
use vep_builtins::annotation_store::AnnotationStore;
use vep_builtins::binary_store::BinaryAnnotator;
use vep_builtins::tabix::{TabixAnnotator, TabixConfig, TabixRecord};

use crate::TabixConvertConfig;

type Query<'a> = (&'a str, u64, u64);

#[derive(Default)]
pub(crate) struct VerifyResult {
    pub checked: usize,
    pub probes: usize,
    pub mismatches: usize,
    pub mismatch_details: Vec<String>,
}

pub(crate) fn verify_concordance(
    vpd_path: &Path,
    input_path: &Path,
    cfg: &TabixConvertConfig,
    source_records: &[(String, Vec<TabixRecord>)],
    verify_all: bool,
) -> Result<VerifyResult> {
    let binary = BinaryAnnotator::open(vpd_path)
        .with_context(|| format!("failed to open {} for verification", vpd_path.display()))?;
    let source = TabixAnnotator::open(TabixConfig {
        file_path: input_path.to_path_buf(),
        chr_col: cfg.chr_col,
        start_col: cfg.start_col,
        end_col: cfg.end_col,
        ref_col: None,
        alt_col: None,
        zero_based: cfg.zero_based,
    })
    .with_context(|| format!("failed to open {} for verification", input_path.display()))?;
    let queries = query_probes(source_records, verify_all);
    verify_queries(&binary, &source, &queries, cfg.chr_col)
}

fn query_probes(source_records: &[(String, Vec<TabixRecord>)], verify_all: bool) -> Vec<Query<'_>> {
    let total: usize = source_records
        .iter()
        .map(|(_, records)| records.len())
        .sum();
    let step = if verify_all {
        1
    } else {
        total.div_ceil(1000).max(1)
    };
    let mut ordinal = 0;
    let mut queries = BTreeSet::new();
    for (chr, records) in source_records {
        let mut reach = 0u64;
        for record in records {
            if ordinal % step == 0 {
                let (start, end) = (record.start, record.end);
                for pos in [start, start + (end - start) / 2, end] {
                    queries.insert((chr.as_str(), pos, pos));
                }
                queries.insert((chr.as_str(), start, end));
                if start > 1 {
                    queries.insert((chr.as_str(), start - 1, start - 1));
                }
                if let Some(after) = end.checked_add(1) {
                    queries.insert((chr.as_str(), after, after));
                }
                if start.saturating_sub(reach) > 1 {
                    queries.insert((chr.as_str(), reach + 1, start - 1));
                }
            }
            reach = reach.max(record.end);
            ordinal += 1;
        }
        if !records.is_empty() {
            if let Some(after) = reach.checked_add(1) {
                queries.insert((chr.as_str(), after, after));
            }
        }
    }
    queries.into_iter().collect()
}

fn verify_queries(
    binary: &dyn AnnotationStore,
    source: &dyn AnnotationStore,
    queries: &[Query<'_>],
    chr_col: usize,
) -> Result<VerifyResult> {
    let mut result = VerifyResult {
        probes: queries.len(),
        ..VerifyResult::default()
    };
    for chunk in queries.chunks(64) {
        let mut expected = Vec::with_capacity(chunk.len());
        for &(chr, start, end) in chunk {
            let query = (chr, start, end);
            let records = source
                .query(chr, start, end)
                .with_context(|| format!("tabix query {chr}:{start}-{end} failed"))?;
            let actual = binary
                .query(chr, start, end)
                .with_context(|| format!("binary single query {chr}:{start}-{end} failed"))?;
            result.compare(
                "single",
                query,
                &records,
                &actual.iter().collect::<Vec<_>>(),
                chr_col,
            );
            let isolated = binary
                .query_batch(&[query])
                .with_context(|| format!("binary isolated batch {chr}:{start}-{end} failed"))?;
            result.compare(
                "isolated batch",
                query,
                &records,
                &isolated.lookup(chr, start, end),
                chr_col,
            );
            expected.push(records);
        }
        let forward: Vec<usize> = (0..chunk.len()).collect();
        let mut reverse = forward.iter().rev().copied().collect::<Vec<_>>();
        reverse.push(0);
        let sparse: Vec<usize> = forward.iter().step_by(2).copied().collect();
        for (mode, indices) in [
            ("batch", forward),
            ("reordered batch with repeated query", reverse),
            ("sparse batch", sparse),
        ] {
            let batch_queries: Vec<_> = indices.iter().map(|&i| chunk[i]).collect();
            let batch = binary
                .query_batch(&batch_queries)
                .with_context(|| format!("binary {mode} failed for {batch_queries:?}"))?;
            for i in indices {
                let (chr, start, end) = chunk[i];
                result.compare(
                    mode,
                    chunk[i],
                    &expected[i],
                    &batch.lookup(chr, start, end),
                    chr_col,
                );
            }
        }
    }
    Ok(result)
}

impl VerifyResult {
    fn compare(
        &mut self,
        mode: &str,
        (chr, start, end): Query<'_>,
        expected: &[TabixRecord],
        actual: &[&TabixRecord],
        chr_col: usize,
    ) {
        self.checked += 1;
        if let Some(difference) = first_difference(expected, actual, chr_col) {
            self.mismatches += 1;
            if self.mismatch_details.len() < 10 {
                self.mismatch_details
                    .push(format!("{mode} {chr}:{start}-{end}: {difference}"));
            }
        }
    }
}

fn first_difference(
    expected: &[TabixRecord],
    actual: &[&TabixRecord],
    chr_col: usize,
) -> Option<String> {
    if expected.len() != actual.len() {
        return Some(format!(
            "record count mismatch: tabix={}, binary={}",
            expected.len(),
            actual.len()
        ));
    }
    for (i, (src, bin)) in expected.iter().zip(actual).enumerate() {
        if (src.chr.as_str(), src.start, src.end) != (bin.chr.as_str(), bin.start, bin.end) {
            return Some(format!(
                "record[{i}] bounds: tabix={}:{}-{}, binary={}:{}-{}",
                src.chr, src.start, src.end, bin.chr, bin.start, bin.end
            ));
        }
        if src.columns.len() != bin.columns.len() {
            return Some(format!(
                "record[{i}] column count: tabix={}, binary={}",
                src.columns.len(),
                bin.columns.len()
            ));
        }
        for (col, (src_value, bin_value)) in src.columns.iter().zip(&bin.columns).enumerate() {
            let src_value = if col == chr_col { &src.chr } else { src_value };
            if src_value != bin_value {
                let name = src.header.as_ref().and_then(|header| header.get(col));
                return Some(format!(
                    "record[{i}] column[{col}] {}: tabix={src_value:?}, binary={bin_value:?}",
                    name.map(String::as_str).unwrap_or("?")
                ));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeMap, HashMap};
    use vep_builtins::tabix::BatchQueryResult;
    use vep_builtins::PluginError;

    struct Records {
        rows: Vec<TabixRecord>,
        starts_only: bool,
        reorder_batches: bool,
    }

    impl AnnotationStore for Records {
        fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError> {
            Ok(self
                .rows
                .iter()
                .filter(|r| {
                    r.chr == chr
                        && r.start <= end
                        && r.end >= start
                        && (!self.starts_only || r.start >= start)
                })
                .cloned()
                .collect())
        }

        fn query_batch(&self, queries: &[Query<'_>]) -> Result<BatchQueryResult, PluginError> {
            let mut rows: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();
            let mut selected: Vec<_> = self
                .rows
                .iter()
                .filter(|r| {
                    queries.iter().any(|&(chr, start, end)| {
                        r.chr == chr
                            && r.start <= end
                            && r.end >= start
                            && (!self.starts_only || r.start >= start)
                    })
                })
                .cloned()
                .collect();
            if self.reorder_batches && queries.len() > 1 {
                selected.reverse();
            }
            for row in selected {
                rows.entry(row.chr.clone())
                    .or_default()
                    .entry(row.start)
                    .or_default()
                    .push(row);
            }
            Ok(BatchQueryResult::from_records(rows))
        }

        fn header(&self) -> Option<&[String]> {
            None
        }
    }

    fn rows() -> Vec<TabixRecord> {
        [(150, 260, "long"), (150, 200, "short"), (300, 300, "point")]
            .map(|(start, end, name)| TabixRecord {
                chr: "21".into(),
                start,
                end,
                columns: vec!["21".into(), name.into()],
                header: None,
            })
            .to_vec()
    }

    #[test]
    fn probes_cover_interiors_endpoints_and_true_gaps() {
        let source = vec![("21".into(), rows())];
        let probes = query_probes(&source, true);
        for query in [
            ("21", 150, 150),
            ("21", 205, 205),
            ("21", 260, 260),
            ("21", 261, 261),
            ("21", 261, 299),
            ("21", 150, 260),
            ("21", 1, 149),
            ("21", 301, 301),
        ] {
            assert!(probes.contains(&query), "{query:?}");
        }
        assert_eq!(probes, query_probes(&source, false));
    }

    #[test]
    fn verification_rejects_start_only_and_batch_order_errors() {
        let source = Records {
            rows: rows(),
            starts_only: false,
            reorder_batches: false,
        };
        let queries = [("21", 150, 150), ("21", 205, 205), ("21", 261, 299)];
        for (starts_only, reorder_batches, diagnostic) in [
            (true, false, "single 21:205-205: record count mismatch"),
            (false, true, "batch 21:150-150: record[0] bounds"),
        ] {
            let binary = Records {
                rows: rows(),
                starts_only,
                reorder_batches,
            };
            let result = verify_queries(&binary, &source, &queries, 0).unwrap();
            assert!(result.mismatches > 0);
            assert!(
                result
                    .mismatch_details
                    .iter()
                    .any(|d| d.contains(diagnostic)),
                "{:?}",
                result.mismatch_details
            );
        }
        assert_eq!(
            verify_queries(&source, &source, &queries, 0)
                .unwrap()
                .mismatches,
            0
        );
    }

    #[test]
    fn verification_reports_duplicate_loss_and_wrong_columns() {
        let source = Records {
            rows: vec![rows()[0].clone(); 2],
            starts_only: false,
            reorder_batches: false,
        };
        let mut binary = Records {
            rows: vec![rows()[0].clone()],
            starts_only: false,
            reorder_batches: false,
        };
        let queries = [("21", 150, 150)];
        let result = verify_queries(&binary, &source, &queries, 0).unwrap();
        assert!(result.mismatch_details[0].contains("tabix=2, binary=1"));
        binary.rows.push(rows()[0].clone());
        binary.rows[1].columns[1] = "wrong".into();
        let result = verify_queries(&binary, &source, &queries, 0).unwrap();
        assert!(result.mismatch_details[0].contains("record[1] column[1]"));
        assert!(result.mismatch_details[0].contains("binary=\"wrong\""));
    }
}
