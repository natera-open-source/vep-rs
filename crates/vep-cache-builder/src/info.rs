// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cache's metadata files: `info.json`, read by the runtime's cache
//! resolution and printed in the output headers, and `chr_synonyms.txt`, the
//! chromosome synonyms the runtime resolves input names through.
//!
//! `info.json` carries what Ensembl VEP's `info.txt` carries for the cache
//! (`Pipeline/DumpVEP/Dumper/Core.pm` `dump_info`: species, assembly, the
//! `sift` and `polyphen` markers, the source versions) in the form the
//! converter writes, plus `contig_order`, and `variation_cols` with `var_type
//! tabix` only when the cache holds variation records
//! (`Pipeline/DumpVEP/Dumper/Variation.pm:140, 177-197`). `chr_synonyms.txt`
//! is the file the Ensembl VEP dump pipeline derives from the core
//! `seq_region_synonym` table (`Pipeline/DumpVEP/CreateDumpJobs.pm:187-212`,
//! `BaseVEP.pm:596-621`).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;
use tracing::info;

use crate::builder::BuildStats;
use crate::core_dumps::{id_to_name, read_table};

/// The `info.json` object the builder writes.
#[derive(Serialize)]
pub struct CacheInfoJson {
    pub species: String,
    pub assembly: String,
    pub cache_version: u32,
    /// `b` when the cache carries SIFT matrices.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sift: Option<String>,
    /// `b` when the cache carries PolyPhen matrices.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polyphen: Option<String>,
    pub regulatory: bool,
    /// The columns of the variation files, present only when the files were
    /// written.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub variation_cols: Option<Vec<String>>,
    /// `tabix` when the variation files are written: one indexed `all_vars.gz`
    /// per chromosome, the layout `CacheDir.pm:184` selects by this key.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub var_type: Option<String>,
    pub source_versions: BTreeMap<String, String>,
    /// The chromosome directories written, in [`contig_order`].
    pub contig_order: Vec<String>,
}

impl CacheInfoJson {
    /// The metadata of a build: `sift` and `polyphen` read `b` only when a
    /// matrix of that tool was written, `variation_cols` and `var_type` are
    /// present only when variation files were (`variation_cols` names their
    /// columns), and `contig_order` lists the chromosomes the build wrote.
    pub fn new(
        species: &str,
        assembly: &str,
        cache_version: u32,
        stats: &BuildStats,
        variation_cols: Option<Vec<String>>,
        source_versions: BTreeMap<String, String>,
    ) -> Self {
        let marker = |written: usize| (written > 0).then(|| "b".to_string());
        Self {
            species: species.to_string(),
            assembly: assembly.to_string(),
            cache_version,
            sift: marker(stats.sift_matrices),
            polyphen: marker(stats.polyphen_matrices),
            regulatory: false,
            var_type: variation_cols.is_some().then(|| "tabix".to_string()),
            variation_cols,
            source_versions,
            contig_order: stats.chromosomes.clone(),
        }
    }

    /// Writes `info.json` into the cache directory.
    pub fn write(&self, output_dir: &Path) -> Result<()> {
        let path = output_dir.join("info.json");
        let json = serde_json::to_string_pretty(self).context("Failed to serialize info.json")?;
        std::fs::write(&path, json)
            .with_context(|| format!("Failed to write {}", path.display()))?;
        info!("Wrote {}", path.display());
        Ok(())
    }
}

/// Chromosome names in the order a cache lists them: autosomes by number, then
/// X, Y and MT, then every other contig by name, the converter's order.
pub fn contig_order(names: impl IntoIterator<Item = String>) -> Vec<String> {
    let rank = |name: &str| -> (u8, u64, String) {
        if let Ok(number) = name.parse::<u64>() {
            return (0, number, String::new());
        }
        let tier = match name {
            "X" => 1,
            "Y" => 2,
            "MT" => 3,
            _ => 4,
        };
        (tier, 0, name.to_string())
    };
    let mut names: Vec<String> = names.into_iter().collect();
    names.sort_by_cached_key(|name| rank(name));
    names
}

/// The synonym pairs the Ensembl VEP dump pipeline writes to
/// `chr_synonyms.txt` from the `seq_region_synonym` table
/// (`CreateDumpJobs.pm:200-209` over the tree `BaseVEP.pm:596-621` builds):
/// each row, taken in `seq_region_synonym_id` order (the `ORDER BY` of
/// `SeqRegionSynonymAdaptor::_final_clause`, which the dump file's row order
/// does not follow), links a sequence region's name to its synonym; when the
/// table holds at most 40,000 rows, each new pair is also linked to every
/// name already linked to either side (one level of indirection, `:606-618`,
/// whose result depends on the row order); the output is every pair of
/// distinct linked names, ordered within the pair and here across the lines.
/// `seq_region.txt` and `seq_region_synonym.txt` (or their `.gz`) are read
/// from `dumps_dir`.
pub fn chr_synonym_pairs(dumps_dir: &Path) -> Result<Vec<(String, String)>> {
    let region_names = id_to_name(dumps_dir, "seq_region", 4)?;
    let mut rows: Vec<(u64, String, String)> = read_table(dumps_dir, "seq_region_synonym", 4)?
        .map(|row| {
            let row = row?;
            Ok((
                row.number(0)?,
                row.field(1).to_string(),
                row.field(2).to_string(),
            ))
        })
        .collect::<Result<_>>()?;
    rows.sort_by_key(|(id, _, _)| *id);
    let indirect = rows.len() <= 40_000;
    let mut tree: HashMap<String, BTreeSet<String>> = HashMap::new();
    for (_, region_id, synonym) in rows {
        let Some(region_name) = region_names.get(&region_id) else {
            continue;
        };
        let (a, b) = ordered(region_name.clone(), synonym);
        tree.entry(a.clone()).or_default().insert(b.clone());
        if indirect {
            let linked_to_a: BTreeSet<String> = tree.get(&a).cloned().unwrap_or_default();
            for name in linked_to_a {
                tree.entry(name).or_default().insert(b.clone());
            }
            let linked_to_b: BTreeSet<String> = tree.get(&b).cloned().unwrap_or_default();
            for name in linked_to_b {
                tree.entry(name).or_default().insert(a.clone());
            }
        }
    }
    let pairs: BTreeSet<(String, String)> = tree
        .iter()
        .flat_map(|(a, linked)| {
            linked
                .iter()
                .filter(move |b| *b != a)
                .map(move |b| ordered(a.clone(), b.clone()))
        })
        .collect();
    Ok(pairs.into_iter().collect())
}

/// A pair in string order, as Perl's `sort` places it.
fn ordered(a: String, b: String) -> (String, String) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Writes `chr_synonyms.txt` beside `info.json`, one `a<TAB>b` pair per line,
/// and returns the number of lines.
pub fn write_chr_synonyms(output_dir: &Path, dumps_dir: &Path) -> Result<usize> {
    let pairs = chr_synonym_pairs(dumps_dir)?;
    let text: String = pairs.iter().map(|(a, b)| format!("{a}\t{b}\n")).collect();
    let path = output_dir.join("chr_synonyms.txt");
    std::fs::write(&path, text).with_context(|| format!("Failed to write {}", path.display()))?;
    info!("Wrote {} ({} synonym pairs)", path.display(), pairs.len());
    Ok(pairs.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_dir() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
    }

    /// The pairs the release 116 GRCh38 cache's `chr_synonyms.txt` lists among
    /// the names of chromosomes 21 and MT: three synonym rows each in the
    /// dump, linked to one another through the one level of indirection, so
    /// every pair of the four names appears, twelve lines in all.
    #[test]
    fn chr_synonyms_match_the_ensembl_dumper() {
        let pairs = chr_synonym_pairs(&fixture_dir()).unwrap();
        let lines: Vec<String> = pairs.iter().map(|(a, b)| format!("{a}\t{b}")).collect();
        assert_eq!(
            lines,
            [
                "21\tCM000683.2",
                "21\tNC_000021.9",
                "21\tchr21",
                "CM000683.2\tNC_000021.9",
                "CM000683.2\tchr21",
                "J01415.2\tMT",
                "J01415.2\tNC_012920.1",
                "J01415.2\tchrM",
                "MT\tNC_012920.1",
                "MT\tchrM",
                "NC_000021.9\tchr21",
                "NC_012920.1\tchrM",
            ]
        );
        let dir = tempfile::tempdir().unwrap();
        assert_eq!(write_chr_synonyms(dir.path(), &fixture_dir()).unwrap(), 12);
        let written = std::fs::read_to_string(dir.path().join("chr_synonyms.txt")).unwrap();
        assert_eq!(written.lines().count(), 12);
        assert!(written.starts_with("21\tCM000683.2\n"));
    }

    /// Rows are taken in `seq_region_synonym_id` order whatever the file's
    /// order, and the indirection is one level, not a closure: the chain
    /// A-B, B-C, C-D yields A-C and B-D but no A-D, which the file's order
    /// (B-C first) would have added.
    #[test]
    fn chr_synonym_rows_follow_the_id_order_not_the_file_order() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(
            dir.path().join("seq_region.txt"),
            "1\tA\t4\t100\n2\tB\t4\t100\n3\tC\t4\t100\n",
        )
        .unwrap();
        std::fs::write(
            dir.path().join("seq_region_synonym.txt"),
            "20\t2\tC\t\\N\n10\t1\tB\t\\N\n30\t3\tD\t\\N\n",
        )
        .unwrap();
        let pairs = chr_synonym_pairs(dir.path()).unwrap();
        let lines: Vec<String> = pairs.iter().map(|(a, b)| format!("{a}\t{b}")).collect();
        assert_eq!(lines, ["A\tB", "A\tC", "B\tC", "B\tD", "C\tD"]);
    }

    /// Autosomes by number, X, Y, MT, then the rest by name.
    #[test]
    fn contig_order_follows_the_converter() {
        let names = [
            "KI270728.1",
            "X",
            "2",
            "MT",
            "10",
            "1",
            "Y",
            "GL000009.2",
            "21",
        ]
        .map(str::to_string);
        assert_eq!(
            contig_order(names),
            [
                "1",
                "2",
                "10",
                "21",
                "X",
                "Y",
                "MT",
                "GL000009.2",
                "KI270728.1"
            ]
        );
    }

    /// `sift` and `polyphen` read `b` only for a tool whose matrices were
    /// written, `variation_cols` and `var_type tabix` are present only with
    /// variation files, and the species, assembly, version, source versions
    /// and contig order are written as given.
    #[test]
    fn info_json_markers_follow_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let versions = BTreeMap::from([
            ("assembly".to_string(), "GRCh38.p14".to_string()),
            ("gencode".to_string(), "GENCODE 50".to_string()),
        ]);
        let stats = BuildStats {
            transcript_count: 0,
            chromosome_count: 2,
            region_file_count: 0,
            chromosomes: vec!["21".to_string(), "MT".to_string()],
            sift_matrices: 3,
            polyphen_matrices: 0,
        };
        CacheInfoJson::new("homo_sapiens", "GRCh38", 116, &stats, None, versions)
            .write(dir.path())
            .unwrap();
        let parsed: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(dir.path().join("info.json")).unwrap())
                .unwrap();
        assert_eq!(parsed["species"], "homo_sapiens");
        assert_eq!(parsed["assembly"], "GRCh38");
        assert_eq!(parsed["cache_version"], 116);
        assert_eq!(parsed["sift"], "b");
        assert!(parsed.get("polyphen").is_none());
        assert!(parsed.get("variation_cols").is_none());
        assert!(parsed.get("var_type").is_none());
        assert_eq!(parsed["regulatory"], false);
        assert_eq!(parsed["source_versions"]["assembly"], "GRCh38.p14");
        assert_eq!(parsed["source_versions"]["gencode"], "GENCODE 50");
        assert_eq!(parsed["contig_order"], serde_json::json!(["21", "MT"]));

        let stats = BuildStats {
            sift_matrices: 0,
            polyphen_matrices: 2,
            chromosomes: vec![],
            ..stats
        };
        let columns = crate::variation::variation_columns(false);
        let with_variations = CacheInfoJson::new(
            "homo_sapiens",
            "GRCh37",
            116,
            &stats,
            Some(columns.clone()),
            BTreeMap::new(),
        );
        let parsed = serde_json::to_value(&with_variations).unwrap();
        assert!(parsed.get("sift").is_none());
        assert_eq!(parsed["polyphen"], "b");
        assert_eq!(parsed["var_type"], "tabix");
        assert_eq!(parsed["variation_cols"], serde_json::json!(columns));
    }
}
