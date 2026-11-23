// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a run reads beside its outputs for the classes of a plugin: the GWAS catalog, the
//! REVEL table, the seeded UTRAnnotator outputs, the candidate JSON cache and the SpliceAI
//! score files, each read in the order the command line's options are taken, with what the
//! report states of each.

use super::api::{Classes, Keyword, PluginReads, Reads};
use super::interrupt;
use super::paths::{self, path_text, OsFailure};
use super::reader::{digest, Reader};
use super::Stop;
use crate::rules::downstream::DownstreamSequences;
use crate::rules::exon_spans::ExonSpans;
use crate::rules::gwas::GwasCatalogPairs;
use crate::rules::maxentscan::SpliceWindows;
use crate::rules::revel::RevelTable;
use crate::rules::spliceai::{RefusalCounts, SpliceAiTable};
use crate::rules::utrannotator::UtrannotatorSeedSweep;
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// The registry's name of the plugin whose associations file `--gwas-catalog` is.
pub const GWAS_PLUGIN: &str = "GWAS";
/// The registry's name of the plugin whose table `--revel-table` is.
pub const REVEL_PLUGIN: &str = "REVEL";
/// The registry's name of the plugin whose classes read the outputs of `--utrannotator-sweep`.
pub const UTRANNOTATOR_PLUGIN: &str = "UTRAnnotator";
/// The registry's name of the plugin whose classes read the cache of `--json-cache`.
pub const DOWNSTREAM_PLUGIN: &str = "Downstream";
/// The registry's name of the plugin whose class reads that cache's transcripts likewise.
pub const BLOSUM62_PLUGIN: &str = "Blosum62";
/// The registry's name of the plugin whose class reads that cache's introns.
pub const MAXENTSCAN_PLUGIN: &str = "MaxEntScan";
/// The registry's name of the plugin whose class `hgvs_insertion_left_where_the_shift_stops`
/// reads that cache's exon spans.
pub const HGVS_INTRON_OFFSET_PLUGIN: &str = "HGVSIntronOffset";
/// The registry's name of the plugin whose score files `--spliceai-table` names.
pub const SPLICEAI_PLUGIN: &str = "SpliceAI";
/// The bytes of a file read to tell a default-format output.
const HEAD_BYTES: usize = 1 << 16;

/// The files and directories a run names for its classes, each path as the command line
/// gives it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Named {
    pub gwas_catalog: Option<String>,
    pub revel_table: Option<String>,
    pub revel_no_match: bool,
    pub utrannotator_sweep: Vec<String>,
    pub json_cache: Option<String>,
    pub spliceai_table: Vec<String>,
}

/// What a run read for its classes.
pub struct Context {
    pub reads: Reads,
    /// Per registry plugin name, the keywords its classes are handed, in the order read.
    pub plugin_reads: PluginReads,
    /// The report's `gwas_catalog`: the catalog and its sha256, or null.
    pub gwas_catalog: Value,
    /// The report's `revel_table`: the table, its sha256 and whether the run set `no_match`.
    pub revel_table: Value,
    /// The report's `utrannotator_sweep`: each output read with its sha256.
    pub utrannotator_sweep: Value,
    /// The report's `json_cache`: the directory, by path alone.
    pub json_cache: Value,
    /// The report's `spliceai_table`: the files, by path alone.
    pub spliceai_table: Value,
}

fn digested(path: &str) -> Result<Value, Stop> {
    let sha256 = digest(path).map_err(Stop::Os)?;
    Ok(json!({"path": path, "sha256": sha256}))
}

/// Whether the file opens as a default-format output: its first line that is no `##`
/// comment is the column header.
pub fn names_default_columns(path: &str) -> bool {
    let Ok(mut reader) = Reader::open(path) else {
        return false;
    };
    let Ok(head) = reader.head(HEAD_BYTES) else {
        return false;
    };
    head.split(|&byte| byte == b'\n')
        .find(|line| !line.starts_with(b"##"))
        .is_some_and(|line| line.starts_with(b"#Uploaded_variation\t"))
}

/// The outputs `paths` name: a file as it is, and of a directory every file in it that
/// opens as a default-format output, in name order.
pub fn default_format_outputs(paths: &[String]) -> Result<Vec<String>, Stop> {
    let mut outputs = Vec::new();
    for path in paths {
        if !paths::is_dir(path) {
            outputs.push(path.clone());
            continue;
        }
        let listed = fs::read_dir(path)
            .and_then(|entries| entries.collect::<Result<Vec<_>, _>>())
            .map_err(|error| Stop::Os(OsFailure::of(&error, Some(path))))?;
        let mut inside: Vec<String> = listed
            .iter()
            .map(|entry| {
                let name = entry.file_name();
                let joined = Path::new(path).join(name);
                path_text(&joined.to_string_lossy())
            })
            .filter(|entry| paths::is_file(entry) && names_default_columns(entry))
            .collect();
        if inside.is_empty() {
            return Err(Stop::Refused(format!(
                "{path}: no default-format output in the directory"
            )));
        }
        inside.sort();
        outputs.extend(inside);
    }
    Ok(outputs)
}

impl Context {
    /// Reads what `named` names. The classes the table declares for `Blosum62` are handed
    /// the cache's sequences, those it declares for `MaxEntScan` its introns and those it
    /// declares for `HGVSIntronOffset` its exon spans, only where the table declares any for
    /// the plugin.
    pub fn read(named: &Named, classes: &Classes) -> Result<Context, Stop> {
        let stopped = || match interrupt::requested() {
            true => Err(Stop::Interrupted),
            false => Ok(()),
        };
        let mut reads = Reads::default();
        let mut plugin_reads = Vec::new();
        let mut hand = |plugin: &str, keywords: &[Keyword]| {
            plugin_reads.push((plugin.to_string(), keywords.to_vec()));
        };
        let mut gwas_catalog = Value::Null;
        if let Some(path) = &named.gwas_catalog {
            gwas_catalog = digested(path)?;
            reads.gwas_pairs = Some(GwasCatalogPairs::from_catalog(Path::new(path))?);
            hand(GWAS_PLUGIN, &[Keyword::GwasPairs]);
            stopped()?;
        }
        let mut revel_table = Value::Null;
        if let Some(path) = &named.revel_table {
            revel_table = digested(path)?;
            revel_table["no_match"] = Value::Bool(named.revel_no_match);
            reads.revel_table = Some(RevelTable::new(Path::new(path), named.revel_no_match)?);
            hand(REVEL_PLUGIN, &[Keyword::RevelTable]);
            stopped()?;
        }
        let mut utrannotator_sweep = Value::Null;
        if !named.utrannotator_sweep.is_empty() {
            let outputs = default_format_outputs(&named.utrannotator_sweep)?;
            let mut listed = Vec::new();
            for path in &outputs {
                listed.push(digested(path)?);
            }
            utrannotator_sweep = Value::Array(listed);
            let outputs: Vec<PathBuf> = outputs.iter().map(PathBuf::from).collect();
            reads.utrannotator_sweep = Some(UtrannotatorSeedSweep::from_outputs(&outputs)?);
            hand(UTRANNOTATOR_PLUGIN, &[Keyword::UtrannotatorSweep]);
            stopped()?;
        }
        let mut json_cache = Value::Null;
        if let Some(path) = &named.json_cache {
            json_cache = json!({"path": path});
            reads.downstream_sequences = Some(Arc::new(DownstreamSequences::new(Path::new(path))?));
            hand(
                DOWNSTREAM_PLUGIN,
                &[Keyword::DownstreamSequences, Keyword::Declared],
            );
            if !classes.of_plugin(BLOSUM62_PLUGIN).is_empty() {
                hand(BLOSUM62_PLUGIN, &[Keyword::DownstreamSequences]);
            }
            if !classes.of_plugin(MAXENTSCAN_PLUGIN).is_empty() {
                reads.splice_windows = Some(SpliceWindows::new(Path::new(path))?);
                hand(MAXENTSCAN_PLUGIN, &[Keyword::SpliceWindows]);
            }
            if !classes.of_plugin(HGVS_INTRON_OFFSET_PLUGIN).is_empty() {
                reads.exon_spans = Some(ExonSpans::new(Path::new(path))?);
                hand(HGVS_INTRON_OFFSET_PLUGIN, &[Keyword::ExonSpans]);
            }
            stopped()?;
        }
        let mut spliceai_table = Value::Null;
        if !named.spliceai_table.is_empty() {
            spliceai_table = Value::Array(
                named
                    .spliceai_table
                    .iter()
                    .map(|path| json!({"path": path}))
                    .collect(),
            );
            let files: Vec<PathBuf> = named.spliceai_table.iter().map(PathBuf::from).collect();
            reads.spliceai_table = Some(SpliceAiTable::new(
                &files,
                named.json_cache.as_deref().map(Path::new),
            )?);
            hand(SPLICEAI_PLUGIN, &[Keyword::SpliceaiTable]);
            stopped()?;
        }
        Ok(Context {
            reads,
            plugin_reads: PluginReads(plugin_reads),
            gwas_catalog,
            revel_table,
            utrannotator_sweep,
            json_cache,
            spliceai_table,
        })
    }
}

/// Forgets the rows the score files were asked about, before an engine is scored.
pub fn clear_refusals(reads: &Reads) {
    if let Some(table) = &reads.spliceai_table {
        table.clear_refusals();
    }
}

/// The rows the score files were asked about since they were last forgotten; `None`
/// without score files.
pub fn refusals(reads: &Reads) -> Option<RefusalCounts> {
    reads
        .spliceai_table
        .as_ref()
        .map(SpliceAiTable::refusal_counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir().join(format!(
                "vep-plugin-pairs-reads-{}-{tag}",
                std::process::id()
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Scratch(path)
        }

        fn file(&self, name: &str, content: &[u8]) -> String {
            let path = self.0.join(name);
            fs::write(&path, content).unwrap();
            path.to_str().unwrap().to_string()
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn gzip(content: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(content).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn directory_gives_its_default_format_outputs_in_name_order() {
        let scratch = Scratch::new("outputs");
        let headed = format!("## a comment\n{HEADER}");
        fs::create_dir(scratch.0.join("sweep")).unwrap();
        scratch.file("sweep/seed2.txt", headed.as_bytes());
        scratch.file("sweep/seed10.txt.gz", &gzip(headed.as_bytes()));
        scratch.file("sweep/notes.txt", b"not an output\n");
        scratch.file("sweep/blank-first.txt", format!("\n{HEADER}").as_bytes());
        let packed = gzip(headed.as_bytes());
        scratch.file("sweep/cut.txt.gz", &packed[..packed.len() - 4]);
        fs::create_dir(scratch.0.join("sweep/inner")).unwrap();
        let single = scratch.file("single.txt", b"anything\n");
        let sweep = scratch.0.join("sweep").to_str().unwrap().to_string();
        assert_eq!(
            default_format_outputs(&[single.clone(), sweep.clone()]).unwrap_or_default(),
            [
                single,
                format!("{sweep}/seed10.txt.gz"),
                format!("{sweep}/seed2.txt")
            ]
        );
        fs::create_dir(scratch.0.join("none")).unwrap();
        scratch.file("none/notes.txt", b"not an output\n");
        let none = scratch.0.join("none").to_str().unwrap().to_string();
        let Err(Stop::Refused(message)) = default_format_outputs(std::slice::from_ref(&none))
        else {
            panic!("a directory without an output was taken");
        };
        assert_eq!(
            message,
            format!("{none}: no default-format output in the directory")
        );
    }
}
