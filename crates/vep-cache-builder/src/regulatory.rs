// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The regulatory tree of a native cache, baked from the regulatory GFF files
//! the Ensembl site distributes (`--regulatory-gff`, `--emars-gff`,
//! `--motif-gff`): `regulatory/<chr>/<start>-<end>.json` shards in the shape
//! `vep-cli` reads, one feature per record, and `regulatory/binding_matrices.json`
//! holding a stub per motif matrix (its identifier, the motifs' span and their
//! factors, no counts), so a run with `--regulatory` over the cache gives the
//! rows a run with `--regulatory_gff` over the same files gives. The records
//! read as Ensembl VEP's GFF regulatory source reads them
//! (`AnnotationSource/File/RegFeat.pm` at release/116.2), through the
//! vocabulary of [`vep_core::regulatory_gff`].

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use tracing::{info, warn};
use vep_core::coordinate::is_cache_path_component;
use vep_core::regulatory_gff::{
    decode_attribute, decode_type, duplicate_identifier_warning, feature_id, ignored_bound_warning,
    transcription_factors, unsupported_type_warning, EXTENDED_BOUNDS_TYPES, MOTIF_FEATURE_TYPES,
    REGULATORY_FEATURE_TYPES,
};

use crate::builder::{region_start, REGION_SIZE};
use crate::core_dumps::BoundedLines;
use crate::gff3::normalize_chr;

/// The longest GFF line read, as the transcript GFF3 reader bounds it.
const MAX_LINE_BYTES: usize = 16 << 20;

/// The regulatory inputs of a build.
#[derive(Debug, Clone, Default)]
pub struct RegulatoryInputs {
    /// The regulatory features (`regulatory_features.gff3.gz`).
    pub regulatory_gff: Option<PathBuf>,
    /// The EMAR features (`EMARs.gff.gz`), regulatory features of their own
    /// identifier space.
    pub emars_gff: Option<PathBuf>,
    /// The motif features (`motif_features.gff3.gz`).
    pub motif_gff: Option<PathBuf>,
    /// Write each promoter at its `extended_start` and `extended_end`.
    pub extended_promoters: bool,
}

impl RegulatoryInputs {
    pub fn is_empty(&self) -> bool {
        self.regulatory_gff.is_none() && self.emars_gff.is_none() && self.motif_gff.is_none()
    }
}

/// What the regulatory build wrote.
#[derive(Debug, Default)]
pub struct RegulatoryStats {
    pub regulatory_features: usize,
    pub motif_features: usize,
    pub matrices: usize,
    pub shards: usize,
    pub chromosomes: BTreeSet<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Regulatory,
    Motif,
}

/// A shard element in the shape `vep-cli`'s regulatory loader reads
/// (`regulatory/source.rs`), the shape `info.json` declares as
/// [`vep_core::cache_info::REGULATORY_SCHEMA`]: a regulatory feature with an
/// empty activity string, the cache listing no cell types; a motif naming its
/// matrix by identifier, with neither a database id, a score nor a sequence.
/// Both carry the file they came from (`source`, the input's slot: the
/// regulatory file 0, the EMAR file 1, the motif file 2) and the record's
/// place in that file (`ordinal`), the two keys Ensembl VEP's GFF source
/// gives a record (one identifier set per file, `RegFeat.pm:485`; records
/// met in file order, `:312-328`), which the shard's own order cannot carry
/// once a widened promoter is binned by its widened start.
#[derive(Serialize)]
#[serde(tag = "feature_type")]
enum ShardFeature {
    RegulatoryFeature {
        stable_id: String,
        start: u64,
        end: u64,
        biotype: String,
        activity: String,
        source: u8,
        ordinal: u32,
    },
    MotifFeature {
        stable_id: String,
        start: u64,
        end: u64,
        strand: i8,
        matrix: String,
        source: u8,
        ordinal: u32,
    },
}

impl ShardFeature {
    fn start(&self) -> u64 {
        match self {
            ShardFeature::RegulatoryFeature { start, .. }
            | ShardFeature::MotifFeature { start, .. } => *start,
        }
    }
}

/// One entry of `binding_matrices.json`: a matrix known by its length and
/// factors alone, with no `frequencies`, the stub Ensembl VEP builds for a GFF
/// motif (`RegFeat.pm:440-444`).
#[derive(Debug, Serialize, PartialEq, Eq)]
struct MatrixStub {
    length: u64,
    transcription_factors: Vec<String>,
}

/// The features read so far, by chromosome, and the matrices they name.
#[derive(Default)]
struct Features {
    by_chromosome: BTreeMap<String, Vec<ShardFeature>>,
    matrices: BTreeMap<String, MatrixStub>,
    regulatory_features: usize,
    motif_features: usize,
}

/// Reads the files `inputs` names and writes the cache's regulatory tree under
/// `output_dir`: the shards, by chromosome and 1 Mb region of each feature's
/// start, and the matrices table. A motif whose matrix identifier an earlier
/// motif gave another span or factor list refuses the build, since the table
/// holds one entry per identifier.
pub fn build_regulatory(inputs: &RegulatoryInputs, output_dir: &Path) -> Result<RegulatoryStats> {
    let mut features = Features::default();
    for (source, (kind, path)) in [
        (Kind::Regulatory, &inputs.regulatory_gff),
        (Kind::Regulatory, &inputs.emars_gff),
        (Kind::Motif, &inputs.motif_gff),
    ]
    .into_iter()
    .enumerate()
    {
        if let Some(path) = path {
            info!("Reading regulatory GFF: {}", path.display());
            read_file(
                path,
                kind,
                source as u8,
                inputs.extended_promoters,
                &mut features,
            )
            .with_context(|| format!("Failed to read {}", path.display()))?;
        }
    }

    let root = output_dir.join("regulatory");
    let mut stats = RegulatoryStats {
        regulatory_features: features.regulatory_features,
        motif_features: features.motif_features,
        matrices: features.matrices.len(),
        ..Default::default()
    };
    for (chr, chromosome_features) in &features.by_chromosome {
        let chr_dir = root.join(chr);
        if !is_cache_path_component(chr) || chr_dir.parent() != Some(root.as_path()) {
            bail!("invalid chromosome name '{chr}'");
        }
        std::fs::create_dir_all(&chr_dir)
            .with_context(|| format!("Failed to create {}", chr_dir.display()))?;
        let mut by_region: BTreeMap<u64, Vec<&ShardFeature>> = BTreeMap::new();
        for feature in chromosome_features {
            by_region
                .entry(region_start(feature.start()))
                .or_default()
                .push(feature);
        }
        for (region_start, shard) in &by_region {
            let path = chr_dir.join(format!(
                "{}-{}.json",
                region_start,
                region_start + REGION_SIZE - 1
            ));
            let json = serde_json::to_string_pretty(shard)
                .context("Failed to serialize regulatory features")?;
            std::fs::write(&path, json)
                .with_context(|| format!("Failed to write {}", path.display()))?;
            stats.shards += 1;
        }
        stats.chromosomes.insert(chr.clone());
        info!(
            "chr {}: {} regulatory and motif features in {} regions",
            chr,
            chromosome_features.len(),
            by_region.len()
        );
    }
    std::fs::create_dir_all(&root)
        .with_context(|| format!("Failed to create {}", root.display()))?;
    let table = root.join("binding_matrices.json");
    let json = serde_json::to_string_pretty(&features.matrices)
        .context("Failed to serialize the binding matrices table")?;
    std::fs::write(&table, json).with_context(|| format!("Failed to write {}", table.display()))?;
    info!(
        "Regulatory tree written: {} regulatory features, {} motif features naming {} matrices, {} shards",
        stats.regulatory_features, stats.motif_features, stats.matrices, stats.shards
    );
    Ok(stats)
}

/// Reads one file's records into `features` as `_record_to_regfeat` builds
/// them (`RegFeat.pm:349-385`): a type the kind does not accept is warned about
/// once and skipped (`:314-323`), its start and end read only as
/// [`crate::gff3::check_skipped_coordinates`] reads a skipped line's; a
/// promoter is widened to its extended bounds
/// under `extended_promoters` (`:360-365`), a bound that is not a whole number
/// ignored with a warning; a motif's strand is the column's, `0` when it
/// holds none (`:369-370`), and a regulatory feature carries none, as the
/// shard shape has none; the identifier is `ID`, else `Name`, else the
/// chromosome as the file spells it with the coordinates and type
/// (`:575-584`); a record the file holds twice is kept once (`:628-640`), and
/// a distinct record under a used identifier is kept with the collision warned
/// about once per identifier (`:478-497`). Each feature carries `source` and
/// its `ordinal` in the file. The file is gzip or bgzip compressed, or plain.
/// A byte-order mark at the start of a line is taken off before the line is
/// read ([`crate::gff3::without_byte_order_marks`]).
fn read_file(
    path: &Path,
    kind: Kind,
    source: u8,
    extended_promoters: bool,
    features: &mut Features,
) -> Result<()> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let reader: Box<dyn std::io::BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let accepted: &[&str] = match kind {
        Kind::Regulatory => &REGULATORY_FEATURE_TYPES,
        Kind::Motif => &MOTIF_FEATURE_TYPES,
    };
    let file_name = path.display().to_string();
    let mut warned_types: HashSet<String> = HashSet::new();
    let mut seen = SeenRecords::default();
    // The identifiers this file has yielded, and the ones warned about: a
    // distinct record under a used identifier is a collision (`RegFeat.pm:487`),
    // warned about once per identifier and file.
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut warned_ids: HashSet<String> = HashSet::new();
    let mut ordinal: u32 = 0;
    for (index, line) in BoundedLines::new(reader, MAX_LINE_BYTES, file_name.clone()).enumerate() {
        let line = crate::gff3::without_byte_order_marks(line?);
        let text = line.strip_suffix(b"\r").unwrap_or(&line);
        if text.is_empty() || text[0] == b'#' {
            continue;
        }
        let line_number = index + 1;
        let text = std::str::from_utf8(text)
            .with_context(|| format!("line {line_number} of {file_name} is not UTF-8"))?;
        let columns: Vec<&str> = text.split('\t').collect();
        if columns.len() < 9 {
            bail!(
                "line {line_number} of {file_name} has {} columns where a GFF record has 9",
                columns.len()
            );
        }
        let coordinate = |name: &str, value: &str| -> Result<u64> {
            value.parse().with_context(|| {
                format!(
                    "the {name} of line {line_number} of {file_name} is not a whole number: {value}"
                )
            })
        };
        let feature_type = decode_type(columns[2].as_bytes());
        let feature_type = feature_type.as_ref();
        if !accepted.contains(&feature_type) {
            // A skipped line's start and end are held to the one rule the
            // transcript readers hold them to.
            crate::gff3::check_skipped_coordinates(text, coordinate)?;
            if warned_types.insert(feature_type.to_string()) {
                warn!("{}", unsupported_type_warning(feature_type, &file_name));
            }
            continue;
        }
        let chr = normalize_chr(columns[0]);
        if !is_cache_path_component(chr) {
            bail!(
                "invalid chromosome name '{}' at line {line_number}",
                columns[0]
            );
        }
        let mut start = coordinate("start", columns[3])?;
        let mut end = coordinate("end", columns[4])?;
        if start == 0 || end < start {
            bail!(
                "GFF coordinates {start}-{end} at line {line_number} are not 1-based and ordered"
            );
        }
        // Checked before any widening, on the start column the file is
        // sorted by.
        if !seen.insert(chr, start, text) {
            continue;
        }
        // The hash holds each value as written, the text Ensembl VEP's tests
        // read; a value is decoded where it is read (`decode_attribute`).
        let attributes = crate::gff3::parse_attributes(columns[8]);
        crate::gff3::check_attributes(columns[8], &attributes)
            .with_context(|| format!("line {line_number} of {file_name}"))?;
        let attribute = |key: &str| attributes.get(key).map(String::as_str);
        let decoded = |key: &str| attribute(key).map(|value| decode_attribute(value.as_bytes()));
        let mut ignored_bounds = Vec::new();
        if extended_promoters && EXTENDED_BOUNDS_TYPES.contains(&feature_type) {
            let bound = |key: &str| decoded(key).map(|value| (value.parse::<u64>(), value));
            match bound("extended_start") {
                Some((Ok(extended), _)) => start = start.min(extended),
                Some((Err(_), value)) => ignored_bounds.push(("extended_start", value)),
                None => {}
            }
            match bound("extended_end") {
                Some((Ok(extended), _)) => end = end.max(extended),
                Some((Err(_), value)) => ignored_bounds.push(("extended_end", value)),
                None => {}
            }
        }
        // A coordinate-derived identifier carries the file's own spelling of
        // the chromosome, the name the run-time reader gives the record
        // (`AnnotationSource.pm:179,198`, `RegFeat.pm:372,583`); the shard
        // directory alone takes the cache's. The identifier is chosen on the
        // values as written, and the attribute value chosen is decoded after;
        // the coordinate form is the file's text as the run-time reader
        // writes it.
        let stable_id = match feature_id(
            attribute("ID"),
            attribute("Name"),
            columns[0],
            start,
            end,
            feature_type,
        ) {
            Cow::Borrowed(value) => decode_attribute(value.as_bytes()).into_owned(),
            Cow::Owned(coordinates) => coordinates,
        };
        for (key, value) in ignored_bounds {
            warn!(
                "{}",
                ignored_bound_warning(key, &value, feature_type, &stable_id, &file_name)
            );
        }
        if !seen_ids.insert(stable_id.clone()) && warned_ids.insert(stable_id.clone()) {
            warn!("{}", duplicate_identifier_warning(&stable_id, &file_name));
        }
        let feature = match kind {
            Kind::Regulatory => {
                features.regulatory_features += 1;
                ShardFeature::RegulatoryFeature {
                    stable_id,
                    start,
                    end,
                    biotype: feature_type.to_string(),
                    activity: String::new(),
                    source,
                    ordinal,
                }
            }
            Kind::Motif => {
                let matrix = decoded("binding_matrix_id")
                    .unwrap_or_default()
                    .into_owned();
                let stub = MatrixStub {
                    length: end - start + 1,
                    transcription_factors: transcription_factors(
                        attribute("transcription_factor").map(str::as_bytes),
                    )
                    .into_iter()
                    .map(Cow::into_owned)
                    .collect(),
                };
                match features.matrices.get(&matrix) {
                    Some(earlier) if *earlier != stub => bail!(
                        "binding matrix {matrix} spans {} bases with factors {:?} at line {line_number} of {file_name} where an earlier motif gave it {} bases with factors {:?}; the table holds one entry per matrix",
                        stub.length,
                        stub.transcription_factors,
                        earlier.length,
                        earlier.transcription_factors
                    ),
                    Some(_) => {}
                    None => {
                        features.matrices.insert(matrix.clone(), stub);
                    }
                }
                features.motif_features += 1;
                ShardFeature::MotifFeature {
                    stable_id,
                    start,
                    end,
                    strand: match columns[6] {
                        "+" => 1,
                        "-" => -1,
                        _ => 0,
                    },
                    matrix,
                    source,
                    ordinal,
                }
            }
        };
        features
            .by_chromosome
            .entry(chr.to_string())
            .or_default()
            .push(feature);
        ordinal += 1;
    }
    Ok(())
}

/// The raw records of one chromosome and start position read so far: a
/// tabix-indexed file is sorted by the start column within a chromosome, so a
/// record the file holds twice meets its earlier copy here.
#[derive(Default)]
struct SeenRecords {
    chr: String,
    start: u64,
    records: HashSet<String>,
}

impl SeenRecords {
    /// Whether `record` at `start` on `chr` is new.
    fn insert(&mut self, chr: &str, start: u64, record: &str) -> bool {
        if (self.chr.as_str(), self.start) != (chr, start) {
            self.chr = chr.to_string();
            self.start = start;
            self.records.clear();
        }
        self.records.insert(record.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const FIXTURES: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../vep-cli/tests/fixtures/regulatory-gff"
    );

    fn write_gff(dir: &Path, name: &str, lines: &[impl AsRef<[u8]>]) -> PathBuf {
        let path = dir.join(name);
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        for line in lines {
            encoder.write_all(line.as_ref()).unwrap();
            encoder.write_all(b"\n").unwrap();
        }
        encoder.finish().unwrap();
        path
    }

    fn shard(dir: &Path, chr: &str, name: &str) -> Vec<serde_json::Value> {
        serde_json::from_slice(&std::fs::read(dir.join("regulatory").join(chr).join(name)).unwrap())
            .unwrap()
    }

    /// The three fixture files write the shard layout the runtime reads: one
    /// shard per 1 Mb region of chromosome 21 holding the features starting
    /// there, each regulatory feature with its type as `biotype` and an empty
    /// activity string, each motif naming its matrix, every feature with its
    /// input's slot as `source` and its place in the file as `ordinal`, and a
    /// table of 50 stubs without counts; `--extended-promoters` widens the
    /// promoters' spans.
    #[test]
    fn gff_inputs_write_the_shard_layout() {
        let fixtures = Path::new(FIXTURES);
        let inputs = RegulatoryInputs {
            regulatory_gff: Some(fixtures.join("regulatory_features.gff3.gz")),
            emars_gff: Some(fixtures.join("EMARs.gff.gz")),
            motif_gff: Some(fixtures.join("motif_features.gff3.gz")),
            extended_promoters: false,
        };
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(&inputs, out.path()).unwrap();
        assert_eq!(
            (
                stats.regulatory_features,
                stats.motif_features,
                stats.matrices
            ),
            (516 + 279, 2200, 50)
        );
        assert_eq!(stats.chromosomes.iter().collect::<Vec<_>>(), ["21"]);
        let mut shards: Vec<String> = std::fs::read_dir(out.path().join("regulatory/21"))
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        shards.sort();
        assert_eq!(
            shards,
            [
                "10000001-11000000.json",
                "5000001-6000000.json",
                "6000001-7000000.json",
                "7000001-8000000.json",
                "8000001-9000000.json",
                "9000001-10000000.json"
            ]
        );
        assert_eq!(stats.shards, 6);
        let first = shard(out.path(), "21", "5000001-6000000.json");
        assert!(first.iter().all(|feature| {
            let start = feature["start"].as_u64().unwrap();
            (5_000_001..=6_000_000).contains(&start)
        }));
        let enhancer = first
            .iter()
            .find(|feature| feature["stable_id"] == "ENSR21_85ZF5")
            .unwrap();
        assert_eq!(
            enhancer,
            &serde_json::json!({
                "feature_type": "RegulatoryFeature",
                "stable_id": "ENSR21_85ZF5",
                "start": 5011320,
                "end": 5011550,
                "biotype": "enhancer",
                "activity": "",
                "source": 0,
                "ordinal": 0
            })
        );
        let emar = first
            .iter()
            .find(|feature| feature["stable_id"] == "ENSR00001793149")
            .unwrap();
        assert_eq!(
            (emar["biotype"].as_str(), emar["source"].as_u64()),
            (Some("EMAR"), Some(1))
        );
        let motif = first
            .iter()
            .find(|feature| feature["stable_id"] == "ENSM00000574652")
            .unwrap();
        assert_eq!(
            motif,
            &serde_json::json!({
                "feature_type": "MotifFeature",
                "stable_id": "ENSM00000574652",
                "start": 5010250,
                "end": 5010263,
                "strand": -1,
                "matrix": "ENSPFM0085",
                "source": 2,
                "ordinal": 0
            })
        );
        let table: BTreeMap<String, serde_json::Value> = serde_json::from_slice(
            &std::fs::read(out.path().join("regulatory/binding_matrices.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(table.len(), 50);
        assert_eq!(
            table["ENSPFM0085"],
            serde_json::json!({"length": 14, "transcription_factors": ["EBF1"]})
        );
        assert_eq!(
            table["ENSPFM0290"]["transcription_factors"],
            serde_json::json!(["HNF4A", "NR2F1", "RXRA"])
        );
        assert!(table
            .values()
            .all(|entry| entry.get("frequencies").is_none()));

        let widened = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                extended_promoters: true,
                ..inputs
            },
            widened.path(),
        )
        .unwrap();
        assert_eq!(stats.regulatory_features, 516 + 279);
        let promoter = |dir: &Path| {
            shard(dir, "21", "5000001-6000000.json")
                .into_iter()
                .find(|feature| feature["stable_id"] == "ENSR21_B34QG")
                .map(|feature| {
                    (
                        feature["start"].as_u64().unwrap(),
                        feature["end"].as_u64().unwrap(),
                    )
                })
                .unwrap()
        };
        assert_eq!(promoter(out.path()), (5_553_169, 5_553_746));
        assert_eq!(promoter(widened.path()), (5_553_154, 5_553_746));
    }

    /// A record of an unsupported type is skipped, a `chr`-prefixed name is
    /// written under the cache's spelling, a record without `ID` or `Name` is
    /// named by the file's own spelling and its coordinates, as the run-time
    /// reader names it, a duplicate record is written once while two distinct
    /// records sharing an `ID` are both written in file order (the runtime
    /// collapses them per variant), an `ID` holding a second `=` is cut
    /// there, and a motif naming a matrix with another span refuses the
    /// build.
    #[test]
    fn records_read_as_the_runtime_reads_them() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "reg.gff3.gz",
            &[
                "##gff-version 3",
                "chr21\tEnsembl\texon\t500\t600\t.\t.\t.\tID=skipped",
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tcolor=#faca00",
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tcolor=#faca00",
                "chr21\tEnsembl\tCTCF_binding_site\t2000\t2020\t.\t-\t.\tName=CTCF_1",
                "chr21\tEnsembl\tenhancer\t4000\t4200\t.\t.\t.\tID=TWICE",
                "chr21\tEnsembl\tpromoter\t4150\t4300\t.\t.\t.\tID=TWICE",
                "chr21\tEnsembl\tenhancer\t5000\t5200\t.\t.\t.\tID=CUT=HERE",
            ],
        );
        let motifs = write_gff(
            dir.path(),
            "motifs.gff3.gz",
            &[
                "chr21\tEnsembl\tTF_binding_site\t3000\t3013\t5.1\t+\t.\tID=M1;binding_matrix_id=PFM;transcription_factor=A,B",
                "chr21\tEnsembl\tTF_binding_site\t3100\t3113\t5.1\t-\t.\tID=M2;binding_matrix_id=PFM;transcription_factor=A,B",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                emars_gff: None,
                motif_gff: Some(motifs),
                extended_promoters: false,
            },
            out.path(),
        )
        .unwrap();
        assert_eq!(
            (
                stats.regulatory_features,
                stats.motif_features,
                stats.shards
            ),
            (5, 2, 1)
        );
        let ids: Vec<String> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        // `ID=CUT=HERE` names the feature `CUT`, the text between the first
        // and the second `=`, as the run-time reader and GFF3.pm:194 cut it.
        assert_eq!(
            ids,
            [
                "chr21:1000-1200:enhancer",
                "CTCF_1",
                "TWICE",
                "TWICE",
                "CUT",
                "M1",
                "M2"
            ]
        );

        // The duplicate check keys on the start column the file is sorted by,
        // so a promoter widened under --extended-promoters still meets its copy.
        let widened = write_gff(
            dir.path(),
            "widened.gff3.gz",
            &[
                "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=P;extended_start=900",
                "21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=E",
                "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=P;extended_start=900",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(widened),
                extended_promoters: true,
                ..Default::default()
            },
            out.path(),
        )
        .unwrap();
        assert_eq!(stats.regulatory_features, 2);
        let spans: Vec<(String, u64)> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| {
                (
                    feature["stable_id"].as_str().unwrap().to_string(),
                    feature["start"].as_u64().unwrap(),
                )
            })
            .collect();
        assert_eq!(spans, [("P".to_string(), 900), ("E".to_string(), 1000)]);

        let conflicting = write_gff(
            dir.path(),
            "conflicting.gff3.gz",
            &[
                "21\tEnsembl\tTF_binding_site\t3000\t3013\t5.1\t+\t.\tID=M1;binding_matrix_id=PFM;transcription_factor=A",
                "21\tEnsembl\tTF_binding_site\t3100\t3114\t5.1\t-\t.\tID=M2;binding_matrix_id=PFM;transcription_factor=A",
            ],
        );
        let refused = tempfile::tempdir().unwrap();
        let error = build_regulatory(
            &RegulatoryInputs {
                motif_gff: Some(conflicting),
                ..Default::default()
            },
            refused.path(),
        )
        .unwrap_err();
        assert!(
            format!("{error:#}")
                .contains("binding matrix PFM spans 15 bases with factors [\"A\"] at line 2"),
            "{error:#}"
        );
    }

    /// The type column is read as `get_type` reads it (`GFF3.pm:161-164`,
    /// `uri_unescape`), so a `%XX` escape in it is decoded before the
    /// include-list check (`RegFeat.pm:312-314`): `CTCF%5Fbinding%5Fsite` is
    /// a `CTCF_binding_site` feature with that biotype, `TF%5Fbinding%5Fsite`
    /// a motif, a `promot%65r` is widened and named by the decoded type, and
    /// a type unsupported once decoded is skipped.
    #[test]
    fn an_escaped_type_column_is_decoded_before_the_include_list_check() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "escaped.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=PLAIN_TYPE",
                "21\tEnsembl\tex%6Fn\t2000\t2100\t.\t.\t.\tID=exon1",
                "21\tEnsembl\tCTCF%5Fbinding%5Fsite\t3000\t3020\t.\t+\t.\tID=ESCAPED_TYPE",
                "21\tEnsembl\tpromot%65r\t4000\t4200\t.\t.\t.\textended_start=3900",
                "21\tEnsembl\tpromoter%5Fflanking%5Fregion\t5000\t5100\t.\t.\t.\tcolor=#faca00",
            ],
        );
        let motifs = write_gff(
            dir.path(),
            "escaped_motifs.gff3.gz",
            &["21\tEnsembl\tTF%5Fbinding%5Fsite\t6000\t6013\t.\t-\t.\tID=ESCAPED_MOTIF;binding_matrix_id=M;transcription_factor=TF"],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                emars_gff: None,
                motif_gff: Some(motifs),
                extended_promoters: true,
            },
            out.path(),
        )
        .unwrap();
        assert_eq!((stats.regulatory_features, stats.motif_features), (4, 1));
        let decoded: Vec<(String, u64, Option<String>)> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| {
                (
                    feature["stable_id"].as_str().unwrap().to_string(),
                    feature["start"].as_u64().unwrap(),
                    feature["biotype"].as_str().map(str::to_string),
                )
            })
            .collect();
        assert_eq!(
            decoded,
            [
                ("PLAIN_TYPE".to_string(), 1000, Some("enhancer".to_string())),
                (
                    "ESCAPED_TYPE".to_string(),
                    3000,
                    Some("CTCF_binding_site".to_string())
                ),
                (
                    "21:3900-4200:promoter".to_string(),
                    3900,
                    Some("promoter".to_string())
                ),
                (
                    "21:5000-5100:promoter_flanking_region".to_string(),
                    5000,
                    Some("promoter_flanking_region".to_string())
                ),
                ("ESCAPED_MOTIF".to_string(), 6000, None),
            ]
        );
    }

    /// A type column whose escapes decode to a byte that is no UTF-8 is a type
    /// like any other once decoded: it matches no accepted type, so the record
    /// is skipped and the build goes on, as Ensembl VEP skips the record with
    /// its warning (`RegFeat.pm:314-322`); an escape decoding to text that is
    /// UTF-8 is skipped the same way. A line holding such a byte raw refuses
    /// the build, naming the line, before any column is read.
    #[test]
    fn a_type_column_that_is_not_utf8_is_skipped() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "escaped.gff3.gz",
            &[
                "21\tEnsembl\tcaf%FF\t1000\t1100\t.\t.\t.\tID=ESCAPED",
                "21\tEnsembl\tcaf%CA%80\t2000\t2100\t.\t.\t.\tID=UNSUPPORTED",
                "21\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tID=KEPT",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                ..Default::default()
            },
            out.path(),
        )
        .unwrap();
        assert_eq!(stats.regulatory_features, 1);
        let ids: Vec<String> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        assert_eq!(ids, ["KEPT"]);

        for (name, line) in [
            (
                "raw.gff3.gz",
                b"21\tEnsembl\tcaf\xFF\t1000\t1100\t.\t.\t.\tID=RAW".as_slice(),
            ),
            (
                "mixed.gff3.gz",
                b"21\tEnsembl\tmix\xCA%80\t1000\t1100\t.\t.\t.\tID=MIXED",
            ),
            (
                "attribute.gff3.gz",
                b"21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=caf\xCA%80",
            ),
        ] {
            let regulatory = write_gff(dir.path(), name, &[line]);
            let out = tempfile::tempdir().unwrap();
            let error = build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(regulatory),
                    ..Default::default()
                },
                out.path(),
            )
            .unwrap_err();
            let message = format!("{error:#}");
            assert!(
                message.contains(&format!(
                    "line 1 of {} is not UTF-8",
                    dir.path().join(name).display()
                )),
                "{message}"
            );
        }
    }

    /// Every `%XX` escape of a read attribute value is decoded, in either
    /// letter case and whether or not GFF3 reserves the character, after the
    /// tests on the value as written: `ID=ESC%3bAPE%2c1%3E2` names the feature
    /// `ESC;APE,1>2`, an empty `ID` falls to `Name=a%2Gb%2`, which stays
    /// written, a motif's `binding_matrix_id=ENSPFM%30` is the matrix
    /// `ENSPFM0` and its `transcription_factor=A%2cB,C` is split on the comma
    /// written, so the table lists `A,B` and `C`, the shape the run-time reader
    /// builds. A record without `ID` or `Name` is named by the file's own
    /// spelling of the chromosome as written, `chr%41:3000-3200:enhancer`,
    /// since the run-time reader decodes no coordinate-derived identifier. A
    /// value whose escapes decode to bytes that are not UTF-8 refuses the
    /// build, naming the line and the file, on a record of an accepted type
    /// alone: the skipped `exon` line's `ID=caf%FF` refuses nothing.
    #[test]
    fn every_escape_of_a_read_attribute_value_is_decoded_after_the_written_tests() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "escaped.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ESC%3bAPE%2c1%3E2",
                "21\tEnsembl\tpromoter\t2000\t2200\t.\t.\t.\tID=;Name=a%2Gb%2",
                "chr%41\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tcolor=#faca00",
                "21\tEnsembl\texon\t4000\t4100\t.\t.\t.\tID=caf%FF",
            ],
        );
        let motifs = write_gff(
            dir.path(),
            "motifs.gff3.gz",
            &["21\tEnsembl\tTF_binding_site\t5000\t5013\t.\t+\t.\tID=m;binding_matrix_id=ENSPFM%30;transcription_factor=A%2cB,C"],
        );
        let out = tempfile::tempdir().unwrap();
        build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                motif_gff: Some(motifs),
                ..Default::default()
            },
            out.path(),
        )
        .unwrap();
        let ids: Vec<String> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        assert_eq!(ids, ["ESC;APE,1>2", "a%2Gb%2", "m"]);
        assert_eq!(
            shard(out.path(), "%41", "1-1000000.json")[0]["stable_id"],
            "chr%41:3000-3200:enhancer"
        );
        let table: serde_json::Value = serde_json::from_slice(
            &std::fs::read(out.path().join("regulatory/binding_matrices.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(
            table["ENSPFM0"],
            serde_json::json!({"length": 14, "transcription_factors": ["A,B", "C"]})
        );

        let refused = write_gff(
            dir.path(),
            "refused.gff3.gz",
            &["21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=caf%FF"],
        );
        let out = tempfile::tempdir().unwrap();
        let error = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(refused),
                ..Default::default()
            },
            out.path(),
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(
            message.ends_with(&format!(
                "line 1 of {}: GFF3 attribute ID \"caf%FF\" decodes to text that is not UTF-8",
                dir.path().join("refused.gff3.gz").display()
            )),
            "{message}"
        );
    }

    /// A byte-order mark at the start of a line changes nothing to the bake
    /// ([`crate::gff3::without_byte_order_marks`]): a file with a mark before
    /// each of its lines, its comment among them, writes the files the same
    /// lines write without the marks.
    #[test]
    fn a_marked_file_bakes_as_the_file_without_its_marks() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "##gff-version 3",
            "20\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A20",
            "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=A21",
        ];
        let marked = lines.map(|line| format!("\u{feff}{line}"));
        let bake = |name: &str, lines: &[String]| {
            let out = dir.path().join(name);
            let inputs = RegulatoryInputs {
                regulatory_gff: Some(write_gff(dir.path(), &format!("{name}.gff3.gz"), lines)),
                ..Default::default()
            };
            build_regulatory(&inputs, &out).unwrap();
            [
                "20/1-1000000.json",
                "21/1-1000000.json",
                "binding_matrices.json",
            ]
            .map(|file| std::fs::read(out.join("regulatory").join(file)).unwrap())
        };
        assert_eq!(
            bake("marked", &marked),
            bake("plain", &lines.map(String::from))
        );
    }

    /// A start or end no number leads, or an end below zero, refuses the bake
    /// on a line of a type the kind skips, on a file's first sequence and on a
    /// later one, with the text a record of an accepted type draws, the rule
    /// the transcript readers follow
    /// ([`crate::gff3::check_skipped_coordinates`]). Nothing else of a
    /// skipped line is read: a start of 0, an end before its start that is
    /// not below zero, a seqid that could name no cache directory and a start
    /// a number leads without being one, an empty start before its end among
    /// them, refuse nothing there, where each refuses an `enhancer` line.
    #[test]
    fn a_skipped_line_is_refused_for_a_coordinate_the_index_build_fails_alone() {
        let dir = tempfile::tempdir().unwrap();
        // Two sequences, each an enhancer and one more line: line 2 of the
        // file is on the first sequence and line 4 on the later one.
        let bake = |second: &str, fourth: &str| {
            let regulatory = write_gff(
                dir.path(),
                "t.gff3.gz",
                &[
                    "20\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A20",
                    second,
                    "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A21",
                    fourth,
                ],
            );
            let out = tempfile::tempdir().unwrap();
            build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(regulatory),
                    ..Default::default()
                },
                out.path(),
            )
            .map(|stats| stats.regulatory_features)
            .map_err(|error| format!("{error:#}"))
        };
        let line = |chr: &str, feature_type: &str, start: &str, end: &str| {
            format!("{chr}\tEnsembl\t{feature_type}\t{start}\t{end}\t.\t.\t.\tID=S{chr}")
        };
        let sound_first = line("20", "exon", "2000", "2100");
        let sound_later = line("21", "exon", "2000", "2100");
        let file = dir.path().join("t.gff3.gz").display().to_string();
        for (start, end, column, value) in [
            ("abc", "2100", "start", "abc"),
            ("2000", "xyz", "end", "xyz"),
            ("", "xyz", "start", ""),
            ("2000", "", "end", ""),
            ("2000", "-9", "end", "-9"),
        ] {
            let message = bake(&line("20", "exon", start, end), &sound_later).unwrap_err();
            let refusal =
                format!("the {column} of line 2 of {file} is not a whole number: {value}");
            assert!(message.contains(&refusal), "{message}");
            let message = bake(&sound_first, &line("21", "exon", start, end)).unwrap_err();
            let refusal =
                format!("the {column} of line 4 of {file} is not a whole number: {value}");
            assert!(message.contains(&refusal), "{message}");
        }
        for (chr, start, end) in [
            ("20", "0", "2100"),
            ("20", "2100", "2000"),
            ("20", "2000x", "2100"),
            ("20", "-5", "2100"),
            ("20", " 2000", "2100"),
            ("20", "0x7D0", "2100"),
            ("20", "", "2100"),
            ("20", "2000", "-0"),
            ("20/x", "2000", "2100"),
        ] {
            let on_later = chr.replace("20", "21");
            for (second, fourth) in [
                (line(chr, "exon", start, end), sound_later.clone()),
                (sound_first.clone(), line(&on_later, "exon", start, end)),
            ] {
                assert_eq!(bake(&second, &fourth), Ok(2), "{chr:?} {start:?} {end:?}");
            }
        }
        for (chr, start, end, refusal) in [
            (
                "21",
                "2000x",
                "2100",
                format!("the start of line 4 of {file} is not a whole number: 2000x"),
            ),
            (
                "21",
                "",
                "2100",
                format!(
                    "the start of line 4 of {file} is not a whole number: : cannot parse integer \
                     from empty string"
                ),
            ),
            (
                "21",
                "0",
                "2100",
                "GFF coordinates 0-2100 at line 4 are not 1-based and ordered".to_string(),
            ),
            (
                "21",
                "2100",
                "2000",
                "GFF coordinates 2100-2000 at line 4 are not 1-based and ordered".to_string(),
            ),
            (
                "21/x",
                "2000",
                "2100",
                "invalid chromosome name '21/x' at line 4".to_string(),
            ),
        ] {
            let message = bake(&sound_first, &line(chr, "enhancer", start, end)).unwrap_err();
            assert!(message.contains(&refusal), "{message}");
        }
    }
}
