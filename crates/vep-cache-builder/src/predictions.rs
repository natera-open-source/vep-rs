// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SIFT/PolyPhen prediction matrices, keyed by the MD5 hash of the protein
//! sequence, from Ensembl's variation database: its public MySQL server
//! ([`fetch_predictions`]) or the table dumps a release publishes of it
//! ([`read_prediction_dumps`]).
//!
//! ## The server
//!
//! - Server: `ensembldb.ensembl.org:3306` (anonymous, read-only)
//! - US mirror: `useastdb.ensembl.org:3306`
//! - Database: `homo_sapiens_variation_{release}_{assembly_number}`
//! - Tables: `translation_md5`, `protein_function_predictions`, `attrib`, `attrib_type`
//!
//! The lookup serves `homo_sapiens` on GRCh37 and GRCh38 alone
//! ([`prediction_database`]).
//!
//! ## The table dumps
//!
//! A release publishes the same four tables of a species' variation database
//! as `<table>.txt.gz` or `<table>.txt` ([`PREDICTION_TABLES`]), the files an
//! `--auto` installation reads the matrices of any species from. A matrix
//! value is written with the escapes of the other values, so a tab, a line
//! end, a NUL or a backslash inside the blob follows a backslash and a row
//! spans as many lines as its blob holds line ends.
//!
//! ## Matrix format
//!
//! Each matrix blob is a gzip-compressed binary, the `prediction_matrix`
//! column as Ensembl stores it and as `Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix`
//! holds it in its `matrix` key. See [`vep_core::prediction`] for the layout.
//! The cache stores each blob base64-encoded under the object's other keys
//! (`analysis`, `sub_analysis`, `translation_md5`), so a converted Ensembl
//! cache and a native one read through one decoder.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::thread;
use std::time::Duration;

use anyhow::{bail, Context, Result};
use base64::Engine;
use flate2::read::MultiGzDecoder;
use mysql::prelude::Queryable;
use mysql::{Conn, Opts};
use tracing::{debug, info, warn};

pub use crate::builder::{PeptidePredictions, PredictionBlob};
use crate::core_dumps::{read_table, table_path};
use crate::count::{counted, counted_as};
use crate::fasta::Peptide;

/// Default Ensembl public MySQL host.
pub const DEFAULT_MYSQL_HOST: &str = "ensembldb.ensembl.org";

/// Default MySQL port.
pub const DEFAULT_MYSQL_PORT: u16 = 3306;

/// The tables of a variation database the matrices are read from, each a
/// table dump of the release: the analysis names (`attrib` under the
/// `attrib_type` code [`ANALYSIS_TYPE`]), the peptide MD5 of each matrix
/// (`translation_md5`) and the matrices.
pub const PREDICTION_TABLES: [&str; 4] = [
    "attrib_type",
    "attrib",
    "translation_md5",
    "protein_function_predictions",
];

/// The `attrib_type` code of the analyses a matrix is stored under.
const ANALYSIS_TYPE: &str = "prot_func_analysis";

/// The analyses whose matrices a cache carries.
const ANALYSES: [&str; 3] = ["sift", "polyphen_humvar", "polyphen_humdiv"];

/// The longest row of the matrices table read, 64 MiB: the matrix column is a
/// `mediumblob`, under 16 MiB, and its escapes at most double it.
const MAX_ROW_BYTES: usize = 64 << 20;

/// A fetched prediction matrix with its analysis type.
#[derive(Debug)]
struct FetchedMatrix {
    analysis: String,
    /// The raw prediction matrix blob (gzip-compressed) from the database.
    blob: Vec<u8>,
}

/// The variation database the matrices of `species` on `assembly` are fetched
/// from, or the refusal of a pair the lookup does not serve. The name ends in
/// the number of the human assembly name (`homo_sapiens_variation_116_38`),
/// compared without regard to case; another species' database ends in a number
/// its assembly name does not give, so the lookup serves `homo_sapiens` on
/// GRCh37 and GRCh38 alone and the refusal names the route that serves any
/// species. No server is contacted.
pub fn prediction_database(species: &str, release: u32, assembly: &str) -> Result<String> {
    let assembly_number = match (species, assembly.to_uppercase().as_str()) {
        ("homo_sapiens", "GRCH38") => "38",
        ("homo_sapiens", "GRCH37") => "37",
        _ => bail!(
            "--include-predictions fetches the SIFT and PolyPhen matrices of homo_sapiens on \
             GRCh37 or GRCh38 alone, not of {species} on {assembly}; --auto without the flag \
             reads the matrices of any species from the release's table dumps, and the \
             conversion carries those of its Ensembl cache (docs/cache-conversion.md)"
        ),
    };
    Ok(format!(
        "homo_sapiens_variation_{release}_{assembly_number}"
    ))
}

/// Fetch prediction matrices for a set of protein sequences from Ensembl MySQL.
///
/// # Arguments
/// - `peptide_by_transcript`: The protein FASTA's records, each under its
///   header's `transcript:` token.
/// - `species`: Ensembl species name; [`prediction_database`] refuses a
///   species or assembly the lookup does not serve before any connection.
/// - `release`: Ensembl release number (e.g., 115).
/// - `assembly`: Genome assembly (e.g., "GRCh38").
/// - `host`: MySQL host (default: `ensembldb.ensembl.org`).
/// - `port`: MySQL port (default: 3306).
///
/// # Returns
/// The prediction data of each record, under the record's key.
pub fn fetch_predictions(
    peptide_by_transcript: &HashMap<String, Peptide>,
    species: &str,
    release: u32,
    assembly: &str,
    host: &str,
    port: u16,
) -> Result<HashMap<String, PeptidePredictions>> {
    let db_name = prediction_database(species, release, assembly)?;
    if peptide_by_transcript.is_empty() {
        return Ok(HashMap::new());
    }

    let md5_to_transcripts = transcripts_by_md5(peptide_by_transcript);

    info!("Connecting to {host}:{port}/{db_name}");
    let url = format!("mysql://anonymous@{host}:{port}/{db_name}");
    let opts =
        Opts::from_url(&url).with_context(|| format!("Not a MySQL address: {host}:{port}"))?;
    let mut conn =
        Conn::new(opts.clone()).with_context(|| format!("Failed to connect to {host}:{port}"))?;

    let attrib_ids = resolve_analysis_attribs(&mut conn)?;
    info!(
        "Resolved {}: {:?}",
        counted(attrib_ids.len(), "analysis attrib ID"),
        attrib_ids
    );

    if attrib_ids.is_empty() {
        warn!("No SIFT/PolyPhen analysis attribs found in database; predictions will be empty");
        return Ok(HashMap::new());
    }

    let md5_list: Vec<&str> = md5_to_transcripts.keys().map(|s| s.as_str()).collect();
    let matrices = batch_fetch_matrices(&opts, &mut conn, &md5_list, &attrib_ids)?;

    info!(
        "Fetched {} for {}",
        counted_as(matrices.len(), "prediction matrix", "prediction matrices"),
        counted(md5_list.len(), "unique MD5")
    );

    Ok(predictions_by_transcript(&md5_to_transcripts, matrices))
}

/// The transcripts of each distinct peptide, under the MD5 of its sequence,
/// the key Ensembl stores a matrix under.
fn transcripts_by_md5(
    peptide_by_transcript: &HashMap<String, Peptide>,
) -> HashMap<String, Vec<String>> {
    let mut md5_to_transcripts: HashMap<String, Vec<String>> = HashMap::new();
    for (tx_id, peptide) in peptide_by_transcript {
        let md5 = format!("{:x}", md5::compute(&peptide.sequence));
        md5_to_transcripts
            .entry(md5)
            .or_default()
            .push(tx_id.clone());
    }
    info!(
        "Computing predictions for {} ({})",
        counted(md5_to_transcripts.len(), "unique peptide MD5"),
        counted(peptide_by_transcript.len(), "transcript")
    );
    md5_to_transcripts
}

/// The prediction data of each transcript: the matrices found under its
/// peptide's MD5, none for a peptide the database holds no matrix of.
fn predictions_by_transcript(
    md5_to_transcripts: &HashMap<String, Vec<String>>,
    matrices: Vec<(String, FetchedMatrix)>,
) -> HashMap<String, PeptidePredictions> {
    let mut by_md5: HashMap<String, Vec<FetchedMatrix>> = HashMap::new();
    for (md5, matrix) in matrices {
        by_md5.entry(md5).or_default().push(matrix);
    }

    let mut result: HashMap<String, PeptidePredictions> = HashMap::new();
    for (md5, tx_ids) in md5_to_transcripts {
        let matrices = by_md5.get(md5).map(|v| v.as_slice()).unwrap_or(&[]);
        let predictions = build_peptide_predictions(matrices, md5);
        for tx_id in tx_ids {
            result.insert(tx_id.clone(), predictions.clone());
        }
    }
    result
}

/// The prediction data of each record of `peptide_by_transcript`, read from
/// the dumps of [`PREDICTION_TABLES`] in `dir`: the matrices the variation
/// database stores under the MD5 of a record's peptide, for the analyses
/// [`fetch_predictions`] asks the server for. A database that names none of
/// those analyses, or holds no matrix of these peptides, gives records
/// without a matrix.
pub fn read_prediction_dumps(
    peptide_by_transcript: &HashMap<String, Peptide>,
    dir: &Path,
) -> Result<HashMap<String, PeptidePredictions>> {
    if peptide_by_transcript.is_empty() {
        return Ok(HashMap::new());
    }
    let md5_to_transcripts = transcripts_by_md5(peptide_by_transcript);

    let mut analysis_type = None;
    for row in read_table(dir, "attrib_type", 2)? {
        let row = row?;
        if row.field(1) == ANALYSIS_TYPE {
            analysis_type = Some(row.field(0).to_string());
        }
    }
    let mut analyses: HashMap<Vec<u8>, String> = HashMap::new();
    if let Some(analysis_type) = &analysis_type {
        for row in read_table(dir, "attrib", 3)? {
            let row = row?;
            if row.field(1) == analysis_type && ANALYSES.contains(&row.field(2)) {
                analyses.insert(row.field(0).as_bytes().to_vec(), row.field(2).to_string());
            }
        }
    }
    if analyses.is_empty() {
        warn!(
            "No SIFT/PolyPhen analysis is named in the attrib table of {}; predictions will be empty",
            dir.display()
        );
        return Ok(predictions_by_transcript(&md5_to_transcripts, Vec::new()));
    }

    let mut md5_of: HashMap<Vec<u8>, String> = HashMap::new();
    for row in read_table(dir, "translation_md5", 2)? {
        let row = row?;
        if md5_to_transcripts.contains_key(row.field(1)) {
            md5_of.insert(row.field(0).as_bytes().to_vec(), row.field(1).to_string());
        }
    }

    let table = "protein_function_predictions";
    let path = table_path(dir, table).with_context(|| {
        format!(
            "Table dump {table}.txt.gz (or {table}.txt) not found under {}",
            dir.display()
        )
    })?;
    info!("Reading prediction matrices: {}", path.display());
    let mut matrices = Vec::new();
    read_matrix_rows(
        &path,
        |md5_id, analysis_id| md5_of.contains_key(md5_id) && analyses.contains_key(analysis_id),
        |md5_id, analysis_id, blob| {
            matrices.push((
                md5_of[md5_id].clone(),
                FetchedMatrix {
                    analysis: analyses[analysis_id].clone(),
                    blob,
                },
            ));
        },
    )?;
    info!(
        "Read {} for {}",
        counted_as(matrices.len(), "prediction matrix", "prediction matrices"),
        counted(md5_to_transcripts.len(), "unique MD5")
    );
    Ok(predictions_by_transcript(&md5_to_transcripts, matrices))
}

/// Reads the dump of the matrices table at `path`, gzip when its name ends
/// `.gz`: three columns a row, the MD5's identifier, the analysis's
/// identifier and the matrix. `wanted` is asked once a row's two identifiers
/// are read, and `matrix` is given each wanted row whose matrix is not NULL;
/// the matrix of any other row is passed over without being kept. The values
/// are read as bytes with the escapes `LOAD DATA` reads: a backslash before
/// `0`, `b`, `n`, `r`, `t` or `Z` is that control byte, before any other byte
/// that byte itself (a tab, a line end and a backslash inside a matrix are
/// written so), and a value that is `\N` alone is NULL.
fn read_matrix_rows(
    path: &Path,
    wanted: impl Fn(&[u8], &[u8]) -> bool,
    mut matrix: impl FnMut(&[u8], &[u8], Vec<u8>),
) -> Result<()> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };

    // The row under way: its first two values, the matrix, and which column
    // the next byte belongs to.
    let mut ids: [Vec<u8>; 2] = [Vec::new(), Vec::new()];
    let mut blob: Vec<u8> = Vec::new();
    let mut column = 0usize;
    let mut keep = false;
    // A value that so far is the NULL marker and nothing else.
    let mut null = false;
    let mut value_started = false;
    let mut escaped = false;
    let mut row_bytes = 0usize;
    let mut row_number = 1u64;

    loop {
        let buffer = reader
            .fill_buf()
            .with_context(|| format!("Failed to read {}", path.display()))?;
        if buffer.is_empty() {
            break;
        }
        let read = buffer.len();
        let mut at = 0;
        while at < read {
            // A run of plain bytes is taken whole.
            let plain = if escaped {
                0
            } else {
                buffer[at..]
                    .iter()
                    .position(|byte| matches!(byte, b'\\' | b'\t' | b'\n'))
                    .unwrap_or(read - at)
            };
            if plain > 0 {
                if null {
                    // `\N` followed by more is the letter, not the marker.
                    null = false;
                    push_value(&mut ids, &mut blob, column, keep, b"N");
                }
                push_value(&mut ids, &mut blob, column, keep, &buffer[at..at + plain]);
                value_started = true;
                row_bytes += plain;
                at += plain;
                continue;
            }
            let byte = buffer[at];
            at += 1;
            row_bytes += 1;
            if escaped {
                escaped = false;
                if byte == b'N' && !value_started {
                    null = true;
                } else {
                    let decoded = match byte {
                        b'0' => 0,
                        b'b' => 8,
                        b'n' => b'\n',
                        b'r' => b'\r',
                        b't' => b'\t',
                        b'Z' => 26,
                        other => other,
                    };
                    if null {
                        null = false;
                        push_value(&mut ids, &mut blob, column, keep, b"N");
                    }
                    push_value(&mut ids, &mut blob, column, keep, &[decoded]);
                }
                value_started = true;
            } else if byte == b'\\' {
                escaped = true;
            } else {
                // The end of a value: a tab, or the row's line end.
                let was_null = null;
                null = false;
                value_started = false;
                if byte == b'\t' {
                    column += 1;
                    if column == 2 {
                        keep = wanted(&ids[0], &ids[1]);
                    }
                } else {
                    if column != 2 {
                        bail!(
                            "{} row {row_number} has {} columns, the schema has 3",
                            path.display(),
                            column + 1
                        );
                    }
                    if keep && !was_null {
                        matrix(&ids[0], &ids[1], std::mem::take(&mut blob));
                    }
                    ids[0].clear();
                    ids[1].clear();
                    blob.clear();
                    column = 0;
                    keep = false;
                    row_bytes = 0;
                    row_number += 1;
                }
            }
        }
        if row_bytes > MAX_ROW_BYTES {
            bail!(
                "{} row {row_number} is longer than the {MAX_ROW_BYTES}-byte limit for one row",
                path.display()
            );
        }
        reader.consume(read);
    }
    if column != 0 || value_started || escaped {
        bail!(
            "{} ends inside row {row_number}, which has no line end",
            path.display()
        );
    }
    Ok(())
}

/// Adds `bytes` to the value of `column` in the row under way: an identifier
/// always, the matrix only when the row is kept, and nothing past the third
/// column.
fn push_value(ids: &mut [Vec<u8>; 2], blob: &mut Vec<u8>, column: usize, keep: bool, bytes: &[u8]) {
    match column {
        0 | 1 => ids[column].extend_from_slice(bytes),
        2 if keep => blob.extend_from_slice(bytes),
        _ => {}
    }
}

/// Resolve analysis_attrib_id values for sift, polyphen_humvar, polyphen_humdiv.
fn resolve_analysis_attribs(conn: &mut Conn) -> Result<HashMap<String, u32>> {
    let rows: Vec<(u32, String)> = conn
        .query(
            "SELECT a.attrib_id, a.value
             FROM attrib a
             JOIN attrib_type at ON a.attrib_type_id = at.attrib_type_id
             WHERE at.code = 'prot_func_analysis'
             AND a.value IN ('sift', 'polyphen_humvar', 'polyphen_humdiv')",
        )
        .context("Failed to query analysis attribs")?;

    let map: HashMap<String, u32> = rows.into_iter().map(|(id, name)| (name, id)).collect();
    Ok(map)
}

/// The attempts one batch gets before the fetch fails: the public server
/// closes a connection that has streamed matrices for long enough, which
/// surfaces as a protocol error in the middle of a response, and a fresh
/// connection carries on from the same batch.
const BATCH_ATTEMPTS: u32 = 3;

/// The pause before a batch runs again grows by this much per failed attempt.
const PAUSE_PER_ATTEMPT: Duration = Duration::from_secs(5);

/// Batch fetch prediction matrix blobs from the database, one query per chunk
/// of `md5s`, each run through [`fetch_batch`] so a failed chunk is retried on
/// a new connection to `opts`.
///
/// Returns a list of (md5, FetchedMatrix) pairs.
fn batch_fetch_matrices(
    opts: &Opts,
    conn: &mut Conn,
    md5s: &[&str],
    attrib_ids: &HashMap<String, u32>,
) -> Result<Vec<(String, FetchedMatrix)>> {
    let mut results = Vec::new();
    let batch_size = 500;

    let id_to_name: HashMap<u32, String> = attrib_ids
        .iter()
        .map(|(name, id)| (*id, name.clone()))
        .collect();

    let attrib_id_list: Vec<u32> = attrib_ids.values().copied().collect();
    let attrib_in = attrib_id_list
        .iter()
        .map(|id| id.to_string())
        .collect::<Vec<_>>()
        .join(",");

    for (batch_idx, chunk) in md5s.chunks(batch_size).enumerate() {
        if batch_idx > 0 && batch_idx % 10 == 0 {
            debug!(
                "Fetching predictions batch {}/{}",
                batch_idx,
                md5s.len().div_ceil(batch_size)
            );
        }

        let placeholders: Vec<&str> = chunk.iter().map(|_| "?").collect();
        let query = format!(
            "SELECT t.translation_md5, p.analysis_attrib_id, p.prediction_matrix
             FROM translation_md5 t
             JOIN protein_function_predictions p ON p.translation_md5_id = t.translation_md5_id
             WHERE t.translation_md5 IN ({})
             AND p.analysis_attrib_id IN ({})",
            placeholders.join(","),
            attrib_in
        );

        let params: Vec<mysql::Value> = chunk.iter().map(|md5| mysql::Value::from(*md5)).collect();

        let rows = fetch_batch(opts, conn, &query, &params, batch_idx)?;

        for (md5, attrib_id, blob) in rows {
            if let Some(analysis_name) = id_to_name.get(&attrib_id) {
                results.push((
                    md5,
                    FetchedMatrix {
                        analysis: analysis_name.clone(),
                        blob,
                    },
                ));
            }
        }
    }

    Ok(results)
}

/// One batch's rows. The first attempt runs the query on `conn`; every later
/// attempt opens a new connection to `opts`, which replaces `conn`, and runs
/// the query there. A failed connection or a failed query costs one attempt
/// and a pause of [`PAUSE_PER_ATTEMPT`] times the attempts so far; the batch
/// gets [`BATCH_ATTEMPTS`] attempts in all.
fn fetch_batch(
    opts: &Opts,
    conn: &mut Conn,
    query: &str,
    params: &[mysql::Value],
    batch_idx: usize,
) -> Result<Vec<(String, u32, Vec<u8>)>> {
    let mut attempt = 1;
    loop {
        let connected = if attempt == 1 {
            Ok(())
        } else {
            Conn::new(opts.clone()).map(|fresh| *conn = fresh)
        };
        match connected.and_then(|()| conn.exec(query, params.to_vec())) {
            Ok(rows) => return Ok(rows),
            Err(error) if attempt < BATCH_ATTEMPTS => {
                warn!(
                    "Prediction batch {batch_idx} failed ({error}); connecting again for attempt {}",
                    attempt + 1
                );
                thread::sleep(PAUSE_PER_ATTEMPT * attempt);
                attempt += 1;
            }
            Err(error) => {
                return Err(error).with_context(|| {
                    format!(
                        "Failed to fetch prediction batch {batch_idx} in {BATCH_ATTEMPTS} attempts"
                    )
                })
            }
        }
    }
}

/// Convert fetched matrices into a PeptidePredictions struct.
fn build_peptide_predictions(
    matrices: &[FetchedMatrix],
    translation_md5: &str,
) -> PeptidePredictions {
    let mut predictions = PeptidePredictions::default();

    for matrix in matrices {
        // The blob is already gzip-compressed; base64 is for JSON storage.
        let base64_data = base64::engine::general_purpose::STANDARD.encode(&matrix.blob);

        let (analysis, sub_analysis) = match matrix.analysis.split_once('_') {
            Some((analysis, sub_analysis)) => {
                (analysis.to_string(), Some(sub_analysis.to_string()))
            }
            None => (matrix.analysis.clone(), None),
        };
        let blob = PredictionBlob {
            analysis,
            sub_analysis,
            translation_md5: translation_md5.to_string(),
            base64_data,
        };

        match matrix.analysis.as_str() {
            "sift" => predictions.sift = Some(blob),
            "polyphen_humvar" => predictions.polyphen_humvar = Some(blob),
            "polyphen_humdiv" => predictions.polyphen_humdiv = Some(blob),
            _ => {}
        }
    }

    predictions
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The database is named for human on its two assemblies, in any case of
    /// the assembly name. Any other species or assembly is refused naming the
    /// species limit and the conversion route, and the refusal is settled
    /// before a connection: a fetch for another species against a host that
    /// does not resolve returns it, not a connection failure.
    #[test]
    fn the_lookup_serves_human_alone_and_refuses_another_species_by_name() {
        assert_eq!(
            prediction_database("homo_sapiens", 116, "GRCh38").unwrap(),
            "homo_sapiens_variation_116_38"
        );
        assert_eq!(
            prediction_database("homo_sapiens", 116, "grch37").unwrap(),
            "homo_sapiens_variation_116_37"
        );
        const REFUSAL: &str = "--include-predictions fetches the SIFT and PolyPhen matrices of \
             homo_sapiens on GRCh37 or GRCh38 alone, not of mus_musculus on GRCm39; --auto without \
             the flag reads the matrices of any species from the release's table dumps, and the \
             conversion carries those of its Ensembl cache (docs/cache-conversion.md)";
        assert_eq!(
            prediction_database("mus_musculus", 116, "GRCm39")
                .unwrap_err()
                .to_string(),
            REFUSAL
        );
        for (species, assembly) in [
            ("mus_musculus", "GRCh38"),
            ("homo_sapiens", "T2T-CHM13v2.0"),
        ] {
            let refusal = prediction_database(species, 116, assembly)
                .unwrap_err()
                .to_string();
            assert!(
                refusal.contains(&format!("not of {species} on {assembly};")),
                "{refusal}"
            );
        }

        let peptides = HashMap::from([(
            "transcript".to_string(),
            Peptide {
                sequence: "M".to_string(),
                identifier: "protein".to_string(),
                version: None,
            },
        )]);
        let refused = fetch_predictions(
            &peptides,
            "mus_musculus",
            116,
            "GRCm39",
            "no-such-host.invalid",
            DEFAULT_MYSQL_PORT,
        )
        .unwrap_err()
        .to_string();
        assert_eq!(refused, REFUSAL);
    }

    /// The matrix object's keys follow the Ensembl adaptor: the analysis
    /// attribute value `polyphen_humvar` splits into the analysis `polyphen`
    /// and the sub-analysis `humvar`; `sift` has no sub-analysis; the blob is
    /// stored as it was fetched, base64-encoded.
    #[test]
    fn matrix_object_keys_follow_the_ensembl_adaptor() {
        let fetched = vec![
            FetchedMatrix {
                analysis: "polyphen_humvar".into(),
                blob: vec![0x1f, 0x8b, 1, 2],
            },
            FetchedMatrix {
                analysis: "sift".into(),
                blob: vec![0x1f, 0x8b, 3],
            },
        ];
        let p = build_peptide_predictions(&fetched, "d41d8cd98f00b204e9800998ecf8427e");
        let humvar = p.polyphen_humvar.as_ref().unwrap();
        assert_eq!(humvar.analysis, "polyphen");
        assert_eq!(humvar.sub_analysis.as_deref(), Some("humvar"));
        assert_eq!(humvar.translation_md5, "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(humvar.base64_data, "H4sBAg==");
        let sift = p.sift.as_ref().unwrap();
        assert_eq!(sift.analysis, "sift");
        assert_eq!(sift.sub_analysis, None);
        assert!(p.polyphen_humdiv.is_none());
    }

    /// `value` as a table dump writes it: a backslash before a backslash, a
    /// tab and a line end, and `\0` for a NUL.
    fn escaped(value: &[u8]) -> Vec<u8> {
        let mut written = Vec::new();
        for &byte in value {
            match byte {
                0 => written.extend_from_slice(b"\\0"),
                b'\\' | b'\t' | b'\n' => written.extend_from_slice(&[b'\\', byte]),
                _ => written.push(byte),
            }
        }
        written
    }

    /// Writes `dump` as `<dir>/<name>`, gzip-compressed when the name ends
    /// `.gz`.
    fn write_dump(dir: &Path, name: &str, dump: &[u8]) -> std::path::PathBuf {
        let path = dir.join(name);
        if name.ends_with(".gz") {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            std::io::Write::write_all(&mut encoder, dump).unwrap();
            std::fs::write(&path, encoder.finish().unwrap()).unwrap();
        } else {
            std::fs::write(&path, dump).unwrap();
        }
        path
    }

    /// One row `read_matrix_rows` hands over: the two identifiers and the
    /// matrix.
    type MatrixRow = (Vec<u8>, Vec<u8>, Vec<u8>);

    /// The rows `read_matrix_rows` hands over from the dump at `path`, every
    /// row wanted but those of the analysis `9`, with the rows it asked
    /// about.
    fn matrix_rows(path: &Path) -> Result<(Vec<MatrixRow>, Vec<String>)> {
        let asked = std::cell::RefCell::new(Vec::new());
        let mut rows = Vec::new();
        read_matrix_rows(
            path,
            |md5, analysis| {
                asked
                    .borrow_mut()
                    .push(String::from_utf8_lossy(md5).into_owned());
                analysis != b"9"
            },
            |md5, analysis, matrix| rows.push((md5.to_vec(), analysis.to_vec(), matrix)),
        )?;
        Ok((rows, asked.into_inner()))
    }

    /// A matrix holding every byte value, twice over, reads back byte for
    /// byte behind the escapes a dump writes, from a `.txt` file as from a
    /// `.txt.gz` one. A row whose matrix is NULL is passed over, an unwanted
    /// row's matrix is not handed over, an empty matrix is one of no byte,
    /// and a value that opens with the NULL marker and goes on is the
    /// letter and what follows it.
    #[test]
    fn a_matrix_of_every_byte_reads_back_from_a_plain_and_a_gzipped_dump() {
        let every_byte: Vec<u8> = (0..=255).chain((0..=255).rev()).collect();
        let mut dump = b"1\t2\t".to_vec();
        dump.extend_from_slice(&escaped(&every_byte));
        dump.extend_from_slice(b"\n3\t4\t\\N\n5\t9\t");
        dump.extend_from_slice(&escaped(&every_byte));
        dump.extend_from_slice(b"\n6\t7\t\n8\t2\t\\Nab\n");
        // The four escaped bytes take a second byte each time they occur.
        assert_eq!(escaped(&every_byte).len(), every_byte.len() + 8);

        let dir = tempfile::tempdir().unwrap();
        for name in [
            "protein_function_predictions.txt",
            "protein_function_predictions.txt.gz",
        ] {
            let path = write_dump(dir.path(), name, &dump);
            let (rows, asked) = matrix_rows(&path).unwrap();
            assert_eq!(
                rows,
                [
                    (b"1".to_vec(), b"2".to_vec(), every_byte.clone()),
                    (b"6".to_vec(), b"7".to_vec(), Vec::new()),
                    (b"8".to_vec(), b"2".to_vec(), b"Nab".to_vec()),
                ],
                "{name}"
            );
            assert_eq!(asked, ["1", "3", "5", "6", "8"], "{name}");
        }
    }

    /// A backslash that ends one read of the file escapes the byte that
    /// opens the next: the reader takes a file a mebibyte a read, and the
    /// line end behind the backslash at the last byte of the first mebibyte
    /// stays in the matrix.
    #[test]
    fn an_escape_split_between_two_reads_is_decoded() {
        let plain = vec![b'a'; (1 << 20) - 5];
        let mut dump = b"1\t2\t".to_vec();
        dump.extend_from_slice(&plain);
        dump.extend_from_slice(b"\\\nb\n");
        assert_eq!(dump[(1 << 20) - 1], b'\\');
        let dir = tempfile::tempdir().unwrap();
        let path = write_dump(dir.path(), "protein_function_predictions.txt", &dump);

        let (rows, _) = matrix_rows(&path).unwrap();
        let mut matrix = plain;
        matrix.extend_from_slice(b"\nb");
        assert_eq!(rows.len(), 1);
        assert!(rows[0].2 == matrix, "the matrix differs");
    }

    /// A row of two columns or of four, and a file that ends inside a row,
    /// after a value or after a backslash, are refused naming the file and
    /// the row.
    #[test]
    fn a_row_of_another_width_and_a_cut_file_are_refused_naming_the_row() {
        let dir = tempfile::tempdir().unwrap();
        let refusal = |dump: &[u8]| {
            let path = write_dump(dir.path(), "protein_function_predictions.txt", dump);
            let refusal = matrix_rows(&path).unwrap_err().to_string();
            refusal.replace(&path.display().to_string(), "<file>")
        };
        assert_eq!(
            refusal(b"1\t2\t\\0\n3\t4\n"),
            "<file> row 2 has 2 columns, the schema has 3"
        );
        assert_eq!(
            refusal(b"1\t2\tab\tcd\n"),
            "<file> row 1 has 4 columns, the schema has 3"
        );
        for cut in [
            &b"1\t2\tab\n3\t4\tcd"[..],
            b"1\t2\tab\n3\t4\tcd\\",
            b"1\t2\tab\n3",
        ] {
            assert_eq!(
                refusal(cut),
                "<file> ends inside row 2, which has no line end"
            );
        }
    }

    fn peptide(sequence: &str) -> Peptide {
        Peptide {
            sequence: sequence.to_string(),
            identifier: "protein".to_string(),
            version: None,
        }
    }

    /// Writes the four tables of a made-up variation database under `dir`:
    /// `attrib_type` holds `analysis_type` beside another type; `attrib`
    /// names the three analyses of a cache and one other under the type 2,
    /// and `sift` once more under the type 1; `translation_md5` holds the
    /// MD5 of the peptide `MKV` and of one no record has; each analysis has
    /// a matrix of the first MD5, and `sift` one of the second.
    fn write_database(dir: &Path, analysis_type: &str) {
        let held = format!("{:x}", md5::compute("MKV"));
        let other = format!("{:x}", md5::compute("no record"));
        write_dump(
            dir,
            "attrib_type.txt",
            format!("1\tother_type\tOther\t\\N\n2\t{analysis_type}\tAnalysis\t\\N\n").as_bytes(),
        );
        write_dump(
            dir,
            "attrib.txt.gz",
            b"10\t2\tsift\n11\t2\tpolyphen_humvar\n12\t2\tpolyphen_humdiv\n13\t2\tcadd\n\
              14\t1\tsift\n",
        );
        write_dump(
            dir,
            "translation_md5.txt",
            format!("100\t{held}\n101\t{other}\n").as_bytes(),
        );
        let mut matrices = Vec::new();
        for (md5, analysis, matrix) in [
            (100, 10, &b"sift\t\\"[..]),
            (100, 11, b"humvar\n\0"),
            (100, 12, b"humdiv"),
            (100, 13, b"another analysis"),
            (100, 14, b"another type"),
            (101, 10, b"another peptide"),
        ] {
            matrices.extend_from_slice(format!("{md5}\t{analysis}\t").as_bytes());
            matrices.extend_from_slice(&escaped(matrix));
            matrices.push(b'\n');
        }
        write_dump(dir, "protein_function_predictions.txt.gz", &matrices);
    }

    /// The matrices the dumps hold under the MD5 of a peptide arrive under
    /// each transcript of that peptide, the analysis split as
    /// `build_peptide_predictions` splits it and the matrix base64-encoded
    /// as it is stored. A matrix of an analysis outside the three, one of an
    /// analysis named under another type and one of an MD5 no peptide has
    /// are left out, and a transcript whose peptide has no matrix gets a
    /// record without one. A database whose `attrib_type` lacks the code of
    /// the analyses gives every transcript a record without a matrix.
    #[test]
    fn the_dumps_give_each_transcript_the_matrices_of_its_peptide() {
        let peptides = HashMap::from([
            ("first".to_string(), peptide("MKV")),
            ("second".to_string(), peptide("MKV")),
            ("third".to_string(), peptide("MQL")),
        ]);
        let dir = tempfile::tempdir().unwrap();
        write_database(dir.path(), "prot_func_analysis");

        let read = read_prediction_dumps(&peptides, dir.path()).unwrap();
        assert_eq!(read.len(), 3);
        let encoded = |matrix: &[u8]| base64::engine::general_purpose::STANDARD.encode(matrix);
        let held = format!("{:x}", md5::compute("MKV"));
        let fields = |matrix: &Option<PredictionBlob>| {
            let matrix = matrix.as_ref().expect("a matrix");
            (
                matrix.analysis.clone(),
                matrix.sub_analysis.clone(),
                matrix.translation_md5.clone(),
                matrix.base64_data.clone(),
            )
        };
        for transcript in ["first", "second"] {
            let matrices = &read[transcript];
            assert_eq!(
                fields(&matrices.sift),
                ("sift".to_string(), None, held.clone(), encoded(b"sift\t\\")),
                "{transcript}"
            );
            assert_eq!(
                fields(&matrices.polyphen_humvar),
                (
                    "polyphen".to_string(),
                    Some("humvar".to_string()),
                    held.clone(),
                    encoded(b"humvar\n\0")
                ),
                "{transcript}"
            );
            assert_eq!(
                fields(&matrices.polyphen_humdiv),
                (
                    "polyphen".to_string(),
                    Some("humdiv".to_string()),
                    held.clone(),
                    encoded(b"humdiv")
                ),
                "{transcript}"
            );
        }
        let without = |matrices: &PeptidePredictions| {
            matrices.sift.is_none()
                && matrices.polyphen_humvar.is_none()
                && matrices.polyphen_humdiv.is_none()
        };
        assert!(without(&read["third"]));

        let unnamed = tempfile::tempdir().unwrap();
        write_database(unnamed.path(), "another_analysis");
        let read = read_prediction_dumps(&peptides, unnamed.path()).unwrap();
        assert_eq!(read.len(), 3);
        assert!(read.values().all(without));
    }
}
