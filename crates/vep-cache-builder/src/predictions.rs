// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Fetch SIFT/PolyPhen prediction matrices from Ensembl's public MySQL.
//!
//! Connects to Ensembl's read-only variation database to retrieve pre-computed
//! prediction matrices keyed by MD5 hash of the protein sequence.
//!
//! ## Data source
//!
//! - Server: `ensembldb.ensembl.org:3306` (anonymous, read-only)
//! - US mirror: `useastdb.ensembl.org:3306`
//! - Database: `homo_sapiens_variation_{release}_{assembly_number}`
//! - Tables: `translation_md5`, `protein_function_predictions`, `attrib`
//!
//! ## Matrix format
//!
//! Each matrix blob is a gzip-compressed binary. See [`vep_core::prediction`]
//! for the full format specification. The blobs are stored in the JSON cache
//! as base64-encoded strings.

use std::collections::HashMap;

use anyhow::{bail, Context, Result};
use base64::Engine;
use mysql::prelude::Queryable;
use mysql::{Pool, PooledConn};
use tracing::{debug, info, warn};

/// Default Ensembl public MySQL host.
pub const DEFAULT_MYSQL_HOST: &str = "ensembldb.ensembl.org";

/// Default MySQL port.
pub const DEFAULT_MYSQL_PORT: u16 = 3306;

/// A fetched prediction matrix with its analysis type.
#[derive(Debug)]
struct FetchedMatrix {
    analysis: String,
    /// The raw prediction matrix blob (gzip-compressed) from the database.
    blob: Vec<u8>,
}

/// Per-peptide prediction data: up to 3 matrices (sift, polyphen_humvar, polyphen_humdiv).
#[derive(Debug, Default)]
pub struct PeptidePredictions {
    pub sift: Option<PredictionBlob>,
    pub polyphen_humvar: Option<PredictionBlob>,
    pub polyphen_humdiv: Option<PredictionBlob>,
}

/// A prediction matrix blob ready for JSON serialization.
#[derive(Debug)]
pub struct PredictionBlob {
    pub analysis: String,
    pub peptide_length: usize,
    /// Base64-encoded gzip-compressed binary matrix.
    pub base64_data: String,
}

/// Fetch prediction matrices for a set of protein sequences from Ensembl MySQL.
///
/// # Arguments
/// - `peptide_by_transcript`: Map of transcript ID to amino acid sequence.
/// - `release`: Ensembl release number (e.g., 115).
/// - `assembly`: Genome assembly (e.g., "GRCh38").
/// - `host`: MySQL host (default: `ensembldb.ensembl.org`).
/// - `port`: MySQL port (default: 3306).
///
/// # Returns
/// Map of transcript ID to prediction data.
pub fn fetch_predictions(
    peptide_by_transcript: &HashMap<String, String>,
    release: u32,
    assembly: &str,
    host: &str,
    port: u16,
) -> Result<HashMap<String, PeptidePredictions>> {
    if peptide_by_transcript.is_empty() {
        return Ok(HashMap::new());
    }

    let mut md5_to_transcripts: HashMap<String, Vec<String>> = HashMap::new();
    let mut md5_to_peptide_length: HashMap<String, usize> = HashMap::new();
    for (tx_id, peptide) in peptide_by_transcript {
        let md5 = format!("{:x}", md5::compute(peptide));
        md5_to_transcripts
            .entry(md5.clone())
            .or_default()
            .push(tx_id.clone());
        md5_to_peptide_length.insert(md5, peptide.len());
    }

    info!(
        "Computing predictions for {} unique peptide MD5s ({} transcripts)",
        md5_to_transcripts.len(),
        peptide_by_transcript.len()
    );

    let assembly_number = match assembly.to_uppercase().as_str() {
        "GRCH38" => "38",
        "GRCH37" => "37",
        _ => bail!("Unsupported assembly for prediction lookup: {assembly}"),
    };
    let db_name = format!("homo_sapiens_variation_{release}_{assembly_number}");

    info!("Connecting to {host}:{port}/{db_name}");
    let url = format!("mysql://anonymous@{host}:{port}/{db_name}");
    let pool =
        Pool::new(url.as_str()).with_context(|| format!("Failed to connect to {host}:{port}"))?;
    let mut conn = pool
        .get_conn()
        .context("Failed to get MySQL connection from pool")?;

    let attrib_ids = resolve_analysis_attribs(&mut conn)?;
    info!(
        "Resolved {} analysis attrib IDs: {:?}",
        attrib_ids.len(),
        attrib_ids
    );

    if attrib_ids.is_empty() {
        warn!("No SIFT/PolyPhen analysis attribs found in database; predictions will be empty");
        return Ok(HashMap::new());
    }

    let md5_list: Vec<&str> = md5_to_transcripts.keys().map(|s| s.as_str()).collect();
    let matrices = batch_fetch_matrices(&mut conn, &md5_list, &attrib_ids)?;

    info!(
        "Fetched {} prediction matrices for {} unique MD5s",
        matrices.len(),
        md5_list.len()
    );

    let mut by_md5: HashMap<String, Vec<FetchedMatrix>> = HashMap::new();
    for (md5, matrix) in matrices {
        by_md5.entry(md5).or_default().push(matrix);
    }

    let mut result: HashMap<String, PeptidePredictions> = HashMap::new();
    for (md5, tx_ids) in &md5_to_transcripts {
        let pep_len = md5_to_peptide_length.get(md5).copied().unwrap_or(0);
        let predictions = build_peptide_predictions(
            by_md5.get(md5).map(|v| v.as_slice()).unwrap_or(&[]),
            pep_len,
        );
        for tx_id in tx_ids {
            result.insert(tx_id.clone(), predictions.clone());
        }
    }

    Ok(result)
}

impl Clone for PeptidePredictions {
    fn clone(&self) -> Self {
        PeptidePredictions {
            sift: self.sift.as_ref().map(|b| PredictionBlob {
                analysis: b.analysis.clone(),
                peptide_length: b.peptide_length,
                base64_data: b.base64_data.clone(),
            }),
            polyphen_humvar: self.polyphen_humvar.as_ref().map(|b| PredictionBlob {
                analysis: b.analysis.clone(),
                peptide_length: b.peptide_length,
                base64_data: b.base64_data.clone(),
            }),
            polyphen_humdiv: self.polyphen_humdiv.as_ref().map(|b| PredictionBlob {
                analysis: b.analysis.clone(),
                peptide_length: b.peptide_length,
                base64_data: b.base64_data.clone(),
            }),
        }
    }
}

/// Resolve analysis_attrib_id values for sift, polyphen_humvar, polyphen_humdiv.
fn resolve_analysis_attribs(conn: &mut PooledConn) -> Result<HashMap<String, u32>> {
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

/// Batch fetch prediction matrix blobs from the database.
///
/// Returns a list of (md5, FetchedMatrix) pairs.
fn batch_fetch_matrices(
    conn: &mut PooledConn,
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

        let rows: Vec<(String, u32, Vec<u8>)> = conn
            .exec(query, params)
            .with_context(|| format!("Failed to fetch prediction batch {batch_idx}"))?;

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

/// Convert fetched matrices into a PeptidePredictions struct.
fn build_peptide_predictions(
    matrices: &[FetchedMatrix],
    peptide_length: usize,
) -> PeptidePredictions {
    let mut predictions = PeptidePredictions::default();

    for matrix in matrices {
        // The blob is already gzip-compressed; base64 is for JSON storage.
        let base64_data = base64::engine::general_purpose::STANDARD.encode(&matrix.blob);

        let blob = PredictionBlob {
            analysis: matrix.analysis.clone(),
            peptide_length,
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
