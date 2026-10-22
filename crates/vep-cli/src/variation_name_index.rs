// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variant identifiers resolved to their genomic mappings from the cache's
//! variation files, the source `--format id` reads.
//!
//! Ensembl VEP's identifier parser asks its variation database for each name
//! (`Parser/ID.pm:162` at release/116.2, `VariationAdaptor::fetch_by_name`,
//! whose `v.name = ?` compares case-insensitively under the database's
//! collation) and takes every variation feature of the variant it finds. The
//! variation files hold those features by name, indexed by position alone, so
//! the run's names are resolved together: every name of the list is held in
//! memory and every variation file is read once, start to end, keeping the
//! records whose name is in the list. Nothing is written beside the cache.
//!
//! The pass costs a read of the whole variation data, about the size of the
//! files, whatever the number of names; the chromosomes are read in parallel.

use std::io::BufReader;
use std::path::Path;

use anyhow::Context;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use vep_core::cache_info::VariationColumns;
use vep_core::cache_manifest::ManifestError;
use vep_io::input::id::NameMapping;

use crate::variation_source::{parse_record, read_variation_line, VariationSource};

/// The mappings of every name in `names`, one list per name in the names'
/// order, from the variation files of `source`: a name's mappings come in
/// karyotype order of their chromosomes and, within a chromosome, in the file's
/// order (by start); a name the files lack has an empty list. Names are
/// matched without regard to ASCII case, as the database compares them. A
/// record the files repeat is one mapping: the dumper's query joins
/// `failed_variation` (`AnnotationSource/Database/Variation.pm:175`), so a
/// variant failed for two reasons is written twice, identical, where the
/// database holds one feature. An empty list reads no file.
pub fn resolve_names(
    source: &VariationSource,
    names: &[&str],
) -> anyhow::Result<Vec<Vec<NameMapping>>> {
    if names.is_empty() {
        return Ok(Vec::new());
    }
    let keys: Vec<Vec<u8>> = names
        .iter()
        .map(|name| name.as_bytes().to_ascii_lowercase())
        .collect();
    let wanted: FxHashSet<&[u8]> = keys.iter().map(Vec::as_slice).collect();
    let columns = source.columns();
    let name_column = columns
        .position("variation_name")
        .context("the cache's variation_cols name no variation_name column")?;
    let delimiter = if source.is_tabix() { b'\t' } else { b' ' };
    let found: Vec<Vec<(Vec<u8>, NameMapping)>> = source
        .files()
        .into_par_iter()
        .map(|(chr, path)| {
            scan_file(source, chr, path, columns, name_column, delimiter, &wanted).map_err(
                |error| {
                    // The refusal names the file itself.
                    if error.is::<ManifestError>() {
                        error
                    } else {
                        error.context(format!("reading {}", path.display()))
                    }
                },
            )
        })
        .collect::<anyhow::Result<_>>()?;
    let mut by_name: FxHashMap<Vec<u8>, Vec<NameMapping>> = FxHashMap::default();
    for (key, mapping) in found.into_iter().flatten() {
        let mappings = by_name.entry(key).or_default();
        if !mappings.contains(&mapping) {
            mappings.push(mapping);
        }
    }
    Ok(keys
        .iter()
        .map(|key| by_name.get(key).cloned().unwrap_or_default())
        .collect())
}

/// The records of one variation file whose lower-cased name is in `wanted`,
/// each with that key, in file order; `chr` is the file's chromosome. The file
/// is read through `source`, which checks it against the cache's manifest, with
/// the length of a line capped ([`read_variation_line`]).
fn scan_file(
    source: &VariationSource,
    chr: &str,
    path: &Path,
    columns: &VariationColumns,
    name_column: usize,
    delimiter: u8,
    wanted: &FxHashSet<&[u8]>,
) -> anyhow::Result<Vec<(Vec<u8>, NameMapping)>> {
    let whole = source.open_whole(path)?;
    let mut reader = BufReader::new(flate2::read::MultiGzDecoder::new(whole));
    let mut line = Vec::new();
    let mut key = Vec::new();
    let mut found = Vec::new();
    let mut line_number = 0;
    loop {
        line.clear();
        line_number += 1;
        if read_variation_line(&mut reader, &mut line, line_number)? == 0 {
            break;
        }
        while line.last().is_some_and(|&b| b == b'\n' || b == b'\r') {
            line.pop();
        }
        let Some(name) = line.split(|&b| b == delimiter).nth(name_column) else {
            continue;
        };
        key.clear();
        key.extend(name.iter().map(u8::to_ascii_lowercase));
        if !wanted.contains(key.as_slice()) {
            continue;
        }
        let text = String::from_utf8_lossy(&line);
        let fields: Vec<&str> = text.split(delimiter as char).collect();
        if let Some(record) = parse_record(columns, &fields, false) {
            found.push((
                key.clone(),
                NameMapping {
                    chr: chr.to_string(),
                    start: record.start,
                    end: record.end,
                    strand: record.strand,
                    allele_string: record.allele_string,
                },
            ));
        }
    }
    reader.into_inner().into_inner().finish()?;
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use std::sync::Arc;
    use vep_core::cache_info::CacheInfo;
    use vep_core::cache_manifest::CacheManifest;

    const COLS: &str = "chr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms,AF";

    /// The variation files under `cache`, opened for an `info.json` whose
    /// `variation_cols` are `cols` and whose `var_type` is `var_type`.
    fn open_source(cache: &Path, cols: &str, var_type: Option<&str>) -> VariationSource {
        let info = CacheInfo {
            variation_cols: cols.split(',').map(str::to_string).collect(),
            var_type: var_type.map(str::to_string),
            ..Default::default()
        };
        VariationSource::open(cache, &info, false, Arc::new(CacheManifest::absent(cache)))
            .unwrap()
            .expect("variation files")
    }

    /// Writes `lines` gzipped as `<cache>/variations/<chr>/<file>`; a pass
    /// reads the stream, so the index of the start-indexed layout is not needed.
    fn write_gz(cache: &Path, chr: &str, file: &str, lines: &[&str]) {
        let dir = cache.join("variations").join(chr);
        fs::create_dir_all(&dir).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(dir.join(file)).unwrap(),
            flate2::Compression::default(),
        );
        for line in lines {
            writeln!(encoder, "{line}").unwrap();
        }
        encoder.finish().unwrap();
    }

    fn mapping(chr: &str, start: u64, end: u64, strand: i8, alleles: &str) -> NameMapping {
        NameMapping {
            chr: chr.to_string(),
            start,
            end,
            strand,
            allele_string: alleles.to_string(),
        }
    }

    /// Two chromosomes of the start-indexed layout: a name maps once, a name
    /// maps on both chromosomes and comes back in karyotype order (1 before 21),
    /// a name written in upper case resolves as the database's case-insensitive
    /// comparison resolves it, an unknown name has no mapping, and a repeated
    /// name is answered each time.
    #[test]
    fn resolves_names_over_every_file_in_karyotype_order() {
        let tmp = tempfile::tempdir().unwrap();
        write_gz(
            tmp.path(),
            "21",
            "all_vars.gz",
            &[
                "21\trs2838519\t.\t.\t44195140\t.\tG/A/C\t.\t.\t1\t.\t.\t.\t.\t.\tA:0.5625",
                "21\trs77\t.\t.\t44195200\t44195199\t-/T\t-1\t.\t.\t.\t.\t.\t.\t.\t.",
            ],
        );
        write_gz(
            tmp.path(),
            "1",
            "all_vars.gz",
            &["1\trs77\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\t."],
        );
        let source = open_source(tmp.path(), COLS, Some("tabix"));
        let resolved =
            resolve_names(&source, &["rs2838519", "RS77", "rs404", "rs2838519"]).unwrap();
        assert_eq!(
            resolved,
            vec![
                vec![mapping("21", 44195140, 44195140, 1, "G/A/C")],
                vec![
                    mapping("1", 100, 100, 1, "A/G"),
                    mapping("21", 44195200, 44195199, -1, "-/T"),
                ],
                vec![],
                vec![mapping("21", 44195140, 44195140, 1, "G/A/C")],
            ]
        );
    }

    /// rs2516718524 of the release 116 GRCh38 cache, a variant failed for two
    /// reasons and so written twice by the dumper, is one mapping, as the
    /// database holds one feature; a second record of the name at another
    /// position stays a mapping of its own.
    #[test]
    fn a_record_the_files_repeat_is_one_mapping() {
        let tmp = tempfile::tempdir().unwrap();
        let repeated =
            "21\trs2516718524\t1\t.\t33248813\t33248814\tTCT/TCTC\t.\t.\t.\t.\t.\t.\t.\t.\t.";
        write_gz(
            tmp.path(),
            "21",
            "all_vars.gz",
            &[
                repeated,
                repeated,
                "21\trs2516718524\t1\t.\t33248900\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\t.",
            ],
        );
        let source = open_source(tmp.path(), COLS, Some("tabix"));
        assert_eq!(
            resolve_names(&source, &["rs2516718524"]).unwrap(),
            vec![vec![
                mapping("21", 33248813, 33248814, 1, "TCT/TCTC"),
                mapping("21", 33248900, 33248900, 1, "A/G"),
            ]]
        );
    }

    /// The region layout: space-delimited records without a chromosome column
    /// (an empty cell is an empty string between two spaces), read through every
    /// region file of the chromosome in start order.
    #[test]
    fn resolves_names_on_the_region_layout() {
        let tmp = tempfile::tempdir().unwrap();
        let cols = COLS.strip_prefix("chr,").unwrap();
        write_gz(
            tmp.path(),
            "21",
            "2000001-3000000_var.gz",
            &["rs9   2500000  A/T"],
        );
        write_gz(
            tmp.path(),
            "21",
            "1-1000000_var.gz",
            &["rs9   500 501 AC/A"],
        );
        let source = open_source(tmp.path(), cols, None);
        let resolved = resolve_names(&source, &["rs9"]).unwrap();
        assert_eq!(
            resolved,
            vec![vec![
                mapping("21", 500, 501, 1, "AC/A"),
                mapping("21", 2500000, 2500000, 1, "A/T"),
            ]]
        );
    }

    /// A variation file whose second line is one byte over the cap fails the
    /// resolution naming the file and the line, with no mapping resolved from
    /// the file; the same file with the line at the cap resolves both names.
    #[test]
    fn line_past_the_cap_fails_the_resolution() {
        let cap = vep_io::input::MAX_LINE_BYTES;
        let tmp = tempfile::tempdir().unwrap();
        let first = "21\trs1\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\t.";
        let mut second = String::from("21\trs2\t.\t.\t200\t.\tA/G\t.\t.\t");
        second.push_str(&"x".repeat(cap + 1 - second.len()));
        write_gz(tmp.path(), "21", "all_vars.gz", &[first, &second]);
        let source = open_source(tmp.path(), COLS, Some("tabix"));
        let error = resolve_names(&source, &["rs1", "rs2"]).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!(
                "reading {}: line 2 is longer than the {cap}-byte limit for one line",
                tmp.path().join("variations/21/all_vars.gz").display()
            )
        );

        second.pop();
        write_gz(tmp.path(), "21", "all_vars.gz", &[first, &second]);
        let source = open_source(tmp.path(), COLS, Some("tabix"));
        assert_eq!(
            resolve_names(&source, &["rs1", "rs2"]).unwrap(),
            vec![
                vec![mapping("21", 100, 100, 1, "A/G")],
                vec![mapping("21", 200, 200, 1, "A/G")],
            ]
        );
    }

    /// The name column may be the last one, so the line end is removed before
    /// the name is read.
    #[test]
    fn resolves_a_name_in_the_last_column() {
        let tmp = tempfile::tempdir().unwrap();
        write_gz(tmp.path(), "21", "1-1000000_var.gz", &["100  A/G 1 rs1\r"]);
        let source = open_source(
            tmp.path(),
            "start,end,allele_string,strand,variation_name",
            None,
        );
        assert_eq!(
            resolve_names(&source, &["rs1"]).unwrap(),
            vec![vec![mapping("21", 100, 100, 1, "A/G")]]
        );
    }

    /// An empty list reads no file, so an unreadable file is not reached.
    #[test]
    fn an_empty_list_reads_nothing() {
        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path().join("variations").join("21");
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("all_vars.gz"), b"not gzip at all").unwrap();
        let source = open_source(tmp.path(), COLS, Some("tabix"));
        assert!(resolve_names(&source, &[]).unwrap().is_empty());
    }

    /// A variation file that cannot be read ends the pass with the file named.
    #[test]
    fn unreadable_file_names_itself() {
        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path().join("variations").join("21");
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("all_vars.gz"), b"not gzip at all").unwrap();
        let source = open_source(tmp.path(), COLS, Some("tabix"));
        let error = resolve_names(&source, &["rs1"]).unwrap_err();
        assert!(format!("{error:#}").contains("all_vars.gz"), "{error:#}");
    }
}
