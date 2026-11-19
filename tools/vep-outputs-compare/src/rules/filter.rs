// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The adjusted F1's pass over `discordant.tsv`: classify the swap pairs, exclude the ones an
//! Ensembl VEP self-contradiction accounts for, run the reference release's one-sided pass
//! and write the rows that stay open.

use super::breakend::{strip_chr, vep_rs_pos};
use super::gvcf::perl_split;
use super::one_sided::{apply_116_2_one_sided_rules, Row};
use super::row::{pair_row, replaced_reference, PairRow};
use super::scope::{ReferenceRecordScope, ScopeFeatures};
use super::{classify, excluded, ExclusionRule, EXCLUSION_REGISTRY};
use crate::input::{for_each_record, for_each_row_at, open_vcf, push_canonical_location};
use memchr::{memchr, memchr_iter};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// The bucket of the 115.2 one-sided pass.
pub const SV_TRANSCRIPT_SELECTION_BUCKET: &str = "sv_transcript_selection";

/// A bucket's statistics as the report writes them, field by field in order.
pub type BucketStats = Vec<(&'static str, u64)>;

/// What the pass returns: the excluded rows per side and every bucket's statistics in
/// report order.
pub struct Exclusions {
    pub excluded_perl: u64,
    pub excluded_rust: u64,
    pub buckets: Vec<(String, BucketStats)>,
}

/// Names the input VCF behind an output file of the given name, or `None` for a file without
/// one.
pub type InputRecords<'a> = &'a dyn Fn(&str) -> Option<PathBuf>;

/// The one-sided collector of the selected reference release.
pub enum Collector<'a> {
    Release115(&'a ScopeFeatures),
    Release116(&'a ReferenceRecordScope),
}

fn read_rows(path: &Path) -> io::Result<(Vec<String>, Vec<Row>)> {
    let mut columns: Option<Vec<String>> = None;
    let mut rows = Vec::new();
    for_each_record(File::open(path)?, |record| {
        let text = String::from_utf8_lossy(record);
        let fields: Vec<String> = text.split('\t').map(str::to_string).collect();
        let Some(columns) = &columns else {
            columns = Some(fields);
            return Ok(());
        };
        let column = |name: &str| -> String {
            columns
                .iter()
                .position(|c| c == name)
                .and_then(|i| fields.get(i))
                .cloned()
                .unwrap_or_default()
        };
        rows.push(Row {
            file_name: column("file_name"),
            source: column("source"),
            location: column("location"),
            allele: column("allele"),
            feature: column("feature"),
            feature_type: column("feature_type"),
            consequence_set: column("consequence_set"),
            fields,
        });
        Ok(())
    })?;
    Ok((columns.unwrap_or_default(), rows))
}

/// The rows of the two outputs behind unmatched pairs, for the row-aware pair rules: per
/// index of `rows`, the distinct rows of its output that carry its tuple. Each output file
/// that carries a tuple of `pairs` is read once more, as the record names of a dropped record
/// are, and nothing else is read; an index one of whose rows has fewer than nine columns is
/// left out, so no row-aware rule reads a pair with such a row.
fn fetch_pair_rows(
    rows: &[Row],
    pairs: &[(usize, usize)],
    perl_dir: &Path,
    rust_dir: &Path,
) -> io::Result<HashMap<usize, Vec<PairRow>>> {
    let mut wanted: HashMap<PathBuf, HashMap<Vec<u8>, usize>> = HashMap::new();
    for (perl_index, rust_index) in pairs {
        for (index, directory) in [(*perl_index, perl_dir), (*rust_index, rust_dir)] {
            let row = &rows[index];
            let key = [
                row.location.as_str(),
                &row.allele,
                &row.feature,
                &row.feature_type,
                &row.consequence_set,
            ]
            .join("\t");
            wanted
                .entry(directory.join(&row.file_name))
                .or_default()
                .insert(key.into_bytes(), index);
        }
    }
    let mut found: HashMap<usize, Vec<PairRow>> = HashMap::new();
    let mut short: HashSet<usize> = HashSet::new();
    for (path, index_of) in &wanted {
        let locations: HashSet<&[u8]> = index_of
            .keys()
            .map(|key| &key[..memchr(b'\t', key).unwrap_or(key.len())])
            .collect();
        for_each_row_at(File::open(path)?, &locations, |key, columns| {
            let Some(index) = index_of.get(key) else {
                return;
            };
            match columns {
                None => {
                    short.insert(*index);
                }
                Some(columns) => {
                    let row = pair_row(&columns);
                    let rows = found.entry(*index).or_default();
                    if !rows.contains(&row) {
                        rows.push(row);
                    }
                }
            }
        })?;
    }
    found.retain(|index, _| !short.contains(index));
    Ok(found)
}

/// Gives each of `rows` that carries no reference bases (`PairRow::reference`) the bases its
/// allele stands in for, read by position from the records of the input VCF `path` (plain or
/// gzipped), as `read_record_references` of `compare_vep_outputs.py` does: a record is asked
/// about a row whose Location starts on one of the record's own positions on its chromosome,
/// a `chr` prefix set aside (`replaced_reference`). A row two records read differently keeps
/// none, and the file is not opened when no row asks.
fn read_record_references(path: &Path, rows: &mut [&mut PairRow]) -> io::Result<()> {
    let mut asking: HashMap<String, HashMap<i128, Vec<usize>>> = HashMap::new();
    let mut location = Vec::new();
    for (at, row) in rows.iter().enumerate() {
        if row.reference.is_some() || !row.location.contains(':') {
            continue;
        }
        location.clear();
        push_canonical_location(&row.location, &mut location);
        let canonical = std::str::from_utf8(&location).expect("utf-8 location");
        let (chrom, coords) = canonical.split_once(':').expect("a Location with `:`");
        let start = coords.split_once('-').map_or(coords, |(start, _)| start);
        if start.is_empty() || !start.bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }
        // A start past `i128` is read as no number.
        if let Ok(start) = start.parse::<i128>() {
            let asked = asking.entry(strip_chr(chrom).to_string()).or_default();
            asked.entry(start).or_default().push(at);
        }
    }
    if asking.is_empty() {
        return Ok(());
    }
    let mut read: Vec<Vec<String>> = vec![Vec::new(); rows.len()];
    for_each_record(open_vcf(path)?, |record| {
        if record.first() == Some(&b'#') {
            return Ok(());
        }
        // The first five columns alone are decoded: the rest of a line names no allele.
        let mut tabs = memchr_iter(b'\t', record);
        if tabs.nth(3).is_none() {
            return Ok(());
        }
        let text = String::from_utf8_lossy(&record[..tabs.next().unwrap_or(record.len())]);
        let mut fields = text.split('\t');
        let mut field = || fields.next().expect("five columns");
        let (chrom, pos, _, reference, alt) = (field(), field(), field(), field(), field());
        let Some(asked) = asking.get(strip_chr(chrom)) else {
            return Ok(());
        };
        let Some(position) = vep_rs_pos(pos) else {
            return Ok(());
        };
        let alternates = perl_split(alt, ',');
        // An insertion's Location starts on the base before the inserted bases, the record's
        // last where the shared bases are dropped.
        let positions = reference.chars().count().max(1) as i128;
        for start in position..position + positions {
            for at in asked.get(&start).into_iter().flatten() {
                let row = &rows[*at];
                let found = replaced_reference(
                    position,
                    reference,
                    &alternates,
                    &row.location,
                    &row.allele,
                );
                if let Some(found) = found.filter(|found| !read[*at].contains(found)) {
                    read[*at].push(found);
                }
            }
        }
        Ok(())
    })?;
    for (row, mut found) in rows.iter_mut().zip(read) {
        if found.len() == 1 {
            row.reference = found.pop();
        }
    }
    Ok(())
}

/// Reads `discordant_path`, applies the pair rules then the collector's one-sided pass, and
/// writes the rows no rule set aside to `open_out` in their original order and layout.
///
/// `row_rules` are the row-aware pair rules of the reference release and `outputs` the
/// reference's and vep-rs's output directory: a pair no other rule classifies or excludes is
/// put to them in their order, the first match taking it, on the rows of the two outputs that
/// carry its tuples, and their buckets follow the pair buckets. Without such a rule, or
/// without `outputs`, no output is read. `input_records` names the input VCF behind an output
/// file of that name (`--input-vcf` or `--input-dir`), or `None` for a file without one: the
/// rows of such a pair whose first column spells no alleles take their record's from it
/// (`read_record_references`).
pub fn filter_intended_divergences(
    discordant_path: &Path,
    open_out: &Path,
    collector: Collector,
    row_rules: &[&ExclusionRule],
    outputs: Option<(&Path, &Path)>,
    input_records: Option<InputRecords>,
) -> io::Result<Exclusions> {
    let (columns, rows) = read_rows(discordant_path)?;
    let mut perl_by_key: HashMap<(&str, &str, &str, &str), Vec<usize>> = HashMap::new();
    let mut rust_by_key: HashMap<(&str, &str, &str, &str), Vec<usize>> = HashMap::new();
    let (mut n_perl_rows, mut n_rust_rows) = (0u64, 0u64);
    for (idx, row) in rows.iter().enumerate() {
        let key = (
            row.location.as_str(),
            row.allele.as_str(),
            row.feature.as_str(),
            row.feature_type.as_str(),
        );
        match row.source.as_str() {
            "missing_in_rust" => {
                perl_by_key.entry(key).or_default().push(idx);
                n_perl_rows += 1;
            }
            "extra_in_rust" => {
                rust_by_key.entry(key).or_default().push(idx);
                n_rust_rows += 1;
            }
            _ => {}
        }
    }
    let (mut excluded_perl, mut excluded_rust) = (0u64, 0u64);
    let mut masked: HashSet<usize> = HashSet::new();
    let mut pair_stats: Vec<(&'static str, u64, u64)> = EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| rule.classifies && !rule.row_aware())
        .map(|rule| (rule.bucket, 0, 0))
        .collect();
    let mut untaken: Vec<(usize, usize)> = Vec::new();
    for (key, perl_idx) in &perl_by_key {
        let Some(rust_idx) = rust_by_key.get(key) else {
            continue;
        };
        for pi in perl_idx {
            let pcsq = &rows[*pi].consequence_set;
            for ri in rust_idx {
                let rcsq = &rows[*ri].consequence_set;
                let bucket = classify(pcsq, rcsq, key.1);
                let is_excluded = excluded(pcsq, rcsq, key.1);
                match bucket {
                    Some(bucket) => {
                        let entry = pair_stats
                            .iter_mut()
                            .find(|(b, _, _)| *b == bucket)
                            .expect("classifying bucket");
                        entry.1 += 1;
                        if is_excluded {
                            entry.2 += 1;
                        }
                    }
                    None if is_excluded => {
                        return Err(io::Error::other(format!(
                            "pair excluded but classified into no bucket: perl={pcsq:?} rust={rcsq:?}"
                        )));
                    }
                    None => untaken.push((*pi, *ri)),
                }
                if is_excluded {
                    excluded_perl += 1;
                    excluded_rust += 1;
                    masked.insert(*pi);
                    masked.insert(*ri);
                }
            }
        }
    }
    // The row-aware rules, on the pairs no rule above took, in the order of the rows in
    // `discordant.tsv`; a row one of them set aside is not paired again.
    let mut row_stats: Vec<(&'static str, u64, u64)> =
        row_rules.iter().map(|rule| (rule.bucket, 0, 0)).collect();
    if !row_rules.is_empty() && !untaken.is_empty() {
        if let Some((perl_dir, rust_dir)) = outputs {
            untaken.sort_unstable();
            let mut fetched = fetch_pair_rows(&rows, &untaken, perl_dir, rust_dir)?;
            if let Some(input_records) = input_records {
                let mut by_file: HashMap<&str, Vec<&mut PairRow>> = HashMap::new();
                for (index, pair_rows) in &mut fetched {
                    let rows_of_file = by_file.entry(&rows[*index].file_name).or_default();
                    rows_of_file.extend(pair_rows.iter_mut());
                }
                for (file_name, rows_of_file) in &mut by_file {
                    if let Some(records) = input_records(file_name) {
                        read_record_references(&records, rows_of_file).map_err(|e| {
                            io::Error::new(e.kind(), format!("{}: {e}", records.display()))
                        })?;
                    }
                }
            }
            for (pi, ri) in &untaken {
                if masked.contains(pi) || masked.contains(ri) {
                    continue;
                }
                let (Some(perl_rows), Some(rust_rows)) = (fetched.get(pi), fetched.get(ri)) else {
                    continue;
                };
                let Some(at) = row_rules
                    .iter()
                    .position(|rule| rule.matches_rows(perl_rows, rust_rows))
                else {
                    continue;
                };
                row_stats[at].1 += 1;
                if row_rules[at].excludes {
                    row_stats[at].2 += 1;
                    excluded_perl += 1;
                    excluded_rust += 1;
                    masked.insert(*pi);
                    masked.insert(*ri);
                }
            }
        }
    }
    pair_stats.extend(row_stats);
    if excluded_perl > n_perl_rows || excluded_rust > n_rust_rows {
        return Err(io::Error::other(format!(
            "intended-divergence exclusions exceed available discordant rows (perl \
             {excluded_perl} > {n_perl_rows} or rust {excluded_rust} > {n_rust_rows}); a key \
             is matching more than one consequence set per side and the cross product is \
             over-counting"
        )));
    }
    let attributed: u64 = pair_stats.iter().map(|(_, _, e)| e).sum();
    if attributed != excluded_perl {
        return Err(io::Error::other(format!(
            "per-bucket excluded counts sum to {attributed} but {excluded_perl} pairs were \
             excluded; the mask is not attributable to the buckets that report it"
        )));
    }
    let mut buckets: Vec<(String, BucketStats)> = pair_stats
        .iter()
        .map(|(b, count, excl)| (b.to_string(), vec![("count", *count), ("excluded", *excl)]))
        .collect();
    match collector {
        Collector::Release115(scope) if !scope.records.is_empty() => {
            let selection = transcript_selection(
                &rows,
                &mut masked,
                scope,
                &mut excluded_perl,
                &mut excluded_rust,
            );
            buckets.push((SV_TRANSCRIPT_SELECTION_BUCKET.to_string(), selection));
        }
        Collector::Release115(_) => {}
        Collector::Release116(scope) => {
            for (bucket, stats) in apply_116_2_one_sided_rules(&rows, &mut masked, scope) {
                excluded_perl += stats.perl;
                excluded_rust += stats.rust;
                buckets.push((
                    bucket.to_string(),
                    vec![
                        ("count", stats.count),
                        ("excluded", stats.excluded),
                        ("records", stats.records),
                        ("perl", stats.perl),
                        ("rust", stats.rust),
                    ],
                ));
            }
        }
    }
    let mut out = BufWriter::new(File::create(open_out)?);
    out.write_all(columns.join("\t").as_bytes())?;
    out.write_all(b"\n")?;
    for (idx, row) in rows.iter().enumerate() {
        if !masked.contains(&idx) {
            let line: Vec<&str> = columns
                .iter()
                .enumerate()
                .map(|(i, _)| row.fields.get(i).map_or("", String::as_str))
                .collect();
            out.write_all(line.join("\t").as_bytes())?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(Exclusions {
        excluded_perl,
        excluded_rust,
        buckets,
    })
}

/// The 115.2 one-sided pass: the structural-variant comparator's transcript-selection class
/// over the records in `batch_dependent_annotation_scope`, under its two guards. A vep-rs-only
/// row on a transcript the reference did not name is set aside only when the reference named
/// at least one transcript for the record and vep-rs named more; a reference-only row on a
/// transcript vep-rs did not name only when vep-rs named at least as many.
fn transcript_selection(
    rows: &[Row],
    masked: &mut HashSet<usize>,
    scope: &ScopeFeatures,
    excluded_perl: &mut u64,
    excluded_rust: &mut u64,
) -> BucketStats {
    let (mut count, mut excluded, mut perl, mut rust) = (0u64, 0u64, 0u64, 0u64);
    for (idx, row) in rows.iter().enumerate() {
        if masked.contains(&idx) {
            continue;
        }
        let Some((perl_feats, rust_feats)) = scope
            .records
            .get(&(row.location.clone(), row.allele.clone()))
        else {
            continue;
        };
        count += 1;
        match row.source.as_str() {
            "extra_in_rust" => {
                if perl_feats.contains(&row.feature)
                    || perl_feats.is_empty()
                    || rust_feats.len() <= perl_feats.len()
                {
                    continue;
                }
                *excluded_rust += 1;
                rust += 1;
            }
            "missing_in_rust" => {
                if rust_feats.contains(&row.feature)
                    || rust_feats.is_empty()
                    || rust_feats.len() < perl_feats.len()
                {
                    continue;
                }
                *excluded_perl += 1;
                perl += 1;
            }
            _ => continue,
        }
        excluded += 1;
        masked.insert(idx);
    }
    vec![
        ("count", count),
        ("excluded", excluded),
        ("perl", perl),
        ("rust", rust),
    ]
}

#[cfg(test)]
mod tests {
    use super::super::Predicate;
    use super::*;
    use crate::TempDir;
    use std::sync::Mutex;

    const COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    const LOCATION: &str = "19:13208046-13208050";
    const FEATURE: &str = "ENST00000635895";

    /// One fabricated row of the default format.
    fn row(name: &str, feature: &str, consequence: &str, cdna: &str, cds: &str) -> String {
        format!(
            "{name}\t{LOCATION}\t-\tENSG1\t{feature}\tTranscript\t{consequence}\t{cdna}\t{cds}\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1\n"
        )
    }

    /// One row of `discordant.tsv` for a tuple of the fabricated outputs.
    fn discordant(source: &str, feature: &str, consequence: &str) -> String {
        format!("output.txt\t{source}\t{LOCATION}\t-\t{feature}\tTranscript\t{consequence}\n")
    }

    struct Outcome {
        exclusions: Exclusions,
        open_rows: Vec<String>,
    }

    impl Outcome {
        fn bucket(&self, name: &str) -> Vec<(&'static str, u64)> {
            let (_, stats) = self
                .exclusions
                .buckets
                .iter()
                .find(|(bucket, _)| bucket == name)
                .expect("bucket");
            stats.clone()
        }
    }

    /// Runs the pass over a fabricated output pair and its discordant rows under the 116.2
    /// collector, no input VCF behind the outputs.
    fn filtered(
        perl: &str,
        rust: &str,
        discordant_rows: &str,
        row_rules: &[&ExclusionRule],
    ) -> Outcome {
        filtered_beside_input(perl, rust, discordant_rows, row_rules, None)
    }

    /// Runs the pass over a fabricated output pair and its discordant rows under the 116.2
    /// collector, the records `input_records` behind the two outputs where given.
    fn filtered_beside_input(
        perl: &str,
        rust: &str,
        discordant_rows: &str,
        row_rules: &[&ExclusionRule],
        input_records: Option<&str>,
    ) -> Outcome {
        let scratch = TempDir::create(&std::env::temp_dir(), "filter_test_").unwrap();
        let input_path = scratch.path.join("input.vcf");
        if let Some(records) = input_records {
            std::fs::write(&input_path, records).unwrap();
        }
        let input_of = |name: &str| (name == "output.txt").then(|| input_path.clone());
        for (side, rows) in [("perl", perl), ("rust", rust)] {
            std::fs::create_dir(scratch.path.join(side)).unwrap();
            std::fs::write(
                scratch.path.join(side).join("output.txt"),
                format!("{COLUMNS}{rows}"),
            )
            .unwrap();
        }
        let discordant_path = scratch.path.join("discordant.tsv");
        std::fs::write(
            &discordant_path,
            format!(
                "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n{discordant_rows}"
            ),
        )
        .unwrap();
        let open_path = scratch.path.join("discordant_open.tsv");
        let scope = ReferenceRecordScope::default();
        let exclusions = filter_intended_divergences(
            &discordant_path,
            &open_path,
            Collector::Release116(&scope),
            row_rules,
            Some((&scratch.path.join("perl"), &scratch.path.join("rust"))),
            input_records.map(|_| &input_of as InputRecords),
        )
        .unwrap();
        let open_rows = std::fs::read_to_string(&open_path)
            .unwrap()
            .lines()
            .skip(1)
            .map(str::to_string)
            .collect();
        Outcome {
            exclusions,
            open_rows,
        }
    }

    fn fabricated(predicate: fn(&PairRow, &PairRow) -> bool) -> ExclusionRule {
        ExclusionRule {
            bucket: "fabricated_row_shape",
            classifies: true,
            excludes: true,
            predicate: Predicate::Rows(predicate),
        }
    }

    #[test]
    fn a_row_aware_rule_is_handed_the_named_columns_of_the_two_rows_of_an_unmatched_pair() {
        static SEEN: Mutex<Vec<(PairRow, PairRow)>> = Mutex::new(Vec::new());
        fn shape(e: &PairRow, a: &PairRow) -> bool {
            SEEN.lock().unwrap().push((e.clone(), a.clone()));
            true
        }
        let name = "19_13208045_GACCTG/G";
        let outcome = filtered(
            &row(name, FEATURE, "stop_retained_variant", "?-7057", "-"),
            &row(name, FEATURE, "stop_lost", "?-7057/7140", "-/7023"),
            &(discordant("missing_in_rust", FEATURE, "stop_retained_variant")
                + &discordant("extra_in_rust", FEATURE, "stop_lost")),
            &[&fabricated(shape)],
        );
        let seen = SEEN.lock().unwrap();
        let [(e, a)] = &seen[..] else {
            panic!("one pair of rows, got {}", seen.len());
        };
        assert_eq!(
            *e,
            PairRow {
                location: LOCATION.to_string(),
                allele: "-".to_string(),
                feature_type: "Transcript".to_string(),
                consequence: "stop_retained_variant".to_string(),
                cdna_position: "?-7057".to_string(),
                cds_position: "-".to_string(),
                indel: Some(true),
                insertion: Some(false),
                reference: Some("ACCTG".to_string()),
            }
        );
        assert_eq!(
            (
                a.consequence.as_str(),
                a.cdna_position.as_str(),
                a.cds_position.as_str(),
                a.indel
            ),
            ("stop_lost", "?-7057/7140", "-/7023", Some(true))
        );
        assert_eq!(
            outcome.bucket("fabricated_row_shape"),
            vec![("count", 1), ("excluded", 1)]
        );
        assert_eq!(
            outcome.exclusions.buckets[5].0, "fabricated_row_shape",
            "after the five classifying buckets, before the one-sided ones"
        );
        assert_eq!(
            (
                outcome.exclusions.excluded_perl,
                outcome.exclusions.excluded_rust
            ),
            (1, 1)
        );
        assert!(outcome.open_rows.is_empty());
    }

    #[test]
    fn only_an_unmatched_pair_no_earlier_rule_takes_reaches_a_row_aware_rule() {
        static CALLS: Mutex<u32> = Mutex::new(0);
        fn shape(_: &PairRow, _: &PairRow) -> bool {
            *CALLS.lock().unwrap() += 1;
            true
        }
        let matched = row(
            "1993671",
            "ENST00000000001",
            "missense_variant",
            "10-12",
            "4-6",
        );
        let outcome = filtered(
            &(matched.clone()
                + &row(
                    "1993671",
                    "ENST00000000002",
                    "stop_retained_variant",
                    "20-22",
                    "14-16",
                )
                + &row("1993671", "ENST00000000004", "intron_variant", "-", "-")),
            &(matched
                + &row("1993671", "ENST00000000003", "stop_lost", "20-22", "14-16")
                + &row(
                    "1993671",
                    "ENST00000000004",
                    "intron_variant,splice_region_variant",
                    "-",
                    "-",
                )),
            &(discordant(
                "missing_in_rust",
                "ENST00000000002",
                "stop_retained_variant",
            ) + &discordant("extra_in_rust", "ENST00000000003", "stop_lost")
                + &discordant("missing_in_rust", "ENST00000000004", "intron_variant")
                + &discordant(
                    "extra_in_rust",
                    "ENST00000000004",
                    "intron_variant,splice_region_variant",
                )),
            &[&fabricated(shape)],
        );
        assert_eq!(
            *CALLS.lock().unwrap(),
            0,
            "a tuple one engine alone writes and a pair an earlier rule takes"
        );
        assert_eq!(
            outcome.bucket("fabricated_row_shape"),
            vec![("count", 0), ("excluded", 0)]
        );
        assert_eq!(
            outcome.bucket("intronic_overcall_swap"),
            vec![("count", 1), ("excluded", 1)]
        );
        assert_eq!(outcome.open_rows.len(), 2);
    }

    #[test]
    fn the_rows_of_one_key_stand_or_fall_together() {
        static SEEN: Mutex<Vec<String>> = Mutex::new(Vec::new());
        fn every(e: &PairRow, _: &PairRow) -> bool {
            SEEN.lock().unwrap().push(e.cdna_position.clone());
            true
        }
        fn all_but_one(e: &PairRow, _: &PairRow) -> bool {
            e.cdna_position != "20-22"
        }
        let perl = row("first", FEATURE, "stop_retained_variant", "10-12", "4-6")
            + &row("second", FEATURE, "stop_retained_variant", "20-22", "14-16")
            + &row("third", FEATURE, "stop_retained_variant", "10-12", "4-6");
        let rust = row("first", FEATURE, "stop_lost", "10-12", "4-6");
        let rows = discordant("missing_in_rust", FEATURE, "stop_retained_variant")
            + &discordant("extra_in_rust", FEATURE, "stop_lost");
        let taken = filtered(&perl, &rust, &rows, &[&fabricated(every)]);
        assert_eq!(
            *SEEN.lock().unwrap(),
            ["10-12", "20-22"],
            "each distinct row of the reference's tuple, once"
        );
        assert_eq!(
            (taken.bucket("fabricated_row_shape"), taken.open_rows.len()),
            (vec![("count", 1), ("excluded", 1)], 0)
        );
        let refused = filtered(&perl, &rust, &rows, &[&fabricated(all_but_one)]);
        assert_eq!(
            (
                refused.bucket("fabricated_row_shape"),
                refused.open_rows.len()
            ),
            (vec![("count", 0), ("excluded", 0)], 2)
        );
    }

    #[test]
    fn a_row_of_fewer_than_nine_columns_is_put_to_no_row_aware_rule() {
        static CALLS: Mutex<u32> = Mutex::new(0);
        fn shape(_: &PairRow, _: &PairRow) -> bool {
            *CALLS.lock().unwrap() += 1;
            true
        }
        let short = |consequence: &str| {
            format!("1993671\t{LOCATION}\t-\tENSG1\t{FEATURE}\tTranscript\t{consequence}\t?-7057\n")
        };
        let outcome = filtered(
            &short("stop_retained_variant"),
            &short("stop_lost"),
            &(discordant("missing_in_rust", FEATURE, "stop_retained_variant")
                + &discordant("extra_in_rust", FEATURE, "stop_lost")),
            &[&fabricated(shape)],
        );
        assert_eq!(*CALLS.lock().unwrap(), 0);
        assert_eq!(
            (
                outcome.bucket("fabricated_row_shape"),
                outcome.open_rows.len()
            ),
            (vec![("count", 0), ("excluded", 0)], 2)
        );
    }

    /// One fabricated row at `location` under `allele`, its first column `name`.
    fn pair_row_at(name: &str, location: &str, allele: &str) -> PairRow {
        pair_row(&[
            name,
            location,
            allele,
            "ENSG1",
            FEATURE,
            "Transcript",
            "intron_variant",
            "-",
            "-",
        ])
    }

    /// The reference bases each of `rows` carries once the records of an input VCF named
    /// `file_name`, of the bytes `content`, are read for them.
    fn references_read(
        file_name: &str,
        content: &[u8],
        mut rows: Vec<PairRow>,
    ) -> Vec<Option<String>> {
        let scratch = TempDir::create(&std::env::temp_dir(), "filter_test_").unwrap();
        let path = scratch.path.join(file_name);
        std::fs::write(&path, content).unwrap();
        let mut asked: Vec<&mut PairRow> = rows.iter_mut().collect();
        read_record_references(&path, &mut asked).unwrap();
        rows.into_iter().map(|row| row.reference).collect()
    }

    fn bases(sequence: &str) -> Option<String> {
        Some(sequence.to_string())
    }

    #[test]
    fn a_row_under_an_identifier_takes_its_reference_bases_from_the_input_record() {
        let records = "##fileformat=VCFv4.2\n\
                       #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                       7\t100999\tid_1\tTA\tGTGTG\t.\t.\t.\n\
                       7\t200998\tid_2\tCTA\tC\t.\t.\t.\n\
                       7\t300999\tid_3\tT\tTGTG\t.\t.\t.\n\
                       7\t400999\tid_4\tTA\tGTGTG\n\
                       7\t500999\tid_5\tTA\tGTG,C\r\n\
                       7\tnot_a_position\tid_6\tTA\tGTGTG\t.\t.\t.\n";
        let rows = vec![
            pair_row_at("id_1", "7:100999-101000", "GTGTG"),
            pair_row_at("id_2", "7:200999-201000", "-"),
            pair_row_at("id_3", "7:300999-301000", "GTG"),
            pair_row_at("id_4", "7:400999-401000", "GTGTG"),
            pair_row_at("id_5", "7:500999-501000", "GTG"),
            pair_row_at("id_5", "7:501000", "GTG"),
            pair_row_at("id_7", "7:600999-601000", "GTGTG"),
            pair_row_at("7_100999_T/TGTGTG", "7:100999-101000", "GTGTG"),
        ];
        assert_eq!(
            references_read("input.vcf", records.as_bytes(), rows),
            [
                bases("TA"),
                bases("TA"),
                bases(""),
                bases("TA"),
                bases("TA"),
                None,
                None,
                bases(""),
            ],
            "a replacement, a deletion behind the record's first base, a pure insertion, a \
             line of five columns, a line that ends in a carriage return and a line feed, a \
             Location no alternate writes, a Location without a record, and a row whose first \
             column spells its own alleles"
        );
    }

    #[test]
    fn a_chr_prefix_is_set_aside_on_the_record_and_on_the_row() {
        let records = "chr7\t100999\tid_1\tTA\tGTGTG\t.\t.\t.\n\
                       7\t200999\tid_2\tTA\tGTGTG\t.\t.\t.\n\
                       17\t300999\tid_3\tTA\tGTGTG\t.\t.\t.\n";
        let rows = vec![
            pair_row_at("id_1", "7:100999-101000", "GTGTG"),
            pair_row_at("id_2", "CHR7:200999-201000", "GTGTG"),
            pair_row_at("id_3", "7:300999-301000", "GTGTG"),
        ];
        assert_eq!(
            references_read("input.vcf", records.as_bytes(), rows),
            [bases("TA"), bases("TA"), None]
        );
    }

    #[test]
    fn a_gzipped_input_is_read_as_a_plain_one() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
        gzip.write_all(b"7\t100999\tid_1\tTA\tGTGTG\t.\t.\t.\n")
            .unwrap();
        let gzipped = gzip.finish().unwrap();
        for file_name in ["input.vcf.gz", "input.vcf.bgz"] {
            let rows = vec![pair_row_at("id_1", "7:100999-101000", "GTGTG")];
            assert_eq!(
                references_read(file_name, &gzipped, rows),
                [bases("TA")],
                "{file_name}"
            );
        }
    }

    #[test]
    fn a_row_two_records_read_differently_keeps_no_reference_bases() {
        let replaced = "7\t100999\tid_1\tTA\tGTG\t.\t.\t.\n";
        let inserted = "7\t100999\tid_2\tT\tTGTG\t.\t.\t.\n";
        let read = |records: String| {
            let rows = vec![pair_row_at("id_1", "7:100999-101000", "GTG")];
            references_read("input.vcf", records.as_bytes(), rows)
        };
        assert_eq!(read(format!("{replaced}{inserted}")), [None]);
        assert_eq!(read(format!("{replaced}{replaced}")), [bases("TA")]);
        assert_eq!(read(inserted.to_string()), [bases("")]);
    }

    #[test]
    fn the_input_is_not_opened_when_no_row_asks() {
        let absent = std::env::temp_dir().join("filter_test_absent_input.vcf");
        let mut spelt = pair_row_at("7_100999_TA/GTGTG", "7:100999-101000", "GTGTG");
        let mut not_a_position = pair_row_at("id_1", "7:x-101000", "GTGTG");
        let mut no_chromosome = pair_row_at("id_1", "100999-101000", "GTGTG");
        let mut rows = vec![&mut spelt, &mut not_a_position, &mut no_chromosome];
        read_record_references(&absent, &mut rows).expect("no file is opened");
        let mut asking = pair_row_at("id_1", "7:100999-101000", "GTGTG");
        assert!(read_record_references(&absent, &mut [&mut asking]).is_err());
    }

    #[test]
    fn the_rows_of_an_unmatched_pair_take_their_reference_bases_from_the_input_of_their_file() {
        static SEEN: Mutex<Vec<(Option<String>, Option<String>)>> = Mutex::new(Vec::new());
        fn shape(e: &PairRow, a: &PairRow) -> bool {
            let seen = (e.reference.clone(), a.reference.clone());
            SEEN.lock().unwrap().push(seen);
            true
        }
        let perl = row("1993671", FEATURE, "stop_retained_variant", "?-7057", "-");
        let rust = row("1993671", FEATURE, "stop_lost", "?-7057", "-");
        let rows = discordant("missing_in_rust", FEATURE, "stop_retained_variant")
            + &discordant("extra_in_rust", FEATURE, "stop_lost");
        let record = "19\t13208045\t1993671\tGACCTG\tG\t.\t.\t.\n";
        filtered_beside_input(&perl, &rust, &rows, &[&fabricated(shape)], Some(record));
        filtered_beside_input(&perl, &rust, &rows, &[&fabricated(shape)], Some(""));
        filtered(&perl, &rust, &rows, &[&fabricated(shape)]);
        assert_eq!(
            *SEEN.lock().unwrap(),
            [(bases("ACCTG"), bases("ACCTG")), (None, None), (None, None)],
            "beside the record, beside an input without it, and without an input"
        );
    }

    #[test]
    fn without_a_row_aware_rule_no_output_is_read_again() {
        fn shape(_: &PairRow, _: &PairRow) -> bool {
            true
        }
        let scratch = TempDir::create(&std::env::temp_dir(), "filter_test_").unwrap();
        let discordant_path = scratch.path.join("discordant.tsv");
        std::fs::write(
            &discordant_path,
            format!(
                "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n{}{}",
                discordant("missing_in_rust", FEATURE, "stop_retained_variant"),
                discordant("extra_in_rust", FEATURE, "stop_lost")
            ),
        )
        .unwrap();
        let open_path = scratch.path.join("discordant_open.tsv");
        let absent = scratch.path.join("absent");
        let scope = ReferenceRecordScope::default();
        let pass = |row_rules: &[&ExclusionRule]| {
            filter_intended_divergences(
                &discordant_path,
                &open_path,
                Collector::Release116(&scope),
                row_rules,
                Some((&absent, &absent)),
                None,
            )
        };
        let exclusions = pass(&[]).expect("no output is opened");
        assert_eq!((exclusions.excluded_perl, exclusions.excluded_rust), (0, 0));
        assert!(exclusions
            .buckets
            .iter()
            .all(|(bucket, _)| bucket != "fabricated_row_shape"));
        assert!(
            pass(&[&fabricated(shape)]).is_err(),
            "a rule opens the outputs"
        );
    }

    #[test]
    fn the_registry_holds_no_row_aware_rule_outside_its_row_pair_rules() {
        let buckets: Vec<&str> = super::super::row_pair_rules()
            .iter()
            .map(|rule| rule.bucket)
            .collect();
        let row_aware: Vec<&str> = EXCLUSION_REGISTRY
            .iter()
            .filter(|rule| rule.row_aware())
            .map(|rule| rule.bucket)
            .collect();
        assert_eq!(buckets, row_aware);
        for rule in EXCLUSION_REGISTRY.iter().filter(|rule| rule.row_aware()) {
            assert!(!rule.matches("stop_retained_variant", "stop_lost", "-"));
            assert!(!rule.matches_rows(&[], &[]));
        }
    }
}
