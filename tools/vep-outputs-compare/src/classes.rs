// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per-consequence-class F1: one Sequence Ontology term's concordance derived from archived
//! reports (the reference's per-term totals, the discordance report's one-sided counts), its
//! pooling across datasets and the floor over terms clearing the tuple threshold.

use crate::csv;
use crate::file_name_of;
use crate::input::KeyPrep;
use crate::pyfmt::py_round_half_up;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

/// The fixed tuple threshold below which a term is reported but kept out of the floor.
pub const PER_CLASS_MIN_TUPLES: u64 = 1000;

/// One term's concordance for one engine: the three independent counts, the rest derived.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PerClassCounts {
    pub term: String,
    pub vep_tuples: u64,
    pub vep_only: u64,
    pub engine_only: u64,
}

impl PerClassCounts {
    /// Refuses more VEP-only tuples than VEP tuples: the discordance report then does not
    /// belong to the reference output it was paired with.
    pub fn new(
        term: &str,
        vep_tuples: u64,
        vep_only: u64,
        engine_only: u64,
    ) -> Result<PerClassCounts, String> {
        if vep_only > vep_tuples {
            return Err(format!(
                "term {term:?}: {vep_only} VEP-only tuples against {vep_tuples} VEP tuples. The \
                 discordance report does not correspond to this Ensembl VEP output: check the run \
                 id, the assembly and the canonical-contigs setting on both."
            ));
        }
        Ok(PerClassCounts {
            term: term.to_string(),
            vep_tuples,
            vep_only,
            engine_only,
        })
    }

    pub fn intersection(&self) -> u64 {
        self.vep_tuples - self.vep_only
    }

    pub fn engine_tuples(&self) -> u64 {
        self.intersection() + self.engine_only
    }

    pub fn precision(&self) -> f64 {
        ratio(self.intersection(), self.engine_tuples())
    }

    /// Intersection over the VEP-side total, 0.0 when that total is empty.
    pub fn recall(&self) -> f64 {
        ratio(self.intersection(), self.vep_tuples)
    }

    /// F1 as `2I / (V + E)`, the form a consumer re-derives from the CSV's integers.
    pub fn f1(&self) -> f64 {
        ratio(
            2 * self.intersection(),
            self.vep_tuples + self.engine_tuples(),
        )
    }

    pub fn f1_str(&self) -> String {
        py_round_half_up(self.f1(), 6)
    }
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

/// A counter that remembers first-insertion order, as Python's `Counter` does.
#[derive(Default, Clone)]
pub struct Counter {
    pub entries: Vec<(String, u64)>,
    index: HashMap<String, usize>,
}

impl Counter {
    /// The count of `key`, entered at zero when the key is new.
    fn slot(&mut self, key: &str) -> &mut u64 {
        let i = match self.index.get(key) {
            Some(&i) => i,
            None => {
                self.index.insert(key.to_string(), self.entries.len());
                self.entries.push((key.to_string(), 0));
                self.entries.len() - 1
            }
        };
        &mut self.entries[i].1
    }

    pub fn add(&mut self, key: &str, n: u64) {
        *self.slot(key) += n;
    }

    /// `counter[key] = n`.
    pub fn set(&mut self, key: &str, n: u64) {
        *self.slot(key) = n;
    }

    pub fn get(&self, key: &str) -> u64 {
        self.index.get(key).map_or(0, |&i| self.entries[i].1)
    }

    pub fn total(&self) -> u64 {
        self.entries.iter().map(|(_, n)| n).sum()
    }

    /// Adds every count of `other`, new keys in `other`'s order.
    pub fn update(&mut self, other: &Counter) {
        for (key, n) in &other.entries {
            self.add(key, *n);
        }
    }

    /// `most_common(limit)`: by count descending, ties in first-insertion order.
    pub fn most_common(&self, limit: usize) -> Vec<(String, u64)> {
        let mut sorted = self.entries.clone();
        sorted.sort_by_key(|entry| std::cmp::Reverse(entry.1));
        sorted.truncate(limit);
        sorted
    }
}

/// Adds the terms of one key line to `counter`; an empty line counts nothing, and a line
/// without five tab-separated fields is refused with its field count.
pub fn add_key_terms(counter: &mut Counter, line: &[u8]) -> Result<(), usize> {
    if line.is_empty() {
        return Ok(());
    }
    let text = String::from_utf8_lossy(line);
    let fields: Vec<&str> = text.split('\t').collect();
    if fields.len() != 5 {
        return Err(fields.len());
    }
    for term in fields[4].split(',') {
        if !term.is_empty() {
            counter.add(term, 1);
        }
    }
    Ok(())
}

/// The refusal of a key line that [`add_key_terms`] counted `fields` fields in.
pub fn key_fields_error(path: &Path, fields: usize) -> io::Error {
    io::Error::other(format!(
        "{}: expected 5 tab-separated key fields, got {fields}",
        path.display()
    ))
}

/// Counts, per term, the tuples of a sorted-unique key file.
pub fn count_terms_in_sorted_keys(path: &Path) -> io::Result<Counter> {
    let mut counter = Counter::default();
    let mut reader = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        add_key_terms(&mut counter, &line).map_err(|fields| key_fields_error(path, fields))?;
    }
    Ok(counter)
}

/// The columns of `needed` that `header` lacks, as `sorted(missing)` lists them.
fn missing_columns<'a>(header: &[String], needed: &[&'a str]) -> Vec<&'a str> {
    let mut missing: Vec<&str> = needed
        .iter()
        .copied()
        .filter(|name| !header.iter().any(|column| column == name))
        .collect();
    missing.sort_unstable();
    missing
}

/// `int(text)` as a count: Python's `ValueError` message where it is not a non-negative
/// integer.
fn parse_count(text: &str) -> io::Result<u64> {
    crate::rules::py_int(text)
        .and_then(|n| u64::try_from(n).ok())
        .ok_or_else(|| {
            io::Error::other(format!("invalid literal for int() with base 10: {text:?}"))
        })
}

/// Counts per-term one-sided tuples in a `discordant.tsv`: (VEP-only, engine-only), keyed on
/// the `source` column; any other source is a schema change and is refused.
pub fn count_terms_in_discordance(path: &Path) -> io::Result<(Counter, Counter)> {
    let (header, rows) = csv::read_table(path, b'\t')?;
    let column = |name: &str| header.iter().position(|c| c == name);
    let (Some(source_col), Some(csq_col)) = (column("source"), column("consequence_set")) else {
        return Err(io::Error::other(format!(
            "{}: discordance report lacks column(s) {}",
            path.display(),
            csv::py_list_repr(&missing_columns(&header, &["source", "consequence_set"]))
        )));
    };
    let (mut vep_only, mut engine_only) = (Counter::default(), Counter::default());
    let mut prep = KeyPrep::default();
    let mut set = Vec::new();
    for row in &rows {
        let source = row.get(source_col).map_or("", String::as_str);
        let target = match source {
            "missing_in_rust" => &mut vep_only,
            "extra_in_rust" => &mut engine_only,
            _ => {
                return Err(io::Error::other(format!(
                "{}: unrecognised source {source:?}; expected 'missing_in_rust' or 'extra_in_rust'",
                path.display()
            )))
            }
        };
        let raw = row.get(csq_col).map_or("", String::as_str);
        set.clear();
        prep.normalized_set(raw, &mut set);
        for term in std::str::from_utf8(&set).expect("utf-8 terms").split(',') {
            if !term.is_empty() {
                target.add(term, 1);
            }
        }
    }
    Ok((vep_only, engine_only))
}

/// One `PerClassCounts` per term seen on any of the three inputs, by descending VEP tuples
/// then term; a term only the engine emits is kept with `vep_tuples` 0.
pub fn derive_per_class_counts(
    vep_totals: &Counter,
    vep_only: &Counter,
    engine_only: &Counter,
) -> Result<Vec<PerClassCounts>, String> {
    let mut terms: Vec<&str> = vep_totals
        .entries
        .iter()
        .chain(&vep_only.entries)
        .chain(&engine_only.entries)
        .map(|(t, _)| t.as_str())
        .collect();
    terms.sort_unstable();
    terms.dedup();
    let mut rows = terms
        .into_iter()
        .map(|term| {
            PerClassCounts::new(
                term,
                vep_totals.get(term),
                vep_only.get(term),
                engine_only.get(term),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    sort_by_tuples(&mut rows);
    Ok(rows)
}

fn sort_by_tuples(rows: &mut [PerClassCounts]) {
    rows.sort_by(|a, b| {
        b.vep_tuples
            .cmp(&a.vep_tuples)
            .then_with(|| a.term.cmp(&b.term))
    });
}

/// Sums the three counts per term across datasets, then re-derives.
pub fn pool_per_class_counts(rows: &[PerClassCounts]) -> Result<Vec<PerClassCounts>, String> {
    let mut totals: Vec<(String, [u64; 3])> = Vec::new();
    let mut index: HashMap<&str, usize> = HashMap::new();
    for row in rows {
        let i = *index.entry(row.term.as_str()).or_insert_with(|| {
            totals.push((row.term.clone(), [0; 3]));
            totals.len() - 1
        });
        totals[i].1[0] += row.vep_tuples;
        totals[i].1[1] += row.vep_only;
        totals[i].1[2] += row.engine_only;
    }
    let mut pooled = totals
        .iter()
        .map(|(term, [v, vo, eo])| PerClassCounts::new(term, *v, *vo, *eo))
        .collect::<Result<Vec<_>, _>>()?;
    sort_by_tuples(&mut pooled);
    Ok(pooled)
}

/// The published floor and the two facts that qualify it.
pub struct PerClassFloor {
    pub f1: Option<String>,
    pub term: Option<String>,
    pub n_eligible: usize,
    pub n_below_threshold: usize,
}

/// The lowest F1 among terms clearing `min_tuples` and the term attaining it, ties broken
/// on the term name.
pub fn per_class_floor(rows: &[PerClassCounts], min_tuples: u64) -> PerClassFloor {
    let eligible: Vec<&PerClassCounts> =
        rows.iter().filter(|r| r.vep_tuples >= min_tuples).collect();
    let below = rows.len() - eligible.len();
    let worst = eligible.iter().min_by(|a, b| {
        let fa: f64 = a.f1_str().parse().expect("six-decimal text");
        let fb: f64 = b.f1_str().parse().expect("six-decimal text");
        fa.total_cmp(&fb).then_with(|| a.term.cmp(&b.term))
    });
    match worst {
        None => PerClassFloor {
            f1: None,
            term: None,
            n_eligible: 0,
            n_below_threshold: below,
        },
        Some(worst) => PerClassFloor {
            f1: Some(worst.f1_str()),
            term: Some(worst.term.clone()),
            n_eligible: eligible.len(),
            n_below_threshold: below,
        },
    }
}

/// F1 with `excluded` one-sided pairs removed from both sides, rendered like `f1`.
pub fn adjusted_f1_str(counts: &PerClassCounts, excluded: (u64, u64)) -> Result<String, String> {
    let (ex_v, ex_e) = excluded;
    if ex_v > counts.vep_only || ex_e > counts.engine_only {
        return Err(format!(
            "term {:?}: per_term.tsv excludes {ex_v}/{ex_e} pairs against {}/{} one-sided \
             tuples; the per-term file does not belong to this run",
            counts.term, counts.vep_only, counts.engine_only
        ));
    }
    let inter = counts.intersection();
    let denom = 2 * inter + (counts.vep_only - ex_v) + (counts.engine_only - ex_e);
    Ok(py_round_half_up(ratio(2 * inter, denom), 6))
}

const VEP_TOTALS_CACHE_PREFIX: &str = "# vep_totals_fingerprint=";

/// Identity of a per-term VEP total: the source files' names, sizes and modification times
/// in nanoseconds plus the contig filter, the only inputs the count depends on.
pub fn vep_totals_fingerprint(
    perl_files: &[std::path::PathBuf],
    canonical_contigs: Option<&[String]>,
) -> io::Result<String> {
    let mut files = Vec::new();
    for path in perl_files {
        let meta = std::fs::metadata(path)?;
        let mtime = meta
            .modified()?
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        files.push(format!("{}:{}:{mtime}", file_name_of(path), meta.len()));
    }
    let contigs = match canonical_contigs {
        Some(contigs) => {
            let mut sorted = contigs.to_vec();
            sorted.sort_unstable();
            sorted.join(",")
        }
        None => "none".to_string(),
    };
    Ok(format!("{}|canonical={contigs}", files.join(";")))
}

/// The cached per-term VEP totals when `path` records exactly `fingerprint`; None on a
/// missing, stale or unreadable cache.
pub fn read_vep_totals_cache(path: &Path, fingerprint: &str) -> Option<Counter> {
    let text = std::fs::read(path).ok()?;
    let text = String::from_utf8_lossy(&text);
    let (first, rest) = text.split_once('\n').unwrap_or((&text, ""));
    if first.trim_end_matches('\r') != format!("{VEP_TOTALS_CACHE_PREFIX}{fingerprint}") {
        return None;
    }
    let mut records: Vec<Vec<String>> = rest
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| csv::split_record(l, ','))
        .collect();
    if records.is_empty() {
        return Some(Counter::default());
    }
    let header = records.remove(0);
    let term = header.iter().position(|c| c == "term")?;
    let tuples = header.iter().position(|c| c == "vep_tuples")?;
    let mut totals = Counter::default();
    for row in records {
        let n = u64::try_from(crate::rules::py_int(row.get(tuples)?)?).ok()?;
        totals.set(row.get(term)?, n);
    }
    Some(totals)
}

/// Writes the cache atomically, terms sorted.
pub fn write_vep_totals_cache(path: &Path, fingerprint: &str, totals: &Counter) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut out = format!("{VEP_TOTALS_CACHE_PREFIX}{fingerprint}\n").into_bytes();
    csv::push_row(&mut out, &["term", "vep_tuples"], ',', "\n");
    let mut terms: Vec<&(String, u64)> = totals.entries.iter().collect();
    terms.sort_by(|a, b| a.0.cmp(&b.0));
    for (term, n) in terms {
        csv::push_row(&mut out, &[term.as_str(), &n.to_string()], ',', "\n");
    }
    let tmp = path.with_file_name(format!("{}.tmp", file_name_of(path)));
    std::fs::write(&tmp, &out)?;
    std::fs::rename(&tmp, path)
}

pub const PER_CLASS_BY_DATASET_COLUMNS: [&str; 14] = [
    "dataset",
    "assembly",
    "engine",
    "sweep_id",
    "binary_md5",
    "term",
    "vep_tuples",
    "vep_only",
    "engine_only",
    "intersection",
    "engine_tuples",
    "precision",
    "recall",
    "f1",
];

pub const PER_CLASS_POOLED_COLUMNS: [&str; 15] = [
    "engine",
    "sweep_id",
    "binary_md5",
    "term",
    "vep_tuples",
    "vep_only",
    "engine_only",
    "intersection",
    "engine_tuples",
    "precision",
    "recall",
    "f1",
    "vep_only_excluded",
    "engine_only_excluded",
    "adj_f1",
];

/// The cells `term` to `recall` of a per-term row, the same in both CSV schemas.
fn per_class_cells(counts: &PerClassCounts) -> [String; 8] {
    [
        counts.term.clone(),
        counts.vep_tuples.to_string(),
        counts.vep_only.to_string(),
        counts.engine_only.to_string(),
        counts.intersection().to_string(),
        counts.engine_tuples().to_string(),
        py_round_half_up(counts.precision(), 6),
        py_round_half_up(counts.recall(), 6),
    ]
}

/// The labels every row of one (dataset, engine) carries.
pub struct DatasetLabels<'a> {
    pub dataset: &'a str,
    pub assembly: &'a str,
    pub engine: &'a str,
    pub sweep_id: &'a str,
    pub binary_md5: &'a str,
}

/// Writes (or appends) one (dataset, engine)'s per-term rows; an append into a file whose
/// header does not match is refused. Returns the rows written.
pub fn write_per_class_by_dataset_csv(
    path: &Path,
    counts: &[PerClassCounts],
    labels: &DatasetLabels,
    append: bool,
) -> io::Result<usize> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let exists = std::fs::metadata(path).is_ok_and(|m| m.is_file() && m.len() > 0);
    let appending = append && exists;
    let mut out = Vec::new();
    if appending {
        let (header, _) = csv::read_table(path, b',')?;
        if header != PER_CLASS_BY_DATASET_COLUMNS {
            return Err(io::Error::other(format!(
                "{}: cannot append, existing header {} does not match {}",
                path.display(),
                csv::py_list_repr(&header.iter().map(String::as_str).collect::<Vec<_>>()),
                csv::py_list_repr(&PER_CLASS_BY_DATASET_COLUMNS)
            )));
        }
    } else {
        csv::push_row(&mut out, &PER_CLASS_BY_DATASET_COLUMNS, ',', csv::CRLF);
    }
    for row in counts {
        let cells = per_class_cells(row);
        let mut fields: Vec<&str> = vec![
            labels.dataset,
            labels.assembly,
            labels.engine,
            labels.sweep_id,
            labels.binary_md5,
        ];
        fields.extend(cells.iter().map(String::as_str));
        let f1 = row.f1_str();
        fields.push(&f1);
        csv::push_row(&mut out, &fields, ',', csv::CRLF);
    }
    if appending {
        let mut file = std::fs::OpenOptions::new().append(true).open(path)?;
        file.write_all(&out)?;
        file.flush()?;
    } else {
        std::fs::write(path, &out)?;
    }
    Ok(counts.len())
}

/// The rows of the per-(dataset, term) CSV grouped by engine, in first-seen engine order,
/// with each engine's (sweep_id, binary_md5) pin; an engine carrying two pins is refused.
#[allow(clippy::type_complexity)]
pub fn read_per_class_by_dataset_csv(
    path: &Path,
) -> io::Result<(
    Vec<(String, Vec<PerClassCounts>)>,
    HashMap<String, (String, String)>,
)> {
    let (header, rows) = csv::read_table(path, b',')?;
    let needed = [
        "engine",
        "term",
        "vep_tuples",
        "vep_only",
        "engine_only",
        "sweep_id",
        "binary_md5",
    ];
    let missing = missing_columns(&header, &needed);
    if !missing.is_empty() {
        return Err(io::Error::other(format!(
            "{}: lacks column(s) {}",
            path.display(),
            csv::py_list_repr(&missing)
        )));
    }
    let col = |row: &Vec<String>, name: &str| -> String {
        header
            .iter()
            .position(|c| c == name)
            .and_then(|i| row.get(i))
            .cloned()
            .unwrap_or_default()
    };
    let mut by_engine: Vec<(String, Vec<PerClassCounts>)> = Vec::new();
    let mut pins: HashMap<String, Vec<(String, String)>> = HashMap::new();
    for row in &rows {
        let engine = col(row, "engine");
        let counts = PerClassCounts::new(
            &col(row, "term"),
            parse_count(&col(row, "vep_tuples"))?,
            parse_count(&col(row, "vep_only"))?,
            parse_count(&col(row, "engine_only"))?,
        )
        .map_err(io::Error::other)?;
        match by_engine.iter_mut().find(|(e, _)| *e == engine) {
            Some((_, list)) => list.push(counts),
            None => by_engine.push((engine.clone(), vec![counts])),
        }
        let pin = (col(row, "sweep_id"), col(row, "binary_md5"));
        let seen = pins.entry(engine).or_default();
        if !seen.contains(&pin) {
            seen.push(pin);
        }
    }
    let mut mixed: Vec<(&String, &Vec<(String, String)>)> =
        pins.iter().filter(|(_, v)| v.len() != 1).collect();
    if !mixed.is_empty() {
        mixed.sort();
        return Err(io::Error::other(format!(
            "{}: rows of one engine carry more than one (sweep_id, binary_md5) pin, so they are \
             not one measurement: {mixed:?}",
            path.display()
        )));
    }
    let pins = pins
        .into_iter()
        .map(|(e, v)| (e, v.into_iter().next().expect("one pin")))
        .collect();
    Ok((by_engine, pins))
}

/// Engines with vep-rs first, then by name.
pub fn engine_order<'a>(engines: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let mut order: Vec<&str> = engines.collect();
    order.sort_by(|a, b| (*a != "vep-rs", *a).cmp(&(*b != "vep-rs", *b)));
    order
}

/// Writes the pooled CSV, one row per (engine, term), vep-rs first and terms by descending
/// VEP tuple count; the rows of the engine `adjusted` names carry the three adjusted columns
/// and every other engine's are empty. Returns the rows written.
pub fn write_per_class_pooled_csv(
    path: &Path,
    pooled: &[(String, Vec<PerClassCounts>)],
    pins: &HashMap<String, (String, String)>,
    adjusted: Option<(&str, &Exclusions)>,
) -> io::Result<usize> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut out = Vec::new();
    csv::push_row(&mut out, &PER_CLASS_POOLED_COLUMNS, ',', csv::CRLF);
    let mut written = 0;
    for engine in engine_order(pooled.iter().map(|(e, _)| e.as_str())) {
        let (sweep_id, binary_md5) = pins.get(engine).expect("pinned engine");
        let rows = &pooled
            .iter()
            .find(|(e, _)| e == engine)
            .expect("pooled engine")
            .1;
        let exclusions = adjusted.filter(|(e, _)| *e == engine).map(|(_, ex)| ex);
        for row in rows {
            let cells = per_class_cells(row);
            let f1 = row.f1_str();
            let (ex_v, ex_e, adj_f1) = match exclusions {
                Some(ex) => {
                    let (ex_v, ex_e) = ex.get(&row.term).copied().unwrap_or((0, 0));
                    (
                        ex_v.to_string(),
                        ex_e.to_string(),
                        adjusted_f1_str(row, (ex_v, ex_e)).map_err(io::Error::other)?,
                    )
                }
                None => (String::new(), String::new(), String::new()),
            };
            let mut fields: Vec<&str> = vec![engine, sweep_id, binary_md5];
            fields.extend(cells.iter().map(String::as_str));
            fields.extend([f1.as_str(), &ex_v, &ex_e, &adj_f1]);
            csv::push_row(&mut out, &fields, ',', csv::CRLF);
            written += 1;
        }
    }
    std::fs::write(path, &out)?;
    Ok(written)
}

/// Per term, the (VEP-only, engine-only) one-sided counts the mask excludes.
pub type Exclusions = HashMap<String, (u64, u64)>;

/// Per term, the (VEP-only, engine-only) one-sided counts the mask excludes, from a per-term
/// TSV with the columns term, vep_only_excluded and engine_only_excluded.
pub fn read_per_term_exclusions(path: &Path) -> io::Result<Exclusions> {
    let (header, rows) = csv::read_table(path, b'\t')?;
    let needed = ["term", "vep_only_excluded", "engine_only_excluded"];
    let missing = missing_columns(&header, &needed);
    if !missing.is_empty() {
        return Err(io::Error::other(format!(
            "{} lacks column(s) {}",
            path.display(),
            csv::py_list_repr(&missing)
        )));
    }
    let mut out = HashMap::new();
    for row in &rows {
        let cell = |name: &str| -> String {
            header
                .iter()
                .position(|c| c == name)
                .and_then(|i| row.get(i))
                .cloned()
                .unwrap_or_default()
        };
        out.insert(
            cell("term"),
            (
                parse_count(&cell("vep_only_excluded"))?,
                parse_count(&cell("engine_only_excluded"))?,
            ),
        );
    }
    Ok(out)
}

#[cfg(test)]
#[path = "classes_per_class_counts_tests.rs"]
mod per_class_counts_tests;
