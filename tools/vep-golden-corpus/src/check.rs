// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `check` subcommand: a corpus that holds an identifier under restricted redistribution
//! terms is refused, and the combinations, strata and classes the reference outputs show are
//! held against the exemplars the manifest records.

use crate::classes::{genotype_class_names, parse_family_columns, record_genotype_classes};
use crate::cli::CheckArgs;
use crate::constants::FIXED_COLUMNS;
use crate::error::{complain, failure, say, Failure, TOOL};
use crate::identifiers::stored::{corpus_files, restricted_identifier_errors};
use crate::json::{Object, Value};
use crate::manifest::{
    covers_lists, matrix_items, names_outside_matrix, names_under, read_manifest,
};
use crate::names::record_name_index;
use crate::pytext::{path_join, path_str};
use crate::rows::blocks::{lines, read_blocks, thread_count, RowReader};
use crate::rows::{
    kept_samples, parse_dataset_arg, record_class_pass, substitution_alts, vcf_lines, Dataset,
};
use std::collections::{BTreeSet, HashMap, HashSet};

/// The dataset a record of the manifest names.
enum RecordDataset<'a> {
    Name(&'a str),
    /// A value that is no text, which names no dataset.
    NoName,
    /// No `dataset` key at all.
    Absent,
}

/// A record of the manifest as `check` reads it: its dataset and its `vcf` line.
struct CorpusRecord<'a> {
    dataset: RecordDataset<'a>,
    vcf: &'a str,
}

impl CorpusRecord<'_> {
    /// The record's chromosome: the first column of its line.
    fn chromosome(&self) -> &str {
        self.vcf.split('\t').next().unwrap_or_default()
    }
}

/// The records of a manifest: a list of objects that each hold a `vcf` line.
fn corpus_records(manifest: &Object) -> Result<Vec<CorpusRecord<'_>>, Failure> {
    let records = manifest.get("records").and_then(Value::as_list);
    let records = records.ok_or_else(|| Failure::own("the manifest has no list of records"))?;
    records
        .iter()
        .enumerate()
        .map(|(place, record)| {
            let vcf = record.get("vcf").and_then(Value::as_str).ok_or_else(|| {
                Failure::own(format!("record {place} of the manifest has no vcf line"))
            })?;
            let dataset = match record.get("dataset") {
                Some(Value::Str(name)) => RecordDataset::Name(name),
                Some(_) => RecordDataset::NoName,
                None => RecordDataset::Absent,
            };
            Ok(CorpusRecord { dataset, vcf })
        })
        .collect()
}

/// The corpus records by the names they can carry in an output.
struct Corpus<'a> {
    records: Vec<CorpusRecord<'a>>,
    by_name: HashMap<String, Vec<usize>>,
}

impl Corpus<'_> {
    /// Whether a record of the corpus carries `name` in the dataset `dataset` on
    /// `chromosome`. The records of the name are asked in their order; one without a
    /// `dataset` key that is reached is the error, as its place.
    fn holds(&self, dataset: &str, name: &str, chromosome: &str) -> Result<bool, usize> {
        for place in self.by_name.get(name).into_iter().flatten() {
            let record = &self.records[*place];
            match record.dataset {
                RecordDataset::Absent => return Err(*place),
                RecordDataset::Name(held)
                    if held == dataset && record.chromosome() == chromosome =>
                {
                    return Ok(true)
                }
                _ => {}
            }
        }
        Ok(false)
    }
}

fn no_dataset(place: usize) -> Failure {
    Failure::own(format!("record {place} of the manifest has no dataset"))
}

/// The chromosome of a Location: the text before its first `:`.
fn chromosome(location: &str) -> &str {
    location.split_once(':').map_or(location, |(head, _)| head)
}

/// What a worker reads of one block of a reference: the consequence sets its rows show,
/// and with family specs the strata they show and the strata rows of corpus records show.
#[derive(Default)]
struct BlockSets {
    sets: HashSet<String>,
    strata: HashSet<String>,
    covered: HashSet<String>,
    /// The place of a record without a `dataset` key that a row reached.
    without_dataset: Option<usize>,
}

fn block_sets(
    block: &str,
    reader: &RowReader<'_>,
    dataset: &str,
    corpus: Option<&Corpus<'_>>,
) -> BlockSets {
    let mut found = BlockSets::default();
    // The rows of one record follow each other, so the answer for a name on a chromosome
    // is asked once a run of rows.
    let mut last: Option<(&str, &str, Result<bool, usize>)> = None;
    for line in lines(block) {
        let Some(view) = reader.view(line) else {
            continue;
        };
        if !found.sets.contains(view.consequence_set.as_ref()) {
            found.sets.insert(view.consequence_set.to_string());
        }
        let Some(corpus) = corpus else {
            continue;
        };
        let (name, on) = (view.uploaded_variation, chromosome(view.location));
        let held = match last {
            Some((last_name, last_on, held)) if last_name == name && last_on == on => held,
            _ => corpus.holds(dataset, name, on),
        };
        last = Some((name, on, held));
        let stratum = view.family_stratum();
        match held {
            Ok(true) if !found.covered.contains(&stratum) => {
                found.covered.insert(stratum.clone());
            }
            Ok(_) => {}
            Err(place) => {
                found.without_dataset.get_or_insert(place);
            }
        }
        found.strata.insert(stratum);
    }
    found
}

/// What iterating a record's `combinations` yields as names: the texts of a list, the
/// characters of a text, the keys of an object.
fn combinations(record: &Value, place: usize) -> Result<Vec<String>, Failure> {
    let unreadable = || {
        Failure::own(format!(
            "record {place} of the manifest has no combinations that can be read"
        ))
    };
    match record.get("combinations").ok_or_else(unreadable)? {
        Value::List(items) => items
            .iter()
            .filter(|item| {
                !matches!(
                    item,
                    Value::Null | Value::Bool(_) | Value::Int(_) | Value::Float(_)
                )
            })
            .map(|item| item.as_str().map(str::to_string).ok_or_else(unreadable))
            .collect(),
        Value::Str(text) => Ok(text.chars().map(String::from).collect()),
        Value::Object(object) => Ok(object.keys().map(str::to_string).collect()),
        _ => Err(unreadable()),
    }
}

/// One line a way `variants.vcf` fails to carry the manifest's kept samples
/// (`sample_column_errors`): a `#CHROM` line whose columns are not the eight fixed ones,
/// `FORMAT` and the manifest's `samples`, or a record whose column count differs from that
/// header's.
fn sample_column_errors(
    corpus: &str,
    samples: &[String],
    keep: usize,
) -> Result<Vec<String>, Failure> {
    let mut expected: Vec<&str> = FIXED_COLUMNS.to_vec();
    expected.push("FORMAT");
    expected.extend(samples.iter().map(String::as_str));
    let mut errors = Vec::new();
    let mut has_header = false;
    let mut lines = vcf_lines(&path_join(corpus, "variants.vcf"), None);
    while let Some(line) = lines.next_columns() {
        let (is_header, f) = line?;
        if is_header {
            has_header = true;
            if f != expected {
                let named = f.get(8..).unwrap_or_default().join(" ");
                errors.push(format!(
                    "variants.vcf header names {}, the manifest keeps FORMAT {}",
                    if named.is_empty() {
                        "no sample columns"
                    } else {
                        &named
                    },
                    samples.join(" ")
                ));
            }
        } else if f.len() != 9 + keep {
            errors.push(format!(
                "variants.vcf record {} has {} columns, the header {}",
                f[..5].join(" "),
                f.len(),
                9 + keep
            ));
        }
    }
    if !has_header {
        errors.push("variants.vcf has no #CHROM header line".to_string());
    }
    Ok(errors)
}

/// The genotype classes the inputs show among the manifest's kept samples, the ones the
/// corpus records show, and those of the first set or of `required` the corpus lacks
/// (`genotype_coverage`).
fn genotype_coverage(
    datasets: &[Dataset],
    records: &[CorpusRecord<'_>],
    samples: &[String],
    keep: usize,
    required: &[String],
) -> Result<(usize, usize, Vec<String>), Failure> {
    let mut observed: BTreeSet<&str> = BTreeSet::new();
    for (_, _, input) in datasets {
        let mut named = false;
        let mut lines = vcf_lines(input, Some(9 + keep));
        while let Some(line) = lines.next_columns() {
            let (is_header, f) = line?;
            if is_header {
                named = !kept_samples(&f, keep, input, samples, "the manifest")?.is_empty();
            } else {
                observed.extend(record_genotype_classes(&f, keep)?);
            }
        }
        if !named {
            return Err(failure(format!(
                "{input} has no #CHROM header line naming its samples"
            )));
        }
    }
    let mut covered: BTreeSet<&str> = BTreeSet::new();
    for record in records {
        let columns: Vec<&str> = record.vcf.split('\t').collect();
        covered.extend(record_genotype_classes(&columns, keep)?);
    }
    let mut wanted: BTreeSet<&str> = observed.clone();
    wanted.extend(required.iter().map(String::as_str));
    let missing = wanted
        .into_iter()
        .filter(|class| !covered.contains(class))
        .map(str::to_string)
        .collect();
    Ok((observed.len(), covered.len(), missing))
}

/// Writes the errors found so far when the run ends on a refusal of the program's own; a
/// refusal of the tool's own stands alone.
fn after_errors(errors: &[String], refusal: Failure) -> Failure {
    if !refusal.line.starts_with(&format!("ERROR: [{TOOL}]")) {
        for error in errors {
            complain(error);
        }
    }
    refusal
}

/// Runs `check`.
pub fn run(args: &CheckArgs) -> Result<i32, Failure> {
    let corpus_path = path_str(&args.corpus);
    let manifest = read_manifest(&corpus_path)?;
    // The read for identifiers takes a bound of no shape as one thread, so the bound is
    // refused before it.
    let threads = thread_count()?;
    // A corpus that holds a restricted identifier is refused before anything else is read.
    let refused = restricted_identifier_errors(&corpus_files(&corpus_path)?, Some(&corpus_path))?;
    for error in &refused {
        complain(error);
    }
    if !refused.is_empty() {
        return Ok(1);
    }
    let mut errors: Vec<String> = Vec::new();
    if let Some(matrix) = &args.matrix {
        let matrix = path_str(matrix);
        let items = matrix_items(&matrix)?;
        for (owner, names) in covers_lists(&manifest) {
            errors.extend(names_outside_matrix(&names, &items, &matrix, &owner));
        }
    }
    let datasets: Vec<Dataset> = args
        .dataset
        .iter()
        .map(|dataset| parse_dataset_arg(dataset))
        .collect::<Result<_, _>>()?;
    let family_specs = parse_family_columns(args.family_columns.as_deref())?;
    let required = genotype_class_names(&args.require_genotype_class)?;
    let keep_samples = manifest
        .get("keep_samples")
        .filter(|keep| keep.truthy())
        .and_then(Value::as_i64)
        .and_then(|keep| usize::try_from(keep).ok());
    if !required.is_empty() && keep_samples.is_none() {
        return Err(failure(
            "--require-genotype-class wants a corpus built with --keep-samples",
        ));
    }
    let records = corpus_records(&manifest)?;
    let by_name = record_name_index(records.iter().map(|record| record.vcf))?;
    let corpus = Corpus { records, by_name };
    let with_strata = !family_specs.is_empty();

    let mut observed: HashSet<String> = HashSet::new();
    let mut strata_observed: HashSet<String> = HashSet::new();
    let mut strata_covered: HashSet<String> = HashSet::new();
    for (sid, reference, input) in &datasets {
        // A family stratum ends in the row's variant class, which is read with the input.
        let substitutions = match with_strata.then(|| substitution_alts(input)).transpose() {
            Ok(substitutions) => substitutions,
            // The refusal that ends the run comes after the errors found before it.
            Err(refusal) => return Err(after_errors(&errors, refusal)),
        };
        let reader = RowReader::new(sid, &family_specs, substitutions.as_ref());
        let read = |block: &str| block_sets(block, &reader, sid, with_strata.then_some(&corpus));
        let take = |_: &str, found: BlockSets| {
            if let Some(place) = found.without_dataset {
                return Err(no_dataset(place));
            }
            observed.extend(found.sets);
            strata_observed.extend(found.strata);
            strata_covered.extend(found.covered);
            Ok(())
        };
        read_blocks(reference, threads, read, take)?;
    }
    let mut covered: HashSet<String> = HashSet::new();
    let manifest_records = manifest.get("records").and_then(Value::as_list);
    for (place, record) in manifest_records.unwrap_or_default().iter().enumerate() {
        covered.extend(combinations(record, place)?);
    }
    let mut missing: Vec<&String> = observed.difference(&covered).collect();
    missing.sort();
    say(&format!(
        "observed {}, covered by corpus {}, missing {}",
        observed.len(),
        observed.intersection(&covered).count(),
        missing.len()
    ));
    for name in missing.iter().take(50) {
        say(&format!("  MISSING {name}"));
    }
    let mut missing_strata: Vec<&String> = strata_observed.difference(&strata_covered).collect();
    missing_strata.sort();
    if with_strata {
        say(&format!(
            "family strata observed {}, covered by corpus {}, missing {}",
            strata_observed.len(),
            strata_covered.len(),
            missing_strata.len()
        ));
        for name in missing_strata.iter().take(50) {
            say(&format!("  MISSING family stratum {name}"));
        }
    }
    let mut missing_classes: Vec<String> = Vec::new();
    if let Some(keep) = keep_samples {
        let samples = names_under(&manifest, "samples");
        let mut wanted = names_under(&manifest, "required_genotype_classes");
        wanted.extend(required);
        let columns = sample_column_errors(&corpus_path, &samples, keep);
        let coverage = columns.and_then(|found| {
            errors.extend(found);
            genotype_coverage(&datasets, &corpus.records, &samples, keep, &wanted)
        });
        let (observed_classes, covered_classes, lacking) = match coverage {
            Ok(coverage) => coverage,
            // The refusal that ends the run comes after the errors found before it.
            Err(refusal) => return Err(after_errors(&errors, refusal)),
        };
        say(&format!(
            "genotype classes observed {observed_classes}, covered by corpus {covered_classes}, missing {}",
            lacking.len()
        ));
        for class in &lacking {
            say(&format!("  MISSING genotype class {class}"));
        }
        missing_classes = lacking;
    }
    let mut missing_record_classes: Vec<&str> = Vec::new();
    let class_key = manifest.get("record_class_key").and_then(Value::as_str);
    if let Some(key) = class_key.filter(|key| !key.is_empty()) {
        if matches!(
            manifest.get("seed"),
            Some(Value::List(_) | Value::Object(_))
        ) {
            return Err(Failure::own(
                "the manifest's seed is a list or an object, which seeds no generator",
            ));
        }
        // The samples of the pass are not read here, so no draw is made for them.
        let pass = record_class_pass(&datasets, key, &mut |_, _, _, _| {})?;
        let mut covered_classes: BTreeSet<&str> = BTreeSet::new();
        for ((sid, name, location), classes) in &pass.classes_of {
            if corpus
                .holds(sid, name, chromosome(location))
                .map_err(no_dataset)?
            {
                covered_classes.extend(classes);
            }
        }
        let seen: BTreeSet<&str> = pass.seen.keys().copied().collect();
        missing_record_classes = seen.difference(&covered_classes).copied().collect();
        say(&format!(
            "record classes observed {}, covered by corpus {}, missing {}",
            seen.len(),
            seen.intersection(&covered_classes).count(),
            missing_record_classes.len()
        ));
        for class in &missing_record_classes {
            say(&format!("  MISSING record class {class}"));
        }
    }
    for error in &errors {
        complain(error);
    }
    let complete = missing.is_empty()
        && missing_strata.is_empty()
        && missing_classes.is_empty()
        && missing_record_classes.is_empty()
        && errors.is_empty();
    Ok(if complete { 0 } else { 1 })
}

#[cfg(test)]
#[path = "check_tests.rs"]
mod tests;
