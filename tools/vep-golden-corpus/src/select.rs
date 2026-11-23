// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `select` subcommand: exemplar records for every consequence set the reference
//! outputs hold, written as `variants.vcf` and the manifest.

use crate::cache_cost::CacheCost;
use crate::classes::{
    genotype_class_names, parse_family_columns, record_class_key, record_genotype_classes,
};
use crate::cli::SelectArgs;
use crate::constants::{CLASSIFY_KEYS, DEFAULT_SPECIES, FIXED_COLUMNS, GENERATOR_VERSION};
use crate::error::{complain, failure, say, Failure};
use crate::identifiers::{fabricated_identifiers, holds_identifier_bytes};
use crate::json::{Object, Value};
use crate::manifest::{
    covers_names, covers_owner, flag_set_name, is_species_name, manifest_path, matrix_items,
    names_outside_matrix, names_under, read_manifest, reference_statements_to_discard,
    report_discarded_statements, run_options, without_repeats, write_manifest,
};
use crate::names::{
    is_symbolic_or_breakend, multi_allelic_stratum, names_of_record, record_names, record_pos,
    vep_auto_names,
};
use crate::pytext::{is_file, os_error_text, path_join, path_str, py_repr, py_repr_str};
use crate::rows::blocks::thread_count;
use crate::rows::{
    kept_samples, parse_dataset_arg, record_class_pass, vcf_lines, ByRecordName, Dataset,
    RecordClassPass, RecordRows, Row, EARLIER_INPUT,
};
use crate::sampling::{
    choose_exemplars, ensure_class_exemplars, ensure_family_exemplars, exemplar_order,
    reservoir_pass, reservoir_step, reservoir_step_with, OrderedMap, Random,
};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::Write;

/// What one scan of the input VCFs yields (`resolve_records`).
#[derive(Debug, Default)]
pub struct ResolvedInputs {
    /// Per dataset id, the record lines whose Uploaded_variation (the ID, or a name Ensembl
    /// VEP writes) was wanted, by chromosome and that name.
    pub found: HashMap<String, ByRecordName<String>>,
    /// The `#CHROM` header line the corpus writes: eight columns, or `FORMAT` and the kept
    /// samples after them.
    pub header: String,
    /// Per dataset id and multi-allelic stratum, a uniform sample of records.
    pub multi: BTreeMap<String, BTreeMap<&'static str, Vec<String>>>,
    /// With kept samples, their names.
    pub samples: Vec<String>,
    /// With kept samples, how many records show each genotype class.
    pub genotype_seen: BTreeMap<&'static str, usize>,
    /// With kept samples, a uniform sample of records per genotype class, each with its
    /// dataset id.
    pub genotype_pool: BTreeMap<&'static str, Vec<(String, String)>>,
}

/// Scans each input VCF once (`resolve_records`). A record line is its first eight columns,
/// or with `keep_samples` its first eight, `FORMAT` and the first `keep_samples` sample
/// columns, passed through `fabricated_identifiers` before it is kept. A record that step
/// would rename (its ID is such an identifier) is never sampled, and is refused when a
/// reference row names it, since the reference's rows carry the name as read. The
/// multi-allelic sample (`multi_per_dataset` sequence records per dataset and stratum) and the
/// genotype-class sample (`genotype_per_class` records per class the kept samples show) are
/// records no consequence set can select for; their draws are those of `random.Random(seed)`
/// over the records in file order.
pub fn resolve_records(
    inputs: &[(String, String)],
    wanted: &HashMap<String, HashSet<String>>,
    multi_per_dataset: i64,
    seed: i64,
    keep_samples: usize,
    genotype_per_class: usize,
) -> Result<ResolvedInputs, Failure> {
    let mut rng = Random::new(seed);
    let mut resolved = ResolvedInputs::default();
    let mut multi_seen: HashMap<(String, &'static str), usize> = HashMap::new();
    let mut header: Option<String> = None;
    let kept_columns = if keep_samples > 0 {
        9 + keep_samples
    } else {
        8
    };
    let per_class = i64::try_from(genotype_per_class).unwrap_or(i64::MAX);
    for (sid, path) in inputs {
        let want = wanted.get(sid).filter(|names| !names.is_empty());
        if want.is_none() && multi_per_dataset == 0 && keep_samples == 0 {
            continue;
        }
        let mut named = false;
        let mut renamed = 0usize;
        let mut lines = vcf_lines(path, Some(kept_columns));
        while let Some(line) = lines.next_columns() {
            let (is_header, f) = line?;
            let kept = &f[..kept_columns.min(f.len())];
            if is_header {
                if keep_samples > 0 {
                    resolved.samples =
                        kept_samples(&f, keep_samples, path, &resolved.samples, EARLIER_INPUT)?;
                    named = true;
                }
                if header.is_none() {
                    header = Some(kept.join("\t"));
                }
                continue;
            }
            if keep_samples > 0 && !named {
                return Err(failure(format!(
                    "{path} has no #CHROM header line naming its samples"
                )));
            }
            let (chrom, ref_allele) = (f[0], f[3]);
            let pos = record_pos(f[0], f[1], f[2])?;
            let alts: Vec<&str> = f[4].split(',').collect();
            // An identifier lies inside one column, so a line whose columns hold none is kept
            // as read and is joined only where a sample takes it.
            let holds = |column: &&str| holds_identifier_bytes(column.as_bytes());
            let fabricated = if kept.iter().any(holds) {
                let as_read = kept.join("\t");
                let record = fabricated_identifiers(&as_read, &BTreeSet::new()).into_owned();
                let names_as_read = record_names(&as_read)?;
                if record_names(&record)? != names_as_read {
                    let auto = vep_auto_names(chrom, pos, ref_allele, &alts);
                    let mut names = names_as_read.iter().chain(&auto);
                    let chosen = want.is_some_and(|want| names.any(|name| want.contains(name)));
                    renamed += usize::from(chosen);
                    continue;
                }
                Some(record)
            } else {
                None
            };
            let record = || fabricated.clone().unwrap_or_else(|| kept.join("\t"));
            if keep_samples > 0 {
                // In sorted order, so that the draws fall on the classes in one order.
                for class in record_genotype_classes(&f, keep_samples)? {
                    let seen = resolved.genotype_seen.entry(class).or_default();
                    *seen += 1;
                    if genotype_per_class > 0 {
                        let pool = resolved.genotype_pool.entry(class).or_default();
                        let item = || (sid.clone(), record());
                        reservoir_step_with(&mut rng, pool, *seen, per_class, item);
                    }
                }
            }
            let sequence_alleles = !alts
                .iter()
                .any(|alt| is_symbolic_or_breakend(alt) || *alt == "*");
            if multi_per_dataset != 0 && alts.len() > 1 && sequence_alleles {
                let stratum = multi_allelic_stratum(ref_allele, &alts)?;
                let seen = multi_seen.entry((sid.clone(), stratum)).or_default();
                *seen += 1;
                let strata = resolved.multi.entry(sid.clone()).or_default();
                let sample = strata.entry(stratum).or_default();
                reservoir_step_with(&mut rng, sample, *seen, multi_per_dataset, record);
            }
            let Some(want) = want else {
                continue;
            };
            let mut names = names_of_record(f[0], f[1], f[2], f[3], f[4])?;
            names.extend(vep_auto_names(chrom, pos, ref_allele, &alts));
            if let Some(name) = names.iter().find(|name| want.contains(*name)) {
                let found = resolved.found.entry(sid.clone()).or_default();
                found.extend(chrom, name, [record()]);
            }
        }
        if renamed > 0 {
            let plural = if renamed == 1 { "" } else { "s" };
            return Err(failure(format!(
                "{path} names {renamed} record{plural} the references chose by an identifier \
                 under restricted redistribution terms in the ID column; a corpus carries no \
                 such record, since its reference rows are named by the identifier"
            )));
        }
    }
    if keep_samples > 0 && resolved.samples.is_empty() {
        return Err(failure(
            "no input has a #CHROM header line naming its samples",
        ));
    }
    resolved.header = header.unwrap_or_else(|| FIXED_COLUMNS.join("\t"));
    Ok(resolved)
}

/// A record the corpus carries outside the sampling: its dataset, the label an error names
/// it by (the `--include` argument, a kept record's first five columns, or a carried
/// `included_records` name) and the Uploaded_variation names its reference rows can carry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NamedRecord {
    pub dataset: String,
    pub label: String,
    pub names: Vec<String>,
}

/// A record of the manifest being regenerated, as `--keep-records` reads it: an object
/// with a `vcf` line and a `dataset`.
fn kept_record(record: &Value, place: usize) -> Result<(&Value, &str), Failure> {
    let read = record
        .get("dataset")
        .zip(record.get("vcf").and_then(Value::as_str));
    read.ok_or_else(|| {
        Failure::own(format!(
            "record {place} of the manifest being regenerated has no dataset and vcf line"
        ))
    })
}

/// The records `--include` names and, under `--keep-records`, every record of the manifest
/// being regenerated and every name of its `included_records`, so that a carried name is
/// held to the checks an `--include` gets; a record of a dataset that is not among
/// `--dataset` is refused.
pub fn named_records(
    args: &SelectArgs,
    existing: &Object,
    dataset_ids: &HashSet<&str>,
) -> Result<Vec<NamedRecord>, Failure> {
    let outside = |label: &str, dataset: &str| {
        failure(format!(
            "{label} names dataset {dataset}, which is not among --dataset"
        ))
    };
    let mut named: Vec<NamedRecord> = Vec::new();
    // A record whose dataset is no text is among no `--dataset`: the first of them in the
    // order the records are named, with the place it holds there.
    let mut no_text: Option<(usize, Failure)> = None;
    for arg in &args.include {
        let Some((sid, name)) = arg
            .split_once('=')
            .filter(|(sid, name)| !sid.is_empty() && !name.is_empty())
        else {
            return Err(failure(format!(
                "--include wants sid=name, got {}",
                py_repr_str(arg)
            )));
        };
        named.push(NamedRecord {
            dataset: sid.to_string(),
            label: format!("--include {arg}"),
            names: vec![name.to_string()],
        });
    }
    if args.keep_records {
        if existing.is_empty() {
            return Err(failure(format!(
                "--keep-records wants a manifest to regenerate in {}",
                args.out
            )));
        }
        let records = match existing.get("records").filter(|records| records.truthy()) {
            None => &[][..],
            Some(Value::List(records)) => records.as_slice(),
            Some(_) => {
                return Err(Failure::own(
                    "the records of the manifest being regenerated are not a list",
                ))
            }
        };
        for (place, record) in records.iter().enumerate() {
            let (dataset, vcf) = kept_record(record, place)?;
            let columns: Vec<&str> = vcf.split('\t').take(5).collect();
            let label = format!("kept record {}", columns.join(" "));
            let names = record_names(vcf)?;
            match dataset {
                Value::Str(dataset) => named.push(NamedRecord {
                    dataset: dataset.clone(),
                    label,
                    names,
                }),
                Value::List(_) | Value::Object(_) => {
                    return Err(Failure::own(format!(
                    "record {place} of the manifest being regenerated has a dataset that is no name"
                )))
                }
                other if no_text.is_none() => {
                    no_text = Some((named.len(), outside(&label, &py_repr(other))));
                }
                _ => {}
            }
        }
        let included = existing.get("included_records").and_then(Value::as_object);
        for (sid, names) in included.into_iter().flat_map(Object::iter) {
            for name in names.as_list().unwrap_or_default() {
                let name = name.as_str().unwrap_or_default();
                named.push(NamedRecord {
                    dataset: sid.to_string(),
                    label: format!("carried include {sid}={name}"),
                    names: vec![name.to_string()],
                });
            }
        }
    }
    for (place, record) in named.iter().enumerate() {
        if let Some((_, refusal)) = no_text.take_if(|(before, _)| *before == place) {
            return Err(refusal);
        }
        if !dataset_ids.contains(record.dataset.as_str()) {
            return Err(outside(&record.label, &py_repr_str(&record.dataset)));
        }
    }
    match no_text {
        Some((_, refusal)) => Err(refusal),
        None => Ok(named),
    }
}

/// The `--include` arguments as dataset id to sorted names, merged with the
/// `included_records` a regeneration under `--keep-records` carries from the manifest it
/// keeps: the manifest's `included_records` (`included_by_dataset`).
pub fn included_by_dataset(includes: &[String], carried: Option<&Object>) -> Object {
    let mut by_dataset: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (sid, names) in carried.into_iter().flat_map(Object::iter) {
        let names = names.as_list().unwrap_or_default();
        let names = names.iter().filter_map(Value::as_str).map(str::to_string);
        by_dataset.entry(sid.to_string()).or_default().extend(names);
    }
    for arg in includes {
        let (sid, name) = arg.split_once('=').unwrap_or((arg, ""));
        by_dataset
            .entry(sid.to_string())
            .or_default()
            .insert(name.to_string());
    }
    by_dataset
        .into_iter()
        .map(|(sid, names)| (sid, Value::strings(&names.into_iter().collect::<Vec<_>>())))
        .collect()
}

/// One record of the corpus as `select` gathers it.
#[derive(Debug, Default)]
struct Record {
    dataset: String,
    vcf: String,
    combinations: Vec<String>,
    /// Per exemplar row: its Uploaded_variation, and with family specs its family stratum.
    exemplar_rows: Vec<(String, Option<String>)>,
    strata: Vec<&'static str>,
    genotype_classes: Vec<&'static str>,
    record_classes: Vec<&'static str>,
}

/// The dataset id and Uploaded_variation and Location of a row's record.
fn record_of(row: &Row) -> (String, String, String) {
    (
        row.dataset.clone(),
        row.uploaded_variation.clone(),
        row.location.clone(),
    )
}

/// The chromosome of a Location: the text before its first `:`.
fn chromosome(location: &str) -> &str {
    location.split_once(':').map_or(location, |(head, _)| head)
}

/// `Path(path).name` of a path as `str(Path(...))` writes it.
fn file_name(path: &str) -> &str {
    if path == "." {
        return "";
    }
    path.rsplit_once('/').map_or(path, |(_, name)| name)
}

/// The place of a chromosome among the records of a corpus: a name of ASCII digits by its
/// number (the count of its digits, then the digits), then X, Y and MT, then every other.
type ChromosomeRank = (u8, usize, String, u8);

/// The key the records of a corpus are ordered by: chromosomes named by a number in
/// number order, then X, Y and MT, then every other; then POS, REF and ALT.
fn record_order(vcf: &str) -> Result<(ChromosomeRank, i64, String, String), Failure> {
    let f: Vec<&str> = vcf.split('\t').collect();
    let chrom = f[0];
    let numbered = !chrom.is_empty() && chrom.bytes().all(|byte| byte.is_ascii_digit());
    let rank = if numbered {
        let digits = chrom.trim_start_matches('0');
        (0, digits.len(), digits.to_string(), 0)
    } else {
        let named = match chrom {
            "X" => 1,
            "Y" => 2,
            "MT" => 3,
            _ => 9,
        };
        (1, 0, String::new(), named)
    };
    let pos = record_pos(f[0], f[1], f[2])?;
    Ok((rank, pos, f[3].to_string(), f[4].to_string()))
}

/// The counts of `names` by name, in the names' sorted order: `dict(sorted(Counter(...).items()))`.
fn sorted_counts<S: AsRef<str>>(names: impl IntoIterator<Item = S>) -> Object {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for name in names {
        *counts.entry(name.as_ref().to_string()).or_default() += 1;
    }
    counts
        .into_iter()
        .map(|(name, count)| (name, Value::from(count)))
        .collect()
}

/// A failure of the tool's own for a file of the corpus that cannot be written.
fn not_writable(path: &str, error: &std::io::Error) -> Failure {
    Failure::own(format!(
        "{path} is not writable ({})",
        os_error_text(error, path)
    ))
}

/// Runs `select`.
pub fn run(args: &SelectArgs) -> Result<i32, Failure> {
    let datasets: Vec<Dataset> = args
        .dataset
        .iter()
        .map(|dataset| parse_dataset_arg(dataset))
        .collect::<Result<_, _>>()?;
    let out = path_str(&args.out);
    let existing = if is_file(&manifest_path(&out)) {
        read_manifest(&out)?
    } else {
        Object::new()
    };
    let carried_species = existing.get("species").and_then(Value::as_str);
    let species = match &args.species {
        Some(species) => species.as_str(),
        None => carried_species
            .filter(|species| !species.is_empty())
            .unwrap_or(DEFAULT_SPECIES),
    };
    if !is_species_name(species) {
        complain(&format!(
            "--species {} is not an Ensembl species name (lower-case letters, digits and underscores)",
            py_repr_str(species)
        ));
        return Ok(1);
    }
    let dataset_ids: HashSet<&str> = datasets.iter().map(|(sid, _, _)| sid.as_str()).collect();
    let named = named_records(args, &existing, &dataset_ids)?;
    let mut covers = names_under(&existing, "covers");
    covers.extend(covers_names(&args.covers)?);
    let covers = without_repeats(covers);
    let carried = |key: &str| existing.get(key).filter(|value| value.truthy()).cloned();
    let existing_sets = existing
        .get("flag_sets")
        .and_then(Value::as_list)
        .unwrap_or_default();
    let set_runs = existing_sets.iter().map(|set| {
        let run = set.as_object().expect("a flag set of the shape held");
        (format!("flag set {}", flag_set_name(set)), run)
    });
    let discarded = reference_statements_to_discard(
        std::iter::once(("the manifest's".to_string(), &existing)).chain(set_runs),
        "select",
        args.drop_reference_statements,
    )?;
    // A set's divergence keys name the records of the manifest being regenerated by
    // ordinal, so they go with the manifest's own.
    let flag_sets: Vec<Value> = existing_sets
        .iter()
        .map(|set| {
            let mut set = set.clone();
            if let Some(entry) = set.as_object_mut() {
                entry.retain(|key, _| !CLASSIFY_KEYS.contains(&key));
            }
            set
        })
        .collect();
    if covers.is_empty() {
        complain(
            "select writes no manifest without --covers: name the parity-matrix rows the corpus exercises",
        );
        return Ok(1);
    }
    if let Some(matrix) = &args.matrix {
        let matrix = path_str(matrix);
        let items = matrix_items(&matrix)?;
        let errors = names_outside_matrix(&covers, &items, &matrix, &covers_owner(None));
        for error in &errors {
            complain(error);
        }
        if !errors.is_empty() {
            return Ok(1);
        }
    }
    let family_specs = parse_family_columns(args.family_columns.as_deref())?;
    let required = genotype_class_names(&args.require_genotype_class)?;
    if args.keep_samples < 0 {
        return Err(failure(format!(
            "--keep-samples wants a count of sample columns, got {}",
            args.keep_samples
        )));
    }
    if args.genotype_exemplars < 0 {
        return Err(failure(format!(
            "--genotype-exemplars wants a count of records per class, got {}",
            args.genotype_exemplars
        )));
    }
    let count = |number: i64, option: &str| {
        usize::try_from(number)
            .map_err(|_| Failure::own(format!("{option} {number} is past what the host counts")))
    };
    let keep_samples = count(args.keep_samples, "--keep-samples")?;
    let genotype_exemplars = count(args.genotype_exemplars, "--genotype-exemplars")?;
    if !required.is_empty() && keep_samples == 0 {
        return Err(failure(
            "--require-genotype-class wants --keep-samples: a corpus without sample columns shows no genotype",
        ));
    }
    let class_key = record_class_key(args.record_classes.as_deref())?;
    if args.record_class_exemplars < 0 {
        return Err(failure(format!(
            "--record-class-exemplars wants a count of records per class, got {}",
            args.record_class_exemplars
        )));
    }
    let record_class_exemplars = count(args.record_class_exemplars, "--record-class-exemplars")?;
    let cache = args.cache.as_deref().filter(|cache| !cache.is_empty());
    if let Some(cache) = cache {
        let transcripts = path_join(&path_str(cache), "transcripts");
        if !std::fs::metadata(&transcripts).is_ok_and(|found| found.is_dir()) {
            return Err(failure(format!("{cache} has no transcripts/ directory")));
        }
    }
    let threads = thread_count()?;
    let mut names_by_dataset: HashMap<String, HashSet<String>> = HashMap::new();
    for record in &named {
        let names = names_by_dataset.entry(record.dataset.clone()).or_default();
        names.extend(record.names.iter().cloned());
    }
    let mut rng = Random::new(args.seed);
    let pass = reservoir_pass(
        &datasets,
        args.seed,
        &family_specs,
        &names_by_dataset,
        threads,
    )?;
    let rowed = |record: &NamedRecord| {
        record.names.iter().any(|name| {
            pass.named_rows
                .contains_key(&(record.dataset.clone(), name.clone()))
        })
    };
    let unrowed: Vec<&NamedRecord> = named.iter().filter(|record| !rowed(record)).collect();
    for record in &unrowed {
        complain(&format!(
            "{} has no row in the {} reference",
            record.label, record.dataset
        ));
    }
    if !unrowed.is_empty() {
        return Ok(1);
    }
    let cost = match cache {
        Some(cache) => Some(CacheCost::new(&path_str(cache), args.flank)?),
        None => None,
    };
    let cost = cost.as_ref();
    if let Some(cost) = cost {
        let sampled = pass.reservoir.values().chain(pass.families.values());
        cost.read_ahead(sampled.flat_map(|sample| &sample.rows), threads);
    }
    let mut chosen = choose_exemplars(
        &pass.reservoir,
        args.k,
        &mut rng,
        args.max_span,
        &args.focus_term,
        args.focus_k,
        cost,
    )?;
    ensure_class_exemplars(&mut chosen, &pass.reservoir, args.max_span, cost)?;
    if !family_specs.is_empty() {
        ensure_family_exemplars(&mut chosen, &pass.families, args.max_span, cost)?;
    }
    for rows in pass.named_rows.values() {
        for row in rows {
            let bucket = chosen.entry(row.consequence_set.as_str());
            if !bucket.contains(row) {
                bucket.push(row.clone());
            }
        }
    }
    let mut class_pass = RecordClassPass::default();
    if let Some(key) = class_key {
        // Every record class the references show gets its exemplars: a class short of them
        // takes, from its own sample, the records that sort first in `exemplar_order`, each
        // counted towards every class it shows.
        let mut class_rng = Random::new(args.seed);
        class_pass = record_class_pass(&datasets, key, &mut |bucket, item, offered, capacity| {
            reservoir_step(&mut class_rng, bucket, item, offered, capacity)
        })?;
        if let Some(cost) = cost {
            let sampled = class_pass.samples.values().flatten();
            cost.read_ahead(sampled.flat_map(|rows| rows.iter()), threads);
        }
        let mut chosen_records: HashSet<(String, String, String)> =
            chosen.values().flatten().map(record_of).collect();
        let mut showing: HashMap<&str, usize> = HashMap::new();
        for record in &chosen_records {
            for class in class_pass.classes_of.get(record).into_iter().flatten() {
                *showing.entry(*class).or_default() += 1;
            }
        }
        for class in class_pass.seen.keys() {
            let empty = Vec::new();
            let sample = class_pass.samples.get(class).unwrap_or(&empty);
            let mut candidates: Vec<(_, &RecordRows)> = Vec::with_capacity(sample.len());
            for rows in sample {
                let mut first = None;
                for row in rows.iter() {
                    let order = exemplar_order(row, args.max_span, cost)?;
                    if first.is_none_or(|held| order < held) {
                        first = Some(order);
                    }
                }
                candidates.push((first, rows));
            }
            candidates.sort_by_key(|(order, _)| *order);
            for (_, rows) in candidates {
                if showing.get(class).copied().unwrap_or(0) >= record_class_exemplars {
                    break;
                }
                let record = record_of(&rows[0]);
                if chosen_records.contains(&record) {
                    continue;
                }
                for shown in class_pass.classes_of.get(&record).into_iter().flatten() {
                    *showing.entry(*shown).or_default() += 1;
                }
                chosen_records.insert(record);
                for row in rows.iter() {
                    chosen.entry(row.consequence_set.as_str()).push(row.clone());
                }
            }
        }
    }
    let mut wanted: HashMap<String, HashSet<String>> = HashMap::new();
    for row in chosen.values().flatten() {
        let names = wanted.entry(row.dataset.clone()).or_default();
        if !names.contains(&row.uploaded_variation) {
            names.insert(row.uploaded_variation.clone());
        }
    }
    let inputs: Vec<(String, String)> = datasets
        .iter()
        .map(|(sid, _, input)| (sid.clone(), input.clone()))
        .collect();
    let resolved = resolve_records(
        &inputs,
        &wanted,
        args.multi_allelic_per_dataset,
        args.seed,
        keep_samples,
        genotype_exemplars,
    )?;
    let mut records: OrderedMap<String, Record> = OrderedMap::new();
    let new_record = |dataset: &str, vcf: &str| Record {
        dataset: dataset.to_string(),
        vcf: vcf.to_string(),
        ..Record::default()
    };
    let mut missing: Vec<String> = Vec::new();
    for (combo, rows) in chosen.iter() {
        for row in rows {
            let lines = resolved
                .found
                .get(&row.dataset)
                .and_then(|found| found.get(chromosome(&row.location), &row.uploaded_variation))
                .unwrap_or_default();
            if lines.is_empty() {
                missing.push(format!(
                    "{}\t{}\t{}",
                    row.dataset, row.uploaded_variation, row.location
                ));
                continue;
            }
            for line in lines {
                let key = format!("{}\t{line}", row.dataset);
                let place =
                    records.place_or_insert_with(key.as_str(), || new_record(&row.dataset, line));
                let record = records.at_mut(place);
                if !record.combinations.contains(combo) {
                    record.combinations.push(combo.clone());
                }
                let stratum = (!family_specs.is_empty()).then(|| row.family_stratum());
                record
                    .exemplar_rows
                    .push((row.uploaded_variation.clone(), stratum));
            }
        }
    }
    for (sid, strata) in &resolved.multi {
        for (stratum, lines) in strata {
            for line in lines {
                let key = format!("{sid}\t{line}");
                let place = records.place_or_insert_with(key.as_str(), || new_record(sid, line));
                records.at_mut(place).strata.push(*stratum);
            }
        }
    }
    let genotype_classes_of = |vcf: &str| -> Result<Vec<&'static str>, Failure> {
        let columns: Vec<&str> = vcf.split('\t').collect();
        Ok(record_genotype_classes(&columns, keep_samples)?
            .into_iter()
            .collect())
    };
    if keep_samples > 0 {
        let mut showing: HashMap<&str, usize> = HashMap::new();
        for record in records.values_mut() {
            record.genotype_classes = genotype_classes_of(&record.vcf)?;
            for class in &record.genotype_classes {
                *showing.entry(*class).or_default() += 1;
            }
        }
        for class in resolved.genotype_seen.keys() {
            let pool = resolved.genotype_pool.get(class);
            for (sid, line) in pool.into_iter().flatten() {
                if showing.get(class).copied().unwrap_or(0) >= genotype_exemplars {
                    break;
                }
                let key = format!("{sid}\t{line}");
                if !records.contains_key(key.as_str()) {
                    let mut record = new_record(sid, line);
                    record.genotype_classes = genotype_classes_of(line)?;
                    for shown in &record.genotype_classes {
                        *showing.entry(*shown).or_default() += 1;
                    }
                    records.insert(key, record);
                }
            }
        }
    }
    if class_key.is_some() {
        // By the record's names and not its exemplar rows: a record the multi-allelic or
        // the genotype sampling carried has none, and shows its classes all the same.
        let mut classes_by_name: HashMap<(&str, &str, &str), BTreeSet<&'static str>> =
            HashMap::new();
        for ((sid, name, location), classes) in &class_pass.classes_of {
            let key = (sid.as_str(), name.as_str(), chromosome(location));
            classes_by_name.entry(key).or_default().extend(classes);
        }
        for record in records.values_mut() {
            let chrom = record.vcf.split('\t').next().unwrap_or_default();
            let mut shown: BTreeSet<&'static str> = BTreeSet::new();
            for name in record_names(&record.vcf)? {
                let key = (record.dataset.as_str(), name.as_str(), chrom);
                shown.extend(classes_by_name.get(&key).into_iter().flatten());
            }
            record.record_classes = shown.into_iter().collect();
        }
    }
    if i64::try_from(records.len()).is_ok_and(|held| held > args.max_records) {
        // The first record of every combination, family stratum, genotype class and record
        // class and every named record survive the cap, which those alone may exceed; the
        // rest fill in sorted record order, so that the cut is the same in every run.
        let mut keep: HashSet<usize> = HashSet::new();
        let mut grouped: HashSet<String> = HashSet::new();
        for (place, (_, record)) in records.iter().enumerate() {
            let combinations = record.combinations.iter();
            let strata = record.exemplar_rows.iter().filter_map(|(_, s)| s.as_ref());
            let groups = combinations
                .map(|combo| format!("combination {combo}"))
                .chain(strata.map(|stratum| format!("stratum {stratum}")))
                .chain(
                    record
                        .genotype_classes
                        .iter()
                        .map(|class| format!("genotype class {class}")),
                )
                .chain(
                    record
                        .record_classes
                        .iter()
                        .map(|class| format!("record class {class}")),
                );
            for group in groups {
                if grouped.insert(group) {
                    keep.insert(place);
                }
            }
            let named_row = names_by_dataset.get(&record.dataset).is_some_and(|names| {
                record
                    .exemplar_rows
                    .iter()
                    .any(|(name, _)| names.contains(name))
            });
            if named_row {
                keep.insert(place);
            }
        }
        let mut by_key: Vec<(&String, usize)> = records.keys().zip(0..).collect();
        by_key.sort();
        for (_, place) in by_key {
            if i64::try_from(keep.len()).is_ok_and(|kept| kept >= args.max_records) {
                break;
            }
            keep.insert(place);
        }
        records.retain_places(|place| keep.contains(&place));
    }
    let mut ordered: Vec<(_, &Record)> = Vec::with_capacity(records.len());
    for record in records.values() {
        ordered.push((record_order(&record.vcf)?, record));
    }
    ordered.sort_by(|(a, _), (b, _)| a.cmp(b));
    let ordered: Vec<&Record> = ordered.into_iter().map(|(_, record)| record).collect();
    std::fs::create_dir_all(&out).map_err(|error| not_writable(&out, &error))?;
    let variants = path_join(&out, "variants.vcf");
    let mut text = format!("##fileformat=VCFv4.2\n{}\n", resolved.header);
    for record in &ordered {
        text.push_str(&record.vcf);
        text.push('\n');
    }
    std::fs::File::create(&variants)
        .and_then(|mut file| file.write_all(text.as_bytes()))
        .map_err(|error| not_writable(&variants, &error))?;
    drop(text);

    let covered: HashSet<&str> = ordered
        .iter()
        .flat_map(|record| record.combinations.iter().map(String::as_str))
        .collect();
    let genotype_classes = sorted_counts(ordered.iter().flat_map(|r| &r.genotype_classes));
    let mut uncovered_genotype_classes: BTreeSet<&str> =
        resolved.genotype_seen.keys().copied().collect();
    uncovered_genotype_classes.extend(required.iter().map(String::as_str));
    uncovered_genotype_classes.retain(|class| !genotype_classes.contains_key(class));
    let record_classes = sorted_counts(ordered.iter().flat_map(|r| &r.record_classes));
    let uncovered_record_classes: Vec<&str> = class_pass
        .seen
        .keys()
        .copied()
        .filter(|class| !record_classes.contains_key(class))
        .collect();
    let included_carried = existing
        .get("included_records")
        .and_then(Value::as_object)
        .filter(|_| args.keep_records);
    let included = included_by_dataset(&args.include, included_carried);
    let mut combos: Vec<(&String, usize)> = pass
        .reservoir
        .iter()
        .map(|(combo, sample)| (combo, sample.seen))
        .collect();
    combos.sort();
    let uncovered_combinations: Vec<&str> = combos
        .iter()
        .map(|(combo, _)| combo.as_str())
        .filter(|combo| !covered.contains(combo))
        .collect();
    let chosen_rows = || chosen.values().flatten();

    let mut manifest = Object::new();
    manifest.insert("release", args.release.as_str());
    manifest.insert("species", species);
    manifest.insert("assembly", args.assembly.as_str());
    manifest.insert("generator_version", GENERATOR_VERSION);
    manifest.insert("seed", args.seed);
    manifest.insert("exemplars_per_combination", args.k);
    if !args.focus_term.is_empty() {
        let mut terms = args.focus_term.clone();
        terms.sort();
        manifest.insert("focus_terms", Value::strings(&terms));
        manifest.insert(
            "focus_exemplars_per_combination",
            args.focus_k.map_or(Value::Null, Value::from),
        );
    }
    if !included.is_empty() {
        manifest.insert("included_records", included);
    }
    if args.keep_records {
        let kept = existing.get("records").and_then(Value::as_list);
        manifest.insert("kept_records", kept.map_or(0, <[Value]>::len));
    }
    manifest.insert("max_span", args.max_span);
    if cache.is_some() {
        manifest.insert("cache_cost_flank", args.flank);
    }
    if keep_samples > 0 {
        manifest.insert("keep_samples", keep_samples);
        manifest.insert("samples", Value::strings(&resolved.samples));
        manifest.insert("genotype_exemplars_per_class", genotype_exemplars);
        manifest.insert("required_genotype_classes", Value::strings(&required));
    }
    if let Some(key) = class_key {
        manifest.insert("record_class_key", key);
        manifest.insert("record_class_exemplars", record_class_exemplars);
    }
    manifest.update(run_options(&args.flag, args.fasta_name.as_deref()));
    manifest.insert("covers", Value::strings(&covers));
    let carried_keys = [
        "masks",
        "perl_mode",
        "reference_warnings",
        "plugin_files",
        "generated",
    ];
    for key in carried_keys {
        if let Some(value) = carried(key) {
            manifest.insert(key, value);
        }
    }
    if !flag_sets.is_empty() {
        manifest.insert("flag_sets", flag_sets);
    }
    // File names only: the directory each dataset was staged in belongs to the host that
    // built the corpus, not to the corpus.
    let mut dataset_entries = Object::new();
    for (sid, reference, input) in &datasets {
        let count = pass.per_dataset.get(sid);
        let mut entry = Object::new();
        entry.insert("reference", file_name(reference));
        entry.insert("input", file_name(input));
        entry.insert("rows", count.map_or(0, |count| count.rows));
        entry.insert(
            "combinations",
            count.map_or(0, |count| count.combinations.len()),
        );
        dataset_entries.insert(sid.as_str(), entry);
    }
    manifest.insert("datasets", dataset_entries);
    let mut combinations = Object::new();
    for (combo, seen) in &combos {
        let exemplars = ordered
            .iter()
            .filter(|record| record.combinations.contains(combo))
            .count();
        let mut entry = Object::new();
        entry.insert("observed_rows", *seen);
        entry.insert("exemplar_records", exemplars);
        combinations.insert(combo.as_str(), entry);
    }
    manifest.insert("combinations", combinations);
    let variant_classes = sorted_counts(chosen_rows().map(Row::klass));
    manifest.insert("variant_classes", variant_classes.clone());
    if !family_specs.is_empty() {
        manifest.insert("family_columns", Value::strings(&family_specs));
        manifest.insert(
            "family_strata",
            sorted_counts(chosen_rows().map(Row::family_stratum)),
        );
    }
    if keep_samples > 0 {
        manifest.insert("genotype_classes", genotype_classes.clone());
    }
    if class_key.is_some() {
        manifest.insert("record_classes", record_classes.clone());
    }
    let record_entries: Vec<Value> = ordered
        .iter()
        .map(|record| {
            let mut combinations = record.combinations.clone();
            combinations.sort();
            let mut entry = Object::new();
            entry.insert("dataset", record.dataset.as_str());
            entry.insert("vcf", record.vcf.as_str());
            entry.insert("combinations", Value::strings(&combinations));
            if !record.strata.is_empty() {
                entry.insert("strata", Value::strings(&record.strata));
            }
            if keep_samples > 0 {
                entry.insert("genotype_classes", Value::strings(&record.genotype_classes));
            }
            if class_key.is_some() {
                entry.insert("record_classes", Value::strings(&record.record_classes));
            }
            Value::Object(entry)
        })
        .collect();
    manifest.insert("records", record_entries);
    manifest.insert("multi_allelic_per_dataset", args.multi_allelic_per_dataset);
    manifest.insert(
        "uncovered_combinations",
        Value::strings(&uncovered_combinations),
    );
    if keep_samples > 0 {
        let uncovered: Vec<&str> = uncovered_genotype_classes.iter().copied().collect();
        manifest.insert("uncovered_genotype_classes", Value::strings(&uncovered));
    }
    if class_key.is_some() {
        manifest.insert(
            "uncovered_record_classes",
            Value::strings(&uncovered_record_classes),
        );
    }
    manifest.insert("unresolved_exemplars", Value::strings(&missing));
    write_manifest(&out, &manifest)?;
    report_discarded_statements("select", &discarded);

    let class_counts: Vec<String> = variant_classes
        .iter()
        .map(|(class, count)| format!("{class} {}", count.as_i64().unwrap_or_default()))
        .collect();
    let mut summary = format!(
        "combinations observed {}, covered {}, records {}, unresolved exemplars {}, variant classes {} ({})",
        combos.len(),
        covered.len(),
        ordered.len(),
        missing.len(),
        variant_classes.len(),
        class_counts.join(", ")
    );
    if keep_samples > 0 {
        summary.push_str(&format!(
            ", genotype classes observed {}, covered {}",
            resolved.genotype_seen.len(),
            genotype_classes.len()
        ));
    }
    if class_key.is_some() {
        summary.push_str(&format!(
            ", record classes observed {}, covered {}",
            class_pass.seen.len(),
            record_classes.len()
        ));
    }
    // A record both included and kept is two named records and one key of the named rows.
    if !named.is_empty() {
        summary.push_str(&format!(", named records {}", pass.named_rows.len()));
    }
    say(&summary);
    let complete = missing.is_empty()
        && uncovered_combinations.is_empty()
        && uncovered_genotype_classes.is_empty()
        && uncovered_record_classes.is_empty();
    Ok(if complete { 0 } else { 1 })
}

#[cfg(test)]
#[path = "select_tests.rs"]
mod tests;
