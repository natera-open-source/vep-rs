// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-consequence-class modes: one (dataset, engine)'s per-term counts from a run's saved
//! reports, and the pooling of a per-class CSV across datasets.

use crate::classes::{
    add_key_terms, count_terms_in_discordance, count_terms_in_sorted_keys, derive_per_class_counts,
    engine_order, key_fields_error, per_class_floor, pool_per_class_counts,
    read_per_class_by_dataset_csv, read_per_term_exclusions, read_vep_totals_cache,
    vep_totals_fingerprint, write_per_class_by_dataset_csv, write_per_class_pooled_csv,
    write_vep_totals_cache, Counter, DatasetLabels, Exclusions,
};
use crate::compare::{ended_without_contig_filter, resolve_canonical_contigs, sort_keys};
use crate::partition::{self, PartitionFiles, ScanSpec, Tuning};
use crate::{file_name_of, glob_sorted, py_resolve, Args, TempDir};
use std::collections::HashSet;
use std::io;
use std::path::Path;

/// The per-term tuples of one reference output file, its keys sorted whole on one thread.
fn file_term_totals(
    perl_file: &Path,
    canonical: Option<&HashSet<String>>,
    scratch: &Path,
) -> Result<Counter, String> {
    let name = file_name_of(perl_file);
    let sorted_keys = scratch.join(format!("{name}.perl.keys.sorted.txt"));
    sort_keys(perl_file, false, canonical, scratch, &sorted_keys)
        .map_err(|e| format!("{name}: {e}"))?;
    let totals = count_terms_in_sorted_keys(&sorted_keys).map_err(|e| e.to_string())?;
    let _ = std::fs::remove_file(&sorted_keys);
    Ok(totals)
}

/// The per-term tuples of one reference output file across the worker threads: its keys
/// routed into partitions, each sorted and counted on its own, the counts joined in
/// partition order.
fn file_term_totals_partitioned(
    perl_file: &Path,
    canonical: Option<&HashSet<String>>,
    scratch: &Path,
    tuning: &Tuning,
) -> io::Result<Counter> {
    let name = file_name_of(perl_file);
    let splitters = partition::splitters_for(&[perl_file], canonical, tuning)?;
    let keys = PartitionFiles::new(scratch, &format!("{name}.perl"), splitters.partitions());
    let spec = ScanSpec {
        canonical_contigs: canonical,
        bare_breakpoints: false,
        gvcf_sites: false,
    };
    partition::scan(perl_file, &spec, &splitters, &keys, tuning)?;
    let counted =
        crate::threads::run_indexed(tuning.threads, splitters.partitions(), |partition| {
            let sorted = keys.take_sorted(partition)?;
            let mut counter = Counter::default();
            for index in 0..sorted.len() {
                add_key_terms(&mut counter, sorted.get(index))
                    .map_err(|fields| key_fields_error(&keys.path(partition), fields))?;
            }
            Ok(counter)
        })?;
    let mut totals = Counter::default();
    for counter in &counted {
        totals.update(counter);
    }
    Ok(totals)
}

/// Counts one (dataset, engine)'s per-term concordance from a run's saved reports.
pub fn run_by_consequence_class(args: &Args) -> Result<u8, String> {
    let perl_dir = py_resolve(args.perl_dir.as_deref().expect("required"));
    let discordant_in = py_resolve(args.discordant_in.as_deref().expect("required"));
    let class_csv = py_resolve(args.class_csv.as_deref().expect("required"));
    if !discordant_in.is_file() {
        return Err(format!(
            "discordance report not found: {}",
            discordant_in.display()
        ));
    }
    let canonical = match resolve_canonical_contigs(args) {
        Ok(canonical) => canonical,
        Err(message) => return ended_without_contig_filter(args, message),
    };
    let perl_files = glob_sorted(&perl_dir, &args.glob)?;
    if perl_files.is_empty() {
        return Err(format!(
            "No files found in {} matching {}",
            perl_dir.display(),
            args.glob
        ));
    }
    let cache = args.vep_totals_cache.as_deref().map(py_resolve);
    let contig_list: Option<Vec<String>> =
        canonical.as_ref().map(|set| set.iter().cloned().collect());
    let fingerprint = vep_totals_fingerprint(&perl_files, contig_list.as_deref())
        .map_err(|e| format!("{}: {e}", perl_dir.display()))?;
    let cached = cache
        .as_deref()
        .and_then(|c| read_vep_totals_cache(c, &fingerprint));
    let vep_totals = match cached {
        Some(totals) => {
            println!(
                "VEP per-term totals: cache hit {}",
                cache.as_deref().expect("cache path").display()
            );
            totals
        }
        None => {
            let mut totals = Counter::default();
            let tmp = TempDir::create(&std::env::temp_dir(), "per_class_tmp_")
                .map_err(|e| format!("temporary directory: {e}"))?;
            let tuning = Tuning::new(args.threads);
            for perl_file in &perl_files {
                let file_totals = if args.threads > 1 {
                    file_term_totals_partitioned(perl_file, canonical.as_ref(), &tmp.path, &tuning)
                        .map_err(|e| format!("{}: {e}", file_name_of(perl_file)))?
                } else {
                    file_term_totals(perl_file, canonical.as_ref(), &tmp.path)?
                };
                totals.update(&file_totals);
            }
            if let Some(cache) = &cache {
                write_vep_totals_cache(cache, &fingerprint, &totals)
                    .map_err(|e| format!("{}: {e}", cache.display()))?;
                println!("VEP per-term totals: cache miss, wrote {}", cache.display());
            }
            totals
        }
    };
    let (vep_only, engine_only) =
        count_terms_in_discordance(&discordant_in).map_err(|e| e.to_string())?;
    let counts = derive_per_class_counts(&vep_totals, &vep_only, &engine_only)?;
    let assembly = if args.assembly_label.is_empty() {
        args.assembly.clone().unwrap_or_default()
    } else {
        args.assembly_label.clone()
    };
    let dataset = args.dataset.as_deref().expect("required");
    let labels = DatasetLabels {
        dataset,
        assembly: &assembly,
        engine: &args.engine,
        sweep_id: args.sweep_id.as_deref().expect("required"),
        binary_md5: args.binary_md5.as_deref().expect("required"),
    };
    let written = write_per_class_by_dataset_csv(&class_csv, &counts, &labels, args.append)
        .map_err(|e| e.to_string())?;
    let floor = per_class_floor(&counts, args.min_tuples);
    println!(
        "wrote per-class counts: {} ({written} rows {}/{dataset})",
        class_csv.display(),
        args.engine
    );
    println!(
        "terms: {}  VEP tuples: {}  VEP-only: {}  engine-only: {}",
        counts.len(),
        vep_totals.total(),
        vep_only.total(),
        engine_only.total()
    );
    match (floor.f1, floor.term) {
        (Some(f1), Some(term)) => println!(
            "this dataset's floor: {f1} on {term} ({} terms at or above {} tuples, {} below). \
             NOT the reported floor, which is pooled across datasets -- see \
             --pool-consequence-class.",
            floor.n_eligible, args.min_tuples, floor.n_below_threshold
        ),
        _ => println!(
            "no term reaches {} VEP tuples on this dataset, so it contributes no per-dataset floor",
            args.min_tuples
        ),
    }
    Ok(0)
}

/// Pools a per-(dataset, term) CSV into the per-term artefact.
pub fn run_pool_consequence_class(args: &Args) -> Result<u8, String> {
    let class_csv = py_resolve(args.class_csv.as_deref().expect("required"));
    let pooled_csv = py_resolve(args.pooled_csv.as_deref().expect("required"));
    if !class_csv.is_file() {
        return Err(format!("per-class CSV not found: {}", class_csv.display()));
    }
    let (by_engine, pins) = read_per_class_by_dataset_csv(&class_csv).map_err(|e| e.to_string())?;
    if by_engine.is_empty() {
        return Err(format!("{} carries no rows to pool", class_csv.display()));
    }
    let pooled: Vec<(String, Vec<_>)> = by_engine
        .iter()
        .map(|(engine, rows)| pool_per_class_counts(rows).map(|p| (engine.clone(), p)))
        .collect::<Result<_, _>>()?;
    let mut adjusted: Option<Exclusions> = None;
    if let Some(per_term) = &args.per_term {
        if !pooled.iter().any(|(e, _)| *e == args.per_term_engine) {
            return Err(format!(
                "--per-term-engine '{}' has no rows in {}",
                args.per_term_engine,
                class_csv.display()
            ));
        }
        adjusted =
            Some(read_per_term_exclusions(&py_resolve(per_term)).map_err(|e| e.to_string())?);
    }
    let written = write_per_class_pooled_csv(
        &pooled_csv,
        &pooled,
        &pins,
        adjusted
            .as_ref()
            .map(|ex| (args.per_term_engine.as_str(), ex)),
    )
    .map_err(|e| e.to_string())?;
    let tail = match adjusted {
        Some(_) => format!("; adj_f1 filled for {})", args.per_term_engine),
        None => ")".to_string(),
    };
    println!(
        "wrote pooled per-class F1: {} ({written} rows{tail}",
        pooled_csv.display()
    );
    for engine in engine_order(pooled.iter().map(|(e, _)| e.as_str())) {
        let rows = &pooled
            .iter()
            .find(|(e, _)| e == engine)
            .expect("pooled engine")
            .1;
        let floor = per_class_floor(rows, args.min_tuples);
        match (floor.f1, floor.term) {
            (Some(f1), Some(term)) => println!(
                "{engine}: floor {f1} on {term}; {} of {} terms at or above {} pooled VEP \
                 tuples, {} below",
                floor.n_eligible,
                rows.len(),
                args.min_tuples,
                floor.n_below_threshold
            ),
            _ => println!(
                "{engine}: no term clears {} pooled VEP tuples, so no floor can be reported \
                 under the fixed threshold",
                args.min_tuples
            ),
        }
    }
    Ok(0)
}
