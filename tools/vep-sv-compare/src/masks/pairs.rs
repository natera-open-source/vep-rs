// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The registry's pair rules on this comparator's pairs: the consequence-set rules, then the
//! row-aware rules on the rows of the two outputs that carry a pair's tuples.
//!
//! A pair is a `(Location, Allele, Feature, Feature_type)` that each side holds exactly one
//! one-sided tuple for, the two differing in their consequence sets.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use memchr::memchr_iter;

use crate::failure::Failure;
use crate::input::{for_each_record, push_canonical_location, KeyPrep};
use crate::masks::NamedCounts;
use crate::normalize::{is_fastvep, normalize_allele_of, push_normalized_location};
use crate::registry::excluding_rules;
use crate::rules::breakend::{strip_chr, vep_rs_pos};
use crate::rules::coding_end::{AtACodingEnd, CodingEndSources};
use crate::rules::gvcf::perl_split;
use crate::rules::row::{
    pair_row, record_names_its_rows, replaced_reference, start_codon_stands, PairRow,
    RecordReadings,
};
use crate::rules::{ExclusionRule, START_COOCCURRENCE_SWAP};
use crate::strings::{Sides, Strings, Tuple, TupleSet};
use crate::text::{open_records, os_failure};
use crate::threads::crew_indexed_prefix;

/// One pair rule as this mask asks it.
pub trait PairRule {
    /// The bucket the rule's pairs are counted under.
    fn bucket(&self) -> &str;

    /// Whether both members of a pair the rule takes are set aside.
    fn excludes(&self) -> bool;

    /// The rule on the two consequence sets and the allele of a pair; false for a rule that
    /// reads rows.
    fn matches(&self, perl_csq: &str, rust_csq: &str, allele: &str) -> bool;

    /// A row-aware rule on the rows of the two outputs that carry the pair's tuples, the
    /// reference's first: every row of one against every row of the other, a rule that asks
    /// the coding ends of the pair's transcript with `at`. False for any other rule and for
    /// a side without a row.
    fn matches_rows_at(&self, perl: &[PairRow], rust: &[PairRow], at: &AtACodingEnd<'_>) -> bool;
}

impl PairRule for ExclusionRule {
    fn bucket(&self) -> &str {
        self.bucket
    }

    fn excludes(&self) -> bool {
        self.excludes
    }

    fn matches(&self, perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
        ExclusionRule::matches(self, perl_csq, rust_csq, allele)
    }

    fn matches_rows_at(&self, perl: &[PairRow], rust: &[PairRow], at: &AtACodingEnd<'_>) -> bool {
        ExclusionRule::matches_rows_at(self, perl, rust, at)
    }
}

fn io_failure(path: &Path, error: io::Error) -> Failure {
    os_failure(&error, Some(path))
}

/// Gives each of `rows` that carries no reference bases (`PairRow::reference`) the bases its
/// allele stands in for, read by position from the records of the input VCF `path` (plain or
/// gzipped): a record is asked about a row whose Location starts on one of the record's own
/// positions on its chromosome, a `chr` prefix set aside (`replaced_reference`). Where two
/// or more records read a row differently, the record whose identifier is the row's own
/// first column (`Uploaded_variation`) decides; where no record carries that identifier, or
/// more than one does and they read the row differently, the row keeps none. The file is
/// not opened when no row asks.
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
    let mut read: Vec<RecordReadings> = vec![RecordReadings::default(); rows.len()];
    for_each_record(open_records(path)?, |record| {
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
        let (chrom, pos, ids, reference, alt) = (field(), field(), field(), field(), field());
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
                if let Some(found) = found {
                    read[*at].add(found, record_names_its_rows(ids, &row.name));
                }
            }
        }
        Ok(())
    })?;
    for (row, found) in rows.iter_mut().zip(read) {
        if let Some(bases) = found.decided() {
            row.reference = Some(bases);
        }
    }
    Ok(())
}

/// The rows of one output behind the wanted tuples.
#[derive(Default)]
struct OutputRows {
    /// Per tuple, in the order its first row stands in the file: its distinct rows of nine
    /// columns or more, in file order, each read to its Codons column, the twelfth, where it
    /// holds one.
    rows: Vec<(Tuple, Vec<PairRow>)>,
    /// The tuples one of whose rows has fewer than nine columns.
    short: HashSet<Tuple>,
}

/// Reads the output `path` as a stream and keeps the rows that carry a tuple of `wanted`,
/// given by the text of its Location; a path that does not exist holds no row. A row is
/// keyed as the comparison keys it: the Location and the Allele normalised (`fastvep` maps a
/// symbolic token), the Feature and the Feature_type as written, the consequence set
/// normalised.
fn read_output_rows(
    strings: &Strings,
    path: &Path,
    wanted: &HashMap<&[u8], &[Tuple]>,
    fastvep: bool,
) -> Result<OutputRows, Failure> {
    let mut output = OutputRows::default();
    if !path.exists() {
        return Ok(output);
    }
    let file = File::open(path).map_err(|error| io_failure(path, error))?;
    let mut index: HashMap<Tuple, usize> = HashMap::new();
    let mut prep = KeyPrep::default();
    let (mut location, mut consequence) = (Vec::new(), Vec::new());
    let mut line_number = 0u64;
    let mut not_utf8: Option<u64> = None;
    let read = for_each_record(file, |record| {
        line_number += 1;
        let Ok(line) = std::str::from_utf8(record) else {
            not_utf8 = Some(line_number);
            return Err(io::Error::other("a line that is not UTF-8"));
        };
        if line.is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let mut fields = line.split('\t');
        let mut columns = [""; 9];
        let mut held = 0;
        for column in fields.by_ref().take(9) {
            columns[held] = column;
            held += 1;
        }
        if held < 7 {
            return Ok(());
        }
        location.clear();
        push_normalized_location(columns[1], &mut location);
        let Some(at_location) = wanted.get(location.as_slice()) else {
            return Ok(());
        };
        let allele = normalize_allele_of(columns[2], fastvep);
        consequence.clear();
        prep.normalized_set(columns[6], &mut consequence);
        let carried = at_location.iter().copied().find(|tuple| {
            let [_, wanted_allele, feature, feature_type, set] = strings.resolve(*tuple);
            wanted_allele == allele
                && feature == columns[4]
                && feature_type == columns[5]
                && set.as_bytes() == consequence.as_slice()
        });
        let Some(tuple) = carried else {
            return Ok(());
        };
        if held < 9 {
            output.short.insert(tuple);
            return Ok(());
        }
        let mut cells = [""; 3];
        let mut further = 0;
        for (cell, field) in cells.iter_mut().zip(&mut fields) {
            *cell = field;
            further += 1;
        }
        let row = pair_row(&columns);
        let row = if further == cells.len() {
            row.with_cells(cells)
        } else {
            row
        };
        let at = *index.entry(tuple).or_insert_with(|| {
            output.rows.push((tuple, Vec::new()));
            output.rows.len() - 1
        });
        let rows = &mut output.rows[at].1;
        if !rows.contains(&row) {
            rows.push(row);
        }
        Ok(())
    });
    match (not_utf8, read) {
        (Some(line), _) => Err(Failure::decode(format!(
            "{}: line {line} is not UTF-8",
            path.display()
        ))),
        (None, Err(error)) => Err(io_failure(path, error)),
        (None, Ok(())) => Ok(output),
    }
}

/// The rows of the outputs `paths` behind each tuple of `wanted`, for the row-aware pair
/// rules: every row the comparison reads into the tuple, its first nine columns as
/// `pair_row` reads them and its Protein_position, Amino_acids and Codons where it has
/// twelve, each distinct row once. A tuple one of whose rows has fewer than nine columns,
/// and so no position column, is left out. No output is read when `wanted` is empty.
///
/// `input_by_output` names the input VCF an output was written from. A row of such an
/// output whose first column spells no alleles is given the reference bases its allele
/// stands in for, read by position from the records of that input, which is opened only
/// where a row asks. The rows of each output are read against its own input before the rows
/// of all outputs are pooled, so a row that two outputs carry stands once for each distinct
/// reading, whatever the order of `paths`.
///
/// The outputs are read side by side; a failure is the one the outputs and inputs read one
/// after the other in the order of `paths` would end in first.
pub fn pair_rows_by_tuple(
    strings: &Strings,
    paths: &[PathBuf],
    wanted: &TupleSet,
    engine: &str,
    input_by_output: &BTreeMap<PathBuf, PathBuf>,
) -> Result<HashMap<Tuple, Vec<PairRow>>, Failure> {
    let mut found: HashMap<Tuple, Vec<PairRow>> = HashMap::new();
    if wanted.is_empty() {
        return Ok(found);
    }
    let mut by_location: HashMap<&[u8], &[Tuple]> = HashMap::new();
    let mut rest = wanted.as_slice();
    while let Some(first) = rest.first() {
        let run = rest
            .iter()
            .take_while(|tuple| tuple.location() == first.location())
            .count();
        by_location.insert(strings.text(first.location()).as_bytes(), &rest[..run]);
        rest = &rest[run..];
    }
    let fastvep = is_fastvep(engine);
    let (outputs, failure) = crew_indexed_prefix(paths.len(), |at| {
        read_output_rows(strings, &paths[at], &by_location, fastvep)
    });
    let mut short: HashSet<Tuple> = HashSet::new();
    for (path, mut output) in paths.iter().zip(outputs) {
        short.extend(output.short);
        if let Some(input) = input_by_output.get(path) {
            let mut asking: Vec<&mut PairRow> = output
                .rows
                .iter_mut()
                .filter(|(tuple, _)| !short.contains(tuple))
                .flat_map(|(_, rows)| rows.iter_mut())
                .collect();
            read_record_references(input, &mut asking).map_err(|error| io_failure(input, error))?;
        }
        for (tuple, rows) in output.rows {
            let pooled = found.entry(tuple).or_default();
            for row in rows {
                if !pooled.contains(&row) {
                    pooled.push(row);
                }
            }
        }
    }
    if let Some(failure) = failure {
        return Err(failure);
    }
    found.retain(|tuple, _| !short.contains(tuple));
    Ok(found)
}

fn count_under(by_bucket: &mut NamedCounts, bucket: &str) {
    if let Some((_, n)) = by_bucket.iter_mut().find(|(name, _)| name == bucket) {
        *n += 1;
    }
}

/// The pair mask over any rules: `excluding` are asked first on the two consequence sets of
/// a pair, the first match taking it and setting both members aside. A pair none of them
/// takes is put to `row_pair_rules` in their order, the first match taking it, on the rows
/// of `perl_paths` and `rust_paths` that carry its two tuples ([`pair_rows_by_tuple`]); such
/// a rule counts the pair and sets it aside where it excludes. Without a row rule, or
/// without a pair left for one, no output is read.
///
/// `sources` are the candidate's cache and the reference sequence the comparison was given.
/// A row rule that asks the coding ends of a pair's transcript is asked with the reading of
/// the transcript's cache record, looked up by the pair's Location and Feature, the
/// transcripts of a sequence read when a pair on it first asks. Beside a cache the rule of
/// `start_cooccurrence_swap` takes a pair only on a transcript whose record shows a start
/// codon (`start_codon_stands`); a pair it leaves is put to the row rules. A cache that
/// cannot be read is a failure.
///
/// Returns `(excluded_rust, excluded_perl, matches_by_bucket)`, the buckets in the order of
/// the two lists, each once.
#[allow(clippy::too_many_arguments)]
pub fn filter_swap_pairs(
    sides: &Sides<'_>,
    excluding: &[&dyn PairRule],
    row_pair_rules: &[&dyn PairRule],
    perl_paths: &[PathBuf],
    rust_paths: &[PathBuf],
    scored_engine: &str,
    input_by_output: &BTreeMap<PathBuf, PathBuf>,
    sources: CodingEndSources<'_>,
) -> Result<(TupleSet, TupleSet, NamedCounts), Failure> {
    let mut by_bucket: NamedCounts = Vec::new();
    for rule in excluding.iter().chain(row_pair_rules) {
        if !by_bucket.iter().any(|(name, _)| name == rule.bucket()) {
            by_bucket.push((rule.bucket().to_string(), 0));
        }
    }
    let sets = |perl: Tuple, rust: Tuple| {
        let [_, allele, _, _, perl_csq] = sides.resolve(perl);
        (perl_csq, sides.strings.text(rust.consequence()), allele)
    };
    // What the sources hold for the pair of the reference's tuple.
    let read_for = |perl: Tuple| {
        let [location, _, feature, feature_type, _] = sides.resolve(perl);
        sources
            .of_pair(location, feature, feature_type)
            .map_err(|error| Failure::value(error.to_string()))
    };
    let (mut excluded_perl, mut excluded_rust) = (Vec::new(), Vec::new());
    let mut untaken: Vec<(Tuple, Tuple)> = Vec::new();
    let mut rest = sides.perl_only.as_slice();
    while let Some(&perl) = rest.first() {
        let key = perl.key4();
        let run = rest.iter().take_while(|tuple| tuple.key4() == key).count();
        rest = &rest[run..];
        let rust_side = sides.rust_only.with_prefix(&key);
        if run != 1 || rust_side.len() != 1 {
            continue;
        }
        let rust = rust_side[0];
        let (perl_csq, rust_csq, allele) = sets(perl, rust);
        // Beside a cache the rule of the start codon is asked the transcript's.
        let start_rule = excluding.iter().any(|rule| {
            rule.bucket() == START_COOCCURRENCE_SWAP && rule.matches(perl_csq, rust_csq, allele)
        });
        let start_stands = match sources.given() && start_rule {
            true => read_for(perl)?.ask(start_codon_stands),
            false => true,
        };
        match excluding.iter().find(|rule| {
            rule.matches(perl_csq, rust_csq, allele)
                && (rule.bucket() != START_COOCCURRENCE_SWAP || start_stands)
        }) {
            Some(rule) => {
                excluded_perl.push(perl);
                excluded_rust.push(rust);
                count_under(&mut by_bucket, rule.bucket());
            }
            None => untaken.push((perl, rust)),
        }
    }
    if !row_pair_rules.is_empty() && !untaken.is_empty() {
        let wanted_perl: TupleSet = untaken.iter().map(|(perl, _)| *perl).collect();
        let wanted_rust: TupleSet = untaken.iter().map(|(_, rust)| *rust).collect();
        let perl_rows = pair_rows_by_tuple(
            sides.strings,
            perl_paths,
            &wanted_perl,
            "vep-rs",
            input_by_output,
        )?;
        let rust_rows = pair_rows_by_tuple(
            sides.strings,
            rust_paths,
            &wanted_rust,
            scored_engine,
            input_by_output,
        )?;
        for (perl, rust) in untaken {
            let (Some(of_perl), Some(of_rust)) = (perl_rows.get(&perl), rust_rows.get(&rust))
            else {
                continue;
            };
            let (perl_csq, rust_csq, allele) = sets(perl, rust);
            let taken = read_for(perl)?.ask(|at| {
                row_pair_rules.iter().find(|rule| {
                    rule.matches_rows_at(of_perl, of_rust, at)
                        || rule.matches(perl_csq, rust_csq, allele)
                })
            });
            let Some(rule) = taken else {
                continue;
            };
            count_under(&mut by_bucket, rule.bucket());
            if rule.excludes() {
                excluded_perl.push(perl);
                excluded_rust.push(rust);
            }
        }
    }
    Ok((
        TupleSet::from_unsorted(excluded_rust),
        TupleSet::from_unsorted(excluded_perl),
        by_bucket,
    ))
}

/// Both members of every pair a registry rule takes, and the pairs matched per bucket in
/// registry order: `(excluded_rust, excluded_perl, matches_by_bucket)`. `sources` are the
/// candidate's cache and the reference sequence the comparison was given
/// ([`filter_swap_pairs`]).
pub fn filter_registry_swap_pairs(
    sides: &Sides<'_>,
    row_pair_rules: &[&'static ExclusionRule],
    perl_paths: &[PathBuf],
    rust_paths: &[PathBuf],
    scored_engine: &str,
    input_by_output: &BTreeMap<PathBuf, PathBuf>,
    sources: CodingEndSources<'_>,
) -> Result<(TupleSet, TupleSet, NamedCounts), Failure> {
    let as_asked = |rule: &&'static ExclusionRule| *rule as &dyn PairRule;
    let excluding: Vec<&dyn PairRule> = excluding_rules().iter().map(as_asked).collect();
    let row_pair_rules: Vec<&dyn PairRule> = row_pair_rules.iter().map(as_asked).collect();
    filter_swap_pairs(
        sides,
        &excluding,
        &row_pair_rules,
        perl_paths,
        rust_paths,
        scored_engine,
        input_by_output,
        sources,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::CacheAuthority;
    use crate::failure::FailureKind;
    use crate::masks::{apply_masks, MaskContext, MaskResult, MAX_SV_SIZE};
    use crate::normalize::{normalize_allele, normalize_consequence_set, normalize_location};
    use crate::rules::cache_transcripts::CacheTranscripts;
    use crate::rules::row_pair_rules;
    use crate::strings::OwnedSides;
    use crate::tempdir::TempDir;
    use std::cell::RefCell;
    use std::collections::BTreeSet;
    use std::fmt;

    type Texts = [String; 5];

    /// One annotation tuple on a transcript, as the comparison normalises it.
    fn t(loc: &str, allele: &str, feature: &str, csq: &str) -> Texts {
        [
            normalize_location(loc),
            normalize_allele(allele, "vep-rs"),
            feature.to_string(),
            "Transcript".to_string(),
            normalize_consequence_set(csq),
        ]
    }

    fn rows_of(texts: &[Texts]) -> Vec<[&str; 5]> {
        texts
            .iter()
            .map(|row| row.each_ref().map(String::as_str))
            .collect()
    }

    fn sides_of(perl: &[Texts], rust: &[Texts]) -> OwnedSides {
        OwnedSides::from_rows(&rows_of(perl), &rows_of(rust))
    }

    fn no_inputs() -> BTreeMap<PathBuf, PathBuf> {
        BTreeMap::new()
    }

    /// The registry's mask on one pair of tuple sets with no output behind them.
    fn pairs_of(perl: Texts, rust: Texts) -> (OwnedSides, TupleSet, TupleSet, NamedCounts) {
        let owned = sides_of(&[perl], &[rust]);
        let (excluded_rust, excluded_perl, buckets) = filter_registry_swap_pairs(
            &owned.sides(),
            &[],
            &[],
            &[],
            "vep-rs",
            &no_inputs(),
            CodingEndSources::default(),
        )
        .unwrap();
        (owned, excluded_rust, excluded_perl, buckets)
    }

    fn bucket(buckets: &NamedCounts, name: &str) -> u64 {
        buckets
            .iter()
            .find(|(held, _)| held == name)
            .unwrap_or_else(|| panic!("no bucket {name}"))
            .1
    }

    fn names(buckets: &NamedCounts) -> Vec<&str> {
        buckets.iter().map(|(name, _)| name.as_str()).collect()
    }

    #[test]
    fn test_structural_start_co_emission_pair_is_excluded_both_sides() {
        let loc = "21:26965975-26966647";
        let (owned, excluded_rust, excluded_perl, buckets) = pairs_of(
            t(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost,start_retained_variant"),
            t(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost"),
        );
        assert_eq!((&excluded_rust, &excluded_perl), (&owned.rust, &owned.perl));
        assert_eq!(bucket(&buckets, "start_cooccurrence_swap"), 1);
        assert_eq!(bucket(&buckets, "covered_splice_region_swap"), 0);
    }

    #[test]
    fn test_sequence_start_co_emission_pair_drops_the_other_member() {
        let (owned, excluded_rust, excluded_perl, _) = pairs_of(
            t(
                "21:36265246",
                "-",
                "ENST00000486278",
                "frameshift_variant,start_lost,start_retained_variant",
            ),
            t(
                "21:36265246",
                "-",
                "ENST00000486278",
                "frameshift_variant,start_retained_variant",
            ),
        );
        assert_eq!((&excluded_rust, &excluded_perl), (&owned.rust, &owned.perl));
    }

    /// A structural allele on which the scored engine dropped `start_lost` instead is its
    /// own divergence, not the co-emission shape.
    #[test]
    fn test_the_wrong_member_for_the_allele_kind_is_not_excluded() {
        let loc = "21:26965975-26966647";
        let (_, excluded_rust, excluded_perl, _) = pairs_of(
            t(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost,start_retained_variant"),
            t(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_retained_variant"),
        );
        assert!(excluded_rust.is_empty() && excluded_perl.is_empty());
    }

    #[test]
    fn test_both_start_terms_missing_on_vep_rs_is_not_the_shape() {
        let loc = "21:39660415-39660694";
        let (_, excluded_rust, excluded_perl, _) = pairs_of(
            t(loc, "Alu_insertion", "ENST00000852645", "5_prime_UTR_variant,coding_sequence_variant,feature_elongation,start_lost,start_retained_variant"),
            t(loc, "Alu_insertion", "ENST00000852645", "5_prime_UTR_variant,coding_sequence_variant,feature_elongation"),
        );
        assert!(excluded_rust.is_empty() && excluded_perl.is_empty());
    }

    #[test]
    fn test_covered_splice_region_pair_is_excluded_both_sides() {
        let loc = "21:31663791-31663794";
        let (owned, excluded_rust, excluded_perl, buckets) = pairs_of(
            t(loc, "CA", "ENST00000270142", "frameshift_variant"),
            t(
                loc,
                "CA",
                "ENST00000270142",
                "frameshift_variant,splice_region_variant",
            ),
        );
        assert_eq!((&excluded_rust, &excluded_perl), (&owned.rust, &owned.perl));
        assert_eq!(bucket(&buckets, "covered_splice_region_swap"), 1);
    }

    #[test]
    fn test_a_polypyrimidine_tract_pair_is_not_the_covered_shape() {
        let loc = "21:29079809-29085808";
        let (_, excluded_rust, excluded_perl, _) = pairs_of(
            t(
                loc,
                "deletion",
                "ENST00000341618",
                "intron_variant,splice_polypyrimidine_tract_variant",
            ),
            t(loc, "deletion", "ENST00000341618", "intron_variant"),
        );
        assert!(excluded_rust.is_empty() && excluded_perl.is_empty());
    }

    #[test]
    fn test_an_extra_vep_rs_term_beside_the_co_emission_rejects_the_pair() {
        let loc = "21:26965975-26966647";
        let (_, excluded_rust, excluded_perl, _) = pairs_of(
            t(
                loc,
                "deletion",
                "ENST00000284987",
                "feature_truncation,start_lost,start_retained_variant",
            ),
            t(
                loc,
                "deletion",
                "ENST00000284987",
                "feature_truncation,start_lost,intron_variant",
            ),
        );
        assert!(excluded_rust.is_empty() && excluded_perl.is_empty());
    }

    #[test]
    fn test_one_sided_tuples_are_never_pairs() {
        let (_, excluded_rust, excluded_perl, _) = pairs_of(
            t(
                "21:1-2",
                "deletion",
                "ENST_A",
                "start_lost,start_retained_variant",
            ),
            t("21:1-2", "deletion", "ENST_B", "start_lost"),
        );
        assert!(excluded_rust.is_empty() && excluded_perl.is_empty());
    }

    const ROW_COLS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

    /// One fabricated row of the default format; it prints as its line.
    #[derive(Clone)]
    struct Row {
        allele: String,
        consequence: String,
        cdna: String,
        cds: String,
        location: String,
        name: String,
        feature: String,
        amino_acids: String,
        codons: String,
    }

    fn stop_codon_row(allele: &str, consequence: &str, cdna: &str, cds: &str) -> Row {
        Row {
            allele: allele.to_string(),
            consequence: consequence.to_string(),
            cdna: cdna.to_string(),
            cds: cds.to_string(),
            location: "19:13208046-13208050".to_string(),
            name: "1993671".to_string(),
            feature: "ENST00000635895".to_string(),
            amino_acids: "-".to_string(),
            codons: "-".to_string(),
        }
    }

    impl Row {
        fn at(mut self, location: &str) -> Row {
            self.location = location.to_string();
            self
        }

        fn named(mut self, name: &str) -> Row {
            self.name = name.to_string();
            self
        }

        fn on(mut self, feature: &str) -> Row {
            self.feature = feature.to_string();
            self
        }

        /// The row with its Amino_acids and its Codons.
        fn coding(mut self, amino_acids: &str, codons: &str) -> Row {
            self.amino_acids = amino_acids.to_string();
            self.codons = codons.to_string();
            self
        }
    }

    impl fmt::Display for Row {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            writeln!(
                f,
                "{}\t{}\t{}\tENSG1\t{}\tTranscript\t{}\t{}\t{}\t-\t{}\t{}\t-\tIMPACT=LOW;STRAND=-1",
                self.name,
                self.location,
                self.allele,
                self.feature,
                self.consequence,
                self.cdna,
                self.cds,
                self.amino_acids,
                self.codons
            )
        }
    }

    /// The tuples of the rows of an output's text, as the comparison reads them.
    fn tuples_of(text: &str) -> Vec<Texts> {
        text.lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                [
                    normalize_location(columns[1]),
                    normalize_allele(columns[2], "vep-rs"),
                    columns[4].to_string(),
                    columns[5].to_string(),
                    normalize_consequence_set(columns[6]),
                ]
            })
            .collect()
    }

    /// One fabricated output pair on disk: both paths and both tuple sets.
    struct Outputs {
        folder: TempDir,
        perl_path: PathBuf,
        rust_path: PathBuf,
        perl: Vec<Texts>,
        rust: Vec<Texts>,
        owned: OwnedSides,
    }

    fn row_outputs_under(perl: &str, rust: &str, cols: &str) -> Outputs {
        let folder = TempDir::new("vep_sv_compare_pairs_").unwrap();
        let write = |side: &str, rows: &str| {
            let path = folder.path.join(format!("{side}.txt"));
            std::fs::write(&path, format!("{cols}{rows}")).unwrap();
            path
        };
        let (perl_path, rust_path) = (write("perl", perl), write("rust", rust));
        let (perl, rust) = (tuples_of(perl), tuples_of(rust));
        let owned = sides_of(&perl, &rust);
        Outputs {
            folder,
            perl_path,
            rust_path,
            perl,
            rust,
            owned,
        }
    }

    fn row_outputs(perl: &str, rust: &str) -> Outputs {
        row_outputs_under(perl, rust, ROW_COLS)
    }

    impl Outputs {
        /// Both outputs written from the input `input`.
        fn inputs(&self, input: &Path) -> BTreeMap<PathBuf, PathBuf> {
            [&self.perl_path, &self.rust_path]
                .into_iter()
                .map(|output| (output.clone(), input.to_path_buf()))
                .collect()
        }

        /// The registry's mask with `rules` as its row-aware rules.
        fn registry_mask(
            &self,
            rules: &[&'static ExclusionRule],
            inputs: &BTreeMap<PathBuf, PathBuf>,
        ) -> Result<(TupleSet, TupleSet, NamedCounts), Failure> {
            filter_registry_swap_pairs(
                &self.owned.sides(),
                rules,
                std::slice::from_ref(&self.perl_path),
                std::slice::from_ref(&self.rust_path),
                "vep-rs",
                inputs,
                CodingEndSources::default(),
            )
        }

        /// The mask with one fabricated rule after the registry's excluding rules.
        fn fabricated_mask(&self, rule: &dyn PairRule) -> (TupleSet, TupleSet, NamedCounts) {
            let excluding = excluding_rules();
            let excluding: Vec<&dyn PairRule> = excluding
                .iter()
                .map(|rule| *rule as &dyn PairRule)
                .collect();
            filter_swap_pairs(
                &self.owned.sides(),
                &excluding,
                &[rule],
                std::slice::from_ref(&self.perl_path),
                std::slice::from_ref(&self.rust_path),
                "vep-rs",
                &no_inputs(),
                CodingEndSources::default(),
            )
            .unwrap()
        }
    }

    /// A fabricated row-aware rule around a predicate over the two rows of a pair.
    struct Fabricated<F: Fn(&PairRow, &PairRow) -> bool> {
        predicate: F,
    }

    impl<F: Fn(&PairRow, &PairRow) -> bool> PairRule for Fabricated<F> {
        fn bucket(&self) -> &str {
            "fabricated_row_shape"
        }

        fn excludes(&self) -> bool {
            true
        }

        fn matches(&self, _perl_csq: &str, _rust_csq: &str, _allele: &str) -> bool {
            false
        }

        fn matches_rows_at(
            &self,
            perl: &[PairRow],
            rust: &[PairRow],
            _at: &AtACodingEnd<'_>,
        ) -> bool {
            !perl.is_empty()
                && !rust.is_empty()
                && perl
                    .iter()
                    .all(|e| rust.iter().all(|a| (self.predicate)(e, a)))
        }
    }

    fn masks_of(
        release: &str,
        owned: &OwnedSides,
        perl_paths: &[PathBuf],
        rust_paths: &[PathBuf],
        inputs: &BTreeMap<PathBuf, PathBuf>,
    ) -> MaskResult {
        apply_masks(
            &owned.sides(),
            &MaskContext {
                reference_release: release,
                cache: &CacheAuthority::default(),
                scored_engine: "vep-rs",
                input_variants: Some(&[]),
                max_sv_size: MAX_SV_SIZE,
                perl_paths,
                rust_paths,
                input_by_output: inputs,
                coding_ends: None,
            },
        )
        .unwrap()
    }

    fn features(owned: &OwnedSides, set: &TupleSet) -> BTreeSet<String> {
        set.iter()
            .map(|tuple| owned.strings.text(tuple.feature()).to_string())
            .collect()
    }

    fn locations(owned: &OwnedSides, set: &TupleSet) -> Vec<String> {
        set.iter()
            .map(|tuple| owned.strings.text(tuple.location()).to_string())
            .collect()
    }

    #[test]
    fn test_a_row_aware_rule_is_handed_the_named_columns_of_the_two_rows_of_an_unmatched_pair() {
        let seen: RefCell<Vec<(PairRow, PairRow)>> = RefCell::new(Vec::new());
        let name = "19_13208045_GACCTG/G";
        let outputs = row_outputs(
            &stop_codon_row("-", "stop_retained_variant", "?-7057", "-")
                .named(name)
                .to_string(),
            &stop_codon_row("-", "stop_lost", "?-7057/7140", "-/7023")
                .named(name)
                .to_string(),
        );
        let rule = Fabricated {
            predicate: |e: &PairRow, a: &PairRow| {
                seen.borrow_mut().push((e.clone(), a.clone()));
                true
            },
        };
        let (excluded_rust, excluded_perl, buckets) = outputs.fabricated_mask(&rule);
        let seen = seen.into_inner();
        let [(e, a)] = &seen[..] else {
            panic!("the rule was asked {} times", seen.len());
        };
        assert_eq!(
            (
                e.location.as_str(),
                e.allele.as_str(),
                e.feature_type.as_str(),
                e.consequence.as_str(),
                e.cdna_position.as_str(),
                e.cds_position.as_str(),
                e.indel
            ),
            (
                "19:13208046-13208050",
                "-",
                "Transcript",
                "stop_retained_variant",
                "?-7057",
                "-",
                Some(true)
            )
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
            (&excluded_rust, &excluded_perl),
            (&outputs.owned.rust, &outputs.owned.perl)
        );
        assert_eq!(bucket(&buckets, "fabricated_row_shape"), 1);
        assert_eq!(names(&buckets).last(), Some(&"fabricated_row_shape"));
    }

    #[test]
    fn test_a_row_aware_rule_is_handed_the_three_cells_after_the_positions_of_a_row_that_has_them()
    {
        type Cells = [Option<String>; 3];
        let seen: RefCell<Vec<(Cells, Cells)>> = RefCell::new(Vec::new());
        let row_of = |consequence: &str, cells: &str| {
            format!(
                "1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{consequence}\t?-7057\t-{cells}\n"
            )
        };
        let outputs = row_outputs(
            &row_of(
                "stop_retained_variant",
                "\t2341-?\tKL/X\taagCTg/aag\t-\tIMPACT=LOW",
            ),
            &row_of("stop_lost", "\t2341-?\tKL/X"),
        );
        let cells = |row: &PairRow| {
            [
                row.protein_position.clone(),
                row.amino_acids.clone(),
                row.codons.clone(),
            ]
        };
        let rule = Fabricated {
            predicate: |e: &PairRow, a: &PairRow| {
                seen.borrow_mut().push((cells(e), cells(a)));
                true
            },
        };
        let (_, _, buckets) = outputs.fabricated_mask(&rule);
        let text = |cell: &str| Some(cell.to_string());
        assert_eq!(
            seen.into_inner(),
            [(
                [text("2341-?"), text("KL/X"), text("aagCTg/aag")],
                [None, None, None]
            )],
            "none of a row that ends before its twelfth column"
        );
        assert_eq!(bucket(&buckets, "fabricated_row_shape"), 1);
    }

    #[test]
    fn test_only_an_unmatched_pair_no_earlier_rule_takes_reaches_a_row_aware_rule() {
        let asked = RefCell::new(0);
        let matched = stop_codon_row("-", "missense_variant", "10-12", "4-6").on("ENST00000000001");
        let outputs = row_outputs(
            &format!(
                "{matched}{}{}",
                stop_codon_row("-", "stop_retained_variant", "20-22", "14-16")
                    .on("ENST00000000002"),
                stop_codon_row("-", "intron_variant", "-", "-").on("ENST00000000004"),
            ),
            &format!(
                "{matched}{}{}",
                stop_codon_row("-", "stop_lost", "20-22", "14-16").on("ENST00000000003"),
                stop_codon_row("-", "intron_variant,splice_region_variant", "-", "-")
                    .on("ENST00000000004"),
            ),
        );
        let rule = Fabricated {
            predicate: |_: &PairRow, _: &PairRow| {
                *asked.borrow_mut() += 1;
                true
            },
        };
        let (excluded_rust, excluded_perl, buckets) = outputs.fabricated_mask(&rule);
        assert_eq!(
            asked.into_inner(),
            0,
            "a matched tuple, a tuple one engine alone writes and a pair an earlier rule takes"
        );
        assert_eq!(
            (
                bucket(&buckets, "fabricated_row_shape"),
                bucket(&buckets, "covered_splice_region_swap")
            ),
            (0, 1)
        );
        assert_eq!(
            features(&outputs.owned, &excluded_rust.union(&excluded_perl)),
            BTreeSet::from(["ENST00000000004".to_string()])
        );
    }

    #[test]
    fn test_the_rows_of_one_tuple_stand_or_fall_together() {
        let outputs = row_outputs(
            &format!(
                "{}{}{}{}",
                stop_codon_row("-", "stop_retained_variant", "10-12", "4-6").named("first"),
                stop_codon_row("-", "stop_retained_variant", "20-22", "14-16").named("second"),
                stop_codon_row("-", "stop_retained_variant", "10-12", "4-6").named("third"),
                stop_codon_row("-", "stop_retained_variant", "10-12", "4-6").named("first"),
            ),
            &stop_codon_row("-", "stop_lost", "10-12", "4-6")
                .named("first")
                .to_string(),
        );
        for (refused, expected, asked) in [(None, 1, 3), (Some("20-22"), 0, 2)] {
            let seen: RefCell<Vec<String>> = RefCell::new(Vec::new());
            let rule = Fabricated {
                predicate: |e: &PairRow, _: &PairRow| {
                    seen.borrow_mut().push(e.cdna_position.clone());
                    Some(e.cdna_position.as_str()) != refused
                },
            };
            let (excluded_rust, excluded_perl, buckets) = outputs.fabricated_mask(&rule);
            assert_eq!(
                seen.into_inner(),
                ["10-12", "20-22", "10-12"][..asked],
                "each distinct row of the reference's tuple, its first column among its cells, \
                 once, up to the row the rule refuses"
            );
            assert_eq!(
                (
                    bucket(&buckets, "fabricated_row_shape"),
                    excluded_rust.len(),
                    excluded_perl.len()
                ),
                (expected, expected as usize, expected as usize),
                "{refused:?}"
            );
        }
    }

    #[test]
    fn test_a_row_of_fewer_than_nine_columns_is_put_to_no_row_aware_rule() {
        let asked = RefCell::new(0);
        let short = |consequence: &str| {
            format!("1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{consequence}\t?-7057\n")
        };
        let cols: Vec<&str> = ROW_COLS.split('\t').take(8).collect();
        let outputs = row_outputs_under(
            &short("stop_retained_variant"),
            &short("stop_lost"),
            &format!("{}\n", cols.join("\t")),
        );
        let rule = Fabricated {
            predicate: |_: &PairRow, _: &PairRow| {
                *asked.borrow_mut() += 1;
                true
            },
        };
        let (excluded_rust, excluded_perl, buckets) = outputs.fabricated_mask(&rule);
        assert_eq!(
            (
                asked.into_inner(),
                excluded_rust.len(),
                excluded_perl.len(),
                bucket(&buckets, "fabricated_row_shape")
            ),
            (0, 0, 0, 0)
        );
    }

    /// Without a row-aware rule the outputs are not opened: here neither can be read, and a
    /// row-aware rule that asks for their rows ends the comparison.
    #[test]
    fn test_without_a_row_aware_rule_no_output_is_read_again() {
        let outputs = row_outputs(
            &stop_codon_row("-", "stop_retained_variant", "?-7057", "-").to_string(),
            &stop_codon_row("-", "stop_lost", "?-7057", "-").to_string(),
        );
        for path in [&outputs.perl_path, &outputs.rust_path] {
            std::fs::write(path, b"\xff\n").unwrap();
        }
        let (excluded_rust, excluded_perl, buckets) =
            outputs.registry_mask(&[], &no_inputs()).unwrap();
        let without_paths = filter_registry_swap_pairs(
            &outputs.owned.sides(),
            &[],
            &[],
            &[],
            "vep-rs",
            &no_inputs(),
            CodingEndSources::default(),
        )
        .unwrap();
        assert_eq!(without_paths.2, buckets);
        assert!(excluded_rust.is_empty() && excluded_perl.is_empty());
        assert_eq!(
            names(&buckets),
            excluding_rules()
                .iter()
                .map(|rule| rule.bucket)
                .collect::<Vec<_>>()
        );
        let read = outputs
            .registry_mask(&row_pair_rules(), &no_inputs())
            .unwrap_err();
        assert_eq!(read.kind, FailureKind::Decode);
    }

    /// The masks of 116.2 hand the registry's row-aware rules to the pair mask and those of
    /// 115.2 none; a rule handed to the pair mask is asked once for the pair.
    #[test]
    fn test_the_masks_of_a_release_ask_its_row_aware_rules() {
        let outputs = row_outputs(
            &stop_codon_row("-", "stop_retained_variant", "?-7057", "-").to_string(),
            &stop_codon_row("-", "stop_lost", "?-7057", "-").to_string(),
        );
        let (perl_paths, rust_paths) = (
            std::slice::from_ref(&outputs.perl_path),
            std::slice::from_ref(&outputs.rust_path),
        );
        let masks = masks_of(
            "116.2",
            &outputs.owned,
            perl_paths,
            rust_paths,
            &no_inputs(),
        );
        assert_eq!(
            (
                &masks.swap_rust,
                &masks.swap_perl,
                bucket(
                    &masks.swap_by_bucket,
                    "stop_codon_removed_without_stop_codon_test"
                )
            ),
            (&outputs.owned.rust, &outputs.owned.perl, 1)
        );
        let asked = RefCell::new(0);
        let rule = Fabricated {
            predicate: |_: &PairRow, _: &PairRow| {
                *asked.borrow_mut() += 1;
                true
            },
        };
        assert_eq!(
            bucket(&outputs.fabricated_mask(&rule).2, "fabricated_row_shape"),
            1
        );
        assert_eq!(asked.into_inner(), 1);
        let published = masks_of(
            "115.2",
            &outputs.owned,
            perl_paths,
            rust_paths,
            &no_inputs(),
        );
        assert!(published.swap_rust.is_empty() && published.swap_perl.is_empty());
        assert_eq!(
            names(&published.swap_by_bucket),
            ["start_cooccurrence_swap", "covered_splice_region_swap"],
            "no bucket of a row-aware rule"
        );
        let with_rows: Vec<&str> = row_pair_rules().iter().map(|rule| rule.bucket).collect();
        assert_eq!(names(&masks.swap_by_bucket)[2..], with_rows[..]);
    }

    /// One row-aware rule of the 116.2 registry on fabricated rows: the reference's row and
    /// the scored engine's row of the rule's pair, the scored row with one further differing
    /// term, and the pair with a position column outside the rule's condition.
    struct RuleCase {
        bucket: &'static str,
        reference: Row,
        candidate: Row,
        further_term: Row,
        other_column: (Row, Row),
    }

    /// The registry pair mask of the 116.2 release over one fabricated output pair: the
    /// rule's count, and the Locations of the tuples of each side the mask leaves charged.
    fn masked(
        bucket_name: &str,
        perl_rows: &str,
        rust_rows: &str,
    ) -> (u64, Vec<String>, Vec<String>) {
        let matched = stop_codon_row("-", "missense_variant", "10-12", "4-6")
            .at("21:5000-5002")
            .on("ENST00000000001");
        let outputs = row_outputs(
            &format!("{matched}{perl_rows}"),
            &format!("{matched}{rust_rows}"),
        );
        let (swap_rust, swap_perl, buckets) = outputs
            .registry_mask(&row_pair_rules(), &no_inputs())
            .unwrap();
        let sides = outputs.owned.sides();
        (
            bucket(&buckets, bucket_name),
            locations(&outputs.owned, &sides.perl_only.difference(&swap_perl)),
            locations(&outputs.owned, &sides.rust_only.difference(&swap_rust)),
        )
    }

    fn the_pair_of_rows_is_set_aside(case: &RuleCase) {
        assert_eq!(
            masked(
                case.bucket,
                &case.reference.to_string(),
                &case.candidate.to_string()
            ),
            (1, Vec::new(), Vec::new())
        );
    }

    fn one_further_differing_term_stays_a_difference(case: &RuleCase) {
        let (count, perl_open, rust_open) = masked(
            case.bucket,
            &case.reference.to_string(),
            &case.further_term.to_string(),
        );
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
    }

    fn a_column_outside_the_condition_stays_a_difference(case: &RuleCase) {
        let (count, perl_open, rust_open) = masked(
            case.bucket,
            &case.other_column.0.to_string(),
            &case.other_column.1.to_string(),
        );
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
    }

    fn a_missing_tuple_and_a_wrong_tuple_stay_differences(case: &RuleCase) {
        let missing = stop_codon_row("-", "intron_variant", "-", "-")
            .at("21:5000-5002")
            .on("ENST00000000002");
        let wrong = "21:7000-7002";
        let (count, perl_open, rust_open) = masked(
            case.bucket,
            &format!(
                "{}{missing}{}",
                case.reference,
                stop_codon_row("-", "missense_variant", "30-32", "24-26").at(wrong)
            ),
            &format!(
                "{}{}",
                case.candidate,
                stop_codon_row("-", "synonymous_variant", "30-32", "24-26").at(wrong)
            ),
        );
        assert_eq!(count, 1);
        assert_eq!(perl_open, ["21:5000-5002", wrong]);
        assert_eq!(rust_open, [wrong]);
    }

    macro_rules! row_pair_rule_checks {
        ($class:ident, $case:expr) => {
            mod $class {
                #[test]
                fn test_the_pair_of_rows_is_set_aside() {
                    super::the_pair_of_rows_is_set_aside(&$case);
                }

                #[test]
                fn test_one_further_differing_term_stays_a_difference() {
                    super::one_further_differing_term_stays_a_difference(&$case);
                }

                #[test]
                fn test_a_column_outside_the_condition_stays_a_difference() {
                    super::a_column_outside_the_condition_stays_a_difference(&$case);
                }

                #[test]
                fn test_a_missing_tuple_and_a_wrong_tuple_stay_differences() {
                    super::a_missing_tuple_and_a_wrong_tuple_stay_differences(&$case);
                }
            }
        };
    }

    /// The reference's coding_sequence_variant against the scored engine's
    /// stop_retained_variant beside splice_region_variant, on a row whose stop codon test
    /// is made.
    fn stop_codon_overlap_genomic_span() -> RuleCase {
        let row = |term: &str, cdna: &str, cds: &str| {
            stop_codon_row(
                "-",
                &format!("splice_region_variant,{term},NMD_transcript_variant"),
                cdna,
                cds,
            )
            .at("3:36993663")
        };
        RuleCase {
            bucket: "stop_codon_overlap_genomic_span",
            reference: row("coding_sequence_variant", "1705", "1630"),
            candidate: row("stop_retained_variant", "1705", "1630"),
            further_term: row("stop_retained_variant,stop_lost", "1705", "1630"),
            other_column: (
                row("coding_sequence_variant", "1705-?", "1630-?"),
                row("stop_retained_variant", "1705-?", "1630-?"),
            ),
        }
    }

    /// The reference's stop_retained_variant against nothing, coding_sequence_variant or
    /// frameshift_variant, on rows whose stop codon test is not made.
    fn stop_retained_without_stop_codon_test() -> RuleCase {
        let row = |term: &str, cdna: &str, cds: &str| {
            stop_codon_row("-", &format!("splice_acceptor_variant,{term}"), cdna, cds)
        };
        let reference = row("stop_retained_variant", "?-7057", "-");
        RuleCase {
            bucket: "stop_retained_without_stop_codon_test",
            reference: reference.clone(),
            candidate: row("coding_sequence_variant", "?-7057", "-"),
            further_term: row("coding_sequence_variant,frameshift_variant", "?-7057", "-"),
            other_column: (
                reference,
                row("coding_sequence_variant", "7055-7059", "6978-6982"),
            ),
        }
    }

    /// The reference's stop_retained_variant against coding_sequence_variant, on rows whose
    /// stop codon test is made.
    fn stop_retained_three_base_stop_codon_test() -> RuleCase {
        let row = |term: &str, cdna: &str| {
            stop_codon_row("-", &format!("{term},3_prime_UTR_variant"), cdna, "2249-?")
                .at("2:190970701-190970707")
        };
        let reference = row("stop_retained_variant", "2558-2564");
        RuleCase {
            bucket: "stop_retained_three_base_stop_codon_test",
            reference: reference.clone(),
            candidate: row("coding_sequence_variant", "2558-2564"),
            further_term: row("coding_sequence_variant,frameshift_variant", "2558-2564"),
            other_column: (reference, row("coding_sequence_variant", "2558-?")),
        }
    }

    /// The reference's stop_retained_variant against the scored engine's stop_lost, on rows
    /// without a cDNA position at one end or a CDS start.
    fn stop_codon_removed_without_stop_codon_test() -> RuleCase {
        let row = |term: &str, cdna: &str, cds: &str| {
            stop_codon_row("-", &format!("splice_acceptor_variant,{term}"), cdna, cds)
        };
        RuleCase {
            bucket: "stop_codon_removed_without_stop_codon_test",
            reference: row("stop_retained_variant", "?-7057", "-"),
            candidate: row("stop_lost", "?-7057", "-"),
            further_term: row("stop_lost,frameshift_variant", "?-7057", "-"),
            other_column: (
                row("stop_retained_variant", "7055-7059", "6978-6982"),
                row("stop_lost", "7055-7059", "6978-6982"),
            ),
        }
    }

    /// The reference's coding_sequence_variant against the scored engine's
    /// stop_retained_variant on an insertion's row with both cDNA positions and a CDS start.
    fn stop_codon_overlap_insertion_range() -> RuleCase {
        let row = |terms: &str, cdna: &str, cds: &str| {
            stop_codon_row("TTGA", terms, cdna, cds).at("19:13208046-13208047")
        };
        RuleCase {
            bucket: "stop_codon_overlap_insertion_range",
            reference: row(
                "inframe_insertion,coding_sequence_variant",
                "7054-7055",
                "6977-6978",
            ),
            candidate: row(
                "inframe_insertion,stop_retained_variant",
                "7054-7055",
                "6977-6978",
            ),
            further_term: row(
                "inframe_insertion,stop_retained_variant,stop_lost",
                "7054-7055",
                "6977-6978",
            ),
            other_column: (
                row(
                    "inframe_insertion,coding_sequence_variant",
                    "7054-?",
                    "6977-?",
                ),
                row(
                    "inframe_insertion,stop_retained_variant",
                    "7054-?",
                    "6977-?",
                ),
            ),
        }
    }

    /// The scored engine's stop_gained beside the frameshift_variant of both rows, on an
    /// insertion's row whose Codons cell settles a stop codon before its unknown bases: the
    /// reference's row writes no residues and the scored engine's residues hold their stop
    /// at that codon's place. The other column is a Codons cell with no stop codon before
    /// its unknown bases.
    fn stop_codon_settled_before_unknown_bases() -> RuleCase {
        let row = |terms: &str, amino_acids: &str, codons: &str| {
            stop_codon_row("GAATAGNNNN", terms, "775-776", "423-424")
                .at("21:10567700-10567701")
                .on("ENST00000622113")
                .coding(amino_acids, codons)
        };
        let (settled, unsettled) = ("-/GAATAGNNNN", "-/GAAGCTNNNN");
        let gained = "stop_gained,frameshift_variant";
        RuleCase {
            bucket: "stop_codon_settled_before_unknown_bases",
            reference: row("frameshift_variant", "-", settled),
            candidate: row(gained, "-/E*XX", settled),
            further_term: row(
                "stop_gained,frameshift_variant,stop_lost",
                "-/E*XX",
                settled,
            ),
            other_column: (
                row("frameshift_variant", "-", unsettled),
                row(gained, "-/E*XX", unsettled),
            ),
        }
    }

    /// The pair of `stop_codon_settled_before_unknown_bases` with residues on the scored
    /// engine's row whose first stop stands one place after the stop codon the Codons cell
    /// settles stays a difference.
    #[test]
    fn test_residues_with_their_stop_at_another_place_stay_a_difference() {
        let case = stop_codon_settled_before_unknown_bases();
        let other_place = case.candidate.clone().coding("-/EE*X", "-/GAATAGNNNN");
        let (count, perl_open, rust_open) = masked(
            case.bucket,
            &case.reference.to_string(),
            &other_place.to_string(),
        );
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
    }

    /// The scored engine's frameshift_variant in the place of the reference's
    /// coding_sequence_variant, on a deletion's row whose reference residues open with a
    /// stop: both rows hold the same residues, the unknown residue in their alternate half,
    /// and the same Codons cell, whose halves differ by one base. The other column is a
    /// Codons cell whose halves hold as many bases.
    fn frameshift_inside_a_held_stop_codon() -> RuleCase {
        let row = |terms: &str, codons: &str| {
            stop_codon_row("-", terms, "36", "30")
                .at("1:136")
                .on("ENST00000000010")
                .coding("*/X", codons)
        };
        RuleCase {
            bucket: "frameshift_inside_a_held_stop_codon",
            reference: row("coding_sequence_variant", "tgA/tg"),
            candidate: row("frameshift_variant", "tgA/tg"),
            further_term: row("stop_gained,frameshift_variant", "tgA/tg"),
            other_column: (
                row("coding_sequence_variant", "tGa/tAa"),
                row("frameshift_variant", "tGa/tAa"),
            ),
        }
    }

    /// The pair of `frameshift_inside_a_held_stop_codon` whose reference residues open with
    /// another residue than a stop stays a difference.
    #[test]
    fn test_residues_that_open_with_no_stop_stay_a_difference() {
        let case = frameshift_inside_a_held_stop_codon();
        let (count, perl_open, rust_open) = masked(
            case.bucket,
            &case.reference.clone().coding("L/X", "tgA/tg").to_string(),
            &case.candidate.clone().coding("L/X", "tgA/tg").to_string(),
        );
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
    }

    /// The reference's splice_donor_variant and intron_variant against the scored engine's
    /// splice_region_variant, on a row under an identifier, which carries no record: its
    /// three positions lie in the transcript's sequence under an allele of five.
    fn splice_term_past_the_record_span() -> RuleCase {
        let row = |terms: &str, cdna: &str, cds: &str| {
            stop_codon_row("GTTCA", terms, cdna, cds).at("7:100998-101000")
        };
        RuleCase {
            bucket: "splice_term_past_the_record_span",
            reference: row(
                "splice_donor_variant,frameshift_variant,intron_variant",
                "998-1000",
                "898-900",
            ),
            candidate: row(
                "frameshift_variant,splice_region_variant",
                "998-1000",
                "898-900",
            ),
            further_term: row(
                "stop_gained,frameshift_variant,splice_region_variant",
                "998-1000",
                "898-900",
            ),
            other_column: (
                row(
                    "splice_donor_variant,frameshift_variant,intron_variant",
                    "999-?",
                    "899-?",
                ),
                row("frameshift_variant,splice_region_variant", "999-?", "899-?"),
            ),
        }
    }

    /// The reference's synonymous_variant against the scored engine's missense_variant on a
    /// codon of a transcript on genetic code 26, told by the rows' Amino_acids and Codons;
    /// the other column is a Codons cell that holds no codon the code reads as alanine.
    fn pachysolen_nuclear_code_read_as_standard() -> RuleCase {
        let row = |term: &str, amino_acids: &str, codons: &str| {
            stop_codon_row("A", term, "114", "114")
                .at("7:1000")
                .coding(amino_acids, codons)
        };
        RuleCase {
            bucket: "pachysolen_nuclear_code_read_as_standard",
            reference: row("synonymous_variant", "L", "ctG/ctA"),
            candidate: row("missense_variant", "A/L", "ctG/ctA"),
            further_term: row("missense_variant,splice_region_variant", "A/L", "ctG/ctA"),
            other_column: (
                row("synonymous_variant", "L", "ctA/ctT"),
                row("missense_variant", "A/L", "ctA/ctT"),
            ),
        }
    }

    /// The scored engine's start_lost beside the inframe_deletion of both rows, on a row
    /// whose cDNA and CDS positions open at 1 and whose Amino_acids writes no alternate
    /// residue. The other column is a cDNA position that opens after a 5' UTR.
    fn start_codon_removed_without_five_prime_utr() -> RuleCase {
        let row = |terms: &str, cdna: &str| {
            stop_codon_row("-", terms, cdna, "1-3")
                .at("21:44817966-44817968")
                .named("21_44817965_ACAT/A")
                .on("ENST00000397893")
                .coding("M/-", "ATG/-")
        };
        let lost = "start_lost,inframe_deletion";
        RuleCase {
            bucket: "start_codon_removed_without_five_prime_utr",
            reference: row("inframe_deletion", "1-3"),
            candidate: row(lost, "1-3"),
            further_term: row("start_lost,inframe_deletion,stop_lost", "1-3"),
            other_column: (row("inframe_deletion", "101-103"), row(lost, "101-103")),
        }
    }

    /// The pair of `start_codon_removed_without_five_prime_utr` with an alternate residue
    /// in the Amino_acids of both rows stays a difference.
    #[test]
    fn test_residues_with_an_alternate_half_stay_a_difference() {
        let case = start_codon_removed_without_five_prime_utr();
        let kept = |row: &Row| row.clone().coding("MA/M", "ATGGCG/ATG").to_string();
        let (count, perl_open, rust_open) =
            masked(case.bucket, &kept(&case.reference), &kept(&case.candidate));
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
    }

    /// The scored engine's start_retained_variant against the reference's start_lost beside
    /// the frameshift_variant of both rows, on a row whose cDNA and CDS positions are one
    /// span that opens on the third coding base. The other column is a cDNA position that
    /// opens after a 5' UTR.
    fn start_codon_kept_without_five_prime_utr() -> RuleCase {
        let row = |terms: &str, cdna: &str| {
            stop_codon_row("-", terms, cdna, "3-12")
                .at("21:30001-30010")
                .on("ENST90000000002")
        };
        let (lost, kept) = (
            "frameshift_variant,start_lost",
            "frameshift_variant,start_retained_variant",
        );
        RuleCase {
            bucket: "start_codon_kept_without_five_prime_utr",
            reference: row(lost, "3-12"),
            candidate: row(kept, "3-12"),
            further_term: row(
                "frameshift_variant,start_retained_variant,stop_lost",
                "3-12",
            ),
            other_column: (row(lost, "104-113"), row(kept, "104-113")),
        }
    }

    /// The pair of `start_codon_kept_without_five_prime_utr` in its second form, an
    /// insertion of whole codons on which the reference's row holds neither start term, is
    /// set aside; the same rows on a span that opens on the fourth coding base, and a
    /// scored row that keeps start_lost beside the retained start, stay differences.
    #[test]
    fn test_a_kept_start_is_set_aside_in_its_second_form_and_on_the_start_codon_alone() {
        let bucket = "start_codon_kept_without_five_prime_utr";
        let inserted = |terms: &str, span: &str| {
            stop_codon_row("GCA", terms, span, span)
                .at("21:30000-30001")
                .on("ENST90000000002")
                .to_string()
        };
        let kept = "inframe_insertion,start_retained_variant";
        assert_eq!(
            masked(
                bucket,
                &inserted("inframe_insertion", "2-3"),
                &inserted(kept, "2-3")
            ),
            (1, Vec::new(), Vec::new())
        );
        let (count, perl_open, rust_open) = masked(
            bucket,
            &inserted("inframe_insertion", "4-5"),
            &inserted(kept, "4-5"),
        );
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
        let case = start_codon_kept_without_five_prime_utr();
        let both_start_terms = case.candidate.to_string().replace(
            "start_retained_variant",
            "start_lost,start_retained_variant",
        );
        let (count, perl_open, rust_open) =
            masked(bucket, &case.reference.to_string(), &both_start_terms);
        assert_eq!((count, perl_open.len(), rust_open.len()), (0, 1, 1));
    }

    row_pair_rule_checks!(
        stop_codon_overlap_genomic_span_rule_tests,
        super::stop_codon_overlap_genomic_span()
    );
    row_pair_rule_checks!(
        stop_retained_without_stop_codon_test_rule_tests,
        super::stop_retained_without_stop_codon_test()
    );
    row_pair_rule_checks!(
        stop_retained_three_base_stop_codon_test_rule_tests,
        super::stop_retained_three_base_stop_codon_test()
    );
    row_pair_rule_checks!(
        stop_codon_removed_without_stop_codon_test_rule_tests,
        super::stop_codon_removed_without_stop_codon_test()
    );
    row_pair_rule_checks!(
        stop_codon_overlap_insertion_range_rule_tests,
        super::stop_codon_overlap_insertion_range()
    );
    row_pair_rule_checks!(
        stop_codon_settled_before_unknown_bases_rule_tests,
        super::stop_codon_settled_before_unknown_bases()
    );
    row_pair_rule_checks!(
        frameshift_inside_a_held_stop_codon_rule_tests,
        super::frameshift_inside_a_held_stop_codon()
    );
    row_pair_rule_checks!(
        splice_term_past_the_record_span_rule_tests,
        super::splice_term_past_the_record_span()
    );
    row_pair_rule_checks!(
        pachysolen_nuclear_code_read_as_standard_rule_tests,
        super::pachysolen_nuclear_code_read_as_standard()
    );
    row_pair_rule_checks!(
        start_codon_removed_without_five_prime_utr_rule_tests,
        super::start_codon_removed_without_five_prime_utr()
    );
    row_pair_rule_checks!(
        start_codon_kept_without_five_prime_utr_rule_tests,
        super::start_codon_kept_without_five_prime_utr()
    );

    /// A row that ends at its ninth, tenth or eleventh column is read without its Codons
    /// and without its Amino_acids. A rule that reads positions takes its pair as on whole
    /// rows, and each of the three rules that read the Codons, and the rule that reads the
    /// Amino_acids alone, takes its pair only where both rows hold the twelfth column.
    #[test]
    fn test_a_row_that_ends_before_its_codons_column_is_read_without_them() {
        let cut = |row: &Row, columns: usize| {
            let line = row.to_string();
            let cells: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
            cells[..columns].join("\t") + "\n"
        };
        let masked_at = |case: &RuleCase, columns: usize| {
            masked(
                case.bucket,
                &cut(&case.reference, columns),
                &cut(&case.candidate, columns),
            )
        };
        let of_positions = [
            stop_codon_overlap_genomic_span(),
            stop_retained_without_stop_codon_test(),
            stop_retained_three_base_stop_codon_test(),
            stop_codon_removed_without_stop_codon_test(),
            stop_codon_overlap_insertion_range(),
            splice_term_past_the_record_span(),
            start_codon_kept_without_five_prime_utr(),
        ];
        let of_codons = [
            stop_codon_settled_before_unknown_bases(),
            frameshift_inside_a_held_stop_codon(),
            pachysolen_nuclear_code_read_as_standard(),
        ];
        for columns in [9, 10, 11] {
            for case in &of_positions {
                assert_eq!(
                    masked_at(case, columns),
                    (1, Vec::new(), Vec::new()),
                    "{} at {columns} columns",
                    case.bucket
                );
            }
            for case in &of_codons {
                let (count, perl_open, rust_open) = masked_at(case, columns);
                assert_eq!(
                    (count, perl_open.len(), rust_open.len()),
                    (0, 1, 1),
                    "{} at {columns} columns",
                    case.bucket
                );
            }
        }
        for case in &of_codons {
            assert_eq!(
                masked_at(case, 12),
                (1, Vec::new(), Vec::new()),
                "{}",
                case.bucket
            );
        }
        let of_amino_acids = start_codon_removed_without_five_prime_utr();
        for columns in [9, 10, 11] {
            let (count, perl_open, rust_open) = masked_at(&of_amino_acids, columns);
            assert_eq!(
                (count, perl_open.len(), rust_open.len()),
                (0, 1, 1),
                "{columns} columns"
            );
        }
        assert_eq!(
            masked_at(&of_amino_acids, 12),
            (1, Vec::new(), Vec::new()),
            "12 columns"
        );
    }

    /// The transcripts of the made-up cache `tests/fixtures/coding_ends/cache` of the tuple
    /// comparison, on sequence `7` from base 1001: `FLAGGED` shows no stop codon and carries
    /// `cds_end_NF`, `CLOSED` shows a start codon and a stop codon, and `START_OPEN` shows
    /// no start codon.
    fn made_up_cache() -> CacheTranscripts {
        let cache = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../vep-outputs-compare/tests/fixtures/coding_ends/cache"
        );
        CacheTranscripts::new(Path::new(cache), 2)
    }

    /// One row of a structural allele on `feature` of the made-up cache.
    fn structural_row(location: &str, feature: &str, consequence: &str) -> String {
        format!(
            "sv1\t{location}\tdeletion\tENSG1\t{feature}\tTranscript\t{consequence}\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH\n"
        )
    }

    /// The pairs the registry's rules take of one reference row and one vep-rs row on
    /// `feature`, by bucket, and how many tuples the mask sets aside, in a comparison that
    /// was given `cache`.
    fn taken_beside(
        cache: Option<&CacheTranscripts>,
        location: &str,
        feature: &str,
        terms: (&str, &str),
    ) -> (Vec<(String, u64)>, usize) {
        let outputs = row_outputs(
            &structural_row(location, feature, terms.0),
            &structural_row(location, feature, terms.1),
        );
        let (excluded_rust, excluded_perl, buckets) = filter_registry_swap_pairs(
            &outputs.owned.sides(),
            &row_pair_rules(),
            std::slice::from_ref(&outputs.perl_path),
            std::slice::from_ref(&outputs.rust_path),
            "vep-rs",
            &no_inputs(),
            CodingEndSources {
                transcripts: cache,
                reference: None,
            },
        )
        .unwrap();
        let taken = buckets
            .into_iter()
            .filter(|(_, pairs)| *pairs > 0)
            .collect();
        (taken, excluded_rust.len() + excluded_perl.len())
    }

    #[test]
    fn test_a_structural_deletion_over_a_coding_end_without_a_stop_codon_is_set_aside_beside_the_cache(
    ) {
        let cache = made_up_cache();
        let lost = "stop_lost,coding_sequence_variant,3_prime_UTR_variant,feature_truncation";
        let kept = "coding_sequence_variant,3_prime_UTR_variant,feature_truncation";
        let class = "structural_deletion_stop_lost_without_stop_codon".to_string();
        let over_the_end = "7:1010-1040";
        for feature in ["FLAGGED", "FLAGGED_0", "OPEN_WHOLE"] {
            assert_eq!(
                taken_beside(Some(&cache), over_the_end, feature, (lost, kept)),
                (vec![(class.clone(), 1)], 2),
                "{feature}"
            );
            // A comparison that was given no cache sets nothing aside.
            assert_eq!(
                taken_beside(None, over_the_end, feature, (lost, kept)),
                (Vec::new(), 0),
                "{feature}"
            );
        }
        // A stop codon that stands, a record that is not read, a transcript the cache does
        // not hold and a Location no record of the transcript holds.
        for (location, feature) in [
            (over_the_end, "CLOSED"),
            (over_the_end, "UNREAD"),
            (over_the_end, "ENST_ABSENT"),
            ("7:3010-3040", "FLAGGED"),
            ("8:1010-1040", "FLAGGED"),
        ] {
            assert_eq!(
                taken_beside(Some(&cache), location, feature, (lost, kept)),
                (Vec::new(), 0),
                "{feature} at {location}"
            );
        }
        // The masks of a release hand the cache they were given to the pair rules.
        let outputs = row_outputs(
            &structural_row(over_the_end, "FLAGGED", lost),
            &structural_row(over_the_end, "FLAGGED", kept),
        );
        let masked = |release: &str, coding_ends: Option<&CacheTranscripts>| {
            let masks = apply_masks(
                &outputs.owned.sides(),
                &MaskContext {
                    reference_release: release,
                    cache: &CacheAuthority::default(),
                    scored_engine: "vep-rs",
                    input_variants: Some(&[]),
                    max_sv_size: MAX_SV_SIZE,
                    perl_paths: std::slice::from_ref(&outputs.perl_path),
                    rust_paths: std::slice::from_ref(&outputs.rust_path),
                    input_by_output: &no_inputs(),
                    coding_ends,
                },
            )
            .unwrap();
            (
                bucket(&masks.swap_by_bucket, &class),
                masks.swap_rust.len() + masks.swap_perl.len(),
            )
        };
        assert_eq!(masked("116.2", Some(&cache)), (1, 2));
        assert_eq!(masked("116.3", Some(&cache)), (1, 2));
        assert_eq!(masked("116.2", None), (0, 0));
        // A row of vep-rs that drops a further term, or keeps a stop term in the lost
        // one's place.
        for wrong in [
            "coding_sequence_variant,3_prime_UTR_variant",
            "stop_retained_variant,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
        ] {
            assert_eq!(
                taken_beside(Some(&cache), over_the_end, "FLAGGED", (lost, wrong)),
                (Vec::new(), 0),
                "{wrong}"
            );
        }
    }

    #[test]
    fn test_the_start_co_emission_pair_beside_a_cache_is_set_aside_only_where_a_start_codon_stands()
    {
        let cache = made_up_cache();
        let both = "5_prime_UTR_variant,coding_sequence_variant,feature_truncation,start_lost,start_retained_variant";
        let one = "5_prime_UTR_variant,coding_sequence_variant,feature_truncation,start_lost";
        let class = "start_cooccurrence_swap".to_string();
        let over_the_start = "7:990-1010";
        let set_aside = (vec![(class, 1)], 2);
        assert_eq!(
            taken_beside(Some(&cache), over_the_start, "CLOSED", (both, one)),
            set_aside
        );
        // No start codon is shown, with the flag and without; the record is not read; the
        // cache holds no record of the transcript.
        for feature in ["START_OPEN", "START_FLAGGED", "UNREAD", "ENST_ABSENT"] {
            assert_eq!(
                taken_beside(Some(&cache), over_the_start, feature, (both, one)),
                (Vec::new(), 0),
                "{feature}"
            );
            // Without a cache the class stands as it stood.
            assert_eq!(
                taken_beside(None, over_the_start, feature, (both, one)),
                set_aside,
                "{feature}"
            );
        }
        // A record of the transcript that does not hold the pair's Location, on its own
        // sequence or on another, is no record of the pair.
        for elsewhere in ["7:2990-3010", "8:990-1010"] {
            assert_eq!(
                taken_beside(Some(&cache), elsewhere, "CLOSED", (both, one)),
                (Vec::new(), 0),
                "{elsewhere}"
            );
        }
        // The row vep-rs owes where no start codon is shown is the start codon overlap's.
        let neither = "5_prime_UTR_variant,coding_sequence_variant,feature_truncation";
        assert_eq!(
            taken_beside(Some(&cache), over_the_start, "START_OPEN", (both, neither)),
            (
                vec![("start_codon_overlap_without_start_codon".to_string(), 1)],
                2
            )
        );
    }

    const SPLICE_BUCKET: &str = "splice_term_past_the_record_span";
    /// The terms of the pair two positions hold: the reference's and the scored engine's.
    const SPLICE_TERMS: [&str; 2] = [
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant",
        "frameshift_variant,splice_region_variant",
    ];
    const SPLICE_LOCATION: &str = "7:100999-101000";

    #[test]
    fn test_the_pair_is_set_aside_through_the_registry_s_rule_alone() {
        let rules = row_pair_rules();
        assert!(rules.iter().any(|rule| rule.bucket == SPLICE_BUCKET));
        let others: Vec<&'static ExclusionRule> = rules
            .into_iter()
            .filter(|rule| rule.bucket != SPLICE_BUCKET)
            .collect();
        let case = splice_term_past_the_record_span();
        let outputs = row_outputs(&case.reference.to_string(), &case.candidate.to_string());
        let (swap_rust, swap_perl, buckets) = outputs.registry_mask(&others, &no_inputs()).unwrap();
        assert!(!names(&buckets).contains(&SPLICE_BUCKET));
        assert!(swap_perl.is_empty() && swap_rust.is_empty());
    }

    /// Two bases at an exon's end replaced by five: the first column gives the row its
    /// record's two reference bases, and the same two rows under an identifier have an
    /// insertion's columns.
    #[test]
    fn test_a_row_whose_first_column_spells_the_alleles_is_set_aside_at_two_positions() {
        for (name, count) in [("7_100999_AC/GTTCA", 1), ("1993671", 0)] {
            let [reference, candidate] = SPLICE_TERMS.map(|terms| {
                stop_codon_row("GTTCA", terms, "999-1000", "899-900")
                    .at(SPLICE_LOCATION)
                    .named(name)
                    .to_string()
            });
            let outputs = row_outputs(&reference, &candidate);
            let (swap_rust, swap_perl, buckets) = outputs
                .registry_mask(&row_pair_rules(), &no_inputs())
                .unwrap();
            assert_eq!(bucket(&buckets, SPLICE_BUCKET), count, "{name}");
            if count == 1 {
                assert_eq!(
                    (&swap_perl, &swap_rust),
                    (&outputs.owned.perl, &outputs.owned.rust),
                    "{name}"
                );
            } else {
                assert!(swap_perl.is_empty() && swap_rust.is_empty(), "{name}");
            }
        }
    }

    const VCF_HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// Two bases replaced by five, and five bases inserted after the first of the two: both
    /// records give their rows the same Location and the same allele.
    const REPLACED: &str = "7\t100999\t1993671\tAC\tGTTCA\t.\tPASS\t.\n";
    const INSERTED: &str = "7\t100999\t1993671\tA\tAGTTCA\t.\tPASS\t.\n";

    /// The reference's row and the scored engine's row of the pair on one transcript.
    fn hand_over_rows(feature: &str) -> [String; 2] {
        SPLICE_TERMS.map(|terms| {
            stop_codon_row("GTTCA", terms, "999-1000", "899-900")
                .at(SPLICE_LOCATION)
                .on(feature)
                .to_string()
        })
    }

    fn input_file(folder: &TempDir, record: &str) -> PathBuf {
        let path = folder.path.join("in.vcf");
        std::fs::write(&path, format!("{VCF_HEAD}{record}")).unwrap();
        path
    }

    #[test]
    fn test_a_row_under_an_identifier_is_set_aside_where_its_input_is_handed() {
        let [reference, candidate] = hand_over_rows("ENST00000635895");
        let outputs = row_outputs(&reference, &candidate);
        let inputs = outputs.inputs(&input_file(&outputs.folder, REPLACED));
        let (swap_rust, swap_perl, buckets) =
            outputs.registry_mask(&row_pair_rules(), &inputs).unwrap();
        assert_eq!(
            (bucket(&buckets, SPLICE_BUCKET), &swap_perl, &swap_rust),
            (1, &outputs.owned.perl, &outputs.owned.rust)
        );
        let masks = masks_of(
            "116.2",
            &outputs.owned,
            std::slice::from_ref(&outputs.perl_path),
            std::slice::from_ref(&outputs.rust_path),
            &inputs,
        );
        assert_eq!(
            (
                bucket(&masks.swap_by_bucket, SPLICE_BUCKET),
                &masks.swap_perl,
                &masks.swap_rust
            ),
            (1, &outputs.owned.perl, &outputs.owned.rust)
        );
    }

    /// The first file's record replaces two bases and the second file's inserts at the same
    /// Location. Each file's outputs carry the pair on a transcript of its own, and both
    /// carry it on a third, whose row then has two readings.
    #[test]
    fn test_each_row_reads_the_input_of_its_own_output_in_either_order_of_the_files() {
        let own = ["ENST00000000001", "ENST00000000002"];
        let files: Vec<(Outputs, PathBuf)> = [REPLACED, INSERTED]
            .into_iter()
            .zip(own)
            .map(|(record, feature)| {
                let ([reference, candidate], [shared_reference, shared_candidate]) =
                    (hand_over_rows(feature), hand_over_rows("ENST00000000003"));
                let outputs = row_outputs(
                    &format!("{reference}{shared_reference}"),
                    &format!("{candidate}{shared_candidate}"),
                );
                let input = input_file(&outputs.folder, record);
                (outputs, input)
            })
            .collect();
        let both = |side: fn(&Outputs) -> &Vec<Texts>| -> Vec<Texts> {
            files
                .iter()
                .flat_map(|(outputs, _)| side(outputs).iter().cloned())
                .collect()
        };
        let pooled = sides_of(
            &both(|outputs| &outputs.perl),
            &both(|outputs| &outputs.rust),
        );
        let inputs: BTreeMap<PathBuf, PathBuf> = files
            .iter()
            .flat_map(|(outputs, input)| outputs.inputs(input))
            .collect();
        for order in [[0, 1], [1, 0]] {
            let paths = |side: fn(&Outputs) -> &PathBuf| -> Vec<PathBuf> {
                order.iter().map(|at| side(&files[*at].0).clone()).collect()
            };
            let masks = masks_of(
                "116.2",
                &pooled,
                &paths(|outputs| &outputs.perl_path),
                &paths(|outputs| &outputs.rust_path),
                &inputs,
            );
            assert_eq!(bucket(&masks.swap_by_bucket, SPLICE_BUCKET), 1);
            let of_the_first = BTreeSet::from([own[0].to_string()]);
            assert_eq!(
                features(&pooled, &masks.swap_perl),
                of_the_first,
                "the pair of the file whose record replaces"
            );
            assert_eq!(features(&pooled, &masks.swap_rust), of_the_first);
        }
        for (at, (outputs, _)) in files.iter().enumerate() {
            let alone = masks_of(
                "116.2",
                &outputs.owned,
                std::slice::from_ref(&outputs.perl_path),
                std::slice::from_ref(&outputs.rust_path),
                &inputs,
            );
            assert_eq!(
                bucket(&alone.swap_by_bucket, SPLICE_BUCKET),
                if at == 0 { 2 } else { 0 },
                "one file reads its own input alone"
            );
        }
    }

    /// Several records of one input at a row's Location and allele: records that give the
    /// same bases are one reading; where they differ the record whose identifier is the row's
    /// first column decides, and the row takes none where no record carries that identifier
    /// or two that do differ.
    #[test]
    fn test_a_row_several_records_read_differently_takes_the_reading_of_its_own_record() {
        let folder = TempDir::new("vep_sv_compare_pairs_").unwrap();
        let records = [
            "7\t100998\tmade_up_0\tACG\tGTTCA",
            "7\t100998\tmade_up_1\tACG\tTTT,GTTCA",
            "7\t500999\tmade_up_2\tAC\tGTT",
            "7\t500999\tmade_up_3\tA\tAGTT",
            "7\t700998\tmade_up_4\tACG\tGTTCA",
            "7\t900999\tmade_up_5\tA\tAA",
            "7\t900999\tmade_up_6\tAC\tA",
            "7\t950999\tmade_up_7\tAC\tGTT",
            "7\t950999\tmade_up_7\tA\tAGTT",
        ]
        .map(|record| format!("{record}\t.\tPASS\t.\n"));
        let path = input_file(&folder, &records.concat());
        let row = |name: &str, location: &str, allele: &str| {
            pair_row(&[
                name,
                location,
                allele,
                "ENSG1",
                "ENST1",
                "Transcript",
                "frameshift_variant",
                "-",
                "-",
            ])
        };
        let mut rows = [
            row("1993671", "7:100998-101000", "GTTCA"),
            row("1993671", "7:500999-501000", "GTT"),
            row("made_up_3", "7:500999-501000", "GTT"),
            row("made_up_2", "7:500999-501000", "GTT"),
            row("1993671", "7:700998-701000", "GTTCA"),
            row("made_up_5", "7:900999-901000", "A"),
            row("made_up_7", "7:950999-951000", "GTT"),
        ];
        read_record_references(&path, &mut rows.iter_mut().collect::<Vec<_>>()).unwrap();
        let read: Vec<Option<&str>> = rows.iter().map(|row| row.reference.as_deref()).collect();
        assert_eq!(
            read,
            [
                Some("ACG"),
                None,
                Some(""),
                Some("AC"),
                Some("ACG"),
                Some(""),
                None
            ],
            "alike; differently under neither identifier, under the inserting record's, \
             under the replacing record's; one record; an insertion beside a deletion whose \
             whole alleles fit its row; two records of the row's own identifier that differ"
        );
    }

    #[test]
    fn test_the_input_is_opened_only_where_a_row_asks() {
        let masked = |perl_rows: &str,
                      rust_rows: &str,
                      rules: &[&'static ExclusionRule]|
         -> Result<NamedCounts, Failure> {
            let outputs = row_outputs(perl_rows, rust_rows);
            let missing = outputs.folder.path.join("absent.vcf");
            outputs
                .registry_mask(rules, &outputs.inputs(&missing))
                .map(|(_, _, buckets)| buckets)
        };
        let rules = row_pair_rules();
        let [reference, candidate] = hand_over_rows("ENST00000635895");
        let asked = masked(&reference, &candidate, &rules).unwrap_err();
        assert_eq!(
            asked.kind,
            FailureKind::Io,
            "a row under an identifier asks"
        );
        assert!(
            !names(&masked(&reference, &candidate, &[]).unwrap()).contains(&SPLICE_BUCKET),
            "no row-aware rule"
        );
        assert_eq!(
            bucket(
                &masked(&reference, &reference, &rules).unwrap(),
                SPLICE_BUCKET
            ),
            0,
            "no unmatched pair"
        );
        let covered = [
            stop_codon_row("-", "intron_variant", "-", "-").to_string(),
            stop_codon_row("-", "intron_variant,splice_region_variant", "-", "-").to_string(),
        ];
        assert_eq!(
            bucket(
                &masked(&covered[0], &covered[1], &rules).unwrap(),
                "covered_splice_region_swap"
            ),
            1,
            "a pair an earlier rule takes"
        );
        let spelled =
            [&reference, &candidate].map(|row| row.replace("1993671", "7_100999_AC/GTTCA"));
        assert_eq!(
            bucket(
                &masked(&spelled[0], &spelled[1], &rules).unwrap(),
                SPLICE_BUCKET
            ),
            1,
            "a first column that spells the alleles"
        );
    }
}

#[cfg(test)]
mod pair_rows_tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::tempdir::TempDir;

    /// A scratch folder of outputs and the table of the tuples asked for.
    struct Folder {
        folder: TempDir,
        strings: Strings,
    }

    const DELETED: [&str; 5] = [
        "21:5-9",
        "deletion",
        "ENST1",
        "Transcript",
        "feature_truncation",
    ];
    const FIRST: [&str; 5] = ["21:5-9", "deletion", "ENST1", "Transcript", "a,b"];
    const SECOND: [&str; 5] = ["MT:7", "N[22:5[", "ENST2", "Transcript", "c"];

    impl Folder {
        fn new() -> Folder {
            Folder {
                folder: TempDir::new("vep_sv_compare_pair_rows_").unwrap(),
                strings: Strings::from_texts([DELETED, FIRST, SECOND].into_iter().flatten()),
            }
        }

        fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
            let path = self.folder.path.join(name);
            std::fs::write(&path, bytes).unwrap();
            path
        }

        fn output(&self, name: &str, rows: &[String]) -> PathBuf {
            let text = format!("#Uploaded_variation\tLocation\n{}", rows.concat());
            self.write(name, text.as_bytes())
        }

        fn tuple(&self, texts: [&str; 5]) -> Tuple {
            self.strings.tuple(texts).unwrap()
        }

        /// What each wanted tuple's rows read, by the tuple's Location: the Location and the
        /// Allele as written and the cDNA position.
        fn rows(
            &self,
            paths: &[PathBuf],
            wanted: &[[&str; 5]],
            engine: &str,
        ) -> Vec<(String, Vec<[String; 3]>)> {
            let wanted: TupleSet = wanted.iter().map(|texts| self.tuple(*texts)).collect();
            let found = pair_rows_by_tuple(&self.strings, paths, &wanted, engine, &BTreeMap::new())
                .unwrap();
            let mut read: Vec<(String, Vec<[String; 3]>)> = found
                .into_iter()
                .map(|(tuple, rows)| {
                    (
                        self.strings.text(tuple.location()).to_string(),
                        rows.into_iter()
                            .map(|row| [row.location, row.allele, row.cdna_position])
                            .collect(),
                    )
                })
                .collect();
            read.sort();
            read
        }
    }

    fn row(
        name: &str,
        location: &str,
        allele: &str,
        feature: &str,
        consequence: &str,
        cdna: &str,
        columns: usize,
    ) -> String {
        let all = [
            name,
            location,
            allele,
            "ENSG1",
            feature,
            "Transcript",
            consequence,
            cdna,
            "4-6",
        ];
        all[..columns].join("\t") + "\n"
    }

    fn read(location: &str, rows: &[[&str; 3]]) -> (String, Vec<[String; 3]>) {
        (
            location.to_string(),
            rows.iter().map(|row| row.map(str::to_string)).collect(),
        )
    }

    #[test]
    fn a_row_of_the_scored_engine_is_keyed_under_that_engine_s_allele() {
        let folder = Folder::new();
        let scored = folder.output(
            "scored.txt",
            &[row(
                "sv1",
                "21:5-9",
                "<DEL>",
                "ENST1",
                "feature_truncation",
                "10-12",
                9,
            )],
        );
        let paths = [scored];
        let mapped = [read("21:5-9", &[["21:5-9", "<DEL>", "10-12"]])];
        assert_eq!(folder.rows(&paths, &[DELETED], "fastvep"), mapped);
        assert_eq!(folder.rows(&paths, &[DELETED], " FastVEP "), mapped);
        assert_eq!(folder.rows(&paths, &[DELETED], "vep-rs"), []);
    }

    #[test]
    fn the_rows_of_all_outputs_are_pooled_in_their_order_and_a_short_row_drops_its_tuple() {
        let folder = Folder::new();
        let a = folder.output(
            "a.txt",
            &[
                row("sv1", "21:5-9", "deletion", "ENST1", "b, a,b", "10-12", 9),
                row("sv2", "chrM:07", "N[chr22:5[", "ENST2", "c", "1-2", 9),
                row("sv2", "chrM:07", "N[chr22:5[", "ENST2", "c", "1-2", 9),
                "\n".to_string(),
                "short\trow\n".to_string(),
            ],
        );
        let b = folder.output(
            "b.txt",
            &[
                row("sv1", "21:5-9", "deletion", "ENST1", "a,b", "10-12", 8),
                row("sv2", "M:7", "N[22:5[", "ENST2", "c", "3-4", 9),
                row("sv2", "chrM:07", "N[chr22:5[", "ENST2", "c", "1-2", 9),
            ],
        );
        let absent = folder.folder.path.join("absent.txt");
        let wanted = [FIRST, SECOND];
        let as_written = ["chrM:07", "N[chr22:5[", "1-2"];
        let bare = ["M:7", "N[22:5[", "3-4"];
        assert_eq!(
            folder.rows(&[a.clone(), b.clone()], &wanted, "vep-rs"),
            [read("MT:7", &[as_written, bare])]
        );
        assert_eq!(
            folder.rows(&[b, a.clone()], &wanted, "vep-rs"),
            [read("MT:7", &[bare, as_written])]
        );
        let alone = [
            read("21:5-9", &[["21:5-9", "deletion", "10-12"]]),
            read("MT:7", &[as_written]),
        ];
        assert_eq!(
            folder.rows(std::slice::from_ref(&a), &wanted, "vep-rs"),
            alone
        );
        assert_eq!(folder.rows(&[absent, a], &wanted, "vep-rs"), alone);
    }

    #[test]
    fn a_line_ends_at_a_line_feed_a_carriage_return_or_both() {
        let folder = Folder::new();
        let dos =
            row("sv1", "21:5-9", "deletion", "ENST1", "a,b", "10-12", 9).replace('\n', "\r\n");
        let mac = row("sv1", "21:5-9", "deletion", "ENST1", "a,b", "7", 9).replace('\n', "\r");
        let output = folder.write("c.txt", format!("#h\r\n{dos}{mac}").as_bytes());
        assert_eq!(
            folder.rows(&[output], &[FIRST], "vep-rs"),
            [read(
                "21:5-9",
                &[["21:5-9", "deletion", "10-12"], ["21:5-9", "deletion", "7"]]
            )]
        );
    }

    #[test]
    fn no_output_is_read_for_no_tuple_and_an_output_that_is_not_utf_8_ends_the_comparison() {
        let folder = Folder::new();
        let broken = folder.write("broken.txt", b"#h\nsv1\t21:5-9\t\xff\n");
        let paths = [broken];
        assert_eq!(folder.rows(&paths, &[], "vep-rs"), []);
        let wanted: TupleSet = [folder.tuple(FIRST)].into_iter().collect();
        let refused =
            pair_rows_by_tuple(&folder.strings, &paths, &wanted, "vep-rs", &BTreeMap::new())
                .unwrap_err();
        assert_eq!(refused.kind, FailureKind::Decode);
        assert!(
            refused.message.ends_with("line 2 is not UTF-8"),
            "{refused}"
        );
    }
}
