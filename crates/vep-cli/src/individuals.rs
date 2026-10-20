// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-individual reading of a VCF record (`--individual`) and the
//! reference-allele row of `--process_ref_homs`.
//!
//! Ensembl VEP builds one variation feature per record, then copies it once per
//! named individual with that individual's genotype (`Parser/VCF.pm:355-376` and
//! `create_individual_VariationFeatures`, `:595-665`, at release/116.2): the
//! genotype is read by ensembl-io (`BaseVCF4.pm:928-955` at release/116), each
//! called allele mapped from the record's written alleles to the chopped ones,
//! the non-reference alleles the individual carries become the copy's allele
//! string, and an individual carrying none is a reference-homozygous copy, the
//! reference alone unless `--process_ref_homs` makes it `REF/REF`. A copy is a
//! feature of its own from there on: validated, annotated, picked and written
//! apart from the record's other copies. Here a copy is the record's alleles the
//! individual carries, each marked with the individual's name and genotype, laid
//! out one copy after another in the individuals' order.
//!
//! `--process_ref_homs` adds the reference allele's row to every feature of every
//! sequence record (`OutputFactory.pm:509-516`), first among the feature's
//! alleles and numbered 0 (`VariationFeatureOverlap.pm:469-474`,
//! `BaseVariationFeatureOverlap.pm:363-373`); a structural record is read through
//! another method and gains none (`:534-560`).

use std::collections::HashMap;
use std::sync::Arc;

use vep_core::sv_routing::perl_split;
use vep_core::variant::{keeps_uploaded_spelling, InputVariant};
use vep_core::variation::perl_numeric;

use crate::vcf_parser::{SampleColumns, VcfColumns};

/// The individuals `--individual` names: every sample of the header, or a list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IndividualRequest {
    /// `all`: every sample column, in header order (`Parser/VCF.pm:609`).
    All,
    /// The named samples, in the order given; every name must be a header sample
    /// (`:612-617`).
    Named(Vec<String>),
}

impl IndividualRequest {
    /// Ensembl VEP reads `all` from the list's first element alone, in any case
    /// (`Parser/VCF.pm:609`): `--individual all,S1` is `all`.
    pub fn from_list(names: &[String]) -> Self {
        match names.first() {
            Some(first) if first.eq_ignore_ascii_case("all") => IndividualRequest::All,
            _ => IndividualRequest::Named(names.to_vec()),
        }
    }

    /// The individuals to read for a header's sample `columns`: every sample, or the
    /// named list in its order, a repeated name included (`Parser/VCF.pm:622`
    /// copies the record once per listed name). A named list with a name the
    /// header lacks is Ensembl VEP's `die` (`:612-617`), both lists in the message.
    /// ensembl-io keys the samples by name, so a name the header repeats is one
    /// sample reading its last column, and `all` lists the names by that column
    /// (BaseVCF4.pm:785-808).
    pub fn resolve(&self, columns: &SampleColumns) -> Result<Vec<(Arc<str>, usize)>, String> {
        let mut last_column: HashMap<&str, usize> = HashMap::new();
        for (column, name) in columns.names.iter().enumerate() {
            last_column.insert(name.as_str(), column);
        }
        match self {
            IndividualRequest::All => {
                let mut by_column: Vec<(Arc<str>, usize)> = last_column
                    .iter()
                    .map(|(name, column)| (Arc::from(*name), *column))
                    .collect();
                by_column.sort_by_key(|(_, column)| *column);
                Ok(by_column)
            }
            IndividualRequest::Named(names) => {
                let mut individuals = Vec::with_capacity(names.len());
                for name in names {
                    let Some(column) = last_column.get(name.as_str()).copied() else {
                        return Err(format!(
                            "Sample IDs given ({}) do not match samples from VCF ({})",
                            names.join(","),
                            columns.names.join(",")
                        ));
                    };
                    individuals.push((Arc::from(name.as_str()), column));
                }
                Ok(individuals)
            }
        }
    }
}

/// How a record's lines are split per individual: the individuals to read, as
/// [`IndividualRequest::resolve`] resolved them against the input's header, with
/// their sample column indices, and whether the records carry genotypes at all.
#[derive(Debug, Clone)]
pub struct IndividualSplit {
    /// The individuals in their output order, each with its 0-based sample
    /// column (the column after FORMAT is 0).
    individuals: Vec<(Arc<str>, usize)>,
    /// Whether the header's ninth column is `FORMAT`; without it no record yields a
    /// genotype and no copy is written (BaseVCF4.pm:98, :614-618).
    genotyped: bool,
}

impl IndividualSplit {
    /// The split for `request` over a header's sample `columns`; `Err` is Ensembl
    /// VEP's message for a list naming an unknown sample.
    pub fn new(request: &IndividualRequest, columns: &SampleColumns) -> Result<Self, String> {
        Ok(Self {
            individuals: request.resolve(columns)?,
            genotyped: columns.genotyped,
        })
    }

    /// The number of copies a record yields at most, one per individual listed.
    pub fn len(&self) -> usize {
        self.individuals.len()
    }

    /// Whether no individual is listed.
    pub fn is_empty(&self) -> bool {
        self.individuals.is_empty()
    }
}

/// How a run reads each sequence record beyond its alleles: per individual, or
/// with the reference allele's row.
#[derive(Debug, Clone, Copy, Default)]
pub struct RecordReading<'a> {
    pub individuals: Option<&'a IndividualSplit>,
    pub process_ref_homs: bool,
}

impl RecordReading<'_> {
    /// The alleles the record writes: its individuals' copies, or the record
    /// itself with the reference row when `--process_ref_homs` asks for one.
    pub fn finish_record(
        &self,
        mut record: Vec<InputVariant>,
        alleles: &RecordAlleles<'_>,
        columns: &VcfColumns<'_>,
    ) -> Vec<InputVariant> {
        match self.individuals {
            Some(split) => {
                split_by_individual(record, alleles, columns, split, self.process_ref_homs)
            }
            None => {
                if self.process_ref_homs {
                    add_reference_allele_row(&mut record);
                }
                record
            }
        }
    }
}

/// The alleles of one record as the per-individual reading needs them: every
/// allele as written (the reference first, then each ALT), and each as the
/// chopped allele string Ensembl VEP maps it to (`Parser/VCF.pm:357-358`), in
/// the written case, since the copy's alleles are compared before `validate_vf`
/// upper-cases them (`Parser.pm:554`).
pub struct RecordAlleles<'a> {
    pub written: Vec<&'a str>,
    pub chopped: Vec<String>,
}

/// The GT index of a record's FORMAT, or `None` when the record carries no FORMAT
/// or FORMAT has no `GT` key: ensembl-io then reads no genotype from the record
/// (`BaseVCF4.pm:884-920`).
fn genotype_index(columns: &VcfColumns<'_>) -> Option<usize> {
    columns.format?.split(':').position(|key| key == "GT")
}

/// The GT field of the sample column `sample_column` among the record's
/// `samples`, or `None` when the column is absent or stops before the GT key.
fn genotype_field<'a>(
    samples: &[&'a str],
    gt_index: usize,
    sample_column: usize,
) -> Option<&'a str> {
    samples.get(sample_column)?.split(':').nth(gt_index)
}

/// Perl's reading of a GT bit as an index into the record's `allele_count`
/// alleles (`$alleles[$_]`, `BaseVCF4.pm:947`): the bit numified (leading blanks,
/// a sign, digits, a fraction and an exponent, so `+1` and ` 1` are 1, `1e1` is
/// 10 and a digit-led value past the unsigned word, `1e999` or a run of digits,
/// reads as -1, while a bit spelled as infinity or not-a-number is 0, the integer
/// slot Perl's conversion leaves untouched), truncated, a negative index
/// counting back from the last allele; an index past either end names no allele
/// (`None`, an undefined bit).
fn allele_index(bit: &str, allele_count: usize) -> Option<usize> {
    let value = perl_numeric(bit);
    let spelled = bit
        .trim_start()
        .trim_start_matches(['+', '-'])
        .starts_with(|c: char| c.is_ascii_alphabetic());
    let index: i64 = if value.is_nan() || (value.is_infinite() && spelled) {
        0
    } else if value >= 18_446_744_073_709_551_616.0 {
        -1
    } else if value >= 9_223_372_036_854_775_808.0 {
        (value as u64) as i64
    } else {
        value.trunc() as i64
    };
    let count = i64::try_from(allele_count).unwrap_or(i64::MAX);
    let resolved = if index < 0 {
        index.checked_add(count)?
    } else {
        index
    };
    (0..count).contains(&resolved).then_some(resolved as usize)
}

/// The alleles a GT field calls, as indices into the record's `allele_count`
/// written alleles, or `None` for a missing genotype (`BaseVCF4.pm:944-951`,
/// `Parser/VCF.pm:626-627`): the field is split on `|` when it has one, else on
/// `/`, a trailing empty bit being no bit; `.` bits are dropped; each bit is read
/// as a Perl array index ([`allele_index`]) and the alleles are joined again, so an
/// undefined bit is an empty allele (`None` in the list) and a trailing undefined
/// bit vanishes with the second split; a genotype that reads as the empty string
/// is missing.
pub fn called_alleles(gt: &str, allele_count: usize) -> Option<Vec<Option<usize>>> {
    let separator = if gt.contains('|') { '|' } else { '/' };
    let mut bits: Vec<Option<usize>> = perl_split(gt, &[separator])
        .filter(|bit| *bit != ".")
        .map(|bit| allele_index(bit, allele_count))
        .collect();
    let joined_is_empty = matches!(bits.as_slice(), [] | [None]);
    if joined_is_empty {
        return None;
    }
    while bits.last() == Some(&None) {
        bits.pop();
    }
    Some(bits)
}

/// One individual's copy of a record (`Parser/VCF.pm:627-662`), `None` for an
/// empty `record`: the alleles of `record` the individual carries, in the
/// record's order, marked with the individual's name and genotype, led by the
/// reference allele's row under `--process_ref_homs`; a reference-homozygous
/// individual (`:641-650`) takes the reference alone as a feature with no
/// allele, and an individual whose only alleles write no row (a `<NON_REF>`
/// beside real alleles, or the `.` of a non-variant record) a feature with no
/// allele carrying its allele string, each the reference row alone under
/// `--process_ref_homs` ([`add_reference_allele_row`]).
fn individual_copy(
    record: &[InputVariant],
    alleles: &RecordAlleles<'_>,
    name: &Arc<str>,
    ordinal: u32,
    called: &[Option<usize>],
    process_ref_homs: bool,
) -> Option<Vec<InputVariant>> {
    let base = record.first()?;
    let genotype: Vec<String> = called
        .iter()
        .map(|bit| bit.map_or_else(String::new, |index| alleles.chopped[index].clone()))
        .collect();
    // `grep {$_ ne $ref && $_ !~ /\*/}` over the mapped alleles, each distinct
    // string once (`:635`), here in the record's allele order; an undefined bit
    // is an empty allele, not the reference, which writes no row.
    let reference = &alleles.chopped[0];
    let mut non_ref: Vec<usize> = Vec::new();
    for index in 1..alleles.chopped.len() {
        let allele = &alleles.chopped[index];
        if !called.contains(&Some(index)) || allele == reference || allele.contains('*') {
            continue;
        }
        if !non_ref.iter().any(|&seen| alleles.chopped[seen] == *allele) {
            non_ref.push(index);
        }
    }
    let empty_allele = called.contains(&None);
    // `UPLOADED_ALLELE` (OutputFactory.pm:1360) is Ensembl VEP's
    // `original_allele_string`, the raw string of a record whose alleles differ
    // in length, which every copy inherits, else the copy's own string.
    let inherits_raw_string = base
        .uploaded_allele_string
        .as_deref()
        .is_some_and(|raw| keeps_uploaded_spelling(raw, false));
    let mark = |variant: &mut InputVariant, copy_string: &str, hom_ref: bool| {
        variant.individual = Some(Arc::clone(name));
        variant.copy_ordinal = ordinal;
        variant.genotype.clone_from(&genotype);
        variant.hom_ref = hom_ref;
        if !inherits_raw_string {
            variant.uploaded_allele_string = Some(copy_string.into());
        }
    };
    if non_ref.is_empty() && !empty_allele {
        let mut copy = feature_without_allele(base);
        // The copy is the record's working reference (`:604`, `:642`, `$ref` of
        // the post-processed feature): the minimised one of a bi-allelic indel
        // `post_process_vfs` trimmed past the chop, so `CTA GA`, read as `CT/G`
        // over two bases, is `CT` to a reference-homozygous individual and its
        // reference row fits the span. Ensembl VEP reaches this branch for that
        // record only when a chopped bit equals the minimised reference, which
        // `CTA` and `GA` never do, so its `0/0` individual carries `CTA` as an
        // alternate (`:635-639`); vep-rs compares chopped with chopped, the
        // divergence `docs/intended-divergences.md` lists under the trimmable
        // indel.
        let working_reference = String::from_utf8_lossy(&base.ref_allele).to_ascii_uppercase();
        let copy_string = if process_ref_homs {
            format!("{working_reference}/{working_reference}")
        } else {
            working_reference
        };
        mark(&mut copy, &copy_string, true);
        copy.record_allele_string_multi = Some(copy_string.into());
        return Some(if process_ref_homs {
            vec![copy.reference_allele_row()]
        } else {
            vec![copy]
        });
    }
    let mut copy_string = reference.to_ascii_uppercase();
    for &index in &non_ref {
        copy_string.push('/');
        copy_string.push_str(&alleles.chopped[index].to_ascii_uppercase());
    }
    if empty_allele {
        copy_string.push('/');
    }
    let mut copy: Vec<InputVariant> = non_ref
        .iter()
        .enumerate()
        .filter_map(|(number, &index)| {
            // The record's row for the allele: a repeated ALT has one row, the
            // first ALT of that string (`perl_allele_numbers`).
            let allele = &alleles.chopped[index];
            let mut variant = record
                .iter()
                .find(|v| !v.no_allele && alleles.chopped.get(v.allele_index + 1) == Some(allele))?
                .clone();
            // The copy's alleles are numbered along its own string
            // (`VariationFeatureOverlap.pm:477-490`).
            variant.allele_num = Some(number as u32 + 1);
            Some(variant)
        })
        .collect();
    if copy.is_empty() {
        copy.push(feature_without_allele(base));
    }
    // A copy whose string is not the record's chopped string (a subset of its
    // alleles, or every allele but a `*`) works with its own: the JSON
    // `allele_string`, `VARIANT_CLASS` and the generated name read it. A copy of
    // the whole record keeps the record's, which `post_process_vfs` may have
    // minimised past the chop (`CTA GA` is `CT/G`, the chop kept as
    // `original_allele_string`), so it reads as the record does without
    // `--individual`; Ensembl VEP rebuilds such a copy from the chopped bits
    // against the minimised reference (`CT/CTA/GA`), the divergence
    // `docs/intended-divergences.md` lists under the trimmable indel.
    let record_string = match &base.original_allele_string {
        Some(original) => original.to_ascii_uppercase(),
        None => base.record_allele_string().to_ascii_uppercase(),
    };
    let subset = copy_string != record_string;
    for variant in &mut copy {
        if subset {
            variant.record_allele_string_multi = Some(copy_string.as_str().into());
        }
        mark(variant, &copy_string, false);
    }
    if process_ref_homs {
        add_reference_allele_row(&mut copy);
    }
    Some(copy)
}

/// A feature over `base`'s span with no allele to annotate, carrying the
/// record-level fields of `base`.
fn feature_without_allele(base: &InputVariant) -> InputVariant {
    let mut copy = base.clone();
    copy.alt_alleles = vec![b"*".to_vec()];
    copy.allele_string = String::from_utf8_lossy(&base.ref_allele).into_owned();
    copy.allele_index = 0;
    copy.allele_num = None;
    copy.no_allele = true;
    copy.minimised = false;
    copy.original_allele_string = None;
    copy.transcript_consequences.clear();
    copy.most_severe_consequence = None;
    copy
}

/// The record's alleles per individual, one copy after another in the split's
/// order: the per-individual reading of `Parser/VCF.pm:355-376`. An individual
/// whose genotype is missing has no copy, and a record every individual misses
/// yields nothing, as Ensembl VEP reads the next record in its place (`:371-374`);
/// a header without `FORMAT` yields no genotype and so no copy.
pub fn split_by_individual(
    record: Vec<InputVariant>,
    alleles: &RecordAlleles<'_>,
    columns: &VcfColumns<'_>,
    split: &IndividualSplit,
    process_ref_homs: bool,
) -> Vec<InputVariant> {
    if !split.genotyped {
        return Vec::new();
    }
    let Some(gt_index) = genotype_index(columns) else {
        return Vec::new();
    };
    let samples: Vec<&str> = columns.samples().collect();
    let mut out = Vec::new();
    for (ordinal, (name, column)) in split.individuals.iter().enumerate() {
        let Some(gt) = genotype_field(&samples, gt_index, *column) else {
            continue;
        };
        let Some(called) = called_alleles(gt, alleles.written.len()) else {
            continue;
        };
        if let Some(copy) = individual_copy(
            &record,
            alleles,
            name,
            ordinal as u32,
            &called,
            process_ref_homs,
        ) {
            out.extend(copy);
        }
    }
    out
}

/// `--process_ref_homs` on a record, or on an individual's copy of one: the
/// reference allele's row before its alleles, or in place of a feature with no
/// allele (an allele-less or non-variant record, or a copy whose alleles all
/// write no row, `VariationFeatureOverlap.pm:479-491`), since the VCF line of a
/// record holding such a feature carries no CSQ; a structural record gains none.
pub fn add_reference_allele_row(record: &mut Vec<InputVariant>) {
    let Some(first) = record.first() else {
        return;
    };
    if first.is_structural || first.reference_allele_row {
        return;
    }
    let row = first.reference_allele_row();
    if first.no_allele {
        record.clear();
    }
    record.insert(0, row);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `BaseVCF4.pm:944-951` and `Parser/VCF.pm:626-627`: `|` splits a phased
    /// genotype, `/` an unphased one, `.` bits are dropped before anything else is
    /// read, and each bit is a Perl array index: a leading empty bit numifies to the
    /// reference, a trailing one is no bit, `-1` is the last allele, `+1`, ` 1` and
    /// `1a` are 1, `1e1` is 10, a digit-led value past the unsigned word (`1e999`,
    /// a run of digits) is the last allele, a spelled `inf` or `nan` is the
    /// reference, a negative overflow (`-1e999`) an undefined bit, as is an index
    /// past the alleles, an empty allele
    /// that counts (`3/0`) unless it trails (`0/3`); a genotype reading as the
    /// empty string is missing, while two undefined bits (`3/3`) are a copy with
    /// no allele called.
    #[test]
    fn called_alleles_follow_ensembl_ios_genotype_reading() {
        let some = |bits: &[usize]| Some(bits.iter().map(|&b| Some(b)).collect::<Vec<_>>());
        assert_eq!(called_alleles("0/1", 2), some(&[0, 1]));
        assert_eq!(called_alleles("1|0", 2), some(&[1, 0]));
        assert_eq!(called_alleles("./1", 2), some(&[1]));
        assert_eq!(called_alleles("./.", 2), None);
        assert_eq!(called_alleles("1", 2), some(&[1]));
        assert_eq!(called_alleles(".", 2), None);
        assert_eq!(called_alleles("0/0/0", 2), some(&[0, 0, 0]));
        assert_eq!(called_alleles("1/2", 3), some(&[1, 2]));
        assert_eq!(called_alleles("0/1a", 2), some(&[0, 1]));
        assert_eq!(called_alleles("1/", 2), some(&[1]));
        assert_eq!(called_alleles("/1", 2), some(&[0, 1]));
        assert_eq!(called_alleles("", 2), None);
        assert_eq!(called_alleles("/", 2), None);
        assert_eq!(called_alleles("-1/0", 2), some(&[1, 0]));
        assert_eq!(called_alleles("1/-1", 2), some(&[1, 1]));
        assert_eq!(called_alleles("+1/0", 2), some(&[1, 0]));
        assert_eq!(called_alleles(" 1", 2), some(&[1]));
        assert_eq!(
            called_alleles("0/99999999999999999999999", 2),
            some(&[0, 1])
        );
        assert_eq!(called_alleles("1e1", 2), None);
        assert_eq!(called_alleles("0/3", 2), some(&[0]));
        assert_eq!(called_alleles("3/0", 2), Some(vec![None, Some(0)]));
        assert_eq!(called_alleles("3/3", 2), Some(Vec::new()));
        assert_eq!(called_alleles("-3/0", 2), Some(vec![None, Some(0)]));
        assert_eq!(called_alleles("inf/0", 2), some(&[0, 0]));
        assert_eq!(called_alleles("-inf", 2), some(&[0]));
        assert_eq!(called_alleles("nan/1", 2), some(&[0, 1]));
        assert_eq!(called_alleles("1e999/0", 2), some(&[1, 0]));
        assert_eq!(called_alleles("-1e999/0", 2), Some(vec![None, Some(0)]));
    }

    /// BaseVCF4.pm:785-808: the samples are keyed by name, so a repeated header
    /// name is one sample reading its last column, `all` lists the names by that
    /// column, and a list may repeat a name (`Parser/VCF.pm:622`).
    #[test]
    fn a_repeated_sample_name_reads_its_last_column() {
        let columns = SampleColumns {
            names: vec!["S1".into(), "S3".into(), "S1".into()],
            genotyped: true,
        };
        let names = |split: &IndividualSplit| -> Vec<(String, usize)> {
            split
                .individuals
                .iter()
                .map(|(name, column)| (name.to_string(), *column))
                .collect()
        };
        let all = IndividualSplit::new(&IndividualRequest::All, &columns).unwrap();
        assert_eq!(names(&all), [("S3".to_string(), 1), ("S1".to_string(), 2)]);
        let twice = IndividualRequest::Named(vec!["S1".into(), "S1".into(), "S3".into()]);
        let split = IndividualSplit::new(&twice, &columns).unwrap();
        assert_eq!(
            names(&split),
            [
                ("S1".to_string(), 2),
                ("S1".to_string(), 2),
                ("S3".to_string(), 1)
            ]
        );
        assert_eq!(split.len(), 3);
    }

    /// `Parser/VCF.pm:609`: the list's first element decides `all`, in any case.
    #[test]
    fn all_is_read_from_the_first_name_in_any_case() {
        let all = |names: &[&str]| {
            IndividualRequest::from_list(&names.iter().map(|s| s.to_string()).collect::<Vec<_>>())
        };
        assert_eq!(all(&["all"]), IndividualRequest::All);
        assert_eq!(all(&["ALL", "S1"]), IndividualRequest::All);
        assert_eq!(
            all(&["S1", "all"]),
            IndividualRequest::Named(vec!["S1".into(), "all".into()])
        );
    }
}
