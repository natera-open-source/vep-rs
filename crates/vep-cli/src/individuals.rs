// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-individual readings of a VCF record (`--individual` and
//! `--individual_zyg`) and the reference-allele row of `--process_ref_homs`.
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
//! Under `--individual_zyg` the record stays one feature and every named
//! individual's genotype is read onto it (`create_individuals_zyg_VariationFeature`,
//! `:683-743`): the feature's `ZYG` lists each individual's zygosity, and when
//! the list names one individual carrying no alternate allele the feature is the
//! reference allele alone (`:727-730`, `OutputFactory.pm:506-508`).
//!
//! `--process_ref_homs` adds the reference allele's row to every feature of every
//! sequence record (`OutputFactory.pm:509-516`), first among the feature's
//! alleles and numbered 0 (`VariationFeatureOverlap.pm:469-474`,
//! `BaseVariationFeatureOverlap.pm:363-373`); a structural record is read through
//! another method and gains none (`:534-560`).

use std::collections::HashMap;
use std::sync::Arc;

use vep_core::sv_routing::perl_split;
use vep_core::variant::{
    keeps_uploaded_spelling, InputVariant, SampleZygosity, Zygosity, ZygosityList,
};
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
    /// header lacks is Ensembl VEP's `die` (`:612-617`), both lists in the message,
    /// the header's as `get_samples` lists it. ensembl-io keys the samples by name,
    /// so a name the header repeats is one sample reading its last column, and
    /// `get_samples` lists the names by that column (BaseVCF4.pm:785-808): `all`
    /// reads them so, and the message names `S1 S1 S3` as `S1,S3`.
    pub fn resolve(&self, columns: &SampleColumns) -> Result<Vec<(Arc<str>, usize)>, String> {
        let mut last_column: HashMap<&str, usize> = HashMap::new();
        for (column, name) in columns.names.iter().enumerate() {
            last_column.insert(name.as_str(), column);
        }
        let mut samples: Vec<(Arc<str>, usize)> = last_column
            .iter()
            .map(|(name, column)| (Arc::from(*name), *column))
            .collect();
        samples.sort_by_key(|(_, column)| *column);
        match self {
            IndividualRequest::All => Ok(samples),
            IndividualRequest::Named(names) => {
                let mut individuals = Vec::with_capacity(names.len());
                for name in names {
                    let Some(column) = last_column.get(name.as_str()).copied() else {
                        let header: Vec<&str> =
                            samples.iter().map(|(name, _)| name.as_ref()).collect();
                        return Err(format!(
                            "Sample IDs given ({}) do not match samples from VCF ({})",
                            names.join(","),
                            header.join(",")
                        ));
                    };
                    individuals.push((Arc::from(name.as_str()), column));
                }
                Ok(individuals)
            }
        }
    }
}

/// What a run reads the named individuals' genotypes for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndividualMode {
    /// `--individual`: one copy of the record per individual
    /// (`Parser/VCF.pm:595-665`).
    Copies,
    /// `--individual_zyg`: the record once, every individual's zygosity in its
    /// `ZYG` list (`:683-743`).
    Zygosities,
}

/// How a record's lines are read per individual: the individuals to read, as
/// [`IndividualRequest::resolve`] resolved them against the input's header, with
/// their sample column indices, what their genotypes are read for, and whether
/// the records carry genotypes at all.
#[derive(Debug, Clone)]
pub struct IndividualSplit {
    /// The individuals in their output order, each with its 0-based sample
    /// column (the column after FORMAT is 0).
    individuals: Vec<(Arc<str>, usize)>,
    /// Whether the header's ninth column is `FORMAT`; without it no record yields a
    /// genotype and no copy is written (BaseVCF4.pm:98, :614-618).
    genotyped: bool,
    /// Whether a name is listed twice, so a `ZYG` list has to drop the repeat
    /// (Ensembl VEP keys the entries by name, `Parser/VCF.pm:737`); `all` never
    /// repeats one.
    repeats_a_name: bool,
    mode: IndividualMode,
}

impl IndividualSplit {
    /// The split for `request` over a header's sample `columns`; `Err` is Ensembl
    /// VEP's message for a list naming an unknown sample.
    pub fn new(
        request: &IndividualRequest,
        columns: &SampleColumns,
        mode: IndividualMode,
    ) -> Result<Self, String> {
        let individuals = request.resolve(columns)?;
        let mut names: Vec<&str> = individuals.iter().map(|(name, _)| &**name).collect();
        names.sort_unstable();
        let repeats_a_name = names.windows(2).any(|pair| pair[0] == pair[1]);
        Ok(Self {
            individuals,
            genotyped: columns.genotyped,
            repeats_a_name,
            mode,
        })
    }

    /// The number of individuals listed, a repeated name counted each time: the
    /// copies a record yields at most under `--individual`, and Ensembl VEP's
    /// `$n_individuals` under `--individual_zyg` (`Parser/VCF.pm:710`).
    pub fn len(&self) -> usize {
        self.individuals.len()
    }

    /// Whether no individual is listed.
    pub fn is_empty(&self) -> bool {
        self.individuals.is_empty()
    }

    /// The features a line yields at most, which Ensembl VEP's buffer counts
    /// (InputBuffer.pm:212): one copy per individual listed under `--individual`,
    /// the one record under `--individual_zyg`.
    pub fn features_per_line(&self) -> usize {
        match self.mode {
            IndividualMode::Copies => self.len().max(1),
            IndividualMode::Zygosities => 1,
        }
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
    /// Whether the record is read as one copy per individual (`--individual`),
    /// each copy's alleles the ones its genotype names.
    pub fn reads_copies(&self) -> bool {
        self.individuals
            .is_some_and(|split| split.mode == IndividualMode::Copies)
    }

    /// The alleles the record writes: its individuals' copies, the record with its
    /// individuals' zygosities, or the record itself with the reference row when
    /// `--process_ref_homs` asks for one. `minimal` is the run's `--minimal`, which
    /// decides whether a reading prints the record's uploaded string as
    /// `UPLOADED_ALLELE` ([`keeps_uploaded_spelling`]).
    pub fn finish_record(
        &self,
        mut record: Vec<InputVariant>,
        alleles: &RecordAlleles<'_>,
        columns: &VcfColumns<'_>,
        minimal: bool,
    ) -> Vec<InputVariant> {
        let Some(split) = self.individuals else {
            if self.process_ref_homs {
                add_reference_allele_row(&mut record);
            }
            return record;
        };
        let run = RunFlags {
            process_ref_homs: self.process_ref_homs,
            minimal,
        };
        match split.mode {
            IndividualMode::Copies => split_by_individual(record, alleles, columns, split, run),
            IndividualMode::Zygosities => zygosity_record(record, alleles, columns, split, run),
        }
    }
}

/// The run's flags a per-individual reading reads beyond the split.
#[derive(Debug, Clone, Copy)]
struct RunFlags {
    process_ref_homs: bool,
    /// `--minimal`, which decides whether a reading prints the record's uploaded
    /// string as `UPLOADED_ALLELE` ([`inherits_raw_string`]).
    minimal: bool,
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
/// as a Perl array index (`allele_index`) and the alleles are joined again, so an
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

/// One individual's genotype read against the record's alleles.
struct CalledGenotype<'a> {
    /// Each called bit as the chopped allele it names, an undefined bit the empty
    /// string (`Parser/VCF.pm:627`, `:716`): what `ZYG` counts.
    alleles: Vec<&'a str>,
    /// The record's alleles the individual carries beyond the reference, `grep
    /// {$_ ne $ref && $_ !~ /\*/}` over the mapped bits (`:635`, `:720`), each
    /// distinct string once, as indices into the record's alleles in the record's
    /// order.
    non_ref: Vec<usize>,
    /// A bit past the record's alleles: an undefined allele, which Ensembl VEP's
    /// `grep` keeps (it is neither the reference nor a `*`) though it writes no
    /// row, so the individual is not reference-homozygous.
    empty_allele: bool,
}

impl<'a> CalledGenotype<'a> {
    fn read(called: &[Option<usize>], alleles: &'a RecordAlleles<'_>) -> Self {
        let genotype = called
            .iter()
            .map(|bit| bit.map_or("", |index| alleles.chopped[index].as_str()))
            .collect();
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
        Self {
            alleles: genotype,
            non_ref,
            empty_allele: called.contains(&None),
        }
    }

    /// Whether the individual carries no alternate allele of the record
    /// (`Parser/VCF.pm:641`, `:723`).
    fn hom_ref(&self) -> bool {
        self.non_ref.is_empty() && !self.empty_allele
    }
}

/// The genotype each listed individual calls on the record, in the split's order,
/// each with the individual's position in the list and its name: an individual
/// whose genotype is missing has none ([`called_alleles`]), and a header without
/// `FORMAT` or a record whose FORMAT has no GT key calls none ([`genotype_index`]).
fn called_genotypes<'a>(
    split: &'a IndividualSplit,
    columns: &VcfColumns<'_>,
    alleles: &'a RecordAlleles<'_>,
) -> Vec<(usize, &'a Arc<str>, CalledGenotype<'a>)> {
    if !split.genotyped {
        return Vec::new();
    }
    let Some(gt_index) = genotype_index(columns) else {
        return Vec::new();
    };
    let samples: Vec<&str> = columns.samples().collect();
    let mut genotypes = Vec::new();
    for (ordinal, (name, column)) in split.individuals.iter().enumerate() {
        let Some(gt) = genotype_field(&samples, gt_index, *column) else {
            continue;
        };
        let Some(called) = called_alleles(gt, alleles.written.len()) else {
            continue;
        };
        genotypes.push((ordinal, name, CalledGenotype::read(&called, alleles)));
    }
    genotypes
}

/// Whether every per-individual reading of the record prints the record's raw
/// string as `UPLOADED_ALLELE` (OutputFactory.pm:1360): Ensembl VEP's
/// `original_allele_string`, set for a record whose alleles differ in length and,
/// under `--minimal`, for every record with an alternate allele (`minimise_alleles`
/// passes a string without one through untouched, Parser.pm:974-976), which every
/// copy inherits; the copy's own string otherwise.
fn inherits_raw_string(base: &InputVariant, minimal: bool) -> bool {
    base.uploaded_allele_string
        .as_deref()
        .is_some_and(|raw| keeps_uploaded_spelling(raw, minimal && raw.contains('/')))
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
    genotype: &CalledGenotype<'_>,
    run: RunFlags,
) -> Option<Vec<InputVariant>> {
    let base = record.first()?;
    let inherits_raw_string = inherits_raw_string(base, run.minimal);
    let mark = |variant: &mut InputVariant, copy_string: &str, hom_ref: bool| {
        variant.individual = Some(Arc::clone(name));
        variant.copy_ordinal = ordinal;
        variant.genotype = genotype
            .alleles
            .iter()
            .map(|allele| allele.to_string())
            .collect();
        variant.hom_ref = hom_ref;
        if !inherits_raw_string {
            variant.uploaded_allele_string = Some(copy_string.into());
        }
    };
    if genotype.hom_ref() {
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
        let written_reference = &alleles.chopped[0];
        let (copy_string, written_string) = if run.process_ref_homs {
            (
                format!("{working_reference}/{working_reference}"),
                format!("{written_reference}/{written_reference}"),
            )
        } else {
            (working_reference, written_reference.to_string())
        };
        mark(&mut copy, &copy_string, true);
        copy.record_allele_string_multi = Some(copy_string.into());
        copy.written_case_allele_string = written_case(&written_string);
        return Some(if run.process_ref_homs {
            vec![copy.reference_allele_row()]
        } else {
            vec![copy]
        });
    }
    // The copy's string as written beside its working, upper-cased form: Ensembl
    // VEP names a copy with a Perl-false ID from its own string before the
    // upper-casing (Parser.pm:462-467, then :554), so `0 c t,g` with `0/1` is
    // `21_<pos>_c/t`.
    let mut written_string = alleles.chopped[0].to_string();
    for &index in &genotype.non_ref {
        written_string.push('/');
        written_string.push_str(&alleles.chopped[index]);
    }
    if genotype.empty_allele {
        written_string.push('/');
    }
    let copy_string = written_string.to_ascii_uppercase();
    let mut copy: Vec<InputVariant> = genotype
        .non_ref
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
            variant.written_case_allele_string = written_case(&written_string);
        }
        mark(variant, &copy_string, false);
    }
    if run.process_ref_homs {
        add_reference_allele_row(&mut copy);
    }
    Some(copy)
}

/// `string` as the copy's written-case allele string when a base of it is lower
/// case ([`InputVariant::written_case_allele_string`]); `None` for an upper-case
/// string.
fn written_case(string: &str) -> Option<Arc<str>> {
    string
        .bytes()
        .any(|b| b.is_ascii_lowercase())
        .then(|| string.into())
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
fn split_by_individual(
    record: Vec<InputVariant>,
    alleles: &RecordAlleles<'_>,
    columns: &VcfColumns<'_>,
    split: &IndividualSplit,
    run: RunFlags,
) -> Vec<InputVariant> {
    let mut out = Vec::new();
    for (ordinal, name, genotype) in called_genotypes(split, columns, alleles) {
        if let Some(copy) = individual_copy(&record, alleles, name, ordinal as u32, &genotype, run)
        {
            out.extend(copy);
        }
    }
    out
}

/// The record once, with every named individual's genotype read onto it: the
/// `--individual_zyg` reading of `Parser/VCF.pm:683-743`. One zygosity per
/// individual whose genotype is called, each name once (`:737` keys them by
/// name) in the split's order; an individual whose genotype is missing has none
/// (`:715`), and a header without `FORMAT` or a record whose FORMAT has no GT
/// yields none, so the record prints no `ZYG`. When the list names one
/// individual and that individual carries no alternate allele, the record is the
/// reference allele's row alone, its allele string `REF/REF` (`:727-730`;
/// `OutputFactory.pm:506-508` writes the reference allele of every feature);
/// otherwise it is the record as read, with the reference row under
/// `--process_ref_homs` (`:509-516`).
fn zygosity_record(
    mut record: Vec<InputVariant>,
    alleles: &RecordAlleles<'_>,
    columns: &VcfColumns<'_>,
    split: &IndividualSplit,
    run: RunFlags,
) -> Vec<InputVariant> {
    let mut zygosities: Vec<SampleZygosity> = Vec::new();
    let mut unique_hom_ref = false;
    for (_, name, genotype) in called_genotypes(split, columns, alleles) {
        let hom_ref = genotype.hom_ref();
        unique_hom_ref |= hom_ref && split.len() == 1;
        if split.repeats_a_name && zygosities.iter().any(|z| z.sample == *name) {
            continue;
        }
        zygosities.push(SampleZygosity {
            sample: Arc::clone(name),
            zygosity: Zygosity::of(&genotype.alleles, hom_ref),
        });
    }
    let zygosities = Arc::new(ZygosityList::new(zygosities));
    let Some(base) = record.first() else {
        return record;
    };
    if unique_hom_ref {
        // The record's working reference twice (`$ref` of the post-processed
        // feature, `:691-692`, `:728`), as `individual_copy` reads a
        // reference-homozygous individual's string.
        let working_reference = String::from_utf8_lossy(&base.ref_allele).to_ascii_uppercase();
        let allele_string = format!("{working_reference}/{working_reference}");
        let written_reference = &alleles.chopped[0];
        let mut feature = feature_without_allele(base);
        feature.hom_ref = true;
        if !inherits_raw_string(base, run.minimal) {
            feature.uploaded_allele_string = Some(allele_string.as_str().into());
        }
        feature.record_allele_string_multi = Some(allele_string.into());
        feature.written_case_allele_string =
            written_case(&format!("{written_reference}/{written_reference}"));
        feature.sample_zygosities = Some(zygosities);
        return vec![feature.reference_allele_row()];
    }
    for variant in &mut record {
        variant.sample_zygosities = Some(Arc::clone(&zygosities));
    }
    if run.process_ref_homs {
        add_reference_allele_row(&mut record);
    }
    record
}

/// `--process_ref_homs` on a record, or on an individual's copy of one: the
/// reference allele's row before its alleles, or in place of a feature with no
/// allele (an allele-less or non-variant record, or a copy whose alleles all
/// write no row, `VariationFeatureOverlap.pm:479-491`), since the VCF line of a
/// record holding such a feature carries no CSQ; a structural record gains none.
/// The row follows the shift of the record's last alternate allele
/// ([`InputVariant::reference_row_shift_allele`]).
pub fn add_reference_allele_row(record: &mut Vec<InputVariant>) {
    let Some(first) = record.first() else {
        return;
    };
    if first.is_structural || first.reference_allele_row {
        return;
    }
    let mut row = first.reference_allele_row();
    if let Some(last) = record.last().filter(|v| !v.no_allele) {
        row.reference_row_shift_allele = Some(Arc::from(
            String::from_utf8_lossy(last.alt_allele()).as_ref(),
        ));
    }
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
    /// column, a list may repeat a name (`Parser/VCF.pm:622`), and the refusal of
    /// an unknown name lists the header's samples the same way (`:616` joins
    /// `get_samples`): Ensembl VEP 116.2 names a header `S1 S3 S1` as `(S3,S1)`.
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
        let all = IndividualSplit::new(&IndividualRequest::All, &columns, IndividualMode::Copies)
            .unwrap();
        assert_eq!(names(&all), [("S3".to_string(), 1), ("S1".to_string(), 2)]);
        let twice = IndividualRequest::Named(vec!["S1".into(), "S1".into(), "S3".into()]);
        let split = IndividualSplit::new(&twice, &columns, IndividualMode::Copies).unwrap();
        assert_eq!(
            names(&split),
            [
                ("S1".to_string(), 2),
                ("S1".to_string(), 2),
                ("S3".to_string(), 1)
            ]
        );
        assert_eq!(split.len(), 3);
        let unknown = IndividualRequest::Named(vec!["S9".into(), "S1".into()]);
        let refusal = unknown.resolve(&columns).unwrap_err();
        assert!(
            refusal.ends_with(" given (S9,S1) do not match samples from VCF (S3,S1)"),
            "{refusal}"
        );
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
