// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `--check_frequency` filter: Ensembl VEP's
//! `AnnotationSource/Cache/BaseCacheVariation.pm` at release/116.2, which
//! drops input alleles by the frequency of the known variants co-located with
//! their record before any transcript is annotated (`check_frequency_filter`
//! :271-289 validates the population, `get_frequency_data` :179-256 reads the
//! frequencies, `_add_check_freq_data_to_vf` :394-446 judges each alternate
//! allele, `check_pass` :460-471 is the threshold test and
//! `frequency_check_buffer` :308-375 rewrites the buffer).
//!
//! The filter is a stage of the binary's pipeline, between co-location and
//! annotation as Ensembl VEP's variation source runs before its transcript
//! source (`AnnotationSourceAdaptor.pm:98` sorts a source that can filter
//! first, `AnnotationType/Variation.pm:96` filters after matching). A record
//! whose every alternate allele is dropped leaves the batch; a record that
//! loses some of them keeps the others, with its allele string, `ALLELE_NUM`
//! and `FREQS` as Ensembl VEP writes them. The library [`crate::Annotator`]
//! annotates a batch in place with every allele kept and never filters.
//!
//! `--freq_pop 1KG_ALL`, the default and what `--filter_common` reads, is
//! resolved to the cache's `AF` column: Ensembl VEP 116.2 reads the
//! `1KG_ALL` population from a record's `minor_allele` (:216-218) or its
//! `1KG_ALL` or `ALL` keys (:232), none of which a release-116 cache carries,
//! so it filters nothing under that population; `docs/intended-divergences.md`
//! states the difference.

use std::collections::HashMap;
use std::sync::Arc;

use vep_core::cache_info::CacheInfo;
use vep_core::sv_routing::{perl_split, perl_true};
use vep_core::variant::{ColocatedVariant, InputVariant};
use vep_core::variation::{perl_number_string, perl_numeric, reverse_comp};

use crate::config::{CliError, Config};

/// The frequency of an allele as the filter holds it: the cache's text, or the
/// number `1 - total` of the one allele a record names no frequency for
/// (`BaseCacheVariation.pm:249-252`). Ensembl VEP keeps both in one hash and
/// reads them with Perl's truth, number and string rules.
#[derive(Debug, Clone, PartialEq)]
enum Frequency {
    Text(String),
    Number(f64),
}

impl Frequency {
    /// Perl's truth of the value (`:419`, `:426`): `0`, `0.0` as a number and the
    /// empty text are false, so a frequency of zero reads as no frequency.
    fn is_true(&self) -> bool {
        match self {
            Frequency::Text(text) => perl_true(text),
            Frequency::Number(value) => *value != 0.0,
        }
    }

    fn numeric(&self) -> f64 {
        match self {
            Frequency::Text(text) => perl_numeric(text),
            Frequency::Number(value) => *value,
        }
    }

    /// The value as `FREQS` prints it: the cache's text verbatim, a number as
    /// Perl stringifies it.
    fn text(&self) -> String {
        match self {
            Frequency::Text(text) => text.clone(),
            Frequency::Number(value) => perl_number_string(*value),
        }
    }
}

/// `freq_pop` upper-cased, then each `(from, to)` of `respellings` applied to
/// its first occurrence, in order (Ensembl VEP's `uc` followed by a run of
/// `s///` substitutions).
fn respelled(freq_pop: &str, respellings: &[(&str, &str)]) -> String {
    let mut full = freq_pop.to_ascii_uppercase();
    for (from, to) in respellings {
        full = full.replacen(from, to, 1);
    }
    full
}

/// `--freq_pop` with Ensembl VEP's case rules applied for the lookup
/// (`BaseCacheVariation.pm:183-187`): upper-cased, with `ExAC`, `gnomADg` and
/// `gnomADe` respelled, then split on `_` into the group and the population.
/// Returns the respelled name and the population, the second field of the
/// split; `None` for a name without an `_`.
fn lookup_names(freq_pop: &str) -> (String, Option<String>) {
    let full = respelled(
        freq_pop,
        &[
            ("EXAC", "ExAC"),
            ("GNOMADG", "gnomADg"),
            ("GNOMADE", "gnomADe"),
        ],
    );
    let population = perl_split(&full, &['_']).nth(1).map(str::to_string);
    (full, population)
}

/// Ensembl VEP's check of `--freq_pop` against the cache's declared variation
/// columns (`BaseCacheVariation.pm:271-289`), run when `--check_frequency` is
/// on and the cache declares columns, as its variation source is built
/// (`CacheDir.pm:179-207`): the population must be `1KG_ALL`, a column name, or
/// `1KG_` followed by a column name, once upper-cased with `ExAC`, `Adj` and
/// `gnomADe` respelled. The `gnomADg` respelling of :280 looks for `GNOMADg`,
/// which no upper-cased name contains, so every gnomAD genome population is
/// refused in every spelling, as Ensembl VEP refuses it.
pub fn check_frequency_filter(config: &Config, info: &CacheInfo) -> Result<(), CliError> {
    if !config.check_frequency {
        return Ok(());
    }
    let Some(columns) = info.variation_columns() else {
        return Ok(());
    };
    let freq_pop = config.freq_pop.as_str();
    if !perl_true(freq_pop) {
        return Err(CliError::new("cache", "No population specified"));
    }
    let full = respelled(
        freq_pop,
        &[
            ("EXAC", "ExAC"),
            ("ADJ", "Adj"),
            ("GNOMADE", "gnomADe"),
            ("GNOMADg", "gnomADg"),
        ],
    );
    let mut fields = perl_split(&full, &['_']);
    let group = fields.next().unwrap_or("");
    let population = fields.next().unwrap_or("");
    let valid = full == "1KG_ALL"
        || columns.position(&full).is_some()
        || (group == "1KG" && columns.position(population).is_some());
    if valid {
        Ok(())
    } else {
        Err(CliError::new(
            "cache",
            format!("Invalid population {freq_pop} specified"),
        ))
    }
}

/// The filter a run's `--check_frequency` flags describe.
#[derive(Debug, Clone, PartialEq)]
pub struct FrequencyFilter {
    /// `--freq_pop` upper-cased: the population `FREQS` names
    /// (`BaseCacheVariation.pm:402`, `OutputFactory.pm:1156`).
    freqs_key: String,
    /// The cache column read first on every co-located variant (:232), the
    /// respelled `--freq_pop`, or `AF` for `1KG_ALL`.
    column: String,
    /// The column read when the first holds nothing: the population half of
    /// `--freq_pop` (:232 reads `$ex->{$freq_pop}` after `$ex->{$freq_pop_full}`,
    /// so `gnomADe_AFR` falls back to the 1000 Genomes `AFR` column).
    fallback: Option<String>,
    /// `--freq_freq`.
    threshold: f64,
    /// `--freq_gt_lt gt`: an allele passes at or above the threshold; `lt`, at
    /// or below (:467-468).
    at_or_above: bool,
    /// `--freq_filter include`: a record is kept only when every alternate
    /// allele passes; `exclude` drops the alleles that pass (:319-359).
    include: bool,
}

impl FrequencyFilter {
    /// The filter of a configuration with `--check_frequency` on; `None`
    /// without it.
    pub fn from_config(config: &Config) -> Option<Self> {
        if !config.check_frequency {
            return None;
        }
        let (full, population) = lookup_names(&config.freq_pop);
        let (column, fallback) = if full == "1KG_ALL" {
            ("AF".to_string(), None)
        } else {
            (full, population)
        };
        Some(Self {
            freqs_key: config.freq_pop.to_ascii_uppercase(),
            column,
            fallback,
            threshold: config.freq_freq,
            at_or_above: config.freq_gt_lt == "gt",
            include: config.freq_filter == "include",
        })
    }

    /// Whether a record with no frequency data is dropped (`--freq_filter
    /// include`, `BaseCacheVariation.pm:322`); under `exclude` it is kept.
    pub fn drops_records_without_data(&self) -> bool {
        self.include
    }

    /// Whether a record with no alternate allele is dropped: the feature of a
    /// non-variant line kept by `--allow_non_variant`, the reference alone
    /// (`Parser/VCF.pm:342`). Under `include` it is dropped without frequency
    /// data (`BaseCacheVariation.pm:322`) and kept with it, since
    /// `_add_check_freq_data_to_vf` marks a record passed when its pass count
    /// equals its allele count (:444) and no allele can fail; under `exclude`
    /// it is kept either way. The data is its co-located variants' as
    /// `get_frequency_data` reads them: a null-allele record at its span
    /// contributes the remainder `1` for its one allele (:249-252).
    pub fn drops_record_without_alleles(&self, record: &[InputVariant]) -> bool {
        self.include && !record.is_empty() && self.frequency_data(record).0.is_empty()
    }

    /// `check_pass` (`BaseCacheVariation.pm:460-471`): inclusive on both sides.
    fn passes(&self, frequency: &Frequency) -> bool {
        let value = frequency.numeric();
        if self.at_or_above {
            value >= self.threshold
        } else {
            value <= self.threshold
        }
    }

    /// `get_frequency_data` (`BaseCacheVariation.pm:179-256`) over the
    /// record's co-located variants in window order (the order of Ensembl VEP's
    /// `existing` list, which the per-allele lists of a record whose alleles
    /// span different positions interleave by start): the frequency of every
    /// allele named by the population's column of any variant, the last variant
    /// naming an allele winning, each variant's one unnamed allele given the
    /// remainder `1 - total`; and the record allele each input allele was
    /// matched to, the last variant matching it winning. A variant on the other
    /// strand has its alleles reverse-complemented first. Empty when no variant
    /// names a frequency and none has a single unnamed allele.
    fn frequency_data(
        &self,
        record: &[InputVariant],
    ) -> (HashMap<String, Frequency>, HashMap<String, String>) {
        let mut frequencies: HashMap<String, Frequency> = HashMap::new();
        let mut matched: HashMap<String, String> = HashMap::new();
        let record_strand = record[0].strand.as_i8();
        let mut colocated: Vec<&ColocatedVariant> = Vec::new();
        for cv in record.iter().flat_map(|v| &v.colocated_variants) {
            if !colocated.iter().any(|seen| seen.id == cv.id) {
                colocated.push(cv);
            }
        }
        colocated.sort_by_key(|cv| cv.start);
        for cv in colocated {
            let same_strand = cv.strand == record_strand;
            let oriented = |allele: &str| -> String {
                if same_strand {
                    allele.to_string()
                } else {
                    reverse_comp(allele)
                }
            };
            // `%ex_alleles` (:199-203): the record's alleles as a set.
            let mut record_alleles: Vec<String> = Vec::new();
            for allele in perl_split(cv.allele_string.as_deref().unwrap_or(""), &['/']) {
                let allele = oriented(allele);
                if !record_alleles.contains(&allele) {
                    record_alleles.push(allele);
                }
            }
            for pair in &cv.matched_alleles {
                matched.insert(pair.input_allele.clone(), pair.record_allele.clone());
            }
            let mut total = 0.0;
            let column_text = [Some(self.column.as_str()), self.fallback.as_deref()]
                .into_iter()
                .flatten()
                .filter_map(|column| cv.frequency_strings.get(column))
                .find(|text| perl_true(text));
            if let Some(text) = column_text {
                for pair in perl_split(text, &[',']) {
                    let mut fields = pair.splitn(3, ':');
                    let allele = oriented(fields.next().unwrap_or(""));
                    let value = fields.next().unwrap_or("").to_string();
                    total += perl_numeric(&value);
                    record_alleles.retain(|a| *a != allele);
                    frequencies.insert(allele, Frequency::Text(value));
                }
            }
            if let [remaining] = record_alleles.as_slice() {
                frequencies.insert(remaining.clone(), Frequency::Number(1.0 - total));
            }
        }
        (frequencies, matched)
    }

    /// Judges one input record, the alleles split from one line, and returns
    /// which of them stay (`frequency_check_buffer`,
    /// `BaseCacheVariation.pm:308-359`, with `_add_check_freq_data_to_vf`,
    /// :394-446). Every allele of a record with frequency data gets the
    /// record's `FREQS` entries; the survivors of a record that lost alleles
    /// get the rewritten allele string and their new allele numbers. An empty
    /// record has an empty verdict.
    ///
    /// Under `include` a record without frequency data is dropped, and so is
    /// one with a failing allele: Ensembl VEP rewrites its allele string to the
    /// passing alleles and then leaves it out of the buffer (:325-334). Under
    /// `exclude` a record without data, or whose every allele fails, is kept
    /// whole; one whose every allele passes is dropped; one with both keeps
    /// its failing alleles (:341-358).
    ///
    /// The alleles judged are those of the record's allele string (:404,
    /// `$vf->alt_alleles`), so a `*` or `.` ALT, which the parser reads as no
    /// variant, has no frequency and fails. When it alone stays, the record
    /// stays as Ensembl VEP keeps its feature: with no allele to annotate, the
    /// rewritten allele string, and its co-located variants.
    pub fn filter_record(&self, record: &mut [InputVariant]) -> Vec<bool> {
        if record.is_empty() {
            return Vec::new();
        }
        let (frequencies, matched) = self.frequency_data(record);
        if frequencies.is_empty() {
            return vec![!self.include; record.len()];
        }
        let alleles = record_alleles(record);
        let mut passes = Vec::with_capacity(alleles.len());
        let mut freqs = Vec::with_capacity(alleles.len());
        for (allele, _) in &alleles {
            // `$freq_data->{$alt} || 'NA'`, then the matched record allele's
            // value when the allele's own is not Perl-true (:419-429).
            let mut frequency = frequencies.get(allele).filter(|f| f.is_true());
            if frequency.is_none() {
                frequency = matched
                    .get(allele)
                    .filter(|record_allele| perl_true(record_allele))
                    .and_then(|record_allele| frequencies.get(record_allele))
                    .filter(|f| f.is_true());
            }
            let text = frequency.map_or_else(|| "NA".to_string(), Frequency::text);
            freqs.push(format!("{}:{allele}:{text}", self.freqs_key));
            passes.push(frequency.is_some_and(|f| self.passes(f)));
        }
        let passed = passes.iter().filter(|&&p| p).count();
        let all_passed = passed == alleles.len();
        let all_failed = passed == 0;
        let stays: Vec<bool> = if self.include {
            vec![all_passed; alleles.len()]
        } else if all_failed {
            vec![true; alleles.len()]
        } else {
            passes.iter().map(|&p| !p).collect()
        };
        for variant in record.iter_mut() {
            variant.freqs = freqs.clone();
        }
        let mut keep = vec![false; record.len()];
        for ((_, variant), &stay) in alleles.iter().zip(&stays) {
            if let Some(index) = variant {
                keep[*index] = stay;
            }
        }
        if stays.iter().any(|&s| s) && !stays.iter().all(|&s| s) {
            rewrite_allele_string(record, &alleles, &stays, &mut keep);
        }
        keep
    }
}

/// The alternate alleles of the record's allele string in order, each with the
/// index of the variant that carries it; `None` for an allele no variant
/// carries (`*` and `.`, which the parser keeps in the string and reads as no
/// allele). A spelling shared by several alleles is carried by its variants in
/// order, and a variant the string does not spell follows, under its own
/// allele, so every variant is judged.
fn record_alleles(record: &[InputVariant]) -> Vec<(String, Option<usize>)> {
    let displayed: Vec<String> = record.iter().map(InputVariant::display_allele).collect();
    let mut taken = vec![false; record.len()];
    let allele_string = record[0].record_allele_string();
    let mut alleles: Vec<(String, Option<usize>)> = perl_split(&allele_string, &['/'])
        .skip(1)
        .map(|allele| {
            let index = (0..record.len()).find(|&i| !taken[i] && displayed[i] == allele);
            if let Some(i) = index {
                taken[i] = true;
            }
            (allele.to_string(), index)
        })
        .collect();
    for (index, allele) in displayed.into_iter().enumerate() {
        if !taken[index] {
            alleles.push((allele, Some(index)));
        }
    }
    alleles
}

/// The record's allele string after `frequency_check_buffer` dropped alleles
/// (`BaseCacheVariation.pm:355`, `join('/', ref_allele_string, @alts)`): the
/// reference and the alleles that stay, set on every kept variant with its new
/// `ALLELE_NUM`, which Ensembl VEP numbers from the rewritten string
/// (`VariationFeatureOverlapAllele::allele_number` of ensembl-variation). When
/// no variant stays, the first variant stays in its place marked as having no
/// allele, spelled as the first staying allele so that nothing reads it as the
/// dropped one.
fn rewrite_allele_string(
    record: &mut [InputVariant],
    alleles: &[(String, Option<usize>)],
    stays: &[bool],
    keep: &mut [bool],
) {
    let reference = record[0]
        .record_allele_string()
        .split('/')
        .next()
        .unwrap_or("")
        .to_string();
    let mut parts = vec![reference];
    let mut number = 0u32;
    let mut numbers: Vec<(usize, u32)> = Vec::new();
    for ((allele, variant), _) in alleles.iter().zip(stays).filter(|(_, &stay)| stay) {
        parts.push(allele.clone());
        number += 1;
        if let Some(index) = variant {
            numbers.push((*index, number));
        }
    }
    let rewritten: Arc<str> = Arc::from(parts.join("/"));
    if numbers.is_empty() {
        let placeholder = &mut record[0];
        placeholder.no_allele = true;
        placeholder.alt_alleles = vec![parts[1].as_bytes().to_vec()];
        placeholder.allele_num = None;
        placeholder.filtered_allele_string = Some(rewritten);
        keep[0] = true;
        return;
    }
    for (index, number) in numbers {
        record[index].allele_num = Some(number);
        record[index].filtered_allele_string = Some(Arc::clone(&rewritten));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::cache_info::CacheInfo;
    use vep_core::variant::MatchedAllele;

    const COLUMNS: &[&str] = &[
        "chr",
        "variation_name",
        "failed",
        "somatic",
        "start",
        "end",
        "allele_string",
        "strand",
        "clin_sig",
        "phenotype_or_disease",
        "clin_sig_allele",
        "clin_sig_ref_allele",
        "clinical_impact",
        "pubmed",
        "var_synonyms",
        "AF",
        "AFR",
        "AMR",
        "EAS",
        "EUR",
        "SAS",
        "gnomADe",
        "gnomADe_AFR",
        "gnomADe_NFE",
        "gnomADg",
        "gnomADg_AFR",
    ];

    fn info() -> CacheInfo {
        CacheInfo {
            variation_cols: COLUMNS.iter().map(|c| c.to_string()).collect(),
            ..CacheInfo::default()
        }
    }

    fn config(args: &[&str]) -> Config {
        let mut argv = vec!["vep"];
        argv.extend_from_slice(args);
        Config::from_argv(argv).expect("a valid command line")
    }

    fn filter(args: &[&str]) -> FrequencyFilter {
        FrequencyFilter::from_config(&config(args)).expect("--check_frequency is on")
    }

    /// A known variant `name` with `allele_string` and the given frequency
    /// columns, matched to the input alleles `matched` as `(input, record)` pairs.
    fn known(
        name: &str,
        allele_string: &str,
        frequencies: &[(&str, &str)],
        matched: &[(&str, &str)],
    ) -> ColocatedVariant {
        ColocatedVariant {
            id: name.to_string(),
            start: 100,
            end: 100,
            allele_string: Some(allele_string.to_string()),
            strand: 1,
            frequency_strings: frequencies
                .iter()
                .map(|(column, text)| (column.to_string(), text.to_string()))
                .collect(),
            matched_alleles: matched
                .iter()
                .enumerate()
                .map(|(i, (input, record))| MatchedAllele {
                    input_allele: input.to_string(),
                    input_index: i,
                    record_allele: record.to_string(),
                    record_index: i,
                })
                .collect(),
            ..ColocatedVariant::default()
        }
    }

    /// The alleles of one record, `reference` against each of `alts`, each
    /// carrying the record's co-located variants.
    fn alleles_of(
        reference: &str,
        alts: &[&str],
        colocated: Vec<ColocatedVariant>,
    ) -> Vec<InputVariant> {
        let multi: Arc<str> = Arc::from(format!("{reference}/{}", alts.join("/")));
        alts.iter()
            .enumerate()
            .map(|(index, alt)| {
                let mut v = InputVariant::new(
                    "21".into(),
                    100,
                    100,
                    reference.as_bytes().to_vec(),
                    alt.as_bytes().to_vec(),
                );
                v.allele_index = index;
                v.input_record = 1;
                if alts.len() > 1 {
                    v.record_allele_string_multi = Some(Arc::clone(&multi));
                }
                v.colocated_variants = colocated.clone();
                v
            })
            .collect()
    }

    fn one(reference: &str, alt: &str, colocated: Vec<ColocatedVariant>) -> Vec<InputVariant> {
        alleles_of(reference, &[alt], colocated)
    }

    /// Ensembl VEP 116.2 on the release-116 GRCh38 cache, `21 5131441 T C` with
    /// `rs1555830109` (`AF` `C:0.006789`), `21 7759128 T A,G` with `rs1555870136`
    /// (`AF` `G:1`) and `21 5030192 A G` with `rs1172282560` (gnomAD genome
    /// columns alone) under `--check_frequency --freq_pop AF --freq_freq 0.01
    /// --freq_gt_lt gt --freq_filter exclude`: the first prints
    /// `FREQS=AF:C:0.006789`, the second keeps allele `A` alone with
    /// `FREQS=AF:A:NA,AF:G:1` and `ALLELE_NUM=1`, the third prints no `FREQS`
    /// and keeps its row.
    #[test]
    fn exclude_keeps_variants_without_data() {
        let f = filter(&["--check_frequency", "--freq_pop", "AF"]);
        let mut no_data = one(
            "A",
            "G",
            vec![known(
                "rs1172282560",
                "A/G",
                &[("gnomADg", "G:0")],
                &[("G", "G")],
            )],
        );
        assert_eq!(f.filter_record(&mut no_data), [true]);
        assert!(no_data[0].freqs.is_empty());
        assert!(no_data[0].filtered_allele_string.is_none());

        let mut rare = one(
            "T",
            "C",
            vec![known(
                "rs1555830109",
                "T/C",
                &[("AF", "C:0.006789")],
                &[("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut rare), [true]);
        assert_eq!(rare[0].freqs, ["AF:C:0.006789"]);

        let mut common = one(
            "T",
            "G",
            vec![known(
                "rs1555870136",
                "T/A/G",
                &[("AF", "G:1")],
                &[("G", "G")],
            )],
        );
        assert_eq!(f.filter_record(&mut common), [false]);
    }

    /// `21 7759128 T A,G` against `rs1555870136` (`AF` `G:1`) under the exclude
    /// run above: Ensembl VEP keeps the failing `A` alone, rewrites the record
    /// to `T/A` (JSON `allele_string`, `UPLOADED_ALLELE`) and numbers it 1; its
    /// `FREQS` lists both alleles. `21 44195140 G A,C` against `rs2838519`
    /// (`AF` `A:0.5625`) keeps `C` as allele 1 of `G/C` with
    /// `FREQS=AF:C:NA,AF:A:0.5625`.
    #[test]
    fn partial_pass_rewrites_allele_string() {
        let f = filter(&["--check_frequency", "--freq_pop", "AF"]);
        let mut record = alleles_of(
            "T",
            &["A", "G"],
            vec![known(
                "rs1555870136",
                "T/A/G",
                &[("AF", "G:1")],
                &[("A", "A"), ("G", "G")],
            )],
        );
        assert_eq!(f.filter_record(&mut record), [true, false]);
        assert_eq!(record[0].freqs, ["AF:A:NA", "AF:G:1"]);
        assert_eq!(record[0].filtered_allele_string.as_deref(), Some("T/A"));
        assert_eq!(record[0].allele_num, Some(1));

        let mut record = alleles_of(
            "G",
            &["A", "C"],
            vec![known(
                "rs2838519",
                "G/A/C",
                &[("AF", "A:0.5625")],
                &[("A", "A"), ("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut record), [false, true]);
        assert_eq!(record[1].freqs, ["AF:A:0.5625", "AF:C:NA"]);
        assert_eq!(record[1].filtered_allele_string.as_deref(), Some("G/C"));
        assert_eq!(record[1].allele_num, Some(1));
        assert_eq!(record[1].display_allele(), "C");
    }

    /// A library caller may hand the filter an empty record; the binary's
    /// record units are never empty.
    #[test]
    fn empty_record_has_an_empty_verdict() {
        let f = filter(&["--check_frequency", "--freq_pop", "AF"]);
        assert_eq!(f.filter_record(&mut []), Vec::<bool>::new());
    }

    /// A `*` ALT is an allele of the record's string with no variant behind it.
    /// Ensembl VEP 116.2 judges it with the others (`$vf->alt_alleles`, :404): on
    /// `21 8529128 T C,*` against `rs66525445` (`AF` `C:0.1939,G:0.5665`) it
    /// prints `FREQS=AF:C:0.1939,AF:*:NA`; under `--freq_filter include` the
    /// failing `*` drops the record; under `exclude` at `--freq_freq 0.01` the
    /// passing `C` goes and the record stays as the feature `T/*` with no
    /// allele to annotate (the JSON object has `allele_string` `T/*` and
    /// `most_severe_consequence` `?`, the VCF line carries no CSQ); at
    /// `--freq_freq 0.3` both fail and the record stays whole. A record with a
    /// passing, a failing and a `*` allele keeps the failing one as allele 1 of
    /// the string that still names the `*`.
    #[test]
    fn star_allele_is_a_failing_allele_of_the_record() {
        let star_record = |alts: &str| {
            let mut record = one(
                "T",
                "C",
                vec![known(
                    "rs66525445",
                    "T/A/C/G",
                    &[("AF", "C:0.1939,G:0.5665")],
                    &[("C", "C")],
                )],
            );
            record[0].record_allele_string_multi = Some(Arc::from(alts));
            record
        };
        let exclude = filter(&["--check_frequency", "--freq_pop", "AF"]);
        let mut record = star_record("T/C/*");
        assert_eq!(exclude.filter_record(&mut record), [true]);
        assert_eq!(record[0].freqs, ["AF:C:0.1939", "AF:*:NA"]);
        assert!(record[0].no_allele);
        assert_eq!(record[0].display_allele(), "*");
        assert_eq!(record[0].filtered_allele_string.as_deref(), Some("T/*"));
        assert_eq!(record[0].allele_num, None);

        let include = filter(&[
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_filter",
            "include",
        ]);
        let mut record = star_record("T/C/*");
        assert_eq!(include.filter_record(&mut record), [false]);
        assert_eq!(record[0].freqs, ["AF:C:0.1939", "AF:*:NA"]);
        assert!(!record[0].no_allele);

        let rare = filter(&[
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_freq",
            "0.3",
        ]);
        let mut record = star_record("T/C/*");
        assert_eq!(rare.filter_record(&mut record), [true]);
        assert_eq!(record[0].freqs, ["AF:C:0.1939", "AF:*:NA"]);
        assert!(!record[0].no_allele);
        assert_eq!(record[0].display_allele(), "C");
        assert!(record[0].filtered_allele_string.is_none());

        let mut mixed = alleles_of(
            "G",
            &["A", "C"],
            vec![known(
                "rs2838519",
                "G/A/C",
                &[("AF", "A:0.5625")],
                &[("A", "A"), ("C", "C")],
            )],
        );
        for variant in &mut mixed {
            variant.record_allele_string_multi = Some(Arc::from("G/A/C/*"));
        }
        assert_eq!(exclude.filter_record(&mut mixed), [false, true]);
        assert_eq!(mixed[1].freqs, ["AF:A:0.5625", "AF:C:NA", "AF:*:NA"]);
        assert_eq!(mixed[1].filtered_allele_string.as_deref(), Some("G/C/*"));
        assert_eq!(mixed[1].allele_num, Some(1));
        assert!(!mixed[1].no_allele);
    }

    /// Under `--freq_filter include` Ensembl VEP 116.2 drops a record with no
    /// frequency data (`BaseCacheVariation.pm:322`) and a record with a failing
    /// allele beside a passing one (:325-334 rewrites it and leaves it out of
    /// the buffer): under `--freq_pop AF --freq_freq 0.01 --freq_gt_lt gt
    /// --freq_filter include` it writes no row for `21 44195140 G A,C` (`A`
    /// 0.5625 passes, `C` has no value). A record whose every allele passes is
    /// kept whole.
    #[test]
    fn include_drops_variants_without_frequency_data() {
        let f = filter(&[
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_filter",
            "include",
        ]);
        let mut no_data = one("A", "G", vec![]);
        assert_eq!(f.filter_record(&mut no_data), [false]);

        let mut partial = alleles_of(
            "G",
            &["A", "C"],
            vec![known(
                "rs2838519",
                "G/A/C",
                &[("AF", "A:0.5625")],
                &[("A", "A"), ("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut partial), [false, false]);

        let mut common = one(
            "T",
            "G",
            vec![known(
                "rs1555870136",
                "T/A/G",
                &[("AF", "G:1")],
                &[("G", "G")],
            )],
        );
        assert_eq!(f.filter_record(&mut common), [true]);
        assert_eq!(common[0].freqs, ["AF:G:1"]);
        assert!(common[0].filtered_allele_string.is_none());
    }

    /// A record with no alternate allele, the feature of a non-variant line
    /// (`Parser/VCF.pm:342`), is dropped under `include` only without frequency
    /// data: a null-allele known variant at its span supplies the remainder `1`
    /// (:249-252), and `_add_check_freq_data_to_vf` marks a record with no
    /// allele to fail passed (:444). Ensembl VEP 116.2 under
    /// `--allow_non_variant --vcf --check_frequency --freq_pop AF --freq_filter
    /// include` writes `21 10521608 . C .` and `21 10521659 . A .`, each at a
    /// `COSMIC_MUTATION` record, and drops `21 5131442 . C .`, where no known
    /// variant lies; `exclude` writes all three.
    #[test]
    fn record_without_alleles_is_dropped_under_include_only_without_data() {
        let include = filter(&[
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_filter",
            "include",
        ]);
        let exclude = filter(&["--check_frequency", "--freq_pop", "AF"]);
        let cosmic = known("COSV107315996", "COSMIC_MUTATION", &[], &[]);
        let feature = |colocated: Vec<ColocatedVariant>| {
            let mut v =
                InputVariant::without_allele("21".into(), 10521608, "C", std::iter::empty());
            v.colocated_variants = colocated;
            vec![v]
        };
        assert!(include.drops_record_without_alleles(&feature(vec![])));
        assert!(!include.drops_record_without_alleles(&feature(vec![cosmic.clone()])));
        assert!(!exclude.drops_record_without_alleles(&feature(vec![])));
        assert!(!exclude.drops_record_without_alleles(&feature(vec![cosmic])));
        assert!(!include.drops_record_without_alleles(&[]));
    }

    /// `check_pass` (`BaseCacheVariation.pm:467-468`) passes at the threshold
    /// on either side: `$f >= $freq_freq` under `gt`, `$f <= $freq_freq` under
    /// `lt`.
    #[test]
    fn threshold_is_inclusive_both_ways() {
        let at = |gt_lt: &str| {
            let f = filter(&[
                "--check_frequency",
                "--freq_pop",
                "AF",
                "--freq_freq",
                "0.05",
                "--freq_gt_lt",
                gt_lt,
                "--freq_filter",
                "include",
            ]);
            let mut record = one(
                "T",
                "C",
                vec![known("rs1", "T/C", &[("AF", "C:0.05")], &[("C", "C")])],
            );
            f.filter_record(&mut record)
        };
        assert_eq!(at("gt"), [true]);
        assert_eq!(at("lt"), [true]);
    }

    /// Under `lt` an allele passes at or below the threshold and fails above
    /// it; a frequency of `0` is no frequency (`$freq_data->{$alt} || 'NA'`,
    /// :419), so an allele the cache lists at `0` prints `NA` and fails. Ensembl
    /// VEP 116.2 under `--freq_pop gnomADe_AFR --freq_freq 0.05 --freq_gt_lt lt
    /// --freq_filter exclude` prints `FREQS=GNOMADE_AFR:C:NA` for `21 5131441 T C`
    /// (`rs1555830109`, `gnomADe_AFR` `C:0`) and keeps it.
    #[test]
    fn lt_passes_at_or_below_threshold() {
        let f = filter(&[
            "--check_frequency",
            "--freq_pop",
            "gnomADe_AFR",
            "--freq_freq",
            "0.05",
            "--freq_gt_lt",
            "lt",
            "--freq_filter",
            "include",
        ]);
        let mut below = one(
            "T",
            "C",
            vec![known(
                "rs1",
                "T/C",
                &[("gnomADe_AFR", "C:0.02")],
                &[("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut below), [true]);
        assert_eq!(below[0].freqs, ["GNOMADE_AFR:C:0.02"]);

        let mut above = one(
            "T",
            "C",
            vec![known(
                "rs1",
                "T/C",
                &[("gnomADe_AFR", "C:0.7")],
                &[("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut above), [false]);

        let mut zero = one(
            "T",
            "C",
            vec![known(
                "rs1555830109",
                "T/C",
                &[("gnomADe_AFR", "C:0")],
                &[("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut zero), [false]);
        assert_eq!(zero[0].freqs, ["GNOMADE_AFR:C:NA"]);
    }

    /// `check_frequency_filter` (`BaseCacheVariation.pm:271-289`) accepts
    /// `1KG_ALL`, any declared column in any case (`af`, `AFR`, `gnomade_nfe`
    /// respelled `gnomADe_NFE`) and `1KG_` before a column (`1KG_AFR`,
    /// `1kg_all`), each of which Ensembl VEP 116.2 accepts on the release-116
    /// GRCh38 cache.
    #[test]
    fn freq_pop_accepts_cache_columns_and_1kg_prefix() {
        for pop in [
            "1KG_ALL",
            "1kg_all",
            "AF",
            "af",
            "AFR",
            "1KG_AFR",
            "1kg_afr",
            "gnomADe_AFR",
            "gnomade_afr",
            "GNOMADE_NFE",
            "gnomADe",
            "1KG_gnomADe_AFR",
        ] {
            let config = config(&["--check_frequency", "--freq_pop", pop]);
            assert_eq!(check_frequency_filter(&config, &info()), Ok(()), "{pop}");
        }
        let off = config(&["--check_existing", "--freq_pop", "nowhere"]);
        assert_eq!(check_frequency_filter(&off, &info()), Ok(()));
        let no_columns = config(&["--check_frequency", "--freq_pop", "nowhere"]);
        assert_eq!(
            check_frequency_filter(&no_columns, &CacheInfo::default()),
            Ok(())
        );
    }

    /// `ERROR: Invalid population <pop> specified` (:285) for a name that is no
    /// column once upper-cased (`chr`, `AFR_AF`, `nowhere`) and for every
    /// gnomAD genome population (`gnomADg_AFR`, `gnomadg_afr`, `gnomadg`), since
    /// the respelling at :280 matches nothing after `uc`, each of which Ensembl
    /// VEP 116.2 refuses on the release-116 GRCh38 cache. An empty or `0`
    /// population is `ERROR: No population specified` (:274).
    #[test]
    fn freq_pop_rejects_unknown_population_with_ensembl_message() {
        for pop in [
            "gnomADg_AFR",
            "gnomadg_afr",
            "gnomadg",
            "chr",
            "AFR_AF",
            "nowhere",
        ] {
            let config = config(&["--check_frequency", "--freq_pop", pop]);
            let error = check_frequency_filter(&config, &info()).unwrap_err();
            assert_eq!(
                error.to_string(),
                format!("ERROR: [cache] Invalid population {pop} specified")
            );
        }
        for pop in ["", "0"] {
            let config = config(&["--check_frequency", "--freq_pop", pop]);
            let error = check_frequency_filter(&config, &info()).unwrap_err();
            assert_eq!(error.to_string(), "ERROR: [cache] No population specified");
        }
    }

    /// `--filter_common` is `--check_frequency` under the defaults of
    /// `Config.pm:319-322` (`freq_pop 1KG_ALL`, `freq_freq 0.01`, `freq_gt_lt gt`,
    /// `freq_filter exclude`), and `1KG_ALL` reads the cache's `AF` column with
    /// `FREQS` naming `1KG_ALL`; Ensembl VEP 116.2 reads no column for it.
    #[test]
    fn filter_common_equals_check_frequency_defaults_with_af_column() {
        let common = filter(&["--filter_common"]);
        let explicit = filter(&[
            "--check_frequency",
            "--freq_pop",
            "1KG_ALL",
            "--freq_freq",
            "0.01",
            "--freq_gt_lt",
            "gt",
            "--freq_filter",
            "exclude",
        ]);
        assert_eq!(common, explicit);
        assert_eq!(common.column, "AF");
        assert_eq!(common.fallback, None);
        assert_eq!(common.freqs_key, "1KG_ALL");
        assert!(config(&["--filter_common"]).check_existing);

        let mut record = one(
            "T",
            "C",
            vec![known(
                "rs1555830109",
                "T/C",
                &[("AF", "C:0.006789")],
                &[("C", "C")],
            )],
        );
        assert_eq!(common.filter_record(&mut record), [true]);
        assert_eq!(record[0].freqs, ["1KG_ALL:C:0.006789"]);
        let mut common_record = one(
            "T",
            "G",
            vec![known(
                "rs1555870136",
                "T/A/G",
                &[("AF", "G:1")],
                &[("G", "G")],
            )],
        );
        assert_eq!(common.filter_record(&mut common_record), [false]);
    }

    /// The lookup reads `$ex->{$freq_pop_full} || $ex->{$freq_pop}` (:232): a
    /// record with no `gnomADe_AFR` value is judged on its 1000 Genomes `AFR`
    /// column. Ensembl VEP 116.2 under `--freq_pop gnomADe_AFR --freq_freq 0.05
    /// --freq_gt_lt lt --freq_filter exclude` prints
    /// `FREQS=GNOMADE_AFR:A:0.7148,GNOMADE_AFR:C:NA` for `21 44195140 G A,C`
    /// against `rs2838519` (`AFR` `A:0.7148`, no gnomAD exome columns) and keeps
    /// both alleles, since neither passes.
    #[test]
    fn population_half_of_the_name_is_the_fallback_column() {
        let f = filter(&[
            "--check_frequency",
            "--freq_pop",
            "gnomADe_AFR",
            "--freq_freq",
            "0.05",
            "--freq_gt_lt",
            "lt",
        ]);
        assert_eq!(f.column, "gnomADe_AFR");
        assert_eq!(f.fallback.as_deref(), Some("AFR"));
        let mut record = alleles_of(
            "G",
            &["A", "C"],
            vec![known(
                "rs2838519",
                "G/A/C",
                &[("AFR", "A:0.7148"), ("gnomADg_AFR", "A:0.001167,C:0")],
                &[("A", "A"), ("C", "C")],
            )],
        );
        assert_eq!(f.filter_record(&mut record), [true, true]);
        assert_eq!(
            record[0].freqs,
            ["GNOMADE_AFR:A:0.7148", "GNOMADE_AFR:C:NA"]
        );
        assert!(record[0].filtered_allele_string.is_none());
    }

    /// The one allele a record names no frequency for takes `1 - total`
    /// (:249-252), printed as Perl prints the number, and a record allele
    /// matched to an input allele spelled otherwise supplies its value through
    /// the match (:425-428): a deletion `GG/G` read as `-` finds the record's
    /// `G`. A record on the other strand is read through its reverse
    /// complement (:201, :237).
    #[test]
    fn remainder_and_matched_allele_supply_a_frequency() {
        let f = filter(&[
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_freq",
            "0.5",
            "--freq_filter",
            "include",
        ]);
        let mut remainder = one(
            "T",
            "C",
            vec![known("rs1", "T/C", &[("AF", "T:0.3")], &[("C", "C")])],
        );
        assert_eq!(f.filter_record(&mut remainder), [true]);
        assert_eq!(remainder[0].freqs, ["AF:C:0.7"]);

        let mut deletion = InputVariant::new("21".into(), 101, 101, b"G".to_vec(), b"-".to_vec());
        deletion.input_record = 1;
        deletion.colocated_variants = vec![known("rs2", "GG/G", &[("AF", "G:0.9")], &[("-", "G")])];
        let mut deletion = vec![deletion];
        assert_eq!(f.filter_record(&mut deletion), [true]);
        assert_eq!(deletion[0].freqs, ["AF:-:0.9"]);

        let mut reverse = known("rs3", "A/G", &[("AF", "G:0.8")], &[("C", "G")]);
        reverse.strand = -1;
        let mut other_strand = one("T", "C", vec![reverse]);
        assert_eq!(f.filter_record(&mut other_strand), [true]);
        assert_eq!(other_strand[0].freqs, ["AF:C:0.8"]);
    }

    /// `%freq_data` is filled from the record's `existing` list in window order
    /// (:195, :238), so of two known variants naming one allele the later in the
    /// window wins, whichever of the record's alleles each was found on.
    #[test]
    fn later_variant_in_the_window_wins_across_alleles() {
        let f = filter(&["--check_frequency", "--freq_pop", "AF"]);
        let mut later = known("rs2", "T/C", &[("AF", "C:0.2")], &[("C", "C")]);
        later.start = 101;
        later.end = 101;
        let earlier = known("rs1", "T/C/G", &[("AF", "C:0.1")], &[("C", "C")]);
        let mut record = alleles_of("T", &["C", "G"], vec![]);
        record[0].colocated_variants = vec![later.clone()];
        record[1].colocated_variants = vec![earlier, later];
        f.filter_record(&mut record);
        assert_eq!(record[0].freqs, ["AF:C:0.2", "AF:G:NA"]);
    }
}
