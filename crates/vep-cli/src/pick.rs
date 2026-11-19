// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pick family and the summary modes, applied to one input record at a time
//! once every allele of the record is annotated: `--pick`, `--pick_allele`,
//! `--per_gene`, `--pick_allele_gene`, `--flag_pick`, `--flag_pick_allele`,
//! `--flag_pick_allele_gene` with `--pick_order`, then `--most_severe` and
//! `--summary`.
//!
//! Ensembl VEP filters the overlap alleles of a whole variation feature, that is
//! of a whole input record with every alternate allele it carries
//! (`OutputFactory.pm:578-676`), ranks them with
//! `pick_worst_VariationFeatureOverlapAllele` (`:703-811`), groups them per gene
//! (`:830-859`) and squashes what is left into one summary row (`:426-458`).
//! vep-rs annotates one allele at a time, so this module reassembles a record's
//! rows from its alleles in Ensembl VEP's order and applies the same steps.
//!
//! Perl citations name modules of ensembl-vep release/116.2 (`Bio/EnsEMBL/VEP/...`).

use smallvec::SmallVec;
use vep_core::consequence::{
    Consequence, ConsequenceList, FeatureType, SourceCache, TranscriptConsequence,
};
use vep_core::pick::{pick_by_categories, Rank};
use vep_core::variant::InputVariant;

use crate::row_filter::{self, RowFilters};

/// One `--pick_order` criterion, in Ensembl VEP's spelling (`Config.pm:601`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PickCriterion {
    ManeSelect,
    ManePlusClinical,
    Canonical,
    Appris,
    Tsl,
    Biotype,
    Ccds,
    Rank,
    Length,
    Ensembl,
    RefSeq,
}

impl PickCriterion {
    /// The default `--pick_order` (`Config.pm:311`).
    pub const DEFAULT_ORDER: [PickCriterion; 11] = [
        PickCriterion::ManeSelect,
        PickCriterion::ManePlusClinical,
        PickCriterion::Canonical,
        PickCriterion::Appris,
        PickCriterion::Tsl,
        PickCriterion::Biotype,
        PickCriterion::Ccds,
        PickCriterion::Rank,
        PickCriterion::Length,
        PickCriterion::Ensembl,
        PickCriterion::RefSeq,
    ];

    /// The criterion a `--pick_order` token names, `None` for a token outside
    /// `%VALID{pick_order}` (`Config.pm:601`).
    pub fn parse(token: &str) -> Option<Self> {
        Some(match token {
            "mane_select" => PickCriterion::ManeSelect,
            "mane_plus_clinical" => PickCriterion::ManePlusClinical,
            "canonical" => PickCriterion::Canonical,
            "appris" => PickCriterion::Appris,
            "tsl" => PickCriterion::Tsl,
            "biotype" => PickCriterion::Biotype,
            "ccds" => PickCriterion::Ccds,
            "rank" => PickCriterion::Rank,
            "length" => PickCriterion::Length,
            "ensembl" => PickCriterion::Ensembl,
            "refseq" => PickCriterion::RefSeq,
            _ => return None,
        })
    }

    /// The token Ensembl VEP names the criterion by.
    pub fn name(self) -> &'static str {
        match self {
            PickCriterion::ManeSelect => "mane_select",
            PickCriterion::ManePlusClinical => "mane_plus_clinical",
            PickCriterion::Canonical => "canonical",
            PickCriterion::Appris => "appris",
            PickCriterion::Tsl => "tsl",
            PickCriterion::Biotype => "biotype",
            PickCriterion::Ccds => "ccds",
            PickCriterion::Rank => "rank",
            PickCriterion::Length => "length",
            PickCriterion::Ensembl => "ensembl",
            PickCriterion::RefSeq => "refseq",
        }
    }
}

/// The `--pick_order` list as criteria: the default order when the flag is
/// absent, else its tokens in the given order. A token outside Ensembl VEP's
/// list is the error `check_config` raises for it (`Config.pm:945-957`), with
/// Perl's own `ERROR: ` prefix stripped.
pub fn pick_order(tokens: Option<&[String]>) -> Result<Vec<PickCriterion>, String> {
    let Some(tokens) = tokens else {
        return Ok(PickCriterion::DEFAULT_ORDER.to_vec());
    };
    tokens
        .iter()
        .map(|t| {
            PickCriterion::parse(t).ok_or_else(|| {
                format!(
                    "\"{t}\" is not a valid value for --pick_order\nValid values: {}",
                    PickCriterion::DEFAULT_ORDER
                        .iter()
                        .map(|c| c.name())
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            })
        })
        .collect()
}

/// The row filters, the pick family and the summary flags of a run, each a
/// command-line flag, and the criteria `pick_worst` compares in order.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// `--coding_only` and `--no_intergenic`, applied before the pick family.
    pub row_filters: RowFilters,
    pub pick: bool,
    pub pick_allele: bool,
    pub per_gene: bool,
    pub pick_allele_gene: bool,
    pub most_severe: bool,
    pub summary: bool,
    pub flag_pick: bool,
    pub flag_pick_allele: bool,
    pub flag_pick_allele_gene: bool,
    /// The `--pick_order` criteria in comparison order; a shorter list compares
    /// its criteria and then falls through to the first row.
    pub pick_order: Vec<PickCriterion>,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            row_filters: RowFilters::default(),
            pick: false,
            pick_allele: false,
            per_gene: false,
            pick_allele_gene: false,
            most_severe: false,
            summary: false,
            flag_pick: false,
            flag_pick_allele: false,
            flag_pick_allele_gene: false,
            pick_order: PickCriterion::DEFAULT_ORDER.to_vec(),
        }
    }
}

impl FilterConfig {
    /// Whether any record-level step, a row filter included, changes a record.
    pub fn any_active(&self) -> bool {
        self.row_filters.any_active()
            || self.pick
            || self.pick_allele
            || self.per_gene
            || self.pick_allele_gene
            || self.most_severe
            || self.summary
            || self.flag_pick
            || self.flag_pick_allele
            || self.flag_pick_allele_gene
    }
}

/// Applies the record-level steps to every record of an annotated batch laid
/// out in record order (every allele of a record adjacent to the others, as the
/// pipeline lays them out). A warning a record draws goes to the run's sink.
pub fn apply_to_batch(batch: &mut [InputVariant], config: &FilterConfig) {
    apply_to_batch_reporting(batch, config, &mut |_, text| {
        crate::warnings::global().warn(&text)
    });
}

/// [`apply_to_batch`] with the warning a record draws handed to `warn` with the
/// index in `batch` of the record's first variant, for the caller to report in
/// record order.
pub(crate) fn apply_to_batch_reporting(
    batch: &mut [InputVariant],
    config: &FilterConfig,
    warn: &mut dyn FnMut(usize, String),
) {
    if !config.any_active() {
        return;
    }
    let mut first = 0usize;
    for record in batch.chunk_by_mut(InputVariant::same_record) {
        apply_to_record_reporting(record, config, &mut |text| warn(first, text));
        first += record.len();
    }
}

/// Applies the record-level steps to the alleles of one input record in Ensembl
/// VEP's order (`OutputFactory.pm:477-516`, then `:387`): the row filters
/// (`:482`, `:490-502`), the pick family over the rows they leave (`:514`),
/// then the summary of the rows left. A warning the record draws goes to the
/// run's sink.
pub fn apply_to_record(record: &mut [InputVariant], config: &FilterConfig) {
    apply_to_record_reporting(record, config, &mut |text| {
        crate::warnings::global().warn(&text)
    });
}

/// [`apply_to_record`] with the warning the record draws handed to `warn`.
pub(crate) fn apply_to_record_reporting(
    record: &mut [InputVariant],
    config: &FilterConfig,
    warn: &mut dyn FnMut(String),
) {
    if !config.any_active() {
        return;
    }
    row_filter::apply_to_record(record, config.row_filters);
    pick_family(record, config);
    if config.most_severe || config.summary {
        // `summary_only` tests `summary` before `most_severe` (`:436`, `:446`), so
        // the pair `--safe` lets through lists every term.
        summarise(record, config.most_severe && !config.summary, warn);
    }
}

/// One row of a record as Ensembl VEP's `VariationFeatureOverlapAllele` list
/// holds it: an allele against a feature, a transcript, a regulatory feature or
/// a motif feature alike (`Transcript`, after the consequence list that holds
/// all three kinds), or an allele's intergenic row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Row {
    Transcript { allele: usize, tc: usize },
    Intergenic { allele: usize },
}

impl Row {
    fn allele(self) -> usize {
        match self {
            Row::Transcript { allele, .. } | Row::Intergenic { allele } => allele,
        }
    }
}

/// The record's rows in Ensembl VEP's order (`VariationFeature.pm:855-864` of
/// ensembl-variation): the transcript rows, then the regulatory rows, then the
/// motif rows, each kind feature-major in the string order of the keys its
/// variations are stored under (`:666`, `:710`, `:750` sort them), allele-minor
/// (each feature lists its reference allele's row under `--process_ref_homs`
/// and then its alternate allele); then the intergenic row of every allele that
/// overlaps no transcript, whether it carries regulatory and motif rows or
/// none. A feature's key is a transcript's stable id with `_` and its database
/// id, a regulatory feature's stable id or a motif's database id
/// ([`TranscriptConsequence::ensembl_vep_key`]); the rejoin of a `--minimal`
/// split stores every allele after the first under `<its allele string>_<key>`
/// (`OutputFactory.pm:2316`), so such an allele's feature sorts by its allele
/// string against the first allele's bare key, before it when the string starts
/// with a byte below the key's first. A structural variant's regulatory and
/// motif rows stand in the order the features were met instead
/// ([`TranscriptConsequence::structural_feature_order`]):
/// `pick_worst_VariationFeatureOverlapAllele` ranks the overlap list
/// `get_all_StructuralVariationOverlaps` returns as pushed
/// (`OutputFactory.pm:556-560`, `StructuralVariationFeature.pm:495-509`), so
/// a full tie falls to the first feature met. An allele whose chromosome has
/// no transcripts has no row.
fn record_rows(record: &[InputVariant]) -> Vec<Row> {
    let mut feature_rows: Vec<(usize, usize)> = Vec::new();
    for (allele, variant) in record.iter().enumerate() {
        for tc in 0..variant.transcript_consequences.len() {
            feature_rows.push((allele, tc));
        }
    }
    let feature = |&(allele, tc): &(usize, usize)| &record[allele].transcript_consequences[tc];
    let feature_key = |row: &(usize, usize)| {
        record[row.0]
            .rejoin_key_prefix
            .as_deref()
            .unwrap_or("")
            .bytes()
            .chain(feature(row).ensembl_vep_key())
    };
    let same_kind_order = |a: &(usize, usize), b: &(usize, usize)| {
        if record[a.0].is_structural && feature(a).feature_type != FeatureType::Transcript {
            feature(a).structural_feature_order(feature(b))
        } else {
            feature_key(a).cmp(feature_key(b))
        }
    };
    feature_rows.sort_by(|a, b| {
        feature(a)
            .feature_type
            .cmp(&feature(b).feature_type)
            .then_with(|| same_kind_order(a, b))
            .then_with(|| {
                record[b.0]
                    .reference_allele_row
                    .cmp(&record[a.0].reference_allele_row)
            })
            .then_with(|| a.0.cmp(&b.0))
    });
    let mut rows: Vec<Row> = feature_rows
        .into_iter()
        .map(|(allele, tc)| Row::Transcript { allele, tc })
        .collect();
    for (allele, variant) in record.iter().enumerate() {
        if variant.intergenic_overlap
            || (variant.transcript_consequences.is_empty()
                && variant.most_severe_consequence.is_some())
        {
            rows.push(Row::Intergenic { allele });
        }
    }
    rows
}

/// The consequence terms of an allele's intergenic row: `intergenic_variant`
/// beside the regulatory and motif rows of an allele that overlaps no
/// transcript ([`InputVariant::intergenic_overlap`]), else the verdict of an
/// allele with no feature row.
fn intergenic_row_terms(variant: &InputVariant) -> &[Consequence] {
    if variant.intergenic_overlap {
        &[Consequence::IntergenicVariant]
    } else {
        variant.most_severe_consequence.as_slice()
    }
}

/// The values `pick_worst_VariationFeatureOverlapAllele` ranks one row by, each
/// lower when better (`OutputFactory.pm:712-768`): the defaults for a row that is
/// not against a transcript, the transcript's attributes otherwise. Every value
/// is a number but the TSL, the attribute's digit run as text (`:753`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PickInfo<'a> {
    mane_select: i64,
    mane_plus_clinical: i64,
    canonical: i64,
    appris: i64,
    tsl: Rank<'a>,
    biotype: i64,
    ccds: i64,
    rank: i64,
    length: i64,
    ensembl: i64,
    refseq: i64,
}

impl<'a> PickInfo<'a> {
    /// The defaults of `:712-727`, with the rank of the row's consequences
    /// (`:780-786`; 1000 for a row with none).
    fn defaults(consequences: &[Consequence]) -> Self {
        PickInfo {
            mane_select: 1,
            mane_plus_clinical: 1,
            canonical: 1,
            appris: 100,
            tsl: Rank::Number(100),
            biotype: 1,
            ccds: 1,
            rank: consequences
                .iter()
                .map(|c| i64::from(c.rank()))
                .min()
                .unwrap_or(1000),
            length: 0,
            ensembl: 1,
            refseq: 1,
        }
    }

    /// A transcript row's values (`:729-768`).
    fn of_transcript(tc: &'a TranscriptConsequence) -> Self {
        let mut info = PickInfo::defaults(&tc.consequences);
        info.mane_select = if tc.mane_select().is_some() { 0 } else { 1 };
        info.mane_plus_clinical = if tc.mane_plus_clinical().is_some() {
            0
        } else {
            1
        };
        info.canonical = if tc.canonical { 0 } else { 1 };
        info.biotype = if tc.biotype() == Some("protein_coding") {
            0
        } else {
            1
        };
        // `$tr->{_ccds} && $tr->{_ccds} ne '-'` (`:740`).
        info.ccds = match tc.ccds() {
            Some(s) if !s.is_empty() && s != "-" => 0,
            _ => 1,
        };
        // `$info->{lc($tr->{_source_cache})} = 0` (`:741`): the source's name,
        // lower-cased, is the criterion it scores, so a custom source counts
        // only when it is named `ensembl` or `refseq` in any case.
        if let Some(source) = tc.source_cache().map(SourceCache::column_text) {
            if source.eq_ignore_ascii_case("ensembl") {
                info.ensembl = 0;
            } else if source.eq_ignore_ascii_case("refseq") {
                info.refseq = 0;
            }
        }
        // Longer is better, so the length is negated (`:743-748`): the
        // translateable sequence of a translated transcript, else the spliced
        // length.
        info.length = -i64::from(tc.cds_length().unwrap_or_else(|| tc.cdna_length()));
        // `m/tsl(\d+)/` then `$1 if $1` (`:751-755`): the level's digit run
        // as text, which the loader leaves out for a zero.
        info.tsl = Rank::tsl(tc.tsl());
        info.appris = i64::from(appris_rank(tc.appris()));
        info
    }

    fn value(&self, criterion: PickCriterion) -> Rank<'a> {
        match criterion {
            PickCriterion::ManeSelect => Rank::Number(self.mane_select),
            PickCriterion::ManePlusClinical => Rank::Number(self.mane_plus_clinical),
            PickCriterion::Canonical => Rank::Number(self.canonical),
            PickCriterion::Appris => Rank::Number(self.appris),
            PickCriterion::Tsl => self.tsl,
            PickCriterion::Biotype => Rank::Number(self.biotype),
            PickCriterion::Ccds => Rank::Number(self.ccds),
            PickCriterion::Rank => Rank::Number(self.rank),
            PickCriterion::Length => Rank::Number(self.length),
            PickCriterion::Ensembl => Rank::Number(self.ensembl),
            PickCriterion::RefSeq => Rank::Number(self.refseq),
        }
    }
}

/// The APPRIS attribute as `pick_worst_VariationFeatureOverlapAllele` scores it
/// (`OutputFactory.pm:758-768`): `m/([A-Za-z]).+(\d+)/` reads the first letter
/// and, `.+` being greedy, the last digit alone, which must come at least two
/// characters after the letter; `principal1` scores 1 and `alternative1` 11 (the
/// grade plus 10 when the letter is a lower-case `a`); a score of 0 and anything
/// the pattern does not match keep the default 100.
fn appris_rank(appris: Option<&str>) -> u32 {
    let Some(s) = appris else { return 100 };
    let bytes = s.as_bytes();
    let Some(letter_pos) = bytes.iter().position(|b| b.is_ascii_alphabetic()) else {
        return 100;
    };
    let Some(digit_pos) = bytes.iter().rposition(|b| b.is_ascii_digit()) else {
        return 100;
    };
    if digit_pos < letter_pos + 2 {
        return 100;
    }
    // `$grade += 10 if substr($type, 0, 1) eq 'a'; $info->{appris} = $grade if $grade;`
    let mut grade = u32::from(bytes[digit_pos] - b'0');
    if bytes[letter_pos] == b'a' {
        grade += 10;
    }
    if grade == 0 {
        return 100;
    }
    grade
}

/// The row `pick_worst_VariationFeatureOverlapAllele` returns from
/// `candidates` (`OutputFactory.pm:703-811`): the per-category stable sort and
/// shrink of `:789-806`, on a full tie the first candidate (`:808`).
fn pick_worst(candidates: &[usize], infos: &[PickInfo], order: &[PickCriterion]) -> Option<usize> {
    pick_by_categories(candidates, order, |candidate, criterion| {
        infos[candidate].value(criterion)
    })
}

/// The candidates grouped by `key`, groups in order of first appearance and
/// each group in candidate order.
fn group_by<K: PartialEq>(candidates: &[usize], key: impl Fn(usize) -> K) -> Vec<Vec<usize>> {
    let mut keys: Vec<K> = Vec::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for &c in candidates {
        let k = key(c);
        match keys.iter().position(|existing| *existing == k) {
            Some(i) => groups[i].push(c),
            None => {
                keys.push(k);
                groups.push(vec![c]);
            }
        }
    }
    groups
}

/// `pick_VariationFeatureOverlapAllele_per_gene` (`OutputFactory.pm:830-859`):
/// every row not against a transcript passes through, and each gene keeps the
/// worst of its transcript rows.
fn pick_per_gene(
    candidates: &[usize],
    rows: &[Row],
    record: &[InputVariant],
    infos: &[PickInfo],
    order: &[PickCriterion],
) -> Vec<usize> {
    let mut picked: Vec<usize> = Vec::new();
    let mut transcript_rows: Vec<usize> = Vec::new();
    for &c in candidates {
        match rows[c] {
            Row::Transcript { allele, tc }
                if record[allele].transcript_consequences[tc].feature_type
                    == FeatureType::Transcript =>
            {
                transcript_rows.push(c)
            }
            _ => picked.push(c),
        }
    }
    let gene = |c: usize| match rows[c] {
        Row::Transcript { allele, tc } => record[allele].transcript_consequences[tc].gene_id(),
        Row::Intergenic { .. } => unreachable!("only transcript rows are grouped by gene"),
    };
    for group in group_by(&transcript_rows, gene) {
        picked.extend(pick_worst(&group, infos, order));
    }
    picked
}

/// `filter_VariationFeatureOverlapAlleles` (`OutputFactory.pm:578-627`) for a
/// sequence variant record and `filter_StructuralVariationOverlapAlleles`
/// (`:653-673`) for a structural one, where the per-allele variants collapse to
/// their record-level forms. The removing variants keep the picked rows alone;
/// the flagging variants mark them.
fn pick_family(record: &mut [InputVariant], config: &FilterConfig) {
    // The variant in force is the first set flag of the chain at `:585-626`:
    // whether the rows are ranked per allele or over the whole record, whether
    // each gene keeps a row, and whether the picks are flagged or kept alone.
    let (per_allele, per_gene, flag) = if config.pick {
        (false, false, false)
    } else if config.pick_allele {
        (true, false, false)
    } else if config.per_gene {
        (false, true, false)
    } else if config.pick_allele_gene {
        (true, true, false)
    } else if config.flag_pick {
        (false, false, true)
    } else if config.flag_pick_allele {
        (true, false, true)
    } else if config.flag_pick_allele_gene {
        (true, true, true)
    } else {
        return;
    };

    let rows = record_rows(record);
    if rows.is_empty() {
        return;
    }
    let infos: Vec<PickInfo> = rows
        .iter()
        .map(|row| match *row {
            Row::Transcript { allele, tc } => {
                PickInfo::of_transcript(&record[allele].transcript_consequences[tc])
            }
            Row::Intergenic { allele } => PickInfo::defaults(intergenic_row_terms(&record[allele])),
        })
        .collect();
    let order = config.pick_order.as_slice();
    let all: Vec<usize> = (0..rows.len()).collect();
    let structural = record.iter().any(|v| v.is_structural);
    // `%by_allele` groups on `variation_feature_seq`, the Allele column (`:592`);
    // a structural record ranks every allele's rows together (`:653-673`).
    let groups: Vec<Vec<usize>> = if per_allele && !structural {
        group_by(&all, |c| record[rows[c].allele()].display_allele())
    } else {
        vec![all]
    };
    let picked: Vec<usize> = if per_gene {
        groups
            .iter()
            .flat_map(|group| pick_per_gene(group, &rows, record, &infos, order))
            .collect()
    } else {
        groups
            .iter()
            .filter_map(|group| pick_worst(group, &infos, order))
            .collect()
    };

    if flag {
        for &c in &picked {
            match rows[c] {
                Row::Transcript { allele, tc } => {
                    record[allele].transcript_consequences[tc].picked = true
                }
                Row::Intergenic { allele } => record[allele].intergenic_pick = true,
            }
        }
    } else {
        retain_rows(record, &rows, &picked);
        // The output walks the list `pick_VariationFeatureOverlapAllele_per_gene`
        // returns (`:838-856`): the rows not against a transcript first, then the
        // gene picks; the writers read the mark to list the record's rows so.
        if per_gene {
            for variant in record.iter_mut() {
                variant.per_gene_pick = true;
            }
        }
    }
}

/// Keeps the rows at `picked` and removes every other row of the record: an
/// allele keeps its picked transcript consequences in their order, and an allele
/// left with no row loses its consequence verdict so no writer prints one for it.
fn retain_rows(record: &mut [InputVariant], rows: &[Row], picked: &[usize]) {
    for (allele, variant) in record.iter_mut().enumerate() {
        let keeps_intergenic = picked
            .iter()
            .any(|&c| rows[c] == Row::Intergenic { allele });
        variant.intergenic_overlap &= keeps_intergenic;
        if variant.transcript_consequences.is_empty() {
            if !keeps_intergenic {
                variant.most_severe_consequence = None;
            }
            continue;
        }
        let kept_tcs: Vec<usize> = picked
            .iter()
            .filter_map(|&c| match rows[c] {
                Row::Transcript { allele: a, tc } if a == allele => Some(tc),
                _ => None,
            })
            .collect();
        let mut index = 0usize;
        variant.transcript_consequences.retain(|_| {
            let keep = kept_tcs.contains(&index);
            index += 1;
            keep
        });
        row_filter::settle_verdict(variant);
    }
}

/// `summary_only` (`OutputFactory.pm:426-458`): the consequence terms of the rows
/// left on the record, sorted by rank, each once (`--summary`) or the first alone
/// (`--most_severe`), carried by the record's first allele as its one row while
/// every other allele keeps none. A record with no row keeps its one summary row
/// with no term and draws `Unable to assign consequence type`, handed to `warn`
/// (`:454-458` warn and return the record's hash with no `Consequence`): a
/// reference-homozygous individual's copy, or a non-variant record.
fn summarise(record: &mut [InputVariant], most_severe_only: bool, warn: &mut dyn FnMut(String)) {
    let mut terms: Vec<Consequence> = Vec::new();
    for row in record_rows(record) {
        let row_terms: &[Consequence] = match row {
            Row::Transcript { allele, tc } => {
                &record[allele].transcript_consequences[tc].consequences
            }
            Row::Intergenic { allele } => intergenic_row_terms(&record[allele]),
        };
        terms.extend_from_slice(row_terms);
    }
    terms.sort_by_key(|c| c.rank());
    terms.dedup();
    if most_severe_only {
        terms.truncate(1);
    }
    if terms.is_empty() {
        warn("Unable to assign consequence type".to_string());
    }
    for variant in record.iter_mut() {
        variant.transcript_consequences.clear();
        variant.intergenic_overlap = false;
        variant.most_severe_consequence = None;
        variant.intergenic_pick = false;
    }
    let Some(leader) = record.first_mut() else {
        return;
    };
    let summary: ConsequenceList = SmallVec::from_vec(terms);
    leader.most_severe_consequence = summary.first().copied();
    leader.summary_terms = Some(summary);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runner::{edit_transcript, set_lengths};
    use smallvec::smallvec;
    use std::sync::Arc;
    use vep_core::consequence::{Consequence, FeatureType, Impact};
    use vep_core::transcript::Transcript;
    use vep_core::variant::InputVariant;

    fn tc(
        transcript_id: &str,
        gene_id: &str,
        consequences: ConsequenceList,
    ) -> TranscriptConsequence {
        let impact = consequences
            .first()
            .map(|c| c.impact())
            .unwrap_or(Impact::MODIFIER);
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: transcript_id.into(),
                gene_stable_id: gene_id.into(),
                biotype: "protein_coding".into(),
                ..Transcript::default()
            })),
            consequences,
            impact,
            feature_type: FeatureType::Transcript,
            strand: 1,
            ..Default::default()
        }
    }

    /// A regulatory or motif row of `feature_type` carrying the one term
    /// `consequence`.
    fn feature(
        id: &str,
        feature_type: FeatureType,
        consequence: Consequence,
    ) -> TranscriptConsequence {
        TranscriptConsequence {
            feature_id: Some(id.into()),
            feature_type,
            consequences: smallvec![consequence],
            ..Default::default()
        }
    }

    fn allele(alt: &[u8], index: usize, tcs: Vec<TranscriptConsequence>) -> InputVariant {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), alt.to_vec());
        v.allele_index = index;
        v.input_record = 1;
        v.transcript_consequences = tcs;
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        v
    }

    fn config(set: impl Fn(&mut FilterConfig)) -> FilterConfig {
        let mut c = FilterConfig::default();
        set(&mut c);
        c
    }

    /// `(allele index, transcript id)` of every row left, in writer order.
    fn kept(record: &[InputVariant]) -> Vec<(usize, String)> {
        let mut out: Vec<(usize, String)> = Vec::new();
        for v in record {
            for t in &v.transcript_consequences {
                out.push((v.allele_index, t.transcript_id().to_string()));
            }
            if v.transcript_consequences.is_empty() && v.most_severe_consequence.is_some() {
                out.push((v.allele_index, "-".to_string()));
            }
        }
        out
    }

    /// `(allele index, transcript id)` of every row flagged PICK.
    fn flagged(record: &[InputVariant]) -> Vec<(usize, String)> {
        let mut out: Vec<(usize, String)> = Vec::new();
        for v in record {
            for t in &v.transcript_consequences {
                if t.picked {
                    out.push((v.allele_index, t.transcript_id().to_string()));
                }
            }
            if v.intergenic_pick {
                out.push((v.allele_index, "-".to_string()));
            }
        }
        out
    }

    /// A bi-allelic record over two genes: ENST1 (gene G1, canonical) and ENST2
    /// (gene G2), the G allele missense on ENST1 and the T allele synonymous.
    fn two_gene_record() -> Vec<InputVariant> {
        let mut t1 = tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]);
        t1.canonical = true;
        let t2 = tc("ENST2", "G2", smallvec![Consequence::IntronVariant]);
        let mut t1_t = tc("ENST1", "G1", smallvec![Consequence::SynonymousVariant]);
        t1_t.canonical = true;
        let t2_t = tc("ENST2", "G2", smallvec![Consequence::IntronVariant]);
        vec![
            allele(b"G", 0, vec![t1, t2]),
            allele(b"T", 1, vec![t1_t, t2_t]),
        ]
    }

    fn structural_two_gene_record() -> Vec<InputVariant> {
        let mut record = two_gene_record();
        for v in record.iter_mut() {
            v.is_structural = true;
        }
        record
    }

    /// `--pick` keeps one row for the whole record, whatever the allele count
    /// (`filter_VariationFeatureOverlapAlleles` ranks every allele's overlap
    /// alleles together, `OutputFactory.pm:585-588`): the canonical ENST1 of the
    /// first allele, whose missense term outranks the second allele's synonymous
    /// one on the same transcript.
    #[test]
    fn pick_keeps_one_row_per_record_across_alts() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENST1".to_string())]);
        assert_eq!(
            record[1].most_severe_consequence, None,
            "an allele left without a row prints none"
        );
        assert_eq!(
            record[0].most_severe_consequence,
            Some(Consequence::MissenseVariant)
        );
    }

    /// Each individual's copy of a record is a feature of its own in Ensembl VEP
    /// (`Parser/VCF.pm:629-664`), picked apart from the record's other copies
    /// (`filter_VariationFeatureOverlapAlleles` runs per feature,
    /// `OutputFactory.pm:514`): `--pick` over a batch of two copies keeps one row
    /// per copy, and the reference row of `--process_ref_homs` ranks with the
    /// copy's other rows.
    #[test]
    fn pick_is_per_individual_variant() {
        let mut batch = two_gene_record();
        let mut second = two_gene_record();
        for v in batch.iter_mut() {
            v.individual = Some("S1".into());
        }
        for v in second.iter_mut() {
            v.individual = Some("S2".into());
        }
        batch.extend(second);
        assert!(!batch[1].same_record(&batch[2]) && batch[0].same_record(&batch[1]));
        apply_to_batch(&mut batch, &config(|c| c.pick = true));
        assert_eq!(
            kept(&batch),
            vec![(0, "ENST1".to_string()), (0, "ENST1".to_string())]
        );
        assert_eq!(batch[0].individual.as_deref(), Some("S1"));
        assert_eq!(batch[2].individual.as_deref(), Some("S2"));

        let mut record = two_gene_record();
        let mut reference = record[0].reference_allele_row();
        let mut t1 = tc("ENST1", "G1", smallvec![Consequence::SynonymousVariant]);
        t1.canonical = true;
        reference.transcript_consequences = vec![t1];
        record.insert(0, reference);
        apply_to_batch(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENST1".to_string())]);
        assert!(
            record[0].transcript_consequences.is_empty(),
            "the reference row's synonymous ENST1 loses to the G allele's missense ENST1"
        );
    }

    /// `--pick_allele` ranks each allele's rows on their own (`:590-595`, grouped
    /// on the Allele column), so a bi-allelic record keeps one row per allele.
    #[test]
    fn pick_allele_keeps_one_row_per_alt() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.pick_allele = true));
        assert_eq!(
            kept(&record),
            vec![(0, "ENST1".to_string()), (1, "ENST1".to_string())]
        );
    }

    /// On a structural record `--pick_allele` is `--pick`, `--flag_pick_allele` is
    /// `--flag_pick`, `--pick_allele_gene` is `--per_gene` and
    /// `--flag_pick_allele_gene` flags one row per gene over the whole record
    /// (`filter_StructuralVariationOverlapAlleles`, `:653-673`).
    #[test]
    fn structural_record_collapses_the_per_allele_flags() {
        let mut record = structural_two_gene_record();
        apply_to_record(&mut record, &config(|c| c.pick_allele = true));
        assert_eq!(kept(&record), vec![(0, "ENST1".to_string())]);

        let mut record = structural_two_gene_record();
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele = true));
        assert_eq!(flagged(&record), vec![(0, "ENST1".to_string())]);
        assert_eq!(kept(&record).len(), 4, "flagging removes nothing");

        let per_gene_rows = vec![(0, "ENST1".to_string()), (0, "ENST2".to_string())];
        let mut record = structural_two_gene_record();
        apply_to_record(&mut record, &config(|c| c.pick_allele_gene = true));
        assert_eq!(kept(&record), per_gene_rows);

        let mut record = structural_two_gene_record();
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele_gene = true));
        assert_eq!(flagged(&record), per_gene_rows);
        assert_eq!(kept(&record).len(), 4);
    }

    /// `--per_gene` keeps the worst row of each gene over the whole record
    /// (`pick_VariationFeatureOverlapAllele_per_gene`, `:830-859`, keyed on the
    /// gene stable id at `:851`): ENST1 of the first allele for G1 and ENST2 of
    /// the first allele for G2, each a full tie resolved to the first row.
    #[test]
    fn per_gene_keeps_one_row_per_gene_per_record() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.per_gene = true));
        assert_eq!(
            kept(&record),
            vec![(0, "ENST1".to_string()), (0, "ENST2".to_string())]
        );
    }

    /// Rows not against a transcript pass through the per-gene picks untouched
    /// (`:838-844`): the regulatory row, the motif row and the intergenic row
    /// survive, and only the two G1 transcripts are reduced to one.
    #[test]
    fn per_gene_passes_regulatory_and_intergenic_rows_through() {
        let build = |with_transcripts: bool| {
            let mut rows = vec![
                feature(
                    "ENSR1",
                    FeatureType::RegulatoryFeature,
                    Consequence::RegulatoryRegionVariant,
                ),
                feature(
                    "ENSM1",
                    FeatureType::MotifFeature,
                    Consequence::TfBindingSiteVariant,
                ),
            ];
            if with_transcripts {
                rows.push(tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]));
                rows.push(tc("ENST2", "G1", smallvec![Consequence::SynonymousVariant]));
            }
            let mut v = allele(b"G", 0, rows);
            v.intergenic_overlap = !with_transcripts;
            vec![v]
        };
        for set in [
            |c: &mut FilterConfig| c.per_gene = true,
            |c: &mut FilterConfig| c.pick_allele_gene = true,
        ] {
            let mut record = build(false);
            apply_to_record(&mut record, &config(set));
            assert_eq!(record[0].transcript_consequences.len(), 2);
            assert!(record[0].intergenic_overlap);

            let mut record = build(true);
            apply_to_record(&mut record, &config(set));
            let ids: Vec<&str> = record[0]
                .transcript_consequences
                .iter()
                .map(|t| t.transcript_id())
                .collect();
            assert_eq!(ids, ["ENSR1", "ENSM1", "ENST1"]);
        }

        let mut record = build(false);
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele_gene = true));
        assert!(record[0].transcript_consequences.iter().all(|t| t.picked));
        assert!(record[0].intergenic_pick);
    }

    /// `(Feature_type, Feature)` of the record's rows in the order the writers
    /// list them.
    fn written(record: &[InputVariant]) -> Vec<(String, String)> {
        let refs: Vec<&InputVariant> = record.iter().collect();
        vep_io::output::fields::record_rows(&refs)
            .into_iter()
            .map(|(_, tc)| match tc {
                Some(tc) => (tc.feature_type.to_string(), tc.transcript_id().to_string()),
                None => ("-".to_string(), "-".to_string()),
            })
            .collect()
    }

    /// `pick_VariationFeatureOverlapAllele_per_gene` returns the alleles that are
    /// not against a transcript first and one transcript allele per gene after
    /// them (`:838-856`), and the output walks that list, so a `--per_gene` or
    /// `--pick_allele_gene` run writes a record's regulatory row, motif row and
    /// intergenic row before its gene picks, where the unpicked record lists its
    /// transcript rows first (`VariationFeature.pm:855-864`); a flagging run
    /// removes nothing and keeps the overlap-list order.
    #[test]
    fn per_gene_passes_regulatory_motif_and_intergenic_rows_through() {
        let build = || {
            let mut v = allele(
                b"G",
                0,
                vec![
                    tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]),
                    tc("ENST2", "G1", smallvec![Consequence::SynonymousVariant]),
                    tc("ENST3", "G2", smallvec![Consequence::IntronVariant]),
                    feature(
                        "ENSR1",
                        FeatureType::RegulatoryFeature,
                        Consequence::RegulatoryRegionVariant,
                    ),
                    feature(
                        "ENSM1",
                        FeatureType::MotifFeature,
                        Consequence::TfBindingSiteVariant,
                    ),
                ],
            );
            v.most_severe_consequence = Some(Consequence::MissenseVariant);
            vec![v]
        };
        let row = |kind: &str, id: &str| (kind.to_string(), id.to_string());
        let record = build();
        assert_eq!(
            written(&record),
            [
                row("Transcript", "ENST1"),
                row("Transcript", "ENST2"),
                row("Transcript", "ENST3"),
                row("RegulatoryFeature", "ENSR1"),
                row("MotifFeature", "ENSM1"),
            ]
        );
        for set in [
            |c: &mut FilterConfig| c.per_gene = true,
            |c: &mut FilterConfig| c.pick_allele_gene = true,
        ] {
            let mut record = build();
            apply_to_record(&mut record, &config(set));
            assert!(record[0].per_gene_pick);
            assert_eq!(
                written(&record),
                [
                    row("RegulatoryFeature", "ENSR1"),
                    row("MotifFeature", "ENSM1"),
                    row("Transcript", "ENST1"),
                    row("Transcript", "ENST3"),
                ]
            );
        }
        let mut record = build();
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele_gene = true));
        assert!(!record[0].per_gene_pick);
        assert_eq!(written(&record).len(), 5);
        assert_eq!(written(&record)[0], row("Transcript", "ENST1"));

        let mut record = build();
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert!(!record[0].per_gene_pick);
        assert_eq!(written(&record), [row("Transcript", "ENST1")]);

        let mut no_transcript = build();
        no_transcript[0].transcript_consequences.drain(..3);
        no_transcript[0].intergenic_overlap = true;
        apply_to_record(&mut no_transcript, &config(|c| c.per_gene = true));
        assert_eq!(
            written(&no_transcript),
            [
                row("RegulatoryFeature", "ENSR1"),
                row("MotifFeature", "ENSM1"),
                row("-", "-"),
            ]
        );
    }

    /// `--flag_pick_allele_gene` flags everything the per-gene pick returns
    /// (`:623-626`), the pass-through regulatory and motif rows included, beside
    /// one transcript row per gene; `--flag_pick` and `--flag_pick_allele` flag the
    /// one worst row of the record, which a canonical transcript's row wins.
    #[test]
    fn flag_pick_allele_gene_flags_regulatory_rows() {
        let build = || {
            let mut canonical = tc("ENST1", "G1", smallvec![Consequence::IntronVariant]);
            canonical.canonical = true;
            vec![allele(
                b"G",
                0,
                vec![
                    canonical,
                    tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]),
                    feature(
                        "ENSR1",
                        FeatureType::RegulatoryFeature,
                        Consequence::RegulatoryRegionVariant,
                    ),
                    feature(
                        "ENSM1",
                        FeatureType::MotifFeature,
                        Consequence::TfBindingSiteVariant,
                    ),
                ],
            )]
        };
        let mut record = build();
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele_gene = true));
        assert_eq!(
            flagged(&record),
            [
                (0, "ENST1".to_string()),
                (0, "ENSR1".to_string()),
                (0, "ENSM1".to_string())
            ]
        );
        assert_eq!(kept(&record).len(), 4, "flagging removes nothing");
        for set in [
            |c: &mut FilterConfig| c.flag_pick = true,
            |c: &mut FilterConfig| c.flag_pick_allele = true,
        ] {
            let mut record = build();
            apply_to_record(&mut record, &config(set));
            assert_eq!(flagged(&record), [(0, "ENST1".to_string())]);
        }
    }

    /// `--pick_allele_gene` groups on the allele first (`:601-605`), so each
    /// allele keeps one row per gene.
    #[test]
    fn pick_allele_gene_keeps_one_row_per_alt_and_gene() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.pick_allele_gene = true));
        assert_eq!(
            kept(&record),
            vec![
                (0, "ENST1".to_string()),
                (0, "ENST2".to_string()),
                (1, "ENST1".to_string()),
                (1, "ENST2".to_string())
            ]
        );
    }

    /// `--flag_pick` marks the one row `--pick` would keep and removes nothing
    /// (`:606-611`).
    #[test]
    fn flag_pick_marks_one_row_per_record() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.flag_pick = true));
        assert_eq!(flagged(&record), vec![(0, "ENST1".to_string())]);
        assert_eq!(kept(&record).len(), 4);
    }

    /// `--flag_pick_allele` marks one row per allele (`:614-619`).
    #[test]
    fn flag_pick_allele_marks_one_row_per_alt() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele = true));
        assert_eq!(
            flagged(&record),
            vec![(0, "ENST1".to_string()), (1, "ENST1".to_string())]
        );
    }

    /// `--flag_pick_allele_gene` marks one row per allele and gene (`:622-626`).
    #[test]
    fn flag_pick_allele_gene_marks_one_row_per_alt_and_gene() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele_gene = true));
        assert_eq!(
            flagged(&record),
            vec![
                (0, "ENST1".to_string()),
                (0, "ENST2".to_string()),
                (1, "ENST1".to_string()),
                (1, "ENST2".to_string())
            ]
        );
    }

    /// A record overlapping no transcript has one intergenic row per allele; the
    /// flagging variants mark it as they mark a transcript row (the intergenic
    /// allele passes through `VariationFeatureOverlapAllele_to_output_hash`,
    /// `:1363`), and `--pick` keeps the first allele's alone.
    #[test]
    fn intergenic_rows_are_picked_and_flagged_like_transcript_rows() {
        let mut record = vec![allele(b"G", 0, vec![]), allele(b"T", 1, vec![])];
        apply_to_record(&mut record, &config(|c| c.flag_pick = true));
        assert_eq!(flagged(&record), vec![(0, "-".to_string())]);

        let mut record = vec![allele(b"G", 0, vec![]), allele(b"T", 1, vec![])];
        apply_to_record(&mut record, &config(|c| c.flag_pick_allele_gene = true));
        assert_eq!(
            flagged(&record),
            vec![(0, "-".to_string()), (1, "-".to_string())],
            "per-gene picking passes a non-transcript row through to the flagging map (:626, :836-844)"
        );

        let mut record = vec![allele(b"G", 0, vec![]), allele(b"T", 1, vec![])];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "-".to_string())]);
    }

    /// With `rank` first in `--pick_order` the most severe row wins over the
    /// canonical one (`:789-806` compares the categories in the given order).
    #[test]
    fn pick_order_rank_first_changes_the_pick() {
        let mut canonical = tc("ENST1", "G1", smallvec![Consequence::IntronVariant]);
        canonical.canonical = true;
        let severe = tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]);
        let mut record = vec![allele(b"G", 0, vec![canonical.clone(), severe.clone()])];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENST1".to_string())]);

        let mut record = vec![allele(b"G", 0, vec![canonical, severe])];
        let rank_first = config(|c| {
            c.pick = true;
            c.pick_order =
                pick_order(Some(&["rank".to_string(), "canonical".to_string()])).unwrap();
        });
        apply_to_record(&mut record, &rank_first);
        assert_eq!(kept(&record), vec![(0, "ENST2".to_string())]);
    }

    /// `length` is the negated translateable-sequence length of a translated
    /// transcript, else the negated spliced length (`:743-748`), so the longer
    /// CDS wins between coding transcripts and the longer cDNA between
    /// non-coding ones.
    #[test]
    fn length_prefers_the_longer_cds_then_cdna() {
        let mut short_cds = tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]);
        edit_transcript(&mut short_cds, |t| set_lengths(t, Some(300), 5000));
        let mut long_cds = tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]);
        edit_transcript(&mut long_cds, |t| set_lengths(t, Some(600), 1000));
        let mut record = vec![allele(b"G", 0, vec![short_cds, long_cds])];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(
            kept(&record),
            vec![(0, "ENST2".to_string())],
            "the longer CDS wins although its cDNA is shorter"
        );

        let mut short = tc("ENST1", "G1", smallvec![Consequence::IntronVariant]);
        edit_transcript(&mut short, |transcript| {
            transcript.biotype = "lncRNA".into();
            set_lengths(transcript, None, 800);
        });
        let mut long = tc("ENST2", "G1", smallvec![Consequence::IntronVariant]);
        edit_transcript(&mut long, |transcript| {
            transcript.biotype = "lncRNA".into();
            set_lengths(transcript, None, 2400);
        });
        let mut record = vec![allele(b"G", 0, vec![short, long])];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENST2".to_string())]);
    }

    /// On a single-source cache no transcript carries `_source_cache`, so the
    /// `ensembl` and `refseq` categories score 1 on every row (`:741` sets a
    /// category to 0 only when the key exists) and decide nothing.
    #[test]
    fn ensembl_and_refseq_tie_on_a_single_source_cache() {
        let a = tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]);
        let b = tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]);
        assert_eq!(a.source_cache(), None);
        let info = PickInfo::of_transcript(&a);
        assert_eq!((info.ensembl, info.refseq), (1, 1));
        let mut record = vec![allele(b"G", 0, vec![a, b])];
        let source_only = config(|c| {
            c.pick = true;
            c.pick_order = vec![PickCriterion::Ensembl, PickCriterion::RefSeq];
        });
        apply_to_record(&mut record, &source_only);
        assert_eq!(
            kept(&record),
            vec![(0, "ENST1".to_string())],
            "a full tie falls to the first row in stable-id order"
        );
        let mut refseq = tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]);
        edit_transcript(&mut refseq, |t| t.source_cache = Some(SourceCache::RefSeq));
        let info = PickInfo::of_transcript(&refseq);
        assert_eq!((info.ensembl, info.refseq), (1, 0));
    }

    /// The `tsl` criterion is the attribute's digit run as text (`:751-753`),
    /// sorted as a number (`<=>`, `:789`): `300` loses to `2` and to a row
    /// without a level (the default 100), and `00` beats `1`. The tie test
    /// reads the text (`eq`, `:796`): `007` sorts with `7` but is no tie of
    /// it, so it ends the tied prefix, and only the rows tied as text go on to
    /// the next criterion.
    #[test]
    fn tsl_ranks_the_digit_run_as_a_number_and_ties_it_as_text() {
        let with_tsl = |id: &str, level: Option<&str>, cdna: u32| {
            let mut row = tc(id, "G1", smallvec![Consequence::MissenseVariant]);
            edit_transcript(&mut row, |t| {
                t.tsl = level.map(str::to_string);
                set_lengths(t, None, cdna);
            });
            row
        };
        let tsl_then_length = config(|c| {
            c.pick = true;
            c.pick_order = vec![PickCriterion::Tsl, PickCriterion::Length];
        });
        let picked = |rows: Vec<TranscriptConsequence>| {
            let mut record = vec![allele(b"G", 0, rows)];
            apply_to_record(&mut record, &tsl_then_length);
            kept(&record)
        };
        let one = |id: &str| vec![(0, id.to_string())];
        assert_eq!(
            picked(vec![
                with_tsl("ENST1", Some("300"), 900),
                with_tsl("ENST2", Some("2"), 100),
            ]),
            one("ENST2")
        );
        assert_eq!(
            picked(vec![
                with_tsl("ENST1", Some("300"), 900),
                with_tsl("ENST2", None, 100),
            ]),
            one("ENST2")
        );
        assert_eq!(
            picked(vec![
                with_tsl("ENST1", Some("1"), 900),
                with_tsl("ENST2", Some("00"), 100),
            ]),
            one("ENST2")
        );
        // Rows are ranked in stable-id order: the `007` second in that order
        // ends the tie at the first `7`, which is the pick though the shortest;
        // third, it leaves the two `7`s to the length criterion.
        assert_eq!(
            picked(vec![
                with_tsl("ENST1", Some("7"), 100),
                with_tsl("ENST2", Some("007"), 900),
                with_tsl("ENST3", Some("7"), 500),
            ]),
            one("ENST1")
        );
        assert_eq!(
            picked(vec![
                with_tsl("ENST1", Some("7"), 100),
                with_tsl("ENST2", Some("7"), 500),
                with_tsl("ENST3", Some("007"), 900),
            ]),
            one("ENST2")
        );
    }

    /// Two rows equal on every criterion: the pick is the first in Ensembl VEP's
    /// overlap-allele order (`:808` returns the first of what the stable sorts
    /// left), the string order of the keys the transcript variations are stored
    /// under (`VariationFeature.pm:666` of ensembl-variation), each the stable id
    /// with `_` and the database id (`BaseVariationFeature.pm:310-320`). `_`
    /// sorts after every digit, so on the TAIR10 record `4 7204233 CA CT,CAAGG`
    /// under `--pick --pick_order rank`, where the eleven `AAGG` rows of
    /// AT4G12020 tie on `inframe_insertion`, Ensembl VEP 116.2 keeps
    /// `AT4G12020.10` (key `AT4G12020.10_118783`), not the bare-id-first
    /// `AT4G12020.1` (`AT4G12020.1_118785`); `--flag_pick` marks that same row.
    /// A transcript without a database id keys by its stable id alone.
    #[test]
    fn full_tie_falls_to_the_first_transcript_key_with_its_database_id() {
        let later = tc(
            "ENST00000400000",
            "G1",
            smallvec![Consequence::MissenseVariant],
        );
        let earlier = tc(
            "ENST00000300000",
            "G1",
            smallvec![Consequence::MissenseVariant],
        );
        let mut record = vec![allele(b"G", 0, vec![later, earlier])];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENST00000300000".to_string())]);

        let isoforms: [(&str, u64); 11] = [
            ("AT4G12020.1", 118785),
            ("AT4G12020.2", 118784),
            ("AT4G12020.3", 118779),
            ("AT4G12020.4", 118777),
            ("AT4G12020.5", 118782),
            ("AT4G12020.6", 118786),
            ("AT4G12020.7", 118776),
            ("AT4G12020.8", 118780),
            ("AT4G12020.9", 118781),
            ("AT4G12020.10", 118783),
            ("AT4G12020.11", 118778),
        ];
        let rows = |terms: &ConsequenceList| -> Vec<TranscriptConsequence> {
            isoforms
                .iter()
                .map(|&(id, db_id)| {
                    let mut row = tc(id, "AT4G12020", terms.clone());
                    edit_transcript(&mut row, |t| t.db_id = Some(db_id));
                    row
                })
                .collect()
        };
        let tair_record = || {
            vec![
                allele(
                    b"T",
                    0,
                    rows(&smallvec![
                        Consequence::SpliceRegionVariant,
                        Consequence::SynonymousVariant
                    ]),
                ),
                allele(
                    b"AAGG",
                    1,
                    rows(&smallvec![
                        Consequence::InframeInsertion,
                        Consequence::SpliceRegionVariant
                    ]),
                ),
            ]
        };
        let by_rank = |flag: bool| {
            config(|c| {
                c.pick = !flag;
                c.flag_pick = flag;
                c.pick_order = vec![PickCriterion::Rank];
            })
        };
        let mut record = tair_record();
        let ids: Vec<(usize, String)> = record_rows(&record)
            .iter()
            .map(|r| match *r {
                Row::Transcript { allele, tc } => (
                    allele,
                    record[allele].transcript_consequences[tc]
                        .transcript_id()
                        .to_string(),
                ),
                Row::Intergenic { allele } => (allele, "-".to_string()),
            })
            .collect();
        assert_eq!(
            &ids[..6],
            [
                (0, "AT4G12020.10".to_string()),
                (1, "AT4G12020.10".to_string()),
                (0, "AT4G12020.11".to_string()),
                (1, "AT4G12020.11".to_string()),
                (0, "AT4G12020.1".to_string()),
                (1, "AT4G12020.1".to_string()),
            ]
        );
        apply_to_record(&mut record, &by_rank(false));
        assert_eq!(kept(&record), vec![(1, "AT4G12020.10".to_string())]);
        let mut record = tair_record();
        apply_to_record(&mut record, &by_rank(true));
        assert_eq!(flagged(&record), vec![(1, "AT4G12020.10".to_string())]);
    }

    /// `pick_VariationFeatureOverlapAllele_per_gene` (`OutputFactory.pm:830-859`)
    /// ranks each gene's rows in the overlap list's order, so a full tie within
    /// a gene falls to the first transcript key (`:808`): on the TAIR10 record
    /// `ENSVATH10224189` (`4 9557579 T TAAATT`), whose nineteen AT4G16990
    /// isoforms are all `upstream_gene_variant` and `protein_coding`, Ensembl
    /// VEP 116.2 keeps `AT4G16990.10` (key `AT4G16990.10_108373`) under
    /// `--per_gene --pick_order biotype` and under `--per_gene --pick_order
    /// rank`, where the bare ids would put `AT4G16990.1` (`AT4G16990.1_108382`)
    /// first; the default order decides on `canonical`, `AT4G16990.7`.
    #[test]
    fn per_gene_full_tie_falls_to_the_first_transcript_key() {
        let isoforms: [(&str, u64); 19] = [
            ("AT4G16990.1", 108382),
            ("AT4G16990.2", 108378),
            ("AT4G16990.3", 108379),
            ("AT4G16990.4", 108384),
            ("AT4G16990.5", 108369),
            ("AT4G16990.6", 108376),
            ("AT4G16990.7", 108372),
            ("AT4G16990.8", 108371),
            ("AT4G16990.9", 108377),
            ("AT4G16990.10", 108373),
            ("AT4G16990.11", 108375),
            ("AT4G16990.12", 108386),
            ("AT4G16990.13", 108383),
            ("AT4G16990.14", 108381),
            ("AT4G16990.15", 108380),
            ("AT4G16990.16", 108387),
            ("AT4G16990.17", 108385),
            ("AT4G16990.18", 108370),
            ("AT4G16990.19", 108374),
        ];
        let gene_record = || {
            let rows = isoforms
                .iter()
                .map(|&(id, db_id)| {
                    let mut row = tc(id, "AT4G16990", smallvec![Consequence::UpstreamGeneVariant]);
                    row.canonical = id == "AT4G16990.7";
                    edit_transcript(&mut row, |t| t.db_id = Some(db_id));
                    row
                })
                .collect();
            vec![allele(b"AAATT", 0, rows)]
        };
        let per_gene = |order: &[PickCriterion]| {
            config(|c| {
                c.per_gene = true;
                c.pick_order = order.to_vec();
            })
        };
        for criterion in [PickCriterion::Biotype, PickCriterion::Rank] {
            let mut record = gene_record();
            apply_to_record(&mut record, &per_gene(&[criterion]));
            assert_eq!(
                kept(&record),
                vec![(0, "AT4G16990.10".to_string())],
                "{}",
                criterion.name()
            );
        }
        let mut record = gene_record();
        apply_to_record(&mut record, &per_gene(&PickCriterion::DEFAULT_ORDER));
        assert_eq!(kept(&record), vec![(0, "AT4G16990.7".to_string())]);
    }

    /// Ensembl VEP keys a motif row by the motif's database id
    /// (`MotifFeatureVariation.pm:70-76`, `VariationFeature.pm:937`) and the
    /// rejoin of a `--minimal` split stores a later allele's rows under its
    /// allele string, `_` and that key (`OutputFactory.pm:2316`), so `sort keys`
    /// (`VariationFeature.pm:750`) lists the first allele's digit-led keys
    /// before any letter-led prefix and after a `-/` one. On
    /// `21 6562008 CGA TGA,CGC` under `--regulatory --minimal --pick`, the first
    /// allele's `C/T` row and the second's `A/C_` row of ENSM00000575394 tie on
    /// every criterion, and Ensembl VEP 116.2 keeps the first allele's (`575394`
    /// before `A/C_575394`); with the second allele an insertion keyed `-/G_`
    /// the second allele's row leads and is kept.
    #[test]
    fn prefixed_motif_rows_tie_falls_to_the_first_key_by_database_id() {
        let motif = || TranscriptConsequence {
            motif_db_id: Some(575394),
            ..feature(
                "ENSM00000575394",
                FeatureType::MotifFeature,
                Consequence::TfBindingSiteVariant,
            )
        };
        let split_record = |prefix: &str| {
            let mut first = allele(b"T", 0, vec![motif()]);
            first.intergenic_overlap = true;
            let mut second = allele(b"C", 1, vec![motif()]);
            second.intergenic_overlap = true;
            second.rejoin_key_prefix = Some(prefix.into());
            vec![first, second]
        };
        let mut record = split_record("A/C_");
        let rows: Vec<usize> = record_rows(&record).iter().map(|r| r.allele()).collect();
        assert_eq!(rows, [0, 1, 0, 1]);
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENSM00000575394".to_string())]);

        let mut record = split_record("-/G_");
        let rows: Vec<usize> = record_rows(&record).iter().map(|r| r.allele()).collect();
        assert_eq!(rows, [1, 0, 0, 1]);
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(1, "ENSM00000575394".to_string())]);
    }

    /// A structural variant's regulatory rows stand in the order their
    /// features were met (`id_rank`, the order the source lists them), not in
    /// stable-id order, so a full tie between a regulatory file's record and
    /// the EMAR of the same span falls to the record met first
    /// (`OutputFactory.pm:556-560` ranks `get_all_StructuralVariationOverlaps`
    /// as pushed, `:808` returns the first); rows sharing a rank fall back to
    /// the feature's start and then its stable id, and a sequence variant's
    /// rows keep the stable-id order.
    #[test]
    fn structural_tie_falls_to_the_first_feature_met() {
        let regulatory = |id: &str, id_rank: u32, start: u64| TranscriptConsequence {
            id_rank,
            feature_start: start,
            feature_end: start + 230,
            ..feature(
                id,
                FeatureType::RegulatoryFeature,
                Consequence::RegulatoryRegionAblation,
            )
        };
        let structural = |rows: Vec<TranscriptConsequence>| {
            let mut v = allele(b"<DEL>", 0, rows);
            v.is_structural = true;
            v.intergenic_overlap = true;
            vec![v]
        };
        let ids = |record: &[InputVariant]| -> Vec<String> {
            record_rows(record)
                .into_iter()
                .filter_map(|row| match row {
                    Row::Transcript { allele, tc } => Some(
                        record[allele].transcript_consequences[tc]
                            .transcript_id()
                            .to_string(),
                    ),
                    Row::Intergenic { .. } => None,
                })
                .collect()
        };
        let pick = config(|c| c.pick = true);

        // The EMAR (`ENSR0...`) sorts before the enhancer (`ENSR2...`) by id;
        // the enhancer was met first.
        let mut met = structural(vec![
            regulatory("ENSR00001793149", 517, 5_011_320),
            regulatory("ENSR21_85ZF5", 1, 5_011_320),
        ]);
        assert_eq!(ids(&met), ["ENSR21_85ZF5", "ENSR00001793149"]);
        apply_to_record(&mut met, &pick);
        assert_eq!(kept(&met), vec![(0, "ENSR21_85ZF5".to_string())]);

        // One shared rank: the earlier start wins, then the stable id.
        let mut by_start = structural(vec![
            regulatory("ENSR21_B", 0, 200),
            regulatory("ENSR21_A", 0, 300),
            regulatory("ENSR00000001", 0, 300),
        ]);
        assert_eq!(ids(&by_start), ["ENSR21_B", "ENSR00000001", "ENSR21_A"]);
        apply_to_record(&mut by_start, &pick);
        assert_eq!(kept(&by_start), vec![(0, "ENSR21_B".to_string())]);

        // The order is total: a lower rank sorts before a higher one whatever
        // the starts.
        let mixed = structural(vec![
            regulatory("ENSR21_A", 0, 300),
            regulatory("ENSR21_B", 7, 100),
        ]);
        assert_eq!(ids(&mixed), ["ENSR21_A", "ENSR21_B"]);

        // A sequence variant's rows keep the stable-id order whatever their rank.
        let mut sequence = vec![allele(
            b"C",
            0,
            vec![
                regulatory("ENSR21_85ZF5", 1, 5_011_320),
                regulatory("ENSR00001793149", 517, 5_011_320),
            ],
        )];
        sequence[0].intergenic_overlap = true;
        assert_eq!(ids(&sequence), ["ENSR00001793149", "ENSR21_85ZF5"]);
        apply_to_record(&mut sequence, &pick);
        assert_eq!(kept(&sequence), vec![(0, "ENSR00001793149".to_string())]);
    }

    /// A full tie falls to the first row of Ensembl VEP's overlap list (`:808`),
    /// which holds the transcript rows, then the regulatory rows, then the motif
    /// rows, then the intergenic row (`VariationFeature.pm:855-864`), whatever
    /// the ids: under `--pick_order canonical` every row here scores 1, so the
    /// transcript row wins over the regulatory and motif rows whose ids sort
    /// before its own, and without a transcript the regulatory row wins over
    /// the motif row and over the allele's intergenic row.
    #[test]
    fn full_tie_falls_to_the_first_row_of_the_overlap_list() {
        let canonical_only = config(|c| {
            c.pick = true;
            c.pick_order = vec![PickCriterion::Canonical];
        });
        let regulatory = || {
            feature(
                "ENSR1",
                FeatureType::RegulatoryFeature,
                Consequence::RegulatoryRegionVariant,
            )
        };
        let motif = || {
            feature(
                "ENSM1",
                FeatureType::MotifFeature,
                Consequence::TfBindingSiteVariant,
            )
        };
        let transcript = || tc("ENST1", "G1", smallvec![Consequence::IntronVariant]);

        let mut record = vec![allele(b"G", 0, vec![motif(), regulatory(), transcript()])];
        assert_eq!(
            record_rows(&record),
            [
                Row::Transcript { allele: 0, tc: 2 },
                Row::Transcript { allele: 0, tc: 1 },
                Row::Transcript { allele: 0, tc: 0 }
            ]
        );
        apply_to_record(&mut record, &canonical_only);
        assert_eq!(kept(&record), vec![(0, "ENST1".to_string())]);

        let mut v = allele(b"G", 0, vec![motif(), regulatory()]);
        v.intergenic_overlap = true;
        let mut record = vec![v];
        assert_eq!(
            record_rows(&record),
            [
                Row::Transcript { allele: 0, tc: 1 },
                Row::Transcript { allele: 0, tc: 0 },
                Row::Intergenic { allele: 0 }
            ]
        );
        apply_to_record(&mut record, &canonical_only);
        assert_eq!(kept(&record), vec![(0, "ENSR1".to_string())]);
        assert!(!record[0].intergenic_overlap);

        let mut v = allele(b"G", 0, vec![regulatory()]);
        v.intergenic_overlap = true;
        let mut record = vec![v];
        apply_to_record(&mut record, &canonical_only);
        assert_eq!(kept(&record), vec![(0, "ENSR1".to_string())]);
        assert!(!record[0].intergenic_overlap);
    }

    /// The alleles of a `--minimal` split sort by the keys the rejoin stores
    /// their transcript variations under (`OutputFactory.pm:2316`: the allele
    /// string before the stable id for every allele after the first;
    /// `VariationFeature.pm:666` sorts the keys as strings): on `AAA AGA,AGA`,
    /// two `A/G` features tied on every criterion, `A/G_ENST...` sorts before
    /// the bare `ENST...`, so `--pick` keeps the second allele, the ALLELE_NUM 2
    /// row Ensembl VEP 116.2 writes for the minimal corpus's
    /// `multi_dup_trimmable`; on `A T,T,AG` the second `A/T` outranks the
    /// `-/G` insertion's frameshift row on rank and precedes the first allele's
    /// equal row (`multi_dup_plus_ins`, ALLELE_NUM 2). A prefix sorting after the
    /// id (`G/A_`) leaves the first allele's row in front, and each feature's
    /// reference row of `--process_ref_homs` stays ahead of its alternate. The
    /// intergenic rows follow in allele order whatever the prefixes: the rejoin
    /// appends a later allele's intergenic alleles to the record's one intergenic
    /// variation (`:2323-2338`) instead of rekeying them.
    #[test]
    fn minimal_split_alleles_sort_by_the_rejoin_key() {
        let same = || {
            tc(
                "ENST00000618007",
                "G1",
                smallvec![Consequence::MissenseVariant],
            )
        };
        let mut repeated = allele(b"G", 1, vec![same()]);
        repeated.rejoin_key_prefix = Some("A/G_".into());
        let mut record = vec![allele(b"G", 0, vec![same()]), repeated];
        let rows: Vec<usize> = record_rows(&record).iter().map(|r| r.allele()).collect();
        assert_eq!(rows, [1, 0]);
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(1, "ENST00000618007".to_string())]);

        let stop = || tc("ENST00000618007", "G1", smallvec![Consequence::StopGained]);
        let mut second = allele(b"T", 1, vec![stop()]);
        second.rejoin_key_prefix = Some("A/T_".into());
        let mut insertion = allele(
            b"G",
            2,
            vec![tc(
                "ENST00000618007",
                "G1",
                smallvec![Consequence::FrameshiftVariant],
            )],
        );
        insertion.rejoin_key_prefix = Some("-/G_".into());
        let mut record = vec![allele(b"T", 0, vec![stop()]), second, insertion];
        let rows: Vec<usize> = record_rows(&record).iter().map(|r| r.allele()).collect();
        assert_eq!(rows, [2, 1, 0]);
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(1, "ENST00000618007".to_string())]);

        let mut later = allele(b"A", 1, vec![same()]);
        later.rejoin_key_prefix = Some("G/A_".into());
        let mut record = vec![allele(b"G", 0, vec![same()]), later];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(kept(&record), vec![(0, "ENST00000618007".to_string())]);

        let mut repeated = allele(b"G", 1, vec![same()]);
        repeated.rejoin_key_prefix = Some("A/G_".into());
        let mut repeated_reference = repeated.reference_allele_row();
        repeated_reference.transcript_consequences = vec![same()];
        let first = allele(b"G", 0, vec![same()]);
        let mut first_reference = first.reference_allele_row();
        first_reference.transcript_consequences = vec![same()];
        let record = vec![first_reference, first, repeated_reference, repeated];
        let rows: Vec<(usize, bool)> = record_rows(&record)
            .iter()
            .map(|r| (r.allele(), record[r.allele()].reference_allele_row))
            .collect();
        assert_eq!(rows, [(2, true), (3, false), (0, true), (1, false)]);

        let regulatory = || {
            feature(
                "ENSR1",
                FeatureType::RegulatoryFeature,
                Consequence::RegulatoryRegionVariant,
            )
        };
        let mut second = allele(b"T", 1, vec![regulatory()]);
        second.rejoin_key_prefix = Some("C/T_".into());
        second.intergenic_overlap = true;
        let mut third = allele(b"-", 2, vec![regulatory()]);
        third.rejoin_key_prefix = Some("A/-_".into());
        third.intergenic_overlap = true;
        let record = vec![allele(b"G", 0, vec![]), second, third];
        assert_eq!(
            record_rows(&record),
            [
                Row::Transcript { allele: 2, tc: 0 },
                Row::Transcript { allele: 1, tc: 0 },
                Row::Intergenic { allele: 0 },
                Row::Intergenic { allele: 1 },
                Row::Intergenic { allele: 2 }
            ]
        );
    }

    /// `--most_severe` collapses the record to one row carrying the lowest-rank
    /// term of every allele's rows (`summary_only`, `:447`): missense over the
    /// second allele's synonymous and both alleles' intron terms.
    #[test]
    fn most_severe_collapses_the_record() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.most_severe = true));
        assert_eq!(kept(&record), vec![(0, "-".to_string())]);
        assert_eq!(
            record[0].summary_terms.as_deref(),
            Some(&[Consequence::MissenseVariant][..])
        );
        assert_eq!(
            record[0].most_severe_consequence,
            Some(Consequence::MissenseVariant)
        );
        assert_eq!(record[1].summary_terms, None);
        assert_eq!(record[1].most_severe_consequence, None);
    }

    /// `--summary` lists every term once, in rank order (`:436-443`).
    #[test]
    fn summary_lists_unique_terms_of_the_record_in_rank_order() {
        let mut record = two_gene_record();
        apply_to_record(&mut record, &config(|c| c.summary = true));
        assert_eq!(kept(&record), vec![(0, "-".to_string())]);
        assert_eq!(
            record[0].summary_terms.as_deref(),
            Some(
                &[
                    Consequence::MissenseVariant,
                    Consequence::SynonymousVariant,
                    Consequence::IntronVariant
                ][..]
            )
        );
    }

    /// The summary modes read every row left on the record (`summary_only`,
    /// `OutputFactory.pm:426-458`), regulatory and intergenic rows included:
    /// `--summary` lists a regulatory row's term after a transcript's missense
    /// term and, on an allele overlapping no transcript, before the allele's
    /// `intergenic_variant`; `--most_severe` keeps the first term of each list.
    #[test]
    fn regulatory_and_intergenic_rows_participate_in_summaries() {
        for most_severe in [true, false] {
            for with_transcript in [false, true] {
                let mut rows = vec![feature(
                    "ENSR1",
                    FeatureType::RegulatoryFeature,
                    Consequence::RegulatoryRegionVariant,
                )];
                if with_transcript {
                    rows.push(tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]));
                }
                let mut variant = allele(b"G", 0, rows);
                variant.intergenic_overlap = !with_transcript;
                variant.most_severe_consequence = Some(if with_transcript {
                    Consequence::MissenseVariant
                } else {
                    Consequence::RegulatoryRegionVariant
                });
                let mut record = vec![variant];
                apply_to_record(
                    &mut record,
                    &config(|c| {
                        c.most_severe = most_severe;
                        c.summary = !most_severe;
                    }),
                );
                let mut expected: ConsequenceList = if with_transcript {
                    smallvec![
                        Consequence::MissenseVariant,
                        Consequence::RegulatoryRegionVariant
                    ]
                } else {
                    smallvec![
                        Consequence::RegulatoryRegionVariant,
                        Consequence::IntergenicVariant
                    ]
                };
                if most_severe {
                    expected.truncate(1);
                }
                assert_eq!(record[0].most_severe_consequence, expected.first().copied());
                assert_eq!(record[0].summary_terms.as_ref(), Some(&expected));
                assert!(record[0].transcript_consequences.is_empty());
                assert!(!record[0].intergenic_overlap);
            }
        }
    }

    /// `--summary` beside `--most_severe`, a pair `--safe` lets through
    /// (`Config.pm:961`), lists every term: `summary_only` routes on either flag
    /// (`:387`) and tests `summary` first (`:436`), `most_severe` in its `else`
    /// (`:446`).
    #[test]
    fn summary_wins_over_most_severe() {
        let mut record = two_gene_record();
        apply_to_record(
            &mut record,
            &config(|c| {
                c.summary = true;
                c.most_severe = true;
            }),
        );
        assert_eq!(kept(&record), vec![(0, "-".to_string())]);
        assert_eq!(
            record[0].summary_terms.as_deref(),
            Some(
                &[
                    Consequence::MissenseVariant,
                    Consequence::SynonymousVariant,
                    Consequence::IntronVariant
                ][..]
            )
        );
    }

    /// `--most_severe --per_gene` summarises the per-gene picks (`:514` filters
    /// before `:387` summarises): G1's canonical intron row beats its
    /// non-canonical missense row, so the record's most severe term is the
    /// intron variant, not the missense the unpicked set holds.
    #[test]
    fn most_severe_follows_the_per_gene_picks() {
        let mut canonical = tc("ENST1", "G1", smallvec![Consequence::IntronVariant]);
        canonical.canonical = true;
        let severe = tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]);
        let mut record = vec![allele(b"G", 0, vec![canonical, severe])];
        apply_to_record(
            &mut record,
            &config(|c| {
                c.most_severe = true;
                c.per_gene = true;
            }),
        );
        assert_eq!(
            record[0].summary_terms.as_deref(),
            Some(&[Consequence::IntronVariant][..])
        );
    }

    /// A record whose alleles have no row (a reference-homozygous copy, a
    /// non-variant record) keeps the one summary row `summary_only` returns for
    /// it (`OutputFactory.pm:454-458`), with no term.
    #[test]
    fn summary_of_a_record_without_rows_is_a_row_with_no_term() {
        let mut v = allele(b"G", 0, vec![]);
        v.most_severe_consequence = None;
        let mut record = vec![v];
        apply_to_record(&mut record, &config(|c| c.summary = true));
        assert_eq!(record[0].summary_terms, Some(SmallVec::new()));
        assert_eq!(record[0].most_severe_consequence, None);
        assert_eq!(kept(&record), Vec::<(usize, String)>::new());
    }

    /// The default `--pick_order` is Ensembl VEP's (`Config.pm:311`); a token
    /// outside the list is `check_config`'s message (`:945-957`).
    #[test]
    fn pick_order_tokens_resolve_to_the_criteria() {
        assert_eq!(
            pick_order(None).unwrap(),
            PickCriterion::DEFAULT_ORDER.to_vec()
        );
        let names: Vec<&str> = PickCriterion::DEFAULT_ORDER
            .iter()
            .map(|c| c.name())
            .collect();
        assert_eq!(
            names,
            [
                "mane_select",
                "mane_plus_clinical",
                "canonical",
                "appris",
                "tsl",
                "biotype",
                "ccds",
                "rank",
                "length",
                "ensembl",
                "refseq"
            ]
        );
        for name in &names {
            assert_eq!(PickCriterion::parse(name).map(|c| c.name()), Some(*name));
        }
        assert_eq!(
            pick_order(Some(&["canonical".to_string(), "rank".to_string()])).unwrap(),
            vec![PickCriterion::Canonical, PickCriterion::Rank]
        );
        assert_eq!(
            pick_order(Some(&["x".to_string()])).unwrap_err(),
            "\"x\" is not a valid value for --pick_order\nValid values: mane_select, mane_plus_clinical, canonical, appris, tsl, biotype, ccds, rank, length, ensembl, refseq"
        );
    }

    /// `m/([A-Za-z]).+(\d+)/` with `+10` for an alternative (`:758-768`).
    #[test]
    fn appris_rank_follows_the_perl_pattern() {
        assert_eq!(appris_rank(Some("principal1")), 1);
        assert_eq!(appris_rank(Some("principal9")), 9);
        assert_eq!(appris_rank(Some("alternative1")), 11);
        assert_eq!(appris_rank(Some("alternative2")), 12);
        assert_eq!(appris_rank(None), 100);
        assert_eq!(appris_rank(Some("")), 100);
        assert_eq!(appris_rank(Some("principal")), 100);
        // `.+` wants a character between the letter and the digits.
        assert_eq!(appris_rank(Some("P1")), 100);
        assert_eq!(appris_rank(Some("p_1")), 1);
        // `$grade if $grade` runs after the +10, so only a zero principal grade
        // keeps the default.
        assert_eq!(appris_rank(Some("principal0")), 100);
        assert_eq!(appris_rank(Some("alternative0")), 10);
        // The greedy `.+` leaves `(\d+)` the last digit alone, and a trailing
        // non-digit does not stop the match.
        assert_eq!(appris_rank(Some("principal12")), 2);
        assert_eq!(appris_rank(Some("principal10")), 100);
        assert_eq!(appris_rank(Some("principal1x")), 1);
        // `substr($type, 0, 1) eq 'a'` is case-sensitive.
        assert_eq!(appris_rank(Some("Alternative1")), 1);
        // The minor-isoform spellings carry no digit.
        assert_eq!(appris_rank(Some("principalM")), 100);
        assert_eq!(appris_rank(Some("alternativeM")), 100);
    }

    /// Each criterion reads the field Ensembl VEP reads, missing values at
    /// their defaults (`:712-727`, `:735-768`).
    #[test]
    fn pick_info_reads_every_criterion() {
        let mut t = tc("ENST1", "G1", smallvec![Consequence::MissenseVariant]);
        t.canonical = true;
        edit_transcript(&mut t, |transcript| {
            transcript.mane_select = Some("NM_1.1".into());
            transcript.mane_plus_clinical = Some("NM_2.1".into());
            transcript.appris = Some("alternative2".into());
            transcript.tsl = Some("3".into());
            transcript.ccds = Some("CCDS1.1".into());
            set_lengths(transcript, Some(900), 2000);
        });
        let info = PickInfo::of_transcript(&t);
        assert_eq!(
            info,
            PickInfo {
                mane_select: 0,
                mane_plus_clinical: 0,
                canonical: 0,
                appris: 12,
                tsl: Rank::Text("3"),
                biotype: 0,
                ccds: 0,
                rank: i64::from(Consequence::MissenseVariant.rank()),
                length: -900,
                ensembl: 1,
                refseq: 1,
            }
        );
        let mut bare = tc("ENST2", "G1", smallvec![Consequence::IntronVariant]);
        edit_transcript(&mut bare, |transcript| {
            transcript.biotype = "lncRNA".into();
            transcript.ccds = Some("-".into());
            transcript.tsl = None;
            set_lengths(transcript, None, 700);
        });
        let info = PickInfo::of_transcript(&bare);
        assert_eq!(
            info,
            PickInfo {
                mane_select: 1,
                mane_plus_clinical: 1,
                canonical: 1,
                appris: 100,
                tsl: Rank::Number(100),
                biotype: 1,
                ccds: 1,
                rank: i64::from(Consequence::IntronVariant.rank()),
                length: -700,
                ensembl: 1,
                refseq: 1,
            }
        );
        assert_eq!(
            PickInfo::defaults(&[]),
            PickInfo {
                mane_select: 1,
                mane_plus_clinical: 1,
                canonical: 1,
                appris: 100,
                tsl: Rank::Number(100),
                biotype: 1,
                ccds: 1,
                rank: 1000,
                length: 0,
                ensembl: 1,
                refseq: 1,
            }
        );
    }

    /// A batch of two records is filtered record by record.
    #[test]
    fn apply_to_batch_treats_each_record_on_its_own() {
        let mut batch = two_gene_record();
        let mut second = two_gene_record();
        for v in second.iter_mut() {
            v.input_record = 2;
            v.start = 200;
            v.end = 200;
        }
        batch.extend(second);
        apply_to_batch(&mut batch, &config(|c| c.pick = true));
        assert_eq!(
            kept(&batch),
            vec![(0, "ENST1".to_string()), (0, "ENST1".to_string())]
        );
        assert_eq!(batch[0].input_record, 1);
        assert_eq!(batch[2].input_record, 2);
    }
}
