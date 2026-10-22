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

use std::cmp::Ordering;

use smallvec::SmallVec;
use vep_core::consequence::{Consequence, ConsequenceList, SourceCache, TranscriptConsequence};
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
/// pipeline lays them out).
pub fn apply_to_batch(batch: &mut [InputVariant], config: &FilterConfig) {
    if !config.any_active() {
        return;
    }
    for record in batch.chunk_by_mut(InputVariant::same_record) {
        apply_to_record(record, config);
    }
}

/// Applies the record-level steps to the alleles of one input record in Ensembl
/// VEP's order (`OutputFactory.pm:477-516`, then `:387`): the row filters
/// (`:482`, `:490-502`), the pick family over the rows they leave (`:514`),
/// then the summary of the rows left.
pub fn apply_to_record(record: &mut [InputVariant], config: &FilterConfig) {
    if !config.any_active() {
        return;
    }
    row_filter::apply_to_record(record, config.row_filters);
    pick_family(record, config);
    if config.most_severe || config.summary {
        // `summary_only` tests `summary` before `most_severe` (`:436`, `:446`), so
        // the pair `--safe` lets through lists every term.
        summarise(record, config.most_severe && !config.summary);
    }
}

/// One row of a record as Ensembl VEP's `VariationFeatureOverlapAllele` list
/// holds it: an allele against a transcript, or an allele's intergenic row.
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

/// The record's rows in Ensembl VEP's order: feature-major in the string order
/// of the keys the record's transcript variations are stored under
/// (`VariationFeature.pm:666` of ensembl-variation sorts them), allele-minor
/// (each feature lists its reference allele's row under `--process_ref_homs`
/// and then its alternate allele), then the intergenic row of every allele that
/// overlaps no transcript. A feature's key is its transcript's stable id; the
/// rejoin of a `--minimal` split stores every allele after the first under
/// `<its allele string>_<stable id>` (`OutputFactory.pm:2316`), so such an
/// allele's feature sorts by its allele string against the first allele's bare
/// stable id, before it when the string starts with a `-`, a digit or a letter
/// below the id's first. An allele whose chromosome has no transcripts has no
/// row.
fn record_rows(record: &[InputVariant]) -> Vec<Row> {
    let mut transcript_rows: Vec<(usize, usize)> = Vec::new();
    for (allele, variant) in record.iter().enumerate() {
        for tc in 0..variant.transcript_consequences.len() {
            transcript_rows.push((allele, tc));
        }
    }
    let feature_key = |&(allele, tc): &(usize, usize)| {
        let variant = &record[allele];
        variant
            .rejoin_key_prefix
            .as_deref()
            .unwrap_or("")
            .bytes()
            .chain(variant.transcript_consequences[tc].transcript_id.bytes())
    };
    transcript_rows.sort_by(|a, b| {
        feature_key(a)
            .cmp(feature_key(b))
            .then_with(|| {
                record[b.0]
                    .reference_allele_row
                    .cmp(&record[a.0].reference_allele_row)
            })
            .then_with(|| a.0.cmp(&b.0))
    });
    let mut rows: Vec<Row> = transcript_rows
        .into_iter()
        .map(|(allele, tc)| Row::Transcript { allele, tc })
        .collect();
    for (allele, variant) in record.iter().enumerate() {
        if variant.transcript_consequences.is_empty() && variant.most_severe_consequence.is_some() {
            rows.push(Row::Intergenic { allele });
        }
    }
    rows
}

/// The values `pick_worst_VariationFeatureOverlapAllele` ranks one row by, each
/// lower when better (`OutputFactory.pm:712-768`): the defaults for a row that is
/// not against a transcript, the transcript's attributes otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PickInfo {
    mane_select: i64,
    mane_plus_clinical: i64,
    canonical: i64,
    appris: i64,
    tsl: i64,
    biotype: i64,
    ccds: i64,
    rank: i64,
    length: i64,
    ensembl: i64,
    refseq: i64,
}

impl PickInfo {
    /// The defaults of `:712-727`, with the rank of the row's consequences
    /// (`:780-786`; 1000 for a row with none).
    fn defaults(consequences: &[Consequence]) -> Self {
        PickInfo {
            mane_select: 1,
            mane_plus_clinical: 1,
            canonical: 1,
            appris: 100,
            tsl: 100,
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
    fn of_transcript(tc: &TranscriptConsequence) -> Self {
        let mut info = PickInfo::defaults(&tc.consequences);
        info.mane_select = if tc.mane_select.is_some() { 0 } else { 1 };
        info.mane_plus_clinical = if tc.mane_plus_clinical.is_some() {
            0
        } else {
            1
        };
        info.canonical = if tc.canonical { 0 } else { 1 };
        info.biotype = if tc.biotype.as_deref() == Some("protein_coding") {
            0
        } else {
            1
        };
        // `$tr->{_ccds} && $tr->{_ccds} ne '-'` (`:740`).
        info.ccds = match tc.ccds.as_deref() {
            Some(s) if !s.is_empty() && s != "-" => 0,
            _ => 1,
        };
        match tc.source_cache {
            Some(SourceCache::Ensembl) => info.ensembl = 0,
            Some(SourceCache::RefSeq) => info.refseq = 0,
            None => {}
        }
        // Longer is better, so the length is negated (`:743-748`): the
        // translateable sequence of a translated transcript, else the spliced
        // length.
        info.length = -i64::from(tc.cds_length.unwrap_or(tc.cdna_length));
        // `m/tsl(\d+)/` then `$1 if $1` (`:751-755`): a zero stays at the default.
        if let Some(tsl) = tc.tsl.filter(|t| *t != 0) {
            info.tsl = i64::from(tsl);
        }
        info.appris = i64::from(appris_rank(tc.appris.as_deref()));
        info
    }

    fn value(&self, criterion: PickCriterion) -> i64 {
        match criterion {
            PickCriterion::ManeSelect => self.mane_select,
            PickCriterion::ManePlusClinical => self.mane_plus_clinical,
            PickCriterion::Canonical => self.canonical,
            PickCriterion::Appris => self.appris,
            PickCriterion::Tsl => self.tsl,
            PickCriterion::Biotype => self.biotype,
            PickCriterion::Ccds => self.ccds,
            PickCriterion::Rank => self.rank,
            PickCriterion::Length => self.length,
            PickCriterion::Ensembl => self.ensembl,
            PickCriterion::RefSeq => self.refseq,
        }
    }

    /// Lexicographic order over the criteria in `order`.
    fn compare(&self, other: &Self, order: &[PickCriterion]) -> Ordering {
        for &criterion in order {
            match self.value(criterion).cmp(&other.value(criterion)) {
                Ordering::Equal => {}
                unequal => return unequal,
            }
        }
        Ordering::Equal
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
/// `candidates` (`OutputFactory.pm:703-811`): the lexicographic minimum over the
/// criteria in order, which is what the per-category stable sort and shrink of
/// `:789-806` leaves; on a full tie the first candidate (`:808`).
fn pick_worst(candidates: &[usize], infos: &[PickInfo], order: &[PickCriterion]) -> usize {
    let mut best = candidates[0];
    for &candidate in &candidates[1..] {
        if infos[candidate].compare(&infos[best], order) == Ordering::Less {
            best = candidate;
        }
    }
    best
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
            Row::Transcript { .. } => transcript_rows.push(c),
            Row::Intergenic { .. } => picked.push(c),
        }
    }
    let gene = |c: usize| match rows[c] {
        Row::Transcript { allele, tc } => {
            record[allele].transcript_consequences[tc].gene_id.clone()
        }
        Row::Intergenic { .. } => unreachable!("only transcript rows are grouped by gene"),
    };
    for group in group_by(&transcript_rows, gene) {
        picked.push(pick_worst(&group, infos, order));
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
            Row::Intergenic { allele } => {
                PickInfo::defaults(record[allele].most_severe_consequence.as_slice())
            }
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
            .map(|group| pick_worst(group, &infos, order))
            .collect()
    };

    if flag {
        for &c in &picked {
            match rows[c] {
                Row::Transcript { allele, tc } => {
                    push_flag(&mut record[allele].transcript_consequences[tc], "PICK")
                }
                Row::Intergenic { allele } => record[allele].intergenic_pick = true,
            }
        }
    } else {
        retain_rows(record, &rows, &picked);
    }
}

/// Keeps the rows at `picked` and removes every other row of the record: an
/// allele keeps its picked transcript consequences in their order, and an allele
/// left with no row loses its consequence verdict so no writer prints one for it.
fn retain_rows(record: &mut [InputVariant], rows: &[Row], picked: &[usize]) {
    for (allele, variant) in record.iter_mut().enumerate() {
        if variant.transcript_consequences.is_empty() {
            let keeps_intergenic = picked
                .iter()
                .any(|&c| rows[c] == Row::Intergenic { allele });
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

/// Appends a flag to a transcript consequence's flags.
fn push_flag(tc: &mut TranscriptConsequence, flag: &str) {
    let mut v: Vec<String> = tc.flags.to_vec();
    v.push(flag.to_string());
    tc.flags = v.into();
}

/// `summary_only` (`OutputFactory.pm:426-458`): the consequence terms of the rows
/// left on the record, sorted by rank, each once (`--summary`) or the first alone
/// (`--most_severe`), carried by the record's first allele as its one row while
/// every other allele keeps none. A record with no row keeps its one summary row
/// with no term and warns `Unable to assign consequence type` (`:454-458` warn
/// and return the record's hash with no `Consequence`): a reference-homozygous
/// individual's copy, or a non-variant record.
fn summarise(record: &mut [InputVariant], most_severe_only: bool) {
    let mut terms: Vec<Consequence> = Vec::new();
    for row in record_rows(record) {
        let row_terms: &[Consequence] = match row {
            Row::Transcript { allele, tc } => {
                &record[allele].transcript_consequences[tc].consequences
            }
            Row::Intergenic { allele } => record[allele].most_severe_consequence.as_slice(),
        };
        terms.extend_from_slice(row_terms);
    }
    terms.sort_by_key(|c| c.rank());
    terms.dedup();
    if most_severe_only {
        terms.truncate(1);
    }
    if terms.is_empty() {
        crate::warnings::global().warn("Unable to assign consequence type");
    }
    for variant in record.iter_mut() {
        variant.transcript_consequences.clear();
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
    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, FeatureType, Impact};
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
            transcript_id: transcript_id.into(),
            gene_id: gene_id.into(),
            biotype: Some("protein_coding".into()),
            consequences,
            impact,
            feature_type: FeatureType::Transcript,
            strand: 1,
            cdna_length: 1000,
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
                out.push((v.allele_index, t.transcript_id.to_string()));
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
                if t.flags.iter().any(|f| f == "PICK") {
                    out.push((v.allele_index, t.transcript_id.to_string()));
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
        short_cds.cds_length = Some(300);
        short_cds.cdna_length = 5000;
        let mut long_cds = tc("ENST2", "G1", smallvec![Consequence::MissenseVariant]);
        long_cds.cds_length = Some(600);
        long_cds.cdna_length = 1000;
        let mut record = vec![allele(b"G", 0, vec![short_cds, long_cds])];
        apply_to_record(&mut record, &config(|c| c.pick = true));
        assert_eq!(
            kept(&record),
            vec![(0, "ENST2".to_string())],
            "the longer CDS wins although its cDNA is shorter"
        );

        let mut short = tc("ENST1", "G1", smallvec![Consequence::IntronVariant]);
        short.biotype = Some("lncRNA".into());
        short.cdna_length = 800;
        let mut long = tc("ENST2", "G1", smallvec![Consequence::IntronVariant]);
        long.biotype = Some("lncRNA".into());
        long.cdna_length = 2400;
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
        assert_eq!(a.source_cache, None);
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
        refseq.source_cache = Some(SourceCache::RefSeq);
        let info = PickInfo::of_transcript(&refseq);
        assert_eq!((info.ensembl, info.refseq), (1, 0));
    }

    /// Two rows equal on every criterion: the pick is the first in Ensembl VEP's
    /// overlap-allele order, the lexicographically first stable id (`:808`
    /// returns the first of what the stable sorts left; `VariationFeature.pm:666`
    /// of ensembl-variation sorts the transcript keys).
    #[test]
    fn full_tie_falls_to_the_lexicographically_first_stable_id() {
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
    /// reference row of `--process_ref_homs` stays ahead of its alternate.
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

    /// The summary modes read transcript and intergenic rows alone: the
    /// regulatory and motif rows Ensembl VEP's `summary_only` would also fold in
    /// (`:431`, over every overlap allele) do not exist in vep-rs, which is why
    /// the `most-severe-with-regulatory` corpus row documents a subset.
    #[test]
    fn most_severe_documents_the_regulatory_exclusion() {
        let mut record = vec![allele(b"G", 0, vec![])];
        apply_to_record(&mut record, &config(|c| c.most_severe = true));
        assert_eq!(
            record[0].summary_terms.as_deref(),
            Some(&[Consequence::IntergenicVariant][..]),
            "an intergenic record summarises to its intergenic term"
        );
        let regulatory = [
            Consequence::TfbsAblation,
            Consequence::TfbsAmplification,
            Consequence::TfBindingSiteVariant,
            Consequence::RegulatoryRegionAblation,
            Consequence::RegulatoryRegionAmplification,
            Consequence::RegulatoryRegionVariant,
        ];
        assert!(
            regulatory
                .iter()
                .all(|c| c.rank() > Consequence::DownstreamGeneVariant.rank()),
            "every regulatory term ranks below the transcript terms, so it would change a summary only when no transcript row remains"
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
        t.mane_select = Some("NM_1.1".into());
        t.mane_plus_clinical = Some("NM_2.1".into());
        t.canonical = true;
        t.appris = Some("alternative2".into());
        t.tsl = Some(3);
        t.ccds = Some("CCDS1.1".into());
        t.cds_length = Some(900);
        t.cdna_length = 2000;
        let info = PickInfo::of_transcript(&t);
        assert_eq!(
            info,
            PickInfo {
                mane_select: 0,
                mane_plus_clinical: 0,
                canonical: 0,
                appris: 12,
                tsl: 3,
                biotype: 0,
                ccds: 0,
                rank: i64::from(Consequence::MissenseVariant.rank()),
                length: -900,
                ensembl: 1,
                refseq: 1,
            }
        );
        let mut bare = tc("ENST2", "G1", smallvec![Consequence::IntronVariant]);
        bare.biotype = Some("lncRNA".into());
        bare.ccds = Some("-".into());
        bare.tsl = Some(0);
        bare.cdna_length = 700;
        let info = PickInfo::of_transcript(&bare);
        assert_eq!(
            info,
            PickInfo {
                mane_select: 1,
                mane_plus_clinical: 1,
                canonical: 1,
                appris: 100,
                tsl: 100,
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
                tsl: 100,
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
