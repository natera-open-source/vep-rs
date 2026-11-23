// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! UTRAnnotator plugin: the effect of a 5' UTR variant on upstream open reading
//! frames (uORFs), as the Ensembl VEP plugin `UTRAnnotator.pm` computes it.
//!
//! **Usage:** `--plugin UTRAnnotator,file=<uORF evidence file>[,max_overlap=<percent>]`
//!
//! A row carrying `5_prime_UTR_variant` (`:144`) gets five keys: `5UTR_consequence`
//! (one flag), `5UTR_annotation` (that flag's effects), and the counts
//! `Existing_uORFs`, `Existing_OutOfFrame_oORFs` and `Existing_InFrame_oORFs` of
//! the reference 5' UTR (`:1029-1070`). The five analyses are `uAUG_gained`
//! (`:307-435`), `uSTOP_gained` (`:437-574`, substitutions alone), `uSTOP_lost`
//! (`:576-715`), `uAUG_lost` (`:717-848`) and `uFrameshift` (`:850-1018`, indels
//! alone), each over the 5' UTR sequence with the allele spliced in
//! (`:1189-1200`).
//!
//! The 5' UTR and coding sequences are read from the reference FASTA through the
//! transcript's exon model, as Ensembl VEP reads them through `spliced_seq`
//! (`Transcript.pm:823-870`, `:903-926`, `:1746-1762`): every UTR exon segment in
//! transcript order, the coding sequence from the coding start to the coding end
//! with one `N` per base of the start exon's phase. The plugin cannot annotate
//! without that reference, so a run that has none refuses it when it loads,
//! where Ensembl VEP's module prints a notice and runs on (`:88`). Each
//! transcript's sequences and reference codon scans are read once per batch and
//! shared by every row and allele of that transcript.
//!
//! The uORF evidence file (`file=`) has six tab-separated columns `chr pos gene
//! strand type stop_pos` with a header (`:98-99`); its chromosome names open with
//! `ch` (`:101`; the key lookup of `:1321` expects the `chr` prefix) and its
//! positions are numbers (`:102`). A path that does not exist is an initialization
//! error, where Ensembl VEP prints a warning and runs without evidence
//! (`:109-112`); without `file=` every `Evidence` value is `NA` (`:551`, `:693`,
//! `:826`, `:990`). A token without `=` is not read, as `params_to_hash` leaves it
//! out (`:90`), and is logged. `max_overlap` (`:115`, default 100) drops a UTR exon
//! the variant covers at or above that percentage (`:184-188`); its value is read
//! as that numeric comparison reads it, the text's leading number and 0 when none
//! leads it (`80%` is 80, `abc` is 0), a value that is not wholly a number logged.
//!
//! Two renderings are fixed where Ensembl VEP's depend on hash order: the
//! `key=value` pairs inside `5UTR_annotation` are in byte order of their keys,
//! capitalised keys first (`:282-284` writes them in hash order), and when more
//! than one analysis fires the flag is the first of `uAUG_gained`, `uSTOP_gained`,
//! `uSTOP_lost`, `uAUG_lost`, `uFrameshift` (`:262-267` keeps the last in hash
//! order). The window searched for a gained start codon is sized by the alternate
//! allele (`atg_window`); `:339` sizes it by the reference allele, which misses a
//! start codon inside an inserted sequence and reports an existing one beside a
//! deletion.
//!
//! An insertion is read from every index at which its mutated sequence can be
//! written (`insertion_writings`): inserted bases that repeat the bases beside
//! them give one sequence at several indices. A start codon stands when its
//! three bases stand whole, at its own index before every writing or moved on
//! by the allele's length after one of them.
//! It is reported lost only when every writing parts it, where `:764-769` pass
//! an insertion at the codon's three indices and `:776-778` report the codon
//! lost from the three bases at its old index. A start codon that stands at
//! its own index is not reported gained, where `:346-353` name the first start
//! codon of the window; of two equal copies of a repeated start codon the
//! gained one is the copy nearer the cap.
//!
//! Where a writing lies at a start codon's three indices or directly after it,
//! the stop-codon and frameshift analyses read all the writings, so that one
//! sequence has one row however the record writes it. A uORF is
//! shifted when a writing lies inside it, after a start codon that stands at
//! its own index, unless it keeps its length: one that still ends on a stop
//! codon at its own index, every codon before its last unchanged, is left out,
//! where `:921` lists it for a record written before that index. A stop codon
//! is lost when a writing lies inside it and the mutated sequence reads no
//! stop codon at its index. When another writing lies at or before its first
//! base, the stop codon stands whole after the inserted bases and is lost only
//! at a length that is no multiple of three; `:628-633` read the stop codon's
//! index only when the record itself writes the insertion inside it. Any other
//! allele is read at the record's own index.
//!
//! A symbolic deletion whose stated span lies inside one 5' UTR exon segment is
//! annotated as the sequence deletion of that span: the analyses read the count
//! of removed bases and the reference UTR alone (`:1189-1200`), so the row
//! carries the five keys of the same deletion given with explicit alleles. A
//! record with several alternates states that deletion for the alternate whose
//! own text is a deletion's (`<DEL>`, `<CN0>`) when it writes `END` as the
//! span's last base; its other alternates state no sequence. An inversion
//! whose stated span lies inside one segment is annotated as that span
//! replaced by its reverse complement, an allele of equal length that the
//! reference holds. A plain symbolic insertion (`<INS>`) that states its length
//! gives no bases: where the length is no multiple of three and the point lies
//! inside a uORF or an overlapping ORF, the row carries the frameshift with the
//! keys that read no inserted base, and `alt_type=uORF` where a stop codon past
//! the point lands in frame whatever is inserted (`:889`, `:921-923`,
//! `:973-974`); `alt_type_length` is never written. Every other structural
//! variant's 5' UTR rows carry the counts and `-` alone: a deletion or an
//! inversion that ends outside the UTR or in another exon segment, a record
//! with the `IMPRECISE` flag or with a confidence interval around its position
//! or its end, an inversion or an insertion of a record with several
//! alternates, an insertion that writes another end than its position or its
//! position plus its length, a mobile-element insertion, a duplication.
//! Ensembl VEP reads the symbolic allele's text (`<DEL>`) as a deleted sequence
//! of its own length (`:151-153`) and analyses that.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use indexmap::IndexMap;
use tracing::warn;

use crate::params::{perl_numeric, warn_unread_parameters, PluginParams};
use crate::tabix::BoundedLines;
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::codon::reverse_complement;
use vep_core::consequence::{Consequence, LofteeContext, TranscriptConsequence};
use vep_core::sv_routing::classify_symbolic_alt;
use vep_core::transcript::Exon;
use vep_core::variant::{InputVariant, VariantClass, WrittenEnd};

const NAME: &str = "UTRAnnotator";

/// The output keys (`:131-135`).
const CONSEQUENCE_KEY: &str = "5UTR_consequence";
const ANNOTATION_KEY: &str = "5UTR_annotation";
const EXISTING_UORFS_KEY: &str = "Existing_uORFs";
const EXISTING_OUT_OF_FRAME_KEY: &str = "Existing_OutOfFrame_oORFs";
const EXISTING_IN_FRAME_KEY: &str = "Existing_InFrame_oORFs";

/// The flag each analysis writes (`:426`, `:562`, `:704`, `:838`, `:1003`).
const UAUG_GAINED_FLAG: &str = "5_prime_UTR_premature_start_codon_gain_variant";
const USTOP_GAINED_FLAG: &str = "5_prime_UTR_uORF_stop_codon_gain_variant";
const USTOP_LOST_FLAG: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
const UAUG_LOST_FLAG: &str = "5_prime_UTR_premature_start_codon_loss_variant";
const UFRAMESHIFT_FLAG: &str = "5_prime_UTR_uORF_frameshift_variant";

/// The value of a key with nothing to report (`:291-292`).
const NONE: &str = "-";

/// A value Ensembl VEP writes where a distance has no stop codon or the
/// evidence file is absent.
const NA: &str = "NA";

/// The effects of one analysis: one `key=value` set per uORF it fired on.
type Effects = Vec<Vec<(&'static str, String)>>;

/// One analysis's outcome: its flag and the effects behind it.
struct Analysis {
    flag: &'static str,
    effects: Effects,
}

/// UTRAnnotator plugin.
pub struct UtrAnnotatorPlugin {
    /// The `chr<N>:<pos>` keys of the evidence file; `None` without `file=`.
    evidence: Option<HashSet<String>>,
    /// The `max_overlap` percentage (`:115`), read as the comparison
    /// `$UTR_overlap_perc >= $self->{max_overlap}` (`:188`) reads the value: its
    /// leading number, 0 when none leads it.
    max_overlap: f64,
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
}

impl UtrAnnotatorPlugin {
    /// Create a new uninitialized UTRAnnotator plugin.
    pub fn new() -> Self {
        Self {
            evidence: None,
            max_overlap: 100.0,
            reference_fasta: None,
        }
    }

    /// Reads the evidence file into its `chr:pos` key set (`:92-108`): the header
    /// line `chr pos ...` is skipped (`:99`); a chromosome that does not open with `ch`
    /// (`:101`) or a position that is not a number (`:102`) is refused. The lines
    /// are split as the bytes the module reads, with the length of a line capped
    /// as every flat file's is ([`BoundedLines`]); a byte that is not UTF-8 in
    /// the four columns past the position changes nothing.
    fn read_evidence(path: &Path) -> Result<HashSet<String>, PluginError> {
        let file = std::fs::File::open(path).map_err(|e| {
            PluginError::Init(format!("{NAME}: failed to open {}: {e}", path.display()))
        })?;
        let mut keys = HashSet::new();
        let mut lines = BoundedLines::new(BufReader::new(file), path);
        let mut line_number = 0u64;
        while let Some(line) = lines.next_line() {
            line_number += 1;
            let line = line?;
            let mut columns = line.split(|byte| *byte == b'\t');
            let chr = columns.next().unwrap_or(b"");
            let pos = columns.next().unwrap_or(b"");
            if chr == b"chr" && pos == b"pos" {
                continue;
            }
            let (chr, pos) = (String::from_utf8_lossy(chr), String::from_utf8_lossy(pos));
            if !chr.starts_with("ch") {
                return Err(PluginError::Init(format!(
                    "{NAME}: chromosome '{chr}' on line {} of {} does not carry the chr prefix",
                    line_number,
                    path.display()
                )));
            }
            if pos.parse::<f64>().is_err() {
                return Err(PluginError::Init(format!(
                    "{NAME}: position '{pos}' on line {} of {} is not a number",
                    line_number,
                    path.display()
                )));
            }
            keys.insert(format!("{chr}:{pos}"));
        }
        Ok(keys)
    }

    /// Whether the uORF starting at `start_pos` of the reference 5' UTR has
    /// evidence (`:1315-1325`): the key is `chr<chr>:<genomic position of the
    /// start codon>`, the prefix added when the chromosome name lacks it (`:1321`),
    /// `NA` without a file (`:551`).
    fn evidence(&self, chr: &str, positions: &UtrPositions, start_pos: usize) -> String {
        let Some(evidence) = &self.evidence else {
            return NA.to_string();
        };
        let prefix = if chr.to_ascii_lowercase().contains("chr") {
            ""
        } else {
            "chr"
        };
        let key = match positions.genomic_of(start_pos) {
            Some(genomic) => format!("{prefix}{chr}:{genomic}"),
            // `$utr_pos{$start_pos}` is undefined past the kept segments and the
            // query `chr:` matches nothing (`:1319-1322`).
            None => format!("{prefix}{chr}:"),
        };
        if evidence.contains(&key) {
            "True".to_string()
        } else {
            "False".to_string()
        }
    }

    /// Annotates one 5' UTR row (`:140-300`).
    fn annotate(
        &self,
        model: &TranscriptUtr,
        chr: &str,
        variant: &InputVariant,
    ) -> IndexMap<String, String> {
        // Every UTR exon dropped: no key at all (`:193`).
        let Some(positions) = model.kept_positions(variant.start, variant.end, self.max_overlap)
        else {
            return IndexMap::new();
        };

        // A structural variant is analysed as the deletion, the inversion or the
        // insertion its record states ([`stated_deletion`], [`stated_inversion`],
        // [`stated_insertion`]), or not at all.
        let chosen = if !variant.is_structural {
            let ref_len = if variant.ref_allele == b"-" {
                0
            } else {
                variant.ref_allele.len()
            };
            let alt_coding = oriented_alt(variant.alt_allele(), model.strand);
            positions
                .allele_span(model.strand, variant.start, ref_len)
                .and_then(|span| {
                    self.first_analysis(model, chr, &positions, span, ref_len, &alt_coding)
                })
        } else if let Some(length) = stated_deletion(variant) {
            positions
                .allele_span(model.strand, variant.start, length)
                // Both ends lie in one kept segment: the spliced distance
                // between them is then the stated span, whose bases the
                // reference UTR holds.
                .filter(|&(mut_pos, end_pos)| end_pos.checked_sub(mut_pos) == Some(length - 1))
                .and_then(|span| self.first_analysis(model, chr, &positions, span, length, &[]))
        } else if let Some(length) = stated_inversion(variant) {
            positions
                .allele_span(model.strand, variant.start, length)
                .filter(|&(mut_pos, end_pos)| end_pos.checked_sub(mut_pos) == Some(length - 1))
                .and_then(|(mut_pos, end_pos)| {
                    // The inverted span read along the transcript is the reverse
                    // complement of the bases the UTR holds there, on either strand.
                    let reference = model.utr.get(mut_pos..=end_pos)?;
                    let inverted = reverse_complement(reference);
                    // A base the inversion leaves as it was at either end of the span
                    // is no part of the change; a span equal to its own reverse
                    // complement changes nothing.
                    let same = |pair: (&u8, &u8)| pair.0 == pair.1;
                    let before = reference
                        .iter()
                        .zip(&inverted)
                        .take_while(|p| same(*p))
                        .count();
                    if before == length {
                        return None;
                    }
                    let after = reference
                        .iter()
                        .rev()
                        .zip(inverted.iter().rev())
                        .take_while(|p| same(*p))
                        .count();
                    let changed = before..length - after;
                    let span = (mut_pos + changed.start, mut_pos + changed.end - 1);
                    self.first_analysis(
                        model,
                        chr,
                        &positions,
                        span,
                        changed.len(),
                        &inverted[changed],
                    )
                })
        } else if let Some(length) = stated_insertion(variant) {
            positions
                .allele_span(model.strand, variant.start, 0)
                .and_then(|(index, _)| {
                    self.frameshift_by_length(model, chr, &positions, index, length)
                })
        } else {
            None
        };

        let (flag, annotation) = match chosen {
            Some(analysis) => (analysis.flag.to_string(), render_effects(&analysis.effects)),
            None => (NONE.to_string(), NONE.to_string()),
        };
        IndexMap::from([
            (CONSEQUENCE_KEY.to_string(), flag),
            (ANNOTATION_KEY.to_string(), annotation),
            (
                EXISTING_UORFS_KEY.to_string(),
                model.counts.uorfs.to_string(),
            ),
            (
                EXISTING_OUT_OF_FRAME_KEY.to_string(),
                model.counts.out_of_frame_oorfs.to_string(),
            ),
            (
                EXISTING_IN_FRAME_KEY.to_string(),
                model.counts.in_frame_oorfs.to_string(),
            ),
        ])
    }

    /// The first analysis that fires on an allele that puts `alt_coding` in the
    /// place of the `ref_len` UTR bases between the two indices of `span`, in
    /// the fixed order `uAUG_gained`, `uSTOP_gained`, `uSTOP_lost`,
    /// `uAUG_lost`, `uFrameshift`, where `:262-267` keeps the last fired flag
    /// in hash order.
    fn first_analysis(
        &self,
        model: &TranscriptUtr,
        chr: &str,
        positions: &UtrPositions,
        (mut_pos, end_pos): (usize, usize),
        ref_len: usize,
        alt_coding: &[u8],
    ) -> Option<Analysis> {
        let mutated = mutated_utr(&model.utr, mut_pos, ref_len, alt_coding);
        let scan = CodonScan::new(&mutated);
        let writings = (ref_len == 0 && !alt_coding.is_empty())
            .then(|| insertion_writings(&model.utr, &mutated));
        let beside_start_codon = writings.filter(|&(first, last)| {
            model
                .scan
                .atgs
                .iter()
                .any(|&start_pos| first <= start_pos + 3 && start_pos <= last)
        });
        let row = MutatedUtr {
            plugin: self,
            chr,
            positions,
            model,
            mutated,
            scan,
            mut_pos,
            end_pos,
            ref_len,
            alt_len: alt_coding.len(),
            writings,
            beside_start_codon,
        };
        let analyses = [
            row.uaug_gained(),
            row.ustop_gained(),
            row.ustop_lost(),
            row.uaug_lost(),
            row.uframeshift(),
        ];
        analyses.into_iter().flatten().next()
    }

    /// The frameshift an insertion of `length` bases the record does not give
    /// causes before UTR index `index`. `uFrameshift` fires on the length and
    /// the reference UTR alone: the length is no multiple of three (`:889`,
    /// `:923`) and the index lies inside a uORF or an overlapping ORF
    /// (`:904-921`). Each effect carries the keys that read no inserted base
    /// ([`reference_frame_keys`]) and `alt_type=uORF` where a stop codon past
    /// the index lands in the start codon's frame once moved by `length`,
    /// which makes the ORF a uORF whatever is inserted (`:973-974`). An
    /// `alt_type` the inserted bases decide, `alt_type_length` (`:982-987`)
    /// and the four other analyses, which read the inserted bases, are left
    /// out.
    fn frameshift_by_length(
        &self,
        model: &TranscriptUtr,
        chr: &str,
        positions: &UtrPositions,
        index: usize,
        length: usize,
    ) -> Option<Analysis> {
        if length.is_multiple_of(3) {
            return None;
        }
        let mut reference_stops: Option<Vec<usize>> = None;
        let mut effects = Effects::new();
        for &start_pos in &model.scan.atgs {
            if index < start_pos + 3 || index > check_point(model, start_pos) {
                continue;
            }
            let mut effect = reference_frame_keys(self, model, chr, positions, start_pos);
            let stops = reference_stops.get_or_insert_with(|| stop_positions(&model.utr));
            let shifted_stop_in_frame = stops
                .iter()
                .any(|&stop| stop >= index && (stop + length - start_pos).is_multiple_of(3));
            if shifted_stop_in_frame {
                effect.push(("alt_type", "uORF".to_string()));
            }
            effects.push(effect);
        }
        (!effects.is_empty()).then_some(Analysis {
            flag: UFRAMESHIFT_FLAG,
            effects,
        })
    }
}

impl Default for UtrAnnotatorPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl crate::traits::sealed::Sealed for UtrAnnotatorPlugin {}

impl BuiltinPlugin for UtrAnnotatorPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['VariationFeature',
    /// 'StructuralVariationFeature']` (UTRAnnotator.pm:80-82).
    fn annotates_structural_variants(&self) -> bool {
        true
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                CONSEQUENCE_KEY.to_string(),
                "Variant consequence from UTRAnnotator".to_string(),
            ),
            (
                ANNOTATION_KEY.to_string(),
                "Variant annotation from UTRAnnotator".to_string(),
            ),
            (
                EXISTING_UORFS_KEY.to_string(),
                "The number of existing uORFs with a stop codon within the 5 prime UTR".to_string(),
            ),
            (
                EXISTING_OUT_OF_FRAME_KEY.to_string(),
                "The number of existing out-of-frame overlapping ORFs (OutOfFrame oORF) at the 5 prime UTR".to_string(),
            ),
            (
                EXISTING_IN_FRAME_KEY.to_string(),
                "The number of existing inFrame overlapping ORFs (inFrame oORF) at the 5 prime UTR".to_string(),
            ),
        ]
    }

    /// `file=` and `max_overlap=` through `params_to_hash` (`:90`); a token
    /// without `=` leaves that hash empty (BaseVepPlugin.pm:123), so it is not
    /// read and the plugin runs without the file, as Ensembl VEP's does.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file", "max_overlap"], params)?;
        warn_unread_parameters(
            NAME,
            params.positional(),
            "file=<path> and max_overlap=<percent>",
        );
        match params.get("file") {
            Some(file) => {
                let path = Path::new(file);
                if !path.is_file() {
                    return Err(PluginError::Init(format!(
                        "{NAME}: uORF file {} not found",
                        path.display()
                    )));
                }
                self.evidence = Some(Self::read_evidence(path)?);
            }
            None => warn!(
                "{NAME}: no file= given; every Evidence value is NA (the curated uORF files are \
                 uORF_5UTR_GRCh37_PUBLIC.txt and uORF_5UTR_GRCh38_PUBLIC.txt)"
            ),
        }
        if let Some(value) = params.get("max_overlap") {
            self.max_overlap = perl_numeric(NAME, "max_overlap", value);
        }
        Ok(())
    }

    /// The 5' UTR and coding sequences come from the reference alone, so the
    /// registry refuses the plugin in a run without one. Ensembl VEP's module
    /// prints a notice and runs on (`:88`).
    fn requires_reference_fasta(&self) -> bool {
        true
    }

    fn set_reference_fasta(&mut self, fasta: Option<Arc<vep_fasta::IndexedFasta>>) {
        self.reference_fasta = fasta;
    }

    /// Every 5' UTR transcript row of the batch, the transcript's sequences and
    /// reference codon scans read once per batch.
    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        let Some(fasta) = self.reference_fasta.as_deref() else {
            return Ok(());
        };
        let mut models: HashMap<Arc<str>, Option<TranscriptUtr>> = HashMap::new();
        for variant in variants.iter_mut() {
            for tc_idx in 0..variant.transcript_consequences.len() {
                let tc = &variant.transcript_consequences[tc_idx];
                if !self.annotates_row(tc)
                    || !tc.consequences.contains(&Consequence::FivePrimeUtrVariant)
                {
                    continue;
                }
                let Some(ctx) = tc.loftee_ctx.as_deref() else {
                    continue;
                };
                let transcript_id = match tc.transcript() {
                    Some(transcript) => Arc::clone(&transcript.stable_id),
                    None => Arc::from(tc.transcript_id()),
                };
                let Some(model) = models
                    .entry(transcript_id)
                    .or_insert_with(|| TranscriptUtr::read(ctx, &variant.chr, fasta))
                    .as_ref()
                else {
                    continue;
                };
                let result = self.annotate(model, &variant.chr, variant);
                let row = &mut variant.transcript_consequences[tc_idx];
                for (key, value) in result {
                    row.set_plugin_value(key, value, false);
                }
            }
        }
        Ok(())
    }

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // `run_batch` annotates every row itself.
        Ok(IndexMap::new())
    }
}

/// The alternate allele as the transcript reads it (`:1167-1177`): `-` is the
/// empty sequence and a reverse-strand transcript reads its reverse complement.
fn oriented_alt(alt: &[u8], strand: i8) -> Vec<u8> {
    if alt == b"-" {
        return Vec::new();
    }
    if strand < 0 {
        reverse_complement(alt)
    } else {
        alt.to_vec()
    }
}

/// Whether a structural variant's record has several alternates.
fn several_alternates(variant: &InputVariant) -> bool {
    variant
        .uploaded_allele_string
        .as_deref()
        .is_some_and(|alleles| alleles.contains('/'))
}

/// Whether a structural variant's record marks its place certain, with
/// neither the `IMPRECISE` flag nor a confidence interval around its position
/// or its end (`CIPOS`, `CIEND` with a bound other than 0).
fn states_its_place(variant: &InputVariant) -> bool {
    !variant.imprecise
        && [variant.ci_pos, variant.ci_end]
            .into_iter()
            .flatten()
            .all(|(low, high)| low == 0 && high == 0)
}

/// Whether a structural variant's record states one allele at one place: it
/// has a single alternate and marks its place certain ([`states_its_place`]).
fn states_one_place(variant: &InputVariant) -> bool {
    !several_alternates(variant) && states_its_place(variant)
}

/// Whether the alternate of `variant`, one of several of its record, states
/// the loss of the record's whole span: its own text is a deletion's as a
/// record with that one alternate is classed (`<DEL>`, `<CN0>`), and the
/// record writes `END` as the span's last base. A span read from a list of
/// lengths is the first alternate's, which need not be this one's.
fn loss_alternate_of_several(variant: &InputVariant) -> bool {
    let own_text_is_a_deletion = std::str::from_utf8(variant.alt_allele()).is_ok_and(|text| {
        matches!(
            classify_symbolic_alt(text, None),
            VariantClass::StructuralDeletion | VariantClass::MobileElementDeletion
        )
    });
    own_text_is_a_deletion && variant.written_end == WrittenEnd::SpanEnd
}

/// The number of reference bases a structural variant's record states as
/// deleted: the length of its stated span when it marks its place certain
/// ([`states_its_place`]) and either its class is a deletion or a
/// mobile-element deletion, or the row's alternate is the loss alternate of a
/// record with several ([`loss_alternate_of_several`]). `None` for a record
/// that states no span and for every other class and alternate.
fn stated_deletion(variant: &InputVariant) -> Option<usize> {
    let deletion = if several_alternates(variant) {
        loss_alternate_of_several(variant)
    } else {
        matches!(
            variant.variant_class,
            VariantClass::StructuralDeletion | VariantClass::MobileElementDeletion
        )
    };
    if !deletion || !states_its_place(variant) || variant.end < variant.start {
        return None;
    }
    usize::try_from(variant.end - variant.start + 1).ok()
}

/// The number of reference bases a structural variant's record states as
/// inverted: the length of its stated span when its class is an inversion
/// and it states one place ([`states_one_place`]). The record then gives the
/// mutated sequence whole, the span replaced by its reverse complement.
/// `None` for a record that states no span and for every other class.
fn stated_inversion(variant: &InputVariant) -> Option<usize> {
    if variant.variant_class != VariantClass::Inversion
        || !states_one_place(variant)
        || variant.end < variant.start
    {
        return None;
    }
    usize::try_from(variant.end - variant.start + 1).ok()
}

/// The number of bases a structural variant's record states as inserted after
/// its position: its stated length when its ALT is the plain `<INS>` and it
/// states one place ([`states_one_place`]). The record's span is then its
/// length in bases from the base after its position, whatever end it writes,
/// so the written end is read apart ([`WrittenEnd`]): none, the position, or
/// the span's own end. `None` for any other written end, which states an
/// interval beside the length, for a record without a length, for a mobile
/// element, whose length need not count the target-site duplication beside
/// it, and for every other insertion text.
fn stated_insertion(variant: &InputVariant) -> Option<usize> {
    let plain = variant.variant_class == VariantClass::StructuralInsertion
        && variant.alt_allele().eq_ignore_ascii_case(b"<INS>");
    if !plain || !states_one_place(variant) || variant.written_end == WrittenEnd::Other {
        return None;
    }
    // A length below zero states no count of inserted bases.
    let length = u64::try_from(variant.sv_len?).ok()?;
    let position = variant.start.checked_sub(1)?;
    if length == 0 || variant.end.checked_sub(position) != Some(length) {
        return None;
    }
    usize::try_from(length).ok()
}

/// The counts of `count_number_ATG` (`:1029-1070`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExistingCounts {
    /// Start codons with an in-frame stop codon inside the 5' UTR.
    uorfs: u32,
    /// Start codons without one whose frame differs from the coding sequence's.
    out_of_frame_oorfs: u32,
    /// Start codons without one in the coding sequence's frame.
    in_frame_oorfs: u32,
}

/// One transcript's 5' UTR as the plugin reads it: the genomic UTR segments in
/// ascending order, the UTR and coding sequences in transcript orientation, and
/// the reference scans every row of the transcript shares (`:227`, `:232`).
struct TranscriptUtr {
    strand: i8,
    /// The 5' UTR exon segments, ascending genomic order, every one of them.
    segments: Vec<(u64, u64)>,
    utr: Vec<u8>,
    /// The translateable sequence, `N`-padded by the start exon's phase.
    cds: Vec<u8>,
    scan: CodonScan,
    /// The first in-frame stop codon of each start codon over the UTR followed
    /// by the coding sequence (`:406-407`, `:815-816`, `:950-951`).
    utr_cds_first_stop: BTreeMap<usize, Option<usize>>,
    counts: ExistingCounts,
}

impl TranscriptUtr {
    /// Reads the transcript's UTR and coding sequences from the FASTA through its
    /// exon model; `None` when the transcript has no 5' UTR or a sequence is
    /// outside the reference.
    fn read(ctx: &LofteeContext, chr: &str, fasta: &vep_fasta::IndexedFasta) -> Option<Self> {
        let cds_start = usize::try_from(ctx.cdna_coding_start?).ok()?;
        let cds_end = usize::try_from(ctx.cdna_coding_end?).ok()?;
        if cds_start < 2 || cds_end < cds_start || ctx.exons.is_empty() {
            return None;
        }
        let strand = ctx.strand;
        let mut segments = Vec::new();
        let mut spliced = Vec::with_capacity(cds_end);
        let mut start_phase: i8 = 0;
        let mut offset = 0usize;
        for exon in &ctx.exons {
            let length = usize::try_from(exon.end.checked_sub(exon.start)? + 1).ok()?;
            let first = offset + 1;
            let last = offset + length;
            if first < cds_start {
                segments.push(exon_prefix_span(
                    exon,
                    strand,
                    last.min(cds_start - 1) - offset,
                ));
            }
            if first <= cds_start && cds_start <= last {
                start_phase = exon.phase;
            }
            if first <= cds_end {
                let need = last.min(cds_end) - offset;
                spliced.extend(exon_prefix(exon, strand, need, chr, fasta)?);
            }
            offset = last;
            if offset >= cds_end {
                break;
            }
        }
        if spliced.len() < cds_end {
            return None;
        }
        segments.sort_unstable();
        let utr = spliced[..cds_start - 1].to_vec();
        let mut cds = Vec::with_capacity(cds_end - cds_start + 1 + 2);
        if start_phase > 0 {
            cds.extend(std::iter::repeat_n(b'N', start_phase as usize));
        }
        cds.extend_from_slice(&spliced[cds_start - 1..cds_end]);
        let scan = CodonScan::new(&utr);
        let mut utr_cds = utr.clone();
        utr_cds.extend_from_slice(&cds);
        let utr_cds_stops = stop_positions(&utr_cds);
        let utr_cds_first_stop = scan
            .atgs
            .iter()
            .map(|&atg| (atg, first_in_frame_stop(&utr_cds_stops, atg)))
            .collect();
        let counts = existing_counts(&scan, utr.len());
        Some(Self {
            strand,
            segments,
            utr,
            cds,
            scan,
            utr_cds_first_stop,
            counts,
        })
    }

    /// The UTR segments the variant leaves under `max_overlap` (`:179-193`) with
    /// their position maps; `None` when none remains.
    fn kept_positions(&self, pos: u64, pos_end: u64, max_overlap: f64) -> Option<UtrPositions> {
        let kept: Vec<(u64, u64)> = self
            .segments
            .iter()
            .copied()
            .filter(|&(start, end)| {
                let length = (end - start + 1) as f64;
                let overlap_start = start.max(pos) as f64;
                let overlap_end = end.min(pos_end) as f64;
                let overlap = 100.0 * (overlap_end - overlap_start) / length;
                // `next if $UTR_overlap_perc >= $self->{max_overlap}` (`:188`): a
                // `nan` value compares false and skips no segment.
                overlap
                    .partial_cmp(&max_overlap)
                    .is_none_or(|order| order.is_lt())
            })
            .collect();
        (!kept.is_empty()).then_some(UtrPositions {
            strand: self.strand,
            segments: kept,
        })
    }
}

/// The genomic span of the first `need` bases of `exon` in transcript orientation.
fn exon_prefix_span(exon: &Exon, strand: i8, need: usize) -> (u64, u64) {
    let last = need as u64 - 1;
    if strand < 0 {
        (exon.end - last, exon.end)
    } else {
        (exon.start, exon.start + last)
    }
}

/// The first `need` bases of `exon` in transcript orientation, upper-cased as
/// Ensembl VEP reads a FASTA (`FastaSequence.pm:659`).
fn exon_prefix(
    exon: &Exon,
    strand: i8,
    need: usize,
    chr: &str,
    fasta: &vep_fasta::IndexedFasta,
) -> Option<Vec<u8>> {
    let need = need as u64;
    if strand < 0 {
        let seq = fasta.sequence(chr, exon.end + 1 - need, exon.end)?;
        Some(reverse_complement(&seq))
    } else {
        fasta.sequence(chr, exon.start, exon.start + need - 1)
    }
}

/// The kept UTR segments' two position maps (`utr_exon_position` `:1217-1249`
/// and `chr_position` `:1251-1283`): UTR indices count from 0 along the segments
/// in transcript order.
struct UtrPositions {
    strand: i8,
    /// Ascending genomic order.
    segments: Vec<(u64, u64)>,
}

impl UtrPositions {
    /// The UTR index of a genomic position inside a kept segment.
    fn index_of(&self, genomic: u64) -> Option<usize> {
        let mut offset = 0u64;
        if self.strand < 0 {
            for &(start, end) in self.segments.iter().rev() {
                if (start..=end).contains(&genomic) {
                    return usize::try_from(offset + (end - genomic)).ok();
                }
                offset += end - start + 1;
            }
        } else {
            for &(start, end) in &self.segments {
                if (start..=end).contains(&genomic) {
                    return usize::try_from(offset + (genomic - start)).ok();
                }
                offset += end - start + 1;
            }
        }
        None
    }

    /// The genomic position of a UTR index.
    fn genomic_of(&self, index: usize) -> Option<u64> {
        let mut index = index as u64;
        if self.strand < 0 {
            for &(start, end) in self.segments.iter().rev() {
                let length = end - start + 1;
                if index < length {
                    return Some(end - index);
                }
                index -= length;
            }
        } else {
            for &(start, end) in &self.segments {
                let length = end - start + 1;
                if index < length {
                    return Some(start + index);
                }
                index -= length;
            }
        }
        None
    }

    /// The UTR indices of the reference allele's 5' and 3' ends in transcript
    /// orientation (`get_allele_exon_pos` `:1285-1313`): on the forward strand
    /// the allele's start and end, on the reverse strand its end and start; an
    /// insertion's are both the index of the base after it in transcript order.
    fn allele_span(&self, strand: i8, pos: u64, ref_len: usize) -> Option<(usize, usize)> {
        let ref_len = ref_len as u64;
        if strand < 0 {
            // `$pos + length($ref_coding) - 1` is `$pos - 1` for an insertion.
            let start = self.index_of((pos + ref_len).checked_sub(1)?)?;
            let end = if ref_len > 0 {
                self.index_of(pos)?
            } else {
                start
            };
            Some((start, end))
        } else {
            let start = self.index_of(pos)?;
            let end = if ref_len > 0 {
                self.index_of(pos + ref_len - 1)?
            } else {
                start
            };
            Some((start, end))
        }
    }
}

/// The 5' UTR with the allele spliced in (`mut_utr_sequence` `:1189-1200`): the
/// bases before `mut_pos`, the alternate allele, the bases from
/// `mut_pos + ref_len` on.
fn mutated_utr(utr: &[u8], mut_pos: usize, ref_len: usize, alt_coding: &[u8]) -> Vec<u8> {
    let head = mut_pos.min(utr.len());
    let tail = (mut_pos + ref_len).min(utr.len());
    let mut out = Vec::with_capacity(head + alt_coding.len() + utr.len() - tail);
    out.extend_from_slice(&utr[..head]);
    out.extend_from_slice(alt_coding);
    out.extend_from_slice(&utr[tail..]);
    out
}

/// Every position where `ATG` starts (`get_ATG_pos` `:1113-1122`).
fn atg_positions(seq: &[u8]) -> Vec<usize> {
    seq.windows(3)
        .enumerate()
        .filter(|(_, codon)| *codon == b"ATG")
        .map(|(i, _)| i)
        .collect()
}

/// Every position where a stop codon starts (`get_stopcodon_pos` `:1135-1153`).
fn stop_positions(seq: &[u8]) -> Vec<usize> {
    seq.windows(3)
        .enumerate()
        .filter(|(_, codon)| matches!(*codon, b"TAA" | b"TAG" | b"TGA"))
        .map(|(i, _)| i)
        .collect()
}

/// The first stop codon at or after `atg` in its frame.
fn first_in_frame_stop(stops: &[usize], atg: usize) -> Option<usize> {
    stops
        .iter()
        .copied()
        .find(|&stop| stop >= atg && (stop - atg).is_multiple_of(3))
}

/// The start and stop codon positions of a sequence and the uORF map over them
/// (`existing_uORF` `:1083-1100`): each start codon with an in-frame stop codon
/// downstream, keyed to the first such stop.
struct CodonScan {
    atgs: Vec<usize>,
    /// Start codon to the first in-frame stop codon, for the start codons that
    /// have one.
    uorfs: BTreeMap<usize, usize>,
}

impl CodonScan {
    fn new(seq: &[u8]) -> Self {
        let atgs = atg_positions(seq);
        let stops = stop_positions(seq);
        let uorfs = atgs
            .iter()
            .filter_map(|&atg| first_in_frame_stop(&stops, atg).map(|stop| (atg, stop)))
            .collect();
        Self { atgs, uorfs }
    }
}

/// `count_number_ATG` (`:1029-1070`) over the reference UTR's scan.
fn existing_counts(scan: &CodonScan, utr_len: usize) -> ExistingCounts {
    let mut counts = ExistingCounts {
        uorfs: 0,
        out_of_frame_oorfs: 0,
        in_frame_oorfs: 0,
    };
    for &atg in &scan.atgs {
        if scan.uorfs.contains_key(&atg) {
            counts.uorfs += 1;
        } else if !(utr_len - atg).is_multiple_of(3) {
            counts.out_of_frame_oorfs += 1;
        } else {
            counts.in_frame_oorfs += 1;
        }
    }
    counts
}

/// The Kozak context and strength of the start codon at `pos_a` (`:368-387`): the
/// three bases before it, `ATG`, and the base after; `-` and an empty strength
/// when the sequence does not reach that far. Strong when position -3 is `A` or
/// `G` and position +4 is `G`, Moderate on one of the two, Weak otherwise
/// (`:375-384`).
fn kozak(seq: &[u8], pos_a: usize) -> (String, String) {
    if pos_a < 3 || pos_a + 3 >= seq.len() {
        return (NONE.to_string(), String::new());
    }
    let mut context = Vec::with_capacity(7);
    context.extend_from_slice(&seq[pos_a - 3..pos_a]);
    context.extend_from_slice(b"ATG");
    context.push(seq[pos_a + 3]);
    let minus_three = matches!(context[0], b'A' | b'G');
    let plus_four = context[6] == b'G';
    let strength = match (minus_three, plus_four) {
        (true, true) => "Strong",
        (false, false) => "Weak",
        _ => "Moderate",
    };
    (
        String::from_utf8_lossy(&context).into_owned(),
        strength.to_string(),
    )
}

/// The mutated UTR's window in which a gained start codon can sit: from two
/// bases before the allele to two after it, sized by the alternate allele
/// (`:345` sizes it by the reference allele through `:339`). Returns the window
/// and the index of its first base, which never wraps below 0 where `:345`'s
/// negative index would read the sequence's end.
fn atg_window(mutated: &[u8], mut_pos: usize, alt_len: usize) -> (&[u8], usize) {
    let first = mut_pos.saturating_sub(2);
    let last = (mut_pos + alt_len + 2).min(mutated.len());
    (&mutated[first.min(last)..last], first)
}

/// The first and the last UTR index at which an insertion that turns `utr`
/// into `mutated` can be written. Inserted bases that repeat the bases beside
/// them give the same mutated sequence at every index between the two: the
/// mutated sequence equals the reference up to the last index, and from the
/// first index plus the allele's length on it equals the reference from the
/// first index on.
fn insertion_writings(utr: &[u8], mutated: &[u8]) -> (usize, usize) {
    let same_before = utr
        .iter()
        .zip(mutated)
        .take_while(|(reference, changed)| reference == changed)
        .count();
    let same_after = utr
        .iter()
        .rev()
        .zip(mutated.iter().rev())
        .take_while(|(reference, changed)| reference == changed)
        .count();
    (utr.len() - same_after, same_before)
}

/// The last UTR index a variant may touch to count as inside the uORF at
/// `start_pos` (`:489-500`, `:904-917`): the base before its first in-frame stop
/// codon, or the UTR's last base for an overlapping ORF.
fn check_point(model: &TranscriptUtr, start_pos: usize) -> usize {
    match model.scan.uorfs.get(&start_pos) {
        Some(&stop) => stop - 1,
        None => model.utr.len() - 1,
    }
}

/// The uORF type of the start codon at `start_pos` in a UTR of `utr_len` bases
/// whose codons `scan` holds: a uORF when an in-frame stop codon follows inside
/// the UTR, else an overlapping ORF out of or in the coding sequence's frame.
fn orf_type(scan: &CodonScan, utr_len: usize, start_pos: usize) -> &'static str {
    if scan.uorfs.contains_key(&start_pos) {
        "uORF"
    } else if !(utr_len - start_pos).is_multiple_of(3) {
        "OutOfFrame_oORF"
    } else {
        "InFrame_oORF"
    }
}

/// The reference uORF type of the start codon at `start_pos` (`:539-545`,
/// `:805-811`, `:953-959`).
fn reference_type(model: &TranscriptUtr, start_pos: usize) -> &'static str {
    orf_type(&model.scan, model.utr.len(), start_pos)
}

/// The keys of a frameshift's effect on the ORF at `start_pos` that read
/// nothing of the allele (`:928-968`, `:990`): the ORF's type, its length
/// over the reference UTR and the coding sequence after it, its start's
/// distance to the coding sequence, its Kozak context and strength, and its
/// evidence from the table.
fn reference_frame_keys(
    plugin: &UtrAnnotatorPlugin,
    model: &TranscriptUtr,
    chr: &str,
    positions: &UtrPositions,
    start_pos: usize,
) -> Vec<(&'static str, String)> {
    let (context, strength) = kozak(&model.utr, start_pos);
    let ref_type_length = model
        .utr_cds_first_stop
        .get(&start_pos)
        .copied()
        .flatten()
        .map_or(NA.to_string(), |stop| (stop - start_pos + 3).to_string());
    vec![
        ("ref_type", reference_type(model, start_pos).to_string()),
        ("ref_type_length", ref_type_length),
        (
            "ref_StartDistanceToCDS",
            (model.utr.len() - start_pos).to_string(),
        ),
        ("KozakContext", context),
        ("KozakStrength", strength),
        ("Evidence", plugin.evidence(chr, positions, start_pos)),
    ]
}

/// Whether a codon read from the mutated UTR is a stop codon; a codon the
/// sequence does not reach (`:639` reads undefined bases) is not.
fn is_stop(codon: Option<&[u8]>) -> bool {
    matches!(codon, Some(b"TAA" | b"TAG" | b"TGA"))
}

/// One allele spliced into a transcript's 5' UTR: what the five analyses read.
struct MutatedUtr<'a> {
    plugin: &'a UtrAnnotatorPlugin,
    chr: &'a str,
    positions: &'a UtrPositions,
    model: &'a TranscriptUtr,
    /// The 5' UTR with the allele spliced in ([`mutated_utr`]).
    mutated: Vec<u8>,
    /// The start and stop codons of `mutated`.
    scan: CodonScan,
    /// The UTR indices of the reference allele's 5' and 3' ends
    /// ([`UtrPositions::allele_span`]).
    mut_pos: usize,
    end_pos: usize,
    /// The lengths of the reference and alternate alleles, `-` counting as 0.
    ref_len: usize,
    alt_len: usize,
    /// The first and the last UTR index at which an insertion's mutated
    /// sequence can be written ([`insertion_writings`]); `None` for any other
    /// allele.
    writings: Option<(usize, usize)>,
    /// `writings` when one of them lies at one of a reference start codon's
    /// three indices or directly after the codon. The stop-codon and
    /// frameshift analyses read such an insertion from all its writings; any
    /// other allele they read at the record's own index.
    beside_start_codon: Option<(usize, usize)>,
}

impl MutatedUtr<'_> {
    /// The evidence for the uORF starting at `start_pos` of the reference UTR
    /// ([`UtrAnnotatorPlugin::evidence`]).
    fn evidence(&self, start_pos: usize) -> String {
        self.plugin.evidence(self.chr, self.positions, start_pos)
    }

    /// Every position where a stop codon starts in the mutated UTR followed by
    /// the coding sequence.
    fn utr_cds_stops(&self) -> Vec<usize> {
        let mut utr_cds = self.mutated.clone();
        utr_cds.extend_from_slice(&self.model.cds);
        stop_positions(&utr_cds)
    }

    /// Whether an insertion with the writings `first..=last` loses the stop
    /// codon at `stop_pos`. One writing lies inside the stop codon and the
    /// mutated UTR reads no stop codon at that index (`:628-642` read those
    /// three bases only for a record written inside the stop codon). When a
    /// writing also lies at or before the stop codon's first base, the stop
    /// codon stands whole as many bases on as the allele is long, so it is
    /// lost only at a length that is no multiple of three.
    fn loses_stop_codon(&self, stop_pos: usize, (first, last): (usize, usize)) -> bool {
        last > stop_pos
            && first <= stop_pos + 2
            && !is_stop(self.mutated.get(stop_pos..stop_pos + 3))
            && (first > stop_pos || !self.alt_len.is_multiple_of(3))
    }

    /// Whether the uORF at `start_pos` keeps its length under an insertion
    /// whose last writing is `last`: the mutated UTR reads a stop codon at the
    /// index of the uORF's own, and the insertion can be written at or after
    /// that index, or directly before a start codon at whose second base that
    /// stop codon begins. Every codon of the uORF before its last then stands
    /// unchanged and the uORF ends where it ended: it is not shifted. `:921`
    /// holds the record's own index against the check point and reads no base
    /// of the mutated sequence, so it lists such a uORF for a record written
    /// before the stop codon.
    fn keeps_its_length(&self, start_pos: usize, last: usize) -> bool {
        let Some(&stop_pos) = self.model.scan.uorfs.get(&start_pos) else {
            return false;
        };
        is_stop(self.mutated.get(stop_pos..stop_pos + 3))
            && (last >= stop_pos
                || (last + 1 == stop_pos
                    && self.model.utr.get(last..last + 3) == Some(b"ATG".as_slice())))
    }

    /// `uAUG_gained` (`:307-435`): a start codon created by the allele, the
    /// first of the window. An insertion's window is that of its first
    /// writing, so the codon named does not depend on the index the record
    /// gives: of two equal copies of a repeated start codon it is the one
    /// nearer the cap. A reference start codon that stands at its own index
    /// before every writing is not created and is passed over, where
    /// `:346-353` name the first start codon of the window whatever the
    /// reference holds.
    fn uaug_gained(&self) -> Option<Analysis> {
        let window_at = self.writings.map_or(self.mut_pos, |(first, _)| first);
        let (window, window_first) = atg_window(&self.mutated, window_at, self.alt_len);
        let pos_a = atg_positions(window)
            .into_iter()
            .map(|at| window_first + at)
            .find(|&at| {
                !self
                    .writings
                    .is_some_and(|(first, last)| at < first && at + 3 <= last)
            })?;
        let (context, strength) = kozak(&self.mutated, pos_a);
        // `inFrame_oORF` here, where the other analyses write `InFrame_oORF`.
        let kind = if self.scan.uorfs.contains_key(&pos_a) {
            "uORF"
        } else if !(self.mutated.len() - pos_a).is_multiple_of(3) {
            "OutOfFrame_oORF"
        } else {
            "inFrame_oORF"
        };
        let distance_to_stop = first_in_frame_stop(&self.utr_cds_stops(), pos_a)
            .map_or(NA.to_string(), |stop| (stop - pos_a).to_string());
        Some(Analysis {
            flag: UAUG_GAINED_FLAG,
            effects: vec![vec![
                ("KozakContext", context),
                ("KozakStrength", strength),
                ("DistanceToCDS", (self.mutated.len() - pos_a).to_string()),
                ("type", kind.to_string()),
                ("DistanceToStop", distance_to_stop),
                ("CapDistanceToStart", pos_a.to_string()),
            ]],
        })
    }

    /// `uSTOP_gained` (`:437-574`): a substitution inside an existing uORF or
    /// overlapping ORF creating a stop codon before its reference stop.
    fn ustop_gained(&self) -> Option<Analysis> {
        if self.ref_len != self.alt_len {
            return None;
        }
        let mut effects = Effects::new();
        for &start_pos in &self.model.scan.atgs {
            let check_point = check_point(self.model, start_pos);
            if self.mut_pos < start_pos + 3 || self.end_pos > check_point {
                continue;
            }
            let Some(&mut_stop) = self.scan.uorfs.get(&start_pos) else {
                continue;
            };
            if mut_stop >= check_point {
                continue;
            }
            let (context, strength) = kozak(&self.model.utr, start_pos);
            effects.push(vec![
                (
                    "ref_type",
                    reference_type(self.model, start_pos).to_string(),
                ),
                (
                    "ref_StartDistanceToCDS",
                    (self.model.utr.len() - start_pos).to_string(),
                ),
                (
                    "newSTOPDistanceToCDS",
                    (self.mutated.len() - mut_stop).to_string(),
                ),
                ("KozakContext", context),
                ("KozakStrength", strength),
                ("Evidence", self.evidence(start_pos)),
            ]);
        }
        (!effects.is_empty()).then_some(Analysis {
            flag: USTOP_GAINED_FLAG,
            effects,
        })
    }

    /// `uSTOP_lost` (`:576-715`): the allele removes the stop codon of an existing
    /// uORF. An insertion at or beside a start codon is read from all its
    /// writings ([`Self::loses_stop_codon`]), whichever of them the record
    /// gives.
    fn ustop_lost(&self) -> Option<Analysis> {
        let mut effects = Effects::new();
        for (&start_pos, &stop_pos) in &self.model.scan.uorfs {
            let lost = match self.beside_start_codon {
                Some(writings) => self.loses_stop_codon(stop_pos, writings),
                // `$mut_pos + length($ref_coding) - 1 < $stop_pos` (`:631`); the
                // insertion test at `:633` can never hold.
                None => {
                    self.mut_pos <= stop_pos + 2
                        && self.mut_pos + self.ref_len > stop_pos
                        && !(self.alt_len != 0 && is_stop(self.mutated.get(stop_pos..stop_pos + 3)))
                }
            };
            if !lost {
                continue;
            }
            let (context, strength) = kozak(&self.model.utr, start_pos);
            let (alt_stop, alt_stop_distance) = match self.scan.uorfs.get(&start_pos) {
                Some(&stop) => ("True", (self.mutated.len() - stop).to_string()),
                None => ("False", NA.to_string()),
            };
            // `:685`: a deletion reaching back over the start codon leaves the
            // start past the mutated UTR's end, and Perl's `%` of that negative
            // distance is its remainder in `0..3`.
            let distance = self.mutated.len() as i64 - start_pos as i64;
            let frame = if distance.rem_euclid(3) != 0 {
                "outOfFrame"
            } else {
                "inFrame"
            };
            effects.push(vec![
                ("AltStop", alt_stop.to_string()),
                ("AltStopDistanceToCDS", alt_stop_distance),
                ("FrameWithCDS", frame.to_string()),
                ("KozakContext", context),
                ("KozakStrength", strength),
                ("Evidence", self.evidence(start_pos)),
            ]);
        }
        (!effects.is_empty()).then_some(Analysis {
            flag: USTOP_LOST_FLAG,
            effects,
        })
    }

    /// `uAUG_lost` (`:717-848`): the allele starts inside an existing start codon
    /// and removes it. An insertion removes the codon only when every one of
    /// its writings parts the codon's three bases. One that can be written
    /// directly before the codon moves it on whole, and one that can be
    /// written directly after it leaves it at its own index: neither is
    /// reported, where `:764-769` pass an insertion at the codon's three
    /// indices and `:776-778` read the three bases at the codon's old index
    /// alone.
    fn uaug_lost(&self) -> Option<Analysis> {
        let mut effects = Effects::new();
        for &start_pos in &self.model.scan.atgs {
            let lost = match self.writings {
                Some((first, last)) => first > start_pos && last < start_pos + 3,
                // `:764-769`: the allele starts at one of the codon's three bases (an
                // allele starting upstream of the codon is skipped even when it covers
                // it, and `:767` can then never hold).
                None => {
                    (start_pos..=start_pos + 2).contains(&self.mut_pos)
                        && !(self.alt_len != 0
                            && self.mutated.get(start_pos..start_pos + 3)
                                == Some(b"ATG".as_slice()))
                }
            };
            if !lost {
                continue;
            }
            let (context, strength) = kozak(&self.model.utr, start_pos);
            let distance_to_stop = self
                .model
                .utr_cds_first_stop
                .get(&start_pos)
                .copied()
                .flatten()
                .map_or(NA.to_string(), |stop| (stop - start_pos).to_string());
            effects.push(vec![
                ("type", reference_type(self.model, start_pos).to_string()),
                ("CapDistanceToStart", start_pos.to_string()),
                (
                    "DistanceToCDS",
                    (self.model.utr.len() - start_pos).to_string(),
                ),
                ("DistanceToStop", distance_to_stop),
                ("KozakContext", context),
                ("KozakStrength", strength),
                ("Evidence", self.evidence(start_pos)),
            ]);
        }
        (!effects.is_empty()).then_some(Analysis {
            flag: UAUG_LOST_FLAG,
            effects,
        })
    }

    /// `uFrameshift` (`:850-1018`): an indel inside an existing uORF or overlapping
    /// ORF whose length is not a multiple of three. An insertion at or beside
    /// a start codon shifts an ORF when one of its writings lies inside the
    /// ORF, after a start codon that stands at its own index and at or before
    /// the check point, unless the ORF keeps its length
    /// ([`Self::keeps_its_length`]). An ORF whose start codon the insertion
    /// moves on whole stands unchanged after the inserted bases. `:921` reads
    /// the record's own index: it leaves out the ORF of a start codon the
    /// record is written inside, and an ORF whose stop codon begins directly
    /// after the record's index.
    fn uframeshift(&self) -> Option<Analysis> {
        if self.ref_len == self.alt_len || self.ref_len.abs_diff(self.alt_len).is_multiple_of(3) {
            return None;
        }
        let mut alt_overlapping_stops: Option<Vec<usize>> = None;
        let mut effects = Effects::new();
        for &start_pos in &self.model.scan.atgs {
            let check_point = check_point(self.model, start_pos);
            let shifted = match self.beside_start_codon {
                Some((first, last)) => {
                    first > start_pos
                        && last >= start_pos + 3
                        && start_pos + 3 <= check_point
                        && first <= check_point
                        && !self.keeps_its_length(start_pos, last)
                }
                None => self.mut_pos >= start_pos + 3 && self.end_pos <= check_point,
            };
            if !shifted {
                continue;
            }
            let alt_type = orf_type(&self.scan, self.mutated.len(), start_pos);
            let alt_stops = alt_overlapping_stops.get_or_insert_with(|| self.utr_cds_stops());
            let alt_type_length = first_in_frame_stop(alt_stops, start_pos)
                .map_or(NA.to_string(), |stop| (stop - start_pos + 3).to_string());
            let mut effect =
                reference_frame_keys(self.plugin, self.model, self.chr, self.positions, start_pos);
            effect.push(("alt_type", alt_type.to_string()));
            effect.push(("alt_type_length", alt_type_length));
            effects.push(effect);
        }
        (!effects.is_empty()).then_some(Analysis {
            flag: UFRAMESHIFT_FLAG,
            effects,
        })
    }
}

/// The text form of an analysis's effects (`:280-288`): each effect's
/// `key=value` pairs in byte order of their keys (capitalised keys first) joined
/// by `:`, the effects joined by `,`.
fn render_effects(effects: &Effects) -> String {
    effects
        .iter()
        .map(|effect| {
            let mut pairs: Vec<&(&str, String)> = effect.iter().collect();
            pairs.sort_unstable_by(|a, b| a.0.cmp(b.0));
            pairs
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<_>>()
                .join(":")
        })
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::ConsequenceList;
    use vep_core::transcript::Transcript;
    use vep_core::variant::VariantClass;

    /// A 50-base 5' UTR with a uORF (`ATG` at 5, `TGA` at 20), an out-of-frame
    /// overlapping ORF (`ATG` at 30) and an in-frame one (`ATG` at 38); `ACG` at
    /// 24 and `TGC` at 14 are one substitution away from a start and a stop codon.
    const UTR: &str = "CCGCCATGGCCCCCTGCCCCTGACACGCCCATGCCACCATGTCCCCCCCC";
    /// A 90-base coding sequence: `ATG`, 28 `GCC` codons, `TAA`.
    const CDS: &str = "ATGGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTAA";

    /// A 50-base 5' UTR whose start codons overlap stop codons. The uORF of the
    /// `ATG` at 5 ends with the `TGA` at 17, which begins inside the `ATG` at 16;
    /// the uORF of that `ATG` ends with the `TGA` at 25, whose last base opens the
    /// `ATG` at 27.
    const UTR_OVERLAPS: &str = "CCGCCATGGCCCCCCCATGACCCCCTGATGCCCCCCCCCCCCCCCCCCCC";

    /// Where the fixture places the transcript on a 400-base contig `21`. On the
    /// forward strand a 40-base UTR exon at 101-140, a 20-base intron, then an
    /// exon at 161-260 holding the last 10 UTR bases and the coding sequence; on
    /// the reverse strand the same transcript read downwards, its first exon at
    /// 261-300 and its second at 141-240.
    const UTR_SPLIT: usize = 40;

    /// The two exons in transcript order (rank 1 first) and the coding region's
    /// genomic bounds.
    fn layout(strand: i8) -> ([(u64, u64); 2], (u64, u64)) {
        if strand < 0 {
            ([(261, 300), (141, 240)], (141, 230))
        } else {
            ([(101, 140), (161, 260)], (171, 260))
        }
    }

    /// The UTR exon segments of [`layout`], ascending.
    fn utr_segments(strand: i8) -> Vec<(u64, u64)> {
        if strand < 0 {
            vec![(231, 240), (261, 300)]
        } else {
            vec![(101, 140), (161, 170)]
        }
    }

    struct Fixture {
        _dir: tempfile::TempDir,
        plugin: UtrAnnotatorPlugin,
        strand: i8,
    }

    impl Fixture {
        fn new(strand: i8, params: &[&str]) -> Self {
            Self::with_utr(strand, UTR, params)
        }

        /// The fixture with `utr` as the transcript's 5' UTR, as long as [`UTR`].
        fn with_utr(strand: i8, utr: &str, params: &[&str]) -> Self {
            assert_eq!(utr.len(), UTR.len());
            let mut genome = vec![b'C'; 400];
            let mut transcript: Vec<u8> = Vec::new();
            transcript.extend_from_slice(utr.as_bytes());
            transcript.extend_from_slice(CDS.as_bytes());
            let ([rank1, rank2], _) = layout(strand);
            let place = |genome: &mut Vec<u8>, span: (u64, u64), bases: &[u8]| {
                let start = span.0 as usize - 1;
                assert_eq!(bases.len(), (span.1 - span.0 + 1) as usize);
                genome[start..start + bases.len()].copy_from_slice(bases);
            };
            if strand < 0 {
                place(
                    &mut genome,
                    rank1,
                    &reverse_complement(&transcript[..UTR_SPLIT]),
                );
                place(
                    &mut genome,
                    rank2,
                    &reverse_complement(&transcript[UTR_SPLIT..]),
                );
            } else {
                place(&mut genome, rank1, &transcript[..UTR_SPLIT]);
                place(&mut genome, rank2, &transcript[UTR_SPLIT..]);
            }
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("ref.fa");
            let mut contents = b">21\n".to_vec();
            contents.extend_from_slice(&genome);
            contents.push(b'\n');
            std::fs::write(&fa, contents).unwrap();
            std::fs::write(dir.path().join("ref.fa.fai"), "21\t400\t4\t400\t401\n").unwrap();
            let mut plugin = UtrAnnotatorPlugin::new();
            let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
            plugin.init(&params).unwrap();
            plugin.set_reference_fasta(Some(Arc::new(
                vep_fasta::IndexedFasta::from_path(&fa).unwrap(),
            )));
            Self {
                _dir: dir,
                plugin,
                strand,
            }
        }

        fn context(&self) -> LofteeContext {
            let ([rank1, rank2], coding) = layout(self.strand);
            let exon = |span: (u64, u64), rank: u32, end_phase: i8| Exon {
                start: span.0,
                end: span.1,
                rank,
                phase: -1,
                end_phase,
            };
            let mut context =
                LofteeContext::from_transcript(&vep_core::transcript::Transcript::default());
            context.exons = vec![exon(rank1, 1, -1), exon(rank2, 2, 0)];
            context.introns = Vec::new();
            context.strand = self.strand;
            context.coding_region_start = Some(coding.0);
            context.coding_region_end = Some(coding.1);
            context.translation_start = None;
            context.translation_end = None;
            context.cdna_coding_start = Some(UTR.len() as u64 + 1);
            context.cdna_coding_end = Some((UTR.len() + CDS.len()) as u64);
            context.cds_start_nf = false;
            context.cds_end_nf = false;
            context
        }

        /// The genomic position of UTR index `index` in this fixture.
        fn genomic(&self, index: usize) -> u64 {
            let positions = UtrPositions {
                strand: self.strand,
                segments: utr_segments(self.strand),
            };
            positions.genomic_of(index).unwrap()
        }

        /// The genomic span of UTR indices `first..=last`, ascending.
        fn span(&self, first: usize, last: usize) -> (u64, u64) {
            let (a, b) = (self.genomic(first), self.genomic(last));
            (a.min(b), a.max(b))
        }

        /// Runs the plugin on one 5' UTR row of the transcript for the allele
        /// `reference/alternate` (transcript orientation, `-` for none) whose
        /// first base in transcript order sits at UTR index `index`; an
        /// insertion sits before that index.
        fn annotate(
            &self,
            index: usize,
            reference: &str,
            alternate: &str,
            consequences: &[Consequence],
        ) -> IndexMap<String, String> {
            let mut variant = self.variant(index, reference, alternate, consequences);
            self.plugin
                .run_batch(std::slice::from_mut(&mut variant))
                .unwrap();
            variant.transcript_consequences.remove(0).plugin_data
        }

        /// The record of the allele `reference/alternate` with one 5' UTR row
        /// of the transcript, placed as [`Self::annotate`] says.
        fn variant(
            &self,
            index: usize,
            reference: &str,
            alternate: &str,
            consequences: &[Consequence],
        ) -> InputVariant {
            let ref_len = (if reference == "-" { 0 } else { reference.len() }) as u64;
            let (start, end, ref_allele, alt_allele) = if self.strand < 0 {
                let five_prime = self.genomic(index);
                let rc = |s: &str| {
                    if s == "-" {
                        b"-".to_vec()
                    } else {
                        reverse_complement(s.as_bytes())
                    }
                };
                if ref_len == 0 {
                    (five_prime + 1, five_prime, rc(reference), rc(alternate))
                } else {
                    (
                        five_prime + 1 - ref_len,
                        five_prime,
                        rc(reference),
                        rc(alternate),
                    )
                }
            } else {
                let start = self.genomic(index);
                if ref_len == 0 {
                    (
                        start,
                        start - 1,
                        b"-".to_vec(),
                        alternate.as_bytes().to_vec(),
                    )
                } else {
                    (
                        start,
                        start + ref_len - 1,
                        reference.as_bytes().to_vec(),
                        alternate.as_bytes().to_vec(),
                    )
                }
            };
            let mut variant = InputVariant::new("21".into(), start, end, ref_allele, alt_allele);
            variant.transcript_consequences.push({
                let mut consequence = TranscriptConsequence::default();
                consequence.transcript = Some(Arc::new({
                    let mut transcript = Transcript::default();
                    transcript.stable_id = Arc::from("ENST00000000001");
                    transcript
                }));
                consequence.consequences = ConsequenceList::from_slice(consequences);
                consequence.loftee_ctx = Some(Box::new(self.context()));
                consequence
            });
            variant
        }

        /// The 5' UTR rows of many insertions, each the bases written before
        /// its UTR index, annotated as one batch.
        fn insertion_rows(&self, insertions: &[(usize, String)]) -> Vec<IndexMap<String, String>> {
            let mut variants: Vec<InputVariant> = insertions
                .iter()
                .map(|(index, bases)| {
                    self.variant(*index, "-", bases, &[Consequence::FivePrimeUtrVariant])
                })
                .collect();
            self.plugin.run_batch(&mut variants).unwrap();
            variants
                .into_iter()
                .map(|mut variant| variant.transcript_consequences.remove(0).plugin_data)
                .collect()
        }
    }

    fn utr_row(
        fixture: &Fixture,
        index: usize,
        reference: &str,
        alternate: &str,
    ) -> IndexMap<String, String> {
        fixture.annotate(
            index,
            reference,
            alternate,
            &[Consequence::FivePrimeUtrVariant],
        )
    }

    /// One 5' UTR row of the transcript for a symbolic allele of `class` whose
    /// record states the genomic span `span`; `record` sets what else the
    /// record states.
    fn symbolic_row(
        fixture: &Fixture,
        class: VariantClass,
        span: (u64, u64),
        record: impl FnOnce(&mut InputVariant),
    ) -> IndexMap<String, String> {
        let mut variant = InputVariant::new(
            "21".into(),
            span.0,
            span.1,
            b"N".to_vec(),
            b"<DEL>".to_vec(),
        );
        variant.is_structural = true;
        variant.variant_class = class;
        variant.sv_end = Some(span.1);
        record(&mut variant);
        variant.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.consequences =
                ConsequenceList::from_slice(&[Consequence::FivePrimeUtrVariant]);
            consequence.loftee_ctx = Some(Box::new(fixture.context()));
            consequence
        });
        fixture
            .plugin
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();
        variant.transcript_consequences.remove(0).plugin_data
    }

    /// The row of a symbolic deletion of UTR indices `first..=last`.
    fn symbolic_deletion(fixture: &Fixture, first: usize, last: usize) -> IndexMap<String, String> {
        symbolic_row(
            fixture,
            VariantClass::StructuralDeletion,
            fixture.span(first, last),
            |_| {},
        )
    }

    /// The row of a plain symbolic insertion of `length` bases before UTR index
    /// `index`, handed as a structural record is: a span of `length` bases from
    /// the base after the record's position. `record` sets what else the record
    /// states.
    fn symbolic_insertion(
        fixture: &Fixture,
        index: usize,
        length: u64,
        record: impl FnOnce(&mut InputVariant),
    ) -> IndexMap<String, String> {
        // The base before the insertion in genomic order.
        let position = if fixture.strand < 0 {
            fixture.genomic(index)
        } else {
            fixture.genomic(index) - 1
        };
        symbolic_row(
            fixture,
            VariantClass::StructuralInsertion,
            (position + 1, position + length),
            |variant| {
                variant.alt_alleles = vec![b"<INS>".to_vec()];
                variant.sv_len = Some(length as i64);
                record(variant);
            },
        )
    }

    /// [`UTR`] with a stop codon at index 43, in no start codon's frame: moved
    /// by one base it ends the frames of the start codons at 5 and 38.
    fn utr_with_late_stop() -> String {
        let mut utr = UTR.to_string();
        utr.replace_range(43..46, "TAA");
        utr
    }

    /// `annotation` without the pairs of `keys`.
    fn without_keys(annotation: &str, keys: &[&str]) -> String {
        annotation
            .split(',')
            .map(|effect| {
                effect
                    .split(':')
                    .filter(|pair| !keys.iter().any(|key| pair.split('=').next() == Some(*key)))
                    .collect::<Vec<_>>()
                    .join(":")
            })
            .collect::<Vec<_>>()
            .join(",")
    }

    /// `utr` with `bases` in the place of the bases at `at`.
    fn utr_with(utr: &str, at: std::ops::Range<usize>, bases: &str) -> String {
        let mut out = utr.to_string();
        out.replace_range(at, bases);
        out
    }

    /// A whole row: the term, its annotation and the reference UTR's counts
    /// of uORFs, out-of-frame and in-frame overlapping ORFs.
    fn whole_row(
        flag: &str,
        annotation: &str,
        (uorfs, out_of_frame, in_frame): (&str, &str, &str),
    ) -> IndexMap<String, String> {
        IndexMap::from([
            (CONSEQUENCE_KEY.to_string(), flag.to_string()),
            (ANNOTATION_KEY.to_string(), annotation.to_string()),
            (EXISTING_UORFS_KEY.to_string(), uorfs.to_string()),
            (
                EXISTING_OUT_OF_FRAME_KEY.to_string(),
                out_of_frame.to_string(),
            ),
            (EXISTING_IN_FRAME_KEY.to_string(), in_frame.to_string()),
        ])
    }

    fn counts_of(out: &IndexMap<String, String>) -> (&str, &str, &str) {
        (
            out[EXISTING_UORFS_KEY].as_str(),
            out[EXISTING_OUT_OF_FRAME_KEY].as_str(),
            out[EXISTING_IN_FRAME_KEY].as_str(),
        )
    }

    #[test]
    fn only_five_prime_utr_rows_are_evaluated() {
        let fixture = Fixture::new(1, &[]);
        let utr = fixture.annotate(12, "C", "T", &[Consequence::FivePrimeUtrVariant]);
        assert_eq!(utr.len(), 5, "{utr:?}");
        let coding = fixture.annotate(12, "C", "T", &[Consequence::MissenseVariant]);
        assert!(coding.is_empty(), "{coding:?}");
        let both = fixture.annotate(
            12,
            "C",
            "T",
            &[
                Consequence::SpliceRegionVariant,
                Consequence::FivePrimeUtrVariant,
            ],
        );
        assert_eq!(both.len(), 5);
    }

    /// The fixture UTR by hand: `ATG` at 5 with the in-frame `TGA` at 20 is a
    /// uORF; `ATG` at 30 has 20 bases to the coding sequence, out of frame;
    /// `ATG` at 38 has 12, in frame.
    #[test]
    fn existing_counts_match_a_hand_counted_utr() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let out = utr_row(&fixture, 12, "C", "T");
            assert_eq!(counts_of(&out), ("1", "1", "1"), "strand {strand}: {out:?}");
            assert_eq!(out[CONSEQUENCE_KEY], NONE);
            assert_eq!(out[ANNOTATION_KEY], NONE);
        }
    }

    /// `ACG` at 24 becomes `ATG`: 24 bases from the cap, 26 to the coding
    /// sequence, no in-frame stop in the UTR or the coding sequence, context
    /// `GAC ATG C` (`-3` is `G`, `+4` is `C`: Moderate). `TCC` at 41 becomes
    /// `ATG`: 9 bases to the coding sequence, in its frame, so the new ORF runs
    /// to the coding stop codon 96 bases on and its type is spelled
    /// `inFrame_oORF`, the one place the plugin writes the lowercase form.
    #[test]
    fn uaug_gained_reports_kozak_context_and_strength() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let out = utr_row(&fixture, 25, "C", "T");
            assert_eq!(out[CONSEQUENCE_KEY], UAUG_GAINED_FLAG, "strand {strand}");
            assert_eq!(
                out[ANNOTATION_KEY],
                "CapDistanceToStart=24:DistanceToCDS=26:DistanceToStop=NA:KozakContext=GACATGC:KozakStrength=Moderate:type=OutOfFrame_oORF",
                "strand {strand}"
            );
            let in_frame = utr_row(&fixture, 41, "TCC", "ATG");
            assert_eq!(
                in_frame[CONSEQUENCE_KEY], UAUG_GAINED_FLAG,
                "strand {strand}"
            );
            assert_eq!(
                in_frame[ANNOTATION_KEY],
                "CapDistanceToStart=41:DistanceToCDS=9:DistanceToStop=96:KozakContext=ATGATGC:KozakStrength=Moderate:type=inFrame_oORF",
                "strand {strand}"
            );
        }
    }

    #[test]
    fn kozak_strength_classes_weak_moderate_strong() {
        // -3 and +4 both favourable.
        assert_eq!(
            kozak(b"GCCATGGCC", 3),
            ("GCCATGG".to_string(), "Strong".to_string())
        );
        assert_eq!(
            kozak(b"ACCATGGCC", 3),
            ("ACCATGG".to_string(), "Strong".to_string())
        );
        // One of the two.
        assert_eq!(
            kozak(b"GCCATGCCC", 3),
            ("GCCATGC".to_string(), "Moderate".to_string())
        );
        assert_eq!(
            kozak(b"CCCATGGCC", 3),
            ("CCCATGG".to_string(), "Moderate".to_string())
        );
        // Neither.
        assert_eq!(
            kozak(b"CCCATGCCC", 3),
            ("CCCATGC".to_string(), "Weak".to_string())
        );
        assert_eq!(
            kozak(b"TCCATGACC", 3),
            ("TCCATGA".to_string(), "Weak".to_string())
        );
        // Too close to either end of the sequence: `-` and an empty strength.
        assert_eq!(kozak(b"CCATGGCC", 2), (NONE.to_string(), String::new()));
        assert_eq!(kozak(b"GCCATG", 3), (NONE.to_string(), String::new()));
    }

    /// `TGC` at 14 becomes `TGA`, in frame with the uORF at 5 (Strong context
    /// `GCC ATG G`) and before its reference stop at 20: 45 bases from the start
    /// to the coding sequence, 36 from the new stop. A deletion of the same base
    /// is a frameshift, never a gained stop.
    #[test]
    fn ustop_gained_is_snv_and_mnv_only() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let snv = utr_row(&fixture, 16, "C", "A");
            assert_eq!(snv[CONSEQUENCE_KEY], USTOP_GAINED_FLAG, "strand {strand}");
            assert_eq!(
                snv[ANNOTATION_KEY],
                "Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:newSTOPDistanceToCDS=36:ref_StartDistanceToCDS=45:ref_type=uORF"
            );
            let mnv = utr_row(&fixture, 15, "GC", "AA");
            assert_eq!(mnv[CONSEQUENCE_KEY], USTOP_GAINED_FLAG, "strand {strand}");
            let deletion = utr_row(&fixture, 16, "C", "-");
            assert_eq!(
                deletion[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                "strand {strand}"
            );
        }
    }

    /// `TGA` at 20 becomes `TGC`: the uORF at 5 loses its stop with no other
    /// in-frame stop in the UTR (45 bases to the coding sequence, in frame). With
    /// a second in-frame stop downstream, the alternative stop is reported with
    /// its distance.
    #[test]
    fn ustop_lost_reports_alt_stop_and_frame() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let out = utr_row(&fixture, 22, "A", "C");
            assert_eq!(out[CONSEQUENCE_KEY], USTOP_LOST_FLAG, "strand {strand}");
            assert_eq!(
                out[ANNOTATION_KEY],
                "AltStop=False:AltStopDistanceToCDS=NA:Evidence=NA:FrameWithCDS=inFrame:KozakContext=GCCATGG:KozakStrength=Strong"
            );
        }
        // A UTR of 30 bases: `ATG` at 3, `TAA` at 12 and `TGA` at 24, both in
        // frame; `TAA` to `TAC` leaves the second stop, 6 bases from the coding
        // sequence; 27 bases from the start to the coding sequence is in frame.
        let utr = b"CCCATGCCCCCCTAACCCCCCCCCTGACCC";
        assert_eq!(utr.len(), 30);
        let fixture = Fixture::new(1, &[]);
        let model = TranscriptUtr {
            strand: 1,
            segments: vec![(101, 130)],
            utr: utr.to_vec(),
            cds: CDS.as_bytes().to_vec(),
            scan: CodonScan::new(utr),
            utr_cds_first_stop: BTreeMap::new(),
            counts: existing_counts(&CodonScan::new(utr), utr.len()),
        };
        let variant = InputVariant::new("21".into(), 115, 115, b"A".to_vec(), b"C".to_vec());
        let out = fixture.plugin.annotate(&model, "21", &variant);
        assert_eq!(out[CONSEQUENCE_KEY], USTOP_LOST_FLAG);
        assert_eq!(
            out[ANNOTATION_KEY],
            "AltStop=True:AltStopDistanceToCDS=6:Evidence=NA:FrameWithCDS=inFrame:KozakContext=CCCATGC:KozakStrength=Weak"
        );
        // A deletion of the first 28 bases removes the start codon and both
        // stops: 2 bases remain and the start sits one base past the mutated
        // UTR's end, a distance of -1 whose remainder by 3 is 2 in Perl, out of
        // frame; no alternative stop remains.
        let variant = InputVariant::new("21".into(), 101, 128, utr[..28].to_vec(), b"-".to_vec());
        let out = fixture.plugin.annotate(&model, "21", &variant);
        assert_eq!(out[CONSEQUENCE_KEY], USTOP_LOST_FLAG);
        assert_eq!(
            out[ANNOTATION_KEY],
            "AltStop=False:AltStopDistanceToCDS=NA:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=CCCATGC:KozakStrength=Weak"
        );
    }

    /// An allele that begins one base before the `TGA` at 20 and makes the UTR
    /// three bases longer is read at index 20 of the mutated UTR, where `:639`
    /// reads it, on either strand. `CTG` at 19 replaced by `GGCCTA` leaves
    /// `GCC` there: the stop codon is lost though `TAA` stands at 23, the
    /// alternative stop, 30 bases from the coding sequence. `CT` at 19
    /// replaced by `GTGAC` writes `TGA` there: nothing is lost though index
    /// 23, three bases on, reads no stop codon.
    #[test]
    fn ustop_lost_reads_the_stop_codons_index_after_a_length_change() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let lost = utr_row(&fixture, 19, "CTG", "GGCCTA");
            assert_eq!(lost[CONSEQUENCE_KEY], USTOP_LOST_FLAG, "strand {strand}");
            assert_eq!(
                lost[ANNOTATION_KEY],
                "AltStop=True:AltStopDistanceToCDS=30:Evidence=NA:FrameWithCDS=inFrame:KozakContext=GCCATGG:KozakStrength=Strong",
                "strand {strand}"
            );
            let kept = utr_row(&fixture, 19, "CT", "GTGAC");
            assert_eq!(
                (
                    kept[CONSEQUENCE_KEY].as_str(),
                    kept[ANNOTATION_KEY].as_str()
                ),
                (NONE, NONE),
                "strand {strand}"
            );
        }
    }

    /// The `G` of `ATG` at 38 deleted: an in-frame overlapping ORF (12 bases to
    /// the coding sequence, 99 to the coding sequence's stop codon, since the
    /// frame runs into it) with a Moderate context `ACC ATG T`. The `T` of `ATG`
    /// at 5 substituted: the uORF (45 to the coding sequence, 15 to its stop).
    #[test]
    fn uaug_lost_reports_type_and_distances() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let oorf = utr_row(&fixture, 40, "G", "-");
            assert_eq!(oorf[CONSEQUENCE_KEY], UAUG_LOST_FLAG, "strand {strand}");
            assert_eq!(
                oorf[ANNOTATION_KEY],
                "CapDistanceToStart=38:DistanceToCDS=12:DistanceToStop=99:Evidence=NA:KozakContext=ACCATGT:KozakStrength=Moderate:type=InFrame_oORF"
            );
            let uorf = utr_row(&fixture, 6, "T", "C");
            assert_eq!(uorf[CONSEQUENCE_KEY], UAUG_LOST_FLAG, "strand {strand}");
            assert_eq!(
                uorf[ANNOTATION_KEY],
                "CapDistanceToStart=5:DistanceToCDS=45:DistanceToStop=15:Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:type=uORF"
            );
        }
    }

    /// An insertion directly before the first base of a start codon leaves the
    /// codon whole, moved on by the allele's length, so it is not lost. Three
    /// bases before the `ATG` at 5 fire no analysis. One base before the `ATG`
    /// at 38 lies inside the overlapping ORF at 30 and shifts its frame: 21
    /// bases to the coding sequence, in frame, and on to the coding stop.
    /// `:769` passes an insertion at the codon's own first index and
    /// `:776-778` read the inserted bases at the codon's old index.
    #[test]
    fn insertion_before_a_start_codon_leaves_it_whole() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let whole = utr_row(&fixture, 5, "-", "CCC");
            assert_eq!(whole[CONSEQUENCE_KEY], NONE, "strand {strand}");
            assert_eq!(whole[ANNOTATION_KEY], NONE, "strand {strand}");
            let shifted = utr_row(&fixture, 38, "-", "C");
            assert_eq!(
                shifted[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                "strand {strand}"
            );
            assert_eq!(
                shifted[ANNOTATION_KEY],
                "Evidence=NA:KozakContext=CCCATGC:KozakStrength=Weak:alt_type=InFrame_oORF:alt_type_length=111:ref_StartDistanceToCDS=20:ref_type=OutOfFrame_oORF:ref_type_length=NA",
                "strand {strand}"
            );
        }
    }

    /// Three bases after the first or the second base of the `ATG` at 38
    /// that neither end with the codon's first bases nor open with its last
    /// part the codon, which is lost; so is the codon whose first base a
    /// longer allele replaces. An insertion directly after the codon's third
    /// base touches no start codon.
    #[test]
    fn insertion_inside_a_start_codon_loses_it() {
        const LOST: &str = "CapDistanceToStart=38:DistanceToCDS=12:DistanceToStop=99:Evidence=NA:KozakContext=ACCATGT:KozakStrength=Moderate:type=InFrame_oORF";
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for index in [39, 40] {
                let parted = utr_row(&fixture, index, "-", "CCC");
                assert_eq!(
                    parted[CONSEQUENCE_KEY], UAUG_LOST_FLAG,
                    "strand {strand}, index {index}"
                );
                assert_eq!(
                    parted[ANNOTATION_KEY], LOST,
                    "strand {strand}, index {index}"
                );
            }
            let replaced = utr_row(&fixture, 38, "A", "CC");
            assert_eq!(replaced[CONSEQUENCE_KEY], UAUG_LOST_FLAG, "strand {strand}");
            assert_eq!(replaced[ANNOTATION_KEY], LOST, "strand {strand}");
            let after = utr_row(&fixture, 41, "-", "CCC");
            assert_eq!(after[CONSEQUENCE_KEY], NONE, "strand {strand}");
            assert_eq!(after[ANNOTATION_KEY], NONE, "strand {strand}");
        }
    }

    /// An allele inserted directly before the `ATG` at 38 that opens with a
    /// start codon, or holds one from its third base, is a gained start at 38
    /// or at 40.
    #[test]
    fn inserted_start_codon_before_a_start_codon_is_a_gain() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (allele, gained) in [("ATGCCC", "38"), ("CCATGC", "40")] {
                let out = utr_row(&fixture, 38, "-", allele);
                assert_eq!(
                    out[CONSEQUENCE_KEY], UAUG_GAINED_FLAG,
                    "strand {strand}, {allele}"
                );
                assert!(
                    out[ANNOTATION_KEY].starts_with(&format!("CapDistanceToStart={gained}:")),
                    "strand {strand}, {allele}: {}",
                    out[ANNOTATION_KEY]
                );
            }
        }
    }

    /// An insertion written after the first or the second base of the `ATG` at
    /// 30 whose inserted bases end with that base, or with those two, leaves the
    /// codon whole three bases on with every base after it unchanged: the
    /// mutated sequence of the same three bases written directly before the
    /// codon. No start codon is created or removed and no other analysis fires,
    /// so the row reports nothing, as the row of that other record does.
    #[test]
    fn insertion_inside_a_start_codon_that_only_moves_it_on_reports_nothing() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele, before) in [(31, "CCA", "ACC"), (32, "CAT", "ATC")] {
                let inside = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    inside[CONSEQUENCE_KEY], NONE,
                    "strand {strand}, index {index}"
                );
                assert_eq!(
                    inside[ANNOTATION_KEY], NONE,
                    "strand {strand}, index {index}"
                );
                assert_eq!(
                    inside,
                    utr_row(&fixture, 30, "-", before),
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// The `ATG` at 38 lies inside the overlapping ORF at 30. One base after
    /// the codon's first base that repeats it, or two after its second that
    /// repeat its first two, move the codon on and shift that ORF's frame: the
    /// row carries that ORF's frameshift alone, the row of the same bases
    /// written directly before the codon. Three bases shift no frame.
    #[test]
    fn insertion_that_moves_a_start_codon_on_inside_another_orf_reports_its_frameshift() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele) in [(39, "A"), (40, "AT")] {
                let inside = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    inside[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                    "strand {strand}, index {index}"
                );
                assert!(
                    inside[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=20:")
                        && !inside[ANNOTATION_KEY].contains(','),
                    "strand {strand}, index {index}: {}",
                    inside[ANNOTATION_KEY]
                );
                assert_eq!(
                    inside,
                    utr_row(&fixture, 38, "-", allele),
                    "strand {strand}, index {index}"
                );
            }
            let in_frame = utr_row(&fixture, 39, "-", "CCA");
            assert_eq!(in_frame[CONSEQUENCE_KEY], NONE, "strand {strand}");
            assert_eq!(in_frame[ANNOTATION_KEY], NONE, "strand {strand}");
        }
    }

    /// In [`UTR_OVERLAPS`] the uORF of the `ATG` at 5 ends with a stop codon
    /// that begins inside the `ATG` at 16. An insertion that only moves that
    /// `ATG` on stands before that stop codon, inside the uORF, wherever the
    /// record writes it. One base after the codon's first base cannot be
    /// written inside the stop codon and reports the uORF's frameshift. Two
    /// bases after its second base stand inside the stop codon, which reads
    /// no stop codon at its index: the row reports the lost stop, first in
    /// the fixed order, the uORF running on into the coding sequence. Three
    /// report nothing. Each is the row of the same mutated sequence written
    /// the other ways.
    #[test]
    fn insertion_that_moves_a_start_codon_on_is_read_directly_before_it() {
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            let shifted = utr_row(&fixture, 17, "-", "A");
            assert_eq!(
                shifted[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                "strand {strand}"
            );
            assert!(
                shifted[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=45:")
                    && !shifted[ANNOTATION_KEY].contains(','),
                "strand {strand}: {}",
                shifted[ANNOTATION_KEY]
            );
            assert_eq!(shifted, utr_row(&fixture, 16, "-", "A"), "strand {strand}");
            let lost = utr_row(&fixture, 18, "-", "AT");
            assert_eq!(lost[CONSEQUENCE_KEY], USTOP_LOST_FLAG, "strand {strand}");
            assert_eq!(
                lost[ANNOTATION_KEY],
                "AltStop=False:AltStopDistanceToCDS=NA:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong",
                "strand {strand}"
            );
            for (index, allele) in [(17, "TA"), (16, "AT")] {
                assert_eq!(
                    lost,
                    utr_row(&fixture, index, "-", allele),
                    "strand {strand}, index {index}"
                );
            }
            let in_frame = utr_row(&fixture, 18, "-", "CAT");
            assert_eq!(in_frame[CONSEQUENCE_KEY], NONE, "strand {strand}");
            assert_eq!(
                in_frame,
                utr_row(&fixture, 16, "-", "ATC"),
                "strand {strand}"
            );
        }
    }

    /// In [`UTR_OVERLAPS`] four or eight bases written after the second base
    /// of the `ATG` at 16 that move it on and complete no stop codon with
    /// that base stand inside the stop codon of the uORF at 5. That stop
    /// codon is lost and the row reports it, with the uORF's next stop codon
    /// in the mutated sequence: the reference's own, four bases on, 25 before
    /// the coding sequence; one the eight inserted bases hold, 38 before it;
    /// or none, the uORF running on into the coding sequence. The row is the
    /// same on the three ways of writing the sequence, and beside a
    /// read-through ORF, whose frameshift comes later in the fixed order.
    #[test]
    fn insertion_written_inside_the_stop_codon_at_a_moved_start_codon_loses_that_stop() {
        let read_through = UTR_OVERLAPS.replacen("CCG", "ATG", 1);
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            let beside = Fixture::with_utr(strand, &read_through, &[]);
            for (second, first, before, alt_stop) in [
                (
                    "CCAT",
                    "TCCA",
                    "ATCC",
                    "AltStop=True:AltStopDistanceToCDS=25",
                ),
                (
                    "CCTAGGAT",
                    "TCCTAGGA",
                    "ATCCTAGG",
                    "AltStop=True:AltStopDistanceToCDS=38",
                ),
                (
                    "CCCCCCAT",
                    "TCCCCCCA",
                    "ATCCCCCC",
                    "AltStop=False:AltStopDistanceToCDS=NA",
                ),
            ] {
                let lost = utr_row(&fixture, 18, "-", second);
                assert_eq!(
                    lost[CONSEQUENCE_KEY], USTOP_LOST_FLAG,
                    "strand {strand}, {second}: {}",
                    lost[ANNOTATION_KEY]
                );
                assert_eq!(
                    lost[ANNOTATION_KEY],
                    format!("{alt_stop}:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong"),
                    "strand {strand}, {second}"
                );
                for (index, allele) in [(17, first), (16, before)] {
                    assert_eq!(
                        lost,
                        utr_row(&fixture, index, "-", allele),
                        "strand {strand}, index {index}, {allele}"
                    );
                }
                for (index, allele) in [(18, second), (17, first), (16, before)] {
                    let row = utr_row(&beside, index, "-", allele);
                    assert_eq!(
                        row[CONSEQUENCE_KEY], USTOP_LOST_FLAG,
                        "strand {strand}, index {index}, {allele}"
                    );
                    assert_eq!(
                        row[ANNOTATION_KEY], lost[ANNOTATION_KEY],
                        "strand {strand}, index {index}, {allele}"
                    );
                }
            }
        }
    }

    /// In [`UTR_OVERLAPS`] the stop codon of the uORF at 5 begins at the second
    /// base of the `ATG` at 16. Four bases after that base whose first two
    /// complete a stop codon with it move the `ATG` on and leave the uORF
    /// ending on a stop codon at its own index, every codon before it
    /// unchanged: the uORF is not shifted and the row leaves it out, whichever
    /// of the three ways the record writes the sequence. With a start codon
    /// put at 0, whose overlapping ORF reads through the `ATG` at 16, the row
    /// lists that ORF alone; without one it reports nothing.
    #[test]
    fn insertion_that_moves_a_start_codon_on_leaves_out_the_orf_that_still_ends_inside_it() {
        let read_through = UTR_OVERLAPS.replacen("CCG", "ATG", 1);
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, &read_through, &[]);
            let alone = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (second, first, before) in [("AAAT", "TAAA", "ATAA"), ("AGAT", "TAGA", "ATAG")] {
                let inside = utr_row(&fixture, 18, "-", second);
                assert_eq!(
                    inside[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                    "strand {strand}, {second}"
                );
                assert!(
                    inside[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=50:")
                        && !inside[ANNOTATION_KEY].contains(','),
                    "strand {strand}, {second}: {}",
                    inside[ANNOTATION_KEY]
                );
                for (index, allele) in [(17, first), (16, before)] {
                    assert_eq!(
                        inside,
                        utr_row(&fixture, index, "-", allele),
                        "strand {strand}, index {index}"
                    );
                }
                for (index, allele) in [(18, second), (17, first), (16, before)] {
                    let kept = utr_row(&alone, index, "-", allele);
                    assert_eq!(
                        kept[CONSEQUENCE_KEY], NONE,
                        "strand {strand}, index {index}: {}",
                        kept[ANNOTATION_KEY]
                    );
                    assert_eq!(kept[ANNOTATION_KEY], NONE, "strand {strand}, index {index}");
                }
            }
        }
    }

    /// In [`UTR_OVERLAPS`] an insertion that moves the `ATG` at 16 on and cannot
    /// be written after that codon's second base stands before the first base
    /// of the stop codon that ends the uORF at 5. At a length that is no
    /// multiple of three the row reports that uORF's frameshift, written after
    /// the codon's first base or directly before the codon.
    #[test]
    fn insertion_before_the_stop_codon_inside_a_moved_start_codon_shifts_its_orf() {
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (index, allele) in [
                (17, "A"),
                (16, "A"),
                (17, "CCCA"),
                (16, "ACCC"),
                (16, "C"),
                (16, "CC"),
            ] {
                let out = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    out[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                    "strand {strand}, index {index}, {allele}"
                );
                assert!(
                    out[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=45:")
                        && !out[ANNOTATION_KEY].contains(','),
                    "strand {strand}, index {index}, {allele}: {}",
                    out[ANNOTATION_KEY]
                );
            }
        }
    }

    /// Six bases that move the `ATG` at 16 of [`UTR_OVERLAPS`] on shift no
    /// frame: the stop codon of the uORF at 5 stands whole six bases on, in
    /// frame, and the row reports nothing however the record writes them.
    #[test]
    fn in_frame_insertion_that_moves_a_start_codon_on_inside_a_stop_codon_reports_nothing() {
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (index, allele) in [(18, "CCCCAT"), (17, "CCCCCA"), (16, "CCCCCC")] {
                let out = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    out[CONSEQUENCE_KEY], NONE,
                    "strand {strand}, index {index}: {}",
                    out[ANNOTATION_KEY]
                );
                assert_eq!(out[ANNOTATION_KEY], NONE, "strand {strand}, index {index}");
            }
        }
    }

    /// At the `ATG` at 16 of [`UTR_OVERLAPS`], inside which the stop codon of
    /// the uORF at 5 begins, an allele that does not only move the codon on
    /// keeps its own term: four bases after its first base that part it lose
    /// the start codon, four directly before it that open with a start codon
    /// gain one, one base after its second that repeats its third loses the
    /// stop codon, and so do a deleted or substituted second base; a
    /// substituted first base loses the start codon.
    #[test]
    fn allele_at_a_start_codon_inside_a_stop_codon_that_does_not_only_move_it_keeps_its_term() {
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (index, reference, alternate, flag) in [
                (17, "-", "TAAC", UAUG_LOST_FLAG),
                (16, "-", "ATGA", UAUG_GAINED_FLAG),
                (18, "-", "G", USTOP_LOST_FLAG),
                (17, "T", "-", USTOP_LOST_FLAG),
                (17, "T", "C", USTOP_LOST_FLAG),
                (16, "A", "C", UAUG_LOST_FLAG),
            ] {
                let out = utr_row(&fixture, index, reference, alternate);
                assert_eq!(
                    out[CONSEQUENCE_KEY], flag,
                    "strand {strand}, index {index}, {reference}/{alternate}: {}",
                    out[ANNOTATION_KEY]
                );
            }
        }
    }

    /// Two bases that repeat the two before the `ATG` at 16 give one mutated
    /// sequence at five indices: 14 and 15 before the codon, the codon's
    /// first index, and 17 and 18 inside it, where 18 lies inside the stop
    /// codon of the uORF at 5. That uORF reads `TAT` where its stop codon
    /// stood and meets no other in the mutated UTR. Each of the five records
    /// carries its lost stop, where `:628-633` read the stop codon for the
    /// record at 18 alone and `:921` lists a frameshift for those at 14, 15
    /// and 16.
    #[test]
    fn repeated_bases_before_a_start_codon_lose_the_stop_codon_inside_it_at_every_index() {
        let utr = utr_with(UTR_OVERLAPS, 14..16, "AT");
        let lost = whole_row(
            USTOP_LOST_FLAG,
            "AltStop=False:AltStopDistanceToCDS=NA:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong",
            ("2", "1", "0"),
        );
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, &utr, &[]);
            for (index, allele) in [(14, "AT"), (15, "TA"), (16, "AT"), (17, "TA"), (18, "AT")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    lost,
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// Four bases that repeat the base before the `ATG` at 16 and leave a stop
    /// codon at 17 give one mutated sequence at four indices: 15 before the
    /// codon, its first index, and 17 and 18 inside it. The uORF at 5 reads
    /// every base before 17 unchanged and ends on the stop codon there: it
    /// keeps its length and no record reports it, where `:921` lists it for
    /// the records at 15 and 16. With a start codon put at 0, whose ORF reads
    /// through the inserted bases, each record lists that ORF alone.
    #[test]
    fn repeated_bases_before_a_start_codon_leave_out_the_orf_that_keeps_its_length_at_every_index()
    {
        let alone = utr_with(UTR_OVERLAPS, 15..16, "A");
        let read_through = alone.replacen("CCG", "ATG", 1);
        let nothing = whole_row(NONE, NONE, ("2", "1", "0"));
        let shifted = whole_row(
            UFRAMESHIFT_FLAG,
            "Evidence=NA:KozakContext=-:KozakStrength=:alt_type=uORF:alt_type_length=24:ref_StartDistanceToCDS=50:ref_type=OutOfFrame_oORF:ref_type_length=NA",
            ("2", "2", "0"),
        );
        for strand in [1, -1] {
            for (utr, row) in [(&alone, &nothing), (&read_through, &shifted)] {
                let fixture = Fixture::with_utr(strand, utr, &[]);
                for (index, allele) in [(15, "AATA"), (16, "ATAA"), (17, "TAAA"), (18, "AAAT")] {
                    assert_eq!(
                        &utr_row(&fixture, index, "-", allele),
                        row,
                        "strand {strand}, index {index}"
                    );
                }
            }
        }
    }

    /// In [`UTR_OVERLAPS`] the stop codon at 25, which ends the uORF at 16,
    /// lends its last base to the `ATG` at 27. Two bases that repeat its
    /// first two give one mutated sequence at 25, directly before the stop
    /// codon, at 26 inside it, and at 27, the start codon's first index. The
    /// uORF reads `TGT` where its stop codon stood and runs into the coding
    /// sequence: each record carries its lost stop, where `:631` skips the
    /// record at 25.
    #[test]
    fn repeated_bases_of_a_stop_codon_that_ends_at_a_start_codon_lose_it_at_every_index() {
        let lost = whole_row(
            USTOP_LOST_FLAG,
            "AltStop=False:AltStopDistanceToCDS=NA:Evidence=NA:FrameWithCDS=inFrame:KozakContext=CCCATGA:KozakStrength=Weak",
            ("2", "1", "0"),
        );
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (index, allele) in [(25, "TG"), (26, "GT"), (27, "TG")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    lost,
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// Bases written after the first base of the `ATG` at 16 of
    /// [`UTR_OVERLAPS`] that open with its last two leave that codon at its
    /// own index and can be written at 17, the first base of the stop codon
    /// of the uORF at 5, at 18 inside that stop codon, and at 19 directly
    /// after the start codon. The stop codon stands whole after the inserted
    /// bases. Three bases leave it in frame: the uORF gains a codon and ends
    /// on it, and no record reports anything, where `:639-642` report the
    /// stop lost for the records at 18 and 19. Four bases put it out of
    /// frame: each record carries the lost stop, with the next stop codon the
    /// uORF meets, where `:631` skips the record at 17.
    #[test]
    fn insertion_writable_before_and_inside_a_stop_codon_loses_it_only_out_of_frame() {
        let nothing = whole_row(NONE, NONE, ("2", "1", "0"));
        let lost = whole_row(
            USTOP_LOST_FLAG,
            "AltStop=True:AltStopDistanceToCDS=25:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong",
            ("2", "1", "0"),
        );
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (index, allele) in [(17, "TGC"), (18, "GCT"), (19, "CTG")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    nothing,
                    "strand {strand}, index {index}"
                );
            }
            for (index, allele) in [(17, "TGCC"), (18, "GCCT"), (19, "CCTG")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    lost,
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// Bases that part the `ATG` at 16 of [`UTR_OVERLAPS`] and repeat its
    /// second base can be written after its first base, at the first base of
    /// the stop codon of the uORF at 5, and after its second, inside that
    /// stop codon. The start codon is lost either way. One base puts the stop
    /// codon, whole one base on, out of frame: both records carry the lost
    /// stop, first in the fixed order, where `:631` skips the record at 17
    /// and the lost start is its only term. Three bases leave the stop codon
    /// in frame: both records carry the lost start, where `:639-642` report
    /// the stop lost for the record at 18.
    #[test]
    fn insertion_that_parts_a_start_codon_inside_a_stop_codon_has_one_term_at_either_index() {
        let lost_stop = whole_row(
            USTOP_LOST_FLAG,
            "AltStop=True:AltStopDistanceToCDS=25:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong",
            ("2", "1", "0"),
        );
        let lost_start = whole_row(
            UAUG_LOST_FLAG,
            "CapDistanceToStart=16:DistanceToCDS=34:DistanceToStop=9:Evidence=NA:KozakContext=CCCATGA:KozakStrength=Weak:type=uORF",
            ("2", "1", "0"),
        );
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for (index, allele) in [(17, "T"), (18, "T")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    lost_stop,
                    "strand {strand}, index {index}"
                );
            }
            for (index, allele) in [(17, "TCC"), (18, "CCT")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    lost_start,
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// A start codon the inserted bases make is gained whether the reference's
    /// codon stands whole behind it or is parted. After the first base of the
    /// `ATG` at 38, `CATGA` makes a start codon at 40 and leaves the
    /// reference's whole at 43; after its second base, `CATGAT` makes one at 41
    /// and leaves the reference's whole at 44; `CATGC` after its first base
    /// makes one at 40 and parts the reference's. The row names the made codon.
    #[test]
    fn start_codon_made_inside_a_start_codon_is_a_gain() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele, gained) in [
                (39, "CATGA", "40"),
                (40, "CATGAT", "41"),
                (39, "CATGC", "40"),
            ] {
                let out = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    out[CONSEQUENCE_KEY], UAUG_GAINED_FLAG,
                    "strand {strand}, {allele}"
                );
                assert!(
                    out[ANNOTATION_KEY].starts_with(&format!("CapDistanceToStart={gained}:")),
                    "strand {strand}, {allele}: {}",
                    out[ANNOTATION_KEY]
                );
            }
        }
    }

    /// A longer allele that replaces the second base of the `ATG` at 38 parts
    /// the codon, which is lost: only an insertion can move a codon on whole.
    #[test]
    fn replacement_that_begins_inside_a_start_codon_loses_it() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let replaced = utr_row(&fixture, 39, "T", "CA");
            assert_eq!(replaced[CONSEQUENCE_KEY], UAUG_LOST_FLAG, "strand {strand}");
            assert!(
                replaced[ANNOTATION_KEY].starts_with("CapDistanceToStart=38:"),
                "strand {strand}: {}",
                replaced[ANNOTATION_KEY]
            );
        }
    }

    /// In [`UTR_OVERLAPS`] an insertion directly before the `ATG` at 27 lies
    /// after the second base of the stop codon at 25, which ends the uORF of
    /// the `ATG` at 16: the lost stop comes first in the fixed order.
    #[test]
    fn insertion_before_a_start_codon_inside_a_stop_codon_reports_the_lost_stop() {
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            for allele in ["C", "CCC"] {
                let out = utr_row(&fixture, 27, "-", allele);
                assert_eq!(
                    out[CONSEQUENCE_KEY], USTOP_LOST_FLAG,
                    "strand {strand}, {allele}"
                );
            }
        }
    }

    /// Three bases written after the first or the second base of the `ATG` at
    /// 30 that open with the codon's remaining bases leave the codon at 30: the
    /// mutated sequence of three bases written directly after the codon. No
    /// start codon is created, three bases shift no frame, and the row reports
    /// nothing, as the row of that other record does.
    #[test]
    fn insertion_that_leaves_a_start_codon_at_its_own_index_gains_nothing() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele, after) in [(31, "TGC", "CTG"), (32, "GCC", "CCG")] {
                let inside = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    inside[CONSEQUENCE_KEY], NONE,
                    "strand {strand}, index {index}"
                );
                assert_eq!(
                    inside[ANNOTATION_KEY], NONE,
                    "strand {strand}, index {index}"
                );
                assert_eq!(
                    inside,
                    utr_row(&fixture, 33, "-", after),
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// One base after the second base of a start codon that repeats its third,
    /// or two after its first that repeat its last two, stand directly after
    /// the codon and shift its own ORF. At the `ATG` at 30 the row reports that
    /// ORF's frameshift alone; at the `ATG` at 38, which lies inside the ORF of
    /// the codon at 30, both ORFs' frameshifts. Each is the row of the same
    /// bases written directly after the codon.
    #[test]
    fn insertion_that_leaves_a_start_codon_at_its_own_index_shifts_its_orf() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele) in [(31, "TG"), (32, "G")] {
                let inside = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    inside[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                    "strand {strand}, index {index}"
                );
                assert!(
                    inside[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=20:")
                        && !inside[ANNOTATION_KEY].contains(','),
                    "strand {strand}, index {index}: {}",
                    inside[ANNOTATION_KEY]
                );
                assert_eq!(
                    inside,
                    utr_row(&fixture, 33, "-", allele),
                    "strand {strand}, index {index}"
                );
            }
            let inside = utr_row(&fixture, 40, "-", "G");
            assert_eq!(inside[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG, "strand {strand}");
            assert!(
                inside[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=20:")
                    && inside[ANNOTATION_KEY].contains(":ref_StartDistanceToCDS=12:"),
                "strand {strand}: {}",
                inside[ANNOTATION_KEY]
            );
            assert_eq!(inside, utr_row(&fixture, 41, "-", "G"), "strand {strand}");
        }
    }

    /// Inserted bases that leave the `ATG` at 38 in place and hold a start
    /// codon of their own make a start codon at 41: the row names that codon,
    /// not the reference's at 38, as the row of the same bases written directly
    /// after the codon does.
    #[test]
    fn start_codon_made_after_one_left_at_its_own_index_is_the_gain() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele, after) in [(39, "TGATGC", "ATGCTG"), (40, "GATGC", "ATGCG")] {
                let inside = utr_row(&fixture, index, "-", allele);
                assert_eq!(
                    inside[CONSEQUENCE_KEY], UAUG_GAINED_FLAG,
                    "strand {strand}, index {index}"
                );
                assert!(
                    inside[ANNOTATION_KEY].starts_with("CapDistanceToStart=41:"),
                    "strand {strand}, index {index}: {}",
                    inside[ANNOTATION_KEY]
                );
                assert_eq!(
                    inside,
                    utr_row(&fixture, 41, "-", after),
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// In [`UTR_OVERLAPS`] one base after the second base of the `ATG` at 16
    /// that repeats its third leaves that codon in place and stands inside the
    /// stop codon at 17, which the mutated sequence reads no more: the row
    /// reports the lost stop of the uORF of the `ATG` at 5, first in the fixed
    /// order, as the row of the same base written directly after the codon
    /// does.
    #[test]
    fn insertion_that_leaves_a_start_codon_at_its_own_index_inside_a_stop_codon_loses_the_stop() {
        for strand in [1, -1] {
            let fixture = Fixture::with_utr(strand, UTR_OVERLAPS, &[]);
            let inside = utr_row(&fixture, 18, "-", "G");
            assert_eq!(inside[CONSEQUENCE_KEY], USTOP_LOST_FLAG, "strand {strand}");
            assert_eq!(inside, utr_row(&fixture, 19, "-", "G"), "strand {strand}");
        }
    }

    /// A duplication of the `ATG` at 38 gives one mutated sequence at four
    /// indices: directly before the codon, after its first and its second
    /// base, and directly after it. The mutated UTR holds two equal start
    /// codons, at 38 and 41, where the reference holds one, and nothing in it
    /// says which of the two is the added one. Each record names the copy
    /// nearer the cap, the first start codon of the changed bases, as
    /// `:346-353` name the first start codon of their window: the gained
    /// start at 38 with that copy's own context and distances.
    #[test]
    fn a_duplicated_start_codon_names_one_copy_for_each_way_of_writing_it() {
        let gained = whole_row(
            UAUG_GAINED_FLAG,
            "CapDistanceToStart=38:DistanceToCDS=15:DistanceToStop=102:KozakContext=ACCATGA:KozakStrength=Moderate:type=inFrame_oORF",
            ("1", "1", "1"),
        );
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (index, allele) in [(38, "ATG"), (39, "TGA"), (40, "GAT"), (41, "ATG")] {
                assert_eq!(
                    utr_row(&fixture, index, "-", allele),
                    gained,
                    "strand {strand}, index {index}"
                );
            }
        }
    }

    /// One base inserted at 44 sits inside both overlapping ORFs (30 and 38) and
    /// shifts their frames: the out-of-frame one becomes in frame and runs to
    /// the coding sequence's stop, the in-frame one becomes out of frame with no
    /// stop. A three-base insertion keeps the frame and a substitution is not
    /// an indel.
    #[test]
    fn uframeshift_is_indel_only() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let ins = utr_row(&fixture, 44, "-", "C");
            assert_eq!(ins[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG, "strand {strand}");
            assert_eq!(
                ins[ANNOTATION_KEY],
                "Evidence=NA:KozakContext=CCCATGC:KozakStrength=Weak:alt_type=InFrame_oORF:alt_type_length=111:ref_StartDistanceToCDS=20:ref_type=OutOfFrame_oORF:ref_type_length=NA,\
                 Evidence=NA:KozakContext=ACCATGT:KozakStrength=Moderate:alt_type=OutOfFrame_oORF:alt_type_length=NA:ref_StartDistanceToCDS=12:ref_type=InFrame_oORF:ref_type_length=102"
            );
            let in_frame = utr_row(&fixture, 44, "-", "CCC");
            assert_eq!(in_frame[CONSEQUENCE_KEY], NONE, "strand {strand}");
            let snv = utr_row(&fixture, 44, "C", "G");
            assert_eq!(snv[CONSEQUENCE_KEY], NONE, "strand {strand}");
        }
    }

    #[test]
    fn text_rendering_joins_effects_with_colon_and_comma() {
        let effects: Effects = vec![
            vec![
                ("type", "uORF".to_string()),
                ("DistanceToCDS", "5".to_string()),
            ],
            vec![
                ("KozakStrength", String::new()),
                ("KozakContext", NONE.to_string()),
            ],
        ];
        assert_eq!(
            render_effects(&effects),
            "DistanceToCDS=5:type=uORF,KozakContext=-:KozakStrength="
        );
        assert_eq!(render_effects(&Effects::new()), "");
    }

    /// With a file the uORF at 5 (genomic 106 on the forward fixture) is `True`
    /// when the file names `chr21:106`, `False` when it names another position,
    /// and `NA` without a file.
    #[test]
    fn evidence_true_false_and_na_without_a_file() {
        let header = "chr\tpos\tgene\tstrand\ttype\tstop_pos\n";
        let dir = tempfile::tempdir().unwrap();
        let listed = dir.path().join("listed.txt");
        std::fs::write(
            &listed,
            format!("{header}chr21\t106\tGENE\tforward\tfive_prime_utr\t121\n"),
        )
        .unwrap();
        let other = dir.path().join("other.txt");
        std::fs::write(
            &other,
            format!("{header}chr21\t107\tGENE\tforward\tfive_prime_utr\t121\n"),
        )
        .unwrap();
        let evidence = |fixture: &Fixture| {
            let out = utr_row(fixture, 16, "C", "A");
            assert_eq!(out[CONSEQUENCE_KEY], USTOP_GAINED_FLAG);
            out[ANNOTATION_KEY]
                .split(':')
                .find_map(|pair| pair.strip_prefix("Evidence="))
                .unwrap()
                .to_string()
        };
        assert_eq!(
            evidence(&Fixture::new(1, &[&format!("file={}", listed.display())])),
            "True"
        );
        assert_eq!(
            evidence(&Fixture::new(1, &[&format!("file={}", other.display())])),
            "False"
        );
        assert_eq!(evidence(&Fixture::new(1, &[])), "NA");
        // The reverse fixture's uORF start sits at genomic 295.
        let reverse = dir.path().join("reverse.txt");
        std::fs::write(
            &reverse,
            format!("{header}chr21\t295\tGENE\treverse\tfive_prime_utr\t280\n"),
        )
        .unwrap();
        assert_eq!(
            evidence(&Fixture::new(-1, &[&format!("file={}", reverse.display())])),
            "True"
        );
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) is a byte of its cell, as
    /// `<$fh>` and `split /\t/` (UTRAnnotator.pm:96-98) read the bytes of a
    /// line: in the gene cell, which the key `chr:pos` never reads, the uORF at
    /// genomic 106 keeps its evidence; in the chromosome cell it names another
    /// chromosome, so the evidence is `False`; in the position cell the
    /// position is no number, refused as `looks_like_number` refuses it
    /// (`:102`); and a `#` line is refused for its chromosome, as `:101` dies
    /// on any line not opening `ch`.
    #[test]
    fn a_byte_that_is_not_utf8_is_a_byte_of_its_cell() {
        let header = b"chr\tpos\tgene\tstrand\ttype\tstop_pos\n";
        let dir = tempfile::tempdir().unwrap();
        let file = |name: &str, line: &[u8]| {
            let path = dir.path().join(name);
            std::fs::write(&path, [header, line].concat()).unwrap();
            format!("file={}", path.display())
        };
        let evidence = |param: String| {
            let out = utr_row(&Fixture::new(1, &[&param]), 16, "C", "A");
            out[ANNOTATION_KEY]
                .split(':')
                .find_map(|pair| pair.strip_prefix("Evidence="))
                .unwrap()
                .to_string()
        };
        assert_eq!(
            evidence(file(
                "gene.txt",
                b"chr21\t106\tGEN\xE9\tforward\tfive_prime_utr\t121\n"
            )),
            "True"
        );
        assert_eq!(
            evidence(file(
                "chr.txt",
                b"chr2\xE91\t106\tGENE\tforward\tfive_prime_utr\t121\n"
            )),
            "False"
        );
        let refused = |param: String| {
            UtrAnnotatorPlugin::new()
                .init(&[param])
                .unwrap_err()
                .to_string()
        };
        let err = refused(file(
            "pos.txt",
            b"chr21\t10\xE96\tGENE\tforward\tfive_prime_utr\t121\n",
        ));
        assert!(err.contains("position '10\u{FFFD}6' on line 2"), "{err}");
        let err = refused(file(
            "comment.txt",
            b"# caf\xE9\t106\tGENE\tforward\tfive_prime_utr\t121\n",
        ));
        assert!(
            err.contains("chromosome '# caf\u{FFFD}' on line 2"),
            "{err}"
        );
    }

    /// An evidence line one byte over the cap is refused as the plugin opens,
    /// naming the file and the line.
    #[test]
    fn evidence_line_past_the_cap_is_refused() {
        let cap = crate::tabix::MAX_LINE_BYTES;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("uORF.txt");
        let mut text = b"chr\tpos\tgene\tstrand\ttype\tstop_pos\nchr21\t106\t".to_vec();
        text.extend(std::iter::repeat_n(b'x', cap + 1));
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        let err = UtrAnnotatorPlugin::new()
            .init(&[format!("file={}", path.display())])
            .unwrap_err()
            .to_string();
        assert_eq!(
            err,
            format!(
                "plugin initialization failed: {}: line 2 is longer than the {cap}-byte limit \
                 for one line",
                path.display()
            )
        );
    }

    #[test]
    fn missing_file_is_an_init_error() {
        let mut plugin = UtrAnnotatorPlugin::new();
        let err = plugin
            .init(&["file=/no/such/uORF_5UTR_GRCh38_PUBLIC.txt".to_string()])
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("UTRAnnotator: uORF file /no/such/uORF_5UTR_GRCh38_PUBLIC.txt not found"),
            "{err}"
        );
        let dir = tempfile::tempdir().unwrap();
        let bare = dir.path().join("bare.txt");
        std::fs::write(
            &bare,
            "chr\tpos\tgene\tstrand\ttype\tstop_pos\n21\t106\tG\tforward\tt\t121\n",
        )
        .unwrap();
        let err = UtrAnnotatorPlugin::new()
            .init(&[format!("file={}", bare.display())])
            .unwrap_err()
            .to_string();
        assert!(err.contains("chromosome '21' on line 2"), "{err}");
        let text = dir.path().join("text.txt");
        std::fs::write(&text, "chr21\tabc\tG\tforward\tt\t121\n").unwrap();
        let err = UtrAnnotatorPlugin::new()
            .init(&[format!("file={}", text.display())])
            .unwrap_err()
            .to_string();
        assert!(err.contains("position 'abc' on line 1"), "{err}");
        // A token without `=` leaves Ensembl VEP's parameter hash empty
        // (BaseVepPlugin.pm:123): the plugin loads without the file.
        let mut bare_token = UtrAnnotatorPlugin::new();
        bare_token
            .init(&["/a/positional/path".to_string()])
            .unwrap();
        assert!(bare_token.evidence.is_none());
        let err = UtrAnnotatorPlugin::new()
            .init(&["cover=1".to_string()])
            .unwrap_err()
            .to_string();
        assert!(err.contains("unknown parameter 'cover'"), "{err}");
    }

    /// A ten-base deletion of the second UTR exon covers 90 percent of it by
    /// `:187`'s measure: under `max_overlap=80` the exon is dropped, the variant's
    /// position maps nowhere and only the counts and `-` are written; under the
    /// default 100 nothing is ever dropped and the deletion is read (it removes
    /// the start codon at 38 and shifts the frame of the one at 30; the lost
    /// start codon wins by precedence). A deletion of both exons drops them
    /// both and writes no key.
    #[test]
    fn max_overlap_drops_utr_exons_the_variant_covers() {
        let exon2_utr = &UTR[UTR_SPLIT..];
        let dropped = Fixture::new(1, &["max_overlap=80"]);
        let out = utr_row(&dropped, UTR_SPLIT, exon2_utr, "-");
        assert_eq!(out[CONSEQUENCE_KEY], NONE);
        assert_eq!(out[ANNOTATION_KEY], NONE);
        assert_eq!(counts_of(&out), ("1", "1", "1"));

        let kept = Fixture::new(1, &[]);
        let out = utr_row(&kept, UTR_SPLIT, exon2_utr, "-");
        assert_eq!(out[CONSEQUENCE_KEY], UAUG_LOST_FLAG);

        let [(exon1_start, exon1_end), (exon2_start, _)] = layout(1).0;
        let whole: String = format!(
            "{}{}{}",
            &UTR[..UTR_SPLIT],
            "C".repeat((exon2_start - exon1_end - 1) as usize),
            exon2_utr
        );
        let mut variant = InputVariant::new(
            "21".into(),
            exon1_start,
            exon2_start + 9,
            whole.into_bytes(),
            b"-".to_vec(),
        );
        variant.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.transcript = Some(Arc::new({
                let mut transcript = Transcript::default();
                transcript.stable_id = Arc::from("ENST00000000001");
                transcript
            }));
            consequence.consequences =
                ConsequenceList::from_slice(&[Consequence::FivePrimeUtrVariant]);
            consequence.loftee_ctx = Some(Box::new(dropped.context()));
            consequence
        });
        dropped
            .plugin
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();
        assert!(variant.transcript_consequences[0].plugin_data.is_empty());
    }

    /// `next if $UTR_overlap_perc >= $self->{max_overlap}` (`:188`) numifies the
    /// value `:115` stores: `80%` is 80 and `1e2` is 100, so the ten-base
    /// deletion covering 90 percent of the second exon is dropped under the first
    /// and read under the second; `abc` and the empty value are 0, so the exon
    /// under a substitution (overlap 0, and `0 >= 0`) is dropped and the row
    /// carries the counts and `-` alone; `nan` compares false with every
    /// percentage and drops nothing.
    #[test]
    fn max_overlap_is_read_as_perl_numifies_it() {
        let exon2_utr = &UTR[UTR_SPLIT..];
        let percent = Fixture::new(1, &["max_overlap=80%"]);
        assert_eq!(percent.plugin.max_overlap, 80.0);
        let out = utr_row(&percent, UTR_SPLIT, exon2_utr, "-");
        assert_eq!(out[CONSEQUENCE_KEY], NONE);
        assert_eq!(out[ANNOTATION_KEY], NONE);
        assert_eq!(counts_of(&out), ("1", "1", "1"));

        let exponent = Fixture::new(1, &["max_overlap=1e2"]);
        assert_eq!(exponent.plugin.max_overlap, 100.0);
        let out = utr_row(&exponent, UTR_SPLIT, exon2_utr, "-");
        assert_eq!(out[CONSEQUENCE_KEY], UAUG_LOST_FLAG);

        for value in ["max_overlap=", "max_overlap=abc", "max_overlap=0"] {
            let zero = Fixture::new(1, &[value]);
            assert_eq!(zero.plugin.max_overlap, 0.0, "{value}");
            let out = utr_row(&zero, 25, "C", "T");
            assert_eq!(out[CONSEQUENCE_KEY], NONE, "{value}");
            assert_eq!(out[ANNOTATION_KEY], NONE, "{value}");
            assert_eq!(counts_of(&out), ("1", "1", "1"), "{value}");
        }
        let default = Fixture::new(1, &[]);
        let out = utr_row(&default, 25, "C", "T");
        assert_eq!(out[CONSEQUENCE_KEY], UAUG_GAINED_FLAG);

        let nan = Fixture::new(1, &["max_overlap=nan"]);
        assert!(nan.plugin.max_overlap.is_nan());
        let out = utr_row(&nan, UTR_SPLIT, exon2_utr, "-");
        assert_eq!(out[CONSEQUENCE_KEY], UAUG_LOST_FLAG);
    }

    /// `CATGC` inserted before index 24 carries a start codon three bases past
    /// the insertion point: a window sized by the alternate allele reaches it
    /// (`CapDistanceToStart=25`), one sized by the reference allele's zero
    /// length (`:339`) reads four bases and misses it. A three-base deletion at
    /// 26 brings no new start codon to the junction, while `:339`'s window,
    /// three bases longer, reads the existing `ATG` at 30 as gained.
    #[test]
    fn atg_window_is_sized_by_the_alternate_allele() {
        let fixture = Fixture::new(1, &[]);
        let ins = utr_row(&fixture, 24, "-", "CATGC");
        assert_eq!(ins[CONSEQUENCE_KEY], UAUG_GAINED_FLAG);
        assert!(
            ins[ANNOTATION_KEY].starts_with("CapDistanceToStart=25:"),
            "{}",
            ins[ANNOTATION_KEY]
        );
        let mutated = mutated_utr(UTR.as_bytes(), 24, 0, b"CATGC");
        assert_eq!(atg_window(&mutated, 24, 5).0, b"ACCATGCAC");
        assert_eq!(atg_window(&mutated, 24, 0).0, b"ACCA");

        let del = utr_row(&fixture, 26, "GCC", "-");
        assert_eq!(del[CONSEQUENCE_KEY], NONE);
        let mutated = mutated_utr(UTR.as_bytes(), 26, 3, b"");
        assert_eq!(atg_window(&mutated, 26, 0).0, b"ACCA");
        assert_eq!(atg_window(&mutated, 26, 3).0, b"ACCATGC");

        // The window never wraps below the sequence's start.
        assert_eq!(atg_window(b"ATGCCC", 1, 1), (b"ATGC".as_slice(), 0));
    }

    /// When two analyses fire the flag is the first of `uAUG_gained`,
    /// `uSTOP_gained`, `uSTOP_lost`, `uAUG_lost`, `uFrameshift`: the ten-base
    /// deletion of the second UTR exon is both a lost start codon (38) and a
    /// frameshift of the overlapping ORF at 30.
    #[test]
    fn two_analyses_choose_the_fixed_precedence() {
        let fixture = Fixture::new(1, &[]);
        let out = utr_row(&fixture, UTR_SPLIT, &UTR[UTR_SPLIT..], "-");
        assert_eq!(out[CONSEQUENCE_KEY], UAUG_LOST_FLAG);
        assert!(out[ANNOTATION_KEY].contains("CapDistanceToStart=38:"));
    }

    /// A direct caller that hands the plugin no reference gets no row
    /// annotated; the registry refuses the plugin before any batch
    /// (`missing_fasta_is_an_init_error`).
    #[test]
    fn no_fasta_annotates_nothing() {
        let fixture = Fixture::new(1, &[]);
        let mut plugin = UtrAnnotatorPlugin::new();
        plugin.init(&[]).unwrap();
        plugin.set_reference_fasta(None);
        let mut variant = InputVariant::new("21".into(), 126, 126, b"C".to_vec(), b"T".to_vec());
        variant.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.consequences =
                ConsequenceList::from_slice(&[Consequence::FivePrimeUtrVariant]);
            consequence.loftee_ctx = Some(Box::new(fixture.context()));
            consequence
        });
        plugin
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();
        assert!(variant.transcript_consequences[0].plugin_data.is_empty());
    }

    /// The registry refuses UTRAnnotator in a run without a reference and
    /// loads it in a run with one. Ensembl VEP's module prints a notice and
    /// runs on (UTRAnnotator.pm:88).
    #[test]
    fn missing_fasta_is_an_init_error() {
        let args = ["UTRAnnotator".to_string()];
        let mut registry = crate::BuiltinRegistry::new();
        let outcome = registry.load_from_args_with_context(&args, None, None);
        let [(name, error)] = &outcome.failed[..] else {
            panic!("one failure, got {:?}", outcome.failed);
        };
        assert_eq!(name, "UTRAnnotator");
        assert_eq!(
            error.to_string(),
            "plugin initialization failed: UTRAnnotator: a reference FASTA is needed; pass --fasta"
        );
        assert!(registry.is_empty());

        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let mut contents = b">21\n".to_vec();
        contents.extend_from_slice(&[b'C'; 400]);
        contents.push(b'\n');
        std::fs::write(&fa, contents).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t400\t4\t400\t401\n").unwrap();
        let reference = Arc::new(vep_fasta::IndexedFasta::from_path(&fa).unwrap());
        let mut registry = crate::BuiltinRegistry::new();
        let outcome = registry.load_from_args_with_context(&args, Some(reference), None);
        assert!(outcome.failed.is_empty(), "{:?}", outcome.failed);
        assert_eq!(registry.plugin_names(), ["UTRAnnotator"]);
    }

    /// A symbolic deletion whose stated span lies inside one UTR exon segment
    /// carries the five keys of the same deletion given with explicit alleles,
    /// on both strands and for each analysis a deletion can fire. Index 8 lies
    /// inside the uORF at 5 and shifts its frame. Indices 36 to 38 bring the
    /// `A` at 35 and the `TG` at 39 together. Index 5 opens the start codon at 5
    /// and index 20 the uORF's stop codon. Index 38 opens the start codon at 38
    /// and shifts the overlapping ORF at 30: the lost start codon is the first
    /// of the two in the fixed order. Index 41, in the second exon's segment,
    /// shifts both overlapping ORFs. Index 0 lies before every start codon.
    #[test]
    fn a_symbolic_deletion_inside_one_utr_exon_is_the_sequence_deletion_of_its_span() {
        let cases: [(usize, usize, &str, &str); 7] = [
            (8, 8, UFRAMESHIFT_FLAG, "Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:alt_type=OutOfFrame_oORF:alt_type_length=NA:ref_StartDistanceToCDS=45:ref_type=uORF:ref_type_length=18"),
            (36, 38, UAUG_GAINED_FLAG, "CapDistanceToStart=35:DistanceToCDS=12:DistanceToStop=99:KozakContext=GCCATGT:KozakStrength=Moderate:type=inFrame_oORF"),
            (5, 5, UAUG_LOST_FLAG, "CapDistanceToStart=5:DistanceToCDS=45:DistanceToStop=15:Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:type=uORF"),
            (20, 20, USTOP_LOST_FLAG, "AltStop=False:AltStopDistanceToCDS=NA:Evidence=NA:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong"),
            (38, 38, UAUG_LOST_FLAG, "CapDistanceToStart=38:DistanceToCDS=12:DistanceToStop=99:Evidence=NA:KozakContext=ACCATGT:KozakStrength=Moderate:type=InFrame_oORF"),
            (41, 41, UFRAMESHIFT_FLAG, "Evidence=NA:KozakContext=CCCATGC:KozakStrength=Weak:alt_type=OutOfFrame_oORF:alt_type_length=NA:ref_StartDistanceToCDS=20:ref_type=OutOfFrame_oORF:ref_type_length=NA,Evidence=NA:KozakContext=ACCATGT:KozakStrength=Moderate:alt_type=OutOfFrame_oORF:alt_type_length=NA:ref_StartDistanceToCDS=12:ref_type=InFrame_oORF:ref_type_length=102"),
            (0, 0, NONE, NONE),
        ];
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (first, last, flag, annotation) in cases {
                let symbolic = symbolic_deletion(&fixture, first, last);
                assert_eq!(
                    symbolic[CONSEQUENCE_KEY], flag,
                    "strand {strand}, {first}..={last}"
                );
                assert_eq!(
                    symbolic[ANNOTATION_KEY], annotation,
                    "strand {strand}, {first}..={last}"
                );
                assert_eq!(counts_of(&symbolic), ("1", "1", "1"));
                let explicit = utr_row(&fixture, first, &UTR[first..=last], "-");
                assert_eq!(symbolic, explicit, "strand {strand}, {first}..={last}");
            }
        }
    }

    /// A mobile-element deletion states its span as a deletion does, and a
    /// confidence interval of no width leaves the span as stated.
    #[test]
    fn a_mobile_element_deletion_and_an_interval_of_no_width_are_annotated() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let explicit = utr_row(&fixture, 8, &UTR[8..=9], "-");
            assert_eq!(explicit[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG);
            let mobile = symbolic_row(
                &fixture,
                VariantClass::MobileElementDeletion,
                fixture.span(8, 9),
                |_| {},
            );
            assert_eq!(mobile, explicit, "strand {strand}");
            let exact = symbolic_row(
                &fixture,
                VariantClass::StructuralDeletion,
                fixture.span(8, 9),
                |record| {
                    record.ci_pos = Some((0, 0));
                    record.ci_end = Some((0, 0));
                },
            );
            assert_eq!(exact, explicit, "strand {strand}");
        }
    }

    /// Under `max_overlap=50` a symbolic deletion of two bases of the first UTR
    /// exon keeps that exon and is annotated; one of the whole exon drops it, so
    /// no end of the span has an index; one that covers both UTR exons gets no
    /// key at all, as the sequence deletion of each span does.
    #[test]
    fn max_overlap_drops_utr_exons_a_symbolic_deletion_covers() {
        let fixture = Fixture::new(1, &["max_overlap=50"]);
        let kept = symbolic_deletion(&fixture, 8, 9);
        assert_eq!(kept[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG);
        assert_eq!(kept, utr_row(&fixture, 8, &UTR[8..=9], "-"));
        let dropped = symbolic_deletion(&fixture, 0, UTR_SPLIT - 1);
        assert_eq!(dropped[CONSEQUENCE_KEY], NONE);
        assert_eq!(dropped[ANNOTATION_KEY], NONE);
        assert_eq!(counts_of(&dropped), ("1", "1", "1"));
        assert_eq!(dropped, utr_row(&fixture, 0, &UTR[..UTR_SPLIT], "-"));
        let both = symbolic_row(
            &fixture,
            VariantClass::StructuralDeletion,
            (101, 170),
            |_| {},
        );
        assert!(both.is_empty(), "{both:?}");
    }

    /// A plain symbolic insertion gives a length and no bases, so of the five
    /// analyses the frameshift alone follows: the length is no multiple of
    /// three and the point lies inside a uORF or an overlapping ORF
    /// (UTRAnnotator.pm:889, :921-923). Each effect carries the six keys the
    /// reference UTR decides, and `alt_type=uORF` where a stop codon past the
    /// point lands in the start codon's frame once moved by the length, which
    /// holds whatever is inserted (UTRAnnotator.pm:973-974). `alt_type_length`
    /// is read from the mutated UTR followed by the coding sequence
    /// (UTRAnnotator.pm:982-987), which the inserted bases change, and is
    /// never written. The same insertion given with bases agrees on every
    /// pair the symbolic row carries.
    #[test]
    fn a_symbolic_insertion_gets_the_frameshift_its_length_decides() {
        for strand in [1, -1] {
            // The UTR's one stop codon ends the uORF: moved by four bases it
            // leaves the frame of the start codon at 5.
            let fixture = Fixture::new(strand, &[]);
            let out = symbolic_insertion(&fixture, 8, 4, |_| {});
            assert_eq!(out[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG, "strand {strand}");
            assert_eq!(
                out[ANNOTATION_KEY],
                "Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:\
                 ref_StartDistanceToCDS=45:ref_type=uORF:ref_type_length=18",
                "strand {strand}"
            );
            assert_eq!(counts_of(&out), ("1", "1", "1"), "strand {strand}");
            let explicit = utr_row(&fixture, 8, "-", "CCCC");
            assert_eq!(explicit[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG);
            assert_eq!(
                out[ANNOTATION_KEY],
                without_keys(&explicit[ANNOTATION_KEY], &["alt_type", "alt_type_length"]),
                "strand {strand}"
            );
            // A record that writes its position, or its position plus its
            // length, as its end states the same insertion.
            let position_as_end = symbolic_insertion(&fixture, 8, 4, |record| {
                record.written_end = WrittenEnd::Position;
            });
            let span_as_end = symbolic_insertion(&fixture, 8, 4, |record| {
                record.written_end = WrittenEnd::SpanEnd;
            });
            assert_eq!(position_as_end, out, "strand {strand}");
            assert_eq!(span_as_end, out, "strand {strand}");

            // The stop codon at 43, moved by one base, ends the frame of the
            // start codon at 5.
            let late_stop = Fixture::with_utr(strand, &utr_with_late_stop(), &[]);
            let out = symbolic_insertion(&late_stop, 8, 1, |_| {});
            assert_eq!(out[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG, "strand {strand}");
            assert_eq!(
                out[ANNOTATION_KEY],
                "Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:alt_type=uORF:\
                 ref_StartDistanceToCDS=45:ref_type=uORF:ref_type_length=18",
                "strand {strand}"
            );
            let explicit = utr_row(&late_stop, 8, "-", "C");
            assert_eq!(
                out[ANNOTATION_KEY],
                without_keys(&explicit[ANNOTATION_KEY], &["alt_type_length"]),
                "strand {strand}"
            );

            // Two start codons' frames hold the point before index 41: the
            // moved stop codon ends the frame of the one at 38 and not of the
            // one at 30.
            let out = symbolic_insertion(&late_stop, 41, 1, |_| {});
            assert_eq!(out[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG, "strand {strand}");
            let explicit = utr_row(&late_stop, 41, "-", "C");
            let explicit: Vec<&str> = explicit[ANNOTATION_KEY].split(',').collect();
            assert_eq!(explicit.len(), 2, "strand {strand}");
            assert!(
                explicit[0].contains(":alt_type=InFrame_oORF:"),
                "{explicit:?}"
            );
            assert!(explicit[1].contains(":alt_type=uORF:"), "{explicit:?}");
            assert_eq!(
                out[ANNOTATION_KEY],
                format!(
                    "{},{}",
                    without_keys(explicit[0], &["alt_type", "alt_type_length"]),
                    without_keys(explicit[1], &["alt_type_length"])
                ),
                "strand {strand}"
            );
            assert!(out[ANNOTATION_KEY].contains("ref_type=OutOfFrame_oORF:"));
            assert!(out[ANNOTATION_KEY].contains(":alt_type=uORF:"));
            assert!(!out[ANNOTATION_KEY].contains("alt_type_length"));
        }
    }

    /// A structural deletion that starts inside the start codon at 5 and runs to
    /// the UTR's end crosses from the first UTR exon into the second, so with
    /// neither exon dropped under the default `max_overlap` its row gets the
    /// reference UTR's counts and `-` alone.
    #[test]
    fn a_structural_variant_row_gets_the_counts_alone() {
        let fixture = Fixture::new(1, &[]);
        let out = symbolic_row(
            &fixture,
            VariantClass::StructuralDeletion,
            (107, 170),
            |_| {},
        );
        assert_eq!(out[CONSEQUENCE_KEY], NONE);
        assert_eq!(out[ANNOTATION_KEY], NONE);
        assert_eq!(counts_of(&out), ("1", "1", "1"));
    }

    /// Every other symbolic allele gets the reference UTR's counts and `-`
    /// alone. A deletion's span must lie inside one UTR exon segment: one that
    /// runs from the first segment into the second, or ends in the intron, in
    /// the coding sequence or past the transcript, or starts before the
    /// transcript, has no sequence deletion of its span in the spliced UTR. A
    /// record that states no end, the `IMPRECISE` flag, a confidence interval
    /// around its position or its end, or several alternates states no single
    /// span. An insertion gets no term when its length is a multiple of three,
    /// when its point lies in no reading frame or outside the UTR, when its
    /// record marks its place uncertain, states no length, writes another end
    /// or has several alternates, and when its text names a mobile element or
    /// another kind. A loss alternate of several gets no term when its record
    /// writes no `END`, and the other alternates of its record never do. A
    /// duplication other than a tandem one does not say where its copy lies,
    /// and no duplication or other class is analysed. An inversion gets no
    /// term when it marks its place uncertain, has several alternates or
    /// leaves the segment. The plugin is handed a structural variant at all because
    /// UTRAnnotator.pm:80-82 lists `StructuralVariationFeature` among its
    /// `variant_feature_types`.
    #[test]
    fn every_other_symbolic_allele_gets_the_counts_alone() {
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            assert!(fixture.plugin.annotates_structural_variants());
            let deletion = VariantClass::StructuralDeletion;
            let inside = fixture.span(8, 9);
            let point = (inside.1 + 1, inside.1);
            // The last two bases of the start codon at 5.
            let start_codon = fixture.span(6, 7);
            // Each span starts in the first UTR exon unless it says otherwise.
            let (into_intron, into_coding, past_transcript, before_transcript) = if strand < 0 {
                ((250, 265), (220, 265), (100, 265), (290, 310))
            } else {
                ((136, 150), (136, 180), (136, 300), (95, 110))
            };
            // A plain insertion of four bases, as its record's span.
            let insertion = VariantClass::StructuralInsertion;
            let (in_intron, in_coding) = if strand < 0 {
                ((251, 254), (201, 204))
            } else {
                ((146, 149), (181, 184))
            };
            let plain_insertion = |record: &mut InputVariant| {
                record.alt_alleles = vec![b"<INS>".to_vec()];
                record.sv_len = Some(4);
            };
            let rows = [
                (
                    "a deletion from one segment into the next",
                    symbolic_row(&fixture, deletion, fixture.span(6, 49), |_| {}),
                ),
                (
                    "a deletion that ends in the intron",
                    symbolic_row(&fixture, deletion, into_intron, |_| {}),
                ),
                (
                    "a deletion that ends in the coding sequence",
                    symbolic_row(&fixture, deletion, into_coding, |_| {}),
                ),
                (
                    "a deletion that ends past the transcript",
                    symbolic_row(&fixture, deletion, past_transcript, |_| {}),
                ),
                (
                    "a deletion that starts before the transcript",
                    symbolic_row(&fixture, deletion, before_transcript, |_| {}),
                ),
                (
                    "a deletion with no stated end",
                    symbolic_row(&fixture, deletion, point, |_| {}),
                ),
                (
                    "a deletion with a confidence interval around its position",
                    symbolic_row(&fixture, deletion, inside, |record| {
                        record.ci_pos = Some((-5, 5));
                    }),
                ),
                (
                    "a deletion with a confidence interval around its end",
                    symbolic_row(&fixture, deletion, inside, |record| {
                        record.ci_end = Some((0, 3));
                    }),
                ),
                (
                    "a deletion marked imprecise with no confidence interval",
                    symbolic_row(&fixture, deletion, inside, |record| {
                        record.imprecise = true;
                    }),
                ),
                (
                    "a loss alternate of several whose record writes no end",
                    symbolic_row(
                        &fixture,
                        VariantClass::CopyNumberVariation,
                        inside,
                        |record| {
                            record.alt_alleles = vec![b"<CN0>".to_vec()];
                            record.uploaded_allele_string = Some("<CN0>/<CN2>".into());
                        },
                    ),
                ),
                (
                    "a loss alternate of several whose record writes another end",
                    symbolic_row(
                        &fixture,
                        VariantClass::CopyNumberVariation,
                        inside,
                        |record| {
                            record.alt_alleles = vec![b"<CN0>".to_vec()];
                            record.uploaded_allele_string = Some("<CN0>/<CN2>".into());
                            record.written_end = WrittenEnd::Other;
                        },
                    ),
                ),
                (
                    "a loss alternate of several marked imprecise",
                    symbolic_row(
                        &fixture,
                        VariantClass::CopyNumberVariation,
                        inside,
                        |record| {
                            record.alt_alleles = vec![b"<CN0>".to_vec()];
                            record.uploaded_allele_string = Some("<CN0>/<CN2>".into());
                            record.written_end = WrittenEnd::SpanEnd;
                            record.imprecise = true;
                        },
                    ),
                ),
                (
                    "the gain alternate beside a loss alternate",
                    symbolic_row(
                        &fixture,
                        VariantClass::CopyNumberVariation,
                        inside,
                        |record| {
                            record.alt_alleles = vec![b"<CN2>".to_vec()];
                            record.uploaded_allele_string = Some("<CN0>/<CN2>".into());
                            record.written_end = WrittenEnd::SpanEnd;
                        },
                    ),
                ),
                (
                    "a copy-number alternate that is neither a loss nor a gain",
                    symbolic_row(
                        &fixture,
                        VariantClass::CopyNumberVariation,
                        inside,
                        |record| {
                            record.alt_alleles = vec![b"<CN1>".to_vec()];
                            record.uploaded_allele_string = Some("<CN0>/<CN1>".into());
                            record.written_end = WrittenEnd::SpanEnd;
                        },
                    ),
                ),
                (
                    "an insertion whose length is a multiple of three",
                    symbolic_insertion(&fixture, 10, 3, |_| {}),
                ),
                (
                    "an insertion inside a start codon",
                    symbolic_insertion(&fixture, 7, 4, |_| {}),
                ),
                (
                    "an insertion before a uORF's stop codon",
                    symbolic_insertion(&fixture, 20, 4, |_| {}),
                ),
                (
                    "an insertion between two reading frames",
                    symbolic_insertion(&fixture, 25, 4, |_| {}),
                ),
                (
                    "an insertion in the intron",
                    symbolic_row(&fixture, insertion, in_intron, plain_insertion),
                ),
                (
                    "an insertion in the coding sequence",
                    symbolic_row(&fixture, insertion, in_coding, plain_insertion),
                ),
                (
                    "an insertion with a confidence interval around its position",
                    symbolic_insertion(&fixture, 10, 4, |record| {
                        record.ci_pos = Some((-5, 5));
                    }),
                ),
                (
                    "an insertion marked imprecise",
                    symbolic_insertion(&fixture, 10, 4, |record| {
                        record.imprecise = true;
                    }),
                ),
                (
                    "an insertion without a stated length",
                    symbolic_row(&fixture, insertion, point, |record| {
                        record.alt_alleles = vec![b"<INS>".to_vec()];
                    }),
                ),
                (
                    "an insertion whose written end is neither its position nor its span's end",
                    symbolic_insertion(&fixture, 10, 4, |record| {
                        record.written_end = WrittenEnd::Other;
                    }),
                ),
                (
                    "an insertion of a record with several alternates",
                    symbolic_insertion(&fixture, 10, 4, |record| {
                        record.uploaded_allele_string = Some("<INS>/<DUP>".into());
                    }),
                ),
                (
                    "a mobile-element insertion",
                    symbolic_insertion(&fixture, 10, 4, |record| {
                        record.variant_class = VariantClass::MobileElementInsertion;
                        record.alt_alleles = vec![b"<INS:ME:ALU>".to_vec()];
                    }),
                ),
                (
                    "an insertion of another text",
                    symbolic_insertion(&fixture, 10, 4, |record| {
                        record.alt_alleles = vec![b"<INS:MT>".to_vec()];
                    }),
                ),
                (
                    "a duplication",
                    symbolic_row(&fixture, VariantClass::Duplication, inside, |record| {
                        record.sv_len = Some(2);
                    }),
                ),
                (
                    "a tandem duplication",
                    symbolic_row(&fixture, VariantClass::TandemDuplication, inside, |_| {}),
                ),
                (
                    "an inversion marked imprecise",
                    symbolic_row(&fixture, VariantClass::Inversion, start_codon, |record| {
                        record.imprecise = true;
                    }),
                ),
                (
                    "an inversion with a confidence interval around its end",
                    symbolic_row(&fixture, VariantClass::Inversion, start_codon, |record| {
                        record.ci_end = Some((-2, 2));
                    }),
                ),
                (
                    "an inversion of a record with several alternates",
                    symbolic_row(&fixture, VariantClass::Inversion, start_codon, |record| {
                        record.uploaded_allele_string = Some("<INV>/<DUP>".into());
                    }),
                ),
                (
                    "an inversion from one segment into the next",
                    symbolic_row(
                        &fixture,
                        VariantClass::Inversion,
                        fixture.span(6, 49),
                        |_| {},
                    ),
                ),
                (
                    "an inversion that ends in the intron",
                    symbolic_row(&fixture, VariantClass::Inversion, into_intron, |_| {}),
                ),
                (
                    "an inversion with no stated end",
                    symbolic_row(&fixture, VariantClass::Inversion, point, |_| {}),
                ),
                (
                    "a copy-number variation",
                    symbolic_row(&fixture, VariantClass::CopyNumberVariation, inside, |_| {}),
                ),
                (
                    "a tandem repeat",
                    symbolic_row(&fixture, VariantClass::TandemRepeat, inside, |_| {}),
                ),
                (
                    "a breakpoint",
                    symbolic_row(&fixture, VariantClass::Translocation, inside, |_| {}),
                ),
                (
                    "a complex alteration",
                    symbolic_row(&fixture, VariantClass::ComplexStructural, inside, |_| {}),
                ),
            ];
            for (shape, out) in rows {
                assert_eq!(out[CONSEQUENCE_KEY], NONE, "strand {strand}: {shape}");
                assert_eq!(out[ANNOTATION_KEY], NONE, "strand {strand}: {shape}");
                assert_eq!(counts_of(&out), ("1", "1", "1"), "strand {strand}: {shape}");
            }
            // The sequence alleles at the same place keep their rows.
            let deleted = utr_row(&fixture, 8, &UTR[8..=9], "-");
            assert_eq!(
                deleted[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                "strand {strand}"
            );
            let inserted = utr_row(&fixture, 8, "-", "CCCC");
            assert_eq!(
                inserted[CONSEQUENCE_KEY], UFRAMESHIFT_FLAG,
                "strand {strand}"
            );
            assert!(
                inserted[ANNOTATION_KEY].contains(":alt_type="),
                "strand {strand}: {}",
                inserted[ANNOTATION_KEY]
            );
        }
    }

    /// Every way of writing the insertion that turns `utr` into `mutated`: the
    /// index and the inserted bases of each record that gives that sequence.
    fn writings(utr: &[u8], mutated: &[u8]) -> Vec<(usize, String)> {
        let length = mutated.len() - utr.len();
        (0..utr.len())
            .filter(|&index| {
                utr[..index] == mutated[..index] && utr[index..] == mutated[index + length..]
            })
            .map(|index| {
                (
                    index,
                    String::from_utf8(mutated[index..index + length].to_vec()).unwrap(),
                )
            })
            .collect()
    }

    /// Every sequence of `length` bases.
    fn base_strings(length: usize) -> Vec<String> {
        let mut out = vec![String::new()];
        for _ in 0..length {
            out = out
                .into_iter()
                .flat_map(|prefix| ["A", "C", "G", "T"].map(|base| format!("{prefix}{base}")))
                .collect();
        }
        out
    }

    /// One mutated sequence has one row however a record writes it. For every
    /// start codon of the fixture UTRs and every insertion of one to six
    /// bases at one of its three indices or directly after it, every record
    /// that gives the same mutated sequence ([`writings`]) carries the same
    /// five keys. Lengths of up to three bases are read on both strands.
    #[test]
    fn one_changed_sequence_at_a_start_codon_has_one_row_however_it_is_written() {
        let utrs = [
            UTR.to_string(),
            utr_with_late_stop(),
            UTR_OVERLAPS.to_string(),
            UTR_OVERLAPS.replacen("CCG", "ATG", 1),
            utr_with(UTR_OVERLAPS, 14..16, "AT"),
            utr_with(UTR_OVERLAPS, 15..16, "A"),
        ];
        let mut groups_read = 0usize;
        let mut several_rows: Vec<String> = Vec::new();
        for (number, utr) in utrs.iter().enumerate() {
            let bytes = utr.as_bytes();
            for (strand, lengths) in [(1, 1..=6), (-1, 1..=3)] {
                let fixture = Fixture::with_utr(strand, utr, &[]);
                let mut seen: HashSet<Vec<u8>> = HashSet::new();
                let mut groups: Vec<Vec<(usize, String)>> = Vec::new();
                for start_pos in atg_positions(bytes) {
                    for index in start_pos..=start_pos + 3 {
                        for length in lengths.clone() {
                            for bases in base_strings(length) {
                                let mutated = mutated_utr(bytes, index, 0, bases.as_bytes());
                                if !seen.insert(mutated.clone()) {
                                    continue;
                                }
                                let ways = writings(bytes, &mutated);
                                if ways.len() > 1 {
                                    groups.push(ways);
                                }
                            }
                        }
                    }
                }
                let records: Vec<(usize, String)> = groups.iter().flatten().cloned().collect();
                let rows = fixture.insertion_rows(&records);
                let mut next = 0usize;
                for ways in &groups {
                    let group = &rows[next..next + ways.len()];
                    next += ways.len();
                    if group.iter().any(|row| row != &group[0]) {
                        let indices: Vec<usize> = ways.iter().map(|(index, _)| *index).collect();
                        several_rows.push(format!(
                            "UTR {number}, strand {strand}, {} bases at {indices:?}",
                            ways[0].1.len()
                        ));
                    }
                }
                groups_read += groups.len();
            }
        }
        // The fixtures hold over a hundred thousand sequences that two or more
        // records give; a count far below it means the enumeration lost them.
        assert!(groups_read > 100_000, "{groups_read} groups");
        assert!(
            several_rows.is_empty(),
            "{} of {groups_read} mutated sequences have several rows, the first: {:?}",
            several_rows.len(),
            &several_rows[..several_rows.len().min(10)]
        );
    }

    /// An inversion that states one place inside one UTR exon segment gives
    /// the mutated sequence whole: its span replaced by the reverse complement
    /// of the bases the UTR holds there. Indices 6 and 7, the last two bases
    /// of the start codon at 5, read `CA` inverted: the codon is lost. Indices
    /// 23 and 24 read `TG` inverted, which the `A` at 22 completes to a start
    /// codon: one is gained there. Each of the two rows is that of the same
    /// bases given with explicit alleles. Indices 4 to 7 read the same four
    /// bases inverted: nothing changes, and the row reports nothing.
    #[test]
    fn an_inversion_inside_one_utr_exon_is_read_as_its_changed_sequence() {
        let lost_start = whole_row(
            UAUG_LOST_FLAG,
            "CapDistanceToStart=5:DistanceToCDS=45:DistanceToStop=15:Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:type=uORF",
            ("1", "1", "1"),
        );
        let gained_start = whole_row(
            UAUG_GAINED_FLAG,
            "CapDistanceToStart=22:DistanceToCDS=28:DistanceToStop=NA:KozakContext=CTGATGC:KozakStrength=Weak:type=OutOfFrame_oORF",
            ("1", "1", "1"),
        );
        let nothing = whole_row(NONE, NONE, ("1", "1", "1"));
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            let inversion = |first: usize, last: usize| {
                symbolic_row(
                    &fixture,
                    VariantClass::Inversion,
                    fixture.span(first, last),
                    |record| record.alt_alleles = vec![b"<INV>".to_vec()],
                )
            };
            for ((first, last), inverted, row) in
                [((6, 7), "CA", &lost_start), ((23, 24), "TG", &gained_start)]
            {
                let out = inversion(first, last);
                assert_eq!(&out, row, "strand {strand}, indices {first} to {last}");
                assert_eq!(
                    out,
                    utr_row(&fixture, first, &UTR[first..=last], inverted),
                    "strand {strand}, indices {first} to {last}"
                );
            }
            assert_eq!(inversion(4, 7), nothing, "strand {strand}");
        }
    }

    /// A record with several alternates is handed to the plugin once for each,
    /// under the record's class. The alternate whose own text is a deletion's
    /// states the loss of the record's span when the record writes `END` as
    /// that span's last base: its row is that of the symbolic deletion and of
    /// the same two bases deleted with explicit alleles, the frameshift of the
    /// uORF at 5, which meets no stop codon in the mutated UTR.
    #[test]
    fn a_loss_alternate_of_several_is_read_as_the_deletion_of_its_span() {
        let shifted = whole_row(
            UFRAMESHIFT_FLAG,
            "Evidence=NA:KozakContext=GCCATGG:KozakStrength=Strong:alt_type=OutOfFrame_oORF:alt_type_length=NA:ref_StartDistanceToCDS=45:ref_type=uORF:ref_type_length=18",
            ("1", "1", "1"),
        );
        for strand in [1, -1] {
            let fixture = Fixture::new(strand, &[]);
            for (own_text, alternates) in [
                ("<CN0>", "<CN0>/<CN2>"),
                ("<CN0>", "<CN2>/<CN0>"),
                ("<DEL>", "<DEL>/<DUP>"),
                ("<CN=0>", "<CN=0>/<CN=2>"),
            ] {
                let out = symbolic_row(
                    &fixture,
                    VariantClass::CopyNumberVariation,
                    fixture.span(8, 9),
                    |record| {
                        record.alt_alleles = vec![own_text.as_bytes().to_vec()];
                        record.uploaded_allele_string = Some(alternates.into());
                        record.written_end = WrittenEnd::SpanEnd;
                    },
                );
                assert_eq!(out, shifted, "strand {strand}, {own_text} of {alternates}");
            }
            assert_eq!(
                symbolic_deletion(&fixture, 8, 9),
                shifted,
                "strand {strand}"
            );
            assert_eq!(
                utr_row(&fixture, 8, &UTR[8..=9], "-"),
                shifted,
                "strand {strand}"
            );
        }
    }

    #[test]
    fn header_info_names_the_five_fields() {
        let names: Vec<String> = UtrAnnotatorPlugin::new()
            .header_info()
            .into_iter()
            .map(|(name, _)| name)
            .collect();
        assert_eq!(
            names,
            [
                CONSEQUENCE_KEY,
                ANNOTATION_KEY,
                EXISTING_UORFS_KEY,
                EXISTING_OUT_OF_FRAME_KEY,
                EXISTING_IN_FRAME_KEY
            ]
        );
    }
}
