// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript source behind HGVS `c.`, `n.` and `p.` input: a stable-id,
//! translation-id and gene-name index over the run's transcript cache,
//! standing in for the core database `fetch_by_hgvs_notation` queries
//! (ensembl-variation release/116, VariationFeatureAdaptor.pm:1829-1849,
//! 1908-1941, 1993-2012).
//!
//! The cache is laid out by chromosome and carries no index by identifier, so
//! the first lookup reads the reference sequences of the run's transcript set
//! (the assembled chromosomes and the unplaced scaffolds) shard by shard,
//! keeps each transcript's identifiers and drops the rest, and indexes them;
//! the alternate loci, patches and LRG regions are read the same way only
//! when a stable id is found on none of the reference sequences. A hit then
//! materializes its chromosome alone, through the same lazy index the
//! annotation reads, so a run with transcript notations parses the cache once
//! more than a run on genomic coordinates and holds the identifiers of every
//! transcript beside the chromosomes its records touch. The lookup reads the
//! run's transcript set as the annotation does (`--gencode_basic` and
//! `--gencode_primary` narrow both), where Ensembl VEP's resolver reads its
//! whole database.
//!
//! A stable id present on several sequences (a pseudoautosomal transcript the
//! cache stores under X and Y) resolves on the first by the cache's contig
//! order, as the database holds the one row the X copy is projected from; a
//! translation id resolves the same way to the transcript it translates. A
//! gene name resolves among the reference sequences alone, as the reference
//! gene is the one `_get_gene_transcripts` keeps (`:1997`, `is_reference`). Of
//! several reference genes carrying one name the first by database order is
//! taken, read here as the one whose transcripts have the lowest database id,
//! and a gene's transcripts come in database-id order, the order
//! `get_all_Transcripts` lists them in.

use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::sync::{Arc, OnceLock};

use rayon::prelude::*;
use tracing::info;
use vep_core::transcript::Transcript;
use vep_io::input::TranscriptSource;

use crate::config::CliError;
use crate::transcript_index::LazyTranscriptIndexes;

/// What the index keeps of a transcript: the identifiers a notation resolves
/// through and the fields the gene pick orders by, with the position of its
/// chromosome in [`TranscriptLookup::chromosomes`].
#[derive(Debug, Clone, PartialEq, Eq)]
struct Entry {
    chromosome: usize,
    stable_id: Arc<str>,
    version: Option<u32>,
    /// The stable id and version of the transcript's translation, the cache's
    /// `protein` id with no version for a transcript whose record carries the
    /// id without its translation; `None` for a transcript that translates
    /// nothing.
    translation_id: Option<Arc<str>>,
    translation_version: Option<u32>,
    gene_stable_id: Arc<str>,
    gene_symbol: Option<Arc<str>>,
    db_id: Option<u64>,
    /// Whether the transcript is a RefSeq model, the cache's `source` for one.
    refseq: bool,
}

impl Entry {
    fn of(transcript: &Transcript, chromosome: usize) -> Self {
        let (translation_id, translation_version) = match &transcript.translation {
            Some(translation) => (
                Some(Arc::from(translation.stable_id.as_str())),
                translation.version,
            ),
            None => (transcript.protein_id.clone(), None),
        };
        Self {
            chromosome,
            stable_id: Arc::clone(&transcript.stable_id),
            version: transcript.version,
            translation_id,
            translation_version,
            gene_stable_id: Arc::clone(&transcript.gene_stable_id),
            gene_symbol: transcript.gene_symbol.clone(),
            db_id: transcript.db_id,
            refseq: transcript.source.eq_ignore_ascii_case("refseq"),
        }
    }

    /// The key the shard loader de-duplicates a chromosome's transcripts by: a
    /// gene's transcripts sit in every shard the gene overlaps.
    fn loader_key(&self) -> (Option<u64>, Arc<str>, Option<u32>) {
        match self.db_id {
            Some(id) => (Some(id), Arc::from(""), None),
            None => (None, Arc::clone(&self.stable_id), self.version),
        }
    }

    /// The entry's place in database order: its database id, a transcript
    /// without one after every one with one, then its stable id and its
    /// gene's, so the order is total and two genes carrying one symbol are
    /// told apart.
    fn database_order(&self) -> (u64, &str, &str) {
        (
            self.db_id.unwrap_or(u64::MAX),
            &self.stable_id,
            &self.gene_stable_id,
        )
    }
}

/// The identifier index over one tier of chromosomes, built once.
struct Index {
    /// Every indexed transcript, in chromosome order and then the cache's.
    entries: Vec<Entry>,
    /// The positions in `entries` of the entries carrying each stable id.
    by_stable_id: HashMap<Arc<str>, Vec<usize>>,
    /// The positions in `entries` of the entries whose translation carries
    /// each stable id.
    by_translation_id: HashMap<Arc<str>, Vec<usize>>,
    /// The positions in `entries` of the entries whose gene carries each symbol.
    by_gene_symbol: HashMap<Arc<str>, Vec<usize>>,
}

impl Index {
    /// `fetch_by_stable_id` within this tier: the entry carrying `stable_id`,
    /// else, when the id ends in `.<version>`, the one carrying the id before
    /// the dot with that version (TranscriptAdaptor.pm:164-181).
    fn find_stable_id(&self, stable_id: &str) -> Option<&Entry> {
        find_versioned(&self.by_stable_id, &self.entries, stable_id, |entry| {
            entry.version
        })
    }

    /// `fetch_by_translation_stable_id` within this tier
    /// (TranscriptAdaptor.pm:307-332): as [`Index::find_stable_id`], over the
    /// translations' ids and versions.
    fn find_translation_id(&self, stable_id: &str) -> Option<&Entry> {
        find_versioned(&self.by_translation_id, &self.entries, stable_id, |entry| {
            entry.translation_version
        })
    }
}

/// The entry `by_id` lists first under `stable_id`, else, when the id ends in
/// `.<version>`, the first under the id before the dot whose `version_of` is
/// that version: the two-step lookup `fetch_by_stable_id` and
/// `fetch_by_translation_stable_id` share (TranscriptAdaptor.pm:164-181,
/// 307-332).
fn find_versioned<'a>(
    by_id: &HashMap<Arc<str>, Vec<usize>>,
    entries: &'a [Entry],
    stable_id: &str,
    version_of: impl Fn(&Entry) -> Option<u32>,
) -> Option<&'a Entry> {
    if let Some(&at) = by_id.get(stable_id).and_then(|all| all.first()) {
        return Some(&entries[at]);
    }
    let (id, version) = stable_id.rsplit_once('.')?;
    if id.is_empty() || !version.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let version: u32 = version.parse().ok()?;
    by_id
        .get(id)?
        .iter()
        .map(|&at| &entries[at])
        .find(|entry| version_of(entry) == Some(version))
}

/// A [`TranscriptSource`] over the run's lazily loaded transcript indexes.
pub struct TranscriptLookup {
    indexes: Arc<LazyTranscriptIndexes>,
    /// The cache's chromosomes in the order identifiers resolve: the reference
    /// sequences first and the alternate ones after, each group in the cache's
    /// contig order and then by name.
    chromosomes: Vec<String>,
    /// How many of `chromosomes` are reference sequences.
    reference_count: usize,
    /// Whether the run's first core group is the RefSeq one (`--refseq`), so
    /// a gene name resolves among the RefSeq transcripts before the Ensembl
    /// ones.
    refseq_first: bool,
    reference: OnceLock<Result<Index, String>>,
    alternate: OnceLock<Result<Index, String>>,
}

impl TranscriptLookup {
    /// A lookup over `indexes`, resolving a stable id present on several
    /// chromosomes on the first by `contig_order` (the cache's, which may be
    /// empty) and then by name, the reference sequences before the others; a
    /// gene name resolves among the RefSeq transcripts first under
    /// `refseq_first`, else among the Ensembl ones first, the other group's
    /// gene listed after it.
    pub fn new(
        indexes: Arc<LazyTranscriptIndexes>,
        contig_order: &[String],
        refseq_first: bool,
    ) -> Self {
        let mut ordered: Vec<String> = contig_order
            .iter()
            .filter(|chr| indexes.contains_key(chr))
            .cloned()
            .collect();
        let mut rest: Vec<String> = indexes
            .chromosomes()
            .filter(|chr| !contig_order.iter().any(|named| named == chr))
            .map(str::to_string)
            .collect();
        rest.sort();
        ordered.extend(rest);
        let (reference, alternate): (Vec<String>, Vec<String>) = ordered
            .into_iter()
            .partition(|chr| is_reference_sequence(chr));
        let reference_count = reference.len();
        let mut chromosomes = reference;
        chromosomes.extend(alternate);
        Self {
            indexes,
            chromosomes,
            reference_count,
            refseq_first,
            reference: OnceLock::new(),
            alternate: OnceLock::new(),
        }
    }

    /// The index over the reference sequences, built on the first call.
    fn reference_index(&self) -> Result<&Index, String> {
        self.reference
            .get_or_init(|| self.build_index(0..self.reference_count, "reference sequences"))
            .as_ref()
            .map_err(Clone::clone)
    }

    /// The index over the alternate loci, patches and LRG regions, built on
    /// the first call.
    fn alternate_index(&self) -> Result<&Index, String> {
        self.alternate
            .get_or_init(|| {
                self.build_index(
                    self.reference_count..self.chromosomes.len(),
                    "alternate sequences",
                )
            })
            .as_ref()
            .map_err(Clone::clone)
    }

    /// The index over `chromosomes[tier]`: a chromosome the annotation has
    /// materialized is read from its index, any other shard by shard, each
    /// shard's transcripts dropped once their identifiers are kept; a shard
    /// that cannot be read fails every lookup through this index with its
    /// reason.
    fn build_index(&self, tier: Range<usize>, what: &str) -> Result<Index, String> {
        info!(
            "HGVS transcript notation: indexing the {what} of the transcript cache ({} chromosomes)",
            tier.len()
        );
        let mut entries = Vec::new();
        let mut by_stable_id: HashMap<Arc<str>, Vec<usize>> = HashMap::new();
        let mut by_translation_id: HashMap<Arc<str>, Vec<usize>> = HashMap::new();
        let mut by_gene_symbol: HashMap<Arc<str>, Vec<usize>> = HashMap::new();
        for chromosome in tier {
            for entry in self.entries_of(chromosome)? {
                let at = entries.len();
                by_stable_id
                    .entry(Arc::clone(&entry.stable_id))
                    .or_default()
                    .push(at);
                if let Some(translation) = &entry.translation_id {
                    by_translation_id
                        .entry(Arc::clone(translation))
                        .or_default()
                        .push(at);
                }
                if let Some(symbol) = &entry.gene_symbol {
                    by_gene_symbol
                        .entry(Arc::clone(symbol))
                        .or_default()
                        .push(at);
                }
                entries.push(entry);
            }
        }
        Ok(Index {
            entries,
            by_stable_id,
            by_translation_id,
            by_gene_symbol,
        })
    }

    /// The entries of `chromosomes[chromosome]`, in the cache's order, each
    /// transcript once (the loader's rule: by database id, else by stable id
    /// and version) and only those the run's transcript set admits.
    fn entries_of(&self, chromosome: usize) -> Result<Vec<Entry>, String> {
        let chr = &self.chromosomes[chromosome];
        if let Some(index) = self.indexes.get(chr) {
            return Ok(index
                .transcripts()
                .iter()
                .map(|transcript| Entry::of(transcript, chromosome))
                .collect());
        }
        let shards = self
            .indexes
            .shards_of(chr)
            .ok_or_else(|| format!("chr{chr} has neither a shard list nor a materialized index"))?;
        let per_shard: Vec<Result<Vec<Entry>, String>> = shards
            .par_iter()
            .map(|shard| {
                let transcripts = self.indexes.load_shard(chr, shard).map_err(error_message)?;
                Ok(transcripts
                    .iter()
                    .filter(|transcript| self.indexes.accepts(transcript))
                    .map(|transcript| Entry::of(transcript, chromosome))
                    .collect())
            })
            .collect();
        let mut seen = HashSet::new();
        let mut entries = Vec::new();
        for shard in per_shard {
            for entry in shard? {
                if seen.insert(entry.loader_key()) {
                    entries.push(entry);
                }
            }
        }
        Ok(entries)
    }

    /// The transcript `entry` stands for, read from its chromosome's index,
    /// materialized for the purpose when the annotation has not.
    fn materialize(&self, entry: &Entry) -> Result<Transcript, String> {
        let chr = &self.chromosomes[entry.chromosome];
        self.indexes.prewarm(chr).map_err(error_message)?;
        self.indexes
            .get(chr)
            .and_then(|index| {
                index.transcripts().iter().find(|transcript| {
                    transcript.stable_id == entry.stable_id
                        && transcript.version == entry.version
                        && transcript.gene_stable_id == entry.gene_stable_id
                })
            })
            .cloned()
            .ok_or_else(|| {
                format!(
                    "transcript {} is indexed on chr{chr} but its materialized index lacks it",
                    entry.stable_id
                )
            })
    }

    /// The entry `find` locates in the reference tier, else, when the cache
    /// has alternate sequences, in the alternate tier, materialized: the
    /// two-tier read every identifier lookup makes.
    fn in_either_tier(
        &self,
        find: impl Fn(&Index) -> Option<&Entry>,
    ) -> Result<Option<Transcript>, String> {
        let mut entry = find(self.reference_index()?);
        let has_alternate_sequences = self.reference_count < self.chromosomes.len();
        if entry.is_none() && has_alternate_sequences {
            entry = find(self.alternate_index()?);
        }
        entry.map(|entry| self.materialize(entry)).transpose()
    }

    /// `_get_gene_transcripts` within one core group, the RefSeq transcripts
    /// under `refseq` and the Ensembl ones otherwise: of the entries at
    /// `carrying_name`, those of the group's gene first in database order, each
    /// transcript once, materialized; empty when the group has no gene of the
    /// name.
    fn gene_in_group(
        &self,
        index: &Index,
        carrying_name: &[usize],
        refseq: bool,
    ) -> Result<Vec<Transcript>, String> {
        let mut entries: Vec<&Entry> = carrying_name
            .iter()
            .map(|&at| &index.entries[at])
            .filter(|entry| entry.refseq == refseq)
            .collect();
        entries.sort_by(|a, b| a.database_order().cmp(&b.database_order()));
        let Some(&lowest) = entries.first() else {
            return Ok(Vec::new());
        };
        // A pseudoautosomal gene is stored under X and Y; each transcript is
        // listed once, as the database holds the one row.
        entries.retain(|entry| entry.gene_stable_id == lowest.gene_stable_id);
        entries.dedup_by(|a, b| a.stable_id == b.stable_id);
        entries
            .into_iter()
            .map(|entry| self.materialize(entry))
            .collect()
    }
}

/// A cache read failure as the parser's `Cannot read the transcript cache`
/// error quotes it: a refusal of the error contract by its message alone,
/// since the parser supplies the line's one `ERROR:` prefix; any other failure
/// with its cause chain.
fn error_message(error: anyhow::Error) -> String {
    match error.downcast_ref::<CliError>() {
        Some(refusal) => refusal.message.clone(),
        None => format!("{error:#}"),
    }
}

/// Whether `chr` names a reference sequence of the assembly: an assembled
/// chromosome or an unplaced scaffold, not an alternate locus (`HSCHR…`,
/// `CHR_…`, `…_CTG…`), a patch (`…_PATCH`, `…_NOVEL…`) or an LRG region.
fn is_reference_sequence(chr: &str) -> bool {
    !(chr.starts_with("HSCHR")
        || chr.starts_with("CHR_")
        || chr.starts_with("LRG_")
        || chr.contains("_CTG")
        || chr.contains("_PATCH")
        || chr.contains("_NOVEL"))
}

impl TranscriptSource for TranscriptLookup {
    fn by_stable_id(&self, stable_id: &str) -> Result<Option<Transcript>, String> {
        self.in_either_tier(|index| index.find_stable_id(stable_id))
    }

    fn by_translation_id(&self, stable_id: &str) -> Result<Option<Transcript>, String> {
        self.in_either_tier(|index| index.find_translation_id(stable_id))
    }

    fn gene_transcripts(&self, name: &str) -> Result<Vec<Vec<Transcript>>, String> {
        let index = self.reference_index()?;
        let carrying_name: &[usize] = index.by_gene_symbol.get(name).map_or(&[], Vec::as_slice);
        // Ensembl VEP tries its core groups in the run's order (Parser/HGVS.pm:164),
        // the RefSeq group first only under --refseq (Config.pm:409).
        [self.refseq_first, !self.refseq_first]
            .into_iter()
            .map(|refseq| self.gene_in_group(index, carrying_name, refseq))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runner::make_runner_test_transcript;
    use crate::transcript_index::TranscriptIndex;

    /// The index keeps identifiers alone: over the corpus cache (chromosomes
    /// 21 and MT) a stable id's lookup leaves every chromosome unmaterialized
    /// but the one the hit is on, and a miss on every reference sequence reads
    /// the alternate tier without materializing anything.
    #[test]
    fn a_lookup_materializes_the_hit_chromosome_alone() {
        let cache = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/golden/116/GRCh38-input-hgvs/json_cache"
        );
        let indexes = Arc::new(LazyTranscriptIndexes::new(cache, 5_000, 5_000).unwrap());
        let lookup = TranscriptLookup::new(Arc::clone(&indexes), &[], false);
        assert!(lookup.by_stable_id("ENST00009999999").unwrap().is_none());
        assert!(lookup.reference.get().is_some());
        assert!(indexes.get("21").is_none() && indexes.get("MT").is_none());
        let found = lookup.by_stable_id("ENST00000291554.6").unwrap().unwrap();
        assert_eq!((found.chr.as_str(), found.version), ("21", Some(6)));
        assert!(indexes.get("21").is_some() && indexes.get("MT").is_none());
        let translated = lookup
            .by_translation_id("ENSP00000291554")
            .unwrap()
            .unwrap();
        assert_eq!(translated.stable_id.as_ref(), "ENST00000291554");
        let gene = lookup.gene_transcripts("TFF3").unwrap();
        assert!(gene[0]
            .iter()
            .all(|t| t.gene_symbol.as_deref() == Some("TFF3")));
        assert!(gene[1].is_empty());
    }

    /// A transcript the cache refuses when its chromosome loads (its mapper's
    /// coding start lies a UTR longer than the span away) fails every lookup
    /// through its tier with the refusal's message alone, no `ERROR: [cache]`
    /// line inside it, so the HGVS parser's wrapping of that message prints
    /// under the one `cache` component.
    #[test]
    fn a_cache_refusal_fails_the_lookup_with_its_message_alone() {
        use vep_io::input::hgvs::HgvsParser;
        use vep_io::input::{InputParser, ParserContext, ReferenceSequence};

        struct NoGenome;
        impl ReferenceSequence for NoGenome {
            fn sequence(&self, _chr: &str, _start: u64, _end: u64) -> Option<Vec<u8>> {
                None
            }
            fn record_name(&self, _chr: &str) -> Option<String> {
                None
            }
            fn sequence_length(&self, _chr: &str) -> Option<u64> {
                None
            }
        }

        let dir = tempfile::Builder::new()
            .prefix("vep_tl_refusal_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let far = u64::MAX;
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!(
                r#"[{{"stable_id":"ENST1","gene_stable_id":"ENSG1","chr":"21","start":1000,"end":1029,
                "strand":1,"biotype":"protein_coding","source":"ensembl","is_canonical":1,
                "exons":[{{"start":1000,"end":1029,"rank":1,"phase":0,"end_phase":0}}],
                "translation":{{"stable_id":"ENSP1","start":1,"end":30}},
                "variation_effect_feature_cache":{{"codon_table":1,
                  "translateable_seq":"ATGGCTGCTGCTGCTGCTGCTGCTGCTTGA","peptide":"MAAAAAAAA",
                  "mapper":{{"start_phase":0,"cdna_coding_start":{far},"cdna_coding_end":30,
                    "pairs":[{{"from_start":1,"from_end":30,"to_start":1000,"to_end":1029,"ori":1}}]}}}}}}]"#
            ),
        )
        .unwrap();
        let indexes =
            LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
        let lookup = Arc::new(TranscriptLookup::new(Arc::new(indexes), &[], false));
        let refusal = format!(
            "transcript ENST1: the mapper's coding start {far} and end 30 on a cDNA of 30 bases \
             leave a UTR longer than the transcript's 30-base span"
        );
        assert_eq!(lookup.by_stable_id("ENST1").unwrap_err(), refusal);
        assert_eq!(lookup.by_translation_id("ENSP1").unwrap_err(), refusal);

        let context = ParserContext {
            fasta: Some(Arc::new(NoGenome)),
            transcripts: Some(Arc::clone(&lookup) as Arc<dyn TranscriptSource>),
            ..Default::default()
        };
        let mut parser = HgvsParser::new(std::io::Cursor::new("ENST1:c.1A>G\n"), &context, false);
        let Some(Err(vep_io::IoError::CacheUnreadable(message))) = parser.next_variant() else {
            panic!("the notation must fail with the cache's refusal");
        };
        assert_eq!(
            CliError::new("cache", message).to_string(),
            format!("ERROR: [cache] Cannot read the transcript cache for HGVS input: {refusal}")
        );
    }

    fn transcript(stable_id: &str, version: u32, db_id: u64, chr: &str, gene: &str) -> Transcript {
        Transcript {
            version: Some(version),
            db_id: Some(db_id),
            gene_symbol: Some(Arc::from("GENE")),
            ..make_runner_test_transcript(stable_id, gene, chr, 1_000)
        }
    }

    fn lookup(by_chr: Vec<(&str, Vec<Transcript>)>, contig_order: &[&str]) -> TranscriptLookup {
        lookup_in_group(by_chr, contig_order, false)
    }

    fn lookup_in_group(
        by_chr: Vec<(&str, Vec<Transcript>)>,
        contig_order: &[&str],
        refseq_first: bool,
    ) -> TranscriptLookup {
        let indexes = by_chr
            .into_iter()
            .map(|(chr, txs)| (chr.to_string(), TranscriptIndex::new(txs, 5_000, 5_000)))
            .collect();
        let order: Vec<String> = contig_order.iter().map(|c| c.to_string()).collect();
        TranscriptLookup::new(
            Arc::new(LazyTranscriptIndexes::from_indexes(indexes)),
            &order,
            refseq_first,
        )
    }

    fn refseq(stable_id: &str, db_id: u64, gene: &str) -> Transcript {
        Transcript {
            source: "RefSeq".to_string(),
            ..transcript(stable_id, 1, db_id, "21", gene)
        }
    }

    /// A stable id finds its transcript whatever the chromosome; one written
    /// with a version finds the transcript of that version after the id as
    /// written finds nothing, and no other version; an id on two reference
    /// chromosomes resolves on the first by contig order, and one on a
    /// reference chromosome and a patch on the reference chromosome whatever
    /// the order; an id on a patch alone is found there, after the reference
    /// sequences miss it.
    #[test]
    fn stable_id_resolves_as_fetch_by_stable_id_does() {
        let lookup = lookup(
            vec![
                ("21", vec![transcript("ENST1", 3, 10, "21", "ENSG1")]),
                ("X", vec![transcript("ENSTPAR", 1, 20, "X", "ENSGP")]),
                ("Y", vec![transcript("ENSTPAR", 1, 20, "Y", "ENSGP")]),
                (
                    "HG1_PATCH",
                    vec![
                        transcript("ENST1", 3, 10, "HG1_PATCH", "ENSG1"),
                        transcript("ENSTPATCH", 1, 40, "HG1_PATCH", "ENSG4"),
                    ],
                ),
            ],
            &["HG1_PATCH", "Y", "X", "21"],
        );
        assert_eq!(lookup.by_stable_id("ENST1").unwrap().unwrap().chr, "21");
        assert_eq!(lookup.by_stable_id("ENST1.3").unwrap().unwrap().chr, "21");
        assert_eq!(lookup.by_stable_id("ENSTPAR").unwrap().unwrap().chr, "Y");
        assert!(lookup.alternate.get().is_none());
        assert!(lookup.by_stable_id("ENST1.4").unwrap().is_none());
        assert!(lookup.alternate.get().is_some());
        assert!(lookup.by_stable_id("ENST1.x").unwrap().is_none());
        assert!(lookup.by_stable_id("ENST1.99999999999").unwrap().is_none());
        assert!(lookup.by_stable_id("ENST9").unwrap().is_none());
        assert_eq!(
            lookup.by_stable_id("ENSTPATCH").unwrap().unwrap().chr,
            "HG1_PATCH"
        );
    }

    /// A translation id finds the transcript it translates, by the translation's
    /// stable id as written or, after that finds nothing, by the id before a
    /// `.<version>` with that version, as `fetch_by_translation_stable_id`
    /// does (TranscriptAdaptor.pm:307-332); a transcript whose record carries
    /// the protein id without its translation is found by that id and by no
    /// version of it; a transcript's own stable id finds nothing here, and a
    /// transcript that translates nothing is never found.
    #[test]
    fn translation_id_resolves_as_fetch_by_translation_stable_id_does() {
        use vep_core::transcript::Translation;
        let translated =
            |stable_id: &str, protein: &str, version: Option<u32>, db_id: u64| Transcript {
                translation: Some(Translation {
                    stable_id: protein.to_string(),
                    version,
                    db_id: None,
                    start: 1,
                    end: 30,
                    start_exon_index: 0,
                    end_exon_index: 0,
                    seq: None,
                }),
                ..transcript(stable_id, 1, db_id, "21", "ENSG1")
            };
        let lookup = lookup(
            vec![(
                "21",
                vec![
                    translated("ENST1", "ENSP1", Some(2), 10),
                    translated("ENST2", "ENSP2", None, 11),
                    Transcript {
                        protein_id: Some(Arc::from("ENSP3")),
                        ..transcript("ENST3", 1, 12, "21", "ENSG1")
                    },
                    transcript("ENST4", 1, 13, "21", "ENSG1"),
                ],
            )],
            &["21"],
        );
        let found = |id: &str| {
            lookup
                .by_translation_id(id)
                .unwrap()
                .map(|t| t.stable_id.to_string())
        };
        assert_eq!(found("ENSP1").as_deref(), Some("ENST1"));
        assert_eq!(found("ENSP1.2").as_deref(), Some("ENST1"));
        assert_eq!(found("ENSP1.3"), None);
        assert_eq!(found("ENSP2").as_deref(), Some("ENST2"));
        assert_eq!(found("ENSP2.1"), None);
        assert_eq!(found("ENSP3").as_deref(), Some("ENST3"));
        assert_eq!(found("ENSP3.1"), None);
        assert_eq!(found("ENST1"), None);
        assert_eq!(found("ENSP4"), None);
    }

    /// A gene name finds the transcripts of the reference gene so named, in
    /// database-id order, a copy on a patch or an alternate locus left out;
    /// two reference genes of one name resolve to the one with the lower
    /// database ids.
    #[test]
    fn gene_name_resolves_among_the_reference_sequences() {
        let lookup = lookup(
            vec![
                (
                    "21",
                    vec![
                        transcript("ENST2", 1, 12, "21", "ENSG1"),
                        transcript("ENST1", 1, 11, "21", "ENSG1"),
                    ],
                ),
                (
                    "HG1_PATCH",
                    vec![transcript("ENST3", 1, 5, "HG1_PATCH", "ENSG3")],
                ),
                (
                    "HSCHR21_1_CTG1",
                    vec![transcript("ENST4", 1, 6, "HSCHR21_1_CTG1", "ENSG4")],
                ),
                (
                    "KI270728.1",
                    vec![transcript("ENST5", 1, 30, "KI270728.1", "ENSG5")],
                ),
            ],
            &["21"],
        );
        let ids: Vec<String> = lookup
            .gene_transcripts("GENE")
            .unwrap()
            .remove(0)
            .iter()
            .map(|t| t.stable_id.to_string())
            .collect();
        assert_eq!(ids, ["ENST1", "ENST2"]);
        assert!(lookup
            .gene_transcripts("OTHER")
            .unwrap()
            .iter()
            .all(Vec::is_empty));
        assert!(is_reference_sequence("GL000008.2"));
        assert!(!is_reference_sequence("LRG_8"));
        assert!(!is_reference_sequence("HG142_HG150_NOVEL_TEST"));
    }

    /// A gene name lists the gene of the run's first core group first, the
    /// Ensembl transcripts in an Ensembl run and the RefSeq ones in a RefSeq
    /// run, then the other group's gene, a group without the gene an empty
    /// list; a gene on both X and Y lists each transcript once.
    #[test]
    fn gene_name_resolves_in_the_first_core_group_and_once_per_transcript() {
        let by_chr = || {
            vec![
                (
                    "21",
                    vec![
                        refseq("NM_1", 1, "GENE_REFSEQ"),
                        transcript("ENST1", 1, 2, "21", "ENSG1"),
                        refseq("NM_2", 3, "GENE_REFSEQ"),
                    ],
                ),
                ("X", vec![transcript("ENSTPAR", 1, 20, "X", "ENSGP")]),
                ("Y", vec![transcript("ENSTPAR", 1, 20, "Y", "ENSGP")]),
            ]
        };
        let ids = |lookup: &TranscriptLookup, name: &str| -> Vec<Vec<String>> {
            lookup
                .gene_transcripts(name)
                .unwrap()
                .iter()
                .map(|group| group.iter().map(|t| t.stable_id.to_string()).collect())
                .collect()
        };
        let ensembl_first = lookup_in_group(by_chr(), &["21"], false);
        assert_eq!(
            ids(&ensembl_first, "GENE"),
            [vec!["ENST1"], vec!["NM_1", "NM_2"]]
        );
        let refseq_first = lookup_in_group(by_chr(), &["21"], true);
        assert_eq!(
            ids(&refseq_first, "GENE"),
            [vec!["NM_1", "NM_2"], vec!["ENST1"]]
        );
        let mut par = by_chr();
        for (_, transcripts) in &mut par {
            for transcript in transcripts.iter_mut() {
                if transcript.gene_stable_id.as_ref() == "ENSGP" {
                    transcript.gene_symbol = Some(Arc::from("PAR"));
                }
            }
        }
        let par = lookup_in_group(par, &["21"], true);
        assert_eq!(ids(&par, "PAR"), [vec![], vec!["ENSTPAR"]]);
    }
}
