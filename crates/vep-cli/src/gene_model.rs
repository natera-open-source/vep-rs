// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Gene-model files read beside the cache: `--gff <file>`, `--gtf <file>` and
//! `--custom file=<file>,format=gff|gtf` in transcript mode, any number of them
//! in one run.
//!
//! Ensembl VEP builds one annotation source per file after the cache's
//! (`AnnotationSourceAdaptor.pm:94-105`); each source annotates the input with
//! its own transcripts (`AnnotationType/Transcript.pm:92-190`) and stamps them
//! with its short name (`AnnotationSource/File/BaseGXF.pm:635`), which the
//! `SOURCE` column prints (`OutputFactory.pm:1592`). Here each file is read
//! once, whole, by the cache builder's reader
//! ([`vep_cache_builder::gene_model`]); its transcripts join the cache's in the
//! run's transcript index and nothing is written to the cache.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use tracing::{info, warn};
use vep_builtins::custom::{parse_custom, CustomFormat, CustomSpec};
use vep_cache_builder::gene_model::{self as reader, Format};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::transcript::{SourceCache, Transcript};

use crate::config::{CliError, GFF_GTF_RECIPE};
use crate::transcript_index::TranscriptSet;

/// The refusal of a GFF or GTF source in a run with no sequence to read. It
/// names the FASTA, the one remedy here, and then Ensembl VEP's text
/// (`BaseGXF.pm:105-106`, less its `ERROR: ` prefix). That text also offers
/// database access, which is no remedy in this program: `--database` is
/// refused and `--cache` reads a local cache without a database.
pub(crate) const FASTA_REQUIRED: &str =
    "a GFF or GTF transcript source needs a FASTA file (--fasta) for its sequences: vep-rs never \
     connects to an Ensembl database (Ensembl VEP: GXF annotation requires either database \
     access (--database or --cache) or a FASTA file (--fasta))";

/// Ensembl VEP's warning for `--nearest` beside a GFF or GTF source
/// (`BaseGXF.pm:113-115`), less the `WARNING: ` prefix the warnings file adds.
pub(crate) const NEAREST_WARNING: &str =
    "--nearest is not currently compatible with GXF annotation sources";

/// One gene-model file of the run: its `--custom` specification and its place
/// among all the run's custom sources, counted from 0 in argument order, where
/// the output headers list it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeneModelFile {
    pub(crate) spec: CustomSpec,
    pub(crate) position: usize,
}

/// `specs` in argument order, split into the sources that annotate records
/// (the overlay sources, opened as custom annotation) and the gene-model
/// files, each keeping its place in `specs`.
pub(crate) fn split(specs: Vec<CustomSpec>) -> (Vec<CustomSpec>, Vec<GeneModelFile>) {
    let mut overlays = Vec::new();
    let mut files = Vec::new();
    for (position, spec) in specs.into_iter().enumerate() {
        if spec.is_transcript_source() {
            files.push(GeneModelFile { spec, position });
        } else {
            overlays.push(spec);
        }
    }
    (overlays, files)
}

/// Whether a `--custom` string of the run names a gene-model file. A string
/// the grammar refuses names none: its own refusal is the one to report.
pub(crate) fn named_in(custom: &[String]) -> bool {
    custom
        .iter()
        .any(|argument| parse_custom(argument).is_ok_and(|spec| spec.is_transcript_source()))
}

/// The refusal of a run that names a gene-model file where no cache stands at
/// the place the run names: the file is read beside a cache, so the text
/// carries the cache's own message and both routes, naming the cache or
/// building one from the file.
pub(crate) fn without_a_cache(cache: &CliError) -> CliError {
    CliError::new(
        "custom",
        format!(
            "a GFF or GTF transcript source is read beside a cache, and this run has none ({}); \
             name the cache of the run with --json_cache <dir>, or with --cache and --dir_cache, \
             or {GFF_GTF_RECIPE}",
            cache.message
        ),
    )
}

/// The refusals a gene-model file draws before it is read, in argument order:
/// the one for a run with no sequence to read ([`FASTA_REQUIRED`], raised
/// where Ensembl VEP raises its own, as each source is built,
/// `BaseGXF.pm:105-106`), then a file that is not there. `fasta_named` is
/// whether the run names or finds a FASTA.
pub(crate) fn check(files: &[GeneModelFile], fasta_named: bool) -> Result<(), CliError> {
    for file in files {
        if !fasta_named {
            return Err(CliError::new("custom", FASTA_REQUIRED));
        }
        std::fs::metadata(&file.spec.file).map_err(|error| {
            CliError::new(
                "custom",
                format!(
                    "cannot open the --custom source {}: {error}",
                    file.spec.file
                ),
            )
        })?;
    }
    Ok(())
}

/// What a run reads a gene-model file's sequences from: its FASTA, where it
/// names or finds one, and the assembly it annotates.
#[derive(Clone, Copy)]
pub(crate) struct Reference<'a> {
    pub(crate) fasta: Option<&'a vep_fasta::IndexedFasta>,
    pub(crate) assembly: Option<&'a str>,
}

/// Whether `transcript` was read from a gene-model file of the run.
pub(crate) fn is_from_file(transcript: &Transcript) -> bool {
    transcript
        .source_cache
        .as_ref()
        .is_some_and(|source| source.file_order().is_some())
}

/// The transcripts of the run's gene-model files, by chromosome: each file
/// read whole and in argument order, its transcripts stamped with the source's
/// short name and the file's place among the run's gene-model files
/// ([`SourceCache::GeneModelFile`]). A sequence region is keyed by the name the
/// cache holds for it (`cache_holds`, through the cache's synonyms as
/// `BaseVEP.pm:697-734` resolves a source's names), else by the file's own.
///
/// `set` holds the run's GENCODE flags, which keep a file's tagged
/// transcripts alone ([`TranscriptSet::accepts_file_transcript`]).
///
/// Refused as Ensembl VEP refuses it: any file without a FASTA
/// ([`FASTA_REQUIRED`]). `--nearest` earns each file Ensembl VEP's warning
/// ([`NEAREST_WARNING`]), and each feature type the reader skips its warning.
/// A file the reader refuses ends the run naming the file, and the line where
/// the refusal is one line's.
///
/// On GRCh37 ([`Reference::assembly`]) a file's `M` or `chrM` is left
/// out, with one warning a region ([`reader::GeneModel::set_apart`]),
/// where the run's FASTA gives `MT` a record of another length than the
/// assembly's mitochondrial sequence, or gives none and names an
/// organelle genome: the cache's synonyms would key it to `MT`, which
/// is not the sequence its models count on. A
/// mitochondrial region whose coding transcripts do not read as coding
/// sequences on the FASTA record they are cut from is left out the same way,
/// with one warning a region ([`reader::GeneModel::off_record`]): a record of
/// the assembly's length can hold the other sequence's bases.
pub(crate) fn load(
    files: &[GeneModelFile],
    reference: Reference<'_>,
    synonyms: &ChromosomeSynonyms,
    cache_holds: impl Fn(&str) -> bool,
    set: &TranscriptSet,
    keep_exon_ids: bool,
    nearest: bool,
) -> Result<HashMap<String, Vec<Transcript>>, CliError> {
    let mut by_chr: HashMap<String, Vec<Transcript>> = HashMap::new();
    let Reference { fasta, assembly } = reference;
    check(files, fasta.is_some())?;
    for (index, file) in files.iter().enumerate() {
        let spec = &file.spec;
        let Some(fasta) = fasta else {
            return Err(CliError::new("custom", FASTA_REQUIRED));
        };
        if nearest {
            warn!("{NEAREST_WARNING}");
        }
        let path = Path::new(&spec.file);
        let format = match spec.format {
            CustomFormat::Gtf => Format::Gtf,
            _ => Format::Gff3,
        };
        let unreadable = |reason: String| {
            CliError::new(
                "custom",
                format!("cannot read the --custom source {}: {reason}", spec.file),
            )
        };
        let model = reader::read(path, format, fasta, assembly)
            .map_err(|error| unreadable(format!("{error:#}")))?;
        for (region, count) in &model.set_apart {
            let transcripts = if *count == 1 {
                "transcript"
            } else {
                "transcripts"
            };
            warn!(
                "Leaving out {count} {transcripts} of {} on {region}: {}, and this run's FASTA \
                 does neither",
                spec.file,
                reader::GRCH37_MITOCHONDRIAL_NAMING
            );
        }
        for (region, found) in &model.off_record {
            let transcripts = if found.transcripts == 1 {
                "transcript"
            } else {
                "transcripts"
            };
            warn!(
                "Leaving out {} {transcripts} of {} on {region}: {}, so this run's FASTA does \
                 not hold the sequence the file's models of that region count on",
                found.transcripts,
                spec.file,
                found.finding()
            );
        }
        // Ensembl VEP warns of each line of a type it passes over as it meets
        // the line (`BaseGXF.pm:239-243`); the file is read whole here, so each
        // such type is warned of once.
        for feature_type in &model.skipped_types {
            warn!(
                "Ignoring '{feature_type}' feature_type from {} GFF/GTF file. This feature_type \
                 is not supported in VEP.",
                spec.file
            );
        }
        let source = SourceCache::GeneModelFile {
            short_name: Arc::from(spec.short_name.as_str()),
            order: u32::try_from(index + 1).unwrap_or(u32::MAX),
        };
        let mut count = 0usize;
        for (region, shard) in model.shards {
            let chr = synonyms.source_chr_name(&region, &cache_holds);
            let mut transcripts =
                crate::json_cache::transcripts_of_shard_bytes(&chr, &shard, keep_exon_ids)
                    .map_err(&unreadable)?;
            transcripts.retain(|transcript| set.accepts_file_transcript(transcript));
            for transcript in &mut transcripts {
                transcript.source_cache = Some(source.clone());
            }
            count += transcripts.len();
            by_chr.entry(chr).or_default().extend(transcripts);
        }
        info!(
            "Gene-model source {}: {} ({}), {count} transcripts",
            spec.short_name,
            spec.file,
            spec.format.as_str()
        );
    }
    Ok(by_chr)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A file of the committed corpus whose reference is Ensembl VEP reading
    /// the GTF: the GTF of eight chromosome 21 genes and the FASTA holding
    /// their sequence.
    fn corpus(file: &str) -> String {
        format!(
            "{}/../../tests/golden/116/GRCh38-gtf-builder/{file}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn gtf(short_name: &str) -> CustomSpec {
        parse_custom(&format!(
            "file={},format=gtf,short_name={short_name}",
            corpus("chr21_genes.gtf.gz")
        ))
        .unwrap()
    }

    fn files(specs: Vec<CustomSpec>) -> Vec<GeneModelFile> {
        split(specs).1
    }

    /// The sources that annotate records keep their order, and each gene-model
    /// file its place among all the run's custom sources, where
    /// `get_all_custom` lists it (`AnnotationSourceAdaptor.pm:327-459`).
    #[test]
    fn specs_split_into_overlays_and_gene_model_files_at_their_places() {
        let specs = vec![
            parse_custom("file=a.vcf.gz,format=vcf,short_name=sites").unwrap(),
            parse_custom("file=genes.gff3.gz,format=gff").unwrap(),
            parse_custom("file=promoters.gff.gz,format=gff,gff_type=gencode_promoter").unwrap(),
            parse_custom("file=genes.gtf.gz,format=gtf").unwrap(),
        ];
        let (overlays, gene_models) = split(specs);
        let names = |specs: &[CustomSpec]| -> Vec<String> {
            specs.iter().map(|s| s.short_name.clone()).collect()
        };
        assert_eq!(names(&overlays), ["sites", "GENCODE_PROMOTER"]);
        let places: Vec<(usize, &str)> = gene_models
            .iter()
            .map(|file| (file.position, file.spec.short_name.as_str()))
            .collect();
        assert_eq!(places, [(1, "genes.gff3.gz"), (3, "genes.gtf.gz")]);
    }

    /// A run naming a gene-model file where no cache stands at the place it
    /// names is told both routes, with the cache's own message; a `--custom`
    /// string of another kind, or one the grammar refuses, names no such file.
    #[test]
    fn a_run_without_a_cache_is_told_both_routes() {
        let named = |arguments: &[&str]| {
            named_in(&arguments.iter().map(|a| a.to_string()).collect::<Vec<_>>())
        };
        assert!(named(&[
            "file=a.vcf.gz,format=vcf",
            "file=a.gtf.gz,format=gtf"
        ]));
        assert!(named(&["file=a.gff3.gz,format=gff"]));
        assert!(!named(&["file=a.vcf.gz,format=vcf"]));
        assert!(!named(&[
            "file=a.gff.gz,format=gff,gff_type=gencode_promoter"
        ]));
        assert!(!named(&["file=a.gtf.gz,format=gtf,type=exact"]));
        assert!(!named(&[]));

        let refusal = without_a_cache(&CliError::new("cache", "Cache directory /d not found"));
        assert_eq!(refusal.component, "custom");
        assert_eq!(
            refusal.message,
            format!(
                "a GFF or GTF transcript source is read beside a cache, and this run has none \
                 (Cache directory /d not found); name the cache of the run with --json_cache \
                 <dir>, or with --cache and --dir_cache, or {GFF_GTF_RECIPE}"
            )
        );
    }

    /// `BaseGXF.pm:105-106`: a GFF or GTF source needs a sequence to read, and
    /// an offline run has the FASTA alone to give it. The refusal names the
    /// FASTA and then Ensembl VEP's text, whose two database routes are no
    /// remedy in a run that reads no database.
    #[test]
    fn a_gene_model_file_without_a_fasta_is_refused_naming_the_fasta() {
        let refusal = load(
            &files(vec![gtf("genes")]),
            Reference {
                fasta: None,
                assembly: None,
            },
            &ChromosomeSynonyms::default(),
            |_| true,
            &TranscriptSet::default(),
            false,
            false,
        )
        .unwrap_err();
        assert_eq!(
            refusal.to_string(),
            "ERROR: [custom] a GFF or GTF transcript source needs a FASTA file (--fasta) for its \
             sequences: vep-rs never connects to an Ensembl database (Ensembl VEP: GXF \
             annotation requires either database access (--database or --cache) or a FASTA \
             file (--fasta))"
        );
    }

    /// Every transcript of a file carries the source's short name and the
    /// file's place among the run's gene-model files (`BaseGXF.pm:635` stamps
    /// the short name), under the name the cache holds its chromosome by; no
    /// transcript of one file is merged with another's.
    #[test]
    fn file_transcripts_carry_their_source_and_their_files_place() {
        let fasta = vep_fasta::IndexedFasta::from_path(corpus("reference.fa.gz")).unwrap();
        let by_chr = load(
            &files(vec![gtf("first"), gtf("second")]),
            Reference {
                fasta: Some(&fasta),
                assembly: None,
            },
            &ChromosomeSynonyms::default(),
            |name| name == "chr21",
            &TranscriptSet::default(),
            false,
            false,
        )
        .unwrap();
        assert_eq!(by_chr.keys().collect::<Vec<_>>(), ["chr21"]);
        let transcripts = &by_chr["chr21"];
        assert_eq!(transcripts.len(), 2 * 59);
        for (order, short_name) in [(1, "first"), (2, "second")] {
            let of_file: Vec<&Transcript> = transcripts
                .iter()
                .filter(|t| {
                    t.source_cache
                        == Some(SourceCache::GeneModelFile {
                            short_name: Arc::from(short_name),
                            order,
                        })
                })
                .collect();
            assert_eq!(of_file.len(), 59, "{short_name}");
            assert!(of_file.iter().all(|t| t.chr == "chr21" && is_from_file(t)));
        }
        assert!(!is_from_file(&Transcript::default()));
    }

    /// The refusals come in argument order, the FASTA's before a missing
    /// file's, and a run that names and holds every file draws none.
    #[test]
    fn a_files_refusals_come_before_it_is_read() {
        let present = gtf("genes");
        let missing = parse_custom("file=/nonexistent/genes.gff3,format=gff").unwrap();
        let both = files(vec![present.clone(), missing]);
        assert_eq!(check(&both, false).unwrap_err().message, FASTA_REQUIRED);
        assert!(check(&both, true)
            .unwrap_err()
            .message
            .starts_with("cannot open the --custom source /nonexistent/genes.gff3: "));
        assert!(check(&files(vec![present]), true).is_ok());
        assert!(check(&[], false).is_ok());
    }

    /// `--gencode_basic` and `--gencode_primary` keep a file's tagged
    /// transcripts alone, as they keep the cache's: the tag the reader stores
    /// is the attribute the flag tests. Ensembl VEP's file sources annotate
    /// every transcript under either flag (`file_cache_honours_gencode_basic`).
    #[test]
    fn the_gencode_flags_keep_a_files_tagged_transcripts() {
        let fasta = vep_fasta::IndexedFasta::from_path(corpus("reference.fa.gz")).unwrap();
        let under = |set: TranscriptSet| -> Vec<Transcript> {
            load(
                &files(vec![gtf("genes")]),
                Reference {
                    fasta: Some(&fasta),
                    assembly: None,
                },
                &ChromosomeSynonyms::default(),
                |_| true,
                &set,
                false,
                false,
            )
            .unwrap()
            .remove("21")
            .unwrap_or_default()
        };
        let all = under(TranscriptSet::default());
        assert_eq!(all.len(), 59);
        let basic = under(TranscriptSet {
            gencode_basic: true,
            ..TranscriptSet::default()
        });
        assert_eq!(basic.len(), all.iter().filter(|t| t.gencode_basic).count());
        assert!(!basic.is_empty() && basic.len() < all.len());
        assert!(basic.iter().all(|t| t.gencode_basic));
        let primary = under(TranscriptSet {
            gencode_primary: true,
            ..TranscriptSet::default()
        });
        assert_eq!(
            primary.len(),
            all.iter().filter(|t| t.gencode_primary).count()
        );
        assert!(primary.len() < all.len());
    }

    /// A file that is not there, and one the reader refuses, end the run under
    /// `custom` with the file named.
    #[test]
    fn an_unreadable_gene_model_file_is_refused_by_name() {
        let fasta = vep_fasta::IndexedFasta::from_path(corpus("reference.fa.gz")).unwrap();
        let load_one = |spec: CustomSpec| {
            load(
                &files(vec![spec]),
                Reference {
                    fasta: Some(&fasta),
                    assembly: None,
                },
                &ChromosomeSynonyms::default(),
                |_| true,
                &TranscriptSet::default(),
                false,
                false,
            )
        };
        let missing = load_one(parse_custom("file=/nonexistent/genes.gtf.gz,format=gtf").unwrap())
            .unwrap_err();
        assert!(
            missing.to_string().starts_with(
                "ERROR: [custom] cannot open the --custom source /nonexistent/genes.gtf.gz: "
            ),
            "{missing}"
        );
        let dir = tempfile::tempdir().unwrap();
        let malformed = dir.path().join("malformed.gtf");
        std::fs::write(
            &malformed,
            "21\tsource\tgene\tx\t20\t.\t+\t.\tgene_id \"g\";\n",
        )
        .unwrap();
        let refused =
            load_one(parse_custom(&format!("file={},format=gtf", malformed.display())).unwrap())
                .unwrap_err()
                .to_string();
        assert!(
            refused.starts_with(&format!(
                "ERROR: [custom] cannot read the --custom source {}: ",
                malformed.display()
            )) && refused.contains("line 1"),
            "{refused}"
        );
    }
}
