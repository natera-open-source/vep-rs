// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A GFF3 or GTF read for a run that names it beside a cache (`vep --gff`,
//! `--gtf`, or a `--custom` source of format `gff` or `gtf` in transcript
//! mode): the file's transcript records, built in memory and never written.
//!
//! The reading is the one a cache build from the file and a genome alone
//! takes (`gff3::Reading::GffSource`, `builder::build_transcripts` with the
//! genome as its only source), so a file annotates the same way read beside a
//! cache or built into one.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;

use anyhow::{Context, Result};

use crate::builder::{
    build_transcripts, grch37_mitochondrial_regions_off_their_record, BuildSources,
};
use crate::gff3::{self, GeneRecord, Reading, SequenceNames};
use crate::gtf;

// `builder` is hidden from the documentation, so the type that
// `GeneModel::off_record` holds is offered under this module's path.
pub use crate::builder::OffRecord;

/// The format of a gene-model file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Format {
    Gff3,
    Gtf,
}

/// The codon table of a sequence region of a gene-model file read as a run's
/// transcript source: 2 for a region named `MT`, 1 for every other, the rule
/// of Ensembl VEP's file sources (`BaseGXF.pm` `_add_translation`).
fn source_codon_table(region: &str) -> u8 {
    if region == "MT" {
        2
    } else {
        1
    }
}

/// What the warning of a region kept apart from GRCh37's `MT` says of the
/// naming, in a run ([`GeneModel::set_apart`]) and in a build.
pub const GRCH37_MITOCHONDRIAL_NAMING: &str =
    "on GRCh37 a gene-model file's M or chrM is read as MT only where the FASTA gives for MT the \
     assembly's mitochondrial sequence, 16,569 bases long, or gives no record for MT and names \
     no organelle genome";

/// A gene-model file as a run reads it.
#[derive(Debug)]
#[non_exhaustive]
pub struct GeneModel {
    /// The file's transcripts by sequence region as the reader names it, each
    /// region's as the bytes of one cache shard: a JSON array of transcript
    /// records in position order.
    pub shards: BTreeMap<String, Vec<u8>>,
    /// The feature types of the lines the reader skipped, in name order.
    pub skipped_types: Vec<String>,
    /// The transcripts left out of `shards`, counted by the region the file
    /// writes: those on an `M` or a `chrM` that a GRCh37 run cannot read as
    /// the assembly's `MT` (`gff3::SequenceNames::keeps_apart_from_mt`). Their
    /// sequences would be cut from a record that is not the sequence the
    /// run's `MT` names.
    pub set_apart: BTreeMap<String, usize>,
    /// The mitochondrial regions of a GRCh37 run left out of `shards`
    /// because none of their coding transcripts reads as a coding sequence
    /// on the FASTA record it is cut from, by the region the reader names
    /// (`builder::grch37_mitochondrial_regions_off_their_record`): the record is
    /// not the sequence the file's models count on.
    pub off_record: BTreeMap<String, OffRecord>,
}

/// The gene-model file at `path`, read whole for a run on `assembly`.
/// `genome` supplies every sequence and the name each record's sequence is
/// read under (`gff3::SequenceNames::of_fasta`), as a build from the file and a
/// genome alone takes them. A file the reader refuses, and a transcript whose
/// coding sequence does not sit in its exons, fail the call with the file
/// named, and the line where the refusal is one line's. On GRCh37 a
/// mitochondrial region is left out where the naming keeps it apart from
/// `MT` ([`GeneModel::set_apart`]) and where its coding transcripts do not
/// read on the record they are cut from ([`GeneModel::off_record`]).
pub fn read(
    path: &Path,
    format: Format,
    genome: &vep_fasta::IndexedFasta,
    assembly: Option<&str>,
) -> Result<GeneModel> {
    let names = SequenceNames::of_fasta(genome, assembly);
    let (mut genes, skipped_types): (HashMap<String, GeneRecord>, Vec<String>) = match format {
        Format::Gff3 => gff3::parse_gff3_listing_skipped(path, Reading::GffSource, &names)?,
        Format::Gtf => {
            let model = gtf::parse_gtf_model(path, Reading::GffSource, &names)?;
            (model.genes, model.skipped_types)
        }
    };
    // A transcript's region is the sequence its own line names, the one its
    // shard is keyed by, whatever sequence its gene line names.
    let mut set_apart: BTreeMap<String, usize> = BTreeMap::new();
    for gene in genes.values_mut() {
        gene.transcripts.retain(|transcript| {
            let apart = names.keeps_apart_from_mt(&transcript.chr);
            if apart {
                *set_apart.entry(transcript.chr.clone()).or_default() += 1;
            }
            !apart
        });
    }
    let codon_tables: HashMap<String, u8> = genes
        .values()
        .flat_map(|gene| &gene.transcripts)
        .map(|transcript| (transcript.chr.clone(), source_codon_table(&transcript.chr)))
        .collect();
    let sources = BuildSources {
        protein_sequences: None,
        prediction_data: None,
        genome_fasta: Some(genome),
        gtf_tags: None,
        core_dumps: None,
    };
    // Tested before any record is built, so no transcript of a region cut
    // from another sequence than its models count on reaches a shard.
    let off_record =
        grch37_mitochondrial_regions_off_their_record(&genes, &sources, &codon_tables, assembly);
    for gene in genes.values_mut() {
        gene.transcripts
            .retain(|transcript| !off_record.contains_key(&transcript.chr));
    }
    let shards = build_transcripts(&genes, path, &sources, &codon_tables)?
        .into_iter()
        .map(|(region, transcripts)| {
            let shard = serde_json::to_vec(&transcripts)
                .with_context(|| format!("Failed to serialize the transcripts of {region}"))?;
            Ok((region, shard))
        })
        .collect::<Result<_>>()?;
    Ok(GeneModel {
        shards,
        skipped_types,
        set_apart,
        off_record,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_mitochondrion_alone_takes_table_two() {
        assert_eq!(source_codon_table("MT"), 2);
        for region in ["1", "X", "Mito", "chrM", "mt"] {
            assert_eq!(source_codon_table(region), 1, "{region}");
        }
    }

    /// A file of the committed corpus built from a GTF and a genome alone.
    fn corpus(file: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh38-gtf-builder")
            .join(file)
    }

    /// A file read for a run holds the records a cache build from the same
    /// file and genome writes: the corpus's committed cache, built by the
    /// binary from this GTF and this FASTA, shard for shard. The feature types
    /// the reader passes over are listed once each.
    #[test]
    fn a_file_read_for_a_run_holds_the_records_a_cache_build_writes() {
        let genome = vep_fasta::IndexedFasta::from_path(corpus("reference.fa.gz")).unwrap();
        let model = read(&corpus("chr21_genes.gtf.gz"), Format::Gtf, &genome, None).unwrap();
        assert_eq!(
            model.skipped_types,
            ["five_prime_utr", "start_codon", "three_prime_utr"]
        );
        assert_eq!(model.shards.keys().collect::<Vec<_>>(), ["21"]);
        let read: Vec<serde_json::Value> = serde_json::from_slice(&model.shards["21"]).unwrap();
        assert_eq!(read.len(), 59);

        let mut built: Vec<serde_json::Value> = Vec::new();
        let shards = corpus("json_cache/transcripts/21");
        for shard in std::fs::read_dir(shards).unwrap() {
            let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
            built.extend(serde_json::from_str::<Vec<serde_json::Value>>(&text).unwrap());
        }
        let position = |record: &serde_json::Value| {
            (
                record["start"].as_u64().unwrap(),
                record["end"].as_u64().unwrap(),
                record["stable_id"].as_str().unwrap().to_string(),
            )
        };
        built.sort_by_key(position);
        assert_eq!(read, built);
    }

    /// A record's sequence is named as a build that reads the same genome
    /// names it ([`SequenceNames::of_fasta`]): a name written with a leading
    /// `chr` that the genome holds without it is read under the genome's
    /// name, where a reading without a genome keeps it as written.
    #[test]
    fn a_files_sequence_is_named_as_the_genome_holds_it() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let residues = "N".repeat(60);
        std::fs::write(&fasta, format!(">scaffold_1\n{residues}\n")).unwrap();
        std::fs::write(
            dir.path().join("genome.fa.fai"),
            "scaffold_1\t60\t12\t60\t61\n",
        )
        .unwrap();
        let genome = vep_fasta::IndexedFasta::from_path(&fasta).unwrap();
        let model = dir.path().join("model.gff3");
        std::fs::write(
            &model,
            "##gff-version 3\n\
chrscaffold_1\tsrc\tncRNA_gene\t5\t50\t.\t+\t.\tID=gene:G1;biotype=lncRNA;gene_id=G1\n\
chrscaffold_1\tsrc\tlnc_RNA\t5\t50\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA;transcript_id=T1\n\
chrscaffold_1\tsrc\texon\t5\t50\t.\t+\t.\tParent=transcript:T1;exon_id=X1;rank=1\n",
        )
        .unwrap();

        let read = read(&model, Format::Gff3, &genome, None).unwrap();
        assert_eq!(read.shards.keys().collect::<Vec<_>>(), ["scaffold_1"]);

        let unnamed =
            gff3::parse_gff3(&model, Reading::GffSource, &SequenceNames::default()).unwrap();
        assert_eq!(unnamed["G1"].chr, "chrscaffold_1");
    }

    /// A transcript's region is the sequence its own line names. On GRCh37,
    /// over a `chrM` record two bases longer than the assembly's
    /// mitochondrial sequence that the synonym table gives for `MT`, a
    /// transcript on `chrM` is set apart under a gene line on `chrMT`, and a
    /// transcript on `chrMT` is read as `MT` under a gene line on `chrM`.
    #[test]
    fn a_transcript_is_set_apart_by_the_sequence_its_own_line_names() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let length = vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH + 2;
        let header = ">chrM\n";
        let mut text = header.to_string();
        let mut left = length;
        while left > 0 {
            let line = left.min(60);
            text.push_str(&"N".repeat(line as usize));
            text.push('\n');
            left -= line;
        }
        std::fs::write(&fasta, text).unwrap();
        std::fs::write(
            dir.path().join("genome.fa.fai"),
            format!("chrM\t{length}\t{}\t60\t61\n", header.len()),
        )
        .unwrap();
        let mut genome = vep_fasta::IndexedFasta::from_path(&fasta).unwrap();
        genome.add_synonyms(&vep_core::synonyms::ChromosomeSynonyms::parse("MT\tchrM\n"));
        genome.set_assembly("GRCh37");
        let model = dir.path().join("model.gff3");
        std::fs::write(
            &model,
            "##gff-version 3\n\
chrMT\tsrc\tncRNA_gene\t5\t50\t.\t+\t.\tID=gene:G1;biotype=lncRNA;gene_id=G1\n\
chrM\tsrc\tlnc_RNA\t5\t50\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA;transcript_id=T1\n\
chrM\tsrc\texon\t5\t50\t.\t+\t.\tParent=transcript:T1;exon_id=X1;rank=1\n\
chrM\tsrc\tncRNA_gene\t60\t90\t.\t+\t.\tID=gene:G2;biotype=lncRNA;gene_id=G2\n\
chrMT\tsrc\tlnc_RNA\t60\t90\t.\t+\t.\tID=transcript:T2;Parent=gene:G2;biotype=lncRNA;transcript_id=T2\n\
chrMT\tsrc\texon\t60\t90\t.\t+\t.\tParent=transcript:T2;exon_id=X2;rank=1\n",
        )
        .unwrap();

        let read = read(&model, Format::Gff3, &genome, Some("GRCh37")).unwrap();
        assert_eq!(read.set_apart, BTreeMap::from([("chrM".to_string(), 1)]));
        assert_eq!(read.shards.keys().collect::<Vec<_>>(), ["MT"]);
        let kept: Vec<serde_json::Value> = serde_json::from_slice(&read.shards["MT"]).unwrap();
        let kept: Vec<&str> = kept
            .iter()
            .map(|record| record["stable_id"].as_str().unwrap())
            .collect();
        assert_eq!(kept, ["T2"]);
    }

    /// A GRCh37 run tests a mitochondrial region's coding transcripts on the
    /// record they are cut from. Over a `chrM` record of the assembly's
    /// length that holds a model's coding bases two positions on, the
    /// region is read as `MT` and left out whole, its non-coding transcript
    /// with it, and the reader states what it found; the same file read for
    /// another assembly keeps the region.
    #[test]
    fn a_region_whose_coding_transcripts_do_not_read_on_their_record_is_left_out() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let length = vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH as usize;
        let mut bases: Vec<u8> = b"ACGT".iter().copied().cycle().take(length).collect();
        bases[12..30].copy_from_slice(b"ATGCTTTGAATAAAATAA");
        let header = ">chrM\n";
        let mut text = header.to_string();
        for line in bases.chunks(60) {
            text.push_str(std::str::from_utf8(line).unwrap());
            text.push('\n');
        }
        std::fs::write(&fasta, text).unwrap();
        std::fs::write(
            dir.path().join("genome.fa.fai"),
            format!("chrM\t{length}\t{}\t60\t61\n", header.len()),
        )
        .unwrap();
        let genome = vep_fasta::IndexedFasta::from_path(&fasta).unwrap();
        let model = dir.path().join("model.gff3");
        std::fs::write(
            &model,
            "##gff-version 3\n\
chrM\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene:G1;biotype=protein_coding;gene_id=G1\n\
chrM\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding;transcript_id=T1\n\
chrM\tsrc\texon\t1\t30\t.\t+\t.\tParent=transcript:T1;exon_id=X1;rank=1\n\
chrM\tsrc\tCDS\t11\t28\t.\t+\t0\tID=CDS:PEPTIDE1;Parent=transcript:T1;protein_id=PEPTIDE1\n\
chrM\tsrc\tncRNA_gene\t60\t90\t.\t+\t.\tID=gene:G2;biotype=lncRNA;gene_id=G2\n\
chrM\tsrc\tlnc_RNA\t60\t90\t.\t+\t.\tID=transcript:T2;Parent=gene:G2;biotype=lncRNA;transcript_id=T2\n\
chrM\tsrc\texon\t60\t90\t.\t+\t.\tParent=transcript:T2;exon_id=X2;rank=1\n",
        )
        .unwrap();

        let on_grch37 = read(&model, Format::Gff3, &genome, Some("GRCh37")).unwrap();
        assert_eq!(on_grch37.set_apart, BTreeMap::new());
        assert_eq!(
            on_grch37.off_record,
            BTreeMap::from([(
                "MT".to_string(),
                OffRecord {
                    transcripts: 2,
                    coding: 1,
                    record: "chrM".to_string(),
                    record_length: vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH,
                    codon_table: 2,
                }
            )])
        );
        assert_eq!(on_grch37.shards.keys().count(), 0);

        let on_grch38 = read(&model, Format::Gff3, &genome, Some("GRCh38")).unwrap();
        assert_eq!(on_grch38.off_record, BTreeMap::new());
        assert_eq!(on_grch38.shards.keys().collect::<Vec<_>>(), ["MT"]);
        let kept: Vec<serde_json::Value> = serde_json::from_slice(&on_grch38.shards["MT"]).unwrap();
        assert_eq!(kept.len(), 2);
    }
}
