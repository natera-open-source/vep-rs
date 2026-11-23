// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A cache built from a GTF or a GFF3 given as the transcript model takes
//! `--transcript_filter` and `--exclude_predicted` as any cache does. So the
//! two-step route keeps the transcripts Ensembl VEP's `--gtf` and `--gff`
//! sources keep: `AnnotationSourceAdaptor.pm` hands those sources the filter
//! and `AnnotationType/Transcript.pm` evaluates it per transcript.
//!
//! One gene model is written in the three shapes the builder reads: an
//! Ensembl GFF3, the shape of the cache's own transcript set; a GTF; and a
//! GFF3 with no Ensembl identifier prefix, read by Ensembl VEP's generic GFF
//! rules. Each build is annotated under the same filters. The transcript
//! sets the three caches yield are compared with each other and with the
//! sets the model dictates, record by record. The Ensembl GFF3 alone adds a
//! TEC transcript on an `unconfirmed_transcript` line, a type Ensembl VEP's
//! `--gff` source skips, so the sets agree only while a build from the model
//! and the genome reads a line by its type.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// One contig, `21`, of 420 bases. Gene 1 (plus strand, 11 to 160) holds the
/// coding transcript ENST00000000001 and the non-coding ENST00000000002 (one
/// exon, 11 to 160). ENST00000000001 has exons 11 to 70 and 101 to 160 and a
/// CDS from 21 to 131, its stop codon at 129 to 131 and one codon split
/// between the exons, two bases in the first. Gene 2 (minus strand, 201 to
/// 320) holds the lncRNA ENST00000000003 (exons 281 to 320 and 201 to 250).
/// Gene 3 (plus strand, 341 to 400) holds the coding XM_000000004 (one exon,
/// CDS 351 to 389 with its stop codon at 387 to 389), named like a predicted
/// model so that `--exclude_predicted` drops it.
const GENOME: &str = "GCTAAAGACAATTACATAACATGGGAGGAGCCTTTAATGCAGCCCGCAATATGAAGGATTGCAGACAGCGCAGCACGAAACTTGTTGGCCCAGTGTGAATATACCCCAGCTTGCTTCGGAAGGGGACATAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCATGGTACACCTAGAAGAGTCCCTGGAATTTATTGTATAGGCTCTATTGACTACGACGCGCTCATTCCCTT";

/// The model as the Ensembl core API writes a GFF3: `gene:` and
/// `transcript:` identifier prefixes, `biotype` on every gene and transcript
/// line, exon ranks and Ensembl phases, the CDS through the stop codon; plus
/// the TEC ENST00000000005 under gene 1 (one exon, 11 to 160), whose
/// `unconfirmed_transcript` line Ensembl VEP's `--gff` source skips.
const ENSEMBL_GFF3: &str = "##gff-version 3\n\
21\tensembl_havana\tgene\t11\t160\t.\t+\t.\tID=gene:ENSG00000000001;Name=GENE1;biotype=protein_coding;gene_id=ENSG00000000001;version=3\n\
21\tensembl_havana\tmRNA\t11\t160\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000000001;Name=GENE1-201;biotype=protein_coding;tag=gencode_basic,Ensembl_canonical;transcript_id=ENST00000000001;transcript_support_level=1;version=2\n\
21\tensembl_havana\texon\t11\t70\t.\t+\t.\tParent=transcript:ENST00000000001;Name=ENSE00000000001;constitutive=0;ensembl_end_phase=2;ensembl_phase=-1;exon_id=ENSE00000000001;rank=1;version=1\n\
21\tensembl_havana\tCDS\t21\t70\t.\t+\t0\tID=CDS:ENSP00000000001;Parent=transcript:ENST00000000001;protein_id=ENSP00000000001;version=1\n\
21\tensembl_havana\texon\t101\t160\t.\t+\t.\tParent=transcript:ENST00000000001;Name=ENSE00000000002;constitutive=0;ensembl_end_phase=-1;ensembl_phase=2;exon_id=ENSE00000000002;rank=2;version=1\n\
21\tensembl_havana\tCDS\t101\t131\t.\t+\t1\tID=CDS:ENSP00000000001;Parent=transcript:ENST00000000001;protein_id=ENSP00000000001;version=1\n\
21\thavana\ttranscript\t11\t160\t.\t+\t.\tID=transcript:ENST00000000002;Parent=gene:ENSG00000000001;Name=GENE1-202;biotype=retained_intron;transcript_id=ENST00000000002;version=1\n\
21\thavana\texon\t11\t160\t.\t+\t.\tParent=transcript:ENST00000000002;Name=ENSE00000000003;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00000000003;rank=1;version=1\n\
21\thavana\tunconfirmed_transcript\t11\t160\t.\t+\t.\tID=transcript:ENST00000000005;Parent=gene:ENSG00000000001;Name=GENE1-203;biotype=TEC;tag=gencode_basic;transcript_id=ENST00000000005;transcript_support_level=NA;version=1\n\
21\thavana\texon\t11\t160\t.\t+\t.\tParent=transcript:ENST00000000005;Name=ENSE00000000007;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00000000007;rank=1;version=1\n\
21\thavana\tncRNA_gene\t201\t320\t.\t-\t.\tID=gene:ENSG00000000002;Name=LINC1;biotype=lncRNA;gene_id=ENSG00000000002;version=1\n\
21\thavana\tlnc_RNA\t201\t320\t.\t-\t.\tID=transcript:ENST00000000003;Parent=gene:ENSG00000000002;Name=LINC1-201;biotype=lncRNA;tag=Ensembl_canonical;transcript_id=ENST00000000003;version=1\n\
21\thavana\texon\t281\t320\t.\t-\t.\tParent=transcript:ENST00000000003;Name=ENSE00000000004;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00000000004;rank=1;version=1\n\
21\thavana\texon\t201\t250\t.\t-\t.\tParent=transcript:ENST00000000003;Name=ENSE00000000005;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00000000005;rank=2;version=1\n\
21\tensembl_havana\tgene\t341\t400\t.\t+\t.\tID=gene:ENSG00000000003;Name=GENE3;biotype=protein_coding;gene_id=ENSG00000000003;version=1\n\
21\tensembl_havana\tmRNA\t341\t400\t.\t+\t.\tID=transcript:XM_000000004;Parent=gene:ENSG00000000003;Name=GENE3-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=XM_000000004;version=1\n\
21\tensembl_havana\texon\t341\t400\t.\t+\t.\tParent=transcript:XM_000000004;Name=ENSE00000000006;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00000000006;rank=1;version=1\n\
21\tensembl_havana\tCDS\t351\t389\t.\t+\t0\tID=CDS:XP_000000004;Parent=transcript:XM_000000004;protein_id=XP_000000004;version=1\n";

/// The same model as a GTF: the stop codon on its own line, exons with no
/// `rank` attribute, tags one per `tag` attribute, the biotype in
/// `transcript_biotype`.
const GTF: &str = "#!genome-build GRCh38.p14\n\
21\tensembl_havana\tgene\t11\t160\t.\t+\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; gene_name \"GENE1\"; gene_biotype \"protein_coding\";\n\
21\tensembl_havana\ttranscript\t11\t160\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; transcript_version \"2\"; gene_name \"GENE1\"; transcript_biotype \"protein_coding\"; tag \"gencode_basic\"; tag \"Ensembl_canonical\"; transcript_support_level \"1\";\n\
21\tensembl_havana\texon\t11\t70\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"1\"; exon_id \"ENSE00000000001\";\n\
21\tensembl_havana\tCDS\t21\t70\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"1\"; protein_id \"ENSP00000000001\";\n\
21\tensembl_havana\texon\t101\t160\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; exon_id \"ENSE00000000002\";\n\
21\tensembl_havana\tCDS\t101\t128\t.\t+\t1\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; protein_id \"ENSP00000000001\";\n\
21\tensembl_havana\tstop_codon\t129\t131\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; protein_id \"ENSP00000000001\";\n\
21\thavana\ttranscript\t11\t160\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000002\"; transcript_version \"1\"; gene_name \"GENE1\"; transcript_biotype \"retained_intron\";\n\
21\thavana\texon\t11\t160\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000002\"; exon_number \"1\"; exon_id \"ENSE00000000003\";\n\
21\thavana\tgene\t201\t320\t.\t-\t.\tgene_id \"ENSG00000000002\"; gene_version \"1\"; gene_name \"LINC1\"; gene_biotype \"lncRNA\";\n\
21\thavana\ttranscript\t201\t320\t.\t-\t.\tgene_id \"ENSG00000000002\"; transcript_id \"ENST00000000003\"; transcript_version \"1\"; gene_name \"LINC1\"; transcript_biotype \"lncRNA\"; tag \"Ensembl_canonical\";\n\
21\thavana\texon\t281\t320\t.\t-\t.\tgene_id \"ENSG00000000002\"; transcript_id \"ENST00000000003\"; exon_number \"1\"; exon_id \"ENSE00000000004\";\n\
21\thavana\texon\t201\t250\t.\t-\t.\tgene_id \"ENSG00000000002\"; transcript_id \"ENST00000000003\"; exon_number \"2\"; exon_id \"ENSE00000000005\";\n\
21\tensembl_havana\tgene\t341\t400\t.\t+\t.\tgene_id \"ENSG00000000003\"; gene_version \"1\"; gene_name \"GENE3\"; gene_biotype \"protein_coding\";\n\
21\tensembl_havana\ttranscript\t341\t400\t.\t+\t.\tgene_id \"ENSG00000000003\"; transcript_id \"XM_000000004\"; transcript_version \"1\"; gene_name \"GENE3\"; transcript_biotype \"protein_coding\"; tag \"Ensembl_canonical\";\n\
21\tensembl_havana\texon\t341\t400\t.\t+\t.\tgene_id \"ENSG00000000003\"; transcript_id \"XM_000000004\"; exon_number \"1\"; exon_id \"ENSE00000000006\";\n\
21\tensembl_havana\tCDS\t351\t386\t.\t+\t0\tgene_id \"ENSG00000000003\"; transcript_id \"XM_000000004\"; exon_number \"1\"; protein_id \"XP_000000004\";\n\
21\tensembl_havana\tstop_codon\t387\t389\t.\t+\t0\tgene_id \"ENSG00000000003\"; transcript_id \"XM_000000004\"; exon_number \"1\"; protein_id \"XP_000000004\";\n";

/// The same model as a GFF3 without Ensembl's identifier prefixes. Its
/// `mRNA` lines carry no biotype and are `protein_coding` by Ensembl VEP's
/// generic GFF rule; the other transcript lines name theirs. Exons carry no
/// rank, the CDS runs through the stop codon, and the symbol comes from the
/// gene line's `Name`.
const GENERIC_GFF3: &str = "##gff-version 3\n\
21\tensembl_havana\tgene\t11\t160\t.\t+\t.\tID=gene-GENE1;Dbxref=GeneID:1001;Name=GENE1;gbkey=Gene;gene_biotype=protein_coding\n\
21\tensembl_havana\tmRNA\t11\t160\t.\t+\t.\tID=rna-ENST00000000001;Parent=gene-GENE1;gene_id=ENSG00000000001;transcript_id=ENST00000000001;gbkey=mRNA;tag=gencode_basic,Ensembl_canonical;transcript_support_level=1;version=2\n\
21\tensembl_havana\texon\t11\t70\t.\t+\t.\tID=exon-ENST00000000001-1;Parent=rna-ENST00000000001;gbkey=mRNA\n\
21\tensembl_havana\texon\t101\t160\t.\t+\t.\tID=exon-ENST00000000001-2;Parent=rna-ENST00000000001;gbkey=mRNA\n\
21\tensembl_havana\tCDS\t21\t70\t.\t+\t0\tID=cds-ENSP00000000001;Parent=rna-ENST00000000001;gbkey=CDS;protein_id=ENSP00000000001\n\
21\tensembl_havana\tCDS\t101\t131\t.\t+\t1\tID=cds-ENSP00000000001;Parent=rna-ENST00000000001;gbkey=CDS;protein_id=ENSP00000000001\n\
21\thavana\ttranscript\t11\t160\t.\t+\t.\tID=rna-ENST00000000002;Parent=gene-GENE1;gene_id=ENSG00000000001;transcript_id=ENST00000000002;biotype=retained_intron;gbkey=misc_RNA;version=1\n\
21\thavana\texon\t11\t160\t.\t+\t.\tID=exon-ENST00000000002-1;Parent=rna-ENST00000000002;gbkey=misc_RNA\n\
21\thavana\tgene\t201\t320\t.\t-\t.\tID=gene-LINC1;Dbxref=GeneID:1002;Name=LINC1;gbkey=Gene;gene_biotype=lncRNA\n\
21\thavana\tlnc_RNA\t201\t320\t.\t-\t.\tID=rna-ENST00000000003;Parent=gene-LINC1;gene_id=ENSG00000000002;transcript_id=ENST00000000003;biotype=lncRNA;gbkey=ncRNA;tag=Ensembl_canonical;version=1\n\
21\thavana\texon\t281\t320\t.\t-\t.\tID=exon-ENST00000000003-1;Parent=rna-ENST00000000003;gbkey=ncRNA\n\
21\thavana\texon\t201\t250\t.\t-\t.\tID=exon-ENST00000000003-2;Parent=rna-ENST00000000003;gbkey=ncRNA\n\
21\tensembl_havana\tgene\t341\t400\t.\t+\t.\tID=gene-GENE3;Dbxref=GeneID:1003;Name=GENE3;gbkey=Gene;gene_biotype=protein_coding\n\
21\tensembl_havana\tmRNA\t341\t400\t.\t+\t.\tID=rna-XM_000000004;Parent=gene-GENE3;gene_id=ENSG00000000003;transcript_id=XM_000000004;gbkey=mRNA;tag=Ensembl_canonical;version=1\n\
21\tensembl_havana\texon\t341\t400\t.\t+\t.\tID=exon-XM_000000004-1;Parent=rna-XM_000000004;gbkey=mRNA\n\
21\tensembl_havana\tCDS\t351\t389\t.\t+\t0\tID=cds-XP_000000004;Parent=rna-XM_000000004;gbkey=CDS;protein_id=XP_000000004\n";

/// Four substitutions, one in each gene and one between genes 1 and 2; on a
/// 420-base contig every transcript lies within the upstream and downstream
/// distance of each, so each unfiltered record writes a row per transcript.
const VARIANTS: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
21\t30\tv1\tG\tA\t.\t.\t.\n\
21\t180\tv2\tA\tG\t.\t.\t.\n\
21\t220\tv3\tG\tT\t.\t.\t.\n\
21\t360\tv4\tC\tT\t.\t.\t.\n";

const ALL: [&str; 4] = [
    "ENST00000000001",
    "ENST00000000002",
    "ENST00000000003",
    "XM_000000004",
];

/// The clause the builder's two warnings on parents in no line read, or in a
/// line of a type not read, state after the count: what the tree does with
/// the records naming such a parent.
const ORPHAN_RULE: &str = "a record whose parent list opens with one of them is placed as its other parents place it, else skipped with its children";

fn set(ids: &[&str]) -> BTreeSet<String> {
    ids.iter().map(|id| id.to_string()).collect()
}

/// Writes the genome and its index beside each other and returns the FASTA.
fn write_genome(dir: &Path) -> PathBuf {
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">21\n{GENOME}\n")).unwrap();
    std::fs::write(
        dir.join("ref.fa.fai"),
        format!(
            "21\t{}\t4\t{}\t{}\n",
            GENOME.len(),
            GENOME.len(),
            GENOME.len() + 1
        ),
    )
    .unwrap();
    fa
}

/// Runs the builder on one model file under `flag` (`--gff3` or `--gtf`), the
/// genome and the `extra` arguments (none for a build with no protein FASTA)
/// into `dir/<name>`; returns the output directory, the model file written as
/// `dir/<name>.model` and the process output.
fn run_build(
    dir: &Path,
    genome: &Path,
    name: &str,
    flag: &str,
    model: &str,
    extra: &[&str],
) -> (PathBuf, PathBuf, Output) {
    let model_path = dir.join(format!("{name}.model"));
    std::fs::write(&model_path, model).unwrap();
    let out = dir.join(name);
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", flag])
        .arg(&model_path)
        .arg("--genome-fasta")
        .arg(genome)
        .arg("--output-dir")
        .arg(&out)
        .args(extra)
        .output()
        .expect("vep-cache-builder runs");
    (out, model_path, output)
}

/// Builds a cache as [`run_build`] does and returns it with the builder's
/// log; a failed build fails the test.
fn build_with_log(
    dir: &Path,
    genome: &Path,
    name: &str,
    flag: &str,
    model: &str,
    extra: &[&str],
) -> (PathBuf, String) {
    let (out, _, output) = run_build(dir, genome, name, flag, model, extra);
    let log = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{name}: {log}");
    (out, log)
}

fn build(dir: &Path, genome: &Path, name: &str, flag: &str, model: &str) -> PathBuf {
    build_with_log(dir, genome, name, flag, model, &[]).0
}

/// The Ensembl GFF3 with the `version` taken off its three CDS lines, the
/// shape of the GRCh37 release-87 file; every other line keeps its own.
fn ensembl_gff3_without_cds_versions() -> String {
    let unversioned: String = ENSEMBL_GFF3
        .lines()
        .map(|line| match line.contains("\tCDS\t") {
            true => line.replace(";version=1", "") + "\n",
            false => line.to_string() + "\n",
        })
        .collect();
    assert_eq!(
        unversioned.matches("version=").count(),
        ENSEMBL_GFF3.matches("version=").count() - 3
    );
    unversioned
}

/// A protein FASTA of `cache`'s translations: one record per translated
/// transcript, its header opening with the protein identifier, `.version`
/// appended when one is given, then the `transcript:` token; the sequence is
/// the cache's own peptide, so a build from the file differs from `cache` in
/// the version alone.
fn protein_fasta_of(cache: &Path, version: Option<u32>) -> String {
    let mut text = String::new();
    for shard in std::fs::read_dir(cache.join("transcripts/21")).unwrap() {
        let json = std::fs::read_to_string(shard.unwrap().path()).unwrap();
        let records: Vec<serde_json::Value> = serde_json::from_str(&json).unwrap();
        for record in records {
            let Some(protein) = record["translation"]["stable_id"].as_str() else {
                continue;
            };
            let suffix = version.map_or(String::new(), |v| format!(".{v}"));
            text.push_str(&format!(
                ">{protein}{suffix} pep transcript:{}\n{}\n",
                record["stable_id"].as_str().unwrap(),
                record["variation_effect_feature_cache"]["peptide"]
                    .as_str()
                    .unwrap()
            ));
        }
    }
    text
}

/// The translation version of every translated transcript the cache holds,
/// by transcript identifier; `None` for a translation without one.
fn translation_versions(cache: &Path) -> BTreeMap<String, Option<u64>> {
    let mut versions = BTreeMap::new();
    for shard in std::fs::read_dir(cache.join("transcripts/21")).unwrap() {
        let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
        let records: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        for record in records {
            if record["translation"].is_object() {
                versions.insert(
                    record["stable_id"].as_str().unwrap().to_string(),
                    record["translation"]["version"].as_u64(),
                );
            }
        }
    }
    versions
}

/// The stable identifiers of every transcript the cache holds.
fn transcript_ids(cache: &Path) -> BTreeSet<String> {
    let mut ids = BTreeSet::new();
    for shard in std::fs::read_dir(cache.join("transcripts/21")).unwrap() {
        let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
        let records: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        for record in records {
            ids.insert(record["stable_id"].as_str().unwrap().to_string());
        }
    }
    ids
}

/// The transcripts each record's rows name, by record identifier; a record
/// whose rows are all `intergenic_variant` maps to the empty set.
type Features = BTreeMap<String, BTreeSet<String>>;

/// A build a test expects refused: its name, the model flag and text, the
/// extra arguments, the parse the error names, the file it names (none for
/// the model file) and the line.
type RefusedBuild<'a> = (
    &'a str,
    &'a str,
    String,
    &'a [&'a str],
    &'a str,
    Option<&'a Path>,
    usize,
);

/// [`RefusedBuild`] whose error names the model file: no file field.
type RefusedModel<'a> = (&'a str, &'a str, String, &'a [&'a str], &'a str, usize);

/// Annotates the four records against `cache` through the vep-cli library
/// under `flags` and returns the output file.
fn annotate(cache: &Path, genome: &Path, input: &Path, flags: &[&str]) -> PathBuf {
    let out = cache.join(format!(
        "out-{}.txt",
        flags.join("_").replace(['/', ' '], "_")
    ));
    let cache_s = cache.display().to_string();
    let genome_s = genome.display().to_string();
    let input_s = input.display().to_string();
    let out_s = out.display().to_string();
    let mut argv = vec![
        "vep",
        "--offline",
        "--json_cache",
        &cache_s,
        "--fasta",
        &genome_s,
        "--no_stats",
        "--quiet",
        "--force_overwrite",
        "-i",
        &input_s,
        "-o",
        &out_s,
    ];
    argv.extend_from_slice(flags);
    let config = vep_cli::Config::from_argv(argv).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
    out
}

/// The transcripts of the rows `annotate` writes, per record.
fn features(cache: &Path, genome: &Path, input: &Path, flags: &[&str]) -> Features {
    let out = annotate(cache, genome, input, flags);
    let mut features = Features::new();
    for line in std::fs::read_to_string(&out).unwrap().lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut columns = line.split('\t');
        let record = columns.next().unwrap().to_string();
        let feature = columns.nth(3).unwrap();
        let transcripts = features.entry(record).or_default();
        if feature != "-" {
            transcripts.insert(feature.to_string());
        }
    }
    features
}

/// The `Extra` column of the one row `record` writes on `feature` in `out`,
/// as its key-value pairs.
fn extra(out: &Path, record: &str, feature: &str) -> BTreeMap<String, String> {
    let text = std::fs::read_to_string(out).unwrap();
    let rows: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            columns[0] == record && columns[4] == feature
        })
        .collect();
    assert_eq!(rows.len(), 1, "{record} on {feature}: {rows:?}");
    rows[0]
        .rsplit('\t')
        .next()
        .unwrap()
        .split(';')
        .filter_map(|pair| pair.split_once('='))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

/// Every record of the input annotated to one transcript set.
fn every_record(ids: &[&str]) -> Features {
    ["v1", "v2", "v3", "v4"]
        .iter()
        .map(|record| (record.to_string(), set(ids)))
        .collect()
}

struct Caches {
    _dir: tempfile::TempDir,
    genome: PathBuf,
    input: PathBuf,
    /// The Ensembl GFF3 build, then the GTF and the generic GFF3 builds.
    built: [(&'static str, PathBuf); 3],
}

fn caches() -> Caches {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let input = dir.path().join("variants.vcf");
    std::fs::write(&input, VARIANTS).unwrap();
    let built = [
        (
            "ensembl",
            build(dir.path(), &genome, "ensembl", "--gff3", ENSEMBL_GFF3),
        ),
        ("gtf", build(dir.path(), &genome, "gtf", "--gtf", GTF)),
        (
            "generic",
            build(dir.path(), &genome, "generic", "--gff3", GENERIC_GFF3),
        ),
    ];
    Caches {
        _dir: dir,
        genome,
        input,
        built,
    }
}

/// Under every filter, for every record, the GTF build and the generic GFF3
/// build yield the transcripts the Ensembl GFF3 build yields, and those are
/// the transcripts the model dictates. The biotype filter keeps the two
/// coding transcripts. `--exclude_predicted` (`not stable_id match ^X._`,
/// `Config.pm`'s option set) drops the one named like a predicted model.
/// Repeated filters are joined by `and`. A bare field tests its presence
/// (the canonical tag). An underscore-prefixed key reads the gene symbol.
#[test]
fn transcript_filter_applies_to_gff_and_gtf_transcripts() {
    let caches = caches();
    let cases: [(&[&str], Features); 6] = [
        (&[], every_record(&ALL)),
        (
            &["--transcript_filter", "biotype is protein_coding"],
            every_record(&["ENST00000000001", "XM_000000004"]),
        ),
        (
            &["--exclude_predicted"],
            every_record(&["ENST00000000001", "ENST00000000002", "ENST00000000003"]),
        ),
        (
            &[
                "--transcript_filter",
                "biotype is protein_coding",
                "--transcript_filter",
                "not stable_id match ^X._",
            ],
            every_record(&["ENST00000000001"]),
        ),
        (
            &["--transcript_filter", "is_canonical"],
            every_record(&["ENST00000000001", "ENST00000000003", "XM_000000004"]),
        ),
        (
            &["--transcript_filter", "_gene_symbol is GENE1"],
            every_record(&["ENST00000000001", "ENST00000000002"]),
        ),
    ];
    for (flags, expected) in &cases {
        for (name, cache) in &caches.built {
            let yielded = features(cache, &caches.genome, &caches.input, flags);
            assert_eq!(&yielded, expected, "{name} build under {flags:?}");
        }
    }
}

/// The built cache holds the cache's keys, where Ensembl VEP's `--gtf` and
/// `--gff` transcripts hold their own. `source` is read from the file's
/// source column, a key those transcripts never hold. `_gene_hgnc_id`, a key
/// no file carries, is undefined on both. `ccds` reads the cache's `-`
/// placeholder, where those transcripts have no key. `_source_cache`, the
/// file's name on those transcripts, is undefined on the built cache.
#[test]
fn keys_a_gtf_or_gff3_cannot_fill_read_as_the_cache_holds_them() {
    let caches = caches();
    let cases: [(&[&str], Features); 4] = [
        (
            &["--transcript_filter", "source is havana"],
            every_record(&["ENST00000000002", "ENST00000000003"]),
        ),
        (
            &["--transcript_filter", "_gene_hgnc_id exists"],
            every_record(&[]),
        ),
        (&["--transcript_filter", "ccds exists"], every_record(&ALL)),
        (
            &["--transcript_filter", "_source_cache exists"],
            every_record(&[]),
        ),
    ];
    for (flags, expected) in &cases {
        for (name, cache) in &caches.built[1..] {
            let yielded = features(cache, &caches.genome, &caches.input, flags);
            assert_eq!(&yielded, expected, "{name} build under {flags:?}");
        }
    }
}

/// A build from the Ensembl GFF3 and the genome alone reads a line by its
/// type as Ensembl VEP's `--gff` source does (`BaseGXF.pm`
/// `_get_records_by_coords`): the TEC `unconfirmed_transcript` line, a type
/// `GFF.pm`'s list lacks, is skipped and counted, its exon names a parent
/// whose line was skipped, and the cache holds the other four transcripts.
#[test]
fn an_unlisted_ensembl_type_is_absent_from_a_genome_only_build() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let (cache, log) = build_with_log(dir.path(), &genome, "ensembl", "--gff3", ENSEMBL_GFF3, &[]);
    assert!(
        log.contains("Skipped the feature types unconfirmed_transcript (1) of"),
        "{log}"
    );
    assert!(
        log.contains("1 parent identifier names a line of")
            && log.contains(&format!(
                "of a type not read; {ORPHAN_RULE}: transcript:ENST00000000005"
            )),
        "{log}"
    );
    assert!(!log.contains("name no line of"), "{log}");
    assert_eq!(transcript_ids(&cache), set(&ALL));

    // A skipped line named by `Name` alone is remembered under that name,
    // the identifier a read line gets (`GFF.pm` `_record_get_id`).
    let (cache, log) = build_with_log(dir.path(), &genome, "named", "--gff3", NAMED_TEC_GFF3, &[]);
    assert!(
        log.contains(&format!("of a type not read; {ORPHAN_RULE}: LINC1-TEC")),
        "{log}"
    );
    assert!(!log.contains("name no line of"), "{log}");
    assert_eq!(transcript_ids(&cache), set(&["ENST00000000003"]));
}

/// A lncRNA gene, its transcript and exon, then a TEC line of the
/// `unconfirmed_transcript` type, named by `Name` alone, and an exon under it.
const NAMED_TEC_GFF3: &str = "##gff-version 3\n\
21\thavana\tgene\t201\t320\t.\t-\t.\tID=gene:ENSG00000000002;biotype=lncRNA;gene_id=ENSG00000000002;version=1\n\
21\thavana\tlnc_RNA\t201\t320\t.\t-\t.\tID=transcript:ENST00000000003;Parent=gene:ENSG00000000002;biotype=lncRNA;transcript_id=ENST00000000003;version=1\n\
21\thavana\texon\t201\t320\t.\t-\t.\tParent=transcript:ENST00000000003;Name=ENSE00000000009;rank=1\n\
21\thavana\tunconfirmed_transcript\t201\t320\t.\t-\t.\tName=LINC1-TEC;Parent=gene:ENSG00000000002;biotype=TEC\n\
21\thavana\texon\t201\t320\t.\t-\t.\tParent=LINC1-TEC;Name=ENSE00000000010;rank=1\n";

/// A build of the cache's content, one given a protein FASTA, remembers a
/// line it skips by its type under the name its attributes give it, as the
/// genome-only build does: the exon under the TEC line is told its parent's
/// line was of a type not read, never that no line names it, and the cache
/// holds the lncRNA transcript alone. The protein FASTA holds no record, the
/// model no translation.
#[test]
fn a_skipped_line_keeps_its_name_in_a_build_of_the_cache_content() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let proteins = dir.path().join("proteins.fa");
    std::fs::write(&proteins, "").unwrap();
    let proteins = proteins.display().to_string();
    let (cache, log) = build_with_log(
        dir.path(),
        &genome,
        "named-cache",
        "--gff3",
        NAMED_TEC_GFF3,
        &["--fasta", &proteins],
    );
    assert!(
        log.contains("Skipped the feature types unconfirmed_transcript (1) of"),
        "{log}"
    );
    assert!(
        log.contains(&format!("of a type not read; {ORPHAN_RULE}: LINC1-TEC")),
        "{log}"
    );
    assert!(!log.contains("name no line of"), "{log}");
    assert_eq!(transcript_ids(&cache), set(&["ENST00000000003"]));
}

/// The warning on parents in no line states what the tree does with the
/// records naming them (`BaseGXF.pm` `_get_parent_child_structure`): a record
/// whose list opens with such a parent is placed as its other parents place
/// it, else dropped with its children, and a record whose list opens with an
/// empty field is no such record: a record of its own, or the transcript of a
/// gene a later field names.
/// `ENST00000000006`, `Parent=,gene:MISSING1`, is in the cache as a
/// transcript with no gene line; `ENST00000000007`, `Parent=gene:MISSING2`,
/// and its exon are not; `ENST00000000008`, `Parent=gene:MISSING2,gene:ENSG00000000002`,
/// is the second gene's transcript. The warning names the two parents once.
#[test]
fn a_parent_in_no_line_is_warned_as_the_tree_treats_its_records() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let orphans = "##gff-version 3\n\
21\thavana\tncRNA_gene\t201\t320\t.\t-\t.\tID=gene:ENSG00000000002;Name=LINC1;biotype=lncRNA;gene_id=ENSG00000000002;version=1\n\
21\thavana\tlnc_RNA\t201\t320\t.\t-\t.\tID=transcript:ENST00000000003;Parent=gene:ENSG00000000002;biotype=lncRNA;transcript_id=ENST00000000003;version=1\n\
21\thavana\texon\t201\t320\t.\t-\t.\tParent=transcript:ENST00000000003;Name=ENSE00000000004;rank=1\n\
21\thavana\tlnc_RNA\t11\t160\t.\t+\t.\tID=transcript:ENST00000000006;Parent=,gene:MISSING1;biotype=lncRNA;transcript_id=ENST00000000006;version=1\n\
21\thavana\texon\t11\t160\t.\t+\t.\tParent=transcript:ENST00000000006;Name=ENSE00000000011;rank=1\n\
21\thavana\tlnc_RNA\t341\t400\t.\t+\t.\tID=transcript:ENST00000000007;Parent=gene:MISSING2;biotype=lncRNA;transcript_id=ENST00000000007;version=1\n\
21\thavana\texon\t341\t400\t.\t+\t.\tParent=transcript:ENST00000000007;Name=ENSE00000000012;rank=1\n\
21\thavana\tlnc_RNA\t281\t320\t.\t-\t.\tID=transcript:ENST00000000008;Parent=gene:MISSING2,gene:ENSG00000000002;biotype=lncRNA;transcript_id=ENST00000000008;version=1\n\
21\thavana\texon\t281\t320\t.\t-\t.\tParent=transcript:ENST00000000008;Name=ENSE00000000013;rank=1\n";
    let (cache, log) = build_with_log(dir.path(), &genome, "orphans", "--gff3", orphans, &[]);
    assert!(
        log.contains("2 parent identifiers name no line of")
            && log.contains(&format!("; {ORPHAN_RULE}: gene:MISSING1, gene:MISSING2")),
        "{log}"
    );
    assert!(!log.contains("so their children are skipped"), "{log}");
    assert!(!log.contains("of a type not read"), "{log}");
    assert_eq!(
        transcript_ids(&cache),
        set(&["ENST00000000003", "ENST00000000006", "ENST00000000008"])
    );
}

/// The GENCODE basic tag is `tag "basic"` in the GRCh37 release-87 and
/// GENCODE GTFs and `tag "gencode_basic"` in the GRCh38 Ensembl one; a cache
/// built from either spelling keeps the tagged transcript, `ENST00000000001`,
/// under `--gencode_basic` and drops the three untagged ones, as the same flag
/// on Ensembl VEP's cache keeps the transcripts carrying the attribute. The
/// records lie within the distance of every transcript, so no record is
/// intergenic under the flag.
#[test]
fn a_basic_tag_filters_a_gtf_build_under_gencode_basic() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let input = dir.path().join("variants.vcf");
    std::fs::write(&input, VARIANTS).unwrap();
    let grch37 = GTF.replace("tag \"gencode_basic\";", "tag \"basic\";");
    assert!(!grch37.contains("gencode_basic"));
    for (name, model) in [("grch38", GTF), ("grch37", grch37.as_str())] {
        let cache = build(dir.path(), &genome, name, "--gtf", model);
        assert_eq!(
            features(&cache, &genome, &input, &["--gencode_basic"]),
            every_record(&["ENST00000000001"]),
            "{name} spelling under --gencode_basic"
        );
        assert_eq!(
            features(&cache, &genome, &input, &[]),
            every_record(&ALL),
            "{name} spelling without the flag"
        );
    }
}

/// A GTF transcript line's `ccds_id` is the cache's `ccds`, the value as the
/// line carries it (the Ensembl GTF writes the identifier without its
/// version), so a filter on the key reads it: `ENST00000000001`, given
/// `ccds_id "CCDS00001"`, is the one transcript `ccds is CCDS00001` keeps,
/// and the others read the cache's `-` placeholder.
#[test]
fn a_gtf_ccds_id_is_the_cache_ccds() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let input = dir.path().join("variants.vcf");
    std::fs::write(&input, VARIANTS).unwrap();
    let with_ccds = GTF.replace(
        "tag \"gencode_basic\";",
        "tag \"gencode_basic\"; tag \"CCDS\"; ccds_id \"CCDS00001\";",
    );
    let cache = build(dir.path(), &genome, "ccds", "--gtf", &with_ccds);
    assert_eq!(
        features(
            &cache,
            &genome,
            &input,
            &["--transcript_filter", "ccds is CCDS00001"]
        ),
        every_record(&["ENST00000000001"])
    );
    assert_eq!(
        features(
            &cache,
            &genome,
            &input,
            &["--transcript_filter", "ccds is -"]
        ),
        every_record(&["ENST00000000002", "ENST00000000003", "XM_000000004"])
    );
}

/// Ensembl VEP's `--gff` source gives a translation the `version` of its
/// first CDS line, else 1 (`BaseGXF.pm` `_add_translation`), and prints that
/// version after the protein identifier under `--protein_version`
/// (`OutputFactory.pm` lines 1569-1574) and in `HGVSp`
/// (`TranscriptVariationAllele.pm` lines 1759-1760). A cache built from the
/// genome and an Ensembl GFF3 whose CDS lines carry no `version`, the shape
/// of the GRCh37 release-87 file, prints the same `.1` as one built from the
/// file with `version=1`. v1, `21:30 G>A`, is the first base of codon 4 of
/// ENST00000000001 (CDS from 21), `GCC` to `ACC`, p.Ala4Thr; v4, `21:360
/// C>T`, is the first base of codon 4 of XM_000000004 (CDS from 351), `CTA`
/// to `TTA`, both leucine.
#[test]
fn a_cds_line_without_a_version_gives_the_protein_version_1() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let input = dir.path().join("variants.vcf");
    std::fs::write(&input, VARIANTS).unwrap();
    let unversioned = ensembl_gff3_without_cds_versions();
    let flags = ["--hgvs", "--protein", "--protein_version"];
    for (name, model) in [("versioned", ENSEMBL_GFF3), ("unversioned", &unversioned)] {
        let cache = build(dir.path(), &genome, name, "--gff3", model);
        let out = annotate(&cache, &genome, &input, &flags);
        let coding = extra(&out, "v1", "ENST00000000001");
        assert_eq!(coding["ENSP"], "ENSP00000000001.1", "{name}");
        assert_eq!(coding["HGVSp"], "ENSP00000000001.1:p.Ala4Thr", "{name}");
        let predicted = extra(&out, "v4", "XM_000000004");
        assert_eq!(predicted["ENSP"], "XP_000000004.1", "{name}");
        assert_eq!(predicted["HGVSp"], "XP_000000004.1:p.Leu4%3D", "{name}");
    }
}

/// A build from an Ensembl GFF3 or a GTF, a protein FASTA and the genome is
/// the cache's content, whose translation version is the release's. When the
/// CDS lines carry no version, the shape of the GRCh37 release-87 GFF3, the
/// version is the one the FASTA header's protein identifier carries (2 from
/// `ENSP00000000001.2`), printed under `--protein_version` and in `HGVSp`; a
/// header written without one, the shape of a release whose identifiers
/// carry no version, leaves it absent, as Ensembl VEP's cache prints such a
/// translation bare; a GFF3 CDS line's `version` and a GTF CDS line's
/// `protein_version`, the attribute every Ensembl GTF carries, stand over the
/// header's. The peptides are the genome-only build's own, so the five caches
/// differ from it in the version alone.
#[test]
fn a_protein_fasta_header_gives_the_version_a_cds_line_lacks() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let input = dir.path().join("variants.vcf");
    std::fs::write(&input, VARIANTS).unwrap();
    let translated = build(dir.path(), &genome, "translated", "--gff3", ENSEMBL_GFF3);
    let unversioned = ensembl_gff3_without_cds_versions();
    let gtf_with_protein_version: String = GTF
        .lines()
        .map(|line| match line.contains("\tCDS\t") {
            true => line.replacen("\";", "\"; protein_version \"3\";", 1) + "\n",
            false => line.to_string() + "\n",
        })
        .collect();
    assert_eq!(
        gtf_with_protein_version.matches("protein_version").count(),
        3
    );
    let flags = ["--hgvs", "--protein", "--protein_version"];
    let cases: [(&str, &str, &str, Option<u32>, &str); 5] = [
        ("header-versioned", "--gff3", &unversioned, Some(2), ".2"),
        ("header-bare", "--gff3", &unversioned, None, ""),
        ("cds-versioned", "--gff3", ENSEMBL_GFF3, Some(2), ".1"),
        ("gtf-header-versioned", "--gtf", GTF, Some(2), ".2"),
        (
            "gtf-protein-version",
            "--gtf",
            &gtf_with_protein_version,
            Some(2),
            ".3",
        ),
    ];
    for (name, flag, model, header_version, suffix) in cases {
        let proteins = dir.path().join(format!("{name}.pep.fa"));
        std::fs::write(&proteins, protein_fasta_of(&translated, header_version)).unwrap();
        let proteins = proteins.display().to_string();
        let (cache, _) = build_with_log(
            dir.path(),
            &genome,
            name,
            flag,
            model,
            &["--fasta", &proteins],
        );
        let out = annotate(&cache, &genome, &input, &flags);
        let coding = extra(&out, "v1", "ENST00000000001");
        assert_eq!(coding["ENSP"], format!("ENSP00000000001{suffix}"), "{name}");
        assert_eq!(
            coding["HGVSp"],
            format!("ENSP00000000001{suffix}:p.Ala4Thr"),
            "{name}"
        );
        let predicted = extra(&out, "v4", "XM_000000004");
        assert_eq!(predicted["ENSP"], format!("XP_000000004{suffix}"), "{name}");
        assert_eq!(
            predicted["HGVSp"],
            format!("XP_000000004{suffix}:p.Leu4%3D"),
            "{name}"
        );
    }
}

/// A protein FASTA that supplies no record for a translated transcript is
/// named in the build's log, since the build is the cache's content whether
/// or not the file matches the model: a file holding one of the two
/// translations is counted against them, and a file with no record naming a
/// transcript (RefSeq-shaped headers, no `transcript:` token) is called out
/// as supplying no peptide; a file holding both records draws neither line.
#[test]
fn a_protein_fasta_short_of_a_transcript_is_named_in_the_log() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let translated = build(dir.path(), &genome, "translated", "--gff3", ENSEMBL_GFF3);
    let both = protein_fasta_of(&translated, Some(1));
    let one: String = both
        .split('>')
        .filter(|record| record.contains("transcript:ENST00000000001\n"))
        .map(|record| format!(">{record}"))
        .collect();
    assert_eq!(one.matches('>').count(), 1);
    let none = ">NP_000000004.1 some protein\nMKV\n";
    let short = "no record in the protein FASTA";
    let no_record = "holds no record naming a transcript (`transcript:` token) with a sequence";
    let cases: [(&str, &str, Option<&str>, bool); 3] = [
        ("both", both.as_str(), None, false),
        (
            "one",
            one.as_str(),
            Some("1 of 2 translated transcripts has"),
            false,
        ),
        (
            "none",
            none,
            Some("2 of 2 translated transcripts have"),
            true,
        ),
    ];
    for (name, proteins, short_count, no_record_expected) in cases {
        let path = dir.path().join(format!("{name}.pep.fa"));
        std::fs::write(&path, proteins).unwrap();
        let path = path.display().to_string();
        let (_, log) = build_with_log(
            dir.path(),
            &genome,
            name,
            "--gff3",
            ENSEMBL_GFF3,
            &["--fasta", &path],
        );
        match short_count {
            Some(count) => assert!(log.contains(&format!("{count} {short}")), "{name}: {log}"),
            None => assert!(!log.contains(short), "{name}: {log}"),
        }
        assert_eq!(
            log.contains(&format!("{path} {no_record}")),
            no_record_expected,
            "{name}: {log}"
        );
    }
}

/// `text` with `from` replaced by `to` on its line `line`, counted from 1.
fn with_line(text: &str, line: usize, from: &str, to: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_string).collect();
    assert!(lines[line - 1].contains(from), "line {line}");
    lines[line - 1] = lines[line - 1].replacen(from, to, 1);
    lines.join("\n") + "\n"
}

/// A version whose digits do not fit the cache's 32-bit field is refused:
/// the build exits 1 with one `ERROR: [cache-builder]` line naming the file,
/// the line and the limit, and writes nothing, since a version read as none
/// would lose what the file states. A gene, transcript or CDS `version` of a
/// GFF3, a GTF CDS line's `version`, a GTF transcript line's
/// `transcript_version`, on a build of the cache's content (here with a
/// protein FASTA) a GTF CDS line's `protein_version`, and the protein version
/// of a protein FASTA header are each refused at their line, counted with the
/// header. A GTF gene line's `gene_version`, which no reading takes, and a
/// CDS line's `protein_version` on a build from the genome alone, which
/// Ensembl VEP's `--gtf` source never reads, build as without the attribute;
/// 4294967295 is held as the translation's version.
#[test]
fn a_version_above_the_cache_limit_is_refused_and_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let translated = build(dir.path(), &genome, "translated", "--gff3", ENSEMBL_GFF3);
    let proteins = dir.path().join("translated.pep.fa");
    std::fs::write(&proteins, protein_fasta_of(&translated, Some(2))).unwrap();
    let above_proteins = dir.path().join("above.pep.fa");
    std::fs::write(
        &above_proteins,
        protein_fasta_of(&translated, Some(2)).replacen(".2 pep", ".4294967296 pep", 1),
    )
    .unwrap();
    let proteins = proteins.display().to_string();
    let above_proteins_arg = above_proteins.display().to_string();
    let above = "4294967296";
    let limit = "version 4294967296 is above 4294967295, the largest version the cache holds";
    let cds_version = format!("protein_id \"ENSP00000000001\"; version \"{above}\";");
    let protein_version = format!("protein_id \"ENSP00000000001\"; protein_version \"{above}\";");
    let refused: [RefusedBuild<'_>; 7] = [
        (
            "gff3-gene",
            "--gff3",
            with_line(ENSEMBL_GFF3, 2, "version=3", &format!("version={above}")),
            &[],
            "GFF3",
            None,
            2,
        ),
        (
            "gff3-transcript",
            "--gff3",
            with_line(ENSEMBL_GFF3, 3, "version=2", &format!("version={above}")),
            &[],
            "GFF3",
            None,
            3,
        ),
        (
            "gff3-cds",
            "--gff3",
            with_line(ENSEMBL_GFF3, 5, "version=1", &format!("version={above}")),
            &[],
            "GFF3",
            None,
            5,
        ),
        (
            "gtf-cds",
            "--gtf",
            with_line(GTF, 5, "protein_id \"ENSP00000000001\";", &cds_version),
            &[],
            "GTF",
            None,
            5,
        ),
        (
            "gtf-transcript",
            "--gtf",
            with_line(
                GTF,
                3,
                "transcript_version \"2\"",
                &format!("transcript_version \"{above}\""),
            ),
            &[],
            "GTF",
            None,
            3,
        ),
        (
            "gtf-protein-version",
            "--gtf",
            with_line(GTF, 5, "protein_id \"ENSP00000000001\";", &protein_version),
            &["--fasta", &proteins],
            "GTF",
            None,
            5,
        ),
        (
            "fasta-header",
            "--gtf",
            GTF.to_string(),
            &["--fasta", &above_proteins_arg],
            "protein FASTA",
            Some(&above_proteins),
            1,
        ),
    ];
    for (name, flag, model, extra, parse, named, line) in &refused {
        let (out, model_path, output) = run_build(dir.path(), &genome, name, flag, model, extra);
        let err = String::from_utf8_lossy(&output.stderr).into_owned();
        assert_eq!(output.status.code(), Some(1), "{name}: {err}");
        let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
        let file = named.unwrap_or(&model_path).display();
        assert_eq!(
            errors,
            [format!(
                "ERROR: [cache-builder] Failed to parse {parse}: {file} line {line}: {limit}"
            )],
            "{name}: {err}"
        );
        assert!(
            !out.exists(),
            "{name}: {:?}",
            std::fs::read_dir(&out).map(|d| d.count())
        );
    }

    let unread_gene_version = with_line(
        GTF,
        2,
        "gene_version \"3\"",
        &format!("gene_version \"{above}\""),
    );
    let gtf = build(
        dir.path(),
        &genome,
        "gtf-gene-version",
        "--gtf",
        &unread_gene_version,
    );
    assert_eq!(transcript_ids(&gtf), set(&ALL));
    let unread_protein_version =
        with_line(GTF, 5, "protein_id \"ENSP00000000001\";", &protein_version);
    let gtf = build(
        dir.path(),
        &genome,
        "gtf-unread-protein-version",
        "--gtf",
        &unread_protein_version,
    );
    assert_eq!(translation_versions(&gtf)["ENST00000000001"], Some(1));
    let held = with_line(ENSEMBL_GFF3, 5, "version=1", "version=4294967295");
    let gff3 = build(dir.path(), &genome, "gff3-held", "--gff3", &held);
    assert_eq!(
        translation_versions(&gff3)["ENST00000000001"],
        Some(4294967295)
    );
}

/// A version whose text is not a whole number is refused as a version above
/// the cache's limit is: the build exits 1 with one `ERROR: [cache-builder]`
/// line naming the file, the line and the text, and writes nothing, where
/// Ensembl VEP's `--gff` and `--gtf` sources print the text after the
/// identifier. A gene, transcript or CDS `version` of a GFF3, a GTF CDS
/// line's `version`, a GTF transcript line's `transcript_version` and, on a
/// build of the cache's content, a GTF CDS line's `protein_version` are each
/// refused at their line. A protein FASTA header whose identifier ends in a
/// `.` followed by anything but digits states no version, the shape other
/// genome projects' protein identifiers take, so the build goes on with the
/// translation unversioned. A whole number written with leading zeros is that
/// number, so a CDS `version=01` is the translation version 1.
#[test]
fn a_version_that_is_not_a_whole_number_is_refused_and_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let translated = build(dir.path(), &genome, "translated", "--gff3", ENSEMBL_GFF3);
    let proteins = dir.path().join("translated.pep.fa");
    std::fs::write(&proteins, protein_fasta_of(&translated, Some(2))).unwrap();
    let cds_proteins = dir.path().join("cds.pep.fa");
    std::fs::write(
        &cds_proteins,
        protein_fasta_of(&translated, None).replace(" pep", ".cds1 pep"),
    )
    .unwrap();
    let proteins = proteins.display().to_string();
    let cds_proteins = cds_proteins.display().to_string();
    let refusal = "version \"v2\" is not a whole number";
    let cds_version = "protein_id \"ENSP00000000001\"; version \"v2\";";
    let protein_version = "protein_id \"ENSP00000000001\"; protein_version \"v2\";";
    let refused: [RefusedModel<'_>; 6] = [
        (
            "gff3-gene-text",
            "--gff3",
            with_line(ENSEMBL_GFF3, 2, "version=3", "version=v2"),
            &[],
            "GFF3",
            2,
        ),
        (
            "gff3-transcript-text",
            "--gff3",
            with_line(ENSEMBL_GFF3, 3, "version=2", "version=v2"),
            &[],
            "GFF3",
            3,
        ),
        (
            "gff3-cds-text",
            "--gff3",
            with_line(ENSEMBL_GFF3, 5, "version=1", "version=v2"),
            &[],
            "GFF3",
            5,
        ),
        (
            "gtf-cds-text",
            "--gtf",
            with_line(GTF, 5, "protein_id \"ENSP00000000001\";", cds_version),
            &[],
            "GTF",
            5,
        ),
        (
            "gtf-transcript-text",
            "--gtf",
            with_line(
                GTF,
                3,
                "transcript_version \"2\"",
                "transcript_version \"v2\"",
            ),
            &[],
            "GTF",
            3,
        ),
        (
            "gtf-protein-version-text",
            "--gtf",
            with_line(GTF, 5, "protein_id \"ENSP00000000001\";", protein_version),
            &["--fasta", &proteins],
            "GTF",
            5,
        ),
    ];
    for (name, flag, model, extra, parse, line) in &refused {
        let (out, model_path, output) = run_build(dir.path(), &genome, name, flag, model, extra);
        let err = String::from_utf8_lossy(&output.stderr).into_owned();
        assert_eq!(output.status.code(), Some(1), "{name}: {err}");
        let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
        assert_eq!(
            errors,
            [format!(
                "ERROR: [cache-builder] Failed to parse {parse}: {} line {line}: {refusal}",
                model_path.display()
            )],
            "{name}: {err}"
        );
        assert!(!out.exists(), "{name}");
    }

    let (gtf, _) = build_with_log(
        dir.path(),
        &genome,
        "gtf-cds-headers",
        "--gtf",
        GTF,
        &["--fasta", &cds_proteins],
    );
    assert_eq!(translation_versions(&gtf)["ENST00000000001"], None);
    let zeros = with_line(ENSEMBL_GFF3, 5, "version=1", "version=01");
    let gff3 = build(dir.path(), &genome, "gff3-zeros", "--gff3", &zeros);
    assert_eq!(translation_versions(&gff3)["ENST00000000001"], Some(1));
}
