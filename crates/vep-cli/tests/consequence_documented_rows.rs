// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes on the inputs of eleven documented
//! differences whose rows the consequence, HGVS and structural-variant code of
//! `vep-effects` decides. Each test annotates the input of its entry's record
//! under `scripts/divergences/reproductions`, over the cache and with the
//! flags that record names, and holds every row of the output; the test of
//! the stop codon holds the rows that carry its entry's term, and the test of
//! the intron span reads the committed corpus that holds its record's model
//! and holds that model's row.
//! Beside those eleven, six tests hold the `code` row of a record to the
//! functions that implement its entry, one holds the rows the record of the
//! intron span compares to that corpus, and one holds the sentence of its
//! page that names the compared cells to the record.

use std::collections::{BTreeMap, HashMap};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The lines of an entry's reproduction record that carry one of `ids` in
/// their ID column, in the record's order: lines of the `text` of the
/// `[input]` table of `scripts/divergences/reproductions/<entry>.toml`, a VCF
/// without its header. The records are read here so that each input line is
/// spelt in one file.
fn record_lines(entry: &str, ids: &[&str]) -> String {
    let path = committed(&format!("scripts/divergences/reproductions/{entry}.toml"));
    let record = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = record
        .lines()
        .skip_while(|line| *line != "[input]")
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .filter(|line| line.split('\t').nth(2).is_some_and(|id| ids.contains(&id)))
        .collect();
    assert_eq!(lines.len(), ids.len(), "{}: {ids:?}", path.display());
    lines.iter().map(|line| format!("{line}\n")).collect()
}

/// One run of `vep` on the input of an entry's record.
struct Run<'a> {
    entry: &'a str,
    /// The records of the entry's input the run annotates, by their IDs; none
    /// for a region input, whose lines carry no ID and are all annotated.
    records: &'a [&'a str],
    json_cache: PathBuf,
    species: &'a str,
    assembly: &'a str,
    fasta: Option<PathBuf>,
    /// The record's flags, given after the arguments every run takes.
    flags: &'a [&'a str],
}

/// The lines of the default output of `run` on `input`, a file of the input
/// format `format`, from its column line on: the `##` lines dropped. The run
/// ends with exit code 0.
fn output_rows(dir: &Path, run: &Run, input: &Path, format: &str) -> Vec<String> {
    let out = dir.join(format!("{}.txt", run.entry));
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(&run.json_cache)
        .args(["--species", run.species, "--assembly", run.assembly])
        .args(["--format", format])
        .args(["--force_overwrite", "--quiet", "--no_stats"]);
    if let Some(fasta) = &run.fasta {
        cmd.arg("--fasta").arg(fasta);
    }
    let output = cmd.args(run.flags).stdin(Stdio::null()).output().unwrap();
    assert!(
        output.status.success(),
        "{} {:?}: {}",
        run.entry,
        run.flags,
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with("##"))
        .map(str::to_string)
        .collect()
}

/// The lines of the default output of `run`, whose record names VCF lines, from
/// its column line on.
fn rows(dir: &Path, run: &Run) -> Vec<String> {
    let input = dir.join(format!("{}.vcf", run.entry));
    let records = record_lines(run.entry, run.records);
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    output_rows(dir, run, &input, "vcf")
}

/// The region lines of an entry's reproduction record: the `text` of its
/// `[input]` table, whose `format` is `region`, whole.
fn record_region_lines(entry: &str) -> String {
    let path = committed(&format!("scripts/divergences/reproductions/{entry}.toml"));
    let record = std::fs::read_to_string(&path).unwrap();
    let mut table = record.lines().skip_while(|line| *line != "[input]").skip(1);
    assert_eq!(
        table.next(),
        Some("format = \"region\""),
        "{}",
        path.display()
    );
    assert_eq!(table.next(), Some("text = \"\"\""), "{}", path.display());
    let lines: Vec<&str> = table.take_while(|line| *line != "\"\"\"").collect();
    assert!(!lines.is_empty(), "{}", path.display());
    lines.iter().map(|line| format!("{line}\n")).collect()
}

/// The lines of the default output of `run`, whose record names region lines,
/// from its column line on.
fn region_rows(dir: &Path, run: &Run) -> Vec<String> {
    let input = dir.join(format!("{}.region.txt", run.entry));
    std::fs::write(&input, record_region_lines(run.entry)).unwrap();
    output_rows(dir, run, &input, "region")
}

/// Under `--shift_3prime 1` the mature miRNA test reads the span the shift
/// gives. The record deletes four bases whose written span begins on the last
/// base of the mature sequence of ENST00000584048; the shift carries the
/// deletion three bases on, past that sequence, and the row is
/// `non_coding_transcript_exon_variant`. Ensembl VEP 116.3 reads the written
/// span there (`VariationEffect.pm:545-549`) and writes `mature_miRNA_variant`.
#[test]
fn a_deletion_shifted_off_a_mature_mirna_is_no_mature_mirna_variant() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "mature_mirna_unshifted_span",
        records: &["rs1324558451"],
        json_cache: committed("tests/golden/116/GRCh38-release116/json_cache"),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(committed(
            "tests/golden/116/GRCh38-plugin-spliceai-fixture/reference.fa.gz",
        )),
        flags: &["--hgvs", "--shift_3prime", "1", "--shift_length"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000154736\tENST00000284987\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000284987.6:c.1237+751_1237+754del;SHIFT_LENGTH=0",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000266133\tENST00000584048\tTranscript\tnon_coding_transcript_exon_variant\t28-31\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000584048.1:n.28_31del;HGVS_OFFSET=3;SHIFT_LENGTH=3",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000154736\tENST00000652031\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;FLAGS=cds_start_NF;HGVSc=ENST00000652031.1:c.419+751_419+754del;SHIFT_LENGTH=0",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000154736\tENST00000970346\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000970346.2:c.1237+751_1237+754del;SHIFT_LENGTH=0",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000154736\tENST00001108182\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00001108182.1:c.1237+751_1237+754del;SHIFT_LENGTH=0",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000154736\tENST00001108183\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00001108183.1:c.1237+751_1237+754del;SHIFT_LENGTH=0",
            "rs1324558451\t21:26953985-26953988\t-\tENSG00000154736\tENST00001131497\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00001131497.1:c.1237+751_1237+754del;SHIFT_LENGTH=0",
        ]
    );
}

/// A record keeps the rows of the transcripts it overlaps whatever the other
/// records of its batch hold. The duplication's two rows and the breakend's two
/// stand beside a breakend whose mate position, with the 5,000-base flank
/// added, reaches 2^31. Ensembl VEP 116.3 chooses the cache regions and the
/// bounds of an input buffer from every record's coordinates and every mate
/// position (`AnnotationSource.pm:190`), and writes every record of that buffer
/// as `intergenic_variant`.
#[test]
fn a_record_beside_a_breakend_whose_mate_reaches_the_32_bit_bound_keeps_its_transcript_rows() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "bnd_mate_32bit_buffer_intergenic",
        records: &["control_dup", "end2_2147478648"],
        json_cache: committed("tests/golden/116/GRCh38-release116/json_cache"),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: None,
        flags: &[],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "control_dup\t21:30425651-30425899\tduplication\tENSG00000186971\tENST00000334068\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=4321;STRAND=1",
            "control_dup\t21:30425651-30425899\tduplication\tENSG00000240432\tENST00000390690\tTranscript\tfeature_elongation,coding_sequence_variant\t28-276\t14-262\t5-88\t-\t-\t-\tIMPACT=HIGH;STRAND=-1;OverlapBP=249;OverlapPC=45.52",
            "end2_2147478648\t21:30425651\tN.\tENSG00000186971\tENST00000334068\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=4569;STRAND=1",
            "end2_2147478648\t21:30425651\tN.\tENSG00000240432\tENST00000390690\tTranscript\tfeature_truncation,coding_sequence_variant\t276\t262\t88\t-\t-\t-\tIMPACT=HIGH;STRAND=-1",
        ]
    );
}

/// Each allele of a record is read at its own positions. The record holds two
/// insertion alleles at the stop codon of ENSMUST00000405524, each shifted one
/// base under `--hgvs`; the rows of the later allele on the two transcripts
/// whose stop codon it stands in hold their own protein position and codons
/// under `inframe_insertion,stop_retained_variant`. Ensembl VEP 116.3 reads a
/// later allele at the protein positions the allele before it left on the
/// object the alleles share (`TranscriptVariationAllele.pm:1724-1733`,
/// `BaseTranscriptVariation.pm:365-380`) and writes that allele's row on
/// ENSMUST00000405524 as `frameshift_variant` with no protein cell.
#[test]
fn a_later_insertion_allele_is_read_at_its_own_positions() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "hgvs_shift_shared_across_alleles",
        records: &["rs265482355"],
        json_cache: committed("tests/golden/116/GRCm39-mouse/json_cache"),
        species: "mus_musculus",
        assembly: "GRCm39",
        fasta: Some(committed("tests/golden/116/GRCm39-mouse/reference.fa.gz")),
        flags: &["--hgvs"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000081714\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000081714.6:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000081714\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000081714.6:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000167381\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000167381.2:n.261-769_261-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000167381\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000167381.2:n.261-769_261-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000172244\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000172244.9:c.194-4073_194-4072insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000172244\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000172244.9:c.194-4073_194-4072insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000405524\tTranscript\tinframe_insertion,stop_retained_variant,NMD_transcript_variant\t597-598\t206-207\t69\t*/*KG\ttaa/taAAAGGGa\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000405524.1:c.207_*1insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000405524\tTranscript\tinframe_insertion,stop_retained_variant,NMD_transcript_variant\t597-598\t206-207\t69\t*/*KX\ttaa/taAAAGa\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000405524.1:c.207_*1insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000405525\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405525.1:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000405525\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405525.1:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000405526\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405526.1:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000405526\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405526.1:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000405527\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405527.1:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000405527\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405527.1:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000405528\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405528.1:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000405528\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405528.1:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000405529\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405529.1:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000405529\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000405529.1:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000577435\tTranscript\tinframe_insertion,stop_retained_variant,NMD_transcript_variant\t240-241\t203-204\t68\t*/*KG\ttaa/taAAAGGGa\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000577435.1:c.204_*1insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000577435\tTranscript\tinframe_insertion,stop_retained_variant,NMD_transcript_variant\t240-241\t203-204\t68\t*/*KX\ttaa/taAAAGa\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000577435.1:c.204_*1insAAGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAGGG\tENSMUSG00000025176\tENSMUST00000577436\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000577436.1:c.194-769_194-768insAAGGGA;HGVS_OFFSET=1",
            "rs265482355\t19:42047557-42047558\tAAAG\tENSMUSG00000025176\tENSMUST00000577436\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000577436.1:c.194-769_194-768insAAGA;HGVS_OFFSET=1",
        ]
    );
}

/// Under `--shift_3prime 1` the intron test reads the span the shift gives. Each
/// of the two records deletes one base of a run before an acceptor site: the
/// shift carries the first onto the first base of the site and the second past
/// the site into the exon, and neither row holds `intron_variant`. Ensembl VEP
/// 116.3 reads the written span for that term
/// (`BaseTranscriptVariationAllele.pm:123-124`) and adds it to both rows.
#[test]
fn a_deletion_shifted_out_of_an_intron_is_no_intron_variant() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "shift_intron_term_unshifted_span",
        records: &["rs1315375439", "rs544638589"],
        json_cache: committed("tests/golden/116/GRCh38-release116/json_cache"),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(committed("tests/golden/116/GRCh38-merged/reference.fa.gz")),
        flags: &["--hgvs", "--shift_3prime", "1", "--shift_length"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "rs1315375439\t21:22883199\t-\tENSG00000225906\tENST00000438328\tTranscript\tsplice_acceptor_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000438328.1:n.308-2del;HGVS_OFFSET=2;SHIFT_LENGTH=2",
            "rs1315375439\t21:22883199\t-\tENSG00000299538\tENST00000764465\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=220;STRAND=-1;SHIFT_LENGTH=0",
            "rs544638589\t21:13065551\t-\tENSG00000224309\tENST00000435744\tTranscript\tsplice_region_variant,non_coding_transcript_exon_variant\t1288\t-\t-\t-\t-\t-\tIMPACT=LOW;STRAND=1;HGVSc=ENST00000435744.5:n.1288del;HGVS_OFFSET=8;SHIFT_LENGTH=8",
        ]
    );
}

/// A FASTA of one sequence written into `dir` for the record `id` of
/// `entry`, a deletion: `N` up to the record, the record's reference allele at
/// its place, after it the first `repeat` deleted bases once more, and one
/// `N`.
fn fasta_of_the_record_and_its_repeat(dir: &Path, entry: &str, id: &str, repeat: usize) -> PathBuf {
    let line = record_lines(entry, &[id]);
    let cells: Vec<&str> = line.trim_end().split('\t').collect();
    let (name, position, reference) = (cells[0], cells[1].parse::<usize>().unwrap(), cells[3]);
    let sequence = format!(
        "{}{reference}{}N",
        "N".repeat(position - 1),
        &reference[1..=repeat]
    );
    let fasta = dir.join("reference.fa");
    std::fs::write(&fasta, format!(">{name}\n{sequence}\n")).unwrap();
    let (bases, first_base) = (sequence.len(), name.len() + 2);
    std::fs::write(
        dir.join("reference.fa.fai"),
        format!("{name}\t{bases}\t{first_base}\t{bases}\t{}\n", bases + 1),
    )
    .unwrap();
    fasta
}

/// Under `--shift_3prime 1` the stop codon test reads the span the shift gives.
/// The record deletes 3,804 bases that end six bases before ENST00001097586;
/// the 513 bases after the deleted span repeat its first 513, so the shift
/// carries the deletion over that transcript's stop codon, and six rows, that
/// transcript's among them, hold `stop_retained_variant`. Ensembl VEP 116.3
/// reads the written span for the overlap (`VariationEffect.pm:1403`) and
/// writes `coding_sequence_variant` on ENST00001097586. The tree holds no
/// GRCh38 reference of chromosome 16, so the run reads one written from the
/// record itself.
#[test]
fn a_deletion_shifted_onto_the_stop_codon_retains_the_stop() {
    let dir = tempfile::tempdir().unwrap();
    let entry = "stop_codon_unshifted_span";
    let run = Run {
        entry,
        records: &["38636"],
        json_cache: committed("tests/golden/116/GRCh38/json_cache"),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(fasta_of_the_record_and_its_repeat(
            dir.path(),
            entry,
            "38636",
            513,
        )),
        flags: &["--hgvs", "--shift_3prime", "1", "--shift_length"],
    };
    let rows = rows(dir.path(), &run);
    assert_eq!(rows.len(), 70, "the column line and 69 rows");
    let retained: Vec<&str> = rows
        .iter()
        .map(String::as_str)
        .filter(|row| row.contains("stop_retained_variant"))
        .collect();
    assert_eq!(
        retained,
        [
            "38636\t16:172871-176674\t-\tENSG00000206172\tENST00001097574\tTranscript\tsplice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00001097574.1:c.-45_*98+54del;HGVS_OFFSET=513;SHIFT_LENGTH=513",
            "38636\t16:172871-176674\t-\tENSG00000206172\tENST00001097576\tTranscript\tsplice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00001097576.1:c.-42_*152+54del;HGVS_OFFSET=513;SHIFT_LENGTH=513",
            "38636\t16:172871-176674\t-\tENSG00000206172\tENST00001097577\tTranscript\tsplice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00001097577.1:c.-38_*98+54del;HGVS_OFFSET=513;SHIFT_LENGTH=513",
            "38636\t16:172871-176674\t-\tENSG00000206172\tENST00001097580\tTranscript\tsplice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00001097580.1:c.-42_*152+54del;HGVS_OFFSET=513;SHIFT_LENGTH=513",
            "38636\t16:172871-176674\t-\tENSG00000206172\tENST00001097583\tTranscript\tsplice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00001097583.1:c.-38_*130+54del;HGVS_OFFSET=513;SHIFT_LENGTH=513",
            "38636\t16:172871-176674\t-\tENSG00000206172\tENST00001097586\tTranscript\tsplice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant\t?-342\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00001097586.1:c.-37_*206del;HGVS_OFFSET=513;SHIFT_LENGTH=513",
        ]
    );
}

/// The items the `code` row of an entry's record names, in the record's order:
/// the strings of the `code` list of the `[implementation]` table of
/// `scripts/divergences/reproductions/<entry>.toml`.
fn code_row(entry: &str) -> Vec<String> {
    let path = committed(&format!("scripts/divergences/reproductions/{entry}.toml"));
    let record = std::fs::read_to_string(&path).unwrap();
    let list = record
        .split_once("\n[implementation]\n")
        .and_then(|(_, table)| table.split_once("code = ["))
        .and_then(|(_, rest)| rest.split_once(']'))
        .unwrap_or_else(|| panic!("{}: no code row", path.display()))
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The record of the stop codon overlap names the three functions that read
/// the stop codon at the span the shift gives. `consequences_of_the_record`
/// hands every predicate the shifted allele. `PerlCodingEval::stop_codon_cil_span`
/// is the span the overlap is read at where the stop codon edit is tested: a
/// deletion with a cDNA coordinate at both ends and a CDS start.
/// `PerlCodingEval::given_stop_codon` decides the term where that edit is not
/// tested, as on the record's own rows, whose shifted span begins before the
/// transcript: it asks that span to hold the three bases of the stop codon and
/// reads the sequence the record gives.
#[test]
fn the_stop_codon_overlap_record_names_the_functions_that_read_the_shifted_span() {
    assert_eq!(
        code_row("stop_codon_unshifted_span"),
        [
            "crates/vep-effects/src/coding.rs::PerlCodingEval::given_stop_codon",
            "crates/vep-effects/src/coding.rs::PerlCodingEval::stop_codon_cil_span",
            "crates/vep-effects/src/consequences.rs::consequences_of_the_record",
        ]
    );
}

/// Under `--shift_3prime 1` the start codon test reads the cDNA span the shift
/// gives. The record deletes the last base of an intron and the first three
/// bases of the exon after it; the shift carries the deletion one base into the
/// exon, far from the start codon, and both rows are
/// `frameshift_variant,splice_region_variant` with a frameshift HGVSp. Ensembl
/// VEP 116.3 adds the shift to the cDNA ends of the written span, one of them
/// undefined (`VariationEffect.pm:977-981`), and on ENST00000360731, which has
/// no 5' UTR, adds `start_lost` and writes a `?` in the HGVSp.
#[test]
fn a_deletion_shifted_into_an_exon_is_tested_against_the_start_codon_where_it_lands() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "shift_start_codon_test",
        records: &["del_exon_boundary_coding"],
        json_cache: committed("tests/golden/116/GRCh37/json_cache"),
        species: "homo_sapiens",
        assembly: "GRCh37",
        fasta: Some(committed("tests/golden/116/GRCh37-hgvs/reference.fa.gz")),
        flags: &["--shift_3prime", "1", "--hgvs"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "del_exon_boundary_coding\t21:36088618-36088621\t-\tENSG00000159212\tENST00000349499\tTranscript\tframeshift_variant,splice_region_variant\t1900-1903\t1900-1903\t634-635\tIV/X\tATTGtg/tg\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000349499.2:c.1900_1903del;HGVSp=ENSP00000290332.4:p.Ile634TrpfsTer14;HGVS_OFFSET=1",
            "del_exon_boundary_coding\t21:36088618-36088621\t-\tENSG00000159212\tENST00000360731\tTranscript\tframeshift_variant,splice_region_variant\t1954-1957\t1954-1957\t652-653\tIV/X\tATTGtg/tg\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000360731.3:c.1954_1957del;HGVSp=ENSP00000353959.3:p.Ile652TrpfsTer14;HGVS_OFFSET=1",
        ]
    );
}

/// The 3' shift of `--hgvs` walks the reference as far as the repeat runs, also
/// within 1,000 bases of a sequence's first base. The record deletes three
/// bases 291 bases into the mitochondrial genome, on a reverse-strand
/// transcript; the repeat admits a shift of four bases, which the HGVSc and
/// HGVS_OFFSET of the row carry. Ensembl VEP 116.3 walks the first 1,000 bases
/// of its clipped slice in place of the bases beside the allele
/// (`TranscriptVariationAllele.pm:427-437`) and shifts two.
#[test]
fn a_deletion_near_a_sequence_edge_is_shifted_to_the_end_of_its_repeat() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "hgvs_shift_sequence_edge_window",
        records: &["mt_del3_290_TCTT_T"],
        json_cache: committed("tests/golden/116/TAIR10-arabidopsis/json_cache"),
        species: "arabidopsis_thaliana",
        assembly: "TAIR10",
        fasta: Some(committed(
            "tests/golden/116/TAIR10-arabidopsis/reference.fa.gz",
        )),
        flags: &["--hgvs"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "mt_del3_290_TCTT_T\tMt:291-293\t-\tATMG00010\tATMG00010.1\tTranscript\tinframe_deletion\t442-444\t442-444\t148\tK/-\tAAG/-\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ATMG00010.1:c.446_448del;HGVSp=ATMG00010.1:p.Lys149del;HGVS_OFFSET=-4",
        ]
    );
}

/// The HGVSg of a row describes the row's own allele. One record names the
/// reference base among its alternate alleles and the other repeats an
/// alternate allele; under `--hgvsg --allele_number` every row of both carries
/// the substitution of its own allele. Ensembl VEP 116.3 takes the allele it
/// describes from the allele string by the row's allele number
/// (`OutputFactory.pm:1386-1393`), so it writes a duplication for the rows of
/// the first record and, for the rows of the last allele of the second, the
/// notation of the repeated allele. Each row prints the number of its allele's
/// place in the record, where Ensembl VEP 116.3 counts without the repeated and
/// the reference-equal entry and prints 1 on the first record's rows and 2 on
/// the rows of the second record's last allele.
#[test]
fn a_repeated_or_reference_equal_alternate_allele_gets_its_own_hgvsg() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "hgvsg_repeated_or_reference_equal_alt",
        records: &["refeq", "dupalt"],
        json_cache: committed("tests/golden/116/GRCh37-hgvs/json_cache"),
        species: "homo_sapiens",
        assembly: "GRCh37",
        fasta: Some(committed("tests/golden/116/GRCh37-hgvs/reference.fa.gz")),
        flags: &["--hgvsg", "--allele_number"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "refeq\t21:34775850\tC\tENSG00000159128\tENST00000290219\tTranscript\tstart_lost\t649\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=2;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "refeq\t21:34775850\tC\tENSG00000159128\tENST00000381995\tTranscript\tstart_lost\t113\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=2;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "refeq\t21:34775850\tC\tENSG00000159128\tENST00000405436\tTranscript\t5_prime_UTR_variant\t79\t-\t-\t-\t-\t-\tALLELE_NUM=2;IMPACT=MODIFIER;STRAND=1;HGVSg=21:g.34775850A>C",
            "refeq\t21:34775850\tC\tENSG00000159128\tENST00000439213\tTranscript\tstart_lost,NMD_transcript_variant\t118\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=2;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "refeq\t21:34775850\tC\tENSG00000159128\tENST00000545369\tTranscript\tstart_lost,NMD_transcript_variant\t104\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=2;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "dupalt\t21:34775850\tC\tENSG00000159128\tENST00000290219\tTranscript\tstart_lost\t649\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=1;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "dupalt\t21:34775850\tG\tENSG00000159128\tENST00000290219\tTranscript\tstart_lost\t649\t1\t1\tM/V\tAtg/Gtg\t-\tALLELE_NUM=3;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>G",
            "dupalt\t21:34775850\tC\tENSG00000159128\tENST00000381995\tTranscript\tstart_lost\t113\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=1;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "dupalt\t21:34775850\tG\tENSG00000159128\tENST00000381995\tTranscript\tstart_lost\t113\t1\t1\tM/V\tAtg/Gtg\t-\tALLELE_NUM=3;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>G",
            "dupalt\t21:34775850\tC\tENSG00000159128\tENST00000405436\tTranscript\t5_prime_UTR_variant\t79\t-\t-\t-\t-\t-\tALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1;HGVSg=21:g.34775850A>C",
            "dupalt\t21:34775850\tG\tENSG00000159128\tENST00000405436\tTranscript\t5_prime_UTR_variant\t79\t-\t-\t-\t-\t-\tALLELE_NUM=3;IMPACT=MODIFIER;STRAND=1;HGVSg=21:g.34775850A>G",
            "dupalt\t21:34775850\tC\tENSG00000159128\tENST00000439213\tTranscript\tstart_lost,NMD_transcript_variant\t118\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=1;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "dupalt\t21:34775850\tG\tENSG00000159128\tENST00000439213\tTranscript\tstart_lost,NMD_transcript_variant\t118\t1\t1\tM/V\tAtg/Gtg\t-\tALLELE_NUM=3;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>G",
            "dupalt\t21:34775850\tC\tENSG00000159128\tENST00000545369\tTranscript\tstart_lost,NMD_transcript_variant\t104\t1\t1\tM/L\tAtg/Ctg\t-\tALLELE_NUM=1;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>C",
            "dupalt\t21:34775850\tG\tENSG00000159128\tENST00000545369\tTranscript\tstart_lost,NMD_transcript_variant\t104\t1\t1\tM/V\tAtg/Gtg\t-\tALLELE_NUM=3;IMPACT=HIGH;STRAND=1;HGVSg=21:g.34775850A>G",
        ]
    );
}

/// The reference codons of a deletion across bases an edited model lacks are
/// read from the model's own coding sequence. NM_030582.4 and NM_130444.3 lack
/// nine genomic bases; the record's deletion of seventeen bases holds eight
/// coding positions of each, and the reference codons of both rows end in the
/// model's base. Ensembl VEP 116.3 writes the genome's allele over the coding
/// span before it reads the codons (`TranscriptVariationAllele.pm:872-874`), so
/// its reference codons end in the genome's base.
#[test]
fn a_deletion_across_bases_an_edited_model_lacks_reads_the_models_codons() {
    let dir = tempfile::tempdir().unwrap();
    let run = Run {
        entry: "refseq_edited_model_deleted_span_codons",
        records: &["rs1179352929"],
        json_cache: committed(
            "tests/golden/116/GRCh38-refseq/cache_root/homo_sapiens_refseq/116_GRCh38",
        ),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(committed("tests/golden/116/GRCh38-refseq/reference.fa.gz")),
        flags: &["--refseq"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "rs1179352929\t21:45504519-45504535\t-\t80781\tNM_001379500.1\tTranscript\tframeshift_variant\t2897-2913\t2831-2847\t944-949\tPGPPGP/X\tcCAGGCCCCCCAGGCCCA/c\t-\tIMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA",
            "rs1179352929\t21:45504519-45504535\t-\t80781\tNM_030582.4\tTranscript\tframeshift_variant\t3380-3387\t3371-3378\t1124-1126\tPGP/X\tcCAGGCCCA/c\t-\tIMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK",
            "rs1179352929\t21:45504519-45504535\t-\t80781\tNM_130444.3\tTranscript\tframeshift_variant\t4085-4092\t4076-4083\t1359-1361\tPGP/X\tcCAGGCCCA/c\t-\tIMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK",
            "rs1179352929\t21:45504519-45504535\t-\t6573\tXM_047440957.1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA",
            "rs1179352929\t21:45504519-45504535\t-\t6573\tXM_047440963.1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA",
            "rs1179352929\t21:45504519-45504535\t-\t6573\tXM_047440964.1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA",
        ]
    );
}

/// The record of the deleted-span codons names the three functions that read
/// the reference codons and their residues from the model, without the splice
/// of `TranscriptVariationAllele.pm:872-874`. `codon` of `display.rs` reads
/// the reference window from the model's coding sequence and so gives the
/// Codons cell its reference half. `peptide` of `display.rs` translates that
/// window and, on a model whose translation carries no SeqEdit, writes the
/// cached protein's residues over the translation where the two differ. The
/// reference half of Amino_acids on the record's two in-frame deletions is
/// therefore the model's residues whichever window is read, and takes Ensembl
/// VEP's value only with both functions changed. `PerlCodingEval::codon` reads
/// the same window for the predicates, so a term that follows the reference
/// codons follows the model's bases.
#[test]
fn the_deleted_span_record_names_the_functions_that_read_the_reference_codons() {
    assert_eq!(
        code_row("refseq_edited_model_deleted_span_codons"),
        [
            "crates/vep-effects/src/coding.rs::PerlCodingEval::codon",
            "crates/vep-effects/src/display.rs::codon",
            "crates/vep-effects/src/display.rs::peptide",
        ]
    );
}

/// The record of the inserted-span codons names the two functions of
/// `display.rs` that write the reference halves of the Codons and Amino_acids
/// cells it compares. `codon` reads the reference window from the model's
/// coding sequence, without the write of `TranscriptVariationAllele.pm:872-873`,
/// and so gives the Codons cell its reference half. `peptide` translates that
/// window and, on a model whose translation carries no SeqEdit, writes the
/// cached protein's residues over the translation where the two differ. The
/// reference half of Amino_acids is therefore the model's residues whichever
/// window is read, and takes Ensembl VEP's value only with both functions
/// changed.
#[test]
fn the_inserted_span_record_names_the_functions_that_write_the_reference_halves() {
    assert_eq!(
        code_row("refseq_edited_model_inserted_span_codons"),
        [
            "crates/vep-effects/src/display.rs::codon",
            "crates/vep-effects/src/display.rs::peptide",
        ]
    );
}

/// The record of the own-bases insertion term names the three methods that
/// give its `synonymous_variant`. `PerlCodingEval::codon` reads the reference
/// window from the model's coding sequence, and `PerlCodingEval::peptide`, on
/// a model whose translation carries no SeqEdit, writes the cached protein's
/// residues over the window's translation where the two differ: each keeps
/// the reference peptide equal to the alternate's with the other changed, and
/// the term leaves the row only with both changed.
/// `PerlCodingEval::synonymous_variant` accepts the row on that equality.
#[test]
fn the_own_bases_record_names_the_methods_that_read_the_reference_peptide() {
    assert_eq!(
        code_row("refseq_edited_model_own_bases_insertion_term"),
        [
            "crates/vep-effects/src/coding.rs::PerlCodingEval::codon",
            "crates/vep-effects/src/coding.rs::PerlCodingEval::peptide",
            "crates/vep-effects/src/coding.rs::PerlCodingEval::synonymous_variant",
        ]
    );
}

/// The record of the shorter insertion term names the three methods that give
/// its `inframe_deletion`. `PerlCodingEval::stands_on_more_model_bases` is the
/// length test, an allele longer than its genomic span and shorter than the
/// CDS span it stands on, and `PerlCodingEval::terms` asks `inframe_deletion`
/// of an allele that passes it, where Ensembl VEP asks the predicate of a
/// deletion alone (`Utils/Constants.pm:866-869`): the two hold the term of the
/// record that keeps both bases beside the point. `PerlCodingEval::codon`
/// reads the reference window from the model's coding sequence, without the
/// write of `TranscriptVariationAllele.pm:872-874`, which cuts the span out
/// under a plain insertion's empty reference allele. The predicate holds on
/// the model's window and fails on the shortened one, so that method holds
/// the term of the plain insertion.
#[test]
fn the_shorter_insertion_record_names_the_methods_that_give_its_inframe_deletion() {
    assert_eq!(
        code_row("refseq_edited_model_shorter_insertion_term"),
        [
            "crates/vep-effects/src/coding.rs::PerlCodingEval::codon",
            "crates/vep-effects/src/coding.rs::PerlCodingEval::stands_on_more_model_bases",
            "crates/vep-effects/src/coding.rs::PerlCodingEval::terms",
        ]
    );
}

/// The reference codons of a deletion that holds a whole intron of an edited
/// model are read from the model's own coding sequence. The record's region
/// line deletes 693 bases of NM_000555.3: the last 12 bases of exon 1, the 333
/// bases of intron 1 and the first 348 bases of exon 2, which is 120 codons in
/// frame. The row holds `inframe_deletion`, the model's CDS 208 to 570 and its
/// residues 70 to 190 before the slashes of Codons and Amino_acids, and a
/// plain deletion as its HGVSp. Ensembl VEP 116.3 writes the genome's allele,
/// intron bases included, over the coding span before it reads the reference
/// codons (`TranscriptVariationAllele.pm:872-874`), translates two stops
/// inside the intron and writes `stop_lost,protein_altering_variant` with a
/// stop-lost extension on the HGVSp. The run reads the cache and the reference
/// of the corpus `GRCh38-refseq-intron-span`, which holds the model.
#[test]
fn a_deletion_across_a_whole_intron_of_an_edited_model_reads_the_models_codons() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = "tests/golden/116/GRCh38-refseq-intron-span";
    let run = Run {
        entry: "refseq_edited_model_intron_span_coding_reference",
        records: &[],
        json_cache: committed(&format!("{corpus}/json_cache")),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(committed(&format!("{corpus}/reference.fa.gz"))),
        flags: &["--refseq", "--hgvs"],
    };
    let rows = region_rows(dir.path(), &run);
    assert_eq!(rows.len(), 12, "the column line and 11 rows");
    let of_the_model: Vec<&str> = rows
        .iter()
        .map(String::as_str)
        .filter(|row| row.split('\t').nth(4) == Some("NM_000555.3"))
        .collect();
    assert_eq!(
        of_the_model,
        [
            "X:111410073-111410765:1/-\tX:111410073-111410765\t-\t1641\tNM_000555.3\tTranscript\tsplice_acceptor_variant,splice_donor_variant,inframe_deletion,splice_donor_5th_base_variant,intron_variant\t382-741\t210-569\t70-190\tRSHCQSLRFHQNMELDFGHFDERDKTSRNMRGSRMNGLPSPTHSAHCSFYRTRTLQALSNEKKAKKVRFYRNGDRYFKGIVYAVSSDRFRSFDALLADLTRSLSDNINLPQGVRYIYTIDG/R\tagGAGCCACTGTCAGTCTCTGAGGTTCCACCAAAATATGGAACTTGATTTTGGACACTTTGACGAAAGAGATAAGACATCCAGGAACATGCGAGGCTCCCGGATGAATGGGTTGCCTAGCCCCACTCACAGCGCCCACTGTAGCTTCTACCGAACCAGAACCTTGCAGGCACTGAGTAATGAGAAGAAAGCCAAGAAGGTACGTTTCTACCGCAATGGGGACCGCTACTTCAAGGGGATTGTGTACGCTGTGTCCTCTGACCGTTTTCGCAGCTTTGACGCCTTGCTGGCTGACCTGACGCGATCTCTGTCTGACAACATCAACCTGCCTCAGGGAGTGCGTTACATTTACACCATTGATGGA/aga\t-\tIMPACT=HIGH;STRAND=-1;GIVEN_REF=CCATCAATGGTGTAAATGTAACGCACTCCCTGAGGCAGGTTGATGTTGTCAGACAGAGATCGCGTCAGGTCAGCCAGCAAGGCGTCAAAGCTGCGAAAACGGTCAGAGGACACAGCGTACACAATCCCCTTGAAGTAGCGGTCCCCATTGCGGTAGAAACGTACCTTCTTGGCTTTCTTCTCATTACTCAGTGCCTGCAAGGTTCTGGTTCGGTAGAAGCTACAGTGGGCGCTGTGAGTGGGGCTAGGCAACCCATTCATCCGGGAGCCTCGCATGTTCCTGGATGTCTTATCTCTTTCGTCAAAGTGTCCAAAATCAAGTTCCATATTTTGGTGGAACCTCAGAGACCTGAGCGTGGGAGAAAGGGATGGGGGTGAAGAGAGGCAAAAACAAAAAGGGACAAGGAGAAGGAAAAAAGAAGGGATTTTAAATATTAGCCAGATCCAGATCTGCAATCTGCTGCTTGTGAAAAGAACTGGTTGAAAAAAGCCTGTGACCCATCTGCTGCTTGCCCCTGACCTGCAGGAATATTGATCTTAATAGAGAAGAAGGAGGGGCTGGGCGGGGGTGCTGGTGGTTGTGGGGGTTGGGAGTAAGAGATAGAGAGGGAGGCACTCTTGGCACTGTTCCAGACAGAGGAGAAGGGGAGATTTTGAAATAGCTGAAAATCCTGAGACTTACTGACAGTGGCTC;USED_REF=CCATCAATGGTGTAAATGTAACGCACTCCCTGAGGCAGGTTGATGTTGTCAGACAGAGATCGCGTCAGGTCAGCCAGCAAGGCGTCAAAGCTGCGAAAACGGTCAGAGGACACAGCGTACACAATCCCCTTGAAGTAGCGGTCCCCATTGCGGTAGAAACGTACCTTCTTGGCTTTCTTCTCATTACTCAGTGCCTGCAAGGTTCTGGTTCGGTAGAAGCTACAGTGGGCGCTGTGAGTGGGGCTAGGCAACCCATTCATCCGGGAGCCTCGCATGTTCCTGGATGTCTTATCTCTTTCGTCAAAGTGTCCAAAATCAAGTTCCATATTTTGGTGGAACCTCAGAGACCTGAGCGTGGGAGAAAGGGATGGGGGTGAAGAGAGGCAAAAACAAAAAGGGACAAGGAGAAGGAAAAAAGAAGGGATTTTAAATATTAGCCAGATCCAGATCTGCAATCTGCTGCTTGTGAAAAGAACTGGTTGAAAAAAGCCTGTGACCCATCTGCTGCTTGCCCCTGACCTGCAGGAATATTGATCTTAATAGAGAAGAAGGAGGGGCTGGGCGGGGGTGCTGGTGGTTGTGGGGGTTGGGAGTAAGAGATAGAGAGGGAGGCACTCTTGGCACTGTTCCAGACAGAGGAGAAGGGGAGATTTTGAAATAGCTGAAAATCCTGAGACTTACTGACAGTGGCTC;BAM_EDIT=OK;HGVSc=NM_000555.3:c.210_569del;HGVSp=NP_000546.2:p.His72_Ser191del",
        ]
    );
}

/// The record of the intron-span coding reference names the four functions
/// that read the reference codons and residues from the model's coding
/// sequence, without the write of `TranscriptVariationAllele.pm:872-874`.
/// `PerlCodingEval::codon` reads the window the predicates compare, and
/// `PerlCodingEval::peptide`, on a model whose translation carries no SeqEdit,
/// writes the cached protein's residues over the window's translation where
/// the two differ. With the first changed alone the row loses
/// `inframe_deletion`; it takes `stop_lost` and the stop-lost extension of its
/// HGVSp only with both changed. `codon` of `display.rs` writes the reference
/// half of Codons, and `peptide` of `display.rs` keeps the reference half of
/// Amino_acids on the cached residues: that half takes the translation of the
/// genome's allele only with both changed.
#[test]
fn the_intron_span_record_names_the_functions_that_read_the_models_codons_and_residues() {
    assert_eq!(
        code_row("refseq_edited_model_intron_span_coding_reference"),
        [
            "crates/vep-effects/src/coding.rs::PerlCodingEval::codon",
            "crates/vep-effects/src/coding.rs::PerlCodingEval::peptide",
            "crates/vep-effects/src/display.rs::codon",
            "crates/vep-effects/src/display.rs::peptide",
        ]
    );
}

/// The strings of the list `key` of the table `table` of an entry's record,
/// in the record's order: a list of
/// `scripts/divergences/reproductions/<entry>.toml` none of whose strings
/// holds a bracket or a double quote.
fn record_list(entry: &str, table: &str, key: &str) -> Vec<String> {
    let path = committed(&format!("scripts/divergences/reproductions/{entry}.toml"));
    let record = std::fs::read_to_string(&path).unwrap();
    let list = record
        .split_once(&format!("\n[{table}]"))
        .and_then(|(_, rest)| rest.split_once(&format!("\n{key} = [")))
        .and_then(|(_, rest)| rest.split_once(']'))
        .unwrap_or_else(|| panic!("{}: no {key} list under [{table}]", path.display()))
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The cells of `line`, a row of a default output whose column line is
/// `header`, under `columns`, joined as the `rows` comparison of a record
/// joins them: the cell of a column, or the value of the `Extra` item of that
/// name, and `-` where the row has neither.
fn compared_cells(header: &str, line: &str, columns: &[String]) -> String {
    let names: Vec<&str> = header.trim_start_matches('#').split('\t').collect();
    let cells: Vec<&str> = line.split('\t').collect();
    let cell = |name: &str| {
        names
            .iter()
            .position(|column| *column == name)
            .and_then(|at| cells.get(at).copied())
    };
    let extra = cell("Extra").unwrap_or_default();
    columns
        .iter()
        .map(|column| {
            cell(column)
                .or_else(|| {
                    extra
                        .split(';')
                        .map(|item| item.split_once('=').unwrap_or((item, "")))
                        .find(|(key, _)| key == column)
                        .map(|(_, value)| value)
                })
                .unwrap_or("-")
        })
        .collect::<Vec<_>>()
        .join(" | ")
}

/// The lines of a stored default output, a gzip file, from its column line
/// on: the `##` lines dropped.
fn stored_rows(path: &Path) -> Vec<String> {
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .filter(|line| !line.starts_with("##"))
        .map(str::to_string)
        .collect()
}

/// The record of the intron-span coding reference compares, on the row of
/// `NM_000555.3`, the two cells its page gives as the entry's own difference,
/// `Amino_acids` and `Codons`, beside the terms, the three positions and the
/// protein notation. Under those columns its `vep_rs` row is the row the
/// binary writes on the record's region line over the corpus
/// `GRCh38-refseq-intron-span`, and its `ensembl_vep` row is the row of that
/// corpus's stored Ensembl VEP 116.3 output, which annotates the same line.
#[test]
fn the_intron_span_record_compares_the_codons_and_residues_of_both_engines() {
    let entry = "refseq_edited_model_intron_span_coding_reference";
    let columns = record_list(entry, "compare", "columns");
    assert_eq!(
        columns,
        [
            "Feature",
            "Consequence",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "HGVSp",
        ]
    );
    let dir = tempfile::tempdir().unwrap();
    let corpus = "tests/golden/116/GRCh38-refseq-intron-span";
    assert_eq!(
        std::fs::read_to_string(committed(&format!("{corpus}/variants.region.txt"))).unwrap(),
        record_region_lines(entry),
        "the corpus annotates the record's line"
    );
    let run = Run {
        entry,
        records: &[],
        json_cache: committed(&format!("{corpus}/json_cache")),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(committed(&format!("{corpus}/reference.fa.gz"))),
        flags: &["--refseq", "--hgvs"],
    };
    let written = region_rows(dir.path(), &run);
    let stored = stored_rows(&committed(&format!(
        "{corpus}/expected_region/default.txt.gz"
    )));
    for (engine, rows) in [("vep_rs", &written), ("ensembl_vep", &stored)] {
        let of_the_model: Vec<String> = rows[1..]
            .iter()
            .filter(|row| row.split('\t').nth(4) == Some("NM_000555.3"))
            .map(|row| compared_cells(&rows[0], row, &columns))
            .collect();
        assert_eq!(of_the_model.len(), 1, "{engine}");
        assert_eq!(
            record_list(entry, "expected", engine),
            of_the_model,
            "{engine}"
        );
    }
}

/// The page of the intron-span coding reference names, in the sentence on
/// what its reproduction compares, the cells the record compares on the row
/// of `NM_000555.3` after the row's `Feature`, in the record's order, the
/// three position cells as one item, and says that the record holds the two
/// cells of codons and residues whole.
#[test]
fn the_intron_span_page_names_the_cells_its_reproduction_compares() {
    let page = std::fs::read_to_string(committed(
        "docs/divergences/refseq-edited-model-intron-span-coding-reference.md",
    ))
    .unwrap();
    let cells = page
        .lines()
        .find_map(|line| {
            line.strip_prefix("On the row of `NM_000555.3` the reproduction compares ")
        })
        .and_then(|rest| rest.split_once(", and it expects from each engine"))
        .expect("the sentence on what the reproduction compares")
        .0;
    let named: Vec<&str> = cells
        .split(", ")
        .flat_map(|items| items.split(" and "))
        .flat_map(|item| match item {
            "the three position cells" => {
                vec!["cDNA_position", "CDS_position", "Protein_position"]
            }
            cell => vec![cell.trim_matches('`')],
        })
        .collect();
    let columns = record_list(
        "refseq_edited_model_intron_span_coding_reference",
        "compare",
        "columns",
    );
    assert_eq!(columns[0], "Feature");
    assert_eq!(named, columns[1..]);
    assert!(page.contains("It holds each engine's `Amino_acids` and `Codons` cells whole."));
    assert!(!page.contains("It does not compare the reference halves"));
}

/// The GFF3 text the first `prepare` step of an entry's record writes: the
/// body of its here-document.
fn record_gff3(entry: &str) -> String {
    let path = committed(&format!("scripts/divergences/reproductions/{entry}.toml"));
    let text: String = std::fs::read_to_string(&path)
        .unwrap()
        .lines()
        .skip_while(|line| !line.ends_with("<<'GFF'"))
        .skip(1)
        .take_while(|line| !line.starts_with("GFF'''"))
        .map(|line| format!("{line}\n"))
        .collect();
    assert!(text.starts_with("##gff-version 3\n"), "{}", path.display());
    text
}

/// The human GRCh38 cache of release 116 that the cache builder writes into
/// `dir` from a GFF3 text and a genome FASTA alone, by the library calls its
/// binary makes for such a build.
fn cache_built_from_gff3(dir: &Path, gff3: &str, genome: &Path) -> PathBuf {
    use vep_cache_builder::{builder, gff3 as model, info};

    let file = dir.join("model.gff3");
    std::fs::write(&file, gff3).unwrap();
    let genome = vep_fasta::IndexedFasta::from_path(genome).unwrap();
    let names = model::SequenceNames::of_fasta(&genome, Some("GRCh38"));
    let genes = model::parse_gff3(&file, model::Reading::GffSource, &names).unwrap();
    let codon_tables = builder::codon_tables(&genes, None, &HashMap::new(), &names).unwrap();
    let cache = dir.join("json_cache");
    let sources = builder::BuildSources {
        protein_sequences: None,
        prediction_data: None,
        genome_fasta: Some(&genome),
        gtf_tags: None,
        core_dumps: None,
    };
    let stats = builder::build_cache(&genes, &file, &cache, &sources, &codon_tables).unwrap();
    info::CacheInfoJson::new(
        "homo_sapiens",
        "GRCh38",
        116,
        &stats,
        None,
        BTreeMap::new(),
        None,
    )
    .write(&cache)
    .unwrap();
    vep_core::cache_manifest::CacheManifest::write(&cache).unwrap();
    cache
}

/// A coding transcript whose translation has no identifier gets no HGVSp and no
/// ENSP. The cache is built from the GFF3 lines of the record, in which the
/// first CDS line of ENST00001113783 in transcript order carries no
/// `protein_id`; the row of that transcript holds neither cell and the row of
/// the other transcript holds both. Ensembl VEP 116.3 reads the identifier from
/// that line alone (`BaseGXF.pm:723`) and names the reference of the HGVSp by
/// the empty identifier and the version
/// (`TranscriptVariationAllele.pm:1748`, `:1759`).
#[test]
fn a_coding_transcript_without_a_protein_identifier_gets_no_hgvsp() {
    let dir = tempfile::tempdir().unwrap();
    let entry = "gff_cds_without_protein_id_hgvsp";
    let genome = committed("tests/golden/116/GRCh38-gtf-builder/reference.fa.gz");
    let run = Run {
        entry,
        records: &["tmprss3_codon6", "pwp2_codon3"],
        json_cache: cache_built_from_gff3(dir.path(), &record_gff3(entry), &genome),
        species: "homo_sapiens",
        assembly: "GRCh38",
        fasta: Some(genome),
        flags: &["--hgvs", "--protein", "--protein_version", "--symbol"],
    };
    assert_eq!(
        rows(dir.path(), &run),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "tmprss3_codon6\t21:42395401\tA\tENSG00000160183\tENST00001113783\tTranscript\tmissense_variant\t249\t17\t6\tP/L\tcCg/cTg\t-\tIMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;HGVSc=ENST00001113783.1:c.17C>T",
            "pwp2_codon3\t21:44107424\tC\tENSG00000241945\tENST00000456705\tTranscript\tmissense_variant\t52\t8\t3\tF/S\ttTc/tCc\t-\tIMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser",
        ]
    );
}
