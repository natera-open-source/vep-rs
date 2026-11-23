// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The names a build gives its sequences: the names its genome FASTA holds.
//! Each build reads the synthetic one-gene model of
//! `tests/fixtures/mito_codon_table` with its sequence renamed in the GFF3,
//! the FASTA or both.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The fixture's GFF3 and genome FASTA under `dir`, the GFF3's sequence
/// written as `annotated` and the FASTA's record named `held`.
fn renamed_inputs(dir: &Path, annotated: &str, held: &str) -> (PathBuf, PathBuf) {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mito_codon_table");
    let model = std::fs::read_to_string(fixture.join("mito.gff3"))
        .unwrap()
        .replace("Mito", annotated);
    let genome = std::fs::read_to_string(fixture.join("mito.fa")).unwrap();
    let (header, bases) = genome.split_once('\n').unwrap();
    assert_eq!(header, ">Mito");
    let index = std::fs::read_to_string(fixture.join("mito.fa.fai")).unwrap();
    let rest = index.strip_prefix("Mito\t90\t6\t").unwrap();
    let gff3 = dir.join("model.gff3");
    let fasta = dir.join("genome.fa");
    std::fs::write(&gff3, model).unwrap();
    std::fs::write(&fasta, format!(">{held}\n{bases}")).unwrap();
    std::fs::write(
        dir.join("genome.fa.fai"),
        format!("{held}\t90\t{}\t{rest}", held.len() + 2),
    )
    .unwrap();
    (gff3, fasta)
}

/// Builds a cache under `dir` from the renamed inputs and `flags`.
fn build(dir: &Path, annotated: &str, held: &str, flags: &[&str]) -> Output {
    let (gff3, fasta) = renamed_inputs(dir, annotated, held);
    Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "made_up", "--gff3"])
        .arg(&gff3)
        .arg("--genome-fasta")
        .arg(&fasta)
        .arg("--output-dir")
        .arg(dir.join("cache"))
        .args(flags)
        .output()
        .expect("vep-cache-builder runs")
}

/// The cache's chromosome directories, and the one transcript's
/// `codon_table` and `translateable_seq`.
fn built(dir: &Path) -> (Vec<String>, serde_json::Value, serde_json::Value) {
    let mut chromosomes: Vec<String> = std::fs::read_dir(dir.join("cache/transcripts"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    chromosomes.sort();
    let shard = dir
        .join("cache/transcripts")
        .join(&chromosomes[0])
        .join("1-1000000.json");
    let records: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(shard).unwrap()).unwrap();
    let cache = &records[0]["variation_effect_feature_cache"];
    (
        chromosomes,
        cache["codon_table"].clone(),
        cache["translateable_seq"].clone(),
    )
}

/// The bases of the fixture's coding sequence, as the FASTA holds them.
const CODING_BASES: &str = "ATGCTTTGAATAAAATAA";

/// A sequence whose own name begins with `chr` keeps it: a genome that calls
/// its sequence `chromosome_1` in the GFF3 and in the FASTA is built under
/// that name, with its coding sequence read from the FASTA, and `info.json`
/// lists the name.
#[test]
fn a_sequence_named_chromosome_1_keeps_its_name() {
    let dir = tempfile::tempdir().unwrap();
    let output = build(dir.path(), "chromosome_1", "chromosome_1", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let (chromosomes, codon_table, coding) = built(dir.path());
    assert_eq!(chromosomes, ["chromosome_1"]);
    assert_eq!(codon_table, 1);
    assert_eq!(coding, CODING_BASES);
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.path().join("cache/info.json")).unwrap())
            .unwrap();
    assert_eq!(info["contig_order"], serde_json::json!(["chromosome_1"]));
}

/// A mitochondrion the GFF3 writes as `chrM` is built as `MT` where the
/// genome FASTA holds one under `chrM` or `MT`, with its coding sequence
/// read from that record, and where the FASTA names no organelle genome at
/// all (`21` alone). It is kept as `chrM` where the FASTA names its
/// mitochondrion `Mito`: no sequence of that genome is called `MT`. Either
/// way the region's table is the one `--codon-table chrM=<n>` names, and
/// without the flag the build stops, naming the region as the cache would
/// hold it.
#[test]
fn a_mitochondrion_written_chrm_is_kept_only_where_the_fasta_names_its_own() {
    for (held, named) in [("chrM", "MT"), ("MT", "MT"), ("21", "MT"), ("Mito", "chrM")] {
        let dir = tempfile::tempdir().unwrap();
        let output = build(dir.path(), "chrM", held, &["--codon-table", "chrM=3"]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{held}: {stderr}");
        let (chromosomes, codon_table, coding) = built(dir.path());
        assert_eq!(chromosomes, [named], "{held}");
        assert_eq!(codon_table, 3, "{held}");
        if held == "Mito" || held == "21" {
            assert!(coding.is_null(), "{held}: {coding}");
        } else {
            assert_eq!(coding, CODING_BASES, "{held}");
        }

        let dir = tempfile::tempdir().unwrap();
        let output = build(dir.path(), "chrM", held, &[]);
        assert!(!output.status.success(), "{held}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        let refusal = format!(
            "ERROR: [cache-builder] codon table of sequence region {named} unknown: pass \
             --core-dumps or --codon-table {named}=<n>"
        );
        assert!(
            stderr.lines().any(|line| line == refusal),
            "{held}: {stderr}"
        );
    }
}

/// A region the genome FASTA alone shortens takes the table the flag names
/// for it as the GFF3 writes it: a sequence written `chrx` under a FASTA
/// that holds `x` is built as `x`, under the table of `--codon-table chrx=4`
/// and of `--codon-table x=4` alike.
#[test]
fn a_region_the_fasta_shortens_takes_the_flags_table_under_either_name() {
    for flag in ["chrx=4", "x=4"] {
        let dir = tempfile::tempdir().unwrap();
        let output = build(dir.path(), "chrx", "x", &["--codon-table", flag]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{flag}: {stderr}");
        assert!(
            !stderr.contains("a sequence region no transcript sits on"),
            "{flag}: {stderr}"
        );
        let (chromosomes, codon_table, coding) = built(dir.path());
        assert_eq!(chromosomes, ["x"], "{flag}");
        assert_eq!(codon_table, 4, "{flag}");
        assert_eq!(coding, CODING_BASES, "{flag}");
    }
}

/// The fixture's GFF3 under `dir` with its sequence written as `annotated`,
/// and a genome FASTA of one record named `held` of `length` bases: the
/// fixture's 90 bases, then a repeating four, on lines of 60.
fn sized_inputs(dir: &Path, annotated: &str, held: &str, length: usize) -> (PathBuf, PathBuf) {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mito_codon_table");
    let model = std::fs::read_to_string(fixture.join("mito.gff3"))
        .unwrap()
        .replace("Mito", annotated);
    let genome = std::fs::read_to_string(fixture.join("mito.fa")).unwrap();
    let mut bases: Vec<u8> = genome
        .lines()
        .skip(1)
        .flat_map(|line| line.bytes())
        .collect();
    assert_eq!(bases.len(), 90);
    bases.extend(b"ACGT".iter().copied().cycle().take(length - 90));
    let header = format!(">{held}\n");
    let mut text = header.clone();
    for line in bases.chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    let gff3 = dir.join("model.gff3");
    let fasta = dir.join("genome.fa");
    std::fs::write(&gff3, model).unwrap();
    std::fs::write(&fasta, text).unwrap();
    std::fs::write(
        dir.join("genome.fa.fai"),
        format!("{held}\t{length}\t{}\t60\t61\n", header.len()),
    )
    .unwrap();
    (gff3, fasta)
}

/// On GRCh37 `MT` names the assembly's mitochondrial sequence, of 16,569
/// bases, and `chrM` and `M` also name the one of the UCSC hg19 naming, of
/// 16,571. A build for GRCh37 names a mitochondrion the GFF3 writes `chrM`
/// `MT` only where the genome FASTA holds a record of the assembly's length
/// for it, under any of the four names, or gives no record for `MT` and
/// names no organelle genome (the test after this one); over a record of
/// another length the region keeps its written name, its coding sequence
/// read from the FASTA's record of that naming, takes the table
/// `--codon-table` names for it, and the build says once that it kept the
/// name. A build for another assembly names it `MT` over either record.
#[test]
fn on_grch37_a_mitochondrion_written_chrm_is_mt_only_over_a_record_of_the_assemblys_length() {
    for (assembly, held, length, named) in [
        ("GRCh37", "chrM", 16_571, "chrM"),
        ("GRCh37", "M", 16_571, "chrM"),
        ("GRCh37", "chrM", 16_569, "MT"),
        ("GRCh37", "M", 16_569, "MT"),
        ("GRCh37", "MT", 16_569, "MT"),
        ("GRCh37", "chrMT", 16_569, "MT"),
        ("GRCh38", "chrM", 16_571, "MT"),
        ("made_up", "chrM", 16_571, "MT"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let (gff3, fasta) = sized_inputs(dir.path(), "chrM", held, length);
        let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", assembly, "--gff3"])
            .arg(&gff3)
            .arg("--genome-fasta")
            .arg(&fasta)
            .arg("--output-dir")
            .arg(dir.path().join("cache"))
            .args(["--codon-table", "MT=2"])
            .output()
            .expect("vep-cache-builder runs");
        let stderr = String::from_utf8_lossy(&output.stderr);
        let case = format!("{assembly}, {held} of {length} bases");
        assert!(output.status.success(), "{case}: {stderr}");
        let (chromosomes, codon_table, coding) = built(dir.path());
        assert_eq!(chromosomes, [named], "{case}");
        assert_eq!(codon_table, 2, "{case}");
        assert_eq!(coding, CODING_BASES, "{case}");
        let notices: Vec<&str> = stderr
            .lines()
            .filter(|line| line.contains("under its written name"))
            .collect();
        if named == "chrM" {
            assert_eq!(notices.len(), 1, "{case}: {stderr}");
            assert!(
                notices[0].ends_with(&kept_as_written(&gff3)),
                "{case}: {}",
                notices[0]
            );
        } else {
            assert_eq!(notices, [] as [&str; 0], "{case}");
        }
    }
}

/// The end of the line a build for GRCh37 writes for the `chrM` of `gff3`
/// that it keeps under its written name.
fn kept_as_written(gff3: &Path) -> String {
    format!(
        " WARN vep_cache_builder: Keeping chrM of {} under its written name: on GRCh37 a \
         gene-model file's M or chrM is read as MT only where the FASTA gives for MT the \
         assembly's mitochondrial sequence, 16,569 bases long, or gives no record for MT and \
         names no organelle genome, and the genome FASTA does neither",
        gff3.display()
    )
}

/// A build for GRCh37 over a genome FASTA that gives no record for `MT`
/// follows the names the FASTA holds, whatever the record's length. Where it
/// names no organelle genome (`21` alone, a part of a genome) a mitochondrion
/// the GFF3 writes `chrM` is built as `MT` and the build says nothing of the
/// name; where it names an organelle genome of its own (`Mito`) the region
/// keeps its written name and the build says so once. Neither FASTA holds
/// the region, so its coding sequence is not read.
#[test]
fn on_grch37_a_fasta_with_no_record_for_mt_gives_mt_only_where_it_names_no_organelle_genome() {
    for (held, length, named) in [
        ("21", 16_569, "MT"),
        ("21", 16_571, "MT"),
        ("Mito", 16_569, "chrM"),
        ("Mito", 16_571, "chrM"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let (gff3, fasta) = sized_inputs(dir.path(), "chrM", held, length);
        let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", "GRCh37", "--gff3"])
            .arg(&gff3)
            .arg("--genome-fasta")
            .arg(&fasta)
            .arg("--output-dir")
            .arg(dir.path().join("cache"))
            .args(["--codon-table", "MT=2"])
            .output()
            .expect("vep-cache-builder runs");
        let stderr = String::from_utf8_lossy(&output.stderr);
        let case = format!("{held} of {length} bases");
        assert!(output.status.success(), "{case}: {stderr}");
        let (chromosomes, codon_table, coding) = built(dir.path());
        assert_eq!(chromosomes, [named], "{case}");
        assert_eq!(codon_table, 2, "{case}");
        assert!(coding.is_null(), "{case}: {coding}");
        let notices: Vec<&str> = stderr
            .lines()
            .filter(|line| line.contains("under its written name"))
            .collect();
        if named == "chrM" {
            assert_eq!(notices.len(), 1, "{case}: {stderr}");
            assert!(
                notices[0].ends_with(&kept_as_written(&gff3)),
                "{case}: {}",
                notices[0]
            );
        } else {
            assert_eq!(notices, [] as [&str; 0], "{case}");
        }
    }
}

/// A region is the sequence a transcript's own line names. In a build for
/// GRCh37 over a `chrM` record of 16,571 bases, a transcript written on
/// `chrM` under a gene line on `chrMT` is built as `chrM` and the build says
/// that it kept the name; a transcript written on `chrMT` under a gene line
/// on `chrM` is built as `MT` and the build names no kept region.
#[test]
fn on_grch37_the_region_kept_as_written_is_the_one_a_transcripts_own_line_names() {
    for (gene_line, transcript_lines, named) in [("chrMT", "chrM", "chrM"), ("chrM", "chrMT", "MT")]
    {
        let dir = tempfile::tempdir().unwrap();
        let (gff3, fasta) = sized_inputs(dir.path(), transcript_lines, "chrM", 16_571);
        let model = std::fs::read_to_string(&gff3).unwrap().replacen(
            &format!("{transcript_lines}\tsynthetic\tgene\t"),
            &format!("{gene_line}\tsynthetic\tgene\t"),
            1,
        );
        assert!(model.contains(&format!("\n{gene_line}\tsynthetic\tgene\t")));
        std::fs::write(&gff3, model).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", "GRCh37", "--gff3"])
            .arg(&gff3)
            .arg("--genome-fasta")
            .arg(&fasta)
            .arg("--output-dir")
            .arg(dir.path().join("cache"))
            .args(["--codon-table", "MT=2"])
            .output()
            .expect("vep-cache-builder runs");
        let stderr = String::from_utf8_lossy(&output.stderr);
        let case = format!("a gene line on {gene_line}, its transcript on {transcript_lines}");
        assert!(output.status.success(), "{case}: {stderr}");
        let (chromosomes, _, _) = built(dir.path());
        assert_eq!(chromosomes, [named], "{case}");
        let kept = format!(
            "Keeping chrM of {} under its written name: ",
            gff3.display()
        );
        let notices = stderr.lines().filter(|line| line.contains(&kept)).count();
        assert_eq!(notices, usize::from(named == "chrM"), "{case}: {stderr}");
    }
}

/// The inputs of [`sized_inputs`] with the record's bases moved two positions
/// on: two bases stand before the fixture's 90, so the bases a model counts
/// on at one position stand two positions further in the record.
fn moved_inputs(dir: &Path, annotated: &str, held: &str, length: usize) -> (PathBuf, PathBuf) {
    let (gff3, fasta) = sized_inputs(dir, annotated, held, length);
    let text = std::fs::read_to_string(&fasta).unwrap();
    let (header, lines) = text.split_once('\n').unwrap();
    let mut bases: Vec<u8> = b"GG".to_vec();
    bases.extend(lines.bytes().filter(|byte| *byte != b'\n'));
    bases.truncate(length);
    let mut moved = format!("{header}\n");
    for line in bases.chunks(60) {
        moved.push_str(std::str::from_utf8(line).unwrap());
        moved.push('\n');
    }
    std::fs::write(&fasta, moved).unwrap();
    (gff3, fasta)
}

/// The bases a record moved two positions on holds at the fixture's coding
/// span: under codon table 2 the first codon, `CCA`, is no start codon and
/// the fifth, `TAA`, is a stop codon.
const MOVED_CODING_BASES: &str = "CCATGCTTTGAATAAAAT";

/// Builds a cache under `dir` for `assembly` from `gff3` and `fasta`, the
/// mitochondrion under codon table 2.
fn build_for(dir: &Path, assembly: &str, gff3: &Path, fasta: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", assembly, "--gff3"])
        .arg(gff3)
        .arg("--genome-fasta")
        .arg(fasta)
        .arg("--output-dir")
        .arg(dir.join("cache"))
        .args(["--codon-table", "MT=2"])
        .output()
        .expect("vep-cache-builder runs")
}

/// On GRCh37 neither a mitochondrion's name nor the length of its FASTA
/// record says which of the two sequences a model counts on, so a build
/// tests the region's coding transcripts on the record they are cut from. A
/// build whose one coding transcript opens on no start codon and holds a
/// stop codon before its last codon there ends with one error that names the
/// region, the model, the record and its length, and writes no cache: over
/// a `chrM` record of 16,571 bases, where the region keeps its written name;
/// over a `chrM` record of the assembly's 16,569 bases, which is read as
/// `MT`; and for a model written `MT` over either record.
#[test]
fn on_grch37_a_mitochondrion_whose_coding_transcripts_do_not_read_on_their_record_is_refused() {
    for (annotated, held, length, region, bases_long) in [
        ("chrM", "chrM", 16_571, "chrM", "16,571"),
        ("chrM", "chrM", 16_569, "MT", "16,569"),
        ("chrM", "M", 16_569, "MT", "16,569"),
        ("MT", "chrM", 16_571, "MT", "16,571"),
        ("MT", "MT", 16_569, "MT", "16,569"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let (gff3, fasta) = moved_inputs(dir.path(), annotated, held, length);
        let output = build_for(dir.path(), "GRCh37", &gff3, &fasta);
        let stderr = String::from_utf8_lossy(&output.stderr);
        let case = format!("a model on {annotated} over {held} of {length} bases");
        assert_eq!(output.status.code(), Some(1), "{case}: {stderr}");
        let refusal = format!(
            "ERROR: [cache-builder] sequence region {region} of {} does not fit the genome FASTA: \
             of 1 coding transcript read from the FASTA record {held}, {bases_long} bases long, \
             under codon table 2, none opens on a start codon and runs to its last codon without \
             a stop codon. On GRCh37 two mitochondrial sequences are in use, the assembly's of \
             16,569 bases and the 16,571-base one of the UCSC hg19 naming: build from a genome \
             FASTA that holds the sequence this model counts on, or from a model without the \
             region's lines",
            gff3.display()
        );
        let errors: Vec<&str> = stderr
            .lines()
            .filter(|line| line.starts_with("ERROR:"))
            .collect();
        assert_eq!(errors, [refusal], "{case}: {stderr}");
        assert!(!dir.path().join("cache").exists(), "{case}");
    }
}

/// The test is GRCh37's alone: on another assembly the names are one
/// sequence, and a build over the same moved record writes the region as
/// `MT` with the coding sequence the record gives.
#[test]
fn on_another_assembly_a_mitochondrion_is_built_from_the_record_as_it_stands() {
    for assembly in ["GRCh38", "made_up"] {
        let dir = tempfile::tempdir().unwrap();
        let (gff3, fasta) = moved_inputs(dir.path(), "chrM", "chrM", 16_569);
        let output = build_for(dir.path(), assembly, &gff3, &fasta);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{assembly}: {stderr}");
        let (chromosomes, codon_table, coding) = built(dir.path());
        assert_eq!(chromosomes, ["MT"], "{assembly}");
        assert_eq!(codon_table, 2, "{assembly}");
        assert_eq!(coding, MOVED_CODING_BASES, "{assembly}");
    }
}

/// A transcript flagged `cds_start_NF` opens on no start codon by its own
/// statement, so the test does not read it: the same GRCh37 build, given a
/// GTF whose `transcript` line tags the model's one coding transcript so,
/// tests no transcript and writes the region from the record as it stands.
#[test]
fn on_grch37_a_transcript_flagged_as_incomplete_at_its_start_is_not_tested() {
    let dir = tempfile::tempdir().unwrap();
    let (gff3, fasta) = moved_inputs(dir.path(), "chrM", "chrM", 16_569);
    let tags = dir.path().join("tags.gtf");
    std::fs::write(
        &tags,
        "chrM\tsynthetic\ttranscript\t1\t30\t.\t+\t.\tgene_id \"MITOG1\"; \
         transcript_id \"MITOT1\"; transcript_version \"1\"; tag \"cds_start_NF\";\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh37", "--gff3"])
        .arg(&gff3)
        .arg("--gtf")
        .arg(&tags)
        .arg("--genome-fasta")
        .arg(&fasta)
        .arg("--output-dir")
        .arg(dir.path().join("cache"))
        .args(["--codon-table", "MT=2"])
        .output()
        .expect("vep-cache-builder runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let (chromosomes, codon_table, coding) = built(dir.path());
    assert_eq!(chromosomes, ["MT"]);
    assert_eq!(codon_table, 2);
    assert_eq!(coding, MOVED_CODING_BASES);
    let shard = dir.path().join("cache/transcripts/MT/1-1000000.json");
    let records: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(shard).unwrap()).unwrap();
    let codes: Vec<&str> = records[0]["attributes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|attribute| attribute["code"].as_str().unwrap())
        .collect();
    assert!(codes.contains(&"cds_start_NF"), "{codes:?}");
}

/// The tag of the model's own transcript line passes the transcript over
/// the same way: a GRCh37 build from a GFF3 whose `mRNA` line carries
/// `tag=cds_start_NF`, over the moved record, tests no transcript and
/// writes the region from the record as it stands. The tag gives the
/// record no attribute.
#[test]
fn on_grch37_a_transcript_its_own_line_flags_as_incomplete_at_its_start_is_not_tested() {
    let dir = tempfile::tempdir().unwrap();
    let (gff3, fasta) = moved_inputs(dir.path(), "chrM", "chrM", 16_569);
    let model = std::fs::read_to_string(&gff3).unwrap();
    let flagged = model.replace(
        "Parent=gene:MITOG1;biotype=protein_coding;version=1\n",
        "Parent=gene:MITOG1;biotype=protein_coding;version=1;tag=cds_start_NF\n",
    );
    assert_ne!(flagged, model);
    std::fs::write(&gff3, flagged).unwrap();
    let output = build_for(dir.path(), "GRCh37", &gff3, &fasta);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let (chromosomes, codon_table, coding) = built(dir.path());
    assert_eq!(chromosomes, ["MT"]);
    assert_eq!(codon_table, 2);
    assert_eq!(coding, MOVED_CODING_BASES);
    let shard = dir.path().join("cache/transcripts/MT/1-1000000.json");
    let records: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(shard).unwrap()).unwrap();
    assert!(records[0].get("attributes").is_none(), "{}", records[0]);
}
