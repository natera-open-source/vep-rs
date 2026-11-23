// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! An installed cache starts and ends a one-exon transcript's coding sequence
//! where its row of the release's `translation` table does.
//!
//! A release's GFF3 can give a transcript a `CDS` line over its whole exon
//! beside a `five_prime_UTR` line that overlaps it, where the `translation`
//! row of the same release starts the coding sequence after that untranslated
//! stretch (`seq_start`, a position within the exon `start_exon_id` names).
//! Ensembl VEP's cache is dumped from the database and holds the row's
//! positions. `--auto` reads the core dumps, so the cache it installs holds
//! them too: the transcript's `cdna_coding_start`, its translation's `start`
//! and `end`, its untranslated sequence, its translateable sequence and the
//! protein FASTA record that is that sequence's peptide. The dumps name an
//! exon by its database identifier alone, so a transcript of several exons
//! whose row differs from its CDS lines is built as the lines place it and
//! counted in one warning.
//!
//! Each test lays a made-up genome out as Ensembl Genomes release 63
//! publishes a species of the plants division and installs it from that
//! directory.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Output};

use serde_json::Value;

/// The BSD `sum` of `bytes`: the 16-bit checksum and the 1 KiB block count a
/// `CHECKSUMS` line carries.
fn bsd_sum(bytes: &[u8]) -> (u16, u64) {
    let sum = bytes.iter().fold(0u16, |sum, &b| {
        sum.rotate_right(1).wrapping_add(u16::from(b))
    });
    (sum, (bytes.len() as u64).div_ceil(1024))
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// Writes `dir/name` under `site` and its line of the directory's
/// `CHECKSUMS`.
fn publish(site: &Path, dir: &str, name: &str, bytes: &[u8]) {
    let dir = site.join(dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join(name), bytes).unwrap();
    let (sum, blocks) = bsd_sum(bytes);
    let mut checksums = std::fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(dir.join("CHECKSUMS"))
        .unwrap();
    writeln!(checksums, "{sum:05} {blocks} {name}").unwrap();
}

const PLANT: &str = "arabidopsis_thaliana";
const DATABASE: &str = "arabidopsis_thaliana_core_63_116_11";
const PLANTS: &str = "pub/release-63/plants";

/// The reverse complement of `bases`.
fn reverse_complement(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => panic!("no complement of {}", other as char),
        })
        .collect()
}

/// The made-up bases of sequence region `1` from 1-based `start` to `end`,
/// as [`genome`] lays them where no text is placed.
fn filler(start: usize, end: usize) -> String {
    (start - 1..end)
        .map(|offset| b"ACGT"[(offset * 7 + offset / 5) % 4] as char)
        .collect()
}

/// The 1,000 bases of sequence region `1`: [`filler`] with `placed` written
/// over it, each text at its 1-based start.
fn genome(placed: &[(usize, String)]) -> Vec<u8> {
    let mut bases = filler(1, 1_000).into_bytes();
    for (start, text) in placed {
        bases[start - 1..start - 1 + text.len()].copy_from_slice(text.as_bytes());
    }
    let mut fasta = b">1 dna:chromosome\n".to_vec();
    for line in bases.chunks(60) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    fasta
}

/// One made-up gene of one transcript on sequence region `1`.
struct Gene {
    /// The number its identifiers end in and its rows of the tables are
    /// keyed by: gene `MADEUPG<n>`, transcript `MADEUPT<n>`, protein
    /// `MADEUPP<n>`.
    number: u32,
    strand: i8,
    /// The exons in transcript order, each `(start, end)`.
    exons: Vec<(u64, u64)>,
    /// The `CDS` lines, each `(start, end)`.
    cds: Vec<(u64, u64)>,
    /// The `five_prime_UTR` line, when the GFF3 writes one.
    five_prime_utr: Option<(u64, u64)>,
    /// The translation row's `seq_start` and `seq_end`; its start exon is the
    /// first exon and its end exon the last.
    row: (u64, u64),
    /// The protein FASTA record's residues.
    residues: String,
}

impl Gene {
    fn span(&self) -> (u64, u64) {
        let start = self.exons.iter().map(|e| e.0).min().unwrap();
        let end = self.exons.iter().map(|e| e.1).max().unwrap();
        (start, end)
    }

    fn sign(&self) -> char {
        if self.strand < 0 {
            '-'
        } else {
            '+'
        }
    }

    /// The gene's lines as an Ensembl GFF3 writes them: the gene, the
    /// transcript, the untranslated stretch, then each exon with its rank and
    /// its phases, the `ensembl_phase` of an exon that opens untranslated
    /// being -1, and the `CDS` lines.
    fn gff3(&self) -> String {
        let n = self.number;
        let (start, end) = self.span();
        let sign = self.sign();
        let mut lines = format!(
            "1\tmade_up\tgene\t{start}\t{end}\t.\t{sign}\t.\tID=gene:MADEUPG{n};biotype=protein_coding;gene_id=MADEUPG{n}\n\
             1\tmade_up\tmRNA\t{start}\t{end}\t.\t{sign}\t.\tID=transcript:MADEUPT{n};Parent=gene:MADEUPG{n};biotype=protein_coding;transcript_id=MADEUPT{n}\n"
        );
        if let Some((start, end)) = self.five_prime_utr {
            lines.push_str(&format!(
                "1\tmade_up\tfive_prime_UTR\t{start}\t{end}\t.\t{sign}\t.\tParent=transcript:MADEUPT{n}\n"
            ));
        }
        for (rank, (start, end)) in self.exons.iter().enumerate() {
            let rank = rank + 1;
            let phase = if rank == 1 { -1 } else { 0 };
            lines.push_str(&format!(
                "1\tmade_up\texon\t{start}\t{end}\t.\t{sign}\t.\tParent=transcript:MADEUPT{n};ensembl_end_phase=0;ensembl_phase={phase};exon_id=MADEUPT{n}.exon{rank};rank={rank}\n"
            ));
        }
        for (start, end) in &self.cds {
            lines.push_str(&format!(
                "1\tmade_up\tCDS\t{start}\t{end}\t.\t{sign}\t0\tID=CDS:MADEUPP{n};Parent=transcript:MADEUPT{n};protein_id=MADEUPP{n}\n"
            ));
        }
        lines
    }

    fn record(&self) -> String {
        let n = self.number;
        let (start, end) = self.span();
        format!(
            ">MADEUPP{n} pep chromosome:TAIR10:1:{start}:{end}:{} gene:MADEUPG{n} transcript:MADEUPT{n} gene_biotype:protein_coding transcript_biotype:protein_coding\n{}\n",
            self.strand, self.residues
        )
    }
}

/// Lays `genes` out as a species of the plants division: the species list,
/// the GFF3, a GTF of header lines, the protein FASTA, the genome and the
/// core tables, whose `translation` rows carry each gene's `row` and whose
/// other tables hold the rows an installation reads.
fn site(genes: &[Gene], placed: &[(usize, String)]) -> tempfile::TempDir {
    let root = tempfile::tempdir().unwrap();
    let site = root.path();
    publish(
        site,
        "pub/release-63",
        "species.txt",
        format!(
            "#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\tmicroarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\tspecies_id\n\
             A genome\t{PLANT}\tEnsemblPlants\t1\tasm\tGCA_000000000.1\t2020-01\tN\tN\tN\tN\tN\tN\t{DATABASE}\t1\n"
        )
        .as_bytes(),
    );
    let gff3: String = std::iter::once("##gff-version 3\n".to_string())
        .chain(genes.iter().map(Gene::gff3))
        .collect();
    publish(
        site,
        &format!("{PLANTS}/gff3/{PLANT}"),
        "Arabidopsis_thaliana.TAIR10.63.gff3.gz",
        &gzip(gff3.as_bytes()),
    );
    publish(
        site,
        &format!("{PLANTS}/gtf/{PLANT}"),
        "Arabidopsis_thaliana.TAIR10.63.gtf.gz",
        &gzip(b"#!genome-build made_up\n#!genome-version made_up\n"),
    );
    let records: String = genes.iter().map(Gene::record).collect();
    publish(
        site,
        &format!("{PLANTS}/fasta/{PLANT}/pep"),
        "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz",
        &gzip(records.as_bytes()),
    );
    publish(
        site,
        &format!("{PLANTS}/fasta/{PLANT}/dna"),
        "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
        &gzip(&genome(placed)),
    );

    let mut gene_rows = String::new();
    let mut transcript_rows = String::new();
    let mut translation_rows = String::new();
    for gene in genes {
        let n = gene.number;
        let (start, end) = gene.span();
        let strand = gene.strand;
        // gene: its canonical transcript; transcript: its canonical
        // translation; translation: seq_start, start_exon_id, seq_end,
        // end_exon_id, the exon identifiers numbers no other table of the
        // installation carries.
        gene_rows.push_str(&format!(
            "{n}\tprotein_coding\t1\t1\t{start}\t{end}\t{strand}\t\\N\tmade_up\t\\N\t1\t{n}1\tMADEUPG{n}\t\\N\t\\N\t\\N\n"
        ));
        transcript_rows.push_str(&format!(
            "{n}1\t{n}\t1\t1\t{start}\t{end}\t{strand}\t\\N\tmade_up\tprotein_coding\t\\N\t1\t{n}2\tMADEUPT{n}\t\\N\t\\N\t\\N\n"
        ));
        let (seq_start, seq_end) = gene.row;
        let start_exon = format!("{n}01");
        let end_exon = format!("{n}{:02}", gene.exons.len());
        translation_rows.push_str(&format!(
            "{n}2\t{n}1\t{seq_start}\t{start_exon}\t{seq_end}\t{end_exon}\tMADEUPP{n}\t\\N\t\\N\t\\N\n"
        ));
    }
    let tables: [(&str, &str); 15] = [
        ("external_db", ""),
        (
            "attrib_type",
            "1\tcodon_table\tCodon table\t\\N\n2\tkaryotype_rank\tRank\t\\N\n",
        ),
        ("seq_region", "1\t1\t1\t1000\n"),
        ("seq_region_attrib", "1\t2\t1\n"),
        ("seq_region_synonym", ""),
        ("meta", "1\t1\tassembly.name\tTAIR10\n"),
        ("gene", &gene_rows),
        ("transcript", &transcript_rows),
        ("xref", ""),
        ("object_xref", ""),
        ("translation", &translation_rows),
        ("translation_attrib", ""),
        ("transcript_attrib", ""),
        ("analysis_description", ""),
        ("protein_feature", ""),
    ];
    for (table, rows) in tables {
        publish(
            site,
            &format!("{PLANTS}/mysql/{DATABASE}"),
            &format!("{table}.txt.gz"),
            &gzip(rows.as_bytes()),
        );
    }
    root
}

/// Installs the species from `site` into `dir_cache` and returns the run.
fn install(site: &Path, dir_cache: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--auto", PLANT, "TAIR10", "--release", "63"])
        .arg("--dir_cache")
        .arg(dir_cache)
        .arg("--site")
        .arg(site)
        .output()
        .expect("vep-cache-builder runs")
}

/// The text of each line `module` logs at `level` in `log`, the line's time,
/// level and module taken off.
fn logged<'a>(log: &'a str, level: &str, module: &str) -> Vec<&'a str> {
    let opening = format!(" {level} {module}: ");
    log.lines()
        .filter_map(|line| line.split_once(&opening))
        .map(|(_, text)| text)
        .collect()
}

/// Every transcript of an installed cache by stable identifier.
fn transcripts(installed: &Path) -> std::collections::BTreeMap<String, Value> {
    let mut transcripts = std::collections::BTreeMap::new();
    for region in std::fs::read_dir(installed.join("transcripts")).unwrap() {
        for shard in std::fs::read_dir(region.unwrap().path()).unwrap() {
            let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
            let records: Vec<Value> = serde_json::from_str(&text).unwrap();
            for record in records {
                transcripts.insert(record["stable_id"].as_str().unwrap().to_string(), record);
            }
        }
    }
    transcripts
}

/// Where a cached transcript places its coding sequence and what it reads
/// there: `cdna_coding_start` and `cdna_coding_end`, on the record and in its
/// mapper, the translation's `start` and `end`, the 5' untranslated
/// sequence, the translateable sequence and the peptide.
fn coding(transcript: &Value) -> Value {
    let cache = &transcript["variation_effect_feature_cache"];
    serde_json::json!({
        "cdna_coding_start": transcript["cdna_coding_start"],
        "cdna_coding_end": transcript["cdna_coding_end"],
        "mapper_cdna_coding_start": cache["mapper"]["cdna_coding_start"],
        "mapper_cdna_coding_end": cache["mapper"]["cdna_coding_end"],
        "translation_start": transcript["translation"]["start"],
        "translation_end": transcript["translation"]["end"],
        "five_prime_utr": cache["five_prime_utr"],
        "three_prime_utr": cache["three_prime_utr"],
        "translateable_seq": cache["translateable_seq"],
        "peptide": cache["peptide"],
    })
}

/// The rows of `feature` that `vep` writes for `records`, VCF data lines, from
/// the cache installed under `dir_cache`, annotated in-process as the `vep`
/// binary would.
fn rows_of(feature: &str, records: &str, dir_cache: &Path, work: &Path) -> Vec<String> {
    let input = work.join("variants.vcf");
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"),
    )
    .unwrap();
    let output = work.join("annotated.txt");
    let (input, output_s, dir_cache) = (
        input.display().to_string(),
        output.display().to_string(),
        dir_cache.display().to_string(),
    );
    let argv = vec![
        "vep",
        "--offline",
        "--cache",
        "--no_stats",
        "--quiet",
        "--force_overwrite",
        "-i",
        &input,
        "-o",
        &output_s,
        "--dir_cache",
        &dir_cache,
        "--species",
        PLANT,
        "--cache_version",
        "63",
        "--assembly",
        "TAIR10",
    ];
    let config = vep_cli::Config::from_argv(argv).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
    std::fs::read_to_string(output)
        .unwrap()
        .lines()
        .filter(|line| line.split('\t').nth(4) == Some(feature))
        .map(str::to_string)
        .collect()
}

/// Gene 1, on the forward strand: one exon at 101 to 160 whose first 18 bases
/// are untranslated, a `CDS` line over the whole exon beside the
/// `five_prime_UTR` line, and a translation row that starts at base 19 of the
/// exon. Bases 119 to 160 are a start codon, twelve alanine codons and a stop
/// codon.
fn forward_gene() -> (Gene, (usize, String)) {
    let exon = format!("{}ATG{}TAA", "CCT".repeat(6), "GCT".repeat(12));
    (
        Gene {
            number: 1,
            strand: 1,
            exons: vec![(101, 160)],
            cds: vec![(101, 160)],
            five_prime_utr: Some((101, 118)),
            row: (19, 60),
            residues: format!("M{}", "A".repeat(12)),
        },
        (101, exon),
    )
}

/// Gene 2, on the reverse strand: one exon at 301 to 360 whose first 15
/// bases in transcript order, 346 to 360, are untranslated, a `CDS` line
/// over the whole exon and a translation row that starts at base 16 of the
/// exon. In transcript order the coding bases are a start codon, thirteen
/// lysine codons and a stop codon.
fn reverse_gene() -> (Gene, (usize, String)) {
    let exon = format!("{}ATG{}TAG", "GGA".repeat(5), "AAA".repeat(13));
    (
        Gene {
            number: 2,
            strand: -1,
            exons: vec![(301, 360)],
            cds: vec![(301, 360)],
            five_prime_utr: Some((346, 360)),
            row: (16, 60),
            residues: format!("M{}", "K".repeat(13)),
        },
        (301, reverse_complement(&exon)),
    )
}

/// Gene 3, whose translation row and `CDS` line agree: one exon at 501 to
/// 560, the `CDS` line at 507 to 548 and the row at bases 7 to 48 of the exon.
/// The coding bases are a start codon, twelve glutamate codons and a stop
/// codon.
fn agreeing_gene() -> (Gene, (usize, String)) {
    let coding = format!("ATG{}TGA", "GAA".repeat(12));
    (
        Gene {
            number: 3,
            strand: 1,
            exons: vec![(501, 560)],
            cds: vec![(507, 548)],
            five_prime_utr: Some((501, 506)),
            row: (7, 48),
            residues: format!("M{}", "E".repeat(12)),
        },
        (507, coding),
    )
}

/// The installation of a one-exon transcript on each strand whose `CDS` line
/// covers the exon while the translation row starts after the untranslated
/// bases, beside a transcript whose line and row agree.
///
/// Each of the two is cached as its row places it: the coding sequence from
/// base 19 (forward) and base 16 (reverse) of the cDNA to its end, the
/// translation's `start` the row's `seq_start`, the bases before it the 5'
/// untranslated sequence, and the peptide the protein FASTA record, whose
/// length is the row's coding sequence's. No record is set aside, one line
/// counts the two transcripts, and the agreeing transcript is cached as its
/// line places it. `vep` then writes a record in the untranslated bases as
/// a 5' UTR variant and counts a coding record's positions from the row's
/// start, on both strands.
#[test]
fn a_one_exon_transcript_is_cached_as_its_translation_row_places_its_coding_sequence() {
    let (forward, forward_bases) = forward_gene();
    let (reverse, reverse_bases) = reverse_gene();
    let (agreeing, agreeing_bases) = agreeing_gene();
    let site = site(
        &[forward, reverse, agreeing],
        &[forward_bases, reverse_bases, agreeing_bases],
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = install(site.path(), &dir_cache);
    let log = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{log}");
    let cached = transcripts(&dir_cache.join("arabidopsis_thaliana/63_TAIR10"));
    assert_eq!(cached.len(), 3);
    assert_eq!(
        coding(&cached["MADEUPT1"]),
        serde_json::json!({
            "cdna_coding_start": 19,
            "cdna_coding_end": 60,
            "mapper_cdna_coding_start": 19,
            "mapper_cdna_coding_end": 60,
            "translation_start": 19,
            "translation_end": 60,
            "five_prime_utr": "CCT".repeat(6),
            "three_prime_utr": null,
            "translateable_seq": format!("ATG{}TAA", "GCT".repeat(12)),
            "peptide": format!("M{}", "A".repeat(12)),
        })
    );
    assert_eq!(
        coding(&cached["MADEUPT2"]),
        serde_json::json!({
            "cdna_coding_start": 16,
            "cdna_coding_end": 60,
            "mapper_cdna_coding_start": 16,
            "mapper_cdna_coding_end": 60,
            "translation_start": 16,
            "translation_end": 60,
            "five_prime_utr": "GGA".repeat(5),
            "three_prime_utr": null,
            "translateable_seq": format!("ATG{}TAG", "AAA".repeat(13)),
            "peptide": format!("M{}", "K".repeat(13)),
        })
    );
    assert_eq!(
        coding(&cached["MADEUPT3"]),
        serde_json::json!({
            "cdna_coding_start": 7,
            "cdna_coding_end": 48,
            "mapper_cdna_coding_start": 7,
            "mapper_cdna_coding_end": 48,
            "translation_start": 7,
            "translation_end": 48,
            "five_prime_utr": filler(501, 506),
            "three_prime_utr": filler(549, 560),
            "translateable_seq": format!("ATG{}TGA", "GAA".repeat(12)),
            "peptide": format!("M{}", "E".repeat(12)),
        })
    );

    // Base 110 is the tenth of the forward exon, untranslated; base 125 is
    // the seventh coding base, the first of the third codon.
    assert_eq!(
        rows_of(
            "MADEUPT1",
            "1\t110\t.\tC\tA\t.\t.\t.\n1\t125\t.\tG\tA\t.\t.\t.\n",
            &dir_cache,
            root.path()
        ),
        [
            "1_110_C/A\t1:110\tA\tMADEUPG1\tMADEUPT1\tTranscript\t5_prime_UTR_variant\t10\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1",
            "1_125_G/A\t1:125\tA\tMADEUPG1\tMADEUPT1\tTranscript\tmissense_variant\t25\t7\t3\tA/T\tGct/Act\t-\tIMPACT=MODERATE;STRAND=1",
        ]
    );
    // Base 345 is the sixteenth of the reverse exon in transcript order, the
    // first base of its start codon.
    assert_eq!(
        rows_of(
            "MADEUPT2",
            "1\t345\t.\tT\tC\t.\t.\t.\n",
            &dir_cache,
            root.path()
        ),
        ["1_345_T/C\t1:345\tC\tMADEUPG2\tMADEUPT2\tTranscript\tstart_lost\t16\t1\t1\tM/V\tAtg/Gtg\t-\tIMPACT=HIGH;STRAND=-1"]
    );

    // The build's lines: the two transcripts counted once, no record set aside.
    assert_eq!(
        logged(&log, "INFO", "vep_cache_builder::builder")
            .into_iter()
            .filter(|text| text.contains("translation table"))
            .collect::<Vec<_>>(),
        [
            "2 one-exon transcripts have a coding sequence that the translation table of the core \
          dumps places otherwise than the model's CDS lines: each is built as the table places it"
        ],
        "{log}"
    );
    assert!(!log.contains("set aside for"), "{log}");
    assert!(!log.contains("no record in the protein FASTA"), "{log}");
}

/// A transcript of two exons whose `CDS` lines cover both while the
/// translation row starts at base 10 of its start exon: the dumps name that
/// exon by a database identifier the model does not carry, so the transcript
/// is cached as its lines place the coding sequence, from the first base of
/// the cDNA, and one warning counts it. The agreeing transcript beside it
/// draws no line.
#[test]
fn a_transcript_of_several_exons_is_cached_as_its_cds_lines_place_it_and_counted() {
    let (agreeing, agreeing_bases) = agreeing_gene();
    let spliced = Gene {
        number: 4,
        strand: 1,
        exons: vec![(701, 730), (761, 790)],
        cds: vec![(701, 730), (761, 790)],
        five_prime_utr: Some((701, 709)),
        row: (10, 30),
        residues: format!("M{}", "G".repeat(15)),
    };
    let site = site(&[agreeing, spliced], &[agreeing_bases]);
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = install(site.path(), &dir_cache);
    let log = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{log}");
    assert_eq!(
        logged(&log, "WARN", "vep_cache_builder::builder")
            .into_iter()
            .filter(|text| text.contains("translation table"))
            .collect::<Vec<_>>(),
        ["1 transcript has a coding sequence that the translation table of the core dumps places \
          otherwise than the model's CDS lines, at a place its exon identifiers alone do not give: \
          it is built as the CDS lines place it"],
        "{log}"
    );
    assert!(
        logged(&log, "INFO", "vep_cache_builder::builder")
            .iter()
            .all(|text| !text.contains("translation table")),
        "{log}"
    );

    let cached = transcripts(&dir_cache.join("arabidopsis_thaliana/63_TAIR10"));
    let spliced = &cached["MADEUPT4"];
    assert_eq!(spliced["cdna_coding_start"], 1, "{spliced}");
    assert_eq!(spliced["cdna_coding_end"], 60, "{spliced}");
    assert_eq!(spliced["translation"]["start"], 1, "{spliced}");
    assert_eq!(spliced["translation"]["end"], 30, "{spliced}");
    assert_eq!(cached["MADEUPT3"]["cdna_coding_start"], 7);
}
