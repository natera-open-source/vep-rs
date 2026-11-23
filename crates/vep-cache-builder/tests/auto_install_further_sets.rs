// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--auto <species> <assembly>` reaches the genomes published outside the
//! Ensembl release tree, each from a local directory laid out as its site.
//!
//! A species of an Ensembl Genomes division installs from
//! `pub/release-<N>/<division>/`, `N` the Ensembl Genomes release: the
//! release's species list names its division and its core database, every
//! input is checked against its directory's `CHECKSUMS`, and the cache lands
//! at `<dir_cache>/<species>/<N>_<assembly>`, where `vep` resolves it through
//! `--dir_cache`. A genome of a collection and an assembly of the rapid
//! release are refused before anything is fetched, each with one line naming
//! what its tree lacks and the routes that build its cache.
//!
//! The Ensembl Genomes tree holds the GFF3 lines and peptides of
//! `tests/fixtures/tair10_slice/` with made-up table dumps, GTF and genome
//! sequence; the release tree holds the chromosome 21 fixture of
//! `tests/fixtures/grch38_chr21/` under another species' names.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles::bgzf;
use serde_json::Value;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn golden(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{relative}"))
}

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

/// A directory laid out as a site, each directory carrying the `CHECKSUMS`
/// its files sum to.
struct Site {
    root: tempfile::TempDir,
}

impl Site {
    fn new() -> Self {
        Self {
            root: tempfile::tempdir().unwrap(),
        }
    }

    fn path(&self) -> &Path {
        self.root.path()
    }

    /// Writes `dir/name` and its `CHECKSUMS` line; `sum` overrides the line's
    /// checksum so a corrupt listing can be laid out.
    fn add_with_sum(&self, dir: &str, name: &str, bytes: &[u8], sum: Option<(u16, u64)>) {
        let dir = self.path().join(dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(name), bytes).unwrap();
        let (sum, blocks) = sum.unwrap_or_else(|| bsd_sum(bytes));
        let mut checksums = std::fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(dir.join("CHECKSUMS"))
            .unwrap();
        writeln!(checksums, "{sum:05} {blocks} {name}").unwrap();
    }

    fn add(&self, dir: &str, name: &str, bytes: &[u8]) {
        self.add_with_sum(dir, name, bytes, None);
    }
}

/// The header line of a release's species list.
const SPECIES_LIST_HEADER: &str = "#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\tmicroarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\tspecies_id\n";

/// A species list row: the species, its division and its core database.
fn species_row(species: &str, division: &str, core_db: &str) -> String {
    format!("A genome\t{species}\t{division}\t1\tasm\tGCA_000000000.1\t2020-01\tN\tN\tN\tN\tN\tN\t{core_db}\t1\n")
}

/// The base at 0-based `offset` of every made-up sequence region.
fn base(offset: usize) -> u8 {
    b"ACGT"[(offset * 7 + offset / 5) % 4]
}

/// A FASTA of made-up sequence, each record of the length named, in 60-base
/// lines.
fn genome_fasta(records: &[(&str, usize)]) -> Vec<u8> {
    let mut fasta = Vec::new();
    for (name, length) in records {
        writeln!(fasta, ">{name} dna:chromosome").unwrap();
        let bases: Vec<u8> = (0..*length).map(base).collect();
        for line in bases.chunks(60) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
    }
    fasta
}

const PLANT: &str = "arabidopsis_thaliana";
const PLANT_DATABASE: &str = "arabidopsis_thaliana_core_63_116_11";
const PLANTS: &str = "pub/release-63/plants";

/// The made-up core tables of the plant: the two genes of the slice, the
/// first with a display xref and both with a canonical transcript; sequence
/// region `1` with a karyotype rank and a synonym; `Pt` with no codon table;
/// and the assembly's name in `meta`. The other tables are empty.
fn plant_tables() -> Vec<(&'static str, &'static str)> {
    vec![
        (
            "external_db",
            "50\tTAIR_SYMBOL\t\\N\tXREF\t50\tTAIR symbol\tMISC\t\\N\t\\N\t\\N\n",
        ),
        (
            "attrib_type",
            "1\tcodon_table\tCodon table\t\\N\n2\tkaryotype_rank\tRank\t\\N\n",
        ),
        ("seq_region", "1\t1\t1\t30427671\n2\tPt\t1\t154478\n"),
        ("seq_region_attrib", "1\t2\t1\n"),
        ("seq_region_synonym", "1\t1\tChr1\t\\N\n"),
        ("meta", "1\t1\tassembly.name\tTAIR10\n"),
        (
            "gene",
            "1\tprotein_coding\t1\t1\t3631\t5899\t1\t900\ttair\t\\N\t1\t11\tAT1G01010\t\\N\t\\N\t\\N\n\
             2\tprotein_coding\t1\t1\t6788\t9130\t-1\t\\N\ttair\t\\N\t1\t21\tAT1G01020\t\\N\t\\N\t\\N\n",
        ),
        (
            "transcript",
            "11\t1\t1\t1\t3631\t5899\t1\t\\N\ttair\tprotein_coding\t\\N\t1\t\\N\tAT1G01010.1\t\\N\t\\N\t\\N\n\
             21\t2\t1\t1\t6788\t9130\t-1\t\\N\ttair\tprotein_coding\t\\N\t1\t\\N\tAT1G01020.1\t\\N\t\\N\t\\N\n",
        ),
        ("xref", "900\t50\tAT1G01010\tNAC001\t0\t\\N\tDIRECT\t\\N\n"),
        ("object_xref", ""),
        ("translation", ""),
        ("translation_attrib", ""),
        ("transcript_attrib", ""),
        ("analysis_description", ""),
        ("protein_feature", ""),
    ]
}

/// A made-up gene on the plastid, appended to the slice's GFF3: one coding
/// exon on sequence region `Pt`.
const PLASTID_GENE: &str = "Pt\tmade_up\tgene\t100\t400\t.\t+\t.\tID=gene:MADEUPG1;biotype=protein_coding\n\
    Pt\tmade_up\tmRNA\t100\t400\t.\t+\t.\tID=transcript:MADEUPT1;Parent=gene:MADEUPG1;biotype=protein_coding;transcript_id=MADEUPT1\n\
    Pt\tmade_up\texon\t100\t400\t.\t+\t.\tParent=transcript:MADEUPT1\n\
    Pt\tmade_up\tCDS\t100\t399\t.\t+\t0\tParent=transcript:MADEUPT1;protein_id=MADEUPP1\n";

/// Lays the plant out as Ensembl Genomes release 63 publishes a species of the
/// plants division: the species list at the release root, the GFF3 (`gff3`),
/// a GTF of header lines, the slice's peptides, the genome under `dna/` and
/// the core database the list names. The GTF's header names another genome
/// build, so the `meta` table alone gives the cache the assembly's name.
fn plant_site(gff3: &[u8]) -> Site {
    let site = Site::new();
    site.add(
        "pub/release-63",
        "species.txt",
        format!(
            "{SPECIES_LIST_HEADER}{}",
            species_row(PLANT, "EnsemblPlants", PLANT_DATABASE)
        )
        .as_bytes(),
    );
    site.add(
        &format!("{PLANTS}/gff3/{PLANT}"),
        "Arabidopsis_thaliana.TAIR10.63.gff3.gz",
        &gzip(gff3),
    );
    site.add(
        &format!("{PLANTS}/gtf/{PLANT}"),
        "Arabidopsis_thaliana.TAIR10.63.gtf.gz",
        &gzip(b"#!genome-build made_up\n#!genome-version made_up\n"),
    );
    site.add(
        &format!("{PLANTS}/fasta/{PLANT}/pep"),
        "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz",
        &gzip(&std::fs::read(fixture("tair10_slice/peptides.fa")).unwrap()),
    );
    site.add(
        &format!("{PLANTS}/fasta/{PLANT}/dna"),
        "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
        &gzip(&genome_fasta(&[("1", 10_000), ("Pt", 1_000)])),
    );
    for (table, rows) in plant_tables() {
        site.add(
            &format!("{PLANTS}/mysql/{PLANT_DATABASE}"),
            &format!("{table}.txt.gz"),
            &gzip(rows.as_bytes()),
        );
    }
    site
}

fn slice_gff3() -> Vec<u8> {
    std::fs::read(fixture("tair10_slice/slice.gff3")).unwrap()
}

/// Runs `--auto <species> <assembly> --release <release>` against `site` into
/// `dir_cache`, with `extra` flags.
fn auto_install(
    site: &Site,
    species: &str,
    assembly: &str,
    release: &str,
    dir_cache: &Path,
    extra: &[&str],
) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--auto", species, assembly, "--release", release])
        .arg("--dir_cache")
        .arg(dir_cache)
        .arg("--site")
        .arg(site.path())
        .args(extra)
        .output()
        .expect("vep-cache-builder runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The one `ERROR:` line of a failed run.
fn error_line(output: &Output) -> String {
    assert!(!output.status.success(), "the run succeeded");
    let text = stderr(output);
    let lines: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("ERROR: "))
        .collect();
    assert_eq!(lines.len(), 1, "{text}");
    lines[0].to_string()
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

/// Annotates `input` in-process with `flags`, as the `vep` binary would, and
/// returns the output's data lines (every line not opening with `#`).
fn annotate(flags: &[&str], input: &Path, out: &Path) -> Vec<String> {
    let input = input.display().to_string();
    let out_s = out.display().to_string();
    let mut argv = vec![
        "vep",
        "--offline",
        "--cache",
        "--no_stats",
        "--quiet",
        "--force_overwrite",
        "-i",
        &input,
        "-o",
        &out_s,
    ];
    argv.extend_from_slice(flags);
    let config = vep_cli::Config::from_argv(argv).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// A species of the plants division installs from the Ensembl Genomes layout
/// under its Ensembl Genomes release. The species list names its core
/// database, `arabidopsis_thaliana_core_63_116_11`, a number the assembly
/// name `TAIR10` does not end in, and the tables of that directory reach the
/// cache: the display xref's symbol and source, the canonical flag, the
/// assembly's name and the chromosome synonym. The cache lands at
/// `<dir_cache>/arabidopsis_thaliana/63_TAIR10` with the `dna/` genome
/// recompressed and indexed inside it, the work directory is gone, and `vep`
/// resolves the cache and its FASTA through
/// `--dir_cache --species --cache_version --assembly`.
#[test]
fn an_ensembl_genomes_species_installs_from_its_division_tree() {
    let site = plant_site(&slice_gff3());
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, PLANT, "TAIR10", "63", &dir_cache, &[]);
    let log = stderr(&output);
    assert!(output.status.success(), "{log}");
    assert!(
        log.contains(&format!("{PLANTS}/mysql/{PLANT_DATABASE}")),
        "{log}"
    );

    let installed = dir_cache.join("arabidopsis_thaliana/63_TAIR10");
    let fasta = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz";
    for file in [
        "info.json".to_string(),
        "chr_synonyms.txt".to_string(),
        "CHECKSUMS.sha256".to_string(),
        fasta.to_string(),
        format!("{fasta}.fai"),
        format!("{fasta}.gzi"),
    ] {
        assert!(installed.join(&file).is_file(), "{file}");
    }
    assert!(installed.join("transcripts/1").is_dir());
    assert!(
        !dir_cache.join("tmp").exists(),
        "the work directory is removed"
    );
    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(installed.join("info.json")).unwrap())
            .unwrap();
    assert_eq!(info["species"], PLANT);
    assert_eq!(info["assembly"], "TAIR10");
    assert_eq!(info["cache_version"], 63);
    assert_eq!(info["source_versions"]["assembly"], "TAIR10", "{info}");
    assert_eq!(
        std::fs::read_to_string(installed.join("chr_synonyms.txt")).unwrap(),
        "1\tChr1\n"
    );

    let transcripts = transcripts(&installed);
    assert_eq!(transcripts.len(), 7);
    let first = &transcripts["AT1G01010.1"];
    assert_eq!(first["gene_symbol"], "NAC001", "{first}");
    assert_eq!(first["gene_symbol_source"], "TAIR_SYMBOL", "{first}");
    assert_eq!(first["is_canonical"], 1, "{first}");
    assert_eq!(transcripts["AT1G01020.1"]["is_canonical"], 1);
    assert!(transcripts["AT1G01020.2"]["is_canonical"].is_null());

    let position = 4_000;
    let reference = base(position - 1) as char;
    let alternate = if reference == 'A' { 'C' } else { 'A' };
    let variants = root.path().join("variants.vcf");
    std::fs::write(
        &variants,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t{position}\tmade_up\t{reference}\t{alternate}\t.\t.\t.\n"
        ),
    )
    .unwrap();
    let dir_cache_s = dir_cache.display().to_string();
    let rows = annotate(
        &[
            "--dir_cache",
            &dir_cache_s,
            "--species",
            PLANT,
            "--cache_version",
            "63",
            "--assembly",
            "TAIR10",
            "--hgvs",
        ],
        &variants,
        &root.path().join("annotated.txt"),
    );
    assert!(
        rows.iter()
            .any(|row| row.contains("\tAT1G01010.1\tTranscript\t")),
        "{rows:?}"
    );
}

/// A table dump whose BSD sum differs from its `CHECKSUMS` line stops the
/// installation with one error line naming the file and both checksums: no
/// species directory exists, the file is not kept, and the inputs verified
/// before it stay in place for the next run.
#[test]
fn checksum_mismatch_in_the_division_tree_is_refused_naming_the_file() {
    let site = plant_site(&slice_gff3());
    let core = format!("{PLANTS}/mysql/{PLANT_DATABASE}");
    let listing = site.path().join(&core).join("CHECKSUMS");
    let kept: String = std::fs::read_to_string(&listing)
        .unwrap()
        .lines()
        .filter(|line| !line.ends_with(" gene.txt.gz"))
        .map(|line| format!("{line}\n"))
        .collect();
    std::fs::write(&listing, kept).unwrap();
    let (_, rows) = plant_tables()
        .into_iter()
        .find(|(table, _)| *table == "gene")
        .unwrap();
    let dump = gzip(rows.as_bytes());
    let (sum, blocks) = bsd_sum(&dump);
    site.add_with_sum(&core, "gene.txt.gz", &dump, Some((1, 1)));
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, PLANT, "TAIR10", "63", &dir_cache, &[]);
    let error = error_line(&output);
    assert!(
        error.starts_with("ERROR: [cache-builder] Failed to install the release: checksum for ")
            && error.contains(&format!(
                "gene.txt.gz ({sum:05} {blocks}) does not match the CHECKSUMS line for gene.txt.gz (00001 1)"
            )),
        "{error}"
    );
    assert!(!dir_cache.join(PLANT).exists());
    let work = dir_cache.join("tmp/arabidopsis_thaliana_63_TAIR10");
    assert!(work
        .join("inputs/Arabidopsis_thaliana.TAIR10.63.gff3.gz")
        .is_file());
    assert!(!work.join("core_dumps/gene.txt.gz").exists());
}

/// A tree that lacks an input is refused before anything is fetched: without
/// the GTF's directory the run ends on one line naming the genome, the
/// missing listing and the routes that build the cache, and the cache root is
/// never created.
#[test]
fn a_tree_lacking_an_input_is_refused_before_any_download() {
    let site = plant_site(&slice_gff3());
    std::fs::remove_dir_all(site.path().join(format!("{PLANTS}/gtf"))).unwrap();
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, PLANT, "TAIR10", "63", &dir_cache, &[]);
    assert_eq!(
        error_line(&output),
        format!(
            "ERROR: [cache-builder] Failed to install the release: arabidopsis_thaliana TAIR10 cannot be installed from {}/: pub/release-63/plants/gtf/arabidopsis_thaliana/CHECKSUMS is not on the site (is the release published?); convert Ensembl VEP's cache of it instead (docs/cache-conversion.md), or build from local copies of its files with --gff3, --fasta, --genome-fasta and --codon-table",
            site.path().display()
        )
    );
    assert!(!dir_cache.exists());
}

/// A division's tree that publishes no variation VCF directory for the
/// species installs without the known variants, under `--include-variations`
/// too: the run ends on one warning a part the tree publishes no file of,
/// the first naming the missing listing, and the cache carries no variation
/// file.
#[test]
fn a_division_tree_without_variation_vcfs_installs_without_known_variants() {
    let site = plant_site(&slice_gff3());
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(
        &site,
        PLANT,
        "TAIR10",
        "63",
        &dir_cache,
        &["--include-variations"],
    );
    let log = stderr(&output);
    assert!(output.status.success(), "{log}");
    let lines: Vec<&str> = log.lines().collect();
    let closing: Vec<&str> = lines[lines.len() - 3..]
        .iter()
        .filter_map(|line| line.split_once(" WARN "))
        .map(|(_, text)| text)
        .collect();
    assert_eq!(
        closing,
        [
            "vep_cache_builder::install: No known variants in the cache of arabidopsis_thaliana TAIR10: pub/release-63/plants/variation/vcf/arabidopsis_thaliana/CHECKSUMS is not on the site",
            "vep_cache_builder::install: No SIFT or PolyPhen predictions in the cache of arabidopsis_thaliana TAIR10: no variation database is published for this release and assembly",
            "vep_cache_builder::install: No regulatory features or motif features in the cache of arabidopsis_thaliana TAIR10: pub/release-63/plants/regulation/arabidopsis_thaliana/TAIR10/annotation/CHECKSUMS is not on the site",
        ],
        "{log}"
    );
    let installed = dir_cache.join("arabidopsis_thaliana/63_TAIR10");
    assert!(installed.join("transcripts/1").is_dir());
    assert!(!installed.join("variations").exists());
}

/// A BGZF VCF holding one made-up record at each `(sequence, position)` of
/// `records`, in the layout of Ensembl's variation files, with no index.
fn variation_vcf(records: &[(&str, u32)]) -> Vec<u8> {
    let mut vcf = bgzf::io::Writer::new(Vec::new());
    write!(
        vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    .unwrap();
    for (sequence, position) in records {
        writeln!(
            vcf,
            "{sequence}\t{position}\tmade_up_{position}\tA\tG\t.\t.\tTSA=SNV"
        )
        .unwrap();
    }
    vcf.finish().unwrap()
}

/// With `--include-variations`, a species of an Ensembl Genomes division is
/// installed from the one variation VCF its tree publishes for it
/// (`variation/vcf/<species>/<species>.vcf.gz`, with no per-chromosome file
/// and no `--chromosomes`): the installed cache carries an `all_vars.gz` for
/// each sequence the file holds records on and for no other.
#[test]
fn an_ensembl_genomes_species_reads_the_one_variation_file_its_tree_publishes() {
    let site = plant_site(&slice_gff3());
    site.add(
        &format!("{PLANTS}/variation/vcf/{PLANT}"),
        "arabidopsis_thaliana.vcf.gz",
        &variation_vcf(&[("1", 5_000), ("Pt", 200)]),
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(
        &site,
        PLANT,
        "TAIR10",
        "63",
        &dir_cache,
        &["--include-variations"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let installed = dir_cache.join("arabidopsis_thaliana/63_TAIR10");
    let mut written: Vec<String> = std::fs::read_dir(installed.join("variations"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    written.sort();
    assert_eq!(written, ["1", "Pt"]);
    for sequence in written {
        let all_vars = installed.join(format!("variations/{sequence}/all_vars.gz"));
        assert!(all_vars.is_file(), "{sequence}");
        assert!(
            installed
                .join(format!("variations/{sequence}/all_vars.gz.csi"))
                .is_file(),
            "{sequence}"
        );
    }
    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(installed.join("info.json")).unwrap())
            .unwrap();
    assert!(info["variation_cols"].is_array(), "{info}");
    assert!(!dir_cache.join("tmp").exists());
}

/// A transcript on a sequence region named like an organelle genome whose
/// codon table the core tables do not state stops the installation as it
/// stops a build from local files, naming `--codon-table`; with the flag the
/// same command installs, reusing the inputs it had verified, and the
/// plastid's transcript carries the table named.
#[test]
fn a_region_without_a_codon_table_takes_it_from_the_flag() {
    let mut gff3 = slice_gff3();
    gff3.extend_from_slice(PLASTID_GENE.as_bytes());
    let site = plant_site(&gff3);
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, PLANT, "TAIR10", "63", &dir_cache, &[]);
    assert_eq!(
        error_line(&output),
        format!(
            "ERROR: [cache-builder] Failed to install the release: codon table of sequence region Pt unknown: the core dumps carry no codon_table attribute for it; pass --codon-table Pt=<n>; the work directory {} is kept for the next run",
            dir_cache.join(format!("tmp/{PLANT}_63_TAIR10")).display()
        )
    );
    assert!(!dir_cache.join(PLANT).exists());

    let output = auto_install(
        &site,
        PLANT,
        "TAIR10",
        "63",
        &dir_cache,
        &["--codon-table", "Pt=11"],
    );
    let log = stderr(&output);
    assert!(output.status.success(), "{log}");
    assert!(
        log.contains("Verified ") && !log.contains("Copying "),
        "{log}"
    );
    let transcripts = transcripts(&dir_cache.join("arabidopsis_thaliana/63_TAIR10"));
    assert_eq!(transcripts.len(), 8);
    assert_eq!(
        transcripts["MADEUPT1"]["variation_effect_feature_cache"]["codon_table"],
        11
    );
}

/// In the release tree a species' core database is the one its species list
/// names: `rattus_norvegicus_core_116_1` for `GRCr8`, an assembly whose name
/// ends in another number. The chromosome 21 fixture installs under that
/// species' names from that directory and no other.
#[test]
fn a_database_numbered_apart_from_its_assembly_is_found() {
    let species = "rattus_norvegicus";
    let database = "rattus_norvegicus_core_116_1";
    let release = "pub/release-116";
    let fixture = fixture("grch38_chr21");
    let site = Site::new();
    site.add(
        release,
        "species_EnsemblVertebrates.txt",
        format!(
            "{SPECIES_LIST_HEADER}{}",
            species_row(species, "EnsemblVertebrates", database)
        )
        .as_bytes(),
    );
    let add_gzipped = |dir: String, name: &str, source: &str| {
        site.add(
            &dir,
            name,
            &gzip(&std::fs::read(fixture.join(source)).unwrap()),
        );
    };
    add_gzipped(
        format!("{release}/gff3/{species}"),
        "Rattus_norvegicus.GRCr8.116.gff3.gz",
        "slice.gff3",
    );
    add_gzipped(
        format!("{release}/gtf/{species}"),
        "Rattus_norvegicus.GRCr8.116.gtf.gz",
        "slice.gtf",
    );
    add_gzipped(
        format!("{release}/fasta/{species}/pep"),
        "Rattus_norvegicus.GRCr8.pep.all.fa.gz",
        "peptides.fa",
    );
    for table in plant_tables().iter().map(|(table, _)| *table) {
        add_gzipped(
            format!("{release}/mysql/{database}"),
            &format!("{table}.txt.gz"),
            &format!("{table}.txt"),
        );
    }
    for suffix in ["", ".fai", ".gzi"] {
        site.add(
            &format!("{release}/fasta/{species}/dna_index"),
            &format!("Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz{suffix}"),
            &std::fs::read(golden(&format!("GRCh37-hgvs/reference.fa.gz{suffix}"))).unwrap(),
        );
    }
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, species, "GRCr8", "116", &dir_cache, &[]);
    let log = stderr(&output);
    assert!(output.status.success(), "{log}");
    assert!(
        log.contains(&format!("{release}/mysql/{database}")),
        "{log}"
    );
    let installed = dir_cache.join("rattus_norvegicus/116_GRCr8");
    assert!(installed.join("transcripts/21").is_dir());
    assert!(installed.join("chr_synonyms.txt").is_file());
}

/// A genome of a collection is refused before anything is fetched: its core
/// database holds every genome of the collection. The one error line names
/// the genome, that database and the routes that build the cache, and the
/// cache root is never created.
#[test]
fn a_collection_genome_is_refused_naming_its_shared_database() {
    let species = "bacterium_fictum_gca_000000001";
    let site = Site::new();
    site.add(
        "pub/release-63",
        "species.txt",
        format!(
            "{SPECIES_LIST_HEADER}{}",
            species_row(
                species,
                "EnsemblBacteria",
                "bacteria_7_collection_core_63_116_1"
            )
        )
        .as_bytes(),
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, species, "ASM1v1", "63", &dir_cache, &[]);
    assert_eq!(
        error_line(&output),
        format!(
            "ERROR: [cache-builder] Failed to install the release: bacterium_fictum_gca_000000001 ASM1v1 cannot be installed from {}/: it has no core database of its own (bacteria_7_collection_core_63_116_1 holds every genome of bacteria_7_collection, and a shared database is not split by genome); convert Ensembl VEP's cache of it instead (docs/cache-conversion.md), or build from local copies of its files with --gff3, --fasta, --genome-fasta and --codon-table",
            site.path().display()
        )
    );
    assert!(!dir_cache.exists());
}

/// An assembly of the rapid release, named for its accession as Ensembl VEP's
/// cache of it is, is refused before anything is fetched: its tree publishes
/// no table dumps. The one error line names the assembly, the tree and the
/// routes that build the cache, and the cache root is never created.
#[test]
fn a_rapid_release_assembly_is_refused_naming_its_tree() {
    let species = "homo_sapiens_gca000000001v1";
    let tree = "pub/rapid-release/species/Homo_sapiens/GCA_000000001.1/ensembl";
    let site = Site::new();
    std::fs::create_dir_all(site.path().join(tree).join("geneset/2022_07")).unwrap();
    std::fs::create_dir_all(site.path().join(tree).join("genome")).unwrap();
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, species, "Made-up1.0", "107", &dir_cache, &[]);
    assert_eq!(
        error_line(&output),
        format!(
            "ERROR: [cache-builder] Failed to install the release: homo_sapiens_gca000000001v1 Made-up1.0 cannot be installed from {}/: {tree}/ publishes no core table dumps; convert Ensembl VEP's cache of it instead (docs/cache-conversion.md), or build from local copies of its files with --gff3, --fasta, --genome-fasta and --codon-table",
            site.path().display()
        )
    );
    assert!(!dir_cache.exists());
}
