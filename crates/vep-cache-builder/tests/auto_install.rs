// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--auto <species> <assembly>` installs a release from a directory laid out
//! as the Ensembl site, every file checked against its directory's `CHECKSUMS`:
//! the cache lands at `<dir_cache>/<species>/<release>_<assembly>` with the
//! indexed genome FASTA inside it, and `vep` resolves both through
//! `--dir_cache`. The release's inputs are the chromosome 21 fixture of
//! `tests/fixtures/grch38_chr21/` under the names the site gives them, and the
//! genome FASTA is the GRCh37 chromosome 21 reference committed with the
//! `GRCh37-hgvs` golden corpus.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles::bgzf;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
}

fn golden(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{relative}"))
}

/// The core and variation table dumps of the fixture, every `.txt` file but
/// the expected-transcripts JSON; the two phenotype tables are the variation
/// database's.
const CORE_TABLES: [&str; 15] = [
    "analysis_description",
    "attrib_type",
    "external_db",
    "gene",
    "meta",
    "object_xref",
    "protein_feature",
    "seq_region",
    "seq_region_attrib",
    "seq_region_synonym",
    "transcript",
    "transcript_attrib",
    "translation",
    "translation_attrib",
    "xref",
];
const VARIATION_TABLES: [&str; 2] = ["phenotype", "phenotype_feature"];

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

/// The species as Ensembl spells it in a GFF3, GTF or FASTA file name
/// (`homo_sapiens` is `Homo_sapiens`).
fn capitalized(species: &str) -> String {
    format!("{}{}", species[..1].to_uppercase(), &species[1..])
}

/// A BGZF VCF holding one dbSNP record, at `pos` on `chrom`, in the layout of
/// Ensembl's per-chromosome variation files.
fn variation_vcf(chrom: &str, pos: u32) -> Vec<u8> {
    let mut vcf = bgzf::io::Writer::new(Vec::new());
    write!(
        vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=dbSNP_156,Number=0,Type=Flag,Description=\"dbSNP\">\n\
         ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         {chrom}\t{pos}\trs{pos}\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n"
    )
    .unwrap();
    vcf.finish().unwrap()
}

/// A directory laid out as the Ensembl site, each directory carrying the
/// `CHECKSUMS` its files sum to.
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

    fn root(&self) -> &str {
        self.path().to_str().unwrap()
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

    fn add_copy(&self, dir: &str, name: &str, source: &Path) {
        self.add(dir, name, &std::fs::read(source).unwrap());
    }

    fn add_gzipped(&self, dir: &str, name: &str, source: &Path) {
        self.add(dir, name, &gzip(&std::fs::read(source).unwrap()));
    }

    /// The release directory of `assembly`: Ensembl keeps GRCh37 in a tree of
    /// its own, with the release-87 gene set and uncompressed variation tables.
    fn release_dir(assembly: &str) -> &'static str {
        if assembly == "GRCh37" {
            "pub/grch37/release-116"
        } else {
            "pub/release-116"
        }
    }

    /// The MySQL directory of one database of release 116:
    /// `mysql/<species>_<kind>_116_<N>` under the release directory, `N` the
    /// assembly's trailing digits (`GRCh38` is 38, `GRCz11` 11).
    fn mysql_dir(species: &str, assembly: &str, kind: &str) -> String {
        let number = &assembly[assembly
            .trim_end_matches(|c: char| c.is_ascii_digit())
            .len()..];
        format!(
            "{}/mysql/{species}_{kind}_116_{number}",
            Self::release_dir(assembly)
        )
    }

    /// Places `gff3` and the fixture's GTF, peptides and core table dumps under
    /// the names the site gives release 116 of `species` and `assembly`.
    fn with_core_inputs(&self, species: &str, assembly: &str, gff3: &[u8]) -> &Self {
        let fixture = fixture_dir();
        let root = Self::release_dir(assembly);
        let gene_set = if assembly == "GRCh37" { 87 } else { 116 };
        let prefix = capitalized(species);
        self.add(
            &format!("{root}/gff3/{species}"),
            &format!("{prefix}.{assembly}.{gene_set}.gff3.gz"),
            &gzip(gff3),
        );
        self.add_gzipped(
            &format!("{root}/gtf/{species}"),
            &format!("{prefix}.{assembly}.{gene_set}.gtf.gz"),
            &fixture.join("slice.gtf"),
        );
        self.add_gzipped(
            &format!("{root}/fasta/{species}/pep"),
            &format!("{prefix}.{assembly}.pep.all.fa.gz"),
            &fixture.join("peptides.fa"),
        );
        let core = Self::mysql_dir(species, assembly, "core");
        for table in CORE_TABLES {
            self.add_gzipped(
                &core,
                &format!("{table}.txt.gz"),
                &fixture.join(format!("{table}.txt")),
            );
        }
        self
    }

    /// Places the fixture's inputs as `homo_sapiens` release 116 of `assembly`,
    /// the two phenotype tables as the variation database's.
    fn with_fixture_inputs(&self, assembly: &str) -> &Self {
        let fixture = fixture_dir();
        let gff3 = std::fs::read(fixture.join("slice.gff3")).unwrap();
        self.with_core_inputs("homo_sapiens", assembly, &gff3);
        let variation = Self::mysql_dir("homo_sapiens", assembly, "variation");
        for table in VARIATION_TABLES {
            if assembly == "GRCh37" {
                self.add_copy(
                    &variation,
                    &format!("{table}.txt"),
                    &fixture.join(format!("{table}.txt")),
                );
            } else {
                self.add_gzipped(
                    &variation,
                    &format!("{table}.txt.gz"),
                    &fixture.join(format!("{table}.txt")),
                );
            }
        }
        self
    }

    /// Places the committed GRCh37 chromosome 21 reference with its indexes as
    /// the `dna_index/` toplevel file of `species` and `assembly`.
    fn with_indexed_genome(&self, species: &str, assembly: &str) -> &Self {
        let dir = format!("{}/fasta/{species}/dna_index", Self::release_dir(assembly));
        let prefix = capitalized(species);
        for suffix in ["", ".fai", ".gzi"] {
            self.add_copy(
                &dir,
                &format!("{prefix}.{assembly}.dna.toplevel.fa.gz{suffix}"),
                &golden(&format!("GRCh37-hgvs/reference.fa.gz{suffix}")),
            );
        }
        self
    }
}

fn builder() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
}

/// Runs `--auto <species> <assembly> --release 116` against `site` into
/// `dir_cache`, with `extra` flags.
fn auto_install(
    site: &Site,
    species: &str,
    assembly: &str,
    dir_cache: &Path,
    extra: &[&str],
) -> Output {
    builder()
        .args(["--auto", species, assembly, "--release", "116"])
        .arg("--dir_cache")
        .arg(dir_cache)
        .args(["--site", site.root()])
        .args(extra)
        .output()
        .expect("vep-cache-builder runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Every regular file under `dir`, relative to it, sorted.
fn files_under(dir: &Path) -> Vec<PathBuf> {
    fn walk(root: &Path, dir: &Path, out: &mut Vec<PathBuf>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                walk(root, &path, out);
            } else {
                out.push(path.strip_prefix(root).unwrap().to_path_buf());
            }
        }
    }
    let mut out = Vec::new();
    walk(dir, dir, &mut out);
    out.sort();
    out
}

/// Annotates `input` in-process with `flags`, as the `vep` binary would, and
/// returns the output's data lines (every line not opening with `##`).
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
        .filter(|line| !line.starts_with("##"))
        .map(str::to_string)
        .collect()
}

/// The installed GRCh38 release sits at `<dir_cache>/homo_sapiens/116_GRCh38`
/// with the `dna_index/` FASTA and its two indexes copied byte for byte beside
/// `info.json`, `chr_synonyms.txt` and the shards, all listed in the checksum
/// manifest; the work directory is gone; the shards, `info.json` and
/// `chr_synonyms.txt` are byte for byte those of a build from the same inputs
/// named one flag at a time; and `vep` resolves the installed cache through
/// `--dir_cache --species --cache_version --assembly`, detects the FASTA in it,
/// and writes the same rows for the GRCh38 golden corpus's records as it writes
/// over the flag-named build and `--fasta`.
#[test]
fn auto_install_builds_the_layout_and_matches_the_explicit_build() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, "homo_sapiens", "GRCh38", &dir_cache, &[]);
    assert!(output.status.success(), "{}", stderr(&output));

    let installed = dir_cache.join("homo_sapiens/116_GRCh38");
    let fasta = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz";
    for file in ["info.json", "chr_synonyms.txt", "CHECKSUMS.sha256"] {
        assert!(installed.join(file).is_file(), "{file}");
    }
    assert!(installed.join("transcripts/21").is_dir());
    for suffix in ["", ".fai", ".gzi"] {
        assert_eq!(
            std::fs::read(installed.join(format!("{fasta}{suffix}"))).unwrap(),
            std::fs::read(golden(&format!("GRCh37-hgvs/reference.fa.gz{suffix}"))).unwrap(),
            "{fasta}{suffix}"
        );
    }
    let manifest = std::fs::read_to_string(installed.join("CHECKSUMS.sha256")).unwrap();
    for suffix in ["", ".fai", ".gzi"] {
        assert!(
            manifest.contains(&format!("  {fasta}{suffix}\n")),
            "the manifest lists {fasta}{suffix}"
        );
    }
    assert!(
        !dir_cache.join("tmp").exists(),
        "the work directory is removed"
    );
    let release_files: Vec<PathBuf> = files_under(&installed)
        .into_iter()
        .map(|p| Path::new("116_GRCh38").join(p))
        .collect();
    assert_eq!(
        files_under(&dir_cache.join("homo_sapiens")),
        release_files,
        "nothing else sits under the species directory"
    );
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(installed.join("info.json")).unwrap())
            .unwrap();
    assert_eq!(info["species"], "homo_sapiens");
    assert_eq!(info["assembly"], "GRCh38");
    assert_eq!(info["cache_version"], 116);

    let explicit = tempfile::tempdir().unwrap();
    let fixture = fixture_dir();
    let output = builder()
        .args(["--release", "116", "--assembly", "GRCh38"])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--gtf")
        .arg(fixture.join("slice.gtf"))
        .arg("--core-dumps")
        .arg(&fixture)
        .arg("--genome-fasta")
        .arg(installed.join(fasta))
        .arg("--output-dir")
        .arg(explicit.path())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let shards: Vec<PathBuf> = files_under(explicit.path())
        .into_iter()
        .filter(|p| p.starts_with("transcripts"))
        .collect();
    assert!(!shards.is_empty());
    for relative in shards.iter().chain([
        &PathBuf::from("info.json"),
        &PathBuf::from("chr_synonyms.txt"),
    ]) {
        assert_eq!(
            std::fs::read(installed.join(relative)).unwrap(),
            std::fs::read(explicit.path().join(relative)).unwrap(),
            "{}",
            relative.display()
        );
    }

    let variants = golden("GRCh38/variants.vcf");
    let dir_cache_s = dir_cache.display().to_string();
    let through_layout = annotate(
        &[
            "--dir_cache",
            &dir_cache_s,
            "--species",
            "homo_sapiens",
            "--cache_version",
            "116",
            "--assembly",
            "GRCh38",
        ],
        &variants,
        &root.path().join("layout.txt"),
    );
    let explicit_s = explicit.path().display().to_string();
    let fasta_s = installed.join(fasta).display().to_string();
    let through_json_cache = annotate(
        &[
            "--json_cache",
            &explicit_s,
            "--assembly",
            "GRCh38",
            "--fasta",
            &fasta_s,
        ],
        &variants,
        &root.path().join("explicit.txt"),
    );
    // The corpus's chromosome 21 records: one row per allele, more for the
    // records on the fixture's transcripts; its records on other chromosomes
    // are skipped, the fixture cache holding chromosome 21 alone.
    assert!(through_layout.len() > 300, "{}", through_layout.len());
    assert!(
        through_layout
            .iter()
            .any(|row| row.contains("\tTranscript\t") && !row.contains("intergenic_variant")),
        "no row on a fixture transcript"
    );
    assert_eq!(through_layout, through_json_cache);
}

/// On GRCh37 the inputs come from the `grch37/` tree under the release-87 gene
/// set names, the variation tables uncompressed, and the genome FASTA from
/// `dna/`, there being no `dna_index/`: a source that is already BGZF is kept
/// as it is and indexed, the `.fai` and `.gzi` written for it byte for byte the
/// committed `samtools faidx` indexes of that file. `vep --hgvs` then runs over
/// the installed release with no `--fasta`, reading the FASTA it finds in the
/// cache directory.
#[test]
fn grch37_dna_source_is_kept_as_bgzf_and_indexed() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh37");
    let fasta = "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz";
    site.add_copy(
        "pub/grch37/release-116/fasta/homo_sapiens/dna",
        fasta,
        &golden("GRCh37-hgvs/reference.fa.gz"),
    );
    site.add(
        "pub/grch37/release-116/fasta/homo_sapiens/dna",
        "Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz",
        b"never chosen",
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, "homo_sapiens", "GRCh37", &dir_cache, &[]);
    assert!(output.status.success(), "{}", stderr(&output));

    let installed = dir_cache.join("homo_sapiens/116_GRCh37");
    assert_eq!(
        std::fs::read(installed.join(fasta)).unwrap(),
        std::fs::read(golden("GRCh37-hgvs/reference.fa.gz")).unwrap()
    );
    for suffix in ["fai", "gzi"] {
        assert_eq!(
            std::fs::read(installed.join(format!("{fasta}.{suffix}"))).unwrap(),
            std::fs::read(golden(&format!("GRCh37-hgvs/reference.fa.gz.{suffix}"))).unwrap(),
            "{suffix}"
        );
    }
    assert!(installed.join("transcripts/21").is_dir());
    assert!(!dir_cache.join("tmp").exists());

    let dir_cache_s = dir_cache.display().to_string();
    let rows = annotate(
        &[
            "--dir_cache",
            &dir_cache_s,
            "--species",
            "homo_sapiens",
            "--cache_version",
            "116",
            "--assembly",
            "GRCh37",
            "--hgvs",
        ],
        &golden("GRCh37-hgvs/variants.vcf"),
        &root.path().join("hgvs.txt"),
    );
    assert!(rows.len() > 293, "{}", rows.len());
}

/// With `--include-variations`, the chromosome's variation VCF (here without
/// a published index) and the variation database's ClinVar tables are fetched
/// and verified like every other input, and the installed cache carries the
/// chromosome's `all_vars.gz` with `variation_cols` in `info.json`. The two
/// ClinVar tables the gene phenotype count also reads (`phenotype`,
/// `phenotype_feature`) are read from the site once.
#[test]
fn auto_install_fetches_the_variation_files() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    site.add(
        "pub/release-116/variation/vcf/homo_sapiens",
        "homo_sapiens-chr21.vcf.gz",
        &variation_vcf("21", 18_000_300),
    );
    let variation = "pub/release-116/mysql/homo_sapiens_variation_116_38";
    for (table, rows) in [
        (
            "source",
            "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n",
        ),
        ("seq_region", "131537\t21\t4\n"),
        (
            "attrib_type",
            "480\tclinvar_clin_sig\tClinVar clinical significance\t\\N\n",
        ),
        ("phenotype_feature_attrib", ""),
    ] {
        site.add(
            variation,
            &format!("{table}.txt.gz"),
            &gzip(rows.as_bytes()),
        );
    }
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(
        &site,
        "homo_sapiens",
        "GRCh38",
        &dir_cache,
        &["--include-variations", "--chromosomes", "21"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    for table in ["phenotype", "phenotype_feature"] {
        let site_file = site
            .path()
            .join(variation)
            .join(format!("{table}.txt.gz"))
            .display()
            .to_string();
        assert_eq!(
            stderr(&output)
                .lines()
                .filter(|line| line.contains("Copying ") && line.contains(&site_file))
                .count(),
            1,
            "{table}.txt.gz is copied from the site once"
        );
    }
    let installed = dir_cache.join("homo_sapiens/116_GRCh38");
    assert!(installed.join("variations/21/all_vars.gz").is_file());
    assert!(installed.join("variations/21/all_vars.gz.csi").is_file());
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(installed.join("info.json")).unwrap())
            .unwrap();
    assert!(info["variation_cols"].is_array(), "{info}");
    assert!(!dir_cache.join("tmp").exists());
}

/// Without `--chromosomes`, a species other than `homo_sapiens` takes the
/// chromosomes whose variation VCFs are fetched from the `##sequence-region`
/// pragmas of its GFF3, the rule a build named one flag at a time follows: the
/// installed cache carries an `all_vars.gz` for each declared region and for
/// no other, on a site publishing none of the 25 human files. The release has
/// no variation database, so the ClinVar tables are passed over and every
/// transcript carries `gene_phenotype` 0, the count the Ensembl VEP dumper
/// stores without a phenotype adaptor (`Database/Transcript.pm` lines 202,
/// 243-248 and 284).
#[test]
fn another_species_takes_its_chromosomes_from_the_gff3() {
    let regions = ["21", "MT"];
    let fixture = std::fs::read_to_string(fixture_dir().join("slice.gff3")).unwrap();
    let (version, features) = fixture.split_once('\n').unwrap();
    let mut gff3 = format!("{version}\n");
    for (region, length) in regions.iter().zip([45_935_719, 16_596]) {
        gff3.push_str(&format!("##sequence-region   {region} 1 {length}\n"));
    }
    gff3.push_str(features);
    let site = Site::new();
    site.with_core_inputs("danio_rerio", "GRCz11", gff3.as_bytes())
        .with_indexed_genome("danio_rerio", "GRCz11");
    for (region, pos) in regions.iter().zip([18_000_300, 100]) {
        site.add(
            "pub/release-116/variation/vcf/danio_rerio",
            &format!("danio_rerio-chr{region}.vcf.gz"),
            &variation_vcf(region, pos),
        );
    }
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(
        &site,
        "danio_rerio",
        "GRCz11",
        &dir_cache,
        &["--include-variations"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let installed = dir_cache.join("danio_rerio/116_GRCz11");
    let mut written: Vec<String> = std::fs::read_dir(installed.join("variations"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    written.sort();
    assert_eq!(written, regions);
    for region in regions {
        assert!(
            installed
                .join(format!("variations/{region}/all_vars.gz"))
                .is_file(),
            "{region}"
        );
    }
    let mut transcripts = 0;
    for shard in std::fs::read_dir(installed.join("transcripts/21")).unwrap() {
        let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
        let records: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        for record in records {
            assert_eq!(record["gene_phenotype"], 0, "{}", record["stable_id"]);
            transcripts += 1;
        }
    }
    assert!(transcripts > 0);
    assert!(!dir_cache.join("tmp").exists());
}

/// A file whose BSD sum differs from its `CHECKSUMS` line stops the
/// installation with one error line naming the file and both checksums, before
/// anything is built: no species directory exists, the file is not kept, and
/// the inputs verified before it stay in place for the next run.
#[test]
fn checksum_mismatch_is_refused_naming_the_file() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let gtf_dir = "pub/release-116/gtf/homo_sapiens";
    let gtf = "Homo_sapiens.GRCh38.116.gtf.gz";
    std::fs::remove_file(site.path().join(gtf_dir).join("CHECKSUMS")).unwrap();
    std::fs::remove_file(site.path().join(gtf_dir).join(gtf)).unwrap();
    site.add_with_sum(
        gtf_dir,
        gtf,
        &gzip(&std::fs::read(fixture_dir().join("slice.gtf")).unwrap()),
        Some((1, 1)),
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, "homo_sapiens", "GRCh38", &dir_cache, &[]);
    assert!(!output.status.success());
    let text = stderr(&output);
    let error = text
        .lines()
        .find(|line| line.starts_with("ERROR: [cache-builder] "))
        .unwrap_or_else(|| panic!("no error line in {text}"));
    assert!(
        error.contains(gtf) && error.contains("00001 1") && error.contains("CHECKSUMS"),
        "{error}"
    );
    assert!(!dir_cache.join("homo_sapiens").exists());
    let kept = dir_cache.join("tmp/homo_sapiens_116_GRCh38/inputs");
    assert!(kept.join("Homo_sapiens.GRCh38.116.gff3.gz").is_file());
    assert!(!kept.join(gtf).exists());
}

/// A release directory that already exists is left untouched: the
/// installation refuses it by name before fetching anything.
#[test]
fn an_existing_release_directory_is_never_written_into() {
    let site = Site::new();
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");
    let existing = dir_cache.join("homo_sapiens/116_GRCh38");
    std::fs::create_dir_all(&existing).unwrap();
    std::fs::write(existing.join("info.json"), "{}").unwrap();

    let output = auto_install(&site, "homo_sapiens", "GRCh38", &dir_cache, &[]);
    assert!(!output.status.success());
    let text = stderr(&output);
    assert!(
        text.contains(&format!(
            "ERROR: [cache-builder] Failed to install the release: {} exists",
            existing.display()
        )),
        "{text}"
    );
    assert_eq!(
        files_under(&dir_cache),
        [PathBuf::from("homo_sapiens/116_GRCh38/info.json")]
    );
    assert_eq!(
        std::fs::read_to_string(existing.join("info.json")).unwrap(),
        "{}"
    );
}
