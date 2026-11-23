// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--auto <species> <assembly>` installs a release from a directory laid out
//! as the Ensembl site, every file checked against its directory's `CHECKSUMS`:
//! the cache lands at `<dir_cache>/<species>/<release>_<assembly>` with the
//! indexed genome FASTA inside it, and `vep` resolves both through
//! `--dir_cache`. The release's inputs are the chromosome 21 fixture of
//! `tests/fixtures/grch38_chr21/` under the names the site gives them, and the
//! genome FASTA is the GRCh37 chromosome 21 reference committed with the
//! `GRCh37-hgvs` golden corpus. The same tree is served by a loopback HTTP
//! server where a test turns on a server's answer for a file it does not
//! hold, 403 or 404.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
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
    variation_vcf_of(&[(chrom, pos)])
}

/// A BGZF VCF holding one dbSNP record at each `(chrom, pos)` of `records`,
/// in the layout of Ensembl's variation files.
fn variation_vcf_of(records: &[(&str, u32)]) -> Vec<u8> {
    let mut vcf = bgzf::io::Writer::new(Vec::new());
    write!(
        vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=dbSNP_156,Number=0,Type=Flag,Description=\"dbSNP\">\n\
         ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    .unwrap();
    for (chrom, pos) in records {
        writeln!(
            vcf,
            "{chrom}\t{pos}\trs{pos}\tA\tG\t.\t.\tdbSNP_156;TSA=SNV"
        )
        .unwrap();
    }
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
    auto_install_from(site.root(), species, assembly, dir_cache, extra)
}

/// As [`auto_install`], against the site `--site` names: a directory or the
/// address of a [`Server`], which is reached past any proxy of the
/// environment.
fn auto_install_from(
    site: &str,
    species: &str,
    assembly: &str,
    dir_cache: &Path,
    extra: &[&str],
) -> Output {
    builder()
        .args(["--auto", species, assembly, "--release", "116"])
        .arg("--dir_cache")
        .arg(dir_cache)
        .args(["--site", site])
        .args(extra)
        .env("NO_PROXY", "127.0.0.1")
        .env("no_proxy", "127.0.0.1")
        .output()
        .expect("vep-cache-builder runs")
}

/// A loopback HTTP server over a directory laid out as a site: a request for
/// a file of the directory is answered 200 with its bytes, and a request for
/// anything else (a file the directory lacks, a directory's index page) with
/// the status line `absent`.
struct Server {
    address: String,
}

impl Server {
    fn over(root: &Path, absent: &'static str) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = format!("http://{}", listener.local_addr().unwrap());
        let root = root.to_path_buf();
        std::thread::spawn(move || {
            for stream in listener.incoming().flatten() {
                let root = root.clone();
                std::thread::spawn(move || Self::answer(stream, &root, absent));
            }
        });
        Self { address }
    }

    /// Reads one request and answers it, closing the connection.
    fn answer(mut stream: TcpStream, root: &Path, absent: &str) {
        let mut request = BufReader::new(stream.try_clone().unwrap());
        let mut first = String::new();
        if request.read_line(&mut first).is_err() {
            return;
        }
        let mut header = String::new();
        while matches!(request.read_line(&mut header), Ok(length) if length > 2) {
            header.clear();
        }
        let path = first.split_whitespace().nth(1).unwrap_or("/");
        let file = root.join(path.trim_start_matches('/'));
        let (status, body) = if file.is_file() {
            ("200 OK", std::fs::read(&file).unwrap())
        } else {
            (absent, Vec::new())
        };
        let _ = write!(
            stream,
            "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        );
        let _ = stream.write_all(&body);
    }
}

/// The one `ERROR:` line of a failed run.
fn error_line(output: &Output) -> String {
    assert_eq!(output.status.code(), Some(1), "{}", stderr(output));
    let text = stderr(output);
    let lines: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("ERROR: "))
        .collect();
    assert_eq!(lines.len(), 1, "{text}");
    lines[0].to_string()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The lines a run logged at the warning level, each from its module's name on.
fn warnings(output: &Output) -> Vec<String> {
    stderr(output)
        .lines()
        .filter_map(|line| line.split_once(" WARN "))
        .map(|(_, text)| text.to_string())
        .collect()
}

/// The lines a run logged after the line naming the installed directory, each
/// from its level on.
fn closing_lines(output: &Output) -> Vec<String> {
    let log = stderr(output);
    let logged: Vec<&str> = log
        .lines()
        .filter_map(|line| line.split_once(' '))
        .map(|(_, logged)| logged.trim_start())
        .collect();
    let installed = logged
        .iter()
        .position(|line| line.starts_with("INFO vep_cache_builder::install: Installed "))
        .unwrap_or_else(|| panic!("no line names the installed directory: {log}"));
    logged[installed + 1..]
        .iter()
        .map(|line| line.to_string())
        .collect()
}

/// Where the BGZF block that opens `file` ends, in the file and in the data it
/// inflates to: the `BC` subfield's value plus one, and the four bytes that
/// close the block.
fn first_block_end(file: &[u8]) -> (u64, u64) {
    let length = usize::from(u16::from_le_bytes([file[16], file[17]])) + 1;
    let inflated = u32::from_le_bytes(file[length - 4..length].try_into().unwrap());
    (length as u64, u64::from(inflated))
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

/// A `dna_index/` file whose published `.gzi` is the index of another
/// compression of the same sequence is installed all the same. The installer
/// names, in the first warning of the run, the file and the first block the
/// index places elsewhere, and the run logs no other warning but the three it
/// ends on, one a part the tree publishes no file of. The installer keeps the
/// published file byte for byte and writes both indexes from it: the `.fai`
/// is the committed `samtools faidx` index of that sequence, the FASTA opens
/// through the two, and `vep --hgvs` over the installed release writes the
/// rows it writes with the committed reference named by `--fasta`.
#[test]
fn a_published_index_that_does_not_describe_its_file_is_written_by_the_installer() {
    let mut sequence = Vec::new();
    bgzf::io::Reader::new(std::fs::File::open(golden("GRCh37-hgvs/reference.fa.gz")).unwrap())
        .read_to_end(&mut sequence)
        .unwrap();
    let mut writer = bgzf::io::writer::Builder::default()
        .set_compression_level(bgzf::io::writer::CompressionLevel::FAST)
        .build_from_writer(Vec::new());
    writer.write_all(&sequence).unwrap();
    let published = writer.finish().unwrap();
    let committed_fai = std::fs::read(golden("GRCh37-hgvs/reference.fa.gz.fai")).unwrap();
    let stale_gzi = std::fs::read(golden("GRCh37-hgvs/reference.fa.gz.gzi")).unwrap();
    // The first entry of the committed index, after its eight-byte count: where
    // the committed file's first block ends, in the file and in the sequence.
    let entry = |at: usize| u64::from_le_bytes(stale_gzi[at..at + 8].try_into().unwrap());
    let (indexed_end, indexed_sequence) = (entry(8), entry(16));
    let (block_end, block_sequence) = first_block_end(&published);
    assert_ne!(
        (block_end, block_sequence),
        (indexed_end, indexed_sequence),
        "the second compression's first block ends where the committed one's does"
    );

    let site = Site::new();
    site.with_fixture_inputs("GRCh38");
    let fasta = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz";
    let dna_index = "pub/release-116/fasta/homo_sapiens/dna_index";
    site.add(dna_index, fasta, &published);
    site.add(dna_index, &format!("{fasta}.fai"), &committed_fai);
    site.add(dna_index, &format!("{fasta}.gzi"), &stale_gzi);
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, "homo_sapiens", "GRCh38", &dir_cache, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    let placed = dir_cache
        .join("tmp/homo_sapiens_116_GRCh38/cache")
        .join(fasta);
    assert_eq!(
        warnings(&output),
        [
            format!(
                "vep_cache_builder::install: The .fai and .gzi published beside {fasta} do not \
                 open it, so both are written here from the file itself: BGZF index \
                 {placed}.gzi does not describe {placed} (the block at byte 0 ends at byte \
                 {block_end} with the sequence at {block_sequence} bytes, the index says \
                 {indexed_end} and {indexed_sequence}): index the file again with samtools faidx",
                placed = placed.display()
            )
            .as_str(),
            "vep_cache_builder::install: No known variants or their frequencies in the cache of \
             homo_sapiens GRCh38: pub/release-116/variation/vcf/homo_sapiens/CHECKSUMS is not \
             on the site",
            "vep_cache_builder::install: No SIFT or PolyPhen predictions in the cache of \
             homo_sapiens GRCh38: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS \
             lists no attrib_type table; --include-predictions fetches them from Ensembl's \
             MySQL server (docs/cache-predictions.md)",
            "vep_cache_builder::install: No regulatory features or motif features in the cache \
             of homo_sapiens GRCh38: \
             pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS is not on the \
             site",
        ]
    );

    let installed = dir_cache.join("homo_sapiens/116_GRCh38");
    assert_eq!(std::fs::read(installed.join(fasta)).unwrap(), published);
    assert_eq!(
        std::fs::read(installed.join(format!("{fasta}.fai"))).unwrap(),
        committed_fai
    );
    assert_ne!(
        std::fs::read(installed.join(format!("{fasta}.gzi"))).unwrap(),
        stale_gzi
    );
    vep_fasta::IndexedFasta::from_path(installed.join(fasta)).unwrap();
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

    let variants = golden("GRCh38/variants.vcf");
    let dir_cache_s = dir_cache.display().to_string();
    let layout = [
        "--dir_cache",
        &dir_cache_s,
        "--species",
        "homo_sapiens",
        "--cache_version",
        "116",
        "--assembly",
        "GRCh38",
        "--hgvs",
    ];
    let through_installed_fasta = annotate(&layout, &variants, &root.path().join("installed.txt"));
    let reference_s = golden("GRCh37-hgvs/reference.fa.gz").display().to_string();
    let through_committed_fasta = annotate(
        &[&layout[..], &["--fasta", &reference_s]].concat(),
        &variants,
        &root.path().join("committed.txt"),
    );
    assert!(
        through_installed_fasta
            .iter()
            .any(|row| row.contains("HGVSc=")),
        "no row carries an HGVSc"
    );
    assert_eq!(through_installed_fasta, through_committed_fasta);
}

/// An unknown `VEP_FASTA_STORAGE` is the environment's refusal, not the
/// published indexes': the installation ends on the one error line that names
/// the value and the work directory, and logs no warning.
#[test]
fn an_unknown_fasta_storage_mode_is_refused_without_a_warning() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = builder()
        .args(["--auto", "homo_sapiens", "GRCh38", "--release", "116"])
        .arg("--dir_cache")
        .arg(&dir_cache)
        .args(["--site", site.root()])
        .env("VEP_FASTA_STORAGE", "neither")
        .output()
        .expect("vep-cache-builder runs");
    let work_dir = dir_cache.join("tmp/homo_sapiens_116_GRCh38");
    assert_eq!(
        error_line(&output),
        format!(
            "ERROR: [cache-builder] Failed to install the release: The indexed FASTA {} does \
             not open: invalid VEP_FASTA_STORAGE=\"neither\": expected 'mmap' or 'memory'; the \
             work directory {} is kept for the next run",
            work_dir
                .join("cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
                .display(),
            work_dir.display()
        )
    );
    assert_eq!(warnings(&output), Vec::<String>::new());
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

/// A species other than `homo_sapiens` is installed from the one genome-wide
/// variation VCF the site publishes for it (`<species>.vcf.gz`), the only
/// form Ensembl publishes for a species but human: with no `--chromosomes`,
/// on a site publishing no per-chromosome file, the installed cache carries
/// an `all_vars.gz` for each sequence the file holds records on and for no
/// other. The release has no variation database, so the ClinVar tables are
/// passed over and every transcript carries `gene_phenotype` 0, the count the
/// Ensembl VEP dumper stores without a phenotype adaptor
/// (`Database/Transcript.pm` lines 202, 243-248 and 284).
#[test]
fn another_species_is_installed_from_its_one_genome_wide_variation_file() {
    let regions = ["21", "MT"];
    let gff3 = std::fs::read(fixture_dir().join("slice.gff3")).unwrap();
    let site = Site::new();
    site.with_core_inputs("danio_rerio", "GRCz11", &gff3)
        .with_indexed_genome("danio_rerio", "GRCz11");
    site.add(
        "pub/release-116/variation/vcf/danio_rerio",
        "danio_rerio.vcf.gz",
        &variation_vcf_of(&[("21", 18_000_300), ("MT", 100)]),
    );
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

/// `--chromosomes` beside `--auto` narrows the read of the one genome-wide
/// variation VCF of a species other than `homo_sapiens` to the sequences it
/// names: no per-chromosome file is looked for, the layout of human alone,
/// and the installed cache carries an `all_vars.gz` for the named sequence
/// and for no other the file holds records on.
#[test]
fn chromosomes_narrow_the_one_genome_wide_file_of_another_species() {
    let gff3 = std::fs::read(fixture_dir().join("slice.gff3")).unwrap();
    let site = Site::new();
    site.with_core_inputs("danio_rerio", "GRCz11", &gff3)
        .with_indexed_genome("danio_rerio", "GRCz11");
    site.add(
        "pub/release-116/variation/vcf/danio_rerio",
        "danio_rerio.vcf.gz",
        &variation_vcf_of(&[("21", 18_000_300), ("MT", 100)]),
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(
        &site,
        "danio_rerio",
        "GRCz11",
        &dir_cache,
        &["--include-variations", "--chromosomes", "21"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let installed = dir_cache.join("danio_rerio/116_GRCz11");
    let written: Vec<String> = std::fs::read_dir(installed.join("variations"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    assert_eq!(written, ["21"]);
    assert!(installed.join("variations/21/all_vars.gz").is_file());
    assert!(!dir_cache.join("tmp").exists());
}

/// A site whose variation directory lists no germline file of the species
/// installs without the known variants, under `--include-variations` too:
/// the run ends on one warning a part the site publishes no file of, the
/// first naming the file looked for, the species' one genome-wide file, and
/// the cache carries no variation file.
#[test]
fn a_site_without_the_species_variation_file_installs_without_known_variants() {
    let gff3 = std::fs::read(fixture_dir().join("slice.gff3")).unwrap();
    let site = Site::new();
    site.with_core_inputs("danio_rerio", "GRCz11", &gff3)
        .with_indexed_genome("danio_rerio", "GRCz11");
    site.add(
        "pub/release-116/variation/vcf/danio_rerio",
        "README",
        b"no variation file\n",
    );
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
    assert_eq!(
        closing_lines(&output),
        [
            "WARN vep_cache_builder::install: No known variants in the cache of danio_rerio \
             GRCz11: pub/release-116/variation/vcf/danio_rerio/CHECKSUMS lists no \
             danio_rerio.vcf.gz",
            "WARN vep_cache_builder::install: No SIFT or PolyPhen predictions in the cache of \
             danio_rerio GRCz11: no variation database is published for this release and \
             assembly",
            "WARN vep_cache_builder::install: No regulatory features or motif features in the \
             cache of danio_rerio GRCz11: \
             pub/release-116/regulation/danio_rerio/GRCz11/annotation/CHECKSUMS is not on the \
             site",
        ]
    );
    let installed = dir_cache.join("danio_rerio/116_GRCz11");
    assert!(installed.join("transcripts/21").is_dir());
    assert!(!installed.join("variations").exists());
    assert!(!dir_cache.join("tmp").exists());
}

/// Under `homo_sapiens` a per-chromosome variation VCF the site does not list
/// is refused before anything is fetched, naming the file and
/// `--chromosomes`, and the cache root is never created.
#[test]
fn a_per_chromosome_file_the_site_lacks_is_refused_before_any_download() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    site.add(
        "pub/release-116/variation/vcf/homo_sapiens",
        "homo_sapiens-chr21.vcf.gz",
        &variation_vcf("21", 18_000_300),
    );
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(
        &site,
        "homo_sapiens",
        "GRCh38",
        &dir_cache,
        &["--include-variations", "--chromosomes", "21,22"],
    );
    assert_eq!(
        error_line(&output),
        "ERROR: [cache-builder] Failed to install the release: homo_sapiens-chr22.vcf.gz is \
         not published for this release and assembly; name the chromosomes with --chromosomes"
    );
    assert!(!dir_cache.exists());
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
    let work_dir = dir_cache.join("tmp/homo_sapiens_116_GRCh38");
    assert!(
        error.ends_with(&format!(
            "; the work directory {} is kept for the next run",
            work_dir.display()
        )),
        "{error}"
    );
    let kept = work_dir.join("inputs");
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

/// A failure that leaves no file in the work directory names no work
/// directory: the first input fetched fails its checksum, so nothing is kept
/// for the next run, the error line ends on the mismatch and the directory
/// is gone with the `tmp` it left empty.
#[test]
fn a_failure_that_keeps_no_file_names_no_work_directory() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let gff3_dir = "pub/release-116/gff3/homo_sapiens";
    let gff3 = "Homo_sapiens.GRCh38.116.gff3.gz";
    let bytes = std::fs::read(site.path().join(gff3_dir).join(gff3)).unwrap();
    std::fs::remove_file(site.path().join(gff3_dir).join("CHECKSUMS")).unwrap();
    site.add_with_sum(gff3_dir, gff3, &bytes, Some((1, 1)));
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install(&site, "homo_sapiens", "GRCh38", &dir_cache, &[]);
    let error = error_line(&output);
    assert!(
        error.starts_with("ERROR: [cache-builder] Failed to install the release: checksum for ")
            && error.ends_with(&format!(
                "does not match the CHECKSUMS line for {gff3} (00001 1)"
            )),
        "{error}"
    );
    assert!(!dir_cache.join("tmp").exists());
    assert!(!dir_cache.join("homo_sapiens").exists());
}

/// Every file under the installed release directory with its bytes.
fn installed_files(installed: &Path) -> Vec<(PathBuf, Vec<u8>)> {
    files_under(installed)
        .into_iter()
        .map(|relative| {
            let bytes = std::fs::read(installed.join(&relative)).unwrap();
            (relative, bytes)
        })
        .collect()
}

/// A server that answers 403 for a file it does not hold, on a tree with no
/// species list and with no index pages, installs as a server that answers
/// 404 does: the release root's `CHECKSUMS` and the `mysql/` index page read
/// as absent, the core database is the one named for the assembly's number,
/// and the two installed caches are the same files byte for byte.
#[test]
fn a_server_answering_403_for_an_absent_file_installs_as_one_answering_404() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let root = tempfile::tempdir().unwrap();
    let mut installed = Vec::new();
    for (answer, absent) in [("403", "403 Forbidden"), ("404", "404 Not Found")] {
        let server = Server::over(site.path(), absent);
        let dir_cache = root.path().join(answer);
        let output = auto_install_from(&server.address, "homo_sapiens", "GRCh38", &dir_cache, &[]);
        assert!(output.status.success(), "{answer}: {}", stderr(&output));
        let release = dir_cache.join("homo_sapiens/116_GRCh38");
        assert!(release.join("info.json").is_file(), "{answer}");
        assert!(release.join("transcripts/21").is_dir(), "{answer}");
        assert!(!dir_cache.join("tmp").exists(), "{answer}");
        installed.push(installed_files(&release));
    }
    assert!(installed[0].len() > 5);
    assert!(installed[0] == installed[1], "the two caches differ");
}

/// On a server that answers 403 for what it does not hold, a directory a
/// tree may lack reads as absent: the GRCh37 tree has no `dna_index/`, so the
/// genome comes from `dna/`, and without the variation database the release
/// installs with the warning that none is published.
#[test]
fn a_403_for_a_directory_a_tree_may_lack_reads_as_absent() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh37");
    let fasta = "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz";
    site.add_copy(
        "pub/grch37/release-116/fasta/homo_sapiens/dna",
        fasta,
        &golden("GRCh37-hgvs/reference.fa.gz"),
    );
    let variation = Site::mysql_dir("homo_sapiens", "GRCh37", "variation");
    std::fs::remove_dir_all(site.path().join(variation)).unwrap();
    let server = Server::over(site.path(), "403 Forbidden");
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install_from(&server.address, "homo_sapiens", "GRCh37", &dir_cache, &[]);
    let log = stderr(&output);
    assert!(output.status.success(), "{log}");
    assert!(
        log.contains("No variation database is published for this release and assembly"),
        "{log}"
    );
    let installed = dir_cache.join("homo_sapiens/116_GRCh37");
    for suffix in ["", ".fai", ".gzi"] {
        assert!(
            installed.join(format!("{fasta}{suffix}")).is_file(),
            "{suffix}"
        );
    }
    assert!(installed.join("transcripts/21").is_dir());
    assert!(!dir_cache.join("tmp").exists());
}

/// On a server with no species list and no index pages, a tree without the
/// core database named for the assembly's number is refused before anything
/// is fetched, with the same error line whether the server answers 403 or
/// 404 for that database's `CHECKSUMS`.
#[test]
fn a_403_for_the_numbered_core_database_is_refused_as_a_404_is() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let core = Site::mysql_dir("homo_sapiens", "GRCh38", "core");
    std::fs::remove_dir_all(site.path().join(core)).unwrap();
    let root = tempfile::tempdir().unwrap();
    for (answer, absent) in [("403", "403 Forbidden"), ("404", "404 Not Found")] {
        let server = Server::over(site.path(), absent);
        let dir_cache = root.path().join(answer);
        let output = auto_install_from(&server.address, "homo_sapiens", "GRCh38", &dir_cache, &[]);
        assert_eq!(
            error_line(&output),
            "ERROR: [cache-builder] Failed to install the release: pub/release-116/mysql/ is \
             not on the site: is the release published?",
            "{answer}"
        );
        assert!(!dir_cache.exists(), "{answer}");
    }
}

/// A 403 for a file the installation cannot do without refuses it before
/// anything is fetched: the one error line names the file, the GTF
/// directory's `CHECKSUMS`, and the server's answer, and the cache root is
/// never created.
#[test]
fn a_403_for_a_listing_the_installation_needs_is_refused_naming_the_answer() {
    let site = Site::new();
    site.with_fixture_inputs("GRCh38")
        .with_indexed_genome("homo_sapiens", "GRCh38");
    let listing = "pub/release-116/gtf/homo_sapiens/CHECKSUMS";
    std::fs::remove_file(site.path().join(listing)).unwrap();
    let server = Server::over(site.path(), "403 Forbidden");
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");

    let output = auto_install_from(&server.address, "homo_sapiens", "GRCh38", &dir_cache, &[]);
    assert_eq!(
        error_line(&output),
        format!(
            "ERROR: [cache-builder] Failed to install the release: Failed to read {listing}: \
             Download failed: {}/{listing}: HTTP status 403 Forbidden",
            server.address
        )
    );
    assert!(!dir_cache.exists());
}
