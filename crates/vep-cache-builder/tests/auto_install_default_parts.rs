// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--auto <species> <assembly>` brings, with no further flag, each part of
//! Ensembl VEP's own cache that a file of the release carries: the known
//! variants with their 1000 Genomes frequencies from the variation VCFs, the
//! SIFT and PolyPhen matrices from the variation database's table dumps and
//! the regulatory and motif features from the regulatory annotation GFF3
//! files. A part the tree publishes no file of is left out. The installation
//! ends, after the line naming the installed directory, on one line saying
//! the cache holds each part or on one warning a part left out, then on one
//! warning for the sequences of the cache it read no variation record on and
//! on one for each reason a cell of a brought part holds no value: the gnomAD
//! frequencies, which no file of the release carries, and each cell whose
//! file or table the tree does not publish, named with the listing that
//! lacks it. A cell the installation brought draws no warning, and the
//! cache's `info.json` names under `left_out` the cells it did not bring.
//!
//! Each tree is a local directory laid out as the Ensembl site, every file
//! with its `CHECKSUMS` line: the chromosome 21 fixture of
//! `tests/fixtures/grch38_chr21/` with the GRCh37 chromosome 21 reference of
//! the `GRCh37-hgvs` golden corpus as its genome, made-up variation VCFs and
//! variation table dumps, and the regulatory GFF files of
//! `crates/vep-cli/tests/fixtures/regulatory-gff/`. A tree may publish beside
//! those GFF files a made-up activity file, and a funcgen database of made-up
//! tables for the matrices the motif fixture names: the cache installed from
//! it carries each regulatory feature's activity by cell type and each
//! matrix's counts, as a conversion of Ensembl VEP's own cache does. A tree
//! may also publish, in its variation database, made-up tables of the
//! publications that cite a variant, of its synonyms and of the failed
//! variants, and beside its variation VCFs the GVF file of its failed
//! variants: the cache installed from it carries the `pubmed`,
//! `var_synonyms` and `failed` cells of its known variants and a record for
//! each failed variant the VCF leaves out.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use base64::Engine;
use noodles::bgzf;
use serde_json::Value;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
}

fn golden(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{relative}"))
}

/// A file of the regulatory GFF fixture.
fn regulatory_gff(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../vep-cli/tests/fixtures/regulatory-gff")
        .join(name)
}

/// The core table dumps of the fixture.
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

/// The species as Ensembl spells it in a file name (`homo_sapiens` is
/// `Homo_sapiens`).
fn capitalized(species: &str) -> String {
    format!("{}{}", species[..1].to_uppercase(), &species[1..])
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

    /// Writes `dir/name` and its `CHECKSUMS` line.
    fn add(&self, dir: &str, name: &str, bytes: &[u8]) {
        let dir = self.path().join(dir);
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

    fn add_copy(&self, dir: &str, name: &str, source: &Path) {
        self.add(dir, name, &std::fs::read(source).unwrap());
    }

    fn add_gzipped(&self, dir: &str, name: &str, source: &Path) {
        self.add(dir, name, &gzip(&std::fs::read(source).unwrap()));
    }
}

/// A genome as release 116 publishes it.
#[derive(Clone, Copy)]
struct Genome {
    species: &'static str,
    assembly: &'static str,
    /// The number the names of its core and variation databases end in.
    database: u32,
    /// The directory of its regulatory annotation: the assembly's name, or
    /// the last dot-separated part of it.
    regulation: &'static str,
}

const HUMAN: Genome = Genome {
    species: "homo_sapiens",
    assembly: "GRCh38",
    database: 38,
    regulation: "GRCh38",
};

const ZEBRAFISH: Genome = Genome {
    species: "danio_rerio",
    assembly: "GRCz11",
    database: 11,
    regulation: "GRCz11",
};

const CHICKEN: Genome = Genome {
    species: "gallus_gallus",
    assembly: "bGalGal1.mat.broiler.GRCg7b",
    database: 7,
    regulation: "GRCg7b",
};

/// What a tree publishes of each part an installation brings.
#[derive(Clone, Copy)]
struct Parts {
    variation: Variation,
    database: Database,
    regulation: Regulation,
}

/// The directory of the variation VCFs.
#[derive(Clone, Copy, PartialEq)]
enum Variation {
    /// The germline file that holds the known variant: chromosome 21's for
    /// `homo_sapiens`, with the 1000 Genomes file, and the one genome-wide
    /// file for any other species.
    Files,
    /// The germline file with its one record on `MT`, a sequence the model
    /// names no transcript on: the file of `MT` for `homo_sapiens`, with the
    /// 1000 Genomes file, and the one genome-wide file for any other species.
    OffTheModel,
    /// The one genome-wide file of `OffTheModel` and, in the GVF directory,
    /// the GVF file of the failed variants.
    OffTheModelAndFailedGvf,
    /// A `README` and no VCF.
    ReadmeAlone,
    /// The files of `Files` for `homo_sapiens`, chromosome 21's holding a
    /// second variant, one the made-up tables and the made-up GVF mark as
    /// failed.
    TwoRecords,
    /// The files of `TwoRecords` and, in the GVF directory, the GVF file of
    /// the failed variants.
    TwoRecordsAndFailedGvf,
    Absent,
}

/// The directory of the variation database.
#[derive(Clone, Copy, PartialEq)]
enum Database {
    /// The tables of the ClinVar columns and the four tables of the
    /// matrices, `attrib` naming SIFT and both PolyPhen analyses.
    EveryTable,
    /// Every table, `attrib` naming SIFT and neither PolyPhen analysis.
    SiftAlone,
    /// Every table but `translation_md5`.
    NoTranslationMd5,
    /// Every table, and the tables the PubMed identifiers, the synonyms and
    /// the failed flag of a known variant are read from.
    VariantTables,
    /// The tables of `VariantTables` but `variation`, the table that names
    /// the variant of each row of the others.
    NoVariationTable,
    /// Every table, and `variation_feature`, the table the known variants
    /// of a sequence are read from where no variation VCF holds a record on
    /// it: `seq_region` names `MT` beside chromosome 21 and `attrib` the
    /// evidence of the phenotype-or-disease flag.
    FeatureTable,
    Absent,
}

/// The directory of the regulatory annotation.
#[derive(Clone, Copy, PartialEq)]
enum Regulation {
    /// The regulatory-feature, motif-feature and EMAR files.
    EveryFile,
    /// The regulatory-feature and EMAR files.
    NoMotifFile,
    /// Every file, the activity file beside them and the funcgen database.
    ActivityAndTables,
    /// Every file and the activity file beside them; no funcgen database.
    ActivityAlone,
    /// Every file and the funcgen database; no activity file.
    TablesAlone,
    Absent,
}

/// The cell types of the made-up activity file, in the order and the
/// spelling of its header.
const CELL_TYPES: [&str; 4] = [
    "tissue two (m, 4 y)",
    "cell line one",
    "tissue one (f, 30 y)",
    "organ three (f, 7 dpf)",
];

/// The cell types as the installed cache lists them: each space an
/// underscore, in byte order.
const LISTED_CELL_TYPES: [&str; 4] = [
    "cell_line_one",
    "organ_three_(f,_7_dpf)",
    "tissue_one_(f,_30_y)",
    "tissue_two_(m,_4_y)",
];

/// The enhancer of the regulatory GFF fixture the record `in_enhancer` lies
/// in, the one feature whose row of the made-up activity file differs from
/// the others'.
const ENHANCER: &str = "ENSR21_85ZF5";

/// One record of a GFF fixture file.
struct GffRecord {
    feature_type: String,
    start: u64,
    end: u64,
    attributes: BTreeMap<String, String>,
}

/// The records of a file of the regulatory GFF fixture, in file order.
fn gff_records(name: &str) -> Vec<GffRecord> {
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(regulatory_gff(name)).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            GffRecord {
                feature_type: columns[2].to_string(),
                start: columns[3].parse().unwrap(),
                end: columns[4].parse().unwrap(),
                attributes: columns[8]
                    .split(';')
                    .filter_map(|attribute| attribute.split_once('='))
                    .map(|(key, value)| (key.to_string(), value.to_string()))
                    .collect(),
            }
        })
        .collect()
}

/// The made-up activity file, gzip-compressed: a row for every record of the
/// regulatory GFF fixture but its CTCF binding sites, as Ensembl's files hold
/// no row of one. [`ENHANCER`] is active in the header's first cell type,
/// inactive in the second, poised in the third and `NA` in the fourth; every
/// other feature is inactive in the first three and repressed in the fourth.
fn activity_file() -> Vec<u8> {
    let mut text = format!("feature_id\tfeature_type\t{}\n", CELL_TYPES.join("\t"));
    for record in gff_records("regulatory_features.gff3.gz") {
        if record.feature_type == "CTCF_binding_site" {
            continue;
        }
        let id = &record.attributes["ID"];
        let activity = if id == ENHANCER {
            ["ACTIVE", "INACTIVE", "POISED", "NA"]
        } else {
            ["INACTIVE", "INACTIVE", "INACTIVE", "REPRESSED"]
        };
        writeln!(
            text,
            "{id}\t{}\t{}",
            record.feature_type,
            activity.join("\t")
        )
        .unwrap();
    }
    gzip(text.as_bytes())
}

/// A motif of the fixture whose span the trees' genome holds as bases, which
/// no other motif overlaps, and whose matrix the made-up tables write to
/// follow the motif's own sequence ([`made_up_counts`]).
struct ScoredMotif {
    id: &'static str,
    matrix: &'static str,
    start: u64,
    end: u64,
    reverse: bool,
}

/// A motif of 15 bases on the forward strand.
const FORWARD_MOTIF: ScoredMotif = ScoredMotif {
    id: "ENSM00000576694",
    matrix: "ENSPFM0001",
    start: 9_568_308,
    end: 9_568_322,
    reverse: false,
};

/// A motif of 13 bases on the reverse strand.
const REVERSE_MOTIF: ScoredMotif = ScoredMotif {
    id: "ENSM00000576738",
    matrix: "ENSPFM0422",
    start: 9_828_884,
    end: 9_828_896,
    reverse: true,
};

/// The four nucleotides in the order of a matrix position's counts.
const NUCLEOTIDES: [u8; 4] = *b"ACGT";

/// The base that pairs with `base`.
fn paired(base: u8) -> u8 {
    NUCLEOTIDES[3 - NUCLEOTIDES.iter().position(|n| *n == base).unwrap()]
}

/// The genome every tree publishes.
fn tree_genome() -> vep_fasta::IndexedFasta {
    vep_fasta::IndexedFasta::from_path(golden("GRCh37-hgvs/reference.fa.gz")).unwrap()
}

impl ScoredMotif {
    /// The motif's bases in the trees' genome as its own strand reads them:
    /// the reverse-strand motif's from its last base to its first, each
    /// paired.
    fn sequence(&self) -> Vec<u8> {
        let bases = tree_genome().sequence("21", self.start, self.end).unwrap();
        if self.reverse {
            bases.iter().rev().map(|base| paired(*base)).collect()
        } else {
            bases
        }
    }

    /// The position on the genome of the motif's 1-based `position`.
    fn on_genome(&self, position: u64) -> u64 {
        if self.reverse {
            self.end + 1 - position
        } else {
            self.start + position - 1
        }
    }

    /// A VCF record named `id` that replaces the base at the motif's
    /// `position` with the next nucleotide of [`NUCLEOTIDES`].
    fn record(&self, position: u64, id: &str) -> String {
        let at = self.on_genome(position);
        let reference = tree_genome().base("21", at).unwrap();
        let held = NUCLEOTIDES.iter().position(|n| *n == reference).unwrap();
        let alternate = NUCLEOTIDES[(held + 1) % 4];
        format!(
            "21\t{at}\t{id}\t{}\t{}\t.\t.\t.\n",
            char::from(reference),
            char::from(alternate)
        )
    }
}

/// The counts of a made-up matrix of `length` positions, one `[A, C, G, T]`
/// row a position: an odd position holds 1000 for the nucleotide `favoured`
/// names for it and none for the others, an even position 250 for each.
fn made_up_counts(length: u64, favoured: impl Fn(u64) -> u8) -> Vec<[u64; 4]> {
    (1..=length)
        .map(|position| {
            if position % 2 == 0 {
                return [250; 4];
            }
            let favoured = favoured(position);
            NUCLEOTIDES.map(|nucleotide| if nucleotide == favoured { 1000 } else { 0 })
        })
        .collect()
}

/// The matrices the motif fixture names, each with the made-up counts of as
/// many positions as its motifs span. The matrix of [`FORWARD_MOTIF`] and of
/// [`REVERSE_MOTIF`] favours, at each odd position, the base that motif
/// holds there; any other matrix favours the nucleotides in turn.
fn made_up_matrices() -> BTreeMap<String, Vec<[u64; 4]>> {
    let mut lengths: BTreeMap<String, u64> = BTreeMap::new();
    for record in gff_records("motif_features.gff3.gz") {
        lengths.insert(
            record.attributes["binding_matrix_id"].clone(),
            record.end - record.start + 1,
        );
    }
    lengths
        .into_iter()
        .map(|(matrix, length)| {
            let scored = [&FORWARD_MOTIF, &REVERSE_MOTIF]
                .into_iter()
                .find(|motif| motif.matrix == matrix)
                .map(ScoredMotif::sequence);
            let counts = made_up_counts(length, |position| match &scored {
                Some(sequence) => sequence[position as usize - 1],
                None => NUCLEOTIDES[position as usize % 4],
            });
            (matrix, counts)
        })
        .collect()
}

/// The version the made-up `regulatory_build` table states.
const REGULATORY_BUILD: &str = "7.5";

/// The made-up tables of a funcgen database. `binding_matrix` holds a row
/// for every matrix of [`made_up_matrices`], numbered from 101 in the order
/// of their stable identifiers, and one for a matrix no motif names;
/// `binding_matrix_frequencies` holds their counts as a dump orders them,
/// every position of A, then of C, G and T. The matrix of [`FORWARD_MOTIF`]
/// is linked to complex 12 and then to complex 5. `regulatory_build` holds
/// one row, marked current.
fn funcgen_tables() -> Vec<(&'static str, String)> {
    let (mut matrices, mut frequencies, mut links) = (String::new(), String::new(), String::new());
    let mut row = 0;
    for (number, (stable_id, counts)) in made_up_matrices().iter().enumerate() {
        let id = 101 + number;
        writeln!(
            matrices,
            "{id}\tmade-up matrix {id}\t4.4\tSELEX\t{stable_id}"
        )
        .unwrap();
        for (nucleotide, letter) in NUCLEOTIDES.iter().enumerate() {
            for (position, counts) in counts.iter().enumerate() {
                row += 1;
                writeln!(
                    frequencies,
                    "{row}\t{id}\t{}\t{}\t{}",
                    position + 1,
                    char::from(*letter),
                    counts[nucleotide]
                )
                .unwrap();
            }
        }
        if stable_id == FORWARD_MOTIF.matrix {
            writeln!(links, "1\t{id}\t12\n2\t{id}\t5").unwrap();
        }
    }
    matrices.push_str("999\tmatrix of no motif\t4.4\tSELEX\tENSPFM9999\n");
    frequencies.push_str("999999\t999\t1\tA\t1\n");
    vec![
        ("binding_matrix", matrices),
        ("binding_matrix_frequencies", frequencies),
        ("binding_matrix_transcription_factor_complex", links),
        (
            "transcription_factor_complex",
            "5\tfactor_b_factor_c\tFactorB::FactorC\n12\tfactor_a\tFactorA\n".to_string(),
        ),
        (
            "regulatory_build",
            format!(
                "1\ta build\t116\ta build\t{REGULATORY_BUILD}\t2024-7\t2024-7\t19\t16\t1\t\\N\n"
            ),
        ),
    ]
}

const EVERY_PART: Parts = Parts {
    variation: Variation::Files,
    database: Database::EveryTable,
    regulation: Regulation::EveryFile,
};

/// The INFO keys the 1000 Genomes file of a tree adds to the known variant:
/// one frequency a continental population.
const FREQUENCIES: &str = ";AFR=0.25;AMR=0.5;EAS=0.125;EUR=0.75;SAS=0.0625";

/// A BGZF VCF holding the known variant of every tree, a made-up dbSNP record
/// on chromosome 21, in the layout of Ensembl's variation files; `info` ends
/// its INFO column.
fn variation_vcf(info: &str) -> Vec<u8> {
    variation_vcf_on("21", info)
}

/// The file of [`variation_vcf`] with its record on `sequence`.
fn variation_vcf_on(sequence: &str, info: &str) -> Vec<u8> {
    let mut vcf = bgzf::io::Writer::new(Vec::new());
    write!(
        vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=dbSNP_156,Number=0,Type=Flag,Description=\"dbSNP\">\n\
         ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
         ##INFO=<ID=AFR,Number=A,Type=Float,Description=\"AFR\">\n\
         ##INFO=<ID=AMR,Number=A,Type=Float,Description=\"AMR\">\n\
         ##INFO=<ID=EAS,Number=A,Type=Float,Description=\"EAS\">\n\
         ##INFO=<ID=EUR,Number=A,Type=Float,Description=\"EUR\">\n\
         ##INFO=<ID=SAS,Number=A,Type=Float,Description=\"SAS\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         {sequence}\t18000300\trs18000300\tT\tC\t.\t.\tdbSNP_156;TSA=SNV{info}\n"
    )
    .unwrap();
    vcf.finish().unwrap()
}

/// The second variant of a tree's chromosome 21 file, which the made-up
/// tables and the made-up GVF mark as failed.
const FAILED_IN_THE_VCF: &str = "rs18000350";

/// A failed variant the made-up GVF holds and no VCF does, a substitution.
const FAILED_SUBSTITUTION: &str = "rs18000400";

/// A failed variant the made-up GVF holds and no VCF does, an insertion.
const FAILED_INSERTION: &str = "rs18000450";

/// A BGZF VCF holding the known variant of every tree and
/// [`FAILED_IN_THE_VCF`], both made-up dbSNP records on chromosome 21.
fn two_record_vcf() -> Vec<u8> {
    let mut vcf = bgzf::io::Writer::new(Vec::new());
    write!(
        vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=dbSNP_156,Number=0,Type=Flag,Description=\"dbSNP\">\n\
         ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t18000300\trs18000300\tT\tC\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t18000350\t{FAILED_IN_THE_VCF}\tG\tA\t.\t.\tdbSNP_156;TSA=SNV\n"
    )
    .unwrap();
    vcf.finish().unwrap()
}

/// The made-up GVF file of the failed variants, gzip-compressed, in the
/// layout of Ensembl's: [`FAILED_IN_THE_VCF`], which the chromosome 21 file
/// of a tree with two records holds; [`FAILED_SUBSTITUTION`] and
/// [`FAILED_INSERTION`], which no VCF holds, the insertion with its start
/// written on its end; and a variant on chromosome 20, of which a tree
/// publishes no VCF.
fn failed_gvf() -> Vec<u8> {
    let feature = |sequence: &str, class: &str, at: u64, name: &str, alleles: &str| {
        format!(
            "{sequence}\tdbSNP\t{class}\t{at}\t{at}\t.\t+\t.\t\
             ID=1;Variant_seq={};Dbxref=dbSNP_156:{name};Reference_seq={};\
             ensembl_failure_reason=Variant maps to more than one genomic location\n",
            alleles.split_once('/').unwrap().1,
            alleles.split_once('/').unwrap().0
        )
    };
    let text = [
        "##gff-version 3\n##gvf-version 1.07\n".to_string(),
        feature("20", "SNV", 18_000_400, "rs20000400", "A/G"),
        feature("21", "SNV", 18_000_350, FAILED_IN_THE_VCF, "G/A"),
        feature("21", "SNV", 18_000_400, FAILED_SUBSTITUTION, "A/G"),
        feature("21", "insertion", 18_000_450, FAILED_INSERTION, "-/TT"),
    ]
    .concat();
    gzip(text.as_bytes())
}

/// The PubMed identifiers of the known variant of a tree whose database
/// publishes the made-up tables, in the order of its publications'
/// identifiers.
const KNOWN_PUBMED: &str = "30000021,30000023";

/// The synonyms of the known variant of such a tree: the source that byte
/// order puts second first, and the names of the other from the table's
/// last row back.
const KNOWN_SYNONYMS: &str = "another_db::ad9--Made_up_db::mud2,mud1";

/// The made-up tables the three cells are read from. `variation` names the
/// known variant of every tree, [`FAILED_IN_THE_VCF`], [`FAILED_SUBSTITUTION`]
/// and a variant no file places, and is left out without `with_variation`.
/// Three publications, one without a PubMed identifier, cite the known
/// variant, in no order of their identifiers, and one of them cites
/// [`FAILED_SUBSTITUTION`]. The known variant has one synonym under a source
/// named for dbSNP and three under two other sources. `failed_variation`
/// names the two failed variants. The rows of `variation` are in no order of
/// their identifiers either.
fn variant_cell_tables(with_variation: bool) -> Vec<(&'static str, Vec<u8>)> {
    let variation: String = [
        (9, FAILED_SUBSTITUTION),
        (7, "rs18000300"),
        (10, "rs10000010"),
        (8, FAILED_IN_THE_VCF),
    ]
    .iter()
    .map(|(id, name)| format!("{id}\t1\t{name}\t\\N\t2\t0\t\\N\t\\N\t\\N\t\\N\t\\N\t1\n"))
    .collect();
    let mut publication = b"23\t".to_vec();
    publication.extend(escaped(b"a title\twith a tab\nand a line end"));
    publication.extend(b"\tAn Author\t30000023\t\\N\t2020\t\\N\t\\N\n");
    publication.extend(b"21\ta title\tAn Author\t30000021\t\\N\t2019\t\\N\t\\N\n");
    publication.extend(b"22\ta title\tAn Author\t\\N\t\\N\t2018\t\\N\t\\N\n");
    let source = "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n\
                  2\tArchive dbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation_synonym\n\
                  3\tMade_up_db\t1\t\\N\t\\N\t\\N\tgermline\tvariation_synonym\n\
                  4\tanother_db\t1\t\\N\t\\N\t\\N\tgermline\tvariation_synonym\n";
    let mut tables = vec![
        ("publication", publication),
        (
            "variation_citation",
            b"7\t23\t\\N\n9\t21\t\\N\n7\t22\t\\N\n7\t21\t\\N\n".to_vec(),
        ),
        ("source", source.as_bytes().to_vec()),
        (
            "variation_synonym",
            b"1\t7\t\\N\t2\trs17000300\n2\t7\t\\N\t3\tmud1\n\
              3\t7\t\\N\t4\tad9\n4\t7\t\\N\t3\tmud2\n"
                .to_vec(),
        ),
        ("failed_variation", b"1\t9\t5\n2\t8\t2\n".to_vec()),
    ];
    if with_variation {
        tables.push(("variation", variation.into_bytes()));
    }
    tables
}

/// The transcript whose peptide the made-up matrices are stored for: the one
/// transcript of a single-exon gene of the fixture.
const MATRIX_TRANSCRIPT: &str = "ENST00000334046";

/// The peptide of [`MATRIX_TRANSCRIPT`] in the fixture's protein FASTA.
fn matrix_peptide() -> String {
    let fasta = std::fs::read_to_string(fixture_dir().join("peptides.fa")).unwrap();
    let named = format!("transcript:{MATRIX_TRANSCRIPT}.");
    let record = fasta
        .split('>')
        .find(|record| record.lines().next().is_some_and(|h| h.contains(&named)))
        .expect("the fixture holds the transcript's peptide");
    record.lines().skip(1).map(str::trim).collect()
}

/// The analyses of the made-up `attrib` table, each with its identifier and
/// the category every cell of its made-up matrix holds: `deleterious` for
/// SIFT, `benign` for PolyPhen HumVar. The last is no analysis a cache
/// carries.
static ANALYSES: [(u32, &str, u16); 4] = [
    (267, "sift", 1),
    (268, "polyphen_humvar", 2),
    (269, "polyphen_humdiv", 0),
    (270, "cadd", 3),
];

/// The amino acids of a matrix position, in the order of its cells.
const AMINO_ACIDS: &str = "ACDEFGHIKLMNPQRSTVWY";

/// A made-up prediction matrix.
struct Matrix {
    /// `VEP` and 40 bytes a residue: 20 cells, each a little-endian `u16`
    /// with the category in its two high bits and the score, in thousandths,
    /// in its ten low bits.
    cells: Vec<u8>,
    /// The cells as the variation database stores them, gzip-compressed.
    blob: Vec<u8>,
}

impl Matrix {
    /// The matrix of `analysis` for the peptide of [`MATRIX_TRANSCRIPT`]:
    /// every cell holds the category [`ANALYSES`] names and a score that
    /// differs from cell to cell. The compressed bytes hold a NUL, a tab, a
    /// line end and a backslash, the four bytes a table dump writes behind a
    /// backslash: the matrix is the first of a run of such matrices whose
    /// compressed bytes do.
    fn of(analysis: &str) -> Self {
        let residues = matrix_peptide().len();
        let category = ANALYSES
            .iter()
            .find_map(|(_, name, category)| (*name == analysis).then_some(*category))
            .expect("a made-up analysis");
        (0..64usize)
            .find_map(|salt| {
                let mut cells = b"VEP".to_vec();
                for cell in 0..residues * AMINO_ACIDS.len() {
                    let score = (cell * (7919 + 2 * salt) + usize::from(category) * 331) % 1001;
                    cells.extend_from_slice(&((category << 14) | score as u16).to_le_bytes());
                }
                let blob = gzip(&cells);
                [0, b'\t', b'\n', b'\\']
                    .iter()
                    .all(|byte| blob.contains(byte))
                    .then_some(Self { cells, blob })
            })
            .expect("a matrix whose compressed bytes hold the four escaped bytes")
    }

    /// The score of `amino_acid` at the 1-based `position`.
    fn score(&self, position: usize, amino_acid: char) -> f64 {
        let at = 3 + (position - 1) * 40 + AMINO_ACIDS.find(amino_acid).unwrap() * 2;
        let cell = u16::from_le_bytes([self.cells[at], self.cells[at + 1]]);
        f64::from(cell & 0x3ff) / 1000.0
    }
}

/// `value` as a table dump writes it: a backslash before a backslash, a tab
/// and a line end, and `\0` for a NUL.
fn escaped(value: &[u8]) -> Vec<u8> {
    let mut written = Vec::new();
    for &byte in value {
        match byte {
            0 => written.extend_from_slice(b"\\0"),
            b'\\' | b'\t' | b'\n' => written.extend_from_slice(&[b'\\', byte]),
            _ => written.push(byte),
        }
    }
    written
}

/// The made-up tables of a variation database beside the fixture's two
/// phenotype tables: those the ClinVar columns are read from, and those the
/// matrices are read from as `database` publishes them. `attrib_type` holds
/// the type of the analyses beside the ClinVar significance; `attrib` names
/// the analyses under that type; `translation_md5` holds the MD5 of the
/// peptide of [`MATRIX_TRANSCRIPT`] and one no peptide of the fixture has;
/// each of the two MD5s has a SIFT matrix, and the first one a matrix of
/// every other analysis of [`ANALYSES`].
fn variation_tables(database: Database) -> Vec<(&'static str, Vec<u8>)> {
    let held = format!("{:x}", md5::compute(matrix_peptide()));
    let other = format!("{:x}", md5::compute("no peptide of the fixture"));
    let attrib: String = ANALYSES
        .iter()
        .filter(|(_, name, _)| database != Database::SiftAlone || !name.starts_with("polyphen"))
        .map(|(id, name, _)| format!("{id}\t497\t{name}\n"))
        .collect();
    let mut matrices = Vec::new();
    let rows = ANALYSES.iter().map(|(id, name, _)| (7, id, name));
    for (md5, analysis, name) in rows.chain([(8, &ANALYSES[0].0, &ANALYSES[0].1)]) {
        matrices.extend_from_slice(format!("{md5}\t{analysis}\t").as_bytes());
        matrices.extend_from_slice(&escaped(&Matrix::of(name).blob));
        matrices.push(b'\n');
    }
    let mut tables = vec![
        (
            "source",
            b"1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n".to_vec(),
        ),
        ("seq_region", b"131537\t21\t4\n".to_vec()),
        (
            "attrib_type",
            b"480\tclinvar_clin_sig\tClinVar clinical significance\t\\N\n\
              497\tprot_func_analysis\tProtein function analysis\t\\N\n"
                .to_vec(),
        ),
        ("phenotype_feature_attrib", Vec::new()),
        ("attrib", attrib.into_bytes()),
        ("protein_function_predictions", matrices),
    ];
    if database != Database::NoTranslationMd5 {
        tables.push((
            "translation_md5",
            format!("7\t{held}\n8\t{other}\n").into_bytes(),
        ));
    }
    if matches!(
        database,
        Database::VariantTables | Database::NoVariationTable
    ) {
        tables.retain(|(table, _)| *table != "source");
        tables.extend(variant_cell_tables(database == Database::VariantTables));
    }
    if database == Database::FeatureTable {
        for (table, rows) in &mut tables {
            match *table {
                "seq_region" => rows.extend_from_slice(b"131550\tMT\t4\n"),
                "attrib" => rows.extend_from_slice(
                    b"367\t498\tMultiple_observations\n418\t498\tPhenotype_or_Disease\n",
                ),
                _ => {}
            }
        }
        tables.push(("variation_feature", feature_table().into_bytes()));
    }
    tables
}

/// The known variant the made-up `variation_feature` table holds where the
/// records of an annotated run name one, on chromosome 21.
const TABLE_VARIANT: &str = "rs21000300";

/// The insertion the made-up `variation_feature` table holds on chromosome
/// 21.
const TABLE_INSERTION: &str = "rs21000451";

/// The made-up `variation_feature` table, its rows in no order of their
/// starts. On chromosome 21: [`TABLE_INSERTION`], its start past its end;
/// [`TABLE_VARIANT`], with two clinical significance terms and, among its
/// evidence, the one `attrib` names for the phenotype-or-disease flag; and a
/// variant that is not on display. On `MT`: one variant on display.
fn feature_table() -> String {
    let row = |id: u32, region: u32, span: (u64, u64), name: &str, alleles: &str, rest: &str| {
        format!(
            "{id}\t{region}\t{}\t{}\t1\t{id}\t{alleles}\t\\N\t{name}\t1\t\\N\t1\t\
             intergenic_variant\t\t2\t0\t\\N\t\\N\t\\N\t\\N\t{rest}\n",
            span.0, span.1
        )
    };
    [
        row(
            1,
            131_537,
            (18_000_451, 18_000_450),
            TABLE_INSERTION,
            "-/TT",
            "\\N\t\\N\t1",
        ),
        row(
            2,
            131_537,
            (18_000_300, 18_000_300),
            TABLE_VARIANT,
            "T/C",
            "367,418\tlikely pathogenic,pathogenic\t1",
        ),
        row(
            3,
            131_537,
            (18_000_500, 18_000_500),
            "rs21000500",
            "G/A",
            "\\N\t\\N\t0",
        ),
        row(4, 131_550, (5000, 5000), "rs25000", "A/G", "\\N\t\\N\t1"),
    ]
    .concat()
}

/// Lays out release 116 of `genome` as the Ensembl site publishes it: the
/// fixture's GFF3, GTF, peptides and core tables, the indexed genome, and of
/// each part what `parts` names.
fn tree(genome: Genome, parts: Parts) -> Site {
    let Genome {
        species, assembly, ..
    } = genome;
    let site = Site::new();
    let fixture = fixture_dir();
    let release = "pub/release-116";
    let named = format!("{}.{assembly}", capitalized(species));
    site.add_gzipped(
        &format!("{release}/gff3/{species}"),
        &format!("{named}.116.gff3.gz"),
        &fixture.join("slice.gff3"),
    );
    site.add_gzipped(
        &format!("{release}/gtf/{species}"),
        &format!("{named}.116.gtf.gz"),
        &fixture.join("slice.gtf"),
    );
    site.add_gzipped(
        &format!("{release}/fasta/{species}/pep"),
        &format!("{named}.pep.all.fa.gz"),
        &fixture.join("peptides.fa"),
    );
    let core = format!("{release}/mysql/{species}_core_116_{}", genome.database);
    for table in CORE_TABLES {
        site.add_gzipped(
            &core,
            &format!("{table}.txt.gz"),
            &fixture.join(format!("{table}.txt")),
        );
    }
    for suffix in ["", ".fai", ".gzi"] {
        site.add_copy(
            &format!("{release}/fasta/{species}/dna_index"),
            &format!("{named}.dna.toplevel.fa.gz{suffix}"),
            &golden(&format!("GRCh37-hgvs/reference.fa.gz{suffix}")),
        );
    }

    let vcfs = format!("{release}/variation/vcf/{species}");
    match parts.variation {
        Variation::Files if species == "homo_sapiens" => {
            site.add(&vcfs, "homo_sapiens-chr21.vcf.gz", &variation_vcf(""));
            site.add(
                &vcfs,
                "1000GENOMES-phase_3.vcf.gz",
                &variation_vcf(FREQUENCIES),
            );
        }
        Variation::Files => site.add(&vcfs, &format!("{species}.vcf.gz"), &variation_vcf("")),
        Variation::OffTheModel if species == "homo_sapiens" => {
            site.add(
                &vcfs,
                "homo_sapiens-chrMT.vcf.gz",
                &variation_vcf_on("MT", ""),
            );
            site.add(
                &vcfs,
                "1000GENOMES-phase_3.vcf.gz",
                &variation_vcf(FREQUENCIES),
            );
        }
        Variation::OffTheModel | Variation::OffTheModelAndFailedGvf => {
            site.add(
                &vcfs,
                &format!("{species}.vcf.gz"),
                &variation_vcf_on("MT", ""),
            );
            if parts.variation == Variation::OffTheModelAndFailedGvf {
                site.add(
                    &format!("{release}/variation/gvf/{species}"),
                    &format!("{species}_failed.gvf.gz"),
                    &failed_gvf(),
                );
            }
        }
        Variation::ReadmeAlone => site.add(&vcfs, "README", b"no variation file\n"),
        Variation::TwoRecords | Variation::TwoRecordsAndFailedGvf => {
            site.add(&vcfs, "homo_sapiens-chr21.vcf.gz", &two_record_vcf());
            site.add(
                &vcfs,
                "1000GENOMES-phase_3.vcf.gz",
                &variation_vcf(FREQUENCIES),
            );
            if parts.variation == Variation::TwoRecordsAndFailedGvf {
                site.add(
                    &format!("{release}/variation/gvf/{species}"),
                    &format!("{species}_failed.gvf.gz"),
                    &failed_gvf(),
                );
            }
        }
        Variation::Absent => {}
    }

    if parts.database != Database::Absent {
        let variation = format!(
            "{release}/mysql/{species}_variation_116_{}",
            genome.database
        );
        for table in ["phenotype", "phenotype_feature"] {
            site.add_gzipped(
                &variation,
                &format!("{table}.txt.gz"),
                &fixture.join(format!("{table}.txt")),
            );
        }
        for (table, rows) in variation_tables(parts.database) {
            site.add(&variation, &format!("{table}.txt.gz"), &gzip(&rows));
        }
    }

    if parts.regulation != Regulation::Absent {
        let annotation = format!(
            "{release}/regulation/{species}/{}/annotation",
            genome.regulation
        );
        let named = format!("{}.{}", capitalized(species), genome.regulation);
        let mut files = vec![
            (
                "regulatory_features.v116.gff3.gz",
                "regulatory_features.gff3.gz",
            ),
            ("EMARs.v116.gff.gz", "EMARs.gff.gz"),
        ];
        if parts.regulation != Regulation::NoMotifFile {
            files.push(("motif_features.v116.gff3.gz", "motif_features.gff3.gz"));
        }
        for (name, source) in files {
            site.add_copy(
                &annotation,
                &format!("{named}.{name}"),
                &regulatory_gff(source),
            );
        }
        if matches!(
            parts.regulation,
            Regulation::ActivityAndTables | Regulation::ActivityAlone
        ) {
            site.add(
                &annotation,
                &format!("{named}.regulatory_activity.v116.tsv.gz"),
                &activity_file(),
            );
        }
        if matches!(
            parts.regulation,
            Regulation::ActivityAndTables | Regulation::TablesAlone
        ) {
            let funcgen = format!("{release}/mysql/{species}_funcgen_116_{}", genome.database);
            for (table, rows) in funcgen_tables() {
                site.add(&funcgen, &format!("{table}.txt.gz"), &gzip(rows.as_bytes()));
            }
        }
    }
    site
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// A release installed from a tree.
struct Installation {
    genome: Genome,
    output: Output,
    /// The cache root the run was given as `--dir_cache`.
    dir_cache: PathBuf,
    root: tempfile::TempDir,
}

/// Installs release 116 of `genome` from the tree of `parts` with
/// `--auto <species> <assembly>` and `flags`; the run exits 0.
fn install(genome: Genome, parts: Parts, flags: &[&str]) -> Installation {
    let site = tree(genome, parts);
    let root = tempfile::tempdir().unwrap();
    let dir_cache = root.path().join("vep");
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--auto", genome.species, genome.assembly])
        .args(["--release", "116"])
        .arg("--dir_cache")
        .arg(&dir_cache)
        .arg("--site")
        .arg(site.path())
        .args(flags)
        .output()
        .expect("vep-cache-builder runs");
    assert!(output.status.success(), "{}", stderr(&output));
    Installation {
        genome,
        output,
        dir_cache,
        root,
    }
}

/// The flags that name chromosome 21, the one chromosome a human tree
/// publishes a variation VCF of.
const CHROMOSOME_21: [&str; 2] = ["--chromosomes", "21"];

/// The closing warning of the gnomAD frequencies, which no file of the
/// release carries: every `homo_sapiens` installation with known variants
/// ends on it.
const NO_GNOMAD: &str = "WARN vep_cache_builder::install: No gnomAD frequencies in the cache of homo_sapiens GRCh38: no file of the release carries them; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The closing warning of the PubMed identifiers and the variant synonyms
/// of a `homo_sapiens` tree whose variation database lists no `variation`
/// table.
const NO_TABLE_CELLS: &str = "WARN vep_cache_builder::install: No PubMed identifiers or variant synonyms in the cache of homo_sapiens GRCh38: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no variation table; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The closing warning of the failed variants of a `homo_sapiens` tree with
/// no GVF directory and no `failed_variation` table to flag them from.
const NO_FAILED_VARIANTS: &str = "WARN vep_cache_builder::install: No failed variants in the cache of homo_sapiens GRCh38: pub/release-116/variation/gvf/homo_sapiens/CHECKSUMS is not on the site; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The closing warning of the activity of the regulatory features of a
/// `homo_sapiens` tree whose annotation directory lists no activity file.
const NO_ACTIVITY: &str = "WARN vep_cache_builder::install: No cell-type activity of regulatory features in the cache of homo_sapiens GRCh38: pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS lists no regulatory_activity file; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The closing warning of the matrix counts of the motif features of a
/// `homo_sapiens` tree with no funcgen database.
const NO_MOTIF_SCORES: &str = "WARN vep_cache_builder::install: No HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features in the cache of homo_sapiens GRCh38: pub/release-116/mysql/homo_sapiens_funcgen_116_38/CHECKSUMS is not on the site; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The line of an installation that left no part out, for `homo_sapiens`.
const EVERY_PART_LINE: &str = "INFO vep_cache_builder::install: The cache of homo_sapiens GRCh38 holds each part of Ensembl VEP's own cache that a file of the release carries: known variants with their 1000 Genomes frequencies, SIFT and PolyPhen predictions, regulatory features and motif features";

/// The closing warning of an installation that read the 1000 Genomes file:
/// the cache's `AF` is recovered from the continental frequencies.
const RECOVERED_AF: &str = "WARN vep_cache_builder::install: The 1000 Genomes AF in the cache of homo_sapiens GRCh38 is each allele's count over the panel, recovered from the continental frequencies, and can differ from the AF of Ensembl VEP's own cache, which a conversion carries (docs/cache-variation-data.md#the-1000-genomes-frequencies)";

/// The variable that carries the flags of a run to
/// [`annotation_run_of_a_child_process`], one a line.
const RUN_FLAGS: &str = "VEP_TEST_ANNOTATION_FLAGS";

/// The records of an annotated run: one inside two motif features and one
/// inside an enhancer of the regulatory GFF files, the known variant of the
/// trees with its alleles, and a substitution in the codon of residue 35 of
/// [`MATRIX_TRANSCRIPT`].
const RECORDS: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    21\t5010255\tin_motif\tA\tG\t.\t.\t.\n\
    21\t5011400\tin_enhancer\tA\tG\t.\t.\t.\n\
    21\t18000300\tknown\tT\tC\t.\t.\t.\n\
    21\t30541731\tsubstitution\tT\tC\t.\t.\t.\n";

impl Installation {
    /// The release directory, `<dir_cache>/<species>/116_<assembly>`.
    fn release(&self) -> PathBuf {
        self.dir_cache
            .join(self.genome.species)
            .join(format!("116_{}", self.genome.assembly))
    }

    fn info(&self) -> Value {
        let info = std::fs::read_to_string(self.release().join("info.json")).unwrap();
        serde_json::from_str(&info).unwrap()
    }

    /// The line of the run that names the installed directory, from its
    /// level on.
    fn installed_line(&self) -> String {
        format!(
            "INFO vep_cache_builder::install: Installed {}: run vep with --dir_cache {} --species {} --cache_version 116 --assembly {}",
            self.release().display(),
            self.dir_cache.display(),
            self.genome.species,
            self.genome.assembly
        )
    }

    /// The lines the run logged from the one naming the installed directory
    /// on, each from its level on.
    fn closing_lines(&self) -> Vec<String> {
        let log = stderr(&self.output);
        let logged: Vec<&str> = log
            .lines()
            .filter_map(|line| line.split_once(' '))
            .map(|(_, logged)| logged.trim_start())
            .collect();
        let installed = logged
            .iter()
            .position(|line| line.starts_with("INFO vep_cache_builder::install: Installed "))
            .unwrap_or_else(|| panic!("no line names the installed directory: {log}"));
        logged[installed..]
            .iter()
            .map(|line| line.to_string())
            .collect()
    }

    /// The warnings the run logged, each from its module's name on.
    fn warnings(&self) -> Vec<String> {
        stderr(&self.output)
            .lines()
            .filter_map(|line| line.split_once(" WARN "))
            .map(|(_, text)| text.to_string())
            .collect()
    }

    /// Every file of the release directory, relative to it, with its bytes.
    fn files(&self) -> BTreeMap<PathBuf, Vec<u8>> {
        fn walk(root: &Path, dir: &Path, files: &mut BTreeMap<PathBuf, Vec<u8>>) {
            for entry in std::fs::read_dir(dir).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    walk(root, &path, files);
                } else {
                    let relative = path.strip_prefix(root).unwrap().to_path_buf();
                    files.insert(relative, std::fs::read(&path).unwrap());
                }
            }
        }
        let mut files = BTreeMap::new();
        walk(&self.release(), &self.release(), &mut files);
        files
    }

    /// The records of every JSON shard under `tree` (`transcripts`,
    /// `regulatory`) of the release directory.
    fn records(&self, tree: &str) -> Vec<Value> {
        self.files()
            .iter()
            .filter(|(file, _)| file.starts_with(tree) && file.components().count() == 3)
            .flat_map(|(_, shard)| serde_json::from_slice::<Vec<Value>>(shard).unwrap())
            .collect()
    }

    /// Annotates [`RECORDS`] over the installed release with `--everything`,
    /// as the `vep` binary would, and returns the rows of the output, each
    /// split into its columns, and what the run wrote to standard error.
    fn annotate(&self) -> (Vec<Vec<String>>, String) {
        let run = self.run_vep(RECORDS, &["--everything"]);
        (run.rows, run.notices)
    }

    /// Annotates `records` over the installed release with `flags`, as the
    /// `vep` binary would. The run is this test binary started again for
    /// [`annotation_run_of_a_child_process`] alone, so the standard error
    /// read is the run's own: the notices of a run are written there.
    fn run_vep(&self, records: &str, flags: &[&str]) -> Run {
        let input = self.root.path().join("records.vcf");
        let annotated = self.root.path().join("annotated.txt");
        std::fs::write(&input, records).unwrap();
        let dir_cache = self.dir_cache.display().to_string();
        let (input, annotated_s) = (input.display().to_string(), annotated.display().to_string());
        let mut flags = flags.to_vec();
        flags.splice(
            0..0,
            [
                "vep",
                "--offline",
                "--cache",
                "--no_stats",
                "--force_overwrite",
            ],
        );
        flags.extend([
            "-i",
            &input,
            "-o",
            &annotated_s,
            "--dir_cache",
            &dir_cache,
            "--species",
            self.genome.species,
            "--cache_version",
            "116",
            "--assembly",
            self.genome.assembly,
        ]);
        let run = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "annotation_run_of_a_child_process"])
            .arg("--nocapture")
            .env(RUN_FLAGS, flags.join("\n"))
            .output()
            .expect("the test binary runs");
        let notices = stderr(&run);
        assert!(run.status.success(), "{notices}");
        let output = std::fs::read_to_string(&annotated).unwrap();
        let (header, rows): (Vec<&str>, Vec<&str>) =
            output.lines().partition(|line| line.starts_with('#'));
        Run {
            header: header.into_iter().map(str::to_string).collect(),
            rows: rows
                .into_iter()
                .map(|line| line.split('\t').map(str::to_string).collect())
                .collect(),
            notices,
        }
    }
}

/// An annotated run over an installed release.
struct Run {
    /// The lines of the output that open with `#`.
    header: Vec<String>,
    /// The rows of the output, each split into its columns.
    rows: Vec<Vec<String>>,
    /// What the run wrote to standard error.
    notices: String,
}

/// The columns of an output row a test reads.
const UPLOADED_VARIATION: usize = 0;
const FEATURE: usize = 4;
const FEATURE_TYPE: usize = 5;
const CONSEQUENCE: usize = 6;
const PROTEIN_POSITION: usize = 9;
const AMINO_ACIDS_CELL: usize = 10;
const EXISTING_VARIATION: usize = 12;
const EXTRA: usize = 13;

/// Runs `vep` in-process with the flags the variable [`RUN_FLAGS`] carries,
/// and does nothing without the variable: the run of the child process
/// `Installation::annotate` starts.
#[test]
fn annotation_run_of_a_child_process() {
    let Ok(flags) = std::env::var(RUN_FLAGS) else {
        return;
    };
    let config = vep_cli::Config::from_argv(flags.lines()).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
}

/// A GRCh38 tree that publishes every part installs every part with no flag
/// but `--chromosomes`. The cache carries the variation file of chromosome 21
/// with the 1000 Genomes columns, the three matrices of the one peptide the
/// variation database holds matrices for, byte for byte as its dump writes
/// them behind the escapes, and the regulatory tree of the two GFF3 files
/// with no record of the EMAR file beside them. The tree publishes no file
/// of a cell of the known variants or of the regulatory tree: the run ends
/// on the line that the cache holds each part, the warnings of those cells,
/// each with the listing that lacks its file or table, and the warning that
/// its `AF` is recovered from the continental frequencies, and `info.json`
/// names each of them under `left_out`.
/// `vep --everything` over the installed release
/// then switches off neither SIFT nor PolyPhen, finds the regulatory and the
/// variation data, prints one line each for the recovered `AF`, for the
/// columns of the cells the installation named and for `MAX_AF`, and fills
/// the known variant's identifier and frequency,
/// the SIFT and PolyPhen cells of the two matrices and a row for each
/// regulatory and motif feature.
#[test]
fn every_part_the_tree_publishes_is_installed() {
    let installed = install(HUMAN, EVERY_PART, &CHROMOSOME_21);
    let release = installed.release();
    let info = installed.info();

    assert!(release.join("variations/21/all_vars.gz").is_file());
    assert!(release.join("variations/21/all_vars.gz.csi").is_file());
    let columns = info["variation_cols"].as_array();
    assert!(
        columns.is_some_and(|columns| columns.iter().any(|column| column == "AF")),
        "{info}"
    );

    assert_eq!(
        info["left_out"],
        serde_json::json!([
            "panel_frequency",
            "gnomad_frequencies",
            "pubmed",
            "var_synonyms",
            "failed_variants",
            "cell_type_activity",
            "motif_scores"
        ])
    );

    assert_eq!(info["sift"], "b", "{info}");
    assert_eq!(info["polyphen"], "b", "{info}");
    let with_matrices: Vec<Value> = installed
        .records("transcripts")
        .into_iter()
        .filter(|record| {
            record["variation_effect_feature_cache"]["protein_function_predictions"]
                .as_object()
                .is_some_and(|matrices| !matrices.is_empty())
        })
        .collect();
    assert_eq!(with_matrices.len(), 1);
    assert_eq!(with_matrices[0]["stable_id"], MATRIX_TRANSCRIPT);
    let matrices =
        &with_matrices[0]["variation_effect_feature_cache"]["protein_function_predictions"];
    assert_eq!(matrices.as_object().unwrap().len(), 3, "{matrices}");
    for analysis in ["sift", "polyphen_humvar", "polyphen_humdiv"] {
        assert_eq!(
            matrices[analysis]["predictions_data"],
            base64::engine::general_purpose::STANDARD.encode(Matrix::of(analysis).blob),
            "{analysis}"
        );
    }

    assert_eq!(info["regulatory"], true, "{info}");
    assert!(info["regulatory_schema"].is_u64(), "{info}");
    let features = installed.records("regulatory");
    let of_type = |feature_type: &str| {
        features
            .iter()
            .filter(|feature| feature["feature_type"] == feature_type)
            .count()
    };
    assert_eq!(
        (of_type("RegulatoryFeature"), of_type("MotifFeature")),
        (516, 2200)
    );
    assert!(features.iter().all(|feature| feature["biotype"] != "EMAR"));
    let binding_matrices: Value = serde_json::from_slice(
        &std::fs::read(release.join("regulatory/binding_matrices.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(binding_matrices.as_object().unwrap().len(), 50);

    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );

    let (rows, notices) = installed.annotate();
    for notice in [
        "disabling SIFT",
        "disabling PolyPhen",
        "this cache holds no regulatory data",
        "no variation data",
    ] {
        assert!(!notices.contains(notice), "{notice}: {notices}");
    }
    // What the run asks for and the cache does not carry, once each: the
    // recovered `AF`, the cells the installation named, and `MAX_AF`.
    for line in [
        "INFO: --everything reads the AF of this cache, which is each allele's count over the 1000 Genomes panel, recovered from the continental frequencies, and can differ from the AF of Ensembl VEP's own cache (docs/cache-variation-data.md#the-1000-genomes-frequencies)",
        "INFO: --everything asks for the gnomAD exome columns, the gnomAD genome columns, PUBMED, HIGH_INF_POS and MOTIF_SCORE_CHANGE, which this cache does not carry: the columns stay empty; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
        "INFO: MAX_AF and MAX_AF_POPS (--everything) range over the 1000 Genomes populations alone on this cache, whose variation data holds no gnomAD exome or gnomAD genome frequency; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
    ] {
        assert_eq!(
            notices.lines().filter(|written| *written == line).count(),
            1,
            "{line}: {notices}"
        );
    }
    let of_record = |record: &str| -> Vec<&Vec<String>> {
        rows.iter()
            .filter(|row| row[UPLOADED_VARIATION] == record)
            .collect()
    };
    let known = of_record("known");
    assert!(!known.is_empty());
    for row in known {
        assert_eq!(row[EXISTING_VARIATION], "rs18000300", "{row:?}");
        assert!(
            row[EXTRA].contains(";AF=") && row[EXTRA].contains(";AFR_AF=0.25;"),
            "{row:?}"
        );
    }
    let substitution = of_record("substitution");
    assert_eq!(substitution.len(), 1, "{substitution:?}");
    let row = substitution[0];
    assert_eq!(
        [
            row[FEATURE].as_str(),
            row[CONSEQUENCE].as_str(),
            row[PROTEIN_POSITION].as_str(),
            row[AMINO_ACIDS_CELL].as_str()
        ],
        [MATRIX_TRANSCRIPT, "missense_variant", "35", "G/V"]
    );
    let sift = format!("SIFT=deleterious({})", Matrix::of("sift").score(35, 'V'));
    assert!(
        row[EXTRA].split(';').any(|cell| cell == sift),
        "{sift}: {row:?}"
    );
    let polyphen = format!(
        "PolyPhen=benign({})",
        Matrix::of("polyphen_humvar").score(35, 'V')
    );
    assert!(
        row[EXTRA].split(';').any(|cell| cell == polyphen),
        "{polyphen}: {row:?}"
    );
    let feature_types = |record: &str| -> Vec<&str> {
        of_record(record)
            .iter()
            .map(|row| row[FEATURE_TYPE].as_str())
            .collect()
    };
    assert_eq!(feature_types("in_enhancer"), ["RegulatoryFeature", "-"]);
    assert_eq!(
        feature_types("in_motif"),
        ["MotifFeature", "MotifFeature", "-"]
    );
}

/// `--include-variations` beside `--auto` changes nothing: the tree that
/// publishes every part installs the same files byte for byte with the flag
/// and without it, the variation file among them, and the two runs end on
/// the same lines.
#[test]
fn include_variations_beside_auto_changes_nothing() {
    let without = install(HUMAN, EVERY_PART, &CHROMOSOME_21);
    let with = install(
        HUMAN,
        EVERY_PART,
        &["--include-variations", "--chromosomes", "21"],
    );
    let files = without.files();
    assert!(files.contains_key(Path::new("variations/21/all_vars.gz")));
    assert!(files == with.files(), "the two caches differ");
    for run in [&without, &with] {
        assert_eq!(run.closing_lines()[0], run.installed_line());
    }
    assert_eq!(without.closing_lines().len(), 8);
    assert_eq!(without.closing_lines()[1..], with.closing_lines()[1..]);
}

/// `--chromosomes` beside `--auto` names the per-chromosome files an
/// installation reads with no `--include-variations`: the run logs no warning
/// that the flag is read only with that one, and the cache carries the
/// variation file of the chromosome named and of no other.
#[test]
fn chromosomes_beside_auto_are_read_without_include_variations() {
    let installed = install(HUMAN, EVERY_PART, &CHROMOSOME_21);
    let warnings = installed.warnings();
    assert!(
        !warnings
            .iter()
            .any(|warning| warning.contains("read only with --include-variations")),
        "{warnings:#?}"
    );
    let written: Vec<String> = std::fs::read_dir(installed.release().join("variations"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    assert_eq!(written, ["21"]);
}

/// A tree with no variation directory installs without the known variants:
/// the run ends on the warning that names the listing the site lacks, then
/// on the two warnings of the cells of the regulatory tree, and names no
/// cell of the known variants; `vep --everything` over the release says the
/// cache holds no variation data.
#[test]
fn a_tree_without_a_variation_directory_installs_without_known_variants() {
    let parts = Parts {
        variation: Variation::Absent,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No known variants or their frequencies in the cache of homo_sapiens GRCh38: pub/release-116/variation/vcf/homo_sapiens/CHECKSUMS is not on the site",
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
        ]
    );
    assert!(!installed.release().join("variations").exists());
    let (_, notices) = installed.annotate();
    assert!(notices.contains("no variation data"), "{notices}");
}

/// A variation directory whose listing names no germline file of the genome
/// installs without the known variants, with no `--chromosomes`: the run
/// ends on the warning that names the files the listing lacks.
#[test]
fn a_variation_listing_without_a_germline_file_installs_without_known_variants() {
    let parts = Parts {
        variation: Variation::ReadmeAlone,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &[]);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No known variants or their frequencies in the cache of homo_sapiens GRCh38: pub/release-116/variation/vcf/homo_sapiens/CHECKSUMS lists none of the per-chromosome files homo_sapiens-chr<N>.vcf.gz",
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
        ]
    );
    assert!(!installed.release().join("variations").exists());
}

/// A tree with no variation database installs without the matrices and with
/// the known variants, whose ClinVar columns stay empty: the run ends on one
/// warning for SIFT and PolyPhen together, which names
/// `--include-predictions`, and `vep --everything` over the release switches
/// both tools off.
#[test]
fn a_tree_without_a_variation_database_installs_without_matrices() {
    let parts = Parts {
        database: Database::Absent,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No SIFT or PolyPhen predictions in the cache of homo_sapiens GRCh38: no variation database is published for this release and assembly; --include-predictions fetches them from Ensembl's MySQL server (docs/cache-predictions.md)",
            NO_GNOMAD,
            "WARN vep_cache_builder::install: No PubMed identifiers or variant synonyms in the cache of homo_sapiens GRCh38: no variation database is published for this release and assembly; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    assert!(installed
        .release()
        .join("variations/21/all_vars.gz")
        .is_file());
    assert!(installed.warnings().contains(
        &"vep_cache_builder::install: No variation database is published: the ClinVar allele and somatic classification columns are empty"
            .to_string()
    ));
    let (_, notices) = installed.annotate();
    for notice in ["INFO: disabling SIFT", "INFO: disabling PolyPhen"] {
        assert!(notices.contains(notice), "{notice}: {notices}");
    }
}

/// A variation database that lists the ClinVar tables and no
/// `translation_md5` installs without the matrices: the run ends on one
/// warning for SIFT and PolyPhen together, which names the table the listing
/// lacks and `--include-predictions`.
#[test]
fn a_variation_database_without_the_md5_table_installs_without_matrices() {
    let parts = Parts {
        database: Database::NoTranslationMd5,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No SIFT or PolyPhen predictions in the cache of homo_sapiens GRCh38: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no translation_md5 table; --include-predictions fetches them from Ensembl's MySQL server (docs/cache-predictions.md)",
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    let info = installed.info();
    assert!(
        info["sift"].is_null() && info["polyphen"].is_null(),
        "{info}"
    );
}

/// Tables whose `attrib` names SIFT and neither PolyPhen analysis install the
/// SIFT matrices alone: the run ends on the warning for PolyPhen, which names
/// the database that holds no such matrix.
#[test]
fn a_variation_database_naming_sift_alone_installs_without_polyphen() {
    let parts = Parts {
        database: Database::SiftAlone,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No PolyPhen predictions in the cache of homo_sapiens GRCh38: pub/release-116/mysql/homo_sapiens_variation_116_38 holds no PolyPhen matrix of a protein of this genome",
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    let info = installed.info();
    assert!(info["sift"] == "b" && info["polyphen"].is_null(), "{info}");
}

/// A tree with no regulation directory installs without the regulatory tree:
/// the run ends on one warning for the regulatory and the motif features
/// together, which names the listing the site lacks, then on the warnings
/// of the cells of the known variants, and names no cell of a regulatory
/// tree; `vep --everything` over the release says the cache holds no
/// regulatory data.
#[test]
fn a_tree_without_a_regulation_directory_installs_without_regulatory_features() {
    let parts = Parts {
        regulation: Regulation::Absent,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No regulatory features or motif features in the cache of homo_sapiens GRCh38: pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS is not on the site",
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            RECOVERED_AF,
        ]
    );
    assert!(!installed.release().join("regulatory").exists());
    let (_, notices) = installed.annotate();
    assert!(
        notices.contains("this cache holds no regulatory data"),
        "{notices}"
    );
}

/// An annotation directory that lists the regulatory-feature file and no
/// motif-feature file installs the regulatory features alone: the run ends on
/// the warning for the motif features, and of the cells of the regulatory
/// tree it names the activity of the regulatory features alone.
#[test]
fn an_annotation_listing_without_a_motif_file_installs_without_motif_features() {
    let parts = Parts {
        regulation: Regulation::NoMotifFile,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No motif features in the cache of homo_sapiens GRCh38: pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS lists no motif_features GFF3",
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            RECOVERED_AF,
        ]
    );
    let features = installed.records("regulatory");
    assert_eq!(features.len(), 516);
    assert!(features
        .iter()
        .all(|feature| feature["feature_type"] == "RegulatoryFeature"));
}

/// A species other than `homo_sapiens` whose tree publishes its one
/// genome-wide variation VCF, matrices for SIFT alone and no regulation
/// directory installs the variants and the SIFT matrices: the run ends on the
/// warnings for PolyPhen and for the regulatory part, then on the warnings
/// of the cells of the known variants, which name no gnomAD frequency:
/// Ensembl VEP's cache holds those for `homo_sapiens` alone.
#[test]
fn another_species_installs_the_parts_its_tree_publishes() {
    let parts = Parts {
        variation: Variation::Files,
        database: Database::SiftAlone,
        regulation: Regulation::Absent,
    };
    let installed = install(ZEBRAFISH, parts, &[]);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No PolyPhen predictions in the cache of danio_rerio GRCz11: pub/release-116/mysql/danio_rerio_variation_116_11 holds no PolyPhen matrix of a protein of this genome",
            "WARN vep_cache_builder::install: No regulatory features or motif features in the cache of danio_rerio GRCz11: pub/release-116/regulation/danio_rerio/GRCz11/annotation/CHECKSUMS is not on the site",
            "WARN vep_cache_builder::install: No PubMed identifiers or variant synonyms in the cache of danio_rerio GRCz11: pub/release-116/mysql/danio_rerio_variation_116_11/CHECKSUMS lists no variation table; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
            "WARN vep_cache_builder::install: No failed variants in the cache of danio_rerio GRCz11: pub/release-116/variation/gvf/danio_rerio/CHECKSUMS is not on the site; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
        ]
    );
    let release = installed.release();
    assert!(release.join("variations/21/all_vars.gz").is_file());
    let info = installed.info();
    assert!(info["sift"] == "b" && info["polyphen"].is_null(), "{info}");
    assert_eq!(info["regulatory"], false, "{info}");
    assert_eq!(
        info["left_out"],
        serde_json::json!(["pubmed", "var_synonyms", "failed_variants"])
    );
}

/// A sequence of the cache that the installation read no variation record on
/// is named in the run's last lines. The one genome-wide variation VCF of a
/// species other than `homo_sapiens` holds its one record on `MT`, and the
/// model names transcripts on chromosome 21 alone: the cache carries the
/// variation file of `MT` and none of chromosome 21, the run ends on a
/// warning that names chromosome 21, before the warnings of the cells of
/// the known variants, and `vep --everything` on the release finds the
/// variation data and fills no known-variant cell of a record on
/// chromosome 21.
#[test]
fn a_sequence_the_variation_file_holds_no_record_on_is_named_in_the_last_lines() {
    let parts = Parts {
        variation: Variation::OffTheModel,
        database: Database::SiftAlone,
        regulation: Regulation::Absent,
    };
    let installed = install(ZEBRAFISH, parts, &[]);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No PolyPhen predictions in the cache of danio_rerio GRCz11: pub/release-116/mysql/danio_rerio_variation_116_11 holds no PolyPhen matrix of a protein of this genome",
            "WARN vep_cache_builder::install: No regulatory features or motif features in the cache of danio_rerio GRCz11: pub/release-116/regulation/danio_rerio/GRCz11/annotation/CHECKSUMS is not on the site",
            "WARN vep_cache_builder::install: No known variants on sequence 21 in the cache of danio_rerio GRCz11: the installation read no variation record there, in a variation VCF or in the variation_feature table; a conversion of Ensembl VEP's cache carries any known variant that cache holds there (docs/cache-conversion.md)",
            "WARN vep_cache_builder::install: No PubMed identifiers or variant synonyms in the cache of danio_rerio GRCz11: pub/release-116/mysql/danio_rerio_variation_116_11/CHECKSUMS lists no variation table; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
            "WARN vep_cache_builder::install: No failed variants in the cache of danio_rerio GRCz11: pub/release-116/variation/gvf/danio_rerio/CHECKSUMS is not on the site; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
        ]
    );
    let release = installed.release();
    assert!(release.join("variations/MT/all_vars.gz").is_file());
    assert!(!release.join("variations/21").exists());
    assert_eq!(installed.info()["contig_order"], serde_json::json!(["21"]));

    let (rows, notices) = installed.annotate();
    assert!(!notices.contains("no variation data"), "{notices}");
    let known: Vec<&Vec<String>> = rows
        .iter()
        .filter(|row| row[UPLOADED_VARIATION] == "known")
        .collect();
    assert!(!known.is_empty());
    for row in known {
        assert_eq!(row[EXISTING_VARIATION], "-", "{row:?}");
    }
}

/// Under `homo_sapiens` the sequences read are those whose per-chromosome
/// files the installation fetched: with `--chromosomes MT` on a model of
/// chromosome 21 the cache carries the variation file of `MT` alone, and the
/// run ends on the line that the cache holds each part, the warning that
/// names chromosome 21, the warnings of the cells and the warning of the
/// recovered `AF`.
#[test]
fn a_chromosome_whose_file_was_not_fetched_is_named_in_the_last_lines() {
    let parts = Parts {
        variation: Variation::OffTheModel,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &["--chromosomes", "MT"]);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            "WARN vep_cache_builder::install: No known variants on sequence 21 in the cache of homo_sapiens GRCh38: the installation read no variation record there, in a variation VCF or in the variation_feature table; a conversion of Ensembl VEP's cache carries any known variant that cache holds there (docs/cache-conversion.md)",
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    let written: Vec<String> = std::fs::read_dir(installed.release().join("variations"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    assert_eq!(written, ["MT"]);
}

/// An assembly whose regulation directory is named for the last
/// dot-separated part of its name installs the regulatory tree from that
/// directory: chicken's `bGalGal1.mat.broiler.GRCg7b` stands under `GRCg7b`,
/// with files named for `GRCg7b`.
#[test]
fn a_regulation_directory_named_for_the_last_part_of_the_assembly_is_installed() {
    let parts = Parts {
        variation: Variation::Absent,
        database: Database::Absent,
        regulation: Regulation::EveryFile,
    };
    let installed = install(CHICKEN, parts, &[]);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            "WARN vep_cache_builder::install: No known variants in the cache of gallus_gallus bGalGal1.mat.broiler.GRCg7b: pub/release-116/variation/vcf/gallus_gallus/CHECKSUMS is not on the site",
            "WARN vep_cache_builder::install: No SIFT or PolyPhen predictions in the cache of gallus_gallus bGalGal1.mat.broiler.GRCg7b: no variation database is published for this release and assembly",
            "WARN vep_cache_builder::install: No cell-type activity of regulatory features in the cache of gallus_gallus bGalGal1.mat.broiler.GRCg7b: pub/release-116/regulation/gallus_gallus/GRCg7b/annotation/CHECKSUMS lists no regulatory_activity file; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
            "WARN vep_cache_builder::install: No HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features in the cache of gallus_gallus bGalGal1.mat.broiler.GRCg7b: pub/release-116/mysql/gallus_gallus_funcgen_116_7/CHECKSUMS is not on the site; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
        ]
    );
    let info = installed.info();
    assert_eq!(info["regulatory"], true, "{info}");
    assert_eq!(installed.records("regulatory").len(), 516 + 2200);
}

/// The table of matrices of an installed release, by stable identifier.
fn binding_matrices(release: &Path) -> BTreeMap<String, Value> {
    serde_json::from_slice(
        &std::fs::read(release.join("regulatory/binding_matrices.json")).unwrap(),
    )
    .unwrap()
}

/// A GRCh38 tree that publishes the activity file beside its regulatory GFF3
/// files, and the funcgen database beside its core database, installs what a
/// conversion of Ensembl VEP's own cache carries of them. The record lists
/// the file's cell types, each space an underscore, in byte order, and the
/// regulatory build's version. Every regulatory feature carries one code a
/// cell type in that order: the codes of its row, and `-` for each where the
/// file holds no row of it, as for the 193 CTCF binding sites. Every motif
/// carries the bases of its span in the genome on its own strand. Every
/// matrix a motif names carries the table's counts, name, source and
/// threshold, and the names of its complexes in the order of the complexes'
/// identifiers; the matrix no motif names is not written. The run's last
/// lines name no cell of the regulatory tree, and the record names neither
/// under `left_out`.
///
/// `vep --cell_type` over the installed release then prints the enhancer's
/// activity in the cell types asked for, in the order asked, and the
/// regulatory build's header line. The made-up matrix of each scored motif
/// holds, at each odd position, 1000 for the base the motif has there and
/// none for another, and 250 for every base at each even position. A
/// substitution at an odd position is therefore at an informative position
/// and lowers the motif's score, which is the best the matrix gives, by one
/// part in as many as the matrix has odd positions: an eighth for the motif
/// of 15 bases, a seventh for the motif of 13. A substitution at an even
/// position is at no informative position and changes the score by nothing.
/// The reverse-strand motif counts its positions from its last base on the
/// genome.
#[test]
fn a_tree_with_the_activity_file_and_the_funcgen_database_installs_their_cells() {
    let parts = Parts {
        regulation: Regulation::ActivityAndTables,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    let info = installed.info();
    assert_eq!(
        info["cell_types"],
        serde_json::json!(LISTED_CELL_TYPES),
        "{info}"
    );
    assert_eq!(
        info["source_versions"]["regbuild"], REGULATORY_BUILD,
        "{info}"
    );
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            RECOVERED_AF,
        ]
    );
    assert_eq!(
        info["left_out"],
        serde_json::json!([
            "panel_frequency",
            "gnomad_frequencies",
            "pubmed",
            "var_synonyms",
            "failed_variants"
        ])
    );

    let features = installed.records("regulatory");
    let with_activity = |codes: &str| {
        features
            .iter()
            .filter(|feature| feature["feature_type"] == "RegulatoryFeature")
            .filter(|feature| feature["activity"] == codes)
            .count()
    };
    assert_eq!(
        (
            with_activity("INPA"),
            with_activity("----"),
            with_activity("IRII")
        ),
        (1, 193, 516 - 194)
    );
    let enhancer = features
        .iter()
        .find(|feature| feature["stable_id"] == ENHANCER)
        .unwrap();
    assert_eq!(enhancer["activity"], "INPA");

    let motifs: Vec<&Value> = features
        .iter()
        .filter(|feature| feature["feature_type"] == "MotifFeature")
        .collect();
    assert_eq!(motifs.len(), 2200);
    let span =
        |motif: &Value| motif["end"].as_u64().unwrap() - motif["start"].as_u64().unwrap() + 1;
    assert!(motifs.iter().all(|motif| motif["seq"]
        .as_str()
        .is_some_and(|seq| seq.len() as u64 == span(motif))));
    for scored in [&FORWARD_MOTIF, &REVERSE_MOTIF] {
        let motif = motifs
            .iter()
            .find(|motif| motif["stable_id"] == scored.id)
            .unwrap();
        // Compared without printing either: a sequence is no part of a log.
        assert!(
            motif["seq"].as_str().map(str::as_bytes) == Some(scored.sequence().as_slice()),
            "the sequence of {} is not that of its span on its own strand",
            scored.id
        );
    }

    let table = binding_matrices(&installed.release());
    let made_up = made_up_matrices();
    assert_eq!((table.len(), made_up.len()), (50, 50));
    for (number, (stable_id, counts)) in made_up.iter().enumerate() {
        let factors = if stable_id == FORWARD_MOTIF.matrix {
            serde_json::json!(["FactorB::FactorC", "FactorA"])
        } else {
            serde_json::json!([])
        };
        assert_eq!(
            table[stable_id],
            serde_json::json!({
                "frequencies": counts,
                "length": counts.len(),
                "name": format!("made-up matrix {}", 101 + number),
                "source": "SELEX",
                "threshold": 4.4,
                "transcription_factors": factors
            }),
            "{stable_id}"
        );
    }

    let records = [
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string(),
        "21\t5011400\tin_enhancer\tA\tG\t.\t.\t.\n".to_string(),
        FORWARD_MOTIF.record(1, "forward_odd"),
        FORWARD_MOTIF.record(2, "forward_even"),
        REVERSE_MOTIF.record(13, "reverse_odd"),
        REVERSE_MOTIF.record(12, "reverse_even"),
    ]
    .concat();
    let asked = "tissue_two_(m,_4_y),cell_line_one,organ_three_(f,_7_dpf)";
    let run = installed.run_vep(&records, &["--regulatory", "--cell_type", asked]);
    let regbuild = format!("## regbuild version {REGULATORY_BUILD}");
    assert!(run.header.contains(&regbuild), "{:?}", run.header);
    let cells_of = |record: &str, feature: &str| -> Vec<String> {
        let row = run
            .rows
            .iter()
            .find(|row| row[UPLOADED_VARIATION] == record && row[FEATURE] == feature)
            .unwrap_or_else(|| panic!("no row of {feature} for {record}: {}", run.notices));
        row[EXTRA].split(';').map(str::to_string).collect()
    };
    let holds = |cells: &[String], cell: &str| cells.iter().any(|held| held == cell);
    let enhancer = cells_of("in_enhancer", ENHANCER);
    assert!(
        holds(
            &enhancer,
            "CELL_TYPE=tissue_two_(m,_4_y):ACTIVE,cell_line_one:INACTIVE,organ_three_(f,_7_dpf):NA"
        ),
        "{enhancer:?}"
    );
    for (record, motif, expected) in [
        (
            "forward_odd",
            &FORWARD_MOTIF,
            ["HIGH_INF_POS=Y", "MOTIF_SCORE_CHANGE=-0.125"],
        ),
        (
            "forward_even",
            &FORWARD_MOTIF,
            ["HIGH_INF_POS=N", "MOTIF_SCORE_CHANGE=0.000"],
        ),
        (
            "reverse_odd",
            &REVERSE_MOTIF,
            ["HIGH_INF_POS=Y", "MOTIF_SCORE_CHANGE=-0.143"],
        ),
        (
            "reverse_even",
            &REVERSE_MOTIF,
            ["HIGH_INF_POS=N", "MOTIF_SCORE_CHANGE=0.000"],
        ),
    ] {
        let cells = cells_of(record, motif.id);
        for cell in expected {
            assert!(holds(&cells, cell), "{record}: {cell}: {cells:?}");
        }
        assert!(
            holds(&cells, &format!("MOTIF_NAME={}", motif.matrix)),
            "{cells:?}"
        );
    }
    let factors = cells_of("forward_odd", FORWARD_MOTIF.id);
    assert!(
        holds(&factors, "TRANSCRIPTION_FACTORS=FactorB::FactorC,FactorA"),
        "{factors:?}"
    );
}

/// A tree that publishes one of the two and not the other installs, with no
/// error, the part it publishes and leaves the other as a tree with neither
/// installs it. Without the activity file the shards are those of the tree
/// with neither, each regulatory feature with an empty activity string, and
/// the record lists no cell type, while the matrices carry their counts.
/// Without the funcgen database the table of matrices is that of the tree
/// with neither, each matrix its stub, and the record names no regulatory
/// build, while the record lists the cell types. Each run says in one line
/// what it did not find, ends on the warning of the one cell it did not
/// bring, with the listing that lacks its file or table, and names that
/// cell alone of the two in its record.
#[test]
fn a_tree_without_the_activity_file_or_the_funcgen_database_installs_the_other() {
    let neither = install(HUMAN, EVERY_PART, &CHROMOSOME_21);
    let of_tree =
        |files: &BTreeMap<PathBuf, Vec<u8>>, shards: bool| -> BTreeMap<PathBuf, Vec<u8>> {
            files
                .iter()
                .filter(|(file, _)| file.starts_with("regulatory"))
                .filter(|(file, _)| (file.components().count() == 3) == shards)
                .map(|(file, bytes)| (file.clone(), bytes.clone()))
                .collect()
        };
    let no_activity = "No cell-type activity is read for the regulatory tree: pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS lists no regulatory_activity file";
    let no_funcgen = "No table of the funcgen database is read for the regulatory tree: pub/release-116/mysql/homo_sapiens_funcgen_116_38/CHECKSUMS is not on the site";
    for line in [no_activity, no_funcgen] {
        assert!(stderr(&neither.output).contains(line), "{line}");
    }
    assert!(binding_matrices(&neither.release())
        .values()
        .all(|matrix| matrix.as_object().unwrap().len() == 2));

    let tables_alone = install(
        HUMAN,
        Parts {
            regulation: Regulation::TablesAlone,
            ..EVERY_PART
        },
        &CHROMOSOME_21,
    );
    let log = stderr(&tables_alone.output);
    assert!(
        log.contains(no_activity) && !log.contains(no_funcgen),
        "{log}"
    );
    let info = tables_alone.info();
    assert!(info.get("cell_types").is_none(), "{info}");
    assert_eq!(
        info["source_versions"]["regbuild"], REGULATORY_BUILD,
        "{info}"
    );
    assert_eq!(of_tree(&tables_alone.files(), true).len(), 6);
    assert!(of_tree(&tables_alone.files(), true) == of_tree(&neither.files(), true));
    assert!(tables_alone
        .records("regulatory")
        .iter()
        .filter(|feature| feature["feature_type"] == "RegulatoryFeature")
        .all(|feature| feature["activity"] == ""));
    let table = binding_matrices(&tables_alone.release());
    assert_eq!(table.len(), 50);
    assert!(table.values().all(|matrix| {
        matrix["frequencies"].as_array().map(Vec::len)
            == matrix["length"].as_u64().map(|n| n as usize)
    }));
    assert_eq!(
        tables_alone.closing_lines(),
        [
            tables_alone.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_ACTIVITY,
            RECOVERED_AF,
        ]
    );
    assert_eq!(
        regulatory_cells_named(&tables_alone),
        ["cell_type_activity"]
    );

    let activity_alone = install(
        HUMAN,
        Parts {
            regulation: Regulation::ActivityAlone,
            ..EVERY_PART
        },
        &CHROMOSOME_21,
    );
    let log = stderr(&activity_alone.output);
    assert!(
        log.contains(no_funcgen) && !log.contains(no_activity),
        "{log}"
    );
    let info = activity_alone.info();
    assert_eq!(
        info["cell_types"],
        serde_json::json!(LISTED_CELL_TYPES),
        "{info}"
    );
    assert!(info["source_versions"].get("regbuild").is_none(), "{info}");
    assert!(of_tree(&activity_alone.files(), false) == of_tree(&neither.files(), false));
    assert_eq!(of_tree(&activity_alone.files(), false).len(), 1);
    assert!(activity_alone
        .records("regulatory")
        .iter()
        .filter(|feature| feature["feature_type"] == "RegulatoryFeature")
        .all(|feature| feature["activity"]
            .as_str()
            .is_some_and(|codes| codes.len() == 4)));
    assert_eq!(
        activity_alone.closing_lines(),
        [
            activity_alone.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_FAILED_VARIANTS,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    assert_eq!(regulatory_cells_named(&activity_alone), ["motif_scores"]);
    assert_eq!(
        regulatory_cells_named(&neither),
        ["cell_type_activity", "motif_scores"]
    );
}

/// The names of `names` that the record of an installed release lists under
/// `left_out`, in the record's order.
fn named_among(installed: &Installation, names: &[&str]) -> Vec<String> {
    installed.info()["left_out"]
        .as_array()
        .map(|listed| {
            listed
                .iter()
                .filter_map(Value::as_str)
                .filter(|name| names.contains(name))
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

/// The cells of the regulatory tree that the record of an installed release
/// names under `left_out`.
fn regulatory_cells_named(installed: &Installation) -> Vec<String> {
    named_among(installed, &["cell_type_activity", "motif_scores"])
}

/// The cells of the known variants that the record of an installed release
/// names under `left_out`.
fn variant_cells_named(installed: &Installation) -> Vec<String> {
    named_among(installed, &["pubmed", "var_synonyms", "failed_variants"])
}

/// A tree whose variation database publishes the made-up tables of the
/// three cells and whose GVF directory publishes the failed variants.
const TABLES_AND_FAILED_GVF: Parts = Parts {
    variation: Variation::TwoRecordsAndFailedGvf,
    database: Database::VariantTables,
    ..EVERY_PART
};

/// The records of the variation file of chromosome 21 of an installed
/// release, each as its cells by column name.
fn known_variants(installed: &Installation) -> Vec<BTreeMap<String, String>> {
    known_variants_on(installed, "21")
}

/// The records of the variation file of `sequence` of an installed release,
/// each as its cells by column name.
fn known_variants_on(installed: &Installation, sequence: &str) -> Vec<BTreeMap<String, String>> {
    let columns: Vec<String> = installed.info()["variation_cols"]
        .as_array()
        .unwrap()
        .iter()
        .map(|column| column.as_str().unwrap().to_string())
        .collect();
    let mut text = String::new();
    bgzf::io::Reader::new(
        std::fs::File::open(
            installed
                .release()
                .join("variations")
                .join(sequence)
                .join("all_vars.gz"),
        )
        .unwrap(),
    )
    .read_to_string(&mut text)
    .unwrap();
    text.lines()
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), columns.len(), "{line}");
            columns
                .iter()
                .cloned()
                .zip(cells.into_iter().map(str::to_string))
                .collect()
        })
        .collect()
}

/// The name, the `failed`, `start`, `end` and `allele_string` cells and the
/// `pubmed` and `var_synonyms` cells of each record of `records`.
fn cells_of(records: &[BTreeMap<String, String>]) -> Vec<[&str; 7]> {
    records
        .iter()
        .map(|record| {
            [
                "variation_name",
                "failed",
                "start",
                "end",
                "allele_string",
                "pubmed",
                "var_synonyms",
            ]
            .map(|column| record[column].as_str())
        })
        .collect()
}

/// The records an annotated run reads the three cells on: the known variant
/// of every tree, [`FAILED_IN_THE_VCF`] and [`FAILED_SUBSTITUTION`], each
/// with the alleles its file gives it.
const KNOWN_AND_FAILED: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    21\t18000300\tknown\tT\tC\t.\t.\t.\n\
    21\t18000350\tin_the_vcf\tG\tA\t.\t.\t.\n\
    21\t18000400\tin_the_gvf\tA\tG\t.\t.\t.\n";

/// The line an installation logs of what it wrote of the three cells, from
/// its module's name on.
fn cells_line(pubmed: &str, synonyms: &str, failed: &str, added: u64) -> String {
    format!(
        "vep_cache_builder::install: The known variants carry PubMed identifiers on {pubmed}, \
         variant synonyms on {synonyms} and the failed flag on {failed}, {added} of them added \
         from the failed-variants GVF"
    )
}

/// A GRCh38 tree whose variation database publishes the tables of the
/// citations, the synonyms and the failed variants, and whose GVF directory
/// publishes the failed variants, installs what a conversion of Ensembl
/// VEP's own cache carries of them. The known variant holds the PubMed
/// identifiers of the two publications that have one, in the order of the
/// publications' identifiers, and its synonyms of the two sources not named
/// for dbSNP, the sources in the order of their names without regard to
/// case and the names of one source from the last row back. The variant the
/// VCF holds and the tables mark as failed is flagged. The two failed
/// variants the VCF leaves out are records of their own, in start order,
/// the insertion one base after the end its feature writes, and the one a
/// publication cites holds that publication's identifier; the feature on a
/// chromosome the installation builds no file of adds nothing. The run's
/// last lines name, of the cells of the known variants, the gnomAD
/// frequencies alone, and the record names none of the three.
#[test]
fn a_tree_with_the_variation_tables_and_the_failed_gvf_installs_their_cells() {
    let installed = install(HUMAN, TABLES_AND_FAILED_GVF, &CHROMOSOME_21);
    let records = known_variants(&installed);
    assert_eq!(
        cells_of(&records),
        [
            [
                "rs18000300",
                ".",
                "18000300",
                ".",
                "T/C",
                KNOWN_PUBMED,
                KNOWN_SYNONYMS
            ],
            [FAILED_IN_THE_VCF, "1", "18000350", ".", "G/A", ".", "."],
            [
                FAILED_SUBSTITUTION,
                "1",
                "18000400",
                ".",
                "A/G",
                "30000021",
                "."
            ],
            [
                FAILED_INSERTION,
                "1",
                "18000451",
                "18000450",
                "-/TT",
                ".",
                "."
            ],
        ]
    );
    assert_eq!(records[0]["AFR"], "C:0.25");
    let log = stderr(&installed.output);
    let line = cells_line("2 records", "1 record", "3 records", 2);
    assert!(log.contains(&line), "{line}: {log}");
    assert!(!log.contains("is read for the known variants"), "{log}");
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    assert_eq!(
        installed.info()["left_out"],
        serde_json::json!([
            "panel_frequency",
            "gnomad_frequencies",
            "cell_type_activity",
            "motif_scores"
        ])
    );
}

/// `vep --everything --var_synonyms` over a release installed from the tree
/// of the made-up tables and the made-up GVF prints the known variant's
/// `PUBMED` and `VAR_SYNONYMS` cells and names neither failed variant, as
/// Ensembl VEP names none without `--failed 1`; with that flag the run names
/// both, the one the GVF brought with its `PUBMED` cell.
#[test]
fn a_run_over_the_installed_cells_prints_them_and_keeps_a_failed_variant_on_request() {
    let installed = install(HUMAN, TABLES_AND_FAILED_GVF, &CHROMOSOME_21);
    let existing = |run: &Run, record: &str| -> Vec<(String, Vec<String>)> {
        let rows: Vec<&Vec<String>> = run
            .rows
            .iter()
            .filter(|row| row[UPLOADED_VARIATION] == record)
            .collect();
        assert!(!rows.is_empty(), "no row of {record}: {}", run.notices);
        rows.into_iter()
            .map(|row| {
                let cells = row[EXTRA]
                    .split(';')
                    .filter(|cell| cell.starts_with("PUBMED=") || cell.starts_with("VAR_SYNONYMS="))
                    .map(str::to_string)
                    .collect();
                (row[EXISTING_VARIATION].clone(), cells)
            })
            .collect()
    };
    let known_cells = vec![
        format!("PUBMED={KNOWN_PUBMED}"),
        format!("VAR_SYNONYMS={KNOWN_SYNONYMS}"),
    ];
    let named = |name: &str, cells: &[String]| (name.to_string(), cells.to_vec());

    let default = installed.run_vep(KNOWN_AND_FAILED, &["--everything", "--var_synonyms"]);
    for row in existing(&default, "known") {
        assert_eq!(row, named("rs18000300", &known_cells));
    }
    for record in ["in_the_vcf", "in_the_gvf"] {
        for row in existing(&default, record) {
            assert_eq!(row, named("-", &[]), "{record}");
        }
    }

    let with_failed = installed.run_vep(
        KNOWN_AND_FAILED,
        &["--everything", "--var_synonyms", "--failed", "1"],
    );
    for row in existing(&with_failed, "known") {
        assert_eq!(row, named("rs18000300", &known_cells));
    }
    for row in existing(&with_failed, "in_the_vcf") {
        assert_eq!(row, named(FAILED_IN_THE_VCF, &[]));
    }
    for row in existing(&with_failed, "in_the_gvf") {
        assert_eq!(
            row,
            named(FAILED_SUBSTITUTION, &["PUBMED=30000021".to_string()])
        );
    }
}

/// A tree whose variation database publishes every table of the cells but
/// `variation`, as the GRCh37 tree does, installs: no row of the other
/// tables can be given its variant's name, so the known variant holds no
/// PubMed identifier and no synonym, and the run says so in one line a cell,
/// naming the table the listing lacks. The failed variants still come from
/// the GVF: the one the VCF holds is flagged and the two it leaves out are
/// records. The run's last lines name the PubMed identifiers and the
/// synonyms in one warning, with the listing that lacks the table, and no
/// failed variant; the record names those two cells and not the third.
#[test]
fn a_tree_without_the_variation_table_installs_the_failed_variants_alone() {
    let parts = Parts {
        variation: Variation::TwoRecordsAndFailedGvf,
        database: Database::NoVariationTable,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        cells_of(&known_variants(&installed)),
        [
            ["rs18000300", ".", "18000300", ".", "T/C", ".", "."],
            [FAILED_IN_THE_VCF, "1", "18000350", ".", "G/A", ".", "."],
            [FAILED_SUBSTITUTION, "1", "18000400", ".", "A/G", ".", "."],
            [
                FAILED_INSERTION,
                "1",
                "18000451",
                "18000450",
                "-/TT",
                ".",
                "."
            ],
        ]
    );
    let log = stderr(&installed.output);
    let no_table = "is read for the known variants: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no variation table";
    for what in [
        "PubMed identifier",
        "variant synonym",
        "failed_variation table",
    ] {
        let line = format!("No {what} {no_table}");
        assert!(log.contains(&line), "{line}: {log}");
    }
    assert!(!log.contains("No failed-variants GVF"), "{log}");
    let line = cells_line("0 records", "0 records", "3 records", 2);
    assert!(log.contains(&line), "{line}: {log}");
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_TABLE_CELLS,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    assert_eq!(
        installed.info()["left_out"],
        serde_json::json!([
            "panel_frequency",
            "gnomad_frequencies",
            "pubmed",
            "var_synonyms",
            "cell_type_activity",
            "motif_scores"
        ])
    );
}

/// A tree that publishes the tables and no GVF directory installs: the
/// known variant holds its PubMed identifiers and its synonyms, the variant
/// the VCF holds is flagged as failed from `failed_variation`, and no record
/// is added for a failed variant the VCF leaves out, of which the run says
/// in one line that the GVF is not on the site. The flag was read, so the
/// run's last lines and the record name none of the three cells.
#[test]
fn a_tree_without_the_failed_gvf_installs_the_cells_of_the_tables() {
    let parts = Parts {
        variation: Variation::TwoRecords,
        database: Database::VariantTables,
        ..EVERY_PART
    };
    let installed = install(HUMAN, parts, &CHROMOSOME_21);
    assert_eq!(
        cells_of(&known_variants(&installed)),
        [
            [
                "rs18000300",
                ".",
                "18000300",
                ".",
                "T/C",
                KNOWN_PUBMED,
                KNOWN_SYNONYMS
            ],
            [FAILED_IN_THE_VCF, "1", "18000350", ".", "G/A", ".", "."],
        ]
    );
    let log = stderr(&installed.output);
    let no_gvf = "No failed-variants GVF is read for the known variants: pub/release-116/variation/gvf/homo_sapiens/CHECKSUMS is not on the site";
    assert!(log.contains(no_gvf), "{log}");
    assert!(!log.contains("lists no variation table"), "{log}");
    let line = cells_line("1 record", "1 record", "1 record", 0);
    assert!(log.contains(&line), "{line}: {log}");
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            NO_ACTIVITY,
            NO_MOTIF_SCORES,
            RECOVERED_AF,
        ]
    );
    assert!(variant_cells_named(&installed).is_empty());
}

/// A tree that publishes none of the tables and no GVF directory installs
/// its known variants as their VCFs give them, each of the three cells
/// empty, says of each source in one line that the site lacks it, names the
/// three cells in its last lines and lists them in its record.
#[test]
fn a_tree_without_the_tables_or_the_failed_gvf_installs_the_three_cells_empty() {
    let installed = install(HUMAN, EVERY_PART, &CHROMOSOME_21);
    assert_eq!(
        cells_of(&known_variants(&installed)),
        [["rs18000300", ".", "18000300", ".", "T/C", ".", "."]]
    );
    let log = stderr(&installed.output);
    for line in [
        "No PubMed identifier is read for the known variants: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no variation table",
        "No variant synonym is read for the known variants: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no variation table",
        "No failed_variation table is read for the known variants: pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no variation table",
        "No failed-variants GVF is read for the known variants: pub/release-116/variation/gvf/homo_sapiens/CHECKSUMS is not on the site",
    ] {
        assert!(log.contains(line), "{line}: {log}");
    }
    assert!(!log.contains("The known variants carry"), "{log}");
    let closing = installed.closing_lines();
    for line in [NO_TABLE_CELLS, NO_FAILED_VARIANTS] {
        assert_eq!(
            closing.iter().filter(|closed| *closed == line).count(),
            1,
            "{line}: {closing:#?}"
        );
    }
    assert_eq!(
        variant_cells_named(&installed),
        ["pubmed", "var_synonyms", "failed_variants"]
    );
}

/// A tree that publishes every file a cell is read from: the tables of the
/// citations, the synonyms and the failed variants, the GVF file of the
/// failed variants, the activity file and the funcgen database.
const EVERY_FILE_OF_THE_CELLS: Parts = Parts {
    variation: Variation::TwoRecordsAndFailedGvf,
    database: Database::VariantTables,
    regulation: Regulation::ActivityAndTables,
};

/// The five cells an installation brings where the tree publishes their
/// files, as a record names them under `left_out`.
const FIVE_CELLS: [&str; 5] = [
    "pubmed",
    "var_synonyms",
    "failed_variants",
    "cell_type_activity",
    "motif_scores",
];

/// The record of an installation names under `left_out` the cells it did
/// not bring. From the tree that publishes every file of the five cells it
/// names none of them: what is left is the recovered `AF` and the gnomAD
/// frequencies, which no file of the release carries. From the tree that
/// publishes none of those files it names all five.
#[test]
fn the_record_names_the_cells_an_installation_did_not_bring() {
    let every = install(HUMAN, EVERY_FILE_OF_THE_CELLS, &CHROMOSOME_21);
    assert_eq!(
        every.info()["left_out"],
        serde_json::json!(["panel_frequency", "gnomad_frequencies"])
    );
    assert!(named_among(&every, &FIVE_CELLS).is_empty());

    let none = install(HUMAN, EVERY_PART, &CHROMOSOME_21);
    assert_eq!(named_among(&none, &FIVE_CELLS), FIVE_CELLS);
}

/// An installation from the tree that publishes every file of the five
/// cells ends on the line that the cache holds each part, the warning of
/// the gnomAD frequencies and the warning of the recovered `AF`: it names no
/// cell it brought. `vep --everything --var_synonyms` over the release then
/// prints one line for the gnomAD columns and none that names `PUBMED`,
/// `VAR_SYNONYMS`, `HIGH_INF_POS` or `MOTIF_SCORE_CHANGE`, and fills each of
/// those cells; the same run over a release installed from the tree that
/// publishes none of those files names the four columns in that line.
#[test]
fn a_run_over_a_cache_that_carries_every_cell_prints_no_line_for_them() {
    let installed = install(HUMAN, EVERY_FILE_OF_THE_CELLS, &CHROMOSOME_21);
    assert_eq!(
        installed.closing_lines(),
        [
            installed.installed_line().as_str(),
            EVERY_PART_LINE,
            NO_GNOMAD,
            RECOVERED_AF,
        ]
    );

    let records = [
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string(),
        FORWARD_MOTIF.record(1, "forward_odd"),
        "21\t18000300\tknown\tT\tC\t.\t.\t.\n".to_string(),
    ]
    .concat();
    let flags = ["--everything", "--var_synonyms"];
    let not_carried = |run: &Run| -> Vec<String> {
        run.notices
            .lines()
            .filter(|line| line.contains("which this cache does not carry"))
            .map(str::to_string)
            .collect()
    };
    let columns = [
        "PUBMED",
        "VAR_SYNONYMS",
        "HIGH_INF_POS",
        "MOTIF_SCORE_CHANGE",
    ];

    let run = installed.run_vep(&records, &flags);
    assert_eq!(
        not_carried(&run),
        ["INFO: --everything asks for the gnomAD exome columns and the gnomAD genome columns, which this cache does not carry: the columns stay empty; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)"],
        "{}",
        run.notices
    );
    for column in columns {
        assert!(
            !run.notices.lines().any(|line| line.contains(column)),
            "{column}: {}",
            run.notices
        );
    }
    let cells_of = |record: &str, feature_type: &str| -> Vec<String> {
        let row = run
            .rows
            .iter()
            .find(|row| row[UPLOADED_VARIATION] == record && row[FEATURE_TYPE] == feature_type)
            .unwrap_or_else(|| panic!("no {feature_type} row of {record}: {}", run.notices));
        row[EXTRA].split(';').map(str::to_string).collect()
    };
    let known = cells_of("known", "Transcript");
    for cell in [
        format!("PUBMED={KNOWN_PUBMED}"),
        format!("VAR_SYNONYMS={KNOWN_SYNONYMS}"),
    ] {
        assert!(known.contains(&cell), "{cell}: {known:?}");
    }
    let motif = cells_of("forward_odd", "MotifFeature");
    for cell in ["HIGH_INF_POS=Y", "MOTIF_SCORE_CHANGE=-0.125"] {
        assert!(motif.iter().any(|held| held == cell), "{cell}: {motif:?}");
    }

    let without = install(HUMAN, EVERY_PART, &CHROMOSOME_21).run_vep(&records, &flags);
    assert_eq!(
        not_carried(&without),
        ["INFO: --everything and --var_synonyms ask for the gnomAD exome columns, the gnomAD genome columns, PUBMED, VAR_SYNONYMS, HIGH_INF_POS and MOTIF_SCORE_CHANGE, which this cache does not carry: the columns stay empty; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)"],
        "{}",
        without.notices
    );
}

/// A tree of a species other than `homo_sapiens` whose one genome-wide
/// variation VCF holds its one record on `MT`, whose model names transcripts
/// on chromosome 21 alone, and whose variation database publishes the
/// made-up `variation_feature` table.
const FEATURE_TABLE_TREE: Parts = Parts {
    variation: Variation::OffTheModel,
    database: Database::FeatureTable,
    regulation: Regulation::Absent,
};

/// The line an installation logs of what the `variation_feature` table gave
/// it, from its module's name on.
fn table_line(variants: &str, sequences: &str) -> String {
    format!(
        "vep_cache_builder::install: The variation_feature table gave {variants} on \
         {sequences}, which the variation VCF holds no record on"
    )
}

/// The line an installation logs where it fetches no `variation_feature`
/// table, from its module's name on; `named` where `--chromosomes` names
/// the sequences read.
fn table_not_fetched_line(named: bool) -> String {
    format!(
        "vep_cache_builder::install: The variation_feature table is not fetched: the variation \
         VCF holds a record on every sequence of the cache{}",
        if named {
            " that --chromosomes names"
        } else {
            ""
        }
    )
}

/// A sequence of the cache that the one genome-wide variation VCF holds no
/// record on takes its known variants from the variation database's
/// `variation_feature` table. The cache carries a variation file of
/// chromosome 21 with the table's two rows on display there, in start
/// order: the variant with its clinical significance, each space an
/// underscore, and its phenotype-or-disease flag, and the insertion with
/// its start past its end; the row that is not on display is no record.
/// `MT`, which the VCF holds a record on, keeps that record alone: the
/// table's row there is not read. The run's last lines name no sequence,
/// and `vep --everything` on the release names the table's variant on the
/// record at its place.
#[test]
fn a_sequence_the_vcf_holds_no_record_on_takes_its_known_variants_from_the_feature_table() {
    let installed = install(ZEBRAFISH, FEATURE_TABLE_TREE, &[]);
    let on_21 = known_variants(&installed);
    assert_eq!(
        cells_of(&on_21),
        [
            [TABLE_VARIANT, ".", "18000300", ".", "T/C", ".", "."],
            [
                TABLE_INSERTION,
                ".",
                "18000451",
                "18000450",
                "-/TT",
                ".",
                "."
            ],
        ]
    );
    assert_eq!(
        (
            on_21[0]["clin_sig"].as_str(),
            on_21[0]["phenotype_or_disease"].as_str(),
            on_21[1]["clin_sig"].as_str(),
            on_21[1]["phenotype_or_disease"].as_str()
        ),
        ("likely_pathogenic,pathogenic", "1", ".", ".")
    );
    assert_eq!(
        cells_of(&known_variants_on(&installed, "MT")),
        [["rs18000300", ".", "18000300", ".", "T/C", ".", "."]]
    );

    let log = stderr(&installed.output);
    let line = table_line("2 known variants", "sequence 21");
    assert!(log.contains(&line), "{line}: {log}");
    let closing = installed.closing_lines();
    assert!(
        !closing
            .iter()
            .any(|line| line.contains("No known variants on")),
        "{closing:#?}"
    );

    let (rows, notices) = installed.annotate();
    assert!(!notices.contains("no variation data"), "{notices}");
    let known: Vec<&Vec<String>> = rows
        .iter()
        .filter(|row| row[UPLOADED_VARIATION] == "known")
        .collect();
    assert!(!known.is_empty());
    for row in known {
        assert_eq!(row[EXISTING_VARIATION], TABLE_VARIANT, "{row:?}");
    }
}

/// Under `--chromosomes` only the sequences named are in question. With
/// `MT` named, chromosome 21 is not read for, the table is not fetched, and
/// the run's last lines name chromosome 21. With chromosome 21 named, which
/// the VCF holds no record on, its variation file is written from the table
/// and no file of `MT` is written.
#[test]
fn only_a_named_sequence_takes_its_known_variants_from_the_feature_table() {
    let mt = install(ZEBRAFISH, FEATURE_TABLE_TREE, &["--chromosomes", "MT"]);
    let log = stderr(&mt.output);
    let line = table_not_fetched_line(true);
    assert!(log.contains(&line), "{line}: {log}");
    assert!(!log.contains("variation_feature.txt.gz"), "{log}");
    assert!(!mt.release().join("variations/21").exists());
    let named = "WARN vep_cache_builder::install: No known variants on sequence 21 in the cache of danio_rerio GRCz11: the installation read no variation record there, in a variation VCF or in the variation_feature table; a conversion of Ensembl VEP's cache carries any known variant that cache holds there (docs/cache-conversion.md)";
    assert!(
        mt.closing_lines().iter().any(|line| line == named),
        "{:#?}",
        mt.closing_lines()
    );

    let chromosome = install(ZEBRAFISH, FEATURE_TABLE_TREE, &["--chromosomes", "21"]);
    assert_eq!(
        cells_of(&known_variants(&chromosome)),
        [
            [TABLE_VARIANT, ".", "18000300", ".", "T/C", ".", "."],
            [
                TABLE_INSERTION,
                ".",
                "18000451",
                "18000450",
                "-/TT",
                ".",
                "."
            ],
        ]
    );
    assert!(!chromosome.release().join("variations/MT").exists());
    let log = stderr(&chromosome.output);
    let line = table_line("2 known variants", "sequence 21");
    assert!(log.contains(&line), "{line}: {log}");
    assert!(
        !chromosome
            .closing_lines()
            .iter()
            .any(|line| line.contains("No known variants on")),
        "{:#?}",
        chromosome.closing_lines()
    );
}

/// The table is fetched only for a sequence in question, and read only
/// where the variation files are one genome-wide file. A tree whose VCF
/// holds a record on the cache's one sequence installs that record alone,
/// though the table holds two more there, and does not fetch the table. A
/// tree of `homo_sapiens`, whose files are one a chromosome, reads nothing
/// of the table for a chromosome whose file it did not fetch: the run's
/// last lines name that chromosome, and the run neither fetches the table
/// nor says anything of reading it.
/// A tree whose database publishes no such table says so in one line and
/// names the sequence in its last lines.
#[test]
fn the_feature_table_is_fetched_only_for_a_sequence_one_genome_wide_vcf_leaves_without_a_record() {
    let covered = install(
        ZEBRAFISH,
        Parts {
            variation: Variation::Files,
            ..FEATURE_TABLE_TREE
        },
        &[],
    );
    assert_eq!(
        cells_of(&known_variants(&covered)),
        [["rs18000300", ".", "18000300", ".", "T/C", ".", "."]]
    );
    let log = stderr(&covered.output);
    let line = table_not_fetched_line(false);
    assert!(log.contains(&line), "{line}: {log}");
    assert!(!log.contains("variation_feature.txt.gz"), "{log}");

    let human = install(
        HUMAN,
        Parts {
            variation: Variation::OffTheModel,
            database: Database::FeatureTable,
            ..EVERY_PART
        },
        &["--chromosomes", "MT"],
    );
    let log = stderr(&human.output);
    for words in [
        "variation_feature.txt",
        "The variation_feature table",
        "No variation_feature table is read",
    ] {
        assert!(!log.contains(words), "{words}: {log}");
    }
    assert!(!human.release().join("variations/21").exists());
    assert!(
        human
            .closing_lines()
            .iter()
            .any(|line| line.contains("No known variants on sequence 21 in the cache of homo_sapiens GRCh38: the installation read no variation record there, in a variation VCF or in the variation_feature table;")),
        "{:#?}",
        human.closing_lines()
    );

    let not_listed = install(
        ZEBRAFISH,
        Parts {
            database: Database::EveryTable,
            ..FEATURE_TABLE_TREE
        },
        &[],
    );
    let log = stderr(&not_listed.output);
    let reason = "No variation_feature table is read for a sequence the variation VCF holds no record on: pub/release-116/mysql/danio_rerio_variation_116_11/CHECKSUMS lists no variation_feature table";
    assert!(log.contains(reason), "{log}");
    assert!(!not_listed.release().join("variations/21").exists());
    assert!(
        not_listed
            .closing_lines()
            .iter()
            .any(|line| line.contains("No known variants on sequence 21 in the cache of danio_rerio GRCz11: the installation read no variation record there, in a variation VCF or in the variation_feature table;")),
        "{:#?}",
        not_listed.closing_lines()
    );
}

/// The failed variants of a sequence the VCF holds no record on join the
/// records the table gives it. With the GVF file of the failed variants
/// published beside the table, chromosome 21's variation file holds the
/// table's two records and the GVF's three features there, each flagged
/// failed, in start order and, at one start, in name order: the failed
/// insertion before the table's at the start they share.
#[test]
fn the_failed_variants_of_a_sequence_join_the_records_the_feature_table_gives_it() {
    let installed = install(
        ZEBRAFISH,
        Parts {
            variation: Variation::OffTheModelAndFailedGvf,
            ..FEATURE_TABLE_TREE
        },
        &[],
    );
    assert_eq!(
        cells_of(&known_variants(&installed)),
        [
            [TABLE_VARIANT, ".", "18000300", ".", "T/C", ".", "."],
            [FAILED_IN_THE_VCF, "1", "18000350", ".", "G/A", ".", "."],
            [FAILED_SUBSTITUTION, "1", "18000400", ".", "A/G", ".", "."],
            [
                FAILED_INSERTION,
                "1",
                "18000451",
                "18000450",
                "-/TT",
                ".",
                "."
            ],
            [
                TABLE_INSERTION,
                ".",
                "18000451",
                "18000450",
                "-/TT",
                ".",
                "."
            ],
        ]
    );
    let log = stderr(&installed.output);
    let line = cells_line("0 records", "0 records", "3 records", 3);
    assert!(log.contains(&line), "{line}: {log}");
    let line = table_line("5 known variants", "sequence 21");
    assert!(log.contains(&line), "{line}: {log}");
    // The GVF brought the failed variants: of the three cells the last
    // lines and the record name the two the tables would give.
    let closing = installed.closing_lines();
    assert_eq!(
        closing
            .iter()
            .filter(|line| line.contains("a conversion of Ensembl VEP's cache carries them"))
            .collect::<Vec<_>>(),
        ["WARN vep_cache_builder::install: No PubMed identifiers or variant synonyms in the cache of danio_rerio GRCz11: pub/release-116/mysql/danio_rerio_variation_116_11/CHECKSUMS lists no variation table; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)"],
        "{closing:#?}"
    );
    assert_eq!(variant_cells_named(&installed), ["pubmed", "var_synonyms"]);
}
