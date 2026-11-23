// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Records at the first and the last base of a sequence, past its last base,
//! at the origin of a sequence its cache declares circular and on a transcript
//! stored across that origin, through the `vep` binary. The cache is written
//! to a temporary directory in the layout `tests/fixtures/species/README.md`
//! describes: made-up sequences that each hold one single-exon coding
//! transcript whose coding sequence is the whole transcript, the shape of a
//! bacterial gene model.
//!
//! Every expected row but those of a transcript stored across an origin is
//! the row Ensembl VEP 116.3 writes for the same record, offline under default
//! flags, on a Storable cache holding the same gene models; it writes the same
//! rows with and without `--fasta`, and no warning.
//! Ensembl VEP reads a sequence as a straight line from base 1. `validate_vf`
//! tests a record's start against its end and never against a length
//! (`Parser.pm:483-488`); the transcripts a record reaches are those within
//! 5,000 bases of its span (`AnnotationSource.pm:190`, `:306-309`); the
//! Location column is the span as the record gives it, an insertion's two
//! flanking positions in ascending order (`Utils.pm:141-163`); and
//! `Slice::is_circular` answers 0 for a slice without an adaptor, which every
//! slice of an offline run is (ensembl core `Slice.pm:637-638`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::json;
use vep_core::cache_manifest::CacheManifest;

mod common;

/// A record as its input line, and the rows written for it.
type Case = (&'static str, &'static [&'static str]);

/// The codons of a made-up coding sequence between its start and stop codons,
/// taken in turn.
const BODY: [&str; 19] = [
    "GCT", "GAA", "AAA", "CTG", "ACC", "GGT", "CGT", "TCT", "GTT", "GAT", "CAG", "AAC", "CCG",
    "ATC", "TTC", "CAT", "TAC", "TGG", "TGT",
];

/// The residue each codon of [`BODY`] codes for.
const RESIDUES: &[u8; 19] = b"AEKLTGRSVDQNPIFHYWC";

/// The bases between and around the genes, repeated.
const FILLER: &[u8] = b"ACGTTGCA";

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A coding sequence of `codons` codons: `ATG`, the codons of [`BODY`] in turn
/// and `TAA`.
fn coding_sequence(codons: usize) -> String {
    let body: String = (0..codons - 2).map(|i| BODY[i % BODY.len()]).collect();
    format!("ATG{body}TAA")
}

/// The peptide [`coding_sequence`] translates to, without its stop.
fn peptide(codons: usize) -> String {
    let body: String = (0..codons - 2)
        .map(|i| char::from(RESIDUES[i % RESIDUES.len()]))
        .collect();
    format!("M{body}")
}

fn reverse_complement(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// A made-up sequence and its one gene, a single-exon coding transcript from
/// `start` to `end` on the forward or the reverse strand. A gene whose `start`
/// is above its `end` lies across the origin of a circular sequence, from
/// `start` to the last base and from base 1 to `end`.
struct Sequence {
    name: &'static str,
    length: u64,
    circular: bool,
    transcript: &'static str,
    gene: &'static str,
    start: u64,
    end: u64,
    forward: bool,
}

impl Sequence {
    fn across_origin(&self) -> bool {
        self.start > self.end
    }

    fn codons(&self) -> usize {
        let bases = if self.across_origin() {
            self.length - self.start + 1 + self.end
        } else {
            self.end - self.start + 1
        };
        (bases / 3) as usize
    }

    /// The sequence's bases: [`FILLER`] repeated, with the gene's coding
    /// sequence on its strand over the gene's span.
    fn bases(&self) -> String {
        let mut bases: Vec<u8> = FILLER
            .iter()
            .copied()
            .cycle()
            .take(self.length as usize)
            .collect();
        let coding = coding_sequence(self.codons());
        let genomic = if self.forward {
            coding
        } else {
            reverse_complement(&coding)
        };
        let from = (self.start - 1) as usize;
        for (offset, base) in genomic.bytes().enumerate() {
            let at = (from + offset) % bases.len();
            bases[at] = base;
        }
        String::from_utf8(bases).unwrap()
    }

    /// The transcript as the cache converter writes a dumped bacterial
    /// transcript: one exon, one mapper pair, the coding sequence from the
    /// first cDNA base to the last. A transcript across the origin is stored
    /// as Ensembl's dump stores it: the transcript and its exon start below 1,
    /// at `start` less the sequence's length (ensembl core
    /// `DBSQL/TranscriptAdaptor.pm:1941-1947`), and the mapper pair keeps
    /// `start` and `end` themselves (`TranscriptMapper.pm:160-161`).
    fn transcript_record(&self, db_id: u64) -> serde_json::Value {
        let codons = self.codons();
        let cdna_end = 3 * codons as u64;
        let ori = if self.forward { 1 } else { -1 };
        let protein = format!("{}P", self.transcript);
        let stored_start = if self.across_origin() {
            self.start as i64 - self.length as i64
        } else {
            self.start as i64
        };
        let exon = json!({
            "stable_id": format!("{}E", self.transcript),
            "start": stored_start,
            "end": self.end,
            "phase": 0,
            "end_phase": 0,
        });
        json!({
            "stable_id": self.transcript,
            "dbID": db_id,
            "gene_stable_id": self.gene,
            "gene_symbol": format!("gene{db_id}"),
            "gene_symbol_source": "made_up",
            "start": stored_start,
            "end": self.end,
            "strand": ori,
            "biotype": "protein_coding",
            "source": "ensembl",
            "is_canonical": 1,
            "protein": protein,
            "cdna_coding_start": 1,
            "cdna_coding_end": cdna_end,
            "exons": [exon.clone()],
            "translation": {"stable_id": protein, "dbID": db_id, "start": 1, "end": cdna_end},
            "attributes": [],
            "variation_effect_feature_cache": {
                "codon_table": 11,
                "introns": [],
                "mapper": {
                    "start_phase": 0,
                    "cdna_coding_start": 1,
                    "cdna_coding_end": cdna_end,
                    "pair_count": 1,
                    "pairs": [{
                        "from_start": 1,
                        "from_end": cdna_end,
                        "to_start": self.start,
                        "to_end": self.end,
                        "ori": ori,
                    }],
                },
                "peptide": peptide(codons),
                "protein_features": [],
                "seq_edits": [],
                "sorted_exons": [exon],
                "translateable_seq": coding_sequence(codons),
            },
        })
    }
}

/// The sequences of the cache:
///
/// - `fwd` and `rev`, 99 bases each, are one gene from the first base to the
///   last, on the forward and on the reverse strand;
/// - `near`, 4,000 bases, holds its gene at 2001-2099, within 5,000 bases of
///   both ends;
/// - `far`, 12,000 bases, holds its gene at 6001-6099, more than 5,000 bases
///   from both ends;
/// - `loop`, 3,000 bases and circular, holds its gene at 101-199: its last
///   base is 2,801 bases after the gene on the straight line and 101 bases
///   before it round the origin;
/// - `ring` and `ring2`, 12,000 bases and circular, each hold a gene of 150
///   bases across the origin, from 11911 to the last base and from base 1 to
///   60, on the forward and on the reverse strand.
const SEQUENCES: [Sequence; 7] = [
    Sequence {
        name: "fwd",
        length: 99,
        circular: false,
        transcript: "SEQT01",
        gene: "SEQG01",
        start: 1,
        end: 99,
        forward: true,
    },
    Sequence {
        name: "rev",
        length: 99,
        circular: false,
        transcript: "SEQT02",
        gene: "SEQG02",
        start: 1,
        end: 99,
        forward: false,
    },
    Sequence {
        name: "near",
        length: 4000,
        circular: false,
        transcript: "SEQT03",
        gene: "SEQG03",
        start: 2001,
        end: 2099,
        forward: true,
    },
    Sequence {
        name: "far",
        length: 12000,
        circular: false,
        transcript: "SEQT04",
        gene: "SEQG04",
        start: 6001,
        end: 6099,
        forward: true,
    },
    Sequence {
        name: "loop",
        length: 3000,
        circular: true,
        transcript: "SEQT07",
        gene: "SEQG07",
        start: 101,
        end: 199,
        forward: true,
    },
    Sequence {
        name: "ring",
        length: 12000,
        circular: true,
        transcript: "SEQT05",
        gene: "SEQG05",
        start: 11911,
        end: 60,
        forward: true,
    },
    Sequence {
        name: "ring2",
        length: 12000,
        circular: true,
        transcript: "SEQT06",
        gene: "SEQG06",
        start: 11911,
        end: 60,
        forward: false,
    },
];

/// Whether the cache carries the optional side file that states each
/// sequence's length and whether it is circular.
#[derive(Clone, Copy, PartialEq)]
enum SideFile {
    Absent,
    Present,
}

/// Writes the cache of [`SEQUENCES`] under `dir/cache` with its checksum
/// manifest, so a run reads it without a warning.
fn write_cache(dir: &Path, side_file: SideFile) -> PathBuf {
    let cache = dir.join("cache");
    let names: Vec<&str> = SEQUENCES.iter().map(|sequence| sequence.name).collect();
    for (index, sequence) in SEQUENCES.iter().enumerate() {
        let shard_dir = cache.join("transcripts").join(sequence.name);
        fs::create_dir_all(&shard_dir).unwrap();
        let records = json!([sequence.transcript_record(index as u64 + 1)]);
        fs::write(shard_dir.join("1-1000000.json"), records.to_string()).unwrap();
    }
    let info = json!({
        "species": "synthetic_species",
        "assembly": "SYN1",
        "cache_version": 116,
        "contig_order": names,
        "source_type": "ensembl",
        "source_versions": {"assembly": "SYN1", "genebuild": "synthetic_1"},
    });
    fs::write(cache.join("info.json"), info.to_string()).unwrap();
    if side_file == SideFile::Present {
        let mut table = String::from("name\tlength\tcircular\n");
        for sequence in &SEQUENCES {
            let circular = u8::from(sequence.circular);
            table.push_str(&format!(
                "{}\t{}\t{circular}\n",
                sequence.name, sequence.length
            ));
        }
        fs::write(cache.join("sequence_regions.tsv"), table).unwrap();
    }
    CacheManifest::write(&cache).unwrap();
    cache
}

/// Writes [`SEQUENCES`] as an uncompressed FASTA with its index, one line a
/// sequence, and returns its path.
fn write_fasta(dir: &Path) -> PathBuf {
    let fasta = dir.join("genome.fa");
    let mut text = String::new();
    let mut index = String::new();
    for sequence in &SEQUENCES {
        text.push_str(&format!(">{}\n", sequence.name));
        let offset = text.len();
        let length = sequence.length;
        index.push_str(&format!(
            "{}\t{length}\t{offset}\t{length}\t{}\n",
            sequence.name,
            length + 1
        ));
        text.push_str(&sequence.bases());
        text.push('\n');
    }
    fs::write(&fasta, text).unwrap();
    fs::write(dir.join("genome.fa.fai"), index).unwrap();
    fasta
}

/// Runs `vep` on the records of `cases`, read as `format`, against `cache` and
/// returns the output's data lines; the run writes no warnings file.
fn annotate(
    dir: &Path,
    cache: &Path,
    name: &str,
    format: &str,
    cases: &[Case],
    fasta: Option<&Path>,
) -> Vec<String> {
    let mut text = String::new();
    if format == "vcf" {
        text.push_str(VCF_HEADER);
    }
    for (record, _) in cases {
        text.push_str(record);
        text.push('\n');
    }
    let input = dir.join(format!("{name}.in"));
    fs::write(&input, text).unwrap();
    let out = dir.join(format!("{name}.txt"));
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(cache)
        .args(["--species", "synthetic_species", "--assembly", "SYN1"])
        .args(["--format", format])
        .args(["--force_overwrite", "--no_stats", "--quiet"]);
    if let Some(fasta) = fasta {
        command.arg("--fasta").arg(fasta);
    }
    let output = command.stdin(Stdio::null()).output().unwrap();
    assert!(
        output.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let warnings = dir.join(format!("{name}.txt_warnings.txt"));
    assert!(
        !warnings.exists(),
        "{name}: {}",
        fs::read_to_string(&warnings).unwrap_or_default()
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// Asserts that the records of `cases`, read as `format`, get the rows the
/// table gives them, without a FASTA, where any coordinate is taken, and with
/// one, where the sequence's length is known.
fn assert_rows(format: &str, cases: &[Case], side_file: SideFile) {
    let dir = tempfile::tempdir().unwrap();
    let cache = write_cache(dir.path(), side_file);
    let fasta = write_fasta(dir.path());
    let expected: Vec<&str> = cases
        .iter()
        .flat_map(|(_, rows)| rows.iter().copied())
        .collect();
    for (name, fasta) in [
        ("without-fasta", None),
        ("with-fasta", Some(fasta.as_path())),
    ] {
        let rows = annotate(dir.path(), &cache, name, format, cases, fasta);
        assert_eq!(rows, expected, "{format} {name}");
    }
}

/// A substitution at base 1: the start codon's first base on the forward gene
/// that begins there, the stop codon's last base on the reverse one. On `near`
/// the gene is 2,000 bases downstream of the record, and on `far` no
/// transcript is within 5,000 bases.
const FIRST_BASE_SUBSTITUTIONS: &[Case] = &[
    (
        "fwd\t1\tfwd_first_snv\tA\tC\t.\t.\t.",
        &[
            "fwd_first_snv\tfwd:1\tC\tSEQG01\tSEQT01\tTranscript\tstart_lost\t1\t1\t1\tM/L\tAtg/Ctg\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "rev\t1\trev_first_snv\tT\tG\t.\t.\t.",
        &[
            "rev_first_snv\trev:1\tG\tSEQG02\tSEQT02\tTranscript\tstop_lost\t99\t99\t33\t*/Y\ttaA/taC\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "near\t1\tnear_first_snv\tA\tC\t.\t.\t.",
        &[
            "near_first_snv\tnear:1\tC\tSEQG03\tSEQT03\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2000;STRAND=1",
        ],
    ),
    (
        "far\t1\tfar_first_snv\tA\tC\t.\t.\t.",
        &[
            "far_first_snv\tfar:1\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// An insertion before base 1 and a deletion of base 1, in Ensembl's default
/// input format. The insertion has start 1 and end 0, which `validate_vf`
/// admits (`Parser.pm:483`: a start is refused only past the end plus one);
/// its Location is `0-1` and it lies outside a gene that begins at base 1, at
/// distance 0.
const FIRST_BASE_INDELS: &[Case] = &[
    (
        "fwd\t1\t0\t-/G\t+\tfwd_before_first_ins",
        &[
            "fwd_before_first_ins\tfwd:0-1\tG\tSEQG01\tSEQT01\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=0;STRAND=1",
        ],
    ),
    (
        "fwd\t1\t1\tA/-\t+\tfwd_first_del",
        &[
            "fwd_first_del\tfwd:1\t-\tSEQG01\tSEQT01\tTranscript\tframeshift_variant,start_lost\t1\t1\t1\tM/X\tAtg/tg\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "rev\t1\t0\t-/G\t+\trev_before_first_ins",
        &[
            "rev_before_first_ins\trev:0-1\tG\tSEQG02\tSEQT02\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=0;STRAND=-1",
        ],
    ),
    (
        "rev\t1\t1\tT/-\t+\trev_first_del",
        &[
            "rev_first_del\trev:1\t-\tSEQG02\tSEQT02\tTranscript\tstop_lost\t99\t99\t33\t*/X\ttaA/ta\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "near\t1\t0\t-/G\t+\tnear_before_first_ins",
        &[
            "near_before_first_ins\tnear:0-1\tG\tSEQG03\tSEQT03\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2000;STRAND=1",
        ],
    ),
    (
        "near\t1\t1\tA/-\t+\tnear_first_del",
        &[
            "near_first_del\tnear:1\t-\tSEQG03\tSEQT03\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2000;STRAND=1",
        ],
    ),
    (
        "far\t1\t0\t-/G\t+\tfar_before_first_ins",
        &[
            "far_before_first_ins\tfar:0-1\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
    (
        "far\t1\t1\tA/-\t+\tfar_first_del",
        &[
            "far_first_del\tfar:1\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// A VCF record at `POS` 1 whose two-base reference allele keeps its second
/// base: the alleles share no first base, so the record is the deletion of
/// base 1.
const FIRST_BASE_DELETIONS_IN_VCF: &[Case] = &[
    (
        "fwd\t1\tfwd_first_del_vcf\tAT\tT\t.\t.\t.",
        &[
            "fwd_first_del_vcf\tfwd:1\t-\tSEQG01\tSEQT01\tTranscript\tframeshift_variant,start_lost\t1\t1\t1\tM/X\tAtg/tg\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "near\t1\tnear_first_del_vcf\tAC\tC\t.\t.\t.",
        &[
            "near_first_del_vcf\tnear:1\t-\tSEQG03\tSEQT03\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2000;STRAND=1",
        ],
    ),
    (
        "far\t1\tfar_first_del_vcf\tAC\tC\t.\t.\t.",
        &[
            "far_first_del_vcf\tfar:1\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// A substitution at the last base: the stop codon's last base on the forward
/// gene that ends there, the start codon's first base on the reverse one.
const LAST_BASE_SUBSTITUTIONS: &[Case] = &[
    (
        "fwd\t99\tfwd_last_snv\tA\tC\t.\t.\t.",
        &[
            "fwd_last_snv\tfwd:99\tC\tSEQG01\tSEQT01\tTranscript\tstop_lost\t99\t99\t33\t*/Y\ttaA/taC\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "rev\t99\trev_last_snv\tT\tG\t.\t.\t.",
        &[
            "rev_last_snv\trev:99\tG\tSEQG02\tSEQT02\tTranscript\tstart_lost\t1\t1\t1\tM/L\tAtg/Ctg\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "near\t4000\tnear_last_snv\tA\tC\t.\t.\t.",
        &[
            "near_last_snv\tnear:4000\tC\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1901;STRAND=1",
        ],
    ),
    (
        "far\t12000\tfar_last_snv\tA\tC\t.\t.\t.",
        &[
            "far_last_snv\tfar:12000\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// Deletions of the last two bases and of the last base.
const LAST_BASE_DELETIONS: &[Case] = &[
    (
        "fwd\t97\tfwd_last_del2\tTAA\tT\t.\t.\t.",
        &[
            "fwd_last_del2\tfwd:98-99\t-\tSEQG01\tSEQT01\tTranscript\tstop_lost\t98-99\t98-99\t33\t*/X\ttAA/t\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "fwd\t98\tfwd_last_del1\tAA\tA\t.\t.\t.",
        &[
            "fwd_last_del1\tfwd:99\t-\tSEQG01\tSEQT01\tTranscript\tstop_lost\t99\t99\t33\t*/X\ttaA/ta\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "rev\t97\trev_last_del2\tCAT\tC\t.\t.\t.",
        &[
            "rev_last_del2\trev:98-99\t-\tSEQG02\tSEQT02\tTranscript\tframeshift_variant,start_lost\t1-2\t1-2\t1\tM/X\tATg/g\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "rev\t98\trev_last_del1\tAT\tA\t.\t.\t.",
        &[
            "rev_last_del1\trev:99\t-\tSEQG02\tSEQT02\tTranscript\tframeshift_variant,start_lost\t1\t1\t1\tM/X\tAtg/tg\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "near\t3998\tnear_last_del2\tGCA\tG\t.\t.\t.",
        &[
            "near_last_del2\tnear:3999-4000\t-\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1900;STRAND=1",
        ],
    ),
    (
        "near\t3999\tnear_last_del1\tCA\tC\t.\t.\t.",
        &[
            "near_last_del1\tnear:4000\t-\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1901;STRAND=1",
        ],
    ),
    (
        "far\t11998\tfar_last_del2\tGCA\tG\t.\t.\t.",
        &[
            "far_last_del2\tfar:11999-12000\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
    (
        "far\t11999\tfar_last_del1\tCA\tC\t.\t.\t.",
        &[
            "far_last_del1\tfar:12000\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// An insertion after the last base, a VCF record at the last base whose
/// alternate allele adds a base: its Location names the last base and the
/// position after it.
const AFTER_LAST_BASE_INSERTIONS: &[Case] = &[
    (
        "fwd\t99\tfwd_after_last_ins\tA\tAG\t.\t.\t.",
        &[
            "fwd_after_last_ins\tfwd:99-100\tG\tSEQG01\tSEQT01\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=0;STRAND=1",
        ],
    ),
    (
        "rev\t99\trev_after_last_ins\tT\tTG\t.\t.\t.",
        &[
            "rev_after_last_ins\trev:99-100\tG\tSEQG02\tSEQT02\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=0;STRAND=-1",
        ],
    ),
    (
        "near\t4000\tnear_after_last_ins\tA\tAG\t.\t.\t.",
        &[
            "near_after_last_ins\tnear:4000-4001\tG\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1901;STRAND=1",
        ],
    ),
    (
        "far\t12000\tfar_after_last_ins\tA\tAG\t.\t.\t.",
        &[
            "far_after_last_ins\tfar:12000-12001\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// VCF records past the last base: a deletion at the base before the last
/// whose four-base reference allele runs two bases past the end, and a
/// substitution ten bases past the end. Each is annotated where it is written.
const PAST_LAST_BASE: &[Case] = &[
    (
        "rev\t98\trev_past_last_del\tATAC\tA\t.\t.\t.",
        &[
            "rev_past_last_del\trev:99-101\t-\tSEQG02\tSEQT02\tTranscript\tcoding_sequence_variant,5_prime_UTR_variant\t?-1\t?-1\t?-1\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1",
        ],
    ),
    (
        "near\t3999\tnear_past_last_del\tCAAC\tC\t.\t.\t.",
        &[
            "near_past_last_del\tnear:4000-4002\t-\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1901;STRAND=1",
        ],
    ),
    (
        "far\t11999\tfar_past_last_del\tCAAC\tC\t.\t.\t.",
        &[
            "far_past_last_del\tfar:12000-12002\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
    (
        "fwd\t109\tfwd_past_last_snv\tA\tG\t.\t.\t.",
        &[
            "fwd_past_last_snv\tfwd:109\tG\tSEQG01\tSEQT01\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=10;STRAND=1",
        ],
    ),
    (
        "rev\t109\trev_past_last_snv\tA\tG\t.\t.\t.",
        &[
            "rev_past_last_snv\trev:109\tG\tSEQG02\tSEQT02\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=10;STRAND=-1",
        ],
    ),
    (
        "near\t4010\tnear_past_last_snv\tA\tG\t.\t.\t.",
        &[
            "near_past_last_snv\tnear:4010\tG\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1911;STRAND=1",
        ],
    ),
    (
        "far\t12010\tfar_past_last_snv\tA\tG\t.\t.\t.",
        &[
            "far_past_last_snv\tfar:12010\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// The deletions of [`PAST_LAST_BASE`] in Ensembl's default input format.
const PAST_LAST_BASE_IN_DEFAULT_FORMAT: &[Case] = &[
    (
        "rev\t99\t101\tTAC/-\t+\trev_past_last_del_ens",
        &[
            "rev_past_last_del_ens\trev:99-101\t-\tSEQG02\tSEQT02\tTranscript\tcoding_sequence_variant,5_prime_UTR_variant\t?-1\t?-1\t?-1\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1",
        ],
    ),
    (
        "near\t4000\t4002\tAAC/-\t+\tnear_past_last_del_ens",
        &[
            "near_past_last_del_ens\tnear:4000-4002\t-\tSEQG03\tSEQT03\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=1901;STRAND=1",
        ],
    ),
    (
        "far\t12000\t12002\tAAC/-\t+\tfar_past_last_del_ens",
        &[
            "far_past_last_del_ens\tfar:12000-12002\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
];

/// The deletion past the last base on `fwd`, where it takes the last base of
/// the stop codon and has no cDNA end. Ensembl VEP writes
/// `stop_retained_variant` in place of `coding_sequence_variant` on these two
/// rows, the documented difference `stop_retained_without_stop_codon_test`
/// (`docs/divergences/stop-retained-without-stop-codon-test.md`).
const PAST_LAST_BASE_OVER_A_STOP_CODON: &[Case] = &[
    (
        "fwd\t98\tfwd_past_last_del\tAAAC\tA\t.\t.\t.",
        &[
            "fwd_past_last_del\tfwd:99-101\t-\tSEQG01\tSEQT01\tTranscript\tcoding_sequence_variant,3_prime_UTR_variant\t99-?\t99-?\t33-?\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1",
        ],
    ),
];

/// The deletion of [`PAST_LAST_BASE_OVER_A_STOP_CODON`] in Ensembl's default
/// input format.
const PAST_LAST_BASE_OVER_A_STOP_CODON_IN_DEFAULT_FORMAT: &[Case] = &[
    (
        "fwd\t99\t101\tAAC/-\t+\tfwd_past_last_del_ens",
        &[
            "fwd_past_last_del_ens\tfwd:99-101\t-\tSEQG01\tSEQT01\tTranscript\tcoding_sequence_variant,3_prime_UTR_variant\t99-?\t99-?\t33-?\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1",
        ],
    ),
];

/// VCF records at the two ends of the circular sequence `loop`. The record at
/// its last base is `downstream_gene_variant` at the straight-line distance
/// 2,801, and the deletion whose reference allele is the last two bases
/// followed by the first two is a deletion of the last base and the two
/// positions past it.
const ORIGIN: &[Case] = &[
    (
        "loop\t1\tloop_first_snv\tA\tC\t.\t.\t.",
        &[
            "loop_first_snv\tloop:1\tC\tSEQG07\tSEQT07\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=100;STRAND=1",
        ],
    ),
    (
        "loop\t3000\tloop_last_snv\tA\tC\t.\t.\t.",
        &[
            "loop_last_snv\tloop:3000\tC\tSEQG07\tSEQT07\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2801;STRAND=1",
        ],
    ),
    (
        "loop\t2999\tloop_last_del1\tCA\tC\t.\t.\t.",
        &[
            "loop_last_del1\tloop:3000\t-\tSEQG07\tSEQT07\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2801;STRAND=1",
        ],
    ),
    (
        "loop\t2999\tloop_past_last_del\tCAAC\tC\t.\t.\t.",
        &[
            "loop_past_last_del\tloop:3000-3002\t-\tSEQG07\tSEQT07\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2801;STRAND=1",
        ],
    ),
    (
        "loop\t3000\tloop_after_last_ins\tA\tAG\t.\t.\t.",
        &[
            "loop_after_last_ins\tloop:3000-3001\tG\tSEQG07\tSEQT07\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2801;STRAND=1",
        ],
    ),
    (
        "loop\t3010\tloop_past_last_snv\tA\tG\t.\t.\t.",
        &[
            "loop_past_last_snv\tloop:3010\tG\tSEQG07\tSEQT07\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2811;STRAND=1",
        ],
    ),
];

/// Default-format records at the two ends of `loop`.
const ORIGIN_IN_DEFAULT_FORMAT: &[Case] = &[
    (
        "loop\t1\t0\t-/G\t+\tloop_before_first_ins",
        &[
            "loop_before_first_ins\tloop:0-1\tG\tSEQG07\tSEQT07\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=100;STRAND=1",
        ],
    ),
    (
        "loop\t1\t1\tA/-\t+\tloop_first_del",
        &[
            "loop_first_del\tloop:1\t-\tSEQG07\tSEQT07\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=100;STRAND=1",
        ],
    ),
    (
        "loop\t3000\t3002\tAAC/-\t+\tloop_past_last_del_ens",
        &[
            "loop_past_last_del_ens\tloop:3000-3002\t-\tSEQG07\tSEQT07\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2801;STRAND=1",
        ],
    ),
];

/// VCF records on `ring` and `ring2`, whose gene is stored across the origin:
/// inside the gene among the first bases (1 and 30) and among the last (11950
/// and the last base), deletions and an insertion at the last base, the
/// deletion whose reference allele runs from the last bases into the first,
/// and a record 5,940 bases from the gene. Each row on the gene is the row
/// Ensembl VEP 116.3 writes for the same bases and the same edit where the
/// gene lies unbroken on a straight sequence, at 11911-12060, with the record's
/// own Location: the coding positions count from the gene's first base on
/// both sides of the origin. On the circular sequence itself Ensembl VEP
/// writes `intergenic_variant` on the transcript for a record among the first
/// bases, where its mapper pair, which runs from the database start down to
/// the database end, maps no base (ensembl core `Mapper.pm:316`,
/// ensembl-variation `BaseVariationFeatureOverlapAllele.pm:285`), and no row
/// of the transcript for a record among the last bases, which the stored span
/// does not reach (`AnnotationSource.pm:306-309`).
const ACROSS_ORIGIN: &[Case] = &[
    (
        "ring\t1\tring_first_snv\tC\tA\t.\t.\t.",
        &[
            "ring_first_snv\tring:1\tA\tSEQG05\tSEQT05\tTranscript\tmissense_variant\t91\t91\t31\tQ/K\tCag/Aag\t-\tIMPACT=MODERATE;STRAND=1",
        ],
    ),
    (
        "ring\t1\tring_first_del_vcf\tCA\tA\t.\t.\t.",
        &[
            "ring_first_del_vcf\tring:1\t-\tSEQG05\tSEQT05\tTranscript\tframeshift_variant\t91\t91\t31\tQ/X\tCag/ag\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "ring\t30\tring_low_snv\tT\tG\t.\t.\t.",
        &[
            "ring_low_snv\tring:30\tG\tSEQG05\tSEQT05\tTranscript\tsynonymous_variant\t120\t120\t40\tA\tgcT/gcG\t-\tIMPACT=LOW;STRAND=1",
        ],
    ),
    (
        "ring\t6000\tring_mid_snv\tA\tC\t.\t.\t.",
        &[
            "ring_mid_snv\tring:6000\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
    (
        "ring\t11950\tring_high_snv\tC\tA\t.\t.\t.",
        &[
            "ring_high_snv\tring:11950\tA\tSEQG05\tSEQT05\tTranscript\tmissense_variant\t40\t40\t14\tP/T\tCcg/Acg\t-\tIMPACT=MODERATE;STRAND=1",
        ],
    ),
    (
        "ring\t11998\tring_last_del2\tGAT\tG\t.\t.\t.",
        &[
            "ring_last_del2\tring:11999-12000\t-\tSEQG05\tSEQT05\tTranscript\tframeshift_variant\t89-90\t89-90\t30\tD/X\tgAT/g\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "ring\t11999\tring_last_del1\tAT\tA\t.\t.\t.",
        &[
            "ring_last_del1\tring:12000\t-\tSEQG05\tSEQT05\tTranscript\tframeshift_variant\t90\t90\t30\tD/X\tgaT/ga\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "ring\t11999\tring_past_last_del\tATCA\tA\t.\t.\t.",
        &[
            "ring_past_last_del\tring:12000-12002\t-\tSEQG05\tSEQT05\tTranscript\tinframe_deletion\t90-92\t90-92\t30-31\tDQ/E\tgaTCAg/gag\t-\tIMPACT=MODERATE;STRAND=1",
        ],
    ),
    (
        "ring\t12000\tring_last_snv\tT\tG\t.\t.\t.",
        &[
            "ring_last_snv\tring:12000\tG\tSEQG05\tSEQT05\tTranscript\tmissense_variant\t90\t90\t30\tD/E\tgaT/gaG\t-\tIMPACT=MODERATE;STRAND=1",
        ],
    ),
    (
        "ring\t12000\tring_after_last_ins\tT\tTG\t.\t.\t.",
        &[
            "ring_after_last_ins\tring:12000-12001\tG\tSEQG05\tSEQT05\tTranscript\tframeshift_variant\t90-91\t90-91\t30-31\t-/X\t-/G\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "ring2\t1\tring2_first_snv\tA\tC\t.\t.\t.",
        &[
            "ring2_first_snv\tring2:1\tC\tSEQG06\tSEQT06\tTranscript\tmissense_variant\t60\t60\t20\tC/W\ttgT/tgG\t-\tIMPACT=MODERATE;STRAND=-1",
        ],
    ),
    (
        "ring2\t1\tring2_first_del_vcf\tAC\tC\t.\t.\t.",
        &[
            "ring2_first_del_vcf\tring2:1\t-\tSEQG06\tSEQT06\tTranscript\tframeshift_variant\t60\t60\t20\tC/X\ttgT/tg\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "ring2\t30\tring2_low_snv\tC\tA\t.\t.\t.",
        &[
            "ring2_low_snv\tring2:30\tA\tSEQG06\tSEQT06\tTranscript\tmissense_variant\t31\t31\t11\tD/Y\tGat/Tat\t-\tIMPACT=MODERATE;STRAND=-1",
        ],
    ),
    (
        "ring2\t6000\tring2_mid_snv\tA\tC\t.\t.\t.",
        &[
            "ring2_mid_snv\tring2:6000\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
        ],
    ),
    (
        "ring2\t11950\tring2_high_snv\tG\tT\t.\t.\t.",
        &[
            "ring2_high_snv\tring2:11950\tT\tSEQG06\tSEQT06\tTranscript\tstop_gained\t111\t111\t37\tY/*\ttaC/taA\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "ring2\t11998\tring2_last_del2\tAGC\tA\t.\t.\t.",
        &[
            "ring2_last_del2\tring2:11999-12000\t-\tSEQG06\tSEQT06\tTranscript\tframeshift_variant\t61-62\t61-62\t21\tA/X\tGCt/t\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "ring2\t11999\tring2_last_del1\tGC\tG\t.\t.\t.",
        &[
            "ring2_last_del1\tring2:12000\t-\tSEQG06\tSEQT06\tTranscript\tframeshift_variant\t61\t61\t21\tA/X\tGct/ct\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "ring2\t11999\tring2_past_last_del\tGCAC\tG\t.\t.\t.",
        &[
            "ring2_past_last_del\tring2:12000-12002\t-\tSEQG06\tSEQT06\tTranscript\tinframe_deletion\t59-61\t59-61\t20-21\tCA/S\ttGTGct/tct\t-\tIMPACT=MODERATE;STRAND=-1",
        ],
    ),
    (
        "ring2\t12000\tring2_last_snv\tC\tA\t.\t.\t.",
        &[
            "ring2_last_snv\tring2:12000\tA\tSEQG06\tSEQT06\tTranscript\tmissense_variant\t61\t61\t21\tA/S\tGct/Tct\t-\tIMPACT=MODERATE;STRAND=-1",
        ],
    ),
    (
        "ring2\t12000\tring2_after_last_ins\tC\tCG\t.\t.\t.",
        &[
            "ring2_after_last_ins\tring2:12000-12001\tG\tSEQG06\tSEQT06\tTranscript\tframeshift_variant\t60-61\t60-61\t20-21\t-/X\t-/C\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
];

/// Default-format records on `ring` and `ring2`: an insertion before base 1,
/// which is the insertion after the last base, a deletion of base 1, and the
/// deletion of the last base and the first two.
const ACROSS_ORIGIN_IN_DEFAULT_FORMAT: &[Case] = &[
    (
        "ring\t1\t0\t-/G\t+\tring_before_first_ins",
        &[
            "ring_before_first_ins\tring:0-1\tG\tSEQG05\tSEQT05\tTranscript\tframeshift_variant\t90-91\t90-91\t30-31\t-/X\t-/G\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "ring\t1\t1\tC/-\t+\tring_first_del",
        &[
            "ring_first_del\tring:1\t-\tSEQG05\tSEQT05\tTranscript\tframeshift_variant\t91\t91\t31\tQ/X\tCag/ag\t-\tIMPACT=HIGH;STRAND=1",
        ],
    ),
    (
        "ring\t12000\t12002\tTCA/-\t+\tring_past_last_del_ens",
        &[
            "ring_past_last_del_ens\tring:12000-12002\t-\tSEQG05\tSEQT05\tTranscript\tinframe_deletion\t90-92\t90-92\t30-31\tDQ/E\tgaTCAg/gag\t-\tIMPACT=MODERATE;STRAND=1",
        ],
    ),
    (
        "ring2\t1\t0\t-/G\t+\tring2_before_first_ins",
        &[
            "ring2_before_first_ins\tring2:0-1\tG\tSEQG06\tSEQT06\tTranscript\tframeshift_variant\t60-61\t60-61\t20-21\t-/X\t-/C\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "ring2\t1\t1\tA/-\t+\tring2_first_del",
        &[
            "ring2_first_del\tring2:1\t-\tSEQG06\tSEQT06\tTranscript\tframeshift_variant\t60\t60\t20\tC/X\ttgT/tg\t-\tIMPACT=HIGH;STRAND=-1",
        ],
    ),
    (
        "ring2\t12000\t12002\tCAC/-\t+\tring2_past_last_del_ens",
        &[
            "ring2_past_last_del_ens\tring2:12000-12002\t-\tSEQG06\tSEQT06\tTranscript\tinframe_deletion\t59-61\t59-61\t20-21\tCA/S\ttGTGct/tct\t-\tIMPACT=MODERATE;STRAND=-1",
        ],
    ),
];

#[test]
fn a_substitution_at_the_first_base_is_annotated() {
    assert_rows("vcf", FIRST_BASE_SUBSTITUTIONS, SideFile::Absent);
}

#[test]
fn an_insertion_before_the_first_base_and_a_deletion_of_it_are_annotated() {
    assert_rows("ensembl", FIRST_BASE_INDELS, SideFile::Absent);
    assert_rows("vcf", FIRST_BASE_DELETIONS_IN_VCF, SideFile::Absent);
}

#[test]
fn a_substitution_at_the_last_base_is_annotated() {
    assert_rows("vcf", LAST_BASE_SUBSTITUTIONS, SideFile::Absent);
}

#[test]
fn a_deletion_ending_at_the_last_base_is_annotated() {
    assert_rows("vcf", LAST_BASE_DELETIONS, SideFile::Absent);
}

#[test]
fn an_insertion_after_the_last_base_is_annotated() {
    assert_rows("vcf", AFTER_LAST_BASE_INSERTIONS, SideFile::Absent);
}

#[test]
fn a_reference_allele_running_past_the_last_base_is_treated_as_ensembl_vep_treats_it() {
    assert_rows("vcf", PAST_LAST_BASE, SideFile::Absent);
    assert_rows(
        "ensembl",
        PAST_LAST_BASE_IN_DEFAULT_FORMAT,
        SideFile::Absent,
    );
    assert_rows("vcf", PAST_LAST_BASE_OVER_A_STOP_CODON, SideFile::Absent);
    assert_rows(
        "ensembl",
        PAST_LAST_BASE_OVER_A_STOP_CODON_IN_DEFAULT_FORMAT,
        SideFile::Absent,
    );
}

/// The rows are those of a straight line whether or not the cache's side file
/// declares the sequence circular.
#[test]
fn a_record_at_the_origin_of_a_circular_sequence_matches_ensembl_vep() {
    for side_file in [SideFile::Present, SideFile::Absent] {
        assert_rows("vcf", ORIGIN, side_file);
        assert_rows("ensembl", ORIGIN_IN_DEFAULT_FORMAT, side_file);
    }
}

/// The loader reads the transcript into coordinates that run on past the
/// last base, and a record on either side of the origin gets the transcript's
/// rows, with or without the cache's side file.
#[test]
fn a_transcript_stored_across_the_origin_is_loaded_and_annotated() {
    for side_file in [SideFile::Absent, SideFile::Present] {
        assert_rows("vcf", ACROSS_ORIGIN, side_file);
        assert_rows("ensembl", ACROSS_ORIGIN_IN_DEFAULT_FORMAT, side_file);
    }
}
