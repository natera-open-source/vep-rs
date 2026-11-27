// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Downstream, NMD and BLOSUM62 classes and the classes of a plugin field that follows
//! a pair of Consequence sets, on small outputs and a small cache. The expected values are
//! those `tools/vep-formats-compare` gives on the same inputs.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use vep_plugin_pairs::pytext::Json;
use vep_plugin_pairs::rules::downstream::classes::{
    downstream_class, hgvs_rotated_allele, later_allele_shared_shift, nmd_reads_shifted_pair,
};
use vep_plugin_pairs::rules::downstream::codes::{blosum62, genetic_code};
use vep_plugin_pairs::rules::downstream::readings::shifted_pair_readings;
use vep_plugin_pairs::rules::downstream::transcript::{codon_table_id, Sequence};
use vep_plugin_pairs::rules::downstream::DownstreamSequences;
use vep_plugin_pairs::rules::entry::{parse_default, Entry, Row, ABSENT};
use vep_plugin_pairs::rules::term_pair::{
    registry_class, term_pair_following_cell, terms_impact, START_COOCCURRENCE_SWAP, TERMS_BY_RANK,
};
use vep_plugin_pairs::rules::{BuiltIn, Consult, RecordRows, Rules};

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

const UTR_FRAME: &str = "downstream_utr_frame_past_coding_end";
const START: &str = "downstream_start_past_own_codon";
const SHIFTED: &str = "plugin_hgvs_shifted_cds_coordinates";
const ROTATION: &str = "downstream_reverse_strand_rotation";
const BLOSUM62_SHARED_SHIFT: &str = "blosum62_hgvs_shift_shared_across_alleles";
const SWAP: &str = "splice_family_swap";
const FOLLOWS_STOP: &str = "plugin_value_follows_stop_codon_term_pair";
const FOLLOWS_SPLICE: &str = "plugin_value_follows_splice_term_past_the_record_span";

// A coding sequence of eight codons and a stop codon whose 3' UTR, read from its first
// base, reaches a stop after three codons; and one of four codons and a stop codon whose
// second codon, less its third base, spells a stop with the base after it.
const TOY_CODING: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";
const TOY_UTR: &str = "GGGCCCAAATAGCCC";
const SHORT_CODING: &str = "ATGTACAAAGGTTAA";
const SHORT_UTR: &str = "GTGACC";
// A 3' UTR the vertebrate mitochondrial code reads unlike the standard one in both frames
// the classes read.
const TABLE_2_UTR: &str = "GTGACCAAGATAGCC";
const RUN_CODING: &str = "ATGCCAAAAGGTTAA";

const TOY: &str = "ENST00000000101";
const TOY_NO_UTR: &str = "ENST00000000102";
const TOY_TABLE_2: &str = "ENST00000000103";
const TOY_UTR_STOP: &str = "ENST00000000104";
const SHORT: &str = "ENST00000000105";
const TOY_TABLE_7: &str = "ENST00000000106";
const TOY_FORWARD: &str = "ENST00000000201";
const TOY_REVERSE: &str = "ENST00000000202";
const TOY_TWO_EXONS: &str = "ENST00000000203";
const TOY_SEAM: &str = "ENST00000000204";
const RUN: &str = "ENST00000000205";
const THREE_EXONS: &str = "ENST00000000206";

const NMD_VALUE: &str = "NMD_escaping_variant";
const ACCEPTOR_TERMS: &str = "splice_acceptor_variant,coding_sequence_variant,intron_variant";

/// An empty directory of one test's own, removed with what it holds when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-downstream-{}-{name}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).expect("a scratch directory");
        Scratch(directory)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A cache transcript without a model: its sequences as the members `sequences` spell them.
fn transcript(stable_id: &str, sequences: &str) -> String {
    format!(r#"{{"stable_id": "{stable_id}", "variation_effect_feature_cache": {{{sequences}}}}}"#)
}

/// A cache transcript with a model: each of `pairs` a run of cDNA coordinates and the
/// genomic bases it maps to, the exons those genomic runs unless `exons` gives them, the
/// coding sequence from cDNA base 1; `text` writes every number as a string.
#[allow(clippy::too_many_arguments)]
fn modelled(
    stable_id: &str,
    strand: i64,
    pairs: &[(i64, i64, i64, i64)],
    coding: &str,
    utr: &str,
    exons: &[(i64, i64)],
    text: bool,
    more: &str,
) -> String {
    let number = |value: i64| {
        if text {
            format!("\"{value}\"")
        } else {
            value.to_string()
        }
    };
    let mapper: Vec<String> = pairs
        .iter()
        .map(|(a, b, c, d)| {
            format!(
                r#"{{"from_start": {}, "from_end": {}, "to_start": {}, "to_end": {}, "ori": {}}}"#,
                number(*a),
                number(*b),
                number(*c),
                number(*d),
                number(strand)
            )
        })
        .collect();
    let mut spans: Vec<(i64, i64)> = if exons.is_empty() {
        pairs.iter().map(|(_, _, c, d)| (*c, *d)).collect()
    } else {
        exons.to_vec()
    };
    spans.sort_unstable();
    if strand == -1 {
        spans.reverse();
    }
    let spans: Vec<String> = spans
        .iter()
        .map(|(start, end)| {
            format!(
                r#"{{"start": {}, "end": {}}}"#,
                number(*start),
                number(*end)
            )
        })
        .collect();
    format!(
        r#"{{"stable_id": "{stable_id}", "strand": {}, "cdna_coding_start": {}, "cdna_coding_end": {}, "exons": [{}], "variation_effect_feature_cache": {{"translateable_seq": "{coding}", "three_prime_utr": "{utr}", "codon_table": 1, "mapper": {{"pairs": [{}]}}{more}}}}}"#,
        number(strand),
        number(1),
        number(coding.len() as i64),
        spans.join(", "),
        mapper.join(", "),
    )
}

fn write_shard(root: &Path, chromosome: &str, transcripts: &[String]) {
    let directory = root.join("transcripts").join(chromosome);
    fs::create_dir_all(&directory).expect("a chromosome directory");
    fs::write(
        directory.join("1-1000000.json"),
        format!("[{}]", transcripts.join(", ")),
    )
    .expect("a shard");
}

/// A candidate JSON cache of twelve coding transcripts. Chromosome 21 holds the toy
/// transcript (its `codon_table` the string a converted cache writes), a copy without a 3'
/// UTR (the number a native cache writes), a copy under codon table 2, a copy whose 3' UTR
/// opens with a stop codon (no `codon_table`) and a copy under table 7, which the codes do
/// not hold, none with a model, and the six transcripts with one (the reverse-strand one's
/// numbers the strings a converted cache writes, the seam's with its `peptide`); chromosome
/// 22 holds the short transcript (table 0, read as table 1).
fn downstream_cache(name: &str) -> Scratch {
    let root = Scratch::new(name);
    let three_exons_coding = format!("ATG{}TAA", "GCT".repeat(38));
    write_shard(
        &root.0,
        "21",
        &[
            transcript(
                TOY,
                &format!(
                    r#""translateable_seq": "{TOY_CODING}", "three_prime_utr": "{TOY_UTR}", "codon_table": "1""#
                ),
            ),
            transcript(
                TOY_NO_UTR,
                &format!(r#""translateable_seq": "{TOY_CODING}", "codon_table": 1"#),
            ),
            transcript(
                TOY_TABLE_2,
                &format!(
                    r#""translateable_seq": "{TOY_CODING}", "three_prime_utr": "{TABLE_2_UTR}", "codon_table": "2""#
                ),
            ),
            transcript(
                TOY_UTR_STOP,
                &format!(
                    r#""translateable_seq": "{TOY_CODING}", "three_prime_utr": "TGA{}""#,
                    &TOY_UTR[3..]
                ),
            ),
            transcript(
                TOY_TABLE_7,
                &format!(
                    r#""translateable_seq": "{TOY_CODING}", "three_prime_utr": "{TOY_UTR}", "codon_table": 7"#
                ),
            ),
            modelled(
                TOY_FORWARD,
                1,
                &[(1, 42, 1001, 1042)],
                TOY_CODING,
                TOY_UTR,
                &[],
                false,
                "",
            ),
            modelled(
                TOY_REVERSE,
                -1,
                &[(1, 42, 2001, 2042)],
                TOY_CODING,
                TOY_UTR,
                &[],
                true,
                "",
            ),
            modelled(
                TOY_TWO_EXONS,
                1,
                &[(1, 9, 3001, 3009), (10, 42, 3014, 3046)],
                TOY_CODING,
                TOY_UTR,
                &[],
                false,
                "",
            ),
            modelled(
                TOY_SEAM,
                1,
                &[(1, 9, 4001, 4009), (13, 42, 4010, 4039)],
                TOY_CODING,
                TOY_UTR,
                &[(4001, 4039)],
                false,
                r#", "peptide": "MKTLAVRG""#,
            ),
            modelled(
                RUN,
                1,
                &[(1, 21, 6001, 6021)],
                RUN_CODING,
                SHORT_UTR,
                &[],
                false,
                "",
            ),
            modelled(
                THREE_EXONS,
                1,
                &[
                    (1, 110, 5001, 5110),
                    (111, 120, 5201, 5210),
                    (121, 130, 5301, 5310),
                ],
                &three_exons_coding,
                "GGGCCCAAAT",
                &[],
                false,
                "",
            ),
        ],
    );
    write_shard(
        &root.0,
        "22",
        &[transcript(
            SHORT,
            &format!(
                r#""translateable_seq": "{SHORT_CODING}", "three_prime_utr": "{SHORT_UTR}", "codon_table": 0"#
            ),
        )],
    );
    root
}

/// The cache's directory and its sequences, which read a chromosome's shards from the
/// directory at the first look-up there.
fn sequences_of(name: &str) -> (Scratch, DownstreamSequences) {
    let cache = downstream_cache(name);
    let sequences = DownstreamSequences::new(&cache.0).expect("a cache with transcripts");
    (cache, sequences)
}

/// The rows of a default-format text of `rows`.
fn parsed(rows: &[String]) -> Vec<Row> {
    parse_default(&format!("{HEADER}{}\n", rows.join("\n"))).expect("rows of fourteen columns")
}

/// A default-format row of one transcript: its coding columns and its Extra text.
fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    coding: [&str; 5],
    extra: &str,
) -> String {
    let [cdna, cds, protein, amino_acids, codons] = coding;
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{consequence}\t{cdna}\t{cds}\t{protein}\t{amino_acids}\t{codons}\t-\t{extra}")
}

/// A frameshift row of the Downstream plugin: its CDS_position, its STRAND and an
/// HGVS_OFFSET where it carries one.
#[derive(Clone)]
struct Frameshift {
    location: &'static str,
    allele: &'static str,
    feature: &'static str,
    cds: &'static str,
    consequence: &'static str,
    strand: i64,
    offset: Option<i64>,
}

fn frameshift(
    location: &'static str,
    allele: &'static str,
    feature: &'static str,
    cds: &'static str,
) -> Frameshift {
    Frameshift {
        location,
        allele,
        feature,
        cds,
        consequence: "frameshift_variant",
        strand: 1,
        offset: None,
    }
}

impl Frameshift {
    fn offset(mut self, offset: i64) -> Frameshift {
        self.offset = Some(offset);
        self
    }

    fn strand(mut self, strand: i64) -> Frameshift {
        self.strand = strand;
        self
    }

    fn terms(mut self, consequence: &'static str) -> Frameshift {
        self.consequence = consequence;
        self
    }

    /// The row with the plugin's two fields, a `None` leaving the field out.
    fn row(&self, protein: Option<&str>, change: Option<&str>) -> Row {
        let mut extra = format!("IMPACT=HIGH;STRAND={}", self.strand);
        if let Some(offset) = self.offset {
            extra += &format!(";HGVS_OFFSET={offset}");
        }
        if let Some(protein) = protein {
            extra += &format!(";DownstreamProtein={protein}");
        }
        if let Some(change) = change {
            extra += &format!(";ProteinLengthChange={change}");
        }
        let row = default_row(
            "v",
            self.location,
            self.allele,
            self.feature,
            self.consequence,
            ["-", self.cds, "-", "-", "-"],
            &extra,
        );
        Row::clone(&parsed(&[row])[0])
    }
}

/// Whether the class `klass` accepts the field `name` of the reference's row `e` and the
/// candidate's `a`, each side's value the field's own or absent.
fn asks(
    klass: &str,
    e: &Row,
    a: &Row,
    name: &str,
    sequences: Option<&DownstreamSequences>,
    record_rows: Option<&[Row]>,
    declared: &[&str],
) -> bool {
    let declared: BTreeSet<String> = declared.iter().map(|name| name.to_string()).collect();
    BuiltIn
        .accepts(&Consult {
            klass,
            e,
            name,
            expected: e.fields.get_or(name, ABSENT),
            actual: a.fields.get_or(name, ABSENT),
            a: Some(a),
            record_rows,
            candidate_rows: None,
            declared: &declared,
            assembly: None,
            source: None,
            gwas_pairs: None,
            revel_table: None,
            utrannotator_sweep: None,
            downstream_sequences: sequences,
            splice_windows: None,
            spliceai_table: None,
            exon_spans: None,
        })
        .expect("no refusal")
}

/// Whether `klass` takes the DownstreamProtein of the row with the reference's value `vep`
/// and the candidate's `vep_rs`, both with the ProteinLengthChange `change`.
fn protein_taken(
    klass: &str,
    row: &Frameshift,
    vep: &str,
    vep_rs: &str,
    change: &str,
    sequences: Option<&DownstreamSequences>,
) -> bool {
    let (e, a) = (
        row.row(Some(vep), Some(change)),
        row.row(Some(vep_rs), Some(change)),
    );
    asks(klass, &e, &a, "DownstreamProtein", sequences, None, &[])
}

#[test]
fn the_genetic_codes_are_those_of_the_codon_module() {
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../crates/vep-core/src/codon.rs");
    let source = fs::read_to_string(source).expect("the codon module");
    let codes = &source[source.find("const CODES:").expect("the tables")..];
    let codes = &codes[..codes.find("\n];").expect("the tables' end")];
    let mut table = 0u8;
    let mut held = 0;
    let mut lines = codes.lines().peekable();
    while let Some(line) = lines.next() {
        if line.starts_with("    None") {
            assert!(genetic_code(table).is_none(), "table {table}");
            table += 1;
        } else if line.starts_with("    code(") {
            let spelt = match line.split_once("b\"") {
                Some((_, rest)) => rest.to_string(),
                None => lines
                    .next()
                    .and_then(|next| next.split_once("b\""))
                    .map(|(_, rest)| rest.to_string())
                    .expect("the table's amino acids"),
            };
            let amino_acids = &spelt[..64];
            assert_eq!(
                genetic_code(table).map(|code| code.as_slice()),
                Some(amino_acids.as_bytes()),
                "table {table}"
            );
            held += 1;
            table += 1;
        }
    }
    assert_eq!(held, 18);
    assert!((table..=u8::MAX).all(|table| genetic_code(table).is_none()));
    // The codons run TTT, TTC, TTA, TTG, TCT ... GGG.
    let standard = genetic_code(1).expect("the standard code");
    assert_eq!(&standard[..5], b"FFLLS");
    assert_eq!(standard[63], b'G');
    let ids = |values: &[Option<Json>]| -> Vec<Option<u8>> {
        values
            .iter()
            .map(|value| codon_table_id(value.as_ref()))
            .collect()
    };
    let (int, text) = (
        |value: i64| Some(Json::int(value)),
        |value: &str| Some(Json::str(value)),
    );
    assert_eq!(
        ids(&[
            None,
            int(0),
            text("0"),
            int(1),
            text("1"),
            int(2),
            text("2"),
            int(24)
        ]),
        [1, 1, 1, 1, 1, 2, 2, 24].map(Some)
    );
    assert_eq!(
        ids(&[
            int(7),
            text("8"),
            int(17),
            int(20),
            int(25),
            text("x"),
            Some(Json::Bool(true))
        ]),
        [None; 7]
    );
    assert_eq!(
        codon_table_id(Some(&Json::Null)),
        Some(1),
        "a null is no key"
    );
}

#[test]
fn a_sequence_is_sliced_as_python_slices_a_str() {
    let sequence = Sequence::new("ACGTAC".to_string());
    assert_eq!(sequence.len(), 6);
    assert_eq!(sequence.slice(None, Some(2)), "AC");
    assert_eq!(sequence.slice(Some(4), None), "AC");
    assert_eq!(
        sequence.slice(None, Some(-1)),
        "ACGTA",
        "an end counted from the end"
    );
    assert_eq!(
        sequence.slice(Some(-3), Some(0)),
        "",
        "a start past the end"
    );
    assert_eq!(sequence.slice(Some(-6), Some(-3)), "ACG");
    assert_eq!(sequence.slice(Some(9), None), "");
    assert_eq!(sequence.slice(None, Some(9)), "ACGTAC");
    let accented = Sequence::new("AÇGTAC".to_string());
    assert_eq!(accented.len(), 6, "characters, not bytes");
    assert_eq!(accented.slice(Some(1), Some(3)), "ÇG");
    assert_eq!(accented.slice(Some(2), None), "GTAC");
}

#[test]
fn downstream_utr_frame_past_coding_end_accepts_a_deletion_running_to_the_last_coding_base() {
    let (_cache, sequences) = sequences_of("utr-frame");
    let versioned: &'static str = "ENST00000000101.3";
    let rows = [
        (
            frameshift("21:124-127", "-", TOY, "24-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            frameshift("21:124-127", "-", versioned, "24-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            frameshift("21:1024-1027", "-", TOY_FORWARD, "24-27").offset(2),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            frameshift("21:2016-2019", "-", TOY_REVERSE, "24-27")
                .offset(-2)
                .strand(-1),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            frameshift("21:224-227", "-", TOY_UTR_STOP, "24-27"),
            "",
            "GDPNS",
            "4",
        ),
        (
            frameshift("21:324-327", "-", TOY_TABLE_2, "24-27"),
            "VTKMA",
            "GWP",
            "2",
        ),
    ];
    for (row, vep, vep_rs, change) in &rows {
        let taken =
            |klass: &str, sequences| protein_taken(klass, row, vep, vep_rs, change, sequences);
        assert!(
            taken(UTR_FRAME, Some(&sequences)),
            "{} {}",
            row.location,
            row.feature
        );
        assert!(!taken(UTR_FRAME, None), "without the cache");
        assert!(
            !taken(START, Some(&sequences)),
            "the other Downstream class"
        );
        assert!(
            !taken(SHIFTED, Some(&sequences)),
            "the class of the shifted pair"
        );
    }
}

#[test]
fn downstream_utr_frame_past_coding_end_leaves_a_row_the_mechanism_does_not_reach_a_difference() {
    let (_cache, sequences) = sequences_of("utr-frame-refused");
    let own = frameshift("21:124-127", "-", TOY, "24-27");
    let taken = |row: &Frameshift, vep: &str, vep_rs: &str, change: &str| {
        protein_taken(UTR_FRAME, row, vep, vep_rs, change, Some(&sequences))
    };
    assert!(taken(&own, "GPK", "GGPNS", "4"));
    let refused = [
        (
            "ends before the last coding base",
            frameshift("21:121-124", "-", TOY, "21-24"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "starts on a codon's first base",
            frameshift("21:125-127", "-", TOY, "25-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "starts on a codon's second base",
            frameshift("21:126-127", "-", TOY, "26-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "the column two bases before without an offset",
            frameshift("21:122-125", "-", TOY, "22-25"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "the column two bases before with the offset",
            frameshift("21:1022-1025", "-", TOY_FORWARD, "22-25").offset(2),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "inserted bases",
            frameshift("21:124-127", "GC", TOY, "24-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "no 3' UTR",
            frameshift("21:124-127", "-", TOY_NO_UTR, "24-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "the standard code's pair under table 2",
            frameshift("21:124-127", "-", TOY_TABLE_2, "24-27"),
            "VTKIA",
            "G",
            "0",
        ),
        (
            "a table the module lacks",
            frameshift("21:124-127", "-", TOY_TABLE_7, "24-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "absent from the cache",
            frameshift("21:124-127", "-", "ENST00000000199", "24-27"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "a planted candidate value",
            own.clone(),
            "GPK",
            "GGPNT",
            "4",
        ),
        (
            "a planted reference value",
            own.clone(),
            "GPR",
            "GGPNS",
            "4",
        ),
        (
            "no frameshift_variant",
            own.clone().terms("inframe_deletion"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "a splice term",
            own.clone()
                .terms("frameshift_variant,splice_region_variant"),
            "GPK",
            "GGPNS",
            "4",
        ),
        (
            "a deletion ending inside the coding sequence",
            frameshift("21:106-107", "-", TOY, "6-7"),
            "SCCPWLRAQIA",
            "NSCCPWLRAQIA",
            "5",
        ),
    ];
    for (label, row, vep, vep_rs, change) in &refused {
        assert!(!taken(row, vep, vep_rs, change), "{label}");
    }
    let ask = |e: &Row, a: &Row| {
        asks(
            UTR_FRAME,
            e,
            a,
            "DownstreamProtein",
            Some(&sequences),
            None,
            &[],
        )
    };
    let (e, a) = (
        own.row(Some("GPK"), Some("4")),
        own.row(Some("GGPNS"), Some("4")),
    );
    assert!(
        !ask(&own.row(None, Some("4")), &a),
        "the reference lacks the field"
    );
    assert!(
        !ask(&e, &own.row(None, Some("4"))),
        "the candidate lacks the field"
    );
    assert!(
        !ask(&e, &own.row(Some("GGPNS"), Some("5"))),
        "a differing ProteinLengthChange"
    );
    let other_span = frameshift("21:124-127", "-", TOY, "23-27");
    assert!(
        !ask(&e, &other_span.row(Some("GGPNS"), Some("4"))),
        "a differing CDS_position"
    );
    let more_terms = own.clone().terms("frameshift_variant,stop_lost");
    assert!(
        !ask(&e, &more_terms.row(Some("GGPNS"), Some("4"))),
        "differing consequence sets"
    );
}

#[test]
fn a_transcript_filed_under_two_chromosomes_is_read_by_an_entry_without_a_location() {
    let filed = |root: &Path, chromosome: &str, start: i64, coding: &str| {
        write_shard(
            root,
            chromosome,
            &[modelled(
                TOY,
                1,
                &[(1, 42, start, start + 41)],
                coding,
                TOY_UTR,
                &[],
                false,
                "",
            )],
        );
    };
    let (two_models, two_sequences) = (Scratch::new("two-models"), Scratch::new("two-sequences"));
    filed(&two_models.0, "X", 1001, TOY_CODING);
    filed(&two_models.0, "Y", 501, TOY_CODING);
    filed(&two_sequences.0, "X", 1001, TOY_CODING);
    filed(
        &two_sequences.0,
        "Y",
        501,
        &format!("{}TGA", &TOY_CODING[..TOY_CODING.len() - 3]),
    );
    let located = frameshift("X:1024-1027", "-", TOY, "24-27").row(Some("GPK"), Some("4"));
    let unlocated: Row = Row::new(Entry {
        location: None,
        ..Entry::clone(&located)
    });
    let two_models = DownstreamSequences::new(&two_models.0).expect("a cache");
    let model_of = |row: &Row| {
        let transcript = two_models.transcript(row).expect("shards that read");
        transcript.expect("the transcript").model.is_some()
    };
    assert!(model_of(&located), "the chromosome's own model");
    assert!(!model_of(&unlocated), "one reading under two models");
    // The row is the class's either way: the span runs to the last coding base.
    let candidate = frameshift("X:1024-1027", "-", TOY, "24-27").row(Some("GGPNS"), Some("4"));
    for row in [&located, &unlocated] {
        assert!(asks(
            UTR_FRAME,
            row,
            &candidate,
            "DownstreamProtein",
            Some(&two_models),
            None,
            &[]
        ));
    }
    let two_sequences = DownstreamSequences::new(&two_sequences.0).expect("a cache");
    let held = |row: &Row| {
        two_sequences
            .transcript(row)
            .expect("shards that read")
            .is_some()
    };
    assert!(held(&located));
    assert!(
        !held(&unlocated),
        "two coding sequences under one stable id"
    );
}

#[test]
fn a_cache_without_a_transcripts_directory_is_refused() {
    let nowhere = Scratch::new("nowhere");
    let refusal = DownstreamSequences::new(&nowhere.0)
        .err()
        .expect("a refusal");
    assert_eq!(
        refusal.0,
        format!(
            "--json-cache: {} has no transcripts/ directory",
            nowhere.0.display()
        )
    );
}

#[test]
fn downstream_start_past_own_codon_accepts_a_translation_started_one_codon_late() {
    let (_cache, sequences) = sequences_of("start");
    let rows = [
        (
            frameshift("21:106", "-", TOY, "6"),
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            "5",
        ),
        (
            frameshift("21:106-107", "-", TOY, "6-7"),
            "SCCPWLRAQIA",
            "NSCCPWLRAQIA",
            "5",
        ),
        (
            frameshift("21:106", "GC", TOY, "6"),
            "HSCCPWLRAQIA",
            "KHSCCPWLRAQIA",
            "6",
        ),
        (frameshift("22:206", "-", SHORT, "6"), "KVK", "", "-3"),
        (
            frameshift("21:206", "-", TOY_NO_UTR, "6"),
            "LLLSVV",
            "KLLLSVV",
            "0",
        ),
        (
            frameshift("21:306", "AC", TOY, "6").strand(-1),
            "YSCCPWLRAQIA",
            "KYSCCPWLRAQIA",
            "6",
        ),
        (
            frameshift("21:1006", "-", TOY_FORWARD, "6").offset(1),
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            "5",
        ),
        (
            frameshift("21:2037", "-", TOY_REVERSE, "6")
                .offset(-1)
                .strand(-1),
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            "5",
        ),
        (
            frameshift("21:406", "-", TOY_TABLE_2, "6"),
            "LLLSVVKWP",
            "KLLLSVVKWP",
            "3",
        ),
    ];
    for (row, vep, vep_rs, change) in &rows {
        let taken =
            |klass: &str, sequences| protein_taken(klass, row, vep, vep_rs, change, sequences);
        assert!(
            taken(START, Some(&sequences)),
            "{} {}",
            row.location,
            row.feature
        );
        assert!(!taken(START, None), "without the cache");
        assert!(
            !taken(UTR_FRAME, Some(&sequences)),
            "the other Downstream class"
        );
    }
}

#[test]
fn downstream_start_past_own_codon_leaves_a_row_the_mechanism_does_not_reach_a_difference() {
    let (_cache, sequences) = sequences_of("start-refused");
    let own = frameshift("21:106", "-", TOY, "6");
    let short = frameshift("22:206", "-", SHORT, "6");
    let taken = |row: &Frameshift, vep: &str, vep_rs: &str, change: &str| {
        protein_taken(START, row, vep, vep_rs, change, Some(&sequences))
    };
    assert!(taken(&own, "LLLSVVKGPNS", "KLLLSVVKGPNS", "5"));
    assert!(taken(&short, "KVK", "", "-3"));
    let pair = ("LLLSVVKGPNS", "KLLLSVVKGPNS", "5");
    let refused = [
        (
            "starts on a codon's first base",
            frameshift("21:107", "-", TOY, "7"),
            pair,
        ),
        (
            "starts on a codon's second base",
            frameshift("21:105", "-", TOY, "5"),
            pair,
        ),
        (
            "the column two bases before without an offset",
            frameshift("21:104", "-", TOY, "4"),
            pair,
        ),
        (
            "the column two bases before with the offset",
            frameshift("21:1004", "-", TOY_FORWARD, "4").offset(2),
            pair,
        ),
        (
            "an insertion carrying a replacement's pair",
            frameshift("21:106-107", "GC", TOY, "6-7"),
            ("HSCCPWLRAQIA", "KHSCCPWLRAQIA", "6"),
        ),
        (
            "an insertion's value against a planted one",
            frameshift("21:106-107", "GC", TOY, "6-7"),
            ("ALLLSVVKGPNS", "KALLLSVVKGPNS", "6"),
        ),
        (
            "a reverse-strand row's pair on a forward-strand row",
            frameshift("21:106", "AC", TOY, "6"),
            ("YSCCPWLRAQIA", "KYSCCPWLRAQIA", "6"),
        ),
        (
            "the standard code's pair under table 2",
            frameshift("21:106", "-", TOY_TABLE_2, "6"),
            ("LLLSVVK", "KLLLSVVK", "1"),
        ),
        (
            "a table the module lacks",
            frameshift("21:106", "-", TOY_TABLE_7, "6"),
            pair,
        ),
        (
            "absent from the cache",
            frameshift("21:106", "-", "ENST00000000199", "6"),
            pair,
        ),
        (
            "a planted candidate value",
            own.clone(),
            ("LLLSVVKGPNS", "KLLLSVVKGPNT", "5"),
        ),
        (
            "a planted reference value",
            own.clone(),
            ("LLLSVVKGPNT", "KLLLSVVKGPNS", "5"),
        ),
        (
            "no frameshift_variant",
            own.clone().terms("inframe_deletion"),
            pair,
        ),
        (
            "a splice term",
            own.clone().terms("frameshift_variant,splice_donor_variant"),
            pair,
        ),
        (
            "a deletion running to the last coding base",
            frameshift("21:124-127", "-", TOY, "24-27"),
            ("GPK", "GGPNS", "4"),
        ),
    ];
    for (label, row, (vep, vep_rs, change)) in &refused {
        assert!(!taken(row, vep, vep_rs, change), "{label}");
    }
    let ask = |e: &Row, a: &Row| {
        asks(
            START,
            e,
            a,
            "DownstreamProtein",
            Some(&sequences),
            None,
            &[],
        )
    };
    let (e, a) = (
        own.row(Some("LLLSVVKGPNS"), Some("5")),
        own.row(Some("KLLLSVVKGPNS"), Some("5")),
    );
    assert!(
        !ask(&own.row(None, Some("5")), &a),
        "the reference lacks the field"
    );
    assert!(
        !ask(&e, &own.row(None, Some("5"))),
        "the candidate lacks the field"
    );
    assert!(
        !ask(
            &short.row(Some("KVK"), Some("-3")),
            &short.row(None, Some("-3"))
        ),
        "the candidate lacks the field where its value is empty"
    );
    assert!(
        !ask(&e, &own.row(Some("KLLLSVVKGPNS"), Some("4"))),
        "a differing ProteinLengthChange"
    );
}

#[test]
fn downstream_classes_take_the_start_term_difference_the_registry_sets_aside() {
    let (_cache, sequences) = sequences_of("start-terms");
    let both_start_terms = "frameshift_variant,start_lost,start_retained_variant";
    let without_start_lost = "frameshift_variant,start_retained_variant";
    let short = (frameshift("22:206", "-", SHORT, "6"), ("KVK", ""), "-3");
    let coding_end = (
        frameshift("21:124-127", "-", TOY, "24-27"),
        ("GPK", "GGPNS"),
        "4",
    );
    let shifted = (
        frameshift("21:1005", "-", TOY_FORWARD, "5").offset(2),
        ("LLLSVVKGPNS", "KLLLSVVKGPNS"),
        "5",
    );
    let judged = |reference_terms: &'static str,
                  candidate_terms: &'static str,
                  (row, values, change): &(Frameshift, (&str, &str), &str)| {
        let e = row
            .clone()
            .terms(reference_terms)
            .row(Some(values.0), Some(change));
        let a = row
            .clone()
            .terms(candidate_terms)
            .row(Some(values.1), Some(change));
        let accepting: Vec<&str> = [UTR_FRAME, START, SHIFTED]
            .into_iter()
            .filter(|klass| {
                asks(
                    klass,
                    &e,
                    &a,
                    "DownstreamProtein",
                    Some(&sequences),
                    None,
                    &[],
                )
            })
            .collect();
        let found = downstream_class(
            &e,
            &a,
            "DownstreamProtein",
            values.0,
            values.1,
            &sequences,
            None,
        )
        .expect("no refusal");
        // A row one cause changes: the class alone, none to declare beside it.
        assert_eq!(
            accepting,
            found
                .map(|(klass, _)| klass)
                .into_iter()
                .collect::<Vec<_>>()
        );
        assert!(found.is_none_or(|(_, beside)| beside.is_none()));
        found.map(|(klass, _)| klass)
    };
    assert_eq!(
        judged("frameshift_variant", "frameshift_variant", &short),
        Some(START)
    );
    assert_eq!(
        judged(both_start_terms, without_start_lost, &short),
        Some(START)
    );
    assert_eq!(
        judged(both_start_terms, without_start_lost, &coding_end),
        Some(UTR_FRAME)
    );
    assert_eq!(
        judged("frameshift_variant", "frameshift_variant", &shifted),
        Some(SHIFTED)
    );
    assert_eq!(
        judged(both_start_terms, without_start_lost, &shifted),
        Some(SHIFTED)
    );
    for (reference_terms, candidate_terms) in [
        (both_start_terms, "frameshift_variant,start_lost"),
        (both_start_terms, "frameshift_variant"),
        (
            both_start_terms,
            "frameshift_variant,start_retained_variant,stop_lost",
        ),
        ("frameshift_variant,start_lost", "frameshift_variant"),
        (without_start_lost, both_start_terms),
        (
            "frameshift_variant,start_lost,start_retained_variant,splice_region_variant",
            "frameshift_variant,start_retained_variant,splice_region_variant",
        ),
    ] {
        assert_eq!(
            judged(reference_terms, candidate_terms, &short),
            None,
            "{candidate_terms}"
        );
        assert_eq!(
            judged(reference_terms, candidate_terms, &shifted),
            None,
            "{candidate_terms}"
        );
    }
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_leaves_a_downstream_protein_the_shift_does_not_explain_a_difference(
) {
    let (_cache, sequences) = sequences_of("shifted-protein");
    let own = frameshift("21:1005", "-", TOY_FORWARD, "5").offset(2);
    let pair = ("LLLSVVKGPNS", "KLLLSVVKGPNS", "5");
    let taken = |row: &Frameshift, (vep, vep_rs, change): (&str, &str, &str)| {
        protein_taken(SHIFTED, row, vep, vep_rs, change, Some(&sequences))
    };
    assert!(taken(&own, pair));
    assert!(
        !protein_taken(SHIFTED, &own, pair.0, pair.1, pair.2, None),
        "without the cache"
    );
    let no_offset = frameshift("21:1005", "-", TOY_FORWARD, "5");
    let coding_end = ("GPK", "GGPNS", "4");
    let refused = [
        ("no HGVS_OFFSET", no_offset.clone(), pair),
        ("an HGVS_OFFSET of 0", no_offset.clone().offset(0), pair),
        (
            "a shift past the last coding base",
            no_offset.clone().offset(23),
            pair,
        ),
        (
            "a shift off the transcript",
            no_offset.clone().offset(60),
            pair,
        ),
        (
            "the offset of the other strand",
            no_offset.clone().offset(-2),
            pair,
        ),
        (
            "a shifted pair read as the own pair is",
            frameshift("21:1004", "-", TOY_FORWARD, "4").offset(1),
            pair,
        ),
        (
            "a Location that does not map to the CDS_position",
            frameshift("21:1006", "-", TOY_FORWARD, "5").offset(2),
            pair,
        ),
        (
            "a transcript without a model",
            frameshift("21:105", "-", TOY, "5").offset(2),
            pair,
        ),
        (
            "absent from the cache",
            frameshift("21:1005", "-", "ENST00000000199", "5").offset(2),
            pair,
        ),
        (
            "a planted candidate value",
            own.clone(),
            ("LLLSVVKGPNS", "KLLLSVVKGPNT", "5"),
        ),
        (
            "a planted reference value",
            own.clone(),
            ("LLLSVVKGPNT", "KLLLSVVKGPNS", "5"),
        ),
        (
            "no frameshift_variant",
            own.clone().terms("inframe_deletion"),
            pair,
        ),
        (
            "a splice term",
            own.clone()
                .terms("frameshift_variant,splice_region_variant"),
            pair,
        ),
        (
            "the start class's row",
            frameshift("21:106", "-", TOY, "6"),
            pair,
        ),
        (
            "the start class's row beside an offset",
            frameshift("21:1006", "-", TOY_FORWARD, "6").offset(1),
            pair,
        ),
        (
            "the coding end class's row",
            frameshift("21:124-127", "-", TOY, "24-27"),
            coding_end,
        ),
        (
            "the coding end class's row beside an offset",
            frameshift("21:1024-1027", "-", TOY_FORWARD, "24-27").offset(2),
            coding_end,
        ),
    ];
    for (label, row, values) in &refused {
        assert!(!taken(row, *values), "{label}");
    }
    let ask = |e: &Row, a: &Row| {
        asks(
            SHIFTED,
            e,
            a,
            "DownstreamProtein",
            Some(&sequences),
            None,
            &[],
        )
    };
    let (e, a) = (
        own.row(Some(pair.0), Some("5")),
        own.row(Some(pair.1), Some("5")),
    );
    assert!(
        !ask(&own.row(None, Some("5")), &a),
        "the reference lacks the field"
    );
    assert!(
        !ask(&e, &own.row(None, Some("5"))),
        "the candidate lacks the field"
    );
    let other_span = frameshift("21:1005", "-", TOY_FORWARD, "4").offset(2);
    assert!(
        !ask(&e, &other_span.row(Some(pair.1), Some("5"))),
        "a differing CDS_position"
    );
    let more_terms = own.clone().terms("frameshift_variant,stop_lost");
    assert!(
        !ask(&e, &more_terms.row(Some(pair.1), Some("5"))),
        "differing consequence sets"
    );
    // The insertion after base 9, the first exon's last, moved five bases to the point
    // after base 10: the row is placed by its Location.
    let seam = frameshift("21:3009-3010", "A", TOY_TWO_EXONS, "9-10").offset(5);
    assert!(taken(&seam, ("HCCPWLRAQIA", "TCCPWLRAQIA", "6")));
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_accepts_a_protein_length_change_read_at_the_shifted_coordinates(
) {
    let (_cache, sequences) = sequences_of("shifted-length");
    let seam = frameshift("21:4006-4007", "A", TOY_SEAM, "6-7").offset(3);
    let e = seam.row(Some("SCPWLRAQIA"), Some("5"));
    let a = seam.row(Some("NSCCPWLRAQIA"), Some("6"));
    for name in ["DownstreamProtein", "ProteinLengthChange"] {
        assert!(
            asks(SHIFTED, &e, &a, name, Some(&sequences), None, &[]),
            "{name}"
        );
        assert!(
            !asks(SHIFTED, &e, &a, name, None, None, &[]),
            "{name} without the cache"
        );
        for other in [UTR_FRAME, START] {
            assert!(
                !asks(other, &e, &a, name, Some(&sequences), None, &[]),
                "{other} {name}"
            );
        }
    }
    let agreed = seam.row(Some("NSCCPWLRAQIA"), Some("5"));
    assert!(asks(
        SHIFTED,
        &agreed,
        &a,
        "ProteinLengthChange",
        Some(&sequences),
        None,
        &[]
    ));
    // Each field is judged on its own values: beside a planted ProteinLengthChange the
    // DownstreamProtein stays the class's.
    let planted = seam.row(Some("NSCCPWLRAQIA"), Some("7"));
    assert!(asks(
        SHIFTED,
        &e,
        &planted,
        "DownstreamProtein",
        Some(&sequences),
        None,
        &[]
    ));
    assert!(!asks(
        SHIFTED,
        &e,
        &planted,
        "ProteinLengthChange",
        Some(&sequences),
        None,
        &[]
    ));
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_leaves_a_protein_length_change_the_shift_does_not_explain_a_difference(
) {
    let (_cache, sequences) = sequences_of("shifted-length-refused");
    let seam = frameshift("21:4006-4007", "A", TOY_SEAM, "6-7");
    let protein = "NSCCPWLRAQIA";
    let taken = |row: &Frameshift, protein: &str, vep: Option<&str>, vep_rs: Option<&str>| {
        let (e, a) = (row.row(Some(protein), vep), row.row(Some(protein), vep_rs));
        asks(
            SHIFTED,
            &e,
            &a,
            "ProteinLengthChange",
            Some(&sequences),
            None,
            &[],
        )
    };
    let own = seam.clone().offset(3);
    assert!(taken(&own, protein, Some("5"), Some("6")));
    let refused = [
        (
            "no HGVS_OFFSET",
            seam.clone(),
            protein,
            Some("5"),
            Some("6"),
        ),
        (
            "a shift past the last coding base",
            seam.clone().offset(22),
            protein,
            Some("5"),
            Some("6"),
        ),
        (
            "a shifted point read as the own point is",
            seam.clone().offset(1),
            protein,
            Some("5"),
            Some("6"),
        ),
        (
            "a planted candidate value",
            own.clone(),
            protein,
            Some("5"),
            Some("7"),
        ),
        (
            "a planted reference value",
            own.clone(),
            protein,
            Some("4"),
            Some("6"),
        ),
        (
            "a value that is no integer",
            own.clone(),
            protein,
            Some("five"),
            Some("6"),
        ),
        (
            "the reference lacks the field",
            own.clone(),
            protein,
            None,
            Some("6"),
        ),
        (
            "the candidate lacks the field",
            own.clone(),
            protein,
            Some("5"),
            None,
        ),
        (
            "absent from the cache",
            frameshift("21:4006-4007", "A", "ENST00000000199", "6-7").offset(3),
            protein,
            Some("5"),
            Some("6"),
        ),
        (
            "no frameshift_variant",
            own.clone().terms("inframe_insertion"),
            protein,
            Some("5"),
            Some("6"),
        ),
        (
            "a splice term",
            own.clone()
                .terms("frameshift_variant,splice_region_variant"),
            protein,
            Some("5"),
            Some("6"),
        ),
        (
            "the start class's row",
            frameshift("21:106", "-", TOY, "6"),
            "KLLLSVVKGPNS",
            Some("5"),
            Some("4"),
        ),
    ];
    for (label, row, protein, vep, vep_rs) in &refused {
        assert!(!taken(row, protein, *vep, *vep_rs), "{label}");
    }
    let e = own.row(Some(protein), Some("5"));
    let more_terms = own.clone().terms("frameshift_variant,stop_lost");
    let a = more_terms.row(Some(protein), Some("6"));
    assert!(
        !asks(
            SHIFTED,
            &e,
            &a,
            "ProteinLengthChange",
            Some(&sequences),
            None,
            &[]
        ),
        "differing consequence sets"
    );
}

/// A row of the NMD plugin: its CDS_position, an HGVS_OFFSET where `offset` gives one and
/// the plugin's field where `nmd` gives its value.
fn nmd_row(
    row: (&str, &str, &str, &str),
    nmd: Option<&str>,
    consequence: &str,
    offset: Option<i64>,
    cdna: &str,
) -> Row {
    let (location, allele, feature, cds) = row;
    let mut extra = "IMPACT=HIGH;STRAND=1".to_string();
    if let Some(offset) = offset {
        extra += &format!(";HGVS_OFFSET={offset}");
    }
    if let Some(nmd) = nmd {
        extra += &format!(";NMD={nmd}");
    }
    let text = default_row(
        "v",
        location,
        allele,
        feature,
        consequence,
        [cdna, cds, "-", "-", "-"],
        &extra,
    );
    Row::clone(&parsed(&[text])[0])
}

const NMD_ACCEPTOR_ROW: (&str, &str, &str, &str) = ("21:3013-3014", "-", TOY_TWO_EXONS, "?-10");
const NMD_PENULTIMATE_ROW: (&str, &str, &str, &str) = ("21:5200-5201", "-", THREE_EXONS, "?-111");
const NMD_LATE_ROW: (&str, &str, &str, &str) = ("21:5100", "-", THREE_EXONS, "100");

/// Whether `plugin_hgvs_shifted_cds_coordinates` takes the NMD of `row` that the
/// reference's entry writes as `vep` under `reference_terms` and the candidate's as
/// `vep_rs` under `candidate_terms`: the class's verdict, held to the rule's own.
#[allow(clippy::too_many_arguments)]
fn nmd_taken(
    reference_terms: &str,
    candidate_terms: &str,
    sequences: Option<&DownstreamSequences>,
    row: (&str, &str, &str, &str),
    offset: Option<i64>,
    vep: Option<&str>,
    vep_rs: Option<&str>,
    cdna: &str,
) -> bool {
    let e = nmd_row(row, vep, reference_terms, offset, cdna);
    let a = nmd_row(row, vep_rs, candidate_terms, offset, cdna);
    let taken = asks(SHIFTED, &e, &a, "NMD", sequences, None, &[]);
    let read = sequences.is_some_and(|sequences| {
        nmd_reads_shifted_pair(
            &e,
            &a,
            vep.unwrap_or(ABSENT),
            vep_rs.unwrap_or(ABSENT),
            sequences,
            None,
        )
        .expect("no refusal")
    });
    assert_eq!(taken, read);
    taken
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_accepts_an_nmd_read_at_the_shifted_coordinates() {
    let (_cache, sequences) = sequences_of("nmd");
    let frameshift = "frameshift_variant";
    for (row, terms, offset, vep, vep_rs) in [
        (NMD_ACCEPTOR_ROW, ACCEPTOR_TERMS, 1, Some(NMD_VALUE), None),
        (
            NMD_PENULTIMATE_ROW,
            ACCEPTOR_TERMS,
            1,
            Some(NMD_VALUE),
            None,
        ),
        (NMD_LATE_ROW, frameshift, 3, None, Some(NMD_VALUE)),
    ] {
        assert!(
            nmd_taken(
                terms,
                terms,
                Some(&sequences),
                row,
                Some(offset),
                vep,
                vep_rs,
                "-"
            ),
            "{}",
            row.0
        );
        assert!(
            !nmd_taken(terms, terms, None, row, Some(offset), vep, vep_rs, "-"),
            "without the cache"
        );
    }
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_leaves_an_nmd_the_shift_does_not_explain_a_difference() {
    let (_cache, sequences) = sequences_of("nmd-refused");
    let frameshift = "frameshift_variant";
    let taken = |row, terms: &str, offset: Option<i64>, vep: Option<&str>, vep_rs: Option<&str>| {
        nmd_taken(
            terms,
            terms,
            Some(&sequences),
            row,
            offset,
            vep,
            vep_rs,
            "-",
        )
    };
    assert!(taken(
        NMD_LATE_ROW,
        frameshift,
        Some(3),
        None,
        Some(NMD_VALUE)
    ));
    let refused = [
        (
            "no HGVS_OFFSET",
            NMD_LATE_ROW,
            frameshift,
            None,
            None,
            Some(NMD_VALUE),
        ),
        (
            "a shift that leaves the exon",
            NMD_LATE_ROW,
            frameshift,
            Some(30),
            None,
            Some(NMD_VALUE),
        ),
        (
            "a shifted pair the rule reads as the own pair",
            NMD_LATE_ROW,
            frameshift,
            Some(1),
            None,
            Some(NMD_VALUE),
        ),
        (
            "the value on the reference's side",
            NMD_LATE_ROW,
            frameshift,
            Some(3),
            Some(NMD_VALUE),
            None,
        ),
        (
            "another text on the reference's side",
            NMD_LATE_ROW,
            frameshift,
            Some(3),
            Some("NMD"),
            Some(NMD_VALUE),
        ),
        (
            "an empty value on the reference's side",
            NMD_LATE_ROW,
            frameshift,
            Some(3),
            Some(""),
            Some(NMD_VALUE),
        ),
        (
            "a text that is not the plugin's",
            NMD_LATE_ROW,
            frameshift,
            Some(3),
            None,
            Some("escaping"),
        ),
        (
            "absent from the cache",
            ("21:5100", "-", "ENST00000000199", "100"),
            frameshift,
            Some(3),
            None,
            Some(NMD_VALUE),
        ),
        (
            "a transcript without a model",
            ("21:105", "-", TOY, "5"),
            frameshift,
            Some(3),
            None,
            Some(NMD_VALUE),
        ),
        (
            "no term of the plugin's",
            NMD_LATE_ROW,
            "inframe_deletion",
            Some(3),
            None,
            Some(NMD_VALUE),
        ),
        (
            "the acceptor row's value on the candidate's side",
            NMD_ACCEPTOR_ROW,
            ACCEPTOR_TERMS,
            Some(1),
            None,
            Some(NMD_VALUE),
        ),
    ];
    for (label, row, terms, offset, vep, vep_rs) in refused {
        assert!(!taken(row, terms, offset, vep, vep_rs), "{label}");
    }
    let e = nmd_row(NMD_LATE_ROW, None, frameshift, Some(3), "-");
    let other_span = nmd_row(
        ("21:5100", "-", THREE_EXONS, "99"),
        Some(NMD_VALUE),
        frameshift,
        Some(3),
        "-",
    );
    assert!(
        !asks(SHIFTED, &e, &other_span, "NMD", Some(&sequences), None, &[]),
        "a differing CDS_position"
    );
    assert!(
        !nmd_taken(
            frameshift,
            "frameshift_variant,stop_lost",
            Some(&sequences),
            NMD_LATE_ROW,
            Some(3),
            None,
            Some(NMD_VALUE),
            "-"
        ),
        "differing consequence sets"
    );
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_takes_an_nmd_across_the_retained_stop_term_pair() {
    let (_cache, sequences) = sequences_of("nmd-retained");
    let pairs = [
        (
            "intron_variant,splice_acceptor_variant,splice_donor_5th_base_variant,splice_donor_variant,stop_retained_variant",
            "coding_sequence_variant,intron_variant,splice_acceptor_variant,splice_donor_5th_base_variant,splice_donor_variant",
        ),
        (
            "splice_acceptor_variant,stop_retained_variant",
            "coding_sequence_variant,splice_acceptor_variant",
        ),
    ];
    let taken = |reference_terms: &str, candidate_terms: &str, sequences, row| {
        nmd_taken(
            reference_terms,
            candidate_terms,
            sequences,
            row,
            Some(1),
            Some(NMD_VALUE),
            None,
            "-",
        )
    };
    for (reference_terms, candidate_terms) in pairs {
        for row in [NMD_ACCEPTOR_ROW, NMD_PENULTIMATE_ROW] {
            assert!(
                taken(reference_terms, candidate_terms, Some(&sequences), row),
                "{}",
                row.0
            );
            assert!(
                !taken(reference_terms, candidate_terms, None, row),
                "{}",
                row.0
            );
        }
    }
    let (retained, coding) = pairs[1];
    let on_the_acceptor_row = |reference_terms: &str, candidate_terms: &str| {
        taken(
            reference_terms,
            candidate_terms,
            Some(&sequences),
            NMD_ACCEPTOR_ROW,
        )
    };
    assert!(on_the_acceptor_row(retained, coding));
    assert!(on_the_acceptor_row(
        retained,
        "frameshift_variant,splice_acceptor_variant"
    ));
    let start_pair = (
        "splice_acceptor_variant,start_lost,start_retained_variant",
        "splice_acceptor_variant,start_retained_variant",
    );
    for (label, reference_terms, candidate_terms) in [
        ("the two terms the other way round", coding, retained),
        (
            "stop_lost in the term's place",
            retained,
            "splice_acceptor_variant,stop_lost",
        ),
        (
            "a further term on the candidate's side",
            retained,
            "coding_sequence_variant,splice_acceptor_variant,intron_variant",
        ),
        (
            "a further term on the reference's side",
            "splice_acceptor_variant,stop_retained_variant,intron_variant",
            coding,
        ),
        ("the registry's start pair", start_pair.0, start_pair.1),
        (
            "no term of the plugin's in the reference's set",
            "stop_retained_variant",
            "frameshift_variant",
        ),
        (
            "no term of the plugin's in either set",
            "intron_variant,stop_retained_variant",
            "coding_sequence_variant,intron_variant",
        ),
    ] {
        assert!(
            !on_the_acceptor_row(reference_terms, candidate_terms),
            "{label}"
        );
    }
    assert_eq!(
        registry_class(start_pair.0, start_pair.1, "-", &[START_COOCCURRENCE_SWAP]),
        Some(START_COOCCURRENCE_SWAP)
    );
    // The deletion of the three-exon transcript's coding base 100, whose value the
    // candidate alone writes: without a cDNA position the row shows a stop codon test that
    // is not made, and with one beside its CDS start a test both engines make.
    let late = |cdna: &str| {
        nmd_taken(
            "frameshift_variant,stop_retained_variant",
            "frameshift_variant",
            Some(&sequences),
            NMD_LATE_ROW,
            Some(3),
            None,
            Some(NMD_VALUE),
            cdna,
        )
    };
    assert!(late("-"));
    assert!(!late("300"));
}

/// A frameshift row of a made-up record on the toy forward transcript that inserts after
/// coding base 4: its coding columns, an HGVS_OFFSET where `offset` gives one and the
/// Downstream plugin's two fields where `downstream` gives the protein.
fn later_allele_row(
    allele: &str,
    cds: &str,
    protein: &str,
    amino_acids: &str,
    downstream: Option<&str>,
    offset: Option<i64>,
) -> String {
    let mut extra = "IMPACT=HIGH;STRAND=1".to_string();
    if let Some(offset) = offset {
        extra += &format!(";HGVS_OFFSET={offset}");
    }
    if let Some(downstream) = downstream {
        extra += &format!(";DownstreamProtein={downstream};ProteinLengthChange=6");
    }
    default_row(
        "v",
        "21:1004-1005",
        allele,
        TOY_FORWARD,
        "frameshift_variant",
        ["-", cds, protein, amino_acids, "-"],
        &extra,
    )
}

/// Whether `klass` takes the field `name` of the last row of the reference's rows against
/// the last of the candidate's, the reference's rows read as those of the input records
/// `inputs`.
fn last_row_taken(
    klass: &str,
    name: &str,
    reference: &[String],
    candidate: &[String],
    sequences: Option<&DownstreamSequences>,
    inputs: &[&str],
) -> bool {
    row_taken(
        klass,
        name,
        reference.len() - 1,
        reference,
        candidate,
        sequences,
        inputs,
    )
}

fn row_taken(
    klass: &str,
    name: &str,
    at: usize,
    reference: &[String],
    candidate: &[String],
    sequences: Option<&DownstreamSequences>,
    inputs: &[&str],
) -> bool {
    let (expected, actual) = (parsed(reference), parsed(candidate));
    let records = RecordRows::new(&expected, inputs);
    let record_rows = records.of(&expected[at]).expect("no refusal");
    asks(
        klass,
        &expected[at],
        &actual[at],
        name,
        sequences,
        record_rows.as_deref(),
        &[],
    )
}

const LATER_ALLELE_RECORD: &str = "21\t1004\tv\tA\tAA,AAA\t.\t.\t.";

#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_the_records_own_pair_on_a_later_alleles_row() {
    let (_cache, sequences) = sequences_of("later-allele");
    // The first allele's rows print one pair; the later allele's reference row prints the
    // record's span moved by the offset and the candidate's the record's own.
    let rows = |reference_later: (&str, &str, Option<i64>),
                candidate_later: (&str, &str, Option<i64>)| {
        let first =
            |downstream| later_allele_row("A", "4-5", "2", "K/KX", Some(downstream), Some(3));
        let later = |(cds, downstream, offset): (&str, &str, Option<i64>), protein, amino_acids| {
            later_allele_row("AA", cds, protein, amino_acids, Some(downstream), offset)
        };
        (
            vec![first("NSCCPWLRAQIA"), later(reference_later, "3", "T/KX")],
            vec![first("KNSCCPWLRAQIA"), later(candidate_later, "2", "K/KX")],
        )
    };
    let (moved, own) = (
        ("7-8", "KLLLSVVKGPNS", Some(3)),
        ("4-5", "KKLLLSVVKGPNS", Some(3)),
    );
    let taken = |reference_later, candidate_later, sequences, inputs: &[&str]| {
        let (reference, candidate) = rows(reference_later, candidate_later);
        [0, 1].map(|at| {
            row_taken(
                SHIFTED,
                "DownstreamProtein",
                at,
                &reference,
                &candidate,
                sequences,
                inputs,
            )
        })
    };
    let record = [LATER_ALLELE_RECORD];
    assert_eq!(taken(moved, own, Some(&sequences), &record), [true, true]);
    assert_eq!(
        taken(moved, own, Some(&sequences), &[]),
        [true, false],
        "the record's alleles in no entry"
    );
    assert_eq!(
        taken(moved, own, None, &record),
        [false, false],
        "without the cache"
    );
    for (label, reference_later, candidate_later) in [
        (
            "a pair moved by another size",
            ("6-7", "KLLLSVVKGPNS", Some(3)),
            own,
        ),
        (
            "a candidate pair that is not the record's own",
            moved,
            ("5-6", "KKLLLSVVKGPNS", Some(3)),
        ),
        (
            "a value that is not the rule's reading",
            ("7-8", "LLLSVVKGPNS", Some(3)),
            own,
        ),
        (
            "a value planted on the candidate's side",
            moved,
            ("4-5", "KKLLLSVVKGPN", Some(3)),
        ),
    ] {
        assert_eq!(
            taken(reference_later, candidate_later, Some(&sequences), &record),
            [true, false],
            "{label}"
        );
    }
    // A later allele's row without an offset reads the pair the first allele's notation
    // left, with the allele as written, which this record's reference value is.
    let no_offset =
        |(cds, downstream, _): (&'static str, &'static str, Option<i64>)| (cds, downstream, None);
    assert_eq!(
        taken(no_offset(moved), no_offset(own), Some(&sequences), &record),
        [true, true]
    );
}

/// The DownstreamProtein of an insertion of `allele` after coding base 4 of the toy
/// forward transcript, read at its span moved by `moved_by` bases with the allele turned as
/// a shift of `turned_by` turns it: Ensembl VEP's rule, or the candidate's.
fn earlier_shift_reading(
    sequences: &DownstreamSequences,
    allele: &str,
    moved_by: i128,
    turned_by: i128,
    ensembl: bool,
) -> String {
    let e = parsed(&[later_allele_row(allele, "4-5", "2", "K/KX", None, Some(3))]);
    let transcript = sequences
        .transcript(&e[0])
        .expect("shards that read")
        .expect("the transcript");
    let mut pair = (5, 4);
    if moved_by != 0 {
        let model = transcript.model.as_ref().expect("a model");
        let readings = shifted_pair_readings(&e[0], model, moved_by, &[(Some(5), Some(4))]);
        assert_eq!(readings.len(), 1);
        pair = readings[0].1.expect("a pair inside the coding sequence");
    }
    transcript
        .downstream_value(
            pair.0,
            pair.1,
            &hgvs_rotated_allele(allele, turned_by, 1),
            false,
            ensembl,
        )
        .expect("a value")
}

/// The fabricated grid of a later allele whose printed CDS pair an earlier allele's
/// notation left: per cell its name, the record's alternates in order, each with the
/// HGVS_OFFSET both engines print on its row, the pair the reference prints on the last
/// allele's row, the reading its DownstreamProtein is (the span moved by so many bases, the
/// allele turned by so many), and whether the class takes the value. The record inserts
/// after coding base 4 of the toy forward transcript; the last alternate is the row asked.
type EarlierShiftCell = (
    &'static str,
    &'static [(&'static str, Option<i64>)],
    &'static str,
    (i128, i128),
    bool,
);
const EARLIER_SHIFT_GRID: [EarlierShiftCell; 12] = [
    (
        "the pair an earlier allele left, the value at the row's own offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "7-8",
        (2, 2),
        true,
    ),
    (
        "the pair moved by the row's own offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "6-7",
        (2, 2),
        true,
    ),
    (
        "an allele without an offset between the two",
        &[("A", Some(3)), ("AG", None), ("AC", Some(2))],
        "7-8",
        (2, 2),
        true,
    ),
    (
        "the pair of the last earlier allele that has an offset",
        &[("A", Some(1)), ("AG", Some(3)), ("AC", Some(2))],
        "7-8",
        (2, 2),
        true,
    ),
    (
        "the pair of the first allele where a later one has another",
        &[("A", Some(1)), ("AG", Some(3)), ("AC", Some(2))],
        "5-6",
        (2, 2),
        false,
    ),
    (
        "a pair moved by neither offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "5-6",
        (2, 2),
        false,
    ),
    (
        "the value read at the printed pair",
        &[("A", Some(3)), ("AC", Some(2))],
        "7-8",
        (3, 0),
        false,
    ),
    (
        "the first allele's row",
        &[("AC", Some(2))],
        "7-8",
        (2, 2),
        false,
    ),
    (
        "no offset on the row, the pair an earlier allele left, the allele as written",
        &[("A", Some(3)), ("AC", None)],
        "7-8",
        (3, 0),
        true,
    ),
    (
        "no offset on the row and none on an earlier allele's",
        &[("A", None), ("AC", None)],
        "7-8",
        (3, 0),
        false,
    ),
    (
        "no offset on the row, another value",
        &[("A", Some(3)), ("AC", None)],
        "7-8",
        (2, 2),
        false,
    ),
    (
        "no offset on the row, the record's own pair printed",
        &[("A", Some(3)), ("AC", None)],
        "4-5",
        (3, 0),
        false,
    ),
];

#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_the_pair_an_earlier_allele_left_and_nothing_wider() {
    let (_cache, sequences) = sequences_of("earlier-shift");
    let readings: BTreeSet<String> = [(0, 0), (2, 2), (3, 0)]
        .into_iter()
        .map(|(moved_by, turned_by)| {
            earlier_shift_reading(&sequences, "AC", moved_by, turned_by, true)
        })
        .collect();
    assert_eq!(
        readings.len(),
        3,
        "the readings the grid tells apart differ in value"
    );
    for (label, alternates, printed, (moved_by, turned_by), taken) in EARLIER_SHIFT_GRID {
        let (allele, offset) = alternates[alternates.len() - 1];
        let earlier: Vec<String> = alternates[..alternates.len() - 1]
            .iter()
            .map(|(other, shift)| later_allele_row(other, "4-5", "2", "K/KX", Some("KNS"), *shift))
            .collect();
        let asked = |cds: &str, protein: &str, amino_acids: &str, downstream: String| {
            let mut rows = earlier.clone();
            rows.push(later_allele_row(
                allele,
                cds,
                protein,
                amino_acids,
                Some(&downstream),
                offset,
            ));
            rows
        };
        let reference = asked(
            printed,
            "3",
            "T/KX",
            earlier_shift_reading(&sequences, allele, moved_by, turned_by, true),
        );
        let candidate = asked(
            "4-5",
            "2",
            "K/KX",
            earlier_shift_reading(&sequences, allele, 0, 0, false),
        );
        let spelt: Vec<String> = alternates
            .iter()
            .map(|(other, _)| format!("A{other}"))
            .collect();
        let record = format!("21\t1004\tv\tA\t{}\t.\t.\t.", spelt.join(","));
        let ask = |sequences, inputs: &[&str]| {
            last_row_taken(
                SHIFTED,
                "DownstreamProtein",
                &reference,
                &candidate,
                sequences,
                inputs,
            )
        };
        assert_eq!(ask(Some(&sequences), &[&record]), taken, "{label}");
        assert!(
            !ask(Some(&sequences), &[]),
            "{label}: without the record's alleles in order"
        );
        assert!(!ask(None, &[&record]), "{label}: without the cache");
    }
}

/// A row of a made-up record on the toy forward transcript that inserts after coding base
/// 6, the last base of codon 2.
fn blosum62_row(
    allele: &str,
    (cds, protein, amino_acids, blosum): (&str, &str, &str, Option<&str>),
    offset: Option<i64>,
) -> String {
    let mut extra = "IMPACT=MODERATE;STRAND=1".to_string();
    if let Some(offset) = offset {
        extra += &format!(";HGVS_OFFSET={offset}");
    }
    if let Some(blosum) = blosum {
        extra += &format!(";BLOSUM62={blosum}");
    }
    default_row(
        "b",
        "21:1006-1007",
        allele,
        TOY_FORWARD,
        "inframe_insertion",
        ["-", cds, protein, amino_acids, "-"],
        &extra,
    )
}

type Blosum62Cells = (
    &'static str,
    &'static str,
    &'static str,
    Option<&'static str>,
);
const BLOSUM62_OWN: Blosum62Cells = ("6-7", "2-3", "-/K", None);
const BLOSUM62_MOVED: Blosum62Cells = ("7-8", "3", "T/K", Some("-1"));
/// The fabricated grid of the class: per cell its name, the record's alternates in order,
/// each with the HGVS_OFFSET on its row, the cells on the reference's side (CDS position,
/// protein position, Amino_acids, BLOSUM62) and on the candidate's, and whether the class
/// takes the value. The last alternate is the row asked: three bases between two codons,
/// whose own pair has no reference residue.
type Blosum62Cell = (
    &'static str,
    &'static [(&'static str, Option<i64>)],
    Blosum62Cells,
    Blosum62Cells,
    bool,
);
const BLOSUM62_GRID: [Blosum62Cell; 10] = [
    (
        "a later allele read one base on, the first allele's offset",
        &[("C", Some(1)), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        true,
    ),
    (
        "a third allele, the first allele's offset",
        &[("C", Some(1)), ("G", None), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        true,
    ),
    (
        "the first allele's row",
        &[("AAA", Some(1))],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        false,
    ),
    (
        "a first allele without an offset",
        &[("C", None), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        false,
    ),
    (
        "a value that is not the matrix cell of its pair",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/K", Some("5")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "a reference residue that is not the moved codon's",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "S/K", Some("0")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "alternate residues that differ",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/R", Some("-1")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "a candidate pair the plugin scores",
        &[("C", Some(1)), ("AAA", None)],
        BLOSUM62_MOVED,
        ("6-7", "2-3", "N/K", None),
        false,
    ),
    (
        "a pair moved by another size",
        &[("C", Some(1)), ("AAA", None)],
        ("8-9", "3", "T/K", Some("-1")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "the value on the candidate's side alone",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/K", None),
        ("6-7", "2-3", "-/K", Some("-1")),
        false,
    ),
];

#[test]
fn blosum62_hgvs_shift_shared_across_alleles_takes_the_score_read_at_the_moved_span_and_nothing_wider(
) {
    let (_cache, sequences) = sequences_of("blosum62");
    assert_eq!(
        (blosum62('T', 'K'), blosum62('S', 'K'), blosum62('T', 'R')),
        (Some(-1), Some(0), Some(-1))
    );
    for (label, alternates, reference_cells, candidate_cells, taken) in BLOSUM62_GRID {
        let (allele, offset) = alternates[alternates.len() - 1];
        let earlier: Vec<String> = alternates[..alternates.len() - 1]
            .iter()
            .map(|(other, shift)| blosum62_row(other, ("6-7", "2-3", "-/X", None), *shift))
            .collect();
        let asked = |cells| {
            let mut rows = earlier.clone();
            rows.push(blosum62_row(allele, cells, offset));
            rows
        };
        let (reference, candidate) = (asked(reference_cells), asked(candidate_cells));
        let spelt: Vec<String> = alternates
            .iter()
            .map(|(other, _)| format!("A{other}"))
            .collect();
        let record = format!("21\t1006\tb\tA\t{}\t.\t.\t.", spelt.join(","));
        let ask = |sequences, inputs: &[&str]| {
            last_row_taken(
                BLOSUM62_SHARED_SHIFT,
                "BLOSUM62",
                &reference,
                &candidate,
                sequences,
                inputs,
            )
        };
        assert_eq!(ask(Some(&sequences), &[&record]), taken, "{label}");
        assert!(!ask(None, &[&record]), "{label}: without the cache");
        assert!(
            !ask(Some(&sequences), &[]),
            "{label}: without the record's alleles"
        );
    }
}

/// The row of 21 30166352 G>GAG on a reverse-strand transcript under --hgvs.
fn rotation_row(
    protein: &str,
    strand: &str,
    offset: &str,
    notation: &str,
    consequence: &str,
    change: &str,
) -> Row {
    let extra = format!("IMPACT=HIGH;STRAND={strand};HGVSc=ENST00000286808:{notation};HGVS_OFFSET={offset};DownstreamProtein={protein};ProteinLengthChange={change}");
    let row = default_row(
        "cldn17_ins_AG",
        "21:30166352-30166353",
        "AG",
        "ENST00000286808",
        consequence,
        ["-", "-", "-", "-", "-"],
        &extra,
    );
    Row::clone(&parsed(&[row])[0])
}

#[test]
fn downstream_reverse_strand_rotation_accepts_the_protein_of_an_allele_shifted_past_its_length() {
    let row = |protein: &str, offset: &str, notation: &str, change: &str| {
        rotation_row(
            protein,
            "-1",
            offset,
            notation,
            "frameshift_variant",
            change,
        )
    };
    let taken = |e: &Row, a: &Row, name: &str| asks(ROTATION, e, a, name, None, None, &[]);
    let (e, a) = (
        row("SS", "-7", "c.271_272dup", "-132"),
        row("SP", "-7", "c.271_272dup", "-132"),
    );
    assert!(taken(&e, &a, "DownstreamProtein"));
    let other_length = row("SP", "-7", "c.271_272dup", "-131");
    assert!(taken(&e, &other_length, "DownstreamProtein"));
    assert!(taken(&e, &other_length, "ProteinLengthChange"));
    let deletion = |protein: &str| row(protein, "-4", "c.100_102del", "-132");
    assert!(taken(&deletion("SS"), &deletion("SP"), "DownstreamProtein"));
}

#[test]
fn downstream_reverse_strand_rotation_refuses_a_shift_the_rotation_reads_alike() {
    let taken = |strand: &str, offset: &str, notation: &str, consequence: &str| {
        let row =
            |protein: &str| rotation_row(protein, strand, offset, notation, consequence, "-132");
        asks(
            ROTATION,
            &row("SS"),
            &row("SP"),
            "DownstreamProtein",
            None,
            None,
            &[],
        )
    };
    let frameshift = "frameshift_variant";
    assert!(taken("-1", "-7", "c.271_272dup", frameshift));
    for (strand, offset, notation, consequence) in [
        ("1", "-7", "c.271_272dup", frameshift),
        ("-1", "-2", "c.271_272dup", frameshift),
        ("-1", "-1", "c.271_272dup", frameshift),
        ("-1", "-6", "c.271_272dup", frameshift),
        ("-1", "-", "c.271_272dup", frameshift),
        ("-1", "-7", "c.271A>G", frameshift),
        ("-1", "-7", "c.10_5del", frameshift),
        ("-1", "-7", "c.271_272dup", "inframe_insertion"),
    ] {
        assert!(
            !taken(strand, offset, notation, consequence),
            "{strand} {offset} {notation} {consequence}"
        );
    }
    // A value one engine lacks stays a difference: the class is two values that differ.
    let row = |protein: &str| rotation_row(protein, "-1", "-7", "c.271_272dup", frameshift, "-132");
    assert!(!asks(
        ROTATION,
        &row("SS"),
        &row("-"),
        "DownstreamProtein",
        None,
        None,
        &[]
    ));
    assert!(!asks(
        ROTATION,
        &row("-"),
        &row("SP"),
        "DownstreamProtein",
        None,
        None,
        &[]
    ));
}

/// The rows of a multi-base delins on one transcript: the reference's with the Downstream
/// plugin's two fields, the candidate's with `candidate_fields`.
fn swap_rows(reference_terms: &str, candidate_terms: &str, candidate_fields: &str) -> (Row, Row) {
    let row = |terms: &str, fields: &str| {
        let text = default_row(
            "delins_swap",
            "12:5000103-5000106",
            "ACACACAC",
            "ENST00000000051",
            terms,
            ["-", "-", "-", "-", "-"],
            &format!("IMPACT=HIGH;STRAND=-1{fields}"),
        );
        Row::clone(&parsed(&[text])[0])
    };
    (
        row(
            reference_terms,
            ";DownstreamProtein=KL;ProteinLengthChange=-120",
        ),
        row(candidate_terms, candidate_fields),
    )
}

/// One question to a class with the two values as given.
fn asked_with(
    klass: &str,
    e: &Row,
    name: &str,
    expected: &str,
    actual: &str,
    a: Option<&Row>,
) -> bool {
    BuiltIn
        .accepts(&Consult {
            klass,
            e,
            name,
            expected,
            actual,
            a,
            record_rows: None,
            candidate_rows: None,
            declared: &BTreeSet::new(),
            assembly: None,
            source: None,
            gwas_pairs: None,
            revel_table: None,
            utrannotator_sweep: None,
            downstream_sequences: None,
            splice_windows: None,
            spliceai_table: None,
            exon_spans: None,
        })
        .expect("no refusal")
}

#[test]
fn splice_family_swap_accepts_the_downstream_fields_the_reference_alone_writes() {
    let (e, a) = swap_rows(
        "frameshift_variant",
        "frameshift_variant,splice_region_variant",
        "",
    );
    for (field, value) in [("DownstreamProtein", "KL"), ("ProteinLengthChange", "-120")] {
        assert!(
            asked_with(SWAP, &e, field, value, ABSENT, Some(&a)),
            "{field}"
        );
    }
}

#[test]
fn splice_family_swap_refuses_a_downstream_field_on_any_other_pair_of_rows() {
    for candidate_terms in [
        "frameshift_variant",
        "frameshift_variant,splice_donor_variant",
    ] {
        let (e, a) = swap_rows("frameshift_variant", candidate_terms, "");
        assert!(
            !asked_with(SWAP, &e, "DownstreamProtein", "KL", ABSENT, Some(&a)),
            "{candidate_terms}"
        );
    }
    // The plugin writes on a frameshift row none of whose terms names a splice site: a
    // reference row without the frameshift, or with a splice term of its own, wrote nothing.
    for reference_terms in [
        "intron_variant",
        "frameshift_variant,splice_donor_5th_base_variant",
    ] {
        let candidate_terms = format!("{reference_terms},splice_region_variant");
        let (e, a) = swap_rows(reference_terms, &candidate_terms, "");
        assert!(
            !asked_with(SWAP, &e, "DownstreamProtein", "KL", ABSENT, Some(&a)),
            "{reference_terms}"
        );
    }
    let (e, a) = swap_rows(
        "frameshift_variant",
        "frameshift_variant,splice_region_variant",
        "",
    );
    assert!(!asked_with(
        SWAP,
        &e,
        "DownstreamProtein",
        "KL",
        "KR",
        Some(&a)
    ));
    assert!(!asked_with(
        SWAP,
        &a,
        "DownstreamProtein",
        ABSENT,
        "KL",
        Some(&e)
    ));
    assert!(!asked_with(SWAP, &e, "REVEL", "0.5", ABSENT, Some(&a)));
    assert!(!asked_with(
        SWAP,
        &e,
        "DownstreamProtein",
        "KL",
        ABSENT,
        None
    ));
}

/// A deletion's Transcript row with these terms and this Extra text.
fn following_entry(consequence: &str, extra: &str) -> Row {
    let text = format!("1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{consequence}\t7055-7059\t6978-?\t-\t-\t-\t-\t{extra}");
    Row::clone(&parsed(&[text])[0])
}

/// `term_pair_following_cell` on the field `name` of two rows of one key, a value of `None`
/// being a field the row leaves out.
fn following(
    name: &str,
    reference_terms: &str,
    reference_value: Option<&str>,
    candidate_terms: &str,
    candidate_value: Option<&str>,
) -> bool {
    let entry = |terms: &str, value: Option<&str>| match value {
        None => following_entry(terms, "STRAND=-1"),
        Some(value) => following_entry(terms, &format!("{name}={value};STRAND=-1")),
    };
    let (e, a) = (
        entry(reference_terms, reference_value),
        entry(candidate_terms, candidate_value),
    );
    term_pair_following_cell(
        &e,
        &a,
        name,
        e.fields.get_or(name, ABSENT),
        a.fields.get_or(name, ABSENT),
    )
}

const RETAINED: &str = "stop_retained_variant,3_prime_UTR_variant";
const UNKNOWN: &str = "coding_sequence_variant,3_prime_UTR_variant";
const FRAMESHIFT: &str = "frameshift_variant,3_prime_UTR_variant";
const STOP_GAINED: &str = "stop_gained,frameshift_variant";

fn spliced(term: &str) -> String {
    format!("splice_acceptor_variant,{term},intron_variant")
}

#[test]
fn term_impact_names_every_ranked_term() {
    let impacts: BTreeSet<&str> = TERMS_BY_RANK.iter().map(|(_, impact)| *impact).collect();
    assert_eq!(
        impacts,
        BTreeSet::from(["HIGH", "MODERATE", "LOW", "MODIFIER"])
    );
    let terms: BTreeSet<&str> = TERMS_BY_RANK.iter().map(|(term, _)| *term).collect();
    assert_eq!(terms.len(), 41, "each term once");
    let impact = |terms: &[&str]| terms_impact(&terms.iter().copied().collect());
    assert_eq!(
        impact(&["3_prime_UTR_variant", "stop_retained_variant"]),
        Some("LOW")
    );
    assert_eq!(
        impact(&[
            "NMD_transcript_variant",
            "frameshift_variant",
            "stop_gained"
        ]),
        Some("HIGH")
    );
    assert_eq!(impact(&["stop_retained_variant", "not_a_term"]), None);
    assert_eq!(impact(&[]), None);
}

#[test]
fn term_pair_following_cell_takes_the_impact_of_each_sides_terms() {
    assert!(following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        UNKNOWN,
        Some("MODIFIER")
    ));
    assert!(following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        FRAMESHIFT,
        Some("HIGH")
    ));
    assert!(following(
        "IMPACT",
        "inframe_deletion,stop_retained_variant",
        Some("MODERATE"),
        "frameshift_variant",
        Some("HIGH")
    ));
    assert!(following(
        "IMPACT",
        UNKNOWN,
        Some("MODIFIER"),
        RETAINED,
        Some("LOW")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        UNKNOWN,
        Some("HIGH")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("MODERATE"),
        UNKNOWN,
        Some("MODIFIER")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        RETAINED,
        Some("MODIFIER")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        Some("LOW"),
        "not_a_term,3_prime_UTR_variant",
        Some("MODIFIER")
    ));
    assert!(!following(
        "IMPACT",
        RETAINED,
        None,
        UNKNOWN,
        Some("MODIFIER")
    ));
}

#[test]
fn term_pair_following_cell_takes_downstream_protein_where_frameshift_comes() {
    let name = "DownstreamProtein";
    assert!(following(name, RETAINED, None, FRAMESHIFT, Some("MKV")));
    assert!(following(name, FRAMESHIFT, Some("MKV"), RETAINED, None));
    assert!(!following(name, RETAINED, None, UNKNOWN, Some("MKV")));
    assert!(!following(name, RETAINED, Some("MKV"), FRAMESHIFT, None));
    assert!(!following(
        name,
        RETAINED,
        Some("MKL"),
        FRAMESHIFT,
        Some("MKV")
    ));
    assert!(!following(
        name,
        &spliced("stop_retained_variant"),
        None,
        &spliced("frameshift_variant"),
        Some("MKV")
    ));
    assert!(!following(
        name,
        FRAMESHIFT,
        Some("MKL"),
        FRAMESHIFT,
        Some("MKV")
    ));
}

#[test]
fn term_pair_following_cell_takes_protein_length_change_where_frameshift_comes() {
    let name = "ProteinLengthChange";
    assert!(following(name, RETAINED, None, FRAMESHIFT, Some("-12")));
    assert!(following(name, FRAMESHIFT, Some("-12"), RETAINED, None));
    assert!(!following(name, RETAINED, None, UNKNOWN, Some("-12")));
    assert!(!following(
        name,
        RETAINED,
        Some("4"),
        FRAMESHIFT,
        Some("-12")
    ));
    assert!(!following(
        name,
        &spliced("stop_retained_variant"),
        None,
        &spliced("frameshift_variant"),
        Some("-12")
    ));
}

#[test]
fn term_pair_following_cell_takes_nmd_where_a_truncating_term_comes() {
    assert!(following(
        "NMD",
        RETAINED,
        None,
        FRAMESHIFT,
        Some(NMD_VALUE)
    ));
    assert!(following(
        "NMD",
        RETAINED,
        None,
        STOP_GAINED,
        Some(NMD_VALUE)
    ));
    assert!(following(
        "NMD",
        FRAMESHIFT,
        Some(NMD_VALUE),
        RETAINED,
        None
    ));
    let (retained, coding) = (
        spliced("stop_retained_variant"),
        spliced("coding_sequence_variant"),
    );
    assert!(!following("NMD", &retained, Some(NMD_VALUE), &coding, None));
    assert!(!following("NMD", &retained, None, &coding, Some(NMD_VALUE)));
    assert!(!following("NMD", RETAINED, None, UNKNOWN, Some(NMD_VALUE)));
    assert!(!following(
        "NMD",
        RETAINED,
        Some(NMD_VALUE),
        FRAMESHIFT,
        None
    ));
}

#[test]
fn term_pair_following_cell_takes_the_loftee_fields_where_a_truncating_term_comes() {
    assert!(following("LoF", RETAINED, None, FRAMESHIFT, Some("HC")));
    assert!(following("LoF", FRAMESHIFT, Some("LC"), RETAINED, None));
    assert!(!following(
        "LoF",
        FRAMESHIFT,
        Some("HC"),
        STOP_GAINED,
        Some("LC")
    ));
    assert!(!following("LoF", RETAINED, None, UNKNOWN, Some("HC")));
    assert!(!following("LoF", RETAINED, Some("HC"), FRAMESHIFT, None));
    let info = "PERCENTILE:0.97";
    assert!(following(
        "LoF_info",
        RETAINED,
        None,
        FRAMESHIFT,
        Some("PERCENTILE:0.97,GERP_DIST:12.5")
    ));
    assert!(!following(
        "LoF_info",
        FRAMESHIFT,
        Some(info),
        STOP_GAINED,
        Some("PERCENTILE:0.96")
    ));
    assert!(!following("LoF_info", RETAINED, None, UNKNOWN, Some(info)));
    assert!(!following("LoF_info", RETAINED, Some(info), UNKNOWN, None));
    assert!(following(
        "LoF_filter",
        RETAINED,
        None,
        FRAMESHIFT,
        Some("END_TRUNC")
    ));
    assert!(!following(
        "LoF_filter",
        FRAMESHIFT,
        Some("END_TRUNC"),
        STOP_GAINED,
        Some("END_TRUNC,GC_TO_GT_DONOR")
    ));
    assert!(!following(
        "LoF_filter",
        RETAINED,
        None,
        UNKNOWN,
        Some("END_TRUNC")
    ));
    let flags = "SINGLE_EXON,PHYLOCSF_WEAK";
    assert!(following(
        "LoF_flags",
        RETAINED,
        None,
        FRAMESHIFT,
        Some("SINGLE_EXON")
    ));
    assert!(following(
        "LoF_flags",
        STOP_GAINED,
        Some(flags),
        RETAINED,
        None
    ));
    assert!(!following(
        "LoF_flags",
        RETAINED,
        Some("SINGLE_EXON"),
        FRAMESHIFT,
        None
    ));
    assert!(!following(
        "LoF_flags",
        RETAINED,
        None,
        UNKNOWN,
        Some("SINGLE_EXON")
    ));
    assert!(!following(
        "LoF_flags",
        FRAMESHIFT,
        Some("SINGLE_EXON"),
        STOP_GAINED,
        Some(flags)
    ));
}

#[test]
fn term_pair_following_cell_takes_a_frameshift_notation_beside_a_lost_retained_stop() {
    let protein = "ENSP00000490161.1:p.";
    let hgvsp = |reference_terms: &str, reference: &str, candidate_terms: &str, candidate: &str| {
        let (reference, candidate) = (
            format!("{protein}{reference}"),
            format!("{protein}{candidate}"),
        );
        following(
            "HGVSp",
            reference_terms,
            Some(&reference),
            candidate_terms,
            Some(&candidate),
        )
    };
    assert!(hgvsp(
        RETAINED,
        "Tyr22_Ter23delinsCys",
        UNKNOWN,
        "Tyr22CysfsTer4"
    ));
    assert!(hgvsp(RETAINED, "Ter23=", UNKNOWN, "Ter23TyrextTer4"));
    assert!(hgvsp(RETAINED, "Ter23%3D", FRAMESHIFT, "Ter23TyrextTer?"));
    assert!(hgvsp(
        UNKNOWN,
        "Tyr22CysfsTer4",
        RETAINED,
        "Tyr22_Ter23delinsCys"
    ));
    assert!(!hgvsp(
        RETAINED,
        "Tyr22_Ter23delinsCys",
        UNKNOWN,
        "Tyr22_Ter23delinsSer"
    ));
    assert!(!hgvsp(
        RETAINED,
        "Tyr22CysfsTer4",
        UNKNOWN,
        "Tyr22CysfsTer9"
    ));
    assert!(!following(
        "HGVSp",
        RETAINED,
        Some(&format!("{protein}Tyr22_Ter23delinsCys")),
        UNKNOWN,
        Some("ENSP00000490162.1:p.Tyr22CysfsTer4")
    ));
    assert!(!following(
        "HGVSp",
        RETAINED,
        None,
        UNKNOWN,
        Some(&format!("{protein}Tyr22CysfsTer4"))
    ));
    assert!(!hgvsp(
        UNKNOWN,
        "Tyr22_Ter23delinsCys",
        FRAMESHIFT,
        "Tyr22CysfsTer4"
    ));
    assert!(!following(
        "HGVSc",
        RETAINED,
        Some("ENST00000635895.1:c.6978_*2del"),
        UNKNOWN,
        Some("ENST00000635895.1:c.6979_*3del")
    ));
}

#[test]
fn plugin_value_follows_stop_codon_term_pair_takes_a_field_one_side_alone_writes() {
    // The deletion's row lacks a CDS end, so the stop codon test is made on neither side.
    let row = |terms: &str, extra: &str| {
        let text = format!("1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{terms}\t?-7057\t-\t-\t-\t-\t-\t{extra}");
        Row::clone(&parsed(&[text])[0])
    };
    let retained = row(
        "splice_acceptor_variant,stop_retained_variant",
        "STRAND=-1;NMD=NMD_escaping_variant",
    );
    let shifted = row(
        "splice_acceptor_variant,frameshift_variant",
        "STRAND=-1;NMD=NMD_escaping_variant;DownstreamProtein=MKV",
    );
    assert!(
        !asks(
            FOLLOWS_STOP,
            &retained,
            &shifted,
            "DownstreamProtein",
            None,
            None,
            &[]
        ),
        "a splice term keeps the Downstream plugin out on both sides"
    );
    let plain = |terms: &str, extra: &str| {
        let text = format!("1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{terms}\t?-7057\t-\t-\t-\t-\t-\t{extra}");
        Row::clone(&parsed(&[text])[0])
    };
    let kept = plain("stop_retained_variant,3_prime_UTR_variant", "STRAND=-1");
    let frameshift = plain(
        "frameshift_variant,3_prime_UTR_variant",
        "STRAND=-1;DownstreamProtein=MKV;NMD=NMD_escaping_variant;LoF=HC",
    );
    for name in ["DownstreamProtein", "NMD", "LoF"] {
        assert!(
            asks(FOLLOWS_STOP, &kept, &frameshift, name, None, None, &[]),
            "{name}"
        );
        assert!(
            !asks(FOLLOWS_STOP, &frameshift, &kept, name, None, None, &[]),
            "{name}: the other way round"
        );
        assert!(
            !asks(FOLLOWS_SPLICE, &kept, &frameshift, name, None, None, &[]),
            "{name}: the splice class"
        );
    }
    assert!(
        !asks(FOLLOWS_STOP, &kept, &frameshift, "IMPACT", None, None, &[]),
        "a field the class does not read"
    );
    let stop_lost = plain(
        "stop_lost,3_prime_UTR_variant",
        "STRAND=-1;NMD=NMD_escaping_variant",
    );
    assert!(
        !asks(FOLLOWS_STOP, &kept, &stop_lost, "NMD", None, None, &[]),
        "a pair of another stop codon class"
    );
}

#[test]
fn plugin_value_follows_splice_term_past_the_record_span_takes_a_field_one_side_alone_writes() {
    // Three positions inside an exon replaced by five letters, under an identifier.
    let row = |name: &str, location: &str, allele: &str, terms: &str, cdna: &str, extra: &str| {
        let text = format!("{name}\t{location}\t{allele}\tENSG1\tENST00000635895\tTranscript\t{terms}\t{cdna}\t-\t-\t-\t-\t-\t{extra}");
        Row::clone(&parsed(&[text])[0])
    };
    let three = |terms: &str, extra: &str| {
        row("rec1", "7:100998-101000", "GTGTG", terms, "998-1000", extra)
    };
    let reference = three(
        "splice_donor_variant,protein_altering_variant",
        "STRAND=1;NMD=NMD_escaping_variant;LoF=HC",
    );
    for candidate_terms in [
        "protein_altering_variant,splice_region_variant",
        "protein_altering_variant",
    ] {
        let candidate = three(candidate_terms, "STRAND=1");
        for name in ["NMD", "LoF"] {
            assert!(
                asks(
                    FOLLOWS_SPLICE,
                    &reference,
                    &candidate,
                    name,
                    None,
                    None,
                    &[]
                ),
                "{name}"
            );
            assert!(
                !asks(FOLLOWS_STOP, &reference, &candidate, name, None, None, &[]),
                "{name}"
            );
        }
    }
    let downstream = three(
        "frameshift_variant",
        "STRAND=1;DownstreamProtein=MKV;ProteinLengthChange=-12",
    );
    let beside_a_donor = three("splice_donor_variant,frameshift_variant", "STRAND=1");
    for name in ["DownstreamProtein", "ProteinLengthChange"] {
        assert!(
            asks(
                FOLLOWS_SPLICE,
                &beside_a_donor,
                &downstream,
                name,
                None,
                None,
                &[]
            ),
            "{name}"
        );
    }
    // Two positions among the record's alleles, which the first column spells; under an
    // identifier the same two positions are an insertion's Location too.
    let two = |name: &str, terms: &str, extra: &str| {
        row(name, "7:100999-101000", "GTGTG", terms, "999-1000", extra)
    };
    let donor_side =
        "splice_donor_variant,protein_altering_variant,splice_donor_region_variant,intron_variant";
    let region_back = "protein_altering_variant,splice_region_variant";
    for (name, taken) in [("7_100999_TA/GTGTG", true), ("rec2", false)] {
        let reference = two(name, donor_side, "STRAND=1;NMD=NMD_escaping_variant");
        let candidate = two(name, region_back, "STRAND=1");
        assert_eq!(
            asks(
                FOLLOWS_SPLICE,
                &reference,
                &candidate,
                "NMD",
                None,
                None,
                &[]
            ),
            taken,
            "{name}"
        );
    }
    // The reference bases a caller reads from the record's input line stand in for the
    // first column's alleles.
    let reference = two("rec2", donor_side, "STRAND=1;NMD=NMD_escaping_variant");
    let candidate = two("rec2", region_back, "STRAND=1");
    *reference.reference.borrow_mut() = Some("TA".to_string());
    assert!(asks(
        FOLLOWS_SPLICE,
        &reference,
        &candidate,
        "NMD",
        None,
        None,
        &[]
    ));
}

#[test]
fn plugin_hgvs_shifted_cds_coordinates_takes_a_row_the_start_rule_changes_too_only_beside_its_class(
) {
    let (_cache, sequences) = sequences_of("two-causes");
    // On a run of four bases from a codon's third the deletion of base 6 beside an offset
    // of 3 is read at base 9: the reference translates from codon 4, its reading of the
    // record's own pair would start at codon 3, and the candidate starts at codon 2.
    let both_causes = frameshift("21:6006", "-", RUN, "6").offset(3);
    let (e, a) = (
        both_causes.row(Some("VK"), Some("1")),
        both_causes.row(Some("PKVK"), Some("1")),
    );
    let found = downstream_class(&e, &a, "DownstreamProtein", "VK", "PKVK", &sequences, None);
    assert_eq!(found.expect("no refusal"), Some((SHIFTED, Some(START))));
    let taken = |klass: &str, declared: &[&str]| {
        asks(
            klass,
            &e,
            &a,
            "DownstreamProtein",
            Some(&sequences),
            None,
            declared,
        )
    };
    assert!(taken(SHIFTED, &[SHIFTED, START]), "beside the start class");
    assert!(
        !taken(START, &[SHIFTED, START]),
        "the row is the shifted pair's class's"
    );
    assert!(!taken(SHIFTED, &[SHIFTED]), "alone");
    assert!(!taken(START, &[START]), "the start class alone");
    assert!(
        !taken(SHIFTED, &[SHIFTED, UTR_FRAME]),
        "beside the other Downstream class"
    );
    assert!(!taken(SHIFTED, &[]), "no class declared");
    // The own reading's value on the same row without the offset is the start class's alone.
    let own = frameshift("21:6006", "-", RUN, "6");
    assert!(protein_taken(
        START,
        &own,
        "KVK",
        "PKVK",
        "1",
        Some(&sequences)
    ));
    assert!(!protein_taken(
        SHIFTED,
        &own,
        "KVK",
        "PKVK",
        "1",
        Some(&sequences)
    ));
    // A row the flag alone changes is set aside by this class alone.
    let flag = frameshift("21:1005", "-", TOY_FORWARD, "5").offset(2);
    assert!(protein_taken(
        SHIFTED,
        &flag,
        "LLLSVVKGPNS",
        "KLLLSVVKGPNS",
        "5",
        Some(&sequences)
    ));
}

/// The coding cells of a row of the made-up record: its CDS position, its protein
/// position and its `Amino_acids`.
type CodingCells = (&'static str, &'static str, &'static str);

#[test]
fn the_shift_an_earlier_allele_left_takes_three_columns_of_a_later_alleles_row() {
    let (_cache, sequences) = sequences_of("shared-shift");
    let columns = ["CDS_position", "Protein_position", "Amino_acids"];
    let (moved, own): (CodingCells, CodingCells) = (("7-8", "3", "T/KX"), ("4-5", "2", "K/KX"));
    // The columns the rule takes on the later allele's row: the reference's first and later
    // rows, the candidate's later row and the offset on the first allele's rows as given.
    let taken = |reference_first: CodingCells,
                 reference_later: CodingCells,
                 candidate_later: CodingCells,
                 first_offset: Option<i64>,
                 inputs: &[&str]| {
        let row = |allele: &str, (cds, protein, amino_acids): CodingCells, offset| {
            later_allele_row(allele, cds, protein, amino_acids, None, offset)
        };
        let reference = parsed(&[
            row("A", reference_first, first_offset),
            row("AA", reference_later, Some(3)),
        ]);
        let candidate = parsed(&[
            row("A", own, first_offset),
            row("AA", candidate_later, Some(3)),
        ]);
        let records = RecordRows::new(&reference, inputs);
        let rows_of = |at: usize| records.of(&reference[at]).expect("no refusal");
        let asked = |at: usize, name: &str| {
            later_allele_shared_shift(
                &reference[at],
                &candidate[at],
                name,
                rows_of(at).as_deref(),
                &sequences,
            )
            .expect("no refusal")
        };
        assert!(
            !columns.iter().any(|name| asked(0, name)),
            "the first allele's row"
        );
        columns
            .into_iter()
            .filter(|name| asked(1, name))
            .collect::<Vec<_>>()
    };
    let record = [LATER_ALLELE_RECORD];
    assert_eq!(taken(own, moved, own, Some(3), &record), columns);
    assert!(
        taken(own, moved, own, Some(3), &[]).is_empty(),
        "without the record's alleles"
    );
    for (label, reference_first, reference_later, candidate_later, first_offset) in [
        (
            "the moved values on a first allele's row",
            moved,
            own,
            own,
            Some(3),
        ),
        ("a first allele without an offset", own, moved, own, None),
        (
            "a first allele that moved two bases",
            own,
            moved,
            own,
            Some(2),
        ),
        (
            "a move of another size",
            own,
            ("6-7", "3", "T/KX"),
            own,
            Some(3),
        ),
        (
            "a candidate that does not print the own pair",
            own,
            moved,
            ("5-6", "2", "K/KX"),
            Some(3),
        ),
    ] {
        assert!(
            taken(
                reference_first,
                reference_later,
                candidate_later,
                first_offset,
                &record
            )
            .is_empty(),
            "{label}"
        );
    }
    // A peptide with another alternate part, or whose reference part is not its codon's
    // residue, stays open alone; a protein position that is not the moved pair's codon
    // leaves the peptide open as well.
    for (label, reference_later, candidate_later, held) in [
        ("another alternate part", ("7-8", "3", "T/TX"), own, 2),
        (
            "a reference part that is not the moved codon's",
            ("7-8", "3", "L/KX"),
            own,
            2,
        ),
        (
            "a candidate reference part that is not the own codon's",
            moved,
            ("4-5", "2", "L/KX"),
            2,
        ),
        (
            "a protein position that is not the moved pair's codon",
            ("7-8", "4", "T/KX"),
            own,
            1,
        ),
    ] {
        assert_eq!(
            taken(own, reference_later, candidate_later, Some(3), &record),
            columns[..held],
            "{label}"
        );
    }
}

/// A later allele's row of an in-frame insertion prints, at the pair the first allele's
/// offset moved, the record's own allele inserted into the moved codon's window
/// (TranscriptVariationAllele.pm:841-843, :877), upper-cased from the own pair's codon
/// position (:911-923), and the own pair's reference codons. The rule takes the
/// Amino_acids and Codons each side writes so, at a codon boundary and inside a codon, on
/// either strand, as `tools/vep-formats-compare` takes them; a reference part that is not
/// the moved codon's residue and a window cased otherwise stay differences.
#[test]
fn the_shift_an_earlier_allele_left_takes_the_codons_of_an_insertion_read_at_the_moved_pair() {
    let (_cache, sequences) = sequences_of("shared-shift-codons");
    // (Location, Feature, strand, first and later alleles, record, the later allele's
    // reference and candidate cells (CDS, protein, Amino_acids, Codons), the first offset)
    type Cells<'a> = [&'a str; 4];
    type Case<'a> = (
        &'a str,
        &'a str,
        &'a str,
        [&'a str; 2],
        &'a str,
        Cells<'a>,
        Cells<'a>,
        i64,
    );
    let boundary_own: Cells = ["6-7", "2-3", "-/PP", "-/CCCCCC"];
    let cases: [Case; 3] = [
        (
            "21:1006-1007",
            TOY_FORWARD,
            "1",
            ["CCC", "CCCCCC"],
            "21\t1006\tw\tA\tACCC,ACCCCCC\t.\t.\t.",
            ["10-11", "4", "L/PPL", "-/CCCCCCctt"],
            boundary_own,
            4,
        ),
        (
            "21:1007-1008",
            TOY_FORWARD,
            "1",
            ["CCC", "CCCCCC"],
            "21\t1007\tw\tA\tACCC,ACCCCCC\t.\t.\t.",
            ["11-12", "4", "L/LPP", "act/cTCCCCCct"],
            ["7-8", "3", "T/TPP", "act/aCCCCCCct"],
            4,
        ),
        (
            "21:2036-2037",
            TOY_REVERSE,
            "-1",
            ["GGG", "GGGGGG"],
            "21\t2036\tw\tT\tTGGG,TGGGGGG\t.\t.\t.",
            ["10-11", "4", "L/PPL", "-/CCCCCCctt"],
            boundary_own,
            -4,
        ),
    ];
    for (location, feature, strand, alleles, record, reference, candidate, offset) in cases {
        let row = |allele: &str, [cds, protein, amino_acids, codons]: Cells| {
            default_row(
                "w",
                location,
                allele,
                feature,
                "inframe_insertion",
                ["-", cds, protein, amino_acids, codons],
                &format!("IMPACT=MODERATE;STRAND={strand};HGVS_OFFSET={offset}"),
            )
        };
        let first: Cells = [candidate[0], candidate[1], "-/X", "-/X"];
        let taken = |later: Cells| {
            let reference = parsed(&[row(alleles[0], first), row(alleles[1], later)]);
            let candidate = parsed(&[row(alleles[0], first), row(alleles[1], candidate)]);
            let records = RecordRows::new(&reference, &[record]);
            let rows = records.of(&reference[1]).expect("no refusal");
            ["Amino_acids", "Codons"]
                .into_iter()
                .filter(|name| {
                    later_allele_shared_shift(
                        &reference[1],
                        &candidate[1],
                        name,
                        rows.as_deref(),
                        &sequences,
                    )
                    .expect("no refusal")
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(taken(reference), ["Amino_acids", "Codons"], "{record}");
        let other_residue = format!("K{}", &reference[2][1..]);
        let shifted = [
            reference[0],
            reference[1],
            other_residue.as_str(),
            reference[3],
        ];
        assert_eq!(taken(shifted), ["Codons"], "{record}");
        let upper = reference[3].to_ascii_uppercase();
        let cased = [reference[0], reference[1], reference[2], upper.as_str()];
        assert_eq!(taken(cased), ["Amino_acids"], "{record}");
    }
    // A two-base replacement prints the pair an insertion between its bases prints, and is
    // no insertion: no window is read for it. The row name spells the record's alleles.
    let row = |allele: &str, [cds, protein, amino_acids, codons]: Cells| {
        default_row(
            "21_1006_AA/CCC/CCCCCC",
            "21:1006-1007",
            allele,
            TOY_FORWARD,
            "inframe_insertion",
            ["-", cds, protein, amino_acids, codons],
            "IMPACT=MODERATE;STRAND=1;HGVS_OFFSET=4",
        )
    };
    let first: Cells = [boundary_own[0], boundary_own[1], "-/X", "-/X"];
    let reference = parsed(&[
        row("CCC", first),
        row("CCCCCC", ["10-11", "4", "L/PPL", "-/CCCCCCctt"]),
    ]);
    let candidate = parsed(&[row("CCC", first), row("CCCCCC", boundary_own)]);
    let records = RecordRows::new(&reference, &["21\t1006\tw\tAA\tCCC,CCCCCC\t.\t.\t."]);
    let rows = records.of(&reference[1]).expect("no refusal");
    for name in ["Amino_acids", "Codons"] {
        let taken = later_allele_shared_shift(
            &reference[1],
            &candidate[1],
            name,
            rows.as_deref(),
            &sequences,
        )
        .expect("no refusal");
        assert!(!taken, "{name}");
    }
}
