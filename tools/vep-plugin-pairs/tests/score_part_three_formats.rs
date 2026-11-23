// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` program under the tool's own registry and classes over two fabricated grids
//! of records of several alternates on one toy transcript: the class of the shifted CDS
//! coordinates on the DownstreamProtein of a later allele, and the class of the BLOSUM62
//! score read at the span an earlier allele's shift moved. Each cell's two outputs, its
//! input of one record and the candidate cache are written by the test.

mod support;

use serde_json::{json, Value};
use support::{counts, entry, message, residual_lines, try_score, Run, Scratch};
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::rules::{
    BuiltIn, BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES, PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
};

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
// A coding sequence that holds a run of one base from a codon's third base.
const RUN_CODING: &str = "ATGCCAAAAGGTTAA";
const THREE_EXONS_UTR: &str = "GGGCCCAAAT";
const SEAM_PEPTIDE: &str = "MKTLAVRG";

/// The transcript the records of both grids stand on: one exon on the forward strand,
/// coding base c at 21:1000+c.
const TOY_FORWARD: &str = "ENST00000000201";

/// A cache transcript without a model: its sequences alone.
fn transcript(stable_id: &str, sequences: Value) -> Value {
    json!({"stable_id": stable_id, "variation_effect_feature_cache": sequences})
}

/// A cache transcript with a model: each of `pairs` a run of cDNA coordinates and the
/// genomic bases it maps to, the exons those genomic runs unless `exons` gives them, the
/// coding sequence from cDNA base 1.
struct Modelled<'a> {
    stable_id: &'a str,
    strand: i64,
    pairs: &'a [[i64; 4]],
    coding: &'a str,
    utr: &'a str,
    exons: &'a [(i64, i64)],
    /// Whether every number is written as a string, as a converted cache writes it.
    text: bool,
    peptide: Option<&'a str>,
}

impl<'a> Modelled<'a> {
    /// A transcript of the toy sequences.
    fn new(stable_id: &'a str, strand: i64, pairs: &'a [[i64; 4]]) -> Modelled<'a> {
        Modelled {
            stable_id,
            strand,
            pairs,
            coding: TOY_CODING,
            utr: TOY_UTR,
            exons: &[],
            text: false,
            peptide: None,
        }
    }

    fn written(&self) -> Value {
        let number = |value: i64| match self.text {
            true => json!(value.to_string()),
            false => json!(value),
        };
        let pairs: Vec<Value> = self
            .pairs
            .iter()
            .map(|[from_start, from_end, to_start, to_end]| {
                json!({
                    "from_start": number(*from_start),
                    "from_end": number(*from_end),
                    "to_start": number(*to_start),
                    "to_end": number(*to_end),
                    "ori": number(self.strand),
                })
            })
            .collect();
        let mut spans: Vec<(i64, i64)> = match self.exons {
            [] => self.pairs.iter().map(|pair| (pair[2], pair[3])).collect(),
            exons => exons.to_vec(),
        };
        spans.sort_unstable();
        if self.strand == -1 {
            spans.reverse();
        }
        let exons: Vec<Value> = spans
            .iter()
            .map(|(start, end)| json!({"start": number(*start), "end": number(*end)}))
            .collect();
        let mut sequences = json!({
            "translateable_seq": self.coding,
            "three_prime_utr": self.utr,
            "codon_table": 1,
            "mapper": {"pairs": pairs},
        });
        if let Some(peptide) = self.peptide {
            sequences["peptide"] = json!(peptide);
        }
        json!({
            "stable_id": self.stable_id,
            "strand": number(self.strand),
            "cdna_coding_start": number(1),
            "cdna_coding_end": number(i64::try_from(self.coding.len()).unwrap()),
            "exons": exons,
            "variation_effect_feature_cache": sequences,
        })
    }
}

/// A candidate JSON cache of twelve coding transcripts. Chromosome 21 holds the toy
/// transcript (its `codon_table` the string a converted cache writes), a copy without a 3'
/// UTR (the number a native cache writes), a copy under codon table 2, a copy whose 3' UTR
/// opens with a stop codon (no `codon_table`) and a copy under table 7, none with a model,
/// and six transcripts with one: one exon on the forward strand, one on the reverse strand
/// (its numbers the strings a converted cache writes), two exons, one exon whose mapper
/// pairs touch on the genome while the model carries three bases between them (with its
/// `peptide`), the run of one base, and a coding sequence of 120 bases over three exons.
/// Chromosome 22 holds the short transcript under table 0.
fn downstream_cache(scratch: &Scratch) -> String {
    let utr_stop = format!("TGA{}", &TOY_UTR[3..]);
    let three_exons_coding = format!("ATG{}TAA", "GCT".repeat(38));
    let chromosome_21 = Value::Array(vec![
        transcript(
            "ENST00000000101",
            json!({"translateable_seq": TOY_CODING, "three_prime_utr": TOY_UTR, "codon_table": "1"}),
        ),
        transcript(
            "ENST00000000102",
            json!({"translateable_seq": TOY_CODING, "codon_table": 1}),
        ),
        transcript(
            "ENST00000000103",
            json!({"translateable_seq": TOY_CODING, "three_prime_utr": TABLE_2_UTR, "codon_table": "2"}),
        ),
        transcript(
            "ENST00000000104",
            json!({"translateable_seq": TOY_CODING, "three_prime_utr": utr_stop}),
        ),
        transcript(
            "ENST00000000106",
            json!({"translateable_seq": TOY_CODING, "three_prime_utr": TOY_UTR, "codon_table": 7}),
        ),
        Modelled::new(TOY_FORWARD, 1, &[[1, 42, 1001, 1042]]).written(),
        Modelled {
            text: true,
            ..Modelled::new("ENST00000000202", -1, &[[1, 42, 2001, 2042]])
        }
        .written(),
        Modelled::new(
            "ENST00000000203",
            1,
            &[[1, 9, 3001, 3009], [10, 42, 3014, 3046]],
        )
        .written(),
        Modelled {
            exons: &[(4001, 4039)],
            peptide: Some(SEAM_PEPTIDE),
            ..Modelled::new(
                "ENST00000000204",
                1,
                &[[1, 9, 4001, 4009], [13, 42, 4010, 4039]],
            )
        }
        .written(),
        Modelled {
            coding: RUN_CODING,
            utr: SHORT_UTR,
            ..Modelled::new("ENST00000000205", 1, &[[1, 21, 6001, 6021]])
        }
        .written(),
        Modelled {
            coding: &three_exons_coding,
            utr: THREE_EXONS_UTR,
            ..Modelled::new(
                "ENST00000000206",
                1,
                &[
                    [1, 110, 5001, 5110],
                    [111, 120, 5201, 5210],
                    [121, 130, 5301, 5310],
                ],
            )
        }
        .written(),
    ]);
    scratch.file(
        "cache/transcripts/21/1-1000000.json",
        chromosome_21.to_string(),
    );
    let chromosome_22 = Value::Array(vec![transcript(
        "ENST00000000105",
        json!({"translateable_seq": SHORT_CODING, "three_prime_utr": SHORT_UTR, "codon_table": 0}),
    )]);
    scratch.file(
        "cache/transcripts/22/1-1000000.json",
        chromosome_22.to_string(),
    );
    scratch.path("cache")
}

/// A row of the toy forward transcript: its CDS position, protein position and amino
/// acids, and its Extra text.
fn coding_row(
    name: &str,
    location: &str,
    allele: &str,
    consequence: &str,
    coding: [&str; 3],
    extra: &str,
) -> String {
    let [cds, protein, amino_acids] = coding;
    format!(
        "{name}\t{location}\t{allele}\tENSG1\t{TOY_FORWARD}\tTranscript\t{consequence}\t-\t{cds}\t{protein}\t{amino_acids}\t-\t-\t{extra}\n"
    )
}

/// The alternates of a record after the reference base `A`, each with the HGVS_OFFSET both
/// engines print on its row.
type Alternates = &'static [(&'static str, Option<i64>)];

/// The ALT column of the record of `alternates`.
fn spelt(alternates: Alternates) -> String {
    let spelt: Vec<String> = alternates
        .iter()
        .map(|(other, _)| format!("A{other}"))
        .collect();
    spelt.join(",")
}

/// One grid scored cell by cell under one class of one plugin's fields.
struct Grid<'a> {
    scratch: &'a Scratch,
    cache: String,
    klass: &'a str,
    /// The class table's one row.
    table: String,
    plugins: &'a [&'a str],
    /// The source the grid's field is scored under.
    source: &'a str,
    /// The residual rows of a cell the class does not take.
    left: i64,
}

impl Grid<'_> {
    /// Scores one cell's outputs with the cache and the cell's record as the input, and
    /// holds the source's counts to `counted`: the class sets the value aside where the
    /// cell is `taken`, and the value stays in the residual file where it is not.
    fn holds(
        &self,
        number: usize,
        label: &str,
        rows: (&[String], &[String]),
        record: [&str; 5],
        taken: bool,
        counted: [i64; 6],
    ) {
        let name = |file: &str| format!("cell{number}/{file}");
        let reference = self.scratch.output(&name("reference.txt"), rows.0);
        let candidate = self.scratch.output(&name("candidate.txt"), rows.1);
        let input = self.scratch.vcf(&name("input.vcf"), &[record]);
        let residuals = self.scratch.path(&name("residuals"));
        let run = Run::new()
            .table(&self.table)
            .applied(&[self.klass])
            .plugins(self.plugins)
            .cache(&self.cache)
            .run_input(&input)
            .residuals(&residuals);
        let report = try_score(
            self.scratch,
            &reference,
            &[("candidate", &candidate)],
            &run,
            &BuiltInRegistry,
            &BuiltIn,
        )
        .unwrap_or_else(|stop| panic!("{label}: the run was stopped: {}", message(&stop)));
        let klass = self.klass;
        for projection in ["feature", "allele"] {
            let item = entry(&report, self.source, "candidate", projection);
            assert_eq!(counts(item), counted, "{label}: {projection}");
            assert_eq!(
                item["masked"],
                json!({klass: i64::from(taken)}),
                "{label}: {projection}"
            );
            assert_eq!(item["not_applied"], json!({}), "{label}: {projection}");
            assert_eq!(
                item["residual_rows"],
                if taken { 0 } else { self.left },
                "{label}: {projection}"
            );
        }
        assert_eq!(
            residual_lines(&residuals, self.source, "candidate", "feature").is_empty(),
            taken,
            "{label}"
        );
    }
}

// The DownstreamProtein of an insertion of `AC` after coding base 4 of the toy forward
// transcript as Ensembl VEP reads it at the span moved by two bases with the allele turned
// by two, and at the span moved by three bases with the allele as written; and as the
// candidate reads it at the record's own span.
const TWO_ON: &str = "TLLLSVVKGPNS";
const THREE_ON: &str = "NLLLSVVKGPNS";
const OWN_SPAN: &str = "NKLLLSVVKGPNS";

/// A frameshift row of a record that inserts after coding base 4: its coding columns, an
/// HGVS_OFFSET where `offset` gives one, and the Downstream plugin's two fields.
fn later_allele_row(
    allele: &str,
    coding: [&str; 3],
    downstream: &str,
    offset: Option<i64>,
) -> String {
    let mut extra = "IMPACT=HIGH;STRAND=1".to_string();
    if let Some(offset) = offset {
        extra += &format!(";HGVS_OFFSET={offset}");
    }
    extra += &format!(";DownstreamProtein={downstream};ProteinLengthChange=6");
    coding_row(
        "v",
        "21:1004-1005",
        allele,
        "frameshift_variant",
        coding,
        &extra,
    )
}

/// One cell of the grid of a later allele whose printed CDS pair an earlier allele's
/// notation left: its name, the record's alternates in order, the pair the reference
/// prints on the last allele's row, the reading its DownstreamProtein is, whether the
/// class takes the value, and the counts of the Downstream source. The last alternate is
/// the row asked.
type EarlierShiftCell = (
    &'static str,
    Alternates,
    &'static str,
    &'static str,
    bool,
    [i64; 6],
);

const EARLIER_SHIFT_GRID: [EarlierShiftCell; 12] = [
    (
        "the pair an earlier allele left, the value at the row's own offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "7-8",
        TWO_ON,
        true,
        [4, 4, 3, 3, 3, 3],
    ),
    (
        "the pair moved by the row's own offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "6-7",
        TWO_ON,
        true,
        [4, 4, 3, 3, 3, 3],
    ),
    (
        "an allele without an offset between the two",
        &[("A", Some(3)), ("AG", None), ("AC", Some(2))],
        "7-8",
        TWO_ON,
        true,
        [6, 6, 5, 5, 5, 5],
    ),
    (
        "the pair of the last earlier allele that has an offset",
        &[("A", Some(1)), ("AG", Some(3)), ("AC", Some(2))],
        "7-8",
        TWO_ON,
        true,
        [6, 6, 5, 5, 5, 5],
    ),
    (
        "the pair of the first allele where a later one has another",
        &[("A", Some(1)), ("AG", Some(3)), ("AC", Some(2))],
        "5-6",
        TWO_ON,
        false,
        [6, 6, 5, 6, 6, 5],
    ),
    (
        "a pair moved by neither offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "5-6",
        TWO_ON,
        false,
        [4, 4, 3, 4, 4, 3],
    ),
    (
        "the value read at the printed pair",
        &[("A", Some(3)), ("AC", Some(2))],
        "7-8",
        THREE_ON,
        false,
        [4, 4, 3, 4, 4, 3],
    ),
    (
        "the first allele's row",
        &[("AC", Some(2))],
        "7-8",
        TWO_ON,
        false,
        [2, 2, 1, 2, 2, 1],
    ),
    (
        "no offset on the row, the pair an earlier allele left, the allele as written",
        &[("A", Some(3)), ("AC", None)],
        "7-8",
        THREE_ON,
        true,
        [4, 4, 3, 3, 3, 3],
    ),
    (
        "no offset on the row and none on an earlier allele's",
        &[("A", None), ("AC", None)],
        "7-8",
        THREE_ON,
        false,
        [4, 4, 3, 4, 4, 3],
    ),
    (
        "no offset on the row, another value",
        &[("A", Some(3)), ("AC", None)],
        "7-8",
        TWO_ON,
        false,
        [4, 4, 3, 4, 4, 3],
    ),
    (
        "no offset on the row, the record's own pair printed",
        &[("A", Some(3)), ("AC", None)],
        "4-5",
        THREE_ON,
        false,
        [4, 4, 3, 4, 4, 3],
    ),
];

/// The rows of the earlier alleles are alike on both sides; the last allele's row prints
/// the cell's pair and value on the reference's side and the record's own pair with the
/// candidate's reading on its own. Given the candidate cache and the record's input line,
/// the class of the shifted CDS coordinates sets the DownstreamProtein aside in the cells
/// the grid marks and leaves it a residual in the others.
#[test]
fn score_sets_aside_the_cells_of_the_earlier_shift_grid_that_the_shifted_cds_class_takes() {
    let scratch = Scratch::new("earlier-shift-grid");
    let klass = PLUGIN_HGVS_SHIFTED_CDS_COORDINATES;
    let grid = Grid {
        scratch: &scratch,
        cache: downstream_cache(&scratch),
        klass,
        table: format!(
            "{klass}\tDownstream;NMD\tDownstreamProtein,ProteinLengthChange;NMD\tmet\t-\n"
        ),
        plugins: &["Downstream", "NMD"],
        source: "Downstream",
        left: 2,
    };
    for (number, (label, alternates, printed, reading, taken, counted)) in
        EARLIER_SHIFT_GRID.into_iter().enumerate()
    {
        let (earlier, last) = alternates.split_at(alternates.len() - 1);
        let (allele, offset) = last[0];
        assert_eq!(
            allele, "AC",
            "{label}: the readings are those of one allele"
        );
        let earlier: Vec<String> = earlier
            .iter()
            .map(|(other, shift)| later_allele_row(other, ["4-5", "2", "K/KX"], "KNS", *shift))
            .collect();
        let with = |asked: String| [earlier.clone(), vec![asked]].concat();
        let reference = with(later_allele_row(
            allele,
            [printed, "3", "T/KX"],
            reading,
            offset,
        ));
        let candidate = with(later_allele_row(
            allele,
            ["4-5", "2", "K/KX"],
            OWN_SPAN,
            offset,
        ));
        grid.holds(
            number,
            label,
            (&reference, &candidate),
            ["21", "1004", "v", "A", &spelt(alternates)],
            taken,
            counted,
        );
    }
}

/// The cells of a row of a record that inserts after coding base 6, the last base of
/// codon 2: its CDS position, protein position, amino acids and BLOSUM62 score.
type Blosum62Cells = (
    &'static str,
    &'static str,
    &'static str,
    Option<&'static str>,
);

/// The record's own pair, which has no reference residue, and the pair one base on with
/// the score of its residues.
const BLOSUM62_OWN: Blosum62Cells = ("6-7", "2-3", "-/K", None);
const BLOSUM62_MOVED: Blosum62Cells = ("7-8", "3", "T/K", Some("-1"));

fn blosum62_row(allele: &str, cells: Blosum62Cells, offset: Option<i64>) -> String {
    let (cds, protein, amino_acids, blosum) = cells;
    let mut extra = "IMPACT=MODERATE;STRAND=1".to_string();
    if let Some(offset) = offset {
        extra += &format!(";HGVS_OFFSET={offset}");
    }
    if let Some(blosum) = blosum {
        extra += &format!(";BLOSUM62={blosum}");
    }
    coding_row(
        "b",
        "21:1006-1007",
        allele,
        "inframe_insertion",
        [cds, protein, amino_acids],
        &extra,
    )
}

/// One cell of the grid of the BLOSUM62 score: its name, the record's alternates in
/// order, the last allele's row on the reference's side and on the candidate's, whether
/// the class takes the value, and the counts of the Blosum62 source. The last alternate is
/// the row asked: three bases between two codons.
type Blosum62Cell = (
    &'static str,
    Alternates,
    Blosum62Cells,
    Blosum62Cells,
    bool,
    [i64; 6],
);

const BLOSUM62_GRID: [Blosum62Cell; 10] = [
    (
        "a later allele read one base on, the first allele's offset",
        &[("C", Some(1)), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        true,
        [1, 0, 0, 0, 0, 0],
    ),
    (
        "a third allele, the first allele's offset",
        &[("C", Some(1)), ("G", None), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        true,
        [1, 0, 0, 0, 0, 0],
    ),
    (
        "the first allele's row",
        &[("AAA", Some(1))],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "a first allele without an offset",
        &[("C", None), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "a value that is not the matrix cell of its pair",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/K", Some("5")),
        BLOSUM62_OWN,
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "a reference residue that is not the moved codon's",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "S/K", Some("0")),
        BLOSUM62_OWN,
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "alternate residues that differ",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/R", Some("-1")),
        BLOSUM62_OWN,
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "a candidate pair the plugin scores",
        &[("C", Some(1)), ("AAA", None)],
        BLOSUM62_MOVED,
        ("6-7", "2-3", "N/K", None),
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "a pair moved by another size",
        &[("C", Some(1)), ("AAA", None)],
        ("8-9", "3", "T/K", Some("-1")),
        BLOSUM62_OWN,
        false,
        [1, 0, 0, 1, 0, 0],
    ),
    (
        "the value on the candidate's side alone",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/K", None),
        ("6-7", "2-3", "-/K", Some("-1")),
        false,
        [0, 1, 0, 0, 1, 0],
    ),
];

/// The rows of the earlier alleles are alike on both sides. Given the candidate cache and
/// the record's input line, the class of the shift shared across alleles sets the BLOSUM62
/// of the last allele's row aside in the cells the grid marks and leaves it a residual in
/// the others.
#[test]
fn score_sets_aside_the_cells_of_the_blosum62_grid_that_the_shared_shift_class_takes() {
    let scratch = Scratch::new("blosum62-grid");
    let klass = BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES;
    let grid = Grid {
        scratch: &scratch,
        cache: downstream_cache(&scratch),
        klass,
        table: format!("{klass}\tBlosum62\tBLOSUM62\tmet\t-\n"),
        plugins: &["Blosum62"],
        source: "Blosum62",
        left: 1,
    };
    for (number, (label, alternates, reference, candidate, taken, counted)) in
        BLOSUM62_GRID.into_iter().enumerate()
    {
        let (earlier, last) = alternates.split_at(alternates.len() - 1);
        let (allele, offset) = last[0];
        let earlier: Vec<String> = earlier
            .iter()
            .map(|(other, shift)| blosum62_row(other, ("6-7", "2-3", "-/X", None), *shift))
            .collect();
        let with = |cells: Blosum62Cells| {
            [earlier.clone(), vec![blosum62_row(allele, cells, offset)]].concat()
        };
        grid.holds(
            number,
            label,
            (&with(reference), &with(candidate)),
            ["21", "1006", "b", "A", &spelt(alternates)],
            taken,
            counted,
        );
    }
}
