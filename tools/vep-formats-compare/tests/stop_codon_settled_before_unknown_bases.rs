// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `stop_codon_settled_before_unknown_bases`: the rows of an insertion whose alternate
//! allele holds unknown bases after bases that settle a stop codon. The reference writes no
//! residues, the term of an unknown change or `frameshift_variant` alone, and no protein
//! notation; the candidate writes the residues the codons give, `stop_gained` and the
//! notation. Where the reference writes `stop_retained_variant` for a stop codon test it
//! does not make, the row is the class's as well. Where the settled stop codon stands at the
//! transcript's own, the reference writes `stop_retained_variant`, `frameshift_variant` or
//! the term of an unknown change, and the candidate `inframe_insertion` beside
//! `stop_retained_variant`, the residues and the notation of a stop codon left unchanged.
//! Two default-format outputs of one row are compared under the class, without it and off
//! its shape, and the rows of a kept stop codon in the tab, VCF and JSON formats as well.
//!
//! An insertion written in an intron that the HGVS 3' shift carries into the coding
//! sequence has a row without a `Codons` cell: the reference writes no notation, and the
//! candidate the protein notation read where the shift ends and the HGVS_OFFSET beside it.
//! The class reads such a row beside a made-up candidate cache, on both strands, and leaves
//! it a difference without the cache and off its shape.

mod compare_common;
mod support;
mod term_pairs_common;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use compare_common::*;
use support::{args_of, compact, member, report, run};
use term_pairs_common::{csq_header, most_severe, open_info, quoted, quoted_list};
use vep_formats_compare::classes;
use vep_formats_compare::compare::{compare_entries, CompareOptions};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

const CLASS: &str = "stop_codon_settled_before_unknown_bases";

/// The class of an insertion that keeps the stop codon, whose rule stands before this
/// class's among the row classes.
const OVERLAP: &str = "stop_codon_overlap_insertion_range";

/// Four unknown bases.
const UNKNOWN: &str = "NNNN";

/// Where a row stands: its Location, its gene, transcript and protein, its three position
/// columns and its strand.
#[derive(Clone, Copy)]
struct Site {
    location: &'static str,
    gene: &'static str,
    transcript: &'static str,
    protein: &'static str,
    cdna: &'static str,
    cds: &'static str,
    protein_position: &'static str,
    strand: &'static str,
}

/// A codon boundary inside the coding sequence of the forward-strand ENST00000622113.
const INSIDE: Site = Site {
    location: "21:10567700-10567701",
    gene: "ENSG00000274391",
    transcript: "ENST00000622113",
    protein: "ENSP00000482040.1",
    cdna: "775-776",
    cds: "423-424",
    protein_position: "141-142",
    strand: "1",
};

/// Between the last sense codon of that transcript and its stop codon.
const BEFORE_THE_STOP: Site = Site {
    location: "21:10605549-10605550",
    cdna: "1951-1952",
    cds: "1599-1600",
    protein_position: "533-534",
    ..INSIDE
};

/// Between the stop codon of the reverse-strand ENST00000681601 and its last sense codon.
const REVERSE_STOP: Site = Site {
    location: "21:14108995-14108996",
    gene: "ENSG00000188992",
    transcript: "ENST00000681601",
    protein: "ENSP00000505323.1",
    cdna: "1490-1491",
    cds: "1380-1381",
    protein_position: "460-461",
    strand: "-1",
};

/// One Transcript row of an insertion.
#[derive(Clone)]
struct Row {
    /// The first column: an identifier, or the record's position and alleles.
    name: &'static str,
    site: Site,
    /// The inserted bases as the row writes them.
    allele: String,
    terms: &'static str,
    residues: &'static str,
    /// The Codons cell; no reference codon and the allele's own bases where `None`.
    codons: Option<String>,
    /// The Extra text.
    extra: String,
}

impl Row {
    fn codons(&self) -> String {
        self.codons
            .clone()
            .unwrap_or_else(|| format!("-/{}", self.allele))
    }

    /// The row's line of a default-format output.
    fn line(&self) -> String {
        format!("{}\t{}", self.columns(), self.extra)
    }

    /// The columns before the Extra text.
    fn columns(&self) -> String {
        let site = &self.site;
        format!(
            "{}\t{}\t{}\t{}\t{}\tTranscript\t{}\t{}\t{}\t{}\t{}\t{}\t-",
            self.name,
            site.location,
            self.allele,
            site.gene,
            site.transcript,
            self.terms,
            site.cdna,
            site.cds,
            site.protein_position,
            self.residues,
            self.codons()
        )
    }

    /// The value the Extra text holds under `name`, empty where it holds none.
    fn extra_value(&self, name: &str) -> &str {
        self.extra
            .split(';')
            .filter_map(|pair| pair.split_once('='))
            .find(|(key, _)| *key == name)
            .map_or("", |(_, value)| value)
    }

    /// The input line of the row's record: the base before the insertion point on the
    /// forward-strand transcript, then that base and the inserted bases.
    fn input_line(&self) -> String {
        format!("21\t10605549\t.\tT\tT{}\t.\t.\t.", self.allele)
    }

    /// The row's entries in the tab, VCF or JSON format: a row before the forward-strand
    /// stop codon, each format writing a cell the row leaves out in its own way.
    fn in_format(&self, fmt: &str) -> Vec<Entry> {
        let site = &self.site;
        let inputs = [self.input_line()];
        let or_dash = |value: &str| {
            if value.is_empty() {
                "-".to_string()
            } else {
                value.to_string()
            }
        };
        let or_empty = |value: &str| if value == "-" { "" } else { value }.to_string();
        let (impact, notation, offset) = (
            self.extra_value("IMPACT"),
            self.extra_value("HGVSp"),
            self.extra_value("HGVS_OFFSET"),
        );
        match fmt {
            "tab" => {
                let header =
                    DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSp\tHGVS_OFFSET");
                let row = [
                    self.columns(),
                    impact.to_string(),
                    site.strand.to_string(),
                    or_dash(notation),
                    or_dash(offset),
                ]
                .join("\t");
                parse_tab(&format!("{header}\n{row}\n")).expect("the output parses")
            }
            "vcf" => {
                let csq = [
                    self.allele.clone(),
                    self.terms.replace(',', "&"),
                    impact.to_string(),
                    site.transcript.to_string(),
                    "Transcript".to_string(),
                    site.cdna.to_string(),
                    site.cds.to_string(),
                    site.protein_position.to_string(),
                    or_empty(self.residues),
                    self.codons(),
                    notation.to_string(),
                    offset.to_string(),
                ]
                .join("|");
                let fields = "Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|HGVSp|HGVS_OFFSET";
                let input_line = self.input_line();
                let text = format!(
                    "{}{}CSQ={csq}\n",
                    csq_header(fields),
                    open_info(&input_line)
                );
                parse_vcf(&text, &inputs).expect("the output parses")
            }
            _ => {
                let held: Vec<&str> = self.terms.split(',').collect();
                let mut entry = format!(
                    "{{\"variant_allele\": {}, \"transcript_id\": {}, \"consequence_terms\": {}, \"impact\": {}, \"strand\": {}, \"codons\": {}",
                    quoted(&self.allele),
                    quoted(site.transcript),
                    quoted_list(&held),
                    quoted(impact),
                    site.strand,
                    quoted(&self.codons()),
                );
                for (key, column) in [
                    ("cdna", site.cdna),
                    ("cds", site.cds),
                    ("protein", site.protein_position),
                ] {
                    let (first, last) = column.split_once('-').unwrap_or((column, column));
                    entry.push_str(&format!(
                        ", \"{key}_start\": {first}, \"{key}_end\": {last}"
                    ));
                }
                if self.residues != "-" {
                    entry.push_str(&format!(", \"amino_acids\": {}", quoted(self.residues)));
                }
                if !notation.is_empty() {
                    // The JSON format writes the notation's `=` as it stands.
                    entry.push_str(&format!(
                        ", \"hgvsp\": {}",
                        quoted(&notation.replace("%3D", "="))
                    ));
                }
                if !offset.is_empty() {
                    entry.push_str(&format!(", \"hgvs_offset\": {offset}"));
                }
                entry.push('}');
                let text = format!(
                    "{{\"input\": {}, \"seq_region_name\": \"21\", \"start\": 10605550, \"end\": 10605549, \"allele_string\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
                    quoted(&self.input_line()),
                    quoted(&format!("-/{}", self.allele)),
                    quoted(most_severe(&held)),
                );
                parse_json(&text, &inputs).expect("the output parses")
            }
        }
    }
}

/// The bases of the other strand, read in its direction.
fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// The reference's row of the bases `read` inserted at `site`, as the transcript reads
/// them: no residues, `terms` and their `impact`. The row of a reverse-strand transcript
/// writes the other strand's bases as its allele.
fn at(site: Site, name: &'static str, read: &str, terms: &'static str, impact: &str) -> Row {
    let allele = if site.strand == "-1" {
        reverse_complement(read)
    } else {
        read.to_string()
    };
    Row {
        name,
        site,
        allele,
        terms,
        residues: "-",
        codons: Some(format!("-/{read}")),
        extra: format!("IMPACT={impact};STRAND={}", site.strand),
    }
}

/// The reference's row of the inserted bases `known`, the unknown bases and `tail` inside
/// the coding sequence: no residues, `terms` and their `impact`.
fn reference(
    name: &'static str,
    known: &str,
    tail: &str,
    terms: &'static str,
    impact: &str,
) -> Row {
    Row {
        codons: None,
        ..at(
            INSIDE,
            name,
            &format!("{known}{UNKNOWN}{tail}"),
            terms,
            impact,
        )
    }
}

/// The candidate's row of the same record: `terms`, their `impact`, `residues` and, where
/// `notation` is given, the protein notation, with `offset` where the notation is read past
/// the record's own placement.
fn written(
    of: &Row,
    terms: &'static str,
    impact: &str,
    residues: &'static str,
    notation: &str,
    offset: &str,
) -> Row {
    let site = &of.site;
    let notation = if notation.is_empty() {
        String::new()
    } else {
        format!(";HGVSp={}:p.{notation}", site.protein)
    };
    let offset = if offset.is_empty() {
        String::new()
    } else {
        format!(";HGVS_OFFSET={offset}")
    };
    Row {
        terms,
        residues,
        extra: format!("IMPACT={impact};STRAND={}{notation}{offset}", site.strand),
        ..of.clone()
    }
}

/// The candidate's row of a gained stop: `residues`, `terms`, IMPACT `HIGH` and the protein
/// notation `notation`, with `offset` where the notation is read past the record's own
/// placement.
fn candidate(
    of: &Row,
    terms: &'static str,
    residues: &'static str,
    notation: &str,
    offset: &str,
) -> Row {
    written(of, terms, "HIGH", residues, notation, offset)
}

/// The verdict of the two outputs of one row under `classes` and what was applied.
fn compared(reference: &Row, candidate: &Row, classes: &[&str]) -> (String, BTreeMap<String, i64>) {
    let (report, documented) = compare_default(
        &default_text(&[reference.line()]),
        &default_text(&[candidate.line()]),
        classes,
    );
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// A sense codon and a stop codon before the unknown bases, ten inserted bases.
fn frameshift() -> (Row, Row) {
    let reference = reference(
        "stop_then_unknown_frameshift",
        "GAATAG",
        "",
        "frameshift_variant",
        "HIGH",
    );
    let candidate = candidate(
        &reference,
        "stop_gained,frameshift_variant",
        "-/E*XX",
        "Ile142GlufsTer2",
        "",
    );
    (reference, candidate)
}

/// The same with two known bases after the unknown ones, twelve inserted bases.
fn inframe() -> (Row, Row) {
    let reference = reference(
        "stop_then_unknown_inframe",
        "GAATAG",
        "CA",
        "coding_sequence_variant",
        "MODIFIER",
    );
    let candidate = candidate(
        &reference,
        "stop_gained,inframe_insertion",
        "-/E*XX",
        "Asn141_Ile142insGluTer",
        "",
    );
    (reference, candidate)
}

/// An insertion whose notation is read four bases 3' of the record's placement.
fn shifted() -> (Row, Row) {
    let reference = reference(
        "stop_then_unknown_shifted",
        "ATTTAG",
        "CA",
        "coding_sequence_variant",
        "MODIFIER",
    );
    let candidate = candidate(
        &reference,
        "stop_gained,inframe_insertion",
        "-/I*XX",
        "Ile142_Leu143insTer",
        "4",
    );
    (reference, candidate)
}

/// A stop codon and the unknown bases as the transcript reads them, seven inserted bases.
const KEPT_SHIFTING: &str = "TAGNNNN";

/// The same with two known bases after the unknown ones, nine inserted bases.
const KEPT_WHOLE_CODONS: &str = "TAGNNNNCA";

/// The candidate's terms on a row whose inserted bases keep the stop codon.
const KEPT: &str = "inframe_insertion,stop_retained_variant";

/// The candidate's row of a kept stop codon: `KEPT`, IMPACT `MODERATE`, `residues` and the
/// notation of the stop codon left unchanged at the residue `stop`, read `offset` bases
/// past the record's own placement.
fn kept(of: &Row, residues: &'static str, stop: &str, offset: &str) -> Row {
    let notation = if stop.is_empty() {
        String::new()
    } else {
        format!("Ter{stop}%3D")
    };
    written(of, KEPT, "MODERATE", residues, &notation, offset)
}

/// The inserted bases before the stop codon of the forward-strand transcript, where the
/// reference writes `terms` and their `impact`.
fn kept_forward(name: &'static str, read: &str, terms: &'static str, impact: &str) -> (Row, Row) {
    let reference = at(BEFORE_THE_STOP, name, read, terms, impact);
    let candidate = kept(&reference, "-/*XX", "534", "2");
    (reference, candidate)
}

/// The inserted bases before the stop codon of the reverse-strand transcript, where the
/// reference writes `stop_retained_variant`.
fn kept_reverse(name: &'static str, read: &str) -> (Row, Row) {
    let reference = at(REVERSE_STOP, name, read, "stop_retained_variant", "LOW");
    let candidate = kept(&reference, "-/*XX", "461", "-2");
    (reference, candidate)
}

/// Seven bases before the forward-strand stop codon: the reference's `frameshift_variant`.
fn kept_frameshift() -> (Row, Row) {
    kept_forward(
        "kept_stop_frameshift",
        KEPT_SHIFTING,
        "frameshift_variant",
        "HIGH",
    )
}

/// Nine bases there: the reference's term of an unknown change.
fn kept_whole_codons() -> (Row, Row) {
    kept_forward(
        "kept_stop_whole_codons",
        KEPT_WHOLE_CODONS,
        "coding_sequence_variant",
        "MODIFIER",
    )
}

/// The rows of a kept stop codon whose reference writes `frameshift_variant`, the term the
/// earlier class's rule reads as well.
fn kept_against_a_frameshift() -> Vec<(Row, Row)> {
    // A term both rows hold, and no protein notation.
    let beside_a_shared_term = {
        let reference = at(
            BEFORE_THE_STOP,
            "kept_stop_beside_a_shared_term",
            KEPT_SHIFTING,
            "frameshift_variant,NMD_transcript_variant",
            "HIGH",
        );
        let candidate = written(
            &reference,
            "inframe_insertion,stop_retained_variant,NMD_transcript_variant",
            "MODERATE",
            "-/*XX",
            "",
            "",
        );
        (reference, candidate)
    };
    // Eight bases inside the last sense codon, which the alternate codons spell again
    // before the stop codon.
    let inside_the_last_codon = {
        let site = Site {
            location: "21:10605548-10605549",
            cdna: "1950-1951",
            cds: "1598-1599",
            protein_position: "533",
            ..BEFORE_THE_STOP
        };
        let reference = Row {
            codons: Some(format!("gat/gaTTAG{UNKNOWN}t")),
            ..at(
                site,
                "kept_stop_inside_the_last_codon",
                &format!("TTAG{UNKNOWN}"),
                "frameshift_variant",
                "HIGH",
            )
        };
        let candidate = kept(&reference, "D/D*XX", "", "");
        (reference, candidate)
    };
    // The record's own alleles in the first column.
    let spelt = kept_forward(
        "21_10605549_T/TTAGNNNN",
        KEPT_SHIFTING,
        "frameshift_variant",
        "HIGH",
    );
    vec![
        kept_frameshift(),
        beside_a_shared_term,
        inside_the_last_codon,
        spelt,
    ]
}

/// The rows of a kept stop codon no earlier class's rule reads: the reference's term of an
/// unknown change and its `stop_retained_variant`.
fn kept_against_another_term() -> Vec<(Row, Row)> {
    // Six bases inside the stop codon, whose first base opens the alternate codons.
    let inside_the_stop_codon = {
        let site = Site {
            location: "21:10605550-10605551",
            cdna: "1952-1953",
            cds: "1600-1601",
            protein_position: "534",
            ..BEFORE_THE_STOP
        };
        let reference = Row {
            codons: Some(format!("taa/tAG{UNKNOWN}aa")),
            ..at(
                site,
                "kept_stop_inside_the_stop_codon",
                &format!("AG{UNKNOWN}"),
                "coding_sequence_variant",
                "MODIFIER",
            )
        };
        let candidate = kept(&reference, "*/*XX", "", "");
        (reference, candidate)
    };
    // Three bases, the last unknown, inside a stop codon an intron splits: the first of
    // them is the letter of the exon base that ends the codon, and both rows hold
    // `splice_region_variant`.
    let beside_a_splice_term = {
        let site = Site {
            location: "21:35201927-35201928",
            gene: "ENSG1",
            transcript: "ENST00000399326",
            protein: "ENSP1.1",
            cdna: "3048-3049",
            cds: "3047-3048",
            protein_position: "1016",
            strand: "1",
        };
        let reference = Row {
            codons: Some("taa/taAANa".to_string()),
            ..at(
                site,
                "21_35201927_G/GAAN",
                "AAN",
                "splice_region_variant,coding_sequence_variant",
                "LOW",
            )
        };
        let candidate = written(
            &reference,
            "inframe_insertion,splice_region_variant,stop_retained_variant",
            "MODERATE",
            "*/*X",
            "",
            "",
        );
        (reference, candidate)
    };
    vec![
        kept_whole_codons(),
        inside_the_stop_codon,
        beside_a_splice_term,
        kept_reverse("kept_stop_reverse_frameshift", KEPT_SHIFTING),
        kept_reverse("kept_stop_reverse_whole_codons", KEPT_WHOLE_CODONS),
    ]
}

fn kept_rows() -> Vec<(Row, Row)> {
    let mut rows = kept_against_a_frameshift();
    rows.extend(kept_against_another_term());
    rows
}

fn nothing() -> BTreeMap<String, i64> {
    count_map(&[])
}

/// The verdict of the two rows in one of the tab, VCF and JSON formats under `classes` and
/// what was applied.
fn compared_in(
    fmt: &str,
    reference: &Row,
    candidate: &Row,
    classes: &[&str],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare(
        &reference.in_format(fmt),
        &candidate.in_format(fmt),
        &mut documented,
    );
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// Each row is set aside under the class alone, its terms, residues, IMPACT, protein
/// notation and offset together; without the class it is a difference.
#[test]
fn a_row_of_the_class_is_set_aside_with_every_cell_the_class_reads() {
    for (reference, candidate) in [frameshift(), inframe(), shifted()] {
        let (verdict, applied) = compared(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict, "green", "{}", reference.name);
        assert_eq!(
            applied.keys().map(String::as_str).collect::<Vec<_>>(),
            [CLASS],
            "{}",
            reference.name
        );
        assert_eq!(
            compared(&reference, &candidate, &[]),
            ("red".to_string(), nothing()),
            "{}",
            reference.name
        );
    }
}

/// The sides exchanged are no row of the class, on a row it takes as written.
#[test]
fn the_sides_exchanged_stay_a_difference() {
    for (reference, candidate) in [frameshift(), inframe(), shifted()] {
        assert_eq!(
            compared(&reference, &candidate, &[CLASS]).0,
            "green",
            "{}",
            reference.name
        );
        assert_eq!(
            compared(&candidate, &reference, &[CLASS]),
            ("red".to_string(), nothing()),
            "{}",
            reference.name
        );
    }
}

/// A row off the shape stays a difference under the class and has no cell set aside; a
/// cell the class does not read stays a difference on a row it takes.
#[test]
fn a_row_off_the_shape_stays_a_difference() {
    let (reference, candidate) = inframe();
    // Each case ends with whether the class takes the row.
    let off: Vec<(&str, Row, Row, bool)> = vec![
        (
            "residues that are not the codons' translation",
            reference.clone(),
            Row {
                residues: "-/E*XK",
                ..candidate.clone()
            },
            false,
        ),
        (
            "known bases that settle no stop codon before the unknown ones",
            Row {
                allele: format!("GAAGCT{UNKNOWN}CA"),
                ..reference.clone()
            },
            Row {
                allele: format!("GAAGCT{UNKNOWN}CA"),
                residues: "-/EAXX",
                ..candidate.clone()
            },
            false,
        ),
        (
            "codons that differ between the sides",
            reference.clone(),
            Row {
                codons: Some(format!("-/GAATAA{UNKNOWN}CA")),
                ..candidate.clone()
            },
            false,
        ),
        (
            "a term other than the gained stop and an in-frame term",
            reference.clone(),
            Row {
                terms: "stop_gained,missense_variant",
                ..candidate.clone()
            },
            false,
        ),
        (
            "no gained stop",
            reference.clone(),
            Row {
                terms: "inframe_insertion",
                ..candidate.clone()
            },
            false,
        ),
        (
            "a protein notation that names no stop",
            reference.clone(),
            Row {
                extra: format!(
                    "IMPACT=HIGH;STRAND=1;HGVSp={}:p.Asn141_Ile142insGlu",
                    INSIDE.protein
                ),
                ..candidate.clone()
            },
            true,
        ),
        (
            "a coding notation the candidate alone writes",
            reference.clone(),
            Row {
                extra: format!(
                    "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000622113.4:c.423_424insGAATAG;HGVSp={}:p.Asn141_Ile142insGluTer",
                    INSIDE.protein
                ),
                ..candidate.clone()
            },
            true,
        ),
        (
            "an IMPACT that is not the terms'",
            reference.clone(),
            Row {
                extra: format!(
                    "IMPACT=MODERATE;STRAND=1;HGVSp={}:p.Asn141_Ile142insGluTer",
                    INSIDE.protein
                ),
                ..candidate.clone()
            },
            true,
        ),
        (
            "residues in the reference",
            Row {
                residues: "-/E*XX",
                ..reference.clone()
            },
            candidate.clone(),
            false,
        ),
    ];
    for (what, reference, candidate, taken) in off {
        let (verdict, applied) = compared(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict, "red", "{what}");
        assert_eq!(applied.contains_key(CLASS), taken, "{what}");
    }
}

/// A row whose inserted bases keep the stop codon is set aside under the class alone, its
/// terms, residues, IMPACT, protein notation and offset together: the candidate's
/// `inframe_insertion` beside the `stop_retained_variant` both rows hold, or the two terms
/// against the reference's `frameshift_variant` or term of an unknown change on an insertion
/// between two neighbouring positions. The residues keep the reference residue of a codon
/// the insertion stands inside. Without the class each row is a difference.
#[test]
fn a_kept_stop_row_is_set_aside_with_every_cell_the_class_reads() {
    let rows = kept_rows();
    assert_eq!(rows.len(), 9);
    for (reference, candidate) in rows {
        assert_eq!(
            compared(&reference, &candidate, &[CLASS]),
            ("green".to_string(), count_map(&[(CLASS, 1)])),
            "{}",
            reference.name
        );
        assert_eq!(
            compared(&reference, &candidate, &[]),
            ("red".to_string(), nothing()),
            "{}",
            reference.name
        );
        assert_eq!(
            compared(&candidate, &reference, &[CLASS]),
            ("red".to_string(), nothing()),
            "{} with its sides exchanged",
            reference.name
        );
    }
}

/// Beside the class of an insertion that keeps the stop codon, whose rule reads the
/// reference's `frameshift_variant` against the candidate's two terms, such a row is that
/// class's by its terms and this class's by its residues, its protein notation and its
/// offset; the other kept stop rows are this class's alone. Every row class of the table
/// together sets each row aside in the same way.
#[test]
fn a_kept_stop_row_is_set_aside_beside_the_class_of_the_kept_stop_codon() {
    let both = [OVERLAP, CLASS];
    let every_row_class = classes::tables().term_pair_classes().to_vec();
    assert!(every_row_class.contains(&OVERLAP) && every_row_class.contains(&CLASS));
    for (reference, candidate) in kept_rows() {
        let with_both = compared(&reference, &candidate, &both);
        assert_eq!(with_both.0, "green", "{}", reference.name);
        assert_eq!(
            compared(&reference, &candidate, &every_row_class),
            with_both,
            "{}",
            reference.name
        );
    }
    for (reference, candidate) in kept_against_a_frameshift() {
        assert_eq!(
            compared(&reference, &candidate, &both),
            ("green".to_string(), count_map(&[(CLASS, 1), (OVERLAP, 1)])),
            "{}",
            reference.name
        );
        // The earlier class reads no residues, so the row stays a difference under it.
        assert_eq!(
            compared(&reference, &candidate, &[OVERLAP]).0,
            "red",
            "{}",
            reference.name
        );
    }
    for (reference, candidate) in kept_against_another_term() {
        assert_eq!(
            compared(&reference, &candidate, &both),
            ("green".to_string(), count_map(&[(CLASS, 1)])),
            "{}",
            reference.name
        );
    }
    // A gained stop is no row of the earlier class.
    for (reference, candidate) in [frameshift(), inframe(), shifted()] {
        assert_eq!(
            compared(&reference, &candidate, &both),
            ("green".to_string(), count_map(&[(CLASS, 1)])),
            "{}",
            reference.name
        );
    }
}

/// In the tab, VCF and JSON formats a kept stop row before the forward-strand stop codon is
/// set aside as in the default format: under the class alone, and beside the earlier class,
/// which takes the terms of the row whose reference writes `frameshift_variant`. Without
/// the class each row is a difference.
#[test]
fn a_kept_stop_row_is_set_aside_in_the_tab_vcf_and_json_formats() {
    let alone = count_map(&[(CLASS, 1)]);
    let beside = count_map(&[(CLASS, 1), (OVERLAP, 1)]);
    let rows = [(kept_frameshift(), &beside), (kept_whole_codons(), &alone)];
    for fmt in ["tab", "vcf", "json"] {
        for ((reference, candidate), with_both) in &rows {
            let label = format!("{fmt} {}", reference.name);
            assert_eq!(
                compared_in(fmt, reference, candidate, &[CLASS]),
                ("green".to_string(), alone.clone()),
                "{label}"
            );
            assert_eq!(
                compared_in(fmt, reference, candidate, &[OVERLAP, CLASS]),
                ("green".to_string(), (*with_both).clone()),
                "{label}"
            );
            assert_eq!(
                compared_in(fmt, reference, candidate, &[]),
                ("red".to_string(), nothing()),
                "{label}"
            );
            assert_eq!(
                compared_in(fmt, candidate, reference, &[OVERLAP, CLASS]),
                ("red".to_string(), nothing()),
                "{label} with its sides exchanged"
            );
        }
    }
}

/// A kept stop row off the shape stays a difference under the class and has no cell set
/// aside; a cell the class does not read stays a difference on a row it takes.
#[test]
fn a_kept_stop_row_off_the_shape_stays_a_difference() {
    let (shifting, kept_shifting) = kept_frameshift();
    let (whole, kept_whole) = kept_whole_codons();
    let (retained, kept_retained) = kept_reverse("kept_stop_reverse_frameshift", KEPT_SHIFTING);
    for (reference, candidate) in [
        (&shifting, &kept_shifting),
        (&whole, &kept_whole),
        (&retained, &kept_retained),
    ] {
        assert_eq!(
            compared(reference, candidate, &[CLASS]).0,
            "green",
            "{}",
            reference.name
        );
    }
    let moved = |site: Site, row: &Row| Row {
        site,
        ..row.clone()
    };
    let three_positions = Site {
        location: "21:10605549-10605551",
        ..BEFORE_THE_STOP
    };
    let no_cdna_end = Site {
        cdna: "1951-?",
        ..BEFORE_THE_STOP
    };
    let no_cds_start = Site {
        cds: "?-1600",
        ..BEFORE_THE_STOP
    };
    let replacement = |row: &Row| Row {
        name: "21_10605549_TA/TTAGNNNNCA",
        ..row.clone()
    };
    // Each case ends with whether the class takes the row.
    let off: Vec<(&str, Row, Row, bool)> = vec![
        (
            "the reference's missense_variant",
            Row {
                terms: "missense_variant",
                ..whole.clone()
            },
            kept_whole.clone(),
            false,
        ),
        (
            "the reference's inframe_insertion",
            Row {
                terms: "inframe_insertion",
                ..retained.clone()
            },
            kept_retained.clone(),
            false,
        ),
        (
            "the reference's inframe_insertion beside its frameshift_variant",
            Row {
                terms: "inframe_insertion,frameshift_variant",
                ..shifting.clone()
            },
            kept_shifting.clone(),
            false,
        ),
        (
            "a term the reference alone holds beside its stop_retained_variant",
            Row {
                terms: "stop_retained_variant,splice_region_variant",
                ..retained.clone()
            },
            kept_retained.clone(),
            false,
        ),
        (
            "a second term the reference alone holds",
            Row {
                terms: "stop_gained,frameshift_variant",
                ..shifting.clone()
            },
            kept_shifting.clone(),
            false,
        ),
        (
            "a third term the candidate alone holds",
            whole.clone(),
            Row {
                terms: "inframe_insertion,stop_retained_variant,splice_region_variant",
                ..kept_whole.clone()
            },
            false,
        ),
        (
            "no stop_retained_variant in the candidate",
            whole.clone(),
            Row {
                terms: "inframe_insertion",
                ..kept_whole.clone()
            },
            false,
        ),
        (
            "no inframe_insertion in the candidate",
            whole.clone(),
            Row {
                terms: "stop_retained_variant",
                ..kept_whole.clone()
            },
            false,
        ),
        (
            "a Location of three positions",
            moved(three_positions, &whole),
            moved(three_positions, &kept_whole),
            false,
        ),
        (
            "a record whose reference allele keeps a base",
            replacement(&whole),
            replacement(&kept_whole),
            false,
        ),
        (
            "no cDNA end",
            moved(no_cdna_end, &shifting),
            moved(no_cdna_end, &kept_shifting),
            false,
        ),
        (
            "no CDS start",
            moved(no_cds_start, &whole),
            moved(no_cds_start, &kept_whole),
            false,
        ),
        (
            "residues that are not the codons' translation",
            retained.clone(),
            Row {
                residues: "-/*XK",
                ..kept_retained.clone()
            },
            false,
        ),
        (
            "residues without a reference half",
            whole.clone(),
            Row {
                residues: "*XX",
                ..kept_whole.clone()
            },
            false,
        ),
        (
            "residues in the reference",
            Row {
                residues: "-/*XX",
                ..retained.clone()
            },
            kept_retained.clone(),
            false,
        ),
        (
            "an unknown base in the stop codon",
            Row {
                allele: format!("TA{UNKNOWN}GCA"),
                codons: None,
                ..whole.clone()
            },
            Row {
                allele: format!("TA{UNKNOWN}GCA"),
                codons: None,
                residues: "-/XXA",
                ..kept_whole.clone()
            },
            false,
        ),
        (
            "codons that differ between the sides",
            whole.clone(),
            Row {
                codons: Some(format!("-/TAA{UNKNOWN}CA")),
                ..kept_whole.clone()
            },
            false,
        ),
        (
            "an IMPACT that is not the terms'",
            whole.clone(),
            Row {
                extra: format!(
                    "IMPACT=LOW;STRAND=1;HGVSp={}:p.Ter534%3D;HGVS_OFFSET=2",
                    BEFORE_THE_STOP.protein
                ),
                ..kept_whole.clone()
            },
            true,
        ),
        (
            "a protein notation that names no stop",
            whole.clone(),
            Row {
                extra: format!(
                    "IMPACT=MODERATE;STRAND=1;HGVSp={}:p.Asp533dup",
                    BEFORE_THE_STOP.protein
                ),
                ..kept_whole.clone()
            },
            true,
        ),
        (
            "an offset beside no protein notation",
            retained.clone(),
            Row {
                extra: "IMPACT=MODERATE;STRAND=-1;HGVS_OFFSET=-2".to_string(),
                ..kept_retained.clone()
            },
            true,
        ),
        (
            "a coding notation the candidate alone writes",
            shifting.clone(),
            Row {
                extra: format!(
                    "IMPACT=MODERATE;STRAND=1;HGVSc=ENST00000622113.4:c.1599_1600insTAG;HGVSp={}:p.Ter534%3D;HGVS_OFFSET=2",
                    BEFORE_THE_STOP.protein
                ),
                ..kept_shifting.clone()
            },
            true,
        ),
    ];
    for (what, reference, candidate, taken) in off {
        let (verdict, applied) = compared(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict, "red", "{what}");
        assert_eq!(applied.contains_key(CLASS), taken, "{what}");
    }
}

/// The class of a `stop_retained_variant` written for a stop codon test that is not made.
const UNTESTED: &str = "stop_retained_without_stop_codon_test";

/// The reference's row of the same record where it writes `stop_retained_variant` for a
/// stop codon test it does not make: that term in the place of its own, and its IMPACT.
fn untested(of: &Row) -> Row {
    Row {
        terms: "stop_retained_variant",
        extra: "IMPACT=LOW;STRAND=1".to_string(),
        ..of.clone()
    }
}

/// The gained stop with no term beside it, as on an insertion inside a codon.
fn alone() -> (Row, Row) {
    let (reference, candidate) = inframe();
    let name = "stop_then_unknown_inside_a_codon";
    let reference = Row { name, ..reference };
    let candidate = Row {
        name,
        terms: "stop_gained",
        ..candidate
    };
    (reference, candidate)
}

/// A row on which the reference writes `stop_retained_variant` without its stop codon test
/// is set aside under the class with every cell the class reads, whether the candidate
/// writes `frameshift_variant`, an in-frame term or nothing beside the gained stop. The
/// class of the untested term does not take it, alone or named first.
#[test]
fn a_reference_row_of_the_untested_stop_retained_term_is_set_aside_under_the_class() {
    for (reference, candidate) in [frameshift(), inframe(), shifted(), alone()] {
        let reference = untested(&reference);
        let name = reference.name;
        for classes in [&[CLASS][..], &[UNTESTED, CLASS][..]] {
            let (verdict, applied) = compared(&reference, &candidate, classes);
            assert_eq!(verdict, "green", "{name} {classes:?}");
            assert_eq!(
                applied.keys().map(String::as_str).collect::<Vec<_>>(),
                [CLASS],
                "{name} {classes:?}"
            );
        }
        for classes in [&[][..], &[UNTESTED][..]] {
            assert_eq!(
                compared(&reference, &candidate, classes),
                ("red".to_string(), nothing()),
                "{name} {classes:?}"
            );
        }
        assert_eq!(
            compared(&candidate, &reference, &[CLASS]),
            ("red".to_string(), nothing()),
            "{name} exchanged"
        );
    }
}

/// Beside the reference's untested `stop_retained_variant` the class takes the gained stop
/// with one term at most, on an insertion's or a deletion's row whose allele shows that the
/// test is not made.
#[test]
fn a_reference_row_of_stop_retained_off_the_shape_stays_a_difference() {
    let (reference, candidate) = frameshift();
    let reference = untested(&reference);
    let known = "GAATAGACGT".to_string();
    let off: Vec<(&str, Row, Row)> = vec![
        (
            "stop_retained_variant kept beside the gained stop",
            reference.clone(),
            Row {
                terms: "stop_gained,frameshift_variant,stop_retained_variant",
                ..candidate.clone()
            },
        ),
        (
            "two terms beside the gained stop",
            reference.clone(),
            Row {
                terms: "stop_gained,frameshift_variant,inframe_insertion",
                ..candidate.clone()
            },
        ),
        (
            "coding_sequence_variant beside the gained stop",
            reference.clone(),
            Row {
                terms: "stop_gained,coding_sequence_variant",
                ..candidate.clone()
            },
        ),
        (
            "a further term the reference alone holds",
            Row {
                terms: "splice_region_variant,stop_retained_variant",
                ..reference.clone()
            },
            candidate.clone(),
        ),
        (
            "an allele of known bases, whose stop codon test is made",
            Row {
                allele: known.clone(),
                codons: Some(reference.allele.clone()),
                ..reference.clone()
            },
            Row {
                allele: known,
                codons: Some(reference.allele.clone()),
                ..candidate.clone()
            },
        ),
    ];
    for (what, reference, candidate) in off {
        assert_eq!(
            compared(&reference, &candidate, &[CLASS]),
            ("red".to_string(), nothing()),
            "{what}"
        );
    }
    // Ten positions under an allele of ten letters: no insertion's and no deletion's row.
    let (reference, candidate) = frameshift();
    let same_length = |row: &Row| {
        row.line()
            .replace("21:10567700-10567701", "21:10567700-10567709")
    };
    let (report, documented) = compare_default(
        &default_text(&[same_length(&untested(&reference))]),
        &default_text(&[same_length(&candidate)]),
        &[CLASS],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(counts(&documented.applied), nothing());
}

/// A forward-strand transcript of two exons, 1:101-130 and 1:201-248, whose coding sequence
/// runs from cDNA base 7 to cDNA base 69: the second exon opens with coding base 25, the
/// first base of codon 9.
const FORWARD: &str = "ENST00000000001";

/// The same transcript on the reverse strand, its exons 1:471-500 and 1:353-400: the second
/// exon opens at 1:400 with coding base 25.
const REVERSE: &str = "ENST00000000002";

/// A transcript of the forward-strand exons whose coding sequence opens at cDNA base 37,
/// six bases inside the second exon.
const LATE_START: &str = "ENST00000000003";

/// The forward-strand transcript under the vertebrate mitochondrial code, which reads
/// `AGA` as a stop codon and `TGA` as a tryptophan codon.
const OTHER_CODE: &str = "ENST00000000004";

/// A transcript of the forward-strand exons whose coding sequence opens at cDNA base 28,
/// with the last three bases of the first exon.
const START_AT_AN_EXON_END: &str = "ENST00000000005";

/// The made-up coding sequence of 21 codons: a start codon and seven alanine codons, which
/// the first exon holds; then an alanine, a proline, a leucine and a serine codon (codons
/// 9 to 12), eight alanine codons and a stop codon.
fn made_coding() -> String {
    format!(
        "ATG{}{}{}TAA",
        "GCT".repeat(7),
        ["GCT", "CCA", "CTT", "AGC"].concat(),
        "GCT".repeat(8)
    )
}

/// One cache record of a transcript of two exons, of 30 and of 48 bases in the
/// transcript's order, whose coding sequence `coding` opens at cDNA base `coding_start`.
fn two_exon_record(
    name: &str,
    strand: i64,
    exons: [(i64, i64); 2],
    coding_start: i64,
    coding: &str,
    table: i64,
) -> String {
    let coding_end = coding_start + i64::try_from(coding.len()).expect("a length") - 1;
    let pairs = [(1, 30, exons[0]), (31, 78, exons[1])]
        .map(|(from_start, from_end, (to_start, to_end))| {
            format!(
                r#"{{"from_start": {from_start}, "from_end": {from_end}, "to_start": {to_start}, "to_end": {to_end}, "ori": {strand}}}"#
            )
        })
        .join(", ");
    let spans = exons
        .map(|(start, end)| format!(r#"{{"start": {start}, "end": {end}}}"#))
        .join(", ");
    let cached = format!(
        r#"{{"translateable_seq": "{coding}", "codon_table": {table}, "mapper": {{"pairs": [{pairs}]}}}}"#
    );
    format!(
        r#"{{"stable_id": "{name}", "version": 1, "strand": {strand}, "biotype": "protein_coding", "exons": [{spans}], "cdna_coding_start": {coding_start}, "cdna_coding_end": {coding_end}, "variation_effect_feature_cache": {cached}}}"#
    )
}

/// A candidate cache of chromosome 1 that holds the five made-up transcripts, in a
/// directory of the test `name`.
fn carried_cache(name: &str) -> PathBuf {
    let cache = scratch(name).join("cache");
    let forward = [(101, 130), (201, 248)];
    let late = format!("ATG{}TAA", "GCT".repeat(12));
    let early = format!("ATG{}TAA", "GCT".repeat(15));
    let records = [
        two_exon_record(FORWARD, 1, forward, 7, &made_coding(), 1),
        two_exon_record(REVERSE, -1, [(471, 500), (353, 400)], 7, &made_coding(), 1),
        two_exon_record(LATE_START, 1, forward, 37, &late, 1),
        two_exon_record(OTHER_CODE, 1, forward, 7, &made_coding(), 2),
        two_exon_record(START_AT_AN_EXON_END, 1, forward, 28, &early, 1),
    ];
    let shard = cache.join("transcripts").join("1");
    std::fs::create_dir_all(&shard).expect("the chromosome's directory is made");
    planted(
        &shard,
        "1-1000000.json",
        &format!("[{}]", records.join(", ")),
    );
    cache
}

/// The side tables of a comparison that reads the transcripts of `cache`.
fn beside(cache: &Path) -> Side {
    Side {
        downstream_sequences: Some(DownstreamSequences::new(cache).expect("the cache reads")),
        ..Side::default()
    }
}

/// The terms both engines write on a row beside an intron's last bases.
const INTRON_TERMS: &str = "splice_region_variant,intron_variant";

/// One Transcript row of an insertion written in an intron, as both engines write it, and
/// the notations each engine adds to it: a coding notation, the description of a protein
/// notation and an HGVS_OFFSET, each empty for none.
#[derive(Clone, Copy)]
struct Carried<'a> {
    name: &'a str,
    location: &'a str,
    /// The inserted bases on the forward strand.
    allele: &'a str,
    feature: &'a str,
    strand: &'a str,
    terms: &'a str,
    /// The cDNA, the CDS and the protein position.
    positions: [&'a str; 3],
    codons: &'a str,
    reference: [&'a str; 3],
    candidate: [&'a str; 3],
}

/// Thirteen bases written two bases before the second exon of the forward-strand
/// transcript and carried five bases on, onto the boundary of codons 9 and 10. The allele
/// opens with the two intron bases the shift passes and with codon 9; read from the shifted
/// point, its next bases spell a glutamate codon where the proline codon stands and a stop
/// codon, before its two unknown bases.
const ON_A_BOUNDARY: Carried = Carried {
    name: "intron_insertion",
    location: "1:198-199",
    allele: "AGGCTGAATAGNN",
    feature: FORWARD,
    strand: "1",
    terms: INTRON_TERMS,
    positions: ["-"; 3],
    codons: "-",
    reference: [""; 3],
    candidate: ["", "Pro10GlufsTer2", "5"],
};

/// Thirteen bases carried six bases on, behind the first base of codon 10: that base and
/// the allele's next two spell a glutamine codon, and the three after them a stop codon.
const INSIDE_A_CODON: Carried = Carried {
    allele: "AGGCTCAATGANN",
    candidate: ["", "Pro10GlnfsTer2", "6"],
    ..ON_A_BOUNDARY
};

/// Sixteen bases carried seven bases on, behind the second base of codon 10: those two
/// bases and the allele's next one spell a proline codon again, the next three a glutamate
/// codon where the leucine codon of residue 11 stands, and the three after them a stop
/// codon.
const BEHIND_A_KEPT_RESIDUE: Carried = Carried {
    allele: "AGGCTCCGGAATAGNN",
    candidate: ["", "Leu11GlufsTer2", "7"],
    ..ON_A_BOUNDARY
};

/// The bases of `ON_A_BOUNDARY` as the reverse-strand transcript reads them, written two
/// bases before its second exon and carried five bases on in its direction: the row's
/// allele is the forward strand's, and its offset is negative.
const ON_THE_REVERSE_STRAND: Carried = Carried {
    location: "1:402-403",
    allele: "NNCTATTCAGCCT",
    feature: REVERSE,
    strand: "-1",
    candidate: ["", "Pro10GlufsTer2", "-5"],
    ..ON_A_BOUNDARY
};

/// Twelve bases carried onto the boundary of codons 9 and 10, whole codons: a glutamate
/// codon and a stop codon stand before the codon that holds the unknown base.
const IN_FRAME: Carried = Carried {
    allele: "AGGCTGAATAGN",
    candidate: ["", "Ala9_Pro10insGluTer", "5"],
    ..ON_A_BOUNDARY
};

impl Carried<'_> {
    /// The coding notation, the protein notation and the offset `notes` spell on the row,
    /// each empty for none.
    fn notations(&self, [coding, protein, offset]: [&str; 3]) -> [String; 3] {
        let named = |name: &str, description: &str| {
            if description.is_empty() {
                String::new()
            } else {
                format!("{name}{description}")
            }
        };
        [
            named(&format!("{}.1:", self.feature), coding),
            named("ENSP00000000001.1:p.", protein),
            offset.to_string(),
        ]
    }

    /// The columns of the row before the Extra text.
    fn columns(&self) -> String {
        let [cdna, cds, protein] = self.positions;
        format!(
            "{}\t{}\t{}\tENSG1\t{}\tTranscript\t{}\t{cdna}\t{cds}\t{protein}\t-\t{}\t-",
            self.name, self.location, self.allele, self.feature, self.terms, self.codons
        )
    }

    /// The row's line of a default-format output with the notations `notes`.
    fn line(&self, notes: [&str; 3]) -> String {
        let mut extra = format!("IMPACT=LOW;STRAND={}", self.strand);
        let named = ["HGVSc", "HGVSp", "HGVS_OFFSET"];
        for (name, value) in named.iter().zip(self.notations(notes)) {
            if !value.is_empty() {
                extra.push_str(&format!(";{name}={value}"));
            }
        }
        format!("{}\t{extra}", self.columns())
    }

    /// A tab-format output of the row with the notations `notes`.
    fn tab_text(&self, notes: [&str; 3]) -> String {
        let header =
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSc\tHGVSp\tHGVS_OFFSET");
        let or_dash = |value: String| {
            if value.is_empty() {
                "-".to_string()
            } else {
                value
            }
        };
        let [coding, protein, offset] = self.notations(notes).map(or_dash);
        format!(
            "{header}\n{}\tLOW\t{}\t{coding}\t{protein}\t{offset}\n",
            self.columns(),
            self.strand
        )
    }

    /// The two positions of the row's Location and its chromosome.
    fn placed(&self) -> (&str, &str, &str) {
        let (chromosome, positions) = self.location.split_once(':').expect("a Location");
        let (first, last) = positions.split_once('-').expect("two positions");
        (chromosome, first, last)
    }

    /// The input line of the row's record: a base at the first position of the Location,
    /// then that base and the inserted bases.
    fn input_line(&self) -> String {
        let (chromosome, first, _) = self.placed();
        format!("{chromosome}\t{first}\t.\tG\tG{}\t.\t.\t.", self.allele)
    }

    /// The entries of the row, which holds no position, in the tab, VCF or JSON format with
    /// the notations `notes`.
    fn in_format(&self, fmt: &str, notes: [&str; 3]) -> Vec<Entry> {
        let inputs = [self.input_line()];
        let [coding, protein, offset] = self.notations(notes);
        match fmt {
            "tab" => parse_tab(&self.tab_text(notes)).expect("the output parses"),
            "vcf" => {
                let csq = [
                    self.allele,
                    &self.terms.replace(',', "&"),
                    "LOW",
                    self.feature,
                    "Transcript",
                    "",
                    "",
                    "",
                    "",
                    "",
                    &coding,
                    &protein,
                    &offset,
                ]
                .join("|");
                let fields = "Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|HGVSc|HGVSp|HGVS_OFFSET";
                let text = format!("{}{}CSQ={csq}\n", csq_header(fields), open_info(&inputs[0]));
                parse_vcf(&text, &inputs).expect("the output parses")
            }
            _ => {
                let (chromosome, first, last) = self.placed();
                let held: Vec<&str> = self.terms.split(',').collect();
                let mut entry = format!(
                    "{{\"variant_allele\": {}, \"transcript_id\": {}, \"consequence_terms\": {}, \"impact\": \"LOW\", \"strand\": {}",
                    quoted(self.allele),
                    quoted(self.feature),
                    quoted_list(&held),
                    self.strand
                );
                for (key, value) in [("hgvsc", &coding), ("hgvsp", &protein)] {
                    if !value.is_empty() {
                        entry.push_str(&format!(", \"{key}\": {}", quoted(value)));
                    }
                }
                if !offset.is_empty() {
                    entry.push_str(&format!(", \"hgvs_offset\": {offset}"));
                }
                entry.push('}');
                let text = format!(
                    "{{\"input\": {}, \"seq_region_name\": {}, \"start\": {last}, \"end\": {first}, \"allele_string\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
                    quoted(&inputs[0]),
                    quoted(chromosome),
                    quoted(&format!("-/{}", self.allele)),
                    quoted(most_severe(&held)),
                );
                parse_json(&text, &inputs).expect("the output parses")
            }
        }
    }
}

/// The verdict of two entry lists under `classes` beside `side`, and what was applied.
fn read_beside(
    reference: &[Entry],
    candidate: &[Entry],
    classes: &[&str],
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare_entries(
        reference,
        candidate,
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The verdict of the default-format row of `reference` with the reference's notations
/// against that of `candidate` with the candidate's.
fn carried_pair(
    reference: &Carried,
    candidate: &Carried,
    classes: &[&str],
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let entries = |line: String| parse_default(&default_text(&[line])).expect("the row parses");
    read_beside(
        &entries(reference.line(reference.reference)),
        &entries(candidate.line(candidate.candidate)),
        classes,
        side,
    )
}

/// The verdict of the two default-format rows of `row`.
fn carried(row: &Carried, classes: &[&str], side: &Side) -> (String, BTreeMap<String, i64>) {
    carried_pair(row, row, classes, side)
}

fn set_aside() -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), count_map(&[(CLASS, 1)]))
}

fn left() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), nothing())
}

/// A row the class takes beside the cache: its protein notation and its offset are set
/// aside together, with the description in a pair of prediction parentheses as well.
/// Without the class named, and with its sides exchanged, the row is a difference.
fn assert_set_aside(row: &Carried, side: &Side) {
    let description = row.candidate[1];
    assert_eq!(carried(row, &[CLASS], side), set_aside(), "{description}");
    let predicted = format!("({description})");
    let in_parentheses = Carried {
        candidate: [row.candidate[0], &predicted, row.candidate[2]],
        ..*row
    };
    assert_eq!(
        carried(&in_parentheses, &[CLASS], side),
        set_aside(),
        "{predicted}"
    );
    assert_eq!(carried(row, &[], side), left(), "{description} unnamed");
    let exchanged = Carried {
        reference: row.candidate,
        candidate: row.reference,
        ..*row
    };
    assert_eq!(
        carried(&exchanged, &[CLASS], side),
        left(),
        "{description} exchanged"
    );
}

/// Each named row stays a difference beside the cache, with no cell set aside.
fn assert_left(rows: &[(&str, Carried)], side: &Side) {
    for (what, row) in rows {
        assert_eq!(carried(row, &[CLASS], side), left(), "{what}");
    }
}

/// The shifted point stands on a codon boundary of the forward-strand transcript. The codon
/// that opens there reads a glutamate for the proline of residue 10 and the next one is a
/// stop codon, both of known bases, so the description `Pro10GlufsTer2` is the one every
/// completion of the unknown bases shares. A row whose first column spells the record's
/// alleles is read the same way.
#[test]
fn an_insertion_carried_onto_a_codon_boundary_is_set_aside_beside_the_cache() {
    let side = beside(&carried_cache("carried-codon-boundary"));
    assert_set_aside(&ON_A_BOUNDARY, &side);
    let spelt = Carried {
        name: "1_198_G/GAGGCTGAATAGNN",
        ..ON_A_BOUNDARY
    };
    assert_set_aside(&spelt, &side);
}

/// The shifted point stands inside codon 10, which opens with its reference bases: behind
/// its first base the changed codon is the first changed residue, and behind its second
/// base, where the changed codon reads the reference's proline again, the first changed
/// residue is the next one and the count runs from it.
#[test]
fn an_insertion_carried_inside_a_codon_is_set_aside_beside_the_cache() {
    let side = beside(&carried_cache("carried-inside-a-codon"));
    assert_set_aside(&INSIDE_A_CODON, &side);
    assert_set_aside(&BEHIND_A_KEPT_RESIDUE, &side);
}

/// On the reverse-strand transcript the row's allele is the other strand of the bases the
/// transcript reads, the shift runs towards lower positions and the offset is negative.
#[test]
fn an_insertion_carried_on_the_reverse_strand_is_set_aside_beside_the_cache() {
    let side = beside(&carried_cache("carried-reverse-strand"));
    assert_set_aside(&ON_THE_REVERSE_STRAND, &side);
}

/// An insertion of whole codons whose notation names the inserted stop is no frameshift:
/// the class takes it where a stop codon stands before the first unknown codon.
#[test]
fn an_inframe_insertion_that_names_an_inserted_stop_is_set_aside_beside_the_cache() {
    let side = beside(&carried_cache("carried-in-frame"));
    assert_set_aside(&IN_FRAME, &side);
}

/// In the tab, VCF and JSON formats the row is read as in the default format, on either
/// strand: a VCF or JSON entry is placed by the Location its record gives its allele.
#[test]
fn a_carried_insertion_is_set_aside_in_the_tab_vcf_and_json_formats() {
    let side = beside(&carried_cache("carried-formats"));
    for row in [ON_A_BOUNDARY, ON_THE_REVERSE_STRAND] {
        for fmt in ["tab", "vcf", "json"] {
            let label = format!("{fmt} strand {}", row.strand);
            let reference = row.in_format(fmt, row.reference);
            let candidate = row.in_format(fmt, row.candidate);
            assert_eq!(
                read_beside(&reference, &candidate, &[CLASS], &side),
                set_aside(),
                "{label}"
            );
            assert_eq!(
                read_beside(&reference, &candidate, &[CLASS], &Side::default()),
                left(),
                "{label} without the cache"
            );
            assert_eq!(
                read_beside(&reference, &candidate, &[], &side),
                left(),
                "{label} unnamed"
            );
            assert_eq!(
                read_beside(&candidate, &reference, &[CLASS], &side),
                left(),
                "{label} exchanged"
            );
        }
    }
}

/// The command line reads the candidate's cache for the class named alone: two tab-format
/// outputs of the row compare red without `--json-cache`, nothing applied, and green beside
/// it, the class applied once.
#[test]
fn a_comparison_reads_the_cache_the_command_line_names_for_a_carried_insertion() {
    let cache = carried_cache("carried-command-line");
    let directory = cache
        .parent()
        .expect("the cache stands in a directory")
        .to_path_buf();
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let output = |name: &str, notes: [&str; 3]| {
        text(&planted(&directory, name, &ON_A_BOUNDARY.tab_text(notes)))
    };
    let reference = output("reference.tab", ON_A_BOUNDARY.reference);
    let candidate = output("candidate.tab", ON_A_BOUNDARY.candidate);
    let compared = |more: &[&str], name: &str| {
        let path = directory.join(name);
        let mut argv = args_of(&[
            "--format",
            "tab",
            "--expected",
            &reference,
            "--actual",
            &candidate,
            "--masks",
            CLASS,
            "--report",
            &text(&path),
        ]);
        argv.extend(args_of(more));
        (run(&directory, &argv), path)
    };
    let (ran, path) = compared(&[], "without.json");
    assert_eq!(ran.code, 1, "{}", ran.err());
    assert_eq!(compact(member(&report(&path), &["masks", "applied"])), "{}");
    let (ran, path) = compared(&["--json-cache", &text(&cache)], "beside.json");
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        compact(member(&report(&path), &["masks", "applied"])),
        format!(r#"{{"{CLASS}": 1}}"#)
    );
}

/// Without the candidate's cache the class reads no row that lacks a `Codons` cell.
#[test]
fn without_the_cache_a_carried_insertion_stays_a_difference() {
    let side = beside(&carried_cache("carried-without-the-cache"));
    for row in [
        ON_A_BOUNDARY,
        INSIDE_A_CODON,
        ON_THE_REVERSE_STRAND,
        IN_FRAME,
    ] {
        assert_eq!(carried(&row, &[CLASS], &side), set_aside());
        assert_eq!(
            carried(&row, &[CLASS], &Side::default()),
            left(),
            "{}",
            row.candidate[1]
        );
    }
}

/// The codons are read in the transcript's frame from the codon that holds the shifted
/// point. An allele whose second codon there holds the unknown bases settles no stop,
/// although the frame one base on reads a stop codon of known bases first and the
/// transcript's frame reads one after the unknown codon.
#[test]
fn a_stop_codon_behind_the_first_unknown_codon_stays_a_difference() {
    let side = beside(&carried_cache("carried-stop-behind-unknown"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    assert_left(
        &[(
            "a stop codon in another frame and one behind the unknown codon",
            Carried {
                allele: "AGGCTGTAGNNTAA",
                candidate: ["", "Pro10ValfsTer3", "5"],
                ..ON_A_BOUNDARY
            },
        )],
        &side,
    );
}

/// The shift stops at a known base: an offset as large as the count of known bases the
/// allele opens with, or larger, is one the bases do not give, although under the larger
/// one the bases passed in the exon are the exon's own and the allele turned by it opens
/// with a stop codon.
#[test]
fn a_shift_that_reaches_an_unknown_base_stays_a_difference() {
    let side = beside(&carried_cache("carried-shift-reaches-unknown"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    assert_left(
        &[
            (
                "five known bases open the allele and the offset is five",
                Carried {
                    allele: "AGGCTNTAGTA",
                    candidate: ["", "Pro10Ter", "5"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "one known base opens the allele and the offset is five",
                Carried {
                    allele: "ANGCTTAGTA",
                    candidate: ["", "Pro10Ter", "5"],
                    ..ON_A_BOUNDARY
                },
            ),
        ],
        &side,
    );
}

/// The reference writes no notation on a row of the class and the candidate no coding
/// notation: a protein notation or an offset of the reference's own, and a coding notation
/// on either side, leave the row a difference.
#[test]
fn a_notation_of_the_reference_or_a_coding_notation_stays_a_difference() {
    let side = beside(&carried_cache("carried-other-notations"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    let coding = "c.27_28insGAATAG";
    assert_left(
        &[
            (
                "another protein notation in the reference",
                Carried {
                    reference: ["", "Pro10GlufsTer9", ""],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "the protein notation in the reference, the offset in the candidate alone",
                Carried {
                    reference: ["", "Pro10GlufsTer2", ""],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "the offset in the reference, the protein notation in the candidate alone",
                Carried {
                    reference: ["", "", "5"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "another offset in the reference",
                Carried {
                    reference: ["", "", "4"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a coding notation the candidate alone writes",
                Carried {
                    candidate: [coding, "Pro10GlufsTer2", "5"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a coding notation both write",
                Carried {
                    reference: [coding, "", ""],
                    candidate: [coding, "Pro10GlufsTer2", "5"],
                    ..ON_A_BOUNDARY
                },
            ),
        ],
        &side,
    );
}

/// The row of the class is annotated in an intron and its allele holds an unknown base: an
/// allele of known bases alone, a row that carries a cDNA, a CDS or a protein position and
/// a row with a `Codons` cell, which the class reads by that cell alone, stay differences.
#[test]
fn an_allele_of_known_bases_or_a_row_with_a_coding_cell_stays_a_difference() {
    let side = beside(&carried_cache("carried-coding-cells"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    assert_left(
        &[
            (
                "an allele of known bases alone",
                Carried {
                    allele: "AGGCTGAATAGAA",
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a cDNA position",
                Carried {
                    positions: ["33-34", "-", "-"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a CDS position",
                Carried {
                    positions: ["-", "27-28", "-"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a protein position",
                Carried {
                    positions: ["-", "-", "9-10"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a Codons cell",
                Carried {
                    codons: "-/GAATAGNNAGGCT",
                    ..ON_A_BOUNDARY
                },
            ),
        ],
        &side,
    );
}

/// The insertion is written outside every exon and carried to a point with a coding base
/// of one exon on either side. A point still in the intron, one in the exon's bases before
/// the coding sequence, one between the intron's last base and the exon's first, an
/// insertion written inside the exon, among its coding bases or among its bases before the
/// coding sequence, and a transcript the cache does not hold stay differences, the last
/// four under the description the same reading gives a point it takes.
#[test]
fn a_shifted_point_outside_one_coding_exon_stays_a_difference() {
    let side = beside(&carried_cache("carried-placement"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    assert_left(
        &[
            (
                "the shifted point lies in the intron",
                Carried {
                    location: "1:150-151",
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "the shifted point lies before the coding sequence",
                Carried {
                    feature: LATE_START,
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "the shifted point lies beside the exon's first base",
                Carried {
                    candidate: ["", "Pro10GlufsTer2", "2"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "the insertion is written inside the exon",
                Carried {
                    location: "1:201-202",
                    allele: "CTGAATAGNNAGG",
                    candidate: ["", "Pro10GlufsTer2", "2"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "the insertion is written inside the exon, before the coding sequence",
                Carried {
                    location: "1:204-205",
                    allele: "CCATGTAGNN",
                    feature: LATE_START,
                    candidate: ["", "Ala2Ter", "5"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a transcript the cache does not hold",
                Carried {
                    feature: "ENST00000000009",
                    ..ON_A_BOUNDARY
                },
            ),
        ],
        &side,
    );
}

/// The row `row` with the candidate's protein notation described as `description`.
fn described<'a>(row: Carried<'a>, description: &'a str) -> Carried<'a> {
    Carried {
        candidate: ["", description, row.candidate[2]],
        ..row
    }
}

/// A frameshift's description is held to the sequence: the number of the first codon from
/// the point's codon on whose residue differs from the reference's, and the count from it
/// to the stop codon, the first counted. Another count, an unknown count, another residue
/// number beside the residue names the sequence gives, and a count from the point's codon
/// where that codon keeps its residue stay differences.
#[test]
fn a_frameshift_description_the_sequence_does_not_give_stays_a_difference() {
    let side = beside(&carried_cache("carried-frameshift-description"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    assert_eq!(
        carried(&BEHIND_A_KEPT_RESIDUE, &[CLASS], &side),
        set_aside()
    );
    assert_left(
        &[
            (
                "a count one too large",
                described(ON_A_BOUNDARY, "Pro10GlufsTer3"),
            ),
            (
                "a count one too small",
                described(ON_A_BOUNDARY, "Pro10GlufsTer1"),
            ),
            (
                "an unknown count",
                described(ON_A_BOUNDARY, "Pro10GlufsTer?"),
            ),
            (
                "the residue number before the first changed one",
                described(ON_A_BOUNDARY, "Pro9GlufsTer2"),
            ),
            (
                "the residue number behind the first changed one",
                described(ON_A_BOUNDARY, "Pro11GlufsTer2"),
            ),
            (
                "the point's codon named where it keeps its residue",
                described(BEHIND_A_KEPT_RESIDUE, "Pro10GlufsTer3"),
            ),
            (
                "a count from the point's codon behind the first changed residue",
                described(BEHIND_A_KEPT_RESIDUE, "Leu11GlufsTer3"),
            ),
        ],
        &side,
    );
}

/// The offset's sign is the transcript's strand. An offset of the other sign stays a
/// difference, also where the point it gives on the genome lies in the coding sequence of
/// the exon before the intron and the allele read there settles a stop codon, and where
/// the bases it passes there, the first two of a coding sequence that opens with the exon's
/// last three bases, are the allele's own and the third is another than the allele's next.
#[test]
fn an_offset_whose_sign_is_against_the_strand_stays_a_difference() {
    let side = beside(&carried_cache("carried-offset-sign"));
    assert_eq!(carried(&ON_A_BOUNDARY, &[CLASS], &side), set_aside());
    assert_eq!(
        carried(&ON_THE_REVERSE_STRAND, &[CLASS], &side),
        set_aside()
    );
    assert_left(
        &[
            (
                "a negative offset on the forward strand",
                Carried {
                    candidate: ["", "Pro10GlufsTer2", "-5"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a positive offset on the reverse strand",
                Carried {
                    candidate: ["", "Pro10GlufsTer2", "5"],
                    ..ON_THE_REVERSE_STRAND
                },
            ),
            (
                "a negative offset into the first exon of the forward strand",
                Carried {
                    location: "1:132-133",
                    candidate: ["", "Ala8GlufsTer2", "-5"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "a positive offset into the first exon of the reverse strand",
                Carried {
                    location: "1:468-469",
                    candidate: ["", "Ala8GlufsTer2", "5"],
                    ..ON_THE_REVERSE_STRAND
                },
            ),
            (
                "a negative offset onto the start codon that ends the first exon",
                Carried {
                    location: "1:132-133",
                    allele: "CATATAGNN",
                    feature: START_AT_AN_EXON_END,
                    candidate: ["", "Met1delinsIleTer", "-3"],
                    ..ON_A_BOUNDARY
                },
            ),
        ],
        &side,
    );
}

/// The offset is held to the exon bases the shift passes: the coding bases of the point's
/// exon that stand before the point, counted from the exon's first base or, in an exon
/// that opens before the coding sequence, from the first coding base, are the last of the
/// allele's bases the shift passes. An offset one base longer than the bases give stays a
/// difference under the description the allele turned by it would have, on either strand,
/// and so does an allele whose passed bases end in another base than the start codon's
/// last. An offset one base shorter than the bases give, which leaves behind the point a
/// base the allele's next base equals, stays a difference under the description both
/// offsets have, on either strand and inside the start codon.
#[test]
fn an_offset_the_exon_bases_do_not_give_stays_a_difference() {
    let side = beside(&carried_cache("carried-exon-bases"));
    // The shift passes two bases of the intron, the exon's six bases before the coding
    // sequence and the start codon; the allele's next three bases spell a stop codon in
    // the place of codon 2.
    let behind_a_late_start = Carried {
        allele: "CCGGTTCCATGTAGNN",
        feature: LATE_START,
        candidate: ["", "Ala2Ter", "11"],
        ..ON_A_BOUNDARY
    };
    // The bases of `INSIDE_A_CODON` as the reverse-strand transcript reads them.
    let inside_a_codon_on_the_reverse_strand = Carried {
        allele: "NNTCATTGAGCCT",
        candidate: ["", "Pro10GlnfsTer2", "-6"],
        ..ON_THE_REVERSE_STRAND
    };
    for row in [
        ON_A_BOUNDARY,
        ON_THE_REVERSE_STRAND,
        behind_a_late_start,
        INSIDE_A_CODON,
        inside_a_codon_on_the_reverse_strand,
    ] {
        assert_eq!(
            carried(&row, &[CLASS], &side),
            set_aside(),
            "{} on {}",
            row.candidate[1],
            row.feature
        );
    }
    assert_left(
        &[
            (
                "six bases passed where the allele's sixth base is not the exon's fourth",
                Carried {
                    candidate: ["", "Pro10GlnfsTer2", "6"],
                    ..ON_A_BOUNDARY
                },
            ),
            (
                "six bases passed on the reverse strand",
                Carried {
                    candidate: ["", "Pro10GlnfsTer2", "-6"],
                    ..ON_THE_REVERSE_STRAND
                },
            ),
            (
                "an eleventh base that is not the start codon's last",
                Carried {
                    allele: "CCGGTTCCATCTAGNN",
                    ..behind_a_late_start
                },
            ),
            (
                "five bases passed where the allele's sixth base is the exon's fourth",
                Carried {
                    candidate: ["", "Pro10GlnfsTer2", "5"],
                    ..INSIDE_A_CODON
                },
            ),
            (
                "five bases passed on the reverse strand where the sixth is the exon's fourth",
                Carried {
                    candidate: ["", "Pro10GlnfsTer2", "-5"],
                    ..inside_a_codon_on_the_reverse_strand
                },
            ),
            (
                "ten bases passed where the allele's eleventh base is the start codon's last",
                Carried {
                    candidate: ["", "Ala2Ter", "10"],
                    ..behind_a_late_start
                },
            ),
        ],
        &side,
    );
}

/// A frameshift's description names the reference's residue at the first changed codon and
/// the residue the changed codon reads there. Another name in either place stays a
/// difference under the number and the count the sequence gives, on either strand.
#[test]
fn a_frameshift_description_with_another_residue_name_stays_a_difference() {
    let side = beside(&carried_cache("carried-frameshift-residue-names"));
    for row in [ON_A_BOUNDARY, BEHIND_A_KEPT_RESIDUE, ON_THE_REVERSE_STRAND] {
        assert_eq!(
            carried(&row, &[CLASS], &side),
            set_aside(),
            "{} on {}",
            row.candidate[1],
            row.feature
        );
    }
    assert_left(
        &[
            (
                "an alanine named where the changed codon reads a glutamate",
                described(ON_A_BOUNDARY, "Pro10AlafsTer2"),
            ),
            (
                "an alanine named where the reference reads a proline",
                described(ON_A_BOUNDARY, "Ala10GlufsTer2"),
            ),
            (
                "the proline the point's codon keeps named for the leucine that changes",
                described(BEHIND_A_KEPT_RESIDUE, "Pro11GlufsTer2"),
            ),
            (
                "a glutamine named on the reverse strand where the changed codon reads a glutamate",
                described(ON_THE_REVERSE_STRAND, "Pro10GlnfsTer2"),
            ),
        ],
        &side,
    );
}

/// The stop codon is one of the transcript's own genetic code. Under the vertebrate
/// mitochondrial code `AGA` settles the stop and `TGA` does not; under the standard code
/// the same alleles read the other way.
#[test]
fn a_stop_codon_is_read_under_the_genetic_code_of_the_transcript() {
    let side = beside(&carried_cache("carried-genetic-code"));
    let arginine_codon = Carried {
        allele: "AGGCTGAAAGANN",
        ..ON_A_BOUNDARY
    };
    assert_set_aside(&INSIDE_A_CODON, &side);
    assert_set_aside(
        &Carried {
            feature: OTHER_CODE,
            ..arginine_codon
        },
        &side,
    );
    assert_left(
        &[
            (
                "TGA under the vertebrate mitochondrial code",
                Carried {
                    feature: OTHER_CODE,
                    ..INSIDE_A_CODON
                },
            ),
            ("AGA under the standard code", arginine_codon),
        ],
        &side,
    );
}

/// A row off the shape stays a difference: a protein notation that names no stop, an
/// offset beside no protein notation, an offset of zero, a Location of three positions, a
/// record that replaces the two bases of the Location, a term that places the row outside
/// the transcript, and two sets of terms that differ.
#[test]
fn a_carried_row_off_the_shape_stays_a_difference() {
    let side = beside(&carried_cache("carried-off-the-shape"));
    assert_eq!(carried(&IN_FRAME, &[CLASS], &side), set_aside());
    assert_left(
        &[
            (
                "a protein notation that names no stop",
                Carried {
                    candidate: ["", "Ala9_Pro10insGlu", "5"],
                    ..IN_FRAME
                },
            ),
            (
                "an offset beside no protein notation",
                Carried {
                    candidate: ["", "", "5"],
                    ..IN_FRAME
                },
            ),
            (
                "an offset of zero",
                Carried {
                    candidate: ["", "Ala9_Pro10insGluTer", "0"],
                    ..IN_FRAME
                },
            ),
            (
                "a Location of three positions",
                Carried {
                    location: "1:197-199",
                    ..IN_FRAME
                },
            ),
            (
                "a record whose reference allele keeps two bases",
                Carried {
                    name: "1_197_CGA/CAGGCTGAATAGN",
                    ..IN_FRAME
                },
            ),
            (
                "a term that places the row outside the transcript",
                Carried {
                    terms: "upstream_gene_variant",
                    ..IN_FRAME
                },
            ),
        ],
        &side,
    );
    let other_terms = Carried {
        terms: "intron_variant",
        ..IN_FRAME
    };
    assert_eq!(
        carried_pair(&other_terms, &IN_FRAME, &[CLASS], &side),
        left(),
        "two sets of terms that differ"
    );
}

/// The `--masks` help states the row the class reads beside the cache, and the
/// `--json-cache` help what the class reads from the cache and what it takes without it.
#[test]
fn the_help_states_what_the_class_reads_beside_the_cache() {
    let help = include_str!("../data/help.txt");
    let help = help.split_whitespace().collect::<Vec<_>>().join(" ");
    assert!(help.contains(
        "most_severe_consequence follow; with --json-cache also the protein notation that \
         names a stop and the HGVS_OFFSET beside it, which the candidate alone writes on the \
         row of an insertion written in an intron, a row whose two sets of terms are equal \
         and that holds no position and no Codons cell, where the cache's transcript places \
         the insertion, moved by the offset, between two coding bases of one exon and the \
         allele's known bases settle a stop codon there before its first unknown base),"
    ));
    assert!(help.contains(&format!(
        "{CLASS} reads its transcripts' strands, exons, coding bounds, mapper pairs, coding \
         sequences and codon tables for the row of an insertion written in an intron and \
         carried into the coding sequence (without it the class takes its rows with a \
         Codons cell alone)"
    )));
}

/// Asked directly, the class's rule takes the two cells of the row beside the transcripts
/// of the cache and neither cell without them; beside them a row whose two sets of terms
/// differ is none of its rows.
#[test]
fn the_rule_asked_directly_reads_a_carried_insertion_beside_the_transcripts_alone() {
    let side = beside(&carried_cache("carried-asked-directly"));
    let entry = |notes: [&str; 3]| {
        parse_default(&default_text(&[ON_A_BOUNDARY.line(notes)]))
            .expect("the row parses")
            .into_iter()
            .next()
            .expect("the output holds the row")
    };
    let e = entry(ON_A_BOUNDARY.reference);
    let a = entry(ON_A_BOUNDARY.candidate);
    let asked = |e: &Entry, name: &str, side: &Side| {
        let context = RuleContext {
            record_rows: None,
            candidate_rows: None,
            declared: &[CLASS],
            assembly: None,
            dbnsfp_columns: &[],
            side,
        };
        let value = a.fields.get(name).expect("the candidate writes the cell");
        field_class_accepts(CLASS, e, name, "-", value, Some(&a), &context)
            .expect("the rule answers")
    };
    let other_terms = Entry {
        consequence: "intron_variant".to_string(),
        ..e.clone()
    };
    for name in ["HGVSp", "HGVS_OFFSET"] {
        assert!(asked(&e, name, &side), "{name}");
        assert!(
            !asked(&e, name, &Side::default()),
            "{name} without the transcripts"
        );
        assert!(
            !asked(&other_terms, name, &side),
            "{name} beside other terms"
        );
    }
}
