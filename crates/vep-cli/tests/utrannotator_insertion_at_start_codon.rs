// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--plugin UTRAnnotator` through the `vep` binary on insertions around a
//! start codon of a 5' UTR. An insertion that only moves the codon on leaves
//! it whole, as many bases further from the cap as the allele is long, with
//! every base after it unchanged: the insertion directly before the codon's
//! first base, and the insertion after its first or second base whose
//! inserted bases end with that base or those two, which is the same mutated
//! sequence. The row reports no start codon gained or lost: nothing where no
//! other analysis fires, the frameshift where the insertion lies inside
//! another upstream ORF, and the same two fields whichever way the record
//! writes it. An upstream ORF whose stop codon begins at the moved codon's
//! second base is left out of the frameshift when the mutated sequence still
//! reads a stop codon at that index; when none stands there, the length is
//! no multiple of three and the inserted bases can be written inside that
//! stop codon, the row reports its lost stop. A start codon the inserted
//! bases make is reported gained. An
//! insertion after the codon's first or second base that opens with the
//! codon's remaining bases, and does not also end with the codon's bases
//! before it, leaves the codon at its own index, the inserted bases directly
//! after it: the row reports no gained start for that codon, and the
//! frameshift analysis reads the insertion directly after the codon; the
//! other analyses read it where the record writes it. Inserted bases that do
//! both repeat the codon: the insertion is read as moving the codon on, and
//! the row reports a gained start at the codon's own index. After any other
//! insertion that follows the codon's first or second base the lost-start
//! test fires where the codon's old index reads no start codon.
//! Ensembl VEP reports the moved codon lost: it passes an insertion at the
//! codon's three indices (UTRAnnotator.pm:764-769) and reads the three bases
//! at the codon's old index alone (UTRAnnotator.pm:776-778).

use std::path::{Path, PathBuf};

mod common;

const GAINED_START: &str = "5_prime_UTR_premature_start_codon_gain_variant";
const LOST_START: &str = "5_prime_UTR_premature_start_codon_loss_variant";
const FRAMESHIFT: &str = "5_prime_UTR_uORF_frameshift_variant";
const LOST_STOP: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";

/// The corpus whose cache carries ENST00000427445 (forward strand; a start
/// codon at 5' UTR indices 145 and 281, cDNA 146 and 282) and ENST00000400211
/// (reverse strand; a start codon at index 48, cDNA 49), with the reference
/// and the uORF evidence table the plugin reads. ENST00001091162 (forward
/// strand; a start codon at index 360) and ENST00000318948 (reverse strand; a
/// start codon at index 493) each hold an upstream ORF whose stop codon
/// begins at that start codon's second base.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-plugin-utr-annotator")
}

/// Runs the binary on `records` with the plugin and returns, per record
/// name, the `5UTR_consequence` and `5UTR_annotation` of its row on `feature`.
fn annotate(records: &[&str], feature: &str) -> Vec<(String, String, String)> {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push_str("\t.\t.\t.\n");
    }
    std::fs::write(&input, text).unwrap();
    let output = tmp.path().join("out.txt");
    let corpus = corpus();
    let run = common::isolated_vep_anywhere()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--plugin")
        .arg(format!(
            "UTRAnnotator,file={}",
            corpus
                .join("data/uORF_5UTR_GRCh38_PUBLIC.chr21.txt")
                .display()
        ))
        .output()
        .expect("vep binary runs");
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let field = |extra: &str, key: &str| {
        extra
            .split(';')
            .find_map(|pair| pair.strip_prefix(key).map(str::to_owned))
            .unwrap_or_default()
    };
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|columns| columns[4] == feature)
        .map(|columns| {
            (
                columns[0].to_owned(),
                field(columns[13], "5UTR_consequence="),
                field(columns[13], "5UTR_annotation="),
            )
        })
        .collect()
}

/// The two fields of the row of `name` among `rows`.
fn row<'a>(rows: &'a [(String, String, String)], name: &str) -> (&'a str, &'a str) {
    let (_, consequence, annotation) = rows
        .iter()
        .find(|(record, _, _)| record == name)
        .unwrap_or_else(|| panic!("no row of {name}"));
    (consequence, annotation)
}

/// Three bases inserted directly before a start codon on which no other
/// analysis fires: the two fields report nothing on either strand.
#[test]
fn insertion_directly_before_a_start_codon_reports_no_lost_start() {
    let forward = annotate(&["21\t10524582\tbefore_fwd\tA\tACCC"], "ENST00000427445");
    assert_eq!(row(&forward, "before_fwd"), ("-", "-"), "forward strand");
    let reverse = annotate(&["21\t14210907\tbefore_rev\tT\tTGGG"], "ENST00000400211");
    assert_eq!(row(&reverse, "before_rev"), ("-", "-"), "reverse strand");
}

/// One base inserted directly before a start codon that lies inside another
/// upstream ORF: the row reports that ORF's frameshift, on either strand.
#[test]
fn insertion_directly_before_a_start_codon_inside_an_orf_reports_the_frameshift() {
    let forward = annotate(&["21\t10527384\tshift_fwd\tT\tTC"], "ENST00000427445");
    let (consequence, annotation) = row(&forward, "shift_fwd");
    assert_eq!(consequence, FRAMESHIFT, "forward strand: {annotation}");
    assert!(
        !annotation.contains("CapDistanceToStart="),
        "forward strand: {annotation}"
    );
    let reverse = annotate(&["21\t14210907\tshift_rev\tT\tTG"], "ENST00000400211");
    let (consequence, annotation) = row(&reverse, "shift_rev");
    assert_eq!(consequence, FRAMESHIFT, "reverse strand: {annotation}");
    assert!(
        !annotation.contains("CapDistanceToStart="),
        "reverse strand: {annotation}"
    );
}

/// Three bases inserted after a start codon's first base and after its
/// second that neither end with the codon's first bases nor open with its
/// last part the codon, which the row reports lost with the codon's own
/// index; inserted after its third base they touch no start codon.
#[test]
fn insertion_inside_a_start_codon_reports_it_lost() {
    let forward = annotate(
        &[
            "21\t10524583\tafter_first_fwd\tA\tACCC",
            "21\t10524584\tafter_second_fwd\tT\tTCCC",
            "21\t10524585\tafter_third_fwd\tG\tGCCC",
        ],
        "ENST00000427445",
    );
    for name in ["after_first_fwd", "after_second_fwd"] {
        let (consequence, annotation) = row(&forward, name);
        assert_eq!(consequence, LOST_START, "{name}");
        assert!(
            annotation.starts_with("CapDistanceToStart=145:"),
            "{name}: {annotation}"
        );
    }
    assert_eq!(row(&forward, "after_third_fwd"), ("-", "-"));
    let reverse = annotate(
        &[
            "21\t14210904\tafter_third_rev\tT\tTGGG",
            "21\t14210905\tafter_second_rev\tC\tCGGG",
            "21\t14210906\tafter_first_rev\tA\tAGGG",
        ],
        "ENST00000400211",
    );
    for name in ["after_first_rev", "after_second_rev"] {
        let (consequence, annotation) = row(&reverse, name);
        assert_eq!(consequence, LOST_START, "{name}");
        assert!(
            annotation.starts_with("CapDistanceToStart=48:"),
            "{name}: {annotation}"
        );
    }
    assert_eq!(row(&reverse, "after_third_rev"), ("-", "-"));
}

/// Three bases inserted after a start codon's first base that end with that
/// base, or after its second base that end with its first two, leave the
/// codon whole three bases on. No other analysis fires on these codons, so
/// the two fields report nothing, as they do for the same mutated sequence
/// written directly before the codon.
#[test]
fn insertion_inside_a_start_codon_that_only_moves_it_on_reports_nothing() {
    let forward = annotate(
        &[
            "21\t10524582\tfirst_fwd_before\tA\tAACC",
            "21\t10524582\tsecond_fwd_before\tA\tAATC",
            "21\t10524583\tfirst_fwd\tA\tACCA",
            "21\t10524584\tsecond_fwd\tT\tTCAT",
        ],
        "ENST00000427445",
    );
    for name in ["first_fwd", "second_fwd"] {
        assert_eq!(row(&forward, name), ("-", "-"), "{name}");
        assert_eq!(
            row(&forward, name),
            row(&forward, &format!("{name}_before")),
            "{name}"
        );
    }
    let reverse = annotate(
        &[
            "21\t14210905\tsecond_rev\tC\tCATG",
            "21\t14210906\tfirst_rev\tA\tATGG",
            "21\t14210907\tfirst_rev_before\tT\tTGGT",
            "21\t14210907\tsecond_rev_before\tT\tTGAT",
        ],
        "ENST00000400211",
    );
    for name in ["first_rev", "second_rev"] {
        assert_eq!(row(&reverse, name), ("-", "-"), "{name}");
        assert_eq!(
            row(&reverse, name),
            row(&reverse, &format!("{name}_before")),
            "{name}"
        );
    }
}

/// One base that repeats a start codon's first base, inserted after it, and
/// two that repeat its first two, inserted after its second, move the codon
/// on inside another upstream ORF: the row reports that ORF's frameshift,
/// field for field the row of the same bases written directly before the
/// codon.
#[test]
fn insertion_that_moves_a_start_codon_on_inside_an_orf_reports_the_frameshift() {
    let forward = annotate(
        &[
            "21\t10527384\tfirst_fwd_before\tT\tTA",
            "21\t10527384\tsecond_fwd_before\tT\tTAT",
            "21\t10527385\tfirst_fwd\tA\tAA",
            "21\t10527386\tsecond_fwd\tT\tTAT",
        ],
        "ENST00000427445",
    );
    let reverse = annotate(
        &[
            "21\t14210905\tsecond_rev\tC\tCAT",
            "21\t14210906\tfirst_rev\tA\tAT",
            "21\t14210907\tfirst_rev_before\tT\tTT",
            "21\t14210907\tsecond_rev_before\tT\tTAT",
        ],
        "ENST00000400211",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        for name in [format!("first_{strand}"), format!("second_{strand}")] {
            let (consequence, annotation) = row(rows, &name);
            assert_eq!(consequence, FRAMESHIFT, "{name}: {annotation}");
            assert_eq!(
                row(rows, &name),
                row(rows, &format!("{name}_before")),
                "{name}"
            );
        }
    }
}

/// An upstream ORF's stop codon begins at the second base of the start
/// codon these insertions move on, so the inserted bases stand before that
/// stop codon wherever the record writes them. One base cannot be written
/// inside the stop codon: the row reports that ORF's frameshift. Two bases
/// written after the codon's second base stand inside the stop codon, whose
/// index reads no stop codon: the row reports that ORF's lost stop. Each is
/// field for field the row of the same mutated sequence written directly
/// before the codon.
#[test]
fn insertion_that_moves_a_start_codon_on_shifts_the_orf_that_ends_inside_it() {
    let forward = annotate(
        &[
            "21\t10535884\tfirst_fwd_before\tC\tCA",
            "21\t10535884\tsecond_fwd_before\tC\tCAT",
            "21\t10535885\tfirst_fwd\tA\tAA",
            "21\t10535886\tsecond_fwd\tT\tTAT",
        ],
        "ENST00001091162",
    );
    let reverse = annotate(
        &[
            "21\t14968490\tsecond_rev\tC\tCAT",
            "21\t14968491\tfirst_rev\tA\tAT",
            "21\t14968492\tfirst_rev_before\tT\tTT",
            "21\t14968492\tsecond_rev_before\tT\tTAT",
        ],
        "ENST00000318948",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        for (name, term) in [
            (format!("first_{strand}"), FRAMESHIFT),
            (format!("second_{strand}"), LOST_STOP),
        ] {
            let (consequence, annotation) = row(rows, &name);
            assert_eq!(consequence, term, "{name}: {annotation}");
            assert_eq!(
                row(rows, &name),
                row(rows, &format!("{name}_before")),
                "{name}"
            );
        }
    }
}

/// Eight bases written after the second base of a start codon they move on
/// that complete no stop codon with that base and hold one a codon further
/// on: the stop codon of the upstream ORF that begins at the codon's second
/// base is lost, and the row reports it with the stop codon the inserted
/// bases hold, the same two fields on the three ways of writing the
/// sequence. Ensembl VEP reads the stop codon's index only for the record
/// written inside it (UTRAnnotator.pm:628-633).
#[test]
fn insertion_written_inside_the_stop_codon_at_a_moved_start_codon_reports_the_lost_stop() {
    let forward = annotate(
        &[
            "21\t10535884\tlost_fwd_before\tC\tCATCCTAGG",
            "21\t10535885\tlost_fwd_first\tA\tATCCTAGGA",
            "21\t10535886\tlost_fwd_second\tT\tTCCTAGGAT",
        ],
        "ENST00001091162",
    );
    let reverse = annotate(
        &[
            "21\t14968490\tlost_rev_second\tC\tCATCCTAGG",
            "21\t14968491\tlost_rev_first\tA\tATCCTAGGA",
            "21\t14968492\tlost_rev_before\tT\tTCCTAGGAT",
        ],
        "ENST00000318948",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        let second = row(rows, &format!("lost_{strand}_second"));
        assert_eq!(second.0, LOST_STOP, "{strand}: {}", second.1);
        assert!(
            second.1.starts_with("AltStop=True:"),
            "{strand}: {}",
            second.1
        );
        for way in ["before", "first"] {
            let name = format!("lost_{strand}_{way}");
            assert_eq!(row(rows, &name), second, "{name}");
        }
    }
}

/// Four bases written after the second base of a start codon they move on,
/// the first two completing a stop codon with that base: the upstream ORF
/// whose stop codon begins there still ends on a stop codon at that index,
/// with the length it had. No ORF is shifted and the row reports nothing,
/// whichever of the three ways the record writes the sequence.
#[test]
fn insertion_that_moves_a_start_codon_on_and_keeps_a_stop_codon_inside_it_shifts_nothing() {
    let forward = annotate(
        &[
            "21\t10535884\tkept_fwd_before\tC\tCATAA",
            "21\t10535885\tkept_fwd_first\tA\tATAAA",
            "21\t10535886\tkept_fwd_second\tT\tTAAAT",
        ],
        "ENST00001091162",
    );
    let reverse = annotate(
        &[
            "21\t14968490\tkept_rev_second\tC\tCATTT",
            "21\t14968491\tkept_rev_first\tA\tATTTA",
            "21\t14968492\tkept_rev_before\tT\tTTTAT",
        ],
        "ENST00000318948",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        for way in ["before", "first", "second"] {
            let name = format!("kept_{strand}_{way}");
            assert_eq!(row(rows, &name), ("-", "-"), "{name}");
        }
    }
}

/// Five bases inserted after a start codon's first base that hold a start
/// codon from their second base: the row reports that codon gained, two
/// bases past the reference codon's index, whether the inserted bases end
/// with the reference codon's first base, leaving it whole five bases on, or
/// part it.
#[test]
fn start_codon_made_inside_a_start_codon_is_reported_gained() {
    let forward = annotate(
        &[
            "21\t10524583\tmoved_fwd\tA\tACATGA",
            "21\t10524583\tparted_fwd\tA\tACATGC",
        ],
        "ENST00000427445",
    );
    for name in ["moved_fwd", "parted_fwd"] {
        let (consequence, annotation) = row(&forward, name);
        assert_eq!(consequence, GAINED_START, "{name}");
        assert!(
            annotation.starts_with("CapDistanceToStart=147:"),
            "{name}: {annotation}"
        );
    }
    let reverse = annotate(
        &[
            "21\t14210906\tmoved_rev\tA\tATCATG",
            "21\t14210906\tparted_rev\tA\tAGCATG",
        ],
        "ENST00000400211",
    );
    for name in ["moved_rev", "parted_rev"] {
        let (consequence, annotation) = row(&reverse, name);
        assert_eq!(consequence, GAINED_START, "{name}");
        assert!(
            annotation.starts_with("CapDistanceToStart=50:"),
            "{name}: {annotation}"
        );
    }
}

/// An insertion directly before a transcript's first base is no 5' UTR
/// variant of that transcript, so its row carries no field of the plugin.
#[test]
fn insertion_before_a_transcript_reports_no_plugin_field() {
    let forward = annotate(&["21\t10521552\tedge_fwd\tG\tGCCC"], "ENST00000427445");
    assert_eq!(row(&forward, "edge_fwd"), ("", ""), "forward strand");
    let reverse = annotate(&["21\t14210955\tedge_rev\tT\tTGGG"], "ENST00000400211");
    assert_eq!(row(&reverse, "edge_rev"), ("", ""), "reverse strand");
}

/// Three bases inserted after a start codon's first base that open with its
/// last two, or after its second base that open with its last, leave the
/// codon at its own index with three bases directly after it. Nothing is
/// gained and no frame shifts: the two fields report nothing, as they do for
/// the same mutated sequence written directly after the codon.
#[test]
fn insertion_that_leaves_a_start_codon_at_its_own_index_reports_nothing() {
    let forward = annotate(
        &[
            "21\t10524583\tfirst_fwd\tA\tATGC",
            "21\t10524584\tsecond_fwd\tT\tTGCC",
            "21\t10524585\tfirst_fwd_after\tG\tGCTG",
            "21\t10524585\tsecond_fwd_after\tG\tGCCG",
        ],
        "ENST00000427445",
    );
    let reverse = annotate(
        &[
            "21\t14210904\tfirst_rev_after\tT\tTCAG",
            "21\t14210904\tsecond_rev_after\tT\tTCGG",
            "21\t14210905\tsecond_rev\tC\tCGGC",
            "21\t14210906\tfirst_rev\tA\tAGCA",
        ],
        "ENST00000400211",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        for name in [format!("first_{strand}"), format!("second_{strand}")] {
            assert_eq!(row(rows, &name), ("-", "-"), "{name}");
            assert_eq!(
                row(rows, &name),
                row(rows, &format!("{name}_after")),
                "{name}"
            );
        }
    }
}

/// Two bases after a start codon's first base that repeat its last two, and
/// one after its second that repeats its last, stand directly after the
/// codon and shift its ORF: the row reports the frameshift, field for field
/// the row of the same bases written directly after the codon, and names no
/// gained start.
#[test]
fn insertion_that_leaves_a_start_codon_at_its_own_index_reports_its_frameshift() {
    let forward = annotate(
        &[
            "21\t10524583\tfirst_fwd\tA\tATG",
            "21\t10524584\tsecond_fwd\tT\tTG",
            "21\t10524585\tfirst_fwd_after\tG\tGTG",
            "21\t10524585\tsecond_fwd_after\tG\tGG",
        ],
        "ENST00000427445",
    );
    let reverse = annotate(
        &[
            "21\t14210904\tfirst_rev_after\tT\tTCA",
            "21\t14210904\tsecond_rev_after\tT\tTC",
            "21\t14210905\tsecond_rev\tC\tCC",
            "21\t14210906\tfirst_rev\tA\tACA",
        ],
        "ENST00000400211",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        for name in [format!("first_{strand}"), format!("second_{strand}")] {
            let (consequence, annotation) = row(rows, &name);
            assert_eq!(consequence, FRAMESHIFT, "{name}: {annotation}");
            assert!(
                !annotation.contains("CapDistanceToStart="),
                "{name}: {annotation}"
            );
            assert_eq!(
                row(rows, &name),
                row(rows, &format!("{name}_after")),
                "{name}"
            );
        }
    }
}

/// Inserted bases that leave a start codon at its own index and hold a start
/// codon of their own: the row reports the codon they make, three bases past
/// the reference codon's index, as the row of the same bases written
/// directly after the codon does.
#[test]
fn start_codon_made_after_one_left_at_its_own_index_is_reported_gained() {
    let forward = annotate(
        &[
            "21\t10524583\tfirst_fwd\tA\tATGATGC",
            "21\t10524584\tsecond_fwd\tT\tTGATGC",
            "21\t10524585\tfirst_fwd_after\tG\tGATGCTG",
            "21\t10524585\tsecond_fwd_after\tG\tGATGCG",
        ],
        "ENST00000427445",
    );
    let reverse = annotate(
        &[
            "21\t14210904\tfirst_rev_after\tT\tTCAGCAT",
            "21\t14210904\tsecond_rev_after\tT\tTCGCAT",
            "21\t14210905\tsecond_rev\tC\tCGCATC",
            "21\t14210906\tfirst_rev\tA\tAGCATCA",
        ],
        "ENST00000400211",
    );
    for (rows, strand, gained) in [(&forward, "fwd", "148"), (&reverse, "rev", "51")] {
        for name in [format!("first_{strand}"), format!("second_{strand}")] {
            let (consequence, annotation) = row(rows, &name);
            assert_eq!(consequence, GAINED_START, "{name}");
            assert!(
                annotation.starts_with(&format!("CapDistanceToStart={gained}:")),
                "{name}: {annotation}"
            );
            assert_eq!(
                row(rows, &name),
                row(rows, &format!("{name}_after")),
                "{name}"
            );
        }
    }
}

/// One base after a start codon's second base that repeats its third leaves
/// the codon at its own index and stands inside the stop codon of an
/// upstream ORF, which begins at the codon's second base: the row reports
/// that ORF's lost stop, as the row of the same base written directly after
/// the codon does.
#[test]
fn insertion_that_leaves_a_start_codon_at_its_own_index_inside_a_stop_codon_reports_the_lost_stop()
{
    let forward = annotate(
        &[
            "21\t10535886\tstop_fwd\tT\tTG",
            "21\t10535887\tstop_fwd_after\tG\tGG",
        ],
        "ENST00001091162",
    );
    let reverse = annotate(
        &[
            "21\t14968489\tstop_rev_after\tT\tTC",
            "21\t14968490\tstop_rev\tC\tCC",
        ],
        "ENST00000318948",
    );
    for (rows, strand) in [(&forward, "fwd"), (&reverse, "rev")] {
        let name = format!("stop_{strand}");
        assert_eq!(row(rows, &name).0, LOST_STOP, "{name}");
        assert_eq!(
            row(rows, &name),
            row(rows, &format!("{name}_after")),
            "{name}"
        );
    }
}

/// Three bases that leave a start codon at the reference codon's index and
/// end with that codon's first bases repeat the codon: one start codon more
/// stands in the UTR, and the row reports the gained start at the reference
/// codon's index, as the row of a start codon written directly before the
/// codon does.
#[test]
fn insertion_that_repeats_a_start_codon_reports_the_gain_at_its_index() {
    let forward = annotate(
        &[
            "21\t10524582\tbefore_fwd\tA\tAATG",
            "21\t10524583\tfirst_fwd\tA\tATGA",
            "21\t10524584\tsecond_fwd\tT\tTGAT",
        ],
        "ENST00000427445",
    );
    let reverse = annotate(
        &[
            "21\t14210905\tsecond_rev\tC\tCATC",
            "21\t14210906\tfirst_rev\tA\tATCA",
            "21\t14210907\tbefore_rev\tT\tTCAT",
        ],
        "ENST00000400211",
    );
    for (rows, strand, gained) in [(&forward, "fwd", "145"), (&reverse, "rev", "48")] {
        let before = row(rows, &format!("before_{strand}"));
        assert_eq!(before.0, GAINED_START, "{strand}");
        assert!(
            before
                .1
                .starts_with(&format!("CapDistanceToStart={gained}:")),
            "{strand}: {}",
            before.1
        );
        for name in [format!("first_{strand}"), format!("second_{strand}")] {
            assert_eq!(row(rows, &name), before, "{name}");
        }
    }
}
