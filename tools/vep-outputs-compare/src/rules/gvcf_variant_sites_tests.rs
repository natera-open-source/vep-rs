// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! gVCF variant sites: which records are one, the rows the scored engine writes for each
//! and the Locations an output or an input VCF names them at.

use super::*;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::PathBuf;

const OUTPUT_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";
const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The 300 generated sites `scripts/parity/make_gvcf_site_fixture.py` writes, with the
/// (Location, Allele) rows the scored engine wrote for each.
const SITE_ROWS_FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../scripts/parity/gvcf_site_rows_fixture.tsv"
);

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "gvcf_site_test_").unwrap()
}

fn written(scratch: &TempDir, name: &str, content: &str) -> PathBuf {
    let path = scratch.path.join(name);
    std::fs::write(&path, content).unwrap();
    path
}

fn locations(items: &[&str]) -> HashSet<String> {
    items.iter().map(|item| item.to_string()).collect()
}

/// Rows as `gvcf_site_rows` returns them: (allele, start, end).
fn rows(items: &[(&str, i128, i128)]) -> Vec<(String, i128, i128)> {
    items
        .iter()
        .map(|(allele, start, end)| (allele.to_string(), *start, *end))
        .collect()
}

/// The (Location, allele) of every row of a site at 21:30425650, sorted.
fn placed(reference: &str, alts: &str) -> Vec<(String, String)> {
    let alts: Vec<&str> = alts.split(',').collect();
    let mut placed: Vec<(String, String)> = gvcf_site_rows(30_425_650, reference, &alts)
        .into_iter()
        .map(|(allele, start, end)| (sequence_location("21", start, end), allele))
        .collect();
    placed.sort_unstable();
    placed
}

fn pairs(items: &[(&str, &str)]) -> Vec<(String, String)> {
    items
        .iter()
        .map(|(location, allele)| (location.to_string(), allele.to_string()))
        .collect()
}

#[test]
fn a_block_beside_an_allele_that_names_a_variant_is_a_site() {
    assert!(is_gvcf_variant_site(&["G", "<NON_REF>"], None));
    assert!(is_gvcf_variant_site(&["<non_ref>", "AT"], None));
    assert!(!is_gvcf_variant_site(&["<NON_REF>"], None));
    assert!(!is_gvcf_variant_site(&["<NON_REF>", "<NON_REF>"], None));
    assert!(!is_gvcf_variant_site(&["*", "<NON_REF>"], None));
    assert!(!is_gvcf_variant_site(&["G", "T"], None));
    assert!(
        !is_gvcf_variant_site(&["<DEL>", "<NON_REF>"], None),
        "a block beside a symbolic allele is an unsupported structural list"
    );
    assert!(!is_gvcf_variant_site(&["G", "<DEL>", "<NON_REF>"], None));
    assert!(
        !is_gvcf_variant_site(&["G", "<NON_REF>"], Some("DEL")),
        "under a true SVTYPE the block is an allele of a structural record"
    );
    assert!(is_gvcf_variant_site(&["G", "<NON_REF>"], Some("0")));
}

#[test]
fn an_input_names_a_site_by_its_shape_whatever_its_id_plain_or_gzipped() {
    let scratch = scratch();
    let named = "rs1\t21:30425651\t-\tG\tENST00000390690\tTranscript\tframeshift_variant\n";
    let blocks = "21_50000000_N/<NON_REF>/<NON_REF>\t21:50000001-50000100\t<NON_REF>\t-\t-\t-\tintergenic_variant\n";
    let output = written(
        &scratch,
        "s.txt",
        &format!("{OUTPUT_COLUMNS}{named}{blocks}"),
    );
    assert_eq!(
        gvcf_site_locations(&output, None).unwrap(),
        locations(&[]),
        "a list of blocks and a row under an identifier name no site"
    );
    let records = format!(
        "{VCF_HEADER}\
         21\t30425650\trs1\tATG\tAG,<NON_REF>\t.\t.\t.\n\
         21\t50000000\t.\tN\t<NON_REF>,<NON_REF>\t.\t.\tEND=50000100\n\
         21\t60000000\t.\tN\t<NON_REF>\t.\t.\tEND=60000100\n\
         chr22\t100\t.\tA\tG,<NON_REF>\t.\t.\t.\n"
    );
    let plain = written(&scratch, "s.vcf", &records);
    let gzipped = scratch.path.join("z.vcf.gz");
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(records.as_bytes()).unwrap();
    std::fs::write(&gzipped, gzip.finish().unwrap()).unwrap();
    for path in [&plain, &gzipped] {
        assert_eq!(
            gvcf_site_locations_from_vcf(path, None).unwrap(),
            locations(&["21:30425651", "chr22:100"]),
            "{}",
            path.display()
        );
    }
    let contigs = locations(&["21"]);
    assert_eq!(
        gvcf_site_locations_from_vcf(&plain, Some(&contigs)).unwrap(),
        locations(&["21:30425651"])
    );
}

#[test]
fn a_site_is_placed_at_the_pos_as_the_scored_engine_reads_it_and_a_refused_pos_names_none() {
    let read = [
        "3.042565e7",
        "30425650.0",
        "+30425650",
        "30425650.",
        " 30425650 ",
    ];
    let refused = [
        "30_425_650",
        "-30425650",
        "abc",
        "３０４２５６５０",
        "0x1D04A1A",
        "inf",
        "1e400",
    ];
    let sites = |positions: &[&str]| -> String {
        positions
            .iter()
            .enumerate()
            .map(|(at, pos)| format!("21\t{pos}\tp{at}\tATG\tAG,<NON_REF>\t.\t.\t.\n"))
            .collect()
    };
    let scratch = scratch();
    let every_spelling = [read.as_slice(), refused.as_slice()].concat();
    let path = written(
        &scratch,
        "pos.vcf",
        &format!("{VCF_HEADER}{}", sites(&every_spelling)),
    );
    assert_eq!(
        gvcf_site_locations_from_vcf(&path, None).unwrap(),
        locations(&["21:30425651"])
    );
    let path = written(
        &scratch,
        "refused.vcf",
        &format!("{VCF_HEADER}{}", sites(&refused)),
    );
    assert_eq!(
        gvcf_site_locations_from_vcf(&path, None).unwrap(),
        locations(&[])
    );
}

#[test]
fn a_sites_allele_is_written_in_upper_case() {
    assert_eq!(
        gvcf_site_rows(30_425_650, "A", &["g", "<NON_REF>"]),
        rows(&[("G", 30_425_650, 30_425_650)])
    );
    assert_eq!(
        gvcf_site_rows(30_425_650, "a", &["g", "<non_ref>"]),
        rows(&[("G", 30_425_650, 30_425_650)])
    );
    assert_eq!(
        gvcf_site_rows(30_425_650, "atg", &["ag", "<NON_REF>"]),
        rows(&[("-", 30_425_651, 30_425_651)])
    );
}

#[test]
fn a_star_or_a_dot_beside_a_sites_indel_keeps_the_anchor_trimmed_alleles() {
    let pos = 30_425_650;
    assert_eq!(
        gvcf_site_rows(pos, "ATT", &["AT", "*", "<NON_REF>"]),
        rows(&[("T", 30_425_651, 30_425_652)])
    );
    assert_eq!(
        gvcf_site_rows(pos, "AT", &["ATT", "*", "<NON_REF>"]),
        rows(&[("TT", 30_425_651, 30_425_651)])
    );
    assert_eq!(
        gvcf_site_rows(pos, "ATTG", &["ATC", "*", "<NON_REF>"]),
        rows(&[("TC", 30_425_651, 30_425_653)])
    );
    assert_eq!(
        gvcf_site_rows(pos, "ATT", &["AT", ".", "<NON_REF>"]),
        rows(&[("AT", 30_425_650, 30_425_652)]),
        "a `.` is compared and never shares the first base, so the anchor stays"
    );
    assert_eq!(
        gvcf_site_rows(pos, "ATT", &["AT", "<NON_REF>"]),
        rows(&[("-", 30_425_652, 30_425_652)]),
        "the one allele beside the block is minimised from both ends"
    );
    assert_eq!(
        gvcf_site_rows(pos, "AT", &["ATT", "<NON_REF>"]),
        rows(&[("T", 30_425_652, 30_425_651)])
    );
    assert_eq!(
        gvcf_site_rows(pos, "ATTG", &["ATC", "<NON_REF>"]),
        rows(&[("C", 30_425_652, 30_425_653)])
    );
}

#[test]
fn a_star_counts_as_an_allele_of_one_base_in_a_sites_anchor_trim() {
    for (reference, alts, expected) in [
        ("A", "G,*,<NON_REF>", vec![("21:30425650", "G")]),
        ("AG", "AC,*,<NON_REF>", vec![("21:30425651", "C")]),
        ("AG", "TC,*,<NON_REF>", vec![("21:30425650-30425651", "TC")]),
        (
            "AGC",
            "ATT,*,<NON_REF>",
            vec![("21:30425651-30425652", "TT")],
        ),
        ("AG", "A,*,<NON_REF>", vec![("21:30425651", "-")]),
        ("A", "AG,*,<NON_REF>", vec![("21:30425650-30425651", "G")]),
        ("AG", "*,AC,<NON_REF>", vec![("21:30425651", "C")]),
        (
            "AG",
            "AC,AT,*,<NON_REF>",
            vec![("21:30425651", "C"), ("21:30425651", "T")],
        ),
        ("A", "*,G,<NON_REF>", vec![("21:30425650", "G")]),
        ("AG", "*,A,<NON_REF>", vec![("21:30425651", "-")]),
        (
            "A",
            "G,T,*,<NON_REF>",
            vec![("21:30425650", "G"), ("21:30425650", "T")],
        ),
        (
            "AG",
            "AC,AT,<NON_REF>",
            vec![
                ("21:30425650-30425651", "AC"),
                ("21:30425650-30425651", "AT"),
            ],
        ),
        ("AG", "AC,.,<NON_REF>", vec![("21:30425650-30425651", "AC")]),
        (
            "AGC",
            "ATT,.,<NON_REF>",
            vec![("21:30425650-30425652", "ATT")],
        ),
        (
            "A",
            "AG,<*>,<NON_REF>",
            vec![
                ("21:30425650-30425651", "<*>"),
                ("21:30425650-30425651", "G"),
            ],
        ),
    ] {
        assert_eq!(
            placed(reference, alts),
            pairs(&expected),
            "{reference} {alts}"
        );
    }
}

#[test]
fn the_rows_of_a_site_follow_the_sequence_path() {
    for (reference, alts, expected) in [
        ("A", vec!["G", "<NON_REF>"], vec![("G", 100, 100)]),
        ("AT", vec!["GC", "<NON_REF>"], vec![("GC", 100, 101)]),
        ("A", vec!["AT", "<NON_REF>"], vec![("T", 101, 100)]),
        ("N", vec!["NA", "<NON_REF>"], vec![("A", 101, 100)]),
        ("ATG", vec!["A", "<NON_REF>"], vec![("-", 101, 102)]),
        ("ATG", vec!["AG", "<NON_REF>"], vec![("-", 101, 101)]),
        (
            "A",
            vec!["AT", "ATT", "<NON_REF>"],
            vec![("T", 101, 100), ("TT", 101, 100)],
        ),
        (
            "A",
            vec!["G", "T", "<NON_REF>"],
            vec![("G", 100, 100), ("T", 100, 100)],
        ),
        ("A", vec!["<NON_REF>", "<NON_REF>"], vec![]),
        ("A", vec!["*", "G", "<NON_REF>"], vec![("G", 100, 100)]),
    ] {
        assert_eq!(
            gvcf_site_rows(100, reference, &alts),
            rows(&expected),
            "{reference} {alts:?}"
        );
    }
    assert_eq!(
        sequence_location("21", 101, 100),
        "21:100-101",
        "an insertion ends one base before its start"
    );
    assert_eq!(sequence_location("21", 101, 102), "21:101-102");
    assert_eq!(sequence_location("chr21", 100, 100), "chr21:100");
}

#[test]
fn a_trailing_empty_alt_field_is_no_field() {
    assert_eq!(perl_split("<DEL>,", ','), ["<DEL>"]);
    assert_eq!(perl_split("G,,T,,", ','), ["G", "", "T"]);
    assert!(perl_split("", ',').is_empty());
    assert!(perl_split(",", ',').is_empty());
    assert_eq!(perl_split("N/G/<NON_REF>/", '/'), ["N", "G", "<NON_REF>"]);
    let mut fields = vec!["G", "<NON_REF>", "", ""];
    drop_trailing_empty(&mut fields);
    assert_eq!(fields, ["G", "<NON_REF>"]);
    assert!(!is_gvcf_variant_site(&["G", ""], None));
    assert!(is_gvcf_variant_site(&["G", "<NON_REF>", ""], None));
    assert!(
        !is_gvcf_variant_site(&["<NON_REF>", ""], None),
        "a block alone is no site"
    );
    assert_eq!(
        gvcf_site_rows(100, "A", &["G", ""]),
        rows(&[("G", 100, 100)])
    );
    assert_eq!(
        gvcf_site_rows(100, "A", &["G", "<NON_REF>", ""]),
        rows(&[("G", 100, 100)])
    );
    let scratch = scratch();
    let output = written(
        &scratch,
        "s.txt",
        &format!(
            "{OUTPUT_COLUMNS}\
             21_30425650_A/G/<NON_REF>/\t21:30425650\tG\tG\tENST_1\tTranscript\tmissense_variant\n\
             21_40000000_A/G/\t21:40000000\tG\tG\tENST_2\tTranscript\tmissense_variant\n"
        ),
    );
    let input = written(
        &scratch,
        "s.vcf",
        &format!(
            "{VCF_HEADER}\
             21\t30425650\trs1\tA\tG,<NON_REF>,\t.\t.\t.\n\
             21\t40000000\trs2\tA\tG,\t.\t.\t.\n\
             21\t50000000\trs3\tA\t<NON_REF>,\t.\t.\tEND=50000100\n"
        ),
    );
    assert_eq!(
        gvcf_site_locations(&output, None).unwrap(),
        locations(&["21:30425650"])
    );
    assert_eq!(
        gvcf_site_locations_from_vcf(&input, None).unwrap(),
        locations(&["21:30425650"])
    );
}

#[test]
fn an_info_value_is_read_up_to_a_second_equals_sign() {
    assert_eq!(info_written("SVTYPE==DEL;END=30425949", "SVTYPE"), Some(""));
    assert_eq!(
        info_written("SVTYPE==DEL;END=30425949", "END"),
        Some("30425949")
    );
    assert_eq!(
        info_written("SVTYPE=BND;CHR2=21;END2=30425700=5", "END2"),
        Some("30425700")
    );
    assert_eq!(info_written("END2=30425700=5=7", "END2"), Some("30425700"));
    assert_eq!(info_written("SVTYPE==DEL;END", "END"), None);
    assert!(is_gvcf_variant_site(
        &["G", "<NON_REF>"],
        info_written("SVTYPE==DEL", "SVTYPE")
    ));
    let scratch = scratch();
    let input = written(
        &scratch,
        "s.vcf",
        &format!("{VCF_HEADER}21\t30425650\trs1\tA\tG,<NON_REF>\t.\t.\tSVTYPE==DEL;END=30425949\n"),
    );
    assert_eq!(
        gvcf_site_locations_from_vcf(&input, None).unwrap(),
        locations(&["21:30425650"]),
        "an empty SVTYPE is false, so the record takes the sequence path"
    );
}

#[test]
fn the_rows_predicted_for_the_generated_sites_are_the_rows_the_scored_engine_wrote() {
    let text = std::fs::read_to_string(SITE_ROWS_FIXTURE)
        .unwrap_or_else(|e| panic!("{SITE_ROWS_FIXTURE}: {e}"));
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("pos\tref\talt_column\tallele\tlocation"));
    // A site is its (POS, REF, ALT column); its rows are (Location, allele) pairs.
    type Rows = BTreeSet<(String, String)>;
    let mut sites: BTreeMap<(i128, &str, &str), Rows> = BTreeMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let [pos, reference, alt_column, allele, location] = fields[..] else {
            panic!("a row of {} columns", fields.len());
        };
        let pos: i128 = pos.parse().expect("a position");
        sites
            .entry((pos, reference, alt_column))
            .or_default()
            .insert((location.to_string(), allele.to_string()));
    }
    assert_eq!(sites.len(), 300);
    assert!(sites
        .keys()
        .all(|(pos, _, _)| (25_000_000..=46_000_000).contains(pos)));
    assert!(sites
        .keys()
        .all(|(_, _, alt_column)| alt_column.split(',').any(|alt| alt == "<NON_REF>")));
    let mut differing = Vec::new();
    for ((pos, reference, alt_column), written) in &sites {
        let alts: Vec<&str> = alt_column.split(',').collect();
        let predicted: Rows = gvcf_site_rows(*pos, reference, &alts)
            .into_iter()
            .map(|(allele, start, end)| {
                let mut location = Vec::new();
                push_canonical_location(&sequence_location("21", start, end), &mut location);
                (String::from_utf8(location).unwrap(), allele)
            })
            .collect();
        if predicted != *written {
            differing.push(*pos);
        }
    }
    assert!(
        differing.is_empty(),
        "the positions of the sites whose predicted rows differ from the written ones: \
         {differing:?}"
    );
}
