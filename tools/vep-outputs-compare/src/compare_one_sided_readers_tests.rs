// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the one-sided classes of the 116.2 reference read of an output pair beside its keys:
//! the record names behind a Location the reference wrote nothing at, the reference's
//! duplicate rows of a bare breakpoint, the gVCF variant sites an output or an input VCF
//! names, and the input VCF behind an output; each followed through the merge into the pass.

use super::*;
use crate::rules::filter::{BucketStats, Exclusions};
use crate::rules::scope::Group;

const OUTPUT_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "one_sided_reader_test_").unwrap()
}

/// Writes an output file of the given rows under the column header.
fn output(scratch: &TempDir, name: &str, rows: &str) -> PathBuf {
    let path = scratch.path.join(name);
    std::fs::write(&path, format!("{OUTPUT_COLUMNS}{rows}")).unwrap();
    path
}

/// One row of an output: the record name, the Location, the allele and the feature.
fn line(name: &str, location: &str, allele: &str, feature: &str) -> String {
    format!("{name}\t{location}\t{allele}\tg1\t{feature}\tTranscript\tfeature_truncation\n")
}

fn locations(items: &[&str]) -> HashSet<String> {
    items.iter().map(|item| item.to_string()).collect()
}

fn groups(file_name: &str, items: &[&str]) -> HashSet<Group> {
    items
        .iter()
        .map(|location| (file_name.to_string(), location.to_string()))
        .collect()
}

/// The sorted unique keys of the two outputs, each in a file of its own.
fn sorted_keys(scratch: &TempDir, perl: &Path, rust: &Path) -> (PathBuf, PathBuf) {
    let sorted = |output: &Path, name: &str| {
        let keys = scratch.path.join(name);
        sort_keys(output, false, None, &scratch.path, &keys).unwrap();
        keys
    };
    (sorted(perl, "perl.keys"), sorted(rust, "rust.keys"))
}

/// Merges the two key files as the file `f`, feeding the 116.2 collector.
fn merged(perl_keys: &Path, rust_keys: &Path, scope: &mut ReferenceRecordScope) -> MergeCounts {
    let collector = Collector::Release116(scope);
    compare_sorted_key_files(perl_keys, rust_keys, "f", None, false, collector).unwrap()
}

/// One discordant row of the file `f`; a row on no feature is an intergenic one.
fn discordant(location: &str, allele: &str, feature: &str, consequence: &str) -> String {
    let feature_type = if feature == "-" { "-" } else { "Transcript" };
    format!("f\textra_in_rust\t{location}\t{allele}\t{feature}\t{feature_type}\t{consequence}\n")
}

/// Runs the 116.2 pass over discordant rows and returns its statistics with the Location of
/// every row that stays open.
fn record_pass(rows: &[String], scope: &ReferenceRecordScope) -> (Exclusions, Vec<String>) {
    let scratch = scratch();
    let discordant = scratch.path.join("discordant.tsv");
    let open = scratch.path.join("discordant_open.tsv");
    let columns = "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set";
    std::fs::write(&discordant, format!("{columns}\n{}", rows.concat())).unwrap();
    let exclusions = filter_intended_divergences(
        &discordant,
        &open,
        FilterCollector::Release116(scope, &CacheContigNames::Human),
        &row_pair_rules(),
        None,
        None,
    )
    .unwrap();
    let open_locations = std::fs::read_to_string(&open)
        .unwrap()
        .lines()
        .skip(1)
        .map(|row| row.split('\t').nth(2).expect("a Location").to_string())
        .collect();
    (exclusions, open_locations)
}

fn class<'a>(exclusions: &'a Exclusions, name: &str) -> &'a BucketStats {
    let (_, stats) = exclusions
        .buckets
        .iter()
        .find(|(bucket, _)| bucket == name)
        .expect("the class");
    stats
}

fn one_sided(count: u64, excluded: u64, records: u64, perl: u64, rust: u64) -> BucketStats {
    vec![
        ("count", count),
        ("excluded", excluded),
        ("records", records),
        ("perl", perl),
        ("rust", rust),
    ]
}

#[test]
fn outputs_with_no_candidate_location_are_not_read_again() {
    let scratch = scratch();
    let (perl, rust) = (scratch.path.join("perl.txt"), scratch.path.join("rust.txt"));
    assert_eq!(
        reference_named_locations(&perl, &rust, &locations(&[]), None).unwrap(),
        locations(&[]),
        "neither output exists, and neither is opened"
    );
    assert!(
        reference_named_locations(&perl, &rust, &locations(&["21:100"]), None).is_err(),
        "a candidate Location opens the outputs"
    );
}

#[test]
fn the_named_locations_are_those_holding_a_name_the_reference_wrote() {
    let scratch = scratch();
    let perl = output(
        &scratch,
        "perl.txt",
        &line("sv_kept", "21:5001-5300", "deletion", "T1"),
    );
    let rust = output(
        &scratch,
        "rust.txt",
        &[
            line("sv_kept", "21:5000-5300", "<DEL>", "T1"),
            line("sv_dropped", "21:7000-27000100", "<DEL>", "T1"),
            line("sv_kept", "21:1-5", "<DEL>", "T9"),
        ]
        .concat(),
    );
    let candidates = ["21:5000-5300", "21:7000-27000100"];
    assert_eq!(
        reference_named_locations(&perl, &rust, &locations(&candidates), None).unwrap(),
        locations(&["21:5000-5300"])
    );
    let mut scope = ReferenceRecordScope::default();
    scope.absent = groups("f", &candidates);
    scope.absent.extend(groups("g", &["21:5000-5300"]));
    assert_eq!(scope.absent_locations("f"), locations(&candidates));
    scope.record_reference_named("f", &locations(&["21:5000-5300"]));
    let mut left = groups("f", &["21:7000-27000100"]);
    left.extend(groups("g", &["21:5000-5300"]));
    assert_eq!(scope.absent, left);
}

#[test]
fn the_reference_allele_rows_of_a_bare_breakpoint_are_set_aside_and_counted() {
    let own = "bare_bnd\t21:30425651\tchromosome_breakpoint\tG\tENST00000390690\tTranscript\tcoding_sequence_variant\n";
    let up = "bare_bnd\t21:30425651\tchromosome_breakpoint\tG\tENST00000334068\tTranscript\tupstream_gene_variant\n";
    let mate = "paired\t21:30425651\tN[21:30425700[\tG\tENST00000390690\tTranscript\tcoding_sequence_variant\n";
    let scratch = scratch();
    let perl = output(
        &scratch,
        "perl.txt",
        &[own, own, up, up, mate, mate].concat(),
    );
    let rust = output(&scratch, "rust.txt", &[own, up, mate].concat());
    let (perl_keys, rust_keys) = sorted_keys(&scratch, &perl, &rust);
    let mut scope = ReferenceRecordScope::default();
    let counts = merged(&perl_keys, &rust_keys, &mut scope);
    scope.record_bare_reference_rows(
        "f",
        &count_bare_breakpoint_reference_rows(&perl, None).unwrap(),
    );
    assert_eq!(
        (
            counts.perl_count,
            counts.rust_count,
            counts.intersection,
            counts.missing,
            counts.extra
        ),
        (3, 3, 3, 0, 0),
        "the reference-allele rows enter neither figure"
    );
    assert_eq!(
        scope.bare_reference_rows,
        HashMap::from([(("f".to_string(), "21:30425651".to_string()), 2)]),
        "the bracket record's duplicate mate rows are not the shape"
    );
    let (exclusions, open_locations) = record_pass(&[], &scope);
    assert_eq!((exclusions.excluded_perl, exclusions.excluded_rust), (0, 0));
    assert!(open_locations.is_empty(), "{open_locations:?}");
    assert_eq!(
        class(&exclusions, "bnd_bare_reference_allele_row"),
        &one_sided(2, 2, 1, 0, 0)
    );
}

#[test]
fn a_gvcf_variant_site_the_reference_dropped_is_set_aside_from_its_own_name() {
    let site = "21_30425650_N/G/<NON_REF>\t21:30425650\tG\tG\tENST00000390690\tTranscript\tcoding_sequence_variant\n\
                21_30425650_N/G/<NON_REF>\t21:30425650\tG\tG\tENST00000334068\tTranscript\tupstream_gene_variant\n";
    let named = "rs1\t21:40000000\tT\tG\tENST_X\tTranscript\tmissense_variant\n";
    let block =
        "21_50000000_N/<NON_REF>\t21:50000001-50000100\t<NON_REF>\t-\t-\t-\tintergenic_variant\n";
    let scratch = scratch();
    let perl = output(&scratch, "perl.txt", "");
    let rust = output(&scratch, "rust.txt", &[site, named, block].concat());
    let sites = gvcf_site_locations(&rust, None).unwrap();
    assert_eq!(sites, locations(&["21:30425650"]));
    assert_eq!(
        gvcf_site_locations(&scratch.path.join("absent.txt"), None).unwrap(),
        locations(&[])
    );
    let (perl_keys, rust_keys) = sorted_keys(&scratch, &perl, &rust);
    let mut scope = ReferenceRecordScope::default();
    scope.record_gvcf_sites("f", &sites);
    merged(&perl_keys, &rust_keys, &mut scope);
    assert_eq!(scope.gvcf_absent, groups("f", &["21:30425650"]));
    assert_eq!(
        scope.absent,
        groups("f", &["21:50000001-50000100"]),
        "the lone block is the skipped record"
    );
    let rows = [
        discordant(
            "21:30425650",
            "G",
            "ENST00000390690",
            "coding_sequence_variant",
        ),
        discordant(
            "21:30425650",
            "G",
            "ENST00000334068",
            "upstream_gene_variant",
        ),
        discordant("21:40000000", "G", "ENST_X", "missense_variant"),
        discordant(
            "21:50000001-50000100",
            "<NON_REF>",
            "-",
            "intergenic_variant",
        ),
    ];
    let (exclusions, open_locations) = record_pass(&rows, &scope);
    assert_eq!((exclusions.excluded_perl, exclusions.excluded_rust), (0, 3));
    assert_eq!(
        class(&exclusions, "gvcf_variant_site_dropped"),
        &one_sided(2, 2, 1, 0, 2)
    );
    assert_eq!(
        class(&exclusions, "reference_skipped_record")[1],
        ("excluded", 1)
    );
    assert_eq!(open_locations, ["21:40000000"]);
}

#[test]
fn a_site_under_an_identifier_is_set_aside_once_the_input_names_it() {
    let named = "rs1\t21:30425651\t-\tG\tENST00000390690\tTranscript\tframeshift_variant\n";
    let blocks = "21_50000000_N/<NON_REF>/<NON_REF>\t21:50000001-50000100\t<NON_REF>\t-\t-\t-\tintergenic_variant\n";
    let scratch = scratch();
    let perl = output(&scratch, "perl.txt", "");
    let rust = output(&scratch, "rust.txt", &[named, blocks].concat());
    let input = scratch.path.join("s.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30425650\trs1\tATG\tAG,<NON_REF>\t.\t.\t.\n\
         21\t50000000\t.\tN\t<NON_REF>,<NON_REF>\t.\t.\tEND=50000100\n\
         21\t60000000\t.\tN\t<NON_REF>\t.\t.\tEND=60000100\n\
         chr22\t100\t.\tA\tG,<NON_REF>\t.\t.\t.\n",
    )
    .unwrap();
    let (perl_keys, rust_keys) = sorted_keys(&scratch, &perl, &rust);
    let rows = [
        discordant("21:30425651", "-", "ENST00000390690", "frameshift_variant"),
        discordant(
            "21:50000001-50000100",
            "<NON_REF>",
            "-",
            "intergenic_variant",
        ),
    ];
    for with_input in [false, true] {
        let mut sites = gvcf_site_locations(&rust, None).unwrap();
        if with_input {
            sites.extend(gvcf_site_locations_from_vcf(&input, None).unwrap());
        }
        let mut scope = ReferenceRecordScope::default();
        scope.record_gvcf_sites("f", &sites);
        merged(&perl_keys, &rust_keys, &mut scope);
        let (exclusions, open_locations) = record_pass(&rows, &scope);
        if with_input {
            assert_eq!(scope.gvcf_absent, groups("f", &["21:30425651"]));
            assert_eq!((exclusions.excluded_perl, exclusions.excluded_rust), (0, 2));
            assert!(open_locations.is_empty(), "{open_locations:?}");
            assert_eq!(
                class(&exclusions, "gvcf_variant_site_dropped")[4],
                ("rust", 1)
            );
        } else {
            assert_eq!(scope.gvcf_absent, groups("f", &[]));
            assert_eq!(
                open_locations,
                ["21:30425651"],
                "without the input the row under an identifier stays charged"
            );
        }
        assert_eq!(
            class(&exclusions, "reference_skipped_record")[1],
            ("excluded", 1),
            "the block list is a skipped record with the input and without it"
        );
    }
}

#[test]
fn the_input_behind_an_output_is_the_file_named_or_the_one_of_the_outputs_stem() {
    let scratch = scratch();
    let inputs = scratch.path.join("inputs");
    std::fs::create_dir(&inputs).unwrap();
    let (plain, gzipped) = (inputs.join("s.vcf"), inputs.join("z.vcf.gz"));
    std::fs::write(&plain, "").unwrap();
    std::fs::write(&gzipped, "").unwrap();
    assert_eq!(
        input_vcf_for("s.txt", Some(&inputs), None),
        Some(plain.clone())
    );
    assert_eq!(
        input_vcf_for("z.txt", Some(&inputs), None),
        Some(gzipped.clone())
    );
    assert_eq!(input_vcf_for("none.txt", Some(&inputs), None), None);
    assert_eq!(
        input_vcf_for("s.txt", Some(&inputs), Some(&gzipped)),
        Some(gzipped.clone()),
        "a file named for every output comes before the directory"
    );
    assert_eq!(input_vcf_for("s.txt", None, None), None);
}
