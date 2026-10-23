// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `VARIANT_CLASS` column as the `vep` binary prints it: the `Allele` and
//! class of a structural record whose type the Sequence Ontology term table
//! lacks, and the class of a record `--check_frequency` rewrote, whole or
//! split under `--minimal`.

use std::path::Path;
use std::process::Stdio;

mod common;

/// The rows of a default-format run on `input` against the named corpus's cache
/// with `extra` flags, each row's identifier, allele and Extra column.
fn rows(
    corpus_name: &str,
    dir: &Path,
    input: &Path,
    extra: &[&str],
) -> Vec<(String, String, String)> {
    let corpus = common::corpora()
        .into_iter()
        .find(|corpus| corpus.name == corpus_name)
        .unwrap_or_else(|| panic!("the {corpus_name} corpus is committed"));
    let out = dir.join("out.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .args(corpus.cache_args())
        .args([
            "--assembly",
            &corpus.assembly,
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            (
                fields[0].to_string(),
                fields[2].to_string(),
                fields[13].to_string(),
            )
        })
        .collect()
}

/// A symbolic ALT with no Sequence Ontology term keeps Perl's type string as
/// its `class_SO_term` (Parser/VCF.pm:469-481: the ALTs joined by `/`, replaced
/// by INFO/SVTYPE when that names the type, kept whole when `get_SO_term` has no
/// term for it), printed in `Allele` and `VARIANT_CLASS` (OutputFactory.pm:967):
/// Ensembl VEP 116.2 under `--dont_skip --variant_class` prints `<CPX>` for
/// `21 19744457 cpx_bare N <CPX> . . END=19744524`, `CPX` for the same record
/// with `SVTYPE=CPX`, and `FOO` with `SVTYPE=FOO`; `<CTX>` reads like `<CPX>`.
#[test]
fn typeless_symbolic_alt_prints_perls_type_string() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("cpx.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t19744457\tcpx_bare\tN\t<CPX>\t.\t.\tEND=19744524\n\
         21\t19744457\tcpx_svtype\tN\t<CPX>\t.\t.\tSVTYPE=CPX;END=19744524\n\
         21\t19744457\tcpx_foo\tN\t<CPX>\t.\t.\tSVTYPE=FOO;END=19744524\n\
         21\t19744457\tctx_bare\tN\t<CTX>\t.\t.\tEND=19744524\n",
    )
    .unwrap();
    let rows = rows(
        "GRCh37-output-flags",
        dir.path(),
        &input,
        &["--dont_skip", "--variant_class"],
    );
    // One triple per record: vep-rs writes a kept record's row per overlapping
    // transcript, the coverage divergence docs/divergences/defect-c.md states.
    assert_eq!(
        classes(&rows),
        [
            ("cpx_bare", "<CPX>", "<CPX>"),
            ("cpx_svtype", "CPX", "CPX"),
            ("cpx_foo", "FOO", "FOO"),
            ("ctx_bare", "<CTX>", "<CTX>"),
        ]
    );
}

/// The `(identifier, allele, VARIANT_CLASS)` triples of `rows`, consecutive
/// duplicates collapsed, so a record's rows read as one triple.
fn classes(rows: &[(String, String, String)]) -> Vec<(&str, &str, &str)> {
    let mut classes: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|(id, allele, extra)| {
            let class = extra
                .split(';')
                .find_map(|kv| kv.strip_prefix("VARIANT_CLASS="))
                .unwrap_or_else(|| panic!("no VARIANT_CLASS on {id}: {extra}"));
            (id.as_str(), allele.as_str(), class)
        })
        .collect();
    classes.dedup();
    classes
}

/// A record `--check_frequency` rewrote prints the class of its kept alleles:
/// the filtering source annotates before the transcript source
/// (AnnotationSourceAdaptor.pm:98 sorts `can_filter_vfs` first), rewrites the
/// feature's allele string to the kept alleles (BaseCacheVariation.pm:355 under
/// `exclude`) and `class_SO_term` reads that string (VariationFeature.pm:1118).
/// Ensembl VEP 116.2 prints `VARIANT_CLASS=SNV` on all 16 rows of
/// `21 31148973 del_plus_snv GG G,GC` under the GRCh38-freqfilter corpus's
/// exclude set, whose common `-` allele (AF 0.6705) the filter drops, and
/// `sequence_alteration` for the same record without the filter.
#[test]
fn frequency_filtered_record_prints_the_class_of_its_kept_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("freq.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t31148973\tdel_plus_snv\tGG\tG,GC\t.\t.\t.\n",
    )
    .unwrap();
    let exclude = [
        "--check_frequency",
        "--freq_pop",
        "AF",
        "--freq_freq",
        "0.01",
        "--freq_gt_lt",
        "gt",
        "--freq_filter",
        "exclude",
        "--variant_class",
    ];
    let filtered = rows("GRCh38-freqfilter", dir.path(), &input, &exclude);
    assert_eq!(filtered.len(), 16);
    assert_eq!(classes(&filtered), [("del_plus_snv", "C", "SNV")]);
    assert!(filtered
        .iter()
        .all(|(_, _, extra)| extra.contains("FREQS=")));
    let unfiltered_rows = rows(
        "GRCh38-freqfilter",
        dir.path(),
        &input,
        &["--variant_class"],
    );
    let mut unfiltered = classes(&unfiltered_rows);
    unfiltered.sort_unstable();
    unfiltered.dedup();
    assert_eq!(
        unfiltered,
        [
            ("del_plus_snv", "-", "sequence_alteration"),
            ("del_plus_snv", "C", "sequence_alteration"),
        ]
    );
}

/// A split feature whose only row is its reference row prints the class of
/// its own allele string: `--minimal` splits the record per ALT
/// (InputBuffer.pm:638-718), the frequency check judges each feature
/// (BaseCacheVariation.pm:308-375), `--process_ref_homs` adds the reference
/// allele's rows (OutputFactory.pm:511) and `class_SO_term` reads the
/// feature's string (VariationFeature.pm:1118). Ensembl VEP 116.2 prints 48
/// rows for `21 31148973 . GG *,G,TG` under the exclude set below: the `GG/G`
/// feature's common `-` allele (AF 0.6705) is dropped whole, the `GG/*`
/// feature, whose `*` writes no row and has no frequency data, keeps 16 `GG`
/// rows with `VARIANT_CLASS=sequence_alteration` and `MINIMISED=1`, and the
/// `G/T` feature writes 16 `G` and 16 `T` rows with `SNV` and
/// `FREQS=AF:T:0.1396`.
#[test]
fn a_split_feature_left_with_its_reference_row_prints_its_own_class() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("star-first.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t31148973\t.\tGG\t*,G,TG\t.\t.\t.\n",
    )
    .unwrap();
    let exclude = [
        "--minimal",
        "--check_frequency",
        "--freq_pop",
        "AF",
        "--freq_freq",
        "0.5",
        "--freq_gt_lt",
        "gt",
        "--freq_filter",
        "exclude",
        "--process_ref_homs",
        "--variant_class",
        "--allele_number",
    ];
    let filtered = rows("GRCh38-freqfilter", dir.path(), &input, &exclude);
    assert_eq!(filtered.len(), 48);
    for (id, _, _) in &filtered {
        assert_eq!(id, "21_31148973_GG/*/G/TG");
    }
    let shapes: Vec<_> = filtered
        .iter()
        .map(|(_, allele, extra)| {
            let value = |key: &str| extra.split(';').find_map(|kv| kv.strip_prefix(key));
            (
                allele.as_str(),
                value("ALLELE_NUM="),
                value("VARIANT_CLASS="),
                value("MINIMISED="),
                value("FREQS="),
            )
        })
        .collect();
    for shape in [
        (
            "GG",
            Some("0"),
            Some("sequence_alteration"),
            Some("1"),
            None,
        ),
        ("G", Some("0"), Some("SNV"), None, Some("AF:T:0.1396")),
        ("T", Some("3"), Some("SNV"), None, Some("AF:T:0.1396")),
    ] {
        assert_eq!(
            shapes.iter().filter(|row| **row == shape).count(),
            16,
            "{shape:?}"
        );
    }
}
