// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How far a `--check_frequency` drop reaches under `--minimal`: a dropped
//! record or a dropped split feature cancels the rejoin of the split records
//! around it, which the `vep` binary then writes as one record per surviving
//! feature, and the unit the drop resets is the dropped record's run of one
//! chromosome inside a `--buffer_size` batch, where Ensembl VEP 116.2 resets
//! its input buffer.

use std::path::Path;
use std::process::Stdio;

mod common;

/// A data row's `(Uploaded_variation, Location, Allele, Extra)`.
type Row = (String, String, String, String);

/// A row's `(Location, Allele, ALLELE_NUM, MINIMISED, FREQS)`, each Extra key
/// `None` where the column lacks it.
type Shape<'a> = (
    &'a str,
    &'a str,
    Option<&'a str>,
    Option<&'a str>,
    Option<&'a str>,
);

/// The rows of a default-format run on `input` against the GRCh38-freqfilter
/// corpus's cache with `extra` flags.
fn rows(dir: &Path, input: &Path, extra: &[&str]) -> Vec<Row> {
    let corpus = common::corpora()
        .into_iter()
        .find(|corpus| corpus.name == "GRCh38-freqfilter")
        .expect("the GRCh38-freqfilter corpus is committed");
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
                fields[1].to_string(),
                fields[2].to_string(),
                fields[13].to_string(),
            )
        })
        .collect()
}

/// The shapes of the rows record `id` wrote, in the output's order.
fn shapes<'a>(rows: &'a [Row], id: &str) -> Vec<Shape<'a>> {
    rows.iter()
        .filter(|(row_id, _, _, _)| row_id == id)
        .map(|(_, location, allele, extra)| {
            let value = |key: &str| extra.split(';').find_map(|kv| kv.strip_prefix(key));
            (
                location.as_str(),
                allele.as_str(),
                value("ALLELE_NUM="),
                value("MINIMISED="),
                value("FREQS="),
            )
        })
        .collect()
}

/// Asserts record `id` wrote each `expected` shape as many times as stated and
/// no row of another shape.
fn assert_shapes(rows: &[Row], id: &str, expected: &[(Shape<'_>, usize)]) {
    let found = shapes(rows, id);
    for (shape, count) in expected {
        assert_eq!(
            found.iter().filter(|row| *row == shape).count(),
            *count,
            "{id} {shape:?}: {found:?}"
        );
    }
    assert_eq!(
        found.len(),
        expected.iter().map(|(_, count)| count).sum::<usize>(),
        "{id}: {found:?}"
    );
}

/// A `--check_frequency` drop under `--minimal` resets the rejoin of the split
/// records in the dropped record's run of one chromosome inside a
/// `--buffer_size` batch (`runs_reset_by_a_drop` in `pipeline.rs`), and a batch
/// is `--buffer_size` input lines counted from the file's first data line
/// (`read_vcf_lines`), a line later skipped as not found included. Ensembl VEP
/// 116.2 resets its input buffer (BaseCacheVariation.pm:365-368,
/// `if($filtered_count) { $buffer->reset_buffer();`, InputBuffer.pm:535-537,
/// `sub reset_buffer { $_[0]->{temp} = {}; }`), which holds one chromosome's
/// records up to `--buffer_size` (InputBuffer.pm:191-199, `while(@$pre_buffer
/// && @$buffer < $buffer_size) { ... if($prev_chr && $vf->{chr} ne $prev_chr) {
/// $self->split_variants() if $self->{minimal}; return $buffer; }`) and never a
/// record its parser skipped (Parser.pm:189); under `--fork N` the unit is a
/// chunk of at most `--buffer_size` over `2N` records (Runner.pm:486, `my
/// $maxForkSize = int($buffer_size / (2 * $fork_number)) || 1;`). Where the
/// units differ, one engine rejoins what the other resets
/// (`docs/output/known-variants.md`).
///
/// Eleven records under the GRCh38-freqfilter corpus's exclude set at 0.5: a
/// chromosome-22 line the cache lacks, `GG G,TG` at 21:31148973, whose `GG/G`
/// feature's common `-` allele (AF 0.6705) is dropped, eight `A T` records at
/// 21:32288094 (AF 0.0003994) and `GTG G,GAG` at 21:42029386, whose `GTG/G`
/// feature's `-` has AF 0.002995 and whose `GTG/GAG` feature has no data. Under
/// `--buffer_size 10` the first batch is lines 1 to 10, the skipped line
/// counted, and `GTG G,GAG` is line 11, a batch with no drop: vep-rs writes
/// `GG G,TG` reset, its `T` feature alone with `ALLELE_NUM=2`,
/// `FREQS=AF:T:0.1396` and no `MINIMISED`, and `GTG G,GAG` rejoined, 34 `A` and
/// 34 `-` rows at the record's `21:42029387-42029388`, every one `MINIMISED=1`
/// with the first feature's `FREQS=AF:-:0.002995`. Without `--buffer_size` the
/// eleven lines are one batch and the drop resets `GTG G,GAG` too: its `A` rows
/// stand at their own `21:42029387` with `ALLELE_NUM=2`, no `MINIMISED` and no
/// `FREQS`. Ensembl VEP 116.2's buffer under the `--buffer_size 10` flags holds
/// the ten chromosome-21 records, the skipped line never having entered it, so
/// its drop resets `GTG G,GAG` as well and it prints the `A` rows at
/// `21:42029387` with `ALLELE_NUM=2` and no `FREQS` where vep-rs rejoins them,
/// the shape the corpus's `expected_minimal-freqfilter-exclude` reference
/// prints for both records in one buffer; with `--fork 1` its chunk arithmetic
/// can place the record outside the drop's chunk, and it then rejoins.
#[test]
fn a_drops_reset_reaches_the_batch_counted_from_the_first_line() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reach.vcf");
    let snv = "21\t32288094\t.\tA\tT\t.\t.\t.\n";
    std::fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             22\t18110861\t.\tCAGCGATGGTG\tC\t.\t.\t.\n\
             21\t31148973\t.\tGG\tG,TG\t.\t.\t.\n\
             {}21\t42029386\t.\tGTG\tG,GAG\t.\t.\t.\n",
            snv.repeat(8)
        ),
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
        "--allele_number",
    ];
    let reset_t = ("21:31148973", "T", Some("2"), None, Some("AF:T:0.1396"));
    let whole_snv = (
        "21:32288094",
        "T",
        Some("1"),
        Some("1"),
        Some("AF:T:0.0003994"),
    );
    let span = "21:42029387-42029388";
    let rejoined_a = (span, "A", Some("2"), Some("1"), Some("AF:-:0.002995"));
    let reset_a = ("21:42029387", "A", Some("2"), None, None);
    let deleted_tg = (span, "-", Some("1"), Some("1"), Some("AF:-:0.002995"));

    let mut batched_flags = exclude.to_vec();
    batched_flags.extend(["--buffer_size", "10"]);
    let batched = rows(dir.path(), &input, &batched_flags);
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(
        warnings.contains("Chromosome 22 not found in annotation sources or synonyms"),
        "{warnings}"
    );
    assert_eq!(batched.len(), 124);
    assert_shapes(&batched, "21_31148973_GG/G/TG", &[(reset_t, 16)]);
    assert_shapes(&batched, "21_32288094_A/T", &[(whole_snv, 40)]);
    assert_shapes(
        &batched,
        "21_42029386_GTG/G/GAG",
        &[(rejoined_a, 34), (deleted_tg, 34)],
    );

    let one_batch = rows(dir.path(), &input, &exclude);
    assert_eq!(one_batch.len(), 124);
    assert_shapes(&one_batch, "21_31148973_GG/G/TG", &[(reset_t, 16)]);
    assert_shapes(&one_batch, "21_32288094_A/T", &[(whole_snv, 40)]);
    assert_shapes(
        &one_batch,
        "21_42029386_GTG/G/GAG",
        &[(reset_a, 34), (deleted_tg, 34)],
    );
}
