// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How far a `--check_frequency` drop reaches under `--minimal`: a dropped
//! record or a dropped split feature cancels the rejoin of the split records
//! around it, which the `vep` binary then writes as one record per surviving
//! feature, and the unit the drop resets is the dropped record's run of one
//! chromosome inside a `--buffer_size` batch, where Ensembl VEP 116.3 resets
//! its input buffer.

use std::path::{Path, PathBuf};
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

/// The committed GRCh38-freqfilter corpus.
fn freqfilter_corpus() -> common::Corpus {
    common::corpora()
        .into_iter()
        .find(|corpus| corpus.name == "GRCh38-freqfilter")
        .expect("the GRCh38-freqfilter corpus is committed")
}

/// A default-format run on `input` against the GRCh38-freqfilter corpus's
/// cache with `stats`, the run's statistics flag, and `extra` flags; returns
/// its output file, `out.txt` in `dir`.
fn run(dir: &Path, input: &Path, stats: &str, extra: &[&str]) -> PathBuf {
    let corpus = freqfilter_corpus();
    run_on(dir, &corpus.cache_args(), input, stats, extra)
}

/// [`run`] against the cache `cache_args` name, on the corpus's assembly.
fn run_on(dir: &Path, cache_args: &[String], input: &Path, stats: &str, extra: &[&str]) -> PathBuf {
    let out = dir.join("out.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .args(cache_args)
        .args([
            "--assembly",
            &freqfilter_corpus().assembly,
            "--force_overwrite",
            stats,
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
    out
}

/// The rows of a default-format run on `input` against the GRCh38-freqfilter
/// corpus's cache with `extra` flags.
fn rows(dir: &Path, input: &Path, extra: &[&str]) -> Vec<Row> {
    let out = run(dir, input, "--no_stats", extra);
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
/// 116.3 resets its input buffer (BaseCacheVariation.pm:366-369,
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
/// `FREQS`. Ensembl VEP 116.3's buffer under the `--buffer_size 10` flags holds
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

/// The `Variants filtered out` value and the data row count of a `--stats_text`
/// run on one VCF record with `extra` flags.
fn filtered_out(dir: &Path, record: &str, extra: &[&str]) -> (String, usize) {
    let input = dir.join("one.vcf");
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    let out = run(dir, &input, "--stats_text", extra);
    let rows = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .count();
    let summary = std::fs::read_to_string(dir.join("out.txt_summary.txt")).unwrap();
    let filtered = summary
        .lines()
        .find_map(|line| line.strip_prefix("Variants filtered out\t"))
        .unwrap_or_else(|| panic!("no Variants filtered out row: {summary}"))
        .to_string();
    (filtered, rows)
}

/// Ensembl VEP 116.3 counts every VariationFeature `frequency_check_buffer` leaves
/// out of its buffer (BaseCacheVariation.pm:363, `$filtered_count` is the buffer's
/// features less the ones passed, :374 `increment_filtered_variants`), which the
/// report prints as `Variants filtered out` (Stats.pm:723). Under
/// `--minimal` a multi-allelic record whose alleles trim is one feature per
/// alternate allele (InputBuffer.pm:638-718), a `*` or an allele equal to the
/// reference once trimmed among them, and such a feature, which writes no row,
/// matches no known allele, so `include` leaves it out
/// (BaseCacheVariation.pm:323) and counts it. On the release-116 GRCh38 cache
/// under `--check_frequency --freq_pop AF
/// --freq_freq 0.0001 --freq_gt_lt gt`: `21 31148973 . GG *,TG` prints `1` under
/// `include` (the `GG/*` feature; the `G/T` feature, AF 0.1396, stays and writes
/// 16 rows) and `1` under `exclude` (the `G/T` feature passes and goes, the
/// `GG/*` feature stays with no row); `GG TG,*` prints `1` under either;
/// `GG G,*,GG`, whose `*` and `GG` features write no row, prints `2` under
/// `include` (the `G/-` feature, AF 0.6705, stays) and `1` under `exclude` (it
/// goes); `GG *,GG`, whose two features are the `*` and `GG` features of
/// `GG G,*,GG` (the same chop, Parser/VCF.pm:306-320, and the same trim), prints
/// `2` under `include` and `0` under `exclude`, with no row under either.
/// `--process_ref_homs` adds a reference row per feature and changes no count.
/// Without `--minimal` the record is one feature and prints `1` with no row
/// under `include`.
#[test]
fn a_dropped_feature_no_row_stands_for_is_a_filtered_variant() {
    let dir = tempfile::tempdir().unwrap();
    let flags = |filter: &'static str, more: &[&'static str]| -> Vec<&'static str> {
        let mut flags = vec![
            "--minimal",
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_freq",
            "0.0001",
            "--freq_gt_lt",
            "gt",
            "--freq_filter",
            filter,
        ];
        flags.extend_from_slice(more);
        flags
    };
    let star_first = "21\t31148973\t.\tGG\t*,TG\t.\t.\t.";
    let star_last = "21\t31148973\t.\tGG\tTG,*\t.\t.\t.";
    let two_rowless = "21\t31148973\t.\tGG\tG,*,GG\t.\t.\t.";
    let all_rowless = "21\t31148973\t.\tGG\t*,GG\t.\t.\t.";

    let include = flags("include", &[]);
    assert_eq!(
        filtered_out(dir.path(), star_first, &include),
        ("1".to_string(), 16)
    );
    assert_eq!(
        filtered_out(dir.path(), star_last, &include),
        ("1".to_string(), 16)
    );
    assert_eq!(
        filtered_out(dir.path(), two_rowless, &include),
        ("2".to_string(), 16)
    );
    assert_eq!(
        filtered_out(dir.path(), all_rowless, &include),
        ("2".to_string(), 0)
    );
    let with_reference_rows = flags("include", &["--process_ref_homs"]);
    assert_eq!(
        filtered_out(dir.path(), star_first, &with_reference_rows),
        ("1".to_string(), 32)
    );
    assert_eq!(
        filtered_out(dir.path(), star_last, &with_reference_rows),
        ("1".to_string(), 32)
    );
    assert_eq!(
        filtered_out(dir.path(), two_rowless, &with_reference_rows),
        ("2".to_string(), 32)
    );
    assert_eq!(
        filtered_out(dir.path(), all_rowless, &with_reference_rows),
        ("2".to_string(), 0)
    );

    let exclude = flags("exclude", &[]);
    assert_eq!(
        filtered_out(dir.path(), star_first, &exclude),
        ("1".to_string(), 0)
    );
    assert_eq!(
        filtered_out(dir.path(), star_last, &exclude),
        ("1".to_string(), 0)
    );
    assert_eq!(
        filtered_out(dir.path(), two_rowless, &exclude),
        ("1".to_string(), 0)
    );
    assert_eq!(
        filtered_out(dir.path(), all_rowless, &exclude),
        ("0".to_string(), 0)
    );

    let without_minimal: Vec<&str> = include
        .iter()
        .copied()
        .filter(|flag| *flag != "--minimal")
        .collect();
    assert_eq!(
        filtered_out(dir.path(), star_first, &without_minimal),
        ("1".to_string(), 0)
    );
}

/// The data lines of a `--stats_text` run and the rows of its report.
struct Written {
    lines: Vec<String>,
    report: String,
}

impl Written {
    /// The value of a general statistics row, `None` for a label printed alone.
    fn general(&self, label: &str) -> Option<&str> {
        self.report
            .lines()
            .find_map(|line| match line.strip_prefix(label) {
                Some("") => Some(None),
                Some(rest) => rest.strip_prefix('\t').map(Some),
                None => None,
            })
            .unwrap_or_else(|| panic!("no {label} row: {}", self.report))
    }

    /// The rows of the chart headed `[title]`, in the report's order.
    fn chart(&self, title: &str) -> Vec<&str> {
        self.report
            .lines()
            .skip_while(|line| *line != format!("[{title}]"))
            .skip(1)
            .take_while(|line| !line.is_empty())
            .collect()
    }
}

/// A `--stats_text` run on the VCF `records` against the cache `cache_args`
/// name, with `extra` flags.
fn written(dir: &Path, cache_args: &[String], records: &[&str], extra: &[&str]) -> Written {
    let input = dir.join("records.vcf");
    std::fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
            records.join("\n")
        ),
    )
    .unwrap();
    let out = run_on(dir, cache_args, &input, "--stats_text", extra);
    Written {
        lines: std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_string)
            .collect(),
        report: std::fs::read_to_string(dir.join("out.txt_summary.txt")).unwrap(),
    }
}

/// The `--minimal --check_frequency` flags on the cache's `AF` column, `gt`,
/// at `threshold` under `filter`, then `more`.
fn frequency_flags<'a>(filter: &'a str, threshold: &'a str, more: &[&'a str]) -> Vec<&'a str> {
    let mut flags = vec![
        "--minimal",
        "--check_frequency",
        "--freq_pop",
        "AF",
        "--freq_freq",
        threshold,
        "--freq_gt_lt",
        "gt",
        "--freq_filter",
        filter,
    ];
    flags.extend_from_slice(more);
    flags
}

/// Ensembl VEP 116.3 writes every feature a reset buffer still holds
/// (OutputFactory.pm:301 maps `get_all_output_hashes_by_VariationFeature` over
/// the buffer, :342-349) and logs each one (:1022, Stats.pm:235 `var_count`,
/// :241 the most severe term, :249 the class), a row following or not: a split
/// feature whose allele is a `*`, or equal to the reference once trimmed
/// (VariationFeatureOverlap.pm:481, :484), is a record of its own with no allele
/// to annotate, so its VCF line carries no `CSQ`, its JSON object no
/// consequence, and the report counts it among `Variants processed` under the
/// class of its own pair with `intergenic_variant` as its most severe term
/// (VariationFeature.pm:1008-1013).
///
/// On the release-116 GRCh38 cache under `--check_frequency --freq_pop AF
/// --freq_gt_lt gt --freq_filter exclude`, traced on Ensembl VEP 116.3:
/// `21 31148973 . GG *,TG` at 0.0001 loses its `G/T` feature (AF 0.1396) and
/// keeps `GG/*`, one VCF line equal to the input, `Variants processed 1` and
/// `Variants filtered out 1`; `GG TG,*` reads the same. `GG G,*,GG` loses `G/-`
/// (AF 0.6705) and keeps `G/*` and `-/-`, two lines and two processed variants.
/// `GG *,G,TG` at 0.5 loses `G/-` alone: the `GG/*` line comes first, then the
/// line of `G/T` with its sixteen `T` rows, and `GG G,TG,*` writes the two in
/// the other order, the features standing in the alleles' order. At 0.5
/// `GG *,TG` loses nothing, so the rejoin makes one record of it, one line and
/// one processed variant. The default format writes no row for the feature and
/// counts it all the same.
#[test]
fn a_kept_feature_no_row_stands_for_is_a_record_of_its_own_after_a_drop() {
    let dir = tempfile::tempdir().unwrap();
    let cache = freqfilter_corpus().cache_args();
    let star_first = "21\t31148973\t.\tGG\t*,TG\t.\t.\t.";
    let star_last = "21\t31148973\t.\tGG\tTG,*\t.\t.\t.";
    let two_rowless = "21\t31148973\t.\tGG\tG,*,GG\t.\t.\t.";
    let star_before_kept = "21\t31148973\t.\tGG\t*,G,TG\t.\t.\t.";
    let star_after_kept = "21\t31148973\t.\tGG\tG,TG,*\t.\t.\t.";
    let vcf = |filter: &'static str, threshold: &'static str| {
        frequency_flags(filter, threshold, &["--vcf"])
    };

    for record in [star_first, star_last] {
        let run = written(dir.path(), &cache, &[record], &vcf("exclude", "0.0001"));
        assert_eq!(run.lines, [record]);
        assert_eq!(run.general("Variants processed"), Some("1"));
        assert_eq!(run.general("Variants filtered out"), Some("1"));
        assert_eq!(run.chart("Variant classes"), ["sequence_alteration\t1"]);
        assert_eq!(
            run.chart("Consequences (most severe)"),
            ["intergenic_variant\t1"]
        );
        assert_eq!(run.chart("Variants by chromosome"), ["21\t1"]);
    }

    let run = written(
        dir.path(),
        &cache,
        &[two_rowless],
        &vcf("exclude", "0.0001"),
    );
    assert_eq!(run.lines, [two_rowless, two_rowless]);
    assert_eq!(run.general("Variants processed"), Some("2"));
    assert_eq!(run.general("Variants filtered out"), Some("1"));
    assert_eq!(run.chart("Variant classes"), ["sequence_alteration\t2"]);

    let annotated = |line: &str, record: &str| {
        let fixed = &record[..record.len() - 1];
        line.strip_prefix(fixed)
            .is_some_and(|info| info.starts_with("CSQ=T|intron_variant|"))
    };
    let run = written(
        dir.path(),
        &cache,
        &[star_before_kept],
        &vcf("exclude", "0.5"),
    );
    assert_eq!(run.lines.len(), 2, "{:?}", run.lines);
    assert_eq!(run.lines[0], star_before_kept);
    assert!(
        annotated(&run.lines[1], star_before_kept),
        "{:?}",
        run.lines
    );
    assert_eq!(run.general("Variants processed"), Some("2"));
    assert_eq!(run.general("Variants filtered out"), Some("1"));
    let run = written(
        dir.path(),
        &cache,
        &[star_after_kept],
        &vcf("exclude", "0.5"),
    );
    assert_eq!(run.lines.len(), 2, "{:?}", run.lines);
    assert!(annotated(&run.lines[0], star_after_kept), "{:?}", run.lines);
    assert_eq!(run.lines[1], star_after_kept);
    assert_eq!(run.general("Variants processed"), Some("2"));

    let run = written(dir.path(), &cache, &[star_first], &vcf("exclude", "0.5"));
    assert_eq!(run.lines.len(), 1, "{:?}", run.lines);
    assert!(annotated(&run.lines[0], star_first), "{:?}", run.lines);
    assert_eq!(run.general("Variants processed"), Some("1"));
    assert_eq!(run.general("Variants filtered out"), Some("0"));

    let run = written(
        dir.path(),
        &cache,
        &[star_first],
        &frequency_flags("exclude", "0.0001", &[]),
    );
    assert!(run.lines.is_empty(), "{:?}", run.lines);
    assert_eq!(run.general("Variants processed"), Some("1"));

    let objects = |record: &str| -> Vec<serde_json::Value> {
        written(
            dir.path(),
            &cache,
            &[record],
            &frequency_flags("exclude", "0.0001", &["--json"]),
        )
        .lines
        .iter()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect()
    };
    let shape = |object: &serde_json::Value| {
        (
            object["allele_string"].as_str().unwrap().to_string(),
            object["start"].as_u64().unwrap(),
            object["end"].as_u64().unwrap(),
            object["minimised"].as_u64(),
            object["most_severe_consequence"]
                .as_str()
                .unwrap()
                .to_string(),
            object.get("transcript_consequences").is_some(),
        )
    };
    let unknown = "?".to_string();
    assert_eq!(
        objects(star_first).iter().map(shape).collect::<Vec<_>>(),
        [(
            "GG/*".to_string(),
            31148973,
            31148974,
            Some(1),
            unknown.clone(),
            false
        )]
    );
    assert_eq!(
        objects(two_rowless).iter().map(shape).collect::<Vec<_>>(),
        [
            (
                "G/*".to_string(),
                31148974,
                31148974,
                None,
                unknown.clone(),
                false
            ),
            ("-/-".to_string(), 31148975, 31148974, None, unknown, false),
        ]
    );
}

/// A split record none of whose alleles writes a row is one record with no
/// allele to annotate, which a reset leaves as it is: one VCF line and one of the
/// `Variants processed`, counted once in every chart. Ensembl VEP 116.3 writes
/// each of its two features and logs each (OutputFactory.pm:301, :1022,
/// Stats.pm:235), two lines and `Variants processed 2`
/// (`docs/output/known-variants.md`). `21 31148973 . GG *,GG` beside
/// `21 32288094 . A T`, whose AF 0.0003994 `exclude` drops at 0.0001.
#[test]
fn a_record_no_allele_of_which_writes_a_row_is_one_processed_variant_after_a_drop() {
    let dir = tempfile::tempdir().unwrap();
    let cache = freqfilter_corpus().cache_args();
    let all_rowless = "21\t31148973\t.\tGG\t*,GG\t.\t.\t.";
    let dropped = "21\t32288094\t.\tA\tT\t.\t.\t.";
    let run = written(
        dir.path(),
        &cache,
        &[all_rowless, dropped],
        &frequency_flags("exclude", "0.0001", &["--vcf"]),
    );
    assert_eq!(run.lines, [all_rowless]);
    assert_eq!(run.general("Variants processed"), Some("1"));
    assert_eq!(run.general("Variants filtered out"), Some("1"));
    assert_eq!(run.chart("Variant classes"), ["sequence_alteration\t1"]);
    assert_eq!(
        run.chart("Consequences (most severe)"),
        ["intergenic_variant\t1"]
    );
    assert_eq!(run.chart("Variants by chromosome"), ["21\t1"]);
}

/// A copy of the GRCh38-freqfilter corpus's cache under `dir` whose chromosome
/// 21 known-variant file holds `known` alone, each `(name, start, end, allele
/// string, AF)` with every other column at its default, block-compressed with a
/// tabix index beside it; returns the arguments naming it.
fn cache_with_known_variants(dir: &Path, known: &[(&str, u64, u64, &str, &str)]) -> Vec<String> {
    use std::io::Write;

    use noodles::bgzf;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

    fn copy_dir(from: &Path, to: &Path) {
        std::fs::create_dir_all(to).unwrap();
        for entry in std::fs::read_dir(from).unwrap() {
            let entry = entry.unwrap();
            let target = to.join(entry.file_name());
            if entry.file_type().unwrap().is_dir() {
                copy_dir(&entry.path(), &target);
            } else {
                std::fs::copy(entry.path(), target).unwrap();
            }
        }
    }

    let cache = dir.join("json_cache");
    copy_dir(&freqfilter_corpus().cache_dir(), &cache);
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(cache.join("info.json")).unwrap()).unwrap();
    let columns: Vec<&str> = info["variation_cols"]
        .as_array()
        .unwrap()
        .iter()
        .map(|column| column.as_str().unwrap())
        .collect();
    let column = |name: &str| columns.iter().position(|c| *c == name).unwrap();
    let variations = cache.join("variations").join("21");
    std::fs::remove_file(variations.join("all_vars.gz.csi")).unwrap();
    let path = variations.join("all_vars.gz");
    let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(
        HeaderBuilder::default()
            .set_reference_sequence_name_index(column("chr"))
            .set_start_position_index(column("start"))
            .set_end_position_index(Some(column("end")))
            .build(),
    );
    for (name, start, end, alleles, af) in known {
        let mut cells = vec![".".to_string(); columns.len()];
        cells[column("chr")] = "21".to_string();
        cells[column("variation_name")] = name.to_string();
        cells[column("start")] = start.to_string();
        cells[column("end")] = end.to_string();
        cells[column("allele_string")] = alleles.to_string();
        cells[column("AF")] = af.to_string();
        let from = writer.virtual_position();
        writeln!(writer, "{}", cells.join("\t")).unwrap();
        indexer
            .add_record(
                "21",
                Position::try_from(*start.min(end) as usize).unwrap(),
                Position::try_from(*start.max(end) as usize).unwrap(),
                Chunk::new(from, writer.virtual_position()),
            )
            .unwrap();
    }
    writer.finish().unwrap();
    let mut index_path = path.into_os_string();
    index_path.push(".tbi");
    noodles::tabix::fs::write(index_path, &indexer.build()).unwrap();
    vep_core::cache_manifest::CacheManifest::write(&cache).unwrap();
    vec!["--json_cache".to_string(), cache.display().to_string()]
}

/// A split feature that writes no row is judged like any other, on the known
/// variants of its own span (BaseCacheVariation.pm:316 at release/116.3). An
/// allele equal to the reference trims to `-/-` at the point after the record's
/// last base (InputBuffer.pm:665), where a known insertion lies: under
/// `--no_check_alleles` that insertion is listed on its start and end alone
/// (AnnotationType/Variation.pm:159-161), its reference `-`, which its frequency
/// text does not name, takes the remainder `1 - total`
/// (BaseCacheVariation.pm:249-252), and the feature's `-`, matched to nothing,
/// reads its own key (:431). Without the flag the pair `-/-` matches no allele
/// of the insertion and the feature has no frequency data.
///
/// The known-variant file here is made up: `rs200166967` as the corpus holds it
/// (`G/T` at 21:31148973, `AF` `T:0.1396`) and an insertion `-/T` after that
/// base with `AF` `T:0.25`, so `-` reads 0.75. `21 31148973 . G T,G` at 0.5:
/// under `exclude` with `--no_check_alleles` the `-/-` feature passes and goes,
/// which the report counts and which resets the buffer, so `G/T` is written on
/// its own under its own pair; under `include` `-/-` passes and stays, the one
/// record written, with the insertion as its known variant, while `G/T`
/// (0.1396) goes. Without `--no_check_alleles` the feature stays under
/// `exclude`, where the rejoined record keeps its string `G/T/G`, and goes
/// under `include` with `G/T`. Ensembl VEP 116.3 writes the four shapes on the
/// release-116 cache for a record built the same way at `rs545643089`, an
/// insertion after 21:30058548 whose `AF` holds 0.0001997.
#[test]
fn a_later_feature_no_row_stands_for_reads_the_known_insertion_at_its_point() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_known_variants(
        dir.path(),
        &[
            ("rs200166967", 31148973, 31148973, "G/T", "T:0.1396"),
            ("made_up_insertion", 31148974, 31148973, "-/T", "T:0.25"),
        ],
    );
    let record = "21\t31148973\t.\tG\tT,G\t.\t.\t.";
    let run = |filter: &'static str, more: &[&'static str]| {
        let mut flags = frequency_flags(filter, "0.5", &["--json"]);
        flags.extend_from_slice(more);
        let run = written(dir.path(), &cache, &[record], &flags);
        let objects: Vec<(String, u64, u64, Vec<String>)> = run
            .lines
            .iter()
            .map(|line| {
                let object: serde_json::Value = serde_json::from_str(line).unwrap();
                let known = object["colocated_variants"]
                    .as_array()
                    .map(|known| {
                        known
                            .iter()
                            .map(|variant| variant["id"].as_str().unwrap().to_string())
                            .collect()
                    })
                    .unwrap_or_default();
                (
                    object["allele_string"].as_str().unwrap().to_string(),
                    object["start"].as_u64().unwrap(),
                    object["end"].as_u64().unwrap(),
                    known,
                )
            })
            .collect();
        (objects, run)
    };
    let snv = vec!["rs200166967".to_string()];

    let (objects, report) = run("exclude", &["--no_check_alleles"]);
    assert_eq!(
        objects,
        [("G/T".to_string(), 31148973, 31148973, snv.clone())]
    );
    assert_eq!(report.general("Variants filtered out"), Some("1"));
    assert_eq!(report.general("Variants processed"), Some("1"));

    let (objects, report) = run("include", &["--no_check_alleles"]);
    assert_eq!(
        objects,
        [(
            "-/-".to_string(),
            31148974,
            31148973,
            vec!["made_up_insertion".to_string()]
        )]
    );
    assert_eq!(report.general("Variants filtered out"), Some("1"));
    assert_eq!(report.general("Variants processed"), Some("1"));
    assert_eq!(report.chart("Variant classes"), ["sequence_alteration\t1"]);

    let (objects, report) = run("exclude", &[]);
    assert_eq!(objects, [("G/T/G".to_string(), 31148973, 31148973, snv)]);
    assert_eq!(report.general("Variants filtered out"), Some("0"));
    assert_eq!(report.general("Variants processed"), Some("1"));

    let (objects, report) = run("include", &[]);
    assert!(objects.is_empty(), "{objects:?}");
    assert_eq!(report.general("Variants filtered out"), Some("2"));
    assert_eq!(report.general("Variants processed"), None);
}
