// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads SPDI input: a notation annotates as its VCF twin
//! does under the chromosome name the cache uses, the reference allele is
//! checked against or read from `--fasta`, which the format needs, a notation
//! Ensembl VEP's adaptor refuses is warned about and skipped, and `--format
//! guess` detects the format.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The identifiers corpus: a chromosome 21 cache slice, its synonyms naming
/// `NC_000021.9` and `chr21`, and a reference FASTA over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Five changes as VCF records, each named by the SPDI notation that spells it.
const VCF_RECORDS: &str = "\
21\t10538885\tNC_000021.9:10538884:A:C\tA\tC\t.\t.\t.
21\t10521699\tchr21:10521699:C:\tGC\tG\t.\t.\t.
21\t10451779\t21:10451779::TA\tA\tATA\t.\t.\t.
21\t10542415\tNC_000021.9:10542414:2:GG\tAA\tGG\t.\t.\t.
21\t13630971\tNC_000021.9:13630970:AATCCTGTGGGATTACCAGAAAACC:A\tAATCCTGTGGGATTACCAGAAAACC\tA\t.\t.\t.
";

/// The same changes as SPDI: a substitution on the RefSeq accession, a
/// deletion by its base on `chr21` with an empty insertion, an insertion by
/// its bases with an empty deletion, a substitution by length, and an
/// un-minimised deletion written with its anchor base as the VCF writes it.
const SPDI_RECORDS: &str = "\
NC_000021.9:10538884:A:C
chr21:10521699:C:
21:10451779::TA
NC_000021.9:10542414:2:GG
NC_000021.9:13630970:AATCCTGTGGGATTACCAGAAAACC:A
";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// The binary's arguments after the input: the corpus cache and the output
/// path, with `extra` last. The corpus FASTA is not among them, so a run can go
/// without one; `vep` adds it.
fn args(out: &Path, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        corpus().join("json_cache").display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
        "--variant_class".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

/// A run on `input` with the corpus FASTA and the [`args`] list for `out` and
/// `extra`.
fn vep(dir: &Path, input: &Path, out: &Path, extra: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args(args(out, extra))
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The output's lines that are not `#` headers.
fn rows(out: &Path) -> Vec<String> {
    String::from_utf8_lossy(&std::fs::read(out).unwrap())
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The text of the warnings file a run left beside its output `out`.
fn warnings(out: &Path) -> String {
    std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap()
}

/// An output path no earlier run of the process wrote.
fn output_path(dir: &Path) -> PathBuf {
    static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    dir.join(format!("out-{n}.txt"))
}

/// The rows of a successful run on `input` with `extra` flags.
fn annotate(dir: &Path, input: &Path, extra: &[&str]) -> Vec<String> {
    let out = output_path(dir);
    let output = vep(dir, input, &out, extra);
    assert!(
        output.status.success(),
        "{} {extra:?}: {}",
        input.display(),
        stderr(&output)
    );
    rows(&out)
}

/// The rows of the records as a VCF: what every SPDI spelling of them is
/// compared with. The VCF names each record by the notation, so the rows
/// agree column for column.
fn vcf_rows(dir: &Path, extra: &[&str]) -> Vec<String> {
    let vcf = write(dir, "records.vcf", &format!("{VCF_HEADER}{VCF_RECORDS}"));
    let rows = annotate(dir, &vcf, extra);
    assert!(rows.len() >= 5, "the VCF annotates: {rows:?}");
    rows
}

/// An SPDI record annotates as its VCF twin does in the default and tab
/// formats, row for row: the position is 0-based (VariationFeatureAdaptor.pm:
/// 2524), a deletion written with its anchor base is minimised as Ensembl VEP's
/// `post_process_vfs` minimises it (Parser.pm:862-881), a deletion or
/// substitution by length reads its bases from the FASTA (`:2563-2564,
/// 2589-2590`), and a record on `NC_000021.9` or `chr21` prints under the
/// cache's name `21` as Ensembl VEP prints the slice's (Parser/SPDI.pm:196).
/// The JSON object carries the notation as its input and its id, and the VCF
/// output rebuilds the record line from the variant.
#[test]
fn spdi_records_annotate_as_their_vcf_twins() {
    let dir = tempfile::tempdir().unwrap();
    let spdi = write(dir.path(), "records.spdi.txt", SPDI_RECORDS);
    let default_rows = annotate(dir.path(), &spdi, &["--format", "spdi"]);
    assert_eq!(default_rows, vcf_rows(dir.path(), &[]));
    assert_eq!(
        annotate(dir.path(), &spdi, &["--format", "spdi", "--tab"]),
        vcf_rows(dir.path(), &["--tab"])
    );
    assert!(
        default_rows
            .iter()
            .all(|row| row.split('\t').nth(1).is_some_and(|l| l.starts_with("21:"))),
        "every Location prints the cache's name: {default_rows:?}"
    );
    assert!(
        default_rows.iter().any(|row| row.starts_with(
            "NC_000021.9:13630970:AATCCTGTGGGATTACCAGAAAACC:A\t21:13630972-13630995\t-\t"
        )),
        "the anchored deletion is minimised: {default_rows:?}"
    );
    let json = annotate(dir.path(), &spdi, &["--format", "spdi", "--json"]);
    assert_eq!(json.len(), vcf_rows(dir.path(), &["--json"]).len());
    assert!(
        json[0].contains(r#""input":"NC_000021.9:10538884:A:C""#)
            && json[0].contains(r#""id":"NC_000021.9:10538884:A:C""#)
            && json[0].contains(r#""seq_region_name":"21""#),
        "the JSON object carries the notation: {}",
        json[0]
    );
    let vcf = annotate(dir.path(), &spdi, &["--format", "spdi", "--vcf"]);
    assert_eq!(vcf.len(), vcf_rows(dir.path(), &["--vcf"]).len());
    assert!(
        vcf[0].starts_with("21\t10538885\tNC_000021.9:10538884:A:C\tA\tC\t.\t.\tCSQ="),
        "{}",
        vcf[0]
    );
}

/// Without `--fasta` the format is refused before any line is read, as the
/// region format is: every notation's reference allele is checked against or
/// read from the sequence (VariationFeatureAdaptor.pm:2557-2590), which
/// Ensembl VEP takes from its database.
#[test]
fn spdi_without_fasta_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let spdi = write(dir.path(), "records.spdi.txt", SPDI_RECORDS);
    let out = output_path(dir.path());
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&spdi)
        .args(args(&out, &["--format", "spdi", "--no_fasta"]))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
    assert_eq!(
        stderr(&output).trim_end(),
        "ERROR: [input] Cannot use SPDI format in offline mode without a FASTA file"
    );
}

/// A notation the adaptor refuses is warned about with Ensembl VEP's own text
/// (Parser/SPDI.pm:183 and the `throw` it quotes) and skipped, the run going
/// on to the next line; the warning quotes the record's line, as every skip of
/// a parser-fed record does. The refusals: too few elements, a reference
/// allele the FASTA does not carry, and a sequence name nothing resolves. A
/// `#` line is such a record, as ensembl-io's list parser has no comment lines.
#[test]
fn unparsable_notation_is_warned_and_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "mixed.spdi.txt",
        "# header\nNC_000021.9:10538884:A\nNC_000021.9:10538884:G:C\nNC_000099.1:100:A:C\nNC_000021.9:10538884:A:C\n",
    );
    let out = output_path(dir.path());
    let output = vep(dir.path(), &input, &out, &["--format", "spdi"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let rows = rows(&out);
    assert!(
        !rows.is_empty()
            && rows
                .iter()
                .all(|r| r.starts_with("NC_000021.9:10538884:A:C\t")),
        "{rows:?}"
    );
    let warnings = warnings(&out);
    for expected in [
        "WARNING: line 1 skipped (# header): Unable to parse SPDI notation '#header': Could not parse the SPDI notation #header. Too few elements present",
        "WARNING: line 2 skipped (NC_000021.9:10538884:A): Unable to parse SPDI notation 'NC_000021.9:10538884:A': Could not parse the SPDI notation NC_000021.9:10538884:A. Too few elements present",
        "WARNING: line 3 skipped (NC_000021.9:10538884:G:C): Unable to parse SPDI notation 'NC_000021.9:10538884:G:C': Reference allele extracted from NC_000021.9:10538885-10538885 (A) does not match reference allele given by SPDI notation NC_000021.9:10538884:G:C (G)",
        "WARNING: line 4 skipped (NC_000099.1:100:A:C): Unable to parse SPDI notation 'NC_000099.1:100:A:C': Sequence name NC_000099.1 not valid",
    ] {
        assert!(warnings.contains(expected), "{expected}\n--\n{warnings}");
    }
}

/// `--format guess` detects a single `sequence:position:deletion:insertion`
/// field as SPDI (Utils/VEP.pm:407), so the file reads with no `--format`.
#[test]
fn format_guess_detects_spdi() {
    let dir = tempfile::tempdir().unwrap();
    let spdi = write(dir.path(), "records.spdi.txt", SPDI_RECORDS);
    assert_eq!(
        annotate(dir.path(), &spdi, &[]),
        annotate(dir.path(), &spdi, &["--format", "spdi"])
    );
}

/// Position 0 names the first base (VariationFeatureAdaptor.pm:2524), so
/// `NC_000021.9:0::G` is an insertion before it: start 1, end 0, printed
/// `21:0-1` with no transcript to overlap, as Ensembl VEP 116.2 prints it; a
/// deletion there is checked against base 1.
#[test]
fn insertion_before_the_first_base_prints_0_1() {
    let dir = tempfile::tempdir().unwrap();
    let spdi = write(
        dir.path(),
        "first.spdi.txt",
        "NC_000021.9:0::G\nNC_000021.9:0:N:G\n",
    );
    let rows = annotate(dir.path(), &spdi, &["--format", "spdi"]);
    assert_eq!(rows.len(), 2, "{rows:?}");
    assert!(
        rows[0].starts_with("NC_000021.9:0::G\t21:0-1\tG\t-\t-\t-\tintergenic_variant\t"),
        "{}",
        rows[0]
    );
    assert!(
        rows[1].starts_with("NC_000021.9:0:N:G\t21:1\tG\t-\t-\t-\tintergenic_variant\t"),
        "{}",
        rows[1]
    );
}

/// A span past the end of the sequence is skipped with a warning naming the
/// span and the sequence's length, where Ensembl VEP's adaptor builds a slice
/// past the end with no bounds test (VariationFeatureAdaptor.pm:2533,
/// 2545-2549, 2563-2564), reads `N` for every base past it
/// (Utils/FastaSequence.pm:619-629) and annotates the record: Ensembl VEP 116.2
/// prints `21:46709984-46709985 A intergenic_variant` and `21:46709984 -
/// intergenic_variant` for the VCF twins of `NC_000021.9:46709984::A` and
/// `NC_000021.9:46709983:1:` on the 46,709,983-base chromosome 21. The
/// insertion after the last base, `NC_000021.9:46709983::A`, reads no base past
/// the end and is annotated as Ensembl VEP annotates it; a name no sequence
/// carries, the version-less `NC_000021` among them, keeps `Sequence name
/// <name> not valid` (VariationFeatureAdaptor.pm:2536).
#[test]
fn span_past_the_sequence_end_is_warned_and_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "end.spdi.txt",
        "NC_000021.9:46709983::A\nNC_000021.9:46709984::A\nNC_000021.9:46709983:1:\nNC_000021.9:46709980:10:\nNC_000021:10538884:A:C\nNC_000099.1:100:A:C\n",
    );
    let out = output_path(dir.path());
    let output = vep(dir.path(), &input, &out, &["--format", "spdi"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let rows = rows(&out);
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert!(
        rows[0].starts_with(
            "NC_000021.9:46709983::A\t21:46709983-46709984\tA\t-\t-\t-\tintergenic_variant\t"
        ),
        "{}",
        rows[0]
    );
    assert_eq!(
        warnings(&out).lines().collect::<Vec<_>>(),
        [
            "WARNING: line 2 skipped (NC_000021.9:46709984::A): Unable to parse SPDI notation 'NC_000021.9:46709984::A': Span 46709984-46709985 lies past the end of sequence NC_000021.9 (46709983 bases)",
            "WARNING: line 3 skipped (NC_000021.9:46709983:1:): Unable to parse SPDI notation 'NC_000021.9:46709983:1:': Span 46709984-46709984 lies past the end of sequence NC_000021.9 (46709983 bases)",
            "WARNING: line 4 skipped (NC_000021.9:46709980:10:): Unable to parse SPDI notation 'NC_000021.9:46709980:10:': Span 46709981-46709990 lies past the end of sequence NC_000021.9 (46709983 bases)",
            "WARNING: line 5 skipped (NC_000021:10538884:A:C): Unable to parse SPDI notation 'NC_000021:10538884:A:C': Sequence name NC_000021 not valid",
            "WARNING: line 6 skipped (NC_000099.1:100:A:C): Unable to parse SPDI notation 'NC_000099.1:100:A:C': Sequence name NC_000099.1 not valid",
        ]
    );
}

/// A deletion of length `0` is skipped with a warning naming the empty span,
/// where Ensembl VEP's adaptor reads it as a deletion ending at the position
/// (VariationFeatureAdaptor.pm:2552, :2563), takes the empty slice's sequence as
/// its reference allele (:2564) and keeps the feature `/-` (:2566), which
/// writes no consequence row: Ensembl VEP 116.2 writes no default or tab row
/// for `21:10538884:0:`, the VCF line `21 10538885 21:10538884:0: <empty> . . .`
/// without consequences and a JSON object with `allele_string` `/-` and
/// `most_severe_consequence` `?`. The notation with `0` for both parts is
/// refused by both engines (:2540-2542); the run goes on after each.
#[test]
fn a_zero_length_deletion_is_warned_and_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "zero.spdi.txt",
        "21:10538884:0:\nNC_000021.9:46709983::A\n21:10538884:00:\n21:10538884:0:0\n",
    );
    let out = output_path(dir.path());
    let output = vep(dir.path(), &input, &out, &["--format", "spdi"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let rows = rows(&out);
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert!(
        rows[0].starts_with("NC_000021.9:46709983::A\t21:46709983-46709984\tA\t"),
        "{}",
        rows[0]
    );
    assert_eq!(
        warnings(&out).lines().collect::<Vec<_>>(),
        [
            "WARNING: line 1 skipped (21:10538884:0:): Unable to parse SPDI notation '21:10538884:0:': Deleted sequence length (0) leaves nothing to delete: span 10538885-10538884 is empty",
            "WARNING: line 3 skipped (21:10538884:00:): Unable to parse SPDI notation '21:10538884:00:': Deleted sequence length (00) leaves nothing to delete: span 10538885-10538884 is empty",
            "WARNING: line 4 skipped (21:10538884:0:0): Unable to parse SPDI notation '21:10538884:0:0': Could not parse the SPDI notation 21:10538884:0:0. SPDI notation not supported.",
        ]
    );
}

/// A record whose span exceeds `--max_sv_size` is skipped with the text
/// Ensembl VEP's `validate_vf` skips it with (Parser.pm:491-499), the span being
/// end minus start: under a limit of 1000 a deletion of 1002 bases (span 1001)
/// is skipped and one of 1001 bases (span 1000) annotated, so a deletion
/// written by its length never annotates more than the limit allows (its bases
/// are read before the span is judged, as Ensembl VEP reads its slice first).
#[test]
fn oversize_deletion_is_skipped_with_the_max_sv_size_warning() {
    let dir = tempfile::tempdir().unwrap();
    let spdi = write(
        dir.path(),
        "long.spdi.txt",
        "NC_000021.9:10538884:1002:\nNC_000021.9:10538884:1001:\n",
    );
    let out = output_path(dir.path());
    let output = vep(
        dir.path(),
        &spdi,
        &out,
        &["--format", "spdi", "--max_sv_size", "1000"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let rows = rows(&out);
    assert!(
        !rows.is_empty()
            && rows
                .iter()
                .all(|r| r.starts_with("NC_000021.9:10538884:1001:\t21:10538885-10539885\t-\t")),
        "{rows:?}"
    );
    let warnings = warnings(&out);
    assert!(
        warnings.contains(
            "WARNING: line 1 skipped (NC_000021.9:10538884:1002:): variant size (1001) is bigger than --max_sv_size (1000)"
        ),
        "{warnings}"
    );
}
