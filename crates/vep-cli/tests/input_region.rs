// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads region input: a region annotates as its VCF twin
//! does once `--fasta` fills its reference allele, a reverse-strand region reads
//! the reverse complement, a structural word names its class, a line outside
//! Ensembl VEP's region expression or over `--max_sv_size` is warned about and
//! skipped, `--format guess` detects the format, and the record is named by its
//! line as written.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The identifiers corpus: a chromosome 21 cache slice, its synonyms naming
/// `chr21`, and a reference FASTA over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Seven changes as VCF records, each named by the region that spells it: an
/// SNV, an SNV written with a colon before its allele, one with a strand token
/// and a colon, one on the `+1` strand with a lower-case allele, a two-base
/// substitution, a deletion written with its anchor base and an insertion.
const VCF_RECORDS: &str = "\
21\t10538885\t21:10538885-10538885:C\tA\tC\t.\t.\t.
21\t10542415\t21:10542415-10542416/GG\tAA\tGG\t.\t.\t.
21\t13615721\t21:13615721-13615721:1:C\tA\tC\t.\t.\t.
21\t13620283\t21:13620284-13620283/TCA\tG\tGTCA\t.\t.\t.
21\t13630971\t21:13630972-13630995/-\tAATCCTGTGGGATTACCAGAAAACC\tA\t.\t.\t.
21\t14221107\t21:14221107-14221107/A\tC\tA\t.\t.\t.
21\t15826987\t21:15826987-15826987:+1/g\tC\tG\t.\t.\t.
";

/// The same changes as regions.
const REGION_RECORDS: &str = "\
21:10538885-10538885:C
21:10542415-10542416/GG
21:13615721-13615721:1:C
21:13620284-13620283/TCA
21:13630972-13630995/-
21:14221107-14221107/A
21:15826987-15826987:+1/g
";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// The binary's arguments after the input: the corpus cache, the corpus FASTA
/// and the output path, with `extra` last.
fn args(out: &Path, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        corpus().join("json_cache").display().to_string(),
        "--fasta".into(),
        corpus().join("reference.fa.gz").display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--offline".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
        "--variant_class".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

fn vep(dir: &Path, input: &Path, out: &Path, extra: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
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

/// The text of the warnings file a run left beside its output `out`, empty
/// when it left none.
fn warnings(out: &Path) -> String {
    std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default()
}

/// An output path no earlier run of the process wrote.
fn output_path(dir: &Path) -> PathBuf {
    static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    dir.join(format!("out-{n}.txt"))
}

/// The rows of a successful run on `input` with `extra` flags, and the text of
/// the warnings file it left.
fn annotate_with_warnings(dir: &Path, input: &Path, extra: &[&str]) -> (Vec<String>, String) {
    let out = output_path(dir);
    let output = vep(dir, input, &out, extra);
    assert!(
        output.status.success(),
        "{} {extra:?}: {}",
        input.display(),
        stderr(&output)
    );
    (rows(&out), warnings(&out))
}

/// The rows of a successful run on `input` with `extra` flags.
fn annotate(dir: &Path, input: &Path, extra: &[&str]) -> Vec<String> {
    let (rows, _warnings) = annotate_with_warnings(dir, input, extra);
    rows
}

/// The rows of the records as a VCF, each named by its region, so the rows
/// agree column for column with the region run's.
fn vcf_rows(dir: &Path, extra: &[&str]) -> Vec<String> {
    let vcf = write(dir, "records.vcf", &format!("{VCF_HEADER}{VCF_RECORDS}"));
    let rows = annotate(dir, &vcf, extra);
    assert!(rows.len() >= 7, "the VCF annotates: {rows:?}");
    rows
}

/// The `(Uploaded_variation, UPLOADED_ALLELE)` pair of each record of a
/// default-format output, in order: the rows of one record repeat the pair,
/// so consecutive equal pairs collapse to one.
fn uploaded_alleles(rows: &[String]) -> Vec<(String, String)> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    for row in rows {
        let columns: Vec<&str> = row.split('\t').collect();
        let name = columns[0].to_string();
        let allele = columns
            .get(13)
            .and_then(|extra| {
                extra
                    .split(';')
                    .find_map(|field| field.strip_prefix("UPLOADED_ALLELE="))
            })
            .unwrap_or("")
            .to_string();
        let pair = (name, allele);
        if pairs.last() != Some(&pair) {
            pairs.push(pair);
        }
    }
    pairs
}

/// Whether `name` is the one record with rows: there is at least one row and
/// every row is its.
fn only_record_has_rows(rows: &[String], name: &str) -> bool {
    let prefix = format!("{name}\t");
    !rows.is_empty() && rows.iter().all(|row| row.starts_with(&prefix))
}

/// A region record annotates as its VCF twin does in the default and tab
/// formats, row for row: the reference allele is read from the FASTA over the
/// span (Parser/Region.pm:121-122, 189), an allele follows a `:` as it follows a
/// `/`, a strand token `:1` or `:+1` is forward, and a lower-case allele is
/// upper-cased (:176). The JSON object carries the region as its input and its
/// id.
#[test]
fn region_records_annotate_as_their_vcf_twins() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(dir.path(), "records.region.txt", REGION_RECORDS);
    let default_rows = annotate(dir.path(), &region, &["--format", "region"]);
    assert_eq!(default_rows, vcf_rows(dir.path(), &[]));
    assert_eq!(
        annotate(dir.path(), &region, &["--format", "region", "--tab"]),
        vcf_rows(dir.path(), &["--tab"])
    );
    let json = annotate(dir.path(), &region, &["--format", "region", "--json"]);
    assert_eq!(json.len(), 7, "{json:?}");
    assert!(
        json[0].contains(r#""input":"21:10538885-10538885:C""#)
            && json[0].contains(r#""id":"21:10538885-10538885:C""#)
            && json[0].contains(r#""allele_string":"A/C""#),
        "the JSON object carries the region and the looked-up allele string: {}",
        json[0]
    );
}

/// `UPLOADED_ALLELE` prints the allele string the reference lookup leaves,
/// `C/A` for `21:14221107-14221107/A`, except for a record whose alleles differ
/// in length with no `-` among them, where Ensembl VEP's minimisation keeps the
/// placeholder as parsed: Ensembl VEP 116.2 prints `UPLOADED_ALLELE=NN/GTA` for
/// `21:10552677-10552678/GTA`, `ATCCTGTGGGATTACCAGAAAACC/-` for
/// `21:13630972-13630995/-` and `-/TCA` for `21:13620284-13620283/TCA`.
#[test]
fn uploaded_allele_is_the_looked_up_string_except_for_an_unminimised_indel() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "alleles.region.txt",
        "21:10552677-10552678/GTA\n21:13620284-13620283/TCA\n21:13630972-13630995/-\n21:14221107-14221107/A\n",
    );
    let rows = annotate(
        dir.path(),
        &region,
        &["--format", "region", "--uploaded_allele"],
    );
    assert_eq!(
        uploaded_alleles(&rows),
        [
            ("21:10552677-10552678/GTA", "NN/GTA"),
            ("21:13620284-13620283/TCA", "-/TCA"),
            ("21:13630972-13630995/-", "ATCCTGTGGGATTACCAGAAAACC/-"),
            ("21:14221107-14221107/A", "C/A"),
        ]
        .map(|(name, allele)| (name.to_string(), allele.to_string()))
    );
}

/// Under `--minimal` Ensembl VEP minimises every bi-allelic record before the
/// reference lookup (Parser.pm:852), keeping the parser's `N` placeholder as the
/// string uploaded whether or not a base is trimmed (:1008), so `UPLOADED_ALLELE`
/// prints the placeholder for every sequence record (OutputFactory.pm:1360):
/// Ensembl VEP 116.2 prints `N/A` for `21:14221107-14221107/A`, `N/-` for
/// `21:10521700-10521700/-`, `NN/CA` for the reverse-strand
/// `21:10542415-10542416:-1/CA`, `NNNNNNNNNNNNNNNNNNNNNNNN/-` for
/// `21:13630972-13630995/-`, and the `NN/GTA` and `-/TCA` it prints without the
/// flag, with `MINIMISED=1` on every row. `21:14221105-14221107/NA` prints its
/// placeholder `NNN/NA` here; Ensembl VEP 116.2 prints `NN/A`, the string its
/// second pass over the once-trimmed record leaves (`:862-881` hand `NN/A`, whose
/// alleles differ in length, to `minimise_alleles` again, and `:1008` resets
/// `original_allele_string` to it), which is not the allele uploaded.
#[test]
fn minimal_keeps_the_placeholder_as_the_uploaded_allele() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "minimal.region.txt",
        "21:10552677-10552678/GTA\n21:13620284-13620283/TCA\n21:13630972-13630995/-\n21:14221107-14221107/A\n21:10521700-10521700/-\n21:10542415-10542416:-1/CA\n21:14221105-14221107/NA\n",
    );
    let rows = annotate(
        dir.path(),
        &region,
        &["--format", "region", "--uploaded_allele", "--minimal"],
    );
    assert_eq!(
        uploaded_alleles(&rows),
        [
            ("21:10552677-10552678/GTA", "NN/GTA"),
            ("21:13620284-13620283/TCA", "-/TCA"),
            ("21:13630972-13630995/-", "NNNNNNNNNNNNNNNNNNNNNNNN/-"),
            ("21:14221107-14221107/A", "N/A"),
            ("21:10521700-10521700/-", "N/-"),
            ("21:10542415-10542416:-1/CA", "NN/CA"),
            ("21:14221105-14221107/NA", "NNN/NA"),
        ]
        .map(|(name, allele)| (name.to_string(), allele.to_string()))
    );
    assert!(
        rows.iter().all(|row| row.ends_with(";MINIMISED=1")),
        "{rows:?}"
    );
}

/// `--max_sv_size` measures the span `post_process_vfs` leaves (Parser.pm:852
/// before :492), and under `--minimal` every bi-allelic record is trimmed first,
/// its length and `-` tests set aside: `21:14221104-14221107/NNAA` is `NN/AA` over
/// 14221106-14221107, the span 1 that `--max_sv_size 2` admits, and Ensembl VEP
/// 116.2 annotates it at `21:14221106-14221107` with Allele `AA` and
/// `UPLOADED_ALLELE=NNNN/NNAA`. Without `--minimal` the same-length alleles are
/// not trimmed and the span 3 is refused.
#[test]
fn minimal_trims_a_same_length_allele_before_the_size_test() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(dir.path(), "same.region.txt", "21:14221104-14221107/NNAA\n");
    let flags = [
        "--format",
        "region",
        "--uploaded_allele",
        "--max_sv_size",
        "2",
    ];
    let (rows, warnings) =
        annotate_with_warnings(dir.path(), &region, &[&flags[..], &["--minimal"]].concat());
    assert!(!rows.is_empty(), "{warnings}");
    assert!(warnings.is_empty(), "{warnings}");
    for row in &rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns[1], "21:14221106-14221107", "{row}");
        assert_eq!(columns[2], "AA", "{row}");
        assert!(columns[13].contains("UPLOADED_ALLELE=NNNN/NNAA;"), "{row}");
    }

    let (rows, warnings) = annotate_with_warnings(dir.path(), &region, &flags);
    assert!(rows.is_empty(), "{rows:?}");
    assert!(
        warnings.contains(
            "WARNING: line 1 skipped (21:14221104-14221107/NNAA): variant size (3) is bigger than --max_sv_size (2)"
        ),
        "{warnings}"
    );
}

/// The region expression admits a `-` among an allele's bases, and
/// `validate_vf` removes it from every allele (Parser.pm:578-580) before the
/// `--lookup_ref` rejoin the format runs (`:645`), so `21:14221107-14221107/A-`
/// annotates as `21:14221107-14221107/A` does, row for row: Ensembl VEP 116.2
/// prints 17 rows with Allele `A`, `UPLOADED_ALLELE=C/A` and
/// `VARIANT_CLASS=SNV`, the ENST00000400577 row a `synonymous_variant` at
/// cDNA 300, CDS 270, protein 90, `R`, `cgC/cgA`.
#[test]
fn dash_inside_an_allele_is_removed_at_the_lookup_ref_rejoin() {
    let dir = tempfile::tempdir().unwrap();
    let dashed = write(dir.path(), "dashed.region.txt", "21:14221107-14221107/A-\n");
    let plain = write(dir.path(), "plain.region.txt", "21:14221107-14221107/A\n");
    let flags = ["--format", "region", "--uploaded_allele"];
    let rows = annotate(dir.path(), &dashed, &flags);
    assert_eq!(rows.len(), 17, "{rows:?}");
    assert!(
        rows.contains(
            &"21:14221107-14221107/A-\t21:14221107\tA\tENSG00000185272\tENST00000400577\tTranscript\tsynonymous_variant\t300\t270\t90\tR\tcgC/cgA\t-\tUPLOADED_ALLELE=C/A;IMPACT=LOW;STRAND=1;VARIANT_CLASS=SNV".to_string()
        ),
        "{rows:?}"
    );
    let renamed: Vec<String> = annotate(dir.path(), &plain, &flags)
        .iter()
        .map(|row| row.replacen("21:14221107-14221107/A\t", "21:14221107-14221107/A-\t", 1))
        .collect();
    assert_eq!(rows, renamed);
}

/// A region on the reverse strand (`:-1`) reads its reference allele as the
/// reverse complement of the FASTA's bases (`sub_Slice` on the record's strand,
/// Parser.pm:751) and its allele is read against a forward-strand transcript as
/// its reverse complement, so `21:14964727-14964727:-1/C` annotates as the
/// forward `C>G` at that base does with the Allele column `C` in place of `G`:
/// Ensembl VEP 116.2 prints `missense_variant`, `E/Q`, `Gaa/Caa` and
/// `UPLOADED_ALLELE=G/C;IMPACT=MODERATE;STRAND=-1` for it.
#[test]
fn reverse_strand_region_reads_the_reverse_complement() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "reverse.region.txt",
        "21:14964727-14964727:-1/C\n",
    );
    let rows = annotate(
        dir.path(),
        &region,
        &["--format", "region", "--uploaded_allele"],
    );
    let vcf = write(
        dir.path(),
        "forward.vcf",
        &format!("{VCF_HEADER}21\t14964727\t21:14964727-14964727:-1/C\tC\tG\t.\t.\t.\n"),
    );
    let forward = annotate(dir.path(), &vcf, &["--uploaded_allele"]);
    assert_eq!(rows.len(), forward.len(), "{rows:?}\n--\n{forward:?}");
    for (row, twin) in rows.iter().zip(&forward) {
        let mut columns: Vec<String> = row.split('\t').map(str::to_string).collect();
        let allele = &mut columns[2];
        assert_eq!(*allele, "C", "{row}");
        *allele = "G".to_string();
        let as_forward = columns
            .join("\t")
            .replace("UPLOADED_ALLELE=G/C", "UPLOADED_ALLELE=C/G");
        assert_eq!(as_forward, *twin);
    }
    assert!(
        rows.iter().any(|row| row.contains(
            "\tENST00000318948\tTranscript\tmissense_variant\t4259\t3466\t1156\tE/Q\tGaa/Caa\t"
        ) && row
            .contains("UPLOADED_ALLELE=G/C;IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV")),
        "{rows:?}"
    );
    let json = annotate(dir.path(), &region, &["--format", "region", "--json"]);
    assert!(
        json[0].contains(r#""allele_string":"G/C""#) && json[0].contains(r#""strand":-1"#),
        "{}",
        json[0]
    );
}

/// An allele outside `ACGTN-` names a structural class through `get_SO_term`
/// (Parser/Region.pm:205) over the span as written: `DUP` is a duplication,
/// `del` a deletion and `INS:ME:ALU` an Alu insertion, each printing its class
/// term as the Allele (Ensembl VEP 116.2 prints `21:14221000-14222000/DUP
/// 21:14221000-14222000 duplication ... feature_elongation,coding_sequence_variant,
/// intron_variant`). A word `%SO_TERMS` lacks is warned about with `get_SO_term`'s
/// text (Parser.pm:725) and skipped, the run going on.
#[test]
fn structural_words_annotate_as_their_class() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "structural.region.txt",
        "21:13609000-13610500/del\n21:14221000-14222000/DUP\n21:14221000-14222000/FOO\n21:14221107-14221107/INS:ME:ALU\n",
    );
    let (rows, warnings) = annotate_with_warnings(dir.path(), &region, &["--format", "region"]);
    let first = |name: &str| {
        rows.iter()
            .find(|row| row.starts_with(&format!("{name}\t")))
            .unwrap_or_else(|| panic!("{name} has a row: {rows:?}"))
            .clone()
    };
    let duplication = first("21:14221000-14222000/DUP");
    assert!(
        duplication.starts_with(
            "21:14221000-14222000/DUP\t21:14221000-14222000\tduplication\tENSG00000185272\tENST00000400577\tTranscript\tfeature_elongation,coding_sequence_variant,intron_variant\t"
        ),
        "{duplication}"
    );
    let deletion = first("21:13609000-13610500/del");
    assert!(
        deletion.contains("\t21:13609000-13610500\tdeletion\t")
            && deletion.contains("VARIANT_CLASS=deletion"),
        "{deletion}"
    );
    let alu_insertion = first("21:14221107-14221107/INS:ME:ALU");
    assert!(
        alu_insertion
            .contains("\t21:14221107\tAlu_insertion\tENSG00000185272\tENST00000400577\tTranscript\tfeature_elongation,coding_sequence_variant\t"),
        "{alu_insertion}"
    );
    assert!(
        !rows.iter().any(|row| row.contains("/FOO\t")),
        "the unsupported word has no row: {rows:?}"
    );
    assert!(
        warnings.contains(
            "WARNING: line 3 skipped (21:14221000-14222000/FOO): FOO is not a supported structural variant type"
        ),
        "{warnings}"
    );
}

/// A structural word `%SO_TERMS` lacks on the first line of the file is warned
/// about with `get_SO_term`'s text (Parser.pm:725) and skipped
/// (`RegionParser::parse_record` fails the line with `IoError::SkippedLine`
/// carrying that text, and the caller warns and skips it), and the lines after
/// it annotate: `21:14221107-14221107/A` writes its 17 rows, each named by its
/// line as written, exit 0. Ensembl VEP 116.2 writes the same warning and then
/// 0 rows for this input, exit 0: `create_VariationFeatures` builds no feature
/// for the word (Parser/Region.pm:205-206, `my $so_term =
/// $self->get_SO_term($allele); return [] unless defined $so_term;`), the
/// parser's `next` returns the empty list's undefined record (Parser.pm:186-187,
/// `my $vf = shift @$cache; return $vf unless $vf;`), the buffer fill ends on it
/// (InputBuffer.pm:212, `while(@$buffer < $buffer_size && (my $vf =
/// $parser->next))`) and the run ends (Runner.pm:424; :477 under `--fork`), the
/// second line unread.
#[test]
fn an_unsupported_structural_word_on_the_first_line_is_skipped_and_the_run_goes_on() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "first.region.txt",
        "21:14221000-14222000/FOO\n21:14221107-14221107/A\n",
    );
    let (rows, warnings) = annotate_with_warnings(dir.path(), &region, &["--format", "region"]);
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [
            "WARNING: line 1 skipped (21:14221000-14222000/FOO): FOO is not a supported \
             structural variant type",
        ],
        "{warnings}"
    );
    assert_eq!(rows.len(), 17, "{rows:?}");
    assert!(
        only_record_has_rows(&rows, "21:14221107-14221107/A"),
        "{rows:?}"
    );
}

/// A line outside Ensembl VEP's region expression (Utils/VEP.pm:391) builds no
/// feature in Ensembl VEP, which ends its buffer fill there and, when the line
/// opens a fill as this file's first line does, its run (InputBuffer.pm:212,
/// Runner.pm:424: Ensembl VEP 116.2 writes no row and no warning for this
/// file); here it is warned about and skipped, the run going on: a position
/// without its end, a trailing blank, a strand other than `1`, `+1` or `-1`,
/// and a multi-allelic string. A `#` line is a comment, counted into the line
/// numbers that follow as Ensembl VEP counts the record it reads there
/// (Parser/Region.pm:167, over ensembl-io's list parser, which has no metadata
/// lines, ListBasedParser.pm:60), so the line after it is line 2; an empty line
/// is warned about and counted too, as `skip_empty_lines` counts it
/// (Parser.pm:266-293, the record after one empty line being line 4). The same
/// lines are skipped under `--dont_skip`, which has no feature to keep.
#[test]
fn lines_outside_the_expression_are_warned_and_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "mixed.region.txt",
        "# a comment\n21:14221107/A\n\n21:14221107-14221107/A \n21:14221107-14221107:2/A\n21:14221107-14221107/A/G\n21:14221107-14221107/A\n",
    );
    for extra in [
        &["--format", "region"][..],
        &["--format", "region", "--dont_skip"][..],
    ] {
        let (rows, warnings) = annotate_with_warnings(dir.path(), &region, extra);
        assert!(
            only_record_has_rows(&rows, "21:14221107-14221107/A"),
            "{extra:?}: {rows:?}"
        );
        assert_eq!(
            warnings.lines().collect::<Vec<_>>(),
            [
                "WARNING: line 2 skipped (21:14221107/A): Unable to parse region '21:14221107/A': expected chr:start-end[:strand]/allele",
                "WARNING: Skipped 1 empty line from line 3",
                "WARNING: line 4 skipped (21:14221107-14221107/A ): Unable to parse region '21:14221107-14221107/A ': expected chr:start-end[:strand]/allele",
                "WARNING: line 5 skipped (21:14221107-14221107:2/A): Unable to parse region '21:14221107-14221107:2/A': expected chr:start-end[:strand]/allele",
                "WARNING: line 6 skipped (21:14221107-14221107/A/G): Unable to parse region '21:14221107-14221107/A/G': expected chr:start-end[:strand]/allele",
            ],
            "{extra:?}: {warnings}"
        );
    }
}

/// A `#` line with a word character counts toward the line numbers that
/// follow, as Ensembl VEP counts it: ensembl-io's list parser has no metadata
/// lines (ListBasedParser.pm:60, `sub is_metadata { return; }`), so the line is
/// a record `create_VariationFeatures` numbers (Parser/Region.pm:167,
/// `$self->line_number($self->line_number + 1);`) and then builds no feature for
/// (`:171`, `return [] unless $region =~ &_valid_region_regex();`), which ends
/// the buffer fill there and reads on with the next fill (InputBuffer.pm:212)
/// unless the line opened the fill; here it is passed over without a word. For
/// `21:14221107-14221107/A`, `# c`, `21:1-20000000/A`, Ensembl VEP 116.2 writes
/// `WARNING: line 3 skipped (21:1-20000000/A): variant size (19999999) is
/// bigger than --max_sv_size (10000000)` and the first record's 17 rows, as
/// here.
#[test]
fn a_hash_line_counts_toward_the_line_numbers() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "hash.region.txt",
        "21:14221107-14221107/A\n# c\n21:1-20000000/A\n",
    );
    let (rows, warnings) = annotate_with_warnings(dir.path(), &region, &["--format", "region"]);
    assert_eq!(rows.len(), 17, "{rows:?}");
    assert!(
        only_record_has_rows(&rows, "21:14221107-14221107/A"),
        "{rows:?}"
    );
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        ["WARNING: line 3 skipped (21:1-20000000/A): variant size (19999999) is bigger than --max_sv_size (10000000)"],
        "{warnings}"
    );
}

/// A `#` line with no word character is an empty line to `skip_empty_lines`
/// (Parser.pm:273, `$parser->{record} !~ /\w+/`), warned about and counted as
/// one: for `21:14221107-14221107/A`, `#`, `21:1-20000000/A`, Ensembl VEP 116.2
/// writes `WARNING: Skipped 1 empty line from line 2` and `WARNING: line 3
/// skipped (21:1-20000000/A): variant size (19999999) is bigger than
/// --max_sv_size (10000000)`. A `#` line with a word character before such a
/// run is the record that ends the run before it and starts the count of the
/// next: for `21:14221107-14221107/A`, `# c`, `#`, `  `, `21:1-20000000/A`,
/// Ensembl VEP 116.2 writes `WARNING: Skipped 2 empty lines from line 3` and
/// `WARNING: line 5 skipped (21:1-20000000/A): ...`. 17 rows each.
#[test]
fn a_bare_hash_line_is_an_empty_line() {
    let dir = tempfile::tempdir().unwrap();
    for (name, text, expected_warnings) in [
        (
            "bare.region.txt",
            "21:14221107-14221107/A\n#\n21:1-20000000/A\n",
            [
                "WARNING: Skipped 1 empty line from line 2",
                "WARNING: line 3 skipped (21:1-20000000/A): variant size (19999999) is bigger than --max_sv_size (10000000)",
            ],
        ),
        (
            "after-comment.region.txt",
            "21:14221107-14221107/A\n# c\n#\n  \n21:1-20000000/A\n",
            [
                "WARNING: Skipped 2 empty lines from line 3",
                "WARNING: line 5 skipped (21:1-20000000/A): variant size (19999999) is bigger than --max_sv_size (10000000)",
            ],
        ),
    ] {
        let region = write(dir.path(), name, text);
        let (rows, warnings) =
            annotate_with_warnings(dir.path(), &region, &["--format", "region"]);
        assert_eq!(rows.len(), 17, "{name}: {rows:?}");
        assert_eq!(
            warnings.lines().collect::<Vec<_>>(),
            expected_warnings,
            "{name}: {warnings}"
        );
    }
}

/// A refused line that kept a `\r` is quoted without it: a line ending `\r\r\n`
/// loses its `\r\n` and keeps one `\r` (ensembl-io TextParser.pm:138 replaces
/// one `\r\n`, then `chomp`), which fails the expression; the skip warning
/// quotes the line with every run of white space one space (BaseVEP.pm:849),
/// the parser's refusal quotes it the same way, and no control byte reaches the
/// warnings file. The record after it, ending `\r\n`, is annotated.
#[test]
fn a_kept_carriage_return_is_quoted_as_a_space() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "carriage.region.txt",
        "21:14221107-14221107/A\r\r\n21:14221107-14221107/A\r\n",
    );
    let (rows, warnings) = annotate_with_warnings(dir.path(), &region, &["--format", "region"]);
    assert!(
        only_record_has_rows(&rows, "21:14221107-14221107/A"),
        "{rows:?}"
    );
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        ["WARNING: line 1 skipped (21:14221107-14221107/A ): Unable to parse region '21:14221107-14221107/A ': expected chr:start-end[:strand]/allele"],
        "{warnings}"
    );
    assert!(!warnings.contains('\r'), "{warnings:?}");
}

/// A span over `--max_sv_size` is skipped with the text Ensembl VEP's
/// `validate_vf` skips it with (Parser.pm:491-499), the span being end minus
/// start, before any reference is read for it: under a limit of 1000 the
/// duplication over 1001 bases (span 1000) is annotated and the one over 1002
/// bases (span 1001) and a 10,001-base deletion are skipped.
#[test]
fn oversize_span_is_skipped_with_the_max_sv_size_warning() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "long.region.txt",
        "21:10560000-10570000/-\n21:14221000-14222000/DUP\n21:14221000-14222001/DUP\n",
    );
    let (rows, warnings) = annotate_with_warnings(
        dir.path(),
        &region,
        &["--format", "region", "--max_sv_size", "1000"],
    );
    assert!(
        only_record_has_rows(&rows, "21:14221000-14222000/DUP"),
        "{rows:?}"
    );
    for expected in [
        "WARNING: line 1 skipped (21:10560000-10570000/-): variant size (10000) is bigger than --max_sv_size (1000)",
        "WARNING: line 3 skipped (21:14221000-14222001/DUP): variant size (1001) is bigger than --max_sv_size (1000)",
    ] {
        assert!(warnings.contains(expected), "{expected}\n--\n{warnings}");
    }
}

/// The `--chr` list silences `validate_vf`'s refusals of a record on another
/// chromosome (Parser.pm:470-472 runs before the size check at :491-500), so
/// an oversize region on `Z` is dropped without a word under `--chr 21`; a line
/// outside the region expression and a class word the table lacks are refused
/// before `validate_vf` runs (Parser/Region.pm:171, 205-206), so their
/// warnings are written whatever the list says.
#[test]
fn chr_list_silences_the_size_check_alone() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "chr.region.txt",
        "Z:1/A\nZ:1-100/FOO\nZ:1-20000000/A\n21:14221107-14221107/A\n",
    );
    let (rows, warnings) =
        annotate_with_warnings(dir.path(), &region, &["--format", "region", "--chr", "21"]);
    assert!(
        only_record_has_rows(&rows, "21:14221107-14221107/A"),
        "{rows:?}"
    );
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        [
            "WARNING: line 1 skipped (Z:1/A): Unable to parse region 'Z:1/A': expected chr:start-end[:strand]/allele",
            "WARNING: line 2 skipped (Z:1-100/FOO): FOO is not a supported structural variant type",
        ],
        "{warnings}"
    );
}

/// `--format guess` detects a single `chr:start-end[:strand]/allele` field as a
/// region (Utils/VEP.pm:399), so the file reads with no `--format`.
#[test]
fn format_guess_detects_region() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(dir.path(), "records.region.txt", REGION_RECORDS);
    assert_eq!(
        annotate(dir.path(), &region, &[]),
        annotate(dir.path(), &region, &["--format", "region"])
    );
}

/// The record is named by its line as written (Parser/Region.pm:198): a
/// lower-case allele and a `chr` prefix stay in `Uploaded_variation` and in the
/// JSON `id` and `input`, while the chromosome prints as the input spells it
/// (Ensembl VEP 116.2 prints `chr21:14144642-14144643/CT chr21:14144642-14144643
/// CT`).
#[test]
fn record_is_named_by_its_line_as_written() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(
        dir.path(),
        "named.region.txt",
        "chr21:14144642-14144643/CT\n21:15826987-15826987:+1/g\n",
    );
    let rows = annotate(dir.path(), &region, &["--format", "region"]);
    assert!(
        rows.iter().any(|row| row.starts_with("chr21:14144642-14144643/CT\tchr21:14144642-14144643\tCT\t")),
        "{rows:?}"
    );
    assert!(
        rows.iter()
            .any(|row| row.starts_with("21:15826987-15826987:+1/g\t21:15826987\tG\t")),
        "{rows:?}"
    );
    let json = annotate(dir.path(), &region, &["--format", "region", "--json"]);
    assert!(
        json[0].contains(r#""id":"chr21:14144642-14144643/CT""#)
            && json[0].contains(r#""input":"chr21:14144642-14144643/CT""#)
            && json[0].contains(r#""seq_region_name":"chr21""#),
        "{}",
        json[0]
    );
}
