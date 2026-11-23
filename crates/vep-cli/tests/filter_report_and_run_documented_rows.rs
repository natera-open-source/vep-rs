// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `vep` binary writes on the inputs of ten documented differences
//! whose code stands in this crate: four of the transcript filter, three of the
//! statistics report, one of `vep recode`, one of `--buffer_size` and one of the
//! VCF reader. Each test runs the input of its entry's record under
//! `scripts/divergences/reproductions` over the committed corpus that record
//! names and holds the rows, the report section or the refusal the entry is
//! about, whole.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A committed corpus under `tests/golden/116`, by its directory name.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// A VCF of one record line under the header.
fn vcf(dir: &Path, name: &str, record: &str) -> PathBuf {
    write(dir, name, &format!("{VCF_HEADER}{record}\n"))
}

/// One `vep` run on `input` over the JSON cache of the corpus `cache`, writing
/// `out`, with `flags` after the arguments every run of this file takes. The
/// run's working directory is `dir`, so a flag names a file of `dir` by its
/// name alone.
fn vep(dir: &Path, input: &Path, out: &Path, cache: &str, flags: &[&str]) -> Output {
    common::isolated_vep(dir)
        .current_dir(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--json_cache")
        .arg(corpus(cache).join("json_cache"))
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--quiet",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The data rows of a run that ends with exit code 0 and writes no statistics
/// report: the lines of `<dir>/<name>` that are not `#` lines.
fn rows(dir: &Path, name: &str, input: &Path, cache: &str, flags: &[&str]) -> Vec<String> {
    let out = dir.join(name);
    let mut all = vec!["--no_stats"];
    all.extend_from_slice(flags);
    let output = vep(dir, input, &out, cache, &all);
    assert!(
        output.status.success(),
        "{name} {flags:?}: {}",
        stderr(&output)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The lines of the section `title` in the text statistics report of a run that
/// ends with exit code 0: the lines after `[<title>]` up to the next empty line.
fn report_section(
    dir: &Path,
    name: &str,
    input: &Path,
    cache: &str,
    flags: &[&str],
    title: &str,
) -> Vec<String> {
    let out = dir.join(name);
    let mut all = vec!["--stats_text"];
    all.extend_from_slice(flags);
    let output = vep(dir, input, &out, cache, &all);
    assert!(
        output.status.success(),
        "{name} {flags:?}: {}",
        stderr(&output)
    );
    let report = std::fs::read_to_string(format!("{}_summary.txt", out.display())).unwrap();
    let heading = format!("[{title}]");
    report
        .lines()
        .skip_while(|line| *line != heading)
        .skip(1)
        .take_while(|line| !line.is_empty())
        .map(str::to_string)
        .collect()
}

/// The record of the transcript-filter corpus that 41 transcripts reach: 34 of
/// the gene `HGNC:1132`, 5 of `HGNC:53145` and 2 of `HGNC:2559`.
const RECORD_OF_THREE_GENES: &str = "21\t17608184\ttf_17608184\tA\tG\t.\t.\t.";

/// The record of the transcript-filter corpus that three transcripts reach,
/// `ENST00000624153` among them.
const RECORD_OF_THREE_TRANSCRIPTS: &str = "21\t7130369\ttf_7130369\tA\tG\t.\t.\t.";

/// The row of [`RECORD_OF_THREE_TRANSCRIPTS`] on `ENST00000624153`.
const ROW_OF_THE_LISTED_TRANSCRIPT: &str = "tf_7130369\t21:7130369\tG\tENSG00000280330\t\
    ENST00000624153\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\t\
    IMPACT=MODIFIER;DISTANCE=82;STRAND=1";

/// An `in` list that mixes a number with other text is judged by its first
/// written entry: `HGNC:2559,1132` compares each transcript's `HGNC:<number>`
/// whole and keeps the two transcripts of `HGNC:2559`, and `1132,HGNC:2559`
/// reads each value through its number and keeps the 34 of `HGNC:1132`.
/// Ensembl VEP 116.3 judges the list by one key of its hash, `(keys
/// %{$value})[0]` (`FilterSet.pm:597-599`), which the hash seed of the process
/// picks: a numeric key reads the number (lines 600-601), any other compares
/// the value whole (lines 603-604), so either written order gives either set.
#[test]
fn an_in_list_of_a_number_and_other_text_is_judged_by_its_first_written_entry() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(dir.path(), "three-genes.vcf", RECORD_OF_THREE_GENES);
    let filtered = |name: &str, list: &str| {
        rows(
            dir.path(),
            name,
            &input,
            "GRCh38-transcript-filter",
            &[
                "--symbol",
                "--transcript_filter",
                &format!("gene_hgnc_id in {list}"),
            ],
        )
    };
    assert_eq!(
        filtered("text-first.txt", "HGNC:2559,1132"),
        [
            "tf_17608184\t21:17608184\tG\tENSG00000154639\tENST00001092433\tTranscript\t\
             downstream_gene_variant\t-\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;DISTANCE=426;STRAND=1;SYMBOL=CXADR;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:2559",
            "tf_17608184\t21:17608184\tG\tENSG00000154639\tENST00001128578\tTranscript\t\
             intron_variant\t-\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;STRAND=1;SYMBOL=CXADR;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:2559",
        ]
    );
    let number_first = filtered("number-first.txt", "1132,HGNC:2559");
    assert_eq!(number_first.len(), 34, "{number_first:?}");
    for row in &number_first {
        assert!(
            row.ends_with(";SYMBOL=BTG3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:1132"),
            "{row}"
        );
    }
}

/// The file an `in` list names is opened under its written name, whatever
/// character the name ends or begins with. A file named `true|`, written by its
/// path, and one named `<ids.txt`, written by its name alone beside an
/// `ids.txt` that lists another transcript, each list `ENST00000624153`, and
/// each run keeps that transcript's row alone. Ensembl VEP 116.3 finds the file
/// with `-e` (`FilterSet.pm:396`) and opens it with the two-argument `open`
/// (line 397), which reads a `|` that ends the written string as a pipe from
/// the command the rest of it spells and a `<` that begins it as the mode of
/// the name after it.
#[test]
fn an_in_list_file_is_opened_under_its_written_name() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "three-transcripts.vcf",
        RECORD_OF_THREE_TRANSCRIPTS,
    );
    let piped = write(dir.path(), "true|", "ENST00000624153\n");
    write(dir.path(), "<ids.txt", "ENST00000624153\n");
    write(dir.path(), "ids.txt", "ENST00001020056\n");
    for (out, list) in [
        ("pipe.txt", piped.display().to_string()),
        ("mode.txt", "<ids.txt".to_string()),
    ] {
        let kept = rows(
            dir.path(),
            out,
            &input,
            "GRCh38-transcript-filter",
            &["--transcript_filter", &format!("stable_id in {list}")],
        );
        assert_eq!(kept, [ROW_OF_THE_LISTED_TRANSCRIPT], "{list}");
    }
}

/// A filter unit that opens with the word `0` has `0` as its field, a name no
/// transcript key matches: `0 is 1` fails on every transcript, so the record is
/// one `intergenic_variant` row and the run ends with exit code 0 and without
/// `No such operator` on standard error, and `0 ne` and `0 is`, whose field is
/// undefined, keep all 41 rows. Ensembl VEP 116.3 tests the unit's field slot
/// with Perl's truth test (`FilterSet.pm:260`), to which `0` is unset, so the
/// next word takes the slot and `0 is 1` is refused with `No such operator "1"`
/// (line 276).
#[test]
fn a_filter_unit_that_opens_with_the_word_zero_has_zero_as_its_field() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(dir.path(), "three-genes.vcf", RECORD_OF_THREE_GENES);
    let filtered = |name: &str, flags: &[&str]| {
        rows(dir.path(), name, &input, "GRCh38-transcript-filter", flags)
    };
    let output = vep(
        dir.path(),
        &input,
        &dir.path().join("zero-is-one.streams.txt"),
        "GRCh38-transcript-filter",
        &["--no_stats", "--transcript_filter", "0 is 1"],
    );
    assert!(
        !stderr(&output).contains("No such operator"),
        "{}",
        stderr(&output)
    );
    assert_eq!(
        filtered("zero-is-one.txt", &["--transcript_filter", "0 is 1"]),
        ["tf_17608184\t21:17608184\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER"]
    );
    let unfiltered = filtered("unfiltered.txt", &[]);
    assert_eq!(unfiltered.len(), 41, "{unfiltered:?}");
    assert_eq!(
        filtered("zero-ne.txt", &["--transcript_filter", "0 ne"]),
        unfiltered
    );
    assert_eq!(
        filtered("zero-is.txt", &["--transcript_filter", "0 is"]),
        unfiltered
    );
}

/// A directory named as an `in` list is refused: exit code 1, the line `ERROR:
/// [filter] Could not read from file <path>` on standard error and no output
/// file. Ensembl VEP 116.3 finds the directory with `-e` and opens it
/// (`FilterSet.pm:396-397`), its first read fails, which the loop of line 398
/// takes for the end of the file, and the empty list drops every transcript, so
/// the record is written as `intergenic_variant` with nothing on standard
/// error.
#[test]
fn a_directory_named_as_an_in_list_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "three-transcripts.vcf",
        RECORD_OF_THREE_TRANSCRIPTS,
    );
    let list = dir.path().join("list_dir");
    std::fs::create_dir(&list).unwrap();
    let out = dir.path().join("directory.txt");
    let output = vep(
        dir.path(),
        &input,
        &out,
        "GRCh38-transcript-filter",
        &[
            "--no_stats",
            "--transcript_filter",
            &format!("stable_id in {}", list.display()),
        ],
    );
    assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
    let refusal = format!(
        "ERROR: [filter] Could not read from file {}",
        list.display()
    );
    assert!(
        stderr(&output).lines().any(|line| line == refusal),
        "{}",
        stderr(&output)
    );
    assert!(!out.exists(), "the refused run writes no output");
}

/// A coding row whose position lies past the end of its peptide is counted
/// under `90-100%` in the report's `Position in protein` section: the record's
/// row on `ENST00000451131`, a model of ten residues whose coding sequence has
/// no stop, lies at residue 11. Ensembl VEP 116.3 keys the count by `int(10 *
/// ($protein_pos / $protein_length))` (`Stats.pm:359-360`), 11 here, and folds
/// the key 10 alone into the last tenth (lines 461-463), so its section holds
/// an eleventh row, `110-120%`.
#[test]
fn a_row_past_the_end_of_its_peptide_is_counted_in_the_last_tenth_of_the_protein() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "short-peptide.vcf",
        "21\t39428463\ttiny_cds_end\tT\tC\t.\t.\t.",
    );
    assert_eq!(
        report_section(
            dir.path(),
            "short-peptide.txt",
            &input,
            "GRCh38-individual",
            &[],
            "Position in protein",
        ),
        [
            "00-10%\t9",
            "10-20%\t6",
            "20-30%\t0",
            "30-40%\t1",
            "40-50%\t0",
            "50-60%\t0",
            "60-70%\t0",
            "70-80%\t0",
            "80-90%\t0",
            "90-100%\t1",
        ]
    );
}

/// The report's `Variant classes` section counts a record that `--minimal`
/// splits under the class of its first alternate allele's trimmed pair: the
/// record of a one-base substitution and an insertion, which no transcript
/// reaches, is one `SNV`. Ensembl VEP 116.3 logs the class of the feature as it
/// stands when its logger meets it (`Stats.pm:246-247`), the allele string of
/// the whole record on a run without `--variant_class` and without a
/// transcript at the record, so its section reads `indel 1`.
#[test]
fn a_record_minimal_splits_is_counted_under_the_class_of_its_first_allele() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "two-alleles.vcf",
        "21\t5030602\tmulti_no_transcript\tA\tG,AT\t.\t.\t.",
    );
    assert_eq!(
        report_section(
            dir.path(),
            "two-alleles.txt",
            &input,
            "GRCh38-individual",
            &["--minimal"],
            "Variant classes",
        ),
        ["SNV\t1"]
    );
}

/// A record its reader marks for skipping and `--dont_skip` keeps, a `<DEL>`
/// with neither `END` nor `SVLEN`, is one of the report's `Variants processed`
/// in the VCF format as in every other. Ensembl VEP 116.3 writes such a record
/// through as its input line without building its output hashes
/// (`OutputFactory/VCF.pm:342-345`), short of the logger that counts a record,
/// so its report's `Variants processed` row carries no count.
#[test]
fn a_record_dont_skip_keeps_is_counted_as_processed_in_the_vcf_format() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "kept.vcf",
        "21\t30165850\tdel_no_end_no_svlen\tN\t<DEL>\t.\t.\tSVTYPE=DEL",
    );
    assert_eq!(
        report_section(
            dir.path(),
            "kept.vcf.out",
            &input,
            "GRCh38-individual-zyg",
            &["--dont_skip", "--vcf"],
            "General statistics",
        ),
        [
            "Lines of input read\t1",
            "Variants processed\t1",
            "Variants filtered out\t0",
            "Novel / existing variants\t-",
            "Overlapped genes\t0",
            "Overlapped transcripts\t0",
            "Overlapped regulatory features\t-",
        ]
    );
}

/// `vep recode` keeps the warnings of every buffer fill that held an input
/// given twice: one identifier written after each of two identifiers no source
/// holds is one object that carries both warnings, in input order. Ensembl
/// VEP's `variant_recoder` 116.3 reuses the object an earlier fill made for the
/// same text (`VariantRecoder.pm:468`) and assigns it the list of the fill at
/// hand (line 485), so its object carries the second warning alone.
#[test]
fn recode_keeps_the_warnings_of_every_fill_of_an_input_given_twice() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "ids.txt",
        "rs_nope1\nrs886039732\nrs_nope2\nrs886039732\n",
    );
    let recoder = corpus("GRCh38-recoder");
    let output = common::isolated_vep(dir.path())
        .arg("recode")
        .arg("-i")
        .arg(&input)
        .args(["--format", "id", "--fields", "hgvsg", "--pretty"])
        .arg("--json_cache")
        .arg(recoder.join("json_cache"))
        .arg("--fasta")
        .arg(recoder.join("reference.fa.gz"))
        .arg("--synonyms")
        .arg(recoder.join("chr_synonyms.txt"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        r#"[
   {
      "-" : {
         "hgvsg" : [
            "NC_000021.9:g.32694272_32694276del"
         ],
         "input" : "rs886039732"
      },
      "warnings" : [
         "No variant found with ID 'rs_nope1'",
         "No variant found with ID 'rs_nope2'"
      ]
   }
]
"#
    );
}

/// `--buffer_size 0` writes every row: the record's five rows, the rows of
/// `--buffer_size 1`. Ensembl VEP 116.3 fills its buffer while it holds fewer
/// records than the size (`InputBuffer.pm:191`, `:212`), so at 0 its first fill
/// is empty, which the run reads as the end of the input, and it writes the
/// header alone.
#[test]
fn a_buffer_size_of_zero_writes_every_row() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "one-record.vcf",
        "21\t45573402\tutr_and_intron\tC\tT\t.\t.\t.",
    );
    let sized = |name: &str, size: &str| {
        rows(
            dir.path(),
            name,
            &input,
            "GRCh38-pick-summary",
            &["--buffer_size", size],
        )
    };
    let at_zero = sized("size-zero.txt", "0");
    assert_eq!(
        at_zero,
        [
            "utr_and_intron\t21:45573402\tT\tENSG00000173638\tENST00000651099\tTranscript\t\
             upstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=37;STRAND=-1",
            "utr_and_intron\t21:45573402\tT\tENSG00000299414\tENST00000763267\tTranscript\t\
             upstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=350;STRAND=1",
            "utr_and_intron\t21:45573402\tT\tENSG00000299414\tENST00000763268\tTranscript\t\
             intron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;STRAND=1",
            "utr_and_intron\t21:45573402\tT\tENSG00000299414\tENST00000763269\tTranscript\t\
             upstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=312;STRAND=1",
            "utr_and_intron\t21:45573402\tT\tENSG00000173638\tENST00000859469\tTranscript\t\
             5_prime_UTR_variant\t163\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1",
        ]
    );
    assert_eq!(at_zero, sized("size-one.txt", "1"));
}

/// A sequence record whose first ALT is empty keeps the allele written after
/// it: `,G` on a one-base reference is two alleles, the empty one numbered 1
/// and `G` numbered 2, each with its own row on `ENST00000416044`, and both
/// rows name the record's working string `C//G` and the class
/// `sequence_alteration`. Ensembl VEP 116.3 minimises that string as the one
/// deletion `C/-` and writes the first row alone, under the class `deletion`.
#[test]
fn a_leading_empty_alt_keeps_the_allele_written_after_it() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf(
        dir.path(),
        "leading-empty.vcf",
        "21\t30289183\tlead_empty\tC\t,G\t.\t.\t.",
    );
    assert_eq!(
        rows(
            dir.path(),
            "leading-empty.txt",
            &input,
            "GRCh38-vcf-edge-cases",
            &["--allele_number", "--uploaded_allele", "--variant_class"],
        ),
        [
            "lead_empty\t21:30289183\t-\tENSG00000232263\tENST00000416044\tTranscript\t\
             stop_retained_variant\t332\t308\t103\t*/X\ttGa/ta\t-\t\
             ALLELE_NUM=1;UPLOADED_ALLELE=C//G;IMPACT=LOW;STRAND=-1;\
             VARIANT_CLASS=sequence_alteration",
            "lead_empty\t21:30289183\tG\tENSG00000232263\tENST00000416044\tTranscript\t\
             stop_lost\t332\t308\t103\t*/S\ttGa/tCa\t-\t\
             ALLELE_NUM=2;UPLOADED_ALLELE=C//G;IMPACT=HIGH;STRAND=-1;\
             VARIANT_CLASS=sequence_alteration",
        ]
    );
}
