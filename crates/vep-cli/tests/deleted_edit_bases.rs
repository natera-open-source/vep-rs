// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Alleles inside the bases a RefSeq model's alignment deleted, through the
//! `vep` binary on the RefSeq golden corpora. `NM_198689.3` (GRCh37) lacks
//! fifteen genomic bases inside its coding sequence, `21:46020656-46020670`,
//! and two inside its 3' UTR, `21:46021957-46021958`; the COL18A1 models
//! `NM_030582.4` and `NM_130444.3` (GRCh38) lack nine inside their shared coding
//! exon, `21:45504526-45504534`. A substitution, an insertion and a deletion
//! wholly inside such bases take the term of the transcript region whose span
//! holds them, with no position, where Ensembl VEP 116.3 writes
//! `intergenic_variant`, the default it assigns when no predicate fires on a
//! mapper gap; an allele reaching a mapped base reads as both engines read it;
//! a deletion across such bases with a mapped base on each side reads the
//! model's own codons, their residues and the term those codons give; and the
//! warnings file names each such model once, in its record's place among the
//! run's warnings and for no record beyond the batch a run ends on.
//! An allele over the first or last base of the COL18A1 models has an HGVSc
//! of the bases inside the model. Its row at the last base, which follows the
//! nine deleted bases, carries the model's `REFSEQ_OFFSET`; its row at the
//! first base, which precedes them, carries none.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

mod common;

/// The header of every input the tests write.
const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A RefSeq golden corpus of release 116 under `tests/golden/116/`.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The binary on `input` against the corpus's cache and reference under
/// `--refseq`, writing to `out`.
fn refseq_run(dir: &Path, corpus_name: &str, assembly: &str, input: &Path, out: &Path) -> Command {
    let corpus = corpus(corpus_name);
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .args(["--offline", "--cache", "--dir_cache"])
        .arg(corpus.join("cache_root"))
        .args(["--cache_version", "116", "--species", "homo_sapiens"])
        .args(["--assembly", assembly, "--refseq", "--fasta"])
        .arg(corpus.join("reference.fa.gz"))
        .args(["--force_overwrite", "--no_stats", "--quiet"]);
    command
}

/// Runs the binary on `records` against the corpus's cache and reference under
/// `--refseq`, and returns the default-format rows of the named models as
/// sorted `id model consequence cDNA CDS protein IMPACT` lines, and the
/// warnings file's lines.
fn annotate(
    dir: &Path,
    corpus_name: &str,
    assembly: &str,
    records: &str,
    models: &[&str],
) -> (Vec<String>, Vec<String>) {
    let input = dir.join("in.vcf");
    fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join("out.txt");
    let output = refseq_run(dir, corpus_name, assembly, &input, &out)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cells| models.contains(&cells[4]))
        .map(|cells| {
            let impact = cells[13]
                .split(';')
                .find_map(|item| item.strip_prefix("IMPACT="))
                .unwrap_or("");
            format!(
                "{} {} {} {} {} {} {impact}",
                cells[0], cells[4], cells[6], cells[7], cells[8], cells[9]
            )
        })
        .collect();
    rows.sort();
    let warnings = fs::read_to_string(format!("{}_warnings.txt", out.display()))
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    (rows, warnings)
}

/// The warning one model with deleted bases draws.
fn deleted_bases_warning(model: &str) -> String {
    format!(
        "WARNING: Bases deleted by the alignment of {model} carry no transcript coordinate; \
         a variant inside them takes the term of the transcript region whose span holds it"
    )
}

/// On `NM_198689.3` every allele kind wholly inside the fifteen deleted coding
/// bases is `coding_sequence_variant` and every one inside the two deleted 3'
/// UTR bases is `3_prime_UTR_variant`, none with a position. A deletion from
/// the deleted bases onto the next mapped base keeps the open-ended position
/// both engines write (`?-160`, `?-1446`), and an insertion between the last
/// mapped coding base and the first deleted one is read at the mapped codon.
/// The warnings file names the model once.
#[test]
fn alleles_inside_the_deleted_bases_of_a_single_exon_model() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t46020662\tcds_substitution\tC\tT\t.\t.\t.
21\t46020662\tcds_insertion\tC\tCA\t.\t.\t.
21\t46020660\tcds_deletion\tTGCG\tT\t.\t.\t.
21\t46020655\tcds_whole_span\tCCTGCTGCGCCCCCAG\tC\t.\t.\t.
21\t46020656\tcds_onto_mapped\tCCTGCTGCGCCCCCAG\tC\t.\t.\t.
21\t46020655\tcds_after_mapped\tC\tCA\t.\t.\t.
21\t46021957\tutr_substitution\tT\tC\t.\t.\t.
21\t46021958\tutr_substitution_2\tG\tA\t.\t.\t.
21\t46021957\tutr_insertion\tT\tTA\t.\t.\t.
21\t46021956\tutr_deletion\tGTG\tG\t.\t.\t.
21\t46021956\tutr_onto_mapped\tGTGG\tG\t.\t.\t.
";
    let (rows, warnings) = annotate(
        dir.path(),
        "GRCh37-refseq",
        "GRCh37",
        records,
        &["NM_198689.3"],
    );
    assert_eq!(
        rows,
        [
            "cds_after_mapped NM_198689.3 frameshift_variant 159-160 134-135 45 HIGH",
            "cds_deletion NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "cds_insertion NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "cds_onto_mapped NM_198689.3 coding_sequence_variant ?-160 ?-135 ?-45 MODIFIER",
            "cds_substitution NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "cds_whole_span NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "utr_deletion NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
            "utr_insertion NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
            "utr_onto_mapped NM_198689.3 3_prime_UTR_variant ?-1446 - - MODIFIER",
            "utr_substitution NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
            "utr_substitution_2 NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
        ]
    );
    assert_eq!(warnings, [deleted_bases_warning("NM_198689.3")]);
}

/// On the two COL18A1 models every allele kind wholly inside the nine deleted
/// coding bases is `coding_sequence_variant` with no position, while the
/// unedited sibling `NM_001379500.1`, which has the bases, reads each at its
/// codon; a deletion from two mapped bases into the deleted ones keeps the
/// open-ended position both engines write. The warnings file names each edited
/// model once.
#[test]
fn alleles_inside_the_deleted_bases_of_a_multi_exon_model() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t45504526\tsubstitution\tC\tA\t.\t.\t.
21\t45504528\tinsertion\tC\tCT\t.\t.\t.
21\t45504527\tdeletion\tCCCC\tC\t.\t.\t.
21\t45504525\twhole_span\tCCCCCCCCCA\tC\t.\t.\t.
21\t45504523\tfrom_mapped\tCCCCCCA\tC\t.\t.\t.
";
    let models = ["NM_030582.4", "NM_130444.3", "NM_001379500.1"];
    let (rows, warnings) = annotate(dir.path(), "GRCh38-refseq", "GRCh38", records, &models);
    assert_eq!(
        rows,
        [
            "deletion NM_001379500.1 inframe_deletion 2906-2908 2840-2842 947-948 MODERATE",
            "deletion NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "deletion NM_130444.3 coding_sequence_variant - - - MODIFIER",
            "from_mapped NM_001379500.1 inframe_deletion 2902-2907 2836-2841 946-947 MODERATE",
            "from_mapped NM_030582.4 coding_sequence_variant 3385-? 3376-? 1126-? MODIFIER",
            "from_mapped NM_130444.3 coding_sequence_variant 4090-? 4081-? 1361-? MODIFIER",
            "insertion NM_001379500.1 frameshift_variant 2906-2907 2840-2841 947 HIGH",
            "insertion NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "insertion NM_130444.3 coding_sequence_variant - - - MODIFIER",
            "substitution NM_001379500.1 synonymous_variant 2904 2838 946 LOW",
            "substitution NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "substitution NM_130444.3 coding_sequence_variant - - - MODIFIER",
            "whole_span NM_001379500.1 inframe_deletion 2904-2912 2838-2846 946-949 MODERATE",
            "whole_span NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "whole_span NM_130444.3 coding_sequence_variant - - - MODIFIER",
        ]
    );
    assert_eq!(
        warnings.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            deleted_bases_warning("NM_030582.4"),
            deleted_bases_warning("NM_130444.3"),
        ])
    );
}

/// The `id model Consequence cDNA_position HGVSc REFSEQ_OFFSET` lines of a
/// `--refseq --hgvs` run on `records` over the GRCh38 corpus, for the named
/// models, sorted, `-` for a key the row does not carry.
fn hgvs_rows(dir: &Path, records: &str, models: &[&str]) -> Vec<String> {
    let input = dir.join("in.vcf");
    fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join("out.txt");
    let output = refseq_run(dir, "GRCh38-refseq", "GRCh38", &input, &out)
        .arg("--hgvs")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cells| models.contains(&cells[4]))
        .map(|cells| {
            let extra = |key: &str| {
                cells[13]
                    .split(';')
                    .find_map(|item| item.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
                    .unwrap_or("-")
            };
            format!(
                "{} {} {} {} {} {}",
                cells[0],
                cells[4],
                cells[6],
                cells[7],
                extra("HGVSc"),
                extra("REFSEQ_OFFSET")
            )
        })
        .collect();
    rows.sort();
    rows
}

/// The COL18A1 models `NM_030582.4` and `NM_130444.3` and their unedited
/// sibling `NM_001379500.1` end on `21:45513720`, c.*1322 of all three (cDNA
/// 5882, 6587 and 5408, their coding sequences ending at 4560, 5265 and 4086),
/// and the two edited models start on `21:45455522`, cDNA 1 and c.-9, their
/// coding sequence starting at cDNA 10. A substitution of the last base and
/// the base after it, `AC` to `CG`, changes c.*1322 from `A` to `C` on each
/// model, and one of the base before the first base and the first base, `CA`
/// to `GC`, changes c.-9 from `A` to `C`: an edited model's notation names
/// the base inside alone, as the unedited model's does. The row of each
/// edited model at its last base carries `REFSEQ_OFFSET=-9`, the nine bases
/// its alignment deleted before that base (`_rna_edit` `3387 3395` and
/// `4092 4100`). Ensembl VEP 116.3 writes `c.*1322delinsCG` and
/// `c.-9delinsGC`, every alternate base over the base inside
/// (`hgvs_clipped_reference_against_whole_alternate`); the sibling's intronic
/// row of the second record is the one it writes.
#[test]
fn a_substitution_over_an_end_of_an_edited_model_names_the_base_inside() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t45455521\tfirst_both\tCA\tGC\t.\t.\t.
21\t45513720\tover_both\tAC\tCG\t.\t.\t.
";
    let models = ["NM_030582.4", "NM_130444.3", "NM_001379500.1"];
    assert_eq!(
        hgvs_rows(dir.path(), records, &models),
        [
            "first_both NM_001379500.1 intron_variant - NM_001379500.1:c.107-12721_107-12720delinsGC -",
            "first_both NM_030582.4 5_prime_UTR_variant ?-1 NM_030582.4:c.-9A>C -",
            "first_both NM_130444.3 5_prime_UTR_variant ?-1 NM_130444.3:c.-9A>C -",
            "over_both NM_001379500.1 3_prime_UTR_variant 5408-? NM_001379500.1:c.*1322A>C -",
            "over_both NM_030582.4 3_prime_UTR_variant 5882-? NM_030582.4:c.*1322A>C -9",
            "over_both NM_130444.3 3_prime_UTR_variant 6587-? NM_130444.3:c.*1322A>C -9",
        ]
    );
}

/// `21 45513720 AC AG` keeps the models' last base and changes the base
/// after it, so no base of a model changes and no row carries an HGVSc. The
/// row of each edited model lies inside the model, at cDNA 5882 and 6587,
/// after the nine bases its alignment deleted (`_rna_edit` `3387 3395` and
/// `4092 4100`), so it carries `REFSEQ_OFFSET=-9`, as every row of the model
/// after those bases does under `--hgvs`; the unedited sibling carries none.
#[test]
fn a_change_past_the_last_base_of_an_edited_model_keeps_the_models_offset() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t45513720\tpast_only\tAC\tAG\t.\t.\t.\n";
    let models = ["NM_030582.4", "NM_130444.3", "NM_001379500.1"];
    assert_eq!(
        hgvs_rows(dir.path(), records, &models),
        [
            "past_only NM_001379500.1 3_prime_UTR_variant 5408-? - -",
            "past_only NM_030582.4 3_prime_UTR_variant 5882-? - -9",
            "past_only NM_130444.3 3_prime_UTR_variant 6587-? - -9",
        ]
    );
}

/// The line of the corpus's `variants.vcf` whose ID column is `id`.
fn corpus_record(corpus_name: &str, id: &str) -> String {
    fs::read_to_string(corpus(corpus_name).join("variants.vcf"))
        .unwrap()
        .lines()
        .find(|line| line.split('\t').nth(2) == Some(id))
        .unwrap()
        .to_string()
}

/// The warning the line `quoted`, record `line` of the input, draws for
/// naming chromosome `Z`, which the cache lacks.
fn unknown_chromosome_warning(line: u64, quoted: &str) -> String {
    format!(
        "WARNING: line {line} skipped ({quoted}): Chromosome Z not found in annotation sources \
         or synonyms; chromosome Z does not overlap any features"
    )
}

/// The warning a model's deleted bases draw stands in its record's place among
/// the run's warnings: after the skips of its batch and after the warning an
/// earlier record draws, in both pipeline modes, at any worker count, with and
/// without a plugin in the run. Record 1 is the corpus record whose model and
/// assembly disagree, record 2 names a chromosome the cache lacks and record 3
/// lies inside the deleted coding bases of `NM_198689.3`.
#[test]
fn the_deleted_bases_warning_stands_in_its_records_place() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(
        &input,
        format!(
            "{VCF_HEADER}{}\n\
             Z\t1000\tunknown_chromosome\tA\tG\t.\t.\t.\n\
             21\t46020662\tinside_deleted\tC\tT\t.\t.\t.\n",
            corpus_record("GRCh37-refseq", "rs1191732671")
        ),
    )
    .unwrap();
    let out = dir.path().join("out.txt");
    let warnings = dir.path().join("out.txt_warnings.txt");
    let expected = format!(
        "{}\nWARNING: Transcript-assembly mismatch in rs1191732671\n{}\n",
        unknown_chromosome_warning(2, "Z 1000 unknown_chromosome A G . . ."),
        deleted_bases_warning("NM_198689.3")
    );
    for plugin in [&[][..], &["--plugin", "TSSDistance"]] {
        for mode in ["streaming", "batched"] {
            for fork in ["1", "2", "8"] {
                let _ = fs::remove_file(&warnings);
                let output = refseq_run(dir.path(), "GRCh37-refseq", "GRCh37", &input, &out)
                    .args(["--pipeline_mode", mode, "--fork", fork])
                    .args(plugin)
                    .output()
                    .unwrap();
                let run = format!("{mode}, --fork {fork}, {plugin:?}");
                assert!(
                    output.status.success(),
                    "{run}: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
                assert_eq!(fs::read_to_string(&warnings).unwrap(), expected, "{run}");
            }
        }
    }
}

/// A run that ends early leaves no deleted-bases warning of a record beyond
/// the batch it ended on, in either pipeline mode at any worker count. Of 400
/// records written as Parquet in batches of 100, record 30 names a chromosome
/// the cache lacks, record 150, in the second batch, is out of order and
/// record 230, in the third, lies inside the deleted coding bases of
/// `NM_198689.3`: every run is refused at record 150 and its warnings file
/// holds the skip of record 30 alone. A streaming run at the built-in
/// `--pipeline_depth` starts no batch beyond the third when its writer ends
/// in the second, so the third is the later batch both modes annotate.
#[cfg(feature = "parquet")]
#[test]
fn a_run_that_ends_early_leaves_no_deleted_bases_warning_of_a_later_batch() {
    let dir = tempfile::tempdir().unwrap();
    let mut text = String::from(VCF_HEADER);
    for line in 1..=400u64 {
        text.push_str(&match line {
            30 => format!("Z\t1030\tr{line}\tA\tG\t.\t.\t.\n"),
            150 => format!("21\t5\tr{line}\tA\tG\t.\t.\t.\n"),
            230 => format!("21\t46020662\tr{line}\tC\tT\t.\t.\t.\n"),
            _ if line < 230 => format!("21\t{}\tr{line}\tA\tG\t.\t.\t.\n", 10_000_000 + line * 100),
            _ => format!("21\t{}\tr{line}\tA\tG\t.\t.\t.\n", 46_100_000 + line),
        });
    }
    let input = dir.path().join("in.vcf");
    fs::write(&input, text).unwrap();
    let out = dir.path().join("out.parquet");
    let warnings = dir.path().join("out.parquet_warnings.txt");
    let expected = format!(
        "{}\n",
        unknown_chromosome_warning(30, "Z 1030 r30 A G . . .")
    );
    for mode in ["streaming", "batched"] {
        for fork in ["1", "2", "8", "32"] {
            for repeat in 0..3 {
                let _ = fs::remove_file(&warnings);
                let output = refseq_run(dir.path(), "GRCh37-refseq", "GRCh37", &input, &out)
                    .args(["--output_format", "parquet", "--buffer_size", "100"])
                    .args(["--pipeline_mode", mode, "--fork", fork])
                    .output()
                    .unwrap();
                let run = format!("{mode}, --fork {fork}, run {repeat}");
                assert_eq!(output.status.code(), Some(1), "{run}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr).trim_end(),
                    "ERROR: [parquet] input is not sorted by chromosome and position at line 150 \
                     (21:5 after 21:10014900); sort the input or pass --no_check_variants_order",
                    "{run}"
                );
                assert!(!out.exists(), "{run}: the refused run leaves no output");
                assert_eq!(fs::read_to_string(&warnings).unwrap(), expected, "{run}");
            }
        }
    }
}

/// Runs the binary as `annotate` does and returns the whole default-format
/// rows of the named models, sorted.
fn whole_rows(
    dir: &Path,
    corpus_name: &str,
    assembly: &str,
    records: &str,
    models: &[&str],
) -> Vec<String> {
    let input = dir.join("in.vcf");
    fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join("out.txt");
    let output = refseq_run(dir, corpus_name, assembly, &input, &out)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| {
            line.split('\t')
                .nth(4)
                .is_some_and(|model| models.contains(&model))
        })
        .map(str::to_string)
        .collect();
    rows.sort();
    rows
}

/// A deletion that removes the nine deleted bases of the two COL18A1 models
/// with a mapped base on each side has more genomic bases than coding
/// positions. Its reference codons are the models' own bases at the row's
/// protein positions and its reference amino acids their translation, where
/// Ensembl VEP reads both from the genomic allele written over the shorter
/// span (TranscriptVariationAllele.pm:872-877, :703). Seventeen bases hold
/// eight coding positions, a frameshift; eighteen hold nine, three whole
/// codons; twelve hold three inside two codons, an `inframe_deletion`, where
/// the allele's bases fail Ensembl VEP's in-frame test
/// (Utils/VariationEffect.pm:1178-1185) and leave `protein_altering_variant`.
#[test]
fn a_deletion_across_the_deleted_bases_reads_the_models_codons_and_residues() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t45504518\trs1179352929\tCCAGGCCCCCCAGGCCCA\tC\t.\t.\t.
21\t45504520\t1338712\tAGGCCCCCCAGGCCCACGT\tA\t.\t.\t.
21\t45504524\ttwelve_bases\tCCCCCAGGCCCAC\tC\t.\t.\t.
";
    let rows = whole_rows(
        dir.path(),
        "GRCh38-refseq",
        "GRCh38",
        records,
        &["NM_030582.4", "NM_130444.3"],
    );
    assert_eq!(
        rows,
        [
            "1338712\t21:45504521-45504538\t-\t80781\tNM_030582.4\tTranscript\tinframe_deletion\t3382-3390\t3373-3381\t1125-1127\tGPR/-\tGGCCCACGT/-\t-\tIMPACT=MODERATE;STRAND=1;GIVEN_REF=GGCCCCCCAGGCCCACGT;USED_REF=GGCCCCCCAGGCCCACGT;BAM_EDIT=OK",
            "1338712\t21:45504521-45504538\t-\t80781\tNM_130444.3\tTranscript\tinframe_deletion\t4087-4095\t4078-4086\t1360-1362\tGPR/-\tGGCCCACGT/-\t-\tIMPACT=MODERATE;STRAND=1;GIVEN_REF=GGCCCCCCAGGCCCACGT;USED_REF=GGCCCCCCAGGCCCACGT;BAM_EDIT=OK",
            "rs1179352929\t21:45504519-45504535\t-\t80781\tNM_030582.4\tTranscript\tframeshift_variant\t3380-3387\t3371-3378\t1124-1126\tPGP/X\tcCAGGCCCA/c\t-\tIMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK",
            "rs1179352929\t21:45504519-45504535\t-\t80781\tNM_130444.3\tTranscript\tframeshift_variant\t4085-4092\t4076-4083\t1359-1361\tPGP/X\tcCAGGCCCA/c\t-\tIMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK",
            "twelve_bases\t21:45504525-45504536\t-\t80781\tNM_030582.4\tTranscript\tinframe_deletion\t3386-3388\t3377-3379\t1126-1127\tPR/R\tcCACGT/cgt\t-\tIMPACT=MODERATE;STRAND=1;GIVEN_REF=CCCCAGGCCCAC;USED_REF=CCCCAGGCCCAC;BAM_EDIT=OK",
            "twelve_bases\t21:45504525-45504536\t-\t80781\tNM_130444.3\tTranscript\tinframe_deletion\t4091-4093\t4082-4084\t1361-1362\tPR/R\tcCACGT/cgt\t-\tIMPACT=MODERATE;STRAND=1;GIVEN_REF=CCCCAGGCCCAC;USED_REF=CCCCAGGCCCAC;BAM_EDIT=OK",
        ]
    );
}
