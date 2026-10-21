// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl-format records whose alleles are given on the reverse strand annotate
//! as Ensembl VEP 116.2 annotates them: the consequences, codons and HGVS
//! notations read the alleles on the forward strand, the `Allele` column and the
//! JSON keep them as written, and a shifted insertion's HGVSp and shifted codons
//! turn the allele as the record wrote it. Every expected row is Ensembl VEP
//! 116.2's output for the same record on the release-116 GRCh38 cache.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38/json_cache")
}

/// A chromosome 21 FASTA holding the two windows the records below read, `N`
/// elsewhere: 26021751-26022050 around the APP insertion site and
/// 46302401-46302700 around the ENST00000291691 one (GRCh38).
fn write_fasta(dir: &Path) -> PathBuf {
    const WINDOWS: [(usize, &str); 2] = [
        (
            26_021_751,
            "GCAGTGGTGCTAGGGTGGATATTTCCAACTCTGGTATTACGCTTTCTTTTCCTTGGGGGTGGGGGGAATCCAAGCAAATGGTGGATTACCTCGAACCACC\
             TCTTCCACAGACTCTGTGGTGGTGGTGGTGGTGGTGGCAATGCTGGTGGTTCTCTCTGTGGCTTCTTCGTAGGGTTCCTCAGCCTCTTCCTCTACCTCAT\
             CACCATCCTCATCGTCCTCGTCATCATCGGCTTCTTCTTCTTCCACCTCAGCCACTTCTTCCTCCTCTGCTACTTCTACTACTTTGTCTTCACTGTGAAG",
        ),
        (
            46_302_401,
            "CTACACAGATGCAGAAATTTCCAGAAGAAAAACCACCCAGGCCCTGTTTCCTTGGCCTAGGGTTCTGCTGGGGGCTAAGCCTACCTGCAGGGCTGGGGGC\
             AGGGCCTTGAGCATCCAGTTCTGCAGGACCAACTGGTGCACCTGTGCGTTCTGCAGCAGCAGCAGCTCCACCATGTCTGCAGGGAAGAAGCAGGGGGACC\
             CTGAATAAAGTTTCCGTTTTTCCTATTTGTTAAAGTGATAGAGCATTATAGGACCAGAGAACAGGTGTGTCTGTACACTGTGCAGGTCCCCGGGGCAGGC",
        ),
    ];
    let length = WINDOWS
        .iter()
        .map(|(start, window)| start + window.len() - 1)
        .max()
        .unwrap();
    let mut sequence = vec![b'N'; length];
    for (start, window) in WINDOWS {
        sequence[start - 1..start - 1 + window.len()].copy_from_slice(window.as_bytes());
    }
    let fasta = dir.join("21.fa");
    let mut contents = b">21\n".to_vec();
    contents.extend_from_slice(&sequence);
    contents.push(b'\n');
    std::fs::write(&fasta, contents).unwrap();
    std::fs::write(
        dir.join("21.fa.fai"),
        format!("21\t{length}\t4\t{length}\t{}\n", length + 1),
    )
    .unwrap();
    fasta
}

/// The data lines of a run on Ensembl-format `records` with `extra` flags.
fn annotate(dir: &Path, name: &str, records: &str, extra: &[&str]) -> Vec<String> {
    let input = dir.join(format!("{name}.txt"));
    std::fs::write(&input, records).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache())
        .args(["--format", "ensembl", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name} {extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The rows naming one of `transcripts`, in output order.
fn rows_on<'a>(rows: &'a [String], transcripts: &[&str]) -> Vec<&'a str> {
    rows.iter()
        .filter(|row| transcripts.iter().any(|t| row.contains(t)))
        .map(String::as_str)
        .collect()
}

/// Three strand -1 SNVs and the forward-strand spelling of the third: the
/// alleles are read on the forward strand against every transcript, so `A/T -1`
/// at 13639588 (forward `T/A`, whose alternate `A` is the genome's own base) is
/// synonymous on ENST00000299443 and `C/A -1` at 14109002 annotates as its twin
/// `G/T 1` does, the `Allele` column keeping each record's spelling.
const SNV_RECORDS: &str = "\
21\t13639588\t13639588\tA/T\t-1\tcv3308991
21\t13610310\t13610310\tT/C\t-1\tcv3216832
21\t14109002\t14109002\tC/A\t-1\tcv4047524
21\t14109002\t14109002\tG/T\t1\tcv4047524_fwd
";

#[test]
fn reverse_strand_snvs_annotate_on_their_forward_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(dir.path(), "snvs", SNV_RECORDS, &[]);
    assert_eq!(
        rows_on(
            &rows,
            &["ENST00000299443", "ENST00000400211", "ENST00000536861"]
        ),
        [
            "cv3308991\t21:13639588\tT\tENSG00000166351\tENST00000299443\tTranscript\tsynonymous_variant\t1935\t1483\t495\tT\tAct/Act\t-\tIMPACT=LOW;STRAND=1",
            "cv3216832\t21:13610310\tC\tENSG00000166351\tENST00000299443\tTranscript\tmissense_variant\t534\t82\t28\tC/G\tTgc/Ggc\t-\tIMPACT=MODERATE;STRAND=1",
            "cv4047524\t21:14109002\tA\tENSG00000188992\tENST00000400211\tTranscript\tmissense_variant\t1195\t1085\t362\tR/K\taGa/aAa\t-\tIMPACT=MODERATE;STRAND=-1",
            "cv4047524\t21:14109002\tA\tENSG00000188992\tENST00000536861\tTranscript\tsynonymous_variant\t1466\t1356\t452\tK\taaG/aaA\t-\tIMPACT=LOW;STRAND=-1",
            "cv4047524_fwd\t21:14109002\tT\tENSG00000188992\tENST00000400211\tTranscript\tmissense_variant\t1195\t1085\t362\tR/K\taGa/aAa\t-\tIMPACT=MODERATE;STRAND=-1",
            "cv4047524_fwd\t21:14109002\tT\tENSG00000188992\tENST00000536861\tTranscript\tsynonymous_variant\t1466\t1356\t452\tK\taaG/aaA\t-\tIMPACT=LOW;STRAND=-1",
        ]
    );
}

/// The JSON document keeps the record as written (`allele_string` `A/T`,
/// `strand` -1, each consequence's `variant_allele` `T`) while its codons read
/// the forward-strand allele, as Ensembl VEP's does.
#[test]
fn reverse_strand_record_json_keeps_the_written_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(dir.path(), "json", SNV_RECORDS, &["--json"]);
    let document: serde_json::Value = serde_json::from_str(&rows[0]).unwrap();
    let consequence = document["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .find(|tc| tc["transcript_id"] == "ENST00000299443")
        .unwrap();
    assert_eq!(document["id"], "cv3308991");
    assert_eq!(document["allele_string"], "A/T");
    assert_eq!(document["strand"], -1);
    assert_eq!(consequence["variant_allele"], "T");
    assert_eq!(consequence["codons"], "Act/Act");
    assert_eq!(
        consequence["consequence_terms"],
        serde_json::json!(["synonymous_variant"])
    );
}

/// Records on APP (ENST00000346798, reverse strand) and ENST00000291691 (reverse
/// strand) under `--hgvs`. The insertion `-/CTT -1` at 26021981 moves one base
/// 3' like its forward twin `-/AAG 1` and takes the same HGVSc, but its HGVSp
/// inserts Leu where the twin's inserts Ser: `shift_feature_seqs` turns the
/// allele as the record wrote it (TranscriptVariationAllele.pm:1330-1347). The
/// insertion `-/CAG 1` at 46302557 moves four bases along a reverse-strand
/// transcript, more than its three, so the allele is left as written for HGVSp:
/// `p.Gln249delinsProGlu`, as for its reverse-strand twin `-/CTG -1`.
const HGVS_RECORDS: &str = "\
21\t26021896\t26021896\tG/A\t-1\tcv1915683
21\t26021943\t26021943\tT/C\t-1\tcv721913
21\t26021866\t26021871\tGTGGTG/-\t-1\tcv664060
21\t26021981\t26021980\t-/CTT\t-1\tcv2041960
21\t26021981\t26021980\t-/AAG\t1\tcv2041960_fwd
21\t46302557\t46302556\t-/CAG\t1\tins_cag_fwd
21\t46302557\t46302556\t-/CTG\t-1\tins_ctg_rev
";

#[test]
fn reverse_strand_records_take_ensembl_veps_hgvs() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path());
    let rows = annotate(
        dir.path(),
        "hgvs",
        HGVS_RECORDS,
        &["--hgvs", "--fasta", fasta.to_str().unwrap()],
    );
    assert_eq!(
        rows_on(&rows, &["ENST00000346798", "ENST00000291691"]),
        [
            "cv1915683\t21:26021896\tA\tENSG00000142192\tENST00000346798\tTranscript\tmissense_variant\t959\t809\t270\tT/N\taCc/aAc\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000346798.8:c.809C>A;HGVSp=ENSP00000284981.4:p.Thr270Asn",
            "cv721913\t21:26021943\tC\tENSG00000142192\tENST00000346798\tTranscript\tsynonymous_variant\t912\t762\t254\tV\tgtA/gtC\t-\tIMPACT=LOW;STRAND=-1;HGVSc=ENST00000346798.8:c.762A>C;HGVSp=ENSP00000284981.4:p.Val254%3D",
            "cv664060\t21:26021866-26021871\t-\tENSG00000142192\tENST00000346798\tTranscript\tinframe_deletion\t984-989\t834-839\t278-280\tTTT/T\tacCACCACa/aca\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000346798.8:c.834_839del;HGVSp=ENSP00000284981.4:p.Thr279_Thr280del",
            "cv2041960\t21:26021980-26021981\tCTT\tENSG00000142192\tENST00000346798\tTranscript\tinframe_insertion\t874-875\t724-725\t242\tA/AS\tgcc/gCTTcc\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000346798.8:c.725_726insTTC;HGVSp=ENSP00000284981.4:p.Ala242_Asp243insLeu;HGVS_OFFSET=-1",
            "cv2041960_fwd\t21:26021980-26021981\tAAG\tENSG00000142192\tENST00000346798\tTranscript\tinframe_insertion\t874-875\t724-725\t242\tA/AS\tgcc/gCTTcc\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000346798.8:c.725_726insTTC;HGVSp=ENSP00000284981.4:p.Ala242_Asp243insSer;HGVS_OFFSET=-1",
            "ins_cag_fwd\t21:46302556-46302557\tCAG\tENSG00000160298\tENST00000291691\tTranscript\tinframe_insertion\t1873-1874\t741-742\t247-248\t-/L\t-/CTG\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000291691.12:c.743_745dup;HGVSp=ENSP00000291691.8:p.Gln249delinsProGlu;HGVS_OFFSET=-4",
            "ins_ctg_rev\t21:46302556-46302557\tCTG\tENSG00000160298\tENST00000291691\tTranscript\tinframe_insertion\t1873-1874\t741-742\t247-248\t-/L\t-/CTG\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000291691.12:c.743_745dup;HGVSp=ENSP00000291691.8:p.Gln249delinsProGlu;HGVS_OFFSET=-4",
        ]
    );
}

/// Under `--shift_3prime 1` the codons read the allele `shift_feature_seqs`
/// leaves (`codon`, TranscriptVariationAllele.pm:820): the APP insertion given on
/// the reverse strand inserts `TCT` on the transcript (`A/AL`) where its forward
/// twin inserts `TTC` (`A/AS`), and the ENST00000291691 insertion, moved further
/// than its length, keeps `CTG` on the transcript and so alters the protein
/// (`Q/PE`, `protein_altering_variant`) instead of duplicating a leucine. The
/// columns from `Consequence` on are compared, then `Allele`: a row whose codon
/// was read prints the turned allele (OutputFactory.pm:1347), `TCT` and `GAA`
/// for the APP insertions, where the non-coding ENST00000474136 row keeps the
/// written `CTT` and `AAG`, and the ENST00000291691 rows keep theirs, the shift
/// of four exceeding the allele's three bases on the reverse strand.
#[test]
fn shifted_codons_turn_the_allele_as_the_record_wrote_it() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path());
    let rows = annotate(
        dir.path(),
        "shift",
        HGVS_RECORDS,
        &[
            "--hgvs",
            "--shift_3prime",
            "1",
            "--shift_length",
            "--fasta",
            fasta.to_str().unwrap(),
        ],
    );
    let from_consequence = |row: &str| -> String {
        let (head, tail) = row.split_once("\tTranscript\t").unwrap();
        format!("{}\t{tail}", head.split('\t').next().unwrap())
    };
    let shifted: Vec<String> = rows_on(&rows, &["ENST00000346798", "ENST00000291691"])
        .into_iter()
        .filter(|row| row.contains("ins"))
        .map(from_consequence)
        .collect();
    assert_eq!(
        shifted,
        [
            "cv2041960\tinframe_insertion\t875-876\t725-726\t242\tA/AL\tgcc/gcTCTc\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000346798.8:c.725_726insTTC;HGVSp=ENSP00000284981.4:p.Ala242_Asp243insLeu;HGVS_OFFSET=-1;SHIFT_LENGTH=1",
            "cv2041960_fwd\tinframe_insertion\t875-876\t725-726\t242\tA/AS\tgcc/gcTTCc\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000346798.8:c.725_726insTTC;HGVSp=ENSP00000284981.4:p.Ala242_Asp243insSer;HGVS_OFFSET=-1;SHIFT_LENGTH=1",
            "ins_cag_fwd\tprotein_altering_variant\t1877-1878\t745-746\t249\tQ/PE\tcag/cCTGag\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000291691.12:c.743_745dup;HGVSp=ENSP00000291691.8:p.Gln249delinsProGlu;HGVS_OFFSET=-4;SHIFT_LENGTH=4",
            "ins_ctg_rev\tprotein_altering_variant\t1877-1878\t745-746\t249\tQ/PE\tcag/cCTGag\t-\tIMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000291691.12:c.743_745dup;HGVSp=ENSP00000291691.8:p.Gln249delinsProGlu;HGVS_OFFSET=-4;SHIFT_LENGTH=4",
        ]
    );
    let allele = |row: &str| {
        let c: Vec<&str> = row.split('\t').collect();
        format!("{} {} {}", c[0], c[4], c[2])
    };
    let alleles: Vec<String> = rows_on(
        &rows,
        &["ENST00000346798", "ENST00000474136", "ENST00000291691"],
    )
    .into_iter()
    .filter(|row| row.contains("ins"))
    .map(allele)
    .collect();
    assert_eq!(
        alleles,
        [
            "cv2041960 ENST00000346798 TCT",
            "cv2041960 ENST00000474136 CTT",
            "cv2041960_fwd ENST00000346798 GAA",
            "cv2041960_fwd ENST00000474136 AAG",
            "ins_cag_fwd ENST00000291691 CAG",
            "ins_ctg_rev ENST00000291691 CTG",
        ]
    );
}
