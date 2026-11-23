// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSc for a plain insertion whose run of
//! placements, the points at which its bases give the same sequence, is
//! longer than Ensembl VEP's roll reads.
//! Both engines roll an insertion 3' along the chromosome on the transcript's
//! strand while each base it passes equals the next of its bases, cyclically.
//! Ensembl VEP 116.3 reads 1,000 bases of the chromosome either side of the
//! record for that roll (TranscriptVariationAllele.pm:427-437) and stops
//! after at most 1,001 steps less the insertion's length
//! (TranscriptVariationAllele.pm:303-306), so an insertion of more than 500
//! bases that copies the reference bases after its point stops short of the
//! copy's end and reads as an insertion there.
//! The HGVS nomenclature assigns a change to the most 3' position possible,
//! so the notation stands at the insertion's last placement however far the
//! run goes: the duplication of the bases before that point where they equal
//! the inserted bases as they stand there, the insertion of those bases
//! otherwise. A last placement directly past the transcript's last base is
//! the duplication of the transcript's last bases where a whole copy stands
//! before it; without one, or further past the last base, the row has no
//! HGVSc.
//! The HGVSp of such an insertion is read at the same placement as its HGVSc:
//! the row has none where that placement lies in an intron, in a UTR or past
//! the transcript's last base, and the insertion's own where the roll
//! carries it from an intron in between two coding bases. Ensembl VEP 116.3
//! reads the HGVSp at the span its roll reaches
//! (TranscriptVariationAllele.pm:1711-1738).
//! `--shift_3prime 1` annotates the whole row where Ensembl VEP's roll stops
//! and describes it there, and under `--shift_hgvs 0` nothing rolls.
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache; the rows that hold an HGVSp, or
//! the lack of one, are read on made-up coding transcripts of a made-up
//! chromosome.

use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::json;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference the roll walks.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The forward-strand lncRNA of one exon whose bases `n.1` to `n.1800` are
/// 21:32279200-32280999.
const LNCRNA: &str = "ENST00000783173";

/// The forward-strand lncRNA of two exons, `n.1` to `n.413` at
/// 21:14961309-14961721 and `n.414` to `n.562` at 21:14964085-14964233, with
/// an intron of 2,363 bases between them.
const SPLICED_LNCRNA: &str = "ENST00000446301";

/// The reverse-strand lncRNA of eight exons whose last exon, `n.634` to
/// `n.752`, is 21:13943744-13943862, with an intron of 1,176 bases before it.
const SHORT_END_LNCRNA: &str = "ENST00000429521";

/// The reverse-strand lncRNA of four exons whose second exon ends on `n.254`
/// at 21:14857542, with an intron of 31,917 bases after it, and whose last
/// base, `n.5929`, is 21:14818843.
const REVERSE_LNCRNA: &str = "ENST00000432230";

/// The reverse-strand coding transcript whose 3' UTR, `c.*1` to `c.*2502`,
/// is 21:14371115-14373616 on the chromosome, inside its last exon.
const REVERSE_CODING: &str = "ENST00000285667";

/// The forward-strand coding transcript whose last base, `c.*1589`, is
/// 21:15880064.
const FORWARD_CODING: &str = "ENST00000285679";

/// The Consequence of a row inside an exon of a lncRNA.
const EXON: &str = "non_coding_transcript_exon_variant";

/// The Consequence of a row inside an intron of a lncRNA, away from its
/// splice sites.
const INTRON: &str = "intron_variant,non_coding_transcript_variant";

/// The Consequence of a row in a coding transcript's 3' UTR.
const UTR: &str = "3_prime_UTR_variant";

/// `n.301` to `n.820` of [`LNCRNA`], the 520 bases at 21:32279500-32280019.
const LNCRNA_COPY: &str = concat!(
    "TTTTTCTATGTTGCCGGTCTCGAACTCCTAGCCTTAAGCGATCTTCCCTACTCGGCCTCCCAAAATGCCAGGATCG",
    "CAAGCGTAAGCCACTGTGCCAGACCACCACCCTTTTTCTTCTTTTTTGTAAGTGGACTAATTATAGTTGTATATAA",
    "GCAGCTCCACTTTTGATGACCACTGCTCTATTTTGTCGGCTCCTGGGGACCTCCTTTCCATGAAGTGGAGAAGGAA",
    "GGGGATAGAGGAGTCCCTCATCACAGGGTCCTTTTGCGCACTGCTATCCCTGGAATGAGACCCTGGGAAGGGCCCC",
    "GGATCACTCAGGCTGTCACTGCAAGGAAAAAGTTGAGGAAGCAGCAGTGCCCCTGGAAAGCCTGGAAGTGTCCCTG",
    "CTGCTTGGATCAAGGAACAATTCGTCCTTGCAAGGAAGAAGCATCCCTTGTTGATGATAAGGACAGATTTCCTCAG",
    "CTGCCCCAGAAAGCCAGCCTCTGCCTAAAACACCGTAGGTCAGGAGGGAAAGCTGGCGGGGAAC",
);

/// `c.*174` to `c.*733` of [`REVERSE_CODING`] as the transcript reads them,
/// the reverse complement of the 560 bases at 21:14372884-14373443.
const UTR_COPY: &str = concat!(
    "ATTGATCCTGTTTGATTTTGGAGAGATCCTATTCTAACAAATACTCTAAAATGATAAAATTGAGGTACAACTCTCT",
    "TAAAAGAGTATGGATAACTATATTTTCTGGATTCTGGAGGTTGATAACCATATGCACTTAACATTATATTCTATAA",
    "ACATTAAGTAGTGCCAGTTATGAGATTCCCAGTTCTTACTAAATTGTATTAGCAGGAGCTGGTAATTACTTGTATT",
    "ATCACATGTAACTAATAATTTGAACTATACTTGAAGGACCGTGTTGATGTCAGGTATTTACAGTGGTTGGAAGATA",
    "GCAGTATTATTAGCATAAGCTGCATACGTAATATTCAGTAACTGCCATATTATATAACAAATTTACATTCACAAAT",
    "TCAGTATCCTGTTAAGTGTCATATTCTTGTAATCTGCATTCTCCAGGAGTTTTATGTGTTTAATAGATGAATTTAT",
    "TTTATTTCTAAAGGTATTCAAATGTTTCAGCACCATATAATAGAAATACCCAATTATATTCTAGTTCCTTTATGTC",
    "CTGTACATCATTCTCTGCTTGGATTTCC",
);

/// The last 100 bases of the first exon of [`SPLICED_LNCRNA`], `n.314` to
/// `n.413`, and the first 420 bases of the intron after it, the 520 bases at
/// 21:14961622-14962141.
const SPLICED_COPY: &str = concat!(
    "AAACCACAGACGTCATGTCCAGAAATGGAATACTGTCATGTTAGGATGAACATCTGATCTGATTTGTTAAGACATT",
    "TTGCCAAATGCACAGATGCACAGGGTCCCCCATGCACTTTTCATGCAGCAGGTAATCACAACTTAAATAAAAGGTT",
    "TATTCTATACATTTGTCCAGACTTGCAACTGTATACATACATGCACAACTTTTAAACCTGTCTGATTATATTTACA",
    "CTTATACATGGAATATACAGGAACCAAAGAGATTAAAAGGCTTTTCTGTTGCATTAGAACAATACAAAATATGTAT",
    "TTTTCATTAAGGAAATCACTATTTACATCACCTTTAAAAACCGATTTTAACATCTTCATGTAACAAGTCAATATAT",
    "ATATATATACATATTATATATATATATATATATATATATATAAAATATATACTCTCACCAGTTTACTCTTCTGTAA",
    "GATGCAGCAGAGAATAAACAGGTAATGAATGCTACCTTACTGATGTCAATAAGCTGTAGTACCC",
);

/// The 201st to the 740th base of the intron after `n.254` of
/// [`REVERSE_LNCRNA`] as the transcript reads them, the reverse complement of
/// the 540 bases at 21:14856802-14857341.
const INTRON_COPY: &str = concat!(
    "TACCAAACTCCCCAGATTTTAGAGTTGACAACTCTACTTATAAGTTATATAAAAGTCAGAGCACCTAATTATTGAC",
    "AAACTTCTCAATCCTTACCAGCTAGTTCCAGATGATAATGTGTAAACAATATTTCCTGGAATTATTTTGATTAGTA",
    "TTACACTTAAGGTTAAAGATTTTGTCAAGCCTAAGCATTACAGAAATTGTCATTCATGTGAAGTTATTCTCTTGTG",
    "CAATAGAAAAGCTTCATTTAGTGTGCTGAGTTTGTAGCTGCATTTGATGCTTAATGGGGTTAGCCCCGTGAAGCAT",
    "ATTGCAAGTGTCCTCTAGAGACTGGACTGGTTTCCAATTGATTTCAAGTTGCAGTTCAAATTGTTGTCTTTCACAA",
    "ATAGCAACTACATTGTGTCAATTGACGACCTCGCAAAAACTACTGTAAATAAAATCTTAAACTTTGAAACTACAGT",
    "AAAGGTTACTCTGTAATGTACGTGAGACAAGCAATCATGAATTAGGAGGGTTTCCTTTCAATTATGAGTCAAGCTT",
    "CAGGATCC",
);

/// The last 520 bases of [`LNCRNA`], `n.1281` to `n.1800`, at
/// 21:32280480-32280999.
const LNCRNA_END: &str = concat!(
    "CATAGAATTCCCTTCCTTCCTTCTGCTTTTTTTTTTCTTTTCTTTATTTTTTCAGAGATGATGTTTCAGTATGTTG",
    "CCCAGGCTGGTCTCGAACTCCTGGCCTCAAGCGATCCTCCTGCCTCAGCCTCCCAAGTAGCTGGAACTACAGGCAT",
    "GCCTCACCATACTTTGCCCCTTCTACTTTTGTAGAGCCAAGCCTGACCTTCCTGGACTGTTGTTTTCAAGCATAAT",
    "ATCTCTGGTACAACATTGATACAATCCACATGAGTTCTTTAAGTAAGTGGCTAGTTACTTGGGGATTCAAAGGATC",
    "ATCTTGTCTCAGACTTTGGGTTCCTATTGTCAATACTTAAATTCTATGATAGGATGCCCCAGGAAAGAAGGAATAC",
    "AAGGCCGAGCTGCGTGGTGACTTTCAGGGGCCCTATTTTTGCCTTAGTAGGCCTTTCTCTCCATAAAACAGTTAAT",
    "AAAAATTTCATTTTATAGTCACATTGCTATGAAGATGATTATATTAAAATGATATCATAAAACA",
);

/// The last 640 bases of [`REVERSE_CODING`], `c.*1863` to `c.*2502`, as the
/// transcript reads them, the reverse complement of 21:14371115-14371754.
const CODING_END: &str = concat!(
    "TGAATTTCTCAAATGTTTGAGATGTTCAAGACAGTCTCTTAATCCGTTAATGCTTTTGGAAACAATTGACAAAATA",
    "GGGCAGGCAGCTCATCTCATGTCCTGAAGTTGGAATTTAAATAATTCCTATTTGCAAATTAGAATGACAGTGTTGG",
    "AATTTGGAGGCAGTAGTTGAGCATATTCTCTAGTATATAGCTACACCTTTAATAAAATGAAGGAATGTCTTCAATC",
    "ATATTTTAGTGGGCTATTTATAAATAGTCTTGAAGTCAATTTAGTTTATTTATTTAAAAGATAATGCATCCTGAAA",
    "GGGATCATTTATGAATAACAATCTGAAGTCTTTTCATAAAAAAAATTAATAAACTTTAGTTGTACATTTAGCCAGT",
    "GTTATTTGAAGTATGTAACTTTTAAAATATTAAGTGTCTTGTATGATTAGAATATGTGAATGAGTAACTTATTTTG",
    "TATCAGGAATGTTTTGGTACTGTGTTTTCACTCAAACCACTGACTTAACAGATACTGCTGTGTATAACATGTACTA",
    "AATATTACAGTTATTGTGCATAACAGATTGTTCCTCTTATATTTGTGTGTATACAGGCAATTCATGTTTTAATGTA",
    "AATAAATACCATTTTGCAGTTTGTTTTTTAAA",
);

/// The last 520 bases of [`SPLICED_LNCRNA`], the last 371 bases of its intron
/// and its second exon, `n.414` to `n.562`, at 21:14963714-14964233.
const SPLICED_END: &str = concat!(
    "ATGCTGCCTATGAATCTGCATTTTAAAACAAGTGCTCCAGATGGTTCTAATGCAAGTGGTCTGAGTTCAGGTTTAA",
    "GGAATCCTGGTCATAACACTCGCAGGTACCACAGGAATGGTGAGAATTAACTTTTACAAGAGATGGTGAAAACTGT",
    "CCAAAGAAAAGCCATTTTGCTTCTGTTAAAAACAATTTCTTAAATATTCAACTTTTCAGGAAACCACCTAAAGATA",
    "AATGTTAACTTTTCTAATTCTGTAAAAGGAGGAGTATTAATTTAGTAAAGTAGTTAAAATACAAAAAATGGAAATC",
    "TGCCTTTTTTTTAAAGGATCTATCAAACTTCCCTCAATTGTAATCAATTTTTACTAAGTAGCACAAGGTTTAAAGA",
    "AATATATATTGAAAAGGATCAGGGTGGGACAGACACATTGGGAACCACAAAACCTCTTCCAAGAGTTAAAAATAAA",
    "CACCTTTTTTTTTTCTTTTTAGCAAAAGTGAAAAATGTAAGCATATATTTTTCTAATTTATTTT",
);

/// The last 520 bases of [`SHORT_END_LNCRNA`], the last 401 bases of its last
/// intron and its last exon, `n.634` to `n.752`, as the transcript reads
/// them, the reverse complement of 21:13943744-13944263.
const SHORT_END: &str = concat!(
    "AATTTTGCTTCTCTTTCTTATTGTTGTTCAGTTTGGCTCTTAATAATTTAGTGTTTGCCTAGTCTCTAGTTAATCT",
    "TCAGAAATTTACATGCACTGTAGGGGCTCACTCTTTCTGGTATGCCGAGGTAAACTCTTTGTAAGAGAGGAAGCTT",
    "TTATAATACTACCTATCATCTTTGAATTCATTTCTGGTAGACTTTACACATAATGCATTAAGTTTAGTCCAAACAA",
    "ACACTGAGAGTTCAGCTTGCCGGTTTATGTTTCTGTCCCGTGTTAAGCCAAGGCAAATTATTTTTCACTTTTTAGT",
    "TACAATCCCATAATTTAAGAGGGGCAACACACAGATTAAGTTTCACAGTTAAATTTTAATTATTTTCTAATATTTC",
    "TTTGTTTATACTTGATTAAAGCTAATTTTAGAACATGCACTCTGACAGAAAAGACATCTGAGAAACAAAACAAGCA",
    "AATTTGTTTTCCTTTTTGCACCTGCCAAAAAAAAAAAAAGAAAAGCCTCAAGAACCAGAACTGG",
);

/// The last 499 bases of [`FORWARD_CODING`], `c.*1091` to `c.*1589`, and the
/// 101 bases after the transcript, the 600 bases at 21:15879566-15880165.
const CODING_PAST_END: &str = concat!(
    "GTTTAAGATTAATGAATTAATATTGTGAATGAAAATCAAAATCCATACTTTAAAGGTAATCATGTTACTAACAACC",
    "TATTTTTGAATTCATAAAAATTTCTTTATAAATGATGTTTTGTGAACATAGTAAAATAGACCATTATACTATGTGT",
    "ATGTTTGATACAGCGTCGCCAAAACTAGTGTTCTTTATTAGTGCCTCTCACAAAAGATCCTGGATGGAGGAGTAAG",
    "ATGAAATATTATGCTATTATATGATGCTGTTTGTAAAGGTATTAATGTACTAGTAAGGTGTTAATGACAAGGAATT",
    "AGTACTATTCCTGTTGTAAAGTTAGATTTTGCATATTGTATCTATCAAAATATGTTTGGGTTTAGATTTTAAGTTG",
    "TCTACTGAGCAGATTTCTGCATTGGTTTTCCAGTCCTGTTAAAAGTTTAGAAACTTCATATGTGTCATCACAGCTT",
    "TTGTAAAGAAAGTATCCTTAATATTTTATGACATTCTACCACAGTGGTAAAGTTGTTTGTGTTTGTGTATAATTGA",
    "CTGCTGAGATACAAAGTGCCCAAACACTATATGGAAATTATTCCTTCAAAGTCACATCCAGAGAATCA",
);

/// The last 515 bases of [`REVERSE_LNCRNA`], `n.5415` to `n.5929`, and the 5
/// bases after the transcript as the transcript's strand reads them, the
/// reverse complement of the 520 bases at 21:14818838-14819357.
const LNCRNA_PAST_END: &str = concat!(
    "CCATTAATGTCTTTGCTCCAAAGTCTTCATTTGATTGTTTTTTTCTAACTTGTTATTTTTTTGAGCAAACATTCTT",
    "GTGAACTTTCATATTTTTAATCCTTCTTTATTGAGGCTTAAAAATATTTTTAAAAAATAGTAAAATGTATTTGCTA",
    "CAAAACTGAGCACTAACTTTTGAGAAGTTGTTCTAAGCTTTCATTTGTTTCTGGTCTGAATCATCAGCCACGAATT",
    "CTGAATAAGGCTGATCTTATATTATTGCCACGTTCAGGTCTCTAGAAATAACCATTATAGGACATGACCTATCATT",
    "CTTATGGATTTTATTGAATGATGAGCTCAGCTGCTTGTACAAATCATTTTATTTGCCTTCTCCTTTAATATATTTA",
    "TCAGCAAATTAGGAAGCAGTATTTCGTTTTGTGTTGTTACTCAAAATAAGAGGATTTTACATTTGTAAAATATGTT",
    "TTTATATCTTAGCCCATTATCTATTAAAAACAAAACCAATAAAATTTTGTCATTTATAACTATT",
);

/// The data lines of a `--hgvs --plugin HGVSIntronOffset` run on `records`
/// with `extra` flags, read on the corpus's cache and reference.
fn annotate(dir: &Path, name: &str, records: &str, extra: &[&str]) -> Vec<String> {
    let corpus = corpus();
    annotate_on(
        dir,
        name,
        records,
        &corpus.join("json_cache"),
        &corpus.join("reference.fa.gz"),
        extra,
    )
}

/// The data lines of a `--hgvs --plugin HGVSIntronOffset` run on `records`
/// with `extra` flags, read on `cache` and the reference `fasta`.
fn annotate_on(
    dir: &Path,
    name: &str,
    records: &str,
    cache: &Path,
    fasta: &Path,
    extra: &[&str],
) -> Vec<String> {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .arg("--fasta")
        .arg(fasta)
        .args(["--assembly", "GRCh38", "--offline", "--hgvs"])
        .args(["--plugin", "HGVSIntronOffset"])
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

/// The value of `key` in a default-format row's Extra column, `None` when the
/// column does not carry it.
fn extra_value<'a>(row: &'a str, key: &str) -> Option<&'a str> {
    let extra = row.split('\t').nth(13)?;
    extra
        .split(';')
        .find_map(|kv| kv.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
}

/// `[Consequence, HGVSc, HGVS_OFFSET, HGVS_IntronStartOffset,
/// HGVS_IntronEndOffset]` of `record`'s row on `transcript`, a field the row
/// does not carry left empty.
fn notation(rows: &[String], record: &str, transcript: &str) -> [String; 5] {
    let row = rows
        .iter()
        .find(|row| {
            let c: Vec<&str> = row.split('\t').collect();
            c[0] == record && c[4] == transcript
        })
        .unwrap_or_else(|| panic!("no row of {record} on {transcript}"));
    let extra = |key: &str| extra_value(row, key).unwrap_or_default();
    [
        row.split('\t').nth(6).unwrap_or_default(),
        extra("HGVSc"),
        extra("HGVS_OFFSET"),
        extra("HGVS_IntronStartOffset"),
        extra("HGVS_IntronEndOffset"),
    ]
    .map(str::to_string)
}

/// One expected row: the record's ID, the transcript and its [`notation`]
/// fields.
type Expected<'a> = (&'a str, &'a str, [&'a str; 5]);

/// Asserts each `expected` row of a run on `records` with `flags`.
fn assert_rows(records: &str, flags: &[&str], expected: &[Expected<'_>]) {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(dir.path(), "records", records, flags);
    for (record, transcript, fields) in expected {
        assert_eq!(
            notation(&rows, record, transcript),
            *fields,
            "{record} under {flags:?}"
        );
    }
}

/// The VCF line of the plain insertion `id`: the bases `added`, as the
/// chromosome's forward strand reads them, after the base `anchor` at
/// 21:`pos`.
fn insertion(pos: u32, id: &str, anchor: &str, added: &str) -> String {
    format!("21\t{pos}\t{id}\t{anchor}\t{anchor}{added}\t.\t.\t.\n")
}

/// The bases of a reverse-strand transcript as the chromosome's forward
/// strand reads them: `bases` reversed and complemented.
fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// `bases` as they stand after rolling `steps` bases 3': each step passes
/// their first base and moves it to their end.
fn rolled(bases: &str, steps: usize) -> String {
    let (passed, rest) = bases.split_at(steps % bases.len());
    format!("{rest}{passed}")
}

/// `bases` with the base at the 0-based `index` replaced by `base`.
fn with_base(bases: &str, index: usize, base: &str) -> String {
    format!("{}{base}{}", &bases[..index], &bases[index + 1..])
}

/// A record that adds after a base a copy of the bases that follow it has a
/// run of as many steps as it has bases, and one more for each base after the
/// copy that goes on repeating it.
/// `21 32279499` adds after `n.300` of the lncRNA a copy of `n.301` to
/// `n.820`, 520 bases inside its one exon, and `n.821` differs from `n.301`:
/// the last placement is after `n.820`, where the 520 bases before it are
/// the copy, so the row is `n.301_820dup` with `HGVS_OFFSET` 520, the bases
/// the notation moves. Ensembl VEP 116.3 stops after 1,001 less 520 steps
/// (TranscriptVariationAllele.pm:303-306) and writes the insertion between
/// `n.781` and `n.782` with `HGVS_OFFSET` 481.
/// On the reverse-strand coding transcript `21 14373443` adds after `c.*173`
/// a copy of `c.*174` to `c.*733`, 560 bases of the 3' UTR, and `c.*734` to
/// `c.*736` repeat `c.*174` to `c.*176`, so the run has 563 steps:
/// `c.*177_*736dup` with `HGVS_OFFSET` -563, negative on a reverse-strand
/// transcript, where Ensembl VEP 116.3 stops after 441 steps, between
/// `c.*614` and `c.*615`.
/// `--shift_3prime 1` annotates each record where that roll stops and writes
/// the insertion there. Under `--shift_hgvs 0` nothing rolls and the row is
/// the duplication of the bases the record copies, which follow its point:
/// Ensembl VEP 116.3 calls an insertion a duplication where the bases after
/// its point equal it, as it does where the bases before its point do.
#[test]
fn a_long_copy_inside_an_exon_is_a_duplication_at_its_last_placement() {
    let records = [
        insertion(32279499, "copy_exon_fwd", "G", LNCRNA_COPY),
        insertion(
            14373443,
            "copy_exon_rev",
            "T",
            &reverse_complement(UTR_COPY),
        ),
    ]
    .concat();
    assert_rows(
        &records,
        &[],
        &[
            (
                "copy_exon_fwd",
                LNCRNA,
                [EXON, "ENST00000783173.1:n.301_820dup", "520", "0", "0"],
            ),
            (
                "copy_exon_rev",
                REVERSE_CODING,
                [UTR, "ENST00000285667.4:c.*177_*736dup", "-563", "0", "0"],
            ),
        ],
    );
    let stopped_fwd = format!("ENST00000783173.1:n.781_782ins{}", rolled(LNCRNA_COPY, 481));
    let stopped_rev = format!("ENST00000285667.4:c.*614_*615ins{}", rolled(UTR_COPY, 441));
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "copy_exon_fwd",
                LNCRNA,
                [EXON, stopped_fwd.as_str(), "481", "0", "0"],
            ),
            (
                "copy_exon_rev",
                REVERSE_CODING,
                [UTR, stopped_rev.as_str(), "-441", "0", "0"],
            ),
        ],
    );
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "copy_exon_fwd",
                LNCRNA,
                [EXON, "ENST00000783173.1:n.301_820dup", "", "0", "0"],
            ),
            (
                "copy_exon_rev",
                REVERSE_CODING,
                [UTR, "ENST00000285667.4:c.*174_*733dup", "", "0", "0"],
            ),
        ],
    );
}

/// The duplicated bases need not lie in one exon, and the two fields of the
/// HGVSIntronOffset plugin carry the intron parts of the first and the last
/// duplicated base.
/// `21 14961621` adds after `n.313` of the spliced lncRNA a copy of the 520
/// bases that follow, the last 100 bases of its first exon, `n.314` to
/// `n.413`, and the first 420 bases of its intron: `n.314_413+420dup` with
/// `HGVS_OFFSET` 520, `HGVS_IntronStartOffset` 0 and `HGVS_IntronEndOffset`
/// 420. Ensembl VEP 116.3 stops after 481 steps, inside the intron, between
/// `n.413+381` and `n.413+382`; `--shift_3prime 1` annotates the record
/// there, so the row's Consequence is the intron's.
/// `21 14857341` adds, 200 bases into the intron after `n.254` of the
/// reverse-strand lncRNA, a copy of the 540 intron bases that follow:
/// `n.254+201_254+740dup` with `HGVS_OFFSET` -540 and the plugin's fields 201
/// and 740, where Ensembl VEP 116.3 stops after 461 steps, between
/// `n.254+661` and `n.254+662`.
/// Under `--shift_hgvs 0` each row is the duplication of the bases the record
/// copies, with no `HGVS_OFFSET`.
#[test]
fn a_long_copy_that_runs_through_an_intron_carries_the_intron_parts_of_its_ends() {
    let records = [
        insertion(14961621, "copy_across_fwd", "T", SPLICED_COPY),
        insertion(
            14857341,
            "copy_intron_rev",
            "A",
            &reverse_complement(INTRON_COPY),
        ),
    ]
    .concat();
    let across = "ENST00000446301.1:n.314_413+420dup";
    let within = "ENST00000432230.6:n.254+201_254+740dup";
    assert_rows(
        &records,
        &[],
        &[
            (
                "copy_across_fwd",
                SPLICED_LNCRNA,
                [EXON, across, "520", "0", "420"],
            ),
            (
                "copy_intron_rev",
                REVERSE_LNCRNA,
                [INTRON, within, "-540", "201", "740"],
            ),
        ],
    );
    let stopped_fwd = format!(
        "ENST00000446301.1:n.413+381_413+382ins{}",
        rolled(SPLICED_COPY, 481)
    );
    let stopped_rev = format!(
        "ENST00000432230.6:n.254+661_254+662ins{}",
        rolled(INTRON_COPY, 461)
    );
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "copy_across_fwd",
                SPLICED_LNCRNA,
                [INTRON, stopped_fwd.as_str(), "481", "381", "382"],
            ),
            (
                "copy_intron_rev",
                REVERSE_LNCRNA,
                [INTRON, stopped_rev.as_str(), "-461", "661", "662"],
            ),
        ],
    );
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "copy_across_fwd",
                SPLICED_LNCRNA,
                [EXON, across, "", "0", "420"],
            ),
            (
                "copy_intron_rev",
                REVERSE_LNCRNA,
                [INTRON, within, "", "201", "740"],
            ),
        ],
    );
}

/// A copy of a transcript's last bases added before them has its last
/// placement directly past the last base, where the bases before it are the
/// transcript's last bases: the row is their duplication.
/// `21 32280479` adds after `n.1280` of the lncRNA a copy of its last 520
/// bases, `n.1281` to `n.1800`, and the base after the transcript differs
/// from `n.1281`, so the run ends directly past `n.1800`: `n.1281_1800dup`
/// with `HGVS_OFFSET` 520. Ensembl VEP 116.3 stops after 481 steps, between
/// `n.1761` and `n.1762`.
/// `21 14371754` adds after `c.*1862` of the reverse-strand coding transcript
/// a copy of its last 640 bases, `c.*1863` to `c.*2502`: `c.*1863_*2502dup`
/// with `HGVS_OFFSET` -640, where Ensembl VEP 116.3 stops after 1,001 less
/// 640 steps, between `c.*2223` and `c.*2224`.
/// `--shift_3prime 1` writes the insertion where that roll stops. Under
/// `--shift_hgvs 0` the row is the same duplication with no `HGVS_OFFSET`:
/// the bases after the record's point are the ones it copies.
#[test]
fn a_long_copy_of_the_last_bases_is_a_duplication_of_the_last_bases() {
    let records = [
        insertion(32280479, "last_fwd", "C", LNCRNA_END),
        insertion(14371754, "last_rev", "A", &reverse_complement(CODING_END)),
    ]
    .concat();
    let last_fwd = "ENST00000783173.1:n.1281_1800dup";
    let last_rev = "ENST00000285667.4:c.*1863_*2502dup";
    assert_rows(
        &records,
        &[],
        &[
            ("last_fwd", LNCRNA, [EXON, last_fwd, "520", "0", "0"]),
            (
                "last_rev",
                REVERSE_CODING,
                [UTR, last_rev, "-640", "0", "0"],
            ),
        ],
    );
    let stopped_fwd = format!(
        "ENST00000783173.1:n.1761_1762ins{}",
        rolled(LNCRNA_END, 481)
    );
    let stopped_rev = format!(
        "ENST00000285667.4:c.*2223_*2224ins{}",
        rolled(CODING_END, 361)
    );
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "last_fwd",
                LNCRNA,
                [EXON, stopped_fwd.as_str(), "481", "0", "0"],
            ),
            (
                "last_rev",
                REVERSE_CODING,
                [UTR, stopped_rev.as_str(), "-361", "0", "0"],
            ),
        ],
    );
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            ("last_fwd", LNCRNA, [EXON, last_fwd, "", "0", "0"]),
            ("last_rev", REVERSE_CODING, [UTR, last_rev, "", "0", "0"]),
        ],
    );
}

/// Where the last exon is shorter than the copy, the duplication of the
/// transcript's last bases starts in the last intron, and
/// `HGVS_IntronStartOffset` carries the intron part of its first base.
/// `21 14963713` adds after `n.414-372` of the spliced lncRNA a copy of the
/// transcript's last 520 bases, 371 intron bases and the 149 bases of its
/// last exon, `n.414` to `n.562`: `n.414-371_562dup` with `HGVS_OFFSET` 520
/// and the plugin's fields -371 and 0. Ensembl VEP 116.3 stops after 481
/// steps, inside the last exon, between `n.523` and `n.524`;
/// `--shift_3prime 1` annotates the record there, so the row's Consequence
/// is the exon's.
/// `21 13944263` adds after `n.634-402` of the reverse-strand lncRNA whose
/// last exon is `n.634` to `n.752` a copy of its last 520 bases, 401 intron
/// bases and that exon: `n.634-401_752dup` with `HGVS_OFFSET` -520 and the
/// plugin's fields -401 and 0, where Ensembl VEP 116.3 stops between `n.713`
/// and `n.714`.
/// Under `--shift_hgvs 0` each row is the same duplication with no
/// `HGVS_OFFSET`.
#[test]
fn a_long_copy_of_the_last_bases_that_starts_in_the_last_intron_carries_its_intron_part() {
    let records = [
        insertion(14963713, "last_intron_fwd", "G", SPLICED_END),
        insertion(
            13944263,
            "last_intron_rev",
            "T",
            &reverse_complement(SHORT_END),
        ),
    ]
    .concat();
    let last_fwd = "ENST00000446301.1:n.414-371_562dup";
    let last_rev = "ENST00000429521.5:n.634-401_752dup";
    assert_rows(
        &records,
        &[],
        &[
            (
                "last_intron_fwd",
                SPLICED_LNCRNA,
                [INTRON, last_fwd, "520", "-371", "0"],
            ),
            (
                "last_intron_rev",
                SHORT_END_LNCRNA,
                [INTRON, last_rev, "-520", "-401", "0"],
            ),
        ],
    );
    let stopped_fwd = format!("ENST00000446301.1:n.523_524ins{}", rolled(SPLICED_END, 481));
    let stopped_rev = format!("ENST00000429521.5:n.713_714ins{}", rolled(SHORT_END, 481));
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "last_intron_fwd",
                SPLICED_LNCRNA,
                [EXON, stopped_fwd.as_str(), "481", "0", "0"],
            ),
            (
                "last_intron_rev",
                SHORT_END_LNCRNA,
                [EXON, stopped_rev.as_str(), "-481", "0", "0"],
            ),
        ],
    );
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "last_intron_fwd",
                SPLICED_LNCRNA,
                [INTRON, last_fwd, "", "-371", "0"],
            ),
            (
                "last_intron_rev",
                SHORT_END_LNCRNA,
                [INTRON, last_rev, "", "-401", "0"],
            ),
        ],
    );
}

/// A copy that takes in bases after the transcript has its last placement
/// after those bases, outside the transcript: the row carries no HGVSc, no
/// `HGVS_OFFSET` and no intron offsets, while its Consequence stays that of
/// the record's point.
/// `21 15879565` adds after `c.*1090` of the forward-strand coding
/// transcript a copy of the 600 bases that follow, its last 499 bases
/// `c.*1091` to `c.*1589` and the 101 bases after it, so the last placement
/// is 101 bases past the last base. Ensembl VEP 116.3 stops after 401 steps,
/// inside the transcript, and writes the insertion between `c.*1491` and
/// `c.*1492` with `HGVS_OFFSET` 401, the row under `--shift_3prime 1`.
/// `21 14819357` adds after `n.5414` of the reverse-strand lncRNA a copy of
/// its last 515 bases, `n.5415` to `n.5929`, and the 5 bases after it, where
/// Ensembl VEP 116.3 stops after 481 steps, between `n.5895` and `n.5896`.
/// Under `--shift_hgvs 0` the record stays at its point, where the bases
/// after it inside the transcript are fewer than the copy, so the row is the
/// insertion there.
#[test]
fn a_long_copy_whose_run_passes_the_last_base_has_no_hgvsc() {
    let records = [
        insertion(15879565, "past_fwd", "T", CODING_PAST_END),
        insertion(
            14819357,
            "past_rev",
            "G",
            &reverse_complement(LNCRNA_PAST_END),
        ),
    ]
    .concat();
    assert_rows(
        &records,
        &[],
        &[
            ("past_fwd", FORWARD_CODING, [UTR, "", "", "", ""]),
            ("past_rev", REVERSE_LNCRNA, [EXON, "", "", "", ""]),
        ],
    );
    let stopped_fwd = format!(
        "ENST00000285679.11:c.*1491_*1492ins{}",
        rolled(CODING_PAST_END, 401)
    );
    let stopped_rev = format!(
        "ENST00000432230.6:n.5895_5896ins{}",
        rolled(LNCRNA_PAST_END, 481)
    );
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "past_fwd",
                FORWARD_CODING,
                [UTR, stopped_fwd.as_str(), "401", "0", "0"],
            ),
            (
                "past_rev",
                REVERSE_LNCRNA,
                [EXON, stopped_rev.as_str(), "-481", "0", "0"],
            ),
        ],
    );
    let written_fwd = format!("ENST00000285679.11:c.*1090_*1091ins{CODING_PAST_END}");
    let written_rev = format!("ENST00000432230.6:n.5414_5415ins{LNCRNA_PAST_END}");
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "past_fwd",
                FORWARD_CODING,
                [UTR, written_fwd.as_str(), "", "0", "0"],
            ),
            (
                "past_rev",
                REVERSE_LNCRNA,
                [EXON, written_rev.as_str(), "", "0", "0"],
            ),
        ],
    );
}

/// A run that ends directly past the last base with less than a whole copy
/// before it duplicates nothing, and an insertion there has no base of the
/// transcript after it: the row carries no HGVSc, no `HGVS_OFFSET` and no
/// intron offsets.
/// `21 15879565` adds after `c.*1090` of the forward-strand coding
/// transcript its last 499 bases and 101 bases of which the first differs
/// from the base after the transcript, so the roll passes `c.*1091` to
/// `c.*1589`, 499 steps, and stops directly past the last base, where the
/// 600 bases before it are not the added bases. Ensembl VEP 116.3 stops
/// after 401 steps and writes the insertion between `c.*1491` and
/// `c.*1492`.
/// `21 14819357` adds after `n.5414` of the reverse-strand lncRNA its last
/// 515 bases and 5 bases of which the first differs from the base after the
/// transcript: 515 steps, where Ensembl VEP 116.3 stops after 481, between
/// `n.5895` and `n.5896`.
/// `--shift_3prime 1` writes the insertion where that roll stops, and
/// `--shift_hgvs 0` the insertion of the bases as the record writes them, at
/// its point.
#[test]
fn a_long_run_that_ends_at_the_last_base_short_of_a_whole_copy_has_no_hgvsc() {
    let edge_fwd = with_base(CODING_PAST_END, 499, "T");
    let edge_rev = with_base(LNCRNA_PAST_END, 515, "A");
    let records = [
        insertion(15879565, "edge_fwd", "T", &edge_fwd),
        insertion(14819357, "edge_rev", "G", &reverse_complement(&edge_rev)),
    ]
    .concat();
    assert_rows(
        &records,
        &[],
        &[
            ("edge_fwd", FORWARD_CODING, [UTR, "", "", "", ""]),
            ("edge_rev", REVERSE_LNCRNA, [EXON, "", "", "", ""]),
        ],
    );
    let stopped_fwd = format!(
        "ENST00000285679.11:c.*1491_*1492ins{}",
        rolled(&edge_fwd, 401)
    );
    let stopped_rev = format!("ENST00000432230.6:n.5895_5896ins{}", rolled(&edge_rev, 481));
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "edge_fwd",
                FORWARD_CODING,
                [UTR, stopped_fwd.as_str(), "401", "0", "0"],
            ),
            (
                "edge_rev",
                REVERSE_LNCRNA,
                [EXON, stopped_rev.as_str(), "-481", "0", "0"],
            ),
        ],
    );
    let written_fwd = format!("ENST00000285679.11:c.*1090_*1091ins{edge_fwd}");
    let written_rev = format!("ENST00000432230.6:n.5414_5415ins{edge_rev}");
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "edge_fwd",
                FORWARD_CODING,
                [UTR, written_fwd.as_str(), "", "0", "0"],
            ),
            (
                "edge_rev",
                REVERSE_LNCRNA,
                [EXON, written_rev.as_str(), "", "0", "0"],
            ),
        ],
    );
}

/// A copy with one base changed near its end rolls up to the changed base,
/// further than Ensembl VEP's roll reads and short of a whole copy, so
/// the bases before its last placement do not equal it: the row is the
/// insertion of the bases as they stand there.
/// `21 32279499` adds after `n.300` of the lncRNA the copy of `n.301` to
/// `n.820` with its 501st base changed: the roll passes `n.301` to `n.800`,
/// and the row is the insertion between `n.800` and `n.801` of the bases
/// rolled 500 steps, with `HGVS_OFFSET` 500, where Ensembl VEP 116.3 stops
/// after 481 steps, between `n.781` and `n.782`.
/// On the reverse-strand coding transcript `21 14373443` adds after `c.*173`
/// the copy of `c.*174` to `c.*733` with its 541st base changed: the
/// insertion between `c.*713` and `c.*714` with `HGVS_OFFSET` -540, where
/// Ensembl VEP 116.3 stops after 441 steps, between `c.*614` and `c.*615`.
/// `--shift_3prime 1` writes the insertion where that roll stops, and
/// `--shift_hgvs 0` the insertion of the bases as the record writes them, at
/// its point.
#[test]
fn a_long_copy_changed_near_its_end_is_an_insertion_at_its_last_placement() {
    let late_fwd = with_base(LNCRNA_COPY, 500, "C");
    let late_rev = with_base(UTR_COPY, 540, "A");
    let records = [
        insertion(32279499, "late_fwd", "G", &late_fwd),
        insertion(14373443, "late_rev", "T", &reverse_complement(&late_rev)),
    ]
    .concat();
    let last_fwd = format!("ENST00000783173.1:n.800_801ins{}", rolled(&late_fwd, 500));
    let last_rev = format!("ENST00000285667.4:c.*713_*714ins{}", rolled(&late_rev, 540));
    assert_rows(
        &records,
        &[],
        &[
            (
                "late_fwd",
                LNCRNA,
                [EXON, last_fwd.as_str(), "500", "0", "0"],
            ),
            (
                "late_rev",
                REVERSE_CODING,
                [UTR, last_rev.as_str(), "-540", "0", "0"],
            ),
        ],
    );
    let stopped_fwd = format!("ENST00000783173.1:n.781_782ins{}", rolled(&late_fwd, 481));
    let stopped_rev = format!("ENST00000285667.4:c.*614_*615ins{}", rolled(&late_rev, 441));
    assert_rows(
        &records,
        &["--shift_3prime", "1"],
        &[
            (
                "late_fwd",
                LNCRNA,
                [EXON, stopped_fwd.as_str(), "481", "0", "0"],
            ),
            (
                "late_rev",
                REVERSE_CODING,
                [UTR, stopped_rev.as_str(), "-441", "0", "0"],
            ),
        ],
    );
    let written_fwd = format!("ENST00000783173.1:n.300_301ins{late_fwd}");
    let written_rev = format!("ENST00000285667.4:c.*173_*174ins{late_rev}");
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "late_fwd",
                LNCRNA,
                [EXON, written_fwd.as_str(), "", "0", "0"],
            ),
            (
                "late_rev",
                REVERSE_CODING,
                [UTR, written_rev.as_str(), "", "0", "0"],
            ),
        ],
    );
}

/// A copy with one base changed halfway rolls up to the changed base, fewer
/// steps than Ensembl VEP's roll reads, so its last placement is where
/// Ensembl VEP 116.3 stops and the row is the insertion it writes there.
/// `21 32279499` adds after `n.300` of the lncRNA the copy of `n.301` to
/// `n.820` with its 261st base changed: the insertion between `n.560` and
/// `n.561` of the bases rolled 260 steps, with `HGVS_OFFSET` 260.
/// On the reverse-strand coding transcript `21 14373443` adds after `c.*173`
/// the copy of `c.*174` to `c.*733` with its 281st base changed: the
/// insertion between `c.*453` and `c.*454` with `HGVS_OFFSET` -280.
/// `--shift_3prime 1` writes the same rows, and `--shift_hgvs 0` the
/// insertion of the bases as the record writes them, at its point.
#[test]
fn a_long_copy_changed_halfway_is_an_insertion_where_its_shorter_run_ends() {
    let half_fwd = with_base(LNCRNA_COPY, 260, "G");
    let half_rev = with_base(UTR_COPY, 280, "T");
    let records = [
        insertion(32279499, "half_fwd", "G", &half_fwd),
        insertion(14373443, "half_rev", "T", &reverse_complement(&half_rev)),
    ]
    .concat();
    let last_fwd = format!("ENST00000783173.1:n.560_561ins{}", rolled(&half_fwd, 260));
    let last_rev = format!("ENST00000285667.4:c.*453_*454ins{}", rolled(&half_rev, 280));
    let last_placement: [Expected<'_>; 2] = [
        (
            "half_fwd",
            LNCRNA,
            [EXON, last_fwd.as_str(), "260", "0", "0"],
        ),
        (
            "half_rev",
            REVERSE_CODING,
            [UTR, last_rev.as_str(), "-280", "0", "0"],
        ),
    ];
    assert_rows(&records, &[], &last_placement);
    assert_rows(&records, &["--shift_3prime", "1"], &last_placement);
    let written_fwd = format!("ENST00000783173.1:n.300_301ins{half_fwd}");
    let written_rev = format!("ENST00000285667.4:c.*173_*174ins{half_rev}");
    assert_rows(
        &records,
        &["--shift_hgvs", "0"],
        &[
            (
                "half_fwd",
                LNCRNA,
                [EXON, written_fwd.as_str(), "", "0", "0"],
            ),
            (
                "half_rev",
                REVERSE_CODING,
                [UTR, written_rev.as_str(), "", "0", "0"],
            ),
        ],
    );
}

/// A short insertion in a short run rolls to the run's end, far fewer steps
/// than Ensembl VEP's roll reads, and is the duplication of the bases
/// before its last placement with the steps in `HGVS_OFFSET`.
/// `21 32280127` adds after `n.928` of the lncRNA the base that `n.929` to
/// `n.932` each are: `n.932dup`, offset 4. `21 32280216` adds after `n.1017`
/// a copy of `n.1018` to `n.1020`, whose first two bases `n.1021` and
/// `n.1022` repeat: `n.1020_1022dup`, offset 5.
/// On the reverse-strand coding transcript `21 14372422` adds after
/// `c.*1194` the base that `c.*1195` to `c.*1199` each are (`c.*1199dup`,
/// offset -5), and `21 14372372` adds after `c.*1244` a copy of `c.*1245`
/// and `c.*1246`, which `c.*1247` to `c.*1249` go on repeating
/// (`c.*1248_*1249dup`, offset -5).
/// `--shift_3prime 1` writes the same rows. Under `--shift_hgvs 0` each row
/// is the duplication of the bases directly after the record's point, with no
/// `HGVS_OFFSET`.
#[test]
fn a_short_insertion_in_a_short_run_is_a_duplication_at_the_run_end() {
    let records = concat!(
        "21\t32280127\tshort_one_fwd\tA\tAT\t.\t.\t.\n",
        "21\t32280216\tshort_three_fwd\tA\tAGTG\t.\t.\t.\n",
        "21\t14372422\tshort_one_rev\tA\tAA\t.\t.\t.\n",
        "21\t14372372\tshort_two_rev\tT\tTAT\t.\t.\t.\n",
    );
    let run_end: [Expected<'_>; 4] = [
        (
            "short_one_fwd",
            LNCRNA,
            [EXON, "ENST00000783173.1:n.932dup", "4", "0", "0"],
        ),
        (
            "short_three_fwd",
            LNCRNA,
            [EXON, "ENST00000783173.1:n.1020_1022dup", "5", "0", "0"],
        ),
        (
            "short_one_rev",
            REVERSE_CODING,
            [UTR, "ENST00000285667.4:c.*1199dup", "-5", "0", "0"],
        ),
        (
            "short_two_rev",
            REVERSE_CODING,
            [UTR, "ENST00000285667.4:c.*1248_*1249dup", "-5", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &run_end);
    assert_rows(records, &["--shift_3prime", "1"], &run_end);
    let as_written: [Expected<'_>; 4] = [
        (
            "short_one_fwd",
            LNCRNA,
            [EXON, "ENST00000783173.1:n.929dup", "", "0", "0"],
        ),
        (
            "short_three_fwd",
            LNCRNA,
            [EXON, "ENST00000783173.1:n.1018_1020dup", "", "0", "0"],
        ),
        (
            "short_one_rev",
            REVERSE_CODING,
            [UTR, "ENST00000285667.4:c.*1195dup", "", "0", "0"],
        ),
        (
            "short_two_rev",
            REVERSE_CODING,
            [UTR, "ENST00000285667.4:c.*1245_*1246dup", "", "0", "0"],
        ),
    ];
    assert_rows(records, &["--shift_hgvs", "0"], &as_written);
}

/// The bases before a made-up transcript on its chromosome: the transcript's
/// first base is base 1,001, so `c.N` of a first exon is base 1,000 + N, and
/// Ensembl VEP's roll has 1,000 bases to read on either side of a record.
const BEFORE: usize = 1000;

/// The first 30 bases of a made-up intron.
const INTRON_START: &str = "GTAAGTCCTTCCTTCCTTCCTTCCTTCCTT";

/// A made-up 3' UTR of 10 bases.
const UTR_BASES: &str = "CCTCTCCTCC";

/// 510 copies of `AG`, 1,020 bases. Read from its first base it is the 340
/// codons `AGA` and `GAG` in turn, arginine and glutamate 170 times.
fn ag_repeat() -> String {
    "AG".repeat(510)
}

/// A made-up coding transcript, `ENST00000000001.1` of the gene
/// `ENSG00000000001` with the protein `ENSP00000000001.1`, on the forward
/// strand of a made-up chromosome `1`.
struct MadeUp {
    /// The chromosome's bases, from base 1.
    chromosome: String,
    /// The transcript's record in the cache.
    model: serde_json::Value,
}

impl MadeUp {
    /// The transcript whose exons and the introns between them are `pieces`
    /// in turn, laid after 1,000 `C` and followed by `CCCCA` and 995 `C`.
    /// Its coding sequence opens its first exon and ends `utr` bases before
    /// its last base, and `peptide` is its protein.
    fn new(pieces: &[&str], utr: usize, peptide: &str) -> Self {
        let chromosome = format!(
            "{}{}CCCCA{}",
            "C".repeat(BEFORE),
            pieces.concat(),
            "C".repeat(995)
        );
        let spliced: String = pieces.iter().step_by(2).copied().collect();
        let (coding, three_prime_utr) = spliced.split_at(spliced.len() - utr);
        let (mut exons, mut introns, mut pairs) = (Vec::new(), Vec::new(), Vec::new());
        let (mut start, mut cdna) = (BEFORE + 1, 1);
        for (index, piece) in pieces.iter().enumerate() {
            let end = start + piece.len() - 1;
            if index % 2 == 0 {
                let last = cdna + piece.len() - 1;
                exons.push(json!({"start": start, "end": end}));
                pairs.push(json!({"from_start": cdna, "from_end": last,
                    "to_start": start, "to_end": end, "ori": 1}));
                cdna = last + 1;
            } else {
                introns.push(json!({"start": start, "end": end}));
            }
            start = end + 1;
        }
        let model = json!({
            "stable_id": "ENST00000000001",
            "version": 1,
            "gene_stable_id": "ENSG00000000001",
            "chr": "1",
            "start": BEFORE + 1,
            "end": start - 1,
            "strand": 1,
            "biotype": "protein_coding",
            "source": "ensembl",
            "protein": "ENSP00000000001",
            "exons": exons,
            "cdna_coding_start": 1,
            "cdna_coding_end": coding.len(),
            "translation": {"stable_id": "ENSP00000000001", "version": 1},
            "variation_effect_feature_cache": {
                "translateable_seq": coding,
                "peptide": peptide,
                "three_prime_utr": three_prime_utr,
                "introns": introns,
                "mapper": {"start_phase": 0, "cdna_coding_start": 1,
                    "cdna_coding_end": coding.len(), "pairs": pairs},
            },
        });
        Self { chromosome, model }
    }

    /// The VCF line of the plain insertion `id`: the bases `added` after the
    /// chromosome's base `after`.
    fn insertion(&self, after: usize, id: &str, added: &str) -> String {
        let anchor = &self.chromosome[after - 1..after];
        format!("1\t{after}\t{id}\t{anchor}\t{anchor}{added}\t.\t.\t.\n")
    }

    /// The data lines of a run of [`annotate_on`] on `records` with `extra`
    /// flags, read on a cache that holds the transcript alone and on the
    /// chromosome.
    fn annotate(&self, records: &str, extra: &[&str]) -> Vec<String> {
        let dir = tempfile::tempdir().unwrap();
        let cache = dir.path().join("cache");
        let shard = cache.join("transcripts").join("1").join("1-1000000.json");
        std::fs::create_dir_all(shard.parent().unwrap()).unwrap();
        std::fs::write(&shard, json!([self.model]).to_string()).unwrap();
        std::fs::write(
            cache.join("info.json"),
            json!({"assembly": "GRCh38"}).to_string(),
        )
        .unwrap();
        let fasta = dir.path().join("reference.fa");
        let length = self.chromosome.len();
        std::fs::write(&fasta, format!(">1\n{}\n", self.chromosome)).unwrap();
        std::fs::write(
            dir.path().join("reference.fa.fai"),
            format!("1\t{length}\t3\t{length}\t{}\n", length + 1),
        )
        .unwrap();
        annotate_on(dir.path(), "records", records, &cache, &fasta, extra)
    }
}

/// The default-format row of a record on the made-up transcript: `record` is
/// its ID, `location` its Location on chromosome `1`, `allele` its Allele,
/// `cells` its six cells from Consequence to Codons and `extra` its Extra
/// column.
fn made_up_row(record: &str, location: &str, allele: &str, cells: &str, extra: &str) -> String {
    format!(
        "{record}\t1:{location}\t{allele}\tENSG00000000001\tENST00000000001\tTranscript\t\
         {cells}\t-\t{extra}"
    )
}

/// The made-up transcript whose first exon ends in the repeat. Its first
/// exon is `ATG`, [`ag_repeat`] and an `A`, `c.1` to `c.1024` at
/// 1:1001-2024; its intron of 36 bases opens `GT`; its second exon is `CT`,
/// `CTG`, `AAA`, the stop codon `TAA` and the 3' UTR. The protein is
/// methionine, arginine and glutamate 170 times, then threonine 342 (`ACT`,
/// the first exon's last base and the second's first two), leucine 343 and
/// lysine 344.
fn repeat_to_the_end_of_a_coding_exon() -> MadeUp {
    let first_exon = format!("ATG{}A", ag_repeat());
    let intron = format!("{INTRON_START}TTTCAG");
    let second_exon = format!("CTCTGAAATAA{UTR_BASES}");
    MadeUp::new(
        &[&first_exon, &intron, &second_exon],
        UTR_BASES.len(),
        &format!("M{}TLK", "RE".repeat(170)),
    )
}

/// An insertion written in a coding exon whose run ends in the intron after
/// it is described there, and a change described in an intron has no
/// protein notation.
/// `1 1005` adds `AG` after `c.5` of [`repeat_to_the_end_of_a_coding_exon`]:
/// the roll passes `c.6` to `c.1024`, 1,019 bases, and the intron's `G`,
/// 1,020 steps, and stops before the intron's `T`. The two bases before
/// that point, `c.1024` and `c.1024+1`, are `AG`: the row is
/// `c.1024_1024+1dup` with `HGVS_OFFSET` 1020, `HGVS_IntronEndOffset` 1 and
/// no HGVSp.
/// Ensembl VEP 116.3 stops after 1,001 less 2 steps
/// (TranscriptVariationAllele.pm:303-306), between `c.1004` and `c.1005`,
/// and reads its HGVSp at that span (TranscriptVariationAllele.pm:1711-1738).
/// The other cells are those of the record's point, between the second and
/// the third base of codon 2, `AGA`.
#[test]
fn a_long_run_from_a_coding_exon_into_its_intron_has_no_hgvsp() {
    let transcript = repeat_to_the_end_of_a_coding_exon();
    let record = transcript.insertion(1005, "exon_to_intron", "AG");
    assert_eq!(
        transcript.annotate(&record, &[]),
        [made_up_row(
            "exon_to_intron",
            "1005-1006",
            "AG",
            "frameshift_variant\t5-6\t5-6\t2\tR/RX\taga/agAGa",
            "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.1024_1024+1dup;HGVS_OFFSET=1020;\
             HGVS_IntronEndOffset=1;HGVS_IntronStartOffset=0",
        )]
    );
}

/// Under `--shift_3prime 1` the record of
/// [`a_long_run_from_a_coding_exon_into_its_intron_has_no_hgvsp`] is
/// annotated where Ensembl VEP's roll stops and both notations are read
/// there, so the row keeps its HGVSp.
/// The roll stops after 999 steps, between `c.1004` and `c.1005`, the second
/// and the third base of codon 335, `GAG`, where the added bases stand as
/// `GA`. The two bases after that point, `c.1005` and `c.1006`, are `GA`:
/// `c.1005_1006dup` with `HGVS_OFFSET` 999.
/// With one more `AG` the repeat ends on base 1,025 of the coding sequence:
/// codons 2 to 341 read as before, codon 342 is `AGA`, arginine in place of
/// threonine, codon 343 is `CTC`, leucine, and codon 344 is `TGA`, the third
/// codon from the changed residue: `p.Thr342ArgfsTer3`.
#[test]
fn shift_3prime_reads_both_notations_of_a_long_run_where_the_roll_stops() {
    let transcript = repeat_to_the_end_of_a_coding_exon();
    let record = transcript.insertion(1005, "exon_to_intron", "AG");
    assert_eq!(
        transcript.annotate(&record, &["--shift_3prime", "1"]),
        [made_up_row(
            "exon_to_intron",
            "1005-1006",
            "GA",
            "frameshift_variant\t1004-1005\t1004-1005\t335\tE/EX\tgag/gaGAg",
            "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.1005_1006dup;\
             HGVSp=ENSP00000000001.1:p.Thr342ArgfsTer3;HGVS_OFFSET=999;\
             HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0",
        )]
    );
}

/// A frameshift whose run ends inside the coding sequence carries one HGVSp
/// wherever in the run it is read, since every placement gives the same
/// coding sequence: this row holds the notation of an insertion that the
/// roll leaves between two coding bases.
/// The made-up transcript's one exon is `ATG`, [`ag_repeat`], `CTG`, `AAA`,
/// the stop codon `TAA` and the 3' UTR: methionine, arginine and glutamate
/// 170 times, leucine 342 and lysine 343. `1 1005` adds `AG` after `c.5`:
/// the roll passes `c.6` to `c.1023`, 1,018 steps, and stops before the `C`
/// of `c.1024`, where the two bases before it are `AG`: `c.1022_1023dup`
/// with `HGVS_OFFSET` 1018.
/// With one more `AG` the repeat ends on base 1,025 of the coding sequence:
/// codons 2 to 341 read as before, codon 342 is `AGC`, serine in place of
/// leucine, and codon 343 is `TGA`, the second codon from the changed
/// residue: `p.Leu342SerfsTer2`.
/// Ensembl VEP 116.3 stops after 999 steps, between `c.1004` and `c.1005`,
/// 19 steps short, and the notation read there is the same.
#[test]
fn a_long_run_inside_the_coding_sequence_keeps_its_frameshift_hgvsp() {
    let exon = format!("ATG{}CTGAAATAA{UTR_BASES}", ag_repeat());
    let transcript = MadeUp::new(
        &[&exon],
        UTR_BASES.len(),
        &format!("M{}LK", "RE".repeat(170)),
    );
    let record = transcript.insertion(1005, "inside_the_coding_sequence", "AG");
    assert_eq!(
        transcript.annotate(&record, &[]),
        [made_up_row(
            "inside_the_coding_sequence",
            "1005-1006",
            "AG",
            "frameshift_variant\t5-6\t5-6\t2\tR/RX\taga/agAGa",
            "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.1022_1023dup;\
             HGVSp=ENSP00000000001.1:p.Leu342SerfsTer2;HGVS_OFFSET=1018;\
             HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0",
        )]
    );
}

/// An insertion written in an intron whose run ends between two coding bases
/// of the next exon is described there and carries the protein notation of
/// the insertion at that point, as one whose shorter run ends there does.
/// The made-up transcript's first exon is `ATG` and nine copies of `GCT`,
/// `c.1` to `c.30` at 1:1001-1030; its intron is [`INTRON_START`] and
/// [`ag_repeat`], 1:1031-2080; its second exon is `AGT`, `CTA`, `AAA`,
/// `GCT`, the stop codon `TGA` and the 3' UTR, `c.31` at 1:2081. The protein
/// is methionine, nine alanines, serine 11, leucine 12, lysine 13 and
/// alanine 14.
/// `1 1060` adds `AG` after the intron's 30th base: the roll passes the
/// 1,020 bases of the repeat, `c.31` and `c.32`, 1,022 steps, and stops
/// before the `T` of `c.33`. The two bases before that point are `c.31` and
/// `c.32`: `c.31_32dup` with `HGVS_OFFSET` 1022. With `AG` between `c.32`
/// and `c.33` the coding sequence reads `AGA`, `GTC`, `TAA` from codon 11:
/// arginine in place of serine 11, valine, and a stop codon, the third codon
/// from the changed residue: `p.Ser11ArgfsTer3`.
/// The same point with `AGAGAG` added has the same run. The six bases before
/// its end are the intron's last four, `c.31` and `c.32`: `c.31-4_32dup`
/// with `HGVS_IntronStartOffset` -4. The coding sequence reads `AGA`, `GAG`,
/// `AGT`, `CTA`, `AAA` from codon 11: arginine and glutamate, then serine,
/// leucine and lysine as before, two residues between alanine 10 and serine
/// 11: `p.Ala10_Ser11insArgGlu`.
/// Ensembl VEP 116.3 stops after 999 and after 995 steps, 21 and 25 bases
/// before the intron's end, and reads its HGVSp at that span
/// (TranscriptVariationAllele.pm:1711-1738). Each row's other cells are
/// those of the record's point, in the intron.
#[test]
fn a_long_run_from_an_intron_into_the_next_exon_carries_the_hgvsp_of_its_last_placement() {
    let first_exon = format!("ATG{}", "GCT".repeat(9));
    let intron = format!("{INTRON_START}{}", ag_repeat());
    let second_exon = format!("AGTCTAAAAGCTTGA{UTR_BASES}");
    let transcript = MadeUp::new(
        &[&first_exon, &intron, &second_exon],
        UTR_BASES.len(),
        &format!("M{}SLKA", "A".repeat(9)),
    );
    let records = [
        transcript.insertion(1060, "intron_to_exon_frameshift", "AG"),
        transcript.insertion(1060, "intron_to_exon_in_frame", "AGAGAG"),
    ]
    .concat();
    let intron_cells = "intron_variant\t-\t-\t-\t-\t-";
    assert_eq!(
        transcript.annotate(&records, &[]),
        [
            made_up_row(
                "intron_to_exon_frameshift",
                "1060-1061",
                "AG",
                intron_cells,
                "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000001.1:c.31_32dup;\
                 HGVSp=ENSP00000000001.1:p.Ser11ArgfsTer3;HGVS_OFFSET=1022;\
                 HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0",
            ),
            made_up_row(
                "intron_to_exon_in_frame",
                "1060-1061",
                "AGAGAG",
                intron_cells,
                "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000001.1:c.31-4_32dup;\
                 HGVSp=ENSP00000000001.1:p.Ala10_Ser11insArgGlu;HGVS_OFFSET=1022;\
                 HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=-4",
            ),
        ]
    );
}

/// The made-up transcript of one exon whose coding sequence is `ATG`, 180
/// copies of `GCT` and the stop codon `TAA`, `c.1` to `c.546` at
/// 1:1001-1546, with a 3' UTR of 20 `C`, `c.*1` to `c.*20` at 1:1547-1566.
/// Its protein is methionine and 180 alanines, and the stop codon is codon
/// 182. The chromosome goes on `CCCCA` after the transcript's last base.
fn alanine_codons_to_a_short_utr() -> MadeUp {
    let exon = format!("ATG{}TAA{}", "GCT".repeat(180), "C".repeat(20));
    MadeUp::new(&[&exon], 20, &format!("M{}", "A".repeat(180)))
}

/// The cells from Consequence to Codons of a record that adds `copy` after
/// coding base `after` of [`alanine_codons_to_a_short_utr`], the second base
/// of the alanine codon `codon`: the codon with `copy` inside it reads as
/// `residues` and one base more.
fn cells_of_a_copy_in_an_alanine_codon(
    after: usize,
    codon: usize,
    residues: &str,
    copy: &str,
) -> String {
    let span = format!("{after}-{}", after + 1);
    format!(
        "inframe_insertion,stop_retained_variant\t{span}\t{span}\t{codon}\tA/{residues}X\t\
         gct/gc{copy}t"
    )
}

/// An insertion written in the coding sequence whose run ends in the 3' UTR
/// is described there, and a change described in a UTR has no protein
/// notation.
/// `1 1038` adds after `c.38` of [`alanine_codons_to_a_short_utr`], the
/// second base of codon 13, a copy of the 520 bases that follow: `c.39`, the
/// 168 codons `c.40` to `c.543`, the stop codon and `c.*1` to `c.*12`.
/// `c.*13` is `C` where `c.39` is `T`, so the run has 520 steps and ends
/// between `c.*12` and `c.*13`, where the 520 bases before it are the copy:
/// `c.39_*12dup` with `HGVS_OFFSET` 520 and no HGVSp.
/// Ensembl VEP 116.3 stops after 1,001 less 520 steps
/// (TranscriptVariationAllele.pm:303-306), between `c.519` and `c.520`, 24
/// bases before the stop codon, and reads its HGVSp at that span
/// (TranscriptVariationAllele.pm:1711-1738).
/// The other cells are those of the record's point: codon 13 with the copy
/// inside it reads `GCT` 169 times, the copy's stop codon, its 12 `C` as
/// four proline codons and the codon's own `T`.
#[test]
fn a_long_copy_from_the_coding_sequence_into_the_3_prime_utr_has_no_hgvsp() {
    let transcript = alanine_codons_to_a_short_utr();
    let copy = &transcript.chromosome[1038..1558];
    let record = transcript.insertion(1038, "coding_to_utr", copy);
    let residues = format!("{}*PPPP", "A".repeat(169));
    assert_eq!(
        transcript.annotate(&record, &[]),
        [made_up_row(
            "coding_to_utr",
            "1038-1039",
            copy,
            &cells_of_a_copy_in_an_alanine_codon(38, 13, &residues, copy),
            "IMPACT=MODERATE;STRAND=1;HGVSc=ENST00000000001.1:c.39_*12dup;HGVS_OFFSET=520;\
             HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0",
        )]
    );
}

/// An insertion written in the coding sequence whose run goes on past the
/// transcript's last base has no description on the transcript: the row
/// carries no HGVSc, no HGVSp, no `HGVS_OFFSET` and no intron offsets.
/// `1 1050` adds after `c.50` of [`alanine_codons_to_a_short_utr`], the
/// second base of codon 17, a copy of the 520 bases that follow: `c.51`, the
/// 164 codons `c.52` to `c.543`, the stop codon, the 20 bases of the 3' UTR
/// and the 4 bases after the transcript. The fifth base after the transcript
/// is `A` where `c.51` is `T`, so the run has 520 steps and its last
/// placement lies 4 bases past the last base.
/// Ensembl VEP 116.3 stops after 481 steps, between `c.531` and `c.532`, 12
/// bases before the stop codon, and reads both notations at that span
/// (TranscriptVariationAllele.pm:1441-1483, 1711-1738).
/// The other cells are those of the record's point: codon 17 with the copy
/// inside it reads `GCT` 165 times, the copy's stop codon, its 24 `C` as
/// eight proline codons and the codon's own `T`.
#[test]
fn a_long_copy_from_the_coding_sequence_past_the_last_base_has_no_hgvs_field() {
    let transcript = alanine_codons_to_a_short_utr();
    let copy = &transcript.chromosome[1050..1570];
    let record = transcript.insertion(1050, "coding_past_the_last_base", copy);
    let residues = format!("{}*PPPPPPPP", "A".repeat(165));
    assert_eq!(
        transcript.annotate(&record, &[]),
        [made_up_row(
            "coding_past_the_last_base",
            "1050-1051",
            copy,
            &cells_of_a_copy_in_an_alanine_codon(50, 17, &residues, copy),
            "IMPACT=MODERATE;STRAND=1",
        )]
    );
}
