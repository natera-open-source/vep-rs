// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads HGVS protein notation: a notation resolves to the
//! nucleotide change Ensembl VEP's protein branch resolves it to, on the
//! transcript's strand; a gene name as the reference is tried over the gene's
//! transcripts with the gene-identifier warning; a notation the resolver
//! refuses is skipped with its message.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-input-hgvs")
}

/// The rows of a successful `--format hgvs` run on `input` against the corpus
/// cache and reference, written to `out`.
fn annotate(dir: &Path, input: &Path, out: &Path) -> Vec<String> {
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args(["--assembly", "GRCh38", "--format", "hgvs"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}: {}",
        input.display(),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&std::fs::read(out).unwrap())
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The rows grouped by `Uploaded_variation`, each row past that column.
fn rows_by_name(rows: &[String]) -> BTreeMap<&str, Vec<&str>> {
    let mut grouped: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for row in rows {
        let (name, rest) = row.split_once('\t').unwrap();
        grouped.entry(name).or_default().push(rest);
    }
    grouped
}

fn warnings_beside(out: &Path) -> String {
    std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default()
}

/// The corpus's protein notations annotate at the nucleotide change Ensembl
/// VEP 116.2 resolves them to in database mode, every row's `Location` and
/// `Allele` as it writes them: a one-base change for a residue a base away
/// from the written one, a two-base change for tryptophan (`Phe10Trp`), the
/// whole codon when the first and third base change (`Arg12Trp`, `TGG` carried
/// with the unchanged middle base), the synonymous base (`Phe10=`), a stop
/// gained and lost, a `delins` with and without the written residue, the
/// codons of a deleted residue and of a deleted range; on the reverse-strand
/// TFF3 transcript the allele is the transcript's base (`Ala2Gly` writes `G`,
/// its VCF twin `C`), the codon spanning its first intron places `Ser28Ala` at
/// the exon's last base, and a translation whose coding sequence opens in
/// phase 2 (`ENSP00000381462`) places its second residue with that shift.
#[test]
fn protein_notations_resolve_to_ensembl_veps_nucleotide_changes() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(
        dir.path(),
        &corpus().join("variants.hgvs.txt"),
        &dir.path().join("corpus.txt"),
    );
    let by_name = rows_by_name(&rows);
    let placed = [
        ("ENSP00000291554:p.Thr4Ala", "21:43169109", "G"),
        ("ENSP00000291554.2:p.Ile5Val", "21:43169112", "G"),
        ("ENSP00000291554.9:p.Q6E", "21:43169115", "G"),
        ("ENST00000291554:p.His7Tyr", "21:43169118", "T"),
        ("ENSP00000291554:p.Phe10Trp", "21:43169128-43169129", "GG"),
        ("ENSP00000291554:p.Arg12Trp", "21:43169133-43169135", "TGG"),
        ("ENSP00000291554:p.Phe10=", "21:43169129", "T"),
        ("ENSP00000291554:p.Gln6Ter", "21:43169115", "T"),
        ("ENSP00000291554:p.Lys11*", "21:43169130", "T"),
        ("ENSP00000291554:p.Ter174Gln", "21:43172278", "C"),
        ("ENSP00000291554:p.4delinsAla", "21:43169109", "G"),
        ("ENSP00000291554:p.Thr13delinsPro", "21:43169136", "C"),
        ("ENSP00000291554:p.Leu14del", "21:43169139-43169141", "-"),
        (
            "ENSP00000291554:p.Gly15_Pro16del",
            "21:43169142-43169147",
            "-",
        ),
        ("ENSP00000430690:p.Ala2Gly", "21:42315370", "G"),
        ("ENSP00000430690:p.Arg4Ter", "21:42315365", "T"),
        ("ENSP00000430690:p.Cys7=", "21:42315354", "T"),
        (
            "ENSP00000430690:p.Met8_Leu9del",
            "21:42315348-42315353",
            "-",
        ),
        ("ENSP00000430690:p.Ser28Ala", "21:42315293", "G"),
        ("ENSP00000381462:p.Ala2Gly", "21:42315334", "G"),
    ];
    for (notation, location, allele) in placed {
        let notation_rows = by_name
            .get(notation)
            .unwrap_or_else(|| panic!("{notation} has no rows"));
        let location_and_allele = format!("{location}\t{allele}\t");
        assert!(
            notation_rows
                .iter()
                .all(|row| row.starts_with(&location_and_allele)),
            "{notation}: {notation_rows:?}"
        );
    }
    assert!(by_name["ENSP00000291554:p.Ter174Gln"]
        .iter()
        .any(|row| row.contains("\tENST00000291554\tTranscript\tstop_lost\t")));
    assert!(by_name["ENSP00000291554:p.Gln6Ter"]
        .iter()
        .any(|row| row.contains("\tENST00000291554\tTranscript\tstop_gained\t")));
    assert!(by_name["ENSP00000291554:p.Phe10="]
        .iter()
        .any(|row| row.contains("\tENST00000291554\tTranscript\tsynonymous_variant\t")));
    assert!(by_name["ENSP00000430690:p.Met8_Leu9del"]
        .iter()
        .any(|row| row.contains("\tENST00000518498\tTranscript\tinframe_deletion\t")));
    assert!(warnings_beside(&dir.path().join("corpus.txt")).is_empty());
}

/// A gene name as the reference tries every transcript of the gene in turn and
/// takes the first that yields a variant (`VariationFeatureAdaptor.pm:1933`,
/// 1946-1985), with the gene-identifier warning (Parser/HGVS.pm:219-229):
/// Ensembl VEP 116.2 in database mode writes `CRYAA:p.Thr4Ala` at 21:43169109
/// with the allele `G`, `CRYAA:p.Met1Ala` at 21:43169100-43169101 with `GC`
/// and `TFF3:p.Ala2Gly` at 21:42315370 with `G`, the rows of the translations'
/// own notations; the warning is written once per gene name, as Ensembl VEP's
/// sink writes each distinct warning once. A qualifier before the colon and a
/// `(p.=)` suffix are removed, a wrongly written residue on a deletion is not
/// checked, and a lower-case three-letter new residue is read.
#[test]
fn gene_names_and_qualifiers_resolve_as_ensembl_vep_resolves_them() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("genes.hgvs");
    std::fs::write(
        &input,
        "CRYAA:p.Thr4Ala\nCRYAA:p.Met1Ala\nTFF3:p.Ala2Gly\nENSP00000291554(CRYAA):p.Thr4Ala\n\
         ENSP00000291554:p.Thr4Ala(p.=)\nENSP00000291554:p.Ala4del\nENSP00000291554:p.Thr4ala\n\
         ENSP00000291554:p.Thr4Ala\nENSP00000430690:p.Ala2Gly\n",
    )
    .unwrap();
    let out = dir.path().join("genes.txt");
    let rows = annotate(dir.path(), &input, &out);
    let by_name = rows_by_name(&rows);
    assert_eq!(by_name.len(), 9, "{rows:?}");
    for (notation, same_as) in [
        ("CRYAA:p.Thr4Ala", "ENSP00000291554:p.Thr4Ala"),
        ("TFF3:p.Ala2Gly", "ENSP00000430690:p.Ala2Gly"),
        (
            "ENSP00000291554(CRYAA):p.Thr4Ala",
            "ENSP00000291554:p.Thr4Ala",
        ),
        (
            "ENSP00000291554:p.Thr4Ala(p.=)",
            "ENSP00000291554:p.Thr4Ala",
        ),
        ("ENSP00000291554:p.Thr4ala", "ENSP00000291554:p.Thr4Ala"),
    ] {
        assert_eq!(by_name[notation], by_name[same_as], "{notation}");
    }
    assert!(by_name["CRYAA:p.Met1Ala"]
        .iter()
        .all(|row| row.starts_with("21:43169100-43169101\tGC\t")));
    assert!(by_name["CRYAA:p.Met1Ala"]
        .iter()
        .any(|row| row.contains("\tENST00000291554\tTranscript\tstart_lost\t")));
    assert!(by_name["ENSP00000291554:p.Ala4del"]
        .iter()
        .all(|row| row.starts_with("21:43169109-43169111\t-\t")));
    let warnings = warnings_beside(&out);
    let gene_warning = |gene: &str| {
        format!(
            "WARNING: Possible invalid use of gene or protein identifier '{gene}' as HGVS reference; \
             most likely transcript will be selected\n"
        )
    };
    assert_eq!(
        warnings,
        format!("{}{}", gene_warning("CRYAA"), gene_warning("TFF3"))
    );
}

/// A protein notation Ensembl VEP's resolver refuses is skipped with its
/// message, the text Ensembl VEP 116.2 writes in database mode for the same
/// line: a residue change with two equally short codon paths (`Asp2Glu`, a
/// base from `GAA` and from `GAG`), a frameshift, a residue no codon codes
/// (`fs`, `?`, `X`), a synonymous change whose only codon is the reference's,
/// a translation the cache lacks, a deletion whose codon crosses an exon
/// boundary, a description without a position or with position 0, a residue
/// written in lower case, a range written as a `delins`. A position past the
/// transcript, where Ensembl VEP's resolver dies on the mapper's gap, is skipped
/// as a sequence that cannot be found, and a deletion range written end
/// first, for which Ensembl VEP writes nothing, likewise. A RefSeq protein
/// or transcript id and a UniProt accession, which Ensembl VEP reads through
/// its external-name tables, are refused as references the cache does not
/// carry. The notations around them annotate.
#[test]
fn unresolvable_protein_notations_are_skipped_with_the_warning() {
    let dir = tempfile::tempdir().unwrap();
    let refused = [
        (
            "ENSP00000291554:p.Asp2Glu",
            "Could not uniquely determine nucleotide change from ENSP00000291554:p.Asp2Glu",
        ),
        (
            "ENSP00000291554:p.Arg12Glyfs*3",
            "Frameshifts are not supported for HGVS protein input",
        ),
        (
            "ENSP00000291554:p.Arg12fs",
            "Could not determine nucleotide change from peptide change R -> fs",
        ),
        (
            "ENSP00000291554:p.Met1?",
            "Could not determine nucleotide change from peptide change M -> ?",
        ),
        (
            "ENSP00000291554:p.Thr4Xaa",
            "Could not determine nucleotide change from peptide change T -> X",
        ),
        (
            "ENSP00000291554:p.Met1=",
            "Reference allele extracted from ENSP00000291554:43169100-43169102 (ATG) matches alt \
             allele given by HGVS notation ENSP00000291554:p.Met1= (ATG)",
        ),
        (
            "ENSP00009999999:p.Thr4Ala",
            "Could not get a Transcript object for 'ENSP00009999999'",
        ),
        (
            "ENSP00000430690:p.Ser28del",
            "Reference allele extracted from ENSP00000430690:42315293-42315293 (T) does not match \
             reference allele given by HGVS notation ENSP00000430690:p.Ser28del (TCT)",
        ),
        ("ENSP00000291554:p.Thr4>Ala", "Must call with start and end"),
        (
            "ENSP00000291554:p.Thr0Ala",
            "Could not parse HGVS protein notation ENSP00000291554:p.Thr0Ala",
        ),
        (
            "ENSP00000291554:p.1Ala",
            "Could not parse HGVS protein notation ENSP00000291554:p.1Ala",
        ),
        (
            "ENSP00000291554:p.thr4ala",
            "Sequence translated from reference (ACC -> T) does not match input sequence (thr)",
        ),
        (
            "ENSP00000291554:p.Thr4_Ile5delinsAla",
            "Could not determine nucleotide change from peptide change T -> _Ile5delinsAla",
        ),
        (
            "ENSP00000291554:p.Thr999Ala",
            "Unable to find the genomic reference sequence for protein ENSP00000291554",
        ),
        (
            "ENSP00000291554:p.Ile5_Thr4del",
            "Unable to find the genomic reference sequence for protein ENSP00000291554",
        ),
        (
            "NP_000385.1:p.Thr4Ala",
            "Could not get a Transcript object for 'NP_000385.1'",
        ),
        (
            "NM_000394.4:p.Thr4Ala",
            "Could not get a Transcript object for 'NM_000394.4'",
        ),
        (
            "P02489:p.Thr4Ala",
            "Could not get a Transcript object for 'P02489'",
        ),
    ];
    let mut text = String::from("ENSP00000291554:p.Thr4Ala\n");
    for (notation, _) in refused {
        text.push_str(notation);
        text.push('\n');
    }
    text.push_str("ENSP00000430690:p.Ala2Gly\n");
    let input = dir.path().join("refused.hgvs");
    std::fs::write(&input, &text).unwrap();
    let out = dir.path().join("refused.txt");
    let rows = annotate(dir.path(), &input, &out);
    let names: Vec<&str> = rows_by_name(&rows).into_keys().collect();
    assert_eq!(
        names,
        ["ENSP00000291554:p.Thr4Ala", "ENSP00000430690:p.Ala2Gly"]
    );
    let expected: String = refused
        .iter()
        .enumerate()
        .map(|(i, (notation, message))| {
            format!(
                "WARNING: line {} skipped ({notation}): Unable to parse HGVS notation '{notation}': {message}\n",
                i + 2
            )
        })
        .collect();
    assert_eq!(warnings_beside(&out), expected);
}
