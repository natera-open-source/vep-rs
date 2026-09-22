// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Real Ensembl prediction matrices decoded under the Perl packing rule.
//!
//! Two `protein_function_predictions.prediction_matrix` blobs from the public Ensembl
//! variation database (`homo_sapiens_variation_115_38`), chosen for small peptides that
//! carry several prediction categories, are embedded here gzip-compressed and
//! base64-encoded. Every cell they hold was independently unpacked by the rule of
//! `ProteinFunctionPredictionMatrix.pm` (ensembl-variation release/115; `unpack 'v'`,
//! category `>> 14`, score `& 0x3FF`, then `/ 1000`) to produce the expected counts
//! and cells below, so a decoder that read the category from any other bit position
//! fails here on real data.

use std::collections::BTreeMap;

use vep_core::prediction::{decode_prediction, format_prediction, gunzip_matrix, AnalysisType};

/// SIFT matrix for translation md5 `685580feaef9fc71e213830bdd64ceab` (39 residues).
const SIFT_B64: &str = "H4sIAAAAAAAA/51UoW4DMQx1ekUBQwFF0TQwBR8YOhCy4qKBA0MFw0OFU7dPmAqG9iVT0bQvGCoYvG8oyRw7uUvSu2qarFwiO7afn527u12BTUXhkta5UgeJTmYe7+C/S1C2ptMKokUXO0eJduek9dKAwQwNSLsDjSdp1xhrCwuxQV1N3rpAcE4aiJi7aqwOlgdx6nkuvrdtQGe6odI8h3NNZkm9TJ/HTGQyASfH6CqwN+JbDPFrYP/XmaRecJyXStoP4O552cGPmKqnxMfSVUvwnhp9IcRpUaPo9Dnz/UrvbyDy6Nw6iacRmd/r7N694I6nEba9Vzlh/5XTCJfI36MYu3WEBvny0lW+akMz6Jwm3ZswyIMkq0GkOvgxH75K5nkuLgTzFetd9n0reVbJN8pceB5SXRv6ewDOwVqPJeaU9EKvxBdMVT3OjS5ebir+1eXTzZi4bq7DZJap+T2PxvsuaDJVQDON6W+xU931LNebPnoLQH0d6pXUI35nA4aWZkD1fT5CVzEvfBvCX47feY4952/AleuZubrQ8bzks9FSJrOXe4ULcMlwlntN++Gpe0YvssGJIIIR7V/tvzlxgCQbBgAA";

/// PolyPhen (HumVar) matrix for translation md5 `77d96fc8e5c080038b043ead02dadfc3` (2 residues).
const POLYPHEN_B64: &str = "H4sIAAAAAAAA/wtzDXjG/AIInzLLOD9nfsJs7HySua/h/3+Q2CPmk8zPmR8z3wOybgJZz4Fiz5kvAMknzHnO15hvOt0Hi94HwhfMz5j//xd2vskMALKohqFTAAAA";

const AMINO_ACIDS: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Minimal RFC 4648 decoder, so vep-core needs no base64 dev-dependency.
fn b64(s: &str) -> Vec<u8> {
    const T: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::with_capacity(s.len() * 3 / 4);
    let mut buf = 0u32;
    let mut bits = 0;
    for c in s.bytes().filter(|c| *c != b'=') {
        let v = T.iter().position(|t| *t == c).expect("base64 alphabet") as u32;
        buf = (buf << 6) | v;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((buf >> bits) as u8);
            buf &= (1 << bits) - 1;
        }
    }
    out
}

fn count_by_label(data: &[u8], residues: usize, analysis: AnalysisType) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for pos in 1..=residues {
        for aa in AMINO_ACIDS {
            if let Some(p) = decode_prediction(data, pos, *aa, analysis) {
                *counts.entry(p.label).or_insert(0) += 1;
            }
        }
    }
    counts
}

fn expected(pairs: &[(&str, usize)]) -> BTreeMap<String, usize> {
    pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
}

#[test]
fn sift_matrix_decodes_every_category_perl_unpacks() {
    let data = gunzip_matrix(&b64(SIFT_B64)).unwrap();
    assert_eq!((data.len() - 3) / 40, 39);
    let counts = count_by_label(&data, 39, AnalysisType::Sift);
    assert_eq!(
        counts,
        expected(&[
            ("tolerated", 78),
            ("deleterious", 606),
            ("tolerated_low_confidence", 2),
            ("deleterious_low_confidence", 55)
        ]),
        "per-category counts under Perl's >> 14 rule"
    );
    assert_eq!(counts.values().sum::<usize>(), 741);
    let p = decode_prediction(&data, 2, b'P', AnalysisType::Sift).unwrap();
    assert_eq!(format_prediction(&p, "b"), "tolerated(0.16)");
    let p = decode_prediction(&data, 1, b'A', AnalysisType::Sift).unwrap();
    assert_eq!(format_prediction(&p, "b"), "deleterious(0)");
    let p = decode_prediction(&data, 37, b'R', AnalysisType::Sift).unwrap();
    assert_eq!(format_prediction(&p, "b"), "tolerated_low_confidence(0.22)");
    let p = decode_prediction(&data, 37, b'A', AnalysisType::Sift).unwrap();
    assert_eq!(
        format_prediction(&p, "b"),
        "deleterious_low_confidence(0.04)"
    );
}

#[test]
fn polyphen_matrix_decodes_every_category_perl_unpacks() {
    let data = gunzip_matrix(&b64(POLYPHEN_B64)).unwrap();
    assert_eq!((data.len() - 3) / 40, 2);
    let counts = count_by_label(&data, 2, AnalysisType::PolyPhenHumVar);
    assert_eq!(
        counts,
        expected(&[
            ("probably_damaging", 32),
            ("possibly_damaging", 5),
            ("benign", 1)
        ]),
        "per-category counts under Perl's >> 14 rule"
    );
    assert_eq!(counts.values().sum::<usize>(), 38);
    let p = decode_prediction(&data, 1, b'A', AnalysisType::PolyPhenHumVar).unwrap();
    assert_eq!(format_prediction(&p, "b"), "probably_damaging(0.998)");
    let p = decode_prediction(&data, 1, b'F', AnalysisType::PolyPhenHumVar).unwrap();
    assert_eq!(format_prediction(&p, "b"), "possibly_damaging(0.796)");
    let p = decode_prediction(&data, 1, b'L', AnalysisType::PolyPhenHumVar).unwrap();
    assert_eq!(format_prediction(&p, "b"), "benign(0.142)");
}
