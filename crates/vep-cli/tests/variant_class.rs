// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `Allele` and `VARIANT_CLASS` columns of a structural record whose type
//! the Sequence Ontology term table lacks, as the `vep` binary prints them.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The GRCh37-output-flags corpus, whose cache covers 21:19744457-19744524.
fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-output-flags")
}

/// The rows of a default-format run on `input` against the corpus's cache with
/// `extra` flags, each row's identifier, allele and Extra column.
fn rows(dir: &Path, input: &Path, extra: &[&str]) -> Vec<(String, String, String)> {
    let corpus = common::corpora()
        .into_iter()
        .find(|corpus| corpus.dir.canonicalize().ok() == corpus_dir().canonicalize().ok())
        .expect("the GRCh37-output-flags corpus is committed");
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
            "GRCh37",
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
                fields[2].to_string(),
                fields[13].to_string(),
            )
        })
        .collect()
}

/// A symbolic ALT with no Sequence Ontology term keeps Perl's type string as
/// its `class_SO_term` (Parser/VCF.pm:469-481: the ALTs joined by `/`, replaced
/// by INFO/SVTYPE when that names the type, kept whole when `get_SO_term` has no
/// term for it), printed in `Allele` and `VARIANT_CLASS` (OutputFactory.pm:967):
/// Ensembl VEP 116.2 under `--dont_skip --variant_class` prints `<CPX>` for
/// `21 19744457 cpx_bare N <CPX> . . END=19744524`, `CPX` for the same record
/// with `SVTYPE=CPX`, and `FOO` with `SVTYPE=FOO`; `<CTX>` reads like `<CPX>`.
#[test]
fn typeless_symbolic_alt_prints_perls_type_string() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("cpx.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t19744457\tcpx_bare\tN\t<CPX>\t.\t.\tEND=19744524\n\
         21\t19744457\tcpx_svtype\tN\t<CPX>\t.\t.\tSVTYPE=CPX;END=19744524\n\
         21\t19744457\tcpx_foo\tN\t<CPX>\t.\t.\tSVTYPE=FOO;END=19744524\n\
         21\t19744457\tctx_bare\tN\t<CTX>\t.\t.\tEND=19744524\n",
    )
    .unwrap();
    let rows = rows(dir.path(), &input, &["--dont_skip", "--variant_class"]);
    // One triple per record: vep-rs writes a kept record's row per overlapping
    // transcript, the coverage divergence docs/divergences/defect-c.md states.
    let mut classes: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|(id, allele, extra)| {
            let class = extra
                .split(';')
                .find_map(|kv| kv.strip_prefix("VARIANT_CLASS="))
                .unwrap_or_else(|| panic!("no VARIANT_CLASS on {id}: {extra}"));
            (id.as_str(), allele.as_str(), class)
        })
        .collect();
    classes.dedup();
    assert_eq!(
        classes,
        [
            ("cpx_bare", "<CPX>", "<CPX>"),
            ("cpx_svtype", "CPX", "CPX"),
            ("cpx_foo", "FOO", "FOO"),
            ("ctx_bare", "<CTX>", "<CTX>"),
        ]
    );
}
