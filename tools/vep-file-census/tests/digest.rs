// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `digest` held to `sha256sum`, `md5sum` and `xxh128sum`: each golden under
//! `tests/fixtures/expected/` is the standard output of the pipeline its case names, run
//! inside `tests/fixtures/` under `set -o pipefail` with the case's fixture as its
//! standard input, and the case's status is that pipeline's.

mod common;

use common::{hold, Case};

const CASES: &[Case] = &[
    Case {
        name: "digest_sha256",
        rule: "sha256sum records.vcf",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "records.vcf"],
        status: 0,
    },
    Case {
        name: "digest_sha256_several",
        rule: "sha256sum empty.vcf records_no_newline.vcf 'spaced name.vcf' missing.vcf records.vcf.gz",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "empty.vcf", "records_no_newline.vcf", "spaced name.vcf", "missing.vcf", "records.vcf.gz"],
        status: 1,
    },
    Case {
        name: "digest_md5",
        rule: "md5sum records.vcf 'spaced name.vcf' empty.vcf",
        stdin: None,
        args: &["digest", "--algorithm", "md5", "records.vcf", "spaced name.vcf", "empty.vcf"],
        status: 0,
    },
    Case {
        name: "digest_xxh128",
        rule: "xxh128sum records.vcf 'spaced name.vcf' empty.vcf",
        stdin: None,
        args: &["digest", "--algorithm", "xxh128", "records.vcf", "spaced name.vcf", "empty.vcf"],
        status: 0,
    },
    Case {
        name: "digest_two_algorithms",
        rule: "sha256sum records.vcf; md5sum records.vcf",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--algorithm", "md5", "records.vcf"],
        status: 0,
    },
    Case {
        name: "digest_standard_input",
        rule: "sha256sum -",
        stdin: Some("records.vcf"),
        args: &["digest", "--algorithm", "sha256", "-"],
        status: 0,
    },
    Case {
        name: "digest_directory",
        rule: "sha256sum tree",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "tree"],
        status: 1,
    },
    Case {
        name: "digest_names_nul",
        rule: "xargs -0 sha256sum",
        stdin: Some("names_tree.nul"),
        args: &["digest", "--algorithm", "sha256", "--names", "nul"],
        status: 0,
    },
    Case {
        name: "digest_names_nul_none",
        rule: "xargs -0 sha256sum",
        stdin: Some("names_none.txt"),
        args: &["digest", "--algorithm", "sha256", "--names", "nul"],
        status: 0,
    },
    Case {
        name: "digest_names_words",
        rule: "xargs sha256sum",
        stdin: Some("names_words.txt"),
        args: &["digest", "--algorithm", "sha256", "--names", "words"],
        status: 0,
    },
    Case {
        name: "digest_names_words_spaced",
        rule: "xargs sha256sum",
        stdin: Some("names_tree.txt"),
        args: &["digest", "--algorithm", "sha256", "--names", "words"],
        status: 123,
    },
    Case {
        name: "digest_names_words_none",
        rule: "xargs md5sum",
        stdin: Some("names_none.txt"),
        args: &["digest", "--algorithm", "md5", "--names", "words"],
        status: 0,
    },
    Case {
        name: "digest_names_read",
        rule: r#"while read -r f; do md5sum "$f"; done"#,
        stdin: Some("names_json.txt"),
        args: &["digest", "--algorithm", "md5", "--names", "read"],
        status: 0,
    },
    Case {
        name: "digest_names_read_blanks",
        rule: r#"while read -r f; do md5sum "$f"; done"#,
        stdin: Some("names_blanks.txt"),
        args: &["digest", "--algorithm", "md5", "--names", "read"],
        status: 0,
    },
    Case {
        name: "digest_names_read_missing_last",
        rule: r#"while read -r f; do md5sum "$f"; done"#,
        stdin: Some("names_missing_last.txt"),
        args: &["digest", "--algorithm", "md5", "--names", "read"],
        status: 1,
    },
    Case {
        name: "digest_names_lines",
        rule: r#"while IFS= read -r f; do sha256sum "$f"; done"#,
        stdin: Some("names_tree.txt"),
        args: &["digest", "--algorithm", "sha256", "--names", "lines"],
        status: 0,
    },
    Case {
        name: "digest_names_lines_blanks",
        rule: r#"while IFS= read -r f; do sha256sum "$f"; done"#,
        stdin: Some("names_blanks.txt"),
        args: &["digest", "--algorithm", "sha256", "--names", "lines"],
        status: 1,
    },
    Case {
        name: "digest_names_lines_none",
        rule: r#"while IFS= read -r f; do sha256sum "$f"; done"#,
        stdin: Some("names_none.txt"),
        args: &["digest", "--algorithm", "sha256", "--names", "lines"],
        status: 0,
    },
    Case {
        name: "digest_gunzip",
        rule: "gzip -dc records.vcf.gz | sha256sum",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--gunzip", "records.vcf.gz"],
        status: 0,
    },
    Case {
        name: "digest_gunzip_blocks",
        rule: "gzip -dc chr_blocks.vcf.gz | md5sum",
        stdin: None,
        args: &["digest", "--algorithm", "md5", "--gunzip", "chr_blocks.vcf.gz"],
        status: 0,
    },
    Case {
        name: "digest_gunzip_garbage",
        rule: "gzip -dc garbage.vcf.gz | sha256sum",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--gunzip", "garbage.vcf.gz"],
        status: 2,
    },
    Case {
        name: "digest_gunzip_missing",
        rule: "gzip -dc missing.vcf.gz | sha256sum",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--gunzip", "missing.vcf.gz"],
        status: 1,
    },
    Case {
        name: "digest_without_stamp",
        rule: "grep -a -v '^## Output produced at' output.txt | sha256sum",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--without-lines-starting", "## Output produced at", "output.txt"],
        status: 0,
    },
    Case {
        name: "digest_without_header",
        rule: "grep -a -v '^##' output_no_newline.txt | xxh128sum",
        stdin: None,
        args: &["digest", "--algorithm", "xxh128", "--without-lines-starting", "##", "output_no_newline.txt"],
        status: 0,
    },
    Case {
        name: "digest_without_header_zero_bytes",
        rule: "grep -a -v '^#' zero_bytes.bin | md5sum",
        stdin: None,
        args: &["digest", "--algorithm", "md5", "--without-lines-starting", "#", "zero_bytes.bin"],
        status: 0,
    },
    Case {
        name: "digest_without_header_empty",
        rule: "grep -a -v '^##' empty.vcf | sha256sum",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--without-lines-starting", "##", "empty.vcf"],
        status: 1,
    },
    Case {
        name: "digest_gunzip_without_header",
        rule: "gzip -dc chr_blocks.vcf.gz | grep -a -v '^##' | sha256sum",
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--gunzip", "--without-lines-starting", "##", "chr_blocks.vcf.gz"],
        status: 0,
    },
    Case {
        name: "digest_after_header",
        rule: r#"tail -c +$(($(sed -n "/^##/!q;p" output.txt | wc -c) + 1)) output.txt | sha256sum"#,
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--after-leading-lines-starting", "##", "output.txt"],
        status: 0,
    },
    Case {
        name: "digest_after_header_none",
        rule: r#"tail -c +$(($(sed -n "/^##/!q;p" output_no_header.txt | wc -c) + 1)) output_no_header.txt | sha256sum"#,
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--after-leading-lines-starting", "##", "output_no_header.txt"],
        status: 0,
    },
    Case {
        name: "digest_after_header_all",
        rule: r#"tail -c +$(($(sed -n "/^##/!q;p" empty.vcf | wc -c) + 1)) empty.vcf | xxh128sum"#,
        stdin: None,
        args: &["digest", "--algorithm", "xxh128", "--after-leading-lines-starting", "##", "empty.vcf"],
        status: 0,
    },
    Case {
        name: "digest_chunks_xxh128",
        rule: r#"tail -c +$(($(sed -n "/^##/!q;p" output.txt | wc -c) + 1)) output.txt | split -b 100 --filter=xxh128sum | awk '{print $1}' | xxh128sum"#,
        stdin: None,
        args: &["digest", "--algorithm", "xxh128", "--after-leading-lines-starting", "##", "--chunk-bytes", "100", "output.txt"],
        status: 0,
    },
    Case {
        name: "digest_chunks_sha256",
        rule: r#"tail -c +$(($(sed -n "/^##/!q;p" output_no_newline.txt | wc -c) + 1)) output_no_newline.txt | split -b 64 --filter=sha256sum | awk '{print $1}' | sha256sum"#,
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--after-leading-lines-starting", "##", "--chunk-bytes", "64", "output_no_newline.txt"],
        status: 0,
    },
    Case {
        name: "digest_chunks_one",
        rule: r#"tail -c +$(($(sed -n "/^##/!q;p" output.txt | wc -c) + 1)) output.txt | split -b 1000000 --filter=sha256sum | awk '{print $1}' | sha256sum"#,
        stdin: None,
        args: &["digest", "--algorithm", "sha256", "--after-leading-lines-starting", "##", "--chunk-bytes", "1000000", "output.txt"],
        status: 0,
    },
    Case {
        name: "digest_chunks_whole_file",
        rule: "split -b 300 --filter=md5sum records.vcf | awk '{print $1}' | md5sum",
        stdin: None,
        args: &["digest", "--algorithm", "md5", "--chunk-bytes", "300", "records.vcf"],
        status: 0,
    },
];

#[test]
fn digest_lines_are_what_the_pipeline_prints_and_the_status_is_the_pipelines() {
    for case in CASES {
        hold(case);
    }
}
