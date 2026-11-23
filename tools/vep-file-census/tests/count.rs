// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `count` held to `grep -c` and to the row rules of two `awk` passes: each golden under
//! `tests/fixtures/expected/` is the standard output of the pipeline its case names, run
//! inside `tests/fixtures/` under `set -o pipefail`, and the case's status is that
//! pipeline's.

mod common;

use common::{hold, Case};

const CASES: &[Case] = &[
    Case {
        name: "count_records",
        rule: "grep -vc '^#' records.vcf",
        stdin: None,
        args: &["count", "--not-starting", "#", "records.vcf"],
        status: 0,
    },
    Case {
        name: "count_no_final_newline",
        rule: "grep -vc '^#' records_no_newline.vcf",
        stdin: None,
        args: &["count", "--not-starting", "#", "records_no_newline.vcf"],
        status: 0,
    },
    Case {
        name: "count_empty",
        rule: "grep -vc '^#' empty.vcf",
        stdin: None,
        args: &["count", "--not-starting", "#", "empty.vcf"],
        status: 1,
    },
    Case {
        name: "count_missing",
        rule: "grep -vc '^#' missing.vcf",
        stdin: None,
        args: &["count", "--not-starting", "#", "missing.vcf"],
        status: 2,
    },
    Case {
        name: "count_spaced_name",
        rule: "grep -vc '^#' 'spaced name.vcf'",
        stdin: None,
        args: &["count", "--not-starting", "#", "spaced name.vcf"],
        status: 0,
    },
    Case {
        name: "count_directory",
        rule: "grep -vc '^#' tree",
        stdin: None,
        args: &["count", "--not-starting", "#", "tree"],
        status: 2,
    },
    Case {
        name: "count_gzip",
        rule: "gzip -dc records.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &["count", "--gunzip", "--not-starting", "#", "records.vcf.gz"],
        status: 0,
    },
    Case {
        name: "count_zcat",
        rule: "zcat chr.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &["count", "--gunzip", "--not-starting", "#", "chr.vcf.gz"],
        status: 0,
    },
    Case {
        name: "count_block_gzip",
        rule: "gzip -dc chr_blocks.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &[
            "count",
            "--gunzip",
            "--not-starting",
            "#",
            "chr_blocks.vcf.gz",
        ],
        status: 0,
    },
    Case {
        name: "count_two_members",
        rule: "gzip -dc two_members.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &[
            "count",
            "--gunzip",
            "--not-starting",
            "#",
            "two_members.vcf.gz",
        ],
        status: 0,
    },
    Case {
        name: "count_empty_gzip",
        rule: "gzip -dc empty.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &["count", "--gunzip", "--not-starting", "#", "empty.vcf.gz"],
        status: 1,
    },
    Case {
        name: "count_missing_gzip",
        rule: "gzip -dc missing.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &["count", "--gunzip", "--not-starting", "#", "missing.vcf.gz"],
        status: 1,
    },
    Case {
        name: "count_garbage_after_member",
        rule: "gzip -dc garbage.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &["count", "--gunzip", "--not-starting", "#", "garbage.vcf.gz"],
        status: 2,
    },
    Case {
        name: "count_zero_padded",
        rule: "gzip -dc zero_padded.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &[
            "count",
            "--gunzip",
            "--not-starting",
            "#",
            "zero_padded.vcf.gz",
        ],
        status: 0,
    },
    Case {
        name: "count_truncated",
        rule: "gzip -dc truncated.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &[
            "count",
            "--gunzip",
            "--not-starting",
            "#",
            "truncated.vcf.gz",
        ],
        status: 1,
    },
    Case {
        name: "count_truncated_blocks",
        rule: "gzip -dc truncated_blocks.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &[
            "count",
            "--gunzip",
            "--not-starting",
            "#",
            "truncated_blocks.vcf.gz",
        ],
        status: 1,
    },
    Case {
        name: "count_not_gzip",
        rule: "gzip -dc not_gzip.vcf.gz | grep -vc '^#'",
        stdin: None,
        args: &[
            "count",
            "--gunzip",
            "--not-starting",
            "#",
            "not_gzip.vcf.gz",
        ],
        status: 1,
    },
    Case {
        name: "count_header_lines",
        rule: "grep -c '^#' output.txt",
        stdin: None,
        args: &["count", "--starting", "#", "output.txt"],
        status: 0,
    },
    Case {
        name: "count_no_header_lines",
        rule: "grep -c '^#' output_no_header.txt",
        stdin: None,
        args: &["count", "--starting", "#", "output_no_header.txt"],
        status: 1,
    },
    Case {
        name: "count_body_lines",
        rule: "grep -vc '^##' output_no_newline.txt",
        stdin: None,
        args: &["count", "--not-starting", "##", "output_no_newline.txt"],
        status: 0,
    },
    Case {
        name: "count_containing",
        rule: r#"grep -c '_consequences":' output.jsonl"#,
        stdin: None,
        args: &[
            "count",
            "--containing",
            r#"_consequences":"#,
            "output.jsonl",
        ],
        status: 0,
    },
    Case {
        name: "count_containing_none",
        rule: "grep -c 'no such text' output.jsonl",
        stdin: None,
        args: &["count", "--containing", "no such text", "output.jsonl"],
        status: 1,
    },
    Case {
        name: "count_zero_bytes",
        rule: "grep -vc '^#' zero_bytes.bin",
        stdin: None,
        args: &["count", "--not-starting", "#", "zero_bytes.bin"],
        status: 0,
    },
    Case {
        name: "count_zero_bytes_starting",
        rule: "grep -c '^#' zero_bytes.bin",
        stdin: None,
        args: &["count", "--starting", "#", "zero_bytes.bin"],
        status: 0,
    },
    Case {
        name: "count_zero_bytes_containing",
        rule: "grep -c 'y' zero_bytes.bin",
        stdin: None,
        args: &["count", "--containing", "y", "zero_bytes.bin"],
        status: 0,
    },
    Case {
        name: "count_info_field",
        rule: r#"FIELD=CSQ awk -F"\t" '!/^#/ && index(";" $8, ";" ENVIRON["FIELD"] "=") {c++} END {print c + 0}' output.vcf"#,
        stdin: None,
        args: &["count", "--info-field", "CSQ", "output.vcf"],
        status: 0,
    },
    Case {
        name: "count_info_field_other",
        rule: r#"FIELD=ANN awk -F"\t" '!/^#/ && index(";" $8, ";" ENVIRON["FIELD"] "=") {c++} END {print c + 0}' output.vcf"#,
        stdin: None,
        args: &["count", "--info-field", "ANN", "output.vcf"],
        status: 0,
    },
    Case {
        name: "count_info_field_none",
        rule: r#"FIELD=NONE awk -F"\t" '!/^#/ && index(";" $8, ";" ENVIRON["FIELD"] "=") {c++} END {print c + 0}' output.vcf"#,
        stdin: None,
        args: &["count", "--info-field", "NONE", "output.vcf"],
        status: 0,
    },
    Case {
        name: "count_info_field_records",
        rule: r#"FIELD=CSQ awk -F"\t" '!/^#/ && index(";" $8, ";" ENVIRON["FIELD"] "=") {c++} END {print c + 0}' records_no_newline.vcf"#,
        stdin: None,
        args: &["count", "--info-field", "CSQ", "records_no_newline.vcf"],
        status: 0,
    },
    Case {
        name: "count_extra_key",
        rule: r#"awk -F"\t" -v k=CLIN_SIG '!/^#/ && index(";" $NF ";", ";" k "=") {c++} END {print c + 0}' output.txt"#,
        stdin: None,
        args: &["count", "--extra-key", "CLIN_SIG", "output.txt"],
        status: 0,
    },
    Case {
        name: "count_extra_key_no_newline",
        rule: r#"awk -F"\t" -v k=HGVSc '!/^#/ && index(";" $NF ";", ";" k "=") {c++} END {print c + 0}' output_no_newline.txt"#,
        stdin: None,
        args: &["count", "--extra-key", "HGVSc", "output_no_newline.txt"],
        status: 0,
    },
    Case {
        name: "count_several_rules",
        rule: r#"grep -vc '^#' output.txt; awk -F"\t" -v k=IMPACT '!/^#/ && index(";" $NF ";", ";" k "=") {c++} END {print c + 0}' output.txt; awk -F"\t" -v k=HGVSp '!/^#/ && index(";" $NF ";", ";" k "=") {c++} END {print c + 0}' output.txt; awk -F"\t" -v k=AF '!/^#/ && index(";" $NF ";", ";" k "=") {c++} END {print c + 0}' output.txt"#,
        stdin: None,
        args: &[
            "count",
            "--not-starting",
            "#",
            "--extra-key",
            "IMPACT",
            "--extra-key",
            "HGVSp",
            "--extra-key",
            "AF",
            "output.txt",
        ],
        status: 0,
    },
];

#[test]
fn a_count_is_what_the_pipeline_prints_and_the_status_is_the_pipelines() {
    for case in CASES {
        hold(case);
    }
}
