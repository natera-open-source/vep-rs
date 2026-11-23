// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `strip-chr`, `contigs` and `slice` held to the `awk` passes they stand for: each golden
//! under `tests/fixtures/expected/` is the standard output of the pipeline its case names,
//! run inside `tests/fixtures/` under `set -o pipefail`, and the case's status is that
//! pipeline's. A slice's `.table` golden is the file its pipeline appended to.

mod common;

use common::{hold, run, Case, Scratch, THREADS};

const STRIP: &[Case] = &[
    Case {
        name: "strip_plain",
        rule: r#"zcat -f chr.vcf | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "chr.vcf"],
        status: 0,
    },
    Case {
        name: "strip_gzip",
        rule: r#"zcat -f chr.vcf.gz | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "chr.vcf.gz"],
        status: 0,
    },
    Case {
        name: "strip_block_gzip",
        rule: r#"zcat -f chr_blocks.vcf.gz | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "chr_blocks.vcf.gz"],
        status: 0,
    },
    Case {
        name: "strip_no_prefix_no_final_newline",
        rule: r#"zcat -f records_no_newline.vcf | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "records_no_newline.vcf"],
        status: 0,
    },
    Case {
        name: "strip_empty",
        rule: r#"zcat -f empty.vcf | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "empty.vcf"],
        status: 0,
    },
    Case {
        name: "strip_missing",
        rule: r#"zcat -f missing.vcf.gz | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "missing.vcf.gz"],
        status: 1,
    },
    Case {
        name: "strip_garbage_after_member",
        rule: r#"zcat -f garbage.vcf.gz | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "garbage.vcf.gz"],
        status: 0,
    },
    Case {
        name: "strip_spaced_name",
        rule: r#"zcat -f 'spaced name.vcf' | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next} /^#/ {print; next} {sub(/^chr/,"",$1); print}'"#,
        stdin: None,
        args: &["strip-chr", "spaced name.vcf"],
        status: 0,
    },
];

#[test]
fn contig_names_lose_their_prefix_as_the_awk_program_writes_them() {
    for case in STRIP {
        hold(case);
    }
}

const CONTIGS: &[Case] = &[
    Case {
        name: "contigs_plain",
        rule: r#"cat chr.vcf | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "chr.vcf"],
        status: 0,
    },
    Case {
        name: "contigs_block_gzip",
        rule: r#"gzip -dc chr_blocks.vcf.gz | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "--gunzip", "chr_blocks.vcf.gz"],
        status: 0,
    },
    Case {
        name: "contigs_no_final_newline",
        rule: r#"cat records_no_newline.vcf | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "records_no_newline.vcf"],
        status: 0,
    },
    Case {
        name: "contigs_empty",
        rule: r#"cat empty.vcf | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "empty.vcf"],
        status: 0,
    },
    Case {
        name: "contigs_missing",
        rule: r#"cat missing.vcf | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "missing.vcf"],
        status: 1,
    },
    Case {
        name: "contigs_gzip",
        rule: r#"gzip -dc chr.vcf.gz | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "--gunzip", "chr.vcf.gz"],
        status: 0,
    },
    Case {
        name: "contigs_spaced_name",
        rule: r#"cat 'spaced name.vcf' | awk -F"\t" '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'"#,
        stdin: None,
        args: &["contigs", "spaced name.vcf"],
        status: 0,
    },
];

#[test]
fn a_tally_is_what_the_awk_pass_prints() {
    for case in CONTIGS {
        hold(case);
    }
}

const SLICES: &[Case] = &[
    Case {
        name: "slice_plain",
        rule: r#"cat chr.vcf | awk -F"\t" -v n=12 -v k=$(grep -c . chr.tally.tsv) -v label=inputs/chr.vcf -v tsv=expected/slice_plain.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' chr.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--records",
            "12",
            "--tally",
            "chr.tally.tsv",
            "--label",
            "inputs/chr.vcf",
            "chr.vcf",
        ],
        status: 0,
    },
    Case {
        name: "slice_block_gzip",
        rule: r#"gzip -dc chr_blocks.vcf.gz | awk -F"\t" -v n=50 -v k=$(grep -c . chr.tally.tsv) -v label=chr -v tsv=expected/slice_block_gzip.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' chr.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--gunzip",
            "--records",
            "50",
            "--tally",
            "chr.tally.tsv",
            "--label",
            "chr",
            "chr_blocks.vcf.gz",
        ],
        status: 0,
    },
    Case {
        name: "slice_more_than_held",
        rule: r#"cat chr.vcf | awk -F"\t" -v n=1000 -v k=$(grep -c . chr.tally.tsv) -v label=all -v tsv=expected/slice_more_than_held.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' chr.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--records",
            "1000",
            "--tally",
            "chr.tally.tsv",
            "--label",
            "all",
            "chr.vcf",
        ],
        status: 0,
    },
    Case {
        name: "slice_one_a_contig",
        rule: r#"cat records_no_newline.vcf | awk -F"\t" -v n=5 -v k=$(grep -c . records.tally.tsv) -v label=five -v tsv=expected/slice_one_a_contig.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' records.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--records",
            "5",
            "--tally",
            "records.tally.tsv",
            "--label",
            "five",
            "records_no_newline.vcf",
        ],
        status: 0,
    },
    Case {
        name: "slice_gzip",
        rule: r#"gzip -dc chr.vcf.gz | awk -F"\t" -v n=9 -v k=$(grep -c . chr.tally.tsv) -v label=gz -v tsv=expected/slice_gzip.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' chr.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--gunzip",
            "--records",
            "9",
            "--tally",
            "chr.tally.tsv",
            "--label",
            "gz",
            "chr.vcf.gz",
        ],
        status: 0,
    },
    Case {
        name: "slice_spaced_name",
        rule: r#"cat 'spaced name.vcf' | awk -F"\t" -v n=3 -v k=$(grep -c . records.tally.tsv) -v label=spaced -v tsv=expected/slice_spaced_name.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' records.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--records",
            "3",
            "--tally",
            "records.tally.tsv",
            "--label",
            "spaced",
            "spaced name.vcf",
        ],
        status: 0,
    },
    Case {
        name: "slice_empty",
        rule: r#"cat empty.vcf | awk -F"\t" -v n=3 -v k=$(grep -c . records.tally.tsv) -v label=none -v tsv=expected/slice_empty.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' records.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--records",
            "3",
            "--tally",
            "records.tally.tsv",
            "--label",
            "none",
            "empty.vcf",
        ],
        status: 0,
    },
    Case {
        name: "slice_missing",
        rule: r#"cat missing.vcf | awk -F"\t" -v n=3 -v k=$(grep -c . records.tally.tsv) -v label=none -v tsv=expected/slice_missing.table '
	FNR == NR {idx[$1] = FNR; c[$1] = $2; next}
	/^#/ {print; next}
	{
		if (!($1 in q)) {
			q[$1] = int(n / k) + (idx[$1] <= n % k)
			if (q[$1] > c[$1]) q[$1] = c[$1]
			order[++m] = $1
		}
		i = ++seen[$1]
		if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print
	}
	END {for (j = 1; j <= m; j++) print label "\t" order[j] "\t" c[order[j]] "\t" q[order[j]] >>tsv}' records.tally.tsv -"#,
        stdin: None,
        args: &[
            "slice",
            "--records",
            "3",
            "--tally",
            "records.tally.tsv",
            "--label",
            "none",
            "missing.vcf",
        ],
        status: 1,
    },
];

#[test]
fn a_slice_and_its_table_are_what_the_awk_pass_writes() {
    for case in SLICES {
        hold(case);
        // A slice without a record has no table golden: its pass appends nothing.
        let golden = std::fs::read(
            common::fixtures()
                .join("expected")
                .join(format!("{}.table", case.name)),
        )
        .unwrap_or_default();
        for (at, threads) in THREADS.into_iter().enumerate() {
            let scratch = Scratch::new(&format!("{}-{at}", case.name));
            let table = scratch.path("inputs.tsv");
            // The pass appends: what the table held before stays.
            std::fs::write(&table, "input\tcontig\trecords\ttaken\n").unwrap();
            let mut args: Vec<String> = case.args.iter().map(|arg| (*arg).to_string()).collect();
            let file = args.pop().unwrap();
            args.extend([
                "--table".to_string(),
                table.to_str().unwrap().to_string(),
                file,
            ]);
            let output = run(&args, None, threads);
            assert_eq!(output.status.code(), Some(case.status), "{}", case.name);
            let mut whole = b"input\tcontig\trecords\ttaken\n".to_vec();
            whole.extend_from_slice(&golden);
            assert_eq!(
                std::fs::read(&table).unwrap(),
                whole,
                "{}: the table",
                case.name
            );
        }
    }
}

#[test]
fn a_slice_without_a_record_leaves_no_table() {
    let scratch = Scratch::new("no-table");
    let table = scratch.path("inputs.tsv");
    let output = run(
        &[
            "slice",
            "--records",
            "3",
            "--tally",
            "chr.tally.tsv",
            "--table",
            table.to_str().unwrap(),
            "empty.vcf",
        ],
        None,
        Some("2"),
    );
    assert_eq!(output.status.code(), Some(0));
    assert!(output.stdout.is_empty());
    assert!(!table.exists());
}
