// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in plugins through the `vep` binary: the behaviours a golden corpus
//! cannot carry because they vary the input or the arguments against one
//! committed fixture; the bundled score files LoFtool and pLI read from
//! `--dir_plugins` when given no file, compared with the files named outright;
//! and the parameter grammars, fields and refusals of the plugins read on the
//! corpus's own synthetic data file; REVEL's and AlphaMissense's grammars and
//! options over their real-data corpora; the shape SpliceAI's fields take in
//! each output format; and LoFTEE's calls on the inputs its corpus does not
//! carry: a one-base acceptor deletion, an exon-spanning deletion and the
//! spellings of `check_complete_cds`.

mod common;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use common::{corpora, read_expected};

fn vep() -> Command {
    common::isolated_vep_anywhere()
}

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// The `gnomAD_*` fields of a default-format output, keyed by the record's
/// `Location` without a `chr` prefix and then by field name; a record no row
/// annotates has an empty map.
fn gnomad_fields(output: &str) -> BTreeMap<String, BTreeMap<String, String>> {
    let mut fields: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for row in output.lines().filter(|line| !line.starts_with('#')) {
        let columns: Vec<&str> = row.split('\t').collect();
        let location = columns[1].trim_start_matches("chr").to_string();
        let record = fields.entry(location).or_default();
        for pair in columns[13].split(';') {
            if let Some((key, value)) = pair.split_once('=') {
                if key.starts_with("gnomAD_") {
                    record.insert(key.to_string(), value.to_string());
                }
            }
        }
    }
    fields
}

/// `variants.vcf` of `dir` with every record's chromosome written `chr21`.
fn chr_prefixed_input(dir: &Path, out: &Path) {
    let text = std::fs::read_to_string(dir.join("variants.vcf")).unwrap();
    let prefixed: String = text
        .lines()
        .map(|line| {
            if line.starts_with('#') {
                line.to_string()
            } else {
                format!("chr{line}")
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(out, prefixed + "\n").unwrap();
}

/// gnomADc reads an input chromosome spelled `chr21` against a coverage file
/// naming its contigs `chr21` (the GRCh38 corpus) and against one naming them
/// `21` (the GRCh37 corpus), as Ensembl VEP's plugin strips the prefix and its
/// tabix base class tries the file's own spelling (`gnomADc.pm:158`,
/// `BaseVepTabixPlugin.pm:693-708`): every record annotated in the corpus
/// reference carries the same coverage means, and the record the file does not
/// cover carries none.
#[test]
fn gnomadc_chr_prefixed_input_matches() {
    let dir = tempfile::tempdir().unwrap();
    for (name, assembly) in [
        ("GRCh38-plugin-gnomadc", "GRCh38"),
        ("GRCh37-plugin-gnomadc", "GRCh37"),
    ] {
        let corpus = corpus(name);
        let input = dir.path().join(format!("{name}.vcf"));
        chr_prefixed_input(&corpus, &input);
        let out = dir.path().join(format!("{name}.txt"));
        let coverage = corpus.join("plugin_data/gnomad.coverage.tsv.gz");
        let run = vep()
            .args([
                "-i",
                input.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--offline",
                "--json_cache",
                corpus.join("json_cache").to_str().unwrap(),
                "--assembly",
                assembly,
                "--force_overwrite",
                "--no_stats",
                "--dir_plugins",
                corpus.join("plugin_data").to_str().unwrap(),
                "--plugin",
                &format!("gnomADc,{}", coverage.display()),
            ])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            run.status.code(),
            Some(0),
            "{name}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let output = std::fs::read_to_string(&out).unwrap();
        assert!(
            output.contains("\n## gnomAD_mean_cov : mean coverage\n"),
            "{name}: no gnomAD_mean_cov header line"
        );
        let expected = gnomad_fields(&read_expected(&corpus, "default.txt"));
        let actual = gnomad_fields(&output);
        assert_eq!(actual.len(), expected.len(), "{name}: record count");
        for (record, fields) in &expected {
            assert_eq!(
                actual.get(record),
                Some(fields),
                "{name}: the coverage fields of {record} under a chr-prefixed input"
            );
        }
        assert!(
            expected.values().any(|fields| !fields.is_empty()),
            "{name}: no record carries coverage fields"
        );
        assert!(
            output
                .lines()
                .filter(|line| !line.starts_with('#'))
                .all(|row| row
                    .split('\t')
                    .nth(1)
                    .is_some_and(|loc| loc.starts_with("chr21:"))),
            "{name}: every row keeps the input's chr21 spelling"
        );
    }
}

/// Runs `vep` on a CADD corpus's `variants.vcf` in the default format without
/// `--quiet`, so a warning reaches standard error, with `--dir_plugins` on the
/// corpus's `plugin_data` as the corpus manifests pass it.
fn vep_on_cadd_corpus(corpus: &Path, out: &Path, plugin: &str) -> std::process::Output {
    vep()
        .current_dir(corpus)
        .env("VEP_DIR_PLUGINS", corpus)
        .args([
            "-i",
            corpus.join("variants.vcf").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--offline",
            "--json_cache",
            corpus.join("json_cache").to_str().unwrap(),
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--dir_plugins",
            corpus.join("plugin_data").to_str().unwrap(),
            "--plugin",
            plugin,
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

/// The `(CADD_PHRED, CADD_RAW)` fields of every row of a default-format
/// output, keyed by the row's `Uploaded_variation`.
fn cadd_fields(output: &str) -> BTreeMap<String, (Option<String>, Option<String>)> {
    output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            let field = |name: &str| {
                columns[13].split(';').find_map(|pair| {
                    pair.strip_prefix(name)?
                        .strip_prefix('=')
                        .map(str::to_owned)
                })
            };
            (
                columns[0].to_string(),
                (field("CADD_PHRED"), field("CADD_RAW")),
            )
        })
        .collect()
}

/// The `(PHRED, raw)` scores of a bgzipped CADD file's rows, keyed by
/// `chrom:pos:ref:alt` for an SNV or indel file and `chrom:start:type` for a
/// CADD-SV file, the columns resolved from the `#Chrom` header as the plugin
/// resolves them (CADD.pm:265-268).
fn cadd_file_scores(path: &Path) -> BTreeMap<String, (Option<String>, Option<String>)> {
    let mut reader = noodles::bgzf::io::Reader::new(std::fs::File::open(path).unwrap());
    let mut text = String::new();
    std::io::Read::read_to_string(&mut reader, &mut text).unwrap();
    let mut header: Vec<&str> = Vec::new();
    let mut scores = BTreeMap::new();
    for line in text.lines() {
        if let Some(columns) = line
            .strip_prefix('#')
            .filter(|_| line.starts_with("#Chrom"))
        {
            header = columns.split('\t').collect();
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let cell = |name: &str| {
            header
                .iter()
                .position(|column| *column == name)
                .map(|index| cells[index].to_string())
        };
        let key = match (cell("Pos"), cell("Start")) {
            (Some(pos), _) => format!(
                "{}:{pos}:{}:{}",
                cells[0],
                cell("Ref").unwrap(),
                cell("Alt").unwrap()
            ),
            (None, Some(start)) => format!("{}:{start}:{}", cells[0], cell("Type").unwrap()),
            (None, None) => panic!("{}: no position column", path.display()),
        };
        let phred = cell("PHRED").or_else(|| cell("CADD-SV_PHRED-score"));
        let raw = cell("RawScore").or_else(|| cell("CADD-SV_Raw-score"));
        scores.insert(key, (phred, raw));
    }
    scores
}

/// CADD reads the SNV and indel files under `snv=` and `indels=`, matching an
/// allele after the common bases are trimmed from both ends (CADD.pm:239-251),
/// so the substitution `TCC/TTC` takes the `C/T` row inside it; a structural
/// variant is skipped while a sequence file is loaded (`:218`) and annotated
/// from a CADD-SV file under `sv=` by its class word and the base before its
/// start (`:220-224`); a trailing `1` of the positional list is `force_annotate`
/// and no file (`:119-120`); a file given under any other key is opened and
/// marks no sequence file (`:124-129`);
/// `force_annotate=1` prints Ensembl VEP's warning (`:137-139`); and the header
/// carries the non-commercial clause (`:100, 174`). The corpora's score files
/// are made up; `whole_genome_SNVs_synthetic.tsv.gz` and
/// `1000G_phase3_SVs_synthetic.tsv.gz` carry the file names of the plugin's
/// usage example (`:31-32`) with `synthetic` added before the extension, and
/// `gnomad.genomes.r4.0.indel.synthetic.tsv.gz` is the corpus's own name for its
/// indel file.
#[test]
fn cadd_snv_indels_and_sv_fixtures() {
    let tmp = tempfile::tempdir().unwrap();
    let cadd = corpus("GRCh38-plugin-cadd-fixture");
    let snv = cadd.join("plugin_data/whole_genome_SNVs_synthetic.tsv.gz");
    let indels = cadd.join("plugin_data/gnomad.genomes.r4.0.indel.synthetic.tsv.gz");
    let plugin = format!("CADD,snv={},indels={}", snv.display(), indels.display());
    let out = tmp.path().join("snv_indels.txt");
    let run = vep_on_cadd_corpus(&cadd, &out, &plugin);
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(0), "{stderr}");
    assert!(
        !stderr.contains("WARNING: Using snv"),
        "no force_annotate warning without the flag: {stderr}"
    );
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        output.contains(
            "\n## CADD_PHRED : PHRED-like scaled CADD score. CADD is only available here for \
             non-commercial use. See CADD website for more information.\n"
        ),
        "no CADD_PHRED header line with the clause"
    );
    let snv_scores = cadd_file_scores(&snv);
    let indel_scores = cadd_file_scores(&indels);
    let found = cadd_fields(&output);
    assert_eq!(found["intergenic_snv"], snv_scores["21:19268555:C:T"]);
    assert_eq!(found["missense_fwd"], snv_scores["21:30396091:C:T"]);
    assert_eq!(
        found["shared_bases_substitution"], snv_scores["21:30396091:C:T"],
        "TCC/TTC at 30396090 is the C/T change at 30396091"
    );
    assert_eq!(found["deletion_rev"], indel_scores["21:30881493:GA:G"]);
    assert_eq!(found["insertion_rev"], indel_scores["21:30881500:T:TA"]);
    assert_eq!(
        found["coding_mnv"],
        (None, None),
        "no row for a two-base change"
    );
    assert_eq!(
        found["unscored_snv"],
        (None, None),
        "the site the file leaves out"
    );
    assert_eq!(
        found["deletion_sv"],
        (None, None),
        "a structural variant beside a sequence file"
    );

    let forced = tmp.path().join("forced.txt");
    let run = vep_on_cadd_corpus(&cadd, &forced, &format!("{plugin},force_annotate=1"));
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(0), "{stderr}");
    assert!(
        stderr.contains(
            "WARNING: Using snv and/or indels CADD annotation file with structural variant can \
             increase run time exponentially."
        ),
        "{stderr}"
    );
    assert_eq!(
        cadd_fields(&std::fs::read_to_string(&forced).unwrap()),
        found,
        "the sequence files score no structural variant under force_annotate"
    );

    let positional = tmp.path().join("positional.txt");
    let run = vep_on_cadd_corpus(
        &cadd,
        &positional,
        &format!("CADD,{},{},1", snv.display(), indels.display()),
    );
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    assert_eq!(
        cadd_fields(&std::fs::read_to_string(&positional).unwrap()),
        found,
        "the trailing 1 is no file"
    );

    // A key the plugin does not name opens its value as a data file, which marks
    // no sequence file (CADD.pm:124-129), so the SNV file under `file=` scores
    // the SNV sites and nothing else.
    let other_key_out = tmp.path().join("other_key.txt");
    let other_key = vep_on_cadd_corpus(
        &cadd,
        &other_key_out,
        &format!("CADD,file={}", snv.display()),
    );
    let stderr = String::from_utf8_lossy(&other_key.stderr);
    assert_eq!(other_key.status.code(), Some(0), "{stderr}");
    assert!(
        !stderr.contains("Failed to instantiate plugin CADD"),
        "{stderr}"
    );
    let other_key_fields = cadd_fields(&std::fs::read_to_string(&other_key_out).unwrap());
    assert_eq!(
        other_key_fields["intergenic_snv"],
        snv_scores["21:19268555:C:T"]
    );
    assert_eq!(
        other_key_fields["shared_bases_substitution"],
        snv_scores["21:30396091:C:T"]
    );
    assert_eq!(
        other_key_fields["deletion_rev"],
        (None, None),
        "no indel file under file="
    );
    assert_eq!(other_key_fields["deletion_sv"], (None, None));

    let sv_corpus = corpus("GRCh38-plugin-cadd-sv-fixture");
    let sv = sv_corpus.join("plugin_data/1000G_phase3_SVs_synthetic.tsv.gz");
    let sv_scores = cadd_file_scores(&sv);
    let sv_out = tmp.path().join("sv.txt");
    let run = vep_on_cadd_corpus(&sv_corpus, &sv_out, &format!("CADD,sv={}", sv.display()));
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let found = cadd_fields(&std::fs::read_to_string(&sv_out).unwrap());
    assert_eq!(found["deletion_sv"], sv_scores["21:5029000:DEL"]);
    assert_eq!(found["duplication_sv"], sv_scores["21:5030100:DUP"]);
    assert_eq!(found["insertion_sv"], sv_scores["21:5031000:INS"]);
    assert_eq!(
        found["transcript_deletion_sv"],
        sv_scores["21:30396000:DEL"]
    );
    assert_eq!(
        found["unscored_deletion_sv"],
        (None, None),
        "no row at its start"
    );
    assert_eq!(
        found["inversion_sv"],
        (None, None),
        "a class the file does not score"
    );
    assert_eq!(
        found["snv_inside_sv"],
        (None, None),
        "a sequence variant against a CADD-SV file alone"
    );
}

/// CADD and AlphaMissense refuse a score file of the other build by the run's
/// assembly, which the cache names where the command gives no `--assembly`
/// (`resolve` in `cache_resolve.rs`). The plugin is handed the assembly
/// and not where it was read, so the warning names the run's assembly and no
/// flag value, and it is the same line with the flag and without it. The run
/// goes on without the plugin.
#[test]
fn a_score_file_of_the_other_build_is_refused_by_the_runs_assembly_with_and_without_the_flag() {
    let tmp = tempfile::tempdir().unwrap();
    let cadd_file =
        corpus("GRCh37-plugin-cadd-fixture").join("plugin_data/whole_genome_SNVs_synthetic.tsv.gz");
    let alphamissense_file =
        corpus("GRCh37-plugin-alphamissense").join("plugin_data/alphamissense.tsv.gz");
    let cases = [
        (
            "GRCh38-plugin-cadd-fixture",
            format!("CADD,snv={}", cadd_file.display()),
            format!(
                "WARNING: Failed to instantiate plugin CADD: the run's assembly is GRCh38 \
                 but the data file '{}' does not name GRCh38 in its header. Supply the GRCh38 \
                 CADD file.",
                cadd_file.display()
            ),
            "## CADD_",
        ),
        (
            "GRCh38-plugin-alphamissense",
            format!("AlphaMissense,file={}", alphamissense_file.display()),
            format!(
                "WARNING: Failed to instantiate plugin AlphaMissense: the run's \
                 assembly is GRCh38 but '{}' declares genome 'hg19'. AlphaMissense publishes \
                 one file per build; point --plugin AlphaMissense at the GRCh38 file.",
                alphamissense_file.display()
            ),
            "## am_",
        ),
    ];
    for (cache_corpus, plugin, warning, header) in &cases {
        let dir = corpus(cache_corpus);
        for assembly_flag in [&[][..], &["--assembly", "GRCh38"][..]] {
            let out = tmp.path().join("refused.txt");
            let run = vep()
                .current_dir(&dir)
                .args([
                    "-i",
                    dir.join("variants.vcf").to_str().unwrap(),
                    "-o",
                    out.to_str().unwrap(),
                    "--offline",
                    "--json_cache",
                    dir.join("json_cache").to_str().unwrap(),
                    "--force_overwrite",
                    "--no_stats",
                    "--dir_plugins",
                    dir.join("plugin_data").to_str().unwrap(),
                    "--plugin",
                    plugin,
                ])
                .args(assembly_flag)
                .stdin(Stdio::null())
                .output()
                .unwrap();
            let stderr = String::from_utf8_lossy(&run.stderr);
            assert_eq!(
                run.status.code(),
                Some(0),
                "{plugin} {assembly_flag:?}: {stderr}"
            );
            assert_eq!(
                stderr
                    .lines()
                    .filter(|line| line.starts_with("WARNING:"))
                    .collect::<Vec<_>>(),
                [warning.as_str()],
                "{plugin} {assembly_flag:?}"
            );
            let output = std::fs::read_to_string(&out).unwrap();
            assert!(
                !output.lines().any(|line| line.starts_with(header)),
                "{plugin} {assembly_flag:?}: the plugin dropped at load writes no header line"
            );
            assert!(
                output.lines().any(|line| !line.starts_with('#')),
                "{plugin} {assembly_flag:?}: the run goes on without the plugin"
            );
        }
    }
}

/// A corpus's own arguments before the output path: its input, its cache and the
/// run's assembly, with the plugins directory the corpus carries.
fn args_for(dir: &Path, assembly: &str) -> Vec<String> {
    let corpus = corpora()
        .into_iter()
        .find(|corpus| corpus.dir.canonicalize().ok() == dir.canonicalize().ok())
        .unwrap_or_else(|| panic!("{} is not a corpus under tests/golden", dir.display()));
    let mut args = vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "--offline".to_string(),
    ];
    args.extend(corpus.cache_args());
    args.extend(
        [
            "--assembly",
            assembly,
            "--force_overwrite",
            "--no_stats",
            "--dir_plugins",
            dir.join("plugin_data").to_str().unwrap(),
        ]
        .map(str::to_string),
    );
    args
}

fn run(dir: &Path, assembly: &str, out_dir: &Path, extra: &[&str]) -> (Output, String) {
    let out = out_dir.join("out.txt");
    let mut cmd: Command = common::isolated_vep(out_dir);
    cmd.args(args_for(dir, assembly));
    cmd.arg("-o").arg(&out);
    cmd.args(extra);
    let output = cmd.output().unwrap();
    let text = std::fs::read_to_string(&out).unwrap_or_default();
    (output, text)
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).to_string()
}

fn rows_of<'a>(output: &'a str, record: &str) -> Vec<&'a str> {
    output
        .lines()
        .filter(|line| !line.starts_with('#') && line.starts_with(&format!("{record}\t")))
        .collect()
}

fn extra_field<'a>(row: &'a str, key: &str) -> Option<&'a str> {
    row.rsplit('\t')
        .next()?
        .split(';')
        .find_map(|pair| pair.strip_prefix(&format!("{key}=")))
}

/// `output_prefix=1` and `transcript_match=1` given as `key=value` parameters
/// beside the positional file and columns: the header names and the values
/// carry the `dbNSFP_` prefix, a transcript-specific column keeps the row's
/// transcript's entry alone, a column the file lists per transcript but the
/// row's transcript is absent from is dropped while a whole-row column stays,
/// and a row outside the default consequence set gets nothing. A parameter
/// name the Ensembl VEP plugin does not know is refused under `--safe`.
#[test]
fn dbnsfp_named_params() {
    let dir = corpus("GRCh38-plugin-dbnsfp-fixture");
    let file = dir.join("plugin_data/dbNSFP4.9a_grch38_synthetic.gz");
    let tmp = tempfile::tempdir().unwrap();
    let plugin = format!(
        "dbNSFP,{},output_prefix=1,transcript_match=1,SIFT_score,GERP++_RS",
        file.display()
    );
    let (output, body) = run(&dir, "GRCh38", tmp.path(), &["--plugin", &plugin]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        body.contains("## dbNSFP_SIFT_score : SIFT_score from dbNSFP file\n"),
        "{body}"
    );
    assert!(body.contains("## dbNSFP_GERP++_RS : GERP++_RS from dbNSFP file\n"));

    let rows = rows_of(&body, "missense_fwd_two_changes");
    let own_change = rows
        .iter()
        .find(|row| row.contains("\tENST00001126521\t"))
        .expect("the transcript with its own amino acid change");
    let sift = extra_field(own_change, "dbNSFP_SIFT_score").expect("a SIFT score");
    assert!(
        !sift.contains(',') && sift.parse::<f64>().is_ok(),
        "one entry for one transcript, got {sift}"
    );
    assert!(extra_field(own_change, "dbNSFP_GERP++_RS").is_some());
    let shared_change: Vec<&&str> = rows
        .iter()
        .filter(|row| row.contains("\tmissense_variant\t"))
        .filter(|row| !row.contains("\tENST00001126521\t"))
        .collect();
    assert_eq!(shared_change.len(), 6, "{rows:?}");
    for row in &shared_change {
        let sift = extra_field(row, "dbNSFP_SIFT_score").unwrap();
        assert!(!sift.contains(','), "{row}");
    }
    let intronic = rows
        .iter()
        .find(|row| row.contains("\tintron_variant\t"))
        .unwrap();
    assert!(!intronic.contains("dbNSFP_"), "{intronic}");

    let refused = format!("dbNSFP,{},pepmatch=0,SIFT_score", file.display());
    let (output, _) = run(
        &dir,
        "GRCh38",
        tmp.path(),
        &["--safe", "--plugin", &refused],
    );
    assert_eq!(output.status.code(), Some(1));
    let stderr = text(&output.stderr);
    assert!(
        stderr.lines().any(|line| line
            == "ERROR: [plugin] Failed to instantiate plugin dbNSFP: unknown parameter 'pepmatch'; dbNSFP accepts consequence, output_prefix, pep_match, transcript_match"),
        "{stderr}"
    );
}

/// `tests/fixtures/dbnsfp-cells`: a row whose `SIFT_score` list is `.;0.297`
/// and whose `SIFT_pred` cell is empty, under `transcript_match=1`. The first
/// transcript's `.` entry is left out and the empty cell is a value the row has
/// (dbNSFP.pm:388, :407-408): `SIFT_pred=` in the default format, an empty cell
/// in the tab format (Tab.pm:128) and `""` in JSON (JSON.pm:327), the rows
/// Ensembl VEP 116.3 wrote on the fixture; a score the row lacks is `-` in the
/// tab format and no JSON key.
#[test]
fn dbnsfp_empty_cell_and_dot_entry_in_every_format() {
    let corpus_dir = corpus("GRCh38-plugin-dbnsfp-fixture");
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dbnsfp-cells");
    let plugin = format!(
        "dbNSFP,{},transcript_match=1,SIFT_score,SIFT_pred",
        fixture
            .join("dbNSFP4.9a_grch38_cells_synthetic.gz")
            .display()
    );
    let tmp = tempfile::tempdir().unwrap();
    let run_format = |format: &[&str]| -> String {
        let out = tmp.path().join("cells.out");
        let mut cmd = common::isolated_vep(tmp.path());
        cmd.arg("-i")
            .arg(fixture.join("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(corpus_dir.join("json_cache"))
            .args([
                "--offline",
                "--species",
                "homo_sapiens",
                "--assembly",
                "GRCh38",
            ])
            .args([
                "--force_overwrite",
                "--no_stats",
                "--quiet",
                "--plugin",
                &plugin,
            ])
            .args(format);
        let output = cmd.output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        std::fs::read_to_string(&out).unwrap()
    };
    let transcript_rows = |body: &str| -> Vec<String> {
        body.lines()
            .filter(|line| {
                line.contains("\tENST00000290200\t") || line.contains("\tENST00000896210\t")
            })
            .map(str::to_string)
            .collect()
    };

    let default = transcript_rows(&run_format(&[]));
    assert_eq!(default.len(), 2, "{default:?}");
    assert!(
        default[0].ends_with("\tIMPACT=MODERATE;STRAND=1;SIFT_pred="),
        "{}",
        default[0]
    );
    assert!(
        default[1].ends_with("\tIMPACT=MODERATE;STRAND=1;SIFT_pred=;SIFT_score=0.297"),
        "{}",
        default[1]
    );

    let tab = run_format(&["--tab"]);
    let header: Vec<&str> = tab
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap()
        .split('\t')
        .collect();
    let score = header.iter().position(|c| *c == "SIFT_score").unwrap();
    let pred = header.iter().position(|c| *c == "SIFT_pred").unwrap();
    let cells: Vec<(String, String)> = transcript_rows(&tab)
        .iter()
        .map(|row| {
            let cells: Vec<&str> = row.split('\t').collect();
            (cells[score].to_string(), cells[pred].to_string())
        })
        .collect();
    assert_eq!(
        cells,
        [
            ("-".to_string(), String::new()),
            ("0.297".to_string(), String::new())
        ]
    );

    let json = run_format(&["--json"]);
    let record: serde_json::Value = serde_json::from_str(json.lines().next().unwrap()).unwrap();
    let mut sift: Vec<(String, Option<String>, Option<String>)> = record["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|tc| {
            matches!(
                tc["transcript_id"].as_str(),
                Some("ENST00000290200" | "ENST00000896210")
            )
        })
        .map(|tc| {
            (
                tc["transcript_id"].as_str().unwrap().to_string(),
                tc.get("sift_score").map(|v| v.to_string()),
                tc.get("sift_pred").map(|v| v.to_string()),
            )
        })
        .collect();
    sift.sort();
    assert_eq!(
        sift,
        [
            (
                "ENST00000290200".to_string(),
                None,
                Some("\"\"".to_string())
            ),
            (
                "ENST00000896210".to_string(),
                Some("0.297".to_string()),
                Some("\"\"".to_string())
            ),
        ]
    );
}

/// The cells a plugin's field `field` takes on the rows of each record of a
/// `*-cells` fixture, in the four formats, keyed by the record's id: the
/// default format's item (`KEY=value`, or `no item`), the tab cell, the CSQ
/// field and the JSON value (its JSON text, or `no key`). Every row of a record
/// is collected, so a set of one cell says all its rows agree.
fn field_cells_in_every_format(
    fixture: &Path,
    corpus_dir: &Path,
    assembly: &str,
    plugin: &str,
    field: &str,
) -> [BTreeMap<String, BTreeSet<String>>; 4] {
    let tmp = tempfile::tempdir().unwrap();
    let run_format = |format: &[&str]| -> String {
        let out = tmp.path().join("cells.out");
        let mut cmd = common::isolated_vep(tmp.path());
        cmd.arg("-i")
            .arg(fixture.join("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(corpus_dir.join("json_cache"))
            .args([
                "--offline",
                "--species",
                "homo_sapiens",
                "--assembly",
                assembly,
            ])
            .args([
                "--force_overwrite",
                "--no_stats",
                "--quiet",
                "--plugin",
                plugin,
            ])
            .args(format);
        let output = cmd.output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        std::fs::read_to_string(&out).unwrap()
    };
    let mut cells: [BTreeMap<String, BTreeSet<String>>; 4] = Default::default();
    let item = format!("{field}=");
    for line in run_format(&[]).lines().filter(|l| !l.starts_with('#')) {
        let columns: Vec<&str> = line.split('\t').collect();
        let cell = columns[13]
            .split(';')
            .find(|entry| entry.starts_with(&item))
            .map_or_else(|| "no item".to_string(), str::to_string);
        cells[0]
            .entry(columns[0].to_string())
            .or_default()
            .insert(cell);
    }
    let tab = run_format(&["--tab"]);
    let header: Vec<&str> = tab
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap()
        .split('\t')
        .collect();
    let column = header.iter().position(|c| *c == field).unwrap();
    for line in tab.lines().filter(|l| !l.starts_with('#')) {
        let columns: Vec<&str> = line.split('\t').collect();
        cells[1]
            .entry(columns[0].to_string())
            .or_default()
            .insert(columns[column].to_string());
    }
    let vcf = run_format(&["--vcf"]);
    let csq_fields: Vec<&str> = vcf
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ"))
        .and_then(|line| line.split("Format: ").nth(1))
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .collect();
    let csq_field = csq_fields.iter().position(|c| *c == field).unwrap();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let columns: Vec<&str> = line.split('\t').collect();
        let csq = columns[7]
            .split(';')
            .find_map(|entry| entry.strip_prefix("CSQ="))
            .unwrap();
        for row in csq.split(',') {
            cells[2]
                .entry(columns[2].to_string())
                .or_default()
                .insert(row.split('|').nth(csq_field).unwrap().to_string());
        }
    }
    let key = field.to_ascii_lowercase();
    for line in run_format(&["--json"]).lines() {
        let record: serde_json::Value = serde_json::from_str(line).unwrap();
        let id = record["id"].as_str().unwrap().to_string();
        for kind in [
            "transcript_consequences",
            "regulatory_feature_consequences",
            "motif_feature_consequences",
            "intergenic_consequences",
        ] {
            for row in record[kind].as_array().into_iter().flatten() {
                // A SpliceVault field nests under `splicevault` with its own name.
                let cell = row
                    .get(&key)
                    .or_else(|| row.get("splicevault").and_then(|object| object.get(field)))
                    .map_or_else(|| "no key".to_string(), serde_json::Value::to_string);
                cells[3].entry(id.clone()).or_default().insert(cell);
            }
        }
    }
    cells
}

/// The four cells each record of a `*-cells` fixture takes in the four formats,
/// in the order default item, tab cell, CSQ field, JSON text, as one set each.
fn expected_cells(rows: &[(&str, [&str; 4])]) -> [BTreeMap<String, BTreeSet<String>>; 4] {
    let mut cells: [BTreeMap<String, BTreeSet<String>>; 4] = Default::default();
    for (id, values) in rows {
        for (format, value) in values.iter().enumerate() {
            cells[format]
                .entry(id.to_string())
                .or_default()
                .insert(value.to_string());
        }
    }
    cells
}

/// `tests/fixtures/clinpred-cells`: two rows whose fifth cell is empty, a row
/// of four cells and a scored row. `ClinPred.pm:127` splits a row into five
/// variables, so Ensembl VEP keeps a trailing empty cell as the empty string and
/// leaves a missing cell undefined, and `:116-122` returns the matched row's
/// result either way: an empty cell is `ClinPred=`, an empty tab cell, an
/// empty CSQ field and `""` in JSON, and the four-cell row is `ClinPred=`, `-`,
/// an empty CSQ field and no JSON key (`VEP_output.pm:117-119`, `Tab.pm:128`,
/// `JSON.pm:327`), the rows Ensembl VEP 116.3 wrote on the fixture.
#[test]
fn clinpred_short_row_and_empty_cell_in_every_format() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/clinpred-cells");
    let plugin = format!(
        "ClinPred,file={}",
        fixture.join("ClinPred_cells_synthetic.tsv.gz").display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-clinpred-fixture"),
        "GRCh38",
        &plugin,
        "ClinPred",
    );
    assert_eq!(
        cells,
        expected_cells(&[
            ("empty_fifth_of_six", ["ClinPred=", "", "", "\"\""]),
            ("empty_cell", ["ClinPred=", "", "", "\"\""]),
            ("short_row", ["ClinPred=", "-", "", "no key"]),
            (
                "scored",
                ["ClinPred=0.411650", "0.411650", "0.411650", "0.41165"]
            ),
        ])
    );
}

/// `tests/fixtures/bayesdel-cells`: the same four row shapes for BayesDel,
/// whose `BayesDel.pm:151` splits into five variables the same way and whose
/// `:146` returns the matched row's result either way; the result prints on
/// every row of the record (`:117`), so the intergenic record's row carries
/// the empty string too. The cells are the rows Ensembl VEP 116.3 wrote on the
/// fixture, as for ClinPred.
#[test]
fn bayesdel_short_row_and_empty_cell_in_every_format() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/bayesdel-cells");
    let plugin = format!(
        "BayesDel,file={}",
        fixture.join("BayesDel_cells_synthetic.txt.gz").display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh37-plugin-bayesdel-fixture"),
        "GRCh37",
        &plugin,
        "BayesDel",
    );
    assert_eq!(
        cells,
        expected_cells(&[
            ("empty_fifth_of_six", ["BayesDel=", "", "", "\"\""]),
            ("empty_cell", ["BayesDel=", "", "", "\"\""]),
            ("short_row", ["BayesDel=", "-", "", "no key"]),
            (
                "scored",
                ["BayesDel=-0.55399", "-0.55399", "-0.55399", "-0.55399"]
            ),
        ])
    );
}

/// `tests/fixtures/primateai-cells`: a table whose header names a column after
/// the score, with two rows whose score cell is empty before that filled cell
/// and a scored row. `PrimateAI.pm:188` splits a row with no limit, so the
/// empty score is defined, passes `:198` and is returned (`:204-207`):
/// `PrimateAI=`, an empty tab cell, an empty CSQ field and `""` in JSON, on a
/// transcript row and on the intergenic row alike, the cells Ensembl VEP 116.3
/// wrote on the fixture.
#[test]
fn primateai_empty_score_before_a_filled_cell_in_every_format() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/primateai-cells");
    let plugin = format!(
        "PrimateAI,{}",
        fixture
            .join("PrimateAI_scores_cells_GRCh38_synthetic.tsv.bgz")
            .display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-primateai-fixture"),
        "GRCh38",
        &plugin,
        "PrimateAI",
    );
    assert_eq!(
        cells,
        expected_cells(&[
            ("empty_score_intergenic", ["PrimateAI=", "", "", "\"\""]),
            ("empty_score", ["PrimateAI=", "", "", "\"\""]),
            ("scored", ["PrimateAI=0.6125", "0.6125", "0.6125", "0.6125"]),
        ])
    );
}

/// `tests/fixtures/bayesdel-mitochondrial-contig`: two records on `MT` of
/// GRCh37, one with a row in the score file and one without, against a file
/// naming the contig `M` and a file naming it `MT`. `MT` and `chrMT` name one
/// mitochondrial sequence there and `chrM` and `M` another, and
/// `BayesDel.pm:144` hands the record's own name to
/// `BaseVepTabixPlugin.pm:693-708`, which read the file under that name with
/// or without `chr`: the record with a row takes its score from the file
/// naming the contig `MT` and none from the file naming it `M`, whose rows
/// the run, given no `--fasta`, cannot compare with a reference sequence. The
/// two headers spell or name their second and fifth columns otherwise than the
/// generator's header for the shape, and `BayesDel.pm:151` takes a row's
/// cells by place. Each output is the one Ensembl VEP 116.3 wrote on the
/// fixture under `--no_headers`, byte for byte.
#[test]
fn bayesdel_reads_a_mitochondrial_contig_of_the_records_own_pair_on_grch37() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bayesdel-mitochondrial-contig");
    let tmp = tempfile::tempdir().unwrap();
    for (pair, rows_with_a_score) in [("other_pair", 0), ("own_pair", 22)] {
        let expected =
            std::fs::read_to_string(fixture.join(format!("ensembl_vep.{pair}.txt"))).unwrap();
        let scored = expected
            .lines()
            .filter(|row| row.contains("BayesDel="))
            .count();
        assert_eq!(scored, rows_with_a_score, "{pair}: the reference's rows");
        let out = tmp.path().join(format!("{pair}.txt"));
        let plugin = format!(
            "BayesDel,file={}",
            fixture
                .join(format!("scores_{pair}_synthetic.txt.gz"))
                .display()
        );
        let output = common::isolated_vep(tmp.path())
            .arg("-i")
            .arg(fixture.join("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(fixture.join("json_cache"))
            .args([
                "--offline",
                "--species",
                "homo_sapiens",
                "--assembly",
                "GRCh37",
                "--buffer_size",
                "5000",
                "--force_overwrite",
                "--no_stats",
                "--quiet",
                "--no_headers",
                "--plugin",
                &plugin,
            ])
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{pair}: {}",
            text(&output.stderr)
        );
        let actual = std::fs::read_to_string(&out).unwrap();
        assert!(
            actual == expected,
            "{pair}: {}",
            common::compare_headerless(&expected, &actual)
        );
    }
}

/// `tests/fixtures/mitochondrial-frame-check`: two made-up records on `MT`
/// and three made-up BayesDel-shaped tables over a made-up reference. On
/// GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence and
/// `chrM` and `M` another one, so a table that names its contig `M` is read
/// for a record on `MT` only where its rows state the bases `--fasta` holds
/// at their own positions. Against the table whose rows state them the run
/// writes the rows Ensembl VEP 116.3 wrote against the same rows named `MT`,
/// which `BayesDel.pm` line 144 and `BaseVepTabixPlugin.pm` lines 693 to 708
/// reach; against the table whose rows state the bases one position on, and
/// in a run without `--fasta`, it writes no `BayesDel` item, as Ensembl VEP
/// wrote none against either table named `M`. The comparison is logged once a
/// table with its counts, and the run that cannot make it warns once, naming
/// the table.
#[test]
fn on_grch37_a_table_of_the_other_mitochondrial_pair_is_read_where_its_rows_state_the_reference() {
    let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let fixture = fixtures.join("mitochondrial-frame-check");
    let cache = fixtures.join("bayesdel-mitochondrial-contig/json_cache");
    let tmp = tempfile::tempdir().unwrap();
    for (table, with_fasta, expected, rows_with_a_score, logged) in [
        (
            "stating_the_reference",
            true,
            "own_name",
            44,
            Some(
                "contig=\"M\" rows_compared=128 rows_stating_the_reference=128 positions=48 \
                 outcome=\"read across the mitochondrial names\"",
            ),
        ),
        (
            "one_position_on",
            true,
            "one_position_on",
            0,
            Some(
                "contig=\"M\" rows_compared=128 rows_stating_the_reference=39 positions=48 \
                 outcome=\"read under its own pair of names alone\"",
            ),
        ),
        (
            "stating_the_reference",
            false,
            "stating_the_reference",
            0,
            None,
        ),
    ] {
        let case = format!("{table}, --fasta {with_fasta}");
        let expected =
            std::fs::read_to_string(fixture.join(format!("ensembl_vep.{expected}.txt"))).unwrap();
        let scored = expected
            .lines()
            .filter(|row| row.contains("BayesDel="))
            .count();
        assert_eq!(scored, rows_with_a_score, "{case}: the reference's rows");
        let out = tmp.path().join(format!("{table}.{with_fasta}.txt"));
        let score_file = fixture.join(format!("scores_{table}_synthetic.txt.gz"));
        let plugin = format!("BayesDel,file={}", score_file.display());
        let mut command = common::isolated_vep(tmp.path());
        command
            .arg("-i")
            .arg(fixture.join("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(&cache)
            .args([
                "--offline",
                "--species",
                "homo_sapiens",
                "--assembly",
                "GRCh37",
                "--buffer_size",
                "5000",
                "--force_overwrite",
                "--no_stats",
                "--no_headers",
                "--plugin",
                &plugin,
            ]);
        if with_fasta {
            command
                .arg("--fasta")
                .arg(fixture.join("reference_synthetic.fa"));
        }
        let output = command.output().unwrap();
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{case}: {stderr}");
        let actual = std::fs::read_to_string(&out).unwrap();
        assert!(
            actual == expected,
            "{case}: {}",
            common::compare_headerless(&expected, &actual)
        );
        let compared: Vec<&str> = stderr
            .lines()
            .filter(|line| line.contains("mitochondrial rows compared with the reference sequence"))
            .collect();
        let warned: Vec<&str> = stderr
            .lines()
            .filter(|line| line.contains("could not be compared with the reference sequence"))
            .collect();
        match logged {
            Some(counts) => {
                assert_eq!(compared.len(), 1, "{case}: {stderr}");
                assert!(
                    compared[0].ends_with(&format!("table={} {counts}", score_file.display())),
                    "{case}: {}",
                    compared[0]
                );
                assert!(warned.is_empty(), "{case}: {stderr}");
            }
            None => {
                assert_eq!(
                    warned,
                    [format!(
                        "WARNING: {}: its contig M is not read for a mitochondrial record \
                         written MT or chrMT: on GRCh37 the names MT and chrMT are one sequence \
                         and chrM and M another, and the table's rows could not be compared \
                         with the reference sequence (no reference sequence was given; pass \
                         --fasta)",
                        score_file.display()
                    )],
                    "{case}: {stderr}"
                );
            }
        }
    }
}

/// `tests/fixtures/avi-cells`: two rows whose sixth cell is empty, a row of
/// five cells and a scored row, under `raw=1,links=1`. `AVI.pm:161` splits a
/// row into six variables with a limit of seven, so Ensembl VEP keeps a trailing
/// empty cell as the empty string, and `:163` drops a row lacking its sixth
/// cell before any match: an empty cell is `AVI_PHRED=`, an empty tab cell, an
/// empty CSQ field and `""` in JSON, and the five-cell row takes no field, so
/// `-` in the tab format, an empty CSQ field and no JSON key, with no default
/// item (`VEP_output.pm:117-119`, `Tab.pm:128`, `JSON.pm:327`); the rows
/// Ensembl VEP 116.3 wrote on the fixture, which print on the intergenic row.
#[test]
fn avi_short_row_and_empty_cell_in_every_format() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/avi-cells");
    let plugin = format!(
        "AVI,file={},raw=1,links=1",
        fixture.join("AVI_cells.synthetic.tsv.gz").display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-avi-fixture"),
        "GRCh38",
        &plugin,
        "AVI_PHRED",
    );
    assert_eq!(
        cells,
        expected_cells(&[
            ("empty_sixth_of_seven", ["AVI_PHRED=", "", "", "\"\""]),
            ("empty_cell", ["AVI_PHRED=", "", "", "\"\""]),
            ("short_row", ["no item", "-", "", "no key"]),
            ("scored", ["AVI_PHRED=1.446", "1.446", "1.446", "1.446"]),
        ])
    );
}

/// `tests/fixtures/cadd-lower-case-alleles`: a made-up CADD table whose row
/// for the first record writes its alleles in lower case and whose row for the
/// second writes them in capitals. Ensembl VEP compares the cells as written
/// (`CADD.pm:239-251`, `Utils/Sequence.pm:1207`, `:1232`) and scores the
/// second record alone; vep-rs reads a row in either letter case and scores
/// both (the documented class `tabix_allele_case_insensitive_match`).
#[test]
fn cadd_scores_a_row_whose_alleles_are_in_lower_case() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cadd-lower-case-alleles");
    let plugin = format!(
        "CADD,snv={}",
        fixture
            .join("whole_genome_SNVs_lower_case.synthetic.tsv.gz")
            .display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-cadd-fixture"),
        "GRCh38",
        &plugin,
        "CADD_PHRED",
    );
    assert_eq!(
        cells,
        expected_cells(&[
            (
                "lower_case_row",
                ["CADD_PHRED=7.654", "7.654", "7.654", "7.654"]
            ),
            (
                "upper_case_row",
                ["CADD_PHRED=6.543", "6.543", "6.543", "6.543"]
            ),
        ])
    );
}

/// `tests/fixtures/opentargets-lower-case-alleles`: a made-up export whose
/// record for the first variant writes its alleles in lower case. Ensembl VEP
/// compares the cells as written (`OpenTargets.pm:376-388`,
/// `Utils/Sequence.pm:1207`, `:1232`) and writes its fields on the second
/// variant alone; vep-rs writes them on both. The JSON format nests the
/// fields, so the default, tab and VCF cells are read.
#[test]
fn opentargets_matches_a_record_whose_alleles_are_in_lower_case() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/opentargets-lower-case-alleles");
    let plugin = format!(
        "OpenTargets,file={}",
        fixture
            .join("open_targets_vep_lower_case.synthetic.tsv.gz")
            .display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-opentargets"),
        "GRCh38",
        &plugin,
        "OpenTargets_gwasGeneId",
    );
    let expected = expected_cells(&[
        (
            "lower_case_row",
            [
                "OpenTargets_gwasGeneId=ENSG00000000008",
                "ENSG00000000008",
                "ENSG00000000008",
                "",
            ],
        ),
        (
            "upper_case_row",
            [
                "OpenTargets_gwasGeneId=ENSG00000000009",
                "ENSG00000000009",
                "ENSG00000000009",
                "",
            ],
        ),
    ]);
    assert_eq!(cells[..3], expected[..3]);
}

/// `tests/fixtures/gwas-lower-case-alleles`: two rows of the GWAS corpus's
/// catalog slice, the first with its risk allele in lower case. Ensembl VEP
/// compares the allele as written (`GWAS.pm:182-194`, `Utils/Sequence.pm:1207`,
/// `:1232`) and writes its fields on the second record alone; vep-rs writes
/// them on both, the risk allele printed as the file writes it.
#[test]
fn gwas_matches_a_risk_allele_in_lower_case() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/gwas-lower-case-alleles");
    let plugin = format!(
        "GWAS,file={}",
        fixture.join("gwas_catalog_lower_case.tsv").display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-gwas"),
        "GRCh38",
        &plugin,
        "GWAS_risk_allele",
    );
    assert_eq!(
        cells,
        expected_cells(&[
            ("lower_case_row", ["GWAS_risk_allele=t", "t", "t", "\"t\""]),
            ("upper_case_row", ["GWAS_risk_allele=G", "G", "G", "\"G\""]),
        ])
    );
}

/// `tests/fixtures/loeuf-trailing-space`: a made-up table whose last column is
/// the LOEUF column, the first row's value followed by one space. Ensembl VEP
/// hands the plugin the line with that space (`BaseVepTabixPlugin.pm:605-606`)
/// and the plugin prints the cell (`LOEUF.pm:181`, `:215`), so under
/// `match_by=transcript` the first row's transcript reads the value with its
/// space and the second its own; under `match_by=gene` the cell is a number
/// (`:192`) and both read the minimum as its cell is written (`:204`). The
/// default-format items are the ones Ensembl VEP 116.3 wrote on the fixture.
#[test]
fn loeuf_prints_a_space_after_a_last_cell() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/loeuf-trailing-space");
    let table = fixture.join("loeuf_narrow_synthetic.tsv.gz");
    let item = |match_by: &str, transcript: &str| -> String {
        let cells = feature_cells_in_every_format(
            &fixture,
            &corpus("GRCh38-plugin-loeuf"),
            &[format!(
                "LOEUF,file={},match_by={match_by}",
                table.display()
            )],
            "LOEUF",
        );
        cells[transcript][0].clone()
    };
    assert_eq!(item("transcript", "ENST00000346798"), "LOEUF=5.0000e-01 ");
    assert_eq!(item("transcript", "ENST00000348990"), "LOEUF=7.0000e-01");
    assert_eq!(item("transcript", "ENST00000354192"), "no item");
    assert_eq!(item("gene", "ENST00000346798"), "LOEUF=5.0000e-01 ");
    assert_eq!(item("gene", "ENST00000348990"), "LOEUF=5.0000e-01 ");
}

/// The cells the field `field` takes on each row of a one-record fixture under
/// the `--plugin` arguments `plugins`, in their order, keyed by the row's
/// feature id: the default format's item (`KEY=value`, or `no item`), the tab
/// cell, the CSQ field and the JSON value (its JSON text, or `no key`); `no row`
/// where a format wrote no row for the feature.
fn feature_cells_in_every_format(
    fixture: &Path,
    corpus_dir: &Path,
    plugins: &[String],
    field: &str,
) -> BTreeMap<String, [String; 4]> {
    let tmp = tempfile::tempdir().unwrap();
    let run_format = |format: &[&str]| -> String {
        let out = tmp.path().join("cells.out");
        let mut cmd = common::isolated_vep(tmp.path());
        cmd.arg("-i")
            .arg(fixture.join("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(corpus_dir.join("json_cache"))
            .args([
                "--offline",
                "--species",
                "homo_sapiens",
                "--assembly",
                "GRCh38",
            ])
            .args(["--force_overwrite", "--no_stats", "--quiet"]);
        for plugin in plugins {
            cmd.args(["--plugin", plugin]);
        }
        let output = cmd.args(format).output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        std::fs::read_to_string(&out).unwrap()
    };
    let mut cells: BTreeMap<String, [String; 4]> = BTreeMap::new();
    let mut set = |feature: &str, format: usize, cell: String| {
        cells
            .entry(feature.to_string())
            .or_insert_with(|| std::array::from_fn(|_| "no row".to_string()))[format] = cell;
    };
    let item = format!("{field}=");
    for line in run_format(&[]).lines().filter(|l| !l.starts_with('#')) {
        let columns: Vec<&str> = line.split('\t').collect();
        let cell = columns[13]
            .split(';')
            .find(|entry| entry.starts_with(&item))
            .map_or_else(|| "no item".to_string(), str::to_string);
        set(columns[4], 0, cell);
    }
    let tab = run_format(&["--tab"]);
    let header: Vec<&str> = tab
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap()
        .split('\t')
        .collect();
    let column = |name: &str| header.iter().position(|c| *c == name).unwrap();
    for line in tab.lines().filter(|l| !l.starts_with('#')) {
        let columns: Vec<&str> = line.split('\t').collect();
        set(
            columns[column("Feature")],
            1,
            columns[column(field)].to_string(),
        );
    }
    let vcf = run_format(&["--vcf"]);
    let csq_fields: Vec<&str> = vcf
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ"))
        .and_then(|line| line.split("Format: ").nth(1))
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .collect();
    let csq_field = |name: &str| csq_fields.iter().position(|c| *c == name).unwrap();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let csq = line
            .split('\t')
            .nth(7)
            .unwrap()
            .split(';')
            .find_map(|entry| entry.strip_prefix("CSQ="))
            .unwrap();
        for row in csq.split(',') {
            let fields: Vec<&str> = row.split('|').collect();
            set(
                fields[csq_field("Feature")],
                2,
                fields[csq_field(field)].to_string(),
            );
        }
    }
    let key = field.to_ascii_lowercase();
    for line in run_format(&["--json"]).lines() {
        let record: serde_json::Value = serde_json::from_str(line).unwrap();
        for row in record["transcript_consequences"]
            .as_array()
            .into_iter()
            .flatten()
        {
            let cell = row
                .get(&key)
                .map_or_else(|| "no key".to_string(), serde_json::Value::to_string);
            set(row["transcript_id"].as_str().unwrap(), 3, cell);
        }
    }
    cells
}

/// The transcript the tables of `tests/fixtures/plugin-shared-key` name, one
/// of the nine the fixture's record overlaps.
const SHARED_KEY_TRANSCRIPT: &str = "ENST00000618007";

/// The cells of `feature_cells_in_every_format` as two groups: the row of
/// [`SHARED_KEY_TRANSCRIPT`], and the one set of cells the other rows share.
fn named_row_and_the_rest(mut cells: BTreeMap<String, [String; 4]>) -> ([String; 4], [String; 4]) {
    let named = cells.remove(SHARED_KEY_TRANSCRIPT).unwrap();
    let rest: BTreeSet<[String; 4]> = cells.into_values().collect();
    assert_eq!(rest.len(), 1, "the other rows agree: {rest:?}");
    (named, rest.into_iter().next().unwrap())
}

/// `tests/fixtures/plugin-shared-key`, BayesDel and then TranscriptAnnotator
/// under `prefix=0` with a column headed `BayesDel`: two plugins that write
/// one key. BayesDel's row for the record has its score cell present and empty,
/// the empty string on every row (`BayesDel.pm:151`); the table's row ends
/// after its five leading cells, so TranscriptAnnotator returns the key
/// undefined for the transcript it names (`TranscriptAnnotator.pm:245-253`) and
/// that later value stands on that row (`OutputFactory.pm:2239-2241`): `-` in
/// the tab format and no JSON key there, the empty string on the other rows,
/// the cells Ensembl VEP 116.3 wrote on the fixture.
#[test]
fn a_later_plugins_undefined_value_replaces_an_earlier_empty_string() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/plugin-shared-key");
    let plugins = [
        format!(
            "BayesDel,file={}",
            fixture
                .join("BayesDel_shared_key_synthetic.txt.gz")
                .display()
        ),
        format!(
            "TranscriptAnnotator,file={},prefix=0",
            fixture
                .join("transcript_annotator_BayesDel.tsv.gz")
                .display()
        ),
    ];
    let (named, rest) = named_row_and_the_rest(feature_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-transcript-annotator"),
        &plugins,
        "BayesDel",
    ));
    assert_eq!(named, ["BayesDel=", "-", "", "no key"]);
    assert_eq!(rest, ["BayesDel=", "", "", "\"\""]);
}

/// `tests/fixtures/plugin-shared-key`, TranscriptAnnotator under `prefix=CADD`
/// with a column headed `PHRED` and the CADD plugin: both write `CADD_PHRED`,
/// and the plugin loaded later stands on a row (`OutputFactory.pm:2239-2241`).
/// With CADD second its score prints on every row, the one the table names
/// included; with CADD first that row prints the table's cell and the others
/// the score, the cells Ensembl VEP 116.3 wrote on the fixture.
#[test]
fn the_plugin_loaded_later_stands_under_a_key_two_plugins_write() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/plugin-shared-key");
    let annotator = format!(
        "TranscriptAnnotator,file={},prefix=CADD",
        fixture.join("transcript_annotator_PHRED.tsv.gz").display()
    );
    let cadd = format!(
        "CADD,{}",
        fixture
            .join("CADD_shared_key_GRCh38_synthetic.tsv.gz")
            .display()
    );
    let cells = |plugins: [&String; 2]| {
        named_row_and_the_rest(feature_cells_in_every_format(
            &fixture,
            &corpus("GRCh38-plugin-transcript-annotator"),
            &plugins.map(String::clone),
            "CADD_PHRED",
        ))
    };
    let score = ["CADD_PHRED=23.1", "23.1", "23.1", "23.1"].map(str::to_string);
    let table_cell = ["CADD_PHRED=7.5", "7.5", "7.5", "7.5"].map(str::to_string);
    assert_eq!(cells([&annotator, &cadd]), (score.clone(), score.clone()));
    assert_eq!(cells([&cadd, &annotator]), (table_cell, score));
}

/// `tests/fixtures/splicevault-cells`: two rows whose twelfth cell is empty, a
/// row of eleven cells, a scored row, a row whose events cell is empty, a row
/// whose count and depth cells carry a blank at one end, a row whose first
/// event ends in its separator and a row whose second event carries a `,`
/// inside a part. `SpliceVault.pm:201` splits a row into twelve
/// variables, so Ensembl VEP keeps a trailing empty cell as the empty string and
/// leaves a missing cell undefined, and `:206-220` write both: an empty cell is
/// `SpliceVault_site_max_depth=`, an empty tab cell, an empty CSQ field and `""`
/// in JSON, and the eleven-cell row is `SpliceVault_site_max_depth=`, `-`, an
/// empty CSQ field and `null` inside the `splicevault` object, which the JSON
/// format encodes as the plugin's hash stands (`VEP_output.pm:117-119`,
/// `Tab.pm:128`; `JSON.pm:327` drops top-level keys alone). The events cell is
/// a list either way (`:216`), empty when the cell is: `SpliceVault_top_events=`,
/// an empty tab cell, an empty CSQ field and `{}`. A padded cell prints with its
/// blank in the default and tab formats, with `_` for it in the VCF format
/// (`VCF.pm:403`) and as the number in JSON, which `numberify` (`Utils.pm:231`)
/// converts through `looks_like_number`, blanks at either end admitted. An
/// event's trailing empty part is the empty string (`:184` splits an event into
/// five variables), so its `frame` is `""` in JSON where the text formats print
/// the event ending in `:`. That split reads `:` alone, so a `,` inside a part
/// is text: the text formats print it inside the `,`-joined list (`&` for it in
/// the VCF format, `VCF.pm:401`), and in JSON the event is one, its percent
/// `"1,5%"`. The rows Ensembl VEP 116.3 wrote on the fixture.
#[test]
fn splicevault_short_row_and_empty_cell_in_every_format() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/splicevault-cells");
    let plugin = format!(
        "SpliceVault,file={}",
        fixture.join("SpliceVault_cells_synthetic.tsv.gz").display()
    );
    let cells = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-splicevault-fixture"),
        "GRCh38",
        &plugin,
        "SpliceVault_site_max_depth",
    );
    // The record whose events cell is empty has two more transcript rows, and
    // the record whose event carries a `,` one more, that the table does not
    // annotate, which carry the absent-field cell.
    let unannotated = ["no item", "-", "", "no key"];
    let with_unannotated_rows = |mut expected: [BTreeMap<String, BTreeSet<String>>; 4]| {
        for (format, cell) in unannotated.into_iter().enumerate() {
            for id in ["empty_events", "comma_in_part"] {
                expected[format]
                    .get_mut(id)
                    .unwrap()
                    .insert(cell.to_string());
            }
        }
        expected
    };
    assert_eq!(
        cells,
        with_unannotated_rows(expected_cells(&[
            (
                "empty_twelfth_of_thirteen",
                ["SpliceVault_site_max_depth=", "", "", "\"\""]
            ),
            (
                "comma_in_part",
                ["SpliceVault_site_max_depth=22", "22", "22", "22"]
            ),
            (
                "empty_cell",
                ["SpliceVault_site_max_depth=", "", "", "\"\""]
            ),
            (
                "short_row",
                ["SpliceVault_site_max_depth=", "-", "", "null"]
            ),
            (
                "scored",
                ["SpliceVault_site_max_depth=93", "93", "93", "93"]
            ),
            (
                "empty_events",
                ["SpliceVault_site_max_depth=191", "191", "191", "191"]
            ),
            (
                "padded_cells",
                ["SpliceVault_site_max_depth=26 ", "26 ", "26_", "26"]
            ),
            (
                "empty_event_part",
                ["SpliceVault_site_max_depth=28", "28", "28", "28"]
            ),
        ]))
    );
    let counts = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-splicevault-fixture"),
        "GRCh38",
        &plugin,
        "SpliceVault_site_sample_count",
    );
    let expected = expected_cells(&[(
        "padded_cells",
        [
            "SpliceVault_site_sample_count= 33837",
            " 33837",
            "_33837",
            "33837",
        ],
    )]);
    for (format, cells) in counts.iter().enumerate() {
        assert_eq!(
            cells["padded_cells"], expected[format]["padded_cells"],
            "format {format}"
        );
    }
    let events = field_cells_in_every_format(
        &fixture,
        &corpus("GRCh38-plugin-splicevault-fixture"),
        "GRCh38",
        &plugin,
        "SpliceVault_top_events",
    );
    let trailing_empty_part = serde_json::json!({
        "Top1": {"site_type": "CA", "transcript_impact": -480,
                 "percent_of_supporting_samples": "7%", "frame": ""},
        "Top2": {"site_type": "CA", "transcript_impact": -76,
                 "percent_of_supporting_samples": "39%", "frame": "Frameshift"},
    })
    .to_string();
    let comma_in_part = serde_json::json!({
        "Top1": {"site_type": "CD", "transcript_impact": -330,
                 "percent_of_supporting_samples": "73_314%", "frame": "inFrame"},
        "Top2": {"site_type": "CA", "transcript_impact": 112,
                 "percent_of_supporting_samples": "1,5%", "frame": "inFrame"},
    })
    .to_string();
    let expected = with_unannotated_rows(expected_cells(&[
        ("empty_events", ["SpliceVault_top_events=", "", "", "{}"]),
        (
            "empty_event_part",
            [
                "SpliceVault_top_events=Top1:CA:-480:7%:,Top2:CA:-76:39%:Frameshift",
                "Top1:CA:-480:7%:,Top2:CA:-76:39%:Frameshift",
                "Top1:CA:-480:7%:&Top2:CA:-76:39%:Frameshift",
                &trailing_empty_part,
            ],
        ),
        (
            "comma_in_part",
            [
                "SpliceVault_top_events=Top1:CD:-330:73_314%:inFrame,Top2:CA:+112:1,5%:inFrame",
                "Top1:CD:-330:73_314%:inFrame,Top2:CA:+112:1,5%:inFrame",
                "Top1:CD:-330:73_314%:inFrame&Top2:CA:+112:1&5%:inFrame",
                &comma_in_part,
            ],
        ),
    ]));
    for (format, cells) in events.iter().enumerate() {
        for id in ["empty_events", "empty_event_part", "comma_in_part"] {
            assert_eq!(cells[id], expected[format][id], "{id}, format {format}");
        }
    }
}

/// `tests/fixtures/dosage-sensitivity-cells`: a two-cell line, a line whose
/// `pHaplo` cell is empty, a line whose `pTriplo` cell is empty and a one-cell
/// line, one corpus record on each gene. `DosageSensitivity.pm:138` assigns
/// the split line to three variables, so Ensembl VEP keeps an empty cell, trailing or
/// not, as the empty string and leaves a missing cell undefined, and `:141-142`
/// store both keys for the gene: a missing cell is `pTriplo=`, `-`, an empty
/// CSQ field and no JSON key, an empty cell `pHaplo=`, an empty tab cell, an
/// empty CSQ field and `""` in JSON (`VEP_output.pm:117-119`, `Tab.pm:128`,
/// `JSON.pm:327`), the rows Ensembl VEP 116.3 wrote on the fixture.
#[test]
fn dosage_sensitivity_short_line_and_empty_cell_in_every_format() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dosage-sensitivity-cells");
    let plugin = format!(
        "DosageSensitivity,file={}",
        fixture
            .join("dosage_sensitivity_cells_synthetic.tsv.gz")
            .display()
    );
    // The SOD1 and KRTAP20-4 records overlap transcripts of neighbouring genes
    // the table lacks, whose rows carry the absent-field cell.
    let unannotated = ["no item", "-", "", "no key"];
    let with_unannotated_rows = |mut expected: [BTreeMap<String, BTreeSet<String>>; 4]| {
        for (format, cell) in unannotated.into_iter().enumerate() {
            for id in ["sod1_exon2", "krtap20_4_snv"] {
                expected[format]
                    .get_mut(id)
                    .unwrap()
                    .insert(cell.to_string());
            }
        }
        expected
    };
    let cells = |field: &str| {
        field_cells_in_every_format(
            &fixture,
            &corpus("GRCh38-plugin-dosage-sensitivity"),
            "GRCh38",
            &plugin,
            field,
        )
    };
    assert_eq!(
        cells("pHaplo"),
        with_unannotated_rows(expected_cells(&[
            ("app_exon", ["pHaplo=0.5", "0.5", "0.5", "0.5"]),
            ("bach1_exon", ["pHaplo=", "", "", "\"\""]),
            ("sod1_exon2", ["pHaplo=0.1", "0.1", "0.1", "0.1"]),
            ("krtap20_4_snv", ["pHaplo=", "-", "", "no key"]),
        ]))
    );
    assert_eq!(
        cells("pTriplo"),
        with_unannotated_rows(expected_cells(&[
            ("app_exon", ["pTriplo=", "-", "", "no key"]),
            ("bach1_exon", ["pTriplo=0.3", "0.3", "0.3", "0.3"]),
            ("sod1_exon2", ["pTriplo=", "", "", "\"\""]),
            ("krtap20_4_snv", ["pTriplo=", "-", "", "no key"]),
        ]))
    );
}

/// The corpus's made-up GRCh38 file, prepared as the Ensembl VEP plugin documents
/// (dbscSNV.pm:53-55: sorted and indexed on columns 5 and 6, its `hg38_chr` and
/// `hg38_pos`), is read under
/// `--assembly GRCh38`: `ada_score` and `rf_score` print on the transcript rows
/// of a single-base substitution, with Ensembl VEP's header descriptions, and
/// an insertion or deletion gets nothing. The positional assembly override that
/// names the other assembly is refused under `--safe`, naming the indexed column.
#[test]
fn dbscsnv_grch38_prepared_as_documented() {
    let dir = corpus("GRCh38-plugin-dbscsnv-fixture");
    let file = dir.join("plugin_data/dbscSNV1.1_synthetic_GRCh38.txt.gz");
    let tmp = tempfile::tempdir().unwrap();
    let plugin = format!("dbscSNV,{}", file.display());
    let (output, body) = run(&dir, "GRCh38", tmp.path(), &["--plugin", &plugin]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        body.contains("## ada_score : dbscSNV ADA score\n"),
        "{body}"
    );
    assert!(body.contains("## rf_score : dbscSNV RF score\n"));
    let donor = rows_of(&body, "splice_donor");
    assert!(!donor.is_empty());
    for row in &donor {
        let ada = extra_field(row, "ada_score").unwrap_or_else(|| panic!("{row}"));
        assert!(ada.parse::<f64>().is_ok(), "{ada}");
    }
    let deletion = rows_of(&body, "deletion_rev");
    assert!(!deletion.is_empty());
    assert!(
        deletion.iter().all(|row| !row.contains("_score=")),
        "{deletion:?}"
    );

    let other = format!("dbscSNV,{},GRCh37", file.display());
    let (output, _) = run(&dir, "GRCh38", tmp.path(), &["--safe", "--plugin", &other]);
    assert_eq!(output.status.code(), Some(1));
    let stderr = text(&output.stderr);
    assert!(
        stderr.contains("a GRCh37 run reads 'pos' (column 2)")
            && stderr.contains("is indexed on 'hg38_pos' (column 6)"),
        "{stderr}"
    );
}

/// `dbscsnv_file_without_gene_column`: the GRCh37 corpus's made-up file cut to
/// the six columns the Ensembl VEP module reads (`chr pos ref alt ada_score
/// rf_score`), indexed as the module documents (dbscSNV.pm:58-59), is read.
/// Ensembl VEP's module dies on a header without `Ensembl_gene`
/// (dbscSNV.pm:111-113), a column whose one test (dbscSNV.pm:169-170) is
/// commented out. Every row equals the row Ensembl VEP writes from the whole
/// file, the corpus's stored reference of its `dbscsnv` flag set, under
/// `--safe` too, and no warning names the plugin.
#[test]
fn dbscsnv_reads_a_file_whose_header_lacks_the_gene_column() {
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;

    let dir = corpus("GRCh37-plugin-dbscsnv-fixture");
    let whole = dir.join("plugin_data/dbscSNV1.1_synthetic_GRCh37.txt.gz");
    let mut reader = noodles::bgzf::io::Reader::new(std::fs::File::open(whole).unwrap());
    let mut table = String::new();
    std::io::Read::read_to_string(&mut reader, &mut table).unwrap();
    // The fixture's lines end CRLF, and the cut keeps each line's end.
    let cut: Vec<String> = table
        .split('\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 18, "{line}");
            [cells[0], cells[1], cells[2], cells[3], cells[16], cells[17]].join("\t")
        })
        .collect();
    assert_eq!(cut[0], "chr\tpos\tref\talt\tada_score\trf_score\r");
    assert_eq!(cut.len(), 19);
    let records: Vec<&str> = cut[1..].iter().map(String::as_str).collect();
    let tmp = tempfile::tempdir().unwrap();
    let six = common::indexed_lines(
        tmp.path(),
        "dbscSNV1.1_synthetic_six_columns_GRCh37.txt.gz",
        HeaderBuilder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'c')
            .build(),
        &[cut[0].as_str()],
        &records,
        |fields| {
            let position = fields[1].parse().unwrap();
            (position, position)
        },
    );

    let corpus_record = corpora()
        .into_iter()
        .find(|corpus| corpus.name == "GRCh37-plugin-dbscsnv-fixture")
        .unwrap();
    let reference = common::runs(&corpus_record)
        .into_iter()
        .find(|run| run.name == "dbscsnv")
        .unwrap()
        .expected("default.txt");
    let data_rows = |output: &str| -> Vec<String> {
        output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect()
    };
    let reference_rows = data_rows(&reference);
    assert_eq!(reference_rows.len(), 67);
    assert_eq!(
        reference_rows
            .iter()
            .filter(|row| row.contains("_score="))
            .count(),
        57
    );

    let plugin = format!("dbscSNV,{}", six.display());
    let plugin = plugin.as_str();
    for extra in [&["--plugin", plugin][..], &["--safe", "--plugin", plugin]] {
        let (output, body) = run(&dir, "GRCh37", tmp.path(), extra);
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{extra:?}: {stderr}");
        assert!(
            !stderr.contains("WARNING") && !tmp.path().join("out.txt_warnings.txt").exists(),
            "{extra:?}: {stderr}"
        );
        assert!(
            body.contains("## ada_score : dbscSNV ADA score\n")
                && body.contains("## rf_score : dbscSNV RF score\n"),
            "{extra:?}: {body}"
        );
        assert_eq!(data_rows(&body), reference_rows, "{extra:?}");
        let donor: Vec<&str> = rows_of(&body, "splice_donor_a")
            .into_iter()
            .filter(|row| row.contains("\tENST00000298232\t"))
            .collect();
        assert_eq!(donor.len(), 1, "{extra:?}: {donor:?}");
        assert_eq!(extra_field(donor[0], "ada_score"), Some("0.644219"));
        assert_eq!(extra_field(donor[0], "rf_score"), Some("0.058437"));
    }
}

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/plugins")
}

/// Runs `vep` on the GRCh38 golden corpus in the default output format without
/// `--quiet`, so warnings reach stderr, with `--symbol`, so each row names the
/// gene symbol the plugins key on (they read it from the cache whether or not
/// it prints).
fn vep_on_corpus(out: &Path, extra: &[&str]) -> Output {
    let corpus = corpus("GRCh38");
    vep()
        .arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--symbol"])
        .args(extra)
        .output()
        .expect("vep binary runs")
}

fn lines(out: &Path, header: bool) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("##") == header)
        .map(str::to_owned)
        .collect()
}

/// The value of `key` in a default-format row's Extra column.
fn extra(row: &str, key: &str) -> Option<String> {
    row.split('\t')
        .nth(13)?
        .split(';')
        .find_map(|pair| pair.strip_prefix(key)?.strip_prefix('=').map(str::to_owned))
}

/// Whether a default-format row is a structural variant's: its Allele column
/// carries the class word (`deletion`, `Alu_insertion`) or a breakend form
/// (`N.`, `N[chr21:10601732[`) where a sequence variant's carries bases or `-`;
/// the unspecified allele `<*>` alone is a sequence allele at POS
/// (Parser/VCF.pm:237-245). Ensembl VEP hands a plugin with the default
/// `variant_feature_types` no structural variant (BaseVepPlugin.pm:100,
/// OutputFactory.pm:2229), so such a row carries no gene-score field.
fn structural_row(row: &str) -> bool {
    let allele = row.split('\t').nth(2).unwrap_or("");
    allele != "<*>"
        && !allele
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N' | b'-'))
}

/// Two whitespace-separated columns of a fixture, keyed lower-cased, as the
/// plugins file them.
fn fixture_table(name: &str) -> HashMap<String, f64> {
    std::fs::read_to_string(fixtures_dir().join(name))
        .unwrap()
        .lines()
        .filter_map(|line| {
            let mut tokens = line.split_ascii_whitespace();
            let key = tokens.next()?.to_ascii_lowercase();
            Some((key, tokens.next()?.parse().ok()?))
        })
        .collect()
}

fn success(output: &Output) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(!stderr.contains("Failed to"), "{stderr}");
}

/// `--plugin LoFtool` with no file reads `LoFtool_scores.txt` from the plugins
/// directory (LoFtool.pm:74-78) and writes the rows the same file named as a
/// parameter does: the field `LoFtool`, its header description, and each
/// scored symbol's percentile as Perl's `%g` prints it (`:89`, `:107`, `:118`).
#[test]
fn loftool_without_path_uses_plugins_dir() {
    let tmp = tempfile::tempdir().unwrap();
    let by_dir = tmp.path().join("by_dir.txt");
    let fixtures = fixtures_dir();
    success(&vep_on_corpus(
        &by_dir,
        &[
            "--dir_plugins",
            fixtures.to_str().unwrap(),
            "--plugin",
            "LoFtool",
        ],
    ));
    assert!(
        lines(&by_dir, true).contains(&"## LoFtool : LoFtool score for gene".to_string()),
        "{:?}",
        lines(&by_dir, true)
    );

    let scores = fixture_table("LoFtool_scores.txt");
    let rows = lines(&by_dir, false);
    let mut scored = 0;
    let mut structural_scorable = 0;
    for row in &rows {
        let symbol = extra(row, "SYMBOL");
        let scorable = symbol
            .as_deref()
            .and_then(|symbol| scores.get(&symbol.to_ascii_lowercase()))
            .map(|score| vep_builtins::perl_format::sprintf_g(*score));
        let expected = if structural_row(row) {
            structural_scorable += usize::from(scorable.is_some());
            None
        } else {
            scorable
        };
        assert_eq!(extra(row, "LoFtool"), expected, "{row}");
        scored += usize::from(expected.is_some());
    }
    assert!(scored > 100, "{scored} scored rows of {}", rows.len());
    assert!(
        structural_scorable > 0,
        "no structural row of a scored gene"
    );

    let by_path = tmp.path().join("by_path.txt");
    let plugin = format!("LoFtool,{}", fixtures.join("LoFtool_scores.txt").display());
    success(&vep_on_corpus(
        &by_path,
        &[
            "--dir_plugins",
            fixtures.to_str().unwrap(),
            "--plugin",
            &plugin,
        ],
    ));
    assert_eq!(lines(&by_path, false), rows);
}

/// pLI's `transcript` mode files the fixture's scores under transcript stable
/// ids and writes `pLI_transcript_value` on each row whose Feature has one,
/// matched case-insensitively (pLI.pm:174-185, 224-228); `sigfig` prints three
/// significant figures (`:105-110`), and the default gene mode without a file
/// reads `pLI_values.txt` from the plugins directory under two decimals
/// (`:115-119`).
#[test]
fn pli_transcript_mode_and_sigfig() {
    let tmp = tempfile::tempdir().unwrap();
    let fixtures = fixtures_dir();
    let transcript = tmp.path().join("transcript.txt");
    let plugin = format!(
        "pLI,{},transcript,sigfig",
        fixtures.join("pLI_transcript_values.txt").display()
    );
    success(&vep_on_corpus(
        &transcript,
        &[
            "--dir_plugins",
            fixtures.to_str().unwrap(),
            "--plugin",
            &plugin,
        ],
    ));
    assert!(
        lines(&transcript, true)
            .contains(&"## pLI_transcript_value : pLI value by transcript".to_string()),
        "{:?}",
        lines(&transcript, true)
    );
    let scores = fixture_table("pLI_transcript_values.txt");
    let rows = lines(&transcript, false);
    let mut scored = 0;
    let mut structural_scorable = 0;
    for row in &rows {
        let feature = row.split('\t').nth(4).unwrap();
        let scorable = scores
            .get(&feature.to_ascii_lowercase())
            .map(|score| vep_builtins::perl_format::sig3(*score));
        let expected = if structural_row(row) {
            structural_scorable += usize::from(scorable.is_some());
            None
        } else {
            scorable
        };
        assert_eq!(extra(row, "pLI_transcript_value"), expected, "{row}");
        assert_eq!(extra(row, "pLI_gene_value"), None, "{row}");
        scored += usize::from(expected.is_some());
    }
    assert!(scored > 10, "{scored} scored rows of {}", rows.len());
    assert!(
        structural_scorable > 0,
        "no structural row of a scored transcript"
    );

    let gene = tmp.path().join("gene.txt");
    success(&vep_on_corpus(
        &gene,
        &[
            "--dir_plugins",
            fixtures.to_str().unwrap(),
            "--plugin",
            "pLI",
        ],
    ));
    assert!(
        lines(&gene, true).contains(&"## pLI_gene_value : pLI value by gene".to_string()),
        "{:?}",
        lines(&gene, true)
    );
    let scores = fixture_table("pLI_values.txt");
    let mut scored = 0;
    for row in lines(&gene, false) {
        let expected = extra(&row, "SYMBOL")
            .filter(|_| !structural_row(&row))
            .and_then(|symbol| scores.get(&symbol.to_ascii_lowercase()).copied())
            .map(|score| vep_builtins::perl_format::fixed(score, 2));
        assert_eq!(extra(&row, "pLI_gene_value"), expected, "{row}");
        scored += usize::from(expected.is_some());
    }
    assert!(scored > 100, "{scored} scored rows");
}

/// Runs `vep` over `corpus` with `extra`, the output at `out`.
fn run_plugin(corpus: &Path, out: &Path, extra: &[&str]) -> Output {
    let mut cmd: Command = vep();
    cmd.args(args_for(corpus, "GRCh38"))
        .arg("-o")
        .arg(out)
        .args(extra)
        .current_dir(corpus)
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

/// Data lines of an output, the header dropped.
fn data_lines(output: &str) -> Vec<&str> {
    output.lines().filter(|l| !l.starts_with('#')).collect()
}

/// The rows of `output` whose Location column is `21:<pos>`.
fn rows_at(output: &str, pos: u64) -> Vec<&str> {
    let location = format!("21:{pos}");
    data_lines(output)
        .into_iter()
        .filter(|row| row.split('\t').nth(1) == Some(location.as_str()))
        .collect()
}

/// GWAS reads the curated catalog and a harmonised summary-statistics file
/// under its `key=value` grammar (GWAS.pm:94-104 at VEP_plugins release/116),
/// matches the risk allele against the record's alleles (`:182-194`) and
/// writes its eight fields on every row of a matching record, an empty field
/// included; the file is read compressed or not, and a file given by position
/// is refused.
#[test]
fn gwas_curated_and_sstate() {
    let corpus = corpus("GRCh38-plugin-gwas");
    let dir = tempfile::tempdir().unwrap();
    // `args_for` names the plugins directory; `--dir` names the directory of
    // the processed copies Ensembl VEP's plugin reads, inert for vep-rs.
    let plugin_dirs = ["--dir", "plugin_data/vepdir"];

    let out = dir.path().join("curated.txt");
    let run = run_plugin(
        &corpus,
        &out,
        &[
            &plugin_dirs[..],
            &["--plugin", "GWAS,file=plugin_data/gwas_catalog_slice.tsv"],
        ]
        .concat(),
    );
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let curated = std::fs::read_to_string(&out).unwrap();
    assert!(
        curated.contains("\n## GWAS_risk_allele : Allele associated with the variant that is strongly associated with the trait\n"),
        "no GWAS_risk_allele header line"
    );
    // The paired site of a 25-SNP haplotype row carries the row's association
    // with the SNP's own risk allele.
    let haplotype = rows_at(&curated, 34300618);
    assert!(!haplotype.is_empty());
    for row in &haplotype {
        assert!(
            row.contains(";GWAS_associated_gene=KCNE2,RCAN1;GWAS_beta_coef=;GWAS_odds_ratio=0.41;GWAS_p_value=5E-7;GWAS_pmid=29299148;GWAS_risk_allele=G;GWAS_study=GCST005275"),
            "{row}"
        );
    }
    // A multi-allelic record whose second allele is the risk allele is
    // annotated on the rows of both alleles, as the plugin reads every
    // alternate of the record.
    let multi_allelic = rows_at(&curated, 14992435);
    for allele in ["C", "G"] {
        assert!(
            multi_allelic
                .iter()
                .any(|row| row.split('\t').nth(2) == Some(allele)
                    && row.contains(";GWAS_risk_allele=G;")),
            "allele {allele} of the multi-allelic record carries no GWAS fields"
        );
    }
    // The complement of the risk allele is another allele.
    let complement = rows_at(&curated, 27428688);
    assert!(!complement.is_empty());
    assert!(complement.iter().all(|row| !row.contains("GWAS_")));

    let out = dir.path().join("curated.tab");
    let run = run_plugin(
        &corpus,
        &out,
        &[
            &plugin_dirs[..],
            &[
                "--plugin",
                "GWAS,file=plugin_data/gwas_catalog_slice.tsv",
                "--tab",
            ],
        ]
        .concat(),
    );
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let tab = std::fs::read_to_string(&out).unwrap();
    let header: Vec<&str> = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap()
        .split('\t')
        .collect();
    let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let row: Vec<&str> = rows_at(&tab, 34300618)[0].split('\t').collect();
    // An empty field is an empty cell, Ensembl VEP's empty string, not `-`.
    assert_eq!(row[column("GWAS_beta_coef")], "");
    assert_eq!(row[column("GWAS_odds_ratio")], "0.41");

    let out = dir.path().join("curated.json");
    let run = run_plugin(
        &corpus,
        &out,
        &[
            &plugin_dirs[..],
            &[
                "--plugin",
                "GWAS,file=plugin_data/gwas_catalog_slice.tsv",
                "--json",
            ],
        ]
        .concat(),
    );
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let json = std::fs::read_to_string(&out).unwrap();
    assert!(
        json.contains("\"gwas_beta_coef\":\"\"") && json.contains("\"gwas_odds_ratio\":0.41"),
        "the JSON format keeps an empty field as the empty string"
    );

    // The summary-statistics file, uncompressed in the corpus and compressed
    // beside it, reads the same.
    let gz = dir.path().join("17463246-GCST000028-EFO_0001360.h.tsv.gz");
    {
        let source = corpus.join("plugin_data/17463246-GCST000028-EFO_0001360.h.tsv");
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&gz).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(&std::fs::read(source).unwrap()).unwrap();
        encoder.finish().unwrap();
    }
    let mut outputs = Vec::new();
    for file in [
        "plugin_data/17463246-GCST000028-EFO_0001360.h.tsv".to_string(),
        gz.display().to_string(),
    ] {
        let out = dir.path().join("sstate.txt");
        let run = run_plugin(
            &corpus,
            &out,
            &[
                &plugin_dirs[..],
                &["--plugin", &format!("GWAS,type=sstate,file={file}")],
            ]
            .concat(),
        );
        assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
        let output = std::fs::read_to_string(&out).unwrap();
        let rows = rows_at(&output, 34036428);
        assert!(!rows.is_empty());
        for row in &rows {
            assert!(
                row.contains(";GWAS_accessions=EFO_0001360;GWAS_associated_gene=;GWAS_beta_coef=NA;GWAS_odds_ratio=NA;GWAS_p_value=0.1965;GWAS_pmid=17463246;GWAS_risk_allele=G;GWAS_study=GCST000028"),
                "{row}"
            );
        }
        outputs.push(data_lines(&output).join("\n"));
    }
    assert_eq!(outputs[0], outputs[1]);

    let out = dir.path().join("refused.txt");
    let refused = run_plugin(
        &corpus,
        &out,
        &[
            &plugin_dirs[..],
            &[
                "--plugin",
                "GWAS,plugin_data/gwas_catalog_slice.tsv",
                "--safe",
            ],
        ]
        .concat(),
    );
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr).lines().any(|l| l
            == "ERROR: [plugin] Failed to instantiate plugin GWAS: no file specified; use file=<path>"),
        "{}",
        text(&refused.stderr)
    );
    assert!(!out.exists());
}

/// The data lines of `output` carrying `key=` in their Extra column.
fn rows_with(output: &str, key: &str) -> usize {
    let needle = format!("{key}=");
    data_lines(output)
        .iter()
        .filter(|row| {
            row.rsplit('\t')
                .next()
                .is_some_and(|extra| extra.split(';').any(|kv| kv.starts_with(&needle)))
        })
        .count()
}

/// The `## <field> : <text>` header lines of a default-format output whose
/// field name opens with `prefix`: a plugin's lines among the columns'.
fn description_lines<'a>(output: &'a str, prefix: &str) -> Vec<&'a str> {
    let opening = format!("## {prefix}");
    output
        .lines()
        .filter(|l| l.starts_with(&opening) && l.contains(" : "))
        .collect()
}

/// The output's lines but the ones that differ between two runs of one record
/// set by construction: the time and the command line.
fn comparable(output: &str) -> String {
    output
        .lines()
        .filter(|l| {
            !l.starts_with("## Output produced at ") && !l.starts_with("## VEP command-line:")
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// `REVEL.pm:105-110` reads the file and `no_match` by position or by key:
/// `--plugin REVEL,<file>,1` is `--plugin REVEL,file=<file>,no_match=1`, and
/// either writes the score on rows the default run leaves unscored because the
/// record's transcripts do not name theirs (`:191-197`). The one field is
/// `REVEL` (`:156`), and the two grammars do not mix.
#[test]
fn revel_positional_no_match_flag() {
    let corpus = corpus("GRCh38-plugin-revel");
    let file = corpus.join("plugin_data/revel.tsv.gz");
    let file = file.to_str().unwrap();
    let dir = tempfile::tempdir().unwrap();

    let out = dir.path().join("default.txt");
    let output = run_plugin(&corpus, &out, &["--plugin", &format!("REVEL,file={file}")]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let default_text = std::fs::read_to_string(&out).unwrap();
    assert_eq!(
        description_lines(&default_text, "REVEL"),
        ["## REVEL : Rare Exome Variant Ensemble Learner"]
    );
    let default_scored = rows_with(&default_text, "REVEL");
    assert!(default_scored > 0, "no row carries a REVEL score");

    let mut no_match_texts = Vec::new();
    for (name, plugin) in [
        ("positional", format!("REVEL,{file},1")),
        ("keyed", format!("REVEL,file={file},no_match=1")),
    ] {
        let out = dir.path().join(format!("{name}.txt"));
        let output = run_plugin(&corpus, &out, &["--plugin", &plugin]);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{name}: {}",
            text(&output.stderr)
        );
        let no_match_text = std::fs::read_to_string(&out).unwrap();
        assert!(
            rows_with(&no_match_text, "REVEL") > default_scored,
            "{name}: no_match must score the rows whose transcript the record does not name"
        );
        no_match_texts.push(comparable(&no_match_text));
    }
    assert_eq!(no_match_texts[0], no_match_texts[1]);

    // A plugin that fails to load is left out of the run with Ensembl VEP's
    // warning; `--safe` makes the failure the run's.
    let out = dir.path().join("mixed.txt");
    let mixed = run_plugin(
        &corpus,
        &out,
        &["--plugin", &format!("REVEL,file={file},1"), "--safe"],
    );
    assert_eq!(mixed.status.code(), Some(1));
    assert!(
        text(&mixed.stderr).lines().any(|l| l
            == format!(
                "ERROR: [plugin] Failed to instantiate plugin REVEL: parameter '1' is \
                 not key=value while 'file={file}' is; give every parameter as key=value, or \
                 every one by position"
            )),
        "{}",
        text(&mixed.stderr)
    );
    assert!(!out.exists());
}

/// A tabix-indexed REVEL table named `name` under `dir` holding `header` and
/// `rows`, every row on chromosome 21 at `position`, indexed on the 1-based
/// column `position_column` as the table of that assembly is.
fn revel_tabix_tsv(
    dir: &Path,
    name: &str,
    header: &str,
    rows: &[String],
    position: usize,
    position_column: usize,
) -> String {
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;

    let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
    common::indexed_lines(
        dir,
        name,
        HeaderBuilder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(position_column - 1)
            .set_end_position_index(None)
            .build(),
        &[header],
        &rows,
        |_| (position, position),
    )
    .display()
    .to_string()
}

/// A REVEL record scores the substitution from its `aaref` to its `aaalt`
/// (`REVEL.pm:210` and `:229` read both cells; `:193` and `:197` test the
/// second alone): the row whose transcript a record names beside its alternate
/// amino acid and another reference amino acid takes no `REVEL`, in the
/// nine-column layout on either assembly's position column and in the
/// seven-column layout, while the row of the record's own pair keeps its
/// score; under `no_match` every row takes the first record of its own pair.
/// On a record of two alternates each alternate's row is put to its own
/// records.
#[test]
fn revel_requires_the_reference_amino_acid() {
    const NINE: &str =
        "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid";
    const SEVEN: &str = "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL";
    let dir = tempfile::tempdir().unwrap();
    let revel = |output: &str, record: &str, allele: &str, transcript: &str| -> Option<String> {
        let rows: Vec<&str> = rows_of(output, record)
            .into_iter()
            .filter(|row| {
                let columns: Vec<&str> = row.split('\t').collect();
                columns[2] == allele && columns[4] == transcript
            })
            .collect();
        assert_eq!(rows.len(), 1, "{record} {allele} {transcript}: {output}");
        assert!(rows[0].contains("missense_variant"), "{}", rows[0]);
        extra(rows[0], "REVEL")
    };
    let run_on = |corpus: &Path, assembly: &str, input: &Path, plugin: &str, tag: &str| {
        let out = dir.path().join(format!("{tag}.txt"));
        let output = common::isolated_vep(dir.path())
            .arg("-i")
            .arg(input)
            .arg("-o")
            .arg(&out)
            .arg("--offline")
            .arg("--json_cache")
            .arg(corpus.join("json_cache"))
            .args(["--species", "homo_sapiens", "--assembly", assembly])
            .args(["--force_overwrite", "--no_stats"])
            .args(["--plugin", plugin])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{tag}: {}",
            text(&output.stderr)
        );
        std::fs::read_to_string(&out).unwrap()
    };

    // GRCh37: `21 34640722 G A` is `E/K` on both transcripts.
    let grch37 = corpus("GRCh37-plugin-revel");
    let (named, other) = ("ENST00000290200", "ENST00000422891");
    let nine = revel_tabix_tsv(
        dir.path(),
        "nine37.tsv.gz",
        NINE,
        &[
            format!("21\t34640722\t999\tG\tA\tD\tK\t0.901\t{named}"),
            format!("21\t34640722\t999\tG\tA\tE\tK\t0.402\t{other}"),
        ],
        34640722,
        2,
    );
    let seven = revel_tabix_tsv(
        dir.path(),
        "seven37.tsv.gz",
        SEVEN,
        &[
            "21\t34640722\tG\tA\tD\tK\t0.901".to_string(),
            "21\t34640722\tG\tA\tE\tK\t0.402".to_string(),
        ],
        34640722,
        2,
    );
    std::fs::create_dir(dir.path().join("one37")).unwrap();
    let one = one_record_vcf(
        &dir.path().join("one37"),
        "21\t34640722\tone_alt\tG\tA\t.\t.\t.",
    );
    let body = run_on(
        &grch37,
        "GRCh37",
        &one,
        &format!("REVEL,file={nine}"),
        "nine37",
    );
    assert_eq!(revel(&body, "one_alt", "A", named), None);
    assert_eq!(
        revel(&body, "one_alt", "A", other).as_deref(),
        Some("0.402")
    );
    let body = run_on(
        &grch37,
        "GRCh37",
        &one,
        &format!("REVEL,file={nine},no_match=1"),
        "nine37-no-match",
    );
    assert_eq!(
        revel(&body, "one_alt", "A", named).as_deref(),
        Some("0.402")
    );
    assert_eq!(
        revel(&body, "one_alt", "A", other).as_deref(),
        Some("0.402")
    );
    let body = run_on(
        &grch37,
        "GRCh37",
        &one,
        &format!("REVEL,file={seven}"),
        "seven37",
    );
    assert_eq!(
        revel(&body, "one_alt", "A", named).as_deref(),
        Some("0.402")
    );
    assert_eq!(
        revel(&body, "one_alt", "A", other).as_deref(),
        Some("0.402")
    );

    // The second alternate, `G C`, is `E/Q` on the same transcripts.
    let two_rows = revel_tabix_tsv(
        dir.path(),
        "two37.tsv.gz",
        NINE,
        &[
            format!("21\t34640722\t999\tG\tA\tE\tK\t0.277\t{named}"),
            format!("21\t34640722\t999\tG\tC\tD\tQ\t0.902\t{named}"),
        ],
        34640722,
        2,
    );
    std::fs::create_dir(dir.path().join("two37")).unwrap();
    let two = one_record_vcf(
        &dir.path().join("two37"),
        "21\t34640722\ttwo_alts\tG\tA,C\t.\t.\t.",
    );
    let body = run_on(
        &grch37,
        "GRCh37",
        &two,
        &format!("REVEL,file={two_rows}"),
        "two37",
    );
    assert_eq!(
        revel(&body, "two_alts", "A", named).as_deref(),
        Some("0.277")
    );
    assert_eq!(revel(&body, "two_alts", "C", named), None);

    // GRCh38: `21 26965413 T C` is `S/G` on the transcript.
    let grch38 = corpus("GRCh38-plugin-revel");
    let named = "ENST00000284987";
    let nine = revel_tabix_tsv(
        dir.path(),
        "nine38.tsv.gz",
        NINE,
        &[
            format!("21\t999\t26965413\tT\tC\tN\tG\t0.901\t{named}"),
            "21\t999\t26965413\tT\tC\tS\tG\t0.402\tENST00000970346".to_string(),
        ],
        26965413,
        3,
    );
    std::fs::create_dir(dir.path().join("one38")).unwrap();
    let one = one_record_vcf(
        &dir.path().join("one38"),
        "21\t26965413\tone_alt\tT\tC\t.\t.\t.",
    );
    let body = run_on(
        &grch38,
        "GRCh38",
        &one,
        &format!("REVEL,file={nine}"),
        "nine38",
    );
    assert_eq!(revel(&body, "one_alt", "C", named), None);
    assert_eq!(
        revel(&body, "one_alt", "C", "ENST00000970346").as_deref(),
        Some("0.402")
    );
    let body = run_on(
        &grch38,
        "GRCh38",
        &one,
        &format!("REVEL,file={nine},no_match=1"),
        "nine38-no-match",
    );
    assert_eq!(
        revel(&body, "one_alt", "C", named).as_deref(),
        Some("0.402")
    );
}

/// `AlphaMissense.pm` reads `key=value` pairs alone (`:175-181`): `cols=all`
/// writes every column of the file under its `am_` name (`:146-147`), a column
/// the header lacks is refused with the valid names (`:161-163`),
/// `transcript_match=1` keeps a prediction on the row of its own transcript
/// alone (`:256-261`), a file given by position is the file given without its
/// key, and the JSON format nests the fields under `alphamissense`
/// (`:281-282`).
#[test]
fn alphamissense_cols_and_transcript_match() {
    let corpus = corpus("GRCh38-plugin-alphamissense");
    let file = corpus.join("plugin_data/alphamissense.tsv.gz");
    let file = file.to_str().unwrap();
    let dir = tempfile::tempdir().unwrap();

    let out = dir.path().join("default.txt");
    let output = run_plugin(
        &corpus,
        &out,
        &["--plugin", &format!("AlphaMissense,file={file}")],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let default_text = std::fs::read_to_string(&out).unwrap();
    let names: Vec<&str> = description_lines(&default_text, "am_")
        .iter()
        .map(|l| l.split(" : ").next().unwrap().trim_start_matches("## "))
        .collect();
    assert_eq!(names, ["am_class", "am_pathogenicity"]);
    let default_rows = rows_with(&default_text, "am_class");
    assert!(default_rows > 0, "no row carries a prediction");
    assert_eq!(rows_with(&default_text, "am_genome"), 0);

    let out = dir.path().join("all.txt");
    let output = run_plugin(
        &corpus,
        &out,
        &["--plugin", &format!("AlphaMissense,file={file},cols=all")],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let all_text = std::fs::read_to_string(&out).unwrap();
    let names: Vec<&str> = description_lines(&all_text, "am_")
        .iter()
        .map(|l| l.split(" : ").next().unwrap().trim_start_matches("## "))
        .collect();
    assert_eq!(
        names,
        [
            "am_class",
            "am_genome",
            "am_pathogenicity",
            "am_protein_variant",
            "am_transcript_id",
            "am_uniprot_id"
        ]
    );
    assert!(all_text.contains(&format!("\n## am_genome : column from {file}\n")));
    assert_eq!(rows_with(&all_text, "am_genome"), default_rows);
    assert_eq!(rows_with(&all_text, "am_class"), default_rows);

    let out = dir.path().join("matched.txt");
    let output = run_plugin(
        &corpus,
        &out,
        &[
            "--plugin",
            &format!("AlphaMissense,file={file},transcript_match=1"),
        ],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let matched_text = std::fs::read_to_string(&out).unwrap();
    let matched_rows = rows_with(&matched_text, "am_class");
    assert!(
        matched_rows > 0 && matched_rows < default_rows,
        "transcript_match=1 keeps the rows of the predicted transcripts alone: {matched_rows} of {default_rows}"
    );

    let out = dir.path().join("json.jsonl");
    let output = run_plugin(
        &corpus,
        &out,
        &["--json", "--plugin", &format!("AlphaMissense,file={file}")],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let mut nested = 0usize;
    for line in std::fs::read_to_string(&out).unwrap().lines() {
        let record: serde_json::Value = serde_json::from_str(line).unwrap();
        for row in record["transcript_consequences"]
            .as_array()
            .map(Vec::as_slice)
            .unwrap_or_default()
        {
            let keys = row.as_object().unwrap();
            assert!(
                !keys.keys().any(|k| k.starts_with("am_")),
                "a flat am_ key in {row}"
            );
            if let Some(group) = row.get("alphamissense") {
                assert!(group["am_pathogenicity"].is_number(), "{group}");
                assert!(group["am_class"].is_string(), "{group}");
                nested += 1;
            }
        }
    }
    assert_eq!(nested, default_rows);

    for (params, message) in [
        (
            format!("AlphaMissense,file={file},cols=am_class:score"),
            "the following columns were not found in file header: am_score; valid columns \
             are: am_genome, am_uniprot_id, am_transcript_id, am_protein_variant, \
             am_pathogenicity, am_class"
                .to_string(),
        ),
        (
            format!("AlphaMissense,{file}"),
            "no file specified; use file=<path>".to_string(),
        ),
    ] {
        let out = dir.path().join("refused.txt");
        let refused = run_plugin(&corpus, &out, &["--plugin", &params, "--safe"]);
        assert_eq!(refused.status.code(), Some(1), "{params}");
        let expected =
            format!("ERROR: [plugin] Failed to instantiate plugin AlphaMissense: {message}");
        assert!(
            text(&refused.stderr).lines().any(|l| l == expected),
            "{params}: {}",
            text(&refused.stderr)
        );
        assert!(!out.exists());
    }
}

/// Runs `vep` over `corpus` as `run_plugin` does, with `--symbol` and the
/// corpus's `reference.fa.gz` as `--fasta`.
fn run_corpus(corpus: &Path, extra: &[&str], out: &Path) -> Output {
    let fasta = corpus.join("reference.fa.gz");
    let mut args = vec!["--fasta", fasta.to_str().unwrap(), "--symbol"];
    args.extend_from_slice(extra);
    run_plugin(corpus, out, &args)
}

/// The keys of a default-format row's Extra column, in the column's order.
fn extra_keys(row: &str) -> Vec<&str> {
    row.split('\t')
        .nth(13)
        .unwrap_or("")
        .split(';')
        .filter_map(|pair| pair.split_once('='))
        .map(|(key, _)| key)
        .collect()
}

/// SpliceAI's fields take Ensembl VEP's shape per format (SpliceAI.pm:173-179,
/// 208-224, 275-292, 316 at VEP_plugins release/116): the default and tab
/// formats print one `SpliceAI_pred` field, `SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|
/// DP_AG|DP_AL|DP_DG|DP_DL`, under one header line naming the file's tool
/// version; the VCF format and `split_output=1` print the nine
/// `SpliceAI_pred_*` fields under nine header lines; the JSON format nests the
/// nine parts under `spliceai` with their prefix dropped and the scores as
/// numbers; `cutoff=0.5` adds `SpliceAI_cutoff`, `PASS` or `FAIL` by the largest
/// delta score of the row's entries, and keeps every record; a row of a gene the file names no
/// record for carries no field, and a record of a gene the transcript does not
/// belong to is not taken. The corpus's score files are made up, under the file
/// names of the plugin's usage example (SpliceAI.pm:127) with `synthetic` added
/// before the extension.
#[test]
fn spliceai_output_shapes_per_format() {
    let corpus = corpus("GRCh38-plugin-spliceai-fixture");
    let plugin = "SpliceAI,snv=plugin_data/spliceai_scores.raw.snv.hg38.synthetic.vcf.gz,\
                  indel=plugin_data/spliceai_scores.raw.indel.hg38.synthetic.vcf.gz";
    let dir = tempfile::tempdir().unwrap();

    let out = dir.path().join("default.txt");
    let run = run_corpus(&corpus, &["--plugin", plugin], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        output.contains(
            "\n## SpliceAI_pred : SpliceAI (v1.3) predicted effect on splicing. These include \
             delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss \
             (AL), donor gain (DG), and donor loss (DL). Format: \
             SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\n"
        ),
        "{output}"
    );
    assert!(!output.contains("## SpliceAI_pred_"), "{output}");
    let missense_rev = rows_of(&output, "snv_missense_rev");
    assert_eq!(missense_rev.len(), 1);
    let pred =
        extra(missense_rev[0], "SpliceAI_pred").expect("the missense row carries SpliceAI_pred");
    let parts: Vec<&str> = pred.split('|').collect();
    assert_eq!(parts.len(), 9, "{pred}");
    assert_eq!(parts[0], "TMPRSS15");
    assert!(
        parts[1..5].iter().all(|score| score.parse::<f64>().is_ok()),
        "{pred}"
    );
    assert!(
        parts[5..]
            .iter()
            .all(|position| position.parse::<i64>().is_ok()),
        "{pred}"
    );
    let keys = extra_keys(missense_rev[0]);
    assert!(
        !keys
            .iter()
            .any(|key| key.starts_with("SpliceAI_pred_") || *key == "SpliceAI_cutoff"),
        "{keys:?}"
    );
    // A symbol the file does not name, and a site the file does not score, take nothing.
    for name in ["snv_symbol_mismatch", "snv_unscored", "snv_intergenic"] {
        let bare = rows_of(&output, name);
        assert!(!bare.is_empty(), "{name} has rows");
        assert!(
            bare.iter().all(|row| !row.contains("SpliceAI")),
            "{name} carries a SpliceAI field: {bare:?}"
        );
    }
    // Two genes at one site: each transcript takes its own gene's record.
    let two_genes: Vec<(String, String)> = rows_of(&output, "snv_two_genes")
        .iter()
        .map(|row| {
            let symbol = extra(row, "SYMBOL").unwrap_or_else(|| "-".to_string());
            let pred_symbol = extra(row, "SpliceAI_pred")
                .map(|pred| pred.split('|').next().unwrap().to_string())
                .unwrap_or_else(|| "-".to_string());
            (symbol, pred_symbol)
        })
        .collect();
    assert!(two_genes.len() >= 2, "{two_genes:?}");
    assert!(
        two_genes
            .iter()
            .all(|(symbol, pred_symbol)| symbol == pred_symbol),
        "{two_genes:?}"
    );
    assert!(
        two_genes.iter().any(|(symbol, _)| symbol == "CXADR"),
        "{two_genes:?}"
    );
    assert!(
        two_genes.iter().any(|(symbol, _)| symbol == "BTG3-AS1"),
        "{two_genes:?}"
    );

    let out = dir.path().join("tab.txt");
    let run = run_corpus(&corpus, &["--plugin", plugin, "--tab"], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    let header = output
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap();
    let columns: Vec<&str> = header.split('\t').collect();
    assert!(columns.contains(&"SpliceAI_pred"), "{header}");
    assert!(
        !columns.iter().any(|c| c.starts_with("SpliceAI_pred_")),
        "{header}"
    );

    let out = dir.path().join("out.vcf");
    let run = run_corpus(&corpus, &["--plugin", plugin, "--vcf"], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    for part in [
        "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    ] {
        assert!(
            output.contains(&format!("\n##SpliceAI_pred_{part}=SpliceAI (v1.3) ")),
            "no header line for SpliceAI_pred_{part}"
        );
    }
    let csq = output
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ"))
        .unwrap();
    assert!(
        csq.contains(
            "|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL\
             |SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL\
             |SpliceAI_pred_SYMBOL"
        ),
        "{csq}"
    );
    assert!(!csq.contains("|SpliceAI_pred|"), "{csq}");
    let missense_rev = data_lines(&output)
        .into_iter()
        .find(|row| row.split('\t').nth(2) == Some("snv_missense_rev"))
        .unwrap_or_else(|| panic!("{output}"));
    assert!(
        missense_rev.ends_with("|-46|-17|-28|-43|0.71|1.00|0.55|0.26|TMPRSS15"),
        "{missense_rev}"
    );

    let out = dir.path().join("split.txt");
    let run = run_corpus(
        &corpus,
        &["--plugin", &format!("{plugin},split_output=1")],
        &out,
    );
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        output.contains("\n## SpliceAI_pred_SYMBOL : SpliceAI (v1.3) gene symbol\n"),
        "{output}"
    );
    let row = rows_of(&output, "snv_missense_rev")[0];
    let keys: Vec<&str> = extra_keys(row)
        .into_iter()
        .filter(|key| key.starts_with("SpliceAI"))
        .collect();
    assert_eq!(
        keys,
        [
            "SpliceAI_pred_DP_AG",
            "SpliceAI_pred_DP_AL",
            "SpliceAI_pred_DP_DG",
            "SpliceAI_pred_DP_DL",
            "SpliceAI_pred_DS_AG",
            "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL",
            "SpliceAI_pred_SYMBOL",
        ]
    );

    let out = dir.path().join("out.json");
    let run = run_corpus(&corpus, &["--plugin", plugin, "--json"], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    let record: serde_json::Value = output
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .find(|record| record["id"] == "snv_missense_rev")
        .unwrap();
    let consequence = &record["transcript_consequences"][0];
    let spliceai = consequence["spliceai"]
        .as_object()
        .unwrap_or_else(|| panic!("no spliceai object: {consequence}"));
    assert_eq!(spliceai.len(), 9, "{spliceai:?}");
    assert_eq!(spliceai["SYMBOL"], "TMPRSS15");
    for part in [
        "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    ] {
        assert!(spliceai[part].is_number(), "{part}: {spliceai:?}");
    }
    assert!(
        consequence
            .as_object()
            .unwrap()
            .keys()
            .all(|key| !key.starts_with("spliceai_")),
        "{consequence}"
    );

    let out = dir.path().join("cutoff.txt");
    let run = run_corpus(
        &corpus,
        &["--plugin", &format!("{plugin},cutoff=0.5")],
        &out,
    );
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        output.contains(
            "\n## SpliceAI_cutoff : Flag if delta score pass the cutoff (PASS) or if it does not (FAIL)\n"
        ),
        "{output}"
    );
    let mut verdicts = std::collections::BTreeSet::new();
    for row in data_lines(&output) {
        let pred = extra(row, "SpliceAI_pred");
        let flag = extra(row, "SpliceAI_cutoff");
        assert_eq!(pred.is_some(), flag.is_some(), "{row}");
        if let (Some(pred), Some(flag)) = (pred, flag) {
            // Every entry of the value counts: the first stands behind its
            // symbol and each later one behind its allele and symbol.
            let largest = pred
                .split(',')
                .enumerate()
                .flat_map(|(place, entry)| {
                    entry.split('|').skip(1 + usize::from(place > 0)).take(4)
                })
                .map(|score| score.parse::<f64>().unwrap())
                .fold(0.0, f64::max);
            assert_eq!(flag, if largest >= 0.5 { "PASS" } else { "FAIL" }, "{row}");
            verdicts.insert(flag);
        }
    }
    assert_eq!(verdicts.len(), 2, "both verdicts appear: {verdicts:?}");
}

/// The rows of a record in a VCF output's CSQ, each the row's Feature and its
/// fields by the names of the CSQ header line.
fn csq_rows(output: &str, record: &str) -> Vec<(String, BTreeMap<String, String>)> {
    let names: Vec<&str> = output
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ"))
        .and_then(|line| line.split("Format: ").nth(1))
        .map(|format| format.trim_end_matches("\">").split('|').collect())
        .unwrap_or_default();
    let line = data_lines(output)
        .into_iter()
        .find(|line| line.split('\t').nth(2) == Some(record))
        .unwrap_or_else(|| panic!("no line of {record}"));
    let csq = line
        .split('\t')
        .nth(7)
        .and_then(|info| info.split(';').find_map(|part| part.strip_prefix("CSQ=")))
        .unwrap_or_default();
    csq.split(',')
        .map(|row| {
            let fields: BTreeMap<String, String> = names
                .iter()
                .zip(row.split('|'))
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect();
            (fields["Feature"].clone(), fields)
        })
        .collect()
}

/// A row whose gene several entries of one score record name carries every
/// one of them, in the record's order, in each output form, where
/// `SpliceAI.pm:328-332` files the record's whole value under its first
/// entry's gene, `:276-287` cut nine members from that value and `:336-338`
/// read its first four scores. The default and tab formats hold the entries
/// as the record spells them, each later one with its allele part; the VCF
/// format holds each of the nine fields as its values joined by `&`, the CSQ
/// field's spelling of the comma, and the JSON format as one text with the
/// commas; `cutoff=0.5` gives a row one verdict, `PASS` where any of its
/// gene's entries reaches the cutoff. The six records and the score rows they
/// read are made up: the corpus's score files are synthetic.
#[test]
fn spliceai_rows_carry_every_entry_of_their_gene_in_each_form() {
    let corpus = corpus("GRCh38-plugin-spliceai-fixture");
    let plugin = "SpliceAI,snv=plugin_data/spliceai_scores.raw.snv.hg38.synthetic.vcf.gz,\
                  indel=plugin_data/spliceai_scores.raw.indel.hg38.synthetic.vcf.gz,cutoff=0.5";
    let dir = tempfile::tempdir().unwrap();
    // BTG3-AS1 has three transcripts at each site and CXADR one.
    let gene_rows = ["ENST00000626994", "ENST00000701920", "ENST00000807308"];
    let other_row = "ENST00001128578";
    // A record, the value and the verdict of its BTG3-AS1 rows, and those of
    // its CXADR row where an entry names CXADR.
    type Kind = (
        &'static str,
        &'static str,
        &'static str,
        Option<(&'static str, &'static str)>,
    );
    let kinds: [Kind; 6] = [
        (
            "snv_gene_two_entries_equal",
            "BTG3-AS1|0.21|0.76|0.43|0.59|-14|7|-6|24,C|BTG3-AS1|0.21|0.76|0.43|0.59|-14|7|-6|24",
            "PASS",
            None,
        ),
        (
            "snv_gene_two_entries_differ",
            "BTG3-AS1|0.18|0.25|0.09|0.49|16|-32|-5|-2,T|BTG3-AS1|0.76|0.84|0.00|0.38|-5|38|-13|-2",
            "PASS",
            None,
        ),
        (
            "snv_gene_four_entries",
            "BTG3-AS1|0.07|0.05|0.09|0.02|-6|-33|16|-20,C|BTG3-AS1|0.09|0.15|0.07|0.01|5|7|-33|-9,\
             C|BTG3-AS1|0.55|0.59|0.54|0.62|39|26|-39|48,C|BTG3-AS1|0.05|0.19|0.11|0.13|19|46|22|-5",
            "PASS",
            None,
        ),
        (
            "snv_gene_twice_then_other_gene",
            "BTG3-AS1|0.19|0.81|0.32|0.68|-36|-31|-2|-36,G|BTG3-AS1|0.16|0.09|0.13|0.01|-33|-7|-24|42",
            "PASS",
            Some(("CXADR|0.15|0.00|0.09|0.14|-37|-2|44|-39", "FAIL")),
        ),
        (
            "snv_gene_around_other_gene",
            "BTG3-AS1|0.09|0.20|0.26|0.04|37|6|-32|31,A|BTG3-AS1|0.19|0.06|0.14|0.13|-36|-1|-26|28",
            "FAIL",
            Some(("CXADR|0.48|0.25|0.72|0.29|-2|33|26|-29", "PASS")),
        ),
        (
            "snv_other_gene_then_gene_twice",
            "BTG3-AS1|0.16|0.01|0.06|0.10|-46|26|-3|30,T|BTG3-AS1|0.63|0.28|0.74|0.83|-17|-16|50|-6",
            "PASS",
            Some(("CXADR|0.08|0.07|0.19|0.01|-25|-42|-20|-43", "FAIL")),
        ),
    ];
    // The nine members of a value in the order of its parts: each member's
    // values over the value's entries, joined by `separator`.
    let members = |value: &str, separator: &str| -> Vec<String> {
        let entries: Vec<Vec<&str>> = value
            .split(',')
            .enumerate()
            .map(|(place, entry)| entry.split('|').skip(usize::from(place > 0)).collect())
            .collect();
        (0..9)
            .map(|at| {
                let values: Vec<&str> = entries.iter().map(|entry| entry[at]).collect();
                values.join(separator)
            })
            .collect()
    };
    let parts = [
        "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    ];

    let out = dir.path().join("default.txt");
    let run = run_corpus(&corpus, &["--plugin", plugin], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let default_output = std::fs::read_to_string(&out).unwrap();
    let out = dir.path().join("tab.txt");
    let run = run_corpus(&corpus, &["--plugin", plugin, "--tab"], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let tab_output = std::fs::read_to_string(&out).unwrap();
    let tab_columns: Vec<&str> = tab_output
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap()
        .split('\t')
        .collect();
    let tab_column = |name: &str| tab_columns.iter().position(|c| *c == name).unwrap();
    let out = dir.path().join("out.vcf");
    let run = run_corpus(&corpus, &["--plugin", plugin, "--vcf"], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let vcf_output = std::fs::read_to_string(&out).unwrap();
    let out = dir.path().join("out.json");
    let run = run_corpus(&corpus, &["--plugin", plugin, "--json"], &out);
    assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
    let json_records: Vec<serde_json::Value> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();

    for (record, value, verdict, other) in kinds {
        let expected = |feature: &str| -> Option<(&str, &str)> {
            if gene_rows.contains(&feature) {
                Some((value, verdict))
            } else {
                assert_eq!(feature, other_row, "{record}");
                other
            }
        };

        let rows = rows_of(&default_output, record);
        assert_eq!(rows.len(), 4, "{record}");
        for row in rows {
            let feature = row.split('\t').nth(4).unwrap();
            let written = extra(row, "SpliceAI_pred").zip(extra(row, "SpliceAI_cutoff"));
            let wanted =
                expected(feature).map(|(value, verdict)| (value.to_string(), verdict.to_string()));
            assert_eq!(written, wanted, "{record} {feature}: default");
        }

        let rows = rows_of(&tab_output, record);
        assert_eq!(rows.len(), 4, "{record}");
        for row in rows {
            let cells: Vec<&str> = row.split('\t').collect();
            let feature = cells[tab_column("Feature")];
            let written = (
                cells[tab_column("SpliceAI_pred")],
                cells[tab_column("SpliceAI_cutoff")],
            );
            assert_eq!(
                written,
                expected(feature).unwrap_or(("-", "-")),
                "{record} {feature}: tab"
            );
        }

        let rows = csq_rows(&vcf_output, record);
        assert_eq!(rows.len(), 4, "{record}");
        for (feature, fields) in rows {
            let written: Vec<&str> = parts
                .iter()
                .map(|part| fields[&format!("SpliceAI_pred_{part}")].as_str())
                .collect();
            match expected(&feature) {
                Some((value, verdict)) => {
                    assert_eq!(written, members(value, "&"), "{record} {feature}: vcf");
                    assert_eq!(fields["SpliceAI_cutoff"], verdict, "{record} {feature}");
                }
                None => {
                    assert_eq!(written, [""; 9], "{record} {feature}: vcf");
                    assert_eq!(fields["SpliceAI_cutoff"], "", "{record} {feature}");
                }
            }
        }

        let json = json_records
            .iter()
            .find(|json| json["id"] == record)
            .unwrap_or_else(|| panic!("no object of {record}"));
        let rows = json["transcript_consequences"].as_array().unwrap();
        assert_eq!(rows.len(), 4, "{record}");
        for row in rows {
            let feature = row["transcript_id"].as_str().unwrap();
            let Some((value, verdict)) = expected(feature) else {
                assert!(row.get("spliceai").is_none(), "{record} {feature}: {row}");
                continue;
            };
            let object = &row["spliceai"];
            let several = value.contains(',');
            for (part, member) in parts.iter().zip(members(value, ",")) {
                if several || *part == "SYMBOL" {
                    assert_eq!(object[part], member.as_str(), "{record} {feature} {part}");
                } else {
                    let number: f64 = member.parse().unwrap();
                    assert_eq!(
                        object[part].as_f64(),
                        Some(number),
                        "{record} {feature} {part}"
                    );
                }
            }
            assert_eq!(object["SpliceAI_cutoff"], verdict, "{record} {feature}");
            assert_eq!(object.as_object().unwrap().len(), 10, "{record} {feature}");
        }
    }
}

/// The made-up score files of the SpliceAI corpora hold some record lines more
/// than once at an allele. On the corpora's records whose line of two entries stands twice
/// (`snv_two_genes_line_twice`) and whose line of three entries stands three
/// times (`snv_two_genes_line_thrice`), on the record that two lines of other
/// scores score (`snv_two_genes_lines_differ`) and on the record that one line
/// a gene scores (`snv_two_genes_record_per_gene`), every transcript row
/// carries the entry of its own gene once: one `SpliceAI_pred` of nine parts
/// that opens with the row's symbol, the same on every row of a gene.
#[test]
fn spliceai_writes_each_gene_its_entry_once_where_a_record_line_stands_again() {
    for (name, assembly, build) in [
        ("GRCh38-plugin-spliceai-fixture", "GRCh38", "hg38"),
        ("GRCh37-plugin-spliceai-fixture", "GRCh37", "hg19"),
    ] {
        let corpus = corpus(name);
        let plugin = format!(
            "SpliceAI,snv=plugin_data/spliceai_scores.raw.snv.{build}.synthetic.vcf.gz,\
             indel=plugin_data/spliceai_scores.raw.indel.{build}.synthetic.vcf.gz"
        );
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("default.txt");
        let fasta = corpus.join("reference.fa.gz");
        let run = vep()
            .args(args_for(&corpus, assembly))
            .arg("-o")
            .arg(&out)
            .args(["--fasta", fasta.to_str().unwrap(), "--symbol"])
            .args(["--plugin", &plugin])
            .current_dir(&corpus)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
        let output = std::fs::read_to_string(&out).unwrap();
        for record in [
            "snv_two_genes_line_twice",
            "snv_two_genes_line_thrice",
            "snv_two_genes_lines_differ",
            "snv_two_genes_record_per_gene",
        ] {
            let mut by_gene: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
            for row in rows_of(&output, record) {
                let symbol = extra(row, "SYMBOL")
                    .unwrap_or_else(|| panic!("{name} {record}: a row without a symbol: {row}"));
                let pred = extra(row, "SpliceAI_pred")
                    .unwrap_or_else(|| panic!("{name} {record}: no SpliceAI_pred: {row}"));
                let parts: Vec<&str> = pred.split('|').collect();
                assert_eq!(parts.len(), 9, "{name} {record}: {pred}");
                assert!(!pred.contains(','), "{name} {record}: {pred}");
                assert_eq!(parts[0], symbol, "{name} {record}: {pred}");
                by_gene.entry(symbol).or_default().insert(pred);
            }
            assert_eq!(by_gene.len(), 2, "{name} {record}: {by_gene:?}");
            assert!(
                by_gene.values().all(|values| values.len() == 1),
                "{name} {record}: {by_gene:?}"
            );
        }
    }
}

/// A score row whose REF and ALT share a trailing base states the record's
/// change beside an unchanged reference base: the row `21 18353100 TG AG` of the
/// GRCh38 corpus's made-up indel file is the change `T` to `A` at 18353100, the `G` at
/// 18353101 kept. vep-rs trims the shared base and writes the row's values on
/// the substitution `T A` and on the record spelled as the row is alike
/// (`spliceai_trimmed_table_row_match`); Ensembl VEP 116.3 compares the row as
/// written (SpliceAI.pm:272 at VEP_plugins release/116) and scores the
/// substitution nothing.
#[test]
fn spliceai_scores_a_row_spelled_with_a_shared_trailing_base() {
    let corpus = corpus("GRCh38-plugin-spliceai-fixture");
    let dir = tempfile::tempdir().unwrap();
    let input = one_record_vcf(
        dir.path(),
        "21\t18353100\tshared_trailing_base\tT\tA\t.\t.\t.\n\
         21\t18353100\trow_spelling\tTG\tAG\t.\t.\t.",
    );
    let out = dir.path().join("out.txt");
    let plugin = format!(
        "SpliceAI,snv={data}/spliceai_scores.raw.snv.hg38.synthetic.vcf.gz,\
         indel={data}/spliceai_scores.raw.indel.hg38.synthetic.vcf.gz",
        data = corpus.join("plugin_data").display()
    );
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats"])
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args(["--plugin", &plugin])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let body = std::fs::read_to_string(&out).unwrap();
    let row_values = "TMPRSS15|0.82|0.64|0.86|0.67|-25|-36|-29|-20";
    for (record, location, allele) in [
        ("shared_trailing_base", "21:18353100", "A"),
        ("row_spelling", "21:18353100-18353101", "AG"),
    ] {
        let rows = rows_of(&body, record);
        assert_eq!(rows.len(), 1, "{record}: {body}");
        let columns: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(
            &columns[1..5],
            [location, allele, "ENSG00000154646", "ENST00000284885"]
        );
        assert_eq!(
            extra(rows[0], "SpliceAI_pred").as_deref(),
            Some(row_values),
            "{record}: {}",
            rows[0]
        );
    }
}

/// `tests/fixtures/spliceai-two-streams`: seven single-base records at sites
/// of the GRCh38 SpliceAI corpus against a score file written as two bgzip
/// streams joined, so the first stream's end marker, a block that holds no
/// bytes, stands inside the file. `tabix` gives the row that opens the second
/// stream the empty block's offset as its chunk start and reads it from the
/// next block (htslib `bgzf_read_block` skips empty blocks), so the record at
/// that row is scored like the records on either side of it. Three sites hold
/// records no single-base variant takes: an INFO column of `.`, which names no
/// gene (`SpliceAI.pm:328-332`), an ALT column listing three bases and a
/// reference allele of `N`, neither equal to a variant's allele (`:272`). The
/// output is the one Ensembl VEP 116.3 wrote on the fixture under
/// `--no_headers`, byte for byte.
#[test]
fn spliceai_reads_a_score_file_of_two_joined_streams() {
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/spliceai-two-streams");
    let corpus = corpus("GRCh38-plugin-spliceai-fixture");
    let expected = std::fs::read_to_string(fixture.join("ensembl_vep.txt")).unwrap();
    let scored = |record: &str| {
        expected
            .lines()
            .filter(|row| row.starts_with(record) && row.contains("SpliceAI_pred="))
            .count()
    };
    assert_eq!(scored("first_stream_two_genes\t"), 4);
    assert_eq!(scored("opens_second_stream\t"), 1);
    assert_eq!(scored("second_stream_two_alternates\t"), 12);
    for unscored in [
        "no_information\t",
        "alternates_joined\t",
        "any_base_reference\t",
    ] {
        assert_eq!(scored(unscored), 0, "{unscored}");
    }
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.txt");
    let plugin = format!(
        "SpliceAI,snv={},indel={}",
        fixture
            .join("scores_two_streams.synthetic.vcf.gz")
            .display(),
        corpus
            .join("plugin_data/spliceai_scores.raw.indel.hg38.synthetic.vcf.gz")
            .display()
    );
    let output = common::isolated_vep(tmp.path())
        .arg("-i")
        .arg(fixture.join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--no_headers",
            "--plugin",
            &plugin,
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let actual = std::fs::read_to_string(&out).unwrap();
    assert!(
        actual == expected,
        "{}",
        common::compare_headerless(&expected, &actual)
    );
}

/// Runs `vep` on `input` against the LoFTEE corpus's cache, reference and plugin
/// data under the corpus's own `LoF` argument string with `params` appended, in
/// the default format; the output text follows the process.
fn run_loftee(input: &Path, out_dir: &Path, params: &str) -> (Output, String) {
    let corpus = corpus("GRCh38-plugin-loftee");
    let plugin_data = corpus.join("plugin_data");
    let plugin = format!(
        "LoF,loftee_path:{data},human_ancestor_fa:{data}/human_ancestor.fa.gz,gerp_bigwig:{data}/gerp_conservation_scores.homo_sapiens.GRCh38.bw{params}",
        data = plugin_data.display()
    );
    let out = out_dir.join("out.txt");
    let output = common::isolated_vep(out_dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats"])
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .arg("--dir_plugins")
        .arg(&plugin_data)
        .args(["--plugin", &plugin])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    (output, std::fs::read_to_string(&out).unwrap_or_default())
}

/// A one-record VCF of `record` under `dir`.
fn one_record_vcf(dir: &Path, record: &str) -> PathBuf {
    let path = dir.join("in.vcf");
    std::fs::write(
        &path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    path
}

/// The Extra column of the row of `output` on `transcript` for `record`.
fn loftee_extra<'a>(output: &'a str, record: &str, transcript: &str) -> &'a str {
    let column = format!("\t{transcript}\t");
    let row = rows_of(output, record)
        .into_iter()
        .find(|row| row.contains(&column))
        .unwrap_or_else(|| panic!("no row of {record} on {transcript}:\n{output}"));
    row.rsplit('\t').next().unwrap()
}

/// LoF.pm:497 tests the length of the variant's span expanded by four bases
/// each side, nine for every one-base reference span, so Ensembl VEP 116.3
/// flags the one-base deletion of JAM2's exon 7 acceptor A (`21 25709431 TA T`,
/// the feature `A/-` at 25709432, context `TTGTAGAAG`) `NAGNAG_SITE` on each of
/// its five protein-coding acceptor rows, as it flags the substitution there.
#[test]
fn loftee_flags_nagnag_site_on_a_one_base_acceptor_deletion() {
    let dir = tempfile::tempdir().unwrap();
    let input = one_record_vcf(dir.path(), "21\t25709431\tnagnag_del\tTA\tT\t.\t.\t.");
    let (output, body) = run_loftee(&input, dir.path(), "");
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let rows = rows_of(&body, "nagnag_del");
    assert_eq!(rows.len(), 15, "{body}");
    for transcript in [
        "ENST00000312957",
        "ENST00000400532",
        "ENST00000480456",
        "ENST00000948521",
        "ENST00001117607",
    ] {
        assert_eq!(
            loftee_extra(&body, "nagnag_del", transcript),
            "IMPACT=HIGH;STRAND=1;LoF=HC;LoF_flags=NAGNAG_SITE;LoF_info=INTRON_SIZE:3347",
            "{transcript}"
        );
    }
    assert_eq!(
        rows.iter().filter(|row| row.contains(";LoF=")).count(),
        5,
        "{rows:?}"
    );
    assert!(rows
        .iter()
        .all(|row| row.split('\t').nth(1) == Some("21:25709432")
            && row.split('\t').nth(2) == Some("-")));
}

/// `tests/fixtures/loftee-acceptor-on-both-strands`: a fabricated model whose
/// forward-strand transcript `FWDT1` and reverse-strand transcript `REVT1`
/// share their acceptor bases, with one record on each shared acceptor. The
/// nine bases around `toy:400` hold the `NAGNAG_SITE` motif on the forward
/// strand alone and those around `toy:699` on the reverse strand alone.
/// LoF.pm:489-494 reads the nine bases once per record on the strand of the
/// first transcript whose acceptor row it tests, `FWDT1` here, and keeps them
/// for the record's other transcripts, so Ensembl VEP 116.3 wrote the flag on
/// both rows at `toy:400` and on neither at `toy:699`. vep-rs reads them on
/// each row's own strand: its rows of `REVT1` carry the flag at `toy:699`
/// alone, and every other byte of its output is the one Ensembl VEP wrote
/// under `--no_headers`.
#[test]
fn loftee_nagnag_site_follows_the_strand_of_each_rows_transcript() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/loftee-acceptor-on-both-strands");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let plugin = format!(
        "LoF,loftee_path:{},use_gerp_end_trunc:0,run_splice_predictions:0",
        fixture.display()
    );
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(fixture.join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus("GRCh38-plugin-loftee").join("json_cache"))
        .arg("--gff")
        .arg(fixture.join("two_strand_acceptor.gff3.gz"))
        .arg("--fasta")
        .arg(fixture.join("two_strand_acceptor.fa"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--buffer_size", "5000", "--force_overwrite", "--no_stats"])
        .args(["--quiet", "--no_headers"])
        .arg("--dir_plugins")
        .arg(&fixture)
        .args(["--plugin", &plugin])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let actual = std::fs::read_to_string(&out).unwrap();
    let flags = |body: &str| -> Vec<(String, String, String)> {
        body.lines()
            .map(|row| {
                let columns: Vec<&str> = row.split('\t').collect();
                (
                    columns[1].to_string(),
                    columns[4].to_string(),
                    extra_field(row, "LoF_flags").unwrap_or("-").to_string(),
                )
            })
            .collect()
    };
    let row = |location: &str, transcript: &str, flags: &str| {
        (
            location.to_string(),
            transcript.to_string(),
            flags.to_string(),
        )
    };
    assert_eq!(
        flags(&actual),
        vec![
            row("toy:400", "FWDT1", "NAGNAG_SITE"),
            row("toy:400", "REVT1", "NON_CAN_SPLICE"),
            row("toy:699", "FWDT1", "NON_CAN_SPLICE"),
            row("toy:699", "REVT1", "NAGNAG_SITE"),
        ]
    );
    let expected = std::fs::read_to_string(fixture.join("ensembl_vep.txt")).unwrap();
    assert_eq!(
        flags(&expected),
        vec![
            row("toy:400", "FWDT1", "NAGNAG_SITE"),
            row("toy:400", "REVT1", "NON_CAN_SPLICE,NAGNAG_SITE"),
            row("toy:699", "FWDT1", "NON_CAN_SPLICE"),
            row("toy:699", "REVT1", "-"),
        ]
    );
    let as_vep_rs_flags_them = expected
        .replace(
            "LoF_flags=NON_CAN_SPLICE,NAGNAG_SITE;LoF_info=INTRON_SIZE:221",
            "LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:221",
        )
        .replace(
            "LoF=HC;LoF_info=INTRON_SIZE:151",
            "LoF=HC;LoF_flags=NAGNAG_SITE;LoF_info=INTRON_SIZE:151",
        );
    assert!(
        actual == as_vep_rs_flags_them,
        "{}",
        common::compare_headerless(&as_vep_rs_flags_them, &actual)
    );
}

/// LoF.pm:268 tests `lc($self->{check_complete_cds}) eq 'true'`: under
/// `check_complete_cds:1` Ensembl VEP 116.3 writes the corpus's rows as with no
/// parameter, and under `check_complete_cds:TRUE` as under `true`, where the
/// seven frameshift rows on a `cds_end_NF` transcript that carry an exon number
/// read `LoF=LC;LoF_filter=INCOMPLETE_CDS` in place of `LoF=HC` and nothing
/// else moves: LoF.pm:262 asks for an exon number before line 268 is reached,
/// so the corpus's three such rows without one keep `LoF=HC`, and the corpus
/// holds no stop-gained row on such a transcript. (Ensembl VEP's plugin also
/// filters every multi-exon transcript under `true`, a documented divergence
/// this test leaves aside.)
#[test]
fn loftee_check_complete_cds_reads_any_case_of_true_only() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus("GRCh38-plugin-loftee");
    let input = corpus.join("variants.vcf");
    let body_of = |params: &str| {
        let (output, body) = run_loftee(&input, dir.path(), params);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        body
    };
    let off = body_of("");
    let expected = read_expected(&corpus, "default.txt");
    assert_eq!(data_lines(&off), data_lines(&expected));
    assert_eq!(
        data_lines(&body_of(",check_complete_cds:1")),
        data_lines(&off)
    );

    let on = body_of(",check_complete_cds:TRUE");
    assert_eq!(
        data_lines(&body_of(",check_complete_cds:true")),
        data_lines(&on)
    );
    let filtered: Vec<&str> = data_lines(&on)
        .into_iter()
        .filter(|row| row.contains("INCOMPLETE_CDS"))
        .collect();
    assert_eq!(filtered.len(), 7, "{filtered:?}");
    assert!(
        filtered
            .iter()
            .all(|row| row.contains("FLAGS=cds_end_NF;LoF=LC;LoF_filter=INCOMPLETE_CDS;")),
        "{filtered:?}"
    );
    assert_eq!(
        loftee_extra(&on, "2815871", "ENST00000416754"),
        "IMPACT=HIGH;STRAND=-1;FLAGS=cds_end_NF;LoF=LC;LoF_filter=INCOMPLETE_CDS;LoF_info=PERCENTILE:0.556701030927835,GERP_DIST:-10.5988001339138,BP_DIST:42,DIST_FROM_LAST_EXON:4,50_BP_RULE:FAIL"
    );
    let restored: Vec<String> = data_lines(&on)
        .into_iter()
        .map(|row| row.replace(";LoF=LC;LoF_filter=INCOMPLETE_CDS;", ";LoF=HC;"))
        .collect();
    assert_eq!(restored, data_lines(&off));
}

/// LoF.pm:262 asks for an exon number before the `check_complete_cds` test of
/// line 268, and BaseTranscriptVariation.pm:692 gives none to an insertion
/// written just before an exon's first base. Of the corpus's frameshift rows on
/// a `cds_end_NF` transcript, Ensembl VEP therefore filters those that carry an
/// exon number; run at 116.3 under `check_complete_cds:true`, it writes the
/// three of record `977916` whose exon starts at 37478208 as the corpus stores
/// them. `docs/plugins/loftee.md` states those counts.
#[test]
fn loftee_check_complete_cds_leaves_a_row_without_an_exon_number() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus("GRCh38-plugin-loftee");
    let (output, on) = run_loftee(
        &corpus.join("variants.vcf"),
        dir.path(),
        ",check_complete_cds:true",
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let incomplete: Vec<&str> = data_lines(&on)
        .into_iter()
        .filter(|row| {
            let frameshift = row
                .split('\t')
                .nth(6)
                .is_some_and(|terms| terms.split(',').any(|term| term == "frameshift_variant"));
            frameshift
                && extra_field(row, "FLAGS")
                    .is_some_and(|flags| flags.split(',').any(|flag| flag == "cds_end_NF"))
        })
        .collect();
    let (filtered, kept): (Vec<&str>, Vec<&str>) = incomplete
        .iter()
        .partition(|row| extra_field(row, "LoF_filter") == Some("INCOMPLETE_CDS"));
    assert_eq!(
        (incomplete.len(), filtered.len(), kept.len()),
        (10, 7, 3),
        "{incomplete:?}"
    );
    assert_eq!(
        data_lines(&on)
            .into_iter()
            .filter(|row| row.contains("INCOMPLETE_CDS"))
            .count(),
        filtered.len()
    );

    let expected = read_expected(&corpus, "default.txt");
    for transcript in ["ENST00000642309", "ENST00000643551", "ENST00000643854"] {
        let column = format!("\t{transcript}\t");
        let row_in = |rows: Vec<&str>| -> String {
            rows.into_iter()
                .find(|row| row.contains(&column))
                .unwrap_or_else(|| panic!("no row of 977916 on {transcript}"))
                .to_string()
        };
        let row = row_in(kept.clone());
        assert_eq!(row, row_in(rows_of(&expected, "977916")), "{transcript}");
        assert!(
            loftee_extra(&on, "977916", transcript).starts_with(
                "IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;LoF=HC;LoF_flags=NO_EXON_NUMBER;LoF_info=PERCENTILE:"
            ),
            "{row}"
        );
        assert!(!row.contains("50_BP_RULE"), "{row}");
    }

    let page = std::fs::read_to_string(format!(
        "{}/../../docs/plugins/loftee.md",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
    .split_whitespace()
    .collect::<Vec<_>>()
    .join(" ");
    for stated in [
        format!(
            "The {} are among the corpus's {} frameshift rows of `cds_end_NF` transcripts, where the engines agree:",
            filtered.len(),
            incomplete.len()
        ),
        format!(
            "| {} | carried | `LoF=LC;LoF_filter=INCOMPLETE_CDS` |",
            filtered.len()
        ),
        format!(
            "| {}, all of record `977916` | none: an insertion just before an exon's first base, a span that overlaps no exon (`BaseTranscriptVariation.pm`, lines 690-693) |",
            kept.len()
        ),
        "| `LoF=HC;LoF_flags=NO_EXON_NUMBER`: line 262 of `LoF.pm` asks for an exon number before the check of line 268 |".to_string(),
    ] {
        assert!(
            page.contains(&stated),
            "docs/plugins/loftee.md lacks {stated:?}"
        );
    }
}

/// A variant overlapping two exons carries the exon number `first-last/total`
/// (BaseTranscriptVariation.pm:700), which LoF.pm:565-566 numifies to the
/// first exon's rank and LoF.pm:457-458 reads the total of, so Ensembl VEP
/// 116.3 applies the 50-bp rule to a frameshift deletion spanning JAM2's exon
/// 6, intron 6 and exon 7 (`21 25706085 <3350 reference bases> A`, exon
/// `6-7/9` on ENST00000312957) and writes `DIST_FROM_LAST_EXON` and
/// `50_BP_RULE` with no `NO_EXON_NUMBER` flag.
#[test]
fn loftee_reads_an_exon_range_as_its_first_exon() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus("GRCh38-plugin-loftee");
    let reference = vep_fasta::IndexedFasta::from_path(corpus.join("reference.fa.gz")).unwrap();
    let deleted = reference.sequence("21", 25706085, 25709434).unwrap();
    assert_eq!(deleted.len(), 3350);
    assert_eq!(deleted[0], b'A');
    let record = format!(
        "21\t25706085\tmultiexon_del\t{}\tA\t.\t.\t.",
        std::str::from_utf8(&deleted).unwrap()
    );
    let input = one_record_vcf(dir.path(), &record);
    let (output, body) = run_loftee(&input, dir.path(), "");
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(!body.contains("NO_EXON_NUMBER"), "{body}");
    for (transcript, info) in [
        (
            "ENST00000312957",
            "PERCENTILE:0.884664131812421,GERP_DIST:109.19599955529,BP_DIST:89,DIST_FROM_LAST_EXON:57,50_BP_RULE:PASS,INTRON_SIZE:3347",
        ),
        (
            "ENST00000400532",
            "PERCENTILE:0.858359957401491,GERP_DIST:28.8360001593828,BP_DIST:131,DIST_FROM_LAST_EXON:57,50_BP_RULE:PASS,INTRON_SIZE:3347",
        ),
        (
            "ENST00000480456",
            "PERCENTILE:0.898550724637681,GERP_DIST:109.19599955529,BP_DIST:89,DIST_FROM_LAST_EXON:57,50_BP_RULE:PASS,INTRON_SIZE:3347",
        ),
        (
            "ENST00000948521",
            "PERCENTILE:0.903397027600849,GERP_DIST:109.19599955529,BP_DIST:89,DIST_FROM_LAST_EXON:57,50_BP_RULE:PASS,INTRON_SIZE:3347",
        ),
        (
            "ENST00001117607",
            "PERCENTILE:0.829218106995885,GERP_DIST:98.0269994065166,BP_DIST:163,DIST_FROM_LAST_EXON:117,50_BP_RULE:PASS,INTRON_SIZE:3347",
        ),
    ] {
        assert_eq!(
            loftee_extra(&body, "multiexon_del", transcript),
            format!("IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info={info}"),
            "{transcript}"
        );
    }
}
