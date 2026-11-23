// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every loaded built-in plugin's header descriptions reach the three headed
//! formats as Ensembl VEP prints them (`OutputFactory.pm` `get_plugin_headers`,
//! :2416-2430): the default format's `## <KEY> : <description>` lines after the
//! flag fields (`OutputFactory/VEP_output.pm:167`), the tab format's column
//! descriptions and columns (`OutputFactory/Tab.pm:149-156`) and the VCF
//! format's `##<KEY>=<description>` metas and CSQ Format entries
//! (`OutputFactory/VCF.pm:272-273`). The text printed is the plugin's own
//! `header_info`, the same registry loaded in process; the descriptions that
//! quote their data's terms or the data file itself (the CADD non-commercial
//! clause, a dbNSFP README, a TranscriptAnnotator column, a gnomADMt `##INFO`
//! definition) print as their sources state them.

use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use vep_builtins::registry::BuiltinRegistry;
use vep_builtins::traits::{OutputFormat, RunContext};
use vep_core::assembly::Assembly;

mod common;

/// The GRCh38 golden corpora whose committed data files the stack loads.
fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The non-commercial clause CADD appends to each field's description
/// (CADD.pm:100, 174).
const CADD_CLAUSE: &str =
    "CADD is only available here for non-commercial use. See CADD website for more information.";

/// The README written beside the dbNSFP fixture: the numbered `column:
/// description` lines dbNSFP reads for its header (dbNSFP.pm:485-543).
const DBNSFP_README: &str = "A made-up list of columns:\n\
1\tchr: A made-up text for the first column.\n\
11\tSIFT_score: A made-up number for each transcript. The numbers are separated\n\
\t\tby \";\", in the order of the transcript list.\n\
\n\
16\tGERP++_RS: A made-up number for the site,  with a doubled blank in its text.\n";

/// The README's two descriptions as dbNSFP prints them: prefixed with the data
/// file's base name, continuation lines joined, whitespace collapsed.
const DBNSFP_GERP_DESCRIPTION: &str = "(from dbNSFP4.9a_grch38_synthetic) A made-up number for \
                                       the site, with a doubled blank in its text.";
const DBNSFP_SIFT_DESCRIPTION: &str = "(from dbNSFP4.9a_grch38_synthetic) A made-up number for \
                                       each transcript. The numbers are separated by \";\", in \
                                       the order of the transcript list.";

/// The fields gnomADMt prints without `fields=` (gnomADMt.pm:159).
const GNOMADMT_DEFAULT_FIELDS: [&str; 5] = [
    "hap_AC_hom",
    "hap_AC_het",
    "hap_AF_hom",
    "hap_AF_het",
    "hap_AN",
];

/// The `--plugin` arguments of a thirty-plugin stack: every built-in but
/// MaxEntScan, whose matrices a generator writes, each on the data files a
/// golden corpus commits, the file paths absolute. The dbNSFP fixture is copied
/// beside `DBNSFP_README` under `dir`. The stack needs `--hgvs` (HGVSIntronOffset)
/// and a reference (UTRAnnotator).
fn plugin_args(dir: &Path) -> Vec<String> {
    let dbnsfp = golden("GRCh38-plugin-dbnsfp-fixture").join("plugin_data");
    for name in [
        "dbNSFP4.9a_grch38_synthetic.gz",
        "dbNSFP4.9a_grch38_synthetic.gz.tbi",
    ] {
        std::fs::copy(dbnsfp.join(name), dir.join(name)).unwrap();
    }
    std::fs::write(dir.join("dbNSFP4.9a.readme.txt"), DBNSFP_README).unwrap();
    let file = |corpus: &str, path: &str| golden(corpus).join(path).display().to_string();
    vec![
        format!(
            "CADD,snv={},indels={}",
            file(
                "GRCh38-plugin-cadd-fixture",
                "plugin_data/whole_genome_SNVs_synthetic.tsv.gz"
            ),
            file(
                "GRCh38-plugin-cadd-fixture",
                "plugin_data/gnomad.genomes.r4.0.indel.synthetic.tsv.gz"
            )
        ),
        format!(
            "dbNSFP,{},SIFT_score,GERP++_RS",
            dir.join("dbNSFP4.9a_grch38_synthetic.gz").display()
        ),
        format!("TranscriptAnnotator,file={}", transcript_annotator_file()),
        format!("gnomADMt,file={}", gnomadmt_file().display()),
        format!(
            "pLI,{}",
            file("GRCh38-plugin-pli", "plugin_data/pLI_values.txt")
        ),
        "NMD".to_string(),
        format!(
            "LOEUF,file={},match_by=gene",
            file("GRCh38-plugin-loeuf", "data/loeuf_dataset_grch38.tsv.gz")
        ),
        format!(
            "DosageSensitivity,file={}",
            file(
                "GRCh38-plugin-dosage-sensitivity",
                "data/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"
            )
        ),
        format!(
            "PrimateAI,{}",
            file(
                "GRCh38-plugin-primateai-fixture",
                "data/PrimateAI_scores_v0.2_GRCh38_sorted_synthetic.tsv.bgz"
            )
        ),
        format!(
            "ClinPred,file={}",
            file(
                "GRCh38-plugin-clinpred-fixture",
                "data/ClinPred_hg38_sorted_tabbed_synthetic.tsv.gz"
            )
        ),
        format!(
            "BayesDel,file={}",
            file(
                "GRCh38-plugin-bayesdel-fixture",
                "data/BayesDel_170824_addAF_all_scores_synthetic.txt.gz"
            )
        ),
        format!(
            "SpliceVault,file={}",
            file(
                "GRCh38-plugin-splicevault-fixture",
                "data/SpliceVault_data_GRCh38_synthetic.tsv.gz"
            )
        ),
        format!(
            "REVEL,file={}",
            file("GRCh38-plugin-revel", "plugin_data/revel.tsv.gz")
        ),
        format!(
            "AlphaMissense,file={}",
            file(
                "GRCh38-plugin-alphamissense",
                "plugin_data/alphamissense.tsv.gz"
            )
        ),
        format!(
            "dbscSNV,{}",
            file(
                "GRCh38-plugin-dbscsnv-fixture",
                "plugin_data/dbscSNV1.1_synthetic_GRCh38.txt.gz"
            )
        ),
        format!(
            "gnomADc,{}",
            file(
                "GRCh38-plugin-gnomadc",
                "plugin_data/gnomad.coverage.tsv.gz"
            )
        ),
        format!(
            "LoFtool,{}",
            file("GRCh38-plugin-loftool", "plugin_data/LoFtool_scores.txt")
        ),
        format!(
            "SpliceAI,snv={},indel={}",
            file(
                "GRCh38-plugin-spliceai-fixture",
                "plugin_data/spliceai_scores.raw.snv.hg38.synthetic.vcf.gz"
            ),
            file(
                "GRCh38-plugin-spliceai-fixture",
                "plugin_data/spliceai_scores.raw.indel.hg38.synthetic.vcf.gz"
            )
        ),
        format!(
            "Conservation,{}",
            file("GRCh38-plugin-conservation", "plugin_data/gerp.bw")
        ),
        format!(
            "UTRAnnotator,file={}",
            file(
                "GRCh38-plugin-utr-annotator",
                "data/uORF_5UTR_GRCh38_PUBLIC.chr21.txt"
            )
        ),
        format!(
            "AncestralAllele,{}",
            file(
                "GRCh38-plugin-ancestral",
                "ancestral/homo_sapiens_ancestor_21.slice.fa.gz"
            )
        ),
        format!(
            "AVI,file={},raw=1,links=1",
            file("GRCh38-plugin-avi-fixture", "avi_scores.synthetic.tsv.gz")
        ),
        format!(
            "OpenTargets,file={}",
            file("GRCh38-plugin-opentargets", "open_targets_vep.slice.tsv.gz")
        ),
        format!(
            "LoF,loftee_path:{},human_ancestor_fa:{},gerp_bigwig:{}",
            file("GRCh38-plugin-loftee", "plugin_data"),
            file("GRCh38-plugin-loftee", "plugin_data/human_ancestor.fa.gz"),
            file(
                "GRCh38-plugin-loftee",
                "plugin_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
            )
        ),
        format!(
            "GWAS,file={}",
            file("GRCh38-plugin-gwas", "plugin_data/gwas_catalog_slice.tsv")
        ),
        "Blosum62".to_string(),
        "Downstream".to_string(),
        "HGVSIntronOffset".to_string(),
        "NearestExonJB".to_string(),
        "TSSDistance".to_string(),
    ]
}

/// The TranscriptAnnotator data file, as the `file=` parameter types it.
fn transcript_annotator_file() -> String {
    golden("GRCh38-plugin-transcript-annotator")
        .join("data/transcript_annotator.tsv.gz")
        .display()
        .to_string()
}

/// The gnomADMt data file, a bgzipped VCF whose `##INFO` lines define its fields.
fn gnomadmt_file() -> PathBuf {
    golden("GRCh38-plugin-gnomadmt").join("data/gnomad.v3.1.chrM.slice.vcf.bgz")
}

/// Runs `vep` with the stack loaded on the UTRAnnotator corpus, whose committed
/// reference serves `--fasta`, under `--hgvs`, in the format `format_flag`
/// names (the default without one), and returns the output's header lines,
/// `##` and `#` alike.
fn header_lines(out: &Path, format_flag: Option<&str>, plugin_args: &[String]) -> Vec<String> {
    let corpus = golden("GRCh38-plugin-utr-annotator");
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--hgvs"])
        .args(format_flag);
    for arg in plugin_args {
        cmd.args(["--plugin", arg]);
    }
    let output = cmd.output().expect("vep binary runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("WARNING"), "{stderr}");
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// The data lines of `out`.
fn data_lines(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// The stack's header info as the binary builds it: the same arguments loaded
/// into a registry under the run's output format, assembly, `--hgvs` and
/// reference, every plugin loaded, plugin by plugin in load order with each
/// plugin's keys sorted.
fn expected_header_info(plugin_args: &[String], format: OutputFormat) -> Vec<(String, String)> {
    let mut registry = BuiltinRegistry::new();
    let mut run = RunContext::default();
    run.output_format = format;
    run.refseq = false;
    run.merged = false;
    run.hgvs = true;
    run.plugin_dirs = &[];
    let reference = vep_fasta::IndexedFasta::from_path(
        golden("GRCh38-plugin-utr-annotator").join("reference.fa.gz"),
    )
    .unwrap();
    let outcome = registry.load_from_args_for_run(
        plugin_args,
        Some(std::sync::Arc::new(reference)),
        Some(Assembly::Grch38),
        run,
    );
    assert!(outcome.failed.is_empty(), "{:?}", outcome.failed);
    assert!(outcome.unresolved.is_empty(), "{:?}", outcome.unresolved);
    assert_eq!(registry.plugin_count(), plugin_args.len());
    let info = registry.all_header_info();
    assert!(info.len() > plugin_args.len(), "{info:?}");
    info
}

/// The `## <KEY> : <description>` lines between `## Extra column keys:` and the
/// command line of a default-format header: the flag fields, then the plugins'.
fn extra_column_lines(header: &[String]) -> Vec<String> {
    header
        .iter()
        .skip_while(|line| *line != "## Extra column keys:")
        .skip(1)
        .take_while(|line| !line.starts_with("## VEP command-line"))
        .cloned()
        .collect()
}

/// The fields of the CSQ `##INFO` definition's `Format:`.
fn csq_fields(header: &[String]) -> Vec<String> {
    let info = header
        .iter()
        .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
        .unwrap();
    info.split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .map(str::to_owned)
        .collect()
}

/// The default format prints one `## <KEY> : <description>` line per plugin
/// field after the flag fields, the tab format describes every plugin column
/// and adds the columns, and the VCF format writes one `##<KEY>=<description>`
/// meta per field and the CSQ Format entries, each text the plugin's own and
/// in the registry's order.
#[test]
fn every_loaded_plugin_description_reaches_default_tab_and_vcf_headers() {
    let tmp = tempfile::tempdir().unwrap();
    let args = plugin_args(tmp.path());

    let out = tmp.path().join("default.txt");
    let header = header_lines(&out, None, &args);
    let expected = expected_header_info(&args, OutputFormat::Vep);
    let expected_lines: Vec<String> = expected
        .iter()
        .map(|(key, description)| format!("## {key} : {description}"))
        .collect();
    let block = extra_column_lines(&header);
    assert!(block.len() > expected_lines.len(), "{block:#?}");
    let (flag_lines, plugin_lines) = block.split_at(block.len() - expected_lines.len());
    assert_eq!(plugin_lines, expected_lines.as_slice());
    for line in flag_lines {
        let key = line.trim_start_matches("## ").split(" : ").next().unwrap();
        assert!(
            !expected.iter().any(|(k, _)| k == key),
            "plugin key {key} printed among the flag fields: {block:#?}"
        );
    }

    let out = tmp.path().join("tab.txt");
    let header = header_lines(&out, Some("--tab"), &args);
    let expected = expected_header_info(&args, OutputFormat::Tab);
    for (key, description) in &expected {
        let line = format!("## {key} : {description}");
        assert!(header.contains(&line), "{line} missing from {header:#?}");
    }
    let columns: Vec<&str> = header.last().unwrap().split('\t').collect();
    assert_eq!(columns[0], "#Uploaded_variation");
    let keys: Vec<&str> = expected.iter().map(|(key, _)| key.as_str()).collect();
    assert_eq!(&columns[columns.len() - keys.len()..], keys.as_slice());
    let rows = data_lines(&out);
    assert!(!rows.is_empty());
    for row in &rows {
        assert_eq!(row.split('\t').count(), columns.len(), "{row}");
    }

    let out = tmp.path().join("out.vcf");
    let header = header_lines(&out, Some("--vcf"), &args);
    let expected = expected_header_info(&args, OutputFormat::Vcf);
    let keys: Vec<&str> = expected.iter().map(|(key, _)| key.as_str()).collect();
    let metas: Vec<&String> = header
        .iter()
        .filter(|line| {
            keys.iter().any(|key| {
                line.strip_prefix("##")
                    .and_then(|l| l.strip_prefix(*key))
                    .is_some_and(|rest| rest.starts_with('='))
            })
        })
        .collect();
    let expected_metas: Vec<String> = expected
        .iter()
        .map(|(key, description)| format!("##{key}={description}"))
        .collect();
    assert_eq!(
        metas.iter().map(|line| line.as_str()).collect::<Vec<_>>(),
        expected_metas
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
    );
    let csq = csq_fields(&header);
    assert_eq!(&csq[csq.len() - keys.len()..], keys.as_slice());
}

/// The `Description` of a `##INFO` definition in the gnomADMt data file.
fn info_description(id: &str) -> String {
    let file = std::fs::File::open(gnomadmt_file()).unwrap();
    let prefix = format!("##INFO=<ID={id},");
    BufReader::new(MultiGzDecoder::new(file))
        .lines()
        .map(Result::unwrap)
        .take_while(|line| line.starts_with("##"))
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("no {prefix} line in the gnomADMt data file"))
        .split_once("Description=\"")
        .unwrap()
        .1
        .trim_end_matches("\">")
        .to_string()
}

/// A description that quotes its source prints as the source states it, in
/// every headed format: CADD's fields end with the non-commercial clause
/// (CADD.pm:174); a dbNSFP column is described from the README beside the data
/// file, prefixed with the file's base name (dbNSFP.pm:263-284, 485-543); a
/// TranscriptAnnotator column is `column from <file>` with the path as typed
/// (TranscriptAnnotator.pm:196); a gnomADMt field carries the data file's own
/// `##INFO` description (gnomADMt.pm:149, 245).
#[test]
fn licence_and_data_derived_descriptions_print_as_their_sources_state_them() {
    let tmp = tempfile::tempdir().unwrap();
    let args = plugin_args(tmp.path());
    let mut expected = vec![
        (
            "CADD_PHRED".to_string(),
            format!("PHRED-like scaled CADD score. {CADD_CLAUSE}"),
        ),
        (
            "CADD_RAW".to_string(),
            format!("Raw CADD score. {CADD_CLAUSE}"),
        ),
        ("GERP++_RS".to_string(), DBNSFP_GERP_DESCRIPTION.to_string()),
        (
            "SIFT_score".to_string(),
            DBNSFP_SIFT_DESCRIPTION.to_string(),
        ),
    ];
    for column in ["Comment", "SIFT_pred", "SIFT_score"] {
        expected.push((
            format!("transcript_annotator_{column}"),
            format!("column from {}", transcript_annotator_file()),
        ));
    }
    for field in GNOMADMT_DEFAULT_FIELDS {
        expected.push((format!("gnomAD_{field}"), info_description(field)));
    }

    let default_header = header_lines(&tmp.path().join("default.txt"), None, &args);
    let tab_header = header_lines(&tmp.path().join("tab.txt"), Some("--tab"), &args);
    let vcf_header = header_lines(&tmp.path().join("out.vcf"), Some("--vcf"), &args);
    for (key, description) in &expected {
        let line = format!("## {key} : {description}");
        assert!(
            default_header.contains(&line),
            "{line} missing from {default_header:#?}"
        );
        assert!(
            tab_header.contains(&line),
            "{line} missing from {tab_header:#?}"
        );
        let meta = format!("##{key}={description}");
        assert!(
            vcf_header.contains(&meta),
            "{meta} missing from {vcf_header:#?}"
        );
    }
}
