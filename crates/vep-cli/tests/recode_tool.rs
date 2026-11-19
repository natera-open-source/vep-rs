// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep recode`, through the binary: the behaviour of Ensembl VEP's
//! `variant_recoder` (release/116.2, `VariantRecoder.pm`, whose line numbers
//! the tests cite) on the inputs of the `tests/golden/116/GRCh38-recoder`
//! corpus over its own cache, and the corpus's
//! `tool_cases`, each naming the arguments of one `variant_recoder` run and the
//! file its standard output was recorded in.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

use serde_json::Value;

mod common;

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-recoder")
}

fn manifest() -> Value {
    serde_json::from_str(&fs::read_to_string(corpus_dir().join("manifest.json")).unwrap()).unwrap()
}

/// The cache, FASTA and synonyms flags the corpus manifest names, as paths
/// below the corpus directory.
fn resources() -> Vec<String> {
    let mut flags = cache_and_fasta();
    flags.push("--synonyms".to_string());
    flags.push(resource_path("synonyms"));
    flags
}

/// The cache and FASTA flags alone: without the synonyms file every notation
/// stays on the chromosome's own name.
fn cache_and_fasta() -> Vec<String> {
    vec![
        "--json_cache".to_string(),
        resource_path("cache"),
        "--fasta".to_string(),
        resource_path("fasta"),
    ]
}

/// The path the manifest names under `key`, below the corpus directory.
fn resource_path(key: &str) -> String {
    corpus_dir()
        .join(
            manifest()[key]
                .as_str()
                .unwrap_or_else(|| panic!("manifest has no {key}")),
        )
        .to_string_lossy()
        .into_owned()
}

/// Runs `vep recode` with `args` and the corpus's resources, feeding `stdin`
/// when given.
fn recode_in(args: &[&str], stdin: Option<&str>) -> Output {
    recode_in_with(args, &resources(), stdin)
}

/// Runs `vep recode` with `args` and the resource flags given, feeding `stdin`
/// when given.
fn recode_in_with(args: &[&str], resources: &[String], stdin: Option<&str>) -> Output {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("recode")
        .args(args)
        .args(resources)
        .current_dir(corpus_dir());
    match stdin {
        Some(text) => {
            let mut child = cmd
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            child
                .stdin
                .take()
                .unwrap()
                .write_all(text.as_bytes())
                .unwrap();
            child.wait_with_output().unwrap()
        }
        None => cmd.stdin(Stdio::null()).output().unwrap(),
    }
}

/// Runs `vep recode` with `args` and the corpus's resources under `home` as
/// the home directory and `extra` added to the environment, for a run that
/// reads `<home>/.vep/vep.ini` or a `VEP_*` variable.
fn recode_in_home(args: &[&str], home: &Path, extra: &[(&str, &str)]) -> Output {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.env("HOME", home);
    for (key, value) in extra {
        cmd.env(key, value);
    }
    cmd.arg("recode")
        .args(args)
        .args(resources())
        .current_dir(corpus_dir())
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The parsed output of a run that must succeed quietly.
fn recoded(args: &[&str]) -> Vec<Value> {
    let output = recode_in(args, None);
    assert!(
        output.status.success(),
        "{args:?} exited {:?}: {}",
        output.status.code(),
        stderr(&output)
    );
    assert_eq!(stderr(&output), "", "{args:?}");
    parse(&stdout(&output), args)
}

fn parse(text: &str, context: impl std::fmt::Debug) -> Vec<Value> {
    serde_json::from_str::<Value>(text)
        .unwrap_or_else(|e| panic!("{context:?}: {e}: {text}"))
        .as_array()
        .unwrap_or_else(|| panic!("{context:?}: not an array: {text}"))
        .clone()
}

/// The one `ERROR: [recode]` line of a refused run, exit 1.
fn refused(args: &[&str]) -> String {
    let output = recode_in(args, None);
    assert_eq!(
        output.status.code(),
        Some(1),
        "{args:?}: {}",
        stderr(&output)
    );
    let text = stderr(&output);
    assert!(
        text.starts_with("ERROR: [recode] ") && text.matches('\n').count() == 1,
        "{args:?}: {text:?}"
    );
    text.trim_end().to_string()
}

/// The string values of an allele object's key, in order.
fn strings(object: &Value, key: &str) -> Vec<String> {
    object[key]
        .as_array()
        .unwrap_or_else(|| panic!("{key} is not an array in {object}"))
        .iter()
        .map(|v| {
            v.as_str()
                .unwrap_or_else(|| panic!("{key} holds {v}"))
                .to_string()
        })
        .collect()
}

/// The one object of a one-input run, as (allele, object) pairs.
fn alleles(results: &[Value]) -> Vec<(String, Value)> {
    assert_eq!(results.len(), 1, "{results:?}");
    results[0]
        .as_object()
        .unwrap()
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect()
}

fn no_refseq_entry(values: &[String]) {
    for value in values {
        assert!(
            !["NM_", "NR_", "XM_", "XR_", "NP_", "XP_", "LRG_"]
                .iter()
                .any(|prefix| value.starts_with(prefix)),
            "{value} names a RefSeq or LRG model, which an Ensembl cache has none of"
        );
    }
}

fn unique(values: &[String]) {
    let set: BTreeSet<&String> = values.iter().collect();
    assert_eq!(set.len(), values.len(), "{values:?} repeats a value");
}

/// An identifier with two alternate alleles becomes one object keyed by its
/// alleles (308-312); each allele's object carries the input, the identifier,
/// one HGVSg and one SPDI on the chromosome's `NC_` accession, the HGVSc and
/// HGVSp of every overlapping transcript as unique arrays (`find_in_ref`,
/// `Utils.pm:337-368`) and the allele's VCF form (318-375). The HGVSc set is
/// the set `vep --json` writes for that allele.
#[test]
fn rsid_recodes_to_every_notation_per_allele() {
    let results = recoded(&["--id", "rs779189159", "--vcf_string"]);
    let alleles = alleles(&results);
    assert_eq!(
        alleles.iter().map(|(a, _)| a.as_str()).collect::<Vec<_>>(),
        ["A", "G"]
    );
    for (allele, object) in &alleles {
        assert_eq!(object["input"], "rs779189159");
        assert_eq!(strings(object, "id"), ["rs779189159"]);
        assert_eq!(
            strings(object, "hgvsg"),
            [format!("NC_000021.9:g.33553711C>{allele}")]
        );
        assert_eq!(
            strings(object, "spdi"),
            [format!("NC_000021.9:33553710:C:{allele}")]
        );
        assert_eq!(
            strings(object, "vcf_string"),
            [format!("21-33553711-C-{allele}")]
        );
        let hgvsc = strings(object, "hgvsc");
        let hgvsp = strings(object, "hgvsp");
        assert!(
            hgvsc.contains(&format!("ENST00000356577.10:c.4480C>{allele}")),
            "{hgvsc:?}"
        );
        unique(&hgvsc);
        unique(&hgvsp);
        no_refseq_entry(&hgvsc);
        no_refseq_entry(&hgvsp);
        let keys: Vec<&str> = object
            .as_object()
            .unwrap()
            .keys()
            .map(String::as_str)
            .collect();
        assert_eq!(
            keys,
            [
                "hgvsc",
                "hgvsg",
                "hgvsp",
                "id",
                "input",
                "spdi",
                "vcf_string"
            ]
        );
    }
    let (_, a) = &alleles[0];
    assert!(strings(a, "hgvsp").contains(&"ENSP00000348984.4:p.Gln1494Lys".to_string()));
    let (_, g) = &alleles[1];
    assert!(strings(g, "hgvsp").contains(&"ENSP00000348984.4:p.Gln1494Glu".to_string()));

    // The regrouping reads the annotation run's own JSON: the HGVSc set of an
    // allele is the set of `hgvsc` over the run's entries for that allele.
    let dir = tempfile::tempdir().unwrap();
    let json_out = dir.path().join("vep.json");
    let mut vep = common::isolated_vep_anywhere();
    vep.args([
        "--id",
        "rs779189159",
        "--format",
        "id",
        "--json",
        "--hgvs",
        "--no_stats",
        "--quiet",
        "-o",
    ])
    .arg(&json_out)
    .args(resources());
    let output = vep.output().unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let record: Value = serde_json::from_str(
        fs::read_to_string(&json_out)
            .unwrap()
            .lines()
            .next()
            .unwrap(),
    )
    .unwrap();
    for (allele, object) in &alleles {
        let from_run: BTreeSet<String> = record["transcript_consequences"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|tc| tc["variant_allele"] == *allele)
            .filter_map(|tc| tc.get("hgvsc").and_then(Value::as_str).map(str::to_string))
            .collect();
        let from_recode: BTreeSet<String> = strings(object, "hgvsc").into_iter().collect();
        assert_eq!(from_recode, from_run, "allele {allele}");
    }
}

/// A transcript (`c.`) notation and a protein (`p.`) notation of one variant
/// each become one object whose allele carries the same genomic notation, SPDI
/// and identifier, and the object is keyed by the input as written.
#[test]
fn hgvs_c_and_p_inputs_recode_to_genomic_and_spdi() {
    let results = recoded(&[
        "--input_data",
        "ENST00000480456.6:c.685C>T\nENSP00000420419.1:p.Arg229Ter\n21:g.25702257C>T",
    ]);
    assert_eq!(results.len(), 3);
    for (result, input) in results.iter().zip([
        "ENST00000480456.6:c.685C>T",
        "ENSP00000420419.1:p.Arg229Ter",
        "21:g.25702257C>T",
    ]) {
        let object = &result["T"];
        assert_eq!(object["input"], input, "{result}");
        assert_eq!(strings(object, "hgvsg"), ["NC_000021.9:g.25702257C>T"]);
        assert_eq!(strings(object, "spdi"), ["NC_000021.9:25702256:C:T"]);
        assert_eq!(strings(object, "id"), ["rs781261918"]);
        assert!(strings(object, "hgvsc").contains(&"ENST00000480456.6:c.685C>T".to_string()));
        assert!(strings(object, "hgvsp").contains(&"ENSP00000420419.1:p.Arg229Ter".to_string()));
    }
}

/// An SPDI input, on the accession or on the chromosome's own name, recodes to
/// its notations; a deletion's VCF form carries the anchor base read from the
/// reference and its allele is `-`.
#[test]
fn spdi_input_recodes() {
    let results = recoded(&[
        "--input_data",
        "NC_000021.9:32694271:TTTGT:\n21:25702256:C:T",
        "--vcf_string",
    ]);
    assert_eq!(results.len(), 2);
    let deletion = &results[0]["-"];
    assert_eq!(deletion["input"], "NC_000021.9:32694271:TTTGT:");
    assert_eq!(
        strings(deletion, "hgvsg"),
        ["NC_000021.9:g.32694272_32694276del"]
    );
    assert_eq!(strings(deletion, "spdi"), ["NC_000021.9:32694271:TTTGT:"]);
    assert_eq!(strings(deletion, "vcf_string"), ["21-32694271-ATTTGT-A"]);
    assert_eq!(strings(deletion, "id"), ["rs886039732"]);
    assert!(strings(deletion, "hgvsc").contains(&"ENST00000674351.1:c.741_745del".to_string()));
    let substitution = &results[1]["T"];
    assert_eq!(substitution["input"], "21:25702256:C:T");
    assert_eq!(strings(substitution, "spdi"), ["NC_000021.9:25702256:C:T"]);
    assert_eq!(strings(substitution, "vcf_string"), ["21-25702257-C-T"]);
}

/// `vcf_string` is matched to its allele as the recoder matches it (318-375):
/// one form per alternate allele of a multi-allelic record, every form built
/// from the record's alleles by `to_VCF_record` (a VCF input's columns are not
/// copied, since no input parser keeps them) with the reference's anchor base,
/// so an indel written from a repeat's first base is written from its
/// minimised position, and an alternate allele whose minimised and
/// VCF-converted pairs match no form gets none (502-531, 578-612).
#[test]
fn vcf_string_minimises_and_matches_alleles_as_perl() {
    let vcf = "21\t33553711\trs779189159\tC\tA,G\t.\t.\t.\n\
               21\t33543115\trs1601245742\tTTTTT\tTTTTTT\t.\t.\t.\n\
               21\t45510105\trs1160505542\tGCGGGGCATCCGCGGGGC\tGCGGGGC,GCGGGGCATCCGCGGGGCATCCGCGGGGC\t.\t.\t.\n";
    let results = recoded(&["--input_data", vcf, "--fields", "spdi", "--vcf_string"]);
    assert_eq!(results.len(), 3);
    let multi = &results[0];
    assert_eq!(strings(&multi["A"], "vcf_string"), ["21-33553711-C-A"]);
    assert_eq!(strings(&multi["G"], "vcf_string"), ["21-33553711-C-G"]);
    assert_eq!(
        multi["A"]["input"],
        "21\t33553711\trs779189159\tC\tA,G\t.\t.\t."
    );
    // The record's alleles are minimised to `-/T` at the repeat's last base,
    // and the VCF form is rebuilt there with that base as its anchor.
    let insertion = &results[1];
    assert_eq!(strings(&insertion["T"], "vcf_string"), ["21-33543119-T-TT"]);
    assert_eq!(
        strings(&insertion["T"], "spdi"),
        ["NC_000021.9:33543119::T"]
    );
    // A multi-allelic indel's alleles lose the VCF anchor base; the forms are
    // rebuilt from the trimmed alleles, and the insertion allele's pair
    // matches no form.
    let indel = &results[2];
    assert_eq!(
        strings(&indel["CGGGGC"], "vcf_string"),
        ["21-45510106-CGGGGCATCCGC-C"]
    );
    assert!(
        indel["CGGGGCATCCGCGGGGCATCCGCGGGGC"]
            .get("vcf_string")
            .is_none(),
        "{indel}"
    );

    // A deletion and an insertion from an HGVS notation take the base before
    // the record from the reference as their anchor.
    let results = recoded(&[
        "--input_data",
        "NC_000021.9:g.34886862del\nNC_000021.9:g.37505341_37505342insT",
        "--fields",
        "spdi",
        "--vcf_string",
    ]);
    assert_eq!(
        strings(&results[0]["-"], "vcf_string"),
        ["21-34886861-GG-G"]
    );
    assert_eq!(
        strings(&results[1]["T"], "vcf_string"),
        ["21-37505341-G-GT"]
    );
}

/// The co-located variants are the cache rows a database-backed run finds
/// (`AnnotationType/Variation.pm:91-109`), whatever the input's spelling of
/// its own variant: a row is kept when its span overlaps the record and an
/// allele matches, so an identifier in another case finds its row under the
/// cache's spelling and on the alleles the row spells alone, a deletion
/// minimised to a repeat's last base finds the row dbSNP spells from the
/// first, and so does an SPDI of it; an insertion overlaps only a row
/// spanning both of its flanking bases, so an identifier whose row is spelled
/// from the repeat's first base gets no `id` at all, as the recoder gives
/// none.
#[test]
fn colocated_variants_are_found_by_span_and_allele_as_the_database_finds_them() {
    let results = recoded(&["--id", "RS781261918", "--fields", "id"]);
    assert_eq!(results[0]["T"]["input"], "RS781261918");
    assert_eq!(strings(&results[0]["T"], "id"), ["rs781261918"]);

    let results = recoded(&[
        "--input_data",
        "21\t25702257\trs781261918\tC\tT,A\t.\t.\t.",
        "--fields",
        "id",
    ]);
    assert_eq!(strings(&results[0]["T"], "id"), ["rs781261918"]);
    assert!(results[0]["A"].get("id").is_none(), "{}", results[0]);

    let results = recoded(&["--id", "rs2057997486", "--fields", "id"]);
    assert_eq!(strings(&results[0]["-"], "id"), ["rs2057997486"]);
    let results = recoded(&[
        "--input_data",
        "NC_000021.9:34886861:G:\nNC_000021.9:44293155:TG:",
        "--fields",
        "id",
    ]);
    assert_eq!(strings(&results[0]["-"], "id"), ["rs2057997486"]);
    assert_eq!(strings(&results[1]["-"], "id"), ["rs759281545"]);

    let results = recoded(&["--id", "rs1601245742", "--fields", "id,spdi"]);
    assert_eq!(
        strings(&results[0]["T"], "spdi"),
        ["NC_000021.9:33543119::T"]
    );
    assert!(results[0]["T"].get("id").is_none(), "{}", results[0]);
}

/// `--mane_select` adds one `hgvsg`, `hgvsc`, `hgvsp` triple per MANE Select
/// transcript entry, once per input (439-463), and switches the three HGVS
/// notations on for the triple while `--fields` keeps them out of the allele's
/// own keys.
#[test]
fn mane_select_triples_are_unique_per_input() {
    let results = recoded(&["--id", "rs781261918", "--mane_select"]);
    let object = &results[0]["T"];
    assert_eq!(
        object["mane_select"],
        serde_json::json!([{
            "hgvsc": "ENST00000480456.6:c.685C>T",
            "hgvsg": "NC_000021.9:g.25702257C>T",
            "hgvsp": "ENSP00000420419.1:p.Arg229Ter"
        }])
    );
    assert!(strings(object, "hgvsc").len() > 1);
    let results = recoded(&["--id", "rs781261918", "--mane_select", "--fields", "spdi"]);
    let object = &results[0]["T"];
    let keys: Vec<&str> = object
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(keys, ["input", "mane_select", "spdi"]);
    assert_eq!(object["mane_select"].as_array().unwrap().len(), 1);
    // A variant with no MANE Select transcript gets no key.
    let results = recoded(&["--id", "rs2838519", "--mane_select"]);
    assert!(
        results[0]["A"].get("mane_select").is_none(),
        "{}",
        results[0]
    );
}

/// `--fields` selects the keys, any key of the JSON output among them
/// (`find_in_ref` walks every entry); `--pretty` indents three spaces per
/// level with ` : ` after each key and ends with a newline, as `JSON->pretty`
/// does; the compact form is one line.
#[test]
fn fields_selects_keys_and_pretty_indents() {
    let results = recoded(&["--id", "rs781261918", "--fields", "hgvsg"]);
    let keys: Vec<&str> = results[0]["T"]
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(keys, ["hgvsg", "input"]);
    let results = recoded(&[
        "--id",
        "rs781261918",
        "--fields",
        "transcript_id,consequence_terms,strand",
    ]);
    let object = &results[0]["T"];
    assert!(strings(object, "transcript_id").contains(&"ENST00000480456".to_string()));
    assert!(strings(object, "consequence_terms").contains(&"stop_gained".to_string()));
    assert_eq!(object["strand"], serde_json::json!([1, -1]));

    let output = recode_in(
        &["--id", "rs781261918", "--fields", "hgvsg", "--pretty"],
        None,
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(
        stdout(&output),
        "[\n   {\n      \"T\" : {\n         \"hgvsg\" : [\n            \"NC_000021.9:g.25702257C>T\"\n         ],\n         \"input\" : \"rs781261918\"\n      }\n   }\n]\n"
    );
    let output = recode_in(&["--id", "rs781261918", "--fields", "hgvsg"], None);
    assert_eq!(
        stdout(&output),
        "[{\"T\":{\"hgvsg\":[\"NC_000021.9:g.25702257C>T\"],\"input\":\"rs781261918\"}}]"
    );
}

/// A repeated option runs with its last value, as `variant_recoder`'s
/// `GetOptions` (36-95) assigns each occurrence of a scalar or switch over the
/// last: the run with every option doubled writes what the single options
/// write. A repeated `--transcript_filter` keeps every filter (66).
#[test]
fn repeated_options_run_with_the_last_value() {
    let single = recode_in(
        &["--id", "rs781261918", "--fields", "hgvsg", "--pretty"],
        None,
    );
    assert!(single.status.success(), "{}", stderr(&single));
    let doubled = recode_in(
        &[
            "--id",
            "rs267606642",
            "--id",
            "rs781261918",
            "--fields",
            "transcript_id",
            "--fields",
            "hgvsg",
            "--pretty",
            "--pretty",
        ],
        None,
    );
    assert!(doubled.status.success(), "{}", stderr(&doubled));
    assert_eq!(stderr(&doubled), "");
    assert_eq!(stdout(&doubled), stdout(&single));
    assert!(stdout(&doubled).contains("NC_000021.9:g.25702257C>T"));
    let every = recoded(&["--id", "rs781261918", "--fields", "transcript_id"]);
    assert!(strings(&every[0]["T"], "transcript_id").len() > 1);
    let filtered = recoded(&[
        "--id",
        "rs781261918",
        "--fields",
        "transcript_id",
        "--transcript_filter",
        "stable_id match ENST00000480456",
        "--transcript_filter",
        "stable_id match ENST",
    ]);
    assert_eq!(
        strings(&filtered[0]["T"], "transcript_id"),
        ["ENST00000480456"]
    );
}

/// `--config` (41) names a file `Config->new` reads as it reads `vep`'s
/// (`Config.pm:701`, `:1068-1110`): `<name> <value>` lines, a `#` line and a
/// blank skipped, one dash before a name dropped (`:1086`), a flag the command
/// line also gives keeping the command line's value (`:1101`), a repeatable
/// flag appending, and a key `vep recode` has no option for dropped. The
/// options `variant_recoder` presets before the file is read, `fields` among
/// them (VariantRecoder.pm:111, `:1101`'s `||=`), the ones that add a field
/// before it (`var_synonyms`, `:135`) and the ones `is_valid_param` (`:1046`)
/// does not know (`pretty`) are left as the command line set them. `--debug`
/// (94) is accepted and changes nothing in the output. A file that cannot be
/// read ends the run with the `ERROR: [config]` line (`:1075`), unless `--help`
/// is given, which `variant_recoder:97` answers before any file is read.
#[test]
fn config_file_and_debug_are_read_as_variant_recoder_reads_them() {
    let dir = tempfile::tempdir().unwrap();
    let ini = dir.path().join("recode.ini");
    std::fs::write(
        &ini,
        "# the fields\n-fields hgvsg\n\npretty 1\ninput_data rs781261918\nno_such_option 1\n\
         ambiguous_hgvs 1\ntranscript_filter \"stable_id match ENST00000480456\"\nconfig other.ini\n",
    )
    .unwrap();
    let ini = ini.display().to_string();
    let default_fields = recode_in(
        &[
            "--id",
            "rs781261918",
            "--transcript_filter",
            "stable_id match ENST00000480456",
        ],
        None,
    );
    assert!(
        default_fields.status.success(),
        "{}",
        stderr(&default_fields)
    );
    let pretty = recode_in(
        &["--id", "rs781261918", "--fields", "hgvsg", "--pretty"],
        None,
    );
    assert!(pretty.status.success(), "{}", stderr(&pretty));
    assert_ne!(stdout(&default_fields), stdout(&pretty));

    // The file's `input_data` and `transcript_filter` apply; its `fields` and
    // `pretty` are the recoder's own, so the output is the default field set,
    // unindented.
    let from_file = recode_in(&["--config", &ini], None);
    assert!(from_file.status.success(), "{}", stderr(&from_file));
    assert_eq!(stderr(&from_file), "");
    assert_eq!(stdout(&from_file), stdout(&default_fields));

    let command_line_wins = recode_in(&["--config", &ini, "--fields", "transcript_id"], None);
    assert!(
        command_line_wins.status.success(),
        "{}",
        stderr(&command_line_wins)
    );
    let results = parse(&stdout(&command_line_wins), "--config with --fields");
    assert_eq!(
        strings(&results[0]["T"], "transcript_id"),
        ["ENST00000480456"],
        "the file's --transcript_filter applied"
    );

    let debug = recode_in(
        &[
            "--id",
            "rs781261918",
            "--fields",
            "hgvsg",
            "--pretty",
            "--debug",
        ],
        None,
    );
    assert!(debug.status.success(), "{}", stderr(&debug));
    assert_eq!(stderr(&debug), "");
    assert_eq!(stdout(&debug), stdout(&pretty));

    let missing = recode_in(&["--config", "/no/such/recode.ini"], None);
    assert_eq!(missing.status.code(), Some(1));
    assert_eq!(
        stderr(&missing),
        "ERROR: [config] Could not open config file \"/no/such/recode.ini\"\n"
    );

    // `variant_recoder:97` prints the usage before the file is read.
    let help = recode_in(&["--help", "--config", "/no/such/recode.ini"], None);
    assert!(help.status.success(), "{}", stderr(&help));
    assert!(stdout(&help).contains("--config"), "{}", stdout(&help));

    // A `~` in the path opens the run's home; `database`, a flag of the
    // annotation tool alone, is a name no option of `vep recode` has, dropped;
    // `var_synonyms 1` adds no field, the field list being composed before the
    // file is read (VariantRecoder.pm:135, `:152`).
    let home = Path::new(env!("CARGO_TARGET_TMPDIR")).join("home");
    fs::create_dir_all(&home).unwrap();
    fs::write(
        home.join("recode_config_tilde.ini"),
        "input_data rs781261918\ndatabase 1\nvar_synonyms 1\n",
    )
    .unwrap();
    let tilde = recode_in(
        &["--config", "~/recode_config_tilde.ini", "--fields", "hgvsg"],
        None,
    );
    assert!(tilde.status.success(), "{}", stderr(&tilde));
    let results = parse(&stdout(&tilde), "--config ~/recode_config_tilde.ini");
    let keys: Vec<&str> = results[0]["T"]
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(keys, ["hgvsg", "input"]);
}

/// Behind `--config`, `Config->new` reads `<dir>/vep.ini` and the `VEP_*`
/// variables for `variant_recoder` as it does for `vep` (`Config.pm:704-708`),
/// `<dir>` being `~/.vep` or the `--config` file's `dir`, each source yielding
/// to the ones above it (`:1101`, `:1158`): an `input_data` from `vep.ini`
/// names the variant recoded, `VEP_INPUT_DATA` gives way to it, a narrower
/// `transcript_filter` from the environment applies beside the ini's, and the
/// file's `input_data` beats the ini's. The recoder's presets stay its own: a
/// `fields` line in `vep.ini` or `VEP_FIELDS`, written for `vep`'s output
/// columns, changes nothing (VariantRecoder.pm:111, `:152`). A directory named
/// by `--config` reads as nothing (`:1075-1077`) and `--config 0` is refused
/// with the empty name (`:1073`), as `vep` reads them.
#[test]
fn vep_ini_and_the_environment_are_read_behind_the_config_file() {
    let dir = tempfile::tempdir().unwrap();
    let home = dir.path();
    fs::create_dir_all(home.join(".vep")).unwrap();
    fs::write(
        home.join(".vep").join("vep.ini"),
        "input_data rs781261918\nfields Uploaded_variation,Location,Allele,Gene\n\
         transcript_filter \"stable_id match ENST\"\n",
    )
    .unwrap();
    let plain = recode_in(
        &[
            "--id",
            "rs781261918",
            "--transcript_filter",
            "stable_id match ENST00000480456",
        ],
        None,
    );
    assert!(plain.status.success(), "{}", stderr(&plain));
    let results = parse(&stdout(&plain), "plain");
    assert!(
        results[0]["T"].get("hgvsg").is_some(),
        "the default fields: {}",
        results[0]
    );

    // The ini's filter keeps every Ensembl transcript, so its run names more
    // than one; the variable's narrower filter, applied beside it, leaves one.
    let from_ini = recode_in_home(&[], home, &[]);
    assert!(from_ini.status.success(), "{}", stderr(&from_ini));
    assert_eq!(stderr(&from_ini), "");
    let results = parse(&stdout(&from_ini), "vep.ini alone");
    assert!(
        results[0]["T"].get("hgvsg").is_some()
            && results[0]["T"]["hgvsc"].as_array().unwrap().len() > 1,
        "vep.ini alone: its input and filter apply, its fields line does not: {}",
        results[0]
    );

    let env_applies = recode_in_home(
        &[],
        home,
        &[
            ("VEP_INPUT_DATA", "rs2838519"),
            ("VEP_FIELDS", "hgvsg"),
            ("VEP_TRANSCRIPT_FILTER", "stable_id match ENST00000480456"),
        ],
    );
    assert!(env_applies.status.success(), "{}", stderr(&env_applies));
    assert_eq!(
        stdout(&env_applies),
        stdout(&plain),
        "the ini's input over the variable's; the variable's filter applied beside the ini's; \
         the variable's fields sealed"
    );

    let config = home.join("recode.ini");
    fs::write(&config, "input_data rs267606642\n").unwrap();
    let config = config.display().to_string();
    let file_wins = recode_in_home(
        &["--config", &config],
        home,
        &[("VEP_INPUT_DATA", "rs2838519")],
    );
    assert!(file_wins.status.success(), "{}", stderr(&file_wins));
    let results = parse(&stdout(&file_wins), "--config over vep.ini");
    let input: BTreeSet<String> = results[0]
        .as_object()
        .unwrap()
        .values()
        .filter_map(|allele| allele.get("input").and_then(Value::as_str))
        .map(str::to_string)
        .collect();
    assert_eq!(
        input,
        BTreeSet::from(["rs267606642".to_string()]),
        "{}",
        results[0]
    );

    let other = home.join("other");
    fs::create_dir_all(&other).unwrap();
    fs::write(other.join("vep.ini"), "input_data rs2838519\n").unwrap();
    let placed = home.join("placed.ini");
    fs::write(&placed, format!("dir {}\n", other.display())).unwrap();
    let placed = placed.display().to_string();
    let moved_ini = recode_in_home(&["--config", &placed], home, &[]);
    assert!(moved_ini.status.success(), "{}", stderr(&moved_ini));
    let results = parse(&stdout(&moved_ini), "--config naming dir");
    let input: BTreeSet<String> = results[0]
        .as_object()
        .unwrap()
        .values()
        .filter_map(|allele| allele.get("input").and_then(Value::as_str))
        .map(str::to_string)
        .collect();
    assert_eq!(
        input,
        BTreeSet::from(["rs2838519".to_string()]),
        "{}",
        results[0]
    );

    let as_dir = recode_in_home(&["--config", &other.display().to_string()], home, &[]);
    assert!(as_dir.status.success(), "{}", stderr(&as_dir));
    assert_eq!(
        stdout(&as_dir),
        stdout(&from_ini),
        "a directory reads as nothing, so the ini alone applies"
    );
    let zero = recode_in_home(&["--config", "0"], home, &[]);
    assert_eq!(zero.status.code(), Some(1));
    assert_eq!(
        stderr(&zero),
        "ERROR: [config] Could not open config file \"\"\n"
    );
}

/// The known synonyms of the variant come from the cache's `var_synonyms`
/// column, one `source: id,id` string per source (425-437), for an identifier
/// input and for a notation that co-locates the variant; a variant without
/// synonyms gets no key.
#[test]
fn var_synonyms_come_from_the_cache_column() {
    let results = recoded(&["--id", "rs267606642", "--var_synonyms", "--fields", "id"]);
    assert_eq!(
        strings(&results[0]["G"], "var_synonyms"),
        ["ClinVar: RCV000003481,VCV000003317", "OMIM: 607358.0011"]
    );
    let results = recoded(&[
        "--id",
        "NC_000021.9:g.44286663T>G",
        "--var_synonyms",
        "--fields",
        "id",
    ]);
    assert_eq!(
        strings(&results[0]["G"], "var_synonyms"),
        ["ClinVar: RCV000003481,VCV000003317", "OMIM: 607358.0011"]
    );
    let results = recoded(&["--id", "rs2838519", "--var_synonyms", "--fields", "id"]);
    assert!(
        results[0]["A"].get("var_synonyms").is_none(),
        "{}",
        results[0]
    );
}

/// `--ga4gh_vrs` adds one GA4GH VRS Allele per SPDI (474-482), its interval
/// start the SPDI's position as text and its end a number, and switches
/// `spdi` on beside the fields named.
#[test]
fn ga4gh_vrs_follows_each_spdi() {
    let results = recoded(&["--id", "rs886039732", "--fields", "hgvsg", "--ga4gh_vrs"]);
    let object = &results[0]["-"];
    assert_eq!(strings(object, "spdi"), ["NC_000021.9:32694271:TTTGT:"]);
    assert_eq!(
        object["ga4gh_vrs"],
        serde_json::json!([{
            "location": {
                "type": "SequenceLocation",
                "interval": {"type": "SimpleInterval", "start": "32694271", "end": 32694276},
                "sequence_id": "refseq:NC_000021.9"
            },
            "state": {"type": "SequenceState", "sequence": ""},
            "type": "Allele"
        }])
    );
}

/// On a chromosome whose synonyms name no `NC_` accession the SPDI stays on the
/// chromosome's own name, and `--ga4gh_vrs` writes the key holding an empty
/// list: `ga4gh_vrs_from_spdi` returns the empty list for such an SPDI
/// (`Variation/Utils/Sequence.pm:1433-1434`), and the push of 479 makes the
/// list and appends nothing to it.
#[test]
fn ga4gh_vrs_is_an_empty_list_without_an_nc_accession() {
    let output = recode_in_with(
        &["--id", "rs779189159", "--fields", "spdi", "--ga4gh_vrs"],
        &cache_and_fasta(),
        None,
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(stderr(&output), "");
    assert_eq!(
        stdout(&output),
        "[{\"A\":{\"ga4gh_vrs\":[],\"input\":\"rs779189159\",\"spdi\":[\"21:33553710:C:A\"]},\"G\":{\"ga4gh_vrs\":[],\"input\":\"rs779189159\",\"spdi\":[\"21:33553710:C:G\"]}}]"
    );
}

/// On a sequence whose RefSeq name is no `NC_` accession, a scaffold or a
/// patch, the SPDI and the HGVSg both take that name, as the recoder's
/// database slice names it (`VariationFeature.pm:1912-1916`, `:2150-2151`),
/// and `--ga4gh_vrs` writes an empty list, the SPDI being off an `NC_` name
/// (`Variation/Utils/Sequence.pm:1434`). The table is built here: one line
/// pairing the corpus's chromosome with a name of that form.
#[test]
fn spdi_takes_the_refseq_name_of_a_sequence_without_an_nc_accession() {
    let dir = tempfile::tempdir().unwrap();
    let synonyms = dir.path().join("chr_synonyms.txt");
    fs::write(&synonyms, "21\tNW_000000021.1\n").unwrap();
    let mut resources = cache_and_fasta();
    resources.push("--synonyms".to_string());
    resources.push(synonyms.to_string_lossy().into_owned());
    let output = recode_in_with(
        &[
            "--id",
            "rs779189159",
            "--fields",
            "hgvsg,spdi",
            "--ga4gh_vrs",
        ],
        &resources,
        None,
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(stderr(&output), "");
    let results: Vec<Value> = serde_json::from_str(&stdout(&output)).unwrap();
    let alleles = results[0].as_object().unwrap();
    assert!(!alleles.is_empty());
    for (allele, entry) in alleles {
        for key in ["hgvsg", "spdi"] {
            let values = entry[key].as_array().unwrap();
            assert!(!values.is_empty(), "{allele} {key}");
            for value in values {
                let text = value.as_str().unwrap();
                assert!(
                    text.starts_with("NW_000000021.1:"),
                    "{allele} {key}: {text}"
                );
            }
        }
        assert_eq!(entry["ga4gh_vrs"], serde_json::json!([]), "{allele}");
    }
}

/// On a row that overlaps no transcript the SPDI names the sequence as the
/// cache names it, whatever the FASTA calls its record: the database slice the
/// recoder reads names itself (`VariationFeature.pm:2151`), and
/// `variant_recoder` writes the sequence's own name there for a record under
/// that name, a table synonym or a prefixed spelling. The FASTA is built here:
/// one record named with a `chr` prefix, the reference base at the record's
/// position.
#[test]
fn spdi_off_a_transcript_takes_the_cache_name_not_the_fasta_record_name() {
    let dir = tempfile::tempdir().unwrap();
    let position = 1000;
    let sequence = format!("{}A{}", "N".repeat(position - 1), "N".repeat(60));
    let fasta = dir.path().join("prefixed.fa");
    fs::write(&fasta, format!(">chr21\n{sequence}\n")).unwrap();
    fs::write(
        dir.path().join("prefixed.fa.fai"),
        format!(
            "chr21\t{0}\t{1}\t{0}\t{2}\n",
            sequence.len(),
            ">chr21\n".len(),
            sequence.len() + 1
        ),
    )
    .unwrap();
    let resources = [
        "--json_cache".to_string(),
        resource_path("cache"),
        "--fasta".to_string(),
        fasta.to_string_lossy().into_owned(),
    ];
    let record = format!("21\t{position}\toff_transcript\tA\tC\t.\t.\t.");
    let output = recode_in_with(
        &[
            "--input_data",
            &record,
            "--format",
            "vcf",
            "--fields",
            "hgvsg,spdi",
        ],
        &resources,
        None,
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(stderr(&output), "");
    let results: Vec<Value> = serde_json::from_str(&stdout(&output)).unwrap();
    let alleles = results[0].as_object().unwrap();
    assert!(!alleles.is_empty());
    for (allele, entry) in alleles {
        for (key, start) in [("hgvsg", "21:g."), ("spdi", "21:")] {
            let values = entry[key].as_array().unwrap();
            assert!(!values.is_empty(), "{allele} {key}");
            for value in values {
                let text = value.as_str().unwrap();
                assert!(text.starts_with(start), "{allele} {key}: {text}");
            }
        }
    }
}

/// `--fields hgvs_offset` names a key no entry carries: the field sets the
/// configuration key of its own name (120), which the output never reads, so
/// no HGVS notation is computed and each allele's object holds the input
/// alone; beside `hgvsc`, whose key switches HGVSc on (`OutputFactory.pm:1750`),
/// the insertion's offset is written under it, the five bases its reverse-strand
/// transcripts shift it by times their strand (`OutputFactory.pm:1755`); beside
/// `hgvsp` alone it is not, since the protein branch (`:1758-1767`) reads an
/// offset only the transcript notation sets (`TranscriptVariationAllele.pm:1447`).
/// The insertion is given as an SPDI, which the recoder shifts; given as its
/// identifier, whose stored class is `indel`, it is left unshifted and carries
/// no offset (`TranscriptVariationAllele.pm:123`).
#[test]
fn hgvs_offset_field_switches_no_notation_on() {
    let insertion = "NC_000021.9:33543119::T";
    let results = recoded(&["--input_data", insertion, "--fields", "hgvs_offset"]);
    assert_eq!(results, [serde_json::json!({"T": {"input": insertion}})]);
    let results = recoded(&["--input_data", insertion, "--fields", "hgvsc,hgvs_offset"]);
    let object = &results[0]["T"];
    assert!(!strings(object, "hgvsc").is_empty(), "{object}");
    assert_eq!(object["hgvs_offset"], serde_json::json!([-5]), "{object}");
    let results = recoded(&["--input_data", insertion, "--fields", "hgvsp,hgvs_offset"]);
    let object = &results[0]["T"];
    let keys: Vec<&str> = object
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(keys, ["hgvsp", "input"], "{object}");
}

/// The warnings raised while an input is read go to the next object written
/// (485-486, with the recoder's buffer of one input), a run that writes no
/// object is one object of its warnings alone (`variant_recoder:116`), an
/// input whose record has no consequence entry is `null`, and standard input
/// is read like a file.
#[test]
fn warnings_are_attached_to_their_input() {
    let results = recoded(&[
        "--input_data",
        "rs781261918\nrs_not_a_thing\nrs2123538435",
        "--fields",
        "hgvsg",
    ]);
    assert_eq!(results.len(), 2);
    assert!(results[0].get("warnings").is_none(), "{}", results[0]);
    assert_eq!(
        results[1]["warnings"],
        serde_json::json!(["No variant found with ID 'rs_not_a_thing'"])
    );
    assert_eq!(results[1]["A"]["input"], "rs2123538435");

    let results = recoded(&["--id", "rs_nope", "--fields", "hgvsg"]);
    assert_eq!(
        results,
        [serde_json::json!({"warnings": ["No variant found with ID 'rs_nope'"]})]
    );

    // A failed variant with no alleles has no consequence entry: `null`.
    let results = recoded(&["--id", "rs1979841086", "--fields", "hgvsg"]);
    assert_eq!(results, [Value::Null]);

    let output = recode_in(
        &["--format", "id", "--fields", "hgvsg"],
        Some("rs_not_a_thing\nrs2123538435\n"),
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let results = parse(&stdout(&output), "stdin");
    assert_eq!(results.len(), 1);
    assert_eq!(
        results[0]["warnings"],
        serde_json::json!(["No variant found with ID 'rs_not_a_thing'"])
    );
    assert_eq!(
        strings(&results[0]["A"], "hgvsg"),
        ["NC_000021.9:g.31693357G>A"]
    );

    // `--warning_file` receives the warnings the output carries, written as
    // the warnings file of `vep` is, at the first warning: a run without one
    // leaves no file.
    let dir = tempfile::tempdir().unwrap();
    let warnings = dir.path().join("recode_warnings.txt");
    let warnings_arg = warnings.to_str().unwrap();
    recoded(&[
        "--id",
        "rs_nope",
        "--fields",
        "hgvsg",
        "--warning_file",
        warnings_arg,
    ]);
    assert_eq!(
        fs::read_to_string(&warnings).unwrap(),
        "WARNING: No variant found with ID 'rs_nope'\n"
    );
    fs::remove_file(&warnings).unwrap();
    recoded(&[
        "--id",
        "rs781261918",
        "--fields",
        "hgvsg",
        "--warning_file",
        warnings_arg,
    ]);
    assert!(
        !warnings.exists(),
        "a run without a warning wrote the warnings file"
    );
}

/// The flags that name an Ensembl database connection, LRG mapping, ambiguous
/// protein notations and an alignment file are refused naming the flag, as
/// Getopt::Long spells it; `--mane_select` is refused for a species other than
/// human (`variant_recoder:104-106`); `--refseq` is refused with the recoder's
/// configuration text, since `VariantRecoder.pm:98` forces `merged` on and
/// `Config.pm:964` checks the `refseq => [gencode_basic merged]` pair (`:621`)
/// in that order.
#[test]
fn database_flags_are_refused_naming_the_flag() {
    assert!(
        refused(&["--id", "rs781261918", "--host", "x"]).starts_with(
            "ERROR: [recode] --host is not supported: it names the Ensembl database connection"
        )
    );
    assert!(refused(&["--id", "rs781261918", "--db_version", "116"])
        .starts_with("ERROR: [recode] --db_version is not supported"));
    assert!(refused(&["--id", "rs781261918", "--lrg"])
        .starts_with("ERROR: [recode] --lrg is not supported"));
    assert!(refused(&["--id", "rs781261918", "--ambiguous_hgvs"])
        .starts_with("ERROR: [recode] --ambiguous_hgvs is not supported"));
    assert!(refused(&["--id", "rs781261918", "--bam", "a.bam"])
        .starts_with("ERROR: [recode] --bam is not supported"));
    assert!(refused(&["--id", "rs781261918", "-u", "anonymous"])
        .starts_with("ERROR: [recode] --user is not supported"));
    assert_eq!(
        refused(&[
            "--id",
            "rs781261918",
            "--species",
            "mus_musculus",
            "--mane_select"
        ]),
        "ERROR: [recode] MANE Select is only available for human"
    );
    assert_eq!(
        refused(&["--id", "rs781261918", "--refseq"]),
        "ERROR: [recode] Can't use --refseq and --merged together"
    );
    assert_eq!(
        refused(&["--id", "rs781261918", "--refseq", "--merged"]),
        "ERROR: [recode] Can't use --refseq and --merged together"
    );
    assert_eq!(
        refused(&["--id", "rs781261918", "--refseq", "--gencode_basic"]),
        "ERROR: [recode] Can't use --refseq and --gencode_basic together"
    );
    // The species check (variant_recoder:104-106) runs before the configuration
    // dies (:108).
    assert_eq!(
        refused(&[
            "--id",
            "rs781261918",
            "--species",
            "mus_musculus",
            "--mane_select",
            "--refseq"
        ]),
        "ERROR: [recode] MANE Select is only available for human"
    );
    // No argument prints the usage and exits 0 (variant_recoder:97).
    let mut cmd = common::isolated_vep_anywhere();
    let output = cmd.arg("recode").output().unwrap();
    assert!(output.status.success());
    assert!(stdout(&output).contains("--input_data"));
}

/// Copies the directory `from` to `to`, creating it.
fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// The corpus of a merged cache: its records, its FASTA and its `json_cache`,
/// whose `info.json` declares `source_type` merged.
fn merged_corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-merged")
}

/// The standard output of `vep recode` over the merged corpus's records and
/// FASTA, with `cache` naming the cache and `flags` added; the run succeeds.
fn recode_merged_records(cache: &[String], flags: &[&str]) -> String {
    let corpus = merged_corpus();
    let input = corpus.join("variants.vcf").to_string_lossy().into_owned();
    let mut args = vec!["-i", input.as_str()];
    args.extend(flags);
    let mut resources = cache.to_vec();
    resources.push("--fasta".to_string());
    resources.push(
        corpus
            .join("reference.fa.gz")
            .to_string_lossy()
            .into_owned(),
    );
    let output = recode_in_with(&args, &resources, None);
    assert!(output.status.success(), "{}", stderr(&output));
    stdout(&output)
}

/// The cache flags naming `root` in Ensembl VEP's layout at cache version 116.
fn layout(root: &Path) -> [String; 4] {
    [
        "--dir_cache".to_string(),
        root.to_string_lossy().into_owned(),
        "--cache_version".to_string(),
        "116".to_string(),
    ]
}

/// `variant_recoder` reads the merged transcript set on every run
/// (`VariantRecoder.pm:98` presets `merged`, and `CacheDir.pm:312` suffixes the
/// species directory with it), so a cache that holds the set is read without
/// `--merged`: a `--json_cache` directory whose `info.json` declares it, and a
/// JSON cache in the `_merged` species directory of a cache root, the one
/// `--merged` resolves (a `<version>_<assembly>` name or the bare version)
/// where the root has no directory of the Ensembl set, and the one named as
/// the Ensembl set's cache beside it, also where the `_merged` directory
/// holds another assembly's cache as well or the entry's `info.json` declares
/// no set. Each run writes what the run under `--merged` writes, RefSeq
/// notations among them.
#[test]
fn a_merged_cache_is_read_without_the_flag() {
    let cache = merged_corpus().join("json_cache");
    let directory = [
        "--json_cache".to_string(),
        cache.to_string_lossy().into_owned(),
    ];
    let flagged = recode_merged_records(&directory, &["--merged"]);
    assert!(flagged.contains("\"NM_"), "no RefSeq notation:\n{flagged}");
    assert_eq!(recode_merged_records(&directory, &[]), flagged);

    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("merged_root");
    copy_dir(&cache, &root.join("homo_sapiens_merged/116_GRCh38"));
    assert_eq!(recode_merged_records(&layout(&root), &[]), flagged);
    assert_eq!(
        recode_merged_records(&layout(&root), &["--assembly", "GRCh38"]),
        flagged
    );
    let ensembl = PathBuf::from(resource_path("cache"));
    copy_dir(&ensembl, &root.join("homo_sapiens/116_GRCh38"));
    assert_eq!(recode_merged_records(&layout(&root), &[]), flagged);
    copy_dir(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-merged/json_cache"),
        &root.join("homo_sapiens_merged/116_GRCh37"),
    );
    assert_eq!(recode_merged_records(&layout(&root), &[]), flagged);

    let undeclared = dir.path().join("undeclared_root");
    copy_dir(&ensembl, &undeclared.join("homo_sapiens/116_GRCh38"));
    let unlabelled = undeclared.join("homo_sapiens_merged/116_GRCh38");
    copy_dir(&cache, &unlabelled);
    let mut info: Value =
        serde_json::from_str(&fs::read_to_string(unlabelled.join("info.json")).unwrap()).unwrap();
    info.as_object_mut().unwrap().remove("source_type").unwrap();
    fs::write(unlabelled.join("info.json"), info.to_string()).unwrap();
    fs::remove_file(unlabelled.join("CHECKSUMS.sha256")).unwrap();
    assert_eq!(recode_merged_records(&layout(&undeclared), &[]), flagged);

    let bare = dir.path().join("bare_root");
    copy_dir(&cache, &bare.join("homo_sapiens_merged/116"));
    assert_eq!(recode_merged_records(&layout(&bare), &[]), flagged);
    assert_eq!(
        recode_merged_records(&layout(&bare), &["--assembly", "GRCh38"]),
        flagged
    );
}

/// A cache of the Ensembl set is read as it is unless a JSON merged cache sits
/// under the same directory name: the run keeps its cache version and
/// assembly. An Ensembl VEP Storable cache in the `_merged` directory (neither
/// `info.json` nor `transcripts/`), a merged cache of another cache version,
/// one of another assembly, under its own name or under the same bare version
/// name, and an entry declaring the Ensembl set, with `--assembly` or without,
/// each leave the Ensembl set's output, and a `--json_cache` directory
/// declaring the Ensembl set, or no set, is read without the flag as it is. A
/// root whose directory of the Ensembl set holds no cache of the version is
/// refused as that directory is, its merged cache left to `--merged`, and an
/// empty root's refusal names the Ensembl set's directory.
#[test]
fn only_a_merged_cache_of_the_same_directory_name_replaces_the_ensembl_set() {
    let merged = merged_corpus().join("json_cache");
    let ensembl = PathBuf::from(resource_path("cache"));
    let dir = tempfile::tempdir().unwrap();
    let root_with = |name: &str| {
        let root = dir.path().join(name);
        copy_dir(&ensembl, &root.join("homo_sapiens/116_GRCh38"));
        root
    };
    let both = [&[][..], &["--assembly", "GRCh38"][..]];

    let alone = root_with("ensembl_root");
    let unmerged = recode_merged_records(&layout(&alone), &[]);
    assert!(!unmerged.contains("\"NM_"), "{unmerged}");

    let storable = alone.join("homo_sapiens_merged/116_GRCh38");
    fs::create_dir_all(storable.join("21")).unwrap();
    fs::write(storable.join("info.txt"), "species\thomo_sapiens\n").unwrap();
    for flags in both {
        assert_eq!(recode_merged_records(&layout(&alone), flags), unmerged);
    }

    let other_version = root_with("other_version_root");
    copy_dir(
        &merged,
        &other_version.join("homo_sapiens_merged/115_GRCh38"),
    );
    for flags in both {
        assert_eq!(
            recode_merged_records(&layout(&other_version), flags),
            unmerged
        );
    }

    let other_assembly = root_with("other_assembly_root");
    copy_dir(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-merged/json_cache"),
        &other_assembly.join("homo_sapiens_merged/116_GRCh37"),
    );
    for flags in both {
        assert_eq!(
            recode_merged_records(&layout(&other_assembly), flags),
            unmerged
        );
    }

    let bare = dir.path().join("bare_root");
    copy_dir(&ensembl, &bare.join("homo_sapiens/116"));
    copy_dir(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-merged/json_cache"),
        &bare.join("homo_sapiens_merged/116"),
    );
    for flags in both {
        assert_eq!(recode_merged_records(&layout(&bare), flags), unmerged);
    }

    let declared_ensembl =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-regulatory");
    let mislabelled = root_with("mislabelled_root");
    copy_dir(
        &declared_ensembl.join("json_cache"),
        &mislabelled.join("homo_sapiens_merged/116_GRCh38"),
    );
    for flags in both {
        assert_eq!(
            recode_merged_records(&layout(&mislabelled), flags),
            unmerged
        );
    }

    let input = declared_ensembl
        .join("variants.vcf")
        .to_string_lossy()
        .into_owned();
    let cache = [
        "--json_cache".to_string(),
        declared_ensembl
            .join("json_cache")
            .to_string_lossy()
            .into_owned(),
    ];
    let output = recode_in_with(&["-i", &input, "--fields", "id"], &cache, None);
    assert!(output.status.success(), "{}", stderr(&output));

    let unlabelled = dir.path().join("unlabelled_cache");
    copy_dir(&merged, &unlabelled);
    let mut info: Value =
        serde_json::from_str(&fs::read_to_string(unlabelled.join("info.json")).unwrap()).unwrap();
    info.as_object_mut().unwrap().remove("source_type").unwrap();
    fs::write(unlabelled.join("info.json"), info.to_string()).unwrap();
    fs::remove_file(unlabelled.join("CHECKSUMS.sha256")).unwrap();
    let cache = [
        "--json_cache".to_string(),
        unlabelled.to_string_lossy().into_owned(),
    ];
    assert_ne!(
        recode_merged_records(&cache, &[]),
        recode_merged_records(&cache, &["--merged"])
    );

    let stale = dir.path().join("stale_root");
    copy_dir(&ensembl, &stale.join("homo_sapiens/115_GRCh38"));
    copy_dir(&merged, &stale.join("homo_sapiens_merged/116_GRCh38"));
    let records = merged_corpus()
        .join("variants.vcf")
        .to_string_lossy()
        .into_owned();
    let output = recode_in_with(&["-i", &records], &layout(&stale), None);
    assert_eq!(
        stderr(&output).trim_end(),
        "ERROR: [cache] No cache found for homo_sapiens, version 116"
    );

    let empty = dir.path().join("empty_root");
    fs::create_dir_all(&empty).unwrap();
    let output = recode_in_with(&["-i", &records], &layout(&empty), None);
    assert_eq!(
        stderr(&output).trim_end(),
        format!(
            "ERROR: [cache] Cache directory {}/homo_sapiens not found",
            empty.display()
        )
    );
}

/// Rewrites the `info.json` of the JSON cache at `cache` to name `assembly`,
/// and drops the checksum manifest, which would refuse the rewritten file
/// first.
fn set_assembly(cache: &Path, assembly: &str) {
    fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
    let info = cache.join("info.json");
    let mut document: Value = serde_json::from_str(&fs::read_to_string(&info).unwrap()).unwrap();
    document["assembly"] = Value::from(assembly);
    fs::write(info, document.to_string()).unwrap();
}

/// Damages the JSON cache at `cache` as `damage` names: a shard its checksum
/// manifest lists removed, an `info.json` that is not JSON, one whose assembly
/// is empty, no `info.json`, or no `transcripts/`; the two that rewrite
/// `info.json` drop the manifest, which would refuse the rewritten file first.
fn damage_cache(cache: &Path, damage: &str) {
    let info = cache.join("info.json");
    match damage {
        "shard" => fs::remove_file(cache.join("transcripts/21/10000001-11000000.json.gz")).unwrap(),
        "not_json" => {
            fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
            fs::write(info, "{not json").unwrap();
        }
        "empty_assembly" => set_assembly(cache, ""),
        "no_info" => fs::remove_file(info).unwrap(),
        "no_transcripts" => fs::remove_dir_all(cache.join("transcripts")).unwrap(),
        other => panic!("unknown damage {other}"),
    }
}

/// A merged cache the run reads without the flag is refused for what `--merged`
/// refuses it for, never passed over for the Ensembl set. With a shard its
/// checksum manifest lists removed, an `info.json` that is not JSON, one whose
/// assembly is empty, no `info.json` or no `transcripts/`, the run ends on the
/// line the run under `--merged` ends on, beside the Ensembl set's cache and
/// under a root with the merged set alone.
#[test]
fn a_damaged_merged_cache_is_refused_as_under_the_flag() {
    let dir = tempfile::tempdir().unwrap();
    let input = merged_corpus()
        .join("variants.vcf")
        .to_string_lossy()
        .into_owned();
    let refusal = |root: &Path, flags: &[&str]| {
        let mut args = vec!["-i", input.as_str()];
        args.extend(flags);
        let output = recode_in_with(&args, &layout(root), None);
        assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
        stderr(&output)
    };
    for damage in [
        "shard",
        "not_json",
        "empty_assembly",
        "no_info",
        "no_transcripts",
    ] {
        for beside_ensembl in [true, false] {
            let root = dir.path().join(format!("{damage}_{beside_ensembl}"));
            let merged = root.join("homo_sapiens_merged/116_GRCh38");
            copy_dir(&merged_corpus().join("json_cache"), &merged);
            if beside_ensembl {
                copy_dir(
                    &PathBuf::from(resource_path("cache")),
                    &root.join("homo_sapiens/116_GRCh38"),
                );
            }
            damage_cache(&merged, damage);
            let flagged = refusal(&root, &["--merged"]);
            assert!(flagged.starts_with("ERROR: [cache] "), "{flagged}");
            assert_eq!(
                refusal(&root, &[]),
                flagged,
                "{damage}, beside the Ensembl set: {beside_ensembl}"
            );
        }
    }
}

/// The merged run is given an assembly only where the directory of the Ensembl
/// set's cache has the name the merged run looks up under it,
/// `<cache_version>_<assembly>`; under any other name the run resolves the
/// `_merged` directory beside it as `--merged` does. Under the bare cache
/// version it reads the entry of the same name, also one whose `info.json`
/// holds an empty assembly, and under `116_hg38`, an assembly `--assembly`
/// reads as `GRCh38`, it reads that entry too. Beside a second entry of the
/// cache version, under the bare version or a name the version only begins,
/// it ends on the line the run under `--merged` ends on.
#[test]
fn under_any_other_directory_name_the_merged_directory_resolves_as_under_the_flag() {
    let dir = tempfile::tempdir().unwrap();
    let merged = merged_corpus().join("json_cache");
    let root_with = |name: &str, leaf: &str| {
        let root = dir.path().join(name);
        copy_dir(
            &PathBuf::from(resource_path("cache")),
            &root.join("homo_sapiens").join(leaf),
        );
        copy_dir(&merged, &root.join("homo_sapiens_merged").join(leaf));
        root
    };
    let input = merged_corpus()
        .join("variants.vcf")
        .to_string_lossy()
        .into_owned();
    let refusal = |root: &Path, flags: &[&str]| {
        let mut args = vec!["-i", input.as_str()];
        args.extend(flags);
        let output = recode_in_with(&args, &layout(root), None);
        assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
        stderr(&output)
    };

    let bare = root_with("bare_root", "116");
    let flagged = recode_merged_records(&layout(&bare), &["--merged"]);
    assert!(flagged.contains("\"NM_"), "no RefSeq notation:\n{flagged}");
    assert_eq!(recode_merged_records(&layout(&bare), &[]), flagged);

    let unnamed = root_with("unnamed_assembly_root", "116");
    damage_cache(&unnamed.join("homo_sapiens_merged/116"), "empty_assembly");
    assert_eq!(
        recode_merged_records(&layout(&unnamed), &[]),
        recode_merged_records(&layout(&unnamed), &["--merged"])
    );

    let alias = root_with("alias_root", "116_hg38");
    for set in ["homo_sapiens", "homo_sapiens_merged"] {
        set_assembly(&alias.join(set).join("116_hg38"), "hg38");
    }
    assert_eq!(
        recode_merged_records(&layout(&alias), &[]),
        recode_merged_records(&layout(&alias), &["--merged"])
    );

    let several = root_with("several_root", "116");
    copy_dir(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-merged/json_cache"),
        &several.join("homo_sapiens_merged/116_GRCh37"),
    );
    let longer = root_with("longer_name_root", "1160_GRCh38");
    copy_dir(&merged, &longer.join("homo_sapiens_merged/116_GRCh38"));
    for root in [&several, &longer] {
        let flagged = refusal(root, &["--merged"]);
        assert!(
            flagged.starts_with("ERROR: [cache] Multiple assemblies found"),
            "{flagged}"
        );
        assert_eq!(refusal(root, &[]), flagged);
    }
}

/// The transcript-id prefixes the manifest's masks drop from the reference:
/// RefSeq models, which an Ensembl cache does not carry, and LRG models.
fn masked_prefixes(masks: &[String]) -> Vec<&'static str> {
    let mut prefixes = Vec::new();
    for mask in masks {
        match mask.as_str() {
            "refseq_transcript_entries" => {
                prefixes.extend(["NM_", "NR_", "XM_", "XR_", "NP_", "XP_"])
            }
            "lrg_transcript_entries" => prefixes.push("LRG_"),
            "archive_dbsnp_synonyms" | "identifier_indel_transcript_notations" => {}
            other => panic!("unknown mask {other}"),
        }
    }
    prefixes
}

/// The manifest's mask classes as the comparison applies them.
struct Masks {
    /// The transcript-id prefixes dropped from the notations.
    prefixes: Vec<&'static str>,
    /// `archive_dbsnp_synonyms`: the `Archive dbSNP` synonym entries, the
    /// merged identifiers the cache's synonyms column never holds, dropped.
    archive_dbsnp: bool,
    /// `identifier_indel_transcript_notations`: the `hgvsc`, `hgvsp` and
    /// `mane_select` keys of an identifier input's insertion or deletion
    /// alleles dropped, the keys the recoder leaves unshifted for those
    /// inputs alone.
    identifier_indels: bool,
}

impl Masks {
    fn of(masks: &[String]) -> Masks {
        Masks {
            prefixes: masked_prefixes(masks),
            archive_dbsnp: masks.iter().any(|m| m == "archive_dbsnp_synonyms"),
            identifier_indels: masks
                .iter()
                .any(|m| m == "identifier_indel_transcript_notations"),
        }
    }
}

/// One allele object in comparable form: each key's values as a sorted set of
/// canonical texts.
type CanonicalAllele = BTreeMap<String, BTreeSet<String>>;

/// One output object in comparable form: its input text and each allele's
/// comparable form.
type CanonicalObject = (String, BTreeMap<String, CanonicalAllele>);

/// Whether an allele object of an identifier input (no notation separator, no
/// space) is an insertion or deletion: its allele is `-`, or an `hgvsg` names
/// a deletion, duplication or insertion, or an `spdi` deletes and inserts
/// sequences of different lengths.
fn is_identifier_indel(input: &str, allele: &str, fields: &Value) -> bool {
    if input.contains(':') || input.contains(char::is_whitespace) {
        return false;
    }
    if allele == "-" {
        return true;
    }
    let texts = |key: &str| {
        fields
            .get(key)
            .and_then(Value::as_array)
            .map(|items| items.iter().filter_map(Value::as_str).collect::<Vec<_>>())
            .unwrap_or_default()
    };
    texts("hgvsg")
        .iter()
        .any(|h| h.contains("del") || h.contains("dup") || h.contains("ins"))
        || texts("spdi").iter().any(|s| {
            let parts: Vec<&str> = s.split(':').collect();
            parts.len() == 4 && parts[2].len() != parts[3].len()
        })
}

/// One allele object in comparable form: every key's values as a sorted set of
/// canonical texts, a `var_synonyms` entry with its identifiers sorted, a
/// `mane_select` triple as one text, a masked transcript entry dropped, the
/// keys a mask class drops left out.
fn canonical_allele(object: &Value, masks: &Masks, dropped_keys: &[&str]) -> CanonicalAllele {
    let is_masked = |text: &str| masks.prefixes.iter().any(|prefix| text.starts_with(prefix));
    let mut out: CanonicalAllele = BTreeMap::new();
    for (key, value) in object.as_object().unwrap() {
        if dropped_keys.contains(&key.as_str()) {
            continue;
        }
        let mut set = BTreeSet::new();
        match value {
            Value::Array(items) => {
                for item in items {
                    let text = match (key.as_str(), item) {
                        ("mane_select", Value::Object(triple)) => {
                            if triple
                                .get("hgvsc")
                                .and_then(Value::as_str)
                                .is_some_and(is_masked)
                                || triple
                                    .get("hgvsp")
                                    .and_then(Value::as_str)
                                    .is_some_and(is_masked)
                            {
                                continue;
                            }
                            serde_json::to_string(item).unwrap()
                        }
                        ("var_synonyms", Value::String(entry)) => match entry.split_once(": ") {
                            Some((source, ids)) => {
                                if masks.archive_dbsnp && source == "Archive dbSNP" {
                                    continue;
                                }
                                let mut ids: Vec<&str> = ids.split(',').collect();
                                ids.sort_unstable();
                                format!("{source}: {}", ids.join(","))
                            }
                            None => entry.clone(),
                        },
                        (_, Value::String(text)) => {
                            if is_masked(text) {
                                continue;
                            }
                            text.clone()
                        }
                        (_, other) => serde_json::to_string(other).unwrap(),
                    };
                    set.insert(text);
                }
            }
            other => {
                set.insert(serde_json::to_string(other).unwrap());
            }
        }
        // A key whose every value a mask dropped is compared as absent: the
        // other side writes no key for an empty list.
        if !set.is_empty() {
            out.insert(key.clone(), set);
        }
    }
    out
}

/// A recoder output in comparable form: one entry per object, keyed by the
/// object's input (`null` objects counted under an empty key).
fn canonical(results: &[Value], masks: &Masks) -> Vec<CanonicalObject> {
    results
        .iter()
        .map(|result| match result {
            Value::Object(object) => {
                let input = object
                    .values()
                    .find_map(|v| v.get("input").and_then(Value::as_str))
                    .unwrap_or("")
                    .to_string();
                let alleles = object
                    .iter()
                    .map(|(allele, fields)| {
                        let canonical = if allele == "warnings" {
                            BTreeMap::from([(
                                "warnings".to_string(),
                                fields
                                    .as_array()
                                    .unwrap()
                                    .iter()
                                    .map(|w| w.as_str().unwrap().to_string())
                                    .collect(),
                            )])
                        } else {
                            let dropped: &[&str] = if masks.identifier_indels
                                && is_identifier_indel(&input, allele, fields)
                            {
                                &["hgvsc", "hgvsp", "mane_select"]
                            } else {
                                &[]
                            };
                            canonical_allele(fields, masks, dropped)
                        };
                        (allele.clone(), canonical)
                    })
                    .collect();
                (input, alleles)
            }
            _ => (String::new(), BTreeMap::new()),
        })
        .collect()
}

/// Every `tool_cases` entry of the corpus manifest reproduces the recorded
/// `variant_recoder` output: the same objects in the same order, each allele's
/// keys holding the same value sets, after the manifest's masks drop the
/// RefSeq and LRG transcript entries an Ensembl cache cannot carry, the
/// `Archive dbSNP` synonyms its synonyms column never holds, and the
/// transcript notations the recoder leaves unshifted for an identifier input's
/// insertions and deletions. The expected files are written by the reference
/// generator in database mode.
#[test]
fn expected_outputs_match_variant_recoder_on_the_corpus() {
    let dir = corpus_dir();
    let manifest = manifest();
    let covers: Vec<&str> = manifest["covers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    assert_eq!(covers, ["variant_recoder", "input-data-and-id"]);
    let masks: Vec<String> = manifest["masks"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| m.as_str().unwrap().to_string())
        .collect();
    let masks = Masks::of(&masks);
    let cases = manifest["tool_cases"].as_array().unwrap();
    assert_eq!(cases.len(), 4, "{} cases", cases.len());
    let mut failures = Vec::new();
    for case in cases {
        let name = case["name"].as_str().unwrap();
        let input = case["input"].as_str().unwrap();
        let expected_path = dir.join(case["expected"].as_str().unwrap());
        let Ok(expected_text) = fs::read_to_string(&expected_path) else {
            failures.push(format!(
                "{name}: {} is not there; the reference generator writes it: \
                 scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38-recoder \
                 --tool variant_recoder",
                expected_path.display()
            ));
            continue;
        };
        let args: Vec<&str> = case["args"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| a.as_str().unwrap())
            .collect();
        let mut full = vec!["-i", input];
        full.extend(args);
        let output = recode_in(&full, None);
        if !output.status.success() {
            failures.push(format!(
                "{name}: exited {:?}: {}",
                output.status.code(),
                stderr(&output)
            ));
            continue;
        }
        let actual = canonical(&parse(&stdout(&output), name), &masks);
        let expected = canonical(&parse(&expected_text, name), &masks);
        if actual.len() != expected.len() {
            failures.push(format!(
                "{name}: {} objects, expected {}",
                actual.len(),
                expected.len()
            ));
            continue;
        }
        for ((input, actual_alleles), (expected_input, expected_alleles)) in
            actual.iter().zip(&expected)
        {
            if input != expected_input {
                failures.push(format!(
                    "{name}: object {input:?} where {expected_input:?} was expected"
                ));
                break;
            }
            if actual_alleles != expected_alleles {
                let detail = actual_alleles
                    .iter()
                    .chain(expected_alleles.iter())
                    .map(|(allele, _)| allele.clone())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .filter_map(|allele| {
                        let a = actual_alleles.get(&allele);
                        let e = expected_alleles.get(&allele);
                        (a != e)
                            .then(|| format!("  allele {allele}: actual {a:?}\n    expected {e:?}"))
                    })
                    .collect::<Vec<_>>()
                    .join("\n");
                failures.push(format!("{name}: input {input:?} differs:\n{detail}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
