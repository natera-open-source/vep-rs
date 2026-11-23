// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `vep` binary says when its command asks for a part its cache does
//! not carry. The cache's own record states each part: the names under
//! `left_out` in `info.json`, which the builder and the converter write, the
//! frequency columns `variation_cols` declares, and the absence of variation
//! data. Each line is an `INFO:` notice, written once to standard error and
//! not under `--quiet`; a cache without variation data draws a warning, in
//! the warnings file as well. None changes a row of the output. `vep recode`,
//! whose annotation run is quiet, writes the two lines that run reaches
//! itself, before its own line: that of a cache converted without sequence
//! edits and, where its fields name `mane_select`, that of the column on a
//! cache built without the core table dumps. Its annotation run asks for no
//! known variant, so the tool, which reads the cache's variation files for
//! `id` and `var_synonyms` itself, writes one `INFO:` line more over a cache
//! without variation data, after the run's lines and before its own.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// The committed cache of the corpus `name`.
fn cache_of(name: &str) -> PathBuf {
    corpus(name).join("json_cache")
}

/// Copies the files and directories of `from` into the new directory `to`.
fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// A copy of the committed cache of the corpus `name` at `to` whose
/// `info.json` lists `names` under `left_out`, as the tool that made such a
/// cache writes them. The copy holds no checksum manifest, so the changed
/// metadata is read unchecked.
fn cache_leaving_out(name: &str, to: &Path, names: &[&str]) -> PathBuf {
    copy_dir(&cache_of(name), to);
    std::fs::remove_file(to.join("CHECKSUMS.sha256")).unwrap();
    let path = to.join("info.json");
    let mut info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    info["left_out"] = serde_json::json!(names);
    std::fs::write(&path, serde_json::to_string_pretty(&info).unwrap()).unwrap();
    to.to_path_buf()
}

/// One completed run of the binary.
struct Run {
    /// The lines of standard error.
    written: Vec<String>,
    /// The data rows of the output file.
    rows: Vec<String>,
    /// The lines of the warnings file, none when the run wrote no file.
    warnings: Vec<String>,
}

impl Run {
    /// How often `line` stands on standard error.
    fn count(&self, line: &str) -> usize {
        self.written
            .iter()
            .filter(|written| *written == line)
            .count()
    }

    /// Whether a line of standard error holds `words`.
    fn says(&self, words: &str) -> bool {
        self.written.iter().any(|line| line.contains(words))
    }

    /// The records the rows name, by their first column.
    fn records(&self) -> BTreeSet<&str> {
        self.rows
            .iter()
            .filter_map(|row| row.split('\t').next())
            .collect()
    }
}

/// Runs the binary over the records of the corpus `records` on `cache` with
/// `flags`, writing `<name>.txt` under `dir`; exit status 0 is required.
fn run(dir: &Path, name: &str, records: &str, cache: &Path, flags: &[&str]) -> Run {
    let out = dir.join(format!("{name}.txt"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(corpus(records).join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args(["--offline", "--force_overwrite", "--no_stats"])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(output.status.code(), Some(0), "{flags:?}: {written}");
    let lines = |text: String| text.lines().map(str::to_string).collect::<Vec<_>>();
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display()))
        .map(lines)
        .unwrap_or_default();
    let rows = lines(std::fs::read_to_string(&out).unwrap())
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .collect();
    Run {
        written: lines(written),
        rows,
        warnings,
    }
}

const HUMAN: [&str; 2] = ["--assembly", "GRCh38"];
const ZEBRAFISH: [&str; 4] = ["--species", "danio_rerio", "--assembly", "GRCz11"];

/// `flags` after the flags that name the genome.
fn with<'a>(genome: &[&'a str], flags: &[&'a str]) -> Vec<&'a str> {
    [genome, flags].concat()
}

/// The route the builder's closing lines name.
const ROUTE: &str = "; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// A cache whose record names `core_tables`, as a build without the core
/// table dumps writes it: a run under `--everything` names once the twelve
/// columns that flag asks for and the consequence term, a run under
/// `--xref_refseq` and `--uniprot` its five, and a run under `--quiet` or one
/// that asks for none of the columns prints no such line. The rows are those
/// of the committed cache, whose record names nothing.
#[test]
fn a_cache_built_without_the_core_tables_names_the_columns_a_run_asks_for() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-gtf-builder";
    let cache = cache_leaving_out(records, &dir.path().join("cache"), &["core_tables"]);
    let everything = "INFO: --everything asks for SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, \
        MANE_PLUS_CLINICAL, APPRIS, miRNA, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS \
        and GENE_PHENO, which this cache does not carry: it was built without the core table \
        dumps (vep-cache-builder --core-dumps), so the columns stay empty and \
        mature_miRNA_variant is never assigned";
    let flags = with(&HUMAN, &["--everything"]);
    let stated = run(dir.path(), "stated", records, &cache, &flags);
    assert_eq!(stated.count(everything), 1, "{:?}", stated.written);
    let committed = run(dir.path(), "committed", records, &cache_of(records), &flags);
    assert!(
        !committed.says("core table dumps"),
        "{:?}",
        committed.written
    );
    assert!(!stated.rows.is_empty());
    assert_eq!(stated.rows, committed.rows);

    let cross_references = run(
        dir.path(),
        "cross_references",
        records,
        &cache,
        &with(&HUMAN, &["--xref_refseq", "--uniprot"]),
    );
    assert_eq!(
        cross_references.count(
            "INFO: --xref_refseq and --uniprot ask for RefSeq, SWISSPROT, TREMBL, UNIPARC and \
             UNIPROT_ISOFORM, which this cache does not carry: it was built without the core \
             table dumps (vep-cache-builder --core-dumps), so the columns stay empty"
        ),
        1,
        "{:?}",
        cross_references.written
    );
    for flags in [
        &["--everything", "--quiet"][..],
        &["--biotype", "--canonical"],
    ] {
        let quiet_or_unasked = run(dir.path(), "none", records, &cache, &with(&HUMAN, flags));
        assert!(
            !quiet_or_unasked.says("core table dumps"),
            "{flags:?}: {:?}",
            quiet_or_unasked.written
        );
    }
}

/// A cache whose record names `source`, as a build of a GTF or GFF3 read as
/// the transcript model writes it: a run whose `--fields` list names
/// `SOURCE` says once that the column is empty on the cache's transcripts,
/// and a run that does not write the column says nothing. The rows are those
/// of the committed cache.
#[test]
fn a_cache_built_from_a_gtf_model_says_the_source_column_is_empty() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-gtf-builder";
    let cache = cache_leaving_out(records, &dir.path().join("cache"), &["source"]);
    let line = "INFO: --fields writes the SOURCE column, which is empty on every transcript of \
        this cache: a cache built from a GTF or GFF3 read as the transcript model names no \
        source";
    let fields = ["--tab", "--fields", "Uploaded_variation,Feature,SOURCE"];
    let stated = run(
        dir.path(),
        "stated",
        records,
        &cache,
        &with(&HUMAN, &fields),
    );
    assert_eq!(stated.count(line), 1, "{:?}", stated.written);
    let committed = run(
        dir.path(),
        "committed",
        records,
        &cache_of(records),
        &with(&HUMAN, &fields),
    );
    assert!(!committed.says("SOURCE column"), "{:?}", committed.written);
    assert!(!stated.rows.is_empty());
    assert_eq!(stated.rows, committed.rows);

    let unasked = run(
        dir.path(),
        "unasked",
        records,
        &cache,
        &with(&HUMAN, &["--tab", "--fields", "Uploaded_variation,Feature"]),
    );
    assert!(!unasked.says("SOURCE column"), "{:?}", unasked.written);
}

/// A cache whose record names `panel_frequency`, as a build that read the
/// 1000 Genomes file writes it: a run under `--af` says once how the `AF` is
/// derived, with the page that states it, and a run under `--check_existing`
/// alone, which prints no `AF`, says nothing. The rows are those of the
/// committed cache.
#[test]
fn a_cache_with_a_recovered_frequency_says_so_where_a_run_reads_it() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-native-builder";
    let cache = cache_leaving_out(records, &dir.path().join("cache"), &["panel_frequency"]);
    let line = "INFO: --af reads the AF of this cache, which is each allele's count over the \
        1000 Genomes panel, recovered from the continental frequencies, and can differ from the \
        AF of Ensembl VEP's own cache \
        (docs/cache-variation-data.md#the-1000-genomes-frequencies)";
    let flags = with(&HUMAN, &["--af"]);
    let stated = run(dir.path(), "stated", records, &cache, &flags);
    assert_eq!(stated.count(line), 1, "{:?}", stated.written);
    let committed = run(dir.path(), "committed", records, &cache_of(records), &flags);
    assert!(!committed.says("recovered from"), "{:?}", committed.written);
    assert!(stated.rows.iter().any(|row| row.contains("AF=")));
    assert_eq!(stated.rows, committed.rows);

    let unread = run(
        dir.path(),
        "unread",
        records,
        &cache,
        &with(&HUMAN, &["--check_existing"]),
    );
    assert!(!unread.says("recovered from"), "{:?}", unread.written);
}

/// A cache the builder made holds the 1000 Genomes columns and no gnomAD
/// column, and its record names the cells it holds no value for. A run that
/// asks for them draws one line that names every column that stays empty,
/// with the route, and one for `MAX_AF`, which ranges over the 1000 Genomes
/// populations alone; the same run on the committed cache, whose record names
/// nothing, draws the same two lines without the route for the gnomAD
/// columns its `variation_cols` do not declare, and writes the same rows.
#[test]
fn a_built_cache_names_the_cells_a_run_asks_for_and_what_max_af_ranges_over() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-native-builder";
    let cache = cache_leaving_out(
        records,
        &dir.path().join("cache"),
        &["gnomad_frequencies", "pubmed", "var_synonyms"],
    );
    let flags = with(
        &HUMAN,
        &["--af_gnomade", "--af_gnomadg", "--max_af", "--pubmed"],
    );
    let cells = "INFO: --af_gnomade, --af_gnomadg and --pubmed ask for the gnomAD exome \
        columns, the gnomAD genome columns and PUBMED, which this cache does not carry: the \
        columns stay empty";
    let max_af = "INFO: MAX_AF and MAX_AF_POPS (--max_af) range over the 1000 Genomes \
        populations alone on this cache, whose variation data holds no gnomAD exome or gnomAD \
        genome frequency";
    let stated = run(dir.path(), "stated", records, &cache, &flags);
    for line in [cells, max_af] {
        assert_eq!(
            stated.count(&format!("{line}{ROUTE}")),
            1,
            "{line}: {:?}",
            stated.written
        );
    }
    let committed = run(dir.path(), "committed", records, &cache_of(records), &flags);
    assert_eq!(
        committed.count(
            "INFO: --af_gnomade and --af_gnomadg ask for the gnomAD exome columns and the \
             gnomAD genome columns, which this cache does not carry: the columns stay empty"
        ),
        1,
        "{:?}",
        committed.written
    );
    assert_eq!(committed.count(max_af), 1, "{:?}", committed.written);
    assert!(!stated.rows.is_empty());
    assert_eq!(stated.rows, committed.rows);
}

/// A cache with known variants and no frequency column is a cache that holds
/// no frequencies. On the zebrafish cache `--af` says once that `AF` stays
/// empty; `--filter_common` says once that the frequency filter finds no
/// `AF` and keeps every record, and keeps each record a run without it
/// writes; and a filter that includes by frequency says it drops every
/// record, and writes none.
#[test]
fn a_cache_with_known_variants_and_no_frequency_column_says_what_the_flags_get() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCz11-zebrafish";
    let cache = cache_of(records);
    let known = run(
        dir.path(),
        "known",
        records,
        &cache,
        &with(&ZEBRAFISH, &["--check_existing"]),
    );
    assert!(
        !known.says("which this cache does not carry"),
        "{:?}",
        known.written
    );
    assert!(known.warnings.is_empty(), "{:?}", known.warnings);

    let af = run(
        dir.path(),
        "af",
        records,
        &cache,
        &with(&ZEBRAFISH, &["--af"]),
    );
    assert_eq!(
        af.count("INFO: --af asks for AF, which this cache does not carry: the column stays empty"),
        1,
        "{:?}",
        af.written
    );

    let common = run(
        dir.path(),
        "common",
        records,
        &cache,
        &with(&ZEBRAFISH, &["--filter_common"]),
    );
    assert_eq!(
        common.count(
            "INFO: --filter_common reads AF, which this cache does not carry: the frequency \
             filter finds no AF and keeps every record"
        ),
        1,
        "{:?}",
        common.written
    );
    assert!(!known.records().is_empty());
    assert_eq!(common.records(), known.records());

    let include = run(
        dir.path(),
        "include",
        records,
        &cache,
        &with(
            &ZEBRAFISH,
            &["--check_frequency", "--freq_filter", "include"],
        ),
    );
    assert_eq!(
        include.count(
            "INFO: --check_frequency reads AF, which this cache does not carry: the frequency \
             filter finds no AF and drops every record"
        ),
        1,
        "{:?}",
        include.written
    );
    assert!(include.rows.is_empty(), "{:?}", include.rows);
}

/// A cache whose record declares every frequency column the command asks
/// for, and names nothing as left out, draws none of the lines.
#[test]
fn a_cache_that_declares_the_columns_a_run_asks_for_draws_no_line() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-colocated";
    let full = run(
        dir.path(),
        "full",
        records,
        &cache_of(records),
        &with(
            &HUMAN,
            &[
                "--af",
                "--af_1kg",
                "--af_gnomade",
                "--af_gnomadg",
                "--max_af",
                "--pubmed",
                "--var_synonyms",
                "--filter_common",
            ],
        ),
    );
    for words in [
        "which this cache does not carry",
        "populations alone on this cache",
        "recovered from",
        "no variation data",
    ] {
        assert!(!full.says(words), "{words}: {:?}", full.written);
    }
    assert!(full.warnings.is_empty(), "{:?}", full.warnings);
}

/// A cache without variation data warns once, on standard error and in the
/// warnings file, naming the flag that asks for known variants and the caches
/// that hold them. Under `--filter_common` the warning says that the
/// frequency filter keeps every record, and the run writes each record a run
/// without the flag writes.
#[test]
fn a_cache_without_variation_data_names_the_flags_and_what_the_filter_does() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-identifiers";
    let cache = cache_of(records);
    let warning = |asked: &str| {
        format!(
            "WARNING: [cache] no variation data under {}: Existing_variation and the \
             frequency columns stay empty under {asked}; a cache from vep-cache-builder --auto \
             or --include-variations, or a conversion of Ensembl VEP's cache, holds it \
             (docs/cache-variation-data.md)",
            cache.display()
        )
    };
    let existing = run(
        dir.path(),
        "existing",
        records,
        &cache,
        &with(&HUMAN, &["--check_existing"]),
    );
    assert_eq!(
        existing.count(&warning("--check_existing")),
        1,
        "{:?}",
        existing.written
    );
    assert_eq!(existing.warnings, [warning("--check_existing")]);

    let common = run(
        dir.path(),
        "common",
        records,
        &cache,
        &with(&HUMAN, &["--filter_common"]),
    );
    let filtered = warning("--filter_common, and the frequency filter keeps every record");
    assert_eq!(common.count(&filtered), 1, "{:?}", common.written);
    assert_eq!(common.warnings, [filtered]);
    assert!(!existing.records().is_empty());
    assert_eq!(common.records(), existing.records());
}

/// A copy of the committed cache of the corpus `name` at `to` without its
/// variation files, the `variation_cols` of its record left as they are. The
/// copy holds no checksum manifest, so the changed tree is read unchecked.
fn cache_without_variation_files(name: &str, to: &Path) -> PathBuf {
    copy_dir(&cache_of(name), to);
    std::fs::remove_file(to.join("CHECKSUMS.sha256")).unwrap();
    std::fs::remove_dir_all(to.join("variations")).unwrap();
    to.to_path_buf()
}

/// A cache whose record declares variation columns and that holds no
/// variation file warns as a cache without variation data does, and its
/// frequency filter runs on the declared columns with no known variant to
/// read. A filter that includes by frequency drops every record there: the
/// warning names that, on standard error and in the warnings file, and the
/// run writes no row. Under `--filter_common` the warning names a filter that
/// keeps every record, and the run writes each record a run without the flag
/// writes.
#[test]
fn a_cache_that_declares_columns_and_holds_no_variation_file_names_what_the_filter_does() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-colocated";
    let cache = cache_without_variation_files(records, &dir.path().join("cache"));
    let warning = |asked: &str| {
        format!(
            "WARNING: [cache] no variation data under {}: Existing_variation and the \
             frequency columns stay empty under {asked}; a cache from vep-cache-builder --auto \
             or --include-variations, or a conversion of Ensembl VEP's cache, holds it \
             (docs/cache-variation-data.md)",
            cache.display()
        )
    };
    let existing = run(
        dir.path(),
        "existing",
        records,
        &cache,
        &with(&HUMAN, &["--check_existing"]),
    );
    assert_eq!(
        existing.count(&warning("--check_existing")),
        1,
        "{:?}",
        existing.written
    );

    let include = run(
        dir.path(),
        "include",
        records,
        &cache,
        &with(&HUMAN, &["--check_frequency", "--freq_filter", "include"]),
    );
    let dropped = warning("--check_frequency, and the frequency filter drops every record");
    assert_eq!(include.count(&dropped), 1, "{:?}", include.written);
    assert!(
        include.warnings.contains(&dropped),
        "{:?}",
        include.warnings
    );
    assert!(!include.says("keeps every record"), "{:?}", include.written);
    assert!(include.rows.is_empty(), "{:?}", include.rows);

    let common = run(
        dir.path(),
        "common",
        records,
        &cache,
        &with(&HUMAN, &["--filter_common"]),
    );
    let kept = warning("--filter_common, and the frequency filter keeps every record");
    assert_eq!(common.count(&kept), 1, "{:?}", common.written);
    assert!(!existing.records().is_empty());
    assert_eq!(common.records(), existing.records());
}

/// The line of a cache converted from transcript files that hold no sequence
/// edits.
const NO_SEQ_EDITS: &str = "INFO: this cache was converted from transcript files that hold no \
    sequence edits (seq_edits): a variant at codon 1 of a transcript with an alternative \
    initiation codon is called from the cached peptide; a conversion of Ensembl VEP's cache \
    with vep-storable-to-json carries the edits (docs/cache-conversion.md)";

/// The line a completed `vep recode` run writes last over a cache read as
/// the Ensembl transcript set.
const RECODE_ENSEMBL_SET: &str = "INFO: this cache is read as the Ensembl transcript set, and \
    RefSeq notation is written from a merged cache; no variant is mapped to an LRG, which needs \
    the Ensembl database";

/// A cache whose record names `seq_edits`, as a conversion of transcript
/// files without the list writes it: every run says so once, whatever its
/// flags, and not under `--quiet`; the rows are those of the committed
/// cache. `vep recode` writes the line first and its own last. The cache
/// holds no variation data, so the tool's line of the unwritten `id` stands
/// between them, and the three are its whole standard error.
#[test]
fn a_cache_converted_without_sequence_edits_says_so_in_every_run() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-identifiers";
    let cache = cache_leaving_out(records, &dir.path().join("cache"), &["seq_edits"]);
    let stated = run(dir.path(), "stated", records, &cache, &HUMAN);
    assert_eq!(stated.count(NO_SEQ_EDITS), 1, "{:?}", stated.written);
    let committed = run(dir.path(), "committed", records, &cache_of(records), &HUMAN);
    assert!(!committed.says("seq_edits"), "{:?}", committed.written);
    assert!(!stated.rows.is_empty());
    assert_eq!(stated.rows, committed.rows);
    let quiet = run(
        dir.path(),
        "quiet",
        records,
        &cache,
        &with(&HUMAN, &["--quiet"]),
    );
    assert!(!quiet.says("seq_edits"), "{:?}", quiet.written);

    let (written, recoded) = recode(records, &cache, &["--fields", "id"]);
    let no_id = recode_without_variation_data(&cache, ID_UNWRITTEN);
    assert_eq!(
        written,
        format!("{NO_SEQ_EDITS}\n{no_id}\n{RECODE_ENSEMBL_SET}\n")
    );
    let (written, recoded_on_the_committed_cache) =
        recode(records, &cache_of(records), &["--fields", "id"]);
    let no_id = recode_without_variation_data(&cache_of(records), ID_UNWRITTEN);
    assert_eq!(written, format!("{no_id}\n{RECODE_ENSEMBL_SET}\n"));
    assert!(recoded == recoded_on_the_committed_cache);
}

/// One completed `vep recode` run over the records of the corpus `records`
/// on `cache` with `flags`, under a home of the run's own: its whole
/// standard error and its output.
fn recode(records: &str, cache: &Path, flags: &[&str]) -> (String, Vec<u8>) {
    let home = tempfile::tempdir().unwrap();
    let output = common::isolated_vep(home.path())
        .arg("recode")
        .arg("-i")
        .arg(corpus(records).join("variants.vcf"))
        .arg("--json_cache")
        .arg(cache)
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(output.status.code(), Some(0), "{flags:?}: {written}");
    (written, output.stdout)
}

/// `vep recode` hands `--mane_select` to its annotation run where its fields
/// name `mane_select`, the one flag of that run that asks for a column a
/// cache built without the core table dumps does not carry. On a cache whose
/// record names `core_tables` the tool writes the line of that column first
/// and with `seq_edits` named too, both lines in the order `vep` prints them.
/// Fields that hand the run no such flag draw no line of that column, and the
/// output is that of the committed cache. The fields name `id` and the cache
/// holds no variation data, so each run ends on the tool's line of the
/// unwritten `id` and then its own, and the lines named are its whole
/// standard error.
#[test]
fn vep_recode_names_the_mane_select_column_of_a_cache_built_without_the_core_tables() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-identifiers";
    let fasta = corpus(records).join("reference.fa.gz");
    let mane = [
        "--fields",
        "id,mane_select",
        "--fasta",
        fasta.to_str().unwrap(),
    ];
    let no_mane_select = "INFO: --mane_select asks for MANE_SELECT, which this cache does not \
        carry: it was built without the core table dumps (vep-cache-builder --core-dumps), so \
        the column stays empty";

    let cache = cache_leaving_out(records, &dir.path().join("tables"), &["core_tables"]);
    let no_id = recode_without_variation_data(&cache, ID_UNWRITTEN);
    let (written, recoded) = recode(records, &cache, &mane);
    assert_eq!(
        written,
        format!("{no_mane_select}\n{no_id}\n{RECODE_ENSEMBL_SET}\n")
    );
    let (written, recoded_on_the_committed_cache) = recode(records, &cache_of(records), &mane);
    let no_id_on_the_committed_cache =
        recode_without_variation_data(&cache_of(records), ID_UNWRITTEN);
    assert_eq!(
        written,
        format!("{no_id_on_the_committed_cache}\n{RECODE_ENSEMBL_SET}\n")
    );
    assert!(recoded == recoded_on_the_committed_cache);
    let (written, _) = recode(records, &cache, &["--fields", "id"]);
    assert_eq!(written, format!("{no_id}\n{RECODE_ENSEMBL_SET}\n"));

    let both = cache_leaving_out(
        records,
        &dir.path().join("both"),
        &["seq_edits", "core_tables"],
    );
    let (written, _) = recode(records, &both, &mane);
    let no_id = recode_without_variation_data(&both, ID_UNWRITTEN);
    assert_eq!(
        written,
        format!("{no_mane_select}\n{NO_SEQ_EDITS}\n{no_id}\n{RECODE_ENSEMBL_SET}\n")
    );
}

/// The line a completed `vep recode` run writes last where it writes
/// `var_synonyms`, over a cache read as the Ensembl transcript set.
const RECODE_ENSEMBL_SET_WITH_SYNONYMS: &str = "INFO: this cache is read as the Ensembl \
    transcript set, and RefSeq notation is written from a merged cache; no variant is mapped to \
    an LRG and no Archive dbSNP synonym is written, which need the Ensembl database";

/// What the line of [`recode_without_variation_data`] says of a run whose
/// fields read known variants for `id` alone.
const ID_UNWRITTEN: &str = "id, read from a cache's known variants, stays unwritten";

/// The line `vep recode` writes over `cache`, which holds no variation data,
/// where `unwritten` names the fields that read known variants and what
/// becomes of them.
fn recode_without_variation_data(cache: &Path, unwritten: &str) -> String {
    format!(
        "INFO: no variation data under {}: {unwritten}; a cache from vep-cache-builder --auto \
         or --include-variations, or a conversion of Ensembl VEP's cache, holds it \
         (docs/cache-variation-data.md)",
        cache.display()
    )
}

/// `vep recode` reads the cache's variation files itself for the fields `id`,
/// one of its default fields, and `var_synonyms`, and its annotation run asks
/// for no known variant, so that run raises no warning over a cache without
/// variation data. The tool writes one line there, after the lines of its
/// annotation run (`a_cache_converted_without_sequence_edits_says_so_in_every_run`)
/// and before its own, the lines its whole standard error: it names the
/// fields of the run that read known variants, says they stay unwritten and
/// names the caches that hold the data. A run whose fields read no known
/// variant, and a run over a cache that holds variation data, draw no such
/// line. A cache whose record declares variation columns and that holds no
/// variation file draws it as a cache without the columns does, and no
/// object of such a run holds either key.
#[test]
fn vep_recode_names_the_fields_a_cache_without_variation_data_leaves_unwritten() {
    let dir = tempfile::tempdir().unwrap();
    let records = "GRCh38-identifiers";
    let cache = cache_of(records);
    let id = recode_without_variation_data(&cache, ID_UNWRITTEN);
    let (written, recoded) = recode(records, &cache, &["--fields", "id"]);
    assert_eq!(written, format!("{id}\n{RECODE_ENSEMBL_SET}\n"));
    let holds = |output: &[u8], key: &str| String::from_utf8_lossy(output).contains(key);
    assert!(!holds(&recoded, "\"id\""));

    let fasta = corpus(records).join("reference.fa.gz");
    let (written, _) = recode(records, &cache, &["--fasta", fasta.to_str().unwrap()]);
    assert_eq!(written, format!("{id}\n{RECODE_ENSEMBL_SET}\n"));

    let both = recode_without_variation_data(
        &cache,
        "id and var_synonyms, read from a cache's known variants, stay unwritten",
    );
    let (written, recoded) = recode(records, &cache, &["--fields", "id", "--var_synonyms"]);
    assert_eq!(
        written,
        format!("{both}\n{RECODE_ENSEMBL_SET_WITH_SYNONYMS}\n")
    );
    assert!(!holds(&recoded, "\"id\"") && !holds(&recoded, "\"var_synonyms\""));

    let synonyms = recode_without_variation_data(
        &cache,
        "var_synonyms, read from a cache's known variants, stays unwritten",
    );
    let (written, _) = recode(records, &cache, &["--fields", "var_synonyms"]);
    assert_eq!(
        written,
        format!("{synonyms}\n{RECODE_ENSEMBL_SET_WITH_SYNONYMS}\n")
    );

    let (written, _) = recode(records, &cache, &["--fields", "spdi"]);
    assert_eq!(written, format!("{RECODE_ENSEMBL_SET}\n"));

    let known = "GRCh38-colocated";
    let (written, recoded) = recode(known, &cache_of(known), &["--fields", "id"]);
    assert_eq!(written, format!("{RECODE_ENSEMBL_SET}\n"));
    assert!(holds(&recoded, "\"id\""));
    let (written, recoded) = recode(
        known,
        &cache_of(known),
        &["--fields", "id", "--var_synonyms"],
    );
    assert_eq!(written, format!("{RECODE_ENSEMBL_SET_WITH_SYNONYMS}\n"));
    assert!(holds(&recoded, "\"id\"") && holds(&recoded, "\"var_synonyms\""));

    let emptied = cache_without_variation_files(known, &dir.path().join("emptied"));
    let both_on_the_copy = recode_without_variation_data(
        &emptied,
        "id and var_synonyms, read from a cache's known variants, stay unwritten",
    );
    let (written, recoded) = recode(known, &emptied, &["--fields", "id", "--var_synonyms"]);
    assert_eq!(
        written,
        format!("{both_on_the_copy}\n{RECODE_ENSEMBL_SET_WITH_SYNONYMS}\n")
    );
    assert!(!holds(&recoded, "\"id\"") && !holds(&recoded, "\"var_synonyms\""));
}

/// Asserts that the page `path` under `docs/` states every phrase, the
/// page's white space read as single blanks.
fn assert_page_states(path: &str, phrases: &[&str]) {
    let page = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs")
        .join(path);
    let text = std::fs::read_to_string(&page)
        .unwrap_or_else(|error| panic!("{}: {error}", page.display()))
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    for phrase in phrases {
        assert!(text.contains(phrase), "docs/{path} lacks {phrase:?}");
    }
}

/// `attribute-mirna`, `uniprot`, `symbol-source-hgnc-id-columns` and
/// `xref_refseq`: the page that lists the columns a build without the core
/// table dumps leaves empty says that the build and the run name them.
#[test]
fn the_setup_page_says_the_build_and_the_run_name_the_columns_without_the_core_tables() {
    assert_page_states(
        "cache-setup.md",
        &[
            "`TSL` is empty on GRCh37 with or without it: that gene set carries no level. The \
             build's closing lines name the columns, as does a run that asks for one.",
        ],
    );
}

/// `gff-gtf-transcript-source`, `custom-output-fields` and `attribute-mirna`:
/// the page of a GTF or GFF3 read as the transcript model says that the
/// build's last lines name the missing source and count the translations
/// that name no protein, that a run says the `SOURCE` column is empty, and
/// that the build and the run name the mature miRNA term.
#[test]
fn the_gtf_page_says_the_build_and_the_run_name_the_source_and_the_mirna_term() {
    assert_page_states(
        "cache-from-gtf-or-gff3.md",
        &[
            "The build's last lines say that the cache names no source and count the \
             transcripts [whose translation names no \
             protein](#the-translation-and-attribute-values), and a run that writes the \
             `SOURCE` column says it is empty there.",
            "The build's closing lines name the term, as does a run that asks for the `miRNA` \
             column.",
        ],
    );
}

/// `af-1kg-columns` and `api-install-options`: the page of a native cache's
/// variation data says that the build and the run say the `AF` is
/// recovered, what `MAX_AF` ranges over, that a run names the columns its
/// command asks for and the cache does not carry, reading the `left_out` key
/// of `info.json`, what a run and `vep recode` say on a cache without
/// variation data, and what a run says on one without a frequency column.
#[test]
fn the_variation_page_says_what_a_run_names_on_a_native_cache() {
    assert_page_states(
        "cache-variation-data.md",
        &[
            "A run that asks such a cache for known variants warns once, naming the flags \
             that ask, and says that a frequency filter then keeps every record.",
            "`vep recode` over such a cache names in one `INFO:` line the fields among `id` \
             and `var_synonyms` that stay unwritten.",
            "The cache's `info.json` lists under `left_out` what its build left out. A run \
             names once, in one line, the columns its command asks for and the cache does not \
             carry: those names, and any frequency column `variation_cols` does not declare. \
             A cache whose `info.json` has no `left_out` key draws no line for a name.",
            "The build's last lines say that `AF` is recovered, as does a run that reads it.",
            "`MAX_AF` and `MAX_AF_POPS` range over the populations a cache holds, on a native \
             cache the 1000 Genomes populations alone, which a run under `--max_af` says in \
             one line.",
            "On a cache with known variants and no frequency column, a native GRCh37 cache \
             among them, the line of the empty columns also says what a frequency filter then \
             does.",
        ],
    );
}

/// `vefc-seq-edits`: the conversion page says that the converter and the run
/// name a cache without sequence edits, that a cache whose record has no
/// such name draws no line, and which native caches hold the edits.
#[test]
fn the_conversion_page_says_the_converter_and_the_run_name_a_cache_without_edits() {
    assert_page_states(
        "cache-conversion.md",
        &[
            "`vep-cache-converter`, copying transcript files without the key, says so in its \
             last line and names `seq_edits` under `left_out` in the cache's `info.json`.",
            "A run on that cache then says so once, as `vep recode` does; a cache whose \
             `info.json` has no such name draws no line.",
        ],
    );
}
