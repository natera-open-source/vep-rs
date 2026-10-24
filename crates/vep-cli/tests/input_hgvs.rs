// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads HGVS genomic, transcript and protein notation:
//! the format needs the reference FASTA, each notation annotates as its VCF
//! twin does and is named by itself, a notation Ensembl VEP cannot resolve is
//! skipped with its warning, a gene name as the reference is warned about, a
//! transcript cache a notation's lookup cannot read stops the run under the
//! `cache` component, and `--max_sv_size` bounds the span a genomic notation
//! may describe. Over the same corpus, the HGVSc `--hgvs` writes for a record
//! on a chromosome the FASTA lacks.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-input-hgvs")
}

/// The binary's arguments after the input: the corpus cache, the output path
/// and `extra`, with the corpus reference under `--fasta` when `fasta` is set.
fn args(out: &Path, fasta: bool, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        corpus().join("json_cache").display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
    ];
    if fasta {
        args.push("--fasta".into());
        args.push(corpus().join("reference.fa.gz").display().to_string());
    }
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

fn vep(dir: &Path, input: &Path, out: &Path, fasta: bool, extra: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .args(args(out, fasta, extra))
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The output's lines that are not `#` headers.
fn rows(out: &Path) -> Vec<String> {
    String::from_utf8_lossy(&std::fs::read(out).unwrap())
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// The rows of a successful run on `input` with the reference and `extra` flags.
fn annotate(dir: &Path, input: &Path, out: &Path, extra: &[&str]) -> Vec<String> {
    let output = vep(dir, input, out, true, extra);
    assert!(
        output.status.success(),
        "{} {extra:?}: {}",
        input.display(),
        stderr(&output)
    );
    rows(out)
}

/// A row's `Uploaded_variation` column.
fn name(row: &str) -> &str {
    row.split_once('\t').unwrap().0
}

/// A row past its `Uploaded_variation` column.
fn after_name(row: &str) -> &str {
    row.split_once('\t').unwrap().1
}

/// HGVS input needs the reference the notations are resolved against: without
/// `--fasta` the run stops before any output, under an explicit `--format` and
/// under the detection of the format from the first line, as Ensembl VEP's own
/// parser refuses the format offline (Parser/HGVS.pm:104).
#[test]
fn hgvs_input_needs_a_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(dir.path(), "records.hgvs", "21:g.30282169C>A\n");
    for extra in [&["--format", "hgvs"][..], &[][..]] {
        let out = dir.path().join("refused.txt");
        let output = vep(dir.path(), &input, &out, false, extra);
        assert_eq!(output.status.code(), Some(1), "{extra:?}");
        assert_eq!(
            stderr(&output),
            "ERROR: [input] Cannot use HGVS format in offline mode without a FASTA file\n",
            "{extra:?}"
        );
        assert!(!out.exists(), "{extra:?}: nothing is written");
    }
}

/// The rows of a run grouped by `Uploaded_variation`, in output order.
fn rows_by_name(rows: &[String]) -> BTreeMap<&str, Vec<&str>> {
    let mut grouped: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for row in rows {
        grouped.entry(name(row)).or_default().push(row);
    }
    grouped
}

/// A row past its name with its `Allele` column reverse-complemented.
fn with_allele_complemented(row: &str) -> String {
    let mut columns: Vec<String> = after_name(row).split('\t').map(str::to_string).collect();
    columns[1] = columns[1]
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => other as char,
        })
        .collect();
    columns.join("\t")
}

/// Every notation of the corpus annotates as its VCF twin does: the twin is
/// the `variants.vcf` record of the same ordinal, named by the notation, and
/// its rows match the notation's past the name column, under `--format hgvs`
/// and under format detection. A notation on a reverse-strand transcript
/// carries its alleles on that strand (VariationFeatureAdaptor.pm:2132), so
/// its `Allele` column is the twin's reverse complement. `Uploaded_variation`
/// is the notation itself, as Ensembl VEP names the record
/// (Parser/HGVS.pm:237): the accession `NC_000021.9` and the `chr21` spelling
/// print their rows at `21:`, the reference FASTA's name for the chromosome,
/// and the `(p.=)` suffix stays in the name.
#[test]
fn hgvs_records_annotate_as_their_vcf_twins_and_are_named_by_the_notation() {
    let dir = tempfile::tempdir().unwrap();
    let vcf_rows = annotate(
        dir.path(),
        &corpus().join("variants.vcf"),
        &dir.path().join("vcf.txt"),
        &[],
    );
    let twins = rows_by_name(&vcf_rows);
    let notations: Vec<String> = std::fs::read_to_string(corpus().join("variants.hgvs.txt"))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(notations.len(), 76);
    assert_eq!(twins.len(), 76);
    for (label, extra) in [("explicit", &["--format", "hgvs"][..]), ("guess", &[][..])] {
        let hgvs_rows = annotate(
            dir.path(),
            &corpus().join("variants.hgvs.txt"),
            &dir.path().join(format!("{label}.txt")),
            extra,
        );
        let by_name = rows_by_name(&hgvs_rows);
        assert_eq!(by_name.len(), notations.len(), "{label}");
        for notation in &notations {
            let twin_rows = &twins[notation.as_str()];
            let hgvs_rows = by_name
                .get(notation.as_str())
                .unwrap_or_else(|| panic!("{label}: {notation} has no rows"));
            let same: Vec<&str> = hgvs_rows.iter().map(|row| after_name(row)).collect();
            let twin: Vec<&str> = twin_rows.iter().map(|row| after_name(row)).collect();
            if same != twin {
                let complemented: Vec<String> = hgvs_rows
                    .iter()
                    .map(|row| with_allele_complemented(row))
                    .collect();
                assert_eq!(complemented, twin, "{label}: {notation}");
            }
        }
        for notation in ["NC_000021.9:g.30282169C>A", "chr21:g.30282169C>A"] {
            let row = by_name[notation][0];
            assert_eq!(
                row.split('\t').nth(1),
                Some("21:30282169"),
                "{label}: {row}"
            );
        }
    }
}

/// A transcript notation resolves through the cache: on the reverse-strand
/// TFF3 transcript `c.5C>T` lands at 21:42315370 with the allele `T`, as
/// written on the transcript (its twin's `A`), and the `n.` notations of a
/// non-coding transcript annotate on its exons and introns (`n.135+5` is the
/// fifth intronic base past a splice donor). A gene name as the reference
/// picks the gene's likely transcript (its MANE Select) with Ensembl VEP's
/// warning (Parser/HGVS.pm:219-229), and Ensembl VEP 116.2 in database mode
/// writes the rows of `CRYAA:c.100T>C` and `TFF3:c.5C>T` as it writes those of
/// the two transcripts' own notations. A transcript the cache lacks, a
/// position before the transcript's first base and a reference allele the
/// genome contradicts are skipped with the resolver's messages, the texts
/// Ensembl VEP 116.2 writes for the same lines.
#[test]
fn transcript_notations_resolve_through_the_cache_and_warn_about_gene_names() {
    let dir = tempfile::tempdir().unwrap();
    let corpus_rows = annotate(
        dir.path(),
        &corpus().join("variants.hgvs.txt"),
        &dir.path().join("corpus.txt"),
        &["--format", "hgvs"],
    );
    let by_name = rows_by_name(&corpus_rows);
    let tff3 = &by_name["ENST00000518498:c.5C>T"];
    assert!(tff3.iter().all(|row| {
        row.split('\t').nth(1) == Some("21:42315370") && row.split('\t').nth(2) == Some("T")
    }));
    assert!(by_name["ENST00000807019:n.135+5A>G"].iter().any(|row| {
        row.contains(
            "ENST00000807019\tTranscript\tsplice_donor_5th_base_variant,intron_variant,non_coding_transcript_variant",
        )
    }));
    let input = write(
        dir.path(),
        "genes.hgvs",
        "CRYAA:c.100T>C\nTFF3:c.5C>T\nENST00000291554:c.1T>G\nENST00000291554:c.-200A>G\n\
         ENST00009999999:c.1A>G\n",
    );
    let out = dir.path().join("genes.txt");
    let rows = annotate(dir.path(), &input, &out, &["--format", "hgvs"]);
    let picked = rows_by_name(&rows);
    assert_eq!(picked.len(), 2, "{rows:?}");
    for (gene, transcript) in [
        ("CRYAA:c.100T>C", "ENST00000291554:c.100T>C"),
        ("TFF3:c.5C>T", "ENST00000518498:c.5C>T"),
    ] {
        assert_eq!(
            picked[gene]
                .iter()
                .map(|row| after_name(row))
                .collect::<Vec<_>>(),
            by_name[transcript]
                .iter()
                .map(|row| after_name(row))
                .collect::<Vec<_>>(),
            "{gene}"
        );
    }
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    let gene_warning = |gene: &str| {
        format!(
            "WARNING: Possible invalid use of gene or protein identifier '{gene}' as HGVS reference; \
             most likely transcript will be selected\n"
        )
    };
    assert!(warnings.starts_with(&gene_warning("CRYAA")), "{warnings}");
    assert!(warnings.contains(&gene_warning("TFF3")), "{warnings}");
    for text in [
        "Unable to parse HGVS notation 'ENST00000291554:c.1T>G': Reference allele extracted from \
         ENST00000291554:43169100-43169100 (A) does not match reference allele given by HGVS notation \
         ENST00000291554:c.1T>G (T)\n",
        "Unable to parse HGVS notation 'ENST00000291554:c.-200A>G': Unable to map the cDNA coordinates \
         -200--200 to genomic coordinates for Transcript ENST00000291554\n",
        "Unable to parse HGVS notation 'ENST00009999999:c.1A>G': Could not get a Transcript object for \
         'ENST00009999999'\n",
    ] {
        assert!(warnings.contains(text), "{warnings}");
    }
}

/// Copies the directory tree at `from` to `to`.
fn copy_tree(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// The cache and the cut shard: the corpus cache copied to `dir/cache` with
/// the chr21 shard holding ENST00000518498 cut to half its bytes, so its gzip
/// stream ends early, and without its checksum manifest, which is checked
/// before a shard is inflated and would refuse the cut shard first.
fn cache_with_a_cut_shard(dir: &Path) -> (PathBuf, PathBuf) {
    let cache = dir.join("cache");
    copy_tree(&corpus().join("json_cache"), &cache);
    std::fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
    let shard = cache.join("transcripts/21/42000001-43000000.json.gz");
    let bytes = std::fs::read(&shard).unwrap();
    std::fs::write(&shard, &bytes[..bytes.len() / 2]).unwrap();
    (cache, shard)
}

/// Runs `input` in the HGVS format against `cache` and asserts the run ended
/// with exit 1, one line of standard error opening with `prefix` and nothing
/// past the header on standard output. The line ends in the inflater's own
/// words for the cut, which are its to choose.
fn assert_run_refused(dir: &Path, input: &Path, cache: &Path, prefix: &str) {
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .args(["-o", "STDOUT", "--format", "hgvs", "--json_cache"])
        .arg(cache)
        .args(["--assembly", "GRCh38", "--no_stats", "--quiet", "--fasta"])
        .arg(corpus().join("reference.fa.gz"))
        .arg("--warning_file")
        .arg(dir.join("warnings.txt"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    assert!(err.starts_with(prefix), "{err}");
    assert_eq!(err.lines().count(), 1, "{err}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.lines().all(|line| line.starts_with('#')), "{stdout}");
}

/// A transcript cache a `c.` notation's lookup cannot read stops the run under
/// the `cache` component, the component of a cache whose content is unusable,
/// with the loader's refusal after the parser's text; the notation itself is
/// sound, so the `input` component does not apply. The lookup reaches the
/// shard holding the named transcript and finds its gzip stream cut.
#[test]
fn an_unreadable_transcript_cache_stops_the_run_under_the_cache_component() {
    let dir = tempfile::tempdir().unwrap();
    let (cache, shard) = cache_with_a_cut_shard(dir.path());
    let input = write(dir.path(), "tff3.hgvs", "ENST00000518498:c.5C>T\n");
    let prefix = format!(
        "ERROR: [cache] Cannot read the transcript cache for HGVS input: Failed to inflate {}: ",
        shard.display()
    );
    assert_run_refused(dir.path(), &input, &cache, &prefix);
}

/// A genomic notation on the chromosome of the cut shard, read before the `c.`
/// notation, resolves without the cache and has its chromosome loaded as the
/// reader hands it on; the load meets the same cut shard, and the run ends
/// under the same `cache` component, with the chromosome in place of the
/// parser's text since no lookup made the refusal.
#[test]
fn a_genomic_notation_before_the_transcript_notation_meets_the_same_cache_refusal() {
    let dir = tempfile::tempdir().unwrap();
    let (cache, shard) = cache_with_a_cut_shard(dir.path());
    let input = write(
        dir.path(),
        "mixed.hgvs",
        "21:g.30282169C>A\nENST00000518498:c.5C>T\n",
    );
    let prefix = format!(
        "ERROR: [cache] Cannot read the transcript cache for chr21: Failed to inflate {}: ",
        shard.display()
    );
    assert_run_refused(dir.path(), &input, &cache, &prefix);
}

/// A notation Ensembl VEP cannot resolve is skipped and warned about with the
/// resolver's message (Parser/HGVS.pm:204-216), the line quoted as every skip
/// of a parser-fed record quotes it: a reference allele the genome contradicts,
/// a repeat, a notation with no type letter. A variant the resolver builds but
/// `validate_vf` refuses (Parser.pm:564-592) is skipped likewise: an insertion
/// at one position, and an upper-case `DEL`, whose listed bases (none) stay its
/// reference allele. The notations around them annotate.
#[test]
fn unresolvable_notations_are_skipped_with_the_warning() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "mixed.hgvs",
        "21:g.30282169G>A\n21:g.30282169C>A\n21:g.30282169C[3]\n21:30282169C>A\n\
         21:g.30282172insG\n21:g.30282172DEL\n21:g.30282172del\n",
    );
    let out = dir.path().join("mixed.txt");
    let rows = annotate(dir.path(), &input, &out, &["--format", "hgvs"]);
    let names: Vec<&str> = rows.iter().map(|row| name(row)).collect();
    assert_eq!(names, ["21:g.30282169C>A", "21:g.30282172del"]);
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (21:g.30282169G>A): Unable to parse HGVS notation '21:g.30282169G>A': \
         Reference allele extracted from 21:30282169-30282169 (C) does not match reference \
         allele given by HGVS notation 21:g.30282169G>A (G)\n\
         WARNING: line 3 skipped (21:g.30282169C[3]): Unable to parse HGVS notation '21:g.30282169C[3]'\n\
         WARNING: line 4 skipped (21:30282169C>A): Unable to parse HGVS notation '21:30282169C>A': \
         HGVS type  not recognised for 21:30282169C>A\n\
         WARNING: line 5 skipped (21:g.30282172insG): Alleles look like an insertion (-/G) \
         but coordinates are not start = end + 1 (START=30282172, END=30282172) \
         for variant (21:g.30282172insG) \n\
         WARNING: line 6 skipped (21:g.30282172DEL): Length of reference allele ( length 0) \
         does not match coordinates 30282172-30282172\n"
    );
}

/// A genomic span past the sequence's last base (chromosome 21 of GRCh38 ends at
/// 46709983) is skipped with the FASTA's refusal, where Ensembl VEP 116.2 reads
/// the span from its slice padded with `N` (`VariationFeatureAdaptor.pm:2060-2062`,
/// `SequenceAdaptor.pm:148-152` and `:222-233`), so `21:g.46709984del` is `N/-` at
/// 46709984 and `21:g.46709983_46709984insA` the insertion `-/A` after the last
/// base, each an intergenic row there. The last base itself is inside the
/// sequence on both engines: `21:g.46709983del` annotates as an intergenic row.
/// `docs/divergences/spdi-input.md` states the difference.
#[test]
fn a_span_past_the_last_base_is_skipped_with_the_fastas_refusal() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "past-end.hgvs",
        "21:g.46709984del\n21:g.46709983_46709984insA\n21:g.46709983del\n",
    );
    let out = dir.path().join("past-end.txt");
    let rows = annotate(dir.path(), &input, &out, &["--format", "hgvs"]);
    let placed: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            (columns[0], columns[1], columns[6])
        })
        .collect();
    assert_eq!(
        placed,
        [("21:g.46709983del", "21:46709983", "intergenic_variant")]
    );
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (21:g.46709984del): Unable to parse HGVS notation \
         '21:g.46709984del': The reference FASTA has no bases at 21:46709984-46709984\n\
         WARNING: line 2 skipped (21:g.46709983_46709984insA): Unable to parse HGVS notation \
         '21:g.46709983_46709984insA': The reference FASTA has no bases at 21:46709983-46709984\n"
    );
}

/// A version-less accession names no FASTA record: the genomic resolver looks the
/// whole name up, in any letter case, among the records and their aliases
/// (`record_name` in vep-fasta, the synonym table's `NC_000021.9` among them), so
/// `NC_000021:g.30282169C>A` is skipped with the FASTA's refusal and writes no
/// row, while `NC_000021.9:g.30282169C>A` and its lower-case spelling
/// `nc_000021.9:g.30282169C>A` each annotate as `21:30282169 A`, named as written.
/// Ensembl VEP 116.2 fetches the slice through `fetch_by_region('chromosome',
/// $reference) || fetch_by_region(undef, $reference)`
/// (`VariationFeatureAdaptor.pm:1889`), whose synonym lookup retries the name as
/// a prefix (`SliceAdaptor.pm:2870-2884`, `"$escaped_seq_region_name%"`), so there
/// the version-less notation annotates as the versioned one does.
/// `docs/divergences/spdi-input.md` states the difference.
#[test]
fn a_version_less_accession_is_skipped_with_the_fastas_refusal() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "accession.hgvs",
        "NC_000021:g.30282169C>A\nNC_000021.9:g.30282169C>A\nnc_000021.9:g.30282169C>A\n",
    );
    let out = dir.path().join("accession.txt");
    let rows = annotate(dir.path(), &input, &out, &["--format", "hgvs"]);
    let placed: BTreeMap<&str, Vec<&str>> = rows_by_name(&rows)
        .into_iter()
        .map(|(notation, rows)| (notation, rows.iter().map(|row| after_name(row)).collect()))
        .collect();
    assert_eq!(
        placed.keys().copied().collect::<Vec<_>>(),
        ["NC_000021.9:g.30282169C>A", "nc_000021.9:g.30282169C>A"]
    );
    assert_eq!(
        placed["NC_000021.9:g.30282169C>A"],
        placed["nc_000021.9:g.30282169C>A"]
    );
    for row in &rows {
        assert!(after_name(row).starts_with("21:30282169\tA\t"), "{row}");
    }
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (NC_000021:g.30282169C>A): Unable to parse HGVS notation \
         'NC_000021:g.30282169C>A': The reference FASTA has no sequence named NC_000021\n"
    );
}

/// A notation whose alternate names no bases builds a feature with the reference
/// alone: `21:g.30282172DUP` under `--lookup_ref` keeps `get_hgvs_alleles`'s `-`
/// (`_hgvs_from_components` tests `/dup/` in lower case only,
/// VariationFeatureAdaptor.pm:2071), which the genome's `T` replaces (`:2096-2097`)
/// beside the empty alternate `dup([A-Z]*)$` captured, and `validate_vf`'s split
/// drops that empty field (Parser.pm:575), so Ensembl VEP 116.2 annotates no allele
/// and writes no row and no warning for it; `21:g.30282172_30282173ins` is the
/// same shape (`-/`). The lower-case `dup` beside them annotates the same rows
/// with and without `--lookup_ref`: its alleles are the `-` and the genome's
/// base either way (`:2076-2077`). Without `--lookup_ref` the `-` fails the
/// reference check and the notation is skipped.
#[test]
fn a_notation_naming_no_alternate_bases_writes_no_row() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "empty-alternate.hgvs",
        "21:g.30282172DUP\n21:g.30282172Dup\n21:g.30282172_30282173ins\n21:g.30282172dup\n",
    );
    let out = dir.path().join("looked-up.txt");
    let looked_up = annotate(
        dir.path(),
        &input,
        &out,
        &["--format", "hgvs", "--lookup_ref"],
    );
    assert!(!looked_up.is_empty());
    assert!(
        looked_up.iter().all(|row| name(row) == "21:g.30282172dup"),
        "{looked_up:?}"
    );
    assert!(
        !std::fs::exists(format!("{}_warnings.txt", out.display())).unwrap(),
        "no warning is written"
    );

    let out = dir.path().join("checked.txt");
    let rows = annotate(dir.path(), &input, &out, &["--format", "hgvs"]);
    assert_eq!(rows, looked_up);
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert!(
        warnings.starts_with(
            "WARNING: line 1 skipped (21:g.30282172DUP): Unable to parse HGVS notation '21:g.30282172DUP': \
             Reference allele extracted from 21:30282172-30282172 (T) does not match reference \
             allele given by HGVS notation 21:g.30282172DUP (-)\n"
        ),
        "{warnings}"
    );
}

/// Each row's `(Location, Allele, UPLOADED_ALLELE)` by its `Uploaded_variation`.
fn placed_by_name(rows: &[String]) -> BTreeMap<&str, (&str, &str, &str)> {
    rows.iter()
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            let uploaded_allele = columns
                .last()
                .unwrap()
                .split(';')
                .find_map(|field| field.strip_prefix("UPLOADED_ALLELE="))
                .unwrap_or_else(|| panic!("no UPLOADED_ALLELE in {row}"));
            (columns[0], (columns[1], columns[2], uploaded_allele))
        })
        .collect()
}

/// Under `--minimal` a notation's alleles are trimmed as written before they are
/// upper-cased (`post_process_vfs`, Parser.pm:852, then `validate_vf`, :554):
/// `21:g.30282168_30282170TCA>TCG` annotates as `A/G` at `21:30282170`, and
/// `21:g.30282168_30282170delinsGGa` keeps its `A`/`a` pair and annotates as `GGA`
/// over `21:30282168-30282170`, as Ensembl VEP 116.2 does; a trim of the
/// upper-cased pair would drop the shared `A` and annotate `GG` over two bases.
/// `UPLOADED_ALLELE` prints the pair as the notation built it, `TCA/TCG`,
/// `TCA/GGa` and `TCA/TCg`, the `original_allele_string` the trim records
/// (:1008) before the upper-casing. Without `--minimal` the same-length
/// `TCA>TCg` prints its upper-cased working pair, `TCA/TCG`, since no trim
/// records the pair as built, and the indels `TCA>TC` and `TCA>TG` print
/// `TCA/TC` and `TCA/TG`, recorded by the trim every indel goes through
/// (:862-881), which places `TCA>TG` as `G` over `21:30282169-30282170` with or
/// without the flag.
///
/// `TCA>TG`, whose first trim leaves `CA/G` with the lengths still unequal, is
/// where the engines part: Ensembl VEP's second pass (:881, `$vf =
/// ${$self->minimise_alleles([$vf])}[0] if $is_non_minimised_indel;`) reads the
/// trimmed string as the uploaded one, since only the VCF parser sets
/// `nontrimmed_allele_string` (Parser/VCF.pm:343) and the HGVS parser hands its
/// records to the same `post_process_vfs` (Parser/HGVS.pm:248), so there
/// `UPLOADED_ALLELE` prints `CA/G` (OutputFactory.pm:1360), the value Ensembl
/// VEP prints for the identical feature built from the Ensembl-format record
/// `21 30282168 30282170 TCA/TG +`; here it prints `TCA/TG`, the string as
/// parsed, and `docs/divergences/region-input.md` states the difference.
#[test]
fn minimal_trims_a_notation_in_its_written_case() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "minimal.hgvs",
        "21:g.30282168_30282170TCA>TCG\n21:g.30282168_30282170delinsGGa\n\
         21:g.30282168_30282170TCA>TCg\n21:g.30282168_30282170TCA>TG\n",
    );
    let out = dir.path().join("minimal.txt");
    let rows = annotate(
        dir.path(),
        &input,
        &out,
        &["--format", "hgvs", "--minimal", "--uploaded_allele"],
    );
    assert_eq!(
        placed_by_name(&rows),
        BTreeMap::from([
            (
                "21:g.30282168_30282170TCA>TCG",
                ("21:30282170", "G", "TCA/TCG")
            ),
            (
                "21:g.30282168_30282170delinsGGa",
                ("21:30282168-30282170", "GGA", "TCA/GGa")
            ),
            (
                "21:g.30282168_30282170TCA>TCg",
                ("21:30282170", "G", "TCA/TCg")
            ),
            (
                "21:g.30282168_30282170TCA>TG",
                ("21:30282169-30282170", "G", "TCA/TG")
            ),
        ])
    );
    assert!(
        rows.iter()
            .filter(|row| name(row) == "21:g.30282168_30282170TCA>TG")
            .all(|row| row.contains("MINIMISED=1")),
        "{rows:?}"
    );
    let input = write(
        dir.path(),
        "plain.hgvs",
        "21:g.30282168_30282170TCA>TCg\n21:g.30282168_30282170TCA>TC\n\
         21:g.30282168_30282170TCA>TG\n",
    );
    let rows = annotate(
        dir.path(),
        &input,
        &out,
        &["--format", "hgvs", "--uploaded_allele"],
    );
    assert_eq!(
        placed_by_name(&rows),
        BTreeMap::from([
            (
                "21:g.30282168_30282170TCA>TCg",
                ("21:30282168-30282170", "TCG", "TCA/TCG")
            ),
            (
                "21:g.30282168_30282170TCA>TC",
                ("21:30282170", "-", "TCA/TC")
            ),
            (
                "21:g.30282168_30282170TCA>TG",
                ("21:30282169-30282170", "G", "TCA/TG")
            ),
        ])
    );
}

/// `docs/divergences/region-input.md` states the `UPLOADED_ALLELE` difference
/// `minimal_trims_a_notation_in_its_written_case` shows on `TCA>TG`: under
/// `--minimal` Ensembl VEP prints the trimmed `CA/G` and vep-rs the parsed
/// `TCA/TG`, in every input format but VCF.
#[test]
fn region_input_page_states_the_minimal_uploaded_allele_difference() {
    let page = std::fs::read_to_string(format!(
        "{}/../../docs/divergences/region-input.md",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
    .split_whitespace()
    .collect::<Vec<_>>()
    .join(" ");
    for sentence in [
        "An Ensembl-format, HGVS, SPDI or identifier record behaves the same way: its parser hands \
         the record to the same pass (`Parser/HGVS.pm` line 248, `Parser/SPDI.pm` line 204, \
         `Parser/ID.pm` line 196, `Parser/VEP_input.pm` line 167) and sets no untrimmed string.",
        "for `21:g.30282168_30282170TCA>TG`, the same feature built from a notation, `CA/G` and \
         `TCA/TG`, both `MINIMISED=1`, and without `--minimal` both engines print `TCA/TG`",
    ] {
        assert!(page.contains(sentence), "region-input.md lacks {sentence:?}");
    }
}

/// The span a genomic notation describes must be shorter than `--max_sv_size`
/// (VariationFeatureAdaptor.pm:1886): a twenty-base deletion annotates under
/// the default limit and is skipped under `--max_sv_size 10` with Ensembl
/// VEP's message and its hint.
#[test]
fn max_sv_size_bounds_the_span_of_a_notation() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(dir.path(), "span.hgvs", "21:g.30282180_30282199del\n");
    let rows = annotate(
        dir.path(),
        &input,
        &dir.path().join("default.txt"),
        &["--format", "hgvs"],
    );
    assert!(!rows.is_empty());
    assert!(rows[0].starts_with("21:g.30282180_30282199del\t21:30282180-30282199\t-\t"));
    let out = dir.path().join("limited.txt");
    let rows = annotate(
        dir.path(),
        &input,
        &out,
        &["--format", "hgvs", "--max_sv_size", "10"],
    );
    assert!(rows.is_empty(), "{rows:?}");
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (21:g.30282180_30282199del): Unable to parse HGVS notation '21:g.30282180_30282199del': \
         Region requested must be smaller than 10. This value can be adjusted with --max_sv_size.\n"
    );
}

/// A `#` line counts toward the line numbers that follow, as Ensembl VEP counts
/// the record it reads there: ensembl-io's list parser has no metadata lines
/// (ListBasedParser.pm:60, `sub is_metadata { return; }`), so
/// `create_VariationFeatures` numbers the line (Parser/HGVS.pm:151,
/// `$self->line_number($self->line_number + 1);`) and hands its text to the
/// resolver, which refuses it (`:204-216`, `WARNING: Unable to parse HGVS
/// notation '#c'` with the resolver's message, and the next record is read);
/// here such a line is a comment, passed over without a word. A `#` line with no
/// word character is an empty line to `skip_empty_lines` (Parser.pm:273,
/// `$parser->{record} !~ /\w+/`), warned about and counted as one. So for
/// `21:g.30282169C>A`, `# c`, `21:g.30282169G>A` the refusal names line 3, and
/// for `21:g.30282169C>A`, `#`, `21:g.30282169G>A` the empty line is line 2 and
/// the refusal names line 3 (Ensembl VEP refuses the format offline,
/// Parser/HGVS.pm:104, so these are its source's numbers, not a run's).
#[test]
fn a_hash_line_counts_toward_the_line_numbers() {
    let dir = tempfile::tempdir().unwrap();
    let refusal = "WARNING: line 3 skipped (21:g.30282169G>A): Unable to parse HGVS notation \
                   '21:g.30282169G>A': Reference allele extracted from 21:30282169-30282169 (C) \
                   does not match reference allele given by HGVS notation 21:g.30282169G>A (G)\n";
    let comment = write(
        dir.path(),
        "comment.hgvs",
        "21:g.30282169C>A\n# c\n21:g.30282169G>A\n",
    );
    let out = dir.path().join("comment.txt");
    let rows = annotate(dir.path(), &comment, &out, &["--format", "hgvs"]);
    assert!(
        rows.iter().all(|row| name(row) == "21:g.30282169C>A"),
        "{rows:?}"
    );
    assert!(!rows.is_empty());
    assert_eq!(
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap(),
        refusal
    );
    let bare = write(
        dir.path(),
        "bare.hgvs",
        "21:g.30282169C>A\n#\n21:g.30282169G>A\n",
    );
    let out = dir.path().join("bare.txt");
    let bare_rows = annotate(dir.path(), &bare, &out, &["--format", "hgvs"]);
    assert_eq!(bare_rows, rows);
    assert_eq!(
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap(),
        format!("WARNING: Skipped 1 empty line from line 2\n{refusal}")
    );
}

/// Under `--hgvs` a record on a chromosome the FASTA has no sequence for gets no
/// HGVSc when its reference bases cannot be read, a substitution and a deletion,
/// while an insertion, which reads no base, keeps its notation and HGVSp is
/// unchanged on every row; no warning is written. Over the corpus reference the
/// same records print `c.764G>T`, `c.307_308insC`, `c.632del` and `c.196_197del`,
/// the notations Ensembl VEP 116.2 prints over the GRCh38 primary assembly.
/// Over a FASTA holding chromosome 22 alone Ensembl VEP 116.2 reads the
/// transcript's slice as `N` per base (`Utils/FastaSequence.pm:597-601`, `return
/// 'N' x (($end - $start) + 1)`), hands it to `hgvs_variant_notation` as the
/// reference (`TranscriptVariationAllele.pm:1474-1476`) and cuts the reference
/// allele from it (`Utils/Sequence.pm:526`), so it prints
/// `ENST00000340345.6:c.764N>T`, `ENST00000416044.1:c.307_308insC`,
/// `ENST00000360542.5:c.632del` and `ENST00000334160.6:c.196_197del` with the
/// same HGVSp and no warning. `docs/cache-setup.md` states the difference.
#[test]
fn a_contig_the_fasta_lacks_gives_no_hgvsc_for_a_substitution_or_deletion() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "stops.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30282169\tsnv_fwd\tC\tA\t.\t.\t.\n\
         21\t30289183\tstop_ins1_lost_rev\tC\tCG\t.\t.\t.\n\
         21\t30319403\tstop_del1_kept_rev\tTC\tT\t.\t.\t.\n\
         21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.\n",
    );
    // Chromosome 22 alone, 120 bases over two lines.
    let other = write(
        dir.path(),
        "other.fa",
        &format!(">22\n{}\n{}\n", "ACGT".repeat(15), "ACGT".repeat(15)),
    );
    write(dir.path(), "other.fa.fai", "22\t120\t4\t60\t61\n");
    let fields = "Uploaded_variation,Feature,HGVSc,HGVSp";
    let notations = |rows: &[String]| -> BTreeMap<String, (String, String)> {
        rows.iter()
            .map(|row| {
                let columns: Vec<&str> = row.split('\t').collect();
                (
                    format!("{} {}", columns[0], columns[1]),
                    (columns[2].to_string(), columns[3].to_string()),
                )
            })
            .collect()
    };
    let protein = [
        (
            "snv_fwd ENST00000340345",
            "ENSP00000339238.4:p.Ter255LeuextTer21",
        ),
        (
            "stop_ins1_lost_rev ENST00000416044",
            "ENSP00000398619.1:p.Ter103delinsSerTer",
        ),
        (
            "stop_del1_kept_rev ENST00000360542",
            "ENSP00000353742.3:p.Ter211%3D",
        ),
        (
            "stop_del2_lost_rev ENST00000334160",
            "ENSP00000346536.3:p.Ter66%3D",
        ),
    ];
    let expected = |coding: [&str; 4]| -> BTreeMap<String, (String, String)> {
        protein
            .iter()
            .zip(coding)
            .map(|((row, hgvsp), hgvsc)| (row.to_string(), (hgvsc.to_string(), hgvsp.to_string())))
            .collect()
    };

    let out = dir.path().join("genome.txt");
    let genome_rows = annotate(
        dir.path(),
        &input,
        &out,
        &["--hgvs", "--tab", "--fields", fields],
    );
    assert_eq!(
        notations(&genome_rows),
        expected([
            "ENST00000340345.6:c.764G>T",
            "ENST00000416044.1:c.307_308insC",
            "ENST00000360542.5:c.632del",
            "ENST00000334160.6:c.196_197del",
        ])
    );

    let out = dir.path().join("other.txt");
    let output = vep(
        dir.path(),
        &input,
        &out,
        false,
        &[
            "--fasta",
            other.to_str().unwrap(),
            "--hgvs",
            "--tab",
            "--fields",
            fields,
        ],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(
        notations(&rows(&out)),
        expected(["-", "ENST00000416044.1:c.307_308insC", "-", "-"])
    );
    assert!(
        !std::fs::exists(format!("{}_warnings.txt", out.display())).unwrap(),
        "no warning is written"
    );
}
