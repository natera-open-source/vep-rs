// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse Ensembl GFF3 files into gene -> transcript -> exon/CDS hierarchy.
//!
//! The attributes read here are the ones the Ensembl core API writes on each
//! line (`Feature.pm`, `Gene.pm`, `Transcript.pm`, `Exon.pm` and `CDS.pm`
//! `summary_as_hash`): the feature `version`, the transcript `tag` list,
//! `ccdsid` and `transcript_support_level`, the exon `ensembl_phase` and
//! `ensembl_end_phase`, and the CDS `protein_id` and translation `version`.
//! A gene is any line whose `ID` is `gene:<id>` and a transcript any line
//! whose `ID` is `transcript:<id>`, whatever the Sequence Ontology term in
//! the type column, because the term follows the biotype (`gene`,
//! `ncRNA_gene`, `pseudogene`; `mRNA`, `lnc_RNA`, `pseudogenic_transcript`
//! and the rest) and the GRCh37 release-87 file uses `pseudogene` for both
//! a gene and a transcript line.
//!
//! Two transcripts never reach the cache, as Ensembl VEP's own dumper drops
//! them (`AnnotationSource/Database/Transcript.pm` `get_features_by_regions_uncached`):
//! a transcript of biotype `artifact` and a readthrough transcript, tagged
//! `readthrough_tra` by the core API's attribute code and
//! `readthrough_transcript` on the release-116 GRCh38 file.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use tracing::{debug, info, warn};
use vep_core::coordinate::is_cache_path_component;

use crate::core_dumps::BoundedLines;

/// The longest GFF3 line read, 16 MiB. A feature line's widest column is its
/// attributes, a few kilobytes even for a gene with every cross-reference
/// listed, so the cap refuses a file with no line end in it long before the
/// line buffer grows to its size.
const MAX_LINE_BYTES: usize = 16 << 20;

/// A parsed gene record with its child transcripts.
#[derive(Debug)]
#[allow(dead_code)]
pub struct GeneRecord {
    pub gene_id: String,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    /// The gene's `Name`, its display xref's label.
    pub name: Option<String>,
    /// The gene line's `version`.
    pub version: Option<u32>,
    pub transcripts: Vec<TranscriptRecord>,
}

/// A parsed transcript record with its child exons and CDS regions.
#[derive(Debug)]
pub struct TranscriptRecord {
    pub transcript_id: String,
    pub gene_id: String,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    pub source: String,
    pub version: Option<u32>,
    /// `Ensembl_canonical` in `tag`.
    pub is_canonical: bool,
    /// `MANE_Select` in `tag`: the membership only. The RefSeq accession is
    /// the `MANE_Select` transcript attribute, which the GFF3 does not carry.
    pub mane_select: bool,
    /// `MANE_Plus_Clinical` in `tag`, the membership only.
    pub mane_plus_clinical: bool,
    /// `gencode_basic` in `tag`, or `basic`, the GRCh37 release-87 spelling.
    pub gencode_basic: bool,
    /// `gencode_primary` in `tag`.
    pub gencode_primary: bool,
    /// `transcript_support_level`: the `TSL` attribute's value with its `tsl`
    /// prefix removed (`Transcript.pm` `tsl`), so `1`, `NA` or
    /// `5 (assigned to previous version 6)`.
    pub tsl: Option<String>,
    /// `ccdsid`, the first CCDS xref's display label.
    pub ccds: Option<String>,
    pub exons: Vec<ExonRecord>,
    pub cds_regions: Vec<CdsRecord>,
}

/// A parsed exon record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ExonRecord {
    pub exon_id: Option<String>,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub rank: Option<u32>,
    /// `ensembl_phase`, the exon's own phase (`Exon.pm` `summary_as_hash`):
    /// -1 outside the coding sequence.
    pub ensembl_phase: Option<i8>,
    /// `ensembl_end_phase`.
    pub ensembl_end_phase: Option<i8>,
}

/// A parsed CDS record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CdsRecord {
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub phase: i8,
    /// `protein_id`, the translation's stable identifier (`CDS.pm`
    /// `summary_as_hash`).
    pub protein_id: Option<String>,
    /// `version`, the translation's version; the GRCh37 release-87 file
    /// carries none.
    pub protein_version: Option<u32>,
}

/// Whether Ensembl VEP's dumper writes a transcript to its cache: it skips
/// one of biotype `artifact` and one carrying a `readthrough_tra` attribute
/// (`AnnotationSource/Database/Transcript.pm`).
fn is_dumped(biotype: &str, tags: &[&str]) -> bool {
    biotype != "artifact"
        && !tags
            .iter()
            .any(|t| *t == "readthrough_tra" || *t == "readthrough_transcript")
}

/// Parse a GFF3 file (optionally gzipped) into a gene map.
///
/// Returns `HashMap<gene_id, GeneRecord>`. Fails when the file cannot be
/// opened or read, on a line longer than `MAX_LINE_BYTES`, naming the file and
/// the line, on a `start` or `end` column that is not a number, naming the
/// column and the value, on a start of 0 or an end before its start, naming
/// the line, and on a seqid that cannot name a chromosome's cache directory
/// ([`is_cache_path_component`]), naming the value and its line.
pub fn parse_gff3(path: &Path) -> Result<HashMap<String, GeneRecord>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut genes: HashMap<String, GeneRecord> = HashMap::new();
    let mut transcript_map: HashMap<String, TranscriptRecord> = HashMap::new();
    let mut exons_by_transcript: HashMap<String, Vec<ExonRecord>> = HashMap::new();
    let mut cds_by_transcript: HashMap<String, Vec<CdsRecord>> = HashMap::new();
    let mut dropped = 0usize;

    for (line_index, line) in
        BoundedLines::new(reader, MAX_LINE_BYTES, path.display().to_string()).enumerate()
    {
        let line = String::from_utf8(line?).context("Failed to read GFF3 line")?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chr = normalize_chr(fields[0]);
        if !is_cache_path_component(chr) {
            bail!(
                "invalid chromosome name '{}' at line {}",
                fields[0],
                line_index + 1
            );
        }
        let feature_type = fields[2];
        let start = parse_coordinate("start", fields[3])?;
        let end = parse_coordinate("end", fields[4])?;
        // Region binning subtracts one from a start and spans run from start to
        // end, so a start of 0 or an end before its start has no place to go.
        if start == 0 || end < start {
            bail!(
                "GFF3 coordinates {start}-{end} at line {} are not 1-based and ordered",
                line_index + 1
            );
        }
        let strand: i8 = match fields[6] {
            "-" => -1,
            _ => 1,
        };
        let phase: i8 = fields[7].parse().unwrap_or(-1);
        let attrs = parse_attributes(fields[8]);
        let id = attrs.get("ID").map(String::as_str).unwrap_or("");

        if let Some(gene_id) = id.strip_prefix("gene:") {
            genes.insert(
                gene_id.to_string(),
                GeneRecord {
                    gene_id: gene_id.to_string(),
                    chr: chr.to_string(),
                    start,
                    end,
                    strand,
                    biotype: attrs.get("biotype").cloned().unwrap_or_default(),
                    name: attrs.get("Name").cloned(),
                    version: parse_version(&attrs),
                    transcripts: Vec::new(),
                },
            );
        } else if let Some(transcript_id) = id.strip_prefix("transcript:") {
            let biotype = attrs.get("biotype").cloned().unwrap_or_default();
            let tags: Vec<&str> = attrs
                .get("tag")
                .map(|t| t.split(',').map(str::trim).collect())
                .unwrap_or_default();
            if !is_dumped(&biotype, &tags) {
                dropped += 1;
                continue;
            }
            let parent = attrs
                .get("Parent")
                .map(|p| strip_prefix(p, "gene:"))
                .unwrap_or_default();
            transcript_map.insert(
                transcript_id.to_string(),
                TranscriptRecord {
                    transcript_id: transcript_id.to_string(),
                    gene_id: parent,
                    chr: chr.to_string(),
                    start,
                    end,
                    strand,
                    biotype,
                    source: fields[1].to_string(),
                    version: parse_version(&attrs),
                    is_canonical: tags.contains(&"Ensembl_canonical"),
                    mane_select: tags.contains(&"MANE_Select"),
                    mane_plus_clinical: tags.contains(&"MANE_Plus_Clinical"),
                    gencode_basic: tags.contains(&"gencode_basic") || tags.contains(&"basic"),
                    gencode_primary: tags.contains(&"gencode_primary"),
                    tsl: attrs.get("transcript_support_level").cloned(),
                    ccds: attrs.get("ccdsid").cloned(),
                    exons: Vec::new(),
                    cds_regions: Vec::new(),
                },
            );
        } else if feature_type == "exon" {
            let exon = ExonRecord {
                exon_id: attrs.get("Name").cloned(),
                start,
                end,
                strand,
                rank: attrs.get("rank").and_then(|r| r.parse::<u32>().ok()),
                ensembl_phase: attrs
                    .get("ensembl_phase")
                    .and_then(|p| p.parse::<i8>().ok()),
                ensembl_end_phase: attrs
                    .get("ensembl_end_phase")
                    .and_then(|p| p.parse::<i8>().ok()),
            };
            for parent in parent_transcripts(&attrs) {
                exons_by_transcript
                    .entry(parent)
                    .or_default()
                    .push(exon.clone());
            }
        } else if feature_type == "CDS" {
            let cds = CdsRecord {
                start,
                end,
                strand,
                phase,
                protein_id: attrs.get("protein_id").cloned(),
                protein_version: parse_version(&attrs),
            };
            for parent in parent_transcripts(&attrs) {
                cds_by_transcript
                    .entry(parent)
                    .or_default()
                    .push(cds.clone());
            }
        } else {
            debug!("Skipping feature type: {}", feature_type);
        }
    }

    for (tid, mut transcript) in transcript_map {
        if let Some(mut exons) = exons_by_transcript.remove(&tid) {
            // A rank compared with a start is no order at all, so the two keys
            // never share one comparator.
            if exons.iter().all(|e| e.rank.is_some()) {
                exons.sort_by_key(|e| e.rank);
            } else {
                if exons.iter().any(|e| e.rank.is_some()) {
                    warn!("Transcript {tid}: exons with and without a rank are ordered by start");
                }
                exons.sort_by_key(|e| e.start);
            }
            transcript.exons = exons;
        }
        if let Some(mut cds) = cds_by_transcript.remove(&tid) {
            cds.sort_by_key(|c| c.start);
            transcript.cds_regions = cds;
        }

        if let Some(gene) = genes.get_mut(&transcript.gene_id) {
            gene.transcripts.push(transcript);
        } else {
            debug!(
                "Transcript {} has no gene line for {}",
                tid, transcript.gene_id
            );
        }
    }

    info!(
        "Dropped {} artifact or readthrough transcripts, as the Ensembl VEP dumper does",
        dropped
    );

    Ok(genes)
}

/// Normalize chromosome names: remove "chr" prefix, "chrM" -> "MT".
fn normalize_chr(chr: &str) -> &str {
    let stripped = chr.strip_prefix("chr").unwrap_or(chr);
    if stripped == "M" {
        "MT"
    } else {
        stripped
    }
}

/// The `start` or `end` column as a number; a value that is not one is an
/// error naming the column and the value.
fn parse_coordinate(column: &str, value: &str) -> Result<u64> {
    value
        .parse()
        .with_context(|| format!("GFF3 {column} {value:?} is not a number"))
}

/// Strip a known prefix (e.g., "gene:", "transcript:") from an ID.
fn strip_prefix(s: &str, prefix: &str) -> String {
    s.strip_prefix(prefix).unwrap_or(s).to_string()
}

/// The `version` attribute as a number.
fn parse_version(attrs: &HashMap<String, String>) -> Option<u32> {
    attrs.get("version").and_then(|v| v.parse::<u32>().ok())
}

/// The transcripts an exon or CDS line belongs to: its comma-separated
/// `Parent` list without the `transcript:` prefix.
fn parent_transcripts(attrs: &HashMap<String, String>) -> Vec<String> {
    attrs
        .get("Parent")
        .map(|p| {
            p.split(',')
                .map(|s| strip_prefix(s, "transcript:"))
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

/// Parse GFF3 attribute column (key=value pairs separated by semicolons).
fn parse_attributes(attrs_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in attrs_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        if let Some((key, value)) = pair.split_once('=') {
            let decoded = value
                .replace("%20", " ")
                .replace("%2C", ",")
                .replace("%3B", ";")
                .replace("%3D", "=")
                .replace("%26", "&")
                .replace("%25", "%");
            map.insert(key.to_string(), decoded);
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parses `lines` (GFF3 body lines, no header) from a temporary file.
    fn try_parse(lines: &str) -> Result<HashMap<String, GeneRecord>> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gff3");
        std::fs::write(&path, lines).unwrap();
        parse_gff3(&path)
    }

    fn parse(lines: &str) -> HashMap<String, GeneRecord> {
        try_parse(lines).unwrap()
    }

    fn transcript<'a>(genes: &'a HashMap<String, GeneRecord>, id: &str) -> &'a TranscriptRecord {
        genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .find(|t| t.transcript_id == id)
            .unwrap_or_else(|| panic!("transcript {id} not parsed"))
    }

    // Lines of Homo_sapiens.GRCh38.116.gff3.gz, chromosome 21 (the transcript `Name`
    // attribute, which is not read, left out).
    const TPTE_GENE: &str = "21\tensembl_havana\tgene\t10521553\t10606140\t.\t+\t.\tID=gene:ENSG00000274391;Name=TPTE;biotype=protein_coding;description=transmembrane phosphatase with tensin homology [Source:HGNC Symbol%3BAcc:HGNC:12023];gene_id=ENSG00000274391;logic_name=ensembl_havana_gene_homo_sapiens;version=6\n";
    const TPTE_205: &str = "21\tensembl_havana\tmRNA\t10521553\t10606140\t.\t+\t.\tID=transcript:ENST00000622113;Parent=gene:ENSG00000274391;biotype=protein_coding;ccdsid=CCDS74772.1;tag=gencode_basic,gencode_primary;transcript_id=ENST00000622113;transcript_support_level=1;version=4\n";
    const TPTE_205_CDS: &str = "21\tensembl_havana\tCDS\t10538724\t10538734\t.\t+\t0\tID=CDS:ENSP00000482040;Parent=transcript:ENST00000622113;protein_id=ENSP00000482040;version=1\n\
21\tensembl_havana\tCDS\t10541112\t10541165\t.\t+\t1\tID=CDS:ENSP00000482040;Parent=transcript:ENST00000622113;protein_id=ENSP00000482040;version=1\n";

    // Lines of Homo_sapiens.GRCh37.87.gff3.gz, chromosome 21: a `pseudogene`
    // gene line and a `pseudogene` transcript line, the transcript tagged `basic`.
    const GRCH37_PSEUDOGENE: &str = "21\thavana\tpseudogene\t10872119\t10872429\t.\t+\t.\tID=gene:ENSG00000235476;Name=AF254982.2;biotype=pseudogene;gene_id=ENSG00000235476;logic_name=havana;version=1\n\
21\thavana\tpseudogene\t10872119\t10872429\t.\t+\t.\tID=transcript:ENST00000411878;Parent=gene:ENSG00000235476;Name=AF254982.2-001;biotype=unprocessed_pseudogene;havana_transcript=OTTHUMT00000157411;havana_version=1;tag=basic;transcript_id=ENST00000411878;version=1\n";

    #[test]
    fn test_normalize_chr() {
        assert_eq!(normalize_chr("chr1"), "1");
        assert_eq!(normalize_chr("chrX"), "X");
        assert_eq!(normalize_chr("chrM"), "MT");
        assert_eq!(normalize_chr("1"), "1");
        assert_eq!(normalize_chr("MT"), "MT");
    }

    #[test]
    fn test_strip_prefix() {
        assert_eq!(strip_prefix("gene:ENSG00000123", "gene:"), "ENSG00000123");
        assert_eq!(
            strip_prefix("transcript:ENST00000456", "transcript:"),
            "ENST00000456"
        );
        assert_eq!(strip_prefix("ENSG00000123", "gene:"), "ENSG00000123");
    }

    #[test]
    fn test_parse_attributes() {
        let attrs = parse_attributes("ID=gene:ENSG00000123;Name=BRCA2;biotype=protein_coding");
        assert_eq!(attrs.get("ID").unwrap(), "gene:ENSG00000123");
        assert_eq!(attrs.get("Name").unwrap(), "BRCA2");
        assert_eq!(attrs.get("biotype").unwrap(), "protein_coding");
    }

    #[test]
    fn test_parse_attributes_url_decode() {
        let attrs = parse_attributes("description=Gene%20name%2C%20variant");
        assert_eq!(attrs.get("description").unwrap(), "Gene name, variant");
    }

    /// `ccdsid` is the first CCDS xref's display label (`Transcript.pm`
    /// `summary_as_hash` through `ccds`), the same entry Ensembl VEP stores
    /// as `_ccds` (`Database/Transcript.pm` `prefetch_transcript_ids`), so
    /// `CCDS74772.1` on ENST00000622113 is the cache's `ccds`.
    #[test]
    fn ccdsid_becomes_ccds() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}"));
        let tr = transcript(&genes, "ENST00000622113");
        assert_eq!(tr.ccds.as_deref(), Some("CCDS74772.1"));
        let no_ccds = TPTE_205.replace("ccdsid=CCDS74772.1;", "");
        let genes = parse(&format!("{TPTE_GENE}{no_ccds}"));
        assert_eq!(transcript(&genes, "ENST00000622113").ccds, None);
    }

    /// The gene line's `version` (`Feature.pm` `summary_as_hash`) is the
    /// `_gene_version` Ensembl VEP appends to `Gene` under `--gene_version`:
    /// `6` for ENSG00000274391.
    #[test]
    fn gene_version_from_gff3() {
        let genes = parse(TPTE_GENE);
        let gene = &genes["ENSG00000274391"];
        assert_eq!(gene.version, Some(6));
        assert_eq!(gene.name.as_deref(), Some("TPTE"));
        assert_eq!(gene.biotype, "protein_coding");
        let unversioned = TPTE_GENE.replace(";version=6", "");
        assert_eq!(parse(&unversioned)["ENSG00000274391"].version, None);
    }

    /// A `start` or `end` column that is not a number fails the whole parse
    /// with an error naming the column and the value, so a malformed line is
    /// never binned at position 0.
    #[test]
    fn non_numeric_coordinate_fails_the_parse() {
        for (lines, message) in [
            (
                TPTE_GENE.replace("\t10521553\t", "\tabc\t"),
                "GFF3 start \"abc\" is not a number",
            ),
            (
                TPTE_GENE.replace("\t10606140\t", "\t1060614x\t"),
                "GFF3 end \"1060614x\" is not a number",
            ),
        ] {
            let err = try_parse(&lines).unwrap_err();
            assert!(format!("{err:#}").contains(message), "{err:#}");
        }
    }

    /// A start of 0 or an end before its start fails the parse, naming the line:
    /// the region of a start is binned from `start - 1`, which a 0 has no room
    /// for, and a span whose end precedes its start has no length.
    #[test]
    fn a_start_of_zero_or_an_end_before_its_start_fails_the_parse() {
        for (lines, message) in [
            (
                TPTE_GENE.replace("\t10521553\t", "\t0\t"),
                "GFF3 coordinates 0-10606140 at line 1 are not 1-based and ordered",
            ),
            (
                TPTE_GENE.replace("\t10606140\t", "\t10521552\t"),
                "GFF3 coordinates 10521553-10521552 at line 1 are not 1-based and ordered",
            ),
        ] {
            let err = try_parse(&lines).unwrap_err();
            assert_eq!(format!("{err:#}"), message);
        }
        let genes = parse(&TPTE_GENE.replace("\t10606140\t", "\t10521553\t"));
        assert_eq!(genes["ENSG00000274391"].end, 10521553);
    }

    /// A seqid that would place the chromosome's directory outside the cache (a
    /// parent-directory walk, a rooted path, a `chr` that normalises to nothing)
    /// fails the parse, naming the seqid as written and its line, counted from
    /// the file's first line with the header included.
    #[test]
    fn seqid_that_is_not_one_path_component_fails_the_parse() {
        for (seqid, message) in [
            (
                "../../escaped",
                "invalid chromosome name '../../escaped' at line 2",
            ),
            (
                "/tmp/escaped",
                "invalid chromosome name '/tmp/escaped' at line 2",
            ),
            ("chr", "invalid chromosome name 'chr' at line 2"),
            ("..", "invalid chromosome name '..' at line 2"),
        ] {
            let lines = format!("##gff-version 3\n{}", TPTE_GENE.replacen("21", seqid, 1));
            let err = try_parse(&lines).unwrap_err();
            assert_eq!(format!("{err:#}"), message);
        }
        let genes = parse(&format!(
            "##gff-version 3\n{}",
            TPTE_GENE.replacen("21", "HSCHR6_MHC_COX", 1)
        ));
        assert_eq!(genes["ENSG00000274391"].chr, "HSCHR6_MHC_COX");
    }

    /// Exons with a rank sort by it; once any exon of the transcript lacks one,
    /// every exon sorts by its start, as a rank compared with a start orders
    /// nothing, and the file is parsed whole.
    #[test]
    fn exons_mixing_ranked_and_unranked_sort_by_start() {
        let exon = |start: u64, rank: Option<u32>| {
            let rank = rank.map(|r| format!(";rank={r}")).unwrap_or_default();
            format!(
                "21\tensembl_havana\texon\t{start}\t{}\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE{start}{rank}\n",
                start + 99
            )
        };
        let exon_starts = |genes: &HashMap<String, GeneRecord>| -> Vec<u64> {
            transcript(genes, "ENST00000622113")
                .exons
                .iter()
                .map(|e| e.start)
                .collect()
        };
        // Sixty-four exons, every other one ranked, listed in a scrambled order.
        let mut lines = format!("{TPTE_GENE}{TPTE_205}");
        for i in (0..64u64).map(|i| (i * 37) % 64) {
            let rank = (i % 2 == 0).then_some(64 - i as u32);
            lines.push_str(&exon(10_521_600 + i * 1_000, rank));
        }
        let expected: Vec<u64> = (0..64).map(|i| 10_521_600 + i * 1_000).collect();
        assert_eq!(exon_starts(&parse(&lines)), expected);

        let ranked = format!(
            "{TPTE_GENE}{TPTE_205}{}{}",
            exon(10_540_000, Some(2)),
            exon(10_600_000, Some(1))
        );
        assert_eq!(exon_starts(&parse(&ranked)), [10_600_000, 10_540_000]);
    }

    /// A line past the cap (the gene line with an attribute value of the cap's
    /// length) stops the parse naming the file and its line, 2.
    #[test]
    fn line_past_the_cap_fails_the_parse() {
        let padding: String = std::iter::repeat_n('x', MAX_LINE_BYTES).collect();
        let lines = format!(
            "{}{}",
            TPTE_205,
            TPTE_GENE.replace("Name=TPTE;", &format!("Name=TPTE;Note={padding};"))
        );
        let err = try_parse(&lines).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "t.gff3 line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
    }

    /// The GENCODE basic tag is `gencode_basic` on the GRCh38 file (the
    /// attribute code, `Transcript.pm` `summary_as_hash`) and `basic` on the
    /// GRCh37 release-87 file; either marks the transcript.
    #[test]
    fn gencode_basic_tag_both_spellings() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}"));
        assert!(transcript(&genes, "ENST00000622113").gencode_basic);

        let genes = parse(GRCH37_PSEUDOGENE);
        let tr = transcript(&genes, "ENST00000411878");
        assert!(tr.gencode_basic);
        assert!(!tr.gencode_primary);

        let untagged = TPTE_205.replace("tag=gencode_basic,gencode_primary;", "");
        let genes = parse(&format!("{TPTE_GENE}{untagged}"));
        assert!(!transcript(&genes, "ENST00000622113").gencode_basic);
    }

    /// `gencode_primary` in `tag` is the GENCODE primary attribute's code
    /// (`Transcript.pm` `summary_as_hash`).
    #[test]
    fn gencode_primary_tag_from_gff3() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}"));
        assert!(transcript(&genes, "ENST00000622113").gencode_primary);
        let basic_only =
            TPTE_205.replace("tag=gencode_basic,gencode_primary;", "tag=gencode_basic;");
        let genes = parse(&format!("{TPTE_GENE}{basic_only}"));
        let tr = transcript(&genes, "ENST00000622113");
        assert!(tr.gencode_basic);
        assert!(!tr.gencode_primary);
    }

    /// `transcript_support_level` is the `TSL` attribute value without its
    /// `tsl` prefix (`Transcript.pm` `tsl`): `1` here, and the `assigned to
    /// previous version` suffix is kept verbatim, as it is in the attribute
    /// Ensembl VEP prints the digit from (`OutputFactory.pm` `m/tsl(\d+)/`).
    #[test]
    fn tsl_from_transcript_support_level() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}"));
        assert_eq!(
            transcript(&genes, "ENST00000622113").tsl.as_deref(),
            Some("1")
        );

        let col6a2 = "21\tensembl_havana\tgene\t46096722\t46132851\t.\t+\t.\tID=gene:ENSG00000142173;Name=COL6A2;biotype=protein_coding;description=collagen type VI alpha 2 chain [Source:HGNC Symbol%3BAcc:HGNC:2212];gene_id=ENSG00000142173;logic_name=ensembl_havana_gene_homo_sapiens;version=19\n\
21\thavana\tmRNA\t46098112\t46129801\t.\t+\t.\tID=transcript:ENST00000397763;Parent=gene:ENSG00000142173;Name=COL6A2-202;biotype=protein_coding;ccdsid=CCDS13729.1;tag=gencode_basic,gencode_primary,MANE_Plus_Clinical;transcript_id=ENST00000397763;transcript_support_level=5 (assigned to previous version 5);version=6\n";
        let genes = parse(col6a2);
        let tr = transcript(&genes, "ENST00000397763");
        assert_eq!(
            tr.tsl.as_deref(),
            Some("5 (assigned to previous version 5)")
        );
        assert!(tr.mane_plus_clinical);
        assert!(!tr.mane_select);
        assert!(!tr.is_canonical);

        let no_tsl = TPTE_205.replace("transcript_support_level=1;", "");
        let genes = parse(&format!("{TPTE_GENE}{no_tsl}"));
        assert_eq!(transcript(&genes, "ENST00000622113").tsl, None);
    }

    /// A CDS line carries `protein_id`, the translation's stable identifier,
    /// and `version`, the translation's version (`CDS.pm` `summary_as_hash`):
    /// ENSP00000482040 version 1 for ENST00000622113, the identifier Ensembl
    /// VEP stores as `_protein` (`Database/Transcript.pm`
    /// `prefetch_translation_ids`). The GRCh37 release-87 file has no
    /// `version` on its CDS lines.
    #[test]
    fn protein_id_and_version_from_cds_line() {
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{TPTE_205_CDS}"));
        let tr = transcript(&genes, "ENST00000622113");
        assert_eq!(tr.cds_regions.len(), 2);
        for cds in &tr.cds_regions {
            assert_eq!(cds.protein_id.as_deref(), Some("ENSP00000482040"));
            assert_eq!(cds.protein_version, Some(1));
        }
        assert_eq!(tr.cds_regions[0].phase, 0);
        assert_eq!(tr.cds_regions[1].phase, 1);

        let grch37_cds = TPTE_205_CDS.replace(";version=1", "");
        let genes = parse(&format!("{TPTE_GENE}{TPTE_205}{grch37_cds}"));
        let cds = &transcript(&genes, "ENST00000622113").cds_regions[0];
        assert_eq!(cds.protein_id.as_deref(), Some("ENSP00000482040"));
        assert_eq!(cds.protein_version, None);
    }

    /// Ensembl VEP's dumper skips a transcript of biotype `artifact` and a
    /// readthrough transcript (`Database/Transcript.pm`
    /// `get_features_by_regions_uncached`, `next if $tr->{biotype} eq
    /// 'artifact'` and `next if @{$tr->get_all_Attributes("readthrough_tra")}`).
    /// The GRCh38 release-116 file spells the readthrough tag
    /// `readthrough_transcript`; the core API's attribute code is
    /// `readthrough_tra`, and both are dropped.
    #[test]
    fn readthrough_and_artifact_transcripts_dropped() {
        let lines = "21\tensembl_havana\tgene\t5022531\t5036771\t.\t+\t.\tID=gene:ENSG00000277117;biotype=artifact;description=novel protein%2C similar to inducible T-cell co-stimulator ligand ICOSLG;gene_id=ENSG00000277117;logic_name=ensembl_havana_gene_homo_sapiens;version=6\n\
21\thavana\ttranscript\t5022531\t5036771\t.\t+\t.\tID=transcript:ENST00000623960;Parent=gene:ENSG00000277117;biotype=artifact;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000623960;transcript_support_level=1 (assigned to previous version 3);version=4\n\
21\thavana\texon\t5022531\t5022693\t.\t+\t.\tParent=transcript:ENST00000623960;Name=ENSE00003755585;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003755585;rank=1;version=1\n\
21\thavana\tncRNA_gene\t10482738\t10605716\t.\t+\t.\tID=gene:ENSG00000273840;biotype=lncRNA;description=novel transcript;gene_id=ENSG00000273840;logic_name=havana_homo_sapiens;version=2\n\
21\thavana\tlnc_RNA\t10482738\t10605716\t.\t+\t.\tID=transcript:ENST00000612267;Parent=gene:ENSG00000273840;biotype=lncRNA;tag=gencode_basic,gencode_primary,Ensembl_canonical,readthrough_transcript;transcript_id=ENST00000612267;transcript_support_level=5;version=1\n\
21\thavana\tlnc_RNA\t10482738\t10500000\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000273840;biotype=lncRNA;tag=gencode_basic,readthrough_tra;transcript_id=ENST00000000001;version=1\n\
21\thavana\tlnc_RNA\t10482738\t10490000\t.\t+\t.\tID=transcript:ENST00000000002;Parent=gene:ENSG00000273840;biotype=lncRNA;tag=gencode_basic;transcript_id=ENST00000000002;version=1\n";
        let genes = parse(lines);
        assert!(genes["ENSG00000277117"].transcripts.is_empty());
        let kept: Vec<&str> = genes["ENSG00000273840"]
            .transcripts
            .iter()
            .map(|t| t.transcript_id.as_str())
            .collect();
        assert_eq!(kept, vec!["ENST00000000002"]);
    }

    /// A gene is any line whose `ID` is `gene:<id>`, so a lncRNA gene's
    /// `ncRNA_gene` line and a `pseudogene` gene line are genes; a
    /// transcript is any `transcript:<id>` line, so a `pseudogenic_transcript`
    /// line is a transcript, as is the GRCh37 release-87 file's `pseudogene`
    /// transcript line.
    #[test]
    fn gene_and_transcript_lines_by_id_prefix() {
        let lines = "21\thavana\tncRNA_gene\t5232668\t5243903\t.\t-\t.\tID=gene:ENSG00000279669;biotype=lncRNA;description=novel transcript;gene_id=ENSG00000279669;logic_name=havana_homo_sapiens;version=2\n\
21\thavana_tagene\tlnc_RNA\t5232668\t5243903\t.\t-\t.\tID=transcript:ENST00000623753;Parent=gene:ENSG00000279669;biotype=lncRNA;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000623753;transcript_support_level=3 (assigned to previous version 1);version=2\n\
21\tensembl_havana\tpseudogene\t9068414\t9073609\t.\t-\t.\tID=gene:ENSG00000188681;Name=TEKT4P2;biotype=transcribed_unprocessed_pseudogene;description=tektin 4 pseudogene 2 [Source:HGNC Symbol%3BAcc:HGNC:40046];gene_id=ENSG00000188681;logic_name=ensembl_havana_gene_homo_sapiens;version=14\n\
21\thavana\tpseudogenic_transcript\t9068414\t9073609\t.\t-\t.\tID=transcript:ENST00000559466;Parent=gene:ENSG00000188681;Name=TEKT4P2-201;biotype=transcribed_unprocessed_pseudogene;tag=gencode_basic,gencode_primary,Ensembl_canonical;transcript_id=ENST00000559466;transcript_support_level=NA;version=1\n\
21\t.\tbiological_region\t5000000\t5000100\t.\t+\t.\tlogic_name=cpg\n";
        let genes = parse(lines);
        assert_eq!(genes.len(), 2);
        let lnc = transcript(&genes, "ENST00000623753");
        assert_eq!(lnc.biotype, "lncRNA");
        assert!(lnc.is_canonical);
        assert_eq!(lnc.source, "havana_tagene");
        assert_eq!(
            lnc.tsl.as_deref(),
            Some("3 (assigned to previous version 1)")
        );
        assert_eq!(genes["ENSG00000279669"].version, Some(2));
        assert_eq!(genes["ENSG00000279669"].name, None);
        let pseudo = transcript(&genes, "ENST00000559466");
        assert_eq!(pseudo.gene_id, "ENSG00000188681");
        assert_eq!(pseudo.tsl.as_deref(), Some("NA"));
        assert_eq!(genes["ENSG00000188681"].name.as_deref(), Some("TEKT4P2"));

        let genes = parse(GRCH37_PSEUDOGENE);
        assert_eq!(genes.len(), 1);
        assert_eq!(
            transcript(&genes, "ENST00000411878").biotype,
            "unprocessed_pseudogene"
        );
    }

    /// An exon line's `ensembl_phase` and `ensembl_end_phase` are the exon's
    /// own phases (`Exon.pm` `summary_as_hash`), the values the Ensembl
    /// transcript cache stores on each exon: -1 on a non-coding exon, and on
    /// ENST00000622113 exon 4 starts in the 5' UTR (-1) and ends in frame 2.
    #[test]
    fn ensembl_phases_read_from_exon_line() {
        let lines = format!(
            "{TPTE_GENE}{TPTE_205}\
21\tensembl_havana\texon\t10521553\t10521694\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00003729301;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003729301;rank=1;version=1\n\
21\tensembl_havana\texon\t10538681\t10538734\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00003735534;constitutive=0;ensembl_end_phase=2;ensembl_phase=-1;exon_id=ENSE00003735534;rank=4;version=1\n\
21\tensembl_havana\texon\t10541112\t10541165\t.\t+\t.\tParent=transcript:ENST00000622113;Name=ENSE00002295075;constitutive=0;ensembl_end_phase=2;ensembl_phase=2;exon_id=ENSE00002295075;rank=5;version=2\n"
        );
        let genes = parse(&lines);
        let exons = &transcript(&genes, "ENST00000622113").exons;
        let phases: Vec<(Option<i8>, Option<i8>)> = exons
            .iter()
            .map(|e| (e.ensembl_phase, e.ensembl_end_phase))
            .collect();
        assert_eq!(
            phases,
            vec![
                (Some(-1), Some(-1)),
                (Some(-1), Some(2)),
                (Some(2), Some(2))
            ]
        );
        assert_eq!(exons[1].exon_id.as_deref(), Some("ENSE00003735534"));
        assert_eq!(exons[1].rank, Some(4));
    }
}
