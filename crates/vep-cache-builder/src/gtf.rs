// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript tags from an Ensembl GTF.
//!
//! Ensembl's GFF3 carries a transcript's `basic`, `Ensembl_canonical`, MANE and
//! support-level tags but not the CDS completeness flags; the GTF of the same
//! release repeats `tag "..."` per transcript line and includes `cds_start_NF`
//! and `cds_end_NF`. VEP's own cache keeps those two plus `gencode_basic` and
//! `gencode_primary` as transcript attributes ([`KEPT_TAGS`]), and its output
//! prints the `cds_` ones in FLAGS. The GENCODE basic tag is spelled `basic` in
//! the GRCh37 (release 87 annotation) GTF and `gencode_basic` in the GRCh38 one.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

/// GTF tags VEP's cache dumper stores as transcript attributes, with the
/// attribute code each becomes.
pub const KEPT_TAGS: [(&str, &str); 5] = [
    ("basic", "gencode_basic"),
    ("gencode_basic", "gencode_basic"),
    ("gencode_primary", "gencode_primary"),
    ("cds_start_NF", "cds_start_NF"),
    ("cds_end_NF", "cds_end_NF"),
];

/// Header fields (`#!key value`) and per-transcript tags of a GTF.
#[derive(Debug, Default)]
pub struct GtfTags {
    /// `#!genome-build`, `#!genebuild-last-updated` and the other header keys.
    pub header: HashMap<String, String>,
    /// Unversioned transcript id to every `tag` value on its `transcript` line.
    pub tags: HashMap<String, Vec<String>>,
}

/// Reads the `transcript` lines of a (possibly gzipped) Ensembl GTF.
pub fn parse_gtf_tags(path: &Path) -> Result<GtfTags> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut out = GtfTags::default();
    for line in reader.lines() {
        let line = line.context("Failed to read GTF line")?;
        if let Some(rest) = line.strip_prefix("#!") {
            if let Some((k, v)) = rest.split_once(' ') {
                out.header.insert(k.to_string(), v.trim().to_string());
            }
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let feature = fields.nth(2);
        if feature != Some("transcript") {
            continue;
        }
        let Some(attrs) = fields.nth(5) else {
            continue;
        };
        let mut transcript_id: Option<String> = None;
        let mut tags = Vec::new();
        for attr in attrs.split(';') {
            let attr = attr.trim();
            let Some((key, value)) = attr.split_once(' ') else {
                continue;
            };
            let value = value.trim().trim_matches('"');
            match key {
                "transcript_id" => {
                    transcript_id = Some(value.split('.').next().unwrap_or(value).to_string());
                }
                "tag" => tags.push(value.to_string()),
                _ => {}
            }
        }
        if let Some(id) = transcript_id {
            out.tags.entry(id).or_default().extend(tags);
        }
    }
    Ok(out)
}

impl GtfTags {
    /// The attribute codes VEP keeps, for one transcript, in [`KEPT_TAGS`] order,
    /// each code once.
    pub fn kept_attribute_codes(&self, transcript_id: &str) -> Vec<&'static str> {
        let unversioned = transcript_id.split('.').next().unwrap_or(transcript_id);
        let Some(tags) = self.tags.get(unversioned) else {
            return Vec::new();
        };
        let mut codes: Vec<&'static str> = Vec::new();
        for (tag, code) in KEPT_TAGS.iter() {
            if tags.iter().any(|t| t == tag) && !codes.contains(code) {
                codes.push(code);
            }
        }
        codes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transcript_tags_and_header_are_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf");
        std::fs::write(
            &path,
            "#!genome-build GRCh38.p14\n#!genebuild-last-updated 2024-01\n\
             1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825\"; gene_version \"2\";\n\
             1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825\"; transcript_id \"ENST00000456328\"; transcript_version \"2\"; tag \"basic\"; tag \"Ensembl_canonical\"; transcript_support_level \"1\";\n\
             1\thavana\ttranscript\t65419\t71585\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; tag \"cds_start_NF\"; tag \"mRNA_start_NF\"; tag \"cds_end_NF\";\n\
             1\thavana\ttranscript\t90000\t91000\t.\t+\t.\tgene_id \"ENSG00000000003\"; transcript_id \"ENST00000000003\"; tag \"gencode_basic\"; tag \"gencode_primary\"; tag \"Ensembl_canonical\";\n\
             1\thavana\texon\t65419\t65433\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; tag \"cds_start_NF\";\n",
        )
        .unwrap();
        let gtf = parse_gtf_tags(&path).unwrap();
        assert_eq!(gtf.header["genome-build"], "GRCh38.p14");
        assert_eq!(gtf.header["genebuild-last-updated"], "2024-01");
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000456328.2"),
            vec!["gencode_basic"]
        );
        // mRNA_start_NF is not an attribute VEP keeps; exon lines do not count.
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000641515"),
            vec!["cds_start_NF", "cds_end_NF"]
        );
        // The GRCh38 GTF spells the basic tag `gencode_basic`.
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000000003"),
            vec!["gencode_basic", "gencode_primary"]
        );
        assert!(gtf.kept_attribute_codes("ENST00000000000").is_empty());
    }
}
