// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A GTF as a transcript model, and the transcript tags of an Ensembl GTF.
//!
//! [`parse_gtf_model`] reads the `gene`, `transcript`, `exon`, `CDS` and
//! `stop_codon` lines into the gene tree of [`crate::gff3`], as Ensembl VEP
//! reads a GTF given to it as an annotation source
//! (`AnnotationSource/File/GTF.pm`): an exon, CDS or stop codon belongs to
//! each transcript its `transcript_id` names, a transcript to the gene its
//! `gene_id` names, a stop codon is one more CDS segment, a transcript's
//! biotype is its `transcript_biotype`, `transcript_type` or `biotype`, else
//! its source column, and its CCDS identifier is its `ccds_id`, as the line
//! carries it.
//!
//! [`parse_gtf_tags`] reads an Ensembl GTF beside an Ensembl GFF3. The GFF3
//! carries a transcript's `basic`, `Ensembl_canonical`, MANE and
//! support-level tags but not the CDS completeness flags; the GTF of the same
//! release repeats `tag "..."` per transcript line and includes `cds_start_NF`
//! and `cds_end_NF`. Ensembl VEP's own cache keeps those two plus
//! `gencode_basic` and `gencode_primary` as transcript attributes
//! ([`KEPT_TAGS`]), and its output prints the `cds_` ones in FLAGS. The
//! GENCODE basic tag is spelled `basic` in the GRCh37 (release 87 annotation)
//! and GENCODE GTFs and `gencode_basic` in the GRCh38 Ensembl one; both
//! readings take a tag through the one table, so `basic` is `gencode_basic`
//! on either.

use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use vep_core::coordinate::is_cache_path_component;
use vep_core::regulatory_gff::{undecodable_attribute_error, uri_unescape};
use vep_core::sv_routing::perl_split;

use crate::core_dumps::BoundedLines;
use crate::gff3::{
    attribute, check_skipped_coordinates, feature_version_keys, is_set, located, names_a_parent,
    normalize_chr, opens_with_a_blank_column, parent_ids, parse_version, without_byte_order_marks,
    CdsRecord, ExonRecord, Feature, GeneRecord, Model, Reading, Records,
};

/// The longest GTF line read, 16 MiB. A transcript line's widest column is its
/// attributes, a few kilobytes even with every tag listed, so the cap refuses
/// a file with no line end in it long before the line buffer grows to its
/// size.
const MAX_LINE_BYTES: usize = 16 << 20;

/// GTF tags Ensembl VEP's cache dumper stores as transcript attributes, with
/// the attribute code each becomes. The one table of both GTF readings: a tag
/// list is spelled by [`attribute_code`] as it is read, so `basic` is
/// `gencode_basic` to the model and to the tags beside a GFF3 alike.
pub const KEPT_TAGS: [(&str, &str); 5] = [
    ("basic", "gencode_basic"),
    ("gencode_basic", "gencode_basic"),
    ("gencode_primary", "gencode_primary"),
    ("cds_start_NF", "cds_start_NF"),
    ("cds_end_NF", "cds_end_NF"),
];

/// The attribute code a GTF tag stands for ([`KEPT_TAGS`]): `gencode_basic`
/// for `basic`; a tag the table does not name is itself.
fn attribute_code(tag: &str) -> &str {
    KEPT_TAGS
        .iter()
        .find(|(kept, _)| *kept == tag)
        .map_or(tag, |(_, code)| code)
}

/// The GTF feature types Ensembl VEP reads from a file given as an annotation
/// source (`GTF.pm` `%INCLUDE_FEATURE_TYPES`); a line of any other type is
/// skipped.
pub const GTF_FEATURE_TYPES: [&str; 6] = ["cds", "CDS", "stop_codon", "exon", "gene", "transcript"];

/// Header fields (`#!key value`) and per-transcript tags of a GTF.
#[derive(Debug, Default)]
pub struct GtfTags {
    /// `#!genome-build`, `#!genebuild-last-updated` and the other header keys.
    pub header: HashMap<String, String>,
    /// A `transcript` line's name ([`tagged_name`]: its `transcript_id` whole,
    /// with its version behind a dot when it carries one) to every `tag`
    /// value on the line, each spelled as its attribute code
    /// ([`attribute_code`]).
    pub tags: HashMap<String, Vec<String>>,
}

/// A GTF read as a transcript model.
#[derive(Debug)]
pub struct GtfModel {
    /// `#!genome-build`, `#!genebuild-last-updated` and the other header keys.
    pub header: HashMap<String, String>,
    pub genes: HashMap<String, GeneRecord>,
}

/// Opens a plain, gzipped or bgzipped GTF for line reading.
fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    Ok(if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// The GTF attribute keys the readers take as one value. A line repeating one
/// is refused: ensembl-io lists a repeated key's values (`GTF.pm` lines
/// 128-133), and Ensembl VEP, reading the list where it expects one value,
/// prints a symbol, biotype or protein identifier as the values joined by
/// commas, a version (and a protein identifier under its version flag) as the
/// list's address, and loses the record whose identifier or parent the list
/// was. `tag` is read as the list it is, and a key no reader takes keeps its
/// last value.
const SINGLE_VALUE_KEYS: [&str; 20] = [
    "gene_id",
    "transcript_id",
    "exon_id",
    "protein_id",
    "version",
    "protein_version",
    "transcript_version",
    "transcript_biotype",
    "transcript_type",
    "biotype",
    "gene_name",
    "gene",
    "name",
    "Name",
    "transcript_support_level",
    "tsl",
    "ccdsid",
    "ccds_id",
    "dbxref",
    "Dbxref",
];

/// Perl's `\s` on a byte string: the six ASCII whitespace characters, the
/// vertical tab among them, and nothing above ASCII.
fn is_perl_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')
}

/// A GTF attribute column as ensembl-io's `Parser/GTF.pm` `get_attributes`
/// reads it (lines 117-135): the column split on `;`; each piece split at its
/// first run of whitespace into key and value, the leading whitespace skipped
/// and the rest of the piece the value, trailing whitespace included (Perl's
/// `split(' ', $attr, 2)`, line 120); a piece with no key, or whose key is
/// `0`, skipped (line 121); every `"` removed from the value wherever it
/// stands (line 123); the value URI-unescaped (line 126,
/// [`vep_core::regulatory_gff::uri_unescape`]), so `"A%2CB"` is `A,B`,
/// `"a"b"` is `ab` and `"%30"` the unset `0`. A piece of one word, a key
/// without a value, sets nothing, as the undefined value it gives Ensembl VEP
/// reads, but repeats a key already set (lines 128-131 list the undefined
/// value after the set one). The `tag` values are the list Ensembl VEP reads
/// (`BaseGXF.pm` lines 683-685: the one value of a `tag` set once split on
/// `,` when it is true to Perl, so `"0"` and `""` list nothing; the values
/// of a repeated `tag` whole), each spelled as the attribute
/// code it stands for ([`attribute_code`]), and built when `with_tags`:
/// without it the list is left unbuilt and empty, the attributes and the
/// errors the same either way. A value whose escapes spell bytes that are
/// not UTF-8, and a line repeating a key of [`SINGLE_VALUE_KEYS`], are
/// errors naming the key.
fn parse_gtf_attributes(
    column: &str,
    with_tags: bool,
) -> Result<(HashMap<String, String>, Vec<String>)> {
    let mut attrs = HashMap::new();
    let mut tag_values = Vec::new();
    let mut tag_repeated = false;
    for piece in column.split(';') {
        let piece = piece.trim_start_matches(is_perl_space);
        let key_end = piece.find(is_perl_space).unwrap_or(piece.len());
        let (key, rest) = piece.split_at(key_end);
        if !is_set(key) {
            continue;
        }
        let repeated = attrs.contains_key(key);
        if repeated && SINGLE_VALUE_KEYS.contains(&key) {
            bail!("GTF attribute {key} is repeated; Ensembl VEP reads a repeated key as a list");
        }
        if repeated && key == "tag" {
            tag_repeated = true;
        }
        if rest.is_empty() {
            continue;
        }
        let unquoted = rest.trim_start_matches(is_perl_space).replace('"', "");
        let value = match uri_unescape(&unquoted) {
            Cow::Borrowed(_) => unquoted,
            Cow::Owned(bytes) => String::from_utf8(bytes)
                .map_err(|_| anyhow!(undecodable_attribute_error("GTF", key, &unquoted)))?,
        };
        if with_tags && key == "tag" {
            tag_values.push(value.clone());
        }
        attrs.insert(key.to_string(), value);
    }
    let tags = match tag_values.as_slice() {
        [one] if !tag_repeated && !is_set(one) => Vec::new(),
        [one] if !tag_repeated => perl_split(one, &[','])
            .map(|tag| attribute_code(tag).to_string())
            .collect(),
        _ => tag_values
            .iter()
            .map(|tag| attribute_code(tag).to_string())
            .collect(),
    };
    Ok((attrs, tags))
}

/// Reads a GTF (plain, gzipped or bgzipped) as a transcript model, with the
/// length of a line capped: a line longer than `MAX_LINE_BYTES` fails the parse, naming
/// the file and the line, as does a line that is not UTF-8. Fails as
/// [`crate::gff3::parse_gff3`] fails on a coordinate that is not a number or
/// not ordered, on the start or end no number leads, or the end below zero,
/// of a line it skips, for its type or for having fewer than nine columns
/// ([`check_skipped_coordinates`]), on a line that is no comment and ends
/// before its start column, on a line that opens with a tab, read at the
/// columns the file writes, or with white space and a tab
/// ([`opens_with_a_blank_column`]), on a seqid that cannot name a
/// chromosome's cache directory,
/// and on a version that is not a whole number or is above the largest the
/// cache holds, read from a gene or transcript line or, under the keys of
/// `reading`, a CDS line, and on an attribute column
/// [`parse_gtf_attributes`] refuses; each error opens with the file and the
/// line. A byte-order mark at the start of a line is taken off before the
/// line is read ([`without_byte_order_marks`]), so a `#!` key behind one is
/// read.
///
/// A record is named by its `<type>_id` attribute (`gene_id`,
/// `transcript_id`, `exon_id`), else by the MD5 of its line (`GTF.pm`
/// `_record_get_id`); an exon, CDS or stop codon's parents are its
/// `transcript_id` list and a transcript's parent its `gene_id` (`GTF.pm`
/// `%PARENTS`). A stop codon is one more CDS segment of its transcript
/// (`BaseGXF.pm` `_create_transcript`, `cds` and `stop_codon` alike), so the
/// coding sequence runs through the stop codon as an Ensembl GFF3's does.
/// `reading` decides a translation's version, and whether a parentless record
/// whose identifier opens with `gene:` is a gene ([`Reading`] names the rule),
/// and nothing else. For the version,
/// Ensembl VEP's `--gtf` source reads the CDS line's `version`, an attribute
/// no Ensembl GTF carries, else 1 (`BaseGXF.pm` `_add_translation`); a build
/// of the cache's content reads `version`, else `protein_version`, the
/// attribute the Ensembl GTF writes the translation's version as, else
/// leaves it to the core dumps or the protein FASTA header. A GTF's lines
/// are read by [`GTF_FEATURE_TYPES`] under either reading.
pub fn parse_gtf_model(path: &Path, reading: Reading) -> Result<GtfModel> {
    let mut header = HashMap::new();
    let mut records = Records::new(reading);
    let version_keys: &[&str] = match reading {
        Reading::GffSource => &["version"],
        Reading::EnsemblCache => &["version", "protein_version"],
    };
    for (line_index, line) in
        BoundedLines::new(open(path)?, MAX_LINE_BYTES, path.display().to_string()).enumerate()
    {
        let line_number = line_index + 1;
        let written = String::from_utf8(without_byte_order_marks(line?))
            .map_err(|_| anyhow!("GTF line is not UTF-8"))
            .with_context(|| located(path, line_number))?;
        if let Some(rest) = written.strip_prefix("#!") {
            if let Some((k, v)) = rest.split_once(' ') {
                header.insert(k.to_string(), v.trim().to_string());
            }
            continue;
        }
        let line = written.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // A record without its start, as the GFF3 reader refuses one: a GTF
        // has no sequence section, so every such line is a record's.
        let columns = written.split('\t').count();
        if columns < 4 {
            return Err(anyhow!(
                "GTF line ends at column {columns}, before its start"
            ))
            .with_context(|| located(path, line_number));
        }
        // A line that opens with a tab is read at the columns the file
        // writes, as the GFF3 reader reads one.
        let line = if written.starts_with('\t') {
            written.trim_end()
        } else {
            line
        };
        let fields: Vec<&str> = line.split('\t').collect();
        let coordinate = |column: &str, value: &str| -> Result<u64> {
            value
                .parse()
                .map_err(|_| anyhow!("GTF {column} {value:?} is not a number"))
                .with_context(|| located(path, line_number))
        };
        if fields.len() < 9 {
            check_skipped_coordinates(&written, coordinate)?;
            continue;
        }
        // The include test comes before any other column is read, as it does
        // in the GFF3 reader: Ensembl VEP's source reads a record's raw start,
        // for its loop's bound, and its type, and nothing else of one it skips
        // (`BaseGXF.pm:217-245`). A skipped line's start and end are held to
        // the one rule the GFF3 reader holds them to.
        let feature_type = fields[2];
        if !GTF_FEATURE_TYPES.contains(&feature_type) {
            check_skipped_coordinates(&written, coordinate)?;
            records.skip(feature_type, "");
            continue;
        }
        // White space and a tab before a record put its type where the index
        // build reads a start, as in the GFF3 reader: the line is held to
        // that reading, as its twin of a skipped type is.
        if opens_with_a_blank_column(&written) {
            check_skipped_coordinates(&written, coordinate)?;
        }
        let chr = normalize_chr(fields[0]);
        if !is_cache_path_component(chr) {
            return Err(anyhow!("invalid chromosome name '{}'", fields[0]))
                .with_context(|| located(path, line_number));
        }
        let start = coordinate("start", fields[3])?;
        let end = coordinate("end", fields[4])?;
        if start == 0 || end < start {
            return Err(anyhow!(
                "GTF coordinates {start}-{end} are not 1-based and ordered"
            ))
            .with_context(|| located(path, line_number));
        }
        let strand: i8 = match fields[6] {
            "-" => -1,
            _ => 1,
        };
        let phase: i8 = fields[7].parse().unwrap_or(-1);
        // The tag list is built in the arm whose record keeps it
        // ([`Feature::tags`]): an exon, CDS or stop codon line repeats its
        // transcript's tags, and its record has no field for them.
        let attributes = |with_tags: bool| {
            parse_gtf_attributes(fields[8], with_tags).with_context(|| located(path, line_number))
        };

        match feature_type {
            "exon" => {
                let (attrs, _) = attributes(false)?;
                records.add_exon(
                    &parent_ids(Model::Gtf, &attrs, &["transcript_id"]),
                    ExonRecord {
                        exon_id: attrs.get("exon_id").cloned(),
                        start,
                        end,
                        strand,
                        rank: None,
                        ensembl_phase: None,
                        ensembl_end_phase: None,
                    },
                );
            }
            "cds" | "CDS" | "stop_codon" => {
                let (attrs, _) = attributes(false)?;
                records.add_cds(
                    &parent_ids(Model::Gtf, &attrs, &["transcript_id"]),
                    CdsRecord {
                        start,
                        end,
                        strand,
                        phase,
                        protein_id: attribute(&attrs, &["protein_id"]),
                        protein_version: parse_version(Model::Gtf, &attrs, version_keys)
                            .with_context(|| located(path, line_number))?,
                        line: line_number,
                    },
                );
            }
            _ => {
                let (attrs, tags) = attributes(true)?;
                let parents = match feature_type {
                    "transcript" => parent_ids(Model::Gtf, &attrs, &["gene_id"]),
                    _ => Vec::new(),
                };
                let id = attribute(&attrs, &[&format!("{feature_type}_id")])
                    .unwrap_or_else(|| format!("{:x}", md5::compute(line)));
                let child = feature_type == "transcript"
                    && names_a_parent(attrs.get("gene_id").map(String::as_str));
                let version = parse_version(
                    Model::Gtf,
                    &attrs,
                    feature_version_keys(reading, child, &id, feature_type),
                )
                .with_context(|| located(path, line_number))?;
                records.add_feature(Feature {
                    id,
                    parents,
                    child,
                    feature_type: feature_type.to_string(),
                    version,
                    chr: chr.to_string(),
                    start,
                    end,
                    strand,
                    source: fields[1].to_string(),
                    attrs,
                    tags,
                    ensembl: false,
                });
            }
        }
    }
    let genes = records.assemble(Model::Gtf, path)?;
    Ok(GtfModel { header, genes })
}

/// Reads the `transcript` lines of an Ensembl GTF (plain, gzipped or
/// bgzipped), with the length of a line capped: a line longer than
/// `MAX_LINE_BYTES` fails the parse, naming the file and the line, as does a
/// line that is not UTF-8 and a line [`parse_gtf_attributes`] refuses. A
/// byte-order mark at the start of a line is taken off first, as
/// [`parse_gtf_model`] takes it.
pub fn parse_gtf_tags(path: &Path) -> Result<GtfTags> {
    let mut out = GtfTags::default();
    for (line_index, line) in
        BoundedLines::new(open(path)?, MAX_LINE_BYTES, path.display().to_string()).enumerate()
    {
        let line = String::from_utf8(without_byte_order_marks(line?))
            .map_err(|_| anyhow!("GTF line is not UTF-8"))
            .with_context(|| located(path, line_index + 1))?;
        if let Some(rest) = line.strip_prefix("#!") {
            if let Some((k, v)) = rest.split_once(' ') {
                out.header.insert(k.to_string(), v.trim().to_string());
            }
            continue;
        }
        // The line's end is trimmed as `parse_gtf_model` trims it, so a value that
        // ends the line without a closing `;` reads the same to both readers.
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
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
        let (attrs, tags) =
            parse_gtf_attributes(attrs, true).with_context(|| located(path, line_index + 1))?;
        if let Some(id) = attrs.get("transcript_id") {
            out.tags
                .entry(tagged_name(id, &attrs))
                .or_default()
                .extend(tags);
        }
    }
    Ok(out)
}

/// The name a `transcript` line's tags are held under: its `transcript_id`
/// whole, every dot kept, and behind one more dot its version when the line
/// carries one (`version`, else `transcript_version`, as
/// [`feature_version_keys`] orders them), read as the number (`01` is 1). It
/// is the form a protein FASTA header writes (`ENST00000456328.2`); a
/// version that is not a whole number the cache holds names nothing.
fn tagged_name(id: &str, attrs: &HashMap<String, String>) -> String {
    match parse_version(Model::Gtf, attrs, &["version", "transcript_version"]) {
        Ok(Some(version)) => format!("{id}.{version}"),
        Ok(None) | Err(_) => id.to_string(),
    }
}

impl GtfTags {
    /// The attribute codes Ensembl VEP keeps, for the model's transcript `id`
    /// at `version`, in [`KEPT_TAGS`] order, each code once; the tag list
    /// holds codes ([`attribute_code`]), so a code is matched, never a file
    /// spelling. The tags are those of the line named `id.version` when the
    /// model states a version and the file holds that name, else of the line
    /// named `id`, each compared whole, as [`crate::fasta::record_of`] finds
    /// a protein FASTA record: the line of another version of `id`, or of a
    /// transcript that shares only its first components, is not its own.
    pub fn kept_attribute_codes(&self, id: &str, version: Option<u32>) -> Vec<&'static str> {
        let tags = version
            .and_then(|version| self.tags.get(format!("{id}.{version}").as_str()))
            .or_else(|| self.tags.get(id));
        let Some(tags) = tags else {
            return Vec::new();
        };
        let mut codes: Vec<&'static str> = Vec::new();
        for (_, code) in KEPT_TAGS.iter() {
            if tags.iter().any(|t| t == code) && !codes.contains(code) {
                codes.push(code);
            }
        }
        codes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gff3::TranscriptRecord;

    /// Parses `text` (a GTF, header lines included) from a temporary file
    /// under `reading`.
    fn try_parse_model_as(text: &str, reading: Reading) -> Result<GtfModel> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf");
        std::fs::write(&path, text).unwrap();
        parse_gtf_model(&path, reading)
    }

    fn parse_model_as(text: &str, reading: Reading) -> GtfModel {
        try_parse_model_as(text, reading).unwrap()
    }

    /// Parses `text` as Ensembl VEP's `--gtf` source reads it.
    fn parse_model(text: &str) -> GtfModel {
        parse_model_as(text, Reading::GffSource)
    }

    fn transcript<'a>(genes: &'a HashMap<String, GeneRecord>, id: &str) -> &'a TranscriptRecord {
        genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .find(|t| t.transcript_id == id)
            .unwrap_or_else(|| panic!("transcript {id} not parsed"))
    }

    /// Seven lines in the shape of Homo_sapiens.GRCh38.116.gtf.gz behind a
    /// `#!genome-build` header: a gene, its transcript, two exons, two CDS lines and
    /// the stop codon, which the Ensembl GTF writes as its own line outside the CDS.
    const SEVEN_LINES: &str = "#!genome-build GRCh38.p14\n\
21\tensembl_havana\tgene\t1000\t2000\t.\t+\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; gene_name \"GENE1\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\";\n\
21\tensembl_havana\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; transcript_id \"ENST00000000001\"; transcript_version \"5\"; gene_name \"GENE1\"; transcript_name \"GENE1-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"Ensembl_canonical\"; tag \"gencode_basic\"; tag \"MANE_Select\"; transcript_support_level \"1\";\n\
21\tensembl_havana\texon\t1000\t1200\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"1\"; exon_id \"ENSE00000000001\"; exon_version \"1\"; transcript_biotype \"protein_coding\"; tag \"Ensembl_canonical\";\n\
21\tensembl_havana\texon\t1500\t2000\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; exon_id \"ENSE00000000002\"; exon_version \"1\"; transcript_biotype \"protein_coding\";\n\
21\tensembl_havana\tCDS\t1100\t1200\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"1\"; protein_id \"ENSP00000000001\"; protein_version \"2\"; transcript_biotype \"protein_coding\";\n\
21\tensembl_havana\tCDS\t1500\t1600\t.\t+\t1\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; protein_id \"ENSP00000000001\"; protein_version \"2\"; transcript_biotype \"protein_coding\";\n\
21\tensembl_havana\tstop_codon\t1601\t1603\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; protein_id \"ENSP00000000001\"; transcript_biotype \"protein_coding\";\n";

    /// A `gene` line is a gene, a `transcript` line a transcript of the gene
    /// its `gene_id` names, and an `exon`, `CDS` or `stop_codon` line belongs
    /// to the transcript its `transcript_id` names (`GTF.pm` lines 71-85,
    /// `%PARENTS` and `%INCLUDE_FEATURE_TYPES`): one gene with one
    /// transcript, two exons and three CDS segments, the stop codon among
    /// them (`BaseGXF.pm` lines 510-511, `cds` and `stop_codon` alike), so
    /// the coding span runs to 1603. The transcript keeps its
    /// `transcript_id` and `gene_id`, `gene_name`, `transcript_biotype`, the
    /// `Ensembl_canonical`, `gencode_basic` and `MANE_Select` tags, the
    /// support level, the `protein_id` and, with no `version` attribute on the
    /// CDS lines, protein version 1 (`_add_translation`, line 730). The gene
    /// line's `gene_version` is not the `version` attribute Ensembl VEP reads,
    /// so the gene has no version, and `start_codon`, `five_prime_utr` and
    /// `three_prime_utr` lines are skipped.
    #[test]
    fn gtf_transcript_model_builds_the_gene_tree_by_transcript_and_gene_id() {
        let model = parse_model(&format!(
            "{SEVEN_LINES}21\tensembl_havana\tstart_codon\t1100\t1102\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\";\n\
21\tensembl_havana\tfive_prime_utr\t1000\t1099\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\";\n"
        ));
        assert_eq!(model.header["genome-build"], "GRCh38.p14");
        assert_eq!(model.genes.len(), 1);
        let gene = &model.genes["ENSG00000000001"];
        assert_eq!(gene.name, None);
        assert_eq!(gene.version, None);
        assert_eq!(gene.biotype, "protein_coding");
        assert_eq!((gene.start, gene.end, gene.strand), (1000, 2000, 1));
        assert_eq!(gene.transcripts.len(), 1);
        let tr = &gene.transcripts[0];
        assert_eq!(tr.transcript_id, "ENST00000000001");
        assert_eq!(tr.gene_id, "ENSG00000000001");
        assert_eq!(tr.symbol.as_deref(), Some("GENE1"));
        assert_eq!(tr.biotype, "protein_coding");
        assert_eq!(tr.source, "ensembl_havana");
        assert_eq!(tr.version, Some(5));
        assert!(tr.is_canonical && tr.gencode_basic && tr.mane_select);
        assert!(!tr.mane_plus_clinical && !tr.gencode_primary);
        assert_eq!(tr.tsl.as_deref(), Some("1"));
        assert_eq!(tr.ccds, None);
        assert_eq!(
            tr.exons
                .iter()
                .map(|e| (e.start, e.end, e.exon_id.as_deref(), e.rank))
                .collect::<Vec<_>>(),
            [
                (1000, 1200, Some("ENSE00000000001"), None),
                (1500, 2000, Some("ENSE00000000002"), None),
            ]
        );
        assert_eq!(
            tr.cds_regions
                .iter()
                .map(|c| (c.start, c.end, c.phase))
                .collect::<Vec<_>>(),
            [(1100, 1200, 0), (1500, 1600, 1), (1601, 1603, 0)]
        );
        for cds in &tr.cds_regions {
            assert_eq!(cds.protein_id.as_deref(), Some("ENSP00000000001"));
            assert_eq!(cds.protein_version, Some(1));
        }
    }

    /// Read for the cache's content, a CDS line's `protein_version`, the
    /// attribute the Ensembl GTF writes the translation's version as, is the
    /// version (2 on the two CDS lines; the stop codon line, which the
    /// Ensembl GTF writes without one, keeps none), where the source reading
    /// leaves it at 1; without the attribute the lines keep no version, so the
    /// build takes the core table dumps' or the protein FASTA header's; a
    /// `version` attribute stands under either reading.
    #[test]
    fn cache_content_reading_takes_a_gtf_cds_line_protein_version_else_none() {
        let versions = |text: &str, reading| {
            parse_model_as(text, reading).genes["ENSG00000000001"].transcripts[0]
                .cds_regions
                .iter()
                .map(|c| c.protein_version)
                .collect::<Vec<_>>()
        };
        assert_eq!(
            versions(SEVEN_LINES, Reading::EnsemblCache),
            [Some(2), Some(2), None]
        );
        let unversioned = SEVEN_LINES.replace(" protein_version \"2\";", "");
        assert!(!unversioned.contains("protein_version"));
        assert_eq!(
            versions(&unversioned, Reading::EnsemblCache),
            [None, None, None]
        );
        assert_eq!(
            versions(&unversioned, Reading::GffSource),
            [Some(1), Some(1), Some(1)]
        );
        let versioned = SEVEN_LINES.replace(
            "protein_id \"ENSP00000000001\";",
            "protein_id \"ENSP00000000001\"; version \"3\";",
        );
        assert_eq!(versioned.matches("\"; version \"3\";").count(), 3);
        for reading in [Reading::GffSource, Reading::EnsemblCache] {
            assert_eq!(
                versions(&versioned, reading),
                [Some(3), Some(3), Some(3)],
                "{reading:?}"
            );
        }
    }

    /// A transcript line without `transcript_biotype`, `transcript_type` or
    /// `biotype` takes the source column as its biotype (`GTF.pm` lines
    /// 176-182), `havana` here; the three attributes are read in that order.
    #[test]
    fn gtf_biotype_falls_back_to_the_source_column() {
        let no_biotype = SEVEN_LINES
            .replace(" transcript_biotype \"protein_coding\";", "")
            .replace("21\tensembl_havana\ttranscript", "21\thavana\ttranscript");
        let model = parse_model(&no_biotype);
        assert_eq!(
            transcript(&model.genes, "ENST00000000001").biotype,
            "havana"
        );

        let typed = no_biotype.replace(
            "transcript_name \"GENE1-201\";",
            "transcript_name \"GENE1-201\"; biotype \"lncRNA\"; transcript_type \"retained_intron\";",
        );
        let model = parse_model(&typed);
        assert_eq!(
            transcript(&model.genes, "ENST00000000001").biotype,
            "retained_intron"
        );
    }

    /// The stop codon's line is one more CDS segment (`BaseGXF.pm` lines
    /// 510-511), so the coding span ends on its last base, 1603, as it does
    /// on an Ensembl GFF3, whose CDS lines include the stop codon; a stop
    /// codon split over two exons adds a segment in each.
    #[test]
    fn stop_codon_lines_extend_the_cds() {
        let model = parse_model(SEVEN_LINES);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert_eq!(tr.cds_regions.iter().map(|c| c.end).max(), Some(1603));

        let without_stop = SEVEN_LINES
            .lines()
            .filter(|l| !l.contains("\tstop_codon\t"))
            .collect::<Vec<_>>()
            .join("\n");
        let model = parse_model(&without_stop);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert_eq!(tr.cds_regions.iter().map(|c| c.end).max(), Some(1600));

        let split = format!(
            "{without_stop}\n21\tensembl_havana\tstop_codon\t1601\t1601\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\";\n\
21\tensembl_havana\texon\t1900\t2000\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_id \"ENSE00000000003\";\n\
21\tensembl_havana\tstop_codon\t1900\t1901\t.\t+\t2\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\";\n"
        );
        let model = parse_model(&split);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert_eq!(
            tr.cds_regions
                .iter()
                .map(|c| (c.start, c.end, c.phase))
                .collect::<Vec<_>>(),
            [
                (1100, 1200, 0),
                (1500, 1600, 1),
                (1601, 1601, 0),
                (1900, 1901, 2)
            ]
        );
    }

    /// A `cds` line, the lowercase spelling [`GTF_FEATURE_TYPES`] lists beside
    /// `CDS`, is a CDS segment of the transcript its `transcript_id` names, as
    /// a `CDS` line is: [`SEVEN_LINES`] with its two `CDS` lines typed `cds`
    /// gives the model of [`SEVEN_LINES`], three segments with the stop codon's,
    /// under either reading.
    #[test]
    fn a_lowercase_cds_line_is_a_cds_segment_as_a_cds_line_is() {
        let lowercase = SEVEN_LINES.replace("\tCDS\t", "\tcds\t");
        assert_eq!(lowercase.matches("\tcds\t").count(), 2);
        for reading in [Reading::GffSource, Reading::EnsemblCache] {
            let model = parse_model_as(&lowercase, reading);
            assert_eq!(
                transcript(&model.genes, "ENST00000000001")
                    .cds_regions
                    .len(),
                3,
                "{reading:?}"
            );
            assert_eq!(
                format!("{:?}", model.genes),
                format!("{:?}", parse_model_as(SEVEN_LINES, reading).genes),
                "{reading:?}"
            );
        }
    }

    /// A transcript's version is its `version` attribute, else its
    /// `transcript_version` (`BaseGXF.pm` line 481), the GTF spelling: 5
    /// here, 7 when a `version` attribute is present as well, and none when
    /// neither is.
    #[test]
    fn transcript_version_attribute_is_the_fallback_version() {
        let model = parse_model(SEVEN_LINES);
        assert_eq!(transcript(&model.genes, "ENST00000000001").version, Some(5));

        let both = SEVEN_LINES.replace(
            "transcript_version \"5\";",
            "transcript_version \"5\"; version \"7\";",
        );
        let model = parse_model(&both);
        assert_eq!(transcript(&model.genes, "ENST00000000001").version, Some(7));

        let neither = SEVEN_LINES.replace(" transcript_version \"5\";", "");
        let model = parse_model(&neither);
        assert_eq!(transcript(&model.genes, "ENST00000000001").version, None);
    }

    /// A `tag "readthrough_transcript"` marks the transcript `readthrough`
    /// and a `transcript_biotype "artifact"` is its biotype, each carried on
    /// the record as Ensembl VEP's `--gtf` source carries them; the model
    /// keeps the transcript either way.
    #[test]
    fn readthrough_tag_and_artifact_biotype_are_carried() {
        let model = parse_model(SEVEN_LINES);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert!(!tr.readthrough);
        assert_eq!(tr.biotype, "protein_coding");

        let tagged = SEVEN_LINES.replace(
            "tag \"MANE_Select\";",
            "tag \"MANE_Select\"; tag \"readthrough_transcript\";",
        );
        let model = parse_model(&tagged);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert!(tr.readthrough && tr.mane_select && tr.is_canonical);

        let artifact = SEVEN_LINES.replacen(
            "transcript_biotype \"protein_coding\"",
            "transcript_biotype \"artifact\"",
            1,
        );
        let model = parse_model(&artifact);
        assert_eq!(
            transcript(&model.genes, "ENST00000000001").biotype,
            "artifact"
        );
    }

    /// The GENCODE basic tag is `tag "basic"` in the GRCh37 release-87 and
    /// GENCODE GTFs and `tag "gencode_basic"` in the GRCh38 Ensembl one; the
    /// model reads both as the `gencode_basic` attribute, the spelling the
    /// tags beside a GFF3 give it ([`KEPT_TAGS`]), so `--gencode_basic` on
    /// the built cache keeps a `basic` transcript as it keeps a
    /// `gencode_basic` one. A transcript tagged neither has no attribute,
    /// and `gencode_primary` is read under its one spelling.
    #[test]
    fn basic_tag_is_the_gencode_basic_attribute_on_the_model() {
        let model = parse_model(SEVEN_LINES);
        assert!(transcript(&model.genes, "ENST00000000001").gencode_basic);

        let grch37 = SEVEN_LINES.replace("tag \"gencode_basic\";", "tag \"basic\";");
        assert!(!grch37.contains("gencode_basic"));
        let model = parse_model(&grch37);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert!(tr.gencode_basic && tr.is_canonical && tr.mane_select);
        assert!(!tr.gencode_primary);

        let untagged = SEVEN_LINES.replace(" tag \"gencode_basic\";", "");
        let model = parse_model(&untagged);
        assert!(!transcript(&model.genes, "ENST00000000001").gencode_basic);

        let primary = SEVEN_LINES.replace("tag \"gencode_basic\";", "tag \"gencode_primary\";");
        let model = parse_model(&primary);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert!(tr.gencode_primary && !tr.gencode_basic);
    }

    /// A transcript line's `ccds_id` is the transcript's CCDS identifier, as
    /// the line carries it: the Ensembl GTF writes the xref's primary
    /// identifier, `CCDS33512` with no version (ensembl-io `GTFSerializer.pm`
    /// lines 459-464), where the GFF3's `ccdsid` is `CCDS33512.1`. A line
    /// without the attribute, or with the value `0`, gives none.
    #[test]
    fn ccds_id_becomes_ccds() {
        let with_ccds = SEVEN_LINES.replace(
            "tag \"gencode_basic\";",
            "tag \"gencode_basic\"; tag \"CCDS\"; ccds_id \"CCDS33512\";",
        );
        let model = parse_model(&with_ccds);
        assert_eq!(
            transcript(&model.genes, "ENST00000000001").ccds.as_deref(),
            Some("CCDS33512")
        );

        let model = parse_model(SEVEN_LINES);
        assert_eq!(transcript(&model.genes, "ENST00000000001").ccds, None);

        let zero = with_ccds.replace("ccds_id \"CCDS33512\";", "ccds_id \"0\";");
        let model = parse_model(&zero);
        assert_eq!(transcript(&model.genes, "ENST00000000001").ccds, None);
    }

    /// A transcript whose `gene_id` names no `gene` line is a child of a
    /// record in no line, which Ensembl VEP drops with its `Parent entries
    /// ... not found` warning, as is an exon whose `transcript_id` names no
    /// `transcript` line; a `transcript_id` list places the exon under each
    /// transcript it names.
    #[test]
    fn gtf_records_without_their_parent_line_are_dropped() {
        let no_gene = SEVEN_LINES
            .lines()
            .filter(|l| !l.contains("\tgene\t"))
            .collect::<Vec<_>>()
            .join("\n");
        let model = parse_model(&no_gene);
        assert!(model.genes.is_empty());

        let shared_exon = format!(
            "{SEVEN_LINES}21\tensembl_havana\ttranscript\t1000\t1200\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000002\"; transcript_biotype \"retained_intron\";\n\
21\tensembl_havana\texon\t1000\t1200\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000002,ENST00000000009\"; exon_id \"ENSE00000000001\";\n"
        );
        let model = parse_model(&shared_exon);
        let second = transcript(&model.genes, "ENST00000000002");
        assert_eq!(second.exons.len(), 1);
        assert!(second.cds_regions.is_empty());
        assert_eq!(model.genes["ENSG00000000001"].transcripts.len(), 2);
    }

    /// A GTF attribute whose value is `0` is absent, as Ensembl VEP's truth
    /// tests read it: `T1`'s `transcript_biotype "0"` and `transcript_type
    /// "0"` yield to `biotype` (`GTF.pm` lines 178-181), and its CDS line's
    /// `protein_id "0"` is no protein (`BaseGXF.pm` line 723) and `version
    /// "0"` version 1 (line 727). A transcript line's `gene_id "0"` is a
    /// parent list naming no set identifier (`GTF.pm` line 137 splits the
    /// value as written): the record is linked to nothing (`BaseGXF.pm`
    /// line 389, `next unless ($parent_id)`) and, having a parent list, is
    /// no record of its own (lines 430-435), so `T0` is dropped. A
    /// `transcript_id "0"` names a record by the MD5 of its line (`GTF.pm`
    /// line 158, `<type>_id || md5`) that no exon's `transcript_id "0"`
    /// reaches, so that transcript of `G2` has no children and is dropped
    /// (`BaseGXF.pm` line 462), and `G3`'s transcript, with no biotype
    /// attribute and a source column of `0`, has no biotype and is dropped
    /// (line 470). A `gene_id ",0"` opens with an empty field, so the
    /// prune leaves `T4` a transcript of its own (line 432, `$parent_ids->[0]
    /// ne ''`), `T5`'s `gene_id ",G1"` places it under `G1` once, and
    /// `T6`'s `gene_id ",X"`, naming no gene line, leaves it a transcript of
    /// its own with that `gene_id`, as Ensembl VEP leaves it.
    #[test]
    fn a_zero_value_is_an_absent_attribute() {
        let model = parse_model(
            "1\thavana\tgene\t1000\t2000\t.\t+\t.\tgene_id \"G1\";\n\
1\thavana\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_biotype \"0\"; transcript_type \"0\"; biotype \"lncRNA\";\n\
1\thavana\texon\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
1\thavana\tCDS\t1100\t1200\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\"; protein_id \"0\"; version \"0\";\n\
1\thavana\ttranscript\t3000\t4000\t.\t+\t.\tgene_id \"0\"; transcript_id \"T0\"; transcript_biotype \"lncRNA\";\n\
1\thavana\texon\t3000\t4000\t.\t+\t.\tgene_id \"0\"; transcript_id \"T0\";\n\
1\thavana\tgene\t5000\t6000\t.\t+\t.\tgene_id \"G2\";\n\
1\thavana\ttranscript\t5000\t6000\t.\t+\t.\tgene_id \"G2\"; transcript_id \"0\"; transcript_biotype \"lncRNA\";\n\
1\thavana\texon\t5000\t6000\t.\t+\t.\tgene_id \"G2\"; transcript_id \"0\";\n\
1\thavana\tgene\t7000\t8000\t.\t+\t.\tgene_id \"G3\";\n\
1\t0\ttranscript\t7000\t8000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";\n\
1\thavana\texon\t7000\t8000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\";\n\
1\thavana\ttranscript\t9000\t9900\t.\t+\t.\tgene_id \",0\"; transcript_id \"T4\"; transcript_biotype \"lncRNA\";\n\
1\thavana\texon\t9000\t9900\t.\t+\t.\tgene_id \",0\"; transcript_id \"T4\";\n\
1\thavana\ttranscript\t1000\t1500\t.\t+\t.\tgene_id \",G1\"; transcript_id \"T5\"; transcript_biotype \"lncRNA\";\n\
1\thavana\texon\t1000\t1500\t.\t+\t.\tgene_id \",G1\"; transcript_id \"T5\";\n\
1\thavana\ttranscript\t12000\t12900\t.\t+\t.\tgene_id \",X\"; transcript_id \"T6\"; transcript_biotype \"lncRNA\";\n\
1\thavana\texon\t12000\t12900\t.\t+\t.\tgene_id \",X\"; transcript_id \"T6\";\n",
        );
        let mut ids: Vec<&str> = model
            .genes
            .values()
            .flat_map(|g| g.transcripts.iter())
            .map(|t| t.transcript_id.as_str())
            .collect();
        ids.sort_unstable();
        assert_eq!(ids, ["T1", "T4", "T5", "T6"]);
        let tr = transcript(&model.genes, "T1");
        assert_eq!(tr.gene_id, "G1");
        assert_eq!(tr.biotype, "lncRNA");
        assert_eq!(tr.exons.len(), 1);
        assert_eq!(tr.cds_regions.len(), 1);
        assert_eq!(tr.cds_regions[0].protein_id, None);
        assert_eq!(tr.cds_regions[0].protein_version, Some(1));
        assert_eq!(model.genes.len(), 5, "{:?}", model.genes.keys());
        assert_eq!(model.genes["G1"].transcripts.len(), 2);
        assert!(model.genes["G2"].transcripts.is_empty());
        assert!(model.genes["G3"].transcripts.is_empty());
        assert_eq!(transcript(&model.genes, "T4").gene_id, ",0");
        assert_eq!(transcript(&model.genes, "T6").gene_id, ",X");
    }

    /// The attribute column read as `GTF.pm` `get_attributes` reads it (lines
    /// 117-135), each expected value the one that code gives the same column:
    /// every `%XX` escape decoded in either letter case (`A%2CB` and `A%2cB`
    /// are `A,B`), every quote removed wherever it stands (`"a"b"` is `ab`),
    /// `"%30"` the `0` the truth tests read as unset, and `"0"`, `0` and the
    /// unquoted `0` alike; the value runs to the `;`, so a space before it
    /// stays in the value (`"X" ` is `X `); a piece of one word sets nothing,
    /// and a key of `0` is skipped; a `;` inside a quoted value ends the value
    /// there and opens a piece of one word; a single `tag` is split on `,`
    /// and a repeated `tag` lists its values whole, a `tag` of one word
    /// after a set one repeating it while one before it does not; the
    /// vertical tab is whitespace to the split and a no-break space is not.
    #[test]
    fn gtf_attributes_are_read_as_the_gtf_parser_reads_them() {
        let read = |column: &str| {
            let (attrs, tags) = parse_gtf_attributes(column, true).unwrap();
            let mut attrs: Vec<(String, String)> = attrs.into_iter().collect();
            attrs.sort();
            (attrs, tags)
        };
        let pair = |k: &str, v: &str| (k.to_string(), v.to_string());
        assert_eq!(
            read("gene_id \"A%2CB\"; gene_name \"a\"b\"; x \"%30\"; y \"0\"; z 0; tag \"basic\"; tag \"x%20y\""),
            (
                vec![
                    pair("gene_id", "A,B"),
                    pair("gene_name", "ab"),
                    pair("tag", "x y"),
                    pair("x", "0"),
                    pair("y", "0"),
                    pair("z", "0"),
                ],
                vec!["gencode_basic".to_string(), "x y".to_string()]
            )
        );
        assert_eq!(read("gene_id \"A%2cB\"").0, [pair("gene_id", "A,B")]);
        assert_eq!(
            read("gene_id \"X\" ; transcript_id   \"T\"   ").0,
            [pair("gene_id", "X "), pair("transcript_id", "T   ")]
        );
        assert_eq!(
            read("gene_name \"foo;bar=baz\"; k").0,
            [pair("gene_name", "foo")]
        );
        assert_eq!(
            read("0 \"x\"; gene_id \"G\"; ;  ").0,
            [pair("gene_id", "G")]
        );
        assert_eq!(read("gene_id; gene_id \"B\"").0, [pair("gene_id", "B")]);
        assert_eq!(
            read("ont \"A\"; ont \"B\"; gene_id \"G\"").0,
            [pair("gene_id", "G"), pair("ont", "B")]
        );
        assert_eq!(read("gene_id \"%C3%A9\"").0, [pair("gene_id", "\u{E9}")]);
        assert_eq!(
            read("tag \"basic,Ensembl_canonical,\"").1,
            ["gencode_basic", "Ensembl_canonical"]
        );
        assert_eq!(
            read("tag \"basic,Ensembl_canonical\"; tag \"c\"").1,
            ["basic,Ensembl_canonical", "c"]
        );
        assert_eq!(
            read("tag \"basic,Ensembl_canonical\"; tag").1,
            ["basic,Ensembl_canonical"]
        );
        assert_eq!(
            read("tag; tag \"basic,Ensembl_canonical\"").1,
            ["gencode_basic", "Ensembl_canonical"]
        );
        assert_eq!(read("gene_id\u{0B}\"G\"").0, [pair("gene_id", "G")]);
        assert!(read("gene_id\u{A0}\"G\"").0.is_empty());
    }

    /// A single `tag` whose value is `"0"`, `""`, the unquoted `0` or `"%30"`
    /// lists nothing: `_add_identifiers` splits the one value only when it is
    /// true to Perl (`BaseGXF.pm` line 683, `if($tr_record->{attributes}->{tag})`),
    /// and `0` and the empty string are not. A repeated `tag` is a list
    /// reference, true whatever it holds, so `"0"` beside `"basic"` is listed.
    #[test]
    fn a_single_tag_that_is_not_true_to_perl_lists_nothing() {
        let tags = |column: &str| parse_gtf_attributes(column, true).unwrap().1;
        for column in [
            "gene_id \"G\"; tag \"0\"",
            "gene_id \"G\"; tag \"\"",
            "gene_id \"G\"; tag 0",
            "gene_id \"G\"; tag \"%30\"",
        ] {
            let listed = tags(column);
            assert!(listed.is_empty(), "{column}: {listed:?}");
        }
        assert_eq!(tags("tag \"0\"; tag \"basic\""), ["0", "gencode_basic"]);
        assert_eq!(tags("tag \"\"; tag \"basic\""), ["", "gencode_basic"]);
        assert_eq!(tags("tag \"basic\"; tag \"0\""), ["gencode_basic", "0"]);
    }

    /// A line read without its tags holds the attributes and draws the
    /// refusals of the same line read with them, its `tag` value among the
    /// attributes; the list alone is left unbuilt.
    #[test]
    fn a_line_read_without_its_tags_keeps_its_attributes_and_refusals() {
        let column = "gene_id \"G\"; tag \"basic\"; tag \"CCDS\"; transcript_id \"T\"";
        let (attrs, tags) = parse_gtf_attributes(column, true).unwrap();
        let (untagged, none) = parse_gtf_attributes(column, false).unwrap();
        assert_eq!(tags, ["gencode_basic", "CCDS"]);
        assert!(none.is_empty());
        assert_eq!(attrs, untagged);
        assert_eq!(untagged["tag"], "CCDS");
        for refused in [
            "tag \"%C3%zz%2\"",
            "tag \"a\"; gene_id \"A\"; gene_id \"B\"",
        ] {
            assert_eq!(
                parse_gtf_attributes(refused, false)
                    .unwrap_err()
                    .to_string(),
                parse_gtf_attributes(refused, true).unwrap_err().to_string(),
                "{refused}"
            );
        }
    }

    /// A transcript's tags are its own line's. The tags of an exon, CDS or
    /// stop codon line reach no record, so a file whose part lines carry
    /// tags their transcript line lacks gives the model of the file without
    /// them.
    #[test]
    fn a_transcripts_tags_are_read_from_its_own_line_alone() {
        let tagged_parts: String = SEVEN_LINES
            .lines()
            .map(|line| {
                let part = ["\texon\t", "\tCDS\t", "\tstop_codon\t"]
                    .iter()
                    .any(|feature_type| line.contains(feature_type));
                if part {
                    format!(
                        "{line} tag \"readthrough_transcript\"; tag \"gencode_primary\"; tag \"MANE_Plus_Clinical\";\n"
                    )
                } else {
                    format!("{line}\n")
                }
            })
            .collect();
        assert_eq!(tagged_parts.matches("readthrough_transcript").count(), 5);
        let model = parse_model(&tagged_parts);
        let tr = transcript(&model.genes, "ENST00000000001");
        assert!(tr.is_canonical && tr.gencode_basic && tr.mane_select);
        assert!(!tr.readthrough && !tr.gencode_primary && !tr.mane_plus_clinical);
        assert_eq!(
            format!("{:?}", model.genes),
            format!("{:?}", parse_model(SEVEN_LINES).genes)
        );
    }

    /// The two real-world GTF shapes that repeat a key are read without a
    /// refusal: a GENCODE transcript line repeats `tag` and `ont`, a RefSeq
    /// line `db_xref`. `tag` is the list `BaseGXF.pm` lines 683-685 read, so
    /// a repeated `tag` lists its values whole and `basic` is the GENCODE
    /// basic code; `ont` and `db_xref` are keys no reader takes, so each keeps
    /// its last value and nothing is refused. The identifiers, the symbol
    /// (`gene_name`, else `gene`, lines 662-666), the biotype
    /// (`transcript_type` or `transcript_biotype`) and the support level are
    /// the one value each line writes, as Ensembl VEP reads them.
    #[test]
    fn gencode_and_refseq_shaped_lines_are_read_without_refusal() {
        let gencode = "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L16\"; level 2; tag \"overlaps_pseudogene\";\n\
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825.2\"; transcript_id \"ENST00000456328.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L16\"; transcript_type \"lncRNA\"; transcript_name \"DDX11L16-264\"; level 2; transcript_support_level \"1\"; tag \"basic\"; tag \"Ensembl_canonical\"; havana_transcript \"OTTHUMT00000362751.1\"; ont \"PGO:0000005\"; ont \"PGO:0000019\";\n\
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id \"ENSG00000290825.2\"; transcript_id \"ENST00000456328.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L16\"; transcript_type \"lncRNA\"; exon_number 1; exon_id \"ENSE00002234944.1\"; level 2; tag \"basic\"; tag \"Ensembl_canonical\"; ont \"PGO:0000005\"; ont \"PGO:0000019\";\n";
        let genes = parse_model(gencode).genes;
        let tr = transcript(&genes, "ENST00000456328.2");
        assert_eq!(tr.gene_id, "ENSG00000290825.2");
        assert_eq!(tr.symbol.as_deref(), Some("DDX11L16"));
        assert_eq!(tr.biotype, "lncRNA");
        assert_eq!(tr.tsl.as_deref(), Some("1"));
        assert!(tr.is_canonical && tr.gencode_basic && !tr.mane_select);
        assert_eq!(tr.exons.len(), 1);
        let (attrs, tags) = parse_gtf_attributes(
            "tag \"basic\"; tag \"Ensembl_canonical\"; ont \"PGO:0000005\"; ont \"PGO:0000019\"",
            true,
        )
        .unwrap();
        assert_eq!(tags, ["gencode_basic", "Ensembl_canonical"]);
        assert_eq!(attrs["ont"], "PGO:0000019");

        let refseq = "NC_000001.11\tBestRefSeq\tgene\t11874\t14409\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"\"; db_xref \"GeneID:100287102\"; db_xref \"HGNC:HGNC:37102\"; description \"DEAD/H-box helicase 11 like 1 (pseudogene)\"; gbkey \"Gene\"; gene \"DDX11L1\"; gene_biotype \"transcribed_pseudogene\"; pseudo \"true\";\n\
NC_000001.11\tBestRefSeq\ttranscript\t11874\t14409\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018.2\"; db_xref \"GeneID:100287102\"; gbkey \"misc_RNA\"; gene \"DDX11L1\"; product \"DEAD/H-box helicase 11 like 1 (pseudogene)\"; pseudo \"true\"; transcript_biotype \"transcript\";\n\
NC_000001.11\tBestRefSeq\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018.2\"; db_xref \"GeneID:100287102\"; gene \"DDX11L1\"; product \"DEAD/H-box helicase 11 like 1 (pseudogene)\"; pseudo \"true\"; transcript_biotype \"transcript\"; exon_number \"1\";\n";
        let genes = parse_model(refseq).genes;
        let tr = transcript(&genes, "NR_046018.2");
        assert_eq!(tr.gene_id, "DDX11L1");
        assert_eq!(tr.symbol.as_deref(), Some("DDX11L1"));
        assert_eq!(tr.biotype, "transcript");
        assert_eq!(tr.chr, "NC_000001.11");
        let (attrs, _) = parse_gtf_attributes(
            "db_xref \"GeneID:100287102\"; db_xref \"HGNC:HGNC:37102\"; gene \"DDX11L1\"",
            true,
        )
        .unwrap();
        assert_eq!(attrs["db_xref"], "HGNC:HGNC:37102");
    }

    /// A line repeating a key the readers take as one value is refused,
    /// naming the key, whatever the values and even when the second
    /// occurrence has none, since Ensembl VEP's parser lists them from the
    /// first defined value (`GTF.pm` lines 128-131) and the source reads the
    /// list as one value; a repeated key no reader takes, GENCODE's `ont`,
    /// passes. An escape spelling a byte that is not UTF-8 is refused as well.
    /// Both refusals name the file and the line through the model, and the
    /// tags reader refuses the same line the same way.
    #[test]
    fn a_repeated_read_key_and_a_non_utf8_escape_are_refused() {
        for column in [
            "gene_id \"A\"; gene_id \"B\"",
            "gene_id \"A\"; gene_id \"A\"",
            "gene_id \"A\"; gene_id",
            "transcript_id \"T\"; gene_name \"n\"; gene_name \"m\"",
        ] {
            let error = parse_gtf_attributes(column, true).unwrap_err().to_string();
            assert!(error.contains("is repeated"), "{column}: {error}");
        }
        assert!(parse_gtf_attributes("ont \"A\"; ont \"B\"; tag \"a\"; tag \"a\"", true).is_ok());
        let error = parse_gtf_attributes("gene_id \"%C3%zz%2\"", true)
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "GTF attribute gene_id \"%C3%zz%2\" decodes to text that is not UTF-8"
        );
        for (text, refusal) in [
            (
                "1\thavana\tgene\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; gene_id \"G1\";\n",
                "line 1: GTF attribute gene_id is repeated; Ensembl VEP reads a repeated key as a list",
            ),
            (
                "1\thavana\tgene\t1000\t2000\t.\t+\t.\tgene_id \"G1\";\n\
1\thavana\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"%FF\";\n",
                "line 2: GTF attribute gene_name \"%FF\" decodes to text that is not UTF-8",
            ),
        ] {
            let err = try_parse_model_as(text, Reading::GffSource).unwrap_err();
            let message = format!("{err:#}");
            assert!(message.contains(refusal), "{message}");
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf");
        std::fs::write(
            &path,
            "1\thavana\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_id \"T2\"; tag \"basic\";\n",
        )
        .unwrap();
        let message = format!("{:#}", parse_gtf_tags(&path).unwrap_err());
        assert!(
            message.contains("line 1: GTF attribute transcript_id is repeated"),
            "{message}"
        );
    }

    /// Through the model: a `gene_name` of `"GENE%2C1"` is the symbol
    /// `GENE,1`, one of `"G"E"NE"` the symbol `GENE`, a
    /// `transcript_support_level` of `"%30"` the unset `0`, so no support
    /// level, and a `gene_name` of `"foo;bar=baz"` the symbol `foo`, the
    /// text before the `;` that ends the value for both parsers.
    #[test]
    fn escaped_and_requoted_gtf_values_reach_the_transcript_decoded() {
        let symbol_and_tsl = |attributes: &str| {
            let model = parse_model(&format!(
                "1\thavana\tgene\t1000\t2000\t.\t+\t.\tgene_id \"G1\";\n\
1\thavana\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_biotype \"lncRNA\"; {attributes}\n\
1\thavana\texon\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
            ));
            let tr = transcript(&model.genes, "T1");
            (tr.symbol.clone(), tr.tsl.clone())
        };
        assert_eq!(
            symbol_and_tsl("gene_name \"GENE%2C1\"; transcript_support_level \"%30\";"),
            (Some("GENE,1".to_string()), None)
        );
        assert_eq!(
            symbol_and_tsl("gene_name \"G\"E\"NE\"; transcript_support_level \"%31\";"),
            (Some("GENE".to_string()), Some("1".to_string()))
        );
        assert_eq!(
            symbol_and_tsl("gene_name \"foo;bar=baz\";"),
            (Some("foo".to_string()), None)
        );
    }

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
             1\thavana\texon\t65419\t65433\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; tag \"cds_start_NF\";\n\
             \x20\x20#\tx\ttranscript\t1\t2\t.\t+\t.\ttranscript_id \"T1\"; tag \"basic\";\n\
             1\thavana\ttranscript\t95000\t96000\t.\t+\t.\tgene_id \"ENSG00000000004\"; transcript_id \"ENST00000000004\"; tag \"basic\" \n",
        )
        .unwrap();
        let gtf = parse_gtf_tags(&path).unwrap();
        assert_eq!(gtf.header["genome-build"], "GRCh38.p14");
        assert_eq!(gtf.header["genebuild-last-updated"], "2024-01");
        assert!(
            !gtf.tags.contains_key("T1"),
            "a line of blanks then `#` is a comment, as it is to the model reader"
        );
        // A value that ends the line with a blank and no `;` loses the blank
        // to the line trim, as it does in the model reader.
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000000004", None),
            vec!["gencode_basic"]
        );
        // A line is named by its identifier and its version: a model that
        // writes the version inside its identifier and one that states it
        // apart both find the line, and a model of another version does not.
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000456328.2", None),
            vec!["gencode_basic"]
        );
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000456328", Some(2)),
            vec!["gencode_basic"]
        );
        assert!(gtf
            .kept_attribute_codes("ENST00000456328", Some(3))
            .is_empty());
        // mRNA_start_NF is not an attribute Ensembl VEP keeps; exon lines do
        // not count.
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000641515", None),
            vec!["cds_start_NF", "cds_end_NF"]
        );
        // The GRCh38 GTF spells the basic tag `gencode_basic`.
        assert_eq!(
            gtf.kept_attribute_codes("ENST00000000003", Some(7)),
            vec!["gencode_basic", "gencode_primary"]
        );
        assert!(gtf.kept_attribute_codes("ENST00000000000", None).is_empty());
    }

    /// Reads the tags of a GTF with one `transcript` line of gene `G1` per
    /// entry of `lines`, each entry the attributes behind the line's
    /// `gene_id`.
    fn tags_of_transcript_lines(lines: &[&str]) -> GtfTags {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf");
        let text: String = lines
            .iter()
            .map(|attributes| {
                format!("1\thavana\ttranscript\t100\t200\t.\t+\t.\tgene_id \"G1\"; {attributes}\n")
            })
            .collect();
        std::fs::write(&path, text).unwrap();
        parse_gtf_tags(&path).unwrap()
    }

    /// A `transcript` line's tags are held under its `transcript_id` whole,
    /// every dot kept, with its version behind one more dot when the line
    /// carries one (`02` is 2): transcripts that share a first component
    /// (`T1.1`, `T1.2`) each hold their own tags. A version that is not a
    /// whole number, or is `0`, names nothing.
    #[test]
    fn tags_are_held_under_each_lines_whole_name() {
        let gtf = tags_of_transcript_lines(&[
            "transcript_id \"T1.1\"; tag \"cds_start_NF\";",
            "transcript_id \"T1.2\"; tag \"cds_end_NF\";",
            "transcript_id \"ENST1\"; transcript_version \"02\"; tag \"basic\";",
            "transcript_id \"ENST2\"; transcript_version \"v3\"; tag \"basic\";",
            "transcript_id \"ENST3\"; transcript_version \"0\"; tag \"basic\";",
            "transcript_id \"ENST4\"; version \"5\"; transcript_version \"6\"; tag \"basic\";",
        ]);
        let mut held: Vec<(String, Vec<String>)> = gtf.tags.into_iter().collect();
        held.sort();
        let codes = |name: &str, code: &str| (name.to_string(), vec![code.to_string()]);
        assert_eq!(
            held,
            [
                codes("ENST1.2", "gencode_basic"),
                codes("ENST2", "gencode_basic"),
                codes("ENST3", "gencode_basic"),
                codes("ENST4.5", "gencode_basic"),
                codes("T1.1", "cds_start_NF"),
                codes("T1.2", "cds_end_NF"),
            ]
        );
    }

    /// A model's transcript takes the tags of the line that names it whole:
    /// a transcript whose identifier holds a dot of its own takes its own
    /// line's and never a sibling's; one that states a version takes the line
    /// of that version, written apart or inside the identifier, else the
    /// line that carries no version, and never the line of another version;
    /// one that states none does not take a line that carries one.
    #[test]
    fn a_transcript_takes_the_tags_of_the_line_that_names_it_whole() {
        let gtf = tags_of_transcript_lines(&[
            "transcript_id \"T1.1\"; tag \"cds_start_NF\";",
            "transcript_id \"T1.2\"; tag \"cds_end_NF\";",
            "transcript_id \"ENST1\"; transcript_version \"3\"; tag \"basic\";",
            "transcript_id \"ENST2.3\"; tag \"gencode_primary\";",
            "transcript_id \"ENST3\"; tag \"cds_start_NF\";",
            "transcript_id \"ENST4\"; transcript_version \"3\"; tag \"basic\";",
            "transcript_id \"ENST4\"; tag \"cds_end_NF\";",
        ]);
        let none: Vec<&str> = Vec::new();

        assert_eq!(gtf.kept_attribute_codes("T1.1", None), ["cds_start_NF"]);
        assert_eq!(gtf.kept_attribute_codes("T1.2", None), ["cds_end_NF"]);
        assert_eq!(gtf.kept_attribute_codes("T1", None), none);
        assert_eq!(gtf.kept_attribute_codes("T1", Some(1)), ["cds_start_NF"]);

        assert_eq!(
            gtf.kept_attribute_codes("ENST1", Some(3)),
            ["gencode_basic"]
        );
        assert_eq!(gtf.kept_attribute_codes("ENST1", Some(4)), none);
        assert_eq!(gtf.kept_attribute_codes("ENST1", None), none);

        assert_eq!(
            gtf.kept_attribute_codes("ENST2", Some(3)),
            ["gencode_primary"]
        );
        assert_eq!(gtf.kept_attribute_codes("ENST2", Some(4)), none);
        assert_eq!(
            gtf.kept_attribute_codes("ENST2.3", None),
            ["gencode_primary"]
        );

        assert_eq!(gtf.kept_attribute_codes("ENST3", Some(9)), ["cds_start_NF"]);
        assert_eq!(gtf.kept_attribute_codes("ENST3", None), ["cds_start_NF"]);

        assert_eq!(
            gtf.kept_attribute_codes("ENST4", Some(3)),
            ["gencode_basic"]
        );
        assert_eq!(gtf.kept_attribute_codes("ENST4", Some(4)), ["cds_end_NF"]);
    }

    /// A bgzipped GTF is a run of gzip members, the shape Ensembl VEP's
    /// `--gtf` takes: every member is read, so the model and the tags are
    /// the plain file's, not the header block's alone.
    #[test]
    fn bgzipped_gtf_is_read_past_its_first_member() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf.gz");
        let mut bytes = Vec::new();
        let (header, body) = SEVEN_LINES.split_once('\n').unwrap();
        for member in [format!("{header}\n"), body.to_string()] {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(member.as_bytes()).unwrap();
            bytes.extend(encoder.finish().unwrap());
        }
        std::fs::write(&path, bytes).unwrap();
        let model = parse_gtf_model(&path, Reading::GffSource).unwrap();
        assert_eq!(model.header["genome-build"], "GRCh38.p14");
        assert_eq!(model.genes["ENSG00000000001"].transcripts.len(), 1);
        let tags = parse_gtf_tags(&path).unwrap();
        assert_eq!(
            tags.kept_attribute_codes("ENST00000000001", Some(5)),
            vec!["gencode_basic"]
        );
    }

    /// A line one byte over the cap (a transcript line whose tag value pads it
    /// to that length) fails the parse, naming the file and the line, 2.
    #[test]
    fn line_past_the_cap_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf");
        let mut transcript_line =
            b"1\thavana\ttranscript\t11869\t14409\t.\t+\t.\ttranscript_id \"ENST00000456328\"; tag \""
                .to_vec();
        transcript_line.resize(MAX_LINE_BYTES + 1, b'x');
        let mut text = b"#!genome-build GRCh38.p14\n".to_vec();
        text.extend(transcript_line);
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        let err = parse_gtf_tags(&path).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "t.gtf line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
        let err = parse_gtf_model(&path, Reading::GffSource).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "t.gtf line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
    }

    /// A `version` whose digits do not fit the cache's 32-bit field fails the
    /// parse under either reading, naming the file, the line (counted with
    /// the header) and the limit: a gene line's `version`, a CDS line's
    /// `version`, a transcript line's `transcript_version`, and, read for the
    /// cache's content, a CDS line's `protein_version`. The source reading
    /// never reads `protein_version`, and neither reading reads
    /// `gene_version`, so those values build as without the attribute.
    #[test]
    fn a_version_above_the_cache_limit_fails_the_parse() {
        let limit = "version 4294967296 is above 4294967295, the largest version the cache holds";
        let version_on_gene = SEVEN_LINES.replacen(
            "gene_version \"3\";",
            "gene_version \"3\"; version \"4294967296\";",
            1,
        );
        let version_on_cds = SEVEN_LINES.replacen(
            "protein_id \"ENSP00000000001\";",
            "protein_id \"ENSP00000000001\"; version \"4294967296\";",
            1,
        );
        let transcript_version = SEVEN_LINES.replace(
            "transcript_version \"5\"",
            "transcript_version \"4294967296\"",
        );
        let protein_version =
            SEVEN_LINES.replacen("protein_version \"2\"", "protein_version \"4294967296\"", 1);
        let gene_version =
            SEVEN_LINES.replacen("gene_version \"3\"", "gene_version \"4294967296\"", 1);
        for (text, reading, line) in [
            (&version_on_gene, Reading::GffSource, 2),
            (&version_on_gene, Reading::EnsemblCache, 2),
            (&version_on_cds, Reading::GffSource, 6),
            (&version_on_cds, Reading::EnsemblCache, 6),
            (&transcript_version, Reading::GffSource, 3),
            (&transcript_version, Reading::EnsemblCache, 3),
            (&protein_version, Reading::EnsemblCache, 6),
        ] {
            let err = try_parse_model_as(text, reading).unwrap_err();
            let message = format!("{err:#}");
            assert!(
                message.ends_with(&format!("t.gtf line {line}: {limit}")),
                "{reading:?}: {message}"
            );
        }
        let unread = parse_model_as(&protein_version, Reading::GffSource);
        assert_eq!(
            unread.genes["ENSG00000000001"].transcripts[0]
                .cds_regions
                .iter()
                .map(|c| c.protein_version)
                .collect::<Vec<_>>(),
            [Some(1), Some(1), Some(1)]
        );
        let unread = parse_model_as(&gene_version, Reading::EnsemblCache);
        assert_eq!(unread.genes["ENSG00000000001"].version, None);
    }

    /// A `version` whose text is not a whole number fails the parse under
    /// either reading, naming the file, the line and the text: a gene line's
    /// `version`, a transcript line's `transcript_version`, a CDS line's
    /// `version` and, read for the cache's content, a CDS line's
    /// `protein_version`, where Ensembl VEP prints the text after the
    /// identifier. The source reading never reads `protein_version`, so that
    /// value builds as without the attribute; `01` is 1.
    #[test]
    fn a_version_that_is_not_a_whole_number_fails_the_parse() {
        let version_on_gene = SEVEN_LINES.replacen(
            "gene_version \"3\";",
            "gene_version \"3\"; version \"v2\";",
            1,
        );
        let version_on_cds = SEVEN_LINES.replacen(
            "protein_id \"ENSP00000000001\";",
            "protein_id \"ENSP00000000001\"; version \"v2\";",
            1,
        );
        let transcript_version =
            SEVEN_LINES.replace("transcript_version \"5\"", "transcript_version \"v2\"");
        let protein_version =
            SEVEN_LINES.replacen("protein_version \"2\"", "protein_version \"v2\"", 1);
        for (text, reading, line) in [
            (&version_on_gene, Reading::GffSource, 2),
            (&version_on_gene, Reading::EnsemblCache, 2),
            (&version_on_cds, Reading::GffSource, 6),
            (&version_on_cds, Reading::EnsemblCache, 6),
            (&transcript_version, Reading::GffSource, 3),
            (&transcript_version, Reading::EnsemblCache, 3),
            (&protein_version, Reading::EnsemblCache, 6),
        ] {
            let err = try_parse_model_as(text, reading).unwrap_err();
            let message = format!("{err:#}");
            assert!(
                message.ends_with(&format!(
                    "t.gtf line {line}: version \"v2\" is not a whole number"
                )),
                "{reading:?}: {message}"
            );
        }
        let unread = parse_model_as(&protein_version, Reading::GffSource);
        assert_eq!(
            unread.genes["ENSG00000000001"].transcripts[0]
                .cds_regions
                .iter()
                .map(|c| c.protein_version)
                .collect::<Vec<_>>(),
            [Some(1), Some(1), Some(1)]
        );
        let zeros = SEVEN_LINES.replace("transcript_version \"5\"", "transcript_version \"05\"");
        assert_eq!(
            parse_model_as(&zeros, Reading::GffSource).genes["ENSG00000000001"].transcripts[0]
                .version,
            Some(5)
        );
    }

    /// [`SEVEN_LINES`] on sequence 20, then `on_first` (line 9 of the file),
    /// the same lines on sequence 21 under their own identifiers, then
    /// `on_later` (line 17).
    fn two_sequences(on_first: &str, on_later: &str) -> String {
        let (header, lines) = SEVEN_LINES.split_once('\n').unwrap();
        let first = format!("\n{lines}")
            .replace("\n21\t", "\n20\t")
            .replace("ENSG00000000001", "ENSG00000000002")
            .replace("ENST00000000001", "ENST00000000002");
        format!("{header}{first}{on_first}{lines}{on_later}")
    }

    /// A line of type `feature_type` on sequence `chr` under the model of
    /// that sequence, with the start and end as written.
    fn line_of(chr: &str, feature_type: &str, start: &str, end: &str) -> String {
        let id = if chr.starts_with("20") { 2 } else { 1 };
        format!(
            "{chr}\tensembl_havana\t{feature_type}\t{start}\t{end}\t.\t+\t.\tgene_id \"ENSG0000000000{id}\"; transcript_id \"ENST0000000000{id}\";\n"
        )
    }

    /// A start or end no number leads, or an end below zero, refuses the
    /// parse on a line of a type the model skips, on a file's first sequence
    /// and on a later one, naming the file, the line, the column and the
    /// value, the rule the GFF3 reader follows
    /// ([`check_skipped_coordinates`]).
    #[test]
    fn a_coordinate_the_index_build_fails_is_refused_on_a_skipped_line_of_any_sequence() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        for (start, end, message) in [
            ("abc", "1099", "GTF start \"abc\" is not a number"),
            ("1000", "xyz", "GTF end \"xyz\" is not a number"),
            ("", "xyz", "GTF start \"\" is not a number"),
            ("1000", "", "GTF end \"\" is not a number"),
            ("1000", "-9", "GTF end \"-9\" is not a number"),
        ] {
            let on_first = line_of("20", "five_prime_utr", start, end);
            let err =
                try_parse_model_as(&two_sequences(&on_first, &sound_later), Reading::GffSource)
                    .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 9: {message}")),
                "{err:#}"
            );
            let on_later = line_of("21", "five_prime_utr", start, end);
            let err =
                try_parse_model_as(&two_sequences(&sound_first, &on_later), Reading::GffSource)
                    .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 17: {message}")),
                "{err:#}"
            );
        }
        // The columns are the line's as written: one that opens with a tab
        // holds its type where the index build reads a start.
        let shifted = format!("\t{sound_first}");
        let err = try_parse_model_as(&two_sequences(&shifted, &sound_later), Reading::GffSource)
            .unwrap_err();
        assert!(
            format!("{err:#}")
                .ends_with("t.gtf line 9: GTF start \"five_prime_utr\" is not a number"),
            "{err:#}"
        );
    }

    /// Nothing else of a skipped line is read, as Ensembl VEP's source reads a
    /// record's raw start, for its loop's bound, and its type, and nothing
    /// else of one it skips (`BaseGXF.pm:217-245`): a `five_prime_utr` line
    /// whose seqid could name no cache directory, whose start is 0, whose end
    /// precedes its start, or whose start or end a number leads without being
    /// one, an empty start before its end among them, the end not below
    /// zero, refuses nothing on either sequence, where the same column on an
    /// exon line refuses the parse, naming the file and the line.
    #[test]
    fn a_skipped_lines_seqid_and_number_led_coordinates_are_not_read() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        for (chr, start, end) in [
            ("20/x", "1000", "1099"),
            ("20", "0", "1099"),
            ("20", "1099", "1000"),
            ("20", "1000x", "1099"),
            ("20", "-5", "1099"),
            ("20", " 1000", "1099"),
            ("20", "0x3E8", "1099"),
            ("20", "1000", "1099x"),
            ("20", "", "1099"),
            ("20", "1000", "-0"),
        ] {
            for (first, later) in [
                (
                    line_of(chr, "five_prime_utr", start, end),
                    sound_later.clone(),
                ),
                (
                    sound_first.clone(),
                    line_of(&chr.replace("20", "21"), "five_prime_utr", start, end),
                ),
            ] {
                let model = try_parse_model_as(&two_sequences(&first, &later), Reading::GffSource)
                    .unwrap_or_else(|err| panic!("{chr:?} {start:?} {end:?}: {err:#}"));
                for gene in ["ENSG00000000001", "ENSG00000000002"] {
                    assert_eq!(model.genes[gene].transcripts[0].exons.len(), 2);
                }
            }
        }
        // An empty end is read on into a score that is a number.
        let scored = "20\tensembl_havana\tfive_prime_utr\t1000\t\t7\t+\t.\tgene_id \"ENSG00000000002\"; transcript_id \"ENST00000000002\";\n";
        assert!(
            try_parse_model_as(&two_sequences(scored, &sound_later), Reading::GffSource).is_ok()
        );
        for (chr, start, end, message) in [
            ("21/x", "1000", "1099", "invalid chromosome name '21/x'"),
            ("21", "abc", "1099", "GTF start \"abc\" is not a number"),
            ("21", "", "1099", "GTF start \"\" is not a number"),
            ("21", "1000x", "1099", "GTF start \"1000x\" is not a number"),
            (
                "21",
                "0",
                "1099",
                "GTF coordinates 0-1099 are not 1-based and ordered",
            ),
            (
                "21",
                "1099",
                "1000",
                "GTF coordinates 1099-1000 are not 1-based and ordered",
            ),
        ] {
            let exon = line_of(chr, "exon", start, end);
            let err = try_parse_model_as(&two_sequences(&sound_first, &exon), Reading::GffSource)
                .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 17: {message}")),
                "{err:#}"
            );
        }
    }

    /// A line of four to eight columns is skipped whatever its type and is
    /// read no further than its start and its end: one the index build fails,
    /// a blank one ending the line included, refuses the parse on either
    /// sequence, as on a line of nine columns, and any other spelling
    /// refuses nothing. The columns are the line's as written, so a line of
    /// a skipped type whose seqid is empty is read at its own start and end.
    #[test]
    fn a_short_line_is_refused_for_a_coordinate_the_index_build_fails_alone() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        for (columns, message) in [
            (
                "\tx\texon\tabc\t1099\t.\t+\t.\n",
                "GTF start \"abc\" is not a number",
            ),
            (
                "\tx\tfive_prime_utr\t1000\txyz\n",
                "GTF end \"xyz\" is not a number",
            ),
            (
                "\tx\tfive_prime_utr\t-\n",
                "GTF start \"-\" is not a number",
            ),
            ("\tx\tfive_prime_utr\t\n", "GTF start \"\" is not a number"),
            (
                "\tx\tfive_prime_utr\t1000\t \t\n",
                "GTF end \" \" is not a number",
            ),
            (
                "\tx\tfive_prime_utr\t1000\t-9\n",
                "GTF end \"-9\" is not a number",
            ),
        ] {
            let on_first = format!("20{columns}");
            let err =
                try_parse_model_as(&two_sequences(&on_first, &sound_later), Reading::GffSource)
                    .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 9: {message}")),
                "{err:#}"
            );
            let on_later = format!("21{columns}");
            let err =
                try_parse_model_as(&two_sequences(&sound_first, &on_later), Reading::GffSource)
                    .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 17: {message}")),
                "{err:#}"
            );
        }
        for columns in [
            "\tx\texon\t0\t1099\t.\t+\t.\n",
            "\tx\texon\t1099\t1000\t.\t+\t.\n",
            "\tx\texon\t1000x\t1099\n",
            "\tx\texon\t1000\n",
        ] {
            for (first, later) in [
                (format!("20{columns}"), sound_later.clone()),
                (sound_first.clone(), format!("21{columns}")),
            ] {
                let model = try_parse_model_as(&two_sequences(&first, &later), Reading::GffSource)
                    .unwrap_or_else(|err| panic!("{columns:?}: {err:#}"));
                for gene in ["ENSG00000000001", "ENSG00000000002"] {
                    assert_eq!(model.genes[gene].transcripts[0].exons.len(), 2);
                }
            }
        }
        // A line of a skipped type whose seqid is empty is read at its columns
        // as written.
        let unnamed = line_of("", "five_prime_utr", "1000", "1099");
        assert!(
            try_parse_model_as(&two_sequences(&unnamed, &sound_later), Reading::GffSource).is_ok()
        );
        let unnamed = line_of("", "five_prime_utr", "abc", "1099");
        let err = try_parse_model_as(&two_sequences(&unnamed, &sound_later), Reading::GffSource)
            .unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gtf line 9: GTF start \"abc\" is not a number"),
            "{err:#}"
        );
    }

    /// A line that is no comment and ends before its fourth column is a
    /// record without its start: the parse refuses it on either sequence,
    /// whatever its type, naming the column it ends at, as the GFF3 reader
    /// does. A GTF has no sequence section, so a line that begins with `>`
    /// is such a line; a blank line, empty or of white space alone, is
    /// skipped.
    #[test]
    fn a_line_that_ends_before_its_start_column_is_refused() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        for (line, message) in [
            (
                "20\tx\texon\n",
                "GTF line ends at column 3, before its start",
            ),
            ("20\tx\n", "GTF line ends at column 2, before its start"),
            (
                ">20 placeholder\n",
                "GTF line ends at column 1, before its start",
            ),
        ] {
            let err = try_parse_model_as(&two_sequences(line, &sound_later), Reading::GffSource)
                .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 9: {message}")),
                "{err:#}"
            );
            let err = try_parse_model_as(&two_sequences(&sound_first, line), Reading::GffSource)
                .unwrap_err();
            assert!(
                format!("{err:#}").ends_with(&format!("t.gtf line 17: {message}")),
                "{err:#}"
            );
        }
        for blank in ["\n", "  \n", "\t\t\t\n", "\r\n"] {
            let model = try_parse_model_as(&two_sequences(blank, &sound_later), Reading::GffSource)
                .unwrap_or_else(|err| panic!("{blank:?}: {err:#}"));
            for gene in ["ENSG00000000001", "ENSG00000000002"] {
                assert_eq!(model.genes[gene].transcripts[0].exons.len(), 2);
            }
        }
    }

    /// A byte-order mark at the start of a line changes nothing, as in the
    /// GFF3 reader ([`without_byte_order_marks`]): a model with a mark before
    /// each of its lines parses to the header and the genes of the same
    /// lines without the marks, the `#!` key behind a mark among them, and a
    /// marked record cut short refuses as it does unmarked.
    #[test]
    fn a_byte_order_mark_at_the_start_of_a_line_changes_nothing() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        let model = two_sequences(&sound_first, &sound_later);
        let marked: String = model.lines().map(|l| format!("\u{feff}{l}\n")).collect();
        let plain = parse_model(&model);
        let marked = parse_model(&marked);
        assert_eq!(marked.header, plain.header);
        assert_eq!(marked.header["genome-build"], "GRCh38.p14");
        assert_eq!(marked.genes.len(), 2);
        for gene in ["ENSG00000000001", "ENSG00000000002"] {
            assert_eq!(
                format!("{:?}", marked.genes[gene]),
                format!("{:?}", plain.genes[gene])
            );
        }
        let cut_short = format!("{model}\u{feff}20\tx\n");
        let err = try_parse_model_as(&cut_short, Reading::GffSource).unwrap_err();
        assert!(
            format!("{err:#}")
                .ends_with("t.gtf line 18: GTF line ends at column 2, before its start"),
            "{err:#}"
        );
    }

    /// The tags reader takes the mark off as the model reader does: a file
    /// with a mark before each of its lines gives the header and the tags of
    /// the same lines without the marks.
    #[test]
    fn the_tags_reader_reads_a_marked_file_as_the_file_without_its_marks() {
        let dir = tempfile::tempdir().unwrap();
        let marked: String = SEVEN_LINES
            .lines()
            .map(|l| format!("\u{feff}{l}\n"))
            .collect();
        let read = |name: &str, text: &str| {
            let path = dir.path().join(name);
            std::fs::write(&path, text).unwrap();
            parse_gtf_tags(&path).unwrap()
        };
        let plain = read("plain.gtf", SEVEN_LINES);
        let marked = read("marked.gtf", &marked);
        assert_eq!(marked.header, plain.header);
        assert_eq!(marked.header["genome-build"], "GRCh38.p14");
        assert_eq!(marked.tags, plain.tags);
        assert!(!plain.tags["ENST00000000001.5"].is_empty());
    }

    /// A line that opens with a tab is read at the columns the file writes,
    /// as the GFF3 reader reads one: a stray tab before a record puts its
    /// type in the start column, so the line refuses, the record's type kept
    /// or skipped, and a record whose seqid alone is empty refuses on a kept
    /// type. A line led by spaces alone, no tab before its seqid, is read
    /// trimmed.
    #[test]
    fn a_line_that_opens_with_a_tab_is_read_at_its_columns_as_written() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        let shifted = format!("\t{}", line_of("20", "exon", "1000", "1099"));
        let err = try_parse_model_as(&two_sequences(&shifted, &sound_later), Reading::GffSource)
            .unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gtf line 9: GTF start \"exon\" is not a number"),
            "{err:#}"
        );
        let shifted = format!("\t{}", line_of("21", "exon", "1000", "1099"));
        let err = try_parse_model_as(&two_sequences(&sound_first, &shifted), Reading::GffSource)
            .unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gtf line 17: GTF start \"exon\" is not a number"),
            "{err:#}"
        );
        let unnamed = line_of("", "exon", "1000", "1099");
        let err = try_parse_model_as(&two_sequences(&unnamed, &sound_later), Reading::GffSource)
            .unwrap_err();
        assert!(
            format!("{err:#}").ends_with("t.gtf line 9: invalid chromosome name ''"),
            "{err:#}"
        );
        let spaced = format!("  {sound_first}");
        assert!(
            try_parse_model_as(&two_sequences(&spaced, &sound_later), Reading::GffSource).is_ok()
        );
    }

    /// White space and then a tab before a record put its type where the
    /// index build reads a start, as in the GFF3 reader: the parse refuses
    /// the line on a kept type as on a skipped one, on either sequence.
    #[test]
    fn a_kept_line_led_by_white_space_and_a_tab_is_refused_as_its_skipped_twin_is() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        for lead in [" \t", "   \t"] {
            for feature_type in ["exon", "five_prime_utr"] {
                let message = format!("GTF start \"{feature_type}\" is not a number");
                let on_first = format!("{lead}{}", line_of("20", feature_type, "1000", "1099"));
                let text = two_sequences(&on_first, &sound_later);
                let err = try_parse_model_as(&text, Reading::GffSource).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gtf line 9: {message}")),
                    "{err:#}"
                );
                let on_later = format!("{lead}{}", line_of("21", feature_type, "1000", "1099"));
                let text = two_sequences(&sound_first, &on_later);
                let err = try_parse_model_as(&text, Reading::GffSource).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gtf line 17: {message}")),
                    "{err:#}"
                );
            }
        }
    }

    /// A line led by spaces alone, no tab before its seqid, is read trimmed:
    /// a gene line's seqid is the one behind the spaces, a kept record joins
    /// its transcript and a skipped one is skipped.
    #[test]
    fn a_line_led_by_spaces_alone_is_read_trimmed() {
        let sound_first = line_of("20", "five_prime_utr", "1000", "1099");
        let sound_later = line_of("21", "five_prime_utr", "1000", "1099");
        let model = two_sequences(&sound_first, &sound_later).replacen("\n20\t", "\n  20\t", 1);
        let parsed = parse_model_as(&model, Reading::GffSource);
        assert_eq!(parsed.genes["ENSG00000000002"].chr, "20");
        let exon = format!("  {}", line_of("20", "exon", "1300", "1400"));
        let parsed = parse_model_as(&two_sequences(&exon, &sound_later), Reading::GffSource);
        assert_eq!(
            parsed.genes["ENSG00000000002"].transcripts[0].exons.len(),
            3
        );
        let utr = format!("  {sound_first}");
        let parsed = parse_model_as(&two_sequences(&utr, &sound_later), Reading::GffSource);
        assert_eq!(
            parsed.genes["ENSG00000000002"].transcripts[0].exons.len(),
            2
        );
    }

    /// Both readers refuse an attribute value whose escapes spell bytes that
    /// are not UTF-8 with one text ([`undecodable_attribute_error`]), the file
    /// and the line before it, so the same offending pair (`Name` holding
    /// `caf%FF`) draws the same words from a GFF3 exon line and a GTF exon line
    /// but for the format's name.
    #[test]
    fn both_readers_refuse_an_undecodable_value_with_one_text() {
        let dir = tempfile::tempdir().unwrap();
        let gff3 = dir.path().join("t.gff3");
        std::fs::write(
            &gff3,
            "21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T1;Name=caf%FF\n",
        )
        .unwrap();
        let gtf = dir.path().join("t.gtf");
        std::fs::write(
            &gtf,
            format!(
                "{SEVEN_LINES}21\tensembl_havana\texon\t1000\t1099\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; Name \"caf%FF\";\n"
            ),
        )
        .unwrap();
        let gff3_text = format!(
            "{:#}",
            crate::gff3::parse_gff3(&gff3, Reading::GffSource).unwrap_err()
        );
        let gtf_text = format!(
            "{:#}",
            parse_gtf_model(&gtf, Reading::GffSource).unwrap_err()
        );
        let gff3_refusal = undecodable_attribute_error("GFF3", "Name", "caf%FF");
        let gtf_refusal = undecodable_attribute_error("GTF", "Name", "caf%FF");
        assert_eq!(
            gff3_text,
            format!("{} line 3: {gff3_refusal}", gff3.display())
        );
        assert_eq!(gtf_text, format!("{} line 9: {gtf_refusal}", gtf.display()));
        let after_format = |text: &str| {
            text.split_once(" attribute ")
                .map(|(_, rest)| rest.to_string())
        };
        assert_eq!(after_format(&gff3_text), after_format(&gtf_text));
        assert_eq!(
            gff3_refusal,
            "GFF3 attribute Name \"caf%FF\" decodes to text that is not UTF-8"
        );
    }

    /// A line that is not UTF-8 fails both readers, naming the file and the
    /// line, as every refusal of a line does.
    #[test]
    fn a_line_that_is_not_utf8_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.gtf");
        let (before_name, after_name) = SEVEN_LINES.split_once("GENE1").unwrap();
        let mut bytes = before_name.as_bytes().to_vec();
        bytes.extend_from_slice(b"GEN\xc9");
        bytes.extend_from_slice(after_name.as_bytes());
        std::fs::write(&path, bytes).unwrap();
        for err in [
            parse_gtf_model(&path, Reading::GffSource).unwrap_err(),
            parse_gtf_tags(&path).unwrap_err(),
        ] {
            assert!(
                format!("{err:#}").ends_with("t.gtf line 2: GTF line is not UTF-8"),
                "{err:#}"
            );
        }
    }

    /// A kept line that would draw two refusals draws the one read first, the
    /// start and end columns before the attribute column and the attribute
    /// column before the version, under either reading: an `exon` line whose
    /// `exon_id` decodes to bytes that are not UTF-8 is refused for a start
    /// or an end that is no number and for an end before its start, and a
    /// `CDS` line that repeats `protein_id` around a `version` that is no
    /// whole number is refused for the repeated key.
    #[test]
    fn a_kept_line_draws_the_refusal_of_the_column_read_first() {
        let exon = |start: &str, end: &str| {
            format!(
                "21\tensembl_havana\texon\t{start}\t{end}\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_id \"caf%FF\";\n"
            )
        };
        let cds = "21\tensembl_havana\tCDS\t1100\t1200\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; protein_id \"ENSP00000000001\"; version \"v2\"; protein_id \"ENSP00000000001\";\n";
        for (line, refusal) in [
            (exon("abc", "1099"), "GTF start \"abc\" is not a number"),
            (exon("1000", "xyz"), "GTF end \"xyz\" is not a number"),
            (
                exon("1200", "1100"),
                "GTF coordinates 1200-1100 are not 1-based and ordered",
            ),
            (
                cds.to_string(),
                "GTF attribute protein_id is repeated; Ensembl VEP reads a repeated key as a list",
            ),
        ] {
            for reading in [Reading::GffSource, Reading::EnsemblCache] {
                let err = try_parse_model_as(&format!("{SEVEN_LINES}{line}"), reading).unwrap_err();
                assert!(
                    format!("{err:#}").ends_with(&format!("t.gtf line 9: {refusal}")),
                    "{reading:?}: {err:#}"
                );
            }
        }
    }

    /// A kept line that would draw a refusal for its seqid and one for its
    /// start draws the seqid's, the column read first, under either reading:
    /// an `exon` line whose seqid can name no cache directory and whose start
    /// is no number is refused for the seqid.
    #[test]
    fn a_kept_line_is_refused_for_its_seqid_before_its_start() {
        let exon = line_of("21/x", "exon", "abc", "1099");
        for reading in [Reading::GffSource, Reading::EnsemblCache] {
            let err = try_parse_model_as(&format!("{SEVEN_LINES}{exon}"), reading).unwrap_err();
            assert!(
                format!("{err:#}").ends_with("t.gtf line 9: invalid chromosome name '21/x'"),
                "{reading:?}: {err:#}"
            );
        }
    }
}
