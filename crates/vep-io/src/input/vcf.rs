// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF input parser using noodles.

use std::io::BufRead;

use noodles::vcf;

use vep_core::variant::{InputVariant, VariantClass};

use crate::error::IoError;

use super::InputParser;

/// VCF input parser that reads VCF records and converts them to `InputVariant`s.
///
/// For multi-allelic sites, each ALT allele produces a separate `InputVariant`.
/// VCF coordinates are converted to VEP-style coordinates (trimming shared prefixes).
pub struct VcfParser {
    reader: vcf::io::Reader<Box<dyn BufRead + Send>>,
    header: vcf::Header,
    /// Buffer of variants from the current record (for multi-allelic sites).
    buffer: Vec<InputVariant>,
}

impl VcfParser {
    /// Create a new VCF parser from a buffered reader.
    pub fn new(reader: Box<dyn BufRead + Send>) -> Result<Self, IoError> {
        let mut vcf_reader = vcf::io::Reader::new(reader);
        let header = vcf_reader
            .read_header()
            .map_err(|e| IoError::VcfParse(format!("failed to read VCF header: {e}")))?;
        Ok(Self {
            reader: vcf_reader,
            header,
            buffer: Vec::new(),
        })
    }

    /// Create a new VCF parser from a file path (single-threaded decompression).
    pub fn from_path(path: &std::path::Path) -> Result<Self, IoError> {
        Self::from_path_with_decompression_threads(path, 1)
    }

    /// Create a new VCF parser from a file path with configurable BGZF worker count.
    /// When `decompression_threads > 1` and the file is gzipped, uses
    /// `noodles::bgzf::MultithreadedReader` for parallel block decompression.
    pub fn from_path_with_decompression_threads(
        path: &std::path::Path,
        decompression_threads: usize,
    ) -> Result<Self, IoError> {
        let file = std::fs::File::open(path)?;
        let reader: Box<dyn BufRead + Send> = if path
            .extension()
            .is_some_and(|ext| ext == "gz" || ext == "bgz")
        {
            match std::num::NonZeroUsize::new(decompression_threads) {
                Some(n) if n.get() > 1 => Box::new(std::io::BufReader::new(
                    noodles::bgzf::io::MultithreadedReader::with_worker_count(n, file),
                )),
                _ => Box::new(std::io::BufReader::new(noodles::bgzf::io::Reader::new(
                    file,
                ))),
            }
        } else {
            Box::new(std::io::BufReader::new(file))
        };
        Self::new(reader)
    }

    /// Parse a single VCF record into one or more InputVariants (one per ALT allele).
    fn parse_record(&self, record: &vcf::Record) -> Result<Vec<InputVariant>, IoError> {
        let chrom = record.reference_sequence_name().to_string();

        let pos: u64 = record
            .variant_start()
            .ok_or_else(|| IoError::VcfParse("missing POS".into()))?
            .map_err(|e| IoError::VcfParse(format!("invalid POS: {e}")))?
            .get() as u64;

        let vcf_ref = record.reference_bases();

        let ids = record.ids();
        let ids_str: &str = ids.as_ref();
        let id = if ids_str.is_empty() {
            None
        } else {
            Some(ids_str.to_string())
        };

        let alt_bases = record.alternate_bases();
        use noodles::vcf::variant::record::AlternateBases;
        let alts: Vec<Vec<u8>> = alt_bases
            .iter()
            .filter_map(|r| r.ok())
            .filter(|alt| *alt != "." && *alt != "*" && *alt != "<*>")
            .map(|s| s.as_bytes().to_vec())
            .collect();

        if alts.is_empty() {
            return Ok(Vec::new());
        }

        let mut symbolic_indices = Vec::new();
        let mut explicit_indices = Vec::new();
        for (i, alt) in alts.iter().enumerate() {
            if is_symbolic_allele(alt) {
                symbolic_indices.push(i);
            } else {
                explicit_indices.push(i);
            }
        }

        let multi = alts.len() > 1;
        let mut variants = Vec::with_capacity(alts.len());

        if !explicit_indices.is_empty() {
            let explicit_alts: Vec<Vec<u8>> =
                explicit_indices.iter().map(|&i| alts[i].clone()).collect();
            let vcf_ref_bytes = vcf_ref.as_bytes();
            let (trimmed_ref, trimmed_alts, prefix_trimmed) =
                trim_multi_alleles(vcf_ref_bytes, &explicit_alts);

            for (local_idx, trimmed_alt) in trimmed_alts.into_iter().enumerate() {
                let allele_idx = explicit_indices[local_idx];
                let (start, end) = compute_vep_coords(
                    pos,
                    prefix_trimmed,
                    &trimmed_ref,
                    &trimmed_alt,
                    vcf_ref_bytes.len(),
                );

                let mut variant =
                    InputVariant::new(chrom.clone(), start, end, trimmed_ref.clone(), trimmed_alt);
                variant.id = id.clone();
                variant.allele_index = allele_idx;
                variant.minimised = multi;
                variants.push(variant);
            }
        }

        if !symbolic_indices.is_empty() {
            let info_fields = parse_info_fields(record, &self.header);

            // In a mixed record (explicit and symbolic ALTs) Perl VEP joins the
            // ALTs ("G/<INV>"), which matches no SV pattern, so it classifies by
            // SVTYPE instead of the symbolic allele.
            let is_mixed = !explicit_indices.is_empty() && !symbolic_indices.is_empty();

            for &allele_idx in &symbolic_indices {
                let alt = &alts[allele_idx];
                let alt_str = std::str::from_utf8(alt).unwrap_or("");

                let sv_class = if is_mixed {
                    if let Some(ref svtype) = info_fields.svtype {
                        classify_symbolic_allele(&format!("<{}>", svtype), Some(svtype.as_str()))
                    } else {
                        classify_symbolic_allele(alt_str, info_fields.svtype.as_deref())
                    }
                } else {
                    classify_symbolic_allele(alt_str, info_fields.svtype.as_deref())
                };

                // End coordinate precedence: INFO/END, then POS + abs(SVLEN), then POS.
                let sv_end = info_fields.end.unwrap_or_else(|| {
                    info_fields
                        .svlen
                        .map(|l| pos + l.unsigned_abs())
                        .unwrap_or(pos)
                });

                let ref_bytes = vcf_ref.as_bytes().to_vec();

                let mut variant =
                    InputVariant::new(chrom.clone(), pos, sv_end, ref_bytes, alt.clone());
                variant.variant_class = sv_class;
                variant.id = id.clone();
                variant.allele_index = allele_idx;
                variant.minimised = multi;
                variant.is_structural = true;
                variant.sv_end = Some(sv_end);
                variant.sv_type = info_fields.svtype.clone();
                variant.sv_len = info_fields.svlen;
                variant.mate_id = info_fields.mateid.clone();
                variants.push(variant);
            }
        }

        variants.sort_by_key(|v| v.allele_index);

        Ok(variants)
    }
}

/// Trim shared prefix/suffix across REF and all ALT alleles together.
///
/// This mirrors Perl VEP multi-allelic handling where all alleles at a site
/// share a single normalized coordinate window.
fn trim_multi_alleles(
    ref_allele: &[u8],
    alt_alleles: &[Vec<u8>],
) -> (Vec<u8>, Vec<Vec<u8>>, usize) {
    let mut prefix_len = 0usize;
    loop {
        if prefix_len >= ref_allele.len() {
            break;
        }
        let base = ref_allele[prefix_len].to_ascii_uppercase();
        if alt_alleles
            .iter()
            .any(|a| prefix_len >= a.len() || a[prefix_len].to_ascii_uppercase() != base)
        {
            break;
        }
        prefix_len += 1;
    }

    let ref_after_prefix = &ref_allele[prefix_len..];
    let alts_after_prefix: Vec<&[u8]> = alt_alleles.iter().map(|a| &a[prefix_len..]).collect();

    // Perl VEP trims only a shared prefix, never a shared suffix; suffix
    // trimming shifts MNP and complex coordinates and alleles away from Perl's.
    let suffix_len = 0usize;

    let ref_core_end = ref_after_prefix.len().saturating_sub(suffix_len);
    let ref_core = &ref_after_prefix[..ref_core_end];
    let trimmed_ref = if ref_core.is_empty() {
        b"-".to_vec()
    } else {
        ref_core.to_vec()
    };

    let mut trimmed_alts = Vec::with_capacity(alts_after_prefix.len());
    for alt in alts_after_prefix {
        let alt_core_end = alt.len().saturating_sub(suffix_len);
        let alt_core = &alt[..alt_core_end];
        trimmed_alts.push(if alt_core.is_empty() {
            b"-".to_vec()
        } else {
            alt_core.to_vec()
        });
    }

    (trimmed_ref, trimmed_alts, prefix_len)
}

/// Compute VEP-style start and end coordinates from VCF position after allele trimming.
fn compute_vep_coords(
    vcf_pos: u64,
    prefix_trimmed: usize,
    trimmed_ref: &[u8],
    trimmed_alt: &[u8],
    _original_ref_len: usize,
) -> (u64, u64) {
    let start = vcf_pos + prefix_trimmed as u64;

    let ref_is_dash = trimmed_ref == b"-";
    let _alt_is_dash = trimmed_alt == b"-";

    if ref_is_dash {
        // VEP convention: an insertion at position N has start=N, end=N-1.
        (start, start - 1)
    } else {
        let ref_len = trimmed_ref.len() as u64;
        (start, start + ref_len - 1)
    }
}

/// Returns true if the ALT allele is a symbolic (structural variant) allele.
///
/// Symbolic alleles include:
/// - Angle-bracket notation: `<DEL>`, `<INS>`, `<DUP:TANDEM>`, etc.
/// - BND notation: `A[chr:pos[`, `]chr:pos]A`, etc.
/// - Single breakend: `A.` or `.A` (length > 1)
fn is_symbolic_allele(alt: &[u8]) -> bool {
    if alt.is_empty() {
        return false;
    }
    if alt.first() == Some(&b'<') && alt.last() == Some(&b'>') {
        return true;
    }
    if alt.contains(&b'[') || alt.contains(&b']') {
        return true;
    }
    if alt.len() > 1 && (alt.last() == Some(&b'.') || alt.first() == Some(&b'.')) {
        return true;
    }
    false
}

/// Parsed INFO fields relevant to structural variants.
struct SvInfoFields {
    end: Option<u64>,
    svtype: Option<String>,
    svlen: Option<i64>,
    mateid: Option<String>,
}

/// Parse structural-variant-relevant INFO fields from the record.
///
/// Uses the noodles header-aware API when possible, with a fallback to raw
/// string parsing for VCFs that lack INFO header definitions (common in
/// minimal or synthetic test files).
fn parse_info_fields(record: &vcf::Record, header: &vcf::Header) -> SvInfoFields {
    let info = record.info();
    let info_raw: &str = info.as_ref();

    let end = get_info_integer(record, header, "END")
        .map(|v| v as u64)
        .or_else(|| get_raw_info_value(info_raw, "END").and_then(|s| s.parse::<u64>().ok()));

    let svtype = get_info_string(record, header, "SVTYPE")
        .or_else(|| get_raw_info_value(info_raw, "SVTYPE").map(|s| s.to_string()));

    let svlen = get_info_integer(record, header, "SVLEN")
        .map(|v| v as i64)
        .or_else(|| get_raw_info_value(info_raw, "SVLEN").and_then(|s| s.parse::<i64>().ok()));

    let mateid = get_info_string(record, header, "MATEID")
        .or_else(|| get_raw_info_value(info_raw, "MATEID").map(|s| s.to_string()));

    SvInfoFields {
        end,
        svtype,
        svlen,
        mateid,
    }
}

/// Extract an integer value from a VCF INFO field via the noodles typed API.
fn get_info_integer(record: &vcf::Record, header: &vcf::Header, key: &str) -> Option<i32> {
    use noodles::vcf::variant::record::info::field::Value;
    record.info().get(header, key)?.ok()?.and_then(|v| match v {
        Value::Integer(n) => Some(n),
        _ => None,
    })
}

/// Extract a string value from a VCF INFO field via the noodles typed API.
fn get_info_string(record: &vcf::Record, header: &vcf::Header, key: &str) -> Option<String> {
    use noodles::vcf::variant::record::info::field::Value;
    record.info().get(header, key)?.ok()?.and_then(|v| match v {
        Value::String(s) => Some(s.to_string()),
        _ => None,
    })
}

/// Fallback: parse a single INFO field value from the raw semicolon-delimited
/// INFO string. Returns the value portion of `KEY=VALUE`, or `None` if the
/// key is not found.
fn get_raw_info_value<'a>(info_raw: &'a str, key: &str) -> Option<&'a str> {
    if info_raw == "." || info_raw.is_empty() {
        return None;
    }
    for field in info_raw.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key {
                return Some(v);
            }
        }
    }
    None
}

/// Classify a symbolic ALT allele string into a `VariantClass`.
///
/// Matches the Perl VEP `get_SO_term()` logic for mapping SV type abbreviations
/// to sequence ontology terms.
fn classify_symbolic_allele(alt: &str, svtype: Option<&str>) -> VariantClass {
    let upper = alt.to_ascii_uppercase();

    if upper.contains(":ME") {
        if upper.starts_with("<DEL") {
            return VariantClass::MobileElementDeletion;
        }
        return VariantClass::MobileElementInsertion;
    }

    if upper.contains("DUP:TANDEM") {
        return VariantClass::TandemDuplication;
    }

    if upper.contains("CNV:TR") {
        return VariantClass::TandemRepeat;
    }

    // Perl VEP maps CN0 -> DEL (deletion), CN2 -> DUP (duplication),
    // all others -> CNV (copy_number_variation).
    if upper.starts_with("<CN") && upper.ends_with('>') {
        let inner = &upper[1..upper.len() - 1];
        if inner == "CNV" {
            return VariantClass::CopyNumberVariation;
        }
        let after_cn = &inner[2..]; // strip "CN"
        let num_str = after_cn.strip_prefix('=').unwrap_or(after_cn);
        if let Ok(cn) = num_str.parse::<u32>() {
            return match cn {
                0 => VariantClass::StructuralDeletion,
                2 => VariantClass::Duplication,
                _ => VariantClass::CopyNumberVariation,
            };
        }
        if inner.starts_with("CN") {
            return VariantClass::CopyNumberVariation;
        }
    }

    if alt.contains('[') || alt.contains(']') {
        return VariantClass::Translocation;
    }

    if alt.len() > 1 && (alt.ends_with('.') || alt.starts_with('.')) {
        return VariantClass::Translocation;
    }

    if upper.contains("CPX") {
        return VariantClass::ComplexStructural;
    }

    let stripped = upper
        .trim_start_matches('<')
        .trim_end_matches('>')
        .split(':')
        .next()
        .unwrap_or("");

    // Perl VEP falls back to INFO/SVTYPE when the ALT matches no standard VCF
    // 4.4 pattern.
    let effective_type = if stripped.is_empty() || stripped == "*" {
        svtype.map(|s| s.to_ascii_uppercase()).unwrap_or_default()
    } else {
        stripped.to_string()
    };

    match effective_type.as_str() {
        "DEL" => VariantClass::StructuralDeletion,
        "INS" => VariantClass::StructuralInsertion,
        "DUP" => VariantClass::Duplication,
        "INV" => VariantClass::Inversion,
        "CNV" => VariantClass::CopyNumberVariation,
        "BND" => VariantClass::Translocation,
        _ => {
            if let Some(sv) = svtype {
                let sv_upper = sv.to_ascii_uppercase();
                match sv_upper.as_str() {
                    "DEL" => return VariantClass::StructuralDeletion,
                    "INS" => return VariantClass::StructuralInsertion,
                    "DUP" => return VariantClass::Duplication,
                    "INV" => return VariantClass::Inversion,
                    "CNV" => return VariantClass::CopyNumberVariation,
                    "BND" => return VariantClass::Translocation,
                    _ => {}
                }
            }
            VariantClass::ComplexStructural
        }
    }
}

/// Parse BND notation to extract the remote chromosome and position.
///
/// BND patterns (per VCF 4.3+ spec):
/// - `N[chr:pos[`: forward to forward
/// - `N]chr:pos]`: forward to reverse
/// - `]chr:pos]N`: reverse to forward
/// - `[chr:pos[N`: reverse to reverse
///
/// Returns `(remote_chr, remote_pos)` if parsing succeeds.
pub fn parse_bnd_alt(alt: &str) -> Option<(String, u64)> {
    if let Some(bracket_start) = alt.find(['[', ']']) {
        let bracket_char = alt.as_bytes()[bracket_start] as char;
        let after_bracket = &alt[bracket_start + 1..];
        if let Some(close_pos) = after_bracket.find(bracket_char) {
            let coords = &after_bracket[..close_pos];
            if let Some((chr, pos_str)) = coords.split_once(':') {
                if let Ok(pos) = pos_str.parse::<u64>() {
                    return Some((chr.to_string(), pos));
                }
            }
        }
    }
    None
}

impl InputParser for VcfParser {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        if let Some(variant) = self.buffer.pop() {
            return Some(Ok(variant));
        }

        let mut record = vcf::Record::default();
        loop {
            match self.reader.read_record(&mut record) {
                Ok(0) => return None,
                Ok(_) => {
                    match self.parse_record(&record) {
                        Ok(variants) => {
                            if variants.is_empty() {
                                continue;
                            }
                            // Reversed so `pop` yields alleles in input order.
                            let mut variants = variants;
                            variants.reverse();
                            let first = variants.pop().unwrap();
                            self.buffer = variants;
                            return Some(Ok(first));
                        }
                        Err(e) => return Some(Err(e)),
                    }
                }
                Err(e) => {
                    return Some(Err(IoError::VcfParse(format!(
                        "error reading VCF record: {e}"
                    ))));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_vcf_data(records: &str) -> Vec<u8> {
        format!(
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             {records}"
        )
        .into_bytes()
    }

    /// Build a VCF carrying structural-variant records (END, SVTYPE, SVLEN,
    /// MATEID in INFO) under the same minimal header as `make_vcf_data`.
    ///
    /// noodles needs no INFO header definitions to parse the records; the raw
    /// INFO string fallback extracts the values.
    fn make_sv_vcf_data(records: &str) -> Vec<u8> {
        make_vcf_data(records)
    }

    #[test]
    fn test_parse_snv() {
        let data = make_vcf_data("21\t25000100\trs123\tA\tG\t.\tPASS\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let variant = parser.next_variant().unwrap().unwrap();
        assert_eq!(variant.chr, "21");
        assert_eq!(variant.start, 25000100);
        assert_eq!(variant.end, 25000100);
        assert_eq!(variant.ref_allele, b"A");
        assert_eq!(variant.alt_allele(), b"G");
        assert_eq!(variant.allele_string, "A/G");
        assert_eq!(variant.id, Some("rs123".into()));
        assert_eq!(variant.uploaded_variation(), "rs123");

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_insertion() {
        // VCF insertion: REF=A, ALT=ATT at POS=100
        // VEP: ref="-", alt="TT", start=101, end=100
        let data = make_vcf_data("21\t100\t.\tA\tATT\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let variant = parser.next_variant().unwrap().unwrap();
        assert_eq!(variant.chr, "21");
        assert_eq!(variant.start, 101);
        assert_eq!(variant.end, 100);
        assert_eq!(variant.ref_allele, b"-");
        assert_eq!(variant.alt_allele(), b"TT");
        assert_eq!(variant.allele_string, "-/TT");
        assert_eq!(variant.id, None);
    }

    #[test]
    fn test_parse_deletion() {
        // VCF deletion: REF=ATT, ALT=A at POS=100
        // VEP: ref="TT", alt="-", start=101, end=102
        let data = make_vcf_data("21\t100\t.\tATT\tA\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let variant = parser.next_variant().unwrap().unwrap();
        assert_eq!(variant.chr, "21");
        assert_eq!(variant.start, 101);
        assert_eq!(variant.end, 102);
        assert_eq!(variant.ref_allele, b"TT");
        assert_eq!(variant.alt_allele(), b"-");
        assert_eq!(variant.allele_string, "TT/-");
    }

    #[test]
    fn test_parse_multi_allelic() {
        let data = make_vcf_data("21\t200\t.\tA\tG,T\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.allele_string, "A/G");
        assert_eq!(v1.allele_index, 0);
        assert!(v1.minimised);

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.allele_string, "A/T");
        assert_eq!(v2.allele_index, 1);
        assert!(v2.minimised);

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_multi_allelic_indel_uses_shared_window() {
        // Perl-style normalization for CTTT/C/CTT is:
        // shared prefix C trimmed once => TTT/-/TT at the same locus.
        let data = make_vcf_data("1\t100\t.\tCTTT\tC,CTT\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.start, 101);
        assert_eq!(v1.end, 103);
        assert_eq!(v1.ref_allele, b"TTT");
        assert_eq!(v1.alt_allele(), b"-");

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.start, 101);
        assert_eq!(v2.end, 103);
        assert_eq!(v2.ref_allele, b"TTT");
        assert_eq!(v2.alt_allele(), b"TT");

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_missing_id() {
        let data = make_vcf_data("21\t300\t.\tC\tT\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let variant = parser.next_variant().unwrap().unwrap();
        assert_eq!(variant.id, None);
        assert_eq!(variant.uploaded_variation(), "21_300_C/T");
    }

    #[test]
    fn test_parse_monomorphic_site_skipped() {
        let data = make_vcf_data("21\t400\t.\tA\t.\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_multiple_records() {
        let data = make_vcf_data(
            "21\t100\trs1\tA\tG\t.\t.\t.\n\
             21\t200\trs2\tC\tT\t.\t.\t.\n",
        );
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.id, Some("rs1".into()));
        assert_eq!(v1.start, 100);

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.id, Some("rs2".into()));
        assert_eq!(v2.start, 200);

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_complex_indel() {
        // VCF complex: REF=ACGT, ALT=ATGT at POS=100
        // After prefix-only trim: ref=CGT, alt=TGT, trimmed 1 from start
        // VEP: start=101, end=103 (Perl VEP does not suffix-trim)
        let data = make_vcf_data("21\t100\t.\tACGT\tATGT\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let variant = parser.next_variant().unwrap().unwrap();
        assert_eq!(variant.start, 101);
        assert_eq!(variant.end, 103);
        assert_eq!(variant.ref_allele, b"CGT");
        assert_eq!(variant.alt_allele(), b"TGT");
    }

    #[test]
    fn test_is_symbolic_allele() {
        assert!(is_symbolic_allele(b"<DEL>"));
        assert!(is_symbolic_allele(b"<INS>"));
        assert!(is_symbolic_allele(b"<DUP:TANDEM>"));
        assert!(is_symbolic_allele(b"<CNV>"));
        assert!(is_symbolic_allele(b"<INS:ME:ALU>"));

        assert!(is_symbolic_allele(b"A[21:100["));
        assert!(is_symbolic_allele(b"]21:100]A"));
        assert!(is_symbolic_allele(b"A]21:100]"));
        assert!(is_symbolic_allele(b"[21:100[A"));

        assert!(is_symbolic_allele(b"A."));
        assert!(is_symbolic_allele(b".A"));

        assert!(!is_symbolic_allele(b"A"));
        assert!(!is_symbolic_allele(b"ACGT"));
        assert!(!is_symbolic_allele(b"G"));
        assert!(!is_symbolic_allele(b".")); // single dot is not symbolic (len=1)
        assert!(!is_symbolic_allele(b""));
    }

    #[test]
    fn test_parse_sv_del_with_end() {
        let data = make_sv_vcf_data("21\t1000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1100\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
        assert_eq!(v.start, 1000);
        assert_eq!(v.end, 1100);
        assert_eq!(v.sv_end, Some(1100));
        assert_eq!(v.sv_type, Some("DEL".into()));
        assert!(v.is_structural);
        // REF should be preserved (not trimmed)
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_allele(), b"<DEL>");
        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_sv_ins() {
        let data = make_sv_vcf_data("21\t500\t.\tA\t<INS>\t.\t.\tSVTYPE=INS\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::StructuralInsertion);
        assert!(v.is_structural);
        assert_eq!(v.start, 500);
        // No END or SVLEN -> end = POS
        assert_eq!(v.end, 500);
    }

    #[test]
    fn test_parse_sv_dup_tandem() {
        let data = make_sv_vcf_data("21\t200\t.\tA\t<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;END=500\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::TandemDuplication);
        assert!(v.is_structural);
        assert_eq!(v.sv_end, Some(500));
    }

    #[test]
    fn test_parse_sv_inv_with_svlen() {
        // <INV> with SVLEN=500 at POS=1000 -> end = 1000 + 500 = 1500
        let data = make_sv_vcf_data("21\t1000\t.\tA\t<INV>\t.\t.\tSVTYPE=INV;SVLEN=500\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::Inversion);
        assert!(v.is_structural);
        assert_eq!(v.sv_len, Some(500));
        // No END -> fallback to POS + abs(SVLEN)
        assert_eq!(v.end, 1500);
        assert_eq!(v.sv_end, Some(1500));
    }

    #[test]
    fn test_parse_sv_cnv() {
        let data = make_sv_vcf_data("21\t300\t.\tA\t<CNV>\t.\t.\tSVTYPE=CNV;END=800\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
        assert!(v.is_structural);
    }

    #[test]
    fn test_parse_sv_cn0() {
        // <CN0> -> StructuralDeletion (Perl VEP maps CN0 to DEL)
        let data = make_sv_vcf_data("21\t300\t.\tA\t<CN0>\t.\t.\tSVTYPE=CNV;END=800\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
        assert!(v.is_structural);
    }

    #[test]
    fn test_parse_sv_mobile_element_insertion() {
        let data = make_sv_vcf_data("21\t400\t.\tA\t<INS:ME:ALU>\t.\t.\tSVTYPE=INS\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::MobileElementInsertion);
        assert!(v.is_structural);
    }

    #[test]
    fn test_parse_sv_bnd() {
        let data =
            make_sv_vcf_data("13\t500\tbnd_1\tA\tA[21:100[\t.\t.\tSVTYPE=BND;MATEID=bnd_2\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::Translocation);
        assert!(v.is_structural);
        assert_eq!(v.mate_id, Some("bnd_2".into()));
        assert_eq!(v.sv_type, Some("BND".into()));
    }

    #[test]
    fn test_parse_mixed_multi_allelic_snv_and_del() {
        // Mixed: ALT = G,<DEL>  -> SNV (allele 0) + StructuralDeletion (allele 1)
        let data = make_sv_vcf_data("21\t1000\t.\tA\tG,<DEL>\t.\t.\tSVTYPE=DEL;END=1100\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.allele_index, 0);
        assert_eq!(v1.variant_class, VariantClass::Snv);
        assert!(!v1.is_structural);
        assert_eq!(v1.ref_allele, b"A");
        assert_eq!(v1.alt_allele(), b"G");
        assert!(v1.minimised);

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.allele_index, 1);
        assert_eq!(v2.variant_class, VariantClass::StructuralDeletion);
        assert!(v2.is_structural);
        assert_eq!(v2.sv_end, Some(1100));
        assert!(v2.minimised);

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_mixed_multi_allelic_svtype_overrides_symbolic() {
        // Mixed: ALT = G,<INV> with SVTYPE=DEL
        // Perl VEP uses SVTYPE for the class in mixed multi-allelic records,
        // producing "deletion" instead of "inversion" for the symbolic allele.
        let data =
            make_sv_vcf_data("21\t1000\t.\tC\tG,<INV>\t.\t.\tSVTYPE=DEL;END=2000;SVLEN=-1000\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.allele_index, 0);
        assert_eq!(v1.variant_class, VariantClass::Snv);
        assert!(!v1.is_structural);

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.allele_index, 1);
        assert_eq!(v2.variant_class, VariantClass::StructuralDeletion);
        assert!(v2.is_structural);
        assert_eq!(v2.sv_end, Some(2000));

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_symbolic_allele_not_trimmed() {
        // Symbolic alleles bypass prefix trimming: REF=A must survive.
        let data = make_sv_vcf_data("21\t100\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=200\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_allele(), b"<DEL>");
        assert_eq!(v.start, 100);
        assert_eq!(v.end, 200);
    }

    #[test]
    fn test_classify_symbolic_allele_patterns() {
        assert_eq!(
            classify_symbolic_allele("<DEL>", None),
            VariantClass::StructuralDeletion
        );
        assert_eq!(
            classify_symbolic_allele("<INS>", None),
            VariantClass::StructuralInsertion
        );
        assert_eq!(
            classify_symbolic_allele("<DUP>", None),
            VariantClass::Duplication
        );
        assert_eq!(
            classify_symbolic_allele("<DUP:TANDEM>", None),
            VariantClass::TandemDuplication
        );
        assert_eq!(
            classify_symbolic_allele("<INV>", None),
            VariantClass::Inversion
        );
        assert_eq!(
            classify_symbolic_allele("<CNV>", None),
            VariantClass::CopyNumberVariation
        );
        assert_eq!(
            classify_symbolic_allele("<CN0>", None),
            VariantClass::StructuralDeletion
        );
        assert_eq!(
            classify_symbolic_allele("<CN2>", None),
            VariantClass::Duplication
        );
        assert_eq!(
            classify_symbolic_allele("<CN3>", None),
            VariantClass::CopyNumberVariation
        );
        assert_eq!(
            classify_symbolic_allele("<INS:ME:ALU>", None),
            VariantClass::MobileElementInsertion
        );
        assert_eq!(
            classify_symbolic_allele("<INS:ME:LINE1>", None),
            VariantClass::MobileElementInsertion
        );
        assert_eq!(
            classify_symbolic_allele("<DEL:ME>", None),
            VariantClass::MobileElementDeletion
        );
        assert_eq!(
            classify_symbolic_allele("<DEL:ME:ALU>", None),
            VariantClass::MobileElementDeletion
        );
        assert_eq!(
            classify_symbolic_allele("<CPX>", None),
            VariantClass::ComplexStructural
        );
        assert_eq!(
            classify_symbolic_allele("<CNV:TR>", None),
            VariantClass::TandemRepeat
        );
        assert_eq!(
            classify_symbolic_allele("A[21:100[", None),
            VariantClass::Translocation
        );
        assert_eq!(
            classify_symbolic_allele("]21:100]A", None),
            VariantClass::Translocation
        );
        assert_eq!(
            classify_symbolic_allele("A.", None),
            VariantClass::Translocation
        );
        assert_eq!(
            classify_symbolic_allele(".A", None),
            VariantClass::Translocation
        );
    }

    #[test]
    fn test_classify_symbolic_allele_svtype_fallback() {
        // When ALT is a generic symbolic like <SV>, SVTYPE provides the real type.
        assert_eq!(
            classify_symbolic_allele("<SV>", Some("DEL")),
            VariantClass::StructuralDeletion
        );
        assert_eq!(
            classify_symbolic_allele("<SV>", Some("INS")),
            VariantClass::StructuralInsertion
        );
    }

    #[test]
    fn test_parse_bnd_alt_patterns() {
        assert_eq!(parse_bnd_alt("A[21:100["), Some(("21".to_string(), 100)));
        assert_eq!(parse_bnd_alt("A]21:200]"), Some(("21".to_string(), 200)));
        assert_eq!(parse_bnd_alt("]13:300]A"), Some(("13".to_string(), 300)));
        assert_eq!(parse_bnd_alt("[13:400[A"), Some(("13".to_string(), 400)));
        assert_eq!(parse_bnd_alt("<DEL>"), None);
        assert_eq!(parse_bnd_alt("A"), None);
    }

    #[test]
    fn test_get_raw_info_value() {
        assert_eq!(
            get_raw_info_value("SVTYPE=DEL;END=1100", "SVTYPE"),
            Some("DEL")
        );
        assert_eq!(
            get_raw_info_value("SVTYPE=DEL;END=1100", "END"),
            Some("1100")
        );
        assert_eq!(get_raw_info_value("SVTYPE=DEL;END=1100", "SVLEN"), None);
        assert_eq!(get_raw_info_value(".", "SVTYPE"), None);
        assert_eq!(get_raw_info_value("", "SVTYPE"), None);
        assert_eq!(
            get_raw_info_value("SVTYPE=BND;MATEID=bnd_2", "MATEID"),
            Some("bnd_2")
        );
    }

    #[test]
    fn test_parse_sv_del_with_negative_svlen() {
        // Deletions use negative SVLEN in many VCF callers
        let data = make_sv_vcf_data("21\t1000\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-100\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
        assert_eq!(v.sv_len, Some(-100));
        // end = POS + abs(SVLEN) = 1000 + 100 = 1100
        assert_eq!(v.end, 1100);
        assert_eq!(v.sv_end, Some(1100));
    }

    #[test]
    fn test_star_symbolic_filtered() {
        let data = make_vcf_data("21\t100\t.\tA\t<*>\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();
        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_non_ref_symbolic_passes_filter() {
        // <NON_REF> is not filtered: Perl VEP annotates it as an SV.
        let data = make_vcf_data("21\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();
        let variant = parser.next_variant();
        assert!(variant.is_some(), "<NON_REF> should not be filtered");
    }
}
