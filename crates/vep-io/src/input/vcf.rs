// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF input parser using noodles.

use std::io::BufRead;
use std::sync::Arc;

use noodles::vcf;

use vep_core::coordinate::parse_pos;
use vep_core::sv_record::{structural_variants, StructuralRecordError, DEFAULT_MAX_SV_SIZE};
use vep_core::sv_routing::{perl_true, routes_to_structural, variant_alts, written_alts};
use vep_core::variant::{InputVariant, VariantClass};
use vep_core::vcf_info::info_value;

use crate::error::IoError;

use super::InputParser;

/// VCF input parser that reads VCF records and converts them to `InputVariant`s.
///
/// For multi-allelic sites, each ALT allele produces a separate `InputVariant`.
/// VCF coordinates are converted to VEP-style coordinates (trimming shared prefixes).
pub struct VcfParser {
    reader: vcf::io::Reader<Box<dyn BufRead + Send>>,
    /// Buffer of variants from the current record (for multi-allelic sites).
    buffer: Vec<InputVariant>,
    /// The `--max_sv_size` a `<CNV:TR>` allele is expanded under (`None` for no cap),
    /// Ensembl VEP's default unless [`Self::with_max_sv_size`] sets it.
    max_sv_size: Option<u64>,
}

impl VcfParser {
    /// Create a new VCF parser from a buffered reader.
    pub fn new(reader: Box<dyn BufRead + Send>) -> Result<Self, IoError> {
        let mut vcf_reader = vcf::io::Reader::new(reader);
        vcf_reader
            .read_header()
            .map_err(|e| IoError::VcfParse(format!("failed to read VCF header: {e}")))?;
        Ok(Self {
            reader: vcf_reader,
            buffer: Vec::new(),
            max_sv_size: Some(DEFAULT_MAX_SV_SIZE),
        })
    }

    /// Sets the `--max_sv_size` a `<CNV:TR>` allele is expanded under (`None` for `-1`).
    pub fn with_max_sv_size(mut self, max_sv_size: Option<u64>) -> Self {
        self.max_sv_size = max_sv_size;
        self
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
    ///
    /// The record is read as the command-line parser reads it: POS from its own text
    /// (`pos_text`) as `parse_pos` reads it, INFO through `vep_core::vcf_info`, the ALT
    /// list as Perl's `split` returns it with a gVCF site read as its real alleles
    /// (`vep_core::sv_routing`), a record whose every ALT is `*` as the one allele-less
    /// variant Ensembl VEP carries, and a structural record's alleles, spans, breakends
    /// and numbering built by `vep_core::sv_record::structural_variants`; a structural
    /// type with no Sequence Ontology term and no symbolic ALT is skipped at warning
    /// level, as Ensembl VEP skips it. The sequence path trims the shared prefix of
    /// every allele and never minimises a bi-allelic indel from both ends, which the
    /// command-line parser does.
    fn parse_record(
        &self,
        record: &vcf::Record,
        pos_text: &str,
    ) -> Result<Vec<InputVariant>, IoError> {
        let chrom = record.reference_sequence_name().to_string();

        // POS from the line's own text, as the command-line parser reads it: noodles
        // reads `0` as no position, refuses a number past the machine word and never
        // numifies a spelling such as `3e7`.
        let pos = parse_pos(pos_text)
            .ok_or_else(|| IoError::VcfParse(format!("Invalid POS: '{pos_text}'")))?;

        let vcf_ref = record.reference_bases();

        let ids = record.ids();
        let ids_str: &str = ids.as_ref();
        let id = if ids_str.is_empty() {
            None
        } else {
            Some(ids_str.to_string())
        };

        let alt_bases = record.alternate_bases();
        let alt_column: &str = alt_bases.as_ref();
        let info = record.info();
        let info_raw: &str = info.as_ref();
        let svtype = info_value(info_raw, "SVTYPE").filter(|v| perl_true(v));
        let written = written_alts(alt_column);
        let mut alts = written.clone();
        variant_alts(&mut alts, svtype);
        // `*` and `.` name no allele; `<*>` is a sequence allele at POS unless INFO/SVTYPE
        // routes the record.
        let valid_alts: Vec<(usize, &str)> = alts
            .iter()
            .copied()
            .filter(|(_, alt)| *alt != "*" && *alt != ".")
            .collect();

        // OutputFactory.pm `VariationFeature_to_output_hash`, the VCF `_line`
        // branch: the ALT column as written, a filtered `*` allele included.
        let line_name: Arc<str> =
            format!("{chrom}_{pos}_{vcf_ref}/{}", alt_column.replace(',', "/")).into();

        if valid_alts.is_empty() {
            if written.is_empty() || written.iter().any(|(_, alt)| *alt != "*") {
                return Ok(Vec::new());
            }
            let mut variant = InputVariant::without_allele(
                chrom,
                pos,
                vcf_ref,
                written.iter().map(|(_, alt)| *alt),
            );
            variant.id = id;
            variant.input_line_name = Some(line_name);
            return Ok(vec![variant]);
        }

        let multi = valid_alts.len() > 1;
        let mut variants = Vec::with_capacity(valid_alts.len());

        if !routes_to_structural(vcf_ref, &alts, svtype) {
            let explicit_alts: Vec<Vec<u8>> = valid_alts
                .iter()
                .map(|(_, alt)| alt.as_bytes().to_vec())
                .collect();
            let vcf_ref_bytes = vcf_ref.as_bytes();
            let (mut trimmed_ref, mut trimmed_alts, prefix_trimmed) =
                trim_multi_alleles(vcf_ref_bytes, &explicit_alts);
            // `validate_vf` upper-cases the allele string (Parser.pm:554).
            trimmed_ref.make_ascii_uppercase();
            trimmed_alts
                .iter_mut()
                .for_each(|alt| alt.make_ascii_uppercase());

            for (local_idx, trimmed_alt) in trimmed_alts.into_iter().enumerate() {
                let allele_idx = valid_alts[local_idx].0;
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
        } else {
            variants = match structural_variants(
                &chrom,
                pos,
                vcf_ref,
                &written,
                info_raw,
                self.max_sv_size,
                id.as_deref(),
            ) {
                Ok(variants) => variants,
                Err(StructuralRecordError::Span(error)) => {
                    return Err(IoError::VcfParse(error.to_string()));
                }
                Err(error @ StructuralRecordError::UnsupportedType(_)) => {
                    tracing::warn!("{chrom}:{pos}: {error}");
                    return Ok(Vec::new());
                }
            };
            if variants
                .iter()
                .any(|v| v.variant_class == VariantClass::ReferenceBlock)
            {
                tracing::warn!(
                    "{chrom}:{pos} {alt_column}: NON_REF is not a supported structural variant \
                     type; a reference-confidence block is written as one intergenic row"
                );
            }
            for variant in &mut variants {
                variant.id = id.clone();
                variant.minimised = multi;
            }
        }

        for variant in &mut variants {
            variant.input_line_name = Some(line_name.clone());
        }

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

impl InputParser for VcfParser {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        if let Some(variant) = self.buffer.pop() {
            return Some(Ok(variant));
        }

        let mut line = Vec::new();
        loop {
            line.clear();
            match self.reader.get_mut().read_until(b'\n', &mut line) {
                Ok(0) => return None,
                Ok(_) => {
                    // The line terminator alone is removed, as Perl's `chomp` removes it
                    // (ensembl-io ColumnBasedParser.pm:249); a blank data line is skipped.
                    if line.ends_with(b"\n") {
                        line.pop();
                        if line.ends_with(b"\r") {
                            line.pop();
                        }
                    }
                    if line.is_empty() {
                        continue;
                    }
                    let record = match vcf::Record::try_from(&line[..]) {
                        Ok(record) => record,
                        Err(e) => {
                            return Some(Err(IoError::VcfParse(format!(
                                "error reading VCF record: {e}"
                            ))))
                        }
                    };
                    let pos_text = std::str::from_utf8(&line)
                        .ok()
                        .and_then(|text| text.split('\t').nth(1))
                        .unwrap_or("");
                    match self.parse_record(&record, pos_text) {
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
    /// noodles needs no INFO header definitions to parse the records; the reader
    /// reads INFO from the raw column (`vep_core::vcf_info`).
    fn make_sv_vcf_data(records: &str) -> Vec<u8> {
        make_vcf_data(records)
    }

    /// Every variant the reader yields for one data line.
    fn read_record(record: &str) -> Vec<InputVariant> {
        let data = make_sv_vcf_data(&format!("{record}\n"));
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();
        let mut variants = Vec::new();
        while let Some(v) = parser.next_variant() {
            variants.push(v.unwrap());
        }
        variants
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
    fn test_parse_sv_del_with_end() {
        let data = make_sv_vcf_data("21\t1000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1100\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
        // The base before the event is the anchor (BaseVCF4.pm get_start).
        assert_eq!(v.start, 1001);
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
        assert_eq!(v.start, 501);
        // No END or SVLEN -> the single base at POS + 1
        assert_eq!(v.end, 501);
    }

    /// The reader shares `vep_core::sv_span::structural_end` with the command-line
    /// parser: SVLEN first, by its first value as Perl numifies it (`abc` is 0, the pair
    /// `start = POS + 1, end = POS`), else END, else the start; a repeated key reads its
    /// last value; an END before POS or not a coordinate fails the line as Ensembl VEP
    /// skips the record.
    #[test]
    fn the_reader_reads_the_span_as_the_command_line_parser_does() {
        let span = |info: &str| {
            let data = make_sv_vcf_data(&format!("21\t1000\t.\tA\t<DUP>\t.\t.\t{info}\n"));
            let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
            let mut parser = VcfParser::new(reader).unwrap();
            parser
                .next_variant()
                .unwrap()
                .map(|v| (v.start, v.end, v.sv_len))
        };
        assert_eq!(
            span("SVTYPE=DUP;SVLEN=-300;END=1100").unwrap(),
            (1001, 1300, Some(-300))
        );
        assert_eq!(
            span("SVTYPE=DUP;SVLEN=abc;END=1100").unwrap(),
            (1001, 1000, Some(0))
        );
        assert_eq!(
            span("SVTYPE=DUP;SVLEN=300abc").unwrap(),
            (1001, 1300, Some(300))
        );
        assert_eq!(
            span("SVTYPE=DUP;SVLEN=-100;SVLEN=-300").unwrap(),
            (1001, 1300, Some(-300))
        );
        assert_eq!(
            span("SVTYPE=DUP;END=1100;END=1200").unwrap(),
            (1001, 1200, None)
        );
        assert_eq!(span("SVTYPE=DUP;END=1000").unwrap(), (1001, 1000, None));
        for info in [
            "SVTYPE=DUP;END=0",
            "SVTYPE=DUP;END=900",
            "SVTYPE=DUP;END=abc",
            "SVTYPE=DUP;END=",
        ] {
            let err = span(info).unwrap_err().to_string();
            assert!(
                err.contains("start > end+1") || err.contains("Invalid END"),
                "{info}: {err}"
            );
        }
        let bnd = make_sv_vcf_data("21\t1000\t.\tA\t<BND>\t.\t.\tSVTYPE=BND;END=1100\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(bnd));
        let v = VcfParser::new(reader)
            .unwrap()
            .next_variant()
            .unwrap()
            .unwrap();
        assert_eq!(
            (v.start, v.end),
            (1001, 1001),
            "a breakpoint never reads END for its span"
        );
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
        // SVLEN gives the end, POS + abs(SVLEN)
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

    /// A mixed record is one structural variant (Parser/VCF.pm:237-245 routes the whole
    /// line; one StructuralVariationFeature carries it, `:560-574`), so the sequence ALT
    /// of `G,<DEL>` with `SVTYPE=DEL` is a deletion allele over END like the symbolic one,
    /// and `G,<INV>` with `SVTYPE=DEL` takes SVTYPE's type, the joined ALTs opening with a
    /// sequence allele (`:472-475`).
    #[test]
    fn a_mixed_record_is_one_structural_variant_under_the_record_type() {
        let v = read_record("21\t1000\t.\tA\tG,<DEL>\t.\t.\tSVTYPE=DEL;END=1100");
        assert_eq!(v.len(), 2);
        for (variant, index) in v.iter().zip(0..) {
            assert_eq!(variant.allele_index, index);
            assert_eq!(variant.variant_class, VariantClass::StructuralDeletion);
            assert!(variant.is_structural && variant.minimised);
            assert_eq!(
                (variant.start, variant.end, variant.sv_end),
                (1001, 1100, Some(1100))
            );
            assert_eq!(variant.record_so_term.as_deref(), Some("deletion"));
        }
        assert_eq!(v[0].alt_allele(), b"G");
        assert_eq!(v[1].alt_allele(), b"<DEL>");
        let v = read_record("21\t1000\t.\tC\tG,<INV>\t.\t.\tSVTYPE=DEL;END=2000;SVLEN=-1000");
        assert_eq!(v.len(), 2);
        for variant in &v {
            assert_eq!(variant.variant_class, VariantClass::StructuralDeletion);
            assert_eq!(variant.sv_end, Some(2000));
        }
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
        assert_eq!(v.start, 101);
        assert_eq!(v.end, 200);
    }

    /// The reader routes, types and classes a record through `vep_core::sv_routing`, as
    /// the command-line parser does: Ensembl VEP 116.2 writes a `deletion` over END for
    /// `N G`, `<*>` and `<SV>` with `SVTYPE=DEL`, the SNV `G` for `A G` with the same INFO
    /// and for `N G` alone, `duplication` for `<del>` with `SVTYPE=DUP`, a `chromosome_breakpoint`
    /// for each part of `<BND>,<DEL>` with `SVTYPE=BND`, `mobile_element_insertion` for
    /// `<INS:ME::ALU>`, and skips `<del>`, `<CN=X>`, `<NON_REF>`, `<NON_REF>,<NON_REF>`
    /// and `N G` with `SVTYPE=FOO` as unsupported types; vep-rs drops the gVCF site's
    /// `<NON_REF>` beside `G`, reads a block list as one reference block per allele and
    /// writes nothing for the unsupported `N G`.
    #[test]
    fn the_reader_routes_a_record_as_the_command_line_parser_does() {
        let summary = |v: &InputVariant| {
            (
                v.is_structural,
                v.variant_class,
                v.vep_skip,
                v.record_so_term.clone(),
                v.start,
                v.end,
            )
        };
        let deletion = (
            true,
            VariantClass::StructuralDeletion,
            false,
            Some("deletion".to_string()),
            1001,
            1100,
        );
        for alt in ["G", "<*>", "<SV>"] {
            let v = read_record(&format!("21\t1000\t.\tN\t{alt}\t.\t.\tSVTYPE=DEL;END=1100"));
            assert_eq!(v.len(), 1, "{alt}");
            assert_eq!(summary(&v[0]), deletion, "{alt}");
        }
        for record in [
            "21\t1000\t.\tA\tG\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\t.\tN\tG\t.\t.\t.",
            "21\t1000\t.\tN\t<*>\t.\t.\t.",
            "21\t1000\t.\tA\tG,<NON_REF>\t.\t.\t.",
        ] {
            let v = read_record(record);
            assert_eq!(v.len(), 1, "{record}");
            assert_eq!(
                (v[0].is_structural, v[0].start, v[0].end, v[0].allele_index),
                (false, 1000, 1000, 0),
                "{record}"
            );
        }
        let v = read_record("21\t1000\t.\tN\t<del>\t.\t.\tSVTYPE=DUP;END=1100");
        assert_eq!(v[0].record_so_term.as_deref(), Some("duplication"));
        assert_eq!(v[0].variant_class, VariantClass::Duplication);
        let v = read_record("21\t1000\t.\tN\t<BND>,<DEL>\t.\t.\tSVTYPE=BND;END=1100");
        assert_eq!(v.len(), 2);
        for variant in &v {
            assert_eq!(variant.variant_class, VariantClass::Translocation);
            assert_eq!((variant.start, variant.end), (1001, 1001));
            assert_eq!(
                variant.uploaded_allele_string.as_deref(),
                Some("<BND>/<DEL>")
            );
        }
        let v = read_record("21\t1000\t.\tN\t<INS:ME::ALU>\t.\t.\tSVTYPE=INS");
        assert_eq!(
            v[0].record_so_term.as_deref(),
            Some("mobile_element_insertion")
        );
        for (record, count, class) in [
            (
                "21\t1000\t.\tN\t<del>\t.\t.\tEND=1100",
                1,
                VariantClass::StructuralDeletion,
            ),
            (
                "21\t1000\t.\tN\t<CN=X>\t.\t.\tEND=1100",
                1,
                VariantClass::CopyNumberVariation,
            ),
            (
                "21\t1000\t.\tN\t<NON_REF>\t.\t.\tEND=1100",
                1,
                VariantClass::ReferenceBlock,
            ),
            (
                "21\t1000\t.\tN\t<NON_REF>,<NON_REF>\t.\t.\tEND=1100",
                2,
                VariantClass::ReferenceBlock,
            ),
        ] {
            let v = read_record(record);
            assert_eq!(v.len(), count, "{record}");
            for variant in &v {
                assert!(
                    variant.vep_skip && variant.record_so_term.is_none(),
                    "{record}"
                );
                assert_eq!(variant.variant_class, class, "{record}");
            }
        }
        assert!(read_record("21\t1000\t.\tN\tG\t.\t.\tSVTYPE=FOO;END=1100").is_empty());
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

    /// `<*>` is outside Perl's structural ALT regex (`[<\[\]][^\*]+[>\]\[]`), so without
    /// INFO/SVTYPE it is a sequence allele at POS, as the command-line parser reads it.
    #[test]
    fn a_star_allele_alone_is_a_sequence_allele_at_pos() {
        let data = make_vcf_data("21\t100\t.\tA\t<*>\t.\t.\t.\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();
        let v = parser.next_variant().unwrap().unwrap();
        assert!(!v.is_structural);
        assert_eq!((v.start, v.end), (100, 100));
        assert_eq!(v.alt_allele(), b"<*>");
        assert!(parser.next_variant().is_none());
    }

    /// The POS column is read from the line's own text, as Ensembl VEP numifies it
    /// (`21 3e7 x N <DEL> . . SVTYPE=DEL;END=30425949` is the deletion from 30000001);
    /// alleles are upper-cased after the anchor trim (Parser.pm:554: `agc a` is `GC/-`);
    /// a record whose every ALT is `*` is the one allele-less feature Ensembl VEP
    /// carries, over the reference bases; and a blank closing the line stays in the
    /// last INFO value, so `END2=1700 ` names the mate `N[21:1700 [`.
    #[test]
    fn pos_spellings_case_star_only_records_and_line_ends_follow_the_command_line_parser() {
        let deletion = read_record("21\t3e7\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425949");
        assert_eq!((deletion[0].start, deletion[0].end), (30000001, 30425949));
        let snv = read_record("21\t30425949.0\tx\tA\tG\t.\t.\t.");
        assert_eq!((snv[0].start, snv[0].end), (30425949, 30425949));
        let upper = read_record("21\t100\tx\tagc\ta\t.\t.\t.");
        assert_eq!(upper[0].allele_string, "GC/-");
        let star = read_record("21\t100\tx\tNAC\t*,*\t.\t.\t.");
        assert_eq!(star.len(), 1);
        assert!(star[0].no_allele);
        assert_eq!(star[0].allele_string, "NAC/*/*");
        assert_eq!((star[0].start, star[0].end), (100, 102));
        let blank = read_record("21\t100\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=1700 ");
        assert_eq!(blank[1].alt_allele(), b"N[21:1700 [");
        assert_eq!(blank[1].mate_pos, Some(1700));
    }

    #[test]
    fn test_non_ref_symbolic_is_a_reference_block() {
        // <NON_REF> is not filtered: the record is one reference-confidence block.
        let data = make_vcf_data("21\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200\n");
        let reader: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(data));
        let mut parser = VcfParser::new(reader).unwrap();
        let variant = parser
            .next_variant()
            .expect("<NON_REF> should not be filtered")
            .unwrap();
        assert_eq!(variant.variant_class, VariantClass::ReferenceBlock);
        assert!(variant.is_structural);
    }
}
