// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The writer of the records and of the summary: a VCF of the eight fixed columns under
//! a header that names the generator, and one tab-separated line a category.

use std::io::{self, Write};

use crate::GENERATOR_VERSION;

/// One input record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// The sequence, as the record names it.
    pub sequence: String,
    /// The position of the reference allele's first base, counted from 1.
    pub position: u64,
    /// The category's name, an underscore and the record's ordinal in its category.
    pub id: String,
    pub reference: String,
    pub alt: String,
}

/// What the header says of a run beside the generator's name and version.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    pub seed: u64,
    /// `light`, or `wide` and how many records a category.
    pub records: String,
    /// Where the reference alleles were read: `cache` or `fasta`.
    pub reference_alleles: String,
}

/// The header line of the summary.
pub const SUMMARY_HEADER: &str = "category\twanted\twritten\tpassed_over\treason";

/// One line of the summary: a category, how many records it wanted and wrote, how many
/// candidates the reference file made it pass over and, where it wrote fewer than it
/// wanted, why; `-` where it wrote them all.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SummaryLine {
    pub category: String,
    pub wanted: usize,
    pub written: usize,
    pub passed_over: u64,
    pub reason: String,
}

/// Writes the header and the records, in the order given.
pub fn write_vcf(out: &mut impl Write, header: &Header, records: &[Record]) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(
        out,
        "##source=vep-seeded-variants generator_version={GENERATOR_VERSION}"
    )?;
    writeln!(out, "##seed={}", header.seed)?;
    writeln!(out, "##records={}", header.records)?;
    writeln!(out, "##reference_alleles={}", header.reference_alleles)?;
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    for record in records {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t.\t.\t.",
            record.sequence, record.position, record.id, record.reference, record.alt
        )?;
    }
    Ok(())
}

/// Writes the summary: its header line, then one line a category.
pub fn write_summary(out: &mut impl Write, lines: &[SummaryLine]) -> io::Result<()> {
    writeln!(out, "{SUMMARY_HEADER}")?;
    for line in lines {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            line.category, line.wanted, line.written, line.passed_over, line.reason
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_header_names_the_generator_and_the_run() {
        let header = Header {
            seed: 7,
            records: "light".to_owned(),
            reference_alleles: "cache".to_owned(),
        };
        let record = Record {
            sequence: "Mito".to_owned(),
            position: 12,
            id: "intron_1".to_owned(),
            reference: "A".to_owned(),
            alt: "C".to_owned(),
        };
        let mut out = Vec::new();
        write_vcf(&mut out, &header, &[record]).expect("a written file");
        let text = String::from_utf8(out).expect("text");
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(
            lines[1],
            format!("##source=vep-seeded-variants generator_version={GENERATOR_VERSION}")
        );
        assert_eq!(lines[2], "##seed=7");
        assert_eq!(lines[5], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        assert_eq!(lines[6], "Mito\t12\tintron_1\tA\tC\t.\t.\t.");
        assert_eq!(lines.len(), 7);
    }

    #[test]
    fn the_summary_is_one_line_a_category() {
        let line = SummaryLine {
            category: "origin".to_owned(),
            wanted: 1,
            written: 0,
            passed_over: 0,
            reason: "no sequence length".to_owned(),
        };
        let mut out = Vec::new();
        write_summary(&mut out, &[line]).expect("a written file");
        assert_eq!(
            String::from_utf8(out).expect("text"),
            format!("{SUMMARY_HEADER}\norigin\t1\t0\t0\tno sequence length\n")
        );
    }
}
