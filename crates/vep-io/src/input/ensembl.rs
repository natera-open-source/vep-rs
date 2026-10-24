// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's default input format, read as `Parser/VEP_input.pm`
//! `create_VariationFeatures` (release/116.2, lines 109-168) reads it over
//! ensembl-io's `VEP_input` parser:
//!
//! ```text
//! chr  start  end  allele_string  [strand]  [identifier]
//! 21   25585733  25585733  A/G  1  rs699
//! 21   25585656  25585660  ACGTG/-  -  .
//! 21   25585700  25585900  DUP
//! ```
//!
//! The columns are split on the run's delimiter ([`Delimiter`], a tab unless
//! `--delimiter` or format detection chose another) and each field loses its
//! trailing whitespace (`:121-128`). The strand is reverse when its token
//! contains a `-` (`-`, `-1`) and forward for any other token or none
//! (`:131`, over ensembl-io's `+`, `.` and `-` conversion). An allele string
//! without a `/` names a structural variant class (`DUP`, `DEL`, `INS`, `INV`,
//! `CNV`, `TDUP`, ...) over the line's start and end (`:136-148`, through
//! `get_SO_term`); one with a `/` is a sequence variant whose alleles are
//! upper-cased as `validate_vf` upper-cases them (Parser.pm:554) and, for a
//! bi-allelic record whose alleles differ in length, minimised as
//! `post_process_vfs` minimises them ([`minimise_alleles`]). A multi-allelic
//! string (`A/G/T`) yields one `InputVariant` per alternate allele.
//!
//! Ensembl VEP dies on a line of fewer than four columns or an empty line
//! (`:136`: its parser has no allele string to read) and reads a `#` line as a
//! record (ensembl-io's VEP_input parser has no metadata lines,
//! VEP_input.pm:61), so a `#` header of four columns is skipped with a
//! coordinate warning and a shorter one dies; either way the line counts
//! toward the line numbers that follow (Parser/VEP_input.pm:119). Here a `#`
//! line is a comment, counted the same way and otherwise passed over without
//! a word (`LineCount`); an empty or whitespace-only line is counted and
//! warned about in `skip_empty_lines`'s words, as the list formats' lines are
//! (Parser.pm:273-290), where Ensembl VEP's own parser dies; and a short line
//! is a parse error the caller warns about and skips.

use std::io::BufRead;
use std::sync::Arc;

use vep_core::coordinate::{parse_coordinate, Strand};
use vep_core::minimal::minimise_bi_allelic;
use vep_core::sv_routing::{perl_true, sv_so_term_and_class, sv_type_abbreviation};
use vep_core::variant::{classify_variant, InputVariant, VariantClass};

use crate::error::IoError;

use super::{is_perl_space, minimise_alleles, Delimiter, InputParser, LineCount, ParserContext};

/// Parser for Ensembl VEP's default input format.
pub struct EnsemblParser<R: BufRead> {
    reader: R,
    delimiter: Delimiter,
    /// `--minimal` ([`ParserContext::minimal`]).
    minimal: bool,
    line_buf: Vec<u8>,
    /// Buffer for multi-allelic variants from a single line.
    buffer: Vec<InputVariant>,
    /// Ensembl VEP's line number over the lines read, and the empty-line
    /// warnings it raises.
    lines: LineCount,
}

impl<R: BufRead> EnsemblParser<R> {
    pub fn new(reader: R, context: &ParserContext) -> Self {
        Self {
            reader,
            delimiter: context.delimiter.clone(),
            minimal: context.minimal,
            line_buf: Vec::new(),
            buffer: Vec::new(),
            lines: LineCount::default(),
        }
    }

    /// Parse one record line, its line end removed, into one `InputVariant`
    /// per alternate allele.
    fn parse_record(&self, line: &str) -> Result<Vec<InputVariant>, IoError> {
        let record = self.delimiter.split(line);
        let fields: Vec<&str> = record
            .iter()
            .map(|field| field.trim_end_matches(is_perl_space))
            .collect();
        if fields.len() < 4 {
            return Err(IoError::VcfParse(format!(
                "Ensembl format requires at least 4 columns (chr, start, end, allele string), got {}",
                fields.len()
            )));
        }

        let chr = fields[0].to_string();
        let start = parse_coordinate(fields[1])
            .ok_or_else(|| IoError::VcfParse(format!("Invalid start position: '{}'", fields[1])))?;
        let end = parse_coordinate(fields[2])
            .ok_or_else(|| IoError::VcfParse(format!("Invalid end position: '{}'", fields[2])))?;
        let allele_string = fields[3];
        let strand = strand_from_token(fields.get(4).copied());
        // The ID column whole, `.` included, as Parser/VEP_input.pm:127 takes it; a
        // `;` is part of the name here, where a VCF ID field is split on it.
        let id = fields
            .get(5)
            .filter(|id| !id.is_empty())
            .map(|id| id.to_string());
        // The line as Ensembl VEP re-emits it (`OutputFactory/JSON.pm:247`): the
        // fields as split, before the trim, joined by the delimiter.
        let raw_input: Arc<str> = record.join(self.delimiter.join_text()).into();
        // OutputFactory.pm `VariationFeature_to_output_hash`, the ensembl `_line`
        // branch: the lower coordinate, so an insertion written `100 99` is named
        // at 99, and the chromosome and allele columns as split, before the trim.
        let line_name: Arc<str> = format!("{}_{}_{}", record[0], start.min(end), record[3]).into();

        let mut variants = if allele_string.contains('/') {
            sequence_variants(chr, start, end, allele_string, id, self.minimal)?
        } else {
            let mut variant = structural_variant(chr, start, end, allele_string);
            variant.id = id;
            vec![variant]
        };
        for variant in &mut variants {
            variant.strand = strand;
            variant.raw_input = Some(raw_input.clone());
            variant.input_line_name = Some(line_name.clone());
        }
        Ok(variants)
    }
}

/// The strand of an Ensembl-format record (Parser/VEP_input.pm:131, over
/// ensembl-io `Parser.pm:69` and `VEP_input.pm:130-135`): the column is read
/// as `1` when absent, `+` as `1`, `.` as `0` and `-` as `-1`, then the record
/// is reverse when the result contains a `-` and forward otherwise, so `-` and
/// `-1` are reverse and `+`, `1`, `+1`, `.`, `0`, an absent column and any word
/// without a `-` are forward.
fn strand_from_token(token: Option<&str>) -> Strand {
    match token {
        Some(token) if token.contains('-') => Strand::Reverse,
        _ => Strand::Forward,
    }
}

/// The alleles of a record whose allele string has a `/`, one variant per
/// alternate allele, each carrying the record's `id`. Ensembl VEP keeps the
/// whole column as the record's allele string (Parser/VEP_input.pm:156),
/// minimises a bi-allelic indel in `post_process_vfs` on the bases as written,
/// and under `minimal` (`--minimal`) every bi-allelic record (Parser.pm:852),
/// then upper-cases the string in `validate_vf` (Parser.pm:554), in that order
/// here too: `AAAG/AGGg` under the flag is `AAG/GGg` over the last three bases,
/// the `G` and `g` left unmatched by the written-case comparison, and then
/// `AAG/GGG`. An empty allele (`A/`, `/G`) fails the line: Ensembl VEP builds a
/// feature with one allele from such a string and writes nothing for it.
///
/// A record whose ID Perl reads as false is named `chr_start_alleles` at
/// Parser.pm:462, between the minimisation and the upper-casing: for alleles
/// written in lower case that name differs from the upper-cased one
/// `VariationFeature_to_output_hash` compares it with, so it is what
/// `Uploaded_variation` and the JSON `id` print, and it is stored here as the
/// record's id. A record written in upper case takes the generated name through
/// [`InputVariant::variation_name`] as Ensembl VEP builds it. The identifier
/// parser builds a mapping's alleles through this too, since `post_process_vfs`
/// runs on its features alike (Parser/ID.pm:196). The string as written is kept
/// beside the upper-cased alleles ([`InputVariant::written_case_allele_string`])
/// for a record `--dont_skip` keeps after a check that returns before the
/// upper-casing refused it (Parser.pm:470-472, `:542-546`), which Ensembl VEP
/// then annotates in the written case.
pub(super) fn sequence_variants(
    chr: String,
    start: u64,
    end: u64,
    allele_string: &str,
    id: Option<String>,
    minimal: bool,
) -> Result<Vec<InputVariant>, IoError> {
    if allele_string.split('/').any(str::is_empty) {
        return Err(IoError::VcfParse(format!(
            "Invalid allele string: '{allele_string}' (an allele is empty)"
        )));
    }
    let mut alleles = allele_string.split('/');
    let ref_allele = alleles.next().unwrap_or_default();
    let alts: Vec<&str> = alleles.collect();
    let multi_allelic = alts.len() > 1;
    let upper = allele_string.to_ascii_uppercase();
    let record_allele_string: Option<Arc<str>> = multi_allelic.then(|| upper.as_str().into());
    // `UPLOADED_ALLELE` reads the line's one VariationFeature (OutputFactory.pm:1360),
    // so every allele keeps the whole column as written, not its own pair.
    let uploaded_allele_string: Arc<str> = allele_string.into();
    let named_from_lower_case_alleles =
        upper != allele_string && !id.as_deref().is_some_and(perl_true);
    // The string Parser.pm:462 reads and `:554` upper-cases is the whole column
    // for a multi-allelic record, shared by its alleles; a bi-allelic record's
    // is its one minimised pair, built in the loop.
    let written_column = if multi_allelic {
        written_case(allele_string)
    } else {
        None
    };

    let mut variants = Vec::with_capacity(alts.len());
    for (i, alt) in alts.iter().enumerate() {
        let mut variant = InputVariant::new(
            chr.clone(),
            start,
            end,
            ref_allele.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant.allele_index = i;
        variant.minimised = multi_allelic;
        variant.record_allele_string_multi = record_allele_string.clone();
        variant.uploaded_allele_string = Some(uploaded_allele_string.clone());
        if !multi_allelic {
            minimise_alleles(&mut variant);
            if minimal {
                minimise_bi_allelic(&mut variant);
            }
        }
        let (working, written) = if multi_allelic {
            (allele_string, written_column.clone())
        } else {
            let pair = variant.allele_string.as_str();
            (pair, written_case(pair))
        };
        variant.id = if named_from_lower_case_alleles {
            Some(format!("{chr}_{}_{working}", variant.start))
        } else {
            id.clone()
        };
        variant.written_case_allele_string = written;
        upper_case_alleles(&mut variant);
        variants.push(variant);
    }
    Ok(variants)
}

/// `text` as its own `Arc` when a base is lower case, the case the input wrote
/// ([`InputVariant::written_case_allele_string`]); `None` for upper-case text.
fn written_case(text: &str) -> Option<Arc<str>> {
    text.bytes()
        .any(|b| b.is_ascii_lowercase())
        .then(|| text.into())
}

/// `validate_vf`'s `tr/[a-z]/[A-Z]/` on the allele string (Parser.pm:554),
/// applied to the alleles the string is built from; the class is read again,
/// since `classify_variant` reads the bases.
pub(super) fn upper_case_alleles(variant: &mut InputVariant) {
    variant.ref_allele.make_ascii_uppercase();
    for alt in &mut variant.alt_alleles {
        alt.make_ascii_uppercase();
    }
    variant.allele_string.make_ascii_uppercase();
    variant.variant_class = classify_variant(&variant.ref_allele, variant.alt_allele());
}

/// The structural variant a record whose allele string has no `/` names
/// (Parser/VEP_input.pm:136-148): the class is the `%SO_TERMS` entry for the
/// word as `get_SO_term` abbreviates it (Parser.pm:663-728), over the line's
/// start and end as written. A word the table lacks is Ensembl VEP's own class
/// term for the record, warned as `<word> is not a supported structural variant
/// type` (Parser.pm:725) and annotated under it in every format: the record
/// carries no `record_so_term`, which is the caller's cue to write the warning.
pub(super) fn structural_variant(chr: String, start: u64, end: u64, word: &str) -> InputVariant {
    let term = sv_so_term_and_class(&sv_type_abbreviation(word, None));
    let mut variant = InputVariant::new(chr, start, end, b"N".to_vec(), word.as_bytes().to_vec());
    variant.is_structural = true;
    variant.variant_class = term.map_or(VariantClass::ComplexStructural, |(_, class)| class);
    variant.record_so_term = term.map(|(term, _)| term.to_string());
    // Ensembl VEP's feature has no allele string, so its class term stands in for
    // it wherever one is printed (`$vf->{allele_string} || $vf->{class_SO_term}`:
    // the generated name, the JSON `allele_string`, `UPLOADED_ALLELE`).
    variant.uploaded_allele_string = Some(term.map_or(word, |(term, _)| term).into());
    variant.sv_end = Some(end);
    variant
}

impl<R: BufRead + Send> InputParser for EnsemblParser<R> {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        if let Some(variant) = self.buffer.pop() {
            return Some(Ok(variant));
        }

        loop {
            self.line_buf.clear();
            match self.reader.read_until(b'\n', &mut self.line_buf) {
                Ok(0) => {
                    self.lines.end_of_input();
                    return None;
                }
                Ok(_) => {}
                Err(e) => return Some(Err(IoError::Io(e))),
            }
            let text = String::from_utf8_lossy(&self.line_buf);
            // ensembl-io's `TextParser` turns `\r\n` into `\n` (:138) and
            // `chomp` removes it; the fields lose any other trailing
            // blank on their own.
            let line = match text.strip_suffix('\n') {
                Some(without_newline) => without_newline
                    .strip_suffix('\r')
                    .unwrap_or(without_newline),
                None => &text,
            };
            if line.trim_matches(is_perl_space).is_empty() {
                self.lines.empty_line();
                continue;
            }
            self.lines.record_line();
            if line.starts_with('#') {
                continue;
            }
            let mut variants = match self.parse_record(line) {
                Ok(variants) => variants,
                Err(e) => return Some(Err(e)),
            };
            variants.reverse();
            let first = variants.pop().expect("a record has an allele");
            self.buffer = variants;
            return Some(Ok(first));
        }
    }

    fn take_warnings(&mut self) -> Vec<String> {
        self.lines.take_warnings()
    }

    fn record_line_number(&self) -> Option<u64> {
        Some(self.lines.line())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use vep_core::coordinate::MAX_COORD;

    fn parser_under(delimiter: Delimiter, data: &str) -> EnsemblParser<Cursor<Vec<u8>>> {
        EnsemblParser::new(
            Cursor::new(data.as_bytes().to_vec()),
            &ParserContext::with_delimiter(delimiter),
        )
    }

    fn parser(data: &str) -> EnsemblParser<Cursor<Vec<u8>>> {
        parser_under(Delimiter::Tab, data)
    }

    fn all(data: &str) -> Vec<InputVariant> {
        let mut parser = parser(data);
        let mut out = Vec::new();
        while let Some(v) = parser.next_variant() {
            out.push(v.unwrap());
        }
        out
    }

    #[test]
    fn test_parse_snv() {
        let mut parser = parser("21\t25585733\t25585733\tA/G\t1\trs699\n");

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.end, 25585733);
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_allele(), b"G");
        assert_eq!(v.strand, Strand::Forward);
        assert_eq!(v.id, Some("rs699".to_string()));
        assert_eq!(
            v.raw_input.as_deref(),
            Some("21\t25585733\t25585733\tA/G\t1\trs699")
        );

        assert!(parser.next_variant().is_none());
    }

    /// Parser/VEP_input.pm:131 over ensembl-io's conversion of `+`, `.` and `-`:
    /// the record is reverse when the strand token contains a `-`, so `-` and
    /// `-1` are reverse and `+`, `1`, `+1`, `.`, `0` and a word such as `minus`
    /// are forward. Ensembl VEP 116.2 reads `-` and `-1` as reverse and `minus`
    /// as forward.
    #[test]
    fn strand_tokens_plus_minus_one() {
        let cases = [
            ("-", Strand::Reverse),
            ("-1", Strand::Reverse),
            ("+", Strand::Forward),
            ("1", Strand::Forward),
            ("+1", Strand::Forward),
            (".", Strand::Forward),
            ("0", Strand::Forward),
            ("minus", Strand::Forward),
        ];
        for (token, strand) in cases {
            let v = all(&format!("21\t25585733\t25585733\tA/G\t{token}\n"));
            assert_eq!(v.len(), 1, "{token:?}");
            assert_eq!(v[0].strand, strand, "{token:?}");
        }
    }

    /// ensembl-io `VEP_input.pm:133` reads an absent strand column as `1`, so a
    /// four-column line is a forward-strand record with no identifier (Ensembl
    /// VEP 116.2 annotates `21 25585733 25585733 A/G` as forward); three columns
    /// leave no allele string to read and fail the line.
    #[test]
    fn four_column_line_defaults_forward() {
        let v = all("21\t25585733\t25585733\tA/G\n");
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].strand, Strand::Forward);
        assert_eq!(v[0].id, None);
        assert_eq!(v[0].uploaded_variation(), "21_25585733_A/G");
        let error = parser("21\t25585733\t25585733\n")
            .next_variant()
            .unwrap()
            .unwrap_err()
            .to_string();
        assert!(error.contains("at least 4 columns"), "{error}");
    }

    /// A refusal is the reason alone, no format named: an empty allele is skipped
    /// with the text `Invalid allele string: 'A/' (an allele is empty)` whether
    /// the line is an Ensembl-format one or an identifier's mapping
    /// (`sequence_variants` serves both).
    #[test]
    fn a_refusal_names_the_field_and_no_format() {
        let error = parser("21\t25585733\t25585733\tA/\t1\n")
            .next_variant()
            .unwrap()
            .unwrap_err()
            .to_string();
        assert_eq!(error, "Invalid allele string: 'A/' (an allele is empty)");
        let error = sequence_variants("21".into(), 100, 100, "/G", None, false)
            .unwrap_err()
            .to_string();
        assert_eq!(error, "Invalid allele string: '/G' (an allele is empty)");
    }

    /// Parser/VEP_input.pm:136-148: an allele string without a `/` is a
    /// structural variant class word read through `get_SO_term`, over the line's
    /// start and end; its class term is the Allele and `VARIANT_CLASS`, stands
    /// in for the allele string in the generated name (Parser.pm:462-467,
    /// `$vf->{allele_string} || $vf->{class_SO_term}`), and the ID-less record
    /// is named from its line. A word `%SO_TERMS` lacks is kept as the class and
    /// the record is annotated, not skipped (Ensembl VEP warns `FOO is not a
    /// supported structural variant type` and annotates it, Parser.pm:725).
    #[test]
    fn sv_allele_word_is_structural() {
        let v = all("21\t25585700\t25585900\tDUP\t1\n21\t25585700\t25585900\tdel\t-\tsv1\n21\t100\t200\tFOO\n");
        assert_eq!(v.len(), 3);
        let dup = &v[0];
        assert!(dup.is_structural);
        assert_eq!((dup.start, dup.end), (25585700, 25585900));
        assert_eq!(dup.variant_class, VariantClass::Duplication);
        assert_eq!(dup.record_so_term.as_deref(), Some("duplication"));
        assert_eq!(dup.display_allele(), "duplication");
        assert_eq!(dup.class_so_term(), "duplication");
        assert_eq!(dup.variation_name(), "21_25585700_duplication");
        assert_eq!(dup.uploaded_variation(), "21_25585700_DUP");
        assert_eq!(dup.sv_end, Some(25585900));
        assert!(!dup.vep_skip);
        // `get_SO_term` keeps a bare word as written (Parser.pm:717), so `del` is
        // no `%SO_TERMS` key, and the strand token `-` is read as it is on a
        // sequence record.
        let del = &v[1];
        assert_eq!(del.strand, Strand::Reverse);
        assert_eq!(del.id.as_deref(), Some("sv1"));
        assert!(!del.vep_skip);
        assert_eq!(del.record_so_term, None);
        let foo = &v[2];
        assert!(foo.is_structural);
        assert!(!foo.vep_skip);
        assert_eq!(foo.record_so_term, None);
        assert_eq!(foo.display_allele(), "FOO");
        assert_eq!(foo.variation_name(), "21_100_FOO");
        assert_eq!(foo.uploaded_variation(), "21_100_FOO");
    }

    /// Parser.pm:862-881 minimises a bi-allelic record whose alleles differ in
    /// length and carry no `-`: `AC/A` at 100-101 becomes `C/-` at 101 and
    /// `A/AT` at 100 becomes `-/T` at 101-100, the written string kept for
    /// `UPLOADED_ALLELE` and the line's name for `Uploaded_variation`; the JSON
    /// `id` is the generated name of the minimised record. `A/G/AT` stays as
    /// written: Ensembl VEP leaves a multi-allelic string to `--minimal`.
    #[test]
    fn ac_a_is_auto_minimised() {
        let v = all("21\t100\t101\tAC/A\t1\n21\t100\t100\tA/AT\t1\n21\t100\t100\tA/G/AT\t1\n");
        assert_eq!(v.len(), 4);
        let del = &v[0];
        assert_eq!((del.start, del.end), (101, 101));
        assert_eq!(del.allele_string, "C/-");
        assert_eq!(del.ref_allele, b"C");
        assert_eq!(del.alt_allele(), b"-");
        assert!(del.minimised);
        assert_eq!(del.original_allele_string.as_deref(), Some("AC/A"));
        assert_eq!(del.uploaded_allele_string.as_deref(), Some("AC/A"));
        assert_eq!(del.uploaded_variation(), "21_100_AC/A");
        assert_eq!(del.variation_name(), "21_101_C/-");
        assert_eq!(del.location(), "21:101");
        let ins = &v[1];
        assert_eq!((ins.start, ins.end), (101, 100));
        assert_eq!(ins.allele_string, "-/T");
        assert_eq!(ins.uploaded_variation(), "21_100_A/AT");
        for multi in &v[2..] {
            assert_eq!((multi.start, multi.end), (100, 100));
            assert_eq!(multi.record_allele_string_multi.as_deref(), Some("A/G/AT"));
            assert_eq!(multi.original_allele_string, None);
        }
        assert_eq!(v[3].allele_string, "A/AT");
    }

    /// Parser.pm:554 upper-cases the allele string before the alleles are read,
    /// so `a/t` is the SNV `A/T` with Allele `T` and the uploaded allele column
    /// keeps the written case. An ID-less record is named at Parser.pm:462 from
    /// the written alleles as they stand after the minimisation, and the name
    /// differs from the upper-cased one, so Ensembl VEP 116.2 prints
    /// `21_30348410_ta/-` for `21 30348409 30348411 cta/c` and `21_30282169_c/a`
    /// for `21 30282169 30282169 c/a` in `Uploaded_variation`; a lower-case
    /// record with an ID is named by it.
    #[test]
    fn lower_case_alleles_are_upper_cased() {
        let v = all("21\t20006769\t20006769\ta/t\t1\n\
             21\t30348409\t30348411\tcta/c\t-\n\
             21\t100\t101\tac/a\t1\trs1\n\
             21\t100\t100\tA/g/T\t1\t0\n");
        assert_eq!(v[0].allele_string, "A/T");
        assert_eq!(v[0].ref_allele, b"A");
        assert_eq!(v[0].display_allele(), "T");
        assert_eq!(v[0].class_so_term(), "SNV");
        assert_eq!(v[0].uploaded_variation(), "21_20006769_a/t");
        assert_eq!(v[0].variation_name(), "21_20006769_a/t");
        assert_eq!(v[0].uploaded_allele_string.as_deref(), Some("a/t"));
        assert_eq!(v[1].allele_string, "TA/-");
        assert_eq!((v[1].start, v[1].end), (30348410, 30348411));
        assert_eq!(v[1].original_allele_string.as_deref(), Some("cta/c"));
        assert_eq!(v[1].uploaded_allele_string.as_deref(), Some("cta/c"));
        assert_eq!(v[1].uploaded_variation(), "21_30348410_ta/-");
        assert_eq!(v[1].variation_name(), "21_30348410_ta/-");
        assert_eq!(v[2].allele_string, "C/-");
        assert_eq!(v[2].uploaded_variation(), "rs1");
        assert_eq!(v[2].variation_name(), "rs1");
        assert_eq!(v[3].uploaded_variation(), "21_100_A/g/T");
        assert_eq!(v[4].variation_name(), "21_100_A/g/T");
        assert_eq!(v[4].record_allele_string_multi.as_deref(), Some("A/G/T"));
        // The bases are compared as written before the upper-casing, so a
        // mixed-case pair shares none and stays over its three bases.
        let mixed = all("21\t100\t101\tAc/a\t1\tmixed\n");
        assert_eq!(mixed[0].allele_string, "AC/A");
        assert_eq!((mixed[0].start, mixed[0].end), (100, 101));
        assert!(!mixed[0].minimised);
    }

    /// The string Parser.pm:554 upper-cases is kept as written beside the
    /// upper-cased alleles when a base is lower case, the minimised pair of a
    /// bi-allelic record and the whole column of a multi-allelic one, so a
    /// record a check ahead of that line refused and `--dont_skip` kept can be
    /// annotated as written; an upper-case record carries none.
    #[test]
    fn lower_case_alleles_keep_their_written_string() {
        let v = all("21\t20006769\t20006769\ta/t\t1\n\
             21\t30348409\t30348411\tcta/c\t-\n\
             21\t100\t100\tA/g/T\t1\t0\n\
             21\t100\t100\tA/G\t1\tupper\n");
        assert_eq!(v[0].written_case_allele_string.as_deref(), Some("a/t"));
        assert_eq!(v[1].written_case_allele_string.as_deref(), Some("ta/-"));
        assert_eq!(v[2].written_case_allele_string.as_deref(), Some("A/g/T"));
        assert!(
            Arc::ptr_eq(
                v[2].written_case_allele_string.as_ref().unwrap(),
                v[3].written_case_allele_string.as_ref().unwrap()
            ),
            "the alleles of one record share the column"
        );
        assert_eq!(v[4].written_case_allele_string, None);
    }

    /// Under `--minimal` the bi-allelic trim runs on the column as written
    /// (Parser.pm:852 before :554), so `AAAG/AGGg` loses its shared `A` and keeps
    /// the `G`/`g` pair, is named from that pair at the trimmed start (:462) and
    /// is upper-cased afterwards; the upper-case column trims on both sides. A
    /// multi-allelic column is left to the split, as without the flag.
    #[test]
    fn minimal_trims_the_written_case_before_the_upper_casing() {
        let context = ParserContext {
            minimal: true,
            ..ParserContext::default()
        };
        let data = "21\t10542414\t10542417\tAAAG/AGGg\t+\n\
                    21\t10542414\t10542417\tAAAG/AGGG\t+\n\
                    21\t100\t100\tA/g/T\t1\n";
        let mut parser = EnsemblParser::new(Cursor::new(data.as_bytes().to_vec()), &context);
        let written = parser.next_variant().unwrap().unwrap();
        assert_eq!((written.start, written.end), (10542415, 10542417));
        assert_eq!(written.allele_string, "AAG/GGG");
        assert_eq!(written.id.as_deref(), Some("21_10542415_AAG/GGg"));
        assert_eq!(
            written.written_case_allele_string.as_deref(),
            Some("AAG/GGg")
        );
        assert_eq!(written.uploaded_allele_string.as_deref(), Some("AAAG/AGGg"));
        assert!(written.minimised);
        let upper = parser.next_variant().unwrap().unwrap();
        assert_eq!((upper.start, upper.end), (10542415, 10542416));
        assert_eq!(upper.allele_string, "AA/GG");
        assert_eq!(upper.id, None);
        let multi = parser.next_variant().unwrap().unwrap();
        assert_eq!(multi.allele_string, "A/G");
        assert_eq!(multi.record_allele_string_multi.as_deref(), Some("A/G/T"));
    }

    /// The columns are split on the context's delimiter (Parser/VEP_input.pm:90
    /// hands the configured one to ensembl-io): a tab inside a line read under
    /// the whitespace pattern still separates, a leading run of spaces under the
    /// space-run pattern is an empty first column, a space-delimited line read
    /// under a tab is one column and fails, and a comma-delimited one splits on
    /// its literal.
    #[test]
    fn columns_split_on_the_context_delimiter() {
        let read = |delimiter: Delimiter, data: &str| {
            parser_under(delimiter, data).next_variant().unwrap()
        };
        let v = read(
            Delimiter::Whitespace,
            "21 25585733  25585733\tA/G 1 rs699\n",
        )
        .unwrap();
        assert_eq!((v.start, v.id.as_deref()), (25585733, Some("rs699")));
        assert_eq!(
            v.raw_input.as_deref(),
            Some("21 25585733 25585733 A/G 1 rs699")
        );
        // A leading run of spaces is a leading empty field under `" +"`, so the
        // columns shift: an empty chromosome, then `21` as the start.
        let v = read(Delimiter::SpaceRun, " 21  100   100 A/G\n").unwrap();
        assert_eq!((v.chr.as_str(), v.start, v.end), ("", 21, 100));
        let v = read(Delimiter::Tab, "21 25585733 25585733 A/G 1\n");
        assert!(v.unwrap_err().to_string().contains("at least 4 columns"));
        let v = read(Delimiter::Literal(",".to_string()), "21,100,100,A/G,-1,x\n").unwrap();
        assert_eq!((v.start, v.strand), (100, Strand::Reverse));
        assert_eq!(v.raw_input.as_deref(), Some("21,100,100,A/G,-1,x"));
    }

    /// Each field loses its trailing whitespace (`map {$_ =~ s/\s+$//}`), so a
    /// `\r` or trailing blanks on the last column do not reach the identifier,
    /// and the `\r` of a `\r\n` line end is removed with it.
    #[test]
    fn fields_lose_their_trailing_whitespace() {
        let v = all("21\t100\t100\tA/G\t1\trs1 \r\n21\t200\t200\tA/G \t-1\trs2\r");
        assert_eq!(v[0].id.as_deref(), Some("rs1"));
        assert_eq!(v[1].allele_string, "A/G");
        assert_eq!(v[1].id.as_deref(), Some("rs2"));
        assert_eq!(v[1].strand, Strand::Reverse);
        // The re-emitted line keeps the fields as split (`_line` is the record,
        // not the trimmed copies): a blank inside a field stays, as does the `\r`
        // of a last line with no `\n` behind it, which no `\r\n` rule removes.
        assert_eq!(
            v[0].raw_input.as_deref(),
            Some("21\t100\t100\tA/G\t1\trs1 ")
        );
        assert_eq!(
            v[1].raw_input.as_deref(),
            Some("21\t200\t200\tA/G \t-1\trs2\r")
        );
        // The line's name reads the same record (OutputFactory.pm:912-917), so an
        // ID-less record keeps the blank too.
        let unnamed = all("21\t200\t200\tA/G \t-1\n");
        assert_eq!(unnamed[0].uploaded_variation(), "21_200_A/G ");
    }

    /// A multi-allelic line keeps the whole allele column as the record's allele
    /// string on every allele (Parser/VEP_input.pm:156), upper-cased as `validate_vf`
    /// leaves it (Parser.pm:554); a bi-allelic line keeps none. Every allele also
    /// keeps the column as written, the string `UPLOADED_ALLELE` is read from.
    #[test]
    fn multi_allelic_line_carries_the_record_allele_string() {
        let data = "21\t38218616\t38218617\ttt/-/t\t1\tmixed\n21\t100\t100\tA/G\t1\tsnv\n";
        let mut parser = parser(data);
        let first = parser.next_variant().unwrap().unwrap();
        let second = parser.next_variant().unwrap().unwrap();
        assert_eq!(first.record_allele_string_multi.as_deref(), Some("TT/-/T"));
        assert_eq!(second.record_allele_string_multi.as_deref(), Some("TT/-/T"));
        assert_eq!(first.uploaded_allele_string.as_deref(), Some("tt/-/t"));
        assert_eq!(second.uploaded_allele_string.as_deref(), Some("tt/-/t"));
        assert_eq!(second.allele_index, 1);
        let snv = parser.next_variant().unwrap().unwrap();
        assert_eq!(snv.record_allele_string_multi, None);
        assert_eq!(snv.uploaded_allele_string.as_deref(), Some("A/G"));
    }

    #[test]
    fn test_parse_deletion() {
        let mut parser = parser("21\t25585656\t25585660\tACGTG/-\t1\t.\n");

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585656);
        assert_eq!(v.end, 25585660);
        assert_eq!(v.ref_allele, b"ACGTG");
        assert_eq!(v.alt_allele(), b"-");
        assert_eq!(v.id, Some(".".to_string()));

        assert!(parser.next_variant().is_none());
    }

    /// OutputFactory.pm `VariationFeature_to_output_hash`, the ensembl `_line`
    /// branch: an ID-less line is named `chr_start_alleles` with the lower of its
    /// two coordinates (99 for the insertion `100 99`) and the allele column as
    /// written (`A/G/T` on every allele of a multi-allelic line); an ID names it.
    #[test]
    fn id_less_line_is_named_by_its_lower_coordinate_and_allele_column() {
        let names: Vec<String> =
            all("21\t100\t99\t-/A\t1\n21\t200\t200\tA/G/T\t1\t.\n21\t300\t300\tC/T\t1\trs1\n")
                .iter()
                .map(InputVariant::uploaded_variation)
                .collect();
        assert_eq!(names, ["21_99_-/A", "21_200_A/G/T", "21_200_A/G/T", "rs1"]);
    }

    /// Parser/VEP_input.pm:127 takes the ID column whole, so `a;b` names the record
    /// `a;b` in `Uploaded_variation` and in the JSON `id` (a VCF ID field is split on
    /// `;`). A `.` stays the JSON `id` (Perl reads it as true) while the record is
    /// named by its line; `0` is replaced by `validate_vf`'s generated name, here
    /// `21_100_-/A` (start 100, allele string `-/A`), which OutputFactory.pm rejects,
    /// so the record is named by its line at the lower coordinate, `21_99_-/A`.
    #[test]
    fn an_ensembl_format_id_column_is_kept_whole() {
        let names: Vec<(String, String)> =
            all("21\t300\t300\tC/T\t1\ta;b\n21\t200\t200\tA/G\t1\t.\n21\t100\t99\t-/A\t1\t0\n")
                .iter()
                .map(|v| (v.uploaded_variation(), v.variation_name()))
                .collect();
        assert_eq!(
            names,
            [
                ("a;b".to_string(), "a;b".to_string()),
                ("21_200_A/G".to_string(), ".".to_string()),
                ("21_99_-/A".to_string(), "21_100_-/A".to_string()),
            ]
        );
    }

    #[test]
    fn test_parse_multi_allelic() {
        let mut parser = parser("21\t25585733\t25585733\tA/G/T\t1\trs699\n");

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.alt_allele(), b"G");
        assert_eq!(v1.allele_index, 0);
        assert!(v1.minimised);

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.alt_allele(), b"T");
        assert_eq!(v2.allele_index, 1);
        assert!(v2.minimised);

        assert!(parser.next_variant().is_none());
    }

    /// A coordinate is read up to `MAX_COORD`; a start or end past it fails the line.
    #[test]
    fn coordinates_past_the_ceiling_fail_the_line() {
        let past = MAX_COORD + 1;
        let data = format!(
            "21\t{MAX_COORD}\t{MAX_COORD}\tA/G\t1\n\
             21\t{past}\t{past}\tA/G\t1\n\
             21\t100\t{}\tA/G\t1\n",
            i64::MAX
        );
        let mut parser = parser(&data);
        let kept = parser.next_variant().unwrap().unwrap();
        assert_eq!((kept.start, kept.end), (MAX_COORD, MAX_COORD));
        let start_err = parser.next_variant().unwrap().unwrap_err().to_string();
        assert!(start_err.contains("Invalid start position"), "{start_err}");
        let end_err = parser.next_variant().unwrap().unwrap_err().to_string();
        assert!(end_err.contains("Invalid end position"), "{end_err}");
    }

    #[test]
    fn test_missing_identifier() {
        let v = all("21\t25585733\t25585733\tA/G\t1\n");
        assert_eq!(v[0].chr, "21");
        assert_eq!(v[0].id, None);
    }

    /// Every line counts toward Ensembl VEP's line number (Parser/VEP_input.pm:119
    /// per record read, a `#` line and an empty line included, as ensembl-io's
    /// VEP_input parser has no metadata lines, VEP_input.pm:61): a `#` line is a
    /// comment here, counted and passed over without a warning (Ensembl VEP
    /// 116.2 numbers `#chr start end allele strand id`, a record, an oversize
    /// record as `line 1 skipped` twice and `line 3 skipped`); an empty or
    /// whitespace-only line, which Ensembl VEP's own parser dies on
    /// (Parser/VEP_input.pm:136), is counted and warned about here in
    /// `skip_empty_lines`'s words (Parser.pm:273-290) when its run ends, at the
    /// next record or the end of the input; a multi-allelic line's alleles share
    /// its number.
    #[test]
    fn every_line_is_counted_and_empty_runs_warned() {
        let mut parser =
            parser("# comment\n\n  \t\n21\t25585733\t25585733\tA/G/T\t1\trs699\n# c\n\n");
        assert!(parser.take_warnings().is_empty());
        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!((v.chr.as_str(), v.allele_index), ("21", 0));
        assert_eq!(parser.record_line_number(), Some(4));
        assert_eq!(
            parser.take_warnings(),
            ["Skipped 2 empty lines from line 2"]
        );
        let second = parser.next_variant().unwrap().unwrap();
        assert_eq!(second.allele_index, 1);
        assert_eq!(parser.record_line_number(), Some(4));
        assert!(parser.take_warnings().is_empty());
        assert!(parser.next_variant().is_none());
        assert_eq!(parser.take_warnings(), ["Skipped 1 empty line from line 6"]);
        assert_eq!(parser.take_empty_lines(), 0);
    }

    #[test]
    fn test_multiple_lines() {
        let v = all(
            "21\t25585733\t25585733\tA/G\t1\trs699\n21\t25587758\t25587758\tG/A\t1\trs4149056\n",
        );
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].start, 25585733);
        assert_eq!(v[0].id, Some("rs699".to_string()));
        assert_eq!(v[1].start, 25587758);
        assert_eq!(v[1].id, Some("rs4149056".to_string()));
    }
}
