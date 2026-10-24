// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Region format input, read as Ensembl VEP's `Parser/Region.pm`
//! `create_VariationFeatures` (release/116.2, lines 157-222) reads it over
//! ensembl-io's list parser:
//!
//! ```text
//! 21:25585733-25585733:1/G
//! 21:25585656-25585660:-1/ACGTG
//! X:100000-100000/T
//! 21:25585700-25585900/DUP
//! ```
//!
//! One region per line, matched whole against `_valid_region_regex`
//! (ensembl-variation release/116, Utils/VEP.pm:391),
//! `^([^:]+):(\d+)-(\d+)(:[-\+]?1)?[\/:]([a-z0-9:]{3,}|[ACGTN-]+)$` without
//! regard to case: the chromosome is the text before the first colon, the
//! start and the end are both written, the strand is `:1`, `:+1` or `:-1`, and
//! the allele follows a `/` or a `:`. The record is reverse when the strand
//! token contains a `-` and forward otherwise (`:183`), and the allele is
//! upper-cased (`:176`). An allele of `ACGTN-` names a sequence variant whose
//! reference allele is `N` per base of the span, or `-` when the start lies
//! past the end (`:188-189`); the format switches `--lookup_ref` on
//! (`:121-122`), so the run replaces that placeholder with the reference
//! sequence over the span, which is why the format needs the run's FASTA
//! (`:115`) and refuses `--check_ref` (`:118`). Any other allele names a
//! structural variant class through `get_SO_term` (`:205`, Parser.pm:663-728).
//! The record is named by its line as written (`:198, 213`), which
//! `Uploaded_variation`, the JSON `id` and the JSON `input` print, and a
//! sequence record goes through the minimisation every parsed record goes
//! through ([`minimise_alleles`]).
//!
//! A line the expression does not match, or whose structural class `%SO_TERMS`
//! lacks, builds no feature in Ensembl VEP (`:171, 206`), which ends the buffer
//! fill at that line (InputBuffer.pm:212) and the run itself when the line
//! opens a fill (Runner.pm:424), the lines after it unread; here such a line
//! fails with [`IoError::SkippedLine`], carrying `get_SO_term`'s warning for
//! the class (Parser.pm:725) or this parser's own text, for the caller to warn
//! about and skip, and the lines after it are read. A span over
//! `--max_sv_size` is refused with the text `validate_vf` refuses it with
//! (Parser.pm:491-500) before its run of `N` is built. The lines are numbered
//! as Ensembl VEP numbers them (`LineCount`): a line with no word character is
//! an empty line, warned about in `skip_empty_lines`'s words and counted
//! (Parser.pm:273), and any other `#` line is a comment here, counted toward
//! the line numbers that follow and otherwise passed over without a word,
//! where ensembl-io's list parser reads it as a record (ListBasedParser.pm:60)
//! that the expression refuses, with the consequences above.

use std::io::BufRead;

use vep_core::coordinate::{parse_coordinate, Strand, MAX_COORD};
use vep_core::sv_routing::{sv_so_term_and_class, sv_type_abbreviation};
use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::ensembl::structural_variant;
use super::{
    has_word_character, is_perl_space, minimise_alleles, InputParser, LineCount, ParserContext,
};

/// Parser for region format input, one region per line.
pub struct RegionParser<R: BufRead> {
    reader: R,
    /// `--max_sv_size`, the span a record is refused beyond; `None` for the
    /// flag's `-1`.
    max_sv_size: Option<u64>,
    /// `--minimal` ([`ParserContext::minimal`]), under which the span is
    /// measured on every sequence record's trimmed alleles.
    minimal: bool,
    line_buf: Vec<u8>,
    /// Ensembl VEP's line number over the lines read, and the empty-line
    /// warnings it raises.
    lines: LineCount,
}

/// A line's parts as `_valid_region_regex` captures them, each as written.
pub(crate) struct Region<'a> {
    pub chr: &'a str,
    /// The start's digits.
    pub start: &'a str,
    /// The end's digits.
    pub end: &'a str,
    /// The strand token with its colon (`:1`, `:+1`, `:-1`), when written.
    pub strand: Option<&'a str>,
    pub allele: &'a str,
}

impl<R: BufRead> RegionParser<R> {
    pub fn new(reader: R, context: &ParserContext) -> Self {
        Self {
            reader,
            max_sv_size: context.max_sv_size,
            minimal: context.minimal,
            line_buf: Vec::new(),
            lines: LineCount::default(),
        }
    }

    /// Parse one line, its line end removed, into its variant.
    fn parse_record(&self, line: &str) -> Result<InputVariant, IoError> {
        // A refusal quotes the line collapsed, not as written: a `\r` the line
        // end left (`\r\r\n`), or a tab inside a chromosome name the expression
        // admits, would carry a control byte into the warning.
        let region = match_region(line).ok_or_else(|| {
            IoError::SkippedLine(format!(
                "Unable to parse region '{}': expected chr:start-end[:strand]/allele",
                collapse_space(line)
            ))
        })?;
        let coordinate = |digits: &str| {
            parse_coordinate(digits).ok_or_else(|| {
                IoError::SkippedLine(format!(
                    "Unable to parse region '{}': a coordinate is over {MAX_COORD}",
                    collapse_space(line)
                ))
            })
        };
        let start = coordinate(region.start)?;
        let end = coordinate(region.end)?;
        let allele = region.allele.to_ascii_uppercase();
        let sequence = is_sequence_allele(&allele);
        if !sequence {
            // `get_SO_term` refuses an unknown class as the feature is built
            // (Parser/Region.pm:205-206), before `validate_vf` measures the span.
            let abbreviation = sv_type_abbreviation(&allele, None);
            if sv_so_term_and_class(&abbreviation).is_none() {
                return Err(IoError::SkippedLine(format!(
                    "{abbreviation} is not a supported structural variant type"
                )));
            }
        }
        // `validate_vf` measures end minus start (Parser.pm:492) on the record
        // `post_process_vfs` leaves, so the span is judged before the run of `N`
        // is built and a span over the limit costs nothing.
        let span = if sequence {
            minimised_span(start, end, allele.as_bytes(), self.minimal)
        } else {
            end.saturating_sub(start)
        };
        if let Some(max) = self.max_sv_size.filter(|&max| span > max) {
            return Err(IoError::SkippedLine(format!(
                "variant size ({span}) is bigger than --max_sv_size ({max})"
            )));
        }
        let chr = region.chr.to_string();
        let mut variant = if sequence {
            let reference = if start <= end {
                "N".repeat((end - start + 1) as usize)
            } else {
                "-".to_string()
            };
            let mut variant =
                InputVariant::new(chr, start, end, reference.into_bytes(), allele.into_bytes());
            // `post_process_vfs` runs `minimise_alleles` before the lookup replaces
            // the placeholder, on every bi-allelic record under `--minimal`
            // (Parser.pm:852) and otherwise on a record whose alleles differ in
            // length and carry no `-` (:862-881), and that keeps the string as
            // parsed in `original_allele_string` (:1008) whether or not a base is
            // trimmed, which `UPLOADED_ALLELE` prints (OutputFactory.pm:1360); the
            // placeholder is kept for every record and `uploaded_allele` makes
            // that choice under the run's `--minimal` setting.
            variant.uploaded_allele_string = Some(variant.allele_string.as_str().into());
            minimise_alleles(&mut variant);
            variant
        } else {
            structural_variant(chr, start, end, &allele)
        };
        variant.strand = if region.strand.is_some_and(|token| token.contains('-')) {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        variant.id = Some(line.to_string());
        variant.raw_input = Some(line.into());
        Ok(variant)
    }
}

/// The captures of `_valid_region_regex` on `line`, `None` when the expression
/// does not match it. The optional strand group is taken whenever an allele
/// follows it, as the expression's engine keeps the group and gives it up only
/// when the rest of the line cannot match without it. Format detection tests a
/// line against the same expression through this ([`super::detect`]).
pub(crate) fn match_region(line: &str) -> Option<Region<'_>> {
    let (chr, rest) = line.split_once(':')?;
    if chr.is_empty() {
        return None;
    }
    let (start, rest) = leading_digits(rest)?;
    let (end, rest) = leading_digits(rest.strip_prefix('-')?)?;
    let with_strand = strand_token(rest)
        .and_then(|(strand, after)| allele_after(after).map(|allele| (Some(strand), allele)));
    let (strand, allele) = match with_strand {
        Some(parts) => parts,
        None => (None, allele_after(rest)?),
    };
    Some(Region {
        chr,
        start,
        end,
        strand,
        allele,
    })
}

/// The span `validate_vf` measures for a sequence record, `end - start`
/// (Parser.pm:492), on the record `post_process_vfs` leaves: a record whose
/// alleles differ in length and carry no `-` (:862-881), and under `minimal`
/// (`--minimal`) every record (:852), goes through `minimise_alleles`, whose
/// `trim_sequences` removes the bases the two alleles share from the left and
/// then from the right, moving the start and the end in by as many bases. The
/// reference is `N` per base, so the only bases it can share with the allele
/// are the allele's own leading and trailing `N`s, and the span is computed
/// from them without building the run. An insertion (start past end) has no
/// run and no span to measure.
fn minimised_span(start: u64, end: u64, allele: &[u8], minimal: bool) -> u64 {
    if start > end {
        return 0;
    }
    let span = end - start;
    let reference_len = span + 1;
    if !minimal && (allele.len() as u64 == reference_len || allele.contains(&b'-')) {
        return span;
    }
    let leading = allele.iter().take_while(|&&b| b == b'N').count();
    let prefix = leading.min(reference_len.min(allele.len() as u64) as usize);
    let rest = &allele[prefix..];
    let trailing = rest.iter().rev().take_while(|&&b| b == b'N').count();
    let suffix = trailing.min((reference_len - prefix as u64).min(rest.len() as u64) as usize);
    span.saturating_sub((prefix + suffix) as u64)
}

/// `line` with every run of white space (Perl's `\s`) one space, as a skipped
/// line is quoted (`s/\s+/ /g`, BaseVEP.pm:849).
fn collapse_space(line: &str) -> String {
    let mut collapsed = String::with_capacity(line.len());
    let mut in_space = false;
    for c in line.chars() {
        if is_perl_space(c) {
            if !in_space {
                collapsed.push(' ');
            }
            in_space = true;
        } else {
            collapsed.push(c);
            in_space = false;
        }
    }
    collapsed
}

/// `(\d+)` opening `text`: the digits and the text after them.
fn leading_digits(text: &str) -> Option<(&str, &str)> {
    let digits = text.bytes().take_while(u8::is_ascii_digit).count();
    (digits > 0).then(|| text.split_at(digits))
}

/// `(:[-\+]?1)` opening `text`: the token and the text after it.
fn strand_token(text: &str) -> Option<(&str, &str)> {
    let after_colon = text.strip_prefix(':')?;
    let unsigned = after_colon.strip_prefix(['-', '+']).unwrap_or(after_colon);
    let after = unsigned.strip_prefix('1')?;
    Some(text.split_at(text.len() - after.len()))
}

/// `[\/:]([a-z0-9:]{3,}|[ACGTN-]+)$` on `text`: the allele after its separator.
fn allele_after(text: &str) -> Option<&str> {
    text.strip_prefix(['/', ':'])
        .filter(|allele| is_allele(allele))
}

/// `([a-z0-9:]{3,}|[ACGTN-]+)` without regard to case, on the whole of `text`.
fn is_allele(text: &str) -> bool {
    let word = text.len() >= 3 && text.bytes().all(|b| b.is_ascii_alphanumeric() || b == b':');
    word || is_sequence_allele(text)
}

/// `^[ACGTN-]+$` without regard to case (`:188`).
fn is_sequence_allele(text: &str) -> bool {
    !text.is_empty()
        && text.bytes().all(|b| {
            matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'N' | b'-'
            )
        })
}

impl<R: BufRead + Send> InputParser for RegionParser<R> {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
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
            // ensembl-io's `TextParser` turns `\r\n` into `\n` (:138) and `chomp`
            // removes it; any other trailing blank stays and fails the expression.
            let line = match text.strip_suffix('\n') {
                Some(without_newline) => without_newline
                    .strip_suffix('\r')
                    .unwrap_or(without_newline),
                None => &text,
            };
            // `skip_empty_lines` reads a `#` line with no word character past as
            // an empty line (Parser.pm:273), before the parser sees it.
            if !has_word_character(line) {
                self.lines.empty_line();
                continue;
            }
            self.lines.record_line();
            if line.starts_with('#') {
                continue;
            }
            return Some(self.parse_record(line));
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
    use vep_core::variant::VariantClass;

    fn parser_in(context: &ParserContext, data: &str) -> RegionParser<Cursor<Vec<u8>>> {
        RegionParser::new(Cursor::new(data.as_bytes().to_vec()), context)
    }

    fn parser(data: &str) -> RegionParser<Cursor<Vec<u8>>> {
        parser_in(&ParserContext::default(), data)
    }

    /// What the parser yields for `region` as the input's one line.
    fn record(region: &str) -> Result<InputVariant, IoError> {
        parser(&format!("{region}\n")).next_variant().unwrap()
    }

    fn one(region: &str) -> InputVariant {
        record(region).unwrap_or_else(|e| panic!("{region}: {e}"))
    }

    fn refused(region: &str) -> String {
        let error = record(region).unwrap_err();
        assert!(
            matches!(error, IoError::SkippedLine(_)),
            "{region}: {error:?}"
        );
        error.to_string()
    }

    /// Parser/Region.pm:189: a sequence record's reference allele is `N` per base
    /// of the span, the placeholder the run's reference lookup replaces, so
    /// `21:100-100/G` is `N/G` at 100 and `21:100-102/A` is `NNN/A` over 100 to
    /// 102; the record is forward without a strand token and named by its line,
    /// which its re-emitted input is too, and the placeholder is kept as the
    /// string uploaded.
    #[test]
    fn ref_is_n_per_base_of_the_span() {
        let v = one("21:100-100/G");
        assert_eq!((v.start, v.end), (100, 100));
        assert_eq!(v.allele_string, "N/G");
        assert_eq!(v.ref_allele, b"N");
        assert_eq!(v.alt_allele(), b"G");
        assert_eq!(v.strand, Strand::Forward);
        assert!(!v.is_structural);
        assert_eq!(v.uploaded_variation(), "21:100-100/G");
        assert_eq!(v.variation_name(), "21:100-100/G");
        assert_eq!(v.raw_input.as_deref(), Some("21:100-100/G"));
        assert_eq!(v.uploaded_allele_string.as_deref(), Some("N/G"));
        let v = one("21:100-102/A");
        assert_eq!((v.start, v.end), (100, 102));
        assert_eq!(v.allele_string, "NNN/A");
        let v = one("21:100-104/-");
        assert_eq!(v.allele_string, "NNNNN/-");
        assert_eq!(v.variant_class, VariantClass::Deletion);
    }

    /// The allele follows a `/` or a `:` (`[\/:]` in the expression), with or
    /// without a strand token before it, so `21:100-100:G`, `21:100-100:1:G` and
    /// `21:100-100:-1:G` read as their `/` forms.
    #[test]
    fn colon_allele_separator() {
        for (region, strand) in [
            ("21:100-100:G", Strand::Forward),
            ("21:100-100:1:G", Strand::Forward),
            ("21:100-100:-1:G", Strand::Reverse),
            ("21:100-100:+1:G", Strand::Forward),
        ] {
            let v = one(region);
            assert_eq!((v.start, v.end), (100, 100), "{region}");
            assert_eq!(v.allele_string, "N/G", "{region}");
            assert_eq!(v.strand, strand, "{region}");
            assert_eq!(v.uploaded_variation(), region);
        }
        // The strand group is taken before the separator is read, so `:1:ACG` is a
        // forward-strand `ACG`, not the structural word `1:ACG`.
        let v = one("21:100-102:1:ACG");
        assert_eq!(v.allele_string, "NNN/ACG");
        assert!(!v.is_structural);
    }

    /// An insertion is written with its start one past its end
    /// (Parser/Region.pm:62), where `'N' x 0` leaves the reference `-`
    /// (`:189`): `21:101-100/G` is `-/G` with start 101 and end 100. A start
    /// further past the end gets the same `-` and is left to `validate_vf`'s
    /// coordinate check.
    #[test]
    fn insertion_start_end_plus_one() {
        let v = one("21:101-100/G");
        assert_eq!((v.start, v.end), (101, 100));
        assert_eq!(v.allele_string, "-/G");
        assert_eq!(v.variant_class, VariantClass::Insertion);
        assert_eq!(v.location(), "21:100-101");
        let v = one("21:105-100/G");
        assert_eq!((v.start, v.end), (105, 100));
        assert_eq!(v.allele_string, "-/G");
    }

    /// An allele outside `ACGTN-` names a structural class through `get_SO_term`
    /// on the upper-cased word (Parser/Region.pm:176, 205): `DUP` and `del` are
    /// the `%SO_TERMS` duplication and deletion over the line's start and end,
    /// `INS:ME:ALU` an Alu insertion, each named by its line with the class term
    /// standing in for the allele string. A word the table lacks is refused with
    /// `get_SO_term`'s warning (Parser.pm:725), as Ensembl VEP builds no feature
    /// for it (`:206`).
    #[test]
    fn sv_word_dup_del() {
        let dup = one("21:25585700-25585900/DUP");
        assert!(dup.is_structural);
        assert_eq!((dup.start, dup.end), (25585700, 25585900));
        assert_eq!(dup.sv_end, Some(25585900));
        assert_eq!(dup.variant_class, VariantClass::Duplication);
        assert_eq!(dup.record_so_term.as_deref(), Some("duplication"));
        assert_eq!(dup.display_allele(), "duplication");
        assert_eq!(dup.uploaded_variation(), "21:25585700-25585900/DUP");
        assert_eq!(dup.variation_name(), "21:25585700-25585900/DUP");
        assert_eq!(dup.strand, Strand::Forward);
        assert!(!dup.vep_skip);
        let del = one("21:25585700-25585900:-1/del");
        assert_eq!(del.record_so_term.as_deref(), Some("deletion"));
        assert_eq!(del.variant_class, VariantClass::StructuralDeletion);
        assert_eq!(del.strand, Strand::Reverse);
        assert_eq!(del.uploaded_variation(), "21:25585700-25585900:-1/del");
        let alu = one("21:25585700-25585900/INS:ME:ALU");
        assert_eq!(alu.record_so_term.as_deref(), Some("Alu_insertion"));
        assert_eq!(
            refused("21:25585700-25585900/FOO"),
            "FOO is not a supported structural variant type"
        );
        assert_eq!(
            refused("21:25585700-25585900/foo:bar"),
            "FOO:BAR is not a supported structural variant type"
        );
    }

    /// Parser/Region.pm:183: the record is reverse when the strand token contains
    /// a `-` and forward for `:1`, `:+1` or no token; the alleles stay as written
    /// on the record's strand for the run to read against the reference and the
    /// features on the other strand.
    #[test]
    fn strand_tokens_one_plus_one_minus_one() {
        for (region, strand) in [
            ("21:25585656-25585660:-1/ACGTG", Strand::Reverse),
            ("21:25585656-25585660:1/ACGTG", Strand::Forward),
            ("21:25585656-25585660:+1/ACGTG", Strand::Forward),
            ("21:25585656-25585660/ACGTG", Strand::Forward),
        ] {
            let v = one(region);
            assert_eq!(v.strand, strand, "{region}");
            assert_eq!(v.allele_string, "NNNNN/ACGTG", "{region}");
            assert_eq!((v.start, v.end), (25585656, 25585660), "{region}");
        }
    }

    /// The expression is read without regard to case (`/i`): the allele is
    /// upper-cased (`:176`) while the record's name and re-emitted line keep the
    /// line as written.
    #[test]
    fn allele_upper_cased_and_line_kept_in_the_name() {
        let v = one("chr21:100-101/ac");
        assert_eq!(v.allele_string, "NN/AC");
        assert_eq!(v.original_chr, "chr21");
        assert_eq!(v.chr, "21");
        assert_eq!(v.uploaded_variation(), "chr21:100-101/ac");
        assert_eq!(v.raw_input.as_deref(), Some("chr21:100-101/ac"));
        assert_eq!(v.id.as_deref(), Some("chr21:100-101/ac"));
        assert_eq!(one("21:100-100/n").allele_string, "N/N");
    }

    /// `post_process_vfs` minimises a bi-allelic record whose alleles differ in
    /// length and carry no `-` (Parser.pm:862-881), keeping the string as parsed
    /// for `UPLOADED_ALLELE` whether or not a base is trimmed: `21:100-101/N` is
    /// `NN/N`, trimmed to `N/-` at 101; `21:100-102/A` shares no base with its
    /// `NNN` and stays, its placeholder kept all the same (Ensembl VEP 116.2
    /// prints `UPLOADED_ALLELE=NN/GTA` for `21:10552677-10552678/GTA`). A
    /// record of one length on each side, or with a `-`, is not minimised and
    /// keeps its placeholder too, for `--minimal`, which minimises every
    /// bi-allelic record (:852) and prints it (Ensembl VEP 116.2 prints
    /// `UPLOADED_ALLELE=N/A` for `21:14221107-14221107/A` under `--minimal` and
    /// the looked-up `C/A` without it).
    #[test]
    fn shared_bases_are_minimised() {
        let v = one("21:100-101/N");
        assert_eq!((v.start, v.end), (101, 101));
        assert_eq!(v.allele_string, "N/-");
        assert!(v.minimised);
        assert_eq!(v.original_allele_string.as_deref(), Some("NN/N"));
        assert_eq!(v.uploaded_allele_string.as_deref(), Some("NN/N"));
        assert_eq!(v.uploaded_variation(), "21:100-101/N");
        let v = one("21:100-102/A");
        assert!(!v.minimised);
        assert_eq!(v.allele_string, "NNN/A");
        assert_eq!(v.uploaded_allele_string.as_deref(), Some("NNN/A"));
        assert_eq!(v.uploaded_allele(false).as_deref(), Some("NNN/A"));
        for (region, placeholder) in [
            ("21:100-101/GG", "NN/GG"),
            ("21:100-102/-", "NNN/-"),
            ("21:101-100/GTA", "-/GTA"),
        ] {
            let v = one(region);
            assert!(!v.minimised, "{region}");
            assert_eq!(v.allele_string, placeholder, "{region}");
            assert_eq!(
                v.uploaded_allele_string.as_deref(),
                Some(placeholder),
                "{region}"
            );
            assert_eq!(
                v.uploaded_allele(true).as_deref(),
                Some(placeholder),
                "{region}"
            );
        }
    }

    /// Every shape outside `_valid_region_regex` is refused: a position without
    /// its end, a missing allele, a multi-allelic string, an allele outside the
    /// two forms, a strand other than `1`, `+1` or `-1`, an empty chromosome, a
    /// trailing blank and a VCF-shaped line. The warning names the line, each
    /// tab of the VCF-shaped one a space.
    #[test]
    fn lines_outside_the_expression_are_refused() {
        for region in [
            "21:100/A",
            "21:100-100",
            "21:100-100/",
            "21:100-100/A/G",
            "21:100-100/XY",
            "21:100-100:11/G",
            "21:100-100:2/G",
            "21:100-100::G",
            ":100-100/G",
            "21:100-100/G ",
            "21:100-100 /G",
            "21:a-100/G",
            "21:100-100-1/G",
            "21\t100\t100\tA/G",
            "21:g.100A>G",
        ] {
            let quoted = region.replace('\t', " ");
            assert_eq!(
                refused(region),
                format!(
                    "Unable to parse region '{quoted}': expected chr:start-end[:strand]/allele"
                )
            );
        }
    }

    /// A refused line is quoted with every run of white space one space, as the
    /// skip warning quotes it (BaseVEP.pm:849): a line ending `\r\r\n` keeps one
    /// `\r` once its line end is removed (ensembl-io TextParser.pm:138 replaces
    /// one `\r\n`, then `chomp`), quoted as a trailing space, and a run of tabs
    /// inside the line is one space; no control byte reaches the warning. The
    /// coordinate refusal quotes the line the same way: `[^:]+` admits a tab in
    /// the chromosome, so a line with one can match the expression and fail there.
    #[test]
    fn a_refused_line_is_quoted_with_its_white_space_collapsed() {
        let past = MAX_COORD + 1;
        let mut parser = parser(&format!(
            "21:100-100/G\r\r\n21:100-100\t\t/G\n2\t1:{past}-{past}/G\n"
        ));
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "Unable to parse region '21:100-100/G ': expected chr:start-end[:strand]/allele"
        );
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "Unable to parse region '21:100-100 /G': expected chr:start-end[:strand]/allele"
        );
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            format!(
                "Unable to parse region '2 1:{past}-{past}/G': a coordinate is over {MAX_COORD}"
            )
        );
        assert!(parser.next_variant().is_none());
    }

    /// `validate_vf`'s size test (Parser.pm:491-500) runs before the run of `N`
    /// is built, so an oversize span, one no memory could hold included, costs
    /// nothing and is refused with Ensembl VEP's text, a structural word too,
    /// and the lines after it are still read; an insertion's span is below
    /// zero and passes; no limit is applied for the flag's `-1`. A class word
    /// the table lacks is refused first, as `get_SO_term` refuses it while the
    /// feature is built (Parser/Region.pm:205-206), whatever its span. The span
    /// is the minimised record's: `21:1-1002/N` shares its first `N` with the
    /// run, so `minimise_alleles` leaves `N` x 1001 over 2 to 1002 and the span
    /// 1000 passes where `21:1-1002/G` is refused.
    #[test]
    fn oversize_span_refused_before_the_reference_is_built() {
        let context = ParserContext {
            max_sv_size: Some(1000),
            ..ParserContext::default()
        };
        let mut parser = parser_in(
            &context,
            &format!(
                "21:1-{MAX_COORD}/A\n21:1-1001/G\n21:1-1002/G\n21:1-5000/DUP\n21:1001-1000/G\n21:1-1000/G\n"
            ),
        );
        let error = parser.next_variant().unwrap().unwrap_err();
        assert!(matches!(error, IoError::SkippedLine(_)), "{error:?}");
        assert_eq!(
            error.to_string(),
            format!(
                "variant size ({}) is bigger than --max_sv_size (1000)",
                MAX_COORD - 1
            )
        );
        assert_eq!(
            parser.next_variant().unwrap().unwrap().allele_string.len(),
            1003
        );
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "variant size (1001) is bigger than --max_sv_size (1000)"
        );
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "variant size (4999) is bigger than --max_sv_size (1000)"
        );
        assert_eq!(parser.next_variant().unwrap().unwrap().allele_string, "-/G");
        assert_eq!(parser.next_variant().unwrap().unwrap().end, 1000);
        assert!(parser.next_variant().is_none());
        let v = one("21:1-2000/G");
        assert_eq!(v.ref_allele.len(), 2000);
        let mut parser = parser_in(&context, "21:1-20000000/FOO\n21:1-1002/N\n21:1-1004/NGN\n");
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "FOO is not a supported structural variant type"
        );
        let trimmed = parser.next_variant().unwrap().unwrap();
        assert_eq!((trimmed.start, trimmed.end), (2, 1002));
        assert_eq!(trimmed.ref_allele.len(), 1001);
        assert_eq!(trimmed.alt_allele(), b"-");
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "variant size (1001) is bigger than --max_sv_size (1000)"
        );
    }

    /// The span `minimise_alleles` leaves, computed from the allele alone: a
    /// shared leading or trailing `N` shortens it, an allele of the run's
    /// length or with a `-` is not minimised, and an insertion has none. Under
    /// `--minimal` every record is minimised (Parser.pm:852, no length or `-`
    /// test), so `21:14221104-14221107/NNAA` leaves `NN/AA` over 14221106-14221107,
    /// the span 1 that `--max_sv_size 2` admits in Ensembl VEP 116.2, and `N-A`
    /// over three bases the span 1 of `NN/-A`.
    #[test]
    fn minimised_span_follows_trim_sequences() {
        assert_eq!(minimised_span(1, 1002, b"N", false), 1000);
        assert_eq!(minimised_span(1, 1002, b"G", false), 1001);
        assert_eq!(minimised_span(1, 1003, b"NGN", false), 1000);
        assert_eq!(minimised_span(1, 3, b"NNN", false), 2);
        assert_eq!(minimised_span(1, 3, b"NN", false), 0);
        assert_eq!(minimised_span(1, 3, b"NNNNN", false), 0);
        assert_eq!(minimised_span(1, 3, b"-", false), 2);
        assert_eq!(minimised_span(5, 4, b"NNN", false), 0);
        assert_eq!(minimised_span(14221104, 14221107, b"NNAA", false), 3);
        assert_eq!(minimised_span(14221104, 14221107, b"NNAA", true), 1);
        assert_eq!(minimised_span(1, 3, b"N-A", true), 1);
        assert_eq!(minimised_span(1, 3, b"NNN", true), 0);
        assert_eq!(minimised_span(1, 1002, b"G", true), 1001);
        assert_eq!(minimised_span(5, 4, b"NNN", true), 0);
    }

    /// Every line counts toward Ensembl VEP's line number (Parser/Region.pm:167
    /// per record, Parser.pm:275 per empty line). A line with no word
    /// character, blank, `-/-` or a bare `#`, is an empty line (Parser.pm:273),
    /// warned about in `skip_empty_lines`'s words when its run ends, at the
    /// next record (`:273-290`, the run's first line named) or at the end of
    /// the input; a `#` line with a word character is a comment, counted and
    /// passed over without a warning, where ensembl-io's list parser reads it
    /// as a record (ListBasedParser.pm:60) the expression refuses
    /// (Parser/Region.pm:171). The `\r` of a `\r\n` line end is removed with
    /// the `\n`, and a last line without a line end is read whole. Ensembl VEP
    /// 116.2 numbers `21:14221107-14221107/A`, `# c`, `#`, `  `,
    /// `21:1-20000000/A` as `Skipped 2 empty lines from line 3` and
    /// `line 5 skipped`.
    #[test]
    fn every_line_is_counted_and_empty_runs_warned() {
        let mut parser =
            parser("\n   \n-/-\n# a comment\n21:100-100/G\r\n\n#\n# c\n\n21:200-200/T");
        let first = parser.next_variant().unwrap().unwrap();
        assert_eq!(first.uploaded_variation(), "21:100-100/G");
        assert_eq!(parser.record_line_number(), Some(5));
        assert_eq!(
            parser.take_warnings(),
            ["Skipped 3 empty lines from line 1"]
        );
        let second = parser.next_variant().unwrap().unwrap();
        assert_eq!(second.uploaded_variation(), "21:200-200/T");
        assert_eq!(parser.record_line_number(), Some(10));
        assert_eq!(
            parser.take_warnings(),
            [
                "Skipped 2 empty lines from line 6",
                "Skipped 1 empty line from line 9"
            ]
        );
        assert!(parser.next_variant().is_none());
        assert!(parser.take_warnings().is_empty());
        assert_eq!(parser.take_empty_lines(), 0);

        // A run the input closes on is warned about at the end of the input,
        // and a `#` line before any record is line 1.
        let mut trailing = parser_in(&ParserContext::default(), "# c\n21:100-100/G\n\n\n");
        let only = trailing.next_variant().unwrap().unwrap();
        assert_eq!(only.uploaded_variation(), "21:100-100/G");
        assert_eq!(trailing.record_line_number(), Some(2));
        assert!(trailing.take_warnings().is_empty());
        assert!(trailing.next_variant().is_none());
        assert_eq!(
            trailing.take_warnings(),
            ["Skipped 2 empty lines from line 3"]
        );
    }

    /// A coordinate is read up to `MAX_COORD`; a start or end past it is a line
    /// the expression matches and the parser still refuses, naming the bound.
    #[test]
    fn coordinates_past_the_ceiling_fail_the_line() {
        let past = MAX_COORD + 1;
        let kept = one(&format!("21:{MAX_COORD}-{MAX_COORD}:1/G"));
        assert_eq!((kept.start, kept.end), (MAX_COORD, MAX_COORD));
        for region in [
            format!("21:{past}-{past}:1/G"),
            format!("21:100-{}:1/G", u64::MAX),
        ] {
            assert_eq!(
                refused(&region),
                format!("Unable to parse region '{region}': a coordinate is over {MAX_COORD}")
            );
        }
    }
}
