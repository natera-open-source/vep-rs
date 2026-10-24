// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Identifier input (`--format id`), read as `Parser/ID.pm`
//! `create_VariationFeatures` (release/116.2, lines 138-197) reads it over
//! ensembl-io's `ListBasedParser`: one variant identifier per line,
//!
//! ```text
//! rs699
//! rs2838519
//! ```
//!
//! Ensembl VEP asks its variation database for each name (`:162`,
//! `fetch_by_name` with failed variants included) and makes one variation
//! feature per genomic mapping of the variant, renamed to the identifier as
//! written (`:185-191`), so a name that maps twice yields two records. A name
//! the database lacks is warned about and skipped (`:175-176`). Here the list is
//! read whole, every name is resolved in one pass over the cache's variation
//! files by the caller, and [`IdParser::next_variant`] yields the mappings in
//! the input's order, each mapping a record of its own with the identifier as
//! its name and its input line. The parser warns about the empty lines it read
//! past itself, through [`InputParser::take_warnings`], since Ensembl VEP numbers
//! them by a line count that takes in every line of the list (Parser.pm:275,
//! Parser/ID.pm:148), a name the cache lacks included. The checks `validate_vf`
//! makes on each feature (Parser.pm:189) are the reader's, as for every format;
//! the parser gives the reader the line a skip warning quotes and the record's
//! number in that same count ([`IdParser::attempted_line`],
//! [`InputParser::record_line_number`]).
//!
//! The variation files hold the variants of the cache's chromosomes by the
//! name of their variation feature, so a structural variant identifier
//! (`nsv`, `esv`; Ensembl VEP's second lookup at `:168-173`), a synonym or a
//! subsnp identifier (`fetch_by_name`'s fallbacks), and a variant with no
//! genomic mapping (`No mappings found`, `:181-182`) are all `No variant found`
//! here.

use std::io::BufRead;
use std::sync::Arc;

use vep_core::coordinate::Strand;
use vep_core::variant::{so_variation_class, InputVariant};

use super::ensembl::sequence_variants;
use super::{
    empty_lines_warning, has_word_character, is_perl_space, read_until_bounded, InputParser,
    LineTooLong,
};
use crate::error::IoError;

/// One genomic mapping of a known variant, as a variation file records it:
/// the chromosome the file belongs to, the span, the strand and the allele
/// string (`BaseCacheVariation.pm::parse_variation`, `:115-145`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameMapping {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub allele_string: String,
}

/// One name of the list, the line it was read from, the empty lines that
/// preceded it and, once supplied, its mappings.
struct Entry {
    name: String,
    /// The line as read without its line end, when whitespace makes it differ
    /// from the name: what a skip warning quotes, as Ensembl VEP quotes the
    /// list parser's `current_block` (BaseVEP.pm:846).
    text: Option<String>,
    /// The 1-based line number, Ensembl VEP's `line_number` for this format
    /// (Parser.pm:275 and Parser/ID.pm:148 count empty lines and names alike).
    line: u64,
    empty_lines_before: u64,
    mappings: Vec<NameMapping>,
}

/// Parser for the identifier list format.
pub struct IdParser {
    entries: Vec<Entry>,
    /// The next entry to yield.
    entry_index: usize,
    /// The next mapping of the current entry to build.
    mapping_index: usize,
    /// The entry whose mapping was last built, set before the build so that a
    /// build that failed names it too.
    attempted: Option<usize>,
    /// The alleles of the mapping being yielded that are still to yield.
    remaining_alleles: std::vec::IntoIter<InputVariant>,
    /// Empty lines after the last identifier, warned about once the names are
    /// exhausted and numbered from `line_count`, the lines of the list.
    trailing_empty_lines: u64,
    line_count: u64,
    warnings: Vec<String>,
}

impl IdParser {
    /// Reads every line of the list: a line without a word character is an
    /// empty line (Parser.pm:273, `skip_empty_lines`), any other line is a name
    /// with its whitespace removed (Parser/ID.pm:153, `s/\s+//g`). The names
    /// wait for [`IdParser::supply_mappings`].
    pub fn read<R: BufRead>(mut reader: R) -> Result<Self, IoError> {
        let mut entries = Vec::new();
        let mut empty_lines_before = 0;
        let mut line_number = 0;
        let mut line = Vec::new();
        loop {
            line.clear();
            // A line past the cap is the list's refusal, with the cap's own
            // sentence as the whole message; any other failure stays an I/O error.
            let read = read_until_bounded(&mut reader, b'\n', &mut line, line_number + 1).map_err(
                |error| match LineTooLong::from_io_error(&error) {
                    Some(too_long) => IoError::Unsupported(too_long.to_string()),
                    None => IoError::Io(error),
                },
            )?;
            if read == 0 {
                break;
            }
            line_number += 1;
            let text = String::from_utf8_lossy(&line);
            if !has_word_character(&text) {
                empty_lines_before += 1;
                continue;
            }
            let name: String = text.chars().filter(|c| !is_perl_space(*c)).collect();
            let text = text.trim_end_matches(['\n', '\r']);
            entries.push(Entry {
                text: (text != name).then(|| text.to_string()),
                name,
                line: line_number,
                empty_lines_before,
                mappings: Vec::new(),
            });
            empty_lines_before = 0;
        }
        Ok(Self {
            entries,
            entry_index: 0,
            mapping_index: 0,
            attempted: None,
            remaining_alleles: Vec::new().into_iter(),
            trailing_empty_lines: empty_lines_before,
            line_count: line_number,
            warnings: Vec::new(),
        })
    }

    /// The names of the list, in order, with their whitespace removed.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|entry| entry.name.as_str())
    }

    /// Gives every name its mappings, one list per name in the order of
    /// [`IdParser::names`]; an empty list is a name the cache lacks.
    ///
    /// # Panics
    ///
    /// When `mappings` has not one entry per name.
    pub fn supply_mappings(&mut self, mappings: Vec<Vec<NameMapping>>) {
        assert_eq!(
            mappings.len(),
            self.entries.len(),
            "one mapping list per name"
        );
        for (entry, mappings) in self.entries.iter_mut().zip(mappings) {
            entry.mappings = mappings;
        }
    }

    /// The line, as read without its line end, of the name whose mapping was
    /// last built, whether the build yielded records or failed: the text a
    /// skip warning on the record quotes, as Ensembl VEP quotes the list
    /// parser's `current_block` (BaseVEP.pm:846). `None` before the first
    /// build.
    pub fn attempted_line(&self) -> Option<&str> {
        let entry = &self.entries[self.attempted?];
        Some(entry.text.as_deref().unwrap_or(&entry.name))
    }
}

/// The alleles of one mapping, as Parser/ID.pm:185-191 shapes the
/// feature: the chromosome of the file, the identifier as the name and as
/// the input line (`_line`), the strand of the mapping. An allele string with
/// a `/` is read as the Ensembl-format parser reads one (a bi-allelic indel
/// minimised by `post_process_vfs`, Parser.pm:862-881, one allele per
/// alternate); a string without one (`HGMD_MUTATION`, `dbSNP_novariation`)
/// names a variant with no alternate allele, which Ensembl VEP annotates
/// against nothing and prints upper-cased, as `validate_vf` upper-cases
/// every allele string (Parser.pm:554). The record's class is the one the
/// database stores for the feature, which Ensembl VEP prints as it is
/// (`VariationFeature.pm:1117`) and which the import derived from the allele
/// string as stored, so it is read from that string before any minimisation:
/// `TTTTT/TTTTTT` stays an `indel` where its minimised `-/T` would read as an
/// insertion.
fn mapping_variants(name: &str, mapping: &NameMapping) -> Result<Vec<InputVariant>, IoError> {
    let id = Some(name.to_string());
    let raw_input: Arc<str> = name.into();
    let upper_allele_string = mapping.allele_string.to_ascii_uppercase();
    let stored_class = so_variation_class(&upper_allele_string);
    let mut variants = if mapping.allele_string.contains('/') {
        // A mapping's alleles are stored upper-case and the record is named by
        // the identifier, so the `--minimal` trim after the read leaves the
        // record the written-case trim would.
        sequence_variants(
            mapping.chr.clone(),
            mapping.start,
            mapping.end,
            &mapping.allele_string,
            id,
            false,
        )?
    } else {
        let mut variant = InputVariant::without_allele(
            mapping.chr.clone(),
            mapping.start,
            &upper_allele_string,
            [],
        );
        variant.end = mapping.end;
        variant.id = id;
        vec![variant]
    };
    let strand = if mapping.strand < 0 {
        Strand::Reverse
    } else {
        Strand::Forward
    };
    for variant in &mut variants {
        variant.strand = strand;
        variant.raw_input = Some(raw_input.clone());
        variant.record_so_term = Some(stored_class.to_string());
    }
    Ok(variants)
}

impl InputParser for IdParser {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        if let Some(variant) = self.remaining_alleles.next() {
            return Some(Ok(variant));
        }
        loop {
            let Some(entry) = self.entries.get(self.entry_index) else {
                let trailing = std::mem::take(&mut self.trailing_empty_lines);
                if trailing > 0 {
                    self.warnings.push(empty_lines_warning(
                        trailing,
                        self.line_count - trailing + 1,
                    ));
                }
                return None;
            };
            if self.mapping_index == 0 {
                // `skip_empty_lines` warns before the name's lookup, numbering the
                // run from its first line (Parser.pm:280-288).
                if entry.empty_lines_before > 0 {
                    self.warnings.push(empty_lines_warning(
                        entry.empty_lines_before,
                        entry.line - entry.empty_lines_before,
                    ));
                }
                if entry.mappings.is_empty() {
                    self.warnings
                        .push(format!("No variant found with ID '{}'", entry.name));
                }
            }
            let Some(mapping) = entry.mappings.get(self.mapping_index) else {
                self.entry_index += 1;
                self.mapping_index = 0;
                continue;
            };
            self.mapping_index += 1;
            self.attempted = Some(self.entry_index);
            let mut alleles = match mapping_variants(&entry.name, mapping) {
                Ok(alleles) => alleles.into_iter(),
                Err(e) => return Some(Err(e)),
            };
            let first = alleles.next().expect("a mapping has an allele");
            self.remaining_alleles = alleles;
            return Some(Ok(first));
        }
    }

    fn take_warnings(&mut self) -> Vec<String> {
        std::mem::take(&mut self.warnings)
    }

    /// The list's line number of the name whose mapping was last built: every
    /// mapping of a name takes the name's line, as Ensembl VEP validates the
    /// features of one name under the one `line_number` (Parser.pm:180-189).
    fn record_line_number(&self) -> Option<u64> {
        Some(self.entries[self.attempted?].line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::MAX_LINE_BYTES;
    use std::io::Cursor;

    fn mapping(chr: &str, start: u64, end: u64, alleles: &str) -> NameMapping {
        NameMapping {
            chr: chr.to_string(),
            start,
            end,
            strand: 1,
            allele_string: alleles.to_string(),
        }
    }

    fn parser(list: &str, mappings: Vec<Vec<NameMapping>>) -> IdParser {
        let mut parser = IdParser::read(Cursor::new(list.as_bytes().to_vec())).unwrap();
        parser.supply_mappings(mappings);
        parser
    }

    fn all(parser: &mut IdParser) -> Vec<InputVariant> {
        let mut out = Vec::new();
        while let Some(v) = parser.next_variant() {
            out.push(v.unwrap());
        }
        out
    }

    /// rs2838519 of the release 116 GRCh38 cache (`21 rs2838519 ... 44195140 .
    /// G/A/C`): one feature at its locus with the identifier as its name and
    /// its line, one allele per alternate (Parser/ID.pm:185-191), so
    /// `Uploaded_variation` is `rs2838519` and the Location `21:44195140`.
    #[test]
    fn rsid_resolves_to_locus() {
        let mut parser = parser(
            "rs2838519\n",
            vec![vec![mapping("21", 44195140, 44195140, "G/A/C")]],
        );
        let variants = all(&mut parser);
        assert_eq!(variants.len(), 2);
        let first = &variants[0];
        assert_eq!(first.chr, "21");
        assert_eq!((first.start, first.end), (44195140, 44195140));
        assert_eq!(first.ref_allele, b"G");
        assert_eq!(first.alt_allele(), b"A");
        assert_eq!(variants[1].alt_allele(), b"C");
        assert_eq!(first.id.as_deref(), Some("rs2838519"));
        assert_eq!(first.uploaded_variation(), "rs2838519");
        assert_eq!(first.location(), "21:44195140");
        assert_eq!(first.raw_input.as_deref(), Some("rs2838519"));
        assert_eq!(first.uploaded_allele_string.as_deref(), Some("G/A/C"));
        assert_eq!(first.allele_index, 0);
        assert_eq!(variants[1].allele_index, 1);
        assert_eq!(first.strand, Strand::Forward);
        assert!(parser.take_warnings().is_empty());
    }

    /// A list whose second line runs past the cap is refused as a whole, with
    /// the cap's sentence naming line 2 as the message, once.
    #[test]
    fn list_line_past_the_cap_is_refused() {
        let mut list = b"rs1\n".to_vec();
        list.extend(std::iter::repeat_n(b'x', MAX_LINE_BYTES + 1));
        list.push(b'\n');
        let error = IdParser::read(Cursor::new(list)).err().expect("refused");
        assert_eq!(
            error.to_string(),
            format!("line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line of input")
        );
        assert!(matches!(error, IoError::Unsupported(_)));
    }

    /// A name with two mappings yields two records (Parser/ID.pm:185: one
    /// feature per mapping), each named by the identifier, each opening at
    /// allele index 0, in the order the mappings were given; a reverse-strand
    /// mapping keeps its strand.
    #[test]
    fn multi_mapping_rsid_yields_several() {
        let mut reverse = mapping("22", 500, 500, "T/C");
        reverse.strand = -1;
        let mut parser = parser("rs1\n", vec![vec![mapping("21", 100, 100, "A/G"), reverse]]);
        let variants = all(&mut parser);
        assert_eq!(variants.len(), 2);
        assert_eq!(
            (
                variants[0].chr.as_str(),
                variants[0].start,
                variants[0].allele_index
            ),
            ("21", 100, 0)
        );
        assert_eq!(
            (
                variants[1].chr.as_str(),
                variants[1].start,
                variants[1].allele_index
            ),
            ("22", 500, 0)
        );
        assert_eq!(variants[1].strand, Strand::Reverse);
        assert!(variants.iter().all(|v| v.id.as_deref() == Some("rs1")));
    }

    /// A name the cache lacks yields nothing and the warning of
    /// Parser/ID.pm:175 (`No variant found with ID '<id>'`), in the input's
    /// order with the names around it; the names after it are still read. An
    /// empty line after it is warned about after it, as `skip_empty_lines` runs
    /// when the next name is read (Parser/ID.pm:144), and takes its line number
    /// from a count the unmapped name advanced too (`:148`).
    #[test]
    fn unknown_id_skipped_with_warning() {
        let mut parser = parser(
            "rs1\nrs404\n\nrs2\n",
            vec![
                vec![mapping("21", 100, 100, "A/G")],
                vec![],
                vec![mapping("21", 200, 200, "C/T")],
            ],
        );
        let first = parser.next_variant().unwrap().unwrap();
        assert_eq!(first.start, 100);
        assert!(parser.take_warnings().is_empty());
        let second = parser.next_variant().unwrap().unwrap();
        assert_eq!(second.start, 200);
        assert_eq!(
            parser.take_warnings(),
            [
                "No variant found with ID 'rs404'".to_string(),
                "Skipped 1 empty line from line 3".to_string(),
            ]
        );
        assert!(parser.next_variant().is_none());
        assert!(parser.take_warnings().is_empty());
    }

    /// A structural variant identifier is not in the variation files, so it
    /// resolves to nothing and takes the `No variant found` warning where
    /// Ensembl VEP's second lookup (Parser/ID.pm:168-173, the structural
    /// variation adaptor) would find it: the documented exclusion.
    #[test]
    fn nsv_esv_documented_unsupported() {
        let mut parser = parser("nsv1\nesv2\n", vec![vec![], vec![]]);
        assert!(parser.next_variant().is_none());
        assert_eq!(
            parser.take_warnings(),
            [
                "No variant found with ID 'nsv1'".to_string(),
                "No variant found with ID 'esv2'".to_string()
            ]
        );
    }

    /// Whitespace inside a line is removed before the lookup (Parser/ID.pm:153)
    /// and a `\r\n` line end is dropped; a line with no word character is an
    /// empty line (Parser.pm:273), warned about in `skip_empty_lines`'s words
    /// before the name that follows it or, after the last name, at the end of
    /// the input, each run numbered from its first line as Ensembl VEP counts
    /// the list's lines (Parser.pm:284-288): lines 1-2, line 4, lines 6-7 here.
    #[test]
    fn whitespace_is_removed_and_empty_lines_warned() {
        let mut parser =
            IdParser::read(Cursor::new(b"\n  \nrs 699\r\n---\nrs2\n\n\n".to_vec())).unwrap();
        assert_eq!(parser.names().collect::<Vec<_>>(), ["rs699", "rs2"]);
        parser.supply_mappings(vec![
            vec![mapping("1", 230710048, 230710048, "A/G")],
            vec![mapping("21", 200, 200, "C/T")],
        ]);
        let warning = |text: &str| [text.to_string()];
        let first = parser.next_variant().unwrap().unwrap();
        assert_eq!(first.id.as_deref(), Some("rs699"));
        assert_eq!(
            parser.take_warnings(),
            warning("Skipped 2 empty lines from line 1")
        );
        let second = parser.next_variant().unwrap().unwrap();
        assert_eq!(second.id.as_deref(), Some("rs2"));
        assert_eq!(
            parser.take_warnings(),
            warning("Skipped 1 empty line from line 4")
        );
        assert!(parser.next_variant().is_none());
        assert_eq!(
            parser.take_warnings(),
            warning("Skipped 2 empty lines from line 6")
        );
        assert!(parser.next_variant().is_none());
        assert!(parser.take_warnings().is_empty());
        assert_eq!(parser.take_empty_lines(), 0);
    }

    /// An un-minimised bi-allelic indel of the files (`TCT/TCTC` over three
    /// bases) is trimmed as `post_process_vfs` trims every parsed feature
    /// (Parser.pm:862-881), to `-/C` after the third base, while its class
    /// stays the stored string's `indel` (`VariationFeature.pm:1117`); an
    /// allele string without a `/` (rs782416306 of the release 116 GRCh38
    /// cache) is a variant with no alternate allele, carried without a row to
    /// annotate, upper-cased as `validate_vf` upper-cases it (Parser.pm:554)
    /// and classed `sequence_alteration`.
    #[test]
    fn mapping_alleles_are_post_processed_as_ensembl_vep_does() {
        let mut parser = parser(
            "rs1\nrs782416306\n",
            vec![
                vec![mapping("21", 100, 102, "TCT/TCTC")],
                vec![mapping("21", 44597428, 44597437, "dbSNP_novariation")],
            ],
        );
        let variants = all(&mut parser);
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].allele_string, "-/C");
        assert_eq!((variants[0].start, variants[0].end), (103, 102));
        assert_eq!(
            variants[0].original_allele_string.as_deref(),
            Some("TCT/TCTC")
        );
        // The stored class of `TCT/TCTC`, not the minimised `-/C`'s insertion.
        assert_eq!(variants[0].class_so_term(), "indel");
        assert!(variants[1].no_allele);
        assert_eq!(variants[1].allele_string, "DBSNP_NOVARIATION");
        assert_eq!(variants[1].class_so_term(), "sequence_alteration");
        assert_eq!((variants[1].start, variants[1].end), (44597428, 44597437));
        assert_eq!(variants[1].id.as_deref(), Some("rs782416306"));
    }

    /// The line and the list line number the reader quotes and numbers a skip
    /// with are the name's whose mapping was built last, whether the build
    /// yielded records or failed on its allele string (`A/`): the first name
    /// here, written `rs1 ` and read at line 2 after an empty line, for each of
    /// its two mappings and the alleles of the second; the second name, at line
    /// 4 past a name the cache lacks, for its failed build. Nothing is quoted
    /// before the first build.
    #[test]
    fn attempted_line_names_the_entry_last_built_on_ok_and_err_alike() {
        fn attempted(parser: &IdParser) -> (Option<&str>, Option<u64>) {
            (parser.attempted_line(), parser.record_line_number())
        }
        let mut parser = IdParser::read(Cursor::new(b"\nrs1 \nrs404\nrs2\n".to_vec())).unwrap();
        parser.supply_mappings(vec![
            vec![
                mapping("21", 100, 100, "A/G"),
                mapping("22", 500, 500, "C/T/G"),
            ],
            vec![],
            vec![mapping("21", 200, 200, "A/")],
        ]);
        assert_eq!(attempted(&parser), (None, None));
        for _ in 0..3 {
            assert!(parser.next_variant().unwrap().is_ok());
            assert_eq!(attempted(&parser), (Some("rs1 "), Some(2)));
        }
        let failed = parser.next_variant().unwrap().unwrap_err();
        assert!(failed.to_string().contains("'A/'"), "{failed}");
        assert_eq!(attempted(&parser), (Some("rs2"), Some(4)));
        assert_eq!(
            parser.take_warnings(),
            [
                "Skipped 1 empty line from line 1".to_string(),
                "No variant found with ID 'rs404'".to_string(),
            ]
        );
        assert!(parser.next_variant().is_none());
        assert_eq!(attempted(&parser), (Some("rs2"), Some(4)));
    }

    /// A name repeated in the list yields its mappings each time.
    #[test]
    fn repeated_name_yields_its_mappings_each_time() {
        let site = mapping("21", 100, 100, "A/G");
        let mut parser = parser("rs1\nrs1\n", vec![vec![site.clone()], vec![site]]);
        let variants = all(&mut parser);
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].start, variants[1].start);
    }

    #[test]
    #[should_panic(expected = "one mapping list per name")]
    fn mappings_must_match_the_names() {
        let mut parser = IdParser::read(Cursor::new(b"rs1\nrs2\n".to_vec())).unwrap();
        parser.supply_mappings(vec![vec![]]);
    }
}
