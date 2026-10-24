// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! HGVS notation input, one notation per line, read as Ensembl VEP's
//! `Parser/HGVS.pm` `create_VariationFeatures` (release/116.2, lines 141-252)
//! reads it over ensembl-io's list parser and resolves it through
//! `VariationFeatureAdaptor::fetch_by_hgvs_notation` (ensembl-variation
//! release/116, lines 1767-1991):
//!
//! ```text
//! 21:g.25585733A>G
//! NC_000021.9:g.25585656_25585660del
//! MT:m.8993T>G
//! ```
//!
//! Whitespace is removed from the line (`:156`), a notation with a `[` or `]`
//! (a repeat) is skipped with the parser's warning (`:158-161`), and the
//! notation is split into its reference, type and description after a
//! `(gene)` qualifier and a `(p.…)` suffix are removed (`Notation`). The
//! genomic types `g.` and `m.` are read by [`genomic`] against the run's
//! reference FASTA, which stands in for the slice Ensembl VEP fetches from its
//! database (its own parser accepts no HGVS input offline, `:104`); the
//! transcript types `c.` and `n.` by [`transcript`] and the protein type `p.`
//! by [`protein`], each against the run's transcript cache and the FASTA,
//! which stand in for its core database, with the parser's warning when the
//! reference reads as a gene name (`:219-229`). An `r.` notation is one
//! Ensembl VEP's resolver throws on (`HGVS type r not recognised`, `:1989`)
//! and its parser then skips with that message; here it fails the line with
//! [`IoError::Unsupported`], which stops the run. A transcript cache a `c.`,
//! `n.` or `p.` lookup cannot read fails the line with
//! [`IoError::CacheUnreadable`] carrying the cache's refusal, which stops the
//! run under the `cache` component.
//!
//! A notation Ensembl VEP cannot resolve is warned about as `Unable to parse
//! HGVS notation '<notation>'` followed by the exception its resolver threw
//! (`:204-216`) and the line is skipped; here such a line fails with
//! [`IoError::SkippedLine`] carrying that text and the exception's message on one
//! line, for the caller to warn about and skip. The record is named by
//! the notation as written without its whitespace (`:237, 242`), which is what
//! `Uploaded_variation`, the JSON `id` and the JSON `input` print. The lines
//! are numbered as Ensembl VEP numbers them (`LineCount`): a line with no
//! word character is an empty line, warned about in `skip_empty_lines`'s words
//! and counted (Parser.pm:273), and any other `#` line is a comment here,
//! counted toward the line numbers that follow and otherwise passed over
//! without a word, where ensembl-io's list parser has no comment lines
//! (ListBasedParser.pm:60) and Ensembl VEP warns about the notation it cannot
//! parse (`:204-216`).

mod components;
pub mod genomic;
pub mod protein;
pub mod transcript;

use std::io::BufRead;
use std::sync::Arc;

use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::{
    has_word_character, is_perl_space, InputParser, LineCount, ParserContext, ReferenceSequence,
    TranscriptSource,
};
use components::ReadFlags;

/// Parser for HGVS notation input.
pub struct HgvsParser<R: BufRead> {
    reader: R,
    reference: Option<Arc<dyn ReferenceSequence>>,
    /// The run's transcripts, which a `c.`, `n.` or `p.` notation is resolved
    /// against.
    transcripts: Option<Arc<dyn TranscriptSource>>,
    /// `--lookup_ref` (`-replace_ref`, Parser/HGVS.pm:192) and `--minimal`
    /// ([`ParserContext::minimal`]), as the resolver reads them.
    flags: ReadFlags,
    /// `--max_sv_size` as given, the span a genomic notation must stay under
    /// (`max_size`, `:169`): `None` is the flag's `-1`, kept and compared as
    /// written (`$end - $start < -1`), so only a span written end first by two
    /// or more bases passes under it; `0` reads as 5000. A transcript
    /// notation's record is bounded by `validate_vf` instead (Parser.pm:491-499,
    /// the span being end minus start), `None` for no bound.
    max_sv_size: Option<u64>,
    line_buf: Vec<u8>,
    /// Ensembl VEP's line number over the lines read, and the empty-line
    /// warnings it raises.
    lines: LineCount,
    /// The notations' own warnings since the last [`InputParser::take_warnings`].
    warnings: Vec<String>,
}

impl<R: BufRead> HgvsParser<R> {
    /// A parser over `reader` resolving its notations against the context's
    /// reference and transcripts under its `--max_sv_size`, replacing each
    /// written reference allele by the genome's when `lookup_ref` is set and
    /// trimming each record's alleles as written under its `--minimal`.
    pub fn new(reader: R, context: &ParserContext, lookup_ref: bool) -> Self {
        Self {
            reader,
            reference: context.fasta.clone(),
            transcripts: context.transcripts.clone(),
            flags: ReadFlags {
                lookup_ref,
                minimal: context.minimal,
            },
            max_sv_size: context.max_sv_size,
            line_buf: Vec::new(),
            lines: LineCount::default(),
            warnings: Vec::new(),
        }
    }

    /// The variant of one notation, its whitespace already removed.
    fn read_notation(&mut self, hgvs: &str) -> Result<InputVariant, IoError> {
        if hgvs.contains(['[', ']']) {
            return Err(unable(hgvs, None));
        }
        let notation = Notation::split(hgvs).map_err(|message| unable(hgvs, Some(&message)))?;
        let no_reference = || {
            IoError::Unsupported(
                "Cannot use HGVS format in offline mode without a FASTA file".into(),
            )
        };
        match notation.kind.map(|kind| kind.to_ascii_lowercase()) {
            Some('g' | 'm') => {
                let reference = self.reference.as_ref().ok_or_else(no_reference)?;
                genomic::variant(
                    hgvs,
                    &notation,
                    reference.as_ref(),
                    self.flags,
                    self.max_sv_size,
                )
            }
            Some(kind @ ('c' | 'n' | 'p')) => {
                let reference = self.reference.as_ref().ok_or_else(no_reference)?;
                let Some(transcripts) = &self.transcripts else {
                    return Err(IoError::Unsupported(
                        "Cannot use HGVS transcript notation without a transcript cache".into(),
                    ));
                };
                let variant = if kind == 'p' {
                    protein::variant(
                        hgvs,
                        &notation,
                        reference.as_ref(),
                        transcripts.as_ref(),
                        self.flags,
                    )?
                } else {
                    transcript::variant(
                        hgvs,
                        &notation,
                        reference.as_ref(),
                        transcripts.as_ref(),
                        self.flags,
                    )?
                };
                // Parser/HGVS.pm:219 warns before the base parser's `validate_vf`
                // (Parser.pm:189) reads the record's span.
                if let Some(warning) = gene_identifier_warning(hgvs) {
                    self.warnings.push(warning);
                }
                let span = variant.end.saturating_sub(variant.start);
                match self.max_sv_size {
                    Some(max) if span > max => Err(IoError::SkippedLine(format!(
                        "variant size ({span}) is bigger than --max_sv_size ({max})"
                    ))),
                    _ => Ok(variant),
                }
            }
            Some('r') => Err(IoError::Unsupported(
                "HGVS r. notation is not supported (Ensembl VEP refuses it too: HGVS type r not recognised)"
                    .into(),
            )),
            other => Err(unable(
                hgvs,
                Some(&format!(
                    "HGVS type {} not recognised for {}",
                    other.map(String::from).unwrap_or_default(),
                    notation.cleaned
                )),
            )),
        }
    }
}

/// The skip Ensembl VEP's parser writes for a notation its resolver refused
/// (Parser/HGVS.pm:159, 214): the notation as written, then the resolver's
/// message where there is one.
pub(super) fn unable(hgvs: &str, message: Option<&str>) -> IoError {
    IoError::SkippedLine(match message {
        Some(message) => format!("Unable to parse HGVS notation '{hgvs}': {message}"),
        None => format!("Unable to parse HGVS notation '{hgvs}'"),
    })
}

/// The failure of a `c.`, `n.` or `p.` lookup whose transcript cache cannot be
/// read: the cache's refusal `message` after the parser's text, as
/// [`IoError::CacheUnreadable`].
pub(super) fn cache_unreadable(message: String) -> IoError {
    IoError::CacheUnreadable(format!(
        "Cannot read the transcript cache for HGVS input: {message}"
    ))
}

/// The warning Ensembl VEP's parser writes once a transcript or protein
/// notation resolves under a reference that is not an Ensembl, RefSeq or LRG
/// identifier (Parser/HGVS.pm:219-229, `/\:[cnp]\./` and not
/// `/^(ENS|[NX][CGMRP]\_|LRG\_)/`): the text before the first colon is
/// most likely a gene name, and one transcript of the gene was selected.
fn gene_identifier_warning(hgvs: &str) -> Option<String> {
    let is_transcript_notation = hgvs
        .match_indices(':')
        .any(|(at, _)| matches!(&hgvs.as_bytes()[at + 1..], [b'c' | b'n' | b'p', b'.', ..]));
    let is_ensembl_refseq_or_lrg = hgvs.starts_with("ENS")
        || hgvs.starts_with("LRG_")
        || matches!(
            hgvs.as_bytes(),
            [b'N' | b'X', b'C' | b'G' | b'M' | b'R' | b'P', b'_', ..]
        );
    if !is_transcript_notation || is_ensembl_refseq_or_lrg {
        return None;
    }
    let reference = hgvs.split(':').next().unwrap_or_default();
    Some(format!(
        "Possible invalid use of gene or protein identifier '{reference}' as HGVS reference; \
         most likely transcript will be selected"
    ))
}

/// A notation split as `fetch_by_hgvs_notation` splits it
/// (VariationFeatureAdaptor.pm:1784-1820): a `(gene)` qualifier before the
/// colon and a `(p.…)` suffix are removed, the reference is the text before
/// the first colon, the type is the letter that follows the colon when it is
/// one of `c`, `g`, `m`, `n`, `r` or `p`, and the description is the rest
/// after that letter's dot. A notation with no colon, an empty reference or no
/// digit, `*` or `-` after the colon cannot be parsed; a ClinVar range with an
/// unknown bound (`-?_`, `(?_`) and an insertion of unspecified bases
/// (`ins56`, `(100_200)insTC`) are refused with Ensembl VEP's messages; any
/// other parenthesised text in the description cannot be interpreted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Notation {
    /// The sequence the notation names before the colon, with the version
    /// removed from an Ensembl or LRG identifier (`:1816-1818`).
    pub reference: String,
    /// The type letter as written, `None` when the description follows the
    /// colon directly.
    pub kind: Option<char>,
    /// The variant description after the type and its dot.
    pub description: String,
    /// The notation with its qualifier and suffix removed, the form Ensembl
    /// VEP's messages quote.
    pub cleaned: String,
}

impl Notation {
    pub(super) fn split(hgvs: &str) -> Result<Self, String> {
        let cleaned = without_qualifiers(hgvs);
        let could_not_parse = || format!("Could not parse the HGVS notation {cleaned}");
        let (reference, rest) = cleaned
            .split_once(':')
            .filter(|(reference, _)| !reference.is_empty())
            .ok_or_else(could_not_parse)?;
        if !rest.contains(|c: char| c.is_ascii_digit() || c == '*' || c == '-') {
            return Err(could_not_parse());
        }
        let kind = rest
            .chars()
            .next()
            .filter(|c| matches!(c.to_ascii_lowercase(), 'c' | 'g' | 'm' | 'n' | 'r' | 'p'));
        let after_kind = if kind.is_some() { &rest[1..] } else { rest };
        let description = after_kind.strip_prefix('.').unwrap_or(after_kind);
        if has_unknown_location(description) {
            return Err(
                "HGVS notation for variation with unknown location is not supported".to_string(),
            );
        }
        if is_imprecise_insertion(description) {
            return Err(format!(
                "HGVS notation for insertion '{description}' is not supported"
            ));
        }
        let description = match parenthesised(description) {
            Some("(p.=)") => &description[..description.find('(').unwrap_or(0)],
            Some(extra) => {
                return Err(format!(
                    "Could not parse the HGVS notation {cleaned} - can't interpret '{extra}'"
                ))
            }
            None => description,
        };
        let reference = if reference.starts_with("ENS") || is_lrg(reference) {
            without_versions(reference)
        } else {
            reference.to_string()
        };
        Ok(Self {
            reference,
            kind,
            description: description.to_string(),
            cleaned,
        })
    }
}

/// `s/\(.+\):/:/` then `s/\(p\..+\)//` (VariationFeatureAdaptor.pm:1785-1786):
/// the text from the first `(` to the last `):` goes, the colon staying, then
/// the text from the first `(p.` to the last `)`.
fn without_qualifiers(hgvs: &str) -> String {
    let mut text = hgvs.to_string();
    if let (Some(open), Some(close)) = (text.find('('), text.rfind("):")) {
        if close > open + 1 {
            text.replace_range(open..=close, "");
        }
    }
    if let (Some(open), Some(close)) = (text.find("(p."), text.rfind(')')) {
        if close > open + 3 {
            text.replace_range(open..=close, "");
        }
    }
    text
}

/// `(\-\?\_.+\?)|(\(\?\_.+\))` (`:1798`): a ClinVar range whose bound is
/// unknown, `88957372-?_89042944+?dup` or `(?_145413072)(145592772_145596903)del`.
fn has_unknown_location(description: &str) -> bool {
    let follows = |marker: &str, closer: char| {
        description.find(marker).is_some_and(|at| {
            description[at + marker.len()..]
                .rfind(closer)
                .is_some_and(|close| close > 0)
        })
    };
    follows("-?_", '?') || follows("(?_", ')')
}

/// `\(.+\_.+\)ins|ins\(?\d+\)?` (`:1805`, lower case only): an insertion
/// between two positions given as a range, or of a number of unspecified bases.
fn is_imprecise_insertion(description: &str) -> bool {
    let ranged = description.match_indices(")ins").any(|(close, _)| {
        description[..close].find('(').is_some_and(|open| {
            let between = &description[open + 1..close];
            between
                .match_indices('_')
                .any(|(underscore, _)| underscore > 0 && underscore + 1 < between.len())
        })
    });
    let counted = description.match_indices("ins").any(|(at, _)| {
        let after = &description[at + 3..];
        let after = after.strip_prefix('(').unwrap_or(after);
        after.starts_with(|c: char| c.is_ascii_digit())
    });
    ranged || counted
}

/// The parenthesised text `(.+?)(\(.+\))` captures when `\(.+\)` matches the
/// description (`:1810-1813`): from the first `(` that has at least one
/// character before it and a `)` two or more characters after it, to the last
/// `)`; the empty string when the description opens with its only such `(`,
/// which leaves the capture undefined; `None` when no `(` encloses a character.
fn parenthesised(description: &str) -> Option<&str> {
    let close = description.rfind(')')?;
    let encloses = |open: usize| close >= open + 2;
    if !description.find('(').is_some_and(encloses) {
        return None;
    }
    let extra = description
        .match_indices('(')
        .map(|(open, _)| open)
        .find(|open| *open > 0)
        .filter(|open| encloses(*open))
        .map_or("", |open| &description[open..=close]);
    Some(extra)
}

/// `^LRG_\d+`: an LRG identifier, whose version is removed as an Ensembl one's is.
fn is_lrg(reference: &str) -> bool {
    reference
        .strip_prefix("LRG_")
        .is_some_and(|rest| rest.starts_with(|c: char| c.is_ascii_digit()))
}

/// `s/\.\d+//g`: every dot followed by digits removed.
fn without_versions(reference: &str) -> String {
    let mut out = String::with_capacity(reference.len());
    let mut chars = reference.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '.' && chars.peek().is_some_and(char::is_ascii_digit) {
            while chars.peek().is_some_and(char::is_ascii_digit) {
                chars.next();
            }
        } else {
            out.push(c);
        }
    }
    out
}

impl<R: BufRead + Send> InputParser for HgvsParser<R> {
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
            // `skip_empty_lines` reads a line with no word character past, a
            // bare `#` included (Parser.pm:273), before the parser sees it.
            if !has_word_character(&text) {
                self.lines.empty_line();
                continue;
            }
            self.lines.record_line();
            if text.starts_with('#') {
                continue;
            }
            // `s/\s+//g` (Parser/HGVS.pm:156): the line end goes with every
            // other blank.
            let hgvs: String = text.chars().filter(|c| !is_perl_space(*c)).collect();
            return Some(self.read_notation(&hgvs));
        }
    }

    fn take_warnings(&mut self) -> Vec<String> {
        let mut warnings = self.lines.take_warnings();
        warnings.append(&mut self.warnings);
        warnings
    }

    fn record_line_number(&self) -> Option<u64> {
        Some(self.lines.line())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::io::Cursor;

    /// A reference of named sequences, each known from an offset, reached by
    /// its own name or an alias.
    struct Genome {
        sequences: HashMap<&'static str, (u64, &'static str)>,
        aliases: HashMap<&'static str, &'static str>,
    }

    impl Genome {
        /// GRCh38 chromosome 21 from 30,282,140 and MT from 8,990, with the
        /// accession and `chr` spellings of chromosome 21 as aliases.
        fn grch38() -> Arc<dyn ReferenceSequence> {
            Arc::new(Self {
                sequences: HashMap::from([
                    ("21", (30_282_140, "TCCAGCTGCTATAAATTCTGGAAGTTGGTCAATAGC")),
                    ("MT", (8_990, "CCCTGGC")),
                ]),
                aliases: HashMap::from([("NC_000021.9", "21"), ("chr21", "21")]),
            })
        }
    }

    impl ReferenceSequence for Genome {
        fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
            let (offset, bases) = self.sequences.get(self.record_name(chr)?.as_str())?;
            let first = usize::try_from(start.checked_sub(*offset)?).ok()?;
            let last = usize::try_from(end.checked_sub(*offset)?).ok()?;
            bases.get(first..=last).map(|s| s.as_bytes().to_vec())
        }

        fn record_name(&self, chr: &str) -> Option<String> {
            if self.sequences.contains_key(chr) {
                return Some(chr.to_string());
            }
            self.aliases.get(chr).map(|name| name.to_string())
        }

        fn sequence_length(&self, chr: &str) -> Option<u64> {
            let (offset, bases) = self.sequences.get(self.record_name(chr)?.as_str())?;
            Some(offset + bases.len() as u64 - 1)
        }
    }

    /// The fake reference under Ensembl VEP's default `--max_sv_size`
    /// (Config.pm:315).
    pub(super) fn context() -> ParserContext {
        context_with(Some(10_000_000))
    }

    pub(super) fn context_with(max_sv_size: Option<u64>) -> ParserContext {
        ParserContext {
            fasta: Some(Genome::grch38()),
            max_sv_size,
            ..ParserContext::default()
        }
    }

    pub(super) fn parser(data: &str) -> HgvsParser<Cursor<Vec<u8>>> {
        HgvsParser::new(Cursor::new(data.as_bytes().to_vec()), &context(), false)
    }

    /// The reference, type letter and description of each notation shape,
    /// split as `fetch_by_hgvs_notation` splits them: a `(gene)` qualifier and
    /// a `(p.=)` suffix go, the version comes off an Ensembl or LRG identifier
    /// and stays on an accession, the type letter keeps its case, and a
    /// notation whose description follows the colon directly has no type.
    #[test]
    fn notation_splits_as_fetch_by_hgvs_notation_does() {
        let split = |hgvs: &str| Notation::split(hgvs).unwrap();
        let n = split("21:g.25585733A>G");
        assert_eq!(
            (n.reference.as_str(), n.kind, n.description.as_str()),
            ("21", Some('g'), "25585733A>G")
        );
        let n = split("NM_000088.3(COL1A1):c.589G>T(p.=)");
        assert_eq!(
            (
                n.reference.as_str(),
                n.kind,
                n.description.as_str(),
                n.cleaned.as_str()
            ),
            ("NM_000088.3", Some('c'), "589G>T", "NM_000088.3:c.589G>T")
        );
        let n = split("ENST00000366667.4:c.803C>T");
        assert_eq!(n.reference, "ENST00000366667");
        let n = split("LRG_8t1.2:c.3891A>T");
        assert_eq!(n.reference, "LRG_8t1");
        let n = split("NC_000021.9:G.25585733A>G");
        assert_eq!((n.reference.as_str(), n.kind), ("NC_000021.9", Some('G')));
        let n = split("21:25585733A>G");
        assert_eq!((n.kind, n.description.as_str()), (None, "25585733A>G"));
        let n = split("MT:m.8993T>G");
        assert_eq!((n.reference.as_str(), n.kind), ("MT", Some('m')));
    }

    /// The notations `fetch_by_hgvs_notation` refuses before resolving them:
    /// no colon or no digit, a ClinVar range with an unknown bound, an
    /// insertion of unspecified bases (lower-case `ins` only, as the regular
    /// expression at `:1805` is), and a parenthesised text other than the
    /// protein suffix: the capture runs from the first `(` with text before it
    /// to the last `)`, is empty for a description opening with its only `(`,
    /// and an empty pair of parentheses is no text to interpret. The capture is
    /// found by character, so a description opening with a character of more
    /// than one byte is refused like any other.
    #[test]
    fn notation_refusals_carry_ensembl_vep_messages() {
        let error = |hgvs: &str| Notation::split(hgvs).unwrap_err();
        assert_eq!(
            error("21.g.123A>G"),
            "Could not parse the HGVS notation 21.g.123A>G"
        );
        assert_eq!(
            error(":g.123A>G"),
            "Could not parse the HGVS notation :g.123A>G"
        );
        assert_eq!(
            error("21:g.A>G"),
            "Could not parse the HGVS notation 21:g.A>G"
        );
        assert_eq!(
            error("13:g.88957372-?_89042944+?dup"),
            "HGVS notation for variation with unknown location is not supported"
        );
        assert_eq!(
            error("13:g.(?_145413072)(145592772_145596903)del"),
            "HGVS notation for variation with unknown location is not supported"
        );
        for description in [
            "(99982241_99982249)insTC",
            "99982241_99982242ins56",
            "99982241_99982242ins(56)",
        ] {
            assert_eq!(
                error(&format!("NC_000013.11:g.{description}")),
                format!("HGVS notation for insertion '{description}' is not supported")
            );
        }
        assert_eq!(
            error("21:g.123A>G(foo)"),
            "Could not parse the HGVS notation 21:g.123A>G(foo) - can't interpret '(foo)'"
        );
        assert_eq!(
            error("21:g.(123)A>G"),
            "Could not parse the HGVS notation 21:g.(123)A>G - can't interpret ''"
        );
        assert_eq!(
            error("21:g.(1)_(2)del"),
            "Could not parse the HGVS notation 21:g.(1)_(2)del - can't interpret '(2)'"
        );
        assert_eq!(
            error("21:g.\u{2013}5(x)A>G"),
            "Could not parse the HGVS notation 21:g.\u{2013}5(x)A>G - can't interpret '(x)'"
        );
        let n = Notation::split("21:g.123()A>G").unwrap();
        assert_eq!(n.description, "123()A>G");
        let n = Notation::split("21:g.123_124INS56").unwrap();
        assert_eq!(n.description, "123_124INS56");
    }

    /// Whitespace inside a notation is removed before anything else and the
    /// record is named by the notation without it, so `21:g.30282169 C>A`
    /// annotates as `21:g.30282169C>A` does and prints that name (Ensembl VEP
    /// 116.2 prints `21:g.30282169C>A` in `Uploaded_variation` for the spaced
    /// line). Every line counts toward Ensembl VEP's line number
    /// (Parser/HGVS.pm:151 per record, Parser.pm:275 per empty line): a `#` line
    /// with a word character is a comment here, counted and passed over without
    /// a warning, where ensembl-io's list parser reads it as a record
    /// (ListBasedParser.pm:60) that `fetch_by_hgvs_notation` refuses
    /// (`:204-216`); a line with no word character, blank or a bare `#`, is an
    /// empty line (Parser.pm:273), warned about in `skip_empty_lines`'s words
    /// when its run ends, at the next record or the end of the input (Ensembl
    /// VEP refuses the format offline, Parser/HGVS.pm:104, so these are its
    /// source's numbers, not a run's).
    #[test]
    fn whitespace_is_removed_and_the_record_named_by_the_notation() {
        let mut parser = parser("# a comment\n\n  \n21:g.30282169 C>A\r\n#\n\n");
        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(parser.record_line_number(), Some(4));
        assert_eq!(
            parser.take_warnings(),
            ["Skipped 2 empty lines from line 2"]
        );
        assert_eq!(v.id.as_deref(), Some("21:g.30282169C>A"));
        assert_eq!(v.raw_input.as_deref(), Some("21:g.30282169C>A"));
        assert_eq!(v.uploaded_variation(), "21:g.30282169C>A");
        assert_eq!(v.variation_name(), "21:g.30282169C>A");
        assert!(v.input_line_name.is_none());
        assert!(parser.next_variant().is_none());
        assert_eq!(
            parser.take_warnings(),
            ["Skipped 2 empty lines from line 5"]
        );
        assert_eq!(parser.take_empty_lines(), 0);
    }

    /// A notation with a `[` or `]` is skipped with the parser's own warning,
    /// which carries no resolver message (Parser/HGVS.pm:158-161); Ensembl VEP
    /// 116.2 warns `Unable to parse HGVS notation '21:g.30282169C[3]'` and
    /// writes no row for it.
    #[test]
    fn repeats_are_skipped_with_the_parser_warning() {
        let mut parser = parser("21:g.30282169C[3]\n21:g.30282169C>A\n");
        let err = parser.next_variant().unwrap().unwrap_err();
        assert!(matches!(err, IoError::SkippedLine(_)), "{err:?}");
        assert_eq!(
            err.to_string(),
            "Unable to parse HGVS notation '21:g.30282169C[3]'"
        );
        assert!(parser.next_variant().unwrap().is_ok());
    }

    /// The types Ensembl VEP's resolver knows (VariationFeatureAdaptor.pm:1829,
    /// 1882, 1900): a `c.`, `n.` or `p.` notation is resolved through the
    /// context's transcripts, and stops the run when the context has none; a
    /// notation with no type letter is skipped with Ensembl VEP's `HGVS type
    /// not recognised` (two spaces, the empty type).
    #[test]
    fn type_letters_dispatch_as_ensembl_vep_dispatches_them() {
        let mut parser = parser(
            "ENST00000366667:c.803C>T\nENST00000366667:n.803C>T\nENSP00000366667:p.Arg97Ter\n\
             21:30282169C>A\n21:g.30282169C>A\n",
        );
        for _ in 0..3 {
            let err = parser.next_variant().unwrap().unwrap_err();
            assert!(matches!(err, IoError::Unsupported(_)), "{err:?}");
            assert_eq!(
                err.to_string(),
                "Cannot use HGVS transcript notation without a transcript cache"
            );
        }
        let err = parser.next_variant().unwrap().unwrap_err();
        assert!(matches!(err, IoError::SkippedLine(_)), "{err:?}");
        assert_eq!(
            err.to_string(),
            "Unable to parse HGVS notation '21:30282169C>A': HGVS type  not recognised for 21:30282169C>A"
        );
        assert!(parser.next_variant().unwrap().is_ok());
    }

    /// An r. notation fails its line with the refusal variant and a message
    /// naming the type, as Ensembl VEP's resolver throws `HGVS type r not
    /// recognised` (DBSQL/VariationFeatureAdaptor.pm:1989), where its parser
    /// warns and skips the line; the parser here reads on, and the caller stops
    /// the run on the variant.
    #[test]
    fn rna_hgvs_is_refused_naming_r() {
        let mut parser = parser("ENST00000366667:r.803C>T\n21:g.30282169C>A\n");
        let err = parser.next_variant().unwrap().unwrap_err();
        assert!(matches!(err, IoError::Unsupported(_)), "{err:?}");
        assert_eq!(
            err.to_string(),
            "HGVS r. notation is not supported (Ensembl VEP refuses it too: HGVS type r not recognised)"
        );
        assert!(parser.next_variant().unwrap().is_ok());
    }

    /// Without a reference the genomic types cannot be resolved: Ensembl VEP
    /// refuses the format offline (Parser/HGVS.pm:104), and the run stops.
    #[test]
    fn genomic_notation_without_a_reference_is_refused() {
        let mut parser = HgvsParser::new(
            Cursor::new(b"21:g.30282169C>A\n".to_vec()),
            &ParserContext::default(),
            false,
        );
        let err = parser.next_variant().unwrap().unwrap_err();
        assert!(matches!(err, IoError::Unsupported(_)), "{err:?}");
        assert_eq!(
            err.to_string(),
            "Cannot use HGVS format in offline mode without a FASTA file"
        );
    }
}
