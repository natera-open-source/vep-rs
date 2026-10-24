// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Input parsers for variant formats, the detection of a format from the
//! first data line ([`detect`]), and what the parsers share: the field
//! delimiter of a column-based input ([`Delimiter`]), the run's reference
//! sequence, chromosome synonyms and chromosome names ([`ParserContext`]), the
//! minimisation Ensembl VEP applies to every parsed record
//! ([`minimise_alleles`]), and the bound on one line's length every reader of
//! an input file keeps ([`read_until_bounded`]).

pub mod detect;
pub mod ensembl;
pub mod hgvs;
pub mod id;
pub mod region;
pub mod spdi;

use std::collections::HashSet;
use std::fmt;
use std::io::{self, BufRead, Read};
use std::sync::Arc;

use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::transcript::Transcript;
use vep_core::variant::{classify_variant, InputVariant};

use crate::error::IoError;

/// The longest line a reader of a variant input accepts, 64 MiB. A record of
/// a few megabytes is legitimate (a long INFO column, thousands of sample
/// columns, a structural variant written out base by base) and sits far below
/// the cap; a file with no line end in it is refused at the cap instead of
/// being read whole into one line buffer.
pub const MAX_LINE_BYTES: usize = 64 << 20;

/// The error a bounded line read returns for a line past [`MAX_LINE_BYTES`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LineTooLong {
    /// The line's number in the input, counted from 1.
    pub line: u64,
}

impl fmt::Display for LineTooLong {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "line {} is longer than the {MAX_LINE_BYTES}-byte limit for one line of input",
            self.line
        )
    }
}

impl std::error::Error for LineTooLong {}

impl LineTooLong {
    /// `self` as the `InvalidData` I/O error a reader returns.
    pub fn into_io_error(self) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, self)
    }

    /// The [`LineTooLong`] an I/O error carries, when it carries one.
    pub fn from_io_error(error: &io::Error) -> Option<&LineTooLong> {
        error.get_ref()?.downcast_ref()
    }
}

/// `BufRead::read_until(delimiter, buf)` with the line's length capped at
/// [`MAX_LINE_BYTES`]: appends line `line_number` through its delimiter, or to
/// the end of the input, and returns the bytes appended (0 at the end of the
/// input). A line still open past the cap is a [`LineTooLong`] I/O error, the
/// reader left inside the line and the bytes read so far left in `buf`.
pub fn read_until_bounded<R: BufRead + ?Sized>(
    reader: &mut R,
    delimiter: u8,
    buf: &mut Vec<u8>,
    line_number: u64,
) -> io::Result<usize> {
    // The cap plus the delimiter: a line of exactly the cap is read whole.
    let limit = MAX_LINE_BYTES as u64 + 1;
    let read = reader.take(limit).read_until(delimiter, buf)?;
    if read as u64 == limit && buf.last() != Some(&delimiter) {
        return Err(LineTooLong { line: line_number }.into_io_error());
    }
    Ok(read)
}

/// Trait for parsing input variants from various formats.
pub trait InputParser: Send {
    /// Read the next variant from the input.
    /// Returns `None` at end of input.
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>>;

    /// The empty and whitespace-only lines read past since the last call, for
    /// the caller to warn about as Ensembl VEP's `skip_empty_lines` does
    /// (Parser.pm:266-293). A parser whose format has no such lines, or that
    /// warns about them itself through [`InputParser::take_warnings`], returns 0.
    fn take_empty_lines(&mut self) -> u64 {
        0
    }

    /// The warnings the parser raised since the last call, in the order they
    /// arose and without the `WARNING: ` prefix, for the caller to write through
    /// its warning sink: each skips no record, as Ensembl VEP's `warning_msg`
    /// skips none. A parser that warns about nothing returns none.
    fn take_warnings(&mut self) -> Vec<String> {
        Vec::new()
    }

    /// Ensembl VEP's line number for the record last yielded, from a parser
    /// whose format counts its lines other than by the records read: the
    /// region, Ensembl and HGVS readers count every line of the input, empty
    /// lines and `#` lines included (`LineCount`), and the identifier parser
    /// counts every line of its list, names the cache lacks included, and
    /// gives every mapping of a name the name's line (Parser/ID.pm:148,
    /// Parser.pm:275). `None` from a parser whose records the caller numbers
    /// by its own count of them.
    fn record_line_number(&self) -> Option<u64> {
        None
    }
}

/// The warning `skip_empty_lines` writes for a run of empty lines
/// (Parser.pm:280-288): `skipped` lines read past, the first at `from_line`.
pub fn empty_lines_warning(skipped: u64, from_line: u64) -> String {
    let plural = if skipped > 1 { "s" } else { "" };
    format!("Skipped {skipped} empty line{plural} from line {from_line}")
}

/// Perl's `/\w+/` on `text`: whether it holds an ASCII letter, digit or `_`.
/// `skip_empty_lines` reads past a record without one (Parser.pm:273).
pub(crate) fn has_word_character(text: &str) -> bool {
    text.chars().any(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// Ensembl VEP's line count over a list-based input, kept by the parser that
/// reads it. `line_number` advances once per record (Parser/Region.pm:167 and
/// its twins in `Parser/HGVS.pm:151` and `Parser/VEP_input.pm:119`) and once
/// per empty line read past (Parser.pm:275), so every line of the input
/// counts, a `#` line included: ensembl-io's list parser has no metadata lines
/// (ListBasedParser.pm:60, VEP_input.pm:61), so such a line is a record that
/// builds no feature. A run of empty lines is warned about in
/// `skip_empty_lines`'s words when the run ends, at the next line that counts
/// as a record or at the end of the input (Parser.pm:273-290), the warning
/// naming the run's first line.
#[derive(Debug, Default)]
pub(crate) struct LineCount {
    /// Lines read so far: the number of the last line read.
    line: u64,
    /// The run of empty lines the last lines read make up.
    empty_run: u64,
    /// The empty-line warnings since the last `take_warnings`, in the order the
    /// runs ended.
    warnings: Vec<String>,
}

impl LineCount {
    /// An empty line read past.
    pub(crate) fn empty_line(&mut self) {
        self.line += 1;
        self.empty_run += 1;
    }

    /// A line read as a record, a `#` line included: ends any run of empty
    /// lines before it.
    pub(crate) fn record_line(&mut self) {
        self.end_empty_run();
        self.line += 1;
    }

    /// The end of the input: ends a run of empty lines the input closes on.
    pub(crate) fn end_of_input(&mut self) {
        self.end_empty_run();
    }

    /// Warns about the run of empty lines the last lines read make up, when
    /// there is one, naming the run's first line.
    fn end_empty_run(&mut self) {
        if self.empty_run == 0 {
            return;
        }
        let run = std::mem::take(&mut self.empty_run);
        self.warnings
            .push(empty_lines_warning(run, self.line - run + 1));
    }

    /// The number of the last line read, 0 before the first.
    pub(crate) fn line(&self) -> u64 {
        self.line
    }

    /// The empty-line warnings since the last call.
    pub(crate) fn take_warnings(&mut self) -> Vec<String> {
        std::mem::take(&mut self.warnings)
    }
}

/// A reference genome a parser reads bases from: the run's `--fasta`, which a
/// notation that leaves its reference allele implicit is completed from.
pub trait ReferenceSequence: Send + Sync {
    /// The bases of `chr` from `start` to `end`, 1-based and inclusive; `None`
    /// when the reference has no such chromosome or span.
    fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>>;

    /// The name the reference gives the sequence `chr` reaches: `chr` itself
    /// when a sequence carries it, else the one its `chr` spelling or a synonym
    /// names; `None` when no sequence is reached. A notation that names a
    /// chromosome by an accession is printed under this name, as Ensembl VEP
    /// prints a slice's `seq_region_name`.
    fn record_name(&self, chr: &str) -> Option<String>;

    /// The number of bases of the sequence `chr` reaches ([`Self::record_name`]),
    /// the last position [`Self::sequence`] can read; `None` when no sequence is
    /// reached.
    fn sequence_length(&self, chr: &str) -> Option<u64>;
}

/// The transcripts an HGVS `c.`, `n.` or `p.` notation is resolved against:
/// the run's transcript cache, standing in for the core database Ensembl VEP's
/// `fetch_by_hgvs_notation` queries (VariationFeatureAdaptor.pm:1829-1849,
/// 1908-1941). A lookup that cannot read the cache reports why.
pub trait TranscriptSource: Send + Sync {
    /// `TranscriptAdaptor::fetch_by_stable_id` (TranscriptAdaptor.pm:164-181):
    /// the transcript whose stable id is `stable_id`, else, when `stable_id`
    /// ends in `.<version>`, the one with the id before the dot and that
    /// version; `None` when the cache has neither.
    fn by_stable_id(&self, stable_id: &str) -> Result<Option<Transcript>, String>;

    /// `TranscriptAdaptor::fetch_by_translation_stable_id`
    /// (TranscriptAdaptor.pm:307-332): the transcript whose translation's stable
    /// id is `stable_id`, else, when `stable_id` ends in `.<version>`, the one
    /// whose translation has the id before the dot and that version; `None`
    /// when the cache has neither.
    fn by_translation_id(&self, stable_id: &str) -> Result<Option<Transcript>, String>;

    /// `_get_gene_transcripts` (VariationFeatureAdaptor.pm:1993-2012) in each
    /// core group the run reads, in the run's group order (Parser/HGVS.pm:164:
    /// the Ensembl transcripts first, the RefSeq ones first under `--refseq`):
    /// one entry per group, holding every transcript of the group's reference
    /// gene named `name` in database order (the order `fetch_all_by_Gene`
    /// returns them, which the pick's stable sort keeps among equals), empty
    /// for a group that carries no gene of the name.
    fn gene_transcripts(&self, name: &str) -> Result<Vec<Vec<Transcript>>, String>;
}

/// What the text-format parsers are built against beyond their reader
/// (Parser.pm:128-135 and :159 hand the format's parser the delimiter; the
/// reference, the synonyms and the valid chromosomes, `Runner.pm:827`, reach it
/// through the run's configuration).
#[derive(Clone, Default)]
pub struct ParserContext {
    /// The pattern a column-based input is split on.
    pub delimiter: Delimiter,
    /// The run's reference FASTA, when one is open.
    pub fasta: Option<Arc<dyn ReferenceSequence>>,
    /// The chromosome synonyms an input name is resolved through.
    pub synonyms: Arc<ChromosomeSynonyms>,
    /// The chromosome names the run's annotation source holds, for a parser
    /// that names a record's chromosome as the source names it (an accession
    /// resolves through the synonyms to one of these).
    pub chromosomes: Arc<HashSet<String>>,
    /// `--max_sv_size`, the span `validate_vf` refuses a record beyond
    /// (Parser.pm:491-499); `None` for the flag's `-1`, no limit.
    pub max_sv_size: Option<u64>,
    /// `--minimal`: `post_process_vfs` minimises every bi-allelic record
    /// (Parser.pm:852) before `validate_vf` measures, names and upper-cases it,
    /// so a parser that measures or names a record as written does so on the
    /// trimmed alleles under the flag.
    pub minimal: bool,
    /// The run's transcripts, for a parser that names a position on a
    /// transcript or its protein (HGVS `c.`, `n.` and `p.`); `None` when the run
    /// annotates without a transcript cache.
    pub transcripts: Option<Arc<dyn TranscriptSource>>,
}

impl ParserContext {
    /// A context with `delimiter` alone.
    pub fn with_delimiter(delimiter: Delimiter) -> Self {
        Self {
            delimiter,
            ..Self::default()
        }
    }
}

/// The pattern a column-based input is split on, as `--delimiter` names it
/// (Config.pm:82, a tab by default at :302) or as format detection chooses it
/// from the first data line (Parser.pm:387-415). Ensembl VEP hands the value to
/// Perl's `split` as a regular expression (ensembl-io ColumnBasedParser.pm:250);
/// the three patterns its detection chooses, `\s+` and text with no
/// regular-expression character are split here directly, and any other value
/// as the regular expression it is, by a splitter the caller compiles
/// ([`Delimiter::from_flag`]).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Delimiter {
    /// A tab: the default, and `detect_format`'s choice for a line with one.
    #[default]
    Tab,
    /// A run of spaces (`" +"`): `detect_format`'s choice for a line with two
    /// spaces in a row and no tab. A leading run yields a leading empty field.
    SpaceRun,
    /// Whitespace in the manner of Perl's `split ' '` (`" "`): `detect_format`'s
    /// choice for a line with a single space and no tab. Leading whitespace is
    /// ignored.
    Whitespace,
    /// A run of whitespace (`\s+`): as [`Delimiter::SpaceRun`] with tabs and the
    /// other blank characters as separators too.
    WhitespaceRun,
    /// Text with no regular-expression character, matched as written.
    Literal(String),
    /// Any other text, the regular expression Perl's `split` reads it as.
    Pattern(PatternDelimiter),
}

/// How a regular-expression delimiter splits a line, as Perl's `split` splits
/// on a pattern: the fields are the text between matches, a leading empty field
/// kept when a match opens the line and trailing empty fields dropped; a match
/// never ends at the field's start, so an empty match there is passed over for
/// the next match one character on; each capture group of a match is a field
/// of its own, an unmatched group an empty one. A splitter of a `str` ends its
/// fields on character boundaries.
pub trait PatternSplit: Send + Sync {
    /// Whether the pattern matches anywhere in `line`.
    fn is_match(&self, line: &[u8]) -> bool;
    /// The fields of a line of text.
    fn split_str<'a>(&self, line: &'a str) -> Vec<&'a str>;
    /// The fields of a line of bytes, which need not be UTF-8.
    fn split_bytes<'a>(&self, line: &'a [u8]) -> Vec<&'a [u8]>;
}

/// A regular-expression delimiter: the `--delimiter` text and the splitter
/// compiled from it. Two are equal when their texts are.
#[derive(Clone)]
pub struct PatternDelimiter {
    text: String,
    split: Arc<dyn PatternSplit>,
}

impl PatternDelimiter {
    /// The delimiter of `text`, split by `split`.
    pub fn new(text: &str, split: Arc<dyn PatternSplit>) -> Self {
        Self {
            text: text.to_string(),
            split,
        }
    }

    /// The `--delimiter` text the splitter was compiled from.
    pub fn text(&self) -> &str {
        &self.text
    }
}

impl std::fmt::Debug for PatternDelimiter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("PatternDelimiter").field(&self.text).finish()
    }
}

impl PartialEq for PatternDelimiter {
    fn eq(&self, other: &Self) -> bool {
        self.text == other.text
    }
}

impl Eq for PatternDelimiter {}

/// The characters Perl's `\s` and `split ' '` separate on: ASCII whitespace
/// plus the vertical tab, which `char::is_ascii_whitespace` leaves out.
pub(crate) fn is_perl_space(c: char) -> bool {
    c.is_ascii_whitespace() || c == '\u{0B}'
}

/// The characters a regular expression reads as other than themselves.
fn is_regex_metacharacter(c: char) -> bool {
    matches!(
        c,
        '\\' | '^' | '$' | '.' | '|' | '?' | '*' | '+' | '(' | ')' | '[' | ']' | '{' | '}'
    )
}

impl Delimiter {
    /// The delimiter a `--delimiter` value names: a tab character or the text
    /// `\t`, a space, a space and a plus, `\s+`, text with no
    /// regular-expression character, or any other text as the regular
    /// expression `compile` makes of it, whose error is this one's. An empty
    /// value is an error.
    pub fn from_flag(
        value: &str,
        compile: impl FnOnce(&str) -> Result<Arc<dyn PatternSplit>, String>,
    ) -> Result<Self, String> {
        match value {
            "\t" | "\\t" => Ok(Delimiter::Tab),
            " +" => Ok(Delimiter::SpaceRun),
            " " => Ok(Delimiter::Whitespace),
            "\\s+" => Ok(Delimiter::WhitespaceRun),
            "" => Err("--delimiter needs at least one character".to_string()),
            pattern if pattern.chars().any(is_regex_metacharacter) => Ok(Delimiter::Pattern(
                PatternDelimiter::new(pattern, compile(pattern)?),
            )),
            literal => Ok(Delimiter::Literal(literal.to_string())),
        }
    }

    /// The delimiter for `line` during format detection (Parser.pm:401-412):
    /// this one when the line contains it, else a run of spaces when the line
    /// has two in a row, else whitespace when it has a space, else this one
    /// still.
    pub fn chosen_for(&self, line: &str) -> Self {
        let present = match self {
            Delimiter::Tab => line.contains('\t'),
            Delimiter::SpaceRun | Delimiter::Whitespace => line.contains(' '),
            Delimiter::WhitespaceRun => line.contains(is_perl_space),
            Delimiter::Literal(text) => line.contains(text.as_str()),
            Delimiter::Pattern(pattern) => pattern.split.is_match(line.as_bytes()),
        };
        if present {
            self.clone()
        } else if line.contains("  ") {
            Delimiter::SpaceRun
        } else if line.contains(' ') {
            Delimiter::Whitespace
        } else {
            self.clone()
        }
    }

    /// The fields of `line` as Perl's `split` yields them: a leading empty
    /// field kept, trailing empty fields dropped, and under
    /// [`Delimiter::Whitespace`] the leading whitespace ignored.
    pub fn split<'a>(&self, line: &'a str) -> Vec<&'a str> {
        let mut fields: Vec<&str> = match self {
            Delimiter::Pattern(pattern) => pattern.split.split_str(line),
            _ => self
                .split_bytes(line.as_bytes())
                .into_iter()
                .map(|field| {
                    // A separator is an ASCII byte or a whole UTF-8 literal, so a
                    // field ends on character bounds.
                    std::str::from_utf8(field).expect("a field of a str is a str")
                })
                .collect(),
        };
        drop_trailing_empty_fields(&mut fields);
        fields
    }

    /// [`Self::split`] on the bytes of a line, for a line that may not be
    /// UTF-8: the fields are byte slices of `line`, and a byte outside ASCII
    /// is a separator only inside a literal delimiter, matched byte for byte,
    /// or where a regular expression matches it.
    pub fn split_bytes<'a>(&self, line: &'a [u8]) -> Vec<&'a [u8]> {
        let mut fields: Vec<&[u8]> = match self {
            Delimiter::Tab => line.split(|b| *b == b'\t').collect(),
            Delimiter::SpaceRun => split_on_runs(line, |b| b == b' '),
            Delimiter::WhitespaceRun => split_on_runs(line, is_perl_space_byte),
            Delimiter::Whitespace => line
                .split(|b| is_perl_space_byte(*b))
                .filter(|field| !field.is_empty())
                .collect(),
            Delimiter::Literal(text) => split_on_literal(line, text.as_bytes()),
            Delimiter::Pattern(pattern) => pattern.split.split_bytes(line),
        };
        drop_trailing_empty_fields(&mut fields);
        fields
    }

    /// The text the fields of an input line are joined with where Ensembl VEP
    /// re-emits the line (`OutputFactory/JSON.pm:172, 247`): the delimiter's
    /// text, or a single space for a delimiter with a `+` in it (`" +"`, `\s+`,
    /// a regular expression with one), which `:172` substitutes before the join.
    pub fn join_text(&self) -> &str {
        match self {
            Delimiter::Tab => "\t",
            Delimiter::SpaceRun | Delimiter::Whitespace | Delimiter::WhitespaceRun => " ",
            Delimiter::Literal(text) => text,
            Delimiter::Pattern(pattern) if pattern.text().contains('+') => " ",
            Delimiter::Pattern(pattern) => pattern.text(),
        }
    }
}

/// [`is_perl_space`] on a byte.
fn is_perl_space_byte(b: u8) -> bool {
    b.is_ascii_whitespace() || b == 0x0B
}

/// Perl's `split / +/` (or `/\s+/`): each run of separator characters is one
/// separator, so a leading run yields a leading empty field.
fn split_on_runs(line: &[u8], is_separator: impl Fn(u8) -> bool) -> Vec<&[u8]> {
    let mut fields: Vec<&[u8]> = Vec::new();
    if line.first().is_some_and(|b| is_separator(*b)) {
        fields.push(b"");
    }
    fields.extend(
        line.split(|b| is_separator(*b))
            .filter(|field| !field.is_empty()),
    );
    fields
}

/// Perl's `split` on literal text: the fields between each occurrence of
/// `text`, empty fields kept.
fn split_on_literal<'a>(line: &'a [u8], text: &[u8]) -> Vec<&'a [u8]> {
    if text.is_empty() {
        return vec![line];
    }
    let mut fields = Vec::new();
    let mut rest = line;
    while let Some(at) = rest.windows(text.len()).position(|window| window == text) {
        fields.push(&rest[..at]);
        rest = &rest[at + text.len()..];
    }
    fields.push(rest);
    fields
}

/// Removes the empty fields at the end of a line, which Perl's `split` drops.
fn drop_trailing_empty_fields<T: AsRef<[u8]>>(fields: &mut Vec<T>) {
    while fields.last().is_some_and(|field| field.as_ref().is_empty()) {
        fields.pop();
    }
}

/// Ensembl VEP's `post_process_vfs` on a parsed sequence variant (Parser.pm:
/// 862-881): a bi-allelic record whose allele string has no `-` and whose
/// alleles differ in length is minimised by `minimise_alleles` (:964-1019)
/// through `trim_sequences` (ensembl-variation Utils/Sequence.pm:965-1035):
/// the bases the two alleles share are trimmed from the left, then from the
/// right, the start and end moving in by the bases trimmed on their side, and
/// an emptied allele is written `-`. The bases are compared as written (`eq`),
/// since `validate_vf` upper-cases the string only afterwards (Parser.pm:554),
/// so a mixed-case pair is not trimmed. The record keeps its written alleles in
/// `original_allele_string` and is marked minimised. A multi-allelic record is
/// left as written: Ensembl VEP splits and trims it only under `--minimal`, in
/// its input buffer. Returns whether the alleles changed.
pub fn minimise_alleles(variant: &mut InputVariant) -> bool {
    if variant.is_structural
        || variant.alt_alleles.len() != 1
        || variant.record_allele_string_multi.is_some()
        || variant.allele_string.contains('-')
        || variant.ref_allele.len() == variant.alt_alleles[0].len()
    {
        return false;
    }
    let (reference, alt) = (&variant.ref_allele, &variant.alt_alleles[0]);
    let prefix_len = reference
        .iter()
        .zip(alt.iter())
        .take_while(|(r, a)| r == a)
        .count();
    let suffix_len = reference[prefix_len..]
        .iter()
        .rev()
        .zip(alt[prefix_len..].iter().rev())
        .take_while(|(r, a)| r == a)
        .count();
    if prefix_len == 0 && suffix_len == 0 {
        return false;
    }
    let trimmed = |allele: &[u8]| match &allele[prefix_len..allele.len() - suffix_len] {
        [] => b"-".to_vec(),
        kept => kept.to_vec(),
    };
    let (min_ref, min_alt) = (trimmed(reference), trimmed(alt));
    variant.original_allele_string = Some(variant.allele_string.clone());
    variant.allele_string = format!(
        "{}/{}",
        String::from_utf8_lossy(&min_ref),
        String::from_utf8_lossy(&min_alt),
    );
    variant.start += prefix_len as u64;
    variant.end = variant.end.saturating_sub(suffix_len as u64);
    variant.variant_class = classify_variant(&min_ref, &min_alt);
    variant.ref_allele = min_ref;
    variant.alt_alleles = vec![min_alt];
    variant.minimised = true;
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A stand-in for a compiled pattern: `,` and `;` separate, every field
    /// kept, as the plumbing around a [`PatternDelimiter`] needs no more.
    struct CommaOrSemicolon;

    impl PatternSplit for CommaOrSemicolon {
        fn is_match(&self, line: &[u8]) -> bool {
            line.iter().any(|b| matches!(b, b',' | b';'))
        }

        fn split_str<'a>(&self, line: &'a str) -> Vec<&'a str> {
            line.split([',', ';']).collect()
        }

        fn split_bytes<'a>(&self, line: &'a [u8]) -> Vec<&'a [u8]> {
            line.split(|b| matches!(b, b',' | b';')).collect()
        }
    }

    fn compiled(text: &str) -> Result<Arc<dyn PatternSplit>, String> {
        if text == "[" {
            return Err(format!("unclosed {text}"));
        }
        Ok(Arc::new(CommaOrSemicolon))
    }

    fn pattern(text: &str) -> Delimiter {
        Delimiter::Pattern(PatternDelimiter::new(text, Arc::new(CommaOrSemicolon)))
    }

    /// The three patterns format detection chooses, `\s+`, a literal, an empty
    /// value and a regular expression, read as `--delimiter` reads them: text
    /// with a regular-expression character is handed to the compiler, whose
    /// error is the flag's.
    #[test]
    fn delimiter_flag_values() {
        let from_flag = |value: &str| Delimiter::from_flag(value, compiled);
        assert_eq!(from_flag("\t"), Ok(Delimiter::Tab));
        assert_eq!(from_flag("\\t"), Ok(Delimiter::Tab));
        assert_eq!(from_flag(" "), Ok(Delimiter::Whitespace));
        assert_eq!(from_flag(" +"), Ok(Delimiter::SpaceRun));
        assert_eq!(from_flag("\\s+"), Ok(Delimiter::WhitespaceRun));
        assert_eq!(from_flag(","), Ok(Delimiter::Literal(",".to_string())));
        assert_eq!(from_flag("::"), Ok(Delimiter::Literal("::".to_string())));
        assert!(from_flag("").unwrap_err().contains("at least one"));
        for regex in ["|", ".", "\\s", "[ ,]", "a+"] {
            assert_eq!(from_flag(regex), Ok(pattern(regex)), "{regex:?}");
        }
        assert_eq!(from_flag("["), Err("unclosed [".to_string()));
    }

    /// A pattern delimiter splits a line of text or of bytes by its splitter,
    /// trailing empty fields dropped; it is present in a line the pattern
    /// matches and chosen over the space patterns there; its join text is the
    /// `--delimiter` text when that has no `+`; and two pattern delimiters are
    /// equal when their texts are, whatever splitter each carries.
    #[test]
    fn pattern_delimiter_splits_by_its_splitter() {
        let delimiter = pattern("[,;]");
        assert_eq!(delimiter.split("a,b;c,,"), ["a", "b", "c"]);
        assert_eq!(
            delimiter.split_bytes(b",caf\xe9;b,,"),
            [&b""[..], b"caf\xe9", b"b"]
        );
        assert_eq!(delimiter.chosen_for("a,b"), delimiter);
        assert_eq!(delimiter.chosen_for("a  b"), Delimiter::SpaceRun);
        assert_eq!(delimiter.chosen_for("a b"), Delimiter::Whitespace);
        assert_eq!(delimiter.chosen_for("ab"), delimiter);
        assert_eq!(delimiter.join_text(), "[,;]");
        assert_eq!(
            format!("{delimiter:?}"),
            r#"Pattern(PatternDelimiter("[,;]"))"#
        );
        assert_ne!(delimiter, pattern("[;,]"));
        assert_ne!(delimiter, Delimiter::Literal("[,;]".to_string()));
    }

    /// Perl's `split` on each pattern: a leading empty field is kept under a tab,
    /// a literal and a run pattern, never under `split ' '`; trailing empty
    /// fields go; `\s+` separates on tabs where `" +"` does not.
    #[test]
    fn delimiter_split_follows_perl() {
        assert_eq!(Delimiter::Tab.split("a\tb c\t\t"), ["a", "b c"]);
        assert_eq!(Delimiter::Tab.split("\ta\tb"), ["", "a", "b"]);
        assert_eq!(Delimiter::SpaceRun.split(" a  b\tc "), ["", "a", "b\tc"]);
        assert_eq!(Delimiter::Whitespace.split("  a \tb \x0B"), ["a", "b"]);
        assert_eq!(Delimiter::WhitespaceRun.split(" a \tb "), ["", "a", "b"]);
        assert_eq!(
            Delimiter::Literal(",".to_string()).split(",a,,b,,"),
            ["", "a", "", "b"]
        );
        assert!(Delimiter::Tab.split("").is_empty());
        assert!(Delimiter::Whitespace.split("   ").is_empty());
    }

    /// The byte split gives the same fields as the text split, and a byte
    /// outside ASCII (`0xE9`) is kept whole in its field under every pattern.
    #[test]
    fn delimiter_split_bytes_keeps_bytes_outside_ascii() {
        let line = b"21 30282169 caf\xe9 C A . . .";
        let fields = Delimiter::Whitespace.split_bytes(line);
        assert_eq!(fields[2], b"caf\xe9");
        assert_eq!(fields.len(), 8);
        assert_eq!(
            Delimiter::SpaceRun.split_bytes(b" a  caf\xe9 "),
            [&b""[..], b"a", b"caf\xe9"]
        );
        assert_eq!(
            Delimiter::Literal(",".to_string()).split_bytes(b",caf\xe9,,b,,"),
            [&b""[..], b"caf\xe9", b"", b"b"]
        );
        for delimiter in [
            Delimiter::Tab,
            Delimiter::SpaceRun,
            Delimiter::Whitespace,
            Delimiter::WhitespaceRun,
            Delimiter::Literal("::".to_string()),
        ] {
            let text = " a::b\tc  d \x0Be ";
            let from_text: Vec<&[u8]> =
                delimiter.split(text).iter().map(|f| f.as_bytes()).collect();
            assert_eq!(
                delimiter.split_bytes(text.as_bytes()),
                from_text,
                "{delimiter:?}"
            );
        }
    }

    /// A line of exactly the cap is read whole with its line end; one byte
    /// longer is the `LineTooLong` error naming the line, with the reader left
    /// inside the line, and the lines before and after it read as usual. The
    /// last line of an input needs no line end.
    #[test]
    fn read_until_bounded_refuses_a_line_past_the_cap() {
        let mut input = b"first\n".to_vec();
        input.extend(std::iter::repeat_n(b'a', MAX_LINE_BYTES));
        input.push(b'\n');
        input.extend(std::iter::repeat_n(b'b', MAX_LINE_BYTES + 1));
        input.extend_from_slice(b"\nlast");
        let mut reader = std::io::Cursor::new(input);
        let mut buf = Vec::new();
        assert_eq!(
            read_until_bounded(&mut reader, b'\n', &mut buf, 1).unwrap(),
            6
        );
        assert_eq!(buf, b"first\n");
        buf.clear();
        assert_eq!(
            read_until_bounded(&mut reader, b'\n', &mut buf, 2).unwrap(),
            MAX_LINE_BYTES + 1
        );
        assert_eq!(buf.len(), MAX_LINE_BYTES + 1);
        buf.clear();
        let error = read_until_bounded(&mut reader, b'\n', &mut buf, 3).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(
            LineTooLong::from_io_error(&error),
            Some(&LineTooLong { line: 3 })
        );
        assert_eq!(
            error.to_string(),
            format!("line 3 is longer than the {MAX_LINE_BYTES}-byte limit for one line of input")
        );
        assert_eq!(buf.len(), MAX_LINE_BYTES + 1);
        // The rest of the refused line, then the last line without a line end.
        buf.clear();
        assert_eq!(
            read_until_bounded(&mut reader, b'\n', &mut buf, 3).unwrap(),
            1
        );
        buf.clear();
        assert_eq!(
            read_until_bounded(&mut reader, b'\n', &mut buf, 4).unwrap(),
            4
        );
        assert_eq!(buf, b"last");
        assert_eq!(
            read_until_bounded(&mut reader, b'\n', &mut buf, 5).unwrap(),
            0
        );
    }

    /// `OutputFactory/JSON.pm:172` joins with a single space under any
    /// delimiter with a `+` in it, so a regular expression with one (`,+`, or
    /// `\+` for a literal plus) joins as the space runs do, and one without
    /// joins with its text.
    #[test]
    fn delimiter_join_text_is_one_space_for_either_space_pattern() {
        assert_eq!(Delimiter::Tab.join_text(), "\t");
        assert_eq!(Delimiter::SpaceRun.join_text(), " ");
        assert_eq!(Delimiter::Whitespace.join_text(), " ");
        assert_eq!(Delimiter::WhitespaceRun.join_text(), " ");
        assert_eq!(Delimiter::Literal(";".to_string()).join_text(), ";");
        assert_eq!(pattern("[,;]").join_text(), "[,;]");
        assert_eq!(pattern(",+").join_text(), " ");
        assert_eq!(pattern("\\+").join_text(), " ");
    }

    fn variant(start: u64, end: u64, alleles: &str) -> InputVariant {
        let (reference, alt) = alleles.split_once('/').unwrap();
        InputVariant::new(
            "21".to_string(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    /// `trim_sequences` with `empty_to_dash`: `AC/A` at 100-101 is `C/-` at 101,
    /// `A/AT` at 100 is `-/T` at 101-100 (an insertion's end precedes its start),
    /// `TAA/GCTCA` at 100-102 trims its shared last base alone, and `ACG/AT`
    /// trims its shared first base and keeps `CG/T` at 101-102; each keeps the
    /// written string and is marked minimised.
    #[test]
    fn minimise_alleles_trims_as_trim_sequences_does() {
        let cases = [
            ((100, 101, "AC/A"), (101, 101, "C/-")),
            ((100, 100, "A/AT"), (101, 100, "-/T")),
            ((100, 102, "TAA/GCTCA"), (100, 101, "TA/GCTC")),
            ((100, 102, "ACG/AT"), (101, 102, "CG/T")),
        ];
        for ((start, end, written), (min_start, min_end, minimised)) in cases {
            let mut v = variant(start, end, written);
            assert!(minimise_alleles(&mut v), "{written}");
            assert_eq!((v.start, v.end), (min_start, min_end), "{written}");
            assert_eq!(v.allele_string, minimised, "{written}");
            assert_eq!(v.original_allele_string.as_deref(), Some(written));
            assert!(v.minimised, "{written}");
        }
    }

    /// Nothing is minimised when the alleles share no base at either end, when
    /// they are of one length, when either is `-`, when the record is
    /// multi-allelic (Ensembl VEP leaves those to `--minimal`) or structural.
    #[test]
    fn minimise_alleles_leaves_the_other_shapes_alone() {
        for (start, end, written) in [
            (100, 101, "AC/G"),
            (100, 100, "A/G"),
            (100, 101, "AC/GT"),
            (100, 101, "AC/-"),
            (101, 100, "-/A"),
        ] {
            let mut v = variant(start, end, written);
            assert!(!minimise_alleles(&mut v), "{written}");
            assert_eq!(v.allele_string, written);
            assert_eq!((v.start, v.end), (start, end));
            assert!(!v.minimised);
        }
        let mut multi = variant(100, 101, "AC/A");
        multi.record_allele_string_multi = Some("AC/A/G".into());
        assert!(!minimise_alleles(&mut multi));
        let mut structural = variant(100, 101, "N/DUP");
        structural.is_structural = true;
        assert!(!minimise_alleles(&mut structural));
    }

    /// `trim_sequences` compares the bases as written and `validate_vf` upper-cases
    /// them only afterwards, so `Ac/a` shares no base at either end and stays as it
    /// is, while `ac/a` is `c/-` at 101; a `-` anywhere in the string, `A-/A`
    /// included, leaves the record alone (Parser.pm:871-873).
    #[test]
    fn minimise_alleles_compares_the_written_case_and_skips_a_dash() {
        let mut mixed = variant(100, 101, "Ac/a");
        assert!(!minimise_alleles(&mut mixed));
        assert_eq!(mixed.allele_string, "Ac/a");
        let mut lower = variant(100, 101, "ac/a");
        assert!(minimise_alleles(&mut lower));
        assert_eq!(
            (lower.start, lower.end, lower.allele_string.as_str()),
            (101, 101, "c/-")
        );
        let mut dashed = variant(100, 101, "A-/A");
        assert!(!minimise_alleles(&mut dashed));
        assert_eq!((dashed.start, dashed.end), (100, 101));
    }
}
