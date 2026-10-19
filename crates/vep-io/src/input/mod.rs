// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Input parsers for variant formats, the detection of a format from the
//! first data line ([`detect`]), and what the parsers share: the field
//! delimiter of a column-based input ([`Delimiter`]), the run's reference
//! sequence, chromosome synonyms and chromosome names ([`ParserContext`]), and
//! the minimisation Ensembl VEP applies to every parsed record
//! ([`minimise_alleles`]).

pub mod detect;
pub mod ensembl;
pub mod hgvs;
pub mod id;
pub mod region;
pub mod spdi;

use std::collections::HashSet;
use std::sync::Arc;

use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{classify_variant, InputVariant};

use crate::error::IoError;

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
    /// identifier parser counts every line of its list, empty lines and names
    /// the cache lacks included, and gives every mapping of a name the name's
    /// line (Parser/ID.pm:148, Parser.pm:275). `None` from a parser whose
    /// records the caller numbers by its own count of them.
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
/// Perl's `split` as a regular expression; the three patterns its detection
/// chooses and `\s+` are read as that, and any other value as literal text, so
/// a value with a regular-expression character in it is refused rather than
/// split differently ([`Delimiter::from_flag`]).
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
    /// Any other text, matched as written.
    Literal(String),
}

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
    /// `\t`, a space, a space and a plus, `\s+`, or literal text. Any other
    /// value with a regular-expression character, and an empty value, is an
    /// error naming what is accepted.
    pub fn from_flag(value: &str) -> Result<Self, String> {
        match value {
            "\t" | "\\t" => Ok(Delimiter::Tab),
            " +" => Ok(Delimiter::SpaceRun),
            " " => Ok(Delimiter::Whitespace),
            "\\s+" => Ok(Delimiter::WhitespaceRun),
            "" => Err("--delimiter needs at least one character".to_string()),
            literal if literal.chars().any(is_regex_metacharacter) => Err(format!(
                "--delimiter {literal:?} is read as a regular expression by Ensembl VEP; \
                 a tab (\\t), a space, a run of spaces ( +), whitespace (\\s+) or literal \
                 text with no regular-expression character is accepted"
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
            Delimiter::Tab => line.split('\t').collect(),
            Delimiter::SpaceRun => split_on_runs(line, |c| c == ' '),
            Delimiter::WhitespaceRun => split_on_runs(line, is_perl_space),
            Delimiter::Whitespace => line
                .split(is_perl_space)
                .filter(|field| !field.is_empty())
                .collect(),
            Delimiter::Literal(text) => line.split(text.as_str()).collect(),
        };
        while fields.last().is_some_and(|field| field.is_empty()) {
            fields.pop();
        }
        fields
    }

    /// The text the fields of an input line are joined with where Ensembl VEP
    /// re-emits the line (`OutputFactory/JSON.pm:172, 247`): the delimiter
    /// itself, a single space for either space pattern.
    pub fn join_text(&self) -> &str {
        match self {
            Delimiter::Tab => "\t",
            Delimiter::SpaceRun | Delimiter::Whitespace | Delimiter::WhitespaceRun => " ",
            Delimiter::Literal(text) => text,
        }
    }
}

/// Perl's `split / +/` (or `/\s+/`): each run of separator characters is one
/// separator, so a leading run yields a leading empty field.
fn split_on_runs(line: &str, is_separator: impl Fn(char) -> bool + Copy) -> Vec<&str> {
    let mut fields = Vec::new();
    if line.starts_with(is_separator) {
        fields.push("");
    }
    fields.extend(line.split(is_separator).filter(|field| !field.is_empty()));
    fields
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

    /// The three patterns format detection chooses, `\s+`, a literal, an empty
    /// value and a regular expression, read as `--delimiter` reads them.
    #[test]
    fn delimiter_flag_values() {
        assert_eq!(Delimiter::from_flag("\t"), Ok(Delimiter::Tab));
        assert_eq!(Delimiter::from_flag("\\t"), Ok(Delimiter::Tab));
        assert_eq!(Delimiter::from_flag(" "), Ok(Delimiter::Whitespace));
        assert_eq!(Delimiter::from_flag(" +"), Ok(Delimiter::SpaceRun));
        assert_eq!(Delimiter::from_flag("\\s+"), Ok(Delimiter::WhitespaceRun));
        assert_eq!(
            Delimiter::from_flag(","),
            Ok(Delimiter::Literal(",".to_string()))
        );
        assert_eq!(
            Delimiter::from_flag("::"),
            Ok(Delimiter::Literal("::".to_string()))
        );
        assert!(Delimiter::from_flag("")
            .unwrap_err()
            .contains("at least one"));
        for regex in ["|", ".", "\\s", "[ ,]", "a+"] {
            let error = Delimiter::from_flag(regex).unwrap_err();
            assert!(error.contains("regular expression"), "{regex:?}: {error}");
        }
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

    #[test]
    fn delimiter_join_text_is_one_space_for_either_space_pattern() {
        assert_eq!(Delimiter::Tab.join_text(), "\t");
        assert_eq!(Delimiter::SpaceRun.join_text(), " ");
        assert_eq!(Delimiter::Whitespace.join_text(), " ");
        assert_eq!(Delimiter::WhitespaceRun.join_text(), " ");
        assert_eq!(Delimiter::Literal(";".to_string()).join_text(), ";");
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
