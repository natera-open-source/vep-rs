// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines `vep filter` reads, as Ensembl VEP's `filter_vep` reads them
//! (`filter_vep` at release/116.3): the header block of an annotated output,
//! the format of its data lines, and each line parsed into the hash a filter
//! is evaluated on.
//!
//! A line starting with `#` is a header line. Among them, a `##` line is a
//! meta line: the `##INFO=<ID=<field>,...>` line of the consequence field
//! names the consequence fields in its `Format:` list, every other
//! `##INFO=<ID=<name>,` line and every `## <name> : <description>` line names
//! a field `--list` prints (`parse_headers`, 351-381). The last `#` line that
//! is not a meta line is the column line. The fields a filter reads are the
//! consequence fields when the header block names them, else the columns
//! (216).
//!
//! A data line is a VCF line when its first word has a word character, its second is
//! a number, its fourth is bases, `N`, `-` or `.` and its fifth is true to Perl,
//! neither empty nor `0` (`detect_format`, 423-439); otherwise it is a tab line. A VCF
//! line's hash holds its columns and the INFO keys, and each entry of every
//! INFO item under the consequence field's name is one hash over those plus
//! the consequence fields (246-271). A tab line's hash holds its columns, with
//! an `Extra` column exploded into its `key=value` pairs (400-405). In every
//! hash an empty value is undefined, and so is `-` except under `Allele`
//! (407-411); every field the header names is a key, undefined when the line
//! lacks it (413-416); the `INFO` key is undefined everywhere (398).
//!
//! A line is read as bytes and written back as the bytes it was read from,
//! as `filter_vep` reads and writes it with no decoding layer (157-158, 184,
//! 319). A byte outside UTF-8 is one character in every comparison, as it is
//! one character of Perl's byte string: [`decode_line`] maps each such byte to
//! a private-use character and [`encode_line`] maps it back. A command-line
//! word is read the same way ([`decode_word`]), so a filter value, a field
//! name or a file name ([`path_of`]) compares byte for byte with the line, as
//! Ensembl VEP compares its `@ARGV` byte strings. The text rules that read a line
//! are Perl's over a byte string: a blank, a word character and a letter's
//! case are ASCII, and no byte above ASCII is any of them.

use std::borrow::Cow;
use std::collections::{BTreeSet, HashMap};
use std::ffi::OsStr;
use std::path::Path;
use std::sync::LazyLock;

use regex::Regex;

use crate::filter_set::FilterData;

/// The header block of an annotated output, parsed as `parse_headers` parses
/// it (351-381).
#[derive(Debug, Clone, Default)]
pub struct Headers {
    /// Every header line, as read.
    pub raw: Vec<String>,
    /// The consequence fields of the `##INFO` line naming the consequence
    /// field, in its `Format:` order; `None` without that line.
    pub consequence_fields: Option<Vec<String>>,
    /// The columns of the column line, in order.
    pub columns: Vec<String>,
    /// Every field name the header block names.
    pub allowed: BTreeSet<String>,
    column_index: HashMap<String, usize>,
    field_index: HashMap<String, usize>,
}

static INFO_ID: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"INFO=<ID=(.+?),").unwrap());
static FORMAT_LIST: LazyLock<Regex> = LazyLock::new(|| Regex::new(r#"Format: (.+?)""#).unwrap());
static DESCRIBED_FIELD: LazyLock<Regex> = LazyLock::new(|| Regex::new(r" (.+?) :").unwrap());

/// The first capture group of `regex`'s first match in `text`.
fn first_group<'a>(regex: &Regex, text: &'a str) -> Option<&'a str> {
    regex
        .captures(text)
        .and_then(|captures| captures.get(1))
        .map(|group| group.as_str())
}

/// The code point a byte of an input line that is not UTF-8 is read as, less
/// the byte: `0x80..=0xFF` land on U+10FF80..=U+10FFFF, the last 128 scalar
/// values of plane 16, private use but for the two noncharacters at its end.
const ESCAPE_BASE: u32 = 0x10FF00;

/// The character that stands for `byte` in a decoded line.
fn escaped(byte: u8) -> char {
    char::from_u32(ESCAPE_BASE + u32::from(byte)).expect("a private-use code point")
}

/// The byte `c` stands for, when it is one of the escape characters.
fn escaped_byte(c: char) -> Option<u8> {
    let code = u32::from(c);
    (code >= ESCAPE_BASE + 0x80).then(|| (code - ESCAPE_BASE) as u8)
}

/// Whether `bytes` holds the first byte of every escape character's UTF-8
/// form, so that a text without it needs no walk.
fn may_hold_escape(bytes: &[u8]) -> bool {
    bytes.contains(&0xF4)
}

/// The text of an input line read as bytes: the line itself when it is UTF-8,
/// else every byte of each sequence that is not UTF-8 read as the escape
/// character of that byte, one character per byte, as each is one character
/// of the byte string `filter_vep` reads. An escape character that stands
/// written in the line as UTF-8 is read as the four escape characters of its
/// own bytes, so that [`encode_line`] returns every line to the bytes it was
/// read from, whatever they are.
pub fn decode_line(bytes: &[u8]) -> Cow<'_, str> {
    match std::str::from_utf8(bytes) {
        Ok(text) if !may_hold_escape(bytes) => Cow::Borrowed(text),
        _ => Cow::Owned(decode_escaping(bytes)),
    }
}

fn decode_escaping(bytes: &[u8]) -> String {
    let mut text = String::with_capacity(bytes.len() + 16);
    for chunk in bytes.utf8_chunks() {
        let valid = chunk.valid();
        if may_hold_escape(valid.as_bytes()) {
            for c in valid.chars() {
                if escaped_byte(c).is_some() {
                    text.extend(c.encode_utf8(&mut [0; 4]).bytes().map(escaped));
                } else {
                    text.push(c);
                }
            }
        } else {
            text.push_str(valid);
        }
        text.extend(chunk.invalid().iter().copied().map(escaped));
    }
    text
}

/// The bytes a decoded line is written as: its UTF-8, with each escape
/// character written as the one byte it stands for, so that
/// `encode_line(&decode_line(bytes)) == bytes` for every `bytes`.
pub fn encode_line(text: &str) -> Cow<'_, [u8]> {
    if !may_hold_escape(text.as_bytes()) {
        return Cow::Borrowed(text.as_bytes());
    }
    let mut bytes = Vec::with_capacity(text.len());
    for c in text.chars() {
        match escaped_byte(c) {
            Some(byte) => bytes.push(byte),
            None => bytes.extend_from_slice(c.encode_utf8(&mut [0; 4]).as_bytes()),
        }
    }
    Cow::Owned(bytes)
}

/// The text of a command-line word read as the bytes the shell passed, as
/// Ensembl VEP reads `@ARGV`, each byte outside UTF-8 one escape character as in a
/// line ([`decode_line`]).
pub fn decode_word(word: &OsStr) -> String {
    #[cfg(unix)]
    let bytes = Cow::Borrowed(std::os::unix::ffi::OsStrExt::as_bytes(word));
    #[cfg(not(unix))]
    let bytes = Cow::Owned(word.to_string_lossy().into_owned().into_bytes());
    decode_line(&bytes).into_owned()
}

/// The path a decoded command-line word names: the bytes the word was read
/// from ([`encode_line`]), so a file name outside UTF-8 opens.
pub fn path_of(text: &str) -> Cow<'_, Path> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        match encode_line(text) {
            Cow::Borrowed(bytes) => Cow::Borrowed(Path::new(OsStr::from_bytes(bytes))),
            Cow::Owned(bytes) => Cow::Owned(std::ffi::OsString::from_vec(bytes).into()),
        }
    }
    #[cfg(not(unix))]
    {
        Cow::Borrowed(Path::new(text))
    }
}

/// Each name's first position in `names`.
fn index_of(names: &[String]) -> HashMap<String, usize> {
    let mut index = HashMap::with_capacity(names.len());
    for (i, name) in names.iter().enumerate() {
        index.entry(name.clone()).or_insert(i);
    }
    index
}

impl Headers {
    /// Parses the header lines; `info_field` is the INFO key the consequence
    /// annotations sit under (`CSQ` unless `--vcf_info_field` says otherwise).
    pub fn parse(raw: Vec<String>, info_field: &str) -> Headers {
        let mut headers = Headers::default();
        for line in &raw {
            // Every `#` of the line goes (358), so `#CHROM` reads as `CHROM`.
            let stripped = line.replace('#', "");
            if line.starts_with("##") {
                match first_group(&INFO_ID, &stripped) {
                    Some(id) if id == info_field => {
                        let fields = first_group(&FORMAT_LIST, &stripped)
                            .map(|list| list.split('|').map(str::to_string).collect())
                            .unwrap_or_default();
                        headers.consequence_fields = Some(fields);
                    }
                    Some(id) => {
                        headers.allowed.insert(id.to_string());
                    }
                    None => {
                        if let Some(name) = first_group(&DESCRIBED_FIELD, &stripped) {
                            headers.allowed.insert(name.to_string());
                        }
                    }
                }
            } else {
                // Perl's `split` without a limit, as Ensembl VEP calls it, drops the
                // trailing empty fields.
                let mut columns: Vec<&str> = stripped.split('\t').collect();
                while columns.last() == Some(&"") {
                    columns.pop();
                }
                headers.columns = columns.into_iter().map(str::to_string).collect();
            }
        }
        let fields = headers.fields().to_vec();
        headers.allowed.extend(fields.iter().cloned());
        headers.allowed.extend(headers.columns.iter().cloned());
        headers.column_index = index_of(&headers.columns);
        headers.field_index = index_of(&fields);
        headers.raw = raw;
        headers
    }

    /// The fields a filter reads: the consequence fields when the header block
    /// names them, else the columns (216).
    pub fn fields(&self) -> &[String] {
        self.consequence_fields.as_deref().unwrap_or(&self.columns)
    }

    /// The field names as `--list` prints them (221): in the order of the
    /// bytes each was read from, as Perl's `sort` orders byte strings.
    pub fn listed(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.allowed.iter().map(String::as_str).collect();
        names.sort_unstable_by(|a, b| encode_line(a).as_ref().cmp(encode_line(b).as_ref()));
        names
    }
}

/// The format of a data line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Vcf,
    Tab,
}

/// The format of `line` as `detect_format` reads it (423-439), from its
/// blank-separated words, leading blanks skipped: Perl's `split /\s+/` (425)
/// keeps a leading blank's empty first field, which fails the word-character
/// test (429) and makes such a line a tab line there. The word character of
/// that test is ASCII, as `\w` is on a byte string.
pub fn detect_format(line: &str) -> Format {
    let words: Vec<&str> = line.split(perl_space).filter(|w| !w.is_empty()).collect();
    let word = |i: usize| words.get(i).copied().unwrap_or("");
    let is_vcf = word(0)
        .bytes()
        .any(|b| b.is_ascii_alphanumeric() || b == b'_')
        && !word(1).is_empty()
        && word(1).bytes().all(|b| b.is_ascii_digit())
        && perl_true(word(3))
        && word(3).bytes().all(|b| {
            matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'N' | b'-' | b'.'
            )
        })
        && perl_true(word(4));
    if is_vcf {
        Format::Vcf
    } else {
        Format::Tab
    }
}

/// Perl's truth of a string: everything but the empty string and `0`.
fn perl_true(s: &str) -> bool {
    !s.is_empty() && s != "0"
}

/// Whether `c` is a blank to Perl's `\s` on the byte string a line is read
/// as (252, 425): the six ASCII blanks, so a non-ASCII space inside an
/// entry is text.
pub(crate) fn perl_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')
}

/// A value as the hash keeps it: undefined when empty, and when `-` under
/// every key but `Allele` (407-411).
fn clean<'a>(key: &str, value: Option<&'a str>) -> Option<&'a str> {
    match value {
        Some("") => None,
        Some("-") if key != "Allele" => None,
        other => other,
    }
}

/// The values of `line` split on `delimiter`, one per name of `names` in
/// order, undefined where the line has no field (393-395, 413-416), each
/// cleaned.
fn split_values<'a>(line: &'a str, names: &[String], delimiter: char) -> Vec<Option<&'a str>> {
    let mut pieces = line.split(delimiter);
    names
        .iter()
        .map(|name| clean(name, pieces.next()))
        .collect()
}

/// What a data line's hashes share: the line's own columns, and the
/// `key=value` pairs of its INFO column (a VCF line) or of its `Extra` column
/// (a tab line). An INFO key without `=` is a flag, defined and empty; an
/// `Extra` word without `=` is undefined; a value is everything after the
/// first `=`.
#[derive(Debug)]
struct LineBase<'a> {
    columns: &'a [String],
    column_index: &'a HashMap<String, usize>,
    column_values: Vec<Option<&'a str>>,
    /// Each pair with whether it is cleaned: an INFO flag is a defined, empty
    /// field everywhere, which the clean never undefines.
    pairs: Vec<Pair<'a>>,
    /// Whether the pairs are cleaned on the line's own hash. The INFO pairs
    /// are set after that hash is cleaned (252-258) and stand as written
    /// there; every entry's hash copies and cleans them (389, 407-411). The
    /// `Extra` pairs are set before the clean (400-411).
    pairs_cleaned: bool,
}

/// One `key=value` pair of an INFO or `Extra` column: the key, the value
/// (empty for a key without `=`) and whether the clean applies to it.
#[derive(Debug, Clone, Copy)]
struct Pair<'a> {
    key: &'a str,
    value: &'a str,
    cleanable: bool,
}

/// One hash a filter is evaluated on: an entry's, with the consequence fields
/// over the line's columns and pairs, or the line's own. A key set later in
/// `parse_line` reads over one set earlier: `INFO` is undefined everywhere
/// (398), then the consequence fields, then the pairs, then the columns.
#[derive(Debug, Clone, Copy)]
pub struct Row<'a, 'b> {
    base: &'b LineBase<'a>,
    fields: Option<(FieldSet<'a>, &'b [Option<&'a str>])>,
}

/// The consequence fields the header names, with each name's position.
#[derive(Debug, Clone, Copy)]
struct FieldSet<'a> {
    names: &'a [String],
    index: &'a HashMap<String, usize>,
}

impl Row<'_, '_> {
    /// The value of the pair under `key`, the last one when the key repeats,
    /// as the hash assignment of 257 and 403 leaves it.
    fn pair(&self, key: &str) -> Option<Option<&str>> {
        let pair = self.base.pairs.iter().rfind(|p| p.key == key)?;
        Some(
            if pair.cleanable && (self.base.pairs_cleaned || self.fields.is_some()) {
                clean(key, Some(pair.value))
            } else {
                Some(pair.value)
            },
        )
    }
}

impl FilterData for Row<'_, '_> {
    fn value(&self, key: &str) -> Option<Cow<'_, str>> {
        if key == "INFO" {
            return None;
        }
        if let Some((fields, values)) = self.fields {
            if let Some(&i) = fields.index.get(key) {
                return values[i].map(Cow::Borrowed);
            }
        }
        if let Some(value) = self.pair(key) {
            return value.map(Cow::Borrowed);
        }
        self.base
            .column_index
            .get(key)
            .and_then(|&i| self.base.column_values[i])
            .map(Cow::Borrowed)
    }

    fn has_key(&self, key: &str) -> bool {
        key == "INFO"
            || self
                .fields
                .is_some_and(|(fields, _)| fields.index.contains_key(key))
            || self.base.pairs.iter().any(|p| p.key == key)
            || self.base.column_index.contains_key(key)
    }

    fn keys(&self) -> Vec<&str> {
        let mut keys: Vec<&str> = self.base.columns.iter().map(String::as_str).collect();
        for pair in &self.base.pairs {
            if !keys.contains(&pair.key) {
                keys.push(pair.key);
            }
        }
        if let Some((fields, _)) = self.fields {
            for field in fields.names {
                if !keys.contains(&field.as_str()) {
                    keys.push(field);
                }
            }
        }
        if !keys.contains(&"INFO") {
            keys.push("INFO");
        }
        keys
    }
}

/// One data line parsed for evaluation: the hashes a filter is evaluated on,
/// one per consequence entry (the line's own hash when it has none), and the
/// consequence field's items as written, for `--only_matched` to rewrite.
#[derive(Debug)]
pub struct ParsedLine<'a> {
    base: LineBase<'a>,
    /// One value list per entry, the items' entries in line order.
    entry_values: Vec<Vec<Option<&'a str>>>,
    fields: Option<FieldSet<'a>>,
    /// The INFO items under the consequence field's name that hold an entry,
    /// in line order; empty when the line has none. An output written with
    /// `--keep_csq` carries two, the earlier annotation first.
    pub items: Vec<ConsequenceItem<'a>>,
}

/// One INFO item under the consequence field's name: the byte range of its
/// value within the line, and the entries that value holds, as written.
#[derive(Debug)]
pub struct ConsequenceItem<'a> {
    pub span: (usize, usize),
    pub entries: Vec<&'a str>,
}

impl ParsedLine<'_> {
    /// The hashes to evaluate, one per consequence entry of every item in
    /// line order, or the line's own hash alone. `filter_vep` evaluates every
    /// item's entries too (261-266): each match of its pattern appends its
    /// entries to one list and parses the whole list again, so the line
    /// passes when any entry of any item passes (282-291).
    pub fn rows(&self) -> Vec<Row<'_, '_>> {
        if self.entry_values.is_empty() {
            return vec![Row {
                base: &self.base,
                fields: None,
            }];
        }
        self.entry_values
            .iter()
            .map(|values| Row {
                base: &self.base,
                fields: self.fields.map(|fields| (fields, values.as_slice())),
            })
            .collect()
    }

    /// Whether a VCF line carried no consequence field (267-271).
    pub fn lacks_consequence_field(&self) -> bool {
        self.items.is_empty()
    }

    /// `line` with each item's value rewritten to its passing entries alone,
    /// joined by `,`; `pass` holds one verdict per row of [`Self::rows`] in
    /// order. `None` when every entry passes, the line standing as written
    /// (294-297). `filter_vep`'s rewrite pairs each passing hash with the
    /// chunk at its index in the one chunk list (290), though with a second
    /// item its hashes hold the first item's entries twice (263-265), and
    /// lands on the first item alone (296); here every item is rewritten to
    /// its own entries.
    pub fn only_matched(&self, line: &str, pass: &[bool]) -> Option<String> {
        let mut rewritten = String::new();
        let mut written_to = 0;
        let mut changed = false;
        let mut verdicts = pass.iter();
        for item in &self.items {
            let kept: Vec<&str> = item
                .entries
                .iter()
                .zip(&mut verdicts)
                .filter_map(|(&entry, &passed)| passed.then_some(entry))
                .collect();
            if kept.len() == item.entries.len() {
                continue;
            }
            let (start, end) = item.span;
            rewritten.push_str(&line[written_to..start]);
            rewritten.push_str(&kept.join(","));
            written_to = end;
            changed = true;
        }
        if !changed {
            return None;
        }
        rewritten.push_str(&line[written_to..]);
        Some(rewritten)
    }
}

/// The pairs of an `Extra` column (400-405). A word without `=` is a key
/// whose value Perl's `split` leaves undefined (402); `key=` holds the
/// empty value, which the clean undefines (407-411). Both read undefined.
fn extra_pairs(column: &str) -> Vec<Pair<'_>> {
    column
        .split(';')
        .map(|item| {
            let (key, value) = item.split_once('=').unwrap_or((item, ""));
            Pair {
                key,
                value,
                cleanable: true,
            }
        })
        .collect()
}

/// One item of an INFO column; a key without `=` is a flag, defined and
/// empty, and its value is never cleaned away.
fn info_pair(item: &str) -> Pair<'_> {
    match item.split_once('=') {
        Some((key, value)) => Pair {
            key,
            value,
            cleanable: true,
        },
        None => Pair {
            key: item,
            value: "",
            cleanable: false,
        },
    }
}

/// A tab line (241-243): one hash over the fields the header names, with its
/// `Extra` column exploded (400-405).
pub fn parse_tab_line<'a>(line: &'a str, headers: &'a Headers) -> ParsedLine<'a> {
    let column_values = split_values(line, headers.fields(), '\t');
    let extra = headers
        .field_index
        .get("Extra")
        .and_then(|&i| column_values[i]);
    ParsedLine {
        base: LineBase {
            columns: headers.fields(),
            column_index: &headers.field_index,
            column_values,
            pairs: extra.map(extra_pairs).unwrap_or_default(),
            pairs_cleaned: true,
        },
        entry_values: Vec::new(),
        fields: None,
        items: Vec::new(),
    }
}

/// A VCF line (246-271): the columns and the INFO keys in one hash, a `CSQ`
/// key always left out of it (254-256), and the entries of every item under
/// the consequence field's name each parsed over a copy; a line without the
/// field, or whose items under its name hold no entry, is evaluated on the
/// columns and INFO keys alone. The INFO column is the eighth
/// blank-separated word (252), as `filter_vep` reads it.
pub fn parse_vcf_line<'a>(line: &'a str, headers: &'a Headers, info_field: &str) -> ParsedLine<'a> {
    let column_values = split_values(line, &headers.columns, '\t');
    let mut info_pairs = Vec::new();
    let mut items = Vec::new();
    if let Some((info_start, info)) = nth_word(line, 7) {
        let mut at = info_start;
        for item in info.split(';') {
            let pair = info_pair(item);
            if pair.key == info_field {
                // Perl's `split` without a limit drops the trailing empty fields (262),
                // so an empty value, or one of commas alone, holds no entry.
                let mut entries: Vec<&str> = pair.value.split(',').collect();
                while entries.last() == Some(&"") {
                    entries.pop();
                }
                if !entries.is_empty() {
                    let span = (at + pair.key.len() + 1, at + item.len());
                    items.push(ConsequenceItem { span, entries });
                }
            }
            // The consequence field is read apart; `CSQ` is always left out (254-256).
            if pair.key != "CSQ" {
                info_pairs.push(pair);
            }
            at += item.len() + 1;
        }
    }
    let entry_values = items
        .iter()
        .flat_map(|item| &item.entries)
        .map(|entry| split_values(entry, headers.fields(), '|'))
        .collect();
    ParsedLine {
        base: LineBase {
            columns: &headers.columns,
            column_index: &headers.column_index,
            column_values,
            pairs: info_pairs,
            pairs_cleaned: false,
        },
        entry_values,
        fields: (!items.is_empty()).then_some(FieldSet {
            names: headers.fields(),
            index: &headers.field_index,
        }),
        items,
    }
}

/// The `n`th (zero-based) blank-separated word of `line` with its byte
/// offset, leading blanks skipped. Perl's `split /\s+/` (252) keeps an
/// empty first field there, so on a line starting with a blank `filter_vep`
/// reads the FILTER column's text as the INFO items.
fn nth_word(line: &str, n: usize) -> Option<(usize, &str)> {
    let mut words = 0;
    let mut start = None;
    for (i, c) in line.char_indices() {
        match (start, perl_space(c)) {
            (None, false) => start = Some(i),
            (Some(s), true) => {
                if words == n {
                    return Some((s, &line[s..i]));
                }
                words += 1;
                start = None;
            }
            _ => {}
        }
    }
    match start {
        Some(s) if words == n => Some((s, &line[s..])),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn headers(lines: &[&str]) -> Headers {
        Headers::parse(lines.iter().map(|l| l.to_string()).collect(), "CSQ")
    }

    const VCF_HEADER: &[&str] = &[
        "##fileformat=VCFv4.2",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ];

    fn value<'a>(row: &Row<'a, '_>, key: &str) -> Option<String> {
        row.value(key).map(|v| v.into_owned())
    }

    /// A line that is UTF-8 is borrowed both ways; each byte of a sequence that
    /// is not UTF-8, a lone Latin-1 byte, a truncated three-byte sequence or
    /// one cut off by the end of the line, is one character that is neither
    /// alphanumeric nor a blank and has no case, and is written back as that
    /// byte; a line that already holds an escape character comes back as
    /// written.
    #[test]
    fn a_line_round_trips_through_decode_and_encode() {
        let ascii = b"21\t100\tr1\tA\tG";
        assert!(matches!(decode_line(ascii), Cow::Borrowed(_)));
        assert!(matches!(encode_line("21\t100\tr1\tA\tG"), Cow::Borrowed(_)));
        let utf8 = "21\t100\tr\u{e9}\tA\tG";
        assert!(matches!(decode_line(utf8.as_bytes()), Cow::Borrowed(text) if text == utf8));
        let latin1 = b"21\t100\tr\xE9\tA\t\xE2\x82B\tG\xC3";
        let text = decode_line(latin1);
        assert_eq!(
            text.as_ref(),
            "21\t100\tr\u{10FFE9}\tA\t\u{10FFE2}\u{10FF82}B\tG\u{10FFC3}"
        );
        let escape_chars: Vec<char> = text
            .chars()
            .filter(|c| escaped_byte(*c).is_some())
            .collect();
        assert_eq!(escape_chars.len(), 4);
        assert!(escape_chars
            .iter()
            .all(|c| !c.is_alphanumeric() && !c.is_whitespace()));
        assert!(escape_chars
            .iter()
            .all(|c| c.to_lowercase().eq([*c]) && c.to_uppercase().eq([*c])));
        assert_eq!(encode_line(&text).as_ref(), latin1);
        let literal = "r\u{10FFE9}\u{100000}";
        let text = decode_line(literal.as_bytes());
        assert_eq!(
            text.chars().filter(|c| escaped_byte(*c).is_some()).count(),
            4
        );
        assert_eq!(encode_line(&text).as_ref(), literal.as_bytes());
        assert_eq!(encode_line(literal).as_ref(), b"r\xE9\xF4\x80\x80\x80");
        assert_eq!(decode_line(b"").as_ref(), "");
        assert_eq!(encode_line("").as_ref(), b"");
    }

    /// A command-line word outside UTF-8 decodes as a line does, and a path
    /// decoded that way opens the file the shell named.
    #[cfg(unix)]
    #[test]
    fn a_word_and_a_path_round_trip_through_their_bytes() {
        use std::os::unix::ffi::OsStrExt;
        let word = OsStr::from_bytes(b"r\xE9");
        assert_eq!(decode_word(word), "r\u{10FFE9}");
        assert_eq!(decode_word(OsStr::new("r\u{e9}")), "r\u{e9}");
        assert_eq!(path_of("r\u{10FFE9}").as_os_str().as_bytes(), b"r\xE9");
        assert!(matches!(path_of("plain.vcf"), Cow::Borrowed(_)));
        assert_eq!(path_of("plain.vcf").as_os_str(), "plain.vcf");
    }

    /// `--list` orders the names by the bytes each was read from, as Perl's
    /// `sort` orders byte strings (221): a raw 0x80 before the UTF-8 `é`
    /// (`C3 A9`), that before a raw 0xE9, that before a four-byte character.
    #[test]
    fn listed_names_are_in_byte_order() {
        let h = Headers::parse(
            [
                b"##INFO=<ID=D\xE9,Number=1,Type=String,Description=\"x\">".as_slice(),
                b"##INFO=<ID=D\xF0\x9F\x98\x80,Number=1,Type=String,Description=\"x\">",
                b"##INFO=<ID=D\xC3\xA9,Number=1,Type=String,Description=\"x\">",
                b"##INFO=<ID=D\x80,Number=1,Type=String,Description=\"x\">",
                b"#D\tDZ",
            ]
            .iter()
            .map(|line| decode_line(line).into_owned())
            .collect(),
            "CSQ",
        );
        let listed: Vec<Vec<u8>> = h
            .listed()
            .iter()
            .map(|name| encode_line(name).into_owned())
            .collect();
        assert_eq!(
            listed,
            [
                b"D".as_slice(),
                b"DZ",
                b"D\x80",
                b"D\xC3\xA9",
                b"D\xE9",
                b"D\xF0\x9F\x98\x80"
            ]
        );
    }

    /// The consequence fields come from the `Format:` list, the columns from
    /// the `#CHROM` line with its `#` removed, and `--list` prints the union of
    /// both with every other `##INFO` id, which `allowed` holds sorted by code
    /// point and `listed` prints in byte order (361-379, 221).
    #[test]
    fn vcf_headers_name_the_consequence_fields_columns_and_info_ids() {
        let h = headers(VCF_HEADER);
        assert_eq!(
            h.consequence_fields.as_deref().unwrap(),
            ["Allele", "Consequence", "IMPACT", "SYMBOL"]
        );
        assert_eq!(
            h.columns,
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        );
        assert_eq!(h.fields(), h.consequence_fields.as_deref().unwrap());
        let listed: Vec<&str> = h.allowed.iter().map(String::as_str).collect();
        assert_eq!(
            listed,
            [
                "ALT",
                "Allele",
                "CHROM",
                "Consequence",
                "DP",
                "FILTER",
                "ID",
                "IMPACT",
                "INFO",
                "POS",
                "QUAL",
                "REF",
                "SYMBOL"
            ]
        );
    }

    /// A tab output's `## <name> : <description>` lines and its column line
    /// name the fields; the fields a filter reads are the columns (371-373, 216).
    #[test]
    fn tab_headers_name_the_described_fields_and_columns() {
        let h = headers(&[
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3",
            "## Output produced at 2026-01-01 09:20:15",
            "## Column descriptions:",
            "## Uploaded_variation : Identifier of uploaded variant",
            "## Extra column keys:",
            "## IMPACT : Subjective impact classification of consequence type",
            "#Uploaded_variation\tLocation\tAllele\tExtra",
        ]);
        assert!(h.consequence_fields.is_none());
        assert_eq!(
            h.fields(),
            ["Uploaded_variation", "Location", "Allele", "Extra"]
        );
        let listed: Vec<&str> = h.allowed.iter().map(String::as_str).collect();
        assert_eq!(
            listed,
            [
                "Allele",
                "Extra",
                "IMPACT",
                "Location",
                "Uploaded_variation"
            ]
        );
    }

    /// The first word's word character, the numeric second word, the allele
    /// fourth word and a fifth that is neither empty nor `0` make a VCF line (428-433);
    /// a tab output's `Location` second word makes a tab line. The words are split
    /// on Perl's `\s` over a byte string (425), so a non-ASCII space is text,
    /// and the word character is ASCII, so a first word of one UTF-8 letter
    /// makes a tab line. A leading blank is skipped.
    #[test]
    fn format_is_detected_from_the_first_data_line() {
        assert_eq!(
            detect_format("\u{e9}\t10541165\tx\tG\tA\t50\t.\tDP=13"),
            Format::Tab
        );
        assert_eq!(
            detect_format("\u{e9}_\t10541165\tx\tG\tA\t50\t.\tDP=13"),
            Format::Vcf
        );
        assert_eq!(
            detect_format("21\t10541165\tx\tG\tA\t50\t.\tDP=13"),
            Format::Vcf
        );
        assert_eq!(
            detect_format(" 21\t10541165\tx\tG\tA\t50\t.\tDP=13"),
            Format::Vcf
        );
        assert_eq!(
            detect_format("chr21\t10541165\t.\tGA\t-\t.\t.\t."),
            Format::Vcf
        );
        assert_eq!(detect_format("21\t10541165\t.\tN\tA,C"), Format::Vcf);
        assert_eq!(
            detect_format("rs1\t21:10541165\tA\tENSG1\tENST1"),
            Format::Tab
        );
        assert_eq!(detect_format("21\t10541165\tx\tG\t0"), Format::Tab);
        assert_eq!(detect_format("21\t10541165\tx\tG"), Format::Tab);
        assert_eq!(detect_format("21\t1054116x\tx\tG\tA"), Format::Tab);
        assert_eq!(detect_format(""), Format::Tab);
        assert_eq!(
            detect_format("21\t10541165\tx\u{a0}y\tG\tA\t50\t.\tDP=13"),
            Format::Vcf
        );
        assert_eq!(detect_format("21\u{a0}10541165\tx\tG\tA"), Format::Tab);
    }

    /// A VCF line's columns and INFO keys reach every entry's hash with the
    /// entry's own fields, which read over them; an INFO flag is defined and
    /// empty, a value holding a second `=` is read whole, and `-` under a
    /// field other than `Allele` is undefined (252-271, 407-411).
    #[test]
    fn vcf_line_is_one_hash_per_consequence_entry() {
        let h = headers(VCF_HEADER);
        let line = "21\t10\tid1\tG\tA,C\t50\tPASS\tDP=13;DB;URL=a=b;SYMBOL=info;CSQ=A|missense_variant|MODERATE|TPTE,C|intron_variant|MODIFIER|-\tGT\t0/1";
        let parsed = parse_vcf_line(line, &h, "CSQ");
        assert_eq!(parsed.items.len(), 1);
        assert_eq!(
            parsed.items[0].entries,
            [
                "A|missense_variant|MODERATE|TPTE",
                "C|intron_variant|MODIFIER|-"
            ]
        );
        let (start, end) = parsed.items[0].span;
        assert_eq!(
            &line[start..end],
            "A|missense_variant|MODERATE|TPTE,C|intron_variant|MODIFIER|-"
        );
        let rows = parsed.rows();
        assert_eq!(rows.len(), 2);
        let first = &rows[0];
        assert_eq!(value(first, "CHROM").as_deref(), Some("21"));
        assert_eq!(value(first, "FILTER").as_deref(), Some("PASS"));
        assert_eq!(value(first, "DP").as_deref(), Some("13"));
        assert_eq!(value(first, "DB").as_deref(), Some(""));
        assert_eq!(value(first, "URL").as_deref(), Some("a=b"));
        assert_eq!(value(first, "SYMBOL").as_deref(), Some("TPTE"));
        assert!(first.has_key("INFO") && first.value("INFO").is_none());
        assert!(!first.has_key("CSQ"));
        let second = &rows[1];
        assert_eq!(
            value(second, "Consequence").as_deref(),
            Some("intron_variant")
        );
        assert!(second.has_key("SYMBOL") && second.value("SYMBOL").is_none());
        assert!(!parsed.lacks_consequence_field());
        assert_eq!(
            first.keys(),
            [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "DP",
                "DB",
                "URL",
                "SYMBOL",
                "Allele",
                "Consequence",
                "IMPACT"
            ]
        );
    }

    /// A line without the consequence field is evaluated on its columns and
    /// INFO keys alone, every consequence field absent from the hash (267-271),
    /// and its INFO values stand as written; an entry's hash cleans them. A
    /// key the INFO column repeats reads its last value (257).
    #[test]
    fn vcf_line_without_the_field_is_the_columns_alone() {
        let h = headers(VCF_HEADER);
        let parsed = parse_vcf_line("21\t10\tid1\tG\tA\t50\t.\tDP=-;AF=", &h, "CSQ");
        assert!(parsed.lacks_consequence_field());
        let rows = parsed.rows();
        assert_eq!(rows.len(), 1);
        assert_eq!(value(&rows[0], "DP").as_deref(), Some("-"));
        assert_eq!(value(&rows[0], "AF").as_deref(), Some(""));
        assert!(!rows[0].has_key("Consequence"));
        let repeated = parse_vcf_line("21\t10\tid1\tG\tA\t50\t.\tDP=1;AF=;DP=7", &h, "CSQ");
        let rows = repeated.rows();
        assert_eq!(value(&rows[0], "DP").as_deref(), Some("7"));
        assert_eq!(rows[0].keys().iter().filter(|k| **k == "DP").count(), 1);
        let parsed = parse_vcf_line(
            "21\t10\tid1\tG\tA\t50\t.\tDP=-;CSQ=A|stop_gained",
            &h,
            "CSQ",
        );
        assert!(parsed.rows()[0].has_key("DP") && parsed.rows()[0].value("DP").is_none());
        let empty = parse_vcf_line("21\t10\tid1\tG\tA\t50\t.\tCSQ=", &h, "CSQ");
        assert!(empty.lacks_consequence_field());
        // A trailing comma adds no entry, as Perl's `split` drops it (262).
        let trailing = parse_vcf_line("21\t10\tid1\tG\tA\t50\t.\tCSQ=A|stop_gained,,", &h, "CSQ");
        assert_eq!(trailing.items[0].entries, ["A|stop_gained"]);
        assert_eq!(trailing.rows().len(), 1);
        let commas = parse_vcf_line("21\t10\tid1\tG\tA\t50\t.\tCSQ=,;DP=1", &h, "CSQ");
        assert!(commas.lacks_consequence_field());
        // Under another INFO key the consequences are read from that key; a
        // `CSQ` beside it is still left out of the hash (254-256).
        let ann = headers(&[
            "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]);
        let ann = Headers::parse(ann.raw, "ANN");
        let parsed = parse_vcf_line(
            "21\t10\tid1\tG\tA\t50\t.\tCSQ=x;ANN=A|stop_gained",
            &ann,
            "ANN",
        );
        assert_eq!(parsed.items.len(), 1);
        assert_eq!(parsed.items[0].entries, ["A|stop_gained"]);
        let rows = parsed.rows();
        assert!(!rows[0].has_key("CSQ"));
        assert_eq!(value(&rows[0], "ANN").as_deref(), Some("A|stop_gained"));
        assert_eq!(
            value(&rows[0], "Consequence").as_deref(),
            Some("stop_gained")
        );
    }

    /// Every INFO item under the field's name with a value is read, in line
    /// order, one hash per entry across the items (261-266); an empty item is
    /// no item. `only_matched` rewrites each item to its own passing entries
    /// and leaves the line alone when every entry passes.
    #[test]
    fn vcf_line_reads_every_item_under_the_field_name() {
        let h = headers(VCF_HEADER);
        let line = "21\t10\tid1\tG\tA\t50\t.\tCSQ=A|intron_variant|MODIFIER|TPTE,A|stop_gained|HIGH|TPTE;DP=4;CSQ=;CSQ=A|missense_variant|MODERATE|TPTE";
        let parsed = parse_vcf_line(line, &h, "CSQ");
        assert!(!parsed.lacks_consequence_field());
        assert_eq!(parsed.items.len(), 2);
        assert_eq!(
            parsed.items[0].entries,
            ["A|intron_variant|MODIFIER|TPTE", "A|stop_gained|HIGH|TPTE"]
        );
        assert_eq!(
            parsed.items[1].entries,
            ["A|missense_variant|MODERATE|TPTE"]
        );
        let (start, end) = parsed.items[1].span;
        assert_eq!(&line[start..end], "A|missense_variant|MODERATE|TPTE");
        let rows = parsed.rows();
        let impacts: Vec<Option<String>> = rows.iter().map(|row| value(row, "IMPACT")).collect();
        assert_eq!(
            impacts,
            [
                Some("MODIFIER".to_string()),
                Some("HIGH".to_string()),
                Some("MODERATE".to_string())
            ]
        );
        assert_eq!(value(&rows[2], "DP").as_deref(), Some("4"));
        assert_eq!(
            parsed.only_matched(line, &[false, true, false]).as_deref(),
            Some("21\t10\tid1\tG\tA\t50\t.\tCSQ=A|stop_gained|HIGH|TPTE;DP=4;CSQ=;CSQ=")
        );
        assert_eq!(
            parsed.only_matched(line, &[true, true, false]).as_deref(),
            Some("21\t10\tid1\tG\tA\t50\t.\tCSQ=A|intron_variant|MODIFIER|TPTE,A|stop_gained|HIGH|TPTE;DP=4;CSQ=;CSQ=")
        );
        assert_eq!(
            parsed.only_matched(line, &[false, false, true]).as_deref(),
            Some("21\t10\tid1\tG\tA\t50\t.\tCSQ=;DP=4;CSQ=;CSQ=A|missense_variant|MODERATE|TPTE")
        );
        assert_eq!(parsed.only_matched(line, &[true, true, true]), None);
        let one_item = "21\t10\tid1\tG\tA\t50\t.\tDP=4;CSQ=A|stop_gained";
        let single = parse_vcf_line(one_item, &h, "CSQ");
        assert_eq!(single.only_matched(one_item, &[true]), None);
        assert_eq!(
            single.only_matched(one_item, &[false]).as_deref(),
            Some("21\t10\tid1\tG\tA\t50\t.\tDP=4;CSQ=")
        );
    }

    /// `nth_word` skips leading blanks, so the eighth word of a line starting
    /// with a blank is still its INFO column; a non-ASCII space is text and a
    /// vertical tab a blank, as Perl's `\s` on a byte string has them.
    #[test]
    fn nth_word_skips_leading_blanks() {
        assert_eq!(nth_word(" a\tb  c", 0), Some((1, "a")));
        assert_eq!(nth_word(" a\tb  c", 2), Some((6, "c")));
        assert_eq!(nth_word(" a\tb  c", 3), None);
        assert_eq!(nth_word("a\tb", 1), Some((2, "b")));
        assert_eq!(nth_word("", 0), None);
        assert_eq!(nth_word("a\u{a0}b\tc", 0), Some((0, "a\u{a0}b")));
        assert_eq!(nth_word("a\u{a0}b\tc", 1), Some((5, "c")));
        assert_eq!(nth_word("a\x0Bb", 1), Some((2, "b")));
    }

    /// A tab line's `Extra` column is exploded into its pairs, which read over
    /// the columns; an empty cell and `-` are undefined while an `Allele` of
    /// `-` stands, and a missing column is a key with no value (241-243,
    /// 400-416).
    #[test]
    fn tab_line_explodes_extra_and_undefines_dashes() {
        let h = headers(&["#Uploaded_variation\tLocation\tAllele\tGene\tConsequence\tExtra"]);
        let parsed = parse_tab_line(
            "rs1\t21:10\t-\t-\tintron_variant\tIMPACT=MODIFIER;STRAND=1;FLAGS=;Gene=ENSG1",
            &h,
        );
        let rows = parsed.rows();
        let row = &rows[0];
        assert_eq!(value(row, "Allele").as_deref(), Some("-"));
        assert_eq!(value(row, "Gene").as_deref(), Some("ENSG1"));
        assert_eq!(value(row, "IMPACT").as_deref(), Some("MODIFIER"));
        assert_eq!(value(row, "STRAND").as_deref(), Some("1"));
        assert!(row.has_key("FLAGS") && row.value("FLAGS").is_none());
        assert!(row.has_key("INFO") && row.value("INFO").is_none());
        assert_eq!(
            row.keys(),
            [
                "Uploaded_variation",
                "Location",
                "Allele",
                "Gene",
                "Consequence",
                "Extra",
                "IMPACT",
                "STRAND",
                "FLAGS",
                "INFO"
            ]
        );
        let short = parse_tab_line("rs1\t21:10", &h);
        let rows = short.rows();
        assert!(rows[0].has_key("Consequence") && rows[0].value("Consequence").is_none());
        assert!(short.items.is_empty());
    }
}
