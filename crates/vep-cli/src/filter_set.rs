// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's filter grammar, the one `--transcript_filter` reads
//! (`Bio::EnsEMBL::VEP::FilterSet`, ensembl-vep release/116.2).
//!
//! A filter string is units of `field operator value` joined by `and` and `or`,
//! grouped by parentheses and inverted by `not`. Perl citations below name
//! `FilterSet.pm` at release/116.2.
//!
//! * The operators are `eq` (`is`, `=`), `ne` (`!=`), `gt` (`>`), `lt` (`<`),
//!   `gte` (`>=`), `lte` (`<=`), `ex` (`exists`, `defined`), `nex`, `re`
//!   (`match`, `matches`, `regex`), `nre`, `in` and `is_child` (the synonym
//!   table at lines 100-117, the predicates at 704-748). A unit without an
//!   operator tests that the field is defined; `ne` without a value is `nex`
//!   (`finish_filter_node`, 370-432).
//! * `and` multiplies and `or` adds the running result, left to right with no
//!   precedence (`evaluate`, 441-515): `a or b and c` is `(a or b) and c`.
//! * A value starting with `#` names another field, whose value is compared.
//! * `in` takes a comma-separated list, or the name of a file with one value
//!   per line (385-419).
//! * `re` and `nre` are case-insensitive, unanchored regular expressions
//!   (`/$pattern/i`, 714-715), compiled here with the `regex` crate: a pattern
//!   it rejects (back-references, look-around) is an error naming the pattern.
//! * `is_child` passes when the value matches the field as `re` does, else when
//!   the field names a descendant of the value in the Sequence Ontology
//!   (722-748; the ancestry is [`vep_core::so_ontology`]).
//! * A field name resolves to a key of the data exactly, then lower-cased,
//!   upper-cased, prefixed with `_`, then as the unique case-insensitive prefix
//!   of one key; a prefix of several keys is an error (`get_input`, 533-613).
//! * A field value of the form `text(number)` or `text:number` is compared
//!   through its number when the filter value is numeric (589-610), and a value
//!   holding several entries joined by `,` or `&` after a digit passes when any
//!   entry passes (487-498).

use std::borrow::Cow;
use std::collections::HashSet;
use std::fmt;
use std::sync::LazyLock;

use regex::Regex;
use vep_core::coordinate::Strand;
use vep_core::so_ontology;
use vep_core::transcript::Transcript;

/// A filter string Ensembl VEP would refuse, with Ensembl VEP's message
/// (its own `ERROR: ` prefix removed), or a regular expression outside the
/// supported dialect.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterError(String);

impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for FilterError {}

/// The data a filter is evaluated against: a value per key, `None` for a key
/// the data does not define.
pub trait FilterData {
    fn value(&self, key: &str) -> Option<Cow<'_, str>>;
}

/// The keys of a transcript as Ensembl VEP's cache stores them, the fields a
/// `--transcript_filter` can name and the set its prefixes resolve over. The
/// first twenty are the scalar fields of `AnnotationSource/Database/Transcript.pm`
/// and of the transcript itself; `dbID`, `cdna_coding_start`, `cdna_coding_end`,
/// `translation` and `_gene_version` are read from the transcript too. `slice`,
/// `attributes`, `_gene`, `_trans_exon_array` and `_variation_effect_feature_cache`
/// hold an object in the cache and `_vep_lazy_loaded` is `1` on every transcript
/// the filter sees, so each is defined and only `exists` and `nex` read it
/// meaningfully. `description`, `coding_region_start` and `coding_region_end`
/// are undefined, as the cache stores them. `_source_cache` is absent,
/// as it is from every cache but a merged one (`Transcript.pm:278`), so `_s`
/// names `_swissprot` alone. `display_xref`, `external_db`,
/// `external_display_name`, `external_name` and `external_status` are left
/// out: the cache holds them, undefined, only on the transcripts without a
/// gene symbol, and with them the prefix `e` would be ambiguous instead of
/// naming `end`, as Ensembl VEP resolves it on a transcript with a symbol.
pub const TRANSCRIPT_KEYS: [&str; 34] = [
    "stable_id",
    "version",
    "biotype",
    "start",
    "end",
    "strand",
    "is_canonical",
    "source",
    "_gene_stable_id",
    "_gene_symbol",
    "_gene_symbol_source",
    "_gene_hgnc_id",
    "_ccds",
    "_refseq",
    "_protein",
    "_swissprot",
    "_trembl",
    "_uniparc",
    "_uniprot_isoform",
    "_gene_phenotype",
    "dbID",
    "cdna_coding_start",
    "cdna_coding_end",
    "translation",
    "description",
    "coding_region_start",
    "coding_region_end",
    "slice",
    "attributes",
    "_gene",
    "_gene_version",
    "_trans_exon_array",
    "_variation_effect_feature_cache",
    "_vep_lazy_loaded",
];

impl FilterData for Transcript {
    /// The cache stores `is_canonical` only on the canonical transcript, so
    /// the key is undefined on every other one.
    fn value(&self, key: &str) -> Option<Cow<'_, str>> {
        match key {
            "stable_id" => Some(Cow::Borrowed(&*self.stable_id)),
            "version" => self.version.map(|v| Cow::Owned(v.to_string())),
            "biotype" => Some(Cow::Borrowed(&*self.biotype)),
            "start" => Some(Cow::Owned(self.start.to_string())),
            "end" => Some(Cow::Owned(self.end.to_string())),
            "strand" => Some(Cow::Borrowed(match self.strand {
                Strand::Forward => "1",
                Strand::Reverse => "-1",
            })),
            "is_canonical" => self.canonical.then_some(Cow::Borrowed("1")),
            "source" => Some(Cow::Borrowed(&self.source)),
            "_gene_stable_id" => Some(Cow::Borrowed(&*self.gene_stable_id)),
            "_gene_version" => self.gene_version.map(|v| Cow::Owned(v.to_string())),
            "_gene_symbol" => self.gene_symbol.as_deref().map(Cow::Borrowed),
            "_gene_symbol_source" => self.gene_symbol_source.as_deref().map(Cow::Borrowed),
            "_gene_hgnc_id" => self.hgnc_id.as_deref().map(Cow::Borrowed),
            "_ccds" => self.ccds.as_deref().map(Cow::Borrowed),
            "_refseq" => self.refseq.as_deref().map(Cow::Borrowed),
            "_protein" => self.protein_id.as_deref().map(Cow::Borrowed),
            "_swissprot" => self.swissprot.as_deref().map(Cow::Borrowed),
            "_trembl" => self.trembl.as_deref().map(Cow::Borrowed),
            "_uniparc" => self.uniparc.as_deref().map(Cow::Borrowed),
            "_uniprot_isoform" => self.uniprot_isoform.as_deref().map(Cow::Borrowed),
            "_gene_phenotype" => self.gene_phenotype.map(|n| Cow::Owned(n.to_string())),
            "dbID" => self.db_id.map(|n| Cow::Owned(n.to_string())),
            "cdna_coding_start" => self.cdna_coding_start.map(|n| Cow::Owned(n.to_string())),
            "cdna_coding_end" => self.cdna_coding_end.map(|n| Cow::Owned(n.to_string())),
            // Perl stringifies a reference as its class and address and numifies
            // it as the address, so a value of that shape keeps `is`, `ne`, a
            // regex and a numeric comparison reading as they do there
            // (`perl_number` reads the address back); the class name alone
            // would equal the written class and compare as 0.
            "translation" => self
                .translation
                .as_ref()
                .map(|t| reference("Bio::EnsEMBL::Translation=HASH", t)),
            "slice" => Some(reference("Bio::EnsEMBL::Slice=HASH", &self.chr)),
            "_gene" => Some(reference("Bio::EnsEMBL::Gene=HASH", &self.gene_stable_id)),
            "attributes" => Some(reference("ARRAY", &self.attributes)),
            "_trans_exon_array" => Some(reference("ARRAY", &self.exons)),
            "_variation_effect_feature_cache" => self.vefc.as_ref().map(|v| reference("HASH", v)),
            "_vep_lazy_loaded" => Some(Cow::Borrowed("1")),
            _ => None,
        }
    }
}

/// Perl's stringification of a reference, `Class=HASH(0xaddr)` or
/// `ARRAY(0xaddr)`, over the address of the Rust value standing in for the
/// object, so two objects of one transcript never compare equal.
fn reference<T>(class: &str, value: &T) -> Cow<'static, str> {
    Cow::Owned(format!("{class}({value:p})"))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Logic {
    And,
    Or,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Operator {
    Eq,
    Ne,
    Gt,
    Lt,
    Gte,
    Lte,
    Ex,
    Nex,
    Re,
    Nre,
    In,
    IsChild,
}

impl Operator {
    /// The operator a word names, directly or through the synonym table
    /// (`%FILTER_SYNONYMS`, 100-117).
    fn parse(word: &str) -> Option<Operator> {
        Some(match word {
            "eq" | "is" | "=" => Operator::Eq,
            "ne" | "!=" => Operator::Ne,
            "gt" | ">" => Operator::Gt,
            "lt" | "<" => Operator::Lt,
            "gte" | ">=" => Operator::Gte,
            "lte" | "<=" => Operator::Lte,
            "ex" | "exists" | "defined" => Operator::Ex,
            "nex" => Operator::Nex,
            "re" | "match" | "matches" | "regex" => Operator::Re,
            "nre" => Operator::Nre,
            "in" => Operator::In,
            "is_child" => Operator::IsChild,
            _ => return None,
        })
    }
}

/// The right-hand side of a unit, bound to the data's keys.
#[derive(Debug, Clone)]
enum Value {
    /// A literal word, with its compiled form under a regex operator.
    Text(String, Option<Regex>),
    /// The `in` operator's set, with its first entry as written.
    List {
        first: Option<String>,
        set: HashSet<String>,
    },
    /// `#name`: the key `name` resolved to, `None` when the data has no such key.
    Field(Option<String>),
}

/// One node of the parsed filter: a unit, or a parenthesised group when
/// `children` is non-empty. `logic` joins the node to the result of the
/// siblings before it.
#[derive(Debug, Clone)]
struct Node {
    logic: Logic,
    not: bool,
    /// The field name as written.
    field: Option<String>,
    /// The key the field name resolved to; `None` leaves the value undefined.
    key: Option<String>,
    operator: Option<Operator>,
    /// The value as written.
    value: Option<String>,
    bound: Option<Value>,
    children: Vec<Node>,
}

impl Node {
    fn new(logic: Logic) -> Node {
        Node {
            logic,
            not: false,
            field: None,
            key: None,
            operator: None,
            value: None,
            bound: None,
            children: Vec::new(),
        }
    }
}

/// A compiled filter, evaluated against one [`FilterData`] at a time.
#[derive(Debug, Clone)]
pub struct FilterSet {
    root: Node,
}

const INCOMPLETE_PARENTHESES: &str = "Error parsing filter string - incomplete parentheses sets?";

impl FilterSet {
    /// Parses `filters`, every string a unit joined to the others by `and`
    /// (`parse_filters`, 209-310), and resolves each field name against
    /// `keys`, the data's key set. No filter string at all is a filter
    /// everything passes.
    pub fn new<S: AsRef<str>>(filters: &[S], keys: &[&str]) -> Result<FilterSet, FilterError> {
        let mut root = Node::new(Logic::And);
        for filter in filters {
            parse_into(&mut root, filter.as_ref(), keys)?;
        }
        Ok(FilterSet { root })
    }

    /// A `--transcript_filter`: the filter strings bound to [`TRANSCRIPT_KEYS`].
    pub fn transcript_filter<S: AsRef<str>>(filters: &[S]) -> Result<FilterSet, FilterError> {
        FilterSet::new(filters, &TRANSCRIPT_KEYS)
    }

    /// Whether `data` passes the filter.
    pub fn evaluate(&self, data: &dyn FilterData) -> bool {
        if self.root.children.is_empty() {
            return true;
        }
        evaluate(&self.root, data).unwrap_or(false)
    }
}

/// Parses one filter string as a unit under `root` (`parse_filters`,
/// 217-307): the words [`next_token`] yields fill a unit's `not`, field,
/// operator and value in turn, `and` and `or` start the next unit, and the
/// parentheses open and close groups.
fn parse_into(root: &mut Node, filter: &str, keys: &[&str]) -> Result<(), FilterError> {
    let chars: Vec<char> = filter.chars().collect();
    // The open groups, innermost last; the root is below them all.
    let mut stack: Vec<Node> = Vec::new();
    let mut current = Node::new(Logic::And);
    let mut pos = 0;
    while let Some((start, end, separator)) = next_token(&chars, pos) {
        let word: String = chars[start..end].iter().collect();
        pos = if separator.is_some() { end + 1 } else { end };

        // The empty match at the end of the string closes the unit under the
        // root, which every opened parenthesis must have returned to (226-236).
        if word.is_empty() && separator.is_none() {
            if !stack.is_empty() {
                return Err(FilterError(INCOMPLETE_PARENTHESES.to_string()));
            }
            root.children.push(finish(current, keys)?);
            return Ok(());
        }

        match word.as_str() {
            "" => {}
            "and" | "or" => {
                let logic = if word == "and" { Logic::And } else { Logic::Or };
                let finished = finish(current, keys)?;
                stack.last_mut().unwrap_or(root).children.push(finished);
                current = Node::new(logic);
            }
            "not" => current.not = true,
            _ if current.field.is_none() => current.field = Some(word),
            _ if current.operator.is_none() => {
                current.operator = Some(
                    Operator::parse(&word)
                        .ok_or_else(|| FilterError(format!("No such operator \"{word}\"")))?,
                );
            }
            _ if current.value.is_none() => current.value = Some(word),
            _ => {}
        }

        match separator {
            Some('(') => {
                stack.push(current);
                current = Node::new(Logic::And);
            }
            Some(')') => {
                let finished = finish(current, keys)?;
                stack.last_mut().unwrap_or(root).children.push(finished);
                // The group just closed is the unit the next word continues.
                current = stack
                    .pop()
                    .ok_or_else(|| FilterError(INCOMPLETE_PARENTHESES.to_string()))?;
            }
            _ => {}
        }
    }
    Ok(())
}

/// Whether `c` ends a word.
fn is_separator(c: char) -> bool {
    c.is_whitespace() || c == '(' || c == ')'
}

/// The token Ensembl VEP's tokenizer, `m/([^\(^\)^\s]*?)(\s|\(|\)|$)/g`, finds
/// from `start`: the word's bounds and the separator after it, which is `None`
/// at the end of the string; `None` as a whole once `start` is past the end.
/// A word runs to a space, a parenthesis or the end, and its character class
/// excludes `^`, so a run of characters ending in `^` is skipped and never
/// becomes a word.
fn next_token(chars: &[char], mut start: usize) -> Option<(usize, usize, Option<char>)> {
    while start <= chars.len() {
        let mut end = start;
        while end < chars.len() && !is_separator(chars[end]) && chars[end] != '^' {
            end += 1;
        }
        if end == chars.len() {
            return Some((start, end, None));
        }
        if is_separator(chars[end]) {
            return Some((start, end, Some(chars[end])));
        }
        start += 1;
    }
    None
}

/// Completes a unit (`finish_filter_node`, 370-424): the default operator,
/// the `in` list, the regular expression, the `is_child` term, and the field
/// name resolved against `keys`. A group's own field and value are never read,
/// so only a unit resolves and compiles them; a group's `in` list is still read
/// as Ensembl VEP reads it.
fn finish(mut node: Node, keys: &[&str]) -> Result<Node, FilterError> {
    let mut operator = node.operator.unwrap_or(Operator::Ex);
    if operator == Operator::Ne && node.value.is_none() {
        operator = Operator::Nex;
    }
    node.operator = Some(operator);
    if operator == Operator::In {
        let (first, set) = read_list(node.value.as_deref())?;
        node.bound = Some(Value::List { first, set });
    }
    if !node.children.is_empty() {
        return Ok(node);
    }
    node.key = resolve(node.field.as_deref().unwrap_or(""), keys)?;
    node.bound = match (&node.value, operator) {
        (_, Operator::In) => node.bound,
        (Some(v), _) if v.starts_with('#') => Some(Value::Field(resolve(&v[1..], keys)?)),
        // An `is_child` value outside the Sequence Ontology is admitted: `filter_is_child`
        // tests the regex first (FilterSet.pm:726) and asks the ontology only on a miss.
        (Some(v), Operator::Re | Operator::Nre | Operator::IsChild) => {
            Some(Value::Text(v.clone(), Some(compile(v)?)))
        }
        (Some(v), _) => Some(Value::Text(v.clone(), None)),
        (None, _) => None,
    };
    Ok(node)
}

/// A case-insensitive, unanchored regular expression (`/$pattern/i`).
fn compile(pattern: &str) -> Result<Regex, FilterError> {
    Regex::new(&format!("(?i){pattern}")).map_err(|e| {
        let message = match &e {
            // The syntax report spans several lines; its last names the fault.
            regex::Error::Syntax(report) => report
                .lines()
                .last()
                .unwrap_or_default()
                .trim_start_matches("error: ")
                .to_string(),
            other => other.to_string(),
        };
        FilterError(format!("pattern \"{pattern}\" is not supported: {message}"))
    })
}

/// The `in` operator's set (385-419): a comma-separated list when the value
/// has a comma with text on both sides, else the lines of the file it names.
fn read_list(value: Option<&str>) -> Result<(Option<String>, HashSet<String>), FilterError> {
    let list = match value {
        Some(v) if !v.is_empty() && v != "0" => v,
        _ => return Err(FilterError("No list/file given".to_string())),
    };
    let bytes = list.as_bytes();
    let inline = bytes
        .iter()
        .enumerate()
        .any(|(i, b)| *b == b',' && i >= 1 && i + 1 < bytes.len());
    let entries: Vec<String> = if inline {
        // Perl's `split` drops the empty fields a trailing comma leaves.
        let mut parts: Vec<&str> = list.split(',').collect();
        while parts.last() == Some(&"") {
            parts.pop();
        }
        parts.into_iter().map(str::to_string).collect()
    } else {
        if !std::path::Path::new(list).exists() {
            return Err(FilterError(format!("Could not find/parse list {list}")));
        }
        let text = std::fs::read_to_string(list)
            .map_err(|_| FilterError(format!("Could not read from file {list}")))?;
        // Perl reads the file line by line, so a final newline ends the last
        // entry and opens no empty one; `split_inclusive` keeps that count.
        text.split_inclusive('\n')
            .map(|line| {
                line.chars()
                    .filter(|c| {
                        !matches!(
                            c,
                            '\n' | '\r' | '\x0B' | '\x0C' | '\u{85}' | '\u{2028}' | '\u{2029}'
                        )
                    })
                    .collect()
            })
            .collect()
    };
    let first = entries.first().cloned();
    Ok((first, entries.into_iter().collect()))
}

/// The key `field` names in `keys` (`get_input`, 538-586): the exact name,
/// then lower-cased, upper-cased, prefixed with `_`, then the one key it is a
/// case-insensitive prefix of; several such keys are an error, none leaves the
/// field undefined.
fn resolve(field: &str, keys: &[&str]) -> Result<Option<String>, FilterError> {
    let exact = |candidate: &str| {
        keys.iter()
            .find(|k| **k == candidate)
            .map(|k| k.to_string())
    };
    if let Some(key) = exact(field)
        .or_else(|| exact(&field.to_lowercase()))
        .or_else(|| exact(&field.to_uppercase()))
        .or_else(|| exact(&format!("_{field}")))
    {
        return Ok(Some(key));
    }
    let lower = field.to_lowercase();
    let matches: Vec<&str> = keys
        .iter()
        .copied()
        .filter(|k| k.to_lowercase().starts_with(&lower))
        .collect();
    match matches.as_slice() {
        [one] => Ok(Some(one.to_string())),
        [] => Ok(None),
        several => Err(FilterError(format!(
            "Multiple data fields found matching query field \"{field}\": {}",
            several.join(", ")
        ))),
    }
}

/// The verdict of `node` (`evaluate`, 441-512). `None` is the early `return 0`
/// of a `#field` value the data leaves undefined (477), which the node's `not`
/// does not invert.
fn evaluate(node: &Node, data: &dyn FilterData) -> Option<bool> {
    let mut result = if node.children.is_empty() {
        unit(node, data)?
    } else {
        let mut acc = true;
        for child in &node.children {
            let verdict = evaluate(child, data).unwrap_or(false);
            acc = match child.logic {
                Logic::And => acc && verdict,
                Logic::Or => acc || verdict,
            };
        }
        acc
    };
    if node.not {
        result = !result;
    }
    Some(result)
}

static EXTRACT: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"([\w\.\-]+)?:?\(?([\-\d\.e]*)\)?").unwrap());
static NUMERIC_VALUE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^[\-\d\.e]+$").unwrap());
static NUMERIC_TEXT: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^\-?\d+\.?\d*(e\-?\d+)?$").unwrap());
static MULTI_VALUE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"[0-9](,|&)").unwrap());

/// The verdict of one unit (464-504).
fn unit(node: &Node, data: &dyn FilterData) -> Option<bool> {
    let operator = node.operator.unwrap_or(Operator::Ex);
    let field = node.field.as_deref().unwrap_or("");
    let raw = node.key.as_deref().and_then(|key| data.value(key));
    // The filter value as `get_input` reads it: a word, a list's first entry,
    // or the `#name` of a field, which stays the written word there (472-476).
    let written = match &node.bound {
        Some(Value::Text(text, _)) => Some(text.as_str()),
        Some(Value::List { first, .. }) => first.as_deref(),
        Some(Value::Field(_)) | None => node.value.as_deref(),
    };
    let input = raw.map(|v| extract(v, written, field));

    let compared: Option<Cow<'_, str>> = match &node.bound {
        Some(Value::Field(key)) => Some(data.value(key.as_deref()?)?),
        Some(Value::Text(text, _)) => Some(Cow::Borrowed(text.as_str())),
        Some(Value::List { .. }) | None => None,
    };
    let value = compared.as_deref();
    let regex = match &node.bound {
        Some(Value::Text(_, regex)) => regex.as_ref(),
        _ => None,
    };
    let list = match &node.bound {
        Some(Value::List { set, .. }) => Some(set),
        _ => None,
    };

    // A unit with a value fails on a field the data does not define (481-483).
    if (value.is_some() || list.is_some()) && input.is_none() {
        return Some(false);
    }
    let input = input.as_deref();
    let predicate = |piece: Option<&str>| match operator {
        Operator::Eq => piece.unwrap_or("") == value.unwrap_or(""),
        Operator::Ne => piece.unwrap_or("") != value.unwrap_or(""),
        Operator::Gt => perl_number(piece) > perl_number(value),
        Operator::Lt => perl_number(piece) < perl_number(value),
        Operator::Gte => perl_number(piece) >= perl_number(value),
        Operator::Lte => perl_number(piece) <= perl_number(value),
        Operator::Ex => piece.is_some(),
        Operator::Nex => piece.is_none(),
        Operator::Re => matches(regex, value, piece),
        Operator::Nre => !matches(regex, value, piece),
        Operator::In => list.is_some_and(|set| set.contains(piece.unwrap_or(""))),
        Operator::IsChild => {
            matches(regex, value, piece)
                || value
                    .and_then(so_ontology::descendants)
                    .is_some_and(|names| {
                        names
                            .iter()
                            .any(|name| name.eq_ignore_ascii_case(piece.unwrap_or("")))
                    })
        }
    };
    // Several values joined by `,` or `&` after a digit pass when one passes (487-498).
    if let Some(caps) = input.and_then(|i| MULTI_VALUE.captures(i)) {
        let delimiter = caps.get(1).map_or(",", |m| m.as_str());
        let mut pieces: Vec<&str> = input.unwrap_or("").split(delimiter).collect();
        while pieces.last() == Some(&"") {
            pieces.pop();
        }
        return Some(pieces.into_iter().any(|piece| predicate(Some(piece))));
    }
    Some(predicate(input))
}

/// Whether `piece` matches the unit's regular expression: the compiled one, or
/// a `#field` value compiled here, which matches nothing when the crate
/// rejects it.
fn matches(regex: Option<&Regex>, value: Option<&str>, piece: Option<&str>) -> bool {
    let piece = piece.unwrap_or("");
    match regex {
        Some(re) => re.is_match(piece),
        None => value
            .and_then(|v| compile(v).ok())
            .is_some_and(|re| re.is_match(piece)),
    }
}

/// `get_input`'s numeric extraction (589-610): a value shaped `text(number)` or
/// `text:number` compared against a numeric filter value is its number (its
/// text when that is itself a number); any other pairing keeps the value whole.
/// `CELL_TYPE` and the HGVS fields are exempt, as Perl's condition names them.
fn extract<'a>(input: Cow<'a, str>, value: Option<&str>, field: &str) -> Cow<'a, str> {
    if field == "CELL_TYPE" || field.contains("HGVS") {
        return input;
    }
    let Some(caps) = EXTRACT.captures(&input) else {
        return input;
    };
    let number = caps.get(2).map_or("", |m| m.as_str());
    if number.is_empty() {
        return input;
    }
    let text = caps.get(1).map_or("", |m| m.as_str());
    match value {
        // Perl's `$value &&` reads "" and "0" as false.
        Some(v) if !v.is_empty() && v != "0" && NUMERIC_VALUE.is_match(v) => {
            if NUMERIC_TEXT.is_match(text) {
                Cow::Owned(text.to_string())
            } else {
                Cow::Owned(number.to_string())
            }
        }
        _ => input,
    }
}

/// A string as Perl's numeric comparison reads it: a stringified reference
/// (`Class=HASH(0xaddr)`, `ARRAY(0xaddr)`) is its address, as Perl numifies
/// the reference itself; anything else is its leading number, else 0.
fn perl_number(s: Option<&str>) -> f64 {
    static LEADING_NUMBER: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)").unwrap());
    static REFERENCE: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"^(?:[\w:]+=)?(?:HASH|ARRAY)\(0x([0-9a-f]+)\)$").unwrap());
    let Some(s) = s else {
        return 0.0;
    };
    if let Some(address) = REFERENCE.captures(s).and_then(|caps| caps.get(1)) {
        return u64::from_str_radix(address.as_str(), 16).map_or(0.0, |n| n as f64);
    }
    LEADING_NUMBER
        .captures(s)
        .and_then(|caps| caps.get(1))
        .and_then(|m| m.as_str().parse::<f64>().ok())
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// A hash of string values, the shape Ensembl VEP's filters read.
    struct Hash(HashMap<&'static str, &'static str>);

    impl FilterData for Hash {
        fn value(&self, key: &str) -> Option<Cow<'_, str>> {
            self.0.get(key).map(|v| Cow::Borrowed(*v))
        }
    }

    fn hash(entries: &[(&'static str, &'static str)]) -> Hash {
        Hash(entries.iter().copied().collect())
    }

    fn keys(entries: &[(&'static str, &'static str)]) -> Vec<&'static str> {
        entries.iter().map(|(k, _)| *k).collect()
    }

    fn passes(filter: &str, entries: &[(&'static str, &'static str)]) -> bool {
        FilterSet::new(&[filter], &keys(entries))
            .unwrap_or_else(|e| panic!("{filter}: {e}"))
            .evaluate(&hash(entries))
    }

    /// The two examples of `FilterSet.pm`'s synopsis (lines 40-43) and the
    /// grouping example of `parse_filters` (171-201): `field1 = foo and
    /// field2 < 10` passes `{foo, 6}` and fails `{bar, 12}`; the group
    /// `(a < 10 or b > 100)` passes when either side does.
    #[test]
    fn parses_the_perl_examples() {
        assert!(passes(
            "field1 = foo and field2 < 10",
            &[("field1", "foo"), ("field2", "6")]
        ));
        assert!(!passes(
            "field1 = foo and field2 < 10",
            &[("field1", "bar"), ("field2", "12")]
        ));
        let data = [("foo", "bar"), ("a", "50"), ("b", "500")];
        assert!(passes("foo is bar and (a < 10 or b > 100)", &data));
        assert!(!passes(
            "foo is bar and (a < 10 or b > 100)",
            &[("foo", "bar"), ("a", "50"), ("b", "5")]
        ));
        assert!(!passes("foo is baz and (a < 10 or b > 100)", &data));
        // No precedence: `and` multiplies and `or` adds, left to right.
        let flags = [("a", "1"), ("b", "0"), ("c", "1")];
        assert!(passes("a is 1 or b is 1 and c is 1", &flags));
        assert!(!passes(
            "a is 1 or b is 1 and c is 1",
            &[("a", "1"), ("b", "0"), ("c", "0")]
        ));
        assert!(passes("b is 1 or a is 1 and (c is 1 or b is 1)", &flags));
        // `not` inverts a unit and a group.
        assert!(passes("not foo is baz", &data));
        assert!(!passes("not (a < 10 or b > 100)", &data));
        // A unit without an operator tests that the field is defined; `ne`
        // without a value tests that it is not.
        assert!(passes("foo", &data));
        assert!(!passes("missing", &data));
        assert!(passes("missing ne", &data));
        assert!(passes("foo exists and missing nex", &data));
        // The synonyms of the operator table.
        assert!(passes(
            "a >= 50 and a <= 50 and a != 51 and b gt 499 and b lte 500",
            &data
        ));
        assert!(passes(
            "foo matches BA and foo regex AR$ and foo nre z",
            &data
        ));
        // Another field's value through `#`.
        assert!(passes("a < #b", &data));
        assert!(!passes("a < #missing", &data));
        assert!(!passes("not a < #missing", &data));
    }

    /// The error texts of `parse_filters` and `finish_filter_node`.
    #[test]
    fn perl_parse_errors_keep_their_text() {
        let error = |filter: &str| {
            FilterSet::new(&[filter], &["a", "b"])
                .unwrap_err()
                .to_string()
        };
        assert_eq!(error("a is b and (a is c"), INCOMPLETE_PARENTHESES);
        assert_eq!(error("a is b)"), INCOMPLETE_PARENTHESES);
        assert_eq!(error("a equals b"), "No such operator \"equals\"");
        assert_eq!(error("a in"), "No list/file given");
        assert_eq!(
            error("a in protein_coding"),
            "Could not find/parse list protein_coding"
        );
    }

    /// Ensembl VEP's tokenizer excludes `^` from a word, so `^X._` reads as
    /// `X._` and `a^b` as `b`; the regex stays unanchored and case-insensitive.
    #[test]
    fn caret_never_reaches_a_pattern() {
        assert!(passes("stable_id match ^X._", &[("stable_id", "NXM_1")]));
        assert!(passes("stable_id match ^X._", &[("stable_id", "xr_001")]));
        assert!(!passes("stable_id match ^X._", &[("stable_id", "ENST1")]));
        assert!(passes("stable_id match a^b", &[("stable_id", "B")]));
    }

    /// `in` reads a comma-separated list, else a file with one value per line
    /// (`finish_filter_node`, 385-419); the file's values are exact, and its
    /// final newline adds no empty value.
    #[test]
    fn in_reads_a_one_value_per_line_file() {
        let data = [("biotype", "lncRNA"), ("stable_id", "ENST3")];
        assert!(passes("biotype in protein_coding,lncRNA", &data));
        assert!(!passes("biotype in protein_coding,miRNA,", &data));
        let dir = tempfile::tempdir().unwrap();
        let list = dir.path().join("ids.txt");
        std::fs::write(&list, "ENST1\r\nENST2\nENST3\n").unwrap();
        let filter = format!("stable_id in {}", list.display());
        assert!(passes(&filter, &data));
        assert!(!passes(
            &filter,
            &[("biotype", "lncRNA"), ("stable_id", "ENST4")]
        ));
        assert!(!passes(&filter, &[("stable_id", "")]));
        assert!(passes(
            &format!("not stable_id in {}", list.display()),
            &[("stable_id", "ENST1 ")]
        ));
        let missing = dir.path().join("none.txt");
        assert_eq!(
            FilterSet::new(
                &[format!("stable_id in {}", missing.display())],
                &["stable_id"]
            )
            .unwrap_err()
            .to_string(),
            format!("Could not find/parse list {}", missing.display())
        );
    }

    /// A field name resolves exactly, lower-cased, upper-cased, with a `_`
    /// prefix, or as the unique prefix of one key; a prefix of several keys is
    /// the error of `get_input` (574-580), and a prefix of none is undefined.
    #[test]
    fn ambiguous_prefix_is_an_error() {
        let data = [
            ("stable_id", "ENST1"),
            ("start", "10"),
            ("strand", "-1"),
            ("_gene_symbol", "ABC"),
            ("SYMBOL", "ABC"),
        ];
        assert!(passes("gene_symbol is ABC", &data));
        assert!(passes("symbol is ABC", &data));
        assert!(passes("stab is ENST1", &data));
        assert!(passes("Stable_ID is ENST1", &data));
        assert!(!passes("nothing", &data));
        let error = FilterSet::new(&["st is 10"], &keys(&data)).unwrap_err();
        assert_eq!(
            error.to_string(),
            "Multiple data fields found matching query field \"st\": stable_id, start, strand"
        );
        // `_s` names `_swissprot` alone, as it does on every cache but a merged one.
        let set = FilterSet::transcript_filter(&["_s exists"]).unwrap();
        assert_eq!(set.root.children[0].key.as_deref(), Some("_swissprot"));
        assert_eq!(
            FilterSet::transcript_filter(&["_u exists"])
                .unwrap_err()
                .to_string(),
            "Multiple data fields found matching query field \"_u\": _uniparc, _uniprot_isoform"
        );
        // `_t` reaches the exon array as well as TrEMBL; `gene` is the `_gene` object.
        assert_eq!(
            FilterSet::transcript_filter(&["_t exists"])
                .unwrap_err()
                .to_string(),
            "Multiple data fields found matching query field \"_t\": _trembl, _trans_exon_array"
        );
        let set = FilterSet::transcript_filter(&["gene exists"]).unwrap();
        assert_eq!(set.root.children[0].key.as_deref(), Some("_gene"));
        assert!(FilterSet::transcript_filter(&[""])
            .unwrap_err()
            .to_string()
            .starts_with(
                "Multiple data fields found matching query field \"\": stable_id, version"
            ));
    }

    /// `is_child` passes on a regex match of the term (`filter_is_child`, 726),
    /// else when the field's term descends from it in the Sequence Ontology
    /// (728-746), which `vep` reaches with no ontology adaptor and dies on
    /// (`No ontology adaptor set`, 736): the walk over the built-in table is
    /// the documented extension. A list of terms splits only after a digit
    /// (487), so two terms compare whole.
    #[test]
    fn is_child_walks_the_so_parent_table() {
        let missense = [("Consequence", "missense_variant")];
        assert!(passes(
            "Consequence is_child coding_sequence_variant",
            &missense
        ));
        assert!(passes("Consequence is_child missense_variant", &missense));
        assert!(passes(
            "Consequence is_child Nonsynonymous_Variant",
            &missense
        ));
        assert!(!passes("Consequence is_child intron_variant", &missense));
        assert!(passes(
            "Consequence is_child splice_region_variant",
            &[("Consequence", "splice_donor_region_variant")]
        ));
        assert!(!passes(
            "Consequence is_child coding_sequence_variant",
            &[("Consequence", "intron_variant,stop_gained")]
        ));
    }

    /// An `is_child` value outside the Sequence Ontology is no parse error:
    /// `filter_is_child` returns on the regex match before it asks the ontology
    /// (726), so `stable_id is_child ENST` passes every Ensembl transcript and
    /// `biotype is_child lncRNA` the lncRNA ones.
    #[test]
    fn is_child_matches_a_non_so_value_as_a_regex() {
        assert!(passes(
            "stable_id is_child ENST",
            &[("stable_id", "ENST00000000001")]
        ));
        assert!(!passes(
            "stable_id is_child ENST",
            &[("stable_id", "NM_000001.1")]
        ));
        assert!(passes("biotype is_child lncRNA", &[("biotype", "lncRNA")]));
        assert!(!passes(
            "biotype is_child lncRNA",
            &[("biotype", "protein_coding")]
        ));
    }

    /// A pattern outside the `regex` crate's dialect is an error naming the
    /// pattern and the construct: Perl's `\Z` and `\K` escapes, a back-reference,
    /// an unclosed class. Parentheses never reach a pattern, since the tokenizer
    /// reads them as grouping, so a pattern with a group opens one instead.
    #[test]
    fn unsupported_regex_syntax_exits_naming_the_pattern() {
        let error = |filter: &str| {
            FilterSet::transcript_filter(&[filter])
                .unwrap_err()
                .to_string()
        };
        assert_eq!(
            error(r"stable_id match ENST\Z"),
            r#"pattern "ENST\Z" is not supported: unrecognized escape sequence"#
        );
        assert_eq!(
            error(r"stable_id nre EN\KST"),
            r#"pattern "EN\KST" is not supported: unrecognized escape sequence"#
        );
        assert_eq!(
            error(r"stable_id match \1"),
            r#"pattern "\1" is not supported: backreferences are not supported"#
        );
        assert_eq!(
            error("stable_id match ENST["),
            r#"pattern "ENST[" is not supported: unclosed character class"#
        );
        assert!(FilterSet::transcript_filter(&["stable_id match (ENST)"]).is_ok());
    }

    /// `get_input`'s numeric extraction (589-610): `HGNC:5` compared with a
    /// number is 5, compared with text is itself (`hgnc_id` names no key: the
    /// `_` rule prefixes one underscore, so the key is `gene_hgnc_id`); a `tolerated(0.05)` score
    /// compares through its number and is never equal to its text; a value
    /// joined by `,` after a digit passes when one entry does (487-498).
    #[test]
    fn mixed_values_compare_through_their_number() {
        let data = [
            ("_gene_hgnc_id", "HGNC:5"),
            ("SIFT", "tolerated(0.05)"),
            ("_refseq", "NM_1.2,NM_3.4"),
        ];
        assert!(passes("gene_hgnc_id > 4 and gene_hgnc_id < 6", &data));
        assert!(!passes("hgnc_id > 4", &data));
        assert!(passes("gene_hgnc_id is HGNC:5", &data));
        assert!(passes("gene_hgnc_id is 5", &data));
        assert!(!passes("gene_hgnc_id is 6", &data));
        assert!(passes("SIFT < 0.1", &data));
        assert!(!passes("SIFT is tolerated", &data));
        assert!(passes("SIFT match tolerated", &data));
        assert!(passes("refseq is NM_3.4", &data));
        assert!(!passes("refseq is NM_3", &data));
        assert!(passes("refseq match NM_3", &data));
        // A numeric comparison reads a non-numeric string as 0.
        assert!(passes("SIFT > -1", &[("SIFT", "-")]));
        assert!(!passes("SIFT > 0", &[("SIFT", "-")]));
    }

    /// A canonical coding transcript of a gene with a symbol, the UniProt
    /// fields as the cache stores them: `-` for TrEMBL, UniParc and the
    /// isoform undefined.
    fn canonical_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST00000380152".into(),
            version: Some(8),
            db_id: Some(4_201_747),
            gene_stable_id: "ENSG00000139618".into(),
            chr: "13".to_string(),
            start: 32315086,
            end: 32400268,
            strand: Strand::Forward,
            biotype: "protein_coding".into(),
            source: "ensembl_havana".to_string(),
            gene_symbol: Some("BRCA2".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:1101".into()),
            gene_phenotype: Some(3),
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: Some("CCDS9344.1".to_string()),
            protein_id: Some("ENSP00000369497".into()),
            refseq: Some("NM_000059.4".to_string()),
            swissprot: Some("P51587.242".to_string()),
            trembl: Some("-".to_string()),
            uniparc: None,
            exons: Vec::new(),
            introns: Vec::new(),
            cdna_coding_start: Some(228),
            cdna_coding_end: Some(10_484),
            coding_region_start: Some(32_316_461),
            coding_region_end: Some(32_398_770),
            translation_start: None,
            translation_end: None,
            translation: Some(vep_core::transcript::Translation {
                stable_id: "ENSP00000369497".to_string(),
                version: None,
                db_id: None,
                start: 1,
                end: 1,
                start_exon_index: 0,
                end_exon_index: 0,
                seq: None,
            }),
            cdna_sequence: None,
            protein_sequence: None,
            flags: Vec::new().into(),
            gencode_primary: false,
            gencode_basic: false,
            gene_version: Some(17),
            uniprot_isoform: None,
            attributes: Vec::new(),
            vefc: Some(vep_core::transcript::TranscriptVEFC::default()),
            derived: Default::default(),
            id_rank: 0,
        }
    }

    /// The transcript keys resolve to the cache's values: `is_canonical` is
    /// defined on the canonical transcript alone, a missing symbol is
    /// undefined, the strand is `1` or `-1`, the object keys are defined and
    /// `translation` only on a coding transcript.
    #[test]
    fn transcript_keys_read_the_transcript() {
        let mut tx = canonical_transcript();
        let filter = |s: &str| FilterSet::transcript_filter(&[s]).unwrap();
        // `gene_version` resolves to `_gene_version` through the `_` prefix
        // (`FilterSet.pm:562-564`), the gene's version as the loader reads it.
        assert!(
            filter("gene_version exists and gene_version is 17 and gene_version > 16")
                .evaluate(&tx)
        );
        assert!(
            filter("dbID is 4201747 and cdna_coding_start < 300 and cdna_coding_end > 10000")
                .evaluate(&tx)
        );
        assert!(filter("translation exists and translation match Translation and gene exists and slice exists and attributes exists and trans_exon_array exists and variation_effect_feature_cache exists and vep_lazy_loaded is 1").evaluate(&tx));
        // The cache leaves the coding region and the description undefined.
        assert!(
            filter("coding_region_start nex and coding_region_end nex and description nex")
                .evaluate(&tx)
        );
        // Perl numifies a reference as its address, and `get_input` reads an
        // unblessed `ARRAY(0x…)` as `0` against a non-zero number alone.
        assert!(
            filter("translation > 0 and slice >= 1 and attributes > 0 and not attributes > 1")
                .evaluate(&tx)
        );
        assert!(!filter("translation is 1 or attributes is ARRAY or gene is Bio::EnsEMBL::Gene=HASH or attributes is #trans_exon_array").evaluate(&tx));
        assert!(filter(
            "gene ne Bio::EnsEMBL::Gene=HASH and gene match Gene=HASH and attributes nre HASH"
        )
        .evaluate(&tx));
        assert!(filter("biotype is protein_coding and is_canonical").evaluate(&tx));
        assert!(filter("strand is 1 and start > 32000000 and end < 33000000").evaluate(&tx));
        assert!(
            filter("gene_symbol is BRCA2 and gene_hgnc_id > 1100 and version is 8").evaluate(&tx)
        );
        assert!(filter("ccds match CCDS and protein exists and trembl is -").evaluate(&tx));
        assert!(filter("gene_phenotype > 2 and source is ensembl_havana").evaluate(&tx));
        assert!(filter("uniparc nex and uniprot_isoform nex").evaluate(&tx));
        assert!(filter("gene_stable_id is ENSG00000139618 and gene_symbol_source is HGNC and swissprot match P51587 and refseq is NM_000059.4").evaluate(&tx));
        tx.canonical = false;
        tx.gene_symbol = None;
        tx.strand = Strand::Reverse;
        tx.translation = None;
        tx.vefc = None;
        assert!(filter("not is_canonical and gene_symbol nex and strand is -1").evaluate(&tx));
        assert!(!filter("gene_symbol ne BRCA2").evaluate(&tx));
        assert!(filter("translation nex and variation_effect_feature_cache nex").evaluate(&tx));
    }

    /// `uniprot_isoform` resolves to the cache's `_uniprot_isoform`
    /// (`FilterSet.pm:562-564`), the UniProt isoform ids of the translation
    /// joined by `,` (`AnnotationSource/Database/Transcript.pm:679-683`, `-`
    /// for a translation without one): a value with several ids passes a unit
    /// any id passes (487-498 split the input at a `,` after a digit), and the
    /// `-` sentinel is a defined value.
    #[test]
    fn uniprot_isoform_reads_the_transcripts_isoform_ids() {
        let mut tx = canonical_transcript();
        tx.uniprot_isoform = Some("P51587-2,P51587-3".to_string());
        let filter = |s: &str| FilterSet::transcript_filter(&[s]).unwrap();
        assert!(filter("uniprot_isoform exists").evaluate(&tx));
        assert!(filter("uniprot_isoform is P51587-3").evaluate(&tx));
        assert!(filter("uniprot_isoform is P51587-2").evaluate(&tx));
        assert!(filter("uniprot_isoform match P51587").evaluate(&tx));
        assert!(filter("uniprot_isoform in P51587-3,Q00001-1").evaluate(&tx));
        assert!(!filter("uniprot_isoform is P51587-4").evaluate(&tx));
        assert!(!filter("uniprot_isoform nex").evaluate(&tx));
        assert!(filter("_uniprot_isoform is P51587-3 and _uniprot_i exists").evaluate(&tx));
        tx.uniprot_isoform = Some("-".to_string());
        assert!(filter("uniprot_isoform exists and uniprot_isoform is -").evaluate(&tx));
        assert!(!filter("uniprot_isoform is P51587-2").evaluate(&tx));
        tx.uniprot_isoform = None;
        assert!(filter("uniprot_isoform nex").evaluate(&tx));
        assert!(!filter("uniprot_isoform is -").evaluate(&tx));
    }

    /// Two filter strings are two units under the root, joined by `and`
    /// (`parse_filters`, 217-221); no filter string at all passes everything.
    #[test]
    fn repeated_filters_are_joined_by_and() {
        let data = [("biotype", "protein_coding"), ("strand", "1")];
        let both =
            FilterSet::new(&["biotype is protein_coding", "strand is 1"], &keys(&data)).unwrap();
        assert!(both.evaluate(&hash(&data)));
        let other =
            FilterSet::new(&["biotype is protein_coding", "strand is -1"], &keys(&data)).unwrap();
        assert!(!other.evaluate(&hash(&data)));
        let none: [&str; 0] = [];
        assert!(FilterSet::new(&none, &keys(&data))
            .unwrap()
            .evaluate(&hash(&data)));
    }

    /// The filter `--exclude_predicted` appends (`Config.pm:566-571`), compiled
    /// as the binary compiles it: the tokenizer drops the `^` (`parse_filters`,
    /// 223), so `/X._/i` drops the RefSeq predicted models `XM_011545924.3` and
    /// `XR_001755219.2` and keeps `ENST00000427445`; a filter before it still
    /// applies, the two being joined by `and`.
    #[test]
    fn exclude_predicted_filter_drops_the_refseq_predicted_models() {
        let predicted = "not stable_id match ^X._";
        let with_id = |id: &str| {
            let mut tx = canonical_transcript();
            tx.stable_id = id.into();
            tx
        };
        let alone = FilterSet::transcript_filter(&[predicted]).unwrap();
        assert!(!alone.evaluate(&with_id("XM_011545924.3")));
        assert!(!alone.evaluate(&with_id("XR_001755219.2")));
        assert!(alone.evaluate(&with_id("ENST00000427445")));
        let both = FilterSet::transcript_filter(&["biotype is protein_coding", predicted]).unwrap();
        assert!(both.evaluate(&with_id("ENST00000427445")));
        assert!(!both.evaluate(&with_id("XM_011545924.3")));
        let mut non_coding = with_id("ENST00000427445");
        non_coding.biotype = "lncRNA".into();
        assert!(!both.evaluate(&non_coding));
    }
}
