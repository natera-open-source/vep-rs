// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's filter grammar, the one `--transcript_filter` reads
//! (`Bio::EnsEMBL::VEP::FilterSet`, ensembl-vep release/116.3).
//!
//! A filter string is units of `field operator value` joined by `and` and `or`,
//! grouped by parentheses and inverted by `not`. Ensembl VEP source citations below name
//! `FilterSet.pm` at release/116.3.
//!
//! * The operators are `eq` (`is`, `=`), `ne` (`!=`), `gt` (`>`), `lt` (`<`),
//!   `gte` (`>=`), `lte` (`<=`), `ex` (`exists`, `defined`), `nex`, `re`
//!   (`match`, `matches`, `regex`), `nre`, `in` and `is_child` (the synonym
//!   table at lines 100-117, the predicates at 704-748). A unit without an
//!   operator tests that the field is defined; `ne` without a value is `nex`
//!   (`finish_filter_node`, 370-424).
//! * `and` multiplies and `or` adds the running result, left to right with no
//!   precedence (`evaluate`, 441-512): `a or b and c` is `(a or b) and c`.
//! * A value starting with `#` names another field, whose value is compared.
//! * `in` takes a comma-separated list, or the name of a file with one value
//!   per line, read as bytes (385-419).
//! * `re` and `nre` are case-insensitive, unanchored regular expressions
//!   (`/$pattern/i`, 714-715), read as Perl reads them (`filter_regex`) and
//!   compiled with the `regex` crate: a construct the crate lacks
//!   (back-references, look-around) is an error naming the pattern. The
//!   pattern runs over the bytes of the value as Perl's runs over a byte
//!   string (`Matcher`): `.` is one byte, a class holds bytes, and `\w`,
//!   `\s`, `\d` and the case fold are ASCII; a `\p{}` class, a `\x{}`,
//!   `\o{}` or octal escape above 255 or a `\N{U+}` escape switches the
//!   pattern to Unicode rules over the bytes read as Latin-1, as Perl's
//!   default rules switch it, and the case fold is then full: `ss` matches
//!   `ß`; a `u`, `a` or `aa` flag switches its own group alone.
//! * Every other text rule is Perl's over a byte string as well: a word splits
//!   on ASCII blanks, a field name's lower and upper case change its ASCII
//!   letters alone, the prefix search folds ASCII case alone, and the numeric
//!   extraction's `\w` and `\d` are ASCII.
//! * `is_child` passes when the value matches the field as `re` does, else when
//!   the field names a descendant of the value in the Sequence Ontology
//!   (722-748; the ancestry is [`vep_core::so_ontology`]).
//! * A field name resolves to a key of the data exactly, then lower-cased,
//!   upper-cased, prefixed with `_`, then as the unique case-insensitive prefix
//!   of one key; a prefix of several keys is an error (`get_input`, 533-613).
//!   [`FilterSet::new`] resolves every name once, when the filter is parsed,
//!   against a fixed key set; [`FilterSet::deferred`] resolves a name on each
//!   evaluation against the keys of the data at hand, through Ensembl VEP's
//!   synonym cache (546-548, 553-583): a name read through a synonym keeps that
//!   synonym for the rest of the run, and under
//!   [`FilterSet::limit_synonym_search`] a name no key matches is recorded as
//!   undefined and never searched again. The prefix branch reads the value it
//!   resolves to on the first evaluation as well, where Ensembl VEP records the
//!   synonym and reads nothing (570-571).
//! * A field value of the form `text(number)` or `text:number` is compared
//!   through its number when the filter value is numeric (589-610), and a value
//!   holding several entries joined by `,` or `&` after a digit passes when any
//!   entry passes (487-498).

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::sync::{LazyLock, Mutex, MutexGuard};

use regex::Regex;
use vep_core::coordinate::Strand;
use vep_core::so_ontology;
use vep_core::transcript::{BamEditStatus, Transcript};

use crate::filter_input::{self, encode_line};
use crate::filter_regex::{self, Refusal, Rules};

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
/// the data does not define. A deferred filter also asks which keys the data
/// holds: a key may be held and undefined, as a hash entry set to `undef` is.
pub trait FilterData {
    fn value(&self, key: &str) -> Option<Cow<'_, str>>;

    /// Whether the data holds `key`, defined or not (`exists`, 538).
    fn has_key(&self, key: &str) -> bool {
        self.value(key).is_some()
    }

    /// Every key the data holds, the set the prefix search runs over (568).
    fn keys(&self) -> Vec<&str> {
        Vec::new()
    }
}

/// The keys of a transcript as Ensembl VEP's cache stores them, the fields a
/// `--transcript_filter` can name and the set its prefixes resolve over. The
/// first twenty-one are the scalar fields of `AnnotationSource/Database/Transcript.pm`
/// and of the transcript itself; `dbID`, `cdna_coding_start`, `cdna_coding_end`,
/// `translation` and `_gene_version` are read from the transcript too. `slice`,
/// `attributes`, `_gene`, `_trans_exon_array` and `_variation_effect_feature_cache`
/// hold an object in the cache and `_vep_lazy_loaded` is `1` on every transcript
/// the filter sees, so each is defined and only `exists` and `nex` read it
/// meaningfully. `description`, `coding_region_start` and `coding_region_end`
/// are undefined, as the cache stores them. `_source_cache` is `Ensembl` or
/// `RefSeq` on a merged cache's transcripts and undefined elsewhere, where
/// Ensembl VEP's key set lacks it (`MergeVEP.pm:103, 110` stamp it on the
/// region files both caches held), so `_s` is a prefix of two keys here on
/// every cache and names neither, as it names neither in Ensembl VEP on a
/// stamped transcript; `_sw` and `_so` name one each. `_gene_symbol_id` and
/// `_bam_edit_status` are a RefSeq or merged cache's alike: the EntrezGene
/// identifier a RefSeq gene's symbol is taken from
/// (`Database/Transcript.pm:563`) and the outcome of editing a model against
/// the cache's alignment file (`AnnotationType/Transcript.pm:601, 621`),
/// undefined on every other transcript, where Ensembl VEP's key set lacks
/// them; so `_gene_symbol_` is a prefix of two keys here on every cache, where
/// Ensembl VEP resolves it to `_gene_symbol_source` on a transcript without
/// the identifier. `display_xref`, `external_db`, `external_display_name`,
/// `external_name` and `external_status` are left out: the cache holds them,
/// undefined, only on the transcripts without a gene symbol, and with them the
/// prefix `e` would be ambiguous instead of naming `end`, as Ensembl VEP
/// resolves it on a transcript with a symbol.
pub const TRANSCRIPT_KEYS: [&str; 37] = [
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
    "_gene_symbol_id",
    "_gene_hgnc_id",
    "_ccds",
    "_refseq",
    "_protein",
    "_swissprot",
    "_trembl",
    "_uniparc",
    "_uniprot_isoform",
    "_gene_phenotype",
    "_source_cache",
    "_bam_edit_status",
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
            "_gene_symbol_id" => self.gene_symbol_id.as_deref().map(Cow::Borrowed),
            "_gene_hgnc_id" => self.hgnc_id.as_deref().map(Cow::Borrowed),
            "_ccds" => self.ccds.as_deref().map(Cow::Borrowed),
            "_refseq" => self.refseq.as_deref().map(Cow::Borrowed),
            "_protein" => self.protein_id.as_deref().map(Cow::Borrowed),
            "_swissprot" => self.swissprot.as_deref().map(Cow::Borrowed),
            "_trembl" => self.trembl.as_deref().map(Cow::Borrowed),
            "_uniparc" => self.uniparc.as_deref().map(Cow::Borrowed),
            "_uniprot_isoform" => self.uniprot_isoform.as_deref().map(Cow::Borrowed),
            "_gene_phenotype" => self.gene_phenotype.map(|n| Cow::Owned(n.to_string())),
            "_source_cache" => self
                .source_cache
                .as_ref()
                .map(|source| Cow::Borrowed(source.column_text())),
            // The cache's spelling, `ok` or `failed`; the `BAM_EDIT` column
            // upper-cases it (OutputFactory.pm:1560-1562), the filter does not.
            "_bam_edit_status" => self.bam_edit_status.map(|status| {
                Cow::Borrowed(match status {
                    BamEditStatus::Ok => "ok",
                    BamEditStatus::Failed => "failed",
                })
            }),
            "dbID" => self.db_id.map(|n| Cow::Owned(n.to_string())),
            "cdna_coding_start" => self.cdna_coding_start.map(|n| Cow::Owned(n.to_string())),
            "cdna_coding_end" => self.cdna_coding_end.map(|n| Cow::Owned(n.to_string())),
            // Perl stringifies a reference as its class and address and numifies
            // it as the address, so a value of that shape keeps `is`, `ne`, a
            // regex and a numeric comparison reading as they do in Ensembl VEP
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
    Text(String, Option<Matcher>),
    /// The `in` operator's set, with its first entry as written.
    List {
        first: Option<String>,
        set: HashSet<String>,
    },
    /// `#name`: the field `name`, with the key it resolved to at parse time
    /// (`None` when the key set has no such key, or under deferred resolution).
    Field { name: String, key: Option<String> },
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
    /// The key the field name resolved to at parse time; `None` leaves the
    /// value undefined, or defers the resolution to each evaluation.
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

/// Ensembl VEP's synonym cache (`synonyms`, 630-634): the key each written
/// field name was resolved to, `None` for a name recorded as undefined under
/// `limit_synonym_search` (652-656).
#[derive(Debug, Clone, Default)]
struct Synonyms {
    map: HashMap<String, Option<String>>,
    limit: bool,
}

/// A compiled filter, evaluated against one [`FilterData`] at a time.
#[derive(Debug)]
pub struct FilterSet {
    root: Node,
    /// The synonym cache of a deferred filter; `None` when every field name
    /// was resolved at parse time.
    synonyms: Option<Mutex<Synonyms>>,
}

impl Clone for FilterSet {
    fn clone(&self) -> Self {
        FilterSet {
            root: self.root.clone(),
            synonyms: self
                .locked_synonyms()
                .map(|synonyms| Mutex::new(Synonyms::clone(&synonyms))),
        }
    }
}

const INCOMPLETE_PARENTHESES: &str = "Error parsing filter string - incomplete parentheses sets?";

impl FilterSet {
    /// Parses `filters`, every string a unit joined to the others by `and`
    /// (`parse_filters`, 209-310), and resolves each field name against
    /// `keys`, the data's key set. No filter string at all is a filter
    /// everything passes.
    pub fn new<S: AsRef<str>>(filters: &[S], keys: &[&str]) -> Result<FilterSet, FilterError> {
        Ok(FilterSet {
            root: parse_all(filters, Some(keys))?,
            synonyms: None,
        })
    }

    /// Parses `filters` as [`FilterSet::new`] does and leaves every field name
    /// to be resolved against the keys of each [`FilterData`] it is evaluated
    /// on, through the synonym cache, as `filter_vep` resolves a name against
    /// each line's hash.
    pub fn deferred<S: AsRef<str>>(filters: &[S]) -> Result<FilterSet, FilterError> {
        Ok(FilterSet {
            root: parse_all(filters, None)?,
            synonyms: Some(Mutex::new(Synonyms::default())),
        })
    }

    /// A `--transcript_filter`: the filter strings bound to [`TRANSCRIPT_KEYS`].
    pub fn transcript_filter<S: AsRef<str>>(filters: &[S]) -> Result<FilterSet, FilterError> {
        FilterSet::new(filters, &TRANSCRIPT_KEYS)
    }

    /// The synonym cache of a deferred filter, locked; `None` for a filter
    /// resolved at parse time.
    fn locked_synonyms(&self) -> Option<MutexGuard<'_, Synonyms>> {
        self.synonyms.as_ref().map(|synonyms| {
            synonyms
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
        })
    }

    /// Whether a deferred filter records a name no key matches as undefined,
    /// so the search never runs again for it (`limit_synonym_search`,
    /// 652-656); off, the search runs on every evaluation. No effect on a
    /// filter resolved at parse time.
    pub fn limit_synonym_search(&self, limit: bool) {
        if let Some(mut synonyms) = self.locked_synonyms() {
            synonyms.limit = limit;
        }
    }

    /// The `is_child` values of every unit, as written.
    pub fn is_child_values(&self) -> Vec<&str> {
        fn walk<'a>(node: &'a Node, out: &mut Vec<&'a str>) {
            if node.children.is_empty() && node.operator == Some(Operator::IsChild) {
                out.extend(node.value.as_deref());
            }
            for child in &node.children {
                walk(child, out);
            }
        }
        let mut out = Vec::new();
        walk(&self.root, &mut out);
        out
    }

    /// Whether a unit reads `key`, as its field or as a `#field` value, in a
    /// filter whose names were resolved at parse time.
    pub(crate) fn names_key(&self, key: &str) -> bool {
        fn walk(node: &Node, key: &str) -> bool {
            let as_value =
                matches!(&node.bound, Some(Value::Field { key: Some(k), .. }) if k == key);
            node.key.as_deref() == Some(key)
                || as_value
                || node.children.iter().any(|child| walk(child, key))
        }
        walk(&self.root, key)
    }

    /// Whether `data` passes the filter; a deferred filter's resolution error
    /// reads as a failure.
    pub fn evaluate(&self, data: &dyn FilterData) -> bool {
        self.try_evaluate(data).unwrap_or(false)
    }

    /// Whether `data` passes the filter, or the error an evaluation meets: a
    /// deferred field name that is a prefix of several of the data's keys
    /// (`get_input`, 574-580), or a `#field` pattern Perl or the dialect
    /// refuses, which Ensembl VEP compiles, and dies on, at the match (714).
    pub fn try_evaluate(&self, data: &dyn FilterData) -> Result<bool, FilterError> {
        if self.root.children.is_empty() {
            return Ok(true);
        }
        let mut synonyms = self.locked_synonyms();
        let mut lookup = Lookup {
            data,
            synonyms: synonyms.as_deref_mut(),
        };
        Ok(evaluate(&self.root, &mut lookup)?.unwrap_or(false))
    }
}

/// Every filter string parsed under one root, each a unit joined to the
/// others by `and` (`parse_filters`, 217-221); the field names resolved
/// against `keys` when given, else left for evaluation.
fn parse_all<S: AsRef<str>>(filters: &[S], keys: Option<&[&str]>) -> Result<Node, FilterError> {
    let mut root = Node::new(Logic::And);
    for filter in filters {
        parse_into(&mut root, filter.as_ref(), keys)?;
    }
    Ok(root)
}

/// Parses one filter string as a unit under `root` (`parse_filters`,
/// 217-307): the words [`next_token`] yields fill a unit's `not`, field,
/// operator and value in turn, `and` and `or` start the next unit, and the
/// parentheses open and close groups.
fn parse_into(root: &mut Node, filter: &str, keys: Option<&[&str]>) -> Result<(), FilterError> {
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

/// Whether `c` ends a word: a blank to Perl's `\s` on a byte string, or a
/// parenthesis.
fn is_separator(c: char) -> bool {
    filter_input::perl_space(c) || c == '(' || c == ')'
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
/// name resolved against `keys` when they are given. A group's own field and
/// value are never read, so only a unit resolves and compiles them; a group's
/// `in` list is still read as Ensembl VEP reads it.
fn finish(mut node: Node, keys: Option<&[&str]>) -> Result<Node, FilterError> {
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
    let resolve_now = |field: &str| match keys {
        Some(keys) => resolve(field, keys),
        None => Ok(None),
    };
    node.key = resolve_now(node.field.as_deref().unwrap_or(""))?;
    node.bound = match (&node.value, operator) {
        (_, Operator::In) => node.bound,
        (Some(v), _) if v.starts_with('#') => Some(Value::Field {
            name: v[1..].to_string(),
            key: resolve_now(&v[1..])?,
        }),
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

/// A compiled pattern, run over the bytes of a value as Perl runs
/// `/$pattern/i` over a byte string. Under Perl's default rules a byte above
/// ASCII has no case and is no `\w`, `\s` or `\d`, so the pattern runs over
/// the bytes with Unicode off. A `\p{}` or `\P{}` class, a `\x{}`, `\o{}` or
/// octal escape above 255 or a `\N{U+}` escape switches Perl to Unicode
/// rules, under which each byte is the Latin-1 character of its value, and a
/// `u`, `a` or `aa` flag switches its own group; such a pattern runs over the
/// UTF-8 of that reading with Unicode on ([`filter_regex`]).
#[derive(Debug, Clone)]
enum Matcher {
    Bytes(regex::bytes::Regex),
    Latin1(regex::bytes::Regex),
}

impl Matcher {
    fn is_match(&self, value: &str) -> bool {
        let bytes = encode_line(value);
        match self {
            Matcher::Bytes(regex) => regex.is_match(&bytes),
            // The Latin-1 reading of ASCII bytes is the bytes themselves.
            Matcher::Latin1(regex) if bytes.is_ascii() => regex.is_match(&bytes),
            Matcher::Latin1(regex) => regex.is_match(latin1(&bytes).as_bytes()),
        }
    }
}

/// Each byte read as the character of its value.
fn latin1(bytes: &[u8]) -> String {
    bytes.iter().map(|&b| char::from(b)).collect()
}

/// A case-insensitive, unanchored regular expression (`/$pattern/i`) over the
/// bytes `pattern` was read from.
fn compile(pattern: &str) -> Result<Matcher, FilterError> {
    let unsupported =
        |message: &str| FilterError(format!("pattern \"{pattern}\" is not supported: {message}"));
    let bytes = encode_line(pattern);
    let translated = filter_regex::translate(&bytes).map_err(|refusal| match refusal {
        Refusal::Perl(message) => FilterError(format!("{message} in regex m/{pattern}/")),
        Refusal::Unsupported(message) => unsupported(&message),
    })?;
    let regex = regex::bytes::RegexBuilder::new(&format!("(?i){}", translated.pattern))
        .unicode(translated.rules == Rules::Latin1)
        .build()
        .map_err(|e| {
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
            unsupported(&message)
        })?;
    Ok(match translated.rules {
        Rules::Bytes => Matcher::Bytes(regex),
        Rules::Latin1 => Matcher::Latin1(regex),
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
        // Perl's `split` in Ensembl VEP drops the empty fields a trailing comma leaves.
        let mut parts: Vec<&str> = list.split(',').collect();
        while parts.last() == Some(&"") {
            parts.pop();
        }
        parts.into_iter().map(str::to_string).collect()
    } else {
        let path = filter_input::path_of(list);
        if !path.exists() {
            return Err(FilterError(format!("Could not find/parse list {list}")));
        }
        let contents = std::fs::read(&path)
            .map_err(|_| FilterError(format!("Could not read from file {list}")))?;
        // Ensembl VEP reads the file line by line, so a final newline ends the last
        // entry and opens no empty one; `split_inclusive` keeps that count.
        // Each line then loses its carriage returns and every byte `\v`
        // matches on a byte string (406): 0x0A to 0x0D and 0x85.
        contents
            .split_inclusive(|b| *b == b'\n')
            .map(|line| {
                let kept: Vec<u8> = line
                    .iter()
                    .copied()
                    .filter(|b| !matches!(b, 0x0A..=0x0D | 0x85))
                    .collect();
                filter_input::decode_line(&kept).into_owned()
            })
            .collect()
    };
    let first = entries.first().cloned();
    Ok((first, entries.into_iter().collect()))
}

/// The data under evaluation with, for a deferred filter, the synonym cache
/// its field names resolve through.
struct Lookup<'a> {
    data: &'a dyn FilterData,
    synonyms: Option<&'a mut Synonyms>,
}

impl<'a> Lookup<'a> {
    /// The value of the field written `name`, resolved at parse time to `key`
    /// (`None` when no key matched, or when the resolution is deferred): the
    /// data's value under `key`, or under the key the deferred search finds.
    fn input(
        &mut self,
        name: &str,
        key: Option<&str>,
    ) -> Result<Option<Cow<'a, str>>, FilterError> {
        match &mut self.synonyms {
            Some(synonyms) => deferred_input(name, self.data, synonyms),
            None => Ok(key.and_then(|key| self.data.value(key))),
        }
    }
}

/// `get_input`'s lookup (538-586) against the data at hand: the exact key
/// reads the value; else a recorded synonym reads through it (an undefined
/// synonym reads nothing); else the lower-cased, upper-cased and `_`-prefixed
/// names each record themselves and read; else the one key the name is a
/// case-insensitive prefix of is recorded and read, several are an error, and
/// none records the name as undefined when the search is limited. The case
/// changes are Perl's `lc` and `uc` on a byte string: ASCII letters alone.
fn deferred_input<'a>(
    field: &str,
    data: &'a dyn FilterData,
    synonyms: &mut Synonyms,
) -> Result<Option<Cow<'a, str>>, FilterError> {
    if data.has_key(field) {
        return Ok(data.value(field));
    }
    if let Some(synonym) = synonyms.map.get(field) {
        return Ok(synonym.as_deref().and_then(|key| data.value(key)));
    }
    for candidate in [
        field.to_ascii_lowercase(),
        field.to_ascii_uppercase(),
        format!("_{field}"),
    ] {
        if data.has_key(&candidate) {
            let value = data.value(&candidate);
            synonyms.map.insert(field.to_string(), Some(candidate));
            return Ok(value);
        }
    }
    match prefixed(field, &data.keys())? {
        Some(key) => {
            let value = data.value(&key);
            synonyms.map.insert(field.to_string(), Some(key));
            Ok(value)
        }
        None => {
            if synonyms.limit {
                synonyms.map.insert(field.to_string(), None);
            }
            Ok(None)
        }
    }
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
        .or_else(|| exact(&field.to_ascii_lowercase()))
        .or_else(|| exact(&field.to_ascii_uppercase()))
        .or_else(|| exact(&format!("_{field}")))
    {
        return Ok(Some(key));
    }
    prefixed(field, keys)
}

/// The one key of `keys` that `field` is a case-insensitive prefix of
/// (567-572); several are the error of 574-580, none is `None`. Ensembl VEP tests
/// each key against `/^$field/i`, so the name is a pattern over the key's
/// bytes with ASCII case folded: a name without a regex metacharacter is
/// compared as the bytes it was read from, any other is compiled.
fn prefixed(field: &str, keys: &[&str]) -> Result<Option<String>, FilterError> {
    let has_metacharacter = field.bytes().any(|b| br"\.+*?()|[]{}^$".contains(&b));
    let matches: Vec<&str> = if has_metacharacter {
        let pattern = compile(&format!("^{field}"))?;
        keys.iter()
            .copied()
            .filter(|k| pattern.is_match(k))
            .collect()
    } else {
        let prefix = encode_line(field);
        keys.iter()
            .copied()
            .filter(|k| {
                let key = encode_line(k);
                key.len() >= prefix.len() && key[..prefix.len()].eq_ignore_ascii_case(&prefix)
            })
            .collect()
    };
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
/// does not invert; `Err` is a deferred resolution's error.
fn evaluate(node: &Node, lookup: &mut Lookup<'_>) -> Result<Option<bool>, FilterError> {
    let mut result = if node.children.is_empty() {
        match unit(node, lookup)? {
            Some(verdict) => verdict,
            None => return Ok(None),
        }
    } else {
        let mut acc = true;
        for child in &node.children {
            let verdict = evaluate(child, lookup)?.unwrap_or(false);
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
    Ok(Some(result))
}

// The classes of `get_input`'s patterns over a byte string: `\w` and `\d` are
// ASCII, so no byte of a UTF-8 letter or digit is one.
static EXTRACT: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?-u)([\w\.\-]+)?:?\(?([\-\d\.e]*)\)?").unwrap());
static NUMERIC_VALUE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?-u)^[\-\d\.e]+$").unwrap());
static NUMERIC_TEXT: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?-u)^\-?\d+\.?\d*(e\-?\d+)?$").unwrap());
static MULTI_VALUE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"[0-9](,|&)").unwrap());

/// The verdict of one unit (464-504); `None` is the early return of 477.
fn unit(node: &Node, lookup: &mut Lookup<'_>) -> Result<Option<bool>, FilterError> {
    let operator = node.operator.unwrap_or(Operator::Ex);
    let field = node.field.as_deref().unwrap_or("");
    let raw = lookup.input(field, node.key.as_deref())?;
    // The filter value as `get_input` reads it: a word, a list's first entry,
    // or the `#name` of a field, which stays the written word there (472-476).
    let written = match &node.bound {
        Some(Value::Text(text, _)) => Some(text.as_str()),
        Some(Value::List { first, .. }) => first.as_deref(),
        Some(Value::Field { .. }) | None => node.value.as_deref(),
    };
    let input = raw.map(|v| extract(v, written, field));

    let compared: Option<Cow<'_, str>> = match &node.bound {
        Some(Value::Field { name, key }) => match lookup.input(name, key.as_deref())? {
            Some(value) => Some(value),
            None => return Ok(None),
        },
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
        return Ok(Some(false));
    }
    let input = input.as_deref();
    let predicate = |piece: Option<&str>| -> Result<bool, FilterError> {
        Ok(match operator {
            Operator::Eq => piece.unwrap_or("") == value.unwrap_or(""),
            Operator::Ne => piece.unwrap_or("") != value.unwrap_or(""),
            Operator::Gt => perl_number(piece) > perl_number(value),
            Operator::Lt => perl_number(piece) < perl_number(value),
            Operator::Gte => perl_number(piece) >= perl_number(value),
            Operator::Lte => perl_number(piece) <= perl_number(value),
            Operator::Ex => piece.is_some(),
            Operator::Nex => piece.is_none(),
            Operator::Re => matches(regex, value, piece)?,
            Operator::Nre => !matches(regex, value, piece)?,
            Operator::In => list.is_some_and(|set| set.contains(piece.unwrap_or(""))),
            Operator::IsChild => {
                matches(regex, value, piece)?
                    || value
                        .and_then(so_ontology::descendants)
                        .is_some_and(|names| {
                            names
                                .iter()
                                .any(|name| name.eq_ignore_ascii_case(piece.unwrap_or("")))
                        })
            }
        })
    };
    // Several values joined by `,` or `&` after a digit pass when one passes
    // (487-498); the pattern runs only on a value holding one of the two.
    if let Some(caps) = input
        .filter(|i| i.contains([',', '&']))
        .and_then(|i| MULTI_VALUE.captures(i))
    {
        let delimiter = caps.get(1).map_or(",", |m| m.as_str());
        let mut pieces: Vec<&str> = input.unwrap_or("").split(delimiter).collect();
        while pieces.last() == Some(&"") {
            pieces.pop();
        }
        for piece in pieces {
            if predicate(Some(piece))? {
                return Ok(Some(true));
            }
        }
        return Ok(Some(false));
    }
    predicate(input).map(Some)
}

/// Whether `piece` matches the unit's regular expression: the compiled one, or
/// a `#field` value compiled here, whose refusal is the unit's error, as Ensembl VEP
/// dies compiling it at the match (714).
fn matches(
    regex: Option<&Matcher>,
    value: Option<&str>,
    piece: Option<&str>,
) -> Result<bool, FilterError> {
    let piece = piece.unwrap_or("");
    match (regex, value) {
        (Some(re), _) => Ok(re.is_match(piece)),
        (None, Some(v)) => Ok(compile(v)?.is_match(piece)),
        (None, None) => Ok(false),
    }
}

/// `get_input`'s numeric extraction (589-610): a value shaped `text(number)` or
/// `text:number` compared against a numeric filter value is its number (its
/// text when that is itself a number); any other pairing keeps the value whole.
/// `CELL_TYPE` and the HGVS fields are exempt, as Ensembl VEP's condition names them.
fn extract<'a>(input: Cow<'a, str>, value: Option<&str>, field: &str) -> Cow<'a, str> {
    if field == "CELL_TYPE" || field.contains("HGVS") {
        return input;
    }
    // Only a numeric filter value reads a number out of the input (Perl's
    // `$value &&` reads "" and "0" as false), and the pattern's number group
    // is empty unless a `:` or `(` precedes it: its first group takes every
    // leading `[\w.-]` character, and the group's class lies inside that one.
    // The pattern is skipped where it would leave the input as it is.
    let numeric = value.is_some_and(|v| !v.is_empty() && v != "0" && NUMERIC_VALUE.is_match(v));
    if !numeric || !input.contains([':', '(']) {
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
    if NUMERIC_TEXT.is_match(text) {
        Cow::Owned(text.to_string())
    } else {
        Cow::Owned(number.to_string())
    }
}

/// A string as Perl's numeric comparison reads it: a stringified reference
/// (`Class=HASH(0xaddr)`, `ARRAY(0xaddr)`) is its address, as Perl numifies
/// the reference itself; anything else is its leading number after ASCII
/// blanks, `inf`, `infinity` or `nan` in any case with an optional sign, else
/// 0.
fn perl_number(s: Option<&str>) -> f64 {
    static LEADING_NUMBER: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(
            r"(?-u)^\s*(?:([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)|([+-]?)((?i:inf|nan)))",
        )
        .unwrap()
    });
    static REFERENCE: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"^(?:[\w:]+=)?(?:HASH|ARRAY)\(0x([0-9a-f]+)\)$").unwrap());
    let Some(s) = s else {
        return 0.0;
    };
    // The common shape, digits with at most one point, is its own number;
    // the patterns below read the rest.
    if let Some(number) = plain_decimal(s) {
        return number;
    }
    if let Some(address) = REFERENCE.captures(s).and_then(|caps| caps.get(1)) {
        return u64::from_str_radix(address.as_str(), 16).map_or(0.0, |n| n as f64);
    }
    let Some(caps) = LEADING_NUMBER.captures(s) else {
        return 0.0;
    };
    if let Some(number) = caps.get(1) {
        return number.as_str().parse::<f64>().unwrap_or(0.0);
    }
    if caps[3].eq_ignore_ascii_case("nan") {
        return f64::NAN;
    }
    if &caps[2] == "-" {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    }
}

/// `s` as a number when it is one or more ASCII digits with at most one
/// `.` among them, which Perl reads whole; `None` for every other string.
fn plain_decimal(s: &str) -> Option<f64> {
    let mut digits = 0;
    let mut points = 0;
    for &b in s.as_bytes() {
        match b {
            b'0'..=b'9' => digits += 1,
            b'.' => points += 1,
            _ => return None,
        }
    }
    if digits == 0 || points > 1 {
        return None;
    }
    s.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use vep_core::transcript::SourceCache;

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

    /// A pattern runs over the bytes of the value as Perl's `/$pattern/i` runs
    /// over a byte string: `.` is one byte, so the UTF-8 `é` is two; a class
    /// holds bytes; the case fold, `\w` and `\s` are ASCII; a backslash before
    /// a byte above ASCII is that byte. A `\p{}` class or a `\x{}` above 255
    /// switches the pattern to Unicode rules over the bytes read as Latin-1,
    /// under which `\w` matches the byte 0xE9 and `É` folds to `é`, as Perl's
    /// default rules switch. A braced `\x{}` is read as Perl reads it: blanks
    /// at either end and an underscore before a digit are skipped, the first
    /// other non-hex character ends the value, so a sign or empty braces give
    /// the NUL byte, and a surrogate or a value past U+10FFFF compiles and
    /// matches no byte. An unbraced `\x` takes the hex digits present, two at
    /// most: `\xE` before `$` or `z` is the byte 0x0E, `\x` before `g`, `$` or
    /// `Z?` is the NUL byte, which the `?` never binds, and `\xE9` stays the
    /// byte 0xE9. A blank never survives the filter string's tokenizer, so a
    /// pattern holding one reaches the compiler through a `#field` value (each
    /// verdict is Perl 5.34's on the same bytes).
    #[test]
    fn patterns_run_over_the_bytes_of_the_value() {
        let utf8 = [("id", "r\u{e9}")];
        let latin1 = [("id", "r\u{10FFE9}")];
        assert!(!passes("id match r.$", &utf8));
        assert!(passes("id match r.{2}$", &utf8));
        assert!(passes("id match r.$", &latin1));
        assert!(!passes("id match r[\u{e9}]$", &utf8));
        assert!(passes(r"id match r[\xC3][\xA9]", &utf8));
        assert!(passes(r"id match r\xE9", &latin1));
        assert!(!passes("id match \u{c9}", &utf8));
        assert!(!passes(r"id match r\xC9", &latin1));
        assert!(passes("id match R", &utf8));
        assert!(!passes(r"id match r\w", &utf8));
        assert!(!passes(r"id match r\w", &latin1));
        assert!(!passes(r"id match r\s", &[("id", "r\u{a0}")]));
        assert!(passes("id match r\\\u{e9}$", &utf8));
        assert!(passes("id match r\u{e9}+$", &utf8));
        assert!(passes(r"id match r\p{L}$", &latin1));
        assert!(!passes(r"id match r\p{L}$", &utf8));
        assert!(passes(r"id match r\p{L}\P{L}$", &utf8));
        assert!(passes(r"id match r\w\p{L}?$", &latin1));
        assert!(passes(r"id match R\xC9\p{L}?", &latin1));
        assert!(passes(r"id match r\x{100}?\xE9", &latin1));
        assert!(passes(r"id match r\x{00E9}$", &latin1));
        assert!(!passes(r"id match R\x{C9}$", &latin1));
        assert!(passes(r"id match r\s\p{L}?", &[("id", "r\u{10FFA0}")]));
        assert!(passes(r"id match r\x{0_00E9}$", &latin1));
        assert!(passes(r"id match r\x{_00E9}$", &latin1));
        assert!(passes(r"id match r\x{00E9_}$", &latin1));
        assert!(passes(r"id match r\x{00E9x}$", &latin1));
        assert!(passes(r"id match r[\x{0_00E9}]$", &latin1));
        assert!(!passes(r"id match r\x{E__9}$", &latin1));
        assert!(!passes(r"id match r\x{+00E9}$", &latin1));
        assert!(!passes(r"id match r\x{}$", &latin1));
        assert!(!passes(r"id match r\x{g}$", &latin1));
        assert!(!passes(r"id match r[\x{+00E9}]$", &latin1));
        assert!(passes(r"id match r\x{1_00}?\w$", &latin1));
        assert!(!passes(r"id match r\x{+100}?\w$", &latin1));
        assert!(passes(r"id match R\x{C9}\x{_100}?$", &latin1));
        assert!(passes(r"id match r\x{D800}?\w$", &latin1));
        assert!(passes(r"id match r\x{110000}?\w$", &latin1));
        assert!(passes(r"id match r[\x{100}-\x{110000}]?\w$", &latin1));
        let via_field = |pattern: &'static str| {
            passes("id match #pat", &[("id", "r\u{10FFE9}"), ("pat", pattern)])
        };
        assert!(via_field(r"r\x{ 00E9}$"));
        assert!(via_field(r"r\x{00E9 }$"));
        assert!(via_field(r"r\x{ 00E9 }$"));
        assert!(via_field("r\\x{\t00E9}$"));
        assert!(via_field(r"r\x{00E9 x}$"));
        assert!(!via_field(r"r\x{ }$"));
        assert!(via_field(r"r\x{ 100}?\xE9"));
        assert!(via_field(r"r\x{100 }?\w$"));
        assert!(!via_field(r"r\x{+100}?\w$"));
        assert!(passes(r"id match r\xE$", &[("id", "r\u{0e}")]));
        assert!(passes(r"id match r\xg", &[("id", "r\u{0}g")]));
        assert!(passes(r"id match r\x$", &[("id", "r\u{0}")]));
        assert!(passes(r"id match r\xEz$", &[("id", "r\u{0e}z")]));
        assert!(!passes(r"id match r\xZ?E", &[("id", "rE")]));
        assert!(passes(r"id match r\xE9$", &latin1));
        assert!(passes(r"id match r\xE9z$", &[("id", "r\u{10FFE9}z")]));
    }

    /// The numeric extraction's pattern yields an empty number group on every
    /// value without a `:` or `(`, whatever its first character, so such a
    /// value is read whole against a numeric filter value, and a non-numeric
    /// filter value reads every value whole: the premises of the paths that
    /// skip the pattern.
    #[test]
    fn extraction_reads_a_value_without_colon_or_parenthesis_whole() {
        for value in [
            "0.5",
            "-5",
            "1e-5",
            "abc-1",
            "x)5",
            "x 5",
            "\u{e9}5",
            "\u{10FFE9}5",
            ".5",
            "5.",
            "e",
            "-",
            "",
            "x,5",
            "1&2",
            "a=5",
            "5x",
            "x5",
        ] {
            let caps = EXTRACT.captures(value).unwrap();
            assert_eq!(caps.get(2).map_or("", |m| m.as_str()), "", "{value:?}");
            assert_eq!(
                extract(Cow::Borrowed(value), Some("1"), "F").as_ref(),
                value
            );
        }
        assert_eq!(
            extract(Cow::Borrowed("x(0.5)"), Some("1"), "F").as_ref(),
            "0.5"
        );
        assert_eq!(
            extract(Cow::Borrowed("HGNC:5"), Some("4"), "F").as_ref(),
            "5"
        );
        assert_eq!(extract(Cow::Borrowed("5:6"), Some("4"), "F").as_ref(), "5");
        assert_eq!(
            extract(Cow::Borrowed("x(0.5)"), Some("y"), "F").as_ref(),
            "x(0.5)"
        );
        assert_eq!(
            extract(Cow::Borrowed("x(0.5)"), Some("0"), "F").as_ref(),
            "x(0.5)"
        );
        assert_eq!(
            extract(Cow::Borrowed("x(0.5)"), None, "F").as_ref(),
            "x(0.5)"
        );
        assert_eq!(
            extract(Cow::Borrowed("x(0.5)"), Some("1"), "HGVSc").as_ref(),
            "x(0.5)"
        );
    }

    /// Perl's numeric reading of a byte string: ASCII blanks lead, a no-break
    /// space does not, and `Inf`, `Infinity`, `-inf` and `NaN` in any case
    /// are the floating-point values, with `NaN` ordered against nothing.
    #[test]
    fn perl_number_reads_blanks_infinities_and_nan_as_perl_does() {
        assert_eq!(perl_number(Some(" 7")), 7.0);
        assert_eq!(perl_number(Some("\x0B7")), 7.0);
        assert_eq!(perl_number(Some("\u{a0}7")), 0.0);
        assert_eq!(perl_number(Some("Inf")), f64::INFINITY);
        assert_eq!(perl_number(Some("Infinity")), f64::INFINITY);
        assert_eq!(perl_number(Some("infinityx")), f64::INFINITY);
        assert_eq!(perl_number(Some("-inf")), f64::NEG_INFINITY);
        assert_eq!(perl_number(Some("+INF")), f64::INFINITY);
        assert!(perl_number(Some("NaN")).is_nan());
        assert!(perl_number(Some("-nan")).is_nan());
        assert_eq!(perl_number(Some("1e400")), f64::INFINITY);
        assert_eq!(perl_number(Some("0x10")), 0.0);
        assert_eq!(perl_number(Some("\u{663}")), 0.0);
        assert_eq!(perl_number(Some("0.01")), 0.01);
        assert_eq!(perl_number(Some("1000")), 1000.0);
        assert_eq!(perl_number(Some("7.")), 7.0);
        assert_eq!(perl_number(Some(".5")), 0.5);
        assert_eq!(perl_number(Some("1.2.3")), 1.2);
        assert_eq!(perl_number(Some(".")), 0.0);
        assert_eq!(perl_number(Some("")), 0.0);
        assert_eq!(perl_number(Some("-3")), -3.0);
        assert_eq!(perl_number(Some("00012")), 12.0);
        let data = [("S", "NaN"), ("T", "Inf")];
        assert!(!passes("S > 1 or S < 1 or S >= 1 or S <= 1", &data));
        assert!(passes("T > 1e308 and not T < 1", &data));
    }

    /// The text rules over a byte string: a field name's lower and upper case
    /// change its ASCII letters alone, the prefix search folds ASCII case
    /// alone and reads the name as a pattern, and a word of the filter string
    /// splits on ASCII blanks alone (each verdict is Perl 5.34's).
    #[test]
    fn text_rules_fold_ascii_alone() {
        let data = [("\u{c9}K", "1"), ("ek", "2"), ("D.P", "3")];
        assert!(!passes("\u{e9}k is 1", &data));
        assert!(passes("\u{c9}k is 1", &data));
        assert!(passes("EK is 2", &data));
        assert!(!passes("\u{e9} is 1", &data));
        assert!(passes("\u{c9} is 1", &data));
        assert!(passes("d.p is 3", &data));
        assert!(passes("e. is 2", &data));
        let tx = canonical_transcript();
        assert!(FilterSet::transcript_filter(&["bio.ype is protein_coding"])
            .unwrap()
            .evaluate(&tx));
        assert!(
            FilterSet::transcript_filter(&["biotype$ is protein_coding"])
                .unwrap()
                .evaluate(&tx)
        );
        assert_eq!(
            FilterSet::new(&["e[ is 2"], &keys(&data))
                .unwrap_err()
                .to_string(),
            "Unmatched [ in regex m/^e[/"
        );
        assert_eq!(
            FilterSet::new(&["ID\u{a0}is rE"], &["ID"])
                .unwrap_err()
                .to_string(),
            "No such operator \"rE\""
        );
        assert!(passes("id\x0Bis 1", &[("id", "1")]));
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
        // `_s` is a prefix of `_swissprot` and `_source_cache`, so it names
        // neither; `_sw` and `_so` name one each.
        assert_eq!(
            FilterSet::transcript_filter(&["_s exists"])
                .unwrap_err()
                .to_string(),
            "Multiple data fields found matching query field \"_s\": _swissprot, _source_cache"
        );
        let set = FilterSet::transcript_filter(&["_sw exists"]).unwrap();
        assert_eq!(set.root.children[0].key.as_deref(), Some("_swissprot"));
        let set = FilterSet::transcript_filter(&["_so exists"]).unwrap();
        assert_eq!(set.root.children[0].key.as_deref(), Some("_source_cache"));
        // `_gene_symbol_` is a prefix of the symbol's source and of its
        // EntrezGene identifier; `_gene_symbol_s` and `_gene_symbol_i` name one
        // each, and `_b` the edit status alone.
        assert_eq!(
            FilterSet::transcript_filter(&["_gene_symbol_ exists"])
                .unwrap_err()
                .to_string(),
            "Multiple data fields found matching query field \"_gene_symbol_\": _gene_symbol_source, _gene_symbol_id"
        );
        let set = FilterSet::transcript_filter(&["_gene_symbol_i exists"]).unwrap();
        assert_eq!(set.root.children[0].key.as_deref(), Some("_gene_symbol_id"));
        let set = FilterSet::transcript_filter(&["_b exists"]).unwrap();
        assert_eq!(
            set.root.children[0].key.as_deref(),
            Some("_bam_edit_status")
        );
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
    /// (728-747), which `vep` reaches with no ontology adaptor and dies on
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

    /// A pattern the `regex` crate cannot express is an error naming the
    /// pattern and the construct, and a pattern Perl refuses is an error with
    /// Perl's message: a reference to a group the pattern lacks, an unclosed
    /// class and a trailing backslash. Perl's `\Z` and `\K` escapes compile.
    /// Parentheses never reach a pattern, since the tokenizer reads them as
    /// grouping, so a pattern with a group opens one instead.
    #[test]
    fn unsupported_regex_syntax_exits_naming_the_pattern() {
        let error = |filter: &str| {
            FilterSet::transcript_filter(&[filter])
                .unwrap_err()
                .to_string()
        };
        assert_eq!(
            error(r"stable_id match ENST\X"),
            r#"pattern "ENST\X" is not supported: unrecognized escape sequence"#
        );
        assert_eq!(
            error("stable_id match s++s$"),
            r#"pattern "s++s$" is not supported: the possessive quantifier ++ is not supported"#
        );
        assert_eq!(
            error(r"stable_id match \1"),
            r"Reference to nonexistent group in regex m/\1/"
        );
        assert_eq!(
            error("stable_id match ENST["),
            "Unmatched [ in regex m/ENST[/"
        );
        assert_eq!(
            error(r"stable_id match ENST\"),
            r"Trailing \ in regex m/ENST\/"
        );
        assert!(FilterSet::transcript_filter(&["stable_id match (ENST)"]).is_ok());
        // A `#field` pattern is compiled at the match: Perl's or the
        // dialect's refusal is the evaluation's error, a run that compiles
        // matches.
        let via_field = FilterSet::new(&["id match #pat"], &["id", "pat"]).unwrap();
        assert_eq!(
            via_field
                .try_evaluate(&hash(&[("id", "r"), ("pat", "r{01}")]))
                .unwrap_err()
                .to_string(),
            "Invalid quantifier in {,} in regex m/r{01}/"
        );
        assert_eq!(
            via_field
                .try_evaluate(&hash(&[("id", "r"), ("pat", "r(?=x)")]))
                .unwrap_err()
                .to_string(),
            "pattern \"r(?=x)\" is not supported: look-around, including look-ahead and look-behind, is not supported"
        );
        assert_eq!(
            via_field
                .try_evaluate(&hash(&[("id", "r"), ("pat", "r{1,2}+")]))
                .unwrap_err()
                .to_string(),
            "pattern \"r{1,2}+\" is not supported: the possessive quantifier {1,2}+ is not supported"
        );
        assert_eq!(
            via_field.try_evaluate(&hash(&[("id", "r"), ("pat", "r{1}")])),
            Ok(true)
        );
        assert!(FilterSet::transcript_filter(&[r"stable_id match ENST\Z"]).is_ok());
        assert!(FilterSet::transcript_filter(&[r"stable_id nre EN\KST"]).is_ok());
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
        let mut model = Transcript::default();
        model.stable_id = "ENST00000380152".into();
        model.version = Some(8);
        model.db_id = Some(4_201_747);
        model.gene_stable_id = "ENSG00000139618".into();
        model.chr = "13".to_string();
        model.start = 32315086;
        model.end = 32400268;
        model.strand = Strand::Forward;
        model.biotype = "protein_coding".into();
        model.source = "ensembl_havana".to_string();
        model.gene_symbol = Some("BRCA2".into());
        model.gene_symbol_source = Some("HGNC".into());
        model.gene_symbol_id = None;
        model.hgnc_id = Some("HGNC:1101".into());
        model.gene_phenotype = Some(3);
        model.canonical = true;
        model.mane_select = None;
        model.mane_plus_clinical = None;
        model.tsl = None;
        model.appris = None;
        model.ccds = Some("CCDS9344.1".to_string());
        model.protein_id = Some("ENSP00000369497".into());
        model.refseq = Some("NM_000059.4".to_string());
        model.swissprot = Some("P51587.242".to_string());
        model.trembl = Some("-".to_string());
        model.uniparc = None;
        model.exons = Vec::new();
        model.exon_ids = None;
        model.introns = Vec::new();
        model.cdna_coding_start = Some(228);
        model.cdna_coding_end = Some(10_484);
        model.coding_region_start = Some(32_316_461);
        model.coding_region_end = Some(32_398_770);
        model.translation_start = None;
        model.translation_end = None;
        model.translation = Some(vep_core::transcript::Translation {
            stable_id: "ENSP00000369497".to_string(),
            version: None,
            db_id: None,
            start: 1,
            end: 1,
            seq: None,
        });
        model.cdna_sequence = None;
        model.protein_sequence = None;
        model.flags = Vec::new().into();
        model.gencode_primary = false;
        model.gencode_basic = false;
        model.gene_version = Some(17);
        model.uniprot_isoform = None;
        model.attributes = Vec::new();
        model.vefc = Some(vep_core::transcript::TranscriptVEFC::default());
        model.source_cache = None;
        model.bam_edit_status = None;
        model.circular_length = None;
        model.region_files = None;
        model.derived = Default::default();
        model.id_rank = 0;
        model
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

    /// `_source_cache` is a transcript key, as it is in Ensembl VEP on a merged
    /// cache's transcripts (`MergeVEP.pm:103, 110` stamp `Ensembl` or `RefSeq`
    /// on the region files both caches held): `--transcript_filter
    /// "_source_cache is RefSeq"`, the Ensembl documentation's own example,
    /// keeps the RefSeq-stamped transcript alone, an unstamped transcript (an
    /// Ensembl cache's, or a one-side region file's) has the key undefined, and
    /// `source_cache` reaches it through the `_` rule (`FilterSet.pm:562-564`).
    /// The prefix `_s` is a prefix of two keys in this fixed key set (the test
    /// `ambiguous_prefix_is_an_error` states the error): Ensembl VEP's key set
    /// holds `_source_cache` on a stamped transcript alone, so there `_s`
    /// names `_swissprot` on every other transcript and errors on a stamped
    /// one; one key set for every cache makes the two-letter prefix ambiguous
    /// on every cache, and `_sw` and `_so` name the keys.
    #[test]
    fn source_cache_is_a_transcript_filter_key() {
        let filter = |s: &str| FilterSet::transcript_filter(&[s]).unwrap();
        let mut refseq = canonical_transcript();
        refseq.stable_id = "NM_000059.4".into();
        refseq.source_cache = Some(SourceCache::RefSeq);
        let mut ensembl = canonical_transcript();
        ensembl.source_cache = Some(SourceCache::Ensembl);
        let unstamped = canonical_transcript();
        let is_refseq = filter("_source_cache is RefSeq");
        assert!(is_refseq.evaluate(&refseq));
        assert!(!is_refseq.evaluate(&ensembl));
        assert!(!is_refseq.evaluate(&unstamped));
        assert!(filter("source_cache is Ensembl").evaluate(&ensembl));
        assert!(filter("_so exists").evaluate(&refseq));
        assert!(filter("_source_cache nex").evaluate(&unstamped));
        assert!(!filter("_source_cache exists").evaluate(&unstamped));
        assert!(filter("not _source_cache is RefSeq").evaluate(&ensembl));
    }

    /// `_gene_symbol_id` and `_bam_edit_status` are transcript keys, as they
    /// are in Ensembl VEP on a RefSeq or merged cache's transcripts: the
    /// EntrezGene identifier the symbol is taken from
    /// (`Database/Transcript.pm:563`) and the outcome of editing the model
    /// against the alignment file (`AnnotationType/Transcript.pm:601, 621`).
    /// A transcript without either (an Ensembl cache's, a RefSeq gene with no
    /// single EntrezGene entry, an unedited model) has the key undefined, as
    /// Ensembl VEP finds a key the hash lacks (`FilterSet.pm:538-586`): `exists`
    /// and a unit with a value drop it (481-483) and `nex` keeps it. `eq`
    /// compares the stored text (`ok`, 704), so `is OK` drops an edited model
    /// `match OK` keeps (714). `_gene_symbol_` is a prefix of two keys here on
    /// every cache, where Ensembl VEP's key set holds `_gene_symbol_id` on a
    /// RefSeq transcript alone (the test `ambiguous_prefix_is_an_error`).
    #[test]
    fn gene_symbol_id_and_bam_edit_status_are_transcript_filter_keys() {
        let filter = |s: &str| FilterSet::transcript_filter(&[s]).unwrap();
        let mut edited = canonical_transcript();
        edited.stable_id = "NR_109925.1".into();
        edited.bam_edit_status = Some(vep_core::transcript::BamEditStatus::Ok);
        let mut failed = canonical_transcript();
        failed.bam_edit_status = Some(vep_core::transcript::BamEditStatus::Failed);
        let unedited = canonical_transcript();
        assert!(filter("_bam_edit_status exists").evaluate(&edited));
        assert!(filter("_bam_edit_status is ok").evaluate(&edited));
        assert!(filter("bam_edit_status match OK and _b exists").evaluate(&edited));
        assert!(!filter("_bam_edit_status is OK").evaluate(&edited));
        assert!(!filter("_bam_edit_status nex").evaluate(&edited));
        assert!(filter("_bam_edit_status is failed").evaluate(&failed));
        assert!(!filter("_bam_edit_status is ok").evaluate(&failed));
        assert!(filter("_bam_edit_status nex").evaluate(&unedited));
        assert!(!filter("_bam_edit_status exists").evaluate(&unedited));
        assert!(!filter("_bam_edit_status is ok").evaluate(&unedited));
        assert!(filter("not _bam_edit_status is ok").evaluate(&unedited));
        let mut refseq = canonical_transcript();
        refseq.stable_id = "NM_003225.3".into();
        refseq.gene_symbol = Some("TFF1".into());
        refseq.gene_symbol_source = Some("EntrezGene".into());
        refseq.gene_symbol_id = Some("7031".into());
        let ensembl = canonical_transcript();
        assert!(filter("_gene_symbol_id exists").evaluate(&refseq));
        assert!(filter("_gene_symbol_id is 7031").evaluate(&refseq));
        assert!(filter("gene_symbol_id > 7000 and gene_symbol_id < 8000").evaluate(&refseq));
        assert!(filter("_gene_symbol_i exists and _gene_symbol_s is EntrezGene").evaluate(&refseq));
        assert!(!filter("_gene_symbol_id is 7157").evaluate(&refseq));
        assert!(!filter("_gene_symbol_id nex").evaluate(&refseq));
        assert!(filter("_gene_symbol_id nex").evaluate(&ensembl));
        assert!(!filter("_gene_symbol_id exists").evaluate(&ensembl));
        assert!(!filter("_gene_symbol_id is 7031").evaluate(&ensembl));
        assert!(filter("not _gene_symbol_id is 7031").evaluate(&ensembl));
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

    /// A measured pattern shape with `filter_vep`'s verdict over the fixture
    /// values (Perl 5.34, the pattern through `#field` or the filter string).
    struct Shape {
        name: &'static str,
        pattern: &'static [u8],
        values: &'static [(u32, &'static [u8])],
        verdict: Verdict,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Verdict {
        /// The positions whose ID the pattern matches.
        Rows(&'static [u32]),
        /// The pattern is refused.
        Refused,
    }

    #[rustfmt::skip]
    const FOLD_VALUES: &[(u32, &[u8])] = &[
        (100, br"ss"),
        (200, br"SS"),
        (300, br"sS"),
        (400, b"\xC3\x9F"),
        (500, b"\xDF"),
        (600, b"\xEF\xAC\x81"),
        (700, br"fi"),
        (800, br"FI"),
        (900, b"\xC5\xBF"),
        (1000, br"s"),
        (1100, br"S"),
        (1200, b"\xE2\x84\xAA"),
        (1300, br"k"),
        (1400, br"K"),
        (1500, b"\xE1\xBA\x9E"),
        (1600, br"st"),
        (1700, b"\xEF\xAC\x86"),
        (1800, br"ass"),
        (1900, b"a\xDFa"),
        (2000, br"assa"),
        (2100, b"\xB5"),
        (2200, br"ff"),
        (2300, b"\xEF\xAC\x80"),
        (2400, br"ffi"),
        (2500, b"\xEF\xAC\x83"),
        (2600, b"a\xDF"),
        (2700, b"\xDFs"),
        (2800, br"sss"),
        (2900, b"s\xDF"),
        (3000, br"ssss"),
        (3100, br"xss"),
        (3200, b"x\xDF"),
        (3300, br"xx"),
        (3400, b"\xDF\xDF"),
        (3500, b"\xC3\x9F\xC3\x9F"),
        (3600, b"\xE9"),
        (3700, b"\xC9"),
        (3800, br"i"),
        (3900, br"I"),
        (4000, br"t"),
        (4100, br"fl"),
        (4200, br"ffl"),
    ];
    #[rustfmt::skip]
    const ESC_VALUES: &[(u32, &[u8])] = &[
        (100, br"r{"),
        (200, br"r{x}"),
        (300, br"r{,3}"),
        (400, br"a"),
        (500, br"aaa"),
        (600, br"aaaa"),
        (700, br"A"),
        (800, b"r\x00"),
        (900, br"r0"),
        (1000, br"rA"),
        (1100, br"r{2}"),
        (1200, br"rr"),
        (1300, br"r{1,2}"),
        (1400, br"r{,}"),
        (1500, br"r{1,"),
        (1600, b"r\xE9"),
        (1700, b"r\x01"),
        (1800, br"{"),
        (1900, br"r}"),
        (2000, br"r"),
        (2100, br"rrr"),
        (2200, b"r\x081"),
        (2300, b"r\x08"),
        (2400, br"ra"),
        (2500, br"r{1}"),
        (2600, br"{x}"),
        (2700, b"r\xC4\x80"),
        (2800, b"r\x100"),
        (2900, br"rA0"),
        (3000, b"r\x001"),
    ];

    #[rustfmt::skip]
    const CHARSET_VALUES: &[(u32, &[u8])] = &[
        (100, br"ab"),
        (200, b"a\xE9"),
        (300, b"\xE9a"),
        (400, b"\xE9\xE9"),
        (500, b"a\xDF"),
        (600, b"\xDFa"),
        (700, b"\xDF\xDF"),
        (800, br"ss"),
        (900, b"ss\xDF"),
        (1000, br"ssss"),
        (1100, b"\xDFss"),
        (1200, br"a1"),
        (1300, br"1a"),
        (1400, b"a\xA0"),
        (1500, b"\xA0a"),
        (1600, b"a\xB2"),
        (1700, b"\xAA"),
        (1800, b"a\xAA"),
        (1900, b"\xC9\xC9"),
        (2000, b"a\xC9"),
        (2100, b"\xE9\xDF"),
        (2200, br"sS"),
        (2300, b"S\xDF"),
        (2400, b"\xDFS"),
        (2500, b"a\x85"),
        (2600, br"k"),
        (2700, br"K"),
        (2800, b"\xB5"),
        (2900, br"a"),
        (3000, b"\xE9"),
        (3100, b"\xDF"),
        (3200, b"s\xDF"),
        (3300, b"\xC9"),
        (3400, b"\xC9a"),
        (3500, b"a\xDFa"),
        (3600, br"assa"),
        (3700, b"\xDF\xE9"),
        (3800, br"xss"),
        (3900, b"\xAA\xAA"),
    ];
    #[rustfmt::skip]
    const NAMED_VALUES: &[(u32, &[u8])] = &[
        (100, br"A"),
        (200, br"AB"),
        (300, br"a"),
        (400, b"r\x00"),
        (500, br"r"),
        (600, br"rA"),
        (700, br"rAB"),
        (800, b"r\x04"),
        (900, b"r\x01"),
        (1000, b"r\xE9"),
        (1100, br"ra"),
        (1200, b"r\xC9"),
        (1300, br"rAA"),
        (1400, b"r\x041"),
        (1500, br"rABAB"),
        (1600, br"rab"),
        (1700, b"r\x00A"),
        (1800, b"rA\x00"),
        (1900, br"rA.B"),
        (2000, br"rA_B"),
        (2100, br"r4_1"),
        (2200, b"r\x04_1"),
        (2300, br"rB"),
        (2400, b"r\xC3\x89"),
    ];

    /// The measured shapes carried through a field: filter_vep's verdicts.
    #[rustfmt::skip]
    const FIELD_SHAPES: &[Shape] = &[
        Shape { name: "fold 1", pattern: br"(?:s)s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 2", pattern: br"(?:ss)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 3", pattern: br"(?:ss|zz)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 4", pattern: br"(?:\x{DF}|zz)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 5", pattern: br"(?u)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 6", pattern: br"(?a)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 7", pattern: br"(?aa)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "fold 8", pattern: br"(?d)ss\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 9", pattern: br"(?l)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 10", pattern: br"(?-i)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 1800, 2800, 3000, 3100]) },
        Shape { name: "fold 11", pattern: br"(?^u)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 1800, 2800, 3000, 3100]) },
        Shape { name: "fold 12", pattern: br"(?iu)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 13", pattern: br"(?n)(ss)$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "fold 14", pattern: br"(ss)$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "fold 15", pattern: br"(ss)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "fold 16", pattern: br"(s)(s)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "fold 17", pattern: br"^ss\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 18", pattern: br"^[^\x{DF}]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1000, 1100, 1300, 1400, 2100, 3600, 3700, 3800, 3900, 4000]) },
        Shape { name: "fold 19", pattern: br"^[^s]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 1300, 1400, 2100, 3600, 3700, 3800, 3900, 4000]) },
        Shape { name: "fold 20", pattern: br"^s{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300]) },
        Shape { name: "fold 21", pattern: br"^(?:ss){2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "fold 22", pattern: br"^(?i:ss)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 23", pattern: br"(?u)^\xDF$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 24", pattern: br"(?u)\xE9$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600, 3700]) },
        Shape { name: "fold 25", pattern: br"^\x{DF}{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "fold 26", pattern: br"^.\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 1000, 1100, 1300, 1400, 2100, 3600, 3700, 3800, 3900, 4000]) },
        Shape { name: "fold 27", pattern: br"^\w\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 1000, 1100, 1300, 1400, 2100, 3600, 3700, 3800, 3900, 4000]) },
        Shape { name: "fold 28", pattern: br"\o{400}?\xe9$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600, 3700]) },
        Shape { name: "fold 29", pattern: br"\400?\xe9$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600, 3700]) },
        Shape { name: "fold 30", pattern: br"^\x{DF}s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 2900]) },
        Shape { name: "fold 31", pattern: br"^(?:\x{DF}|ss)s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800]) },
        Shape { name: "fold 32", pattern: br"^s(?:s)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 33", pattern: br"^[\x{DF}s]{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 2700, 2800, 2900, 3000, 3400]) },
        Shape { name: "fold 34", pattern: br"^[^\x{DF}s]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1300, 1400, 2100, 3600, 3700, 3800, 3900, 4000]) },
        Shape { name: "fold 35", pattern: br"^(?-i)\x{DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "fold 36", pattern: br"^(?-i:\x{DF})\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "fold 37", pattern: br"^\N{U+DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 38", pattern: br"^s\x{DF}?\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1000, 1100, 2800, 2900]) },
        Shape { name: "fold 39", pattern: br"^[\x{DF}]{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "fold 40", pattern: br"^ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300]) },
        Shape { name: "fold 41", pattern: br"^\xDF$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "fold 42", pattern: br"(?u)^ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 43", pattern: br"^(?u)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 44", pattern: br"^(?:\x{DF})\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 45", pattern: br"^\x{DF}(?:)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 46", pattern: br"^s(?:)s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 47", pattern: br"^(?:s|x)s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300]) },
        Shape { name: "fold 48", pattern: br"^a(?:\x{DF}|ss)a\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1900, 2000]) },
        Shape { name: "fold 49", pattern: br"^(?:ss|s)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1000, 1100, 2900]) },
        Shape { name: "fold 50", pattern: br"^(?:s|ss)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1000, 1100, 2900]) },
        Shape { name: "fold 51", pattern: br"^\x{DF}|^zz\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 2700, 2800, 3000, 3400]) },
        Shape { name: "fold 52", pattern: br"^(?:ss)+\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 3000, 3400]) },
        Shape { name: "fold 53", pattern: br"^\x{1E9E}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 54", pattern: br"^x\x{DF}|ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100, 3200]) },
        Shape { name: "fold 55", pattern: br"^(?-i)ss|SS$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 2800, 3000]) },
        Shape { name: "fold 56", pattern: br"^[\x{DF}]\w\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 3400]) },
        Shape { name: "fold 57", pattern: br"^s[s]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 58", pattern: br"^\x{DF}{1,2}s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800]) },
        Shape { name: "fold 59", pattern: br"^ss{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2800]) },
        Shape { name: "fold 60", pattern: br"^a\x{DF}+\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1800, 2600]) },
        Shape { name: "fold 61", pattern: br"^\x{FB01}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[700, 800]) },
        Shape { name: "fold 62", pattern: br"^(?:f)(?:i)\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[700, 800]) },
        Shape { name: "fold 63", pattern: br"^f\x{FB01}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "fold 64", pattern: br"^\x{FB00}i\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "fold 65", pattern: br"^\x{FB03}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "fold 66", pattern: br"^[\x{FB03}]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "fold 67", pattern: br"^[\x{FB03}x]{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3300]) },
        Shape { name: "fold 68", pattern: br"^f\x{FB01}?\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "fold 69", pattern: br"^(?i)ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300]) },
        Shape { name: "fold 70", pattern: br"^(?u:ss)$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 71", pattern: br"^(?u)\xB5$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2100]) },
        Shape { name: "fold 72", pattern: br"^\x{3BC}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2100]) },
        Shape { name: "fold 73", pattern: br"^\p{Lu}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 1000, 1100, 1300, 1400, 2100, 3600, 3700, 3800, 3900, 4000]) },
        Shape { name: "fold 74", pattern: br"^\p{Lu}+$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 700, 800, 1000, 1100, 1300, 1400, 1600, 1800, 1900, 2000, 2100, 2200, 2400, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3600, 3700, 3800, 3900, 4000, 4100, 4200]) },
        Shape { name: "fold 75", pattern: br"^\x{DF}\p{L}$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 3400]) },
        Shape { name: "fold 76", pattern: br"^(?:ss\x{DF})\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "fold 77", pattern: br"^\x{DF}{2,}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "fold 78", pattern: br"^s{2,}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 2800, 3000]) },
        Shape { name: "fold 79", pattern: br"^s\x{17F}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "fold 80", pattern: br"^S\x{17F}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500]) },
        Shape { name: "brace 1", pattern: br"^(?:r){$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 2", pattern: br"^(r)\1$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200]) },
        Shape { name: "brace 3", pattern: br"^(r)\10$", values: ESC_VALUES, verdict: Verdict::Rows(&[2300]) },
        Shape { name: "brace 4", pattern: br"^r{01}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 5", pattern: br"^r{,}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1400]) },
        Shape { name: "brace 6", pattern: br"^r{1}{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 7", pattern: br"^r{1}{2}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 8", pattern: br"^r**$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 9", pattern: br"^r+?{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 10", pattern: br"^(?:r{)$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 11", pattern: br"^r{x,y}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 12", pattern: br"^r\{x,y}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 13", pattern: br"^r{x,$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 14", pattern: br"^r\{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 15", pattern: br"^r\}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1900]) },
        Shape { name: "brace 16", pattern: br"^\{$", values: ESC_VALUES, verdict: Verdict::Rows(&[1800]) },
        Shape { name: "brace 17", pattern: br"^r{1,2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2000]) },
        Shape { name: "brace 18", pattern: br"^r{,3}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2000, 2100]) },
        Shape { name: "brace 19", pattern: br"^r{0}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 20", pattern: br"^r{00}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 21", pattern: br"^r{0,02}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 22", pattern: br"^r{2,1}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 23", pattern: br"^r\o{}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 24", pattern: br"^r\o{x}$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "brace 25", pattern: br"^r\o{_101}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "brace 26", pattern: br"^r\o{101_}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "brace 27", pattern: br"^r\o{1__01}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "brace 28", pattern: br"^r\o{+101}$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "brace 29", pattern: br"^r\o{101", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 30", pattern: br"^r\o101$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 31", pattern: br"^r\o$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 32", pattern: br"^r\x{}$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "brace 33", pattern: br"^r\0{2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 34", pattern: br"^r\101{2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 35", pattern: br"^r\o{101}{2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 36", pattern: br"^r[\101-\132]$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 1200, 2400]) },
        Shape { name: "brace 37", pattern: br"^r[\o{101}-\o{132}]$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 1200, 2400]) },
        Shape { name: "brace 38", pattern: br"^r[\0-\1]$", values: ESC_VALUES, verdict: Verdict::Rows(&[800, 1700]) },
        Shape { name: "brace 39", pattern: br"^(r)(r)(r)(r)(r)(r)(r)(r)(r)(r)\10$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 40", pattern: br"^r\12$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 41", pattern: br"^r\8$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 42", pattern: br"^r\9$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 43", pattern: br"^r\18$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 44", pattern: br"^r\08$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 45", pattern: br"^r{,3}{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100, 1800]) },
        Shape { name: "brace 46", pattern: br"^r{1,2}?$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2000]) },
        Shape { name: "brace 47", pattern: br"^r{1,2}+$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2000]) },
        Shape { name: "brace 48", pattern: br"^r{{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 49", pattern: br"^r{1{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 50", pattern: br"^r{1,{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 51", pattern: br"^r{{1}$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 52", pattern: br"^r{1,2}}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1900]) },
        Shape { name: "brace 53", pattern: br"^r}}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 54", pattern: br"^[{]{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 55", pattern: br"^\p{Lu}{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "brace 56", pattern: br"^\p{Lu}{1}$", values: ESC_VALUES, verdict: Verdict::Rows(&[400, 700, 2000]) },
        Shape { name: "brace 57", pattern: br"^r\x{41}{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 58", pattern: br"^r\o{101}{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 59", pattern: br"^r\101{$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 60", pattern: br"^r{32766}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 61", pattern: br"^r{65534}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "brace 62", pattern: br"^r{65535}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 63", pattern: br"^r{99999}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "brace 64", pattern: br"^r{1,99999}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "u-o400tr", pattern: br"\o{400}?\xe9$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600, 3700]) },
        Shape { name: "u-400tr", pattern: br"\400?\xe9$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600, 3700]) },
        Shape { name: "e-rxcomma", pattern: br"r{x,$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 1", pattern: br"^(?u:\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 2", pattern: br"^\w(?u:\w)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 500, 800, 1200, 1300, 1800, 2000, 2200, 2300, 3200]) },
        Shape { name: "charset 3", pattern: br"^\w\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 4", pattern: br"^(?u:\w)\w\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 500, 600, 700, 800, 1200, 1300, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 3200, 3400, 3700, 3900]) },
        Shape { name: "charset 5", pattern: br"^(?a:\w)\w\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 500, 800, 1200, 1300, 1800, 2000, 2200, 2300, 3200]) },
        Shape { name: "charset 6", pattern: br"^\w(?a:\w)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 7", pattern: br"^(?a:\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 8", pattern: br"^(?aa:ss)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900]) },
        Shape { name: "charset 9", pattern: br"^(?aa:ss)\xdf\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900, 1000]) },
        Shape { name: "charset 10", pattern: br"^ss(?aa:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900]) },
        Shape { name: "charset 11", pattern: br"^(?u:ss)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900]) },
        Shape { name: "charset 12", pattern: br"^\xdf(?u:ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 1100]) },
        Shape { name: "charset 13", pattern: br"^a(?u:\s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1400, 2500]) },
        Shape { name: "charset 14", pattern: br"^a\s(?u:)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 15", pattern: br"^(?u:a)\s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 16", pattern: br"^(?a:a)\s\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1400, 2500]) },
        Shape { name: "charset 17", pattern: br"^a(?a:\s)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 18", pattern: br"^(?u:\xe9)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000, 3300]) },
        Shape { name: "charset 19", pattern: br"^(?u:a)\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200]) },
        Shape { name: "charset 20", pattern: br"^\xe9(?u:a)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[300]) },
        Shape { name: "charset 21", pattern: br"^(?aa:\xe9)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000, 3300]) },
        Shape { name: "charset 22", pattern: br"^(?aa:k)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700]) },
        Shape { name: "charset 23", pattern: br"^(?u:[[:alpha:]])[[:alpha:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 2200, 2400, 3400]) },
        Shape { name: "charset 24", pattern: br"^(?-u:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 25", pattern: br"^(?d:\w)\w\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 500, 600, 700, 800, 1200, 1300, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 3200, 3400, 3700, 3900]) },
        Shape { name: "charset 26", pattern: br"^(?u:\w)(?a:\w)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 27", pattern: br"^(?a:\w)(?u:\w)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 500, 800, 1200, 1300, 1800, 2000, 2200, 2300, 3200]) },
        Shape { name: "charset 28", pattern: br"^(?a)\w(?u:\w)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 500, 800, 1200, 1300, 1800, 2000, 2200, 2300, 3200]) },
        Shape { name: "charset 29", pattern: br"^(?u)\w(?a:\w)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 30", pattern: br"^\w(?u)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 500, 800, 1200, 1300, 1800, 2000, 2200, 2300, 3200]) },
        Shape { name: "charset 31", pattern: br"^\w(?a)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 32", pattern: br"^\w(?aa)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 33", pattern: br"^(?u:.)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 1500, 2200, 2400, 3400]) },
        Shape { name: "charset 34", pattern: br"^(?a:ss)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900]) },
        Shape { name: "charset 35", pattern: br"^(?a:\xdf)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 36", pattern: br"^ss(?a:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900, 1000]) },
        Shape { name: "charset 37", pattern: br"^(?aa:\xdf)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1100]) },
        Shape { name: "charset 38", pattern: br"^(?u:\xdf)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 39", pattern: br"^(?aa:ss)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000]) },
        Shape { name: "charset 40", pattern: br"^(?aa:s)s\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900]) },
        Shape { name: "charset 41", pattern: br"^s(?aa:s)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900]) },
        Shape { name: "charset 42", pattern: br"^(?u:s)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 43", pattern: br"^s(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 44", pattern: br"^(?u:s)(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 45", pattern: br"^(?u:ss|zz)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2300, 3100, 3200]) },
        Shape { name: "charset 46", pattern: br"^(?u:a)(?a:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "charset 47", pattern: br"^(?i:\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 48", pattern: br"^(?u:\W)\W$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 49", pattern: br"^\W(?u:\W)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 50", pattern: br"^(?u:\b)a\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200]) },
        Shape { name: "charset 51", pattern: br"^a\xe9(?u:\b)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200]) },
        Shape { name: "charset 52", pattern: br"^a\xe9\b$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 53", pattern: br"^a\xe9(?u:$)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200]) },
        Shape { name: "charset 54", pattern: br"^(?u:[\w])\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 55", pattern: br"^\w(?u:[\w])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 500, 800, 1200, 1300, 1800, 2000, 2200, 2300, 3200]) },
        Shape { name: "charset 56", pattern: br"^(?u:[^\W])[^\W]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 57", pattern: br"^(?u:\d)\d$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 58", pattern: br"^(?u:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 59", pattern: br"^(?aa:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3100]) },
        Shape { name: "charset 60", pattern: br"^(?a:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 61", pattern: br"^(?u:ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 62", pattern: br"^(?aa:ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 63", pattern: br"^(?u:)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 64", pattern: br"^ss(?u:)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 65", pattern: br"^(?u:x)|^ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3800]) },
        Shape { name: "charset 66", pattern: br"^ss|^(?u:x)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 900, 1000, 2200]) },
        Shape { name: "charset 67", pattern: br"^(?aa:x)|^ss\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100, 3800]) },
        Shape { name: "charset 68", pattern: br"^(?u:\xaa)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "charset 69", pattern: br"^(?u:\xb5)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2800]) },
        Shape { name: "charset 70", pattern: br"^(?u:[[:lower:]])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 71", pattern: br"^[[:lower:]](?u:)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 72", pattern: br"^(?a:[[:alpha:]])\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 73", pattern: br"^(?u:\p{L})$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 74", pattern: br"^(?a:\p{L})$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 75", pattern: br"^\p{L}(?a:\w)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 2200, 2400, 3400]) },
        Shape { name: "charset 76", pattern: br"^\w(?u:\x{100}?)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 77", pattern: br"^\w(?a:\x{100}?)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 78", pattern: br"^(?a:\x{100}?)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 79", pattern: br"^(?ua:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 80", pattern: br"^(?au:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 81", pattern: br"^(?uu:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 82", pattern: br"^(?aaa:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 83", pattern: br"^(?dd:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 84", pattern: br"^(?ad:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 85", pattern: br"^(?da:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 86", pattern: br"^(?-a:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 87", pattern: br"^(?-d:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 88", pattern: br"^(?-l:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 89", pattern: br"^(?aua:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 90", pattern: br"^(?^a:\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 91", pattern: br"^(?^:\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 92", pattern: br"^(?u)(?^i)\w\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200]) },
        Shape { name: "charset 93", pattern: br"^(?u:)[^\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3100, 3300]) },
        Shape { name: "charset 94", pattern: br"^(?u:)[\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 95", pattern: br"^(?u:)[a-z]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 96", pattern: br"^(?u:)[[:lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 97", pattern: br"^(?u:)[[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 98", pattern: br"^(?u:)[^a-z]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 99", pattern: br"^(?u:)[\xc0-\xde]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3300]) },
        Shape { name: "charset 100", pattern: br"^(?u:)[\xe0-\xfe]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 101", pattern: br"^(?u:)[\xdfa]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2900, 3100]) },
        Shape { name: "charset 102", pattern: br"^(?u:)\xe9+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400, 3000]) },
        Shape { name: "charset 103", pattern: br"^(?u:)[\w\xe9]+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 800, 1000, 1200, 1300, 2200, 2600, 2700, 2900, 3000, 3600, 3800]) },
        Shape { name: "charset 104", pattern: br"^(?u:)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3100]) },
        Shape { name: "charset 105", pattern: br"^(?u:)ss(?u:)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 106", pattern: br"^(?u:(?a)\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 107", pattern: br"^(?a:(?u)\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 108", pattern: br"^(?u:(?d)\w)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 109", pattern: br"^(?u)\w(?d)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 300, 600, 800, 1200, 1300, 2200, 2400, 3400]) },
        Shape { name: "charset 110", pattern: br"^(?u)\w(?d)\w\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 500, 600, 700, 800, 1200, 1300, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 3200, 3400, 3700, 3900]) },
        Shape { name: "charset 111", pattern: br"^(?u:k)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700]) },
        Shape { name: "charset 112", pattern: br"^(?aa:[\xdf])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3100]) },
        Shape { name: "charset 113", pattern: br"^(?u:[\xdf])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 114", pattern: br"^(?u:)[\xdf]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3100]) },
        Shape { name: "charset 115", pattern: br"^(?u:)\b\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 116", pattern: br"^(?u:)\B\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 117", pattern: br"^(?u:)[[:alpha:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 118", pattern: br"^(?u:)[[:^alpha:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 119", pattern: br"^(?u:)[^[:alpha:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 120", pattern: br"^(?u:)[^[:lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 121", pattern: br"^(?u:)[[:upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 122", pattern: br"^(?u:)[^\xe9a]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 3100, 3300]) },
        Shape { name: "charset 123", pattern: br"^(?u:)[^a-k]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 124", pattern: br"^(?u:)[\xe9-\xeb]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 125", pattern: br"^(?u:)[\xc9-\xcb]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3300]) },
        Shape { name: "charset 126", pattern: br"^(?u:)[k-\xb5]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800]) },
        Shape { name: "charset 127", pattern: br"^(?u:)[A-\xc9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3300]) },
        Shape { name: "charset 128", pattern: br"^(?u:)[^A-\xc9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000, 3100]) },
        Shape { name: "charset 129", pattern: br"^(?u:)[\xe9A]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2900, 3000]) },
        Shape { name: "charset 130", pattern: br"^(?u:)[^\xe9A-Z]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3100, 3300]) },
        Shape { name: "charset 131", pattern: br"^(?u:)(?i)[K]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700]) },
        Shape { name: "charset 132", pattern: br"^(?u:)(?-i)[k]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600]) },
        Shape { name: "charset 133", pattern: br"^(?u:)(?-i)\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 134", pattern: br"^(?aa:)\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 135", pattern: br"^(?aa:)[\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 136", pattern: br"^(?a:)[\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 137", pattern: br"^(?u:[[:^lower:]])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 138", pattern: br"^(?u:[^[:lower:]])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 139", pattern: br"^(?u:)[[:^upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 140", pattern: br"^(?u:)[\d\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 141", pattern: br"^(?u:)[^\d]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 142", pattern: br"^(?u:)[\xe9]+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400, 3000]) },
        Shape { name: "charset 143", pattern: br"^(?u:)[\xe9]{2}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400]) },
        Shape { name: "charset 144", pattern: br"^(?u:)[\xc9\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000, 3300]) },
        Shape { name: "charset 145", pattern: br"^(?u:)[\xdf\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000, 3100]) },
        Shape { name: "charset 146", pattern: br"^(?u:)[\xa0a]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2900]) },
        Shape { name: "charset 147", pattern: br"^(?u:)\xa0a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "charset 148", pattern: br"^(?u:)[s\xdf]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3100]) },
        Shape { name: "charset 149", pattern: br"^(?u:)[^s\xdf]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3300]) },
        Shape { name: "charset 150", pattern: br"^(?u:)[^\xdf]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3300]) },
        Shape { name: "charset 151", pattern: br"^(?u:)[^\xdf]+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 800, 1000, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2500, 2600, 2700, 2800, 2900, 3000, 3300, 3400, 3600, 3800, 3900]) },
        Shape { name: "charset 152", pattern: br"^(?u:)[s-\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 153", pattern: br"^(?u:)[^k]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 154", pattern: br"^(?u:)[kK]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700]) },
        Shape { name: "charset 155", pattern: br"^(?u:)[^kK]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 156", pattern: br"^(?u:)[^\W]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 157", pattern: br"^(?u:)[\W]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 158", pattern: br"^(?u:)[^\w\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3100, 3300]) },
        Shape { name: "charset 159", pattern: br"^(?u:[^\w])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 160", pattern: br"^(?u:)[^\xc0-\xde]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100]) },
        Shape { name: "charset 161", pattern: br"^(?u:)[^\xe0-\xfe]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3100, 3300]) },
        Shape { name: "charset 162", pattern: br"^(?u:)[^[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 163", pattern: br"^(?u:)[[:lower:][:upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 164", pattern: br"^(?u:)[^[:lower:][:upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 165", pattern: br"^(?u:)[\x{e9}]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 166", pattern: br"^(?u:)\x{e9}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 167", pattern: br"^(?u:)\351$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 168", pattern: br"^(?u:)\o{351}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 169", pattern: br"^(?u:)[^a]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 170", pattern: br"^(?u:)\xdf+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 3100]) },
        Shape { name: "charset 171", pattern: br"^(?u:)\xe9{2}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400]) },
        Shape { name: "charset 172", pattern: br"^(?u:)(?:\xe9)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 173", pattern: br"^(?u:)(?:\xe9|a)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2900, 3000]) },
        Shape { name: "charset 174", pattern: br"^(?u:)(?i:\xe9)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 175", pattern: br"^(?u:)(?s:\xe9)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "charset 176", pattern: br"^(?u:)\xe9|^(?u:)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[300, 400, 2100, 3000, 3100]) },
        Shape { name: "charset 177", pattern: br"^(?u:)(?:\xe9){2}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400]) },
        Shape { name: "charset 178", pattern: br"^(?u:)[\xc0-\xff]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3000, 3100, 3300]) },
        Shape { name: "charset 179", pattern: br"^(?u:)[\x80-\xff]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 180", pattern: br"^(?u:)[^\x80-\xff]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 181", pattern: br"^(?u:)[\x00-\xff]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 182", pattern: br"^(?u:)\xb5$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2800]) },
        Shape { name: "charset 183", pattern: br"^(?u:)\xff$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 184", pattern: br"^(?u:)\xaa$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "charset 185", pattern: br"^(?u:)\xc9\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 186", pattern: br"^(?u:)[[:word:]][[:word:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 1200, 1300, 2200]) },
        Shape { name: "charset 187", pattern: br"^(?u:)[[:punct:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 188", pattern: br"^(?u:)[[:space:]]a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 189", pattern: br"^(?u:)\s(?u:\s)a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 190", pattern: br"^(?u:)\h\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "charset 191", pattern: br"^(?u:)\Ha$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[300, 600, 1300, 3400]) },
        Shape { name: "charset 192", pattern: br"^(?u:)\v(?u:)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 193", pattern: br"^\N{U+DF}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 194", pattern: br"^\N{U+41}\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200, 2000]) },
        Shape { name: "charset 195", pattern: br"^\N{U+73.73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 196", pattern: br"^\N{U+DF.DF}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900, 1000, 1100]) },
        Shape { name: "charset 197", pattern: br"^[\N{U+DF}]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 198", pattern: br"^[[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 199", pattern: br"^[^[:lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 200", pattern: br"^[[:^upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 201", pattern: br"^[[:^alpha:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 202", pattern: br"^[^a-z]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 203", pattern: br"^[^[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 204", pattern: br"^[[:^lower:]]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 205", pattern: br"^[^[:lower:]]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 206", pattern: br"^[[:^upper:]]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 207", pattern: br"^(?u:[[:^upper:]])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 208", pattern: br"^(?aa:[[:^lower:]])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 209", pattern: br"^(?a:[[:lower:]])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 210", pattern: br"^[^[:lower:]]+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400, 700, 1700, 1900, 2100, 2800, 3000, 3100, 3300, 3700, 3900]) },
        Shape { name: "charset 211", pattern: br"^(?u:)[[:^word:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 212", pattern: br"^(?u:)[[:^digit:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 213", pattern: br"^(?u:)[[:^alpha:]]+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400, 700, 1700, 1900, 2100, 2800, 3000, 3100, 3300, 3700, 3900]) },
        Shape { name: "charset 214", pattern: br"^(?u:)\xdf{2}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700]) },
        Shape { name: "charset 215", pattern: br"^(?u:)[\xdf]{2}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700]) },
        Shape { name: "charset 216", pattern: br"^(?u:)(?:ss)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[900]) },
        Shape { name: "charset 217", pattern: br"^(?u:)(?:s)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 218", pattern: br"^(?u:)\N{U+41}\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200, 2000]) },
        Shape { name: "charset 219", pattern: br"^[[:lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 220", pattern: br"^[[:upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 221", pattern: br"^(?-i)[[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 222", pattern: br"^(?-i)[[:^upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 223", pattern: br"^(?-i)[^[:lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 224", pattern: br"^(?u:)(?-i)[[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 225", pattern: br"^(?-i)[[:^lower:]]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2700, 3300]) },
        Shape { name: "charset 226", pattern: br"^[[:^lower:][:^upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 227", pattern: br"^[[:^lower:]\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 228", pattern: br"^[[:^lower:]a]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 229", pattern: br"^[^[:^lower:]a]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700]) },
        Shape { name: "charset 230", pattern: br"^[\W]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 231", pattern: br"^[^\W]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 232", pattern: br"^[\W]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 233", pattern: br"^[^\W]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 234", pattern: br"^[\D]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 235", pattern: br"^[^\S]a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 236", pattern: br"^[^\D]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 237", pattern: br"^[\S\D]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 238", pattern: br"^[^\w]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 239", pattern: br"^[^\w\xe9]\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 240", pattern: br"^(?u:)[^\w]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 241", pattern: br"^(?u:)[^\d\xe9]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3100, 3300]) },
        Shape { name: "charset 242", pattern: br"^(?u:)[^[:digit:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2600, 2700, 2800, 2900, 3000, 3100, 3300]) },
        Shape { name: "charset 243", pattern: br"^(?u:)[^[:^digit:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 244", pattern: br"^(?u:s)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 245", pattern: br"^(?u:)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 246", pattern: br"^(?u:s)(?aa:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 247", pattern: br"^(?aa:s)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 248", pattern: br"^s(?aa:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 249", pattern: br"^(?a:s)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 250", pattern: br"^s(?a:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 251", pattern: br"^(?u:a)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 252", pattern: br"^ss(?u:a)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 253", pattern: br"^(?u:a)(?:ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 254", pattern: br"^(?u:a)s(?u:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 255", pattern: br"^(?u:s)[s]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 256", pattern: br"^(?u:[s])s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 257", pattern: br"^(?u:s)\x73$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 258", pattern: br"^(?u:s)s{1}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 259", pattern: br"^(?u:s)(?:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 260", pattern: br"^(?u:s)(s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 261", pattern: br"^(?u:x)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3800]) },
        Shape { name: "charset 262", pattern: br"^x(?u:)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3800]) },
        Shape { name: "charset 263", pattern: br"^(?u:)x(?u:)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3800]) },
        Shape { name: "charset 264", pattern: br"^(?u:x)s(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3800]) },
        Shape { name: "charset 265", pattern: br"^(?u:s)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 266", pattern: br"^(?u:)s(?u:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 267", pattern: br"^s(?u:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 268", pattern: br"^(?u:)[s]s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 269", pattern: br"^s(?u:)s(?u:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 270", pattern: br"^(?u:)s(?:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 271", pattern: br"^(?u:)s(?i)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 272", pattern: br"^(?u:)s(?#c)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 273", pattern: br"^(?u:)s()s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 274", pattern: br"^(?u:)s\x73$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 275", pattern: br"^(?u:a)s\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 276", pattern: br"^(?u:\xdf)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 277", pattern: br"^(?u:\xdf)s(?u:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 278", pattern: br"^(?u:a)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "charset 279", pattern: br"^(?a:a)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 280", pattern: br"^(?aa:a)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 281", pattern: br"^(?u:ab)|^ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 800, 2200]) },
        Shape { name: "charset 282", pattern: br"^(?u:a)(?:x|ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 283", pattern: br"^(?u:a)(?:ss|x)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 284", pattern: br"^(?u:a)(?:ss)a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3600]) },
        Shape { name: "charset 285", pattern: br"^(?u:a)ss(?u:a)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3600]) },
        Shape { name: "charset 286", pattern: br"^(?u:a)s(?u:s)a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3500, 3600]) },
        Shape { name: "charset 287", pattern: br"^(?u:x)|ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 1000, 1100, 2200, 3800]) },
        Shape { name: "charset 288", pattern: br"^(?u:x)ss|\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[500, 700, 900, 2100, 2300, 3100, 3200, 3800]) },
        Shape { name: "charset 289", pattern: br"^(?u:ss)(?:)ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 290", pattern: br"^(?u:ss)\x73\x73$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 291", pattern: br"^(?u:a)[s][s]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 292", pattern: br"^(?u:a)s[s]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 293", pattern: br"^(?u:a)\xdf(?u:a)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[3500]) },
        Shape { name: "charset 294", pattern: br"^(?u:\xe9)a$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[300, 3400]) },
        Shape { name: "charset 295", pattern: br"^a(?u:\xe9)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[200, 2000]) },
        Shape { name: "charset 296", pattern: br"^(?u:\xe9)\xe9$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[400]) },
        Shape { name: "charset 297", pattern: br"^(?u:\xdf)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900]) },
        Shape { name: "charset 298", pattern: br"^(?u:s)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 3200]) },
        Shape { name: "charset 299", pattern: br"^\xdf(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 300", pattern: br"^(?u:\xdf)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 301", pattern: br"^s(?u:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 3200]) },
        Shape { name: "charset 302", pattern: br"^(?u:a)\N{U+41}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 303", pattern: br"^(?u:\w)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2400]) },
        Shape { name: "charset 304", pattern: br"^(?u:s)\w$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 305", pattern: br"^(?u:s)[sS]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 306", pattern: br"^[sS](?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 307", pattern: br"^(?u:s).$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2300, 3200]) },
        Shape { name: "charset 308", pattern: br"^(?u:s)s?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 309", pattern: br"^(?u:s)s*$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 1000, 2200]) },
        Shape { name: "charset 310", pattern: br"^(?u:s)s+$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 1000, 2200]) },
        Shape { name: "charset 311", pattern: br"^(?u:s)(?:s|t)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 312", pattern: br"^(?u:s)(?:s|)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2300, 3200]) },
        Shape { name: "charset 313", pattern: br"^(?u:s)(?:s|s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2300, 3200]) },
        Shape { name: "charset 314", pattern: br"^(?u:s)(?:s)(?:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 315", pattern: br"^(?u:s)(?:(?:s))$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 316", pattern: br"^(?u:s)(?u:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 317", pattern: br"^(?u:s)(?aa:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 318", pattern: br"^(?u:s)(?a:)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 319", pattern: br"^(?u:s)(?d:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 320", pattern: br"^(?u:s)(?^:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "charset 321", pattern: br"^(?u:s)(?^i:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 322", pattern: br"^(?u:s)(?i:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 323", pattern: br"^(?u:s)(?s:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 324", pattern: br"^(?u:s)(?-i:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "charset 325", pattern: br"^(?u:s)\Ks$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 326", pattern: br"^(?u:s)\bs$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 327", pattern: br"^(?u:s)\x{73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 328", pattern: br"^(?u:s)\o{163}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 329", pattern: br"^(?u:s)\163$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 330", pattern: br"^(?u:s)s\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 331", pattern: br"^(?lu:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 332", pattern: br"^(?ul:\w)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 333", pattern: br"^(?u:ss)(?n)(ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1000, 1100]) },
        Shape { name: "charset 334", pattern: br"^(?u)ss(?n)(ss)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900, 1000, 1100]) },
        Shape { name: "charset 335", pattern: br"^\N{}$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 336", pattern: br"^r\N{u+41}$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 337", pattern: br"^(?aa)[[:lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 338", pattern: br"^(?a)[[:^lower:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[1700, 2800, 3000, 3100, 3300]) },
        Shape { name: "charset 339", pattern: br"^(?u:)(?i)[[:upper:]]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2600, 2700, 2900]) },
        Shape { name: "charset 340", pattern: br"^(?u:\xdf)(?u:\xdf)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900, 1000, 1100]) },
        Shape { name: "charset 341", pattern: br"^(?u:s)(?u:s)\xdf$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900]) },
        Shape { name: "charset 342", pattern: br"^(?u:[s])(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 343", pattern: br"^(?u:s)(?u:[s])$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 344", pattern: br"^(?u:s)(?aa:s)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 345", pattern: br"^(?u:s)(?a:s)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 346", pattern: br"^(?aa:s)(?aa:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 347", pattern: br"^(?a:s)(?a:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 348", pattern: br"^(?a:s)(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 349", pattern: br"^(?u:\w)(?d:\w)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 500, 600, 700, 800, 1200, 1300, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 3200, 3400, 3700, 3900]) },
        Shape { name: "charset 350", pattern: br"^(?u:)(?u:s)s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 351", pattern: br"^(?u:s)(?:)(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 352", pattern: br"^(?u:s)(?#c)(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 353", pattern: br"^(?u:s)(?u:(?u:s))$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 354", pattern: br"^(?u:(?u:s)s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 355", pattern: br"^(?u:s(?d:s))$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 356", pattern: br"^(?u:s(?d)s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 357", pattern: br"^(?d:s)(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 358", pattern: br"^(?i:(?u:s)(?u:s))$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 359", pattern: br"^\N{U+73.73}s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 360", pattern: br"^s\N{U+73.73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 3200]) },
        Shape { name: "charset 361", pattern: br"^\N{U+73}s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 362", pattern: br"^s\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 363", pattern: br"^\N{U+73}\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 364", pattern: br"^\N{U+73.73}s\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 365", pattern: br"^\N{U+73.73}(?:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 366", pattern: br"^(?:\N{U+73.73})s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 367", pattern: br"^\N{U+73.73}\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 368", pattern: br"^\N{U+DF}s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 369", pattern: br"^s\N{U+DF}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 370", pattern: br"^\N{U+DF}\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 371", pattern: br"^[\N{U+61.73}]s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "charset 372", pattern: br"^(?^d:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 373", pattern: br"^(?u-u:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 374", pattern: br"^(?a-u:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 375", pattern: br"^(?d-u:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 376", pattern: br"^(?^id:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 377", pattern: br"^(?^dd:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 378", pattern: br"^(?a-d:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 379", pattern: br"^(?i^:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 380", pattern: br"^(?^-i:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 381", pattern: br"^(?ei:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 382", pattern: br"^\N{U+73.73}\x{100}?s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 383", pattern: br"^(?u:\N{U+73.73})s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 384", pattern: br"^\x{73}\x{73}s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 385", pattern: br"^\N{U+73.73}(?u:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 386", pattern: br"^\N{U+73.73}[s]$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "charset 387", pattern: br"^\N{U+73.73}\N{U+73.73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[700, 900, 1000, 1100]) },
        Shape { name: "charset 388", pattern: br"^\N{U+73}\N{U+73}\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 389", pattern: br"^ss\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 390", pattern: br"^\N{U+73}ss$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 391", pattern: br"^s\N{U+73}s$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 392", pattern: br"^\N{U+DF.73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 393", pattern: br"^\N{U+73.DF}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 394", pattern: br"^\N{U+73}\x{DF}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 395", pattern: br"^\x{DF}\N{U+73}$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[2300, 2400, 3200]) },
        Shape { name: "charset 396", pattern: br"^(?i-m-s:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 397", pattern: br"^(?-m-s:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 398", pattern: br"^(?i-ms:s)$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "charset 399", pattern: br"^(?i--m:s)$", values: CHARSET_VALUES, verdict: Verdict::Refused },
        Shape { name: "charset 400", pattern: br"^s(?s:(?:s))\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 401", pattern: br"^s(?s:s)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 402", pattern: br"^(?s:(?:s))s\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 403", pattern: br"^s(?m:(?:(?:s)))\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 3100]) },
        Shape { name: "charset 404", pattern: br"^s(?s:(?:s|s))\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2300, 3200]) },
        Shape { name: "charset 405", pattern: br"^s(?s:(s))\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200]) },
        Shape { name: "charset 406", pattern: br"^s(?s:(?:s).)\x{100}?$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[600, 700, 900, 2400, 3700]) },
        Shape { name: "charset 407", pattern: br"^s(?s:(?:.))$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[800, 2200, 2300, 3200]) },
        Shape { name: "charset 408", pattern: br"^s(?s:(?:s)(?:s))$", values: CHARSET_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 1", pattern: br"^r\N{U+41}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100]) },
        Shape { name: "named 2", pattern: br"^r\N{U+}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 3", pattern: br"^r\N{U+zz}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 4", pattern: br"^r\N{U+4_1}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100]) },
        Shape { name: "named 5", pattern: br"^r\N{U+41.42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 6", pattern: br"^r\N{U+0041}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100]) },
        Shape { name: "named 7", pattern: br"^r\N{U+41_}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 8", pattern: br"^r\N{U+_41}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 9", pattern: br"^r\N{U+41.}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 10", pattern: br"^r\N{U+41..42}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 11", pattern: br"^r\N{U+41.4_2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 12", pattern: br"^r\N{U+.41}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 13", pattern: br"^r\N{U+4__1}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 14", pattern: br"^r[\N{U+41}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100]) },
        Shape { name: "named 15", pattern: br"^r[\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 16", pattern: br"^r\N{U+41.42}{2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 17", pattern: br"^r\N{U+41}{2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1300]) },
        Shape { name: "named 18", pattern: br"^r\N{U+0}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[400]) },
        Shape { name: "named 19", pattern: br"^r\N{U+00}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[400]) },
        Shape { name: "named 20", pattern: br"^r\N{U+41", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 21", pattern: br"^r\N{U+41}+$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100, 1300]) },
        Shape { name: "named 22", pattern: br"^r\N{U+4}1$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1400]) },
        Shape { name: "named 23", pattern: br"^r\N{U+41.42}+$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1500, 1600]) },
        Shape { name: "named 24", pattern: br"^r\N{U+-41}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 25", pattern: br"^r\N{U+x41}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 26", pattern: br"^r\N{U+41.42.}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 27", pattern: br"^r\N{U+41.42.41.42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 28", pattern: br"^r\N{U+41.42}?$", values: NAMED_VALUES, verdict: Verdict::Rows(&[500, 700, 1600]) },
        Shape { name: "named 29", pattern: br"^r\N{U+41.42}*$", values: NAMED_VALUES, verdict: Verdict::Rows(&[500, 700, 1500, 1600]) },
        Shape { name: "named 30", pattern: br"^r\N{U+41.42}{1,2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1500, 1600]) },
        Shape { name: "named 31", pattern: br"^r(?:\N{U+41.42})$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 32", pattern: br"^r\N{U+41.42}|^r\N{U+42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1500, 1600, 2300]) },
        Shape { name: "named 33", pattern: br"^r\N{U+41.42}B$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 34", pattern: br"^r\N{U+41.100}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 35", pattern: br"^r\N{U+41.42.}+$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 36", pattern: br"^r\N{U+41.1F}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 37", pattern: br"^r\N{U+41.4__2}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 38", pattern: br"^r\N{U+41._42}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 39", pattern: br"^r\N{U+41.42_}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 40", pattern: br"^r\N{U+0041.0042}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 41", pattern: br"^r\N{U+00000041}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100]) },
        Shape { name: "named 42", pattern: br"^r\N{U+FFFFFFFFFF}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 43", pattern: br"^r\N{U+41.}+$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 44", pattern: br"^r\N{U+41.42}\N{U+41.42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 45", pattern: br"^r\N{U+41}\N{U+42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 46", pattern: br"^r(?-i)\N{U+41.42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700]) },
        Shape { name: "named 47", pattern: br"^r(?aa)\N{U+41.42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 48", pattern: br"^r\N{U+A}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 49", pattern: br"^r\N{U+41.42}{0}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[500]) },
        Shape { name: "named 50", pattern: br"^r\N{U+41.42}B+$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 51", pattern: br"^r\N{U+41.42}{2,}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 52", pattern: br"^r\N{U+41}\N{U+42}{2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 53", pattern: br"^r[\N{U+41}\N{U+42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100, 2300]) },
        Shape { name: "named 54", pattern: br"^r[\N{U+41}-\N{U+42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100, 2300]) },
        Shape { name: "named 55", pattern: br"^r[^\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[400, 600, 800, 900, 1000, 1100, 1200, 2300]) },
        Shape { name: "named 56", pattern: br"^r[x\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 57", pattern: br"^r[\N{U+41.42}]+$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1500, 1600]) },
        Shape { name: "named 58", pattern: br"^r[\N{U+41.42}B]{2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 59", pattern: br"^r[\N{U+41.42}-\N{U+43}]$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 60", pattern: br"^r[\N{U+41.42}\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 61", pattern: br"^r(?i)\N{U+61.62}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 62", pattern: br"^r\N{U+00.41}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "named 63", pattern: br"^r\N{U+41.0}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1800]) },
        Shape { name: "named 64", pattern: br"^r\N{U+61}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100]) },
        Shape { name: "named 65", pattern: br"^r\N{U+e9}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1000, 1200]) },
        Shape { name: "named 66", pattern: br"^r\N{U+41.42}B?$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 67", pattern: br"^r(\N{U+41.42})$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 68", pattern: br"^r\N{U+41.42}|B$", values: NAMED_VALUES, verdict: Verdict::Rows(&[200, 700, 1500, 1600, 1900, 2000, 2300]) },
        Shape { name: "named 69", pattern: br"^r\N{U+4_1.4_2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 70", pattern: br"^r\N{U+41.42.41}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 71", pattern: br"^r[^\N{U+41}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[400, 800, 900, 1000, 1200, 2300]) },
        Shape { name: "named 72", pattern: br"^r[\N{U+41.42}x]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 73", pattern: br"^r[^\N{U+41.42}x]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[400, 600, 800, 900, 1000, 1100, 1200, 2300]) },
        Shape { name: "named 74", pattern: br"^r[\N{U+41.42}]{2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 75", pattern: br"^r[\N{U+41.42}]?$", values: NAMED_VALUES, verdict: Verdict::Rows(&[500, 700, 1600]) },
        Shape { name: "named 76", pattern: br"^r[\N{U+41.42}][\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 77", pattern: br"^r[\N{U+41.42}\N{U+42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600, 2300]) },
        Shape { name: "named 78", pattern: br"^r[\N{U+42}\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600, 2300]) },
        Shape { name: "named 79", pattern: br"^r[\N{U+41.42}-B]$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 80", pattern: br"^r[A-\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100, 2300]) },
        Shape { name: "named 81", pattern: br"^r\N{U+041.042}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 82", pattern: br"^r\N{U+41.42}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 83", pattern: br"^r\N{U+0}A$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "named 84", pattern: br"^r\N{U+100000000}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 85", pattern: br"^r\N{U+7FFFFFFFFFFFFFFF}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 86", pattern: br"^r\N{U+FFFFFFFFFFFFFFFFFFFF}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 87", pattern: br"^r\N{U+D800}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 88", pattern: br"^r\N{U+110000}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 89", pattern: br"^r\N{U+41.D800}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 90", pattern: br"^r\N{U+41.42}{,2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[500, 700, 1500, 1600]) },
        Shape { name: "named 91", pattern: br"^r\N{U+41.42}{2}{2}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 92", pattern: br"^r\N{U+41.42}+?$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1500, 1600]) },
        Shape { name: "named 93", pattern: br"^r\N{U+41.42}++$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1500, 1600]) },
        Shape { name: "named 94", pattern: br"^r\N{U+41.42}(?#x){2}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[1500]) },
        Shape { name: "named 95", pattern: br"^r\N{U+41.42}\N{U+41}$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 96", pattern: br"^r[\N{U+41.42}-]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 97", pattern: br"^r[-\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 98", pattern: br"^r[\N{U+41.42}]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[700, 1600]) },
        Shape { name: "named 99", pattern: br"^r\N{U+41.42}B$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 100", pattern: br"^r[\N{U+41.42}-Z]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 101", pattern: br"^r[\N{U+61.62}-c]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[600, 1100, 2300]) },
        Shape { name: "named 102", pattern: br"^r[\N{U+41.42}-N]$", values: NAMED_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "named 103", pattern: br"^r[\N{U+41.42}-M]$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 104", pattern: br"^r\N{U+8_000_000_000_000_000}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 105", pattern: br"^r\N{U+08000000000000000}$", values: NAMED_VALUES, verdict: Verdict::Refused },
        Shape { name: "named 106", pattern: br"^r\N{U+FFFF_FFFF_FFFF_FFFF_FFFF}$", values: NAMED_VALUES, verdict: Verdict::Refused },
    ];

    /// The measured shapes written in the filter string, through its tokenizer.
    #[rustfmt::skip]
    const FILTER_SHAPES: &[Shape] = &[
        Shape { name: "b-ss", pattern: br"^ss$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "b-xdf", pattern: br"^\xDF$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 2600, 2900, 3200, 3400]) },
        Shape { name: "b-bxdf", pattern: br"^\x{DF}$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 2600, 2900, 3200, 3400]) },
        Shape { name: "b-utf8ss", pattern: b"^\xC3\x9F$", values: FOLD_VALUES, verdict: Verdict::Rows(&[400, 3500]) },
        Shape { name: "b-fi", pattern: br"^fi$", values: FOLD_VALUES, verdict: Verdict::Rows(&[700, 800, 2400]) },
        Shape { name: "b-s", pattern: br"^s$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]) },
        Shape { name: "b-k", pattern: br"^k$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1300, 1400]) },
        Shape { name: "b-eacute", pattern: br"^\xE9$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600]) },
        Shape { name: "b-rawdf", pattern: b"^\xDF$", values: FOLD_VALUES, verdict: Verdict::Rows(&[500, 2600, 2900, 3200, 3400]) },
        Shape { name: "u-ss", pattern: br"^ss\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-xdf", pattern: br"^\xDF\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-bxdf", pattern: br"^\x{DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-rawdf", pattern: b"^\xDF\\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-utf8ss", pattern: b"^\xC3\x9F\\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[400, 3500]) },
        Shape { name: "u-fb01", pattern: br"^\x{FB01}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[700, 800, 2400]) },
        Shape { name: "u-fi", pattern: br"^fi\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[700, 800, 2400]) },
        Shape { name: "u-17f", pattern: br"^\x{17F}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]) },
        Shape { name: "u-s", pattern: br"^s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]) },
        Shape { name: "u-212a", pattern: br"^\x{212A}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1300, 1400]) },
        Shape { name: "u-kelvin", pattern: br"^k\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1300, 1400]) },
        Shape { name: "u-1e9e", pattern: br"^\x{1E9E}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-st", pattern: br"^st\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1600]) },
        Shape { name: "u-fb06", pattern: br"^\x{FB06}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1600]) },
        Shape { name: "u-fb05", pattern: br"^\x{FB05}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1600]) },
        Shape { name: "u-fb00", pattern: br"^\x{FB00}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2200]) },
        Shape { name: "u-fb03", pattern: br"^\x{FB03}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2400]) },
        Shape { name: "u-fb02", pattern: br"^\x{FB02}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[4100, 4200]) },
        Shape { name: "u-fb04", pattern: br"^\x{FB04}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[4200]) },
        Shape { name: "u-splus", pattern: br"^s+\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]) },
        Shape { name: "u-s2", pattern: br"^s{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "u-clsdf", pattern: br"^[\x{DF}]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-nclsdf", pattern: br"^[^\x{DF}]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "u-clss", pattern: br"^[s]s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-assa", pattern: br"^assa\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1900, 2000]) },
        Shape { name: "u-adfa", pattern: br"^a\x{DF}a\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1900, 2000]) },
        Shape { name: "u-sdf", pattern: br"^s\x{DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 2900, 3000]) },
        Shape { name: "u-sss", pattern: br"^sss\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 2900, 3000]) },
        Shape { name: "u-dfs", pattern: br"^\x{DF}s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 2900, 3000]) },
        Shape { name: "u-ssss", pattern: br"^ssss\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "u-dfdf", pattern: br"^\x{DF}\x{DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3400]) },
        Shape { name: "u-dfplus", pattern: br"^\x{DF}+\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-clsdfx2", pattern: br"^[\x{DF}x]{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3000, 3100, 3200, 3300, 3400]) },
        Shape { name: "u-clsdfx", pattern: br"^[\x{DF}x]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3300, 3400]) },
        Shape { name: "u-xdfq", pattern: br"^x\x{DF}?\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3100, 3200, 3300]) },
        Shape { name: "u-assplus", pattern: br"^ass+\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1800]) },
        Shape { name: "u-as2", pattern: br"^as{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1800]) },
        Shape { name: "u-eacute", pattern: br"^\xE9\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[3600, 3700]) },
        Shape { name: "u-mu", pattern: br"^\x{3BC}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2100]) },
        Shape { name: "u-micro", pattern: br"^\xB5\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2100]) },
        Shape { name: "u-cap", pattern: br"^\x{39C}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2100]) },
        Shape { name: "u-letter-i", pattern: br"^i\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[700, 800, 2400, 3800, 3900]) },
        Shape { name: "u-130", pattern: br"^\x{130}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "u-131", pattern: br"^\x{131}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "u-17fs", pattern: br"^\x{17F}s\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-17f17f", pattern: br"^\x{17F}\x{17F}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-Sdf", pattern: br"^S\x{DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[2700, 2800, 2900, 3000]) },
        Shape { name: "u-dfdot", pattern: br"^\x{DF}.\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1900, 2000, 2700, 2800, 3000, 3400]) },
        Shape { name: "u-dotdf", pattern: br"^.\x{DF}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-sw", pattern: br"^s\w\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1600, 1800, 2000, 2800, 2900, 3000, 3100]) },
        Shape { name: "u-fl", pattern: br"^fl\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[4100, 4200]) },
        Shape { name: "u-ffl", pattern: br"^ffl\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[4200]) },
        Shape { name: "u-clsdfs", pattern: br"^[\x{DF}s]\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1000, 1100, 1800, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-clsdfsplus", pattern: br"^[\x{DF}s]+\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1000, 1100, 1800, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "u-clsssdf", pattern: br"^[s]{2}\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 1800, 2800, 3000, 3100]) },
        Shape { name: "u-wordss", pattern: br"^\w\w\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 700, 800, 1600, 1800, 1900, 2000, 2200, 2400, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 4100, 4200]) },
        Shape { name: "u-dfz", pattern: br"^\x{DF}z?\x{100}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "f-pss", pattern: br"^ss\p{L}?$", values: FOLD_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 500, 1800, 1900, 2000, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3400]) },
        Shape { name: "e-rbrace", pattern: br"^r{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "e-rbracex", pattern: br"^r{x}$", values: ESC_VALUES, verdict: Verdict::Rows(&[200]) },
        Shape { name: "e-rcomma3", pattern: br"^r{,3}$", values: ESC_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]) },
        Shape { name: "e-acomma3", pattern: br"^a{,3}$", values: ESC_VALUES, verdict: Verdict::Rows(&[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]) },
        Shape { name: "e-o101", pattern: br"^r\o{101}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-101", pattern: br"^r\101$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-0", pattern: br"^r\0$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "e-00", pattern: br"^r\00$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "e-000", pattern: br"^r\000$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "e-0001", pattern: br"^r\0001$", values: ESC_VALUES, verdict: Verdict::Rows(&[3000]) },
        Shape { name: "e-0101", pattern: br"^r\0101$", values: ESC_VALUES, verdict: Verdict::Rows(&[2200]) },
        Shape { name: "e-351", pattern: br"^r\351$", values: ESC_VALUES, verdict: Verdict::Rows(&[1600]) },
        Shape { name: "e-1", pattern: br"^r\1$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "e-brace", pattern: br"^{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100, 1800]) },
        Shape { name: "e-rclose", pattern: br"^r}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1900]) },
        Shape { name: "e-r1comma", pattern: br"^r{1,$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-rcomma", pattern: br"^r{,}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1400]) },
        Shape { name: "e-r12", pattern: br"^r{1,2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2000, 2100]) },
        Shape { name: "e-resc12", pattern: br"^r\{1,2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-r2", pattern: br"^r{2}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2100]) },
        Shape { name: "e-rbx", pattern: br"^r{x$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-r21", pattern: br"^r{2,1}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-oblank", pattern: br"^r\o{ 101 }$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "e-ounder", pattern: br"^r\o{1_01}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-oempty", pattern: br"^r\o{}$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "e-o8", pattern: br"^r\o{8}$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "e-o101x", pattern: br"^r\o{101x}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-o400", pattern: br"^r\o{400}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-400", pattern: br"^r\400$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-o0", pattern: br"^r\o{0}$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "e-8", pattern: br"^r\8$", values: ESC_VALUES, verdict: Verdict::Refused },
        Shape { name: "e-07", pattern: br"^r\07$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-010", pattern: br"^r\010$", values: ESC_VALUES, verdict: Verdict::Rows(&[2300]) },
        Shape { name: "e-10", pattern: br"^r\10$", values: ESC_VALUES, verdict: Verdict::Rows(&[2300]) },
        Shape { name: "e-101z", pattern: br"^r\1010$", values: ESC_VALUES, verdict: Verdict::Rows(&[2900]) },
        Shape { name: "e-0x", pattern: br"^r\0A$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-2000", pattern: br"^r\2000$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-rbrace1", pattern: br"^r{1}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1200, 2000, 2100]) },
        Shape { name: "e-bracex", pattern: br"^{x}$", values: ESC_VALUES, verdict: Verdict::Rows(&[200, 2600]) },
        Shape { name: "e-rbracecls", pattern: br"^r[{]$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "e-brace3", pattern: br"^r{3,2,1}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-plusbrace", pattern: br"^r+{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "e-dotbrace", pattern: br"^.{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "e-lbrace", pattern: br"^r{1}{$", values: ESC_VALUES, verdict: Verdict::Rows(&[100]) },
        Shape { name: "e-ocomma", pattern: br"^r\o{1,01}$", values: ESC_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "e-o101plus", pattern: br"^r\o{101}+$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-o101q", pattern: br"^r\o{101}0?$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400, 2900]) },
        Shape { name: "e-101q", pattern: br"^r\1010?$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400, 2900]) },
        Shape { name: "e-rbracec", pattern: br"^r{,x}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-rbracec2", pattern: br"^r{1,x}$", values: ESC_VALUES, verdict: Verdict::Rows(&[]) },
        Shape { name: "e-cls-o", pattern: br"^r[\o{101}]$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-cls-101", pattern: br"^r[\101]$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 2400]) },
        Shape { name: "e-cls-0", pattern: br"^r[\0]$", values: ESC_VALUES, verdict: Verdict::Rows(&[800]) },
        Shape { name: "e-cls-1", pattern: br"^r[\1]$", values: ESC_VALUES, verdict: Verdict::Rows(&[1700]) },
        Shape { name: "e-o101q2", pattern: br"^r\o{101}?$", values: ESC_VALUES, verdict: Verdict::Rows(&[1000, 1200, 2000, 2100, 2400]) },
        Shape { name: "e-o400u", pattern: br"^r\o{400}?\xe9$", values: ESC_VALUES, verdict: Verdict::Rows(&[1600]) },
        Shape { name: "e-400u", pattern: br"^r\400?\xe9$", values: ESC_VALUES, verdict: Verdict::Rows(&[1600]) },
    ];

    /// The dialect sentence behind a literal run Perl's optimizer merges across
    /// a charset scope boundary and folds whole.
    const MERGED_RUN: &str =
        "a literal run that crosses a charset scope boundary folds whole in Perl alone";
    /// The dialect sentence behind a trie word ending inside the fold of a
    /// value's ß.
    const TRIE_WORD: &str =
        "an alternation of literal words whose word ends inside the fold of a value's ß matches in Perl alone";
    /// The dialect sentence behind a `\N{U+}` sequence of several characters
    /// Perl keeps as a run of its own beside a literal.
    const SEQUENCE_RUN: &str =
        "a \\N{U+} escape of several characters and the literal beside it fold as two runs in Perl alone";
    /// The dialect sentence behind the `\N{U+}` sequence shapes vep-rs refuses.
    const SEQUENCE_SHAPE: &str = "a \\N{U+} escape of several characters in a negated class, before a `-` or as the end of a range is refused";
    /// The dialect sentence behind a possessive quantifier, refused in every
    /// place, the end of a pattern included.
    const POSSESSIVE: &str = "a possessive quantifier is refused";

    /// The measured shapes vep-rs reads otherwise, each with vep-rs's verdict
    /// and the dialect sentence of `docs/divergences/transcript-filter-fields.md`
    /// it rests on.
    const EXCEPTIONS: &[(&str, Verdict, &str)] = &[
        (
            "fold 9",
            Verdict::Refused,
            "the locale flag (?l) is refused",
        ),
        (
            "fold 49",
            Verdict::Rows(&[100, 200, 300, 500, 1000, 1100]),
            TRIE_WORD,
        ),
        (
            "fold 50",
            Verdict::Rows(&[100, 200, 300, 500, 1000, 1100]),
            TRIE_WORD,
        ),
        ("brace 2", Verdict::Refused, "a back-reference is refused"),
        ("brace 39", Verdict::Refused, "a back-reference is refused"),
        ("brace 47", Verdict::Refused, POSSESSIVE),
        ("charset 42", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 244", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 43", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 249", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 250", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 255", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 256", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 257", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 259", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 305", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 306", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 315", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 316", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 317", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 318", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 319", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 321", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 322", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 323", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 327", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 328", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 329", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 350", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 355", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 356", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 357", Verdict::Rows(&[800, 2200]), MERGED_RUN),
        ("charset 286", Verdict::Rows(&[3600]), MERGED_RUN),
        (
            "charset 359",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 360",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 365",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 366",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 367",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 383",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 385",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        (
            "charset 386",
            Verdict::Rows(&[2300, 2400, 3200]),
            SEQUENCE_RUN,
        ),
        ("charset 45", Verdict::Rows(&[800, 2200, 3100]), TRIE_WORD),
        ("charset 312", Verdict::Rows(&[800, 2200]), TRIE_WORD),
        ("charset 313", Verdict::Rows(&[800, 2200]), TRIE_WORD),
        ("charset 404", Verdict::Rows(&[800, 2200]), TRIE_WORD),
        ("charset 371", Verdict::Rows(&[]), TRIE_WORD),
        ("named 55", Verdict::Refused, SEQUENCE_SHAPE),
        ("named 73", Verdict::Refused, SEQUENCE_SHAPE),
        ("named 80", Verdict::Refused, SEQUENCE_SHAPE),
        ("named 93", Verdict::Refused, POSSESSIVE),
        ("named 96", Verdict::Refused, SEQUENCE_SHAPE),
        ("named 100", Verdict::Refused, SEQUENCE_SHAPE),
        ("named 101", Verdict::Refused, SEQUENCE_SHAPE),
        ("named 102", Verdict::Refused, SEQUENCE_SHAPE),
    ];

    fn rows(matcher: &Matcher, values: &[(u32, &[u8])]) -> Vec<u32> {
        values
            .iter()
            .filter(|(_, id)| matcher.is_match(&filter_input::decode_line(id)))
            .map(|(pos, _)| *pos)
            .collect()
    }

    fn expected(shape: &Shape) -> Verdict {
        EXCEPTIONS
            .iter()
            .find(|(name, _, _)| *name == shape.name)
            .map_or(shape.verdict, |(_, verdict, _)| *verdict)
    }

    fn check(shapes: &[Shape], verdict: impl Fn(&Shape) -> Verdict) {
        let mismatches: Vec<String> = shapes
            .iter()
            .filter_map(|shape| {
                let ours = verdict(shape);
                let theirs = expected(shape);
                (ours != theirs).then(|| {
                    format!(
                        "{} {:?}: vep-rs {:?}, filter_vep {:?}",
                        shape.name,
                        String::from_utf8_lossy(shape.pattern),
                        ours,
                        theirs
                    )
                })
            })
            .collect();
        assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
    }

    /// Every measured shape carried through a field value compiles and matches
    /// as `filter_vep` does on the same bytes, the listed exceptions aside.
    #[test]
    fn field_patterns_match_as_filter_vep_does() {
        check(FIELD_SHAPES, |shape| {
            match compile(&filter_input::decode_line(shape.pattern)) {
                Ok(matcher) => {
                    Verdict::Rows(Box::leak(rows(&matcher, shape.values).into_boxed_slice()))
                }
                Err(_) => Verdict::Refused,
            }
        });
    }

    /// Every measured shape written in the filter string, through its
    /// tokenizer, matches as `filter_vep` does on the same bytes.
    #[test]
    fn filter_patterns_match_as_filter_vep_does() {
        struct Owned(String);
        impl FilterData for Owned {
            fn value(&self, key: &str) -> Option<Cow<'_, str>> {
                (key == "ID").then_some(Cow::Borrowed(self.0.as_str()))
            }
        }
        check(FILTER_SHAPES, |shape| {
            let filter = format!("ID match {}", filter_input::decode_line(shape.pattern));
            match FilterSet::new(&[filter], &["ID"]) {
                Ok(set) => {
                    let mut kept = Vec::new();
                    for (pos, id) in shape.values {
                        match set.try_evaluate(&Owned(filter_input::decode_line(id).into_owned())) {
                            Ok(true) => kept.push(*pos),
                            Ok(false) => {}
                            Err(_) => return Verdict::Refused,
                        }
                    }
                    Verdict::Rows(Box::leak(kept.into_boxed_slice()))
                }
                Err(_) => Verdict::Refused,
            }
        });
    }
}
