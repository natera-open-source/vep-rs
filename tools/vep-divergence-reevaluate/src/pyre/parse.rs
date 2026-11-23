// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pattern parser: CPython's `re` syntax for a `str` pattern, read into a tree, with the
//! message and the position `re.PatternError` gives a pattern it refuses.

use crate::pyfmt::{is_py_alpha, py_repr_str};

pub const MAXREPEAT: u32 = u32::MAX;
const MAXGROUPS: u128 = 0x3fff_ffff;
/// The cap of a width, above every repeat count and every text length.
const MAXWIDTH: u128 = 1 << 64;

pub const FLAG_IGNORECASE: u32 = 2;
pub const FLAG_LOCALE: u32 = 4;
pub const FLAG_MULTILINE: u32 = 8;
pub const FLAG_DOTALL: u32 = 16;
pub const FLAG_UNICODE: u32 = 32;
pub const FLAG_VERBOSE: u32 = 64;
pub const FLAG_ASCII: u32 = 256;
const TYPE_FLAGS: u32 = FLAG_ASCII | FLAG_LOCALE | FLAG_UNICODE;

/// A pattern the parser or the compiler refuses, as `str(re.PatternError)` reads.
#[derive(Clone, Debug, PartialEq)]
pub struct PatternError(pub String);

impl std::fmt::Display for PatternError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Category {
    Digit,
    NotDigit,
    Space,
    NotSpace,
    Word,
    NotWord,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum At {
    Beginning,
    BeginningString,
    Boundary,
    NonBoundary,
    End,
    EndString,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SetItem {
    Literal(u32),
    Range(u32, u32),
    Category(Category),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RepeatKind {
    Greedy,
    Lazy,
    Possessive,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Node {
    Literal(u32),
    NotLiteral(u32),
    Any,
    In {
        negate: bool,
        items: Vec<SetItem>,
    },
    At(At),
    Branch(Vec<Vec<Node>>),
    Subpattern {
        group: Option<usize>,
        add_flags: u32,
        del_flags: u32,
        body: Vec<Node>,
    },
    Atomic(Vec<Node>),
    Repeat {
        kind: RepeatKind,
        min: u32,
        max: u32,
        body: Vec<Node>,
    },
    GroupRef(usize),
    GroupRefExists {
        group: usize,
        yes: Vec<Node>,
        no: Option<Vec<Node>>,
    },
    Assert {
        behind: bool,
        negate: bool,
        body: Vec<Node>,
    },
    Failure,
}

/// A parsed pattern: its tree, its flags, the number of groups (the whole match included),
/// the named groups and the widths of the groups.
#[derive(Debug)]
pub struct Parsed {
    pub nodes: Vec<Node>,
    pub flags: u32,
    pub groups: usize,
    pub names: Vec<(String, usize)>,
    pub widths: Vec<Option<(u128, u128)>>,
}

/// One token: a character, or a backslash and the character after it.
type Token = Option<String>;

struct Tokenizer {
    chars: Vec<char>,
    index: usize,
    next: Token,
}

fn token_len(token: &Token) -> usize {
    token.as_ref().map_or(0, |t| t.chars().count())
}

impl Tokenizer {
    fn new(pattern: &str) -> Result<Self, PatternError> {
        let mut tokenizer = Tokenizer {
            chars: pattern.chars().collect(),
            index: 0,
            next: None,
        };
        tokenizer.advance()?;
        Ok(tokenizer)
    }

    fn advance(&mut self) -> Result<(), PatternError> {
        let mut index = self.index;
        let Some(&first) = self.chars.get(index) else {
            self.next = None;
            return Ok(());
        };
        let mut token = String::from(first);
        if first == '\\' {
            index += 1;
            match self.chars.get(index) {
                Some(&second) => token.push(second),
                None => {
                    return Err(self.error_at(
                        "bad escape (end of pattern)",
                        self.chars.len().saturating_sub(1),
                    ))
                }
            }
        }
        self.index = index + 1;
        self.next = Some(token);
        Ok(())
    }

    fn next_is(&self, text: &str) -> bool {
        self.next.as_deref() == Some(text)
    }

    fn matches(&mut self, text: &str) -> Result<bool, PatternError> {
        if self.next_is(text) {
            self.advance()?;
            return Ok(true);
        }
        Ok(false)
    }

    fn get(&mut self) -> Result<Token, PatternError> {
        let this = self.next.take();
        self.advance()?;
        Ok(this)
    }

    /// Up to `n` next tokens that are one character of the set each.
    fn get_while(&mut self, n: usize, set: fn(char) -> bool) -> Result<String, PatternError> {
        let mut result = String::new();
        for _ in 0..n {
            match self.next_char() {
                Some(c) if set(c) => {
                    result.push(c);
                    self.advance()?;
                }
                _ => break,
            }
        }
        Ok(result)
    }

    /// The next token when it is one character.
    fn next_char(&self) -> Option<char> {
        let token = self.next.as_deref()?;
        let mut chars = token.chars();
        let first = chars.next()?;
        chars.next().is_none().then_some(first)
    }

    fn get_until(&mut self, terminator: &str, name: &str) -> Result<String, PatternError> {
        let mut result = String::new();
        loop {
            let c = self.next.take();
            self.advance()?;
            let Some(c) = c else {
                if result.is_empty() {
                    return Err(self.error(&format!("missing {name}"), 0));
                }
                return Err(self.error(
                    &format!("missing {terminator}, unterminated name"),
                    result.chars().count(),
                ));
            };
            if c == terminator {
                if result.is_empty() {
                    return Err(self.error(&format!("missing {name}"), 1));
                }
                break;
            }
            result.push_str(&c);
        }
        Ok(result)
    }

    fn tell(&self) -> usize {
        self.index - token_len(&self.next)
    }

    fn seek(&mut self, index: usize) -> Result<(), PatternError> {
        self.index = index;
        self.advance()
    }

    fn error(&self, msg: &str, offset: usize) -> PatternError {
        // A position before the pattern cannot arise from the offsets the parser passes.
        self.error_at(msg, self.tell().saturating_sub(offset))
    }

    fn error_at(&self, msg: &str, pos: usize) -> PatternError {
        let mut text = format!("{msg} at position {pos}");
        if self.chars.contains(&'\n') {
            let before = &self.chars[..pos.min(self.chars.len())];
            let line = before.iter().filter(|c| **c == '\n').count() + 1;
            let column = match before.iter().rposition(|c| *c == '\n') {
                Some(newline) => pos - newline,
                None => pos + 1,
            };
            text.push_str(&format!(" (line {line}, column {column})"));
        }
        PatternError(text)
    }

    fn check_group_name(&self, name: &str, offset: usize) -> Result<(), PatternError> {
        if !is_identifier(name) {
            let msg = format!("bad character in group name {}", py_repr_str(name));
            return Err(self.error(&msg, name.chars().count() + offset));
        }
        Ok(())
    }
}

/// `str.isidentifier()`.
fn is_identifier(name: &str) -> bool {
    let mut chars = name.chars();
    match chars.next() {
        Some(first) if first == '_' || unicode_ident::is_xid_start(first) => {}
        _ => return false,
    }
    chars.all(unicode_ident::is_xid_continue)
}

struct State {
    flags: u32,
    names: Vec<(String, usize)>,
    widths: Vec<Option<(u128, u128)>>,
    lookbehind_groups: Option<usize>,
    /// The numbered groups conditionals name before they open, with where each is named.
    groupref_pos: Vec<(u128, usize)>,
}

impl State {
    fn groups(&self) -> usize {
        self.widths.len()
    }

    fn group_named(&self, name: &str) -> Option<usize> {
        self.names
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, gid)| *gid)
    }

    /// Opens a group; the error is the message alone, which the caller places.
    fn open_group(&mut self, name: Option<&str>) -> Result<usize, String> {
        let gid = self.groups();
        self.widths.push(None);
        if let Some(name) = name {
            if let Some(earlier) = self.group_named(name) {
                return Err(format!(
                    "redefinition of group name {} as group {gid}; was group {earlier}",
                    py_repr_str(name)
                ));
            }
            self.names.push((name.to_string(), gid));
        }
        Ok(gid)
    }

    fn check_group(&self, gid: u128) -> bool {
        usize::try_from(gid)
            .ok()
            .and_then(|gid| self.widths.get(gid))
            .is_some_and(Option::is_some)
    }

    fn check_lookbehind_group(&self, gid: u128, source: &Tokenizer) -> Result<(), PatternError> {
        if let Some(first) = self.lookbehind_groups {
            if !self.check_group(gid) {
                return Err(source.error("cannot refer to an open group", 0));
            }
            if gid >= first as u128 {
                return Err(source.error(
                    "cannot refer to group defined in the same lookbehind subpattern",
                    0,
                ));
            }
        }
        Ok(())
    }
}

/// The least and the greatest number of characters a sequence takes.
pub fn width(nodes: &[Node], widths: &[Option<(u128, u128)>]) -> (u128, u128) {
    let mut lo: u128 = 0;
    let mut hi: u128 = 0;
    for node in nodes {
        let (i, j) = match node {
            Node::Branch(alternatives) => {
                let mut i = MAXWIDTH;
                let mut j = 0;
                for alternative in alternatives {
                    let (l, h) = width(alternative, widths);
                    i = i.min(l);
                    j = j.max(h);
                }
                (i, j)
            }
            Node::Atomic(body) | Node::Subpattern { body, .. } => width(body, widths),
            Node::Repeat { min, max, body, .. } => {
                let (i, j) = width(body, widths);
                let lo_part = i.saturating_mul(u128::from(*min));
                if *max == MAXREPEAT && j != 0 {
                    lo = lo.saturating_add(lo_part);
                    hi = MAXWIDTH;
                    continue;
                }
                (lo_part, j.saturating_mul(u128::from(*max)))
            }
            Node::Any | Node::In { .. } | Node::Literal(_) | Node::NotLiteral(_) => (1, 1),
            Node::GroupRef(group) => widths[*group].unwrap_or((0, 0)),
            Node::GroupRefExists { yes, no, .. } => {
                let (mut i, mut j) = width(yes, widths);
                match no {
                    Some(no) => {
                        let (l, h) = width(no, widths);
                        i = i.min(l);
                        j = j.max(h);
                    }
                    None => i = 0,
                }
                (i, j)
            }
            Node::At(_) | Node::Assert { .. } | Node::Failure => (0, 0),
        };
        lo = lo.saturating_add(i);
        hi = hi.saturating_add(j);
    }
    (lo.min(MAXWIDTH), hi.min(MAXWIDTH))
}

/// The items of a set, each once, in the order first written.
fn unique(set: Vec<SetItem>) -> Vec<SetItem> {
    let mut unique: Vec<SetItem> = Vec::with_capacity(set.len());
    for item in set {
        if !unique.contains(&item) {
            unique.push(item);
        }
    }
    unique
}

/// Whether two nodes of this kind are the same when their parts are: the kinds that hold
/// no pattern of their own.
fn compares_equal(node: &Node) -> bool {
    matches!(
        node,
        Node::Literal(_)
            | Node::NotLiteral(_)
            | Node::Any
            | Node::In { .. }
            | Node::At(_)
            | Node::GroupRef(_)
            | Node::Failure
    )
}

fn is_digit(c: char) -> bool {
    c.is_ascii_digit()
}

fn is_octal(c: char) -> bool {
    ('0'..='7').contains(&c)
}

fn is_hex(c: char) -> bool {
    c.is_ascii_hexdigit()
}

/// What an escape reads as inside a character set or outside one.
enum Escaped {
    Literal(u32),
    Category(Category),
    At(At),
    GroupRef(usize),
}

fn category_escape(escape: &str) -> Option<Escaped> {
    Some(match escape {
        "\\A" => Escaped::At(At::BeginningString),
        "\\b" => Escaped::At(At::Boundary),
        "\\B" => Escaped::At(At::NonBoundary),
        "\\d" => Escaped::Category(Category::Digit),
        "\\D" => Escaped::Category(Category::NotDigit),
        "\\s" => Escaped::Category(Category::Space),
        "\\S" => Escaped::Category(Category::NotSpace),
        "\\w" => Escaped::Category(Category::Word),
        "\\W" => Escaped::Category(Category::NotWord),
        "\\z" | "\\Z" => Escaped::At(At::EndString),
        _ => return None,
    })
}

fn literal_escape(escape: &str) -> Option<u32> {
    Some(match escape {
        "\\a" => 0x07,
        "\\b" => 0x08,
        "\\f" => 0x0c,
        "\\n" => 0x0a,
        "\\r" => 0x0d,
        "\\t" => 0x09,
        "\\v" => 0x0b,
        "\\\\" => 0x5c,
        _ => return None,
    })
}

struct Parser {
    source: Tokenizer,
    state: State,
}

impl Parser {
    /// The escapes a set and an expression read alike: `\x`, `\u`, `\U` and `\N`. None for
    /// any other escape; the inner error is a code point outside Unicode, which falls to the
    /// bad-escape message with the digits read so far.
    fn common_escape(&mut self, escape: &mut String) -> Result<Option<u32>, PatternError> {
        let kind = escape.chars().nth(1);
        let digits = match kind {
            Some('x') => 2,
            Some('u') => 4,
            Some('U') => 8,
            Some('N') => {
                if !self.source.matches("{")? {
                    return Err(self.source.error("missing {", 0));
                }
                let name = self.source.get_until("}", "character name")?;
                // The tool carries no table of character names.
                let msg = format!("undefined character name {}", py_repr_str(&name));
                return Err(self.source.error(&msg, name.chars().count() + 4));
            }
            _ => return Ok(None),
        };
        escape.push_str(&self.source.get_while(digits, is_hex)?);
        if escape.len() != digits + 2 {
            let msg = format!("incomplete escape {escape}");
            return Err(self.source.error(&msg, escape.len()));
        }
        let code = u32::from_str_radix(&escape[2..], 16).expect("hex digits");
        if code > 0x10ffff {
            return Err(self.bad_escape(escape));
        }
        Ok(Some(code))
    }

    fn bad_escape(&self, escape: &str) -> PatternError {
        self.source
            .error(&format!("bad escape {escape}"), escape.chars().count())
    }

    fn octal_range_error(&self, escape: &str) -> PatternError {
        let msg = format!("octal escape value {escape} outside of range 0-0o377");
        self.source.error(&msg, escape.chars().count())
    }

    /// An escape inside a character set.
    fn class_escape(&mut self, token: &str) -> Result<SetItem, PatternError> {
        if let Some(code) = literal_escape(token) {
            return Ok(SetItem::Literal(code));
        }
        if let Some(Escaped::Category(category)) = category_escape(token) {
            return Ok(SetItem::Category(category));
        }
        let mut escape = token.to_string();
        if let Some(code) = self.common_escape(&mut escape)? {
            return Ok(SetItem::Literal(code));
        }
        let c = escape.chars().nth(1).expect("an escape of two characters");
        if is_octal(c) {
            escape.push_str(&self.source.get_while(2, is_octal)?);
            let code = u32::from_str_radix(&escape[1..], 8).expect("octal digits");
            if code > 0o377 {
                return Err(self.octal_range_error(&escape));
            }
            return Ok(SetItem::Literal(code));
        }
        if is_digit(c) || c.is_ascii_alphabetic() {
            return Err(self.bad_escape(&escape));
        }
        Ok(SetItem::Literal(c as u32))
    }

    /// An escape outside a character set.
    fn escape(&mut self, token: &str) -> Result<Escaped, PatternError> {
        if let Some(escaped) = category_escape(token) {
            return Ok(escaped);
        }
        if let Some(code) = literal_escape(token) {
            return Ok(Escaped::Literal(code));
        }
        let mut escape = token.to_string();
        if let Some(code) = self.common_escape(&mut escape)? {
            return Ok(Escaped::Literal(code));
        }
        let c = escape.chars().nth(1).expect("an escape of two characters");
        if c == '0' {
            escape.push_str(&self.source.get_while(2, is_octal)?);
            let code = u32::from_str_radix(&escape[1..], 8).expect("octal digits");
            return Ok(Escaped::Literal(code));
        }
        if is_digit(c) {
            // An octal escape of three digits, or a group reference of one or two.
            if let Some(second) = self.source.next_char().filter(|c| is_digit(*c)) {
                self.source.advance()?;
                escape.push(second);
                if is_octal(c) && is_octal(second) {
                    if let Some(third) = self.source.next_char().filter(|c| is_octal(*c)) {
                        self.source.advance()?;
                        escape.push(third);
                        let code = u32::from_str_radix(&escape[1..], 8).expect("octal digits");
                        if code > 0o377 {
                            return Err(self.octal_range_error(&escape));
                        }
                        return Ok(Escaped::Literal(code));
                    }
                }
            }
            let group: usize = escape[1..].parse().expect("one or two digits");
            if group < self.state.groups() {
                if !self.state.check_group(group as u128) {
                    return Err(self
                        .source
                        .error("cannot refer to an open group", escape.len()));
                }
                self.state
                    .check_lookbehind_group(group as u128, &self.source)?;
                return Ok(Escaped::GroupRef(group));
            }
            let msg = format!("invalid group reference {group}");
            return Err(self.source.error(&msg, escape.len() - 1));
        }
        if c.is_ascii_alphabetic() {
            return Err(self.bad_escape(&escape));
        }
        Ok(Escaped::Literal(c as u32))
    }

    /// An alternation, `a|b|c`.
    fn parse_sub(&mut self, verbose: bool, nested: usize) -> Result<Vec<Node>, PatternError> {
        let mut verbose = verbose;
        let mut items = Vec::new();
        loop {
            let first = nested == 0 && items.is_empty();
            items.push(self.parse(verbose, nested + 1, first)?);
            if !self.source.matches("|")? {
                break;
            }
            if nested == 0 {
                verbose = self.state.flags & FLAG_VERBOSE != 0;
            }
        }
        if items.len() == 1 {
            return Ok(items.pop().expect("one alternative"));
        }
        // What every alternative starts with moves in front of the alternation.
        let mut sub = Vec::new();
        loop {
            let shared = match items[0].first() {
                Some(first) if compares_equal(first) => {
                    items.iter().all(|item| item.first() == Some(first))
                }
                _ => false,
            };
            if !shared {
                break;
            }
            sub.push(items[0][0].clone());
            for item in &mut items {
                item.remove(0);
            }
        }
        // Alternatives of one character or one set each read as one set.
        let mut set: Vec<SetItem> = Vec::new();
        let all_sets = items.iter().all(|item| match item.as_slice() {
            [Node::Literal(code)] => {
                set.push(SetItem::Literal(*code));
                true
            }
            [Node::In {
                negate: false,
                items,
            }] => {
                set.extend(items.iter().cloned());
                true
            }
            _ => false,
        });
        if all_sets {
            sub.push(Node::In {
                negate: false,
                items: unique(set),
            });
        } else {
            sub.push(Node::Branch(items));
        }
        Ok(sub)
    }

    fn unterminated_group(&self, start: usize) -> PatternError {
        self.source.error(
            "missing ), unterminated subpattern",
            self.source.tell() - start,
        )
    }

    /// A sequence up to the next `|` or `)` of its own level.
    fn parse(
        &mut self,
        verbose: bool,
        nested: usize,
        first: bool,
    ) -> Result<Vec<Node>, PatternError> {
        let mut verbose = verbose;
        let mut sub: Vec<Node> = Vec::new();
        loop {
            let this = match self.source.next.as_deref() {
                None | Some("|") | Some(")") => break,
                Some(token) => token.to_string(),
            };
            self.source.advance()?;
            if verbose {
                if matches!(
                    this.as_str(),
                    " " | "\t" | "\n" | "\r" | "\u{0b}" | "\u{0c}"
                ) {
                    continue;
                }
                if this == "#" {
                    loop {
                        match self.source.get()?.as_deref() {
                            None | Some("\n") => break,
                            Some(_) => {}
                        }
                    }
                    continue;
                }
            }
            if this.starts_with('\\') {
                sub.push(match self.escape(&this)? {
                    Escaped::Literal(code) => Node::Literal(code),
                    Escaped::Category(category) => Node::In {
                        negate: false,
                        items: vec![SetItem::Category(category)],
                    },
                    Escaped::At(at) => Node::At(at),
                    Escaped::GroupRef(group) => Node::GroupRef(group),
                });
                continue;
            }
            let c = this.chars().next().expect("a token of one character");
            match c {
                '[' => sub.push(self.parse_set()?),
                '*' | '+' | '?' | '{' => {
                    let here = self.source.tell();
                    let (min, max) = match c {
                        '?' => (0, 1),
                        '*' => (0, MAXREPEAT),
                        '+' => (1, MAXREPEAT),
                        _ => {
                            if self.source.next_is("}") {
                                sub.push(Node::Literal('{' as u32));
                                continue;
                            }
                            let mut lo = String::new();
                            let mut hi = String::new();
                            while let Some(digit) = self.source.next_char().filter(|c| is_digit(*c))
                            {
                                lo.push(digit);
                                self.source.advance()?;
                            }
                            if self.source.matches(",")? {
                                while let Some(digit) =
                                    self.source.next_char().filter(|c| is_digit(*c))
                                {
                                    hi.push(digit);
                                    self.source.advance()?;
                                }
                            } else {
                                hi = lo.clone();
                            }
                            if !self.source.matches("}")? {
                                sub.push(Node::Literal('{' as u32));
                                self.source.seek(here)?;
                                continue;
                            }
                            let mut min = 0;
                            let mut max = MAXREPEAT;
                            if !lo.is_empty() {
                                min = repeat_count(&lo)?;
                            }
                            if !hi.is_empty() {
                                max = repeat_count(&hi)?;
                                if max < min {
                                    return Err(self.source.error(
                                        "min repeat greater than max repeat",
                                        self.source.tell() - here,
                                    ));
                                }
                            }
                            (min, max)
                        }
                    };
                    let offset = self.source.tell() - here + 1;
                    let body = match sub.last() {
                        None | Some(Node::At(_)) => {
                            return Err(self.source.error("nothing to repeat", offset))
                        }
                        Some(Node::Repeat { .. }) => {
                            return Err(self.source.error("multiple repeat", offset))
                        }
                        Some(Node::Subpattern {
                            group: None,
                            add_flags: 0,
                            del_flags: 0,
                            body,
                        }) => body.clone(),
                        Some(item) => vec![item.clone()],
                    };
                    let kind = if self.source.matches("?")? {
                        RepeatKind::Lazy
                    } else if self.source.matches("+")? {
                        RepeatKind::Possessive
                    } else {
                        RepeatKind::Greedy
                    };
                    *sub.last_mut().expect("an item to repeat") = Node::Repeat {
                        kind,
                        min,
                        max,
                        body,
                    };
                }
                '.' => sub.push(Node::Any),
                '(' => {
                    if let Some(node) = self.parse_group(&mut verbose, nested, first, &sub)? {
                        sub.push(node);
                    }
                }
                '^' => sub.push(Node::At(At::Beginning)),
                '$' => sub.push(Node::At(At::End)),
                _ => sub.push(Node::Literal(c as u32)),
            }
        }
        // A group that neither captures nor sets a flag reads as its content.
        let mut unpacked = Vec::with_capacity(sub.len());
        for node in sub {
            match node {
                Node::Subpattern {
                    group: None,
                    add_flags: 0,
                    del_flags: 0,
                    body,
                } => unpacked.extend(body),
                other => unpacked.push(other),
            }
        }
        Ok(unpacked)
    }

    /// A character set, the `[` already read.
    fn parse_set(&mut self) -> Result<Node, PatternError> {
        let here = self.source.tell() - 1;
        let mut set: Vec<SetItem> = Vec::new();
        let negate = self.source.matches("^")?;
        let unterminated = |parser: &Parser| -> PatternError {
            parser
                .source
                .error("unterminated character set", parser.source.tell() - here)
        };
        loop {
            let Some(this) = self.source.get()? else {
                return Err(unterminated(self));
            };
            if this == "]" && !set.is_empty() {
                break;
            }
            let code1 = if this.starts_with('\\') {
                self.class_escape(&this)?
            } else {
                SetItem::Literal(this.chars().next().expect("a token") as u32)
            };
            if !self.source.matches("-")? {
                set.push(code1);
                continue;
            }
            // A possible range.
            let Some(that) = self.source.get()? else {
                return Err(unterminated(self));
            };
            if that == "]" {
                set.push(code1);
                set.push(SetItem::Literal('-' as u32));
                break;
            }
            let code2 = if that.starts_with('\\') {
                self.class_escape(&that)?
            } else {
                SetItem::Literal(that.chars().next().expect("a token") as u32)
            };
            let bad_range = |parser: &Parser| -> PatternError {
                let msg = format!("bad character range {this}-{that}");
                parser
                    .source
                    .error(&msg, this.chars().count() + 1 + that.chars().count())
            };
            match (&code1, &code2) {
                (SetItem::Literal(lo), SetItem::Literal(hi)) if hi >= lo => {
                    set.push(SetItem::Range(*lo, *hi));
                }
                _ => return Err(bad_range(self)),
            }
        }
        let unique = unique(set);
        Ok(match unique.as_slice() {
            [SetItem::Literal(code)] if negate => Node::NotLiteral(*code),
            [SetItem::Literal(code)] => Node::Literal(*code),
            _ => Node::In {
                negate,
                items: unique,
            },
        })
    }

    /// A group, the `(` already read. None for a construct that adds no node of its own: a
    /// comment or global flags.
    fn parse_group(
        &mut self,
        verbose: &mut bool,
        nested: usize,
        first: bool,
        before: &[Node],
    ) -> Result<Option<Node>, PatternError> {
        let start = self.source.tell() - 1;
        let mut capture = true;
        let mut atomic = false;
        let mut name: Option<String> = None;
        let mut add_flags = 0;
        let mut del_flags = 0;
        if self.source.matches("?")? {
            let Some(char) = self.source.get()? else {
                return Err(self.source.error("unexpected end of pattern", 0));
            };
            match char.as_str() {
                "P" => {
                    if self.source.matches("<")? {
                        let group_name = self.source.get_until(">", "group name")?;
                        self.source.check_group_name(&group_name, 1)?;
                        name = Some(group_name);
                    } else if self.source.matches("=")? {
                        let group_name = self.source.get_until(")", "group name")?;
                        self.source.check_group_name(&group_name, 1)?;
                        let length = group_name.chars().count();
                        let Some(gid) = self.state.group_named(&group_name) else {
                            let msg = format!("unknown group name {}", py_repr_str(&group_name));
                            return Err(self.source.error(&msg, length + 1));
                        };
                        if !self.state.check_group(gid as u128) {
                            return Err(self
                                .source
                                .error("cannot refer to an open group", length + 1));
                        }
                        self.state
                            .check_lookbehind_group(gid as u128, &self.source)?;
                        return Ok(Some(Node::GroupRef(gid)));
                    } else {
                        let Some(char) = self.source.get()? else {
                            return Err(self.source.error("unexpected end of pattern", 0));
                        };
                        let msg = format!("unknown extension ?P{char}");
                        return Err(self.source.error(&msg, char.chars().count() + 2));
                    }
                }
                ":" => capture = false,
                "#" => {
                    loop {
                        if self.source.next.is_none() {
                            return Err(self.source.error(
                                "missing ), unterminated comment",
                                self.source.tell() - start,
                            ));
                        }
                        if self.source.get()?.as_deref() == Some(")") {
                            break;
                        }
                    }
                    return Ok(None);
                }
                "=" | "!" | "<" => {
                    let mut char = char;
                    let mut behind = false;
                    let mut outer_lookbehind = None;
                    if char == "<" {
                        char = match self.source.get()? {
                            Some(char) => char,
                            None => return Err(self.source.error("unexpected end of pattern", 0)),
                        };
                        if char != "=" && char != "!" {
                            let msg = format!("unknown extension ?<{char}");
                            return Err(self.source.error(&msg, char.chars().count() + 2));
                        }
                        behind = true;
                        outer_lookbehind = self.state.lookbehind_groups;
                        if outer_lookbehind.is_none() {
                            self.state.lookbehind_groups = Some(self.state.groups());
                        }
                    }
                    let body = self.parse_sub(*verbose, nested + 1)?;
                    if behind && outer_lookbehind.is_none() {
                        self.state.lookbehind_groups = None;
                    }
                    if !self.source.matches(")")? {
                        return Err(self.unterminated_group(start));
                    }
                    return Ok(Some(if char == "=" {
                        Node::Assert {
                            behind,
                            negate: false,
                            body,
                        }
                    } else if !body.is_empty() {
                        Node::Assert {
                            behind,
                            negate: true,
                            body,
                        }
                    } else {
                        Node::Failure
                    }));
                }
                "(" => {
                    let cond_name = self.source.get_until(")", "group name")?;
                    let length = cond_name.chars().count();
                    let group: u128 = if cond_name.bytes().all(|b| b.is_ascii_digit()) {
                        let digits = cond_name.trim_start_matches('0');
                        if digits.is_empty() {
                            return Err(self.source.error("bad group number", length + 1));
                        }
                        let group: u128 = digits.parse().unwrap_or(u128::MAX);
                        if group >= MAXGROUPS {
                            let msg = format!("invalid group reference {digits}");
                            return Err(self.source.error(&msg, length + 1));
                        }
                        if !self.state.groupref_pos.iter().any(|(g, _)| *g == group) {
                            self.state
                                .groupref_pos
                                .push((group, self.source.tell() - length - 1));
                        }
                        group
                    } else {
                        self.source.check_group_name(&cond_name, 1)?;
                        match self.state.group_named(&cond_name) {
                            Some(gid) => gid as u128,
                            None => {
                                let msg = format!("unknown group name {}", py_repr_str(&cond_name));
                                return Err(self.source.error(&msg, length + 1));
                            }
                        }
                    };
                    self.state.check_lookbehind_group(group, &self.source)?;
                    let yes = self.parse(*verbose, nested + 1, false)?;
                    let no = if self.source.matches("|")? {
                        let no = self.parse(*verbose, nested + 1, false)?;
                        if self.source.next_is("|") {
                            return Err(self
                                .source
                                .error("conditional backref with more than two branches", 0));
                        }
                        Some(no)
                    } else {
                        None
                    };
                    if !self.source.matches(")")? {
                        return Err(self.unterminated_group(start));
                    }
                    return Ok(Some(Node::GroupRefExists {
                        group: group as usize,
                        yes,
                        no,
                    }));
                }
                ">" => {
                    capture = false;
                    atomic = true;
                }
                "i" | "L" | "m" | "s" | "x" | "a" | "u" | "-" => match self.parse_flags(&char)? {
                    None => {
                        if !first || !before.is_empty() {
                            return Err(self.source.error(
                                "global flags not at the start of the expression",
                                self.source.tell() - start,
                            ));
                        }
                        *verbose = self.state.flags & FLAG_VERBOSE != 0;
                        return Ok(None);
                    }
                    Some((add, del)) => {
                        add_flags = add;
                        del_flags = del;
                        capture = false;
                    }
                },
                _ => {
                    let msg = format!("unknown extension ?{char}");
                    return Err(self.source.error(&msg, char.chars().count() + 1));
                }
            }
        }
        let group = if capture {
            match self.state.open_group(name.as_deref()) {
                Ok(gid) => Some(gid),
                Err(msg) => {
                    let length = name.as_deref().map_or(0, |n| n.chars().count());
                    return Err(self.source.error(&msg, length + 1));
                }
            }
        } else {
            None
        };
        let sub_verbose =
            (*verbose || add_flags & FLAG_VERBOSE != 0) && del_flags & FLAG_VERBOSE == 0;
        let body = self.parse_sub(sub_verbose, nested + 1)?;
        if !self.source.matches(")")? {
            return Err(self.unterminated_group(start));
        }
        if let Some(gid) = group {
            self.state.widths[gid] = Some(width(&body, &self.state.widths));
        }
        Ok(Some(if atomic {
            Node::Atomic(body)
        } else {
            Node::Subpattern {
                group,
                add_flags,
                del_flags,
                body,
            }
        }))
    }

    /// Inline flags, the first character already read. None for global flags, which close
    /// on `)`; the flags a group turns on and off otherwise.
    fn parse_flags(&mut self, first: &str) -> Result<Option<(u32, u32)>, PatternError> {
        let flag_of = |token: &str| -> Option<u32> {
            Some(match token {
                "i" => FLAG_IGNORECASE,
                "L" => FLAG_LOCALE,
                "m" => FLAG_MULTILINE,
                "s" => FLAG_DOTALL,
                "x" => FLAG_VERBOSE,
                "a" => FLAG_ASCII,
                "u" => FLAG_UNICODE,
                _ => return None,
            })
        };
        let is_alpha = |token: &str| !token.is_empty() && token.chars().all(is_py_alpha);
        let mut add_flags = 0;
        let mut del_flags = 0;
        let mut char = first.to_string();
        if char != "-" {
            loop {
                let flag = flag_of(&char).expect("a flag letter");
                if char == "L" {
                    return Err(self.source.error(
                        "bad inline flags: cannot use 'L' flag with a str pattern",
                        0,
                    ));
                }
                add_flags |= flag;
                if flag & TYPE_FLAGS != 0 && add_flags & TYPE_FLAGS != flag {
                    return Err(self.source.error(
                        "bad inline flags: flags 'a', 'u' and 'L' are incompatible",
                        0,
                    ));
                }
                char = match self.source.get()? {
                    Some(char) => char,
                    None => return Err(self.source.error("missing -, : or )", 0)),
                };
                if matches!(char.as_str(), ")" | "-" | ":") {
                    break;
                }
                if flag_of(&char).is_none() {
                    let msg = if is_alpha(&char) {
                        "unknown flag"
                    } else {
                        "missing -, : or )"
                    };
                    return Err(self.source.error(msg, char.chars().count()));
                }
            }
        }
        if char == ")" {
            self.state.flags |= add_flags;
            return Ok(None);
        }
        if char == "-" {
            char = match self.source.get()? {
                Some(char) => char,
                None => return Err(self.source.error("missing flag", 0)),
            };
            if flag_of(&char).is_none() {
                let msg = if is_alpha(&char) {
                    "unknown flag"
                } else {
                    "missing flag"
                };
                return Err(self.source.error(msg, char.chars().count()));
            }
            loop {
                let flag = flag_of(&char).expect("a flag letter");
                if flag & TYPE_FLAGS != 0 {
                    return Err(self.source.error(
                        "bad inline flags: cannot turn off flags 'a', 'u' and 'L'",
                        0,
                    ));
                }
                del_flags |= flag;
                char = match self.source.get()? {
                    Some(char) => char,
                    None => return Err(self.source.error("missing :", 0)),
                };
                if char == ":" {
                    break;
                }
                if flag_of(&char).is_none() {
                    let msg = if is_alpha(&char) {
                        "unknown flag"
                    } else {
                        "missing :"
                    };
                    return Err(self.source.error(msg, char.chars().count()));
                }
            }
        }
        if add_flags & del_flags != 0 {
            return Err(self
                .source
                .error("bad inline flags: flag turned on and off", 1));
        }
        Ok(Some((add_flags, del_flags)))
    }
}

/// A repeat count below the engine's bound; a count at or above it is the `OverflowError`
/// CPython raises, stated in its words.
fn repeat_count(digits: &str) -> Result<u32, PatternError> {
    match digits.parse::<u32>() {
        Ok(count) if count < MAXREPEAT => Ok(count),
        _ => Err(PatternError(
            "the repetition number is too large".to_string(),
        )),
    }
}

/// Parses a pattern as `re.compile(pattern)` does, with no flag given.
pub fn parse(pattern: &str) -> Result<Parsed, PatternError> {
    let source = Tokenizer::new(pattern)?;
    let state = State {
        flags: 0,
        names: Vec::new(),
        widths: vec![None],
        lookbehind_groups: None,
        groupref_pos: Vec::new(),
    };
    let mut parser = Parser { source, state };
    let nodes = parser.parse_sub(false, 0)?;
    let mut flags = parser.state.flags;
    if flags & FLAG_ASCII == 0 {
        flags |= FLAG_UNICODE;
    } else if flags & FLAG_UNICODE != 0 {
        return Err(PatternError(
            "ASCII and UNICODE flags are incompatible".to_string(),
        ));
    }
    if parser.source.next.is_some() {
        return Err(parser.source.error("unbalanced parenthesis", 0));
    }
    let groups = parser.state.groups();
    for (group, pos) in &parser.state.groupref_pos {
        if *group >= groups as u128 {
            let msg = format!("invalid group reference {group}");
            return Err(parser.source.error_at(&msg, *pos));
        }
    }
    Ok(Parsed {
        nodes,
        flags,
        groups,
        names: parser.state.names,
        widths: parser.state.widths,
    })
}
