// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pattern compiler: a parsed tree into the instructions the matcher runs, with the
//! flags of each part applied as CPython's `re` applies them (case folding, the ASCII and
//! Unicode categories, multi-line anchors).

use super::parse::{
    width, At, Category, Node, Parsed, PatternError, RepeatKind, SetItem, FLAG_DOTALL,
    FLAG_IGNORECASE, FLAG_MULTILINE, FLAG_UNICODE, MAXREPEAT,
};
use crate::pyfmt::{is_py_alnum, is_py_decimal, is_py_space};

const TYPE_FLAGS: u32 = super::parse::FLAG_ASCII | super::parse::FLAG_LOCALE | FLAG_UNICODE;
/// The widest look-behind the engine's code word holds.
const MAXCODE: u128 = 0xffff_ffff;

/// How a character is folded before it is held to a literal or a set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Fold {
    None,
    Ascii,
    Unicode,
}

#[derive(Clone, Debug)]
pub enum Atom {
    Literal(u32),
    Range(u32, u32),
    /// A range a folded character or its uppercase may fall in.
    RangeFolded(u32, u32),
    Category(Category, bool),
}

#[derive(Clone, Debug)]
pub struct CharSet {
    negate: bool,
    fold: Fold,
    /// The folded code points below U+10000 of a set that ignores case, sorted.
    folded: Vec<u32>,
    atoms: Vec<Atom>,
}

/// One character's worth of pattern.
#[derive(Clone, Debug)]
pub enum Unit {
    Literal(u32),
    NotLiteral(u32),
    LiteralFolded(u32, Fold),
    NotLiteralFolded(u32, Fold),
    Any,
    AnyAll,
    Set(usize),
}

#[derive(Clone, Copy, Debug)]
pub enum AtKind {
    Beginning,
    BeginningLine,
    End,
    EndLine,
    EndString,
    Boundary(bool),
    NonBoundary(bool),
}

#[derive(Clone, Debug)]
pub enum Inst {
    One(Unit),
    At(AtKind),
    /// Try the first target, then the second.
    Split(usize, usize),
    Jump(usize),
    Save(usize),
    RepeatStart(usize),
    /// A repeat of a body of any width, at the instruction the body returns to: `body` runs
    /// once more or `exit` follows.
    Repeat {
        register: usize,
        min: u32,
        max: u32,
        lazy: bool,
        body: usize,
        exit: usize,
    },
    /// A repeat of one character's worth of pattern.
    RepeatOne {
        unit: Unit,
        min: u32,
        max: u32,
        kind: RepeatKind,
    },
    PossessiveRepeat {
        min: u32,
        max: u32,
        body: usize,
        exit: usize,
    },
    Atomic {
        body: usize,
        exit: usize,
    },
    Look {
        behind: Option<usize>,
        negate: bool,
        body: usize,
        exit: usize,
    },
    /// The end of a body that runs as a match of its own.
    SubEnd,
    BackRef(usize, Fold),
    /// Falls through when the group took part in the match, else jumps.
    Conditional {
        group: usize,
        otherwise: usize,
    },
    Fail,
    Match,
}

#[derive(Debug)]
pub struct Program {
    pub insts: Vec<Inst>,
    pub sets: Vec<CharSet>,
    pub groups: usize,
    pub registers: usize,
    pub names: Vec<(String, usize)>,
    /// The set a search holds a match's first character to before it tries the pattern.
    pub first: Option<usize>,
}

pub fn ascii_lower(c: u32) -> u32 {
    if (0x41..=0x5a).contains(&c) {
        c + 0x20
    } else {
        c
    }
}

/// The first code point of the character's full lowercase mapping.
pub fn unicode_lower(c: u32) -> u32 {
    char::from_u32(c)
        .and_then(|ch| ch.to_lowercase().next())
        .map_or(c, |lower| lower as u32)
}

pub fn unicode_upper(c: u32) -> u32 {
    char::from_u32(c)
        .and_then(|ch| ch.to_uppercase().next())
        .map_or(c, |upper| upper as u32)
}

fn ascii_is_cased(c: u32) -> bool {
    c < 128 && (c as u8).is_ascii_alphabetic()
}

fn unicode_is_cased(c: u32) -> bool {
    c != unicode_lower(c) || c != unicode_upper(c)
}

/// The lowercase letters that share an uppercase with another lowercase letter, each with
/// the others: the equivalences a case-insensitive match adds to the lowercase mapping.
fn extra_cases(c: u32) -> &'static [u32] {
    match c {
        0x0069 => &[0x0131],
        0x0073 => &[0x017f],
        0x00b5 => &[0x03bc],
        0x0131 => &[0x0069],
        0x017f => &[0x0073],
        0x0345 => &[0x03b9, 0x1fbe],
        0x0390 => &[0x1fd3],
        0x03b0 => &[0x1fe3],
        0x03b2 => &[0x03d0],
        0x03b5 => &[0x03f5],
        0x03b8 => &[0x03d1],
        0x03b9 => &[0x0345, 0x1fbe],
        0x03ba => &[0x03f0],
        0x03bc => &[0x00b5],
        0x03c0 => &[0x03d6],
        0x03c1 => &[0x03f1],
        0x03c2 => &[0x03c3],
        0x03c3 => &[0x03c2],
        0x03c6 => &[0x03d5],
        0x03d0 => &[0x03b2],
        0x03d1 => &[0x03b8],
        0x03d5 => &[0x03c6],
        0x03d6 => &[0x03c0],
        0x03f0 => &[0x03ba],
        0x03f1 => &[0x03c1],
        0x03f5 => &[0x03b5],
        0x0432 => &[0x1c80],
        0x0434 => &[0x1c81],
        0x043e => &[0x1c82],
        0x0441 => &[0x1c83],
        0x0442 => &[0x1c84, 0x1c85],
        0x044a => &[0x1c86],
        0x0463 => &[0x1c87],
        0x1c80 => &[0x0432],
        0x1c81 => &[0x0434],
        0x1c82 => &[0x043e],
        0x1c83 => &[0x0441],
        0x1c84 => &[0x0442, 0x1c85],
        0x1c85 => &[0x0442, 0x1c84],
        0x1c86 => &[0x044a],
        0x1c87 => &[0x0463],
        0x1c88 => &[0xa64b],
        0x1e61 => &[0x1e9b],
        0x1e9b => &[0x1e61],
        0x1fbe => &[0x0345, 0x03b9],
        0x1fd3 => &[0x0390],
        0x1fe3 => &[0x03b0],
        0xa64b => &[0x1c88],
        0xfb05 => &[0xfb06],
        0xfb06 => &[0xfb05],
        _ => &[],
    }
}

/// Whether a code point is of a category, by the ASCII rule or the Unicode one.
pub fn in_category(category: Category, unicode: bool, c: u32) -> bool {
    let ascii = c < 128;
    let byte = c as u8;
    let ch = char::from_u32(c);
    let digit = if unicode {
        ch.is_some_and(is_py_decimal)
    } else {
        ascii && byte.is_ascii_digit()
    };
    let space = if unicode {
        ch.is_some_and(is_py_space)
    } else {
        ascii && matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
    };
    let word = if unicode {
        ch.is_some_and(|ch| ch == '_' || is_py_alnum(ch))
    } else {
        ascii && (byte.is_ascii_alphanumeric() || byte == b'_')
    };
    match category {
        Category::Digit => digit,
        Category::NotDigit => !digit,
        Category::Space => space,
        Category::NotSpace => !space,
        Category::Word => word,
        Category::NotWord => !word,
    }
}

impl CharSet {
    pub fn contains(&self, c: u32) -> bool {
        let folded = match self.fold {
            Fold::None => c,
            Fold::Ascii => ascii_lower(c),
            Fold::Unicode => unicode_lower(c),
        };
        let hit = self.folded.binary_search(&folded).is_ok()
            || self.atoms.iter().any(|atom| match atom {
                Atom::Literal(code) => folded == *code,
                Atom::Range(lo, hi) => (*lo..=*hi).contains(&folded),
                Atom::RangeFolded(lo, hi) => {
                    (*lo..=*hi).contains(&folded) || (*lo..=*hi).contains(&unicode_upper(folded))
                }
                Atom::Category(category, unicode) => in_category(*category, *unicode, folded),
            });
        hit != self.negate
    }
}

pub fn unit_matches(unit: &Unit, sets: &[CharSet], c: u32) -> bool {
    let fold = |fold: &Fold| match fold {
        Fold::None => c,
        Fold::Ascii => ascii_lower(c),
        Fold::Unicode => unicode_lower(c),
    };
    match unit {
        Unit::Literal(code) => c == *code,
        Unit::NotLiteral(code) => c != *code,
        Unit::LiteralFolded(code, how) => fold(how) == *code,
        Unit::NotLiteralFolded(code, how) => fold(how) != *code,
        Unit::Any => c != '\n' as u32,
        Unit::AnyAll => true,
        Unit::Set(index) => sets[*index].contains(c),
    }
}

struct Compiler<'a> {
    insts: Vec<Inst>,
    sets: Vec<CharSet>,
    registers: usize,
    widths: &'a [Option<(u128, u128)>],
}

fn combine_flags(flags: u32, add_flags: u32, del_flags: u32) -> u32 {
    let flags = if add_flags & TYPE_FLAGS != 0 {
        flags & !TYPE_FLAGS
    } else {
        flags
    };
    (flags | add_flags) & !del_flags
}

fn fold_of(flags: u32) -> Fold {
    if flags & FLAG_IGNORECASE == 0 {
        Fold::None
    } else if flags & FLAG_UNICODE != 0 {
        Fold::Unicode
    } else {
        Fold::Ascii
    }
}

impl Compiler<'_> {
    fn add_set(&mut self, set: CharSet) -> Unit {
        self.sets.push(set);
        Unit::Set(self.sets.len() - 1)
    }

    fn literal_unit(&mut self, code: u32, negate: bool, flags: u32) -> Unit {
        let plain = |code| {
            if negate {
                Unit::NotLiteral(code)
            } else {
                Unit::Literal(code)
            }
        };
        let fold = fold_of(flags);
        let cased = match fold {
            Fold::None => false,
            Fold::Ascii => ascii_is_cased(code),
            Fold::Unicode => unicode_is_cased(code),
        };
        if !cased {
            return plain(code);
        }
        let lower = match fold {
            Fold::Unicode => unicode_lower(code),
            _ => ascii_lower(code),
        };
        let extras = if fold == Fold::Unicode {
            extra_cases(lower)
        } else {
            &[]
        };
        if extras.is_empty() {
            return if negate {
                Unit::NotLiteralFolded(lower, fold)
            } else {
                Unit::LiteralFolded(lower, fold)
            };
        }
        let mut atoms = vec![Atom::Literal(lower)];
        atoms.extend(extras.iter().map(|extra| Atom::Literal(*extra)));
        self.add_set(CharSet {
            negate,
            fold,
            folded: Vec::new(),
            atoms,
        })
    }

    fn set_unit(&mut self, negate: bool, items: &[SetItem], flags: u32) -> Unit {
        let unicode = flags & FLAG_UNICODE != 0;
        let fold = fold_of(flags);
        if fold == Fold::None {
            let atoms = items
                .iter()
                .map(|item| match item {
                    SetItem::Literal(code) => Atom::Literal(*code),
                    SetItem::Range(lo, hi) => Atom::Range(*lo, *hi),
                    SetItem::Category(category) => Atom::Category(*category, unicode),
                })
                .collect();
            return self.add_set(CharSet {
                negate,
                fold,
                folded: Vec::new(),
                atoms,
            });
        }
        let lower = |c: u32| match fold {
            Fold::Unicode => unicode_lower(c),
            _ => ascii_lower(c),
        };
        let is_cased = |c: u32| match fold {
            Fold::Unicode => unicode_is_cased(c),
            _ => ascii_is_cased(c),
        };
        let mut folded: Vec<u32> = Vec::new();
        let mut atoms = Vec::new();
        let mut has_cased = false;
        let add_folded = |folded: &mut Vec<u32>, code: u32| {
            folded.push(code);
            if fold == Fold::Unicode {
                folded.extend_from_slice(extra_cases(code));
            }
        };
        for item in items {
            match item {
                SetItem::Literal(code) => {
                    let code = lower(*code);
                    if code < 0x10000 {
                        add_folded(&mut folded, code);
                    } else {
                        atoms.push(Atom::Literal(code));
                    }
                    has_cased = has_cased || is_cased(code);
                }
                SetItem::Range(lo, hi) => {
                    let mut beyond = false;
                    for code in *lo..=*hi {
                        let code = lower(code);
                        if code >= 0x10000 {
                            beyond = true;
                            break;
                        }
                        add_folded(&mut folded, code);
                    }
                    if beyond {
                        atoms.push(if fold == Fold::Unicode {
                            Atom::RangeFolded(*lo, *hi)
                        } else {
                            Atom::Range(*lo, *hi)
                        });
                        has_cased = true;
                    } else if !has_cased {
                        has_cased = (*lo..=*hi).any(is_cased);
                    }
                }
                SetItem::Category(category) => atoms.push(Atom::Category(*category, unicode)),
            }
        }
        folded.sort_unstable();
        folded.dedup();
        self.add_set(CharSet {
            negate,
            // A set with no cased member is held to the character as it is.
            fold: if has_cased { fold } else { Fold::None },
            folded,
            atoms,
        })
    }

    /// The one character's worth of pattern a body is, when it is that and nothing else.
    fn simple_unit(&mut self, body: &[Node], flags: u32) -> Option<Unit> {
        let [node] = body else {
            return None;
        };
        match node {
            Node::Subpattern {
                group: None,
                add_flags,
                del_flags,
                body,
            } => self.simple_unit(body, combine_flags(flags, *add_flags, *del_flags)),
            Node::Literal(code) => Some(self.literal_unit(*code, false, flags)),
            Node::NotLiteral(code) => Some(self.literal_unit(*code, true, flags)),
            Node::Any => Some(self.any_unit(flags)),
            Node::In { negate, items } => Some(self.set_unit(*negate, items, flags)),
            _ => None,
        }
    }

    fn any_unit(&self, flags: u32) -> Unit {
        if flags & FLAG_DOTALL != 0 {
            Unit::AnyAll
        } else {
            Unit::Any
        }
    }

    fn emit(&mut self, inst: Inst) -> usize {
        self.insts.push(inst);
        self.insts.len() - 1
    }

    fn sequence(&mut self, nodes: &[Node], flags: u32) -> Result<(), PatternError> {
        for node in nodes {
            self.node(node, flags)?;
        }
        Ok(())
    }

    fn node(&mut self, node: &Node, flags: u32) -> Result<(), PatternError> {
        match node {
            Node::Literal(code) => {
                let unit = self.literal_unit(*code, false, flags);
                self.emit(Inst::One(unit));
            }
            Node::NotLiteral(code) => {
                let unit = self.literal_unit(*code, true, flags);
                self.emit(Inst::One(unit));
            }
            Node::Any => {
                let unit = self.any_unit(flags);
                self.emit(Inst::One(unit));
            }
            Node::In { negate, items } => {
                let unit = self.set_unit(*negate, items, flags);
                self.emit(Inst::One(unit));
            }
            Node::At(at) => {
                let multiline = flags & FLAG_MULTILINE != 0;
                let unicode = flags & FLAG_UNICODE != 0;
                self.emit(Inst::At(match at {
                    At::Beginning if multiline => AtKind::BeginningLine,
                    At::Beginning | At::BeginningString => AtKind::Beginning,
                    At::End if multiline => AtKind::EndLine,
                    At::End => AtKind::End,
                    At::EndString => AtKind::EndString,
                    At::Boundary => AtKind::Boundary(unicode),
                    At::NonBoundary => AtKind::NonBoundary(unicode),
                }));
            }
            Node::Branch(alternatives) => {
                let mut jumps = Vec::new();
                for (index, alternative) in alternatives.iter().enumerate() {
                    let last = index + 1 == alternatives.len();
                    let split = (!last).then(|| self.emit(Inst::Split(0, 0)));
                    self.sequence(alternative, flags)?;
                    if !last {
                        jumps.push(self.emit(Inst::Jump(0)));
                    }
                    if let Some(split) = split {
                        self.insts[split] = Inst::Split(split + 1, self.insts.len());
                    }
                }
                let end = self.insts.len();
                for jump in jumps {
                    self.insts[jump] = Inst::Jump(end);
                }
            }
            Node::Subpattern {
                group,
                add_flags,
                del_flags,
                body,
            } => {
                if let Some(group) = group {
                    self.emit(Inst::Save(group * 2));
                }
                self.sequence(body, combine_flags(flags, *add_flags, *del_flags))?;
                if let Some(group) = group {
                    self.emit(Inst::Save(group * 2 + 1));
                }
            }
            Node::Atomic(body) => {
                let head = self.emit(Inst::Fail);
                self.sequence(body, flags)?;
                self.emit(Inst::SubEnd);
                self.insts[head] = Inst::Atomic {
                    body: head + 1,
                    exit: self.insts.len(),
                };
            }
            Node::Repeat {
                kind,
                min,
                max,
                body,
            } => {
                if let Some(unit) = self.simple_unit(body, flags) {
                    self.emit(Inst::RepeatOne {
                        unit,
                        min: *min,
                        max: *max,
                        kind: *kind,
                    });
                } else if *kind == RepeatKind::Possessive {
                    let head = self.emit(Inst::Fail);
                    self.sequence(body, flags)?;
                    self.emit(Inst::SubEnd);
                    self.insts[head] = Inst::PossessiveRepeat {
                        min: *min,
                        max: *max,
                        body: head + 1,
                        exit: self.insts.len(),
                    };
                } else {
                    let register = self.registers;
                    self.registers += 1;
                    self.emit(Inst::RepeatStart(register));
                    let head = self.emit(Inst::Fail);
                    self.sequence(body, flags)?;
                    self.emit(Inst::Jump(head));
                    self.insts[head] = Inst::Repeat {
                        register,
                        min: *min,
                        max: *max,
                        lazy: *kind == RepeatKind::Lazy,
                        body: head + 1,
                        exit: self.insts.len(),
                    };
                }
            }
            Node::GroupRef(group) => {
                self.emit(Inst::BackRef(*group, fold_of(flags)));
            }
            Node::GroupRefExists { group, yes, no } => {
                let head = self.emit(Inst::Fail);
                self.sequence(yes, flags)?;
                match no.as_deref().filter(|no| !no.is_empty()) {
                    Some(no) => {
                        let jump = self.emit(Inst::Jump(0));
                        let otherwise = self.insts.len();
                        self.sequence(no, flags)?;
                        self.insts[jump] = Inst::Jump(self.insts.len());
                        self.insts[head] = Inst::Conditional {
                            group: *group,
                            otherwise,
                        };
                    }
                    None => {
                        self.insts[head] = Inst::Conditional {
                            group: *group,
                            otherwise: self.insts.len(),
                        };
                    }
                }
            }
            Node::Assert {
                behind,
                negate,
                body,
            } => {
                let back = if *behind {
                    let (lo, hi) = width(body, self.widths);
                    if lo > MAXCODE {
                        return Err(PatternError("looks too much behind".to_string()));
                    }
                    if lo != hi {
                        return Err(PatternError(
                            "look-behind requires fixed-width pattern".to_string(),
                        ));
                    }
                    Some(lo as usize)
                } else {
                    None
                };
                let head = self.emit(Inst::Fail);
                self.sequence(body, flags)?;
                self.emit(Inst::SubEnd);
                self.insts[head] = Inst::Look {
                    behind: back,
                    negate: *negate,
                    body: head + 1,
                    exit: self.insts.len(),
                };
            }
            Node::Failure => {
                self.emit(Inst::Fail);
            }
        }
        Ok(())
    }
}

pub fn compile(parsed: &Parsed) -> Result<Program, PatternError> {
    let mut compiler = Compiler {
        insts: Vec::new(),
        sets: Vec::new(),
        registers: 0,
        widths: &parsed.widths,
    };
    compiler.sequence(&parsed.nodes, parsed.flags)?;
    compiler.emit(Inst::Match);
    let first = first_set(parsed).map(|set| {
        compiler.sets.push(set);
        compiler.sets.len() - 1
    });
    Ok(Program {
        insts: compiler.insts,
        sets: compiler.sets,
        groups: parsed.groups,
        registers: compiler.registers,
        names: parsed.names.clone(),
        first,
    })
}

fn is_cased(flags: u32, code: u32) -> bool {
    match fold_of(flags) {
        Fold::None => false,
        Fold::Ascii => ascii_is_cased(code),
        Fold::Unicode => unicode_is_cased(code),
    }
}

/// The literal characters a pattern starts with, and whether the pattern is those alone.
fn literal_prefix(nodes: &[Node], flags: u32) -> (Vec<u32>, bool) {
    let mut prefix = Vec::new();
    for node in nodes {
        match node {
            Node::Literal(code) if !is_cased(flags, *code) => prefix.push(*code),
            Node::Subpattern {
                add_flags,
                del_flags,
                body,
                ..
            } => {
                let (inner, whole) =
                    literal_prefix(body, combine_flags(flags, *add_flags, *del_flags));
                prefix.extend(inner);
                if !whole {
                    return (prefix, false);
                }
            }
            _ => return (prefix, false),
        }
    }
    (prefix, true)
}

/// The set of characters a match can start with, where CPython's search keeps one: for a
/// pattern that takes at least one character, starts with no literal and opens on a set,
/// on a character or on alternatives that each open on a character. The search skips every
/// position whose character is outside the set. Its categories are read by the pattern's
/// own flags, not by those of a group the set stands in, so a set under `(?a:` at the
/// start of a pattern is searched by the Unicode categories: CPython's behaviour, kept.
fn first_set(parsed: &Parsed) -> Option<CharSet> {
    if width(&parsed.nodes, &parsed.widths).0 == 0 {
        return None;
    }
    if !literal_prefix(&parsed.nodes, parsed.flags).0.is_empty() {
        return None;
    }
    let mut nodes = parsed.nodes.as_slice();
    let mut flags = parsed.flags;
    let mut first = nodes.first()?;
    while let Node::Subpattern {
        add_flags,
        del_flags,
        body,
        ..
    } = first
    {
        nodes = body.as_slice();
        flags = combine_flags(flags, *add_flags, *del_flags);
        first = nodes.first()?;
    }
    let (negate, items): (bool, Vec<SetItem>) = match first {
        Node::Literal(code) if !is_cased(flags, *code) => (false, vec![SetItem::Literal(*code)]),
        Node::Branch(alternatives) => {
            let mut items = Vec::new();
            for alternative in alternatives {
                match alternative.first() {
                    Some(Node::Literal(code)) if !is_cased(flags, *code) => {
                        items.push(SetItem::Literal(*code));
                    }
                    _ => return None,
                }
            }
            (false, items)
        }
        Node::In { negate, items } => {
            if fold_of(flags) != Fold::None {
                for item in items {
                    match item {
                        SetItem::Literal(code) if is_cased(flags, *code) => return None,
                        SetItem::Range(_, hi) if *hi > 0xffff => return None,
                        SetItem::Range(lo, hi) if (*lo..=*hi).any(|c| is_cased(flags, c)) => {
                            return None
                        }
                        _ => {}
                    }
                }
            }
            (*negate, items.clone())
        }
        _ => return None,
    };
    let unicode = parsed.flags & FLAG_UNICODE != 0;
    Some(CharSet {
        negate,
        fold: Fold::None,
        folded: Vec::new(),
        atoms: items
            .iter()
            .map(|item| match item {
                SetItem::Literal(code) => Atom::Literal(*code),
                SetItem::Range(lo, hi) => Atom::Range(*lo, *hi),
                SetItem::Category(category) => Atom::Category(*category, unicode),
            })
            .collect(),
    })
}

/// Whether a repeat count is the open bound.
pub fn unbounded(max: u32) -> bool {
    max == MAXREPEAT
}
