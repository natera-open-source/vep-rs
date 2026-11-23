// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The matcher: a backtracking run of a compiled pattern over a text of characters, with the
//! order of attempts CPython's `re` has, so the first match found is the one it finds. The
//! choices left open live on a stack of the heap, and a repeat of one character's worth of
//! pattern keeps one of them however often it repeats.

use super::compile::{
    ascii_lower, in_category, unbounded, unicode_lower, unit_matches, AtKind, Fold, Inst, Program,
};
use super::parse::{Category, RepeatKind};

const UNSET: usize = usize::MAX;

enum Undo {
    Capture(usize, usize),
    Register(usize, (u32, usize)),
}

enum ChoiceKind {
    /// Resume at the instruction as it stands.
    Plain,
    /// A greedy repeat of one unit that took `count` characters from `pos` and may give
    /// them back one at a time down to `min`.
    GiveBack { count: u32, min: u32 },
    /// A lazy repeat of one unit, at the instruction, that took `count` and may take more.
    TakeMore { count: u32 },
    /// A lazy repeat of a body, at its head, that may run the body once more.
    RunMore,
}

struct Choice {
    pc: usize,
    pos: usize,
    undo: usize,
    kind: ChoiceKind,
}

pub struct Vm<'a> {
    program: &'a Program,
    text: &'a [char],
    pub captures: Vec<usize>,
    /// Per repeat: the bodies run so far and where the last one started.
    registers: Vec<(u32, usize)>,
    undo: Vec<Undo>,
    stack: Vec<Choice>,
    start: usize,
    must_advance: bool,
    match_all: bool,
}

impl<'a> Vm<'a> {
    pub fn new(program: &'a Program, text: &'a [char]) -> Self {
        Vm {
            program,
            text,
            captures: vec![UNSET; program.groups * 2],
            registers: vec![(0, UNSET); program.registers],
            undo: Vec::new(),
            stack: Vec::new(),
            start: 0,
            must_advance: false,
            match_all: false,
        }
    }

    /// A match that starts at `start`: its end, with the groups in `captures`. With
    /// `must_advance` an empty match is refused; with `match_all` one that stops short of
    /// the text's end.
    pub fn match_at(&mut self, start: usize, must_advance: bool, match_all: bool) -> Option<usize> {
        self.captures.fill(UNSET);
        self.undo.clear();
        self.stack.clear();
        self.start = start;
        self.must_advance = must_advance;
        self.match_all = match_all;
        let end = self.run(0, start)?;
        self.captures[0] = start;
        self.captures[1] = end;
        Some(end)
    }

    /// Whether the character at a position is of a set of the program.
    pub fn starts_in(&self, set: usize, pos: usize) -> bool {
        self.char_at(pos)
            .is_some_and(|c| self.program.sets[set].contains(c))
    }

    pub fn group(&self, index: usize) -> Option<(usize, usize)> {
        let start = *self.captures.get(index * 2)?;
        let end = self.captures[index * 2 + 1];
        (start != UNSET && end != UNSET && end >= start).then_some((start, end))
    }

    fn unwind(&mut self, mark: usize) {
        while self.undo.len() > mark {
            match self.undo.pop().expect("an entry above the mark") {
                Undo::Capture(slot, old) => self.captures[slot] = old,
                Undo::Register(register, old) => self.registers[register] = old,
            }
        }
    }

    fn set_register(&mut self, register: usize, value: (u32, usize)) {
        self.undo
            .push(Undo::Register(register, self.registers[register]));
        self.registers[register] = value;
    }

    fn char_at(&self, pos: usize) -> Option<u32> {
        self.text.get(pos).map(|c| *c as u32)
    }

    fn is_word(&self, pos: usize, unicode: bool) -> bool {
        self.char_at(pos)
            .is_some_and(|c| in_category(Category::Word, unicode, c))
    }

    fn at(&self, kind: AtKind, pos: usize) -> bool {
        let len = self.text.len();
        let newline = '\n' as u32;
        match kind {
            AtKind::Beginning => pos == 0,
            AtKind::BeginningLine => pos == 0 || self.char_at(pos - 1) == Some(newline),
            AtKind::End => pos == len || (pos + 1 == len && self.char_at(pos) == Some(newline)),
            AtKind::EndLine => pos == len || self.char_at(pos) == Some(newline),
            AtKind::EndString => pos == len,
            AtKind::Boundary(unicode) | AtKind::NonBoundary(unicode) => {
                let before = pos > 0 && self.is_word(pos - 1, unicode);
                let after = self.is_word(pos, unicode);
                (before != after) == matches!(kind, AtKind::Boundary(_))
            }
        }
    }

    /// Runs from an instruction and a position to the end of the pattern or of the body the
    /// instruction opens; the choices the run left open are dropped when it ends.
    fn run(&mut self, mut pc: usize, mut pos: usize) -> Option<usize> {
        let program = self.program;
        let base = self.stack.len();
        loop {
            let advanced = match &program.insts[pc] {
                Inst::One(unit) => match self.char_at(pos) {
                    Some(c) if unit_matches(unit, &program.sets, c) => {
                        pos += 1;
                        pc += 1;
                        true
                    }
                    _ => false,
                },
                Inst::At(kind) => {
                    pc += 1;
                    self.at(*kind, pos)
                }
                Inst::Split(first, second) => {
                    self.stack.push(Choice {
                        pc: *second,
                        pos,
                        undo: self.undo.len(),
                        kind: ChoiceKind::Plain,
                    });
                    pc = *first;
                    true
                }
                Inst::Jump(target) => {
                    pc = *target;
                    true
                }
                Inst::Save(slot) => {
                    self.undo.push(Undo::Capture(*slot, self.captures[*slot]));
                    self.captures[*slot] = pos;
                    pc += 1;
                    true
                }
                Inst::RepeatStart(register) => {
                    self.set_register(*register, (0, UNSET));
                    pc += 1;
                    true
                }
                Inst::Repeat {
                    register,
                    min,
                    max,
                    lazy,
                    body,
                    exit,
                } => {
                    let (count, last) = self.registers[*register];
                    if count < *min {
                        self.set_register(*register, (count + 1, last));
                        pc = *body;
                    } else if *lazy {
                        self.stack.push(Choice {
                            pc,
                            pos,
                            undo: self.undo.len(),
                            kind: ChoiceKind::RunMore,
                        });
                        pc = *exit;
                    } else if (count < *max || unbounded(*max)) && pos != last {
                        // The body runs once more unless its last run took nothing.
                        self.stack.push(Choice {
                            pc: *exit,
                            pos,
                            undo: self.undo.len(),
                            kind: ChoiceKind::Plain,
                        });
                        self.set_register(*register, (count + 1, pos));
                        pc = *body;
                    } else {
                        pc = *exit;
                    }
                    true
                }
                Inst::RepeatOne {
                    unit,
                    min,
                    max,
                    kind,
                } => {
                    let available = (self.text.len() - pos).min(*max as usize);
                    let most = if *kind == RepeatKind::Lazy {
                        *min as usize
                    } else {
                        available
                    };
                    let mut count = 0;
                    while count < most
                        && count < available
                        && unit_matches(unit, &program.sets, self.text[pos + count] as u32)
                    {
                        count += 1;
                    }
                    if count < *min as usize {
                        false
                    } else {
                        match kind {
                            RepeatKind::Greedy => {
                                if count > *min as usize {
                                    self.stack.push(Choice {
                                        pc: pc + 1,
                                        pos,
                                        undo: self.undo.len(),
                                        kind: ChoiceKind::GiveBack {
                                            count: count as u32,
                                            min: *min,
                                        },
                                    });
                                }
                            }
                            RepeatKind::Lazy => self.stack.push(Choice {
                                pc,
                                pos: pos + count,
                                undo: self.undo.len(),
                                kind: ChoiceKind::TakeMore {
                                    count: count as u32,
                                },
                            }),
                            RepeatKind::Possessive => {}
                        }
                        pos += count;
                        pc += 1;
                        true
                    }
                }
                Inst::PossessiveRepeat {
                    min,
                    max,
                    body,
                    exit,
                } => {
                    let mark = self.undo.len();
                    let mut count = 0;
                    let mut failed = false;
                    while count < *min {
                        match self.run(*body, pos) {
                            Some(end) => {
                                pos = end;
                                count += 1;
                            }
                            None => {
                                failed = true;
                                break;
                            }
                        }
                    }
                    if failed {
                        self.unwind(mark);
                        false
                    } else {
                        let mut last = UNSET;
                        while (count < *max || unbounded(*max)) && pos != last {
                            let before = self.undo.len();
                            last = pos;
                            match self.run(*body, pos) {
                                Some(end) => {
                                    pos = end;
                                    count += 1;
                                }
                                None => {
                                    self.unwind(before);
                                    break;
                                }
                            }
                        }
                        pc = *exit;
                        true
                    }
                }
                Inst::Atomic { body, exit } => {
                    let mark = self.undo.len();
                    match self.run(*body, pos) {
                        Some(end) => {
                            pos = end;
                            pc = *exit;
                            true
                        }
                        None => {
                            self.unwind(mark);
                            false
                        }
                    }
                }
                Inst::Look {
                    behind,
                    negate,
                    body,
                    exit,
                } => {
                    let mark = self.undo.len();
                    let from = match behind {
                        None => Some(pos),
                        Some(back) => pos.checked_sub(*back),
                    };
                    let found = from.is_some_and(|from| self.run(*body, from).is_some());
                    if found != *negate {
                        if *negate {
                            self.unwind(mark);
                        }
                        pc = *exit;
                        true
                    } else {
                        self.unwind(mark);
                        false
                    }
                }
                Inst::SubEnd => {
                    self.stack.truncate(base);
                    return Some(pos);
                }
                Inst::BackRef(group, fold) => match self.group(*group) {
                    None => false,
                    Some((start, end)) => {
                        let length = end - start;
                        let same = |a: char, b: char| match fold {
                            Fold::None => a == b,
                            Fold::Ascii => ascii_lower(a as u32) == ascii_lower(b as u32),
                            Fold::Unicode => unicode_lower(a as u32) == unicode_lower(b as u32),
                        };
                        if pos + length <= self.text.len()
                            && (0..length).all(|i| same(self.text[pos + i], self.text[start + i]))
                        {
                            pos += length;
                            pc += 1;
                            true
                        } else {
                            false
                        }
                    }
                },
                Inst::Conditional { group, otherwise } => {
                    pc = if self.group(*group).is_some() {
                        pc + 1
                    } else {
                        *otherwise
                    };
                    true
                }
                Inst::Fail => false,
                Inst::Match => {
                    let refused = (self.match_all && pos != self.text.len())
                        || (self.must_advance && pos == self.start);
                    if !refused {
                        self.stack.truncate(base);
                        return Some(pos);
                    }
                    false
                }
            };
            if advanced {
                continue;
            }
            // Back to the newest open choice of this run.
            loop {
                if self.stack.len() == base {
                    return None;
                }
                let top = self.stack.len() - 1;
                let mark = self.stack[top].undo;
                self.unwind(mark);
                let choice = &mut self.stack[top];
                match &mut choice.kind {
                    ChoiceKind::Plain => {
                        pc = choice.pc;
                        pos = choice.pos;
                        self.stack.pop();
                        break;
                    }
                    ChoiceKind::GiveBack { count, min } => {
                        *count -= 1;
                        pc = choice.pc;
                        pos = choice.pos + *count as usize;
                        if *count == *min {
                            self.stack.pop();
                        }
                        break;
                    }
                    ChoiceKind::TakeMore { count } => {
                        let Inst::RepeatOne { unit, max, .. } = &program.insts[choice.pc] else {
                            unreachable!("a lazy repeat of one unit")
                        };
                        let here = choice.pos;
                        let more = (*count < *max || unbounded(*max))
                            && self
                                .text
                                .get(here)
                                .is_some_and(|c| unit_matches(unit, &program.sets, *c as u32));
                        if more {
                            *count += 1;
                            choice.pos += 1;
                            pc = choice.pc + 1;
                            pos = choice.pos;
                            break;
                        }
                        self.stack.pop();
                    }
                    ChoiceKind::RunMore => {
                        let Inst::Repeat {
                            register,
                            max,
                            body,
                            ..
                        } = &program.insts[choice.pc]
                        else {
                            unreachable!("a lazy repeat of a body")
                        };
                        let here = choice.pos;
                        self.stack.pop();
                        let (count, last) = self.registers[*register];
                        if (count >= *max && !unbounded(*max)) || here == last {
                            continue;
                        }
                        self.set_register(*register, (count + 1, here));
                        pc = *body;
                        pos = here;
                        break;
                    }
                }
            }
        }
    }
}
