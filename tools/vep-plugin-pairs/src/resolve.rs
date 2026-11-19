// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The dialogue with the driver, which owns the plugin registry: the source of a field
//! name and the canonical form of a Location the tool leaves to it, each asked once.
//!
//! A question is one line on standard output, `field <hex>` or `location <hex>`, and its
//! answer one line on standard input: `source <name>`, `core`, `unclaimed` or `refused`
//! for a field, `location <hex>` for a Location.

use crate::order::{cmp_column, hex, unhex};
use crate::{Refusal, Res};
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Write};
use std::sync::{Mutex, RwLock};

/// What a field name is to the scored run.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Source {
    /// A column of Ensembl VEP's own.
    Core,
    /// A key no source claims.
    Unclaimed,
    /// A field of a source, by its number in the order the fields are first met.
    Field(u32),
}

struct Dialogue {
    ask: Box<dyn Write + Send>,
    answer: Box<dyn BufRead + Send>,
}

impl Dialogue {
    fn put(&mut self, question: &str, subject: &[u8]) -> Res<Vec<u8>> {
        let gone = || Refusal("the driver left the dialogue".to_string());
        writeln!(self.ask, "{question} {}", hex(subject)).map_err(|_| gone())?;
        self.ask.flush().map_err(|_| gone())?;
        let mut line = Vec::new();
        if self
            .answer
            .read_until(b'\n', &mut line)
            .map_err(|_| gone())?
            == 0
        {
            return Err(gone());
        }
        while line.last().is_some_and(|byte| b"\r\n".contains(byte)) {
            line.pop();
        }
        Ok(line)
    }
}

#[derive(Default)]
struct Claimed {
    names: Vec<Vec<u8>>,
    source: Vec<u32>,
    sources: Vec<String>,
}

/// The fields met so far and the dialogue that resolves a new one.
pub struct Registry {
    resolved: RwLock<HashMap<Vec<u8>, Source>>,
    claimed: Mutex<Claimed>,
    dialogue: Mutex<Dialogue>,
}

/// The claimed fields at one moment: per field its name, its source and its place among the
/// names in the order tuple lines sort them.
pub struct Table {
    pub names: Vec<Vec<u8>>,
    pub source: Vec<u32>,
    pub rank: Vec<u32>,
    pub sources: Vec<String>,
}

impl Registry {
    pub fn over_standard_streams() -> Registry {
        Registry {
            resolved: RwLock::default(),
            claimed: Mutex::default(),
            dialogue: Mutex::new(Dialogue {
                ask: Box::new(io::stdout()),
                answer: Box::new(BufReader::new(io::stdin())),
            }),
        }
    }

    /// The source of `name`, asked of the driver the first time the name is met.
    pub fn source_of(&self, name: &[u8]) -> Res<Source> {
        if let Some(&source) = self.resolved.read().unwrap().get(name) {
            return Ok(source);
        }
        let mut dialogue = self.dialogue.lock().unwrap();
        if let Some(&source) = self.resolved.read().unwrap().get(name) {
            return Ok(source);
        }
        let answer = dialogue.put("field", name)?;
        let source = match answer.as_slice() {
            b"core" => Source::Core,
            b"unclaimed" => Source::Unclaimed,
            b"refused" => {
                return Err(Refusal(format!(
                    "the registry refuses the field {}",
                    hex(name)
                )))
            }
            _ => match answer.strip_prefix(b"source ") {
                Some(source) => {
                    let source = String::from_utf8_lossy(source).into_owned();
                    let mut claimed = self.claimed.lock().unwrap();
                    let known = claimed.sources.iter().position(|other| *other == source);
                    let number = known.unwrap_or_else(|| {
                        claimed.sources.push(source);
                        claimed.sources.len() - 1
                    });
                    claimed.names.push(name.to_vec());
                    claimed.source.push(number as u32);
                    Source::Field(claimed.names.len() as u32 - 1)
                }
                None => return Err(unreadable(&answer)),
            },
        };
        self.resolved.write().unwrap().insert(name.to_vec(), source);
        Ok(source)
    }

    /// The canonical form of a Location that `canonical_location` leaves to the driver.
    pub fn location(&self, written: &[u8]) -> Res<Vec<u8>> {
        let answer = self.dialogue.lock().unwrap().put("location", written)?;
        answer
            .strip_prefix(b"location ")
            .and_then(unhex)
            .ok_or_else(|| unreadable(&answer))
    }

    pub fn table(&self) -> Table {
        let claimed = self.claimed.lock().unwrap();
        let mut order: Vec<u32> = (0..claimed.names.len() as u32).collect();
        order.sort_by(|&a, &b| cmp_column(&claimed.names[a as usize], &claimed.names[b as usize]));
        let mut rank = vec![0; order.len()];
        for (place, &field) in order.iter().enumerate() {
            rank[field as usize] = place as u32;
        }
        Table {
            names: claimed.names.clone(),
            source: claimed.source.clone(),
            rank,
            sources: claimed.sources.clone(),
        }
    }
}

fn unreadable(answer: &[u8]) -> Refusal {
    Refusal(format!(
        "the driver's answer cannot be read: {}",
        String::from_utf8_lossy(answer)
    ))
}
