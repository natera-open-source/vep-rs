// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The first read of the table: its header and, for every row, the numbers of the four
//! values its class depends on and the engine that wrote it. The text of a row is not
//! kept.
//!
//! The calling thread reads the text and cuts it into records; the values of the four
//! columns go in batches to the threads that number them, each column in row order on one
//! thread, so the numbers are those of one thread at any bound.

use crate::csv::{Reader, Record};
use crate::intern::{Batch, Bits, Interner, MOST_ROWS};
use crate::table::{is_reference_only, start_of, Header};
use crate::text::{self, Text};
use std::io::{self, Read};
use std::sync::mpsc::{sync_channel, SyncSender};

/// The columns numbered: Location, Allele, Feature and the consequence set.
const COLUMNS: usize = 4;
const LOCATION: usize = 0;

/// The rows of one batch.
const BATCH_ROWS: usize = 1 << 15;

/// The batches that may wait for one numbering thread.
const WAITING_BATCHES: usize = 4;

/// Why a read of the table ends early.
#[derive(Debug)]
pub enum Failure {
    Read(io::Error),
    NotUtf8(u64),
    /// The table itself is refused; the text says why.
    Table(String),
}

impl From<text::Failure> for Failure {
    fn from(failure: text::Failure) -> Failure {
        match failure {
            text::Failure::Read(error) => Failure::Read(error),
            text::Failure::NotUtf8(offset) => Failure::NotUtf8(offset),
        }
    }
}

impl From<&'static str> for Failure {
    fn from(detail: &'static str) -> Failure {
        Failure::Table(detail.to_string())
    }
}

/// What the classes of a table are derived from, a number a value.
pub struct Keys {
    pub rows: usize,
    /// The rows that the reference alone wrote.
    pub reference: Bits,
    pub location: Vec<u32>,
    pub allele: Vec<u32>,
    pub feature: Vec<u32>,
    pub consequence: Vec<u32>,
    /// The number of the start of each distinct Location.
    pub start_of_location: Vec<u32>,
}

impl Keys {
    /// The number of the start of a row's Location.
    pub fn start(&self, row: usize) -> u32 {
        self.start_of_location[self.location[row] as usize]
    }
}

/// One column as it is numbered.
struct Column {
    numbers: Vec<u32>,
    values: Interner,
    /// For the Location column: the starts, and the number of the start of each Location.
    starts: Option<(Interner, Vec<u32>)>,
}

impl Column {
    fn new(with_starts: bool) -> Column {
        Column {
            numbers: Vec::new(),
            values: Interner::new(),
            starts: with_starts.then(|| (Interner::new(), Vec::new())),
        }
    }

    fn absorb(&mut self, batch: &Batch) {
        for value in batch.values() {
            let number = self.values.intern_optional(value);
            if let (Some((starts, of_location)), Some(value)) = (&mut self.starts, value) {
                if number as usize == of_location.len() {
                    of_location.push(starts.intern(start_of(value)));
                }
            }
            self.numbers.push(number);
        }
    }
}

/// Where the batches of a run go: to the columns on the calling thread, or to the threads
/// that hold them.
enum Numbering {
    Here(Vec<Column>),
    Threads(Vec<SyncSender<(usize, Batch)>>),
}

struct Rows {
    header: Option<Header>,
    rows: usize,
    reference: Bits,
    batches: [Batch; COLUMNS],
    numbering: Numbering,
}

impl Rows {
    fn record(&mut self, record: Record<'_>) -> Result<(), Failure> {
        let Some(header) = &self.header else {
            self.header = Some(Header::read(&record).map_err(Failure::Table)?);
            return Ok(());
        };
        if record.is_empty() {
            return Ok(());
        }
        let Some(location) = record.get(header.location) else {
            return Err(Failure::Table(format!(
                "row {} ends before its location column",
                self.rows + 1
            )));
        };
        if self.rows == MOST_ROWS {
            return Err(Failure::Table(format!("more than {MOST_ROWS} rows")));
        }
        self.reference
            .push(is_reference_only(record.get(header.source)));
        let values = [
            Some(location),
            record.get(header.allele),
            record.get(header.feature),
            record.get(header.consequence),
        ];
        for (batch, value) in self.batches.iter_mut().zip(values) {
            batch.push(value);
        }
        self.rows += 1;
        if self.rows.is_multiple_of(BATCH_ROWS) {
            self.hand_on();
        }
        Ok(())
    }

    fn hand_on(&mut self) {
        let batches = std::mem::take(&mut self.batches);
        match &mut self.numbering {
            Numbering::Here(columns) => {
                for (column, batch) in columns.iter_mut().zip(&batches) {
                    column.absorb(batch);
                }
            }
            Numbering::Threads(senders) => {
                for (index, batch) in batches.into_iter().enumerate() {
                    // A numbering thread ends only when its sender is dropped.
                    let _ = senders[index % senders.len()].send((index, batch));
                }
            }
        }
    }
}

/// Reads the table once. `threads` bounds the threads at work, the calling one included.
pub fn scan<R: Read>(mut text: Text<R>, threads: usize) -> Result<(Header, Keys), Failure> {
    let columns: Vec<Column> = (0..COLUMNS)
        .map(|index| Column::new(index == LOCATION))
        .collect();
    let numbering_threads = threads.saturating_sub(1).min(COLUMNS);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        let numbering = if numbering_threads == 0 {
            Numbering::Here(columns)
        } else {
            let mut shares: Vec<Vec<(usize, Column)>> =
                (0..numbering_threads).map(|_| Vec::new()).collect();
            for (index, column) in columns.into_iter().enumerate() {
                shares[index % numbering_threads].push((index, column));
            }
            let mut senders = Vec::new();
            for mut share in shares {
                let (sender, receiver) = sync_channel::<(usize, Batch)>(WAITING_BATCHES);
                senders.push(sender);
                handles.push(scope.spawn(move || {
                    for (index, batch) in receiver {
                        if let Some((_, column)) = share.iter_mut().find(|(held, _)| *held == index)
                        {
                            column.absorb(&batch);
                        }
                    }
                    share
                }));
            }
            Numbering::Threads(senders)
        };
        let mut rows = Rows {
            header: None,
            rows: 0,
            reference: Bits::default(),
            batches: Default::default(),
            numbering,
        };
        let read = read_records(&mut text, &mut rows);
        rows.hand_on();
        let mut columns: Vec<(usize, Column)> = match rows.numbering {
            Numbering::Here(columns) => columns.into_iter().enumerate().collect(),
            Numbering::Threads(senders) => {
                drop(senders);
                handles
                    .into_iter()
                    .flat_map(|handle| handle.join().expect("a numbering thread"))
                    .collect()
            }
        };
        read?;
        let header = rows
            .header
            .ok_or_else(|| Failure::Table("no header row".to_string()))?;
        columns.sort_by_key(|(index, _)| *index);
        let mut columns = columns.into_iter().map(|(_, column)| column);
        let mut next = || columns.next().expect("four columns");
        let location = next();
        let (allele, feature, consequence) = (next(), next(), next());
        let (_, start_of_location) = location.starts.expect("the starts of the Location column");
        Ok((
            header,
            Keys {
                rows: rows.rows,
                reference: rows.reference,
                location: location.numbers,
                allele: allele.numbers,
                feature: feature.numbers,
                consequence: consequence.numbers,
                start_of_location,
            },
        ))
    })
}

fn read_records<R: Read>(text: &mut Text<R>, rows: &mut Rows) -> Result<(), Failure> {
    let mut reader = Reader::new();
    let mut run = Vec::new();
    while text.next_run(&mut run)? {
        reader.feed(&run, &mut |record| rows.record(record))?;
    }
    reader.finish(&mut |record| rows.record(record))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scanned(table: &str, threads: usize) -> Result<(Header, Keys), Failure> {
        scan(Text::with_run_bytes(table.as_bytes(), 7), threads)
    }

    const TABLE: &str =
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n\
        o\tmissing_in_rust\t1:100-200\tA\tT1\tTranscript\tx\n\
        \n\
        o\textra_in_rust\t1:100-300\tA\tT1\tTranscript\tx\n\
        o\tonly_perl\t2:100\t\tT2\n\
        o\tother\t1:100-200\n";

    #[test]
    fn equal_values_share_a_number_and_an_absent_one_has_none_at_any_bound() {
        for threads in [1, 2, 3, 4, 5, 64] {
            let (header, keys) = scanned(TABLE, threads).unwrap();
            assert_eq!(header.names.len(), 7);
            assert_eq!(keys.rows, 4, "an empty line is no row");
            let reference: Vec<bool> = (0..4).map(|row| keys.reference.get(row)).collect();
            assert_eq!(reference, [true, false, true, false]);
            assert_eq!(keys.location, [0, 1, 2, 0]);
            assert_eq!(keys.allele, [0, 0, 1, crate::intern::ABSENT]);
            assert_eq!(keys.feature, [0, 0, 1, crate::intern::ABSENT]);
            assert_eq!(
                keys.consequence,
                [0, 0, crate::intern::ABSENT, crate::intern::ABSENT]
            );
            assert_eq!(keys.start_of_location, [0, 0, 0]);
            assert_eq!(
                (0..4).map(|row| keys.start(row)).collect::<Vec<_>>(),
                [0; 4]
            );
        }
    }

    #[test]
    fn a_table_of_many_batches_is_numbered_alike_at_any_bound() {
        let mut table =
            String::from("file_name\tsource\tlocation\tallele\tfeature\tconsequence_set\n");
        let rows = 2 * BATCH_ROWS + 17;
        for row in 0..rows {
            table.push_str(&format!(
                "f\textra_in_rust\t{}:{}-{}\tA{}\tT{}\tc{}\n",
                row % 5,
                row % 1000,
                row,
                row % 3,
                row % 77,
                row % 2
            ));
        }
        let (_, one) = scanned(&table, 1).unwrap();
        assert_eq!(one.rows, rows);
        for threads in [2, 3, 5] {
            let (_, keys) = scanned(&table, threads).unwrap();
            assert_eq!(keys.location, one.location);
            assert_eq!(keys.allele, one.allele);
            assert_eq!(keys.feature, one.feature);
            assert_eq!(keys.consequence, one.consequence);
            assert_eq!(keys.start_of_location, one.start_of_location);
        }
        assert_eq!(one.start_of_location.iter().max(), Some(&999));
    }

    #[test]
    fn a_table_is_refused_by_what_is_wrong_with_it() {
        let detail = |table: &str| match scanned(table, 2) {
            Err(Failure::Table(detail)) => detail,
            other => panic!("{table:?}: {:?}", other.map(|(header, _)| header)),
        };
        assert_eq!(detail(""), "no header row");
        assert_eq!(
            detail("\nfile_name\n"),
            "Unsupported discordant TSV schema; missing required columns: \
             allele, consequence_set, feature, file_name, location, source"
        );
        let header = "file_name\tsource\tlocation\tallele\tfeature\tconsequence_set\n";
        assert_eq!(
            detail(&format!("{header}o\ts\t1:10\nf\ts\n")),
            "row 2 ends before its location column"
        );
        assert_eq!(
            detail(&format!("{header}o\ts\t{}\n", "x".repeat(131_073))),
            "field larger than field limit (131072)"
        );
        assert!(matches!(
            scanned(&format!("{header}o\ts\t1:10\n"), 2),
            Ok((_, Keys { rows: 1, .. }))
        ));
        assert!(matches!(scanned(header, 2), Ok((_, Keys { rows: 0, .. }))));
    }

    #[test]
    fn a_byte_outside_utf_8_is_named_by_its_offset() {
        let table = b"file_name\tsource\tlocation\tallele\tfeature\tconsequence_set\no\ts\t\xff\n";
        match scan(Text::with_run_bytes(&table[..], 5), 3) {
            Err(Failure::NotUtf8(offset)) => assert_eq!(offset, 61),
            other => panic!("{:?}", other.map(|(header, _)| header)),
        }
    }
}
