// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A run over a cache. Reading threads take the sequences in the cache's order, one
//! sequence a thread, and each hands on what it read and then the sequence's rows, a few
//! blocks ahead of the writer. The calling thread takes the sequences in order, writes
//! the rows to the table and passes the transcript records to a thread that writes the
//! database. Each file is written under a temporary name beside its own and renamed
//! when the run is whole, so a failed run leaves no file.

use std::collections::VecDeque;
use std::fs;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::{Mutex, PoisonError};
use std::thread;

use vep_seeded_variants::cache::Cache;

use crate::database::{read_rows, DatabaseRows, DatabaseWriter};
use crate::exons::{Record, SequenceExons};
use crate::table::{holds_name, rows, Block, TableRows, TableWriter, HIGHEST_POSITION};
use crate::{DATABASE_FILE, INDEX_FILE, TABLE_FILE};

/// How many blocks of rows a reading thread holds ready for the writer.
const BLOCKS_AHEAD: usize = 16;

/// How many sequences' transcript records wait for the database's thread.
const RECORDS_AHEAD: usize = 4;

const ENDED_EARLY: &str = "a thread of the run ended early";

/// What a run is asked for.
#[derive(Debug, Clone)]
pub struct Options {
    /// The seed every value derives from.
    pub seed: u64,
    /// The reading threads and the compressing threads; `None` takes one a processor.
    pub threads: Option<usize>,
    /// Leaves out the transcript records the reader gives no exon for, where a run
    /// without the option fails on the first.
    pub skip_transcripts_without_exons: bool,
    /// Overwrites a file that exists.
    pub force: bool,
}

/// The counts of one sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceCounts {
    /// The sequence's name.
    pub name: String,
    /// The transcript records the reader gives.
    pub transcripts: u64,
    /// The exons of the records that have one.
    pub exons: u64,
    /// The bases that lie in an exon, each once: the sequence's rows in the table.
    pub exon_bases: u64,
}

/// What a run wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Summary {
    /// The seed.
    pub seed: u64,
    /// The counts of every sequence, in the table's order.
    pub sequences: Vec<SequenceCounts>,
    /// The records left out, under the option that leaves them out.
    pub transcripts_without_exons: Option<u64>,
    /// The table's rows and their digest.
    pub table: TableRows,
    /// The database's rows and their digest.
    pub database: DatabaseRows,
}

impl Summary {
    fn total(&self, count: impl Fn(&SequenceCounts) -> u64) -> u64 {
        self.sequences.iter().map(count).sum()
    }

    /// The lines a run prints: a key, a tab and the value, then one line a file.
    pub fn lines(&self) -> String {
        let mut lines = format!(
            "seed\t{}\nsequences\t{}\ntranscripts\t{}\n",
            self.seed,
            self.sequences.len(),
            self.total(|sequence| sequence.transcripts)
        );
        if let Some(without) = self.transcripts_without_exons {
            lines.push_str(&format!("transcripts_without_exons\t{without}\n"));
        }
        lines.push_str(&format!(
            "exons\t{}\nexon_bases\t{}\n",
            self.total(|sequence| sequence.exons),
            self.total(|sequence| sequence.exon_bases)
        ));
        lines.push_str(&format!(
            "{TABLE_FILE}\trows={}\trows_sha256={}\n",
            self.table.rows, self.table.sha256
        ));
        lines.push_str(&format!(
            "{DATABASE_FILE}\trows={}\tgerp_exons={}\tphylocsf_data={}\trows_sha256={}\n",
            self.database.gerp_exons + self.database.phylocsf_data,
            self.database.gerp_exons,
            self.database.phylocsf_data,
            self.database.sha256
        ));
        lines
    }

    /// The table of the sequences: the header, then one line a sequence in the table's
    /// order.
    pub fn per_sequence(&self) -> String {
        let mut table = "sequence\ttranscripts\texons\texon_bases\n".to_owned();
        for sequence in &self.sequences {
            table.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                sequence.name, sequence.transcripts, sequence.exons, sequence.exon_bases
            ));
        }
        table
    }
}

/// What a reading thread hands on for one sequence: first what it read, then the
/// sequence's rows.
enum Message {
    Read(Result<SequenceExons, String>),
    Rows(Block),
}

/// The three files under the names they are written under, each removed when the value
/// is dropped unless the run renamed it.
struct Temporary {
    table: PathBuf,
    index: PathBuf,
    database: PathBuf,
}

impl Temporary {
    fn files(&self) -> [&PathBuf; 3] {
        [&self.table, &self.index, &self.database]
    }
}

impl Drop for Temporary {
    fn drop(&mut self) {
        for file in self.files() {
            let _ = fs::remove_file(file);
        }
    }
}

/// The path of `file` under `dir` with `ending` after its name.
fn beside(dir: &Path, file: &str, ending: &str) -> PathBuf {
    dir.join(format!("{file}{ending}"))
}

/// Reads one sequence for a reading thread: `None` when the run stopped meanwhile.
fn read_sequence(
    cache: &Cache,
    sequence: &str,
    stopped: &AtomicBool,
) -> Option<Result<SequenceExons, String>> {
    if !holds_name(sequence) {
        return Some(Err(format!(
            "the sequence name {sequence:?} cannot stand in a row of the table"
        )));
    }
    let read = match SequenceExons::read(cache, sequence, stopped) {
        Ok(read) => read?,
        Err(error) => return Some(Err(error.to_string())),
    };
    match read.spans.last() {
        Some(last) if last.end > HIGHEST_POSITION => Some(Err(format!(
            "sequence {sequence} has an exon base at {}, above the {HIGHEST_POSITION} a tabix \
             index reaches",
            last.end
        ))),
        _ => Some(Ok(read)),
    }
}

/// What the calling thread gathers while it writes the table.
struct Written {
    sequences: Vec<SequenceCounts>,
    without_exons: u64,
    table: TableRows,
}

/// Takes the sequences in order from the reading threads: their rows go to the table
/// and their transcript records to the database's thread. After the first record
/// without an exon, in a run that does not leave such records out, the sequences are
/// read on for their count alone.
fn write_table(
    sequences: &[String],
    receivers: Vec<Receiver<Message>>,
    records: &SyncSender<Vec<Record>>,
    temporary: &Temporary,
    workers: NonZero<usize>,
    options: &Options,
) -> Result<Written, String> {
    let mut table = TableWriter::create(&temporary.table, workers)?;
    let mut counts = Vec::with_capacity(sequences.len());
    let mut without_exons = 0;
    let mut first_without_exons = None;
    for (name, receiver) in sequences.iter().zip(receivers) {
        let read = match receiver.recv() {
            Ok(Message::Read(read)) => read?,
            _ => return Err(ENDED_EARLY.to_owned()),
        };
        let sequence = SequenceCounts {
            name: name.clone(),
            transcripts: read.transcripts(),
            exons: read.exons(),
            exon_bases: read.exon_bases(),
        };
        without_exons += read.without_exons;
        if first_without_exons.is_none() {
            first_without_exons = read.first_without_exons.map(|id| (id, name));
        }
        if without_exons > 0 && !options.skip_transcripts_without_exons {
            continue;
        }
        records
            .send(read.records)
            .map_err(|_| ENDED_EARLY.to_owned())?;
        let mut written = 0;
        for message in receiver {
            let Message::Rows(block) = message else {
                return Err(ENDED_EARLY.to_owned());
            };
            table.write(name, &block)?;
            written += block.rows();
        }
        if written != sequence.exon_bases {
            return Err(ENDED_EARLY.to_owned());
        }
        counts.push(sequence);
    }
    if let Some((id, sequence)) = first_without_exons {
        if !options.skip_transcripts_without_exons {
            return Err(format!(
                "transcript records the cache reader gives no exon for: {without_exons}, the \
                 first {id} on sequence {sequence}; --skip-transcripts-without-exons writes \
                 the files without them"
            ));
        }
    }
    Ok(Written {
        sequences: counts,
        without_exons,
        table: table.finish(&temporary.index)?,
    })
}

/// Writes the three files for `cache` into `out_dir`, which is created where it is
/// absent: `GERP_scores.synthetic.txt.gz`, its `.tbi` index and
/// `phylocsf_gerp.synthetic.sql`.
pub fn generate(cache: &Cache, out_dir: &Path, options: &Options) -> Result<Summary, String> {
    fs::create_dir_all(out_dir).map_err(|error| format!("{}: {error}", out_dir.display()))?;
    let names = [TABLE_FILE, INDEX_FILE, DATABASE_FILE];
    if !options.force {
        for name in names {
            let file = out_dir.join(name);
            if file.exists() {
                return Err(format!("{} exists; --force overwrites it", file.display()));
            }
        }
    }
    // A rename over a directory fails, after the files renamed before it are in place.
    for name in names {
        let file = out_dir.join(name);
        if fs::symlink_metadata(&file).is_ok_and(|entry| entry.is_dir()) {
            return Err(format!("{} is a directory", file.display()));
        }
    }
    let ending = format!(".{}.tmp", std::process::id());
    let temporary = Temporary {
        table: beside(out_dir, TABLE_FILE, &ending),
        index: beside(out_dir, INDEX_FILE, &ending),
        database: beside(out_dir, DATABASE_FILE, &ending),
    };
    // SQLite opens a file a killed run left under this process id as a database to add to.
    for file in temporary.files() {
        let _ = fs::remove_file(file);
    }

    let sequences = cache.sequences();
    let workers = options
        .threads
        .or_else(|| thread::available_parallelism().ok().map(usize::from))
        .and_then(NonZero::new)
        .unwrap_or(NonZero::<usize>::MIN);
    let readers = workers.get().min(sequences.len());
    let stopped = AtomicBool::new(false);
    let (senders, receivers): (Vec<_>, Vec<_>) = sequences
        .iter()
        .map(|_| sync_channel::<Message>(BLOCKS_AHEAD))
        .unzip();
    let queue: Mutex<VecDeque<(&String, SyncSender<Message>)>> =
        Mutex::new(sequences.iter().zip(senders).collect());
    let seed = options.seed;
    let writes_without = options.skip_transcripts_without_exons;

    let read = || loop {
        let next = queue
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .pop_front();
        let Some((name, sender)) = next else {
            return;
        };
        let Some(read) = read_sequence(cache, name, &stopped) else {
            return;
        };
        let spans = match &read {
            Ok(read) if writes_without || read.without_exons == 0 => read.spans.clone(),
            _ => Vec::new(),
        };
        if sender.send(Message::Read(read)).is_ok() {
            rows(seed, name, &spans, |block| {
                sender.send(Message::Rows(block)).is_ok()
            });
        }
    };

    let (written, database) = thread::scope(|scope| {
        for _ in 0..readers {
            scope.spawn(read);
        }
        let (records, queued) = sync_channel::<Vec<Record>>(RECORDS_AHEAD);
        let database = scope.spawn(|| -> Result<Option<DatabaseRows>, String> {
            let wrong =
                |error: rusqlite::Error| format!("{}: {error}", temporary.database.display());
            let mut writer = DatabaseWriter::create(&temporary.database, seed)?;
            for records in queued {
                writer.add(&records).map_err(wrong)?;
            }
            if stopped.load(Ordering::Relaxed) {
                return Ok(None);
            }
            writer.finish().map_err(wrong)?;
            read_rows(&temporary.database).map(Some)
        });
        let written = write_table(sequences, receivers, &records, &temporary, workers, options);
        if written.is_err() {
            stopped.store(true, Ordering::Relaxed);
        }
        drop(records);
        let database = database
            .join()
            .unwrap_or_else(|_| Err(ENDED_EARLY.to_owned()));
        // A database that failed is why the table's thread could hand it no record.
        match (written, database) {
            (_, Err(what)) | (Err(what), _) => Err(what),
            (Ok(_), Ok(None)) => Err(ENDED_EARLY.to_owned()),
            (Ok(written), Ok(Some(database))) => Ok((written, database)),
        }
    })?;

    let exons: u64 = written
        .sequences
        .iter()
        .map(|sequence| sequence.exons)
        .sum();
    let exon_bases: u64 = written
        .sequences
        .iter()
        .map(|sequence| sequence.exon_bases)
        .sum();
    if written.table.rows != exon_bases
        || database.gerp_exons != exons
        || database.phylocsf_data != exons
    {
        return Err(format!(
            "the files hold {} rows, {} totals and {} PhyloCSF rows for {exon_bases} exon bases \
             and {exons} exons",
            written.table.rows, database.gerp_exons, database.phylocsf_data
        ));
    }
    for (from, name) in temporary.files().into_iter().zip(names) {
        let to = out_dir.join(name);
        fs::rename(from, &to).map_err(|error| format!("{}: {error}", to.display()))?;
    }
    Ok(Summary {
        seed,
        sequences: written.sequences,
        transcripts_without_exons: options
            .skip_transcripts_without_exons
            .then_some(written.without_exons),
        table: written.table,
        database,
    })
}
