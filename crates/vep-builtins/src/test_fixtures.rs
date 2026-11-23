// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Test-only helpers for the plugin tests: a real bgzip+tabix fixture builder
//! and a subscriber that collects the messages of `warn!` events.
//!
//! The fixture helpers write genuine bgzip-compressed, tabix-indexed files so tests
//! exercise the bgzf and `.tbi` read path in [`crate::tabix`], the path every tabix
//! plugin runs through. Two defect classes live only there and are invisible to
//! a test built on in-memory strings:
//!
//! - a data file whose header lacks a leading `#` leaves `header = None`, so
//!   every score lookup by name returns `None` and the plugin annotates
//!   nothing while exiting 0 (the dbscSNV plugin documents a header line that
//!   opens with `c`, dbscSNV.pm:61-62);
//! - a file carrying two coordinate columns read on the first one discards
//!   every record when indexed on the second (REVEL on GRCh38).
//!
//! Fixtures are tiny (a handful of rows), so they cost milliseconds and need no
//! network or external `bgzip`/`tabix` binaries.
//!
//! [`Warnings`] collects the `warn!` events a test raises, the texts the run
//! would write to its warnings file.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, Once};

use noodles::bgzf;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::tabix;

/// Collects the message of every `warn!` event raised while it is the
/// thread's default subscriber (`tracing::subscriber::with_default`), so a
/// test reads the warnings a plugin writes at load or match time; the second
/// list holds, for each event, the input record it names through the field
/// `input_record` (0 for none).
pub struct Warnings(pub Mutex<Vec<String>>, pub Mutex<Vec<u64>>);

impl Default for Warnings {
    /// Installs [`Interested`] as the test binary's global default subscriber
    /// before the first `Warnings` exists, so no thread of the binary ever
    /// runs under `tracing`'s no-op default.
    fn default() -> Self {
        static INSTALLED: Once = Once::new();
        INSTALLED.call_once(|| {
            let _ = tracing::subscriber::set_global_default(Interested);
        });
        Warnings(Mutex::new(Vec::new()), Mutex::new(Vec::new()))
    }
}

/// The test binary's global default subscriber: interested in every callsite,
/// enabled for none, with a level hint of off. `tracing` caches a callsite's
/// interest when the callsite is first hit, computed from the hitting thread's
/// subscriber alone while one dispatcher is registered, and its no-op default
/// answers never, which silences that `warn!` for a [`Warnings`] on another
/// thread until the next subscriber registration. Answering sometimes makes
/// every event ask its own thread's subscriber instead, and the hint of off
/// keeps every callsite unregistered until a `Warnings` raises the level.
struct Interested;

impl tracing::Subscriber for Interested {
    fn register_callsite(
        &self,
        _: &'static tracing::Metadata<'static>,
    ) -> tracing::subscriber::Interest {
        tracing::subscriber::Interest::sometimes()
    }

    fn max_level_hint(&self) -> Option<tracing::level_filters::LevelFilter> {
        Some(tracing::level_filters::LevelFilter::OFF)
    }

    fn enabled(&self, _: &tracing::Metadata<'_>) -> bool {
        false
    }

    fn new_span(&self, _: &tracing::span::Attributes<'_>) -> tracing::span::Id {
        tracing::span::Id::from_u64(1)
    }

    fn record(&self, _: &tracing::span::Id, _: &tracing::span::Record<'_>) {}

    fn record_follows_from(&self, _: &tracing::span::Id, _: &tracing::span::Id) {}

    fn event(&self, _: &tracing::Event<'_>) {}

    fn enter(&self, _: &tracing::span::Id) {}

    fn exit(&self, _: &tracing::span::Id) {}
}

impl tracing::Subscriber for Warnings {
    fn enabled(&self, metadata: &tracing::Metadata<'_>) -> bool {
        *metadata.level() == tracing::Level::WARN
    }

    fn new_span(&self, _: &tracing::span::Attributes<'_>) -> tracing::span::Id {
        tracing::span::Id::from_u64(1)
    }

    fn record(&self, _: &tracing::span::Id, _: &tracing::span::Record<'_>) {}

    fn record_follows_from(&self, _: &tracing::span::Id, _: &tracing::span::Id) {}

    fn event(&self, event: &tracing::Event<'_>) {
        struct Message(String, u64);
        impl tracing::field::Visit for Message {
            fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
                if field.name() == "message" {
                    self.0 = format!("{value:?}");
                }
            }
            fn record_u64(&mut self, field: &tracing::field::Field, value: u64) {
                if field.name() == vep_core::warning_hold::RECORD_FIELD {
                    self.1 = value;
                }
            }
        }
        let mut message = Message(String::new(), 0);
        event.record(&mut message);
        self.0.lock().unwrap().push(message.0);
        self.1.lock().unwrap().push(message.1);
    }

    fn enter(&self, _: &tracing::span::Id) {}

    fn exit(&self, _: &tracing::span::Id) {}
}

/// Which column layout a fixture's tabix index describes.
///
/// Column numbers are 1-based here, matching the `tabix -s/-b/-e` flags the
/// upstream data-preparation instructions use.
#[derive(Debug, Clone, Copy)]
pub struct IndexSpec {
    /// 1-based sequence-name column (`tabix -s`).
    pub seq_col: usize,
    /// 1-based start-position column (`tabix -b`).
    pub start_col: usize,
    /// 1-based end-position column (`tabix -e`); `None` indexes a record on
    /// its start alone.
    pub end_col: Option<usize>,
    /// Number of leading lines to skip, matching `tabix -S`.
    pub skip_lines: usize,
    /// Comment-line prefix, matching `tabix -c`. `None` writes no prefix, for a
    /// header line without a `#` as the dbscSNV plugin documents one
    /// (dbscSNV.pm:61-62).
    pub comment_prefix: Option<u8>,
    /// Whether the start column is zero-based with an exclusive end, matching
    /// `tabix -0` (the BED convention; the made-up CADD-SV fixtures use it).
    pub zero_based: bool,
}

impl IndexSpec {
    /// The conventional layout: sequence in column 1, position in column 2.
    pub fn standard() -> Self {
        Self {
            seq_col: 1,
            start_col: 2,
            end_col: None,
            skip_lines: 0,
            comment_prefix: Some(b'#'),
            zero_based: false,
        }
    }

    /// Index zero-based starts with exclusive ends, as `tabix -0` does.
    pub fn zero_based(mut self) -> Self {
        self.zero_based = true;
        self
    }

    /// Index on an arbitrary 1-based position column, as the GRCh38 REVEL
    /// preparation does with `-b 3`.
    #[allow(dead_code)]
    pub fn with_start_col(mut self, start_col: usize) -> Self {
        self.start_col = start_col;
        self
    }

    /// Index on arbitrary 1-based sequence, start and end columns, as the LOEUF
    /// preparation does with `-s 76 -b 77 -e 78`.
    pub fn with_span_cols(mut self, seq_col: usize, start_col: usize, end_col: usize) -> Self {
        self.seq_col = seq_col;
        self.start_col = start_col;
        self.end_col = Some(end_col);
        self
    }

    /// Treat the first line as a skipped header carrying no `#` prefix, the
    /// header form the dbscSNV plugin documents (dbscSNV.pm:61-62).
    pub fn bare_header(mut self) -> Self {
        self.skip_lines = 1;
        self.comment_prefix = None;
        self
    }
}

/// The four bases a made-up sequence is drawn from.
pub const BASES: [u8; 4] = *b"ACGT";

/// A made-up mitochondrial reference of `length` bases under the name
/// `MT`, written into `dir` with its index, and its bases: each drawn from
/// a fixed linear congruential sequence, so the file is the same on every
/// run and copies no published sequence.
pub fn made_up_reference(
    dir: &Path,
    length: usize,
) -> (std::sync::Arc<vep_fasta::IndexedFasta>, Vec<u8>) {
    made_up_reference_named(dir, "MT", length)
}

/// [`made_up_reference`] with its one sequence named `name`.
pub fn made_up_reference_named(
    dir: &Path,
    name: &str,
    length: usize,
) -> (std::sync::Arc<vep_fasta::IndexedFasta>, Vec<u8>) {
    let mut state: u64 = 37;
    let bases: Vec<u8> = (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            BASES[(state >> 33) as usize % 4]
        })
        .collect();
    let path = dir.join(format!("reference_{name}_synthetic.fa"));
    let mut text = format!(">{name}\n").into_bytes();
    let offset = text.len();
    text.extend_from_slice(&bases);
    text.push(b'\n');
    std::fs::write(&path, text).unwrap();
    std::fs::write(
        dir.join(format!("reference_{name}_synthetic.fa.fai")),
        format!("{name}\t{length}\t{offset}\t{length}\t{}\n", length + 1),
    )
    .unwrap();
    (
        std::sync::Arc::new(vep_fasta::IndexedFasta::from_path(&path).unwrap()),
        bases,
    )
}

/// The alternate allele a made-up row gives beside `reference`: the base that
/// follows it in [`BASES`], the first after the last.
pub fn made_up_alternate(reference: u8) -> u8 {
    let place = BASES
        .iter()
        .position(|base| *base == reference)
        .expect("a made-up base");
    BASES[(place + 1) % 4]
}

/// One made-up table row for each of `positions`, 1-based, written by `line`
/// from the contig, the position, the row's reference allele and its
/// alternate ([`made_up_alternate`]). The reference allele of the row at a
/// position is the base `bases` holds `offset` positions on: `offset` 0 is a
/// table that counts on the sequence of `bases`, 1 a table whose rows state
/// the bases one position on.
pub fn made_up_rows(
    contig: &str,
    bases: &[u8],
    positions: impl Iterator<Item = usize>,
    offset: usize,
    line: impl Fn(&str, usize, char, char) -> String,
) -> Vec<String> {
    positions
        .map(|position| {
            let reference = bases[position - 1 + offset];
            line(
                contig,
                position,
                reference as char,
                made_up_alternate(reference) as char,
            )
        })
        .collect()
}

/// The first 1-based position from `from` on whose base `bases` repeats one
/// position on: there a row of a table one position on carries the alleles of
/// a record on the sequence of `bases`.
pub fn repeated_base(bases: &[u8], from: usize) -> usize {
    (from..bases.len())
        .find(|position| bases[position - 1] == bases[*position])
        .expect("a made-up sequence repeats a base")
}

/// Write `body` as a bgzip-compressed file at `dir/file_name`, with no index:
/// the shape a reader of a file's preamble alone opens. The body is bytes, so a
/// test can place a byte that is not UTF-8 where a plugin reads one.
pub fn write_bgzf_bytes(dir: &Path, file_name: &str, body: &[u8]) -> PathBuf {
    let data_path = dir.join(file_name);
    let mut writer =
        bgzf::io::Writer::new(std::fs::File::create(&data_path).expect("create fixture"));
    writer.write_all(body).expect("write fixture body");
    writer.try_finish().expect("finish bgzf");
    data_path
}

/// Write `lines` as a real bgzip-compressed, tabix-indexed file.
///
/// Returns the path to the `.gz`; its `.tbi` sits beside it, which is where
/// [`crate::tabix::TabixAnnotator`] looks. `lines` must already be sorted by the
/// spec's position column, as tabix requires. A line is text or bytes, so a
/// fixture can hold a byte that is not UTF-8; its sequence and position cells
/// are text.
pub fn write_tabix_fixture(
    dir: &Path,
    file_name: &str,
    lines: &[impl AsRef<[u8]>],
    spec: IndexSpec,
) -> PathBuf {
    let data_path = dir.join(file_name);

    // The whole body is one bgzf block, as `bgzip` writes a small file. Calling
    // `virtual_position()` between lines forces a block per record and hides the
    // truncation of records that share a block.
    let body: Vec<u8> = lines
        .iter()
        .flat_map(|l| l.as_ref().iter().copied().chain(*b"\n"))
        .collect();
    let mut writer =
        bgzf::io::Writer::new(std::fs::File::create(&data_path).expect("create fixture"));
    writer.write_all(&body).expect("write fixture body");
    writer.try_finish().expect("finish bgzf");

    // One chunk per contig, from its first record's byte offset in the single
    // block to the end of the body: the shape a real tabix bin has.
    let mut byte_offset: usize = 0;
    let mut chunks: Vec<(String, u64, Chunk)> = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        let line = line.as_ref();
        let line_start = byte_offset;
        byte_offset += line.len() + 1; // + '\n'

        if i < spec.skip_lines || line.starts_with(b"#") {
            continue;
        }
        let fields: Vec<&str> = line
            .split(|&byte| byte == b'\t')
            .map(|cell| std::str::from_utf8(cell).unwrap_or(""))
            .collect();
        let seq = fields[spec.seq_col - 1].to_string();
        // The leading integer, as `tabix` reads a start cell, one below zero
        // placed at the first base and digits past the range at the last.
        let cell = fields[spec.start_col - 1];
        let pos: u64 = crate::tabix::leading_integer(cell)
            .or_else(|| crate::tabix::minus_leads_a_number(cell).then_some(0))
            .or_else(|| {
                cell.starts_with(|digit: char| digit.is_ascii_digit())
                    .then_some(u64::MAX)
            })
            .expect("fixture position column must open with an integer");

        // Virtual position: block offset 0 (single block), byte offset within it.
        let start_vpos = bgzf::VirtualPosition::try_from((0, line_start as u16))
            .expect("fixture body must fit one bgzf block");
        let end_vpos = bgzf::VirtualPosition::try_from((0, byte_offset as u16))
            .expect("fixture body must fit one bgzf block");

        match chunks.iter_mut().find(|(s, _, _)| *s == seq) {
            Some((_, first_pos, span)) => {
                *first_pos = (*first_pos).min(pos);
                *span = Chunk::new(span.start(), end_vpos);
            }
            None => chunks.push((seq, pos, Chunk::new(start_vpos, end_vpos))),
        }
    }

    write_index(&data_path, &chunks, spec);
    data_path
}

/// Build and write the `.tbi` for a fixture.
fn write_index(data_path: &Path, chunks: &[(String, u64, Chunk)], spec: IndexSpec) {
    use indexmap::IndexMap;
    use noodles::csi::binning_index::index::header::format::CoordinateSystem;
    use noodles::csi::binning_index::index::header::Format;
    use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
    use noodles::csi::binning_index::index::reference_sequence::Bin;
    use noodles::csi::binning_index::index::ReferenceSequence;
    use std::collections::HashMap;

    // Preserve first-seen sequence order so reference ids are stable.
    let mut seq_names: Vec<String> = Vec::new();
    for (seq, _, _) in chunks {
        if !seq_names.iter().any(|s| s == seq) {
            seq_names.push(seq.clone());
        }
    }

    let mut per_seq: HashMap<String, Vec<(u64, Chunk)>> = HashMap::new();
    for (seq, pos, chunk) in chunks {
        per_seq.entry(seq.clone()).or_default().push((*pos, *chunk));
    }

    let reference_sequences: Vec<ReferenceSequence<LinearIndex>> = seq_names
        .iter()
        .map(|seq| {
            let entries = per_seq.get(seq).cloned().unwrap_or_default();
            // One bin covering everything: the reader still filters records by the
            // query interval, so a coarse bin exercises the same path as a precise one.
            let mut bins: IndexMap<usize, Bin> = IndexMap::new();
            let all_chunks: Vec<Chunk> = entries.iter().map(|(_, c)| *c).collect();
            let first = all_chunks.first().copied();
            bins.insert(0, Bin::new(all_chunks));

            let linear_index = first.map(|c| vec![c.start()]).unwrap_or_default();
            ReferenceSequence::new(bins, linear_index, None)
        })
        .collect();

    let coordinates = if spec.zero_based {
        CoordinateSystem::Bed
    } else {
        CoordinateSystem::Gff
    };
    let header = noodles::csi::binning_index::index::header::Builder::default()
        .set_format(Format::Generic(coordinates))
        .set_reference_sequence_name_index(spec.seq_col - 1)
        .set_start_position_index(spec.start_col - 1)
        .set_end_position_index(spec.end_col.map(|col| col - 1))
        .set_line_comment_prefix(spec.comment_prefix.unwrap_or(b'#'))
        .set_line_skip_count(spec.skip_lines as u32)
        .set_reference_sequence_names(seq_names.iter().map(|s| s.as_bytes().into()).collect())
        .build();

    let index = tabix::Index::builder()
        .set_header(header)
        .set_reference_sequences(reference_sequences)
        .build();

    let index_path = PathBuf::from(format!("{}.tbi", data_path.display()));
    let mut index_writer =
        tabix::io::Writer::new(std::fs::File::create(&index_path).expect("create tbi"));
    index_writer.write_index(&index).expect("write tbi");
}

mod tests {
    use super::Warnings;

    fn warn_once() {
        tracing::warn!("one shared callsite");
    }

    /// A `warn!` whose callsite another thread reached first with no subscriber
    /// of its own is still collected by this thread's `Warnings`.
    #[test]
    fn a_warning_first_raised_without_a_subscriber_is_still_collected() {
        let warnings = std::sync::Arc::new(Warnings::default());
        tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            std::thread::spawn(warn_once).join().unwrap();
            warn_once();
        });
        assert_eq!(*warnings.0.lock().unwrap(), ["one shared callsite"]);
    }
}
