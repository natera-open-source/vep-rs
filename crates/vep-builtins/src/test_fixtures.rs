// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Real bgzip+tabix fixture builder for plugin tests (test builds only).
//!
//! These helpers write genuine bgzip-compressed, tabix-indexed files so tests
//! exercise the bgzf + CSI read path in [`crate::tabix`], the path every tabix
//! plugin runs through. Two defect classes live only there and are invisible to
//! a test built on in-memory strings:
//!
//! - a data file whose header lacks a leading `#` leaves `header = None`, so
//!   every score lookup by name returns `None` and the plugin annotates
//!   nothing while exiting 0 (dbscSNV ships such headers on both assemblies);
//! - a file carrying two coordinate columns read on the first one discards
//!   every record when indexed on the second (REVEL on GRCh38).
//!
//! Fixtures are tiny (a handful of rows), so they cost milliseconds and need no
//! network or external `bgzip`/`tabix` binaries.

use std::io::Write;
use std::path::{Path, PathBuf};

use noodles::bgzf;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::tabix;

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
    /// Number of leading lines to skip, matching `tabix -S`.
    pub skip_lines: usize,
    /// Comment-line prefix, matching `tabix -c`. `None` writes no prefix, which
    /// is how dbscSNV ships (its header line has no `#`).
    pub comment_prefix: Option<u8>,
}

impl IndexSpec {
    /// The conventional layout: sequence in column 1, position in column 2.
    pub fn standard() -> Self {
        Self {
            seq_col: 1,
            start_col: 2,
            skip_lines: 0,
            comment_prefix: Some(b'#'),
        }
    }

    /// Index on an arbitrary 1-based position column, as the GRCh38 REVEL
    /// preparation does with `-b 3`.
    #[allow(dead_code)]
    pub fn with_start_col(mut self, start_col: usize) -> Self {
        self.start_col = start_col;
        self
    }

    /// Treat the first line as a skipped header carrying no `#` prefix, which is
    /// exactly how dbscSNV1.1 is distributed.
    pub fn bare_header(mut self) -> Self {
        self.skip_lines = 1;
        self.comment_prefix = None;
        self
    }
}

/// Write `lines` as a real bgzip-compressed, tabix-indexed file.
///
/// Returns the path to the `.gz`; its `.tbi` sits beside it, which is where
/// [`crate::tabix::TabixAnnotator`] looks. `lines` must already be sorted by the
/// spec's position column, as tabix requires.
pub fn write_tabix_fixture(
    dir: &Path,
    file_name: &str,
    lines: &[&str],
    spec: IndexSpec,
) -> PathBuf {
    let data_path = dir.join(file_name);

    // The whole body is one bgzf block, as `bgzip` writes a small file. Calling
    // `virtual_position()` between lines forces a block per record and hides the
    // truncation of records that share a block.
    let body: String = lines.iter().map(|l| format!("{l}\n")).collect();
    let mut writer =
        bgzf::io::Writer::new(std::fs::File::create(&data_path).expect("create fixture"));
    write!(writer, "{body}").expect("write fixture body");
    writer.try_finish().expect("finish bgzf");

    // One chunk per contig, from its first record's byte offset in the single
    // block to the end of the body: the shape a real tabix bin has.
    let mut byte_offset: usize = 0;
    let mut chunks: Vec<(String, u64, Chunk)> = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        let line_start = byte_offset;
        byte_offset += line.len() + 1; // + '\n'

        if i < spec.skip_lines || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let seq = fields[spec.seq_col - 1].to_string();
        let pos: u64 = fields[spec.start_col - 1]
            .trim()
            .parse()
            .expect("fixture position column must be numeric");

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

    let header = noodles::csi::binning_index::index::header::Builder::default()
        .set_reference_sequence_name_index(spec.seq_col - 1)
        .set_start_position_index(spec.start_col - 1)
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
