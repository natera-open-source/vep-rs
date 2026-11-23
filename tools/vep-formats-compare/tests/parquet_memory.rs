// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The memory a Parquet comparison holds: the rows both outputs share are never held, so
//! the bytes live at once stay below the bytes of the reference's rows. The one test of
//! this file counts every allocation of its process.

mod compare_common;
mod parquet_writer;

use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::BTreeSet;
use std::sync::atomic::{AtomicUsize, Ordering};

use compare_common::*;
use indexmap::IndexMap;
use parquet_writer::plant_tab_rows;
use vep_formats_compare::cache::THREADS_VARIABLE;
use vep_formats_compare::compare_files::compare_files_through;
use vep_formats_compare::remainder::FrontEnd;
use vep_formats_compare::remainder::STREAM;
use vep_formats_compare::rules::Side;

/// The system's allocator, counting the bytes live at once and their peak.
struct Counting;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

fn grown(by: usize) {
    let live = LIVE.fetch_add(by, Ordering::Relaxed) + by;
    PEAK.fetch_max(live, Ordering::Relaxed);
}

// SAFETY: every call is handed to the system's allocator unchanged; the counters hold no
// pointer.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let block = System.alloc(layout);
        if !block.is_null() {
            grown(layout.size());
        }
        block
    }

    unsafe fn dealloc(&self, block: *mut u8, layout: Layout) {
        System.dealloc(block, layout);
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
    }

    unsafe fn realloc(&self, block: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let moved = System.realloc(block, layout, size);
        if !moved.is_null() {
            if size >= layout.size() {
                grown(size - layout.size());
            } else {
                LIVE.fetch_sub(layout.size() - size, Ordering::Relaxed);
            }
        }
        moved
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT";
/// The records of the planted pair and the transcript rows of each.
const RECORDS: usize = 30_000;
const TRANSCRIPTS: usize = 20;

/// 600,000 rows of 30,000 records, the same in a tab reference and in a Parquet candidate
/// of two part files that hold the reference's halves in the other order. The comparison
/// through the stream of record groups, on one worker, holds at no time as many bytes as
/// the reference's rows alone, beyond what the test itself holds: it reads both outputs in
/// blocks and keeps one digest a record group.
#[test]
fn a_parquet_comparison_holds_no_row_both_outputs_share() {
    let directory = scratch("parquet_memory");
    let rows: Vec<String> = (0..RECORDS)
        .flat_map(|record| {
            (0..TRANSCRIPTS).map(move |transcript| {
                format!(
                    "rec{record}\t21:{}\tA\tENSG00000000001\tENST{transcript:011}\tTranscript\t\
                     intron_variant,non_coding_transcript_variant\tMODIFIER",
                    1000 + record
                )
            })
        })
        .collect();
    let mut text = format!("{COLUMNS}\n");
    for row in &rows {
        text.push_str(row);
        text.push('\n');
    }
    let row_bytes = text.len() - COLUMNS.len() - 1;
    let expected = planted(&directory, "tab.txt", &text);
    drop(text);
    const HALF: usize = RECORDS * TRANSCRIPTS / 2;
    let out = directory.join("out.parquet");
    std::fs::create_dir_all(&out).expect("the output's directory is made");
    for (at, half) in [&rows[HALF..], &rows[..HALF]].into_iter().enumerate() {
        let half: Vec<&str> = half.iter().map(String::as_str).collect();
        plant_tab_rows(&out.join(format!("part-{at:05}.parquet")), COLUMNS, &half);
    }
    drop(rows);
    std::env::set_var(THREADS_VARIABLE, "1");
    let before = LIVE.load(Ordering::Relaxed);
    PEAK.store(before, Ordering::Relaxed);
    let mut front_end = FrontEnd::new(STREAM);
    let report = compare_files_through(
        Some(&mut front_end),
        "parquet",
        &expected,
        &out,
        None,
        &mut no_documented(),
        &BTreeSet::new(),
        &IndexMap::new(),
        &IndexMap::new(),
        4,
        false,
        None,
        None,
        None,
        &mut Side::default(),
    )
    .expect("the outputs compare");
    let held = PEAK.load(Ordering::Relaxed) - before;
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(
        number(&report, &["counts", "compared"]),
        (RECORDS * TRANSCRIPTS) as i64
    );
    assert!(
        held < row_bytes,
        "the comparison held {held} bytes at once, the reference's rows being {row_bytes}; {}",
        front_end.note
    );
}
