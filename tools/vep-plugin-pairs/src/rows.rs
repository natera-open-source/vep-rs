// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows a class of the scorer can read for the remainder, in the scorer's own form:
//! the key, the row's ordinal among the input's data rows on twelve digits, then the row as
//! written, sorted as `sort` sorts the lines under C collation.
//!
//! The scorer reads each side's rows of a whole Location: the rows under one key for the
//! pair in hand, and the rows of the same input record on the Feature for a class that reads
//! a record's other alleles. So both sides' rows are kept by Location.

use crate::extract::{Extracted, Locator, Run};
use crate::input::{lines, Input, Kind};
use crate::order::strip;
use crate::rules::Res;
use crate::score::interrupt;
use crate::slots::Slot;
use memchr::memchr_iter;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Writes the rows of a default-format `input` at one of `locations` to `path` and returns
/// their number. The caller's slot, `own`, is lent while the pass's threads read.
pub fn write(
    run: &Run,
    input: &Input,
    extracted: &Extracted,
    locations: &HashSet<Vec<u8>>,
    path: &Path,
    own: &mut Slot,
) -> Res<u64> {
    let mut file = BufWriter::with_capacity(1 << 20, File::create(path)?);
    if locations.is_empty() || input.kind != Kind::Default {
        return Ok(0);
    }
    let next = AtomicUsize::new(0);
    let lent = own.lent();
    let kept = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..run.threads.min(input.blocks()).max(1))
            .map(|_| {
                scope.spawn(|| -> Res<(Vec<u8>, Vec<usize>)> {
                    let _slot = run.slots.take();
                    let mut locator = Locator::new(run);
                    let (mut kept, mut ends) = (Vec::new(), Vec::new());
                    let (mut block, mut key) = (Vec::new(), Vec::new());
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= input.blocks() || interrupt::requested() {
                            return Ok((kept, ends));
                        }
                        input.read_block(index, &mut block)?;
                        let mut ordinal = extracted.rows_before[index];
                        for row in lines(&block) {
                            if row.is_empty() || row.starts_with(b"#") {
                                continue;
                            }
                            ordinal += 1;
                            let mut tabs = memchr_iter(b'\t', row);
                            let mut tab = || {
                                tabs.next()
                                    .expect("the extraction counted the row's columns")
                            };
                            let at = [tab(), tab(), tab(), tab(), tab()];
                            locator.set(&row[at[0] + 1..at[1]])?;
                            if !locations.contains(&locator.location) {
                                continue;
                            }
                            key.clear();
                            let (allele, feature) =
                                (&row[at[1] + 1..at[2]], &row[at[3] + 1..at[4]]);
                            for column in
                                [locator.location.as_slice(), strip(allele), strip(feature)]
                            {
                                key.extend_from_slice(column);
                                key.push(b'\t');
                            }
                            kept.extend_from_slice(&key);
                            write!(kept, "{:012}\t", ordinal - 1)
                                .expect("a vector takes every write");
                            kept.extend_from_slice(row);
                            ends.push(kept.len());
                        }
                    }
                })
            })
            .collect();
        workers
            .into_iter()
            .map(|worker| worker.join().expect("a row thread panicked"))
            .collect::<Res<Vec<_>>>()
    });
    drop(lent);
    let kept = kept?;
    let mut rows: Vec<&[u8]> = Vec::new();
    for (bytes, ends) in &kept {
        let mut start = 0;
        for &end in ends {
            rows.push(&bytes[start..end]);
            start = end;
        }
    }
    rows.sort_unstable();
    for row in &rows {
        file.write_all(row)?;
        file.write_all(b"\n")?;
    }
    file.flush()?;
    Ok(rows.len() as u64)
}
