// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The ordering precondition of the module documentation of [`super`]: the
//! monitor that reads every row's record against the rows before it, the
//! refusal it ends in, and the warning of a run that skips it.

use std::collections::HashSet;
use std::fmt;

use arrow_array::cast::AsArray;
use arrow_array::types::Int32Type;

use super::rows::RowBatch;

/// Whether the writer checks the order of its input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputOrder {
    /// Every row passes the [`OrderMonitor`], and the file is laid out sorted
    /// by (`contig_index`, `pos`).
    Checked {
        /// The `contig_index` of every chromosome outside the footer's
        /// `contig_order` ([`super::ParquetOutputFormatter::unlisted_rank`]),
        /// whose records are told apart by their spelling.
        unlisted_rank: i32,
    },
    /// Nothing is checked: rows are encoded in input order and no row group
    /// declares a sort.
    Unchecked,
}

/// The warning a run that does not check its input order writes once.
pub const UNCHECKED_ORDER_WARNING: &str = "input order not checked: Parquet row groups are not sorted and parts are not rolled at contig boundaries";

/// A record out of chromosome and position order: the first row that broke
/// the precondition, as the run refuses it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnsortedInput {
    /// The record's number in the input, as the skip warnings count records.
    pub line: u64,
    /// The record's chromosome as the input spells it.
    pub chrom: String,
    /// The record's written position.
    pub position: u64,
    /// The chromosome of the record it should have followed.
    pub previous_chrom: String,
    /// The greatest written position seen on that chromosome.
    pub previous_position: u64,
}

impl fmt::Display for UnsortedInput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "input is not sorted by chromosome and position at line {} ({}:{} after {}:{}); sort the input or pass --no_check_variants_order",
            self.line, self.chrom, self.position, self.previous_chrom, self.previous_position
        )
    }
}

impl std::error::Error for UnsortedInput {}

/// The chromosome the rows are on and how far along it they have come.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Current {
    /// The chromosome's `contig_index`.
    pub rank: i32,
    /// The chromosome as the input spelled it when the run began.
    pub chrom: String,
    /// The greatest written position seen on the chromosome.
    pub max_position: u64,
}

/// Reads rows in input order and refuses the first that is out of order: a
/// chromosome that returns after another, or a record whose written position
/// is below the greatest seen on its chromosome. A chromosome is its
/// `contig_index`, so the spellings and synonyms that resolve to one cache
/// chromosome are one run of records; the chromosomes outside `contig_order`
/// share one rank and are told apart by their spelling.
#[derive(Debug)]
pub struct OrderMonitor {
    unlisted_rank: i32,
    current: Option<Current>,
    /// The chromosomes whose runs of records have ended.
    finished: HashSet<(i32, String)>,
}

impl OrderMonitor {
    pub fn new(unlisted_rank: i32) -> Self {
        Self {
            unlisted_rank,
            current: None,
            finished: HashSet::new(),
        }
    }

    /// The key a row's chromosome is told apart by: its rank, and its spelling
    /// for a chromosome outside `contig_order`.
    fn identity(&self, rank: i32, chrom: &str) -> (i32, String) {
        let spelling = if rank == self.unlisted_rank {
            chrom.to_string()
        } else {
            String::new()
        };
        (rank, spelling)
    }

    /// The chromosome of the last row read and its greatest position; `None`
    /// before the first row.
    pub fn current(&self) -> Option<&Current> {
        self.current.as_ref()
    }

    /// Reads a batch's rows in order.
    pub fn observe(&mut self, batch: &RowBatch) -> Result<(), UnsortedInput> {
        let ranks = batch.rows.column(0).as_primitive::<Int32Type>().values();
        let chroms = batch.rows.column(1).as_string::<i32>();
        for (i, (&rank, &position)) in ranks.iter().zip(&batch.positions).enumerate() {
            let chrom = chroms.value(i);
            let line = batch.lines[i];
            let current_run = self.current.as_mut().filter(|current| {
                current.rank == rank && (rank != self.unlisted_rank || current.chrom == chrom)
            });
            if let Some(current) = current_run {
                if position < current.max_position {
                    return Err(self.refusal(line, chrom, position));
                }
                current.max_position = position;
                continue;
            }
            if self.finished.contains(&self.identity(rank, chrom)) {
                return Err(self.refusal(line, chrom, position));
            }
            if let Some(previous) = self.current.take() {
                let identity = self.identity(previous.rank, &previous.chrom);
                self.finished.insert(identity);
            }
            self.current = Some(Current {
                rank,
                chrom: chrom.to_string(),
                max_position: position,
            });
        }
        Ok(())
    }

    /// The refusal of the record at `line`, written at `chrom:position`, against
    /// the current run.
    fn refusal(&self, line: u64, chrom: &str, position: u64) -> UnsortedInput {
        let current = self.current.as_ref().expect("a row was read before");
        UnsortedInput {
            line,
            chrom: chrom.to_string(),
            position,
            previous_chrom: current.chrom.clone(),
            previous_position: current.max_position,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::fields::FieldOptions;
    use crate::output::parquet::{ParquetOutputFormatter, RowBatchBuilder};
    use std::sync::Arc;
    use vep_core::consequence::Consequence;
    use vep_core::variant::InputVariant;

    /// A batch of one intergenic record per `(chrom, position, line)`, the
    /// chromosomes ranked `1`, `2`, `10` in that order.
    fn batch(records: &[(&str, u64, u64)]) -> RowBatch {
        let formatter = Arc::new(ParquetOutputFormatter::new(
            FieldOptions::default(),
            vec![],
            vec!["1".into(), "2".into(), "10".into()],
        ));
        let mut builder = RowBatchBuilder::new(formatter);
        for &(chrom, position, line) in records {
            let mut v = InputVariant::new(
                chrom.into(),
                position,
                position,
                b"A".to_vec(),
                b"G".to_vec(),
            );
            v.most_severe_consequence = Some(Consequence::IntergenicVariant);
            v.input_position = position;
            v.input_record = line;
            builder.push_record(&[&v], &mut Vec::new());
        }
        builder.finish()
    }

    /// The ranks of the test formatter: `1`, `2`, `10`, then every other
    /// chromosome at 3.
    const UNLISTED: i32 = 3;

    /// A fresh monitor's reading of `records` as one batch.
    fn observe(records: &[(&str, u64, u64)]) -> Result<(), UnsortedInput> {
        OrderMonitor::new(UNLISTED).observe(&batch(records))
    }

    #[test]
    fn positions_may_repeat_but_not_decrease_within_a_chromosome() {
        assert_eq!(
            observe(&[("21", 100, 1), ("21", 100, 2), ("21", 150, 3)]),
            Ok(())
        );
        let error = observe(&[("21", 100, 1), ("21", 50, 2)]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "input is not sorted by chromosome and position at line 2 (21:50 after 21:100); sort the input or pass --no_check_variants_order"
        );
        assert_eq!(
            observe(&[("21", 100, 1), ("21", 200, 2), ("21", 150, 3)])
                .unwrap_err()
                .to_string(),
            "input is not sorted by chromosome and position at line 3 (21:150 after 21:200); sort the input or pass --no_check_variants_order"
        );
    }

    /// Chromosomes may come in any order and reset the position, as long as
    /// each one's records stay together; the refusal names the record that
    /// returned to an earlier chromosome and the last record of the run it
    /// broke.
    #[test]
    fn a_chromosome_may_not_return_after_another() {
        assert_eq!(observe(&[("1", 500, 1), ("10", 5, 2), ("2", 7, 3)]), Ok(()));
        assert_eq!(observe(&[("2", 500, 1), ("1", 5, 2)]), Ok(()));
        let error = observe(&[("1", 500, 1), ("2", 5, 2), ("1", 600, 3)]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "input is not sorted by chromosome and position at line 3 (1:600 after 2:5); sort the input or pass --no_check_variants_order"
        );
    }

    /// Spellings of one ranked chromosome (`2`, `chr2`) are one run, so the
    /// positions must keep rising across them and a spelling may return; two
    /// chromosomes outside `contig_order` share the trailing rank and are
    /// told apart by their spelling, each with positions of its own.
    #[test]
    fn a_run_is_a_ranked_chromosome_or_an_unlisted_spelling() {
        assert_eq!(
            observe(&[("2", 100, 1), ("chr2", 150, 2), ("2", 200, 3)]),
            Ok(())
        );
        let error = observe(&[("2", 100, 1), ("chr2", 50, 2)]).unwrap_err();
        assert_eq!(error.line, 2);
        assert_eq!((error.chrom.as_str(), error.position), ("chr2", 50));
        assert_eq!(
            (error.previous_chrom.as_str(), error.previous_position),
            ("2", 100)
        );
        assert_eq!(observe(&[("Z", 1000, 1), ("Q", 5, 2), ("Q", 6, 3)]), Ok(()));
        let error = observe(&[("Z", 1000, 1), ("Q", 5, 2), ("Z", 2000, 3)]).unwrap_err();
        assert_eq!(error.line, 3);
        assert_eq!(
            (error.previous_chrom.as_str(), error.previous_position),
            ("Q", 5)
        );
    }

    /// The monitor's state carries across batches, and its current chromosome
    /// reports the greatest position read.
    #[test]
    fn the_check_spans_batches() {
        let mut monitor = OrderMonitor::new(UNLISTED);
        monitor
            .observe(&batch(&[("21", 100, 1), ("21", 120, 2)]))
            .unwrap();
        assert_eq!(
            monitor.current(),
            Some(&Current {
                rank: 3,
                chrom: "21".into(),
                max_position: 120,
            })
        );
        let error = monitor.observe(&batch(&[("21", 119, 3)])).unwrap_err();
        assert_eq!(error.line, 3);
        assert_eq!(error.previous_position, 120);
    }
}
