// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A two-column score table as Ensembl VEP's LoFtool and pLI plugins read one:
//! each line split on whitespace into a key and a score (`split`, LoFtool.pm:87,
//! pLI.pm:143), the key lower-cased (`lc`, LoFtool.pm:89, pLI.pm:145) and the
//! score read as Perl numifies the cell and written once, at load, as the
//! module's `sprintf` writes it. A lookup lower-cases the symbol or stable id
//! it is given (LoFtool.pm:118, pLI.pm:221, 227), so a later line with the
//! same key in another case replaces an earlier one, as the Perl hash keeps
//! the last assignment.
//!
//! Both modules read a bundled file beside their own module when the argument
//! names none (`$INC{'<Plugin>.pm'}`, LoFtool.pm:74-78, pLI.pm:115-119):
//! [`bundled_file`] is that default over the run's plugins directories.
//!
//! `lc` folds ASCII letters alone on the byte strings these modules read, so
//! the fold here is ASCII too, and a line is split as those bytes: a byte that
//! is not UTF-8 in a key or a cell is a byte of that token, not a refusal of
//! the file.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use tracing::warn;

use crate::perl_format::numify;
use crate::tabix::BoundedLines;
use crate::traits::PluginError;

/// The longest key [`ScoreTable::get`] lower-cases on the stack: HGNC symbols
/// and Ensembl stable ids run to about 30 bytes.
const STACK_KEY_BYTES: usize = 64;

/// A byte Perl's `split` with no pattern splits on (`\s` on a byte string):
/// ASCII whitespace plus the vertical tab, which [`u8::is_ascii_whitespace`]
/// leaves out.
pub(crate) fn perl_space(byte: u8) -> bool {
    byte.is_ascii_whitespace() || byte == 0x0B
}

/// What the module does with a line whose key has no score cell after it:
/// `my ($gene, $score) = split` leaves `$score` undefined there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingCell {
    /// The line is skipped (`next if !defined($score)`, LoFtool.pm:88).
    Skipped,
    /// The key takes `sprintf`'s `0.00` of an undefined value, which the module
    /// prints (pLI.pm:144-145), so here the key takes no score, as for a cell
    /// with nothing numeric in front.
    NoScore,
}

/// Lower-cased keys to scores as the plugin prints them.
#[derive(Debug, Default)]
pub struct ScoreTable {
    scores: HashMap<String, String>,
}

impl ScoreTable {
    /// Reads `path` for the plugin `name` as `my ($key, $score) = split` reads
    /// each line: the key is the first whitespace-separated token, the score
    /// the second, and a line with blanks after its key but no second token has
    /// the empty cell, which the list assignment keeps (a CRLF file's bare
    /// `SOD1` line), where a line ending at its key has none
    /// (`missing_cell`). An entry is keyed by its key lower-cased, its score
    /// read as Perl numifies the cell ([`numify`]) and written by `format`; the
    /// header, whose score is `header_score` (`LoFtool_percentile`,
    /// LoFtool.pm:88; `pLI`, pLI.pm:144), is skipped, and the later of two
    /// lines for one key wins, as the Perl hash assignment keeps the last.
    ///
    /// A missing file is `<name>: <kind> file <path> not found` (LoFtool.pm:80,
    /// pLI.pm:121), a file with no key at all `<name>: No scores read from
    /// <path>` (LoFtool.pm:94, pLI.pm:191). A cell with nothing numeric in
    /// front gives its key no score, taking an earlier line's score away:
    /// Ensembl VEP stores the `0` its `sprintf` makes of such a cell
    /// (LoFtool.pm:89, pLI.pm:145), which LoFtool's truth test drops and pLI
    /// prints. A cell not read whole is reported once, at load, as a warning
    /// naming the key, the cell and what the key takes from the line.
    pub fn load(
        name: &str,
        kind: &str,
        path: &Path,
        header_score: &str,
        missing_cell: MissingCell,
        format: impl Fn(f64) -> String,
    ) -> Result<Self, PluginError> {
        if !path.exists() {
            return Err(PluginError::Init(format!(
                "{name}: {kind} file {} not found",
                path.display()
            )));
        }
        let file = std::fs::File::open(path).map_err(|e| {
            PluginError::Init(format!("{name}: failed to open {}: {e}", path.display()))
        })?;
        let mut scores = HashMap::new();
        let mut keys_read = false;
        let mut lines = BoundedLines::chomped(BufReader::new(file), path);
        while let Some(line) = lines.next_line() {
            let line = line?;
            let start = line
                .iter()
                .position(|byte| !perl_space(*byte))
                .unwrap_or(line.len());
            let line = &line[start..];
            let (key, after_key) = match line.iter().position(|byte| perl_space(*byte)) {
                Some(end) => (&line[..end], Some(&line[end..])),
                None => (line, None),
            };
            if key.is_empty() {
                // pLI keys a blank line's `0.00` under the empty string, which
                // no symbol matches, so the line counts as read and nothing more.
                if missing_cell == MissingCell::NoScore {
                    keys_read = true;
                }
                continue;
            }
            let cell = after_key.map(|rest| {
                rest.split(|byte| perl_space(*byte))
                    .find(|token| !token.is_empty())
                    .unwrap_or(b"")
            });
            let score = match (cell, missing_cell) {
                (Some(cell), _) => cell,
                (None, MissingCell::Skipped) => continue,
                (None, MissingCell::NoScore) => b"",
            };
            if score == header_score.as_bytes() {
                continue;
            }
            keys_read = true;
            let key = String::from_utf8_lossy(key);
            let score = String::from_utf8_lossy(score);
            match numify(&score) {
                Some((value, rest)) => {
                    let text = format(value);
                    if !rest.is_empty() {
                        warn!(
                            "{name}: the {kind} file {} gives {key} the score {score}, which is \
                             not a number; {key} takes its leading number, {text}",
                            path.display()
                        );
                    }
                    scores.insert(key.to_ascii_lowercase(), text);
                }
                None if score.is_empty() => {
                    warn!(
                        "{name}: the {kind} file {} gives {key} no score cell; {key} takes no \
                         score from the line",
                        path.display()
                    );
                    scores.remove(&key.to_ascii_lowercase());
                }
                None => {
                    warn!(
                        "{name}: the {kind} file {} gives {key} the score {score}, which is not \
                         a number; {key} takes no score from it",
                        path.display()
                    );
                    scores.remove(&key.to_ascii_lowercase());
                }
            }
        }
        if !keys_read {
            return Err(PluginError::Init(format!(
                "{name}: No scores read from {}",
                path.display()
            )));
        }
        Ok(Self { scores })
    }

    /// The score of `key`, matched case-insensitively. The key is lower-cased
    /// on the stack, so a lookup on a transcript row allocates nothing; a key
    /// over [`STACK_KEY_BYTES`] bytes is lower-cased on the heap.
    pub fn get(&self, key: &str) -> Option<&str> {
        let mut lowered = [0u8; STACK_KEY_BYTES];
        let score = match lowered.get_mut(..key.len()) {
            Some(buffer) => {
                buffer.copy_from_slice(key.as_bytes());
                buffer.make_ascii_lowercase();
                // Lower-casing ASCII letters keeps the UTF-8 of `key` valid.
                self.scores
                    .get(std::str::from_utf8(buffer).expect("lower-cased UTF-8"))
            }
            None => self.scores.get(&key.to_ascii_lowercase()),
        };
        score.map(String::as_str)
    }

    /// The number of entries.
    pub fn len(&self) -> usize {
        self.scores.len()
    }

    /// Whether the table has no entry.
    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }
}

/// The bundled data file `file` of the plugin `name`, given no file parameter:
/// the copy in the first of `dirs` (the run's plugins directories in search
/// order) that holds one; when none does, the path in the first directory, so
/// the refusal names where it was looked for. With no directory at all the
/// file is refused by its bare name, `<name>: <kind> file <file> not found`,
/// never read from the working directory.
pub fn bundled_file(
    name: &str,
    kind: &str,
    dirs: &[String],
    file: &str,
) -> Result<PathBuf, PluginError> {
    let candidates = || dirs.iter().map(|dir| Path::new(dir).join(file));
    candidates()
        .find(|path| path.exists())
        .or_else(|| candidates().next())
        .ok_or_else(|| PluginError::Init(format!("{name}: {kind} file {file} not found")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perl_format::sprintf_g;
    use crate::test_fixtures::Warnings;

    /// [`ScoreTable::load`] on a file whose header score is `header`, the
    /// scores written as `%g` writes them.
    fn load(
        name: &str,
        kind: &str,
        path: &Path,
        missing_cell: MissingCell,
    ) -> Result<ScoreTable, PluginError> {
        ScoreTable::load(name, kind, path, "header", missing_cell, sprintf_g)
    }

    fn table(text: &[u8]) -> ScoreTable {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(&path, text).unwrap();
        load("Example", "scores", &path, MissingCell::Skipped).unwrap()
    }

    /// `split` reads whitespace-separated tokens, CRLF line ends and a
    /// vertical tab or form feed included, and `lc` keys the hash, so the
    /// later of two lines spelling one symbol in two cases replaces the
    /// earlier (LoFtool.pm:87-89).
    #[test]
    fn keys_are_lower_cased_and_the_later_line_wins() {
        let table = table(
            b"Gene\theader\r\nC2ORF15\t0.013\r\nC2orf15 \t 0.543\r\nSOD1\x0B1\x0B\r\nTP53\x0C0.25\r\n",
        );
        assert_eq!(table.len(), 3);
        assert_eq!(table.get("C2orf15"), Some("0.543"));
        assert_eq!(table.get("c2ORF15"), Some("0.543"));
        assert_eq!(table.get("SOD1"), Some("1"));
        assert_eq!(table.get("sod1"), Some("1"));
        assert_eq!(table.get("TP53"), Some("0.25"));
        assert_eq!(table.get("SOD2"), None);
    }

    /// A long symbol matches, and a non-ASCII symbol folds its ASCII letters
    /// alone, as `lc` folds the byte strings Ensembl VEP's plugins read.
    #[test]
    fn long_and_non_ascii_keys_match() {
        let long = "G".repeat(65);
        let text = format!("{long}\t0.5\nÄBC\t0.25\n");
        let table = table(text.as_bytes());
        assert_eq!(table.get(&long), Some("0.5"));
        assert_eq!(table.get(&long.to_ascii_lowercase()), Some("0.5"));
        assert_eq!(table.get("Äbc"), Some("0.25"));
        assert_eq!(table.get("äbc"), None);
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) is a byte of its token, as
    /// `split` reads the bytes of a line (LoFtool.pm:87, pLI.pm:143): in a
    /// comment line it is part of a cell with nothing numeric in front, in a
    /// third cell it is never read, in a key it makes another key, which no
    /// symbol matches, and after a score it is the text `sprintf` numifies
    /// past (LoFtool.pm:89, pLI.pm:145). The warnings name the cell with the
    /// byte written as U+FFFD.
    #[test]
    fn a_byte_that_is_not_utf8_is_a_byte_of_its_token() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(
            &path,
            b"# caf\xE9 scores\nGene\theader\nAPP\t0.25\tcaf\xE9\nSOD\xE91\t0.5\nTP53\t0.125\xE9\n",
        )
        .unwrap();
        let warnings = std::sync::Arc::new(Warnings::default());
        let table = tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            load("LoFtool", "scores", &path, MissingCell::Skipped).unwrap()
        });
        assert_eq!(table.len(), 3);
        assert_eq!(table.get("APP"), Some("0.25"));
        assert_eq!(table.get("SOD1"), None);
        assert_eq!(table.get("SOD\u{FFFD}1"), Some("0.5"));
        assert_eq!(table.get("TP53"), Some("0.125"));
        assert_eq!(table.get("#"), None);
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                format!(
                    "LoFtool: the scores file {} gives # the score caf\u{FFFD}, which is not a \
                     number; # takes no score from it",
                    path.display()
                ),
                format!(
                    "LoFtool: the scores file {} gives TP53 the score 0.125\u{FFFD}, which is \
                     not a number; TP53 takes its leading number, 0.125",
                    path.display()
                ),
            ]
        );
    }

    /// A line with one token or a cell with no leading number is no entry, and
    /// a file with no key at all is refused naming the file. Ensembl VEP 116.2
    /// formats such a cell with `sprintf` (pLI.pm:145, LoFtool.pm:89), which
    /// makes `0.00` or `0` of a word or, in pLI, of a missing cell; the cell is
    /// skipped here and reported once, at load, in a warning naming the key and
    /// the cell.
    #[test]
    fn lines_without_a_numeric_score_are_no_entry_and_a_word_is_warned_about() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(
            &path,
            b"Gene\theader\nALONE\nSOD1\tNA\nSOD2\t0.5\nSCAF4\tn/a\n",
        )
        .unwrap();
        let warnings = std::sync::Arc::new(Warnings::default());
        let table = tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            load("pLI", "values", &path, MissingCell::NoScore).unwrap()
        });
        assert_eq!(table.len(), 1);
        assert_eq!(table.get("ALONE"), None);
        assert_eq!(table.get("SOD1"), None);
        assert_eq!(table.get("SCAF4"), None);
        assert_eq!(table.get("SOD2"), Some("0.5"));
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                format!(
                    "pLI: the values file {} gives ALONE no score cell; ALONE takes no score \
                     from the line",
                    path.display()
                ),
                format!(
                    "pLI: the values file {} gives SOD1 the score NA, which is not a number; \
                     SOD1 takes no score from it",
                    path.display()
                ),
                format!(
                    "pLI: the values file {} gives SCAF4 the score n/a, which is not a number; \
                     SCAF4 takes no score from it",
                    path.display()
                ),
            ]
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(&path, b"Gene\theader\n").unwrap();
        let error = load("Example", "scores", &path, MissingCell::Skipped).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: Example: No scores read from {}",
                path.display()
            )
        );
    }

    /// A cell is read as the Perl `sprintf` numifies it: its leading number,
    /// where `1.5abc`, `12.5%`, a decimal comma and a dangling exponent read
    /// the number before the text, `inf` and `nan` read by their first letters
    /// in any case, Visual C's `1.#INF` and `1.#IND` as infinity and not a
    /// number where `1.#QNAN` is `1` (Perl keeps its own reading of a `nan`
    /// only when text follows it, and `strtod` reads `1.`), a second sign as
    /// C's `strtod` reads it after Perl's own, and a `0x` or `0b` form as `0`.
    /// A partially read cell is reported with the number the key takes from it.
    #[test]
    fn a_cell_is_read_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(
            &path,
            b"Gene\theader\r\nAPP\t1.5abc\r\nGENEA2\t12.5%\r\nGENEA3\t1,5\r\nGENEA4\t.5\r\n\
              GENEA5\t5.\r\nGENEA6\t1.5e\r\nGENEA7\t1.5E-3x\r\nGENEA8\t1_000\r\n\
              GENEA9\t0x1A\r\nGENEB1\t0abc\r\nGENEB2\t-0abc\r\nGENEB3\t--1\r\n\
              GENEB4\t+-1\r\nGENEB5\tInf\r\nGENEB6\tinfx\r\nGENEB7\t-Inf\r\n\
              GENEB8\t1.#INF\r\nGENEB9\t1.#IND\r\nGENEC1\t1.#QNAN\r\nGENEC2\tnan(123)\r\n\
              GENEC3\t1e999\r\nGENEC4\t+-inf\r\nGENEC5\t1.26E-01\r\n",
        )
        .unwrap();
        let warnings = std::sync::Arc::new(Warnings::default());
        let table = tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            load("LoFtool", "scores", &path, MissingCell::Skipped).unwrap()
        });
        let expected = [
            ("APP", "1.5"),
            ("GENEA2", "12.5"),
            ("GENEA3", "1"),
            ("GENEA4", "0.5"),
            ("GENEA5", "5"),
            ("GENEA6", "1.5"),
            ("GENEA7", "0.0015"),
            ("GENEA8", "1"),
            ("GENEA9", "0"),
            ("GENEB1", "0"),
            ("GENEB2", "-0"),
            ("GENEB3", "1"),
            ("GENEB4", "-1"),
            ("GENEB5", "Inf"),
            ("GENEB6", "Inf"),
            ("GENEB7", "-Inf"),
            ("GENEB8", "Inf"),
            ("GENEB9", "NaN"),
            ("GENEC1", "1"),
            ("GENEC2", "NaN"),
            ("GENEC3", "Inf"),
            ("GENEC4", "-Inf"),
            ("GENEC5", "0.126"),
        ];
        assert_eq!(table.len(), expected.len());
        for (key, score) in expected {
            assert_eq!(table.get(key), Some(score), "{key}");
        }
        let warnings = warnings.0.lock().unwrap();
        assert_eq!(
            warnings[0],
            format!(
                "LoFtool: the scores file {} gives APP the score 1.5abc, which is not a \
                 number; APP takes its leading number, 1.5",
                path.display()
            )
        );
        let partially_read = [
            "APP", "GENEA2", "GENEA3", "GENEA6", "GENEA7", "GENEA8", "GENEA9", "GENEB1", "GENEB2",
            "GENEB6", "GENEC1", "GENEC2",
        ];
        assert_eq!(warnings.len(), partially_read.len());
        for (warning, key) in warnings.iter().zip(partially_read) {
            assert!(
                warning.contains(&format!(" gives {key} the score ")),
                "{warning}"
            );
        }
    }

    /// The later of two lines for one key wins whatever it carries, as the Perl
    /// hash assignment (LoFtool.pm:89, pLI.pm:145) keeps the last: a later cell
    /// with no leading number, which Perl stores as `0`, takes the earlier score
    /// away, and a later number replaces an earlier such cell. Each line's
    /// warning says what the key takes from that line alone.
    #[test]
    fn the_later_line_wins_whether_or_not_it_carries_a_number() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(
            &path,
            b"Gene\theader\r\nCXADR\t1.38E-04\r\nCXADR\tNA\r\nSCAF4\tNA\r\nSCAF4\t0.85\r\n",
        )
        .unwrap();
        let warnings = std::sync::Arc::new(Warnings::default());
        let table = tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            load("LoFtool", "scores", &path, MissingCell::Skipped).unwrap()
        });
        assert_eq!(table.len(), 1);
        assert_eq!(table.get("CXADR"), None);
        assert_eq!(table.get("SCAF4"), Some("0.85"));
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                format!(
                    "LoFtool: the scores file {} gives CXADR the score NA, which is not a \
                     number; CXADR takes no score from it",
                    path.display()
                ),
                format!(
                    "LoFtool: the scores file {} gives SCAF4 the score NA, which is not a \
                     number; SCAF4 takes no score from it",
                    path.display()
                ),
            ]
        );
    }

    /// A file whose every cell lacks a leading number loads with no score: the
    /// Perl hash holds a `0` per key, so `scalar keys %scores` (LoFtool.pm:94,
    /// pLI.pm:191) is not zero and the plugin runs, writing no field. Only a
    /// file with no key at all is refused.
    #[test]
    fn a_file_of_only_non_numeric_cells_loads_with_no_score() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(&path, b"Gene\theader\r\nAPP\tNA\r\nSOD1\tn/a\r\n").unwrap();
        let table = load("LoFtool", "scores", &path, MissingCell::Skipped).unwrap();
        assert!(table.is_empty());
        assert_eq!(table.get("APP"), None);
    }

    /// `my ($gene, $score) = split` reads a bare key's trailing blank as an
    /// empty score cell (a list assignment keeps the trailing empty field) and
    /// a line ending at its key as no cell. The empty cell is a cell with
    /// nothing numeric in front for both modules; the missing cell is skipped
    /// by LoFtool (`:88`) and scores `0.00` in pLI (`:145`), so there the key
    /// takes no score. A blank line counts as a key read for pLI alone.
    #[test]
    fn a_bare_key_follows_the_modules_split() {
        let text = b"Gene\theader\r\nSOD1\t0.5\r\nSOD1\r\nAPP\t0.25\nAPP\nTP53 \n";
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.txt");
        std::fs::write(&path, text).unwrap();
        let warnings = std::sync::Arc::new(Warnings::default());
        let skipped = tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            load("LoFtool", "scores", &path, MissingCell::Skipped).unwrap()
        });
        assert_eq!(skipped.get("SOD1"), None);
        assert_eq!(skipped.get("APP"), Some("0.25"));
        assert_eq!(skipped.get("TP53"), None);
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                format!(
                    "LoFtool: the scores file {} gives SOD1 no score cell; SOD1 takes no score \
                     from the line",
                    path.display()
                ),
                format!(
                    "LoFtool: the scores file {} gives TP53 no score cell; TP53 takes no score \
                     from the line",
                    path.display()
                ),
            ]
        );
        let no_score = load("pLI", "values", &path, MissingCell::NoScore).unwrap();
        assert_eq!(no_score.get("SOD1"), None);
        assert_eq!(no_score.get("APP"), None);
        assert!(no_score.is_empty());

        std::fs::write(&path, b"Gene\theader\nSOD1\nAPP\n\n").unwrap();
        let error = load("LoFtool", "scores", &path, MissingCell::Skipped).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: LoFtool: No scores read from {}",
                path.display()
            )
        );
        let no_score = load("pLI", "values", &path, MissingCell::NoScore).unwrap();
        assert!(no_score.is_empty());
    }

    #[test]
    fn a_missing_file_is_refused_by_name() {
        let error = load(
            "Example",
            "scores",
            Path::new("/no/such/scores.txt"),
            MissingCell::Skipped,
        )
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "plugin initialization failed: Example: scores file /no/such/scores.txt not found"
        );
    }

    /// The bundled file is read from the first plugins directory that holds
    /// it; with none holding it, the first directory names the path looked
    /// for, and with no directory at all the file is refused by its bare name.
    #[test]
    fn bundled_file_is_the_first_directory_holding_it() {
        let dir = tempfile::tempdir().unwrap();
        let first = dir.path().join("first");
        let second = dir.path().join("second");
        std::fs::create_dir_all(&first).unwrap();
        std::fs::create_dir_all(&second).unwrap();
        std::fs::write(second.join("scores.txt"), b"SOD1\t1\n").unwrap();
        let dirs = [first.display().to_string(), second.display().to_string()];
        assert_eq!(
            bundled_file("Example", "scores", &dirs, "scores.txt").unwrap(),
            second.join("scores.txt")
        );
        assert_eq!(
            bundled_file("Example", "scores", &dirs, "other.txt").unwrap(),
            first.join("other.txt")
        );
        assert_eq!(
            bundled_file("Example", "scores", &[], "scores.txt")
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: Example: scores file scores.txt not found"
        );
    }
}
