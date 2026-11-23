// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A SQLite database file read as an Ensembl VEP plugin reads one through
//! Perl's DBI.
//!
//! A plugin module that opens `dbi:SQLite:dbname=<file>` prepares a statement,
//! binds its placeholders and takes the first row as a hash of column names
//! (`fetchrow_hashref`). [`SqliteFile::first_row`] is that read: the bound
//! values are text, as `DBD::SQLite` binds a value given without a type, so
//! SQLite applies each compared column's affinity to them, and the row's
//! values are [`PerlValue`]s, the scalars the driver hands Perl. The file is
//! opened read-only and never created.

use std::path::Path;
use std::sync::{Mutex, PoisonError};

use rusqlite::types::ValueRef;
use rusqlite::{Connection, OpenFlags};

use crate::perl_format::{numify_arithmetic, perl_number};

/// One column value as `DBD::SQLite` hands it to Perl: `undef` for NULL, an
/// integer, a floating-point number, or text. A blob reads as its bytes'
/// text.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum PerlValue {
    Undef,
    Integer(i64),
    Real(f64),
    Text(String),
}

impl PerlValue {
    /// The value in string context: `undef` is empty, a floating-point number
    /// prints with fifteen significant digits, text stays as stored.
    pub(crate) fn text(&self) -> String {
        match self {
            Self::Undef => String::new(),
            Self::Integer(value) => value.to_string(),
            Self::Real(value) => perl_number(*value),
            Self::Text(value) => value.clone(),
        }
    }

    /// The value in numeric context: `undef` is 0, and text is its leading
    /// number, 0 when none leads it.
    pub(crate) fn number(&self) -> f64 {
        match self {
            Self::Undef => 0.0,
            Self::Integer(value) => *value as f64,
            Self::Real(value) => *value,
            Self::Text(value) => numify_arithmetic(value).unwrap_or(0.0),
        }
    }
}

impl From<ValueRef<'_>> for PerlValue {
    fn from(value: ValueRef<'_>) -> Self {
        match value {
            ValueRef::Null => Self::Undef,
            ValueRef::Integer(value) => Self::Integer(value),
            ValueRef::Real(value) => Self::Real(value),
            ValueRef::Text(bytes) | ValueRef::Blob(bytes) => {
                Self::Text(String::from_utf8_lossy(bytes).into_owned())
            }
        }
    }
}

/// The first row of a query, by column name.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PerlRow(Vec<(String, PerlValue)>);

impl PerlRow {
    /// The value under `column`, `undef` for a column the row lacks, as a
    /// lookup in the hash `fetchrow_hashref` returns reads; of two columns
    /// with one name the hash keeps the later.
    pub(crate) fn get(&self, column: &str) -> &PerlValue {
        self.0
            .iter()
            .rev()
            .find(|(name, _)| name == column)
            .map_or(&PerlValue::Undef, |(_, value)| value)
    }
}

/// A SQLite database file open for reading.
pub(crate) struct SqliteFile {
    /// A connection is not `Sync`, and a plugin is queried through `&self`.
    connection: Mutex<Connection>,
    path: String,
}

impl SqliteFile {
    /// Opens the database file at `path` read-only and reads its schema, so a
    /// file that is no SQLite database is an error here and not at the first
    /// query. The error's text names no path; the caller's message does.
    pub(crate) fn open(path: &str) -> Result<Self, String> {
        if !Path::new(path).is_file() {
            return Err("no such file".to_string());
        }
        let flags = OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX;
        let connection = Connection::open_with_flags(path, flags).map_err(|e| e.to_string())?;
        connection
            .query_row("SELECT count(*) FROM sqlite_master", [], |_| Ok(()))
            .map_err(|e| e.to_string())?;
        Ok(Self {
            connection: Mutex::new(connection),
            path: path.to_string(),
        })
    }

    /// The path the file was opened by.
    pub(crate) fn path(&self) -> &str {
        &self.path
    }

    /// The first row `sql` yields with its placeholders bound to `bindings` as
    /// text, `None` when it yields none. An `Err` is the database's own
    /// message for a statement it cannot prepare or step, a table the file
    /// lacks among them.
    pub(crate) fn first_row(
        &self,
        sql: &str,
        bindings: &[&str],
    ) -> Result<Option<PerlRow>, String> {
        let connection = self
            .connection
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        let mut statement = connection.prepare(sql).map_err(|e| e.to_string())?;
        let names: Vec<String> = statement.column_names().map(str::to_string).collect();
        let mut rows = statement
            .query(rusqlite::params_from_iter(bindings))
            .map_err(|e| e.to_string())?;
        let Some(row) = rows.next().map_err(|e| e.to_string())? else {
            return Ok(None);
        };
        let mut values = Vec::with_capacity(names.len());
        for (index, name) in names.into_iter().enumerate() {
            let value = row.get_ref(index).map_err(|e| e.to_string())?;
            values.push((name, PerlValue::from(value)));
        }
        Ok(Some(PerlRow(values)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A database of one table per column affinity, each holding the same
    /// rows, written through the same crate.
    fn database(dir: &Path) -> String {
        let path = dir.join("values.synthetic.sqlite");
        let connection = Connection::open(&path).unwrap();
        connection
            .execute_batch(
                "CREATE TABLE typed (k INTEGER, v REAL);
                 CREATE TABLE untyped (k, v);
                 INSERT INTO typed VALUES (3, 1.5), (4, -2), (5, NULL), (6, 0.1 + 0.2),
                                          (7, 1e21), (8, 'abc'), (9, 100.0);
                 INSERT INTO untyped SELECT * FROM typed;",
            )
            .unwrap();
        path.to_str().unwrap().to_string()
    }

    /// `DBD::SQLite` binds an untyped value as text, so a column with integer
    /// affinity matches the bound `3` and a column with no declared type,
    /// which holds the integer as stored, matches nothing.
    #[test]
    fn a_bound_value_is_text_compared_under_the_column_s_affinity() {
        let dir = tempfile::tempdir().unwrap();
        let file = SqliteFile::open(&database(dir.path())).unwrap();
        let typed = file
            .first_row("SELECT * FROM typed WHERE k = ?", &["3"])
            .unwrap()
            .expect("the integer column matches the text 3");
        assert_eq!(typed.get("v"), &PerlValue::Real(1.5));
        assert_eq!(
            file.first_row("SELECT * FROM untyped WHERE k = ?", &["3"]),
            Ok(None)
        );
        assert_eq!(
            file.first_row("SELECT * FROM typed WHERE k = ?", &["3-4"]),
            Ok(None)
        );
    }

    /// The values of a row print and count as the Perl scalars the driver
    /// returns: a real through fifteen significant digits, NULL as the empty
    /// string and 0, text as stored and as its leading number.
    #[test]
    fn column_values_read_as_perl_scalars() {
        let dir = tempfile::tempdir().unwrap();
        let file = SqliteFile::open(&database(dir.path())).unwrap();
        let value_under = |key: &str| {
            file.first_row("SELECT * FROM typed WHERE k = ?", &[key])
                .unwrap()
                .unwrap()
                .get("v")
                .clone()
        };
        for (key, text, number) in [
            ("3", "1.5", 1.5),
            ("4", "-2", -2.0),
            ("5", "", 0.0),
            ("6", "0.3", 0.1 + 0.2),
            ("7", "1e+21", 1e21),
            ("8", "abc", 0.0),
            ("9", "100", 100.0),
        ] {
            let value = value_under(key);
            assert_eq!(value.text(), text, "k = {key}");
            assert_eq!(value.number(), number, "k = {key}");
        }
        let row = file
            .first_row("SELECT * FROM typed WHERE k = ?", &["3"])
            .unwrap()
            .unwrap();
        assert_eq!(row.get("no_such_column"), &PerlValue::Undef);
    }

    /// A missing file and a file that is no database are refused when opened;
    /// a table the database lacks is the query's error.
    #[test]
    fn a_missing_file_a_non_database_and_a_missing_table_are_errors() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("absent.sqlite");
        assert_eq!(
            SqliteFile::open(missing.to_str().unwrap()).err().as_deref(),
            Some("no such file")
        );
        assert!(!missing.exists(), "opening creates no file");

        let text = dir.path().join("not-a-database.txt");
        std::fs::write(
            &text,
            "a line of text that is not a database header\n".repeat(40),
        )
        .unwrap();
        let refusal = SqliteFile::open(text.to_str().unwrap())
            .err()
            .expect("a text file is no database");
        assert!(refusal.contains("not a database"), "{refusal}");

        let file = SqliteFile::open(&database(dir.path())).unwrap();
        let refusal = file
            .first_row("SELECT * FROM absent WHERE k = ?", &["3"])
            .expect_err("the table is absent");
        assert!(refusal.contains("no such table: absent"), "{refusal}");
    }
}
