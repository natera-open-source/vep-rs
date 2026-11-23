// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The join of two digest files on the transcript id, its four id lists, its summary
//! and the row of the parity table it extends.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Write;
use std::os::unix::ffi::OsStringExt;
use std::path::PathBuf;

use crate::csv::{write_record, Table};
use crate::py::{
    grouped, is_decimal, os_text, path_join, push_json_str, repr_str, universal_newlines,
};

/// The name in brackets that opens every message of the join.
const TAG: &str = "[compare_storable_cache_versions]";

/// What the command line asks for. Paths are normalised as `pathlib` writes them.
#[derive(Debug, Clone, Default)]
pub struct Request {
    pub version_a: String,
    pub dump_a: Vec<u8>,
    pub version_b: String,
    pub dump_b: Vec<u8>,
    pub out_dir: Vec<u8>,
    pub parity_csv: Option<Vec<u8>>,
    pub sweep_id: Option<String>,
    pub assembly: Option<String>,
}

/// How a run ends: the exit code, with a line for standard error when it failed.
pub struct Failure {
    pub code: i32,
    pub line: String,
}

impl Failure {
    /// A failure CPython reports itself, in one line of the tool's own.
    fn refused(detail: impl std::fmt::Display) -> Failure {
        Failure {
            code: 1,
            line: format!("ERROR: [vep-storable-digest-join] {detail}"),
        }
    }
}

fn os_path(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

fn shown(path: &[u8]) -> String {
    String::from_utf8_lossy(path).into_owned()
}

/// The four columns behind the id: contig, start, end, `md5_full`, `md5_content`.
type Columns = [String; 5];

/// One digest file as a map from id to its columns; a later line of an id replaces an
/// earlier one.
fn load(path: &[u8]) -> Result<HashMap<String, Columns>, Failure> {
    let bytes = std::fs::read(os_path(path))
        .map_err(|error| Failure::refused(format!("{}: {}", shown(path), os_text(&error))))?;
    let text = String::from_utf8(bytes)
        .map_err(|_| Failure::refused(format!("{}: the file is not UTF-8", shown(path))))?;
    let text = universal_newlines(&text);
    let mut ids = HashMap::new();
    for line in text.split_inclusive('\n') {
        let columns: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        let [id, contig, start, end, full, content] = columns.as_slice() else {
            return Err(Failure {
                code: 1,
                line: format!(
                    "ERROR: {TAG} malformed dump line in {}: {}",
                    shown(path),
                    repr_str(line)
                ),
            });
        };
        ids.insert(
            (*id).to_string(),
            [*contig, *start, *end, *full, *content].map(str::to_string),
        );
    }
    Ok(ids)
}

/// `re.fullmatch(r"sweep-\d{8}T\d{6}Z", id)`: the digits are decimal digits of any script.
fn is_sweep_id(id: &str) -> bool {
    let Some(rest) = id.strip_prefix("sweep-") else {
        return false;
    };
    let chars: Vec<char> = rest.chars().collect();
    chars.len() == 16
        && chars[..8].iter().all(|c| is_decimal(*c))
        && chars[8] == 'T'
        && chars[9..15].iter().all(|c| is_decimal(*c))
        && chars[15] == 'Z'
}

/// A JSON object whose keys keep the place of their first writing, as a Python dict's do.
#[derive(Default)]
struct Summary {
    entries: Vec<(String, String)>,
}

impl Summary {
    /// Sets `key` to a value given as its JSON text.
    fn set(&mut self, key: String, json: String) {
        match self.entries.iter_mut().find(|(name, _)| *name == key) {
            Some(entry) => entry.1 = json,
            None => self.entries.push((key, json)),
        }
    }

    fn string(&mut self, key: &str, value: &str) {
        let mut json = String::new();
        push_json_str(&mut json, value);
        self.set(key.to_string(), json);
    }

    fn count(&mut self, key: String, value: usize) {
        self.set(key, value.to_string());
    }

    /// `json.dumps(summary, indent=2)`.
    fn dumps(&self) -> String {
        let mut out = String::from("{");
        for (index, (key, json)) in self.entries.iter().enumerate() {
            out.push_str(if index > 0 { ",\n  " } else { "\n  " });
            push_json_str(&mut out, key);
            out.push_str(": ");
            out.push_str(json);
        }
        out.push_str(if self.entries.is_empty() { "}" } else { "\n}" });
        out
    }
}

fn write_file(path: &[u8], text: &str) -> Result<(), Failure> {
    std::fs::write(os_path(path), text)
        .map_err(|error| Failure::refused(format!("{}: {}", shown(path), os_text(&error))))
}

/// The path with its last suffix replaced by `.csv.tmp`, as
/// `Path.with_suffix(".csv.tmp")` gives it; `None` for a path without a name.
fn with_tmp_suffix(path: &[u8]) -> Option<Vec<u8>> {
    let split = path
        .iter()
        .rposition(|byte| *byte == b'/')
        .map_or(0, |at| at + 1);
    let (parent, name) = path.split_at(split);
    if name.is_empty() || name == b"." {
        return None;
    }
    // The stem is the name up to its last dot, unless nothing but dots stand before it.
    let stem = match name.iter().rposition(|byte| *byte == b'.') {
        Some(at) if name[..at].iter().any(|byte| *byte != b'.') => &name[..at],
        _ => name,
    };
    Some([parent, stem, b".csv.tmp"].concat())
}

/// Runs the join. `stdout` takes the summary line.
pub fn run(request: &Request, stdout: &mut dyn Write) -> Result<i32, Failure> {
    if let Some(id) = request.sweep_id.as_deref().filter(|id| !id.is_empty()) {
        if !is_sweep_id(id) {
            return Err(Failure {
                code: 2,
                line: format!("ERROR: {TAG} sweep id must be sweep-<UTC>: {id}"),
            });
        }
    }
    let assembly = request.assembly.as_deref().filter(|name| !name.is_empty());
    if request.parity_csv.is_some() != assembly.is_some() {
        return Err(Failure {
            code: 2,
            line: format!("ERROR: {TAG} --parity-csv and --assembly go together"),
        });
    }

    // The two files are read side by side; a failure of the first is the one reported.
    let (a, b) = std::thread::scope(|scope| {
        let b = scope.spawn(|| load(&request.dump_b));
        let a = load(&request.dump_a);
        (a, b.join().expect("the reader of the second file"))
    });
    let (a, b) = (a?, b?);

    let mut only_a: Vec<&str> = a
        .keys()
        .filter(|id| !b.contains_key(*id))
        .map(String::as_str)
        .collect();
    let mut only_b: Vec<&str> = b
        .keys()
        .filter(|id| !a.contains_key(*id))
        .map(String::as_str)
        .collect();
    let mut both: Vec<&str> = a
        .keys()
        .filter(|id| b.contains_key(*id))
        .map(String::as_str)
        .collect();
    only_a.sort_unstable();
    only_b.sort_unstable();
    both.sort_unstable();
    let full = |id: &str| (&a[id][3], &b[id][3]);
    let content = |id: &str| (&a[id][4], &b[id][4]);
    let same = |pair: (&String, &String)| pair.0 == pair.1;
    let identical = both.iter().filter(|id| same(full(id))).count();
    let content_only: Vec<&str> = both
        .iter()
        .copied()
        .filter(|id| !same(full(id)) && same(content(id)))
        .collect();
    let differing: Vec<&str> = both
        .iter()
        .copied()
        .filter(|id| !same(content(id)))
        .collect();

    let (va, vb) = (request.version_a.as_str(), request.version_b.as_str());
    std::fs::create_dir_all(os_path(&request.out_dir)).map_err(|error| {
        Failure::refused(format!("{}: {}", shown(&request.out_dir), os_text(&error)))
    })?;
    let lists: [(String, &[&str]); 4] = [
        (format!("only_in_{va}"), &only_a),
        (format!("only_in_{vb}"), &only_b),
        (format!("identical_after_strip_{va}_{vb}"), &content_only),
        (format!("differing_{va}_{vb}"), &differing),
    ];
    for (name, ids) in &lists {
        let mut text = String::from("transcript_id\tcontig\tstart\tend\n");
        for id in *ids {
            let columns = a.get(*id).unwrap_or_else(|| &b[*id]);
            writeln!(text, "{id}\t{}\t{}\t{}", columns[0], columns[1], columns[2])
                .expect("a string takes text");
        }
        write_file(
            &path_join(&request.out_dir, format!("{name}.tsv").as_bytes()),
            &text,
        )?;
    }

    let mut summary = Summary::default();
    summary.string("version_a", va);
    summary.string("version_b", vb);
    summary.count(format!("storable_{va}_transcripts"), a.len());
    summary.count(format!("storable_{vb}_transcripts"), b.len());
    summary.count("in_both".to_string(), both.len());
    summary.count(format!("only_in_{va}"), only_a.len());
    summary.count(format!("only_in_{vb}"), only_b.len());
    summary.count("identical_full_dump".to_string(), identical);
    summary.count("identical_content_only".to_string(), content_only.len());
    summary.count("differing_content".to_string(), differing.len());
    summary.set(
        "content_digest_normalises".to_string(),
        "[\n    \"dbID\",\n    \"adaptor\",\n    \"created_date\",\n    \"modified_date\",\n    \"_gene_version\",\n    \"integer scalars stored as strings\"\n  ]"
            .to_string(),
    );
    let content_identical = only_a.is_empty() && only_b.is_empty() && differing.is_empty();
    summary.set(
        "content_identical".to_string(),
        content_identical.to_string(),
    );
    write_file(
        &path_join(
            &request.out_dir,
            format!("summary_{va}_vs_{vb}.json").as_bytes(),
        ),
        &(summary.dumps() + "\n"),
    )?;
    let line = format!(
        "{TAG} {va}: {} transcripts; {vb}: {}; both {}; only {va} {}; only {vb} {}; identical {}; identical after strip {}; differing {}\n",
        grouped(a.len()),
        grouped(b.len()),
        grouped(both.len()),
        grouped(only_a.len()),
        grouped(only_b.len()),
        grouped(identical),
        grouped(content_only.len()),
        grouped(differing.len()),
    );
    stdout
        .write_all(line.as_bytes())
        .and_then(|()| stdout.flush())
        .map_err(|error| Failure::refused(format!("standard output: {}", os_text(&error))))?;

    let (Some(parity_csv), Some(assembly)) = (&request.parity_csv, assembly) else {
        return Ok(0);
    };
    let bytes = std::fs::read(os_path(parity_csv))
        .map_err(|error| Failure::refused(format!("{}: {}", shown(parity_csv), os_text(&error))))?;
    let text = String::from_utf8(bytes)
        .map_err(|_| Failure::refused(format!("{}: the file is not UTF-8", shown(parity_csv))))?;
    let mut table = Table::read(&text)
        .map_err(|error| Failure::refused(format!("{}: {error}", shown(parity_csv))))?;
    let mut extra: Vec<(String, String)> = Vec::new();
    let mut put = |key: String, value: String| match extra.iter_mut().find(|(name, _)| *name == key)
    {
        Some(entry) => entry.1 = value,
        None => extra.push((key, value)),
    };
    put(format!("storable_{va}_transcripts"), a.len().to_string());
    put(format!("storable_{vb}_transcripts"), b.len().to_string());
    put(
        format!("storable_{va}_{vb}_identical_transcripts"),
        (identical + content_only.len()).to_string(),
    );
    put(
        format!("storable_{va}_{vb}_differing_transcripts"),
        differing.len().to_string(),
    );
    put(format!("storable_{va}_only"), only_a.len().to_string());
    put(format!("storable_{vb}_only"), only_b.len().to_string());
    if let Some(id) = request.sweep_id.as_deref().filter(|id| !id.is_empty()) {
        put(format!("storable_{va}_{vb}_sweep_id"), id.to_string());
    }
    for (column, _) in &extra {
        if !table.header.contains(column) {
            table.header.push(column.clone());
            for row in &mut table.rows {
                // The row's dict has no such key yet: the writer's default is the empty string.
                row.push(Some(String::new()));
            }
        }
    }
    let mut hits = Vec::new();
    for (index, row) in table.rows.iter().enumerate() {
        let value = table.get(row, "assembly").map_err(|_| {
            Failure::refused(format!(
                "{}: the table has no column 'assembly'",
                shown(parity_csv)
            ))
        })?;
        if value == Some(assembly) {
            hits.push(index);
        }
    }
    let [hit] = hits.as_slice() else {
        return Err(Failure {
            code: 1,
            line: format!(
                "ERROR: {TAG} {} has {} rows for {assembly}",
                shown(parity_csv),
                hits.len()
            ),
        });
    };
    for (column, value) in &extra {
        // A dict holds one value a name: every column of the name takes it.
        for (index, name) in table.header.iter().enumerate() {
            if name == column {
                table.rows[*hit][index] = Some(value.clone());
            }
        }
    }
    let mut out = String::new();
    let names: Vec<Option<&str>> = table
        .header
        .iter()
        .map(|name| Some(name.as_str()))
        .collect();
    write_record(&mut out, &names);
    for row in &table.rows {
        let values: Vec<Option<&str>> = table
            .header
            .iter()
            .map(|name| table.get(row, name).unwrap_or(None))
            .collect();
        write_record(&mut out, &values);
    }
    let tmp = with_tmp_suffix(parity_csv)
        .ok_or_else(|| Failure::refused(format!("{} has an empty name", shown(parity_csv))))?;
    write_file(&tmp, &out)?;
    std::fs::rename(os_path(&tmp), os_path(parity_csv))
        .map_err(|error| Failure::refused(format!("{}: {}", shown(parity_csv), os_text(&error))))?;
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::{is_sweep_id, with_tmp_suffix};

    /// `sweep-<date>T<time>Z` with the two runs of digits given.
    fn id(date: &str, time: &str) -> String {
        format!("sweep-{date}T{time}Z")
    }

    #[test]
    fn a_sweep_id_is_eight_digits_a_letter_six_digits_and_a_letter() {
        assert!(is_sweep_id(&id("12345678", "123456")));
        // `\d` takes a decimal digit of any script.
        assert!(is_sweep_id(&id("1234567\u{0668}", "123456")));
        for other in [
            id("12345678", "12345"),
            id("1234567", "123456"),
            id("12345678", "123456").replace('T', "t"),
            id("12345678", "123456").replace("sweep", "Sweep"),
            id("12345678", "123456") + " ",
            id("12345678", "123456") + "\n",
            // A superscript two is a digit and no decimal digit.
            id("1234567\u{00b2}", "123456"),
            String::new(),
        ] {
            assert!(!is_sweep_id(&other), "{other:?}");
        }
    }

    /// Each expectation is `str(pathlib.Path(p).with_suffix(".csv.tmp"))` of CPython 3.14.
    #[test]
    fn the_temporary_file_replaces_the_last_suffix() {
        let tmp = |path: &str| {
            with_tmp_suffix(path.as_bytes()).map(|path| String::from_utf8(path).expect("UTF-8"))
        };
        assert_eq!(tmp("a/parity.csv").as_deref(), Some("a/parity.csv.tmp"));
        assert_eq!(tmp("parity").as_deref(), Some("parity.csv.tmp"));
        assert_eq!(tmp("x.tar.gz").as_deref(), Some("x.tar.csv.tmp"));
        assert_eq!(tmp(".hidden").as_deref(), Some(".hidden.csv.tmp"));
        assert_eq!(tmp("..x").as_deref(), Some("..x.csv.tmp"));
        assert_eq!(tmp("name.").as_deref(), Some("name.csv.tmp"));
        assert_eq!(tmp("..").as_deref(), Some("...csv.tmp"));
        assert_eq!(tmp("/abs/t.csv").as_deref(), Some("/abs/t.csv.tmp"));
        assert_eq!(tmp("."), None);
        assert_eq!(tmp("/"), None);
    }
}
