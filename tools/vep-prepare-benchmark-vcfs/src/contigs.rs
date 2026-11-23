// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The canonical contigs a run keeps: the human set, or the set a `--canonical-contigs`
//! value names.

use crate::pyfmt::{is_py_space, py_path_str, py_repr, py_splitlines, py_str_cmp};
use std::collections::HashSet;
use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

/// The human canonical contigs, Ensembl-style, for GRCh37 and GRCh38 alike; `M` and `MT`
/// both name the mitochondrial genome.
const HUMAN: [&str; 26] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT", "M",
];

pub struct ContigSet {
    /// Every name once, in code point order, in the bytes it was given in.
    names: Vec<Vec<u8>>,
    /// The names a variant line can carry: a line is decoded with replacement characters,
    /// so no name with a byte outside UTF-8 is among them.
    decodable: HashSet<Vec<u8>>,
}

/// Why a `--canonical-contigs` value gives no set.
#[derive(Debug, PartialEq, Eq)]
pub enum NoSet {
    /// The value names no contig: the line for stderr.
    Empty(String),
    /// The file the value names cannot be read as UTF-8 text.
    Unreadable(String),
}

impl ContigSet {
    fn from_names(names: impl IntoIterator<Item = Vec<u8>>) -> ContigSet {
        let mut names: Vec<Vec<u8>> = names.into_iter().collect();
        names.sort_by(|a, b| py_str_cmp(a, b));
        names.dedup();
        let decodable = names
            .iter()
            .filter(|name| std::str::from_utf8(name).is_ok())
            .cloned()
            .collect();
        ContigSet { names, decodable }
    }

    pub fn human() -> ContigSet {
        ContigSet::from_names(HUMAN.iter().map(|name| name.as_bytes().to_vec()))
    }

    /// The set a `--canonical-contigs` value names: the names in the file at `value` where
    /// one exists, on every line that is not blank and does not start with `#` after its
    /// leading blanks, else the names in `value` itself; names are separated by commas and
    /// by the blanks of `str.isspace()`.
    pub fn named(value: &[u8]) -> Result<ContigSet, NoSet> {
        let path = py_path_str(value);
        let path = Path::new(OsStr::from_bytes(&path));
        let names = if path.is_file() {
            let unreadable =
                |detail: String| NoSet::Unreadable(format!("{}: {detail}", path.display()));
            let bytes = std::fs::read(path).map_err(|e| unreadable(e.to_string()))?;
            let text = String::from_utf8(bytes)
                .map_err(|_| unreadable("the contig list is not UTF-8 text".to_string()))?;
            py_splitlines(&text)
                .into_iter()
                .filter(|line| !line.trim_start_matches(is_py_space).starts_with('#'))
                .flat_map(|line| names_in(line.as_bytes()))
                .collect()
        } else {
            names_in(value)
        };
        if names.is_empty() {
            return Err(NoSet::Empty(format!(
                "error: --canonical-contigs {} names no contig",
                py_repr(value)
            )));
        }
        Ok(ContigSet::from_names(names))
    }

    pub fn contains(&self, name: &[u8]) -> bool {
        self.decodable.contains(name)
    }

    /// The names in code point order.
    pub fn sorted(&self) -> &[Vec<u8>] {
        &self.names
    }
}

/// The names in `text`, which commas and blanks separate. A byte outside UTF-8 is part of
/// the name it stands in.
fn names_in(text: &[u8]) -> Vec<Vec<u8>> {
    let mut names = Vec::new();
    let mut name: Vec<u8> = Vec::new();
    for chunk in text.utf8_chunks() {
        for c in chunk.valid().chars() {
            if c == ',' || is_py_space(c) {
                if !name.is_empty() {
                    names.push(std::mem::take(&mut name));
                }
            } else {
                name.extend_from_slice(c.encode_utf8(&mut [0u8; 4]).as_bytes());
            }
        }
        name.extend_from_slice(chunk.invalid());
    }
    if !name.is_empty() {
        names.push(name);
    }
    names
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sorted(set: &ContigSet) -> Vec<String> {
        set.sorted()
            .iter()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect()
    }

    #[test]
    fn the_human_set_is_the_autosomes_the_sex_chromosomes_and_both_mitochondrial_names() {
        let set = ContigSet::human();
        assert_eq!(set.sorted().len(), 26);
        assert_eq!(
            sorted(&set)[..6],
            ["1", "10", "11", "12", "13", "14"].map(String::from)
        );
        assert_eq!(sorted(&set)[22..], ["M", "MT", "X", "Y"].map(String::from));
        for name in ["1", "22", "X", "Y", "MT", "M"] {
            assert!(set.contains(name.as_bytes()), "{name}");
        }
        for name in ["chr1", "23", "x", "", "GL000191.1"] {
            assert!(!set.contains(name.as_bytes()), "{name}");
        }
    }

    #[test]
    fn a_value_names_its_contigs_between_commas_and_blanks() {
        let set = ContigSet::named(b"2L, 2R\t3L,,X \x1f mitochondrion_genome,2L").unwrap();
        assert_eq!(
            sorted(&set),
            ["2L", "2R", "3L", "X", "mitochondrion_genome"]
        );
        let set = ContigSet::named("chrI\u{a0}chrII\u{3000}chrI".as_bytes()).unwrap();
        assert_eq!(sorted(&set), ["chrI", "chrII"]);
    }

    #[test]
    fn a_name_with_a_byte_outside_utf8_is_listed_and_never_matched() {
        let set = ContigSet::named(b"1,sca\xfffold,\xef\xbf\xbd").unwrap();
        assert_eq!(
            set.sorted(),
            [
                b"1".to_vec(),
                b"sca\xfffold".to_vec(),
                b"\xef\xbf\xbd".to_vec()
            ]
        );
        assert!(set.contains(b"1"));
        assert!(set.contains("\u{fffd}".as_bytes()));
        assert!(!set.contains(b"sca\xfffold"));
        assert!(!set.contains("sca\u{fffd}fold".as_bytes()));
    }

    #[test]
    fn a_value_that_names_no_contig_is_refused_in_reprs_words() {
        for (value, line) in [
            (&b""[..], "error: --canonical-contigs '' names no contig"),
            (b",", "error: --canonical-contigs ',' names no contig"),
            (
                b" ,\t,",
                r"error: --canonical-contigs ' ,\t,' names no contig",
            ),
        ] {
            assert_eq!(
                ContigSet::named(value).err(),
                Some(NoSet::Empty(line.to_string()))
            );
        }
    }

    #[test]
    fn a_file_names_its_contigs_on_the_lines_that_are_not_comments() {
        let dir = std::env::temp_dir().join(format!(
            "vep-prepare-benchmark-vcfs-contigs-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("contigs.txt");
        std::fs::write(
            &file,
            "# the chromosome arms\r\n2L,2R\r\n\n  # 3L is left out\n3R 4\rX\u{2028}#Y\u{85}\t#Z\nlast",
        )
        .unwrap();
        let named = |value: &std::path::Path| ContigSet::named(value.as_os_str().as_bytes());
        let set = named(&file).unwrap();
        assert_eq!(sorted(&set), ["2L", "2R", "3R", "4", "X", "last"]);
        // A trailing separator names the same file.
        let slashed = format!("{}/", file.display());
        assert_eq!(
            sorted(&ContigSet::named(slashed.as_bytes()).unwrap()),
            sorted(&set)
        );

        std::fs::write(&file, "# nothing\n\n   \n").unwrap();
        let refusal = format!(
            "error: --canonical-contigs '{}' names no contig",
            file.display()
        );
        assert_eq!(named(&file).err(), Some(NoSet::Empty(refusal)));

        std::fs::write(&file, b"1\n\xff\n").unwrap();
        assert!(matches!(named(&file), Err(NoSet::Unreadable(_))));

        // A directory is no file: its name is read as a list of one name.
        assert_eq!(
            named(&dir).unwrap().sorted(),
            [dir.as_os_str().as_bytes().to_vec()]
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
