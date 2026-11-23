// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse Ensembl protein FASTA files to extract transcript -> protein record mapping.
//!
//! A record is held under its header's `transcript:` token, whole, and found
//! for a model's transcript by `record_of`: a dot precedes the identifier's
//! version in a release that writes one (`ENST00000456.3`) and sits inside
//! the identifier in one that writes none (`AT1G01010.1`,
//! `HORVU.MOREX.r3.7HG0737000.1`), and only the model tells the two apart.
//! `ModelNames` holds what the model states, the identifier and version of
//! each transcript, and orders the names a transcript is looked up under.

use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;

use crate::core_dumps::BoundedLines;
use crate::gff3::{is_whole_number, located, version_number, without_byte_order_marks};

/// The longest protein FASTA line read, 16 MiB. Ensembl wraps sequences at 60
/// residues, and the longest human protein written on one line runs to some
/// 35,000, so the cap refuses a file with no line end in it long before the
/// line buffer grows to its size.
const MAX_LINE_BYTES: usize = 16 << 20;

/// One record of a protein FASTA: the peptide, the protein identifier its
/// header opens with as written (`ENSP00000123456.2`), and the version that
/// identifier ends in, the translation's version in the release the file was
/// written from; `None` when the identifier is written without one, as a
/// release whose identifiers carry no version writes it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Peptide {
    pub sequence: String,
    pub identifier: String,
    pub version: Option<u32>,
}

impl Peptide {
    /// The version the header states for the translation the model names
    /// `protein_id`: none when the header's identifier is `protein_id` whole,
    /// its last component then the identifier's own (`AT1G01010.1`, a
    /// translation Ensembl VEP's cache holds without a version); else the
    /// header's.
    pub(crate) fn version_beside(&self, protein_id: Option<&str>) -> Option<u32> {
        if protein_id == Some(self.identifier.as_str()) {
            None
        } else {
            self.version
        }
    }
}

/// The names the transcripts of a model state, each the pair of an identifier
/// as written and a version (`ENST00000456` at 3, `AT1G01010.1` at none), and
/// the order a transcript is looked up in among the names a file states
/// ([`ModelNames::name_of`]).
#[derive(Debug, Default)]
pub(crate) struct ModelNames<'a> {
    versions: HashMap<&'a str, Vec<Option<u32>>>,
}

impl<'a> ModelNames<'a> {
    /// The names of a model whose transcripts state `pairs`, each an
    /// identifier as written with its version.
    pub(crate) fn of(pairs: impl IntoIterator<Item = (&'a str, Option<u32>)>) -> Self {
        let mut versions: HashMap<&'a str, Vec<Option<u32>>> = HashMap::new();
        for (identifier, version) in pairs {
            let stated = versions.entry(identifier).or_default();
            if !stated.contains(&version) {
                stated.push(version);
            }
        }
        Self { versions }
    }

    /// The versions the transcripts written `identifier` state.
    fn versions_of(&self, identifier: &str) -> &[Option<u32>] {
        self.versions.get(identifier).map_or(&[], Vec::as_slice)
    }

    /// The name a file states for the transcript `id` at `version`, the pair
    /// its line or record is held under; `held` says whether the file states
    /// a pair. The pairs are tried in this order, each compared whole:
    ///
    /// 1. the transcript's own, `id` at `version`;
    /// 2. for a transcript with a version: the dotted form `id.version`
    ///    stating no version (`ENST00000456.3`); the dotted form stating that
    ///    version again, a version written twice; `id` stating no version;
    /// 3. for a transcript without one whose identifier ends in a dot and a
    ///    whole number (`T1.3`): the identifier before that dot stating that
    ///    number; `id` stating that number, a version written twice.
    ///
    /// A pair past the first may be another transcript's: the transcript
    /// that states it as its own takes it, then one whose identifier as
    /// written is the pair's, and one that reaches it by the dotted form
    /// last. A pair another transcript takes ahead of this one is passed
    /// over, so `T1` at version 1 and a sibling written `T1.1` each find
    /// their own line or record. A pair of another version of `id`, or one
    /// sharing only its first components, is never tried.
    pub(crate) fn name_of<'n>(
        &self,
        id: &'n str,
        version: Option<u32>,
        held: &dyn Fn(&str, Option<u32>) -> bool,
    ) -> Option<(Cow<'n, str>, Option<u32>)> {
        if held(id, version) {
            return Some((Cow::Borrowed(id), version));
        }
        match version {
            Some(number) => {
                let dotted = format!("{id}.{number}");
                for stated in [None, version] {
                    if held(&dotted, stated) && !self.found_by_its_own(&dotted, stated, held) {
                        return Some((Cow::Owned(dotted), stated));
                    }
                }
                (held(id, None) && !self.states(id, None)).then_some((Cow::Borrowed(id), None))
            }
            None => {
                let (stem, number) = ending_version(id)?;
                let stated = Some(number);
                if held(stem, stated) && !self.found_by_its_own(stem, stated, held) {
                    return Some((Cow::Borrowed(stem), stated));
                }
                (held(id, stated) && !self.states(id, stated))
                    .then_some((Cow::Borrowed(id), stated))
            }
        }
    }

    /// Whether a transcript of the model states `identifier` at `version` as
    /// its own pair.
    fn states(&self, identifier: &str, version: Option<u32>) -> bool {
        self.versions_of(identifier).contains(&version)
    }

    /// Whether a transcript written `identifier` finds the pair `identifier`
    /// at `stated` for itself ([`ModelNames::name_of`]). The search ends: a
    /// transcript without a version asks only about a held pair, which the
    /// transcript stating it answers at once, and a transcript with one asks
    /// only about a longer identifier.
    fn found_by_its_own(
        &self,
        identifier: &str,
        stated: Option<u32>,
        held: &dyn Fn(&str, Option<u32>) -> bool,
    ) -> bool {
        self.versions_of(identifier).iter().any(|&version| {
            self.name_of(identifier, version, held)
                .is_some_and(|(name, found)| name == identifier && found == stated)
        })
    }
}

/// `identifier` split at its last dot when a version follows it, written as
/// the number prints: `T1.3` is `T1` and 3, and `T1.03`, `T1.` and `T1` are
/// not split.
fn ending_version(identifier: &str) -> Option<(&str, u32)> {
    let (stem, digits) = identifier.rsplit_once('.')?;
    let number: u32 = digits.parse().ok()?;
    (number.to_string() == digits).then_some((stem, number))
}

/// The record of the transcript `id` at `version`, with the key it is held
/// under. A record states its `transcript:` token as an identifier and no
/// version, so the token is looked up in the order of
/// [`ModelNames::name_of`]: `id` itself for a transcript without a version;
/// for one with a version `id.version`, unless that token is the record of
/// another transcript of `names`, then `id`. A token that names another
/// version of `id`, or shares only its first components, is another
/// transcript's record.
pub(crate) fn record_of<'a>(
    records: &'a HashMap<String, Peptide>,
    names: &ModelNames,
    id: &str,
    version: Option<u32>,
) -> Option<(&'a str, &'a Peptide)> {
    let held = |token: &str, stated: Option<u32>| stated.is_none() && records.contains_key(token);
    let (token, _) = names.name_of(id, version, &held)?;
    records
        .get_key_value(token.as_ref())
        .map(|(key, record)| (key.as_str(), record))
}

/// Parse a protein FASTA file (plain, gzipped or bgzipped) into a map of each
/// header's `transcript:` token, whole, to its record (`record_of` finds a
/// model's), with the length of a line capped: a line longer than
/// `MAX_LINE_BYTES` fails the parse, naming the file and the line, as does a
/// line that is not UTF-8 and a header whose protein version is above the
/// largest the cache holds (`version_number`). A byte-order mark at the
/// start of a line is taken off before the line is read
/// (`without_byte_order_marks`): a header behind one is a header, where
/// the mark kept would file it with the sequence before it.
///
/// Ensembl protein FASTA headers look like:
/// ```text
/// >ENSP00000123456.2 pep chromosome:GRCh38:1:100:500:1 gene:ENSG00000123 transcript:ENST00000456 ...
/// MKVLWAALLLLAAMYTISVQTPMGKSRCMKDRHGAFEEHG...
/// ```
///
/// The `transcript:` token (`ENST00000456` here) is mapped to the amino acid
/// sequence and the protein identifier with its version.
pub fn parse_protein_fasta(path: &Path) -> Result<HashMap<String, Peptide>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    // A bgzipped file is a run of gzip members, so a single-member decoder
    // would read its first block alone and stop.
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut sequences: HashMap<String, Peptide> = HashMap::new();
    let mut current: Option<(String, String, Option<u32>)> = None;
    let mut current_seq = String::new();

    for (line_index, line) in
        BoundedLines::new(reader, MAX_LINE_BYTES, path.display().to_string()).enumerate()
    {
        let line_number = line_index + 1;
        let line = String::from_utf8(without_byte_order_marks(line?))
            .map_err(|_| anyhow!("FASTA line is not UTF-8"))
            .with_context(|| located(path, line_number))?;
        let line = line.trim();

        if line.starts_with('>') {
            finish_record(&mut sequences, current.take(), &mut current_seq);
            let version =
                extract_protein_version(line).with_context(|| located(path, line_number))?;
            current = extract_transcript_id(line)
                .map(|tid| (tid, protein_identifier(line).to_string(), version));
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }
    finish_record(&mut sequences, current, &mut current_seq);

    Ok(sequences)
}

/// Stores the record read so far under its transcript when it has one and a
/// sequence, and empties the sequence buffer for the next. The sequence is
/// copied at its length, never moved: a moved buffer keeps its spare capacity
/// in every stored record and leaves the next record to grow a new one.
fn finish_record(
    sequences: &mut HashMap<String, Peptide>,
    record: Option<(String, String, Option<u32>)>,
    seq: &mut String,
) {
    if let Some((tid, identifier, version)) = record {
        if !seq.is_empty() {
            sequences.insert(
                tid,
                Peptide {
                    sequence: seq.clone(),
                    identifier,
                    version,
                },
            );
        }
    }
    seq.clear();
}

/// The transcript a FASTA header names: the text after `transcript:` in the
/// first token that opens with it, whole, every dot kept (`ENST00000456.3`,
/// `AT1G01010.1`); `None` for a header with no such token.
fn extract_transcript_id(header: &str) -> Option<String> {
    header
        .split_whitespace()
        .find_map(|part| part.strip_prefix("transcript:"))
        .map(str::to_string)
}

/// The protein identifier a FASTA header opens with, its first token as
/// written; empty for a header that opens with none.
fn protein_identifier(header: &str) -> &str {
    header
        .trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or_default()
}

/// The version of the protein identifier a FASTA header opens with: the
/// digits after the last `.` of the first token (`>ENSP00000123456.2` is 2);
/// an error for digits above the largest the cache holds
/// ([`version_number`]). A token without a `.` (`>FBpp0091078`), or whose
/// last `.` is followed by anything but digits, states no version: the
/// identifier conventions of other genome projects end a protein identifier
/// in `.cds1`, `.CDS1`, `.p1`, `.1:pep`, `.1-p1`, `.1-mRNA-1` or `.t1`, where the
/// dotted part names the coding sequence or the translation, never a
/// release of the protein.
fn extract_protein_version(header: &str) -> Result<Option<u32>> {
    protein_identifier(header)
        .rsplit_once('.')
        .map(|(_, suffix)| suffix)
        .filter(|suffix| is_whole_number(suffix))
        .map(version_number)
        .transpose()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_transcript_id() {
        let header = ">ENSP00000123456.2 pep chromosome:GRCh38:1:100:500:1 gene:ENSG00000123 transcript:ENST00000456.3 gene_biotype:protein_coding";
        assert_eq!(
            extract_transcript_id(header),
            Some("ENST00000456.3".to_string())
        );
    }

    #[test]
    fn test_extract_transcript_id_no_version() {
        let header = ">ENSP00000123456 pep transcript:ENST00000789";
        assert_eq!(
            extract_transcript_id(header),
            Some("ENST00000789".to_string())
        );
    }

    #[test]
    fn test_extract_transcript_id_missing() {
        let header = ">ENSP00000123456 pep gene:ENSG00000123";
        assert_eq!(extract_transcript_id(header), None);
    }

    /// The `transcript:` token is read whole, every dot kept: an identifier
    /// of several dotted components (`AT1G01020.2`, `Y74C9A.3.1`,
    /// `HORVU.MOREX.r3.7HG0737000.1`) is the token's whole text, as is an
    /// identifier with its version (`ENST00000456.3`).
    #[test]
    fn a_transcript_token_is_read_whole() {
        for token in [
            "AT1G01020.2",
            "Y74C9A.3.1",
            "HORVU.MOREX.r3.7HG0737000.1",
            "SPAC212.11.1",
            "ENST00000456.3",
            "FBtr0091836",
        ] {
            let header =
                format!(">ENSP1 pep gene:G1 transcript:{token} gene_biotype:protein_coding");
            assert_eq!(extract_transcript_id(&header).as_deref(), Some(token));
        }
    }

    /// Records whose transcript identifiers share their first component are
    /// each kept under their own identifier with their own peptide; none
    /// takes another's place.
    #[test]
    fn records_sharing_a_first_component_are_each_kept() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        std::fs::write(
            &path,
            ">AT1G01020.1 pep transcript:AT1G01020.1\nMKV\n\
             >AT1G01020.2 pep transcript:AT1G01020.2\nMKVL\n\
             >AT1G01020.10 pep transcript:AT1G01020.10\nMK\n",
        )
        .unwrap();
        let records = parse_protein_fasta(&path).unwrap();
        let mut held: Vec<(&str, &str)> = records
            .iter()
            .map(|(id, record)| (id.as_str(), record.sequence.as_str()))
            .collect();
        held.sort_unstable();
        assert_eq!(
            held,
            [
                ("AT1G01020.1", "MKV"),
                ("AT1G01020.10", "MK"),
                ("AT1G01020.2", "MKVL")
            ]
        );
    }

    /// A model's transcript finds the record whose token is its identifier
    /// whole: a dotted identifier without a version by itself, never by its
    /// first component; a versioned one by `identifier.version`, where the
    /// record of another version is not its own, or by the identifier alone
    /// when the header writes no version; the versioned token first when the
    /// file holds both.
    #[test]
    fn a_record_is_found_by_the_models_whole_identifier() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        std::fs::write(
            &path,
            ">AT1G01020.2 pep transcript:AT1G01020.2\nMKVL\n\
             >ENSP00000123456.2 pep transcript:ENST00000456.3\nMKV\n\
             >ENSP00000123789 pep transcript:ENST00000789\nMK\n\
             >ENSP00000123999.1 pep transcript:ENST00000999.2\nMKVLW\n\
             >ENSP00000123999 pep transcript:ENST00000999\nMKVLWA\n",
        )
        .unwrap();
        let records = parse_protein_fasta(&path).unwrap();
        let alone = ModelNames::default();
        let found = |id: &str, version| {
            record_of(&records, &alone, id, version)
                .map(|(key, record)| (key, record.sequence.as_str()))
        };
        assert_eq!(found("AT1G01020.2", None), Some(("AT1G01020.2", "MKVL")));
        assert_eq!(found("AT1G01020", None), None);
        assert_eq!(found("AT1G01020", Some(2)), Some(("AT1G01020.2", "MKVL")));

        assert_eq!(
            found("ENST00000456", Some(3)),
            Some(("ENST00000456.3", "MKV"))
        );
        assert_eq!(found("ENST00000456", Some(4)), None);
        assert_eq!(found("ENST00000456", None), None);

        assert_eq!(found("ENST00000789", Some(2)), Some(("ENST00000789", "MK")));
        assert_eq!(found("ENST00000789", None), Some(("ENST00000789", "MK")));

        assert_eq!(
            found("ENST00000999", Some(2)),
            Some(("ENST00000999.2", "MKVLW"))
        );
        assert_eq!(
            found("ENST00000999", Some(3)),
            Some(("ENST00000999", "MKVLWA"))
        );
    }

    /// A token goes to the transcript written as it before one that reaches
    /// it by the dotted form. Beside a sibling written `T1.1`, the transcript
    /// `T1` at version 1 takes the record of `T1`, and none where the file
    /// holds the sibling's alone; a sibling with a version of its own takes
    /// `T1.1` only where the file holds no `T1.1.1`, and leaves it to `T1`
    /// where it does. Without the sibling the dotted token is tried first.
    #[test]
    fn a_token_written_as_a_sibling_is_the_siblings_record() {
        let file = |tokens: &[&str]| -> HashMap<String, Peptide> {
            let record = Peptide {
                sequence: "MK".to_string(),
                identifier: "ENSP1".to_string(),
                version: None,
            };
            tokens
                .iter()
                .map(|token| (token.to_string(), record.clone()))
                .collect()
        };
        let keys = |model: &[(&str, Option<u32>)], tokens: &[&str]| -> Vec<Option<String>> {
            let names = ModelNames::of(model.iter().copied());
            let records = file(tokens);
            model
                .iter()
                .map(|(id, version)| {
                    record_of(&records, &names, id, *version).map(|(key, _)| key.to_string())
                })
                .collect()
        };
        let key = |token: &str| Some(token.to_string());

        let plain_sibling = [("T1", Some(1)), ("T1.1", None)];
        assert_eq!(
            keys(&plain_sibling, &["T1", "T1.1"]),
            [key("T1"), key("T1.1")]
        );
        assert_eq!(keys(&plain_sibling, &["T1.1"]), [None, key("T1.1")]);

        let versioned_sibling = [("T1", Some(1)), ("T1.1", Some(1))];
        assert_eq!(
            keys(&versioned_sibling, &["T1", "T1.1"]),
            [key("T1"), key("T1.1")]
        );
        assert_eq!(
            keys(&versioned_sibling, &["T1.1", "T1.1.1"]),
            [key("T1.1"), key("T1.1.1")]
        );
        assert_eq!(keys(&versioned_sibling, &["T1.1"]), [None, key("T1.1")]);

        let three = [("T1", Some(1)), ("T1.1", Some(1)), ("T1.1.1", None)];
        assert_eq!(
            keys(&three, &["T1", "T1.1", "T1.1.1"]),
            [key("T1"), key("T1.1"), key("T1.1.1")]
        );

        assert_eq!(keys(&[("T1", Some(1))], &["T1", "T1.1"]), [key("T1.1")]);
    }

    /// A header whose protein identifier is the model's own, whole, states no
    /// version for it: the `1` that ends `AT1G01010.1` is the identifier's.
    /// Beside the identifier without its last component the header's digits
    /// are the version, as they are for a model that names no protein.
    #[test]
    fn a_header_identifier_that_is_the_models_states_no_version() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        std::fs::write(
            &path,
            ">AT1G01010.1 pep transcript:AT1G01010.1\nMKV\n\
             >ENSP00000123456.2 pep transcript:ENST00000456.3\nMK\n",
        )
        .unwrap();
        let records = parse_protein_fasta(&path).unwrap();
        let dotted = &records["AT1G01010.1"];
        assert_eq!(
            (dotted.identifier.as_str(), dotted.version),
            ("AT1G01010.1", Some(1))
        );
        assert_eq!(dotted.version_beside(Some("AT1G01010.1")), None);
        assert_eq!(dotted.version_beside(Some("AT1G01010")), Some(1));
        assert_eq!(dotted.version_beside(None), Some(1));

        let versioned = &records["ENST00000456.3"];
        assert_eq!(versioned.version_beside(Some("ENSP00000123456")), Some(2));
        assert_eq!(versioned.version_beside(Some("ENSP00000123456.2")), None);
    }

    /// A record keeps the protein identifier its header opens with and that
    /// identifier's version, the translation's version in the release: 2 from
    /// `>ENSP00000123456.2`, 1 from `>ENSP00000123457.01`, none from an
    /// identifier written without one (`>FBpp0091078`, the shape of a release
    /// whose identifiers carry no version); the transcript's own version in
    /// the `transcript:` token is not the translation's.
    #[test]
    fn a_record_keeps_the_header_protein_version() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        std::fs::write(
            &path,
            ">ENSP00000123456.2 pep transcript:ENST00000456.3\nMKV\n\
             >FBpp0091078 pep transcript:FBtr0091836\nMKVL\n\
             >ENSP00000123457.01 pep transcript:ENST00000457.1\nMK\n",
        )
        .unwrap();
        let peptide = |sequence: &str, identifier: &str, version| Peptide {
            sequence: sequence.to_string(),
            identifier: identifier.to_string(),
            version,
        };
        assert_eq!(
            parse_protein_fasta(&path).unwrap(),
            HashMap::from([
                (
                    "ENST00000456.3".to_string(),
                    peptide("MKV", "ENSP00000123456.2", Some(2))
                ),
                (
                    "FBtr0091836".to_string(),
                    peptide("MKVL", "FBpp0091078", None)
                ),
                (
                    "ENST00000457.1".to_string(),
                    peptide("MK", "ENSP00000123457.01", Some(1))
                ),
            ])
        );
    }

    /// A header whose first token ends in a `.` followed by anything but
    /// digits states no version: the dotted part is the identifier's own,
    /// as the protein identifiers of other genome projects write it
    /// (`.cds1`, `.CDS1`, `.p1`, `.1:pep`, `.1-p1`, `.t1`, a trailing `.`),
    /// so the record keeps its peptide with no version and the build goes on.
    #[test]
    fn a_header_suffix_that_is_not_digits_is_part_of_the_identifier() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        let text: String = [
            ">TraesCS7D02G255100.1.cds1 pep transcript:TraesCS7D02G255100.1\nMKV\n",
            ">HORVU.MOREX.r3.7HG0737000.1.CDS1 pep transcript:HORVU.MOREX.r3.7HG0737000.1\nMKVL\n",
            ">cds.TCONS_00101750.p1 pep transcript:TCONS_00101750\nMK\n",
            ">FCUL_00001.1-mRNA-1 pep transcript:FCUL_00001.1-mRNA-1\nMKVLW\n",
            ">SPAC212.11.1:pep pep transcript:SPAC212.11.1\nMKVLWA\n",
            ">PF3D7_0528000.1-p1 pep transcript:PF3D7_0528000.1\nMKVLWAA\n",
            ">g00297.t1 pep transcript:g00297.t1\nMKVLWAAL\n",
            ">ENSP00000123457.4294967296a pep transcript:ENST00000457.1\nMKVLWAALL\n",
            ">ENSP00000123458. pep transcript:ENST00000458.2\nMKVLWAALLL\n",
        ]
        .concat();
        std::fs::write(&path, text).unwrap();
        let versions: Vec<Option<u32>> = parse_protein_fasta(&path)
            .unwrap()
            .into_values()
            .map(|peptide| peptide.version)
            .collect();
        assert_eq!(versions.len(), 9);
        assert!(versions.iter().all(Option::is_none), "{versions:?}");
    }

    /// A header whose protein version does not fit the cache's 32-bit field
    /// fails the parse, naming the file, the line and the limit, since a
    /// record read without the version would lose what the header states;
    /// 4294967295 is held.
    #[test]
    fn a_header_version_above_the_cache_limit_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        std::fs::write(
            &path,
            ">ENSP00000123456.4294967295 pep transcript:ENST00000456.3\nMKV\n\
             >ENSP00000123457.4294967296 pep transcript:ENST00000457.1\nMK\n",
        )
        .unwrap();
        let err = parse_protein_fasta(&path).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(
                "pep.fa line 3: version 4294967296 is above 4294967295, \
                 the largest version the cache holds"
            ),
            "{err:#}"
        );
        std::fs::write(
            &path,
            ">ENSP00000123456.4294967295 pep transcript:ENST00000456.3\nMKV\n",
        )
        .unwrap();
        assert_eq!(
            parse_protein_fasta(&path).unwrap()["ENST00000456.3"].version,
            Some(4294967295)
        );
    }

    /// A line that is not UTF-8 fails the parse, naming the file and the line,
    /// as every refusal of a line does.
    #[test]
    fn a_line_that_is_not_utf8_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        std::fs::write(
            &path,
            b">ENSP00000123456.2 pep transcript:ENST00000456.3 \xe9\nMKV\n",
        )
        .unwrap();
        let err = parse_protein_fasta(&path).unwrap_err();
        assert!(
            format!("{err:#}").ends_with("pep.fa line 1: FASTA line is not UTF-8"),
            "{err:#}"
        );
    }

    /// A bgzipped protein FASTA, one gzip member per 64 KiB block, is read to
    /// its last record: every peptide of a file spanning several blocks is
    /// mapped, as from the same file written plain.
    #[test]
    fn a_bgzipped_fasta_is_read_past_its_first_block() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let mut text = String::new();
        for i in 0..2_000 {
            text.push_str(&format!(
                ">ENSP{i:011} pep transcript:ENST{i:011}.1\n{}\n",
                "MKV".repeat(20)
            ));
        }
        let plain = dir.path().join("pep.fa");
        std::fs::write(&plain, &text).unwrap();
        let bgzipped = dir.path().join("pep.fa.gz");
        let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(&bgzipped).unwrap());
        writer.write_all(text.as_bytes()).unwrap();
        writer.finish().unwrap();
        assert!(crate::fasta_index::is_bgzf(&bgzipped).unwrap());
        // BGZF blocks hold at most 64 KiB of input, so the file has three or more.
        assert!(text.len() > 2 * 65_536);

        let from_bgzip = parse_protein_fasta(&bgzipped).unwrap();
        assert_eq!(from_bgzip.len(), 2_000);
        assert_eq!(from_bgzip["ENST00000001999.1"].sequence, "MKV".repeat(20));
        assert_eq!(from_bgzip, parse_protein_fasta(&plain).unwrap());
    }

    /// A sequence line one byte over the cap fails the parse, naming the file
    /// and the line, 2; a file within the cap parses.
    #[test]
    fn line_past_the_cap_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        let header = ">ENSP00000123456.2 pep transcript:ENST00000456.3\n";
        std::fs::write(&path, format!("{header}MKV\n")).unwrap();
        assert_eq!(
            parse_protein_fasta(&path).unwrap()["ENST00000456.3"].sequence,
            "MKV"
        );
        let mut text = header.as_bytes().to_vec();
        text.extend(std::iter::repeat_n(b'M', MAX_LINE_BYTES + 1));
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        let err = parse_protein_fasta(&path).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "pep.fa line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
    }
}
