// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader of an engine's default-format output: its tuples, its Locations and the variant type of each printed Location.
//!
//! A file is read as a stream in runs of whole lines, and each run is split and normalised
//! on whichever thread is free. A column text is held once, in an [`Interner`] all the
//! outputs of a run share, and a row is five numbers from the moment it is read; when every
//! output is read the texts are ranked into the run's [`Strings`] and the numbers become
//! the ranks.

use std::collections::{BTreeMap, HashMap};
use std::hash::{BuildHasherDefault, Hasher};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, PoisonError};

use memchr::memchr_iter;

use crate::failure::Failure;
use crate::input::KeyPrep;
use crate::normalize::{
    fastvep_allele_integer_digits, is_fastvep, location_integer_digits, normalize_allele_of,
    push_normalized_location,
};
use crate::strings::{Id, Strings, Tuple, TupleSet};
use crate::text::{
    first_defect, line_count, lines, open_file, os_failure, read_blocks, reading_failure,
    LineDefect, RunDefect, BLOCK_BYTES,
};
use crate::threads::{read_sources, BlockSource};
use crate::vcf::VcfVariant;

/// The number a text has while outputs are still being read: the order texts were first
/// met in, which differs from run to run and is never written anywhere.
pub type Provisional = u32;

const SHARDS: usize = 256;

/// FNV-1a over the bytes hashed.
struct Fnv(u64);

impl Default for Fnv {
    fn default() -> Fnv {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
}

impl Hasher for Fnv {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 = (self.0 ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
}

type Shard = HashMap<Box<str>, Provisional, BuildHasherDefault<Fnv>>;

/// Every distinct column text met so far, each with its [`Provisional`] number, shared by
/// the threads that read the outputs of a run.
pub struct Interner {
    shards: Vec<Mutex<Shard>>,
    next: AtomicU32,
}

impl Default for Interner {
    fn default() -> Interner {
        Interner::new()
    }
}

impl Interner {
    pub fn new() -> Interner {
        Interner {
            shards: (0..SHARDS).map(|_| Mutex::new(Shard::default())).collect(),
            next: AtomicU32::new(0),
        }
    }

    /// The number of a text, the same for every thread that asks.
    pub fn intern(&self, text: &str) -> Provisional {
        let mut hasher = Fnv::default();
        hasher.write(text.as_bytes());
        // Bits the multiplications have mixed and that a shard's own table, which hashes
        // the text again, does not index by.
        let shard = (hasher.finish() >> 40) as usize % SHARDS;
        let mut texts = self.shards[shard]
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        if let Some(number) = texts.get(text) {
            return *number;
        }
        let number = self.next.fetch_add(1, Ordering::Relaxed);
        assert!(
            number < Provisional::MAX,
            "more distinct column texts than an Id counts"
        );
        texts.insert(text.into(), number);
        number
    }

    /// The table of every text met, and the rank of each provisional number in it.
    pub fn finish(self) -> (Strings, Ranks) {
        let mut texts: Vec<(Box<str>, Provisional)> = Vec::new();
        for shard in self.shards {
            texts.extend(shard.into_inner().unwrap_or_else(PoisonError::into_inner));
        }
        texts.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        let mut ranks = vec![0; texts.len()];
        let sorted: Vec<Box<str>> = texts
            .into_iter()
            .enumerate()
            .map(|(rank, (text, number))| {
                ranks[number as usize] = rank as Id;
                text
            })
            .collect();
        (Strings::from_sorted_unique(sorted), Ranks(ranks))
    }
}

/// The rank in the run's [`Strings`] of each [`Provisional`] number.
pub struct Ranks(Vec<Id>);

impl Ranks {
    pub fn id(&self, number: Provisional) -> Id {
        self.0[number as usize]
    }

    /// The set of the rows read, each column as its rank.
    pub fn tuple_set(&self, rows: &[[Provisional; 5]]) -> TupleSet {
        rows.iter()
            .map(|row| Tuple(row.map(|number| self.id(number))))
            .collect()
    }
}

/// The Locations a set of tuples carries, sorted, each once.
pub fn locations_of(tuples: &TupleSet) -> Vec<Id> {
    let mut locations: Vec<Id> = tuples.iter().map(Tuple::location).collect();
    locations.dedup();
    locations
}

/// The variant type of every Location that has one, the Location as its rank in a run's
/// [`Strings`].
#[derive(Clone, Debug, Default)]
pub struct LocationTypes {
    names: Vec<Box<str>>,
    by_location: HashMap<Id, u32>,
}

impl LocationTypes {
    fn name_index(&mut self, variant_type: &str) -> u32 {
        match self.names.iter().position(|name| **name == *variant_type) {
            Some(index) => index as u32,
            None => {
                self.names.push(variant_type.into());
                self.names.len() as u32 - 1
            }
        }
    }

    pub fn len(&self) -> usize {
        self.by_location.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_location.is_empty()
    }

    pub fn get(&self, location: Id) -> Option<&str> {
        self.by_location
            .get(&location)
            .map(|index| &*self.names[*index as usize])
    }

    /// Gives the Location this type, whatever it had.
    pub fn set(&mut self, location: Id, variant_type: &str) {
        let index = self.name_index(variant_type);
        self.by_location.insert(location, index);
    }

    /// Gives the Location this type unless it has one.
    pub fn set_default(&mut self, location: Id, variant_type: &str) {
        if !self.by_location.contains_key(&location) {
            self.set(location, variant_type);
        }
    }

    /// Gives every Location of `other` the type it has there.
    pub fn update(&mut self, other: &LocationTypes) {
        let indices: Vec<u32> = other
            .names
            .iter()
            .map(|name| self.name_index(name))
            .collect();
        self.by_location.reserve(other.by_location.len());
        for (location, index) in &other.by_location {
            self.by_location.insert(*location, indices[*index as usize]);
        }
    }

    /// Every Location with its type, in no order.
    pub fn iter(&self) -> impl Iterator<Item = (Id, &str)> {
        self.by_location
            .iter()
            .map(|(location, index)| (*location, &*self.names[*index as usize]))
    }

    /// The types of the Locations a table holds, from a map by Location text.
    pub fn from_texts(strings: &Strings, types: &BTreeMap<String, String>) -> LocationTypes {
        let mut out = LocationTypes::default();
        for (location, variant_type) in types {
            if let Some(id) = strings.get(location) {
                out.set(id, variant_type);
            }
        }
        out
    }
}

/// The variant type of each record ID of one input file, an ID being the first identifier
/// of the ID column ([`VcfVariant::vid`]): the first record of an ID gives it its type, and
/// a record without an ID (empty or `.`) gives none. The types are indices into `names`.
pub struct TypesById<'a> {
    pub names: Vec<&'a str>,
    by_id: HashMap<&'a str, u32>,
}

impl<'a> TypesById<'a> {
    pub fn of(input_variants: &'a [VcfVariant]) -> TypesById<'a> {
        let mut names: Vec<&'a str> = Vec::new();
        let mut by_id: HashMap<&'a str, u32> = HashMap::new();
        for variant in input_variants {
            if variant.vid.is_empty() || variant.vid == "." {
                continue;
            }
            by_id.entry(variant.vid.as_str()).or_insert_with(|| {
                let name = variant.variant_type.as_str();
                match names.iter().position(|held| *held == name) {
                    Some(index) => index as u32,
                    None => {
                        names.push(name);
                        names.len() as u32 - 1
                    }
                }
            });
        }
        TypesById { names, by_id }
    }
}

/// What one run of whole lines of an output holds.
#[derive(Debug, Default)]
pub struct OutputBlock {
    /// How many lines the run has.
    pub lines: u64,
    /// The rows before the first defective line, in file order.
    pub rows: Vec<[Provisional; 5]>,
    /// Location and type index of rows whose first column is the ID of a typed record, in
    /// file order; within a run of rows at one Location, the first such row alone.
    pub printed: Vec<(Provisional, u32)>,
    pub defect: Option<RunDefect>,
}

/// One engine's output as a source of runs of lines.
pub struct OutputSource<'a> {
    pub path: PathBuf,
    /// Whether fastVEP wrote the output: a raw symbolic token in its Allele column, which
    /// fastVEP through v0.5.0 writes, is then mapped to its class name.
    pub fastvep: bool,
    pub interner: &'a Interner,
    /// Whether the rows are read as tuples: a row of fewer than seven columns is then a
    /// defect. Otherwise rows of three columns or more are read for their Location alone.
    pub tuples: bool,
    /// The record types the first column of a row is looked up in, when the printed
    /// Locations are typed.
    pub types: Option<&'a TypesById<'a>>,
    /// The bytes a run holds, about.
    pub block_bytes: usize,
}

/// The last text of a column and the number it was given.
#[derive(Default)]
struct Last {
    text: String,
    number: Option<Provisional>,
}

impl Last {
    fn holds(&self, text: &str) -> Option<Provisional> {
        self.number.filter(|_| self.text == text)
    }

    fn hold(&mut self, text: &str, number: Provisional) -> Provisional {
        self.text.clear();
        self.text.push_str(text);
        self.number = Some(number);
        number
    }
}

/// What a thread keeps between the lines of one run.
#[derive(Default)]
struct RunState {
    location: Last,
    allele: Last,
    feature: Last,
    feature_type: Last,
    consequence: Last,
    record_id: String,
    record_type: Option<u32>,
    typed_location: Option<Provisional>,
    scratch: Vec<u8>,
    prep: KeyPrep,
}

impl<'a> OutputSource<'a> {
    pub fn new(
        path: PathBuf,
        engine: &str,
        interner: &'a Interner,
        types: Option<&'a TypesById<'a>>,
    ) -> OutputSource<'a> {
        OutputSource {
            path,
            fastvep: is_fastvep(engine),
            interner,
            tuples: true,
            types,
            block_bytes: BLOCK_BYTES,
        }
    }

    fn location(&self, state: &mut RunState, written: &str) -> Result<Provisional, LineDefect> {
        if let Some(number) = state.location.holds(written) {
            return Ok(number);
        }
        if let Some(digits) = location_integer_digits(written) {
            return Err(LineDefect::IntegerDigits(digits));
        }
        state.scratch.clear();
        push_normalized_location(written, &mut state.scratch);
        let compared = std::str::from_utf8(&state.scratch).expect("cut at ASCII bytes of a str");
        Ok(state.location.hold(written, self.interner.intern(compared)))
    }

    fn allele(&self, state: &mut RunState, written: &str) -> Result<Provisional, LineDefect> {
        if let Some(number) = state.allele.holds(written) {
            return Ok(number);
        }
        if self.fastvep {
            if let Some(digits) = fastvep_allele_integer_digits(written) {
                return Err(LineDefect::IntegerDigits(digits));
            }
        }
        let compared = normalize_allele_of(written, self.fastvep);
        Ok(state.allele.hold(written, self.interner.intern(&compared)))
    }

    fn verbatim(&self, last: &mut Last, written: &str) -> Provisional {
        match last.holds(written) {
            Some(number) => number,
            None => last.hold(written, self.interner.intern(written)),
        }
    }

    fn consequence(&self, state: &mut RunState, written: &str) -> Provisional {
        if let Some(number) = state.consequence.holds(written) {
            return number;
        }
        state.scratch.clear();
        state.prep.normalized_set(written, &mut state.scratch);
        let compared = std::str::from_utf8(&state.scratch).expect("the terms of a str");
        state
            .consequence
            .hold(written, self.interner.intern(compared))
    }

    /// Reads one line into `block`; a comment line, an empty line and, where rows are read
    /// for their Location alone, a row of fewer than three columns hold nothing.
    fn line(
        &self,
        state: &mut RunState,
        line: &[u8],
        block: &mut OutputBlock,
    ) -> Result<(), LineDefect> {
        let line = std::str::from_utf8(line).map_err(|_| LineDefect::NotUtf8)?;
        if line.is_empty() || line.starts_with('#') {
            return Ok(());
        }
        // The tabs that end the columns read: seven of a tuple row, whose seventh column
        // may be the last, and the first two of a row read for its Location.
        let (wanted, needed) = if self.tuples { (7, 6) } else { (2, 2) };
        let mut tabs = [line.len(); 7];
        let mut found = 0;
        for tab in memchr_iter(b'\t', line.as_bytes()).take(wanted) {
            tabs[found] = tab;
            found += 1;
        }
        if found < needed {
            return if self.tuples {
                Err(LineDefect::ShortRow(found + 1))
            } else {
                Ok(())
            };
        }
        let column = |index: usize| -> &str {
            let start = if index == 0 { 0 } else { tabs[index - 1] + 1 };
            &line[start..tabs[index]]
        };
        let record_type = self.types.and_then(|types| {
            let record_id = column(0);
            if state.record_id != record_id {
                state.record_id.clear();
                state.record_id.push_str(record_id);
                state.record_type = types.by_id.get(record_id).copied();
            }
            state.record_type
        });
        // The Location of a row read for it alone is read once a record has the row's ID:
        // a Location that cannot be read fails the row of a typed record and no other.
        if !self.tuples && record_type.is_none() {
            return Ok(());
        }
        let location = self.location(state, column(1))?;
        if self.tuples {
            let allele = self.allele(state, column(2))?;
            let feature = self.verbatim(&mut state.feature, column(4));
            let feature_type = self.verbatim(&mut state.feature_type, column(5));
            let consequence = self.consequence(state, column(6));
            let row = [location, allele, feature, feature_type, consequence];
            if block.rows.last() != Some(&row) {
                block.rows.push(row);
            }
        }
        if let Some(record_type) = record_type {
            if state.typed_location != Some(location) {
                state.typed_location = Some(location);
                block.printed.push((location, record_type));
            }
        }
        Ok(())
    }

    /// What the file holds, from its parsed runs: `None` for a file that does not exist,
    /// or the failure its reading ends in when a line of it is defective.
    pub fn finish(
        &self,
        read: Result<bool, Failure>,
        blocks: Vec<OutputBlock>,
    ) -> Result<Option<OutputRows>, Failure> {
        // A line is read, and can fail, before the file is known to end badly.
        let defect = first_defect(blocks.iter().map(|block| (block.lines, block.defect)));
        if let Some((line, defect)) = defect {
            return Err(reading_failure(
                &self.path,
                open_file(&self.path),
                line,
                defect,
            ));
        }
        if !read? {
            return Ok(None);
        }
        let mut rows = Vec::with_capacity(blocks.iter().map(|block| block.rows.len()).sum());
        let mut printed = Vec::with_capacity(blocks.iter().map(|block| block.printed.len()).sum());
        for mut block in blocks {
            rows.append(&mut block.rows);
            printed.append(&mut block.printed);
        }
        Ok(Some(OutputRows { rows, printed }))
    }
}

/// The rows of one output in file order, and the typed Locations among them; see
/// [`OutputBlock`].
#[derive(Debug, Default)]
pub struct OutputRows {
    pub rows: Vec<[Provisional; 5]>,
    pub printed: Vec<(Provisional, u32)>,
}

impl BlockSource for OutputSource<'_> {
    /// Whether the file exists.
    type Read = Result<bool, Failure>;
    type Parsed = OutputBlock;

    fn read(&self, emit: &mut dyn FnMut(Vec<u8>)) -> Result<bool, Failure> {
        if std::fs::metadata(&self.path).is_err() {
            return Ok(false);
        }
        let file = open_file(&self.path).map_err(|error| os_failure(&error, Some(&self.path)))?;
        read_blocks(file, self.block_bytes, emit).map_err(|error| os_failure(&error, None))?;
        Ok(true)
    }

    fn parse(&self, _sequence: u64, block: Vec<u8>) -> OutputBlock {
        let mut parsed = OutputBlock {
            lines: line_count(&block),
            ..OutputBlock::default()
        };
        let mut state = RunState::default();
        for (index, line) in lines(&block).enumerate() {
            if let Err(defect) = self.line(&mut state, line, &mut parsed) {
                parsed.defect = Some((index as u64 + 1, defect));
                break;
            }
        }
        parsed
    }
}

/// One output of several read into one table: its tuples, its Locations and whether the
/// file exists.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct ParsedOutput {
    pub tuples: TupleSet,
    /// The Locations of the tuples, sorted.
    pub locations: Vec<Id>,
    pub exists: bool,
}

/// Reads each `(path, engine)` output into one table of texts, on at most `threads`
/// threads; see [`parse_vep_output`].
pub fn parse_vep_outputs(
    outputs: &[(&Path, &str)],
    threads: usize,
) -> Result<(Strings, Vec<ParsedOutput>), Failure> {
    let interner = Interner::new();
    let sources: Vec<OutputSource> = outputs
        .iter()
        .map(|(path, engine)| OutputSource::new(path.to_path_buf(), engine, &interner, None))
        .collect();
    let mut rows = Vec::with_capacity(sources.len());
    for (source, (read, blocks)) in sources.iter().zip(read_sources(threads, &sources)) {
        rows.push(source.finish(read, blocks)?);
    }
    drop(sources);
    let (strings, ranks) = interner.finish();
    let parsed = rows
        .into_iter()
        .map(|rows| {
            let exists = rows.is_some();
            let tuples = ranks.tuple_set(&rows.unwrap_or_default().rows);
            ParsedOutput {
                locations: locations_of(&tuples),
                tuples,
                exists,
            }
        })
        .collect();
    Ok((strings, parsed))
}

/// The tuples and the Locations of one engine's tab-delimited output, with the table their
/// columns are ranks in. `engine` names the engine that wrote the file; fastVEP alone
/// changes anything, a raw symbolic ALT token, which fastVEP through v0.5.0 writes, being
/// mapped to the class name the other engines write.
///
/// A file that does not exist is two empty sets: the engine wrote no output for the input.
/// A row of fewer than seven columns is a failure, since leaving it out would turn a
/// truncated file into a smaller tuple set that scores as better concordance.
pub fn parse_vep_output(
    path: &Path,
    engine: &str,
) -> Result<(Strings, TupleSet, Vec<Id>), Failure> {
    let (strings, mut parsed) = parse_vep_outputs(&[(path, engine)], 1)?;
    let only = parsed.pop().expect("one output");
    Ok((strings, only.tuples, only.locations))
}

/// The variant type of every Location: of each input record's `chrom:pos`, a later record
/// at a position replacing an earlier one, and of every other Location an output prints
/// for a record with an ID, the first row to print a Location deciding its type, the
/// outputs in the order given and their rows in file order.
///
/// A record's `chrom:pos` equals only a point Location; the engines print a symbolic or
/// breakend record at another one, so the join is through the first column of a row, the
/// record's ID: the first of the identifiers its ID column holds. A record without an ID
/// keeps its `chrom:pos` alone, and the alleles of a multi-allelic record share its ID and
/// its type.
pub fn location_types_from_outputs(
    input_variants: &[VcfVariant],
    output_paths: &[&Path],
) -> Result<BTreeMap<String, String>, Failure> {
    let mut types: BTreeMap<String, String> = BTreeMap::new();
    for variant in input_variants {
        types.insert(variant.location.clone(), variant.variant_type.clone());
    }
    let by_id = TypesById::of(input_variants);
    let interner = Interner::new();
    let sources: Vec<OutputSource> = output_paths
        .iter()
        .map(|path| OutputSource {
            tuples: false,
            ..OutputSource::new(path.to_path_buf(), "", &interner, Some(&by_id))
        })
        .collect();
    let mut printed = Vec::new();
    for (source, (read, blocks)) in sources.iter().zip(read_sources(1, &sources)) {
        printed.extend(source.finish(read, blocks)?.unwrap_or_default().printed);
    }
    drop(sources);
    let (strings, ranks) = interner.finish();
    for (location, record_type) in printed {
        types
            .entry(strings.text(ranks.id(location)).to_string())
            .or_insert_with(|| by_id.names[record_type as usize].to_string());
    }
    Ok(types)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::text::INT_MAX_STR_DIGITS;
    use crate::TempDir;

    /// Writes a 17-column tab output: one `##` line, the column header, then each row's
    /// seven columns and ten dashes.
    pub(crate) fn write_vep(path: &Path, rows: &[[&str; 7]]) {
        let mut text = String::from(
            "## fastVEP output\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\n",
        );
        for row in rows {
            text.push_str(&row.join("\t"));
            text.push_str(&"\t-".repeat(10));
            text.push('\n');
        }
        std::fs::write(path, text).unwrap();
    }

    const DEL_ROW: [&str; 7] = [
        "v1",
        "21:100",
        "<DEL>",
        "G1",
        "ENST1",
        "Transcript",
        "feature_truncation",
    ];
    const DUP_ROW: [&str; 7] = [
        "v2",
        "21:200",
        "<DUP>",
        "G2",
        "ENST2",
        "Transcript",
        "feature_elongation",
    ];

    fn append(path: &Path, text: &str) {
        let mut held = std::fs::read(path).unwrap();
        held.extend_from_slice(text.as_bytes());
        std::fs::write(path, held).unwrap();
    }

    fn texts<'a>(strings: &'a Strings, tuples: &TupleSet) -> Vec<[&'a str; 5]> {
        tuples.iter().map(|tuple| strings.resolve(tuple)).collect()
    }

    #[test]
    fn test_wellformed_file_parses() {
        let dir = TempDir::new("output_ok_").unwrap();
        let path = dir.path.join("ok.txt");
        write_vep(&path, &[DEL_ROW, DUP_ROW]);
        let (strings, tuples, locations) = parse_vep_output(&path, "vep-rs").unwrap();
        assert_eq!(tuples.len(), 2);
        let locations: Vec<&str> = locations.iter().map(|id| strings.text(*id)).collect();
        assert_eq!(locations, ["21:100", "21:200"]);
    }

    #[test]
    fn test_truncated_midline_raises() {
        let dir = TempDir::new("output_trunc_").unwrap();
        let path = dir.path.join("trunc.txt");
        write_vep(&path, &[DEL_ROW]);
        append(&path, "v2\t21:200\t<DUP>\tG2\n");
        let failure = parse_vep_output(&path, "vep-rs").unwrap_err();
        assert_eq!(failure.kind, FailureKind::Value);
        assert!(failure.message.contains("expected at least 7"));
        assert_eq!(
            failure.message,
            format!(
                "{}: line 4 has 4 columns, expected at least 7. A truncated or malformed engine output must fail the comparison rather than silently reduce the tuple set.",
                path.display()
            )
        );
    }

    #[test]
    fn test_missing_file_returns_empty_without_raising() {
        let dir = TempDir::new("output_missing_").unwrap();
        let (strings, tuples, locations) =
            parse_vep_output(&dir.path.join("does_not_exist.txt"), "vep-rs").unwrap();
        assert!(strings.is_empty());
        assert!(tuples.is_empty());
        assert!(locations.is_empty());
    }

    #[test]
    fn test_header_and_blank_lines_skipped() {
        let dir = TempDir::new("output_hdr_").unwrap();
        let path = dir.path.join("hdr.txt");
        write_vep(&path, &[DEL_ROW]);
        append(&path, "\n");
        let (_, tuples, _) = parse_vep_output(&path, "vep-rs").unwrap();
        assert_eq!(tuples.len(), 1);
    }

    fn deletion_row(allele: &str) -> [&str; 7] {
        [
            "v1",
            "21:1000-2000",
            allele,
            "G1",
            "ENST_A",
            "Transcript",
            "feature_truncation",
        ]
    }

    #[test]
    fn test_fastvep_del_matches_perl_deletion_under_the_fastvep_engine() {
        let dir = TempDir::new("output_fastvep_").unwrap();
        let (perl_path, fast_path) = (dir.path.join("perl.txt"), dir.path.join("fastvep.txt"));
        write_vep(&perl_path, &[deletion_row("deletion")]);
        write_vep(&fast_path, &[deletion_row("<DEL>")]);
        let (strings, parsed) =
            parse_vep_outputs(&[(&perl_path, "vep-rs"), (&fast_path, "fastvep")], 2).unwrap();
        let (perl, fast) = (&parsed[0].tuples, &parsed[1].tuples);
        assert_eq!(perl, fast, "the two rows must normalise to one tuple");
        let metrics = crate::metrics::compute_metrics(&strings, &[], perl, &[], fast, &[]);
        assert_eq!(metrics.f1, 1.0);
    }

    #[test]
    fn test_fastvep_del_does_not_match_without_the_scored_engine_flag() {
        let dir = TempDir::new("output_no_flag_").unwrap();
        let (perl_path, fast_path) = (dir.path.join("perl.txt"), dir.path.join("fastvep.txt"));
        write_vep(&perl_path, &[deletion_row("deletion")]);
        write_vep(&fast_path, &[deletion_row("<DEL>")]);
        let (strings, parsed) = parse_vep_outputs(
            &[
                (&perl_path, "vep-rs"),
                (&fast_path, "vep-rs"),
                (&fast_path, "vep-rs"),
            ],
            1,
        )
        .unwrap();
        let (perl, default, as_vep_rs) = (&parsed[0].tuples, &parsed[1].tuples, &parsed[2].tuples);
        assert!(perl.intersection(default).is_empty());
        assert!(perl.intersection(as_vep_rs).is_empty());
        let metrics = crate::metrics::compute_metrics(&strings, &[], perl, &[], default, &[]);
        assert_eq!(metrics.f1, 0.0);
    }

    #[test]
    fn a_row_is_five_columns_none_of_them_stripped() {
        let dir = TempDir::new("output_columns_").unwrap();
        let path = dir.path.join("out.txt");
        std::fs::write(
            &path,
            "#h\n\
             id\tchr21:+0100-0200\tN[chr1:5[\tG\t ENST1 \t\tb, a ,b\n\
             id\tchrM:7\t<DEL>\tG\tENST2\tTranscript\t\textra\tcolumns\r\n\
             \r\
             id\tnoloc\t-\tG\t-\t-\tx,,\rid\tchr21:+0100-0200\tN[chr1:5[\tG\t ENST1 \t\ta,b",
        )
        .unwrap();
        let (strings, tuples, locations) = parse_vep_output(&path, "perl").unwrap();
        assert_eq!(
            texts(&strings, &tuples),
            [
                ["21:100-200", "N[1:5[", " ENST1 ", "", "a,b"],
                ["MT:7", "<DEL>", "ENST2", "Transcript", ""],
                ["noloc", "-", "-", "-", "x"],
            ]
        );
        assert_eq!(locations.len(), 3);
        let (strings, tuples, _) = parse_vep_output(&path, " FastVEP ").unwrap();
        assert!(texts(&strings, &tuples).contains(&[
            "MT:7",
            "deletion",
            "ENST2",
            "Transcript",
            ""
        ]));
    }

    #[test]
    fn the_first_defective_line_is_the_failure_whatever_the_runs_and_threads() {
        let dir = TempDir::new("output_defect_").unwrap();
        let path = dir.path.join("out.txt");
        let mut text = String::from("## header\n");
        for n in 0..300 {
            text.push_str(&format!(
                "id{n}\t21:{n}\tA\tG\tENST{n}\tTranscript\tintron_variant\n"
            ));
        }
        text.push_str("short\trow\n");
        for n in 0..50 {
            text.push_str(&format!("id\t21:{n}\tA\n"));
        }
        std::fs::write(&path, &text).unwrap();
        let expected = format!(
            "{}: line 302 has 2 columns, expected at least 7. A truncated or malformed engine output must fail the comparison rather than silently reduce the tuple set.",
            path.display()
        );
        for block_bytes in [1, 50, 4096, BLOCK_BYTES] {
            for threads in [1, 3, 8] {
                let interner = Interner::new();
                let source = OutputSource {
                    block_bytes,
                    ..OutputSource::new(path.clone(), "vep-rs", &interner, None)
                };
                let (read, blocks) = read_sources(threads, std::slice::from_ref(&source))
                    .pop()
                    .unwrap();
                let failure = source.finish(read, blocks).unwrap_err();
                assert_eq!(failure.message, expected, "{block_bytes} {threads}");
            }
        }
        let mut bytes = text.into_bytes();
        bytes[60] = 0xff;
        std::fs::write(&path, &bytes).unwrap();
        let failure = parse_vep_output(&path, "vep-rs").unwrap_err();
        assert_eq!(failure.kind, FailureKind::Decode);
        assert_eq!(
            failure.message,
            "'utf-8' codec can't decode byte 0xff in position 60: invalid start byte"
        );
        // A byte that is not UTF-8 after the short row is met first when it stands in the
        // 8,192 bytes that hold the end of the short row's line, and after the row when it
        // stands further on.
        bytes[60] = b'A';
        let row_end = bytes
            .windows(10)
            .position(|window| window == b"short\trow\n")
            .unwrap()
            + 10;
        assert_eq!(row_end, 14990);
        bytes.resize(17000, b'x');
        bytes[row_end + 20] = 0xff;
        std::fs::write(&path, &bytes).unwrap();
        let failure = parse_vep_output(&path, "vep-rs").unwrap_err();
        assert_eq!(
            failure.message,
            "'utf-8' codec can't decode byte 0xff in position 6818: invalid start byte"
        );
        bytes[row_end + 20] = b'A';
        bytes[16500] = 0xff;
        std::fs::write(&path, &bytes).unwrap();
        let failure = parse_vep_output(&path, "vep-rs").unwrap_err();
        assert_eq!(failure.kind, FailureKind::Value);
        assert_eq!(failure.message, expected);
    }

    #[test]
    fn runs_of_any_size_on_any_threads_give_the_same_sets() {
        let dir = TempDir::new("output_runs_").unwrap();
        let (perl_path, rust_path) = (dir.path.join("perl.txt"), dir.path.join("rust.txt"));
        let rows = |shift: usize| -> String {
            let mut text = String::from("## header\r\n#Uploaded_variation\tLocation\n");
            for n in 0..400 {
                let m = (n * 7 + shift) % 211;
                text.push_str(&format!(
                    "id{m}\tchr{}:{m}-{}\tA\tG\tENST{}\tTranscript\tz_term,a_term{}\tmore\n",
                    m % 3 + 1,
                    m + shift,
                    n % 5,
                    if n % 2 == 0 { ",a_term" } else { "" },
                ));
            }
            text
        };
        std::fs::write(&perl_path, rows(0)).unwrap();
        std::fs::write(&rust_path, rows(1)).unwrap();
        let outputs = [
            (perl_path.as_path(), "vep-rs"),
            (rust_path.as_path(), "vep-rs"),
        ];
        let (strings, whole) = parse_vep_outputs(&outputs, 1).unwrap();
        let whole: Vec<Vec<[&str; 5]>> = whole.iter().map(|o| texts(&strings, &o.tuples)).collect();
        assert!(whole[0].len() > 100);
        for block_bytes in [1, 33, 500] {
            for threads in [1, 2, 9] {
                let interner = Interner::new();
                let sources: Vec<OutputSource> = outputs
                    .iter()
                    .map(|(path, engine)| OutputSource {
                        block_bytes,
                        ..OutputSource::new(path.to_path_buf(), engine, &interner, None)
                    })
                    .collect();
                let rows: Vec<OutputRows> = sources
                    .iter()
                    .zip(read_sources(threads, &sources))
                    .map(|(source, (read, blocks))| source.finish(read, blocks).unwrap().unwrap())
                    .collect();
                drop(sources);
                let (strings, ranks) = interner.finish();
                let sets: Vec<TupleSet> = rows.iter().map(|r| ranks.tuple_set(&r.rows)).collect();
                let read: Vec<Vec<[&str; 5]>> = sets.iter().map(|s| texts(&strings, s)).collect();
                assert_eq!(read, whole, "{block_bytes} {threads}");
            }
        }
    }

    #[test]
    fn test_per_type_table_keys_a_ranged_record_by_its_printed_location() {
        let variants = [VcfVariant {
            chrom: "21".into(),
            pos: 38_589_780,
            vid: "synth_nonref_0188".into(),
            reference: "A".into(),
            alt: "<NON_REF>".into(),
            variant_type: "NON_REF".into(),
            location: "21:38589780".into(),
            info: String::new(),
            pos_text: None,
        }];
        let loc = "21:38589781-38590008";
        let owned = crate::strings::OwnedSides::from_rows(
            &[[loc, "<NON_REF>", "ENST_A", "Transcript", "intron_variant"]],
            &[[loc, "<NON_REF>", "-", "-", "intergenic_variant"]],
        );
        let dir = TempDir::new("output_types_").unwrap();
        let out = dir.path.join("perl.txt");
        write_vep(
            &out,
            &[[
                "synth_nonref_0188",
                loc,
                "<NON_REF>",
                "-",
                "ENST_A",
                "Transcript",
                "intron_variant",
            ]],
        );
        let types = location_types_from_outputs(&variants, &[&out]).unwrap();
        let expected: BTreeMap<String, String> = [
            ("21:38589780".to_string(), "NON_REF".to_string()),
            (loc.to_string(), "NON_REF".to_string()),
        ]
        .into();
        assert_eq!(types, expected);
        let types = LocationTypes::from_texts(&owned.strings, &types);
        let printed = [owned.strings.get(loc).unwrap()];
        let per_type = |types: Option<&LocationTypes>| {
            let rows = crate::metrics::compute_per_type_metrics(
                &owned.strings,
                &variants,
                &owned.perl,
                &printed,
                &owned.rust,
                &printed,
                types,
            );
            let (name, metrics) = rows.into_iter().next().unwrap();
            assert_eq!(name, "NON_REF");
            metrics
        };
        let without = per_type(None);
        let with_map = per_type(Some(&types));
        assert_eq!((without.perl_tuples, without.rust_tuples), (0, 0));
        assert_eq!(
            (
                with_map.perl_tuples,
                with_map.rust_tuples,
                with_map.intersection
            ),
            (1, 1, 0)
        );
        let typed = |types: Option<&LocationTypes>| -> Vec<String> {
            crate::metrics::collect_discordants(
                &owned.strings,
                "10_special_alleles",
                &variants,
                &owned.perl,
                &owned.rust,
                types,
            )
            .into_iter()
            .map(|record| record.variant_type)
            .collect()
        };
        assert_eq!(typed(Some(&types)), ["NON_REF", "NON_REF"]);
        assert_eq!(typed(None), ["Unknown", "Unknown"]);
    }

    #[test]
    fn test_a_record_without_an_id_keeps_its_position_key() {
        let variants = [VcfVariant {
            chrom: "21".into(),
            pos: 100,
            vid: ".".into(),
            reference: "A".into(),
            alt: "T".into(),
            variant_type: "SNV".into(),
            location: "21:100".into(),
            info: String::new(),
            pos_text: None,
        }];
        let dir = TempDir::new("output_no_id_").unwrap();
        let out = dir.path.join("perl.txt");
        write_vep(
            &out,
            &[[
                ".",
                "21:100",
                "T",
                "-",
                "ENST_A",
                "Transcript",
                "missense_variant",
            ]],
        );
        let types = location_types_from_outputs(&variants, &[&out]).unwrap();
        let expected: BTreeMap<String, String> = [("21:100".to_string(), "SNV".to_string())].into();
        assert_eq!(types, expected);
    }

    #[test]
    fn a_printed_location_takes_the_type_of_the_first_row_with_a_known_id() {
        let record = |vid: &str, variant_type: &str, location: &str| VcfVariant {
            chrom: "21".into(),
            pos: 1,
            vid: vid.into(),
            reference: "N".into(),
            alt: "<DEL>".into(),
            variant_type: variant_type.into(),
            location: location.into(),
            info: String::new(),
            pos_text: None,
        };
        let variants = [
            record("a", "Symbolic_DEL", "21:1"),
            record("b", "Symbolic_DUP", "21:1"),
            record("a", "BND", "21:2"),
            record("", "SNV", "21:3"),
            record(".", "MNP", "21:4"),
        ];
        let dir = TempDir::new("output_first_").unwrap();
        let (perl, rust, absent) = (
            dir.path.join("perl.txt"),
            dir.path.join("rust.txt"),
            dir.path.join("absent.txt"),
        );
        std::fs::write(
            &perl,
            "#h\nunknown\tchr21:5-9\tx\nb\t21:5-9\tx\tmore\na\t21:5-9\nshort\t21:6\n\nb\t21:1\tx\n.\t21:7\tx\n\t21:8\tx\n",
        )
        .unwrap();
        std::fs::write(&rust, "a\t21:5-9\tx\na\t21:10\tx\nb\t21:2\tx\n").unwrap();
        let types = location_types_from_outputs(&variants, &[&perl, &absent, &rust]).unwrap();
        let expected: BTreeMap<String, String> = [
            ("21:1", "Symbolic_DUP"),
            ("21:2", "BND"),
            ("21:3", "SNV"),
            ("21:4", "MNP"),
            ("21:5-9", "Symbolic_DUP"),
            ("21:10", "Symbolic_DEL"),
        ]
        .into_iter()
        .map(|(location, variant_type)| (location.to_string(), variant_type.to_string()))
        .collect();
        assert_eq!(types, expected);
    }

    #[test]
    fn a_record_with_several_identifiers_is_typed_by_its_first() {
        let dir = TempDir::new("output_ids_").unwrap();
        let input = dir.path.join("in.vcf");
        // The second record opens its ID column with the first identifier of the first.
        std::fs::write(
            &input,
            "21\t1000\tsv_one;sv_two\tN\t<DEL>\t.\t.\tEND=2000\n\
             21\t5000\tsv_one;sv_three\tN\t<DUP>\t.\t.\tEND=6000\n",
        )
        .unwrap();
        let variants = crate::vcf::parse_vcf(&input).unwrap();
        let ranged = "21:1001-2000";
        let out = dir.path.join("perl.txt");
        write_vep(
            &out,
            &[[
                "sv_one",
                ranged,
                "deletion",
                "-",
                "ENST_A",
                "Transcript",
                "intron_variant",
            ]],
        );
        let types = location_types_from_outputs(&variants, &[&out]).unwrap();
        let expected: BTreeMap<String, String> = [
            ("21:1000", "Symbolic_DEL"),
            ("21:5000", "Symbolic_DUP"),
            (ranged, "Symbolic_DEL"),
        ]
        .into_iter()
        .map(|(location, variant_type)| (location.to_string(), variant_type.to_string()))
        .collect();
        assert_eq!(types, expected);
        let owned = crate::strings::OwnedSides::from_rows(
            &[[ranged, "deletion", "ENST_A", "Transcript", "intron_variant"]],
            &[[ranged, "deletion", "-", "-", "intergenic_variant"]],
        );
        let types = LocationTypes::from_texts(&owned.strings, &types);
        let typed: Vec<String> = crate::metrics::collect_discordants(
            &owned.strings,
            "in",
            &variants,
            &owned.perl,
            &owned.rust,
            Some(&types),
        )
        .into_iter()
        .map(|record| record.variant_type)
        .collect();
        assert_eq!(typed, ["Symbolic_DEL", "Symbolic_DEL"]);
    }

    #[test]
    fn a_location_is_read_only_on_a_row_whose_id_a_record_has() {
        let variants = [VcfVariant {
            chrom: "21".into(),
            pos: 1,
            vid: "a".into(),
            reference: "N".into(),
            alt: "<DEL>".into(),
            variant_type: "Symbolic_DEL".into(),
            location: "21:1".into(),
            info: String::new(),
            pos_text: None,
        }];
        let long = "1".repeat(INT_MAX_STR_DIGITS + 1);
        let dir = TempDir::new("output_unread_").unwrap();
        let out = dir.path.join("perl.txt");
        let types = |text: String| {
            std::fs::write(&out, text).unwrap();
            location_types_from_outputs(&variants, &[&out]).map(|types| {
                types
                    .into_iter()
                    .map(|(location, variant_type)| format!("{location}={variant_type}"))
                    .collect::<Vec<String>>()
            })
        };
        // No record has the ID of the row: its coordinate, longer than an integer text may
        // be, is not read.
        assert_eq!(
            types(format!("unknown\t21:{long}\tx\na\t21:5-9\tx\n")),
            Ok(vec![
                "21:1=Symbolic_DEL".to_string(),
                "21:5-9=Symbolic_DEL".to_string()
            ])
        );
        assert_eq!(
            types(format!("\t21:{long}\tx\na\t21:5-9\tx\n")),
            Ok(vec![
                "21:1=Symbolic_DEL".to_string(),
                "21:5-9=Symbolic_DEL".to_string()
            ])
        );
        assert_eq!(
            types(format!(".\t21:{long}-7\tx\n")),
            Ok(vec!["21:1=Symbolic_DEL".to_string()])
        );
        // A row of two columns is passed over whatever its ID.
        assert_eq!(
            types(format!("unknown\t21:{long}\na\t21:{long}\n")),
            Ok(vec!["21:1=Symbolic_DEL".to_string()])
        );
        // A record has the ID: the coordinate is read, and fails.
        for text in [
            format!("a\t21:5-9\tx\na\t21:{long}\tx\n"),
            format!("a\t21:{long}\tx\tmore\n"),
        ] {
            let refused = types(text).unwrap_err();
            assert_eq!(refused.kind, FailureKind::Value);
            assert_eq!(
                refused.message,
                "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 digits; use sys.set_int_max_str_digits() to increase the limit"
            );
        }
    }

    #[test]
    fn a_tuple_row_has_its_location_read_whatever_its_id() {
        let long = "1".repeat(INT_MAX_STR_DIGITS + 1);
        let dir = TempDir::new("output_read_").unwrap();
        let out = dir.path.join("perl.txt");
        write_vep(
            &out,
            &[[
                "unknown",
                &format!("21:{long}"),
                "T",
                "-",
                "ENST_A",
                "Transcript",
                "missense_variant",
            ]],
        );
        let refused = parse_vep_output(&out, "vep-rs").unwrap_err();
        assert_eq!(refused.kind, FailureKind::Value);
        assert!(
            refused
                .message
                .starts_with("Exceeds the limit (4300 digits)"),
            "{refused}"
        );
    }

    #[test]
    fn location_types_keep_the_first_or_the_last_writer_as_asked() {
        let mut types = LocationTypes::default();
        types.set(3, "SNV");
        types.set(3, "MNP");
        types.set_default(3, "BND");
        types.set_default(4, "BND");
        assert_eq!(
            (types.get(3), types.get(4), types.get(5)),
            (Some("MNP"), Some("BND"), None)
        );
        let mut later = LocationTypes::default();
        later.set(4, "SNV");
        later.set(9, "Other");
        types.update(&later);
        let mut all: Vec<(Id, &str)> = types.iter().collect();
        all.sort_unstable();
        assert_eq!(all, [(3, "MNP"), (4, "SNV"), (9, "Other")]);
        assert_eq!(types.len(), 3);
    }
}
