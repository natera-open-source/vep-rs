// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a plugin value at a Location of a record the reference dropped.
//!
//! Ensembl VEP writes at least one row for every record it keeps, so a Location at which
//! the engine wrote a structural row and the reference none holds a record the reference
//! skipped whole, and every engine row there goes with it. The Location is given back where
//! an engine row there names no record, or a record whose name a row of the reference's
//! output carries: the reference kept that record and wrote it at another Location. A gVCF
//! variant site, which Ensembl VEP routes to its structural path and drops, is recognised
//! by its own shape and is none of these Locations.

use crate::pytext::{is_ascii_digits, strip};
use crate::rules::record::{perl_fields, perl_split, vep_rs_pos};
use crate::rules::shared::sequence_location;
use crate::rules::{Refusal, Res};
use flate2::bufread::GzDecoder;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::LazyLock;

/// The most digits CPython's `int()` converts (`sys.get_int_max_str_digits()`); a longer
/// number ends the reading.
const INT_MAX_DIGITS: usize = 4300;

const NON_REF_BLOCK: &str = "<NON_REF>";

/// The one symbolic record both engines keep but write at different Locations: Ensembl VEP
/// expands a `<CNV:TR>` to literal sequence (Parser/VCF.pm:382-431), so an engine row
/// under this allele never has a reference row at its Location.
const TANDEM_REPEAT_ALLELE: &str = "tandem_repeat";

/// A first column that names no input record.
const NO_RECORD_NAMES: [&str; 3] = ["", ".", "-"];

const READ_BLOCK: usize = 1 << 20;

/// `CHROM_POS_REF/ALT...`, the name VEP gives a record without an ID; the chromosome may
/// carry underscores, the position is the first `_digits_` group after it.
static GENERATED_NAME: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^.+?_[0-9]+_(.+)$").expect("the pattern compiles"));

/// Parser/VCF.pm:241, the ALT forms that route a line to the structural path. The word
/// class is CPython's `\w` in a text pattern: letters, numbers and the underscore.
static STRUCTURAL_ALT_LIST: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"[<\[\]][^*]+[>\]\[]|^\.[\p{L}\p{N}_]+|[\p{L}\p{N}_]+\.$")
        .expect("the pattern compiles")
});

/// The mate position of a bracket-notation breakend allele.
static BRACKET_MATE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\[\]][^:\[\]]+:([0-9]+)[\[\]]").expect("the pattern compiles"));

fn too_long(what: &str, digits: usize) -> Refusal {
    Refusal(format!(
        "{what} of {digits} digits is longer than the {INT_MAX_DIGITS} digits an integer is read from"
    ))
}

fn short_row() -> Refusal {
    Refusal("a row of fewer than three columns states no Location and no allele".to_string())
}

/// The Location with each coordinate written as its number: Ensembl VEP prints an INFO/END
/// as given, so `21:30425651-030425949` and `21:30425651-+30425949` name the span
/// `21:30425651-30425949`. A coordinate that is not in the digits `0` to `9`, and a
/// Location without `:`, stays as written.
fn canonical_location(location: &str) -> Res<String> {
    let Some((chrom, coords)) = location.split_once(':') else {
        return Ok(location.to_string());
    };
    let mut canonical = String::with_capacity(location.len());
    canonical.push_str(chrom);
    canonical.push(':');
    for (n, coordinate) in coords.split('-').enumerate() {
        if n > 0 {
            canonical.push('-');
        }
        let digits = coordinate.strip_prefix('+').unwrap_or(coordinate);
        if is_ascii_digits(digits) {
            if digits.len() > INT_MAX_DIGITS {
                return Err(too_long("a Location coordinate", digits.len()));
            }
            let bare = digits.trim_start_matches('0');
            canonical.push_str(if bare.is_empty() { "0" } else { bare });
        } else {
            canonical.push_str(coordinate);
        }
    }
    Ok(canonical)
}

/// The lines of rows held one per line: a row's trailing `\n`s end its last line, a `\n`
/// inside a row ends a line and a `\r` ends none.
fn written_lines(rows: &[String]) -> impl Iterator<Item = &str> {
    rows.iter()
        .flat_map(|row| row.trim_end_matches('\n').split('\n'))
}

/// The lines of `text` as CPython's text mode reads a file with universal newlines: `\n`,
/// `\r` and `\r\n` each end a line.
fn universal_lines(text: &str) -> impl Iterator<Item = &str> {
    let mut rest = Some(text);
    std::iter::from_fn(move || {
        let text = rest?;
        let Some(at) = text.find(['\n', '\r']) else {
            rest = None;
            return Some(text);
        };
        let after = &text[at + 1..];
        rest = Some(if text.as_bytes()[at] == b'\r' {
            after.strip_prefix('\n').unwrap_or(after)
        } else {
            after
        });
        Some(&text[..at])
    })
}

fn is_non_ref_block(alt: &str) -> bool {
    alt.eq_ignore_ascii_case(NON_REF_BLOCK)
}

/// Whether an ALT names a variant: not the `<NON_REF>` block, `*` or `.`.
fn names_an_allele(alt: &str) -> bool {
    !is_non_ref_block(alt) && alt != "*" && alt != "."
}

/// Whether Perl reads an INFO value as true: defined, not empty and not `0`.
fn perl_true(value: Option<&str>) -> bool {
    !matches!(value, None | Some("") | Some("0"))
}

/// The value of `key=value` in an INFO column, the last occurrence winning as in the hash
/// Ensembl VEP's parser builds (BaseVCF4.pm:581-583), the text up to a second `=` kept;
/// `None` for an absent key or a bare flag.
fn info_written<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    let mut value = None;
    for item in info.split(';') {
        match item.split_once('=') {
            Some((name, text)) if name == key => {
                value = Some(text.split_once('=').map_or(text, |(kept, _)| kept));
            }
            None if item == key => value = None,
            _ => {}
        }
    }
    value
}

/// Whether a record is a gVCF variant site vep-rs reads on its sequence path: `<NON_REF>`
/// beside alleles that name a variant, no true INFO/SVTYPE, and no remaining allele written
/// in a structural form (Parser/VCF.pm:237-245). A list of `<NON_REF>` alone is a reference
/// block, and a block beside a symbolic allele is an unsupported structural list.
fn is_gvcf_variant_site(alts: &[&str], svtype: Option<&str>) -> bool {
    let alts = perl_fields(alts.to_vec());
    let real: Vec<&str> = alts
        .iter()
        .copied()
        .filter(|alt| names_an_allele(alt))
        .collect();
    if real.is_empty() || !alts.iter().any(|alt| is_non_ref_block(alt)) || perl_true(svtype) {
        return false;
    }
    !STRUCTURAL_ALT_LIST.is_match(&real.join(","))
}

fn same_letter(a: char, b: char) -> bool {
    a.to_uppercase().eq(b.to_uppercase())
}

/// A bi-allelic indel minimised as vep-rs minimises it: the shared prefix trimmed while
/// both alleles have bases, then the shared suffix, an emptied allele reading `-`; returns
/// the two alleles and the bases trimmed in front.
fn minimised(reference: &[char], alt: &[char]) -> (Vec<char>, Vec<char>, usize) {
    let mut shared = 0;
    while shared < reference.len().min(alt.len()) && same_letter(reference[shared], alt[shared]) {
        shared += 1;
    }
    let (mut reference, mut alt) = (&reference[shared..], &alt[shared..]);
    while let (Some(r), Some(a)) = (reference.last(), alt.last()) {
        if !same_letter(*r, *a) {
            break;
        }
        reference = &reference[..reference.len() - 1];
        alt = &alt[..alt.len() - 1];
    }
    let or_dash = |bases: &[char]| {
        if bases.is_empty() {
            vec!['-']
        } else {
            bases.to_vec()
        }
    };
    (or_dash(reference), or_dash(alt), shared)
}

/// The (allele, start, end) vep-rs writes for each real allele of a gVCF variant site, the
/// allele in upper case, on the sequence path of its parser: when any ALT's length differs
/// from REF's, the first base REF shares with every compared ALT is trimmed and the start
/// moves one base (Parser/VCF.pm:296-328, an ALT carrying a `*` kept whole and, beside other
/// ALTs, left out of the comparison); a single remaining allele whose length still differs
/// from REF's is minimised from both ends (Parser.pm:978-984); an insertion's end is its
/// start minus one.
pub(super) fn gvcf_site_rows(
    pos: i128,
    reference: &str,
    alts: &[&str],
) -> Vec<(String, i128, i128)> {
    let alts = perl_fields(alts.to_vec());
    let named: Vec<&str> = alts
        .iter()
        .copied()
        .filter(|alt| !is_non_ref_block(alt))
        .collect();
    let reference: Vec<char> = reference.chars().collect();
    let is_indel = named
        .iter()
        .any(|alt| alt.chars().count() != reference.len());
    let compared: Vec<&str> = if named.len() == 1 {
        named.clone()
    } else {
        named
            .iter()
            .copied()
            .filter(|alt| !alt.contains('*'))
            .collect()
    };
    let anchored = !reference.is_empty()
        && compared.iter().all(|alt| {
            alt.chars()
                .next()
                .is_some_and(|first| same_letter(first, reference[0]))
        });
    let prefix = usize::from(is_indel && anchored);
    let trimmed_reference: Vec<char> = if reference.len() > prefix {
        reference[prefix..].to_vec()
    } else {
        vec!['-']
    };
    let mut rows = Vec::new();
    for alt in named.iter().copied().filter(|alt| names_an_allele(alt)) {
        let alt_bases: Vec<char> = alt.chars().collect();
        let mut allele: Vec<char> = if alt.contains('*') {
            alt_bases
        } else if alt_bases.len() > prefix {
            alt_bases[prefix..].to_vec()
        } else {
            vec!['-']
        };
        let mut start = pos + prefix as i128;
        let mut reference_part = trimmed_reference.clone();
        if named.len() == 1
            && reference_part != ['-']
            && allele != ['-']
            && reference_part.len() != allele.len()
        {
            let shared;
            (reference_part, allele, shared) = minimised(&reference_part, &allele);
            start += shared as i128;
        }
        let end = if reference_part == ['-'] {
            start - 1
        } else {
            start + reference_part.len() as i128 - 1
        };
        let allele: String = allele.into_iter().collect();
        rows.push((allele.to_uppercase(), start, end));
    }
    rows
}

/// The POS column as vep-rs's parsers read it, `None` where vep-rs fails the line; an
/// unpadded integer of more digits than `int()` converts ends the reading.
fn pos_of(text: &str) -> Res<Option<i128>> {
    let digits = text.strip_prefix('+').unwrap_or(text);
    if is_ascii_digits(digits) && digits.len() > INT_MAX_DIGITS {
        return Err(too_long("a POS", digits.len()));
    }
    Ok(vep_rs_pos(text))
}

/// Adds the Locations of the engine's rows whose record is a gVCF variant site by its
/// name: a record without an ID is named `CHROM_POS_REF/ALT...` from its own line, so a
/// name whose ALT list carries `<NON_REF>` beside another allele is the site.
fn add_sites_of_generated_names<'a>(
    lines: impl Iterator<Item = &'a str>,
    sites: &mut HashSet<String>,
) -> Res<()> {
    for line in lines {
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let (Some(name), Some(location)) = (fields.next(), fields.next()) else {
            continue;
        };
        if fields.count() < 5 {
            continue;
        }
        let Some(named) = GENERATED_NAME.captures(name) else {
            continue;
        };
        let alleles = perl_split(&named[1], "/");
        if is_gvcf_variant_site(alleles.get(1..).unwrap_or(&[]), None) {
            sites.insert(canonical_location(strip(location))?);
        }
    }
    Ok(())
}

/// Adds the Locations vep-rs writes for every gVCF variant site among the input's record
/// lines, recognised by the record's shape whatever its ID, at the POS as vep-rs reads it;
/// a record whose POS vep-rs refuses has no rows.
fn add_sites_of_records(input_lines: &[String], sites: &mut HashSet<String>) -> Res<()> {
    for line in input_lines.iter().flat_map(|line| universal_lines(line)) {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let (chrom, reference) = (fields[0], fields[3]);
        let pos = pos_of(fields[1])?;
        let alts = perl_split(fields[4], ",");
        let svtype = fields.get(7).and_then(|info| info_written(info, "SVTYPE"));
        let Some(pos) = pos else {
            continue;
        };
        if !is_gvcf_variant_site(&alts, svtype) {
            continue;
        }
        for (_, start, end) in gvcf_site_rows(pos, reference, &alts) {
            sites.insert(canonical_location(&sequence_location(chrom, start, end))?);
        }
    }
    Ok(())
}

/// Whether an allele column names a structural allele: VEP writes a sequence allele as its
/// bases or `-` and a symbolic one as its Sequence Ontology class, so anything outside `-`
/// and the nucleotide alphabet is structural.
fn is_structural_allele(allele: &str) -> bool {
    allele != "-"
        && !allele
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// The Locations at which the engine wrote a structural row and the reference none, read
/// from the rows in the order they are handed: a run of rows at one Location is one group,
/// which closes when the Location changes. A gVCF variant site is none of them, whatever
/// its rows' alleles.
struct Scope<'a> {
    gvcf_sites: &'a HashSet<String>,
    absent: HashSet<String>,
    group: Option<String>,
    reference_wrote: bool,
    engine_wrote_structural: bool,
}

impl<'a> Scope<'a> {
    fn new(gvcf_sites: &'a HashSet<String>) -> Scope<'a> {
        Scope {
            gvcf_sites,
            absent: HashSet::new(),
            group: None,
            reference_wrote: false,
            engine_wrote_structural: false,
        }
    }

    /// Observes one row by its Location and allele and the side that wrote it.
    fn observe(&mut self, location: &str, allele: &str, by_reference: bool) -> Res<()> {
        if self.group.as_deref() != Some(location) {
            self.close();
            self.group = Some(location.to_string());
        }
        if by_reference {
            self.reference_wrote = true;
        } else if self.gvcf_sites.contains(location)
            || (is_structural_allele(allele) && allele != TANDEM_REPEAT_ALLELE)
        {
            self.engine_wrote_structural = true;
        }
        // The mate position of a breakend allele is read as an integer.
        if let Some(mate) = BRACKET_MATE.captures(allele) {
            if mate[1].len() > INT_MAX_DIGITS {
                return Err(too_long("a breakend mate position", mate[1].len()));
            }
        }
        Ok(())
    }

    /// Closes the open group.
    fn close(&mut self) {
        if let Some(location) = self.group.take() {
            if self.engine_wrote_structural
                && !self.reference_wrote
                && !self.gvcf_sites.contains(&location)
            {
                self.absent.insert(location);
            }
        }
        self.reference_wrote = false;
        self.engine_wrote_structural = false;
    }
}

/// The record name and the Location of a default-format row that states a Location, an
/// allele and a Feature in seven columns or more: the name is the first column, the input
/// record's ID or the name the engine generates for a record without one. With `names`,
/// only a row whose name is among them, asked before its Location is read.
fn named_row<'a>(line: &'a str, names: Option<&HashSet<&str>>) -> Res<Option<(&'a str, String)>> {
    if line.starts_with('#') {
        return Ok(None);
    }
    let mut fields = line.splitn(7, '\t');
    let mut columns = [""; 7];
    for column in &mut columns {
        match fields.next() {
            Some(field) => *column = field,
            None => return Ok(None),
        }
    }
    let name = strip(columns[0]);
    if names.is_some_and(|names| !names.contains(name)) {
        return Ok(None);
    }
    let location = canonical_location(strip(columns[1]))?;
    if location.is_empty() || strip(columns[2]).is_empty() || strip(columns[4]).is_empty() {
        return Ok(None);
    }
    Ok(Some((name, location)))
}

/// The members of a gzip file one after the other; zero bytes after a member are padding,
/// as CPython's `gzip` reads them.
struct Members<R: BufRead> {
    member: Option<GzDecoder<R>>,
}

impl<R: BufRead> Read for Members<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        loop {
            let Some(member) = self.member.as_mut() else {
                return Ok(0);
            };
            let count = member.read(buffer)?;
            if count > 0 || buffer.is_empty() {
                return Ok(count);
            }
            let mut stored = self.member.take().expect("an open member").into_inner();
            loop {
                let padding = stored.fill_buf()?;
                let zeros = padding.iter().take_while(|&&byte| byte == 0).count();
                if zeros == 0 {
                    break;
                }
                stored.consume(zeros);
            }
            if stored.fill_buf()?.is_empty() {
                return Ok(0);
            }
            self.member = Some(GzDecoder::new(stored));
        }
    }
}

/// Calls `take` with each line of `reader` without its end, as CPython's text mode reads a
/// file with universal newlines: `\n`, `\r` and `\r\n` each end a line, and a last line
/// needs no end. A failed read is refused as `unread` words it.
fn for_each_line(
    mut reader: impl Read,
    unread: &dyn Fn(io::Error) -> Refusal,
    mut take: impl FnMut(&[u8]) -> Res<()>,
) -> Res<()> {
    let mut block = vec![0u8; READ_BLOCK];
    let (mut start, mut end, mut ended) = (0, 0, false);
    loop {
        while let Some(found) = memchr::memchr2(b'\n', b'\r', &block[start..end]) {
            let at = start + found;
            // A `\r` that ends the block may be the first half of a `\r\n`.
            if block[at] == b'\r' && at + 1 == end && !ended {
                break;
            }
            take(&block[start..at])?;
            start = at + 1;
            if block[at] == b'\r' && start < end && block[start] == b'\n' {
                start += 1;
            }
        }
        if ended {
            if start < end {
                take(&block[start..end])?;
            }
            return Ok(());
        }
        if start > 0 {
            block.copy_within(start..end, 0);
            end -= start;
            start = 0;
        }
        if end == block.len() {
            block.resize(block.len() * 2, 0);
        }
        match reader.read(&mut block[end..]) {
            Ok(0) => ended = true,
            Ok(count) => end += count,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(unread(error)),
        }
    }
}

/// Those of `wanted` that a row of the stored output `reference`, plain or gzip, carries as
/// its record name; its bytes are read as UTF-8, a malformed sequence as U+FFFD.
fn names_written<'a>(reference: &Path, wanted: &HashSet<&'a str>) -> Res<HashSet<&'a str>> {
    let shown = reference.display();
    let unread = |error: io::Error| Refusal(format!("{shown}: {error}"));
    let broken =
        |error: io::Error| Refusal(format!("{shown}: cannot read the gzip stream: {error}"));
    let stored = File::open(reference).map_err(unread)?;
    let mut magic = [0u8; 2];
    let gzip = stored.read_at(&mut magic, 0).map_err(unread)? == 2 && magic == [0x1f, 0x8b];
    let mut written = HashSet::new();
    let take = |line: &[u8]| -> Res<()> {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let Some(tab) = memchr::memchr(b'\t', line) else {
            return Ok(());
        };
        let Some(&name) = wanted.get(strip(&String::from_utf8_lossy(&line[..tab]))) else {
            return Ok(());
        };
        if named_row(&String::from_utf8_lossy(line), None)?.is_some() {
            written.insert(name);
        }
        Ok(())
    };
    if gzip {
        let members = Members {
            member: Some(GzDecoder::new(BufReader::with_capacity(READ_BLOCK, stored))),
        };
        for_each_line(members, &broken, take)?;
    } else {
        for_each_line(stored, &unread, take)?;
    }
    Ok(written)
}

/// Of `locations`, each a Location at which the engine wrote and the reference nothing,
/// those that are not a record the reference dropped: an engine row there carries the name
/// of a record the reference wrote a row for, or names no record. Both engines print the
/// input record's name in the first column, so the reference's output tells a record it
/// dropped, whose name is on none of its rows, from a record it wrote at another Location.
/// Two records at one Location stand or fall together. The reference's output is read
/// once, and only when a Location holds a name to look for.
fn reference_named_locations(
    reference: &Path,
    engine_rows: &[String],
    locations: &HashSet<String>,
) -> Res<HashSet<String>> {
    if locations.is_empty() {
        return Ok(HashSet::new());
    }
    let mut names_at: HashMap<String, HashSet<&str>> = HashMap::new();
    for line in written_lines(engine_rows) {
        if let Some((name, location)) = named_row(line, None)? {
            if locations.contains(&location) {
                names_at.entry(location).or_default().insert(name);
            }
        }
    }
    let wanted: HashSet<&str> = names_at
        .values()
        .flatten()
        .copied()
        .filter(|name| !NO_RECORD_NAMES.contains(name))
        .collect();
    let written = if wanted.is_empty() {
        HashSet::new()
    } else {
        names_written(reference, &wanted)?
    };
    Ok(names_at
        .into_iter()
        .filter(|(_, names)| {
            names
                .iter()
                .any(|name| written.contains(name) || NO_RECORD_NAMES.contains(name))
        })
        .map(|(location, _)| location)
        .collect())
}

/// Whether the rows the two sides wrote at one Location, each a default-format row as
/// written, are those of a record the reference dropped: the engine wrote a structural row
/// there and the reference none, and no engine row there names no record or one that the
/// reference's stored output `reference` carries on a row. The gVCF variant sites that the
/// engine's generated names and the record lines `input_lines` show are never such a
/// Location. `reference` is read only for a Location the first test took.
pub fn dropped_by_the_reference(
    reference_rows: &[String],
    engine_rows: &[String],
    reference: &Path,
    input_lines: &[String],
) -> Result<bool, Refusal> {
    let mut sites = HashSet::new();
    add_sites_of_generated_names(written_lines(engine_rows), &mut sites)?;
    add_sites_of_records(input_lines, &mut sites)?;
    let mut scope = Scope::new(&sites);
    let mut locations = HashSet::new();
    for (rows, by_reference) in [(reference_rows, true), (engine_rows, false)] {
        for row in rows {
            let mut columns = row.splitn(4, '\t').skip(1);
            let location = canonical_location(strip(columns.next().ok_or_else(short_row)?))?;
            let allele = strip(columns.next().ok_or_else(short_row)?);
            scope.observe(&location, allele, by_reference)?;
            locations.insert(location);
        }
    }
    scope.close();
    let mut absent = scope.absent;
    for named in reference_named_locations(reference, engine_rows, &absent)? {
        absent.remove(&named);
    }
    Ok(!locations.is_empty() && locations == absent)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(pos: i128, reference: &str, alts: &[&str]) -> Vec<(i128, i128)> {
        gvcf_site_rows(pos, reference, alts)
            .into_iter()
            .map(|(_, start, end)| (start, end))
            .collect()
    }

    /// A reader that hands its text over a few bytes at a time.
    struct Chunks<'a>(Vec<&'a [u8]>);

    impl Read for Chunks<'_> {
        fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
            if self.0.is_empty() {
                return Ok(0);
            }
            let chunk = self.0.remove(0);
            buffer[..chunk.len()].copy_from_slice(chunk);
            Ok(chunk.len())
        }
    }

    #[test]
    fn a_location_is_written_with_each_coordinate_as_its_number() {
        for (written, canonical) in [
            ("21:005-+9", "21:5-9"),
            ("21:30425651-030425949", "21:30425651-30425949"),
            ("21:30425651-+30425949", "21:30425651-30425949"),
            ("21:000", "21:0"),
            ("21:-4-295", "21:-4-295"),
            ("21:30425949.5", "21:30425949.5"),
            ("21:10\u{b2}", "21:10\u{b2}"),
            ("21:++5", "21:++5"),
            ("21", "21"),
            ("", ""),
            ("HLA-A*01:01:01:007", "HLA-A*01:01:01:007"),
            ("a:01:5-007", "a:01:5-7"),
        ] {
            assert_eq!(canonical_location(written).unwrap(), canonical, "{written}");
        }
    }

    #[test]
    fn a_coordinate_past_the_digits_an_integer_is_read_from_is_refused() {
        let held = format!("21:{}", "9".repeat(INT_MAX_DIGITS));
        assert_eq!(canonical_location(&held).unwrap(), held);
        let zeros = format!("21:+{}", "0".repeat(INT_MAX_DIGITS));
        assert_eq!(canonical_location(&zeros).unwrap(), "21:0");
        for long in ["9", "0"] {
            let refused = canonical_location(&format!("21:5-{}", long.repeat(INT_MAX_DIGITS + 1)));
            assert_eq!(
                refused,
                Err(Refusal(
                    "a Location coordinate of 4301 digits is longer than the 4300 digits an integer is read from"
                        .to_string()
                ))
            );
        }
        assert_eq!(pos_of(&"1".repeat(INT_MAX_DIGITS)), Ok(None));
        assert_eq!(
            pos_of(&format!("+{}", "0".repeat(INT_MAX_DIGITS - 1))),
            Ok(Some(0))
        );
        assert_eq!(pos_of(" 30425650 "), Ok(Some(30425650)));
        for long in ["1", "0"] {
            let digits = long.repeat(INT_MAX_DIGITS + 1);
            assert!(pos_of(&digits).is_err());
            assert!(pos_of(&format!("+{digits}")).is_err());
            assert_eq!(
                pos_of(&format!(" {digits}")),
                Ok((long == "0").then_some(0))
            );
            assert_eq!(
                pos_of(&format!("{digits}.0")),
                Ok((long == "0").then_some(0))
            );
        }
        let sites = HashSet::new();
        let mut scope = Scope::new(&sites);
        let held = format!("N[chr5:{}[", "1".repeat(INT_MAX_DIGITS));
        assert_eq!(scope.observe("21:5", &held, false), Ok(()));
        let long = format!("N[chr5:{}[", "1".repeat(INT_MAX_DIGITS + 1));
        assert!(scope.observe("21:5", &long, true).is_err());
    }

    #[test]
    fn a_non_ref_block_beside_a_real_allele_is_a_gvcf_variant_site() {
        assert!(is_gvcf_variant_site(&["G", "<NON_REF>"], None));
        assert!(is_gvcf_variant_site(&["<non_ref>", "AT"], None));
        assert!(!is_gvcf_variant_site(&["<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["<NON_REF>", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["*", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&[".", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["G", "T"], None));
        assert!(!is_gvcf_variant_site(&[], None));
        assert!(!is_gvcf_variant_site(&["<DEL>", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["G", "<DEL>", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["G.", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&[".G", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["G]21:5]", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["G", "<NON_REF>"], Some("DEL")));
        assert!(is_gvcf_variant_site(&["G", "<NON_REF>"], None));
        assert!(is_gvcf_variant_site(&["G", "<NON_REF>"], Some("0")));
        assert!(is_gvcf_variant_site(&["G", "<NON_REF>"], Some("")));
    }

    #[test]
    fn a_trailing_empty_alt_is_no_field() {
        assert!(!is_gvcf_variant_site(&["G", ""], None));
        assert!(is_gvcf_variant_site(&["G", "<NON_REF>", ""], None));
        assert!(!is_gvcf_variant_site(&["<NON_REF>", ""], None));
        assert_eq!(spans(100, "A", &["G", ""]), [(100, 100)]);
        assert_eq!(spans(100, "A", &["G", "<NON_REF>", ""]), [(100, 100)]);
    }

    #[test]
    fn the_word_class_of_a_structural_alt_is_letters_numbers_and_the_underscore() {
        // A superscript digit is a number and a combining mark is neither.
        assert!(!is_gvcf_variant_site(&[".\u{b2}", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["\u{e9}.", "<NON_REF>"], None));
        assert!(!is_gvcf_variant_site(&["._", "<NON_REF>"], None));
        assert!(is_gvcf_variant_site(&[".\u{301}", "<NON_REF>"], None));
        assert!(is_gvcf_variant_site(&["\u{203f}.", "<NON_REF>"], None));
    }

    #[test]
    fn an_info_value_is_the_last_written_up_to_a_second_equals() {
        assert_eq!(info_written("SVTYPE==DEL;END=30425949", "SVTYPE"), Some(""));
        assert_eq!(
            info_written("SVTYPE==DEL;END=30425949", "END"),
            Some("30425949")
        );
        assert_eq!(
            info_written("SVTYPE=BND;CHR2=21;END2=30425700=5", "END2"),
            Some("30425700")
        );
        assert_eq!(info_written("END2=30425700=5=7", "END2"), Some("30425700"));
        assert_eq!(info_written("SVTYPE==DEL;END", "END"), None);
        assert_eq!(info_written("END=5;END", "END"), None);
        assert_eq!(info_written("END=5;END=7", "END"), Some("7"));
        assert_eq!(info_written("", "END"), None);
    }

    #[test]
    fn the_spans_of_a_site_follow_the_sequence_path_of_the_parser() {
        assert_eq!(spans(100, "A", &["G", "<NON_REF>"]), [(100, 100)]);
        assert_eq!(spans(100, "a", &["g", "<non_ref>"]), [(100, 100)]);
        assert_eq!(spans(100, "AT", &["GC", "<NON_REF>"]), [(100, 101)]);
        assert_eq!(spans(100, "A", &["AT", "<NON_REF>"]), [(101, 100)]);
        assert_eq!(spans(100, "N", &["NA", "<NON_REF>"]), [(101, 100)]);
        assert_eq!(spans(100, "ATG", &["A", "<NON_REF>"]), [(101, 102)]);
        assert_eq!(spans(100, "ATG", &["AG", "<NON_REF>"]), [(101, 101)]);
        assert_eq!(spans(100, "atg", &["ag", "<NON_REF>"]), [(101, 101)]);
        assert_eq!(
            spans(100, "A", &["AT", "ATT", "<NON_REF>"]),
            [(101, 100), (101, 100)]
        );
        assert_eq!(
            spans(100, "A", &["G", "T", "<NON_REF>"]),
            [(100, 100), (100, 100)]
        );
        assert_eq!(spans(100, "A", &["<NON_REF>", "<NON_REF>"]), []);
        assert_eq!(spans(100, "A", &["*", "G", "<NON_REF>"]), [(100, 100)]);
    }

    #[test]
    fn a_star_or_a_dot_beside_an_indel_keeps_the_anchor_trimmed_alleles() {
        assert_eq!(spans(100, "ATT", &["AT", "*", "<NON_REF>"]), [(101, 102)]);
        assert_eq!(spans(100, "AT", &["ATT", "*", "<NON_REF>"]), [(101, 101)]);
        assert_eq!(spans(100, "ATTG", &["ATC", "*", "<NON_REF>"]), [(101, 103)]);
        assert_eq!(spans(100, "ATT", &["AT", ".", "<NON_REF>"]), [(100, 102)]);
        assert_eq!(spans(100, "ATT", &["AT", "<NON_REF>"]), [(102, 102)]);
        assert_eq!(spans(100, "AT", &["ATT", "<NON_REF>"]), [(102, 101)]);
        assert_eq!(spans(100, "ATTG", &["ATC", "<NON_REF>"]), [(102, 103)]);
    }

    #[test]
    fn a_site_is_read_from_a_generated_name_and_from_a_record_line() {
        let named = |name: &str, location: &str| {
            format!("{name}\t{location}\tG\tENSG1\tENST1\tTranscript\tintron_variant\t-")
        };
        let rows = [
            named("21_30425650_A/G/<NON_REF>", " 21:030425650 "),
            named("21_30425700_A/<NON_REF>", "21:30425700"),
            named("21_30425800_A/<DEL>/<NON_REF>", "21:30425800"),
            named("rs1", "21:30425900"),
            named("chr_Un_5_A/T/<non_ref>", "chr_Un:5"),
            "21_7_A/G/<NON_REF>\t21:7\tG\tENSG1\tENST1\tTranscript".to_string(),
            format!("#{}", named("21_8_A/G/<NON_REF>", "21:8")),
        ];
        let mut sites = HashSet::new();
        add_sites_of_generated_names(written_lines(&rows), &mut sites).unwrap();
        let expected: HashSet<String> = ["21:30425650", "chr_Un:5"].map(String::from).into();
        assert_eq!(sites, expected);

        let records = [
            "#CHROM\tPOS\tID\tREF\tALT".to_string(),
            "21\t30425650\trs1\tA\tAT,<NON_REF>\t.\t.\tEND=5".to_string(),
            "chr22\t100\t.\tA\tG,<NON_REF>".to_string(),
            "21\t200\t.\tA\tG,<NON_REF>\t.\t.\tSVTYPE=DEL".to_string(),
            "21\t300\t.\tA\tG,<NON_REF>\t.\t.\tSVTYPE=0".to_string(),
            "21\t.\t.\tA\tG,<NON_REF>".to_string(),
            "21\t400\t.\tA\t<NON_REF>".to_string(),
            "21\t500\t.\tATG\tAG,<NON_REF>\r".to_string(),
            "21\t600\t.\tA".to_string(),
        ];
        let mut sites = HashSet::new();
        add_sites_of_records(&records, &mut sites).unwrap();
        let expected: HashSet<String> = ["21:30425650-30425651", "chr22:100", "21:300", "21:501"]
            .map(String::from)
            .into();
        assert_eq!(sites, expected);
    }

    #[test]
    fn lines_end_as_their_reader_ends_them() {
        let rows = [
            "a\tb\n".to_string(),
            "c\rd\ne\n\n".to_string(),
            String::new(),
        ];
        assert_eq!(
            written_lines(&rows).collect::<Vec<_>>(),
            ["a\tb", "c\rd", "e", ""]
        );
        assert_eq!(
            universal_lines("a\rb\r\nc\n\nd").collect::<Vec<_>>(),
            ["a", "b", "c", "", "d"]
        );
        for size in [1usize, 2, 3, 7] {
            let mut lines = Vec::new();
            for_each_line(
                Chunks(b"a\r\nb\rc\n\r\rd".chunks(size).collect()),
                &|error| Refusal(error.to_string()),
                |line| {
                    lines.push(String::from_utf8_lossy(line).into_owned());
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(lines, ["a", "b", "c", "", "", "d"], "reads of {size}");
        }
    }

    #[test]
    fn a_structural_allele_is_anything_outside_the_bases_and_the_dash() {
        for allele in [
            "deletion",
            "<NON_REF>",
            "N[chr5:100[",
            "A.",
            "tandem_repeat",
            "AC-",
        ] {
            assert!(is_structural_allele(allele), "{allele}");
        }
        for allele in ["-", "A", "acgtn", "ACGTN", ""] {
            assert!(!is_structural_allele(allele), "{allele}");
        }
    }

    #[test]
    fn a_run_of_rows_at_one_location_is_one_group() {
        let sites: HashSet<String> = ["21:9".to_string()].into();
        let mut scope = Scope::new(&sites);
        for (location, allele, by_reference) in [
            ("21:1", "deletion", true),
            ("21:2", "deletion", true),
            ("21:1", "deletion", false),
            ("21:2", "T", false),
            ("21:2", "duplication", false),
            ("21:3", "tandem_repeat", false),
            ("21:9", "T", false),
            ("21:4", "insertion", false),
        ] {
            scope.observe(location, allele, by_reference).unwrap();
        }
        scope.close();
        let expected: HashSet<String> = ["21:1", "21:2", "21:4"].map(String::from).into();
        assert_eq!(scope.absent, expected);
    }
}
