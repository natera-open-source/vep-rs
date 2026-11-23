// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! List items, order-only fields and declared columns.

use crate::argparse::path_text;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::{Error, Result};
use crate::input::byte_stream;
use crate::parse::json::{json_compact, json_scalar};
use crate::parse::vcf::csq_format_fields;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::report::strerror;
use crate::shared::ABSENT;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// The NearestExonJB plugin lists its exon items in the order of a hash's keys as well
/// (NearestExonJB.pm:267-271). `CLIN_SIG` under `--clin_sig_allele 1` holds each different
/// term list of a row's known variants once, the lists joined by `;` in the order of a
/// hash's keys (OutputFactory.pm:1120, 1143-1144); under `order_only` its lists agree as a
/// multiset, each list compared whole (`clin_sig_list_separator`). JSON has no such cell:
/// each known variant carries its own terms.
pub const ORDER_ONLY_COLUMNS: &[&str] = &["DOMAINS", "FREQS", "ZYG", "NearestExonJB", "CLIN_SIG"];

/// The UTRAnnotator plugin writes its `5UTR_annotation` as groups of `key=value` pairs
/// joined by `:`, the pairs in hash order within each group (UTRAnnotator.pm:282-288); under
/// `order_only` the groups agree as a multiset of groups whose pairs agree as a multiset.
pub const ORDER_ONLY_PAIR_GROUP_COLUMN: &str = "5UTR_annotation";

/// The separator of the term lists of `CLIN_SIG` in the row format `fmt`
/// (`clin_sig_allele_list_order`): the `;` Ensembl VEP joins them with, which the default and
/// VCF formats write as `%3B` (VEP_output.pm:118, VCF.pm:402).
pub fn clin_sig_list_separator(fmt: &str) -> &'static str {
    if fmt == "default" || fmt == "vcf" {
        "%3B"
    } else {
        ";"
    }
}

/// The marker `order_only_fields` gives a JSON array field in place of a list separator, and
/// the prefix it gives the pair-group column's separator.
pub const JSON_ARRAY: &str = "json";

pub const PAIR_GROUPS: &str = "pairs:";

/// The JSON field of each order-only column with its separator: the `domains` array of a
/// consequence entry, the record entry's `zyg` array and the consequence's `nearestexonjb`
/// string, whose items `json_scalar` joins with `,`; FREQS has no JSON key.
pub const ORDER_ONLY_JSON_FIELDS: &[(&str, &str)] =
    &[("domains", "json"), ("zyg", ","), ("nearestexonjb", ",")];

/// `json.loads("[" + text + "]")`: the items, or `None` for a text that is no such list.
fn loaded_items(text: &str) -> Option<Vec<Value>> {
    match pyjson::loads(&format!("[{text}]")) {
        Ok(items) => Some(items.into_list().unwrap_or_default()),
        Err(error) if error.is_value_error() => None,
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// `text.split(separator)`, which refuses an empty separator.
fn parts<'a>(text: &'a str, separator: &'a str) -> impl Iterator<Item = &'a str> {
    if separator.is_empty() {
        crate::uncaught!("empty separator");
    }
    text.split(separator)
}

/// Whether two lists hold the same items the same number of times.
fn same_multiset<T: Ord>(mut expected: Vec<T>, mut actual: Vec<T>) -> bool {
    expected.sort_unstable();
    actual.sort_unstable();
    expected == actual
}

/// The items of a list-valued field: the `separator`-joined parts of a text format's
/// value, or the elements of the JSON array `json_scalar` joined (`JSON_ARRAY`), each
/// as its compact text; under a `PAIR_GROUPS` separator the parts are groups of `:`-joined
/// pairs, each group read with its pairs sorted; `None` when the value is not such a list.
pub fn list_items(value: &str, separator: &str) -> Option<Vec<String>> {
    if let Some(group_separator) = separator.strip_prefix(PAIR_GROUPS) {
        return Some(
            parts(value, group_separator)
                .map(|group| {
                    let mut pairs: Vec<&str> = group.split(':').collect();
                    pairs.sort_unstable();
                    pairs.join(":")
                })
                .collect(),
        );
    }
    if separator != JSON_ARRAY {
        return Some(parts(value, separator).map(str::to_string).collect());
    }
    Some(loaded_items(value)?.iter().map(json_compact).collect())
}

/// Whether two JSON `domains` arrays, each object as its compact text, can be the arrays of
/// one `DOMAINS` list written in two orders (`domains_order_region_file_copy`). Ensembl
/// VEP's JSON writer re-splits each entry of the list with `(\w+):(\w+)` and pushes the two
/// captures (JSON.pm:363-371); an entry the pattern does not match pushes the captures of
/// the last entry that matched once more, and nothing where none has matched, so the order
/// of the list decides which object such an entry repeats. Every run of one object
/// therefore opens with a matched entry, and its further objects are matched entries of the
/// same captures or repeats: one list gives both arrays exactly where each side holds every
/// object at least as many times as the other side has runs of it. Two arrays of the same
/// objects the same number of times are such a pair.
pub fn arrays_of_one_domains_list(expected: &[String], actual: &[String]) -> bool {
    fn runs_and_copies(objects: &[String]) -> BTreeMap<&str, (usize, usize)> {
        let mut seen: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
        let mut previous = None;
        for object in objects {
            let (runs, copies) = seen.entry(object.as_str()).or_default();
            if previous != Some(object) {
                *runs += 1;
            }
            *copies += 1;
            previous = Some(object);
        }
        seen
    }
    let (expected, actual) = (runs_and_copies(expected), runs_and_copies(actual));
    expected.len() == actual.len()
        && expected.iter().all(|(object, (runs, copies))| {
            actual
                .get(object)
                .is_some_and(|(other_runs, other_copies)| {
                    runs <= other_copies && other_runs <= copies
                })
        })
}

/// Whether two list-valued field texts carry the same items the same number of times,
/// in any order, or, under `JSON_ARRAY`, the marker of the `domains` array, are the arrays
/// of one list (`arrays_of_one_domains_list`); an empty or absent value agrees with itself
/// alone.
pub fn same_items(expected: &str, actual: &str, separator: &str) -> bool {
    let empty = |text: &str| text.is_empty() || text == "-" || text == ABSENT;
    if empty(expected) || empty(actual) {
        return expected == actual;
    }
    let expected_items = list_items(expected, separator);
    let actual_items = list_items(actual, separator);
    match (expected_items, actual_items) {
        (Some(expected_items), Some(actual_items)) if separator == JSON_ARRAY => {
            arrays_of_one_domains_list(&expected_items, &actual_items)
        }
        (Some(expected_items), Some(actual_items)) => same_multiset(expected_items, actual_items),
        _ => false,
    }
}

/// `list_parts` without a copy of the items.
fn borrowed_parts<'a>(expected: &'a str, actual: &'a str) -> (Vec<&'a str>, Vec<&'a str>) {
    if [expected, actual]
        .iter()
        .any(|text| text.contains(',') && text.contains('&'))
    {
        return (Vec::new(), Vec::new());
    }
    let separator = if expected.contains('&') { '&' } else { ',' };
    (
        expected.split(separator).collect(),
        actual.split(separator).collect(),
    )
}

/// The items of two list-valued field texts of one format, split at the format's list
/// separator: `&` where the reference's text carries one, as a VCF value does (VCF.pm:391),
/// `,` otherwise; no items for a text carrying both, which neither separator splits into
/// the items its engine wrote.
pub fn list_parts(expected: &str, actual: &str) -> (Vec<String>, Vec<String>) {
    let (expected_items, actual_items) = borrowed_parts(expected, actual);
    let owned = |items: Vec<&str>| items.into_iter().map(str::to_string).collect();
    (owned(expected_items), owned(actual_items))
}

/// Whether two list-valued field texts carry the same items the same number of times in
/// different orders, so two or more of them (`nearest_exon_jb_item_order`).
pub fn reordered_items(expected: &str, actual: &str) -> bool {
    let (expected_items, actual_items) = borrowed_parts(expected, actual);
    expected_items != actual_items && same_multiset(expected_items, actual_items)
}

/// Whether two `5UTR_annotation` texts carry the same groups of the same `key=value`
/// pairs with the groups, or the pairs inside one, in different orders, the candidate's pairs
/// in byte order of their keys inside every group (`utrannotator_annotation_pair_order`).
pub fn reordered_pair_groups(expected: &str, actual: &str) -> bool {
    fn groups(items: Vec<&str>) -> Vec<Vec<&str>> {
        items
            .into_iter()
            .map(|group| group.split(':').collect())
            .collect()
    }
    fn sorted(mut groups: Vec<Vec<&str>>) -> Vec<Vec<&str>> {
        groups.iter_mut().for_each(|group| group.sort_unstable());
        groups
    }
    let (expected_items, actual_items) = borrowed_parts(expected, actual);
    let (expected_groups, actual_groups) = (groups(expected_items), groups(actual_items));
    if expected_groups == actual_groups {
        return false;
    }
    if !expected_groups
        .iter()
        .chain(&actual_groups)
        .flatten()
        .all(|pair| pair.contains('='))
    {
        return false;
    }
    fn key<'a>(pair: &&'a str) -> &'a str {
        pair.split('=').next().unwrap_or("")
    }
    if actual_groups
        .iter()
        .any(|group| !group.is_sorted_by(|pair, next| key(pair) <= key(next)))
    {
        return false;
    }
    same_multiset(sorted(expected_groups), sorted(actual_groups))
}

/// The fields the `order_only` class compares as lists in `fmt`, each with its
/// list separator: `,` in the default and tab formats (`convert_arrayref`,
/// Utils.pm:185-198), `&` in VCF (VCF.pm:391), and in JSON the arrays of
/// `ORDER_ONLY_JSON_FIELDS`; `CLIN_SIG` is split at `clin_sig_list_separator` alone, into
/// its term lists.
pub fn order_only_fields(fmt: &str) -> Fields {
    if fmt == "json" {
        return ORDER_ONLY_JSON_FIELDS
            .iter()
            .map(|(field, separator)| (field.to_string(), separator.to_string()))
            .collect();
    }
    let separator = if fmt == "vcf" { "&" } else { "," };
    let mut fields: Fields = ORDER_ONLY_COLUMNS
        .iter()
        .map(|column| {
            let items = if *column == "CLIN_SIG" {
                clin_sig_list_separator(fmt)
            } else {
                separator
            };
            (column.to_string(), items.to_string())
        })
        .collect();
    fields.insert(
        ORDER_ONLY_PAIR_GROUP_COLUMN.to_string(),
        format!("{PAIR_GROUPS}{separator}"),
    );
    fields
}

/// The objects of a field whose text is one JSON object or an array of them as
/// `json_scalar` renders it, the items' compact texts joined with `,` (a plugin writing one
/// object per association, `opentargets`, reads as `{..},{..}`); empty for a scalar or any
/// other text.
pub fn field_objects(text: &str) -> Vec<Value> {
    if !text.starts_with('{') || !text.ends_with('}') {
        return Vec::new();
    }
    loaded_items(text)
        .unwrap_or_default()
        .into_iter()
        .filter(|item| item.as_dict().is_some())
        .collect()
}

/// Whether `name.lower()` is `lowered`, a text already lower-cased.
pub(crate) fn lowers_to(name: &str, lowered: &str) -> bool {
    if name.is_ascii() {
        name.eq_ignore_ascii_case(lowered) && !lowered.bytes().any(|byte| byte.is_ascii_uppercase())
    } else {
        pyfmt::lower(name) == lowered
    }
}

/// The value of `key` inside the entry's `parent` field when that field is a JSON
/// object, or an array of objects, carrying the key in any case (a plugin that nests its fields
/// under its own name, `alphamissense` holding `am_class`; the JSON factory lower-cases the
/// object's key and leaves the members as the plugin spelt them): the first object's value,
/// else `None`.
pub fn json_object_member(e: &Entry, parent: &str, key: &str) -> Option<String> {
    let wanted = pyfmt::lower(key);
    let text = e.fields.get(parent).map_or("", String::as_str);
    for object in field_objects(text) {
        for (name, value) in object.as_dict().into_iter().flatten() {
            if lowers_to(name, &wanted) {
                return Some(json_scalar(value));
            }
        }
    }
    None
}

/// The consequence entries (never a record's own entry, which repeats `strand` and
/// the coordinates) carrying any of `fields` with a value other than `-`: the JSON
/// factory drops a key whose default-format value is `-` and keeps the allele keys
/// (JSON.pm:327), so a present key is a populated column, an empty `domains` array
/// included (JSON.pm:363-371 keep the key when no item re-splits). An item `parent.key`
/// counts an entry whose `parent` object carries `key` (`json_object_member`).
pub fn populated_under(entries: &[Entry], fields: &[&str]) -> i64 {
    let populated = |e: &Entry, field: &str| -> bool {
        if let Some(value) = e.fields.get(field) {
            return value != "-";
        }
        let Some((parent, key)) = field.split_once('.') else {
            return false;
        };
        json_object_member(e, parent, key).is_some_and(|value| value != "-")
    };
    let count = entries
        .iter()
        .filter(|e| !e.key.ends_with("|record") && fields.iter().any(|field| populated(e, field)))
        .count();
    i64::try_from(count).unwrap_or(i64::MAX)
}

/// The next line of `reader` without its end, a line ending at `\n`, `\r\n` or a lone
/// `\r`; `None` at the end of the text.
fn next_line(reader: &mut dyn BufRead) -> std::io::Result<Option<Vec<u8>>> {
    let mut line = Vec::new();
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok((!line.is_empty()).then_some(line));
        }
        match memchr::memchr2(b'\n', b'\r', buffer) {
            None => {
                let taken = buffer.len();
                line.extend_from_slice(buffer);
                reader.consume(taken);
            }
            Some(at) => {
                let carriage_return = buffer[at] == b'\r';
                line.extend_from_slice(&buffer[..at]);
                reader.consume(at + 1);
                if carriage_return && reader.fill_buf()?.first() == Some(&b'\n') {
                    reader.consume(1);
                }
                return Ok(Some(line));
            }
        }
    }
}

/// The columns a tab or VCF output's header declares, read from the file's leading `#`
/// lines alone: the tab `#` header line's columns, the VCF CSQ `Format:` list
/// (`csq_format_fields`) of the last `##INFO=<ID=CSQ` line, the one `parse_vcf` reads
/// the entries by (`--keep_csq` writes the kept list first); empty for a file without one.
/// Every row of these formats carries every declared column, so a column no header declares
/// appears in no row.
pub fn declared_columns(fmt: &str, path: &Path) -> Result<BTreeSet<String>> {
    let cannot_read =
        |what: String| Error::comparator(format!("cannot read {}: {what}", path_text(path)));
    let gzipped = path.as_os_str().as_encoded_bytes().ends_with(b".gz");
    let mut reader =
        BufReader::new(byte_stream(path, gzipped).map_err(|error| cannot_read(strerror(&error)))?);
    let mut declared = BTreeSet::new();
    while let Some(line) = next_line(&mut reader).map_err(|error| cannot_read(strerror(&error)))? {
        if !line.starts_with(b"#") {
            break;
        }
        let line = pyfmt::decode_utf8(&line).map_err(cannot_read)?;
        if fmt == "vcf" && line.starts_with("##INFO=<ID=CSQ") {
            declared = csq_format_fields(line)?.into_iter().collect();
        } else if fmt == "tab" && !line.starts_with("##") {
            return Ok(line[1..].split('\t').map(str::to_string).collect());
        }
    }
    Ok(declared)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn list_parts_split_at_the_separator_the_reference_carries() {
        let owned = |items: &[&str]| {
            items
                .iter()
                .map(|item| item.to_string())
                .collect::<Vec<_>>()
        };
        assert_eq!(
            list_parts("a,b", "b,a"),
            (owned(&["a", "b"]), owned(&["b", "a"]))
        );
        assert_eq!(
            list_parts("a&b", "b,a"),
            (owned(&["a", "b"]), owned(&["b,a"]))
        );
        assert_eq!(list_parts("a,b&c", "x"), (Vec::new(), Vec::new()));
        assert_eq!(list_parts("x", "a,b&c"), (Vec::new(), Vec::new()));
        assert_eq!(list_parts("", ""), (owned(&[""]), owned(&[""])));
    }

    #[test]
    fn reordered_items_are_the_same_items_in_another_order() {
        assert!(reordered_items("a,b", "b,a"));
        assert!(reordered_items("a&b&a", "a&a&b"));
        assert!(!reordered_items("a,b", "a,b"));
        assert!(!reordered_items("a,b", "b,b"));
        assert!(!reordered_items("a,b&c", "b&c,a"));
    }

    #[test]
    fn reordered_pair_groups_hold_the_candidate_to_the_byte_order_of_its_keys() {
        assert!(reordered_pair_groups("b=2:a=1", "a=1:b=2"));
        assert!(reordered_pair_groups("b=2:a=1,c=3", "c=3,a=1:b=2"));
        assert!(!reordered_pair_groups("a=1:b=2", "b=2:a=1"));
        assert!(!reordered_pair_groups("a=1:b=2", "a=1:b=2"));
        assert!(!reordered_pair_groups("b=2:a=1", "a=1:b=3"));
        assert!(!reordered_pair_groups("b:a=1", "a=1:b"));
        assert!(!reordered_pair_groups("b=2:a=1,x&y", "a=1:b=2,x&y"));
    }

    #[test]
    fn the_order_only_fields_carry_the_separator_of_their_format() {
        let fields =
            |fmt: &str| -> Vec<(String, String)> { order_only_fields(fmt).into_iter().collect() };
        let pairs = |items: &[(&str, &str)]| -> Vec<(String, String)> {
            items
                .iter()
                .map(|(name, separator)| (name.to_string(), separator.to_string()))
                .collect()
        };
        assert_eq!(
            fields("json"),
            pairs(&[("domains", "json"), ("zyg", ","), ("nearestexonjb", ",")])
        );
        assert_eq!(
            fields("vcf"),
            pairs(&[
                ("DOMAINS", "&"),
                ("FREQS", "&"),
                ("ZYG", "&"),
                ("NearestExonJB", "&"),
                ("CLIN_SIG", "%3B"),
                ("5UTR_annotation", "pairs:&")
            ])
        );
        assert_eq!(
            fields("default"),
            pairs(&[
                ("DOMAINS", ","),
                ("FREQS", ","),
                ("ZYG", ","),
                ("NearestExonJB", ","),
                ("CLIN_SIG", "%3B"),
                ("5UTR_annotation", "pairs:,")
            ])
        );
        // The tab format, and a Parquet output read as tab rows, write the `;` as it is.
        for fmt in ["tab", "parquet"] {
            let of_format = fields(fmt);
            assert_eq!(
                of_format[4],
                ("CLIN_SIG".to_string(), ";".to_string()),
                "{fmt}"
            );
            assert_eq!(of_format[..4], fields("default")[..4], "{fmt}");
            assert_eq!(of_format[5], fields("default")[5], "{fmt}");
        }
    }

    #[test]
    fn text_lists_agree_as_multisets_and_an_empty_value_with_itself_alone() {
        assert!(same_items("a,b,a", "a,a,b", ","));
        assert!(!same_items("a,b", "a,a,b", ","));
        assert!(same_items("-", "-", ","));
        assert!(!same_items("-", "", ","));
        assert!(!same_items("a", ABSENT, ","));
        assert!(same_items("b=2:a=1,c=3", "c=3,a=1:b=2", "pairs:,"));
        assert!(!same_items("b=2:a=1,c=3", "c=3,a=1:b=3", "pairs:,"));
        assert_eq!(
            list_items("b=2:a=1&c=3", "pairs:&"),
            Some(vec!["a=1:b=2".to_string(), "c=3".to_string()])
        );
    }

    /// The objects of a made-up `domains` array, `a` to `c`, as `list_items` gives them.
    fn objects(names: &str) -> Vec<String> {
        names
            .chars()
            .map(|name| format!("{{\"db\":\"D\",\"name\":\"{name}\"}}"))
            .collect()
    }

    #[test]
    fn two_domains_arrays_agree_where_one_list_in_two_orders_gives_both() {
        let one_list = |expected: &str, actual: &str| {
            arrays_of_one_domains_list(&objects(expected), &objects(actual))
        };
        // The same objects the same number of times, in any order.
        assert!(one_list("abc", "cab"));
        assert!(one_list("aab", "baa"));
        // An unmatched entry repeats the entry before it, or nothing where it stands first.
        assert!(one_list("abc", "baac"));
        assert!(one_list("aab", "bba"));
        assert!(one_list("aaab", "ba"));
        assert!(one_list("abb", "bbba"));
        // Two runs of one object are two matched entries, which the other side holds too.
        assert!(one_list("aba", "aab"));
        assert!(!one_list("aba", "ab"));
        assert!(!one_list("ab", "babab"));
        // An object one side alone holds is an entry the other list lacks.
        assert!(!one_list("abc", "abb"));
        assert!(!one_list("aa", "ab"));
        assert!(!one_list("a", ""));
    }

    #[test]
    fn the_json_array_marker_reads_objects_and_agrees_as_one_list() {
        let (pfam, panther) = (
            r#"{"db": "Pfam", "name": "PF99001"}"#,
            r#"{"name": "PTHR99001", "db": "PANTHER"}"#,
        );
        let array = |objects: &[&str]| objects.join(",");
        assert_eq!(
            list_items(&array(&[panther, pfam]), JSON_ARRAY),
            Some(vec![
                r#"{"db":"PANTHER","name":"PTHR99001"}"#.to_string(),
                r#"{"db":"Pfam","name":"PF99001"}"#.to_string()
            ])
        );
        let (once, twice) = (array(&[pfam, panther]), array(&[panther, panther, pfam]));
        assert!(same_items(&once, &twice, JSON_ARRAY));
        assert!(same_items(&once, &array(&[panther, pfam]), JSON_ARRAY));
        assert!(!same_items(
            &once,
            &array(&[pfam, panther, pfam]),
            JSON_ARRAY
        ));
        assert!(!same_items(&once, "", JSON_ARRAY));
        assert!(!same_items(&once, "Pfam:PF99001", JSON_ARRAY));
        assert!(same_items("", "", JSON_ARRAY));
        // A text list still agrees as a multiset alone.
        assert!(!same_items("a,b", "b,b,a", ","));
    }

    #[test]
    fn a_line_ends_at_a_line_feed_a_carriage_return_or_both() {
        let mut reader = BufReader::with_capacity(4, &b"#a\r\n#bcd\r#e\n\n#f"[..]);
        let mut lines = Vec::new();
        while let Some(line) = next_line(&mut reader).unwrap() {
            lines.push(String::from_utf8(line).unwrap());
        }
        assert_eq!(lines, ["#a", "#bcd", "#e", "", "#f"]);
    }

    #[test]
    fn a_tab_header_declares_the_columns_of_its_first_single_hash_line() {
        use crate::cache::testing::{plant, scratch};
        let root = scratch("declared-columns");
        let declared = |name: &str, text: &str| -> Vec<String> {
            plant(&root, name, text);
            declared_columns("tab", &root.join(name))
                .unwrap()
                .into_iter()
                .collect()
        };
        assert_eq!(
            declared(
                "plain.txt",
                "## a\n## b\n#Uploaded_variation\tLocation\tSYMBOL\nrs1\t21:5\tG\n"
            ),
            ["Location", "SYMBOL", "Uploaded_variation"]
        );
        assert_eq!(
            declared(
                "returns.txt",
                "## a\r\n#Uploaded_variation\tLocation\r\nrs1\t21:5\r\n"
            ),
            ["Location", "Uploaded_variation"]
        );
        assert_eq!(declared("bare.txt", "#\n"), [""]);
        assert!(declared("data-first.txt", "rs1\t21:5\n#Uploaded_variation\n").is_empty());
        assert!(declared("hashes.txt", "## a\n## b\n").is_empty());
        assert!(declared("empty.txt", "").is_empty());
        assert!(declared("empty.gz", "").is_empty());
        let refused = declared_columns("tab", &root.join("absent.txt"))
            .unwrap_err()
            .to_string();
        assert_eq!(
            refused,
            format!(
                "cannot read {}: No such file or directory",
                root.join("absent.txt").display()
            )
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_damaged_gzip_file_declares_the_columns_read_before_its_damage() {
        use crate::cache::testing::scratch;
        use std::io::Write;
        let root = scratch("declared-damaged");
        let text = "## a\n#Uploaded_variation\tLocation\tExtra\nrs1\t21:5\t-\n";
        let mut zipped = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        zipped.write_all(text.as_bytes()).unwrap();
        let whole = zipped.finish().unwrap();
        let mut checksum = whole.clone();
        let at = checksum.len() - 8;
        checksum[at] ^= 0xff;
        let after = [&whole[..], b"no"].concat();
        let declared = |name: &str, bytes: &[u8], fmt: &str| {
            let path = root.join(name);
            std::fs::write(&path, bytes).unwrap();
            declared_columns(fmt, &path)
                .map(|columns| columns.into_iter().collect::<Vec<String>>())
                .map_err(|error| error.to_string().replace(&path_text(&path), "<path>"))
        };
        // The text before a cut, a wrong checksum or bytes that are no member is read first,
        // and a header ends before any of them.
        for (name, bytes) in [
            ("cut.txt.gz", &whole[..whole.len() - 3]),
            ("checksum.txt.gz", &checksum[..]),
            ("after.txt.gz", &after[..]),
        ] {
            assert_eq!(
                declared(name, bytes, "tab").unwrap(),
                ["Extra", "Location", "Uploaded_variation"],
                "{name}"
            );
            assert!(declared(name, bytes, "vcf").unwrap().is_empty(), "{name}");
        }
        for fmt in ["tab", "vcf"] {
            assert_eq!(
                declared("plain.txt.gz", text.as_bytes(), fmt).unwrap_err(),
                "cannot read <path>: Not a gzipped file (b'##')"
            );
            assert_eq!(
                declared("head.txt.gz", &whole[..10], fmt).unwrap_err(),
                "cannot read <path>: Compressed file ended before the end-of-stream marker was reached"
            );
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_vcf_header_declares_the_format_list_of_its_last_csq_line() {
        use crate::cache::testing::{plant, scratch};
        use std::io::Write;
        let root = scratch("declared-csq");
        let header = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Consequence|IMPACT\">\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Kept. Format: Allele|SYMBOL|Gene\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t5\t.\tA\tT\t.\t.\tCSQ=T|x|y\n";
        plant(&root, "plain.vcf", header);
        let mut zipped = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        zipped.write_all(header.as_bytes()).unwrap();
        std::fs::write(root.join("zipped.vcf.gz"), zipped.finish().unwrap()).unwrap();
        for name in ["plain.vcf", "zipped.vcf.gz"] {
            let declared: Vec<String> = declared_columns("vcf", &root.join(name))
                .unwrap()
                .into_iter()
                .collect();
            assert_eq!(declared, ["Allele", "Gene", "SYMBOL"], "{name}");
            assert!(declared_columns("tab", &root.join(name)).unwrap().len() == 8);
        }
        plant(
            &root,
            "no-list.vcf",
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"none\">\n#CHROM\n",
        );
        assert_eq!(
            declared_columns("vcf", &root.join("no-list.vcf"))
                .unwrap_err()
                .to_string(),
            "CSQ header without a Format list: ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"none\">"
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn the_objects_of_a_field_are_those_of_one_object_or_of_several_joined() {
        let compact = |text: &str| -> Vec<String> {
            field_objects(text)
                .iter()
                .map(pyjson::dumps_compact)
                .collect()
        };
        assert_eq!(
            compact(r#"{"am_class":"likely_benign","am_pathogenicity":0.1261}"#),
            [r#"{"am_class":"likely_benign","am_pathogenicity":0.1261}"#]
        );
        assert_eq!(
            compact(r#"{"db":"Pfam","name":"PF1"},{"db":"Pfam","name":"PF2"}"#),
            [
                r#"{"db":"Pfam","name":"PF1"}"#,
                r#"{"db":"Pfam","name":"PF2"}"#
            ]
        );
        assert!(compact("[1,2]").is_empty());
        assert!(compact("MODERATE").is_empty());
        assert!(compact("{").is_empty());
        assert!(compact("{x},{y}").is_empty());
        assert_eq!(
            compact(r#"{"a":1},2,{"b":2}"#),
            [r#"{"a":1}"#, r#"{"b":2}"#]
        );
    }
}
