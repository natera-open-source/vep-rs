// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The index page, `docs/intended-divergences.md`: its sections, the rows of its entry
//! tables and the rule that holds each entry's input class to its tier and its section.
//! The page is read as `scripts/divergences/check_divergence_order.py` reads it, which
//! checks the same page against the trigger-frequency table.

use crate::pyfmt::{py_repr_str, py_splitlines, py_strip};
use crate::pyre::Regex;

pub const INPUT_CLASSES: [&str; 3] = ["observed", "unobserved", "malformed"];
/// The heading of the section that holds the `malformed` entries, without any `{#id}`.
pub const MALFORMED_SECTION: &str = "Differences only malformed input reaches";
pub const ENTRY_TABLE_CELLS: [&str; 5] = [
    "Entry",
    "In one sentence",
    "Tier and count",
    "Verified against",
    "Page",
];
pub const COUNT_CELL: usize = 2;
pub const VERIFIED_CELL: usize = 3;

/// A rendered `Tier and count` cell: the tier, the figure and what reaches it.
pub const CELL: &str = r"^(?P<tier>[1-3]) · (?P<figure>[^·]+?) · (?P<tail>.+)$";

#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub key: String,
    pub cells: Vec<String>,
    pub line: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Section {
    pub title: String,
    pub line: usize,
    pub rows: Vec<Row>,
    pub subsections: Vec<Section>,
    /// The keys in code font in the opening prose of the malformed-input section, each with
    /// its line and whether it is linked; empty for every other section.
    pub named: Vec<(String, usize, bool)>,
}

/// The patterns the page is read with, compiled once.
pub struct PagePatterns {
    heading: Regex,
    heading_id: Regex,
    entry_cell: Regex,
    named_key: Regex,
    closing_pipe: Regex,
    pipe: Regex,
    pub cell: Regex,
}

impl PagePatterns {
    pub fn new() -> Self {
        let compile = |pattern: &str| Regex::new(pattern).expect("a pattern of the page reader");
        PagePatterns {
            heading: compile(r"^(?P<marks>#{2,4}) (?P<title>.+)$"),
            heading_id: compile(r"\s*\{#[^}]*\}$"),
            entry_cell: compile(r"^\[`(?P<key>[a-z0-9_]+)`\]\([^)]+\)$"),
            named_key: compile(r"(?P<link>\[)?`(?P<key>[a-z0-9_]+)`(?(link)\]\([^)]+\))"),
            closing_pipe: compile(r"(?<!\\)\|$"),
            pipe: compile(r"(?<!\\)\|"),
            cell: compile(CELL),
        }
    }

    /// The cells of a Markdown table row, a `\|` inside a cell kept as a pipe.
    fn split_cells(&self, line: &str) -> Vec<String> {
        let stripped = py_strip(line);
        let inner = self
            .closing_pipe
            .sub("", stripped.strip_prefix('|').unwrap_or(stripped));
        self.pipe
            .split(&inner)
            .iter()
            .map(|cell| py_strip(cell).replace("\\|", "|"))
            .collect()
    }

    /// A heading's words, without a trailing `{#id}` attribute.
    fn heading_text(&self, title: &str) -> String {
        self.heading_id.sub("", title)
    }

    pub fn is_malformed_section(&self, section: &Section) -> bool {
        self.heading_text(&section.title) == MALFORMED_SECTION
    }
}

impl Default for PagePatterns {
    fn default() -> Self {
        Self::new()
    }
}

/// Where the open section stands in the tree.
#[derive(Clone, Copy, PartialEq)]
enum Open {
    None,
    Top,
    Category,
    Group,
}

/// The `##` sections of the page with their entry rows and `###` subsections and, inside
/// the malformed-input section, the keys its opening prose writes in code font and the
/// `####` subsections of each `###` category. An entry table opens on its header row and
/// runs to the next blank line or section heading. A `####` line outside the
/// malformed-input section is body text. The error is one line that names the page's line.
pub fn parse_page(text: &str, patterns: &PagePatterns) -> Result<Vec<Section>, String> {
    let mut sections: Vec<Section> = Vec::new();
    let mut open = Open::None;
    // Whether a `###` section is open for a `####` heading to go under.
    let mut category = false;
    let mut in_table = false;
    for (index, line) in py_splitlines(text).into_iter().enumerate() {
        let n = index + 1;
        let heading = patterns.heading.match_start(line);
        let mut level = heading
            .as_ref()
            .map_or(0, |m| m.name("marks").map_or(0, |marks| marks.len()));
        if level == 4
            && !sections
                .last()
                .is_some_and(|last| patterns.is_malformed_section(last))
        {
            level = 0;
        }
        if level != 0 {
            let title = heading
                .and_then(|m| m.name("title"))
                .expect("a heading's title");
            let section = Section {
                title: py_strip(&title).to_string(),
                line: n,
                rows: Vec::new(),
                subsections: Vec::new(),
                named: Vec::new(),
            };
            match level {
                2 => {
                    sections.push(section);
                    category = false;
                    open = Open::Top;
                }
                3 => {
                    let Some(top) = sections.last_mut() else {
                        return Err(format!("line {n}: a subsection before any section"));
                    };
                    top.subsections.push(section);
                    category = true;
                    open = Open::Category;
                }
                _ => {
                    if !category {
                        return Err(format!(
                            "line {n}: a `####` heading outside a `###` subsection"
                        ));
                    }
                    sections
                        .last_mut()
                        .and_then(|top| top.subsections.last_mut())
                        .expect("an open category")
                        .subsections
                        .push(section);
                    open = Open::Group;
                }
            }
            in_table = false;
            continue;
        }
        let unindented = line.trim_start_matches(' ');
        if !in_table
            && unindented.starts_with('|')
            && line.len() - unindented.len() <= 3
            && patterns.split_cells(line) == ENTRY_TABLE_CELLS
        {
            if open == Open::None {
                return Err(format!("line {n}: an entry table before any heading"));
            }
            in_table = true;
        } else if in_table && !py_strip(line).is_empty() {
            // A table runs to the next blank line; every line until then is a row, however
            // it is indented and with or without its outer pipes.
            let cells = patterns.split_cells(line);
            let is_rule = cells
                .iter()
                .all(|cell| cell.contains('-') && cell.chars().all(|c| c == '-' || c == ':'));
            if is_rule {
                continue;
            }
            if cells.len() != ENTRY_TABLE_CELLS.len() {
                return Err(format!("line {n}: an entry row with {} cells", cells.len()));
            }
            let Some(key) = patterns
                .entry_cell
                .match_start(&cells[0])
                .and_then(|m| m.name("key"))
            else {
                return Err(format!(
                    "line {n}: the first cell is not a linked key in code font: {}",
                    py_repr_str(&cells[0])
                ));
            };
            let current = open_section(&mut sections, open);
            current.rows.push(Row {
                key,
                cells,
                line: n,
            });
        } else {
            in_table = false;
            // The `##` section is open until its first `###` heading, so a subsection that
            // repeats the title names nothing.
            if open == Open::Top {
                let top = sections.last_mut().expect("an open section");
                if patterns.is_malformed_section(top) {
                    for found in patterns.named_key.find_iter(line) {
                        top.named.push((
                            found.name("key").expect("a key"),
                            n,
                            found.name("link").is_some(),
                        ));
                    }
                }
            }
        }
    }
    Ok(sections)
}

fn open_section(sections: &mut [Section], open: Open) -> &mut Section {
    let top = sections.last_mut().expect("an open section");
    match open {
        Open::None | Open::Top => top,
        Open::Category => top.subsections.last_mut().expect("an open category"),
        Open::Group => top
            .subsections
            .last_mut()
            .and_then(|category| category.subsections.last_mut())
            .expect("an open group"),
    }
}

/// A section and every subsection below it, in page order.
pub fn walk(section: &Section) -> Vec<&Section> {
    let mut out = vec![section];
    for sub in &section.subsections {
        out.extend(walk(sub));
    }
    out
}

/// Every section that lists entries, in page order.
pub fn leaves(sections: &[Section]) -> Vec<&Section> {
    sections
        .iter()
        .flat_map(walk)
        .filter(|section| !section.rows.is_empty())
        .collect()
}

/// How many entries name each input class, as a summary line states it.
pub fn class_counts<'a>(classes: impl Iterator<Item = &'a str> + Clone) -> String {
    INPUT_CLASSES
        .iter()
        .map(|input_class| {
            let count = classes.clone().filter(|named| named == input_class).count();
            format!("{count} {input_class}")
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// Every violation of the input classes against the page's sections and tiers, as one line
/// each: a tier 1 or tier 2 entry is `observed`, a `malformed` entry stands in the
/// malformed-input section and no other entry does. A row whose key has no class is not
/// judged, and one whose cell states no tier is held to its section alone.
pub fn check_classes(
    sections: &[Section],
    classes: &[(String, String)],
    patterns: &PagePatterns,
) -> Vec<String> {
    let mut errors = Vec::new();
    for top in sections {
        let inside = patterns.is_malformed_section(top);
        for section in walk(top) {
            for row in &section.rows {
                let Some(input_class) = classes
                    .iter()
                    .find(|(key, _)| *key == row.key)
                    .map(|(_, class)| class.as_str())
                else {
                    continue;
                };
                let tier = patterns
                    .cell
                    .match_start(&row.cells[COUNT_CELL])
                    .and_then(|m| m.name("tier"));
                if let Some(tier) = tier.filter(|tier| tier != "3") {
                    if input_class != "observed" {
                        errors.push(format!(
                            "line {}: {} is tier {tier} with input_class {}: a tier 1 or tier 2 entry is observed",
                            row.line,
                            row.key,
                            py_repr_str(input_class)
                        ));
                    }
                }
                if input_class == "malformed" && !inside {
                    errors.push(format!(
                        "line {}: {} is malformed and stands under {}, outside {}",
                        row.line,
                        row.key,
                        py_repr_str(&section.title),
                        py_repr_str(MALFORMED_SECTION)
                    ));
                }
                if input_class != "malformed" && inside {
                    errors.push(format!(
                        "line {}: {} stands under {} with input_class {}",
                        row.line,
                        row.key,
                        py_repr_str(MALFORMED_SECTION),
                        py_repr_str(input_class)
                    ));
                }
            }
        }
    }
    errors
}

#[cfg(test)]
#[path = "page_tests.rs"]
mod tests;
