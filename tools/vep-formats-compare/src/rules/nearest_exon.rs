// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The NearestExonJB items.

use crate::rules::whole_text;
use std::sync::LazyLock;

/// One item of a `NearestExonJB` value, the exon id, distance, side (`start_end` on a tie)
/// and length: `|` joins the parts and `,` the items in the default, tab and JSON formats,
/// `+` and `&` in the VCF format.
pub const NEAREST_EXON_ITEM: &str = "([^|+&,]*)[|+](\\d+)[|+](start|end|start_end)[|+](\\d+)";

pub const NEAREST_EXON_ITEM_RE_PATTERN: &str =
    r"([^|+&,]*)[|+](\d+)[|+](start|end|start_end)[|+](\d+)";
pub static NEAREST_EXON_ITEM_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(NEAREST_EXON_ITEM_RE_PATTERN).expect("the pattern compiles")
});

pub const NEAREST_EXON_LIST_RE_PATTERN: &str = r"([^|+&,]*)[|+](\d+)[|+](start|end|start_end)[|+](\d+)(?:[,&]([^|+&,]*)[|+](\d+)[|+](start|end|start_end)[|+](\d+))*";
pub static NEAREST_EXON_LIST_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(NEAREST_EXON_LIST_RE_PATTERN).expect("the pattern compiles")
});

/// `NEAREST_EXON_LIST_RE` held to the whole text, as `fullmatch` holds it.
static NEAREST_EXON_LIST_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NEAREST_EXON_LIST_RE_PATTERN));

/// The (exon id, distance, side, length) items of a `NearestExonJB` value in any
/// format's spelling; `None` for a value that is not a list of such items.
pub fn nearest_exon_items(value: &str) -> Option<Vec<(String, String, String, String)>> {
    if !NEAREST_EXON_LIST_WHOLE_RE.is_match(value) {
        return None;
    }
    Some(
        NEAREST_EXON_ITEM_RE
            .captures_iter(value)
            .map(|found| {
                let part = |group: usize| found[group].to_string();
                (part(1), part(2), part(3), part(4))
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn items(value: &str) -> Option<Vec<String>> {
        nearest_exon_items(value).map(|items| {
            items
                .into_iter()
                .map(|(exon, distance, side, length)| format!("{exon} {distance} {side} {length}"))
                .collect()
        })
    }

    #[test]
    fn a_value_is_read_as_its_items_in_any_format() {
        assert_eq!(
            items("ENSE1|3|start|117,ENSE2|3|end|207"),
            Some(vec![
                "ENSE1 3 start 117".to_string(),
                "ENSE2 3 end 207".to_string()
            ])
        );
        assert_eq!(
            items("ENSE1+3+start_end+117&+0+end+9"),
            Some(vec![
                "ENSE1 3 start_end 117".to_string(),
                " 0 end 9".to_string()
            ])
        );
        // An exon without an id is an item with an empty id.
        assert_eq!(
            items("|12|start|40"),
            Some(vec![" 12 start 40".to_string()])
        );
    }

    #[test]
    fn a_text_that_is_no_list_of_items_has_none() {
        for value in [
            "",
            "-",
            "ENSE1|3|start",
            "ENSE1|3|middle|117",
            "ENSE1|x|start|117",
            "ENSE1|3|start|117,",
            "ENSE1|3|start|117;ENSE2|3|end|207",
            "ENSE1|3|start|117\n",
        ] {
            assert_eq!(items(value), None, "{value:?}");
        }
    }
}
