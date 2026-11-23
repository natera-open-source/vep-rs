// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CPython rules the tool's bytes follow, character for character: the width of a
//! `str` under a format specification, `f"{x:5.1f}"`, `str(pathlib.Path(arg).parent / name)`
//! and the `\d` of `re`.

/// The code points of UTF-8 text: every byte that does not continue a sequence.
pub fn char_count(text: &[u8]) -> usize {
    text.iter().filter(|byte| (**byte & 0xc0) != 0x80).count()
}

/// Appends `f"{text:{width}}"` for a `str`: the text whole, then blanks up to `width` code
/// points.
pub fn push_padded(out: &mut Vec<u8>, text: &[u8], width: usize) {
    out.extend_from_slice(text);
    let blanks = width.saturating_sub(char_count(text));
    out.resize(out.len() + blanks, b' ');
}

/// `f"{100 * count / total:5.1f}"`: the quotient of two integers as the nearest double, then
/// one decimal of its exact value, ties to the even digit, right-aligned in five characters.
pub fn percent(count: u64, total: u64) -> String {
    format!("{:5.1}", (100 * count) as f64 / total as f64)
}

/// `str(pathlib.PurePosixPath(input).parent / name)` on the argument's bytes: repeated
/// separators and `.` components dropped, exactly two leading slashes kept, `..` left
/// alone, the last component replaced by `name`, which stands alone where there is none.
pub fn sibling(input: &[u8], name: &[u8]) -> Vec<u8> {
    let leading_slashes = input.iter().take_while(|byte| **byte == b'/').count();
    let root: &[u8] = match leading_slashes {
        0 => b"",
        2 => b"//",
        _ => b"/",
    };
    let mut parts: Vec<&[u8]> = input
        .split(|byte| *byte == b'/')
        .filter(|part| !part.is_empty() && *part != b".")
        .collect();
    parts.pop();
    parts.push(name);
    [root, &parts.join(&b'/')].concat()
}

/// The first code point of each run of decimal digits (general category Nd) and the length
/// of the run, in Unicode 16.0.
const DIGIT_RUNS: [(u32, u32); 71] = [
    (0x0030, 10),
    (0x0660, 10),
    (0x06f0, 10),
    (0x07c0, 10),
    (0x0966, 10),
    (0x09e6, 10),
    (0x0a66, 10),
    (0x0ae6, 10),
    (0x0b66, 10),
    (0x0be6, 10),
    (0x0c66, 10),
    (0x0ce6, 10),
    (0x0d66, 10),
    (0x0de6, 10),
    (0x0e50, 10),
    (0x0ed0, 10),
    (0x0f20, 10),
    (0x1040, 10),
    (0x1090, 10),
    (0x17e0, 10),
    (0x1810, 10),
    (0x1946, 10),
    (0x19d0, 10),
    (0x1a80, 10),
    (0x1a90, 10),
    (0x1b50, 10),
    (0x1bb0, 10),
    (0x1c40, 10),
    (0x1c50, 10),
    (0xa620, 10),
    (0xa8d0, 10),
    (0xa900, 10),
    (0xa9d0, 10),
    (0xa9f0, 10),
    (0xaa50, 10),
    (0xabf0, 10),
    (0xff10, 10),
    (0x104a0, 10),
    (0x10d30, 10),
    (0x10d40, 10),
    (0x11066, 10),
    (0x110f0, 10),
    (0x11136, 10),
    (0x111d0, 10),
    (0x112f0, 10),
    (0x11450, 10),
    (0x114d0, 10),
    (0x11650, 10),
    (0x116c0, 10),
    (0x116d0, 20),
    (0x11730, 10),
    (0x118e0, 10),
    (0x11950, 10),
    (0x11bf0, 10),
    (0x11c50, 10),
    (0x11d50, 10),
    (0x11da0, 10),
    (0x11f50, 10),
    (0x16130, 10),
    (0x16a60, 10),
    (0x16ac0, 10),
    (0x16b50, 10),
    (0x16d70, 10),
    (0x1ccf0, 10),
    (0x1d7ce, 50),
    (0x1e140, 10),
    (0x1e2f0, 10),
    (0x1e4f0, 10),
    (0x1e5f1, 10),
    (0x1e950, 10),
    (0x1fbf0, 10),
];

/// `\d` of `re` on a `str`: a decimal digit of any script.
pub fn is_decimal_digit(c: char) -> bool {
    let code = u32::from(c);
    let run = DIGIT_RUNS.partition_point(|(first, _)| *first <= code);
    run > 0 && code - DIGIT_RUNS[run - 1].0 < DIGIT_RUNS[run - 1].1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_str_is_padded_to_its_width_in_code_points_and_never_cut() {
        let padded = |text: &str, width: usize| {
            let mut out = Vec::new();
            push_padded(&mut out, text.as_bytes(), width);
            String::from_utf8(out).unwrap()
        };
        assert_eq!(padded("true_missing", 25), "true_missing             ");
        assert_eq!(
            padded("\u{e9}\u{e9}\u{e9}", 25),
            format!("\u{e9}\u{e9}\u{e9}{}", " ".repeat(22))
        );
        assert_eq!(padded(&"x".repeat(30), 25), "x".repeat(30));
        assert_eq!(padded("", 12), " ".repeat(12));
    }

    #[test]
    fn a_percentage_is_one_decimal_of_the_exact_quotient() {
        let cases: [(u64, u64, &str); 31] = [
            (1, 400, "  0.2"),
            (3, 400, "  0.8"),
            (5, 400, "  1.2"),
            (7, 400, "  1.8"),
            (1, 800, "  0.1"),
            (3, 800, "  0.4"),
            (5, 800, "  0.6"),
            (7, 800, "  0.9"),
            (1, 8, " 12.5"),
            (1, 16, "  6.2"),
            (3, 16, " 18.8"),
            (1, 3, " 33.3"),
            (2, 3, " 66.7"),
            (1, 7, " 14.3"),
            (99999, 100000, "100.0"),
            (1, 1, "100.0"),
            (1, 1000000, "  0.0"),
            (5, 5, "100.0"),
            (1, 6, " 16.7"),
            (5, 6, " 83.3"),
            (1, 32, "  3.1"),
            (19, 20, " 95.0"),
            (199, 200, " 99.5"),
            (1999, 2000, "100.0"),
            (1, 2, " 50.0"),
            (1, 2000, "  0.1"),
            (3, 2000, "  0.1"),
            (9, 2000, "  0.5"),
            (123456, 7890123, "  1.6"),
            (10934, 301350, "  3.6"),
            (290416, 301350, " 96.4"),
        ];
        for (count, total, text) in cases {
            assert_eq!(percent(count, total), text, "{count} of {total}");
        }
    }

    #[test]
    fn the_default_table_stands_beside_the_input_as_pathlib_names_it() {
        let cases = [
            ("discordant.tsv", "discordant_classified.tsv"),
            ("a/b/discordant.tsv", "a/b/discordant_classified.tsv"),
            ("/discordant.tsv", "/discordant_classified.tsv"),
            ("/", "/discordant_classified.tsv"),
            ("//", "//discordant_classified.tsv"),
            ("///x/y", "/x/discordant_classified.tsv"),
            ("//x/y", "//x/discordant_classified.tsv"),
            ("//x", "//discordant_classified.tsv"),
            (".", "discordant_classified.tsv"),
            ("..", "discordant_classified.tsv"),
            ("a/..", "a/discordant_classified.tsv"),
            ("../a", "../discordant_classified.tsv"),
            ("./a", "discordant_classified.tsv"),
            ("a/./b/", "a/discordant_classified.tsv"),
            ("a//b//", "a/discordant_classified.tsv"),
            ("", "discordant_classified.tsv"),
            ("a", "discordant_classified.tsv"),
            ("/a", "/discordant_classified.tsv"),
            ("a/b/../c", "a/b/../discordant_classified.tsv"),
            ("x/.", "discordant_classified.tsv"),
            ("./.", "discordant_classified.tsv"),
        ];
        for (input, expected) in cases {
            assert_eq!(
                String::from_utf8(sibling(input.as_bytes(), b"discordant_classified.tsv")).unwrap(),
                expected,
                "{input:?}"
            );
        }
    }

    #[test]
    fn a_decimal_digit_of_any_script_is_a_digit_and_other_numerals_are_not() {
        for digit in [
            '0',
            '9',
            '\u{663}',
            '\u{ff15}',
            '\u{116e3}',
            '\u{1d7ff}',
            '\u{1fbf9}',
        ] {
            assert!(is_decimal_digit(digit), "{digit:?}");
        }
        for other in [
            '/',
            ':',
            'a',
            '\u{b2}',
            '\u{bd}',
            '\u{2167}',
            '\u{116ca}',
            '\u{1d7cd}',
            '\u{1fbfa}',
        ] {
            assert!(!is_decimal_digit(other), "{other:?}");
        }
        let digits = (0..=0x10ffff)
            .filter_map(char::from_u32)
            .filter(|c| is_decimal_digit(*c));
        assert_eq!(digits.count(), 760);
    }
}
