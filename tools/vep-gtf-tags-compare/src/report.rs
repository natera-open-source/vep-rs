// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of the two sides, the JSON report and the stdout lines.

use crate::gtf::{CodeSets, CODES};
use crate::pyfmt::{push_json_str, py_grouped};
use std::fmt::Write as _;
use std::io::{self, Write};

/// A transcript on which the two sides disagree about one code.
pub struct Example<'a> {
    pub transcript: &'a [u8],
    pub gtf: bool,
    pub cache: bool,
}

#[derive(Default)]
pub struct CodeRow<'a> {
    pub agree: usize,
    pub disagree: usize,
    pub present_in_gtf: usize,
    pub present_in_cache: usize,
    pub examples: Vec<Example<'a>>,
}

pub struct Report<'a> {
    pub gtf: usize,
    pub cache: usize,
    pub shared: usize,
    pub per_code: [CodeRow<'a>; 4],
}

impl Report<'_> {
    pub fn exact_agreement(&self) -> bool {
        self.per_code.iter().all(|row| row.disagree == 0)
    }
}

/// Compares the code sets of the transcripts both sides hold, in the order of their ids by
/// code point. A code keeps its first disagreements as examples while it holds fewer than
/// `max_mismatches` of them.
pub fn compare<'a>(gtf: &'a CodeSets, cache: &CodeSets, max_mismatches: i64) -> Report<'a> {
    // The first sixteen bytes of an id, zero-padded, order two ids as the ids themselves do
    // wherever they differ, so most comparisons of the sort read no id.
    let leading = |id: &[u8]| {
        let mut bytes = [0u8; 16];
        let len = id.len().min(16);
        bytes[..len].copy_from_slice(&id[..len]);
        u128::from_be_bytes(bytes)
    };
    let mut shared: Vec<(u128, &[u8], u8, u8)> = gtf
        .iter()
        .filter_map(|(id, in_gtf)| {
            cache
                .get(id)
                .map(|in_cache| (leading(id), &**id, *in_gtf, *in_cache))
        })
        .collect();
    // An id both sides hold came out of the GTF and is UTF-8, which orders as its code
    // points do.
    shared.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(b.1)));
    let mut per_code: [CodeRow; 4] = Default::default();
    for &(_, id, in_gtf, in_cache) in &shared {
        for (bit, row) in per_code.iter_mut().enumerate() {
            let gtf = in_gtf >> bit & 1 == 1;
            let cache = in_cache >> bit & 1 == 1;
            row.present_in_gtf += usize::from(gtf);
            row.present_in_cache += usize::from(cache);
            if gtf == cache {
                row.agree += 1;
            } else {
                row.disagree += 1;
                let kept = i64::try_from(row.examples.len()).unwrap_or(i64::MAX);
                if kept < max_mismatches {
                    row.examples.push(Example {
                        transcript: id,
                        gtf,
                        cache,
                    });
                }
            }
        }
    }
    Report {
        gtf: gtf.len(),
        cache: cache.len(),
        shared: shared.len(),
        per_code,
    }
}

/// Writes the report as `json.dump(report, fh, indent=1)` writes it, and the newline
/// behind it.
pub fn write_json(report: &Report, out: &mut impl Write) -> io::Result<()> {
    let mut text = String::new();
    write!(
        text,
        "{{\n \"transcripts\": {{\n  \"gtf\": {},\n  \"cache\": {},\n  \"shared\": {},\n  \
         \"gtf_only\": {},\n  \"cache_only\": {}\n }},\n \"per_code\": {{\n",
        report.gtf,
        report.cache,
        report.shared,
        report.gtf - report.shared,
        report.cache - report.shared,
    )
    .unwrap();
    for (i, (code, row)) in CODES.iter().zip(&report.per_code).enumerate() {
        write!(
            text,
            "  \"{code}\": {{\n   \"agree\": {},\n   \"disagree\": {},\n   \
             \"present_in_gtf\": {},\n   \"present_in_cache\": {},\n   \"examples\": [",
            row.agree, row.disagree, row.present_in_gtf, row.present_in_cache,
        )
        .unwrap();
        for (n, example) in row.examples.iter().enumerate() {
            text.push_str(if n == 0 { "\n    {\n" } else { ",\n    {\n" });
            text.push_str("     \"transcript\": ");
            push_json_str(&mut text, &String::from_utf8_lossy(example.transcript));
            write!(
                text,
                ",\n     \"gtf\": {},\n     \"cache\": {}\n    }}",
                example.gtf, example.cache
            )
            .unwrap();
            if text.len() >= 1 << 16 {
                out.write_all(text.as_bytes())?;
                text.clear();
            }
        }
        text.push_str(if row.examples.is_empty() {
            "]\n  }"
        } else {
            "\n   ]\n  }"
        });
        text.push_str(if i + 1 < CODES.len() { ",\n" } else { "\n" });
    }
    writeln!(
        text,
        " }},\n \"exact_agreement\": {}\n}}",
        report.exact_agreement()
    )
    .unwrap();
    out.write_all(text.as_bytes())
}

/// The stdout lines: the transcript counts, then a line a code with its first five
/// examples under it.
pub fn summary(report: &Report) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(
        format!(
            "transcripts: gtf {} cache {} shared {} gtf_only {} cache_only {}\n",
            py_grouped(report.gtf),
            py_grouped(report.cache),
            py_grouped(report.shared),
            py_grouped(report.gtf - report.shared),
            py_grouped(report.cache - report.shared),
        )
        .as_bytes(),
    );
    for (code, row) in CODES.iter().zip(&report.per_code) {
        out.extend_from_slice(
            format!(
                "  {code:<16} agree {:>9}  disagree {:>7}  present gtf {:>8} cache {:>8}\n",
                py_grouped(row.agree),
                py_grouped(row.disagree),
                py_grouped(row.present_in_gtf),
                py_grouped(row.present_in_cache),
            )
            .as_bytes(),
        );
        for example in row.examples.iter().take(5) {
            out.extend_from_slice(b"    ");
            out.extend_from_slice(example.transcript);
            out.extend_from_slice(
                format!(
                    ": gtf={} cache={}\n",
                    py_bool(example.gtf),
                    py_bool(example.cache)
                )
                .as_bytes(),
            );
        }
    }
    out
}

fn py_bool(value: bool) -> &'static str {
    if value {
        "True"
    } else {
        "False"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sets(entries: &[(&str, u8)]) -> CodeSets {
        entries
            .iter()
            .map(|(id, codes)| (id.as_bytes().into(), *codes))
            .collect()
    }

    fn json(report: &Report) -> String {
        let mut out = Vec::new();
        write_json(report, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn shared_transcripts_are_compared_code_by_code_in_id_order() {
        let gtf = sets(&[("T3", 1), ("T1", 1 | 2), ("T2", 4), ("G", 8)]);
        let cache = sets(&[("T1", 1), ("T2", 4 | 8), ("T3", 0), ("C1", 1), ("C2", 2)]);
        let report = compare(&gtf, &cache, 50);
        assert_eq!((report.gtf, report.cache, report.shared), (4, 5, 3));
        assert!(!report.exact_agreement());
        let rows: Vec<_> = report
            .per_code
            .iter()
            .map(|r| (r.agree, r.disagree, r.present_in_gtf, r.present_in_cache))
            .collect();
        assert_eq!(
            rows,
            [(2, 1, 2, 1), (2, 1, 1, 0), (3, 0, 1, 1), (2, 1, 0, 1)]
        );
        let examples: Vec<_> = report
            .per_code
            .iter()
            .map(|r| {
                r.examples
                    .iter()
                    .map(|e| (e.transcript, e.gtf, e.cache))
                    .collect::<Vec<_>>()
            })
            .collect();
        assert_eq!(examples[0], [(&b"T3"[..], true, false)]);
        assert_eq!(examples[1], [(&b"T1"[..], true, false)]);
        assert!(examples[2].is_empty());
        assert_eq!(examples[3], [(&b"T2"[..], false, true)]);
    }

    #[test]
    fn ids_are_ordered_bytewise_at_any_length() {
        let ids = [
            "",
            "A",
            "AB",
            "AB\0",
            "AB\0\0C",
            "ABC",
            "ENST00000000001X",
            "ENST00000000001XA",
            "ENST00000000001XB",
            "ENST00000000001Y",
            "ENST0000000000",
            "T\u{e9}",
            "T\u{4e2d}",
            "T\u{1F600}",
            "Tz",
        ];
        let gtf: CodeSets = ids.iter().map(|id| (id.as_bytes().into(), 1)).collect();
        let cache: CodeSets = ids.iter().map(|id| (id.as_bytes().into(), 0)).collect();
        let report = compare(&gtf, &cache, 50);
        let ordered: Vec<&[u8]> = report.per_code[0]
            .examples
            .iter()
            .map(|e| e.transcript)
            .collect();
        let mut expected: Vec<&[u8]> = ids.iter().map(|id| id.as_bytes()).collect();
        expected.sort_unstable();
        assert_eq!(ordered, expected);
    }

    #[test]
    fn examples_stop_below_the_bound_and_a_bound_of_zero_or_less_keeps_none() {
        let gtf = sets(&[("A", 1), ("B", 1), ("C", 1)]);
        let cache = sets(&[("A", 0), ("B", 0), ("C", 0)]);
        for (bound, kept) in [
            (50, 3),
            (3, 3),
            (2, 2),
            (1, 1),
            (0, 0),
            (-1, 0),
            (i64::MAX, 3),
        ] {
            let report = compare(&gtf, &cache, bound);
            assert_eq!(report.per_code[0].examples.len(), kept, "bound {bound}");
            assert_eq!(report.per_code[0].disagree, 3);
        }
    }

    #[test]
    fn the_report_is_json_dump_with_an_indent_of_one() {
        let gtf = sets(&[("T1", 1), ("T\u{e9}\"2", 0)]);
        let cache = sets(&[("T1", 1), ("T\u{e9}\"2", 2)]);
        let report = compare(&gtf, &cache, 50);
        let escaped = format!("T{}u00e9{}\"2", '\\', '\\');
        let expected = [
            "{",
            " \"transcripts\": {",
            "  \"gtf\": 2,",
            "  \"cache\": 2,",
            "  \"shared\": 2,",
            "  \"gtf_only\": 0,",
            "  \"cache_only\": 0",
            " },",
            " \"per_code\": {",
            "  \"gencode_basic\": {",
            "   \"agree\": 2,",
            "   \"disagree\": 0,",
            "   \"present_in_gtf\": 1,",
            "   \"present_in_cache\": 1,",
            "   \"examples\": []",
            "  },",
            "  \"gencode_primary\": {",
            "   \"agree\": 1,",
            "   \"disagree\": 1,",
            "   \"present_in_gtf\": 0,",
            "   \"present_in_cache\": 1,",
            "   \"examples\": [",
            "    {",
            &format!("     \"transcript\": \"{escaped}\","),
            "     \"gtf\": false,",
            "     \"cache\": true",
            "    }",
            "   ]",
            "  },",
            "  \"cds_start_NF\": {",
            "   \"agree\": 2,",
            "   \"disagree\": 0,",
            "   \"present_in_gtf\": 0,",
            "   \"present_in_cache\": 0,",
            "   \"examples\": []",
            "  },",
            "  \"cds_end_NF\": {",
            "   \"agree\": 2,",
            "   \"disagree\": 0,",
            "   \"present_in_gtf\": 0,",
            "   \"present_in_cache\": 0,",
            "   \"examples\": []",
            "  }",
            " },",
            " \"exact_agreement\": false",
            "}",
            "",
        ]
        .join("\n");
        assert_eq!(json(&report), expected);
    }

    #[test]
    fn the_stdout_lines_group_digits_and_pad_as_the_format_specs_do() {
        let mut report = Report {
            gtf: 1_234_567,
            cache: 669_092,
            shared: 1_000,
            per_code: Default::default(),
        };
        report.per_code[0].agree = 12_345_678;
        report.per_code[0].disagree = 1_234_567;
        report.per_code[0].present_in_gtf = 12_345_678;
        report.per_code[0].present_in_cache = 999;
        report.per_code[1].agree = 1_000;
        for id in ["T1", "T2", "T3", "T4", "T5", "T6"] {
            report.per_code[1].examples.push(Example {
                transcript: id.as_bytes(),
                gtf: id != "T2",
                cache: id == "T2",
            });
        }
        let expected = "\
transcripts: gtf 1,234,567 cache 669,092 shared 1,000 gtf_only 1,233,567 cache_only 668,092
  gencode_basic    agree 12,345,678  disagree 1,234,567  present gtf 12,345,678 cache      999
  gencode_primary  agree     1,000  disagree       0  present gtf        0 cache        0
    T1: gtf=True cache=False
    T2: gtf=False cache=True
    T3: gtf=True cache=False
    T4: gtf=True cache=False
    T5: gtf=True cache=False
  cds_start_NF     agree         0  disagree       0  present gtf        0 cache        0
  cds_end_NF       agree         0  disagree       0  present gtf        0 cache        0
";
        assert_eq!(String::from_utf8(summary(&report)).unwrap(), expected);
    }
}
