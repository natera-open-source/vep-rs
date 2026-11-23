// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn rules(mode: Mode, contigs: Option<ContigSet>) -> Rules {
    Rules {
        mode,
        contigs: contigs.map(Arc::new),
    }
}

fn human() -> Option<ContigSet> {
    Some(ContigSet::human())
}

/// What `process` writes for `unit`, its counts and whether a cap stopped it.
fn processed(unit: &[u8], rules: &Rules, cap: Option<u64>) -> (Vec<u8>, Counts, bool) {
    let mut counts = Counts::default();
    let (text, stopped) = process(unit, rules, cap, &mut counts);
    let bytes = match text {
        Text::Same(kept) => unit[..kept].to_vec(),
        Text::Rewritten(bytes) => bytes,
    };
    (bytes, counts, stopped)
}

fn counts(headers: u64, variants: u64, dropped: u64, chr_stripped: bool) -> Counts {
    Counts {
        headers,
        variants,
        dropped,
        chr_stripped,
    }
}

#[test]
fn lines_that_stay_as_read_are_not_rewritten() {
    let unit = b"##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t100\tA\n22\t5\nX\n";
    let mut seen = Counts::default();
    let (text, stopped) = process(unit, &rules(Mode::Keep, human()), None, &mut seen);
    assert_eq!(text, Text::Same(unit.len()));
    assert!(!stopped);
    assert_eq!(seen, counts(2, 3, 0, false));
    // A header without a line end at the end of the input stays without one.
    let unit = b"1\t100\n#CHROM";
    let (text, _) = process(unit, &rules(Mode::Keep, None), None, &mut seen);
    assert_eq!(text, Text::Same(unit.len()));
}

#[test]
fn the_first_column_goes_through_the_mode_and_then_the_filter() {
    let unit = b"chr1\t1\nCHR2\t2\n3\t3\nchrUn\t4\nchr\t5\n\t6\nchrchr1\t7\nGL1\t8\n";
    let cases: [(Mode, Option<ContigSet>, &[u8], Counts); 7] = [
        (Mode::Keep, None, unit, counts(0, 8, 0, false)),
        (
            Mode::Strip,
            None,
            b"1\t1\n2\t2\n3\t3\nUn\t4\n\t5\n\t6\nchr1\t7\nGL1\t8\n",
            counts(0, 8, 0, false),
        ),
        (
            Mode::Add,
            None,
            b"chr1\t1\nCHR2\t2\nchr3\t3\nchrUn\t4\nchr\t5\nchr\t6\nchrchr1\t7\nchrGL1\t8\n",
            counts(0, 8, 0, false),
        ),
        // Under the filter a name that still starts with `chr` loses it once, whatever
        // the mode, and what is left is looked up.
        (
            Mode::Keep,
            human(),
            b"1\t1\n2\t2\n3\t3\n",
            counts(0, 3, 5, true),
        ),
        (
            Mode::Strip,
            human(),
            b"1\t1\n2\t2\n3\t3\n1\t7\n",
            counts(0, 4, 4, true),
        ),
        (
            Mode::Add,
            human(),
            b"1\t1\n2\t2\n3\t3\n",
            counts(0, 3, 5, true),
        ),
        // A set that names a prefixed contig keeps that spelling.
        (
            Mode::Keep,
            Some(ContigSet::named(b"chr1,2,chr,chr1\t").unwrap()),
            b"chr1\t1\n2\t2\nchr\t5\nchr1\t7\n",
            counts(0, 4, 4, true),
        ),
    ];
    for (mode, contigs, expected, expected_counts) in cases {
        let (bytes, seen, _) = processed(unit, &rules(mode, contigs), None);
        assert_eq!(
            String::from_utf8_lossy(&bytes),
            String::from_utf8_lossy(expected),
            "{mode:?}"
        );
        assert_eq!(seen, expected_counts, "{mode:?}");
    }
}

#[test]
fn every_line_end_is_written_as_a_line_feed_and_a_variant_line_always_has_one() {
    let unit = b"##a\r\n##b\r#c\n1\t1\r\n1\t2\r2\t3\n\n\r\n3\t4";
    let (bytes, seen, _) = processed(unit, &rules(Mode::Keep, None), None);
    assert_eq!(bytes, b"##a\n##b\n#c\n1\t1\n1\t2\n2\t3\n\n\n3\t4\n");
    assert_eq!(seen, counts(3, 6, 0, false));
    // A line feed after a carriage return is its line end, not a line of its own; a
    // carriage return that ends the run is a line end.
    let (bytes, seen, _) = processed(b"1\r\r\n2\r", &rules(Mode::Keep, None), None);
    assert_eq!(bytes, b"1\n\n2\n");
    assert_eq!(seen.variants, 3);
    // An empty line is a variant line with an empty first column: dropped under a set,
    // named `chr` by the add mode.
    let (bytes, seen, _) = processed(b"\n1\t1\n", &rules(Mode::Keep, human()), None);
    assert_eq!(bytes, b"1\t1\n");
    assert_eq!(seen, counts(0, 1, 1, false));
    let (bytes, _, _) = processed(b"\n", &rules(Mode::Add, None), None);
    assert_eq!(bytes, b"chr\n");
    let (bytes, seen, _) = processed(b"\n", &rules(Mode::Add, human()), None);
    assert_eq!(bytes, b"");
    assert_eq!(seen, counts(0, 0, 1, true));
}

#[test]
fn bytes_outside_utf8_are_written_as_replacement_characters() {
    let unit = b"#h\xff\n1\t\xe2\x82x\xf0\x9f\x98\n\xff\t2\nchr\xc0\x80\t3\ncaf\xc3\xa9\t4\n";
    let (bytes, seen, _) = processed(unit, &rules(Mode::Keep, None), None);
    assert_eq!(
        String::from_utf8(bytes).unwrap(),
        "#h\u{fffd}\n1\t\u{fffd}x\u{fffd}\n\u{fffd}\t2\nchr\u{fffd}\u{fffd}\t3\ncaf\u{e9}\t4\n"
    );
    assert_eq!(seen, counts(1, 4, 0, false));
    // A set can name the replacement character a stray byte is read as.
    let set = ContigSet::named("\u{fffd},caf\u{e9}".as_bytes()).unwrap();
    let (bytes, seen, _) = processed(unit, &rules(Mode::Keep, Some(set)), None);
    assert_eq!(
        String::from_utf8(bytes).unwrap(),
        "#h\u{fffd}\n\u{fffd}\t2\ncaf\u{e9}\t4\n"
    );
    assert_eq!(seen, counts(1, 2, 2, true));
    // A byte order mark is not a `#`: the line it starts is a variant line.
    let (_, seen, _) = processed(b"\xef\xbb\xbf##a\n#b\n", &rules(Mode::Keep, None), None);
    assert_eq!(seen, counts(1, 1, 0, false));
}

#[test]
fn a_cap_stops_before_the_first_variant_line_once_it_is_reached() {
    let unit = b"#a\n1\t1\n#b\nGL1\t2\n2\t3\n#c\n#d\n3\t4\n#e\n";
    let keep = rules(Mode::Keep, human());
    // Headers after the last written variant are written until a variant line is met,
    // and a dropped line does not count towards the cap.
    let cases: [(u64, &[u8], Counts, bool); 4] = [
        (0, b"#a\n", counts(1, 0, 0, false), true),
        (1, b"#a\n1\t1\n#b\n", counts(2, 1, 0, false), true),
        (
            2,
            b"#a\n1\t1\n#b\n2\t3\n#c\n#d\n",
            counts(4, 2, 1, false),
            true,
        ),
        (
            3,
            b"#a\n1\t1\n#b\n2\t3\n#c\n#d\n3\t4\n#e\n",
            counts(5, 3, 1, false),
            false,
        ),
    ];
    for (cap, expected, expected_counts, expected_stop) in cases {
        let (bytes, seen, stopped) = processed(unit, &keep, Some(cap));
        assert_eq!(bytes, expected, "cap {cap}");
        assert_eq!(seen, expected_counts, "cap {cap}");
        assert_eq!(stopped, expected_stop, "cap {cap}");
    }
    // Without a rewrite the kept bytes are the start of the run.
    let mut seen = Counts::default();
    let (text, stopped) = process(
        b"#a\n1\t1\n2\t2\n",
        &rules(Mode::Keep, None),
        Some(1),
        &mut seen,
    );
    assert_eq!((text, stopped), (Text::Same(7), true));
}

/// A deterministic stream of numbers for the cuts below.
struct Sequence(u64);

impl Sequence {
    fn below(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }
}

/// Feeds `pieces` to a collector, through `split` where `splitting`, and returns what it
/// emits and counts.
fn collected(
    pieces: &[&[u8]],
    rules: &Rules,
    cap: Option<u64>,
    splitting: bool,
) -> (Vec<u8>, Counts) {
    let mut written = Vec::new();
    let emit = |chunk: Chunk| -> io::Result<()> {
        assert!(!chunk.bytes().is_empty(), "an empty chunk is not emitted");
        written.extend_from_slice(chunk.bytes());
        Ok(())
    };
    let mut collector = Collector::new(rules, cap, emit);
    let mut stopped = false;
    for piece in pieces {
        if splitting {
            match split(piece.to_vec(), rules) {
                Ok(split) => collector.split(split).unwrap(),
                Err(piece) => stopped = collector.raw(piece).unwrap(),
            }
        } else {
            stopped = collector.raw(piece.to_vec()).unwrap();
        }
        if stopped {
            break;
        }
    }
    if !stopped {
        collector.finish().unwrap();
    }
    let counts = collector.counts;
    drop(collector);
    (written, counts)
}

const INPUTS: [&[u8]; 8] = [
    b"",
    b"##fileformat=VCFv4.2\n#CHROM\tPOS\nchr1\t1\tx\nchr1\t2\ty\nchr2_alt\t3\tz\n2\t4\nCHRX\t5\n",
    b"##a\r\n##b\r#c\n1\t1\r\n1\t2\r2\t3\n\n\r\n3\t4\r\r\n\rchr5\t5\r",
    b"1\t1\n#late header\nGL1\t2\n\t3\nno_tab\n#h\xff\r\n\xff\t4\nchr\xe2\x82\t5\ncaf\xc3\xa9\t6",
    b"\xef\xbb\xbf##a\n#b\nchrM\t1\nchrMT\t2\nM\t3\n",
    b"\r\r\r\n\n\r",
    b"#only\n#headers\r\n#here",
    b"chr1\t1\tA\tC\nchr1\t2\tA\tC\nchr1\t3\tA\tC\nchr1\t4\tA\tC\nchr1\t5\tA\tC\nchr1\t6\tA\tC\n",
];

fn rule_sets() -> Vec<Rules> {
    vec![
        rules(Mode::Keep, None),
        rules(Mode::Keep, human()),
        rules(Mode::Strip, None),
        rules(Mode::Add, None),
        rules(Mode::Add, human()),
        rules(Mode::Keep, Some(ContigSet::named(b"chr1,2,X,chr").unwrap())),
    ]
}

fn cut_in_two(input: &[u8], cut: usize) -> Vec<&[u8]> {
    [&input[..cut], &input[cut..]]
        .into_iter()
        .filter(|piece| !piece.is_empty())
        .collect()
}

#[test]
fn an_input_cut_anywhere_is_written_as_the_whole_input_is() {
    for input in INPUTS {
        for rules in rule_sets() {
            for cap in [None, Some(0), Some(1), Some(2), Some(4), Some(1000)] {
                let whole = collected(&[input], &rules, cap, false);
                for cut in 0..=input.len() {
                    assert_eq!(
                        collected(&cut_in_two(input, cut), &rules, cap, false),
                        whole,
                        "cut at {cut} of {} bytes, cap {cap:?}",
                        input.len()
                    );
                }
                let single: Vec<&[u8]> = input.chunks(1).collect();
                assert_eq!(
                    collected(&single, &rules, cap, false),
                    whole,
                    "byte by byte"
                );
            }
        }
    }
}

#[test]
fn pieces_whose_whole_lines_are_processed_apart_join_to_the_same_bytes() {
    let mut sequence = Sequence(0x9e37_79b9_7f4a_7c15);
    for input in INPUTS {
        for rules in rule_sets() {
            let whole = collected(&[input], &rules, None, false);
            for cut in 0..=input.len() {
                assert_eq!(
                    collected(&cut_in_two(input, cut), &rules, None, true),
                    whole,
                    "cut at {cut} of {} bytes",
                    input.len()
                );
            }
            for _ in 0..40 {
                let mut pieces: Vec<&[u8]> = Vec::new();
                let mut rest = input;
                while !rest.is_empty() {
                    let (piece, after) = rest.split_at(1 + sequence.below(rest.len().min(24)));
                    pieces.push(piece);
                    rest = after;
                }
                assert_eq!(collected(&pieces, &rules, None, true), whole, "{pieces:?}");
                assert_eq!(collected(&pieces, &rules, None, false), whole, "{pieces:?}");
            }
        }
    }
}

#[test]
fn a_piece_without_a_line_feed_is_not_split() {
    let keep = rules(Mode::Keep, None);
    assert!(split(b"chr1\t1\r2\t2".to_vec(), &keep).is_err());
    let piece = split(b"tail of a line\n1\t1\n2\t2\nhead of".to_vec(), &keep).unwrap();
    assert_eq!((piece.head_end, piece.tail_start), (15, 23));
    assert_eq!(piece.middle, None);
    assert_eq!(piece.counts.variants, 2);
    let piece = split(b"x\nchr1\t1\n".to_vec(), &rules(Mode::Strip, None)).unwrap();
    assert_eq!(piece.middle.as_deref(), Some(&b"1\t1\n"[..]));
}
