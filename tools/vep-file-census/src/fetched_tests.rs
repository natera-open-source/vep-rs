// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

#[test]
fn the_part_sizes_tried_are_the_usual_ones_first_and_twelve_at_most() {
    // Three parts of 20 MiB: 7 to 9 MiB fit, the usual 8 first.
    assert_eq!(part_sizes(20 * MIB, 3), [8, 7, 9]);
    // One part: the first size that fits alone.
    assert_eq!(part_sizes(20 * MIB, 1), [32]);
    assert_eq!(part_sizes(3 * MIB, 1), [8]);
    // Two parts of 6 GiB + 1 fit many sizes: twelve are tried.
    let sizes = part_sizes(6144 * MIB + 1, 2);
    assert_eq!(sizes.len(), 12);
    assert_eq!(sizes[..3], [3073, 3074, 3075]);
    assert_eq!(part_sizes(0, 1), [0u128; 0]);
    assert_eq!(part_sizes(20 * MIB, 4), [5, 6]);
    assert_eq!(part_sizes(20 * MIB, 6), [0u128; 0]);
}

// "a\nb\r\nc\rd\x0be\x1cf g".splitlines()
#[test]
fn the_list_is_cut_into_lines_as_splitlines_cuts_it() {
    assert_eq!(
        split_lines("a\nb\r\nc\rd\u{b}e\u{1c}f\u{2028}g"),
        ["a", "b", "c", "d", "e", "f", "g"]
    );
    assert_eq!(split_lines("a\n\nb\n"), ["a", "", "b"]);
    assert_eq!(split_lines(""), [""; 0]);
}

#[test]
fn the_closing_line_names_sha256_only_where_a_key_is_held_so() {
    let words =
        |list: &[&str]| -> Vec<String> { list.iter().map(|word| (*word).to_string()).collect() };
    assert_eq!(
        held_line(&words(&[
            "md5",
            "multipart-md5 8",
            "multipart-md5 16",
            BYTES_ALONE,
            DIFFERS
        ])),
        "held: 1 by md5, 2 by multipart md5, 1 by bytes alone, 1 differ of 5"
    );
    assert_eq!(
        held_line(&words(&["sha256", "md5"])),
        "held: 1 by md5, 0 by multipart md5, 0 by bytes alone, 1 by sha256, 0 differ of 2"
    );
}

// { head -c 1048576 F | md5sum; tail -c +1048577 F | md5sum; } | cut -c1-32 | xxd -r -p | md5sum
#[test]
fn the_etag_of_parts_is_the_md5_of_the_parts_md5s_and_their_count() {
    let dir =
        std::env::temp_dir().join(format!("vep-file-census-unit-{}-parts", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("object");
    let bytes: Vec<u8> = (0..(MIB as u32 + 4096))
        .map(|at| (at % 253) as u8)
        .collect();
    std::fs::write(&path, &bytes).unwrap();
    let (first, second) = bytes.split_at(MIB as usize);
    let joined = [md5::compute(first).0, md5::compute(second).0].concat();
    assert_eq!(
        etag_of_parts(&path, MIB).unwrap(),
        format!("{:x}-2", md5::compute(joined))
    );
    assert_eq!(
        etag_of_parts(&path, 8 * MIB).unwrap(),
        format!("{:x}-1", md5::compute(md5::compute(&bytes).0))
    );
    std::fs::write(&path, b"").unwrap();
    assert_eq!(
        etag_of_parts(&path, MIB).unwrap(),
        format!("{:x}-0", md5::compute(b""))
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
