// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The dumper against Data::Dumper's own text.
//!
//! `tests/fixtures/values/<name>.pst` is a Storable image written by perl 5.34.0 with
//! Storable 3.23, and `<name>.indent0.txt` and `<name>.indent1.txt` are what
//! `Data::Dumper::Dumper(retrieve("<name>.pst"))` returned on the same perl
//! (Data::Dumper 2.179) under `Sortkeys = 1`, `Terse = 1` and `Indent` 0 and 1. The
//! table below names what each image holds.

use std::path::PathBuf;

use vep_storable_digest::dumper::{dump, Indent, Overlay, Val};
use vep_storable_digest::shards::Shard;
use vep_storable_digest::storable::Image;

fn values() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/values")
}

fn check(name: &str) {
    let image = Image::read_file(&values().join(format!("{name}.pst"))).expect("an image");
    let shard = Shard::new(image);
    let root = Val::RefTo(shard.image.root());
    for (indent, suffix) in [(Indent::None, "indent0"), (Indent::Fixed, "indent1")] {
        let expected =
            std::fs::read(values().join(format!("{name}.{suffix}.txt"))).expect("expected text");
        let mut text = Vec::new();
        dump(&shard, &Overlay::default(), &root, indent, &mut text).expect("a text");
        assert!(
            text == expected,
            "{name} {suffix}\n--- expected\n{}\n--- got\n{}",
            String::from_utf8_lossy(&expected),
            String::from_utf8_lossy(&text)
        );
    }
}

macro_rules! value_tests {
    ($($test:ident => $name:literal,)*) => {
        $(#[test] fn $test() { check($name); })*

        /// Every image of the directory has a test above.
        #[test]
        fn every_image_is_tested() {
            let mut on_disk: Vec<String> = std::fs::read_dir(values())
                .expect("the fixtures")
                .filter_map(|entry| entry.ok()?.file_name().into_string().ok())
                .filter_map(|name| name.strip_suffix(".pst").map(str::to_string))
                .collect();
            on_disk.sort();
            let mut tested = vec![$($name.to_string(),)*];
            tested.sort();
            assert_eq!(on_disk, tested);
        }
    };
}

value_tests! {
    // A hash of scalars: integers of up to and past ten characters, strings of digits
    // (`'12345'`, `'007'`, `'-0'`, `' 12'`, `"12\n"`, `'1e3'`), doubles (which `nstore`
    // writes as strings), the empty string, undef, strings with `'`, `\`, `$`, `@`, `"`,
    // control characters and every byte from 128 to 255, and keys of the same kinds.
    scalars_are_quoted_by_the_kind_the_image_stores => "scalars",
    // `store` in the machine's own order: `{ integer => 42, big => 4294967296,
    // decimal => 1.5, whole_double => 3.0, tiny => 1e-7, huge_double => 1e20,
    // minus_zero => -0.0, fifteen => 0.1 + 0.2, list => [1, 2.5, 'x'], ... }`.
    native_order_integers_and_doubles_keep_their_kind => "native",
    // Values with Perl's UTF-8 flag: an ASCII string, `"caf\x{e9}"`, a string with
    // U+263A beside `$`, `@`, `"`, `\` and `'`, U+1F600; keys with U+263A, with a byte
    // 0xE9 and no flag, with U+00E9 under the flag, and an ASCII key under the flag.
    strings_and_keys_with_the_utf8_flag => "utf8",
    // `{ empty_hash => {}, empty_array => [], nested => [[], {}, [[1, [2, [3]]]],
    // { a => { b => { c => [] } } }], holes => <items 0 and 3 of four>,
    // undefs => [undef, 1, undef], one => [1], hash_one => { only => 1 } }`.
    containers_nest_and_stay_empty_or_hold_holes => "containers",
    // `[ 1, 'two', [3], { four => 4 } ]`.
    an_array_at_the_root => "array-root",
    // `\ 'just a string'`.
    a_scalar_at_the_root => "scalar-root",
    // A blessed hash, array and scalar, the classes `It's::Odd`, `main` and `0`,
    // an array blessed into `HASH`, and three classes nested.
    blessed_values_are_wrapped_in_bless => "blessed",
    // References to a string, a number, undef, a reference, a reference to a
    // reference, an array reference and a hash reference; two references to one
    // scalar; and references to two values of the hash itself, one under a key that
    // sorts before its value's key and one after.
    references_to_scalars_and_to_references => "scalar-refs",
    // A transcript whose two exons stand under `_trans_exon_array`, `again` (the same
    // array), `sorted` (the same exons in another array) and the translation's
    // `start_exon` and `end_exon`; the translation holds a weak reference to the
    // transcript, `self` a strong one; `ring` is an array that holds itself twice.
    a_value_reached_twice_is_the_path_of_its_first_place => "shared",
    // `{ kept => $kept, weak_kept => <weak $kept>, weak_lost => <weak $lost>,
    // list => [ <weak $lost>, 1 ] }`: nothing else holds `$lost`.
    a_weak_reference_to_a_freed_value_is_undef => "weak",
    // `sub { \@_ }->(!!1, !!0, undef, !!1, !!0, undef)`: the array holds Perl's
    // immortal true, false and undef themselves.
    an_immortal_met_twice_is_the_path_of_its_first_place => "immortals",
    // `[ <reference to the immortal true>, <a copy of it>, <reference to the immortal
    // undef>, undef ]`.
    references_to_immortals => "immortal-refs",
    // `sub { unshift @_, \$_[0]; \@_ }->(undef, undef)`: a reference to the immortal
    // undef, then that scalar twice.
    the_immortal_undef_after_a_reference_to_it => "immortal-undef-after-its-reference",
    // `{ plain => qr/abc/, flags => qr/a.b/msix, slash => qr{a/b\/c}, class => qr/[\/]x/,
    // reblessed => bless(qr/x+/, 'My::Pattern'), version => v1.22.333,
    // version_ref => \ v9.8 }`.
    compiled_patterns_and_version_strings => "patterns",
    // `{ locked => <a hash of keep, drop and also with its keys locked and drop deleted> }`.
    a_restricted_hash_without_its_deleted_key => "restricted",
    // A hash at the bottom of twelve levels of `{ down => [ ... ] }`, and again
    // beside them.
    a_path_through_many_levels => "deep",
}
