// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader against images Perl wrote. Every fixture under `tests/fixtures/storable/`
//! was written once by `perl -MStorable=store,nstore -e '<the line its test names>'
//! <fixture>` on perl 5.32.1 with Storable 3.37 (aarch64, little-endian, 8-byte
//! integers); the images a perl of this machine cannot write (another byte order,
//! 4-byte integers, the boolean opcodes, damaged images) are built in the tests.

// The reader is compiled into this test the way another tool compiles it.
#[allow(dead_code, unused_imports)]
#[path = "../src/storable/mod.rs"]
mod storable;

use std::path::PathBuf;

use storable::{ByteOrder, ErrorKind, Image, NodeId, Value};

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/storable")
        .join(name)
}

fn fixture(name: &str) -> Image {
    Image::read_file(&fixture_path(name)).unwrap_or_else(|error| panic!("{name}: {error}"))
}

fn fixture_bytes(name: &str) -> Vec<u8> {
    std::fs::read(fixture_path(name)).expect("fixture")
}

fn items(image: &Image, id: NodeId) -> Vec<NodeId> {
    image.items(id).expect("an array").to_vec()
}

fn text(image: &Image, id: NodeId) -> (Vec<u8>, bool) {
    let (bytes, utf8) = image.string(id).expect("a string scalar");
    (bytes.to_vec(), utf8)
}

/// A network-order header (`nstore`, binary 2.11) before `body`.
fn network_image(body: &[u8]) -> Vec<u8> {
    let mut bytes = b"pst0\x05\x0b".to_vec();
    bytes.extend_from_slice(body);
    bytes
}

fn kind_of(bytes: Vec<u8>) -> (ErrorKind, u64) {
    let error = Image::read(bytes).expect_err("an image the reader refuses");
    (error.kind, error.offset)
}

// nstore([5, -3, 300, -70000, 2**40, 1.5, "abc", "x" x 300, undef, ""], shift)
#[test]
fn network_order_scalars_keep_the_kind_they_were_stored_as() {
    let image = fixture("net_scalars.pst");
    let header = image.header();
    assert_eq!((header.major, header.minor), (2, 11));
    assert!(header.network_order);
    assert_eq!(header.byte_order, ByteOrder::Big);
    assert_eq!(header.native, None);
    assert_eq!(header.body_offset, 6);
    assert_eq!(image.root().tag(), 0);
    let values = items(&image, image.root());
    assert_eq!(values.len(), 10);
    // SX_BYTE twice, SX_NETINT twice.
    assert_eq!(image.value(values[0]), Value::Integer(5));
    assert_eq!(image.value(values[1]), Value::Integer(-3));
    assert_eq!(image.value(values[2]), Value::Integer(300));
    assert_eq!(image.value(values[3]), Value::Integer(-70000));
    // An integer past 32 bits and a double are strings in network order.
    assert_eq!(text(&image, values[4]), (b"1099511627776".to_vec(), false));
    assert_eq!(text(&image, values[5]), (b"1.5".to_vec(), false));
    // SX_SCALAR, SX_LSCALAR, SX_UNDEF and the empty SX_SCALAR.
    assert_eq!(text(&image, values[6]), (b"abc".to_vec(), false));
    assert_eq!(text(&image, values[7]), (vec![b'x'; 300], false));
    assert_eq!(image.value(values[8]), Value::Undef);
    assert!(!image.is_defined(values[8]));
    assert_eq!(text(&image, values[9]), (Vec::new(), false));
    assert!(image.is_defined(values[9]));
    assert_eq!(image.node_count(), 11);
    assert_eq!(image.trailing_len(), 0);
}

// store([5, -3, 300, -70000, 2**40, 1.5, "abc", "x" x 300, undef, ""], shift)
#[test]
fn native_order_scalars_are_integers_and_doubles() {
    let image = fixture("native_scalars.pst");
    let header = image.header();
    assert!(!header.network_order);
    assert_eq!(header.byte_order, ByteOrder::Little);
    assert_eq!(header.integer_size, 8);
    let native = header.native.as_ref().expect("a native layout");
    assert_eq!(native.byte_order_string, "12345678");
    assert_eq!(
        (
            native.int_size,
            native.long_size,
            native.pointer_size,
            native.double_size
        ),
        (4, 8, 8, Some(8))
    );
    let values = items(&image, image.root());
    assert_eq!(image.value(values[0]), Value::Integer(5));
    assert_eq!(image.value(values[1]), Value::Integer(-3));
    // SX_INTEGER three times, SX_DOUBLE once.
    assert_eq!(image.value(values[2]), Value::Integer(300));
    assert_eq!(image.value(values[3]), Value::Integer(-70000));
    assert_eq!(image.value(values[4]), Value::Integer(1_099_511_627_776));
    assert_eq!(image.value(values[5]), Value::Double(1.5));
    assert_eq!(text(&image, values[6]), (b"abc".to_vec(), false));
    assert_eq!(text(&image, values[7]), (vec![b'x'; 300], false));
    assert_eq!(image.value(values[8]), Value::Undef);
}

// nstore(["\x{263a}", "\x{263a}" x 100, "caf\xe9"], shift)
#[test]
fn strings_carry_the_utf8_flag_they_were_stored_with() {
    let image = fixture("utf8_strings.pst");
    let values = items(&image, image.root());
    // SX_UTF8STR and SX_LUTF8STR hold the UTF-8 bytes; a Latin-1 byte string is SX_SCALAR.
    assert_eq!(
        text(&image, values[0]),
        ("\u{263a}".as_bytes().to_vec(), true)
    );
    assert_eq!(
        text(&image, values[1]),
        ("\u{263a}".repeat(100).into_bytes(), true)
    );
    assert_eq!(text(&image, values[2]), (b"caf\xe9".to_vec(), false));
}

// sub f { nstore(\@_, $ARGV[0]) } f(undef, !!1, !!0)
#[test]
fn immortals_are_told_from_ordinary_scalars() {
    let image = fixture("immortals.pst");
    let values = items(&image, image.root());
    assert_eq!(image.value(values[0]), Value::ElementUndef);
    assert_eq!(image.value(values[1]), Value::Yes);
    assert_eq!(image.value(values[2]), Value::No);
    assert!(!image.is_defined(values[0]));
    assert!(image.is_defined(values[2]));
}

// my @a; $a[2] = 7; nstore(\@a, shift)
#[test]
fn an_array_element_that_does_not_exist_is_the_immortal_undef() {
    let image = fixture("array_holes.pst");
    let values = items(&image, image.root());
    assert_eq!(image.value(values[0]), Value::ImmortalUndef);
    assert_eq!(image.value(values[1]), Value::ImmortalUndef);
    assert_eq!(image.value(values[2]), Value::Integer(7));
    // Each stored item has a tag of its own.
    assert_ne!(values[0], values[1]);
}

// my $h = {k => 1}; my $s = "v"; nstore([$h, $h, \$s, \\$s, \1], shift)
#[test]
fn a_value_reached_twice_is_one_node() {
    let image = fixture("refs.pst");
    let values = items(&image, image.root());
    assert_eq!(values.len(), 5);
    let hash = image.target(values[0]).expect("a reference");
    assert!(matches!(
        image.value(hash),
        Value::Hash {
            restricted: false,
            ..
        }
    ));
    assert_eq!(
        image.value(image.fetch(hash, b"k").expect("key k")),
        Value::Integer(1)
    );
    // The second reference is a scalar of its own that points at the same hash.
    assert_ne!(values[0], values[1]);
    assert_eq!(image.target(values[1]), Some(hash));
    let scalar = image.target(values[2]).expect("a reference to a scalar");
    assert_eq!(text(&image, scalar), (b"v".to_vec(), false));
    let inner = image.target(values[3]).expect("a reference to a reference");
    assert_eq!(image.target(inner), Some(scalar));
    assert_eq!(
        image.value(image.target(values[4]).expect("a reference to a number")),
        Value::Integer(1)
    );
    assert_eq!(
        image.value(values[0]),
        Value::Ref {
            target: hash,
            weak: false,
            overloaded: false
        }
    );
}

// my $a = []; push @$a, $a; nstore($a, shift)
#[test]
fn a_structure_that_holds_itself_points_back_at_its_own_tag() {
    let image = fixture("cycle.pst");
    let values = items(&image, image.root());
    assert_eq!(values.len(), 1);
    assert_eq!(image.target(values[0]), Some(image.root()));
}

// nstore([bless({a => 1}, "Some::Class"), bless([], "Some::Class"),
//         bless(\(my $x = "s"), "Other"), bless({}, "L" x 200)], shift)
#[test]
fn blessing_is_seen_on_the_value_a_reference_points_at() {
    let image = fixture("blessed.pst");
    let values = items(&image, image.root());
    let targets: Vec<NodeId> = values
        .iter()
        .map(|value| image.target(*value).expect("a reference"))
        .collect();
    // SX_BLESS with a short name, SX_IX_BLESS for its second use, a blessed scalar,
    // and SX_BLESS with a name whose length takes four bytes.
    assert_eq!(image.class(targets[0]), Some(&b"Some::Class"[..]));
    assert!(matches!(image.value(targets[0]), Value::Hash { .. }));
    assert_eq!(image.class(targets[1]), Some(&b"Some::Class"[..]));
    assert!(matches!(image.value(targets[1]), Value::Array(_)));
    assert_eq!(image.class(targets[2]), Some(&b"Other"[..]));
    assert_eq!(text(&image, targets[2]), (b"s".to_vec(), false));
    assert_eq!(image.class(targets[3]), Some(&[b'L'; 200][..]));
    // The references and the unblessed root carry no class.
    assert!(values.iter().all(|value| image.class(*value).is_none()));
    assert_eq!(image.class(image.root()), None);
}

// nstore([map { bless({}, "C$_") } (0 .. 129, 0 .. 129)], shift)
#[test]
fn a_class_index_past_127_is_read_from_four_bytes() {
    let image = fixture("many_classes.pst");
    let values = items(&image, image.root());
    assert_eq!(values.len(), 260);
    for (index, value) in values.iter().enumerate() {
        let target = image.target(*value).expect("a reference");
        let expected = format!("C{}", index % 130);
        assert_eq!(
            image.class(target),
            Some(expected.as_bytes()),
            "item {index}"
        );
    }
}

// my $t = {n => 1}; my $r = [$t, $t]; weaken($r->[1]); nstore($r, shift)   (-MScalar::Util=weaken)
#[test]
fn a_weak_reference_is_marked_weak() {
    let image = fixture("weak.pst");
    let values = items(&image, image.root());
    let target = image.target(values[0]).expect("a reference");
    assert_eq!(
        image.value(values[0]),
        Value::Ref {
            target,
            weak: false,
            overloaded: false
        }
    );
    assert_eq!(
        image.value(values[1]),
        Value::Ref {
            target,
            weak: true,
            overloaded: false
        }
    );
}

// package O; use overload q("") => sub { "o" }; package main; my $o = bless {}, "O";
// my $r = [$o, $o]; weaken($r->[1]); nstore($r, shift)   (-MScalar::Util=weaken)
#[test]
fn a_reference_to_an_overloaded_object_is_marked_overloaded() {
    let image = fixture("overload.pst");
    let values = items(&image, image.root());
    let target = image.target(values[0]).expect("a reference");
    assert_eq!(image.class(target), Some(&b"O"[..]));
    assert_eq!(
        image.value(values[0]),
        Value::Ref {
            target,
            weak: false,
            overloaded: true
        }
    );
    assert_eq!(
        image.value(values[1]),
        Value::Ref {
            target,
            weak: true,
            overloaded: true
        }
    );
}

// my $k = "caf\xe9"; utf8::upgrade($k); nstore({"\x{263a}" => 1, $k => 2, plain => 3}, shift)
#[test]
fn a_flagged_hash_keeps_its_key_flags() {
    let image = fixture("flag_hash.pst");
    let root = image.root();
    assert!(matches!(
        image.value(root),
        Value::Hash {
            restricted: false,
            ..
        }
    ));
    let entries = image.entries(root).expect("a hash");
    assert_eq!(entries.len(), 3);
    let by_key = |key: &[u8]| {
        entries
            .iter()
            .find(|entry| image.key(entry) == key)
            .unwrap_or_else(|| panic!("key {key:?}"))
    };
    let plain = by_key(b"plain");
    assert_eq!(plain.flags.0, 0);
    assert_eq!(image.value(plain.value), Value::Integer(3));
    let wide = by_key("\u{263a}".as_bytes());
    assert!(wide.flags.utf8() && !wide.flags.was_utf8());
    assert_eq!(image.value(wide.value), Value::Integer(1));
    // A key that was UTF-8 and fits one byte a character is stored as those bytes.
    let downgraded = by_key(b"caf\xe9");
    assert!(downgraded.flags.was_utf8() && !downgraded.flags.utf8());
    assert_eq!(image.value(downgraded.value), Value::Integer(2));
    assert_eq!(image.fetch(root, b"plain"), Some(plain.value));
    assert_eq!(image.fetch(root, b"absent"), None);
}

// my %h = (a => 1, b => 2); Internals::SvREADONLY(%h, 1); delete $h{a};
// Internals::SvREADONLY($h{b}, 1); nstore(\%h, shift)
#[test]
fn a_restricted_hash_keeps_its_placeholders_and_locks() {
    let image = fixture("restricted_hash.pst");
    let root = image.root();
    assert!(matches!(
        image.value(root),
        Value::Hash {
            restricted: true,
            ..
        }
    ));
    let entries = image.entries(root).expect("a hash");
    assert_eq!(entries.len(), 2);
    let deleted = entries
        .iter()
        .find(|entry| image.key(entry) == b"a")
        .expect("the deleted key");
    assert!(deleted.flags.placeholder());
    let kept = entries
        .iter()
        .find(|entry| image.key(entry) == b"b")
        .expect("the kept key");
    assert!(kept.flags.locked() && !kept.flags.placeholder());
    assert_eq!(image.value(kept.value), Value::Integer(2));
    // Perl's fetch does not see a placeholder.
    assert_eq!(image.fetch(root, b"a"), None);
    assert_eq!(image.fetch(root, b"b"), Some(kept.value));
}

// nstore([v1.22.333], shift)
#[test]
fn a_version_string_keeps_its_literal() {
    let image = fixture("vstring.pst");
    let values = items(&image, image.root());
    assert_eq!(text(&image, values[0]), (vec![1, 22, 0xc5, 0x8d], true));
    assert_eq!(image.vstring(values[0]), Some(&b"v1.22.333"[..]));
    assert_eq!(image.vstring(image.root()), None);
}

// nstore([qr/ab+c/ix], shift)
#[test]
fn a_compiled_pattern_keeps_its_source_and_flags() {
    let image = fixture("regexp.pst");
    let values = items(&image, image.root());
    let pattern = image.target(values[0]).expect("a reference");
    assert_eq!(image.class(pattern), Some(&b"Regexp"[..]));
    assert_eq!(image.regexp(pattern), Some((&b"ab+c"[..], &b"ix"[..])));
}

// store({list => [1, 2, 3], name => "n"}, shift)
#[test]
fn native_order_lengths_follow_the_byte_order_string() {
    let image = fixture("native_hash.pst");
    let root = image.root();
    let list = image
        .target(image.fetch(root, b"list").expect("key list"))
        .expect("a reference");
    let numbers: Vec<Value> = items(&image, list)
        .iter()
        .map(|item| image.value(*item))
        .collect();
    assert_eq!(
        numbers,
        [Value::Integer(1), Value::Integer(2), Value::Integer(3)]
    );
    let name = image.fetch(root, b"name").expect("key name");
    assert_eq!(text(&image, name), (b"n".to_vec(), false));
}

// nstore({e => [], h => {}}, shift)
#[test]
fn empty_containers_have_no_items() {
    let image = fixture("empty_containers.pst");
    let root = image.root();
    let array = image
        .target(image.fetch(root, b"e").expect("key e"))
        .expect("a reference");
    assert_eq!(image.items(array), Some(&[][..]));
    let hash = image
        .target(image.fetch(root, b"h").expect("key h"))
        .expect("a reference");
    assert_eq!(image.entries(hash).map(<[_]>::len), Some(0));
}

// $Storable::Deparse = 1; nstore([sub { 1 }], shift)
// package TH; sub TIEHASH { bless {}, shift } package main; tie my %h, "TH"; nstore(\%h, shift)
// package TA; sub TIEARRAY { bless [], shift } sub FETCHSIZE { 0 } package main; tie my @a, "TA"; nstore(\@a, shift)
// package TS; sub TIESCALAR { bless [], shift } sub FETCH { 1 } package main; tie my $s, "TS"; nstore(\$s, shift)
// package H; sub STORABLE_freeze { "frozen" } sub STORABLE_thaw { } package main; nstore([bless {}, "H"], shift)
#[test]
fn code_tied_values_and_hooks_are_refused_by_opcode_and_offset() {
    for (name, opcode, opcode_name, offset) in [
        ("code.pst", 26u8, "SX_CODE", 12u64),
        ("tied_hash.pst", 12, "SX_TIED_HASH", 6),
        ("tied_array.pst", 11, "SX_TIED_ARRAY", 6),
        ("tied_scalar.pst", 13, "SX_TIED_SCALAR", 6),
        ("hook.pst", 19, "SX_HOOK", 12),
    ] {
        let error = Image::read_file(&fixture_path(name)).expect_err(name);
        assert_eq!(
            error.kind,
            ErrorKind::UnsupportedOpcode {
                opcode,
                name: opcode_name
            },
            "{name}"
        );
        assert_eq!(error.offset, offset, "{name}");
        assert_eq!(
            error.to_string(),
            format!("unsupported Storable opcode {opcode_name} ({opcode}) at offset {offset}")
        );
    }
}

// gzip -n -c net_scalars.pst > net_scalars.pst.gz
#[test]
fn a_gzip_compressed_image_reads_as_the_plain_one() {
    let plain = fixture("net_scalars.pst");
    let compressed = Image::read_gz_file(&fixture_path("net_scalars.pst.gz")).expect("gzip image");
    assert_eq!(compressed.node_count(), plain.node_count());
    assert_eq!(compressed.image_len(), plain.image_len());
    let values = items(&compressed, compressed.root());
    assert_eq!(compressed.value(values[2]), Value::Integer(300));
    assert_eq!(text(&compressed, values[6]), (b"abc".to_vec(), false));
}

#[test]
fn a_gzip_file_of_two_members_is_read_whole_and_the_first_image_decoded() {
    let member = fixture_bytes("net_scalars.pst.gz");
    let directory =
        std::env::temp_dir().join(format!("vep-storable-reader-{}", std::process::id()));
    std::fs::create_dir_all(&directory).expect("temporary directory");
    let path = directory.join("two-members.gz");
    std::fs::write(&path, [member.clone(), member].concat()).expect("write");
    let image = Image::read_gz_file(&path).expect("gzip image");
    std::fs::remove_dir_all(&directory).ok();
    let plain_len = fixture_bytes("net_scalars.pst").len();
    assert_eq!(image.image_len(), plain_len);
    assert_eq!(image.trailing_len(), plain_len);
}

#[test]
fn a_damaged_gzip_stream_and_a_missing_file_are_told_apart() {
    let directory = std::env::temp_dir().join(format!("vep-storable-gzip-{}", std::process::id()));
    std::fs::create_dir_all(&directory).expect("temporary directory");
    let member = fixture_bytes("net_scalars.pst.gz");
    let short = directory.join("short.gz");
    std::fs::write(&short, &member[..member.len() - 12]).expect("write");
    let error = Image::read_gz_file(&short).expect_err("a short gzip stream");
    assert!(matches!(error.kind, ErrorKind::Gzip(_)), "{error:?}");
    let error = Image::read_gz_file(&directory.join("absent.gz")).expect_err("no file");
    assert!(matches!(error.kind, ErrorKind::Io(_)), "{error:?}");
    std::fs::remove_dir_all(&directory).ok();
}

#[test]
fn bytes_after_the_image_are_left_unread() {
    let mut bytes = fixture_bytes("refs.pst");
    let len = bytes.len();
    bytes.extend_from_slice(b"junk");
    let image = Image::read(bytes).expect("image");
    assert_eq!(image.image_len(), len);
    assert_eq!(image.trailing_len(), 4);
}

#[test]
fn every_image_cut_short_ends_early_and_never_reads_as_whole() {
    for name in [
        "net_scalars.pst",
        "native_scalars.pst",
        "utf8_strings.pst",
        "refs.pst",
        "blessed.pst",
        "flag_hash.pst",
        "restricted_hash.pst",
        "vstring.pst",
        "regexp.pst",
        "native_hash.pst",
        "overload.pst",
    ] {
        let bytes = fixture_bytes(name);
        for cut in 0..bytes.len() {
            let (kind, _) = kind_of(bytes[..cut].to_vec());
            assert_eq!(kind, ErrorKind::Truncated, "{name} cut at {cut}");
        }
    }
}

#[test]
fn a_file_of_another_kind_is_not_a_storable() {
    let error = Image::read(b"not a Storable image".to_vec()).expect_err("text");
    assert_eq!(error.kind, ErrorKind::NotStorable);
    assert_eq!(error.to_string(), "File is not a perl storable");
    assert_eq!(
        kind_of(b"perl-store\x04".to_vec()).0,
        ErrorKind::NotStorable
    );
}

#[test]
fn a_binary_major_other_than_two_is_refused() {
    // Major 1 has no minor byte; major 3 is later than the format read here.
    assert_eq!(
        kind_of(b"pst0\x03\x05".to_vec()),
        (ErrorKind::UnsupportedVersion { major: 1, minor: 0 }, 4)
    );
    assert_eq!(
        kind_of(b"pst0\x07\x00\x05".to_vec()),
        (ErrorKind::UnsupportedVersion { major: 3, minor: 0 }, 4)
    );
}

#[test]
fn a_big_endian_native_image_reads_its_lengths_and_numbers_big_endian() {
    // store([300, 1.5], ...) as a big-endian perl with 8-byte integers writes it.
    let mut bytes = b"pst0\x04\x0b\x0887654321\x04\x08\x08\x08".to_vec();
    bytes.extend_from_slice(&[2, 0, 0, 0, 2]);
    bytes.push(6);
    bytes.extend_from_slice(&300i64.to_be_bytes());
    bytes.push(7);
    bytes.extend_from_slice(&1.5f64.to_be_bytes());
    let image = Image::read(bytes).expect("image");
    assert_eq!(image.header().byte_order, ByteOrder::Big);
    assert!(!image.header().network_order);
    let values = items(&image, image.root());
    assert_eq!(image.value(values[0]), Value::Integer(300));
    assert_eq!(image.value(values[1]), Value::Double(1.5));
}

#[test]
fn a_native_image_of_four_byte_integers_reads_them_as_four_bytes() {
    // store([-70000], ...) as a little-endian perl with 4-byte integers writes it,
    // before binary minor 2 added the size of a double to the header.
    let mut bytes = b"pst0\x04\x01\x041234\x04\x04\x04".to_vec();
    bytes.extend_from_slice(&[2, 1, 0, 0, 0]);
    bytes.push(6);
    bytes.extend_from_slice(&(-70000i32).to_le_bytes());
    let image = Image::read(bytes).expect("image");
    let header = image.header();
    assert_eq!(header.integer_size, 4);
    assert_eq!(header.native.as_ref().expect("native").double_size, None);
    let values = items(&image, image.root());
    assert_eq!(image.value(values[0]), Value::Integer(-70000));
}

#[test]
fn an_unknown_byte_order_string_is_refused() {
    let (kind, offset) = kind_of(b"pst0\x04\x0b\x043412\x04\x08\x08\x08\x05".to_vec());
    assert_eq!(kind, ErrorKind::UnsupportedByteOrder("3412".to_string()));
    assert_eq!(offset, 6);
}

#[test]
fn a_double_that_is_not_eight_bytes_is_refused_where_it_is_met() {
    let mut bytes = b"pst0\x04\x0b\x0812345678\x04\x08\x08\x10".to_vec();
    let offset = bytes.len() as u64;
    bytes.push(7);
    bytes.extend_from_slice(&[0; 16]);
    assert_eq!(
        kind_of(bytes),
        (ErrorKind::BadLength("double that is not 8 bytes"), offset)
    );
}

#[test]
fn a_back_reference_to_a_tag_not_yet_seen_is_refused() {
    let (kind, offset) = kind_of(network_image(&[0, 0, 0, 0, 5]));
    assert_eq!(kind, ErrorKind::UnknownTag(5));
    assert_eq!(offset, 6);
    let (kind, offset) = kind_of(network_image(&[18, 3, 5]));
    assert_eq!(kind, ErrorKind::UnknownClass(3));
    assert_eq!(offset, 6);
}

#[test]
fn an_opcode_the_format_does_not_define_is_refused_by_number() {
    let error = Image::read(network_image(&[2, 0, 0, 0, 1, 99])).expect_err("opcode 99");
    assert_eq!(
        error.kind,
        ErrorKind::UnsupportedOpcode {
            opcode: 99,
            name: "unknown"
        }
    );
    assert_eq!(
        error.to_string(),
        "unsupported Storable opcode unknown (99) at offset 11"
    );
}

#[test]
fn the_boolean_opcodes_of_binary_minor_12_are_scalars_of_their_own() {
    let mut bytes = b"pst0\x05\x0c".to_vec();
    bytes.extend_from_slice(&[2, 0, 0, 0, 2, 34, 35]);
    let image = Image::read(bytes).expect("image");
    let values = items(&image, image.root());
    assert_eq!(image.value(values[0]), Value::Bool(true));
    assert_eq!(image.value(values[1]), Value::Bool(false));
    assert!(image.is_defined(values[1]));
}

#[test]
fn a_class_named_before_a_back_reference_is_dropped() {
    // [ "a", <SX_BLESS "C"> <SX_OBJECT 1> ]: the second item is the first one again,
    // and Storable blesses nothing.
    let image = Image::read(network_image(&[
        2, 0, 0, 0, 2, 10, 1, b'a', 17, 1, b'C', 0, 0, 0, 0, 1,
    ]))
    .expect("image");
    let values = items(&image, image.root());
    assert_eq!(values[0], values[1]);
    assert_eq!(image.class(values[0]), None);
}

#[test]
fn a_key_stored_twice_fetches_the_later_value() {
    // { k => 1, k => 2 } as a writer that repeats a key would store it.
    let image = Image::read(network_image(&[
        3, 0, 0, 0, 2, 8, 129, 0, 0, 0, 1, b'k', 8, 130, 0, 0, 0, 1, b'k',
    ]))
    .expect("image");
    let root = image.root();
    assert_eq!(image.entries(root).map(<[_]>::len), Some(2));
    assert_eq!(
        image.value(image.fetch(root, b"k").expect("key k")),
        Value::Integer(2)
    );
}

#[test]
fn a_large_object_is_read_as_a_string_or_a_tag_and_refused_as_a_container() {
    let mut body = vec![2, 0, 0, 0, 3, 33, 1];
    body.extend_from_slice(&3u64.to_be_bytes());
    body.extend_from_slice(b"abc");
    body.extend_from_slice(&[33, 24]);
    body.extend_from_slice(&3u64.to_be_bytes());
    body.extend_from_slice("\u{263a}".as_bytes());
    body.extend_from_slice(&[33, 0]);
    body.extend_from_slice(&1u64.to_be_bytes());
    let image = Image::read(network_image(&body)).expect("image");
    let values = items(&image, image.root());
    assert_eq!(text(&image, values[0]), (b"abc".to_vec(), false));
    assert_eq!(
        text(&image, values[1]),
        ("\u{263a}".as_bytes().to_vec(), true)
    );
    assert_eq!(values[2], values[0]);

    let mut body = vec![33, 2];
    body.extend_from_slice(&0u64.to_be_bytes());
    assert_eq!(
        kind_of(network_image(&body)),
        (
            ErrorKind::UnsupportedOpcode {
                opcode: 33,
                name: "SX_LOBJECT"
            },
            6
        )
    );
}

#[test]
fn a_hash_key_held_as_a_scalar_is_refused() {
    let (kind, offset) = kind_of(network_image(&[25, 0, 0, 0, 0, 1, 8, 129, 8, 10, 1, b'k']));
    assert_eq!(kind, ErrorKind::ScalarHashKey);
    assert_eq!(offset, 14);
}

#[test]
fn references_nested_two_hundred_thousand_deep_read_without_recursion() {
    let depth = 200_000;
    let mut body = vec![4u8; depth];
    body.push(5);
    let image = Image::read(network_image(&body)).expect("image");
    assert_eq!(image.node_count(), depth + 1);
    let mut node = image.root();
    let mut hops = 0;
    while let Some(target) = image.target(node) {
        node = target;
        hops += 1;
    }
    assert_eq!(hops, depth);
    assert_eq!(image.value(node), Value::Undef);
}

#[test]
fn a_count_past_the_bytes_left_ends_early_before_anything_is_reserved() {
    let (kind, _) = kind_of(network_image(&[2, 0x7f, 0xff, 0xff, 0xff, 5]));
    assert_eq!(kind, ErrorKind::Truncated);
    // A negative count is an empty container, as Storable reads it.
    let image = Image::read(network_image(&[2, 0xff, 0xff, 0xff, 0xff])).expect("image");
    assert_eq!(image.items(image.root()), Some(&[][..]));
}

/// A native-order header of a little-endian perl with 8-byte integers before `body`.
fn little_endian_image(body: &[u8]) -> Vec<u8> {
    let mut bytes = b"pst0\x04\x0b\x0812345678\x04\x08\x08\x08".to_vec();
    bytes.extend_from_slice(body);
    bytes
}

#[test]
fn a_string_length_past_the_bytes_left_ends_early() {
    // SX_SCALAR and SX_UTF8STR carry one length byte, SX_LSCALAR and SX_LUTF8STR four.
    for body in [
        vec![10, 200, b'a', b'b'],
        vec![23, 200, b'a', b'b'],
        vec![1, 0x7f, 0xff, 0xff, 0xff, b'a'],
        vec![1, 0xff, 0xff, 0xff, 0xff, b'a'],
        vec![24, 0xff, 0xff, 0xff, 0xff, b'a'],
    ] {
        let image = network_image(&body);
        let len = image.len() as u64;
        assert_eq!(kind_of(image), (ErrorKind::Truncated, len), "{body:?}");
    }
    let image = little_endian_image(&[1, 0xff, 0xff, 0xff, 0x7f, b'a']);
    let len = image.len() as u64;
    assert_eq!(kind_of(image), (ErrorKind::Truncated, len));
}

#[test]
fn a_hash_count_past_the_bytes_left_ends_early_before_anything_is_reserved() {
    for body in [
        vec![3, 0x7f, 0xff, 0xff, 0xff, 5],
        // SX_FLAG_HASH: the hash flags, then the count.
        vec![25, 0, 0x7f, 0xff, 0xff, 0xff, 5],
        // Two entries are promised and the bytes hold one.
        vec![3, 0, 0, 0, 2, 5, 0, 0, 0, 1, b'k'],
    ] {
        assert_eq!(
            kind_of(network_image(&body)).0,
            ErrorKind::Truncated,
            "{body:?}"
        );
    }
    assert_eq!(
        kind_of(little_endian_image(&[2, 0xff, 0xff, 0xff, 0x7f, 5])).0,
        ErrorKind::Truncated
    );
    // A negative count is an empty hash, as Storable reads it.
    let image = Image::read(network_image(&[3, 0xff, 0xff, 0xff, 0xff])).expect("image");
    assert_eq!(image.entries(image.root()).map(<[_]>::len), Some(0));
}

#[test]
fn a_hash_key_length_past_the_bytes_left_ends_early_and_a_negative_one_is_refused() {
    // { <undef> => key }: the value, then the key's four-byte length and its bytes.
    let entry = |length: [u8; 4]| {
        let mut body = vec![3, 0, 0, 0, 1, 5];
        body.extend_from_slice(&length);
        body.push(b'k');
        network_image(&body)
    };
    assert_eq!(
        kind_of(entry([0x7f, 0xff, 0xff, 0xff])).0,
        ErrorKind::Truncated
    );
    assert_eq!(kind_of(entry([0, 0, 0, 2])).0, ErrorKind::Truncated);
    // The key starts behind the header, the opcode, the count and the value.
    for length in [[0xff, 0xff, 0xff, 0xff], [0x80, 0, 0, 0]] {
        assert_eq!(
            kind_of(entry(length)),
            (ErrorKind::BadLength("negative hash key length"), 12),
            "{length:?}"
        );
    }
}

#[test]
fn a_class_name_length_past_the_bytes_left_ends_early() {
    for body in [
        vec![17, 5, b'C'],
        // A length byte with its top bit set is followed by a four-byte length.
        vec![17, 0x80, 0x7f, 0xff, 0xff, 0xff, b'C'],
        vec![17, 0x80, 0xff, 0xff, 0xff, 0xff, b'C'],
    ] {
        assert_eq!(
            kind_of(network_image(&body)).0,
            ErrorKind::Truncated,
            "{body:?}"
        );
    }
}

#[test]
fn a_class_index_of_four_bytes_that_no_class_carries_is_refused() {
    // [ bless(<undef>, "C"), <SX_IX_BLESS 1> ... ]: one class is known, index 0.
    let (kind, offset) = kind_of(network_image(&[
        2, 0, 0, 0, 2, 17, 1, b'C', 5, 18, 0x80, 0, 0, 0, 1, 5,
    ]));
    assert_eq!(kind, ErrorKind::UnknownClass(1));
    assert_eq!(offset, 15);
    let (kind, _) = kind_of(network_image(&[
        2, 0, 0, 0, 2, 17, 1, b'C', 5, 18, 0x80, 0xff, 0xff, 0xff, 0xff, 5,
    ]));
    assert_eq!(kind, ErrorKind::UnknownClass(u64::from(u32::MAX)));
}

#[test]
fn a_back_reference_names_a_value_met_before_it_and_no_other() {
    // [ <SX_OBJECT tag> ]: the array is tag 0 and the only value met.
    let item = |tag: u32| {
        let mut body = vec![2, 0, 0, 0, 1, 0];
        body.extend_from_slice(&tag.to_be_bytes());
        network_image(&body)
    };
    for tag in [1, 2, u32::MAX] {
        assert_eq!(
            kind_of(item(tag)),
            (ErrorKind::UnknownTag(u64::from(tag)), 11),
            "tag {tag}"
        );
    }
    // Tag 0 is the array itself, still open: the item is the array.
    let image = Image::read(item(0)).expect("image");
    assert_eq!(items(&image, image.root()), [image.root()]);
    assert_eq!(image.node_count(), 1);
}

#[test]
fn a_large_object_tag_past_32_bits_is_refused_whole() {
    // [ <undef>, <SX_LOBJECT SX_OBJECT tag> ]: a tag whose low 32 bits are those of a
    // value met before is still a tag no value carries.
    let item = |tag: u64| {
        let mut body = vec![2, 0, 0, 0, 2, 5, 33, 0];
        body.extend_from_slice(&tag.to_be_bytes());
        network_image(&body)
    };
    for tag in [1 << 32, (1 << 32) + 1, u64::MAX] {
        assert_eq!(
            kind_of(item(tag)),
            (ErrorKind::UnknownTag(tag), 12),
            "tag {tag}"
        );
    }
}

#[test]
fn a_large_object_string_of_4_gib_or_more_is_refused_before_its_bytes_are_looked_for() {
    for kind in [1u8, 24] {
        for length in [1u64 << 32, u64::MAX] {
            let mut body = vec![33, kind];
            body.extend_from_slice(&length.to_be_bytes());
            body.push(b'a');
            assert_eq!(
                kind_of(network_image(&body)),
                (ErrorKind::BadLength("string of 4 GiB or more"), 6),
                "kind {kind}, length {length}"
            );
        }
    }
    // A length that fits 32 bits and passes the bytes left ends early.
    let mut body = vec![33, 1];
    body.extend_from_slice(&u64::from(u32::MAX).to_be_bytes());
    body.push(b'a');
    assert_eq!(kind_of(network_image(&body)).0, ErrorKind::Truncated);
}

#[test]
fn a_version_string_or_pattern_length_past_the_bytes_left_ends_early() {
    for body in [
        // SX_VSTRING and SX_LVSTRING.
        vec![29, 200, b'v'],
        vec![30, 0x7f, 0xff, 0xff, 0xff, b'v'],
        // SX_REGEXP: a pattern length of one byte, of four bytes, and the flags' length.
        vec![32, 0, 200, b'a'],
        vec![32, 1, 0x7f, 0xff, 0xff, 0xff, b'a'],
        vec![32, 0, 1, b'a', 200, b'i'],
    ] {
        assert_eq!(
            kind_of(network_image(&body)).0,
            ErrorKind::Truncated,
            "{body:?}"
        );
    }
}

#[test]
fn a_reference_that_points_at_itself_is_one_node() {
    // <SX_REF> <SX_OBJECT 0>: a walk that follows targets must count its steps.
    let image = Image::read(network_image(&[4, 0, 0, 0, 0, 0])).expect("image");
    assert_eq!(image.node_count(), 1);
    assert_eq!(image.target(image.root()), Some(image.root()));
    assert_eq!(image.trailing_len(), 0);
}

#[test]
fn arrays_and_hashes_nested_two_hundred_thousand_deep_read_without_recursion() {
    let depth = 200_000;
    // [[[ ... <undef> ... ]]]: SX_ARRAY and a count of one, `depth` times.
    let mut body = Vec::with_capacity(depth * 5 + 1);
    for _ in 0..depth {
        body.extend_from_slice(&[2, 0, 0, 0, 1]);
    }
    body.push(5);
    let image = Image::read(network_image(&body)).expect("image");
    assert_eq!(image.node_count(), depth + 1);
    let mut node = image.root();
    for _ in 0..depth {
        let [item] = image.items(node).expect("an array") else {
            panic!("an array of one item");
        };
        node = *item;
    }
    assert_eq!(image.value(node), Value::Undef);

    // {k => {k => ... <undef> ...}}: every key follows its value, so the keys close the
    // image innermost first.
    let mut body = Vec::with_capacity(depth * 10 + 1);
    for _ in 0..depth {
        body.extend_from_slice(&[3, 0, 0, 0, 1]);
    }
    body.push(5);
    for _ in 0..depth {
        body.extend_from_slice(&[0, 0, 0, 1, b'k']);
    }
    let image = Image::read(network_image(&body)).expect("image");
    assert_eq!(image.node_count(), depth + 1);
    assert_eq!(image.trailing_len(), 0);
    let mut node = image.root();
    for _ in 0..depth {
        node = image.fetch(node, b"k").expect("key k");
    }
    assert_eq!(image.value(node), Value::Undef);
}
