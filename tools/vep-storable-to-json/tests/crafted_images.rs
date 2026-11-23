// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The conversion of a transcript on images built here: shapes no cache holds and a
//! damaged file can, where a copied field is a structure that nests past what JSON::XS
//! writes or that holds itself, a biotype held as an object, and a slice whose length
//! is damaged or is a double.

use vep_storable_to_json::perl::Die;
use vep_storable_to_json::regions::{Noted, SequenceRegion};
use vep_storable_to_json::storable::Image;
use vep_storable_to_json::transcripts::{convert, walk};

/// JSON::XS's message for a structure nested past its limit.
const TOO_DEEP: &[u8] =
    b"json text or perl structure exceeds maximum nesting level (max_depth set too low?)";

/// A shard of one transcript whose `biotype` is `value`: the network-order image of
/// `[ { biotype => <value> } ]`.
fn shard(value: &[u8]) -> Image {
    let mut bytes = b"pst0\x05\x0b".to_vec();
    // The list, a reference to the transcript's hash of one entry, the value, the key.
    bytes.extend_from_slice(&[2, 0, 0, 0, 1, 4, 3, 0, 0, 0, 1]);
    bytes.extend_from_slice(value);
    bytes.extend_from_slice(&[0, 0, 0, 7]);
    bytes.extend_from_slice(b"biotype");
    Image::read(bytes).expect("image")
}

/// A reference to arrays nested `depth` deep, each the one item of the array around
/// it and the innermost empty.
fn nested_arrays(depth: usize) -> Vec<u8> {
    let mut value = Vec::with_capacity(depth * 6);
    for _ in 1..depth {
        value.extend_from_slice(&[4, 2, 0, 0, 0, 1]);
    }
    value.extend_from_slice(&[4, 2, 0, 0, 0, 0]);
    value
}

/// A reference to a hash blessed into `Bio::EnsEMBL::Biotype` whose one entry is
/// `name => <name>`.
fn biotype_object(name: &[u8]) -> Vec<u8> {
    let class = b"Bio::EnsEMBL::Biotype";
    let mut value = vec![4, 17, class.len() as u8];
    value.extend_from_slice(class);
    value.extend_from_slice(&[3, 0, 0, 0, 1]);
    value.extend_from_slice(name);
    value.extend_from_slice(&[0, 0, 0, 4]);
    value.extend_from_slice(b"name");
    value
}

/// A reference to a hash blessed into `Bio::EnsEMBL::Slice` whose entries are
/// `seq_region_length => <length>` and `circular => 1`.
fn slice(length: &[u8]) -> Vec<u8> {
    let class = b"Bio::EnsEMBL::Slice";
    let mut value = vec![4, 17, class.len() as u8];
    value.extend_from_slice(class);
    value.extend_from_slice(&[3, 0, 0, 0, 2]);
    value.extend_from_slice(length);
    value.extend_from_slice(&[0, 0, 0, 17]);
    value.extend_from_slice(b"seq_region_length");
    // The integer 1 as one byte.
    value.extend_from_slice(&[8, 129, 0, 0, 0, 8]);
    value.extend_from_slice(b"circular");
    value
}

/// A shard of one transcript: the network-order image of `[ { <entries> } ]`, each entry
/// a key and its value.
fn transcript_shard(entries: &[(&[u8], &[u8])]) -> Image {
    let mut bytes = b"pst0\x05\x0b".to_vec();
    bytes.extend_from_slice(&[2, 0, 0, 0, 1, 4, 3, 0, 0, 0, entries.len() as u8]);
    for (key, value) in entries {
        bytes.extend_from_slice(value);
        bytes.extend_from_slice(&[0, 0, 0, key.len() as u8]);
        bytes.extend_from_slice(key);
    }
    Image::read(bytes).expect("image")
}

/// A shard of one transcript whose slice states the double `length`: the image of
/// `[ { slice => bless({ seq_region_length => <length> }, 'Bio::EnsEMBL::Slice') } ]` in
/// the native order of a little-endian host with 8-byte integers and doubles.
fn native_double_length_shard(length: f64) -> Image {
    let class = b"Bio::EnsEMBL::Slice";
    let mut bytes = b"pst0\x04\x0b\x0812345678\x04\x08\x08\x08".to_vec();
    bytes.extend_from_slice(&[2, 1, 0, 0, 0, 4, 3, 1, 0, 0, 0]);
    bytes.extend_from_slice(&[4, 17, class.len() as u8]);
    bytes.extend_from_slice(class);
    bytes.extend_from_slice(&[3, 1, 0, 0, 0, 7]);
    bytes.extend_from_slice(&length.to_le_bytes());
    bytes.extend_from_slice(&[17, 0, 0, 0]);
    bytes.extend_from_slice(b"seq_region_length");
    bytes.extend_from_slice(&[5, 0, 0, 0]);
    bytes.extend_from_slice(b"slice");
    Image::read(bytes).expect("image")
}

fn encode_error(value: &[u8]) -> Option<Vec<u8>> {
    let converted = convert(&shard(value)).expect("a shard that is walked to its end");
    assert_eq!(converted.transcripts, 1);
    converted.encode_error.map(|error| error.0)
}

#[test]
fn a_copied_field_of_nested_arrays_is_written_as_it_nests() {
    let converted = convert(&shard(&nested_arrays(100))).expect("a shard");
    assert_eq!(converted.transcripts, 1);
    assert!(converted.encode_error.is_none());
    let expected = [
        b"[{\"biotype\":".to_vec(),
        vec![b'['; 100],
        vec![b']'; 100],
        b"}]".to_vec(),
    ]
    .concat();
    assert_eq!(converted.json, expected);
}

#[test]
fn a_copied_field_nested_past_the_encoder_limit_fails_the_shard_at_the_limit() {
    // The deeper image is refused where the limit is met: the walk goes no further
    // down, whatever the depth below.
    for depth in [600, 200_000] {
        assert_eq!(
            encode_error(&nested_arrays(depth)),
            Some(TOO_DEEP.to_vec()),
            "{depth} arrays"
        );
    }
}

#[test]
fn a_copied_field_that_holds_itself_fails_the_shard() {
    // biotype => [ <a reference to that array> ]: the list is tag 0, the reference to
    // the transcript 1, its hash 2, the field's reference 3 and its array 4.
    assert_eq!(
        encode_error(&[4, 2, 0, 0, 0, 1, 4, 0, 0, 0, 0, 4]),
        Some(TOO_DEEP.to_vec())
    );
    // biotype => <a reference to itself>: tag 3 points at tag 3.
    assert_eq!(
        encode_error(&[4, 0, 0, 0, 0, 3]),
        Some(b"cannot encode reference to scalar unless the scalar is 0 or 1".to_vec())
    );
}

#[test]
fn a_biotype_stored_as_an_object_is_written_as_its_name() {
    // name => 'protein_coding', a string of 14 bytes.
    let mut name = vec![10, 14];
    name.extend_from_slice(b"protein_coding");
    let converted = convert(&shard(&biotype_object(&name))).expect("a shard");
    assert_eq!(converted.transcripts, 1);
    assert!(converted.encode_error.is_none());
    assert_eq!(converted.json, b"[{\"biotype\":\"protein_coding\"}]");
}

#[test]
fn a_biotype_object_named_by_a_reference_of_the_class_0_fails_the_shard() {
    // name => bless({}, '0'): `ref` of it is the false string `0`, so Perl takes the
    // name for a plain scalar and writes its address.
    let error = convert(&shard(&biotype_object(&[4, 17, 1, b'0', 3, 0, 0, 0, 0])))
        .err()
        .expect("a shard that fails");
    assert_eq!(
        error.0,
        b"a reference is used as a string: Perl writes its address in memory, which no other run repeats"
    );
}

#[test]
fn a_damaged_length_is_kept_for_the_run_to_judge_and_the_walk_goes_on() {
    // seq_region_length => "5\x80" under the UTF-8 flag: Perl's match on it dies.
    let walked = walk(&transcript_shard(&[(
        b"slice",
        &slice(&[23, 2, b'5', 0x80]),
    )]));
    assert_eq!(
        walked.sequence,
        Noted::Dies(Die::new("Malformed UTF-8 character (fatal)"))
    );
    let shard = walked.shard.expect("a shard that is walked to its end");
    assert_eq!(shard.transcripts, 1);
    assert_eq!(shard.json, b"[{}]");
}

#[test]
fn a_transcript_that_ends_the_walk_has_stated_its_length() {
    // translation => \'x': a reference to a scalar where the conversion reads a hash.
    let walked = walk(&transcript_shard(&[
        (b"slice", &slice(&[8, 128 + 60])),
        (b"translation", &[4, 10, 1, b'x']),
    ]));
    assert_eq!(
        walked.shard.err().expect("a walk that ends"),
        Die::new("Not a HASH reference")
    );
    assert_eq!(
        walked.sequence,
        Noted::Region(SequenceRegion {
            length: b"60".to_vec(),
            utf8: false,
            circular: Some(true),
        })
    );
}

#[test]
fn a_double_is_a_length_when_perl_prints_it_as_digits() {
    assert_eq!(
        walk(&native_double_length_shard(5000.0)).sequence,
        Noted::Region(SequenceRegion {
            length: b"5000".to_vec(),
            utf8: false,
            circular: None,
        })
    );
    // Printed as 1e+15 and as 5000.5.
    for length in [1e15, 5000.5] {
        assert_eq!(
            walk(&native_double_length_shard(length)).sequence,
            Noted::Nothing,
            "{length}"
        );
    }
}
