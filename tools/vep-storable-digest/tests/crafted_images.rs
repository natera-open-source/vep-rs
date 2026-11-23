// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The dumper and the content copy on images built here: values nested to the depth
//! at which `Data::Dumper` gives up and past it, which no cache holds and a damaged
//! file can. Every walk runs on a thread with the stack of a worker thread.

use vep_storable_digest::digest::{digest_shard, CONTENT_RULE};
use vep_storable_digest::dumper::{dump, Indent, Overlay, Val};
use vep_storable_digest::md5::md5_hex;
use vep_storable_digest::shards::Shard;
use vep_storable_digest::storable::Image;

/// The stack of a thread the standard library starts: what a worker of the run has.
const WORKER_STACK: usize = 2 * 1024 * 1024;

/// `Data::Dumper`'s message at its `Maxrecurse`.
const TOO_DEEP: &str = "Recursion limit of 1000 exceeded";

fn on_a_worker_stack<T: Send + 'static>(work: impl FnOnce() -> T + Send + 'static) -> T {
    std::thread::Builder::new()
        .stack_size(WORKER_STACK)
        .spawn(work)
        .expect("a thread")
        .join()
        .expect("a walk that ends")
}

/// The shard of a network-order image whose values are `body`.
fn shard(body: &[u8]) -> Shard {
    let mut bytes = b"pst0\x05\x0b".to_vec();
    bytes.extend_from_slice(body);
    Shard::new(Image::read(bytes).expect("image"))
}

/// The text of the image's root at `Indent = 0`, as `Dumper(retrieve(...))` gives it.
fn dumped(body: Vec<u8>) -> Result<Vec<u8>, String> {
    on_a_worker_stack(move || {
        let shard = shard(&body);
        let mut text = Vec::new();
        dump(
            &shard,
            &Overlay::default(),
            &Val::RefTo(shard.image.root()),
            Indent::None,
            &mut text,
        )
        .map(|()| text)
        .map_err(|error| error.to_string())
    })
}

// my $top = []; my $r = $top; for (1 .. $below) { my $next = []; push @$r, $next; $r = $next }
// Dumper($top) returns a text for 999 and dies for 1000.
#[test]
fn arrays_nested_to_the_recursion_limit_are_dumped_and_one_deeper_is_refused() {
    let nested = |below: usize| {
        let mut body = Vec::with_capacity(below * 6 + 5);
        for _ in 0..below {
            body.extend_from_slice(&[2, 0, 0, 0, 1, 4]);
        }
        body.extend_from_slice(&[2, 0, 0, 0, 0]);
        body
    };
    assert_eq!(
        dumped(nested(999)),
        Ok([vec![b'['; 1000], vec![b']'; 1000]].concat())
    );
    for below in [1000, 1001, 200_000] {
        assert_eq!(
            dumped(nested(below)),
            Err(TOO_DEEP.to_string()),
            "{below} arrays below the top"
        );
    }
}

// my $top = {}; my $r = $top; for (1 .. $below) { my $next = {}; $r->{k} = $next; $r = $next }
// Dumper($top) returns a text for 999 and dies for 1000.
#[test]
fn hashes_nested_to_the_recursion_limit_are_dumped_and_one_deeper_is_refused() {
    let nested = |below: usize| {
        let mut body = Vec::with_capacity(below * 11 + 5);
        for _ in 0..below {
            body.extend_from_slice(&[3, 0, 0, 0, 1, 4]);
        }
        body.extend_from_slice(&[3, 0, 0, 0, 0]);
        // Every key follows its value: the keys close the hashes innermost first.
        for _ in 0..below {
            body.extend_from_slice(&[0, 0, 0, 1, b'k']);
        }
        body
    };
    let expected = ["{'k' => ".repeat(999), "{}".to_string(), "}".repeat(999)].concat();
    assert_eq!(dumped(nested(999)), Ok(expected.into_bytes()));
    for below in [1000, 200_000] {
        assert_eq!(
            dumped(nested(below)),
            Err(TOO_DEEP.to_string()),
            "{below} hashes below the top"
        );
    }
}

// my $v = 5; my $r = \$v; for (2 .. $references) { my $next = $r; $r = \$next }
// Dumper($r) returns a text for 1000 references and dies for 1001.
#[test]
fn references_nested_to_the_recursion_limit_are_dumped_and_one_deeper_is_refused() {
    // `retrieve` returns one reference more than the image holds.
    let nested = |in_image: usize| {
        let mut body = vec![4u8; in_image];
        body.extend_from_slice(&[8, 133]);
        body
    };
    let mut expected = vec![b'\\'; 1000];
    expected.push(b'5');
    assert_eq!(dumped(nested(999)), Ok(expected));
    for in_image in [1000, 200_000] {
        assert_eq!(
            dumped(nested(in_image)),
            Err(TOO_DEEP.to_string()),
            "{in_image} references in the image"
        );
    }
}

/// The values of a shard of one transcript,
/// `{ adaptor => [ \@l0, \@l1, ... ], b => \@last, stable_id => 'T' }`, where `@l0` is
/// empty and every later array of the `links` holds one reference to the array before
/// it. The whole transcript is dumped from `adaptor`, where every array is met after
/// the one it holds, and so stays shallow; the content copy leaves `adaptor` out and
/// meets the arrays from `b`, each inside the one before.
fn chain_behind_a_key_the_copy_leaves_out(links: usize) -> Vec<u8> {
    // Tags in the order the values are written: the list 0, the reference to the
    // transcript 1, its hash 2, the reference under `adaptor` 3 and its array 4.
    let mut body = vec![2, 0, 0, 0, 1, 4, 3, 0, 0, 0, 3, 4, 2];
    body.extend_from_slice(&(links as u32).to_be_bytes());
    let mut next_tag: u32 = 5;
    let mut last_array: u32 = 0;
    for link in 0..links {
        // A reference (tag `next_tag`) to an array (the tag after it).
        body.push(4);
        let array = next_tag + 1;
        if link == 0 {
            body.extend_from_slice(&[2, 0, 0, 0, 0]);
            next_tag += 2;
        } else {
            // Its one item: a reference of its own to the array before, by its tag.
            body.extend_from_slice(&[2, 0, 0, 0, 1, 4, 0]);
            body.extend_from_slice(&last_array.to_be_bytes());
            next_tag += 3;
        }
        last_array = array;
    }
    body.extend_from_slice(&[0, 0, 0, 7]);
    body.extend_from_slice(b"adaptor");
    body.extend_from_slice(&[4, 0]);
    body.extend_from_slice(&last_array.to_be_bytes());
    body.extend_from_slice(&[0, 0, 0, 1, b'b']);
    body.extend_from_slice(&[10, 1, b'T', 0, 0, 0, 9]);
    body.extend_from_slice(b"stable_id");
    body
}

/// The lines and the failure of the one shard, digested on a worker's stack.
fn digested(links: usize) -> (Vec<Vec<u8>>, Option<Vec<u8>>) {
    on_a_worker_stack(move || {
        let shard = shard(&chain_behind_a_key_the_copy_leaves_out(links));
        let digest = digest_shard(Some(&shard), b"1", CONTENT_RULE);
        (
            digest
                .records
                .into_iter()
                .map(|record| record.line)
                .collect(),
            digest.failure.map(|die| die.0),
        )
    })
}

#[test]
fn a_content_copy_nested_to_the_recursion_limit_is_digested() {
    // The copy is { b => [[ ... [] ... ]], stable_id => 'T' }: 999 arrays below the hash.
    let (lines, failure) = digested(999);
    assert_eq!(failure, None);
    let [line] = lines.as_slice() else {
        panic!("one line, got {}", lines.len());
    };
    let copy = format!(
        "{{'b' => {}{},'stable_id' => 'T'}}",
        "[".repeat(999),
        "]".repeat(999)
    );
    let line = String::from_utf8(line.clone()).expect("ASCII");
    let columns: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
    assert_eq!(columns[..4], ["T", "1", "", ""]);
    assert_eq!(columns[5].as_bytes(), md5_hex(copy.as_bytes()));
    assert_eq!(columns.len(), 6);
}

#[test]
fn a_content_copy_nested_past_the_recursion_limit_ends_the_shard_with_the_dumper_message() {
    // One array more is past the limit; two hundred thousand more are past any stack.
    for links in [1000, 200_000] {
        assert_eq!(
            digested(links),
            (Vec::new(), Some(TOO_DEEP.as_bytes().to_vec())),
            "{links} arrays"
        );
    }
}
