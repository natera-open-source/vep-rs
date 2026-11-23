// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The consequence set as the key holds it.

use super::*;

fn normalised(raw: &str) -> String {
    let mut out = Vec::new();
    push_consequence_set(&mut Vec::new(), raw, &mut out);
    String::from_utf8(out).unwrap()
}

#[test]
fn consequence_order_is_normalised_in_the_key() {
    assert_eq!(
        normalised("splice_region_variant,intron_variant"),
        "intron_variant,splice_region_variant"
    );
    assert_eq!(normalised("a,a,b"), "a,b");
}
