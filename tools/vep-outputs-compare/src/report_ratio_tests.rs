// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Precision, recall and F1 where a side holds no tuple.

use super::*;

#[test]
fn a_ratio_over_nothing_is_zero() {
    assert_eq!(prf(0, 0, 0), (0.0, 0.0, 0.0));
    assert_eq!(prf(0, 0, 3), (0.0, 0.0, 0.0));
    assert_eq!(prf(0, 3, 0), (0.0, 0.0, 0.0));
}
