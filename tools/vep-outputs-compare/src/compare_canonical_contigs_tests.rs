// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The contig set a `--canonical-contigs` value names.

use super::*;

#[test]
fn a_value_names_its_contigs_between_commas_and_blanks_and_one_that_names_none_is_refused() {
    let named: HashSet<String> = ["2L", "X", "mitochondrion_genome"]
        .iter()
        .map(|name| name.to_string())
        .collect();
    assert_eq!(
        parse_canonical_contigs("2L, X mitochondrion_genome"),
        Ok(named)
    );
    let scratch = TempDir::create(&std::env::temp_dir(), "contigs_test_").unwrap();
    let empty = scratch.path.join("empty.txt");
    std::fs::write(&empty, "# nothing\n\n").unwrap();
    let empty = empty.to_str().unwrap();
    assert_eq!(
        parse_canonical_contigs(empty),
        Err(format!(
            "error: --canonical-contigs '{empty}' names no contig"
        ))
    );
    assert_eq!(
        parse_canonical_contigs(","),
        Err("error: --canonical-contigs ',' names no contig".to_string())
    );
}
