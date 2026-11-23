// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The plugin comparison of two VEP outputs: the passes that read, sort and pair plugin
//! values across every core, the rules that settle what the passes leave, and the `score`
//! and `compare` programs over both.

/// The shapes of the class `alternate_codon_past_the_coding_end`: a source file of
/// `tools/vep-formats-compare`, which `rules::coding_ends` asks by this name.
#[path = "../../vep-formats-compare/src/rules/coding_end.rs"]
pub mod coding_end;
pub mod compare;
pub mod extract;
/// The genetic codes `coding_end` translates under, a source file of the same tool.
#[path = "../../vep-formats-compare/src/genetic_codes.rs"]
pub mod genetic_codes;
pub mod input;
pub mod order;
pub mod ordered;
pub mod pair;
pub mod pytext;
pub mod registry;
pub mod resolve;
pub mod rows;
pub mod rules;
pub mod score;
pub mod slots;
