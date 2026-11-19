// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--custom` annotation sources: the argument grammar ([`grammar`]), the
//! source abstraction and the rules every format shares ([`source`]), the
//! readers of the tabix-indexed formats ([`vcf`], [`bed`], [`gff_overlay`])
//! and the bigWig reader ([`bigwig`]).
//!
//! A run parses every `--custom` string first ([`parse_all`]), expands a
//! `###CHR###` file over the run's chromosomes ([`expand_all`]), opens the
//! sources ([`CustomSources::open`]) and runs them on each batch after the
//! plugins ([`CustomSources::run_batch`]).

pub mod bed;
pub mod bigwig;
pub mod gff_overlay;
pub mod grammar;
pub mod source;
pub mod vcf;

pub use grammar::{
    default_chromosomes, expand_chr_placeholder, parse_custom, CustomFormat, CustomSpec, GffType,
    GrammarError, MatchType, RecordLimit, SummaryStat, GFF_GTF_RECIPE,
};
pub use source::{CustomError, CustomPrefetch, CustomSource, CustomSources};

/// The run-wide flags a VCF source reads (`File/VCF.pm:114-115`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CustomRunOptions {
    /// `--custom_multi_allelic`: a comma-separated INFO value is never split per ALT.
    pub multi_allelic: bool,
    /// `--custom_suppress_filter`: the FILTER column is left out unless `fields` names it.
    pub suppress_filter: bool,
}

/// Every `--custom` string parsed in the order given; the first refused one
/// ends the parse.
pub fn parse_all(arguments: &[String]) -> Result<Vec<CustomSpec>, GrammarError> {
    arguments
        .iter()
        .map(|argument| parse_custom(argument))
        .collect()
}

/// Every specification with a `###CHR###` file expanded over `chromosomes`
/// ([`expand_chr_placeholder`]), the others as they are, in order.
pub fn expand_all(
    specs: Vec<CustomSpec>,
    chromosomes: &[String],
) -> Result<Vec<CustomSpec>, GrammarError> {
    let mut expanded = Vec::with_capacity(specs.len());
    for spec in specs {
        expanded.extend(expand_chr_placeholder(spec, chromosomes)?);
    }
    Ok(expanded)
}
