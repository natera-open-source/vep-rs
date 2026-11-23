// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two classes of a tabix-backed plugin's field on a record of the mitochondrial contig.
//!
//! BaseVepTabixPlugin.pm:693-708 map a record's chromosome name to a plugin's file by adding
//! or removing `chr` alone, so a record spelt `MT`, `chrM` or `chrMT` finds nothing in a file
//! that names the contig otherwise. dbNSFP.pm:313 names a record on `MT` or `chrMT` `M` to
//! its table first, and :335-338 match a row on its position and alternate allele, so on
//! GRCh37, where `M` and `chrM` name another sequence than `MT` and `chrMT`, the reference
//! writes the value of the row a table naming the contig `M` or `chrM` holds at the record's
//! position, whatever sequence the table counts on.
//!
//! The engine reads a table under all four spellings where they name one sequence, and on
//! GRCh37 under the other pair than the record's only where the table's rows state the bases
//! the run's reference sequence holds at their positions. A rule here reads the two outputs'
//! rows: it tells the record's spelling and which output fills the field, and neither the
//! name the table gives the contig, nor whether its rows state the reference bases, nor
//! whether the run gave the engine a reference sequence.

use super::entry::{blank, field_named, Entry, Row};
use super::shared::record_chromosome;

pub const MITOCHONDRIAL_SPELLINGS: [&str; 3] = ["MT", "chrM", "chrMT"];
/// The record spellings dbNSFP.pm:313 names `M` to the plugin's table.
pub const DBNSFP_NAMED_M_SPELLINGS: [&str; 2] = ["MT", "chrMT"];
/// The registry's name of the dbNSFP plugin.
pub const DBNSFP_PLUGIN: &str = "dbNSFP";
/// The names the engine reads as GRCh37, in lower case.
pub const GRCH37_NAMES: [&str; 3] = ["grch37", "hg19", "37"];
pub const GNOMAD_PLUGIN_PREFIX: &str = "gnomad_";
/// The `gnomAD_` fields are two plugins': gnomADc's end in `_cov`, the others are gnomADMt's.
pub const GNOMADC_FIELD_SUFFIX: &str = "_cov";

/// Whether the field `name` is one the engine alone fills (`ev` blank, `av` not) on a row of
/// a record whose chromosome is one of `MITOCHONDRIAL_SPELLINGS`; a `gnomAD_` field is read
/// only where it is gnomADc's. The rule is the same on every assembly.
pub fn contig_mitochondrial_spelling(
    e: &Entry,
    name: &str,
    ev: &str,
    av: &str,
    record_rows: Option<&[Row]>,
) -> bool {
    if !blank(ev) || blank(av) {
        return false;
    }
    if field_named(name, &[GNOMAD_PLUGIN_PREFIX])
        && !name.to_lowercase().ends_with(GNOMADC_FIELD_SUFFIX)
    {
        return false;
    }
    record_chromosome(e, record_rows)
        .is_some_and(|chromosome| MITOCHONDRIAL_SPELLINGS.contains(&chromosome.as_str()))
}

/// Whether the field is one the reference alone fills (`ev` not blank, `av` blank) on a row
/// of a record whose chromosome is one of `DBNSFP_NAMED_M_SPELLINGS`, in a run whose
/// `assembly` is one of `GRCH37_NAMES`. A run that names no assembly is not known to be on
/// GRCh37, and nothing of it is taken.
pub fn dbnsfp_record_named_m(
    e: &Entry,
    ev: &str,
    av: &str,
    record_rows: Option<&[Row]>,
    assembly: Option<&str>,
) -> bool {
    if blank(ev) || !blank(av) {
        return false;
    }
    let on_grch37 = assembly
        .is_some_and(|assembly| GRCH37_NAMES.contains(&assembly.trim().to_lowercase().as_str()));
    on_grch37
        && record_chromosome(e, record_rows)
            .is_some_and(|chromosome| DBNSFP_NAMED_M_SPELLINGS.contains(&chromosome.as_str()))
}
