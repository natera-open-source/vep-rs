// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! AncestralAllele plugin: the ancestral allele of a variant's span, read from
//! Ensembl's ancestral-sequence FASTA (`AncestralAllele.pm`, VEP_plugins
//! release/116).
//!
//! The plugin takes one parameter, the indexed FASTA (`.fa` with its `.fai`,
//! or the bgzipped `.fa.gz` with its `.fai` and `.gzi`), and writes one field,
//! `AA`, at the record level, so every row of a record carries it, intergenic
//! rows included (`feature_types ['Feature','Intergenic']`, `:105-107`). The
//! value is the FASTA's bases over the record's span as the file spells them:
//! Ensembl's dumps mark a high-confidence call in upper case, a low-confidence
//! call in lower case, a failed call `N`, an insertion in the extant species
//! `-` and no alignment coverage `.`. Two special values stand in for a span
//! the file cannot answer (`:57-60`): `-` for an insertion, whose span is empty
//! (`:124`), and `?` for a chromosome the FASTA has no record for (`:129`).
//! A later positional parameter is not read, as `:90` reads the first alone,
//! and is logged.
//!
//! The FASTA names its records as Ensembl's ancestral dumps do,
//! `ANCESTOR_for_chromosome:GRCh38:21:1:46709983:1`, and the chromosome is the
//! third colon-separated field (`fasta_name_map`, `:273`). A record's
//! chromosome is looked up by name, then under its `chr`-stripped or
//! `chr`-prefixed spelling (`get_sr_name`, `:252-262`). Structural variants
//! get no value: the plugin keeps the default `variant_feature_types` of
//! `['VariationFeature']` (BaseVepPlugin.pm:100), so Ensembl VEP never hands it
//! one.

use std::collections::HashMap;

use indexmap::IndexMap;
use vep_core::codon::reverse_complement;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::coordinate::Strand;
use vep_core::variant::InputVariant;
use vep_fasta::IndexedFasta;

use crate::params::{warn_unread_parameters, PluginParams};
use crate::traits::{BuiltinPlugin, PluginError};

const NAME: &str = "AncestralAllele";

/// The output field (`:110`).
const FIELD: &str = "AA";

/// The value of an insertion (`:124`).
const INSERTION: &str = "-";

/// The value of a chromosome the FASTA has no record for (`:129`).
const UNKNOWN_CHROMOSOME: &str = "?";

/// Ancestral allele plugin.
pub struct AncestralAllelePlugin {
    /// The FASTA path as the parameter gave it, named in the header (`:110`).
    fasta_file: String,
    /// The opened FASTA and its chromosome-to-record map, after `init`.
    source: Option<Source>,
}

/// An opened ancestral FASTA with its chromosome names.
struct Source {
    fasta: IndexedFasta,
    /// Chromosome (the third colon field of a record name) to the record's
    /// name; a chromosome two records claim keeps the later record, as a Perl
    /// hash built from the index's list does (`:273-276`).
    records_by_chromosome: HashMap<String, String>,
}

impl Source {
    fn open(fasta_file: &str) -> Result<Self, PluginError> {
        let fasta = IndexedFasta::from_path(fasta_file)
            .map_err(|error| PluginError::Init(format!("{NAME}: {error}")))?;
        let records_by_chromosome = fasta
            .record_names()
            .filter_map(|name| {
                let chromosome = name.split(':').nth(2)?;
                Some((chromosome.to_string(), name.to_string()))
            })
            .collect();
        Ok(Self {
            fasta,
            records_by_chromosome,
        })
    }

    /// The record `chromosome` reaches (`get_sr_name`, `:244-265`): the
    /// chromosome's own name; a name opening `chr` in any case without those
    /// three characters; any other name with `chr` and with `CHR` in front.
    fn record_for(&self, chromosome: &str) -> Option<&str> {
        let records = &self.records_by_chromosome;
        if let Some(record) = records.get(chromosome) {
            return Some(record);
        }
        let record = match chromosome.get(..3) {
            Some(prefix) if prefix.eq_ignore_ascii_case("chr") => records.get(&chromosome[3..]),
            _ => records
                .get(&format!("chr{chromosome}"))
                .or_else(|| records.get(&format!("CHR{chromosome}"))),
        };
        record.map(String::as_str)
    }

    /// The bases of `record` over `start..=end`, as written, the span cut at
    /// the record's end as htslib's `fai_fetch` cuts a region (`:136-137`), so
    /// a span wholly past the end is empty.
    fn bases(&self, record: &str, start: u64, end: u64) -> Vec<u8> {
        let length = self.fasta.sequence_length(record).unwrap_or(0);
        if start > length {
            return Vec::new();
        }
        self.fasta
            .sequence_as_written(record, start, end.min(length))
            .unwrap_or_default()
    }

    /// The `AA` value of a variant (`run`, `:113-152`).
    fn ancestral_allele(&self, variant: &InputVariant) -> String {
        if variant.start > variant.end {
            return INSERTION.to_string();
        }
        let Some(record) = self.record_for(&variant.original_chr) else {
            return UNKNOWN_CHROMOSOME.to_string();
        };
        let mut bases = self.bases(record, variant.start, variant.end);
        if variant.strand == Strand::Reverse {
            bases = reverse_complement(&bases);
        }
        String::from_utf8_lossy(&bases).into_owned()
    }
}

impl Default for AncestralAllelePlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl AncestralAllelePlugin {
    /// Create an uninitialized plugin.
    pub fn new() -> Self {
        Self {
            fasta_file: String::new(),
            source: None,
        }
    }
}

impl crate::traits::sealed::Sealed for AncestralAllelePlugin {}

impl BuiltinPlugin for AncestralAllelePlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`
    /// (`AncestralAllele.pm:105-106`): every row, the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            FIELD.to_string(),
            format!("Ancestral allele from {}", self.fasta_file),
        )]
    }

    /// The FASTA is the first positional parameter (`:90`) or `file=`; a later
    /// positional parameter is not read, as `:90` reads the first alone. A path
    /// that is not there is refused as Ensembl VEP refuses it (`:91`).
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file"], params)?;
        let (fasta_file, unread) = match (params.get("file"), params.positional()) {
            (Some(file), _) => (file, &[][..]),
            (None, []) => ("", &[][..]),
            (None, [file, unread @ ..]) => (file.as_str(), unread),
        };
        if fasta_file.is_empty() || !std::path::Path::new(fasta_file).exists() {
            return Err(PluginError::Init(format!(
                "{NAME}: Ancestral FASTA file {fasta_file} not provided or found"
            )));
        }
        warn_unread_parameters(NAME, unread, "the ancestral FASTA alone");
        self.source = Some(Source::open(fasta_file)?);
        self.fasta_file = fasta_file.to_string();
        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        let Some(source) = &self.source else {
            return Ok(());
        };
        for variant in variants
            .iter_mut()
            .filter(|variant| !variant.variant_class.is_structural())
        {
            let value = source.ancestral_allele(variant);
            variant.set_plugin_value(FIELD.to_string(), value, false, self.feature_types());
        }
        Ok(())
    }

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Record-level plugin: the value is written in `run_batch`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use vep_core::variant::VariantClass;

    /// A FASTA in Ensembl's ancestral naming: chromosome 21 over two lines
    /// with a low-confidence run, a gap and a no-coverage base, chromosome X
    /// in one line, and a record whose name has no chromosome field.
    fn write_fasta(dir: &Path) -> String {
        let fasta = dir.join("ancestor.fa");
        std::fs::write(
            &fasta,
            b">ANCESTOR_for_chromosome:GRCh38:21:1:12:1\nACGTac\ngt-.NA\n>ANCESTOR_for_chromosome:GRCh38:X:1:4:1\nTTGC\n>scaffold\nAAAA\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("ancestor.fa.fai"),
            b"ANCESTOR_for_chromosome:GRCh38:21:1:12:1\t12\t42\t6\t7\nANCESTOR_for_chromosome:GRCh38:X:1:4:1\t4\t96\t4\t5\nscaffold\t4\t111\t4\t5\n",
        )
        .unwrap();
        fasta.display().to_string()
    }

    fn plugin(fasta: &str) -> AncestralAllelePlugin {
        let mut plugin = AncestralAllelePlugin::new();
        plugin.init(&[fasta.to_string()]).unwrap();
        plugin
    }

    fn variant(chr: &str, start: u64, end: u64, reference: &[u8], alt: &[u8]) -> InputVariant {
        InputVariant::new(
            chr.to_string(),
            start,
            end,
            reference.to_vec(),
            alt.to_vec(),
        )
    }

    fn aa(plugin: &AncestralAllelePlugin, mut variant: InputVariant) -> Option<String> {
        plugin
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();
        variant.plugin_data.get(FIELD).cloned()
    }

    #[test]
    fn reports_the_ancestral_base() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = write_fasta(dir.path());
        let plugin = plugin(&fasta);
        assert_eq!(
            aa(&plugin, variant("21", 1, 1, b"A", b"G")),
            Some("A".into())
        );
        // The case of a low-confidence call is kept.
        assert_eq!(
            aa(&plugin, variant("21", 5, 5, b"A", b"G")),
            Some("a".into())
        );
        // The second line, with its gap and no-coverage characters.
        assert_eq!(
            aa(&plugin, variant("21", 9, 9, b"A", b"G")),
            Some("-".into())
        );
        assert_eq!(
            aa(&plugin, variant("21", 10, 10, b"A", b"G")),
            Some(".".into())
        );
        // The chromosome's `chr` spelling reaches the same record.
        assert_eq!(
            aa(&plugin, variant("chr21", 2, 2, b"C", b"T")),
            Some("C".into())
        );
        assert_eq!(
            aa(&plugin, variant("CHRX", 1, 1, b"T", b"C")),
            Some("T".into())
        );
        assert_eq!(
            plugin.header_info(),
            [(FIELD.to_string(), format!("Ancestral allele from {fasta}"))]
        );
    }

    #[test]
    fn insertion_is_dash() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&write_fasta(dir.path()));
        assert_eq!(
            aa(&plugin, variant("21", 4, 3, b"-", b"TT")),
            Some("-".into())
        );
        // The insertion rule comes before the chromosome lookup.
        assert_eq!(
            aa(&plugin, variant("7", 4, 3, b"-", b"TT")),
            Some("-".into())
        );
    }

    #[test]
    fn unknown_chromosome_is_question_mark() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&write_fasta(dir.path()));
        assert_eq!(
            aa(&plugin, variant("7", 1, 1, b"A", b"G")),
            Some("?".into())
        );
        // A record named without a chromosome field answers no chromosome.
        assert_eq!(
            aa(&plugin, variant("scaffold", 1, 1, b"A", b"G")),
            Some("?".into())
        );
        // A structural variant is not handed to the plugin.
        let mut sv = variant("21", 2, 8, b"N", b"<DEL>");
        sv.variant_class = VariantClass::StructuralDeletion;
        assert_eq!(aa(&plugin, sv), None);
    }

    #[test]
    fn multi_base_span() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&write_fasta(dir.path()));
        // A deletion's span across the line end, as written.
        assert_eq!(
            aa(&plugin, variant("21", 3, 9, b"GTACGGT", b"-")),
            Some("GTacgt-".into())
        );
        // A reverse-strand record reads the reverse complement, case kept.
        let mut reverse = variant("21", 3, 9, b"GTACGGT", b"-");
        reverse.strand = Strand::Reverse;
        assert_eq!(aa(&plugin, reverse), Some("-acgtAC".into()));
        // A span past the record's end is cut there; one wholly past it is empty.
        assert_eq!(
            aa(&plugin, variant("X", 3, 6, b"GCAA", b"-")),
            Some("GC".into())
        );
        assert_eq!(
            aa(&plugin, variant("X", 5, 6, b"AA", b"-")),
            Some(String::new())
        );
    }

    /// A span wholly past the end of an uncompressed FASTA's record is the
    /// empty string, one base past the end or three, as it is from the
    /// bgzipped form; the record's last base reads as itself. Ensembl VEP's
    /// plugin reads a plain `.fa` through `Bio::DB::Fasta`
    /// (AncestralAllele.pm:95, :140), whose `_calc_offset`
    /// (Bio/DB/IndexedBase.pm:779-780) lowers a residue index past the record's
    /// length to the last residue, so such a span reads that one base there.
    #[test]
    fn a_span_past_a_plain_fasta_records_end_is_empty() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&write_fasta(dir.path()));
        // The chromosome 21 record is 12 bases long and ends on `A`.
        assert_eq!(
            aa(&plugin, variant("21", 12, 12, b"A", b"G")),
            Some("A".into())
        );
        assert_eq!(
            aa(&plugin, variant("21", 13, 13, b"C", b"T")),
            Some(String::new())
        );
        assert_eq!(
            aa(&plugin, variant("21", 14, 16, b"TAG", b"-")),
            Some(String::new())
        );
        assert_eq!(
            aa(&plugin, variant("21", 1000, 1000, b"C", b"T")),
            Some(String::new())
        );
        // The X record is 4 bases long and ends on `C`.
        assert_eq!(
            aa(&plugin, variant("X", 5, 5, b"A", b"G")),
            Some(String::new())
        );
    }

    #[test]
    fn the_fasta_is_positional_or_keyed_and_must_exist() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = write_fasta(dir.path());
        let mut keyed = AncestralAllelePlugin::new();
        keyed.init(&[format!("file={fasta}")]).unwrap();
        assert_eq!(
            keyed.header_info()[0].1,
            format!("Ancestral allele from {fasta}")
        );

        let message = |params: &[&str]| {
            let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
            AncestralAllelePlugin::new()
                .init(&params)
                .unwrap_err()
                .to_string()
        };
        assert!(
            message(&[]).contains("AncestralAllele: Ancestral FASTA file  not provided or found")
        );
        assert!(message(&["/no/such.fa"])
            .contains("AncestralAllele: Ancestral FASTA file /no/such.fa not provided or found"));
        assert!(message(&["fasta=x"]).contains("AncestralAllele: unknown parameter 'fasta'"));
        // A token after the FASTA is not read (`AncestralAllele.pm:90` reads the
        // first alone): the plugin loads the FASTA and runs.
        let mut surplus = AncestralAllelePlugin::new();
        surplus.init(&[fasta.clone(), "extra".to_string()]).unwrap();
        assert_eq!(
            surplus.header_info()[0].1,
            format!("Ancestral allele from {fasta}")
        );
        // A FASTA without its index is refused with the index's name.
        let unindexed = dir.path().join("bare.fa");
        std::fs::write(&unindexed, b">x\nACGT\n").unwrap();
        assert!(message(&[&unindexed.display().to_string()]).contains("bare.fa.fai not found"));
    }
}
