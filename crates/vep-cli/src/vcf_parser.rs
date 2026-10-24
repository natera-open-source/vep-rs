// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF line parsing and allele processing.
//!
//! Parses a VCF data line into `InputVariant`s: the line is routed through
//! `vep_core::sv_routing`, a sequence record is split per ALT with its alleles
//! trimmed and minimised here, a `<CNV:TR>` record with the reference sequence to
//! hand is expanded to explicit alleles and read as a sequence record, and any
//! other structural record's alleles are built by
//! `vep_core::sv_record::structural_variants`.
//!
//! Perl citations name modules of ensembl-vep release/116 (`Bio/EnsEMBL/VEP/...`);
//! `Sequence.pm`, `StructuralVariationFeature.pm` and `StructuralVariationOverlap.pm`
//! are in ensembl-variation release/116 and `BaseVCF4.pm` is in ensembl-io release/116.

use std::sync::Arc;

use anyhow::{bail, Context};

use vep_core::coordinate::{parse_pos, MAX_COORD};
use vep_core::minimal::{
    apply_minimal_split, split_changes_record, split_reference_row, trim_sequences, MinimalSplit,
};
use vep_core::sv_record::{structural_variants, tandem_repeat_allele_string, DEFAULT_MAX_SV_SIZE};
use vep_core::sv_routing::{
    joined, perl_split, perl_true, routes_to_structural, sv_so_term_and_class,
    sv_type_abbreviation, sv_type_string, variant_alts, written_alts,
};
use vep_core::sv_span::{structural_end, structural_start, SpanClass, SpanError};
use vep_core::variant::{
    first_vcf_id, id_names_record, so_variation_class, InputVariant, VariantClass,
};
use vep_core::vcf_info::{info_value_written, InfoFields};
use vep_fasta::IndexedFasta;

use crate::individuals::{add_reference_allele_row, RecordAlleles, RecordReading};

/// What a run parses every VCF data line with.
#[derive(Clone, Copy)]
pub struct VcfParse<'a> {
    /// Whether each produced variant keeps the source line in
    /// `InputVariant::raw_input`, one buffer shared by the record's alleles. The
    /// VCF writer writes a VCF input's line with the INFO field appended
    /// (OutputFactory/VCF.pm:315) and, when no line is held or the input is not
    /// a VCF, the line Ensembl VEP builds from the record (:336-337,
    /// `vep_io::output::vcf_output::to_vcf_record`); the JSON writer prints the
    /// line as the `input` key (OutputFactory/JSON.pm:230, 247); the default,
    /// tab and Parquet writers read nothing from it. The Ensembl, HGVS, region,
    /// SPDI and identifier parsers of `vep_io::input` set the field on every
    /// record; this parser alone gates it, on the output format, through
    /// `runner::needs_raw_input_capture`.
    pub capture_raw_input: bool,
    /// The run's `--max_sv_size` (`None` for `-1`), which decides whether a
    /// `<CNV:TR>` record's alternate allele is expanded to an explicit sequence or
    /// kept symbolic, as `_expand_tandem_repeat_allele_string` (`Parser/VCF.pm:403`)
    /// decides under `max_sv_size`.
    pub max_sv_size: Option<u64>,
    /// `--minimal`: every bi-allelic sequence record is trimmed and marked
    /// (`vep_core::minimal`), and a multi-allelic record one of whose alleles trims
    /// keeps every ALT, a repeated one and one equal to the reference included, for
    /// the split `split_minimal_record` applies once the record is validated.
    pub minimal: bool,
    /// The run's reference FASTA, which a `<CNV:TR>` record's expansion reads the
    /// run's reference bases from (`Parser/VCF.pm:516-518`); without one the record
    /// keeps its symbolic reading.
    pub fasta: Option<&'a IndexedFasta>,
    /// The run's per-individual reading and `--process_ref_homs`
    /// ([`RecordReading`]): a sequence record's alleles are copied per individual,
    /// or gain the reference allele's row; a structural record is returned as it
    /// is, since Ensembl VEP builds it before either step (Parser/VCF.pm:237-245
    /// returns ahead of :355), and a `<CNV:TR>` expansion, built there too, is read
    /// once and gains the row alone.
    pub reading: RecordReading<'a>,
}

impl VcfParse<'_> {
    /// Ensembl VEP's default `--max_sv_size`, no `--minimal` and no FASTA.
    pub fn default_run(capture_raw_input: bool) -> Self {
        VcfParse {
            capture_raw_input,
            max_sv_size: Some(DEFAULT_MAX_SV_SIZE),
            minimal: false,
            fasta: None,
            reading: RecordReading::default(),
        }
    }
}

/// Parse a VCF data line into one or more `InputVariant`s.
///
/// Multi-allelic sites are split into separate records, one per ALT allele.
/// This is the `vep` binary's one VCF reader: `vep-io` declares no VCF input
/// module. `capture_raw_input` is [`VcfParse::capture_raw_input`]; the line is read
/// under Ensembl VEP's default `--max_sv_size`, without `--minimal` and without a
/// FASTA ([`parse_vcf_line_with`] takes the run's own settings).
pub fn parse_vcf_line(line: &str, capture_raw_input: bool) -> anyhow::Result<Vec<InputVariant>> {
    parse_vcf_line_with(line, &VcfParse::default_run(capture_raw_input))
}

/// [`parse_vcf_line`] under the run's `--max_sv_size` (`None` for `-1`).
pub fn parse_vcf_line_with_max_sv_size(
    line: &str,
    capture_raw_input: bool,
    max_sv_size: Option<u64>,
) -> anyhow::Result<Vec<InputVariant>> {
    parse_vcf_line_with(
        line,
        &VcfParse {
            max_sv_size,
            ..VcfParse::default_run(capture_raw_input)
        },
    )
}

/// [`parse_vcf_line`] under the run's settings.
pub fn parse_vcf_line_with(line: &str, parse: &VcfParse<'_>) -> anyhow::Result<Vec<InputVariant>> {
    let columns = vcf_columns(line)?;
    parse_columns(&columns, &columns.info(), line, parse)
}

/// What the reader gets from one data line: the record's alleles, or the parse
/// error Ensembl VEP reports as a refused record, the warnings Ensembl VEP's
/// tandem-repeat branch writes as it reads the record (Parser/VCF.pm:502-527),
/// ahead of `validate_vf` and of any refusal, so the reader reports them before
/// the alleles' skips, and the record's INFO column for the lookups the skip
/// rules make on it after the parse (`None` for a line of fewer than five
/// columns, which has no record).
pub(crate) struct ParsedRecord<'a> {
    pub alleles: anyhow::Result<Vec<InputVariant>>,
    pub warnings: Vec<String>,
    pub info: Option<RecordInfo<'a>>,
}

/// One parsed record's INFO column, read after the parse by the skip rules.
pub(crate) enum RecordInfo<'a> {
    /// The view the parse split, borrowing the line.
    Split(InfoFields<'a>),
    /// The column copied out of a line the parser read through a copy with its
    /// invalid UTF-8 replaced ([`ParsedRecord::with_info_copied`]).
    Copied(String),
}

impl RecordInfo<'_> {
    /// The value of `key` when the column writes it `key=value`
    /// ([`InfoFields::value_written`]), a copied column scanned for the one key.
    pub(crate) fn value_written(&self, key: &str) -> Option<&str> {
        match self {
            RecordInfo::Split(fields) => fields.value_written(key),
            RecordInfo::Copied(text) => info_value_written(text, key),
        }
    }
}

impl ParsedRecord<'_> {
    /// The record with its INFO column copied out of `line`, the text it was
    /// parsed from, in place of the view over that text: for a reader whose
    /// text ends before the skip rules read the record.
    pub(crate) fn with_info_copied<'b>(self, line: &str) -> ParsedRecord<'b> {
        ParsedRecord {
            alleles: self.alleles,
            warnings: self.warnings,
            info: vcf_columns(line)
                .ok()
                .map(|columns| RecordInfo::Copied(columns.fixed[7].to_string())),
        }
    }
}

/// One data line read once: its columns split once, its INFO column split by the
/// first lookup and read by every other, the tandem-repeat warnings, then the
/// record as [`non_variant_record`] reads a line whose first ALT is `.`
/// ([`vcf_line_is_non_variant`]) and as [`parse_vcf_line_with`] reads any other,
/// the INFO view carried with the record. A line of fewer than five columns is
/// the parse error alone.
pub(crate) fn parse_data_line<'a>(line: &'a str, parse: &VcfParse<'_>) -> ParsedRecord<'a> {
    let columns = match vcf_columns(line) {
        Ok(columns) => columns,
        Err(error) => {
            return ParsedRecord {
                alleles: Err(error),
                warnings: Vec::new(),
                info: None,
            }
        }
    };
    let info = columns.info();
    let warnings = tandem_repeat_warnings(&columns, &info, parse.fasta.is_some());
    let alleles = if columns_are_non_variant(&columns, &info) {
        non_variant_columns(&columns, &info, line, parse)
    } else {
        parse_columns(&columns, &info, line, parse)
    };
    ParsedRecord {
        alleles,
        warnings,
        info: Some(RecordInfo::Split(info)),
    }
}

/// [`parse_vcf_line_with`] over the line's columns and INFO view.
fn parse_columns(
    columns: &VcfColumns<'_>,
    info: &InfoFields<'_>,
    line: &str,
    parse: &VcfParse<'_>,
) -> anyhow::Result<Vec<InputVariant>> {
    let fields = columns.fixed;

    let chr = fields[0].to_string();
    let pos = parse_pos(fields[1]).with_context(|| format!("Invalid POS: '{}'", fields[1]))?;
    let id = fields[2];
    let ref_allele = fields[3];
    let alt_field = fields[4];

    if alt_field == "." {
        bail!("No alternate allele");
    }

    let svtype = info.value("SVTYPE").filter(|v| perl_true(v));

    // The ALT column as Perl's list-context `split` reads it, a gVCF variant site read
    // as its real alleles (`variant_alts`), each ALT keeping its written position.
    let written = written_alts(alt_field);
    let mut alts = written.clone();
    variant_alts(&mut alts, svtype);
    // `*` and `.` name no allele of their own here, where the record's shape is
    // decided; both stay in the record's allele string, which Ensembl VEP builds from
    // every ALT (Parser/VCF.pm:469). On the sequence path Ensembl VEP numbers a `*`
    // and writes no row for it (VariationFeatureOverlap.pm:483-484,
    // `perl_allele_numbers`) and reads a `.` that is not the first ALT as an allele of
    // its own (a list opening with `.` is a non-variant record, `:259`); on the
    // structural path it writes a row per `/`-part of the string, `*` and `.` included
    // (StructuralVariationOverlap.pm:64-77, `structural_variants`). `<*>` is not
    // filtered: Perl treats it as a non-SV allele at POS unless INFO/SVTYPE routes the
    // record. `<CPX>` is not filtered: absent from `%SO_TERMS`, it is parsed with the
    // `vep_skip` mark (`sv_so_term_and_class`). `<DEL:*>` is routed first: under a true
    // INFO/SVTYPE it is a structural allele of the record's class, and on the sequence
    // path an allele Ensembl VEP numbers and writes no row for (`perl_allele_numbers`).
    let valid_alts: Vec<(usize, &str)> = alts
        .iter()
        .copied()
        .filter(|(_, alt)| *alt != "*" && *alt != ".")
        .collect();

    let structural = routes_to_structural(ref_allele, &alts, svtype);
    let raw_input: Option<Arc<str>> = parse
        .capture_raw_input
        .then(|| stored_line(columns, line, structural));
    let record_id = first_vcf_id(id).map(str::to_string);
    let explicit_id = record_id.as_deref().filter(|id| id_names_record(id));
    let line_name: Arc<str> = input_line_name(&chr, pos, ref_allele, alt_field).into();

    // A sequence record carries its alleles as written and as chopped for the
    // per-individual reading; a structural record carries none, and neither does
    // a `<CNV:TR>` expansion, which Ensembl VEP builds in its structural method.
    let (mut variants, alleles) = if valid_alts.is_empty() {
        if written.is_empty() || written.iter().any(|(_, alt)| *alt != "*") {
            bail!("No alternate allele");
        }
        let chopped = chop_record(ref_allele, &written, pos);
        let alleles = record_alleles(ref_allele, &written, &written, &chopped);
        (
            with_allele_less_fallback(&chr, chopped, Vec::new(), parse.minimal),
            Some(alleles),
        )
    } else if !structural {
        let chopped = chop_record(ref_allele, &alts, pos);
        let variants = sequence_variants(&chr, &chopped, ref_allele, &alts, parse);
        let alleles = record_alleles(ref_allele, &written, &alts, &chopped);
        (
            with_allele_less_fallback(&chr, chopped, variants, parse.minimal),
            Some(alleles),
        )
    } else if let Some(expanded) =
        tandem_repeat_record(&chr, pos, ref_allele, &alts, svtype, info, parse)
    {
        // The expansion's own strings are the record's alleles as given.
        let expanded_alts: Vec<(usize, &str)> = expanded
            .written_case_alts
            .iter()
            .map(String::as_str)
            .enumerate()
            .collect();
        let variants = sequence_variants(
            &chr,
            &expanded,
            &expanded.written_case_reference,
            &expanded_alts,
            parse,
        );
        let mut expansion = with_allele_less_fallback(&chr, expanded, variants, parse.minimal);
        // The structural method's return (Parser/VCF.pm:543, through :244) comes
        // ahead of the per-individual step (:355), so the expansion is read once
        // whatever `--individual` names, while the reference row of
        // `--process_ref_homs` is added at output to every variation feature
        // (OutputFactory.pm:509-516), the expansion included.
        if parse.reading.process_ref_homs {
            add_reference_allele_row(&mut expansion);
        }
        (expansion, None)
    } else {
        (
            structural_variants(&chr, pos, ref_allele, &written, info, explicit_id)?,
            None,
        )
    };

    for variant in &mut variants {
        variant.raw_input = raw_input.clone();
        variant.id = record_id.clone();
        variant.input_line_name = Some(line_name.clone());
    }

    Ok(match alleles {
        Some(alleles) => parse
            .reading
            .finish_record(variants, &alleles, columns, parse.minimal),
        None => variants,
    })
}

/// The record of a VCF line whose first ALT is `.` under `--allow_non_variant`
/// (Parser/VCF.pm:259-267, 342): the reference bases from POS with the reference
/// allele as its whole allele string and no allele to annotate, so no format
/// writes a row for it, the VCF output carries its line and the JSON output its
/// object with `most_severe_consequence` `?`, as Ensembl VEP writes them. A
/// Perl-true INFO/SVTYPE is kept as the record's type, the mark of the line
/// Ensembl VEP routes to its structural path (`:237-245`), whose stored line is
/// filled as that path fills it ([`stored_line`]). Under `--individual`
/// the record is copied per individual like any other (`:355-376`): a `0/0`
/// genotype is a reference-homozygous copy, and the alleles a genotype indexes
/// are the reference and every ALT as written, none chopped (`:300-336` never
/// chops an allele list whose first allele is `.`); a copy whose genotype names
/// an ALT beyond `.` carries that allele (`:627-639` builds the copy's allele
/// string from the alleles named, and nothing reads the record's non-variant
/// mark on it), so such a record is read as a sequence record for its copies.
pub fn non_variant_record(line: &str, parse: &VcfParse<'_>) -> anyhow::Result<Vec<InputVariant>> {
    let columns = vcf_columns(line)?;
    non_variant_columns(&columns, &columns.info(), line, parse)
}

/// [`non_variant_record`] over the line's columns and INFO view.
fn non_variant_columns(
    columns: &VcfColumns<'_>,
    info: &InfoFields<'_>,
    line: &str,
    parse: &VcfParse<'_>,
) -> anyhow::Result<Vec<InputVariant>> {
    let [chr, pos_text, id, reference, alt_field, _, _, _] = columns.fixed;
    let svtype = info.value("SVTYPE").filter(|v| perl_true(v));
    let mut alts = written_alts(alt_field);
    variant_alts(&mut alts, svtype);
    if parse.reading.reads_copies() && alts.iter().any(|(_, alt)| *alt != "." && *alt != "*") {
        return parse_columns(columns, info, line, parse);
    }
    let pos = parse_pos(pos_text).with_context(|| format!("Invalid POS: '{pos_text}'"))?;
    let mut variant =
        InputVariant::without_allele(chr.to_string(), pos, reference, std::iter::empty());
    variant.uploaded_allele_string = Some(variant.allele_string.as_str().into());
    variant.id = first_vcf_id(id).map(str::to_string);
    variant.input_line_name = Some(input_line_name(chr, pos, reference, alt_field).into());
    variant.sv_type = svtype.map(str::to_string);
    variant.raw_input = parse.capture_raw_input.then(|| {
        stored_line(
            columns,
            line,
            routes_to_structural(reference, &alts, svtype),
        )
    });
    let written: Vec<&str> = std::iter::once(reference)
        .chain(perl_split(alt_field, &[',']))
        .collect();
    let alleles = RecordAlleles {
        chopped: written.iter().map(|allele| allele.to_string()).collect(),
        written,
    };
    Ok(parse
        .reading
        .finish_record(vec![variant], &alleles, columns, parse.minimal))
}

/// [`RecordAlleles`] of a sequence record: the reference and every written ALT,
/// each as the chopped allele of `chopped`, found by its written position among
/// `alts`, the ALTs `chop_record` ran over. An ALT outside `alts` (a gVCF
/// `<NON_REF>` beside real alleles, which `variant_alts` leaves out) was never
/// chopped and stays as written.
fn record_alleles<'a>(
    ref_allele: &'a str,
    written: &[(usize, &'a str)],
    alts: &[(usize, &str)],
    chopped: &Chopped,
) -> RecordAlleles<'a> {
    let chopped_alt = |written_index: usize, alt: &str| -> String {
        match alts.iter().position(|&(index, _)| index == written_index) {
            Some(chopped_index) => chopped.written_case_alts[chopped_index].clone(),
            None => alt.to_string(),
        }
    };
    RecordAlleles {
        written: std::iter::once(ref_allele)
            .chain(written.iter().map(|&(_, alt)| alt))
            .collect(),
        chopped: std::iter::once(chopped.written_case_reference.clone())
            .chain(written.iter().map(|&(index, alt)| chopped_alt(index, alt)))
            .collect(),
    }
}

/// `variants`, or, when a sequence record's ALTs write no row, the record as the
/// allele-less feature over its chopped alleles ([`allele_less_record`]). The
/// feature is minimised as `minimise_alleles` (Parser.pm:987-1011) minimises it
/// before the alleles are read: under `--minimal` every bi-allelic record, marked
/// whatever its allele writes, so `A A` at POS is the feature `-/-` at POS + 1 to
/// POS, its class `sequence_alteration`, and `A *` keeps its alleles, both marked
/// minimised; whatever the flags a bi-allelic record whose alleles differ in
/// length (`post_process_vfs`, :876-881), so `C *C` at POS is the feature `-/*`
/// at POS to POS - 1, marked minimised.
fn with_allele_less_fallback(
    chr: &str,
    chopped: Chopped,
    variants: Vec<InputVariant>,
    minimal: bool,
) -> Vec<InputVariant> {
    if !variants.is_empty() {
        return variants;
    }
    let bi_allelic = chopped.written_case_alts.len() == 1;
    let minimises = bi_allelic
        && (minimal
            || perl_minimises(
                &chopped.written_case_reference,
                &chopped.written_case_alts[0],
            ));
    let trimmed = minimises.then(|| {
        trim_sequences(
            chopped.written_case_reference.as_bytes(),
            chopped.written_case_alts[0].as_bytes(),
            chopped.start,
            chopped.end,
        )
    });
    let mut record = allele_less_record(chr.to_string(), chopped);
    record.minimised = bi_allelic && minimal;
    if let Some((reference, alt, start, end, true)) = trimmed {
        record.minimised = true;
        record.written_case_allele_string = written_case_pair(&reference, &alt);
        let reference = reference.to_ascii_uppercase();
        let alt = alt.to_ascii_uppercase();
        record.allele_string = format!(
            "{}/{}",
            String::from_utf8_lossy(&reference),
            String::from_utf8_lossy(&alt)
        );
        record.record_so_term = Some(so_variation_class(&record.allele_string).to_string());
        record.ref_allele = reference;
        record.start = start;
        record.end = end;
    }
    vec![record]
}

/// The variants of a sequence record, one per ALT that writes a row, over its
/// chopped alleles: `raw_reference` and `alts` are the alleles as the record gave
/// them (the VCF columns, or the expansion of a `<CNV:TR>` record), each ALT with
/// its written position, which the `.` test and `post_process_vfs`'s length test
/// read; the `*` test reads the chopped allele.
fn sequence_variants(
    chr: &str,
    chopped: &Chopped,
    raw_reference: &str,
    alts: &[(usize, &str)],
    parse: &VcfParse<'_>,
) -> Vec<InputVariant> {
    let Chopped {
        reference: trimmed_ref,
        alts: trimmed_alts,
        uploaded: uploaded_alleles,
        working: record_string,
        written_case_reference,
        written_case_alts,
        ..
    } = chopped;
    let (vep_start, vep_end) = (chopped.start, chopped.end);
    let bi_allelic = alts.len() == 1;
    // Under `--minimal` a multi-allelic record splits when any of its chopped alleles
    // trims (InputBuffer.pm:665-666, :700); every ALT then takes its written position
    // as its number and a repeated ALT writes its own rows (:684), so the overlap
    // rules of `perl_allele_numbers` give way to the split's, which
    // `split_minimal_record` applies once `validate_vf` has read the record.
    let split_changes = parse.minimal
        && !bi_allelic
        && split_changes_record(trimmed_ref, trimmed_alts.iter().map(String::as_str));
    let numbers = perl_allele_numbers(trimmed_ref, trimmed_alts);
    let written_case_record = chopped.written_case_string();
    let mut variants = Vec::with_capacity(alts.len());
    for (k, (((allele_idx, alt), vep_alt), number)) in
        alts.iter().zip(trimmed_alts).zip(&numbers).enumerate()
    {
        // The `*` test reads the chopped allele, not the written one: the one ALT
        // `C*` of a `C` record is the allele `*` VariationFeatureOverlap.pm:484 skips.
        if vep_alt == "*" || *alt == "." || vep_alt == "<DEL:*>" {
            continue;
        }
        // A dropped allele, or one Ensembl VEP numbers and writes no row for; under
        // the split it writes its rows under its written position. A gVCF block
        // `variant_alts` left out keeps its written position's number, a record
        // Ensembl VEP drops whole (`sv_routing::variant_alts`).
        let number = match *number {
            Some(number) => Some(number + (*allele_idx - k) as u32),
            None if split_changes => None,
            None => continue,
        };
        let mut variant = InputVariant::new(
            chr.to_string(),
            vep_start,
            vep_end,
            trimmed_ref.as_bytes().to_vec(),
            vep_alt.as_bytes().to_vec(),
        );

        // `trim_sequences` (ensembl-variation Sequence.pm:971) reads an end of 0
        // as no end and resets it to `start + length(ref) - 1`, so the one
        // record it reaches with that end, an insertion at POS 0 whose REF is one
        // base (`0 A AT`), leaves the parser as `-/T` at 1-1 and `validate_vf`
        // refuses it as an insertion whose coordinates are not start = end + 1.
        // An empty REF never reaches it: unchopped, its string `/A` fails the
        // `.+/.+` test of `minimise_alleles` (Parser.pm:974), so `1 <empty> A`
        // keeps the end of 0.
        if bi_allelic
            && !raw_reference.is_empty()
            && perl_minimises(raw_reference, alt)
            && variant.end == 0
        {
            variant.end = variant.start + variant.ref_allele.len() as u64 - 1;
        }

        // Secondary minimisation of a bi-allelic record. Ensembl VEP
        // (Parser.pm `post_process_vfs`) minimises a bi-allelic record whose
        // alleles differ in length whatever the flags, and under `--minimal`
        // every bi-allelic record (`minimise_alleles`, :852, :995): `trim_sequences`
        // (Sequence.pm) trims from the left while both alleles are non-empty,
        // then from the right, and maps an emptied allele to `-`, so `T/TT`
        // at pos P becomes `-/T` at P+1 and `CAT/CCT` under the flag `A/C` at P+1.
        // It compares the bases as written (`eq`), before `validate_vf` upper-cases
        // the string (Parser.pm:554), so `TTGA t` stays `TTGA/T`. `minimise_alleles`
        // returns a multi-allelic string (`.+/.+/.+`) unchanged, a `*` part counted
        // (Parser.pm:978-984): without `--minimal` a multi-allelic record keeps its
        // anchor-chopped alleles, and every position, codon and distance is
        // computed on that representation.
        let mut written_case = written_case_record.clone();
        if bi_allelic
            && variant.ref_allele != b"-"
            && variant.alt_alleles[0] != b"-"
            && (parse.minimal || variant.ref_allele.len() != variant.alt_alleles[0].len())
        {
            // The end is the record's own, not derived from the trimmed reference
            // (Sequence.pm:1004-1031): a `<CNV:TR>` expansion over the REF column spans
            // the run, which its reference does not measure, so the span `validate_vf`
            // refuses stays refused.
            let (min_ref, min_alt, new_start, new_end, changed) = trim_sequences(
                written_case_reference.as_bytes(),
                written_case_alts[k].as_bytes(),
                variant.start,
                variant.end,
            );
            if changed {
                written_case = written_case_pair(&min_ref, &min_alt);
                let min_ref = min_ref.to_ascii_uppercase();
                let min_alt = min_alt.to_ascii_uppercase();
                variant.original_allele_string = Some(variant.allele_string.clone());
                variant.allele_string = format!(
                    "{}/{}",
                    String::from_utf8_lossy(&min_ref),
                    String::from_utf8_lossy(&min_alt),
                );
                variant.ref_allele = min_ref;
                variant.alt_alleles = vec![min_alt];
                variant.start = new_start;
                variant.end = new_end;
                variant.variant_class = vep_core::variant::classify_variant(
                    &variant.ref_allele,
                    &variant.alt_alleles[0],
                );
                variant.minimised = true;
            }
        }
        variant.written_case_allele_string = written_case;
        // `minimise_alleles` marks every bi-allelic record it handles (Parser.pm:1011),
        // trimmed or not.
        if bi_allelic && parse.minimal {
            variant.minimised = true;
        }
        // The overlap alleles are built from the minimised string
        // (VariationFeatureOverlap.pm:484 after Parser.pm:881), so a `*` the
        // minimisation leaves (`C *C`, or `CA C*` under `--minimal`) writes no row
        // either: the record is the allele-less feature, minimised the same way.
        if variant.alt_alleles[0] == b"*" {
            continue;
        }

        variant.allele_index = *allele_idx;
        // ALLELE_NUM is Ensembl VEP's count of the alleles it kept, which falls
        // behind the written position once an ALT was dropped.
        variant.allele_num = number.filter(|number| *number as usize != *allele_idx + 1);
        variant.uploaded_allele_string = Some(uploaded_alleles.clone());
        variant.record_allele_string_multi = record_string.clone();
        variants.push(variant);
    }
    variants
}

/// Whether `split_variants` splits a validated record (InputBuffer.pm:648,
/// :665-666, :704): a sequence record whose working string names two or more
/// alleles beyond the reference (`.+/.+/.+`, an empty trailing allele counting
/// as none) and one of whose alleles trims against its reference, the reference
/// as it stands after `validate_vf`. An individual's copy carrying one alternate
/// allele of a multi-allelic record is bi-allelic and is neither split nor
/// minimised: `minimise_alleles` runs on the record before the copy is made
/// (Parser/VCF.pm:362-364), so the copy keeps the record's chopped alleles.
pub(crate) fn record_splits(variant: &InputVariant) -> bool {
    if variant.is_structural || variant.no_allele {
        return false;
    }
    variant
        .record_allele_string_multi
        .as_deref()
        .is_some_and(|working| {
            let mut parts = working.split('/').filter(|part| !part.is_empty());
            let reference = parts.next().unwrap_or("");
            let alts: Vec<&str> = parts.collect();
            alts.len() >= 2 && split_changes_record(reference, alts)
        })
}

/// `split_variants` (InputBuffer.pm:638-718) over the alleles of one validated
/// record, in their written order. A record [`record_splits`] names takes each
/// allele on its own trimmed span (`vep_core::minimal::apply_minimal_split`), a
/// repeated ALT with its own rows (:658) and an allele equal to the reference
/// once trimmed with none (VariationFeatureOverlap.pm:481); a record left with no
/// allele is the feature with no allele to annotate, over the record's span and
/// working string, as [`allele_less_record`] builds it. The alleles are, in order,
/// the working string's alternate alleles that write a row, so each finds its
/// place in that string by its chopped allele. A reference allele's row leading
/// them (`--process_ref_homs`) gives way to one per allele, each its feature's
/// (`vep_core::minimal::split_reference_row`), before the allele's own rows, the
/// allele equal to the reference included, and one for the feature of a `*` or
/// `<DEL:*>` allele, which the alleles leave out. Any other record is returned as
/// it stands. Runs after `validate_vf` read the record, so a reference
/// `--lookup_ref` replaced is the one every allele trims against.
pub(crate) fn split_minimal_record(mut alleles: Vec<InputVariant>) -> Vec<InputVariant> {
    let reference_rows = alleles.first().is_some_and(|v| v.reference_allele_row);
    let Some(first) = alleles.get(usize::from(reference_rows)) else {
        return alleles;
    };
    if !record_splits(first) {
        return alleles;
    }
    let Some(working) = first.record_allele_string_multi.clone() else {
        return alleles;
    };
    let record_row = reference_rows.then(|| alleles.remove(0));
    let working_alts: Vec<&str> = working.split('/').skip(1).collect();
    let mut kept = Vec::with_capacity(alleles.len() * 2);
    let mut first_without_allele = None;
    let mut cursor = 0;
    // A `*` or `<DEL:*>` allele writes no row and the parser carries none for it,
    // but its feature, which no trim changes, builds the reference allele before
    // the allele is passed over (VariationFeatureOverlap.pm:469-474, :484), so
    // under the reference rows it writes that row alone; the record's row is the
    // template, carrying the record's span and fields.
    let rowless_features = |kept: &mut Vec<InputVariant>, from: usize, to: usize| {
        let Some(row) = &record_row else {
            return;
        };
        for (own, alt) in working_alts.iter().enumerate().take(to).skip(from) {
            if !matches!(*alt, "*" | "<DEL:*>") {
                continue;
            }
            let mut feature = row.clone();
            feature.alt_alleles = vec![alt.as_bytes().to_vec()];
            apply_minimal_split(&mut feature, own);
            kept.push(split_reference_row(&feature, own));
        }
    };
    for mut allele in alleles {
        let own = working_alts[cursor.min(working_alts.len())..]
            .iter()
            .position(|alt| alt.as_bytes() == allele.alt_allele())
            .map_or(cursor, |found| cursor + found);
        rowless_features(&mut kept, cursor, own);
        cursor = own + 1;
        let outcome = apply_minimal_split(&mut allele, own);
        if reference_rows {
            kept.push(split_reference_row(&allele, own));
        }
        if outcome == MinimalSplit::NoAllele {
            first_without_allele.get_or_insert(allele);
        } else {
            kept.push(allele);
        }
    }
    rowless_features(&mut kept, cursor, working_alts.len());
    if kept.is_empty() {
        if let Some(mut record) = first_without_allele {
            record.no_allele = true;
            record.alt_alleles = vec![b"*".to_vec()];
            record.allele_string = working.to_string();
            kept.push(record);
        }
    }
    kept
}

/// Whether a VCF record is one Ensembl VEP's tandem-repeat branch reads
/// (`$so_term =~ /tandem/`, Parser/VCF.pm:502): its type abbreviation, from the
/// ALTs joined by `/` with `svtype`, the Perl-true INFO/SVTYPE, standing in for a
/// string outside the VCF 4.4 forms, names a tandem repeat or a tandem duplication.
fn is_tandem_record(alts: &[(usize, &str)], svtype: Option<&str>) -> bool {
    matches!(
        sv_so_term_and_class(&sv_type_abbreviation(&joined(alts, "/"), svtype)),
        Some((
            _,
            VariantClass::TandemRepeat | VariantClass::TandemDuplication
        ))
    )
}

/// The warnings Ensembl VEP's tandem-repeat branch writes for a record
/// (Parser/VCF.pm:502-527, 382-388), in its order, for the columns of a line the
/// reader parsed and their INFO view: that the largest of several SVLEN values
/// is used (`:507-510`, the values compared as written); without a FASTA, once per
/// run through the sink's deduplication, that the calculation is less precise
/// (`:521-524`); and with one, when the record carries RUS and RUC or RB, that
/// CIRUC and CIRB are ignored (`:387-388`). Empty for any other record: one the
/// structural routing (`:237-245`) leaves on the sequence path, INFO/SVTYPE included,
/// or one of another structural class.
fn tandem_repeat_warnings(
    columns: &VcfColumns<'_>,
    info: &InfoFields<'_>,
    fasta_present: bool,
) -> Vec<String> {
    let (ref_allele, alt_field) = (columns.fixed[3], columns.fixed[4]);
    let svtype = info.value("SVTYPE").filter(|v| perl_true(v));
    let mut alts = written_alts(alt_field);
    variant_alts(&mut alts, svtype);
    if !routes_to_structural(ref_allele, &alts, svtype) || !is_tandem_record(&alts, svtype) {
        return Vec::new();
    }
    let mut warnings = Vec::new();
    if let Some(svlen) = info.value_written("SVLEN") {
        // `uniq(split(/,/, $info->{SVLEN}))` (`:505`): a trailing empty value is no
        // value, so `SVLEN=9,` is one reference.
        let mut values = perl_split(svlen, &[',']);
        let first = values.next().unwrap_or("");
        if values.any(|value| value != first) {
            warnings.push(format!(
                "found tandem repeats with different references per alternative allele: SVLEN={svlen}; only using reference with largest size"
            ));
        }
    }
    let has_units = |key: &str| info.value(key).is_some_and(perl_true);
    if !fasta_present {
        warnings.push(
            "could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta"
                .to_string(),
        );
    } else if has_units("RUS")
        && (has_units("RUC") || has_units("RB"))
        && (has_units("CIRUC") || has_units("CIRB"))
    {
        warnings.push(
            "CIRUC and CIRB INFO fields are ignored when calculating alternative alleles in tandem repeats"
                .to_string(),
        );
    }
    warnings
}

/// A `<CNV:TR>` record expanded to explicit alleles as `create_StructuralVariationFeatures`
/// expands it with a reference sequence to hand (Parser/VCF.pm:502-545): the run is
/// the bases from POS + 1 to the end the largest SVLEN or END gives (`:505-511`),
/// its reference read from the FASTA (`:516-518`; the record's REF column when the
/// FASTA has no such span, as `sub_Slice` then yields nothing and `$ref` keeps the
/// column) and each alternate allele the repeat sequence `tandem_repeat_allele_string`
/// builds from INFO. The record is then a sequence record (`:529-543`), an
/// insertion or deletion of repeat units once `post_process_vfs` minimises it.
/// `None`, the symbolic reading, without a FASTA, for a record that is no tandem
/// repeat or tandem duplication (`$so_term =~ /tandem/`, `:502`), for one the
/// expansion leaves symbolic (no RUS, no RUC or RB, a sequence over
/// `--max_sv_size`) and for a span Ensembl VEP refuses. The branch's warnings
/// are [`tandem_repeat_warnings`]'s, written by the reader.
fn tandem_repeat_record(
    chr: &str,
    pos: u64,
    ref_allele: &str,
    alts: &[(usize, &str)],
    svtype: Option<&str>,
    info: &InfoFields<'_>,
    parse: &VcfParse<'_>,
) -> Option<Chopped> {
    let fasta = parse.fasta?;
    if !is_tandem_record(alts, svtype) {
        return None;
    }
    let start = structural_start(pos, VariantClass::TandemRepeat);
    if start > MAX_COORD {
        return None;
    }
    let end = match structural_end(
        pos,
        info.value_written("SVLEN"),
        info.value_written("END"),
        SpanClass::Tandem,
    ) {
        Ok(end) => end,
        // A negative largest SVLEN puts the end before the start (`$start +
        // max(@svlen) - 1`, Parser/VCF.pm:511); with a reference sequence to hand the
        // record is still built, over the REF column (`sub_Slice` yields nothing), and
        // `validate_vf` refuses it as `start > end+1`. An end before the chromosome's
        // first base has no coordinate here and leaves the record to the symbolic
        // path, whose refusal names the same start and end.
        Err(SpanError::StartPastEnd { end, .. }) => u64::try_from(end).ok()?,
        Err(_) => return None,
    };
    let run_bases = fasta
        .sequence(chr, start, end)
        .map(|bases| String::from_utf8_lossy(&bases).into_owned())
        .unwrap_or_else(|| ref_allele.to_string());
    // Without RN the record has one repeat sequence per `CNV:TR` in Perl's type
    // string (`() = $type =~ /CNV:TR/g`, Parser/VCF.pm:397), which is INFO/SVTYPE
    // when the ALTs open outside the VCF 4.4 forms (`:473-474`) and is matched
    // case-sensitively where `get_SO_term` matches the type without case.
    let tandem_alts = sv_type_string(&joined(alts, "/"), svtype)
        .matches("CNV:TR")
        .count();
    let expanded = tandem_repeat_allele_string(info, tandem_alts, parse.max_sv_size, &run_bases)?;
    // The string as Perl's `split` reads it (a trailing empty part is no part), an
    // empty part the `-` allele `_raw_allele_hashes` makes of it
    // (VariationFeatureOverlap.pm:466).
    let mut parts = perl_split(&expanded, &['/']);
    let reference = parts.next().unwrap_or("").to_string();
    let written_case_alts: Vec<String> = parts
        .map(|alt| {
            if alt.is_empty() {
                "-".to_string()
            } else {
                alt.to_string()
            }
        })
        .collect();
    Some(Chopped::unchopped(
        reference,
        written_case_alts,
        expanded.into(),
        start,
        end,
    ))
}

/// `REF/ALT1/ALT2`: the reference and every ALT joined by `/`, one buffer for every
/// allele of the record.
fn allele_string<'a>(reference: &str, alts: impl IntoIterator<Item = &'a str>) -> Arc<str> {
    let mut s = reference.to_string();
    for alt in alts {
        s.push('/');
        s.push_str(alt);
    }
    s.into()
}

/// The name VEP gives a VCF record whose ID column cannot name it
/// (OutputFactory.pm `VariationFeature_to_output_hash`, the VCF `_line` branch):
/// CHROM, POS, REF and ALT as written, the ALT commas as `/`, REF alone when the
/// ALT column is empty.
fn input_line_name(chr: &str, pos: u64, ref_allele: &str, alt_field: &str) -> String {
    if alt_field.is_empty() {
        format!("{chr}_{pos}_{ref_allele}")
    } else {
        format!("{chr}_{pos}_{ref_allele}/{}", alt_field.replace(',', "/"))
    }
}

/// The line Ensembl VEP stores for the record (`_line`, Parser/VCF.pm:349, :541,
/// :573), which its VCF output re-emits (OutputFactory/VCF.pm:315-355) and its JSON
/// output prints as `input` (OutputFactory/JSON.pm:230, :247): the fields Perl's
/// `split` read, joined by tabs, so a trailing empty column is gone. For a record
/// Ensembl VEP routes to its structural path (`structural`, `:237-245`), QUAL, FILTER
/// and INFO are filled with `.` where absent, empty or `0` (`$record->[5] ||= '.'`
/// and the next two columns, `:464-466`, on the array `_line` holds), so a
/// five-column `N <DUP>` is stored with eight columns and a QUAL of `0` as `.`; the
/// sequence path (`:349`) fills nothing.
fn stored_line(columns: &VcfColumns<'_>, line: &str, structural: bool) -> Arc<str> {
    if !structural {
        return Arc::from(line.trim_end_matches('\t'));
    }
    let mut stored: Vec<&str> = columns.fixed.to_vec();
    for column in &mut stored[5..] {
        if !perl_true(column) {
            *column = ".";
        }
    }
    stored.extend(columns.format);
    if !columns.samples.is_empty() {
        stored.push(columns.samples);
    }
    stored.join("\t").into()
}

/// The columns of a VCF data line as Ensembl VEP's parser holds them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfColumns<'a> {
    /// CHROM through INFO. A line that stops after ALT, QUAL or FILTER reads `.` for
    /// the columns it lacks: Ensembl VEP accepts a record of five columns
    /// (Parser/VCF.pm:464-466 fills QUAL, FILTER and INFO; `t/Parser_VCF.t:265-282`),
    /// and a trailing empty column is no column to Perl's `split` (ensembl-io
    /// ColumnBasedParser.pm:250), so `C G\t\t\t` stops after ALT.
    pub fixed: [&'a str; 8],
    /// The FORMAT column, when the line carries one.
    pub format: Option<&'a str>,
    /// The sample columns after FORMAT, tab-joined as written; empty when the line
    /// stops at FORMAT or earlier. [`Self::samples`] splits them.
    samples: &'a str,
}

impl<'a> VcfColumns<'a> {
    /// The INFO column as the view every lookup of the record reads
    /// ([`InfoFields`]), split by the first of them.
    fn info(&self) -> InfoFields<'a> {
        InfoFields::new(self.fixed[7])
    }

    /// The sample columns in their written order (ensembl-io BaseVCF4.pm:928-956 reads
    /// a genotype per sample from them under `--individual`).
    pub fn samples(&self) -> impl Iterator<Item = &'a str> + 'a {
        let samples = self.samples;
        (!samples.is_empty())
            .then(|| samples.split('\t'))
            .into_iter()
            .flatten()
    }
}

/// Splits a VCF data line into its columns as Perl's `split` does, its trailing
/// empty columns dropped. FORMAT and the samples are not split until asked for,
/// since a multisample line can carry thousands of them. A line with fewer than
/// five columns names no ALT and is refused here, where Ensembl VEP reads the
/// missing ALT as an empty string and carries the record as a feature with no
/// allele to annotate.
pub fn vcf_columns(line: &str) -> anyhow::Result<VcfColumns<'_>> {
    let mut parts = line.trim_end_matches('\t').splitn(10, '\t');
    let mut fixed = ["."; 8];
    let mut present = 0;
    for field in fixed.iter_mut() {
        let Some(part) = parts.next() else {
            break;
        };
        *field = part;
        present += 1;
    }
    if present < 5 {
        bail!("VCF line has fewer than 5 fields");
    }
    Ok(VcfColumns {
        fixed,
        format: parts.next(),
        samples: parts.next().unwrap_or(""),
    })
}

/// The sample columns of a `#CHROM` header line, the names Ensembl VEP matches
/// `--individual` against (BaseVCF4.pm:97-98). ensembl-io splits the header on
/// tabs whatever the run's `--delimiter`, which Parser/VCF.pm:129 hands to the
/// data records alone, so the line is read as written, Perl's `split` dropping
/// the empty fields a trailing tab leaves.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SampleColumns {
    /// The sample names in column order: the columns after a ninth column
    /// `FORMAT`, else the ninth column on, as `sample_begin` places them.
    pub names: Vec<String>,
    /// Whether the header's ninth column is `FORMAT`: without it ensembl-io reads
    /// no FORMAT from any record and so no genotype (BaseVCF4.pm:614-618), and a
    /// per-individual run writes no copy.
    pub genotyped: bool,
}

/// The sample columns of a `#CHROM` header line ([`SampleColumns`]).
pub fn sample_columns(header_line: &str) -> SampleColumns {
    let columns: Vec<&str> = perl_split(header_line, &['\t']).collect();
    let genotyped = columns.get(8) == Some(&"FORMAT");
    let first_sample = if genotyped { 9 } else { 8 };
    SampleColumns {
        names: columns
            .into_iter()
            .skip(first_sample)
            .map(str::to_string)
            .collect(),
        genotyped,
    }
}

/// Whether a VCF data line is a non-variant record: its first comma-separated
/// ALT is `.` (Parser/VCF.pm:259 tests `$alts->[0]` of the list ensembl-io
/// BaseVCF4.pm:489 splits on commas), so `.,T` is one and `T,.` is not, and
/// the list is not in a structural shape, which Ensembl VEP tests first
/// (`:237-245`), so `.,<DEL>` is a structural record. A lone `.` under a true
/// INFO/SVTYPE stays a non-variant here ([`non_variant_record`]).
pub(crate) fn vcf_line_is_non_variant(line: &str) -> bool {
    let alt_field = line.split('\t').nth(4).unwrap_or("");
    perl_split(alt_field, &[',']).next() == Some(".")
        && vcf_columns(line).is_ok_and(|columns| columns_are_non_variant(&columns, &columns.info()))
}

/// [`vcf_line_is_non_variant`] over the line's columns and INFO view.
fn columns_are_non_variant(columns: &VcfColumns<'_>, info: &InfoFields<'_>) -> bool {
    let (reference, alt_field) = (columns.fixed[3], columns.fixed[4]);
    if perl_split(alt_field, &[',']).next() != Some(".") {
        return false;
    }
    let svtype = info.value("SVTYPE").filter(|v| perl_true(v));
    let mut alts = written_alts(alt_field);
    variant_alts(&mut alts, svtype);
    alts.len() == 1 || !routes_to_structural(reference, &alts, svtype)
}

/// A sequence record after the anchor chop of Parser/VCF.pm:294-335, upper-cased as
/// `validate_vf` upper-cases the allele string (Parser.pm:554), except
/// `written_case_reference` and `written_case_alts`, which keep the written case the
/// minimisation compares before that upper-casing.
struct Chopped {
    reference: String,
    /// The chopped ALTs in their written order, `*`, `.` and `<DEL:*>` as written.
    alts: Vec<String>,
    start: u64,
    end: u64,
    /// The allele string Ensembl VEP names an ID-less record by: the raw REF and
    /// every ALT joined by `/`. Ensembl VEP keeps the raw string as
    /// `nontrimmed_allele_string` and promotes it to `original_allele_string` for any
    /// record whose alleles differ in length; a record without an indel allele is
    /// never trimmed, so its raw and working strings coincide (Parser/VCF.pm,
    /// Parser.pm `post_process_vfs` and `minimise_alleles`, OutputFactory.pm).
    uploaded: Arc<str>,
    /// Whether the anchor base was chopped (Parser/VCF.pm:309-310, :331-334).
    anchor_chopped: bool,
    /// Ensembl VEP's working string of a record with several ALTs, `*` and `.`
    /// counted (`scalar @$alts > 1`, Parser/VCF.pm:301-322): every allele behind the
    /// chopped anchor, a string `minimise_alleles` leaves as it is; `None` for one ALT.
    working: Option<Arc<str>>,
    /// The chopped reference and ALTs in their written case, which `trim_sequences`
    /// compares with `eq` (Sequence.pm:1006, :1019) before `validate_vf` upper-cases
    /// the string, so `TTGA t` is never minimised to `TGA/-`, and which Parser.pm:462
    /// names an ID-less record by ([`Chopped::written_case_string`]).
    written_case_reference: String,
    written_case_alts: Vec<String>,
}

/// `reference/alt` in the case the input wrote it, when a base is lower case
/// ([`InputVariant::written_case_allele_string`]); `None` for an upper-case pair.
fn written_case_pair(reference: &[u8], alt: &[u8]) -> Option<Arc<str>> {
    reference
        .iter()
        .chain(alt)
        .any(u8::is_ascii_lowercase)
        .then(|| {
            format!(
                "{}/{}",
                String::from_utf8_lossy(reference),
                String::from_utf8_lossy(alt)
            )
            .into()
        })
}

fn chop_record(ref_allele: &str, alts: &[(usize, &str)], pos: u64) -> Chopped {
    let alt_strings: Vec<&str> = alts.iter().map(|(_, alt)| *alt).collect();
    let (written_case_reference, written_case_alts, prefix_trimmed) =
        trim_multi_alleles(ref_allele, &alt_strings);
    let reference = written_case_reference.to_ascii_uppercase();
    let working_alts: Vec<String> = written_case_alts
        .iter()
        .map(|alt| alt.to_ascii_uppercase())
        .collect();
    // Ensembl VEP fixes the end before the chop, `POS + length(REF) - 1`
    // (Parser/VCF.pm:254), and moves the start alone (`:310`), so an empty REF
    // column ends the record one base before POS whether or not it is chopped;
    // here the end stops at 0, where Ensembl VEP carries -1.
    let start = pos + prefix_trimmed as u64;
    let end = (pos + ref_allele.len() as u64).saturating_sub(1);
    let uploaded = allele_string(ref_allele, alt_strings.iter().copied());
    let working = (alts.len() > 1)
        .then(|| allele_string(&reference, working_alts.iter().map(String::as_str)));
    // An ALT the working string keeps empty is the allele `-` to every overlap
    // (VariationFeatureOverlap.pm:462 `$a ||= '-'`).
    let chopped_alts: Vec<String> = working_alts
        .into_iter()
        .map(|alt| if alt.is_empty() { "-".to_string() } else { alt })
        .collect();
    Chopped {
        reference,
        alts: chopped_alts,
        start,
        end,
        uploaded,
        anchor_chopped: prefix_trimmed > 0,
        working,
        written_case_reference,
        written_case_alts,
    }
}

impl Chopped {
    /// A record whose alleles are given over the span `start..=end` with no anchor
    /// to chop: the expansion of a `<CNV:TR>` record, whose feature Ensembl VEP
    /// builds from the run's bases and the repeat sequences as they stand
    /// (Parser/VCF.pm:529-543). The strings are upper-cased as `validate_vf`
    /// upper-cases them, the written case kept for the minimisation's comparison;
    /// `uploaded` is the string as built, which `UPLOADED_ALLELE` prints.
    fn unchopped(
        reference: String,
        alts: Vec<String>,
        uploaded: Arc<str>,
        start: u64,
        end: u64,
    ) -> Chopped {
        let upper_alts: Vec<String> = alts.iter().map(|alt| alt.to_ascii_uppercase()).collect();
        let upper_reference = reference.to_ascii_uppercase();
        let working = (alts.len() > 1)
            .then(|| allele_string(&upper_reference, upper_alts.iter().map(String::as_str)));
        Chopped {
            reference: upper_reference,
            alts: upper_alts,
            start,
            end,
            uploaded,
            anchor_chopped: false,
            working,
            written_case_reference: reference,
            written_case_alts: alts,
        }
    }

    /// The working string in the case the input wrote it, when a base is lower
    /// case ([`InputVariant::written_case_allele_string`]); `None` for an
    /// upper-case record.
    fn written_case_string(&self) -> Option<Arc<str>> {
        let alts = self.written_case_alts.iter().map(String::as_str);
        self.written_case_reference
            .bytes()
            .chain(alts.clone().flat_map(str::bytes))
            .any(|b| b.is_ascii_lowercase())
            .then(|| allele_string(&self.written_case_reference, alts))
    }
}

/// The record of a line none of whose ALTs writes a row, over its chopped (or
/// expanded) alleles: every ALT a `*` or a `<DEL:*>`, equal to the reference or a
/// repeat of an earlier one. Ensembl VEP builds the feature over the chopped alleles
/// and finds no allele to annotate (VariationFeatureOverlap.pm:459-491), so it writes
/// no row in the default and tab formats, carries the VCF line without consequences
/// and writes the JSON object with `most_severe_consequence` `?`:
/// `21 30396602 x TTGA *,*` is the feature `TGA/*/*` at 30396603-30396605 and
/// `21 30289183 x C C,<DEL:*>` the feature `-/-/<DEL:*>` at 30289184-30289183 in
/// Ensembl VEP 116.2's JSON.
fn allele_less_record(chr: String, chopped: Chopped) -> InputVariant {
    let mut variant = InputVariant::new(
        chr,
        chopped.start,
        chopped.end,
        chopped.reference.as_bytes().to_vec(),
        b"*".to_vec(),
    );
    // The feature's string is the chopped one (Parser/VCF.pm:342), which a
    // bi-allelic record's raw string equals unless the anchor was chopped: the
    // one ALT `C*` of a `C` record is the feature `-/*`.
    variant.allele_string = match chopped.working.as_deref() {
        Some(working) => working.to_string(),
        None if chopped.anchor_chopped => {
            allele_string(&chopped.reference, chopped.alts.iter().map(String::as_str)).to_string()
        }
        None => chopped.uploaded.to_ascii_uppercase(),
    };
    variant.no_allele = true;
    variant.written_case_allele_string = chopped.written_case_string();
    variant.uploaded_allele_string = Some(chopped.uploaded);
    variant.record_allele_string_multi = chopped.working;
    variant
}

/// The allele numbers Ensembl VEP gives a sequence record's ALTs when it builds the
/// overlap alleles (VariationFeatureOverlap.pm:459-491), one entry per chopped ALT
/// in written order: the reference is seen first and an allele seen before is dropped
/// without a number (`:461`, an ALT equal to the reference or a repeat of an earlier
/// ALT), every other ALT takes the next number (`:483`), `*` and `<DEL:*>` theirs
/// without a row (`:484`). So `C G,C` numbers `G` 1 alone, `C G,G,T` numbers `G` 1
/// and `T` 2, and `C CG,<DEL:*>` numbers the insertion 1 and the `<DEL:*>` 2, as
/// Ensembl VEP 116.2 prints under `--allele_number`.
pub(crate) fn perl_allele_numbers(reference: &str, alts: &[String]) -> Vec<Option<u32>> {
    let mut seen: Vec<&str> = vec![reference];
    let mut number = 0;
    alts.iter()
        .map(|alt| {
            if seen.contains(&alt.as_str()) {
                return None;
            }
            seen.push(alt);
            number += 1;
            Some(number)
        })
        .collect()
}

/// Whether `post_process_vfs` (Parser.pm:866-881) hands a bi-allelic record to
/// `minimise_alleles`: its alleles as written differ in length and carry no `-`.
fn perl_minimises(ref_allele: &str, alt: &str) -> bool {
    ref_allele.len() != alt.len() && !ref_allele.contains('-') && !alt.contains('-')
}

/// The anchor chop of Parser/VCF.pm:294-335 over every written ALT, `*` and `.` included:
/// the record is an indel when any ALT's length differs from REF's (a `*` of one base
/// against a two-base REF counts, so `AG AC,*` is chopped where `AG AC` is not, :296-298)
/// or when any ALT begins with `D` or `I`, the pileup spelling of a deletion or insertion
/// (`/^[DI]/`, :297), so `DA DG` is chopped to `A/G` at POS + 1 like `DA DC,DG` to `A/C/G`;
/// and the first base is chopped when REF and every ALT other than `*` share it (:306-308
/// with several ALTs, :328 with one), an allele the chop emptied reading `-` (:309, :316,
/// :331-332). With several ALTs one carrying a `*` is kept as written (:315); the one ALT
/// of a record is chopped whatever it carries (:332), so `C C*` reads `-/*` and writes no
/// row. An ALT written empty between two others (`G,,T`) stays empty
/// when nothing is chopped, as :342 joins the parts as they are, so the record's working
/// string reads `C/G//T`; an empty REF reads `-` either way. The comparison is Perl's
/// `eq` on the first characters, case-sensitive (the upper-casing comes later,
/// Parser.pm:554), so `C c,<DEL:*>` is not chopped, and an empty REF is chopped when every
/// ALT carries a `*`, its first character the empty string alone in the vote. Same-length
/// alleles are never trimmed; the secondary minimisation of a bi-allelic indel is
/// `parse_vcf_line`'s (Parser.pm:860-881). Returns the trimmed REF, the trimmed ALTs in
/// their written order and the number of bases chopped.
pub(crate) fn trim_multi_alleles(
    ref_allele: &str,
    alt_alleles: &[&str],
) -> (String, Vec<String>, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let is_indel = alt_alleles
        .iter()
        .any(|a| a.len() != ref_bytes.len() || a.starts_with(['D', 'I']));
    let shares_first_base = |allele: &&str| allele.as_bytes().first() == ref_bytes.first();
    // With one ALT its own first base decides (`<*>` against `N` is not chopped); with
    // several, the alleles carrying a `*` are left out of the comparison.
    let anchored = is_indel
        && match alt_alleles {
            [single] => shares_first_base(single),
            _ => alt_alleles
                .iter()
                .filter(|a| !a.contains('*'))
                .all(shares_first_base),
        };
    let prefix_len = usize::from(anchored);
    // `-` for an allele the chop emptied (:309, :316, :331-332); an ALT written
    // empty stays empty in an unchopped string, which :342 joins as written.
    let chop = |allele: &str, dash_if_empty: bool| {
        let rest: &[u8] = allele.as_bytes().get(prefix_len..).unwrap_or(&[]);
        if rest.is_empty() && dash_if_empty {
            "-".to_string()
        } else {
            String::from_utf8_lossy(rest).into_owned()
        }
    };
    // The `*` exemption is the multi-ALT branch's alone (:315); the one ALT of
    // a record is chopped whatever it carries (:332), so `C C*` reads `-/*`.
    let keep_stars = alt_alleles.len() > 1;
    let trimmed_alts = alt_alleles
        .iter()
        .map(|alt| {
            if keep_stars && alt.contains('*') {
                alt.to_string()
            } else {
                chop(alt, anchored)
            }
        })
        .collect();
    (chop(ref_allele, true), trimmed_alts, prefix_len)
}

/// Convert VCF padded alleles to VEP minimal representation.
///
/// VCF: POS=100, REF=ACGT, ALT=A -> VEP: start=101, end=103, ref=CGT, alt=-
/// VCF: POS=100, REF=A, ALT=ACGT -> VEP: start=101, end=100, ref=-, alt=CGT
/// VCF: POS=100, REF=A, ALT=G -> VEP: start=100, end=100, ref=A, alt=G
#[cfg(test)]
fn vcf_to_vep_alleles(
    ref_allele: &str,
    alt_allele: &str,
    pos: u64,
) -> (Vec<u8>, Vec<u8>, u64, u64) {
    let (trimmed_ref, trimmed_alt, start, end, _) = trim_sequences(
        ref_allele.as_bytes(),
        alt_allele.as_bytes(),
        pos,
        pos + ref_allele.len() as u64 - 1,
    );
    (trimmed_ref, trimmed_alt, start, end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::individuals::{IndividualMode, IndividualRequest, IndividualSplit};
    use vep_core::sv_record::parse_ci_field;
    use vep_core::sv_routing::{classify_symbolic_alt, sv_so_term_and_class, sv_type_abbreviation};

    /// The sample columns of a header whose ninth column is `FORMAT`.
    fn genotyped_samples(names: &[&str]) -> SampleColumns {
        SampleColumns {
            names: names.iter().map(|name| name.to_string()).collect(),
            genotyped: true,
        }
    }

    /// [`parse_vcf_line_with`] under a per-individual `reading` and the default
    /// run's other settings, no line captured.
    fn parse_vcf_record(
        line: &str,
        reading: RecordReading<'_>,
    ) -> anyhow::Result<Vec<InputVariant>> {
        parse_vcf_line_with(
            line,
            &VcfParse {
                reading,
                ..VcfParse::default_run(false)
            },
        )
    }

    /// [`non_variant_record`] under a per-individual `reading` and the default
    /// run's other settings, no line captured.
    fn non_variant_record_read(
        line: &str,
        reading: RecordReading<'_>,
    ) -> anyhow::Result<Vec<InputVariant>> {
        non_variant_record(
            line,
            &VcfParse {
                reading,
                ..VcfParse::default_run(false)
            },
        )
    }

    /// A FASTA of one record `21` reading `bases` from position 1, with its index, in a
    /// fresh temporary directory (kept alive by the returned handle).
    fn fasta_of(bases: &str) -> (tempfile::TempDir, IndexedFasta) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ref.fa");
        std::fs::write(&path, format!(">21\n{bases}\n")).unwrap();
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            format!(
                "21\t{}\t4\t{}\t{}\n",
                bases.len(),
                bases.len(),
                bases.len() + 1
            ),
        )
        .unwrap();
        let fasta = IndexedFasta::from_path(&path).unwrap();
        (dir, fasta)
    }

    /// The run `TCATCATCA` at 11-19 of chromosome 21, behind ten `G`s and before ten
    /// `A`s: the shape of the run at 21:13620284-13620292 of GRCh38, whose padding
    /// base `G` sits at POS 10 here.
    const RUN_FASTA: &str = "GGGGGGGGGGTCATCATCAAAAAAAAAAA";

    /// Parses `line` with `fasta` as the run's reference FASTA under `max_sv_size`.
    fn parse_with_fasta(
        line: &str,
        fasta: Option<&IndexedFasta>,
        max_sv_size: Option<u64>,
    ) -> Vec<InputVariant> {
        parse_vcf_line_with(
            line,
            &VcfParse {
                capture_raw_input: true,
                max_sv_size,
                minimal: false,
                fasta,
                reading: RecordReading::default(),
            },
        )
        .unwrap()
    }

    /// The tandem-repeat warnings the reader writes for `line`, with or without a
    /// FASTA ([`tandem_repeat_warnings`] over the line's columns).
    fn line_warnings(line: &str, fasta_present: bool) -> Vec<String> {
        let columns = vcf_columns(line).unwrap();
        tandem_repeat_warnings(&columns, &columns.info(), fasta_present)
    }

    /// `(location, allele string, variant class, structural)` of each variant.
    fn shapes(variants: &[InputVariant]) -> Vec<(String, String, VariantClass, bool)> {
        variants
            .iter()
            .map(|v| {
                (
                    v.location(),
                    v.allele_string.clone(),
                    v.variant_class,
                    v.is_structural,
                )
            })
            .collect()
    }

    /// Parses `line` under `--minimal` and splits a multi-allelic record as the
    /// reader splits it once validated.
    fn parse_minimal(line: &str) -> Vec<InputVariant> {
        split_minimal_record(
            parse_vcf_line_with(
                line,
                &VcfParse {
                    minimal: true,
                    ..VcfParse::default_run(true)
                },
            )
            .unwrap(),
        )
    }

    #[test]
    fn test_vcf_to_vep_snv() {
        let (r, a, s, e) = vcf_to_vep_alleles("A", "G", 100);
        assert_eq!(r, b"A");
        assert_eq!(a, b"G");
        assert_eq!(s, 100);
        assert_eq!(e, 100);
    }

    #[test]
    fn test_vcf_to_vep_deletion() {
        let (r, a, s, e) = vcf_to_vep_alleles("ACGT", "A", 100);
        assert_eq!(r, b"CGT");
        assert_eq!(a, b"-");
        assert_eq!(s, 101);
        assert_eq!(e, 103);
    }

    #[test]
    fn test_vcf_to_vep_insertion() {
        let (r, a, s, e) = vcf_to_vep_alleles("A", "ACGT", 100);
        assert_eq!(r, b"-");
        assert_eq!(a, b"CGT");
        assert_eq!(s, 101);
        assert_eq!(e, 100);
    }

    #[test]
    fn test_vcf_to_vep_complex_indel() {
        let (r, a, s, e) = vcf_to_vep_alleles("ACG", "AT", 100);
        assert_eq!(r, b"CG");
        assert_eq!(a, b"T");
        assert_eq!(s, 101);
        assert_eq!(e, 102);
    }

    #[test]
    fn test_vcf_to_vep_prefix_and_suffix_trim() {
        let (r, a, s, e) = vcf_to_vep_alleles("ACGT", "ATGT", 100);
        assert_eq!(r, b"C");
        assert_eq!(a, b"T");
        assert_eq!(s, 101);
        assert_eq!(e, 101);
    }

    /// The gate elides the capture: a parser that filled `raw_input` regardless
    /// of it would pass every other test while costing one line-length `String`
    /// allocation per variant.
    #[test]
    fn raw_input_not_captured_when_gated_off() {
        let line = "21\t25000100\trs123\tA\tG\t50\tPASS\t.";
        let variants = parse_vcf_line(line, false).unwrap();
        assert_eq!(variants.len(), 1);
        assert!(
            variants[0].raw_input.is_none(),
            "raw_input must be absent when capture is gated off"
        );
    }

    /// The VCF writer needs the verbatim source line, so capture must retain it
    /// byte-for-byte (it re-emits the record with a CSQ INFO field appended).
    #[test]
    fn raw_input_captured_verbatim_when_gated_on() {
        let line = "21\t25000100\trs123\tA\tG\t50\tPASS\tDP=42";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants[0].raw_input.as_deref(), Some(line));
    }

    /// Multi-allelic records split into one variant per ALT, and the gate must
    /// apply to every split product.
    #[test]
    fn raw_input_gate_applies_to_every_split_allele() {
        let line = "21\t25000100\t.\tA\tG,T\t50\tPASS\t.";
        let off = parse_vcf_line(line, false).unwrap();
        assert_eq!(off.len(), 2);
        assert!(off.iter().all(|v| v.raw_input.is_none()));

        let on = parse_vcf_line(line, true).unwrap();
        assert_eq!(on.len(), 2);
        assert!(on.iter().all(|v| v.raw_input.as_deref() == Some(line)));
    }

    /// Symbolic-allele (SV) records take a separate assignment site in the
    /// parser from explicit alleles, so the gate must hold on both.
    #[test]
    fn raw_input_gate_applies_to_symbolic_alleles() {
        let line = "21\t25000100\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=25001000";
        let off = parse_vcf_line(line, false).unwrap();
        assert!(!off.is_empty());
        assert!(off.iter().all(|v| v.raw_input.is_none()));

        let on = parse_vcf_line(line, true).unwrap();
        assert!(on.iter().all(|v| v.raw_input.as_deref() == Some(line)));
    }

    /// FORMAT and the sample columns name no allele, so a multisample line
    /// must parse exactly as its eight-column prefix does, `raw_input` aside.
    #[test]
    fn multisample_line_parses_like_its_eight_field_prefix() {
        let prefix = "21\t25000100\trs123\tA\tG,T\t50\tPASS\tDP=42;AF=0.5";
        let mut line = format!("{prefix}\tGT:DP");
        for i in 0..3_000 {
            line.push_str(if i % 2 == 0 { "\t0/1:12" } else { "\t1/1:7" });
        }
        let from_prefix = parse_vcf_line(prefix, false).unwrap();
        let from_line = parse_vcf_line(&line, false).unwrap();
        assert_eq!(from_line.len(), 2);
        assert_eq!(
            serde_json::to_string(&from_line).unwrap(),
            serde_json::to_string(&from_prefix).unwrap()
        );
        let captured = parse_vcf_line(&line, true).unwrap();
        assert!(captured
            .iter()
            .all(|v| v.raw_input.as_deref() == Some(line.as_str())));
    }

    /// A line with an INFO column but no FORMAT is the eight-field minimum.
    #[test]
    fn eight_field_line_reads_info() {
        let line = "21\t25000100\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=25001000";
        let variants = parse_vcf_line(line, false).unwrap();
        assert_eq!(variants[0].sv_type.as_deref(), Some("DEL"));
        assert_eq!(variants[0].sv_end, Some(25001000));
    }

    /// Ensembl VEP accepts a VCF record of five columns (Parser/VCF.pm:464-466 fills
    /// QUAL, FILTER and INFO with `.`; its own test `t/Parser_VCF.t:265-282` reads
    /// `21 25587759 test G C` as the feature `G/C` at 25587759): 116.2 annotates
    /// `21 30289183 x C G`, `... G 50` and `... G 50 PASS` as it annotates the
    /// eight-column line, writes `21 30289183 x C G . . CSQ=...` and
    /// `... G 50 PASS CSQ=...` in its VCF output, keeps the short line of a sequence
    /// record as the JSON `input` (`stored_line`), and drops a five-column `N <DEL>`
    /// as `deletion looks incomplete`. A line with fewer than five columns names no
    /// ALT and is refused here (Ensembl VEP carries it as a feature with no allele,
    /// the missing ALT read as an empty string).
    #[test]
    fn five_to_seven_column_records_parse() {
        let full = parse_vcf_line("21\t30289183\tx\tC\tG\t50\tPASS\t.", true).unwrap();
        for short in [
            "21\t30289183\tx\tC\tG",
            "21\t30289183\tx\tC\tG\t50",
            "21\t30289183\tx\tC\tG\t50\tPASS",
        ] {
            let variants = parse_vcf_line(short, true).unwrap();
            assert_eq!(variants.len(), 1, "{short}");
            let (v, f) = (&variants[0], &full[0]);
            assert_eq!(
                (
                    v.chr.as_str(),
                    v.start,
                    v.end,
                    &v.ref_allele,
                    v.alt_allele(),
                    &v.id
                ),
                (
                    f.chr.as_str(),
                    f.start,
                    f.end,
                    &f.ref_allele,
                    f.alt_allele(),
                    &f.id
                ),
                "{short}"
            );
            assert_eq!(v.raw_input.as_deref(), Some(short));
        }
        let columns = vcf_columns("21\t30289183\tx\tC\tG\t50").unwrap();
        assert_eq!(
            columns.fixed,
            ["21", "30289183", "x", "C", "G", "50", ".", "."]
        );
        assert_eq!(columns.format, None);
        assert_eq!(columns.samples().count(), 0);
        let sv = parse_vcf_line("21\t30289183\tx\tN\t<DEL>", true).unwrap();
        assert_eq!(sv.len(), 1);
        assert!(sv[0].is_structural && sv[0].vep_skip);
        assert_eq!((sv[0].start, sv[0].end), (30289184, 30289184));
        for too_short in ["21\t25000100\trs123\tA", "21\t25000100", ""] {
            let err = parse_vcf_line(too_short, false).unwrap_err();
            assert_eq!(
                err.to_string(),
                "VCF line has fewer than 5 fields",
                "{too_short:?}"
            );
        }
    }

    /// The stored line is the one Ensembl VEP 116.2 echoes. Trailing empty columns are
    /// no columns (ensembl-io ColumnBasedParser.pm:250): `21 30289183 x C G,G,T` with
    /// three trailing tabs is a five-column record whose VCF line prints `. . CSQ=...`
    /// and whose JSON `input` ends at the ALT. The structural path fills QUAL, FILTER
    /// and INFO (Parser/VCF.pm:464-466, `||=`) on the array `_line` holds: the JSON
    /// `input` of `21 30165850 x N <DUP>` is `21 30165850 x N <DUP> . . .`, a
    /// structural QUAL of `0` prints `.` in the VCF line and the `input`
    /// (`N <DEL> 0 PASS SVTYPE=DEL;END=30165950` prints `. PASS`), an empty or `0`
    /// FILTER prints `.`, and `0.0`, Perl-true, stays; the sequence path fills nothing,
    /// so `C G 0 PASS .` keeps its `0` and `C G` with an INFO/SVTYPE, a run of `ACGT`
    /// and so a sequence record (`:237-245`), keeps an empty FILTER. A record Ensembl
    /// VEP routes to its structural path by its INFO/SVTYPE alone (`N G`, `N *`) is
    /// filled whether or not it names an allele here.
    #[test]
    fn the_stored_line_folds_trailing_columns_and_fills_a_structural_record() {
        let trailing = "21\t30289183\ttrailing_tabs\tC\tG,G,T\t\t\t";
        let columns = vcf_columns(trailing).unwrap();
        assert_eq!(
            columns.fixed,
            [
                "21",
                "30289183",
                "trailing_tabs",
                "C",
                "G,G,T",
                ".",
                ".",
                "."
            ]
        );
        assert_eq!(columns.format, None);
        let stored = |line: &str| {
            let variants = parse_vcf_line(line, true).unwrap();
            assert!(!variants.is_empty(), "{line}");
            variants[0].raw_input.as_deref().unwrap().to_string()
        };
        assert_eq!(stored(trailing), "21\t30289183\ttrailing_tabs\tC\tG,G,T");
        assert_eq!(
            stored("21\t30289183\tx\tC\tG\t50\t"),
            "21\t30289183\tx\tC\tG\t50"
        );
        for (line, expected) in [
            (
                "21\t30165850\tdup_short\tN\t<DUP>",
                "21\t30165850\tdup_short\tN\t<DUP>\t.\t.\t.",
            ),
            (
                "21\t30165850\tsv_qual0\tN\t<DEL>\t0\tPASS\tSVTYPE=DEL;END=30165950",
                "21\t30165850\tsv_qual0\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30165950",
            ),
            (
                "21\t30165850\tx\tN\t<DUP>\t50\t\tSVTYPE=DUP;END=30165950",
                "21\t30165850\tx\tN\t<DUP>\t50\t.\tSVTYPE=DUP;END=30165950",
            ),
            (
                "21\t30165850\tx\tN\t<DUP>\t50\t0",
                "21\t30165850\tx\tN\t<DUP>\t50\t.\t.",
            ),
            (
                "21\t30165850\tx\tN\t<DUP>\t0.0\t.\t.",
                "21\t30165850\tx\tN\t<DUP>\t0.0\t.\t.",
            ),
            (
                "21\t30165850\tx\tN\t<DUP>\t0\tPASS\tSVTYPE=DUP;END=30165950\tGT\t0/1\t1/1",
                "21\t30165850\tx\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=30165950\tGT\t0/1\t1/1",
            ),
            (
                "21\t30289183\tseq_qual0\tC\tG\t0\tPASS\t.",
                "21\t30289183\tseq_qual0\tC\tG\t0\tPASS\t.",
            ),
            (
                "21\t30289183\tx\tN\tG\t0\t\tSVTYPE=DEL;END=30289190",
                "21\t30289183\tx\tN\tG\t.\t.\tSVTYPE=DEL;END=30289190",
            ),
            (
                "21\t30289183\tx\tC\tG\t0\t\tSVTYPE=DEL;END=30289190",
                "21\t30289183\tx\tC\tG\t0\t\tSVTYPE=DEL;END=30289190",
            ),
            (
                "21\t100\t.\tN\t*\t0\tPASS\tSVTYPE=DEL;END=400",
                "21\t100\t.\tN\t*\t.\tPASS\tSVTYPE=DEL;END=400",
            ),
        ] {
            assert_eq!(stored(line), expected, "{line}");
        }
        let non_variant =
            non_variant_record("21\t30289183\tx\tTG\t.\t\t\t", &VcfParse::default_run(true))
                .unwrap();
        assert_eq!(
            non_variant[0].raw_input.as_deref(),
            Some("21\t30289183\tx\tTG\t.")
        );
        let uncaptured = parse_vcf_line("21\t30165850\tx\tN\t<DUP>", false).unwrap();
        assert!(uncaptured[0].raw_input.is_none());
    }

    /// The parser keeps FORMAT and the sample columns for the per-sample modes and
    /// reads the sample names from the `#CHROM` header, as ensembl-io
    /// BaseVCF4.pm:928-956 reads a genotype per named sample under `--individual`.
    /// The columns are split on demand, so a line with thousands of samples costs
    /// the allele parser nothing.
    #[test]
    fn sample_columns_retained_and_header_names_captured() {
        let line = "21\t30289183\tx\tC\tG\t50\tPASS\tDP=4\tGT:DP\t0/1:3\t1/1:5";
        let columns = vcf_columns(line).unwrap();
        assert_eq!(
            columns.fixed,
            ["21", "30289183", "x", "C", "G", "50", "PASS", "DP=4"]
        );
        assert_eq!(columns.format, Some("GT:DP"));
        assert_eq!(columns.samples().collect::<Vec<_>>(), ["0/1:3", "1/1:5"]);
        let format_only = vcf_columns("21\t30289183\tx\tC\tG\t50\tPASS\tDP=4\tGT").unwrap();
        assert_eq!(format_only.format, Some("GT"));
        assert_eq!(format_only.samples().count(), 0);
        let one_sample = vcf_columns("21\t30289183\tx\tC\tG\t.\t.\t.\tGT\t0/1").unwrap();
        assert_eq!(one_sample.samples().collect::<Vec<_>>(), ["0/1"]);
        let header =
            sample_columns("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2");
        assert_eq!(header.names, ["S1", "S2"]);
        assert!(header.genotyped);
        // Perl's `split("\t")` (BaseVCF4.pm:97) drops the empty fields a trailing tab
        // leaves, so a header ending in a tab names no empty sample.
        let trailing_tab =
            sample_columns("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\t\t");
        assert_eq!(trailing_tab.names, ["S1", "S2"]);
        assert!(
            sample_columns("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
                .names
                .is_empty()
        );
        assert!(
            sample_columns("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
                .names
                .is_empty()
        );
        // BaseVCF4.pm:98: a ninth column that is not FORMAT begins the samples, and
        // no record then carries a FORMAT to read a genotype from.
        let no_format = sample_columns("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2\tS3");
        assert_eq!(no_format.names, ["S1", "S2", "S3"]);
        assert!(!no_format.genotyped);
        // The alleles are read from the fixed columns alone.
        let variants = parse_vcf_line(line, false).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].alt_allele(), b"G");
    }

    /// `<DEL:*>` on the sequence path (Parser/VCF.pm:237-245 route a record there when
    /// no bracket token is free of `*` and INFO/SVTYPE is not true) is an allele Ensembl
    /// VEP numbers and writes no row for (VariationFeatureOverlap.pm:483-484, beside
    /// `*`), left whole by the anchor chop (`:315`). Ensembl VEP 116.2 writes for
    /// `21 30289183 x C CG,<DEL:*>` one row, `21:30289183-30289184 G` with ALLELE_NUM 1
    /// and VARIANT_CLASS `sequence_alteration`, and for `C G,<DEL:*>` one row `G`,
    /// ALLELE_NUM 1; under `SVTYPE=DEL` with an END the token is the record's structural
    /// deletion (`21:30165851-30165950 deletion`).
    #[test]
    fn del_star_kept_as_small_variant_allele() {
        let ins = parse_vcf_line("21\t30289183\tx\tC\tCG,<DEL:*>\t.\t.\t.", true).unwrap();
        assert_eq!(ins.len(), 1);
        let v = &ins[0];
        assert!(!v.is_structural && !v.no_allele);
        assert_eq!((v.start, v.end), (30289184, 30289183));
        assert_eq!(
            (v.ref_allele.as_slice(), v.alt_allele()),
            (&b"-"[..], &b"G"[..])
        );
        assert_eq!((v.allele_index, v.allele_num), (0, None));
        assert_eq!(v.uploaded_allele_string.as_deref(), Some("C/CG/<DEL:*>"));
        assert_eq!(v.record_allele_string(), "-/G/<DEL:*>");
        assert_eq!(v.class_so_term(), "sequence_alteration");

        let snv = parse_vcf_line("21\t30289183\tx\tC\tG,<DEL:*>\t.\t.\t.", true).unwrap();
        assert_eq!(snv.len(), 1);
        assert_eq!(
            (snv[0].start, snv[0].end, snv[0].alt_allele()),
            (30289183, 30289183, &b"G"[..])
        );
        assert_eq!(snv[0].record_allele_string(), "C/G/<DEL:*>");
        assert_eq!(snv[0].class_so_term(), "sequence_alteration");

        let structural = parse_vcf_line(
            "21\t30165850\tx\tN\t<DEL:*>\t.\t.\tSVTYPE=DEL;END=30165950",
            true,
        )
        .unwrap();
        assert_eq!(structural.len(), 1);
        assert!(structural[0].is_structural && !structural[0].vep_skip);
        assert_eq!(
            structural[0].variant_class,
            VariantClass::StructuralDeletion
        );
        assert_eq!(
            (structural[0].start, structural[0].end),
            (30165851, 30165950)
        );
    }

    /// A record whose only ALT is `<DEL:*>` has no allele to annotate: Ensembl VEP
    /// 116.2 writes no row for `21 30289183 x C <DEL:*>` in the default and tab
    /// formats, carries its VCF line without CSQ and writes the JSON object
    /// `allele_string` `C/<DEL:*>`, 30289183-30289183, `most_severe_consequence` `?`,
    /// `variant_class` `deletion`; `C C,<DEL:*>` is chopped first (`:306-316`) and is
    /// the feature `-/-/<DEL:*>` at 30289184-30289183.
    #[test]
    fn lone_del_star_no_rows() {
        let lone = parse_vcf_line("21\t30289183\tx\tC\t<DEL:*>\t.\t.\t.", true).unwrap();
        assert_eq!(lone.len(), 1);
        let v = &lone[0];
        assert!(v.no_allele && !v.is_structural);
        assert_eq!((v.start, v.end), (30289183, 30289183));
        assert_eq!(v.allele_string, "C/<DEL:*>");
        assert_eq!(v.record_allele_string(), "C/<DEL:*>");
        assert_eq!(v.class_so_term(), "deletion");
        assert_eq!(v.uploaded_variation(), "x");

        let chopped = parse_vcf_line("21\t30289183\t.\tC\tC,<DEL:*>\t.\t.\t.", true).unwrap();
        assert_eq!(chopped.len(), 1);
        let v = &chopped[0];
        assert!(v.no_allele);
        assert_eq!((v.start, v.end), (30289184, 30289183));
        assert_eq!(v.record_allele_string(), "-/-/<DEL:*>");
        assert_eq!(v.uploaded_allele_string.as_deref(), Some("C/C/<DEL:*>"));
        assert_eq!(v.uploaded_variation(), "21_30289183_C/C/<DEL:*>");

        // The chop compares the first characters as written (Parser/VCF.pm:306, before
        // the upper-casing at Parser.pm:554): Ensembl VEP 116.2 leaves `C c,<DEL:*>`
        // unchopped, the JSON object `C/C/<DEL:*>` at 30289183-30289183.
        let lower = parse_vcf_line("21\t30289183\t.\tC\tc,<DEL:*>\t.\t.\t.", true).unwrap();
        assert!(lower[0].no_allele);
        assert_eq!((lower[0].start, lower[0].end), (30289183, 30289183));
        assert_eq!(lower[0].record_allele_string(), "C/C/<DEL:*>");
    }

    /// An empty REF column (no base, not valid VCF) keeps Ensembl VEP's end, POS - 1
    /// (Parser/VCF.pm:254), whether or not the record is chopped: with every ALT a `*`
    /// the vote holds the empty string alone and the start moves (`:306-310`), so
    /// `21 30289183 x <empty> *,*` is `-/*/*` at 30289184-30289182, the span
    /// `validate_vf` refuses as `start > end+1 : (START=30289184, END=30289182)` on both
    /// engines; `<empty> G` is `-/G` at 30289183-30289182 on both, and at POS 1 the
    /// insertion `-/G` at 1-0, an end `minimise_alleles` never resets since the
    /// unchopped `/G` fails its `.+/.+` test (Parser.pm:974).
    #[test]
    fn an_empty_ref_keeps_ensembl_veps_end() {
        let stars = parse_vcf_line("21\t30289183\tx\t\t*,*\t.\t.\t.", true).unwrap();
        assert_eq!(stars.len(), 1);
        assert!(stars[0].no_allele);
        assert_eq!((stars[0].start, stars[0].end), (30289184, 30289182));
        assert_eq!(stars[0].record_allele_string(), "-/*/*");
        let insertion = parse_vcf_line("21\t30289183\tx\t\tG\t.\t.\t.", true).unwrap();
        assert_eq!((insertion[0].start, insertion[0].end), (30289183, 30289182));
        assert_eq!(insertion[0].allele_string, "-/G");
        let at_one = parse_vcf_line("21\t1\tx\t\tG\t.\t.\t.", true).unwrap();
        assert_eq!((at_one[0].start, at_one[0].end), (1, 0));
        assert_eq!(at_one[0].allele_string, "-/G");
        assert_eq!(
            trim_multi_alleles("", &["*", "*"]),
            ("-".to_string(), vec!["*".to_string(), "*".to_string()], 1)
        );
        assert_eq!(
            trim_multi_alleles("TTGA", &["T", "ttga"]),
            (
                "TTGA".to_string(),
                vec!["T".to_string(), "ttga".to_string()],
                0
            )
        );
    }

    /// An ALT equal to the reference is dropped when Ensembl VEP builds the overlap
    /// alleles (VariationFeatureOverlap.pm:461 and :481, before the counter at :483), so
    /// it takes no row and no allele number, after the anchor chop and the upper-casing
    /// of `validate_vf` (Parser.pm:554). Ensembl VEP 116.2 writes for `21 30289183 x C C`
    /// no row and the JSON object `C/C`, `?`, `variant_class` `SNV`; for `C G,C` the one
    /// row `G` with ALLELE_NUM 1; for `21 30396602 x TTGA T,TTGA` the one row
    /// `21:30396603-30396605 -` with ALLELE_NUM 1 under the string `TGA/-/TGA`; for
    /// `C CG,C` the one row `21:30289183-30289184 G`; for `TTGA TTGA` no row and the
    /// object `TTGA/TTGA`, `substitution`; and for `C c` no row and the object `C/C`.
    #[test]
    fn alt_equal_to_ref_dropped_no_row() {
        let same = parse_vcf_line("21\t30289183\tx\tC\tC\t.\t.\t.", true).unwrap();
        assert_eq!(same.len(), 1);
        assert!(same[0].no_allele);
        assert_eq!((same[0].start, same[0].end), (30289183, 30289183));
        assert_eq!(same[0].record_allele_string(), "C/C");
        assert_eq!(same[0].class_so_term(), "SNV");

        let second = parse_vcf_line("21\t30289183\tx\tC\tG,C\t.\t.\t.", true).unwrap();
        assert_eq!(second.len(), 1);
        assert_eq!(second[0].alt_allele(), b"G");
        assert_eq!((second[0].allele_index, second[0].allele_num), (0, None));
        assert_eq!(second[0].record_allele_string(), "C/G/C");

        let deletion = parse_vcf_line("21\t30396602\tx\tTTGA\tT,TTGA\t.\t.\t.", true).unwrap();
        assert_eq!(deletion.len(), 1);
        let v = &deletion[0];
        assert_eq!((v.start, v.end), (30396603, 30396605));
        assert_eq!(
            (v.ref_allele.as_slice(), v.alt_allele()),
            (&b"TGA"[..], &b"-"[..])
        );
        assert_eq!((v.allele_index, v.allele_num), (0, None));
        assert_eq!(v.record_allele_string(), "TGA/-/TGA");
        assert_eq!(v.class_so_term(), "sequence_alteration");

        let insertion = parse_vcf_line("21\t30289183\tx\tC\tCG,C\t.\t.\t.", true).unwrap();
        assert_eq!(insertion.len(), 1);
        assert_eq!((insertion[0].start, insertion[0].end), (30289184, 30289183));
        assert_eq!(insertion[0].alt_allele(), b"G");
        assert_eq!(insertion[0].record_allele_string(), "-/G/-");

        let mnp = parse_vcf_line("21\t30396602\tx\tTTGA\tTTGA\t.\t.\t.", true).unwrap();
        assert!(mnp[0].no_allele);
        assert_eq!((mnp[0].start, mnp[0].end), (30396602, 30396605));
        assert_eq!(mnp[0].record_allele_string(), "TTGA/TTGA");
        assert_eq!(mnp[0].class_so_term(), "substitution");

        let lower = parse_vcf_line("21\t30289183\t.\tC\tc\t.\t.\t.", true).unwrap();
        assert!(lower[0].no_allele);
        assert_eq!(lower[0].record_allele_string(), "C/C");
        assert_eq!(lower[0].uploaded_variation(), "21_30289183_C/c");
    }

    /// A repeated ALT is dropped at VariationFeatureOverlap.pm:461 before the counter,
    /// so the alleles after it take the next numbers: Ensembl VEP 116.2 writes for
    /// `21 30289183 x C G,G,T` the rows `G` ALLELE_NUM 1 and `T` ALLELE_NUM 2, for
    /// `C G,T,G` the same two rows, for `C G,G` and `C g,G` the one row `G`, and keeps
    /// the written string (`C/G/G/T`) as the record's allele string; `*` and `<DEL:*>`
    /// take a number without a row (`:483-484`), so `C *,G` numbers `G` 2.
    #[test]
    fn duplicate_alts_deduplicated_allele_num() {
        let repeated = parse_vcf_line("21\t30289183\tx\tC\tG,G,T\t.\t.\t.", true).unwrap();
        assert_eq!(
            repeated
                .iter()
                .map(|v| (v.alt_allele(), v.allele_index, v.allele_num))
                .collect::<Vec<_>>(),
            [(&b"G"[..], 0, None), (&b"T"[..], 2, Some(2))]
        );
        assert!(repeated
            .iter()
            .all(|v| v.record_allele_string() == "C/G/G/T"));

        let apart = parse_vcf_line("21\t30289183\tx\tC\tG,T,G\t.\t.\t.", true).unwrap();
        assert_eq!(
            apart
                .iter()
                .map(|v| (v.alt_allele(), v.allele_index, v.allele_num))
                .collect::<Vec<_>>(),
            [(&b"G"[..], 0, None), (&b"T"[..], 1, None)]
        );
        for line in [
            "21\t30289183\tx\tC\tG,G\t.\t.\t.",
            "21\t30289183\tx\tC\tg,G\t.\t.\t.",
        ] {
            let variants = parse_vcf_line(line, true).unwrap();
            assert_eq!(variants.len(), 1, "{line}");
            assert_eq!(variants[0].alt_allele(), b"G");
            assert_eq!(variants[0].record_allele_string(), "C/G/G");
        }
        let star_first = parse_vcf_line("21\t30289183\tx\tC\t*,G\t.\t.\t.", true).unwrap();
        assert_eq!(star_first.len(), 1);
        assert_eq!(
            (star_first[0].allele_index, star_first[0].allele_num),
            (1, None)
        );
        assert_eq!(
            perl_allele_numbers("C", &["*".into(), "G".into()]),
            [Some(1), Some(2)]
        );
        assert_eq!(
            perl_allele_numbers("-", &["T".into(), "-".into(), "T".into()]),
            [Some(1), None, None]
        );
    }

    /// Parser/VCF.pm:297 reads an ALT beginning with `D` or `I` as an indel whatever its
    /// length, so the anchor chop runs on same-length alleles behind that letter: Ensembl
    /// VEP 116.2 annotates `21 30396602 x DA DG` as the SNV `A/G` at 30396603 (`:328-334`)
    /// and `21 30289182 x DA DC,DG` as the SNVs `A/C/G` at 30289183 (`:306-308`), the
    /// uploaded alleles kept as written; `DA DG` without the letter, `TA TG`, is not
    /// chopped.
    #[test]
    fn an_alt_beginning_with_d_or_i_is_chopped_as_an_indel() {
        assert_eq!(
            trim_multi_alleles("DA", &["DG"]),
            ("A".to_string(), vec!["G".to_string()], 1)
        );
        assert_eq!(
            trim_multi_alleles("DA", &["DC", "DG"]),
            ("A".to_string(), vec!["C".to_string(), "G".to_string()], 1)
        );
        assert_eq!(
            trim_multi_alleles("TA", &["TG"]),
            ("TA".to_string(), vec!["TG".to_string()], 0)
        );

        let single = parse_vcf_line("21\t30396602\tx\tDA\tDG\t.\t.\t.", true).unwrap();
        assert_eq!(single.len(), 1);
        assert_eq!((single[0].start, single[0].end), (30396603, 30396603));
        assert_eq!(single[0].allele_string, "A/G");
        assert_eq!(single[0].variant_class, VariantClass::Snv);
        assert_eq!(single[0].uploaded_allele_string.as_deref(), Some("DA/DG"));

        let multi = parse_vcf_line("21\t30289182\tx\tDA\tDC,DG\t.\t.\t.", true).unwrap();
        assert_eq!(
            multi
                .iter()
                .map(|v| (v.start, v.end, v.alt_allele(), v.variant_class))
                .collect::<Vec<_>>(),
            [
                (30289183, 30289183, &b"C"[..], VariantClass::Snv),
                (30289183, 30289183, &b"G"[..], VariantClass::Snv),
            ]
        );
        assert!(multi.iter().all(|v| v.record_allele_string() == "A/C/G"));
        assert!(multi
            .iter()
            .all(|v| v.uploaded_allele_string.as_deref() == Some("DA/DC/DG")));
    }

    /// Parser/VCF.pm:555-558 marks a deletion with neither END nor SVLEN `vep_skip`
    /// and warns `deletion looks incomplete`; at 116.2 `validate_svf` (Parser.pm:821-826)
    /// then drops the record unless `--dont_skip` keeps it. The test runs on the term
    /// (`deletion` for `<DEL>` and `<CN0>`, `mobile_element_deletion` for
    /// `<DEL:ME:ALU>`) or on the type string of a type `%SO_TERMS` lacks
    /// (`<INS>/<DEL>`), so Ensembl VEP 116.2 drops `21 30165850 x N <DEL>` with
    /// `SVTYPE=DEL`, keeps it with `SVLEN=-100` or `END=30165851`, and never marks a
    /// `<DUP>` or `<INV>` for a missing extent.
    #[test]
    fn del_without_end_or_svlen_flagged() {
        let parse = |alt: &str, info: &str| {
            parse_vcf_line(&format!("21\t30165850\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap()
        };
        let bare = parse("<DEL>", "SVTYPE=DEL");
        assert!(bare[0].vep_skip && bare[0].is_structural);
        assert_eq!((bare[0].start, bare[0].end), (30165851, 30165851));
        assert_eq!(bare[0].record_so_term.as_deref(), Some("deletion"));
        for (alt, info) in [
            ("<DEL>", "."),
            ("<CN0>", "."),
            ("<DEL:ME:ALU>", "SVTYPE=DEL"),
            ("<INS>,<DEL>", "."),
        ] {
            assert!(parse(alt, info)[0].vep_skip, "{alt} {info}");
        }
        for (alt, info) in [
            ("<DEL>", "SVLEN=-100"),
            ("<DEL>", "END=30165851"),
            ("<DEL>", "SVTYPE=DEL;END=30165950"),
            ("<DUP>", "."),
            ("<INV>", "SVTYPE=INV"),
            ("<INS>", "."),
        ] {
            assert!(!parse(alt, info)[0].vep_skip, "{alt} {info}");
        }
        // The unsupported type with an extent keeps the type's mark alone.
        let unsupported = parse("<INS>,<DEL>", "END=30166000");
        assert!(unsupported
            .iter()
            .all(|v| v.vep_skip && v.record_so_term.is_none()));
    }

    /// A deletion whose END is POS has an empty span (`$start > $end`, Parser/VCF.pm:555):
    /// Ensembl VEP 116.2 marks it, warns `deletion looks incomplete` and drops
    /// `21 30165850 x N <DEL> . . END=30165850`, which `--dont_skip` keeps as the pair
    /// 30165851-30165850; an END before POS fails `validate_vf`'s own coordinate check
    /// too (`start > end+1 : (START=30165851, END=30165800)`, Parser.pm:484-486) and is
    /// refused here, as Ensembl VEP writes both warnings; an END at POS + 1 is a
    /// one-base deletion and is kept.
    #[test]
    fn del_end_le_pos_flagged() {
        let at_pos = parse_vcf_line("21\t30165850\tx\tN\t<DEL>\t.\t.\tEND=30165850", true).unwrap();
        assert!(at_pos[0].vep_skip);
        assert_eq!((at_pos[0].start, at_pos[0].end), (30165851, 30165850));
        let after_pos =
            parse_vcf_line("21\t30165850\tx\tN\t<DEL>\t.\t.\tEND=30165851", true).unwrap();
        assert!(!after_pos[0].vep_skip);
        assert_eq!((after_pos[0].start, after_pos[0].end), (30165851, 30165851));
        let before_pos =
            parse_vcf_line("21\t30165850\tx\tN\t<DEL>\t.\t.\tEND=30165800", true).unwrap_err();
        assert_eq!(
            before_pos.to_string(),
            "start > end+1 : (START=30165851, END=30165800)"
        );
    }

    /// POS 0 is the telomere's position in VCF. Ensembl VEP 116.2 annotates `21 0 x A G`
    /// at `21:0` and `21 0 x AT A` at `21:1`, and skips `21 0 x A AT` with `Alleles look
    /// like an insertion (-/T) but coordinates are not start = end + 1 (START=1, END=1)`,
    /// because `trim_sequences` (Sequence.pm:971) reads the chopped record's end of 0 as
    /// no end and resets it to the start; the end an empty REF column gives at POS 0,
    /// -1 in Ensembl VEP, is 0 here, the insertion check refusing the record where
    /// Ensembl VEP writes an intergenic row at `21:-1-0`.
    #[test]
    fn pos_zero_records_follow_ensembl_vep_without_underflow() {
        let snv = parse_vcf_line("21\t0\tx\tA\tG\t.\t.\t.", true).unwrap();
        assert_eq!((snv[0].start, snv[0].end), (0, 0));
        let deletion = parse_vcf_line("21\t0\tx\tAT\tA\t.\t.\t.", true).unwrap();
        assert_eq!((deletion[0].start, deletion[0].end), (1, 1));
        assert_eq!(deletion[0].allele_string, "T/-");
        let insertion = parse_vcf_line("21\t0\tx\tA\tAT\t.\t.\t.", true).unwrap();
        assert_eq!((insertion[0].start, insertion[0].end), (1, 1));
        assert_eq!(insertion[0].allele_string, "-/T");
        let later = parse_vcf_line("21\t5\tx\tA\tAT\t.\t.\t.", true).unwrap();
        assert_eq!((later[0].start, later[0].end), (6, 5));
        let empty_ref = parse_vcf_line("21\t0\tx\t\tA\t.\t.\t.", true).unwrap();
        assert_eq!((empty_ref[0].start, empty_ref[0].end), (0, 0));
        assert_eq!(empty_ref[0].allele_string, "-/A");
        let star = parse_vcf_line("21\t0\tx\tA\t*\t.\t.\t.", true).unwrap();
        assert!(star[0].no_allele);
        assert_eq!((star[0].start, star[0].end), (0, 0));
        let empty_star = parse_vcf_line("21\t0\tx\t\t*\t.\t.\t.", true).unwrap();
        assert!(empty_star[0].no_allele);
        assert_eq!(empty_star[0].allele_string, "/*");
    }

    #[test]
    fn test_parse_vcf_line_snv() {
        let line = "21\t25000100\trs123\tA\tG\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].chr, "21");
        assert_eq!(variants[0].start, 25000100);
        assert_eq!(variants[0].end, 25000100);
        assert_eq!(variants[0].ref_allele, b"A");
        assert_eq!(variants[0].alt_allele(), b"G");
        assert_eq!(variants[0].id.as_deref(), Some("rs123"));
    }

    #[test]
    fn test_parse_vcf_line_multiallelic() {
        let line = "21\t25000100\t.\tA\tG,T\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].display_allele(), "G");
        assert_eq!(variants[1].display_allele(), "T");
        assert_eq!(variants[0].allele_index, 0);
        assert_eq!(variants[1].allele_index, 1);
    }

    /// An empty ALT between two others, `C G,,T`: the first characters `C`, `G`,
    /// none and `T` are not one, so nothing is chopped and Parser/VCF.pm:342 joins
    /// the parts as written, the working string `C/G//T` the JSON `allele_string`
    /// prints; the empty part is the allele `-` to every overlap
    /// (VariationFeatureOverlap.pm:462), numbered 2 between `G` and `T`. Ensembl
    /// VEP 116.2 writes three rows, `G`, `-` and `T`, each `UPLOADED_ALLELE=C/G//T`.
    #[test]
    fn an_empty_alt_between_others_stays_empty_in_the_working_string() {
        let variants = parse_vcf_line("21\t30289183\t.\tC\tG,,T\t.\t.\t.", true).unwrap();
        let alleles: Vec<(String, usize, Option<u32>)> = variants
            .iter()
            .map(|v| (v.display_allele(), v.allele_index, v.allele_num))
            .collect();
        assert_eq!(
            alleles,
            [
                ("G".to_string(), 0, None),
                ("-".to_string(), 1, None),
                ("T".to_string(), 2, None),
            ]
        );
        assert_eq!(variants[1].record_allele_string(), "C/G//T");
        assert_eq!(variants[1].allele_string, "C/-");
        assert_eq!(
            variants[1].uploaded_allele_string.as_deref(),
            Some("C/G//T")
        );
        assert_eq!((variants[1].start, variants[1].end), (30289183, 30289183));

        // The chop keeps `-` for an allele it emptied (:316).
        let chopped = parse_vcf_line("21\t30289183\t.\tCT\tC,CG\t.\t.\t.", true).unwrap();
        assert_eq!(chopped[0].record_allele_string(), "T/-/G");
    }

    #[test]
    fn test_parse_vcf_line_multiallelic_indel_shared_window() {
        let line = "1\t100\t.\tCTTT\tC,CTT\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        // ALT C: prefix C trimmed -> TTT/- (deletion, alt is "-" -> no secondary trim)
        assert_eq!(variants[0].start, 101);
        assert_eq!(variants[0].end, 103);
        assert_eq!(variants[0].ref_allele, b"TTT");
        assert_eq!(variants[0].alt_allele(), b"-");

        // ALT CTT: prefix C chopped -> TTT/TT and kept there: a multi-allelic
        // record is never minimised further, so VEP computes this allele's
        // positions and distances on the three-base window.
        assert_eq!(variants[1].start, 101);
        assert_eq!(variants[1].end, 103);
        assert_eq!(variants[1].ref_allele, b"TTT");
        assert_eq!(variants[1].alt_allele(), b"TT");
        assert!(!variants[1].minimised);
        assert_eq!(variants[1].uploaded_variation(), "1_100_CTTT/C/CTT");
        assert_eq!(variants[1].record_allele_string(), "TTT/-/TT");
    }

    /// `Uploaded_variation` of every allele parsed from `line`.
    fn uploaded_variations(line: &str) -> Vec<String> {
        parse_vcf_line(line, false)
            .unwrap()
            .iter()
            .map(|v| v.uploaded_variation())
            .collect()
    }

    /// OutputFactory.pm `VariationFeature_to_output_hash`: the ID column names
    /// the record (the first of `;`-joined IDs, as Parser/VCF.pm takes it).
    #[test]
    fn uploaded_variation_is_the_id_when_the_record_carries_one() {
        assert_eq!(
            uploaded_variations("21\t9937918\trs867932594\tGGAGGGTCTGAATGTTTGTCCC\tG\t.\t.\t."),
            ["rs867932594"]
        );
        assert_eq!(
            uploaded_variations("21\t34971477\trs1;rs2\tC\tA\t.\t.\t."),
            ["rs1"]
        );
    }

    /// Only the first `;`-joined ID is read (Parser/VCF.pm:347 `$ids->[0]`), and the
    /// record is named by its line when that element is `.` (OutputFactory.pm:900), empty
    /// or `0`: Parser.pm:462 `||=` replaces a Perl-false name with `chr_start_alleles`,
    /// here `21_34836522_GAGACGGA/-`, which OutputFactory.pm:900 rejects as the generated
    /// name. On the GRCh37 corpus record `21 34836521 TGAGACGGA T` the IDs `.;rs2`,
    /// `;rs2`, `0` and `0;rs2` therefore all yield the line name
    /// `21_34836521_TGAGACGGA/T`, the name Ensembl VEP 116.2 writes for the record's
    /// `.` ID, and `rs2;.` yields `rs2`. The JSON `id` (OutputFactory/JSON.pm:236) is the
    /// stored name after the `||=`: `.` for `.;rs2`, `21_34836522_GAGACGGA/-` for `;rs2`,
    /// `0` and `0;rs2`, `rs2` for `rs2;.`.
    #[test]
    fn a_dot_empty_or_zero_first_id_names_the_record_by_its_line() {
        for (id, json_id) in [
            (".;rs2", "."),
            (";rs2", "21_34836522_GAGACGGA/-"),
            ("0", "21_34836522_GAGACGGA/-"),
            ("0;rs2", "21_34836522_GAGACGGA/-"),
        ] {
            let variants =
                parse_vcf_line(&format!("21\t34836521\t{id}\tTGAGACGGA\tT\t.\t.\t."), false)
                    .unwrap();
            assert_eq!(variants[0].start, 34836522, "ID {id}");
            assert_eq!(variants[0].allele_string, "GAGACGGA/-", "ID {id}");
            assert_eq!(
                variants[0].uploaded_variation(),
                "21_34836521_TGAGACGGA/T",
                "ID {id}"
            );
            assert_eq!(variants[0].variation_name(), json_id, "ID {id}");
        }
        let variants = parse_vcf_line("21\t34836521\trs2;.\tTGAGACGGA\tT\t.\t.\t.", false).unwrap();
        assert_eq!(variants[0].uploaded_variation(), "rs2");
        assert_eq!(variants[0].variation_name(), "rs2");
    }

    /// Ensembl VEP 116.2 on the GRCh37-hgvs corpus names `21 34971477 . C A`
    /// `21_34971477_C/A`.
    #[test]
    fn dot_id_names_an_snv_by_its_line() {
        assert_eq!(
            uploaded_variations("21\t34971477\t.\tC\tA\t.\t.\t."),
            ["21_34971477_C/A"]
        );
    }

    /// Ensembl VEP 116.2 names `21 10463009 . T TG` (OutputFactory.pm, the line-name branch)
    /// `21_10463009_T/TG` and `21 30677516 . TCACA T` (GRCh37-hgvs corpus)
    /// `21_30677516_TCACA/T`: the POS as written, not the anchor-trimmed start.
    #[test]
    fn dot_id_names_an_indel_by_its_line() {
        let variants = parse_vcf_line("21\t10463009\t.\tT\tTG\t.\t.\t.", false).unwrap();
        assert_eq!(variants[0].start, 10463010);
        assert_eq!(variants[0].uploaded_variation(), "21_10463009_T/TG");
        let variants = parse_vcf_line("21\t30677516\t.\tTCACA\tT\t.\t.\t.", false).unwrap();
        assert_eq!(variants[0].start, 30677517);
        assert_eq!(variants[0].allele_string, "CACA/-");
        assert_eq!(variants[0].uploaded_variation(), "21_30677516_TCACA/T");
    }

    /// Parser.pm `minimise_alleles` copies the record's line onto the minimised
    /// copy, so `21 100 . ATT AT`, minimised to `T/-` at 102, is named
    /// `21_100_ATT/AT` by OutputFactory.pm.
    #[test]
    fn minimised_copy_is_named_by_its_line() {
        let variants = parse_vcf_line("21\t100\t.\tATT\tAT\t.\t.\t.", false).unwrap();
        assert!(variants[0].minimised);
        assert_eq!(variants[0].start, 102);
        assert_eq!(variants[0].allele_string, "T/-");
        assert_eq!(variants[0].uploaded_variation(), "21_100_ATT/AT");
    }

    /// Ensembl VEP 116.2 on the GRCh37 corpus names both alleles of
    /// `21 9489854 . T TAA,TAC` `21_9489854_T/TAA/TAC`; the ALT column is taken
    /// as written, a `*` allele included (OutputFactory.pm turns its commas to `/`).
    #[test]
    fn dot_id_names_a_multi_allelic_record_by_its_line() {
        assert_eq!(
            uploaded_variations("21\t9489854\t.\tT\tTAA,TAC\t.\t.\t."),
            ["21_9489854_T/TAA/TAC", "21_9489854_T/TAA/TAC"]
        );
        assert_eq!(
            uploaded_variations("21\t100\t.\tA\tG,*\t.\t.\t."),
            ["21_100_A/G/*"]
        );
    }

    /// Ensembl VEP 116.2 on the GRCh37 corpus names `21 33867341 . C <CN2>`
    /// `21_33867341_C/<CN2>`: the REF base and the POS, where the structural
    /// variant itself starts at 33867342.
    #[test]
    fn dot_id_names_a_symbolic_record_by_its_line() {
        let variants =
            parse_vcf_line("21\t33867341\t.\tC\t<CN2>\t.\t.\tEND=34000018", false).unwrap();
        assert_eq!(variants[0].start, 33867342);
        assert_eq!(variants[0].uploaded_variation(), "21_33867341_C/<CN2>");
    }

    /// A `*`-only record is named like every other VCF record: Parser/VCF.pm:347 gives
    /// its feature the ID column as written (`.` for a dot), the `id` of its JSON object
    /// (OutputFactory/JSON.pm:236), and OutputFactory.pm's line branch names it
    /// `CHROM_POS_REF/ALT`; the alleles that build that name are the record's own.
    #[test]
    fn dot_id_names_a_star_only_record_by_its_line() {
        let variants = parse_vcf_line("21\t30425650\t.\tN\t*\t.\t.\t.", false).unwrap();
        assert!(variants[0].no_allele);
        assert_eq!(variants[0].variation_name(), ".");
        assert_eq!(variants[0].uploaded_variation(), "21_30425650_N/*");
        assert_eq!(variants[0].uploaded_allele_string.as_deref(), Some("N/*"));
        let named = parse_vcf_line("21\t30425650\tstar\tN\t*,*\t.\t.\t.", false).unwrap();
        assert_eq!(named[0].variation_name(), "star");
        assert_eq!(named[0].uploaded_allele_string.as_deref(), Some("N/*/*"));
    }

    /// OutputFactory.pm: an ID equal to the `chr_start_alleles` name Parser.pm
    /// `validate_vf` generates for an empty ID column is not an identifier, so the
    /// record is named by its line like an ID-less one.
    #[test]
    fn an_id_equal_to_the_generated_name_yields_the_line_name() {
        assert_eq!(
            uploaded_variations("21\t30677516\t21_30677517_CACA/-\tTCACA\tT\t.\t.\t."),
            ["21_30677516_TCACA/T"]
        );
        assert_eq!(
            uploaded_variations("21\t33867341\t21_33867342_<CN2>\tC\t<CN2>\t.\t.\tEND=34000018"),
            ["21_33867341_C/<CN2>"]
        );
        // An ID that merely resembles the generated name is kept.
        assert_eq!(
            uploaded_variations("21\t30677516\t21_30677516_CACA/-\tTCACA\tT\t.\t.\t."),
            ["21_30677516_CACA/-"]
        );
    }

    /// A bracket breakend record's generated name is built from the one allele
    /// string Parser/VCF.pm `create_StructuralVariationFeatures` gives the record,
    /// REF joined to the bracket ALT (`N/N[21:9412900[`) at the start POS + 1, so
    /// an ID equal to it yields the line name on the record's own breakend row and
    /// on the mate row alike (OutputFactory.pm `VariationFeature_to_output_hash`
    /// names the record once); a native single breakend keeps its bare `.N`.
    #[test]
    fn breakend_id_equal_to_the_generated_name_yields_the_line_name_on_both_rows() {
        let variants = parse_vcf_line(
            "21\t9412500\t21_9412501_N/N[21:9412900[\tN\tN[21:9412900[\t.\t.\tSVTYPE=BND",
            false,
        )
        .unwrap();
        assert_eq!(variants.len(), 2);
        assert!(variants[0].is_single_breakend);
        for v in &variants {
            assert_eq!(v.record_allele_string(), "N/N[21:9412900[");
            assert_eq!(v.uploaded_variation(), "21_9412500_N/N[21:9412900[");
        }
        assert_eq!(
            uploaded_variations("21\t9412700\t21_9412701_.N\tN\t.N\t.\t.\tSVTYPE=BND"),
            ["21_9412700_N/.N"]
        );
        // An ID that merely resembles the generated name (no REF) is kept on both rows.
        assert_eq!(
            uploaded_variations(
                "21\t9412500\t21_9412501_N[21:9412900[\tN\tN[21:9412900[\t.\t.\tSVTYPE=BND"
            ),
            ["21_9412501_N[21:9412900[", "21_9412501_N[21:9412900["]
        );
    }

    /// Without MATEID or PARID a mate's `mate_id` stands in the record's identifier
    /// (`id_names_record`): the first `;`-joined ID when it names the record, else the
    /// breakend's `chr:start`, so `rs1;rs2` gives `rs1` and `.`, `;rs2` and `0` give
    /// `21:1001` on a bracket mate and on a CHR2 mate alike.
    #[test]
    fn the_mate_stand_in_is_the_id_only_when_it_names_the_record() {
        let mate_id = |id: &str, alt: &str, info: &str| {
            let variants =
                parse_vcf_line(&format!("21\t1000\t{id}\tN\t{alt}\t.\t.\t{info}"), false).unwrap();
            assert_eq!(variants.len(), 2, "ID {id}");
            variants[1].mate_id.clone()
        };
        assert_eq!(
            mate_id("rs1;rs2", "N[21:1700[", "SVTYPE=BND").as_deref(),
            Some("rs1")
        );
        for id in [".", ";rs2", "0"] {
            assert_eq!(
                mate_id(id, "N[21:1700[", "SVTYPE=BND").as_deref(),
                Some("21:1001"),
                "ID {id}"
            );
            assert_eq!(
                mate_id(id, "<BND>", "SVTYPE=BND;CHR2=21;END2=1700").as_deref(),
                Some("21:1001"),
                "ID {id}"
            );
        }
    }

    #[test]
    fn test_parse_vcf_line_ignores_star_for_joint_trimming() {
        // * (spanning deletion) is filtered: Ensembl VEP does not annotate it.
        let line = "1\t100\t.\tATAAT\tA,*\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].start, 101);
        assert_eq!(variants[0].end, 104);
        assert_eq!(variants[0].ref_allele, b"TAAT");
        assert_eq!(variants[0].alt_allele(), b"-");
    }

    #[test]
    fn test_parse_vcf_line_no_alt() {
        let line = "21\t25000100\t.\tA\t.\t50\tPASS\t.";
        assert!(parse_vcf_line(line, true).is_err());
    }

    #[test]
    fn test_parse_vcf_line_deletion() {
        let line = "21\t25000100\t.\tACGT\tA\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].start, 25000101);
        assert_eq!(variants[0].end, 25000103);
        assert_eq!(variants[0].ref_allele, b"CGT");
        assert_eq!(variants[0].alt_allele(), b"-");
    }

    /// A sequence ALT beside a symbolic one: VEP 116.2 reads the record's type from
    /// INFO/SVTYPE (the joined ALTs open with a sequence allele) and writes one row
    /// per ALT, each with the class term as its Allele (the `synth_multi_mix` records
    /// of `tests/golden/116/GRCh37`, `deletion` twice per transcript).
    #[test]
    fn test_parse_vcf_line_mixed_symbolic_uses_svtype_override() {
        let line = "21\t1000\t.\tC\tG,<INV>\t.\t.\tSVTYPE=DEL;END=2000;SVLEN=-1000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "one variant per ALT, the sequence ALT included"
        );
        for (variant, idx, alt) in [
            (&variants[0], 0, b"G" as &[u8]),
            (&variants[1], 1, b"<INV>"),
        ] {
            assert_eq!(variant.allele_index, idx);
            assert_eq!(variant.alt_allele(), alt);
            assert_eq!(variant.variant_class, VariantClass::StructuralDeletion);
            assert_eq!(variant.display_allele(), "deletion");
            assert!(variant.is_structural);
            assert!(!variant.vep_skip);
            assert_eq!(variant.start, 1001);
            assert_eq!(variant.sv_end, Some(2000));
            assert_eq!(variant.uploaded_allele_string.as_deref(), Some("G/<INV>"));
        }
    }

    /// A record whose REF and ALTs are a run of `ACGT` is never a structural variant,
    /// whatever INFO/SVTYPE says (Parser/VCF.pm:238, `join("", $ref, @$alts) !~
    /// /^[ACGT]+$/`).
    #[test]
    fn all_sequence_multi_allelic_record_ignores_svtype() {
        let line = "21\t1000\t.\tC\tG,T\t.\t.\tSVTYPE=DEL;END=2000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        for v in &variants {
            assert_eq!(v.variant_class, VariantClass::Snv);
            assert!(!v.is_structural);
        }
    }

    /// A mixed record without INFO/SVTYPE has no type VEP knows (`G/<INV>` is in no
    /// SO table), so the symbolic ALT alone is parsed, with the `vep_skip` mark.
    #[test]
    fn test_parse_vcf_line_mixed_symbolic_without_svtype_is_unsupported() {
        let line = "21\t1000\t.\tC\tG,<INV>\t.\t.\tEND=2000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].alt_allele(), b"<INV>");
        assert_eq!(variants[0].variant_class, VariantClass::Inversion);
        assert!(variants[0].vep_skip);
    }

    #[test]
    fn test_parse_vcf_line_paired_bnd_dual_emission() {
        let line = "13\t500\tbnd_1\tA\tA[21:100[\t.\t.\tSVTYPE=BND;MATEID=bnd_2";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert_eq!(paired.alt_allele(), b"A[21:100[");
        assert_eq!(paired.start, 501);
        assert_eq!(paired.end, 501);
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        assert_eq!(paired.mate_pos, Some(100));
        assert!(!paired.is_single_breakend);

        let single = &variants[0];
        assert_eq!(single.variant_class, VariantClass::Translocation);
        assert_eq!(single.alt_allele(), b"A.");
        assert!(single.is_single_breakend);
        assert!(single.mate_chr.is_none());
        assert!(single.mate_pos.is_none());
    }

    #[test]
    fn test_parse_vcf_line_native_single_breakend_not_duplicated() {
        let line = "13\t500\tbnd_1\tA\tA.\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);

        let variant = &variants[0];
        assert_eq!(variant.variant_class, VariantClass::Translocation);
        assert_eq!(variant.alt_allele(), b"A.");
        assert!(variant.is_single_breakend);
        assert!(variant.mate_chr.is_none());
        assert!(variant.mate_pos.is_none());
        // One allele, no bracket mate: VEP prints the class term as its Allele.
        assert_eq!(variant.display_allele(), "chromosome_breakpoint");
    }

    /// INFO/MATEID on a record whose only ALT is a single breakend names no mate:
    /// Ensembl VEP never reads MATEID, `_parse_breakends`
    /// (`StructuralVariationFeature.pm:1279-1283`) finds no bracket in `.N`, and the
    /// row's Allele falls to `class_SO_term` (`OutputFactory.pm:2029`),
    /// `chromosome_breakpoint`, as for the record without MATEID.
    #[test]
    fn single_breakend_with_mateid_is_parsed_as_without() {
        let with =
            parse_vcf_line("21\t100\t.\tN\t.N\t.\t.\tSVTYPE=BND;MATEID=bnd_x", false).unwrap();
        let without = parse_vcf_line("21\t100\t.\tN\t.N\t.\t.\tSVTYPE=BND", false).unwrap();
        assert_eq!(format!("{with:?}"), format!("{without:?}"));
        assert_eq!(with.len(), 1);
        assert!(with[0].is_single_breakend);
        assert_eq!(with[0].mate_id, None);
        assert_eq!(with[0].display_allele(), "chromosome_breakpoint");
    }

    /// A bracket breakend record names its own breakend `N.` or `.N` whether or not
    /// INFO carries a MATEID (`StructuralVariationFeature.pm` `_parse_breakends`
    /// builds the string from the bracket alleles alone): Ensembl VEP 116.2 writes
    /// `N.` on both rows of `bnd_genic_with_intergenic_mate` (no MATEID) in the
    /// GRCh38-release116 golden corpus.
    #[test]
    fn bracket_breakend_without_mateid_names_its_own_allele_n_dot() {
        let line =
            "21\t30425650\tbnd_genic_with_intergenic_mate\tN\tN[21:11600000[\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        assert!(variants[0].is_single_breakend);
        assert_eq!(variants[0].display_allele(), "N.");
        assert_eq!(variants[1].display_allele(), "N[21:11600000[");

        let left =
            parse_vcf_line("21\t11500000\t.\tN\t]21:11000000]N\t.\t.\tSVTYPE=BND", true).unwrap();
        assert_eq!(left[0].display_allele(), ".N");
    }

    /// A record listing two bracket ALTs has one own breakend and two mates, three
    /// variants: `_parse_breakends` (`StructuralVariationFeature.pm:1274-1321`) pushes
    /// one breakend per bracket part of `N/N[21:11500000[/N]21:11600000]` and names
    /// the record's own string once, `N.` when every mate is right-placed and `.N`
    /// when any is left-placed (`:1298`, `:1321`), and `StructuralVariationOverlap::new`
    /// annotates `($vf, @$breakends)`, the own breakend first (`:90`).
    #[test]
    fn two_bracket_alts_yield_one_own_breakend_and_two_mates() {
        let line = "21\t11000000\tmb2\tN\tN[21:11500000[,N]21:11600000]\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        let alleles: Vec<String> = variants.iter().map(|v| v.display_allele()).collect();
        assert_eq!(alleles, ["N.", "N[21:11500000[", "N]21:11600000]"]);
        assert!(variants[0].is_single_breakend);
        assert_eq!(variants[0].allele_index, 0);
        assert_eq!(variants[1].mate_pos, Some(11500000));
        assert_eq!(variants[2].mate_pos, Some(11600000));
        assert_eq!(variants[2].allele_index, 1);
        // The record string of every variant is the whole record, REF and both ALTs,
        // as `Parser/VCF.pm` joins the ALTs (`:469`) and prepends REF (`:499`); JSON
        // `allele_string` prints it (`OutputFactory/JSON.pm:241`).
        for v in &variants {
            assert_eq!(v.record_allele_string(), "N/N[21:11500000[/N]21:11600000]");
        }

        let mixed = "21\t11000000\tmb_mixed\tN\tN[21:11500000[,]21:11600000]N\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(mixed, true).unwrap();
        let alleles: Vec<String> = variants.iter().map(|v| v.display_allele()).collect();
        assert_eq!(alleles, [".N", "N[21:11500000[", "]21:11600000]N"]);
    }

    /// A `<BND>` with SVLEN and no CHR2 names no mate: Parser/VCF.pm:491-497 builds
    /// a bracket mate only under `defined $info->{CHR2}`, so the record stays
    /// `N/<BND>` (`:499`) spanning POS+1 to POS+|SVLEN| (BaseVCF4.pm:302-305) with
    /// no breakend string, one allele whose Allele is `chromosome_breakpoint`
    /// (OutputFactory.pm:2029); Ensembl VEP 116.2 writes `21:10870052-10922272
    /// chromosome_breakpoint` for it. INFO/END is deleted on a breakpoint record
    /// (`:487-489`) and never ends the span.
    #[test]
    fn symbolic_bnd_with_svlen_and_no_chr2_is_one_single_breakend_over_the_span() {
        let line = "21\t10870051\tgnomAD_BND\tN\t<BND>\t.\t.\tSVTYPE=BND;END=10870052;SVLEN=52221";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::Translocation);
        assert_eq!((v.start, v.end), (10870052, 10922272));
        assert_eq!(v.sv_end, Some(10922272));
        assert!(v.is_single_breakend);
        assert_eq!(v.mate_id, None);
        assert!(v.mate_chr.is_none() && v.mate_pos.is_none());
        assert_eq!(v.alt_allele(), b"<BND>");
        assert_eq!(v.display_allele(), "chromosome_breakpoint");
        assert_eq!(v.record_allele_string(), "N/<BND>");
    }

    /// A breakpoint record's INFO/END never ends its span (Parser/VCF.pm:487-489
    /// deletes it before the end is read again): `<BND>` with `END=30425750` and a
    /// bracket breakend with the same END are the single base at POS+1, as Ensembl
    /// VEP 116.2 writes `21:30425651` for both; with CHR2 the deleted END is the
    /// mate position (`:493`), the span still the point.
    #[test]
    fn a_breakpoint_record_s_end_never_ends_its_span() {
        let symbolic = parse_vcf_line(
            "21\t30425650\tend_span\tN\t<BND>\t.\t.\tSVTYPE=BND;END=30425750",
            true,
        )
        .unwrap();
        assert_eq!(symbolic.len(), 1);
        assert_eq!((symbolic[0].start, symbolic[0].end), (30425651, 30425651));
        assert_eq!(symbolic[0].display_allele(), "chromosome_breakpoint");
        let bracket = parse_vcf_line(
            "21\t30425650\tbracket_end\tN\tN[21:30425700[\t.\t.\tSVTYPE=BND;END=30425750",
            true,
        )
        .unwrap();
        assert_eq!(bracket.len(), 2);
        for v in &bracket {
            assert_eq!((v.start, v.end), (30425651, 30425651));
        }
        let chr2_end = parse_vcf_line(
            "21\t30425650\tchr2_far_end\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=2;END=40000000",
            true,
        )
        .unwrap();
        assert_eq!(chr2_end.len(), 2);
        assert_eq!((chr2_end[1].start, chr2_end[1].end), (30425651, 30425651));
        assert_eq!(chr2_end[1].mate_chr.as_deref(), Some("2"));
        assert_eq!(chr2_end[1].mate_pos, Some(40000000));
        assert_eq!(chr2_end[1].display_allele(), "N[2:40000000[");
    }

    /// A `<BND>` names a mate through CHR2 with END2, else the deleted INFO/END
    /// (Parser/VCF.pm:493, `$info->{END2} || $incorrect_end`); POS2 is never read
    /// and an END2 that is not a coordinate names nothing, so Ensembl VEP 116.2
    /// writes `chromosome_breakpoint` with no mate row for `CHR2=2;POS2=40000000`
    /// and for `CHR2=21;END2=-7`. Perl's `||` falls to END only for an END2 that is
    /// empty or `0`: `END2=-7` beside `END=30425700` still names no mate, while
    /// `END2=0` reads END. A signed `+30425700` as END2 or END names none, as
    /// `_parse_breakends` reads the position as `[0-9]+`
    /// (`StructuralVariationFeature.pm:1276`). A bare `CHR2` names none (undef to
    /// `defined $info->{CHR2}`, `:491`), and `CHR2=` is defined and empty, so it builds
    /// `N[:30425800[`, which the bracket regex refuses: Ensembl VEP 116.2 writes
    /// `chromosome_breakpoint` for `CHR2=;END2=30425800`.
    #[test]
    fn pos2_or_an_unreadable_end2_names_no_mate() {
        for info in [
            "SVTYPE=BND;CHR2=2;POS2=40000000",
            "SVTYPE=BND;CHR2=21;END2=-7",
            "SVTYPE=BND;CHR2=21;END2=-7;END=30425700",
            "SVTYPE=BND;CHR2=21;END2=abc;END=30425700",
            "SVTYPE=BND;CHR2=21;END2=+30425700",
            "SVTYPE=BND;CHR2=21;END=+30425700",
            "SVTYPE=BND;CHR2=;END2=30425800",
            "SVTYPE=BND;CHR2;END2=30425800",
        ] {
            let variants =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<BND>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(variants.len(), 1, "{info}");
            assert!(variants[0].is_single_breakend, "{info}");
            assert_eq!(variants[0].mate_pos, None, "{info}");
            assert_eq!(
                variants[0].display_allele(),
                "chromosome_breakpoint",
                "{info}"
            );
        }
        let with_end2 = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700",
            true,
        )
        .unwrap();
        assert_eq!(with_end2.len(), 2);
        assert_eq!(with_end2[1].mate_pos, Some(30425700));
        let zero_end2 = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=0;END=30425700",
            true,
        )
        .unwrap();
        assert_eq!(zero_end2.len(), 2);
        assert_eq!(zero_end2[1].mate_pos, Some(30425700));
    }

    /// The mate INFO/CHR2 names is written `N[CHR2:pos[` whatever the reference base
    /// (Parser/VCF.pm:495-496, `"N"` for a lone `<BND>`), and the own breakend takes the
    /// reference base (`_parse_breakends`, the `A` part of `A/N[21:30425700[`): Ensembl
    /// VEP 116.2 writes `A.` and `N[21:30425700[` with ALLELE_NUM 0 and 1 for REF `A`,
    /// allele string `A/N[21:30425700[`, over the SVLEN span when there is one. A
    /// `<BND:X>` keeps its part (`N/<BND:X>/N[21:30425700[`) and `<BND>,<BND>` both of
    /// its (`N/<BND>/<BND>/N[21:30425700[`), each with one own breakend and one mate.
    #[test]
    fn a_synthesised_mate_is_a_literal_n_and_the_own_breakend_the_reference_base() {
        let variants = parse_vcf_line(
            "21\t30425650\tx\tA\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700",
            true,
        )
        .unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].display_allele(), "A.");
        assert_eq!(variants[1].display_allele(), "N[21:30425700[");
        assert_eq!(variants[1].mate_pos, Some(30425700));
        assert_eq!(variants[1].mate_chr.as_deref(), Some("21"));
        for v in &variants {
            assert_eq!(v.record_allele_string(), "A/N[21:30425700[");
            assert_eq!((v.start, v.end), (30425651, 30425651));
        }
        let spanned = parse_vcf_line(
            "21\t30425650\tx\tA\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700;SVLEN=100",
            true,
        )
        .unwrap();
        assert_eq!(spanned.len(), 2);
        assert!(spanned
            .iter()
            .all(|v| (v.start, v.end) == (30425651, 30425750)));
        for (alt, record) in [
            ("<BND:X>", "N/<BND:X>/N[21:30425700["),
            ("<BND>,<BND>", "N/<BND>/<BND>/N[21:30425700["),
        ] {
            let variants = parse_vcf_line(
                &format!("21\t30425650\tx\tN\t{alt}\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700"),
                true,
            )
            .unwrap();
            assert_eq!(variants.len(), 2, "{alt}");
            assert_eq!(variants[0].display_allele(), "N.", "{alt}");
            assert_eq!(variants[1].display_allele(), "N[21:30425700[", "{alt}");
            assert_eq!(variants[0].record_allele_string(), record, "{alt}");
        }
    }

    /// The bracket text keeps INFO's values as written and the record's string carries
    /// the text whether or not a mate is parsed from it: Ensembl VEP 116.2 writes the
    /// Allele `N[21:00030425700[` for `END2=00030425700` (the mate at 30425700), the
    /// allele string `N/N[21:-7[`, `N/N[21:+30425700[` and `N/N[21:abc[` with
    /// `chromosome_breakpoint` rows for an END2 the regex refuses, `N/N[:30425800[` for
    /// `CHR2=` (defined and empty, Perl's split keeps the field), and `N/<BND>` for a bare
    /// `CHR2` flag or an END2 written empty with no END.
    #[test]
    fn a_mate_text_follows_the_info_values_as_written() {
        let padded = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=00030425700",
            true,
        )
        .unwrap();
        assert_eq!(padded.len(), 2);
        assert_eq!(padded[1].display_allele(), "N[21:00030425700[");
        assert_eq!(padded[1].mate_pos, Some(30425700));
        assert_eq!(padded[0].record_allele_string(), "N/N[21:00030425700[");
        for (info, record) in [
            ("SVTYPE=BND;CHR2=21;END2=-7", "N/N[21:-7["),
            ("SVTYPE=BND;CHR2=21;END2=+30425700", "N/N[21:+30425700["),
            ("SVTYPE=BND;CHR2=21;END2=abc", "N/N[21:abc["),
            ("SVTYPE=BND;CHR2=;END2=30425800", "N/N[:30425800["),
            ("SVTYPE=BND;CHR2;END2=30425800", "N/<BND>"),
            ("SVTYPE=BND;CHR2=21;END2=", "N/<BND>"),
        ] {
            let variants =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<BND>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(variants.len(), 1, "{info}");
            assert_eq!(
                variants[0].display_allele(),
                "chromosome_breakpoint",
                "{info}"
            );
            assert_eq!(variants[0].record_allele_string(), record, "{info}");
            assert!(
                variants[0].is_single_breakend && variants[0].mate_id.is_none(),
                "{info}"
            );
        }
    }

    /// A dot-bounded ALT beside a bracket ALT lends the own breakend its base and gets no
    /// allele of its own (`_parse_breakends`, StructuralVariationFeature.pm:1280-1283
    /// against :1322): Ensembl VEP 116.2 writes `N.` and `N[21:30425700[` for
    /// `.N,N[21:30425700[` with REF `N` or `A`, and `.N` and `]21:30425700]N` for
    /// `N.,]21:30425700]N` (left-placed mate), whose string is `N/N./]21:30425700]N`.
    #[test]
    fn a_dot_alt_beside_a_bracket_alt_lends_the_own_breakend_its_base() {
        for reference in ["N", "A"] {
            let variants = parse_vcf_line(
                &format!("21\t30425650\tx\t{reference}\t.N,N[21:30425700[\t.\t.\tSVTYPE=BND"),
                true,
            )
            .unwrap();
            assert_eq!(variants.len(), 2, "{reference}");
            assert_eq!(variants[0].display_allele(), "N.", "{reference}");
            assert_eq!(
                variants[1].display_allele(),
                "N[21:30425700[",
                "{reference}"
            );
            assert_eq!(
                variants[0].record_allele_string(),
                ".N/N[21:30425700[",
                "{reference}"
            );
        }
        let left = parse_vcf_line(
            "21\t30425650\tx\tN\tN.,]21:30425700]N\t.\t.\tSVTYPE=BND",
            true,
        )
        .unwrap();
        assert_eq!(left.len(), 2);
        assert_eq!(left[0].display_allele(), ".N");
        assert_eq!(left[1].display_allele(), "]21:30425700]N");
        assert_eq!(left[0].record_allele_string(), "N/N./]21:30425700]N");
        assert_eq!(left[1].mate_pos, Some(30425700));
    }

    /// SVLEN is read as Perl's `abs()` numifies it (BaseVCF4.pm:304-305): Ensembl VEP
    /// 116.2 writes the pair `21:30425650-30425651` for `<DUP>`, `<INV>`, `<CNV>`, `<INS>`,
    /// `<BND>` and `<INS:ME:ALU>` with `SVLEN=0`, `SVLEN=abc` or `SVLEN=` (an END beside
    /// them unread), and `21:30425651-30425950` for `SVLEN=300abc`, `SVLEN= 300` and
    /// `SVLEN=3e2`. A `<DEL>` with the pair is a record Ensembl VEP skips ("deletion looks
    /// incomplete"), parsed here as the pair.
    #[test]
    fn svlen_is_read_as_perl_numifies_it() {
        for alt in [
            "<DUP>",
            "<INV>",
            "<CNV>",
            "<INS>",
            "<BND>",
            "<INS:ME:ALU>",
            "<DEL>",
        ] {
            for svlen in ["0", "abc", ""] {
                let variants = parse_vcf_line(
                    &format!("21\t30425650\tx\tN\t{alt}\t.\t.\tSVLEN={svlen};END=30425899"),
                    true,
                )
                .unwrap();
                assert_eq!(variants.len(), 1, "{alt} {svlen:?}");
                assert_eq!(
                    (variants[0].start, variants[0].end),
                    (30425651, 30425650),
                    "{alt} {svlen:?}"
                );
                assert_eq!(variants[0].sv_len, Some(0), "{alt} {svlen:?}");
            }
        }
        for svlen in ["300abc", " 300", "3e2", "+300", "-300"] {
            let variants = parse_vcf_line(
                &format!("21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;SVLEN={svlen}"),
                true,
            )
            .unwrap();
            assert_eq!(
                (variants[0].start, variants[0].end),
                (30425651, 30425950),
                "{svlen:?}"
            );
        }
    }

    /// `validate_vf` (Parser.pm:483-488) skips a record whose end lies more than one base
    /// before its start, and Ensembl VEP 116.2 does so for `END=0` and `END=30425000` at
    /// POS 30425650 for every class, a `<BND>` included, whose deleted END is re-read only
    /// when Perl reads it as true (Parser/VCF.pm:489), so `END=0` keeps the 0 and `END=`
    /// the empty string ("Invalid start or end coordinate"). `END=30425650` is the pair
    /// and `END=30425651` the single base at POS + 1.
    #[test]
    fn an_end_before_pos_is_refused_and_an_end_at_pos_is_the_pair() {
        for (alt, info) in [
            ("<DUP>", "SVTYPE=DUP;END=0"),
            ("<DUP>", "SVTYPE=DUP;END=30425000"),
            ("<DEL>", "SVTYPE=DEL;END=0"),
            ("<BND>", "SVTYPE=BND;CHR2=21;END=0"),
            ("<BND>", "SVTYPE=BND;END=0"),
        ] {
            let err = parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("start > end+1"), "{alt} {info}: {err}");
        }
        for (alt, info) in [("<BND>", "SVTYPE=BND;END="), ("<DUP>", "SVTYPE=DUP;END=")] {
            let err = parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("Invalid END: ''"), "{alt} {info}: {err}");
        }
        let pair = parse_vcf_line(
            "21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=30425650",
            true,
        )
        .unwrap();
        assert_eq!((pair[0].start, pair[0].end), (30425651, 30425650));
        let point = parse_vcf_line(
            "21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=30425651",
            true,
        )
        .unwrap();
        assert_eq!((point[0].start, point[0].end), (30425651, 30425651));
        let kept_bnd = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;SVLEN=100;END=0",
            true,
        )
        .unwrap();
        assert_eq!((kept_bnd[0].start, kept_bnd[0].end), (30425651, 30425750));
    }

    /// A repeated INFO key reads its last value, as the hash Ensembl VEP's parser builds
    /// does (BaseVCF4.pm:581-583): Ensembl VEP 116.2 writes `21:30425651-30425950` for
    /// `SVLEN=-100;SVLEN=-300`, `21:30425651-30425900` for `END=30425700;END=30425900`,
    /// the mate `N[21:30425800[` for `END2=30425700;END2=30425800` and for
    /// `CHR2=2;CHR2=21`, and `duplication` for `N.` with `SVTYPE=DEL;SVTYPE=DUP`. A bare
    /// key written last is undef, so `CHR2=21;CHR2` names no mate.
    #[test]
    fn a_repeated_info_key_reads_its_last_value() {
        let svlen = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-100;SVLEN=-300",
            true,
        )
        .unwrap();
        assert_eq!((svlen[0].start, svlen[0].end), (30425651, 30425950));
        let end = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425700;END=30425900",
            true,
        )
        .unwrap();
        assert_eq!(end[0].end, 30425900);
        for info in [
            "SVTYPE=BND;CHR2=21;END2=30425700;END2=30425800",
            "SVTYPE=BND;CHR2=2;CHR2=21;END2=30425800",
        ] {
            let variants =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<BND>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(variants.len(), 2, "{info}");
            assert_eq!(variants[1].display_allele(), "N[21:30425800[", "{info}");
            assert_eq!(variants[1].mate_chr.as_deref(), Some("21"), "{info}");
        }
        let svtype = parse_vcf_line(
            "21\t30425650\tx\tN\tN.\t.\t.\tSVTYPE=DEL;SVTYPE=DUP;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(svtype[0].variant_class, VariantClass::Duplication);
        let bare_last = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;CHR2;END2=30425800",
            true,
        )
        .unwrap();
        assert_eq!(bare_last.len(), 1);
        assert_eq!(bare_last[0].record_allele_string(), "N/<BND>");
    }

    /// A signed coordinate is a number to Perl (`looks_like_number`, Parser.pm:475;
    /// `$start++` on the raw POS): Ensembl VEP 116.2 annotates `END=+30425899` as the end
    /// 30425899 and `SVLEN=+300` as 300, as a POS `+30425650` is 30425650; the sign it
    /// prints inside the Location (`21:30425651-+30425899`) is not reproduced.
    #[test]
    fn a_signed_coordinate_reads_as_its_number() {
        let end = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=+30425899",
            true,
        )
        .unwrap();
        assert_eq!((end[0].start, end[0].end), (30425651, 30425899));
        let pos = parse_vcf_line(
            "21\t+30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!((pos[0].start, pos[0].end), (30425651, 30425899));
        let svlen = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=+300",
            true,
        )
        .unwrap();
        assert_eq!(svlen[0].end, 30425950);
    }

    /// A record with a parsed breakend is one breakpoint whatever its other ALTs: Ensembl
    /// VEP 116.2 writes `N.` and `N[21:30425700[` at the point `21:30425651` for
    /// `<DEL>,N[21:30425700[` with `SVTYPE=BND;END=30425899`, the `<DEL>` part lending
    /// nothing (`_parse_breakends` reads no base from it) and the END unread.
    #[test]
    fn a_record_with_a_parsed_breakend_is_one_breakpoint_whatever_its_other_alts() {
        let mixed = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>,N[21:30425700[\t.\t.\tSVTYPE=BND;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(mixed.len(), 2);
        for v in &mixed {
            assert_eq!(v.variant_class, VariantClass::Translocation);
            assert_eq!((v.start, v.end), (30425651, 30425651));
            assert_eq!(v.record_allele_string(), "N/<DEL>/N[21:30425700[");
        }
        assert_eq!(mixed[0].display_allele(), "N.");
        assert_eq!(mixed[1].display_allele(), "N[21:30425700[");
    }

    /// A bare `SVLEN` or `END` flag is undef to Perl (`get_end` tests both with `defined`,
    /// BaseVCF4.pm:302-307), so `<DUP>` with a bare `SVLEN` beside `END=30425899` ends at
    /// 30425899 and one with a bare `END` alone is the single base at POS + 1.
    #[test]
    fn a_bare_svlen_or_end_flag_is_undef() {
        let svlen_flag = parse_vcf_line(
            "21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;SVLEN;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(
            (svlen_flag[0].start, svlen_flag[0].end),
            (30425651, 30425899)
        );
        assert_eq!(svlen_flag[0].sv_len, None);
        let end_flag =
            parse_vcf_line("21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END", true).unwrap();
        assert_eq!((end_flag[0].start, end_flag[0].end), (30425651, 30425651));
    }

    /// Parser/VCF.pm:237-245 routes a line to the structural path when REF and the ALTs
    /// together are not a run of `ACGT` and INFO/SVTYPE is true or an ALT is written in
    /// a structural form: Ensembl VEP 116.2 writes `21:30425651-30425949 deletion` for
    /// `N G` and for `<*>` with `SVTYPE=DEL;END=30425949`, and for `A g` (a lowercase
    /// base is outside `ACGT`), while `A G` with the same INFO is the SNV
    /// `21:30425650 G missense_variant` and `N G` alone the SNV
    /// `coding_sequence_variant`. `SVTYPE=` and `SVTYPE=0` are false to Perl and route
    /// nothing; `<*>` alone is a sequence allele at POS.
    #[test]
    fn svtype_routes_a_record_whose_alleles_are_not_plain_acgt() {
        for (reference, alt) in [("N", "G"), ("N", "<*>"), ("A", "g")] {
            let v = parse_vcf_line(
                &format!("21\t30425650\tx\t{reference}\t{alt}\t.\t.\tSVTYPE=DEL;END=30425949"),
                true,
            )
            .unwrap();
            assert_eq!(v.len(), 1, "{reference} {alt}");
            assert_eq!(
                v[0].variant_class,
                VariantClass::StructuralDeletion,
                "{reference} {alt}"
            );
            assert_eq!(
                (v[0].start, v[0].end),
                (30425651, 30425949),
                "{reference} {alt}"
            );
            assert_eq!(v[0].display_allele(), "deletion", "{reference} {alt}");
        }
        for (reference, alt, info) in [
            ("A", "G", "SVTYPE=DEL;END=30425949"),
            ("N", "G", "."),
            ("N", "G", "SVTYPE=;END=30425949"),
            ("N", "G", "SVTYPE=0;END=30425949"),
            ("N", "<*>", "."),
        ] {
            let v = parse_vcf_line(
                &format!("21\t30425650\tx\t{reference}\t{alt}\t.\t.\t{info}"),
                true,
            )
            .unwrap();
            assert_eq!(v.len(), 1, "{reference} {alt} {info}");
            assert!(!v[0].is_structural, "{reference} {alt} {info}");
            assert_eq!(
                (v[0].start, v[0].end),
                (30425650, 30425650),
                "{reference} {alt} {info}"
            );
            assert_eq!(v[0].display_allele(), alt, "{reference} {alt} {info}");
        }
    }

    /// Parser/VCF.pm:473 tests the VCF 4.4 forms without `/i`, so a lowercase ALT is
    /// outside them and SVTYPE names its type; `get_SO_term`'s fall-through branches
    /// (Parser.pm:712-717) keep the type as written and `%SO_TERMS` has upper-case keys
    /// only. Ensembl VEP 116.2 writes `insertion` for `<ins:me:alu>` with `SVTYPE=INS`,
    /// `duplication` for `<del>` with `SVTYPE=DUP`, `copy_number_variation` for `<cnv:tr>`
    /// with `SVTYPE=CNV` and for `<DEL>,<dup>` with `SVTYPE=DEL`, and skips `<del>` alone
    /// and `<BND>` with `SVTYPE=bnd` as unsupported types.
    #[test]
    fn a_lowercase_symbolic_alt_is_outside_the_vcf44_forms() {
        for (alt, info, term, class) in [
            (
                "<ins:me:alu>",
                "SVTYPE=INS",
                "insertion",
                VariantClass::StructuralInsertion,
            ),
            (
                "<del>",
                "SVTYPE=DUP;END=30425949",
                "duplication",
                VariantClass::Duplication,
            ),
            (
                "<cnv:tr>",
                "SVTYPE=CNV;END=30425949",
                "copy_number_variation",
                VariantClass::CopyNumberVariation,
            ),
            (
                "<DEL>,<dup>",
                "SVTYPE=DEL;END=30425949",
                "copy_number_variation",
                VariantClass::CopyNumberVariation,
            ),
            (
                "<DEL>,<del>",
                "SVTYPE=DEL;END=30425949",
                "deletion",
                VariantClass::StructuralDeletion,
            ),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            for variant in &v {
                assert_eq!(
                    variant.record_so_term.as_deref(),
                    Some(term),
                    "{alt} {info}"
                );
                assert_eq!(variant.variant_class, class, "{alt} {info}");
                assert!(!variant.vep_skip, "{alt} {info}");
            }
        }
        for (alt, info) in [
            ("<del>", "."),
            ("<BND>", "SVTYPE=bnd"),
            ("<del>", "SVTYPE=del;END=30425949"),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            assert_eq!(v.len(), 1, "{alt} {info}");
            assert!(v[0].vep_skip, "{alt} {info}");
            assert_eq!(v[0].record_so_term, None, "{alt} {info}");
        }
        assert_eq!(sv_type_abbreviation("<del>", None), "del");
        assert_eq!(sv_type_abbreviation("<BND>", Some("bnd")), "bnd");
        assert_eq!(
            sv_type_abbreviation("<cn0>", None),
            "CNV",
            "`/^<?CN=?0>?$/` is upper case only"
        );
        assert_eq!(sv_type_abbreviation("<CN0>", None), "DEL");
    }

    /// A gVCF variant site lists `<NON_REF>` beside its real alleles; Ensembl VEP 116.2
    /// drops the line (`G/NON_REF is not a supported structural variant type`), and
    /// vep-rs reads the real alleles, numbered as they stand: `N G,<NON_REF>` is the SNV
    /// `G` at POS. A record whose only ALT is `<NON_REF>` is a reference block.
    #[test]
    fn a_gvcf_variant_site_is_its_real_alleles() {
        let site = parse_vcf_line("21\t30425650\tx\tA\tG,<NON_REF>\t.\t.\t.", true).unwrap();
        assert_eq!(site.len(), 1);
        assert_eq!(
            (
                site[0].display_allele(),
                site[0].allele_index,
                site[0].is_structural
            ),
            ("G".to_string(), 0, false)
        );
        let indel = parse_vcf_line("21\t30425650\tx\tN\tNA,<NON_REF>\t.\t.\t.", true).unwrap();
        assert_eq!(
            (indel.len(), indel[0].display_allele()),
            (1, "A".to_string())
        );
        let block =
            parse_vcf_line("21\t30425650\tx\tN\t<NON_REF>\t.\t.\tEND=30425700", true).unwrap();
        assert_eq!(
            (block.len(), block[0].variant_class),
            (1, VariantClass::ReferenceBlock)
        );
    }

    /// One StructuralVariationFeature carries the record (Parser/VCF.pm:560-574), so a
    /// record whose class is the breakpoint's gives every ALT the breakpoint's class and
    /// point: Ensembl VEP 116.2 writes `21:30425651 chromosome_breakpoint` for each part of
    /// `<BND>,<DEL>` with `SVTYPE=BND;END=30425949`, the END unread, with the allele string
    /// `N/<BND>/<DEL>`. A `*` allele is a part of the string (`*/<DEL>`) and a variant of
    /// the record's class, numbered at its written place (StructuralVariationOverlap.pm:64-77
    /// makes one allele per part): 116.2 writes two `deletion` rows per transcript for
    /// `N *,<DEL>` with `SVTYPE=DEL;END=30425949`. A lone `.` ALT names none whatever
    /// SVTYPE says.
    #[test]
    fn a_breakpoint_record_puts_every_alt_at_the_point_and_keeps_star_in_its_string() {
        let v = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>,<DEL>\t.\t.\tSVTYPE=BND;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(v.len(), 2);
        for variant in &v {
            assert_eq!(variant.variant_class, VariantClass::Translocation);
            assert_eq!((variant.start, variant.end), (30425651, 30425651));
            assert_eq!(variant.display_allele(), "chromosome_breakpoint");
            assert_eq!(variant.record_allele_string(), "N/<BND>/<DEL>");
        }
        let star = parse_vcf_line(
            "21\t30425650\tx\tN\t*,<DEL>\t.\t.\tSVTYPE=DEL;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(star.len(), 2);
        for (variant, (index, alt)) in star.iter().zip([(0, b"*".as_slice()), (1, b"<DEL>")]) {
            assert_eq!(variant.record_allele_string(), "*/<DEL>");
            assert_eq!((variant.allele_index, variant.alt_allele()), (index, alt));
            assert_eq!(variant.variant_class, VariantClass::StructuralDeletion);
            assert_eq!((variant.start, variant.end), (30425651, 30425949));
            assert_eq!(variant.display_allele(), "deletion");
        }
        for info in ["SVTYPE=BND", "SVTYPE=DEL;END=30425949", "."] {
            let err = parse_vcf_line(&format!("21\t30425650\tx\tN\t.\t.\t.\t{info}"), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("No alternate allele"), "{info}: {err}");
        }
        let err = parse_vcf_line("21\t30425650\tx\tN\tG\t.\t.\tSVTYPE=FOO;END=30425949", true)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("FOO is not a supported structural variant type"),
            "{err}"
        );
    }

    /// ALLELE_NUM numbers the parts of Ensembl VEP's allele string from 1
    /// (StructuralVariationOverlap.pm:64-77): a bare `<BND>` is `N/<BND>`, so its row
    /// prints 2 (116.2 writes `ALLELE_NUM=2` beside the reference allele's 1); `.N` is
    /// not prefixed and prints 1; a record with parsed breakends prints 0 for its own
    /// breakend and 1 upward for its mates; every other class prints the ALT's position.
    #[test]
    fn allele_num_follows_ensembl_veps_allele_string() {
        let bare = parse_vcf_line("21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND", true).unwrap();
        assert_eq!(bare[0].allele_num, Some(2));
        let dot = parse_vcf_line("21\t30425650\tx\tN\t.N\t.\t.\tSVTYPE=BND", true).unwrap();
        assert_eq!(dot[0].allele_num, Some(1));
        let paired = parse_vcf_line(
            "21\t30425650\tx\tN\tN[21:30425700[,N]22:100]\t.\t.\tSVTYPE=BND",
            true,
        )
        .unwrap();
        assert_eq!(
            paired.iter().map(|v| v.allele_num).collect::<Vec<_>>(),
            [Some(0), Some(1), Some(2)]
        );
        let del = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(
            del.iter().map(|v| v.allele_num).collect::<Vec<_>>(),
            [None, None]
        );
    }

    /// `get_SO_term`'s regexes decide the type of a shape at their edges: `/CN=?[0-9]/i`
    /// needs a digit, the list normalisation strips one bracket a side, the mobile-element
    /// pattern takes one colon before the element, and `s/\:.+//g` needs a character after
    /// the colon. Ensembl VEP 116.2 skips `<CN=X>`, `<<DEL>>,<DEL>` and `<DEL:>` as
    /// unsupported types (`CN=X`, `DEL/DEL`, `DEL:`), writes `copy_number_variation` for
    /// `<<CN0>>` and `mobile_element_insertion` for `<INS:ME::ALU>` with `SVTYPE=INS`, and
    /// skips `<NON_REF>,<NON_REF>` as `NON_REF`, which vep-rs reads as a reference block
    /// of two alleles, one intergenic row each; under `SVTYPE=DEL` it writes two `deletion`
    /// rows per transcript for `G,<NON_REF>`, the block an allele of the record.
    #[test]
    fn the_type_follows_get_so_terms_regexes_at_their_edges() {
        let parse = |alt: &str, info: &str| {
            parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap()
        };
        for (alt, info, abbrev) in [
            ("<CN=X>", "END=30425949", "CN=X"),
            ("<<DEL>>,<DEL>", "END=30425949", "DEL/DEL"),
            ("<DEL:>", "END=30425949", "DEL:"),
        ] {
            let v = parse(alt, info);
            assert!(v.iter().all(|v| v.vep_skip), "{alt}");
            assert!(v.iter().all(|v| v.record_so_term.is_none()), "{alt}");
            assert_eq!(sv_type_abbreviation(&alt.replace(',', "/"), None), abbrev);
        }
        for (alt, info, term) in [
            ("<CN=X>", "SVTYPE=DEL;END=30425949", "deletion"),
            ("<<DEL>>", "END=30425949", "deletion"),
            ("<<CN0>>", "END=30425949", "copy_number_variation"),
            ("<CN0>", "END=30425949", "deletion"),
            ("<INS:ME::ALU>", "SVTYPE=INS", "mobile_element_insertion"),
            ("<INS:ME:ALU>", "SVTYPE=INS", "Alu_insertion"),
            ("<DEL:X>", "END=30425949", "deletion"),
            ("<CNX>,<CN3>", "END=30425949", "copy_number_variation"),
        ] {
            let v = parse(alt, info);
            assert!(
                !v.is_empty() && v.iter().all(|v| !v.vep_skip),
                "{alt} {info}"
            );
            assert!(
                v.iter().all(|v| v.record_so_term.as_deref() == Some(term)),
                "{alt} {info}: {:?}",
                v[0].record_so_term
            );
        }
        let site = parse("G,<NON_REF>", "SVTYPE=DEL;END=30425949");
        assert_eq!(
            site.len(),
            2,
            "a true SVTYPE keeps the block as a deletion allele"
        );
        for variant in &site {
            assert_eq!(variant.record_so_term.as_deref(), Some("deletion"));
            assert_eq!(variant.record_allele_string(), "G/<NON_REF>");
        }
        let blocks = parse("<NON_REF>,<NON_REF>", "END=30425700");
        assert_eq!(blocks.len(), 2);
        for (variant, index) in blocks.iter().zip(0..) {
            assert_eq!(variant.variant_class, VariantClass::ReferenceBlock);
            assert!(variant.vep_skip);
            assert_eq!(variant.allele_index, index);
            assert_eq!(variant.record_allele_string(), "<NON_REF>/<NON_REF>");
        }
    }

    /// Perl's list-context `split` drops a trailing empty ALT field (BaseVCF4.pm:489) and a
    /// trailing empty type part (Parser.pm:669): Ensembl VEP 116.2 writes the one allele `G`
    /// for `A G,`, a `deletion` over END for `<DEL>,` with `SVTYPE=DEL`, and a `deletion`
    /// for `N.` with `SVTYPE=DEL/`.
    #[test]
    fn a_trailing_empty_alt_field_names_no_allele() {
        let seq = parse_vcf_line("21\t30425650\tx\tA\tG,\t.\t.\t.", true).unwrap();
        assert_eq!(seq.len(), 1);
        assert_eq!(
            (seq[0].display_allele(), seq[0].allele_index),
            ("G".to_string(), 0)
        );
        assert_eq!(seq[0].uploaded_allele_string.as_deref(), Some("A/G"));
        for (alt, info) in [
            ("<DEL>,", "SVTYPE=DEL;END=30425949"),
            ("N.", "SVTYPE=DEL/;END=30425949"),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            assert_eq!(v.len(), 1, "{alt} {info}");
            assert_eq!(
                v[0].record_so_term.as_deref(),
                Some("deletion"),
                "{alt} {info}"
            );
            assert!(!v[0].vep_skip, "{alt} {info}");
            assert_eq!((v[0].start, v[0].end), (30425651, 30425949), "{alt} {info}");
        }
    }

    /// A gVCF site's real allele keeps its written position whichever side of `<NON_REF>`
    /// it stands on, so `<NON_REF>,G` prints ALLELE_NUM 2 for `G`, the number a genotype
    /// `0/2` names.
    #[test]
    fn a_gvcf_sites_allele_keeps_its_written_position() {
        let after = parse_vcf_line("21\t30425650\tx\tA\t<NON_REF>,G\t.\t.\t.", true).unwrap();
        assert_eq!(after.len(), 1);
        assert_eq!(
            (after[0].display_allele(), after[0].allele_index),
            ("G".to_string(), 1)
        );
        let before = parse_vcf_line("21\t30425650\tx\tA\tG,<NON_REF>\t.\t.\t.", true).unwrap();
        assert_eq!(before[0].allele_index, 0);
    }

    /// BaseVCF4.pm:581-583 builds the INFO hash key by key as written (`$info_data{$key} =
    /// $value`), the last occurrence winning: a bare `END` written last leaves END undef, and
    /// `end=` or ` END=` is not `END`, so a deletion carrying only those has no end and is
    /// Ensembl VEP's `deletion looks incomplete` (Parser/VCF.pm:555-558); vep-rs ends such a
    /// record at its start.
    #[test]
    fn an_info_key_matches_exactly_and_a_bare_key_written_last_is_undef() {
        let dup =
            parse_vcf_line("21\t1000\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=1100;END", true).unwrap();
        assert_eq!((dup[0].start, dup[0].end), (1001, 1001));
        for info in [
            "svtype=DEL;end=30425949",
            "SVTYPE=DEL; END=30425949",
            "SVTYPE=DEL;END =30425949",
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<DEL>\t.\t.\t{info}"), true).unwrap();
            assert_eq!((v[0].start, v[0].end), (30425651, 30425651), "{info}");
            assert_eq!(v[0].record_so_term.as_deref(), Some("deletion"), "{info}");
        }
        let exact = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(exact[0].end, 30425949);
    }

    /// `_parse_breakends`' match tolerates whitespace around the brackets, the name and the
    /// colon, and BaseVCF4.pm:582 cuts an INFO value at a second `=`: Ensembl VEP 116.2
    /// writes the own breakend `N.` (ALLELE_NUM 0) and a mate row for `N[ 21 : 30425700 [`,
    /// for `<BND>` with `CHR2= 21;END2=30425700` (Allele `N[ 21:30425700[`) and for `<BND>`
    /// with `CHR2=21;END2=30425700=5` (Allele `N[21:30425700[`).
    #[test]
    fn a_breakend_mate_is_read_through_the_whitespace_perl_allows() {
        for (alt, info, mate_text) in [
            ("N[ 21 : 30425700 [", "SVTYPE=BND", "N[ 21 : 30425700 ["),
            (
                "<BND>",
                "SVTYPE=BND;CHR2= 21;END2=30425700",
                "N[ 21:30425700[",
            ),
            (
                "<BND>",
                "SVTYPE=BND;CHR2=21;END2=30425700=5",
                "N[21:30425700[",
            ),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            assert_eq!(v.len(), 2, "{alt} {info}");
            assert_eq!(
                (v[0].display_allele(), v[0].allele_num),
                ("N.".to_string(), Some(0)),
                "{alt}"
            );
            assert_eq!(
                (v[1].display_allele(), v[1].allele_num),
                (mate_text.to_string(), Some(1)),
                "{alt}"
            );
            assert_eq!(
                (v[1].mate_chr.as_deref(), v[1].mate_pos),
                (Some("21"), Some(30425700)),
                "{alt}"
            );
        }
    }

    /// A unit of no bases under RUC is the empty allele (`'' x 5`, Parser/VCF.pm:419), an
    /// allele string `REF/` that names no alternate allele, so the record is the feature
    /// with no allele to annotate, carrying that string; under RB it divides by zero
    /// (`:414`) and keeps the symbolic reading; a count is read as Perl numifies it
    /// (`RUC=2abc` repeats twice, `TCATCATCA/TCATCA`, the loss of the run's last unit).
    #[test]
    fn a_zero_length_repeat_unit_is_read_as_perl_expands_it() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let record = |info: &str| {
            parse_with_fasta(
                &format!("21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;{info}"),
                Some(&fasta),
                Some(DEFAULT_MAX_SV_SIZE),
            )
        };
        let empty = record("RUS=,AC;RUC=5");
        assert!(empty[0].no_allele, "{empty:?}");
        assert_eq!(empty[0].allele_string, "TCATCATCA/");
        let by_zero = record("RUS=,AC;RB=12");
        assert!(by_zero[0].is_structural);
        assert_eq!(by_zero[0].variant_class, VariantClass::TandemRepeat);
        assert_eq!(
            shapes(&record("RUS=TCA;RUC=2abc")),
            [(
                "21:17-19".to_string(),
                "TCA/-".to_string(),
                VariantClass::Deletion,
                false
            )]
        );
    }

    /// The record-wide strings (the uploaded allele string, the record's working string,
    /// the line name and the raw line) are one buffer per record shared by its alleles.
    #[test]
    fn the_record_strings_are_shared_by_its_alleles() {
        let shared = |variants: &[InputVariant]| {
            let first = &variants[0];
            variants.iter().skip(1).all(|v| {
                Arc::ptr_eq(
                    v.uploaded_allele_string.as_ref().unwrap(),
                    first.uploaded_allele_string.as_ref().unwrap(),
                ) && Arc::ptr_eq(
                    v.input_line_name.as_ref().unwrap(),
                    first.input_line_name.as_ref().unwrap(),
                ) && Arc::ptr_eq(
                    v.raw_input.as_ref().unwrap(),
                    first.raw_input.as_ref().unwrap(),
                ) && match (
                    &v.record_allele_string_multi,
                    &first.record_allele_string_multi,
                ) {
                    (Some(a), Some(b)) => Arc::ptr_eq(a, b),
                    (None, None) => true,
                    _ => false,
                }
            })
        };
        let seq = parse_vcf_line("21\t1000\tx\tA\tG,T,C\t.\t.\t.", true).unwrap();
        assert!(seq.len() == 3 && shared(&seq));
        let sv = parse_vcf_line(
            "21\t1000\tx\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;END=1100",
            true,
        )
        .unwrap();
        assert!(sv.len() == 2 && shared(&sv));
        let bnd = parse_vcf_line(
            "21\t1000\tx\tN\tN[21:1700[,N]22:100]\t.\t.\tSVTYPE=BND",
            true,
        )
        .unwrap();
        assert!(bnd.len() == 3 && shared(&bnd));
    }

    /// A comma-separated SVLEN is read as its first value for the span
    /// (BaseVCF4.pm:304, `(split(',',$info->{SVLEN}))[0]`, so `-300,400` ends at
    /// POS+300) and as its largest for a tandem repeat (Parser/VCF.pm:505-511,
    /// `$end = $start + max(@svlen) - 1`, so `300,400` ends at POS+400), as Ensembl
    /// VEP 116.2 writes `21:30425601-30425900` and `21:30425601-30426000`.
    #[test]
    fn a_comma_separated_svlen_reads_its_first_value_and_a_tandem_repeat_its_largest() {
        let cnv = parse_vcf_line(
            "21\t30425600\tx\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;SVLEN=-300,400",
            true,
        )
        .unwrap();
        assert_eq!(cnv.len(), 2);
        for v in &cnv {
            assert_eq!((v.start, v.end), (30425601, 30425900));
            assert_eq!(v.sv_len, Some(-300));
        }
        let repeat = parse_vcf_line(
            "21\t30425600\tx\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;SVLEN=300,400;RUS=AT,AT;RUC=100,50",
            true,
        )
        .unwrap();
        assert_eq!(repeat.len(), 2);
        for v in &repeat {
            assert_eq!((v.start, v.end), (30425601, 30426000));
        }
    }

    /// A tandem duplication reads SVLEN as a tandem repeat does: Parser/VCF.pm:502
    /// routes every record whose SO term matches `/tandem/` (`tandem_duplication`
    /// included) to `$end = $start + max(@svlen) - 1` (`:505-511`), so Ensembl VEP
    /// 116.2 writes `21:30425601-30426000 tandem_duplication` for `<DUP:TANDEM>` with
    /// `SVLEN=300,400` and for `<DUP:TANDEM>,<DUP:TANDEM>` with the same SVLEN, and
    /// skips `SVLEN=-300` (`start > end+1`) with or without END, while a `<DUP>` with
    /// `SVLEN=300,400` or `SVLEN=-300` is `21:30425601-30425900 duplication`.
    #[test]
    fn a_tandem_duplication_reads_svlen_as_a_tandem_repeat_does() {
        let tandem = parse_vcf_line(
            "21\t30425600\tx\tN\t<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;SVLEN=300,400",
            true,
        )
        .unwrap();
        assert_eq!(tandem.len(), 1);
        assert_eq!(tandem[0].variant_class, VariantClass::TandemDuplication);
        assert_eq!((tandem[0].start, tandem[0].end), (30425601, 30426000));
        assert_eq!(tandem[0].sv_len, Some(300));
        let two = parse_vcf_line(
            "21\t30425600\tx\tN\t<DUP:TANDEM>,<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;SVLEN=300,400",
            true,
        )
        .unwrap();
        assert_eq!(two.len(), 2);
        for v in &two {
            assert_eq!(v.variant_class, VariantClass::TandemDuplication);
            assert_eq!((v.start, v.end), (30425601, 30426000));
        }
        for info in [
            "SVTYPE=DUP;SVLEN=-300",
            "SVTYPE=DUP;SVLEN=-300;END=30425900",
        ] {
            let err = parse_vcf_line(
                &format!("21\t30425600\tx\tN\t<DUP:TANDEM>\t.\t.\t{info}"),
                true,
            )
            .unwrap_err()
            .to_string();
            assert!(
                err.contains("start > end+1 : (START=30425601, END=30425300)"),
                "{info}: {err}"
            );
        }
        for info in ["SVTYPE=DUP;SVLEN=300,400", "SVTYPE=DUP;SVLEN=-300"] {
            let dup =
                parse_vcf_line(&format!("21\t30425600\tx\tN\t<DUP>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(dup[0].variant_class, VariantClass::Duplication, "{info}");
            assert_eq!((dup[0].start, dup[0].end), (30425601, 30425900), "{info}");
        }
    }

    /// `N.` with `SVTYPE=DEL` is a deletion: the ALT is outside the VCF 4.4 forms, so
    /// INFO/SVTYPE is the type (Parser/VCF.pm:472-475) and `get_SO_term` reads
    /// `deletion`; Ensembl VEP 116.2 writes `21:30425601-30425899 deletion` with
    /// `frameshift_variant,feature_truncation`. A `<BND>` without SVTYPE is a
    /// breakpoint from its own ALT.
    #[test]
    fn svtype_names_the_class_of_a_dot_alt_and_a_bare_bnd_is_a_breakpoint() {
        let deletion = parse_vcf_line(
            "21\t30425600\tgnomad_del\tN\tN.\t.\t.\tSVTYPE=DEL;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(deletion.len(), 1);
        assert_eq!(deletion[0].variant_class, VariantClass::StructuralDeletion);
        assert_eq!((deletion[0].start, deletion[0].end), (30425601, 30425899));
        assert!(!deletion[0].is_single_breakend);
        assert_eq!(deletion[0].display_allele(), "deletion");
        let breakpoint =
            parse_vcf_line("21\t30425650\tnosvtype\tN\t<BND>\t.\t.\tEND=30425651", true).unwrap();
        assert_eq!(breakpoint[0].variant_class, VariantClass::Translocation);
        assert_eq!(breakpoint[0].display_allele(), "chromosome_breakpoint");
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_without_svlen_single_variant() {
        // <BND> without SVLEN or CHR2: the single base at POS+1, naming no mate.
        let line = "21\t100\ttest_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND;END=101";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            1,
            "symbolic <BND> without SVLEN should emit 1 variant"
        );
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::Translocation);
        assert_eq!(v.start, 101); // POS+1
        assert_eq!(v.end, 101); // INFO/END never ends a breakpoint's span
    }

    /// A coordinate is read up to `MAX_COORD`, 2^62 - 1: a POS past it fails the line,
    /// as does a negative POS, `-0` included; an END past it (`i64::MAX` included), a
    /// negative END or text fails the line of every class but a breakpoint, as Ensembl
    /// VEP skips the record (Parser.pm `validate_vf`: `start > end+1`, the size test,
    /// "Invalid start or end coordinate"), while a breakpoint record ignores END for
    /// its span and is the single base at POS + 1 and a record with SVLEN never reads
    /// END (BaseVCF4.pm:302-308); a tandem repeat whose SVLEN is negative ends before
    /// its start and fails the line as Ensembl VEP skips it; and a span from SVLEN
    /// saturates at the ceiling, so a POS just below it with SVLEN `i64::MIN` ends
    /// there whatever the class. A breakpoint's deleted END names its mate with CHR2
    /// (Parser/VCF.pm:493): a run of digits past the ceiling is a mate whose text keeps
    /// the digits (Ensembl VEP 116.2 writes `N[21:18446744073709551615[`) and whose
    /// position saturates at the ceiling; `-5` or `abc` names none and the record's
    /// string is `N/N[21:-5[`.
    #[test]
    fn coordinates_are_read_up_to_the_ceiling() {
        let past = [MAX_COORD + 1, i64::MAX as u64, u64::MAX];
        for pos in past.iter().map(|n| n.to_string()).chain(["-0".to_string()]) {
            let err = parse_vcf_line(&format!("21\t{pos}\t.\tA\tG\t.\t.\t."), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("Invalid POS"), "{pos}: {err}");
        }
        for end in past
            .iter()
            .map(|n| n.to_string())
            .chain(["-5".to_string(), "abc".to_string()])
        {
            let err = parse_vcf_line(
                &format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END={end}"),
                true,
            )
            .unwrap_err()
            .to_string();
            assert!(err.contains("Invalid END"), "{end}: {err}");
            let breakpoint = parse_vcf_line(
                &format!("21\t100\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END={end}"),
                true,
            )
            .unwrap();
            let digits = end.bytes().all(|b| b.is_ascii_digit());
            assert_eq!(breakpoint.len(), if digits { 2 } else { 1 }, "{end}");
            assert_eq!(
                (breakpoint[0].start, breakpoint[0].end),
                (101, 101),
                "{end}"
            );
            assert_eq!(
                breakpoint[0].record_allele_string(),
                format!("N/N[21:{end}["),
                "{end}"
            );
            if digits {
                assert_eq!(breakpoint[0].display_allele(), "N.", "{end}");
                assert_eq!(breakpoint[1].display_allele(), format!("N[21:{end}["));
                assert_eq!(breakpoint[1].mate_pos, Some(MAX_COORD), "{end}");
            } else {
                assert_eq!(
                    breakpoint[0].display_allele(),
                    "chromosome_breakpoint",
                    "{end}"
                );
            }
            let with_svlen = parse_vcf_line(
                &format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=100;END={end}"),
                true,
            )
            .unwrap();
            assert_eq!(
                (with_svlen[0].start, with_svlen[0].end),
                (101, 200),
                "{end}"
            );
        }
        let err = parse_vcf_line(
            "21\t100\t.\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;SVLEN=-300;RUS=AT;RUC=1",
            true,
        )
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("start > end+1 : (START=101, END=-200)"),
            "{err}"
        );
        let kept = parse_vcf_line(
            &format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END={MAX_COORD}"),
            true,
        )
        .unwrap();
        assert_eq!(kept[0].end, MAX_COORD);
        let svlen = format!("SVLEN={}", i64::MIN);
        let below = MAX_COORD - 1;
        let widest = parse_vcf_line(
            &format!("21\t{below}\t.\tN\t<DUP>\t.\t.\tSVTYPE=DUP;{svlen}"),
            true,
        )
        .unwrap();
        assert_eq!(widest[0].sv_end, Some(MAX_COORD));
        let insertion = parse_vcf_line(
            &format!("21\t{below}\t.\tN\t<INS>\t.\t.\tSVTYPE=INS;{svlen}"),
            true,
        )
        .unwrap();
        assert_eq!(
            (insertion[0].start, insertion[0].end),
            (MAX_COORD, MAX_COORD)
        );
        let breakend = parse_vcf_line(
            &format!(
                "21\t{}\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;{svlen};CHR2=21;END2={MAX_COORD}",
                MAX_COORD - 2
            ),
            true,
        )
        .unwrap();
        assert_eq!(breakend.len(), 2);
        assert_eq!(breakend[1].end, MAX_COORD);
        assert_eq!(breakend[1].mate_pos, Some(MAX_COORD));
    }

    /// A bare `<BND>`, with neither CHR2 nor a span, names no mate: `Parser/VCF.pm`
    /// leaves the ALT as written without CHR2 (the `:491-498` block that rewrites it
    /// is skipped) and prepends REF (`:499`), `_parse_breakends` finds no bracket in
    /// `N/<BND>` (`StructuralVariationFeature.pm:1279-1283`), and the row's Allele
    /// falls to `class_SO_term` (`OutputFactory.pm:2029`), `chromosome_breakpoint`,
    /// with or without INFO/MATEID.
    #[test]
    fn bare_symbolic_bnd_prints_the_class_term() {
        for info in [
            "SVTYPE=BND",
            "SVTYPE=BND;MATEID=bnd_x",
            "SVTYPE=BND;END=101",
        ] {
            let line = format!("21\t100\tbare\tN\t<BND>\t.\t.\t{info}");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 1, "{info}");
            assert_eq!(variants[0].mate_id, None, "{info}");
            assert_eq!(
                variants[0].display_allele(),
                "chromosome_breakpoint",
                "{info}"
            );
        }
        let with_mate = parse_vcf_line(
            "21\t100\tpaired\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=500",
            true,
        )
        .unwrap();
        assert_eq!(with_mate[0].display_allele(), "N.");
        assert_eq!(with_mate[1].display_allele(), "N[21:500[");
    }

    /// ensembl-io `BaseVCF4::get_end`: SVLEN first (`start + |SVLEN| - 1`), then END,
    /// then `get_start` for a symbolic record with neither. A symbolic insertion whose
    /// END is its position (1000 Genomes `<INS:MT>`) or whose SVLEN is 0 (`<INS:ME:ALU>`)
    /// is therefore the pair `start = POS + 1, end = POS`, printed as the two flanks;
    /// one with neither field is the single base at POS + 1.
    #[test]
    fn test_parse_vcf_line_symbolic_insertion_end_at_pos_keeps_the_pair() {
        let line = "11\t7309883\t.\tT\t<INS:MT>\t.\t.\tCIEND=-1,0;CIPOS=0,1;END=7309883;SVTYPE=INS;IMPRECISE";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!(v.variant_class, VariantClass::StructuralInsertion);
        assert_eq!(
            (v.start, v.end, v.sv_end),
            (7309884, 7309883, Some(7309883))
        );
        assert_eq!(v.location(), "11:7309883-7309884");

        let line = "17\t78056052\t.\tG\t<INS:ME:ALU>\t.\t.\tMEINFO=AluY,0,0,-;SVLEN=0;SVTYPE=ALU;TSD=CAGGCGTG";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!(v.variant_class, VariantClass::MobileElementInsertion);
        assert_eq!((v.start, v.end), (78056053, 78056052));
        assert_eq!(v.location(), "17:78056052-78056053");

        let line = "21\t25880572\tsynth_ins_sva_0001\tT\t<INS:ME:SVA>\t.\tPASS\tSVTYPE=INS;MEINFO=SVA_E,64,2808,-";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!((v.start, v.end), (25880573, 25880573));
        assert_eq!(v.location(), "21:25880573");

        // A ranged insertion: SVLEN gives end = POS + SVLEN.
        let line = "21\t1000\t.\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=300";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!((v.start, v.end), (1001, 1300));
    }

    #[test]
    fn test_parse_vcf_line_bracket_bnd_takes_dual_emission_path() {
        // Bracket BND notation takes the dual-emission path, not the symbolic one.
        let line = "13\t500\tbnd_1\tA\tA[21:100[\t.\t.\tSVTYPE=BND;MATEID=bnd_2;SVLEN=5000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(paired.alt_allele(), b"A[21:100[");
        assert_eq!(paired.start, 501); // POS+1
                                       // Bracket BND: END falls back to POS+abs(SVLEN)=5500, clamped to max(sv_start)
        assert_eq!(paired.end, 5500);
        assert!(!paired.is_single_breakend);

        let single = &variants[0];
        assert_eq!(single.alt_allele(), b"A.");
        assert!(single.is_single_breakend);
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_with_chr2_end2() {
        // gnomAD <BND> with CHR2/END2: Perl VCF.pm uses END2 as the mate position
        // instead of the SVLEN-derived span end. This converts <BND> to a standard
        // paired BND (e.g., N[21:42458165[).
        let line = "21\t15727643\tBND_21_56725\tN\t<BND>\t.\t.\tSVTYPE=BND;END=15727644;SVLEN=26730521;CHR2=21;END2=42458165";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "symbolic <BND> with CHR2/END2 should emit 2 variants"
        );

        let paired = &variants[1];
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert_eq!(paired.start, 15727644); // POS+1
        assert!(!paired.is_single_breakend);
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        assert_eq!(paired.mate_pos, Some(42458165)); // END2, not end+1
        let alt_str = std::str::from_utf8(paired.alt_allele()).unwrap();
        assert!(
            alt_str.contains("42458165"),
            "paired ALT should use END2 position: {alt_str}"
        );

        let single = &variants[0];
        assert!(single.is_single_breakend);
        assert_eq!(single.alt_allele(), b"N.");
        assert!(single.mate_chr.is_none());
    }

    /// With CHR2 but no END2, the mate position is the raw `INFO/END`.
    ///
    /// Perl's breakpoint branch (VCF.pm:483-498) stashes `INFO/END` as
    /// `$incorrect_end`, deletes it so it cannot be read as a span end, and then
    /// resolves `my $breakend_pos = $info->{END2} || $incorrect_end;`. POS2 is
    /// never consulted, so a record carrying POS2 but not END2 still resolves to
    /// `INFO/END`.
    #[test]
    fn test_parse_vcf_line_symbolic_bnd_pos2_no_end2_uses_info_end() {
        let line = "21\t5000000\tBND_test\tN\t<BND>\t.\t.\tSVTYPE=BND;END=5000001;SVLEN=100000;CHR2=21;POS2=5100001";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(
            paired.mate_pos,
            Some(5000001),
            "Perl resolves END2 || INFO/END and never reads POS2 on this branch"
        );
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_inter_chromosomal_chr2() {
        let line =
            "21\t1000000\tBND_inter\tN\t<BND>\t.\t.\tSVTYPE=BND;END=1000001;SVLEN=500;CHR2=chrX;END2=50000000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        // The mate is carried as written; `pipeline::name_source_chromosomes`
        // resolves it against the cache.
        assert_eq!(paired.mate_chr.as_deref(), Some("chrX"));
        assert_eq!(paired.mate_pos, Some(50000000));
        // The allele string spells the mate chromosome as the input did.
        assert_eq!(paired.alt_allele(), b"N[chrX:50000000[");
        assert_eq!(paired.allele_string, "N/N[chrX:50000000[");
    }

    /// CHR2 with neither END2 nor POS2 resolves to the raw `INFO/END`.
    ///
    /// CHR2 and END present, END2 absent: the mate position is the raw INFO/END
    /// (5000001 for the record below), not POS + |SVLEN| (5010001).
    #[test]
    fn test_parse_vcf_line_symbolic_bnd_chr2_without_end2_uses_info_end() {
        let line =
            "21\t5000000\tBND_noend2\tN\t<BND>\t.\t.\tSVTYPE=BND;END=5000001;SVLEN=10000;CHR2=21";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(
            paired.mate_pos,
            Some(5000001),
            "Perl uses the raw INFO/END as the mate position, not POS + |SVLEN|"
        );
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        let alt_str = std::str::from_utf8(paired.alt_allele()).unwrap();
        assert!(
            alt_str.contains("21:5000001"),
            "paired ALT must reference CHR2:END: {alt_str}"
        );
    }

    /// Guard: END2 still wins over `INFO/END` when both are present.
    ///
    /// `$info->{END2} || $incorrect_end` is an ordered fallback, not a swap: a
    /// record that carries END2 must keep using it, so the rule cannot degenerate
    /// into "always use INFO/END".
    #[test]
    fn test_parse_vcf_line_symbolic_bnd_end2_preferred_over_info_end() {
        let line = "21\t5000000\tBND_both\tN\t<BND>\t.\t.\tSVTYPE=BND;END=5000001;SVLEN=500;CHR2=21;END2=5900000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(
            paired.mate_pos,
            Some(5900000),
            "END2 must take precedence over INFO/END"
        );
    }

    /// A `<NON_REF>` record is a reference-confidence block: one structural record
    /// spanning POS+1..END that no consequence path annotates (VEP's own row for it
    /// is `21:101-200 <NON_REF> - - - intergenic_variant`), and whose VCF and JSON
    /// output VEP leaves without consequences (`vep_skip`).
    #[test]
    fn test_parse_vcf_line_non_ref_is_a_reference_block() {
        let line = "21\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::ReferenceBlock);
        assert!(v.is_structural);
        assert!(v.vep_skip);
        assert_eq!((v.start, v.end), (101, 200));
        assert_eq!(v.sv_end, Some(200));
        assert_eq!(v.location(), "21:101-200");
        assert_eq!(v.display_allele(), "<NON_REF>");
    }

    /// Only `<NON_REF>` (in any case) is the reference block; the other symbolic
    /// alleles keep their classes.
    #[test]
    fn classify_symbolic_alt_reference_block() {
        assert_eq!(
            classify_symbolic_alt("<NON_REF>", None),
            VariantClass::ReferenceBlock
        );
        assert_eq!(
            classify_symbolic_alt("<non_ref>", None),
            VariantClass::ReferenceBlock
        );
        assert_eq!(
            classify_symbolic_alt("<NON_REF>", Some("NON_REF")),
            VariantClass::ReferenceBlock
        );
        // A supported SVTYPE names the record's type (Parser/VCF.pm replaces an ALT
        // outside the VCF 4.4 forms with INFO/SVTYPE before the SO lookup).
        assert_eq!(
            classify_symbolic_alt("<NON_REF>", Some("CNV")),
            VariantClass::CopyNumberVariation
        );
        assert_eq!(
            classify_symbolic_alt("<CNV>", None),
            VariantClass::CopyNumberVariation
        );
        assert_eq!(
            classify_symbolic_alt("<CPX>", Some("CPX")),
            VariantClass::ComplexStructural
        );
    }

    #[test]
    fn vep_skip_marks_unsupported_types_only() {
        for (alts, svtype, supported) in [
            ("<DEL>", Some("DEL"), true),
            ("<INS:ME:ALU>", Some("INS"), true),
            ("<DEL:ME:L1>", None, true),
            ("<DUP:TANDEM>", None, true),
            ("<CNV:TR>", None, true),
            ("<CN0>", None, true),
            ("<CN0>/<CN2>", None, true),
            ("<CNV>", None, true),
            ("<INV>", None, true),
            ("<BND>", Some("BND"), true),
            ("N[21:100[", Some("BND"), true),
            ("<DUP:INT>", None, true),
            ("<CPX>", Some("CPX"), false),
            ("<NON_REF>", None, false),
            ("<TRA>", Some("TRA"), false),
        ] {
            assert_eq!(
                sv_so_term_and_class(&sv_type_abbreviation(alts, svtype)).is_some(),
                supported,
                "{alts} {svtype:?}"
            );
        }
        let cpx = parse_vcf_line("21\t100\t.\tA\t<CPX>\t.\t.\tSVTYPE=CPX;END=500", true).unwrap();
        assert!(cpx[0].vep_skip);
        let del = parse_vcf_line("21\t100\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=500", true).unwrap();
        assert!(!del[0].vep_skip);
    }

    #[test]
    fn test_parse_vcf_line_cpx_is_annotated() {
        // <CPX> variants are annotated by Ensembl VEP (despite vep_skip flag) with
        // overlap-only consequences. Allele displayed as "CPX".
        let line = "21\t100\t.\tA\t<CPX>\t.\t.\tSVTYPE=CPX;END=500";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::ComplexStructural);
        assert!(v.is_structural);
        assert_eq!(v.sv_end, Some(500));
        assert_eq!(v.display_allele(), "CPX");
    }

    #[test]
    fn test_parse_vcf_line_copy_number_classification() {
        // CN0 → DEL, CN2 → DUP, CN3 → generic CNV (matching Perl single-allele behavior)
        let mut cn0_variants =
            parse_vcf_line("21\t300\t.\tA\t<CN0>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        let cn0 = cn0_variants.remove(0);
        assert_eq!(cn0.variant_class, VariantClass::StructuralDeletion);

        let mut cn2_variants =
            parse_vcf_line("21\t300\t.\tA\t<CN2>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        let cn2 = cn2_variants.remove(0);
        assert_eq!(cn2.variant_class, VariantClass::Duplication);

        let mut cn3_variants =
            parse_vcf_line("21\t300\t.\tA\t<CN3>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        let cn3 = cn3_variants.remove(0);
        assert_eq!(cn3.variant_class, VariantClass::CopyNumberVariation);
    }

    /// `<CN0>,<CN2>` is one `copy_number_variation` record (Parser.pm `get_SO_term`:
    /// the joined `<CN0>/<CN2>` is neither `<CN0>` nor `<CN2>` alone) written as one
    /// row per ALT, both labelled `copy_number_variation` (record `multi_cn0_cn2` of the
    /// GRCh38-release116 golden corpus: two rows identical except ALLELE_NUM). Each
    /// variant keeps its own ALT bytes and carries the record's term, which is what
    /// gives `<CN0>` the record's copy-number reading rather than a deletion's.
    #[test]
    fn test_parse_vcf_line_multi_allelic_cn_one_variant_per_alt() {
        let line = "21\t14504804\tDUP_gs_CNV\tC\t<CN0>,<CN2>\t.\t.\tSVTYPE=CNV;END=14530597";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2, "one variant per ALT");
        for (idx, v) in variants.iter().enumerate() {
            assert_eq!(v.allele_index, idx);
            assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
            assert!(v.is_structural);
            assert!(!v.vep_skip);
            assert_eq!(v.display_allele(), "copy_number_variation");
            assert_eq!(v.record_so_term.as_deref(), Some("copy_number_variation"));
            assert_eq!(v.uploaded_allele_string.as_deref(), Some("<CN0>/<CN2>"));
            assert_eq!(v.start, 14504805);
            assert_eq!(v.sv_end, Some(14530597));
        }
        assert_eq!(variants[0].alt_allele(), b"<CN0>");
        assert_eq!(variants[1].alt_allele(), b"<CN2>");
    }

    #[test]
    fn test_parse_vcf_line_multi_allelic_cn_four_alleles() {
        let line = "21\t15244513\t.\tT\t<CN0>,<CN2>,<CN3>,<CN4>\t.\t.\tSVTYPE=CNV;END=15253042";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 4);
        for (v, alt) in variants.iter().zip(["<CN0>", "<CN2>", "<CN3>", "<CN4>"]) {
            assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
            assert_eq!(v.alt_allele(), alt.as_bytes());
            assert_eq!(v.record_so_term.as_deref(), Some("copy_number_variation"));
        }
    }

    /// A `<DEL>,<DUP>` record is one `copy_number_variation` (Parser.pm `get_SO_term`:
    /// DEL beside DUP) whose two variants keep the ALT bytes the line wrote, `<DEL>`
    /// and `<DUP>`, under the record's term; its JSON `allele_string` is the ALTs
    /// joined by `/`. Ensembl VEP 116.2 writes two `copy_number_variation` rows with
    /// `3_prime_UTR_variant` and `IMPACT=MODIFIER` for `21 30165700 multi_del_dup N
    /// <DEL>,<DUP>` against ENST00000286808 in the GRCh38-release116 golden corpus,
    /// the reading of neither a deletion nor a duplication.
    #[test]
    fn del_beside_dup_keeps_each_alt_under_the_record_term() {
        for (alts, expected_alts) in [
            ("<DEL>,<DUP>", ["<DEL>", "<DUP>"]),
            ("<DUP>,<DEL>", ["<DUP>", "<DEL>"]),
            ("<DEL>,<DUP:TANDEM>", ["<DEL>", "<DUP:TANDEM>"]),
        ] {
            let line = format!("21\t30165700\tmulti_del_dup\tN\t{alts}\t.\t.\tEND=30165900");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 2, "{alts}");
            for (v, alt) in variants.iter().zip(expected_alts) {
                assert_eq!(v.alt_allele(), alt.as_bytes(), "{alts}");
                assert_eq!(v.variant_class, VariantClass::CopyNumberVariation, "{alts}");
                assert_eq!(
                    v.record_so_term.as_deref(),
                    Some("copy_number_variation"),
                    "{alts}"
                );
                assert_eq!(v.display_allele(), "copy_number_variation", "{alts}");
                assert_eq!(
                    v.record_allele_string(),
                    alts.replace(',', "/"),
                    "{alts}: JSON allele_string"
                );
            }
        }
    }

    /// Every multi-allelic symbolic shape of ensembl-vep PR 1998 with its type and
    /// row count under 115.2 (the joined string looked up as written: one allele,
    /// and unsupported unless it was a copy-number list) and under 116.2 (identical
    /// parts collapse, DEL beside DUP is a copy-number variation, one allele per
    /// ALT, every row's Allele the record's `class_SO_term`). The 116.2 column is
    /// the GRCh38-release116 golden corpus's `multi_*` records; `<INS>,<DEL>` and
    /// `<INV>,<DUP>` are skipped by 116.2 and stay annotated here with the
    /// `vep_skip` mark, as 115.2 left them. The mobile-element shapes follow
    /// Parser.pm `get_SO_term` (the `(INS|DEL):ME` match runs once over the joined
    /// string and before the DEL-and-DUP test) and `%SO_TERMS` (`Utils/Config.pm`):
    /// `<INS:ME:ALU>,<INS:ME:LINE1>` is `INS_ALU`, `Alu_insertion` on both rows
    /// (OutputFactory.pm `BaseStructuralVariationOverlapAllele_to_output_hash`).
    #[test]
    fn multi_allelic_symbolic_shapes_take_the_record_class() {
        use VariantClass::*;
        // (ALTs, 115.2 type and alleles, 116.2 classes per variant, supported,
        // 116.2 Allele column per variant)
        let shapes = [
            (
                "<INS>,<INS>",
                "unsupported INS/INS, 1 allele",
                vec![StructuralInsertion; 2],
                true,
                vec!["insertion"; 2],
            ),
            (
                "<DEL>,<DEL>",
                "unsupported DEL/DEL, 1 allele",
                vec![StructuralDeletion; 2],
                true,
                vec!["deletion"; 2],
            ),
            (
                "<DEL>,<DUP>",
                "unsupported DEL/DUP, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<DUP>,<DEL>",
                "unsupported DUP/DEL, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<DEL>,<DUP:TANDEM>",
                "tandem_duplication, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<CN0>,<CN2>",
                "copy_number_variation, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<INS>,<DEL>",
                "unsupported INS/DEL, 1 allele",
                vec![StructuralInsertion, StructuralDeletion],
                false,
                vec!["insertion", "deletion"],
            ),
            (
                "<INV>,<DUP>",
                "unsupported INV/DUP, 1 allele",
                vec![Inversion, Duplication],
                false,
                vec!["inversion", "duplication"],
            ),
            (
                "G,<INV>",
                "deletion, 1 allele",
                vec![StructuralDeletion; 2],
                true,
                vec!["deletion"; 2],
            ),
            (
                "<INS:ME:ALU>,<INS:ME:ALU>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<INS:ME:ALU>,<INS:ME:LINE1>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<INS:ME:ALU>,<INS>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<INS>,<INS:ME:ALU>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<DEL:ME:ALU>,<DUP>",
                "Alu_deletion, 1 allele",
                vec![MobileElementDeletion; 2],
                true,
                vec!["Alu_deletion"; 2],
            ),
        ];
        for (alts, release_115, classes, supported, labels) in shapes {
            let line = format!("21\t30165700\t.\tN\t{alts}\t.\t.\tSVTYPE=DEL;END=30165900");
            let variants = parse_vcf_line(&line, true).unwrap();
            let got: Vec<VariantClass> = variants.iter().map(|v| v.variant_class).collect();
            assert_eq!(got, classes, "{alts} (115.2: {release_115})");
            assert!(
                variants.iter().all(|v| v.vep_skip != supported),
                "{alts}: vep_skip must be {}",
                !supported
            );
            let alleles: Vec<usize> = variants.iter().map(|v| v.allele_index).collect();
            assert_eq!(alleles, (0..variants.len()).collect::<Vec<_>>(), "{alts}");
            let shown: Vec<String> = variants.iter().map(|v| v.display_allele()).collect();
            assert_eq!(shown, labels, "{alts}");
            let bytes: Vec<&[u8]> = variants.iter().map(|v| v.alt_allele()).collect();
            let written: Vec<&[u8]> = alts.split(',').map(str::as_bytes).collect();
            assert_eq!(bytes, written, "{alts}: each variant keeps its own ALT");
        }
    }

    /// A record carrying both a mobile-element deletion and a mobile-element
    /// insertion takes the leftmost one: Parser.pm `get_SO_term` matches
    /// `/(INS|DEL):(ME):?([A-Z0-9]+)?/i` once over the joined string (`:685`), so
    /// `<DEL:ME:ALU>,<INS:ME:LINE1>` is `DEL_ALU`, `Alu_deletion` on both rows, and
    /// the reversed order `INS_LINE1`, `LINE1_insertion` on both.
    #[test]
    fn mixed_mobile_element_record_takes_the_leftmost_alt() {
        use VariantClass::*;
        for (alts, abbrev, class, label) in [
            (
                "<DEL:ME:ALU>,<INS:ME:LINE1>",
                "DEL_ALU",
                MobileElementDeletion,
                "Alu_deletion",
            ),
            (
                "<INS:ME:LINE1>,<DEL:ME:ALU>",
                "INS_LINE1",
                MobileElementInsertion,
                "LINE1_insertion",
            ),
        ] {
            assert_eq!(
                sv_type_abbreviation(&alts.replace(',', "/"), Some("DEL")),
                abbrev,
                "{alts}"
            );
            let line = format!("21\t30165700\tme_mix\tN\t{alts}\t.\t.\tEND=30165900;SVTYPE=DEL");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 2, "{alts}");
            for v in &variants {
                assert_eq!(v.variant_class, class, "{alts}");
                assert_eq!(v.record_so_term.as_deref(), Some(label), "{alts}");
                assert_eq!(v.display_allele(), label, "{alts}");
            }
        }
    }

    /// A single mobile-element ALT keeps its subtype label from the same table:
    /// `<INS:ME:ALU>` is `Alu_insertion`, `<INS:ME>` `mobile_element_insertion`,
    /// `<INS:ME:L1>` `LINE1_insertion` (Parser.pm `get_SO_term` maps `L1` to
    /// `LINE1`), `<DEL:ME:SVA>` `SVA_deletion`, and an element outside Parser.pm's
    /// list (`<INS:ME:LINE>`) the bare `mobile_element_insertion`. Ensembl VEP
    /// 116.2 prints `Alu_insertion` for `21 10524521 gnomAD-SV_v3_INS_chr21_1d284ae3
    /// N <INS:ME:ALU>` and `mobile_element_insertion` for `21 31659650
    /// synth_ins_me_0117 A <INS:ME>` in the GRCh38 golden corpus, and
    /// `mobile_element_insertion` for `21 27284085 synth_ins_line_0011 G
    /// <INS:ME:LINE>` and `mobile_element_deletion` for `21 33041214
    /// synth_del_me_0080 G <DEL:ME>` in the GRCh37 golden corpus.
    #[test]
    fn single_mobile_element_alts_keep_their_subtype_labels() {
        for (alt, class, label) in [
            (
                "<INS:ME:ALU>",
                VariantClass::MobileElementInsertion,
                "Alu_insertion",
            ),
            (
                "<INS:ME>",
                VariantClass::MobileElementInsertion,
                "mobile_element_insertion",
            ),
            (
                "<INS:ME:L1>",
                VariantClass::MobileElementInsertion,
                "LINE1_insertion",
            ),
            (
                "<INS:ME:LINE>",
                VariantClass::MobileElementInsertion,
                "mobile_element_insertion",
            ),
            (
                "<DEL:ME:SVA>",
                VariantClass::MobileElementDeletion,
                "SVA_deletion",
            ),
            (
                "<DEL:ME>",
                VariantClass::MobileElementDeletion,
                "mobile_element_deletion",
            ),
        ] {
            let line = format!("21\t10524521\t.\tN\t{alt}\t.\t.\tSVTYPE=INS;END=10524572");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 1, "{alt}");
            assert_eq!(variants[0].variant_class, class, "{alt}");
            assert_eq!(variants[0].display_allele(), label, "{alt}");
            assert_eq!(variants[0].record_so_term.as_deref(), Some(label), "{alt}");
        }
    }

    /// `<CN0>,<CN0>` collapses to `<CN0>`, a deletion, under 116.2's normalisation
    /// (Parser.pm `get_SO_term`: identical parts reduce to one).
    #[test]
    fn repeated_cn0_collapses_to_a_deletion() {
        let variants =
            parse_vcf_line("21\t100\t.\tA\t<CN0>,<CN0>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        assert_eq!(variants.len(), 2);
        for v in &variants {
            assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
            assert_eq!(v.alt_allele(), b"<CN0>");
            assert!(!v.vep_skip);
        }
    }

    #[test]
    fn test_parse_vcf_line_single_cn0_still_deletion() {
        let line = "21\t300\t.\tA\t<CN0>\t.\t.\tSVTYPE=CNV;END=800";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_class, VariantClass::StructuralDeletion);
    }

    #[test]
    fn test_parse_vcf_line_single_cn2_still_duplication() {
        let line = "21\t300\t.\tA\t<CN2>\t.\t.\tSVTYPE=CNV;END=800";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_class, VariantClass::Duplication);
    }

    #[test]
    fn test_parse_vcf_line_too_few_fields() {
        let line = "21\t25000100\t.\tA";
        assert!(parse_vcf_line(line, true).is_err());
    }

    #[test]
    fn test_parse_vcf_line_tandem_repeat_uses_svlen_for_end() {
        // VCF.pm:502-512: $end = $start + max(@svlen) - 1 when $so_term =~ /tandem/
        let line = "21\t27253937\tsynth_cnv_tr_0001\tC\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;END=27256025;SVLEN=307";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::TandemRepeat);
        // sv_start = POS+1 = 27253938, sv_end = sv_start + SVLEN - 1 = 27254244
        assert_eq!(v.start, 27253938);
        assert_eq!(v.sv_end, Some(27254244));
        // Must not use the END field value
        assert_ne!(v.sv_end, Some(27256025));
    }

    #[test]
    fn test_parse_vcf_line_cnv_svlen_over_end() {
        // SVLEN=307 and END=27256025: SVLEN wins (BaseVCF4.pm get_end);
        // end = POS + abs(SVLEN) = 27253937 + 307 = 27254244
        let line = "21\t27253937\t.\tC\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=27256025;SVLEN=307";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
        assert_eq!(v.sv_end, Some(27254244));
    }

    /// With the reference sequence to hand a `<CNV:TR>` record is the sequence record
    /// Ensembl VEP 116.2 builds over its run (the run `TCATCATCA` at
    /// 21:13620284-13620292 of GRCh38): `RUC=4` and `RB=12` are
    /// `TCATCATCA/TCATCATCATCA`, minimised to the insertion `-/TCA` after the run
    /// (`21:19-20`); `RUC=2` the deletion `TCA/-` of its last unit (`21:17-19`);
    /// `RUC=3` an allele equal to the reference, a feature with
    /// no allele; `RUC=4.9` four units (`x` truncates); `RUS=TCA;RUC=5` without RN one
    /// sequence of five; `RN=2;RUS=TCA,TC;RUC=3,1` the two sequences joined, `-/TC`;
    /// `RUS=.` the base `N` four times, `TCATCATCA/NNNN` left unminimised; a lower-case
    /// `RUS=tca` never trims against the upper-case run and is upper-cased after; two
    /// ALTs two alleles over the run, left unminimised, with the chopped string as the
    /// record's; `END` without `SVLEN` the same run; and `UPLOADED_ALLELE` prints the
    /// expansion as built.
    #[test]
    fn cnv_tr_expanded_with_fasta() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let record = |alt: &str, info: &str| {
            parse_with_fasta(
                &format!("21\t10\tx\tG\t{alt}\t.\t.\t{info}"),
                Some(&fasta),
                Some(DEFAULT_MAX_SV_SIZE),
            )
        };
        let one = |info: &str| shapes(&record("<CNV:TR>", info));
        let ins = |alt: &str| {
            vec![(
                "21:19-20".to_string(),
                format!("-/{alt}"),
                VariantClass::Insertion,
                false,
            )]
        };
        assert_eq!(one("SVLEN=9;RN=1;RUS=TCA;RUC=4"), ins("TCA"));
        assert_eq!(one("SVLEN=9;RUS=TCA;RB=12"), ins("TCA"));
        assert_eq!(one("SVLEN=9;RN=1;RUS=TCA;RUC=4.9"), ins("TCA"));
        assert_eq!(one("END=19;RN=1;RUS=TCA;RUC=4"), ins("TCA"));
        assert_eq!(one("SVLEN=9;RUS=TCA;RUC=5"), ins("TCATCA"));
        assert_eq!(one("SVLEN=9;RN=2;RUS=TCA,TC;RUC=3,1"), ins("TC"));
        assert_eq!(
            one("SVLEN=9;RN=1;RUS=TCA;RUC=2"),
            [(
                "21:17-19".to_string(),
                "TCA/-".to_string(),
                VariantClass::Deletion,
                false
            )]
        );
        let same = record("<CNV:TR>", "SVLEN=9;RN=1;RUS=TCA;RUC=3");
        assert!(same[0].no_allele);
        assert_eq!(same[0].allele_string, "TCATCATCA/TCATCATCA");
        assert_eq!(
            one("SVLEN=9;RN=1;RUS=.;RUC=4"),
            [(
                "21:11-19".to_string(),
                "TCATCATCA/NNNN".to_string(),
                VariantClass::Indel,
                false
            )]
        );
        let lower = record("<CNV:TR>", "SVLEN=9;RN=1;RUS=tca;RUC=4");
        assert_eq!(lower[0].allele_string, "TCATCATCA/TCATCATCATCA");
        assert_eq!(lower[0].location(), "21:11-19");
        assert_eq!(
            lower[0].uploaded_allele_string.as_deref(),
            Some("TCATCATCA/tcatcatcatca")
        );
        let gain = record("<CNV:TR>", "SVLEN=9;RN=1;RUS=TCA;RUC=4");
        assert_eq!(
            gain[0].uploaded_allele_string.as_deref(),
            Some("TCATCATCA/TCATCATCATCA")
        );
        assert_eq!(gain[0].input_line_name.as_deref(), Some("21_10_G/<CNV:TR>"));
        let multi = record("<CNV:TR>,<CNV:TR>", "SVLEN=9,9;RN=1,1;RUS=TCA,TCA;RUC=4,2");
        assert_eq!(
            shapes(&multi),
            [
                (
                    "21:11-19".to_string(),
                    "TCATCATCA/TCATCATCATCA".to_string(),
                    VariantClass::Indel,
                    false
                ),
                (
                    "21:11-19".to_string(),
                    "TCATCATCA/TCATCA".to_string(),
                    VariantClass::Indel,
                    false
                ),
            ]
        );
        assert_eq!(multi[0].class_so_term(), "indel");
        assert_eq!(
            multi[1].record_allele_string_multi.as_deref(),
            Some("TCATCATCA/TCATCATCATCA/TCATCA")
        );
        assert_eq!(multi[1].allele_index, 1);
    }

    /// Two SVLEN values give the run the largest (`uniq`, `max`, Parser/VCF.pm:505-511):
    /// `SVLEN=9,6` with `RUC=4,1` is the run `TCATCATCA`, so the second allele `TCA` is
    /// the loss of two units; the reader's warning names the values as written.
    #[test]
    fn cnv_tr_svlen_max_of_uniq() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let line = "21\t10\tx\tG\t<CNV:TR>,<CNV:TR>\t.\t.\tSVLEN=9,6;RN=1,1;RUS=TCA,TCA;RUC=4,1";
        let variants = parse_with_fasta(line, Some(&fasta), Some(DEFAULT_MAX_SV_SIZE));
        assert_eq!(
            shapes(&variants),
            [
                (
                    "21:11-19".to_string(),
                    "TCATCATCA/TCATCATCATCA".to_string(),
                    VariantClass::Indel,
                    false
                ),
                (
                    "21:11-19".to_string(),
                    "TCATCATCA/TCA".to_string(),
                    VariantClass::Indel,
                    false
                ),
            ]
        );
        assert_eq!(
            line_warnings(line, true),
            ["found tandem repeats with different references per alternative allele: SVLEN=9,6; only using reference with largest size"]
        );
        let ciruc = "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=4;CIRUC=3,5";
        assert_eq!(
            line_warnings(ciruc, true),
            ["CIRUC and CIRB INFO fields are ignored when calculating alternative alleles in tandem repeats"]
        );
        assert!(line_warnings(ciruc, false)
            .iter()
            .any(|w| w.starts_with("could not fetch sequence for tandem repeats;")));
        assert!(line_warnings("21\t10\tx\tG\t<DEL>\t.\t.\tSVLEN=-9", false).is_empty());
    }

    /// The warnings follow Ensembl VEP's routing (Parser/VCF.pm:237-245, :473-474):
    /// `N N` under `SVTYPE=CNV:TR` is a tandem repeat (the alleles are not all
    /// `ACGT` and the type comes from SVTYPE), as is `<CNV:TR>` under `SVTYPE=DEL`
    /// (the ALT names a VCF 4.4 type, so SVTYPE is not read); `G A` under
    /// `SVTYPE=CNV:TR` is a sequence record and `<DEL>` under `SVTYPE=CNV:TR` a
    /// deletion, and neither warns. Ensembl VEP 116.2 writes the no-FASTA notice
    /// for the first two records and three `tandem_repeat` rows for each.
    #[test]
    fn tandem_warnings_follow_the_structural_routing() {
        let warns = |line: &str| !line_warnings(line, false).is_empty();
        assert!(warns(
            "21\t13620283\tx\tN\tN\t.\t.\tSVTYPE=CNV:TR;SVLEN=9;RN=1;RUS=TCA;RUC=4"
        ));
        assert!(warns(
            "21\t13620283\tx\tG\t<CNV:TR>\t.\t.\tSVTYPE=DEL;SVLEN=9;RN=1;RUS=TCA;RUC=4"
        ));
        assert!(!warns(
            "21\t13620283\tx\tG\tA\t.\t.\tSVTYPE=CNV:TR;SVLEN=9;RN=1;RUS=TCA;RUC=4"
        ));
        assert!(!warns(
            "21\t13620283\tx\tN\t<DEL>\t.\t.\tSVTYPE=CNV:TR;SVLEN=9"
        ));
    }

    /// The SVLEN values are read as Perl's `split` reads them (`:505`): `SVLEN=9,`
    /// is the one value `9`, so it warns of no second reference, where `SVLEN=,9`
    /// and `SVLEN=9,,9` each hold an empty value beside the `9` and warn.
    #[test]
    fn a_trailing_empty_svlen_value_is_one_reference() {
        let svlen_warning = |svlen: &str| {
            let line =
                format!("21\t13620283\tx\tG\t<CNV:TR>\t.\t.\tSVLEN={svlen};RN=1;RUS=TCA;RUC=4");
            line_warnings(&line, true)
                .iter()
                .any(|w| w.starts_with("found tandem repeats with different references"))
        };
        assert!(!svlen_warning("9,"));
        assert!(!svlen_warning("9"));
        assert!(svlen_warning(",9"));
        assert!(svlen_warning("9,,9"));
        assert!(svlen_warning("9,6"));
    }

    /// Without RN the number of repeat sequences is the number of `CNV:TR` in Perl's
    /// type string (`() = $type =~ /CNV:TR/g`, Parser/VCF.pm:397), INFO/SVTYPE when the
    /// ALT opens outside the VCF 4.4 forms (`:473-474`), matched with case where
    /// `get_SO_term` matches without it: `N N` under `SVTYPE=CNV:TR` builds one
    /// sequence, the insertion `-/TCA` Ensembl VEP 116.2 writes for it, and `<cnv:tr>`
    /// is a tandem repeat whose type string holds no `CNV:TR`, so its expansion `REF/`
    /// names no allele.
    #[test]
    fn repeat_sequences_without_rn_are_counted_in_perls_type_string() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let record = |alt: &str, info: &str| {
            parse_with_fasta(
                &format!("21\t10\tx\tN\t{alt}\t.\t.\t{info}"),
                Some(&fasta),
                Some(DEFAULT_MAX_SV_SIZE),
            )
        };
        assert_eq!(
            shapes(&record("N", "SVTYPE=CNV:TR;SVLEN=9;RUS=TCA;RUC=4")),
            [(
                "21:19-20".to_string(),
                "-/TCA".to_string(),
                VariantClass::Insertion,
                false
            )]
        );
        let lower = record("<cnv:tr>", "SVLEN=9;RUS=TCA;RUC=4");
        assert!(lower[0].no_allele, "{lower:?}");
        assert_eq!(lower[0].allele_string, "TCATCATCA/");
        assert_eq!(
            shapes(&record("<CNV:TR>", "SVTYPE=cnv:tr;SVLEN=9;RUS=TCA;RUC=4")),
            [(
                "21:19-20".to_string(),
                "-/TCA".to_string(),
                VariantClass::Insertion,
                false
            )]
        );
    }

    /// Without a FASTA, under `--max_sv_size -1`, under a limit a repeat sequence
    /// exceeds (`RUC=4`, twelve bases, at 10; `RUC=2` expands) and without RUS or a
    /// count list, a `<CNV:TR>` record keeps Ensembl VEP's symbolic reading over the run
    /// SVLEN spans (Parser/VCF.pm:385, :403, :417-418, :428, :521-527); one oversized
    /// allele keeps every allele of its record symbolic; a cap of `0` is Perl's 5000.
    #[test]
    fn cnv_tr_stays_sv_without_fasta() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let gain = "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=4";
        let symbolic = vec![(
            "21:11-19".to_string(),
            "G/<CNV:TR>".to_string(),
            VariantClass::TandemRepeat,
            true,
        )];
        assert_eq!(
            shapes(&parse_with_fasta(gain, None, Some(DEFAULT_MAX_SV_SIZE))),
            symbolic
        );
        assert_eq!(
            shapes(&parse_with_fasta(gain, Some(&fasta), None)),
            symbolic
        );
        assert_eq!(
            shapes(&parse_with_fasta(gain, Some(&fasta), Some(10))),
            symbolic
        );
        assert_eq!(
            shapes(&parse_with_fasta(
                "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=2",
                Some(&fasta),
                Some(10)
            )),
            [(
                "21:17-19".to_string(),
                "TCA/-".to_string(),
                VariantClass::Deletion,
                false
            )]
        );
        for info in [
            "SVLEN=9",
            "SVLEN=9;RUS=TCA",
            "SVLEN=9;RUC=4",
            "SVLEN=9;RUS=TCA;RUC=0",
        ] {
            let line = format!("21\t10\tx\tG\t<CNV:TR>\t.\t.\t{info}");
            assert_eq!(
                shapes(&parse_with_fasta(
                    &line,
                    Some(&fasta),
                    Some(DEFAULT_MAX_SV_SIZE)
                )),
                symbolic,
                "{info}"
            );
        }
        let multi = "21\t10\tx\tG\t<CNV:TR>,<CNV:TR>\t.\t.\tSVLEN=9;RN=1,1;RUS=TCA,TCA;RUC=4,2";
        let both = parse_with_fasta(multi, Some(&fasta), Some(10));
        assert_eq!(both.len(), 2);
        assert!(both.iter().all(|v| v.is_structural));
        let five_thousand = "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RUS=AT;RUC=2500";
        assert!(!parse_with_fasta(five_thousand, Some(&fasta), Some(0))[0].is_structural);
        let over = "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RUS=AT;RUC=2501";
        assert!(parse_with_fasta(over, Some(&fasta), Some(0))[0].is_structural);
    }

    /// `_expand_tandem_repeat_allele_string` builds one allele per RN entry (`:408`), so
    /// a `<CNV:TR>` ALT past a short RN list has no allele: `RN=1;RUS=AT,CAGA;RUC=10,4000000`
    /// on two ALTs is the one allele `AT x 10`, whose 16,000,000-base sibling is never
    /// built and so never oversizes the record; `RN=1,1` builds both, and two alleles
    /// that expand to one sequence are one allele, as a repeated ALT is
    /// (VariationFeatureOverlap.pm:461).
    #[test]
    fn an_alt_past_the_rn_list_builds_no_allele() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let record = |info: &str| {
            let line = format!("21\t10\tx\tG\t<CNV:TR>,<CNV:TR>\t.\t.\tSVLEN=9;{info}");
            parse_with_fasta(&line, Some(&fasta), Some(DEFAULT_MAX_SV_SIZE))
        };
        let one = record("RN=1;RUS=AT,CAGA;RUC=10,4000000");
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].allele_string, "TCATCATCA/ATATATATATATATATATAT");
        assert_eq!(record("RN=1;RUS=AT,CA;RUC=10,10").len(), 1);
        assert_eq!(record("RN=1,1;RUS=AT,CA;RUC=10,10").len(), 2);
        assert_eq!(record("RUS=AT,CA;RUC=10,10").len(), 2);
        assert_eq!(record("RN=1,1;RUS=AT,AT;RUC=10,10").len(), 1);
    }

    /// A span the FASTA lacks leaves the record's REF column as the run's reference, as
    /// `sub_Slice` yields nothing and `$ref` keeps the column (Parser/VCF.pm:516-518):
    /// `G/TCATCATCATCA` over 21:31-39, which `validate_vf` then refuses as a reference
    /// allele whose length is not the span's. A negative SVLEN puts the end before the
    /// start (`:511`), so `SVLEN=-9` at 10 is `G/TCATCATCATCA` over 11-1, which
    /// `validate_vf` refuses as `start > end+1`, as Ensembl VEP 116.2 refuses it with
    /// a FASTA; an end before the chromosome's first base (`SVLEN=-20` at 10) is the
    /// symbolic path's refusal, `start > end+1 : (START=11, END=-10)`.
    #[test]
    fn a_run_past_the_fasta_keeps_the_ref_column() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let past = parse_with_fasta(
            "21\t30\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=4",
            Some(&fasta),
            Some(DEFAULT_MAX_SV_SIZE),
        );
        assert_eq!(
            shapes(&past),
            [(
                "21:31-39".to_string(),
                "G/TCATCATCATCA".to_string(),
                VariantClass::Indel,
                false
            )]
        );
        let negative = parse_with_fasta(
            "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=-9;RN=1;RUS=TCA;RUC=4",
            Some(&fasta),
            Some(DEFAULT_MAX_SV_SIZE),
        );
        assert_eq!(negative.len(), 1);
        assert_eq!((negative[0].start, negative[0].end), (11, 1));
        assert_eq!(negative[0].allele_string, "G/TCATCATCATCA");
        assert!(!negative[0].is_structural);
        let before_the_chromosome = parse_vcf_line_with(
            "21\t10\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=-20;RN=1;RUS=TCA;RUC=4",
            &VcfParse {
                fasta: Some(&fasta),
                ..VcfParse::default_run(true)
            },
        )
        .unwrap_err()
        .to_string();
        assert!(
            before_the_chromosome.contains("start > end+1 : (START=11, END=-10)"),
            "{before_the_chromosome}"
        );
    }

    /// A REF column that shares a base with the expansion is trimmed as any
    /// bi-allelic indel is (`post_process_vfs`, Parser.pm:862-881), the start moving
    /// past the bases trimmed on the left and the end before those trimmed on the
    /// right, the run's own end kept otherwise (Sequence.pm:1004-1031): Ensembl VEP
    /// 116.2 with a FASTA refuses `T <CNV:TR>` with `SVLEN=-9` as `start > end+1`
    /// with the start one past the run's and the end the run's, `A <CNV:TR>` with the
    /// start the run's and the end one before it, and a run past the FASTA from a `T`
    /// as an insertion whose coordinates are not `start = end + 1`; none of the three
    /// is an insertion on a span its alleles measure.
    #[test]
    fn a_ref_column_sharing_a_base_with_the_expansion_keeps_the_runs_end() {
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let span_and_string = |line: &str| {
            let parsed = parse_with_fasta(line, Some(&fasta), Some(DEFAULT_MAX_SV_SIZE));
            assert_eq!(parsed.len(), 1, "{line}");
            (
                parsed[0].start,
                parsed[0].end,
                parsed[0].allele_string.clone(),
            )
        };
        assert_eq!(
            span_and_string("21\t10\tx\tT\t<CNV:TR>\t.\t.\tSVLEN=-9;RN=1;RUS=TCA;RUC=4"),
            (12, 1, "-/CATCATCATCA".to_string())
        );
        assert_eq!(
            span_and_string("21\t10\tx\tA\t<CNV:TR>\t.\t.\tSVLEN=-9;RN=1;RUS=TCA;RUC=4"),
            (11, 0, "-/TCATCATCATC".to_string())
        );
        assert_eq!(
            span_and_string("21\t30\tx\tT\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=4"),
            (32, 39, "-/CATCATCATCA".to_string())
        );
    }

    #[test]
    fn test_parse_vcf_line_tandem_repeat_svlen_only() {
        let line = "21\t27253937\t.\tC\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=500";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::TandemRepeat);
        // sv_start = POS+1 = 27253938, sv_end = 27253938 + 500 - 1 = 27254437
        assert_eq!(v.sv_end, Some(27254437));
    }

    #[test]
    fn test_parse_vcf_line_complex_substitution_suffix_trimmed() {
        // VCF: TAA -> GCTCA (3bp ref, 5bp alt: complex substitution)
        // Ensembl VEP trims trailing A: ref=TA, alt=GCTC, end=27256844
        let line = "21\t27256843\t.\tTAA\tGCTCA\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"TA");
        assert_eq!(v.alt_alleles[0], b"GCTC");
        assert_eq!(v.start, 27256843);
        assert_eq!(v.end, 27256844); // 27256843 + 2 - 1
        assert_eq!(v.allele_string, "TA/GCTC");
        assert!(v.minimised);
    }

    #[test]
    fn test_parse_vcf_line_simple_deletion_no_extra_suffix_trim() {
        // VCF: ACGT -> A (simple deletion: prefix trim removes A, leaving CGT/-)
        // Already minimal after prefix trim, no suffix trimming needed
        let line = "1\t100\t.\tACGT\tA\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"CGT");
        assert_eq!(v.alt_alleles[0], b"-");
        assert_eq!(v.start, 101);
        assert_eq!(v.end, 103);
    }

    #[test]
    fn test_parse_vcf_line_mnp_no_trimming() {
        // Same-length MNP: no trimming at all (ref.len() == alt.len())
        let line = "1\t100\t.\tAC\tGT\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"AC");
        assert_eq!(v.alt_alleles[0], b"GT");
        assert_eq!(v.start, 100);
        assert_eq!(v.end, 101);
        assert!(!v.minimised);
    }

    #[test]
    fn test_parse_vcf_line_insertion_with_shared_suffix() {
        // VCF: AG -> TCCAG (2bp ref, 5bp alt)
        // Multi-allele trim: A vs T no match -> prefix_len=0, ref=AG, alt=TCCAG
        // Secondary minimisation: no prefix match, suffix G matches G -> trim,
        // then suffix A matches A -> trim -> ref=-, alt=TCC
        let line = "1\t100\t.\tAG\tTCCAG\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"-");
        assert_eq!(v.alt_alleles[0], b"TCC");
        assert!(v.minimised);
    }

    /// A bi-allelic record whose REF is a strict prefix of its ALT must be
    /// re-minimised to a pure insertion.
    ///
    /// Perl runs two minimisation stages. `VEP/Parser/VCF.pm:328-335` strips one
    /// shared leading base, then `post_process_vfs` (`VEP/Parser.pm:860-882`) calls
    /// `minimise_alleles` -> `trim_sequences` (`Sequence.pm:1004-1035`), which is
    /// not gated on `--minimal` and trims from the left while both alleles are
    /// non-empty. With `empty_to_dash` set, `T/TT` therefore becomes `-/T`; Perl
    /// does not keep a base on the shorter allele.
    #[test]
    fn test_parse_vcf_line_biallelic_prefix_ref_minimises_to_insertion() {
        // AT -> ATT: stage 1 strips the shared leading A leaving T/TT, then
        // trim_sequences strips the shared leading T, emptying ref.
        let line = "1\t100\t.\tAT\tATT\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(
            v.ref_allele,
            b"-",
            "Perl's trim_sequences empties the ref and empty_to_dash maps it to \
             '-'; got ref={:?} alt={:?}",
            String::from_utf8_lossy(&v.ref_allele),
            String::from_utf8_lossy(&v.alt_alleles[0])
        );
        assert_eq!(v.alt_alleles[0], b"T");
        // VEP insertion convention: end == start - 1.
        assert_eq!(v.start, 102);
        assert_eq!(v.end, 101);
        assert!(v.minimised);
    }

    /// `trim_sequences` compares the bases as written (`eq`, Sequence.pm:1006, :1019)
    /// and `validate_vf` upper-cases the string afterwards (Parser.pm:554), so a
    /// lowercase ALT is minimised only where its letters match the reference's case:
    /// Ensembl VEP 116.2 annotates `21 30396602 x TTGA t` as `TTGA/T` at
    /// 30396602-30396605 (`stop_lost,inframe_deletion`), `TTGA Ttgatc` as `TGA/TGATC`
    /// at 30396603-30396605 (`stop_retained_variant`), and the upper-cased
    /// `TTGA TTGATC` as `-/TC` at 30396606-30396605.
    #[test]
    fn minimisation_compares_the_written_case_as_trim_sequences_does() {
        let lower = parse_vcf_line("21\t30396602\tx\tTTGA\tt\t.\t.\t.", true).unwrap();
        assert_eq!((lower[0].start, lower[0].end), (30396602, 30396605));
        assert_eq!(lower[0].allele_string, "TTGA/T");
        assert!(!lower[0].minimised);
        let mixed = parse_vcf_line("21\t30396602\tx\tTTGA\tTtgatc\t.\t.\t.", true).unwrap();
        assert_eq!((mixed[0].start, mixed[0].end), (30396603, 30396605));
        assert_eq!(mixed[0].allele_string, "TGA/TGATC");
        assert!(!mixed[0].minimised);
        let upper = parse_vcf_line("21\t30396602\tx\tTTGA\tTTGATC\t.\t.\t.", true).unwrap();
        assert_eq!((upper[0].start, upper[0].end), (30396606, 30396605));
        assert_eq!(upper[0].allele_string, "-/TC");
        assert!(upper[0].minimised);
    }

    /// Under `--minimal` every bi-allelic record goes through `minimise_alleles`
    /// (Parser.pm:852, :987-1014): `GAAG/GGGG` at 10542414 is two-side trimmed to
    /// `AA/GG` at 10542415-10542416 (Ensembl VEP 116.2's `21:10542415-10542416 GG`),
    /// `GAA/GGG` and `AAG/GGG` lose one side, `GAAG/GGGg` keeps its right `G` (the
    /// bases compare as written, `GGG` at 15-17), `gaag/gggg` is `AA/GG` and an SNV
    /// keeps its alleles; each is marked minimised, trimmed or not, and
    /// `UPLOADED_ALLELE` keeps the raw string. Without the flag `GAAG/GGGG` is left
    /// as it is.
    #[test]
    fn minimal_two_sided_trims_mnp() {
        let record = |alleles: &str| {
            let (reference, alt) = alleles.split_once('/').unwrap();
            parse_minimal(&format!("21\t10542414\tx\t{reference}\t{alt}\t.\t.\t."))
        };
        let shape = |alleles: &str| {
            let v = record(alleles);
            (v[0].location(), v[0].allele_string.clone(), v[0].minimised)
        };
        assert_eq!(
            shape("GAAG/GGGG"),
            ("21:10542415-10542416".into(), "AA/GG".into(), true)
        );
        assert_eq!(
            shape("GAA/GGG"),
            ("21:10542415-10542416".into(), "AA/GG".into(), true)
        );
        assert_eq!(
            shape("AAG/GGG"),
            ("21:10542414-10542415".into(), "AA/GG".into(), true)
        );
        assert_eq!(
            shape("GAAG/GGGg"),
            ("21:10542415-10542417".into(), "AAG/GGG".into(), true)
        );
        assert_eq!(
            shape("gaag/gggg"),
            ("21:10542415-10542416".into(), "AA/GG".into(), true)
        );
        assert_eq!(shape("A/G"), ("21:10542414".into(), "A/G".into(), true));
        assert_eq!(
            record("gaag/gggg")[0].uploaded_allele(true).as_deref(),
            Some("gaag/gggg")
        );
        let plain = parse_vcf_line("21\t10542414\tx\tGAAG\tGGGG\t.\t.\t.", true).unwrap();
        assert_eq!(plain[0].allele_string, "GAAG/GGGG");
        assert!(!plain[0].minimised);
    }

    /// A bi-allelic record whose one allele equals the reference is minimised under
    /// `--minimal` like any other (Parser.pm:987-1011): `21 10542414 x A A` is the
    /// feature with no allele `-/-` at 10542415-10542414, marked minimised, its class
    /// `sequence_alteration`, as Ensembl VEP 116.2's JSON object carries it; without
    /// the flag it stays `A/A` at 10542414. `A *` trims nothing and is marked too
    /// (`:1011` marks every bi-allelic record).
    #[test]
    fn a_reference_equal_bi_allelic_record_is_minimised_under_the_flag() {
        let line = "21\t10542414\tx\tA\tA\t.\t.\t.";
        let minimal = parse_minimal(line);
        assert_eq!(minimal.len(), 1);
        let record = &minimal[0];
        assert!(record.no_allele);
        assert_eq!(
            (record.start, record.end, record.allele_string.as_str()),
            (10542415, 10542414, "-/-")
        );
        assert!(record.minimised);
        assert_eq!(record.class_so_term(), "sequence_alteration");
        let plain = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            (
                plain[0].start,
                plain[0].end,
                plain[0].allele_string.as_str()
            ),
            (10542414, 10542414, "A/A")
        );
        assert!(!plain[0].minimised);
        let star = parse_minimal("21\t10542414\tx\tA\t*\t.\t.\t.");
        assert!(star[0].no_allele && star[0].minimised);
        assert_eq!(star[0].allele_string, "A/*");
    }

    /// Under `--minimal` a multi-allelic record splits per allele when any allele
    /// trims (InputBuffer.pm:638-718): `G T,GA` at 10542414 keeps `G/T` and trims
    /// `G/GA` to `-/A` at 10542415-10542414, every allele numbered by its written
    /// position, Location printing the record's span (`21:10542414` for both, as
    /// Ensembl VEP 116.2 prints after the rejoin) and the record's class that of the
    /// first trimmed pair; `GAA GGA,GGA` keeps its repeated allele as a second row
    /// (ALLELE_NUM 1 and 2); `A A,G` drops the allele equal to the reference and
    /// numbers `G` 2; `G T,C` trims nothing and stays as it is, unmarked;
    /// `G GA,*` is chopped to `-/A/*`, which trims nothing more; a gVCF block
    /// before the alleles (`G <NON_REF>,T,GA`, a record Ensembl VEP drops whole)
    /// leaves each allele trimmed against its own part of the string, numbered by
    /// its written position.
    #[test]
    fn minimal_per_alt_trim_multi_allelic() {
        let record = |reference: &str, alts: &str| {
            parse_minimal(&format!("21\t10542414\tx\t{reference}\t{alts}\t.\t.\t."))
        };
        let rows = |variants: &[InputVariant]| -> Vec<(String, String, u32, bool, String)> {
            variants
                .iter()
                .map(|v| {
                    (
                        v.location(),
                        v.allele_string.clone(),
                        v.allele_num.unwrap_or(v.allele_index as u32 + 1),
                        v.minimised,
                        v.class_so_term(),
                    )
                })
                .collect()
        };
        let split = record("G", "T,GA");
        assert_eq!(
            rows(&split),
            [
                ("21:10542414".into(), "G/T".into(), 1, true, "SNV".into()),
                ("21:10542414".into(), "-/A".into(), 2, true, "SNV".into()),
            ]
        );
        assert_eq!((split[1].start, split[1].end), (10542415, 10542414));
        assert_eq!(split[1].rejoin_key_prefix.as_deref(), Some("-/A_"));
        assert_eq!(split[0].rejoin_key_prefix, None);
        assert_eq!(
            rows(&record("GAA", "GGA,GGA")),
            [
                (
                    "21:10542414-10542416".into(),
                    "A/G".into(),
                    1,
                    true,
                    "SNV".into()
                ),
                (
                    "21:10542414-10542416".into(),
                    "A/G".into(),
                    2,
                    true,
                    "SNV".into()
                ),
            ]
        );
        assert_eq!(
            rows(&record("A", "A,G")),
            [(
                "21:10542414".into(),
                "A/G".into(),
                2,
                true,
                "sequence_alteration".into()
            )]
        );
        assert_eq!(
            rows(&record("G", "T,C")),
            [
                ("21:10542414".into(), "G/T".into(), 1, false, "SNV".into()),
                ("21:10542414".into(), "G/C".into(), 2, false, "SNV".into()),
            ]
        );
        assert_eq!(
            rows(&record("G", "GA,*")),
            [(
                "21:10542414-10542415".into(),
                "-/A".into(),
                1,
                false,
                "sequence_alteration".into()
            )]
        );
        let gvcf = record("G", "<NON_REF>,T,GA");
        assert_eq!(
            rows(&gvcf),
            [
                ("21:10542414".into(), "G/T".into(), 2, true, "SNV".into()),
                ("21:10542414".into(), "-/A".into(), 3, true, "SNV".into()),
            ]
        );
        assert_eq!(gvcf[0].rejoin_key_prefix, None);
        assert_eq!(gvcf[1].rejoin_key_prefix.as_deref(), Some("-/A_"));
    }

    /// Under `--process_ref_homs` a split record's leading reference row gives way to
    /// one per allele, each its feature's trimmed reference over the feature's span
    /// and under its rejoin key, before the allele's own row, as Ensembl VEP 116.2
    /// prints with `--minimal --process_ref_homs --allele_number` for
    /// `21 10542414 A T,AG` (`A` 0 and `T` 1 of the feature `A/T`; `-` 0 and `G` 2
    /// of the feature `-/G`, on the insertion span), for `21 10542415 A A,G` (`-` 0
    /// of the feature `-/-` alone, then `A` 0 and `G` 2) and for `21 10542414 A T,T,AG`
    /// (a reference row per repeated allele); a record that does not split keeps
    /// its one row in front.
    #[test]
    fn minimal_split_writes_a_reference_row_per_allele_under_process_ref_homs() {
        let record = |line: &str| {
            split_minimal_record(
                parse_vcf_line_with(
                    line,
                    &VcfParse {
                        minimal: true,
                        reading: RecordReading {
                            individuals: None,
                            process_ref_homs: true,
                        },
                        ..VcfParse::default_run(true)
                    },
                )
                .unwrap(),
            )
        };
        let rows =
            |variants: &[InputVariant]| -> Vec<(String, u32, bool, Option<String>, u64, u64)> {
                variants
                    .iter()
                    .map(|v| {
                        (
                            v.allele_string.clone(),
                            v.allele_num.unwrap_or(v.allele_index as u32 + 1),
                            v.reference_allele_row,
                            v.rejoin_key_prefix.as_deref().map(str::to_string),
                            v.start,
                            v.end,
                        )
                    })
                    .collect()
            };
        let row = |alleles: &str, number, reference, prefix: Option<&str>, start, end| {
            (
                alleles.to_string(),
                number,
                reference,
                prefix.map(str::to_string),
                start,
                end,
            )
        };
        let split = record("21\t10542414\tx\tA\tT,AG\t.\t.\t.");
        assert_eq!(
            rows(&split),
            [
                row("A/A", 0, true, None, 10542414, 10542414),
                row("A/T", 1, false, None, 10542414, 10542414),
                row("-/-", 0, true, Some("-/G_"), 10542415, 10542414),
                row("-/G", 2, false, Some("-/G_"), 10542415, 10542414),
            ]
        );
        assert!(split.iter().all(|v| v.minimised
            && v.uploaded_allele_string.as_deref() == Some("A/T/AG")
            && v.location() == "21:10542414"));
        assert_eq!(
            rows(&record("21\t10542415\tx\tA\tA,G\t.\t.\t.")),
            [
                row("-/-", 0, true, None, 10542416, 10542415),
                row("A/A", 0, true, Some("A/G_"), 10542415, 10542415),
                row("A/G", 2, false, Some("A/G_"), 10542415, 10542415),
            ]
        );
        assert_eq!(
            rows(&record("21\t10542414\tx\tA\tT,T,AG\t.\t.\t.")),
            [
                row("A/A", 0, true, None, 10542414, 10542414),
                row("A/T", 1, false, None, 10542414, 10542414),
                row("A/A", 0, true, Some("A/T_"), 10542414, 10542414),
                row("A/T", 2, false, Some("A/T_"), 10542414, 10542414),
                row("-/-", 0, true, Some("-/G_"), 10542415, 10542414),
                row("-/G", 3, false, Some("-/G_"), 10542415, 10542414),
            ]
        );
        assert_eq!(
            rows(&record("21\t10542414\tx\tG\tT,C\t.\t.\t.")),
            [
                row("G/G", 0, true, None, 10542414, 10542414),
                row("G/T", 1, false, None, 10542414, 10542414),
                row("G/C", 2, false, None, 10542414, 10542414),
            ]
        );
        // `21 10542414 AAAG AGAG,*` is chopped to `AAG/GAG/*`; the feature `AAG/*`
        // writes its reference row `AAG` numbered 0 over the record's span and
        // nothing else, as `<DEL:*>` does, and leads the record when written first.
        let star = record("21\t10542414\tx\tAAAG\tAGAG,*\t.\t.\t.");
        assert_eq!(
            rows(&star),
            [
                row("A/A", 0, true, None, 10542415, 10542415),
                row("A/G", 1, false, None, 10542415, 10542415),
                row("AAG/AAG", 0, true, Some("AAG/*_"), 10542415, 10542417),
            ]
        );
        assert!(star.iter().all(|v| v.minimised
            && v.uploaded_allele_string.as_deref() == Some("AAG/GAG/*")
            && v.location() == "21:10542415-10542417"));
        assert_eq!(
            rows(&record("21\t10542414\tx\tAAAG\tAGAG,<DEL:*>\t.\t.\t.")),
            [
                row("A/A", 0, true, None, 10542415, 10542415),
                row("A/G", 1, false, None, 10542415, 10542415),
                row("AAG/AAG", 0, true, Some("AAG/<DEL:*>_"), 10542415, 10542417),
            ]
        );
        assert_eq!(
            rows(&record("21\t10542414\tx\tAAAG\t*,AGAG\t.\t.\t.")),
            [
                row("AAG/AAG", 0, true, None, 10542415, 10542417),
                row("A/A", 0, true, Some("A/G_"), 10542415, 10542415),
                row("A/G", 2, false, Some("A/G_"), 10542415, 10542415),
            ]
        );
    }

    /// A split allele carries the record's span as its original span and the record's
    /// chopped string as its `UPLOADED_ALLELE` (InputBuffer.pm:694-697, the string
    /// `split_variants` writes over the raw one): `GAAG GGAG,GAGG,G` at 10542414 is
    /// chopped to `AAG/GAG/AGG/-` at 15-17, its three alleles `A/G` at 15, `A/G` at 16
    /// and `AAG/-` at 15-17 (Ensembl VEP 116.2's rows `G` 1, `G` 2 and `-` 3, each
    /// `UPLOADED_ALLELE=AAG/GAG/AGG/-` at `21:10542415-10542417`); a bi-allelic record's
    /// `UPLOADED_ALLELE` stays the raw `GAAG/GGGG`.
    #[test]
    fn minimal_sets_minimised_and_original() {
        let three = parse_minimal("21\t10542414\tx\tGAAG\tGGAG,GAGG,G\t.\t.\t.");
        let spans: Vec<(u64, u64, Option<u64>, Option<u64>)> = three
            .iter()
            .map(|v| (v.start, v.end, v.original_start, v.original_end))
            .collect();
        assert_eq!(
            spans,
            [
                (10542415, 10542415, Some(10542415), Some(10542417)),
                (10542416, 10542416, Some(10542415), Some(10542417)),
                (10542415, 10542417, Some(10542415), Some(10542417)),
            ]
        );
        for v in &three {
            assert!(v.minimised);
            assert_eq!(v.location(), "21:10542415-10542417");
            assert_eq!(v.uploaded_allele(true).as_deref(), Some("AAG/GAG/AGG/-"));
            assert_eq!(
                v.record_allele_string_multi.as_deref(),
                Some("AAG/GAG/AGG/-")
            );
        }
        assert_eq!(
            three
                .iter()
                .map(|v| v.rejoin_key_prefix.as_deref())
                .collect::<Vec<_>>(),
            [None, Some("A/G_"), Some("AAG/-_")]
        );
        let two = parse_minimal("21\t10542414\tx\tGAAG\tGGGG\t.\t.\t.");
        assert_eq!(two[0].uploaded_allele(true).as_deref(), Some("GAAG/GGGG"));
        assert_eq!(two[0].original_start, None);
    }

    /// Guard: the same shape with a second ALT must not be re-minimised.
    ///
    /// `minimise_alleles` returns any allele string matching `.+/.+/.+` unchanged
    /// (`VEP/Parser.pm:978-984`), so a multi-allelic record never reaches
    /// `trim_sequences` and keeps its anchor base.
    #[test]
    fn test_parse_vcf_line_multiallelic_prefix_ref_keeps_anchor_base() {
        let line = "1\t100\t.\tAT\tATT,A\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        let ins = &variants[0];
        assert_eq!(
            ins.ref_allele,
            b"T",
            "a 2-ALT record returns early at Parser.pm:978 and keeps its anchor \
             base; got ref={:?}",
            String::from_utf8_lossy(&ins.ref_allele)
        );
        assert_eq!(ins.alt_alleles[0], b"TT");
        assert_eq!(ins.start, 101);
        assert_eq!(ins.end, 101);
    }

    #[test]
    fn test_parse_ci_field() {
        assert_eq!(parse_ci_field(Some("-100,50")), Some((-100, 50)));
        assert_eq!(parse_ci_field(Some("0,0")), Some((0, 0)));
        assert_eq!(parse_ci_field(Some("-500,500")), Some((-500, 500)));
        assert_eq!(parse_ci_field(None), None);
        assert_eq!(parse_ci_field(Some("")), None);
        assert_eq!(parse_ci_field(Some("abc")), None);
        assert_eq!(parse_ci_field(Some("-100")), None);
    }

    #[test]
    fn test_parse_vcf_line_imprecise_sv_cipos_ciend() {
        let line =
            "21\t27253938\tid1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=27260000;CIPOS=-100,50;CIEND=-50,100";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ci_pos, Some((-100, 50)));
        assert_eq!(v.ci_end, Some((-50, 100)));
        assert!(v.is_structural);
        assert_eq!(v.sv_end, Some(27260000));
    }

    #[test]
    fn test_parse_vcf_line_sv_without_cipos_ciend() {
        let line = "21\t27253938\tid1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=27260000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ci_pos, None);
        assert_eq!(v.ci_end, None);
    }

    #[test]
    fn test_parse_vcf_line_star_ref_block_not_filtered() {
        // Ensembl VEP does not filter <*>: the SV regex explicitly excludes it, so without
        // INFO/SVTYPE it is a non-SV allele at position POS.
        let line = "21\t27253474\trefblock1\tT\t<*>\t.\tPASS\tEND=27255568;LEN=2094";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            1,
            "<*> should produce a variant, not be filtered"
        );
        let v = &variants[0];
        assert_eq!(v.start, 27253474);
        assert_eq!(v.end, 27253474); // single position, not SV span
        assert!(!v.is_structural, "<*> should not be structural");
        assert_eq!(v.alt_allele(), b"<*>");
        assert_eq!(v.allele_string, "T/<*>");
    }

    /// The anchor chop counts a `*` as an allele of one base (Parser/VCF.pm:296-298) and
    /// leaves it out of the first-base comparison (:306), so Ensembl VEP 116.2 chops the
    /// shared first base of a record whose other alleles have the REF's length: at 30425650
    /// it writes `AG AC,*` as `21:30425651 C` (allele string `G/C/*`), `AGC ATT,*` as
    /// `21:30425651-30425652 TT` (`GC/TT/*`), `AG *,AC` as `21:30425651 C` (ALLELE_NUM 2,
    /// `G/*/C`) and `AG AC,AT,*` as `21:30425651 C` and `T` (`G/C/T/*`); the chop needs the
    /// shared base (`AG TC,*` stays `21:30425650-30425651 TC`), a one-base REF has nothing
    /// to chop (`A G,*`, `A *,G`, `A G,T,*` stay at 30425650) and without the `*` the
    /// same-length record is not an indel (`AG AC,AT` stays `21:30425650-30425651`); an
    /// indel beside a `*` is chopped as before (`AG A,*` `21:30425651 -`, `A AG,*`
    /// `21:30425650-30425651 G`, `AG *,A` `21:30425651 -` with ALLELE_NUM 2). A `.` beside the
    /// alleles is compared and never shares the base, so `A AT,.` keeps its anchor
    /// (`21:30425650 AT`, as Ensembl VEP 116.2 writes it beside its own rows under `.`).
    #[test]
    fn a_star_counts_as_an_allele_of_one_base_in_the_anchor_chop() {
        for (reference, alts, allele_string, rows) in [
            ("A", "AT,.", "A/AT/.", vec![("21:30425650", "AT", 0)]),
            (
                "ATT",
                "AT,.",
                "ATT/AT/.",
                vec![("21:30425650-30425652", "AT", 0)],
            ),
            ("A", "G,*", "A/G/*", vec![("21:30425650", "G", 0)]),
            ("AG", "AC,*", "G/C/*", vec![("21:30425651", "C", 0)]),
            (
                "AG",
                "TC,*",
                "AG/TC/*",
                vec![("21:30425650-30425651", "TC", 0)],
            ),
            (
                "AGC",
                "ATT,*",
                "GC/TT/*",
                vec![("21:30425651-30425652", "TT", 0)],
            ),
            ("AG", "A,*", "G/-/*", vec![("21:30425651", "-", 0)]),
            ("A", "AG,*", "-/G/*", vec![("21:30425650-30425651", "G", 0)]),
            ("AG", "*,AC", "G/*/C", vec![("21:30425651", "C", 1)]),
            (
                "AG",
                "AC,AT,*",
                "G/C/T/*",
                vec![("21:30425651", "C", 0), ("21:30425651", "T", 1)],
            ),
            ("A", "*,G", "A/*/G", vec![("21:30425650", "G", 1)]),
            ("AG", "*,A", "G/*/-", vec![("21:30425651", "-", 1)]),
            (
                "A",
                "G,T,*",
                "A/G/T/*",
                vec![("21:30425650", "G", 0), ("21:30425650", "T", 1)],
            ),
            (
                "AG",
                "AC,AT",
                "AG/AC/AT",
                vec![
                    ("21:30425650-30425651", "AC", 0),
                    ("21:30425650-30425651", "AT", 1),
                ],
            ),
        ] {
            let record = format!("21\t30425650\t.\t{reference}\t{alts}\t.\t.\t.");
            let variants = parse_vcf_line(&record, false).unwrap();
            let read: Vec<(String, String, usize)> = variants
                .iter()
                .map(|v| (v.location(), v.display_allele(), v.allele_index))
                .collect();
            let expected: Vec<(String, String, usize)> = rows
                .iter()
                .map(|(location, allele, index)| (location.to_string(), allele.to_string(), *index))
                .collect();
            assert_eq!(read, expected, "{reference} {alts}");
            for v in &variants {
                assert_eq!(
                    v.record_allele_string(),
                    allele_string,
                    "{reference} {alts}"
                );
                assert!(!v.minimised, "{reference} {alts}");
            }
        }
    }

    /// The `*` exemption from the anchor chop belongs to the multi-ALT branch alone
    /// (Parser/VCF.pm:315, inside `if(scalar @$alts > 1)` at :301); the one ALT of a
    /// record is chopped whatever it carries (:332 `substr($alts->[0], 1) || '-'`), so
    /// Ensembl VEP 116.2 reads `21 30289183 x C C*` as the feature `-/*` at
    /// 30289184-30289183, whose one allele is a `*` it writes no row for
    /// (VariationFeatureOverlap.pm:484): no default or tab row, the VCF line carried
    /// without consequences, the JSON object with `allele_string` `-/*` and
    /// `most_severe_consequence` `?`. Beside another ALT the same `C*` is kept whole
    /// (`C CT,C*` reads `-/T/C*`) and annotated, as `C*` is not `*` at :484.
    #[test]
    fn the_one_alt_of_a_record_is_chopped_whatever_it_carries() {
        assert_eq!(
            trim_multi_alleles("C", &["C*"]),
            ("-".to_string(), vec!["*".to_string()], 1)
        );
        let single = parse_vcf_line("21\t30289183\tx\tC\tC*\t.\t.\t.", false).unwrap();
        assert_eq!(single.len(), 1, "{single:?}");
        assert!(single[0].no_allele, "{single:?}");
        assert_eq!(single[0].allele_string, "-/*");
        assert_eq!((single[0].start, single[0].end), (30289184, 30289183));
        assert_eq!(single[0].location(), "21:30289183-30289184");

        assert_eq!(
            trim_multi_alleles("C", &["CT", "C*"]),
            ("-".to_string(), vec!["T".to_string(), "C*".to_string()], 1)
        );
        let beside = parse_vcf_line("21\t30289183\tx\tC\tCT,C*\t.\t.\t.", false).unwrap();
        let read: Vec<(String, String, usize)> = beside
            .iter()
            .map(|v| (v.location(), v.display_allele(), v.allele_index))
            .collect();
        assert_eq!(
            read,
            [
                ("21:30289183-30289184".to_string(), "T".to_string(), 0),
                ("21:30289183-30289184".to_string(), "C*".to_string(), 1),
            ]
        );
        assert!(beside.iter().all(|v| v.record_allele_string() == "-/T/C*"));
    }

    /// A `*` the minimisation leaves writes no row either: Ensembl VEP builds the
    /// overlap alleles from the minimised string (VariationFeatureOverlap.pm:484
    /// runs after `post_process_vfs`, Parser.pm:881), so `21 30289183 x C *C`,
    /// whose alleles differ in length and are minimised whatever the flags
    /// (`:876-881`), is the feature `-/*` at 30289183-30289182, and `CA C*` under
    /// `--minimal` the feature `A/*` at 30289184, each marked minimised
    /// (`:1011`), with no row. Without the flag `CA C*` keeps its lengths and
    /// annotates the allele `C*`, as the unminimised string names no `*`.
    #[test]
    fn a_star_left_by_the_minimisation_writes_no_row() {
        let suffix = parse_vcf_line("21\t30289183\tx\tC\t*C\t.\t.\t.", false).unwrap();
        assert_eq!(suffix.len(), 1, "{suffix:?}");
        assert!(suffix[0].no_allele && suffix[0].minimised, "{suffix:?}");
        assert_eq!(suffix[0].allele_string, "-/*");
        assert_eq!((suffix[0].start, suffix[0].end), (30289183, 30289182));

        let minimal = parse_minimal("21\t30289183\tx\tCA\tC*\t.\t.\t.");
        assert_eq!(minimal.len(), 1, "{minimal:?}");
        assert!(minimal[0].no_allele && minimal[0].minimised, "{minimal:?}");
        assert_eq!(minimal[0].allele_string, "A/*");
        assert_eq!((minimal[0].start, minimal[0].end), (30289184, 30289184));

        let plain = parse_vcf_line("21\t30289183\tx\tCA\tC*\t.\t.\t.", false).unwrap();
        assert_eq!(plain.len(), 1, "{plain:?}");
        assert!(!plain[0].no_allele && !plain[0].minimised, "{plain:?}");
        assert_eq!(plain[0].allele_string, "CA/C*");
    }

    /// A `*` beside an indel keeps the record multi-allelic: `minimise_alleles`
    /// (Parser.pm:978-984) leaves a three-part string as it is, so Ensembl VEP 116.2 writes
    /// `21 30425650 x ATT AT,*` at `21:30425651-30425652` with the Allele `T` (the
    /// anchor-chopped `TT/T/*`), `AT ATT,*` at `21:30425651` with `TT`, and `ATTG ATC,*` at
    /// `21:30425651-30425653` with `TC`; without the `*` the bi-allelic indel is minimised
    /// from both ends (`21:30425652 -`, `21:30425651-30425652 T`, `21:30425652-30425653 C`).
    #[test]
    fn a_star_beside_an_indel_keeps_the_anchor_chopped_alleles() {
        let read = |record: &str| {
            let v = parse_vcf_line(record, false).unwrap();
            assert_eq!(v.len(), 1, "{record}");
            (
                v[0].location(),
                v[0].display_allele(),
                v[0].record_allele_string(),
                v[0].minimised,
            )
        };
        for (record, location, allele, string) in [
            (
                "21\t30425650\tx\tATT\tAT,*\t.\t.\t.",
                "21:30425651-30425652",
                "T",
                "TT/T/*",
            ),
            (
                "21\t30425650\tx\tAT\tATT,*\t.\t.\t.",
                "21:30425651",
                "TT",
                "T/TT/*",
            ),
            (
                "21\t30425650\tx\tATTG\tATC,*\t.\t.\t.",
                "21:30425651-30425653",
                "TC",
                "TTG/TC/*",
            ),
        ] {
            let (loc, al, st, minimised) = read(record);
            assert_eq!(
                (loc.as_str(), al.as_str(), st.as_str()),
                (location, allele, string),
                "{record}"
            );
            assert!(!minimised, "{record}");
        }
        for (record, location, allele) in [
            ("21\t30425650\tx\tATT\tAT\t.\t.\t.", "21:30425652", "-"),
            (
                "21\t30425650\tx\tAT\tATT\t.\t.\t.",
                "21:30425651-30425652",
                "T",
            ),
            (
                "21\t30425650\tx\tATTG\tATC\t.\t.\t.",
                "21:30425652-30425653",
                "C",
            ),
        ] {
            let (loc, al, _, minimised) = read(record);
            assert_eq!((loc.as_str(), al.as_str()), (location, allele), "{record}");
            assert!(minimised, "{record}");
        }
    }

    /// A line whose every ALT is `*` names no allele: Ensembl VEP 116.2 writes
    /// `21 30425650 x N * . . .` and `N *,*` with no row, the VCF line without CSQ and a
    /// JSON object (`"allele_string":"N/*"`, start and end 30425650,
    /// `"most_severe_consequence":"?"`). With `SVTYPE=DEL;END=400` Ensembl VEP routes
    /// the record to its structural path (Parser/VCF.pm:237-245) and annotates a
    /// deletion from INFO alone, which vep-rs does not: the record is read allele-less
    /// whatever INFO says (`docs/divergences/defect-m.md`). The anchor chop counts the
    /// REF alone in its first-base vote when every ALT carries a `*`
    /// (Parser/VCF.pm:306), so `21 30396602 x TTGA *,*` is the feature `TGA/*/*` at
    /// 30396603-30396605 and `TTGA *` stays `TTGA/*` at 30396602-30396605. The one
    /// variant carries the record and is never annotated; a `*` beside an allele on
    /// the sequence path stays filtered.
    #[test]
    fn a_star_only_record_is_carried_as_the_allele_less_feature() {
        for (line, allele_string, span) in [
            ("21\t100\t.\tA\t*\t.\tPASS\t.", "A/*", (100, 100)),
            ("21\t100\t.\tNAC\t*,*\t.\tPASS\t.", "AC/*/*", (101, 102)),
            ("21\t100\t.\tNAC\t*\t.\tPASS\t.", "NAC/*", (100, 102)),
            (
                "21\t100\t.\tN\t*\t.\tPASS\tSVTYPE=DEL;END=400",
                "N/*",
                (100, 100),
            ),
        ] {
            let variants = parse_vcf_line(line, true).unwrap();
            assert_eq!(variants.len(), 1, "{line}");
            let v = &variants[0];
            assert!(v.no_allele, "{line}");
            assert_eq!(v.allele_string, allele_string, "{line}");
            assert_eq!((v.start, v.end), span, "{line}");
            assert!(!v.is_structural && v.transcript_consequences.is_empty());
            assert_eq!(v.record_allele_string(), allele_string, "{line}");
        }
        assert_eq!(
            parse_vcf_line("21\t100\t.\tA\t*\t.\tPASS\t.", true).unwrap()[0].uploaded_variation(),
            "21_100_A/*"
        );
        let beside = parse_vcf_line("21\t100\t.\tA\tG,*\t.\tPASS\t.", true).unwrap();
        assert_eq!(beside.len(), 1);
        assert!(!beside[0].no_allele);
        assert_eq!(beside[0].alt_allele(), b"G");
    }

    /// A `*` beside a real allele under a true INFO/SVTYPE is routed with it to the
    /// structural path (`CG*` is no run of `ACGT`, Parser/VCF.pm:237-245), where the
    /// record's allele string `G/*` (`:469`) gives one allele per part with no `*`
    /// filter (StructuralVariationOverlap.pm:64-77): Ensembl VEP 116.2 writes two
    /// `deletion` rows per transcript for `21 30289183 x C G,* . . SVTYPE=DEL;END=30289190`,
    /// `ALLELE_NUM=1` and `ALLELE_NUM=2`, at 21:30289184-30289190, and the same two for
    /// `C *,G`; the JSON `allele_string` is `G/*` and `*/G`. A `.` beside a real allele
    /// is routed the same way (`CG.` is no run of `ACGT`) and is a part like any other:
    /// two `deletion` rows for `C G,.`, `C .,G` and `TC .,<DEL>` under the same INFO.
    #[test]
    fn a_star_beside_an_alt_under_svtype_is_a_structural_allele() {
        for (alts, allele_string) in [
            ("G,*", "G/*"),
            ("*,G", "*/G"),
            ("G,.", "G/."),
            (".,G", "./G"),
        ] {
            let line = format!("21\t30289183\tx\tC\t{alts}\t.\t.\tSVTYPE=DEL;END=30289190");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 2, "{alts}");
            for (index, variant) in variants.iter().enumerate() {
                assert!(variant.is_structural && !variant.no_allele, "{alts}");
                assert_eq!(variant.allele_index, index, "{alts}");
                assert_eq!(variant.display_allele(), "deletion", "{alts}");
                assert_eq!((variant.start, variant.end), (30289184, 30289190), "{alts}");
                assert_eq!(variant.record_allele_string(), allele_string, "{alts}");
                assert_eq!(variant.raw_input.as_deref(), Some(line.as_str()), "{alts}");
            }
            let written: Vec<&[u8]> = variants.iter().map(|v| v.alt_allele()).collect();
            assert_eq!(written.join(&b","[..]), alts.as_bytes(), "{alts}");
        }
    }

    #[test]
    fn test_parse_vcf_line_parid_as_mate_id() {
        // VCF 4.5 uses PARID instead of MATEID for BND partner identification;
        // PARID is read when MATEID is absent.
        let line = "21\t100\tbnd_1\tA\tA[21:200[\t.\t.\tSVTYPE=BND;PARID=bnd_2";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        let paired = &variants[1];
        assert_eq!(paired.mate_id.as_deref(), Some("bnd_2"));
        let single = &variants[0];
        assert_eq!(single.mate_id.as_deref(), Some("bnd_2"));
        assert!(single.is_single_breakend);
        assert_eq!(single.display_allele(), "A.");
    }

    #[test]
    fn test_parse_vcf_line_mateid_preferred_over_parid() {
        // When both MATEID and PARID are present, MATEID takes precedence
        let line = "21\t100\tbnd_1\tA\tA[21:200[\t.\t.\tSVTYPE=BND;MATEID=mate1;PARID=par1";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants[0].mate_id.as_deref(), Some("mate1"));
    }

    #[test]
    fn test_gnomad_grch38_del_positive_svlen() {
        // gnomAD v4.1 GRCh38 DEL: chr21, positive SVLEN, with END and CHR2
        // Real example: chr21:5033802 SVLEN=157 END=5033959
        let line = "chr21\t5033802\tgnomAD-SV_v3_DEL_chr21_test\tN\t<DEL>\t.\t.\tALGORITHMS=manta;CHR2=chr21;END=5033959;SVLEN=157;SVTYPE=DEL";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.chr, "21");
        assert_eq!(v.original_chr, "chr21");
        // sv_start = POS + 1 (DEL is span_type)
        assert_eq!(v.start, 5033803);
        // sv_end = POS + unsigned_abs(SVLEN) = 5033802 + 157 = 5033959
        assert_eq!(v.end, 5033959);
        assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
        assert!(v.is_structural);
        assert_eq!(v.location(), "chr21:5033803-5033959");
        assert_eq!(v.display_allele(), "deletion");
    }

    #[test]
    fn test_gnomad_grch38_dup() {
        let line = "chr21\t10000000\tgnomAD_DUP_test\tN\t<DUP>\t.\t.\tCHR2=chr21;END=10000500;SVLEN=500;SVTYPE=DUP";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.chr, "21");
        assert_eq!(v.original_chr, "chr21");
        assert_eq!(v.start, 10000001);
        assert_eq!(v.end, 10000500);
        assert_eq!(v.variant_class, VariantClass::Duplication);
        assert_eq!(v.location(), "chr21:10000001-10000500");
        assert_eq!(v.display_allele(), "duplication");
    }

    #[test]
    fn test_gnomad_grch38_ins_me_alu() {
        let line = "chr21\t20000000\tgnomAD_ME_ALU_test\tN\t<INS:ME:ALU>\t.\t.\tCHR2=chr21;END=20000001;SVLEN=300;SVTYPE=INS:ME:ALU";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.chr, "21");
        // ME insertions are point-type symbolic alleles; `is_span_type()` is true for every SV class.
        assert_eq!(v.variant_class, VariantClass::MobileElementInsertion);
        assert_eq!(v.display_allele(), "Alu_insertion");
    }

    #[test]
    fn test_gnomad_grch38_bnd_svlen_minus1() {
        // BND with SVLEN=-1 and CHR2+END2: CHR2 with END2 names the mate whatever
        // the span, so the record yields its own breakend and the mate even though
        // SVLEN leaves end == start.
        let line = "chr21\t10838490\tgnomAD_BND_test\tN\t<BND>\t.\t.\tCHR2=chr2;END=10838491;END2=89831875;SVLEN=-1;SVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "BND SVLEN=-1 with CHR2+END2: should produce paired + single-breakend"
        );
        let paired = &variants[1];
        assert_eq!(paired.chr, "21");
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert!(!paired.is_single_breakend);
        assert_eq!(paired.mate_chr.as_deref(), Some("chr2"));
        assert_eq!(paired.mate_pos, Some(89831875));
        let sb = &variants[0];
        assert_eq!(sb.chr, "21");
        assert_eq!(sb.variant_class, VariantClass::Translocation);
        assert!(sb.is_single_breakend);
        assert!(sb.mate_chr.is_none());
    }

    #[test]
    fn test_gnomad_grch38_bnd_large_svlen_with_chr2_end2() {
        let line = "chr21\t10838490\tgnomAD_BND_large_test\tN\t<BND>\t.\t.\tCHR2=chr2;END=10873831;END2=89831875;SVLEN=35340;SVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        // Large SVLEN: sv_start = 10838491, sv_end = 10838491 + 35340 - 1 = 10873830
        // CHR2 with END2 names the mate -> 2 variants
        assert_eq!(variants.len(), 2);
        let paired = &variants[1];
        // The mate is carried as written; `pipeline::name_source_chromosomes`
        // resolves it against the cache.
        assert_eq!(paired.mate_chr.as_deref(), Some("chr2"));
        assert_eq!(paired.mate_pos, Some(89831875));
        let single = &variants[0];
        assert!(single.is_single_breakend);
    }

    #[test]
    fn test_gnomad_grch38_ctx() {
        // <CTX> has no explicit handling in classify_symbolic_alt and SVTYPE=CTX
        // is in no match arm, so it lands in ComplexStructural.
        let line = "chr21\t10000000\tgnomAD_CTX_test\tN\t<CTX>\t.\t.\tCHR2=chr3;END=10000001;END2=50000000;SVLEN=-1;SVTYPE=CTX";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::ComplexStructural);
        // Ensembl VEP treats CTX as a translocation/BND; vep-rs does not.
    }

    /// The per-individual reading of `Parser/VCF.pm:595-665` on `21 100 x A G,T,*`
    /// with eight samples, one genotype shape each, against ensembl-io's reading
    /// of the GT field (`BaseVCF4.pm:944-951`): `0/1` is the one allele G
    /// with `HET` (two distinct bits, A and G); `1/1` is G with `HOM`; `1/2` is the
    /// alleles G and T, each `HET`, numbered 1 and 2 along the copy's string
    /// `A/G/T` (the `*` outside it, `:635`); `./1` drops its `.` bit and is the one
    /// bit G, so `HOM`; `./.` is no genotype and the individual has no copy
    /// (`:626`); a haploid `1` is G with `HOM`; `0/3`, the `*` allele beside the
    /// reference, carries no alternate allele and is a reference-homozygous copy
    /// (`:641`) with no allele to annotate and `HETREF`, the two distinct bits A
    /// and `*` (`OutputFactory.pm:976`); `0/0` is a reference-homozygous copy with
    /// `HOMREF`. The copies come in header order, each carrying its sample's name
    /// (`:660`) and genotype (`:657`).
    #[test]
    fn individual_all_emits_one_variant_set_per_individual() {
        let samples = genotyped_samples(&["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]);
        let split = IndividualSplit::new(&IndividualRequest::All, &samples, IndividualMode::Copies)
            .unwrap();
        let reading = RecordReading {
            individuals: Some(&split),
            process_ref_homs: false,
        };
        let line = "21\t100\tx\tA\tG,T,*\t.\t.\t.\tGT\t0/1\t1/1\t1/2\t./1\t./.\t1\t0/3\t0/0";
        let variants = parse_vcf_record(line, reading).unwrap();
        let summary: Vec<(String, String, String, Option<u32>, bool)> = variants
            .iter()
            .map(|v| {
                (
                    v.individual.as_deref().unwrap().to_string(),
                    if v.no_allele {
                        "-".to_string()
                    } else {
                        v.display_allele()
                    },
                    v.zygosity().unwrap(),
                    v.allele_num,
                    v.no_allele,
                )
            })
            .collect();
        let row = |ind: &str, allele: &str, zyg: &str, num: Option<u32>, none: bool| {
            (
                ind.to_string(),
                allele.to_string(),
                zyg.to_string(),
                num,
                none,
            )
        };
        assert_eq!(
            summary,
            [
                row("S1", "G", "HET", Some(1), false),
                row("S2", "G", "HOM", Some(1), false),
                row("S3", "G", "HET", Some(1), false),
                row("S3", "T", "HET", Some(2), false),
                row("S4", "G", "HOM", Some(1), false),
                row("S6", "G", "HOM", Some(1), false),
                row("S7", "-", "HETREF", None, true),
                row("S8", "-", "HOMREF", None, true),
            ]
        );
        assert_eq!(variants[0].genotype, ["A", "G"]);
        assert_eq!(variants[3].genotype, ["G", "T"]);
        assert_eq!(variants[4].genotype, ["G"]);
        assert_eq!(variants[6].genotype, ["A", "*"]);
        assert!(variants[6].hom_ref && variants[7].hom_ref && !variants[0].hom_ref);
        // The copy's own allele string is what its JSON `allele_string`, its
        // generated name and VARIANT_CLASS read; the record's string keeps the `*`.
        assert_eq!(variants[0].record_allele_string(), "A/G");
        assert_eq!(variants[2].record_allele_string(), "A/G/T");
        assert_eq!(variants[6].record_allele_string(), "A");
        assert_eq!(variants[0].uploaded_variation(), "x");
        // Each copy is a record of its own to the pipeline; the alleles of one
        // copy belong together.
        assert!(variants[2].same_record(&variants[3]));
        assert!(!variants[1].same_record(&variants[2]));
        // Without the flag the record is the one set of alleles, no sample read.
        let plain = parse_vcf_record(line, RecordReading::default()).unwrap();
        assert_eq!(plain.len(), 2);
        assert!(plain
            .iter()
            .all(|v| v.individual.is_none() && v.genotype.is_empty()));
    }

    /// The named list is read in its order (`Parser/VCF.pm:609`, `$include` is the
    /// list as given), and a name the header lacks is the `die` of `:612-617`, both
    /// lists in the message; a VCF without a `#CHROM` line names no sample, so a
    /// list is refused against an empty header and `all` reads nothing.
    #[test]
    fn individual_list_unknown_sample_dies_with_vep_text() {
        let samples = genotyped_samples(&["S1", "S2", "S3"]);
        let named = |names: &[&str]| {
            IndividualRequest::Named(names.iter().map(|s| s.to_string()).collect())
        };
        assert_eq!(
            IndividualSplit::new(&named(&["S9"]), &samples, IndividualMode::Copies).err(),
            Some("Sample IDs given (S9) do not match samples from VCF (S1,S2,S3)".to_string())
        );
        assert_eq!(
            IndividualSplit::new(
                &named(&["S1", "S9"]),
                &SampleColumns::default(),
                IndividualMode::Copies
            )
            .err(),
            Some("Sample IDs given (S1,S9) do not match samples from VCF ()".to_string())
        );
        let split =
            IndividualSplit::new(&named(&["S3", "S1"]), &samples, IndividualMode::Copies).unwrap();
        let reading = RecordReading {
            individuals: Some(&split),
            process_ref_homs: false,
        };
        let line = "21\t100\tx\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t1/1";
        let variants = parse_vcf_record(line, reading).unwrap();
        let order: Vec<&str> = variants
            .iter()
            .map(|v| v.individual.as_deref().unwrap())
            .collect();
        assert_eq!(order, ["S3", "S1"]);
        let all = IndividualSplit::new(
            &IndividualRequest::All,
            &SampleColumns::default(),
            IndividualMode::Copies,
        )
        .unwrap();
        let reading = RecordReading {
            individuals: Some(&all),
            process_ref_homs: false,
        };
        assert!(parse_vcf_record(line, reading).unwrap().is_empty());
    }

    /// `--process_ref_homs` (`Parser/VCF.pm:643-650`, `OutputFactory.pm:505-516`): a
    /// reference-homozygous individual's copy is the string `REF/REF`, whose one
    /// allele is the reference, so the copy is the reference row alone, numbered 0
    /// with `HOMREF`; a carrier's copy leads with the reference row before its
    /// alleles (`BaseVariationFeatureOverlap.pm:363-373`); a record read without
    /// `--individual` gains the row before its alleles, and a non-variant record
    /// (ALT `.`) is the reference row alone.
    #[test]
    fn process_ref_homs_keeps_hom_ref_individuals_as_homref_rows() {
        let samples = genotyped_samples(&["S1", "S2"]);
        let split = IndividualSplit::new(&IndividualRequest::All, &samples, IndividualMode::Copies)
            .unwrap();
        let reading = RecordReading {
            individuals: Some(&split),
            process_ref_homs: true,
        };
        let line = "21\t100\tx\tA\tG\t.\t.\t.\tGT\t0/0\t0/1";
        let variants = parse_vcf_record(line, reading).unwrap();
        assert_eq!(variants.len(), 3);
        let hom_ref = &variants[0];
        assert!(hom_ref.reference_allele_row && hom_ref.hom_ref && !hom_ref.no_allele);
        assert_eq!(hom_ref.individual.as_deref(), Some("S1"));
        assert_eq!(hom_ref.alt_allele(), b"A");
        assert_eq!(hom_ref.allele_num, Some(0));
        assert_eq!(hom_ref.zygosity().as_deref(), Some("HOMREF"));
        assert_eq!(hom_ref.record_allele_string(), "A/A");
        assert_eq!(hom_ref.allele_string, "A/A");
        let (reference, alternate) = (&variants[1], &variants[2]);
        assert!(reference.reference_allele_row && !alternate.reference_allele_row);
        assert_eq!(reference.individual.as_deref(), Some("S2"));
        assert_eq!(reference.allele_num, Some(0));
        assert_eq!(alternate.allele_num, Some(1));
        assert_eq!(reference.zygosity().as_deref(), Some("HET"));
        assert_eq!(reference.record_allele_string(), "A/G");
        assert!(reference.same_record(alternate));

        let reading = RecordReading {
            individuals: None,
            process_ref_homs: true,
        };
        let record = parse_vcf_record("21\t100\tx\tA\tG,T\t.\t.\t.", reading).unwrap();
        let alleles: Vec<&[u8]> = record.iter().map(|v| v.alt_allele()).collect();
        assert_eq!(alleles, [b"A".as_slice(), b"G", b"T"]);
        assert!(record[0].reference_allele_row && record[0].individual.is_none());
        assert_eq!(record[0].allele_num, Some(0));
        assert_eq!(record[0].record_allele_string(), "A/G/T");
        let non_variant = non_variant_record_read("21\t100\tx\tA\t.\t.\t.\t.", reading).unwrap();
        assert_eq!(non_variant.len(), 1);
        assert!(non_variant[0].reference_allele_row && !non_variant[0].no_allele);
        assert_eq!(non_variant[0].record_allele_string(), "A");
        // A structural record is built before either step (Parser/VCF.pm:244).
        let structural = parse_vcf_record(
            "21\t100\tsv\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=200\tGT\t0/1\t0/0",
            RecordReading {
                individuals: Some(&split),
                process_ref_homs: true,
            },
        )
        .unwrap();
        assert_eq!(structural.len(), 1);
        assert!(structural[0].individual.is_none() && !structural[0].reference_allele_row);
        // A `<CNV:TR>` expansion is built there too (`:543`), so it is read once,
        // and gains the reference row at output like any variation feature.
        let (_dir, fasta) = fasta_of(RUN_FASTA);
        let expansion = parse_vcf_line_with(
            "21\t10\ttr\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=4\tGT\t0/1\t0/0",
            &VcfParse {
                fasta: Some(&fasta),
                reading: RecordReading {
                    individuals: Some(&split),
                    process_ref_homs: true,
                },
                ..VcfParse::default_run(false)
            },
        )
        .unwrap();
        assert_eq!(expansion.len(), 2);
        assert!(expansion.iter().all(|v| v.individual.is_none()));
        assert!(expansion[0].reference_allele_row && !expansion[1].reference_allele_row);
        assert_eq!(expansion[0].allele_num, Some(0));
        assert_eq!(expansion[1].alt_allele(), b"TCA");
    }

    /// A non-variant record under `--individual` is copied like any other
    /// (`Parser/VCF.pm:355-376` runs after `:259-262` kept it): a `0/0` sample is a
    /// reference-homozygous copy, a haploid `0` too, and a missing genotype none; a
    /// record every individual misses yields nothing (`:371-374`); a repeated ALT
    /// maps every genotype naming it to the record's one row for that allele
    /// (`:635` keys the alleles by string), numbered along the copy's string.
    #[test]
    fn individual_copies_of_non_variant_and_repeated_alt_records() {
        let samples = genotyped_samples(&["S1", "S2", "S3"]);
        let split = IndividualSplit::new(&IndividualRequest::All, &samples, IndividualMode::Copies)
            .unwrap();
        let reading = RecordReading {
            individuals: Some(&split),
            process_ref_homs: false,
        };
        let copies =
            non_variant_record_read("21\t100\tx\tA\t.\t.\t.\t.\tGT\t0/0\t./.\t0", reading).unwrap();
        let names: Vec<&str> = copies
            .iter()
            .map(|v| v.individual.as_deref().unwrap())
            .collect();
        assert_eq!(names, ["S1", "S3"]);
        assert!(copies.iter().all(|v| v.no_allele && v.hom_ref));
        assert_eq!(copies[0].zygosity().as_deref(), Some("HOMREF"));
        assert_eq!(copies[0].record_allele_string(), "A");
        assert!(
            parse_vcf_record("21\t100\tx\tA\tG\t.\t.\t.\tGT\t./.\t./.\t.", reading)
                .unwrap()
                .is_empty()
        );
        let repeated =
            parse_vcf_record("21\t100\tx\tC\tG,G\t.\t.\t.\tGT\t0/2\t1/2\t0/0", reading).unwrap();
        assert_eq!(repeated.len(), 3);
        assert_eq!(repeated[0].alt_allele(), b"G");
        assert_eq!(repeated[0].allele_num, Some(1));
        assert_eq!(repeated[0].zygosity().as_deref(), Some("HET"));
        assert_eq!(repeated[0].record_allele_string(), "C/G");
        assert_eq!(repeated[1].zygosity().as_deref(), Some("HOM"));
        assert_eq!(repeated[1].genotype, ["G", "G"]);
        assert!(repeated[2].no_allele);
    }

    /// The record read once ([`parse_data_line`]) is the record the line's own
    /// reader builds, a non-variant line's included, with the tandem-repeat
    /// warnings Ensembl VEP writes before `validate_vf` (Parser/VCF.pm:502-527):
    /// written for a tandem record whose POS the parser then refuses, and not for a
    /// sequence record, a non-variant record or a line of fewer than five columns,
    /// whose parse error is the line's whole reading.
    #[test]
    fn the_record_read_once_is_the_readers_record_with_its_warnings() {
        let parse = VcfParse::default_run(true);
        let alleles = |variants: &[InputVariant]| -> Vec<(String, String, bool)> {
            variants
                .iter()
                .map(|v| (v.location(), v.allele_string.clone(), v.no_allele))
                .collect()
        };
        let tandem = "21\t500\tx\tG\t<CNV:TR>,<CNV:TR>\t.\t.\tSVLEN=9,6;RN=1,1;RUS=TCA,TCA;RUC=4,2";
        let record = parse_data_line(tandem, &parse);
        assert_eq!(
            alleles(&record.alleles.unwrap()),
            alleles(&parse_vcf_line_with(tandem, &parse).unwrap())
        );
        assert_eq!(record.warnings, line_warnings(tandem, false));
        assert_eq!(record.warnings.len(), 2);
        let info = record.info.as_ref().unwrap();
        assert_eq!(info.value_written("SVLEN"), Some("9,6"));
        assert_eq!(info.value_written("END"), None);
        let copied = parse_data_line(tandem, &parse).with_info_copied(tandem);
        let info = copied.info.as_ref().unwrap();
        assert!(
            matches!(info, RecordInfo::Copied(column) if column == "SVLEN=9,6;RN=1,1;RUS=TCA,TCA;RUC=4,2")
        );
        assert_eq!(info.value_written("SVLEN"), Some("9,6"));
        assert_eq!(info.value_written("END"), None);
        let refused = parse_data_line(
            "21\tbad\tx\tG\t<CNV:TR>\t.\t.\tSVLEN=9;RN=1;RUS=TCA;RUC=4",
            &parse,
        );
        assert_eq!(
            refused.alleles.unwrap_err().to_string(),
            "Invalid POS: 'bad'"
        );
        assert_eq!(refused.warnings.len(), 1);
        for line in [
            "21\t100\tx\tA\tG,T\t.\t.\tSVTYPE=CNV:TR;RUS=TCA;RUC=4",
            "21\t100\tx\tTG\t.,TA\t.\t.\t.",
            "21\t100\tx\tTG\t.,TA\t.\t.\tSVTYPE=DEL;END=200",
        ] {
            let record = parse_data_line(line, &parse);
            assert!(record.warnings.is_empty(), "{line}");
            let own = if vcf_line_is_non_variant(line) {
                non_variant_record(line, &parse)
            } else {
                parse_vcf_line_with(line, &parse)
            };
            assert_eq!(
                alleles(&record.alleles.unwrap()),
                alleles(&own.unwrap()),
                "{line}"
            );
        }
        let short = "21\t100\tx\tTG";
        let record = parse_data_line(short, &parse);
        assert!(record.warnings.is_empty());
        assert!(record.info.is_none());
        assert_eq!(
            record.alleles.unwrap_err().to_string(),
            "VCF line has fewer than 5 fields"
        );
        assert!(parse_data_line(short, &parse)
            .with_info_copied(short)
            .info
            .is_none());
    }

    /// `Parser/VCF.pm:259` tests the first comma-separated ALT alone (ensembl-io
    /// `BaseVCF4.pm:489` splits the field on commas, Perl's `split` dropping a
    /// trailing empty field), so `.,TA` and `.,` are non-variant records where
    /// `TA,.` is a record with the alleles `TA` and `.`; an empty or missing ALT
    /// is no non-variant. The structural test of `:237-245` comes first, so
    /// `.,<DEL>` and `.,TA` under a true `SVTYPE` are structural records, while
    /// a lone `.` under one keeps its non-variant reading here. The record keeps
    /// the reference as its whole allele string over the reference's span
    /// (`:342`, `:254`), and a genotype indexes the reference and every written
    /// ALT: `0/2` on `.,TA` names `TA`, whose copy carries that allele (`:627-639`
    /// build the copy's allele string from the alleles named), while `0/0` is a
    /// reference-homozygous copy with none; past the list (`2/0` on `.,`) a
    /// genotype names nothing, an undefined allele Ensembl VEP keeps as not the
    /// reference (`:635`).
    #[test]
    fn non_variant_is_the_first_alt_alone() {
        assert!(vcf_line_is_non_variant("21\t100\tx\tTG\t.\t.\t.\t."));
        assert!(vcf_line_is_non_variant("21\t100\tx\tTG\t.,TA\t.\t.\t."));
        assert!(vcf_line_is_non_variant("21\t100\tx\tTG\t.,\t.\t.\t."));
        assert!(vcf_line_is_non_variant(
            "21\t100\tx\tTG\t.\t.\t.\tSVTYPE=DEL;END=200"
        ));
        assert!(!vcf_line_is_non_variant(
            "21\t100\tx\tTG\t.,<DEL>\t.\t.\tEND=200"
        ));
        assert!(!vcf_line_is_non_variant(
            "21\t100\tx\tTG\t.,TA\t.\t.\tSVTYPE=DEL;END=200"
        ));
        assert!(!vcf_line_is_non_variant("21\t100\tx\tTG\tTA,.\t.\t.\t."));
        assert!(!vcf_line_is_non_variant("21\t100\tx\tTG\t,.\t.\t.\t."));
        assert!(!vcf_line_is_non_variant("21\t100\tx\tTG\t\t.\t.\t."));
        assert!(!vcf_line_is_non_variant("21\t100\tx\tTG"));
        let reading = RecordReading {
            individuals: None,
            process_ref_homs: false,
        };
        let record = non_variant_record_read("21\t100\tx\tTG\t.,TA\t.\t.\t.", reading).unwrap();
        assert_eq!(record.len(), 1);
        assert!(record[0].no_allele);
        assert_eq!(record[0].record_allele_string(), "TG");
        assert_eq!((record[0].start, record[0].end), (100, 101));

        let samples = genotyped_samples(&["S1", "S2"]);
        let split = IndividualSplit::new(&IndividualRequest::All, &samples, IndividualMode::Copies)
            .unwrap();
        let per_individual = RecordReading {
            individuals: Some(&split),
            process_ref_homs: false,
        };
        let copies = non_variant_record_read(
            "21\t100\tx\tTG\t.,TA\t.\t.\t.\tGT\t0/0\t0/2",
            per_individual,
        )
        .unwrap();
        let names: Vec<&str> = copies
            .iter()
            .map(|v| v.individual.as_deref().unwrap())
            .collect();
        assert_eq!(names, ["S1", "S2"]);
        assert!(copies[0].hom_ref && copies[0].no_allele && copies[0].genotype == ["TG", "TG"]);
        assert!(!copies[1].hom_ref && !copies[1].no_allele);
        assert_eq!(copies[1].alt_allele(), b"TA");
        assert_eq!(copies[1].genotype, ["TG", "TA"]);
        assert_eq!((copies[1].start, copies[1].end), (100, 101));
        let past_the_list =
            non_variant_record_read("21\t100\tx\tTG\t.,\t.\t.\t.\tGT\t0/0\t2/0", per_individual)
                .unwrap();
        assert_eq!(past_the_list.len(), 2);
        assert!(!past_the_list[1].hom_ref && past_the_list[1].genotype == ["", "TG"]);
    }

    /// A copy of every allele of a bi-allelic record keeps the record's string as
    /// `post_process_vfs` left it (Parser.pm:866-881 minimises an indel whose
    /// alleles differ in length): `CTA GA`, which the chop leaves alone and the
    /// minimisation reads as `CT/G` over two bases, is `CT/G` to its carrier as
    /// to a run without `--individual`, `CT/CT` to a reference-homozygous
    /// individual (`Parser/VCF.pm:604`, `:642`, `$ref` of the post-processed
    /// feature), and every reference row under `--process_ref_homs` is `CT` over
    /// that span, the allele `validate_vf`'s length check admits. Ensembl VEP
    /// itself compares the chopped bits `CTA` and `GA` with the minimised `CT`
    /// (`:635-639`) and gives the carrier `CT/CTA/GA` and the `0/0` individual
    /// `CT/CTA`: the divergence `docs/intended-divergences.md` lists under the
    /// trimmable indel, which this test pins on the vep-rs side.
    #[test]
    fn whole_record_copy_keeps_the_minimised_string() {
        let samples = genotyped_samples(&["S1", "S2"]);
        let split = IndividualSplit::new(&IndividualRequest::All, &samples, IndividualMode::Copies)
            .unwrap();
        let reading = RecordReading {
            individuals: Some(&split),
            process_ref_homs: true,
        };
        let line = "21\t100\tx\tCTA\tGA\t.\t.\t.\tGT\t0/1\t0/0";
        let plain = parse_vcf_record(line, RecordReading::default()).unwrap();
        assert_eq!(plain[0].record_allele_string(), "CT/G");
        let copies = parse_vcf_record(line, reading).unwrap();
        assert_eq!(copies.len(), 3);
        assert_eq!(copies[1].record_allele_string(), "CT/G");
        assert_eq!(copies[1].allele_string, "CT/G");
        assert_eq!(copies[2].record_allele_string(), "CT/CT");
        assert!(copies[2].hom_ref && copies[2].individual.as_deref() == Some("S2"));
        for row in [&copies[0], &copies[2]] {
            assert!(row.reference_allele_row);
            assert_eq!(row.alt_allele(), b"CT");
            assert_eq!(row.ref_allele, b"CT");
            assert_eq!((row.start, row.end), (100, 101));
        }
    }

    /// `(sample, ZYG)` of each entry of the record's zygosity list.
    fn zygosities(variant: &InputVariant) -> Vec<(String, String)> {
        variant
            .sample_zygosities
            .as_deref()
            .expect("a record read under --individual_zyg carries the list")
            .entries
            .iter()
            .map(|z| (z.sample.to_string(), z.zygosity.to_string()))
            .collect()
    }

    fn zyg(entries: &[(&str, &str)]) -> Vec<(String, String)> {
        entries
            .iter()
            .map(|(s, z)| (s.to_string(), z.to_string()))
            .collect()
    }

    /// `--individual_zyg all` (`Parser/VCF.pm:683-743`) on the record of
    /// [`individual_all_emits_one_variant_set_per_individual`]: the record stays
    /// its two alleles, no copy and no `IND`, and every allele carries one `ZYG`
    /// entry per sample whose genotype is called, read as `OutputFactory.pm:981-987`
    /// prints `genotype_ind` and `hom_ref_samples` (`:723-725`, `:737`): `0/1` is
    /// `HET` (A and G), `1/1` `HOM`, `1/2` `HET`, `./1` the one bit G so `HOM`,
    /// `./.` no entry (`:715`), a haploid `1` `HOM`, `0/3` (the `*` beside the
    /// reference) `HETREF` and `0/0` `HOMREF`, in header order. A named list is read
    /// in its order with a repeated name once (`:737` keys the entries by name); a
    /// record every sample misses keeps its alleles with an empty list, where
    /// `--individual` writes nothing; `--process_ref_homs` adds the reference row,
    /// which carries the list; a header without `FORMAT` reads no genotype; a
    /// structural record carries no list (`:244`).
    #[test]
    fn individual_zyg_keeps_one_variant_with_a_zyg_entry_per_sample() {
        let samples = genotyped_samples(&["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]);
        let split = IndividualSplit::new(
            &IndividualRequest::All,
            &samples,
            IndividualMode::Zygosities,
        )
        .unwrap();
        let reading = RecordReading {
            individuals: Some(&split),
            process_ref_homs: false,
        };
        let line = "21\t100\tx\tA\tG,T,*\t.\t.\t.\tGT\t0/1\t1/1\t1/2\t./1\t./.\t1\t0/3\t0/0";
        let variants = parse_vcf_record(line, reading).unwrap();
        let plain = parse_vcf_record(line, RecordReading::default()).unwrap();
        assert_eq!(variants.len(), 2);
        let alleles: Vec<&[u8]> = variants.iter().map(|v| v.alt_allele()).collect();
        assert_eq!(alleles, [b"G".as_slice(), b"T"]);
        assert!(variants
            .iter()
            .all(|v| v.individual.is_none() && v.genotype.is_empty() && !v.hom_ref));
        assert!(variants[0].same_record(&variants[1]));
        let expected = zyg(&[
            ("S1", "HET"),
            ("S2", "HOM"),
            ("S3", "HET"),
            ("S4", "HOM"),
            ("S6", "HOM"),
            ("S7", "HETREF"),
            ("S8", "HOMREF"),
        ]);
        assert_eq!(zygosities(&variants[0]), expected);
        assert_eq!(zygosities(&variants[1]), expected);
        assert_eq!(
            variants[0].record_allele_string(),
            plain[0].record_allele_string()
        );
        assert_eq!(variants[0].allele_num, plain[0].allele_num);
        assert!(plain.iter().all(|v| v.sample_zygosities.is_none()));

        let named = IndividualRequest::Named(vec!["S3".into(), "S1".into(), "S3".into()]);
        let named_split =
            IndividualSplit::new(&named, &samples, IndividualMode::Zygosities).unwrap();
        let listed = parse_vcf_record(
            line,
            RecordReading {
                individuals: Some(&named_split),
                process_ref_homs: true,
            },
        )
        .unwrap();
        assert_eq!(listed.len(), 3);
        assert!(listed[0].reference_allele_row && !listed[1].reference_allele_row);
        for variant in &listed {
            assert_eq!(zygosities(variant), zyg(&[("S3", "HET"), ("S1", "HET")]));
        }

        let missing = parse_vcf_record(
            "21\t100\tx\tA\tG\t.\t.\t.\tGT\t./.\t./.\t.\t./.\t./.\t./.\t./.\t./.",
            reading,
        )
        .unwrap();
        assert_eq!(missing.len(), 1);
        assert_eq!(zygosities(&missing[0]), []);
        let no_format = IndividualSplit::new(
            &IndividualRequest::All,
            &SampleColumns {
                names: samples.names.clone(),
                genotyped: false,
            },
            IndividualMode::Zygosities,
        )
        .unwrap();
        let ungenotyped = parse_vcf_record(
            line,
            RecordReading {
                individuals: Some(&no_format),
                process_ref_homs: false,
            },
        )
        .unwrap();
        assert_eq!(ungenotyped.len(), 2);
        assert_eq!(zygosities(&ungenotyped[0]), []);
        let structural = parse_vcf_record(
            "21\t100\tsv\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=200\tGT\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0",
            reading,
        )
        .unwrap();
        assert_eq!(structural.len(), 1);
        assert!(structural[0].sample_zygosities.is_none());
    }

    /// `Parser/VCF.pm:727-730`: when `--individual_zyg` names one sample
    /// (`$n_individuals`, the list's length, is 1) and that sample carries no
    /// alternate allele (`:720`, the `*` excluded as under `--individual`), the
    /// record's allele string becomes `REF/REF` and `unique_ind` routes the output
    /// through the reference allele's row alone (`OutputFactory.pm:506-508`), with
    /// the sample's `ZYG`; a list of one name repeated twice is two individuals, so
    /// the record keeps its alleles; a carrier or a missing genotype leaves the
    /// record as read; a non-variant record (`A .`) with a `0/0` sample is the
    /// reference row too. Under `--minimal` the row keeps the record's `minimised`
    /// mark and its uploaded string, which `minimise_alleles` left as Ensembl VEP's
    /// `original_allele_string` (Parser.pm:987-1011). The reference is the record's
    /// working one, `CT` for `CTA GA`, as a reference-homozygous copy's is under
    /// `--individual`: Ensembl VEP compares the chopped bits `CTA` and `GA` with the
    /// minimised `CT` (`:691-692`, `:720`) and reads that sample as `HOM` with the
    /// record's alleles kept, the divergence `docs/intended-divergences.md` lists
    /// under the trimmable indel, pinned here on the vep-rs side.
    #[test]
    fn individual_zyg_single_hom_ref_sample_annotates_the_reference_allele() {
        let samples = genotyped_samples(&["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]);
        let one = |name: &str| {
            IndividualSplit::new(
                &IndividualRequest::Named(vec![name.into()]),
                &samples,
                IndividualMode::Zygosities,
            )
            .unwrap()
        };
        let line = "21\t100\tx\tA\tG,T,*\t.\t.\t.\tGT\t0/1\t1/1\t1/2\t./1\t./.\t1\t0/3\t0/0";
        let read = |split: &IndividualSplit| {
            parse_vcf_record(
                line,
                RecordReading {
                    individuals: Some(split),
                    process_ref_homs: false,
                },
            )
            .unwrap()
        };
        let hom_ref = read(&one("S8"));
        assert_eq!(hom_ref.len(), 1);
        let row = &hom_ref[0];
        assert!(row.reference_allele_row && row.hom_ref && !row.no_allele);
        assert!(row.individual.is_none());
        assert_eq!(row.alt_allele(), b"A");
        assert_eq!(row.allele_num, Some(0));
        assert_eq!(row.record_allele_string(), "A/A");
        assert_eq!(row.uploaded_allele_string.as_deref(), Some("A/A"));
        assert_eq!(zygosities(row), zyg(&[("S8", "HOMREF")]));
        let star = read(&one("S7"));
        assert_eq!(star.len(), 1);
        assert!(star[0].reference_allele_row);
        assert_eq!(zygosities(&star[0]), zyg(&[("S7", "HETREF")]));
        let carrier = read(&one("S1"));
        assert_eq!(carrier.len(), 2);
        assert!(carrier.iter().all(|v| !v.reference_allele_row));
        assert_eq!(zygosities(&carrier[0]), zyg(&[("S1", "HET")]));
        let missing = read(&one("S5"));
        assert_eq!(missing.len(), 2);
        assert_eq!(zygosities(&missing[0]), []);
        let twice = IndividualSplit::new(
            &IndividualRequest::Named(vec!["S8".into(), "S8".into()]),
            &samples,
            IndividualMode::Zygosities,
        )
        .unwrap();
        let repeated = read(&twice);
        assert_eq!(repeated.len(), 2);
        assert!(repeated.iter().all(|v| !v.reference_allele_row));
        assert_eq!(zygosities(&repeated[0]), zyg(&[("S8", "HOMREF")]));
        // `all` over a single-sample header names one individual too.
        let single = IndividualSplit::new(
            &IndividualRequest::All,
            &genotyped_samples(&["S1"]),
            IndividualMode::Zygosities,
        )
        .unwrap();
        let reading = RecordReading {
            individuals: Some(&single),
            process_ref_homs: false,
        };
        let all_of_one = parse_vcf_record("21\t100\tx\tA\tG\t.\t.\t.\tGT\t0|0", reading).unwrap();
        assert_eq!(all_of_one.len(), 1);
        assert!(all_of_one[0].reference_allele_row);
        assert_eq!(zygosities(&all_of_one[0]), zyg(&[("S1", "HOMREF")]));
        let non_variant =
            non_variant_record_read("21\t100\tx\tA\t.\t.\t.\t.\tGT\t0/0", reading).unwrap();
        assert_eq!(non_variant.len(), 1);
        assert!(non_variant[0].reference_allele_row && !non_variant[0].no_allele);
        assert_eq!(non_variant[0].record_allele_string(), "A/A");
        assert_eq!(zygosities(&non_variant[0]), zyg(&[("S1", "HOMREF")]));
        let minimised = parse_vcf_record("21\t100\tx\tCTA\tGA\t.\t.\t.\tGT\t0/0", reading).unwrap();
        assert_eq!(minimised.len(), 1);
        assert_eq!(minimised[0].record_allele_string(), "CT/CT");
        assert_eq!(minimised[0].alt_allele(), b"CT");
        assert_eq!((minimised[0].start, minimised[0].end), (100, 101));
        // The row keeps the record's mark (`minimise_alleles` ran on the indel),
        // which only `--minimal` prints.
        assert!(minimised[0].minimised && !all_of_one[0].minimised);
        assert_eq!(all_of_one[0].uploaded_allele(false).as_deref(), Some("A/A"));
        let under_minimal = parse_vcf_line_with(
            "21\t100\tx\tG\tA\t.\t.\t.\tGT\t0|0",
            &VcfParse {
                minimal: true,
                reading,
                ..VcfParse::default_run(false)
            },
        )
        .unwrap();
        assert_eq!(under_minimal.len(), 1);
        assert!(under_minimal[0].reference_allele_row && under_minimal[0].minimised);
        assert_eq!(under_minimal[0].record_allele_string(), "G/G");
        assert_eq!(
            under_minimal[0].uploaded_allele(true).as_deref(),
            Some("G/A")
        );
        // A non-variant record has no alternate allele for `minimise_alleles` to
        // keep a string of (Parser.pm:974-976), so its row prints the reference twice.
        let non_variant_minimal = non_variant_record(
            "21\t100\tx\tC\t.\t.\t.\t.\tGT\t0/0",
            &VcfParse {
                minimal: true,
                reading,
                ..VcfParse::default_run(false)
            },
        )
        .unwrap();
        assert_eq!(
            non_variant_minimal[0].uploaded_allele(true).as_deref(),
            Some("C/C")
        );
        assert!(!non_variant_minimal[0].minimised);
    }

    /// `split_variants` reads only a feature whose allele string names two
    /// alternate alleles (InputBuffer.pm:648, `.+/.+/.+`), and `minimise_alleles`
    /// runs on the record before its per-individual copies are made
    /// (Parser/VCF.pm:362-364), so under `--minimal` a copy of a multi-allelic
    /// record carrying one alternate allele, or two bits naming one allele, is
    /// neither split nor minimised and keeps the record's chopped alleles; a copy
    /// naming two alleles that trim is split, in the record's ALT order, and a copy
    /// carrying the empty allele of a bit past the record's alleles (`3/1`) is not.
    #[test]
    fn a_bi_allelic_copy_of_a_multi_allelic_record_is_not_split_under_minimal() {
        let samples = genotyped_samples(&["S1", "S2", "S3", "S4", "S5"]);
        let split = IndividualSplit::new(&IndividualRequest::All, &samples, IndividualMode::Copies)
            .unwrap();
        let parse = VcfParse {
            minimal: true,
            reading: RecordReading {
                individuals: Some(&split),
                process_ref_homs: false,
            },
            ..VcfParse::default_run(false)
        };
        let line = "21\t100\tx\tAGCAA\tCGCAA,A\t.\t.\t.\tGT\t1/2\t0/1\t2/2\t3/1\t0/0";
        let copies = parse_vcf_line_with(line, &parse).unwrap();
        let of = |name: &str| -> Vec<&InputVariant> {
            copies
                .iter()
                .filter(|v| v.individual.as_deref() == Some(name))
                .collect()
        };
        assert_eq!(of("S1").len(), 2);
        assert!(of("S1").iter().all(|v| record_splits(v)));
        for (name, allele, working) in [
            ("S2", "CGCAA", "AGCAA/CGCAA"),
            ("S3", "A", "AGCAA/A"),
            ("S4", "CGCAA", "AGCAA/CGCAA/"),
        ] {
            let copy = of(name);
            assert_eq!(copy.len(), 1, "{name}");
            assert!(!record_splits(copy[0]), "{name}");
            assert_eq!(copy[0].record_allele_string(), working, "{name}");
            let kept = split_minimal_record(vec![copy[0].clone()]);
            assert_eq!(kept.len(), 1, "{name}");
            assert_eq!(kept[0].display_allele(), allele, "{name}");
            assert_eq!(kept[0].ref_allele, b"AGCAA", "{name}");
            assert_eq!((kept[0].start, kept[0].end), (100, 104), "{name}");
            assert!(!kept[0].minimised, "{name}");
            assert_eq!(kept[0].allele_num, Some(1), "{name}");
        }
        assert!(of("S5")[0].hom_ref);
        let split = split_minimal_record(of("S1").into_iter().cloned().collect());
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].display_allele(), "C");
        assert_eq!((split[0].start, split[0].end), (100, 100));
        assert_eq!(split[1].display_allele(), "-");
        assert_eq!((split[1].start, split[1].end), (101, 104));
        assert!(split.iter().all(|v| v.minimised));
    }
}
