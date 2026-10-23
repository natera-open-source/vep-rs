// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF output format: the original VCF line with a CSQ INFO entry appended.
//!
//! One line per input record: every allele's consequences become comma-joined
//! CSQ entries on that record's line, as VEP writes them. The CSQ subfields are
//! VEP's `@VCF_COLS` followed by the active flag fields and plugin fields, or the
//! `--fields` list alone (`build_csq_field_list`); values inside an entry are
//! escaped as VEP escapes them (`,` and `|` to `&`, `;` to `%3B`, whitespace to
//! `_`), lists join with `&`, and an absent value is empty. A record of another
//! input format gets the VCF line Ensembl VEP builds for it ([`to_vcf_record`]).

use std::borrow::Cow;
use std::collections::HashSet;
use std::io::Write;
use std::sync::Arc;

use vep_core::codon::reverse_complement;
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

use crate::error::IoError;
use crate::input::ReferenceSequence;

use super::fields::{self, FieldOptions, PluginHeader};
use super::OutputFormatter;

/// The field-selection options, shared with every other format.
pub type VcfFieldOptions = FieldOptions;

/// The key of the engine identity meta line, `##vep-rs=<version>`: a
/// `##key=value` line as VCF 4.3 section 1.4 defines meta-information, where
/// the `## vep-rs <version>` form of the text formats would be malformed.
const VCF_IDENTITY_KEY: &str = "##vep-rs=";

/// VCF output formatter that adds a CSQ INFO field to original VCF lines.
pub struct VcfOutputFormatter {
    /// Name of the INFO field (default "CSQ").
    pub info_field_name: String,
    /// Sub-field names within the CSQ value.
    pub fields: Vec<String>,
    /// Per subfield, whether the configuration produces it; a `--fields` name
    /// whose flag is off is an empty subfield in every entry.
    produced: Vec<bool>,
    /// Runtime field-selection options mirroring CLI flags.
    options: FieldOptions,
    /// The plugin fields and their descriptions, one `##<KEY>=<description>`
    /// header line each.
    plugin_headers: Vec<PluginHeader>,
    /// `--keep_csq`: an input's existing annotation stays, header lines and INFO
    /// entry, beside the new one.
    keep_csq: bool,
    /// The `##VEP=` meta line (version, time, cache and source versions).
    vep_meta_line: String,
    /// The invocation for the `##VEP-command-line=` meta line.
    command_line: String,
    /// Whether the records come from a VCF, whose own line each output line
    /// extends (OutputFactory/VCF.pm:314-317); a record of another input format
    /// carries its input line too, and its VCF line is built from the record.
    vcf_input: bool,
    /// The run's reference FASTA, read for the base before an indel of a record
    /// built into a VCF line; without one that base is `N`.
    reference: Option<Arc<dyn ReferenceSequence>>,
}

impl VcfOutputFormatter {
    /// `selected_fields` is the `--fields` list, which replaces the CSQ subfields.
    pub fn new(
        info_field_name: String,
        options: FieldOptions,
        plugin_headers: Vec<PluginHeader>,
        selected_fields: Option<Vec<String>>,
    ) -> Self {
        let fields = build_csq_field_list(
            &options,
            &fields::plugin_field_names(&plugin_headers),
            selected_fields.as_deref(),
        );
        let produced = fields
            .iter()
            .map(|f| fields::field_is_produced(f, &options))
            .collect();
        Self {
            info_field_name,
            fields,
            produced,
            options,
            plugin_headers,
            keep_csq: false,
            vep_meta_line: fields::vcf_vep_meta_line("", None),
            command_line: String::new(),
            vcf_input: true,
            reference: None,
        }
    }

    /// Keeps an input VCF's existing annotation (`--keep_csq`): its `##INFO`,
    /// `##VEP` and `##vep-rs` header lines and its INFO entry under the field's
    /// name stay, and the new entry is appended under the same name.
    pub fn with_keep_csq(mut self, keep_csq: bool) -> Self {
        self.keep_csq = keep_csq;
        self
    }

    /// Sets whether the input is a VCF, whose records' own lines the output
    /// extends; `true` unless set.
    pub fn with_vcf_input(mut self, vcf_input: bool) -> Self {
        self.vcf_input = vcf_input;
        self
    }

    /// Sets the reference the base before an indel is read from when a record
    /// of a non-VCF input is built into a VCF line (`--fasta`); `N` without one.
    pub fn with_reference(mut self, reference: Option<Arc<dyn ReferenceSequence>>) -> Self {
        self.reference = reference;
        self
    }

    /// The record's VCF line, read from any allele `first` of it: a VCF input's
    /// own line as read, the bytes when the line is not UTF-8
    /// (`OutputFactory/VCF.pm:315`), else the line Ensembl VEP builds from the
    /// record (`:336-337`, [`to_vcf_record`]).
    fn record_line<'a>(&self, first: &'a InputVariant) -> Cow<'a, [u8]> {
        if self.vcf_input {
            if let Some(bytes) = &first.raw_input_bytes {
                return Cow::Borrowed(bytes);
            }
            if let Some(line) = &first.raw_input {
                return Cow::Borrowed(line.as_bytes());
            }
        }
        Cow::Owned(to_vcf_record(first, self.reference.as_deref()))
    }

    /// Sets the run description the meta lines carry: the cache directory,
    /// its info and the command line as VEP prints it.
    pub fn with_run_info(
        mut self,
        cache_dir: &str,
        cache_info: Option<&vep_core::cache_info::CacheInfo>,
        command_line: &str,
    ) -> Self {
        self.vep_meta_line = fields::vcf_vep_meta_line(cache_dir, cache_info);
        self.command_line = command_line.to_string();
        self
    }

    /// The header block VEP writes for a VCF output given the input's own
    /// header lines (`##` metas and the `#CHROM` line, in order; empty when the
    /// input was not a VCF): a `##fileformat` line when the input lacks one, the
    /// input metas minus any earlier definition of the INFO field, `##VEP` or
    /// `##vep-rs` line (every input meta under `--keep_csq`), the `##VEP=` line,
    /// the `##vep-rs=<version>` line, the field's `##INFO` definition, one
    /// `##<KEY>=<description>` line per plugin field, the command line, then the
    /// column line (`OutputFactory/VCF.pm` `headers`).
    pub fn header_lines(&self, input_headers: &[String]) -> Vec<String> {
        let mut lines: Vec<String> = Vec::new();
        let mut col_heading = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
        if input_headers.is_empty() {
            lines.push("##fileformat=VCFv4.1".to_string());
        } else {
            if !input_headers
                .iter()
                .any(|h| h.to_ascii_lowercase().starts_with("##fileformat=vcfv4"))
            {
                lines.push("##fileformat=VCFv4.1".to_string());
            }
            let csq_marker = format!("ID={},", self.info_field_name).to_ascii_lowercase();
            let mut metas: Vec<&String> = input_headers.iter().collect();
            if let Some(heading) = metas.pop_if(|h| h.starts_with("#CHROM")) {
                col_heading = heading.clone();
            }
            for h in metas {
                if !self.keep_csq
                    && (h.to_ascii_lowercase().contains(&csq_marker)
                        || h.starts_with("##VEP")
                        || h.starts_with(VCF_IDENTITY_KEY))
                {
                    continue;
                }
                lines.push(h.clone());
            }
        }
        lines.push(self.vep_meta_line.clone());
        lines.push(format!("{VCF_IDENTITY_KEY}{}", fields::VEP_RS_VERSION));
        lines.push(self.info_line());
        for (name, description) in &self.plugin_headers {
            lines.push(format!("##{name}={description}"));
        }
        lines.push(format!("##VEP-command-line='{}'", self.command_line));
        lines.push(col_heading);
        lines
    }

    fn info_line(&self) -> String {
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {}\">",
            self.info_field_name,
            self.fields.join("|"),
        )
    }
}

impl Default for VcfOutputFormatter {
    fn default() -> Self {
        Self::new("CSQ".to_string(), FieldOptions::default(), Vec::new(), None)
    }
}

/// VEP's fixed CSQ prefix (`@VCF_COLS`).
const VCF_COLS: [&str; 18] = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
];

/// The CSQ subfields: `@VCF_COLS`, then the active flag fields not already in
/// it, then plugin fields, each name once; or, under `--fields`, that list in its
/// own order and nothing else (`OutputFactory/VCF.pm` `fields`, where the flag's
/// value pre-empts the list).
pub(crate) fn build_csq_field_list(
    options: &FieldOptions,
    plugin_fields: &[String],
    selected_fields: Option<&[String]>,
) -> Vec<String> {
    if let Some(selected_fields) = selected_fields {
        return selected_fields.to_vec();
    }
    let mut fields: Vec<String> = VCF_COLS.iter().map(|s| s.to_string()).collect();
    for f in fields::flag_fields(options) {
        if !fields.iter().any(|x| x == f) {
            fields.push(f.to_string());
        }
    }
    for f in plugin_fields {
        if !fields.contains(f) {
            fields.push(f.clone());
        }
    }
    fields
}

/// Perl's `\s` on a byte string, as the output modules match it (no
/// `unicode_strings` pragma): the six ASCII whitespace characters, vertical tab
/// included. `char::is_ascii_whitespace` omits the vertical tab, and
/// `char::is_whitespace` adds the Unicode spaces, which Perl leaves as they are.
fn is_perl_whitespace(c: char) -> bool {
    c.is_ascii_whitespace() || c == '\x0B'
}

/// Escape one CSQ subfield value the way VEP does (`OutputFactory/VCF.pm:401-404`):
/// `,` and `|` become `&`, `;` becomes `%3B`, a run of ASCII whitespace becomes
/// one `_`. A control character that is not whitespace, which VEP writes as it
/// is, becomes its `%XX` percent-encoding, as the default and tab formats write
/// it ([`fields::percent_encode_control_bytes`]).
fn vcf_encode(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut in_whitespace = false;
    for ch in value.chars() {
        let whitespace = is_perl_whitespace(ch);
        if whitespace {
            if !in_whitespace {
                out.push('_');
            }
        } else {
            match ch {
                ';' => out.push_str("%3B"),
                ',' | '|' => out.push('&'),
                c if c.is_ascii_control() => {
                    out.push_str(&format!("%{:02X}", c as u8));
                }
                c => out.push(c),
            }
        }
        in_whitespace = whitespace;
    }
    out
}

pub(crate) fn value_for_field(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> String {
    fields::field_value(field, variant, tc, options, "&")
}

/// One CSQ entry. A gene symbol's whitespace is percent-encoded before the
/// subfield escaping turns whitespace into `_`, so it reaches the output as
/// `%20`; `--no_escape` leaves it to the `_` rule (`OutputFactory.pm`,
/// `BaseTranscriptVariationAllele_to_output_hash`, the `symbol` block). The
/// `Allele` of a record on the reverse strand is reverse-complemented, whatever
/// its text (`OutputFactory/VCF.pm:392-393`, `reverse_comp` on the value), so
/// it names the allele as the record's VCF line does.
fn format_csq_entry(
    fields: &[String],
    produced: &[bool],
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> String {
    fields
        .iter()
        .zip(produced)
        .map(|(field, produced)| {
            if !produced {
                return String::new();
            }
            let v = value_for_field(field, variant, tc, options);
            if v == "-" && field != "Allele" {
                String::new()
            } else if field == "Allele" && variant.strand.as_i8() < 0 {
                vcf_encode(&String::from_utf8_lossy(&reverse_complement(v.as_bytes())))
            } else if field == "SYMBOL" && !options.no_escape {
                vcf_encode(&v.replace(is_perl_whitespace, "%20"))
            } else {
                vcf_encode(&v)
            }
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// Whether `allele` is one or more of the bases `A`, `C`, `G`, `T`, `N`
/// (Perl's `/^[ACGTN]+$/`).
fn is_acgtn(allele: &[u8]) -> bool {
    !allele.is_empty() && allele.iter().all(|b| b"ACGTN".contains(b))
}

/// Perl's truth of an allele's text: the empty string and `0` are false, every
/// other text true.
fn perl_true_bytes(text: &[u8]) -> bool {
    !text.is_empty() && text != b"0"
}

/// ensembl-variation `Utils/Sequence.pm` `trim_right` (release 116, lines
/// 1058-1096): while every allele ends in the same base and none is one base
/// long, that base is removed from each; an allele with no base counts as
/// ending in no base, so nothing is trimmed beside one.
fn trim_right(alleles: &mut [Vec<u8>]) {
    if alleles.len() < 2 {
        return;
    }
    loop {
        if alleles.iter().any(|allele| allele.len() == 1) {
            return;
        }
        let Some(last) = alleles[0].last().copied() else {
            return;
        };
        if alleles.iter().any(|allele| allele.last() != Some(&last)) {
            return;
        }
        for allele in alleles.iter_mut() {
            allele.pop();
        }
    }
}

/// The base before `start` on the record's chromosome, upper-cased, read from
/// `reference`; `N` without a reference, for a chromosome it lacks or at the
/// chromosome's first base (`BaseVariationFeature.pm` `_get_prev_base`, release
/// 116, lines 339-353).
fn base_before(
    reference: Option<&dyn ReferenceSequence>,
    variant: &InputVariant,
    start: u64,
) -> Vec<u8> {
    if start <= 1 {
        return b"N".to_vec();
    }
    let position = start - 1;
    reference
        .and_then(|reference| {
            reference
                .sequence(&variant.chr, position, position)
                .or_else(|| reference.sequence(&variant.original_chr, position, position))
        })
        .unwrap_or_else(|| b"N".to_vec())
}

/// The eight fixed VCF columns `chr`, `pos`, `id`, `ref`, `alt`, `.`, `.`,
/// `info`, tab-separated.
fn fixed_columns(
    chr: &str,
    pos: u64,
    id: &str,
    reference: &[u8],
    alt: &[u8],
    info: &[u8],
) -> Vec<u8> {
    let mut line = Vec::new();
    for (i, column) in [
        chr.as_bytes(),
        pos.to_string().as_bytes(),
        id.as_bytes(),
        reference,
        alt,
        b".",
        b".",
        info,
    ]
    .into_iter()
    .enumerate()
    {
        if i > 0 {
            line.push(b'\t');
        }
        let start = line.len();
        line.extend_from_slice(column);
        fields::percent_encode_control_bytes(&mut line, start);
    }
    line
}

/// The VCF line Ensembl VEP builds for a record that did not come from a VCF,
/// as `OutputFactory/VCF.pm:336-337` (release 116.2) builds one for every such
/// record: ensembl-variation `VariationFeature.pm` `to_VCF_record` (release
/// 116, lines 2340-2444) for a sequence variant and
/// `StructuralVariationFeature.pm` `to_VCF_record` (lines 1149-1171) for a
/// structural one. `first` is any allele of the record: the line is built from
/// the record's own allele string, strand, start and name, which every allele
/// carries.
///
/// A sequence variant's alleles are the record's allele string split on `/`;
/// on the reverse strand each allele of `A`, `C`, `G`, `T` and `N` alone is
/// reverse-complemented (`:2376`), and every `-` is removed (`:2381`). Alleles
/// of one length that are all bases or empty make a balanced record: the
/// record's start, the first allele as REF and the rest comma-joined as ALT,
/// `.` when the join is empty or `0` (`:2432-2441`). Otherwise the alleles are
/// trimmed on the right (`trim_right`, `:2390`, unless one holds a character
/// outside `ACGTN.`; an allele Perl reads as false, empty or `0`, is never
/// such a one, `:2383`), a lone allele gets `-` before it (`:2392`), and the
/// alleles get an anchor base: none when every allele free of `*` opens with
/// the same base (`:2399-2402`), else the base before the record's start and
/// the record's position moved one base back (`:2404-2407`, `base_before`).
/// An allele of bases or empty is written after the anchor, any other inside
/// `<` and `>` (`:2409-2419`). The ID column is the record's name (`:2425`).
/// The `COSMIC_MUTATION` and `HGMD_MUTATION` strings (`:2356-2373`) name
/// database records and are not an input shape. Every column's control bytes
/// are percent-encoded as the default and tab formats write them
/// ([`fields::percent_encode_control_bytes`]), so a tab or line end inside a
/// name or a chromosome cannot split the line.
///
/// A structural variant's line is the base before its start and its position
/// moved one base back, `<INS>`, `<DEL>`, `<TDUP>` or `<DUP>` for an
/// insertion, deletion, tandem duplication or duplication and the record's
/// class term inside `<` and `>` for any other class, with `END=` and the
/// record's end as INFO.
pub fn to_vcf_record(first: &InputVariant, reference: Option<&dyn ReferenceSequence>) -> Vec<u8> {
    let chr = &first.original_chr;
    let name = first.variation_name();
    let start = first.original_start.unwrap_or(first.start);

    if first.is_structural {
        let class: Cow<str> = match &first.record_so_term {
            Some(term) => Cow::Borrowed(term),
            None => String::from_utf8_lossy(first.alt_allele()),
        };
        let term = match class.as_ref() {
            "insertion" => "INS",
            "deletion" => "DEL",
            "tandem_duplication" => "TDUP",
            "duplication" => "DUP",
            other => other,
        };
        let end = first.original_end.unwrap_or(first.end);
        return fixed_columns(
            chr,
            start.saturating_sub(1),
            &name,
            &base_before(reference, first, start),
            format!("<{term}>").as_bytes(),
            format!("END={end}").as_bytes(),
        );
    }

    let allele_string = first.working_allele_string();
    let mut alleles: Vec<Vec<u8>> = allele_string
        .split('/')
        .map(|allele| {
            let allele = allele.as_bytes();
            if first.strand.as_i8() < 0 && is_acgtn(allele) {
                reverse_complement(allele)
            } else {
                allele.to_vec()
            }
        })
        .collect();
    for allele in &mut alleles {
        allele.retain(|b| *b != b'-');
    }
    let lengths: HashSet<usize> = alleles.iter().map(Vec::len).collect();
    let non_acgt = alleles
        .iter()
        .any(|allele| perl_true_bytes(allele) && !allele.iter().all(|b| b"ACGTN.".contains(b)));

    let position = if non_acgt || lengths.len() > 1 {
        if !non_acgt {
            trim_right(&mut alleles);
        }
        if alleles.len() == 1 {
            alleles.insert(0, b"-".to_vec());
        }
        let first_bases: HashSet<Option<u8>> = alleles
            .iter()
            .filter(|allele| !allele.contains(&b'*'))
            .map(|allele| allele.first().copied())
            .collect();
        let (anchor, position) = if first_bases.len() == 1 {
            (Vec::new(), start)
        } else {
            (
                base_before(reference, first, start),
                start.saturating_sub(1),
            )
        };
        for allele in &mut alleles {
            allele.retain(|b| *b != b'-');
            if allele.is_empty() || is_acgtn(allele) {
                let mut anchored = anchor.clone();
                anchored.extend_from_slice(allele);
                *allele = anchored;
            } else {
                allele.insert(0, b'<');
                allele.push(b'>');
            }
        }
        position
    } else {
        start
    };

    let reference_allele = alleles.remove(0);
    // The alternates joined with `,`, then Perl's truth test for the `.`
    // fallback (lines 2427 and 2440).
    let mut alt = alleles.join(&b","[..]);
    if !perl_true_bytes(&alt) {
        alt = b".".to_vec();
    }
    fixed_columns(chr, position, &name, &reference_allele, &alt, b".")
}

impl OutputFormatter for VcfOutputFormatter {
    fn write_header(&mut self, writer: &mut dyn Write) -> Result<(), IoError> {
        writeln!(writer, "{}", self.info_line())?;
        Ok(())
    }

    fn format_variant(&self, variant: &InputVariant) -> Result<Vec<String>, IoError> {
        self.format_record(&[variant])
    }

    /// The line's bytes read as text, each invalid sequence replaced; the
    /// pipeline writes the bytes of [`Self::format_record_bytes`].
    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        Ok(self
            .format_record_bytes(record)?
            .into_iter()
            .map(|line| String::from_utf8_lossy(&line).into_owned())
            .collect())
    }

    fn format_record_bytes(&self, record: &[&InputVariant]) -> Result<Vec<Vec<u8>>, IoError> {
        if record.is_empty() {
            return Ok(Vec::new());
        }
        let line = self.record_line(record[0]);
        let mut cols: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
        while cols.len() < 8 {
            cols.push(&b"."[..]);
        }
        // Exactly the first entry under the field's name, whatever its value, is
        // removed, whether or not new consequences follow; --keep_csq keeps it
        // beside the new one. OutputFactory/VCF.pm:327-331 removes it by one
        // substitution, `s/(^|\;)$fieldname\=\S+?(\;\S|$)/$2/`, whose `\S+?` must
        // take one character, so on an empty-valued entry (`ANN=;`) it also removes
        // the field that follows; that accident is not followed here.
        let csq_key = format!("{}=", self.info_field_name);
        let mut info: Vec<u8> = if cols[7] == b"." {
            Vec::new()
        } else if self.keep_csq {
            cols[7].to_vec()
        } else {
            let mut entries: Vec<&[u8]> = cols[7].split(|b| *b == b';').collect();
            if let Some(first) = entries
                .iter()
                .position(|kv| kv.starts_with(csq_key.as_bytes()))
            {
                entries.remove(first);
                // A `;` closing the removed entry at the end of the column goes
                // with it, as the substitution's `\S+?` takes it before `$`
                // matches, so the strip leaves no empty field behind.
                if entries.len() == first + 1 && entries[first].is_empty() {
                    entries.pop();
                }
            }
            entries.join(&b";"[..])
        };
        let carried_without_csq = record
            .iter()
            .any(|v| v.vep_skip || v.oversize_sv || v.no_allele);
        let csq_parts: Vec<String> = if carried_without_csq {
            Vec::new()
        } else {
            fields::record_rows(record)
                .into_iter()
                .map(|(variant, tc)| {
                    format_csq_entry(&self.fields, &self.produced, variant, tc, &self.options)
                })
                .collect()
        };
        // A record with no consequences is carried through without the field, as
        // Ensembl VEP writes a record it annotates nothing for (OutputFactory/VCF.pm:341-353
        // appends the field only when there is a chunk and leaves `.` otherwise).
        if !csq_parts.is_empty() {
            if !info.is_empty() {
                info.push(b';');
            }
            info.extend_from_slice(csq_key.as_bytes());
            info.extend_from_slice(csq_parts.join(",").as_bytes());
        }
        if info.is_empty() {
            info.push(b'.');
        }
        cols[7] = &info;
        Ok(vec![cols.join(&b"\t"[..])])
    }

    fn finish(&mut self, _writer: &mut dyn Write) -> Result<(), IoError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, FeatureType, Impact, TranscriptConsequence};

    #[test]
    fn test_csq_header_generation() {
        let mut formatter = VcfOutputFormatter::default();
        let mut buf = Vec::new();
        formatter.write_header(&mut buf).unwrap();
        let header = String::from_utf8(buf).unwrap();

        assert!(header.starts_with("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: "));
        assert!(header.contains("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation"));
    }

    fn plugin_headers() -> Vec<PluginHeader> {
        vec![
            (
                "CADD_PHRED".to_string(),
                "CADD PHRED-scaled pathogenicity score".to_string(),
            ),
            (
                "REVEL_score".to_string(),
                "REVEL pathogenicity score".to_string(),
            ),
        ]
    }

    /// The identity meta follows the `##VEP=` line in the key-value form of
    /// VCF 4.3 section 1.4, and a re-annotation replaces an earlier one along
    /// with the earlier `##VEP` lines rather than keeping both.
    #[test]
    fn identity_meta_follows_the_vep_line_and_replaces_an_earlier_one() {
        let formatter = VcfOutputFormatter::default();
        let identity = format!("##vep-rs={}", fields::VEP_RS_VERSION);
        let lines = formatter.header_lines(&[]);
        assert_eq!(lines[0], "##fileformat=VCFv4.1");
        assert!(lines[1].starts_with("##VEP=\"v"));
        assert_eq!(lines[2], identity);
        assert!(lines[3].starts_with("##INFO=<ID=CSQ,"));
        assert_eq!(
            lines.iter().filter(|l| l.starts_with("##vep-rs=")).count(),
            1
        );

        let input = [
            "##fileformat=VCFv4.2".to_string(),
            "##VEP=\"v116.2\" API=\"v116\"".to_string(),
            "##vep-rs=0.4.0".to_string(),
            "##source=keep".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let lines = formatter.header_lines(&input);
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(lines[1], "##source=keep");
        assert!(lines[2].starts_with("##VEP=\"v"));
        assert_eq!(lines[3], identity);
        assert!(!lines.contains(&"##vep-rs=0.4.0".to_string()));
    }

    #[test]
    fn test_dynamic_plugin_fields_in_header() {
        let mut formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                allele_number: true,
                variant_class: true,
                ..Default::default()
            },
            plugin_headers(),
            None,
        );
        let mut buf = Vec::new();
        formatter.write_header(&mut buf).unwrap();
        let header = String::from_utf8(buf).unwrap();
        assert!(header.contains("ALLELE_NUM"));
        assert!(header.contains("VARIANT_CLASS"));
        assert!(header.contains("CADD_PHRED"));
        assert!(header.contains("REVEL_score"));
    }

    /// A one-transcript variant whose raw line carries an earlier annotation.
    fn annotated_variant() -> InputVariant {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.raw_input = Some("21\t100\trs1\tA\tG\t.\t.\tDP=3;CSQ=G|old|entry".into());
        v.id = Some("rs1".into());
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            gene_symbol: Some("GENE1".into()),
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            feature_type: FeatureType::Transcript,
            strand: 1,
            distance: None,
            ..Default::default()
        }];
        v
    }

    fn info_of(line: &str) -> &str {
        line.split('\t').nth(7).unwrap()
    }

    /// The options a VCF output runs with: the `vcf` option set switches on
    /// symbol, biotype and numbers, so the fixed CSQ columns are produced.
    fn vcf_options() -> VcfFieldOptions {
        VcfFieldOptions {
            symbol: true,
            biotype: true,
            numbers: true,
            ..Default::default()
        }
    }

    /// `OutputFactory/VCF.pm` `fields` (release 116.2): the `--fields` value stands
    /// in for the whole subfield list, `@VCF_COLS` included, in the order given, so
    /// the `##INFO` Format string and every entry follow it; a name nothing
    /// produces is an empty subfield (`output_hash_to_vcf_info_chunk` pushes `''`
    /// for a key the hash lacks).
    #[test]
    fn fields_replace_the_csq_format() {
        let chosen: Vec<String> = ["Feature", "Allele", "SYMBOL", "STRAND", "NOT_A_FIELD"]
            .map(str::to_string)
            .to_vec();
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                symbol: true,
                ..Default::default()
            },
            plugin_headers(),
            Some(chosen.clone()),
        );
        assert_eq!(formatter.fields, chosen);
        let header = formatter.header_lines(&[]);
        assert!(
            header.contains(
                &"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Feature|Allele|SYMBOL|STRAND|NOT_A_FIELD\">".to_string()
            ),
            "{header:?}"
        );
        let lines = formatter.format_record(&[&annotated_variant()]).unwrap();
        assert_eq!(info_of(&lines[0]), "DP=3;CSQ=ENST1|G|GENE1|1|");
        // A named field whose flag is off is an empty subfield: VEP's row hash
        // carries no SYMBOL key without --symbol.
        let off = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions::default(),
            Vec::new(),
            Some(["Feature", "SYMBOL", "Allele"].map(str::to_string).to_vec()),
        );
        let lines = off.format_record(&[&annotated_variant()]).unwrap();
        assert_eq!(info_of(&lines[0]), "DP=3;CSQ=ENST1||G");
        let default = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions::default(),
            plugin_headers(),
            None,
        );
        assert_eq!(&default.fields[..18], &VCF_COLS);
        assert_eq!(
            &default.fields[18..],
            &["DISTANCE", "STRAND", "FLAGS", "CADD_PHRED", "REVEL_score"]
        );
    }

    /// `OutputFactory/VCF.pm` `headers` 215-223 and `get_all_lines_by_InputBuffer`
    /// 327-331 (release 116.2): under `--keep_csq` every input header line is
    /// kept, the earlier `##INFO=<ID=CSQ` and `##VEP` lines included, and the
    /// line's existing `CSQ=` entry stays while the new one is appended under the
    /// same key, so the INFO column carries two `CSQ=` keys; without the flag the
    /// earlier definition and entry are removed.
    #[test]
    fn keep_csq_keeps_the_old_info_line_and_the_old_entry() {
        let input_headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##VEP=\"v116.2\" API=\"v116\" time=\"2026-01-01 00:00:00\"".to_string(),
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence\">".to_string(),
            "##VEP-command-line='vep --vcf'".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let kept = VcfOutputFormatter::new("CSQ".to_string(), vcf_options(), Vec::new(), None)
            .with_keep_csq(true);
        let header = kept.header_lines(&input_headers);
        assert_eq!(&header[..5], &input_headers[..5]);
        assert_eq!(
            header
                .iter()
                .filter(|l| l.starts_with("##INFO=<ID=CSQ,"))
                .count(),
            2
        );
        assert_eq!(header.iter().filter(|l| l.starts_with("##VEP=")).count(), 2);
        assert_eq!(header.last().unwrap(), &input_headers[5]);
        let info = kept.format_record(&[&annotated_variant()]).unwrap();
        let info = info_of(&info[0]).to_string();
        assert!(
            info.starts_with("DP=3;CSQ=G|old|entry;CSQ=G|missense_variant|MODERATE|GENE1|ENSG1|"),
            "{info}"
        );
        assert_eq!(info.matches("CSQ=").count(), 2);

        let replaced = VcfOutputFormatter::new("CSQ".to_string(), vcf_options(), Vec::new(), None);
        let header = replaced.header_lines(&input_headers);
        assert_eq!(
            &header[..2],
            &input_headers[..1]
                .iter()
                .chain([&input_headers[4]])
                .cloned()
                .collect::<Vec<_>>()[..]
        );
        assert_eq!(
            header
                .iter()
                .filter(|l| l.starts_with("##INFO=<ID=CSQ,"))
                .count(),
            1
        );
        let info = replaced.format_record(&[&annotated_variant()]).unwrap();
        let info = info_of(&info[0]).to_string();
        assert!(info.starts_with("DP=3;CSQ=G|missense_variant|"), "{info}");
        assert_eq!(info.matches("CSQ=").count(), 1);
    }

    /// `OutputFactory/VCF.pm` `headers` 266-273 (release 116.2): after the `##INFO`
    /// definition, whose Format string is the `--fields` list when one is given,
    /// one `##<KEY>=<description>` line per plugin field, then the command line
    /// and the column line.
    #[test]
    fn header_carries_plugin_key_lines_and_the_fields_format() {
        let chosen: Vec<String> = ["Allele", "Consequence", "CADD_PHRED"]
            .map(str::to_string)
            .to_vec();
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions::default(),
            plugin_headers(),
            Some(chosen),
        )
        .with_run_info("/c", None, "vep --vcf");
        let header = formatter.header_lines(&[]);
        let tail: Vec<&str> = header.iter().map(String::as_str).skip(3).collect();
        assert_eq!(
            tail,
            [
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|CADD_PHRED\">",
                "##CADD_PHRED=CADD PHRED-scaled pathogenicity score",
                "##REVEL_score=REVEL pathogenicity score",
                "##VEP-command-line='vep --vcf'",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            ]
        );
        assert_eq!(header[0], "##fileformat=VCFv4.1");
        assert!(header[1].starts_with("##VEP=\"v"));
        assert!(header[2].starts_with(VCF_IDENTITY_KEY));
    }

    /// `--vcf_info_field ANN` (`VCF.pm` 147, 221, 267, 328, 349): the `##INFO`
    /// definition and the INFO key take the name, an input header defining `ANN`
    /// is removed while one defining `CSQ` stays, and an existing `CSQ=` entry on
    /// the line stays beside the new `ANN=` entry, since only the named field's
    /// entry is replaced.
    #[test]
    fn vcf_info_field_ann_names_header_strip_and_key() {
        let formatter = VcfOutputFormatter::new(
            "ANN".to_string(),
            VcfFieldOptions::default(),
            Vec::new(),
            None,
        );
        let input_headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=ANN,Number=.,Type=String,Description=\"old\">".to_string(),
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"other\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let header = formatter.header_lines(&input_headers);
        assert!(!header.iter().any(|l| l.contains("Description=\"old\"")));
        assert!(header.contains(&input_headers[2]));
        assert!(header
            .iter()
            .any(|l| l.starts_with("##INFO=<ID=ANN,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|")));
        let lines = formatter.format_record(&[&annotated_variant()]).unwrap();
        let info = info_of(&lines[0]);
        assert!(
            info.starts_with("DP=3;CSQ=G|old|entry;ANN=G|missense_variant|"),
            "{info}"
        );
    }

    /// An INFO ending in `;` after the entry under the field's name: the `\S+?` of
    /// `OutputFactory/VCF.pm:329` takes that `;` with the value before `$` matches,
    /// so Ensembl VEP 116.2 writes `DP=87;ANN=T|...` for the input
    /// `DP=87;ANN=T|old;` and no empty field is left between `DP=87` and the new
    /// entry. A `;` ending the column elsewhere stays, as it does there
    /// (`ANN=x;DP=87;` gives `DP=87;;ANN=T|...` in both).
    #[test]
    fn a_semicolon_closing_the_removed_entry_goes_with_it() {
        let formatter = VcfOutputFormatter::new(
            "ANN".to_string(),
            VcfFieldOptions::default(),
            Vec::new(),
            None,
        );
        let mut trailing = annotated_variant();
        trailing.raw_input = Some("21\t100\trs1\tA\tG\t.\t.\tDP=87;ANN=T|old;".into());
        let lines = formatter.format_record(&[&trailing]).unwrap();
        let info = info_of(&lines[0]);
        assert!(info.starts_with("DP=87;ANN=G|missense_variant|"), "{info}");

        let mut elsewhere = annotated_variant();
        elsewhere.raw_input = Some("21\t100\trs1\tA\tG\t.\t.\tANN=x;DP=87;".into());
        let lines = formatter.format_record(&[&elsewhere]).unwrap();
        let info = info_of(&lines[0]);
        assert!(info.starts_with("DP=87;;ANN=G|missense_variant|"), "{info}");
    }

    #[test]
    fn test_format_variant_with_one_consequence() {
        let formatter = VcfOutputFormatter::new("CSQ".to_string(), vcf_options(), Vec::new(), None);
        let mut variant = InputVariant::new(
            "21".into(),
            25585733,
            25585733,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25585733\trs123\tA\tG\t100\tPASS\t.".into());
        variant.id = Some("rs123".into());
        variant.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST00000352957".into(),
            gene_id: "ENSG00000154719".into(),
            gene_symbol: Some("MRPL39".into()),
            consequences: smallvec![Consequence::SynonymousVariant],
            impact: Impact::LOW,
            feature_type: FeatureType::Transcript,
            biotype: Some("protein_coding".into()),
            cdna_position: Some("1033".into()),
            cds_position: Some("991".into()),
            protein_position: Some("331".into()),
            amino_acids: Some("A/A".into()),
            codons: Some("Gca/Gca".into()),
            strand: -1,
            ..Default::default()
        }];

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(lines.len(), 1);

        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields[0], "21");
        assert_eq!(fields[1], "25585733");
        let info = fields[7];
        assert!(
            info.starts_with("CSQ="),
            "INFO field should start with CSQ=, got: {}",
            info
        );
        assert!(info.contains("synonymous_variant"));
        assert!(info.contains("MRPL39"));
        assert!(info.contains("ENST00000352957"));
        assert!(info.contains("protein_coding"));
    }

    #[test]
    fn test_format_variant_with_multiple_consequences() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25585733,
            25585733,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25585733\t.\tA\tG\t100\tPASS\t.".into());
        variant.transcript_consequences = vec![
            TranscriptConsequence {
                transcript_id: "ENST00000001".into(),
                gene_id: "ENSG00000001".into(),
                consequences: smallvec![Consequence::MissenseVariant],
                impact: Impact::MODERATE,
                feature_type: FeatureType::Transcript,
                strand: 1,
                ..Default::default()
            },
            TranscriptConsequence {
                transcript_id: "ENST00000002".into(),
                gene_id: "ENSG00000001".into(),
                consequences: smallvec![Consequence::IntronVariant],
                impact: Impact::MODIFIER,
                feature_type: FeatureType::Transcript,
                strand: 1,
                ..Default::default()
            },
        ];

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(
            lines.len(),
            1,
            "VCF should produce exactly one line per variant"
        );

        let fields: Vec<&str> = lines[0].split('\t').collect();
        let info = fields[7];
        let csq_value = info.strip_prefix("CSQ=").unwrap();
        let csq_entries: Vec<&str> = csq_value.split(',').collect();
        assert_eq!(csq_entries.len(), 2, "Should have 2 CSQ entries");
        assert!(csq_entries[0].contains("missense_variant"));
        assert!(csq_entries[1].contains("intron_variant"));
    }

    /// A VCF `SYMBOL` has each ASCII whitespace character written `%20`
    /// (OutputFactory.pm:1527 at release/116.2, before the subfield encoding
    /// turns whitespace into `_`), unless `--no_escape`, under which the
    /// subfield encoding's `_` stands; every other subfield keeps `_`, and a
    /// Unicode space, outside Perl's `\s` on a byte string, stays at both sites.
    #[test]
    fn symbol_whitespace_escapes_percent20() {
        let mut v = annotated_variant();
        v.transcript_consequences[0].gene_symbol = Some("GENE 1\t X\u{A0}Y".into());
        v.transcript_consequences[0]
            .plugin_data
            .insert("Note".to_string(), "a b".to_string());
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            vcf_options(),
            vec![("Note".to_string(), "a plugin field".to_string())],
            None,
        );
        let line = formatter.format_record(&[&v]).unwrap().remove(0);
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[3], "GENE%201%20%20X\u{A0}Y", "{csq}");
        assert_eq!(*fields.last().unwrap(), "a_b", "{csq}");

        let unescaped = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                no_escape: true,
                ..vcf_options()
            },
            Vec::new(),
            None,
        );
        let line = unescaped.format_record(&[&v]).unwrap().remove(0);
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        assert_eq!(csq.split('|').nth(3).unwrap(), "GENE_1_X\u{A0}Y", "{csq}");
    }

    /// The cache's `-` sentinel is an empty subfield, as Ensembl VEP never sets
    /// the key (OutputFactory.pm:1534-1537), and a comma-separated accession
    /// list joins with `&` inside the entry.
    #[test]
    fn dash_sentinel_is_an_empty_subfield_and_lists_join_with_ampersand() {
        let mut v = annotated_variant();
        v.transcript_consequences[0].ccds = Some("-".into());
        v.transcript_consequences[0].refseq = Some("NM_1.1,NM_2.2".into());
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                ccds: true,
                xref_refseq: true,
                ..vcf_options()
            },
            Vec::new(),
            None,
        );
        let line = formatter.format_record(&[&v]).unwrap().remove(0);
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        let fields: Vec<&str> = csq.split('|').collect();
        let at = |name: &str| formatter.fields.iter().position(|f| f == name).unwrap();
        assert_eq!(fields[at("CCDS")], "", "{csq}");
        assert_eq!(fields[at("RefSeq")], "NM_1.1&NM_2.2", "{csq}");
    }

    #[test]
    fn test_vcf_encoding() {
        assert_eq!(vcf_encode("a,b"), "a&b");
        assert_eq!(vcf_encode("a;b"), "a%3Bb");
        assert_eq!(vcf_encode("a|b"), "a&b");
        assert_eq!(vcf_encode("a%b"), "a%b");
        assert_eq!(vcf_encode("p.Pro336="), "p.Pro336=");
        assert_eq!(vcf_encode("a b"), "a_b");
        assert_eq!(vcf_encode("a\tb"), "a_b");
        assert_eq!(vcf_encode("a \t b"), "a_b");
        assert_eq!(vcf_encode(" \na\r\n"), "_a_");
        assert_eq!(vcf_encode("a\x0B\x0Cb"), "a_b");
        assert_eq!(vcf_encode("a\u{A0}b"), "a\u{A0}b");
        assert_eq!(vcf_encode("normal"), "normal");
        // A control character outside Perl's `\s` is percent-encoded; the
        // whitespace rule keeps the tab as `_`.
        assert_eq!(vcf_encode("rs\x001\tx\x1Fy\x7F"), "rs%001_x%1Fy%7F");
        assert_eq!(vcf_encode("a\u{85}b"), "a\u{85}b");
    }

    /// A CSQ value carrying a NUL, a tab and another control character is
    /// encoded: the NUL and the `0x01` as `%00` and `%01`, the tab by the
    /// whitespace rule (`_`, or `%20` in a symbol); the record's own columns,
    /// the ID among them, stay the input line's bytes.
    #[test]
    fn control_bytes_in_a_csq_value_are_percent_encoded() {
        let mut v = annotated_variant();
        v.transcript_consequences[0].gene_symbol = Some("GEN\0E\tX".into());
        v.transcript_consequences[0]
            .plugin_data
            .insert("Note".to_string(), "a\0b\tc\x01".to_string());
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            vcf_options(),
            vec![("Note".to_string(), "a plugin field".to_string())],
            None,
        );
        let line = formatter.format_record(&[&v]).unwrap().remove(0);
        assert!(line.starts_with("21\t100\trs1\tA\tG\t"), "{line}");
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[3], "GEN%00E%20X", "{csq}");
        assert_eq!(*fields.last().unwrap(), "a%00b_c%01", "{csq}");
    }

    /// `OutputFactory.pm` encodes whitespace in a gene symbol as `%20` for VCF
    /// output unless `--no_escape` (`BaseTranscriptVariationAllele_to_output_hash`,
    /// the `symbol` block), before `output_hash_to_vcf_info_chunk` turns the
    /// whitespace of every other value into `_` (`VCF.pm`); under `--no_escape`
    /// the symbol takes that `_` rule like any other value.
    #[test]
    fn symbol_whitespace_is_percent20_in_vcf_unless_no_escape() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            gene_symbol: Some("MY GENE".into()),
            gene_symbol_source: Some("HGNC Symbol".into()),
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            feature_type: FeatureType::Transcript,
            strand: 1,
            ..Default::default()
        };
        let fields = ["SYMBOL".to_string(), "SYMBOL_SOURCE".to_string()];
        let escaped = FieldOptions {
            symbol: true,
            ..Default::default()
        };
        assert_eq!(
            format_csq_entry(&fields, &[true, true], &variant, Some(&tc), &escaped),
            "MY%20GENE|HGNC_Symbol"
        );
        let plain = FieldOptions {
            symbol: true,
            no_escape: true,
            ..Default::default()
        };
        assert_eq!(
            format_csq_entry(&fields, &[true, true], &variant, Some(&tc), &plain),
            "MY_GENE|HGNC_Symbol"
        );
    }

    #[test]
    fn multi_allelic_record_is_one_line_with_every_allele_in_csq() {
        let formatter = VcfOutputFormatter::default();
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.raw_input = Some("21\t100\t.\tA\tG,T\t.\t.\t.".into());
        a.allele_index = 0;
        a.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            ..Default::default()
        }];
        let mut b = a.clone();
        b.alt_alleles = vec![b"T".to_vec()];
        b.allele_string = "A/T".into();
        b.allele_index = 1;
        b.transcript_consequences[0].consequences = smallvec![Consequence::SynonymousVariant];
        b.transcript_consequences[0].impact = Impact::LOW;
        let lines = formatter.format_record(&[&a, &b]).unwrap();
        assert_eq!(lines.len(), 1);
        let info = lines[0].split('\t').nth(7).unwrap();
        let entries: Vec<&str> = info.strip_prefix("CSQ=").unwrap().split(',').collect();
        assert_eq!(entries.len(), 2);
        assert!(entries[0].starts_with("G|missense_variant|MODERATE|"));
        assert!(entries[1].starts_with("T|synonymous_variant|LOW|"));
    }

    #[test]
    fn an_existing_csq_entry_is_replaced() {
        let formatter = VcfOutputFormatter::default();
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tDP=3;CSQ=old|stuff".into());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let lines = formatter.format_record(&[&v]).unwrap();
        let info = lines[0].split('\t').nth(7).unwrap();
        assert!(
            info.starts_with("DP=3;CSQ=G|intergenic_variant|MODIFIER|"),
            "{info}"
        );
        assert_eq!(info.matches("CSQ=").count(), 1);
    }

    /// A record carried through without consequences loses a stale `CSQ=` as an
    /// annotated one does (OutputFactory/VCF.pm:327-331 strips the field before
    /// it looks for chunks): exactly the first entry under the name goes, whatever
    /// its value (Ensembl VEP's substitution takes the next field too when the
    /// value is empty), the rest of INFO survives, and `.` stands in when the
    /// stale field was all of it. A `vep_skip` mark reaches the same line even
    /// with a consequence attached.
    #[test]
    fn a_carried_record_drops_a_stale_csq_entry() {
        let formatter = VcfOutputFormatter::default();
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tDP=3;CSQ=old|stuff".into());
        assert_eq!(
            formatter.format_record(&[&v]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\tDP=3"]
        );

        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tCSQ=old|stuff".into());
        assert_eq!(
            formatter.format_record(&[&v]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\t."]
        );

        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tCSQ=old|stuff;DP=3".into());
        v.vep_skip = true;
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            consequences: smallvec![Consequence::MissenseVariant],
            ..Default::default()
        }];
        assert_eq!(
            formatter.format_record(&[&v]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\tDP=3"]
        );
    }

    /// Each individual's copy is one VCF line (OutputFactory/VCF.pm:310-355 writes
    /// one per feature): the raw record, sample columns included, with that copy's
    /// CSQ, whose entries carry `IND` and `ZYG`; a reference-homozygous copy with
    /// no allele to annotate is the record without consequences, INFO `.` when it
    /// had none.
    #[test]
    fn individual_writes_one_line_per_individual_per_record() {
        let formatter = VcfOutputFormatter::new(
            "CSQ".into(),
            FieldOptions {
                individual: true,
                ..Default::default()
            },
            Vec::new(),
            Some(vec![
                "Allele".into(),
                "Feature".into(),
                "IND".into(),
                "ZYG".into(),
            ]),
        );
        let line = "21\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0";
        let mut het = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        het.raw_input = Some(line.into());
        het.individual = Some("S1".into());
        het.genotype = vec!["A".into(), "G".into()];
        het.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            consequences: smallvec![Consequence::MissenseVariant],
            ..Default::default()
        }];
        let mut hom_ref = InputVariant::without_allele("21".into(), 100, "A", std::iter::empty());
        hom_ref.raw_input = Some(line.into());
        hom_ref.individual = Some("S2".into());
        hom_ref.genotype = vec!["A".into(), "A".into()];
        hom_ref.hom_ref = true;
        assert_eq!(
            formatter.format_record(&[&het]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\tCSQ=G|ENST1|S1|HET\tGT\t0/1\t0/0"]
        );
        assert_eq!(
            formatter.format_record(&[&hom_ref]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0"]
        );
    }

    /// A record the pipeline judged intergenic (`most_severe_consequence` set,
    /// no transcript consequence) gets one intergenic CSQ entry; a record with
    /// neither a consequence nor a verdict is carried through without a CSQ
    /// field, as OutputFactory/VCF.pm:341-353 writes a record it has no chunk for.
    #[test]
    fn test_format_intergenic_variant() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25000100\t.\tA\tG\t.\t.\t.".into());

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(lines, ["21\t25000100\t.\tA\tG\t.\t.\t."]);

        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(lines.len(), 1);

        let fields: Vec<&str> = lines[0].split('\t').collect();
        let info = fields[7];
        assert!(info.starts_with("CSQ="));
        assert!(info.contains("intergenic_variant"));
        assert!(info.contains("MODIFIER"));
    }

    #[test]
    fn test_csq_appended_to_existing_info() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25000100\t.\tA\tG\t100\tPASS\tDP=50".into());
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);

        let lines = formatter.format_variant(&variant).unwrap();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        let info = fields[7];
        assert!(
            info.starts_with("DP=50;CSQ="),
            "Should append CSQ to existing INFO, got: {}",
            info
        );
    }

    /// With `raw_input` present, QUAL/FILTER/INFO from the source record survive.
    ///
    /// This is the direction the capture gate must be set to: `-o vcf` re-emits each
    /// original record with CSQ appended, so columns 6-8 come from the input line.
    /// Paired with `test_format_variant_without_raw_input_drops_qual_filter_info`
    /// below, which pins what is lost when the capture is off; neither test needs
    /// a JSON cache, so unlike the `-o vcf` integration test neither can skip.
    #[test]
    fn test_format_variant_with_raw_input_preserves_qual_filter_info() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25000100\trs99\tA\tG\t99\tPASS\tDP=30".into());
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);

        let lines = formatter.format_variant(&variant).unwrap();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields[2], "rs99", "ID must survive");
        assert_eq!(fields[5], "99", "QUAL must survive");
        assert_eq!(fields[6], "PASS", "FILTER must survive");
        assert!(
            fields[7].starts_with("DP=30;CSQ="),
            "the original INFO must be preserved with CSQ appended, got: {}",
            fields[7]
        );
    }

    /// Without `raw_input`, the reconstructed record loses QUAL, FILTER, and the
    /// original INFO, silently.
    ///
    /// `runner::needs_raw_input_capture` gates the capture on the output format.
    /// Setting that gate wrong for `vcf` produces exactly this line: a well-formed
    /// VCF record, with a valid `##INFO` header and a populated `CSQ=`, that has
    /// quietly replaced the caller's QUAL/FILTER/INFO with `.`. Nothing errors and
    /// the exit code stays 0, so this test pins the fallback shape as the thing the
    /// gate must prevent.
    #[test]
    fn test_format_variant_without_raw_input_drops_qual_filter_info() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.id = Some("rs99".into());
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        assert!(variant.raw_input.is_none(), "capture is off in this case");

        let lines = formatter.format_variant(&variant).unwrap();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields.len(), 8, "reconstructed record has all 8 columns");
        assert_eq!(fields[0], "21");
        assert_eq!(fields[1], "25000100");
        assert_eq!(fields[2], "rs99");
        assert_eq!(fields[5], ".", "QUAL is lost when raw_input is absent");
        assert_eq!(fields[6], ".", "FILTER is lost when raw_input is absent");
        assert!(
            fields[7].starts_with("CSQ="),
            "the original INFO is lost; only CSQ remains, got: {}",
            fields[7]
        );
        // And it fails silently: the line is well-formed and CSQ is populated, which
        // is why header-plus-`CSQ=` assertions cannot detect a wrong gate.
        assert!(fields[7].len() > "CSQ=".len(), "CSQ is populated");
    }

    /// A reference of one run of bases on one chromosome, for the anchor base.
    struct Bases {
        chr: &'static str,
        start: u64,
        bases: &'static [u8],
    }

    impl ReferenceSequence for Bases {
        fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
            if chr != self.chr || start < self.start || end < start {
                return None;
            }
            let from = usize::try_from(start - self.start).ok()?;
            let to = usize::try_from(end - self.start + 1).ok()?;
            self.bases.get(from..to).map(<[u8]>::to_vec)
        }

        fn record_name(&self, chr: &str) -> Option<String> {
            (chr == self.chr).then(|| chr.to_string())
        }

        fn sequence_length(&self, chr: &str) -> Option<u64> {
            (chr == self.chr).then(|| self.start + self.bases.len() as u64 - 1)
        }
    }

    /// A record of a non-VCF input: `alleles` is the record's allele string.
    fn text_record(start: u64, end: u64, alleles: &str, id: Option<&str>) -> InputVariant {
        let (reference, alt) = alleles.split_once('/').unwrap();
        let mut v = InputVariant::new(
            "21".into(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        v.id = id.map(str::to_string);
        v
    }

    /// A structural record of a non-VCF input: `word` is the class word written
    /// and `term` the SO term the parser gave it, `None` for a word the table
    /// lacks.
    fn structural_record(
        start: u64,
        end: u64,
        word: &str,
        term: Option<&str>,
        id: &str,
    ) -> InputVariant {
        let mut v = InputVariant::new("21".into(), start, end, b"N".to_vec(), word.into());
        v.is_structural = true;
        v.record_so_term = term.map(str::to_string);
        v.id = Some(id.to_string());
        v
    }

    fn line_of(variant: &InputVariant, reference: Option<&dyn ReferenceSequence>) -> String {
        String::from_utf8(to_vcf_record(variant, reference)).unwrap()
    }

    /// `VariationFeature.pm` `to_VCF_record` (release 116, lines 2385-2430) on an
    /// indel: the alleles differ in length, so each takes the base before the
    /// record's start and the position moves one base back; without a reference
    /// that base is `N` (`BaseVariationFeature.pm:345`). Ensembl VEP 116.2's VCF
    /// output for the `GRCh38-input-ensembl` corpus (no FASTA) writes
    /// `21 30289183 stop_ins1_lost_rev N NG` for the insertion `-/G` at
    /// 30289184-30289183 and `21 30348409 21_30348410_ta/- NTA N` for the
    /// deletion the parser minimised to `TA/-` at 30348410-30348411; the
    /// `GRCh38-input-hgvs` corpus (FASTA) writes `21 30282168 ... TC TAT` for
    /// `21:g.30282169delCinsAT`, the `C/AT` at 30282169, and `21 30282170 ...
    /// AATA A` for `21:g.30282171_30282173del`, the `ATA/-` at 30282171-30282173.
    #[test]
    fn to_vcf_record_indel_reanchored_pos_minus_one() {
        let insertion = text_record(30289184, 30289183, "-/G", Some("stop_ins1_lost_rev"));
        assert_eq!(
            line_of(&insertion, None),
            "21\t30289183\tstop_ins1_lost_rev\tN\tNG\t.\t.\t."
        );
        let deletion = text_record(30348410, 30348411, "TA/-", Some("21_30348410_ta/-"));
        assert_eq!(
            line_of(&deletion, None),
            "21\t30348409\t21_30348410_ta/-\tNTA\tN\t.\t.\t."
        );
        // Chromosome 21 from 30282168: TCAATA (GRCh38).
        let genome = Bases {
            chr: "21",
            start: 30282168,
            bases: b"TCAATA",
        };
        let delins = text_record(30282169, 30282169, "C/AT", Some("21:g.30282169delCinsAT"));
        assert_eq!(
            line_of(&delins, Some(&genome)),
            "21\t30282168\t21:g.30282169delCinsAT\tTC\tTAT\t.\t.\t."
        );
        let deletion = text_record(
            30282171,
            30282173,
            "ATA/-",
            Some("21:g.30282171_30282173del"),
        );
        assert_eq!(
            line_of(&deletion, Some(&genome)),
            "21\t30282170\t21:g.30282171_30282173del\tAATA\tA\t.\t.\t."
        );
    }

    /// A name and a chromosome carrying a tab, a NUL and a line end reach the
    /// synthesised line percent-encoded, so the line keeps its eight columns.
    #[test]
    fn to_vcf_record_percent_encodes_control_bytes_in_its_columns() {
        let mut record = text_record(100, 100, "C/A", Some("rs\x001\tx\ny"));
        record.original_chr = "2\x7F1".to_string();
        let line = line_of(&record, None);
        assert_eq!(line, "2%7F1\t100\trs%001%09x%0Ay\tC\tA\t.\t.\t.");
        assert_eq!(line.split('\t').count(), 8);
    }

    /// Lines 2399-2402: alleles that all open with the same base need no anchor
    /// and keep the record's position, the shape of a record written with its
    /// anchor base (a multi-allelic `TC/T/TG` left unminimised); lines 2432-2441:
    /// alleles of one length are a balanced record at its position, every
    /// alternate comma-joined. Ensembl VEP 116.2 writes `21 30282169 multi C A,T`
    /// for the `GRCh38-input-ensembl` corpus's `C/A/T` and names an ID-less
    /// record `chr_start_alleles` (`21_30282169_C/A`, Parser.pm:462).
    #[test]
    fn to_vcf_record_shared_first_base_and_balanced_alleles() {
        let mut anchored = text_record(100, 101, "TC/T", None);
        anchored.record_allele_string_multi = Some("TC/T/TG".into());
        assert_eq!(
            line_of(&anchored, None),
            "21\t100\t21_100_TC/T/TG\tTC\tT,TG\t.\t.\t."
        );
        let mut multi = text_record(30282169, 30282169, "C/A", Some("multi"));
        multi.record_allele_string_multi = Some("C/A/T".into());
        assert_eq!(
            line_of(&multi, None),
            "21\t30282169\tmulti\tC\tA,T\t.\t.\t."
        );
        let unnamed = text_record(30282169, 30282169, "C/A", None);
        assert_eq!(
            line_of(&unnamed, None),
            "21\t30282169\t21_30282169_C/A\tC\tA\t.\t.\t."
        );
    }

    /// Line 2376: on the reverse strand every allele of `A`, `C`, `G`, `T` and
    /// `N` alone is reverse-complemented, a `-` left as it is. Ensembl VEP
    /// 116.2 writes `21 11000000 intergenic_rev T C` for the
    /// `GRCh38-input-ensembl` corpus's `A/G` on strand `-1`; an insertion `-/G`
    /// on that strand is `N NC` one base back.
    #[test]
    fn to_vcf_record_strand_minus_one_revcomp() {
        let mut snv = text_record(11000000, 11000000, "A/G", Some("intergenic_rev"));
        snv.strand = vep_core::coordinate::Strand::Reverse;
        assert_eq!(
            line_of(&snv, None),
            "21\t11000000\tintergenic_rev\tT\tC\t.\t.\t."
        );
        let mut insertion = text_record(30289184, 30289183, "-/G", Some("rev_ins"));
        insertion.strand = vep_core::coordinate::Strand::Reverse;
        assert_eq!(
            line_of(&insertion, None),
            "21\t30289183\trev_ins\tN\tNC\t.\t.\t."
        );
        // The CSQ `Allele` of the record is reverse-complemented too
        // (`OutputFactory/VCF.pm:392-393`): Ensembl VEP 116.2's entry for
        // `intergenic_rev` opens `C|intergenic_variant|`, where the default
        // format's Allele column prints `G`.
        snv.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let formatter = VcfOutputFormatter::default().with_vcf_input(false);
        let line = formatter.format_record(&[&snv]).unwrap().remove(0);
        assert!(
            line.starts_with(
                "21\t11000000\tintergenic_rev\tT\tC\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|"
            ),
            "{line}"
        );
        assert_eq!(
            info_of(&line).matches('|').count(),
            formatter.fields.len() - 1
        );
    }

    /// `StructuralVariationFeature.pm` `to_VCF_record` (release 116, lines
    /// 1149-1171): the base before the span and the position one base back,
    /// the class as `<DEL>`, `<DUP>`, `<INS>` or `<TDUP>` and any other class
    /// term inside the brackets, and `END=` as INFO. Ensembl VEP 116.2 writes
    /// `21 30213500 del_from_before_transcript N <DEL> . . END=30214300` and
    /// `21 30215900 dup_past_transcript_end N <DUP> . . END=30216600` for the
    /// `GRCh38-input-ensembl` corpus's `DEL` and `DUP` records (no FASTA).
    #[test]
    fn to_vcf_record_sv_symbolic() {
        let deletion = structural_record(
            30213501,
            30214300,
            "DEL",
            Some("deletion"),
            "del_from_before_transcript",
        );
        assert_eq!(
            line_of(&deletion, None),
            "21\t30213500\tdel_from_before_transcript\tN\t<DEL>\t.\t.\tEND=30214300"
        );
        let duplication = structural_record(
            30215901,
            30216600,
            "DUP",
            Some("duplication"),
            "dup_past_transcript_end",
        );
        assert_eq!(
            line_of(&duplication, None),
            "21\t30215900\tdup_past_transcript_end\tN\t<DUP>\t.\t.\tEND=30216600"
        );
        let inversion = structural_record(1001, 2000, "INV", Some("inversion"), "inv");
        let genome = Bases {
            chr: "21",
            start: 1000,
            bases: b"G",
        };
        assert_eq!(
            line_of(&inversion, Some(&genome)),
            "21\t1000\tinv\tG\t<inversion>\t.\t.\tEND=2000"
        );
        let unsupported = structural_record(1001, 2000, "FOO", None, "foo");
        assert_eq!(
            line_of(&unsupported, None),
            "21\t1000\tfoo\tN\t<FOO>\t.\t.\tEND=2000"
        );
    }

    /// Lines 2409-2419: an allele that is not bases alone is written inside `<`
    /// and `>`; `*` alleles are left out of the shared-first-base test (line
    /// 2399), so `A/*` keeps its position and `A` its one base. An allele `0`
    /// is false to Perl (line 2383), so `A/0` is a balanced record, and its
    /// joined ALT `0` is false too (lines 2427 and 2440), so the ALT column is
    /// `.`.
    #[test]
    fn to_vcf_record_writes_a_symbolic_allele_in_brackets() {
        let spanning = text_record(100, 100, "A/*", Some("star"));
        assert_eq!(line_of(&spanning, None), "21\t100\tstar\tA\t<*>\t.\t.\t.");
        let zero = text_record(100, 100, "A/0", Some("zero"));
        assert_eq!(line_of(&zero, None), "21\t100\tzero\tA\t.\t.\t.\t.");
    }

    /// `trim_right` (ensembl-variation `Utils/Sequence.pm:1058-1096`): the last
    /// base is removed from every allele while they all share it and none is
    /// one base long; an empty allele ends the trimming.
    #[test]
    fn trim_right_follows_perl() {
        let mut shared = vec![b"ACGT".to_vec(), b"AGT".to_vec()];
        trim_right(&mut shared);
        assert_eq!(shared, [b"AC".to_vec(), b"A".to_vec()]);
        let mut different = vec![b"AC".to_vec(), b"AG".to_vec()];
        trim_right(&mut different);
        assert_eq!(different, [b"AC".to_vec(), b"AG".to_vec()]);
        let mut one_base = vec![b"T".to_vec(), b"GT".to_vec()];
        trim_right(&mut one_base);
        assert_eq!(one_base, [b"T".to_vec(), b"GT".to_vec()]);
        let mut empty = vec![Vec::new(), b"TT".to_vec()];
        trim_right(&mut empty);
        assert_eq!(empty, [Vec::new(), b"TT".to_vec()]);
        let mut alone = vec![b"AA".to_vec()];
        trim_right(&mut alone);
        assert_eq!(alone, [b"AA".to_vec()]);
    }

    /// A record built into a VCF line gets its consequences as a VCF record
    /// does: `CSQ=` alone as INFO (`OutputFactory/VCF.pm:338, 345-347`), after
    /// a structural record's `END=` with a `;`.
    #[test]
    fn built_line_carries_csq_after_its_info() {
        let formatter = VcfOutputFormatter::default().with_vcf_input(false);
        let mut v = text_record(30282169, 30282169, "C/A", Some("snv_fwd"));
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let line = formatter.format_record(&[&v]).unwrap().remove(0);
        assert!(
            line.starts_with("21\t30282169\tsnv_fwd\tC\tA\t.\t.\tCSQ=A|intergenic_variant|"),
            "{line}"
        );
        let mut sv = structural_record(30213501, 30214300, "DEL", Some("deletion"), "del");
        sv.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let line = formatter.format_record(&[&sv]).unwrap().remove(0);
        assert!(
            line.starts_with(
                "21\t30213500\tdel\tN\t<DEL>\t.\t.\tEND=30214300;CSQ=deletion|intergenic_variant|"
            ),
            "{line}"
        );
    }

    /// A VCF input's line that is not UTF-8 is written back byte for byte, as
    /// Ensembl VEP writes the bytes it read (`Parser.pm:235-242` opens the input
    /// with no decoding layer): the Latin-1 `é` of an ID stays the byte `0xE9`
    /// in the bytes the pipeline writes, and the text form of the line carries
    /// the replacement character.
    #[test]
    fn a_line_that_is_not_utf8_is_written_as_read() {
        let formatter = VcfOutputFormatter::default();
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let bytes: &[u8] = b"21\t100\tid_\xE9\tA\tG\t.\t.\tDP=3";
        v.raw_input = Some(String::from_utf8_lossy(bytes).into());
        v.raw_input_bytes = Some(Arc::from(bytes));
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let line = formatter.format_record_bytes(&[&v]).unwrap().remove(0);
        assert!(
            line.starts_with(b"21\t100\tid_\xE9\tA\tG\t.\t.\tDP=3;CSQ=G|"),
            "{line:?}"
        );
        let text = formatter.format_record(&[&v]).unwrap().remove(0);
        assert!(
            text.starts_with("21\t100\tid_\u{FFFD}\tA\tG\t.\t.\tDP=3;CSQ=G|"),
            "{text}"
        );
    }
}
