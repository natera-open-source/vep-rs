// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SPDI input (Sequence, Position, Deletion, Insertion: NCBI's variation
//! notation), read as `Parser/SPDI.pm` `create_VariationFeatures`
//! (release/116.2, lines 133-208) reads it over ensembl-variation's
//! `VariationFeatureAdaptor::fetch_by_spdi_notation` (release/116,
//! lines 2477-2611):
//!
//! ```text
//! NC_000021.9:25585732:A:G      substitution of the base after position 25585732
//! 21:25585655:ACGTG:            deletion by sequence
//! 21:25585655:5:                deletion by length
//! chr21:25585732::TT            insertion between two bases
//! 21:25585732:2:GG              substitution by length
//! ```
//!
//! One notation per line, its whitespace removed (`:148`), the position
//! 0-based, so the variant starts at the position plus one (`:2524`). An empty
//! or `0` deletion with inserted bases is an insertion ending at the position
//! (`:2545-2549`); a deletion by letters is checked against the reference and
//! one by length read from it (`:2552-2567`); a substitution or indel checks
//! its letters or reads a run of the length given, which the inserted bases
//! must match in length (`:2570-2596`). Every allele is upper-cased. The
//! sequence name resolves to the annotation source's chromosome through the
//! run's synonyms, as Ensembl's slice adaptor resolves an accession, and the
//! record prints under that name (`:196`) while its name and its re-emitted
//! line are the notation as written (`:193, 198`). A version is stripped from
//! an `LRG_` or `NT_` name alone (`:2500-2502`). A record Ensembl VEP's
//! post-processing minimises (shared bases at either end of a bi-allelic
//! indel) is minimised here too ([`minimise_alleles`]).
//!
//! A notation the adaptor refuses is a line the run warns about and skips
//! (`:177-184`); the warning carries Ensembl VEP's own text for the shape. A
//! line with no word character counts as empty (`Parser.pm:273`); a `#` line is
//! a record, as ensembl-io's list parser has no comment lines. The reference
//! bases come from the run's FASTA, so the format needs one; a span past the
//! end of its sequence is refused with a warning naming the span and the
//! sequence's length, where Ensembl VEP's slice reads `N` for every base past
//! the end (Utils/FastaSequence.pm:619-629) and the record is annotated; and a
//! deletion of length 0 is refused with a warning naming the empty span, where
//! Ensembl VEP keeps the feature `/-` (`:2563-2566`) and writes no row for it.

use std::collections::HashSet;
use std::io::BufRead;
use std::sync::Arc;

use vep_core::coordinate::parse_coordinate;
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::{
    has_word_character, is_perl_space, minimise_alleles, InputParser, ParserContext,
    ReferenceSequence,
};

/// Parser for SPDI notation, one notation per line.
pub struct SpdiParser<R: BufRead> {
    reader: R,
    fasta: Option<Arc<dyn ReferenceSequence>>,
    synonyms: Arc<ChromosomeSynonyms>,
    chromosomes: Arc<HashSet<String>>,
    max_sv_size: Option<u64>,
    line_buf: Vec<u8>,
    /// Empty lines read since the last [`InputParser::take_empty_lines`].
    empty_lines: u64,
}

/// The parts of a notation between its three colons, whitespace removed.
struct Notation<'a> {
    sequence_id: &'a str,
    position: u64,
    deleted: &'a str,
    inserted: &'a str,
}

impl<R: BufRead> SpdiParser<R> {
    pub fn new(reader: R, context: &ParserContext) -> Self {
        Self {
            reader,
            fasta: context.fasta.clone(),
            synonyms: Arc::clone(&context.synonyms),
            chromosomes: Arc::clone(&context.chromosomes),
            max_sv_size: context.max_sv_size,
            line_buf: Vec::new(),
            empty_lines: 0,
        }
    }

    /// Parse one line into its variant, or the warning Ensembl VEP writes for it
    /// (`WARNING: Unable to parse SPDI notation '<notation>'` and the exception's
    /// message, `:183`). A span over `--max_sv_size` is then refused with the
    /// text `validate_vf` refuses it with (Parser.pm:491-499, the span being end
    /// minus start of the minimised record), here under `--dont_skip` too, where
    /// Ensembl VEP carries the feature unannotated.
    fn parse_record(&self, line: &str) -> Result<InputVariant, IoError> {
        let spdi: String = line.chars().filter(|c| !is_perl_space(*c)).collect();
        let variant = self.variant(&spdi).map_err(|reason| {
            IoError::SkippedLine(format!("Unable to parse SPDI notation '{spdi}': {reason}"))
        })?;
        let span = variant.end.saturating_sub(variant.start);
        match self.max_sv_size {
            Some(max) if span > max => Err(IoError::SkippedLine(format!(
                "variant size ({span}) is bigger than --max_sv_size ({max})"
            ))),
            _ => Ok(variant),
        }
    }

    /// `fetch_by_spdi_notation` (`:2488-2608`) over the run's reference, every
    /// `throw` text kept as written.
    fn variant(&self, spdi: &str) -> Result<InputVariant, String> {
        let notation = split_notation(spdi)?;
        let raw_sequence_id = notation.sequence_id;
        let sequence_id = strip_version(raw_sequence_id);
        let chr = self.chromosome_name(&sequence_id);
        let position = notation.position;
        let start = position + 1;
        // The slice is fetched before the shape is read (`:2533-2537`), so a
        // sequence the reference lacks fails every shape with the same text
        // (`:2536`, naming the sequence as written).
        let not_valid = || format!("Sequence name {raw_sequence_id} not valid");
        let fasta = self.fasta.as_ref().ok_or_else(not_valid)?;
        let sequence_length = fasta.sequence_length(&chr).ok_or_else(not_valid)?;
        // Ensembl VEP has no refusal here: its slice reads `N` for every base past
        // the end (Utils/FastaSequence.pm:619-629) and annotates the record.
        let past_end = |from: u64, to: u64| {
            format!(
                "Span {}-{} lies past the end of sequence {raw_sequence_id} \
                 ({sequence_length} bases)",
                from.min(to),
                from.max(to)
            )
        };

        let deleted = notation.deleted;
        let inserted = notation.inserted;
        let no_deletion = deleted.is_empty() || deleted == "0";
        let no_insertion = inserted.is_empty() || inserted == "0";
        let deleted_is_letters = is_bases(deleted);
        let deleted_is_length = is_digits(deleted);
        let inserted_is_letters = is_bases(inserted);
        // The reference bases from `start` to `end`, as `get_reference_allele`
        // reads them (`:2614-2627`).
        let reference = |end: u64| {
            fasta.sequence(&chr, start, end).ok_or_else(|| {
                if end > sequence_length {
                    past_end(start, end)
                } else {
                    not_valid()
                }
            })
        };
        let unsupported =
            || format!("Could not parse the SPDI notation {spdi}. SPDI notation not supported.");

        if deleted == "0" && inserted == "0" {
            return Err(unsupported());
        }
        let (end, ref_allele, alt_allele) = if no_deletion && inserted_is_letters {
            // An insertion ends at the position, before its start (`:2545-2549`),
            // and reads no base: one after the last base is within the sequence,
            // one further on is past its end.
            if position > sequence_length {
                return Err(past_end(position, start));
            }
            let alt_allele = inserted.to_ascii_uppercase().into_bytes();
            (position, b"-".to_vec(), alt_allele)
        } else if (deleted_is_letters || deleted_is_length) && (no_insertion || inserted_is_letters)
        {
            // A deletion (`:2552-2567`), a substitution or an indel (`:2570-2596`):
            // the reference allele is the deleted letters, checked against the
            // reference, or a run of the length given, read from it.
            let (end, ref_allele) = if deleted_is_letters {
                let end = position + deleted.len() as u64;
                let ref_allele = deleted.to_ascii_uppercase().into_bytes();
                check_reference(spdi, &sequence_id, start, end, &ref_allele, reference(end)?)?;
                (end, ref_allele)
            } else {
                let length = parse_coordinate(deleted)
                    .ok_or_else(|| format!("Could not parse the SPDI notation {spdi}"))?;
                if inserted_is_letters && inserted.len() as u64 != length {
                    return Err(format!(
                        "Could not parse the SPDI notation {spdi}. Deleted sequence length \
                         ({deleted}) does not match inserted sequence length ({}).",
                        inserted.len()
                    ));
                }
                // Ensembl VEP keeps a deletion of length 0 as the feature `/-` over
                // the empty span (`:2563-2566`; the slice's sequence is empty when
                // its start is its end plus one), a feature that writes no
                // consequence row. The span is refused by name here, before
                // `reference` reads it as a sequence the FASTA lacks.
                if length == 0 {
                    return Err(format!(
                        "Deleted sequence length ({deleted}) leaves nothing to delete: span \
                         {start}-{position} is empty"
                    ));
                }
                let end = position + length;
                (end, reference(end)?)
            };
            let alt_allele = if no_insertion {
                b"-".to_vec()
            } else {
                let alt_allele = inserted.to_ascii_uppercase().into_bytes();
                if ref_allele == alt_allele {
                    let allele = String::from_utf8_lossy(&ref_allele);
                    return Err(format!(
                        "Reference allele given by SPDI notation {spdi} ({allele}) matches alt \
                         allele given by SPDI notation {spdi} ({allele})"
                    ));
                }
                alt_allele
            };
            (end, ref_allele, alt_allele)
        } else {
            return Err(unsupported());
        };

        let mut variant = InputVariant::new(chr, start, end, ref_allele, alt_allele);
        variant.uploaded_allele_string = Some(variant.allele_string.as_str().into());
        variant.id = Some(spdi.to_string());
        variant.raw_input = Some(spdi.into());
        minimise_alleles(&mut variant);
        Ok(variant)
    }

    /// The chromosome the annotation source knows `sequence_id` by: the name
    /// itself when the source holds it, else a synonym it holds, else its `chr`
    /// spelling (`BaseVEP.pm:689-738`), each tried again regardless of case,
    /// since the database collation Ensembl's slice adaptor resolves a name
    /// under ignores it; a name that resolves nowhere is kept as written. Among
    /// several names differing only in case, the first in byte order wins, so
    /// the choice does not depend on the order the tables are iterated in.
    fn chromosome_name(&self, sequence_id: &str) -> String {
        let holds = |name: &str| self.chromosomes.contains(name);
        let resolved = self.synonyms.source_chr_name(sequence_id, holds);
        if holds(&resolved) {
            return resolved;
        }
        let mut same_letters: Vec<&str> = self
            .chromosomes
            .iter()
            .map(String::as_str)
            .chain(self.synonyms.names())
            .filter(|name| name.eq_ignore_ascii_case(sequence_id))
            .collect();
        same_letters.sort_unstable();
        same_letters
            .into_iter()
            .map(|name| self.synonyms.source_chr_name(name, holds))
            .find(|name| holds(name))
            .unwrap_or(resolved)
    }
}

/// The notation split at its colons (`:2489-2518`): exactly three, a sequence
/// name, a position of digits, and a deletion and an insertion each empty or a
/// run of word characters.
fn split_notation(spdi: &str) -> Result<Notation<'_>, String> {
    let parts: Vec<&str> = spdi.split(':').collect();
    let [sequence_id, position, deleted, inserted] = parts[..] else {
        let reason = if parts.len() > 4 {
            "Too many elements present"
        } else {
            "Too few elements present"
        };
        return Err(format!(
            "Could not parse the SPDI notation {spdi}. {reason}"
        ));
    };
    let well_formed = !sequence_id.is_empty()
        && is_digits(position)
        && (deleted.is_empty() || is_word(deleted))
        && (inserted.is_empty() || is_word(inserted));
    let could_not_parse = || format!("Could not parse the SPDI notation {spdi}");
    if !well_formed {
        return Err(could_not_parse());
    }
    let position = parse_coordinate(position).ok_or_else(could_not_parse)?;
    Ok(Notation {
        sequence_id,
        position,
        deleted,
        inserted,
    })
}

/// The reference check of `:2559-2560` and `:2579-2580`: the deleted letters
/// must be the reference's bases over the span, the message naming the
/// sequence as the lookup used it (version stripped) and the bases of both.
fn check_reference(
    spdi: &str,
    sequence_id: &str,
    start: u64,
    end: u64,
    given: &[u8],
    reference: Vec<u8>,
) -> Result<(), String> {
    if reference == given {
        return Ok(());
    }
    Err(format!(
        "Reference allele extracted from {sequence_id}:{start}-{end} ({}) does not match \
         reference allele given by SPDI notation {spdi} ({})",
        String::from_utf8_lossy(&reference),
        String::from_utf8_lossy(given)
    ))
}

/// `:2498-2502`: every version (`.<digits>`) is removed from a sequence name
/// beginning `LRG` or `NT` in any case, since the database names those without
/// one; an `NC_` accession keeps its version, which names the assembly.
fn strip_version(sequence_id: &str) -> String {
    let lower = sequence_id.to_ascii_lowercase();
    if !(lower.starts_with("lrg") || lower.starts_with("nt")) {
        return sequence_id.to_string();
    }
    let mut parts = sequence_id.split('.');
    let mut stripped = parts.next().unwrap_or_default().to_string();
    for part in parts {
        let digits = part.bytes().take_while(u8::is_ascii_digit).count();
        if digits == 0 {
            stripped.push('.');
        }
        stripped.push_str(&part[digits..]);
    }
    stripped
}

/// Perl's `^[0-9]+\z` on a non-empty string.
fn is_digits(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit())
}

/// Perl's `^[ACGTN]+$`, case-insensitive, on a non-empty string.
fn is_bases(text: &str) -> bool {
    !text.is_empty()
        && text
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// Perl's `^\w+$` on a non-empty string of bytes: ASCII letters, digits and `_`.
fn is_word(text: &str) -> bool {
    !text.is_empty() && text.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
}

impl<R: BufRead + Send> InputParser for SpdiParser<R> {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        loop {
            self.line_buf.clear();
            match self.reader.read_until(b'\n', &mut self.line_buf) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(IoError::Io(e))),
            }
            let line = String::from_utf8_lossy(&self.line_buf);
            // `skip_empty_lines` reads past a record with no word character
            // (Parser.pm:273), so `:::` counts as an empty line too.
            if !has_word_character(&line) {
                self.empty_lines += 1;
                continue;
            }
            return Some(self.parse_record(&line));
        }
    }

    fn take_empty_lines(&mut self) -> u64 {
        std::mem::take(&mut self.empty_lines)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::io::Cursor;
    use vep_core::coordinate::Strand;
    use vep_core::variant::VariantClass;

    /// A reference of named sequences, each a string of bases from position 1.
    struct Bases(HashMap<&'static str, &'static str>);

    impl ReferenceSequence for Bases {
        fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
            let bases = self.0.get(chr)?;
            if start == 0 || end < start || end as usize > bases.len() {
                return None;
            }
            Some(bases.as_bytes()[start as usize - 1..end as usize].to_vec())
        }

        fn record_name(&self, chr: &str) -> Option<String> {
            self.0.contains_key(chr).then(|| chr.to_string())
        }

        fn sequence_length(&self, chr: &str) -> Option<u64> {
            self.0.get(chr).map(|bases| bases.len() as u64)
        }
    }

    /// Chromosome 21 as `TGCATACGTTAGGCATGACC...`, 38 bases, the cache holding
    /// `21` and `MT`, and the synonym lines of a release-116 cache for
    /// chromosome 21.
    fn context() -> ParserContext {
        let mut reference = HashMap::new();
        reference.insert("21", "TGCATACGTTAGGCATGACCTTAGCCGATTACAGGTTC");
        reference.insert("MT", "GATCACAGGTCTATCACCCT");
        ParserContext {
            fasta: Some(Arc::new(Bases(reference))),
            synonyms: Arc::new(ChromosomeSynonyms::parse(
                "21\tNC_000021.9\n21\tchr21\n21\tCM000683.2\nMT\tNC_012920.1\nMT\tchrM\n",
            )),
            chromosomes: Arc::new(["21", "MT"].into_iter().map(String::from).collect()),
            ..ParserContext::default()
        }
    }

    fn parser_in(context: &ParserContext, data: &str) -> SpdiParser<Cursor<Vec<u8>>> {
        SpdiParser::new(Cursor::new(data.as_bytes().to_vec()), context)
    }

    fn parser(data: &str) -> SpdiParser<Cursor<Vec<u8>>> {
        parser_in(&context(), data)
    }

    /// What the parser yields for `spdi` as the input's one line.
    fn record(spdi: &str) -> Result<InputVariant, IoError> {
        parser(&format!("{spdi}\n")).next_variant().unwrap()
    }

    fn one(spdi: &str) -> InputVariant {
        record(spdi).unwrap_or_else(|e| panic!("{spdi}: {e}"))
    }

    fn refused(spdi: &str) -> String {
        let error = record(spdi).unwrap_err();
        assert!(
            matches!(error, IoError::SkippedLine(_)),
            "{spdi}: {error:?}"
        );
        error.to_string()
    }

    /// A deletion written by its bases (`:2554-2561`) is checked against the
    /// reference over `position + 1 .. position + length`, with an empty or `0`
    /// insertion either way; `21:3:ATAC:` deletes bases 4 to 7, `ATAC` of the
    /// reference, and `21:3:ATAG:` fails with Ensembl VEP's mismatch text naming
    /// the span and both alleles. The record is named by its notation.
    #[test]
    fn del_sequence_checked() {
        for spdi in ["21:3:ATAC:", "21:3:ATAC:0", "21:3:atac:"] {
            let v = one(spdi);
            assert_eq!((v.start, v.end), (4, 7), "{spdi}");
            assert_eq!(v.allele_string, "ATAC/-", "{spdi}");
            assert_eq!(v.variant_class, VariantClass::Deletion);
            assert_eq!(v.strand, Strand::Forward);
            assert_eq!(v.uploaded_variation(), spdi);
            assert_eq!(v.variation_name(), spdi);
            assert_eq!(v.raw_input.as_deref(), Some(spdi));
            assert_eq!(v.uploaded_allele_string.as_deref(), Some("ATAC/-"));
            assert!(!v.minimised);
        }
        assert_eq!(
            refused("21:3:ATAG:"),
            "Unable to parse SPDI notation '21:3:ATAG:': Reference allele extracted from \
             21:4-7 (ATAC) does not match reference allele given by SPDI notation 21:3:ATAG: (ATAG)"
        );
        // A substitution's letters are checked the same way (`:2577-2580`).
        assert_eq!(
            refused("21:3:G:A"),
            "Unable to parse SPDI notation '21:3:G:A': Reference allele extracted from \
             21:4-4 (A) does not match reference allele given by SPDI notation 21:3:G:A (G)"
        );
    }

    /// A deletion written by its length (`:2562-2565`) reads its bases from the
    /// reference: `21:3:4:` is `ATAC/-` over 4 to 7, and a substitution by
    /// length (`:2584-2590`) reads them too: `21:3:2:GG` is `AT/GG` over 4 to 5,
    /// `21:3:1:C` the SNV `A/C` at 4. The inserted bases must match the length
    /// given (`:2586-2587`), and a reference equal to the alternate allele is
    /// refused (`:2594-2595`).
    #[test]
    fn del_length_read_from_fasta() {
        let v = one("21:3:4:");
        assert_eq!((v.start, v.end), (4, 7));
        assert_eq!(v.allele_string, "ATAC/-");
        let v = one("21:3:4:0");
        assert_eq!(v.allele_string, "ATAC/-");
        let v = one("21:3:2:GG");
        assert_eq!((v.start, v.end), (4, 5));
        assert_eq!(v.allele_string, "AT/GG");
        assert_eq!(v.class_so_term(), "substitution");
        let v = one("21:3:1:C");
        assert_eq!((v.start, v.end), (4, 4));
        assert_eq!(v.allele_string, "A/C");
        assert_eq!(v.class_so_term(), "SNV");
        assert_eq!(
            refused("21:3:2:G"),
            "Unable to parse SPDI notation '21:3:2:G': Could not parse the SPDI notation \
             21:3:2:G. Deleted sequence length (2) does not match inserted sequence length (1)."
        );
        assert_eq!(
            refused("21:3:1:A"),
            "Unable to parse SPDI notation '21:3:1:A': Reference allele given by SPDI notation \
             21:3:1:A (A) matches alt allele given by SPDI notation 21:3:1:A (A)"
        );
        assert_eq!(
            refused("21:3:ATAC:atac"),
            "Unable to parse SPDI notation '21:3:ATAC:atac': Reference allele given by SPDI \
             notation 21:3:ATAC:atac (ATAC) matches alt allele given by SPDI notation \
             21:3:ATAC:atac (ATAC)"
        );
    }

    /// An empty or `0` insertion makes the record a deletion (`:2552`), as an
    /// empty or `0` deletion with inserted bases makes it an insertion with
    /// reference `-` ending at the position, so start is end plus one
    /// (`:2545-2549`); `0` for both is refused (`:2540-2542`), and so is an
    /// empty pair (`:2600`).
    #[test]
    fn ins_empty_is_deletion() {
        assert_eq!(one("21:3:A:").allele_string, "A/-");
        assert_eq!(one("21:3:A:0").allele_string, "A/-");
        for spdi in ["21:3::TT", "21:3:0:TT", "21:3::tt"] {
            let v = one(spdi);
            assert_eq!((v.start, v.end), (4, 3), "{spdi}");
            assert_eq!(v.allele_string, "-/TT", "{spdi}");
            assert_eq!(v.variant_class, VariantClass::Insertion);
            assert_eq!(v.location(), "21:3-4");
        }
        let not_supported = |spdi: &str| {
            format!(
                "Unable to parse SPDI notation '{spdi}': Could not parse the SPDI notation \
                 {spdi}. SPDI notation not supported."
            )
        };
        assert_eq!(refused("21:3:0:0"), not_supported("21:3:0:0"));
        assert_eq!(refused("21:3::"), not_supported("21:3::"));
    }

    /// The position is 0-based (`:2524`): `21:0:T:G` is the first base, start
    /// and end 1, and `21:0::G` an insertion before it, start 1 and end 0, which
    /// prints as `21:0-1` as Ensembl VEP 116.2 prints it; the reference's first
    /// base is read for the check.
    #[test]
    fn pos0_to_start() {
        let v = one("21:0:T:G");
        assert_eq!((v.start, v.end), (1, 1));
        assert_eq!(v.allele_string, "T/G");
        let v = one("21:0::G");
        assert_eq!((v.start, v.end), (1, 0));
        assert_eq!(v.location(), "21:0-1");
        assert_eq!(
            refused("21:0:A:G"),
            "Unable to parse SPDI notation '21:0:A:G': Reference allele extracted from 21:1-1 \
             (T) does not match reference allele given by SPDI notation 21:0:A:G (A)"
        );
    }

    /// The sequence name resolves through the synonyms to the chromosome the
    /// cache holds, and the record prints under that name as Ensembl VEP prints
    /// the slice's (`:196`; Ensembl VEP 116.2 prints `21:10538885` for
    /// `NC_000021.9:10538884:A:C`, `chr21:...` and `nc_000021.9:...:a:c`), while
    /// its name keeps the notation as written. A name no synonym reaches is
    /// `Sequence name <name> not valid` (`:2536`), and an insertion on it fails
    /// the same way before its shape is read.
    #[test]
    fn nc_accession_via_synonyms() {
        for spdi in [
            "NC_000021.9:3:A:C",
            "chr21:3:A:C",
            "CM000683.2:3:A:C",
            "nc_000021.9:3:a:c",
            "CHR21:3:A:C",
        ] {
            let v = one(spdi);
            assert_eq!(v.original_chr, "21", "{spdi}");
            assert_eq!(v.chr, "21", "{spdi}");
            assert_eq!(v.location(), "21:4", "{spdi}");
            assert_eq!(v.allele_string, "A/C", "{spdi}");
            assert_eq!(v.uploaded_variation(), spdi, "{spdi}");
            assert_eq!(v.display_allele(), "C", "{spdi}");
        }
        let mt = one("NC_012920.1:0:G:A");
        assert_eq!((mt.original_chr.as_str(), mt.start), ("MT", 1));
        assert_eq!(
            refused("NC_000099.1:3:A:C"),
            "Unable to parse SPDI notation 'NC_000099.1:3:A:C': Sequence name NC_000099.1 not valid"
        );
        assert_eq!(
            refused("22:3::A"),
            "Unable to parse SPDI notation '22:3::A': Sequence name 22 not valid"
        );
        // Two names differing only in case: the first in byte order (`FOO`
        // before `foo`) resolves, whichever order the tables yield them in.
        let context = ParserContext {
            synonyms: Arc::new(ChromosomeSynonyms::parse("MT\tfoo\n21\tFOO\n")),
            ..context()
        };
        let mut parser = parser_in(&context, "Foo:3:A:C\n");
        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!((v.original_chr.as_str(), v.start), ("21", 4));
    }

    /// Every refusal of `fetch_by_spdi_notation` by shape, in its order: the
    /// colon count (`:2494-2495`), then a name, a position of digits and word
    /// characters in the deletion and insertion (`:2516-2518`), then a deletion
    /// or insertion outside the four branches (`:2600`). The texts are Ensembl
    /// VEP 116.2's, as its warnings file carries them.
    #[test]
    fn bad_shape_error_text() {
        let could_not_parse = |spdi: &str, tail: &str| {
            format!("Unable to parse SPDI notation '{spdi}': Could not parse the SPDI notation {spdi}{tail}")
        };
        let cases = [
            ("21:3:A", ". Too few elements present"),
            ("21:3", ". Too few elements present"),
            ("#acommentline", ". Too few elements present"),
            ("21:3:A:C:X", ". Too many elements present"),
            ("21:3:A:C:", ". Too many elements present"),
            ("21:abc:A:C", ""),
            ("21:-5:A:C", ""),
            ("21:+5:A:C", ""),
            ("21::A:C", ""),
            (":3:A:C", ""),
            ("21:3:A-:C", ""),
            ("21:3:A:C-", ""),
            ("21:3:R:C", ". SPDI notation not supported."),
            ("21:3:2:2", ". SPDI notation not supported."),
            ("21:3:A:1", ". SPDI notation not supported."),
            ("21:3:AX:", ". SPDI notation not supported."),
            ("21:3::1", ". SPDI notation not supported."),
        ];
        for (spdi, tail) in cases {
            assert_eq!(refused(spdi), could_not_parse(spdi, tail), "{spdi}");
        }
        // The colon count is read before the shape, so a malformed position with
        // four colons is `Too many elements present`.
        assert_eq!(
            refused("21:abc:A:C:D"),
            could_not_parse("21:abc:A:C:D", ". Too many elements present")
        );
    }

    /// Shared bases at either end of a bi-allelic indel are trimmed as
    /// `post_process_vfs` trims them (Parser.pm:862-881): Ensembl VEP 116.2
    /// annotates `NC_000021.9:13630970:AATCCTGTGGGATTACCAGAAAACC:A` at
    /// 13630972-13630995 with Allele `-`; here `21:3:ATAC:A` is `TAC/-` over 5
    /// to 7, the written alleles kept for `UPLOADED_ALLELE` and the name. A
    /// substitution of one length, `21:3:AT:GG`, and an indel sharing no end,
    /// `21:3:AT:GGC`, stay as written.
    #[test]
    fn shared_bases_are_minimised() {
        let v = one("21:3:ATAC:A");
        assert_eq!((v.start, v.end), (5, 7));
        assert_eq!(v.allele_string, "TAC/-");
        assert!(v.minimised);
        assert_eq!(v.original_allele_string.as_deref(), Some("ATAC/A"));
        assert_eq!(v.uploaded_allele_string.as_deref(), Some("ATAC/A"));
        assert_eq!(v.uploaded_variation(), "21:3:ATAC:A");
        assert_eq!(v.class_so_term(), "deletion");
        let v = one("21:3:AT:GG");
        assert_eq!((v.start, v.end, v.allele_string.as_str()), (4, 5, "AT/GG"));
        assert!(!v.minimised);
        let v = one("21:3:AT:GGC");
        assert_eq!((v.start, v.end, v.allele_string.as_str()), (4, 5, "AT/GGC"));
        assert_eq!(v.class_so_term(), "indel");
    }

    /// Whitespace anywhere in the line is removed before the notation is read
    /// (`:148`), the `\r` of a `\r\n` line end with it, and the record's name and
    /// re-emitted line are the stripped notation. Letters keep their case in the
    /// name and are upper-cased in the alleles.
    #[test]
    fn whitespace_stripped_and_case_kept_in_the_name() {
        let v = one("NC_000021.9 : 3 : a : c \t\r");
        assert_eq!(v.uploaded_variation(), "NC_000021.9:3:a:c");
        assert_eq!(v.raw_input.as_deref(), Some("NC_000021.9:3:a:c"));
        assert_eq!(v.allele_string, "A/C");
        assert_eq!(v.id.as_deref(), Some("NC_000021.9:3:a:c"));
    }

    /// A version is stripped from an `LRG_` or `NT_` name, in any case, and from
    /// no other (`:2500-2502`); the mismatch text names the stripped name while
    /// `Sequence name ... not valid` names the one written (`:2536`).
    #[test]
    fn version_stripped_for_lrg_and_nt_only() {
        assert_eq!(strip_version("LRG_1.2"), "LRG_1");
        assert_eq!(strip_version("lrg_1.2"), "lrg_1");
        assert_eq!(strip_version("LRG_1.2.3"), "LRG_1");
        assert_eq!(strip_version("LRG_1.x"), "LRG_1.x");
        assert_eq!(strip_version("NT_187361.1"), "NT_187361");
        assert_eq!(strip_version("NC_000021.9"), "NC_000021.9");
        assert_eq!(strip_version("NW_003315905.1"), "NW_003315905.1");
        assert_eq!(strip_version("LRG_1"), "LRG_1");
        assert_eq!(strip_version("21"), "21");
        assert_eq!(
            refused("NT_187361.1:3:A:C"),
            "Unable to parse SPDI notation 'NT_187361.1:3:A:C': Sequence name NT_187361.1 not valid"
        );
    }

    /// A line without a word character, blank or `:::`, is an empty line counted
    /// for the caller's warning (Parser.pm:273); a `#` line is a record, which
    /// ensembl-io's list parser has no rule to skip, and fails as its text
    /// fails. The parser reads on after a refused line.
    #[test]
    fn empty_lines_counted_and_a_hash_line_is_a_record() {
        let mut parser = parser("\n  \n:::\n# note\n21:3:A:C\n\n");
        assert_eq!(parser.take_empty_lines(), 0);
        let error = parser.next_variant().unwrap().unwrap_err().to_string();
        assert!(
            error.starts_with("Unable to parse SPDI notation '#note'"),
            "{error}"
        );
        assert_eq!(parser.take_empty_lines(), 3);
        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.allele_string, "A/C");
        assert!(parser.next_variant().is_none());
        assert_eq!(parser.take_empty_lines(), 1);
    }

    /// Parser.pm:491-499 drops a record whose end minus start exceeds
    /// `--max_sv_size` (10,000,000 by default, `-1` for none) with
    /// `variant size (N) is bigger than --max_sv_size (M)`: with a limit of 5,
    /// `21:3:6:` (end 9, start 4) is kept and `21:3:7:` refused, by letters as by
    /// length, after the reference check, so a mismatch still names the bases;
    /// without a limit a 30-base deletion is read whole.
    #[test]
    fn oversize_span_is_refused_as_validate_vf_refuses_it() {
        let limited = ParserContext {
            max_sv_size: Some(5),
            ..context()
        };
        let read = |spdi: &str| {
            parser_in(&limited, &format!("{spdi}\n"))
                .next_variant()
                .unwrap()
        };
        assert_eq!(read("21:3:6:").unwrap().allele_string, "ATACGT/-");
        for spdi in ["21:3:7:", "21:3:ATACGTT:", "21:3:7:AAAAAAA"] {
            assert_eq!(
                read(spdi).unwrap_err().to_string(),
                "variant size (6) is bigger than --max_sv_size (5)",
                "{spdi}"
            );
        }
        assert!(read("21:3:GGGGGGG:").unwrap_err().to_string().starts_with(
            "Unable to parse SPDI notation '21:3:GGGGGGG:': Reference allele extracted"
        ));
        assert_eq!(one("21:3:30:").allele_string.len(), 32);
    }

    /// A span past the end of a sequence the reference holds is refused naming
    /// the span and the sequence's length, where a sequence the reference lacks
    /// keeps `Sequence name <name> not valid` (`:2536`): on the 38-base `21`,
    /// `21:38::A` inserts after the last base and is read (Ensembl VEP 116.2
    /// annotates `NC_000021.9:46709983::A` at `21:46709983-46709984` on the
    /// 46,709,983-base chromosome), `21:39::A` falls between bases 39 and 40,
    /// `21:38:1:` deletes base 39, `21:35:5:` runs to base 40 and `21:37:NN:A`
    /// names bases 38 and 39. The name as written is kept in the text.
    #[test]
    fn span_past_the_end_is_refused_naming_the_length() {
        let v = one("21:38::A");
        assert_eq!((v.start, v.end), (39, 38));
        assert_eq!(v.allele_string, "-/A");
        let past_end = |spdi: &str, span: &str| {
            format!(
                "Unable to parse SPDI notation '{spdi}': Span {span} lies past the end of \
                 sequence 21 (38 bases)"
            )
        };
        assert_eq!(refused("21:39::A"), past_end("21:39::A", "39-40"));
        assert_eq!(refused("21:38:1:"), past_end("21:38:1:", "39-39"));
        assert_eq!(refused("21:35:5:"), past_end("21:35:5:", "36-40"));
        assert_eq!(refused("21:37:NN:A"), past_end("21:37:NN:A", "38-39"));
        assert_eq!(
            refused("NC_000021.9:40:1:"),
            "Unable to parse SPDI notation 'NC_000021.9:40:1:': Span 41-41 lies past the end \
             of sequence NC_000021.9 (38 bases)"
        );
        assert_eq!(
            refused("NC_000021:40:1:"),
            "Unable to parse SPDI notation 'NC_000021:40:1:': Sequence name NC_000021 not valid"
        );
    }

    /// A deletion of length `0` (or `00`) with nothing inserted is refused
    /// naming the empty span, start to the position, never the sequence: Ensembl
    /// VEP reads it as a deletion (`:2552`), ends it at the position
    /// (`:2563`), takes the empty slice's sequence as its reference allele
    /// (`:2564`) and keeps the feature `/-` (`:2566`), which writes no
    /// consequence row. `0` for both parts stays `SPDI notation not supported.`
    /// (`:2540-2542`), and a length that matches no inserted bases keeps the
    /// mismatch text (`:2586`).
    #[test]
    fn a_zero_length_deletion_is_refused_naming_the_empty_span() {
        let empty_span = |spdi: &str, length: &str, span: &str| {
            format!(
                "Unable to parse SPDI notation '{spdi}': Deleted sequence length ({length}) \
                 leaves nothing to delete: span {span} is empty"
            )
        };
        assert_eq!(
            refused("21:10538884:0:"),
            empty_span("21:10538884:0:", "0", "10538885-10538884")
        );
        assert_eq!(refused("21:3:00:"), empty_span("21:3:00:", "00", "4-3"));
        assert_eq!(refused("21:3:00:0"), empty_span("21:3:00:0", "00", "4-3"));
        assert_eq!(
            refused("21:3:0:0"),
            "Unable to parse SPDI notation '21:3:0:0': Could not parse the SPDI notation \
             21:3:0:0. SPDI notation not supported."
        );
        assert_eq!(
            refused("21:3:00:G"),
            "Unable to parse SPDI notation '21:3:00:G': Could not parse the SPDI notation \
             21:3:00:G. Deleted sequence length (00) does not match inserted sequence length (1)."
        );
    }

    /// Without a reference every record fails as a sequence the reference lacks:
    /// the format needs the run's FASTA.
    #[test]
    fn no_reference_fails_every_record() {
        let context = ParserContext {
            fasta: None,
            ..context()
        };
        let mut parser = parser_in(&context, "21:3:A:C\n");
        assert_eq!(
            parser.next_variant().unwrap().unwrap_err().to_string(),
            "Unable to parse SPDI notation '21:3:A:C': Sequence name 21 not valid"
        );
    }
}
