// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The breakend classes of the 116.2 reference: own-end rows lost to the mate chromosome,
//! mate rows lost to a chromosome spelling, and the three shapes that are counted.

use std::collections::{BTreeSet, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use regex::Regex;

use crate::cache::CacheAuthority;
use crate::failure::Failure;
use crate::input::{for_each_record, KeyPrep};
use crate::locus::{
    chromosome_of_location, location_span, strip_chr, BREAKEND_FEATURE_DISTANCE,
    INTERGENIC_TUPLE_TAIL,
};
use crate::mate_context::decimal_value;
use crate::perl::info_field;
use crate::record::{perl_skip_route, record_location, span_location};
use crate::rules::breakend::{
    is_own_breakend_allele, reference_mate_chromosome, BARE_BREAKPOINT_ALLELE,
    BND_MATE_32_BIT_LIMIT,
};
use crate::rules::gvcf::perl_split;
use crate::strings::{Sides, Tuple, TupleSet};
use crate::text::os_failure;
use crate::threads::{crew_indexed, crew_map};
use crate::vcf::VcfVariant;

/// The Feature of a row that names no feature.
const NO_FEATURE: &str = "-";

static BRACKET_MATE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]").expect("regex"));

/// The mate a bracket-notation breakend allele names, as its chromosome without a leading
/// `chr` and the digits of its position, of any script and any length; `None` for any
/// other allele.
fn bracket_mate(allele: &str) -> Option<(&str, &str)> {
    let m = BRACKET_MATE.captures(allele)?;
    let (chrom, position) = (m.get(1)?.as_str(), m.get(2)?.as_str());
    Some((strip_chr(chrom), position))
}

/// The tuples of a set by Location: each run of the sorted set at one Location.
fn by_location(set: &TupleSet) -> impl Iterator<Item = &[Tuple]> {
    set.as_slice().chunk_by(|a, b| a.location() == b.location())
}

/// The scored side's own-end rows of a bracket breakend whose own-end rows the reference
/// lost to its mate chromosome.
///
/// A record's slice is named by the first transcript that reaches it
/// (AnnotationType/Transcript.pm 116.2:122, `$vf->{slice} ||= $tr->{slice}`). A buffer
/// holds the records of one chromosome (InputBuffer.pm 116.2:191-207, 223-234), its regions
/// load in the order its records first ask for them, each record's own region before its
/// mates' and a region listed once at its first request (AnnotationSource.pm
/// 116.2:235-258), and the transcripts are visited in that order. When an earlier record of
/// the buffer asked for a region holding a transcript within reach of the record's mate
/// before the record's own region was asked for, that mate-chromosome transcript names the
/// slice, every own-chromosome transcript then fails the slice test (:125-143) and the
/// reference writes nothing under the own-end allele while keeping the mate rows. The
/// result depends on the order of the input; the record alone, with its mate, or under
/// `--buffer_size 1` keeps the rows.
///
/// At a Location whose bracket allele names a mate on another chromosome, where the
/// reference wrote a Transcript tuple under a bracket allele and no tuple under any own-end
/// allele: every scored-only Transcript tuple under an own-end allele whose transcript the
/// cache files on the Location's chromosome within [`BREAKEND_FEATURE_DISTANCE`] of the
/// breakpoint, the reference's own admission rule (StructuralVariationOverlap.pm
/// 116:144-160). A transcript farther away, or on another chromosome, stays charged.
/// Without a cache nothing is set aside.
pub fn filter_bnd_own_rows_lost(sides: &Sides<'_>, cache: &CacheAuthority) -> TupleSet {
    if cache.by_chr.is_empty() || cache.spans.is_empty() {
        return TupleSet::new();
    }
    let allele = |tuple: &Tuple| sides.strings.text(tuple.allele());
    let feature = |tuple: &Tuple| sides.strings.text(tuple.feature());
    let mut excluded: Vec<Tuple> = Vec::new();
    for rust_here in by_location(&sides.rust_only) {
        let location = rust_here[0].location();
        let loc = sides.strings.text(location);
        let (Some(own), Some((_, pos, _))) = (chromosome_of_location(loc), location_span(loc))
        else {
            continue;
        };
        let perl_here = sides.perl.with_prefix(&[location]);
        let names_another_chromosome = perl_here
            .iter()
            .chain(rust_here)
            .filter_map(|tuple| bracket_mate(allele(tuple)))
            .any(|(mate, _)| mate != own);
        let reference_kept_an_own_end = perl_here
            .iter()
            .any(|tuple| is_own_breakend_allele(allele(tuple)));
        let reference_wrote_a_mate_transcript = perl_here
            .iter()
            .any(|tuple| bracket_mate(allele(tuple)).is_some() && feature(tuple) != NO_FEATURE);
        if !names_another_chromosome
            || reference_kept_an_own_end
            || !reference_wrote_a_mate_transcript
        {
            continue;
        }
        let Some(local) = cache.transcripts(own) else {
            continue;
        };
        for tuple in rust_here {
            let transcript = feature(tuple);
            if !is_own_breakend_allele(allele(tuple))
                || transcript == NO_FEATURE
                || !local.contains(transcript)
            {
                continue;
            }
            let Some((start, end)) = cache.span(own, transcript) else {
                continue;
            };
            if start.saturating_sub(BREAKEND_FEATURE_DISTANCE) <= pos
                && pos <= end.saturating_add(BREAKEND_FEATURE_DISTANCE)
            {
                excluded.push(*tuple);
            }
        }
    }
    TupleSet::from_unsorted(excluded)
}

/// The scored side's mate rows of a bracket breakend whose mate chromosome is spelled in a
/// way the reference does not resolve.
///
/// A mate whose chromosome carries a prefix the mate-allele test does not resolve, `Chr21`
/// for `21`, reaches the record's transcripts through the record's own coordinates and
/// passes the slice test; the mate allele is then added only when `_close_to_feature`
/// passes (StructuralVariationOverlap.pm 116:86-103, 144-160), whose name test
/// (Utils/VariationEffect.pm 116:111-124) matches a synonym table key verbatim or the whole
/// name case-folded, and `Chr21` is neither, so the reference writes the own-end rows
/// alone. A mate named by an accession takes the other route: on another chromosome the
/// overlap test resolves it (InputBuffer.pm 116.2:348-357, `get_source_chr_name`) and the
/// slice test then compares the raw name, `chr` stripped, with the cache's
/// (AnnotationType/Transcript.pm 116.2:134), so the reference writes no mate row. The
/// scored engine resolves either name through the synonym table and writes the mate rows.
///
/// At a Location where the reference wrote at least one tuple, for a bracket allele the
/// reference wrote no tuple under and whose mate name, read as the reference reads it
/// ([`reference_mate_chromosome`]), is no chromosome of the cache: every scored-only
/// Transcript tuple under that allele whose transcript the cache files on the chromosome
/// the scored engine resolved the name to (`N[Chr21:30425700[` on chromosome 21) or, for a
/// name the cache files no chromosome under (`NC_000022.11`), whose transcript is not on
/// the Location's own chromosome. Without a cache nothing is set aside.
pub fn filter_bnd_synonym_mate_rows(sides: &Sides<'_>, cache: &CacheAuthority) -> TupleSet {
    if cache.by_chr.is_empty() {
        return TupleSet::new();
    }
    let mut excluded: Vec<Tuple> = Vec::new();
    for rust_here in by_location(&sides.rust_only) {
        let location = rust_here[0].location();
        let Some(own) = chromosome_of_location(sides.strings.text(location)) else {
            continue;
        };
        if sides.perl.with_prefix(&[location]).is_empty() {
            continue;
        }
        let local = cache.transcripts(own);
        for tuple in rust_here {
            let allele = sides.strings.text(tuple.allele());
            let transcript = sides.strings.text(tuple.feature());
            let Some((mate, _)) = bracket_mate(allele) else {
                continue;
            };
            let reference_resolves_the_name = reference_mate_chromosome(allele)
                .is_some_and(|name| cache.by_chr.contains_key(name));
            let reference_wrote_the_allele = !sides
                .perl
                .with_prefix(&[location, tuple.allele()])
                .is_empty();
            if reference_resolves_the_name || reference_wrote_the_allele || transcript == NO_FEATURE
            {
                continue;
            }
            let lost = match cache.transcripts(mate) {
                Some(resolved) => resolved.contains(transcript),
                None => !local.is_some_and(|local| local.contains(transcript)),
            };
            if lost {
                excluded.push(*tuple);
            }
        }
    }
    TupleSet::from_unsorted(excluded)
}

/// The records whose only reference row is intergenic while the scored side wrote
/// Transcript rows and no intergenic row: the shape of a slice named after a FASTA whose
/// sequence names differ from the cache's (Parser.pm 116.2:749 through BaseVEP.pm
/// `get_slice`, every transcript failing AnnotationType/Transcript.pm:141). The condition
/// is the run's configuration, so the records are counted and nothing is set aside.
pub fn count_fasta_named_slice_intergenic(sides: &Sides<'_>) -> u64 {
    let [Some(feature), Some(feature_type), Some(consequence)] =
        INTERGENIC_TUPLE_TAIL.map(|text| sides.strings.get(text))
    else {
        return 0;
    };
    let mut records = 0;
    for perl_here in sides
        .perl
        .as_slice()
        .chunk_by(|a, b| a.variant() == b.variant())
    {
        let &[only] = perl_here else {
            continue;
        };
        if only.0[2..] != [feature, feature_type, consequence] {
            continue;
        }
        let rust_here = sides.rust.with_prefix(&[only.location(), only.allele()]);
        if !rust_here.is_empty() && rust_here.iter().all(|tuple| tuple.feature() != feature) {
            records += 1;
        }
    }
    records
}

/// Whether a text of decimal digits of any script names a coordinate at or past
/// [`BND_MATE_32_BIT_LIMIT`], at any length.
fn reaches_the_bound(digits: &str) -> bool {
    decimal_value(digits).is_some_and(|position| position >= BND_MATE_32_BIT_LIMIT)
}

/// The input records with a coordinate at or past the buffer's 32-bit bound, and the
/// reference's tuples at those records' Locations: `(records, reference_tuples)`.
///
/// The coordinate is a breakend mate in a bracket ALT or through CHR2 with END2 (else END),
/// or, for a record the reference keeps under `max_sv_size` ([`perl_skip_route`]), the
/// record's own start or the end its span derives (`21 30425650 x N <DEL>
/// SVTYPE=DEL;END=2147483647` under `--max_sv_size -1`, `21 2147483648 x N N[21:30425700[
/// SVTYPE=BND`). The reference writes every record of the buffer holding such a coordinate
/// as `intergenic_variant`, so the count names the shape and sets no row aside: the rows
/// reach the buffer's other records, and the outputs do not say which records shared the
/// buffer.
pub fn count_bnd_mate_over_32_bit(
    records: &[VcfVariant],
    sides: &Sides<'_>,
    max_sv_size: i128,
) -> (u64, u64) {
    // Whether a record is counted depends on that record alone; a counted record gives
    // its Location.
    let reaching = crew_map(records, |v| {
        let bracket_mate_reaches = perl_split(&v.alt, ',')
            .into_iter()
            .filter_map(bracket_mate)
            .any(|(_, position)| reaches_the_bound(position));
        let named_mate_reaches = info_field(&v.info, "CHR2").is_some()
            && info_field(&v.info, "END2")
                .filter(|end2| !end2.is_empty())
                .or_else(|| info_field(&v.info, "END"))
                .is_some_and(reaches_the_bound);
        let (chrom, start, end) = record_location(v);
        let own_span_reaches = perl_skip_route(v, max_sv_size).is_none()
            && (start >= BND_MATE_32_BIT_LIMIT || end >= BND_MATE_32_BIT_LIMIT);
        (bracket_mate_reaches || named_mate_reaches || own_span_reaches)
            .then(|| span_location(&chrom, start, end))
    });
    let counted = reaching.iter().flatten().count() as u64;
    let locations: BTreeSet<String> = reaching.into_iter().flatten().collect();
    let tuples: usize = locations
        .iter()
        .filter_map(|loc| sides.strings.get(loc))
        .map(|location| sides.perl.with_prefix(&[location]).len())
        .sum();
    (counted, tuples as u64)
}

/// The reference's duplicate rows under `chromosome_breakpoint` in one raw output; a
/// missing file counts none.
///
/// A breakpoint record that names no mate is held as `N/<BND>` and, with no breakend
/// parsed, gets one allele per slash-separated part, the reference base `N` included
/// (StructuralVariationOverlap.pm:64-77), so every row is written twice, the second for the
/// reference allele, under the same key. The second and later rows of a key are counted:
/// the rows the reader of an output drops as duplicates of a tuple. A row that reader
/// skips (a `#` line, fewer than seven columns, an empty Location, Allele or Feature)
/// counts nothing.
pub fn count_bare_breakpoint_reference_rows(path: &Path) -> Result<u64, Failure> {
    if !path.exists() {
        return Ok(0);
    }
    let unreadable = |error: std::io::Error| os_failure(&error, Some(path));
    let file = File::open(path).map_err(unreadable)?;
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut prep = KeyPrep::default();
    let mut key = Vec::new();
    let mut rows = 0;
    for_each_record(file, |record| {
        key.clear();
        if !prep.prepare(record, false, None, &mut key) {
            return Ok(());
        }
        if key.split(|byte| *byte == b'\t').nth(1) != Some(BARE_BREAKPOINT_ALLELE.as_bytes()) {
            return Ok(());
        }
        if seen.contains(&key) {
            rows += 1;
        } else {
            seen.insert(key.clone());
        }
        Ok(())
    })
    .map_err(unreadable)?;
    Ok(rows)
}

/// The reference-allele rows of bare breakpoint records in the reference's raw outputs,
/// summed over the files; a missing file counts zero.
///
/// The raw rows are read because the reader of an output has already dropped each
/// duplicate from the tuple set: the rows counted here are outside both the raw and the
/// adjusted figure, and the count is what the report states as set aside.
pub fn filter_bnd_bare_reference_allele_rows(perl_paths: &[PathBuf]) -> Result<u64, Failure> {
    let rows = crew_indexed(perl_paths.len(), |at| {
        count_bare_breakpoint_reference_rows(&perl_paths[at])
    })?;
    Ok(rows.into_iter().sum())
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::masks::{apply_masks, MaskContext, MaskResult};
    use crate::record::testing::{intergenic, rec_id, rows, t};
    use crate::strings::OwnedSides;
    use crate::TempDir;

    fn sides_of(perl: &[[String; 5]], rust: &[[String; 5]]) -> OwnedSides {
        OwnedSides::from_rows(&rows(perl), &rows(rust))
    }

    /// The chromosome view of a cache: each chromosome with its transcript ids.
    fn tx_by_chr(filed: &[(&str, &[&str])]) -> CacheAuthority {
        let mut cache = CacheAuthority::default();
        for (chrom, ids) in filed {
            cache.by_chr.insert(
                chrom.to_string(),
                ids.iter().map(|id| id.to_string()).collect(),
            );
        }
        cache
    }

    /// The tuples of a set as rows of texts, in order.
    fn texts(owned: &OwnedSides, set: &TupleSet) -> Vec<[String; 5]> {
        set.iter()
            .map(|tuple| owned.strings.resolve(tuple).map(str::to_string))
            .collect()
    }

    fn sorted(mut expected: Vec<[String; 5]>) -> Vec<[String; 5]> {
        expected.sort_unstable();
        expected
    }

    /// The 116.2 masks on two tuple sets with no cache, as the reference ran under `cap`.
    fn masks_116(
        owned: &OwnedSides,
        records: &[VcfVariant],
        cap: i128,
        perl_paths: &[PathBuf],
    ) -> MaskResult {
        let context = MaskContext {
            reference_release: "116.2",
            cache: &CacheAuthority::default(),
            scored_engine: "vep-rs",
            input_variants: Some(records),
            max_sv_size: cap,
            perl_paths,
            rust_paths: &[],
            input_by_output: &BTreeMap::new(),
            coding_ends: None,
        };
        apply_masks(&owned.sides(), &context).expect("the masks run")
    }

    mod sv_filter_set_tests {
        use super::*;

        /// Ensembl VEP 116.2 writes a bare `<BND>`'s rows twice, ALLELE_NUM 1 for the
        /// reference allele N of `N/<BND>` and 2 for the ALT
        /// (StructuralVariationOverlap.pm:64-77). The mask reads the raw rows and returns
        /// the reference's second row per feature, the rows the reader of an output dropped
        /// from the tuple set, so the set-aside count is 2 while the tuple set carries 3
        /// `chromosome_breakpoint` tuples; a bracket record's own and mate rows, which
        /// differ in Allele, count nothing, and a missing file counts zero.
        #[test]
        fn test_the_bare_breakpoint_reference_allele_rows_are_set_aside_and_counted() {
            let lines = [
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence",
                "bare\t21:30425651\tchromosome_breakpoint\tG1\tENST1\tTranscript\tcoding_sequence_variant",
                "bare\t21:30425651\tchromosome_breakpoint\tG1\tENST1\tTranscript\tcoding_sequence_variant",
                "bare\t21:30425651\tchromosome_breakpoint\tG2\tENST2\tTranscript\tupstream_gene_variant",
                "bare\t21:30425651\tchromosome_breakpoint\tG2\tENST2\tTranscript\tupstream_gene_variant",
                "paired\t21:30425651\tN.\tG1\tENST1\tTranscript\tfeature_truncation,coding_sequence_variant",
                "paired\t21:30425651\tN[21:30425700[\tG1\tENST1\tTranscript\tfeature_truncation,coding_sequence_variant",
                "single\t21:30425651\tchromosome_breakpoint\tG3\tENST3\tTranscript\tintron_variant",
            ];
            let dir = TempDir::new("vep_sv_compare_bare_").expect("a scratch directory");
            let path = dir.path.join("perl.txt");
            std::fs::write(&path, lines.join("\n") + "\n").expect("an output");
            let one = [path.clone()];
            assert_eq!(filter_bnd_bare_reference_allele_rows(&one), Ok(2));
            let tuples: Vec<[String; 5]> = lines[1..]
                .iter()
                .map(|line| {
                    let column: Vec<&str> = line.split('\t').collect();
                    t(column[1], column[2], column[4], column[6])
                })
                .collect();
            let owned = sides_of(&tuples, &tuples);
            let bare = owned
                .perl
                .iter()
                .filter(|tuple| owned.strings.text(tuple.allele()) == "chromosome_breakpoint")
                .count();
            assert_eq!((owned.perl.len(), bare), (5, 3));
            let absent = [dir.path.join("absent.txt")];
            assert_eq!(filter_bnd_bare_reference_allele_rows(&absent), Ok(0));
            let twice = [path.clone(), path];
            assert_eq!(filter_bnd_bare_reference_allele_rows(&twice), Ok(4));
            let masks = masks_116(&owned, &[], crate::masks::MAX_SV_SIZE, &one);
            assert_eq!(masks.bare_reference_rows, 2);
            assert_eq!(
                masks.all_perl(),
                TupleSet::new(),
                "rows, not tuples: nothing leaves the tuple arithmetic"
            );
        }
    }

    /// A bracket breakend whose mate is on another chromosome loses its own-end rows when
    /// a mate-chromosome transcript, loaded for an earlier record of the same buffer, names
    /// its slice first (AnnotationType/Transcript.pm 116.2:122-143 with InputBuffer.pm
    /// 116.2:191-207, 223-234 and AnnotationSource.pm 116.2:235-258).
    mod bnd_own_rows_lost_tests {
        use super::*;

        const LOC: &str = "21:26037668";
        const MATE: &str = "N[22:26037668[";

        /// ENST21_NEAR lies 2,668 bases from the breakpoint, ENST21_FAR 27,668 bases away.
        fn cache() -> CacheAuthority {
            let mut cache = tx_by_chr(&[
                ("21", &["ENST21_NEAR", "ENST21_FAR"]),
                ("22", &["ENST22_M"]),
            ]);
            for (chrom, id, span) in [
                ("21", "ENST21_NEAR", (26030000, 26035000)),
                ("21", "ENST21_FAR", (26000000, 26010000)),
                ("22", "ENST22_M", (26036000, 26040000)),
            ] {
                cache
                    .spans
                    .entry(chrom.to_string())
                    .or_default()
                    .insert(id.to_string(), span);
            }
            cache
        }

        fn mate_row() -> [String; 5] {
            t(LOC, MATE, "ENST22_M", "intron_variant")
        }

        fn rust(own_features: &[&str]) -> Vec<[String; 5]> {
            let mut tuples = vec![mate_row()];
            for feature in own_features {
                tuples.push(t(LOC, "N.", feature, "upstream_gene_variant"));
            }
            tuples
        }

        fn lost(
            perl: &[[String; 5]],
            rust: &[[String; 5]],
            cache: &CacheAuthority,
        ) -> Vec<[String; 5]> {
            let owned = sides_of(perl, rust);
            texts(&owned, &filter_bnd_own_rows_lost(&owned.sides(), cache))
        }

        #[test]
        fn test_own_end_rows_within_the_admission_window_are_excluded() {
            let got = lost(&[mate_row()], &rust(&["ENST21_NEAR"]), &cache());
            assert_eq!(got, [t(LOC, "N.", "ENST21_NEAR", "upstream_gene_variant")]);
        }

        /// The scored engine writes the local end of a bracket breakend as
        /// `chromosome_breakpoint`.
        #[test]
        fn test_vep_rs_own_allele_string_may_be_the_class_term() {
            let rust = [
                mate_row(),
                t(
                    LOC,
                    "chromosome_breakpoint",
                    "ENST21_NEAR",
                    "upstream_gene_variant",
                ),
            ];
            assert_eq!(lost(&[mate_row()], &rust, &cache()).len(), 1);
        }

        /// The reference's own rule admits a breakend within 5,000 bases of the transcript
        /// (StructuralVariationOverlap.pm 116:144-160); a farther row is the scored
        /// engine's alone.
        #[test]
        fn test_a_transcript_beyond_the_window_stays_charged() {
            let got = lost(
                &[mate_row()],
                &rust(&["ENST21_NEAR", "ENST21_FAR"]),
                &cache(),
            );
            let features: Vec<&str> = got.iter().map(|tuple| tuple[2].as_str()).collect();
            assert_eq!(features, ["ENST21_NEAR"]);
        }

        #[test]
        fn test_the_reference_keeping_any_own_end_row_is_not_the_shape() {
            let perl = [mate_row(), intergenic(LOC, "N.")];
            assert_eq!(lost(&perl, &rust(&["ENST21_NEAR"]), &cache()), [[""; 5]; 0]);
        }

        #[test]
        fn test_a_same_chromosome_mate_is_not_the_shape() {
            let (loc, mate) = ("21:24808769-35740566", "N[21:35740567[");
            let perl = [t(loc, mate, "ENST21_NEAR", "intron_variant")];
            let rust = [
                perl[0].clone(),
                t(loc, "N.", "ENST21_NEAR", "intron_variant"),
            ];
            assert_eq!(lost(&perl, &rust, &cache()), [[""; 5]; 0]);
        }

        /// No reference Transcript row under the bracket allele means the record was not
        /// reached through its mate; that is the skipped-record class or a real divergence.
        #[test]
        fn test_the_reference_dropping_the_mate_rows_too_is_not_the_shape() {
            let perl = [intergenic(LOC, MATE)];
            assert_eq!(lost(&perl, &rust(&["ENST21_NEAR"]), &cache()), [[""; 5]; 0]);
        }

        #[test]
        fn test_no_cache_excludes_nothing() {
            let got = lost(
                &[mate_row()],
                &rust(&["ENST21_NEAR"]),
                &CacheAuthority::default(),
            );
            assert_eq!(got, [[""; 5]; 0]);
        }

        #[test]
        fn a_mate_position_of_any_length_and_script_names_its_chromosome() {
            for position in ["9".repeat(60), "\u{0662}\u{0666}".to_string()] {
                let mate = format!("N[22:{position}[");
                let mate_row = t(LOC, &mate, "ENST22_M", "intron_variant");
                let rust = [
                    mate_row.clone(),
                    t(LOC, "N.", "ENST21_NEAR", "upstream_gene_variant"),
                ];
                assert_eq!(lost(&[mate_row], &rust, &cache()).len(), 1, "{position}");
            }
        }
    }

    /// A mate whose chromosome carries a prefix the mate-allele test does not resolve
    /// (`Chr21` for `21`) loses its mate rows in the reference alone
    /// (StructuralVariationOverlap.pm 116:86-103, 144-160; Utils/VariationEffect.pm
    /// 116:111-124); a mate named by an accession on another chromosome passes the overlap
    /// test (InputBuffer.pm 116.2:348-357) and fails the slice test
    /// (AnnotationType/Transcript.pm 116.2:134), while the scored engine resolves it
    /// through the synonym table and writes the mate rows, which are set aside.
    mod bnd_synonym_mate_rows_tests {
        use super::*;

        const LOC: &str = "21:26037668";
        const SYN: &str = "N[NC_000022.11:26037668[";

        fn cache() -> CacheAuthority {
            tx_by_chr(&[("21", &["ENST21_OWN"]), ("22", &["ENST22_M"])])
        }

        fn own_row() -> [String; 5] {
            t(LOC, "N.", "ENST21_OWN", "upstream_gene_variant")
        }

        fn lost(
            perl: &[[String; 5]],
            rust: &[[String; 5]],
            cache: &CacheAuthority,
        ) -> Vec<[String; 5]> {
            let owned = sides_of(perl, rust);
            texts(&owned, &filter_bnd_synonym_mate_rows(&owned.sides(), cache))
        }

        #[test]
        fn test_mate_side_rows_under_a_synonym_named_mate_are_excluded() {
            let mate_row = t(LOC, SYN, "ENST22_M", "intron_variant");
            let got = lost(&[own_row()], &[own_row(), mate_row.clone()], &cache());
            assert_eq!(got, [mate_row]);
        }

        #[test]
        fn test_a_cache_named_mate_is_not_the_shape() {
            let rust = [
                own_row(),
                t(LOC, "N[22:26037668[", "ENST22_M", "intron_variant"),
            ];
            assert_eq!(lost(&[own_row()], &rust, &cache()), [[""; 5]; 0]);
        }

        #[test]
        fn test_an_own_chromosome_transcript_under_the_mate_allele_stays_charged() {
            let rust = [own_row(), t(LOC, SYN, "ENST21_OWN", "intron_variant")];
            assert_eq!(lost(&[own_row()], &rust, &cache()), [[""; 5]; 0]);
        }

        #[test]
        fn test_the_reference_writing_the_mate_allele_is_not_the_shape() {
            let perl = [t(LOC, SYN, "ENST22_M", "intron_variant")];
            let rust = [
                perl[0].clone(),
                t(LOC, SYN, "ENST22_X", "downstream_gene_variant"),
            ];
            assert_eq!(lost(&perl, &rust, &cache()), [[""; 5]; 0]);
        }

        #[test]
        fn test_a_record_the_reference_wrote_nothing_for_is_the_skipped_class_not_this_one() {
            let rust = [t(LOC, SYN, "ENST22_M", "intron_variant")];
            assert_eq!(lost(&[], &rust, &cache()), [[""; 5]; 0]);
        }

        #[test]
        fn test_no_cache_excludes_nothing() {
            let rust = [own_row(), t(LOC, SYN, "ENST22_M", "intron_variant")];
            let got = lost(&[own_row()], &rust, &CacheAuthority::default());
            assert_eq!(got, [[""; 5]; 0]);
        }
    }

    mod fasta_named_slice_intergenic_tests {
        use super::*;

        #[test]
        fn test_counts_records_whose_only_reference_row_is_intergenic_against_vep_rs_transcripts() {
            let perl = [
                intergenic("21:100", "A"),
                intergenic("21:200", "T"),
                t("21:300", "G", "ENST_3", "intron_variant"),
            ];
            let rust = [
                t("21:100", "A", "ENST_1", "missense_variant"),
                t("21:100", "A", "ENST_1b", "intron_variant"),
                intergenic("21:200", "T"),
                t("21:300", "G", "ENST_3", "intron_variant"),
            ];
            let owned = sides_of(&perl, &rust);
            assert_eq!(count_fasta_named_slice_intergenic(&owned.sides()), 1);
        }

        /// Ensembl VEP 116.2 turns the buffer intergenic for a mate at 2147478648 named by
        /// a bracket ALT or by CHR2 with END2 (else END); 2147478647 and a mate on another
        /// route are not the shape. The reference's tuples at the record's Location are
        /// counted with it.
        #[test]
        fn test_a_mate_past_the_buffers_32_bit_bound_is_counted_from_the_input_record() {
            let record = |vid: &str, alt: &str, info: &str| VcfVariant {
                variant_type: "BND".to_string(),
                ..rec_id(30425650, "N", alt, info, vid)
            };
            let records = [
                record("bracket", "N[21:2147478648[", "SVTYPE=BND"),
                record("chr2_end2", "<BND>", "SVTYPE=BND;CHR2=21;END2=2147478648"),
                record("chr2_end", "<BND>", "SVTYPE=BND;CHR2=21;END=2147478648"),
                record("below", "<BND>", "SVTYPE=BND;CHR2=21;END2=2147478647"),
                record("pos2", "<BND>", "SVTYPE=BND;CHR2=21;POS2=2147478648"),
            ];
            let perl = [
                intergenic("21:30425651", "N."),
                intergenic("21:30425651", "N[21:2147478648["),
                t("21:40000001", "N.", "ENST_1", "intron_variant"),
            ];
            let owned = sides_of(&perl, &[]);
            let cap = crate::masks::MAX_SV_SIZE;
            assert_eq!(
                count_bnd_mate_over_32_bit(&records, &owned.sides(), cap),
                (3, 2)
            );
            assert_eq!(
                count_bnd_mate_over_32_bit(&records[3..], &owned.sides(), cap),
                (0, 0)
            );
        }

        #[test]
        fn test_a_vep_rs_intergenic_row_beside_its_transcripts_is_not_the_shape() {
            let perl = [intergenic("21:100", "A")];
            let rust = [
                intergenic("21:100", "A"),
                t("21:100", "A", "ENST_1", "missense_variant"),
            ];
            let owned = sides_of(&perl, &rust);
            assert_eq!(count_fasta_named_slice_intergenic(&owned.sides()), 0);
        }
    }

    mod synonym_mate_prefix_tests {
        use super::*;

        /// `21 30425650 x N N[Chr21:30425700[ . . SVTYPE=BND`: Ensembl VEP 116.2 writes the
        /// own end only, on ENST00000334068 and ENST00000390690; the scored engine adds the
        /// mate rows `21:30425651 N[Chr21:30425700[` on both transcripts, which are set
        /// aside. With `N[chr21:30425700[` both engines write the four rows and nothing is
        /// set aside, and a `N[NC_000022.11:...[` mate keeps its own-chromosome rule.
        #[test]
        fn test_a_mate_prefix_the_reference_does_not_resolve_loses_the_mate_rows_on_the_own_chromosome(
        ) {
            let cache = tx_by_chr(&[
                ("21", &["ENST00000334068", "ENST00000390690"]),
                ("22", &["ENST22"]),
            ]);
            let on_both = |allele: &str| {
                vec![
                    t(
                        "21:30425651",
                        allele,
                        "ENST00000334068",
                        "upstream_gene_variant",
                    ),
                    t(
                        "21:30425651",
                        allele,
                        "ENST00000390690",
                        "feature_truncation,coding_sequence_variant",
                    ),
                ]
            };
            let lost = |perl: &[[String; 5]], rust: &[[String; 5]]| {
                let owned = sides_of(perl, rust);
                texts(
                    &owned,
                    &filter_bnd_synonym_mate_rows(&owned.sides(), &cache),
                )
            };
            let own = on_both("N.");
            let mate = on_both("N[Chr21:30425700[");
            assert_eq!(
                reference_mate_chromosome("N[Chr21:30425700["),
                Some("Chr21")
            );
            assert_eq!(reference_mate_chromosome("N[chr21:30425700["), Some("21"));
            assert_eq!(
                lost(&own, &[own.clone(), mate.clone()].concat()),
                sorted(mate)
            );
            let lower = [own.clone(), on_both("N[chr21:30425700[")].concat();
            assert_eq!(lost(&lower, &lower), [[""; 5]; 0]);
            let synonym_own = t(
                "21:30425651",
                "N[NC_000022.11:26037668[",
                "ENST00000390690",
                "intron_variant",
            );
            let synonym_other = t(
                "21:30425651",
                "N[NC_000022.11:26037668[",
                "ENST22",
                "intron_variant",
            );
            let mut rust = own.clone();
            rust.extend([synonym_own, synonym_other.clone()]);
            assert_eq!(lost(&own, &rust), [synonym_other]);
        }
    }

    mod buffer_coordinate_overflow_tests {
        use super::*;

        /// Under `--max_sv_size -1`, Ensembl VEP 116.2 writes one row,
        /// `21:30425651-2147483647 deletion - intergenic_variant`, for `21 30425650 x N
        /// <DEL> . . SVTYPE=DEL;END=2147483647`, as for `SVLEN=2147483647`
        /// (`21:30425651-2177909297`) and `SVLEN=-2147483648` (`21:30425651-2177909298`),
        /// where the scored engine writes every transcript from the record to the
        /// chromosome's end; `21 2147483648 x N N[21:30425700[ . . SVTYPE=BND` is
        /// `21:2147483649 N.` and `21:2147483649 N[21:30425700[`, both intergenic, as POS
        /// 2147483647 is at 2147483648. The counter counts the records and the reference's
        /// rows at their Locations; under the default cap the span records are skipped and
        /// not counted, and an end at 2147478647 is below the bound.
        #[test]
        fn test_a_span_end_or_own_position_past_the_buffers_32_bit_bound_is_counted() {
            let records = [
                rec_id(30425650, "N", "<DEL>", "SVTYPE=DEL;END=2147483647", "end"),
                rec_id(
                    30425650,
                    "N",
                    "<DEL>",
                    "SVTYPE=DEL;SVLEN=2147483647",
                    "svlen",
                ),
                rec_id(
                    30425650,
                    "N",
                    "<DEL>",
                    "SVTYPE=DEL;SVLEN=-2147483648",
                    "negative",
                ),
                rec_id(2147483648, "N", "N[21:30425700[", "SVTYPE=BND", "pos"),
                rec_id(
                    2147483647,
                    "N",
                    "N[21:30425700[",
                    "SVTYPE=BND",
                    "pos_minus_one",
                ),
                rec_id(30425650, "N", "<DEL>", "SVTYPE=DEL;END=2147478647", "below"),
            ];
            let perl = [
                intergenic("21:30425651-2147483647", "deletion"),
                intergenic("21:30425651-2177909297", "deletion"),
                intergenic("21:30425651-2177909298", "deletion"),
                intergenic("21:2147483649", "N."),
                intergenic("21:2147483649", "N[21:30425700["),
                intergenic("21:2147483648", "N."),
                intergenic("21:2147483648", "N[21:30425700["),
                intergenic("21:30425651-2147478647", "deletion"),
            ];
            let owned = sides_of(&perl, &[]);
            let cap = crate::masks::MAX_SV_SIZE;
            assert_eq!(
                count_bnd_mate_over_32_bit(&records, &owned.sides(), -1),
                (5, 7)
            );
            assert_eq!(
                count_bnd_mate_over_32_bit(&records, &owned.sides(), cap),
                (2, 4)
            );
            let masks = masks_116(&owned, &records, -1, &[]);
            assert_eq!(
                (masks.mate_overflow_records, masks.mate_overflow_perl_tuples),
                (5, 7)
            );
        }
    }

    #[test]
    fn a_mate_coordinate_is_read_in_any_script_and_at_any_length() {
        let record = |alt: &str, info: &str| rec_id(100, "N", alt, info, "r");
        let counted = |v: VcfVariant| {
            let owned = sides_of(&[intergenic("21:101", "N.")], &[]);
            count_bnd_mate_over_32_bit(&[v], &owned.sides(), crate::masks::MAX_SV_SIZE)
        };
        let long = "9".repeat(60);
        for (alt, info) in [
            (format!("N[21:{long}["), "SVTYPE=BND".to_string()),
            (
                "N[21:\u{ff12}147478648[".to_string(),
                "SVTYPE=BND".to_string(),
            ),
            (
                "N[21:5[,N[21:2147478648[".to_string(),
                "SVTYPE=BND".to_string(),
            ),
            (
                "<BND>".to_string(),
                format!("SVTYPE=BND;CHR2=21;END2={long}"),
            ),
            (
                "<BND>".to_string(),
                "SVTYPE=BND;CHR2=21;END2=;END=2147478648".to_string(),
            ),
            (
                "<BND>".to_string(),
                "SVTYPE=BND;CHR2;END=2147478648".to_string(),
            ),
        ] {
            assert_eq!(counted(record(&alt, &info)), (1, 1), "{alt} {info}");
        }
        for (alt, info) in [
            ("N[21:2147478647[", "SVTYPE=BND"),
            ("<BND>", "SVTYPE=BND;END2=2147478648"),
            ("<BND>", "SVTYPE=BND;CHR2=21;END2=5;END=2147478648"),
            ("<BND>", "SVTYPE=BND;CHR2=21;END2=+2147478648"),
            ("<BND>", "SVTYPE=BND;CHR2=21;END2=\u{b2}147478648"),
        ] {
            assert_eq!(counted(record(alt, info)), (0, 0), "{alt} {info}");
        }
    }

    #[test]
    fn the_rows_of_an_output_are_read_as_its_reader_reads_them() {
        let dir = TempDir::new("vep_sv_compare_bare_").expect("a scratch directory");
        let path = dir.path.join("perl.txt");
        let row = "b\t21:5\tchromosome_breakpoint\tG\tENST1\tTranscript\tintron_variant";
        let text = [
            "## a comment".to_string(),
            row.to_string(),
            // The same key: blanks around a column, a coordinate with leading zeros, the
            // terms in another order and twice.
            "b\t 21:05 \tchromosome_breakpoint\tG\tENST1\tTranscript\tintron_variant, intron_variant\textra".to_string(),
            format!("{row}\r"),
            // No row of the reader: a comment, six columns, an empty Feature.
            format!("#{row}"),
            "b\t21:5\tchromosome_breakpoint\tG\tENST1\tTranscript".to_string(),
            "b\t21:5\tchromosome_breakpoint\tG\t\tTranscript\tintron_variant".to_string(),
            // Another key: another Feature_type, another allele.
            "b\t21:5\tchromosome_breakpoint\tG\tENST1\tRegulatoryFeature\tintron_variant".to_string(),
            "b\t21:5\tN.\tG\tENST1\tTranscript\tintron_variant".to_string(),
            "b\t21:5\tN.\tG\tENST1\tTranscript\tintron_variant".to_string(),
        ]
        .join("\n");
        std::fs::write(&path, text).expect("an output");
        assert_eq!(count_bare_breakpoint_reference_rows(&path), Ok(2));
        let failure = count_bare_breakpoint_reference_rows(&dir.path).expect_err("a directory");
        assert_eq!(failure.kind, crate::failure::FailureKind::Io);
    }
}
