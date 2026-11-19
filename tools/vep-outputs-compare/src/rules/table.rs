// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every divergence class of the 116.2 reference registry as the report states it.

/// One class: the one-sentence definition, whether the adjusted F1 sets it aside, its
/// taxonomy class and the Ensembl VEP source it cites.
pub struct DivergenceClass {
    pub bucket: &'static str,
    pub definition: &'static str,
    pub excludes: bool,
    pub taxonomy_class: &'static str,
    pub perl_citation: &'static str,
}

/// The pair rules of the registry come first, in classification order.
pub const PAIR_RULE_COUNT: usize = 12;

/// Every class of the 116.2 registry as the report states it: pair rules in
/// classification order, then the one-sided rules.
pub const DIVERGENCE_CLASSES_116_2: &[DivergenceClass] = &[
    DivergenceClass {
        bucket: "intronic_overcall_swap",
        definition: "Both sides are pure-intronic sets and vep-rs adds splice_region_variant or splice_polypyrimidine_tract_variant, which the reference's last-written region drops.",
        excludes: false,
        taxonomy_class: "splice_family_swap",
        perl_citation: "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    },
    DivergenceClass {
        bucket: "ppt_overcall_swap",
        definition: "Both sides carry intron_variant and vep-rs adds splice_polypyrimidine_tract_variant or splice_region_variant to the same set.",
        excludes: false,
        taxonomy_class: "splice_family_swap",
        perl_citation: "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    },
    DivergenceClass {
        bucket: "frameshift_splice_cascade_swap",
        definition: "An identical frameshift or in-frame coding backbone on both sides, vep-rs adding splice sub-terms the reference's last-written region drops.",
        excludes: false,
        taxonomy_class: "splice_family_swap",
        perl_citation: "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    },
    DivergenceClass {
        bucket: "splice_lastwrite_swap",
        definition: "vep-rs adds splice terms, and nothing else, to the reference's set: the reference stores the splice verdict by assignment and a later differing region overwrote it.",
        excludes: false,
        taxonomy_class: "splice_family_swap",
        perl_citation: "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    },
    DivergenceClass {
        bucket: "start_cooccurrence_swap",
        definition: "The reference writes start_lost and start_retained_variant together for one allele; vep-rs writes the set without the member the allele kind cannot support.",
        excludes: true,
        taxonomy_class: "start_cooccurrence_swap",
        perl_citation: "VariationEffect.pm:851-963, 1028-1075",
    },
    DivergenceClass {
        bucket: "covered_splice_region_swap",
        definition: "The pairs of the splice family whose only difference is a splice_region_variant the reference's assignment order dropped; the exclusion policy of that family.",
        excludes: true,
        taxonomy_class: "splice_family_swap",
        perl_citation: "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    },
    DivergenceClass {
        bucket: "stop_codon_overlap_genomic_span",
        definition: "A Transcript row of an insertion or deletion with both cDNA positions and a CDS start on which the reference writes coding_sequence_variant and vep-rs stop_retained_variant, both beside splice_region_variant: the reference's genomic stop codon overlap reads the last three genomic positions of the coding region and misses a stop base that ends an exon.",
        excludes: true,
        taxonomy_class: "stop_codon_overlap_genomic_span",
        perl_citation: "VariationEffect.pm:1389-1431, 1616-1649",
    },
    DivergenceClass {
        bucket: "stop_retained_without_stop_codon_test",
        definition: "A Transcript row of an insertion or deletion whose stop codon test the reference does not make (an allele holding a letter that names no base, or a span without a cDNA position at one end or a CDS start) and writes stop_retained_variant for; vep-rs writes the same set without the term, with nothing, coding_sequence_variant or frameshift_variant in its place.",
        excludes: true,
        taxonomy_class: "stop_retained_without_stop_codon_test",
        perl_citation: "VariationEffect.pm:1290-1324, 1488-1539",
    },
    DivergenceClass {
        bucket: "stop_retained_three_base_stop_codon_test",
        definition: "A Transcript row of an insertion or deletion whose stop codon test is made, on which the reference writes stop_retained_variant from the three bases at the old stop position alone; vep-rs writes coding_sequence_variant in its place, or, where the reference writes it alone or beside one of inframe_insertion, inframe_deletion and protein_altering_variant, frameshift_variant with or without stop_gained in the place of those terms.",
        excludes: true,
        taxonomy_class: "stop_retained_three_base_stop_codon_test",
        perl_citation: "VariationEffect.pm:1318, 1488-1539",
    },
    DivergenceClass {
        bucket: "stop_codon_removed_without_stop_codon_test",
        definition: "A Transcript row of an insertion or deletion without a cDNA position at one end or a CDS start, on which the reference writes stop_retained_variant with no stop codon test made and vep-rs stop_lost for the stop codon the record removes.",
        excludes: true,
        taxonomy_class: "stop_codon_removed_without_stop_codon_test",
        perl_citation: "VariationEffect.pm:1234-1288, 1290-1324, 1508",
    },
    DivergenceClass {
        bucket: "stop_codon_overlap_insertion_range",
        definition: "A Transcript row of an insertion with both cDNA positions and a CDS start on which vep-rs writes stop_retained_variant, with or without inframe_insertion, for the stop codon the insertion keeps, and the reference coding_sequence_variant, or frameshift_variant with or without stop_gained: the reference's genomic stop codon overlap moves an insertion's start by the inserted length. Under an identifier a row is read as an insertion's by a Location of two neighbouring positions and an Allele of another length than two, so two bases replaced by a sequence of another length are set aside with the insertions.",
        excludes: true,
        taxonomy_class: "stop_codon_overlap_insertion_range",
        perl_citation: "VariationEffect.pm:1389-1431, 1549",
    },
    DivergenceClass {
        bucket: "splice_term_past_the_record_span",
        definition: "On a record whose alternate allele is longer than its reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine tract or intron term from positions past the end of the record's span; where no base of the span lies in a window of that term, vep-rs does not write it. For the intron term the window is the intron. With the record's alleles, from the first column or the input records, a row is set aside where the record is no pure insertion and one intron stands beside the span, at a place the row's cDNA positions allow, at which the record's regions give the reference's terms and the span gives vep-rs's, the donor region and splice_region_variant coming where the terms that withhold them went. Without the record a row is set aside where three or more positions lie in the transcript's sequence and vep-rs's row holds none of the six terms; a Location of two neighbouring positions, which an insertion's row shares, and a span that holds an intron's base then stay differences. A row that differs by the terms of a stop codon class as well is set aside where that class's rule takes the rest of the difference.",
        excludes: true,
        taxonomy_class: "splice_term_past_the_record_span",
        perl_citation: "VariationFeatureOverlapAllele.pm:405-460; BaseTranscriptVariationAllele.pm:117-118, 143-204",
    },
    DivergenceClass {
        bucket: "reference_skipped_record",
        definition: "A structural record the reference dropped before its input buffer, marked vep_skip for a span over --max_sv_size (its SVLEN, else its END), an ALT type with no Sequence Ontology term, an unsupported multi-allelic list or a deletion with neither END nor SVLEN, so it has no reference row in any format; vep-rs's rows for the record are set aside.",
        excludes: true,
        taxonomy_class: "reference_skipped_record",
        perl_citation: "Parser.pm:491-501, 821-827; Parser/VCF.pm:477-481, 555-558, 575",
    },
    DivergenceClass {
        bucket: "bnd_own_chromosome_rows_lost",
        definition: "A bracket breakend whose mate is on another chromosome loses the rows of its own end when a mate-chromosome transcript, loaded for an earlier record of the same buffer, reaches it first and names its slice: the reference keeps the mate rows and writes nothing under the own-end allele, an order-dependent result.",
        excludes: false,
        taxonomy_class: "bnd_own_chromosome_rows_lost",
        perl_citation: "InputBuffer.pm:191-207, 223-234; AnnotationSource.pm:235-258; AnnotationType/Transcript.pm:122-143",
    },
    DivergenceClass {
        bucket: "bnd_synonym_mate_rows_lost",
        definition: "A bracket breakend whose mate's chromosome is spelled with a prefix the mate-allele test does not resolve, Chr21 for 21, loses its mate rows in the reference and keeps them in vep-rs; a mate named by an accession on another chromosome loses them in the reference alone too, vep-rs resolving the accession through the synonym table, and one on the record's own chromosome keeps them in both engines.",
        excludes: false,
        taxonomy_class: "bnd_synonym_mate_rows_lost",
        perl_citation: "StructuralVariationOverlap.pm:86-102, 144-160; Utils/VariationEffect.pm:111-124; InputBuffer.pm:348-357; AnnotationType/Transcript.pm:134",
    },
    DivergenceClass {
        bucket: "fasta_named_slice_intergenic",
        definition: "Under --check_ref or --lookup_ref with a FASTA whose sequence names differ from the cache's, the variant's slice is named after the FASTA and every transcript fails the slice test, so the reference writes intergenic_variant where vep-rs writes transcript rows; a property of the run, counted and never set aside.",
        excludes: false,
        taxonomy_class: "fasta_named_slice_intergenic",
        perl_citation: "Parser.pm:749; BaseVEP.pm get_slice; AnnotationType/Transcript.pm:141",
    },
    DivergenceClass {
        bucket: "bnd_bare_reference_allele_row",
        definition: "A breakpoint record that names no mate (a bare <BND>, or one with SVLEN, POS2 or an unreadable END2 and no bracket) is held as N/<BND> and, with no breakend parsed, one allele is built per slash-separated part, the reference base N included, so the reference writes every row twice under chromosome_breakpoint, ALLELE_NUM 1 for the reference allele and 2 for the ALT; a reference allele has no variant effect, so vep-rs writes the ALT row once and the reference-allele rows are set aside: each duplicates the ALT row's tuple, so the unique-tuple key stream drops it before the raw and the adjusted figure, and both comparators count the rows so set aside per record.",
        excludes: true,
        taxonomy_class: "bnd_bare_reference_allele_row",
        perl_citation: "Parser/VCF.pm:499; StructuralVariationOverlap.pm:64-77; OutputFactory.pm:2029",
    },
    DivergenceClass {
        bucket: "bnd_mate_32bit_buffer_intergenic",
        definition: "A coordinate of the reference's input buffer that, plus the --distance flank, reaches 2^31 (a breakend's mate position, a span's end or the record's own position) turns every record of the buffer intergenic, that record and the records loaded with it (the first such position is 2147478648 under the default flank of 5000 and 2147483647 under --distance 1); the coordinates enter the buffer's bounds with the flank added when the cache regions are chosen. vep-rs annotates the buffer's records against their transcripts and a mate at its coordinate. Counted on the reference's rows at the record's Location and never set aside, since the records it reaches beyond that one depend on the buffer boundaries.",
        excludes: false,
        taxonomy_class: "bnd_mate_32bit_buffer_intergenic",
        perl_citation: "AnnotationSource.pm:190, 251-254, 306-309; InputBuffer.pm:386-390",
    },
    DivergenceClass {
        bucket: "gvcf_variant_site_dropped",
        definition: "A gVCF variant site lists <NON_REF> beside its real allele (`A G,<NON_REF>`); the <NON_REF> token routes the line to the structural path, where `G/<NON_REF>` has no Sequence Ontology term, and the reference drops the site whole, the real allele with it. vep-rs annotates the real allele and reads <NON_REF> as the absence of a variant. vep-rs's rows at the site are set aside: recognised here by the record's own name (`CHROM_POS_REF/ALT/<NON_REF>`) at a Location the reference wrote nothing for, and in the structural-variant comparator from the input record.",
        excludes: true,
        taxonomy_class: "gvcf_variant_site_dropped",
        perl_citation: "Parser/VCF.pm:237-245, 469-481; Parser.pm:663-728",
    },
    DivergenceClass {
        bucket: "dot_breakend_without_svtype_pair",
        definition: "A dot-notation breakend without INFO/SVTYPE (`N.`, `.N`) takes the sequence rules for its coordinates: the start is POS + 1 when the ALT's length differs from the REF's and POS when it does not, the end is POS + len(REF) - 1, so the reference writes the pair `POS-(POS+1)` with the positions of an insertion for a one-base REF, the span `(POS+1)-(POS+len(REF)-1)` for a longer REF beside an ALT of another length and `POS-(POS+len(REF)-1)` beside an ALT of its own length; with `SVTYPE=BND` it writes the single base at POS + 1, as vep-rs writes it for every REF. The reference's span row and vep-rs's point row at the span's start or the base after it on the same feature are set aside together.",
        excludes: true,
        taxonomy_class: "dot_breakend_without_svtype_pair",
        perl_citation: "BaseVCF4.pm:255-275, 297-317",
    },
    DivergenceClass {
        bucket: "missing_alt_annotated",
        definition: "A record whose ALT column is the missing value `.` is routed to the structural path by INFO/SVTYPE and annotated as an allele: a `chromosome_breakpoint` for `SVTYPE=BND`, a deletion over END for `SVTYPE=DEL`; without SVTYPE the reference writes nothing for it, as vep-rs writes nothing in both cases. The reference's rows are set aside by the structural-variant comparator, which reads the input; the SNP/indel comparator cannot tell them from a record with an allele.",
        excludes: false,
        taxonomy_class: "missing_alt_annotated",
        perl_citation: "BaseVCF4.pm:486-490; Parser/VCF.pm:237-245; Parser.pm:710-711",
    },
    DivergenceClass {
        bucket: "non_integer_end_annotated",
        definition: "An INFO/END, SVLEN or POS that is not an integer (`END=30425949.5`, `SVLEN=300.7`, `SVLEN=NaN`, POS `30425949.0` on the sequence path) is taken as written, passes `looks_like_number`, and the record is annotated with the text in its Location (`21:30425651-30425950.7`, `21:30425651-NaN`, `21:30425949.0`) and a fractional DISTANCE; vep-rs refuses an END that is not a coordinate and reads a SVLEN or POS at its numified value (`21:30425651-30425950`, the pair `21:30425650-30425651`, `21:30425949`). The reference's rows, recognised by the Location, are set aside in both comparators, and vep-rs's rows of the same record by the structural-variant comparator, which reads the input.",
        excludes: true,
        taxonomy_class: "non_integer_end_annotated",
        perl_citation: "BaseVCF4.pm:263, 302-308; Parser.pm:475",
    },
    DivergenceClass {
        bucket: "cnv_tr_empty_expansion_dropped",
        definition: "With a FASTA, a <CNV:TR> whose expansion yields no alternate allele, or one equal to the reference run (`RN=.`, a RUC Perl reads as true that gives zero or fewer repeats, a RUC times unit length equal to SVLEN), is expanded to an allele string with no alternate and the reference writes nothing for it, without a warning; vep-rs with a FASTA writes nothing for it either, and without one writes the symbolic rows over the run, as the reference does without one. vep-rs's rows are set aside by the structural-variant comparator when the input record has the shape and the reference wrote nothing at or beside the run; the SNP/indel comparator does not judge it.",
        excludes: false,
        taxonomy_class: "cnv_tr_empty_expansion_dropped",
        perl_citation: "Parser/VCF.pm:382-431, 516-543; Parser.pm:553-592",
    },
    DivergenceClass {
        bucket: "star_alone_svtype_annotated",
        definition: "A record whose every ALT is `*` under a true INFO/SVTYPE is routed to the structural path and annotated from INFO alone, a deletion over END for `SVTYPE=DEL` or a `chromosome_breakpoint` at POS + 1 for `SVTYPE=BND`, although `*` names no allele of the record; vep-rs writes no row for it. The reference's rows are set aside by the structural-variant comparator, which reads the input; the SNP/indel comparator cannot tell them from a record with an allele.",
        excludes: false,
        taxonomy_class: "star_alone_svtype_annotated",
        perl_citation: "Parser/VCF.pm:237-245, 469-475; StructuralVariationOverlap.pm:64-77",
    },
    DivergenceClass {
        bucket: "non_numeric_pos_numified",
        definition: "A POS that is not a non-negative number (`.`, empty, `0x1D04A1D`, fullwidth digits, `30_425_949`, `-5`) is numified by the structural path's increment, `.` to 0 and the record annotated from base 1, `30_425_949` from 31, `-5` from -4, where vep-rs refuses the line; the reference's rows at the record's start are set aside by the structural-variant comparator from the input, and by the SNP/indel comparator when the Location's first coordinate is negative, the one spelling no written POS reaches.",
        excludes: true,
        taxonomy_class: "non_numeric_pos_numified",
        perl_citation: "BaseVCF4.pm:255-264; Parser.pm:475-488",
    },
    DivergenceClass {
        bucket: "cnv_tr_rn_zero_allele_dropped",
        definition: "With a FASTA, a two-ALT <CNV:TR> of which one allele has `RN=0` loses an allele: the empty expansion is an empty part of the allele string, which is not read as multi-allelic, so the first allele alone is kept when the empty part leads and the trailing empty part is dropped by the split; vep-rs writes both alleles when the empty part leads, and the structural-variant comparator sets aside its rows of the dropped allele: the symbolic rows told by their direction term when the real allele is a gain, or unchanged beside a leading empty part (vep-rs reads every symbolic <CNV:TR> allele as a gain), and the literal rows a run with a FASTA writes for the real allele behind a leading empty part.",
        excludes: false,
        taxonomy_class: "cnv_tr_rn_zero_allele_dropped",
        perl_citation: "Parser/VCF.pm:411, 428; Parser.pm:881, 978",
    },
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::filter::SV_TRANSCRIPT_SELECTION_BUCKET;
    use crate::rules::EXCLUSION_REGISTRY;
    use regex::Regex;
    use std::collections::BTreeSet;

    /// The Python comparator's source as it stands.
    fn python_source() -> String {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../scripts/concordance/compare_vep_outputs.py"
        );
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"))
    }

    /// The first group of every match of `pattern`.
    fn captured(pattern: &str, source: &str) -> BTreeSet<String> {
        Regex::new(pattern)
            .expect("regex")
            .captures_iter(source)
            .map(|found| found[1].to_string())
            .collect()
    }

    fn names<'a>(buckets: impl IntoIterator<Item = &'a str>) -> BTreeSet<String> {
        buckets.into_iter().map(str::to_string).collect()
    }

    /// Fails with the classes only one side names.
    fn assert_same_classes(what: &str, python: &BTreeSet<String>, rust: &BTreeSet<String>) {
        let only_python: Vec<&String> = python.difference(rust).collect();
        let only_rust: Vec<&String> = rust.difference(python).collect();
        assert!(
            only_python.is_empty() && only_rust.is_empty(),
            "{what}: only the Python registry names {only_python:?}; only the Rust tool names \
             {only_rust:?}"
        );
    }

    /// The Python registry names a pair rule by a `bucket="..."` literal and every other
    /// class by a module-level string constant, one named `..._BUCKET` or one a rule takes as
    /// `bucket=NAME`; each side's names are compared as sets, so a class one program names
    /// and the other lacks fails here by name.
    #[test]
    fn the_python_registry_and_the_rust_table_name_the_same_classes() {
        let source = python_source();
        let python_pair_rules = captured(r#"\bbucket=["']([a-z0-9_]+)["']"#, &source);
        let referenced = captured(r"\bbucket=([A-Z][A-Z0-9_]*)\b", &source);
        let python_constants: BTreeSet<String> =
            Regex::new(r#"(?m)^([A-Z][A-Z0-9_]*) = ["']([a-z0-9_]+)["']$"#)
                .expect("regex")
                .captures_iter(&source)
                .filter(|found| found[1].ends_with("_BUCKET") || referenced.contains(&found[1]))
                .map(|found| found[2].to_string())
                .collect();
        let (pair_rules, one_sided) = DIVERGENCE_CLASSES_116_2.split_at(PAIR_RULE_COUNT);
        assert_same_classes(
            "pair rules of the registry",
            &python_pair_rules,
            &names(EXCLUSION_REGISTRY.iter().map(|rule| rule.bucket)),
        );
        assert_same_classes(
            "pair rules of the class table",
            &python_pair_rules,
            &names(pair_rules.iter().map(|class| class.bucket)),
        );
        assert_same_classes(
            "one-sided classes",
            &python_constants,
            &names(
                one_sided
                    .iter()
                    .map(|class| class.bucket)
                    .chain([SV_TRANSCRIPT_SELECTION_BUCKET]),
            ),
        );
        assert_eq!(python_pair_rules.len(), PAIR_RULE_COUNT);
        assert!(python_pair_rules.is_disjoint(&python_constants));
        assert!(
            referenced.len() >= one_sided.len(),
            "every one-sided rule takes its class as a constant"
        );
    }
}
