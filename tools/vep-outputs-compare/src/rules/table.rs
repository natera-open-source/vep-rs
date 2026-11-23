// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every divergence class of the 116.2 reference registry as the report states it.

/// One class: the one-sentence definition, whether the adjusted F1 sets it aside, its
/// taxonomy class and the Ensembl VEP source it cites.
#[derive(Clone, Copy)]
pub struct DivergenceClass {
    pub bucket: &'static str,
    pub definition: &'static str,
    pub excludes: bool,
    pub taxonomy_class: &'static str,
    /// The cited files, joined by `; `. A file is followed by a colon and its lines, joined by
    /// `, `: one line names a statement, and a range that cites a whole block runs to the line
    /// of its closing brace. A file followed by a blank and a name cites that sub whole, which
    /// holds where a line differs between the reference releases.
    pub perl_citation: &'static str,
}

/// The pair rules of the registry come first, in classification order.
pub const PAIR_RULE_COUNT: usize = 22;

/// Every class of the 116.2 registry as the report states it: pair rules in
/// classification order, then the one-sided rules. `data/divergence_classes.tsv` states the
/// same classes for the programs that read them, and the test below holds the two together.
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
        bucket: "structural_deletion_stop_lost_without_stop_codon",
        definition: "A Transcript row of a structural deletion on a transcript at whose coding end no stop codon stands, with cds_end_NF or without, or whose stop codon polyadenylation completes: the reference writes stop_lost for a span over the last three genomic positions of the coding region without testing a flag or reading a base, and vep-rs writes the row's other terms without it. The class reads the transcript's record in the candidate's cache and takes no pair without the cache.",
        excludes: true,
        taxonomy_class: "structural_deletion_stop_lost_without_stop_codon",
        perl_citation: "VariationEffect.pm:1268-1280",
    },
    DivergenceClass {
        bucket: "stop_codon_completed_by_polyadenylation",
        definition: "A Transcript row of a sequence allele on a transcript whose coding sequence ends on its last base one or two bases into its stop codon and whose record shows the model complete, by a mitochondrial codon table or a consensus coding sequence id: the reference reads those bases as an incomplete terminal codon and writes incomplete_terminal_codon_variant and no stop term, and vep-rs writes the terms of the row on the stop codon polyadenylation completes. The class derives vep-rs's coding terms from the transcript's record in the candidate's cache, the allele and the positions both rows print equal, and takes no pair without the cache; a deletion that runs past the transcript's last base is derived only beside the reference sequence.",
        excludes: true,
        taxonomy_class: "stop_codon_completed_by_polyadenylation",
        perl_citation: "VariationEffect.pm:1240, 1296, 1548, 1588-1614",
    },
    DivergenceClass {
        bucket: "stop_codon_overlap_without_stop_codon",
        definition: "A Transcript row of a sequence allele on a transcript without cds_end_NF at whose coding end no stop codon stands: the reference reads the last three coding positions as the stop codon and writes stop_lost or stop_retained_variant from its overlap tests, and vep-rs writes frameshift_variant, an in-frame term or coding_sequence_variant as the row's residues, codons and positions give them. The class derives vep-rs's coding terms from the cells both rows print equal, or the reference's where Amino_acids and Codons differ as alternate_codon_past_the_coding_end reads them, and the transcript's record in the candidate's cache, and takes no pair without the cache or, for a row under an identifier, without the input record.",
        excludes: true,
        taxonomy_class: "stop_codon_overlap_without_stop_codon",
        perl_citation: "VariationEffect.pm:1318, 1373-1381, 1400, 1408-1427, 1466-1480, 1508",
    },
    DivergenceClass {
        bucket: "stop_codon_behind_last_residue_without_stop_codon",
        definition: "A Transcript row of a sequence allele on a transcript at whose coding end no stop codon stands, with cds_end_NF or without: the reference reads a stop codon the allele brings directly behind the last residue as the transcript's own kept and writes stop_retained_variant beside the in-frame term, and vep-rs writes stop_gained beside it. The class derives vep-rs's coding terms from the cells both rows print equal, or the reference's where Amino_acids and Codons differ as alternate_codon_past_the_coding_end reads them, and the transcript's record in the candidate's cache, and takes no pair without the cache or, for a row under an identifier, without the input record.",
        excludes: true,
        taxonomy_class: "stop_codon_behind_last_residue_without_stop_codon",
        perl_citation: "VariationEffect.pm:1344, 1357",
    },
    DivergenceClass {
        bucket: "start_codon_overlap_without_start_codon",
        definition: "A Transcript row on a transcript without cds_start_NF at whose coding start no start codon is shown: the reference reads the first three coding positions as the start codon and writes start_lost or start_retained_variant, and vep-rs writes the terms the row's residues, codons and positions give without a start term. The class derives vep-rs's coding terms from the cells both rows print equal, or the reference's where Amino_acids and Codons differ as alternate_codon_past_the_coding_end reads them, and the transcript's record in the candidate's cache, and takes no pair without the cache or, for a row of a sequence allele under an identifier, without the input record.",
        excludes: true,
        taxonomy_class: "start_codon_overlap_without_start_codon",
        perl_citation: "VariationEffect.pm:975-986",
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
        definition: "A Transcript row of an insertion with both cDNA positions and a CDS start on which vep-rs writes stop_retained_variant, with or without inframe_insertion, for the stop codon the insertion keeps, and the reference coding_sequence_variant, or frameshift_variant with or without stop_gained: the reference's genomic stop codon overlap moves an insertion's start by the inserted length. A row at a Location of two neighbouring positions is an insertion's where its record's alleles, from the first column or the input records, write the Allele between the two positions and replace no base. Under an identifier a row whose input record replaces the two bases stays a difference; without the input records a row under an identifier is read as an insertion's by an Allele of another length than two, so two bases replaced by a sequence of another length are set aside with the insertions and an insertion of two bases stays a difference. The input records are read by Location and Allele: an insertion no record places stays a difference where a record replaces the same two bases by the same Allele.",
        excludes: true,
        taxonomy_class: "stop_codon_overlap_insertion_range",
        perl_citation: "VariationEffect.pm:1389-1431, 1549",
    },
    DivergenceClass {
        bucket: "stop_codon_settled_before_unknown_bases",
        definition: "A Transcript row with one Allele, one CDS position that carries both ends and one Codons cell on both sides, on which vep-rs reads the allele's peptide and the reference does not. vep-rs writes stop_gained where the reference writes none: added to the terms of a row the reference calls frameshift_variant, or in place of the reference's coding_sequence_variant, with at most one of inframe_insertion, inframe_deletion and protein_altering_variant. Where the reference writes stop_retained_variant for a stop codon test it does not make, vep-rs writes stop_gained in that term's place, with frameshift_variant, one of those three terms or no other term. Or vep-rs writes inframe_insertion beside stop_retained_variant where the reference writes no inframe_insertion: added to the terms of a row the reference calls stop_retained_variant, or with stop_retained_variant in the place of the reference's frameshift_variant or coding_sequence_variant on a row of an insertion between two neighbouring positions that carries both cDNA positions and a CDS start on both sides. The pair in the place of frameshift_variant is one of stop_codon_overlap_insertion_range as well, whose rule stands first and counts it. The alternate half of the Codons cell holds a letter other than A, C, G and T, and a stop codon of the standard genetic code stands among its whole codons before the first codon with such a letter. The reference translates no allele that holds such a letter and writes no Amino_acids; vep-rs writes the row every completion of the unknown bases shares, and the alternate half of its Amino_acids holds its first stop at that codon's place. A row that ends before its Codons column stays a difference.",
        excludes: true,
        taxonomy_class: "stop_codon_settled_before_unknown_bases",
        perl_citation: "TranscriptVariationAllele.pm:701; VariationEffect.pm:1100-1127, 1211-1232, 1290-1324, 1616-1649",
    },
    DivergenceClass {
        bucket: "frameshift_inside_a_held_stop_codon",
        definition: "A row of an insertion or a deletion whose length change is no multiple of three and whose first base lies inside the codon of a stop that the transcript's peptide holds as a stop and reads on past (a polymorphic pseudogene's), on which vep-rs writes frameshift_variant and the reference does not. Ensembl VEP's frameshift test returns 0 wherever the reference peptide of the allele's codons begins with a stop, which holds for the stop codon no residue follows, so the reference writes coding_sequence_variant alone, inframe_deletion alone, or inframe_insertion with coding_sequence_variant; vep-rs writes frameshift_variant in the place of that one group, and every other term stands on both rows. The two rows hold the same Amino_acids cell, whose reference half opens with a stop and whose alternate half holds the unknown residue, and the same Codons cell, whose two halves differ in length by a number of bases that is no multiple of three. A pair whose reference row holds stop_retained_variant or stop_lost stays a difference, and so does a pair with a row that ends before its Codons column.",
        excludes: true,
        taxonomy_class: "frameshift_inside_a_held_stop_codon",
        perl_citation: "VariationEffect.pm:1554",
    },
    DivergenceClass {
        bucket: "splice_term_past_the_record_span",
        definition: "On a record whose alternate allele is longer than its reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine tract or intron term from positions past the end of the record's span; where no base of the span lies in a window of that term, vep-rs does not write it, unless it is a term of an insertion written after the span. For the intron term the window is the intron. With the record's alleles, from the first column or the input records, a row is set aside where the record is no pure insertion and one intron stands beside the span, at a place the row's cDNA positions allow, at which the record's regions give the reference's terms and the span gives vep-rs's, the donor region and splice_region_variant coming where the terms that withhold them went. A record whose alternate allele opens with its whole reference allele of two or more bases changes no base of its span. vep-rs reads it as the insertion it writes, between the span's last base and the base after it, and writes that insertion's terms. With the record's alleles such a row is set aside where that insertion gives vep-rs's terms at the intron's place, the polypyrimidine tract and intron terms among those vep-rs alone may hold. A pair whose allele opens with its whole reference allele of two or more bases is put to the row-aware rules also where a set rule classifies it without setting it aside, and is counted under the row-aware rule that takes it. Without the record a row is set aside where three or more positions lie in the transcript's sequence and vep-rs's row holds none of the six terms; a Location of two neighbouring positions, which an insertion's row shares, and a span that holds an intron's base then stay differences, and no row is read as that of an allele that opens with its whole reference allele. A row that differs by the terms of a class of a coding end or of a stop codon class as well is set aside where that class's rule takes the rest of the difference.",
        excludes: true,
        taxonomy_class: "splice_term_past_the_record_span",
        perl_citation: "VariationFeatureOverlapAllele.pm:405-460; BaseTranscriptVariationAllele.pm:117-118, 143-206",
    },
    DivergenceClass {
        bucket: "pachysolen_nuclear_code_read_as_standard",
        definition: "A Transcript row of a transcript on genetic code 26, the nuclear code of Pachysolen tannophilus, whose Codons cell holds CTG as its reference or its alternate codon and a leucine codon as the other: the reference translates under the last table of its codon module, which holds no table 26, and writes synonymous_variant with leucine alone in Amino_acids; vep-rs reads CTG as alanine, as the code gives it, and writes missense_variant with alanine on that codon's side and leucine on the other. A pair is set aside where the two sets are otherwise equal and the two rows hold the same Codons cell; a row that ends before its Codons column stays a difference.",
        excludes: true,
        taxonomy_class: "pachysolen_nuclear_code_read_as_standard",
        perl_citation: "Bio/Tools/CodonTable.pm:328-338, 412-460; TranscriptVariationAllele.pm:712, 736",
    },
    DivergenceClass {
        bucket: "start_codon_removed_without_five_prime_utr",
        definition: "A Transcript row of a deletion of whole codons from the first coding base of a transcript with no 5' UTR, on which the reference writes inframe_deletion without start_lost and vep-rs start_lost beside it for the start codon the record removes.",
        excludes: true,
        taxonomy_class: "start_codon_removed_without_five_prime_utr",
        perl_citation: "VariationEffect.pm:851-904, 906-945",
    },
    DivergenceClass {
        bucket: "start_codon_kept_without_five_prime_utr",
        definition: "A Transcript row of an insertion or deletion that overlaps the start codon of a transcript with no 5' UTR and leaves a start codon at the start of the coding sequence, on which the reference writes no start_retained_variant, and start_lost where its length or peptide test holds, and vep-rs start_retained_variant alone.",
        excludes: true,
        taxonomy_class: "start_codon_kept_without_five_prime_utr",
        perl_citation: "VariationEffect.pm:851-904, 947-963, 1028-1075",
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
        perl_citation: "StructuralVariationOverlap.pm:86-103, 144-160; Utils/VariationEffect.pm:111-124; InputBuffer.pm:348-357; AnnotationType/Transcript.pm:134",
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
        perl_citation: "Parser/VCF.pm:499; StructuralVariationOverlap.pm:64-77; OutputFactory.pm BaseStructuralVariationOverlapAllele_to_output_hash",
    },
    DivergenceClass {
        bucket: "bnd_mate_32bit_buffer_intergenic",
        definition: "A coordinate of the reference's input buffer that, plus the --distance flank, reaches 2^31 (a breakend's mate position, a span's end or the record's own position) turns every record of the buffer intergenic, that record and the records loaded with it (the first such position is 2147478648 under the default flank of 5000 and 2147483647 under --distance 1); the coordinates enter the buffer's bounds with the flank added when the cache regions are chosen. vep-rs annotates the buffer's records against their transcripts and a mate at its coordinate. Counted on the reference's rows at the record's Location and never set aside, since the records it reaches beyond that one depend on the buffer boundaries.",
        excludes: false,
        taxonomy_class: "bnd_mate_32bit_buffer_intergenic",
        perl_citation: "AnnotationSource.pm:190, 251-257, 306-309; InputBuffer.pm:386-390",
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
        perl_citation: "BaseVCF4.pm:255-276, 297-317",
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
        perl_citation: "BaseVCF4.pm:263, 302-309; Parser.pm:475",
    },
    DivergenceClass {
        bucket: "cnv_tr_empty_expansion_dropped",
        definition: "With a FASTA, a <CNV:TR> whose expansion yields no alternate allele, or one equal to the reference run (`RN=.`, a RUC Perl reads as true that gives zero or fewer repeats, a RUC times unit length equal to SVLEN), is expanded to an allele string with no alternate and the reference writes nothing for it, without a warning; vep-rs with a FASTA writes nothing for it either, and without one writes the symbolic rows over the run, as the reference does without one. vep-rs's rows are set aside by the structural-variant comparator when the input record has the shape and the reference wrote nothing at or beside the run; the SNP/indel comparator does not judge it.",
        excludes: false,
        taxonomy_class: "cnv_tr_empty_expansion_dropped",
        perl_citation: "Parser/VCF.pm:382-431, 516-544; Parser.pm:553-592",
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
        definition: "With a FASTA, a two-ALT <CNV:TR> of which one allele has `RN=0` loses an allele: the empty expansion is an empty part of the allele string, which is not read as multi-allelic, so the first allele alone is kept when the empty part leads and the trailing empty part is dropped by the split; vep-rs writes both alleles when the empty part leads, and the structural-variant comparator sets aside its rows of the dropped allele: the symbolic rows told by their direction term when the real allele is a gain, or unchanged beside a leading empty part (vep-rs reads a symbolic <CNV:TR> allele as a gain unless its record states it shorter than its run or, with no direction term, as long as it), and the literal rows a run with a FASTA writes for the real allele behind a leading empty part.",
        excludes: false,
        taxonomy_class: "cnv_tr_rn_zero_allele_dropped",
        perl_citation: "Parser/VCF.pm:411, 428; Parser.pm:881, 978",
    },
];

/// The class of a comparison that names the converted cache: it is applied, and stated in
/// the report after the classes of `DIVERGENCE_CLASSES_116_2`, only there, since its rows
/// are told from the cache's transcripts and not from the two outputs.
pub const MERGED_SAME_DBID_COLLAPSE: DivergenceClass = DivergenceClass {
    bucket: "merged_same_dbid_collapse",
    definition: "Under --merged the reference keeps one transcript per dbID among those an input buffer loads, and a merged cache numbers its two transcript sets apart, so a transcript whose dbID a transcript of the other set carries loses its rows wherever the buffer read the other first, and a record left with no transcript is written as intergenic_variant; vep-rs keeps both transcripts, and a comparison that names the converted cache replays the reference's buffers from the reference's rows and that cache and sets aside vep-rs's rows of each transcript the pass dropped and the reference's intergenic row at a record left with none.",
    excludes: true,
    taxonomy_class: "merged_same_dbid_collapse",
    perl_citation: "AnnotationType/Transcript.pm:283-287; AnnotationSource.pm:109-145; Pipeline/DumpVEP/MergeVEP.pm:94-120",
};

/// The second class of a comparison that names the converted cache, stated in the report
/// after `MERGED_SAME_DBID_COLLAPSE`: its rows are told from the same replay of the
/// reference's buffers over the cache's transcripts.
pub const MERGED_SAME_STABLE_ID_COLLAPSE: DivergenceClass = DivergenceClass {
    bucket: "merged_same_stable_id_collapse",
    definition: "Under --refseq or --merged the reference keeps one transcript per stable identifier among those an input buffer loads, the ensembl-sourced one or else the one of the lowest dbID, so a transcript that shares its identifier with another the buffer loads loses its rows, and a record left with no transcript is written as intergenic_variant; vep-rs keeps every transcript the cache tells apart by source and dbID, and a comparison that names the converted cache replays the reference's buffers from the reference's rows and that cache and sets aside vep-rs's rows of each transcript the pass dropped and the reference's intergenic row at a record left with none.",
    excludes: true,
    taxonomy_class: "merged_same_stable_id_collapse",
    perl_citation: "AnnotationType/Transcript.pm:275, 303-311, 314-343; AnnotationSource.pm:109-145",
};

/// The class of a comparison that names the reference's cache, where that cache's synonyms
/// file gives a name to two or more sequences the cache holds: it is applied, and stated in
/// the report after the classes of `DIVERGENCE_CLASSES_116_2`, only there, since its records
/// are told from the cache's files and not from the two outputs.
pub const SYNONYM_OF_SEVERAL_SEQUENCES_BY_HASH_ORDER: DivergenceClass = DivergenceClass {
    bucket: "synonym_of_several_sequences_by_hash_order",
    definition: "A record whose sequence name the cache holds no sequence by and the cache's synonyms file gives to two or more sequences it holds, so that the file does not say which of them the name stands for: the reference reads the record on whichever of them comes first in the order Perl returns the keys of the name's synonyms, which changes between runs of one command, so its rows for the record are not reproducible; vep-rs skips the record with a warning. With the reference's cache named, whose directories and synonyms file say which names are such, the reference's rows of such a record are set aside where the reference output's header states a command line that names no --synonyms, so that the run read the cache's own file, and vep-rs wrote no row at the record's Location; a record vep-rs wrote a row for stays a difference. Without the cache the class is not applied.",
    excludes: true,
    taxonomy_class: "synonym_of_several_sequences_by_hash_order",
    perl_citation: "BaseVEP.pm:646-656, 711-716; AnnotationSource.pm:179; AnnotationSource/Cache.pm:90-92",
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::filter::SV_TRANSCRIPT_SELECTION_BUCKET;
    use crate::rules::{row_pair_rules, ExclusionRule, Predicate, EXCLUSION_REGISTRY};
    use std::collections::BTreeSet;

    /// The columns of the class file, in the order of its header line.
    const COLUMNS: [&str; 11] = [
        "registry",
        "order",
        "bucket",
        "reference_releases",
        "reads",
        "classifies",
        "excludes",
        "excludes_with_cache",
        "taxonomy_class",
        "perl_citation",
        "definition",
    ];

    /// The registries of the class file, in the order of its lines.
    const REGISTRIES: [&str; 5] = [
        "pair",
        "row_pair_116.2",
        "one_sided_115.2",
        "one_sided_116.2",
        "cache_116.2",
    ];

    /// The class file, `data/divergence_classes.tsv` of this tool, by a path that holds from
    /// the root of every crate under `tools/` that compiles this module.
    fn class_file() -> String {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../vep-outputs-compare/data/divergence_classes.tsv"
        );
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"))
    }

    /// One class as a line of the class file states it.
    struct Listed<'a> {
        registry: &'a str,
        order: usize,
        bucket: &'a str,
        reference_releases: &'a str,
        reads: &'a str,
        classifies: Option<bool>,
        excludes: bool,
        excludes_with_cache: Option<bool>,
        taxonomy_class: &'a str,
        perl_citation: &'a str,
        definition: &'a str,
    }

    /// A flag of the class file: `true`, `false`, or empty where the rules of the class's
    /// registry carry no such flag.
    fn flag(cell: &str) -> Option<bool> {
        match cell {
            "" => None,
            "true" => Some(true),
            "false" => Some(false),
            other => panic!("a flag is `true`, `false` or empty, not {other:?}"),
        }
    }

    fn listed(line: &str) -> Listed<'_> {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), COLUMNS.len(), "the fields of the line {line}");
        let cell = |column: &str| {
            let at = COLUMNS.iter().position(|name| *name == column);
            cells[at.unwrap_or_else(|| panic!("no column {column}"))]
        };
        let bucket = cell("bucket");
        Listed {
            registry: cell("registry"),
            order: cell("order")
                .parse()
                .unwrap_or_else(|e| panic!("{bucket}: {e}")),
            bucket,
            reference_releases: cell("reference_releases"),
            reads: cell("reads"),
            classifies: flag(cell("classifies")),
            excludes: flag(cell("excludes")).unwrap_or_else(|| panic!("{bucket}: no `excludes`")),
            excludes_with_cache: flag(cell("excludes_with_cache")),
            taxonomy_class: cell("taxonomy_class"),
            perl_citation: cell("perl_citation"),
            definition: cell("definition"),
        }
    }

    /// What a pair rule is asked on, as the class file's `reads` names it.
    fn reads(rule: &ExclusionRule) -> &'static str {
        match rule.predicate {
            Predicate::Sets(_) => "sets",
            Predicate::SetsAndAllele(_) => "sets_and_allele",
            Predicate::Rows(_) => "rows",
            Predicate::RowsAtACodingEnd(_) => "rows_and_coding_ends",
        }
    }

    fn of_release<'a, 'b>(classes: &'b [Listed<'a>], release: &str) -> Vec<&'b Listed<'a>> {
        let held = |class: &&Listed| class.reference_releases.split(',').any(|r| r == release);
        classes.iter().filter(held).collect()
    }

    /// The class file and the rule tables state the same classes: every registry of the
    /// file in the order of its lines, a class's place counted from 1; the pair rules of
    /// `EXCLUSION_REGISTRY` in classification order with what each reads and whether it
    /// classifies and excludes, the row-aware ones under the 116.2 release alone; the classes
    /// of the 116.2 summary in its order, each with its definition, exclusion flag, taxonomy
    /// class and citation, the two classes of a comparison that names the converted cache and
    /// the class of a comparison that reads the cache's synonyms last, in that order; and the
    /// one class of the 115.2 one-sided pass. A class, a place, a flag, a citation or a definition
    /// that one of them states and the other does not fails here by the class's name.
    #[test]
    fn the_class_file_and_the_rule_tables_state_the_same_classes() {
        let file = class_file();
        assert!(file.ends_with('\n') && !file.contains('\r') && !file.contains("\n\n"));
        let mut lines = file.lines();
        assert_eq!(lines.next(), Some(COLUMNS.join("\t").as_str()));
        let classes: Vec<Listed> = lines.map(listed).collect();

        let mut registries: Vec<&str> = classes.iter().map(|class| class.registry).collect();
        registries.dedup();
        assert_eq!(registries, REGISTRIES);
        for registry in REGISTRIES {
            let of_registry = classes.iter().filter(|class| class.registry == registry);
            let places: Vec<usize> = of_registry.map(|class| class.order).collect();
            assert_eq!(places, (1..=places.len()).collect::<Vec<_>>(), "{registry}");
        }
        let buckets: BTreeSet<&str> = classes.iter().map(|class| class.bucket).collect();
        assert_eq!(buckets.len(), classes.len(), "a class stands on one line");

        let (pair_rules, others): (Vec<&Listed>, Vec<&Listed>) = classes
            .iter()
            .partition(|class| class.registry == "pair" || class.registry == "row_pair_116.2");
        let (one_sided, of_the_cache): (Vec<&Listed>, Vec<&Listed>) = others
            .into_iter()
            .partition(|class| class.registry.starts_with("one_sided_"));
        let rule_buckets: Vec<&str> = EXCLUSION_REGISTRY.iter().map(|rule| rule.bucket).collect();
        let listed_buckets: Vec<&str> = pair_rules.iter().map(|class| class.bucket).collect();
        assert_eq!(
            listed_buckets, rule_buckets,
            "the pair rules, in classification order"
        );
        assert_eq!(pair_rules.len(), PAIR_RULE_COUNT);
        for (class, rule) in pair_rules.iter().zip(&EXCLUSION_REGISTRY) {
            let (registry, releases) = match rule.row_aware() {
                true => ("row_pair_116.2", "116.2"),
                false => ("pair", "115.2,116.2"),
            };
            assert_eq!(class.registry, registry, "{}", rule.bucket);
            assert_eq!(class.reference_releases, releases, "{}", rule.bucket);
            assert_eq!(class.reads, reads(rule), "{}", rule.bucket);
            assert_eq!(class.classifies, Some(rule.classifies), "{}", rule.bucket);
            assert_eq!(class.excludes, rule.excludes, "{}", rule.bucket);
            assert_eq!(class.excludes_with_cache, None, "{}", rule.bucket);
        }
        let row_aware: Vec<&str> = row_pair_rules().iter().map(|rule| rule.bucket).collect();
        let listed_row_aware = classes.iter().filter(|c| c.registry == "row_pair_116.2");
        assert_eq!(
            listed_row_aware.map(|c| c.bucket).collect::<Vec<_>>(),
            row_aware
        );

        for class in &one_sided {
            let release = class.registry.trim_start_matches("one_sided_");
            assert_eq!(class.reference_releases, release, "{}", class.bucket);
            assert_eq!(class.reads, "record", "{}", class.bucket);
            assert_eq!(class.classifies, None, "{}", class.bucket);
            assert!(class.excludes_with_cache.is_some(), "{}", class.bucket);
        }
        assert_eq!(of_the_cache.len(), 3, "the classes that read the cache");
        for class in &of_the_cache {
            assert_eq!(class.registry, "cache_116.2", "{}", class.bucket);
            assert_eq!(class.reference_releases, "116.2", "{}", class.bucket);
            assert_eq!(class.reads, "cache", "{}", class.bucket);
            assert_eq!(class.classifies, None, "{}", class.bucket);
            assert_eq!(class.excludes_with_cache, None, "{}", class.bucket);
        }

        let stated: Vec<DivergenceClass> = DIVERGENCE_CLASSES_116_2
            .iter()
            .chain([
                &MERGED_SAME_DBID_COLLAPSE,
                &MERGED_SAME_STABLE_ID_COLLAPSE,
                &SYNONYM_OF_SEVERAL_SEQUENCES_BY_HASH_ORDER,
            ])
            .copied()
            .collect();
        let of_116_2 = of_release(&classes, "116.2");
        let listed_116_2: Vec<&str> = of_116_2.iter().map(|class| class.bucket).collect();
        let stated_116_2: Vec<&str> = stated.iter().map(|class| class.bucket).collect();
        assert_eq!(
            listed_116_2, stated_116_2,
            "the classes of the 116.2 summary, in its order"
        );
        for (class, stated) in of_116_2.iter().zip(&stated) {
            assert_eq!(class.definition, stated.definition, "{}", stated.bucket);
            assert_eq!(class.excludes, stated.excludes, "{}", stated.bucket);
            assert_eq!(
                class.taxonomy_class, stated.taxonomy_class,
                "{}",
                stated.bucket
            );
            assert_eq!(
                class.perl_citation, stated.perl_citation,
                "{}",
                stated.bucket
            );
        }

        let listed_115_2: Vec<&str> = of_release(&classes, "115.2")
            .iter()
            .map(|class| class.bucket)
            .collect();
        let held_115_2: Vec<&str> = EXCLUSION_REGISTRY
            .iter()
            .filter(|rule| !rule.row_aware())
            .map(|rule| rule.bucket)
            .chain([SV_TRANSCRIPT_SELECTION_BUCKET])
            .collect();
        assert_eq!(listed_115_2, held_115_2, "the classes of the 115.2 release");
    }

    /// What a citation states wrongly, or nothing: its parts are joined by `; `, and a part
    /// is a `.pm` file followed by a colon and lines joined by `, ` (a line, or two lines
    /// joined by a hyphen, the first below the second), or by a blank and the name of a sub.
    fn citation_fault(citation: &str) -> Option<String> {
        let number = |text: &str| -> Option<u32> {
            let plain = !text.is_empty()
                && !text.starts_with('0')
                && text.bytes().all(|b| b.is_ascii_digit());
            if plain {
                text.parse().ok()
            } else {
                None
            }
        };
        for part in citation.split("; ") {
            let Some((file, cited)) = part.split_once(".pm") else {
                return Some(format!("`{part}` names no .pm file"));
            };
            let named = file.starts_with(|c: char| c.is_ascii_uppercase())
                && file.chars().all(|c| c.is_ascii_alphanumeric() || c == '/');
            if !named {
                return Some(format!("`{part}` opens with no file name"));
            }
            if let Some(sub) = cited.strip_prefix(' ') {
                let word = sub.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_')
                    && sub.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
                if !word {
                    return Some(format!("`{part}` names no sub after its file"));
                }
                continue;
            }
            let Some(lines) = cited.strip_prefix(':') else {
                return Some(format!("`{part}` has neither lines nor a sub"));
            };
            for item in lines.split(", ") {
                let whole = match item.split_once('-') {
                    None => number(item).is_some(),
                    Some((first, last)) => {
                        matches!((number(first), number(last)), (Some(a), Some(b)) if a < b)
                    }
                };
                if !whole {
                    return Some(format!("`{item}` of `{part}` is no line and no range"));
                }
            }
        }
        None
    }

    #[test]
    fn every_citation_range_is_written_low_to_high_with_its_file() {
        for class in DIVERGENCE_CLASSES_116_2 {
            assert_eq!(
                citation_fault(class.perl_citation),
                None,
                "{}",
                class.bucket
            );
        }
        let file = class_file();
        for class in file.lines().skip(1).map(listed) {
            assert_eq!(
                citation_fault(class.perl_citation),
                None,
                "{}",
                class.bucket
            );
        }
        for whole in [
            "Parser.pm:749; BaseVEP.pm get_slice; AnnotationType/Transcript.pm:141",
            "Parser/VCF.pm:382-431, 516-544",
        ] {
            assert_eq!(citation_fault(whole), None, "{whole}");
        }
        for faulty in [
            "BaseVCF4.pm:276-255",
            "BaseVCF4.pm:255-255",
            "BaseVCF4.pm:255-",
            "BaseVCF4.pm:255,297",
            "BaseVCF4.pm:0255",
            "BaseVCF4.pm:",
            "BaseVCF4.pm",
            "BaseVCF4.pm get start",
            "lines 255-276",
            "baseVCF4.pm:255",
            "Parser.pm:749;BaseVEP.pm get_slice",
        ] {
            assert!(citation_fault(faulty).is_some(), "{faulty}");
        }
    }
}

#[cfg(test)]
#[path = "table_definitions_tests.rs"]
mod definitions_tests;

#[cfg(test)]
#[path = "table_reference_classes_tests.rs"]
mod reference_classes_tests;
