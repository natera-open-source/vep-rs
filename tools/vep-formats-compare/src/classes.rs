// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class names and the class tables.

use std::sync::LazyLock;

/// The class a vep_rs_only_tuples entry carries when every record it names is one the
/// reference dropped before annotation.
pub const REFERENCE_SKIPPED_RECORD: &str = "reference_skipped_record";
/// The built-in mask class dropping the reference's regulatory rows (`without_regulatory_rows`).
pub const REGULATORY_FEATURE_ROWS: &str = "regulatory_feature_rows";
/// The built-in field class `motif_score_reverse_strand_mnv`: MOTIF_SCORE_CHANGE on a
/// reverse-strand MotifFeature entry whose allele is two or more bases, which Ensembl VEP
/// writes from the motif position of the variant's genomic start, the allele's last position
/// along a reverse-strand motif (MotifFeatureVariationAllele.pm:245), so its score, or its
/// absence, differs from the score of the allele written over the positions it covers. An
/// insertion or a length-changing allele carries no score on either side, so a difference on
/// such an allele is this class alone.
pub const MOTIF_SCORE_REVERSE_STRAND_MNV: &str = "motif_score_reverse_strand_mnv";
/// The key under which a manifest entry names the class of a second cause that changes the value
/// it documents: the entry is applied only where both of its classes are named.
pub const SECOND_EXPECTED_DIVERGENCE: &str = "second_expected_divergence";
/// The recoder format's two other built-in classes: the var_synonyms source the cache's column
/// never carries (the merged dbSNP identifiers), and the keys the recoder writes unshifted for an
/// identifier input of an insertion or deletion (TranscriptVariationAllele.pm `_return_3prime`
/// shifts only a feature whose class is insertion or deletion, and the database row's stored
/// class, `indel` for a non-minimal allele string, survives Parser.pm `minimise_alleles`).
pub const ARCHIVE_DBSNP_SYNONYMS: &str = "archive_dbsnp_synonyms";
/// The built-in class under which a set-valued column agrees when its items agree as a
/// multiset: Ensembl VEP prints DOMAINS in the order of the cache's protein_features array
/// (ensembl-variation BaseTranscriptVariation.pm:634-643 at release/116 copies it from the
/// core API's Translation.pm:772-819, which buckets a translation's features by analysis and
/// flattens the buckets in Perl hash order, the database's rows inside each), FREQS in hash
/// order (OutputFactory.pm:1152-1153) and the ZYG list of --individual_zyg in the hash order
/// of its samples (OutputFactory.pm:983-987), so none of the orders carries meaning.
/// Under the class two JSON `domains` arrays agree where they can be the arrays of one
/// DOMAINS list written in two orders (`arrays_of_one_domains_list`): an entry the pattern
/// of Ensembl VEP's JSON writer does not match adds the object of the last entry that
/// matched once more, or none where none has matched (JSON.pm:363-371), so two orders of
/// one list can give arrays that differ in the object repeated or in the number of objects.
pub const ORDER_ONLY: &str = "order_only";
/// The built-in class under which an individual's alleles agree as a set. Under --individual,
/// Parser/VCF.pm:635-639 at release/116.3 collects the alternate alleles an individual's
/// genotype names in a hash and writes the copy's allele string as the reference allele
/// followed by the hash's keys, so a copy carrying two alternate alleles lists them in hash
/// order, and that order runs through the copy's rows within each feature, its ALLELE_NUM
/// values (the allele's index in the string, BaseVariationFeatureOverlapAllele.pm:217-218 at
/// release/116) and the JSON record's allele_string (JSON.pm:241). UPLOADED_ALLELE follows it
/// only on a record with no original allele string (OutputFactory.pm:1360 prefers that
/// string, which minimise_alleles sets on an indel record and under --minimal, Parser.pm:852,
/// 881, 978-982).
pub const INDIVIDUAL_ALLELE_ORDER: &str = "individual_allele_order";
/// The built-in class spelling the reference's conservation key as the candidate does. Each
/// `--phyloP <N>` expands into a --custom source whose short name is `phlyoP<N>way` while the
/// file it names is spelt `phyloP<N>way` (Config.pm:555 at release/116.3), so the misspelt key
/// heads the reference's header line, tab column, Extra item, CSQ field and JSON
/// `custom_annotations` member where the candidate writes `phyloP<N>way`; the class renames
/// the key throughout the reference text before it is parsed.
pub const CONSERVATION_SHORT_NAME_SPELLING: &str = "conservation_short_name_spelling";
/// The built-in class mapping the reference's two descriptions of a custom source's `PC` field
/// onto the shares the field reports. AnnotationSourceAdaptor.pm:395-397 at release/116.3 names
/// the field under reciprocal=0 as the share of the input variant covered and under reciprocal=1
/// as a minimum reciprocal overlap, while File.pm:556-564 computes the share of the annotation
/// record the variant covers without reciprocal and the share of the variant the record covers
/// with it. The candidate's header names the share computed; the class puts that text into the
/// reference's `## <key>_PC : <text> from <file>` and `##INFO=<ID=<key>_PC,...>` lines.
pub const PC_DESCRIPTION_TEXT: &str = "pc_description_text";
/// The two pair classes of the SNP/indel comparator's registry (`EXCLUSION_REGISTRY` of
/// `tools/vep-outputs-compare`, whose predicates decide them here as well), built in so that a
/// population case, which has no manifest beside its input, sets aside the pairs the
/// concordance figures set aside: a key both engines write whose consequence sets differ by
/// the one shape each rule admits. `start_cooccurrence_swap` is Ensembl VEP's emission of
/// `start_lost` and `start_retained_variant` together on one allele (Utils/VariationEffect.pm
/// :851-881, :947-963 and :1028-1075), `splice_family_swap` its loss of `splice_region_variant`,
/// with `intron_variant` or alone, to the assignment at BaseTranscriptVariationAllele.pm:215.
/// A matched pair is skipped whole, as a manifest's documented key is: its terms, IMPACT and
/// HGVSp follow the term the engines disagree on.
pub const START_COOCCURRENCE_SWAP: &str = "start_cooccurrence_swap";
pub const SPLICE_FAMILY_SWAP: &str = "splice_family_swap";
/// The built-in field classes of the corrected behaviours, each one shape of disagreement in
/// the fields of one feature or plugin, applied without a manifest; every other field of the
/// entry compares exactly. A value is empty when it is `-`, the empty string or absent.
///
/// `gxf_reverse_strand_utr_sequence`: on a transcript Ensembl VEP reads from a GTF or GFF
/// file through --gtf or --gff, BaseTranscriptVariation.pm:1101-1113 (`_generic_utr`)
/// transfers the transcript onto its feature slice and reads the UTR from the copy; Exon.pm:985
/// drops the exons' cached bases on the transfer, and FastaSequence.pm:369-371 maps a
/// reverse-strand slice's relative coordinates from its genomic start, where the core API counts
/// them from its end, so a reverse-strand transcript's 3' UTR is read as the reverse complement
/// of the span mirrored about the slice. A cache transcript carries its UTR from the dumper and
/// never reaches the read. The stop codon a frameshift or extension runs into lies in that UTR,
/// so the count after `fsTer` or `extTer` of the HGVSp differs on a `STRAND` -1 row and
/// nothing else does. The rule takes a row whose `SOURCE` in the reference names the file,
/// the stamp BaseGXF.pm:635 gives such a transcript alone: a cache transcript's row, with no
/// `SOURCE` or with the `Ensembl` or `RefSeq` of a merged cache, stays a difference.
pub const GXF_REVERSE_STRAND_UTR_SEQUENCE: &str = "gxf_reverse_strand_utr_sequence";
/// `alternate_codon_past_the_coding_end`: TranscriptVariationAllele.pm:2421-2433 build the
/// alternate coding sequence, empty it where it is shorter than one codon (`_trim_incomplete_codon`,
/// :2533-2539, whose test assigns) and join the 3' UTR onto it; :877 cuts the alternate codon
/// from that string, and the reference codon from the coding sequence alone (:868). Where the
/// coding sequence ends inside the codon window the alternate codon so holds one or two bases
/// of the UTR the reference codon lacks, and where a deletion leaves one or two coding bases
/// it is UTR bases alone. The candidate cuts both codons from coding bases. The first shape:
/// two Codons values with equal reference halves that are no whole number of codons, the
/// reference's alternate half the candidate's with one or two lower-case bases behind it, no
/// more than the last codon lacks; and on such a row two Amino_acids values that open with
/// the same reference residues, each alternate half the translation of its own alternate
/// codons under the standard genetic code. The second shape: two Codons values with equal
/// reference halves of whole codons on a row whose protein positions open at residue 1, the
/// candidate's alternate half the one or two lower-case bases of the reference half, the
/// reference's as many lower-case bases that spell something else. An upper-case base
/// behind the candidate's half, a longer tail, differing reference halves and residues
/// that are no translation of the row's codons stay differences.
pub const ALTERNATE_CODON_PAST_THE_CODING_END: &str = "alternate_codon_past_the_coding_end";
/// `organelle_peptide_paths_translation_table`: TranscriptVariationAllele.pm:2333, :2450,
/// :2492 and :2555 translate a frameshift's new frame, the reference an inserted peptide is
/// tested against, the sequence a stop distance is counted in and an empty-window deletion's
/// peptide with `translate()` bare, the standard code, against a reference peptide of the
/// transcript's codon table; the candidate translates both sides under that table. The rows do
/// not state the table, so the class reads it from the candidate JSON cache
/// (`OrganelleCodonTables`, `--json-cache`). The shape: an HGVSp on a row whose transcript
/// the cache files under a table that reads some codon unlike the standard one, both values
/// naming one protein, the row's consequence set and the two descriptions placing the
/// difference in one of the paths, and nothing else of the row's coding change differing
/// (`organelle_peptide_path`). Without the cache the class applies from a manifest alone.
pub const ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE: &str =
    "organelle_peptide_paths_translation_table";
/// `conservation_max_text_order`: Conservation.pm:242-243 sorts a span's covered scores as
/// text under MAX and reports the last, so a span whose scores are all negative reports the
/// most negative (`-5.2` sorts after `-2.17`); the candidate reports the numeric maximum.
/// The shape: both values numbers, the reference's below the candidate's, the candidate's below
/// zero.
pub const CONSERVATION_MAX_TEXT_ORDER: &str = "conservation_max_text_order";
/// `gwas_haplotype_cross_product`: GWAS.pm:338-340 hands every location of a haplotype or
/// interaction row (`CHR_ID` `21;21;...`, `SNPS` `rs1; rs2; ...`) to every SNP of the
/// row, so a site carries every SNP's risk allele and study; the candidate pairs SNP i with
/// location i. The reference annotates a superset of sites, so the shape is a `GWAS_` field
/// the reference fills and the candidate leaves empty, or one both fill differently. A field
/// the candidate alone fills is the class's where both rows are annotated under different
/// studies: the reference reports another SNP's catalog row, which leaves the field empty, and
/// the candidate the site's own; on any other row it stays a difference. A field the catalog row
/// leaves empty is a defined empty string on an annotated reference row (`GWAS_beta_coef=`)
/// and no key on a candidate row the plugin never annotated: the pair is the class's where the
/// reference row fills another `GWAS_` field and a difference where it fills none. The class
/// reads the catalog's rows of two or more SNPs, as the interaction class below says. A row of
/// several locations whose SNPs differ from them in count is skipped by the candidate and
/// written whole by the reference, so every SNP's allele at every location of it is the class's.
pub const GWAS_HAPLOTYPE_CROSS_PRODUCT: &str = "gwas_haplotype_cross_product";
/// `gwas_interaction_first_snp_allele_whitespace`: GWAS.pm:343-348 take each risk allele of a
/// `SNP x SNP` row as the split on `x` leaves it, so every allele but the last keeps the
/// space before the ` x ` and matches no variant; the candidate trims it. The shape: a
/// `GWAS_` field the candidate alone fills on a row at the location of such an SNP, the row's
/// allele that SNP's risk allele and the candidate's `GWAS_study` naming the catalog row's
/// study; an empty field of the candidate's annotated row (`GWAS_odds_ratio=`) is the class's
/// with its filled ones. Neither output says which catalog row annotated a site, so this class and
/// `gwas_haplotype_cross_product` read the catalog's rows of two or more SNPs
/// (`GwasCatalogPairs`, `--gwas-catalog`) and accept nothing without them: the cross-product
/// class a row at a location of such a row whose allele is another of its SNPs' risk alleles
/// and whose reference `GWAS_study` names its study.
pub const GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE: &str =
    "gwas_interaction_first_snp_allele_whitespace";
/// `alphamissense_two_digit_transcript_version`: AlphaMissense.pm:259 strips a prediction's
/// transcript version with `s/\.[0-9]//g`, one digit, so under `transcript_match=1` a
/// two-digit version (`ENST00000359568.10` read as `ENST0000035956810`) matches no row and
/// the reference writes no `am_` field there; the candidate strips the whole version. The
/// shape: an `am_` field (the JSON `alphamissense` object) the candidate alone fills.
pub const ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION: &str =
    "alphamissense_two_digit_transcript_version";
/// `cadd_structural_allele_match`: CADD.pm:222-223 gives a structural allele outside
/// `insertion`, `deletion` and `duplication` no allele string, and Sequence.pm:1182 and
/// :1196 reverse-complement a reverse-strand record's class word, so a structural row of an
/// Ensembl-format reverse-strand record matches no CADD-SV row in Ensembl VEP and every row of
/// the other classes raises `Plugin 'CADD' went wrong` (OutputFactory.pm:2253) in its
/// warnings; the candidate matches the row and raises nothing. The shape: a `CADD_` field the
/// candidate alone fills on a structural allele's row, and in --warnings mode the reference's
/// `Plugin 'CADD' went wrong` lines.
pub const CADD_STRUCTURAL_ALLELE_MATCH: &str = "cadd_structural_allele_match";
/// `splicevault_sibling_allele_rows`: SpliceVault.pm:156 matches a table row against every
/// alternate allele of the record (`alt_alleles`), so a row naming one allele of a
/// multi-allelic record fills the `SpliceVault_` fields of every allele's rows; the candidate
/// matches the row's own allele. The shape: a `SpliceVault_` field (the JSON `splicevault`
/// object) the reference fills with the value its row of another allele of the record carries
/// on the same feature, where the candidate leaves it empty or fills it with its own allele's
/// row; a field the candidate alone fills, and any value on a record of one allele, stays a
/// difference. `splicevault_json_later_rows`:
/// under `transcript_match=0` the plugin rewrites its cached row at :190 and dies at :183 on a
/// record's later transcript rows in the JSON format alone (`Plugin 'SpliceVault' went wrong:
/// Not an ARRAY reference`), so the reference's `splicevault` object is missing there; the
/// shape: the JSON `splicevault` object the candidate alone fills.
pub const SPLICEVAULT_SIBLING_ALLELE_ROWS: &str = "splicevault_sibling_allele_rows";
pub const SPLICEVAULT_JSON_LATER_ROWS: &str = "splicevault_json_later_rows";
/// `utrannotator_gained_start_window`: UTRAnnotator.pm:339-345 tests a gained start codon over
/// a window of the reference allele's length, so an insertion completing an `ATG` whose third
/// base lies past the window reports the lost codon alone, and a long deletion's window wraps
/// to the UTR's end through negative indices and reports a pre-existing `ATG` as gained; the
/// candidate sizes the window by the alternate allele. The window opens two bases before the
/// variant, so on a substitution at cDNA position 1 or 2 it wraps the same way. The shape
/// (`gained_start_window_row`): one engine alone reports the gained start, at a
/// `CapDistanceToStart` its own window holds and the other's does not, the candidate on a row
/// whose allele is longer than the reference allele and the reference on a row whose allele is
/// shorter; or the reference reports it at a negative `CapDistanceToStart`, which the wrapped
/// window alone yields. The row's `5UTR_annotation` goes with its `5UTR_consequence`. A
/// sequence allele at a two-base Location under an identifier is an insertion's or a two-base
/// reference allele's, so such a row stays a difference without the record's own alleles. A
/// start codon an insertion only moves on is none the inserted bases make: where the reference
/// reports a start codon lost and the candidate's gained start stands the allele's length past
/// that codon's index, the row is refused. The candidate names the first start codon of the
/// changed bases of the first way of writing an insertion's mutated sequence, so for a record
/// written further on it reports a codon that begins before the base after the insertion,
/// which the reference's window at the record's own index does not hold. On an insertion's
/// row the candidate alone may therefore report the gained start three or more bases before
/// the lower cDNA position, where the `Allele`, read on the transcript's strand (`STRAND`)
/// and on from its last base to its first, spells `ATG` from the base that stands as many
/// bases before its end;
/// a row that states no strand is refused there.
pub const UTRANNOTATOR_GAINED_START_WINDOW: &str = "utrannotator_gained_start_window";
/// `utrannotator_consequence_precedence`: when several of the plugin's five analyses fire on
/// one row, UTRAnnotator.pm:262-267 walks their flags as a hash and overwrites the output with
/// every truthy one, so `5UTR_consequence` and `5UTR_annotation` carry whichever effect the
/// process's hash order listed last; the candidate reports one under a fixed precedence. The
/// shape: a `5UTR_consequence` whose two values are both among the five `5_prime_UTR_`
/// terms (`:426`, `:562`, `:704`, `:838`, `:1003`) and differ, set aside with the
/// row's `5UTR_annotation`; a `-` on either side stays a difference. Two terms on a row do
/// not show that two analyses fire there, so the class reads Ensembl VEP's own runs of the
/// input under other Perl hash seeds (`UtrannotatorSeedSweep`) and accepts nothing without
/// them: both terms are ones those runs write on the row, and each side's `5UTR_annotation`
/// is one they write beside its term. The class follows
/// `utrannotator_gained_start_window` in `FIELD_CLASSES`, so a row that class explains is
/// counted under it alone.
pub const UTRANNOTATOR_CONSEQUENCE_PRECEDENCE: &str = "utrannotator_consequence_precedence";
/// `utrannotator_annotation_pair_order`: UTRAnnotator.pm:282-283 walks the reported effect's
/// groups and each group's `key=value` pairs as hashes, so the groups of `5UTR_annotation`
/// and the pairs inside each come in the process's hash order; the candidate writes each
/// group's pairs in byte order of their keys. The shape: a `5UTR_annotation` whose two values
/// carry the same groups of the same pairs in different orders, the candidate's pairs in key
/// order. The class precedes the plugin's other eight in `FIELD_CLASSES`, so a value differing
/// in order alone is counted under it on every row.
pub const UTRANNOTATOR_ANNOTATION_PAIR_ORDER: &str = "utrannotator_annotation_pair_order";
/// `utrannotator_structural_allele_text`: UTRAnnotator.pm:80-82 accept a structural variant,
/// whose allele string is the record's ALT text (Parser/VCF.pm:572), and `run` (:151-153)
/// reads that text as the reference allele beside an empty alternate, so the plugin removes as
/// many bases of the 5' UTR as the text has characters (:1158-1164, :1197, :1297) and reports
/// the effect of that deletion; the candidate annotates a symbolic deletion whose stated span
/// lies inside one exon of the 5' UTR as the sequence deletion of that span, a symbolic
/// inversion there as that span replaced by its reverse complement and the loss alternate of
/// a record with several alternates as the deletion of the record's span, gives a plain
/// symbolic insertion the frameshift its stated length decides and leaves every other
/// structural allele's two fields empty. Five shapes. Where the candidate fills neither
/// field (`structural_allele_text_row`): on a row whose Allele is a structural class word
/// the reference's `5UTR_consequence` is one of the plugin's five terms; the row's
/// `5UTR_annotation` goes with its `5UTR_consequence`. Where the candidate fills them
/// under a deletion's class word (`structural_deletion_span_row`): a default or tab row
/// whose Location's span is the distance between its two cDNA positions and whose
/// Consequence names `5_prime_UTR_variant` beside no term but `feature_truncation` and
/// `NMD_transcript_variant`, so the span lies inside one exon of the 5' UTR; the
/// candidate's `5UTR_consequence` is a term a deletion can fire beside a `5UTR_annotation`,
/// the frameshift only where the span is no multiple of three.
/// Where the candidate fills them under `inversion` (`structural_inversion_span_row`): a
/// row of the same span whose Consequence names `5_prime_UTR_variant` beside no term but
/// `NMD_transcript_variant`; the candidate's `5UTR_consequence` is any of the plugin's five
/// terms but the frameshift, which two alleles of one length do not fire (`:889`), beside a
/// `5UTR_annotation`; a record handed with the row states `<INV>` as its one alternate.
/// Where the candidate fills them under `copy_number_variation`
/// (`structural_copy_number_loss_row`): a row of the same span and Consequence as an
/// inversion's on which the candidate's `5UTR_consequence` is a term a deletion can fire
/// beside a `5UTR_annotation`, the frameshift only where the span is no multiple of three;
/// a record handed with the row states two or more alternates, and the row is that of one
/// whose own text is `<DEL>`, `<CN0>` or `<CN=0>`, the rows of a record on a feature
/// standing in its alternates' order.
/// Where the candidate fills them under `insertion` (`structural_insertion_length_row`):
/// a default or tab row whose Consequence names `5_prime_UTR_variant` and whose Location's
/// span, the stated length, is no multiple of three; the candidate's `5UTR_consequence` is
/// the frameshift, and each effect of its `5UTR_annotation` holds the six keys read from
/// the reference UTR, with or without `alt_type=uORF`, and no `alt_type_length`. Under
/// all four, the reference's `5UTR_consequence` is one of the plugin's five terms or both
/// of its fields are empty. A sequence allele's row, a breakend's own string, a row the
/// candidate fills under any other class word or outside those four shapes, a VCF or JSON
/// entry the candidate fills, an annotation beside no term and the plugin's `Existing_`
/// counts stay differences. The class follows the pair-order, the gained-start window and
/// the precedence classes in `FIELD_CLASSES`, so a row one of them explains is counted
/// under that one.
pub const UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT: &str = "utrannotator_structural_allele_text";
/// `utrannotator_lost_start_after_insertion`: UTRAnnotator.pm:764-769 pass an insertion whose
/// index lies at one of a start codon's three bases, and `:776-778` read the three bases at
/// the codon's own index of the mutated UTR and report the codon lost unless they spell `ATG`.
/// Where the insertion only moves the codon on, the codon stands whole, as many bases further
/// from the cap as the allele is long, and the candidate reports no lost start: the insertion
/// stands directly before the codon, or after its first or second base with inserted bases that
/// end with that base or those two, which is the same mutated sequence. The shape
/// (`lost_start_after_insertion_row`): the reference's `5UTR_consequence` is the lost-start
/// term and its `5UTR_annotation` names one `CapDistanceToStart`; the row's cDNA positions
/// are an insertion's two flanking bases, the lower of them the 5' UTR index of the base after
/// the insertion; that number is the reported index, or lies one above it on a row whose
/// `Allele`, read on the transcript's strand (`STRAND`), ends with `A`, or two above it on
/// a row whose allele so read ends with `AT`; the candidate's row carries the plugin's counts
/// and its `5UTR_consequence` is not the lost-start term. The row's `5UTR_annotation` goes
/// with its `5UTR_consequence`. A row that states no strand is read at the codon's own index
/// alone. A reference allele's row reports an index below its first position (`:764-769`), so
/// the equality is an insertion's alone; a row whose own record states a reference allele that
/// is not empty is refused all the same. The class follows the gained-start window, the
/// precedence and the pair-order classes in `FIELD_CLASSES`, so a row the gained-start window
/// or the seeded runs explain is counted under that class.
pub const UTRANNOTATOR_LOST_START_AFTER_INSERTION: &str = "utrannotator_lost_start_after_insertion";
/// `utrannotator_gained_start_remade_in_place`: UTRAnnotator.pm:345 opens the gained-start
/// window two bases before an insertion's index and `:346-353` report the first `ATG` in it,
/// whatever the reference holds. An insertion written after the first or second base of a
/// start codon whose inserted bases open with the codon's remaining bases, and do not also end
/// with the codon's bases before them, leaves an `ATG` at the codon's own index: the mutated
/// UTR is the reference with the inserted bases directly after the untouched codon, and that
/// codon is reported gained. The candidate reports no gained start for it and reads the
/// frameshift directly after the codon, so the codon's own ORF is among the frames it lists
/// unless a stop codon or the UTR's end follows the codon directly, where `:921` reads the
/// insertion inside the codon and leaves that ORF out. The shape
/// (`gained_start_remade_in_place_row`): the row's cDNA
/// positions are an insertion's two flanking bases; the reference's `5UTR_consequence` is the
/// gained-start term and its `5UTR_annotation` names one `CapDistanceToStart`, one below
/// the lower cDNA position on a row whose `Allele`, read on the transcript's strand
/// (`STRAND`), opens with `TG`, or two below it on a row whose allele so read opens with
/// `G`; the candidate's row carries the plugin's counts and reports no gained start at that
/// index. Where Ensembl VEP's hash order keeps another of its terms on the row, the gained
/// start is the one its runs under other hash seeds (`UtrannotatorSeedSweep`) write there,
/// at one index. The row's `5UTR_annotation` goes with its `5UTR_consequence`. A row that
/// states no strand, and a row whose own record states a reference allele that is not empty,
/// is refused. The class follows the plugin's other five in `FIELD_CLASSES`, so a row the
/// gained-start window or the seeded runs explain is counted under that class.
pub const UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE: &str =
    "utrannotator_gained_start_remade_in_place";
/// `utrannotator_frameshift_at_stop_boundary`: UTRAnnotator.pm:911-917 set an upstream ORF's
/// check point one index before its stop codon, and `:921` admits the ORF to the frameshift
/// only where the record's own end index is at or before the check point, reading no base of
/// the mutated UTR. An insertion written directly before the first base of a stop codon stands
/// inside the ORF that ends there: at a length that is no multiple of three the ORF is shifted
/// unless the inserted bases put a stop codon at that index; it then runs on and its length
/// changes.
/// `:921` leaves that ORF out (the comment at `:919-920` gives the boundary to the stop
/// codon's analysis), and `:631` skips the lost-stop test for an insertion at the stop
/// codon's first index, so no term reports the ORF; the candidate lists it. Conversely,
/// `:921` admits an ORF for every insertion written at or before the check point,
/// whatever the mutated UTR reads: where the inserted bases leave a stop codon at the stop
/// codon's own index, the ORF ends there with its length unchanged, and the candidate, which
/// reads the mutated sequence at every index that writes it, leaves the ORF out. Two shapes
/// (`frameshift_at_stop_boundary_row`), never both on one row. In each, an entry lists the
/// groups of its `5UTR_annotation` where its `5UTR_consequence` is the frameshift term and
/// none where it fills neither field; at least one entry reports the frameshift, and a
/// candidate's row that fills neither field carries the plugin's counts. The stop codon's
/// first base of a uORF (`ref_type=uORF`) is its `ref_StartDistanceToCDS` less its
/// `ref_type_length` plus three bases before the coding start. Where the candidate lists
/// more: the reference's groups are a proper subset of the candidate's; every group the
/// candidate adds is a uORF whose `alt_type_length` differs from its `ref_type_length` and
/// whose stop codon's first base is the base directly after the insertion as the record
/// writes it (`utr_base_after_insertion`: the row's HGVSc coordinates moved back by the
/// size of `HGVS_OFFSET`); and no group both list is such an ORF. Where the reference
/// lists more: the candidate's groups are a proper subset of the reference's; every group
/// the reference adds is a uORF that the mutated UTR closes too (`alt_type=uORF`) whose
/// `alt_type_length` equals its `ref_type_length` and whose stop codon's first base lies
/// at or after the base after the insertion as the record writes it, and no more than one
/// base beyond the base after the insertion as the HGVSc itself names it, the last index
/// that writes the sequence. A row whose HGVSc states no insertion point in the 5' UTR is
/// refused. The class sets aside `5UTR_annotation`, and `5UTR_consequence` with it where
/// one entry fills neither field, and follows the plugin's other six in `FIELD_CLASSES`,
/// so a row the seeded runs or a start codon left in place explain is counted under that
/// class.
pub const UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY: &str =
    "utrannotator_frameshift_at_stop_boundary";
/// `utrannotator_lost_stop_written_before_stop_codon`: UTRAnnotator.pm:628-633 read the three
/// bases at a stop codon's index of the mutated UTR only where the record's own index lies one
/// or two bases past the stop codon's first (`:631` skips an insertion at or before it). An
/// insertion written directly before a start codon whose second base opens the stop codon of
/// an upstream ORF, with inserted bases that open with the start codon's first two bases, is
/// the mutated sequence of the same bases written after the start codon's second base, inside
/// the stop codon: the stop codon's index reads no stop codon and the stop codon is lost.
/// Ensembl VEP reports the lost stop for the record written inside the stop codon and, for the
/// record written before the start codon, the ORF's frameshift (`:921`); the candidate
/// reports the lost stop on every way of writing the sequence. The candidate reads a stop
/// codon by every index that writes the mutated sequence: it reports a stop codon lost where
/// some writing lies inside it and the frame does not survive, and not lost where a writing
/// at its first base leaves it whole and in frame. Three shapes
/// (`lost_stop_written_before_stop_codon_row`), in each of which one entry alone reports
/// the lost-stop term, the row's cDNA positions are an insertion's two flanking bases, its
/// record, where it states its alleles, has an empty reference allele, and the candidate's
/// row carries the plugin's counts. Where the reference's `5UTR_consequence` is the
/// frameshift term and the candidate's the lost-stop term: each group of the candidate's
/// `5UTR_annotation` carries the
/// `KozakContext` and `KozakStrength` of exactly one uORF (`ref_type=uORF`) the
/// reference lists, a different one for each group, and the alternative stop codon where that
/// uORF's `alt_type_length` ends it: `AltStop=True` with `AltStopDistanceToCDS` equal to
/// the ORF's `ref_StartDistanceToCDS` plus the allele's length less its `alt_type_length`
/// plus three where its `alt_type` is `uORF`, `AltStop=False` with `NA` otherwise.
/// Where the candidate reports a lost stop the reference does not read (`:631`): the
/// reference fills neither field or its `5UTR_consequence` is the lost start, the gained
/// start or the frameshift, no group of a frameshift's list carrying the `KozakContext` and
/// `KozakStrength` of a group of the candidate's, unless that group's stop codon lies
/// outside the listed ORF's frame and whole codons past a stop codon the insertion can be
/// written inside (`lost_stop_of_another_orf`); the `Allele`'s length is no multiple of
/// three; the row's `HGVS_OFFSET` is a whole number other than zero, so the same sequence
/// can be written further on, inside the stop codon; and none of the seeded runs
/// (`UtrannotatorSeedSweep`), where they are given, writes the lost stop on the row. Where
/// the reference reports a lost stop the candidate does not (`:639-642`): the `Allele`'s
/// length is a multiple of three, and each group of the reference's `5UTR_annotation` reads
/// `AltStop=True` with an `AltStopDistanceToCDS` that is one or two bases more than the
/// distance of the base after the insertion as the record writes it
/// (`utr_base_after_insertion`), the ORF's own stop codon moved on by the allele's length,
/// or more by another number below the allele's length that is no multiple of three, a stop
/// codon the inserted bases put in frame before it; the candidate's `5UTR_consequence` is
/// another of the plugin's five terms or both of its fields are empty. The
/// row's `5UTR_annotation` goes with its `5UTR_consequence`. The class follows the
/// plugin's other classes but `utrannotator_without_reference` in `FIELD_CLASSES`, so a row
/// whose two terms the seeded runs write is counted under
/// `utrannotator_consequence_precedence`.
pub const UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON: &str =
    "utrannotator_lost_stop_written_before_stop_codon";
/// `utrannotator_without_reference`: UTRAnnotator.pm:88 prints a notice when the run has no
/// reference sequence and goes on, and `run` (:170, :173) then reads the coding and 5' UTR
/// sequences of a transcript whose exons hold none, so a row carrying `5_prime_UTR_variant`
/// gets `0` in the plugin's three `Existing_` counts (:297-298) and `-` in its two text
/// fields (:291-292), or the gained start where the allele's own bases spell a start codon
/// inside the window `uAUG_gained` reads (:345-354); the candidate refuses the plugin when
/// it loads, so its output names none of the five fields. The shape
/// (`without_reference_row`): the candidate's entry carries none of the five fields; the
/// reference's `5UTR_consequence` and `5UTR_annotation` are empty, or its
/// `5UTR_consequence` is the gained start; and the reference's three counts are each `0`,
/// or all empty on an entry of the tab or VCF format that the plugin left alone, which
/// reports no effect. An entry the candidate fills, any other reference term, an
/// annotation beside no term and a count other than `0` stay differences. The class
/// follows every class of `FIELD_CLASSES` but the two that read a row's term pair, so a
/// row another class explains is counted under that one.
pub const UTRANNOTATOR_WITHOUT_REFERENCE: &str = "utrannotator_without_reference";
/// `nearest_exon_jb_item_order`: NearestExonJB.pm:267-271, and `:245-263` under
/// `intronic=1`, build the value by walking the keys of the hash of the row's exons, so two
/// or more exon items come in the process's hash order; the candidate lists them in transcript
/// order. The shape: a `NearestExonJB` whose two values carry the same items the same number
/// of times in different orders.
pub const NEAREST_EXON_JB_ITEM_ORDER: &str = "nearest_exon_jb_item_order";
/// `loeuf_zero_minimum`: LOEUF.pm:204 drops a gene's minimum whose string is `0` by Perl
/// truthiness; the candidate prints it. The shape: a `LOEUF` the candidate alone fills, with
/// `0`.
pub const LOEUF_ZERO_MINIMUM: &str = "loeuf_zero_minimum";
/// `maxentscan_overlapping_windows`: a substitution in an exon of at most five bases lies in
/// the acceptor window of the intron before it and the donor window of the intron after it;
/// MaxEntScan.pm:270 scores the window of whichever intron Ensembl VEP's interval tree returns
/// first, the candidate the first window in transcript order, the upstream intron's acceptor. A
/// substitution in both windows of one intron, an intron of at most 25 bases, is scored at the
/// donor by both. The shape: on a row whose record is one base in the acceptor window of one
/// intron of the row's transcript and the donor window of the next, read from the candidate's
/// JSON cache (`SpliceWindows`), a `MaxEntScan_` field both engines fill with different
/// numbers, and a window text of a `verbose` run where the reference writes a donor's nine
/// bases and the candidate an acceptor's twenty-three. Without the cache nothing is accepted.
pub const MAXENTSCAN_OVERLAPPING_WINDOWS: &str = "maxentscan_overlapping_windows";
/// `alphamissense_sibling_allele_rows`: AlphaMissense.pm:245 matches a table row against every
/// alternate allele of the record and :250-253 keep the first row whose amino-acid letters equal
/// the output row's, so two alleles of one record giving one amino-acid change both take the
/// table row that comes first; the candidate matches each allele's own row. The shape: an
/// `am_` field (the JSON `alphamissense` object) the reference fills with the value its row
/// of another allele of the record carries on the same feature beside the same `Amino_acids`,
/// where the candidate leaves it empty or fills it differently.
pub const ALPHAMISSENSE_SIBLING_ALLELE_ROWS: &str = "alphamissense_sibling_allele_rows";
/// `revel_sibling_allele_rows`: REVEL.pm:165 matches a score row against every alternate
/// allele of the record and :176-199 return the first row whose alternate amino acid equals
/// the output row's peptide, the widening AlphaMissense's match has. The shape: a `REVEL`
/// the reference fills with the value its row of another allele of the record carries on the
/// same feature beside the same `Amino_acids`, where the candidate leaves it empty or fills
/// it differently.
pub const REVEL_SIBLING_ALLELE_ROWS: &str = "revel_sibling_allele_rows";
/// `revel_other_reference_amino_acid_rows`: a REVEL table row scores one substitution, from
/// its `aaref` to its `aaalt`, on the transcripts it lists. REVEL.pm:210 and :229 read the
/// `aaref` cell and :193 and :197 test the `aaalt` alone, so the reference writes the score
/// of a row whose reference amino acid is not the output row's; the candidate tests both
/// cells. The shape: a `REVEL` the reference fills with the score of the table row its walk
/// reaches (REVEL.pm:170-199), a row whose `aaref` is not the first residue of the row's
/// `Amino_acids`, where the candidate writes what its own walk gives: nothing, or the score
/// of the first row of the row's own alternate that names both residues. The class reads the
/// table the run's plugin read (`RevelTable`, `--revel-table`) and accepts nothing
/// without it.
pub const REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS: &str = "revel_other_reference_amino_acid_rows";
/// `gnomadmt_info_value_second_equals`: gnomADMt.pm:327-330 reads each INFO item of a file
/// record with `my ($key, $value) = split /=/, $field`, which Perl limits to three parts, so
/// a value holding `=` is cut there; the candidate splits the item at its first `=` alone.
/// The shape: a `gnomAD_` field both fill, the reference's value the candidate's up to its
/// first `=`, read with the percent-encoding of `=` undone.
pub const GNOMADMT_INFO_VALUE_SECOND_EQUALS: &str = "gnomadmt_info_value_second_equals";
/// `loftee_anc_allele_multiallelic_record`: LoF.pm:531 returns before the ancestral comparison
/// when the record's allele string holds `-` or is not three characters long, so no allele of
/// a multi-allelic record is filtered `ANC_ALLELE`; the candidate judges each one-base
/// substitution by its own base. The shape: on a one-base substitution's stop-gained or
/// splice-site row of a record that has a row of another allele on the same feature, a
/// `LoF_filter` the candidate
/// writes as the reference's plus `ANC_ALLELE`, set aside with the row's `LoF`, `LC`
/// beside the item and `HC` on a side with no item.
pub const LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD: &str = "loftee_anc_allele_multiallelic_record";
/// `loftee_first_alt_donor_filter`: LoF.pm:363-367 test the reference and the first alternate
/// allele of the record for every allele's row, so a later allele of a multi-allelic donor
/// record is filtered `GC_TO_GT_DONOR`, or spared, by its sibling's base; the candidate tests
/// each allele's own base. The shape: on a `splice_donor_variant` row of a record that has a
/// row of another allele on the same feature, and that is not the row of the record's first
/// alternate, a `LoF_filter` whose two values differ by `GC_TO_GT_DONOR` alone, either
/// side carrying it, set aside with the row's `LoF` as the ancestral class's is. The first
/// alternate is read where an entry states the record's alleles in order (`first_alternate`);
/// where none does the class accepts nothing.
pub const LOFTEE_FIRST_ALT_DONOR_FILTER: &str = "loftee_first_alt_donor_filter";
/// `loftee_gerp_zoom_level_summary`: gerp_dist.pl:70-79 total each exon interval's GERP score
/// through one summary bin of the bigWig, which the library serves from a zoom level once the
/// interval is twice that level's span, 80 bases on the file LoFTEE's GRCh38 bundle ships, so
/// `GERP_DIST` approximates the per-base sum the candidate writes, and LoF.pm:254 pushes
/// `END_TRUNC` by that total. The shape: on a stop-gained or frameshift row whose `BP_DIST`
/// is at least 79 (an interval of N bases walks N - 1, so 79 is the 80-base interval), a
/// `LoF_info` whose two values differ in the `GERP_DIST` item alone, both
/// numbers; and, where the row's `DIST_FROM_LAST_EXON` is at most 50 and the two totals
/// straddle the cutoff -58, a `LoF_filter` differing by `END_TRUNC` alone on the side of
/// the lower total, set aside with the row's `LoF` as the ancestral and donor classes' is. The
/// cutoff is the plugin's default: a run that sets `gerp_end_trunc_cutoff` is not read.
pub const LOFTEE_GERP_ZOOM_LEVEL_SUMMARY: &str = "loftee_gerp_zoom_level_summary";
/// `loftee_gerp_walk_without_conservation_database`: on the `master` branch of the LOFTEE
/// repository, the plugin of a GRCh37 run, gerp_dist.pl:59-60 read the GERP total of each
/// whole exon the walk crosses, an exon after the variant's that does not hold the stop codon,
/// from the conservation database, and :129 call `prepare` on its handle. LoF.pm:119-120
/// leave the text `false` in the handle's place in a run without `conservation_file`, and
/// LoF.pm:258 hands it to the walk of a stop-gained or frameshift row, so the plugin dies on
/// a row whose walk crosses a whole exon and Ensembl VEP writes the row without the plugin's
/// fields; the candidate totals the exon's bases from the per-base GERP file and writes
/// them. The shape: on a row whose Consequence holds `stop_gained` or `frameshift_variant`
/// on both sides, where the reference's row carries none of `LoF`, `LoF_filter`,
/// `LoF_flags` and `LoF_info` and the candidate's row carries `LoF` and a `LoF_info` that
/// holds a `GERP_DIST` item, each of the four fields the candidate fills. Under
/// `skip_lof_info` LoF.pm:243 and :339 write no `LoF_info` on any row, so the candidate's
/// row of such a run carries `LoF` and no `LoF_info`, and the class takes that row as well.
/// A row on which the reference writes one of the four, a row without either term, a row on
/// which the candidate writes no `LoF` and a row whose candidate `LoF_info` holds items and
/// no `GERP_DIST` among them stay differences.
pub const LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE: &str =
    "loftee_gerp_walk_without_conservation_database";
/// `loftee_nagnag_site_first_transcript_strand`: LoF.pm:489-494 read the nine reference bases
/// that decide the `NAGNAG_SITE` item of `LoF_flags` once for an input record, on the strand
/// of the first transcript whose `splice_acceptor_variant` row reaches the test
/// (LoF.pm:306-307), and keep them for the record's other transcripts; the candidate reads
/// them on each row's own strand. Where transcripts of both strands share an acceptor base,
/// a row can carry the item on one side alone. The shape
/// (`nagnag_site_first_transcript_strand_row`): a row whose terms include
/// `splice_acceptor_variant` on both sides; whose two `LoF_flags` hold the same items in the
/// same order but for one `NAGNAG_SITE`, either side carrying it, an empty or absent value
/// holding no item; whose `LoF`, `LoF_filter` and `LoF_info` are each equal on both sides;
/// and whose record has, among the reference's rows, a `splice_acceptor_variant` row on the
/// other strand, each strand read from the row's `STRAND`. A list that differs in another
/// item or in its order, a row without the term, a row whose `LoF`, `LoF_filter` or
/// `LoF_info` differs, a row that states no strand and a record with no such row on the
/// other strand stay differences, and the class sets aside no other field of the row. The
/// record's rows are read as `loftee_first_alt_donor_filter` reads them, a default or tab
/// row under an identifier from the input's lines, and without them the class accepts
/// nothing. The class follows every class of `FIELD_CLASSES` but the two that read a row's
/// term pair.
pub const LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND: &str =
    "loftee_nagnag_site_first_transcript_strand";
/// `nearest_exon_jb_stretch_priming`: NearestExonJB.pm:184 reads `_overlapped_exons`, whose
/// 12-base stretch (BaseTranscriptVariation.pm:858-865) applies only once a record inside the
/// transcript has set its frameshift-intron flag (:1032-1034), so a record within 12 bases
/// outside the first or last exon gets a value or none by the records before it; the candidate
/// decides the stretch from the transcript's introns. The shape: a `NearestExonJB` the
/// candidate alone fills on a row whose `DISTANCE` is at most 12.
pub const NEAREST_EXON_JB_STRETCH_PRIMING: &str = "nearest_exon_jb_stretch_priming";
/// `nearest_exon_jb_shifted_span`: under `--shift_3prime 1` the consequence tests fill the
/// `_overlapped_exons` memo at the shifted span (BaseVariationFeatureOverlapAllele.pm:454-463,
/// :514) and NearestExonJB.pm:184 reads that memo while :225-226 measure from the record's
/// span; the candidate reads the exons at the record's span. The shape: a `NearestExonJB` on
/// a transcript row with a positive `SHIFT_LENGTH` that one side fills and the other leaves
/// empty, or whose two values name different exons.
pub const NEAREST_EXON_JB_SHIFTED_SPAN: &str = "nearest_exon_jb_shifted_span";
/// `nearest_exon_jb_id_less_exon_fold`: NearestExonJB.pm:220-237 file each exon's distance,
/// side and length under its stable id, so the exons of a row that carry none share the empty
/// key and :267-271 (:244-263 under `intronic=1`) print one item from their parts, the
/// nearest exon's distance and side with the length filed last; the candidate prints each
/// exon's own item. The shape: a `NearestExonJB` both fill with items of an empty exon id,
/// the reference's one item whose distance and side are those of one of the candidate's items,
/// on a row the candidate gives two or more items, or one whose Consequence carries both
/// splice-site terms and whose cDNA position names two positions: a variant that begins in one
/// exon and ends in another. A deletion over one whole exon carries both terms and no such
/// position, and overlaps that exon alone.
pub const NEAREST_EXON_JB_ID_LESS_EXON_FOLD: &str = "nearest_exon_jb_id_less_exon_fold";
/// `nearest_start_site_one_base_above`: TranscriptTree.pm:165 at release/116.3 stores a
/// protein-coding start site at `t` as the half-open interval from `t - 1` to `t`, and
/// :233-236 take the nearest site at or below a record's end and the nearest above its
/// start, asking `fetch_nearest_up` with the 1-based start. The library returns only an
/// interval that begins above its query, so the search above the start returns neither a
/// site on the start nor a site one base above it. Such a site is named only where the
/// search at or below the end returns it; on any other record where the two searches
/// return a site but none nearer, a farther site, or one of two sites at one distance, is
/// named in its place. The candidate searches from the start as well. The shape, with
/// `--nearest-sites`: a `NEAREST` on a row whose Location has a site on its first base or
/// one base above it, where the reference's two searches keep one set of sites and the
/// candidate's three another, the reference's values name the first set and the
/// candidate's the second. A row with no such site, a value that names another site and a
/// pair whose two sets are equal stay differences.
pub const NEAREST_START_SITE_ONE_BASE_ABOVE: &str = "nearest_start_site_one_base_above";
/// `downstream_reverse_strand_rotation`: TranscriptVariationAllele.pm:1330-1347 rotate a
/// 3'-shifted allele `seq_length - shift_length` places on a reverse-strand transcript, none
/// once the shift exceeds the allele, and Downstream.pm:110 splices the bases as they stand;
/// the candidate rotates by the shift modulo the allele's length. The shape: a
/// `DownstreamProtein` or `ProteinLengthChange` both fill differently on a `frameshift_variant` row
/// of a `STRAND` -1 transcript whose `HGVS_OFFSET` exceeds the number of bases its HGVSc
/// inserts, duplicates or deletes and is not a multiple of it.
pub const DOWNSTREAM_REVERSE_STRAND_ROTATION: &str = "downstream_reverse_strand_rotation";
/// `spliceai_chopped_multiallelic`: for a record whose allele string holds `-`,
/// SpliceAI.pm:259-266 rebuild every allele with the anchor base before matching, so a
/// substitution beside a deletion is compared as a two-base pair, which equals no single-base
/// row (:272); the candidate compares each allele at its own position. The shape: a `SpliceAI_`
/// field (the JSON `spliceai` object) the candidate alone fills on a substitution allele's
/// row of a record whose rows carry an allele written `-`.
pub const SPLICEAI_CHOPPED_MULTIALLELIC: &str = "spliceai_chopped_multiallelic";
/// `spliceai_multi_entry_cell`: SpliceAI.pm:326-332 remove `SpliceAI=` from a score record's
/// INFO column and split the rest on `|` into the allele, the first symbol and everything
/// after it, and :305 file the record under that symbol alone. A record whose key holds several
/// comma-joined entries so gives the rows of the first entry's symbol the whole remainder as
/// their value (the nine fields of the other formats cut from it, :276-287) and the rows of a
/// later entry's symbol nothing (:314-316); the candidate writes a row every entry of the
/// record that names the row's gene, in the record's order. Where a record line stands more
/// than once at an allele, as in the SpliceAI golden
/// corpora, :305 keep one line a first symbol, the last, so the reference writes a row the
/// whole value of the last line its gene leads, or nothing where it leads none, and the
/// candidate the gene's entries of the last line that names the gene. A cutoff flag is the
/// reference's verdict on the first entry of its line (:336-338) and the candidate's on any
/// of the gene's entries. The shape is read from the score
/// files themselves (`SpliceAiTable`,
/// `--spliceai-table`) and nothing is accepted without them. A row's gene is its `SYMBOL`
/// where the row writes one, else its transcript's symbol in the candidate JSON cache
/// (`TranscriptSymbols`, `--json-cache`); a row whose gene neither names is not accepted.
pub const SPLICEAI_MULTI_ENTRY_CELL: &str = "spliceai_multi_entry_cell";
/// `plugin_contig_mitochondrial_spelling`: BaseVepTabixPlugin.pm:693-708 map a record's
/// chromosome name to a plugin's file by adding or removing `chr` alone, so a record spelt
/// `MT`, `chrM` or `chrMT` finds nothing in a file naming the contig otherwise; the
/// candidate queries the file under all four spellings where they name one sequence, on any
/// assembly but GRCh37, and on GRCh37, where `MT` and `chrMT` name another sequence than
/// `chrM` and `M`, under the other pair than the record's only where the file's rows state
/// the bases the run's reference sequence holds at their positions. The shape, on every
/// assembly: a field of a tabix-backed plugin other than gnomADMt that the candidate alone
/// fills on a row of a record so spelt. The comparison reads the two outputs, so it tells
/// the record's spelling and which output fills the field, and neither the name the file
/// gives the contig nor whether its rows state the reference bases. The fields are those
/// the plugins name themselves; a column a dbNSFP run selects and a TranscriptAnnotator
/// table's fields carry names of the user's and are not read.
pub const PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING: &str = "plugin_contig_mitochondrial_spelling";
/// `tabix_allele_case_insensitive_match`: `get_matched_variant_alleles` compares the
/// `ref_alt_pos` keys of a table row and of a record as text (Utils/Sequence.pm:1207, 1232)
/// and Parser.pm:554 has upper-cased the record's alleles, so a table row whose allele cells
/// are written in lower case matches no record; SpliceAI.pm:272 compares its row's alleles
/// with `eq` to the same effect. The candidate reads such a row in either letter case. The
/// shape: a field of the AlphaMissense, AVI, CADD, gnomADMt, GWAS, OpenTargets, REVEL,
/// SpliceAI or SpliceVault plugin that the candidate alone fills, and a `GWAS_` field the
/// candidate alone writes as the empty string, as the plugin writes a cell its catalog row
/// leaves empty, on a row whose other `GWAS_` fields the candidate alone fills. A row holds
/// no sign of its table's letter case, so the class is one a case names for a table whose
/// row it wrote in lower case. A gnomADc coverage field and a TranscriptAnnotator table's
/// fields, whose names are the user's, are not read.
pub const TABIX_ALLELE_CASE_INSENSITIVE_MATCH: &str = "tabix_allele_case_insensitive_match";
/// `dbnsfp_mitochondrial_record_named_m`: dbNSFP.pm:313 names a record on `MT` or `chrMT`
/// `M` to the plugin's table and :335-338 match a row on its position and alternate allele,
/// so on GRCh37, where `M` and `chrM` name another sequence than `MT` and `chrMT`, the
/// reference writes the value of the row a table naming the contig `M` or `chrM` holds at
/// the record's position, whatever sequence the table counts on; the candidate reads such a
/// table for such a record only where its rows state the bases the run's reference sequence
/// holds at their positions. The shape: beside a manifest whose assembly is GRCh37
/// (`GRCH37_NAMES`, the names `Assembly::parse` reads as that assembly), a dbNSFP field the
/// reference alone fills on a row of a record spelt `MT` or `chrMT`. A field is dbNSFP's
/// where its name carries the plugin's own prefix `dbNSFP_` (dbNSFP.pm:159, :257-261) or is
/// a column the manifest's run hands the plugin without `output_prefix` (`dbnsfp_columns`).
/// The comparison reads the two outputs and the manifest, so it tells neither the columns of
/// a run that selects `ALL`, which the table's header alone names, nor the name the table
/// gives the contig, nor the sequence its rows count on, nor whether the run gave the
/// candidate a reference sequence, without which it reads no table under the other pair of
/// names than the record's; beside no manifest, or one that names no assembly, the class
/// accepts nothing.
pub const DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M: &str = "dbnsfp_mitochondrial_record_named_m";
/// `dosage_sensitivity_regulatory_row_warnings`: DosageSensitivity.pm's `run` reads
/// `$tva->transcript` on a regulatory or motif row under --regulatory and dies inside
/// `run_plugins`' eval, one `Plugin 'DosageSensitivity' went wrong` warning per row
/// (OutputFactory.pm:2253); the candidate writes nothing on those rows and no warning. Applied
/// in --warnings mode alone: the reference's lines are set aside.
pub const DOSAGE_SENSITIVITY_REGULATORY_ROW_WARNINGS: &str =
    "dosage_sensitivity_regulatory_row_warnings";
/// `fork_warning_pid_prefix`: under --fork, Ensembl VEP's parent relays a child's captured
/// standard error as one warning headed by the child's process id (Runner.pm:565), and
/// `warning_msg` (BaseVEP.pm:806) puts `WARNING: ` before it, so the block's first line reads
/// `WARNING: <pid> : WARNING: <text>` where the serial run and the candidate write
/// `WARNING: <text>`. Applied in --warnings mode alone: the prefix is removed from the
/// reference's lines before they are parsed, each line changed counted under the class.
pub const FORK_WARNING_PID_PREFIX: &str = "fork_warning_pid_prefix";
/// `json_mismatch_warning_names_no_record`: in a --json run on a RefSeq or merged cache
/// Ensembl VEP builds each row on an empty hash (OutputFactory/JSON.pm:320 at release/116.3),
/// so the transcript-assembly mismatch warning joins a name the hash lacks
/// (OutputFactory.pm:1727) and reads `WARNING: Transcript-assembly mismatch in ` with no
/// record name, written once a process (BaseVEP.pm:809-813); the candidate writes the named
/// line for each record that draws the warning. Applied in --warnings mode alone, and only
/// to that shape: the reference holds the unnamed line and no named one, and the candidate a
/// named one. The reference's unnamed lines and the candidate's named lines are set aside,
/// with each line of the reference that is the interpreter's `Use of uninitialized value`
/// warning of the statement that wrote an unnamed line: a forked run relays that warning
/// into the file, and the statement is read from the line directly before an unnamed line.
/// Each line set aside is counted under the class. A relayed line that keeps its `<pid> : `
/// prefix is left to `fork_warning_pid_prefix`.
pub const JSON_MISMATCH_WARNING_NAMES_NO_RECORD: &str = "json_mismatch_warning_names_no_record";
/// `refseq_edited_model_hgvsc_offset`: on a transcript model whose alignment to the genome
/// inserted or deleted bases (`BAM_EDIT=OK`), Ensembl VEP 116.3 numbers the HGVSc from the
/// unedited exons: the cache dumper drops the model's edit flag and mapper (Dumper/Core.pm:164,
/// :169), so `_get_cDNA_position` rebuilds an unfolded one, and an intronic allele counts no
/// edit (TranscriptVariationAllele.pm:1608-1610); the candidate numbers every position against
/// the edited model. The shape: an HGVSc on a `BAM_EDIT` `OK` row whose two values name
/// the same reference and differ in their coordinate tokens alone, each pair of coordinates
/// covering as many positions in one value as in the other. One pair of other lengths is the
/// class's too, a deletion across bases the model lacks, whose end the reference numbers
/// through them (`deletion_across_lacked_bases`): on a row whose STRAND is 1, two lone
/// `del` over coding positions that open at one position, the reference's covering as many
/// positions as GIVEN_REF has bases and the candidate's as many as the row's cDNA
/// position, the fewer. Any other two values that cover different numbers of positions,
/// those on bases the model holds and the genome lacks among them, are refused and stay a
/// difference, each value printed as its engine wrote it. The count is read where a pair
/// states it, from two plain positions or two offsets off one exon boundary; a pair that
/// mixes the two, or runs over the stop codon, is not measured.
pub const REFSEQ_EDITED_MODEL_HGVSC_OFFSET: &str = "refseq_edited_model_hgvsc_offset";
/// `refseq_edited_model_hgvs_offset`: on an edited `NM_` or `XM_` model, `_return_3prime`
/// (TranscriptVariationAllele.pm:164-187) copies the first shift hash attached to the variant,
/// whatever transcript or strand made it, because its length test at :182-183 fails for every
/// insertion and deletion; a transcript on the other strand inherits that transcript's shift, a
/// zero that OutputFactory.pm:1761-1764 does not print, while the 3'-first duplication test of
/// Utils/Sequence.pm:578 still names the shifted notation; the candidate shifts each row on its
/// own transcript. The shape: an HGVS_OFFSET on a coding `BAM_EDIT` `OK` row whose HGVSc
/// names a duplication, the reference leaving the offset empty and the candidate filling it with
/// a non-zero integer, the distance of a 3' shift, negative on the reverse strand. The branch is
/// entered only for an allele with CDS coordinates (TranscriptVariationAllele.pm:147-148), and
/// only a duplication keeps its notation whatever the shift: a deletion or a plain insertion
/// that inherits the foreign shift changes its HGVSc too, another shape.
pub const REFSEQ_EDITED_MODEL_HGVS_OFFSET: &str = "refseq_edited_model_hgvs_offset";
/// `refseq_edited_model_shifted_insertion_bases`: on a model carrying an `_rna_edit`
/// (`BAM_EDIT=OK`), `hgvs_transcript` respells the notation's alleles from `feature_seq`
/// (TranscriptVariationAllele.pm:1510-1513) after `hgvs_variant_notation` described an insertion
/// at its 3'-shifted coordinates with the rotated allele (:1444, :1474-1483). `feature_seq` holds
/// the rotated allele only on a coding row under --shift_3prime 1, where `codon` (:813-820) has
/// run `shift_feature_seqs` (:1338-1342) before the notation is built; `hgvs_protein` (:1778)
/// rotates it too, after `hgvs_transcript` has read it (OutputFactory.pm:1760, :1769). On every
/// other row the notation spells the inserted bases as written at the shifted position, a sequence
/// the variant does not produce. The candidate spells the rotated bases on every row. The shape:
/// an HGVSc on a `BAM_EDIT` `OK` insertion row whose HGVS_OFFSET is one non-zero integer on
/// both sides, the two notations insertions of one transcript, numbering and coordinate pair, the
/// reference's bases the row's Allele read on the transcript's strand (reverse-complemented where
/// STRAND is -1; a row without STRAND admits either reading) and the candidate's those bases turned
/// towards the transcript's 3' end by the offset modulo their length; the same bases, another
/// length, a duplication and a row without the offset stay differences.
pub const REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES: &str =
    "refseq_edited_model_shifted_insertion_bases";
/// `refseq_edited_model_seam_insertion_hgvsc`: an insertion standing on the bases an edited
/// model holds and the genome lacks. On a RefSeq
/// model whose alignment to the genome inserted bases (an `_rna_edit`; the row carries
/// `BAM_EDIT` `OK` or `FAILED`) the cached mapper numbers the two genomic bases around
/// those bases more than one apart, so an insertion between the two has the model's bases between
/// them as its span: a `CDS_position` of `75-77` for three bases and `1059` for one, where
/// an ordinary insertion prints two neighbouring positions. `hgvs_transcript` slides the
/// insertion's two transcript positions by the misalignment offset
/// (TranscriptVariationAllele.pm:1517, :1524-1525; `get_misalignment_offset`, :1594-1624) and
/// writes `ins` of the allele between them, at plain coordinates (`c.77_78`) or with intron
/// offsets (`c.2034+2_2034+3`), naming neither the model's bases nor their span. The candidate
/// names the replacement of those bases by the allele, typed as `hgvs_variant_notation` types
/// an allele replacing reference bases (Utils/Sequence.pm:493), clipped as `_clip_alleles`
/// clips one (TranscriptVariationAllele.pm:1492, :2188) and placed by the HGVS 3' rule on the
/// model's own sequence, whose bases beyond the row's codons the row does not carry, so the
/// notation is judged by what it says on those codons. The shape: an HGVSc on an insertion's row
/// carrying `BAM_EDIT` `OK` or `FAILED` on both sides, the reference's notation an `ins`
/// of the row's Allele read on the transcript's strand (reverse-complemented where STRAND is -1;
/// a row without STRAND admits either reading) at any coordinate pair of the `c.` numbering,
/// the candidate's row with a CDS position of one position or a span inside the codons its
/// Protein_position and Codons give, the alternate codons being the reference codons with the
/// span's bases replaced by the allele, and the candidate's notation none where the allele is
/// those bases and otherwise one that, applied to the reference codons, gives the alternate
/// codons: under the reference's sequence name and at plain coding positions, a substitution of
/// the base it names, a `del`, a `delins`, a `dup`, an `inv` or a repeat count of
/// positions inside the codons, or an `ins` between two neighbouring positions, of which one
/// may lie one base outside them. An ordinary insertion's alternate codons are longer than its
/// reference codons by the whole allele, which tells it from a two-base span printing the same
/// two positions. A notation that gives another sequence or names a position outside the codons,
/// an `ins` between positions that are not neighbours, any notation where the allele is the
/// model's bases, a row without `BAM_EDIT` (a mitochondrial model's), a row without a CDS
/// position or Codons (a span in an untranslated region or on a non-coding model, whose row
/// does not carry the model's bases), a reference notation that is not an `ins` of the allele
/// and another reference sequence stay differences
/// (`edited_model_seam_insertion_replacement`).
pub const REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC: &str =
    "refseq_edited_model_seam_insertion_hgvsc";
/// `refseq_edited_model_unminimised_allele_hgvsc`: the HGVSc of an allele that keeps a base it
/// shares with its reference on a model whose `_rna_edit` the alignment applied; the comment
/// above `GIVEN_REF_FIELDS` states the shape `edited_model_clipped_alleles` reads.
pub const REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC: &str =
    "refseq_edited_model_unminimised_allele_hgvsc";
/// `hgvs_shifted_insertion_clamped_end`: under --shift_3prime 1 with --hgvs, an insertion
/// anchored on the base before a transcript's first base and shifted k bases into exon 1 maps to
/// the transcript slice coordinates (1, 0) in `_var2transcript_slice_coords`
/// (TranscriptVariationAllele.pm:2689-2736, either strand at :2704-2714), which :2729 keeps,
/// since only the end is below 1, and :2733 clamps the end to 1; `hgvs_transcript` adds k to
/// both ends (:1474-1483), so `hgvs_variant_notation` reads the one base at n.k+1 against the
/// inserted bases and the row prints a substitution or a one-base delins there (`n.15A>T`,
/// `n.15delinsTT`), or, where the inserted sequence ends with that base, `_clip_alleles`
/// (:2231-2242, :2271-2278) strips it and prints an insertion of the other bases at (n.k, n.k+1).
/// The candidate writes the notation of the shifted span, as Ensembl VEP does for the allele
/// written at the shifted position: a duplication ending at n.k or an insertion between n.k and
/// n.k+1. The shape: an HGVSc differing on an insertion's row whose SHIFT_LENGTH k is positive,
/// whose cDNA position starts at k (the shift carried the insertion from before the first base to
/// after base k) and whose HGVS_OFFSET is set, the reference's notation at one position P
/// (`P X>Y`, `P delins`, or `P-1_P ins`) and the candidate's a `dup` ending at P-1 or an
/// `ins` spanning P-1 and P, both on the same transcript and numbering, the candidate's
/// duplicated or inserted bases as many as the allele's, the reference's as many less the one it
/// stripped, and a candidate `ins` carrying the reference's bases; an empty HGVSc on either side
/// stays a difference. On a model that carries an `_rna_edit`, `hgvs_transcript` then writes
/// both alleles over from the row (:1508-1514), the reference allele an insertion's `-` and
/// the alternate the row's whole allele: the reference's substitution reads `P -><base>`, and
/// its `ins` holds the allele's bases whole, none stripped. The class takes `P -><base>` as it
/// takes `P X>Y`, a form no other notation has, and on a row carrying `BAM_EDIT` `OK` a
/// reference `ins` between P-1 and P whose bases are the row's Allele on the transcript's
/// strand against a candidate `dup` of as many bases ending at P-1. Without `BAM_EDIT` `OK`
/// a reference `ins` of the whole allele is an ordinary shifted insertion's notation and stays
/// a difference, so the same `ins` on a mitochondrial model, which carries an edit and prints
/// no `BAM_EDIT`, is not read.
pub const HGVS_SHIFTED_INSERTION_CLAMPED_END: &str = "hgvs_shifted_insertion_clamped_end";
/// `hgvs_shifted_indel_written_span_gate`: `_var2transcript_slice_coords`
/// (TranscriptVariationAllele.pm:2708-2713) maps the span as written onto the transcript slice and
/// refuses one wholly outside it (:2729), and `hgvs_transcript` refuses a clamped end the shift's
/// offset carries past the last base (:1472) before adding the offset (:1477-1478), so an insertion
/// or deletion the shift carries into a transcript from wholly before it, or over its last base,
/// has no HGVSc and no HGVS_OFFSET on a row annotated inside the transcript. The candidate
/// describes the shifted span, as Ensembl VEP does for the allele written there. The shape: an
/// empty reference HGVSc against a candidate insertion or deletion notation of the row's
/// transcript, an empty reference HGVS_OFFSET against a candidate non-zero integer, and an empty
/// reference HGVSIntronOffset plugin field against a candidate integer, on an indel's row annotated
/// inside the transcript (no upstream or downstream term; one whose allele or Location says so, or,
/// for an allele of two or more bases at a Location of as many, which the text formats print alike
/// for an insertion and a substitution, one whose candidate notation is a `dup` or `ins`) whose
/// SHIFT_LENGTH, where the run prints one, is positive; the candidate row carries the notation and
/// the offset together, so each pair is read beside the row's other fields. The plugin's two
/// fields are the intron offsets of the notation's coordinates, so beside a `del` notation,
/// whose printed coordinates are the ones the plugin reads, they are held to them. On an edited
/// model's row (`BAM_EDIT` `OK` on both sides) the return at :1472 comes before
/// `get_misalignment_offset` (:1517) stores the model's offset (:1621), so the reference prints
/// no REFSEQ_OFFSET either (OutputFactory.pm:1778): an empty reference REFSEQ_OFFSET against a
/// candidate non-zero integer is the class's on a row whose HGVSc and HGVS_OFFSET it reads. The
/// row does not say that the value is the sum of the model's edits before the allele, which
/// the cache holds. A REFSEQ_OFFSET the reference alone writes, two different values, a value
/// on a row whose notation the class does not read and a row without `BAM_EDIT` `OK` on a
/// side stay differences.
pub const HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE: &str = "hgvs_shifted_indel_written_span_gate";
/// `hgvs_intron_offset_abutting_insertion`: HGVSIntronOffset.pm:62 asks for the HGVS notation
/// on every Transcript row, so on an insertion directly beside a transcript's first or last
/// base the reference fills the plugin's two fields from a notation the row never prints, and
/// the candidate writes neither where the HGVS 3' rule places no notation on the transcript.
/// The shape and the three pairs it takes are the predicate's
/// (`intron_offset_abutting_insertion`).
pub const HGVS_INTRON_OFFSET_ABUTTING_INSERTION: &str = "hgvs_intron_offset_abutting_insertion";
/// `hgvsp_shifted_indel_written_span_gate`: `hgvs_protein` shifts the allele 3' and reads its
/// translation coordinates at the shift (TranscriptVariationAllele.pm:1711, 1732-1733), and tests
/// the `coding` pre-consequence predicate of the allele as written (:1731), computed and cached
/// with the row's terms (BaseVariationFeatureOverlapAllele.pm:340, 394, 460-462), where an allele
/// written in an intron inside the coding region reads `non_coding` (:535-537) and one written
/// on the transcript before its coding region, on exon bases of the 5' UTR or in an intron of
/// it, reads `utr` (:523-525). So an insertion or deletion written in either place whose shift
/// carries it into the coding sequence has the shifted HGVSc and HGVS_OFFSET and no HGVSp, except
/// in the first on a transcript with a frameshift intron, whose exons the predicate reads
/// 12 bases wider (BaseTranscriptVariation.pm:861-866). The candidate writes the protein notation
/// of the shifted allele, as Ensembl VEP does for the allele written there. The shape: an empty
/// reference HGVSp against a candidate protein notation, on an indel's row annotated inside the
/// transcript with no CDS or protein position, whose HGVSc and non-zero HGVS_OFFSET both engines
/// write alike, the HGVSc a `c.` insertion or deletion notation of the row's transcript that places
/// the shifted allele in the coding sequence: a `del` whose every coordinate is a coding base; a
/// `dup` ending on a coding base, its first coordinate free to lie back in the intron or before the
/// coding start; an `ins` between two coding bases. The row is one of two kinds. With no cDNA
/// position the allele as written touches no exon and lies in an intron. With a cDNA position the
/// row is the class's where one of its terms is `5_prime_UTR_variant`: the allele is written on
/// exon bases of the 5' UTR, the class reads the row alone, and a `dup` there begins on a coding
/// base or before the coding start (`c.-2_12dup`), with or without an intron offset. An insertion
/// the shift leaves between the intron's last base and the exon's first
/// base is a change at the junction, written as a `dup` that ends on the intron's last base
/// or as an `ins` between that base and the coding base after it: the changed sequence
/// reads as a longer intron before the unchanged exon as well as a longer exon, the candidate
/// leaves the HGVSp empty there as Ensembl VEP does, and a protein notation on one side alone
/// is a difference, whatever bases are inserted and however long the shift. So is an
/// insertion before the stop codon's last base, whose intron base is printed `*-1`,
/// without that base's coordinate (TranscriptVariationAllele.pm:2854-2859). On a row with no cDNA
/// position the notation alone shows an intron inside the coding region where the first base of the
/// exon the allele is shifted into is coding base 2 or later. A notation with a coordinate in the
/// intron names that base (the coordinate counts back from it, or forward from the coding base
/// before the intron, which it follows); one of coding bases alone bounds it from below, since
/// the allele as written lies before it: the coding coordinate the shifted allele starts at,
/// less the bases of the shift, plus 1. A `dup` that ends on a coding base and begins before the
/// coding start (`c.-4+851_271dup`) is taken from the row alone too: its first coordinate is the
/// first base of the copy that stands before the shifted insertion point, and the allele as
/// written stands at or before that base where the shift is as long as the allele, or inside the
/// copy where it is shorter, so it is written in an intron of the 5' UTR or in one inside the
/// coding region that the copy spans. A notation of coding bases alone, or with an intron
/// coordinate counted from a coding base, whose bound is below 2 is read with the exon spans of
/// the candidate's JSON cache (`--json-cache`): the class takes the row where the allele as the
/// row's Location writes it lies wholly in the intron before the exon that holds the shifted
/// allele, whatever that exon's first coding base. Without the cache, or where the cache does not
/// place the row so, such a row stays a difference; the cache is not read for a row with a cDNA
/// position. A row with a CDS or protein position, a row with a cDNA position whose terms lack
/// `5_prime_UTR_variant`, a row whose SHIFT_LENGTH is positive (annotated at the shifted span), a
/// reference notation against an empty candidate and two differing notations stay differences. A
/// candidate notation at the start codon is taken as any other, the start codon lost or written as
/// unchanged (`p.Met1=`). A candidate notation at the stop codon is taken in two forms, each of
/// which keeps every residue before the stop: the stop written as unchanged (`p.Ter84=`) and an
/// extension named at the stop (`p.Ter84GlnextTer17`), which the candidate writes for a shifted
/// change that keeps the whole protein. Any other notation whose first changed residue is the
/// stop codon (`p.Ter84Gln`), one that names the stop codon as a flank of an insertion
/// (`_Ter85ins`) and one that ends a duplicated, deleted or replaced range at it stay
/// differences: bases inserted inside the stop codon that leave a stop codon after the new
/// residues are written like bases inserted before it. A row either engine flags
/// `cds_start_NF` (OutputFactory.pm:1487) is read as a row without a flag: on a coding
/// sequence that lacks its start and ends on its stop codon the candidate writes the notation
/// of the shifted change, its residues numbered on the codons read from the sequence's first
/// base. A row either engine flags `cds_end_NF` is a row of the class with any notation but
/// the two forms at the stop codon: on a coding sequence that does not end on its stop codon
/// the candidate writes the notation of a shifted change that the sequence's whole codons
/// give, a frameshift with an unknown count to the new stop among them, and Ensembl VEP reads
/// no insertion or deletion as overlapping the stop codon of a transcript so flagged
/// (VariationEffect.pm:1373, 1400), a reading the candidate follows. A coding sequence
/// without its closing stop codon shows in no field of a row that carries no flag.
pub const HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE: &str = "hgvsp_shifted_indel_written_span_gate";
/// `hgvsp_shifted_indel_written_span_start_test`: `hgvs_protein` shifts the allele 3' and reads
/// its translation coordinates and its peptides at the shift
/// (TranscriptVariationAllele.pm:1711, 1724-1733, 1774-1793), deletes the shift (:1796), reads
/// the frameshift from the coding coordinates the shift left (:2053) and then asks
/// `start_lost` (:2161), which answers for the allele as the record writes it: from the
/// verdict cached with the row's terms (Utils/VariationEffect.pm:855-862) or from the
/// unshifted cDNA span (:977-980), and with no test where that span has no start (:981). So an
/// insertion or deletion that passes the `coding` test as written
/// (TranscriptVariationAllele.pm:1731), has no start codon test there and reaches the first
/// codon at its 3' place alone is written as a frameshift at residue 1 with the unknown
/// count, or as a deletion of residues from the first. The
/// candidate makes the start codon test on the moved allele and writes the notation the
/// allele gets when its record writes it at the 3' place. The shape: on an insertion's or
/// deletion's row both engines give one non-zero HGVS_OFFSET and the same terms,
/// `5_prime_UTR_variant` among them and neither `start_lost` nor `start_retained_variant`,
/// two HGVSp values of one protein that are one of three pairs. On a row with no CDS position
/// (an insertion written at the last base before the coding start): the reference's a
/// frameshift at residue 1 with the unknown count (`p.Met1AsnfsTer?`) and the candidate's
/// that residue unchanged (`p.Met1=`, the `=` as the text formats escape it too), the coding
/// sequence standing whole after the inserted bases. On a row that holds
/// `coding_sequence_variant` and whose CDS position has an end and no start (a deletion that
/// begins before the first base of a transcript with no 5' UTR and takes its start codon):
/// the same frameshift against the lost start named at that residue, `p.Met1?` or
/// `p.Met1_?3`, its last residue that of the last coding base of the candidate's `c.`
/// deletion; or a deletion from residue 1 to residue n (`p.Met1_Pro3del`) against the lost
/// start over the same n residues, each named (`p.MetValPro1_?3`). A pair of prediction
/// parentheses both values carry is set aside. An empty value on either side, another
/// residue number, a count in the reference's frameshift, a candidate of another form, a row
/// with a term of a start codon test and the sides exchanged stay differences.
pub const HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST: &str =
    "hgvsp_shifted_indel_written_span_start_test";
/// `hgvsp_start_codon_test_past_the_shifted_allele`: under `--shift_3prime 1` the output step
/// has every allele serve its cDNA span with the shift (OutputFactory.pm:297-318,
/// TranscriptVariationAllele.pm:502-527), and `_inv_start_altered` adds the shift length to
/// that span again (Utils/VariationEffect.pm:920, 928-930), so the second start codon test of
/// `start_lost` (:871) edits the sequence the shift length past the allele. An insertion or
/// deletion that the shift moves onto the start codon and that leaves the whole coding sequence
/// standing (`_ins_del_start_altered`, :1071) then carries `start_retained_variant` and, where
/// the start codon still stands at the old coding start after that edit, no `start_lost`.
/// `_get_hgvs_peptides` writes the lost start only under `start_lost`
/// (TranscriptVariationAllele.pm:2161), so the notation is read from the coding sequence edited
/// without the 5' UTR (:2394-2439): a frameshift at residue 1 with the unknown count
/// (:2320-2365, 2003, 2483), or an insertion between a residue 0 and residue 1, both flanks
/// named as the last residue of the translation (:2231-2238, 2134-2144, 2383). The candidate
/// tests the start codon at the moved allele and writes the first residue unchanged, the
/// notation it gives every retained start whose coding sequence stands whole. The shape: on an
/// insertion's or deletion's row both engines give one positive SHIFT_LENGTH, the same terms,
/// `start_retained_variant` among them and no `start_lost`, a CDS position with a start and an
/// end, and two HGVSp values of one protein: the candidate's the first residue unchanged
/// (`p.Met1=`, the `=` as the text formats escape it too), the reference's a frameshift at
/// residue 1 with the unknown count and the same first residue (`p.Met1TrpfsTer?`), or an
/// insertion between residues 0 and 1 whose two flanks carry one name (`p.Ala0_Ala1insSerVal`).
/// A pair of prediction parentheses both values carry is set aside. An empty value on either
/// side, a count in the reference's frameshift, another residue number, unequal flank names, a
/// row that holds `start_lost` or lacks `start_retained_variant`, a missing, zero or unequal
/// SHIFT_LENGTH and the sides exchanged stay differences.
pub const HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE: &str =
    "hgvsp_start_codon_test_past_the_shifted_allele";
/// `hgvs_shifted_deletion_clamped_span`: `_var2transcript_slice_coords` clamps a written span
/// straddling the transcript's first base onto that base (TranscriptVariationAllele.pm:2732-2733)
/// before the offset moves it, so a deletion the shift carries inside is described from the base
/// after the shift where its shifted span starts earlier: two deleted bases print as one. The
/// shape: two `del` notations of one transcript and numbering ending at the same cDNA coordinate,
/// the reference's starting later than the candidate's, on a deletion's row (allele `-`) both
/// engines give a non-zero HGVS_OFFSET and, where the run prints one, a positive SHIFT_LENGTH; the
/// candidate's coordinates agree with the row's cDNA position, which under that run is the shifted
/// span's (an `n.` coordinate equal to the column's, a same-region `c.` span as long as it), so
/// a span off by one base stays a difference. The two notations start at different coordinates,
/// so the HGVSIntronOffset plugin's start offset differs with them on such a row: each engine's
/// is the intron offset of its own notation's start coordinate.
pub const HGVS_SHIFTED_DELETION_CLAMPED_SPAN: &str = "hgvs_shifted_deletion_clamped_span";
/// `hgvs_clipped_reference_against_whole_alternate`: `_var2transcript_slice_coords` clamps a
/// written span running over a transcript's edge onto the bases inside
/// (TranscriptVariationAllele.pm:2732-2733) and `hgvs_transcript` passes that span with the whole
/// alternate allele (:1474-1483), which `hgvs_variant_notation` reads against the reference bases
/// of the span alone (Utils/Sequence.pm:526), so a substitution of several bases running over an
/// edge, or a multi-allelic record's allele still carrying a shared outside base, is described as
/// a `delins` of every alternate base over the bases inside, as the substitution, `ins` or
/// `delins` the clip leaves of them (:2231-2242, :2271-2284), or, where every base inside is
/// unchanged, as an `ins` of the outside bases at the coordinate `_get_cDNA_position` gives
/// the base before the transcript's genomic start (:2800-2827, an intron offset off the first exon
/// in genomic order). The candidate cuts the alternate as the clamp cut the reference and
/// describes the change inside the transcript, none when there is none. The shape: a differing
/// HGVSc on a row of a sequence allele whose cDNA position has an open side (`?`), the
/// reference's notation naming the row's transcript, or empty where the whole alternate equals
/// the bases inside (Utils/Sequence.pm:529-531); the candidate's of the same transcript and
/// numbering, every coordinate within the row's span plus one of every coordinate of the
/// reference's (both engines clip against the bases inside, so both notations lie within the
/// allele's span over the transcript, give or take the base before exon 1 and the base after an
/// insertion site), or empty where the change lies outside; against the empty reference, a
/// `del`; against an empty candidate, an `ins` or a `dup`, what the outside bases read as
/// once the clip has stripped the unchanged bases inside. An allele shorter than its span that
/// deletes bases over the transcript's edge, the deletion written nowhere wholly inside the
/// transcript, has no description on it, a deletion extending beyond a transcript being the
/// genomic description's alone,
/// and the candidate writes none: against that empty candidate the reference's is a `del`, a
/// `delins` or a substitution of the bases inside, what the whole alternate reads as against
/// them (`deletion_over_the_edge`). The coordinate test reads one scale,
/// on which the stop codon or an intron inside a notation's span misstates a distance, and an
/// allele whose other end lies in an intron or past the other edge has no open side at all:
/// neither end maps to a cDNA base and the position is empty on both engines
/// (BaseTranscriptVariation.pm:145-148, Utils.pm:163-165), as it is on every intronic row. On
/// an edited model's row (`BAM_EDIT` `OK`) whose position is open at its end,
/// `hgvs_transcript` adds the model's offset to the clamped span's slice positions (:1517,
/// :1524) and `_get_cDNA_position` writes a position past the last exon as that exon's end,
/// numbered without the model's edits, with the distance as an intron offset (:2798-2829). A
/// reference coordinate whose intron offset is the row's REFSEQ_OFFSET, a positive whole
/// number, therefore stands in the coordinate test for every base from the exon base it counts
/// from to itself: the candidate numbers the model's last base with the bases the edits add
/// before it, between none and the offset, a count the row does not state. A distance that is
/// not the row's REFSEQ_OFFSET, a row without `BAM_EDIT` `OK` and a position open at its
/// start are held to the coordinate test as it stands. On such a model `hgvs_transcript` also
/// writes the row's whole reference allele over the notation's (:1510-1513), so a substitution
/// typed against the one base the clamp left inside names two or more reference bases at one
/// position (`c.*852AC>G`), a base outside the transcript among them; on a `BAM_EDIT` `OK` row
/// the class reads it as the substitution of one base at that position, against the
/// candidate's notation of the base inside. The row does not say which of the bases named
/// lies inside. On a
/// substitution's row the two notations show the clamp against the transcript's exon spans,
/// read from the candidate cache (`whole_alternate_over_fewer_bases`, `TranscriptEdges`):
/// each `delins` writes a run of the row's allele, and the offset between the two runs plus
/// the signed distance between the two start coordinates, and likewise at the end, is the count
/// of the allele's bases that lie outside the transcript at that end; at an end inside the
/// transcript the two notations name one coordinate, and each coordinate of either notation is
/// written as the spans place its base, in an exon or an intron. Such a row is the class whether
/// its position is open or empty; without the cache no row is.
/// A deletion's row (allele `-`) and a closed cDNA position stay differences.
/// The HGVSIntronOffset plugin's two fields are the intron parts of the notation's coordinates
/// (:1530-1533, :1543, :1545), so on a row whose HGVSc the class takes against a candidate that
/// prints none, the reference alone writes them and the candidate neither: they are the class
/// where they equal those intron parts (`clipped_notation_intron_offsets`). At a
/// reverse-strand transcript's first base, its genomic end, the clip leaves Ensembl VEP an
/// insertion at slice position 0, for which `_get_cDNA_position` returns no coordinate
/// (:2834): on such a transcript's row whose cDNA position is open at its start the reference
/// prints no notation where the candidate prints an `ins` of the bases a longer allele adds
/// inside the span, and the candidate alone writes the two fields.
/// The class takes the two fields of such a row, the intron parts of the candidate's notation,
/// and leaves its HGVSc a difference. Where both print an `ins` the class takes as a pair,
/// each side's fields are the intron parts of its own notation.
/// Before a forward-strand transcript's first base the clip leaves Ensembl VEP an insertion
/// that ends at slice position 0, the base before the first exon, which `_get_cDNA_position`
/// writes through its branch for a position between two exons as the first base's coordinate
/// with `-1` (:2800-2826): the end offset stored beside that `ins` is -1 where the
/// transcript has no intron. Where the candidate describes the bases such an allele adds as
/// the duplication the HGVS 3' rule gives, both coordinates exonic and the bases the rule
/// carried them in HGVS_OFFSET, its end offset is 0. The class takes the end offset of such a
/// row where the two rows show the clamp at the first base
/// (`insertion_before_the_first_base_against_duplication`); the row's HGVS_OFFSET is no
/// field of the class and its HGVSc is the class's only by the coordinate test above.
pub const HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE: &str =
    "hgvs_clipped_reference_against_whole_alternate";
/// `hgvs_added_bases_left_at_the_span_end`: `_return_3prime` shifts a record whose class is
/// insertion or deletion alone (TranscriptVariationAllele.pm:120-123), so the bases added by an
/// allele of any other record, one that keeps bases it shares with its reference, stay where
/// `_clip_alleles` leaves them, at the end of the record's span, typed an `ins` with no
/// duplication lookup (:2204-2284). Where that end is a reverse-strand transcript's last base,
/// `_get_cDNA_position` numbers the position past it off the last base with an intron offset
/// of 1 (:2802, :2820-2827): `c.*27+1_*27insAC`, beside which the HGVSIntronOffset plugin
/// writes 1 and 0. The candidate places such bases by the HGVS 3' rule where the run of points
/// at which they give the same sequence reaches the transcript's first or last base: a
/// duplication of the transcript's last bases where they repeat them, the notation at the most
/// 3' point of a run that starts at or before the first base, with the bases the rule moved
/// them in HGVS_OFFSET, and no notation where they have no placement inside the transcript.
/// The shapes are the predicate's (`added_bases_left_at_the_span_end`).
pub const HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END: &str = "hgvs_added_bases_left_at_the_span_end";
/// `hgvs_insertion_left_where_the_shift_stops`: `perform_shift` reads 1,000 bases of the
/// reference either side of an insertion (TranscriptVariationAllele.pm:427-437) and walks them
/// at most 1001 less the allele's length times, 1,000 times for a longer allele (:303-306);
/// going 5', on a reverse-strand transcript, the limit of an allele of 1,001 bases is 0, which
/// :304 leaves, and nothing is walked. An
/// insertion whose run of placements goes further is described where that walk stops, beside
/// an HGVS_OFFSET that is the bound itself, none where the bound is 0. The candidate walks to the
/// end of the run, the most
/// 3' placement, which the HGVS 3' rule assigns the change to: a `dup` or an `ins` further
/// 3' beside a larger HGVS_OFFSET, or no notation where that placement lies past the
/// transcript's last base. The HGVSIntronOffset plugin's two fields follow each notation.
/// A `dup` names a run of as many bases of the chromosome as the allele holds, each end
/// mapped to a coordinate alone (:1524-1525), so two exonic coordinates count those bases
/// less the bases of every intron between them. Given the candidate cache's exon spans, the
/// class holds the `dup` the reference writes for the allele that is not walked to such a
/// run beside the point, and takes a candidate `dup` of fewer bases than the allele holds
/// where it is one. A row whose transcript the cache does not hold on the row's strand
/// stays a difference there, as does a `dup` with one coordinate in the 3' UTR and the
/// other before it, whose distance the notation does not carry. Without the cache the class
/// reads the reference's `dup` by its count alone, from 2 bases up to the allele's own but
/// for two fewer, whatever intron would explain the count, and leaves that candidate `dup`
/// a difference. The HGVSp follows the placement as well: the reference reads it at the span
/// its bounded walk reaches (:1711) and writes none where that span has no translation start
/// or end (:1730-1738), and the candidate reads it where the HGVSc is described. Where the
/// last placement lies outside the coding sequence the reference's notation stands against
/// none, and where a record written outside the coding sequence has its last placement after
/// a coding base the candidate's frameshift, insertion, duplication or deletion-insertion
/// stands against none (`protein_notation_at_the_last_placement`). The shapes are the
/// predicate's (`insertion_left_where_the_shift_stops`).
pub const HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS: &str =
    "hgvs_insertion_left_where_the_shift_stops";
/// `hgvs_added_bases_written_outside_the_transcript`: OutputFactory.pm:1730 fills the HGVSc
/// column only on a row within the feature, and an insertion written directly beside a
/// transcript, or before its first base, is a row outside it
/// (BaseVariationFeatureOverlapAllele.pm:480), so the reference prints no HGVSc and no
/// HGVS_OFFSET there, while HGVSIntronOffset.pm:62 asks for the notation on every Transcript
/// row and writes the pair of a notation the row never prints, or none. The description of a
/// change does not turn on which of its placements the record writes: where the HGVS 3' rule
/// places the inserted bases on the transcript, as a duplication of its last bases or at a
/// point between two of its bases that a run from before its first base ends at, the
/// candidate writes that notation on the row, with the bases the rule moved it in HGVS_OFFSET
/// and the intron offsets of its own coordinates. The shape is the predicate's
/// (`added_bases_written_outside_the_transcript`).
pub const HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT: &str =
    "hgvs_added_bases_written_outside_the_transcript";
/// `hgvsp_stop_codon_protein_comparison`: `hgvs_protein` describes a change from the peptides
/// of the allele's codon window (`_get_hgvs_protein_type`, TranscriptVariationAllele.pm:2047,
/// and `_get_hgvs_protein_format`, :1898), where the HGVS nomenclature describes a protein
/// change by comparing the variant protein with the reference protein. At the stop codon the
/// window's description is not the nomenclature's in eleven shapes, all but the tenth on a row
/// whose two values name one protein. In the first five both name the stop at one residue number above
/// 1; the candidate writes the stop as unchanged where the protein does not change, and an
/// extension named at the stop (the first new residue, the count of residues added or `?`)
/// where the whole protein is kept and residues follow it. In the sixth the candidate writes a
/// frameshift, named at the first changed residue, where a change reaches from before the stop
/// codon into it. In the seventh both write an extension and the candidate counts its residues.
/// In the eighth the reference writes the stop as deleted before an extension's count, and the
/// candidate names the residue that stands in the stop's place. In the ninth the reference
/// writes residues and a stop as inserted before the first of them, the last residue or the
/// residues from an earlier one to the last, where an insertion brings a stop codon of its
/// own to the stop's residue and the protein does not change, and the candidate writes the
/// stop as unchanged. In the tenth the reference writes no value,
/// where an insertion inside the last residue's codon keeps that residue and brings a stop
/// codon to the stop's residue, and the candidate writes the stop as unchanged. In the
/// eleventh both write an extension at a stop that the transcript's peptide holds, or at the
/// stop codon that closes the coding sequence behind one, and the candidate counts its
/// residues from the lost stop; or both write one deletion of whole residues that takes such
/// a stop, with an extension's count after it, and the candidate's count is the place of the
/// first stop after the change less the lost stop's. The fifth to the eleventh shape read the
/// row's terms, which are the reference's; the ninth reads them for one inserted residue and
/// not for more.
///   The stop replaced by a stop (`Ter<N>delinsTer`) against the stop unchanged (`Ter<N>=`).
///   The stop alone replaced by the first residues that follow the protein, with or without a
///   closing `Ter` (`Ter<N>delins<residues>[Ter]`), against an extension with the first of
///   them and `?` or a count above the residues the reference names (their number where the
///   reference closes with no `Ter`).
///   The same where the reference names every added residue and then the stop the inserted
///   bases bring, against an extension whose count is the number of residues named.
///   Residues inserted between the last residue and the stop
///   (`<Aaa><N-1>_Ter<N>ins<residues>`) against an extension with the first of them and
///   their number.
///   The stop written as unchanged (`Ter<N>=`) on a row whose terms hold `stop_lost`
///   against an extension.
///   A range from a residue to the stop replaced by residues, none of them a stop
///   (`<Aaa><M>_Ter<N>delins<residues>`), on a row whose terms hold `stop_lost`, against a
///   frameshift named at that residue with the first of them (`<Aaa><M><Bbb>fsTer<count>`),
///   its count `?`, or above their number and other than `N - M + 1`, the count that leaves
///   the new stop at the stop's own residue.
///   An extension with an unknown count (`Ter<N><Xaa>extTer?`) on a row whose terms hold
///   `stop_lost` against the same extension with a count: `_stop_loss_extra_AA`
///   (TranscriptVariationAllele.pm:2477) ends at the first stop codon of the translated bases, so a
///   stop codon inside the coding sequence that the transcript's peptide goes on past,
///   reading it as a residue or holding it as a stop, leaves its count unknown where a stop
///   follows the lost one.
///   The stop written as deleted with an extension count and no residue in its place
///   (`Ter<N>delextTer<count>`) on a row whose terms hold `stop_lost`, against an
///   extension at the same residue with the same count and the first new residue.
///   One residue and a stop inserted between two neighbouring residues, the residue the right
///   flank's (`<Aaa><M>_<Bbb><M+1>ins<Bbb>Ter`), on a row whose terms hold
///   `stop_retained_variant`, against the stop unchanged at the residue after the right
///   flank (`Ter<M+2>=`): by the reference's value the alternate protein is the reference
///   protein up to the right flank's residue and then a stop, at the residue where the
///   candidate's value has the reference protein's stop.
///   Two or more residues and a stop inserted between two neighbouring residues
///   (`<Aaa><M>_<Bbb><M+1>ins<residues>Ter`), on a row of any terms, against the stop
///   unchanged (`Ter<N>=`), beside the candidate JSON cache (`DownstreamSequences`,
///   `--json-cache`): the coding sequence of the row's transcript, translated under its own
///   genetic code as the candidate reads it, spells the left flank at `M`, the right flank
///   and the inserted residues from `M+1` to the last, and a stop codon after them at `N`
///   (`stop_unchanged_where_inserted_residues_keep_the_protein`).
///   No value, on a row whose terms hold `stop_retained_variant`, whose protein position
///   names one codon `M`, whose coding positions lie inside that codon and whose amino acids
///   are the codon's residue against the same residue and then a stop (`V/V*...`), against
///   the stop unchanged at the residue after the codon (`Ter<M+1>=`), alone or inside
///   prediction parentheses: by the row the alternate protein is the reference protein up to
///   residue `M` and then a stop (`stop_unchanged_where_the_reference_writes_none`).
///   An extension with an unknown count or with a count (`Ter<N><Xaa>extTer<count>`), on a
///   row whose terms hold `stop_lost`, against the same extension with another count, beside
///   the candidate JSON cache (`DownstreamSequences`, `--json-cache`): the peptide of the
///   row's transcript holds a stop and is the translation of its coding sequence, residue `N`
///   is a stop it holds or the stop codon that closes the coding sequence, and the next stop
///   of the frame, sought through the 3' UTR and no farther, stands the candidate's count on.
///   The reference's count is held to the one `_stop_loss_extra_AA` gives on those bases: it
///   subtracts one and the peptide's length from the place of the first stop of the
///   translated bases (TranscriptVariationAllele.pm:2503 and 2513) and makes no count where
///   the figure is not above nought or `N` is 1 (2519 and 2483). The row's FLAGS are not read
///   and `N` may be 1 (`extension_counted_from_a_stop_the_peptide_holds`).
///   One deletion with an extension's count after it, of a range or of one residue
///   (`<Aaa><M>_<Bbb><N>delextTer<count>`, `<Aaa><M>delextTer<count>`), the count unknown or
///   a number, on a row whose terms hold `stop_lost`, against the same deletion with another
///   count, beside the candidate JSON cache: the peptide of the row's transcript holds a stop
///   and is the translation of its coding sequence; the coding bases the row names are three
///   for each residue the description names, and with them out of the coding
///   sequence the first residue that differs is `M`; the first stop of the reference at or
///   after `M` is a stop the peptide holds or the stop codon that closes the coding sequence;
///   and the first stop the alternate bases meet at or after `M`, sought through the 3' UTR
///   and no farther, stands the candidate's count of places after it. The reference's count
///   is held to the one `_stop_loss_extra_AA` gives on the alternate bases, by the same
///   lines, with no count where `M` is 1 (the same function).
/// A row whose FLAGS name `cds_end_NF`, where no stop codon closes the coding sequence, is a
/// row of the eleventh shape alone. The row alone does not say that the whole reference
/// protein is kept, that
/// a change runs to the stop codon's last base or shifts the reading frame, that the coding
/// sequence holds a stop codon its peptide reads as a residue, nor what a count is that the
/// reference does not spell out; the same notation with the sides exchanged, another
/// residue number or first residue, a range ending at the stop on a row whose terms hold no
/// `stop_lost`, and a reference that is empty off the tenth shape or writes the stop as deleted
/// with no extension count stay differences. For the ninth shape the reference's value names
/// the two flanks alone, so the stop's residue number is the candidate's, held to the one
/// after the right flank. More than one residue before the inserted stop is read beside the
/// cache alone, the row not naming the reference residues past the right flank: without the
/// cache, under a genetic code no table is held for, on a coding sequence no stop codon
/// closes and where a named residue is not the protein's at its number, the pair stays a
/// difference. So do one inserted residue that is not the right flank's, flanks that are not
/// two neighbouring residues, and one inserted residue on a row whose terms hold no
/// `stop_retained_variant`. For the tenth shape the row does not number the reference
/// protein's stop either, and the number is the candidate's, held to the one after the row's
/// codon; a row that names two codons, coding positions outside the codon, amino acids that
/// do not read the codon's residue and then a stop, and the same pair with the sides
/// exchanged stay differences.
pub const HGVSP_STOP_CODON_PROTEIN_COMPARISON: &str = "hgvsp_stop_codon_protein_comparison";
/// `hgvsp_stop_named_on_incomplete_coding_sequence`: `hgvs_protein` three-letters the `X` of
/// a partial codon and of an unknown residue as `Xaa` and rewrites it `Ter`
/// (`s/Xaa/Ter/g`), so on a coding sequence that no stop codon closes, and at residue 1 of
/// one that opens inside its first codon, it names a stop the reference protein does not hold.
/// The HGVS nomenclature takes as a reference only a transcript with a complete coding
/// sequence, and the candidate writes no HGVSp there. A row of the class: the reference's
/// description names the stop alone, as unchanged (`Ter<N>=`), deleted (`Ter<N>del`),
/// replaced by one residue (`Ter<N><Xaa>`) or extended (`Ter<N><Xaa>extTer<count>`); the
/// candidate's HGVSp is empty; and the row shows the incomplete coding sequence: its FLAGS name
/// `cds_end_NF`, its terms hold `incomplete_terminal_codon_variant`, or the description is
/// `Ter1=` and its FLAGS name `cds_start_NF`. The row alone does not show a coding sequence
/// that no flag marks incomplete: the candidate's JSON cache does (`PartialCodonModels`,
/// `--json-cache`), where it holds the row's transcript with a coding sequence that is no
/// whole number of codons. A frameshift that the 3' shift carries past the last whole codon
/// is such a row: `_get_fs_peptides` closes the reference translation with a stop whatever
/// the coding sequence's end (`TranscriptVariationAllele.pm:2338`). Without the cache such a
/// row stays a difference, and so do any other reference description, a candidate that
/// writes a notation and the sides exchanged.
pub const HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE: &str =
    "hgvsp_stop_named_on_incomplete_coding_sequence";
/// `hgvsp_frameshift_count_past_a_stop_read_as_a_residue`: `_stop_loss_extra_AA` counts a
/// frameshift's residues to the first stop of the alternate sequence's translation
/// (TranscriptVariationAllele.pm:2492, :2503 and :2508), in which a stop codon that the
/// transcript's peptide goes on past is a stop: one the peptide reads as a residue, a
/// selenocysteine codon among them, or one it holds as a stop before the residues it names
/// after it, as a polymorphic pseudogene's does. After such a codon the count is not positive
/// and Ensembl VEP writes the unknown count (:2519 and :2006), also where the new reading frame
/// meets a stop codon. The candidate counts from the first changed residue to that stop where
/// it is another codon than every stop codon the peptide reads as a residue. A row of the
/// class: the reference's description is a frameshift with the unknown count; the candidate's
/// is the same frameshift, at the same residue and with the same new residue, with a count; and
/// the coding sequence and 3' UTR of the row's transcript in the candidate's cache, with the
/// candidate row's allele at its coding position, give that count to such a stop past a stop
/// codon inside the peptide's residues. A frameshift whose first changed residue is a stop
/// is written by both as an extension named at that stop (:2007 to :2010), with the same
/// count less the stop itself, and is a row of the class too: on a row whose terms hold
/// `frameshift_variant`, the reference's description is that extension with the unknown
/// count and the candidate's the same extension, at the same stop and with the same new
/// residue, with a count; the peptide holds a stop and is the translation of its coding
/// sequence at every codon, its first residue a methionine whatever the first codon; a stop
/// stands before the stop named, which is a stop the peptide holds or the stop codon that
/// closes the coding sequence; and the sequences give that count from the lost stop to the
/// new frame's first stop. The class reads the cache and accepts nothing without one. A
/// count the sequence does not give, a new stop that is the same codon as one the peptide
/// reads as a residue and the sides exchanged stay differences.
pub const HGVSP_FRAMESHIFT_COUNT_PAST_A_STOP_READ_AS_A_RESIDUE: &str =
    "hgvsp_frameshift_count_past_a_stop_read_as_a_residue";
/// `hgvsp_frameshift_with_no_changed_residue`: `_get_fs_peptides` reads the reference and
/// the alternate translation residue by residue from the change on and leaves the walk by
/// its bound where the alternate translation ends before a residue differs
/// (TranscriptVariationAllele.pm:2348). It returns the pair of residues read last, two equal
/// ones, at the position after them (:2351-2352, :2361 and :2364), and Ensembl VEP writes a
/// frameshift that names one residue on both sides of its number, with the unknown count.
/// The walk's own exit is on two different residues (:2360), so no other path writes a
/// frameshift that names one residue twice. The HGVS nomenclature names a frameshift by its
/// first changed residue, which such a translation does not hold, and the candidate writes
/// no HGVSp there. A row of the class: the reference's description is a frameshift with the
/// unknown count whose two residue names are one name and no stop's (`<Xaa><N><Xaa>fsTer?`);
/// the candidate's HGVSp is empty; and the row's terms hold `frameshift_variant` on both
/// sides. The class reads the row alone. The sides exchanged, a candidate that writes a
/// notation, two different residue names, a count that is a number and any other
/// description stay differences.
pub const HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE: &str =
    "hgvsp_frameshift_with_no_changed_residue";
/// `hgvsp_frameshift_count_at_the_first_residue`: under --hgvs, `_stop_loss_extra_AA` makes
/// no count for the position before the first residue (`return undef unless $ref_var_pos`,
/// TranscriptVariationAllele.pm:2483, handed 0 at :2003), so Ensembl VEP writes the unknown
/// count for every frameshift whose first changed residue is residue 1 (:2006), also where
/// the new reading frame meets a stop codon. On a transcript at whose coding start no start
/// codon stands, one that carries `cds_start_NF` or whose record shows none, the candidate
/// counts the residues from residue 1 to the first stop codon of the new frame, sought
/// through the coding sequence and the 3' UTR and no farther; where a start codon stands it
/// keeps the unknown count. A row of the class: neither row's terms hold `start_lost`; both
/// rows name one protein, the reference's description is a frameshift named at residue 1
/// with the unknown count and the candidate's the same frameshift with a number; and the
/// record of the row's transcript in the candidate's cache (`PartialCodonModels`,
/// `--json-cache`), which shows no start codon, gives that residue, that new residue and
/// that number for the row's allele over its `CDS_position`, a span of as many positions as
/// the row's record has reference bases, which tells an insertion between two positions
/// from a replacement of both: the class is one of `RECORD_REFERENCE_CLASSES`. A deletion
/// whose `CDS_position` is open at its start is read on the span its HGVS_OFFSET moves it
/// to, as long as those reference bases. The class reads the cache and accepts nothing
/// without one, and takes the notation whatever other terms the two rows hold. A transcript
/// on which a start codon stands, a lost start on either row, another count, residue or
/// number, a count on the reference's side, a first codon the start exon's phase pads, a new
/// frame that meets no stop codon, a codon before that stop codon that holds a letter that
/// is no base, an edit that leaves fewer than three coding bases, a row without its input
/// record and the sides exchanged stay differences.
pub const HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE: &str =
    "hgvsp_frameshift_count_at_the_first_residue";
/// `hgvsp_intron_end_insertion`: under --hgvs, an insertion the input writes between an intron's
/// last base and the first base of the next exon in the transcript's direction, that base inside
/// the coding region, which the 3' rule leaves in that gap. Ensembl VEP reads the inserted bases
/// as the exon's first bases and writes the HGVSp of that reading. The candidate writes none,
/// because the HGVS nomenclature gives a change at an intron's end no predicted protein
/// consequence. Both engines write the same HGVSc, which names the form, `N` the exon's first
/// coding base: a copy of the intron's last base (`c.N-1dup`); a duplication of two or more
/// intron positions ending on that base (`c.N-K_N-1dup`, or `c.M+J_N-1dup` with `M` the
/// coding base before `N` where the copy starts in the intron's first half); an insertion
/// after that base (`c.N-1_Nins`) of one base, or of bases that do not end in one of the two
/// pairs that end an intron; or an insertion there of two or more bases that end in such a pair.
/// Every `c.N-1_Nins` insertion is therefore of a form. The candidate keeps the notation at the
/// end of an intron of 13 bases or fewer and on a row the 3' rule moved, which carries an
/// HGVS_OFFSET, so an HGVSp that differs there stays a difference. The shape: an HGVSp the
/// reference fills and the candidate leaves empty, on a row whose HGVSc is one notation of the
/// row's transcript on both sides, of one of the four forms, beside no HGVS_OFFSET, the record
/// written in such a gap of that transcript. A record that carries another allele beside the
/// insertion holds it with the reference base it shares, the allele that base and the inserted
/// bases: both engines write the `c.N-1_Nins` notation of the allele clipped of that base, and
/// the class takes the row by the gap the clipped insertion lies in (`held_insertion_gaps`).
/// A text row gives such a record's Location, not the gap, and spells its allele on the record's
/// strand, so the class places it by the Location's span: the gaps its allele's clippings leave
/// on either strand, or, where the allele is the insertion alone (a minimised allele), every gap
/// of the span (`lone_insertion_gaps`).
/// A span of two adjacent positions with a minimised allele is an insertion written alone or a
/// reference allele of two bases, so it names its three gaps.
/// A VCF or JSON entry spells its allele in more ways than a text row does, so the class places it
/// from its input line, a VCF line or an Ensembl-format one, by the notation alone: an ALT that
/// clips against the REF to a pure insertion names a gap, which the entry claims where the
/// inserted bases are the notation's or as many as the positions a `dup` notation duplicates
/// (`notation_insertion_gaps`).
/// Neither output carries the intron or the coding bounds, so the class reads them from the
/// candidate's JSON cache (`IntronEnds`, `--json-cache`) and accepts nothing without it.
pub const HGVSP_INTRON_END_INSERTION: &str = "hgvsp_intron_end_insertion";
/// `hgvsp_intron_start_insertion`: under --hgvs, an insertion the 3' rule leaves between a
/// coding exon's last base and the first base of the intron after it, written in that gap or
/// inside the exon and carried to it, whose bases are no copy of the bases before the gap.
/// Ensembl VEP reads the inserted bases as the exon's last bases and writes the HGVSp of that
/// reading. The candidate writes none, because the HGVS nomenclature gives a change at an
/// intron's start no predicted protein consequence. A copy of the bases before the gap is a
/// duplication of coding positions and keeps its notation on both engines. Both engines write
/// the same HGVSc, `c.N_N+1ins` with `N` the exon's last coding base in both coordinates,
/// and, on a record the 3' rule moved, the same HGVS_OFFSET: the bases moved, positive on a
/// forward-strand transcript and negative on a reverse-strand one, so the gap the insertion is
/// left in is the gap the record is written in moved by the offset as signed. A run that prints
/// the moved location (`--shift_genomic 1`) names the gap the insertion is left in, so the
/// class accepts a Location's gap as named or moved by the offset. The intron's length plays no
/// part. The shape: an HGVSp the reference fills and the candidate leaves empty, on a row whose
/// HGVSc is one notation of the row's transcript on both sides, of that form, beside no
/// HGVS_OFFSET or one non-zero integer on both sides, the gap the insertion is left in one in
/// which an exon of that transcript ends inside the coding region before an intron. A
/// notation on both sides, an HGVSc of another form and offsets that differ stay differences.
/// A record that holds the insertion with a reference base it shares carries no offset and is
/// taken by the gap its clipped allele lies in (`held_insertion_gaps`).
/// A text row of such a record is placed by its Location's span, as in the intron-end class; a
/// minimised allele the 3' rule moved is the notation's bases turned back by the offset
/// (`lone_insertion_gaps`).
/// A run under `--shift_3prime 1` prints a moved record's allele as the rule moved it, the
/// notation's bases as they stand, so on a row that carries the offset the class reads the bases
/// either way, a text row's minimised allele against the notation's with its span's gaps as named
/// or moved and, for a VCF or JSON entry, the bases an ALT of its input line inserts against the
/// notation's turned back by the offset (`notation_insertion_gaps`).
/// The class reads the exons and the coding bounds from the candidate's JSON cache
/// (`IntronEnds`, `--json-cache`) and accepts nothing without it.
pub const HGVSP_INTRON_START_INSERTION: &str = "hgvsp_intron_start_insertion";
/// `shift_intron_term_unshifted_span`: under --shift_3prime 1, `_intron_effects`
/// (BaseTranscriptVariationAllele.pm:99-150) sets `intronic` when the shifted span or the
/// unshifted span (:123-124, tested at :146-147) lies in an intron's interior, while every other
/// key of the method reads the shifted span alone; an insertion or deletion in an intron that
/// the shift carries onto the two bases of the acceptor site, or past them onto the exon's first
/// bases, keeps `intron_variant` beside the shifted span's terms (an intron is tested for a
/// span beginning within three bases past its end, the margin of the intron interval tree,
/// BaseTranscriptVariation.pm:1027). The candidate reads the shifted span alone, as Ensembl VEP
/// does for the allele written at the shifted position. The shape: a Consequence set the
/// reference writes as the candidate's plus `intron_variant`, on a row whose SHIFT_LENGTH is
/// positive and whose HGVSc places the shifted span at the acceptor's two bases or in the exon
/// (every coordinate exonic or at the intron offset -1 or -2, a duplication's first excepted,
/// since its bases run back into the intron; a notation without an offset on a row without a
/// splice-site term, since a deletion over a whole intron is numbered over the exons alone and
/// covers both sites), so that no intron base away from the splice site is touched; a deletion
/// reaching such a base carries the term from its shifted span on both engines. A row without
/// an HGVSc is read by its allele and terms: an insertion's shifted span lies between the two
/// acceptor bases whenever the row carries `splice_acceptor_variant`.
pub const SHIFT_INTRON_TERM_UNSHIFTED_SPAN: &str = "shift_intron_term_unshifted_span";
/// `mature_mirna_unshifted_span`: under --shift_3prime 1, the predicates read the shifted span
/// from `$bvf->{start}` and `{end}`, which BaseVariationFeatureOverlapAllele.pm:269-282 offset
/// for each predicate call; `within_mature_miRNA` (Utils/VariationEffect.pm:532-558) overlaps
/// the mature miRNA's genomic span with `$bvf->seq_region_start` and `seq_region_end`
/// (:545-549) instead, which BaseVariationFeature.pm:279-282 and :304-307 cache at their first
/// read, before any shift, so the test reads the unshifted span. An insertion or deletion the
/// shift carries off the mature sequence keeps `mature_miRNA_variant` where the shifted span
/// is `non_coding_transcript_exon_variant` (`non_coding_exon_variant`, :502-519, reads the
/// shifted span and yields to the mature test), and one carried onto it keeps the exon term. The
/// candidate reads the shifted span alone, as Ensembl VEP does for the allele written at the
/// shifted position. The shape: two Consequence sets equal but for one engine's
/// `mature_miRNA_variant` against the other's `non_coding_transcript_exon_variant`, on a row
/// whose SHIFT_LENGTH is positive, which only an insertion or deletion is given.
pub const MATURE_MIRNA_UNSHIFTED_SPAN: &str = "mature_mirna_unshifted_span";
/// `stop_codon_unshifted_span`: `_overlaps_stop_codon_cil` (Utils/VariationEffect.pm
/// :1389-1430) overlaps the stop codon's genomic span with `$bvf->seq_region_start` and
/// `seq_region_end` (:1403, :1424-1427), the cached unshifted coordinates, where every other
/// test of `stop_retained` (:1290-1323) reads the shifted span. A deletion the shift carries
/// onto the stop codon from the bases before it therefore fails the overlap: `stop_retained` is
/// false, and `coding_unknown` (:1616-1649) writes `coding_sequence_variant` for a span whose
/// shifted cDNA start is undefined. The candidate reads the shifted span, as Ensembl VEP does for
/// the allele written at the shifted position. The shape: two Consequence sets equal but for the
/// reference's `coding_sequence_variant` against the candidate's `stop_retained_variant`, on
/// a row whose SHIFT_LENGTH is positive. The other direction cannot arise: a span the 3' shift
/// carries off the stop codon lies in the 3' UTR, where the term is not a candidate.
pub const STOP_CODON_UNSHIFTED_SPAN: &str = "stop_codon_unshifted_span";
/// The five classes of pair at a coding end where vep-rs writes a stop codon's or a start
/// codon's terms only where that codon stands. Each reads the coding ends of the row's
/// transcript from the candidate's JSON cache (`PartialCodonModels::coding_ends`,
/// `--json-cache`), by the record of the row's Feature whose span holds the row's Location,
/// and takes a pair only where vep-rs's terms of the coding family are the ones derived for
/// the row (`coding_ends::terms::expected_terms`) and every other term is equal on the two
/// rows; `coding_ends::terms::CodingEndClass::takes` states each class's premise and its
/// reading of Ensembl VEP's row. A sequence allele's row is derived from the five cells
/// both rows print (`cDNA_position`, `CDS_position`, `Protein_position`, `Amino_acids`,
/// `Codons`) and the length of the record's reference bases, which a row under an
/// identifier takes from the records of `--input`. Without the cache, without a record of
/// the transcript at the row's Location, on a sequence allele's row of an output without
/// the five cells and on a row that is not derived, a class takes nothing. A row of two
/// classes is the first's in this
/// order. On a row a class took, the IMPACT is each side's own terms', a Downstream, NMD or
/// LOFTEE field stands on the side whose terms let the plugin run, and the JSON record's
/// `most_severe_consequence` is each side's most severe term over its own entries of the
/// record; every other cell is compared.
///
/// `structural_deletion_stop_lost_without_stop_codon`: Ensembl VEP writes `stop_lost` for a
/// structural deletion over the last three genomic positions of the coding region, with no
/// test of what stands there (Utils/VariationEffect.pm:1268-1280). The class takes a
/// structural allele's row on a transcript at whose coding end no stop codon stands, with
/// `cds_end_NF` or without, or whose stop codon polyadenylation completes where the span
/// holds none of the one or two bases that begin it, where vep-rs writes the same terms
/// without `stop_lost`.
pub const STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON: &str =
    "structural_deletion_stop_lost_without_stop_codon";
/// `stop_codon_completed_by_polyadenylation`: on a transcript whose coding sequence ends one
/// or two bases into a stop codon that polyadenylation completes, Ensembl VEP reads the
/// last codon as incomplete (Utils/VariationEffect.pm:1588-1614) and vep-rs reads the row
/// on the coding sequence followed by adenines. The class takes a sequence allele's row
/// whose reference terms hold `incomplete_terminal_codon_variant`, or whose reference
/// `Amino_acids` ends on the unknown residue at the protein position after the whole
/// codons. The candidate's `Amino_acids`, `Codons` and protein notation are the ones
/// derived on the completed coding sequence. A deletion from the stop codon past the
/// transcript's last base is derived with the reference bases of `--reference-fasta`, and
/// at the record's own place without them.
pub const STOP_CODON_COMPLETED_BY_POLYADENYLATION: &str = "stop_codon_completed_by_polyadenylation";
/// `stop_codon_overlap_without_stop_codon`: Ensembl VEP reads the last three coding
/// positions as the stop codon on every transcript without `cds_end_NF`
/// (Utils/VariationEffect.pm:1373-1381 with 1466-1480; 1400, 1408-1427 and 1508 read at
/// 1318). The class takes a sequence allele's row on a transcript without the flag at whose
/// coding end no stop codon stands, or whose stop codon polyadenylation completes, where the
/// reference terms hold `stop_lost` or `stop_retained_variant` by that overlap and the
/// derived terms hold neither. Where the transcript's last base cuts a frameshift's first
/// changed codon short, the reference's protein notation is the deleted residues with an
/// extension of unknown length and the candidate's the first of those residues named as a
/// frameshift.
pub const STOP_CODON_OVERLAP_WITHOUT_STOP_CODON: &str = "stop_codon_overlap_without_stop_codon";
/// `stop_codon_behind_last_residue_without_stop_codon`: Ensembl VEP reads a stop codon that
/// an allele brings behind the last residue as the transcript's own stop codon kept
/// (Utils/VariationEffect.pm:1344 and the first clause of 1357). The class takes a sequence
/// allele's row on a transcript at whose coding end no stop codon stands, with `cds_end_NF`
/// or without, where the reference terms hold `stop_retained_variant` by the peptide
/// comparison and the derived terms hold none.
pub const STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON: &str =
    "stop_codon_behind_last_residue_without_stop_codon";
/// `start_codon_overlap_without_start_codon`: Ensembl VEP reads the first three coding
/// positions as the start codon on every transcript without `cds_start_NF`
/// (Utils/VariationEffect.pm:975-986). The class takes a row on a transcript without the
/// flag that shows no start codon, where the reference terms hold `start_lost` or
/// `start_retained_variant` and the derived terms hold neither. On a substitution of one
/// residue at protein position 1 the reference's protein notation is the lost start and
/// the candidate's names the two residues. On a frameshift under the standard genetic code
/// the reference's is the lost start at the residue the candidate's frameshift notation
/// opens with, written on a row without HGVS_OFFSET as a range to the end of the
/// translation span where that is another residue.
pub const START_CODON_OVERLAP_WITHOUT_START_CODON: &str = "start_codon_overlap_without_start_codon";
/// The classes of the stop codon tests behind `stop_retained` (Utils/VariationEffect.pm
/// :1290-1324); the comment above each one's predicate states its shape.
pub const STOP_RETAINED_WITHOUT_STOP_CODON_TEST: &str = "stop_retained_without_stop_codon_test";
pub const STOP_CODON_OVERLAP_GENOMIC_SPAN: &str = "stop_codon_overlap_genomic_span";
pub const STOP_RETAINED_THREE_BASE_STOP_CODON_TEST: &str =
    "stop_retained_three_base_stop_codon_test";
pub const STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST: &str =
    "stop_codon_removed_without_stop_codon_test";
pub const STOP_CODON_OVERLAP_INSERTION_RANGE: &str = "stop_codon_overlap_insertion_range";
/// `cnv_tr_symbolic_loss_named_a_gain`: `copy_number_gain` (Utils/VariationEffect.pm:246-256)
/// holds for every allele of the class `tandem_repeat`, whatever its counts, and
/// BaseVariationFeatureOverlapAllele.pm:429-436 give such an allele the keys `insertion` and
/// `increase_length`. `feature_elongation` (Utils/VariationEffect.pm:337-348) is therefore
/// the term of every `<CNV:TR>` allele Ensembl VEP reads symbolically, and
/// `feature_truncation` (:350-373), which asks for a loss or a deletion, is never its term.
/// The candidate writes `feature_truncation` for such an allele whose record states fewer
/// bases than the run it spans, and neither term for one whose record states exactly as
/// many. The shape: a Transcript row whose Allele is `tandem_repeat` on both sides and whose
/// two lists of terms differ in the reference's one `feature_elongation` alone, where the
/// candidate writes `feature_truncation` in its place or no term. The class takes the
/// Consequence of such a row, and the JSON record's `most_severe_consequence` is left
/// uncompared on its record, as on the record of every row class. The two terms carry the
/// impact HIGH (Utils/Constants.pm:803-834), so no IMPACT follows a row that names the
/// truncation. A row that names neither carries the IMPACT of the first term left by rank
/// (OutputFactory.pm:2043-2050): the class takes it where each side's IMPACT is the one its
/// terms give, and takes no row whose candidate IMPACT is another. A row of another Allele,
/// the terms the other way round, a second place that differs and every other cell stay
/// differences.
pub const CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN: &str = "cnv_tr_symbolic_loss_named_a_gain";
/// `splice_term_past_the_record_span`, the class of the six terms an intron gives a row: the
/// five splice terms whose windows lie inside it and the intron term, whose window it is
/// (VariationFeatureOverlapAllele.pm:405-460, BaseTranscriptVariationAllele.pm:117-118 and
/// :143-206). On a record whose alternate allele is longer than its reference allele, Ensembl
/// VEP writes a splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine
/// tract or intron term from positions past the end of the record's span; where no base of the
/// span lies in a window of that term, vep-rs does not write it, unless it is a term of an
/// insertion written after the span. A record whose alternate allele opens with its whole
/// reference allele changes no base of its span. vep-rs reads it as the
/// insertion it writes, between the span's last base and the base after it, and writes that
/// insertion's terms. The comment above the predicate
/// (`shared::splice_term_past_the_record_span_pair`) states the shape with the record's
/// alleles and without them.
pub const SPLICE_TERM_PAST_THE_RECORD_SPAN: &str = "splice_term_past_the_record_span";
/// `stop_codon_settled_before_unknown_bases`: `peptide` returns nothing for an allele that
/// holds a letter outside A, C, G and T (TranscriptVariationAllele.pm:701), wherever the letter
/// stands, so Ensembl VEP writes no residues for it, no term that reads the peptide
/// (`stop_gained`, Utils/VariationEffect.pm:1211-1232; `inframe_insertion`, :1100-1159), the
/// term of an unknown change where the length is whole codons (`coding_unknown`, :1616-1649)
/// and no protein notation (TranscriptVariationAllele.pm:1683-1684). Where the whole codons
/// of the allele's codon string hold a stop codon before the first codon that holds such a
/// letter, the protein ends at that stop whatever the unknown bases are, and vep-rs writes
/// the row every completion of them shares: the residues, `stop_gained` beside the term the
/// length gives, and under `--hgvs` the protein notation with its HGVS_OFFSET. Where that
/// stop codon stands at the transcript's own behind the reference's residues, the inserted
/// bases keep the stop codon, and the row every completion shares is `inframe_insertion`,
/// which cuts the alternate peptide after its first stop (Utils/VariationEffect.pm:1126), beside
/// `stop_retained_variant`, with the residues and the notation of a stop codon left
/// unchanged. Ensembl VEP writes `stop_retained_variant` alone on such a row from the second
/// arm of `stop_retained` (:1318), which makes no edit test for a letter that names no base
/// (:1498); where the insertion's moved start fails that arm's overlap (:1408-1415), it
/// writes `frameshift_variant` by the length (:1564) or the term of an unknown change. The
/// shape (`rules::settled_stop::stop_codon_settled_before_unknown_bases_pair`): a Transcript
/// row whose `Codons` cell is the same on both sides, its alternate half holding such a
/// letter after a stop codon the standard genetic code reads among the whole codons before
/// it; no residues in the reference and, in the candidate, the translation of that half; and
/// terms that differ by the gained stop alone or by the kept stop alone. A gained stop is
/// `frameshift_variant` on both sides with `stop_gained` added, or
/// `coding_sequence_variant` in the reference against `stop_gained` and at most one of
/// `inframe_insertion`, `inframe_deletion` and `protein_altering_variant` in the candidate.
/// Where the reference writes `stop_retained_variant` on an insertion's or a deletion's row
/// for a stop codon test it does not make (Utils/VariationEffect.pm:1318, :1498), both rows
/// showing so, the candidate writes `stop_gained` in that term's place, with
/// `frameshift_variant`, one of those three terms or no other term.
/// A kept stop is the candidate's `inframe_insertion` added to terms that hold
/// `stop_retained_variant` on both sides, or the candidate's `inframe_insertion` and
/// `stop_retained_variant` against the reference's `frameshift_variant` or
/// `coding_sequence_variant` on an insertion between two neighbouring positions
/// (`shared::inserted_between_neighbours`) with both cDNA positions and a CDS start on both
/// sides. The class takes the Consequence of such a row, the residues, the
/// cells that follow the two sets of terms as a class of `TERM_PAIR_FOLLOWING_CLASSES`, a
/// protein notation that names a stop where the reference writes none, and the HGVS_OFFSET
/// written beside that notation. Where such an allele is an insertion written in an intron
/// that the HGVS 3' shift carries into the coding sequence, the row holds no `Codons` cell:
/// Ensembl VEP writes no notation and no HGVS_OFFSET (TranscriptVariationAllele.pm:1397,
/// :1684; OutputFactory.pm:1764, :1776), and vep-rs the protein notation read where the
/// shift ends with the offset beside it. Beside the candidate's JSON cache the class takes
/// those two cells of such a row (`rules::settled_stop::carried_insertion_cell`), whose two
/// sets of terms are equal and which holds no position, where the cache's transcript places
/// the insertion, moved by an offset smaller than the count of known bases the allele opens
/// with, between two coding bases of one exon, and the allele read there in the
/// transcript's frame settles a stop codon before its first codon with an unknown base, a
/// frameshift's residue number and count held to that reading; without the cache such a
/// row stays a difference. The rule of `stop_codon_overlap_insertion_range` stands
/// before this class's and takes the terms of a kept stop row whose reference writes
/// `frameshift_variant`; beside that class this one takes the residues, the notation and the
/// offset of the row. A row whose `Codons` cells differ, an allele whose known
/// bases settle no stop codon before its first unknown base, a transcript whose own genetic
/// code reads a codon of that half otherwise than the standard code, residues other than
/// the codons' translation, any other difference of terms and a coding notation one side
/// alone writes stay differences.
pub const STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES: &str = "stop_codon_settled_before_unknown_bases";
/// `frameshift_inside_a_held_stop_codon`: an insertion or a deletion whose length change is no
/// multiple of three and whose first base lies inside the codon of a stop that the
/// transcript's peptide holds as a stop and reads on past (a polymorphic pseudogene's).
/// Ensembl VEP writes no `frameshift_variant` on that row: its `frameshift` returns 0
/// wherever the reference peptide of the allele's codons begins with a stop
/// (Utils/VariationEffect.pm:1554), which holds for the stop codon no residue follows. Its
/// row reads `coding_sequence_variant`, `inframe_deletion`, or `inframe_insertion` with
/// `coding_sequence_variant`, and its HGVSp names the stop unchanged or a deletion-insertion.
/// vep-rs writes `frameshift_variant` and names the first residue that changes (an extension
/// where that residue is the stop, a frameshift where the new frame keeps the stop). The
/// shape (`rules::frameshift_at_a_held_stop::frameshift_inside_a_held_stop_codon_pair`): the
/// candidate's terms hold `frameshift_variant` and the reference's do not; the candidate's
/// other terms are the reference's without one of those three groups, which the reference
/// holds whole, every other term standing on both sides; the reference holds neither
/// `stop_retained_variant` nor `stop_lost`; the `Amino_acids` cell is the same on both
/// sides, opens with `*` and holds `X` after its `/`; and the `Codons` cell is the same on
/// both sides, its two halves differing in length by no multiple of three. The class reads
/// the two rows alone, no cache and no FASTA. It takes the Consequence of such a row and, as
/// a class of `TERM_PAIR_FOLLOWING_CLASSES`, the IMPACT each side's terms give, an HGVSp of
/// any value on either side, and a field of Downstream, NMD or LOFTEE that the candidate
/// alone writes where the terms let the plugin run on its row alone, an empty
/// `DownstreamProtein` included. A row whose reference holds `stop_retained_variant`,
/// `stop_lost` or `frameshift_variant`, a term one side alone holds beside the pair's own,
/// residues or codons that differ between the sides, a length change of whole codons and
/// every other cell stay differences.
pub const FRAMESHIFT_INSIDE_A_HELD_STOP_CODON: &str = "frameshift_inside_a_held_stop_codon";
/// `start_codon_removed_without_five_prime_utr`: for a deletion that removes whole codons
/// from the first coding base of a transcript whose coding sequence opens on the
/// transcript's first base, no line of `start_lost` (Utils/VariationEffect.pm:851-904)
/// finds the start codon kept and none names it lost. The length test, which reads the
/// start as altered (:1069), is set aside for an in-frame deletion (:870);
/// `_inv_start_altered` (:906-945) returns before its edit where the transcript has no
/// 5' UTR (:926); and the peptide test returns on the empty alternate peptide (:875).
/// Ensembl VEP writes `inframe_deletion` without `start_lost`. vep-rs reads the edited
/// coding sequence and writes `start_lost` beside `inframe_deletion` where the first codon
/// is a start codon of the transcript's genetic code that codes methionine and the edited
/// coding sequence opens on no triplet that codes it: the row Ensembl VEP writes for the
/// same removal behind a 5' UTR. The shape
/// (`rules::removed_start::start_codon_removed_without_five_prime_utr_pair`): a Transcript
/// row whose candidate terms are the reference's and `start_lost`, both sets holding
/// `inframe_deletion` and neither `start_retained_variant`; a cDNA position and a CDS
/// position that begin at 1 and end at a number, the same on both sides; and an
/// `Amino_acids` cell, the same on both sides, whose alternate half is `-`. The class
/// takes the Consequence of such a row and, as a class of `TERM_PAIR_FOLLOWING_CLASSES`,
/// the IMPACT each side's terms give and the protein notation of a run under `--hgvs`: a
/// deletion from residue 1 in the reference (TranscriptVariationAllele.pm:2018-2027)
/// against the lost start of the same residues in the candidate (:2161-2165, :2030-2037),
/// both naming one protein. A row whose coding sequence opens past the transcript's first
/// base, a span that begins past the first coding base, an alternate residue, any other
/// difference of terms and every other cell stay differences.
pub const START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR: &str =
    "start_codon_removed_without_five_prime_utr";
/// `start_codon_kept_without_five_prime_utr`: for an insertion or deletion
/// `_ins_del_start_altered` (Utils/VariationEffect.pm:1028-1075) writes the allele over its
/// span of the 5' UTR and the coding sequence (:1058) and reads the start as kept where the
/// UTR is unchanged and the three bases at the start position read `ATG` (:1062-1065). That
/// test stands inside `if ($utr)` (:1061), which is false for a transcript whose coding
/// sequence opens on the transcript's first base: there every deletion reads as an altered
/// start (:1069) and an insertion by the comparison of the sequence's tail (:1071). For an
/// allele that overlaps the start codon of such a transcript and leaves a start codon at the
/// start of the coding sequence, Ensembl VEP writes no `start_retained_variant` (:947-963),
/// and it writes `start_lost` (:851-904) where the allele is a frameshift (:870) or where
/// its peptide test holds (:878-882). vep-rs makes the kept-start test at the transcript's
/// first base and writes `start_retained_variant` and no `start_lost`: the row Ensembl VEP
/// writes for the same allele behind a 5' UTR. The shape
/// (`rules::kept_start::start_codon_kept_without_five_prime_utr_pair`): a Transcript row
/// whose two sets of terms both hold `frameshift_variant`, `inframe_insertion` or
/// `inframe_deletion` and are equal but for the candidate's `start_retained_variant`,
/// against the reference's `start_lost` or against no term of the reference; and a cDNA
/// span that is the row's CDS span, the same on both sides, and opens on one of the first
/// three coding bases. The class takes the Consequence of such a row and, as a class of
/// `TERM_PAIR_FOLLOWING_CLASSES`, the IMPACT each side's terms give and, where the
/// reference holds `start_lost`, the protein notation of a run under `--hgvs`: the
/// reference's holds the `?` of a lost start (TranscriptVariationAllele.pm:2161-2165) and
/// the candidate's, the notation of the change, holds none, both naming one protein. A row
/// behind a 5' UTR, a span that opens past the third coding base, any other difference of
/// terms, a protein notation that differs where the reference holds no `start_lost` and
/// every other cell stay differences.
pub const START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR: &str = "start_codon_kept_without_five_prime_utr";
/// `refseq_edited_model_deleted_span_term`: on a model whose alignment to the genome deleted
/// bases (`BAM_EDIT=OK` with a deletion edit), TranscriptMapper.pm:168-221 pairs no cDNA
/// position with the deleted bases (:213-218), so a sequence allele wholly inside them maps to
/// a Bio::EnsEMBL::Mapper::Gap: `cds_coords` and `cdna_coords` (BaseTranscriptVariation.pm
/// :478-528) hold the gap alone, `within_cds` and `within_cdna` (Utils/VariationEffect.pm
/// :636-693) need a Coordinate or a frameshift intron and find neither, no predicate is true,
/// and BaseVariationFeatureOverlapAllele.pm:285 assigns the default term, `intergenic_variant`
/// (Utils/Constants.pm:636-649, whose `within_feature => 0` condition the default never
/// consults), on the transcript's row although `_bvfo_preds` (:454-558) marked the allele
/// within the feature, in an exon and coding. The candidate writes the term of the transcript
/// region whose genomic span holds the allele: `coding_sequence_variant` inside the coding
/// sequence, the UTR term inside a UTR. The shape: a Transcript row of a sequence allele
/// carrying `BAM_EDIT` `OK` on both sides whose Consequence is `intergenic_variant` in the
/// reference and exactly one of the three region terms in the candidate, neither row with a
/// cDNA position and both with one IMPACT. A structural allele's row (a breakend mate within
/// --distance takes the same default), a row with a cDNA position (an allele reaching a mapped
/// base is `coding_sequence_variant` with a position on both engines) and any other candidate
/// term stay differences.
pub const REFSEQ_EDITED_MODEL_DELETED_SPAN_TERM: &str = "refseq_edited_model_deleted_span_term";
/// `refseq_edited_model_tail_term`: on a model whose edits inserted bases (`BAM_EDIT=OK`
/// with an insertion edit), the cached mapper (TranscriptMapper.pm:168-221, :207-211) numbers
/// every base after an insertion in the edited frame, so the model's last mapped bases carry
/// cDNA coordinates above the sum of its exon lengths; `within_cdna` (Utils/VariationEffect.pm
/// :670-694, line 680) requires the coordinate to be at most `$feat->length`, and
/// `Transcript::length` (the core API's Transcript.pm:1666-1682) adds the edits' lengths only
/// under `edits_enabled`, which the cached transcript object does not carry, so it returns the
/// unedited sum; `within_3_prime_utr` (:744-756) is false for as many final bases as the edits
/// inserted, no predicate is true, and BaseVariationFeatureOverlapAllele.pm:285 assigns the
/// default term, `intergenic_variant`, on a row that prints the cDNA position the same
/// coordinate gave. The candidate writes `3_prime_UTR_variant`. The shape: a Transcript row of
/// a sequence allele carrying `BAM_EDIT` `OK` on both sides whose Consequence is
/// `intergenic_variant` in the reference and `3_prime_UTR_variant` in the candidate, a cDNA
/// position on both and one IMPACT. A row without a cDNA position (the deleted-span class's
/// shape), a structural allele's row and any other candidate term stay differences; a
/// non-coding model's exon term reads the exon's genomic span, not the mapper, so it has no
/// such row.
pub const REFSEQ_EDITED_MODEL_TAIL_TERM: &str = "refseq_edited_model_tail_term";
/// `refseq_edited_model_equal_span_insertion_term`: on a model whose edits inserted bases
/// (`BAM_EDIT=OK` with an insertion edit), an insertion written where the model holds bases
/// the genome lacks stands on those bases, and its CDS span is theirs
/// (BaseTranscriptVariation.pm:258-264). BaseVariationFeatureOverlapAllele.pm:446 takes the
/// reference length from the genomic span, which holds no base for an insertion, so :384-386
/// mark the allele `increase_length`; the `include` rule of `missense_variant`
/// (Utils/Constants.pm:886-890) asks for `increase_length => 0`, and its predicate
/// (Utils/VariationEffect.pm:1085-1098) is not run. For an allele as long as its CDS span,
/// `frameshift` (:1556-1564) finds no length difference, `inframe_insertion` (:1111) finds
/// alternate codons as long as the reference codons, no predicate is true, and
/// BaseVariationFeatureOverlapAllele.pm:285 assigns the default term, `intergenic_variant`,
/// on a row that prints a cDNA, a CDS and a protein position. The candidate asks
/// `missense_variant` of an allele as long as its CDS span and writes it where the residues
/// change. The shape: a Transcript row of an insertion between two neighbouring positions
/// (`shared::inserted_between_neighbours`) carrying `BAM_EDIT` `OK` on both sides, its CDS
/// span one range of numbers, the same on both sides, that holds as many bases as the
/// allele, and its Consequence `intergenic_variant` alone in the reference and
/// `missense_variant` alone in the candidate. The class takes the Consequence of such a row
/// and, as a class of `TERM_PAIR_FOLLOWING_CLASSES`, the IMPACT each side's term gives. A
/// row without `BAM_EDIT` `OK`, an allele of another length than its CDS span, a row without
/// a CDS span, a deletion, a structural allele, any other term on either side and an IMPACT
/// that is not its side's term's stay differences. A default or tab row under an identifier
/// spells no alleles: the class is one of `RECORD_REFERENCE_CLASSES`, and a comparison given
/// the input reads the row's reference bases from its records, taking an insertion of any
/// length and refusing bases a record replaces. Where no record gives a row its bases, with
/// the input or without it, a two-base allele reads as two replaced bases and stays a
/// difference, and two bases replaced by a sequence of another length read as an insertion.
pub const REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM: &str =
    "refseq_edited_model_equal_span_insertion_term";
/// `refseq_edited_model_deleted_span_codons`: TranscriptVariationAllele.pm:872-874 splices the
/// reference allele's bases over the CDS span of a model carrying `_rna_edit` before reading
/// the reference codons (:877), `substr($cds, $cds_start - 1, $vf_nt_len) = $seq`, which
/// assumes the allele has the span's length. For an allele spanning bases the model's alignment
/// deleted, TranscriptVariation.pm:459 refuses the model's own bases as the reference allele
/// (the span maps to three segments) and the genomic bases stand, more of them than the CDS
/// span holds; the splice lengthens the CDS and the codon window reads the genomic bases from
/// the span's start, among them bases the model lacks, while the model's own final bases of the
/// span fall out of the window. The reference allele's amino acids are the translation of
/// that window (:703, :724-736), and `inframe_deletion` reads the same window
/// (Utils/VariationEffect.pm:1178, :1185): where it accepts by neither of its two tests,
/// `protein_altering_variant` (:375-392) stands in its place. The candidate reads the model's
/// codons. The shape: a Transcript row of a sequence allele carrying `BAM_EDIT` `OK` on both
/// sides whose USED_REF is longer than its CDS span. On such a row the class takes three
/// cells: the Codons whose alternate codons are equal and whose reference codons, of one
/// length, agree at every lower-case (unchanged) base and differ in an upper-case one, the
/// reference's upper-case bases spelling, on the transcript's strand, the leading bases of
/// USED_REF or of a GIVEN_REF of its length, which keeps the spliced allele where the
/// protein notation of a shifted deletion puts the shifted bases over the allele USED_REF
/// prints (TranscriptVariationAllele.pm:1785-1786, after OutputFactory.pm:1743 writes the
/// Codons);
/// the Amino_acids beside such Codons whose alternate halves are equal and whose reference
/// halves, of one length, differ, each the translation of its own side's reference codons
/// under the standard genetic code; and the Consequence beside such Codons that the
/// reference writes with `protein_altering_variant` where the candidate writes
/// `inframe_deletion`, every other term shared, where the two tests refuse the reference's
/// codons and accept the candidate's. The `most_severe_consequence` of a JSON record that
/// carries such a Consequence is left uncompared. A row whose allele has the span's length
/// (the splice replaces the model's bases with themselves), a differing alternate codon,
/// upper-case bases that are not the allele's, an amino acid that is not the translation of
/// its side's codons and any other pair of term sets stay differences.
pub const REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS: &str = "refseq_edited_model_deleted_span_codons";
/// `refseq_edited_model_intron_span_coding_reference`: for a span that has both ends on coding
/// bases of a model carrying `_rna_edit` and holds a whole intron, TranscriptVariation.pm:459
/// refuses the model's own bases as the reference allele (the span maps to three or more cDNA
/// segments) and the genomic bases stand, the intron's among them.
/// TranscriptVariationAllele.pm:872-874 write them over the CDS span and :877 reads the
/// reference codons from the span's first codon, so the reference half of Codons holds intron
/// bases and the reference half of Amino_acids is their translation. `stop_lost`
/// (Utils/VariationEffect.pm:1259-1261) accepts on a stop that translation holds,
/// `inframe_deletion` (:1178, :1185) accepts by neither of its two tests, so
/// `protein_altering_variant` stands in its place, and the protein notation of a row that
/// passes `stop_lost` is written with the extension (TranscriptVariationAllele.pm:1916,
/// :1926). The candidate reads the model's codons. The shape: a Transcript row of a sequence
/// allele carrying `BAM_EDIT` `OK` on both sides, its CDS span one range of numbers on both
/// sides, one USED_REF on both sides with more bases than that range, and
/// `splice_acceptor_variant`, `splice_donor_variant` and `intron_variant` among the terms of
/// both sides. On such a row the class takes four cells: a Consequence the reference writes
/// with `stop_lost` and `protein_altering_variant` where the candidate writes
/// `inframe_deletion`, every other term shared; the Codons that
/// `refseq_edited_model_deleted_span_codons` describes; the Amino_acids beside such Codons
/// whose alternate halves are equal and whose reference halves, of one length, differ, each
/// the translation of its own side's reference codons under the standard genetic code; and
/// an HGVSp beside such a Consequence that the reference writes as the candidate's followed
/// by `extTer?`. A row without the shape, any other pair of term sets, any other notation
/// and every other cell stay differences. The class stands before
/// `refseq_edited_model_deleted_span_codons` in `FIELD_CLASSES`, so a comparison that names
/// both counts the Codons of a row of this shape under this class.
pub const REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE: &str =
    "refseq_edited_model_intron_span_coding_reference";
/// `refseq_edited_model_inserted_span_codons`: an insertion written where a model carrying
/// `_rna_edit` holds bases the genome lacks has those bases as its coding span
/// (Mapper.pm:485-558, BaseTranscriptVariation.pm:263-264) and an empty reference allele
/// (TranscriptVariation.pm:451). TranscriptVariationAllele.pm:872-873 write that allele over
/// the span and :877 reads the reference codons from a coding sequence that has lost the
/// span: the other bases of the row's residues, then bases from past the row's last
/// residue. The reference half of Amino_acids is their translation, a last base or two
/// read as `X` (:724-736, :774-776). The candidate reads the model's
/// bases over the row's residues. The shape: a Transcript row of an allele of bases
/// carrying `BAM_EDIT` `OK` on both sides, its USED_REF empty on both, its coding span and
/// its residues the same numbers on both. On such a row the class takes two cells, each
/// derived from the coding sequence the candidate JSON cache holds for the transcript
/// (`DownstreamSequences`, `--json-cache`): the Codons whose alternate halves are equal,
/// the candidate's reference half the sequence's bases over the row's residues and the
/// reference's as many bases read from the sequence without the span; and the Amino_acids
/// beside such Codons whose alternate halves are equal and whose reference halves differ,
/// each the translation of its own side's reference codons under the transcript's genetic
/// code as that engine reads it. Such Amino_acids on a transcript under a codon table no
/// table is held for (`genetic_codes::held_tables`) end the run with an error that names
/// the transcript: they are read under no other code. Without the cache, on a transcript
/// it lacks and on two readings that are equal the class takes nothing; any other pair of
/// values and every other cell stay differences.
pub const REFSEQ_EDITED_MODEL_INSERTED_SPAN_CODONS: &str =
    "refseq_edited_model_inserted_span_codons";
/// `frequency_filter_zero_value`: BaseCacheVariation.pm:419 at release/116.2 reads an allele's
/// frequency with `|| 'NA'`, so a frequency of exactly 0, Perl-false, reads `NA` and fails
/// --check_frequency and --filter_common; the candidate judges the zero, as Ensembl VEP 116.3
/// does. Two shapes: a `FREQS` whose items agree but for a value the reference reads `NA`
/// and the candidate `0`, and under `--freq_filter include` the rows of a record the
/// candidate alone keeps, every `FREQS` item of which carries the value `0`, with the
/// record's VCF line and JSON object. A JSON consequence entry carries no `FREQS`, so in that
/// format the row shape applies from a manifest's candidate-only tuples alone. Ensembl VEP
/// 116.3 corrects the defect, so the class applies under an earlier reference release alone
/// (`corrected_under`).
pub const FREQUENCY_FILTER_ZERO_VALUE: &str = "frequency_filter_zero_value";
/// `merged_same_stable_id_collapse`: under a RefSeq or merged source AnnotationType/Transcript.pm
/// :313-343 at release/116.3 keeps one transcript per stable identifier among the models a
/// buffer's regions load (the `ensembl`-sourced copy, else the lowest dbID), where the method's
/// first pass and its description uniquify by dbID; a model the cache carries at two loci keeps
/// the copy so chosen and loses every row at the other, whether or not a record overlaps the
/// kept copy. The candidate annotates both copies. The shape: a Transcript entry the candidate
/// alone writes on a record the reference wrote a Transcript row for, whose SOURCE, where the
/// run writes one, is RefSeq (a merged cache's Ensembl-sourced models are unique per locus; a
/// RefSeq cache writes no SOURCE), of a model the reference wrote no row for on that record; the
/// JSON record's most_severe_consequence is left uncompared as for a documented key. A record
/// whose every overlapping model is a lost copy gets an intergenic row from the reference, a
/// shape outside the class, as are another allele's row of a kept model and a duplicated row.
pub const MERGED_SAME_STABLE_ID_COLLAPSE: &str = "merged_same_stable_id_collapse";
/// `use_transcript_ref_cache_cdna`: under --use_transcript_ref, which a cache corrected against
/// an alignment switches on, TranscriptVariation.pm:462-463 cut a row's reference allele from
/// the transcript's `spliced_seq` at the row's cDNA start. A model whose correction failed has
/// no cached sequence (AnnotationType/Transcript.pm:617 caches one for a corrected model
/// alone), the cache writer deletes the key that switches a transcript's edits on
/// (Pipeline/DumpVEP/Dumper/Core.pm:164), and Transcript.pm:866 of Ensembl's core API then
/// joins the exons without the model's own `_rna_edit`, while the cached mapper numbers the
/// cDNA with it. Past that edit the cut reads the base one place on along the transcript, or
/// one place back where the edit removes a position, and TranscriptVariationAllele.pm:872-877
/// read the reference codon with that base in it. The candidate reads the cache's cDNA, which
/// holds the edit, and its reference allele is the record's. The shape: a Transcript row of a
/// one-base allele carrying `BAM_EDIT` `FAILED` on both sides, the GIVEN_REF one base and the
/// same on both, the candidate's USED_REF that base and the reference's another. On such a
/// row the class takes seven cells: the USED_REF; the Codons of one codon a side whose
/// alternate halves are equal and whose reference halves differ in one upper-case base, the
/// reference's its USED_REF and the candidate's the GIVEN_REF on the transcript's strand; the
/// Amino_acids each side's Codons translate to under the standard genetic code; a Consequence
/// whose two sets differ in the one term each side's two residues give (`missense_variant`,
/// `synonymous_variant`, `stop_gained`, `stop_lost` or `stop_retained_variant`), every other
/// term shared; the IMPACT of each side's terms beside such a Consequence; an HGVSc of one
/// substitution a side at one position, the base changed the reference's USED_REF on its
/// side and the GIVEN_REF on the candidate's (TranscriptVariationAllele.pm:1510-1513 write
/// the reference allele into the notation); and an HGVSp of one residue a side at one
/// position, each side's the change of the reference residue of its own Codons to the
/// alternate one, the two reference residues differing and neither a stop (:1793 takes the
/// residue from the reference codon). The `most_severe_consequence` of a JSON record that
/// holds a row whose Consequence is taken is left uncompared: OutputFactory/JSON.pm:388
/// writes the lowest-ranked term of the record's rows. A row without the shape, a reference
/// allele of two or more bases, any other term, a notation at another position, an HGVSp
/// beside one reference residue on both sides or a stop on either, and every other cell
/// stay differences. The entry's other rows, written where the FASTA holds other bases than
/// the cache, are not told from the rows alone: they stay differences too, but for one that
/// has this shape.
pub const USE_TRANSCRIPT_REF_CACHE_CDNA: &str = "use_transcript_ref_cache_cdna";
/// `refseq_failed_edit_model_hgvsc_offset`: a model whose alignment correction failed keeps
/// the `_rna_edit` it came with, which marks a base that translation reads twice or passes
/// over, and its transcript is the genome's exons joined. The cached numbering counts the
/// edit; the candidate counts every HGVSc position along the genome's exons. Ensembl VEP
/// 116.3 numbers a one-base substitution that has a CDS position at the row's CDS start,
/// the cached number (TranscriptVariationAllele.pm:1519-1522), and any other notation along
/// the exons moved on by `get_misalignment_offset` (:1516-1517, :1524-1525), which adds the
/// length of the text of an edit of three fields (:1617) and takes off the span of an edit
/// of two (:1614), between the cached coding positions. The shape: an HGVSc on a row
/// carrying `BAM_EDIT` `FAILED` and one pair of cDNA positions on both sides, whose
/// transcript carries, in the candidate's cache, an edit that changes the cDNA's length
/// before the row's first cDNA position; the two values name one sequence in coding
/// numbering with as many coordinates, none with an intron offset. An edit stands before a
/// row whose first position is past the positions the edit adds, or at or past the first
/// position of an edit that takes its span out; a row that holds a position an edit adds,
/// or runs across the bases one takes out, is not the shape. On a row with one CDS
/// position, a substitution of one base stands at that CDS position in the reference and at
/// the count of the row's cDNA position in the candidate. Each coordinate of any other
/// notation is, in the reference, the candidate's count plus the offset summed for the
/// edits before the row, read between the cached coding positions so that a coordinate
/// carried over the stop codon is read too, and the candidate's count lies within the
/// row's `HGVS_OFFSET` and its allele's length of the row's positions. The text after the coordinates is the same on
/// both sides, or the base changed differs as `use_transcript_ref_cache_cdna` reads it
/// where the comparison names that class as well. The class reads the cache and accepts
/// nothing without one. Any other distance, a row that holds the edit's own position, a row
/// before the edit and the sides exchanged stay differences. So does, on a row with one CDS
/// position, a one-base substitution whose record holds an allele that is not one base,
/// where the offset's number is not that CDS position: :1519 tests the record, and Ensembl
/// VEP numbers such an allele by the offset.
pub const REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET: &str = "refseq_failed_edit_model_hgvsc_offset";
/// `hgnc_id_copy`: AnnotationType/Transcript.pm:290 at release/116.3 records the HGNC id of
/// every transcript carrying one under its gene symbol and :305 lends it to every same-symbol
/// transcript among the features fetched for one buffer, so a transcript whose gene is named
/// outside HGNC (EntrezGene, RFAM; every RefSeq transcript of a merged cache) prints an id in a
/// buffer that also loaded the same-named HGNC gene and nothing in one that did not. The
/// candidate lends the id once, as the cache is built. The shape: the HGNC_ID of a Transcript
/// entry whose SYMBOL_SOURCE is not HGNC, one side empty and the other an id (`HGNC:10528`,
/// or `10528` in a cache before the prefix), in either direction. Only a run that writes the
/// field compares it: `--symbol`, or the VCF format, whose option set implies it
/// (`VCF_IMPLIED_COLUMNS`).
pub const HGNC_ID_COPY: &str = "hgnc_id_copy";
/// The two field classes of the Downstream plugin's `DownstreamProtein`. Downstream.pm:108-110
/// write the allele over the change's first to last coding base, :112-113 round the first base's
/// coordinate, counted from 1, down to a multiple of three, :115-119 read the changed coding
/// sequence from that offset, counted from 0, and follow it with the 3' UTR, and :139-142
/// translate the result and cut it at its first stop. For a deletion, with or without inserted
/// bases, whose first base is a codon's third base the offset is the base after it, so the
/// translation starts one codon after the first codon the change touches; the candidate starts
/// at the touched codon. The rows do not state the sequences, so both classes read them from the
/// candidate JSON cache (`DownstreamSequences`, `--json-cache`) and recompute both values
/// (`downstream_class`); without the cache they accept nothing. The shape the two share: a
/// DownstreamProtein both sides carry (an empty value is a value, which the plugin writes when
/// the first codon read is a stop; `-` or no key is none) on a row whose two entries have one
/// consequence set, or two that differ by the reference's `start_lost` as the registry's
/// `start_cooccurrence_swap` admits (`registry_class`), as the row rule of
/// `start_codon_kept_without_five_prime_utr` takes them
/// (`start_codon_kept_without_five_prime_utr_pair`) or as a class of a coding end takes
/// them on the transcript's record in the cache (`coding_end_class`), with
/// `frameshift_variant` and no splice term (:95-99), one CDS_position (or the two of a later
/// allele's row,
/// `later_allele_prints_shifted_pair`) and one ProteinLengthChange; the cache holds the
/// row's transcript with a
/// `translateable_seq`; and on the record's own coordinates, the CDS_position as printed, the
/// allele read on the transcript's strand (`-` as no bases, any other reverse-complemented
/// where the cache's `strand`, or the row's STRAND for a transcript the cache gives none, is
/// -1), the first base is a codon's third base, the reference's value is the changed sequence
/// translated under the transcript's codon table (:130-133, `genetic_codes`) from the offset
/// with the UTR and the candidate's the same translated from three bases before it. Under
/// --hgvs the plugin holds the coordinates the protein notation shifted
/// (`plugin_hgvs_shifted_cds_coordinates` below): a row is one of these two classes' where
/// the plugin's reading of the shifted coordinates gives the value its reading of the record's
/// own gives, and that class's where the two differ. A value one side lacks, another value on
/// either side, a first base that is a codon's first or second, a differing
/// ProteinLengthChange and a transcript the cache lacks stay differences. A row that asks
/// for a transcript under a codon table no table is held for (`genetic_codes::held_tables`)
/// ends the run with an error that names the transcript: it is read under no other code.
///
/// `downstream_utr_frame_past_coding_end`: the offset lies past the changed sequence's end, a
/// deletion without inserted bases (allele `-`) that runs to the last base of the
/// `translateable_seq`. The changed sequence is one base shorter than the offset, the
/// `substr` of :115 returns nothing, and the 3' UTR is translated from its first base, one
/// base out of the changed sequence's frame, where the candidate translates the two coding
/// bases before the deletion with the UTR. The transcript has a 3' UTR: without one both
/// values are empty.
///
/// `downstream_start_past_own_codon`: the offset lies inside the changed sequence or at its
/// end, every other deletion and every deletion with inserted bases of the shared shape; the
/// reference's value is the candidate's less the touched codon's residue, or the residues past
/// a stop the touched codon spells where the candidate's value is empty. A transcript without
/// a 3' UTR is read with an empty one. An insertion between two bases the model keeps adjacent
/// is never the class's: its CDS_position prints the two bases around it, and read as a
/// replacement of those two bases its recomputed values are neither engine's. An insertion on
/// a seam of the model, two mapper pairs that touch on the genome with model bases between
/// them, prints those bases as a span of its own (Mapper.pm:504-558), which is read as any
/// other span. The two classes share no row.
pub const DOWNSTREAM_UTR_FRAME_PAST_CODING_END: &str = "downstream_utr_frame_past_coding_end";
pub const DOWNSTREAM_START_PAST_OWN_CODON: &str = "downstream_start_past_own_codon";
/// `plugin_hgvs_shifted_cds_coordinates`: under --hgvs a plugin reads the CDS coordinates the
/// protein notation left shifted on the row. OutputFactory.pm:1746 prints CDS_position from the
/// record's own span and :1769 calls `hgvs_protein`, whose `codon` deletes the row's kept
/// CDS pair and reads it again at the 3'-shifted position (TranscriptVariationAllele.pm:844-858);
/// `_get_alternate_cds` reads it there once more and leaves it there (:2406-2410), and
/// `shift_feature_seqs` rotates the allele's sequence for the shifted position (:1330-1345,
/// called at :1778). `hgvs_protein` returns before all of it when the shifted span has no
/// translation coordinates (:1730-1738). Nothing clears the pair before OutputFactory.pm:398
/// runs the plugins, so Downstream.pm:108 and :110 read the shifted pair and the rotated
/// allele, and NMD.pm:88, :98 and :103 the shifted pair (:137 reads the record's own genomic
/// end). Without --hgvs the same plugins read the record's own pair, which the candidate reads
/// with and without the flag. The class reads the candidate JSON cache as the two Downstream
/// classes do (`DownstreamSequences`, `--json-cache`), with the transcript's strand, exons,
/// coding bounds and mapper pairs, and accepts nothing without it.
///
/// The two pairs of a row whose reference entry carries a non-zero HGVS_OFFSET: the own pair
/// is the CDS_position, on a later allele's row whose reference entry prints the shifted pair
/// the candidate's (`later_allele_prints_shifted_pair`: the alleles of a record share one
/// TranscriptVariation, and OutputFactory.pm:1746 prints on a later allele's row the pair an
/// earlier allele's notation left there, the span moved by the row's own offset or by the
/// offset of the last earlier allele's row that carries one, `earlier_allele_shift`); the
/// shifted pair is the record's genomic span moved by the offset as
/// printed (signed along the genome) and mapped through the mapper pairs
/// (`TranscriptModel.cds_pair`: BaseTranscriptVariation.pm:252-264, and Mapper.pm:504-558 for
/// an insertion, whose pair runs from the base after its 5' flank to the base before its 3'
/// flank). The record's span is the Location of a text entry, which must map to the row's
/// CDS_position; a VCF or JSON entry carries none, so its span is the one the CDS pair maps
/// back to, and it has none where the pair holds a gap or a base outside the mapper pairs or
/// is an insertion's between two bases the genome does not keep adjacent. The shift reaches
/// the plugin where both ends of the shifted span are CDS coordinates; the allele it reads
/// there is the Allele turned left one base as many times as the offset's absolute value on a
/// forward-strand transcript, as its length less that value on a reverse-strand one (not at
/// all when that is not positive), then read on the transcript's strand. An insertion and a
/// two-base replacement print the same two adjacent numbers in Location and CDS_position, so a
/// sequence allele at such a span is read both ways, a reading counting only where it
/// recomputes both entries' values.
///
/// The Downstream shape: a DownstreamProtein or a ProteinLengthChange both entries carry, on a
/// row of the two Downstream classes' consequence set and CDS_position, where the candidate's
/// value is the candidate's rule on the own pair with the Allele, the reference's is
/// Downstream.pm:104-149 on the shifted pair with the rotated allele, and that differs from
/// Downstream.pm:104-149 on the own pair, the value Ensembl VEP writes without the flag. A row
/// whose two readings of Downstream.pm give one value is left to the two classes above, so the
/// three share no row. Where the candidate's value is not that own-pair value either, the own
/// pair has the shape of one of the two classes above and the start rule changes the value as
/// well as the flag: the row is this class's only where that class is declared beside it, and
/// with either of the two declared alone it stays a difference. The NMD shape: an
/// `NMD_escaping_variant` one entry carries and the
/// other lacks, on a row whose two entries have one consequence set, or two that differ as
/// `stop_retained_without_stop_codon_test_pair` admits, each naming one of the plugin's four
/// terms (NMD.pm:65), and one CDS_position or a later allele's two, where the reference carries it
/// exactly when the rule of
/// NMD.pm:98-109 holds on the shifted pair (both ends coordinates, and the CDS end at most 101
/// or the record's genomic end in the last exon or within 51 bases of the 3' end of the second
/// to last, `variant_exon_check`, :132-176), the candidate exactly when it holds on the own
/// pair, and the two differ. :90-91 count the one array reference `get_all_Introns` returns,
/// so the rule of :109 for a transcript without introns never holds, and a value the candidate
/// writes by it alone is outside the class. A later allele's row without an HGVS_OFFSET whose
/// earlier allele's row carries one is the Downstream shape with the pair that allele's
/// notation left, the pair the reference prints, in the shifted pair's place and the Allele
/// as written in the rotated allele's. The same values on any other row without an HGVS_OFFSET,
/// an NMD on a row without one, with a shift that does not reach the plugin or with a shifted pair
/// that gives the own
/// pair's value, a value planted on either side, a Downstream value one side lacks, an NMD
/// both entries carry or both lack, and a transcript the cache lacks or holds without its
/// mapper pairs stay differences.
pub const PLUGIN_HGVS_SHIFTED_CDS_COORDINATES: &str = "plugin_hgvs_shifted_cds_coordinates";
/// `plugin_value_on_reference_skipped_record`: a field of DosageSensitivity or UTRAnnotator
/// that the candidate alone writes on a row of a record the reference dropped
/// (`reference_skipped_record`). Ensembl VEP writes no row for a record it marks skipped
/// (Parser.pm:491-501, :821-827), so no plugin runs on it; the candidate annotates the record
/// and its plugins write on those rows as on any structural row. The shape is the base class's
/// own: a Location at which the candidate wrote a structural row and the reference none, no
/// candidate row there naming a record the reference wrote (`ReferenceRecordScope` and
/// `reference_named_locations` of `tools/vep-outputs-compare`). A row of such a record has no
/// reference entry, so this comparison sets the row aside whole under the structural-variant
/// class and `field_class_accepts` takes no pair for the class; the plugin tuple scorer,
/// which scores a value and not a row, decides it from the rows of the value's Location. A
/// value on the reference's side, a value at a Location the reference wrote a row for and a
/// value of a tuple-file engine, which states no rows, stay differences.
pub const PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD: &str =
    "plugin_value_on_reference_skipped_record";
/// The three classes of the structural-variant comparator (`tools/vep-sv-compare`), built in
/// so that a structural-variant case sets aside what the
/// adjusted structural-variant figures set aside: each is decided by that comparator's own
/// predicate over the entries read as its tuples (Location, Allele, Feature, Feature_type,
/// Consequence set), a VCF or JSON entry's Location the one its input record gives
/// (`vep_rs_record_rows`). `reference_skipped_record` (`filter_reference_skipped_records`):
/// a row the candidate alone writes on a record Ensembl VEP drops before its input buffer, an
/// unsupported type (Parser/VCF.pm:477-481: `<NON_REF>`, `<CPX>`, `<CTX>`), a span over
/// `--max_sv_size` (Parser.pm:491-501), an incomplete deletion or an invalid span
/// (`perl_skip_route`), and wrote no row for; the record's VCF line and JSON object are
/// counted, not compared, as a manifest's are. It reads the input VCF in every format.
/// `breakend_mate_local_read` (`filter_breakend_mate_local_read_pairs`): a key both engines
/// write under a bracket-mate allele whose consequence sets are the local-coordinate and the
/// mate-coordinate derivations over the cache's transcript model
/// (BaseVariationFeatureOverlapAllele.pm:257,273 evaluate the mate row's predicates on the local
/// variation feature); the pair is skipped whole, as a registry pair is, so its DISTANCE, cDNA,
/// CDS and Protein positions follow. `bnd_own_chromosome_rows_lost`
/// (`filter_bnd_own_rows_lost`): a Transcript row the candidate alone writes under an own-end
/// allele at a Location where the reference wrote mate rows naming another chromosome and no
/// own-end row (AnnotationType/Transcript.pm:122 names the record's slice by the first
/// transcript reached, :125-143 then fail every own-chromosome transcript), the transcript on
/// the Location's chromosome within 5,000 bases of the breakpoint. The two breakend classes
/// read the candidate's JSON cache (`--json-cache`), as the structural-variant comparator does.
pub const BREAKEND_MATE_LOCAL_READ: &str = "breakend_mate_local_read";
pub const BND_OWN_CHROMOSOME_ROWS_LOST: &str = "bnd_own_chromosome_rows_lost";
/// The classes documented on the divergence pages that a manifest alone applies (a corpus's
/// `divergences`, `vep_rs_only_tuples` or `field_divergences` under the class), each
/// a shape no predicate decides from the two outputs: the fields a TranscriptAnnotator table
/// names carry the user's prefix; a later allele of a multi-allelic insertion record annotated
/// at an earlier allele's shifted span (OutputFactory.pm:306-322, TranscriptVariationAllele.pm
/// :502-527) and an indel within 1,000 bases of a sequence region's edge shifted short of the
/// bases the sequence admits (TranscriptVariationAllele.pm:427-437, :462) change the notation
/// of a row whose coordinates the outputs do not carry; and a matched allele's frequency read
/// after the input allele's own key (BaseCacheVariation.pm:419-429 at release/116.2) decides
/// which records --check_frequency keeps by a cache value the outputs do not print. The others
/// are shapes of a flag run:
/// under --minimal a structural record with two symbolic ALTs is trimmed as two sequences
/// (Parser.pm:852, Utils/Sequence.pm:965-1038), so its rows move to a span that starts after
/// the characters the two names share and ends before the `>` they share, marked MINIMISED,
/// and a record with three or more keeps its span and loses the rows of every allele past the
/// second (InputBuffer.pm:638-718, OutputFactory.pm:2320-2327): the reference's keys and the
/// candidate's are both documented, the reference's with no candidate set; under
/// --shift_3prime 1 the alleles of one record share one transcript variation whose span
/// predicates the first allele stores at its own shifted span and whose coding coordinates
/// are read back cached (BaseVariationFeatureOverlapAllele.pm:363, BaseTranscriptVariation.pm
/// :142), so a later allele shifting by another length, or the unshifted reference row of
/// --process_ref_homs, takes another span's terms; and the start-codon test of the same flag
/// adds the shift to the unshifted cDNA start a deletion beginning in an intron lacks
/// (Utils/VariationEffect.pm:977-981), so such a deletion reaching past the start codon gains
/// start_retained_variant, or start_lost with a lost-start HGVSp on a transcript without a 5'
/// UTR (:1069). The outputs carry neither the record's other alleles nor the unshifted span.
/// Two of the classes are decided without a manifest as well, each on the rows its own comment
/// names: `hgvs_shift_shared_across_alleles` (above `SHARED_SHIFT_FIELDS` of `consts`) and
/// `hgvs_shift_sequence_edge_window` (in `rules::sequence_edges`).
pub const TRANSCRIPT_ANNOTATOR_SIBLING_ALLELE_ROWS: &str =
    "transcript_annotator_sibling_allele_rows";
pub const HGVS_SHIFT_SHARED_ACROSS_ALLELES: &str = "hgvs_shift_shared_across_alleles";
/// `blosum62_hgvs_shift_shared_across_alleles`: the BLOSUM62 the reference alone writes on a
/// later allele's row of `hgvs_shift_shared_across_alleles`. Blosum62.pm:110 scores a row
/// whose peptide pair is one residue over one residue, with the matrix cell of the two
/// (:50-71, :73, :84-88). On a later allele's row the reference reads the reference residue at
/// the span an earlier allele's notation moved, so an insertion between two codons, whose own
/// pair has no reference residue and gets no score, carries one residue over one and a score.
/// The shape: the candidate writes no BLOSUM62 and its own `Amino_acids` fails the plugin's
/// test; the row's CDS and protein positions are the base class's
/// (`later_allele_shared_shift`); the reference's `Amino_acids` is one residue over one,
/// the first the translation of the moved pair's codon and the second the candidate's own
/// alternate residue; and the reference's value is the matrix cell of that pair. A value that
/// is not that cell, a value on the candidate's side, a candidate pair the plugin scores, a
/// first allele's row, a record whose first allele carries no offset and a run without the
/// cache stay differences.
pub const BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES: &str =
    "blosum62_hgvs_shift_shared_across_alleles";
pub const HGVS_SHIFT_SEQUENCE_EDGE_WINDOW: &str = "hgvs_shift_sequence_edge_window";
/// `frequency_filter_allele_key_order`: Ensembl VEP releases before 116.3 read the input
/// allele's own key before the matched allele's (BaseCacheVariation.pm:419-429 at
/// release/116.2). Ensembl VEP 116.3 corrects the defect, so the class applies under an
/// earlier reference release alone (`corrected_under`).
pub const FREQUENCY_FILTER_ALLELE_KEY_ORDER: &str = "frequency_filter_allele_key_order";
pub const MINIMAL_SYMBOLIC_ALLELES: &str = "minimal_symbolic_alleles";
pub const SHIFT_3PRIME_SHARED_PREDICATES: &str = "shift_3prime_shared_predicates";
pub const SHIFT_START_CODON_TEST: &str = "shift_start_codon_test";
/// `transcript_annotator_insertion_window`: TranscriptAnnotator.pm:218 reads the table rows
/// from `start - 2` to `end`, a two-base window for an insertion (`end` = `start` - 1)
/// that never reaches a row written at the insertion's start; the candidate reads from two
/// bases before the lower coordinate to the higher. The shape, the table's fields the candidate
/// alone fills on an insertion's row, is documented from a corpus manifest.
pub const TRANSCRIPT_ANNOTATOR_INSERTION_WINDOW: &str = "transcript_annotator_insertion_window";
/// `plugin_value_follows_splice_term_past_the_record_span`: a field of Downstream, NMD or
/// LOFTEE that one engine alone writes on a row `splice_term_past_the_record_span_pair`
/// takes. The three plugins run on a row's terms. LOFTEE sets its splice gate on
/// `splice_acceptor_variant` or `splice_donor_variant` (LoF.pm:167) and returns without a
/// confidence for a row that passes neither that gate, the one on `stop_gained` and
/// `frameshift_variant` nor its own splice prediction (:222-232); NMD runs on one of the
/// same four terms (NMD.pm:65 and :83). Where the donor or acceptor term the reference alone
/// holds is its row's only such term, the fields of both stand on the reference's side alone.
/// Downstream writes no value beside a term that holds `splice` (Downstream.pm:97-99):
/// where the reference holds `frameshift_variant` beside the term and the candidate's row
/// holds no splice term, its fields stand on the candidate's side alone. `intron_variant` is
/// a term of none of those tests, so the intron term of the class moves no field. The shape:
/// the row's two entries are a pair of that predicate, read with the record's alleles where
/// an entry carries them (`Entry.reference`) and without them where none does, and the field
/// is a cell `term_pair_following_cell` reads from the two term sets. A row the predicate
/// refuses, a field on a side whose terms keep the plugin out, a field both sides carry and
/// every field no term decides stay differences.
pub const PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN: &str =
    "plugin_value_follows_splice_term_past_the_record_span";
/// `plugin_value_follows_stop_codon_term_pair`: a field of a plugin that runs by a row's
/// terms, which one engine alone writes on a row whose two Consequence sets differ as a stop
/// codon class or a class of a coding end admits. Downstream, NMD and LOFTEE run on a
/// frameshift or a truncating term (Downstream.pm:97-99, NMD.pm:65 and :83, LoF.pm:166-167),
/// REVEL and AlphaMissense on `missense_variant` (REVEL.pm:162, AlphaMissense.pm:242), and
/// ClinPred and dbNSFP on `missense_variant`, `stop_lost`, `stop_gained` or `start_lost`
/// (ClinPred.pm:69 and :106; dbNSFP.pm:157 and :290-291, its fields read under the prefix
/// `dbNSFP_`), so where the terms differ the field stands on the side whose terms let the
/// plugin run and is left out on the other. The shape: the row's two entries are a pair of
/// one predicate of `STOP_CODON_PLUGIN_VALUE_PAIRS`, of which the five of a coding end take
/// a pair beside the candidate's cache alone, and the field is a cell
/// `term_pair_following_cell` reads from the two term sets. A row of any other pair, a field
/// on a side whose terms keep the plugin out, a field on one side where both sides' terms let
/// the plugin run, a field both sides carry and every field no term decides stay
/// differences.
pub const PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR: &str =
    "plugin_value_follows_stop_codon_term_pair";
/// `total_length_json_range_end`: under `--total_length`, OutputFactory.pm:1737, :1747 and
/// :1751 at release/116.3 append `/` and the sequence's length to a position column after its
/// range is written (`484-487/2674`), and OutputFactory/JSON.pm:335-337 split the column at
/// the hyphen and write the second part under the end key only where it is all digits or holds
/// a `?`, else the first part. The last position of a range that carries a total is never
/// written: the end key repeats the start, and a range whose start is unknown (`?-487/2674`)
/// loses both keys (:340-341). The candidate writes the last position with its total under the
/// end key, the form both engines give a single position under the flag (`"487/2674"`). The
/// shape, on `cdna_end`, `cds_end` or `protein_end` of a JSON consequence object: the
/// candidate's value is digits, `/` and digits; either both entries carry the start key with
/// one value of digits, the reference's end is that value and the candidate's last position
/// is above it, or neither entry carries the start key and the reference carries no end. A
/// start that differs, an end below or at the start, a candidate end without a total and an
/// end the reference writes beside no start stay differences.
pub const TOTAL_LENGTH_JSON_RANGE_END: &str = "total_length_json_range_end";
/// `allele_number_repeated_or_reference_equal_alt`: under `--allele_number`,
/// VariationFeatureOverlap.pm:459-467 at release/116.3 drop from a record's allele list an
/// allele string seen before, the reference allele's being the first seen, and :477-491 number
/// what is left: the reference allele's own entry is passed over and every other entry takes
/// the next number, `*` and `<DEL:*>` among them. An alternate written after a dropped entry
/// is so numbered below its place in the record. The candidate numbers an allele by its place
/// in the record's written allele list, 1 for the first written alternate. The shape, on
/// `ALLELE_NUM` or `allele_num` (`rules::allele_number::allele_number_place_pair`): both
/// entries carry the written list on the row itself, under `UPLOADED_ALLELE` or
/// `uploaded_allele`, as one text of three or more parts between its `/`, an empty part read
/// as `-`; the row's allele is a written alternate, or an alternate as
/// Parser/VCF.pm:294-335 spell it on the rows of an indel record, a record with an alternate
/// of another length than the reference allele's whose reference allele and alternates that
/// hold no `*` begin with one character: without that character, `-` for an allele left
/// empty, an alternate that holds a `*` as written; the candidate's number is the place of
/// the first such alternate among the written alternates, counted from 1; the reference's
/// number is the count of the written alternates up to and including that place that equal
/// neither the reference allele nor an alternate written before them; and the two numbers
/// differ. A row without the list, a list the two entries write differently, an allele the
/// list holds under neither spelling, an allele that is an alternate equal to the reference
/// allele, a row whose allele fits two alternates at different places, one as written and
/// one as the parser spells it, and a number that is neither the place nor the count stay
/// differences.
pub const ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT: &str =
    "allele_number_repeated_or_reference_equal_alt";
/// `colocated_list_cut_at_cache_cap`: the Ensembl cache dump ends a stored `pubmed` or
/// `var_synonyms` value of a known variant at 1,024 bytes, inside an id as readily as after a
/// separator, and Ensembl VEP prints the stored value as it stands (OutputFactory.pm:1090,
/// :1126 and :1141, OutputFactory/JSON.pm:452-457 and :469 at release/116.3), the cut item as
/// an id. The candidate reads a value of exactly that length ending in none of its separators
/// (`,` in `pubmed`; `--`, `::` and `,` in `var_synonyms`) as the bytes before its
/// last separator, a `var_synonyms` group whose only id is the cut one going with its source
/// name and the `--` before that, and writes no field for a value left with no id. The shape:
/// the reference's list is the candidate's plus one last item, or plus one last group of one
/// id, and the reference's own text shows the stored value at exactly 1,024 bytes and ending in
/// no separator. A shorter list, a list ending in a separator, a pair whose earlier items
/// differ and a candidate that keeps the source name of such a group, bare in a row or with an
/// empty array in JSON, stay differences.
///
/// The default, tab and VCF formats join the lists of a row's known variants, `PUBMED` ids
/// with `,` and `VAR_SYNONYMS` stored values with `--` (VCF.pm:401 writes each `,` as
/// `&`), so a stored value is a run of the row's value. The candidate's value is the
/// reference's without the last item of one or more runs, each exactly 1,024 bytes from the
/// value's start or the end of a joining separator to the value's end or the start of one; a
/// run left with no id goes whole with one joining separator, and the field with it when the
/// row joins nothing else. A run of that length spanning the values of two known variants
/// reads as one stored value: the row's text does not tell the two apart.
///
/// The JSON format holds one `colocated_variants` entry per known variant. A `pubmed` array
/// holds the stored value's items, so its items joined with `,` are the stored value. A
/// `var_synonyms` object is rebuilt as its groups, `source::id,id`, joined with `--`; the
/// length does not depend on the groups' order, which the object loses, and the cut item's
/// group is the one whose last id the candidate lacks, or the source the candidate lacks, a
/// name with one id or with no `::`. The object does not show which group the stored value
/// ends in: a candidate lacking the last id of any one group, or any one group of one id, is
/// read as the shape on the rebuilt length alone, where a row's text holds it to the last
/// group of a run. Three shapes do not show the stored length there and stay
/// differences: a second group listing no id (its stored text is the name with or without
/// `::`), a source name the stored value repeats (JSON.pm:457 keeps the ids of its last
/// group alone) and an id whose JSON number is shorter than its stored text (`numberify`,
/// Utils.pm:241, prints `999001.0010` as `999001.001`). One such id is read: a cut item
/// printed as a number with no decimal part stands for its digits alone or for its digits
/// and a closing point, so the rebuilt value is held to 1,024 bytes with the point or
/// without it.
///
/// An accepted field is counted under the class once per entry: every row of the record in the
/// default, tab and VCF formats, the record's object in JSON. The class is no member of
/// `FIELD_CLASSES` or of any other class list of this module: `compare_entries` applies it
/// when `Documented.masks` names it.
pub const COLOCATED_LIST_CUT_AT_CACHE_CAP: &str = "colocated_list_cut_at_cache_cap";

/// Every class constant: its name and its value.
pub const NAMES: &[(&str, &str)] = &[
    ("REFERENCE_SKIPPED_RECORD", REFERENCE_SKIPPED_RECORD),
    ("REGULATORY_FEATURE_ROWS", REGULATORY_FEATURE_ROWS),
    (
        "MOTIF_SCORE_REVERSE_STRAND_MNV",
        MOTIF_SCORE_REVERSE_STRAND_MNV,
    ),
    ("SECOND_EXPECTED_DIVERGENCE", SECOND_EXPECTED_DIVERGENCE),
    ("ARCHIVE_DBSNP_SYNONYMS", ARCHIVE_DBSNP_SYNONYMS),
    ("ORDER_ONLY", ORDER_ONLY),
    ("INDIVIDUAL_ALLELE_ORDER", INDIVIDUAL_ALLELE_ORDER),
    (
        "CONSERVATION_SHORT_NAME_SPELLING",
        CONSERVATION_SHORT_NAME_SPELLING,
    ),
    ("PC_DESCRIPTION_TEXT", PC_DESCRIPTION_TEXT),
    ("START_COOCCURRENCE_SWAP", START_COOCCURRENCE_SWAP),
    ("SPLICE_FAMILY_SWAP", SPLICE_FAMILY_SWAP),
    (
        "GXF_REVERSE_STRAND_UTR_SEQUENCE",
        GXF_REVERSE_STRAND_UTR_SEQUENCE,
    ),
    (
        "ALTERNATE_CODON_PAST_THE_CODING_END",
        ALTERNATE_CODON_PAST_THE_CODING_END,
    ),
    (
        "ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE",
        ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE,
    ),
    ("CONSERVATION_MAX_TEXT_ORDER", CONSERVATION_MAX_TEXT_ORDER),
    ("GWAS_HAPLOTYPE_CROSS_PRODUCT", GWAS_HAPLOTYPE_CROSS_PRODUCT),
    (
        "GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE",
        GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE,
    ),
    (
        "ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION",
        ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION,
    ),
    ("CADD_STRUCTURAL_ALLELE_MATCH", CADD_STRUCTURAL_ALLELE_MATCH),
    (
        "SPLICEVAULT_SIBLING_ALLELE_ROWS",
        SPLICEVAULT_SIBLING_ALLELE_ROWS,
    ),
    ("SPLICEVAULT_JSON_LATER_ROWS", SPLICEVAULT_JSON_LATER_ROWS),
    (
        "UTRANNOTATOR_GAINED_START_WINDOW",
        UTRANNOTATOR_GAINED_START_WINDOW,
    ),
    (
        "UTRANNOTATOR_CONSEQUENCE_PRECEDENCE",
        UTRANNOTATOR_CONSEQUENCE_PRECEDENCE,
    ),
    (
        "UTRANNOTATOR_ANNOTATION_PAIR_ORDER",
        UTRANNOTATOR_ANNOTATION_PAIR_ORDER,
    ),
    (
        "UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT",
        UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT,
    ),
    (
        "UTRANNOTATOR_LOST_START_AFTER_INSERTION",
        UTRANNOTATOR_LOST_START_AFTER_INSERTION,
    ),
    (
        "UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE",
        UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE,
    ),
    (
        "UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY",
        UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY,
    ),
    (
        "UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON",
        UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON,
    ),
    (
        "UTRANNOTATOR_WITHOUT_REFERENCE",
        UTRANNOTATOR_WITHOUT_REFERENCE,
    ),
    ("NEAREST_EXON_JB_ITEM_ORDER", NEAREST_EXON_JB_ITEM_ORDER),
    ("LOEUF_ZERO_MINIMUM", LOEUF_ZERO_MINIMUM),
    (
        "MAXENTSCAN_OVERLAPPING_WINDOWS",
        MAXENTSCAN_OVERLAPPING_WINDOWS,
    ),
    (
        "ALPHAMISSENSE_SIBLING_ALLELE_ROWS",
        ALPHAMISSENSE_SIBLING_ALLELE_ROWS,
    ),
    ("REVEL_SIBLING_ALLELE_ROWS", REVEL_SIBLING_ALLELE_ROWS),
    (
        "REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS",
        REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS,
    ),
    (
        "GNOMADMT_INFO_VALUE_SECOND_EQUALS",
        GNOMADMT_INFO_VALUE_SECOND_EQUALS,
    ),
    (
        "LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD",
        LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD,
    ),
    (
        "LOFTEE_FIRST_ALT_DONOR_FILTER",
        LOFTEE_FIRST_ALT_DONOR_FILTER,
    ),
    (
        "LOFTEE_GERP_ZOOM_LEVEL_SUMMARY",
        LOFTEE_GERP_ZOOM_LEVEL_SUMMARY,
    ),
    (
        "LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE",
        LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE,
    ),
    (
        "LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND",
        LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND,
    ),
    (
        "NEAREST_EXON_JB_STRETCH_PRIMING",
        NEAREST_EXON_JB_STRETCH_PRIMING,
    ),
    ("NEAREST_EXON_JB_SHIFTED_SPAN", NEAREST_EXON_JB_SHIFTED_SPAN),
    (
        "NEAREST_EXON_JB_ID_LESS_EXON_FOLD",
        NEAREST_EXON_JB_ID_LESS_EXON_FOLD,
    ),
    (
        "NEAREST_START_SITE_ONE_BASE_ABOVE",
        NEAREST_START_SITE_ONE_BASE_ABOVE,
    ),
    (
        "DOWNSTREAM_REVERSE_STRAND_ROTATION",
        DOWNSTREAM_REVERSE_STRAND_ROTATION,
    ),
    (
        "SPLICEAI_CHOPPED_MULTIALLELIC",
        SPLICEAI_CHOPPED_MULTIALLELIC,
    ),
    ("SPLICEAI_MULTI_ENTRY_CELL", SPLICEAI_MULTI_ENTRY_CELL),
    (
        "PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING",
        PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING,
    ),
    (
        "TABIX_ALLELE_CASE_INSENSITIVE_MATCH",
        TABIX_ALLELE_CASE_INSENSITIVE_MATCH,
    ),
    (
        "DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M",
        DBNSFP_MITOCHONDRIAL_RECORD_NAMED_M,
    ),
    (
        "DOSAGE_SENSITIVITY_REGULATORY_ROW_WARNINGS",
        DOSAGE_SENSITIVITY_REGULATORY_ROW_WARNINGS,
    ),
    ("FORK_WARNING_PID_PREFIX", FORK_WARNING_PID_PREFIX),
    (
        "JSON_MISMATCH_WARNING_NAMES_NO_RECORD",
        JSON_MISMATCH_WARNING_NAMES_NO_RECORD,
    ),
    (
        "REFSEQ_EDITED_MODEL_HGVSC_OFFSET",
        REFSEQ_EDITED_MODEL_HGVSC_OFFSET,
    ),
    (
        "REFSEQ_EDITED_MODEL_HGVS_OFFSET",
        REFSEQ_EDITED_MODEL_HGVS_OFFSET,
    ),
    (
        "REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES",
        REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES,
    ),
    (
        "REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC",
        REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC,
    ),
    (
        "REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC",
        REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC,
    ),
    (
        "HGVS_SHIFTED_INSERTION_CLAMPED_END",
        HGVS_SHIFTED_INSERTION_CLAMPED_END,
    ),
    (
        "HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE",
        HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE,
    ),
    (
        "HGVS_INTRON_OFFSET_ABUTTING_INSERTION",
        HGVS_INTRON_OFFSET_ABUTTING_INSERTION,
    ),
    (
        "HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE",
        HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE,
    ),
    (
        "HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST",
        HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST,
    ),
    (
        "HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE",
        HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE,
    ),
    (
        "HGVS_SHIFTED_DELETION_CLAMPED_SPAN",
        HGVS_SHIFTED_DELETION_CLAMPED_SPAN,
    ),
    (
        "HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE",
        HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE,
    ),
    (
        "HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END",
        HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END,
    ),
    (
        "HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS",
        HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS,
    ),
    (
        "HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT",
        HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT,
    ),
    (
        "HGVSP_STOP_CODON_PROTEIN_COMPARISON",
        HGVSP_STOP_CODON_PROTEIN_COMPARISON,
    ),
    (
        "HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE",
        HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE,
    ),
    (
        "HGVSP_FRAMESHIFT_COUNT_PAST_A_STOP_READ_AS_A_RESIDUE",
        HGVSP_FRAMESHIFT_COUNT_PAST_A_STOP_READ_AS_A_RESIDUE,
    ),
    (
        "HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE",
        HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE,
    ),
    (
        "HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE",
        HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE,
    ),
    ("HGVSP_INTRON_END_INSERTION", HGVSP_INTRON_END_INSERTION),
    ("HGVSP_INTRON_START_INSERTION", HGVSP_INTRON_START_INSERTION),
    (
        "SHIFT_INTRON_TERM_UNSHIFTED_SPAN",
        SHIFT_INTRON_TERM_UNSHIFTED_SPAN,
    ),
    ("MATURE_MIRNA_UNSHIFTED_SPAN", MATURE_MIRNA_UNSHIFTED_SPAN),
    ("STOP_CODON_UNSHIFTED_SPAN", STOP_CODON_UNSHIFTED_SPAN),
    (
        "STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON",
        STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON,
    ),
    (
        "STOP_CODON_COMPLETED_BY_POLYADENYLATION",
        STOP_CODON_COMPLETED_BY_POLYADENYLATION,
    ),
    (
        "STOP_CODON_OVERLAP_WITHOUT_STOP_CODON",
        STOP_CODON_OVERLAP_WITHOUT_STOP_CODON,
    ),
    (
        "STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON",
        STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON,
    ),
    (
        "START_CODON_OVERLAP_WITHOUT_START_CODON",
        START_CODON_OVERLAP_WITHOUT_START_CODON,
    ),
    (
        "STOP_RETAINED_WITHOUT_STOP_CODON_TEST",
        STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
    ),
    (
        "STOP_CODON_OVERLAP_GENOMIC_SPAN",
        STOP_CODON_OVERLAP_GENOMIC_SPAN,
    ),
    (
        "STOP_RETAINED_THREE_BASE_STOP_CODON_TEST",
        STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
    ),
    (
        "STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST",
        STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
    ),
    (
        "STOP_CODON_OVERLAP_INSERTION_RANGE",
        STOP_CODON_OVERLAP_INSERTION_RANGE,
    ),
    (
        "STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES",
        STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES,
    ),
    (
        "FRAMESHIFT_INSIDE_A_HELD_STOP_CODON",
        FRAMESHIFT_INSIDE_A_HELD_STOP_CODON,
    ),
    (
        "START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR",
        START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR,
    ),
    (
        "START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR",
        START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR,
    ),
    (
        "CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN",
        CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN,
    ),
    (
        "SPLICE_TERM_PAST_THE_RECORD_SPAN",
        SPLICE_TERM_PAST_THE_RECORD_SPAN,
    ),
    (
        "REFSEQ_EDITED_MODEL_DELETED_SPAN_TERM",
        REFSEQ_EDITED_MODEL_DELETED_SPAN_TERM,
    ),
    (
        "REFSEQ_EDITED_MODEL_TAIL_TERM",
        REFSEQ_EDITED_MODEL_TAIL_TERM,
    ),
    (
        "REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM",
        REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM,
    ),
    (
        "REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS",
        REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS,
    ),
    (
        "REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE",
        REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE,
    ),
    (
        "REFSEQ_EDITED_MODEL_INSERTED_SPAN_CODONS",
        REFSEQ_EDITED_MODEL_INSERTED_SPAN_CODONS,
    ),
    ("FREQUENCY_FILTER_ZERO_VALUE", FREQUENCY_FILTER_ZERO_VALUE),
    (
        "MERGED_SAME_STABLE_ID_COLLAPSE",
        MERGED_SAME_STABLE_ID_COLLAPSE,
    ),
    (
        "USE_TRANSCRIPT_REF_CACHE_CDNA",
        USE_TRANSCRIPT_REF_CACHE_CDNA,
    ),
    (
        "REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET",
        REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET,
    ),
    ("HGNC_ID_COPY", HGNC_ID_COPY),
    (
        "DOWNSTREAM_UTR_FRAME_PAST_CODING_END",
        DOWNSTREAM_UTR_FRAME_PAST_CODING_END,
    ),
    (
        "DOWNSTREAM_START_PAST_OWN_CODON",
        DOWNSTREAM_START_PAST_OWN_CODON,
    ),
    (
        "PLUGIN_HGVS_SHIFTED_CDS_COORDINATES",
        PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
    ),
    (
        "PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD",
        PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD,
    ),
    ("BREAKEND_MATE_LOCAL_READ", BREAKEND_MATE_LOCAL_READ),
    ("BND_OWN_CHROMOSOME_ROWS_LOST", BND_OWN_CHROMOSOME_ROWS_LOST),
    (
        "TRANSCRIPT_ANNOTATOR_SIBLING_ALLELE_ROWS",
        TRANSCRIPT_ANNOTATOR_SIBLING_ALLELE_ROWS,
    ),
    (
        "HGVS_SHIFT_SHARED_ACROSS_ALLELES",
        HGVS_SHIFT_SHARED_ACROSS_ALLELES,
    ),
    (
        "BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES",
        BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES,
    ),
    (
        "HGVS_SHIFT_SEQUENCE_EDGE_WINDOW",
        HGVS_SHIFT_SEQUENCE_EDGE_WINDOW,
    ),
    (
        "FREQUENCY_FILTER_ALLELE_KEY_ORDER",
        FREQUENCY_FILTER_ALLELE_KEY_ORDER,
    ),
    ("MINIMAL_SYMBOLIC_ALLELES", MINIMAL_SYMBOLIC_ALLELES),
    (
        "SHIFT_3PRIME_SHARED_PREDICATES",
        SHIFT_3PRIME_SHARED_PREDICATES,
    ),
    ("SHIFT_START_CODON_TEST", SHIFT_START_CODON_TEST),
    (
        "TRANSCRIPT_ANNOTATOR_INSERTION_WINDOW",
        TRANSCRIPT_ANNOTATOR_INSERTION_WINDOW,
    ),
    (
        "PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN",
        PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
    ),
    (
        "PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR",
        PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR,
    ),
    ("TOTAL_LENGTH_JSON_RANGE_END", TOTAL_LENGTH_JSON_RANGE_END),
    (
        "ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT",
        ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT,
    ),
    (
        "COLOCATED_LIST_CUT_AT_CACHE_CAP",
        COLOCATED_LIST_CUT_AT_CACHE_CAP,
    ),
];

const TABLE: &str = include_str!("../data/classes.tsv");

/// The Ensembl VEP releases a reference output is read as written by, oldest first: the
/// choices of `--reference-release`.
pub const REFERENCE_RELEASES: &[&str] = &["115.2", "116.2", "116.3"];

/// The release a reference output is read as written by where the command line names none.
pub const DEFAULT_REFERENCE_RELEASE: &str = "116.3";

/// The column of `data/classes.tsv` whose cell is the first Ensembl VEP release that
/// corrects the defect the row's class sets aside, empty for a class no release corrects.
pub const CORRECTED_IN: &str = "CORRECTED_IN";

/// The two numbers of an Ensembl VEP release, `116.3` as `(116, 3)`; `None` for a text that
/// is not two runs of ASCII digits joined by one dot.
pub fn release_numbers(release: &str) -> Option<(u32, u32)> {
    let number = |digits: &str| -> Option<u32> {
        let whole = digits.bytes().all(|digit| digit.is_ascii_digit());
        digits.parse().ok().filter(|_| whole)
    };
    let (version, patch) = release.split_once('.')?;
    Some((number(version)?, number(patch)?))
}

/// The first release that corrects the defect `class` sets aside, where a comparison whose
/// reference output Ensembl VEP `release` wrote does not apply the class: its `CORRECTED_IN`
/// release is `release` or an earlier one. `None` for a class such a comparison applies,
/// and for a `release` that is no release.
pub fn corrected_under(class: &str, release: &str) -> Option<&'static str> {
    let corrected = tables().corrected_in(class)?;
    (release_numbers(corrected)? <= release_numbers(release)?).then_some(corrected)
}

/// The class tables of `data/classes.tsv`, each in its own order.
#[derive(Debug)]
pub struct ClassTables {
    field_classes: Vec<(&'static str, Vec<&'static str>)>,
    registry_classes: Vec<&'static str>,
    term_pair_classes: Vec<&'static str>,
    recoder_masks: Vec<(&'static str, Vec<&'static str>)>,
    plugin_warning_classes: Vec<(&'static str, &'static str)>,
    warning_text_classes: Vec<&'static str>,
    intron_gap_classes: Vec<&'static str>,
    downstream_classes: Vec<&'static str>,
    transcript_reading_classes: Vec<&'static str>,
    structural_variant_classes: Vec<&'static str>,
    cache_reading_classes: Vec<&'static str>,
    manifest_classes: Vec<&'static str>,
    record_dropping_classes: Vec<&'static str>,
    record_reading_classes: Vec<&'static str>,
    allele_reading_classes: Vec<&'static str>,
    sibling_allele_classes: Vec<&'static str>,
    term_pair_following_classes: Vec<&'static str>,
    record_reference_classes: Vec<&'static str>,
    front_end_whole_classes: Vec<&'static str>,
    corrected_classes: Vec<(&'static str, &'static str)>,
}

/// The class tables, read once.
pub fn tables() -> &'static ClassTables {
    static TABLES: LazyLock<ClassTables> = LazyLock::new(ClassTables::read);
    &TABLES
}

/// The place of the column `name` in the header line of `data/classes.tsv`.
fn column_place(header: &[&str], name: &str) -> usize {
    header
        .iter()
        .position(|cell| *cell == name)
        .unwrap_or_else(|| panic!("data/classes.tsv lacks the column {name}"))
}

/// The classes whose cell of the column `name` is not empty, each with that cell, in the
/// file's order.
fn cells(name: &str) -> Vec<(&'static str, &'static str)> {
    let mut lines = TABLE.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let place = column_place(&header, name);
    lines
        .filter_map(|line| {
            let cells: Vec<&'static str> = line.split('\t').collect();
            let cell = cells.get(place).copied().filter(|cell| !cell.is_empty())?;
            Some((cells[0], cell))
        })
        .collect()
}

/// The classes of the column `table`, each with the cell of the column `values`, in the order
/// the column `table` numbers them.
fn column(table: &str, values: Option<&str>) -> Vec<(&'static str, &'static str)> {
    let mut lines = TABLE.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let find = |name: &str| column_place(&header, name);
    let place = find(table);
    let value = values.map(find);
    let mut found: Vec<(usize, &'static str, &'static str)> = Vec::new();
    for line in lines {
        let cells: Vec<&'static str> = line.split('\t').collect();
        if let Some(Ok(position)) = cells.get(place).map(|cell| cell.parse::<usize>()) {
            let cell = value.and_then(|at| cells.get(at).copied()).unwrap_or("");
            found.push((position, cells[0], cell));
        }
    }
    found.sort();
    found
        .into_iter()
        .map(|(_, class, cell)| (class, cell))
        .collect()
}

impl ClassTables {
    fn read() -> ClassTables {
        ClassTables {
            field_classes: column("FIELD_CLASSES", Some("FIELD_CLASSES.values"))
                .into_iter()
                .map(|(class, cell)| {
                    (
                        class,
                        cell.split(',').filter(|item| !item.is_empty()).collect(),
                    )
                })
                .collect(),
            registry_classes: column("REGISTRY_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            term_pair_classes: column("TERM_PAIR_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            recoder_masks: column("RECODER_MASKS", Some("RECODER_MASKS.values"))
                .into_iter()
                .map(|(class, cell)| {
                    (
                        class,
                        cell.split(',').filter(|item| !item.is_empty()).collect(),
                    )
                })
                .collect(),
            plugin_warning_classes: column(
                "PLUGIN_WARNING_CLASSES",
                Some("PLUGIN_WARNING_CLASSES.value"),
            ),
            warning_text_classes: column("WARNING_TEXT_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            intron_gap_classes: column("INTRON_GAP_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            downstream_classes: column("DOWNSTREAM_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            transcript_reading_classes: column("TRANSCRIPT_READING_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            structural_variant_classes: column("STRUCTURAL_VARIANT_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            cache_reading_classes: column("CACHE_READING_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            manifest_classes: column("MANIFEST_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            record_dropping_classes: column("RECORD_DROPPING_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            record_reading_classes: column("RECORD_READING_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            allele_reading_classes: column("ALLELE_READING_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            sibling_allele_classes: column("SIBLING_ALLELE_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            term_pair_following_classes: column("TERM_PAIR_FOLLOWING_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            record_reference_classes: column("RECORD_REFERENCE_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            front_end_whole_classes: column("FRONT_END_WHOLE_CLASSES", None)
                .into_iter()
                .map(|(class, _)| class)
                .collect(),
            corrected_classes: cells(CORRECTED_IN)
                .into_iter()
                .inspect(|(class, release)| {
                    assert!(
                        release_numbers(release).is_some(),
                        "data/classes.tsv: {CORRECTED_IN} of {class} is no release: {release}"
                    );
                })
                .collect(),
        }
    }

    /// `CORRECTED_IN`: each class whose defect an Ensembl VEP release corrects, with the
    /// first release that corrects it, in the file's order.
    pub fn corrected_classes(&self) -> &[(&'static str, &'static str)] {
        &self.corrected_classes
    }

    /// The first Ensembl VEP release that corrects the defect `class` sets aside; `None` for
    /// a class no release corrects.
    pub fn corrected_in(&self, class: &str) -> Option<&'static str> {
        self.corrected_classes
            .iter()
            .find(|(held, _)| *held == class)
            .map(|(_, release)| *release)
    }

    /// `FIELD_CLASSES`: each class with its items, in the table's order.
    pub fn field_classes(&self) -> &[(&'static str, Vec<&'static str>)] {
        &self.field_classes
    }

    /// The items `FIELD_CLASSES` holds for `class`; `None` for a class it lacks.
    pub fn field_classes_of(&self, class: &str) -> Option<&[&'static str]> {
        self.field_classes
            .iter()
            .find(|(held, _)| *held == class)
            .map(|(_, items)| items.as_slice())
    }

    /// Whether `class` is a key of `FIELD_CLASSES`.
    pub fn in_field_classes(&self, class: &str) -> bool {
        self.field_classes_of(class).is_some()
    }

    /// `REGISTRY_CLASSES`, in its order.
    pub fn registry_classes(&self) -> &[&'static str] {
        &self.registry_classes
    }

    /// Whether `class` is one of `REGISTRY_CLASSES`.
    pub fn in_registry_classes(&self, class: &str) -> bool {
        self.registry_classes.contains(&class)
    }

    /// `TERM_PAIR_CLASSES`, in its order.
    pub fn term_pair_classes(&self) -> &[&'static str] {
        &self.term_pair_classes
    }

    /// Whether `class` is one of `TERM_PAIR_CLASSES`.
    pub fn in_term_pair_classes(&self, class: &str) -> bool {
        self.term_pair_classes.contains(&class)
    }

    /// `RECODER_MASKS`: each class with its items, in the table's order.
    pub fn recoder_masks(&self) -> &[(&'static str, Vec<&'static str>)] {
        &self.recoder_masks
    }

    /// The items `RECODER_MASKS` holds for `class`; `None` for a class it lacks.
    pub fn recoder_masks_of(&self, class: &str) -> Option<&[&'static str]> {
        self.recoder_masks
            .iter()
            .find(|(held, _)| *held == class)
            .map(|(_, items)| items.as_slice())
    }

    /// Whether `class` is a key of `RECODER_MASKS`.
    pub fn in_recoder_masks(&self, class: &str) -> bool {
        self.recoder_masks_of(class).is_some()
    }

    /// `PLUGIN_WARNING_CLASSES`: each class with its text, in the table's order.
    pub fn plugin_warning_classes(&self) -> &[(&'static str, &'static str)] {
        &self.plugin_warning_classes
    }

    /// Whether `class` is a key of `PLUGIN_WARNING_CLASSES`.
    pub fn in_plugin_warning_classes(&self, class: &str) -> bool {
        self.plugin_warning_classes
            .iter()
            .any(|(held, _)| *held == class)
    }

    /// `WARNING_TEXT_CLASSES`, in its order.
    pub fn warning_text_classes(&self) -> &[&'static str] {
        &self.warning_text_classes
    }

    /// Whether `class` is one of `WARNING_TEXT_CLASSES`.
    pub fn in_warning_text_classes(&self, class: &str) -> bool {
        self.warning_text_classes.contains(&class)
    }

    /// `INTRON_GAP_CLASSES`, in its order.
    pub fn intron_gap_classes(&self) -> &[&'static str] {
        &self.intron_gap_classes
    }

    /// Whether `class` is one of `INTRON_GAP_CLASSES`.
    pub fn in_intron_gap_classes(&self, class: &str) -> bool {
        self.intron_gap_classes.contains(&class)
    }

    /// `DOWNSTREAM_CLASSES`, in its order.
    pub fn downstream_classes(&self) -> &[&'static str] {
        &self.downstream_classes
    }

    /// Whether `class` is one of `DOWNSTREAM_CLASSES`.
    pub fn in_downstream_classes(&self, class: &str) -> bool {
        self.downstream_classes.contains(&class)
    }

    /// `TRANSCRIPT_READING_CLASSES`, in its order.
    pub fn transcript_reading_classes(&self) -> &[&'static str] {
        &self.transcript_reading_classes
    }

    /// Whether `class` is one of `TRANSCRIPT_READING_CLASSES`.
    pub fn in_transcript_reading_classes(&self, class: &str) -> bool {
        self.transcript_reading_classes.contains(&class)
    }

    /// `STRUCTURAL_VARIANT_CLASSES`, in its order.
    pub fn structural_variant_classes(&self) -> &[&'static str] {
        &self.structural_variant_classes
    }

    /// Whether `class` is one of `STRUCTURAL_VARIANT_CLASSES`.
    pub fn in_structural_variant_classes(&self, class: &str) -> bool {
        self.structural_variant_classes.contains(&class)
    }

    /// `CACHE_READING_CLASSES`, in its order.
    pub fn cache_reading_classes(&self) -> &[&'static str] {
        &self.cache_reading_classes
    }

    /// Whether `class` is one of `CACHE_READING_CLASSES`.
    pub fn in_cache_reading_classes(&self, class: &str) -> bool {
        self.cache_reading_classes.contains(&class)
    }

    /// `MANIFEST_CLASSES`, in its order.
    pub fn manifest_classes(&self) -> &[&'static str] {
        &self.manifest_classes
    }

    /// Whether `class` is one of `MANIFEST_CLASSES`.
    pub fn in_manifest_classes(&self, class: &str) -> bool {
        self.manifest_classes.contains(&class)
    }

    /// `RECORD_DROPPING_CLASSES`, in its order.
    pub fn record_dropping_classes(&self) -> &[&'static str] {
        &self.record_dropping_classes
    }

    /// Whether `class` is one of `RECORD_DROPPING_CLASSES`.
    pub fn in_record_dropping_classes(&self, class: &str) -> bool {
        self.record_dropping_classes.contains(&class)
    }

    /// `RECORD_READING_CLASSES`, in its order.
    pub fn record_reading_classes(&self) -> &[&'static str] {
        &self.record_reading_classes
    }

    /// Whether `class` is one of `RECORD_READING_CLASSES`.
    pub fn in_record_reading_classes(&self, class: &str) -> bool {
        self.record_reading_classes.contains(&class)
    }

    /// `ALLELE_READING_CLASSES`, in its order.
    pub fn allele_reading_classes(&self) -> &[&'static str] {
        &self.allele_reading_classes
    }

    /// Whether `class` is one of `ALLELE_READING_CLASSES`.
    pub fn in_allele_reading_classes(&self, class: &str) -> bool {
        self.allele_reading_classes.contains(&class)
    }

    /// `SIBLING_ALLELE_CLASSES`, in its order.
    pub fn sibling_allele_classes(&self) -> &[&'static str] {
        &self.sibling_allele_classes
    }

    /// Whether `class` is one of `SIBLING_ALLELE_CLASSES`.
    pub fn in_sibling_allele_classes(&self, class: &str) -> bool {
        self.sibling_allele_classes.contains(&class)
    }

    /// `TERM_PAIR_FOLLOWING_CLASSES`, in its order.
    pub fn term_pair_following_classes(&self) -> &[&'static str] {
        &self.term_pair_following_classes
    }

    /// Whether `class` is one of `TERM_PAIR_FOLLOWING_CLASSES`.
    pub fn in_term_pair_following_classes(&self, class: &str) -> bool {
        self.term_pair_following_classes.contains(&class)
    }

    /// `RECORD_REFERENCE_CLASSES`, in its order.
    pub fn record_reference_classes(&self) -> &[&'static str] {
        &self.record_reference_classes
    }

    /// Whether `class` is one of `RECORD_REFERENCE_CLASSES`.
    pub fn in_record_reference_classes(&self, class: &str) -> bool {
        self.record_reference_classes.contains(&class)
    }

    /// `FRONT_END_WHOLE_CLASSES`, in its order.
    pub fn front_end_whole_classes(&self) -> &[&'static str] {
        &self.front_end_whole_classes
    }

    /// Whether `class` is one of `FRONT_END_WHOLE_CLASSES`.
    pub fn in_front_end_whole_classes(&self, class: &str) -> bool {
        self.front_end_whole_classes.contains(&class)
    }
}
